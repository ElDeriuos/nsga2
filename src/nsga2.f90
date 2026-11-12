module nsga2
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: dominates, nondominated_rank, crowding_distance, crowded_tournament
  public :: sbx_crossover, polynomial_mutation
  public :: igd, hypervolume_2d
  public :: nsga2_optimize, objective_function, constraint_function

  !> Settings for one run; every field has a usable default (Deb et al. 2002).
  type, public :: nsga2_options
    integer  :: population_size = 100
    integer  :: max_generations = 250
    real(dp) :: crossover_rate  = 0.9_dp   !< probability that a parent pair is recombined
    real(dp) :: eta_c           = 20       !< SBX distribution index
    real(dp) :: mutation_rate   = -1       !< per variable; < 0 means 1/(number of variables)
    real(dp) :: eta_m           = 20       !< polynomial mutation distribution index
  end type nsga2_options

  !> The non-dominated (rank-1) solutions of the final population.
  type, public :: nsga2_result
    integer :: generations = 0
    real(dp), allocatable :: x(:,:)   !< (n_var, n_front) decision vectors
    real(dp), allocatable :: f(:,:)   !< (n_obj, n_front) objective values
    real(dp), allocatable :: cv(:)    !< (n_front) constraint violation, 0 = feasible
  end type nsga2_result

  abstract interface
    !> Fill f(1:n_obj) with the objective values of x; all objectives are minimized.
    subroutine objective_function(x, f)
      import :: dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f(:)
    end subroutine objective_function

    !> Fill g(1:n_con) for x; x is feasible when every g <= 0.
    subroutine constraint_function(x, g)
      import :: dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
    end subroutine constraint_function
  end interface

contains

  !> NSGA-II (Deb, Pratap, Agarwal, Meyarivan 2002) on x in [lower, upper] with n_obj objectives,
  !> optionally subject to n_con constraints g(x) <= 0 (Deb's constrained dominance).
  !> Invalid input: sets stat /= 0 and errmsg if stat is present, else error stop.
  subroutine nsga2_optimize(objectives, n_obj, lower, upper, opts, res, constraints, n_con, stat, errmsg)
    procedure(objective_function) :: objectives
    integer, intent(in) :: n_obj
    real(dp), intent(in) :: lower(:), upper(:)
    type(nsga2_options), intent(in) :: opts
    type(nsga2_result), intent(out) :: res
    procedure(constraint_function), optional :: constraints
    integer, intent(in), optional :: n_con
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    real(dp), allocatable :: x(:,:), f(:,:), cv(:), crowd(:), xa(:,:), fa(:,:), cva(:), crowd_all(:)
    integer, allocatable :: rank(:), rank_all(:), keep(:)
    character(:), allocatable :: msg
    real(dp) :: pm
    integer :: n, np, g, i

    msg = invalid_input()
    if (present(stat)) stat = merge(1, 0, len(msg) > 0)
    if (len(msg) > 0) then
      if (.not. present(stat)) error stop 'nsga2: '//msg
      if (present(errmsg)) errmsg = msg
      return
    end if

    n = size(lower)
    np = opts%population_size
    pm = opts%mutation_rate
    if (pm < 0) pm = 1._dp / n
    allocate(x(n, np), f(n_obj, np), cv(np), xa(n, 2*np), fa(n_obj, 2*np), cva(2*np))

    call random_number(x)
    x = spread(lower, 2, np) + x * spread(upper - lower, 2, np)
    call evaluate(x, f, cv)
    rank = nondominated_rank(f, cv)
    crowd = crowding_by_front(f, rank)

    do g = 1, opts%max_generations
      xa(:, :np) = x
      fa(:, :np) = f
      cva(:np) = cv
      xa(:, np+1:) = offspring()
      call evaluate(xa(:, np+1:), fa(:, np+1:), cva(np+1:))
      rank_all = nondominated_rank(fa, cva)
      crowd_all = crowding_by_front(fa, rank_all)
      ! elitist survival: best np of parents + children by rank, then by crowding distance
      keep = sort_index(rank_all + 0.5_dp / (1 + crowd_all))
      keep = keep(:np)
      x = xa(:, keep)
      f = fa(:, keep)
      cv = cva(keep)
      rank = rank_all(keep)
      crowd = crowd_all(keep)
    end do

    keep = pack([(i, i = 1, np)], rank == 1)
    res%generations = opts%max_generations
    res%x = x(:, keep)
    res%f = f(:, keep)
    res%cv = cv(keep)

  contains

    !> Empty string when input is valid, else what is wrong (names the argument or field).
    function invalid_input() result(msg)
      character(:), allocatable :: msg
      msg = ''
      if (n_obj < 1) then
        msg = 'n_obj must be >= 1'
      else if (size(lower) < 1) then
        msg = 'need at least one decision variable (size(lower) >= 1)'
      else if (size(lower) /= size(upper)) then
        msg = 'lower and upper must have the same size'
      else if (any(lower > upper)) then
        msg = 'lower must be <= upper'
      else if (present(constraints) .neqv. present(n_con)) then
        msg = 'constraints and n_con must be given together'
      else if (opts%population_size < 2) then
        msg = 'population_size must be >= 2'
      else if (opts%max_generations < 0) then
        msg = 'max_generations must be >= 0'
      else if (opts%crossover_rate < 0 .or. opts%crossover_rate > 1) then
        msg = 'crossover_rate must be in [0, 1]'
      else if (opts%mutation_rate > 1) then
        msg = 'mutation_rate must be <= 1 (< 0 means 1/variables)'
      else if (opts%eta_c < 0) then
        msg = 'eta_c must be >= 0'
      else if (opts%eta_m < 0) then
        msg = 'eta_m must be >= 0'
      end if
      if (len(msg) == 0 .and. present(n_con)) then
        if (n_con < 1) msg = 'n_con must be >= 1'
      end if
    end function invalid_input

    subroutine evaluate(xs, fs, cvs)
      real(dp), intent(in) :: xs(:,:)
      real(dp), intent(out) :: fs(:,:), cvs(:)
      integer :: j
      do j = 1, size(xs, 2)
        call objectives(xs(:, j), fs(:, j))
      end do
      cvs = 0
      ! passed on, not called here: gfortran 16 wrongly flags a call of an optional dummy procedure
      ! from an internal procedure as an implicit interface (-Wimplicit-interface)
      if (present(constraints)) cvs = violation(constraints, n_con, xs)
    end subroutine evaluate

    !> np children: crowded binary tournaments, SBX on pairs, polynomial mutation.
    function offspring() result(children)
      real(dp) :: children(n, np), r(2, np + 1), u, uc(n), mc(n), sc(n), c1(n), c2(n)
      integer :: parents(np + 1), j
      call random_number(r)
      parents = crowded_tournament(rank, crowd, min(np, 1 + int(r * np)))
      do j = 1, np, 2
        call random_number(u)
        if (u < opts%crossover_rate) then
          call random_number(uc)
          call random_number(mc)
          call random_number(sc)
          call sbx_crossover(x(:, parents(j)), x(:, parents(j+1)), mc < 0.5_dp, uc, lower, upper, &
                             opts%eta_c, c1, c2, swap=sc < 0.5_dp)
        else
          c1 = x(:, parents(j))
          c2 = x(:, parents(j+1))
        end if
        call mutate(c1)
        call mutate(c2)
        children(:, j) = c1
        if (j < np) children(:, j+1) = c2
      end do
    end function offspring

    subroutine mutate(c)
      real(dp), intent(inout) :: c(:)
      real(dp) :: um(n), r(n)
      call random_number(um)
      call random_number(r)
      call polynomial_mutation(c, um < pm, r, lower, upper, opts%eta_m)
    end subroutine mutate

  end subroutine nsga2_optimize

  !> Total constraint violation sum(max(0, g)) of each column of xs.
  function violation(constraints, n_con, xs) result(cv)
    procedure(constraint_function) :: constraints
    integer, intent(in) :: n_con
    real(dp), intent(in) :: xs(:,:)
    real(dp) :: cv(size(xs, 2)), g(n_con)
    integer :: j
    do j = 1, size(xs, 2)
      call constraints(xs(:, j), g)
      cv(j) = sum(max(0._dp, g))
    end do
  end function violation

  !> Crowding distance computed separately within each front.
  pure function crowding_by_front(f, rank) result(crowd)
    real(dp), intent(in) :: f(:,:)
    integer, intent(in) :: rank(:)
    real(dp) :: crowd(size(rank))
    integer, allocatable :: members(:)
    integer :: k, i

    do k = 1, maxval(rank)
      members = pack([(i, i = 1, size(rank))], rank == k)
      if (size(members) > 0) crowd(members) = crowding_distance(f(:, members))
    end do
  end function crowding_by_front

  !> Inverted generational distance: mean Euclidean distance from each reference point
  !> (columns of reference) to the nearest point of front. Lower is better.
  pure real(dp) function igd(front, reference)
    real(dp), intent(in) :: front(:,:), reference(:,:)
    integer :: j

    igd = 0
    do j = 1, size(reference, 2)
      igd = igd + sqrt(minval(sum((front - spread(reference(:, j), 2, size(front, 2)))**2, dim=1)))
    end do
    igd = igd / size(reference, 2)
  end function igd

  !> Area dominated by a 2-objective front (minimization) and bounded by ref_point.
  !> Points not strictly better than ref_point in both objectives are ignored.
  pure real(dp) function hypervolume_2d(front, ref_point)
    real(dp), intent(in) :: front(:,:), ref_point(2)
    integer :: idx(size(front, 2))
    real(dp) :: lowest
    integer :: k

    idx = sort_index(front(1, :))
    hypervolume_2d = 0
    lowest = ref_point(2)
    do k = 1, size(idx)
      associate (p => front(:, idx(k)))
        if (p(1) < ref_point(1) .and. p(2) < lowest) then
          hypervolume_2d = hypervolume_2d + (ref_point(1) - p(1)) * (lowest - p(2))
          lowest = p(2)
        end if
      end associate
    end do
  end function hypervolume_2d

  !> Polynomial mutation, bounded version (Deb). Genes where mask is true move by a step drawn
  !> from r in [0,1): r = 0 reaches the lower bound, r = 0.5 no change, r -> 1 the upper bound.
  !> Larger eta gives smaller steps.
  pure subroutine polynomial_mutation(x, mask, r, lower, upper, eta)
    real(dp), intent(inout) :: x(:)
    logical, intent(in) :: mask(:)
    real(dp), intent(in) :: r(:), lower(:), upper(:), eta
    real(dp) :: span, val, dq
    integer :: i

    do i = 1, size(x)
      span = upper(i) - lower(i)
      if (.not. mask(i) .or. span <= 0) cycle
      if (r(i) < 0.5_dp) then
        val = 2 * r(i) + (1 - 2 * r(i)) * (1 - (x(i) - lower(i)) / span)**(eta + 1)
        dq = val**(1 / (eta + 1)) - 1
      else
        val = 2 * (1 - r(i)) + 2 * (r(i) - 0.5_dp) * (1 - (upper(i) - x(i)) / span)**(eta + 1)
        dq = 1 - val**(1 / (eta + 1))
      end if
      x(i) = min(max(x(i) + dq * span, lower(i)), upper(i))
    end do
  end subroutine polynomial_mutation

  !> Simulated binary crossover, bounded version (Deb & Agrawal 1995; as in Deb's NSGA-II code).
  !> Only genes where mask is true are crossed, with one uniform draw u in [0,1) per gene;
  !> the child near p1 goes to c1, unless swap is true for that gene (Deb's code swaps each
  !> crossed gene with probability 0.5, which is what mixes the two parents).
  !> Larger eta keeps children closer to their parents.
  pure subroutine sbx_crossover(p1, p2, mask, u, lower, upper, eta, c1, c2, swap)
    real(dp), intent(in) :: p1(:), p2(:), u(:), lower(:), upper(:), eta
    logical, intent(in) :: mask(:)
    real(dp), intent(out) :: c1(:), c2(:)
    logical, intent(in), optional :: swap(:)
    real(dp) :: y1, y2, d, lo_child, hi_child
    integer :: i

    c1 = p1
    c2 = p2
    do i = 1, size(p1)
      if (.not. mask(i) .or. abs(p1(i) - p2(i)) <= 1e-14_dp) cycle
      y1 = min(p1(i), p2(i))
      y2 = max(p1(i), p2(i))
      d = y2 - y1
      lo_child = 0.5_dp * ((y1 + y2) - spread_factor(1 + 2 * (y1 - lower(i)) / d) * d)
      hi_child = 0.5_dp * ((y1 + y2) + spread_factor(1 + 2 * (upper(i) - y2) / d) * d)
      lo_child = min(max(lo_child, lower(i)), upper(i))
      hi_child = min(max(hi_child, lower(i)), upper(i))
      if (p1(i) <= p2(i)) then
        c1(i) = lo_child
        c2(i) = hi_child
      else
        c1(i) = hi_child
        c2(i) = lo_child
      end if
      if (present(swap)) then
        if (swap(i)) then
          d = c1(i)
          c1(i) = c2(i)
          c2(i) = d
        end if
      end if
    end do
  contains
    pure real(dp) function spread_factor(beta)
      real(dp), intent(in) :: beta
      real(dp) :: alpha
      alpha = 2 - beta**(-(eta + 1))
      if (u(i) <= 1 / alpha) then
        spread_factor = (u(i) * alpha)**(1 / (eta + 1))
      else
        spread_factor = (1 / (2 - u(i) * alpha))**(1 / (eta + 1))
      end if
    end function spread_factor
  end subroutine sbx_crossover

  !> Binary tournaments: each column of contestants (2 x n_picks) is one tournament.
  !> Lower rank wins; on equal rank the larger crowding distance wins; then the first contestant.
  pure function crowded_tournament(rank, crowd, contestants) result(winners)
    integer, intent(in) :: rank(:), contestants(:,:)
    real(dp), intent(in) :: crowd(:)
    integer :: winners(size(contestants, 2)), j, a, b
    logical :: b_wins

    do j = 1, size(contestants, 2)
      a = contestants(1, j)
      b = contestants(2, j)
      b_wins = rank(b) < rank(a) .or. (rank(b) == rank(a) .and. crowd(b) > crowd(a))
      winners(j) = merge(b, a, b_wins)
    end do
  end function crowded_tournament

  !> Crowding distance of the points of ONE front, f(n_obj, n): per objective, the gap between
  !> each point's two neighbours divided by that objective's range; boundary points get huge().
  pure function crowding_distance(f) result(d)
    real(dp), intent(in) :: f(:,:)
    real(dp) :: d(size(f, 2)), span
    integer :: idx(size(f, 2)), m, i, n

    n = size(f, 2)
    d = 0
    if (n <= 2) then
      d = huge(1._dp)
      return
    end if
    do m = 1, size(f, 1)
      idx = sort_index(f(m, :))
      span = f(m, idx(n)) - f(m, idx(1))
      d(idx(1)) = huge(1._dp)
      d(idx(n)) = huge(1._dp)
      if (span <= 0) cycle
      do i = 2, n - 1
        if (d(idx(i)) < huge(1._dp)) d(idx(i)) = d(idx(i)) + (f(m, idx(i+1)) - f(m, idx(i-1))) / span
      end do
    end do
  end function crowding_distance

  !> Indices that sort a ascending (stable).
  pure function sort_index(a) result(idx)
    real(dp), intent(in) :: a(:)
    integer :: idx(size(a)), i, j, t

    ! shortcut: insertion sort, O(n^2) on fronts of at most 2 x population_size; use merge sort if populations grow large
    idx = [(i, i = 1, size(a))]
    do i = 2, size(a)
      t = idx(i)
      j = i - 1
      do while (j >= 1)
        if (a(idx(j)) <= a(t)) exit
        idx(j+1) = idx(j)
        j = j - 1
      end do
      idx(j+1) = t
    end do
  end function sort_index

  !> Fast non-dominated sort (Deb et al. 2002): rank 1 = not dominated by anyone,
  !> rank k = dominated only by ranks < k. f is (n_obj, n); cv is the constraint violation.
  pure function nondominated_rank(f, cv) result(rank)
    real(dp), intent(in) :: f(:,:)
    real(dp), intent(in), optional :: cv(:)
    integer :: rank(size(f, 2))
    real(dp) :: v(size(f, 2))
    integer :: dominated_by(size(f, 2)), i, j, k
    logical :: beats(size(f, 2), size(f, 2))   ! beats(i, j): i dominates j

    v = 0
    if (present(cv)) v = cv
    do j = 1, size(f, 2)
      do i = 1, size(f, 2)
        beats(i, j) = dominates(f(:, i), f(:, j), v(i), v(j))
      end do
      dominated_by(j) = count(beats(:, j))
    end do
    rank = 0
    k = 0
    do while (any(rank == 0))
      k = k + 1
      where (rank == 0 .and. dominated_by == 0) rank = -k     ! mark the whole current front first
      do i = 1, size(f, 2)
        if (rank(i) == -k) dominated_by = dominated_by - merge(1, 0, beats(i, :))
      end do
      where (rank == -k) rank = k
    end do
  end function nondominated_rank

  !> Does a dominate b (all objectives minimized)? With violations (cv >= 0, 0 = feasible),
  !> Deb's rule applies: feasible beats infeasible, smaller violation beats larger.
  pure logical function dominates(fa, fb, cva, cvb)
    real(dp), intent(in) :: fa(:), fb(:)
    real(dp), intent(in), optional :: cva, cvb
    real(dp) :: va, vb

    va = 0
    vb = 0
    if (present(cva)) va = cva
    if (present(cvb)) vb = cvb
    if (va > 0 .or. vb > 0) then
      dominates = va < vb
    else
      dominates = all(fa <= fb) .and. any(fa < fb)
    end if
  end function dominates

end module nsga2
