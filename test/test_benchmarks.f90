!> NSGA-II on the standard multi-objective benchmarks. Each test: one seeded run with default
!> settings (population 100, 250 generations unless noted), all returned solutions feasible, and
!> IGD to the true front below a limit set at about 2x the worst of 5 calibration seeds.
!> Reference fronts are analytic where known, else brute force over a grid of the decision space.
program test_benchmarks
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use nsga2
  use benchmarks
  use test_check, only: check, seed_rng
  implicit none

  real(dp), parameter :: pi = acos(-1._dp)
  real(dp), parameter :: zero(30) = 0, one(30) = 1

  ! --- result contract (ZDT1)
  block
    type(nsga2_result) :: res
    real(dp) :: f(2)
    integer :: i
    logical :: same
    call seed_rng(1)
    call nsga2_optimize(zdt1, 2, zero, one, nsga2_options(), res)
    call check(all(nondominated_rank(res%f) == 1), 'the returned solutions are mutually non-dominated')
    same = .true.
    do i = 1, size(res%f, 2)
      call zdt1(res%x(:, i), f)
      same = same .and. all(f == res%f(:, i))
    end do
    call check(same .and. size(res%x, 2) == size(res%f, 2), 'res%f holds the objective values of res%x')
  end block

  ! --- classic 2-objective problems (Deb et al. 2002) and Viennet (3 objectives)
  call solve('SCH', sch, 2, [-1000._dp], [1000._dp], front_sch(1000), 0.04_dp, 1)
  call solve('FON', fon, 2, spread(-4._dp, 1, 3), spread(4._dp, 1, 3), front_fon(1000), 0.012_dp, 2)
  call solve('POL', pol, 2, [-pi, -pi], [pi, pi], grid_front(pol, 2, [-pi, -pi], [pi, pi], 500), 0.15_dp, 3)
  call solve('KUR', kur, 2, spread(-5._dp, 1, 3), spread(5._dp, 1, 3), &
             grid_front(kur, 2, spread(-5._dp, 1, 3), spread(5._dp, 1, 3), 100), 0.25_dp, 4)
  call solve('VNT', vnt, 3, [-3._dp, -3._dp], [3._dp, 3._dp], grid_front(vnt, 3, [-3._dp, -3._dp], [3._dp, 3._dp], 100), &
             0.09_dp, 5)

  ! --- ZDT (Zitzler, Deb, Thiele 2000)
  call solve('ZDT1', zdt1, 2, zero, one, front_zdt1(1000), 0.011_dp, 1)
  call solve('ZDT2', zdt2, 2, zero, one, front_zdt2(1000), 0.011_dp, 2)
  call solve('ZDT3', zdt3, 2, zero, one, front_zdt3(10000), 0.07_dp, 3)
  call solve('ZDT4', zdt4, 2, [0._dp, spread(-5._dp, 1, 9)], [1._dp, spread(5._dp, 1, 9)], front_zdt1(1000), 0.02_dp, 4)
  call solve('ZDT6', zdt6, 2, zero(:10), one(:10), front_zdt6(1000), 0.02_dp, 5)

  ! --- DTLZ, 3 objectives (Deb, Thiele, Laumanns, Zitzler); DTLZ1/3/6 need 1000 generations
  call solve('DTLZ1', dtlz1, 3, zero(:7), one(:7), front_dtlz1(40), 0.06_dp, 1, 1000)
  call solve('DTLZ2', dtlz2, 3, zero(:12), one(:12), front_sphere(40), 0.14_dp, 2)
  call solve('DTLZ3', dtlz3, 3, zero(:12), one(:12), front_sphere(40), 0.15_dp, 3, 1000)
  call solve('DTLZ4', dtlz4, 3, zero(:12), one(:12), front_sphere(40), 0.15_dp, 4)
  call solve('DTLZ5', dtlz5, 3, zero(:12), one(:12), front_dtlz5(1000), 0.012_dp, 5)
  call solve('DTLZ6', dtlz6, 3, zero(:12), one(:12), front_dtlz5(1000), 0.21_dp, 1, 1000)
  call solve('DTLZ7', dtlz7, 3, zero(:22), one(:22), front_dtlz7(100), 0.17_dp, 2)

  ! --- constrained (Deb's constrained dominance)
  call solve('BNH', bnh, 2, [0._dp, 0._dp], [5._dp, 3._dp], &
             grid_front(bnh, 2, [0._dp, 0._dp], [5._dp, 3._dp], 500, bnh_g, 2), 1.1_dp, 3, g=bnh_g, nc=2)
  call solve('SRN', srn, 2, [-20._dp, -20._dp], [20._dp, 20._dp], &
             grid_front(srn, 2, [-20._dp, -20._dp], [20._dp, 20._dp], 500, srn_g, 2), 2.5_dp, 4, g=srn_g, nc=2)
  call solve('TNK', tnk, 2, [0._dp, 0._dp], [pi, pi], &
             grid_front(tnk, 2, [0._dp, 0._dp], [pi, pi], 500, tnk_g, 2), 0.011_dp, 5, g=tnk_g, nc=2)
  call solve('CONSTR', constr, 2, [0.1_dp, 0._dp], [1._dp, 5._dp], &
             grid_front(constr, 2, [0.1_dp, 0._dp], [1._dp, 5._dp], 500, constr_g, 2), 0.045_dp, 1, g=constr_g, nc=2)

  ! --- OSY: 6 variables, 6 constraints, too many for a grid; must clearly beat random search
  !     that uses 8x the evaluations (hypervolume, reference point (0, 100))
  block
    real(dp), parameter :: lo(6) = [0, 0, 1, 0, 1, 0], hi(6) = [10, 10, 5, 6, 5, 10]
    real(dp), allocatable :: pts(:,:)
    real(dp) :: x(6), g(6), hv_random
    type(nsga2_result) :: res
    integer :: j, m
    call seed_rng(99)
    allocate(pts(2, 200000))
    m = 0
    do j = 1, 200000
      call random_number(x)
      x = lo + x * (hi - lo)
      call osy_g(x, g)
      if (all(g <= 0)) then
        m = m + 1
        call osy(x, pts(:, m))
      end if
    end do
    hv_random = hypervolume_2d(nondominated_points(pts(:, :m)), [0._dp, 100._dp])
    call seed_rng(2)
    call nsga2_optimize(osy, 2, lo, hi, nsga2_options(), res, constraints=osy_g, n_con=6)
    call check(all(res%cv == 0) .and. hypervolume_2d(res%f, [0._dp, 100._dp]) > 1.2_dp * hv_random, &
               'OSY: all feasible, hypervolume > 1.2x that of random search with 8x the evaluations')
  end block

contains

  subroutine solve(name, f, n_obj, lower, upper, ref, limit, seed, generations, g, nc)
    character(*), intent(in) :: name
    procedure(objective_function) :: f
    integer, intent(in) :: n_obj, seed
    real(dp), intent(in) :: lower(:), upper(:), ref(:,:), limit
    integer, intent(in), optional :: generations, nc
    procedure(constraint_function), optional :: g
    type(nsga2_options) :: opts
    type(nsga2_result) :: res
    character(16) :: txt

    if (present(generations)) opts%max_generations = generations
    call seed_rng(seed)
    call nsga2_optimize(f, n_obj, lower, upper, opts, res, constraints=g, n_con=nc)
    write(txt, '(f6.3)') limit
    call check(all(res%cv == 0) .and. igd(res%f, ref) < limit, &
               name//': all feasible, IGD < '//trim(adjustl(txt))//' to the true front')
  end subroutine solve

end program test_benchmarks
