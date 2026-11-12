!> Standard multi-objective benchmark problems (all objectives minimized, constraints g <= 0
!> feasible) and their reference Pareto fronts. Module procedures, not internal ones: passing
!> internal procedures makes gfortran build a trampoline, which needs an executable stack.
!>
!> Sources: SCH, FON, POL, KUR, VNT and the constrained set as used in Deb et al. 2002 (NSGA-II)
!> and Deb's NSGA-II code; ZDT: Zitzler, Deb & Thiele 2000; DTLZ: Deb, Thiele, Laumanns & Zitzler
!> 2002/2005 (3 objectives). Definitions are cross-checked against pymoo (see README).
module benchmarks
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: sch, fon, pol, kur, vnt
  public :: zdt1, zdt2, zdt3, zdt4, zdt6
  public :: dtlz1, dtlz2, dtlz3, dtlz4, dtlz5, dtlz6, dtlz7
  public :: bnh, bnh_g, srn, srn_g, tnk, tnk_g, constr, constr_g, osy, osy_g
  public :: front_sch, front_fon, front_zdt1, front_zdt2, front_zdt3, front_zdt6
  public :: front_dtlz1, front_sphere, front_dtlz5, front_dtlz7
  public :: grid_front, nondominated_points
  public :: one_min_max, lotz, zdt5, knapsack, knapsack_g, front_lotz, front_zdt5, front_knapsack
  public :: gear, front_gear

  real(dp), parameter :: pi = acos(-1._dp)

  !> Knapsack instance: profit and weight of item j in knapsack i, random integers in [10, 100]
  !> drawn with Python's random.Random(1999) (profits row by row, then weights).
  integer, parameter :: kp_profit(20, 2) = reshape([ &
    24, 55, 28, 82, 73, 42, 89, 82, 35, 31, 86, 36, 92, 16, 61, 60, 21, 89, 16, 28, &
    45, 17, 100, 32, 56, 75, 76, 28, 14, 18, 78, 77, 28, 17, 40, 83, 44, 13, 77, 65], [20, 2])
  integer, parameter :: kp_weight(20, 2) = reshape([ &
    48, 14, 31, 49, 16, 25, 70, 70, 86, 84, 24, 55, 11, 57, 75, 81, 27, 67, 43, 48, &
    19, 65, 29, 100, 36, 45, 89, 88, 37, 95, 69, 36, 19, 72, 48, 80, 47, 84, 26, 85], [20, 2])

  abstract interface
    subroutine vector_function(x, y)
      import :: dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y(:)
    end subroutine vector_function
  end interface

contains

  ! ================================================================ classic 2-objective (unconstrained)

  !> SCH (Schaffer): 1 variable in [-1000, 1000].
  subroutine sch(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = x(1)**2
    f(2) = (x(1) - 2)**2
  end subroutine sch

  !> FON (Fonseca & Fleming): 3 variables in [-4, 4].
  subroutine fon(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = 1 - exp(-sum((x - 1 / sqrt(3._dp))**2))
    f(2) = 1 - exp(-sum((x + 1 / sqrt(3._dp))**2))
  end subroutine fon

  !> POL (Poloni): 2 variables in [-pi, pi]; disconnected front.
  subroutine pol(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    real(dp), parameter :: a1 = 0.5_dp*sin(1._dp) - 2*cos(1._dp) + sin(2._dp) - 1.5_dp*cos(2._dp)
    real(dp), parameter :: a2 = 1.5_dp*sin(1._dp) - cos(1._dp) + 2*sin(2._dp) - 0.5_dp*cos(2._dp)
    real(dp) :: b1, b2
    b1 = 0.5_dp*sin(x(1)) - 2*cos(x(1)) + sin(x(2)) - 1.5_dp*cos(x(2))
    b2 = 1.5_dp*sin(x(1)) - cos(x(1)) + 2*sin(x(2)) - 0.5_dp*cos(x(2))
    f(1) = 1 + (a1 - b1)**2 + (a2 - b2)**2
    f(2) = (x(1) + 3)**2 + (x(2) + 1)**2
  end subroutine pol

  !> KUR (Kursawe): 3 variables in [-5, 5]; disconnected front.
  subroutine kur(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = sum(-10 * exp(-0.2_dp * sqrt(x(:size(x)-1)**2 + x(2:)**2)))
    f(2) = sum(abs(x)**0.8_dp + 5 * sin(x**3))
  end subroutine kur

  !> VNT (Viennet): 2 variables in [-3, 3], 3 objectives.
  subroutine vnt(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    real(dp) :: r2
    r2 = x(1)**2 + x(2)**2
    f(1) = 0.5_dp * r2 + sin(r2)
    f(2) = (3*x(1) - 2*x(2) + 4)**2 / 8 + (x(1) - x(2) + 1)**2 / 27 + 15
    f(3) = 1 / (r2 + 1) - 1.1_dp * exp(-r2)
  end subroutine vnt

  ! ================================================================ ZDT (2 objectives)

  pure real(dp) function zdt_g(x)
    real(dp), intent(in) :: x(:)
    zdt_g = 1 + 9 * sum(x(2:)) / (size(x) - 1)
  end function zdt_g

  !> ZDT1: 30 variables in [0, 1]; convex front.
  subroutine zdt1(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    real(dp) :: g
    g = zdt_g(x)
    f(1) = x(1)
    f(2) = g * (1 - sqrt(x(1) / g))
  end subroutine zdt1

  !> ZDT2: 30 variables in [0, 1]; concave front.
  subroutine zdt2(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    real(dp) :: g
    g = zdt_g(x)
    f(1) = x(1)
    f(2) = g * (1 - (x(1) / g)**2)
  end subroutine zdt2

  !> ZDT3: 30 variables in [0, 1]; disconnected front.
  subroutine zdt3(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    real(dp) :: g
    g = zdt_g(x)
    f(1) = x(1)
    f(2) = g * (1 - sqrt(x(1) / g) - x(1) / g * sin(10 * pi * x(1)))
  end subroutine zdt3

  !> ZDT4: 10 variables, x1 in [0, 1], others in [-5, 5]; 21^9 local fronts.
  subroutine zdt4(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    real(dp) :: g
    g = 1 + 10 * (size(x) - 1) + sum(x(2:)**2 - 10 * cos(4 * pi * x(2:)))
    f(1) = x(1)
    f(2) = g * (1 - sqrt(x(1) / g))
  end subroutine zdt4

  !> ZDT6: 10 variables in [0, 1]; non-uniform, concave front.
  subroutine zdt6(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    real(dp) :: g
    f(1) = 1 - exp(-4 * x(1)) * sin(6 * pi * x(1))**6
    g = 1 + 9 * (sum(x(2:)) / (size(x) - 1))**0.25_dp
    f(2) = g * (1 - (f(1) / g)**2)
  end subroutine zdt6

  ! ================================================================ DTLZ (3 objectives)

  pure real(dp) function g_rastrigin(xm)   ! DTLZ1 / DTLZ3
    real(dp), intent(in) :: xm(:)
    g_rastrigin = 100 * (size(xm) + sum((xm - 0.5_dp)**2 - cos(20 * pi * (xm - 0.5_dp))))
  end function g_rastrigin

  pure function sphere3(x1, x2, g) result(f)   ! DTLZ2-6 shape, angles as fractions of pi/2
    real(dp), intent(in) :: x1, x2, g
    real(dp) :: f(3)
    f(1) = (1 + g) * cos(x1 * pi / 2) * cos(x2 * pi / 2)
    f(2) = (1 + g) * cos(x1 * pi / 2) * sin(x2 * pi / 2)
    f(3) = (1 + g) * sin(x1 * pi / 2)
  end function sphere3

  !> DTLZ1: 7 variables in [0, 1]; linear front sum(f) = 0.5, 11^5 - 1 local fronts.
  subroutine dtlz1(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    real(dp) :: g
    g = g_rastrigin(x(3:))
    f(1) = 0.5_dp * x(1) * x(2) * (1 + g)
    f(2) = 0.5_dp * x(1) * (1 - x(2)) * (1 + g)
    f(3) = 0.5_dp * (1 - x(1)) * (1 + g)
  end subroutine dtlz1

  !> DTLZ2: 12 variables in [0, 1]; spherical front sum(f^2) = 1.
  subroutine dtlz2(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(:3) = sphere3(x(1), x(2), sum((x(3:) - 0.5_dp)**2))
  end subroutine dtlz2

  !> DTLZ3: 12 variables in [0, 1]; spherical front, 3^10 - 1 local fronts.
  subroutine dtlz3(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(:3) = sphere3(x(1), x(2), g_rastrigin(x(3:)))
  end subroutine dtlz3

  !> DTLZ4: 12 variables in [0, 1]; spherical front, biased density (alpha = 100).
  subroutine dtlz4(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(:3) = sphere3(x(1)**100, x(2)**100, sum((x(3:) - 0.5_dp)**2))
  end subroutine dtlz4

  !> DTLZ5: 12 variables in [0, 1]; degenerate front (a curve).
  subroutine dtlz5(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    real(dp) :: g
    g = sum((x(3:) - 0.5_dp)**2)
    f(:3) = sphere3(x(1), (1 + 2 * g * x(2)) / (2 * (1 + g)), g)
  end subroutine dtlz5

  !> DTLZ6: 12 variables in [0, 1]; degenerate front, harder g.
  subroutine dtlz6(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    real(dp) :: g
    g = sum(x(3:)**0.1_dp)
    f(:3) = sphere3(x(1), (1 + 2 * g * x(2)) / (2 * (1 + g)), g)
  end subroutine dtlz6

  !> DTLZ7: 22 variables in [0, 1]; 4 disconnected front regions.
  subroutine dtlz7(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    real(dp) :: g
    g = 1 + 9 * sum(x(3:)) / (size(x) - 2)
    f(1) = x(1)
    f(2) = x(2)
    f(3) = (1 + g) * (3 - sum(f(:2) / (1 + g) * (1 + sin(3 * pi * f(:2)))))
  end subroutine dtlz7

  ! ================================================================ constrained (g <= 0 feasible)

  !> BNH (Binh & Korn): x1 in [0, 5], x2 in [0, 3].
  subroutine bnh(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = 4 * x(1)**2 + 4 * x(2)**2
    f(2) = (x(1) - 5)**2 + (x(2) - 5)**2
  end subroutine bnh

  subroutine bnh_g(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    g(1) = (x(1) - 5)**2 + x(2)**2 - 25
    g(2) = 7.7_dp - (x(1) - 8)**2 - (x(2) + 3)**2
  end subroutine bnh_g

  !> SRN (Srinivas & Deb): 2 variables in [-20, 20].
  subroutine srn(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = 2 + (x(1) - 2)**2 + (x(2) - 1)**2
    f(2) = 9 * x(1) - (x(2) - 1)**2
  end subroutine srn

  subroutine srn_g(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    g(1) = x(1)**2 + x(2)**2 - 225
    g(2) = x(1) - 3 * x(2) + 10
  end subroutine srn_g

  !> TNK (Tanaka): 2 variables in [0, pi]; front on a wavy constraint boundary.
  subroutine tnk(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = x(1)
    f(2) = x(2)
  end subroutine tnk

  subroutine tnk_g(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    g(1) = -x(1)**2 - x(2)**2 + 1 + 0.1_dp * cos(16 * atan2(x(1), x(2)))   ! atan2: defined at x2 = 0
    g(2) = (x(1) - 0.5_dp)**2 + (x(2) - 0.5_dp)**2 - 0.5_dp
  end subroutine tnk_g

  !> CONSTR (Deb): x1 in [0.1, 1], x2 in [0, 5].
  subroutine constr(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = x(1)
    f(2) = (1 + x(2)) / x(1)
  end subroutine constr

  subroutine constr_g(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    g(1) = 6 - (x(2) + 9 * x(1))
    g(2) = 1 - (-x(2) + 9 * x(1))
  end subroutine constr_g

  !> OSY (Osyczka & Kundu): 6 variables, x1,x2,x6 in [0,10], x3,x5 in [1,5], x4 in [0,6]; 6 constraints.
  subroutine osy(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = -(25 * (x(1) - 2)**2 + (x(2) - 2)**2 + (x(3) - 1)**2 + (x(4) - 4)**2 + (x(5) - 1)**2)
    f(2) = sum(x**2)
  end subroutine osy

  subroutine osy_g(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    g(1) = 2 - x(1) - x(2)
    g(2) = x(1) + x(2) - 6
    g(3) = x(2) - x(1) - 2
    g(4) = x(1) - 3 * x(2) - 2
    g(5) = (x(3) - 3)**2 + x(4) - 4
    g(6) = 4 - (x(5) - 3)**2 - x(6)
  end subroutine osy_g

  ! ================================================================ bit strings

  !> OneMinMax (Giel & Lehre 2010): f = (number of ones, number of zeros); every string is Pareto optimal.
  subroutine one_min_max(bits, f)
    logical, intent(in) :: bits(:)
    real(dp), intent(out) :: f(:)
    f(1) = count(bits)
    f(2) = size(bits) - count(bits)
  end subroutine one_min_max

  !> LOTZ, Leading Ones Trailing Zeros (Laumanns, Thiele & Zitzler 2004), as minimization:
  !> f = (n - leading ones, n - trailing zeros). Pareto set: 1...10...0.
  subroutine lotz(bits, f)
    logical, intent(in) :: bits(:)
    real(dp), intent(out) :: f(:)
    integer :: n
    n = size(bits)
    f(1) = n - (findloc([bits, .false.], .false., dim=1) - 1)          ! sentinels: all ones / all zeros
    f(2) = n - (n + 1 - findloc([.true., bits], .true., dim=1, back=.true.))
  end subroutine lotz

  !> ZDT5 (deceptive): 80 bits = x1 (30 bits) + ten 5-bit blocks; u = number of ones.
  !> f1 = 1 + u(x1), f2 = g / f1 with g = sum over the blocks of (2 + u if u < 5, else 1).
  subroutine zdt5(bits, f)
    logical, intent(in) :: bits(:)
    real(dp), intent(out) :: f(:)
    integer :: i, u
    real(dp) :: g
    g = 0
    do i = 31, size(bits), 5
      u = count(bits(i:i+4))
      g = g + merge(1, 2 + u, u == 5)
    end do
    f(1) = 1 + count(bits(:30))
    f(2) = g / f(1)
  end subroutine zdt5

  !> Bi-objective 0/1 knapsack (Zitzler & Thiele 1999), 20 items, 2 knapsacks: maximize the profit
  !> in each knapsack (here: minimize its negative).
  subroutine knapsack(bits, f)
    logical, intent(in) :: bits(:)
    real(dp), intent(out) :: f(:)
    f(:2) = -matmul(merge(1, 0, bits), kp_profit)
  end subroutine knapsack

  !> Capacity of each knapsack: half its total weight (Zitzler & Thiele 1999).
  subroutine knapsack_g(bits, g)
    logical, intent(in) :: bits(:)
    real(dp), intent(out) :: g(:)
    g(:2) = matmul(merge(1, 0, bits), kp_weight) - 0.5_dp * sum(kp_weight, dim=1)
  end subroutine knapsack_g

  ! ================================================================ integers

  !> Gear train (Deb, Pratap & Moitra 2000): teeth k in [12, 60]^4; minimize the squared error to
  !> the gear ratio 1/6.931 and the size of the largest gear.
  subroutine gear(k, f)
    integer, intent(in) :: k(:)
    real(dp), intent(out) :: f(:)
    f(1) = (1 / 6.931_dp - real(k(1) * k(2), dp) / (k(3) * k(4)))**2
    f(2) = maxval(k)
  end subroutine gear

  ! ================================================================ analytic reference fronts

  !> Gear train: the exact front, by enumerating all 49^4 designs (smallest error per largest gear).
  function front_gear() result(r)
    real(dp), allocatable :: r(:,:)
    real(dp) :: best(2, 12:60), f(2)
    integer :: a, b, c, d
    best(1, :) = huge(1._dp)
    best(2, :) = [(real(a, dp), a = 12, 60)]
    do a = 12, 60
      do b = 12, 60
        do c = 12, 60
          do d = 12, 60
            call gear([a, b, c, d], f)
            best(1, nint(f(2))) = min(best(1, nint(f(2))), f(1))
          end do
        end do
      end do
    end do
    r = nondominated_points(best)
  end function front_gear

  !> Knapsack: the exact front, by enumerating all 2^20 item subsets.
  function front_knapsack() result(r)
    real(dp), allocatable :: r(:,:), pts(:,:)
    real(dp) :: g(2)
    logical :: bits(20)
    integer :: mask, m, i
    allocate(pts(2, 2**20))
    m = 0
    do mask = 0, 2**20 - 1
      bits = [(btest(mask, i), i = 0, 19)]
      call knapsack_g(bits, g)
      if (all(g <= 0)) then
        m = m + 1
        call knapsack(bits, pts(:, m))
      end if
    end do
    r = nondominated_points(pts(:, :m))
  end function front_knapsack

  !> ZDT5: g = 10, so f2 = 10 / f1 for f1 = 1 .. 31 (31 points).
  function front_zdt5() result(r)
    real(dp) :: r(2, 31)
    integer :: k
    r(1, :) = [(real(k, dp), k = 1, 31)]
    r(2, :) = 10 / r(1, :)
  end function front_zdt5

  !> OneMinMax and LOTZ on n bits: the n + 1 points (k, n - k).
  function front_lotz(n) result(r)
    integer, intent(in) :: n
    real(dp) :: r(2, n + 1)
    integer :: k
    r(1, :) = [(real(k, dp), k = 0, n)]
    r(2, :) = n - r(1, :)
  end function front_lotz

  !> SCH: x in [0, 2] -> (x^2, (x-2)^2).
  function front_sch(n) result(r)
    integer, intent(in) :: n
    real(dp) :: r(2, n), t(n)
    t = linspace(0._dp, 2._dp, n)
    r(1, :) = t**2
    r(2, :) = (t - 2)**2
  end function front_sch

  !> FON: x1 = x2 = x3 = t in [-1/sqrt3, 1/sqrt3].
  function front_fon(n) result(r)
    integer, intent(in) :: n
    real(dp) :: r(2, n), t(n)
    t = linspace(-1 / sqrt(3._dp), 1 / sqrt(3._dp), n)
    r(1, :) = 1 - exp(-3 * (t - 1 / sqrt(3._dp))**2)
    r(2, :) = 1 - exp(-3 * (t + 1 / sqrt(3._dp))**2)
  end function front_fon

  !> ZDT1 and ZDT4: f2 = 1 - sqrt(f1).
  function front_zdt1(n) result(r)
    integer, intent(in) :: n
    real(dp) :: r(2, n)
    r(1, :) = linspace(0._dp, 1._dp, n)
    r(2, :) = 1 - sqrt(r(1, :))
  end function front_zdt1

  !> ZDT2: f2 = 1 - f1^2.
  function front_zdt2(n) result(r)
    integer, intent(in) :: n
    real(dp) :: r(2, n)
    r(1, :) = linspace(0._dp, 1._dp, n)
    r(2, :) = 1 - r(1, :)**2
  end function front_zdt2

  !> ZDT3: the non-dominated part of f2 = 1 - sqrt(f1) - f1 sin(10 pi f1).
  function front_zdt3(n) result(r)
    integer, intent(in) :: n
    real(dp), allocatable :: r(:,:)
    real(dp) :: c(2, n)
    c(1, :) = linspace(0._dp, 1._dp, n)
    c(2, :) = 1 - sqrt(c(1, :)) - c(1, :) * sin(10 * pi * c(1, :))
    r = nondominated_points(c)
  end function front_zdt3

  !> ZDT6: f2 = 1 - f1^2 for f1 in [0.2807753191, 1] (the smallest f1 reachable).
  function front_zdt6(n) result(r)
    integer, intent(in) :: n
    real(dp) :: r(2, n)
    r(1, :) = linspace(0.2807753191_dp, 1._dp, n)
    r(2, :) = 1 - r(1, :)**2
  end function front_zdt6

  !> DTLZ1: the triangle f >= 0, sum(f) = 0.5, sampled on a grid with k divisions per edge.
  function front_dtlz1(k) result(r)
    integer, intent(in) :: k
    real(dp), allocatable :: r(:,:)
    integer :: i, j, m
    allocate(r(3, (k + 1) * (k + 2) / 2))
    m = 0
    do i = 0, k
      do j = 0, k - i
        m = m + 1
        r(:, m) = 0.5_dp * [real(i, dp), real(j, dp), real(k - i - j, dp)] / k
      end do
    end do
  end function front_dtlz1

  !> DTLZ2-4: the octant of the unit sphere, k x k angles.
  function front_sphere(k) result(r)
    integer, intent(in) :: k
    real(dp) :: r(3, k * k), t(k)
    integer :: i, j
    t = linspace(0._dp, 1._dp, k)
    do i = 1, k
      do j = 1, k
        r(:, (i - 1) * k + j) = sphere3(t(i), t(j), 0._dp)
      end do
    end do
  end function front_sphere

  !> DTLZ5-6: the quarter circle through (0,0,1) and (1/sqrt2, 1/sqrt2, 0).
  function front_dtlz5(n) result(r)
    integer, intent(in) :: n
    real(dp) :: r(3, n), t(n)
    integer :: i
    t = linspace(0._dp, 1._dp, n)
    do i = 1, n
      r(:, i) = sphere3(t(i), 0.5_dp, 0._dp)
    end do
  end function front_dtlz5

  !> DTLZ7: g = 1 surface over a k x k grid of (f1, f2), keeping its non-dominated part.
  function front_dtlz7(k) result(r)
    integer, intent(in) :: k
    real(dp), allocatable :: r(:,:)
    real(dp) :: c(3, k * k), t(k)
    integer :: i, j
    t = linspace(0._dp, 1._dp, k)
    do i = 1, k
      do j = 1, k
        c(:2, (i - 1) * k + j) = [t(i), t(j)]
        c(3, (i - 1) * k + j) = 2 * (3 - sum(c(:2, (i - 1) * k + j) / 2 * (1 + sin(3 * pi * c(:2, (i - 1) * k + j)))))
      end do
    end do
    r = nondominated_points(c)
  end function front_dtlz7

  ! ================================================================ brute-force reference fronts

  !> Evaluate a k-per-variable grid over [lower, upper], keep the feasible points and return the
  !> non-dominated ones (n_obj, m). Independent of the library under test.
  function grid_front(obj, n_obj, lower, upper, k, con, n_con) result(front)
    procedure(vector_function) :: obj
    integer, intent(in) :: n_obj, k
    real(dp), intent(in) :: lower(:), upper(:)
    procedure(vector_function), optional :: con
    integer, intent(in), optional :: n_con
    real(dp), allocatable :: front(:,:), pts(:,:), g(:)
    real(dp) :: x(size(lower))
    integer :: idx(size(lower)), n, m, total, j

    allocate(g(1))
    if (present(n_con)) then
      deallocate(g)
      allocate(g(n_con))
    end if
    n = size(lower)
    total = k**n
    allocate(pts(n_obj, total))
    m = 0
    idx = 0
    do j = 1, total
      x = lower + (upper - lower) * idx / real(k - 1, dp)
      if (present(con)) then
        call con(x, g)
      else
        g = 0
      end if
      if (all(g <= 0)) then
        m = m + 1
        call obj(x, pts(:, m))
      end if
      call next_index(idx, k)
    end do
    front = nondominated_points(pts(:, :m))
  end function grid_front

  !> Odometer increment of a mixed-radix index.
  subroutine next_index(idx, k)
    integer, intent(inout) :: idx(:)
    integer, intent(in) :: k
    integer :: i
    do i = 1, size(idx)
      idx(i) = idx(i) + 1
      if (idx(i) < k) return
      idx(i) = 0
    end do
  end subroutine next_index

  !> Non-dominated subset of points (n_obj, n). Two objectives: sort by f1 and sweep, O(n log n).
  !> More objectives: pairwise check, O(n^2) (only used on small sets).
  function nondominated_points(p) result(front)
    real(dp), intent(in) :: p(:,:)
    real(dp), allocatable :: front(:,:)
    logical :: keep(size(p, 2))
    integer, allocatable :: order(:)
    real(dp) :: best
    integer :: i, j

    keep = .false.
    if (size(p, 1) == 2) then
      order = merge_sort_index(p(1, :), p(2, :))   ! by f1, ties by f2
      best = huge(1._dp)
      do i = 1, size(order)
        if (p(2, order(i)) < best) then
          keep(order(i)) = .true.
          best = p(2, order(i))
        end if
      end do
    else
      do i = 1, size(p, 2)
        keep(i) = .true.
        do j = 1, size(p, 2)
          if (all(p(:, j) <= p(:, i)) .and. any(p(:, j) < p(:, i))) then
            keep(i) = .false.
            exit
          end if
        end do
      end do
    end if
    front = p(:, pack([(i, i = 1, size(p, 2))], keep))
  end function nondominated_points

  !> Indices sorting (a, b) lexicographically ascending; merge sort, O(n log n).
  recursive function merge_sort_index(a, b) result(idx)
    real(dp), intent(in) :: a(:), b(:)
    integer :: idx(size(a))
    integer, allocatable :: l(:), r(:)
    integer :: i, j, k, h

    if (size(a) <= 1) then
      idx = [(i, i = 1, size(a))]
      return
    end if
    h = size(a) / 2
    l = merge_sort_index(a(:h), b(:h))
    r = merge_sort_index(a(h+1:), b(h+1:)) + h
    i = 1
    j = 1
    do k = 1, size(a)
      if (j > size(r)) then
        idx(k) = l(i); i = i + 1
      else if (i > size(l)) then
        idx(k) = r(j); j = j + 1
      else if (a(l(i)) < a(r(j)) .or. (a(l(i)) == a(r(j)) .and. b(l(i)) <= b(r(j)))) then
        idx(k) = l(i); i = i + 1
      else
        idx(k) = r(j); j = j + 1
      end if
    end do
  end function merge_sort_index

  function linspace(a, b, n) result(t)
    real(dp), intent(in) :: a, b
    integer, intent(in) :: n
    real(dp) :: t(n)
    integer :: i
    t = [(a + (b - a) * (i - 1) / real(n - 1, dp), i = 1, n)]
  end function linspace

end module benchmarks
