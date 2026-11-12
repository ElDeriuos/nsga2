!> Invalid input is rejected through stat/errmsg instead of stopping the caller.
program test_errors
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use nsga2
  use benchmarks, only: zdt1, bnh, bnh_g
  use test_check, only: check
  implicit none

  type(nsga2_result) :: res
  integer :: stat
  character(:), allocatable :: msg
  real(dp), parameter :: lo(3) = 0, hi(3) = 1

  call nsga2_optimize(zdt1, 2, lo, hi, nsga2_options(max_generations=1), res, stat=stat, errmsg=msg)
  call check(stat == 0, 'valid input gives stat = 0')

  call nsga2_optimize(zdt1, 2, lo, hi, nsga2_options(population_size=1), res, stat=stat, errmsg=msg)
  call check(rejected('population_size'), 'population_size < 2 is rejected')

  call nsga2_optimize(zdt1, 2, lo, hi, nsga2_options(max_generations=-1), res, stat=stat, errmsg=msg)
  call check(rejected('max_generations'), 'negative max_generations is rejected')

  call nsga2_optimize(zdt1, 0, lo, hi, nsga2_options(), res, stat=stat, errmsg=msg)
  call check(rejected('n_obj'), 'n_obj < 1 is rejected')

  call nsga2_optimize(zdt1, 2, lo(:0), hi(:0), nsga2_options(), res, stat=stat, errmsg=msg)
  call check(rejected('variable'), 'no decision variables is rejected')

  call nsga2_optimize(zdt1, 2, lo, hi(:2), nsga2_options(), res, stat=stat, errmsg=msg)
  call check(rejected('size'), 'lower and upper of different size are rejected')

  call nsga2_optimize(zdt1, 2, hi, lo, nsga2_options(), res, stat=stat, errmsg=msg)
  call check(rejected('lower'), 'lower > upper is rejected')

  call nsga2_optimize(zdt1, 2, lo, hi, nsga2_options(crossover_rate=1.5_dp), res, stat=stat, errmsg=msg)
  call check(rejected('crossover_rate'), 'crossover_rate outside [0,1] is rejected')

  call nsga2_optimize(zdt1, 2, lo, hi, nsga2_options(mutation_rate=2._dp), res, stat=stat, errmsg=msg)
  call check(rejected('mutation_rate'), 'mutation_rate > 1 is rejected')

  call nsga2_optimize(zdt1, 2, lo, hi, nsga2_options(eta_c=-1._dp), res, stat=stat, errmsg=msg)
  call check(rejected('eta_c'), 'negative eta_c is rejected')

  call nsga2_optimize(zdt1, 2, lo, hi, nsga2_options(eta_m=-1._dp), res, stat=stat, errmsg=msg)
  call check(rejected('eta_m'), 'negative eta_m is rejected')

  call nsga2_optimize(bnh, 2, lo(:2), hi(:2), nsga2_options(), res, constraints=bnh_g, stat=stat, errmsg=msg)
  call check(rejected('n_con'), 'constraints without n_con are rejected')

  call nsga2_optimize(bnh, 2, lo(:2), hi(:2), nsga2_options(), res, constraints=bnh_g, n_con=0, stat=stat, errmsg=msg)
  call check(rejected('n_con'), 'n_con < 1 is rejected')

contains

  logical function rejected(field)
    character(*), intent(in) :: field
    rejected = stat /= 0 .and. index(msg, field) > 0
  end function rejected

end program test_errors
