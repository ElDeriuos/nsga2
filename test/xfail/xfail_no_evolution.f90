!> EXPECTED TO FAIL (negative control): with evolution switched off, a random population is far
!> from the ZDT1 front, so this check must fail. If it ever passes, the benchmark assertions are
!> not able to detect a broken NSGA-II.
program xfail_no_evolution
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use nsga2
  use benchmarks, only: zdt1, front_zdt1
  use test_check, only: check, seed_rng
  implicit none
  type(nsga2_result) :: res

  call seed_rng(1)
  call nsga2_optimize(zdt1, 2, spread(0._dp, 1, 30), spread(1._dp, 1, 30), nsga2_options(max_generations=0), res)
  call check(igd(res%f, front_zdt1(1000)) < 0.011_dp, 'ZDT1 solved without any generation')
end program xfail_no_evolution
