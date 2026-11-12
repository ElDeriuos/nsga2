!> EXPECTED TO FAIL: without stat, invalid input must stop the program with a message.
program xfail_invalid_input_without_stat
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use nsga2
  use benchmarks, only: zdt1
  implicit none
  type(nsga2_result) :: res

  call nsga2_optimize(zdt1, 2, [0._dp, 0._dp], [1._dp, 1._dp], nsga2_options(population_size=1), res)
  print '(a)', 'not reached: invalid input was accepted'
end program xfail_invalid_input_without_stat
