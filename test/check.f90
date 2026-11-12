module test_check
  implicit none
  private
  public :: check, seed_rng

contains

  subroutine check(condition, name)
    logical, intent(in) :: condition
    character(*), intent(in) :: name
    if (.not. condition) error stop 'FAILED: '//name
    print '(a)', 'ok    '//name
  end subroutine check

  !> Fix the intrinsic RNG so stochastic tests are reproducible.
  subroutine seed_rng(value)
    integer, intent(in) :: value
    integer :: n, i
    call random_seed(size=n)
    call random_seed(put=[(value + 7919*i, i = 1, n)])
  end subroutine seed_rng

end module test_check
