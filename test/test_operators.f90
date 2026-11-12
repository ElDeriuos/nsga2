program test_operators
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use nsga2
  use test_check, only: check
  implicit none

  logical, parameter :: T = .true., F = .false.

  ! --- dominance (minimization) and Deb's constrained dominance
  call check(dominates([1._dp, 2._dp], [1._dp, 3._dp]), 'better in one objective, equal in the other: dominates')
  call check(.not. dominates([1._dp, 3._dp], [2._dp, 2._dp]), 'a trade-off does not dominate')
  call check(.not. dominates([1._dp, 2._dp], [1._dp, 2._dp]), 'equal solutions do not dominate each other')
  call check(dominates([9._dp, 9._dp], [1._dp, 1._dp], 0._dp, 0.5_dp), 'feasible dominates infeasible, whatever the objectives')
  call check(dominates([9._dp, 9._dp], [1._dp, 1._dp], 0.2_dp, 0.5_dp), 'of two infeasible, smaller violation dominates')
  call check(.not. dominates([1._dp, 3._dp], [2._dp, 2._dp], 0._dp, 0._dp), 'two feasible compare by objectives')

  ! --- non-dominated sorting: A(1,5) B(2,3) C(3,1) D(2,4) E(4,4) F(5,5)
  !     fronts {A,B,C}, then {D} (B dominates it), then {E} (D dominates it), then {F}
  block
    real(dp) :: f(2, 6)
    f = reshape([1,5, 2,3, 3,1, 2,4, 4,4, 5,5], [2, 6])
    call check(all(nondominated_rank(f) == [1, 1, 1, 2, 3, 4]), 'non-dominated sort finds the fronts')
    call check(all(nondominated_rank(f, [0._dp, 0._dp, 1._dp, 0._dp, 0._dp, 0._dp]) == [1, 1, 5, 2, 3, 4]), &
               'an infeasible solution ranks behind every feasible one')
  end block

  ! --- crowding distance: front (1,5) (2,3) (3,1); ranges 2 and 4; middle = 2/2 + 4/4 = 2, ends infinite
  block
    real(dp) :: d(3)
    d = crowding_distance(reshape([1._dp,5._dp, 2._dp,3._dp, 3._dp,1._dp], [2, 3]))
    call check(d(1) == huge(1._dp) .and. d(2) == 2 .and. d(3) == huge(1._dp), &
               'crowding distance: ends infinite, middle sums normalized neighbour gaps')
    d = crowding_distance(reshape([2._dp,3._dp, 1._dp,5._dp, 3._dp,1._dp], [2, 3]))
    call check(d(1) == 2 .and. d(2) == huge(1._dp) .and. d(3) == huge(1._dp), &
               'crowding distance does not depend on the order of the points')
    d = crowding_distance(reshape([1._dp,1._dp, 1._dp,2._dp, 1._dp,3._dp], [2, 3]))
    call check(d(2) == 1, 'an objective with zero range adds nothing (no division by zero)')
  end block

  ! --- crowded tournament: lower rank wins, then larger crowding distance, then the first contestant
  call check(all(crowded_tournament([1, 2, 1, 1], [.5_dp, 9._dp, 2._dp, .5_dp], &
                                    reshape([1,2, 1,3, 4,1, 2,3], [2, 4])) == [1, 3, 4, 3]), &
             'crowded tournament: rank first, then crowding distance')

  ! --- SBX (bounded): expected values worked out independently from Deb's formulas
  block
    real(dp) :: c1(3), c2(3)
    call sbx_crossover([.25_dp, .4_dp, .1_dp], [.75_dp, .4_dp, .9_dp], [T, T, F], [.25_dp, .7_dp, .3_dp], &
                       [0._dp, 0._dp, 0._dp], [1._dp, 1._dp, 1._dp], 1._dp, c1, c2)
    call check(near(c1, [0.334640543058463_dp, .4_dp, .1_dp]) .and. near(c2, [0.665359456941537_dp, .4_dp, .9_dp]), &
               'SBX matches the worked example; equal and unmasked genes are copied')
    call sbx_crossover([.4_dp], [.2_dp], [T], [.8_dp], [0._dp], [1._dp], 2._dp, c1(:1), c2(:1))
    call check(near(c1(:1), [0.435458109832015_dp]) .and. near(c2(:1), [0.167473745192818_dp]), &
               'SBX puts the child near each parent in that parent''s slot')
    call sbx_crossover([.25_dp, .4_dp, .1_dp], [.75_dp, .4_dp, .9_dp], [T, T, F], [.25_dp, .7_dp, .3_dp], &
                       [0._dp, 0._dp, 0._dp], [1._dp, 1._dp, 1._dp], 1._dp, c1, c2, swap=[T, F, T])
    call check(near(c1, [0.665359456941537_dp, .4_dp, .1_dp]) .and. near(c2, [0.334640543058463_dp, .4_dp, .9_dp]), &
               'SBX swap exchanges crossed genes between the children (this mixes the parents)')
  end block

  ! --- polynomial mutation (bounded): r = 0 -> lower bound, r = 0.5 -> no change; worked examples
  block
    real(dp) :: x(3)
    x = [.3_dp, .3_dp, .3_dp]
    call polynomial_mutation(x, [T, T, F], [0._dp, .5_dp, .9_dp], [0._dp, 0._dp, 0._dp], [1._dp, 1._dp, 1._dp], 20._dp)
    call check(near(x, [0._dp, .3_dp, .3_dp]), 'polynomial mutation: r = 0 hits the lower bound, r = 0.5 changes nothing')
    x(:2) = [.3_dp, .3_dp]
    call polynomial_mutation(x(:2), [T, T], [.25_dp, .9_dp], [0._dp, -1._dp], [1._dp, 1._dp], 1._dp)
    call check(near(x(:2), [0.163133825081603_dp, 0.833030334328620_dp]), 'polynomial mutation matches the worked example')
  end block

  ! --- two-point crossover (binary genes, as in Deb's NSGA-II code): genes cut1+1 .. cut2 are exchanged
  block
    real(dp) :: c1(6), c2(6)
    call two_point_crossover([1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp], [0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp], &
                             2, 4, c1, c2)
    call check(all(c1 == [1, 1, 0, 0, 1, 1]) .and. all(c2 == [0, 0, 1, 1, 0, 0]), &
               'two-point crossover exchanges genes cut1+1 .. cut2')
    call two_point_crossover([1._dp, 0._dp, 1._dp], [0._dp, 1._dp, 1._dp], 3, 3, c1(:3), c2(:3))
    call check(all(c1(:3) == [1, 0, 1]) .and. all(c2(:3) == [0, 1, 1]), 'two-point crossover with cut1 = cut2 copies the parents')
  end block

  ! --- quality indicators
  call check(near([igd(reshape([0._dp, 0._dp], [2, 1]), reshape([0._dp,0._dp, 3._dp,4._dp], [2, 2]))], [2.5_dp]), &
             'IGD = mean distance from each reference point to the nearest solution: (0 + 5)/2')
  call check(near([hypervolume_2d(reshape([1._dp,3._dp, 2._dp,2._dp, 3._dp,1._dp], [2, 3]), [4._dp, 4._dp])], [6._dp]), &
             'hypervolume of a 3-point staircase: 3 + 2 + 1')
  call check(near([hypervolume_2d(reshape([3._dp,1._dp, 3._dp,3._dp, 1._dp,3._dp, 2._dp,2._dp, 5._dp,0._dp], [2, 5]), &
                                  [4._dp, 4._dp])], [6._dp]), &
             'hypervolume ignores dominated points and points outside the reference box')

contains

  logical function near(a, b)
    real(dp), intent(in) :: a(:), b(:)
    near = size(a) == size(b)
    if (near) near = all(abs(a - b) < 1e-12_dp)
  end function near

end program test_operators
