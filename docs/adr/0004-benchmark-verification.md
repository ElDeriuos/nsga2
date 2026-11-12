# 4. How the benchmarks are verified

Date: 2026-11-12 · Status: accepted

**Problem definitions** (22 problems) were checked against independent sources: the pymoo source for ZDT1-4, ZDT6, DTLZ1-7, KUR, BNH, SRN, TNK, OSY, and Wikipedia's *Test functions for optimization* for SCH, FON, POL, VNT, CONSTR.

**Reference fronts** are independent of the library: analytic where known (SCH, FON, ZDT, DTLZ), otherwise brute force, i.e. the non-dominated feasible points of a dense grid over the decision space (POL, KUR, VNT, BNH, SRN, TNK, CONSTR). OSY (6 variables, 6 constraints) is too large for a grid; its test requires a hypervolume above 1.2x that of random search with 8x the evaluations.

**Limits** are about 2x the worst IGD of 5 calibration seeds (default settings; DTLZ1, DTLZ3, DTLZ6 at 1000 generations, where 250 is clearly not converged). Each test runs one of those seeds.

**Negative control:** `test/xfail/xfail_no_evolution` must fail, proving the IGD checks can detect a run that does not evolve.
