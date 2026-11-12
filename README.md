# nsga2

NSGA-II, the elitist multi-objective genetic algorithm of Deb et al. (2002), as one Modern Fortran module (F2018, `real64`, no dependencies), usable from any program via `use nsga2`.

- Three encodings:
  - **Real:** simulated binary crossover (SBX) and polynomial mutation, both bounded.
  - **Binary:** two-point crossover and bit-flip mutation, as in Deb's NSGA-II code.
  - **Integer:** SBX and polynomial mutation followed by rounding, as in pymoo.
- Fast non-dominated sorting, crowding distance, crowded binary tournament, elitist (μ+λ) survival
- Constraints `g(x) <= 0` with Deb's constrained dominance, in every encoding
- Quality indicators: IGD and 2-D hypervolume
- Tested on 27 standard benchmarks:
  - **Real:** SCH, FON, POL, KUR, VNT, ZDT1–4, ZDT6, DTLZ1–7, BNH, SRN, TNK, CONSTR, OSY.
  - **Binary:** OneMinMax, LOTZ, ZDT5, the 0/1 knapsack.
  - **Integer:** the gear train.
  - Every front is plotted in [Pareto fronts](#pareto-fronts).

Terms are defined in [GLOSSARY.md](GLOSSARY.md); design decisions in [docs/adr](docs/adr). A single-objective companion package: [genetic_algorithm](https://github.com/ElDeriuos/genetic_algorithm).

## Install

**fpm:** add to your project's `fpm.toml`, then `fpm run` / `fpm build`:

```toml
[dependencies]
nsga2 = { git = "https://github.com/ElDeriuos/nsga2", tag = "v0.2.0" }
```

**CMake:** `add_subdirectory(nsga2)` then `target_link_libraries(your_app PRIVATE nsga2)`.

**Plain compiler:** compile `src/nsga2.f90` before your own files: `gfortran nsga2.f90 your_code.f90 -o your_app`.

## Quick start

```fortran
module my_problem
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
contains
  !> Two objectives, both MINIMIZED.
  subroutine objectives(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = 4 * x(1)**2 + 4 * x(2)**2
    f(2) = (x(1) - 5)**2 + (x(2) - 5)**2
  end subroutine objectives

  !> Constraints: feasible when every g <= 0.
  subroutine constraints(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    g(1) = (x(1) - 5)**2 + x(2)**2 - 25
    g(2) = 7.7_dp - (x(1) - 8)**2 - (x(2) + 3)**2
  end subroutine constraints
end module my_problem

program demo
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use nsga2
  use my_problem
  implicit none
  type(nsga2_result) :: res
  integer :: i, stat
  character(:), allocatable :: msg

  call nsga2_optimize(objectives, 2, [0._dp, 0._dp], [5._dp, 3._dp], nsga2_options(), res, &
                      constraints=constraints, n_con=2, stat=stat, errmsg=msg)
  if (stat /= 0) error stop msg
  print '(i0, a)', size(res%f, 2), ' Pareto-optimal solutions (x1, x2 -> f1, f2):'
  do i = 1, size(res%f, 2), 10
    print '(2f8.3, a, 2f9.3)', res%x(:, i), '  ->', res%f(:, i)
  end do
end program demo
```

## Encodings

`nsga2_optimize` is generic: the type of the bounds picks the encoding, and your subroutines receive the variables in that type.

```fortran
call nsga2_optimize(objectives, n_obj, lower, upper, opts, res ...)   ! real(real64) bounds -> x(:) real,    result res%x
call nsga2_optimize(objectives, n_obj, lower, upper, opts, res ...)   ! integer bounds      -> k(:) integer, result res%ints
call nsga2_optimize(objectives, n_obj, n_bits, opts, res ...)         ! bit-string length   -> bits(:) logical, result res%bits
```

For example, a bi-objective knapsack (bits) and a gear-train design (integers), both minimized:

```fortran
subroutine knapsack(bits, f)          ! binary_objective_function
  logical, intent(in) :: bits(:)
  real(dp), intent(out) :: f(:)
  f(:2) = -matmul(merge(1, 0, bits), profit)       ! profit(n_items, 2): maximize -> minimize the negative
end subroutine knapsack

subroutine gear(k, f)                 ! integer_objective_function
  integer, intent(in) :: k(:)
  real(dp), intent(out) :: f(:)
  f(1) = (1 / 6.931_dp - real(k(1) * k(2), dp) / (k(3) * k(4)))**2
  f(2) = maxval(k)
end subroutine gear

call nsga2_optimize(knapsack, 2, 20, nsga2_options(), res, constraints=capacity, n_con=2)   ! res%bits(20, :)
call nsga2_optimize(gear, 2, [12, 12, 12, 12], [60, 60, 60, 60], nsga2_options(), res)       ! res%ints(4, :)
```

Which to use:

- **Integer** when the variables are counts, sizes or indices.
- **Binary** for yes/no choices: subsets, feature selection, knapsack-type problems.
- **Real** for everything continuous.

A real variable can also be encoded in bits, but then you decode it yourself (for example with `decode` from [genetic_algorithm](https://github.com/ElDeriuos/genetic_algorithm)).

## `nsga2_optimize`

```fortran
call nsga2_optimize(objectives, n_obj, lower, upper, opts, res [, constraints, n_con] [, stat, errmsg])   ! real or integer
call nsga2_optimize(objectives, n_obj, n_bits, opts, res [, constraints, n_con] [, stat, errmsg])         ! binary
```

| Argument | Type | Meaning |
|---|---|---|
| `objectives` | subroutine `(x, f)` | Fills `f(1:n_obj)` for the decision vector, which is `real(real64) x(:)`, `integer k(:)` or `logical bits(:)` depending on the encoding. **All objectives are minimized**; to maximize, return the negative. Called `population_size × (max_generations + 1)` times. |
| `n_obj` | `integer`, in | Number of objectives (≥ 1; NSGA-II is meant for 2–3). |
| `lower`, `upper` | `real(real64)(:)` or `integer(:)`, in | Bounds of each decision variable (inclusive); same size, `lower <= upper`. The number of variables is `size(lower)`. Real bounds: write `5._dp`, not `5.0`. |
| `n_bits` | `integer`, in | Binary encoding only: length of the bit string (≥ 1). |
| `opts` | `type(nsga2_options)`, in | Settings, required. `nsga2_options()` gives the defaults; set only what you need: `nsga2_options(population_size=200)`. |
| `res` | `type(nsga2_result)`, out | The non-dominated front of the final population. Reset on every call. |
| `constraints` | subroutine `(x, g)`, optional | Fills `g(1:n_con)`; the solution is feasible when every `g <= 0`. Its first argument has the same type as in `objectives`. Give it together with `n_con`. |
| `n_con` | `integer`, in, optional | Number of constraints (≥ 1). |
| `stat` | `integer`, out, optional | `0` = success, `1` = invalid input (nothing was run). **If omitted, invalid input stops the program** with `nsga2: <message>`. |
| `errmsg` | `character(:), allocatable`, out, optional | Set only when `stat /= 0`: which input is wrong. |

Both user functions are subroutines with assumed-shape arrays. Their interfaces are in the module:
- `objective_function` / `constraint_function` (real),
- `integer_objective_function` / `integer_constraint_function`,
- `binary_objective_function` / `binary_constraint_function`.

 Write them as **module procedures**: internal procedures also work, but gfortran then needs an executable stack. For extra data (coefficients, measurements), use module variables set before the call.

### `nsga2_options`

| Field | Type | Default | Meaning | Valid |
|---|---|---|---|---|
| `population_size` | integer | `100` | Solutions per generation; also the maximum size of the returned front. | ≥ 2 |
| `max_generations` | integer | `250` | Generations to run. `0` returns the non-dominated part of the random initial population. | ≥ 0 |
| `crossover_rate` | real64 | `0.9` | Probability that a parent pair is recombined; otherwise the children are copies. Real and integer: within a crossed pair each variable is crossed by SBX with probability 0.5, and the two children's values are swapped with probability 0.5 (as in Deb's code). Binary: two-point crossover. | [0, 1] |
| `eta_c` | real64 | `20` | SBX distribution index: larger keeps children closer to their parents (pymoo uses 15). Not used by the binary encoding. | ≥ 0 |
| `mutation_rate` | real64 | `-1` | Probability that each variable (or bit) mutates; negative means `1 / number of variables` (or bits). | ≤ 1 |
| `eta_m` | real64 | `20` | Polynomial mutation distribution index: larger gives smaller steps. Not used by the binary encoding. | ≥ 0 |

### `nsga2_result`

| Field | Type | Meaning |
|---|---|---|
| `x(:,:)` | real64 `(n_var, n_front)` | Real encoding: decision vectors of the non-dominated (rank-1) solutions of the final population, one per column. |
| `ints(:,:)` | integer `(n_var, n_front)` | Integer encoding: the same, as integers. |
| `bits(:,:)` | logical `(n_bits, n_front)` | Binary encoding: the same, as bit strings. |
| `f(:,:)` | real64 `(n_obj, n_front)` | Their objective values. |
| `cv(:)` | real64 `(n_front)` | Their constraint violation `sum(max(0, g))`; 0 = feasible. If no feasible solution was found, the front holds the least-violating ones, so check `cv`. |
| `generations` | integer | Generations run. |

Only the field of the encoding you used is allocated. `n_front = size(res%f, 2)` is at most `population_size`; duplicates are possible, and on discrete problems common.

### Reproducible runs

The library never seeds the random number generator; call `random_init(repeatable=.true., image_distinct=.true.)` or `random_seed(put=...)` yourself before `nsga2_optimize`.

### Building blocks

The operators are public. Each takes its random numbers as arguments, so it is deterministic and testable.

| Procedure | What it does |
|---|---|
| `dominates(fa, fb [, cva, cvb])` | Pareto dominance for minimization; with violations, Deb's rule (feasible beats infeasible, smaller violation beats larger). |
| `nondominated_rank(f [, cv])` | Front number of every column of `f(n_obj, n)`: 1 = non-dominated, 2 = dominated only by front 1, … (O(M·N²)). |
| `crowding_distance(f)` | Crowding distance of the points of one front; boundary points get `huge()`. |
| `crowded_tournament(rank, crowd, contestants)` | Winner of each column of `contestants(2, m)`: lower rank, then larger crowding distance. |
| `sbx_crossover(p1, p2, mask, u, lower, upper, eta, c1, c2 [, swap])` | Bounded SBX on genes where `mask`, one draw `u` in [0,1) per gene; `swap` exchanges the children's values. |
| `polynomial_mutation(x, mask, r, lower, upper, eta)` | Bounded polynomial mutation where `mask`; `r` in [0,1): 0 → lower bound, 0.5 → no change, → 1 upper bound. |
| `two_point_crossover(p1, p2, cut1, cut2, c1, c2)` | The children exchange genes `cut1+1 .. cut2` (`0 <= cut1 <= cut2 <= size(p1)`). |
| `igd(front, reference)` | Inverted generational distance of `front(n_obj, n)` to `reference(n_obj, m)`; lower is better. |
| `hypervolume_2d(front, ref_point)` | Area dominated by a 2-objective front up to `ref_point`; higher is better. |

## Benchmark results

Settings and calibration:
- **Settings:** the defaults (population 100, 250 generations); DTLZ1, DTLZ3 and DTLZ6 get 1000 generations.
- **Calibration:** 5 seeds per problem.
- **Reference fronts:** IGD is measured against the true Pareto front. It is computed without the library: analytic where known, otherwise the non-dominated set of a dense grid over the decision space (brute force), or for the discrete problems the enumeration of every solution.
- **Pass rule:** the tests require every returned solution to be feasible and IGD below the limit shown, about 2× the worst seed. For discrete fronts the limit is a number of exact front points found instead.

| Problem | Variables | Objectives | Constraints | Reference front | IGD (best – worst of 5) | Test limit |
|---|---|---|---|---|---|---|
| SCH | 1 | 2 | – | analytic | 0.018 – 0.020 | 0.04 |
| FON | 3 | 2 | – | analytic | 0.0055 – 0.0060 | 0.012 |
| POL | 2 | 2 | – | grid 500² | 0.066 – 0.070 | 0.15 |
| KUR | 3 | 2 | – | grid 100³ | 0.10 – 0.11 | 0.25 |
| VNT | 2 | 3 | – | grid 100² | 0.041 – 0.045 | 0.09 |
| ZDT1 | 30 | 2 | – | analytic | 0.0048 – 0.0051 | 0.011 |
| ZDT2 | 30 | 2 | – | analytic | 0.0048 – 0.0051 | 0.011 |
| ZDT3 | 30 | 2 | – | analytic (disconnected) | 0.0054 – 0.035 | 0.07 |
| ZDT4 | 10 | 2 | – | analytic | 0.0060 – 0.0091 | 0.02 |
| ZDT6 | 10 | 2 | – | analytic | 0.0056 – 0.0081 | 0.02 |
| DTLZ1 | 7 | 3 | – | analytic (plane) | 0.027 – 0.030 | 0.06 |
| DTLZ2 | 12 | 3 | – | analytic (sphere) | 0.064 – 0.070 | 0.14 |
| DTLZ3 | 12 | 3 | – | analytic (sphere) | 0.063 – 0.071 | 0.15 |
| DTLZ4 | 12 | 3 | – | analytic (sphere) | 0.065 – 0.073 | 0.15 |
| DTLZ5 | 12 | 3 | – | analytic (curve) | 0.0054 – 0.0058 | 0.012 |
| DTLZ6 | 12 | 3 | – | analytic (curve) | 0.025 – 0.10 | 0.21 |
| DTLZ7 | 22 | 3 | – | analytic (disconnected) | 0.073 – 0.082 | 0.17 |
| BNH | 2 | 2 | 2 | grid 500² | 0.43 – 0.50 | 1.1 |
| SRN | 2 | 2 | 2 | grid 500² | 1.03 – 1.21 | 2.5 |
| TNK | 2 | 2 | 2 | grid 500² | 0.0042 – 0.0050 | 0.011 |
| CONSTR | 2 | 2 | 2 | grid 500² | 0.020 – 0.022 | 0.045 |
| OSY | 6 | 2 | 6 | – | hypervolume 1.30–1.34× random search with 8× the evaluations | > 1.2× |
| OneMinMax | 20 bits | 2 | – | analytic (21 points) | 21 of 21 points found, every seed | all 21 |
| LOTZ | 20 bits | 2 | – | analytic (21 points) | 21 of 21 points found, every seed | all 21 |
| ZDT5 | 80 bits | 2 | – | analytic (31 points) | 0.26 – 0.59 | 1.2 |
| Knapsack | 20 bits | 2 | 2 | enumeration, 2²⁰ subsets (13 points) | 0 – 4.9 (12 or 13 of 13 points found) | 10 |
| Gear train | 4 integers | 2 | – | enumeration, 49⁴ designs (28 points) | 23 – 25 of 28 points found | ≥ 18 |

Notes on reading the table:
- **Units:** IGD is in the units of each problem's objectives (SRN's objectives span about 220, for example).
- **3-objective DTLZ:** 100 points cannot cover a surface closely, so an IGD of about 0.06–0.07 is the floor for a fully converged run.
- **ZDT3:** the worst seed misses the front's last, narrow segment (f1 ≈ 0.82–0.85); that is the seed the test and its figure use.
- **ZDT5:** it is deceptive by design. Each 5-bit block pulls towards all zeros, but its optimum is all ones. Runs stall at g = 12–16 (the true front has g = 10), and 1000 generations change nothing; a random population scores IGD 3.7–5.2.

Where the benchmarks come from:
- **Real-coded:** Deb et al. 2002 (SCH, FON, POL, KUR, VNT, BNH, SRN, TNK, CONSTR, OSY), Zitzler, Deb & Thiele 2000 (ZDT) and Deb, Thiele, Laumanns & Zitzler 2002/2005 (DTLZ, 3 objectives). They were checked against the pymoo source (ZDT, DTLZ, KUR, BNH, SRN, TNK, OSY) and Wikipedia's *Test functions for optimization* (SCH, FON, POL, VNT, CONSTR).
- **Binary:**
  - ZDT5 was checked against pymoo; it is used unnormalized here.
  - OneMinMax and LOTZ are the bit-string problems of runtime analysis (Giel & Lehre; Laumanns, Thiele & Zitzler).
  - The knapsack follows Zitzler & Thiele 1999: profits and weights are random integers in [10, 100], and each capacity is half the knapsack's total weight. This instance has 20 items, so its exact front can be enumerated; the data is in `test/benchmarks.f90`.
- **Integer:** the gear train is from Deb, Pratap & Moitra 2000, in the form given in arXiv:2012.00115, Appendix A. Its exact front ends at the known optimum, error 2.7·10⁻¹² with largest gear 49.

## Pareto fronts

Each figure shows a test's own run, with the same seed and settings as the test. The grey points are the reference front, the blue ones what NSGA-II returns. The caption gives the measured value next to the test's limit.

**Classic and ZDT** (real, 2 objectives)

| | |
|---|---|
| ![SCH](docs/images/front_sch.png) | ![FON](docs/images/front_fon.png) |
| ![POL](docs/images/front_pol.png) | ![KUR](docs/images/front_kur.png) |
| ![ZDT1](docs/images/front_zdt1.png) | ![ZDT2](docs/images/front_zdt2.png) |
| ![ZDT3](docs/images/front_zdt3.png) | ![ZDT4](docs/images/front_zdt4.png) |
| ![ZDT6](docs/images/front_zdt6.png) | |

**3 objectives** (real)

| | |
|---|---|
| ![VNT](docs/images/front_vnt.png) | ![DTLZ1](docs/images/front_dtlz1.png) |
| ![DTLZ2](docs/images/front_dtlz2.png) | ![DTLZ3](docs/images/front_dtlz3.png) |
| ![DTLZ4](docs/images/front_dtlz4.png) | ![DTLZ5](docs/images/front_dtlz5.png) |
| ![DTLZ6](docs/images/front_dtlz6.png) | ![DTLZ7](docs/images/front_dtlz7.png) |

**Constrained** (real)

| | |
|---|---|
| ![BNH](docs/images/front_bnh.png) | ![SRN](docs/images/front_srn.png) |
| ![TNK](docs/images/front_tnk.png) | ![CONSTR](docs/images/front_constr.png) |
| ![OSY](docs/images/front_osy.png) | |

**Binary and integer**

| | |
|---|---|
| ![OneMinMax](docs/images/front_oneminmax.png) | ![LOTZ](docs/images/front_lotz.png) |
| ![ZDT5](docs/images/front_zdt5.png) | ![Knapsack](docs/images/front_knapsack.png) |
| ![Gear train](docs/images/front_gear.png) | |

## Tests

```sh
fpm test                                               # operators, benchmarks, invalid input
cmake -S . -B build/cmake && cmake --build build/cmake
ctest --test-dir build/cmake --output-on-failure       # the same, plus expected-failure tests
```

Operator tests compare against worked examples computed independently from the published formulas. Expected-failure tests (in `test/xfail`, CMake only, since fpm cannot mark a test as expected to fail): invalid input without `stat` must stop the program, and a negative control proves the benchmark checks detect an NSGA-II that does not evolve.

## References

- K. Deb, A. Pratap, S. Agarwal, T. Meyarivan, "A fast and elitist multiobjective genetic algorithm: NSGA-II", *IEEE Trans. Evol. Comput.* 6(2), 2002.
- K. Deb, R. B. Agrawal, "Simulated binary crossover for continuous search space", *Complex Systems* 9, 1995.
- E. Zitzler, K. Deb, L. Thiele, "Comparison of multiobjective evolutionary algorithms: empirical results", *Evol. Comput.* 8(2), 2000.
- K. Deb, L. Thiele, M. Laumanns, E. Zitzler, "Scalable test problems for evolutionary multiobjective optimization", 2005 (DTLZ).
- E. Zitzler, L. Thiele, "Multiobjective evolutionary algorithms: a comparative case study and the strength Pareto approach", *IEEE Trans. Evol. Comput.* 3(4), 1999 (0/1 knapsack).
- M. Laumanns, L. Thiele, E. Zitzler, "Running time analysis of multiobjective evolutionary algorithms on pseudo-Boolean functions", *IEEE Trans. Evol. Comput.* 8(2), 2004 (LOTZ).
- O. Giel, P. K. Lehre, "On the effect of populations in evolutionary multi-objective optimisation", *Evol. Comput.* 18(3), 2010 (OneMinMax).
- K. Deb, A. Pratap, S. Moitra, "Mechanical component design for multiple objectives using elitist non-dominated sorting GA", PPSN VI, LNCS 1917, 2000 (gear train).
- J. Blank, K. Deb, "pymoo: Multi-objective optimization in Python", *IEEE Access* 8, 2020 (integer variables by rounding).

## License

[MIT](LICENSE): free to use, modify and redistribute, including commercially; keep the copyright notice.
