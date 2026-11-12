# nsga2

NSGA-II, the elitist multi-objective genetic algorithm of Deb et al. (2002), as one Modern Fortran module (F2018, `real64`, no dependencies), usable from any program via `use nsga2`.

- Real-coded: simulated binary crossover (SBX) and polynomial mutation, both bounded
- Fast non-dominated sorting, crowding distance, crowded binary tournament, elitist (μ+λ) survival
- Constraints `g(x) <= 0` with Deb's constrained dominance
- Quality indicators: IGD and 2-D hypervolume
- Tested on 22 standard benchmarks: SCH, FON, POL, KUR, VNT, ZDT1–4, ZDT6, DTLZ1–7, BNH, SRN, TNK, CONSTR, OSY

Terms are defined in [GLOSSARY.md](GLOSSARY.md); design decisions in [docs/adr](docs/adr). A single-objective companion package: [genetic_algorithm](https://github.com/ElDeriuos/genetic_algorithm).

## Install

**fpm:** add to your project's `fpm.toml`, then `fpm run` / `fpm build`:

```toml
[dependencies]
nsga2 = { git = "https://github.com/ElDeriuos/nsga2", tag = "v0.1.0" }
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

## `nsga2_optimize`

```fortran
call nsga2_optimize(objectives, n_obj, lower, upper, opts, res [, constraints, n_con] [, stat, errmsg])
```

| Argument | Type | Meaning |
|---|---|---|
| `objectives` | subroutine `(x, f)` | Fills `f(1:n_obj)` for the decision vector `x`. **All objectives are minimized**; to maximize, return the negative. Called `population_size × (max_generations + 1)` times. |
| `n_obj` | `integer`, in | Number of objectives (≥ 1; NSGA-II is meant for 2–3). |
| `lower`, `upper` | `real(real64)(:)`, in | Bounds of each decision variable; same size, `lower <= upper`. The number of variables is `size(lower)`. Write `5._dp`, not `5.0`. |
| `opts` | `type(nsga2_options)`, in | Settings, required. `nsga2_options()` gives the defaults; set only what you need: `nsga2_options(population_size=200)`. |
| `res` | `type(nsga2_result)`, out | The non-dominated front of the final population. Reset on every call. |
| `constraints` | subroutine `(x, g)`, optional | Fills `g(1:n_con)`; `x` is feasible when every `g <= 0`. Give it together with `n_con`. |
| `n_con` | `integer`, in, optional | Number of constraints (≥ 1). |
| `stat` | `integer`, out, optional | `0` = success, `1` = invalid input (nothing was run). **If omitted, invalid input stops the program** with `nsga2: <message>`. |
| `errmsg` | `character(:), allocatable`, out, optional | Set only when `stat /= 0`: which input is wrong. |

Both user functions are subroutines with assumed-shape `real(real64)` arrays (interfaces `objective_function` and `constraint_function` in the module). Write them as **module procedures**: internal procedures also work, but gfortran then needs an executable stack. For extra data (coefficients, measurements), use module variables set before the call.

### `nsga2_options`

| Field | Type | Default | Meaning | Valid |
|---|---|---|---|---|
| `population_size` | integer | `100` | Solutions per generation; also the maximum size of the returned front. | ≥ 2 |
| `max_generations` | integer | `250` | Generations to run. `0` returns the non-dominated part of the random initial population. | ≥ 0 |
| `crossover_rate` | real64 | `0.9` | Probability that a parent pair is recombined by SBX; otherwise the children are copies. Within a crossed pair each variable is crossed with probability 0.5 and the two children's values are swapped with probability 0.5 (as in Deb's code). | [0, 1] |
| `eta_c` | real64 | `20` | SBX distribution index: larger keeps children closer to their parents (pymoo uses 15). | ≥ 0 |
| `mutation_rate` | real64 | `-1` | Probability that each variable mutates; negative means `1 / number of variables`. | ≤ 1 |
| `eta_m` | real64 | `20` | Polynomial mutation distribution index: larger gives smaller steps. | ≥ 0 |

### `nsga2_result`

| Field | Type | Meaning |
|---|---|---|
| `x(:,:)` | real64 `(n_var, n_front)` | Decision vectors of the non-dominated (rank-1) solutions of the final population, one per column. |
| `f(:,:)` | real64 `(n_obj, n_front)` | Their objective values. |
| `cv(:)` | real64 `(n_front)` | Their constraint violation `sum(max(0, g))`; 0 = feasible. If no feasible solution was found, the front holds the least-violating ones, so check `cv`. |
| `generations` | integer | Generations run. |

`n_front = size(res%f, 2)` is at most `population_size`; duplicates are possible.

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
| `igd(front, reference)` | Inverted generational distance of `front(n_obj, n)` to `reference(n_obj, m)`; lower is better. |
| `hypervolume_2d(front, ref_point)` | Area dominated by a 2-objective front up to `ref_point`; higher is better. |

## Benchmark results

Default settings (population 100, 250 generations; DTLZ1, DTLZ3 and DTLZ6 1000 generations), 5 seeds. IGD is measured against the true Pareto front: analytic where known, otherwise the non-dominated set of a dense grid over the decision space (brute force, independent of the library). The tests require IGD below the limit shown (about 2× the worst seed) and every returned solution feasible.

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

IGD is in the units of each problem's objectives (SRN's objectives span about 220, for example). For the 3-objective DTLZ problems, 100 points cannot cover a surface closely, so an IGD of about 0.06–0.07 is the floor for a fully converged run.

Benchmark definitions follow Deb et al. 2002 (SCH, FON, POL, KUR, VNT, BNH, SRN, TNK, CONSTR, OSY), Zitzler, Deb & Thiele 2000 (ZDT) and Deb, Thiele, Laumanns & Zitzler 2002/2005 (DTLZ, 3 objectives). They were checked against the pymoo source (ZDT, DTLZ, KUR, BNH, SRN, TNK, OSY) and Wikipedia's *Test functions for optimization* (SCH, FON, POL, VNT, CONSTR).

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

## License

[MIT](LICENSE): free to use, modify and redistribute, including commercially; keep the copyright notice.
