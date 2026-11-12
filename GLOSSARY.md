# Glossary

Terms as used in the code, tests and docs of `nsga2`.

| Term | Meaning here |
|---|---|
| **Objective** | One of `n_obj` values the user's function returns for a solution. All objectives are **minimized**; to maximize `f`, return `-f`. |
| **Constraint** | One of `n_con` values `g(j)` returned by the optional constraint function. A solution is **feasible** when every `g(j) <= 0`. |
| **Constraint violation (cv)** | `sum(max(0, g))`; 0 for feasible solutions. |
| **Dominates** | `a` dominates `b` when `a` is no worse in every objective and strictly better in at least one. With constraints (Deb's rule): a feasible solution dominates an infeasible one; of two infeasible solutions the one with smaller cv dominates; two feasible ones compare by objectives. |
| **Front / rank** | Non-dominated sorting splits a population into fronts: rank 1 is dominated by nobody, rank 2 only by rank 1, and so on. |
| **Pareto front** | The best possible trade-off curve/surface of a problem. The result returns the rank-1 front of the final population as its approximation. |
| **Crowding distance** | Per front: how much empty space surrounds a solution in objective space (sum over objectives of the normalized gap between its neighbours). Boundary solutions get infinity. Prefers spread-out fronts. |
| **Crowded tournament** | Binary tournament: lower rank wins; on equal rank, larger crowding distance wins. |
| **SBX** | Simulated binary crossover (Deb & Agrawal 1995), bounded version. The distribution index `eta_c` controls how close children stay to their parents (larger = closer). |
| **Polynomial mutation** | Bounded mutation (Deb); `eta_m` controls the step size (larger = smaller steps). |
| **Generation** | Make `N` children (tournament → SBX → mutation), merge with the `N` parents, keep the best `N` by rank then crowding distance (elitist). |
| **Encoding** | How decision variables are represented: real (`real64` in bounds), integer (whole numbers in bounds) or binary (a bit string of `n_bits`). Chosen by the arguments of `nsga2_optimize`; see ADR 5. |
| **Two-point crossover** | Binary encoding: the two children exchange the bits between two cut points. |
| **Bit-flip mutation** | Binary encoding: each bit flips with probability `mutation_rate` (default `1/n_bits`). |
| **Rounding (integer encoding)** | Integer children come from SBX and polynomial mutation on the integers as reals, then are rounded to the nearest whole number. |
| **Discrete front** | A Pareto front of isolated points (binary/integer problems); tests count how many of its exact points are found. |
| **IGD** | Inverted generational distance: mean distance from each point of a reference front to the nearest found solution. Lower is better; 0 = reference fully covered. |
| **Hypervolume (2-D)** | Area dominated by the front and bounded by a reference point. Higher is better. |
| **Expected-failure test** | A test program that must exit with a specific error; ctest checks its output for that message. |
