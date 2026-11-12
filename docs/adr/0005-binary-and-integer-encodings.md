# 5. Binary and integer encodings

Date: 2026-11-12 · Status: accepted (supersedes the "rejected for now" note on encodings in ADR 1)

**Decision.**
- **Interface:** `nsga2_optimize` becomes a generic with three specifics, mirroring `genetic_algorithm`.
  - Real bounds give `objective_function(x, f)` and the result in `res%x`.
  - Integer bounds give `integer_objective_function(k, f)` and the result in `res%ints`.
  - `n_bits` gives `binary_objective_function(bits, f)` with `logical` bits and the result in `res%bits`.
  - Constraints are available in every encoding, typed the same way.
- **One engine:** a single loop serves all three. Genes are stored as `real(dp)`: bits as 0/1, integers as whole numbers.
- **Binary operators:** two-point crossover over the whole string and bit-flip mutation (rate 1/n_bits). This is Deb's NSGA-II C code (`bincross`, `bin_mutate_ind`) applied to one bit string instead of per variable.
- **Integer operators:** SBX and polynomial mutation on the integers as reals, then rounding. This is pymoo's standard setup (`SBX` and `PM` with `RoundingRepair`).

**Evidence.**
- **Integer:** on the gear train (4 integers in [12, 60], exact front of 28 points by enumeration), the defaults η_c = η_m = 20 found 23–25 points over 5 seeds; pymoo's integer example uses η = 3, which found 22–24. So the real-coded defaults are kept for integers.
- **Binary:** LOTZ and OneMinMax (20 bits, population 100 ≥ 4(n+1), the size Zheng, Liu & Doerr 2022 analyse) cover all 21 front points on every seed. On the knapsack, 12–13 of the 13 exact front points are found.

**Rejected.**
- **A single `real(dp)` callback for all encodings:** less code, but callers would `nint()` or compare bits themselves, unlike `genetic_algorithm`.
- **Creep mutation for integers** (as in `genetic_algorithm`): SBX/PM with rounding already reaches the exact gear front, so there is no measured need.
- **A `decode` for binary-coded reals:** the real encoding covers continuous variables.
- **Duplicate elimination (pymoo's `eliminate_duplicates`):** not needed by any test front.
