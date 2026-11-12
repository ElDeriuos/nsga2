# 3. Conventions shared with genetic_algorithm

Date: 2026-11-12 · Status: accepted

Same choices as the `genetic_algorithm` package (its ADRs 3-6), for the same reasons:

- **The caller owns the RNG:** the library uses `random_number` and never seeds it.
- **Public operators take their random draws as arguments,** so they are deterministic and tested against literal worked examples.
- **Invalid input:** optional `stat`/`errmsg`; without `stat`, `error stop 'nsga2: <message>'`.
- **Packaging:** `fpm.toml` and `CMakeLists.txt`; tests are plain programs that `error stop` on failure; expected-failure programs in `test/xfail` run under ctest only.
- **No trampolines:** user functions reach the library as procedure arguments; the optional constraint function is passed on to a module procedure rather than called from an internal one (gfortran 16 flags that call as an implicit interface).
