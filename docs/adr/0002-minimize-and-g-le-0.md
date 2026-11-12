# 2. Objectives are minimized; constraints are g(x) <= 0

Date: 2026-11-12 · Status: accepted

**Decision.** All objectives are minimized (the multi-objective convention of Deb and pymoo); a constraint is satisfied when `g <= 0`. The violation of a solution is `sum(max(0, g))`, and dominance follows Deb's rule: feasible beats infeasible, smaller violation beats larger, two feasible solutions compare by Pareto dominance.

**Note.** This is the opposite sign convention to the `genetic_algorithm` package, which maximizes a single fitness. Each follows the convention of its field.
