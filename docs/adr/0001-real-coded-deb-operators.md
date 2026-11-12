# 1. Real-coded NSGA-II with Deb's bounded SBX and polynomial mutation

Date: 2026-11-12 · Status: accepted

**Decision.** Real-valued decision variables only. Offspring by crowded binary tournament, bounded SBX (crossover_rate 0.9; each variable crossed with probability 0.5; the two children's values swapped with probability 0.5) and bounded polynomial mutation (rate 1/n). Distribution indices `eta_c = eta_m = 20`, as in Deb et al. 2002 (pymoo uses `eta_c = 15`). Survival: best N of parents + children by rank, then crowding distance.

**Why the swap matters.** Without it each child is a perturbed copy of one parent. Measured on ZDT1 (seed 1): IGD 0.137 after 250 generations without the swap, 0.005 with it.

**Rejected for now.** Binary/integer encodings (the classic NSGA-II paper also had binary; no request here), NSGA-III for many objectives, parallel evaluation.
