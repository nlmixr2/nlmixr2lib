# Validation patterns for endogenous and mechanistic models

PKNCA-based NCA is the wrong validation target for endogenous turnover, steady-state-balance, and enzyme-kinetic models. There is no dose and no absorption-distribution-elimination profile to integrate. Instead, the validations that catch translation errors in this model class are:

1. **Steady-state check** — does the model hold at the reported baseline indefinitely?
2. **Perturbation-recovery** — if the state is displaced, does it return to the reported baseline?
3. **Mass-balance / flux check** — at steady state, do production and elimination fluxes sum to zero?
4. **Dimensional analysis** — do the units on every term in every ODE line up?

Use one per section of the validation vignette. Skip PKNCA entirely unless the paper reports an NCA-style analysis of one of the model outputs.

## 1. Steady-state check

Solve the model without any intervention for several times the slowest relaxation time, with the initial condition set to the paper's reported baseline. The state should stay at that baseline to machine precision.

```r
mod <- readModelDb("igg_kim_2006")
ev  <- et(amt = 0, time = 0) %>% et(seq(0, 90, by = 1))  # 90 days
s   <- rxSolve(mod, ev)
stopifnot(all.equal(range(s$igg), c(12.1, 12.1), tolerance = 1e-6))
```

Failure modes this catches:

- A rate constant with the wrong sign (production vs. loss).
- A missing term in the ODE (recycling, renal sink, transamination).
- A reference value (baseline, `css`) typed wrong in `ini()`.

## 2. Perturbation-recovery

Displace the state to, say, `0.5 * bl` and `2 * bl` and run forward. The trajectory should monotonically approach the reported baseline (assuming a single stable attractor, as most endogenous turnover models have).

```r
ev_low  <- et(amt = 0, time = 0) %>% et(seq(0, 90, by = 1))
s_low   <- rxSolve(mod, ev_low, inits = c(igg = 0.5 * 12.1))
s_high  <- rxSolve(mod, ev_low, inits = c(igg = 2.0 * 12.1))
# Both should end within tolerance of 12.1
```

Failure modes this catches:

- Dynamic rate terms (e.g., `krmr = Jmax / (V1·(Km + igg))`) implemented as the steady-state value rather than concentration-dependent.
- A fictional equilibrium that is not the paper's baseline.

## 3. Mass-balance / flux check

At steady state, compute each flux term (production, elimination, recycling) and confirm they sum to zero within rounding.

For `igg_kim_2006` at `igg = Css = 12.1`:

```
krmr  = Jmax / (V1 · (Km + Css))
      = 147 / (42 · (21 + 12.1)) = 0.1058 /day
kcat  = kint - krmr = 0.18 - 0.1058 = 0.0742 /day
jpro  = kcat · V1 · Css = 0.0742 · 42 · 12.1 = 37.72 mg/kg/day
d/dt(igg) = jpro / V1 - kcat · igg
         = 37.72 / 42 - 0.0742 · 12.1 = 0.8981 - 0.8981 = 0
```

Do this symbolically as well as numerically. Symbolic cancellation confirms the equations are self-consistent; numerical simulation confirms the solver sees the same steady state.

## 4. Dimensional analysis

Mechanistic models routinely mix per-kg volumes, mass concentrations, and fractional rate constants. A missing `WT`, a `mL/kg` label on an `mg/kg` value, or a stray unit conversion silently gives a numerically plausible but biologically wrong result.

For every ODE line:

1. Write down the units of every symbol on the right-hand side.
2. Multiply them out.
3. The result must equal `d/dt(state)` = `[state units] / [time units]`.
4. Do the same for every derived quantity (intermediate rates, flux variables, augmentation outputs like `daily_phe_intake`).

Past real examples:

- `igg_kim_2006`: `V1` labeled `(mg/kg)` in the file but must be `(mL/kg)` for the equations to be dimensionally consistent (value was correct; label was wrong). Caught only by dimensional analysis.
- `phenylalanine_charbonneau_2021`: an augmentation line `daily_phe_intake = 24 * vd * (...) / f_gut_plasma` carried a stray `vd` factor, reporting `(L/kg)·mg/day` instead of `mg/day`. Caught only by dimensional analysis (all parameter values matched the paper; simulation ran without error).

Dimensional analysis is **not optional** for endogenous and mechanistic models. Include an explicit unit-table in the vignette's "Source trace" section showing units for every ODE term.

## Vignette outline for endogenous models

Replace the PKNCA section of the standard template with these sections, in order:

1. Header and setup.
2. Population / biological context (healthy subjects, disease indicators like PKU, age range).
3. Source trace (tables and equations) — **with an explicit units table**.
4. Parameter table (paper vs. file, side by side).
5. Steady-state check (code chunk, expect no drift).
6. Perturbation-recovery (code chunk, figure of trajectories converging to baseline).
7. Mass-balance / flux check (symbolic, confirm all fluxes sum to zero at SS).
8. Dimensional analysis (narrative or table; call out any non-obvious conversions).
9. Figure replication — if the paper shows a trajectory under perturbation, dietary change, or enzyme-activity reduction, reproduce it.
10. Assumptions and deviations.

## Stop-and-ask triggers for endogenous models

- Paper reports only symbolic equations without point values for every rate constant → ask user for values or whether to leave them as TODOs.
- Paper's reference baseline depends on a covariate (e.g., body weight, disease state) and the file needs a default → propose the value and confirm.
- Dimensional analysis reveals the paper's published equation carries extra or missing unit factors → **do not silently "fix" the equation** (the file must reproduce the paper). Flag the discrepancy to the user and document it in a comment. The Charbonneau `v_renal = phe * cl_renal * vd` term is an example: dimensionally mixed but reproduces the paper's published form exactly.
