# Known vignette-render failure patterns

These are bugs that have shipped repeatedly because the per-paper
render gate (Phase 5 step 2) was skipped, fabricated, or run with
events too simple to exercise the failure. The runner-merge-claude-
branches consolidation runs every vignette in parallel as a HARD
gate, so anything you ship broken at the paper level will be caught
at merge time and will block the consolidation PR. The cost of
fixing it then is much higher (you're un-stacking from a 130-branch
merge, not from your local one-paper worktree).

Each pattern below shows the error you'll see, the root cause, and
the canonical fix. Whenever you write a model or vignette, scan
this list and double-check you're not stepping into one.

## 1. `chol(): decomposition failed` during `rxSolve`

**Cause.** The model's OMEGA (IIV covariance) matrix has a
non-positive eigenvalue — usually a multi-eta block where the
published paper reports "perfect positive (or negative) correlation"
between several etas. Encoding that literally as a 3x3 matrix with
off-diagonals = `sqrt(var_i * var_j)` gives a rank-1 matrix
(determinant 0). rxode2's Cholesky-based sampler cannot decompose
a singular matrix. Tiny rounding noise in the off-diagonals can also
flip the matrix to mildly indefinite.

**Wrong** (singular OMEGA — `chol()` fails):
```r
ini({
  # ...
  etalvc + etalq2 + etalvp2 ~ c(
    0.014297,
    0.044990, 0.141562,   # sqrt(0.014297 * 0.141562) — perfect +1 correlation
    0.029444, 0.092639, 0.060625
  )
})
```
Every pair has correlation exactly 1, so the block is rank-1: its determinant is 0 (and rounding can tip it mildly negative), and the Cholesky sampler cannot decompose it.

**Right** — keep the IIV in `ini()` and nudge the matrix positive definite. Scale the **off-diagonals** by 0.99 (correlation 1.00 -> 0.99); the diagonal variances — the published CVs — are kept exactly, so this is the smallest change that makes OMEGA positive definite. **Do not move the IIV into `model()`.**
```r
ini({
  # ...
  # Paper reports "perfect (+1) correlation" across the V2/Q4/V4 block, which
  # is exactly singular. Nudge each correlation to 0.99 by scaling the
  # off-diagonals by 0.99; variances (the published CVs) are unchanged.
  # (If chol() still fails on a given matrix, drop the factor a little more,
  # e.g. 0.95.)
  etalvc + etalq2 + etalvp2 ~ c(
    0.014297,
    0.99 * 0.044990, 0.141562,
    0.99 * 0.029444, 0.99 * 0.092639, 0.060625
  )
})
model({
  # model() is UNCHANGED — structure only, no eta bookkeeping:
  vc  <- exp(lvc  + etalvc)  * allom_v  * cov_factor
  q2  <- exp(lq2  + etalq2)  * allom_cl * cov_factor
  vp2 <- exp(lvp2 + etalvp2) * allom_v  * cov_factor
})
```
Inline arithmetic in the `c(...)` is valid — rxode2 evaluates it, so the `0.99 *` nudge stays visible in the source trace. The nudge perturbs only the idealized "perfect" correlation, never the reported variances, and leaves `model()` for ODE structure and algebra.

**Do not** instead introduce a shared standardized eta (`eta_v2q4v4 ~ 1.0`) and scale it per parameter inside `model()` (`etalvc <- 0.119571 * eta_v2q4v4`, ...). That is mathematically equivalent but pushes covariance bookkeeping into the model body, which is reserved for structure. `inst/modeldb/specificDrugs/Fanta_2007_ciclosporin.R` currently uses that older `model()`-scaling form and is a candidate to migrate to the `ini()`-nudge above.

## 2. `'cmt' on observation record ... undefined compartment` OR `The following parameter(s) are required for solving: <state>`

**Cause.** The model declares ODE states (`d/dt(central) <- ...`)
plus algebraic observables (`Cc <- central / vc`). The VIGNETTE
event table writes observations with `cmt = "Cc"` (the observable
name, not the ODE state name). When rxode2 / rxUi processes this,
it auto-injects `cmt(Cc)` AFTER the ODE state slots, renumbering
the slot indices. References to ODE states past the inserted slot
(e.g. `peripheral1`, `csf`, `fetus`, `central_mhd`) become
unresolvable and rxode2 reports them as "required parameters"
because their slot is gone.

The error usually surfaces as one of these:
- `'cmt' on observation record or on a undefined compartment`
- `The following parameter(s) are required for solving: <name>`
- `tad(depot)` returning NA on every observation row (downstream:
  PKNCA gets a zero-row filter)

**Wrong (the symptom — bug is in the event table):**
```r
# In the vignette:
obs <- subj |>
  tidyr::crossing(time = seq(0, 24, by = 1)) |>
  mutate(amt = NA_real_, evid = 0,
         cmt = "Cc",                  # <-- bug: observable name, not a compartment
         ...)
```

**Wrong (the "fix" that pollutes the model body — REJECTED):**
```r
# In the model file — DO NOT DO THIS:
model({
  cmt(depot)
  cmt(central)
  cmt(peripheral1)
  cmt(Cc)        # <-- pollutes model body to silence the symptom
  d/dt(depot) <- -ka * depot
  ...
})
```

**Right (fix the event table to use the actual ODE state name):**
```r
# In the vignette:
obs <- subj |>
  tidyr::crossing(time = seq(0, 24, by = 1)) |>
  mutate(amt = NA_real_, evid = 0,
         cmt = "central",             # <-- the ODE state, not the observable
         ...)
```

rxode2 always returns every algebraic observable (like `Cc`) as a
column in the output dataframe regardless of which compartment the
observation row pointed at. The `cmt` on an observation row tells
rxode2 *when* (which compartment's slot to align the timing on),
not *what* to report. The observable lookup is automatic.

For multi-output models where you have several observables backed
by different ODE states (e.g. `Cc <- central/vc`, `Ccsf <- csf/vcsf`),
use the corresponding ODE state name per observation row, or use
`dvid()` in the model body to declare a DV id mapping — `dvid()` is
semantic (declares the endpoint structure) and is fine in `model()`;
`cmt()` declarations are not.

**Process discipline.** When you write the vignette event table, the
`cmt` values for observation rows must be one of the model's declared
ODE state names (i.e. names appearing in a `d/dt(<name>) <- ...`
line). Search the model file for `d/dt(` and use one of those names.
Never use the name of an algebraic intermediate (anything that's just
`X <- ...` in `model({})`, including the observation symbol like
`Cc`). This is the single most common consolidation-blocker.

## 3. `summarise: unique(x) returned >1 value` in dplyr

**Cause.** A typical-value summary inherits per-id variability
from a preceding non-typical simulation, so a column that you
expected to have one value per group actually varies.

**Wrong:**
```r
sim_typical <- sim |> filter(...)
summary <- sim_typical |>
  group_by(cohort, AGE, CRCL) |>
  summarise(CL_model = unique(cl))   # blows up if cl varies within (cohort, AGE, CRCL)
```

**Right** — either tighten the grouping to include the covariate
that varies, or aggregate explicitly:
```r
summary <- sim_typical |>
  group_by(cohort, AGE, CRCL) |>
  summarise(CL_model = mean(cl))     # averages within group
# or:
summary <- sim_typical |>
  group_by(id, cohort, AGE, CRCL) |>  # finer grouping
  summarise(CL_model = first(cl))
```

## 4. `PKNCAconc: data must have at least one row`

**Cause.** A `filter()` chain upstream dropped every row. Usually
caused by referencing an observable that ended up NA at every
observation (often a follow-on from pattern #2: `tad(depot)` is NA,
so `Cc` is NA, so a `filter(Cc > 0)` drops everything).

**Fix.** Trace the filter chain back to the first stage where
`nrow(df)` becomes 0. Often the upstream fix is pattern #2.

## 5. `callr timed out` (parallel build only)

**Cause.** The vignette runs fine in isolation (e.g. 4 min) but
exceeds the per-vignette ceiling (default 900s) when 32 callr
workers contend for CPU. Symptoms: passes in your local single-
vignette render, fails in the consolidation merge's parallel
validator.

**Fix.** Reduce simulation cohort sizes (`n_per_group` from 100
to 50, etc.). A vignette should comfortably complete in under
~5 minutes single-threaded; if it doesn't, the cohort is too
large for the validation use case (which is illustrative
demonstration, not production VPC).

## 5b. `'dvid'->'cmt' on observation record` or `required parameter: <ode_state>` AFTER you've already used named compartments and dvid

If you've already followed pattern 2 (named ODE-state `cmt =` in event tables + `dvid = 1L` on observation rows for multi-output models) and rxSolve STILL errors with the dvid mapping complaint or a "required parameter" for one of the ODE states, you've hit the second rxode2 bug catalogued in `reports/rxode2-tad-state-arg-bug-issue.md`: `rxSolve.rxUi`'s default `useLinCmt = TRUE` performs an automatic ODE→linCmt conversion that corrupts the dvid→cmt mapping for many multi-output / multi-state models.

**Workaround.** Pass `useLinCmt = FALSE` to every `rxode2::rxSolve()` call in the vignette:

```r
sim <- rxode2::rxSolve(
  mod, events = events,
  keep = c("WT", "treatment"),
  useLinCmt = FALSE   # rxode2's ODE->linCmt auto-conversion breaks dvid mapping for this model
) |> as.data.frame()
```

Affects models like Germovsek 2018 meropenem (Cc + Ccsf on [central, csf]), Jelliffe 2014 digoxin (Cc on [depot, central, peripheral1]), Ngamprasertwong 2016 propofol sheep (Cc + Cfetus on [central, peripheral1, fetus]), Rodrigues 2017 oxcarbazepine (Cc + Cc_mhd on [depot, central, peripheral1, central_mhd]), Themans 2019 meropenem (Cc + Celf on [central, peripheral1, peripheral2]). The model body itself is fine — only the simulation call needs the flag.

## 6. `something went wrong in compilation`

**Cause.** rxode2's C compiler failed on the model. Usually a
syntax error in the `model({})` block, an undefined variable, or
a state reference that doesn't match a declared `cmt()` /
`d/dt()`.

**Fix.** Load the model directly with verbose compilation and
read the actual compiler error:

```r
m <- readModelDb("YourModel")
m_built <- m()  # this is where compilation happens
```

The error message at this level is much more specific than the
one in the vignette render.

## Process reminder

The runner-merge-claude-branches skill runs `verify_vignettes_parallel.R`
as a HARD gate before pushing the consolidation branch. Patterns
1-6 above are caught by that gate, but the cost of fixing them at
merge time is high (un-stacking 14 fixes from a 130-branch worktree
takes hours). **Fix them at the paper level, before you push the
per-paper branch.**

When you run Phase 5 step 2 (the per-paper RENDER_GATE), make sure
your vignette events table actually produces each algebraic
observable — at minimum, observation rows **on the relevant ODE
state** (`cmt = "central"`, etc.). rxode2 returns the observable
(e.g. `Cc`) as a column at those rows, so this exercises the
observable-computation path WITHOUT the slot-renumbering bug. Do
NOT write `cmt = "<observable-name>"` (`cmt = "Cc"`) to "exercise"
the pattern — that is the bug itself (pattern above). A render that
only doses (no observations) or that only observes with
`cmt = NA_character_` will pass the gate and still ship a broken
vignette.
