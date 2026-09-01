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
m_built <- rxode2::rxode(m)  # this is where compilation happens
```

The error message at this level is much more specific than the
one in the vignette render.

`readModelDb("YourModel")()` also works and is the idiom used
throughout the existing vignettes — but only because each vignette
calls `library(rxode2)`. Calling the model function relies on
`ini()`/`model()` being **attached**; in a script that merely does
`devtools::load_all()` or `library(nlmixr2lib)` it fails with
`could not find function "ini"`. `rxode2::rxode()` has no such
dependency, so prefer it in helper scripts and diagnostics.

## 7. `object of type 'closure' is not subsettable`, or an empty/NULL result with no error

**Cause.** `readModelDb()` returns the model **function**, not an
`rxUi`. Calling it (`f()`) yields the ui, but treating the function
itself as the ui does not:

```r
# WRONG — errors:
mod <- readModelDb("Some_2024_model")
mod$iniDf          # closure is not subsettable
mod$theta          # closure is not subsettable

# WRONG — SILENT, returns NULL and renders an empty table:
attr(readModelDb(nm), "population")
environment(readModelDb(nm))$population

# RIGHT — either of these:
readModelDb(nm)()$population              # needs library(rxode2) attached
rxode2::rxode(readModelDb(nm))$population # works anywhere
```

**Fix.** Resolve to the ui exactly once, then use it:

```r
ui <- rxode2::rxode(readModelDb("Some_2024_model"))
ui$iniDf; ui$theta; ui$population; ui$reference; ui$description
```

When a vignette holds several models in a list, convert at load
time (`lapply(nms, \(n) rxode2::rxode(readModelDb(n)))`) so every
downstream accessor works.

The silent variants matter more than the erroring ones: `attr()`
returns `NULL` without failing, so `str(pop)` prints ` NULL` and a
table renders empty while the vignette reports success.

## 8. `Column 'id' doesn't exist`, or an all-`NA` prediction from a solve that clearly worked

**Cause.** `rxSolve()` **omits the `id` column entirely** when the
event table holds a single subject. Code that indexes on it then
breaks, or worse, silently yields `NA`:

```r
out <- rxode2::rxSolve(model, ev, returnType = "data.frame")
out$viability[match(seq_along(dox), out$id)]  # match(1L, NULL) -> NA
```

**Fix.** Handle the single-subject case explicitly, and check
rather than assume:

```r
if (is.null(out$id)) {
  stopifnot(nrow(out) == length(dox))
  return(out$viability)
}
idx <- match(seq_along(dox), out$id)
stopifnot(!anyNA(idx))
out$viability[idx]
```

If downstream code (PKNCA grouping, `distinct(id, ...)`) needs the
column, restore it right after the solve:
`if (is.null(sim$id)) sim$id <- 1L`.

## 9. `rep(0, dim(.omega)[1]) : invalid 'times' argument`

**Cause.** `omega = NA` was passed to a model that declares **no**
`eta` terms. rxode2 evaluates `dim(NA)[1]`, which is `NA`, and
calls `rep(0, NA)`.

**Fix.** `omega = NA` suppresses IIV, so it is only meaningful when
there is IIV to suppress. Pass it conditionally when one helper
serves both typical-value and IIV models:

```r
solve_typical <- function(mod, ev, ...) {
  ui <- rxode2::rxode(mod)
  if (any(!is.na(ui$iniDf$neta1))) {
    rxode2::rxSolve(mod, ev, omega = NA, ...)
  } else {
    rxode2::rxSolve(mod, ev, ...)
  }
}
```

## 10. A validation check that passes without testing anything

**Cause.** A lookup that matches no rows returns a zero-length
vector, and `all(logical(0))` is `TRUE`. A conclusions table then
reports "yes" for a claim it never evaluated. The usual trigger is
a label mismatch — including one from `format()`, which is
**vectorised** and pads to a common precision:

```r
format(c(0.5, 4, 6), trim = TRUE)   # "0.5" "4.0" "6.0"  <- 4 became "4.0"
```

so a generated `"1 g q8h (4.0 h)"` silently fails to match a
hand-written `"1 g q8h (4 h)"` elsewhere in the same file.

**Fix.** Format each element independently
(`vapply(x, format, character(1), trim = TRUE)`), derive both label
sets from one source, and make lookups fail loudly:

```r
cell <- function(reg, mic) {
  v <- tab$value[tab$regimen == reg & tab$MIC == mic]
  if (length(v) != 1L) stop("no unique row for '", reg, "' at MIC ", mic)
  v
}
stopifnot(all(wanted_labels %in% tab$regimen))  # guard %in% filters too
```

Whenever a check "passes", confirm it had rows to test. A gate
that cannot go red is worse than no gate.

## 11. A validation assertion that fails while the model is correct

Before changing a model to satisfy an assertion, check that the
assertion measures what it claims. Recurring cases:

* **Terminal half-life.** Time-to-50% measured from the moment of
  withdrawal includes the washout transient and reads long. Fit the
  slope over a window well after washout to recover `log(2)/k`.
* **NCA against Dose/CL.** A time grid that does not resolve `Tmax`
  understates AUC by several percent. Sample the distribution phase
  finely.
* **Concentrations decayed into solver noise** go slightly negative
  in the far tail; PKNCA then takes `log()` of a negative value and
  `aucinf.obs` is `NaN`. Truncate to a window that is still many
  half-lives long, and assert `all(conc >= 0)`.
* **A median across subjects is not the typical-value prediction.**
  With large IIV on a non-linear response, the median of the
  transform sits above the transform of the median. Assert the
  typical-value quantity against a typical-value threshold.

When the model turns out to be right, **tighten** the assertion to
the accuracy actually achieved rather than loosening it — that is
what makes it catch a future regression.

**That rule applies to DETERMINISTIC quantities only.** For anything
computed from a simulated cohort, "the accuracy actually achieved"
is one draw on one machine, and tightening to it produces an
assertion that fails everywhere else. See pattern 12.

## 12. An assertion that passes locally and fails in CI

**Cause.** `rxSetSeed()` makes an rxode2 simulation reproducible for
a **given number of solver threads**, not across thread counts: the
parallel RNG streams are partitioned per thread. A dev box solving
on 16 threads and a CI runner on 2 therefore draw *different
cohorts* from identical source. `set.seed()` does not help — it
seeds R's RNG, not rxode2's.

Any assertion on a cohort-derived quantity that is tight enough to
sit inside that spread will pass where it was written and fail
where it runs. This is the single largest source of CI-only vignette
failures: a 2026-08-31 consolidation shipped **twelve** of them, and
CI reported only four because each render shard aborts on its first
failure.

**Three shapes to never write on a simulated quantity:**

```r
# 1. The SIGN or ORDERING of a near-zero effect.
stopifnot(slope < 0)                    # flat effect -> sign is a coin flip
stopifnot(all(diff(pv) < 0))            # every adjacent pair ordered
stopifnot(all(x > 0))                   # one arm sits on zero
stopifnot(spread_C < spread_A)          # two noisy ranges raced

# 2. EXACT zero or equality.
stopifnot(all(pct_over_limit == 0))     # "no subject exceeds" = one draw

# 3. A bound with no headroom over one observed run.
stopifnot(abs(pct_diff) < 2)            # 2 was simply what that run gave
```

**Write instead:**

```r
# Magnitude, not sign, when the claim is that an effect is small:
stopifnot(abs(slope) < 8)
# Trend, not step-by-step monotonicity:
stopifnot(pv[length(pv)] < pv[1])
# A tolerance that admits the noise, not exact zero:
stopifnot(mean(pct_over_limit) < 2)
# An ABSOLUTE bound the paper itself states, not a race between two
# noisy statistics:
stopifnot(spread_C <= 15)   # paper reports 74/72/65%, a 9-point spread
```

**Choosing a bound.** Never take it from one run. Render the vignette
at more than one thread count and set the bound outside the range you
observe, then record that range in a comment so the next reader does
not "tighten" it back:

```r
# Realised 2.99 / 8.23 / 5.35 at 2 / 4 / 16 threads; 8 sat inside the
# noise. 12 still breaks on a mis-transcribed volume, dose or unit,
# which move Cmax by tens of percent.
stopifnot(max(abs(cmax_pct)) < 12)
```

Check the loosened form can still go red. Widening a bound on a
**proportion** past 100, or on a percentage past its achievable
range, produces a gate that cannot fail — worse than no gate
(pattern 10).

**When the assertion is right and the model disagrees, record the
deviation — do not widen until it passes.** If a claim is
reproducibly outside tolerance rather than flickering at its edge,
flag it and exclude it from the gate, keeping it visible in the
rendered table:

```r
claim <- function(text, value, achieved, deviation = FALSE) { ... }
stopifnot(all(tab$Pass[!tab$Deviation]))
```

and say in prose what is not reproduced, with the measured values
and the likely mechanism. A vignette that documents a known
disagreement is worth more than one whose gate was widened until the
disagreement disappeared.

**Do NOT "fix" this by pinning threads inside the vignette.** The
cohort would then depend on a call the reader has to notice, and the
published numbers would still be one draw. Write assertions that
hold for any cohort the model can produce.

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
