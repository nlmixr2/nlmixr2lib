# SI-units follow-up TODO (2026-06-19 audit, PR 6 partial scope)

## Status

PR 6 of the 2026-06-19 canonical-register standardization audit standardized the
**canonical Units field** to SI for the clinical-chemistry covariates `ALB`,
`TPRO`, `CSF_TPRO`, `IGG`, `IGM`, `TBILI`, `DBIL` (US-convention papers must now
convert their data to SI on ingestion, or apply an inline conversion in their
model file). The enzyme canonicals `AST`, `ALT`, `ALP`, `GGT`, `LDH`, `CPK`
were already in U/L (SI) and the Units field was annotated to confirm.

**Not yet done (follow-up work):** per-model inline unit-conversion lines in
the 94 model files that consume one or more of these covariates. Each affected
model that was calibrated to a US-convention value (e.g., albumin in g/dL) must
add an inline conversion in the `model()` block so the structural coefficients
stay aligned with their original calibration, AND each conversion must be
validated with a before-vs-after simulation equivalence check.

## Operator directive (verbatim from the 2026-06-19 audit)

> Enforce SI units. 100% validate that the model is unchanged with unit
> conversions by running it without/with the conversion. Do the unit conversion
> in the math in model or ini block without hiding it. For example
> `tbili_mgdL <- TBILI*17.1 # unit conversion to conventional units` in the
> model block and then use that tbili_mgdL in the equations.

## Per-model conversion pattern

For each affected model file, the conversion pattern is:

```r
model({
  # ... existing variable definitions ...

  # Inline unit conversion: SI canonical column -> US-convention value
  # used in the original published calibration.
  alb_gdL   <- ALB    * 0.1      # SI g/L     -> US-convention g/dL
  tbili_mgdL <- TBILI  / 17.1    # SI umol/L  -> US-convention mg/dL
  igg_mgdL  <- IGG    * 100      # SI g/L     -> US-convention mg/dL
  # ... etc.

  # Use the converted variable in downstream equations:
  cl <- exp(lcl + etalcl) * (alb_gdL / 4.0)^e_alb_cl
})
```

**Conversion factors:**

| Canonical | SI unit (post-audit) | US-convention | Conversion |
|---|---|---|---|
| ALB | g/L | g/dL | g/dL = g/L * 0.1 |
| TPRO | g/L | g/dL | g/dL = g/L * 0.1 |
| CSF_TPRO | g/L | mg/dL | mg/dL = g/L * 100 |
| IGG | g/L | mg/dL | mg/dL = g/L * 100 |
| IGM | g/L | mg/dL | mg/dL = g/L * 100 |
| TBILI | umol/L | mg/dL | mg/dL = umol/L / 17.1 |
| DBIL | umol/L | mg/dL | mg/dL = umol/L / 17.1 |

## Validation protocol (per affected model)

For EACH affected model:

1. Identify the unit used in the original published calibration (read the model
   file's `covariateData[[X]]$units` field).
2. If already SI, no per-model edit needed — verify the canonical Units field
   is now SI in the register and proceed.
3. If US-convention, add the inline conversion in `model()` block per the pattern
   above and update all downstream equation references to use the converted
   variable.
4. Update `covariateData[[X]]$units` to reflect the canonical SI unit, with a
   `notes` line that says "Input column expected in SI g/L; converted inline to
   g/dL for use in the original Table-N coefficients."
5. **Validation:** run a simulation BEFORE the conversion change and AFTER.
   - BEFORE: use the original input values in the original units (e.g.,
     `ALB = 4.0` g/dL).
   - AFTER: use SI-equivalent input values (e.g., `ALB = 40` g/L).
   - Verify the predicted Cc / output trajectories are identical within
     numerical tolerance (e.g., `max(abs(Cc_before - Cc_after) / Cc_before) <
     1e-6`).

## Affected model files (94 total)

These files contain a `covariateData[[X]] = list(...)` entry for at least one of
the SI-standardized canonicals (ALB, TPRO, CSF_TPRO, IGG, IGM, TBILI, DBIL).
Many also use the enzyme canonicals (AST, ALT, ALP, GGT, LDH, CPK), which are
already in U/L (SI) and need no conversion — but their `notes` field may need a
"confirmed SI" annotation when this follow-up sweep runs.

Files grouped by primary affected canonical:

### Files using ALB

(Most of the 94 files; see `git grep "ALB = list" inst/modeldb/specificDrugs/`
for the authoritative list.)

### Files using TBILI

(Sub-set; see `git grep "TBILI = list" inst/modeldb/specificDrugs/`.)

### Files using IGG

`Zhou_2021_belimumab.R`, `Yang_2021_cemiplimab.R`, `Struemper_2017_belimumab.R`,
`Cheng_2026_immunoglobulin.R`, `Lin_2024_casirivimab.R`, ...

### Files using DBIL

`Chen_2015_voriconazole.R`, plus any future DBIL-consuming papers.

### Files using IGM

`Cheng_2026_immunoglobulin.R`.

### Files using TPRO

`Frey_2010_tocilizumab.R`.

### Files using CSF_TPRO

`Germovsek_2018_meropenem.R`.

## Recommended split for follow-up PRs

Given 94 affected files, a reasonable split would be:

* **PR 6a (this PR)** — Register Units field standardized to SI. No model
  file changes. Documents the follow-up work needed.
* **PR 6b** — Inline conversions for the ~10 highest-impact files (single-source
  papers with explicit US-convention values, full validation per the operator
  directive). Each PR runs the before/after simulation equivalence check.
* **PR 6c+** — Continue model-file batches in groups of ~10-15 files per PR.

Each downstream PR (6b, 6c, ...) MUST:

* Read the per-model `covariateData[[X]]$units` field to identify the
  calibration unit.
* Apply the inline conversion in `model()` block.
* Run a before/after simulation with equivalent input values and verify
  numerical-tolerance equivalence of outputs.
* Document the validation result in the per-PR deliverable report.
