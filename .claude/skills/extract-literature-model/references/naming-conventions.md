# Naming conventions for nlmixr2lib models

Authoritative source: `vignettes/create-model-library.Rmd`. This file collects the conventions relevant to model extraction and adds the standards agreed on when the `extract-literature-model` skill was introduced. When in doubt, prefer this file; if this file conflicts with `create-model-library.Rmd`, raise the conflict with the user rather than silently picking one.

## Compartments

Lower case. Snake case only when combining concepts.

- `depot` — extravascular dosing compartment (oral, SC, IM).
- `central` — IV / central sampling compartment.
- `peripheral1`, `peripheral2` — peripheral compartments for 2- and 3-compartment models.
- `effect` — effect compartment for PK/PD models.
- `transit1`, `transit2`, … — transit-compartment absorption chains.
- `target` — free (unbound) target species in explicit-binding TMDD models.
- `complex` — drug–target complex in explicit-binding TMDD models (Mager & Jusko 2001).
- `total_target` — total (free + bound) target in QSS / MM TMDD approximations where the bound species is not carried as a separate state (Gibiansky et al. 2008).
- Therapeutic-area or mechanism-specific compartments: open a GitHub issue before adding new names.

The `target` / `complex` / `total_target` names follow the convention proposed by @iamstein in review of PR #60 and are standard in the TMDD literature.

### Do not declare compartment order

Do **not** add explicit `cmt(depot)`, `cmt(central)`, … declarations at the top
of `model()`. rxode2 / nlmixr2 infer compartment order from the order the
`d/dt(...)` lines appear, and user data can set the compartment via the `cmt`
column on dose rows. Explicit `cmt()` declarations are almost never needed; they
clutter the model and can mask ordering bugs. Only add them when the source
paper *requires* a specific numeric order that isn't what the ODE declarations
would imply (very rare) and flag the reason in a comment.

## Observation variable

- `Cc` — concentration in the central compartment. Use even with `linCmt()`: `Cc <- linCmt()`.
- Multi-output models add named output variables (e.g., `Cbrain_cerebellum`, `tumorSize`). Apply residual error to each output by name.

## Structural PK parameters

Log-transform any parameter that must be positive. Prefix `l` (lower-case L).

- `lka` — absorption rate
- `lcl` — clearance
- `lvc` — central volume of distribution
- `lvp`, `lvp2` — first / second peripheral volume
- `lq`, `lq2` — inter-compartmental clearance to `peripheral1` / `peripheral2`
- `lfdepot` — log bioavailability for the depot compartment

Derived quantities inside `model()` are un-prefixed:
`ka`, `cl`, `vc`, `vp`, `vp2`, `q`, `q2`, and micro-constants `kel` (= `cl/vc`), `k12`, `k21`, `k13`, `k31`.

## Transform prefixes

- `l` — natural log (`lka`, `lcl`, `lfdepot`)
- `logit` — logit (`logitfr`, `logitemax`)
- `probit` — probit
- Any other transform: spell it out as a prefix

Always label every parameter inside `ini()` with units and a short interpretation.

## Inter-individual variability (IIV)

Prefix `eta` + the **transformed** parameter name. Example: IIV on `lcl` is `etalcl`.

- Single IIV: `etalcl ~ 0.09`
- Correlated IIV (block): `etalcl + etalvc ~ c(var_cl, cov, var_vc)`
- Fixed: wrap in `fixed()` when the source reports a known value.

Do **not** use `iiv_`, `IIV_`, or `bsv_` prefixes for new models. The two most recent existing models (`Clegg_2024_nirsevimab`, `Hu_2026_clesrovimab`) write `etacl` without the `l`; this skill standardizes on `etalcl` going forward. Existing files are not migrated.

Document `omega²` = `log(CV² + 1)` in a comment when translating CV% to the internal variance scale.

## Covariate-effect parameters

Multiplicative additive effect: `e_<cov>_<param>`. Examples: `e_wt_cl`, `e_ada_cl`, `e_black_cl`.

Power-style effect: keep the same naming and apply as `cov^e_cov_param` if the source parameterizes it that way; comment the form.

## Endogenous / mechanistic parameters

For endogenous turnover, steady-state, and enzyme-kinetic models (e.g., `igg_kim_2006`, `phenylalanine_charbonneau_2021`), parameters come from mechanism rather than from a CL/V parameterization. Use the names the paper uses; lower-case snake-case by default. Log-transform only positive-constrained parameters that are being *estimated* — not mechanistic constants that the source paper reports as point values.

Recommended patterns:

- `Vmax`, `Km` — Michaelis–Menten constants for each enzyme / transporter. Name-disambiguate when several coexist: `vmax_pah`, `km_pah`, `vmax_trans`, `km_trans`.
- `kint`, `kcat`, `kpro`, `krmr` — fractional rate constants (1/time). When a rate is recomputed at steady state vs. dynamically, suffix the steady-state value with `_0` (e.g., `kcat_0`) and leave the dynamic one unsuffixed.
- `bl_<species>` — baseline concentration of an endogenous species (e.g., `bl_phe`, `bl_gut`). Use this as the initial condition: `<species>(0) <- bl_<species>`.
- `f_<fraction>` — unitless fractional-activity scalars (e.g., `f_pah` = fraction of healthy PAH activity).
- `vd` — body-weight-normalized volume of distribution (L/kg) when the paper uses it that way. Don't rename to `vc` if the paper's mechanism makes `vd` meaningful.

Constants spelled out in the `model()` block (molecular weights, stoichiometric conversion factors, reference weights) should sit at the **top** of `model()` before any derived quantity, with a unit comment.

Endogenous models typically have:

- **No `eta*` IIV parameters.** The model is a typical-value mechanism.
- **No residual error.** Deterministic simulation is the intended use.
- **No dosing events.** The state starts at biological baseline.

When the paper *does* report variability or residual error for an endogenous model, follow the standard IIV / residual conventions below.

## Residual error

- Proportional: `propSd`
- Additive: `addSd`
- Combined: `Cc ~ add(addSd) + prop(propSd)`

For multi-output models, prefix with the output variable: `CcpropSd`, `CcaddSd`, `tumorSizepropSd`. Note in a comment when the source reports log-additive error (NONMEM `EPS(1)` on log-transformed observation); that maps to proportional error in nlmixr2.

## Covariate column names

Not listed here. See `inst/references/covariate-columns.md` — the authoritative register. Any covariate used in a model must exist in that register; if it doesn't, propose adding an entry (skill workflow step).

## File-level metadata

Every model function body begins with:

```r
description <- "<One-sentence summary used in model listings>"
reference   <- "<Full citation including DOI>"
units       <- list(time = "<unit>", dosing = "<unit>", concentration = "<unit>")
covariateData <- list(
  <NAME> = list(
    description       = "<what it is>",
    units             = "<unit or (binary)>",
    type              = "<continuous | binary | categorical | count>",
    reference_category = "<for categorical/binary; the 0 group>",
    notes             = "<derivation, time-varying vs fixed, etc.>",
    source_name       = "<the column name used in the source paper, if different>"
  )
)
population <- list(
  n_subjects   = <integer>,
  n_studies    = <integer>,
  age_range    = "<e.g., 0-24 months>",
  weight_range = "<e.g., 2-12 kg>",
  sex_female_pct = <numeric>,
  race_ethnicity = c(White = <pct>, Black = <pct>, Asian = <pct>, Other = <pct>),
  disease_state = "<e.g., healthy infants at risk for RSV>",
  dose_range    = "<e.g., 50-200 mg IM>",
  regions       = "<e.g., North America, EU, Japan>",
  notes         = "<free text>"
  # Additional keys permitted — add any other important population details
)
```

After this metadata block come `ini()` and `model()`.

## File naming

- Path: `inst/modeldb/<category>/<FirstAuthor>_<Year>_<drug>.R`
- Function name **must** equal the filename minus `.R`. Enforced by `buildModelDb()`.
- Use the paper's first-author surname (no accents or spaces), four-digit year, and drug INN in lowercase.
