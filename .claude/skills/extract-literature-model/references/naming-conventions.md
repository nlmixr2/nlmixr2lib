# Naming conventions for nlmixr2lib models

Authoritative source: `vignettes/create-model-library.Rmd`. This file collects the conventions relevant to model extraction and adds the standards agreed on when the `extract-literature-model` skill was introduced. When in doubt, prefer this file; if this file conflicts with `create-model-library.Rmd`, raise the conflict with the user rather than silently picking one.

## Compartments

Lower case. Snake case only when combining concepts.

- `depot` — extravascular dosing compartment (oral, SC, IM).
- `central` — IV / central sampling compartment.
- `peripheral1`, `peripheral2` — peripheral compartments for 2- and 3-compartment models.
- `effect` — effect compartment for PK/PD models.
- `transit1`, `transit2`, … — transit-compartment absorption chains.
- `precursor1`, `precursor2`, … — precursor / maturation chains (e.g., platelet maturation in Petrov 2024 romiplostim).
- `lat0`, `lat1`, `lat2`, … — latency chains (e.g., CSF/BBB delay in LeTilly 2021 trastuzumab).
- `target` — free (unbound) target species in explicit-binding TMDD models.
- `complex` — drug–target complex in explicit-binding TMDD models (Mager & Jusko 2001).
- `total_target` — total (free + bound) target in QSS / MM TMDD approximations where the bound species is not carried as a separate state (Gibiansky et al. 2008).
- Therapeutic-area or mechanism-specific compartments: open a GitHub issue before adding new names.

The `target` / `complex` / `total_target` names follow the convention proposed by @iamstein in review of PR #60 and are standard in the TMDD literature.

### Numbered compartments and chain prefixes

Generic numbered compartments (`central1`, `central2`, `cmt1`, `comp_a`) are
**forbidden** — use the canonical `central` / `peripheral1` / `peripheral2`
names. The only numeric suffixes the conventions allow are:

1. The blessed chain prefixes above: `transit<n>`, `effect<n>`,
   `precursor<n>`, `lat<n>` — used for chained / serial compartments.
2. DAR-numbered ADC isoforms: `dar<n>_<base>` (e.g., `dar0_central`,
   `dar3_peripheral1`). Each `dar<n>` represents a biologically distinct
   drug-antibody-ratio isoform of the parent ADC, not arbitrary numbering.
3. Metabolite-suffixed compartments: `<canonical>_<metab>` (see § Parent
   drug + metabolite below).

### Parent drug + metabolite compartments (ADCs and similar)

For models that track a parent drug AND one or more metabolites
(antibody-drug conjugates with payload, mAbs that release a small
molecule, etc.), the **parent uses the canonical names** and the
**metabolite carries a lowercase suffix** drawn from the registered
metabolite list:

```
# Parent ADC
central
peripheral1
peripheral2

# MMAE metabolite (suffix at end)
central_mmae
peripheral1_mmae
```

Registered metabolite suffixes (defined in `R/conventions.R::registeredMetabolites`):
`mmae`, `dxd`, `sn38`, `dm4`, `medm4`, `mcmmaf`, `complex`, `ige`, `tab`,
`nab`, `dar0..dar8`. Add a new metabolite by editing
`R/conventions.R::registeredMetabolites`; do not introduce ad hoc
prefixes (`adc_central`, `mmae_central`).

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
- Multi-output PK metabolite outputs use `Cc_<metab>` (e.g., `Cc_mmae`,
  `Cc_dxd`, `Cc_sn38`, `Cc_tab`, `Cc_dar0`, `Cc_dm4`). Do **not** use the
  deprecated `C<metab>` form (`Cmmae`, `Cdxd`).
- Non-PK multi-output variables (e.g., `tumorSize`, `freeIgE`, `totalIgE`,
  `Cbrain_cerebellum`, `Ccsf`) are paper-named and exempt from the
  `Cc_<metab>` rule.
- Apply residual error to each output by name. Multi-output residual
  error is `<output>propSd` / `<output>addSd` (e.g., `CcpropSd`,
  `Cc_mmaepropSd`, `tumorSizepropSd`).

## Structural PK parameters

Log-transform any parameter that must be positive. Prefix `l` (lower-case L).

- `lka` — absorption rate
- `lcl` — clearance
- `lvc` — central volume of distribution
- `lvp`, `lvp2` — first / second peripheral volume
- `lq`, `lq2` — inter-compartmental clearance to `peripheral1` / `peripheral2`
- `lfdepot` — log bioavailability for the depot compartment
- `lvmax` — Michaelis-Menten Vmax (for nonlinear / saturable elimination)

Derived quantities inside `model()` are un-prefixed:
`ka`, `cl`, `vc`, `vp`, `vp2`, `q`, `q2`, `vmax`, and micro-constants
`kel` (= `cl/vc`), `k12`, `k21`, `k13`, `k31`.

**Volumes**: always `vc` / `vp` / `vp2`. Never `v`, `v1`, `v2`, `vc1`,
`vc2`, etc. The `v1`/`v2`/`v3` numbering used by NONMEM `linCmt()` macros
maps directly: `v1` → `vc`, `v2` → `vp`, `v3` → `vp2`.

**Michaelis-Menten Vmax**: always `lvmax` / `vmax`. Never `lvm` / `vm`.

### Multi-component clearance

Some models report two CL terms (a steady-state component and a
time-varying decay, or a Michaelis-Menten arm plus a linear arm). Use
the registered CL component suffixes:

- `lcl_ss` — steady-state arm of CL.
- `lcl_time` — time-varying decay arm of CL.

Inside `model()` the bare names are `cl_ss`, `cl_time`. Source-paper
terms like `clinf`, `clss`, `clt`, `clD` map onto these:

- `clinf` (terminal / asymptotic CL), `clss` (steady-state CL) → `cl_ss`.
- `clt` (time-varying CL component) → `cl_time`.
- `cld` (distributional clearance, i.e., inter-compartmental) → `lq`.

### Parent drug + metabolite parameters

Same rule as compartments: the **parent uses the canonical name**, and
the **metabolite has the lowercase paper-name suffix** at the end:

```
# Parent ADC
lcl, lvc, lq, lvp                 # canonical, unsuffixed
e_wt_cl, e_alb_cl, e_wt_vc        # parent covariate effects

# MMAE metabolite
lcl_mmae, lvc_mmae, lq_mmae       # _mmae suffix
e_wt_cl_mmae, e_alb_cl_mmae       # metabolite covariate effects
```

Do **not** suffix the parent with `_adc`, `_dato`, `_t-dm1`, etc. — the
parent always wins canonical naming.

Multi-analyte models (e.g., `Sathe_2024_sacituzumab` with SG ADC, free
SN-38, and total antibody tAB) treat the ADC as the parent and each
non-parent analyte as a distinct metabolite suffix:

```
central, peripheral1                              # SG ADC parent
central_sn38, peripheral1_sn38                    # free SN-38 metabolite
central_tab,  peripheral1_tab                     # total antibody analyte
```

### DAR-explicit mechanistic ADCs

For models that carry every drug-antibody-ratio isoform as a separate
species (Bender_2014 T-DM1 mechanistic, Pouzin_2022 tusamitamab),
compartments are `dar<n>_<base>` and outputs are `Cc_dar<n>`:

```
dar0_central, dar0_peripheral1, dar0_peripheral2
dar1_central, dar1_peripheral1, dar1_peripheral2
...
dar7_central, dar7_peripheral1, dar7_peripheral2

Cc_dar0, Cc_dar1, ... Cc_dar7
```

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

Canonical form: `e_<cov>_<param>` (`<cov>` first, `<param>` second). Examples:
`e_wt_cl`, `e_ada_cl`, `e_black_cl`. Never reverse the order
(`e_cl_wt` is wrong).

Power-style effect: keep the same naming and apply as `cov^e_cov_param` if
the source parameterizes it that way; comment the form.

### Three-token covariate effects

When a covariate effect is tied to a specific metabolite, a paired
parameter (shared exponent), or a multi-component CL arm, append the
appropriate token at the end:

- **Metabolite-suffixed**: `e_<cov>_<param>_<metab>`. Examples:
  `e_wt_cl_mmae`, `e_alb_cl_dxd`, `e_creat_cl_sn38`. The `<metab>` token
  must be one of the registered metabolites in
  `R/conventions.R::registeredMetabolites`.
- **Shared exponent (single estimated value applied to two parameters)**:
  `e_<cov>_<param1>_<param2>`. Examples: `e_wt_cl_q` (shared WT exponent on
  CL and Q), `e_wt_vc_vp` (shared WT exponent on Vc and Vp). The
  `<param2>` token must be a bare PK parameter
  (`pkBareParams`). Source-paper terms like `clq`, `vss`, `vcvp` translate
  to: `clq` → `cl_q`; `vss` → `vc_vp` (Vss = Vc + Vp); `vcvp` → `vc_vp`.
- **Multi-component CL**: `e_<cov>_cl_ss` or `e_<cov>_cl_time`. Examples:
  `e_wt_cl_ss`, `e_alb_cl_time`. Source-paper terms `clinf` / `clss` →
  `cl_ss`; `clt` → `cl_time`.

### Compound covariate names

Some canonical covariates contain underscores (e.g., `RACE_BLACK`,
`RACE_ASIAN`, `ADA_TITER`, `PRIOR_GAST`). The covariate-effect name
preserves the underscore: `e_race_black_cl`, `e_race_asian_vc`,
`e_ada_titer_cl`. The `checkModelConventions()` covariate-effect regex
accepts these without complaint.

### Categorical / binary covariate effects

When the same parameter is modified by a set of binary indicators (e.g.,
race indicators), give each indicator its own covariate-effect parameter:

```r
e_asian_cl       # CL effect for Asian (vs reference)
e_black_cl       # CL effect for Black (vs reference)
e_multiracial_cl # CL effect for Multiracial (vs reference)
```

Do not lump them under a single `e_race_cl` (ambiguous); do not omit the
parameter token (`e_asian` is wrong — it must end in `_cl` or `_vc` or
similar to identify what it modifies).

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
