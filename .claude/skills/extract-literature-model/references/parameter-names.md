# Parameter naming for nlmixr2lib models

**Machine-readable source of truth:** `inst/references/parameter-names.md` in the `nlmixr2lib` package. That file is parsed at runtime by `R/conventions.R::.parseTypedNamesMd` and seeds the canonical `pkParams`, `pkBareParams`, and `paperNamedParams` lists used by `checkModelConventions()`. New canonical parameter names MUST be added there; this skill reference is the user-facing summary kept in sync with that file. Companion references: `vignettes/create-model-library.Rmd` (narrative). Covers structural PK parameters, transform prefixes, fixed parameters, IIV, covariate-effect parameters, endogenous / mechanistic parameters, residual error, and file-level metadata.

**Stop-and-ask gate (Phase 1 pre-flight + Phase 3 drafting):** If the model you are extracting needs a parameter name that is NOT in this document, file a sidecar BEFORE writing the model file. Propose: the canonical name (with `l` prefix if log-transformed), its role (one sentence), source paper's local name(s) it would replace, why it isn't an alias of an existing canonical (e.g., why this isn't just `cl_ss` / `cl_time` / `cl_renal` / `cl_nonren` under a different label), and any cross-precedent in existing registered model files. Wait for operator approval before committing. Trivial notation differences (case-only, NONMEM `V1`/`V2`/`V3` → `vc`/`vp`/`vp2`, paper's `Kel` → canonical `kel`) translate silently and do NOT need a sidecar.

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

### K-PD / single-rate-constant parameterisation

Some popPK papers report a single first-order elimination rate constant
(K-PD or one-compartment-without-explicit-V parameterisation) rather than the
`cl`/`vc` decomposition. In these models the elimination rate constant is a
**primary** `ini()` parameter, not a derived micro-constant:

- `lkel` — log first-order elimination rate constant (K-PD / single-rate-constant form).
- Inside `model()` the bare name is `kel`. When `vc` is also estimated, use the
  derived form `kel <- cl/vc` (canonical) — `lkel` as a primary parameter is
  reserved for the K-PD case where no explicit `vc` exists.

Never use `lke` for this role; standardise on `lkel`.

### Multi-component clearance

Some models report two CL terms (a steady-state component and a
time-varying decay, or a Michaelis-Menten arm plus a linear arm). Use
the registered CL component suffixes:

- `lcl_ss` — steady-state arm of CL.
- `lcl_time` — time-varying decay arm of CL.
- `lcl_renal` / `lcl_nonren` — renal vs non-renal arm (e.g., cefepime additive renal-plus-non-renal).

Inside `model()` the bare names are `cl_ss`, `cl_time`, `cl_renal`, `cl_nonren`. Source-paper terms like `clinf`, `clss`, `clt`, `clD` map onto these:

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

## Absorption: transit and lag-time

### Transit-absorption chain (Savic parameterisation)

- `lmtt` — log mean transit time.
- `lktr` — log first-order transit rate constant (= n_transit / MTT for a chain of length n).

Inside `model()` the bare names are `mtt` and `ktr`. Source-paper aliases
(`MAT`, `MTT`, `KTR`) translate silently.

### Lag-time

Canonical lag-time prefix: **`ltlag`** (log absorption-lag time).

Source-paper aliases that translate to `ltlag` without sidecar:
`ALAG1` (NONMEM), `tlag`, `Tlag`, `alag`, `LAG`. The legacy forms
`lalag`, `llag` are deprecated in favour of `ltlag`.

Inside `model()` the bare name is `tlag`. Apply via `alag(depot) <- tlag`
or `alag(<cmt>) <- tlag` (preferred over carrying a separate `lag` compartment).

### Acrophase / circadian peak time

Canonical acrophase prefix: **`ltacro`** (log time of peak in a
circadian-rhythm rate constant). Distinct from absorption-lag time:
`tacro` is a phase shift inside a sinusoidal modulation of `kin` /
`kout` rather than a delay between dose administration and absorption.

Used in `indirect_circ_*` circadian-IDR templates with kinetic forms
such as `kout_t <- kin + amp * sin(2*pi*(t - tacro) / period)` (and
the cos / amplitude-ratio variants). The legacy form `ltz` is
deprecated in favour of `ltacro`.

## Indirect-response (IDR) / turnover parameters

Canonical rate-constant prefixes for indirect-response and turnover-model
families (Dayneka 1993; Jusko & Friberg traditions). Use these as
primary `ini()` parameters when the paper reports the IDR / turnover
rate constants directly:

- `lkin` — log zero-order production rate constant (synthesis into a turnover pool).
- `lkout` — log first-order elimination rate constant of the turnover pool.
- `lkdeg` — log first-order degradation rate constant (synonym for elimination in some papers).
- `lksyn` — log zero-order synthesis rate constant (alternate name for `kin`).
- `lkpin` — log zero-order rate constant for production of a **precursor** pool (used in `indirect_prec_*` precursor-pool models).
- `lkpout` — log first-order rate constant for loss of a **precursor** pool (used in `indirect_prec_*` precursor-pool models).

Inside `model()` the bare names are `kin`, `kout`, `kdeg`, `ksyn`,
`kpin`, `kpout`. Suffix variants `lkin_<analyte>` etc. are permitted
for multi-analyte / combination-therapy IDR models.

### Baseline values

Canonical baseline-value prefix: **`lrbase`** (log baseline value of a
turnover state, in the units of the state). Use this for the steady-state
"R0" value of IDR / turnover / endogenous-cycle models, regardless of
the paper's local terminology.

Source-paper aliases that translate to `lrbase` without sidecar:
`R0`, `Base`, `BASE`, `BL`, `S0`, `TS0`, `R_baseline`. The legacy forms
`lr0`, `lbl`, `lbase`, `lBase`, `ls0`, `lts0` are deprecated in favour
of `lrbase`.

Inside `model()` the bare name is `rbase`. The compartment initial
condition is set as `<state>(0) <- rbase`. For TGI models with an
estimated initial tumour size, use `lrbase` (per-output baseline can be
suffixed: `lrbase_tumor`, `lrbase_anc`).

## Sigmoidal PD shape parameters

Canonical Hill-coefficient prefix: **`lhill`** (log Hill exponent in
sigmoidal Emax / Imax functions).

The canonical pattern is the sigmoidal Emax / Imax form

```
eff = emax * Cc^hill / (ec50^hill + Cc^hill)
```

(or the Imax / inhibitor analogue). Use `lhill` as the primary `ini()`
parameter whenever the source paper reports a Hill coefficient for a
sigmoidal stimulation / inhibition function — regardless of whether the
paper calls it "γ", "gamma", "hill", or "n".

Inside `model()` the bare name is `hill`.

**Role distinction — `lgamma` is NOT a synonym for `lhill`.** Several
mechanistic roles carry the paper-symbol `γ` but are NOT Hill exponents
and should retain their `lgamma` / `gamma` name:

- **Friberg myelosuppression feedback exponent** in `(Circ0/circ)^gamma`
  — the feedback amplification on circulating-cell count. Keep as `lgamma`.
- **TGI power-law / generalised von Bertalanffy / generalised logistic
  growth exponents** (`tumorSize^gamma`, `kge*tumorSize^gamma`,
  `(tumorSize/tsmax)^gamma`). Keep as `lgamma`.
- **Linear amplification factors** (e.g., Tetschke 2018 EPO feedback,
  `Fb <- gamma * (THB_MASS - thb) / THB_MASS`). Keep as `lgamma`.
- **Death-rate / kill-rate constants** named `γ` in the source paper
  (e.g., Mazzocco 2015 TMZ tumour-cell death-rate constant). Keep as
  `lgamma`.
- **Gamma-distribution shape parameters** for distributional / hazard
  models. Keep as `lgamma`.
- **Power coefficients on transit inputs** (e.g., Ait-Oudhia 2012
  CRP-transit amplification). Keep as `lgamma`.

The `lnn` / `nn_fix` form is reserved for Wang & co-authors'
sigmoidicity exponent in specific BDE / morphine-like models and is a
distinct canonical from `lhill`.

## Paper-specific etas and residual-SDs

For models with genuinely paper-mechanistic IIV / residual-error
parameter names that do not pair 1-to-1 with an `lX` fixed-effect
parameter, declare them via the `paper_specific_etas` and
`paper_specific_residual_sds` metadata fields at the top of the
model function body (analogous to the `depends` field for
upstream-imported covariates and the `paper_specific_compartments`
field for paper-mechanistic compartment names):

```r
my_model <- function() {
  description <- "..."
  reference <- "..."
  paper_specific_etas <- c("etalogit", "etap1", "etap2")
  paper_specific_residual_sds <- c("propSd_vact_l1", "propSd_vact_l2")
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")
  ...
}
```

`checkModelConventions()` subtracts these from the parameter-naming
warning set. Use sparingly: when an eta's underlying typical-value
parameter is part of a paper-mechanistic structural equation
(mixture-probability logit transform, transit-rate fractions,
indirect-response baselines, etc.) rather than a 1-to-1 `lX` ini
parameter; or when a residual SD uses a paper-specific multi-token
output suffix (`propSd_vact_l1` for a lesion-stratified PD output)
that the canonical `propSd_<output>` matcher does not recognise.

## Transform prefixes

- `l` — natural log (`lka`, `lcl`, `lfdepot`)
- `logit` — logit (`logitfr`, `logitemax`)
- `probit` — probit
- Any other transform: spell it out as a prefix

Always label every parameter inside `ini()` with units and a short interpretation.

## Fixed parameters

Wrap a parameter value in `fixed()` whenever the source paper holds it constant rather than estimating it. This applies to **every** parameter type: structural THETAs, allometric exponents, IIV variances and covariances, residual-error magnitudes, covariate-effect coefficients, bioavailability anchors. The `fixed()` wrapper is load-bearing provenance — a downstream user must be able to tell which values are estimated point estimates vs structural assumptions.

Source signals that a parameter is fixed:

- Prose: "fixed during estimation", "fixed at <value>", "held fixed at the literature value", "not estimated", "set to 1 (fixed)".
- Allometric exponents without RSE / SE / CI (canonical 0.75 / 1, especially under "fixed exponents" wording).
- NONMEM `$THETA` / `$OMEGA` / `$SIGMA` with a `FIX` flag.
- Bioavailability `F1 = 1` set as structural anchor.
- Inherited parameters from an upstream publication that the current paper re-uses without re-fitting.

Examples:

```r
lcl       <- log(0.225)            ; label("Clearance (L/h)")        # estimated
lcl       <- fixed(log(2))         ; label("Clearance (L/h)")        # fixed log-transformed value: log() goes inside fixed()
e_wt_cl   <- fixed(0.75)           ; label("Allometric exp on CL")   # fixed
lfdepot   <- fixed(log(1))         ; label("Bioavailability")        # fixed anchor (log of 1 = 0)
etalcl    ~ 0.32                                                     # estimated IIV
etalvc    ~ fixed(0.18)                                              # fixed IIV
CcaddSd   <- fixed(0.10)           ; label("Additive SD (LTBS)")     # fixed residual
```

Note the `fixed(log(2))` form: when a log-transformed structural parameter is fixed (e.g. `lcl` for a paper that holds CL constant at 2 L/h from an upstream publication), the `log()` goes **inside** `fixed()`, not the other way around. `lcl <- log(fixed(2))` is wrong — it would fix the linear-scale `2` and then take its log, which is not what `fixed()` is meant to mark. Always: `<lparam> <- fixed(log(<linear_value>))`.

If a parameter is reported without uncertainty but the paper does not explicitly say "fixed", sidecar-ask the operator before guessing. Mis-encoding fixed-vs-estimated is a real downstream error.

## Inter-individual variability (IIV)

Prefix `eta` + the **transformed** parameter name. Example: IIV on `lcl` is `etalcl`.

- Single IIV: `etalcl ~ 0.09`
- Correlated IIV (block): `etalcl + etalvc ~ c(var_cl, cov, var_vc)`
- Fixed: wrap in `fixed()` when the source reports a known value (see "Fixed parameters" above for the broader guidance applying to all parameter types).

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
- **Multi-component CL**: `e_<cov>_cl_ss` or `e_<cov>_cl_time` or `e_<cov>_cl_renal`. Examples: `e_wt_cl_ss`, `e_alb_cl_time`. Source-paper terms `clinf` / `clss` → `cl_ss`; `clt` → `cl_time`.

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

When the paper *does* report variability or residual error for an endogenous model, follow the standard IIV / residual conventions.

## Residual error

- Proportional: `propSd` (used with `~ prop(...)` or `~ prop(...) + add(...)`)
- Additive: `addSd` (used with `~ add(...)`)
- Log-normal / log-scale: `expSd` (used with `~ lnorm(...)`; a distinct error structure from proportional — the SD applies on the exponentiated scale)
- Combined: `Cc ~ add(addSd) + prop(propSd)`

For multi-output models, the per-output form is `<errorname>_<output>` (parameter name first, output suffix second): `propSd_Ccsf`, `addSd_tumorSize`, `expSd_TotalTau`. The parent observation `Cc` uses the bare suffix-free form (`propSd`, `addSd`, `expSd`).

Note in a comment when the source reports log-additive error (NONMEM `EPS(1)` on log-transformed observation); that often maps to either `propSd` (linear-space proportional) or `expSd` (log-normal lnorm), depending on how the residual is interpreted.

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
  species      = "<required: 'human', 'rat (Sprague-Dawley)', 'mouse (HBCx-9 PDX)', 'beagle dog', 'in vitro (SKBR3 cell line)', etc.; for pooled cohorts list each (e.g., 'human + rat')>",
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

- Path: `inst/modeldb/<category>/<FirstAuthor>_<Year>_<drug>.R`.
- Function name **must** equal the filename minus `.R`. Enforced by `buildModelDb()`.
- Use the paper's first-author surname (no accents or spaces), four-digit year, and drug INN in lowercase.
- **Year-letter collision suffix.** When two extractions resolve to the same `<FirstAuthor>_<Year>_<drug>` name (e.g., two same-author/year/drug entries with different scenarios), append a lowercase letter to the year — `Author_2019a_drug.R`, `Author_2019b_drug.R`. Allocate letters in chronological model-development order when known. Never overwrite an existing file silently.
