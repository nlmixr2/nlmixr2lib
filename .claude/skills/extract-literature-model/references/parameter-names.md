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
- `lvmax` — Michaelis-Menten vmax (for nonlinear / saturable elimination)

Derived quantities inside `model()` are un-prefixed:
`ka`, `cl`, `vc`, `vp`, `vp2`, `q`, `q2`, `vmax`, and micro-constants
`kel` (= `cl/vc`), `k12`, `k21`, `k13`, `k31`.

**Volumes**: always `vc` / `vp` / `vp2`. Never `v`, `v1`, `v2`, `vc1`,
`vc2`, etc. The `v1`/`v2`/`v3` numbering used by NONMEM `linCmt()` macros
maps directly: `v1` → `vc`, `v2` → `vp`, `v3` → `vp2`.

**Michaelis-Menten vmax**: always `lvmax` / `vmax`. Never `lvm` / `vm`.

### PD structural parameters

For effect-compartment, indirect-response, and direct-effect PD models, use:

- `lec50` — log of EC50 (sigmoid Emax / Imax models). Bare form `ec50`.
- `lemax` / `limax` — log of Emax / Imax (maximum drug effect on production / loss).
- `lhill` / `lgamma` — log of Hill / sigmoidicity coefficient.
- `lke0` — log of effect-compartment equilibration rate constant (when an effect compartment equilibrates with the central compartment).
- `lthalfrec` — log of an effect-compartment / indirect-response **recovery half-life** (units of time). Use when a paper parameterises a delay-and-recovery process by its half-life rather than by a rate constant. Inside `model()` the bare name is `thalfrec` and the corresponding first-order rate constant is `krec = log(2) / thalfrec`. Founding example: `deVriesSchultink_2018_trastuzumab_LVEF.R` (recovery half-life T1/2rec from de Vries Schultink 2018 Table 2, where dCeff/dt = Ctrastuzumab - log(2)/T1/2rec * Ceff). When the paper parameterises by rate constant instead, use `lkrec` directly and do not introduce a parallel `lthalfrec`.
- `lke_kpd` — DO NOT use; the canonical K-PD elimination rate is `lkel` (same name as the PK-elimination canonical, with the K-PD context conveyed by the surrounding `depot_kpd` compartment).

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

### Fraction metabolised (`fm`, `fm_<pathway>`)

One identified route: bare `fm` (or `logitfm` when the paper holds it on
the logit scale). When the paper splits the metabolic flux across two or
more **named** routes, use the `fm_<pathway>` family:

```
fm_cyp3a4, fm_cyp3a5, fm_other    # enzyme routes + the residual route
fm_m3g, fm_m6g                    # routes named by the metabolite formed
```

- `<pathway>` is **lowercase**, even when the paper capitalises it
  (`fm_h4`, never `fm_H4`), and reuses the canonical metabolite suffix
  from `compartment-names.md` where the route is named by its product.
- The residual / unspecified route is **singular `fm_other`** — never
  `fm_others` — matching the registered `lkp_other` / `kp_other`.
- A trailing `_frac` marks a **nested** share, i.e. a fraction *of* `fm`
  rather than of total clearance (`fm * fm_ko516_frac` is the share of
  parent clearance reaching that metabolite). Use it only for that case.
- Papers that parameterise a nested split **positionally** keep the
  paper's own `fm1` / `fm2` names; those are not `fm_<pathway>` members.

Adding a new pathway that follows these rules needs **no sidecar** — add
it to the `### fm_...` heading in `inst/references/parameter-names.md`.
Anything else is a build error: `buildModelDb()` fails on an `fm_`-prefixed
name bound in `ini()` or `model()` that is not in that heading.

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

## Gastrointestinal transit and enterohepatic rate constants

First-order rate constants for the gut-lumen mass balance in oral PBPK models
that resolve the gastric and intestinal lumen as explicit states (`stomach`,
`gut_lumen`). All three are log-transformed in `ini()`; bare names inside
`model()` drop the `l` prefix. Units 1/time.

- `lkst` -- **gastric emptying rate constant**, draining the `stomach` depot into
  the intestinal lumen. Source-paper aliases: `Kst`, `kst`, `k_st`, `KST`.
  Founding example: `Ai_2024_ractopamine_goat_pbpk.R` (`Kst = 0.0910 /h`); also
  `Yang_2025_matrine_pig_pbpk.R` (`kst = 0.8545 /h`).
- `lkfec` -- **faecal excretion rate constant** for unabsorbed drug leaving
  `gut_lumen` as faeces. Names the *process* (excretion to faeces) rather than
  the compartment, so it cannot be misread as a general gut-transit rate.
  Source-paper aliases: `ke`, `k_e`, `Kgut`, `kgut`, `kint`, `kF`. Founding
  example: `Yang_2025_matrine_pig_pbpk.R` (`ke = 0.007358 /h`).
- `lkbile` -- **biliary excretion rate constant** for parent drug leaving the
  `liver` compartment, either into `gut_lumen` (enterohepatic recirculation) or
  as a terminal biliary elimination route. Source-paper aliases: `kbi`, `k_bi`,
  `Kbile`, `KbileC`. Founding example: `Zhang_2024_f53b_mouse_pbpk.R` (`KbileC`,
  terminal biliary elimination); also `Yang_2025_matrine_pig_pbpk.R`
  (`kbi = 0.05835 /h`, liver into gut lumen).

Three boundaries to respect:

- `lkbile` is **not** `kbm`. `kbm` is registered as the biliary-*metabolite*
  excretion rate constant (a metabolite leaving a plasma / central compartment);
  `lkbile` is parent drug leaving `liver`. Do not widen `kbm` to cover this case.
- `lkgut` appears in `Yang_2023_diclazuril_chicken_pbpk.R` and
  `Ai_2024_ractopamine_goat_pbpk.R` for the same concept as `lkfec`. It predates
  this section and is an unregistered legacy spelling; use `lkfec` in new models.
- The absorption rate constant out of `gut_lumen` remains the ordinary canonical
  `lka` -- do not mint a gut-specific absorption name.

## Time-varying protein binding

Canonical slope for linear time-varying unbound-fraction models:
**`lbfu`** (log slope of unbound fraction with respect to time, 1/time-unit).
Bare name inside `model()`: `bfu`.

Use this for popPK models that encode an evolving free fraction during the
study window, typically of the form

```
fu(t) = fu_ref + bfu * (t - t_ref)
```

where `fu_ref` is the unbound-fraction anchor at reference time `t_ref` (taken
from external literature, not estimated), and `bfu` is the time slope. The
slope is usually fixed at a literature value (`lbfu <- fixed(log(<value>))`)
because the unbound fraction itself is not measured directly in plasma-
concentration popPK studies; only the resulting time-trends in apparent CL
and V are observed. Inter-individual variability on `bfu` may be estimated
when the source paper reports it.

When the source paper imposes a window of validity for the linear-evolution
assumption (e.g. "fu evolves linearly from 0 to 72 h"), clamp `t` to the
declared window inside `model()`:

```
fu <- fu_ref + bfu * (min(t, t_max) - t_ref)
```

Founding example: `LeJouan_2005_quinine.R` (children with falciparum malaria;
fu = b*(t-36) + 0.15, b fixed at 0.001/h, evolution clamped at t = 72 h).

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

## Body-weight evolution (perioperative / fluid-resuscitation models)

Some surgical / critical-care popPK papers fit a transient body-weight curve in
parallel with the drug PK and use the time-varying body weight `BW(t)` as the
allometric size descriptor on `V` (and, less commonly, on `CL`). The Oualha
2018 paediatric-LT enoxaparin paper parameterises BW(t) with a Hill-type
saturation curve:

```
BW(t) = (BWPREOP + PFA / 1000) * (1 - (1 - fbw) * t^hill_bw / (tbw50^hill_bw + t^hill_bw))
```

When the paper estimates this curve jointly with the drug PK and reports IIV on
its parameters, encode the BW(t)-evolution parameters as log-transformed
primary `ini()` entries:

- `lfbw` -- log asymptotic fractional body-weight retention after fluid loss.
  Bare name `fbw`, dimensionless, typically in (0, 1].
- `lhill_bw` -- log Hill steepness of the BW(t) decline curve. Bare name
  `hill_bw`. Suffixed `_bw` to disambiguate from the PD canonical `lhill`
  reserved for sigmoidal Emax / Imax (parameter-role distinction, not a
  notation difference).
- `ltbw50` -- log time to 50% of the BW loss. Bare name `tbw50`, time units
  matching the model time unit.

IIV follows the standard pattern (`etalfbw`, `etalhill_bw`, `etaltbw50`); use
`fixed()` if the source paper holds any of these parameters constant.

Founding example: `Oualha_2018_enoxaparin.R`.

## Lactation: milk-to-plasma transfer

Canonical milk-to-plasma **concentration** ratio: **`lcmpr`** (log-transformed),
bare name `cmpr` inside `model()`. Dimensionless.

Use `cmpr` when a lactation popPK paper estimates a milk:plasma concentration
ratio as a model parameter and the milk data are too sparse to support a
separate milk compartment. Breast-milk concentration is then an algebraic
observable rather than an ODE state:

```
cmpr  <- exp(lcmpr + etalcmpr)
Cmilk <- cmpr * Cc
```

The ratio may carry IIV (`etalcmpr`) and covariate effects — a power function of
time postpartum is common during the colostrum period. Pair it with the
`Cmilk` observable and the `propSd_Cmilk` residual (the milk observable follows
the paper-named multi-output exemption in `compartment-names.md`, alongside the
registered `propSd_Ccsf` / `propSd_Cbrain` precedents).

Source-paper aliases that translate without a sidecar: `MPRcon`, `M/P`.

Do **not** encode an AUC-based milk:plasma ratio as a parameter unless the paper
actually parameterises the model on exposure ratio. Papers commonly report an
`MPRauc` as a *derived simulation output*; compute that in the vignette from the
simulated profiles instead.

Distinct from `kp_milk` below: `cmpr` is **fitted** to paired maternal plasma and
milk observations, whereas a `kp` is a predicted or literature-fixed partition
constant feeding a mechanistic milk compartment. The two coexist.

Founding example: `Li_2023_ornidazole.R` (ornidazole into colostrum;
`lcmpr = log(0.58)`, IIV variance 0.327, power effect of time postpartum).

## Tissue partition coefficients

Canonical per-tissue partition-coefficient family: **`lkp_<tissue>`**
(log-transformed), bare names `kp_<tissue>` inside `model()`. Unitless. Each is
the equilibrium ratio of concentration in one named anatomical tissue to
concentration in plasma or blood, driving that tissue's perfusion- or
permeability-limited distribution term:

```
d/dt(liver) <- q_liver * (Cc - liver / vliver / kp_liver) - ...
```

`<tissue>` is the same lowercase anatomical token as the corresponding PBPK
compartment name in `compartment-names.md`, so `lkp_liver` belongs to the
`liver` compartment. Members in use include `adipose`, `bone`, `brain`, `csf`,
`gut`, `heart`, `intestine`, `kidney`, `liver`, `lung`, `milk`, `muscle`,
`other`, `rest`, `skin`, `spleen`, `trachea`. A new anatomical tissue may be
added without a fresh naming sidecar so long as the token matches the canonical
compartment name. Usually fixed from the paper's physicochemical predictions
(Rodgers-Rowland, Poulin-Theil) or measured tissue:plasma ratios, but may be
estimated when the paper fits them.

`kp_milk` / `lkp_milk` is the lactation member — the milk:plasma partition
coefficient for **mechanistic milk-compartment** models. Use it only when the
paper supplies a partition coefficient as a mechanistic input; use `cmpr` above
when the paper estimates a concentration ratio as a popPK parameter.

Three boundaries to respect:

- `kpu<n>` is a *cluster*-indexed unbound partition coefficient shared by several
  tissues, not one named organ; `sf<n>` scales a predicted Kpu rather than
  replacing it. Neither is a member of this family.
- `lkp_f` and `lkp_hb` (`Li_2015_taspoglutide_mbma.R`) are placebo-response rate
  constants for fasting plasma glucose and HbA1c, and `lkpin` / `lkpout` are
  precursor-pool turnover constants. These share the prefix but are **not**
  partition coefficients.
- Bare `kp_*` names that qualify rather than name a tissue (`kp_free`,
  `kp_bound`, `kp_preg`) or index a loop (`kp_i`) are local `model()`
  intermediates, not canonicals.

## Transform prefixes

- `l` — natural log (`lka`, `lcl`, `lfdepot`)
- `logit` — logit (`logitffo`, `logitemax`)
- `probit` — probit
- Any other transform: spell it out as a prefix

Always label every parameter inside `ini()` with units and a short interpretation.

## Fixed parameters

<!-- Canonical home — referenced from SKILL.md Phase 3 and verification-checklist.md A. -->

Wrap a parameter value in `fixed()` whenever the source paper holds it constant rather than estimating it. This applies to **every** parameter type: structural THETAs, allometric exponents, IIV variances and covariances, residual-error magnitudes, covariate-effect coefficients, bioavailability anchors. The `fixed()` wrapper is load-bearing provenance — a downstream user must be able to tell which values are estimated point estimates vs structural assumptions.

Source signals that a parameter is fixed:

- Prose: "fixed during estimation", "fixed at <value>", "held fixed at the literature value", "not estimated", "set to 1 (fixed)".
- Allometric exponents without RSE / SE / CI (canonical 0.75 / 1, especially under "fixed exponents" wording).
- NONMEM `$THETA` / `$OMEGA` / `$SIGMA` with a `FIX` flag.
- Bioavailability `F1 = 1` set as structural anchor.
- Inherited parameters from an upstream publication that the current paper re-uses without re-fitting.

**Do not repeat `fixed()` in the label.** `fixed()` is the machine-readable
statement that the value was not estimated; writing "fixed" in the label as well
says nothing extra, and `checkModelConventions()` raises an **error** for it.
996 labels across 373 models were cleaned of this in 2026-08.

Keep whatever sits around the word, because that is the part `fixed()` cannot
express -- **where** the value came from, or that *you* assumed it rather than
the paper fixing it:

| Write this | Not this |
|---|---|
| `label("Allometric exponent on CL (unitless)")` | `label("Allometric exponent on CL (unitless; fixed)")` |
| `label("Maternal apparent clearance, from Rizk 2015 (CL, L/h)")` | `label("Maternal apparent clearance, fixed from Rizk 2015 (CL, L/h)")` |
| `label("Additive residual SD (mg/L; placeholder)")` | `label("Additive residual SD (mg/L; FIXED placeholder)")` |
| `label("Vascular reflection coefficient (unitless)")` | `label("Vascular reflection coefficient (unitless; fixed at 0.95)")` |

The last row is the general case: if the label only restates the number, drop the
clause -- `fixed(0.95)` already holds it.

Words that **must stay**, because they distinguish *the paper fixed this* from
*the encoder assumed this* and `fixed()` cannot tell them apart (issue #479):
`assumed`, `not published`, `not paper-derived`, `literature value`,
`taken from`, `transferred from`.

"fixed effect", "fixed dose" and "fixed-dose" are exempt -- there "fixed"
modifies something else and is not a claim about this parameter.

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

### Stratum-suffixed parameters

Some papers estimate the *same quantity twice* -- once in each of two or more
sub-populations -- inside a **single joint fit**. The canonical name supplies
exactly one slot for that quantity, so each stratum's estimate takes an explicit
suffix naming the stratum:

```r
lcl_agegt2     <- log(2.59) ; label("Typical clearance for age > 2 years ... (L/h)")
lcl_agele2     <- log(1.98) ; label("Typical clearance for age <= 2 years ... (L/h)")
e_wt_cl_agegt2 <- 0.38      ; label("Power exponent on (WT/12) for CL, age > 2 years (unitless)")
e_wt_cl_agele2 <- 0.739     ; label("Power exponent on (WT/12) for CL, age <= 2 years (unitless)")
```

Rules:

- **Symmetric.** Every stratum carries a suffix; none keeps the bare canonical.
  A bare `lcl` silently meaning "the age > 2 y value" is exactly the ambiguity
  the suffix exists to remove, and papers of this shape report parallel
  estimates rather than a reference plus an offset.
- **The suffix names the stratum, not the paper's symbol.** `_agele2` /
  `_agegt2` (age <= 2 / age > 2), `_center1` / `_center2`, `_ped` / `_adult`.
  Keep it lower-case and readable; a reader must be able to tell which subgroup
  the number belongs to without opening the paper.
- **Only suffix what is actually stratum-specific.** Anything the joint fit
  shares (in the founding example: the volume, the CLcr exponent, the IIV)
  keeps its bare canonical name. Suffixing a shared parameter falsely implies
  the paper estimated it more than once.
- **Applies on top of any canonical form**, including the covariate-effect
  grammar: `e_<cov>_<param>_<stratum>`. This is *not* one of the three-token
  forms above -- the trailing token is a stratum label, not a metabolite, a
  second PK parameter, or a CL component. Disambiguate by what the token names:
  a registered metabolite (`R/conventions.R::registeredMetabolites`) -> the
  metabolite form; a bare PK parameter (`pkBareParams`) -> the shared-exponent
  form; `ss` / `time` / `renal` -> multi-component CL; anything else -> a
  stratum.
- **One model file, not N.** A stratified joint fit is a single model
  (`references/replicate-author-structure.md`); the suffixes are what let one
  file hold every stratum. Split into separate files only when the paper fit
  separate models.

Founding example: `Shen_2024_vancomycin.R` (Shen 2024 vancomycin popPK -- an
age-cutoff model estimating both a typical CL and a body-weight allometric
exponent in each of two age strata within one NONMEM run, OFV 2272.810).
Precedent predating this entry, in files that use the pattern without having
registered it: `Duong_2016_WHIG_T2DM.R` (`b0_s1` / `b0_s23`),
`Frohlich_2018_mRNA_translation.R` (`lkdeg_egfp` / `lkdeg_d2egfp`),
`Friberg_2012_voriconazole.R` (`expSdStdy1` / `expSdStdy2` / `expSdStdy34`).

Note that `checkModelConventions()` does not validate structural-parameter
names, so a clean lint is not evidence that a stratum suffix is well-formed --
this section is the record.

## Endogenous / mechanistic parameters

For endogenous turnover, steady-state, and enzyme-kinetic models (e.g., `igg_kim_2006`, `phenylalanine_charbonneau_2021`), parameters come from mechanism rather than from a CL/V parameterization. Use the names the paper uses; lower-case snake-case by default. Log-transform only positive-constrained parameters that are being *estimated* — not mechanistic constants that the source paper reports as point values.

Recommended patterns:

- `vmax`, `km` — Michaelis–Menten constants for each enzyme / transporter. Name-disambiguate when several coexist: `vmax_pah`, `km_pah`, `vmax_trans`, `km_trans`.
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
- `<FirstAuthor>` is the paper's first-author surname, normalised to plain ASCII CamelCase per the rules below. `<Year>` is the four-digit publication year. `<drug>` is the drug INN in lowercase (with a species or mechanism-class suffix when applicable — see SKILL.md Phase 1 steps 3 and 3a).

### Author-surname normalization (hard default — no sidecar)

<!-- Canonical home — referenced from SKILL.md Phase 1 Step 3b and pre-flight-checklist.md. -->

CamelCase across separators is the **default behaviour**. Apply silently — do NOT raise a stop-and-ask sidecar to confirm the filename form.

| Source author surname | File basename stem | Notes |
|---|---|---|
| `Lohy Das`     | `LohyDas`       | Drop spaces, CamelCase across the drop |
| `Ait-Oudhia`   | `AitOudhia`     | Drop hyphens, CamelCase across the drop |
| `O'Brien`      | `OBrien`        | Drop apostrophes, CamelCase across the drop |
| `Câmara`       | `Camara`        | Transliterate accents |
| `Müller`       | `Muller`        | Transliterate accents |
| `Łukasz`       | `Lukasz`        | Transliterate accented Latin extensions |
| `van Rongen`   | `vanRongen`     | Lowercase particle preserved when author publishes that way |
| `de Winter`    | `deWinter`      | Lowercase particle preserved |
| `Von Bonin`    | `VonBonin`      | Capitalised particle preserved when author publishes that way |
| `Mac Donald`   | `MacDonald`     | Drop space, CamelCase |
| `Câmara De Souza` | `CamaraDeSouza` | Combination: transliterate accent AND drop space AND preserve published particle case |

Rules:

- **Hyphens, spaces, apostrophes**: drop, CamelCase across the drop.
- **Accents and other diacritics**: transliterate to plain ASCII. The filename must be ASCII-only so the function name (= filename minus `.R`) is a valid R identifier without quoting.
- **Lowercase particles** ("van", "de", "von", "ten", "le", "du", "la", "del", "el", etc.): preserve the lowercase first letter when the author publishes that way (`van Rongen` → `vanRongen`, `de Kock` → `deKock`).
- **Capitalised particles**: follow the published author form (`Von Bonin` → `VonBonin`, `De Kock` → `DeKock`).

Precedents in `inst/modeldb/specificDrugs/`: `AitOudhia_2012_canakinumab.R`. The same normalised stem is used for the vignette basename, the `vignette <- "..."` field, the worktree branch name in Phase 2, and the PR title.

### Year-letter collision suffix

When two extractions resolve to the same `<FirstAuthor>_<Year>_<drug>` name (e.g., two same-author/year/drug entries with different scenarios), append a lowercase letter to the year — `Author_2019a_drug.R`, `Author_2019b_drug.R`. Allocate letters in chronological model-development order when known. Never overwrite an existing file silently.

## Time-varying clearance (issue #481)

If clearance depends on time, use the stem that matches the functional form so
the structure can be found by name. `checkModelConventions()` warns when a
clearance expression references `t` / `time` without one of these.

| Form | Stem | Roles |
|---|---|---|
| Sigmoidal in time: `cl <- cl_base * exp(max * t^g / (t50^g + t^g))` | `cl_time_` | `cl_time_max`, `cl_t50`, `cl_time_hill` |
| Exponential decay to a constant: `cl <- cl_exp_inf + cl_exp_component * exp(-k * t)` | `cl_exp_` | `cl_exp_inf`, `cl_exp_component`, `cl_exp_kdes` |

Prefix `l` for the log scale (`lcl_t50`), `eta` for the IIV partner
(`etacl_time_max`), `e_<cov>_` for a covariate effect
(`e_nhl_cl_exp_kdes`).

Do **not** use `emax`, `imax`, `gamma`, `hill` or `t50` BARE for clearance
time-dependence: all of them are also standard PD parameter names, and several
models carry both a PD `emax` and a clearance one. The `cl_` prefix is what
keeps `cl_t50` and `cl_time_hill` distinct from their PD namesakes.

**The symbol the ODE consumes is the total clearance.** Name components so they
cannot be mistaken for it -- `cl_exp_component`, not `cl_time` or `cl_t_now`. A
decaying component can fall by dozens of orders of magnitude over a treatment
course, which is meaningless in isolation but looks like a clearance value.

Periodic (diurnal / circadian) variation is a different structure and keeps its
own names; do not fold it into `cl_time_` or `cl_exp_`.
