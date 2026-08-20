# Canonical parameter names

This file is the authoritative register of structural and paper-mechanistic parameter names used in nlmixr2lib models. Every parameter referenced inside a model's `ini()` and `model()` blocks is expected to match one of the canonical names below (or, for covariate effects, the documented `e_<cov>_<param>` shape — see `R/conventions.R::covEffectPattern`). The register is seeded from the 2026-05-28 naming audit and extended whenever a new paper introduces a parameter that isn't yet registered.

## How to use this register

1. **Before adding a parameter to a new model**, search this file (by canonical name and by source alias) for the concept you need.
2. **If the canonical name exists**, use it exactly. Document any source-paper rename in a code comment or in the model's `description` block.
3. **If the source paper uses an alias listed under an existing canonical name**, prefer the canonical name. Aliases are documented for cross-reference, not as a free pass to introduce the deprecated form in new models.
4. **If the parameter is not in this register at all**, propose a new entry with a canonical name, type, role, source aliases, and example models. Verify with the user before committing. The addition is part of the model's PR.
5. **Do not modify existing model files when you discover a missing entry**; simply register the canonical here. Retrofitting existing models is a separate effort.
6. **Never add a second entry for a name that already has one at the same `Type`.** Extend the existing block instead -- add your source alias and example model to it. A name may appear twice only when the two entries carry *different* `Type` values. This register is resolved in document order, last one wins, so a same-`Type` repeat silently discards the earlier block along with any alias or example recorded only there. `buildModelDb()` fails the build on one.

## Entry schema

Each canonical entry is an H3 heading whose first whitespace-separated token (before the parenthetical) is the canonical name. Required fields:

```yaml
- name: <CANONICAL_NAME>
  type: log-transformed-pk | bare-pk | paper-named-param
  role: <one-sentence description of mechanistic role>
  source_aliases:
    - <ALIAS_NAME> -- used in <model.R>
  example_models:
    - <model.R>
  notes: <free text>
```

The `Type:` field is the routing tag the runtime parser uses to assign the entry to the appropriate static vector:

- `log-transformed-pk` -> `pkParams` (used by `.isPkParam` and the `l<base>` convention check)
- `bare-pk` -> `pkBareParams` (used by `.isPkBareParam`, the bare-counterpart check, and covariate-effect shared-exponent detection)
- `paper-named-param` -> `paperNamedParams` (paper-mechanistic parameters that fall outside the standard `ka`/`cl`/`vc` shape but recur across published models)

## Regex constants (kept in R, not migrated)

The following pattern constants remain hard-coded in `R/conventions.R::.nlmixr2libConventionsStatic` because they are structural regular expressions rather than name lists. They are documented here for cross-reference:

- `covEffectPattern = "^e_[A-Za-z0-9]+(_[A-Za-z0-9]+){1,5}$"` — canonical shape of a covariate-effect parameter name: a leading `e_`, a covariate token, and one to five additional tokens (final token is the affected PK parameter, optional middle tokens carry metabolite / shared-exponent / CL-component context).
- `clComponents = c("ss", "time", "renal", "nonren")` — suffix tokens permitted on multi-component CL parameters (e.g., `cl_renal + cl_nonren`).
- `transformPrefixes = c("l", "logit", "probit")` — accepted parameter-name transform prefixes.
- Within-subject random-effect levels, matched by `R/checkModelConventions.R::.isIOVEtaSuffix`: `etaiov_<param>_<occ>` for inter-occasion variability and `etabvv_<param>_<visit>` for between-visit variability, where `<param>` is an existing fixed-effect parameter (bare, or with an `l` / `ltv` / `lv` prefix) and the trailing integer is the occasion / visit index. The two prefixes are kept separate because a paper can fit both levels at once on different parameters -- `Abdelgawad_2024_linezolid.R` carries five-occasion IOV on `ka` and `mtt` alongside a two-visit BVV on `vmax`, and folding them into one prefix would lose the published random-effect hierarchy. NONMEM `$OMEGA BLOCK(1) SAME` repeats map to `~ fix(<variance>)` on every index after the first.
- `residualError = c("propSd", "addSd", "expSd")` — canonical residual-error SD names. `propSd` / `addSd` are the proportional and additive SDs for the standard `~ prop(...)`, `~ add(...)`, and combined `~ prop(...) + add(...)` forms. `expSd` is the log-scale residual SD used with `~ lnorm(...)`.
- `deprecatedResidualError`, `deprecatedIivPrefixes`, `deprecatedVolumeNames`, `deprecatedVmaxNames`, `deprecatedParentSuffix` — deprecation lists also remain in R.

## Covariate-expression shape tokens

Most published covariate models are multiplicative and need exactly one
coefficient per (covariate, parameter) pair, which `e_<cov>_<param>` covers.
Machine-learned covariate models (symbolic regression, equation learning,
sparse neural networks) instead produce low-order **polynomials and ratios of
polynomials**, so one pair can carry several coefficients that must be told
apart by name. The tokens below are the canonical way to do that. They are
shape tokens on the existing `e_<cov>_<param>` family and on bare
`<param>_<token>` names -- not new PK parameters -- so no `Type:` entry is
required for them, and `covEffectPattern` already admits the `e_` forms.

- `_quad` -- coefficient of a **squared** normalised covariate, appended as the
  final token: `e_wt_vc_quad` is the coefficient of `(WT/ref)^2` on `vc`. The
  unsuffixed `e_wt_vc` remains the linear coefficient. Precedent:
  `e_pdl1_orr_pd1_quad` / `e_pdl1_orr_pd1_lin` in
  `Franzese_2026_pdl1_nsclc_mbma.R`. Prefer leaving the linear term unsuffixed
  over writing `_lin`, because `_lin` already means "linear clearance
  component" elsewhere in the library.
- `_num` / `_den` -- position within a **rational** covariate expression of the
  form `(a*C1 + b*C2 + c) / (d*C1 + e*C2 + f)`, inserted before the final
  parameter token would be for a bare name and appended for an `e_` name:
  `e_age_k12_num` and `e_age_k12_den` are the two age coefficients of the k12
  ratio. Use these **only** when the numerator and denominator have no
  pharmacological identity of their own. When they do -- a saturable / Hill
  function of a covariate, where the numerator is an asymptote and the
  denominator a half-max -- name the sub-parameter instead and attach the
  covariate effect to it (`e_sexf_lmax` vs `e_pnpla3_cg_bmi50` in
  `Oniki_2018_nafld_risk.R`; `crcl50_cl_renal` / `hill_cl_renal` in
  `Ganesan_2023_tebipenem.R`). That idiom is preferred wherever it applies.
- `e_<cov1>_<cov2>_<param>` -- coefficient of a **covariate-by-covariate
  interaction**, covariates in the order they are printed in the source:
  `e_age_wt_vc` multiplies `(AGE/ref) * (WT/ref)` on `vc`. Distinct from
  covariate-by-stratum effects, which put the stratum in the middle token
  (`e_pdl1_orr_pd1`, `e_fed_ktr_tab`).
- `<param>_int` -- **signed additive intercept** of a covariate expression, for
  the case `l<param>` cannot express: a negative intercept, or a rational
  expression with no single intercept (`vc_int`, `k12_num_int`). When the
  intercept is positive and the expression is `intercept + covariate terms`,
  use the canonical log form instead (`lkel`, `lk21`, `lk31` with
  `kel <- exp(lkel) + e_wt_kel * wtn`), as in
  `McLachlan_1996_fluconazole.R` and `Blair_2004_raltitrexed.R`. Precedent for
  the bare form: `b1_int` / `b2_int` in `Sherer_2012_AAA.R`.

Founding example: `Wahlquist_2024_propofol.R` (symbolic-regression network
covariate model on all four inter-compartmental rate constants, the
elimination rate constant and the central volume).

---

## Log-transformed structural PK parameters

The `l<base>` convention denotes a population mean estimated on the log scale (`ini()` carries `l<base> <- log(typical_value)`; `model()` declares `<base> <- exp(l<base> + eta_<base>)`). Use the log-transformed form whenever the parameter is strictly positive and the source paper reports an exponential typical-value form.

### lka (**canonical log-transformed absorption rate constant**)
- **Type:** log-transformed-pk
- **Role:** First-order absorption rate constant from the depot into the central compartment (1 / time).
- **Source aliases:** none.
- **Example models:** widespread across oral popPK extractions.
- **Notes:** Always paired with a `depot` compartment.

### lcl (**canonical log-transformed total clearance**)
- **Type:** log-transformed-pk
- **Role:** Apparent total drug clearance from the central compartment (volume / time).
- **Source aliases:** none.
- **Example models:** universal in popPK extractions with explicit clearance.
- **Notes:** When clearance is decomposed (steady-state / time-varying / renal / non-renal), use the multi-component forms `lcl_ss`, `lcl_time`, `lcl_renal`, `lcl_nonren`.

### lvc (**canonical log-transformed central volume**)
- **Type:** log-transformed-pk
- **Role:** Apparent volume of distribution of the central compartment (volume).
- **Source aliases:** none.
- **Example models:** universal in popPK extractions.
- **Notes:** Replaces the deprecated bare-volume names `v`, `v1`, `lv`, `lv1` (see `deprecatedVolumeNames`).

### lvp (**canonical log-transformed first peripheral volume**)
- **Type:** log-transformed-pk
- **Role:** Apparent volume of the first peripheral compartment in 2- and 3-compartment models (volume).
- **Source aliases:** none.
- **Example models:** any 2-compartment or higher popPK extraction.
- **Notes:** Pairs with `peripheral1` compartment and inter-compartmental clearance `lq`.

### lvp2 (**canonical log-transformed second peripheral volume**)
- **Type:** log-transformed-pk
- **Role:** Apparent volume of the second peripheral compartment in 3-compartment models (volume).
- **Source aliases:** none.
- **Example models:** any 3-compartment popPK extraction.
- **Notes:** Pairs with `peripheral2` compartment and inter-compartmental clearance `lq2`.

### lvp3 (**canonical log-transformed third peripheral volume**)
- **Type:** log-transformed-pk
- **Role:** Apparent volume of the third peripheral compartment in 4-compartment models (volume).
- **Source aliases:** none.
- **Example models:** `Schmitt_2018_vinflunine.R`, `Li_2017_brentuximab.R`, `Weatherley_2009_maraviroc_iv.R`.
- **Notes:** Pairs with `peripheral3` compartment and inter-compartmental clearance `lq3`.

### lvelf (**canonical log-transformed epithelial-lining-fluid compartment volume**)
- **Type:** log-transformed-pk
- **Role:** Apparent volume of the canonical `elf` compartment, used to convert the ELF drug amount to the ELF concentration `Celf <- elf / velf` (volume). Applies both to a plasma-plus-ELF popPK model in which `elf` is a distribution compartment sampled by bronchoalveolar lavage, and to a single-compartment model fitted directly to ELF concentrations, where the value is an apparent volume that also absorbs bioavailability (the paper's `Vc/F`).
- **Source aliases:**
  - `VELF`, `V ELF`, `V_ELF` -- Abouelhassan 2024 Results and Parmar 2023 Table 2.
  - `Vc/F` -- Abouelhassan 2024 Table 2, for the murine one-compartment ELF fit whose only disposition compartment IS the ELF.
- **Example models:** `Abouelhassan_2024_sulbactam_human.R` (VELF = 2.44 L in a plasma + peripheral + ELF model), `Abouelhassan_2024_sulbactam_mouse.R` (apparent ELF Vc/F = 0.5125 L/kg).
- **Notes:** Member of the `lv<compartment>` family that names a volume after the canonical compartment it scales, following the existing `lvcsf` / `lvcns` precedent in `Luu_2017_nusinersen.R`. Distinct from `lvc`: `lvelf` never scales a plasma concentration. `Parmar_2023_spectinamide_1599_mouse_pbpk.R` carries the same quantity as a hardcoded physiologic constant `v_elf` inside `model()` rather than as an `ini()` parameter, which is correct for a PBPK model where the volume is anatomy rather than a fitted parameter.

### lq (**canonical log-transformed first inter-compartmental clearance**)
- **Type:** log-transformed-pk
- **Role:** Inter-compartmental clearance between central and first peripheral compartment (volume / time).
- **Source aliases:** none.
- **Example models:** universal in 2-compartment+ popPK extractions.

### lq2 (**canonical log-transformed second inter-compartmental clearance**)
- **Type:** log-transformed-pk
- **Role:** Inter-compartmental clearance that feeds the second peripheral compartment (volume / time). In the usual parallel (mammillary) layout the flow originates in `central`.
- **Source aliases:** `Q3` (third inter-compartmental clearance as numbered in Tuffal 2023, where `Q1` is not defined and the numbering follows the compartment index `V3` rather than the clearance index) -- used in `Tuffal_2023_avalglucosidase_alfa.R`.
- **Example models:** 3-compartment popPK extractions; `Tuffal_2023_avalglucosidase_alfa.R` for the series-chain case.
- **Notes:** Pairs with `peripheral2` and `lvp2`. **Series ("concatenated") peripheral chains:** some papers arrange the peripheral compartments in series, so `peripheral2` is fed from `peripheral1` rather than from `central`. `lq2` remains the correct canonical in that case -- it is still "the second inter-compartmental clearance, the one that feeds `peripheral2`", and the source compartment is read off the ODEs rather than inferred from the parameter name. Do not introduce a parallel inter-peripheral canonical (`lq_p1p2` or similar) for this shape. Example: `Tuffal_2023_avalglucosidase_alfa.R`, where the paper's fixed `Q3` = 1.87 L/h carries drug one-way from `peripheral1` into `peripheral2`. When such a chain also returns drug directly to `central`, that return flow is `lqpc`, not `lq2`.

### lqpc (**canonical log-transformed back-redistribution clearance**)
- **Type:** log-transformed-pk
- **Role:** One-way back-redistribution clearance returning drug from the terminal peripheral compartment of a series ("concatenated") peripheral chain directly to `central` (volume / time).
- **Source aliases:** `Qpc` (`Q_pc`, peripheral-to-central) -- used in `Tuffal_2023_avalglucosidase_alfa.R`.
- **Example models:** `Tuffal_2023_avalglucosidase_alfa.R` (founding example; Tuffal 2023 Table 2 `Qpc` = 0.0206 L/h).
- **Notes:** Use when a paper closes the peripheral chain into a cycle -- `central` <-> `peripheral1` -> `peripheral2` -> `central` -- rather than letting each peripheral compartment exchange bidirectionally with `central`. The flow is unidirectional by construction: it appears as `+ qpc * Cp2` in `d/dt(central)` and `- qpc * Cp2` in the terminal peripheral compartment, with **no matching reverse term**. The directionality is load-bearing, not incidental: Tuffal 2023 explicitly tested a bidirectional variant (Supplemental Digital Content 1 Table 2, rows "Q3 one way" vs "Q3 two ways") and rejected it with the comment "Markedly increased OFV". A slow `qpc` returning a marginal drug fraction is what generates a late secondary concentration rebound / second kinetic sequence. Distinct from `lq` / `lq2` (bidirectional `central` <-> peripheral exchange), from `lq3` (`central` <-> `peripheral3` in 4-compartment models, which pairs with `peripheral3` / `lvp3` -- neither of which exists in a 3-compartment series chain), and from `qp` (inter-compartmental clearance of a **target species** in TMDD models, not of the drug). Bare counterpart is `qpc`.

### lq3 (**canonical log-transformed third inter-compartmental clearance**)
- **Type:** log-transformed-pk
- **Role:** Inter-compartmental clearance between central and third peripheral compartment (volume / time) in 4-compartment models.
- **Source aliases:** none.
- **Example models:** `Schmitt_2018_vinflunine.R`, `Li_2017_brentuximab.R`, `Weatherley_2009_maraviroc_iv.R`.
- **Notes:** Pairs with `peripheral3` and `lvp3`.

### lq_p1_p2 (**canonical log-transformed one-way peripheral1-to-peripheral2 inter-compartmental clearance**)
- **Type:** log-transformed-pk
- **Role:** One-way (unidirectional) inter-compartmental clearance carrying drug from `peripheral1` into `peripheral2` (volume / time), with no return flow along the same path.
- **Source aliases:**
  - `Q3` -- used in `Tiraboschi_2023_avalglucosidaseAlfa.R` for the one-way V2 -> V3 clearance.
- **Example models:** `Tiraboschi_2023_avalglucosidaseAlfa.R` (founding example).
- **Notes:** Deliberately NOT `lq2`. `lq` / `lq2` / `lq3` are mamillary canonicals -- bidirectional exchange between `central` and the first / second / third peripheral compartment. `lq_p1_p2` names a flow that neither starts nor ends at `central` and that runs in one direction only, so reusing `lq2` would overload a bidirectional canonical with a directional concept. Directional naming follows the register's existing practice of giving one-way flows their own names when direction is load-bearing (`kbm` biliary, `lcl_csf` one-way CSF-to-plasma). Pairs with `lq_p2_c` to close a catenary / cyclic topology `central <-> peripheral1 -> peripheral2 -> central`. The source-to-destination token pattern `_<from>_<to>` (`p1`, `p2`, `c` for `peripheral1`, `peripheral2`, `central`) is available to future catenary models.

### lq_p2_c (**canonical log-transformed one-way peripheral2-to-central back-redistribution clearance**)
- **Type:** log-transformed-pk
- **Role:** One-way (unidirectional) inter-compartmental clearance returning drug from `peripheral2` to `central` (volume / time) -- the "back-redistribution" limb that closes a catenary / cyclic distribution loop.
- **Source aliases:**
  - `Qpc` / `QPC` -- used in `Tiraboschi_2023_avalglucosidaseAlfa.R` for the low one-way V3 -> V1 clearance.
- **Example models:** `Tiraboschi_2023_avalglucosidaseAlfa.R` (founding example).
- **Notes:** Not an alias of `lq2` (bidirectional `central` <-> `peripheral2`), of `qp` (target inter-compartmental clearance), or of `lqbile` (one-way plasma-to-gut biliary transfer). Typically much smaller than the forward limb, so it governs the terminal phase: in the founding example Qpc = 0.0157 L/h against Q3 = 1.87 L/h. Pairs with `lq_p1_p2`; see that entry for the `_<from>_<to>` token pattern.

### lfdepot (**canonical log-transformed depot fraction**)
- **Type:** log-transformed-pk
- **Role:** Log-scale parameter for the fraction of dose absorbed via the depot (e.g., parallel-absorption oral models with a fast/slow depot split).
- **Source aliases:** none.
- **Example models:** parallel-absorption oral PK models.

### lvmax (**canonical log-transformed Michaelis-Menten Vmax**)
- **Type:** log-transformed-pk
- **Role:** Log-scale maximum metabolic rate for saturable (Michaelis-Menten) elimination (amount / time).
- **Source aliases:** none.
- **Example models:** TMDD models, saturable-elimination popPK extractions.
- **Notes:** Replaces the deprecated `vm` / `lvm` names (see `deprecatedVmaxNames`).

### lbmax (**canonical log-transformed maximum binding capacity**)
- **Type:** log-transformed-pk
- **Role:** Log-scale maximum binding capacity of a saturable binding site or tissue -- the ceiling on bound drug (or bound radioligand), expressed in the units of the state it saturates (amount or concentration). Distinct from `lvmax`, which is a maximum *rate*: `lbmax` is a maximum *amount*.
- **Source aliases:**
  - `BMAX` / `Bmax` -- upper-case source-paper forms; the canonical bare name is `bmax`.
- **Example models:** `Siebinga_2023_lu177psma617.R` (Bmax 40.4 MBq for saturable [177Lu]Lu-PSMA-617 binding in the salivary glands, Table 2), plus the wider `lbmax` / `bmax` population in `inst/modeldb/` (`GonzalezSales_2024_imetelstat.R`, `Nielsen_2011_*.R`, `Svensson_2016_rifampicin.R`, `Sikma_2020_tacrolimus_unbound_plasma.R`, `deWinter_2009_mycophenolic_acid.R`, ...).
- **Notes:** Registered 2026-07-30 to give the long-used `lbmax` / `bmax` names an explicit register entry. Analyte- or site-suffixed variants (`lbmax_c`, `lbmax_p`, `lbmax_rbc`) are permitted on the standard suffix pattern when a model carries more than one binding pool.

### ltmax_abs (**canonical log-transformed saturable absorption Vmax**)
- **Type:** log-transformed-pk
- **Role:** Log-scale maximum rate of saturable Michaelis-Menten ABSORPTION from the depot / absorption compartment into central (amount / time). Distinct from `lvmax`, which is MM ELIMINATION from central.
- **Source aliases:**
  - `Tmax` -- used in `Jansson_2008_eflornithine_rat.R` for the saturable-absorption Vmax. The paper symbol Tmax is not a time-to-Cmax descriptor but the maximum absorption rate.
- **Example models:** `Jansson_2008_eflornithine_rat.R`.
- **Notes:** Paired with `lkt_abs` (half-saturation amount). The `_abs` suffix disambiguates from `lvmax` (MM elimination).

### lkt_abs (**canonical log-transformed saturable absorption half-saturation amount**)
- **Type:** log-transformed-pk
- **Role:** Log-scale amount of drug in the depot / absorption compartment at which the saturable absorption rate equals half its maximum (amount). Functionally a `km` but expressed in amount-in-compartment units rather than concentration.
- **Source aliases:**
  - `Kt` -- used in `Jansson_2008_eflornithine_rat.R`.
- **Example models:** `Jansson_2008_eflornithine_rat.R`.
- **Notes:** Paired with `ltmax_abs`. The amount-in-compartment formulation differs from a concentration `km` because the absorption compartment may not have a well-defined volume.

### lcl_ss (**canonical log-transformed steady-state clearance arm**)
- **Type:** log-transformed-pk
- **Role:** Steady-state component of a time-varying clearance decomposition `CL_total(t) = CL_ss + CL_time(t)`.
- **Source aliases:** none.
- **Example models:** time-varying-clearance popPK extractions (e.g., mAb induction-decay disposition).

### lcl_time (**canonical log-transformed time-varying clearance arm**)
- **Type:** log-transformed-pk
- **Role:** Time-decay component of a time-varying clearance decomposition (paired with `lcl_ss`).
- **Source aliases:** none.
- **Example models:** time-varying-clearance popPK extractions.

### lcl_renal (**canonical log-transformed renal clearance arm**)
- **Type:** log-transformed-pk
- **Role:** Renal (glomerular-filtration / tubular-secretion) component of an additive renal + non-renal clearance decomposition `CL_total = CL_renal + CL_nonren`.
- **Source aliases:** none.
- **Example models:** `Jonckheere_2019_cefepime.R`.

### lcl_nonren (**canonical log-transformed non-renal clearance arm**)
- **Type:** log-transformed-pk
- **Role:** Non-renal (hepatic / metabolic / extra-renal) component of an additive renal + non-renal clearance decomposition.
- **Source aliases:** none.
- **Example models:** `Jonckheere_2019_cefepime.R`.

### lcl_hemodialysis (**canonical log-transformed dialysis-active clearance arm**)
- **Type:** log-transformed-pk
- **Role:** Extracorporeal renal-replacement-therapy clearance arm contributed by an active hemodialysis (or hemofiltration / hemodiafiltration) session. Paired with the canonical body-CL parameter `lcl` (and the per-time-point `HEMODIALYSIS` covariate) to express total apparent clearance as `cl_total <- cl + HEMODIALYSIS * cl_hemodialysis`. The `cl_hemodialysis` arm is added to the body baseline only while a session is running.
- **Source aliases:**
  - `CL_HD` -- Veinstein 2013 paper notation (the hemodialysis-arm clearance estimated as a primary structural THETA with its own IIV; gated by the per-time-point hemodialysis-active indicator).
  - `CLHD`, `CL_HF`, `CL_HDF`, `CL_dialysis` -- variant abbreviations used in adjacent ESRD / CRRT popPK literature.
- **Example models:** `Veinstein_2013_gentamicin.R` (primary `ini()` parameter with IIV; the dialysis arm is estimated as a structural THETA gated by `HEMODIALYSIS`).
- **Notes:** Distinct from `lcl_renal` (= residual renal CL, an intrinsic-body component) and `lcl_nonren` (= non-renal intrinsic-body CL). `Liesenfeld_2013_dabigatran.R` derives an equivalent dialysis-arm quantity from the Michaels equation (a function of blood flow rate, dialysate flow rate, and a hemodialyzer mass-transfer-area coefficient) as a derived `cl_dialysis` expression in `model()` rather than a primary `ini()` parameter, so its file does not include `lcl_hemodialysis`. Covariate-effect names on this arm follow the standard shape `e_<cov>_cl_hemodialysis`.

### lcl_met (**canonical log-transformed metabolic-formation clearance arm**)
- **Type:** log-transformed-pk
- **Role:** Metabolite-formation component of an additive metabolic + non-metabolic clearance decomposition `CL_parent_total = CL_met + CL_nonmet` in parent + metabolite popPK models where the source paper estimates both elimination arms of the parent separately (the formation flux into the metabolite compartment, and the non-formation loss out of the system). Distinct from `lqm` (paper-specific NONMEM `Qm` formation-rate symbol used by Kunarajah 2017) and from `lcl_<metab>` (which is the metabolite's own elimination clearance, not the parent-to-metabolite formation rate). Used in the parent's `d/dt(central)` total-loss term `-(cl_met + cl_nonmet)/vc * central` and the metabolite's `d/dt(central_<metab>)` input term `+cl_met/vc * central`. Parallel to `lcl_renal` / `lcl_nonren`.
- **Source aliases:**
  - `CL_met` -- Lehr 2010 paper notation (the formation clearance of M1 from tesofensine).
- **Example models:** `Lehr_2010_tesofensine.R` (paper Table I: CL_met/F = 0.416 L/h, no IIV; formation of M1 from tesofensine).

### lcl_nonmet (**canonical log-transformed non-metabolic-formation clearance arm**)
- **Type:** log-transformed-pk
- **Role:** Non-formation (renal / non-metabolite-forming hepatic / other) component of an additive metabolic + non-metabolic clearance decomposition (paired with `lcl_met`). Used in the parent's `d/dt(central)` total-loss term `-(cl_met + cl_nonmet)/vc * central` to account for parent elimination not captured by the metabolite-formation flux.
- **Source aliases:**
  - `CL_non-met` -- Lehr 2010 paper notation.
- **Example models:** `Lehr_2010_tesofensine.R` (paper Table I: CL_non-met/F = 1.31 L/h, IIV 42.2% CV; carries the parent CL IIV).

### lcl_form (**canonical log-transformed metabolite formation clearance**)
- **Type:** log-transformed-pk
- **Role:** Formation (parent-to-metabolite) clearance parameter for parent + metabolite popPK models in which each metabolite's formation flux is estimated as a distinct clearance rather than as a fraction of the parent's total CL. Applied with a metabolite suffix as `lcl_form_<metab>` (e.g., `lcl_form_m3g`, `lcl_form_m6g`, `lcl_form_dnef`). The apparent-clearance form absorbs the fraction metabolised and the metabolite bioavailability into a single positive coefficient, so the metabolite ODE is simply `d/dt(central_<metab>) <- cl_form_<metab> * central / vc - cl_<metab> * central_<metab> / vc_<metab>`. Distinct from `lcl_met`: `lcl_met` presupposes a `CL_parent_total = CL_met + CL_nonmet` decomposition of the parent's total elimination (single-metabolite systems), whereas `lcl_form_<metab>` is used when the parent's total CL (`lcl`) is estimated as a separate parameter or when multiple parallel formation clearances co-exist (multi-metabolite systems).
- **Source aliases:**
  - `Qm` -- paper-named formation-rate symbol (Kunarajah 2017 doxorubicin), when the source parameterises the metabolite input as a Q-analogue in mass / volume units.
  - `K13` -- paper-named apparent clearance of metabolisation (Djerada 2014 nefopam Figure 1B; nefopam -> desmethyl-nefopam).
- **Example models:** `Knibbe_2009_morphine.R` (PNA-stratified `lcl_form_m3g_le10` / `lcl_form_m3g_gt10` and `lcl_form_m6g_le10` / `lcl_form_m6g_gt10` for the morphine -> M3G and morphine -> M6G glucuronidation arms), `Franken_2015_morphine.R` (`fm_m3g * cl`, `fm_m6g * cl` derived within `model()` from fixed fractions), `deHoogd_2017_morphine.R`, `Hennig_2015_rifabutin.R`, `Djerada_2014_nefopam.R` (`lcl_form_dnef` = log(K13); the apparent metabolic clearance of nefopam to desmethyl-nefopam).
- **Notes:** The `lcl_form_<metab>` pattern predates the parameter-names.md register (in use since at least Knibbe 2009); formalised on 2026-06-21 alongside the Djerada 2014 nefopam extraction so the convention is discoverable to future extractions. The `lcl_form_<metab>` form does not require a companion `lcl_nonform_<metab>` parameter (mass balance is enforced only when the parent's total CL is independently estimated); use `lcl_met` / `lcl_nonmet` in place of `lcl_form` when the source paper explicitly decomposes total parent CL into formation + non-formation arms with mass balance.

### lcl_2b6 (**canonical log-transformed CYP2B6-mediated clearance arm**)
- **Type:** log-transformed-pk
- **Role:** CYP2B6-mediated metabolic clearance arm of a parent + metabolite popPK decomposition. When the parent drug induces CYP2B6, the time-varying enzyme amount in `enzyme_2b6` multiplies this arm to produce the dynamic CL contribution. Combined with metabolite suffixes (`lcl_2b6_8oh`) to express the same isoenzyme-mediated arm on a metabolite.
- **Source aliases:**
  - `CL-EFV,2B6` -- Heathman 2024 paper notation (parent EFV cleared via CYP2B6 forming 8-OH EFV).
- **Example models:** `Heathman_2024_efavirenz.R`.
- **Notes:** Founding example: `Heathman_2024_efavirenz.R` (parent CL-EFV,2B6 = 3.64 L/h; modulated dynamically by `enzyme_2b6(t)` and by CYP2B6 IM / PM phenotype reductions).

### lcl_2a6 (**canonical log-transformed CYP2A6-mediated clearance arm**)
- **Type:** log-transformed-pk
- **Role:** CYP2A6-mediated metabolic clearance arm of a parent + metabolite popPK decomposition. When the parent drug induces CYP2A6, the time-varying enzyme amount in `enzyme_2a6` multiplies this arm to produce the dynamic CL contribution.
- **Source aliases:**
  - `CL-EFV,2A6` -- Heathman 2024 paper notation (parent EFV cleared via CYP2A6 forming 7-OH EFV).
- **Example models:** `Heathman_2024_efavirenz.R`.
- **Notes:** Founding example: `Heathman_2024_efavirenz.R` (parent CL-EFV,2A6 = 0.0947 L/h; modulated dynamically by `enzyme_2a6(t)`).

### lcl_ugt (**canonical log-transformed UGT-mediated clearance arm**)
- **Type:** log-transformed-pk
- **Role:** UGT (uridine-diphosphate glucuronosyltransferase, typically UGT2B7) -mediated metabolic clearance arm of a parent + metabolite popPK decomposition. Not autoinduced by the drugs in the founding example (EFV does not induce UGT2B7 in Heathman 2024); the arm is constant in time.
- **Source aliases:**
  - `CL-EFV,UGT` -- Heathman 2024 paper notation (parent EFV cleared via UGT2B7).
  - `CL-8OH,UGT` -- Heathman 2024 paper notation when applied to the 8-OH metabolite (in that role the canonical name is `lcl_ugt_8oh`).
- **Example models:** `Heathman_2024_efavirenz.R`.
- **Notes:** Founding example: `Heathman_2024_efavirenz.R` (parent CL-EFV,UGT = 0.0504 L/h; constant in time).

### ld1 (**canonical log-transformed zero-order absorption duration**)
- **Type:** log-transformed-pk
- **Role:** Duration of the zero-order input rate into the depot (or central) compartment for sequential zero+first-order or pure zero-order absorption models. Applied via `dur(depot) <- d1` or `dur(central) <- d1` inside `model()`; the dose mass is delivered uniformly over `[0, d1]`.
- **Source aliases:**
  - `D1` -- common NONMEM / paper notation (Aouri 2017, Goggin 2004, Heathman 2024, Gupta 2016).
  - `lduration` -- Gupta 2016 paper notation; an equivalent log-transformed entry.
- **Example models:** `Aouri_2017_rilpivirine.R` (pure zero-order absorption directly into central; D1 = 4 h), `Goggin_2004_emfilermin.R` (zero-order SC absorption, D1 = 0.84 h), `Heathman_2024_efavirenz.R` (sequential zero-order then first-order; D1 = 1.74 h followed by KA = 0.165/h).

### lkel (**canonical log-transformed elimination rate constant (K-PD)**)
- **Type:** log-transformed-pk
- **Role:** First-order elimination rate constant used when no explicit `vc` is estimated (K-PD or single-rate-constant elimination form).
- **Source aliases:**
  - `lke` -- legacy name; replaced 2026-05-28 by the naming audit.
  - `lkde` -- paper-named (Mazzocco 2015 / Shoji 2017 KDE) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
  - `lkp` -- paper-named (van Hasselt 2015 KP) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
- **Example models:** `Mazzocco_2015_temozolomide.R`, `Shoji_2017_fosdagrocorat_oc.R`, `Shoji_2017_fosdagrocorat_p1np.R`, `vanHasselt_2015_eribulin.R`.
- **Notes:** Canonical `lkel` adopted 2026-05-28 per the naming audit.

### lkel_exp_kdes (**canonical log-transformed decay-rate constant of a time-varying elimination rate constant**)
- **Type:** log-transformed-pk
- **Role:** Log-scale rate constant governing how fast a time-varying elimination rate constant relaxes from its value at t = 0 toward its asymptote, in the exponential form `kel(t) = kel * (1 + kel_exp_famp * (1 - exp(-kel_exp_kdes * t)))` (1 / time). `log(2) / kel_exp_kdes` is the half-life of the change.
- **Source aliases:**
  - `ln(kappa_k)` -- used in `Schreib_2024_busulfan.R` (paper `theta_kappa_k`).
- **Example models:** `Schreib_2024_busulfan.R` (`lkel_exp_kdes = -2.965`, giving 0.0516 1/h and a 13.4 h half-life for the decline in busulfan `k` and `CL` over a conditioning course).
- **Notes:** The `kel` counterpart of the registered `cl_exp_kdes` role in the time-varying-clearance family; reuses the `_kdes` role token deliberately so both families are found by the same stem. Ratified on 2026-08-05 (task `oare_PMC11154452` sidecar question q1, answer A).

### ltlag (**canonical log-transformed absorption lag time**)
- **Type:** log-transformed-pk
- **Role:** Log-scale absorption lag time before drug enters the depot (time).
- **Source aliases:**
  - `lalag` -- legacy.
  - `llag` -- legacy.
- **Example models:** delayed-absorption oral PK models.
- **Notes:** Replaces the legacy `lalag` / `llag` forms per the 2026-05-28 naming audit.

### ltacro (**canonical log-transformed acrophase**)
- **Type:** log-transformed-pk
- **Role:** Log-scale peak-time (acrophase) parameter for circadian-IDR templates (time of day at which the rhythm peaks).
- **Source aliases:**
  - `ltz` -- legacy name for the same circadian peak-time parameter.
- **Example models:** circadian indirect-response templates.
- **Notes:** Replaces the legacy `ltz` form. Semantically distinct from a lag-time even though both denote a time-shift parameter (operator clarification, 2026-05-28).

### lclin (**canonical log-transformed plasma-to-tissue influx clearance**)
- **Type:** log-transformed-pk
- **Role:** Influx clearance from plasma central to a tissue extracellular-fluid compartment (volume / time). Drives unbound-drug delivery into brain / tumor ECF in physiological CNS-distribution popPK models.
- **Source aliases:** none.
- **Example models:** `Campagne_2019_cyclophosphamide_mouse.R` (CNS penetration popPK).
- **Notes:** Paired with `lclef` (efflux) and an `ecf` compartment.

### lclef (**canonical log-transformed tissue-to-plasma efflux clearance**)
- **Type:** log-transformed-pk
- **Role:** Efflux clearance from a tissue extracellular-fluid compartment back to plasma central (volume / time). Paired with `lclin`.
- **Source aliases:** none.
- **Example models:** `Campagne_2019_cyclophosphamide_mouse.R`.

### lkamax (**canonical log-transformed Weibull-absorption asymptotic maximum rate constant**)
- **Type:** log-transformed-pk
- **Role:** Log of the maximum / asymptotic first-order absorption rate constant in a Piotrovskij-style time-dependent (Weibull) absorption model `ka(t) = kamax * (1 - exp(-(ra * tad)^gam1))` (1 / time). The bare counterpart inside `model()` is `kamax`.
- **Source aliases:**
  - `KAMAX` -- NONMEM `$THETA` convention used in NONMEM 7-era Weibull-absorption control streams.
- **Example models:** `Desai_2016_isavuconazole.R` (founding example).
- **Notes:** Distinct from `lka` (the single first-order rate constant of the simple zero-shape absorption model). Used together with `lra` and `lgam1` to parameterize Weibull / Piotrovskij time-dependent absorption. Ratified canonically alongside the Desai 2016 isavuconazole extraction.

### lra (**canonical log-transformed Weibull-absorption rate-scaling parameter**)
- **Type:** log-transformed-pk
- **Role:** Log of the rate-scaling parameter inside a Piotrovskij-style Weibull time-dependent ka (1 / time). The product `(ra * tad)^gam1` drives the time-dependence of ka; larger `ra` shifts the ka rise earlier. The bare counterpart inside `model()` is `ra`.
- **Source aliases:**
  - `RA` -- NONMEM `$THETA` convention used in Weibull-absorption control streams.
- **Example models:** `Desai_2016_isavuconazole.R` (founding example).
- **Notes:** Distinct from `lka` (simple first-order absorption) and from any infusion-rate parameter. Used together with `lkamax` and `lgam1`. Ratified canonically alongside the Desai 2016 isavuconazole extraction.

### lgam1 (**canonical log-transformed Weibull-absorption shape parameter**)
- **Type:** log-transformed-pk
- **Role:** Log of the unitless Weibull shape (sigmoidicity) parameter inside a Piotrovskij-style time-dependent ka. Larger values make the ka rise more abruptly; `gam1 = 1` reduces the Weibull form to a simple saturating exponential. The bare counterpart inside `model()` is `gam1`.
- **Source aliases:**
  - `GAM1` / `GAMMA1` -- NONMEM `$THETA` convention used in Weibull-absorption control streams.
- **Example models:** `Desai_2016_isavuconazole.R` (founding example).
- **Notes:** Distinct from `lhill` (sigmoidal Emax / Imax exponent) and from `lgamma` (Friberg myelosuppression feedback / TGI growth exponents). The `gam1` suffix follows the NONMEM convention for Weibull-absorption sigmoidicity. Ratified canonically alongside the Desai 2016 isavuconazole extraction.

### lbeta_cl (**canonical log-transformed exponential-nonlinear-CL slope**)
- **Type:** log-transformed-pk
- **Role:** Log of the concentration-slope coefficient in an exponential-nonlinear-clearance function of the form `CL(C) = CL_0 * exp(beta_cl * C)`, where `CL_0` is the linear-scale clearance at C = 0 (encoded as the standard `lcl` parameter) and `C` is the observed drug concentration in the same units as the paper's Table. Used when the source paper describes a monotonically-increasing clearance-vs-concentration relationship arising from a saturable protein-binding buffer (e.g., FVIII / von Willebrand factor complex saturation at supraphysiological rFVIII doses). Larger `beta_cl` gives a steeper rise in CL with C; `beta_cl = 0` recovers linear elimination. The bare counterpart inside `model()` is `beta_cl` (units of 1 / concentration).
- **Source aliases:**
  - `β` -- Greek letter used in Larsen 2018 Table 3 and Eq. 1.
- **Example models:** `Larsen_2018_factorviii_rat.R`, `Larsen_2018_factorviii_monkey.R` (founding examples; rat `beta_cl = 0.162 mL/IU`, monkey `beta_cl = 0.0355 mL/IU`, both fitting the exponential-nonlinear-CL form on total FVIII activity).
- **Notes:** Distinct from `lhill` (sigmoidal Emax / Imax exponent), `lgamma` (Friberg / TGI kinetic exponents), and `lvmax` / `km` (Michaelis-Menten saturable elimination — a different mathematical form). Choose `lbeta_cl` when the source paper explicitly parameterises CL as an exponential function of concentration; choose `lvmax` / `km` when the paper fits Michaelis-Menten saturation instead. The `_cl` suffix marks that the slope acts on clearance; parallel `beta_<param>` names are permitted for other parameters when a paper extends the exponential-nonlinear form to `V` or `Q`.

### lcmpr (**canonical log-transformed milk-to-plasma concentration ratio**)
- **Type:** log-transformed-pk
- **Role:** Log of the ratio of drug concentration in breast milk to the concentration in maternal plasma, used to derive a milk observable algebraically from the plasma compartment (`Cmilk <- cmpr * Cc`) when the source data are too sparse to support a separate milk compartment with its own transfer clearances. Unitless. May be estimated with IIV (`etalcmpr`) and may carry covariate effects, e.g. a power function of time postpartum during the colostrum period. The bare counterpart inside `model()` is `cmpr`.
- **Source aliases:**
  - `MPRcon` -- Li 2023 ornidazole notation for the milk-to-plasma concentration ratio.
  - `MP` -- common lactation-literature shorthand for the milk:plasma ratio.
- **Example models:** `Li_2023_ornidazole.R` (founding example; `lcmpr = log(0.58)`, RSE 8.63%, IIV variance 0.327, power effect of time postpartum with exponent 1.37 centred on a back-solved median postpartum sampling time of 54 h).
- **Notes:** Distinct from the `lkp_<tissue>` partition-coefficient family (including `lkp_milk`) in that `cmpr` is an **estimated popPK parameter** fitted to paired maternal plasma and milk observations, whereas a `kp` is a physicochemically predicted or literature-fixed partition constant used inside a mechanistic distribution model. Use `cmpr` when the paper estimates a concentration ratio as a model parameter; use `kp_milk` when the paper supplies a milk:plasma partition coefficient as a mechanistic input to a milk compartment. Use an AUC-based ratio name only if the paper actually parameterises the model on an exposure ratio -- Li 2023 reports `MPRauc` as a derived simulation **output**, not a model parameter, and it is deliberately not encoded. Pair with the `Cmilk` observable and the `propSd_Cmilk` residual. Ratified 2026-08-05 alongside the Li 2023 ornidazole extraction, jointly with the `lkp_<tissue>` family registration below.

### lkp_adipose, lkp_brain, lkp_intestine, lkp_kidney, lkp_liver, lkp_lung, lkp_milk, lkp_muscle, lkp_other, lkp_rest, lkp_trachea (**canonical log-transformed tissue-to-plasma partition coefficients**)
- **Type:** log-transformed-pk
- **Role:** Log of the tissue-to-plasma (or tissue-to-blood) partition coefficient Kp for one named anatomical tissue, i.e. the equilibrium ratio of tissue concentration to plasma (or blood) concentration. Unitless. The family shape is `lkp_<tissue>`, where `<tissue>` is the same lowercase anatomical token used for the corresponding PBPK compartment name, so `lkp_liver` is the partition coefficient of the `liver` compartment. Each drives the perfusion- or permeability-limited distribution term for its tissue (e.g. `d/dt(liver) <- q_liver * (Cc - liver / vliver / kp_liver) - ...`). The bare counterparts inside `model()` are `kp_<tissue>`, registered below. Usually fixed from the source paper's physicochemical predictions (Rodgers-Rowland, Poulin-Theil) or measured tissue:plasma ratios, but may be estimated when the paper fits them.
- **Source aliases:**
  - `Kp` -- near-universal PBPK notation, subscripted by tissue in the source tables.
  - `Kt` / `Ptp` -- tissue:plasma partition notation used by some perfusion-limited PBPK papers.
- **Example models:** `Levitt_2005_propofol_pbpk.R` (`lkp_adipose`, `lkp_brain`, `lkp_liver`, `lkp_intestine`, `lkp_rest`), `Mi_2023_cefquinome_pbpk.R` (`lkp_liver`, `lkp_kidney`, `lkp_lung`, `lkp_muscle`, `lkp_other`), `Kang_2023_artesunate_hamster_pbpk.R` and `Kang_2023_pyronaridine_hamster_pbpk.R` (`lkp_lung`, `lkp_trachea`, `lkp_other`, plus metabolite-suffixed `lkp_lung_dihydroart` etc.).
- **Notes:** This entry formally registers a family that had been de-facto in use across the PBPK model files listed above since before the register existed; the members enumerated in the heading are those currently used, and a new anatomical tissue may be added to the heading without a fresh naming sidecar so long as `<tissue>` matches the canonical compartment name in `compartment-names.md`. `lkp_milk` is the lactation member: the milk:plasma partition coefficient for mechanistic milk-compartment models. Distinct from `lkpu<n>`, which is a *cluster*-indexed unbound partition coefficient shared by several tissues rather than one named organ, and from `sf<n>`, which scales a predicted Kpu rather than replacing it. Distinct from `lcmpr` above, which is an estimated popPK milk:plasma **concentration ratio** rather than a mechanistic partition constant. Two existing names collide on prefix but are NOT members of this family and must not be read as partition coefficients: `lkp_f` and `lkp_hb` in `Li_2015_taspoglutide_mbma.R` are placebo-response rate constants for fasting plasma glucose and HbA1c, and `lkpin` / `lkpout` are precursor-pool turnover rate constants. Family registered 2026-08-05 alongside the Li 2023 ornidazole extraction.

---

## Bare structural PK parameters

The bare counterparts of the log-transformed parameters above. Used when the source paper estimates the parameter directly on the linear scale, or when the parameter appears in the `model()` block as the exponentiated form `<base> <- exp(l<base> + eta_<base>)`.

### ka (**canonical bare absorption rate constant**)
- **Type:** bare-pk
- **Role:** First-order absorption rate constant from depot into central (1 / time).
- **Source aliases:** none.
- **Example models:** universal in oral popPK extractions.

### cl (**canonical bare total clearance**)
- **Type:** bare-pk
- **Role:** Apparent total clearance from central (volume / time).
- **Source aliases:** none.
- **Example models:** universal in popPK extractions.

### vc (**canonical bare central volume**)
- **Type:** bare-pk
- **Role:** Apparent volume of central compartment (volume).
- **Source aliases:** none.
- **Example models:** universal in popPK extractions.
- **Notes:** Replaces the deprecated bare-volume names `v`, `v1` (see `deprecatedVolumeNames`).

### vp (**canonical bare first peripheral volume**)
- **Type:** bare-pk
- **Role:** Apparent first-peripheral volume (volume).
- **Source aliases:** none.
- **Example models:** universal in 2-compartment+ popPK extractions.

### vp2 (**canonical bare second peripheral volume**)
- **Type:** bare-pk
- **Role:** Apparent second-peripheral volume (volume).
- **Source aliases:** none.
- **Example models:** universal in 3-compartment popPK extractions.

### vp3 (**canonical bare third peripheral volume**)
- **Type:** bare-pk
- **Role:** Apparent third-peripheral volume in 4-compartment popPK models (volume).
- **Source aliases:** none.
- **Example models:** `Schmitt_2018_vinflunine.R`, `Li_2017_brentuximab.R`, `Weatherley_2009_maraviroc_iv.R`.

### q (**canonical bare first inter-compartmental clearance**)
- **Type:** bare-pk
- **Role:** Inter-compartmental clearance between central and `peripheral1` (volume / time).
- **Source aliases:** none.
- **Example models:** universal in 2-compartment+ popPK extractions.

### q2 (**canonical bare second inter-compartmental clearance**)
- **Type:** bare-pk
- **Role:** Inter-compartmental clearance that feeds `peripheral2` (volume / time); from `central` in the usual parallel layout, or from `peripheral1` in a series ("concatenated") chain.
- **Source aliases:** `Q3` -- used in `Tuffal_2023_avalglucosidase_alfa.R` (see `lq2` notes on the paper's compartment-indexed numbering).
- **Example models:** universal in 3-compartment popPK extractions.
- **Notes:** See `lq2` for the series-chain rule. Pairs with `qpc` when the chain returns drug directly to `central`.

### qpc (**canonical bare back-redistribution clearance**)
- **Type:** bare-pk
- **Role:** One-way back-redistribution clearance from the terminal peripheral compartment of a series ("concatenated") peripheral chain into `central` (volume / time).
- **Source aliases:** `Qpc` -- used in `Tuffal_2023_avalglucosidase_alfa.R`.
- **Example models:** `Tuffal_2023_avalglucosidase_alfa.R` (founding example).
- **Notes:** Bare counterpart of `lqpc`; see that entry for the full role, the unidirectionality requirement, and the distinctions from `q`, `q2`, `q3`, and `qp`.

### q3 (**canonical bare third inter-compartmental clearance**)
- **Type:** bare-pk
- **Role:** Inter-compartmental clearance between central and `peripheral3` (volume / time) in 4-compartment popPK models.
- **Source aliases:** none.
- **Example models:** `Schmitt_2018_vinflunine.R`, `Li_2017_brentuximab.R`, `Weatherley_2009_maraviroc_iv.R`.

### q_p1_p2 (**canonical bare one-way peripheral1-to-peripheral2 inter-compartmental clearance**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lq_p1_p2`: one-way inter-compartmental clearance from `peripheral1` into `peripheral2` (volume / time).
- **Source aliases:**
  - `Q3` -- used in `Tiraboschi_2023_avalglucosidaseAlfa.R`.
- **Example models:** `Tiraboschi_2023_avalglucosidaseAlfa.R` (founding example).
- **Notes:** See `lq_p1_p2` for why this is not `q2`.

### q_p2_c (**canonical bare one-way peripheral2-to-central back-redistribution clearance**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lq_p2_c`: one-way back-redistribution clearance from `peripheral2` to `central` (volume / time).
- **Source aliases:**
  - `Qpc` / `QPC` -- used in `Tiraboschi_2023_avalglucosidaseAlfa.R`.
- **Example models:** `Tiraboschi_2023_avalglucosidaseAlfa.R` (founding example).
- **Notes:** See `lq_p2_c`.

### q_bone_rbc, lq_bone_rbc (**canonical bare one-way bone-to-red-cell erythropoietic utilisation flow**)
- **Type:** bare-pk
- **Role:** One-way inter-compartmental clearance (volume / time) carrying a haematinic substrate out of the `bone` compartment and into the red-cell pool as it is incorporated into haemoglobin during erythropoiesis. It appears twice in a whole-body PBPK mass balance: as an extra loss term on `bone` (alongside the perfusion-limited efflux `q_bone / kp_bone`) and as the sole input to the red-cell pool. Because it is a flow rather than a rate constant, it multiplies a bone *concentration*, not a bone amount.
- **Source aliases:**
  - `QE` -- used in `Fan_2025_iron_mouse_ironadequate_pbpk.R` and siblings ("production rate of red blood cells in the bone compartment by utilizing iron").
- **Example models:** `Fan_2025_iron_mouse_irondeficient_pbpk.R`, `Fan_2025_iron_mouse_ironadequate_pbpk.R`, `Fan_2025_iron_mouse_ironloaded_pbpk.R`, `Fan_2025_iron_ferriccarboxymaltose_rat_pbpk.R`, `Fan_2025_iron_ferriccarboxymaltose_human_pbpk.R` (founding examples).
- **Notes:** Member of the `q_<from>_<to>` one-way inter-compartmental clearance family established by `q_p1_p2` and `q_p2_c`; the `<from>` / `<to>` tokens are canonical compartment names rather than the `p1` / `p2` / `c` positional abbreviations, because a whole-body PBPK model has no positional peripheral numbering. NOT a `kel`-style rate constant and NOT a `kin` / `kpro` production rate: those have units of 1/time or amount/time, whereas this is a clearance. The paired red-cell drain is `mtt_rbc` (mean red-cell lifespan), so the steady-state red-cell burden is `q_bone_rbc * C_bone * mtt_rbc`.

### kel (**canonical bare elimination rate constant (K-PD)**)
- **Type:** bare-pk
- **Role:** First-order elimination rate constant in K-PD / single-rate-constant elimination models with no explicit `vc`.
- **Source aliases:**
  - `ke` -- legacy.
  - `kde` -- paper-named (Mazzocco 2015 / Shoji 2017 / Xia 2024 KDE) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
  - `kp` -- paper-named (van Hasselt 2015 KP) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
  - `ps_elim`, `pc_elim` -- paper-named (Wilson 2015 p_S / p_C) bare drug-specific K-PD elim rates; replaced 2026-05-30 by `kel_sunitinib` / `kel_irinotecan`.
- **Example models:** `Mazzocco_2015_temozolomide.R`, `Shoji_2017_fosdagrocorat_oc.R`, `Shoji_2017_fosdagrocorat_p1np.R`, `vanHasselt_2015_eribulin.R`, `Wilson_2015_sunitinib_irinotecan_mouse.R` (bare drug-suffixed `kel_<drug>`), `Xia_2024_warfarin.R`.

### kel_exp_kdes (**canonical bare decay-rate constant of a time-varying elimination rate constant**)
- **Type:** bare-pk
- **Role:** Bare (natural-scale) form of `lkel_exp_kdes`: the rate constant at which a time-varying elimination rate constant relaxes toward its asymptote (1 / time).
- **Source aliases:**
  - `kappa_k` -- used in `Schreib_2024_busulfan.R`.
- **Example models:** `Schreib_2024_busulfan.R` (`kel_exp_kdes = exp(lkel_exp_kdes) = 0.0516 1/h`).
- **Notes:** Bare counterpart of `lkel_exp_kdes`; see that entry for the functional form. The sibling `kel_exp_famp` is registered under "Paper-named mechanistic parameters" instead of here, because it is signed and must never be log-transformed. Ratified on 2026-08-05 (task `oare_PMC11154452` sidecar question q1, answer A).

### k12 (**canonical bare central-to-first-peripheral rate constant**)
- **Type:** bare-pk
- **Role:** First-order distribution rate constant from central to `peripheral1` (1 / time). Used in rate-constant-parameterised 2- and 3-compartment models.
- **Source aliases:** none.
- **Example models:** rate-constant-parameterised popPK extractions.

### k21 (**canonical bare first-peripheral-to-central rate constant**)
- **Type:** bare-pk
- **Role:** First-order distribution rate constant from `peripheral1` to central (1 / time).
- **Source aliases:** none.
- **Example models:** rate-constant-parameterised popPK extractions.

### k13 (**canonical bare central-to-second-peripheral rate constant**)
- **Type:** bare-pk
- **Role:** First-order distribution rate constant from central to `peripheral2` (1 / time).
- **Source aliases:** none.
- **Example models:** rate-constant-parameterised 3-compartment popPK extractions.

### k31 (**canonical bare second-peripheral-to-central rate constant**)
- **Type:** bare-pk
- **Role:** First-order distribution rate constant from `peripheral2` to central (1 / time).
- **Source aliases:** none.
- **Example models:** rate-constant-parameterised 3-compartment popPK extractions.

### k_central_elf, k_elf_central (**canonical bare central-to-ELF and ELF-to-central rate constants**)
- **Type:** bare-pk
- **Role:** First-order transfer rate constants between `central` and the canonical `elf` compartment (1 / time). Deliberately a directional PAIR rather than a single inter-compartmental clearance, because plasma-to-ELF distribution is typically asymmetric -- `k_central_elf < k_elf_central` is what generates the sub-unity ELF penetration ratio that lung PK/PD analyses turn on, and a symmetric `q` / `velf` parameterisation cannot express it. The log-transformed `lk_central_elf` / `lk_elf_central` forms are used in `ini()`.
- **Source aliases:**
  - `K12` / `K21` -- Abouelhassan 2024 Results. The paper's prose calls ELF "the third compartment" while its printed subscripts number ELF second; the printed subscripts are the ones that reproduce the paper's own Table 4 PTA, so the ELF pair is the printed `K12` / `K21`. See the source-trace note in `Abouelhassan_2024_sulbactam_human.R`.
- **Example models:** `Abouelhassan_2024_sulbactam_human.R`.
- **Notes:** Member of the established `k_<from>_<to>` directional-transfer family already used for `k_central_milk` / `k_milk_central` (`Wattanakul_2024_primaquine.R`, `Wattanakul_2024_primaquine_motherinfant.R`) and `k_csf_plasma` / `k_csf_brain` / `k_brain_csf` (`Biliouris_2018_nusinersen.R`). Use this family, not `k13` / `k31`, whenever a transfer connects `central` to a named non-numbered compartment: `k13` / `k31` are reserved for `central` <-> `peripheral2` and reusing them for a physiologic matrix hides which compartment the constant actually feeds. Pairs with `lvelf`.

### fdepot (**canonical bare depot fraction**)
- **Type:** bare-pk
- **Role:** Fraction of dose absorbed via the depot (parallel-absorption models).
- **Source aliases:** none.
- **Example models:** parallel-absorption oral PK models.

### vmax (**canonical bare Michaelis-Menten Vmax**)
- **Type:** bare-pk
- **Role:** Maximum metabolic rate for saturable elimination (amount / time).
- **Source aliases:** none.
- **Example models:** TMDD and saturable-elimination popPK extractions.
- **Notes:** Replaces the deprecated `vm` name.

### bmax (**canonical bare maximum binding capacity**)
- **Type:** bare-pk
- **Role:** Maximum binding capacity of a saturable binding site or tissue, in the units of the saturating state (amount or concentration). See `lbmax` for the log-transformed primary form.
- **Source aliases:**
  - `BMAX` / `Bmax` -- upper-case source-paper forms.
- **Example models:** `Siebinga_2023_lu177psma617.R`, `GonzalezSales_2024_imetelstat.R`, `Nielsen_2011_cefuroxime.R`, `Svensson_2016_rifampicin.R`.
- **Notes:** Registered 2026-07-30 alongside `lbmax`.

### tmax_abs (**canonical bare saturable absorption Vmax**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `ltmax_abs`. Maximum rate of saturable Michaelis-Menten absorption from depot / absorption compartment into central (amount / time).
- **Source aliases:** none.
- **Example models:** `Jansson_2008_eflornithine_rat.R`.

### kt_abs (**canonical bare saturable absorption half-saturation amount**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lkt_abs`. Amount of drug in the depot / absorption compartment at which the saturable absorption rate is half-max (amount).
- **Source aliases:** none.
- **Example models:** `Jansson_2008_eflornithine_rat.R`.

### cl_ss (**canonical bare steady-state clearance arm**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lcl_ss`. Steady-state component of a time-varying clearance decomposition.
- **Source aliases:** none.
- **Example models:** time-varying-clearance popPK extractions.

### cl_time (**canonical bare time-varying clearance arm**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lcl_time`. Time-decay component of a time-varying clearance decomposition.
- **Source aliases:** none.
- **Example models:** time-varying-clearance popPK extractions.

### cl_renal (**canonical bare renal clearance arm**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lcl_renal`. Renal component of an additive renal + non-renal clearance decomposition.
- **Source aliases:** none.
- **Example models:** `Jonckheere_2019_cefepime.R`.

### cl_nonren (**canonical bare non-renal clearance arm**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lcl_nonren`. Non-renal component of an additive renal + non-renal clearance decomposition.
- **Source aliases:** none.
- **Example models:** `Jonckheere_2019_cefepime.R`.

### cl_hemodialysis (**canonical bare dialysis-active clearance arm**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lcl_hemodialysis`. Extracorporeal renal-replacement-therapy clearance arm; added to the body baseline `cl` only when the time-varying `HEMODIALYSIS` covariate is 1.
- **Source aliases:** `CL_HD`, `CLHD`, `CL_HF`, `CL_HDF`.
- **Example models:** `Veinstein_2013_gentamicin.R`.

### cl_met (**canonical bare metabolic-formation clearance arm**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lcl_met`. Metabolite-formation component (parent-to-metabolite mass flux) of an additive metabolic + non-metabolic clearance decomposition.
- **Source aliases:** none.
- **Example models:** `Lehr_2010_tesofensine.R`.

### cl_nonmet (**canonical bare non-metabolic-formation clearance arm**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lcl_nonmet`. Non-formation (non-metabolite-producing) elimination component of an additive metabolic + non-metabolic clearance decomposition.
- **Source aliases:** none.
- **Example models:** `Lehr_2010_tesofensine.R`.

### cl_2b6 (**canonical bare CYP2B6-mediated clearance arm**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lcl_2b6`. CYP2B6-mediated metabolic clearance arm; dynamically multiplied by `enzyme_2b6(t)` for the autoinduction contribution.
- **Source aliases:** none.
- **Example models:** `Heathman_2024_efavirenz.R`.

### cl_2a6 (**canonical bare CYP2A6-mediated clearance arm**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lcl_2a6`. CYP2A6-mediated metabolic clearance arm; dynamically multiplied by `enzyme_2a6(t)` for the autoinduction contribution.
- **Source aliases:** none.
- **Example models:** `Heathman_2024_efavirenz.R`.

### cl_ugt (**canonical bare UGT-mediated clearance arm**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lcl_ugt`. UGT (typically UGT2B7) -mediated metabolic clearance arm; constant in time (not autoinduced in the founding example).
- **Source aliases:** none.
- **Example models:** `Heathman_2024_efavirenz.R`.

### d1 (**canonical bare zero-order absorption duration**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `ld1`. Duration of the zero-order input rate into the depot (or central) compartment for sequential zero+first-order or pure zero-order absorption models.
- **Source aliases:**
  - `D1` -- NONMEM / paper notation.
  - `duration` -- alternate paper notation.
- **Example models:** `Aouri_2017_rilpivirine.R`, `Goggin_2004_emfilermin.R`, `Heathman_2024_efavirenz.R`.

### tlag (**canonical bare absorption lag time**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `ltlag`. Absorption lag time before drug enters the depot (time).
- **Source aliases:**
  - `alag` -- legacy.
- **Example models:** delayed-absorption oral PK models.

### tacro (**canonical bare acrophase**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `ltacro`. Peak-time parameter for circadian-IDR templates.
- **Source aliases:**
  - `tz` -- legacy.
- **Example models:** circadian indirect-response templates.

### clin (**canonical bare plasma-to-tissue influx clearance**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lclin`. Influx clearance from plasma central to a tissue ECF compartment.
- **Source aliases:** none.
- **Example models:** `Campagne_2019_cyclophosphamide_mouse.R`.

### clef (**canonical bare tissue-to-plasma efflux clearance**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lclef`. Efflux clearance from a tissue ECF compartment back to plasma central.
- **Source aliases:** none.
- **Example models:** `Campagne_2019_cyclophosphamide_mouse.R`.

### kamax (**canonical bare Weibull-absorption asymptotic maximum rate constant**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lkamax`. Maximum / asymptotic first-order absorption rate constant in a Piotrovskij-style time-dependent (Weibull) absorption model (1 / time).
- **Source aliases:**
  - `KAMAX` -- NONMEM convention.
- **Example models:** `Desai_2016_isavuconazole.R` (founding example).

### ra (**canonical bare Weibull-absorption rate-scaling parameter**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lra`. Rate-scaling parameter inside a Piotrovskij-style Weibull time-dependent ka (1 / time).
- **Source aliases:**
  - `RA` -- NONMEM convention.
- **Example models:** `Desai_2016_isavuconazole.R` (founding example).

### gam1 (**canonical bare Weibull-absorption shape parameter**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lgam1`. Unitless Weibull shape (sigmoidicity) parameter inside a Piotrovskij-style time-dependent ka.
- **Source aliases:**
  - `GAM1` / `GAMMA1` -- NONMEM convention.
- **Example models:** `Desai_2016_isavuconazole.R` (founding example).

### beta_cl (**canonical bare exponential-nonlinear-CL slope**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lbeta_cl`. Concentration-slope coefficient in an exponential-nonlinear-clearance function of the form `CL(C) = CL_0 * exp(beta_cl * C)`, units of 1 / concentration. Used inside `model()` after being exponentiated from `lbeta_cl`.
- **Source aliases:**
  - `β` -- Greek letter used in Larsen 2018 Table 3 and Eq. 1.
- **Example models:** `Larsen_2018_factorviii_rat.R`, `Larsen_2018_factorviii_monkey.R` (founding examples).

### cmpr (**canonical bare milk-to-plasma concentration ratio**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lcmpr`. Unitless ratio of breast-milk to maternal-plasma drug concentration, used inside `model()` after being exponentiated from `lcmpr` (`cmpr <- exp(lcmpr + etalcmpr) * <covariate terms>`) and then applied to derive the milk observable, `Cmilk <- cmpr * Cc`.
- **Source aliases:**
  - `MPRcon` -- Li 2023 ornidazole notation.
- **Example models:** `Li_2023_ornidazole.R` (founding example).
- **Notes:** See `lcmpr` above for the distinction from the `kp_<tissue>` partition-coefficient family and for why an AUC-based ratio is not encoded as a model parameter.

### kp_adipose, kp_bone, kp_brain, kp_cerebellum, kp_choroid_plexus, kp_cortex, kp_csf, kp_gut, kp_heart, kp_hippocampus, kp_intestine, kp_kidney, kp_liver, kp_lnode, kp_lung, kp_milk, kp_muscle, kp_other, kp_rest, kp_skin, kp_spleen, kp_striatum, kp_trachea (**canonical bare tissue-to-plasma partition coefficients**)
- **Type:** bare-pk
- **Role:** Bare counterparts of the `lkp_<tissue>` family. Unitless equilibrium ratio of concentration in one named anatomical tissue to concentration in plasma (or blood), used inside `model()` either as the exponentiated form of an `ini()` entry (`kp_liver <- exp(lkp_liver)`) or as a value computed in place from the source paper's physicochemical prediction equations.
- **Source aliases:**
  - `Kp` -- near-universal PBPK notation, subscripted by tissue.
- **Example models:** `Levitt_2005_propofol_pbpk.R`, `Mi_2023_cefquinome_pbpk.R`, `Gaohua_2012_pregnancy_pbpk_midazolam.R`, `Litjens_2023_linezolid_cns_pbpk.R`, `Yang_2023_diclazuril_chicken_pbpk.R`, `Kang_2023_artesunate_hamster_pbpk.R`.
- **Notes:** Registered 2026-08-05 together with the log-transformed family; see `lkp_adipose, ...` above for the family shape, for the `kp_milk` lactation member, for the boundary against `kpu<n>` / `sf<n>` / `cmpr`, and for the prefix-collision names that are not partition coefficients. Many PBPK models compute additional bare `kp_*` names that are qualifiers rather than tissues (`kp_free`, `kp_bound`, `kp_preg`) or per-tissue loop indices (`kp_i`); those are local `model()` intermediates and are intentionally not registered as canonicals.

---

## Paper-named mechanistic parameters

Parameters that don't fit the standard `ka` / `cl` / `vc` shape but recur across published models. Each entry is treated as a canonical bare name; the log-transformed form (`l<name>`) is acceptable wherever the parameter is strictly positive and the source paper reports an exponential typical-value form. Add to this list rather than introducing a new ad-hoc pattern.

### kel_exp_famp (**canonical fractional amplitude of a time-varying elimination rate constant**)
- **Type:** paper-named-param
- **Role:** Dimensionless fractional amplitude of an exponential change in the elimination rate constant over a course of therapy: `kel(t) = kel * (1 + kel_exp_famp * (1 - exp(-kel_exp_kdes * t)))`, so `kel(0) = kel`, `kel(inf) = kel * (1 + kel_exp_famp)`, and a negative value means `kel` (and hence `CL`) falls over time.
- **Source aliases:**
  - `dk` -- used in `Schreib_2024_busulfan.R` (paper `theta_dk1`).
- **Example models:** `Schreib_2024_busulfan.R` (`kel_exp_famp = -0.167`, a 16.7% average fall in busulfan `k` and `CL` at steady state; carries its own additive IIV `etakel_exp_famp` and an additive covariate effect `e_hlhxlp_kel_exp_famp = -0.145` that takes the amplitude to -0.312).
- **Notes:** `_famp` is the fractional-amplitude role token, introduced because the registered `cl_exp_` time-varying-clearance family parameterises by an absolute decaying component (`cl_exp_inf` + `cl_exp_component`) rather than by a fraction of the baseline. Registered as `paper-named-param` rather than `bare-pk` precisely because it is signed: it is estimated on the natural scale, is typically negative, and must NEVER be log-transformed (`bare-pk` makes `checkModelConventions()` demand an `l`-prefixed form). It is additive in its covariate and random-effect terms for the same reason. Prefer this fractional parameterisation over reparameterising into `cl_exp_inf` / `cl_exp_component` whenever the source estimates `V` and `kel` separately with *different* covariate sets, because the absolute-component form is a product of the two and entangles their covariates -- in the founding model the infusion-duration covariate has opposite signs on `ln(V)` and `ln(k)`, so the product would cancel it. The sibling decay rate is `kel_exp_kdes` / `lkel_exp_kdes` in the bare and log-transformed PK sections. The parallel `cl_exp_famp` was explicitly considered and NOT registered (sidecar option C, declined); add it only when a source paper needs it. Ratified on 2026-08-05 (task `oare_PMC11154452` sidecar question q1, answer A).

### kd (**canonical mechanistic dissociation / dissociation-like rate**)
- **Type:** paper-named-param
- **Role:** Dissociation rate constant or paper-mechanistic "kd" parameter (paper-specific meaning).
- **Source aliases:** none.
- **Example models:** TMDD models.
- **Notes:** Paper-mechanistic; verify the source paper's specific definition before reuse.

### kd0 (**canonical baseline mechanistic dissociation rate**)
- **Type:** paper-named-param
- **Role:** Baseline / initial value of a time-varying dissociation parameter.
- **Source aliases:** none.
- **Example models:** time-varying TMDD / receptor-occupancy models.

### k1 (**canonical association rate constant in reversible binding**)
- **Type:** paper-named-param
- **Role:** Second-order association (forward / on) rate constant for reversible drug-target or drug-drug complex formation, with units of (1 / (concentration * time)). The paired dissociation rate constant is `k2` and the equilibrium dissociation constant is `kd = k2 / k1`. Used by mechanistic 2-drug interaction models that need separate on / off rates rather than a single steady-state `kd`. When `kd` is fixed from an external (in-vitro) source and `k2` is estimated, `k1 = k2 / kd` is derived inside `model()` rather than estimated directly.
- **Source aliases:**
  - `k_on`, `kon` -- general kinetics notation.
- **Example models:** `Kleijn_2011_sugammadex_rocuronium.R` (sugammadex-rocuronium dynamic interaction: `kd` fixed at 0.0559 uM from in-vitro, `log(k2)` estimated, `k1 = k2 / kd = 0.61 1/(min*uM)` derived).
- **Notes:** Distinct from `kss` (the quasi-steady-state binding parameter used in QSS-TMDD approximations) and from `kd` (the equilibrium dissociation constant). Use the `k1` / `k2` / `kd` triple when the source paper reports the full dynamic-interaction kinetics; use `kss` alone for QSS-TMDD shortcuts.

### k2 (**canonical dissociation rate constant in reversible binding**)
- **Type:** paper-named-param
- **Role:** First-order dissociation (reverse / off) rate constant for reversible drug-target or drug-drug complex formation (1 / time). Paired with `k1` (association) and `kd = k2 / k1` (equilibrium dissociation). Inside `model()` the bare name `k2` is the rate constant; the log-transformed `lk2` form is used in `ini()` whenever the source paper reports `log_e(k2)` directly.
- **Source aliases:**
  - `k_off`, `koff` -- general kinetics notation.
- **Example models:** `Kleijn_2011_sugammadex_rocuronium.R` (sugammadex-rocuronium dynamic interaction: `lk2 = -3.38`, giving `k2 = 0.034 1/min`, RSE 16.5%; paired with the fixed `kd` and the derived `k1`).
- **Notes:** Pairs with `k1` (association). The paper-numerical convention `k_d = k_2 / k_1` should be encoded explicitly in `model()` so the dimensionally-correct relationship is visible at the call site.

### ks (**canonical drug-mediated effect-compartment elimination rate**)
- **Type:** paper-named-param
- **Role:** Second-order rate constant for one drug's modulation of another drug's elimination from an effect compartment, with units of (1 / (concentration * time)). Used in two-drug PK-PD interaction models where the modulating drug's plasma concentration multiplies a target drug's effect-compartment amount to drive an additional elimination route (mechanistic abstraction for site-of-action drug-drug interaction such as the sugammadex-mediated rocuronium reversal). Inside `model()` the bare name is `ks`; the log-transformed `lks` form is used in `ini()` when the source paper reports `log_e(ks)`.
- **Source aliases:** none.
- **Example models:** `Kleijn_2011_sugammadex_rocuronium.R` (`lks = -3.43`, giving `ks = 0.033 1/(min*uM)`, RSE 0.222%; modulates rocuronium effect-compartment elimination by sugammadex plasma concentration: `d/dt(effect_roc) = ... - ks * csug * effect_roc`).
- **Notes:** Mechanistically distinct from `kel` / `kdeg` (single-drug elimination rates) and from `kint` (target-mediated internalisation). The second drug's plasma concentration provides the multiplier; the parameter encodes the rate at which the modulating drug consumes the target drug at its site of action.

### ke0 (**canonical effect-compartment equilibration rate constant**)
- **Type:** paper-named-param
- **Role:** First-order rate constant for equilibration between the central plasma concentration and the effect-compartment concentration (1 / time), used by standard hysteresis effect-compartment PK-PD models: `d Ce / dt = ke0 * (Cc - Ce)`. The corresponding equilibration half-life is `log(2) / ke0`. Inside `model()` the bare name is `ke0`; the log-transformed `lke0` form is used in `ini()` when the source paper reports `log_e(ke0)` or uses an exponential typical-value form.
- **Source aliases:**
  - `keo`, `Keo` -- equivalent paper notation; both spellings (`keo` and `ke0`) appear in the literature.
  - `Kcb` -- Ojara 2024 notation for the plasma-to-breast-milk equilibration rate constant.
- **Example models:** `Kleijn_2011_sugammadex_rocuronium.R` (`lke0 = log(0.134) = -2.01`, giving `ke0 = 0.134 1/min` for the rocuronium effect-compartment equilibration; allometric scaling `(BW/70)^-0.25`), `Savic_2010_warfarin.R` (founding registry example -- 1-compartment PK with effect-compartment-driven proportional-odds PD), `Schindler_2016_sunitinib.R` (effect-compartment equilibration with daily AUC at ln(2)/50 1/h), `Ojara_2024_lamivudine.R` (`lke0 = log(0.245)` for plasma-to-breast-milk equilibration; the effect compartment is the physiologic `milk` matrix rather than a latent effect site).
- **Notes:** Distinct from `lke_kpd` (which is K-PD elimination rate and was deprecated in favour of the canonical `lkel`). `lke0` is specifically the effect-compartment equilibration parameter for hysteresis PK-PD modelling.

### ppc (**canonical effect-compartment pseudo-partition coefficient**)
- **Type:** paper-named-param
- **Role:** Steady-state ratio of the effect-compartment concentration to the central (plasma) concentration in a Sheiner-style effect-compartment model, `d Ce / dt = ke0 * (ppc * Cc - Ce)` (unitless fraction). Unlike a plain effect compartment, where the steady-state ratio is 1 by construction, `ppc` scales the driving concentration so the effect compartment equilibrates to a fraction (or multiple) of plasma -- the standard way to describe tissue penetration when only a delayed, partitioned site-of-action concentration is measured. Inside `model()` the bare name is `ppc`; the log-transformed `lppc` form is used in `ini()`, since the coefficient is strictly positive.
- **Source aliases:**
  - `PPC` -- Abdelgawad 2024 paper and control-stream notation (`PPC plasma-CSF`).
  - `Kp,uu`-style "penetration ratio" / "partition coefficient" prose in CSF- and tissue-penetration papers, when the quantity is the steady-state ratio of an effect compartment rather than a physiologic tissue-to-plasma partition.
  - `Rcb` -- Ojara 2024 notation for the lamivudine milk-to-plasma concentration ratio, the gain of a breast-milk partitioned effect compartment (`dA3/dt = Kcb * (Rcb * A2/VC - A3)`).
- **Example models:** `Abdelgawad_2024_linezolid.R` (founding example -- `lppc = log(0.365)`, the maximal CSF-to-plasma pseudo-partition coefficient for linezolid in tuberculous meningitis, modulated by CSF total protein through the covariate-effect parameter `e_csftpro_ppc`), `Ojara_2024_lamivudine.R` (`lppc = log(1.77)`, the milk-to-plasma concentration ratio of a breast-milk partitioned effect compartment; the first `ppc` above 1, exercising the "or multiple" case).
- **Notes:** Distinct from `ke0`, which sets how fast the effect compartment equilibrates, whereas `ppc` sets the level it equilibrates to; a model with a partitioned effect compartment needs both. Distinct from the PBPK `kpu<n>` family (clustered unbound tissue-to-plasma partition coefficients derived from tissue composition), which are physiologic constants of a whole-body model rather than an estimated effect-compartment parameter. Covariate effects on `ppc` follow the standard `e_<cov>_ppc` shape. **Lactation models: `ppc` is the right canonical when a milk-to-plasma ratio is the gain of a partitioned effect compartment**, and is distinct from the two other lactation ratio canonicals -- use `cmpr`/`lcmpr` when the milk observable is derived purely algebraically (`Cmilk <- cmpr * Cc`, no milk state), and `pcmilk` when the paper supplies a mechanistic partition coefficient entering a central-to-milk micro-rate constant (`k = (q_milk / vc) * pcmilk`). All three generate a milk:plasma concentration ratio at pseudo-equilibrium; the discriminator is which structure the paper actually parameterised.

### lec50 (**canonical log-transformed effect-compartment EC50**)
- **Type:** paper-named-param
- **Role:** Log-transformed concentration producing half-maximal effect in sigmoid Emax / Imax PD models (concentration units). Inside `model()` the bare name is `ec50`.
- **Source aliases:**
  - `lEC50`, `lec_50` -- equivalent paper notation.
- **Example models:** `deVriesSchultink_2018_trastuzumab_LVEF.R`, `Kleijn_2011_sugammadex_rocuronium.R`, `Zhang_2022_ormutivimab.R`.
- **Notes:** Pairs with `lhill` (Hill exponent). The bare name `ec50` is for use in `model()` derivations.

### ec50 (**canonical bare effect-compartment EC50**)
- **Type:** paper-named-param
- **Role:** Bare counterpart of `lec50`; the half-maximal effect concentration on the linear scale.
- **Source aliases:** none.
- **Example models:** widespread sigmoid-Emax PD extractions.

### lki50 (**canonical log-transformed half-maximal upstream-biomarker level driving a downstream effect rate**)
- **Type:** paper-named-param
- **Role:** Log-transformed half-maximal point of a sigmoid whose DRIVER is an upstream biomarker response rather than a drug concentration, entering a downstream effect-RATE expression: `rate = rmax * I^hill / (ki50^hill + I^hill)`, where `I` is the biomarker response on its own axis (percent inhibition of a phosphoprotein, percent receptor occupancy, fold-change of a transcript). Carries the units of that biomarker axis, NOT concentration units. Inside `model()` the bare name is `ki50`.
- **Source aliases:**
  - `KI50` -- Moein 2024 notation ("I where Ks is 50% of Kmax", Tables 3 and 4); constrained to (1, 100) in both NONMEM control streams because the driver is a percentage.
- **Example models:** `Moein_2024_apitolisib_mouse.R` (founding example; `ki50 = 67.2` %pAkt inhibition for the 786-O xenograft tumour-shrinkage sigmoid), `Moein_2024_apitolisib_human.R` (`ki50 = 58.0` %pAkt inhibition for the phase 1 RECIST sum-of-longest-diameters fit).
- **Notes:** Ratified 2026-08-05 alongside the Moein 2024 apitolisib extraction (sidecar request-001 q2, option A). Distinct from `lec50` / `lic50`, which are DRUG CONCENTRATIONS producing half-maximal effect -- a model can carry both at once, as the Moein files do (`ic50` on the drug-to-pAkt step, `ki50` on the pAkt-to-tumour step of the same cascade). Distinct also from `kc50` as used for a trimer CONCENTRATION at half-maximum kill rate: `ki50` is deliberately axis-agnostic and is the name to reach for when the sigmoid's x-axis is a biomarker readout. Pairs with `kmax` (the maximum of the driven rate) and `hill` (the sigmoidicity factor). A downstream user must read the axis units off `label()` / `compartmentData`, since they are model-specific by construction.

### ki50 (**canonical bare half-maximal upstream-biomarker level driving a downstream effect rate**)
- **Type:** paper-named-param
- **Role:** Bare counterpart of `lki50`; the half-maximal biomarker-response level on the linear scale, for use inside `model()`.
- **Source aliases:**
  - `KI50` -- Moein 2024 notation.
- **Example models:** `Moein_2024_apitolisib_mouse.R`, `Moein_2024_apitolisib_human.R`.

### lreg_qm (**canonical log-transformed once-monthly regimen EC50 multiplier**)
- **Type:** paper-named-param
- **Role:** Log-transformed multiplicative modifier applied to EC50 when a static Emax-on-AUC exposure-response model parameterised on an integrated-window AUC predictor (rather than an instantaneous concentration) needs to distinguish between dosing regimens whose AUCs encode different time-course profiles. Paired with the binary covariate `REGI_QM`: applied as `ec50_eff = ec50 * reg_qm^REGI_QM`, i.e., `ec50_eff = ec50` for Q2W (REGI_QM = 0) and `ec50_eff = ec50 * exp(lreg_qm)` for QM (REGI_QM = 1). Captures the "lack of kinetics" gap in static Emax-on-AUC models where two regimens that produce similar window-AUC values can still produce different effects because of differences in target-saturation time courses. Inside `model()` the bare name is `reg_qm`.
- **Source aliases:**
  - `REG` -- Kuchimanchi 2018 notation for the linear-scale multiplier (paper estimate 2.30, RSE 10.3%); the linear-scale form maps to `exp(lreg_qm)`.
- **Example models:** `Kuchimanchi_2018_evolocumab_ldlc.R` (founding example; evolocumab exposure-response Emax on LDL-C with `reg_qm = 2.30` distinguishing 420 mg SC QM from 140 mg SC Q2W regimens).
- **Notes:** Paired with the canonical covariate `REGI_QM`. Use only in static Emax-on-AUC exposure-response models that explicitly need a regimen-specific EC50; dynamic ODE-based PD models that resolve the concentration time course naturally do not need this parameter. Distinct from `lec50` (the half-maximal-effect concentration itself, with which `lreg_qm` is multiplicatively composed).

### reg_qm (**canonical bare once-monthly regimen EC50 multiplier**)
- **Type:** paper-named-param
- **Role:** Bare counterpart of `lreg_qm`; the linear-scale multiplicative modifier on EC50 for once-monthly dosing regimens in static Emax-on-AUC exposure-response models. Applied as `ec50_eff = ec50 * reg_qm^REGI_QM` so that `reg_qm` acts as the EC50 multiplier when `REGI_QM = 1` and the EC50 is unchanged when `REGI_QM = 0`.
- **Source aliases:**
  - `REG` -- Kuchimanchi 2018 notation.
- **Example models:** `Kuchimanchi_2018_evolocumab_ldlc.R`.

### le0 (**canonical log-transformed PD baseline parameter**)
- **Type:** paper-named-param
- **Role:** Log-transformed baseline (drug-free) value of a PD response, e.g., baseline T4/T1 twitch ratio in neuromuscular blockade PK-PD models (`E0 = exp(le0)`). Distinct from `lrbase`, which is the log-transformed steady-state baseline for indirect-response / turnover models with explicit `kin` / `kout`. Use `le0` when the source paper reports a non-turnover PD baseline (a typical T4/T1 ratio, a typical pre-treatment biomarker level) that enters the sigmoid Emax expression as an additive baseline plus the maximum-effect-bounded modulation. Inside `model()` the bare name is `e0`.
- **Source aliases:**
  - `lE0`, `lTOF0` -- equivalent paper notation.
- **Example models:** `Kleijn_2011_sugammadex_rocuronium.R` (`le0 = log(104)` for typical baseline T4/T1 x 100; Emax of the sigmoid is forced equal to E0 so the per-subject baseline shape of the readout is preserved).
- **Notes:** Distinct from `lrbase` (turnover-state steady-state baseline). When `Emax = E0` is forced in the sigmoid (the standard NMB parameterisation), the readout decreases monotonically from `E0` toward 0 as the effect-compartment concentration rises.

### e0 (**canonical bare PD baseline parameter**)
- **Type:** paper-named-param
- **Role:** Bare counterpart of `le0`; the baseline (drug-free) PD response value used inside `model()`.
- **Source aliases:**
  - `beta0`, `Intercept` -- the intercept of a binary logistic exposure-response regression `logit(p) = beta0 + beta1 * x`. The intercept IS the baseline (zero-exposure) response on the logit scale, so it takes the registered `logit` transform prefix: `logite0`. When one paper fits several parallel univariate regressions of the same endpoint against different exposure predictors, suffix by endpoint and predictor -- `logite0_<endpoint>_<predictor>` -- in the same way multi-output residual SDs are suffixed by output. The matching slope is a normal `e_<predictor>_<endpoint>` covariate-effect parameter. Used in `Assmus_2025_benznidazole_mouse.R` (`logite0_cure_auc` / `e_auc_cure` and `logite0_cure_tmic` / `e_tmic_cure`, Assmus 2025 Table 4).
- **Example models:** `Kleijn_2011_sugammadex_rocuronium.R`, `Assmus_2025_benznidazole_mouse.R`.

### kdes (**canonical desensitisation rate**)
- **Type:** paper-named-param
- **Role:** Receptor / target desensitisation rate constant (1 / time).
- **Source aliases:** none.
- **Example models:** receptor-desensitisation PD models.

### kdecay (**canonical decay rate**)
- **Type:** paper-named-param
- **Role:** First-order decay rate constant for a paper-mechanistic state.
- **Source aliases:** none.
- **Example models:** paper-mechanistic decay PD models.

### krel (**canonical release rate**)
- **Type:** paper-named-param
- **Role:** Release / liberation rate constant (e.g., drug release from a complex or depot reservoir).
- **Source aliases:** none.
- **Example models:** complex-release PD models.

### kss (**canonical steady-state dissociation parameter**)
- **Type:** paper-named-param
- **Role:** Steady-state binding parameter in quasi-steady-state TMDD approximations (concentration).
- **Source aliases:** none.
- **Example models:** QSS-TMDD popPK extractions.

### kint (**canonical internalisation rate**)
- **Type:** paper-named-param
- **Role:** First-order internalisation rate constant for receptor-mediated endocytosis of a drug-target complex (1 / time).
- **Source aliases:** none.
- **Example models:** TMDD models with explicit complex internalisation.

### kphag (**canonical threshold-gated phagocytic elimination rate**)
- **Type:** paper-named-param
- **Role:** First-order elimination rate constant for phagocytic clearance of antibody-coated target cells (platelets, RBCs, etc.) in TMDD models where phagocytosis is gated by a receptor-occupancy threshold: kphag acts on both the free-target and the drug-target-complex states while receptor occupancy exceeds the paper-fixed `thres` (percent), and is switched off below threshold. Distinct from the constitutive `kdeg` (baseline natural turnover), which is always on, and from `kint` (receptor-mediated complex internalisation, no threshold gate). Used in threshold-switch TMDD models where the paper reports a rate constant explicitly named "phagocytosis" alongside a receptor-occupancy threshold.
- **Source aliases:**
  - `KPH` -- Moc Willeford 2024 NONMEM abbreviation (supplement differential equations).
  - `kphagocytosis` -- Moc Willeford 2024 Table 1 / 2 typeset name.
- **Example models:** `MocWilleford_2024_rlyb212.R` (RLYB212 anti-HPA-1a IgG1 with threshold-gated phagocytic elimination of HPA-1a-positive platelets; founding example).
- **Notes:** Paired with the canonical `thres` receptor-occupancy threshold. When the paper reports a rate constant for phagocytic clearance without a threshold gate (i.e., always on), it is functionally equivalent to `kdeg` and should be recorded as such; `kphag` is reserved for the threshold-switched form.

### thres (**canonical receptor-occupancy threshold**)
- **Type:** paper-named-param
- **Role:** Threshold value of receptor occupancy (typically expressed as a percent) above which a threshold-gated elimination pathway (e.g., `kphag`) is switched on. Below the threshold, the gated pathway contributes zero to the target elimination rate; above the threshold, it contributes its full rate constant. Used in TMDD models where receptor engagement must exceed a minimum coating fraction before macrophage / RES-mediated phagocytic clearance is triggered.
- **Source aliases:**
  - `THRES` -- Moc Willeford 2024 NONMEM / typeset name.
- **Example models:** `MocWilleford_2024_rlyb212.R` (THRES = 10 percent fixed, describing the minimum receptor coating required to drive phagocytic clearance of HPA-1a-positive platelets; founding example).
- **Notes:** Usually held fixed at a value inferred from the paper's model-development results (Moc Willeford 2024: estimation ranged 3-22 percent when co-estimated with other parameters, so THRES was fixed at 10 percent). Distinct from a generic sigmoidal `hill` or `ec50` because the pathway is a hard switch rather than a smooth curve. When expressed as a percent, `thres` is compared to a percent receptor occupancy `ro = 100 * complex / (target + complex)`; when expressed as a fraction, `thres` is compared to the corresponding unitless ratio.

### qp (**canonical target inter-compartmental clearance**)
- **Type:** paper-named-param
- **Role:** Inter-compartmental clearance between the central and peripheral compartments of a target species (platelet, receptor, or other cell population tracked as an ODE state alongside the drug in a TMDD model). Units L / time. Distinct from the canonical drug `lq` / `q` because the two flows are semantically different: `q` moves drug between drug compartments, `qp` moves target between target compartments. Used when the source paper carries a peripheral distribution for the target species with its own inter-compartmental clearance parameter.
- **Source aliases:**
  - `QP` -- Moc Willeford 2024 Table 1 / 2 / supplement notation.
- **Example models:** `MocWilleford_2024_rlyb212.R` (QP = 2.45 L/h: platelet inter-compartmental clearance between the shared central compartment (V1) and the target peripheral compartment (V3); founding example).
- **Notes:** Encoded via `lqp` in `ini()` and `qp` in `model()`; scales allometrically with body weight at exponent 0.75 in the Moc Willeford 2024 simulation.

### vp_target (**canonical target peripheral volume**)
- **Type:** paper-named-param
- **Role:** Peripheral volume of distribution for a target species tracked as an ODE state alongside the drug in a TMDD model (units L). Distinct from the drug's canonical `lvp` / `vp` because the drug and the target may have physically different peripheral distribution volumes even when they share a central compartment. Paired with `qp` as the corresponding target inter-compartmental clearance.
- **Source aliases:**
  - `V3` -- Moc Willeford 2024 Table 1 / 2 / supplement notation for the target peripheral volume (the paper's V1 is the shared central, V2 is the drug peripheral, V3 is the target peripheral).
- **Example models:** `MocWilleford_2024_rlyb212.R` (V3 = 0.523 L: peripheral platelet distribution volume, paired with target central sharing V1; founding example).
- **Notes:** Encoded via `lvp_target` in `ini()` and `vp_target` in `model()`. Do not conflate with the drug's `lvp2` (second peripheral volume of a 3-compartment drug model); those are drug compartments, `vp_target` is a target compartment.

### kbm (**canonical biliary-metabolite excretion rate**)
- **Type:** paper-named-param
- **Role:** First-order rate constant for biliary excretion of a drug or metabolite from a plasma / central compartment into a downstream gut / bile compartment (1 / time). Used in enterohepatic-recirculation and interconversion submodels where the source paper carries biliary transport as a separate ODE flux (distinct from the drug's total clearance terms).
- **Source aliases:**
  - `k_bm` -- Hamren 2008 paper notation.
- **Example models:** `Hamren_2008_tesaglitazar.R` (kbm = 11.7 1/h: rate at which the acyl glucuronide metabolite is biliary-secreted from the metabolite plasma compartment into the gut compartment, where it then undergoes interconversion back to parent tesaglitazar; founding example).
- **Notes:** Ratified canonically on 2026-06-17 alongside the Hamren 2008 tesaglitazar extraction (per operator decision in sidecar request 001). Distinct from a generic intercompartmental clearance `lq` because biliary excretion in interconversion models is one-way (drug leaves the plasma compartment for the gut and does not return via the same route -- the return path is via a separate hydrolysis / reabsorption process parameterised by `kicv`).

### kicv (**canonical interconversion rate constant for metabolite-to-parent reverse-metabolism**)
- **Type:** paper-named-param
- **Role:** First-order rate constant for interconversion of a metabolite back to its parent drug (1 / time). Used for acyl-glucuronide-style reverse-metabolism mechanisms where the gut microbial / brush-border beta-glucuronidases hydrolyse a biliary-excreted conjugated metabolite back to the parent aglycone, which is then reabsorbed into systemic circulation. The single rate constant subsumes both the hydrolysis step and the reabsorption step (the source paper typically reports them as a single first-order process because the underlying gut microflora kinetics and intestinal-uptake kinetics cannot be separated from oral PK data alone).
- **Source aliases:**
  - `k_int` -- Hamren 2008 paper notation (the paper writes `kint` for interconversion; this naming clashes with the canonical `kint` for TMDD complex internalisation, so the canonical here is `kicv` to preserve semantic disambiguation).
- **Example models:** `Hamren_2008_tesaglitazar.R` (kicv = 0.79 1/h: rate at which acyl glucuronide accumulated in the gut compartment is hydrolysed and reabsorbed as parent tesaglitazar into the parent's central plasma compartment; founding example).
- **Notes:** Ratified canonically on 2026-06-17 alongside the Hamren 2008 tesaglitazar extraction (per operator decision in sidecar request 001). The canonical name `kicv` was chosen specifically to disambiguate from the TMDD `kint` canonical -- the two mechanisms (drug-target complex internalisation vs. metabolite-to-parent interconversion via gut hydrolysis) are semantically distinct, and reusing the same name with paper-specific meaning would create a latent reader-confusion failure mode. Paper-numbered alternatives (`k57`, `k72`-style as in deWinter 2009 MPA / MPAG) are discouraged for new extractions because the numbering is non-portable across papers; prefer the role-based `kbm` / `kicv` pair for any future interconversion popPK model.

### frac (**canonical fraction parameter**)
- **Type:** paper-named-param
- **Role:** Generic fraction (mixing weight, fraction-of-arm, etc.) bounded in (0, 1).
- **Source aliases:** none.
- **Example models:** mixture-distribution popPK extractions.

### alfm (**canonical mixing-mixture parameter**)
- **Type:** paper-named-param
- **Role:** Paper-mechanistic mixing parameter; paper-specific meaning.
- **Source aliases:** none.
- **Example models:** mixture models.

### ksyn (**canonical biomarker synthesis rate**)
- **Type:** paper-named-param
- **Role:** Zero-order or first-order synthesis rate constant for an endogenous biomarker turnover state.
- **Source aliases:** none.
- **Example models:** biomarker-turnover PD models.

### p (**canonical proliferation / growth-rate parameter**)
- **Type:** paper-named-param
- **Role:** Proliferation / growth-rate constant of a tumor-growth-inhibition or cell-population state. Paper-mechanistic meaning depends on the growth law. Inside `model()` the bare name is `p`; the log-transformed `lp` form is used in `ini()` because the rate is strictly positive.
- **Source aliases:**
  - `k` -- Yates 2023 notation for the growth-law rate constant of all three laws (exponential, Mayneord, von Bertalanffy).
- **Example models:** `Yates_2023_tgi_exponential.R` (founding example; `p = 0.005` 1/day, the exponential fractional growth rate), `Yates_2023_tgi_mayneord.R` (`p = 0.005` L^(1/3)/day), `Yates_2023_tgi_bertalanffy.R` (`p = 0.005` L^(1/3)/day).
- **Notes:** `p` is the TGI-family *generic*: use it whenever a paper's growth constant is not one of the mechanistically-specific alternatives. Distinct from `kgrow` (life-cycle-anchored multiplication rate constrained by a known cycle time and per-cycle amplification factor) and from `kge` (the growing-fraction exponent of the Stein two-exponential decomposition); both of those entries cross-reference `p` as the family generic. **Dimensions are set by the growth law, not by the name.** In an exponential law `dV/dt = p * V` the constant acts on volume and carries 1 / time; in the sub-exponential Mayneord and von Bertalanffy laws `dV/dt = p * V^(2/3) - ...` it acts on the tumor radius and carries length / time (L^(1/3)/day when V is in L). Record the dimensions in the parameter's `label()` rather than assuming 1 / time. Founding examples added 2026-08-05 with the Yates 2023 extraction; the entry itself predates them. Paired with `kdeath` in the von Bertalanffy law.

### vd (**canonical apparent volume of distribution**)
- **Type:** paper-named-param
- **Role:** Apparent volume of distribution in 1-compartment K-PD parameterisations or paper-mechanistic volume terms where the standard `vc` shape doesn't apply.
- **Source aliases:** none.
- **Example models:** K-PD / paper-mechanistic popPK extractions.

### kcat (**canonical catalytic rate constant**)
- **Type:** paper-named-param
- **Role:** Catalytic rate constant for an enzyme-mediated conversion (1 / time).
- **Source aliases:** none.
- **Example models:** enzyme-mediated metabolism PD models.

### kpro (**canonical production rate constant**)
- **Type:** paper-named-param
- **Role:** Zero-order production rate constant for a paper-mechanistic synthesis term.
- **Source aliases:** none.
- **Example models:** synthesis-rate PD models.

### krmr (**canonical removal-rate parameter**)
- **Type:** paper-named-param
- **Role:** First-order removal / clearance rate constant for a paper-mechanistic state.
- **Source aliases:** none.
- **Example models:** paper-mechanistic removal PD models.

### kgrow (**canonical growth-rate parameter**)
- **Type:** paper-named-param
- **Role:** First-order growth / net-multiplication rate constant of a paper-mechanistic proliferating state (1 / time). Founding use: fixed net asexual-parasite growth rate `kgrow = ln(10) / 48 = 0.0479 /h` for a Plasmodium falciparum life-cycle model in which parasites multiply 10-fold per 48-h intraerythrocytic cycle (Hien 2017 Table 3 row `K_grow (1/h) = 0.0479 fix`). Inside `model()` the bare name is `kgrow`; the log-transformed `lkgrow` form is used in `ini()` when the rate itself is log-parameterised (typical for positivity constraint).
- **Source aliases:**
  - `K_grow`, `Kgrow` -- Hien 2017 notation.
- **Example models:** `Hien_2017_cipargamin.R` (founding example; `kgrow` fixed at ln(10)/48 = 0.0479 /h for the 10-fold per 48-h cycle Plasmodium falciparum multiplication rate).
- **Notes:** Mechanistically distinct from `p` (generic proliferation / growth-rate constant, TGI-family) in that `kgrow` is specifically a life-cycle-anchored multiplication rate constrained by an independently-known cycle time and per-cycle amplification factor -- typically fixed rather than estimated. Also distinct from `kin` / `ksyn` (zero-order production rates into a turnover pool) because `kgrow` is a first-order growth rate proportional to the state itself. Ratified canonically on 2026-07-08 alongside the Hien 2017 cipargamin extraction.

### kge (**canonical Stein bi-exponential growth-rate constant**)
- **Type:** paper-named-param
- **Role:** First-order growth-rate constant of the treatment-resistant fraction in the Stein (2011) bi-exponential tumor-dynamics model, `y(t) = y0 * (exp(-kse * t) + exp(kge * t) - 1)` (1 / time). Drives the `growth` sub-state: `d/dt(growth) = kge * growth` with `growth(0) = y0`. Paired with `kse`. Inside `model()` the bare name is `kge`; the log-transformed `lkge` form is used in `ini()` because the rate is strictly positive. Per-treatment-arm variants take an arm suffix (`lkge_pembro`, `lkge_chemo`, ...); per-endpoint variants in a joint multi-biomarker model take an endpoint suffix (`kge_ctdna`).
- **Source aliases:**
  - `KG`, `TVKG`, `kg`, `k_g`, `kgT` -- Stein-family paper notation for the growth rate (`kgT` when the paper needs to distinguish a tumor-size growth rate from a second endpoint's).
- **Example models:** `Struemper_2025_tumorsize_OS_nsclc.R` (founding example; 12 per-arm `lkge_<arm>` values, 1/week), `Ribba_2022_sld.R` (`kge = 0.0016` 1/day for the OAK sum-of-longest-diameters fit).
- **Notes:** Register entry backfilled 2026-07-28; the name has been in use since Struemper 2025 (2026-06-28) and is described inside the `growth` / `shrink` entries of `compartment-names.md`, but had no entry of its own. Distinct from `kgrow` (life-cycle-anchored multiplication rate constrained by a known cycle time) and from `p` (generic TGI proliferation constant): `kge` is specifically the growing-fraction exponent of the two-exponential Stein decomposition, whose defining feature is that growth and shrinkage act on two independent sub-populations that are summed rather than on a single net-rate state. The related time-to-tumor-growth summary statistic is `TTG = (log(kse) - log(kge)) / (kge + kse)`.

### kse (**canonical Stein bi-exponential shrinkage-rate constant**)
- **Type:** paper-named-param
- **Role:** First-order decay-rate constant of the treatment-sensitive fraction in the Stein (2011) bi-exponential tumor-dynamics model (1 / time). Drives the `shrink` sub-state: `d/dt(shrink) = -kse * shrink` with `shrink(0) = y0`. Paired with `kge`. Inside `model()` the bare name is `kse`; the log-transformed `lkse` form is used in `ini()`. Arm- and endpoint-suffixed variants follow the same pattern as `kge`.
- **Source aliases:**
  - `KS`, `TVKS`, `ks`, `k_s`, `ksT` -- Stein-family paper notation for the shrinkage / decay rate.
- **Example models:** `Struemper_2025_tumorsize_OS_nsclc.R` (founding example; 12 per-arm `lkse_<arm>` values, 1/week), `Ribba_2022_sld.R` (`kse = 0.0014` 1/day for the OAK sum-of-longest-diameters fit).
- **Notes:** Register entry backfilled 2026-07-28 together with `kge`. **Name-collision warning:** the bare symbol `ks` used by much of the Stein literature is already taken in this register by a different canonical -- `ks` is the drug-mediated effect-compartment elimination rate of Kleijn 2011, with units 1 / (concentration * time). Always map a Stein-model paper's `ks` to `kse`, never to `ks`.

### kdeath (**canonical constitutive cell-death / loss rate**)
- **Type:** paper-named-param
- **Role:** Constitutive, drug-independent first-order cell-death / loss rate of a tumor state (1 / time). Founding use: the `- kd * V` loss term of the von Bertalanffy growth law `dV/dt = p * V^(2/3) - kdeath * V`, which is what makes an untreated tumor plateau at the finite carrying volume `(p / kdeath)^3`. Inside `model()` the bare name is `kdeath`; the log-transformed `lkdeath` form is used in `ini()` because the rate is strictly positive.
- **Source aliases:**
  - `kd` -- Yates 2023 notation for the von Bertalanffy natural cell-death rate.
- **Example models:** `Yates_2023_tgi_bertalanffy.R` (founding example; `kdeath = 0.004` 1/day, giving an untreated plateau volume of `(0.005 / 0.004)^3 = 1.95` L), `tgi_sat_VonBertalanffy.R`, `tgi_sat_genVonBertalanffy.R` (both templates, retrofitted from `lkd` 2026-08-06).
- **Notes:** Ratified 2026-08-05 with the Yates 2023 extraction. **Deliberately not named `kd`, and deliberately not merged into `kk`.** The bare `kd` is already this register's canonical for the mechanistic *dissociation* rate constant and is in live use with that meaning (`Kleijn_2011_sugammadex_rocuronium.R`, `Hayashi_2007_omalizumab.R`, `Lowe_2009_omalizumab.R`, `Anbari_2023`); overloading it would give one canonical two incompatible roles. The templates `tgi_sat_VonBertalanffy.R` and `tgi_sat_genVonBertalanffy.R` originally used `lkd` for this concept -- they predated the register and were the collision, not the precedent; both were retrofitted to `lkdeath` / `kdeath` on 2026-08-06, so no model in the library now uses `kd` for a cell-death rate. Note that `lkd` is still correct, and widely used, for the *dissociation* constant (~18 models). Distinct also from Cardilin 2018's `kk`, which is a *drug- or radiation-modulated* kill rate rather than the constitutive loss term encoded here, and from `kdeg` (first-order degradation of a turnover pool, a synthesis/degradation-balance concept rather than a growth-law loss term). Paired with `p`.

### kge_ctdna (**canonical Stein ctDNA growth-rate constant**)
- **Type:** paper-named-param
- **Role:** Endpoint-suffixed `kge` for the circulating-tumor-DNA arm of a Stein bi-exponential biomarker model (1 / time). Drives `d/dt(growth_ctdna) = kge_ctdna * growth_ctdna`. The `_ctdna` suffix is required whenever a single model fits both a tumor-size and a ctDNA Stein pair, so the two growth rates are unambiguous. Log-transformed form `lkge_ctdna` in `ini()`.
- **Source aliases:**
  - `kg` -- Ribba 2022 notation for the ctDNA growth rate in both Eq. 1 and Eq. 2.
- **Example models:** `Ribba_2022_ctdna.R` (founding example; `kge_ctdna = 0.0038` 1/day, RSE 30.1%), `Ribba_2022_ctdna_sld_joint.R` (same value, fixed).
- **Notes:** Registered 2026-07-28. See `kge` for the underlying Stein decomposition and `growth_ctdna` in `compartment-names.md` for the paired state.

### kse_ctdna (**canonical Stein ctDNA decay-rate constant**)
- **Type:** paper-named-param
- **Role:** Endpoint-suffixed `kse` for the circulating-tumor-DNA arm of a Stein bi-exponential biomarker model (1 / time). Drives `d/dt(shrink_ctdna) = -kse_ctdna * shrink_ctdna`. May be estimated directly (log-transformed `lkse_ctdna` in `ini()`) or derived inside `model()` from a tumor-size decay rate via the cross-endpoint link parameter `zeta`.
- **Source aliases:**
  - `ks` -- Ribba 2022 Eq. 1 notation for the freely-estimated ctDNA decay rate. Do NOT map this to the unrelated canonical `ks` (Kleijn 2011 effect-compartment elimination rate).
- **Example models:** `Ribba_2022_ctdna.R` (founding example; freely estimated, `kse_ctdna = 0.0081` 1/day, RSE 27.4%), `Ribba_2022_ctdna_sld_joint.R` (derived as `kse_ctdna <- zeta * kse`, not estimated).
- **Notes:** Registered 2026-07-28. See `kse` for the underlying Stein decomposition and the `ks` name-collision warning.

### zeta (**canonical cross-endpoint decay-rate link**)
- **Type:** paper-named-param
- **Role:** Dimensionless multiplier that ties one endpoint's Stein decay-rate constant to a second, jointly-modeled endpoint's decay-rate constant, so the two biomarkers share a single underlying treatment-response process rather than decaying independently. Founding use: `kse_ctdna = zeta * kse` in Ribba 2022 Eq. 2, coupling the ctDNA decay rate to the tumor-size (sum-of-longest-diameters) decay rate. `zeta > 1` means the ctDNA signal falls faster than tumor size. Inside `model()` the bare name is `zeta`; the log-transformed `lzeta` form is used in `ini()` because the multiplier is strictly positive, with IIV `etalzeta`.
- **Source aliases:**
  - Greek `zeta` -- Ribba 2022 Eq. 2 typeset symbol.
- **Example models:** `Ribba_2022_ctdna_sld_joint.R` (founding example; `zeta = 1.94`, RSE 37.3%, `omega_zeta = 0.86`, RSE 35.0%).
- **Notes:** Ratified canonically on 2026-07-28 alongside the Ribba 2022 joint ctDNA / tumor-size extraction. Distinct from `frac` (a bounded 0-1 fraction of a pool) and from `alfm` (a mixture-model mixing proportion): `zeta` is an unbounded positive ratio between two rate constants of the same dimension, and is the identifiable free parameter of a joint fit whose remaining structural parameters are held at their single-endpoint values. Use this canonical for any joint-biomarker model that couples two endpoints through a shared rate constant scaled by an estimated factor; if a future model needs more than one such link, suffix by the linked endpoint (`zeta_<endpoint>`).

### kact (**canonical activation-rate parameter**)
- **Type:** paper-named-param
- **Role:** First-order activation rate constant for a paper-mechanistic dormant / refractory state to become active (1 / time). Founding use: rate constant `kact` at which refractory (drug-tolerant) Plasmodium falciparum parasites become active and re-enter the drug-sensitive pool in a two-population parasite clearance model (Hien 2017 Table 3 row `K_act (1/h) = 0.0987`). Inside `model()` the bare name is `kact`; the log-transformed `lkact` form is used in `ini()` when the rate itself is log-parameterised.
- **Source aliases:**
  - `K_act`, `Kact` -- Hien 2017 notation.
- **Example models:** `Hien_2017_cipargamin.R` (founding example; `kact = 0.0987 /h` for refractory-to-active first-order transition, with IIV 41.5% CV per Hien 2017 Table 3; drives the awakening term `+ kact * parasite_refractory` in the sensitive-pool ODE and the loss term `- kact * parasite_refractory` in the refractory-pool ODE).
- **Notes:** Mechanistically distinct from `kel` / `kdeg` (single-drug or single-pool elimination), from `kmet` (parent-to-metabolite conversion), and from `kint` (target-mediated internalisation) in that `kact` transfers mass from an inactive pool to an active pool of the same species (no drug binding, no metabolism, no target sequestration). Ratified canonically on 2026-07-08 alongside the Hien 2017 cipargamin extraction.

### fsen (**canonical bare drug-sensitive fraction**)
- **Type:** paper-named-param
- **Role:** Fraction of a paper-mechanistic asexual / cycling pool that is fully drug-sensitive at model initialisation, bounded in `[0, 1]`. Complement `(1 - fsen)` is the drug-refractory subpool. Inside `model()` the bare name is `fsen`; the log-transformed `lfsen` form is used in `ini()`. Founding use: population fraction of asexual Plasmodium falciparum parasites that are drug-sensitive at enrolment in a two-population parasite clearance model (Hien 2017 Table 3 row `F_sen (%) = 99.1`).
- **Source aliases:**
  - `F_sen`, `Fsen` -- Hien 2017 notation.
- **Example models:** `Hien_2017_cipargamin.R` (founding example; `fsen = 0.991` typical; large IIV 81.8% CV per Hien 2017 Table 3; seeds the initial condition of the two ODE pools: `parasite_sensitive(0) = fsen * PARA` and `parasite_refractory(0) = (1 - fsen) * PARA`).
- **Notes:** Encoded on the log scale (`lfsen = log(0.991)`) for consistency with the paper's log-normal `omega^2 = log(1 + CV^2)` reporting formula (Hien 2017 Table 3 footnote a). Individual `fsen_i = exp(lfsen + etalfsen)` can occasionally exceed 1 for extreme etas; this is a documented limitation of the log-normal-on-fraction encoding and is preferable to reinterpreting the paper's reported %CV on a logit scale. Downstream simulations that need strictly-bounded individual fractions can clamp `fsen_i` to `min(fsen_i, 1)` at post-processing. Distinct from `fu` (fraction unbound in plasma; time-invariant physicochemical / binding property), `fr` (Bergstrand-Karlsson mixed-model fraction of MAT in transit delay), and `fm` (fraction metabolised through a specific pathway). Ratified canonically on 2026-07-08 alongside the Hien 2017 cipargamin extraction.

### lfsen (**canonical log-transformed drug-sensitive fraction**)
- **Type:** paper-named-param
- **Role:** Log-transformed counterpart of `fsen` for use in `ini()`. Individual fraction `fsen_i = exp(lfsen + etalfsen)`.
- **Source aliases:** none.
- **Example models:** `Hien_2017_cipargamin.R`.
- **Notes:** See `fsen` for the full role description and the caveat on the log-transform-on-a-bounded-fraction encoding choice.

### lkgrow (**canonical log-transformed growth-rate parameter**)
- **Type:** paper-named-param
- **Role:** Log-transformed counterpart of `kgrow` for use in `ini()`.
- **Source aliases:** none.
- **Example models:** `Hien_2017_cipargamin.R` (`lkgrow = fixed(log(0.0479))`).

### lkact (**canonical log-transformed activation-rate parameter**)
- **Type:** paper-named-param
- **Role:** Log-transformed counterpart of `kact` for use in `ini()`.
- **Source aliases:** none.
- **Example models:** `Hien_2017_cipargamin.R` (`lkact = log(0.0987)`).

### mat (**canonical mean absorption time**)
- **Type:** paper-named-param
- **Role:** Mean absorption time for transit-absorption popPK parameterisations (time).
- **Source aliases:** none.
- **Example models:** `Svensson_2016_bedaquiline.R` (DDMODEL00000219), `Kovalenko_2020_dupilumab.R`.
- **Notes:** Pairs with `mtt`, `fr`, `ktr` in the mean-absorption-time / fraction-of-MAT family.

### mtt (**canonical mean transit time**)
- **Type:** paper-named-param
- **Role:** Mean transit time for a transit-absorption chain (time).
- **Source aliases:** none.
- **Example models:** transit-absorption popPK extractions.
- **Notes:** First-order transit rate constant `ktr = n_transit / mtt`. For a mean *cell* lifespan use the `mtt_<celltype>` forms below rather than bare `mtt`.

### mtt_rbc, lmtt_rbc (**canonical mean red-blood-cell lifespan**)
- **Type:** paper-named-param
- **Role:** Mean lifespan of circulating red blood cells (time). Drains a single-compartment red-cell pool at first-order rate `1 / mtt_rbc`, so the pool's mean residence time equals the lifespan. In iron / haematinic PBPK models this term is the erythrophagocytic recycling limb: senescent red cells are cleared by splenic macrophages and their haemoglobin iron returns to the `spleen` compartment.
- **Source aliases:**
  - `TRBC` -- used in `Fan_2025_iron_mouse_ironadequate_pbpk.R` and siblings ("lifespan of RBC").
  - `MTT` -- used in `Naik_2013_peginesatide.R` for the same quantity ("Mean red-blood-cell lifespan").
- **Example models:** `Fan_2025_iron_mouse_irondeficient_pbpk.R`, `Fan_2025_iron_mouse_ironadequate_pbpk.R`, `Fan_2025_iron_mouse_ironloaded_pbpk.R`, `Fan_2025_iron_ferriccarboxymaltose_rat_pbpk.R`, `Fan_2025_iron_ferriccarboxymaltose_human_pbpk.R`.
- **Notes:** Registers the explicit `mtt_<celltype>` form for red cells, matching the platelet form `mtt_pl` already used in `Boak_2014_*` and making the celltype visible without reading the label. `Naik_2013_peginesatide.R` uses bare `lmtt` for exactly this quantity and is a candidate to migrate. Distinct from `ltr` / `ltp` (`PerezRuixo_2008_epoetinAlfa.R`, `Wang_2010_romiplostim.R`), which are the lifespans of *reticulocytes* and of *precursor* cells in a multi-stage maturation chain rather than of the mature circulating cell. Where a maturation chain is divided into `n` transit compartments the per-compartment rate is `n / mtt_rbc`; in the founding examples the pool is a single compartment, so the rate is `1 / mtt_rbc`.

### fr (**canonical fraction of MAT in transit delay**)
- **Type:** paper-named-param
- **Role:** Fraction of mean absorption time accounted for by the transit-delay chain.
- **Source aliases:** none.
- **Example models:** `Svensson_2016_bedaquiline.R`.

### ktr (**canonical transit-chain rate constant**)
- **Type:** paper-named-param
- **Role:** First-order rate constant for transit-absorption chains; equals `n_transit / mtt` for a chain of length `n_transit`.
- **Source aliases:** none.
- **Example models:** transit-absorption popPK extractions.

### ntr (**canonical transit-chain length**)
- **Type:** paper-named-param
- **Role:** Number of absorption transit compartments in a Savic (2007) transit chain (unitless, and typically non-integer because it is estimated as the shape parameter of the equivalent gamma-density input). Closes the `mtt` / `ktr` / `ntr` family: `ktr = (ntr + 1) / mtt`. Inside `model()` the bare name is `ntr`; the log-transformed `lntr` form is used in `ini()`, since the count is strictly positive.
- **Source aliases:**
  - `NN` -- NONMEM control-stream notation for the transit-chain shape parameter (Abdelgawad 2024 `$THETA 6 NN`, Svensson 2018 rifampicin `THETA(9)`).
  - `N`, `n_transit`, `NTR` -- equivalent paper notation.
- **Example models:** `Abdelgawad_2024_linezolid.R` (founding example -- `lntr = log(5.68)`, an estimated chain length feeding rxode2's analytical `transit()`).
- **Notes:** `ntr` is the canonical name for new extractions. Do NOT use `lnn` / `nn` for a transit-chain length: the `lnn` / `nn_fix` form is reserved for the Wang and co-authors sigmoidicity exponent in specific BDE / morphine-like models, so reusing it here would be a role collision. Legacy files predating this entry (`Svensson_2018_rifampicin.R`, `Smythe_2013_gatifloxacin.R`) still carry `lnn` / `nn` and are not retrofitted. Distinct from `ktr`, which is the per-compartment rate, and from `mtt`, which is the mean time through the whole chain.

### kmet (**canonical metabolite-formation rate constant**)
- **Type:** paper-named-param
- **Role:** First-order rate constant for in-vivo formation of an active metabolite from the parent central compartment, used when the source paper parameterises metabolite formation independently of the parent's total clearance.
- **Source aliases:** none.
- **Example models:** `Krause_2017_selexipag.R` (parent -> ACT-333679 active metabolite).
- **Notes:** Fraction-metabolised is implicit (`kmet * Vc_parent / CL_parent`) rather than estimated as `fm`.

### fm (**canonical fraction metabolised**)
- **Type:** paper-named-param
- **Role:** Fraction of parent clearance routed to an active metabolite, estimated as a structural parameter in parent-metabolite joint popPK models. Unitless, bounded in (0, 1].
- **Source aliases:** none.
- **Example models:** `Danielak_2017_clopidogrel.R` (clopidogrel -> H4 active thiol, doi:10.1007/s00228-017-2334-z).
- **Notes:** Distinct from `kmet` (formation rate constant): `fm` is unitless and bounded; `kmet` has rate units.

### logitfm (**canonical logit-transformed fraction metabolised**)
- **Type:** paper-named-param
- **Role:** Logit-transformed fraction of parent clearance routed to an active metabolite. Used when the source paper's estimation routine holds FM on the logit scale so that FM is bounded in (0, 1) regardless of covariate + eta combinations. Inside `model()` the bare form is `fm = 1 / (1 + exp(-logitfm_ind))` where `logitfm_ind` collects the fixed effect, covariate shifts, and IIV on the logit scale.
- **Source aliases:** none.
- **Example models:** `Mitra_2026_ziftomenib.R` (base `logitfm <- fixed(0.14)` corresponding to `logit^-1(0.14) = 0.535`; additive shift `e_dis_healthy_logitfm = -1.62` on the logit scale for the healthy-volunteer cohort; IIV `etalogitfm ~ 0.280` on the logit scale).
- **Notes:** Follows the `logit`-transform-prefix family (`logitffo`, `logitemax`) applied to the existing `fm` canonical root. Prefer `logitfm` (paper's logit encoding) over `lfm` (log encoding) when the source paper's NONMEM control stream stores FM on the logit scale, because a log-scale encoding of a bounded quantity can leak above 1 under moderate eta or covariate values. Ratified canonically on 2026-07-24 alongside the Mitra 2026 ziftomenib extraction.

### lq_milk (**canonical log-transformed central-to-breast-milk inter-compartmental clearance**)
- **Type:** log-transformed-pk
- **Role:** Apparent inter-compartmental clearance between a central compartment and its paired `milk` compartment in a lactation-transfer model (volume / time). Extends the existing `lq` family with a destination token, on the same principle as the `kin_<compartment>` / `kout_<compartment>` family: the canonical bare `lq` / `lq2` mean exchange with `peripheral1` / `peripheral2` specifically, so reusing them for a breast-milk compartment would silently mislead. The bare form `q_milk` is used inside `model()`.
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine.R`, `Wattanakul_2024_primaquine_motherinfant.R` (`Q/F = 0.400 L/h`, shared between primaquine and carboxyprimaquine per Wattanakul 2024 Table 2 footnote c; the source paper assumed and encoded `Q/F_CPQ = Q/F_PQ`).
- **Notes:** Ratified 2026-08-06 with the library's first lactation model (operator sidecar `oare_PMC11078975` request-001 / response-001, question q2, option A). Metabolite forms take the standard suffix (`lq_milk_<metab>`) when a paper estimates the parent and metabolite transfer clearances separately.

### lq_elf (**canonical log-transformed central-to-lung-ELF inter-compartmental clearance**)
- **Type:** log-transformed-pk
- **Role:** Inter-compartmental clearance between a central compartment and the lung epithelial-lining-fluid space (volume / time). Second member of the `lq_<destination>` family founded by `lq_milk`, for the same reason: the bare `lq` / `lq2` mean exchange with `peripheral1` / `peripheral2` specifically, so a lung-ELF exchange needs its own destination token. The bare form `q_elf` is used inside `model()`. Note that the exchange is generally NOT symmetric in rate -- `q_elf` is a clearance applied to whichever compartment the drug is leaving, so the ELF-to-central micro-constant is `q_elf / v_elf` while the central-to-ELF one is `q_elf / vc`.
- **Source aliases:** `Q1` -- Kurup 2024 Table 1 ("Intercompartmental clearance ELF - central") and supplement eqs. (1)-(4).
- **Example models:** `Kurup_2024_DZIF10c.R` (`Q1 = 3.68e-5 L/h` at WT = 70 kg, allometrically scaled with a fixed exponent of 0.85).
- **Notes:** Ratified 2026-08-14 (operator sidecar `oare_PMC11621205` request-001 / response-001, question q2, option A). Applies whether the ELF space is carried as the bare `elf` compartment alone or as `elf` plus regional `elf_<region>` pools that share one volume; in the founding example one clearance serves all three pools.

### lcl_elf (**canonical log-transformed lung-ELF loss clearance**)
- **Type:** log-transformed-pk
- **Role:** First-order clearance of drug out of the lung epithelial-lining-fluid space by a route that is not distribution into the systemic circulation (volume / time) -- mucociliary escalator transport, local catabolism, or swallowing. It is a clearance FROM the `elf` compartment, so it is deliberately NOT one of the `lcl_<arm>` names (`lcl_renal`, `lcl_nonren`, `lcl_hemodialysis`, ...), all of which are additive arms of the CENTRAL clearance. Reading it as an arm of `cl` would misplace the sink by a whole compartment. The bare form `cl_elf` is used inside `model()`, and the loss micro-constant is `cl_elf / v_elf`.
- **Source aliases:** `CL_L`, `CLL` -- Kurup 2024 Table 1 ("Clearance from ELF (IT/INH)") and supplement eqs. (1)-(2).
- **Example models:** `Kurup_2024_DZIF10c.R` (`CL_L = 0.000412 L/h` at WT = 70 kg; the paper attributes it to mucociliary clearance and it is what gives an inhaled dose its ~2.5-day effective lung half-life against the ~20-day systemic half-life).
- **Notes:** Ratified 2026-08-14 (operator sidecar `oare_PMC11621205` request-001 / response-001, question q2, option A). A model that resolves several distinct non-systemic lung sinks should suffix them further rather than overload this name.

### lv_elf (**canonical log-transformed lung-ELF volume**)
- **Type:** log-transformed-pk
- **Role:** Volume of the lung epithelial-lining-fluid space (volume). Needed because the bare volume canonicals `lvc` / `lvp` / `lvp2` / `lvp3` are reserved for the central and numbered peripheral compartments of a classical-PK model, and an ELF volume is a named physiological space rather than a peripheral compartment. The bare form `v_elf` inside `model()` predates this entry -- it is already used as a hardcoded physiological constant in `Parmar_2023_spectinamide_1599_mouse_pbpk.R` -- so this registers the `l`-prefixed `ini()` form for models that estimate it.
- **Source aliases:** `V_L`, `VL` -- Kurup 2024 Table 1 ("ELF compartment volume") and supplement eqs. (1)-(3), (6).
- **Example models:** `Kurup_2024_DZIF10c.R` (`V_L = 0.0364 L` at WT = 70 kg, estimated; allometrically scaled with a fixed exponent of 1). Bare-form precedent: `Parmar_2023_spectinamide_1599_mouse_pbpk.R` (`v_elf = 1.00e-5 L`, a fixed mouse physiological constant).
- **Notes:** Ratified 2026-08-14 (operator sidecar `oare_PMC11621205` request-001 / response-001, question q2, option A). When a model carries regional `elf_<region>` pools alongside the bare `elf`, a single `v_elf` normally serves all of them and the whole-lung ELF concentration is the summed amount over that one volume; a model that gives each region its own volume should register `lv_elf_<region>` names.

### logitfdepot (**canonical logit-transformed depot fraction**)
- **Type:** log-transformed-pk
- **Role:** Logit-scale encoding of the fraction of a dose that reaches the depot compartment -- the same quantity as `lfdepot`, held on the logit scale so that it is bounded in (0, 1) regardless of the covariate and eta combination. Inside `model()` the bare form is `f = exp(phi) / (1 + exp(phi))` where `phi` collects the fixed effect, covariate shifts, and IIV on the logit scale; IIV is therefore additive on the logit scale, not multiplicative.
- **Source aliases:** `F1` (estimated as `THETA` with `F1 = EXP(PHI)/(1+EXP(PHI))`) -- Kurup 2024 Table 1 and supplement eqs. (8)-(9), (13).
- **Example models:** `Kurup_2024_DZIF10c.R` (inhaled lung deposition, human `F1h = 29.6%`; the species-specific sibling is `logitfdepot_macaque` at `F1m = 1.00%`).
- **Notes:** Ratified 2026-08-14 (operator sidecar `oare_PMC11621205` request-001 / response-001, question q2, option A). Follows the `logit`-transform-prefix family (`logitffo`, `logitemax`, `logitfm`) applied to the existing `fdepot` root. Prefer `logitfdepot` over `lfdepot` when the source paper's control stream estimates the fraction on the logit scale, for the reason given in the `logitfm` entry: a log-scale encoding of a bounded quantity can leak above 1 under moderate eta or covariate values. Species-specific and stratum-specific values take the usual suffix (`logitfdepot_macaque`, per the `lfdepot_macaque` precedent in `Nagy_2017_obiltoxaximab.R`), with the reference species keeping the bare name.

### pcmilk (**canonical milk:plasma partition coefficient**)
- **Type:** paper-named-param
- **Role:** Milk:plasma partition coefficient -- the fraction of the drug in the central compartment that is freely distributed into breast milk. Unitless, normally in (0, 1]. Enters the central-to-milk micro-rate constant as `k = (q_milk / vc) * pcmilk`, so at pseudo-equilibrium the milk:plasma concentration ratio (and hence the milk:plasma AUC ratio) equals `pcmilk` exactly. That identity is the natural falsifier for a lactation extraction: a published milk:plasma AUC ratio that does not equal the published partition coefficient means the model that produced it had an additional milk sink. Metabolite forms take the standard suffix (`pcmilk_<metab>`).
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine.R`, `Wattanakul_2024_primaquine_motherinfant.R` (`PC_PQ = 0.376`, `PC_CPQ = 0.00889` per Wattanakul 2024 Table 2; Table 3 reports the simulated milk:plasma AUC ratio as 0.376 (0.375-0.377), equal to `pcmilk`).
- **Notes:** Ratified 2026-08-06 (operator sidecar `oare_PMC11078975` request-001 / response-001, question q2, option A). Distinct from a "milk-to-plasma ratio" reported as a raw observed concentration ratio: `pcmilk` is a model parameter that generates such a ratio, not the ratio itself.

### cfcap (**canonical capillary:venous conversion factor**)
- **Type:** paper-named-param
- **Role:** Proportional conversion factor between a venous-plasma model prediction and the corresponding capillary (finger-prick) plasma prediction, estimated when a study assays both matrices and fits them simultaneously. Unitless; the capillary observable is `Ccap = Cc * cfcap`. Metabolite forms take the standard suffix (`cfcap_<metab>`).
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine.R`, `Wattanakul_2024_primaquine_motherinfant.R` (`CF_PQ = 0.898`, `CF_CPQ = 1.06` per Wattanakul 2024 Table 2; the associated observation variables are `Ccap` and `Ccap_cpq`).
- **Notes:** Ratified 2026-08-06 as the mechanical consequence of the `cfcap` concept approved in operator sidecar `oare_PMC11078975` question q2 option A -- a capillary:venous conversion factor is only usable if the capillary observable it produces also has a name. Capillary sampling is standard in field malaria and paediatric pharmacokinetics, so this is expected to recur. Do NOT reuse `cfcap` for a plasma:whole-blood or plasma:serum conversion; those are different matrices and should get their own canonical.

### kmilkinf (**canonical breast-milk-to-infant transfer rate constant**)
- **Type:** paper-named-param
- **Role:** First-order rate constant (1/time) transferring drug from a mother's `milk` compartment into the breastfed infant's `infant_depot` compartment during a feeding window in a mother-to-infant dyad model. Normally fixed to a large value so that essentially all of the milk-compartment amount transfers within the window rather than estimated, because the feed is fast relative to every other process in the model.
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine_motherinfant.R` (`MTINF_PQ = MTINF_CPQ = 100 /h`, fixed high so that more than 95% of the amount in the breast-milk compartment is transferred during the 24-minute feeding window per Wattanakul 2024 Results, 'Predicting infant concentrations').
- **Notes:** Ratified 2026-08-06 (operator sidecar `oare_PMC11078975` request-001 / response-001, question q2, option A).

### feed_n, feed_window, feed_first, milk_intake (**canonical breastfeeding-pattern constants**)
- **Type:** paper-named-param
- **Role:** The four design constants that define a breastfeeding schedule in a lactation-transfer model. `feed_n` is the number of feeds per day (giving a feeding cycle of `24 / feed_n` hours); `feed_window` is the duration of each feed (time); `feed_first` is the time of the first feed after time zero (time); `milk_intake` is the average daily breast-milk intake per kg of infant body weight (volume / mass / time), from which the milk-compartment volume is derived as `milk_intake * WT_INFANT / feed_n`. All four are normally `fixed()` -- they describe the simulated feeding scenario, not an estimated property of the drug -- and are the natural handles for a feeding-pattern sensitivity analysis.
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine.R`, `Wattanakul_2024_primaquine_motherinfant.R` (`feed_n = 10` feeds/day, `feed_window = 0.4` h, `feed_first = 1` h, `milk_intake = 0.15` L/kg/day, giving a 2.4 h cycle with a 24-minute feeding window; Wattanakul 2024 Eq. 1 and Results, 'Predicting infant concentrations').
- **Notes:** Ratified 2026-08-06 (operator sidecar `oare_PMC11078975` request-001 / response-001, question q2, option A). Kept as separate scalars rather than folded into one schedule object so that each can be varied independently in a sensitivity analysis, which is exactly what the founding paper does (Supplementary Figs. 2-4 vary feeding frequency and volume per feed).

### lmtt_infant (**canonical log-transformed breastfed-infant mean transit time**)
- **Type:** log-transformed-pk
- **Role:** Mean transit absorption time of the breastfed infant in a mother-to-infant dyad model, when the infant's absorption is not scaled from the mother's but taken from a separate paediatric source. Uses the `_infant` suffix in the same sense as the `infant_<canonical>` compartment namespace and the `WT_INFANT` / `AGE_INFANT` covariates: it marks the quantity as belonging to the dyad PARTNER rather than to the modelled subject, whose own mean transit time is the plain `lmtt`.
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine_motherinfant.R` (`MTT_infant = 0.706 h`, fixed from a paediatric literature model to account for more rapid absorption in infants, with 2 transit compartments so `ktr_infant = 3 / MTT_infant`; Wattanakul 2024 Methods, 'Predicting infant concentrations').
- **Notes:** Ratified 2026-08-06 as the mechanical extension of the `infant_` dyad-partner namespace approved in operator sidecar `oare_PMC11078975` questions q1 and q3 (option A both). Any other canonical parameter that a dyad model must carry separately for the infant partner takes the same `_infant` suffix; parameters that are merely SCALED from the mother's (allometric clearances and volumes) do not get their own `ini()` name at all, because they are derived inside `model()`.

### kin (**canonical indirect-response production rate**)
- **Type:** paper-named-param
- **Role:** Zero-order production rate of an indirect-response / turnover pool (Dayneka 1993; Jusko traditions).
- **Source aliases:** none.
- **Example models:** indirect-response PD models.
- **Notes:** Codified 2026-05-28 per the naming audit. The compartment-suffixed forms `kin_<compartment>` / `kout_<compartment>` (with `lkin_<compartment>` / `lkout_<compartment>` as the log-transformed primaries) are the canonical family for **tissue-exchange rate constants** in semi-physiologic / dosimetry models that parameterise transport as first-order rate constants rather than as clearances: `kin_<tissue>` is central-to-tissue uptake and `kout_<tissue>` is tissue-to-central return, both 1/time. `<compartment>` must be a canonical compartment from `compartment-names.md`. Prefer these role-based names over the source paper's numeric micro-constants (`k12`, `k21`, `k13`, ...) whenever the compartments are anatomically named -- the canonical `k12` / `k21` / `k13` / `k31` entries mean central-to-`peripheral1` / `peripheral2` exchange specifically, so reusing them for organ compartments would silently mislead. Founding examples: `Lindauer_2017_pembrolizumab.R` and `Yamazaki_2008_crizotinib_mouse.R` (`kin_tumor` / `kout_tumor`); `Siebinga_2023_lu177psma617.R` (`kin_salivary_gland`, `kin_kidney`, `kin_liver`, `kin_tumor`, `kin_other` and their `kout_` partners, replacing Siebinga's `k12`/`k21`/`k13`/`k31`/`k14`/`k41`/`k15`/`k51`/`k16`/`k61`).

### kout (**canonical indirect-response elimination rate**)
- **Type:** paper-named-param
- **Role:** First-order elimination rate of an indirect-response turnover pool (1 / time). Also the tissue-to-central return leg of the `kin_<compartment>` / `kout_<compartment>` tissue-exchange family -- see the `kin` entry above.
- **Source aliases:** none.
- **Example models:** indirect-response PD models; `Lindauer_2017_pembrolizumab.R`, `Siebinga_2023_lu177psma617.R` (tissue-exchange form).

### lkinf_rbc, kinf_rbc (**canonical first-order influx rate constant into the red-cell analyte pool**)
- **Type:** paper-named-param
- **Role:** First-order rate constant (1 / time) for LINEAR influx of drug or active metabolite from plasma into the intracellular red-cell pool (`rbc_<analyte>`, see `compartment-names.md`). Enters as `d/dt(rbc_<analyte>) = kinf_rbc * Cc - keff_rbc * rbc_<analyte>`, so the plasma concentration `Cc` must be expressed in the same concentration units as the red-cell state for the term to be dimensionally consistent.
- **Source aliases:**
  - `Kin` -- used in `Gebhard_2023_methotrexate.R` (paper symbol `K_in^MTX`).
- **Example models:** `Gebhard_2023_methotrexate.R` (`K_in^MTX = 0.031 1/day`), `Gebhard_2023_mercaptopurine_anc.R`.
- **Notes:** Registered 2026-07-30 (sidecar `oare_PMC10359452` request-001 q2, option A). The `kinf` stem is deliberately NOT `kin`: the canonical `kin` is a ZERO-ORDER production rate into an indirect-response turnover pool, which is a different quantity with different units, so reusing that stem would actively mislead. Distinct also from `clin` / `clef` (`Campagne_2019_cyclophosphamide_mouse.R`), which are plasma-to-tissue influx / efflux CLEARANCES (volume / time) paired with an `ecf` compartment; `kinf_rbc` / `keff_rbc` are first-order RATE CONSTANTS. The `_rbc` suffix marks the destination pool; no analyte suffix is carried because each model file holds a single drug arm. A future model fitting two red-cell arms jointly would extend to `lkinf_rbc_<analyte>`.

### lkeff_rbc, keff_rbc (**canonical first-order efflux / elimination rate constant out of the red-cell analyte pool**)
- **Type:** paper-named-param
- **Role:** First-order rate constant (1 / time) for loss of drug or active metabolite from the intracellular red-cell pool (`rbc_<analyte>`), combining efflux back to plasma, intracellular catabolism, and red-cell turnover into a single lumped rate. Enters as the `- keff_rbc * rbc_<analyte>` term of the red-cell ODE.
- **Source aliases:**
  - `Keff` -- used in `Gebhard_2023_methotrexate.R` and `Gebhard_2023_mercaptopurine.R` (paper symbols `K_eff^MTX`, `K_eff^6MP`).
- **Example models:** `Gebhard_2023_methotrexate.R` (`K_eff^MTX = 0.018 1/day`), `Gebhard_2023_mercaptopurine.R` (`K_eff^6MP = 0.041 1/day`), `Gebhard_2023_mercaptopurine_anc.R` (`K_eff^6MP = 0.050 1/day`).
- **Notes:** Registered 2026-07-30 (sidecar `oare_PMC10359452` request-001 q2, option A). Because the rate is lumped, its reciprocal half-life is a red-cell residence property rather than a pure membrane-transport property -- Gebhard 2023's Discussion validates `K_eff^MTX = 0.018 1/day` against literature red-cell methotrexate half-lives of 30-40 days (0.017-0.023 1/day).

### lvmax_rbc, vmax_rbc (**canonical saturable-influx maximum rate into the red-cell analyte pool**)
- **Type:** paper-named-param
- **Role:** Maximum rate (concentration / time, e.g. umol/L/day) of SATURABLE Michaelis-Menten influx from plasma into the intracellular red-cell pool: `d/dt(rbc_<analyte>) = vmax_rbc * Cc / (km_rbc + Cc) - keff_rbc * rbc_<analyte>`.
- **Source aliases:**
  - `Vmm` -- used in `Gebhard_2023_mercaptopurine.R` (paper symbol `V_mm^6MP`).
- **Example models:** `Gebhard_2023_mercaptopurine.R` (`V_mm^6MP = 0.096 umol/L/day`), `Gebhard_2023_mercaptopurine_anc.R` (`V_mm^6MP = 0.21 umol/L/day`).
- **Notes:** Registered 2026-07-30 (sidecar `oare_PMC10359452` request-001 q2, option A). Extends the blessed `vmax_<suffix>` / `km_<suffix>` disambiguation pattern to a saturable INFLUX. Distinct from the canonical bare `vmax` / log `lvmax`, which are registered for saturable ELIMINATION and carry amount/time units; `vmax_rbc` is an influx into a concentration state and therefore carries concentration/time units.

### lkm_rbc, km_rbc (**canonical saturable-influx half-saturation concentration for the red-cell analyte pool**)
- **Type:** paper-named-param
- **Role:** Michaelis constant (concentration, e.g. umol/L) of the saturable influx into the intracellular red-cell pool; the plasma concentration at which influx reaches half of `vmax_rbc`.
- **Source aliases:**
  - `Kmm` -- used in `Gebhard_2023_mercaptopurine.R` (paper symbol `K_mm^6MP`).
- **Example models:** `Gebhard_2023_mercaptopurine.R` (`K_mm^6MP = 0.016 umol/L`), `Gebhard_2023_mercaptopurine_anc.R` (`K_mm^6MP = 0.14 umol/L`).
- **Notes:** Registered 2026-07-30 (sidecar `oare_PMC10359452` request-001 q2, option A). Paired with `lvmax_rbc`; both are meaningless alone.

### linieff, inieff (**canonical mid-therapy initialisation fraction for a turnover / maturation chain**)
- **Type:** paper-named-param
- **Role:** Dimensionless fraction of the turnover baseline (`rbase`) at which a turnover or maturation chain is initialised when the observation record starts DURING ongoing therapy rather than at drug-free baseline. The terminal state is initialised at `inieff * rbase` and the upstream chain states at the steady-state values implied by that terminal value, so the chain begins at a treatment steady state carrying an unobserved historical drug effect rather than at the untreated baseline.
- **Source aliases:**
  - `inieff` -- used in `Gebhard_2023_mercaptopurine_anc.R` (Gebhard 2023 main text: "we assume to have reached a treatment steady state with an additional parameter `inieff` describing the drug effect at the initial time point and `X_ma(0) = inieff * base`").
- **Example models:** `Gebhard_2023_mercaptopurine_anc.R` (`inieff = 0.87`, RSE 22%, CV 54%; Friberg chain initialised mid-maintenance-therapy for childhood ALL).
- **Notes:** Registered 2026-07-30 (sidecar `oare_PMC10359452` request-001 q3, option A -- paper-faithful naming, following the register's existing paper-named precedents `le0` and `kmet`). Estimated with IIV in Gebhard 2023, which is what distinguishes it from a fixed initialisation convention: it absorbs the between-patient spread in how suppressed each patient's ANC already was when their record began. A value of 1 recovers the ordinary drug-free-baseline initialisation. Distinct from `rbase`, which is the untreated baseline itself.

### kdeg (**canonical degradation rate**)
- **Type:** paper-named-param
- **Role:** First-order degradation rate constant (paper-synonym for elimination in turnover models).
- **Source aliases:** none.
- **Example models:** turnover PD models.

### kpin (**canonical precursor-pool production rate**)
- **Type:** paper-named-param
- **Role:** Production rate of a precursor pool in `indirect_prec_*` template family.
- **Source aliases:** none.
- **Example models:** indirect-response with precursor templates.

### kpout (**canonical precursor-pool loss rate**)
- **Type:** paper-named-param
- **Role:** Loss rate of a precursor pool in `indirect_prec_*` template family.
- **Source aliases:** none.
- **Example models:** indirect-response with precursor templates.

### hill (**canonical Hill / sigmoid-shape coefficient**)
- **Type:** paper-named-param
- **Role:** Hill / sigmoid-shape coefficient in sigmoidal Emax / Imax functions: `Cc^hill / (ec50^hill + Cc^hill)`.
- **Source aliases:** none.
- **Example models:** sigmoidal Emax / Imax PD templates.
- **Notes:** Codified 2026-05-28 per the naming audit. Distinct from `gamma` for Friberg myelosuppression feedback / TGI power-law growth exponents, which retain `gamma` as a mechanistic-role designator.

### iplac (**canonical constant fractional placebo inhibition of an IDR production rate**)
- **Type:** paper-named-param
- **Role:** Constant fractional inhibition of an indirect-response production rate attributable to the placebo / background effect, entering as `kin * (1 - iplac - <drug effect>)`. Unitless and bounded in [0, 1]; carries no driver (no concentration, no biomarker) and is therefore NOT a maximum. Used in placebo-controlled longitudinal disease-score IDR models where the total inhibition of the input function decomposes additively into a placebo term and a drug term.
- **Source aliases:**
  - `PLAC`, `TVPLAC` -- Zhang 2023 NONMEM control streams (`PLAC = TVPLAC + ETA(9)`; `DADT(4) = KIN*(1-PLAC-EFF) - KOUT*A(4)`).
  - `Inhibitory placebo effect` -- Zhang 2023 Table 2 row label.
- **Example models:** `Zhang_2023_brazikumab_il22.R` (`iplac = 0.209`, i.e. placebo inhibits the CDAI input rate by 20.9%), `Zhang_2023_brazikumab_crp.R` (`iplac = 0.178`).
- **Notes:** Ratified 2026-07-29 alongside the Zhang 2023 brazikumab extraction (sidecar request-002 / response-002, option A). Kept on the bare (untransformed) linear scale rather than log-transformed because the published between-subject variability on it is ADDITIVE (Zhang 2023 Methods: "the placebo effect ... is modeled to have an additive variability (i.e., PLAC_i = PLAC_TV + eta_i), with the assumption that the placebo effect in any individual may be higher or lower than the typical value ... by an equal probability"), so an `etaiplac` enters as `iplac + etaiplac`. Distinct from `imax`, which is the MAXIMUM of a driver-dependent (sigmoid) inhibition; `iplac` is a driver-free constant. Distinct also from `lplac` as used in `Schoemaker_2018_levetiracetam.R`, which is a log-scale placebo effect on a seizure RATE rather than a fractional inhibition of an IDR input function. Later placebo-controlled IDR extractions in the IBD / CDAI family (certolizumab, ustekinumab, mirikizumab, ...) should reuse this canonical.

### rbase (**canonical baseline-value parameter**)
- **Type:** paper-named-param
- **Role:** Baseline-value parameter for IDR / turnover state initial conditions and TGI initial tumour sizes.
- **Source aliases:**
  - `r0` -- legacy.
  - `bl` -- legacy.
  - `base` -- legacy.
  - `s0` -- legacy.
  - `ts0` -- legacy.
- **Example models:** IDR / TGI template family.
- **Notes:** Codified 2026-05-28 per the naming audit (replaces the five legacy names listed above).

### fu (**canonical fraction unbound in plasma**)
- **Type:** paper-named-param
- **Role:** Fraction unbound in plasma; fixed unitless multiplier (in [0, 1]) used in cerebral-microdialysis-style CNS-distribution models where only free drug crosses the BBB.
- **Source aliases:** none.
- **Example models:** `Campagne_2019_cyclophosphamide_mouse.R`.
- **Notes:** Usually held fixed at the in-vitro equilibrium-dialysis-derived value. The BBB transfer term is `CLin * fu * Cp`. When the unbound fraction is itself concentration-dependent, `fu` remains the name of the DERIVED quantity computed inside `model()` and the saturable-binding submodel is parameterised with the `fumin` / `fumax` / `cup50` family below.

### fumin, fumax, cup50 (**canonical saturable plasma-protein-binding parameters**)
- **Type:** paper-named-param
- **Role:** The three parameters of a concentration-dependent (saturable) plasma protein binding submodel, in which the unbound fraction rises with the unbound plasma concentration along an Emax / Michaelis-Menten curve: `fu = fumin + fumax * Cu / (cup50 + Cu)`. `fumin` (unitless) is the limiting unbound fraction as the concentration approaches zero, i.e. the value when every high-affinity binding site is available; `fumax` (unitless) is the maximum INCREMENT added at saturating concentration, so the asymptotic unbound fraction is `fumin + fumax`, not `fumax`; `cup50` (concentration units) is the UNBOUND plasma concentration at which half of that increment has been reached. All three are normally held `fixed()`, because they come from an in-vitro or dedicated-study binding characterisation rather than from the popPK fit. The driver is the unbound concentration, which is what the `cup50` name records: in a free-drug-referenced model (`Cu = central / vc`) the measured TOTAL concentration is recovered as `Cc = Cu / fu`, so a `cu50` named against total drug would misstate the curve by a factor of `fu`.
- **Source aliases:**
  - `fu,min` / `fu, min`, `fu,max` / `fu, max`, `fu50` / `Cu,plasma 50` -- Bian 2024 Equations 16-19 and Supplementary Tables S1-S2; FDA XENLETA Multi-Discipline Review (NDA 211672/211673) section 16.3.2.4, which states the submodel and, uniquely, that the driver is the unbound rather than the total concentration.
- **Example models:** `Bian_2024_lefamulin_original_ppb.R` (founding example; `fumin = 0.0997`, `fumax = 0.259`, `cup50 = 1.35` mg/L, all fixed, from the diluted-plasma binding assay) and `Bian_2024_lefamulin_higher_ppb.R` (`fumin = 0.0260`, `fumax = 0.194`, `cup50 = 0.814` mg/L, from the non-diluted assay).
- **Notes:** Keep the bare lower-case names -- do NOT log-transform, because all three are held fixed and `fumin` / `fumax` are bounded fractions for which a log encoding could leak above 1. Distinct from `fu`, which is the derived unbound fraction this family produces (and which remains the right canonical for a time-invariant fixed unbound fraction). Distinct from `lbfu` / `bfu`, the slope of a linear TIME-varying unbound fraction: that family models drift over a study window, this one models saturation of binding sites with concentration. Distinct from the sigmoidal-PD `emax` / `ec50` canonicals: those describe a drug effect, whereas these describe a physicochemical binding equilibrium, and reusing `ec50` here would collide with a PD role. A future model that fits a Hill exponent on the binding curve should register a companion `fuhill` rather than reuse `hill`. Ratified canonically on 2026-08-09 alongside the Bian 2024 lefamulin extraction (sidecar `oare_PMC11544005` request-001 q1, option A).

### penratio_elf, pwr_penratio_elf (**canonical algebraic ELF penetration-ratio parameters**)
- **Type:** paper-named-param
- **Role:** The two parameters of an ALGEBRAIC (non-ODE) tissue-penetration function that returns the epithelial lining fluid concentration directly from the plasma concentration: `Celf = penratio_elf * (Cfree_plasma)^pwr_penratio_elf`. `penratio_elf` (unitless) is the lung penetration ratio evaluated at unit unbound plasma concentration -- the ELF-to-free-plasma concentration ratio at 1 mg/L. `pwr_penratio_elf` (unitless) is the exponent that lets the penetration ratio change with concentration: 1 recovers a plain proportional-partition model, and values below 1 mean the ratio falls as concentration rises (a saturable penetration process). Bare names are used inside `model()`; the log-transformed `lpenratio_elf` is used in `ini()` because the ratio is strictly positive, while the exponent stays on the linear scale.
- **Source aliases:**
  - `LPR (1 mg/L)` and `power` -- Bian 2024 Equation 21 and Supplementary Table S2; FDA XENLETA Multi-Discipline Review (NDA 211672/211673) section 16.3.2.4.1 "ELF PK Model", which introduces the link function as `C_ELF = LPR(1 mg/L) * [C_P(t) * 0.0379]^power`.
- **Example models:** `Bian_2024_lefamulin_higher_ppb.R` (founding example; `lpenratio_elf = log(3.45)`, `pwr_penratio_elf = 0.51`).
- **Notes:** Deliberately NOT named `lpr`: in this register a leading `l` means log-transformed (`lcl`, `lvc`, `lkdeath`), so `lpr` would read as "log of pr" for a parameter that is not log-transformed -- the same class of trap as the `lkd` / `kdeath` collision. Deliberately NOT folded into `ppc`: `ppc` is precisely defined as a `ke0`-driven effect-compartment pseudo-partition coefficient reached over time, whereas this pair describes an instantaneous algebraic penetration ratio with no dynamics, and widening `ppc` would give one canonical two roles. The `pwr_` prefix on the exponent keeps it from being read as a Hill / sigmoidicity term. Use the ODE form (`kin_elf` / `kout_elf`, from the `kin_<compartment>` / `kout_<compartment>` tissue-exchange family) whenever the paper gives ELF its own compartment with rate constants; use this pair only when the paper's ELF layer is a closed-form function of the plasma concentration. Generalises to other matrices as `penratio_<compartment>` / `pwr_penratio_<compartment>`, where `<compartment>` must be a canonical compartment from `compartment-names.md`. Ratified canonically on 2026-08-09 alongside the Bian 2024 lefamulin extraction (sidecar `oare_PMC11544005` request-001 q2, operator rename ruling).

### mic (**canonical minimum inhibitory concentration**)
- **Type:** paper-named-param
- **Role:** Minimum inhibitory concentration of the drug against the target pathogen, in the model's concentration units (note that mg/L and ug/mL are numerically identical, which is why antimicrobial papers move between them freely). Used as the denominator of the PK/PD indices that drive antimicrobial PK/PD-integration, probability-of-target-attainment and PK/PD-cutoff models: `AUC/MIC`, `fAUC/MIC`, `Cmax/MIC`, and as the threshold for `T>MIC` / `fT>MIC`. Normally held `fixed()`, because it is a measured property of the isolate (CLSI / EUCAST broth microdilution) rather than an estimated model parameter; papers that sweep a range of MICs expose it so the user can override it per simulation.
  Because the same isolate has different MICs by different reference methods, a model whose exposure-response was fitted against one susceptibility methodology must state which one in the `label()`; when a paper fits parallel exposure-response relationships against two methodologies, those are separate model files (see `Beredaki_2023_micafungin_clsi.R` / `..._eucast.R`).
- **Source aliases:**
  - `MIC` -- the near-universal paper symbol.
  - `MIC50`, `MIC90` -- population-distribution percentiles; use `mic` only for a single organism's own MIC, and record which percentile the value is in the `label()` and the source-trace comment.
  - `MEC` -- minimum effective concentration, the echinocandin/mould analogue of the MIC, reported for filamentous fungi.
  - `IC90`, `IC90,plasma`, `IC50` -- the in vitro inhibitory-concentration analogue reported for intracellular protozoa and other organisms for which broth microdilution is not applicable, and which are assayed by image-based percent-inhibition readouts instead. Functionally identical in the model: a measured in vitro potency held `fixed()` and used as the threshold of a `T>threshold` index. Record the percent-inhibition level, the strain, the host-cell line, and any protein-binding correction in the `label()` and the source-trace comment, because an in vitro IC90 read in assay medium is not the same number as the plasma-equivalent threshold. Used in `Assmus_2025_benznidazole_mouse.R` (6.427 ug/mL: in vitro IC90 against *Trypanosoma cruzi* amastigotes, Tulahuen strain, in 3T3 host cells, corrected for protein binding to assay medium and mouse plasma).
- **Example models:** `Chen_2023_tilmicosin.R` (fixed at 0.25 ug/mL, the CLSI 2013 broth-microdilution MIC of *Pasteurella multocida* strain C44-15, serovar D:7; forms the `AUC24h/MIC` index that drives a sigmoidal Emax kill model), `Lallemand_2023_benzylpenicillin_horse.R` (fixed at 0.25 mg/L, the paper's concluded PK/PD cutoff, and overridden via `params = c(mic = ...)` to sweep 0.0625-2 mg/L; used both as the `fAUC/MIC` denominator and as the threshold integrated by the `t_above_mic` state), `Beredaki_2023_micafungin_clsi.R` (0.008 mg/L, CLSI M27 median MIC of *Candida albicans* CA 580) and `Beredaki_2023_micafungin_eucast.R` (0.016 mg/L, EUCAST E.Def 7.3 median MIC of the same isolate), both driving an `fAUC0-24/MIC` sigmoidal Emax.
- **Notes:** Keep the bare lower-case name -- do NOT log-transform (an `lmic` would imply the MIC is being estimated) and do NOT rename to a generic `thres` / `ec50`, which would lose the microbiological meaning and collide with the sigmoidal-PD canonicals. Distinct from `ec50` / `lec50`: `ec50` is the *fitted* index value producing half-maximal effect (and in an `AUC/MIC`-indexed model is unitless, or has units of time), whereas `mic` is a *measured* concentration used to form the index in the first place. Pair with `fu` when the index is defined on free rather than total drug: apply the unbound fraction to the exposure numerator and leave `mic` on the total-drug scale the susceptibility method actually reports, since the free-drug index `fu * Cc` against `mic` is equivalent to total `Cc` against `mic / fu`. For models that carry more than one pathogen or isolate, suffix per analyte in the usual way (`mic_<organism>`). Ratified with the Chen 2023 tilmicosin extraction and extended with the Beredaki 2023 micafungin extraction; the two entries those left behind were merged 2026-08-06, since a repeated name at the same `Type` is resolved last-writes-win and the earlier one was being silently discarded.

### kpu<n> (**canonical clustered unbound tissue:plasma partition coefficient**)
- **Type:** paper-named-param
- **Role:** Estimated unbound tissue-to-plasma partition coefficient (Kpu) shared by a *cluster* of PBPK tissues, numbered `kpu1`, `kpu2`, `kpu3`, `kpu4`, ... (log-transform prefix `lkpu<n>`). Unitless. Used by "steady-state commonality" PBPK simplifications that keep the full whole-body kinetic structure but reduce the number of unknown distribution parameters by assigning one common Kpu to every tissue in a composition-derived cluster. The numeric suffix indexes the cluster, not a compartment; the cluster-to-tissue mapping must be written out explicitly in `model()` (e.g. `kpu_bone <- kpu1`) and cited to the source table footnote. Convert to the tissue:blood coefficient with `kb = kpu * fu_p / BP`.
- **Source aliases:**
  - `KPU1` .. `KPU4` -- NONMEM `$PK` names in Yau 2023 Appendix S1.
  - `Kpu1` .. `Kpu4` -- Yau 2023 Table 2 / Table S6.
- **Example models:** `Yau_2023_diazepam_pbpk_kpu_human.R`, `Yau_2023_diazepam_pbpk_kpu_rat.R`, `Yau_2023_diazepam_pbpk_lumped_human.R`, `Yau_2023_diazepam_pbpk_lumped_rat.R`.
- **Notes:** Distinct from the per-tissue `kp_<tissue>` family (a separate Kp value named for one specific organ, as in the Gaohua 2012 pregnancy PBPK models) -- `kpu<n>` is deliberately cluster-indexed because the whole point of the parameterisation is that several anatomically distinct tissues share one estimate. Also distinct from `sf<n>` below, which scales a *predicted* Kpu rather than replacing it. In kinetically lumped models the same `kpu<n>` names index the lumped compartments (`kpu1` = central lump, `kpu2` = peripheral1 lump, ...). Introduced 2026-07-26 with the Yau 2023 diazepam PBPK extraction.

### sf<n> (**canonical clustered Kpu scaling factor**)
- **Type:** paper-named-param
- **Role:** Estimated multiplicative correction factor applied to a *bottom-up predicted* unbound tissue:plasma partition coefficient, shared by a cluster of PBPK tissues and numbered `sf1`, `sf2`, `sf3`, `sf4`, ... (log-transform prefix `lsf<n>`). Unitless. Implements `Kpu_i = Kpu_pred,i * SF_cluster(i)`, where `Kpu_pred,i` comes from a mechanistic tissue-composition model (Rodgers and Rowland, Poulin and Theil, Berezhkovskiy, ...) evaluated inside `model()`. The scalar therefore *quantifies the systematic bias* of the bottom-up prediction instead of discarding it, which is what makes this parameterisation attractive for interspecies translation (the bias is assumed to be conserved across species while the composition inputs change).
- **Source aliases:**
  - `SF1` .. `SF4` -- Yau 2023 Table 2 / Table S6 and Eq 10.
- **Example models:** `Yau_2023_diazepam_pbpk_scalar_human.R`, `Yau_2023_diazepam_pbpk_scalar_rat.R`.
- **Notes:** A value of 1 means the bottom-up prediction needs no correction, so `sf<n>` estimates are directly interpretable as fold-bias. Do NOT reuse `sf<n>` for a generic "scaling factor" on a non-partition-coefficient parameter -- allometric exponents use `allo_<param>`, and covariate multipliers use the `e_<cov>_<param>` family. Paired with the `kpurr_<tissue>` derived `model()` quantities in the founding examples. Introduced 2026-07-26 with the Yau 2023 diazepam PBPK extraction.

### eh (**canonical hepatic extraction ratio**)
- **Type:** paper-named-param
- **Role:** Hepatic extraction ratio in the well-stirred liver model (`CL_H = FQ * eh`, `F_HEP = 1 - eh`). Unitless, bounded in [0, 1]. Use whenever a paper parameterises hepatic clearance through the extraction-ratio physiology rather than through `lcl` / `lcl_nonren` directly. The companion log-transform prefix is `logiteh` (logit-scale `ini()` typical value) because the linear-scale `eh` is bounded; logit-additive eta on `logiteh` keeps every individual `eh_i` inside the (0, 1) box.
- **Source aliases:** `EH` -- used in `Chan_2008_maraviroc.R` and Brussee 2018 mAb PBPK (the latter as a derived `model()`-block quantity rather than an `ini()` parameter, so the canonical applies to the `ini()` use case introduced by Chan 2008).
- **Example models:** `Chan_2008_maraviroc.R` (estimated `logiteh = logit(0.662)` with logit-additive IIV per Chan 2008 Eq 9; downstream `eh` enters `clh = fq * eh` and `fhep = 1 - eh`).
- **Notes:** Paired with the canonical compartment / pseudo-parameter `fq` (hepatic plasma flow, fixed at a literature value) and the canonical `clr` (renal clearance, often fixed) when the paper decomposes total CL into renal + hepatic with hepatic-extraction physiology. Distinct from `lcl_nonren` (additive renal + non-renal CL decomposition without an explicit extraction-ratio bound): use `eh` only when the source paper writes hepatic clearance as `CL_H = FQ * E_H` and constrains E_H in [0, 1] (e.g., via a logit-form IIV transformation).

### nowsmax (**canonical bare NAS-natural-history maximum-baseline term**)
- **Type:** paper-named-param
- **Role:** Multiplicative baseline of the natural neonatal-abstinence-syndrome (NAS) severity-decay-with-postnatal-age term `NOWST = nowsmax * exp(-nowsm * PNA_days)` used inside the Eudy-Byrne 2021 buprenorphine indirect-response PD model of MOTHER NAS scores. `nowsmax` is unitless (a multiplier on `kin/kout`) and enters both the ODE production term `kin * (1 + nowst)` and the drug-free quasi-steady-state initial condition `nows(0) = kin * (1 + nowst) / kout`. Log-transformed form is `lnowsmax`.
- **Source aliases:**
  - `NOWSMAX` -- Eudy-Byrne 2021 Table S3 notation (all-caps).
- **Example models:** `EudyByrne_2021_buprenorphine.R` (typical `nowsmax = 1.92`, unitless; 95% CI 1.76-2.08; log-normal IIV omega^2 = 1.14).
- **Notes:** Paired with `nowsm` (natural decay rate) and the canonical PD state / observation `nows`. Registered 2026-07-25 alongside the Eudy-Byrne 2021 extraction. Distinct from `rbase` (generic IDR baseline) because `nowsmax` multiplies an exponential decay with postnatal age rather than being a fixed steady-state baseline.

### lnowsmax (**canonical log-transformed NAS-natural-history maximum-baseline term**)
- **Type:** paper-named-param
- **Role:** Log-transformed form of `nowsmax` for `ini()` with log-normal IIV. Inside `model()` the bare name is `nowsmax = exp(lnowsmax + etalnowsmax)`.
- **Source aliases:**
  - `log NOWSMAX`, `lNOWSMAX` -- equivalent paper notation.
- **Example models:** `EudyByrne_2021_buprenorphine.R`.

### nowsm (**canonical bare NAS-natural-history decay rate with postnatal age**)
- **Type:** paper-named-param
- **Role:** First-order decay rate constant (1/day) of the natural NAS-severity term `NOWST = nowsmax * exp(-nowsm * PNA_days)` with chronological postnatal age. As PNA grows, NOWST decays toward zero and the drug-free quasi-steady-state NAS score `nows0 = kin * (1 + nowst) / kout` approaches its long-term floor `kin / kout`. Log-transformed form is `lnowsm`.
- **Source aliases:**
  - `NOWSM` -- Eudy-Byrne 2021 Table S3 notation (all-caps).
- **Example models:** `EudyByrne_2021_buprenorphine.R` (typical `nowsm = 0.107 1/day`; 95% CI 0.102-0.112; log-normal IIV omega^2 = 1.42; correlated with `nowsmax` at corr = 0.778).
- **Notes:** Units are 1/day. When paired with the canonical `PNA` covariate (which is in MONTHS), inside `model()` compute `pna_days = PNA * 30.4375` before use so the exponent stays dimensionless (same pattern as Zhao 2018).

### lnowsm (**canonical log-transformed NAS-natural-history decay rate**)
- **Type:** paper-named-param
- **Role:** Log-transformed form of `nowsm` for `ini()` with log-normal IIV. Inside `model()` the bare name is `nowsm = exp(lnowsm + etalnowsm)`.
- **Source aliases:**
  - `log NOWSM`, `lNOWSM` -- equivalent paper notation.
- **Example models:** `EudyByrne_2021_buprenorphine.R`.

### cal_slope_<assay> (**canonical assay cross-calibration slope**)
- **Type:** paper-named-param
- **Role:** Slope of a linear recalibration mapping one measurement modality's predictions onto another modality's scale, used when a single structural state is observed by two assays with different bias / gain and the paper estimates the mapping as part of the model: `pred_reference = cal_slope_<assay> * pred_<assay> + cal_int_<assay>`. Unitless. `<assay>` is a short lower-case token naming the *non-reference* modality (`spect`, `pet`, `dbs`, `saliva`, ...).
- **Source aliases:**
  - `beta` -- Siebinga 2023 Equation 1 (`Cpred = Cpred_SPECT * beta + alpha`).
- **Example models:** `Siebinga_2023_lu177psma617.R` (`cal_slope_spect` = 0.828, fixed; recalibrates SPECT/CT-derived blood activity onto the blood-sample scale).
- **Notes:** Registered 2026-07-30. Distinct from a covariate effect (`e_<cov>_<param>`): the calibration parameters describe the *measurement* process, not a biological source of variability, and belong with the residual-error block rather than with the structural PK. Typically fixed after estimation in a data-source-only sub-model.

### cal_int_<assay> (**canonical assay cross-calibration intercept**)
- **Type:** paper-named-param
- **Role:** Intercept of the linear recalibration described under `cal_slope_<assay>`, in the units of the observation.
- **Source aliases:**
  - `alpha` -- Siebinga 2023 Equation 1.
- **Example models:** `Siebinga_2023_lu177psma617.R` (`cal_int_spect` = 6.27 MBq/L, fixed).
- **Notes:** Registered 2026-07-30 alongside `cal_slope_<assay>`.

### cal_bias_<matrix> (**canonical structural measurement bias**)
- **Type:** paper-named-param
- **Role:** Additive structural measurement offset applied to the model prediction for a given sampling matrix before the residual-error model, in the units of the observation. Encodes a known / assumed systematic difference between the measured and the true quantity (assay calibration bias, matrix interference, background signal) that the paper estimates as a structural parameter rather than absorbing into the residual error. `<matrix>` names the sampling matrix (`blood`, `plasma`, `urine`, ...).
- **Source aliases:**
  - `gamma` -- Siebinga 2023 Equation 2 (`Cobs = (Cpred + gamma) * (1 + eps_p) + eps_add`).
- **Example models:** `Siebinga_2023_lu177psma617.R` (`cal_bias_blood` = 0.273 MBq/L, fixed; the paper attributes it to calibration uncertainty from extreme calibration ranges for blood samples).
- **Notes:** Registered 2026-07-30. A positive `cal_bias_<matrix>` raises the prediction, so it forces predictions above the drug-free baseline; Siebinga 2023 notes this is the source of the apparent under-prediction of low blood observations in its CWRES plots.

## Nested (multi-level) random effects

Registered 2026-08-06 with the `Qi_2024_vosoritide.R` extraction, the first
model in this library to carry a second hierarchical level of random effects.
This section documents a naming pattern, not a canonical parameter, so it has
no H3 entries and contributes nothing to the runtime name vectors.

When a paper fits a **second hierarchical level** of random effects on top of
the per-subject etas -- a between-study, between-site or between-centre random
effect whose draw is shared by every subject in that group -- name the nested
eta `eta<lparam>_<level>`, where `<lparam>` is the same transformed parameter
name the subject-level eta uses and `<level>` is a short lowercase token naming
the grouping level. Declare it in `ini()` with rxode2's native nesting syntax,
piping the grouping column after the variance:

```r
etalcl       ~ 0.112896            # subject-level IIV on CL/F
etalcl_study ~ 0.066049 | SIDN     # study-level (nested) random effect on CL/F
```

and add it to the same log-scale sum as the subject-level eta in `model()`:

```r
cl <- exp(lcl + etalcl + etalcl_study) * (WT / ref_wt)^e_wt_cl
```

Founding example: `Qi_2024_vosoritide.R` (`etalcl_study` and `etalvc_study`,
nested on the `SIDN` study-identity column; Qi 2024 Sect. 2.4 eta6 / eta7).

**The grouping column is a covariate register entry of its own** (see `SIDN` in
`covariate-columns.md`), not a parameter. Because it is consumed by the `ini()`
nesting syntax rather than by any `model()` expression,
`checkModelConventions()` reports `covariateData[['<level>']] has an entry but
is not referenced in model()`. That warning is structural for a nesting level.

**Distinct from inter-occasion variability.** IOV indexes occasions *within* a
subject (`OCC`, decomposed `etalcl_oc1..N` sharing one variance -- see
`Jonsson_2011_ethambutol.R`, `Chen_2023_nemonoxacin.R`); a nested level groups
*across* subjects, so every record of a subject shares one draw.

**Two rxode2 simulation constraints apply** (verified against rxode2 5.1.3).
The event table must contain at least two distinct grouping-column values --
with only one level `rxSolve()` fails with `The following parameter(s) are
required for solving: THETA[2], THETA[1]`, an error naming thetas the model does
not have and giving no hint that the nesting is the cause. And `omega` must be
passed explicitly (`omega = mod$omega`), because a nested omega is a *list* of
matrices keyed by level rather than a single matrix.

## Unit spellings

Registered 2026-08-03. The machine-readable `units` block wrote the same time
unit three ways -- `"hour"` in 643 models, `"h"` in 208, `"hr"` in 28 -- so a
consumer parsing `units$time` could not canonicalise without its own spelling
table. Same for `"minute"` vs `"min"` and `"microgram"` vs `"ug"`.

Canonical spellings live in `conventions$timeUnitSpellings` and
`conventions$doseUnitSpellings`; `checkModelConventions()` raises an
**error** for a non-canonical spelling and `buildModelDb()` aborts on it.

**Spelling is normalised; dimension is never converted.** `min` and `h` are
both canonical and are never conflated -- rewriting one as the other would
misstate every value in the model. The check only maps alternative spellings
of the *same* unit onto one form.

Generic structural models (`PK_1cmt`, `PK_2cmt`, ...) are dimensionless by
design and declare `"time_unit"` / `"dose_unit"`; `Beal_2001_iv1cmt_bql`
expresses time in half-lives. These are listed in
`conventions$placeholderUnits` and are exempt -- they are correct, not
unnormalised.

### Why `kon` is NOT canonicalised

A consumer asked for a canonical unit for "the kon concept". There isn't one,
because `kon*` is not one concept. Of the 91 parameters whose name starts
`kon` in this library:

- `kon51_col_atcc`, `kon52_col_aru` are an **EC50 in mg/L**, and
  `konca_col_atcc`, `koncp_col_aru` an **Emax in 1/h** -- `KON51` and `KONCA`
  are the source paper's own parameter names, not association rates.
- QSP 2D on-rates carry a **length** dimension:
  `1/(micromolarity*nanometer*second)`.
- Opioid receptor-binding models use `pM^-n s^-1`, with a Hill-like exponent.
- Some are mass-concentration based: `(mg/L)^-1 day^-1`.
- The remainder are ordinary 3D molar rates: `1/(nM*day)`, `1/nM/h`,
  `L/nmol/day`.

At least three distinct dimensionalities share the prefix. A blanket
canonicalisation would silently restate values -- the same trap as `lrbase`
(see #477), where 89 of 102 uses were indirect-response PD baselines rather
than target baselines. What *was* normalised in these labels is spelling and
separators only (`/hour` -> `/h`, `pM.day` -> `pM*day`, `mcmol` -> `umol`),
which cannot change a value.

A consumer that needs kon in one unit has to read the dimensionality per
model. That is a real cost, and it is smaller than the cost of being
confidently wrong.

## Permeability-limited whole-body PBPK parameters (Gaohua 2023)

Vocabulary for permeability-limited whole-body PBPK models, in which each tissue is resolved into residual blood cells, residual plasma, extracellular water and intracellular water (compartment suffixes `_bc` / `_plasma` / `_ew` / `_iw`; see `compartment-names.md`). Two conventions define the family. First, system physiology is carried as *fractions* -- `fvol_<organ>` of body weight and `fq_<organ>` of cardiac output -- so that a published table of percentages transcribes directly and both columns can be checked to sum to 100%. Second, every permeability-surface-area product and every clearance is expressed as a **fold-multiple of the compartment's own plasma flow** (`fold_*`), which is how these papers parameterise their what-if scenarios: one multiplier per mechanism, swept around the tissue blood flow. Registered per the operator naming decision of 2026-08-05.

### bw (**canonical body weight system constant**)
- **Type:** paper-named-param
- **Role:** Reference body weight (kg) of the simulated individual. In whole-body PBPK it is the normalising constant for `fvol_<organ>` and the denominator of the body-weight-normalised volume of distribution `Vdt`.
- **Source aliases:** `BW`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### qc (**canonical cardiac output**)
- **Type:** paper-named-param
- **Role:** Cardiac output (L/h); the normalising constant for `fq_<organ>`. Splits into a total plasma flow `(1 - hct) * qc` and a total blood-cell flow `hct * qc` when a model resolves the two blood phases separately.
- **Source aliases:** `Qc`, `CO`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### hct (**canonical haematocrit**)
- **Type:** paper-named-param
- **Role:** Haematocrit: the volume fraction of blood that is blood cells. Splits both blood volumes and blood flows into their blood-cell and plasma phases, and splits a tissue's residual blood volume into `<organ>_bc` and `<organ>_plasma`.
- **Source aliases:** `Haematocrit`, `HCT`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### density (**canonical tissue density**)
- **Type:** paper-named-param
- **Role:** Tissue and blood density (kg/L), converting an organ's share of body weight into a volume. Commonly assumed to be 1 kg/L, which makes tissue volumes in L numerically equal to their share of body weight in kg.
- **Source aliases:** `Density`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### fvol_<organ> (**canonical fractional organ volume**)
- **Type:** paper-named-param
- **Role:** Volume of the named organ as a fraction of body weight (L per kg), so `v_<organ> = fvol_<organ> * bw / density`. One entry per organ, plus `fvol_blood` for total blood; the set is expected to sum to 1. Distinct from the `vp_<organ>` / `v_<organ>` *absolute* volumes, which are the derived quantities.
- **Source aliases:** `% of total body weight`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### fq_<organ> (**canonical fractional organ blood flow**)
- **Type:** paper-named-param
- **Role:** Blood flow to the named organ as a fraction of cardiac output, so `qplasma_<organ> = fq_<organ> * (1 - hct) * qc`. One entry per organ; the set covering the organs that drain directly to venous blood is expected to sum to 1. Where an organ has a dual blood supply the arterial share is named explicitly (e.g. `fq_liver_arterial`) and the portal share is derived from the splanchnic organs it collects.
- **Source aliases:** `% of total cardiac output`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### frb_<organ> (**canonical residual blood volume fraction**)
- **Type:** paper-named-param
- **Role:** Volume of blood remaining in the named organ after bleeding, as a fraction of that organ's volume. Split by haematocrit into the `<organ>_bc` and `<organ>_plasma` subcompartment volumes. Conventionally sourced from rat tissue-residual-blood measurements.
- **Source aliases:** `residual blood`, `Vre`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### few_<organ> (**canonical extracellular water volume fraction**)
- **Type:** paper-named-param
- **Role:** Extracellular water of the named organ as a fraction of that organ's volume, giving the `<organ>_ew` volume. The intracellular water volume `<organ>_iw` is then usually taken as the remainder of the organ volume once residual blood and extracellular water are removed.
- **Source aliases:** `EW fraction`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### ps_bc_plasma (**canonical blood-cell / plasma permeability-surface area product**)
- **Type:** paper-named-param
- **Role:** Passive permeability-surface area product (L/h) between a compartment's residual blood cells and its residual plasma. Set to a deliberately large value when the source model assumes immediate equilibrium between the two blood phases, which collapses that pair into a single well-stirred space.
- **Source aliases:** `PStc/tp`, `PSrbc2plasma`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### fold_ps_plasma_ew (**canonical vascular-membrane permeability fold multiplier**)
- **Type:** paper-named-param
- **Role:** Passive permeability-surface area product between residual plasma and extracellular water, expressed as a fold-multiple of the compartment's plasma flow (the vascular membrane). Paired with `fold_ps_ew_iw`; when the two differ, the smaller one is rate-limiting, which is what makes a distribution metric saturate as the other is increased.
- **Source aliases:** `FoldPlasmaEw`, `PS/Q`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### fold_ps_ew_iw (**canonical cell-membrane permeability fold multiplier**)
- **Type:** paper-named-param
- **Role:** Passive permeability-surface area product between extracellular and intracellular water, expressed as a fold-multiple of the compartment's plasma flow (the cell membrane). For a typical small molecule the cell membrane is the less permeable of the two, so this is the multiplier a permeability sweep usually varies.
- **Source aliases:** `FoldEwIw`, `PS/Q`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### fold_clint_iw_<compartment> (**canonical intracellular water metabolic-clearance fold multiplier**)
- **Type:** paper-named-param
- **Role:** Metabolic (intrinsic) clearance in the intracellular water of the named compartment, as a fold-multiple of that compartment's plasma flow. Resolved per compartment rather than globally because these papers apply metabolism both to every tissue at once and to a single organ (e.g. the liver) in isolation; setting every compartment to one value recovers the single global multiplier.
- **Source aliases:** `FoldClearanceIW`, `CL/Q`, `CLmet`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### fold_clint_ew_<compartment> (**canonical extracellular water metabolic-clearance fold multiplier**)
- **Type:** paper-named-param
- **Role:** Metabolic (intrinsic) clearance in the extracellular water of the named compartment, as a fold-multiple of that compartment's plasma flow. Resolved per compartment rather than globally because these papers apply metabolism both to every tissue at once and to a single organ (e.g. the liver) in isolation; setting every compartment to one value recovers the single global multiplier.
- **Source aliases:** `FoldClearanceEW`, `CL/Q`, `CLmet`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### fold_clint_plasma_<compartment> (**canonical residual plasma metabolic-clearance fold multiplier**)
- **Type:** paper-named-param
- **Role:** Metabolic (intrinsic) clearance in the residual plasma of the named compartment, as a fold-multiple of that compartment's plasma flow. Resolved per compartment rather than globally because these papers apply metabolism both to every tissue at once and to a single organ (e.g. the liver) in isolation; setting every compartment to one value recovers the single global multiplier.
- **Source aliases:** `FoldClearancePLASMA`, `CL/Q`, `CLmet`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### fold_clint_bc_<compartment> (**canonical residual blood cells metabolic-clearance fold multiplier**)
- **Type:** paper-named-param
- **Role:** Metabolic (intrinsic) clearance in the residual blood cells of the named compartment, as a fold-multiple of that compartment's blood-cell flow. Resolved per compartment rather than globally because these papers apply metabolism both to every tissue at once and to a single organ (e.g. the liver) in isolation; setting every compartment to one value recovers the single global multiplier.
- **Source aliases:** `FoldClearanceRBC`, `CL/Q`, `CLmet`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### fold_uptake (**canonical active uptake transporter fold multiplier**)
- **Type:** paper-named-param
- **Role:** Active uptake transporter clearance carrying drug from extracellular into intracellular water, as a fold-multiple of the compartment's plasma flow. Paired with `fold_efflux`; uptake raises the tissue/plasma partition coefficient and the volume of distribution, efflux lowers both.
- **Source aliases:** `FoldUptake`, `CLew/iw`, `Cluptake`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### fold_efflux (**canonical active efflux transporter fold multiplier**)
- **Type:** paper-named-param
- **Role:** Active efflux transporter clearance carrying drug from intracellular back to extracellular water, as a fold-multiple of the compartment's plasma flow. Paired with `fold_uptake`.
- **Source aliases:** `FoldEfflux`, `CLiw/ew`, `Clefflux`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`

### Kp_<organ> (**canonical dynamic tissue/plasma partition coefficient output**)
- **Type:** paper-named-param
- **Role:** Time-varying tissue/plasma partition coefficient of the named organ: the whole-tissue concentration (total amount across all four subcompartments divided by the total organ volume) divided by the venous plasma concentration. An algebraic **output** of a permeability-limited model, not an input -- which is what distinguishes it from the perfusion-limited `kp_<organ>` input parameter, where the partition coefficient is a fixed constant supplied to the model.
- **Source aliases:** `Kp(t)`, `Kpss`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`
- **Notes:** Capitalised `Kp_` deliberately, to keep the dynamic output visually distinct from the lowercase `kp_<organ>` perfusion-limited input parameter used by e.g. `Gaohua_2012_pregnancy_pbpk_midazolam.R`.

### Vdt (**canonical dynamic volume of distribution output**)
- **Type:** paper-named-param
- **Role:** Time-varying volume of distribution (L/kg): the total drug amount in blood and tissues divided by the venous plasma concentration and body weight. An algebraic output; its pseudo-steady-state value is the reported Vdss. Distinct from the registered `vd` (apparent volume of distribution as a fitted structural parameter) -- `Vdt` is a model prediction that changes over time, `vd` is an input.
- **Source aliases:** `Vd(t)`, `Vdss`
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`
- **Notes:** Named `Vdt` rather than `vd` precisely so that the two never collide in one model.
