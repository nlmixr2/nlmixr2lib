# Canonical parameter names

This file is the authoritative register of structural and paper-mechanistic parameter names used in nlmixr2lib models. Every parameter referenced inside a model's `ini` and `model` blocks is expected to match one of the canonical names below (or, for covariate effects, the documented `e_<cov>_<param>` shape — see `R/conventions.R::covEffectPattern`).The register is extended extended whenever a new paper introduces a parameter that isn't yet registered.

## How to use this register

1. **Before adding a parameter to a new model**, search this file (by canonical name and by source alias) for the concept you need.
2. **If the canonical name exists**, use it exactly. Document any source-paper rename in a code comment or in the model's `description` block.
3. **If the source paper uses an alias listed under an existing canonical name**, prefer the canonical name. Aliases are documented for cross-reference, not as a free pass to introduce the deprecated form in new models.
4. **If the parameter is not in this register at all**, propose a new entry with a canonical name, type, role, source aliases, and example models. Verify with the user before committing. The addition is part of the model's PR.
5. **Do not modify existing model files when you discover a missing entry**; simply register the canonical here. Retrofitting existing models is a separate effort.
6. **Never add a second entry for a name that already has one at the same `Type`.** Extend the existing block instead -- add your source alias and example model to it. A name may appear twice only when the two entries carry *different* `Type` values. This register is resolved in document order, last one wins, so a same-`Type` repeat silently discards the earlier block along with any alias or example recorded only there. `buildModelDb` fails the build on one.

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
  - `e_wt_<param>` is the body-weight member of that family and is the canonical spelling for an **allometric exponent**. The deprecated `allo_<param>` / `allovc` / `dCLdWT` / `dVdWT` spellings are rewritten onto it by the rename map in `R/conventions.R` (`allo_vc` -> `e_wt_vc`, `allo_cl` -> `e_wt_cl`, ...), so `checkModelConventions()` will reject `allo_vc` and name `e_wt_vc` as the fix.
  - When a paper states a **single** exponent that it applies to every rate constant at once, rather than one exponent per named parameter, use the generic token `e_wt_k` for that shared exponent. Founding examples: `Biliouris_2018_nusinersen.R` (`e_wt_k` = -0.08, applied to all rate constants) and `Siebinga_2024_lu177psmaIT.R` (`e_wt_k` = -0.25 across all nine rate constants, alongside `e_wt_vc` = 1). Ratified 2026-08-05 with the Siebinga 2024 extraction. Prefer a specific `e_wt_<param>` whenever the paper does report per-parameter exponents.
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

The `l<base>` convention denotes a population mean estimated on the log scale (`ini` carries `l<base> <- log(typical_value)`; `model` declares `<base> <- exp(l<base> + eta_<base>)`). Use the log-transformed form whenever the parameter is strictly positive and the source paper reports an exponential typical-value form.

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
- **Notes:** Member of the `lv<compartment>` family that names a volume after the canonical compartment it scales, following the existing `lvcsf` / `lvcns` precedent in `Luu_2017_nusinersen.R`. Distinct from `lvc`: `lvelf` never scales a plasma concentration. `Parmar_2023_spectinamide_1599_mouse_pbpk.R` carries the same quantity as a hardcoded physiologic constant `v_elf` inside `model` rather than as an `ini` parameter, which is correct for a PBPK model where the volume is anatomy rather than a fitted parameter.

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

### lffo (**canonical log-transformed first-order absorption-site fraction**)
- **Type:** log-transformed-pk
- **Role:** Log-scale parameter for the fraction of the bioavailable dose routed to the fast (first-order, `ka1`) site of a parallel two-site absorption model, the remainder going to the slow site.
- **Source aliases:**
  - `Frapid` -- used in `Winter_2024_oxytetracycline_cattle.R` (Winter 2024 Table 1 `tvFrapid`).
- **Example models:** `Winter_2024_oxytetracycline_cattle.R`.
- **Notes:** The log-scale sibling of the existing `logitffo` canonical (`Baverel_2015_tralokinumab.R`, `Cirincione_2017_exenatide.R`, `Lacy_2018_cabozantinib.R`, `PerezRuixo_2008_epoetinAlfa.R`), sharing the same `ffo` bare root. Use `logitffo` whenever the source paper estimates the fraction on the logit scale, which is the safer encoding because it cannot leak above 1. Use `lffo` only when the source explicitly parameterises the fraction multiplicatively with exponential IIV, as Winter 2024 does (`stparm(Frapid = tvFrapid * exp(nFrapid))` in its Phoenix control stream); in that case the published parameterisation genuinely admits values above 1 in the tails, and re-encoding it on the logit scale would silently change the model.

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
- **Notes:** Analyte- or site-suffixed variants (`lbmax_c`, `lbmax_p`, `lbmax_rbc`) are permitted on the standard suffix pattern when a model carries more than one binding pool.

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

### ltclchange (**canonical log-transformed clearance-step breakpoint time**)
- **Type:** log-transformed-pk
- **Role:** Log-scale time, measured from treatment initiation, at which a piecewise-constant ("step") time-varying clearance switches from its early arm to its late arm (time).
- **Source aliases:**
  - `TCLchange` -- Park 2025.
  - `tNab` -- Yoneyama 2017 (an onset time gating a covariate effect on CL/F, the same MTIME construction).
- **Example models:** `Park_2025_efineptakin_alfa.R`.
- **Notes:** For the piecewise-constant form `cl <- cl_early * (t < tclchange) + cl_late * (t >= tclchange)`, which is a switch between two mutually exclusive arms and NOT the additive decomposition that `lcl_ss` / `lcl_time` describe. The naming of that form is deliberately asymmetric: **the early arm is carried by the plain canonical `lcl`**, so there is no `lcl_early`, and a covariate acting before the breakpoint takes the ordinary two-token shape (`e_sexf_cl`). Only the post-breakpoint arm needs its own name (`lcl_late`). In NONMEM this structure is written with MTIME. Use `lcl_t50` instead when the time-course is a smooth sigmoid rather than a step.

### lcl_late (**canonical log-transformed post-breakpoint clearance arm**)
- **Type:** log-transformed-pk
- **Role:** Clearance in force after the breakpoint of a piecewise-constant time-varying clearance (volume / time). Paired with `ltclchange`; the pre-breakpoint arm is the plain `lcl`.
- **Source aliases:**
  - `CL>TCLchange` -- Park 2025.
- **Example models:** `Park_2025_efineptakin_alfa.R`.
- **Notes:** Only the LATE arm carries a suffix, because the early arm is an ordinary clearance already covered by `lcl`. Covariates may act on the two arms independently -- in the founding example sex acts on `lcl` only and the late arm is sex-independent.

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

### lcl_tsnet (**canonical log-transformed net-tubular-secretion clearance arm**)
- **Type:** log-transformed-pk
- **Role:** Net renal tubular-secretion component of a renal clearance decomposition `CL_renal_total = GFR + CL_tsnet`, where the filtration arm is carried by the model's glomerular-filtration parameter (in the founding model, the co-administered iohexol probe's `cl`) and this arm carries the transporter-mediated secretory flux *net of tubular reabsorption*. The `tsnet` token keeps the NET qualifier deliberately: the quantity published under this name is secretion minus reabsorption, so a later model that resolves the two separately still has both gross-secretion and reabsorption names available. Applied with an analyte suffix when the model tracks more than one species (`lcl_tsnet_creatinine`).
- **Source aliases:**
  - `nCTS` -- Chen 2025 paper notation ("net creatinine tubular secretion", Table 3 row `nCTS (mL/min)` = 39.7).
  - `CL_SEC` -- Chen 2025 NONMEM control-stream `$PK` symbol (`TVCL_SEC = THETA(2)`).
- **Example models:** `Chen_2025_iohexol_creatinine.R` (founding example; `lcl_tsnet_creatinine = log(2.38177)` L/h, the OCT2/MATE-mediated secretory arm of creatinine clearance, ~31% of total CrCL in healthy adults).
- **Notes:** Joins the established `lcl_*` clearance-arm family (`lcl_ss`, `lcl_time`, `lcl_renal`, `lcl_nonren`, `lcl_hemodialysis`, `lcl_met`, `lcl_nonmet`, `lcl_form`, `lcl_ugt`) rather than being registered under the paper's own `ncts` acronym: net tubular secretion is another arm of the same additive clearance decomposition, and `nCTS` has "creatinine" baked into the acronym, so it could never generalise to the other OCT2/MATE and OAT substrates (metformin, cimetidine, trimethoprim) whose papers report the same quantity. Distinct from `lcl_renal`, which is a *lumped* renal arm opposed to `lcl_nonren`; `lcl_tsnet` decomposes *within* the renal arm and is paired with a filtration term, not with a non-renal term. The token is also registered in `R/conventions.R::clComponents` so covariate effects classify as `e_<cov>_cl_tsnet`.

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
- **Example models:** `Knibbe_2009_morphine.R` (PNA-stratified `lcl_form_m3g_le10` / `lcl_form_m3g_gt10` and `lcl_form_m6g_le10` / `lcl_form_m6g_gt10` for the morphine -> M3G and morphine -> M6G glucuronidation arms), `Franken_2015_morphine.R` (`fm_m3g * cl`, `fm_m6g * cl` derived within `model` from fixed fractions), `deHoogd_2017_morphine.R`, `Hennig_2015_rifabutin.R`, `Djerada_2014_nefopam.R` (`lcl_form_dnef` = log(K13); the apparent metabolic clearance of nefopam to desmethyl-nefopam).
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
- **Role:** Duration of the zero-order input rate into the depot (or central) compartment for sequential zero+first-order or pure zero-order absorption models. Applied via `dur(depot) <- d1` or `dur(central) <- d1` inside `model`; the dose mass is delivered uniformly over `[0, d1]`.
- **Source aliases:**
  - `D1` -- common NONMEM / paper notation (Aouri 2017, Goggin 2004, Heathman 2024, Gupta 2016).
  - `lduration` -- Gupta 2016 paper notation; an equivalent log-transformed entry.
- **Example models:** `Aouri_2017_rilpivirine.R` (pure zero-order absorption directly into central; D1 = 4 h), `Goggin_2004_emfilermin.R` (zero-order SC absorption, D1 = 0.84 h), `Heathman_2024_efavirenz.R` (sequential zero-order then first-order; D1 = 1.74 h followed by KA = 0.165/h).

### lkel (**canonical log-transformed elimination rate constant (K-PD)**)
- **Type:** log-transformed-pk
- **Role:** First-order elimination rate constant, in either of two settings: (a) K-PD or single-rate-constant elimination forms where no explicit `vc` is estimated, and (b) rate-constant-parameterised disposition models that DO estimate an explicit `lvc` alongside it, i.e. the source paper fits `ke` and `V/F` as two separate parameters rather than fitting a clearance. Case (b) must not be silently reparameterised to `lcl` + `lvc`: when the paper places independent random effects on `ke` and on `V`, the algebraically equivalent clearance form would require `etalcl = etalkel + etalvc`, a correlated eta block the authors did not fit, so the reparameterisation would change the model. This mirrors the register's existing acceptance of rate-constant parameterisations at `k12` / `k21`.
- **Source aliases:**
  - `lke` -- legacy name; replaced 2026-05-28 by the naming audit.
  - `lkde` -- paper-named (Mazzocco 2015 / Shoji 2017 KDE) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
  - `lkp` -- paper-named (van Hasselt 2015 KP) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
- **Example models:** `Mazzocco_2015_temozolomide.R`, `Shoji_2017_fosdagrocorat_oc.R`, `Shoji_2017_fosdagrocorat_p1np.R`, `vanHasselt_2015_eribulin.R` (all case (a), no explicit `vc`); `Baklouti_2026_amoxicillin.R` (case (b), `ke = 0.31 1/h` estimated alongside `V/F = 82.38 L`, each with its own eta).
- **Notes:** Canonical `lkel` adopted 2026-05-28 per the naming audit. Role broadened 2026-09-02 to cover rate-constant-parameterised disposition models that also estimate an explicit volume (operator sidecar `oare_PMC13206287` request-001 / response-001, question q2, option A); the entry already carried the paper symbol `ke` as an alias, so registering a separate `lk10` would have duplicated it.

### lkel_exp_kdes (**canonical log-transformed decay-rate constant of a time-varying elimination rate constant**)
- **Type:** log-transformed-pk
- **Role:** Log-scale rate constant governing how fast a time-varying elimination rate constant relaxes from its value at t = 0 toward its asymptote, in the exponential form `kel(t) = kel * (1 + kel_exp_famp * (1 - exp(-kel_exp_kdes * t)))` (1 / time). `log(2) / kel_exp_kdes` is the half-life of the change.
- **Source aliases:**
  - `ln(kappa_k)` -- used in `Schreib_2024_busulfan.R` (paper `theta_kappa_k`).
- **Example models:** `Schreib_2024_busulfan.R` (`lkel_exp_kdes = -2.965`, giving 0.0516 1/h and a 13.4 h half-life for the decline in busulfan `k` and `CL` over a conditioning course).
- **Notes:** The `kel` counterpart of the registered `cl_exp_kdes` role in the time-varying-clearance family; reuses the `_kdes` role token deliberately so both families are found by the same stem.

### ltlag (**canonical log-transformed absorption lag time**)
- **Type:** log-transformed-pk
- **Role:** Log-scale absorption lag time before drug enters the depot (time).
- **Source aliases:**
  - `lalag` -- legacy.
  - `llag` -- legacy.
- **Example models:** delayed-absorption oral PK models.
- **Notes:** Replaces the legacy `lalag` / `llag` forms.

### lka_early (**canonical log-transformed absorption rate constant before a sequential-absorption switch**)
- **Type:** log-transformed-pk
- **Role:** Log-scale first-order absorption rate constant applying from the dose until the sequential-absorption switch time `tkacut` (1/time).
- **Source aliases:**
  - `Ka1` -- used in `vanSchaick_2016_prucalopride_pediatric.R` (van Schaick 2016 Table 3; the paper's "absorption rate at time less than MTIME").
- **Example models:** `vanSchaick_2016_prucalopride_pediatric.R`.
- **Notes:** Paired with `lka_late` and `tkacut` for papers describing absorption as two *sequential* first-order processes. The window suffix names the window, not the paper's ordinal symbol: `lka1` / `lka2` are deliberately NOT canonical here, because a bare ordinal encodes only sequence and not which time window the rate belongs to, and `lka2` already denotes an IV-conversion rate in `Kumpulainen_2010_flurbiprofen.R`. When a paper reports three or more stages whose boundaries are study-design windows rather than a reported parameter, prefer the boundary-labelled form used by `Othman_2007_carvedilol.R` (`lka_cr_0to2` / `lka_cr_2to4` / `lka_cr_gt4`, boundaries as literals in `model()`) and omit `tkacut`.

### lka_late (**canonical log-transformed absorption rate constant after a sequential-absorption switch**)
- **Type:** log-transformed-pk
- **Role:** Log-scale first-order absorption rate constant applying after the sequential-absorption switch time `tkacut` (1/time).
- **Source aliases:**
  - `Ka2` -- used in `vanSchaick_2016_prucalopride_pediatric.R` (van Schaick 2016 Table 3; the paper's "absorption rate at time greater than MTIME").
- **Example models:** `vanSchaick_2016_prucalopride_pediatric.R`.
- **Notes:** See `lka_early` for the family rationale and for when to prefer the boundary-labelled `Othman_2007_carvedilol.R` form instead.

### tkacut (**canonical sequential-absorption switch time**)
- **Type:** paper-named-param
- **Role:** Time after dose at which absorption switches from one first-order rate to the next, in the model's time unit (time).
- **Source aliases:**
  - `MTIME` -- used in `vanSchaick_2016_prucalopride_pediatric.R` (van Schaick 2016 Table 3, "cut-off time between first and second absorption rate"; the symbol comes from the NONMEM model-event-time mechanism used to implement the switch).
- **Example models:** `vanSchaick_2016_prucalopride_pediatric.R` (`tkacut <- fixed(0.734)` h, with `lka_early` = log(0.792 /h) before and `lka_late` = log(3.87 /h) after).
- **Notes:** Carried on the LINEAR scale, so a held-constant value is written `tkacut <- fixed(0.734)`; use `ltkacut` if a future paper estimates the switch time on the log scale. Named `tkacut` rather than `mtime` to avoid colliding with rxode2's `mtime()` model-event-time function. Apply inside `model()` by assigning the early rate first and letting the late window over-write it, switching on time after dose so the sequence resets at every administration: `ka <- exp(lka_early + etalka_early); if (tad(depot) >= tkacut) ka <- exp(lka_late + etalka_late)`. Note that when the two windows carry different IIVs the conditional assignment is not mu-referenced and rxode2 warns accordingly; that is expected and affects re-estimation, not simulation. Reach for `tkacut` only when the switch time is itself a reported parameter -- for study-design window boundaries use literals as in `Othman_2007_carvedilol.R`.

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
- **Role:** Log of the maximum / asymptotic first-order absorption rate constant in a Piotrovskij-style time-dependent (Weibull) absorption model `ka(t) = kamax * (1 - exp(-(ra * tad)^gam1))` (1 / time). The bare counterpart inside `model` is `kamax`.
- **Source aliases:**
  - `KAMAX` -- NONMEM `$THETA` convention used in NONMEM 7-era Weibull-absorption control streams.
- **Example models:** `Desai_2016_isavuconazole.R` (founding example).
- **Notes:** Distinct from `lka` (the single first-order rate constant of the simple zero-shape absorption model). Used together with `lra` and `lgam1` to parameterize Weibull / Piotrovskij time-dependent absorption.

### lra (**canonical log-transformed Weibull-absorption rate-scaling parameter**)
- **Type:** log-transformed-pk
- **Role:** Log of the rate-scaling parameter inside a Weibull absorption / release function (1 / time). Serves two related structures. (1) The Piotrovskij-style time-dependent ka, `ka(t) = kamax * (1 - exp(-(ra * tad)^gam1))`, where the product `(ra * tad)^gam1` drives the time-dependence of ka and larger `ra` shifts the ka rise earlier. (2) A prescribed Weibull *release* input function of the form `S(t) = exp(-(ra * t)^gam1)` giving the fraction of a depot dose still unreleased, whose hazard `h(t) = gam1 * ra * (ra * t)^(gam1 - 1)` empties a depot compartment. Both are the same one-parameter Weibull time scale. The bare counterpart inside `model` is `ra`.
- **Source aliases:**
  - `RA` -- NONMEM `$THETA` convention used in Weibull-absorption control streams.
  - `TD` -- release-function papers commonly report the **time to release 63.2%** of the dose rather than the rate scaler; the two are exact reciprocals, `ra = 1 / TD`, so a source reporting `TD = 117 h` becomes `lra <- log(1/117)`. Always spell the reciprocal out in the in-file comment, because the `ini` value is then not the number printed in the source table.
- **Example models:** `Desai_2016_isavuconazole.R` (founding example, saturating-ka form), `Perlstein_2026_olanzapine_lai.R` (release-input form; `TD = 117 h` released as `ra = 1/117`).
- **Notes:** Distinct from `lka` (simple first-order absorption) and from any infusion-rate parameter. Used together with `lkamax` and `lgam1` in the saturating-ka form, and with `lgam1` alone (plus the second-process partners `lra2` / `lgam2`) in the release-input form.

### lgam1 (**canonical log-transformed Weibull-absorption shape parameter**)
- **Type:** log-transformed-pk
- **Role:** Log of the unitless Weibull shape (sigmoidicity) parameter paired with `lra`, in either the Piotrovskij-style time-dependent ka or the Weibull release-input function (see `lra` for both forms). Larger values make the rise more abrupt; `gam1 = 1` reduces the Weibull form to a simple saturating exponential in the ka form and to first-order release in the release form. The bare counterpart inside `model` is `gam1`.
- **Source aliases:**
  - `GAM1` / `GAMMA1` -- NONMEM `$THETA` convention used in Weibull-absorption control streams.
  - `SS` -- "sigmoidicity factor" in release-function papers (e.g. Perlstein 2026 Table 1). Note that `ss` on its own is a **reserved rxode2 symbol** (steady state) and cannot be used as a parameter or state name -- `rxode2` hard-errors with `'ss' cannot be a state or lhs expression` -- which is one reason a paper's `SS` / `SS1` notation must be renamed to `gam1` / `gam2`.
- **Example models:** `Desai_2016_isavuconazole.R` (founding example, saturating-ka form), `Perlstein_2026_olanzapine_lai.R` (release-input form; `SS = 1.4`).
- **Notes:** Distinct from `lhill` (sigmoidal Emax / Imax exponent) and from `lgamma` (Friberg myelosuppression feedback / TGI growth exponents). The `gam1` suffix follows the NONMEM convention for Weibull-absorption sigmoidicity.

### lra2 (**canonical log-transformed second-process Weibull release rate-scaling parameter**)
- **Type:** log-transformed-pk
- **Role:** Log of the rate-scaling parameter of the **second** Weibull process in a multi-phase release input function (1 / time), the direct partner of `lra`. Used for long-acting-injectable depots whose in vivo release is described as a mixture of two parallel Weibull processes, `r(t) = frel * exp(-(ra * t)^gam1) + (1 - frel) * exp(-(ra2 * t)^gam2)`, where `r(t)` is the fraction of the dose still unreleased. The second process is typically the slower, more sigmoidal, sustained-release phase. The bare counterpart inside `model` is `ra2`.
- **Source aliases:**
  - `TD1` / `TD1_0` -- time to release 63.2% of the second process's dose fraction, the reciprocal of `ra2` (Perlstein 2026 Table 1). As with `lra`, spell the reciprocal out in the in-file comment.
  - `RA2` -- the natural NONMEM `$THETA` spelling when a control stream numbers the two processes.
- **Example models:** `Perlstein_2026_olanzapine_lai.R` (founding example; `TD1_0 = 323 h`, encoded as `lra2 <- log(1/323)`, with a linear dose effect `e_dose_ra2 = 0.0939 h/mg` applied on the reciprocal so the paper's printed additive form `TD1 = TD1_0 + TD1_1 * DOSE` is preserved exactly).
- **Notes:** The operator's ruling was explicit: a double-Weibull release input is the existing Weibull machinery applied twice, so it takes the registered `ra` / `gam1` stems plus numbered second-process partners, rather than a parallel `trel` / `ssrel` family. A third or later process would continue the numbering (`lra3` / `lgam3`). Distinct from `lkamax`, which has no analogue in the release-input form: a release mixture has no asymptotic rate constant, only a fraction split. Prior art with paper-local names: `Bonner_2015_gastric_emptying.R` implements the same double-Weibull mixture as `lgamma1` / `lbeta1` / `lgamma2` / `lbeta2` / `llogit_pr`; it is not retrofitted, per the register's standing rule that existing files are not migrated.

### lgam2 (**canonical log-transformed second-process Weibull release shape parameter**)
- **Type:** log-transformed-pk
- **Role:** Log of the unitless Weibull shape (sigmoidicity) parameter of the **second** release process, the direct partner of `lgam1` in the double-Weibull release input function described under `lra2`. The bare counterpart inside `model` is `gam2`.
- **Source aliases:**
  - `SS1` -- "sigmoidicity factor for the second process" (Perlstein 2026 Table 1). See the `lgam1` entry for why a paper's `SS` / `SS1` notation cannot be carried across literally.
  - `GAM2` / `GAMMA2` -- NONMEM `$THETA` convention.
- **Example models:** `Perlstein_2026_olanzapine_lai.R` (founding example; `SS1 = 3.25`, a markedly more sigmoidal second phase than the first process's `SS = 1.4`).
- **Notes:** Note the deliberate asymmetry in the numbering: the FIRST process keeps the unsuffixed-except-for-`1` NONMEM name `gam1` (which predates this entry) while its rate partner is the unnumbered `ra`, so a two-process release model reads `ra` / `gam1` then `ra2` / `gam2`. This is inherited from the Desai 2016 registration and is not worth churning existing files over.

### logitfrel (**canonical logit-transformed release-process fraction**)
- **Type:** log-transformed-pk
- **Role:** Logit-scale encoding of the fraction of a dose entering the **first** process of a multi-phase release input function (unitless, bounded in (0, 1)); the remaining `1 - frel` enters the second process. Inside `model` the bare form is `frel = 1 / (1 + exp(-logitfrel_ind))` where `logitfrel_ind` collects the fixed effect, covariate shifts, and IIV on the logit scale -- so IIV is additive on the logit, not multiplicative. Collect the fixed effect and eta on their own line (`logitfrel_ind <- logitfrel + etalogitfrel`) so the term stays in a mu-referenced position; rxode2 otherwise warns that the eta defaulted to non-mu-referenced.
- **Source aliases:**
  - `FF` -- "fraction of the available dose released in the first process" (Perlstein 2026 Methods p. 4 and Table 1).
- **Example models:** `Perlstein_2026_olanzapine_lai.R` (founding example; `FF = 0.236`, encoded as `logitfrel <- log(0.236 / (1 - 0.236))`, with the rapid first process carrying 23.6% of the dose and the sustained second process the remaining 76.4%).
- **Notes:** The operator directed the **logit** scale specifically, rejecting a log-scale `lfrel`, for the reason already recorded under `logitfm` and `logitfdepot`: a log-scale encoding of a bounded quantity can leak above 1 under moderate eta or covariate values. Follows the `logit`-transform-prefix family (`logitffo`, `logitemax`, `logitfm`, `logitfdepot`). Distinct from `logitfdepot` / `lfdepot`, which is a **bioavailability** -- the fraction of the administered dose that reaches the depot at all -- whereas `logitfrel` splits an already-bioavailable dose between two parallel release kinetics; a model can legitimately carry both. Also distinct from `fr` ("fraction of MAT in transit delay"). Three or more processes would take `logitfrel2` etc., but a softmax / stick-breaking encoding should be considered at that point.

### lbeta_cl (**canonical log-transformed exponential-nonlinear-CL slope**)
- **Type:** log-transformed-pk
- **Role:** Log of the concentration-slope coefficient in an exponential-nonlinear-clearance function of the form `CL(C) = CL_0 * exp(beta_cl * C)`, where `CL_0` is the linear-scale clearance at C = 0 (encoded as the standard `lcl` parameter) and `C` is the observed drug concentration in the same units as the paper's Table. Used when the source paper describes a monotonically-increasing clearance-vs-concentration relationship arising from a saturable protein-binding buffer (e.g., FVIII / von Willebrand factor complex saturation at supraphysiological rFVIII doses). Larger `beta_cl` gives a steeper rise in CL with C; `beta_cl = 0` recovers linear elimination. The bare counterpart inside `model` is `beta_cl` (units of 1 / concentration).
- **Source aliases:**
  - `β` -- Greek letter used in Larsen 2018 Table 3 and Eq. 1.
- **Example models:** `Larsen_2018_factorviii_rat.R`, `Larsen_2018_factorviii_monkey.R` (founding examples; rat `beta_cl = 0.162 mL/IU`, monkey `beta_cl = 0.0355 mL/IU`, both fitting the exponential-nonlinear-CL form on total FVIII activity).
- **Notes:** Distinct from `lhill` (sigmoidal Emax / Imax exponent), `lgamma` (Friberg / TGI kinetic exponents), and `lvmax` / `km` (Michaelis-Menten saturable elimination — a different mathematical form). Choose `lbeta_cl` when the source paper explicitly parameterises CL as an exponential function of concentration; choose `lvmax` / `km` when the paper fits Michaelis-Menten saturation instead. The `_cl` suffix marks that the slope acts on clearance; parallel `beta_<param>` names are permitted for other parameters when a paper extends the exponential-nonlinear form to `V` or `Q`.

### lcmpr (**canonical log-transformed milk-to-plasma concentration ratio**)
- **Type:** log-transformed-pk
- **Role:** Log of the ratio of drug concentration in breast milk to the concentration in maternal plasma, used to derive a milk observable algebraically from the plasma compartment (`Cmilk <- cmpr * Cc`) when the source data are too sparse to support a separate milk compartment with its own transfer clearances. Unitless. May be estimated with IIV (`etalcmpr`) and may carry covariate effects, e.g. a power function of time postpartum during the colostrum period. The bare counterpart inside `model` is `cmpr`.
- **Source aliases:**
  - `MPRcon` -- Li 2023 ornidazole notation for the milk-to-plasma concentration ratio.
  - `MP` -- common lactation-literature shorthand for the milk:plasma ratio.
- **Example models:** `Li_2023_ornidazole.R` (founding example; `lcmpr = log(0.58)`, RSE 8.63%, IIV variance 0.327, power effect of time postpartum with exponent 1.37 centred on a back-solved median postpartum sampling time of 54 h).
- **Notes:** Distinct from the `lkp_<tissue>` partition-coefficient family (including `lkp_milk`) in that `cmpr` is an **estimated popPK parameter** fitted to paired maternal plasma and milk observations, whereas a `kp` is a physicochemically predicted or literature-fixed partition constant used inside a mechanistic distribution model. Use `cmpr` when the paper estimates a concentration ratio as a model parameter; use `kp_milk` when the paper supplies a milk:plasma partition coefficient as a mechanistic input to a milk compartment. Use an AUC-based ratio name only if the paper actually parameterises the model on an exposure ratio -- Li 2023 reports `MPRauc` as a derived simulation **output**, not a model parameter, and it is deliberately not encoded. Pair with the `Cmilk` observable and the `propSd_Cmilk` residual.

### lkp_adipose, lkp_brain, lkp_intestine, lkp_kidney, lkp_liver, lkp_lung, lkp_milk, lkp_muscle, lkp_other, lkp_rest, lkp_trachea (**canonical log-transformed tissue-to-plasma partition coefficients**)
- **Type:** log-transformed-pk
- **Role:** Log of the tissue-to-plasma (or tissue-to-blood) partition coefficient Kp for one named anatomical tissue, i.e. the equilibrium ratio of tissue concentration to plasma (or blood) concentration. Unitless. The family shape is `lkp_<tissue>`, where `<tissue>` is the same lowercase anatomical token used for the corresponding PBPK compartment name, so `lkp_liver` is the partition coefficient of the `liver` compartment. Each drives the perfusion- or permeability-limited distribution term for its tissue (e.g. `d/dt(liver) <- q_liver * (Cc - liver / vliver / kp_liver) -...`). The bare counterparts inside `model()` are `kp_<tissue>`, registered below. Usually fixed from the source paper's physicochemical predictions (Rodgers-Rowland, Poulin-Theil) or measured tissue:plasma ratios, but may be estimated when the paper fits them.
- **Source aliases:**
  - `Kp` -- near-universal PBPK notation, subscripted by tissue in the source tables.
  - `Kt` / `Ptp` -- tissue:plasma partition notation used by some perfusion-limited PBPK papers.
- **Example models:** `Levitt_2005_propofol_pbpk.R` (`lkp_adipose`, `lkp_brain`, `lkp_liver`, `lkp_intestine`, `lkp_rest`), `Mi_2023_cefquinome_pbpk.R` (`lkp_liver`, `lkp_kidney`, `lkp_lung`, `lkp_muscle`, `lkp_other`), `Kang_2023_artesunate_hamster_pbpk.R` and `Kang_2023_pyronaridine_hamster_pbpk.R` (`lkp_lung`, `lkp_trachea`, `lkp_other`, plus metabolite-suffixed `lkp_lung_dihydroart` etc.).
- **Notes:** This entry formally registers a family that had been de-facto in use across the PBPK model files listed above since before the register existed; the members enumerated in the heading are those currently used, and a new anatomical tissue may be added to the heading without a fresh naming sidecar so long as `<tissue>` matches the canonical compartment name in `compartment-names.md`. `lkp_milk` is the lactation member: the milk:plasma partition coefficient for mechanistic milk-compartment models. Distinct from `lkpu<n>`, which is a *cluster*-indexed unbound partition coefficient shared by several tissues rather than one named organ, and from `sf<n>`, which scales a predicted Kpu rather than replacing it. Distinct from `lcmpr` above, which is an estimated popPK milk:plasma **concentration ratio** rather than a mechanistic partition constant. Two existing names collide on prefix but are NOT members of this family and must not be read as partition coefficients: `lkp_f` and `lkp_hb` in `Li_2015_taspoglutide_mbma.R` are placebo-response rate constants for fasting plasma glucose and HbA1c, and `lkpin` / `lkpout` are precursor-pool turnover rate constants. Family registered 2026-08-05 alongside the Li 2023 ornidazole extraction.

### lfuA, lkfu, lfuB (**canonical log-transformed mono-exponentially decaying free fraction**)
- **Type:** log-transformed-pk
- **Role:** Three-parameter description of a plasma free (unbound) fraction that decays exponentially after each dose towards a non-zero floor, the shape produced by irreversible, time-dependent plasma protein binding. `lfuA` is the log amplitude of the time-decaying component (unitless), `lkfu` the log first-order rate constant of the binding decay (1 / time), and `lfuB` the log time-independent, asymptotic free fraction (unitless). The bare counterparts inside `model` are `fuA`, `kfu` and `fuB`, combined as `fu <- fuA * exp(-kfu * tad(<cmt>)) + fuB`, so the free fraction starts each dosing interval at `fuA + fuB` and decays towards `fuB`. The decay clock is time **after the last dose**, so it must be driven by `tad(<cmt>)` rather than by `t`.
- **Source aliases:**
  - `fp0_A`, `alpha`, `fp0_B` -- used in `Yamada_2025_oxaliplatin.R` (Yamada 2025 Article Equation b, `fp = fp0_A * exp(-alpha * TALD) + fp0_B`, where `fp` is the free fraction of platinum in plasma and `TALD` is time after the start of infusion of the last dose).
- **Example models:** `Yamada_2025_oxaliplatin.R` (simultaneous free + total plasma platinum after oxaliplatin; `fp = 0.194 * exp(-0.0393 * TALD) + 0.0318`, all three estimated, with IIV on `fp0_A` and `fp0_B` but none on `alpha` -- founding example).
- **Notes:** Distinct from `lbfu` (registered in the skill-side naming reference), which is the **linear** time-varying unbound-fraction slope of the form `fu = fu_ref + bfu * (t - t_ref)` anchored to a study-level reference time and usually fixed from literature. Choose between the two by the source paper's own printed equation, not by convenience: `lbfu` for a linear drift with a literature-fixed slope, this family for an estimated exponential decay re-started by each dose. Unlike `lbfu`, these three are normally estimated, because a paper that fits free and total concentrations simultaneously identifies them directly from the free:total ratio. Characteristic of platinum drugs, whose binding to plasma protein is irreversible and time dependent rather than a rapid reversible equilibrium; do not reuse this family for concentration-dependent (saturable) binding, which needs a capacity/affinity parameterisation instead.

### lkst (**canonical log-transformed gastric emptying rate constant**)
- **Type:** log-transformed-pk
- **Role:** First-order rate constant for gastric emptying: transit of drug out of a `stomach` compartment into the first intestinal segment / intestinal lumen, in oral PBPK models that resolve the gastric and intestinal lumen as explicit states (1 / time). The bare counterpart inside `model` is `kst`. Founding member of the lumped gastrointestinal transit family below.
- **Source aliases:**
  - `Kst` -- used in `Ai_2024_ractopamine_goat_pbpk.R`; also Sun 2026 Table 2 ("gastric emptying rate constant", 0.8 /h) and the Supporting Information Berkeley Madonna listing (`RAST = RDOSEoral - Kst * AST`).
  - `kst`, `k_st`, `KST` -- generic acslX / Berkeley-Madonna spellings.
- **Example models:** `Ai_2024_ractopamine_goat_pbpk.R` (`Kst` = 0.0910 1/h, gastric contents into gut lumen in goats; founding example), `Yang_2025_matrine_pig_pbpk.R` (`kst` = 0.8545 1/h, reconstructed gastric depot into the sampled intestinal lumen in pigs), `Sun_2026_tilmicosin_pbpk.R` (`Kst` = 0.8 1/h in swine).
- **Notes:** Distinct from `lka`, which moves drug from a depot into the *systemic* circulation; `lkst` moves drug between two luminal compartments and is not an absorption term. Precedent existed in the repository before ratification -- `Ai_2024_ractopamine_goat_pbpk.R` already carried `lkst` in `ini()`, and `Luo_2024_*_pbpk.R`, `Back_2018_fenofibrate.R` and `Guiastrennec_2016_gastric_emptying.R` carried the same concept as `model()` literals named `kt_sto` / `kg`. Register the constant in `ini()` rather than burying it as a `model()` literal whenever the source paper optimised it, so that half of one optimisation is not hidden from the reader.

### lkt_duodenum (**canonical log-transformed duodenal transit rate constant**)
- **Type:** log-transformed-pk
- **Role:** First-order rate constant for luminal transit out of the `duodenum` into the next small-intestinal segment (1 / time).
- **Source aliases:**
  - `Kd` -- Sun 2026 Table 2 ("duodenal transit rate constant", 2.2 /h) and the Supporting Information listing (`RDI = Kst * AST - Kd * ADI - Ka1 * ADI + Rmet`).
- **Example models:** `Sun_2026_tilmicosin_pbpk.R` (founding example).
- **Notes:** Paired with `lka_duodenum`, which drains the same compartment by absorption rather than by transit; the two compete for the duodenal amount and must not be conflated. Named `lkt_<segment>` (transit out of the named segment) so the family extends to any resolved luminal segment.

### lkt_ileocolic (**canonical log-transformed ileocolic transit rate constant**)
- **Type:** log-transformed-pk
- **Role:** First-order rate constant for luminal transit across the ileocolic junction, i.e. out of the small intestine into the large intestine (1 / time).
- **Source aliases:**
  - `Kint` -- Sun 2026 Table 2 ("ileocolic transit rate constant", 2.0 /h) and the Supporting Information listing (`ROSI = Kdm * ADI - Ka2 * AOSI - Kint * AOSI`).
- **Example models:** `Sun_2026_tilmicosin_pbpk.R` (founding example).
- **Notes:** Named for the anatomical junction it crosses rather than for the donating segment, following the source paper's own term, so the name stays correct whether the donor is a lumped `a_small_intestine` or a resolved `ileum`.

### lka_duodenum (**canonical log-transformed duodenal absorption rate constant**)
- **Type:** log-transformed-pk
- **Role:** First-order rate constant for absorption of drug out of the `duodenum` lumen into the systemic or portal circulation (1 / time). Segment-resolved counterpart of `lka`.
- **Source aliases:**
  - `Ka1` -- Sun 2026 Table 2 ("duodenal absorption rate constant", 0.8 /h) and the Supporting Information listing (`RAO = Ka1 * ADI + Ka2 * AOSI`).
- **Example models:** `Sun_2026_tilmicosin_pbpk.R` (founding example).
- **Notes:** Use the segment-suffixed `lka_<segment>` form only when the model resolves more than one absorbing luminal segment with distinct rate constants; a single-site oral model keeps the plain `lka`.

### lka_small_intestine (**canonical log-transformed lumped small-intestinal absorption rate constant**)
- **Type:** log-transformed-pk
- **Role:** First-order rate constant for absorption out of the lumped post-duodenal small intestine (`a_small_intestine`) into the systemic or portal circulation (1 / time). Sibling of `lka_duodenum`.
- **Source aliases:**
  - `Ka2` -- Sun 2026 Table 2 ("remaining small intestinal absorption rate constant", 0.007 /h) and the Supporting Information listing (`RAO = Ka1 * ADI + Ka2 * AOSI`).
- **Example models:** `Sun_2026_tilmicosin_pbpk.R` (founding example).
- **Notes:** Segment suffix matches the compartment it drains, so it is `_small_intestine` (lumped) rather than `_jejunum` / `_ileum` (resolved).

### lkf (**canonical log-transformed faecal excretion rate constant**)
- **Type:** log-transformed-pk
- **Role:** First-order rate constant for excretion of drug from the terminal luminal compartment into faeces (1 / time). Drains `a_large_intestine` into the `a_feces` sink.
- **Source aliases:**
  - `Kf` -- Sun 2026 Table 2 ("fecal excretion rate", 0.05 /h) and the Supporting Information listing (`Rf = Kf * ALI`).
- **Example models:** `Sun_2026_tilmicosin_pbpk.R` (founding example).
- **Notes:** A transit constant out of the body, not a clearance: it multiplies a luminal *amount*, so it carries units of 1 / time rather than volume / time. Distinct from `lcl_nonren`, which is the hepatobiliary clearance feeding drug *into* the gut.

---

### lkfec (**canonical log-transformed faecal excretion rate constant**)
- **Type:** log-transformed-pk
- **Role:** First-order rate constant for unabsorbed drug leaving the intestinal lumen as faeces, competing with absorption out of the same compartment (1 / time). The bare counterpart inside `model` is `kfec`.
- **Source aliases:**
  - `ke` -- used in `Yang_2025_matrine_pig_pbpk.R`.
  - `Kgut`, `kgut`, `kint`, `kF`, `k_e` -- generic spellings; note that `Kgut` / `kint` name the compartment or the interval, not the process.
- **Example models:** `Yang_2025_matrine_pig_pbpk.R` (`ke` = 0.007358 1/h, unabsorbed matrine leaving the pig intestinal lumen; founding example).
- **Notes:** `lkfec` was chosen over the alternative `lkgut` because it names the **process** -- excretion out of the gut lumen to faeces -- rather than the compartment; `lkgut` would read to a later user as a general gut-transit rate and could be misapplied to the absorption or the emptying step. `lkgut` remains in `Yang_2023_diclazuril_chicken_pbpk.R` and `Ai_2024_ractopamine_goat_pbpk.R` as an unregistered legacy spelling that predates this entry; use `lkfec` in new models. Do **not** mint a gut-specific absorption name for the competing route out of the same compartment -- that remains the ordinary canonical `lka`.

### lkbile (**canonical log-transformed biliary excretion rate constant**)
- **Type:** log-transformed-pk
- **Role:** First-order rate constant for **parent drug** leaving the `liver` compartment by biliary excretion, either into `gut_lumen` (enterohepatic recirculation) or into a terminal faecal / bile sink (1 / time). The bare counterpart inside `model` is `kbile`.
- **Source aliases:**
  - `KbileC` -- used in `Zhang_2024_f53b_mouse_pbpk.R`.
  - `kbi` -- used in `Yang_2025_matrine_pig_pbpk.R`.
  - `Kbile`, `k_bi` -- generic spellings.
- **Example models:** `Zhang_2024_f53b_mouse_pbpk.R` (`KbileC` = 0.00001, allometrically scaled terminal biliary elimination from liver to faeces; founding example), `Yang_2025_matrine_pig_pbpk.R` (`kbi` = 0.05835 1/h, liver back into the intestinal lumen, which is what produces the observed two-phase luminal decay).
- **Notes:** Deliberately **not** `kbm`, and `kbm` must not be widened to cover this case: `kbm` is registered as the biliary-*metabolite* excretion rate constant, moving a **metabolite** out of a **plasma / central** compartment, whereas `lkbile` moves **parent drug** out of **`liver`**. Broadening `kbm` would silently change the meaning of an existing entry that `Hamren_2008_tesaglitazar.R` depends on. The two can coexist in one model: a parent drug excreted in bile via `lkbile` and its glucuronide returned via `kbm` / `kicv`.

### lwdist (**canonical log-transformed amplitude of a transient distribution-phase elimination flux**)
- **Type:** log-transformed-pk
- **Role:** Log-scale amplitude of an additive, dose-proportional elimination flux that is active immediately after a dose and decays exponentially with time after dose, standing in for a distribution phase in a model that carries no peripheral state: `d/dt(central) <- -kel * central - D * wdist * exp(-kdist * tad())`. Units 1 / time, since the flux is `dose * wdist`.
- **Source aliases:**
  - `w` -- Braem 2026 Equation 21 and Table S2.
  - `q` -- Braem 2026 Equation 20 (the same quantity; the paper switches symbol between the ODE and its explicit solution).
- **Example models:** `Bram_2026_biexponential.R` (`wdist = 0.1716832 1/h`, back-solved from the Table S2 macro-parameters).
- **Notes:** Ratified 2026-09-02 with the Braem 2026 extraction (sidecar `oare_PMC13291805` request-001 q2 = A). Deliberately **not** a member of the registered `cl_exp_` / `cl_time_` time-varying-clearance families: those name a decaying component **of a clearance** that multiplies the state, whereas this term is an additive, state-independent flux proportional to the dose. Recording it as a clearance would be a category error visible to anyone reading the register. The `dist` token names what the term means -- Braem 2026's own Discussion reads it as "an increased observed elimination from the central compartment which corresponds to a distribution process until equilibration is reached" -- and keeps it clear of the PD `hill` / `emax` namespace. Pairs 1:1 with `lkdist`. Because the term carries no memory of previous doses, a model using it is valid for multiple dosing only when each dose falls after the previous distribution phase completes; the pseudo-compartment repair is in Braem 2026 Figure S2.

### lkdist (**canonical log-transformed decay rate of a transient distribution-phase elimination flux**)
- **Type:** log-transformed-pk
- **Role:** Log-scale first-order rate constant at which the `wdist` transient distribution-phase flux decays with time after dose (1 / time). In the explicit solution it is the fast macro-exponent, i.e. the `alpha` of the equivalent two-compartment model.
- **Source aliases:**
  - `p` -- Braem 2026 Equations 20 and 21 and Table S2.
- **Example models:** `Bram_2026_biexponential.R` (`kdist = 0.46 1/h`, Table S2).
- **Notes:** Ratified 2026-09-02 with `lwdist` (same sidecar). Not `lkel_exp_kdes`, which is the relaxation rate of a time-varying *elimination rate constant*; `kdist` decays a flux, not a rate constant. Not `lkde` / `lkel` either -- the ordinary first-order elimination of the same model is a separate parameter (`lkel`) and both appear together.

---

### lflnode (**canonical log-transformed lymphatic-absorption fraction**)
- **Type:** log-transformed-pk
- **Role:** Fraction of a subcutaneously absorbed dose that is routed through the `lnode` lymph-node compartment rather than entering `central` directly through the blood capillaries of the SC tissue (unitless, bounded in (0, 1)). The complementary fraction `1 - flnode` reaches plasma directly, so a single absorption rate constant `ka` empties the `depot` and the two destinations split its flux. The bare counterpart inside `model()` is `flnode`.
- **Source aliases:**
  - `Frc` -- Wu 2012 Table II and Eqs. (3), (5); Fig. 1 legend ("The fraction absorbed by the lymphatic compartment is given by Frc, and the fraction taken up by the blood capillaries is 1-Frc").
- **Example models:** `Wu_2012_bevacizumab_mouse.R` (`Frc` = 0.00964, CV 19.6%; founding example -- bevacizumab-IRDye 800CW after a 0.45 mg/kg footpad SC dose in SKH-1 mice).
- **Notes:** Ratified 2026-09-02 (operator sidecar `oare_PMC3326166` request-001 / response-001, question q2, option A). Named by DESTINATION COMPARTMENT, exactly as `fdepot` names the fraction routed to `depot`, so the name inherits the `lnode` token and stays synchronised with `lv_lnode` and `lk_lnode_central`. Deliberately **not** any `flymph` / `f_lymph` spelling: `f_lymph` is already in use in `Lindauer_2017_pembrolizumab.R` and `Michelet_2025_BI754111_mpbpk.R` for a **different** quantity -- lymph flow as a fraction of plasma flow -- and the two would be indistinguishable to a reader. Also deliberately not `lffo`, which is the fraction routed to the FAST site of a parallel two-site absorption model; the lymphatic route is the slow one (into the node at `ka`, out of it at the much smaller `k_lnode_central`), so `lffo` would assert the opposite of the mechanism. Held on the log scale to match the sibling fraction canonicals `lfdepot` / `lffo`; use the `logitflnode` spelling instead if a future source estimates the fraction on the logit scale or carries IIV on it, per the reasoning in the `logitfdepot` entry.

### lv_lnode (**canonical log-transformed lymph-node volume**)
- **Type:** log-transformed-pk
- **Role:** Volume of distribution of the `lnode` lymph-node compartment, used to convert the node amount to the measured node concentration (`Clnode <- lnode / v_lnode`). The bare counterpart inside `model()` is `v_lnode`.
- **Source aliases:**
  - `VLN` -- Wu 2012 Table II, Eq. (7) (`ALN,sc = CLN,sc x VLN`) and Fig. 1 legend.
- **Example models:** `Wu_2012_bevacizumab_mouse.R` (`VLN` = 0.33 mL/kg, FIXED from the measured axillary-node weight at an assumed 1 g/mL specific density over mean SKH-1 body weight; founding example). Bare-form precedent that predates this entry: `Ramachandran_2023_rifampicin_pbpk.R`, `Ramachandran_2023_isoniazid_pbpk.R`, `Ramachandran_2023_ethambutol_pbpk.R`, `Ramachandran_2023_pyrazinamide_pbpk.R` (`v_lnode` = 0.274 L, a hardcoded physiological constant).
- **Notes:** Ratified 2026-09-02 (operator sidecar `oare_PMC3326166` request-001 / response-001, question q3, option A). Member of the `lv_<space>` family founded by `lv_elf`, and registered for the same reason that entry gives: the bare volume canonicals `lvc` / `lvp` / `lvp2` / `lvp3` are reserved for the central and numbered peripheral compartments of a classical-PK model, and a lymph node is a named physiological space rather than a peripheral compartment. Preferred over the paper's own compact `lvln` spelling because `lv_lnode` carries the `lnode` compartment token, which is already established in this register by `kp_lnode` / `lkp_lnode`. Register it in `ini()` -- wrapped in `fixed()` when the source fixes it from anatomy, as Wu 2012 does -- rather than hardcoding it in `model()`, so the value stays visible to anyone refitting the model.

### lk_lnode_central (**canonical log-transformed lymph-node-to-central transfer rate constant**)
- **Type:** log-transformed-pk
- **Role:** First-order rate constant carrying drug out of the `lnode` lymph-node compartment into `central` (1 / time) -- the efferent-lymphatic return limb of a subcutaneous-absorption model. The bare counterpart inside `model()` is `k_lnode_central`.
- **Source aliases:**
  - `ka2` -- Wu 2012 Table II and Eqs. (3), (5); Fig. 1 legend ("ka2 is the first-order rate constant describing transfer from the lymphatic system to plasma").
- **Example models:** `Wu_2012_bevacizumab_mouse.R` (`ka2` = 0.723 /h, CV 9.64%; founding example).
- **Notes:** Ratified 2026-09-02 (operator sidecar `oare_PMC3326166` request-001 / response-001, question q4, option A). Member of the registered `k_<from>_<to>` directional-transfer family (`k_central_elf` / `k_elf_central`, `k_central_milk` / `k_milk_central`, `k_csf_plasma`, `k_presystemic_central`), whose own note directs that the family be used whenever a transfer connects `central` to a named non-numbered compartment. Deliberately **not** `lka2`, even though that is the paper's literal symbol: this register states that `lka1` / `lka2` are not canonical because a bare ordinal encodes only sequence, and `lka2` is already in wide use (`Mauro_2025_nilotinib.R`, `Khwarg_2024_donepezil_im.R`, `Perlstein_2025_risperidone_tv46000.R`, `Qi_2024_vosoritide.R`) for absorption from a SECOND DEPOT. Also not `lka_lnode`: the site-labelled absorption family (`lka_duodenum`, `lka_small_intestine`) describes absorption from a lumen outside the body, whereas the node is a modelled internal state, so this is a transfer and not an absorption.

---

## Bare structural PK parameters

The bare counterparts of the log-transformed parameters above. Used when the source paper estimates the parameter directly on the linear scale, or when the parameter appears in the `model` block as the exponentiated form `<base> <- exp(l<base> + eta_<base>)`.

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
- **Role:** Bare (natural-scale) form of `lkel`: the first-order elimination rate constant, both in K-PD / single-rate-constant elimination models with no explicit `vc` and in rate-constant-parameterised disposition models that estimate an explicit `vc` beside it (see `lkel` for why such a model must not be reparameterised to `cl` / `vc`).
- **Source aliases:**
  - `ke` -- legacy.
  - `kde` -- paper-named (Mazzocco 2015 / Shoji 2017 / Xia 2024 KDE) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
  - `kp` -- paper-named (van Hasselt 2015 KP) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
  - `ps_elim`, `pc_elim` -- paper-named (Wilson 2015 p_S / p_C) bare drug-specific K-PD elim rates; replaced 2026-05-30 by `kel_sunitinib` / `kel_irinotecan`.
- **Example models:** `Mazzocco_2015_temozolomide.R`, `Shoji_2017_fosdagrocorat_oc.R`, `Shoji_2017_fosdagrocorat_p1np.R`, `vanHasselt_2015_eribulin.R`, `Wilson_2015_sunitinib_irinotecan_mouse.R` (bare drug-suffixed `kel_<drug>`), `Xia_2024_warfarin.R`, `Baklouti_2026_amoxicillin.R` (rate-constant disposition with an explicit `vc`).

### kel_exp_kdes (**canonical bare decay-rate constant of a time-varying elimination rate constant**)
- **Type:** bare-pk
- **Role:** Bare (natural-scale) form of `lkel_exp_kdes`: the rate constant at which a time-varying elimination rate constant relaxes toward its asymptote (1 / time).
- **Source aliases:**
  - `kappa_k` -- used in `Schreib_2024_busulfan.R`.
- **Example models:** `Schreib_2024_busulfan.R` (`kel_exp_kdes = exp(lkel_exp_kdes) = 0.0516 1/h`).
- **Notes:** Bare counterpart of `lkel_exp_kdes`; see that entry for the functional form. The sibling `kel_exp_famp` is registered under "Paper-named mechanistic parameters" instead of here, because it is signed and must never be log-transformed.

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
- **Role:** First-order transfer rate constants between `central` and the canonical `elf` compartment (1 / time). Deliberately a directional PAIR rather than a single inter-compartmental clearance, because plasma-to-ELF distribution is typically asymmetric -- `k_central_elf < k_elf_central` is what generates the sub-unity ELF penetration ratio that lung PK/PD analyses turn on, and a symmetric `q` / `velf` parameterisation cannot express it. The log-transformed `lk_central_elf` / `lk_elf_central` forms are used in `ini`.
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


### cl_time_max (**canonical maximum magnitude of a sigmoidal-in-time clearance change**)
- **Type:** bare-pk
- **Role:** Magnitude of the log-clearance change approached as `t >> cl_t50`, in the
  sigmoidal-in-time clearance form
  `td_cl <- exp(cl_time_max_i * t^cl_time_hill / (cl_t50^cl_time_hill + t^cl_time_hill))`.
  Log form `lcl_time_max`; the per-individual value after IIV is written `cl_time_max_i`.
- **Source aliases:**
  - `cl_hill_max` / `lcl_hill_max` -- prior name (pre-2026-08-31 rename; see Notes).
- **Example models:** `Kuchimanchi_2024_dostarlimab.R`, `Masters_2022_avelumab.R`, `PK_2cmt_tdcl_des.R`.
  STRUCTURE (time-varying clearance, `cl_time_*`) rather than only the shape of the curve.
  `hill` now appears solely on the shape coefficient `cl_time_hill`, where it is meaningful.
  Ties the family to the registered `cl_time` / `lcl_time` stem. The trio was used in ~20
  shipped models for a long time without a register entry; documented here as part of that
  rename. Distinct from maturation Hill coefficients such as `e_page_cl_hill`
  (`Tikiso_2021_abacavir.R`), which are covariate effects on CL maturation, not time.

### cl_t50 (**canonical half-time of a sigmoidal-in-time clearance change**)
- **Type:** bare-pk
- **Role:** Time at which half of `cl_time_max` has been reached, in the same sigmoidal
form. Log form `lcl_t50`. Units are the model's time unit.
- **Source aliases:**
  - `cl_hill_t50` / `lcl_hill_t50` -- prior name (pre-2026-08-31 rename).
- **Example models:** `Kuchimanchi_2024_dostarlimab.R`, `Masters_2022_avelumab.R`.
- **Notes:** Named without a `time` infix because `t50` already denotes a time; see
[[cl_time_max]] for the rename rationale.

### cl_time_hill (**canonical Hill coefficient of a sigmoidal-in-time clearance change**)
- **Type:** bare-pk
- **Role:** Sigmoid steepness exponent of the same time-varying clearance form. Log form
`lcl_time_hill`.
- **Source aliases:**
  - `cl_hill_gamma` / `lcl_hill_gamma` -- prior name (pre-2026-08-31 rename).
- **Example models:** `Kuchimanchi_2024_dostarlimab.R`, `Masters_2022_avelumab.R`.
- **Notes:** This is the one member where `hill` is the correct role token -- it names the
shape coefficient itself. See [[cl_time_max]] for the rename rationale.

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

### cl_tsnet (**canonical bare net-tubular-secretion clearance arm**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lcl_tsnet`. Net renal tubular-secretion component (secretion minus reabsorption) of a renal clearance decomposition `CL_renal_total = GFR + CL_tsnet`.
- **Source aliases:** `nCTS`, `CL_SEC`.
- **Example models:** `Chen_2025_iohexol_creatinine.R` (as the analyte-suffixed `cl_tsnet_creatinine` derived in `model`).

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
- **Role:** Bare counterpart of `lgam1`. Unitless Weibull shape (sigmoidicity) parameter inside a Piotrovskij-style time-dependent ka, or of the first process of a Weibull release input function.
- **Source aliases:**
  - `GAM1` / `GAMMA1` -- NONMEM convention.
  - `SS` -- release-function "sigmoidicity factor"; cannot be carried across literally because `ss` is a reserved rxode2 symbol.
- **Example models:** `Desai_2016_isavuconazole.R` (founding example), `Perlstein_2026_olanzapine_lai.R` (release-input form).

### ra2 (**canonical bare second-process Weibull release rate-scaling parameter**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lra2`. Rate-scaling parameter of the second Weibull process in a multi-phase release input function (1 / time). Used inside `model` after being exponentiated from `lra2`, and after any covariate effect on the release TIME has been applied on the reciprocal.
- **Source aliases:**
  - `TD1` / `TD1_0` -- the reciprocal release time; `RA2` -- NONMEM convention.
- **Example models:** `Perlstein_2026_olanzapine_lai.R` (founding example).

### gam2 (**canonical bare second-process Weibull release shape parameter**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lgam2`. Unitless Weibull shape (sigmoidicity) parameter of the second process of a Weibull release input function.
- **Source aliases:**
  - `SS1` -- release-function "sigmoidicity factor for the second process"; `GAM2` / `GAMMA2` -- NONMEM convention.
- **Example models:** `Perlstein_2026_olanzapine_lai.R` (founding example).

### frel (**canonical bare release-process fraction**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `logitfrel`. Fraction of the dose entering the first process of a multi-phase release input function (unitless, in (0, 1)); the second process receives `1 - frel`. Used inside `model` after the inverse-logit back-transform, and normally applied as the bioavailability split across two parallel depot compartments (`f(depot) <- frel`, `f(depot2) <- 1 - frel`).
- **Source aliases:**
  - `FF` -- "fraction of the available dose released in the first process".
- **Example models:** `Perlstein_2026_olanzapine_lai.R` (founding example).

### beta_cl (**canonical bare exponential-nonlinear-CL slope**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lbeta_cl`. Concentration-slope coefficient in an exponential-nonlinear-clearance function of the form `CL(C) = CL_0 * exp(beta_cl * C)`, units of 1 / concentration. Used inside `model` after being exponentiated from `lbeta_cl`.
- **Source aliases:**
  - `β` -- Greek letter used in Larsen 2018 Table 3 and Eq. 1.
- **Example models:** `Larsen_2018_factorviii_rat.R`, `Larsen_2018_factorviii_monkey.R` (founding examples).

### cmpr (**canonical bare milk-to-plasma concentration ratio**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lcmpr`. Unitless ratio of breast-milk to maternal-plasma drug concentration, used inside `model` after being exponentiated from `lcmpr` (`cmpr <- exp(lcmpr + etalcmpr) * <covariate terms>`) and then applied to derive the milk observable, `Cmilk <- cmpr * Cc`.
- **Source aliases:**
  - `MPRcon` -- Li 2023 ornidazole notation.
- **Example models:** `Li_2023_ornidazole.R` (founding example).
- **Notes:** See `lcmpr` above for the distinction from the `kp_<tissue>` partition-coefficient family and for why an AUC-based ratio is not encoded as a model parameter.

### kp_adipose, kp_bone, kp_brain, kp_cerebellum, kp_choroid_plexus, kp_cortex, kp_csf, kp_gut, kp_heart, kp_hippocampus, kp_intestine, kp_kidney, kp_liver, kp_lnode, kp_lung, kp_milk, kp_muscle, kp_other, kp_rest, kp_skin, kp_spleen, kp_striatum, kp_testes, kp_trachea (**canonical bare tissue-to-plasma partition coefficients**)
- **Type:** bare-pk
- **Role:** Bare counterparts of the `lkp_<tissue>` family. Unitless equilibrium ratio of concentration in one named anatomical tissue to concentration in plasma (or blood), used inside `model` either as the exponentiated form of an `ini` entry (`kp_liver <- exp(lkp_liver)`) or as a value computed in place from the source paper's physicochemical prediction equations.
- **Source aliases:**
  - `Kp` -- near-universal PBPK notation, subscripted by tissue.
- **Example models:** `Levitt_2005_propofol_pbpk.R`, `Mi_2023_cefquinome_pbpk.R`, `Gaohua_2012_pregnancy_pbpk_midazolam.R`, `Litjens_2023_linezolid_cns_pbpk.R`, `Yang_2023_diclazuril_chicken_pbpk.R`, `Kang_2023_artesunate_hamster_pbpk.R`.
- **Notes:** Many PBPK models compute additional bare `kp_*` names that are qualifiers rather than tissues (`kp_free`, `kp_bound`, `kp_preg`) or per-tissue loop indices (`kp_i`); those are local `model()` intermediates and are intentionally not registered as canonicals.

### kst, kfec, kbile (**canonical bare gastrointestinal transit and enterohepatic rate constants**)
- **Type:** bare-pk
- **Role:** Bare counterparts of `lkst`, `lkfec` and `lkbile`. `kst` empties a `stomach` depot into the intestinal lumen; `kfec` removes unabsorbed drug from the intestinal lumen as faeces, in competition with `ka`; `kbile` moves parent drug out of `liver` into `gut_lumen` or a terminal bile sink. All are 1 / time.
- **Source aliases:**
  - `Kst`, `kt_sto`, `kg` -- gastric emptying.
  - `ke`, `Kgut`, `kint`, `kF` -- faecal excretion.
  - `KbileC`, `kbi`, `Kbile` -- biliary excretion.
- **Example models:** `Ai_2024_ractopamine_goat_pbpk.R`, `Zhang_2024_f53b_mouse_pbpk.R`, `Yang_2025_matrine_pig_pbpk.R`.
- **Notes:** Some PBPK models compute these in place as `model()` literals rather than exponentiating an `ini()` entry; register them in `ini()` whenever the source paper reports them as optimised quantities.

---

### tclchange (**canonical bare clearance-step breakpoint time**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `ltclchange`. Time from treatment initiation at which a piecewise-constant clearance steps between its early and late arms.
- **Source aliases:**
  - `TCLchange` -- Park 2025.
  - `tNab` -- Yoneyama 2017.
- **Example models:** `Park_2025_efineptakin_alfa.R`.
- **Notes:** The early arm is carried by the plain `cl`; there is no `cl_early`. See `ltclchange` for the full form.

### cl_late (**canonical bare post-breakpoint clearance arm**)
- **Type:** bare-pk
- **Role:** Bare counterpart of `lcl_late`. Clearance in force after the `tclchange` breakpoint.
- **Source aliases:**
  - `CL>TCLchange` -- Park 2025.
- **Example models:** `Park_2025_efineptakin_alfa.R`.
- **Notes:** Paired with `tclchange`; the pre-breakpoint arm is the plain `cl`.

### flnode, v_lnode, k_lnode_central (**canonical bare lymph-node absorption-limb parameters**)
- **Type:** bare-pk
- **Role:** Bare counterparts of `lflnode`, `lv_lnode` and `lk_lnode_central`. `flnode` is the unitless fraction of an absorbed subcutaneous dose routed through the `lnode` compartment rather than straight into `central`; `v_lnode` is the lymph-node volume that converts the node amount to the `Clnode` observable; `k_lnode_central` is the first-order rate constant returning drug from the node to plasma (1 / time). Together they are the lymphatic limb of a subcutaneous-absorption model.
- **Source aliases:**
  - `Frc` -- lymphatic absorption fraction (Wu 2012).
  - `VLN` -- lymph-node volume (Wu 2012).
  - `ka2` -- node-to-plasma return rate constant (Wu 2012).
- **Example models:** `Wu_2012_bevacizumab_mouse.R`. Bare-form `v_lnode` precedent as a hardcoded physiological constant: `Ramachandran_2023_rifampicin_pbpk.R` and its isoniazid / ethambutol / pyrazinamide siblings.
- **Notes:** Registered 2026-09-02 together with the log-transformed family; see `lflnode`, `lv_lnode` and `lk_lnode_central` above for the ratification rationale, for the collision against the existing `f_lymph` (lymph FLOW fraction), and for why `ka2` is not carried as `lka2`. Pair with the `lnode` compartment, the `Clnode` observable, and the `propSd_Clnode` residual.

---

## Paper-named mechanistic parameters

Parameters that don't fit the standard `ka` / `cl` / `vc` shape but recur across published models. Each entry is treated as a canonical bare name; the log-transformed form (`l<name>`) is acceptable wherever the parameter is strictly positive and the source paper reports an exponential typical-value form. Add to this list rather than introducing a new ad-hoc pattern.

### lrbase (**canonical log-transformed baseline level of a PD / turnover state**)
- **Type:** paper-named-param
- **Role:** Log-transformed baseline (pre-treatment steady-state) value of a PD, turnover or endogenous state, expressed in the units of the state itself -- not a concentration and not an increment. Inside `model()` the bare name is `rbase`, and the compartment initial condition is set as `<state>(0) <- rbase`. Per-output forms are suffixed: `lrbase_tumor`, `lrbase_anc`.
- **Source aliases:**
  - `R0`, `Base`, `BASE`, `BL`, `S0`, `TS0`, `R_baseline` -- common source-paper notations; all translate to `lrbase` without a sidecar.
- **Example models:** 187 shipped models, e.g. `Zhang_2022_ormutivimab.R`, `Hansson_2013_sunitinib_dbp.R`, `Ibrahim_2023_ibrutinib_sbp.R`.
- **Notes:** Documented 2026-09-01 as part of the Dings 2026 extraction. This entry records a long-standing convention that had never been written into the register: `lrbase` was already in consistent use across 187 model files (log baseline testosterone, baseline COWA score, baseline endogenous FVIII,...) and was cross-referenced from `le0`'s Notes as though registered. The deprecated legacy forms `lr0`, `lbl`, `lbase`, `lBase`, `ls0`, `lts0` should be migrated to `lrbase`. Distinct from `lemax` (an effect *increment*) and from `lrmax_<output>` (the plateau *level* reached under maximal drug effect); `lrbase` and `lrmax_<output>` are the two ends of the same scale, and `emax = rmax - rbase` relates them.

### lemax, limax (**canonical log-transformed maximal drug-effect increment**)
- **Type:** paper-named-param
- **Role:** Log-transformed maximum achievable drug effect in a sigmoid Emax / Imax PD model, expressed as the *increment* (or fractional increment) added to -- or removed from -- the baseline, never as an absolute attained level. Inside `model()` the bare names are `emax` / `imax`. Per-output forms are suffixed (`lemax_tumor`).
- **Source aliases:**
  - `lEmax`, `lEMAX`, `lImax` -- equivalent paper notation.
- **Example models:** 181 shipped models, e.g. `deVriesSchultink_2018_trastuzumab_LVEF.R`, `Kleijn_2011_sugammadex_rocuronium.R`, `Boucher_2018_naproxen_mbma.R`.
- **Notes:** Documented 2026-09-01 alongside `lrbase`, for the same reason -- an established convention (181 models) that had no register entry. Pairs with `lec50` (half-maximal concentration) and `lhill` (Hill exponent). **Do not use `lemax` when the source paper estimates the attained plateau level rather than the increment** -- that is `lrmax_<output>`; the two differ by the baseline and silently substituting one for the other misstates every prediction whose baseline is not zero.

### lrmax_map, lrmax_sbp, lrmax_hr (**canonical log-transformed maximum attainable plateau level of a PD output**)
- **Type:** paper-named-param
- **Role:** Log-transformed maximum attainable *level* of a PD endpoint under maximal drug effect, in the endpoint's own clinical units (mmHg, beats/min,...) -- the plateau-level counterpart of `lrbase`. Used when the source paper estimates the ceiling of the response rather than the increment, so the Emax term must be recovered as a difference: `emax_map <- rmax_map - MAP_BL`. Inside `model()` the bare names are `rmax_map`, `rmax_sbp`, `rmax_hr`; the per-output suffix follows the registered multi-analyte precedent.
- **Source aliases:**
  - `MAXMAP`, `MAXSBP`, `MAXHR` -- Dings 2026 Table A1 parameter labels.
- **Example models:** `Dings_2026_cafedrine_theodrenaline_ephedrine.R` (`rmax_hr` = 77.8 beats/min, `rmax_map` = 119 mmHg, `rmax_sbp` = 169 mmHg; each enters `endpoint = BL + (rmax - BL) * conc/(conc + ec50) + ...`).
- **Notes:** The plateau-level sibling of `lrbase`; `lemax` is explicitly **not** a substitute because it denotes an increment, so a model that reported `lemax = log(77.8)` for a 77.8 beats/min ceiling above an 84 beats/min baseline would predict a rise where the paper predicts a fall. A signed consequence worth noting: because the ceiling can sit *below* the baseline, the recovered `emax` may be negative (in the founding model `rmax_hr - HR_BL` = 77.8 - 84 = -6.2 beats/min for cafedrine/theodrenaline, and +5.5 for ephedrine after its +15% effect) -- which is exactly the paper's finding that one drug is heart-rate-neutral. Register additional `lrmax_<output>` members as new endpoints appear rather than reusing an existing one.

### lktr, lktr_slow, lktr_fast (**canonical log-transformed first-order transit / effect-delay rate constant**)
- **Type:** paper-named-param
- **Role:** Log-transformed first-order rate constant of a transit or effect-delay cascade (1/time). Inside `model()` the bare names are `ktr`, `ktr_slow`, `ktr_fast`. For a chain of `n` stages the mean transit time is `n / ktr`. `lktr` is the single-chain form; `lktr_slow` / `lktr_fast` are the qualified forms used when a model carries two parallel cascades of different speed, and they pair one-to-one with the `effect_slow<n>` / `effect_fast<n>` compartment families.
- **Source aliases:**
  - `lKTR`, `lktr1`, `lktr2` -- equivalent paper notation; `Ktr1` / `Ktr2` in Dings 2026 Table A1 map to `lktr_slow` / `lktr_fast` respectively (the paper numbers them by equation order, not by speed, so translate by the reported transit time rather than by the digit).
- **Example models:** `Dings_2026_cafedrine_theodrenaline_ephedrine.R` (`ktr_slow` = 0.107/min over three stages, mean transit time 28.0 min; `ktr_fast` = 1.76/min over four stages, mean transit time 2.27 min).
- **Notes:** `lktr` documented and the `_slow` / `_fast` qualified pair registered 2026-09-01 alongside the Dings 2026 extraction. Prefer the qualified names over numeric suffixes whenever the two chains differ in *role* (onset vs persistence), because a paper's numbering of its chains carries no guarantee of speed ordering -- in the founding paper `Ktr1` is the slower constant. The mean-transit-time relation is `n / ktr`, **not** `(n + 1) / ktr`; the latter applies only to rxode2's `transit()` closed form, where the dosing depot counts as a stage. Pairs with `lmtt` when a source parameterises by mean transit time instead.

### lkdel (**canonical log-transformed rate constant of a procedure- or intervention-onset effect**)
- **Type:** paper-named-param
- **Role:** Log-transformed first-order rate constant (1/time) governing how quickly the physiological effect of a non-drug clinical intervention (spinal anaesthesia, a surgical manoeuvre, a device) develops toward its full magnitude. Inside `model()` the bare name is `kdel`. Distinct from `ke0`, which is a plasma-to-effect-site *equilibration* rate for a drug.
- **Source aliases:**
  - `Kdel` -- Dings 2026 Table A1 ("Rate of BP effect of anesthesia").
- **Example models:** `Dings_2026_cafedrine_theodrenaline_ephedrine.R` (`kdel` = 0.171/min; drives the decaying `effect_anaesthesia` state whose complement `exp(-kdel * T_ANAESTHESIA) - effect_anaesthesia` carries the progressive spinal-anaesthesia fall in blood pressure that continues after the antihypotensive bolus).
- **Notes:** The paired magnitude parameter is the additive `eff_anaes_bp` (mmHg). Reuse only for the same intervention-onset role; a drug's effect-compartment equilibration stays `lke0`.

### kel_exp_famp (**canonical fractional amplitude of a time-varying elimination rate constant**)
- **Type:** paper-named-param
- **Role:** Dimensionless fractional amplitude of an exponential change in the elimination rate constant over a course of therapy: `kel(t) = kel * (1 + kel_exp_famp * (1 - exp(-kel_exp_kdes * t)))`, so `kel(0) = kel`, `kel(inf) = kel * (1 + kel_exp_famp)`, and a negative value means `kel` (and hence `CL`) falls over time.
- **Source aliases:**
  - `dk` -- used in `Schreib_2024_busulfan.R` (paper `theta_dk1`).
- **Example models:** `Schreib_2024_busulfan.R` (`kel_exp_famp = -0.167`, a 16.7% average fall in busulfan `k` and `CL` at steady state; carries its own additive IIV `etakel_exp_famp` and an additive covariate effect `e_hlhxlp_kel_exp_famp = -0.145` that takes the amplitude to -0.312).
- **Notes:** `_famp` is the fractional-amplitude role token, introduced because the registered `cl_exp_` time-varying-clearance family parameterises by an absolute decaying component (`cl_exp_inf` + `cl_exp_component`) rather than by a fraction of the baseline. It is additive in its covariate and random-effect terms for the same reason. Prefer this fractional parameterisation over reparameterising into `cl_exp_inf` / `cl_exp_component` whenever the source estimates `V` and `kel` separately with *different* covariate sets, because the absolute-component form is a product of the two and entangles their covariates -- in the founding model the infusion-duration covariate has opposite signs on `ln(V)` and `ln(k)`, so the product would cancel it. The sibling decay rate is `kel_exp_kdes` / `lkel_exp_kdes` in the bare and log-transformed PK sections. The parallel `cl_exp_famp` was explicitly considered and NOT registered (sidecar option C, declined); add it only when a source paper needs it.

### cl_circ_famp_day (**canonical daytime circadian fractional amplitude on clearance**)
- **Type:** paper-named-param
- **Role:** Dimensionless *signed* fractional amplitude of the daytime half of a piecewise-sine circadian rhythm multiplying a clearance term: `cl(t) = cl * (1 + cl_circ_famp_day * sin(tod / interval_day * pi))` over the daytime window, so the rhythm peaks at `cl * (1 + cl_circ_famp_day)` at the midpoint of the daytime window and returns to `cl` at both ends of it. The window lengths (`interval_day`, `interval_night`) are structural constants set by the source paper's day / night split, not estimated parameters.
- **Source aliases:**
  - `Circadian rhythm during daytime (%)` -- Chen 2025 Table 3 row label (3.70, reported as a percentage and as a magnitude only; the sign is fixed by the model code and by the Results statement that GFR and nCTS fluctuate between 104% and 91.6% of the mean).
  - `CIR_DAY` / `THETA(13)` -- Chen 2025 NONMEM control-stream token.
- **Example models:** `Chen_2025_iohexol_creatinine.R` (founding example; `cl_circ_famp_day = 0.0370398`, a +3.70% daytime peak on both the glomerular-filtration clearance and the net-tubular-secretion arm, which share a single amplitude pair).
- **Notes:** Named parameter-first / role-token-last to match the one existing precedent for the `_famp` role token, `kel_exp_famp`. The `circ_` token marks *periodic* structure and is deliberately kept distinct from the monotone time-varying-clearance families `cl_time_*` and `cl_exp_*`. The same rhythm applied to a volume, a production rate, or a biomarker baseline would generalise as `v_circ_famp_day`, `ksyn_circ_famp_day`, etc. Distinct from `lamp` / `ltacro` (Gonzalez-Sales 2015 testosterone circadian IDR), which parameterise an ABSOLUTE amplitude in the state's own units together with an estimated acrophase; here the phase is fixed by the day / night split and the amplitude is a fraction of the typical value. Pairs with [[cl_circ_famp_night]].

### cl_circ_famp_night (**canonical nighttime circadian fractional amplitude on clearance**)
- **Type:** paper-named-param
- **Role:** Dimensionless signed fractional amplitude of the nighttime half of a piecewise-sine circadian rhythm on clearance, the night-window partner of [[cl_circ_famp_day]]: `cl(t) = cl * (1 - cl_circ_famp_night * sin((tod - interval_day) / interval_night * pi))` over the nighttime window, so the rhythm troughs at `cl * (1 - cl_circ_famp_night)` at the midpoint of the night window. The two amplitudes are separate parameters because published diurnal renal-function rhythms are asymmetric -- the nocturnal fall is larger than the daytime rise.
- **Source aliases:**
  - `Circadian rhythm during nighttime (%)` -- Chen 2025 Table 3 row label (8.42, magnitude only).
  - `CIR_NIGHT` / `THETA(14)` -- Chen 2025 NONMEM control-stream token.
- **Example models:** `Chen_2025_iohexol_creatinine.R` (founding example; `cl_circ_famp_night = 0.0842088`, a -8.42% nocturnal trough).
- **Notes:** The minus sign lives in the model equation, not in the stored value, matching the source: Chen 2025 Table 3 prints both amplitudes as positive percentages and the deposited control stream carries the polarity (`CIR = -SIN((TIME-14)/10*PI)*THETA(14)+1` for the night branch versus `CIR = SIN(...)*THETA(13)+1` for the day branch). Storing the magnitude keeps the extracted value identical to the printed table; a reader who wants the signed form reads it off the `model` block. The function is continuous at both window boundaries because `sin(0) = sin(pi) = 0`, so no solver discontinuity is introduced. See [[cl_circ_famp_day]] for the family rationale.

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
- **Role:** Second-order association (forward / on) rate constant for reversible drug-target or drug-drug complex formation, with units of (1 / (concentration * time)). The paired dissociation rate constant is `k2` and the equilibrium dissociation constant is `kd = k2 / k1`. Used by mechanistic 2-drug interaction models that need separate on / off rates rather than a single steady-state `kd`. When `kd` is fixed from an external (in-vitro) source and `k2` is estimated, `k1 = k2 / kd` is derived inside `model` rather than estimated directly.
- **Source aliases:**
  - `k_on`, `kon` -- general kinetics notation.
- **Example models:** `Kleijn_2011_sugammadex_rocuronium.R` (sugammadex-rocuronium dynamic interaction: `kd` fixed at 0.0559 uM from in-vitro, `log(k2)` estimated, `k1 = k2 / kd = 0.61 1/(min*uM)` derived).
- **Notes:** Distinct from `kss` (the quasi-steady-state binding parameter used in QSS-TMDD approximations) and from `kd` (the equilibrium dissociation constant). Use the `k1` / `k2` / `kd` triple when the source paper reports the full dynamic-interaction kinetics; use `kss` alone for QSS-TMDD shortcuts.

### k2 (**canonical dissociation rate constant in reversible binding**)
- **Type:** paper-named-param
- **Role:** First-order dissociation (reverse / off) rate constant for reversible drug-target or drug-drug complex formation (1 / time). Paired with `k1` (association) and `kd = k2 / k1` (equilibrium dissociation). Inside `model` the bare name `k2` is the rate constant; the log-transformed `lk2` form is used in `ini` whenever the source paper reports `log_e(k2)` directly.
- **Source aliases:**
  - `k_off`, `koff` -- general kinetics notation.
- **Example models:** `Kleijn_2011_sugammadex_rocuronium.R` (sugammadex-rocuronium dynamic interaction: `lk2 = -3.38`, giving `k2 = 0.034 1/min`, RSE 16.5%; paired with the fixed `kd` and the derived `k1`).
- **Notes:** Pairs with `k1` (association). The paper-numerical convention `k_d = k_2 / k_1` should be encoded explicitly in `model` so the dimensionally-correct relationship is visible at the call site.

### ks (**canonical drug-mediated effect-compartment elimination rate**)
- **Type:** paper-named-param
- **Role:** Second-order rate constant for one drug's modulation of another drug's elimination from an effect compartment, with units of (1 / (concentration * time)). Used in two-drug PK-PD interaction models where the modulating drug's plasma concentration multiplies a target drug's effect-compartment amount to drive an additional elimination route (mechanistic abstraction for site-of-action drug-drug interaction such as the sugammadex-mediated rocuronium reversal). Inside `model` the bare name is `ks`; the log-transformed `lks` form is used in `ini` when the source paper reports `log_e(ks)`.
- **Source aliases:** none.
- **Example models:** `Kleijn_2011_sugammadex_rocuronium.R` (`lks = -3.43`, giving `ks = 0.033 1/(min*uM)`, RSE 0.222%; modulates rocuronium effect-compartment elimination by sugammadex plasma concentration: `d/dt(effect_roc) =... - ks * csug * effect_roc`).
- **Notes:** Mechanistically distinct from `kel` / `kdeg` (single-drug elimination rates) and from `kint` (target-mediated internalisation). The second drug's plasma concentration provides the multiplier; the parameter encodes the rate at which the modulating drug consumes the target drug at its site of action.

### ke0 (**canonical effect-compartment equilibration rate constant**)
- **Type:** paper-named-param
- **Role:** First-order rate constant for equilibration between the central plasma concentration and the effect-compartment concentration (1 / time), used by standard hysteresis effect-compartment PK-PD models: `d Ce / dt = ke0 * (Cc - Ce)`. The corresponding equilibration half-life is `log(2) / ke0`. Inside `model` the bare name is `ke0`; the log-transformed `lke0` form is used in `ini` when the source paper reports `log_e(ke0)` or uses an exponential typical-value form.
- **Source aliases:**
  - `keo`, `Keo` -- equivalent paper notation; both spellings (`keo` and `ke0`) appear in the literature.
  - `Kcb` -- Ojara 2024 notation for the plasma-to-breast-milk equilibration rate constant.
- **Example models:** `Kleijn_2011_sugammadex_rocuronium.R` (`lke0 = log(0.134) = -2.01`, giving `ke0 = 0.134 1/min` for the rocuronium effect-compartment equilibration; allometric scaling `(BW/70)^-0.25`), `Savic_2010_warfarin.R` (founding registry example -- 1-compartment PK with effect-compartment-driven proportional-odds PD), `Schindler_2016_sunitinib.R` (effect-compartment equilibration with daily AUC at ln(2)/50 1/h), `Ojara_2024_lamivudine.R` (`lke0 = log(0.245)` for plasma-to-breast-milk equilibration; the effect compartment is the physiologic `milk` matrix rather than a latent effect site), `Sun_2025_paclitaxel.R` (`lke0 = log(0.00022)`; see the degenerate-form note below).
- **Notes:** Distinct from `lke_kpd` (which is K-PD elimination rate and was deprecated in favour of the canonical `lkel`). `lke0` is specifically the effect-compartment equilibration parameter for hysteresis PK-PD modelling. **Degenerate pure-integrator form.** A minority of papers print the effect-compartment equation without its loss term, i.e. `d Ce / dt = <driving concentration> * ke0` rather than `ke0 * (Cc - Ce)`, which makes the effect compartment accumulate monotonically and never equilibrate. The name is retained because those papers themselves report the parameter as `ke0` and it plays the same role -- the constant linking plasma to the effect site -- but the parameter is then an input *gain* rather than a first-order rate constant, and its units follow the paper's own equation instead of being 1 / time. `Sun_2025_paclitaxel.R` is the registered instance: `d ce / dt = cact * ke0` with `cact` in ng/mL and `ke0` printed in mL/ng (Sun 2025 Equation 3, Table 2), a combination that is not dimensionally self-consistent as published; the value is applied verbatim and the resulting effect-compartment magnitude is verified against the paper's own reported value in that model's vignette. When encountering this form, keep `ke0`, record the printed units in the `label`, and document the divergence in the model file and the vignette Errata rather than inventing a new canonical.

### ppc (**canonical effect-compartment pseudo-partition coefficient**)
- **Type:** paper-named-param
- **Role:** Steady-state ratio of the effect-compartment concentration to the central (plasma) concentration in a Sheiner-style effect-compartment model, `d Ce / dt = ke0 * (ppc * Cc - Ce)` (unitless fraction). Unlike a plain effect compartment, where the steady-state ratio is 1 by construction, `ppc` scales the driving concentration so the effect compartment equilibrates to a fraction (or multiple) of plasma -- the standard way to describe tissue penetration when only a delayed, partitioned site-of-action concentration is measured. Inside `model` the bare name is `ppc`; the log-transformed `lppc` form is used in `ini`, since the coefficient is strictly positive.
- **Source aliases:**
  - `PPC` -- Abdelgawad 2024 paper and control-stream notation (`PPC plasma-CSF`).
  - `Kp,uu`-style "penetration ratio" / "partition coefficient" prose in CSF- and tissue-penetration papers, when the quantity is the steady-state ratio of an effect compartment rather than a physiologic tissue-to-plasma partition.
  - `Rcb` -- Ojara 2024 notation for the lamivudine milk-to-plasma concentration ratio, the gain of a breast-milk partitioned effect compartment (`dA3/dt = Kcb * (Rcb * A2/VC - A3)`).
- **Example models:** `Abdelgawad_2024_linezolid.R` (founding example -- `lppc = log(0.365)`, the maximal CSF-to-plasma pseudo-partition coefficient for linezolid in tuberculous meningitis, modulated by CSF total protein through the covariate-effect parameter `e_csftpro_ppc`), `Ojara_2024_lamivudine.R` (`lppc = log(1.77)`, the milk-to-plasma concentration ratio of a breast-milk partitioned effect compartment; the first `ppc` above 1, exercising the "or multiple" case).
- **Notes:** Distinct from `ke0`, which sets how fast the effect compartment equilibrates, whereas `ppc` sets the level it equilibrates to; a model with a partitioned effect compartment needs both. Distinct from the PBPK `kpu<n>` family (clustered unbound tissue-to-plasma partition coefficients derived from tissue composition), which are physiologic constants of a whole-body model rather than an estimated effect-compartment parameter. Covariate effects on `ppc` follow the standard `e_<cov>_ppc` shape. **Lactation models: `ppc` is the right canonical when a milk-to-plasma ratio is the gain of a partitioned effect compartment**, and is distinct from the two other lactation ratio canonicals -- use `cmpr`/`lcmpr` when the milk observable is derived purely algebraically (`Cmilk <- cmpr * Cc`, no milk state), and `pcmilk` when the paper supplies a mechanistic partition coefficient entering a central-to-milk micro-rate constant (`k = (q_milk / vc) * pcmilk`). All three generate a milk:plasma concentration ratio at pseudo-equilibrium; the discriminator is which structure the paper actually parameterised.

### lec50 (**canonical log-transformed effect-compartment EC50**)
- **Type:** paper-named-param
- **Role:** Log-transformed concentration producing half-maximal effect in sigmoid Emax / Imax PD models (concentration units). Inside `model` the bare name is `ec50`.
- **Source aliases:**
  - `lEC50`, `lec_50` -- equivalent paper notation.
- **Example models:** `deVriesSchultink_2018_trastuzumab_LVEF.R`, `Kleijn_2011_sugammadex_rocuronium.R`, `Zhang_2022_ormutivimab.R`.
- **Notes:** Pairs with `lhill` (Hill exponent). The bare name `ec50` is for use in `model` derivations.

### ec50 (**canonical bare effect-compartment EC50**)
- **Type:** paper-named-param
- **Role:** Bare counterpart of `lec50`; the half-maximal effect concentration on the linear scale.
- **Source aliases:** none.
- **Example models:** widespread sigmoid-Emax PD extractions.

### tvec50 (**canonical time-zero value of a time-varying EC50**)
- **Type:** paper-named-param
- **Role:** Value of a drifting potency `ec50t` at time zero, i.e. the intercept of the
  sigmoidal-in-time EC50 form
  `ec50t <- tvec50 + ec50_time_max * t^ec50_time_hill / (ec50_t50^ec50_time_hill + t^ec50_time_hill)`.
  Concentration units. Log form `ltvec50`. Use plain [[ec50]] when potency is constant;
  reach for this family only when the source paper estimates a drift of EC50 with
  treatment duration (pharmacodynamic tolerance / attenuation of response).
- **Source aliases:**
  - `tvEC50` -- used in `Zhang_2016_burosumab.R` (Zhang 2016 Table 3, "tvEC50, EC50,t at time 0").
- **Example models:** `Zhang_2016_burosumab.R`.
- **Notes:** The derived time-varying quantity itself is written `ec50t` inside `model()`,
  following the earlier `Jeon_2013_interferonAlfa2a.R` precedent. The `ec50_time_*` / `ec50_t50`
  trio mirrors the registered [[cl_time_max]] / [[cl_t50]] / [[cl_time_hill]] family, on EC50
  instead of clearance: the family name states the STRUCTURE (a sigmoidal-in-time change) and
  `hill` appears only on the shape coefficient. Note the two families differ in how the sigmoid
  enters -- the clearance family is multiplicative on log-clearance (`exp(...)`), whereas
  Zhang 2016 equation 10 is additive on EC50 -- so read the model body, not just the names.

### ec50_time_max (**canonical maximum magnitude of a sigmoidal-in-time EC50 change**)
- **Type:** paper-named-param
- **Role:** Change in EC50 approached as `t >> ec50_t50` in the sigmoidal-in-time EC50 form,
  in the same concentration units as [[tvec50]]. Log form `lec50_time_max`.
- **Source aliases:**
  - `a` -- used in `Zhang_2016_burosumab.R` (Zhang 2016 Table 3 row `a`, tabulated in ng/mL/week
    and glossed "maximum rate of increase of EC50,t"; in the paper's equation 10 it multiplies a
    dimensionless saturating function of time, so it is an increment rather than a rate. The
    vignette gate reproduces all three published EC50,t values under the increment reading).
- **Example models:** `Zhang_2016_burosumab.R`.
- **Notes:** Sibling of [[cl_time_max]]; see [[tvec50]] for the family rationale.

### ec50_t50 (**canonical half-time of a sigmoidal-in-time EC50 change**)
- **Type:** paper-named-param
- **Role:** Time at which half of [[ec50_time_max]] has been reached, in the same sigmoidal
  form. Log form `lec50_t50`. Units are whatever time unit the source equation uses, which
  need not be the model's time unit -- convert inside `model()` and say so.
- **Source aliases:** none.
- **Example models:** `Zhang_2016_burosumab.R` (fixed at 32 weeks, a constant printed only
  inside Zhang 2016 equation 10 and absent from Table 3; the model converts the solver's days
  to weeks before evaluating the sigmoid).
- **Notes:** Named without a `time` infix because `t50` already denotes a time, matching
  [[cl_t50]].

### ec50_time_hill (**canonical Hill coefficient of a sigmoidal-in-time EC50 change**)
- **Type:** paper-named-param
- **Role:** Sigmoid steepness exponent of the same time-varying EC50 form. Log form
  `lec50_time_hill`.
- **Source aliases:**
  - `g` / `gamma` -- used in `Zhang_2016_burosumab.R` (Zhang 2016 Table 3, glossed "Hill
    coefficient"; it is the exponent on TIME in equation 10, not on concentration).
- **Example models:** `Zhang_2016_burosumab.R`.
- **Notes:** Distinct from [[hill]], which is the sigmoidicity exponent on CONCENTRATION in a
  sigmoid Emax / Imax term. A model may carry both. Sibling of [[cl_time_hill]].

### lki50 (**canonical log-transformed half-maximal upstream-biomarker level driving a downstream effect rate**)
- **Type:** paper-named-param
- **Role:** Log-transformed half-maximal point of a sigmoid whose DRIVER is an upstream biomarker response rather than a drug concentration, entering a downstream effect-RATE expression: `rate = rmax * I^hill / (ki50^hill + I^hill)`, where `I` is the biomarker response on its own axis (percent inhibition of a phosphoprotein, percent receptor occupancy, fold-change of a transcript). Carries the units of that biomarker axis, NOT concentration units. Inside `model` the bare name is `ki50`.
- **Source aliases:**
  - `KI50` -- Moein 2024 notation ("I where Ks is 50% of Kmax", Tables 3 and 4); constrained to (1, 100) in both NONMEM control streams because the driver is a percentage.
- **Example models:** `Moein_2024_apitolisib_mouse.R` (founding example; `ki50 = 67.2` %pAkt inhibition for the 786-O xenograft tumour-shrinkage sigmoid), `Moein_2024_apitolisib_human.R` (`ki50 = 58.0` %pAkt inhibition for the phase 1 RECIST sum-of-longest-diameters fit).
- **Notes:** Distinct from `lec50` / `lic50`, which are DRUG CONCENTRATIONS producing half-maximal effect -- a model can carry both at once, as the Moein files do (`ic50` on the drug-to-pAkt step, `ki50` on the pAkt-to-tumour step of the same cascade). Distinct also from `kc50` as used for a trimer CONCENTRATION at half-maximum kill rate: `ki50` is deliberately axis-agnostic and is the name to reach for when the sigmoid's x-axis is a biomarker readout. Pairs with `kmax` (the maximum of the driven rate) and `hill` (the sigmoidicity factor). A downstream user must read the axis units off `label()` / `compartmentData`, since they are model-specific by construction.

### ki50 (**canonical bare half-maximal upstream-biomarker level driving a downstream effect rate**)
- **Type:** paper-named-param
- **Role:** Bare counterpart of `lki50`; the half-maximal biomarker-response level on the linear scale, for use inside `model`.
- **Source aliases:**
  - `KI50` -- Moein 2024 notation.
- **Example models:** `Moein_2024_apitolisib_mouse.R`, `Moein_2024_apitolisib_human.R`.

### lreg_qm (**canonical log-transformed once-monthly regimen EC50 multiplier**)
- **Type:** paper-named-param
- **Role:** Log-transformed multiplicative modifier applied to EC50 when a static Emax-on-AUC exposure-response model parameterised on an integrated-window AUC predictor (rather than an instantaneous concentration) needs to distinguish between dosing regimens whose AUCs encode different time-course profiles. Paired with the binary covariate `REGI_QM`: applied as `ec50_eff = ec50 * reg_qm^REGI_QM`, i.e., `ec50_eff = ec50` for Q2W (REGI_QM = 0) and `ec50_eff = ec50 * exp(lreg_qm)` for QM (REGI_QM = 1). Captures the "lack of kinetics" gap in static Emax-on-AUC models where two regimens that produce similar window-AUC values can still produce different effects because of differences in target-saturation time courses. Inside `model` the bare name is `reg_qm`.
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
- **Role:** Log-transformed baseline (drug-free) value of a PD response, e.g., baseline T4/T1 twitch ratio in neuromuscular blockade PK-PD models (`E0 = exp(le0)`). Distinct from `lrbase`, which is the log-transformed steady-state baseline for indirect-response / turnover models with explicit `kin` / `kout`. Use `le0` when the source paper reports a non-turnover PD baseline (a typical T4/T1 ratio, a typical pre-treatment biomarker level) that enters the sigmoid Emax expression as an additive baseline plus the maximum-effect-bounded modulation. Inside `model` the bare name is `e0`.
- **Source aliases:**
  - `lE0`, `lTOF0` -- equivalent paper notation.
- **Example models:** `Kleijn_2011_sugammadex_rocuronium.R` (`le0 = log(104)` for typical baseline T4/T1 x 100; Emax of the sigmoid is forced equal to E0 so the per-subject baseline shape of the readout is preserved).
- **Notes:** Distinct from `lrbase` (turnover-state steady-state baseline). When `Emax = E0` is forced in the sigmoid (the standard NMB parameterisation), the readout decreases monotonically from `E0` toward 0 as the effect-compartment concentration rises.

### e0 (**canonical bare PD baseline parameter**)
- **Type:** paper-named-param
- **Role:** Bare counterpart of `le0`; the baseline (drug-free) PD response value used inside `model`.
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

### kphag (**canonical phagocytic uptake rate**)
- **Type:** paper-named-param
- **Role:** First-order rate constant for phagocytic uptake of a target species by host phagocytes (1 / time). Covers both shapes the library has met so far. (a) **Threshold-gated elimination**: phagocytic clearance of antibody-coated target cells (platelets, RBCs, etc.) in TMDD models where uptake is gated by a receptor-occupancy threshold -- `kphag` acts on both the free-target and the drug-target-complex states while receptor occupancy exceeds the paper-fixed `thres` (percent), and is switched off below it. (b) **Ungated transfer**: always-on phagocytosis of bacteria by granulocytes in semi-mechanistic infection PKPD models, where the rate moves bacteria from the extracellular states into a still-counted phagocytosed compartment rather than eliminating them, and a separate digestion rate empties that compartment. Distinct from the constitutive `kdeg` (baseline turnover of a synthesis/degradation pool) and from `kint` (receptor-mediated complex internalisation). The log-transformed `lkphag` form is used in `ini()` when the rate is estimated directly.
- **Source aliases:**
  - `KPH` -- Moc Willeford 2024 NONMEM abbreviation (supplement differential equations).
  - `kphagocytosis` -- Moc Willeford 2024 Table 1 / 2 typeset name.
  - `k_phag` -- Saporta 2026 Figure 2 and equations 1-3 typeset name.
- **Example models:** `MocWilleford_2024_rlyb212.R` (RLYB212 anti-HPA-1a IgG1 with threshold-gated phagocytic elimination of HPA-1a-positive platelets; founding example of shape (a), `lkphag` estimated directly in `ini()`), `Saporta_2026_meropenem.R` (founding example of shape (b): ungated granulocyte phagocytosis of *Klebsiella pneumoniae* into a `bact_phagocytosed` state in a mouse lung infection model).
- **Notes:** Paired with the canonical `thres` receptor-occupancy threshold **in shape (a) only** -- the threshold gate is a Moc Willeford-specific feature of that model, not part of this canonical's definition. Broadened from "threshold-gated phagocytic elimination rate" to cover ungated phagocytic uptake on 2026-09-01 (operator sidecar `oare_PMC13041408` request-001 / response-001, question q1, option A); the previous Notes directed the ungated case to `kdeg`, which is wrong whenever the uptake is a transfer into a counted compartment rather than a degradation. A paper may parameterize `kphag` indirectly: `Saporta_2026_meropenem.R` carries it as a sum of additive immune-state increments (`irNeu`, `irInt`, `irCom` -- the paper's `IR_neu` / `IR_int` / `IR_com`) selected by the `IMMUNE_STATE` covariate, and derives `kphag` inside `model()`. When a paper reports a separate phagosomal digestion rate that it constrains to equal `kphag`, write it as the derived local `kdig <- kphag` inside `model()` rather than registering a canonical for a parameter that was never independently estimated (same sidecar, question q2, option A).

### thres (**canonical receptor-occupancy threshold**)
- **Type:** paper-named-param
- **Role:** Threshold value of receptor occupancy (typically expressed as a percent) above which a threshold-gated elimination pathway (e.g., `kphag`) is switched on. Below the threshold, the gated pathway contributes zero to the target elimination rate; above the threshold, it contributes its full rate constant. Used in TMDD models where receptor engagement must exceed a minimum coating fraction before macrophage / RES-mediated phagocytic clearance is triggered.
- **Source aliases:**
  - `THRES` -- Moc Willeford 2024 NONMEM / typeset name.
- **Example models:** `MocWilleford_2024_rlyb212.R` (THRES = 10 percent fixed, describing the minimum receptor coating required to drive phagocytic clearance of HPA-1a-positive platelets; founding example).
- **Notes:** Usually held fixed at a value inferred from the paper's model-development results (Moc Willeford 2024: estimation ranged 3-22 percent when co-estimated with other parameters, so THRES was fixed at 10 percent). Distinct from a generic sigmoidal `hill` or `ec50` because the pathway is a hard switch rather than a smooth curve. When expressed as a percent, `thres` is compared to a percent receptor occupancy `ro = 100 * complex / (target + complex)`; when expressed as a fraction, `thres` is compared to the corresponding unitless ratio. `thres` is deliberately scoped to a receptor-occupancy threshold and is NOT the canonical for a threshold expressed as a drug *concentration* -- see the sibling `cthres` for that. Widening `thres` to "the threshold of whatever quantity gates a hard switch, in that quantity's own units" was considered and rejected when `cthres` was ratified, because it would make one canonical span two different physical dimensions (percent versus ng/mL) and force a reader to open the model to learn which one a given `thres` carries.

### cthres (**canonical gating threshold concentration**)
- **Type:** paper-named-param
- **Role:** Threshold *drug concentration* above which a gated downstream pharmacodynamic input accrues, in the model's concentration units. Below the threshold the gated term contributes exactly zero; at or above it, the gated term is the excess over the threshold. The canonical hard-switch form is `cact = max(Cc - cthres, 0)`, and the gated quantity `cact` then drives the downstream state (an effect compartment, a turnover input, or a hazard). Inside `model` the bare name is `cthres`; the log-transformed `lcthres` form is used in `ini`, since the threshold is strictly positive. Use `cthres` whenever a paper estimates a concentration below which the drug is modelled as having no pharmacodynamic effect at all.
- **Source aliases:**
  - `x`, `x_pop` -- Sun 2025 notation (Equation 1 and Table 2; the paper writes the gate as `f(x)`, so the threshold itself carries the bare symbol `x`).
- **Example models:** `Sun_2025_paclitaxel.R` (`lcthres = log(1735.75)` for the paclitaxel plasma concentration above which chemotherapy-induced peripheral neuropathy accrues; the gated excess drives a pure-integrator effect compartment feeding a CIPN8 turnover pool; founding example).
- **Notes:** The leading `c` marks that the gated quantity is a concentration, distinguishing it from the occupancy-scoped `thres` (a percent or fraction), with which it is otherwise a direct sibling; both are hard switches. Distinct from `mic` on provenance rather than form: `mic` is a *measured* susceptibility property of a pathogen isolate and is normally held `fixed`, whereas `cthres` is an *estimated* patient-level PD parameter reported with a standard error (Sun 2025: SE 1.84, RSE 0.106%). Distinct from `ec50` / `lec50` on functional form: `ec50` parameterises a smooth saturable curve in which effect accrues at every non-zero concentration, whereas `cthres` gates a discontinuous switch that contributes exactly nothing below the threshold. Do not substitute `lec50` for a hard-switch threshold; doing so misrepresents the published structure. A model whose gate is applied to an effect-site rather than a central concentration still uses `cthres` -- the canonical names the threshold, not the compartment it is compared against.

### qp (**canonical target inter-compartmental clearance**)
- **Type:** paper-named-param
- **Role:** Inter-compartmental clearance between the central and peripheral compartments of a target species (platelet, receptor, or other cell population tracked as an ODE state alongside the drug in a TMDD model). Units L / time. Distinct from the canonical drug `lq` / `q` because the two flows are semantically different: `q` moves drug between drug compartments, `qp` moves target between target compartments. Used when the source paper carries a peripheral distribution for the target species with its own inter-compartmental clearance parameter.
- **Source aliases:**
  - `QP` -- Moc Willeford 2024 Table 1 / 2 / supplement notation.
- **Example models:** `MocWilleford_2024_rlyb212.R` (QP = 2.45 L/h: platelet inter-compartmental clearance between the shared central compartment (V1) and the target peripheral compartment (V3); founding example).
- **Notes:** Encoded via `lqp` in `ini` and `qp` in `model`; scales allometrically with body weight at exponent 0.75 in the Moc Willeford 2024 simulation.

### vp_target (**canonical target peripheral volume**)
- **Type:** paper-named-param
- **Role:** Peripheral volume of distribution for a target species tracked as an ODE state alongside the drug in a TMDD model (units L). Distinct from the drug's canonical `lvp` / `vp` because the drug and the target may have physically different peripheral distribution volumes even when they share a central compartment. Paired with `qp` as the corresponding target inter-compartmental clearance.
- **Source aliases:**
  - `V3` -- Moc Willeford 2024 Table 1 / 2 / supplement notation for the target peripheral volume (the paper's V1 is the shared central, V2 is the drug peripheral, V3 is the target peripheral).
- **Example models:** `MocWilleford_2024_rlyb212.R` (V3 = 0.523 L: peripheral platelet distribution volume, paired with target central sharing V1; founding example).
- **Notes:** Encoded via `lvp_target` in `ini` and `vp_target` in `model`. Do not conflate with the drug's `lvp2` (second peripheral volume of a 3-compartment drug model); those are drug compartments, `vp_target` is a target compartment.

### kbm (**canonical biliary-metabolite excretion rate**)
- **Type:** paper-named-param
- **Role:** First-order rate constant for biliary excretion of a drug or metabolite from a plasma / central compartment into a downstream gut / bile compartment (1 / time). Used in enterohepatic-recirculation and interconversion submodels where the source paper carries biliary transport as a separate ODE flux (distinct from the drug's total clearance terms).
- **Source aliases:**
  - `k_bm` -- Hamren 2008 paper notation.
- **Example models:** `Hamren_2008_tesaglitazar.R` (kbm = 11.7 1/h: rate at which the acyl glucuronide metabolite is biliary-secreted from the metabolite plasma compartment into the gut compartment, where it then undergoes interconversion back to parent tesaglitazar; founding example).
- **Notes:** Distinct from a generic intercompartmental clearance `lq` because biliary excretion in interconversion models is one-way (drug leaves the plasma compartment for the gut and does not return via the same route -- the return path is via a separate hydrolysis / reabsorption process parameterised by `kicv`).

### kehc (**canonical gallbladder-emptying rate constant**)
- **Type:** paper-named-param
- **Role:** First-order rate constant (1 / time) at which a `gallbladder` compartment discharges its contents back towards the central or absorption compartment during a post-prandial emptying window, completing an enterohepatic-circulation loop. It is the *return* leg of the pair whose *outbound* leg is `kbm`: `kbm` moves drug central -> gallbladder continuously and independently of food, `kehc` moves it back only while a meal gate is open. Encoded via `lkehc` in `ini()` (strictly positive) and `kehc` in `model()`, and always multiplied by the meal gate rather than applied unconditionally, e.g. `d/dt(gallbladder) <- kbm * central - meal * kehc * gallbladder`.
- **Source aliases:**
  - `k_GE`, `kGE` -- Keunecke 2020 notation ("GE" = gallbladder emptying).
- **Example models:** `Keunecke_2020_regorafenib_phase1.R`, `Keunecke_2020_regorafenib_phase3.R` (`kehc` = 100 1/h, fixed, so that emptying is "complete and (almost) instantaneous" over the 0.5 h window; doi:10.1111/bcp.14334; founding example).
- **Notes:** Named `kehc` rather than `kge` because `kge` is already taken by the unrelated Stein bi-exponential growth-rate constant, and rather than a paper-numbered transfer name (`k84` in `Kim_2018_tacrolimus.R`, `k21` in `Ide_2009_pravastatin.R`) because inter-compartment numbering is not portable across papers -- the same reasoning recorded under `kbm` / `kicv`. A large fixed value paired with a short gate window is the usual parameterisation, which makes the flux a stiff, short-lived discontinuity: an event table used with such a model needs an explicit record at each gate edge or the solver will step over the window. See `tmeal1` / `dge` for the gate itself.

### tmeal1, tmeal2, tmeal3 (**canonical meal clock times gating gallbladder emptying**)
- **Type:** paper-named-param
- **Role:** Clock times, in hours since midnight, of the daily meals that trigger gallbladder emptying in an enterohepatic-circulation model. One member per meal the source paper models (`tmeal1` = breakfast, `tmeal2` = lunch, `tmeal3` = dinner); extend the numbering if a paper models more. Carried on the LINEAR scale, so a held-constant value is written `tmeal1 <- fixed(8)`. Combined with `dge` inside `model()` to build the gate that multiplies `kehc`, reading a *time of day* derived from absolute time: `tod <- t - 24 * floor(t / 24)`, then `meal <- (tod >= tmeal1) * (tod < tmeal1 + dge) +...`.
- **Source aliases:**
  - `Breakfast`, `Lunch`, `Dinner` -- Keunecke 2020 Table 1 / 2 / 3 row labels, reported as clock times (08:00, 12:00, 18:00).
- **Example models:** `Keunecke_2020_regorafenib_phase1.R`, `Keunecke_2020_regorafenib_phase3.R` (`tmeal1` / `tmeal2` / `tmeal3` = 8 / 12 / 18 h, all fixed; doi:10.1111/bcp.14334; founding example).
- **Notes:** Because these are absolute clock times rather than times after dose, **any event table used with such a model must take midnight as its origin** and place doses at the intended clock time; that requirement belongs in the model's `description` and the validation vignette. This is the distinguishing feature against `Kim_2018_tacrolimus.R`, whose `lmtime1` / `lmtime2` carry a single NONMEM-`MTIME`-style emptying window measured *post-dose* for a single-dose study; a multiple-dose model with three meals a day cannot be expressed that way, which is why a time-of-day family is registered separately rather than reusing that name. Do not use `mtime`, which would collide with rxode2's `mtime()` model-event-time function -- the same collision recorded under `tkacut`.

### dge (**canonical duration of gallbladder emptying**)
- **Type:** paper-named-param
- **Role:** Duration (time units) of the post-prandial window during which a `gallbladder` compartment empties, i.e. the width of the gate opened at each `tmeal<n>` and closed at `tmeal<n> + dge`. Carried on the LINEAR scale, so a held-constant value is written `dge <- fixed(0.5)`.
- **Source aliases:**
  - `DGE` -- Keunecke 2020 Table 1 / 2 / 3 row label ("duration of gallbladder emptying").
- **Example models:** `Keunecke_2020_regorafenib_phase1.R`, `Keunecke_2020_regorafenib_phase3.R` (`dge` = 0.5 h, fixed, matching the paper's assumption that "the duration of transfer of regorafenib from the gallbladder compartment to the central compartment was approximately 30 min"; doi:10.1111/bcp.14334; founding example).
- **Notes:** Distinct from a `dur()` / infusion-duration construct: `dge` does not describe a zero-order input of a dose, it parameterises a time gate on an internal recirculation flux. Paired with `tmeal1` / `tmeal2` / `tmeal3` and `kehc`; the three names are only meaningful together. `Kim_2018_tacrolimus.R` carries the same quantity as `lmtime2` on the log scale for a single post-dose window -- see the note under `tmeal1` for why that name was not reused.

### kicv (**canonical interconversion rate constant for metabolite-to-parent reverse-metabolism**)
- **Type:** paper-named-param
- **Role:** First-order rate constant for interconversion of a metabolite back to its parent drug (1 / time). Used for acyl-glucuronide-style reverse-metabolism mechanisms where the gut microbial / brush-border beta-glucuronidases hydrolyse a biliary-excreted conjugated metabolite back to the parent aglycone, which is then reabsorbed into systemic circulation. The single rate constant subsumes both the hydrolysis step and the reabsorption step (the source paper typically reports them as a single first-order process because the underlying gut microflora kinetics and intestinal-uptake kinetics cannot be separated from oral PK data alone).
- **Source aliases:**
  - `k_int` -- Hamren 2008 paper notation (the paper writes `kint` for interconversion; this naming clashes with the canonical `kint` for TMDD complex internalisation, so the canonical here is `kicv` to preserve semantic disambiguation).
- **Example models:** `Hamren_2008_tesaglitazar.R` (kicv = 0.79 1/h: rate at which acyl glucuronide accumulated in the gut compartment is hydrolysed and reabsorbed as parent tesaglitazar into the parent's central plasma compartment; founding example).
- **Notes:** The canonical name `kicv` was chosen specifically to disambiguate from the TMDD `kint` canonical -- the two mechanisms (drug-target complex internalisation vs. metabolite-to-parent interconversion via gut hydrolysis) are semantically distinct, and reusing the same name with paper-specific meaning would create a latent reader-confusion failure mode. Paper-numbered alternatives (`k57`, `k72`-style as in deWinter 2009 MPA / MPAG) are discouraged for new extractions because the numbering is non-portable across papers; prefer the role-based `kbm` / `kicv` pair for any future interconversion popPK model.

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
- **Role:** Proliferation / growth-rate constant of a tumor-growth-inhibition or cell-population state. Paper-mechanistic meaning depends on the growth law. Inside `model` the bare name is `p`; the log-transformed `lp` form is used in `ini` because the rate is strictly positive.
- **Source aliases:**
  - `k` -- Yates 2023 notation for the growth-law rate constant of all three laws (exponential, Mayneord, von Bertalanffy).
- **Example models:** `Yates_2023_tgi_exponential.R` (founding example; `p = 0.005` 1/day, the exponential fractional growth rate), `Yates_2023_tgi_mayneord.R` (`p = 0.005` L^(1/3)/day), `Yates_2023_tgi_bertalanffy.R` (`p = 0.005` L^(1/3)/day).
- **Notes:** `p` is the TGI-family *generic*: use it whenever a paper's growth constant is not one of the mechanistically-specific alternatives. Distinct from `kgrow` (life-cycle-anchored multiplication rate constrained by a known cycle time and per-cycle amplification factor) and from `kge` (the growing-fraction exponent of the Stein two-exponential decomposition); both of those entries cross-reference `p` as the family generic. **Dimensions are set by the growth law, not by the name.** In an exponential law `dV/dt = p * V` the constant acts on volume and carries 1 / time; in the sub-exponential Mayneord and von Bertalanffy laws `dV/dt = p * V^(2/3) -...` it acts on the tumor radius and carries length / time (L^(1/3)/day when V is in L). Record the dimensions in the parameter's `label()` rather than assuming 1 / time. Founding examples added 2026-08-05 with the Yates 2023 extraction; the entry itself predates them. Paired with `kdeath` in the von Bertalanffy law.

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
- **Role:** First-order growth / net-multiplication rate constant of a paper-mechanistic proliferating state (1 / time). Founding use: fixed net asexual-parasite growth rate `kgrow = ln(10) / 48 = 0.0479 /h` for a Plasmodium falciparum life-cycle model in which parasites multiply 10-fold per 48-h intraerythrocytic cycle (Hien 2017 Table 3 row `K_grow (1/h) = 0.0479 fix`). Inside `model` the bare name is `kgrow`; the log-transformed `lkgrow` form is used in `ini` when the rate itself is log-parameterised (typical for positivity constraint).
- **Source aliases:**
  - `K_grow`, `Kgrow` -- Hien 2017 notation.
- **Example models:** `Hien_2017_cipargamin.R` (founding example; `kgrow` fixed at ln(10)/48 = 0.0479 /h for the 10-fold per 48-h cycle Plasmodium falciparum multiplication rate).
- **Notes:** Mechanistically distinct from `p` (generic proliferation / growth-rate constant, TGI-family) in that `kgrow` is specifically a life-cycle-anchored multiplication rate constrained by an independently-known cycle time and per-cycle amplification factor -- typically fixed rather than estimated. Also distinct from `kin` / `ksyn` (zero-order production rates into a turnover pool) because `kgrow` is a first-order growth rate proportional to the state itself.

### kge (**canonical Stein bi-exponential growth-rate constant**)
- **Type:** paper-named-param
- **Role:** First-order growth-rate constant of the treatment-resistant fraction in the Stein (2011) bi-exponential tumor-dynamics model, `y(t) = y0 * (exp(-kse * t) + exp(kge * t) - 1)` (1 / time). Drives the `growth` sub-state: `d/dt(growth) = kge * growth` with `growth(0) = y0`. Paired with `kse`. Inside `model()` the bare name is `kge`; the log-transformed `lkge` form is used in `ini()` because the rate is strictly positive. Per-treatment-arm variants take an arm suffix (`lkge_pembro`, `lkge_chemo`,...); per-endpoint variants in a joint multi-biomarker model take an endpoint suffix (`kge_ctdna`).
- **Source aliases:**
  - `KG`, `TVKG`, `kg`, `k_g`, `kgT` -- Stein-family paper notation for the growth rate (`kgT` when the paper needs to distinguish a tumor-size growth rate from a second endpoint's).
- **Example models:** `Struemper_2025_tumorsize_OS_nsclc.R` (founding example; 12 per-arm `lkge_<arm>` values, 1/week), `Ribba_2022_sld.R` (`kge = 0.0016` 1/day for the OAK sum-of-longest-diameters fit).
- **Notes:** Register entry backfilled 2026-07-28; the name has been in use since Struemper 2025 and is described inside the `growth` / `shrink` entries of `compartment-names.md`, but had no entry of its own. Distinct from `kgrow` (life-cycle-anchored multiplication rate constrained by a known cycle time) and from `p` (generic TGI proliferation constant): `kge` is specifically the growing-fraction exponent of the two-exponential Stein decomposition, whose defining feature is that growth and shrinkage act on two independent sub-populations that are summed rather than on a single net-rate state. The related time-to-tumor-growth summary statistic is `TTG = (log(kse) - log(kge)) / (kge + kse)`.

### kse (**canonical Stein bi-exponential shrinkage-rate constant**)
- **Type:** paper-named-param
- **Role:** First-order decay-rate constant of the treatment-sensitive fraction in the Stein (2011) bi-exponential tumor-dynamics model (1 / time). Drives the `shrink` sub-state: `d/dt(shrink) = -kse * shrink` with `shrink(0) = y0`. Paired with `kge`. Inside `model` the bare name is `kse`; the log-transformed `lkse` form is used in `ini`. Arm- and endpoint-suffixed variants follow the same pattern as `kge`.
- **Source aliases:**
  - `KS`, `TVKS`, `ks`, `k_s`, `ksT` -- Stein-family paper notation for the shrinkage / decay rate.
- **Example models:** `Struemper_2025_tumorsize_OS_nsclc.R` (founding example; 12 per-arm `lkse_<arm>` values, 1/week), `Ribba_2022_sld.R` (`kse = 0.0014` 1/day for the OAK sum-of-longest-diameters fit).
- **Notes:** Register entry backfilled 2026-07-28 together with `kge`. **Name-collision warning:** the bare symbol `ks` used by much of the Stein literature is already taken in this register by a different canonical -- `ks` is the drug-mediated effect-compartment elimination rate of Kleijn 2011, with units 1 / (concentration * time). Always map a Stein-model paper's `ks` to `kse`, never to `ks`.

### kdeath (**canonical constitutive cell-death / loss rate**)
- **Type:** paper-named-param
- **Role:** Constitutive, drug-independent first-order cell-death / loss rate of a tumor state (1 / time). Founding use: the `- kd * V` loss term of the von Bertalanffy growth law `dV/dt = p * V^(2/3) - kdeath * V`, which is what makes an untreated tumor plateau at the finite carrying volume `(p / kdeath)^3`. Inside `model` the bare name is `kdeath`; the log-transformed `lkdeath` form is used in `ini` because the rate is strictly positive.
- **Source aliases:**
  - `kd` -- Yates 2023 notation for the von Bertalanffy natural cell-death rate.
- **Example models:** `Yates_2023_tgi_bertalanffy.R` (founding example; `kdeath = 0.004` 1/day, giving an untreated plateau volume of `(0.005 / 0.004)^3 = 1.95` L), `tgi_sat_VonBertalanffy.R`, `tgi_sat_genVonBertalanffy.R` (both templates, retrofitted from `lkd` 2026-08-06).
- **Notes:** **Deliberately not named `kd`, and deliberately not merged into `kk`.** The bare `kd` is already this register's canonical for the mechanistic *dissociation* rate constant and is in live use with that meaning (`Kleijn_2011_sugammadex_rocuronium.R`, `Hayashi_2007_omalizumab.R`, `Lowe_2009_omalizumab.R`, `Anbari_2023`); overloading it would give one canonical two incompatible roles. The templates `tgi_sat_VonBertalanffy.R` and `tgi_sat_genVonBertalanffy.R` originally used `lkd` for this concept -- they predated the register and were the collision, not the precedent; both were retrofitted to `lkdeath` / `kdeath` on 2026-08-06, so no model in the library now uses `kd` for a cell-death rate. Note that `lkd` is still correct, and widely used, for the *dissociation* constant (~18 models). Distinct also from Cardilin 2018's `kk`, which is a *drug- or radiation-modulated* kill rate rather than the constitutive loss term encoded here, and from `kdeg` (first-order degradation of a turnover pool, a synthesis/degradation-balance concept rather than a growth-law loss term). Paired with `p`.

### kge_ctdna (**canonical Stein ctDNA growth-rate constant**)
- **Type:** paper-named-param
- **Role:** Endpoint-suffixed `kge` for the circulating-tumor-DNA arm of a Stein bi-exponential biomarker model (1 / time). Drives `d/dt(growth_ctdna) = kge_ctdna * growth_ctdna`. The `_ctdna` suffix is required whenever a single model fits both a tumor-size and a ctDNA Stein pair, so the two growth rates are unambiguous. Log-transformed form `lkge_ctdna` in `ini`.
- **Source aliases:**
  - `kg` -- Ribba 2022 notation for the ctDNA growth rate in both Eq. 1 and Eq. 2.
- **Example models:** `Ribba_2022_ctdna.R` (founding example; `kge_ctdna = 0.0038` 1/day, RSE 30.1%), `Ribba_2022_ctdna_sld_joint.R` (same value, fixed).
- **Notes:** See `kge` for the underlying Stein decomposition and `growth_ctdna` in `compartment-names.md` for the paired state.

### kse_ctdna (**canonical Stein ctDNA decay-rate constant**)
- **Type:** paper-named-param
- **Role:** Endpoint-suffixed `kse` for the circulating-tumor-DNA arm of a Stein bi-exponential biomarker model (1 / time). Drives `d/dt(shrink_ctdna) = -kse_ctdna * shrink_ctdna`. May be estimated directly (log-transformed `lkse_ctdna` in `ini`) or derived inside `model` from a tumor-size decay rate via the cross-endpoint link parameter `zeta`.
- **Source aliases:**
  - `ks` -- Ribba 2022 Eq. 1 notation for the freely-estimated ctDNA decay rate. Do NOT map this to the unrelated canonical `ks` (Kleijn 2011 effect-compartment elimination rate).
- **Example models:** `Ribba_2022_ctdna.R` (founding example; freely estimated, `kse_ctdna = 0.0081` 1/day, RSE 27.4%), `Ribba_2022_ctdna_sld_joint.R` (derived as `kse_ctdna <- zeta * kse`, not estimated).
- **Notes:** See `kse` for the underlying Stein decomposition and the `ks` name-collision warning.

### zeta (**canonical cross-endpoint decay-rate link**)
- **Type:** paper-named-param
- **Role:** Dimensionless multiplier that ties one endpoint's Stein decay-rate constant to a second, jointly-modeled endpoint's decay-rate constant, so the two biomarkers share a single underlying treatment-response process rather than decaying independently. Founding use: `kse_ctdna = zeta * kse` in Ribba 2022 Eq. 2, coupling the ctDNA decay rate to the tumor-size (sum-of-longest-diameters) decay rate. `zeta > 1` means the ctDNA signal falls faster than tumor size. Inside `model` the bare name is `zeta`; the log-transformed `lzeta` form is used in `ini` because the multiplier is strictly positive, with IIV `etalzeta`.
- **Source aliases:**
  - Greek `zeta` -- Ribba 2022 Eq. 2 typeset symbol.
- **Example models:** `Ribba_2022_ctdna_sld_joint.R` (founding example; `zeta = 1.94`, RSE 37.3%, `omega_zeta = 0.86`, RSE 35.0%).
- **Notes:** Distinct from `frac` (a bounded 0-1 fraction of a pool) and from `alfm` (a mixture-model mixing proportion): `zeta` is an unbounded positive ratio between two rate constants of the same dimension, and is the identifiable free parameter of a joint fit whose remaining structural parameters are held at their single-endpoint values. Use this canonical for any joint-biomarker model that couples two endpoints through a shared rate constant scaled by an estimated factor; if a future model needs more than one such link, suffix by the linked endpoint (`zeta_<endpoint>`).

### kact (**canonical activation-rate parameter**)
- **Type:** paper-named-param
- **Role:** First-order activation rate constant for a paper-mechanistic dormant / refractory state to become active (1 / time). Founding use: rate constant `kact` at which refractory (drug-tolerant) Plasmodium falciparum parasites become active and re-enter the drug-sensitive pool in a two-population parasite clearance model (Hien 2017 Table 3 row `K_act (1/h) = 0.0987`). Inside `model` the bare name is `kact`; the log-transformed `lkact` form is used in `ini` when the rate itself is log-parameterised.
- **Source aliases:**
  - `K_act`, `Kact` -- Hien 2017 notation.
- **Example models:** `Hien_2017_cipargamin.R` (founding example; `kact = 0.0987 /h` for refractory-to-active first-order transition, with IIV 41.5% CV per Hien 2017 Table 3; drives the awakening term `+ kact * parasite_refractory` in the sensitive-pool ODE and the loss term `- kact * parasite_refractory` in the refractory-pool ODE).
- **Notes:** Mechanistically distinct from `kel` / `kdeg` (single-drug or single-pool elimination), from `kmet` (parent-to-metabolite conversion), and from `kint` (target-mediated internalisation) in that `kact` transfers mass from an inactive pool to an active pool of the same species (no drug binding, no metabolism, no target sequestration).

### fsen (**canonical bare drug-sensitive fraction**)
- **Type:** paper-named-param
- **Role:** Fraction of a paper-mechanistic asexual / cycling pool that is fully drug-sensitive at model initialisation, bounded in `[0, 1]`. Complement `(1 - fsen)` is the drug-refractory subpool. Inside `model` the bare name is `fsen`; the log-transformed `lfsen` form is used in `ini`. Founding use: population fraction of asexual Plasmodium falciparum parasites that are drug-sensitive at enrolment in a two-population parasite clearance model (Hien 2017 Table 3 row `F_sen (%) = 99.1`).
- **Source aliases:**
  - `F_sen`, `Fsen` -- Hien 2017 notation.
- **Example models:** `Hien_2017_cipargamin.R` (founding example; `fsen = 0.991` typical; large IIV 81.8% CV per Hien 2017 Table 3; seeds the initial condition of the two ODE pools: `parasite_sensitive(0) = fsen * PARA` and `parasite_refractory(0) = (1 - fsen) * PARA`).
- **Notes:** Encoded on the log scale (`lfsen = log(0.991)`) for consistency with the paper's log-normal `omega^2 = log(1 + CV^2)` reporting formula (Hien 2017 Table 3 footnote a). Individual `fsen_i = exp(lfsen + etalfsen)` can occasionally exceed 1 for extreme etas; this is a documented limitation of the log-normal-on-fraction encoding and is preferable to reinterpreting the paper's reported %CV on a logit scale. Downstream simulations that need strictly-bounded individual fractions can clamp `fsen_i` to `min(fsen_i, 1)` at post-processing. Distinct from `fu` (fraction unbound in plasma; time-invariant physicochemical / binding property), `fr` (Bergstrand-Karlsson mixed-model fraction of MAT in transit delay), and `fm` (fraction metabolised through a specific pathway).

### lfsen (**canonical log-transformed drug-sensitive fraction**)
- **Type:** paper-named-param
- **Role:** Log-transformed counterpart of `fsen` for use in `ini`. Individual fraction `fsen_i = exp(lfsen + etalfsen)`.
- **Source aliases:** none.
- **Example models:** `Hien_2017_cipargamin.R`.
- **Notes:** See `fsen` for the full role description and the caveat on the log-transform-on-a-bounded-fraction encoding choice.

### lkgrow (**canonical log-transformed growth-rate parameter**)
- **Type:** paper-named-param
- **Role:** Log-transformed counterpart of `kgrow` for use in `ini`.
- **Source aliases:** none.
- **Example models:** `Hien_2017_cipargamin.R` (`lkgrow = fixed(log(0.0479))`).

### lkact (**canonical log-transformed activation-rate parameter**)
- **Type:** paper-named-param
- **Role:** Log-transformed counterpart of `kact` for use in `ini`.
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
- **Role:** Number of absorption transit compartments in a Savic (2007) transit chain (unitless, and typically non-integer because it is estimated as the shape parameter of the equivalent gamma-density input). Closes the `mtt` / `ktr` / `ntr` family: `ktr = (ntr + 1) / mtt`. Inside `model` the bare name is `ntr`; the log-transformed `lntr` form is used in `ini`, since the count is strictly positive.
- **Source aliases:**
  - `NN` -- NONMEM control-stream notation for the transit-chain shape parameter (Abdelgawad 2024 `$THETA 6 NN`, Svensson 2018 rifampicin `THETA(9)`).
  - `N`, `n_transit`, `NTR` -- equivalent paper notation.
- **Example models:** `Abdelgawad_2024_linezolid.R` (founding example -- `lntr = log(5.68)`, an estimated chain length feeding rxode2's analytical `transit`).
- **Notes:** `ntr` is the canonical name for new extractions. Do NOT use `lnn` / `nn` for a transit-chain length: the `lnn` / `nn_fix` form is reserved for the Wang and co-authors sigmoidicity exponent in specific BDE / morphine-like models, so reusing it here would be a role collision. Legacy files predating this entry (`Svensson_2018_rifampicin.R`, `Smythe_2013_gatifloxacin.R`) still carry `lnn` / `nn` and are not retrofitted. Distinct from `ktr`, which is the per-compartment rate, and from `mtt`, which is the mean time through the whole chain.

### kmet (**canonical metabolite-formation rate constant**)
- **Type:** paper-named-param
- **Role:** First-order rate constant for in-vivo formation of an active metabolite from the parent central compartment, used when the source paper parameterises metabolite formation independently of the parent's total clearance.
- **Source aliases:** none.
- **Example models:** `Krause_2017_selexipag.R` (parent -> ACT-333679 active metabolite).
- **Notes:** Fraction-metabolised is implicit (`kmet * Vc_parent / CL_parent`) rather than estimated as `fm`. Distinct from `k_presystemic_central_sa` and its siblings, which describe PRE-systemic (first-pass) metabolite formation out of the `presystemic` compartment rather than systemic conversion out of `central`. Also has a log-transformed primary form `lkmet` (`Krause_2017_selexipag.R`, `Koh_2025_aspirin.R`).

### k_presystemic_central (**canonical pre-systemic-to-parent-central transfer rate constant**)
- **Type:** paper-named-param
- **Role:** First-order rate constant carrying INTACT parent drug out of the `presystemic` compartment into the parent's `central` compartment, in oral models that resolve first-pass metabolism as an explicit structural step. Paired with `k_presystemic_central_<metab>`, the competing route that carries pre-systemically formed metabolite straight into the metabolite's central compartment; the two rate constants together determine the pre-systemic extraction ratio. Log-transformed primary form `lk_presystemic_central`.
- **Source aliases:**
  - `k23` -- Koh 2025 Table 1 / Appendix 1 micro-constant notation (pre-systemic compartment 2 to ASA central compartment 3).
- **Example models:** `Koh_2025_aspirin.R` (doi:10.2147/DDDT.S533428; k23 = 2.32 1/h, RSE 4.11%).
- **Notes:** Named by source-and-destination role rather than by the paper's `k23` / `k24` digits: the registered micro-constant names `k12` / `k21` / `k13` / `k31` are bound to the standard central/peripheral topology, whose numbering is fixed, whereas this paper's numbering indexes its own non-standard pre-systemic compartment, so the digits carry no transferable meaning and would silently break if a state were inserted.

### k_presystemic_central_sa (**canonical pre-systemic metabolite-formation rate constant (salicylic acid)**)
- **Type:** paper-named-param
- **Role:** First-order rate constant carrying pre-systemically formed metabolite out of the `presystemic` compartment directly into the metabolite's `central_<metab>` compartment -- i.e. the drug that is hydrolysed / metabolised in the gut wall before ever reaching the systemic circulation. The `_sa` form is the salicylic-acid member; the general shape is `k_presystemic_central_<metab>` for any registered metabolite suffix. Log-transformed primary form `lk_presystemic_central_sa`.
- **Source aliases:**
  - `k24` -- Koh 2025 Table 1 / Appendix 1 micro-constant notation (pre-systemic compartment 2 to SA central compartment 4).
- **Example models:** `Koh_2025_aspirin.R` (doi:10.2147/DDDT.S533428; k24 = 0.57 1/h, fixed from Dings & Lehr per Koh 2025 reference 14 to aid convergence).
- **Notes:** Distinct from `kmet`, which is SYSTEMIC parent-to-metabolite conversion out of `central`; a model carrying both is precisely how the pre-systemic and systemic contributions are separated.

### fm (**canonical fraction metabolised**)
- **Type:** paper-named-param
- **Role:** Fraction of parent clearance routed to an active metabolite, estimated as a structural parameter in parent-metabolite joint popPK models. Unitless, bounded in (0, 1].
- **Source aliases:** none.
- **Example models:** `Danielak_2017_clopidogrel.R` (clopidogrel -> H4 active thiol, doi:10.1007/s00228-017-2334-z).
- **Notes:** Distinct from `kmet` (formation rate constant): `fm` is unitless and bounded; `kmet` has rate units.

### fm_125d3, fm_1ohm, fm_25d3, fm_aca, fm_cyp3a4, fm_cyp3a5, fm_gx, fm_h4, fm_ko516_frac, fm_m2, fm_m27, fm_m3034, fm_m3g, fm_m5, fm_m6g, fm_megx, fm_morphine, fm_morphine_raw, fm_other (**canonical fraction metabolised by a named pathway**)
- **Type:** paper-named-param
- **Role:** Fraction of a parent drug's clearance routed to one named elimination pathway, used when a paper splits the metabolic flux across two or more identified routes instead of estimating the single `fm` above. Unitless, bounded in (0, 1]; the members of one model's split are expected to sum to 1, either explicitly or with the residual route written inline as `1 - fm_a - fm_b`. The family shape is `fm_<pathway>`, where `<pathway>` is a lowercase token naming the route -- an enzyme (`fm_cyp3a4`, `fm_cyp3a5`) or the metabolite that route forms, reusing the metabolite suffix registered in `compartment-names.md` (`fm_h4`, `fm_m3g`, `fm_m6g`).
- **Source aliases:**
  - `fmCYP3A4` -- Pei 2023 Table 1 notation for the CYP3A4 fraction (`fmCYP3A` is the Discussion's name for the CYP3A4 + CYP3A5 sum, from which `fm_other` is the complement).
  - `fmCYP3A5` -- Pei 2023 Table 1 notation for the CYP3A5 fraction.
  - `fmH4` -- Jung 2024 Table 3 notation for the clopidogrel H4 fraction.
  - `fmcarbo` -- Jung 2024 Table 3 notation for the clopidogrel carboxylic acid fraction; see the `cloca` metabolite suffix in `compartment-names.md`.
  - `fmothers` -- Jung 2024 Table 3 notation for the residual route, and the reason the canonical is fixed as the singular `fm_other` below.
  - `fm CYP3A4`, `fm,CYP3A4`, `fm3A4` -- the usual CYP3A4 table-row spellings.
- **Example models:** `Pei_2023_tacrolimus_pbpk.R` (`fm_cyp3a5` 0.55, `fm_cyp3a4` 0.35, `fm_other` 0.10, summing to 1 over the hepatic routes), `Franken_2015_morphine.R` (`fm_m3g` 0.55 and `fm_m6g` 0.10, with the residual route formed inline as `cl_other <- (1 - fm_m3g - fm_m6g) * cl`), `Pejcic_2024_clopidogrel.R` (`fm_h4` 0.12), `Mitra_2026_ziftomenib.R` (`fm_ko516_frac`; see the `_frac` note below), `Jaiswal_2025_dordaviprone.R` (`fm_cyp3a4` fixed at 0.8, the remaining 0.2 being CYP2D6).
- **Notes:** A new pathway may be added to the heading without a fresh naming sidecar so long as `<pathway>` is lowercase and matches the enzyme or the canonical metabolite suffix it names. Three conventions the family fixes, each of which had already drifted:
  - **The residual route is `fm_other`, singular.** It carries the fraction eliminated by all remaining or unspecified pathways. This matches the singular `_other` member registered for the `lkp_<tissue>` / `kp_<tissue>` families (`lkp_other`, `kp_other`), and the de-facto `cl_other` used in `Smythe_2013_gatifloxacin.R` and `Franken_2015_morphine.R`. A plural `fm_others` was in the tree until 2026-08-20 and is not a permitted alias.
  - **The pathway token is lowercase**, even when the source paper capitalises the metabolite (`fm_h4`, not `fm_H4`). Both spellings of the clopidogrel H4 pathway were in the tree until 2026-08-20.
  - **Positional splits keep the paper's own names.** A paper that parameterises a **nested** split positionally -- `Jung_2024_clopidogrel.R`, whose `fm1` and `fm2` are carried on the logit scale as `logitfm1` / `logitfm2` -- keeps `fm1` / `fm2` rather than being recoded into `fm_<pathway>` members, because `fm1` is a share of the flux that `fm2` already selected, not a share of total clearance, so a `fm_<pathway>` name would misstate the denominator. Those names do not carry the `fm_` prefix and so are not gated.
  - **A trailing `_frac` marks a nested share**, i.e. a fraction *of `fm`* rather than of total clearance. `Mitra_2026_ziftomenib.R` splits its already-estimated `fm` between two metabolites as `k23 <- (cl / vc) * fm * (1 - fm_ko516_frac)` and `k24 <- (cl / vc) * fm * fm_ko516_frac`, so the fraction of *parent clearance* reaching KO-516 is `fm * fm_ko516_frac`, and a bare `fm_ko516` would name the wrong quantity. Use `_frac` only for this nested case; a plain `fm_<pathway>` is always a share of total clearance.

Two classes of name collide on the `fm_` prefix but are **not** members of this family and must not be read as fractions metabolised. Both are outside the gate's scope by construction, because it inspects only names bound inside `ini` and `model`:
  - `fm` also abbreviates **fat mass**, the canonical covariate `FM` in `covariate-columns.md`. `Robarge_2017_efavirenz.R` carries `fm_range = "6.3-42.7 kg (median 18.8)"` in its `population` metadata beside `ffm_range`, and its vignette's simulation chunk derives a local `fm_kg` as total body weight minus fat-free mass. Both are body composition in kg, not pathway fractions, and neither is a model parameter.
  - Source-paper notation quoted inside `label` strings and comments is prose, not a parameter name. `Svensson_2013_bedaquiline.R` labels its apparent metabolite parameters `"...CL_M2/(F*fm_M2)..."`, reproducing the paper's own capitalisation consistently with the `CL_M2` / `V_M2` / `VP1_M2` symbols in the same strings; leave such quotations as the source wrote them.

  - **A trailing `_raw` marks a pre-rescaling value.** `Ashraf_2024_codeine.R` estimates `fm_morphine_raw` as the paper's printed `f_morphine` and then applies the source's own `f/(1+f)` rescaling to obtain `fm_morphine`, so both denote a fraction metabolised to morphine but only the rescaled one is a share of total clearance. Use `_raw` only where the paper itself reports an unrescaled value alongside the rescaled one.

  Two members were added 2026-09-02 alongside the Keunecke 2020 regorafenib extraction: `fm_m2` (the pre-systemically formed fraction of the dose appearing as regorafenib N-oxide) and `fm_m5` (the fraction of parent clearance forming regorafenib N-oxide N-desmethyl). Both pathway tokens are registered metabolite suffixes in `compartment-names.md`, which is the condition this heading sets for adding one without a naming sidecar. Keunecke 2020 holds both on the **logit** scale, so `ini()` carries `logitfm_m2` / `logitfm_m5` and `model()` derives the bare `fm_m2` / `fm_m5` -- see the `logitfm` entry below for that metabolite-suffixed form.

Ten members were added 2026-08-23 during the round-3 consolidation, the first run of this gate over the whole modeldb: `fm_morphine` / `fm_morphine_raw` (Ashraf 2024 codeine), `fm_aca` (Demin 2025), `fm_megx` / `fm_gx` (He 2025 lidocaine), `fm_25d3` / `fm_125d3` (Tuey 2024 cholecalciferol), `fm_m3034` / `fm_m27` (Willmann 2024 elinzanetant) and `fm_1ohm` (Xie 2025 midazolam). Every pathway token was already a registered metabolite suffix in `compartment-names.md`, which is the condition this heading sets for adding one without a naming sidecar. An eleventh name the gate caught, `fm_mass` in `Xie_2025_midazolam.R`, was NOT registered: it is `fm_1ohm` converted to a mass basis rather than a share of clearance, so it was renamed `f_mass_1ohm` and put outside the family.

Enforced mechanically by `.checkFmFamily` in `R/checkModelConventions.R`, which raises an error-severity `fm_family` issue for any `fm_`-prefixed name bound in a model's `ini` or `model` block that is not registered in this heading, so `buildModelDb` fails rather than shipping a new spelling. Adding a pathway means adding it to the heading above -- that is the enumerating step.

### logitfm (**canonical logit-transformed fraction metabolised**)
- **Type:** paper-named-param
- **Role:** Logit-transformed fraction of parent clearance routed to an active metabolite. Used when the source paper's estimation routine holds FM on the logit scale so that FM is bounded in (0, 1) regardless of covariate + eta combinations. Inside `model` the bare form is `fm = 1 / (1 + exp(-logitfm_ind))` where `logitfm_ind` collects the fixed effect, covariate shifts, and IIV on the logit scale.
- **Source aliases:** `FRM2` / `FRM5` -- Keunecke 2020 notation for the two regorafenib metabolite formation fractions, both reported on the logit scale.
- **Example models:** `Mitra_2026_ziftomenib.R` (base `logitfm <- fixed(0.14)` corresponding to `logit^-1(0.14) = 0.535`; additive shift `e_dis_healthy_logitfm = -1.62` on the logit scale for the healthy-volunteer cohort; IIV `etalogitfm ~ 0.280` on the logit scale), `Keunecke_2020_regorafenib_phase1.R` / `Keunecke_2020_regorafenib_phase3.R` (metabolite-suffixed `logitfm_m2 = -0.355` and `logitfm_m5 = -1.09`, whose inverse-logits reproduce the paper's stated 41% and 25%).
- **Metabolite-suffixed form:** when a paper splits the flux across named routes AND holds each on the logit scale, take the metabolite suffix on the `logitfm` root -- `logitfm_<metab>`, with IIV `etalogitfm_<metab>` and the individual logit-scale value collected as `logitfm_<metab>_ind` before the inverse-logit. The bare `fm_<metab>` derived inside `model()` must be registered in the `fm_<pathway>` heading above, because that is the name the `.checkFmFamily()` gate inspects.
- **Notes:** Follows the `logit`-transform-prefix family (`logitffo`, `logitemax`) applied to the existing `fm` canonical root. Keep the `_ind` intermediate on its own line: folding the eta into the inverse-logit expression (`1 / (1 + exp(-(logitfm + etalogitfm)))`) breaks rxode2's mu-referencing detection and raises "some etas defaulted to non-mu referenced". Prefer `logitfm` (paper's logit encoding) over `lfm` (log encoding) when the source paper's NONMEM control stream stores FM on the logit scale, because a log-scale encoding of a bounded quantity can leak above 1 under moderate eta or covariate values.


### lq_milk (**canonical log-transformed central-to-breast-milk inter-compartmental clearance**)
- **Type:** log-transformed-pk
- **Role:** Apparent inter-compartmental clearance between a central compartment and its paired `milk` compartment in a lactation-transfer model (volume / time). Extends the existing `lq` family with a destination token, on the same principle as the `kin_<compartment>` / `kout_<compartment>` family: the canonical bare `lq` / `lq2` mean exchange with `peripheral1` / `peripheral2` specifically, so reusing them for a breast-milk compartment would silently mislead. The bare form `q_milk` is used inside `model`.
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine.R`, `Wattanakul_2024_primaquine_motherinfant.R` (`Q/F = 0.400 L/h`, shared between primaquine and carboxyprimaquine per Wattanakul 2024 Table 2 footnote c; the source paper assumed and encoded `Q/F_CPQ = Q/F_PQ`).
- **Notes:** Metabolite forms take the standard suffix (`lq_milk_<metab>`) when a paper estimates the parent and metabolite transfer clearances separately.

### lq_elf (**canonical log-transformed central-to-lung-ELF inter-compartmental clearance**)
- **Type:** log-transformed-pk
- **Role:** Inter-compartmental clearance between a central compartment and the lung epithelial-lining-fluid space (volume / time). Second member of the `lq_<destination>` family founded by `lq_milk`, for the same reason: the bare `lq` / `lq2` mean exchange with `peripheral1` / `peripheral2` specifically, so a lung-ELF exchange needs its own destination token. The bare form `q_elf` is used inside `model`. Note that the exchange is generally NOT symmetric in rate -- `q_elf` is a clearance applied to whichever compartment the drug is leaving, so the ELF-to-central micro-constant is `q_elf / v_elf` while the central-to-ELF one is `q_elf / vc`.
- **Source aliases:** `Q1` -- Kurup 2024 Table 1 ("Intercompartmental clearance ELF - central") and supplement eqs. (1)-(4).
- **Example models:** `Kurup_2024_DZIF10c.R` (`Q1 = 3.68e-5 L/h` at WT = 70 kg, allometrically scaled with a fixed exponent of 0.85).
- **Notes:** Applies whether the ELF space is carried as the bare `elf` compartment alone or as `elf` plus regional `elf_<region>` pools that share one volume; in the founding example one clearance serves all three pools.

### lcl_elf (**canonical log-transformed lung-ELF loss clearance**)
- **Type:** log-transformed-pk
- **Role:** First-order clearance of drug out of the lung epithelial-lining-fluid space by a route that is not distribution into the systemic circulation (volume / time) -- mucociliary escalator transport, local catabolism, or swallowing. It is a clearance FROM the `elf` compartment, so it is deliberately NOT one of the `lcl_<arm>` names (`lcl_renal`, `lcl_nonren`, `lcl_hemodialysis`,...), all of which are additive arms of the CENTRAL clearance. Reading it as an arm of `cl` would misplace the sink by a whole compartment. The bare form `cl_elf` is used inside `model()`, and the loss micro-constant is `cl_elf / v_elf`.
- **Source aliases:** `CL_L`, `CLL` -- Kurup 2024 Table 1 ("Clearance from ELF (IT/INH)") and supplement eqs. (1)-(2).
- **Example models:** `Kurup_2024_DZIF10c.R` (`CL_L = 0.000412 L/h` at WT = 70 kg; the paper attributes it to mucociliary clearance and it is what gives an inhaled dose its ~2.5-day effective lung half-life against the ~20-day systemic half-life).
- **Notes:** A model that resolves several distinct non-systemic lung sinks should suffix them further rather than overload this name.

### lv_elf (**canonical log-transformed lung-ELF volume**)
- **Type:** log-transformed-pk
- **Role:** Volume of the lung epithelial-lining-fluid space (volume). Needed because the bare volume canonicals `lvc` / `lvp` / `lvp2` / `lvp3` are reserved for the central and numbered peripheral compartments of a classical-PK model, and an ELF volume is a named physiological space rather than a peripheral compartment. The bare form `v_elf` inside `model` predates this entry -- it is already used as a hardcoded physiological constant in `Parmar_2023_spectinamide_1599_mouse_pbpk.R` -- so this registers the `l`-prefixed `ini` form for models that estimate it.
- **Source aliases:** `V_L`, `VL` -- Kurup 2024 Table 1 ("ELF compartment volume") and supplement eqs. (1)-(3), (6).
- **Example models:** `Kurup_2024_DZIF10c.R` (`V_L = 0.0364 L` at WT = 70 kg, estimated; allometrically scaled with a fixed exponent of 1). Bare-form precedent: `Parmar_2023_spectinamide_1599_mouse_pbpk.R` (`v_elf = 1.00e-5 L`, a fixed mouse physiological constant).
- **Notes:** When a model carries regional `elf_<region>` pools alongside the bare `elf`, a single `v_elf` normally serves all of them and the whole-lung ELF concentration is the summed amount over that one volume; a model that gives each region its own volume should register `lv_elf_<region>` names.

### lq_elf_lung (**canonical log-transformed ELF-to-lung inter-compartmental clearance**)
- **Type:** log-transformed-pk
- **Role:** Inter-compartmental clearance between the lung epithelial-lining-fluid space and a second, deeper lung distribution compartment (volume / time). Third member of the `lq_<destination>` family founded by `lq_milk` and continued by `lq_elf`, and shaped like the `lq_p1_p2` two-token form because neither end of the exchange is the central compartment: the bare `lq` / `lq2` mean exchange between `central` and `peripheral1` / `peripheral2` specifically, so an exchange whose SOURCE is `elf` needs both endpoints in the name. Unlike `lq_p1_p2` the flow is bidirectional. The bare form `q_elf_lung` is used inside `model()`; the ELF-to-lung micro-constant is `q_elf_lung / v_elf` and the lung-to-ELF one is `q_elf_lung / v_lung`.
- **Source aliases:** `Q2` -- Saporta 2026 Table 1 ("Apparent intercompartmental clearance (ELF-lung)") and Figure 2.
- **Example models:** `Saporta_2026_meropenem.R` (`Q2 = 2.48 L/(h.kg)`, RSE 13%, in a mouse lung infection model whose ELF limb is a two-compartment system driven by, but not depleting, the plasma compartment).
- **Notes:** Auto-approved 2026-09-01 as a well-formed member of the documented `lq_<destination>` family alongside the Saporta 2026 meropenem extraction (operator sidecar `oare_PMC13041408` request-001, "COMPARTMENT NAMES NEED NO SIDECAR" paragraph). Pairs with `lv_lung`. Distinct from `lq_elf`, which is the CENTRAL-to-ELF clearance; a model resolving both, as Saporta 2026 does, uses the two names together.

### lv_lung (**canonical log-transformed lung-compartment volume**)
- **Type:** log-transformed-pk
- **Role:** Apparent volume of the canonical `lung` compartment (volume). Needed for the same reason as `lv_elf`: the bare volume canonicals `lvc` / `lvp` / `lvp2` / `lvp3` are reserved for the central and numbered peripheral compartments of a classical-PK model, and a lung distribution space is a named physiological compartment rather than a peripheral one. The bare form `v_lung` is used inside `model()` to convert the lung amount to a concentration, `Clung <- lung / v_lung`.
- **Source aliases:** `VL2`, `V L2` -- Saporta 2026 Table 1 ("Apparent volume of distribution of second lung compartment").
- **Example models:** `Saporta_2026_meropenem.R` (`VL2 = 9.27 L/kg`, RSE 26%; an apparent volume, because the lung limb is driven by plasma without mass transfer so the whole limb absorbs the unknown partitioning).
- **Notes:** Auto-approved 2026-09-01 as a well-formed member of the documented `lv_<compartment>` family founded by `lv_elf`, alongside the Saporta 2026 meropenem extraction (operator sidecar `oare_PMC13041408` request-001, "COMPARTMENT NAMES NEED NO SIDECAR" paragraph). Pairs with `lq_elf_lung`. Distinct from the `kp_lung` / `lkp_lung` tissue-to-plasma partition coefficients, which are whole-body-PBPK partitioning constants rather than a fitted compartment volume.

### lq_kidney (**canonical log-transformed central-to-kidney inter-compartmental flow**)
- **Type:** log-transformed-pk
- **Role:** Blood (or plasma) flow perfusing a `kidney` compartment, applied as a bidirectional inter-compartmental clearance between `central` and `kidney` (volume / time). Third member of the `lq_<destination>` family founded by `lq_milk` and continued by `lq_elf`, for the same reason those exist: the bare `lq` / `lq2` mean exchange with `peripheral1` / `peripheral2` specifically, so a perfused organ compartment needs its own destination token. The bare form `q_kidney` is used inside `model()`. In a flow-limited organ the same `q_kidney` appears with opposite sign on both sides (`q_kidney * Cc` leaving central, `q_kidney * Ck` returning), so the organ is effectively well-stirred with a partition coefficient of 1 unless the model states otherwise.
- **Source aliases:** `QR` -- Felmlee 2010 Table I ("Renal blood flow") and eqs. (1), (5).
- **Example models:** `Felmlee_2010_ghb_rat.R` (`QR = 12.5 mL/min`, fixed to the rat physiological value; Table I footnote a).
- **Notes:** Ratified 2026-09-02 (operator sidecar `oare_PMC2895455` request-001 / response-001, question q1, option A). Distinct from every `lcl_<arm>` name (`lcl_renal`, `lcl_nonren`, ...), which are additive arms of the CENTRAL clearance rather than a flow into a separate organ compartment; reading a perfusion flow as a clearance arm would delete the organ compartment entirely. Other perfused organs should take the same `lq_<organ>` shape.

### lv_kidney (**canonical log-transformed kidney compartment volume**)
- **Type:** log-transformed-pk
- **Role:** Volume of a `kidney` compartment in a hybrid-physiological (semi-PBPK) model (volume). Second member of the `lv_<space>` family founded by `lv_elf`, for the reason given there: the bare volume canonicals `lvc` / `lvp` / `lvp2` / `lvp3` are reserved for the central and numbered peripheral compartments of a classical-PK model, and a named organ is not a numbered peripheral compartment. The bare form `v_kidney` is used inside `model()`, and the organ concentration is `kidney / v_kidney`.
- **Source aliases:** `Vkidney` -- Felmlee 2010 Table I ("Kidney volume") and eqs. (5), (6).
- **Example models:** `Felmlee_2010_ghb_rat.R` (`Vkidney = 4.0 mL`, fixed to the rat physiological value; Table I footnote a).
- **Notes:** Ratified 2026-09-02 (operator sidecar `oare_PMC2895455` request-001 / response-001, question q1, option A). Paired with `lq_kidney` in the founding example but independent of it -- a model may carry a kidney volume without an explicit perfusion flow, or vice versa.

### lgfr, gfr (**canonical glomerular filtration rate**)
- **Type:** log-transformed-pk
- **Role:** Glomerular filtration rate as a physiologic volumetric flow OUT of a kidney compartment and INTO a renal ultrafiltrate / tubular compartment (volume / time). It is deliberately NOT `cl_renal`: `lcl_renal` is an additive arm of the CENTRAL clearance that lumps filtration, secretion and reabsorption into one net number applied to plasma, whereas `gfr` is the filtration step alone, applied to the kidney compartment concentration, with reabsorption represented explicitly as a separate returning flux. The registered `lcl_tsnet` entry already presumes this parameter exists ("the filtration arm is carried by the model's glomerular-filtration parameter"). The bare form `gfr` is used inside `model()`.
- **Source aliases:** `GFR` -- Felmlee 2010 Table I ("Glomerular filtration rate") and eqs. (5), (6).
- **Example models:** `Felmlee_2010_ghb_rat.R` (`GFR = 10 mL/min/kg`, entering as `gfr <- exp(lgfr) * WT`; the only weight-scaled term in that model).
- **Notes:** Ratified 2026-09-02 (operator sidecar `oare_PMC2895455` request-001 / response-001, question q1, option A). Physiology tables report GFR per unit body weight, so a source value is usually per-kg even when every other flow in the same table is an absolute per-animal value -- check the units column before wiring it in, and scale by `WT` inside `model()` rather than baking a body weight into the `ini()` value.

### luf, uf (**canonical urine flow through the renal tubule**)
- **Type:** log-transformed-pk
- **Role:** Volumetric flow of ultrafiltrate along the renal tubule (volume / time), carrying drug by convection from one ultrafiltrate compartment to the next and finally into a cumulative `urine` compartment. The bare form `uf` is used inside `model()`.
- **Source aliases:** `UF` -- Felmlee 2010 Table I ("Urine flow") and eqs. (6), (7), (8).
- **Example models:** `Felmlee_2010_ghb_rat.R` (`UF = 0.1 mL/min`, fixed to the rat physiological value with an estimated 114% CV between-subject variability; Table I footnote a).
- **Notes:** Ratified 2026-09-02 (operator sidecar `oare_PMC2895455` request-001 / response-001, question q1, option A, in preference to widening the registered `lurprod` / `urprod`). Kept distinct from `urprod` deliberately: `urprod` is scoped to an ESTIMATED zero-order rate of urine production that integrates the canonical `urine_vol` volume state, whereas `uf` is a FIXED physiological flow that moves drug mass between tubular compartments and never integrates a volume. The two are the same physical dimension and could be numerically equal in a given paper, but they occupy different structural roles, and the register's standing warning against giving one name two meanings applies.

### lvmax_reab, vmax_reab (**canonical maximum rate of saturable renal tubular reabsorption**)
- **Type:** log-transformed-pk
- **Role:** Michaelis-Menten maximum rate for active (transporter-mediated) reabsorption of drug out of a renal ultrafiltrate compartment and back into the kidney (mass / time). Paired with `lkm_reab`; both are meaningless alone. Follows the documented `vmax_<process>` / `km_<process>` disambiguation pattern (`vmax_pah`, `km_pah`, `vmax_trans`, `km_trans`, and the registered `lvmax_rbc` / `lkm_rbc` pair), and is required whenever the bare `lvmax` / `lkm` are already taken by a metabolic arm in the same model.
- **Source aliases:** `Vmax,R` -- Felmlee 2010 Table I ("Maximum renal reabsorption rate") and eqs. (4), (11)-(13).
- **Example models:** `Felmlee_2010_ghb_rat.R` (`Vmax,R = 2.34 mg/min`, estimated; the reabsorption flux is `vmax_reab * Culf1 / (km_reab + Culf1)`).
- **Notes:** Ratified 2026-09-02 (operator sidecar `oare_PMC2895455` request-001 / response-001, question q1, option A). Note the encoding convention in the founding example: the source paper defines a quantity it calls "Reabsorption" that is `Vmax,R / (Km,R + Culf1)` -- a CLEARANCE -- and then multiplies it by `Culf1` in the differential equations, so the mass flux is the usual Michaelis-Menten form.

### lkm_reab, km_reab (**canonical Michaelis constant for saturable renal tubular reabsorption**)
- **Type:** log-transformed-pk
- **Role:** Michaelis-Menten constant (concentration at half-maximal rate) for active reabsorption of drug out of a renal ultrafiltrate compartment, expressed as a concentration in that compartment. Paired with `lvmax_reab`; both are meaningless alone.
- **Source aliases:** `Km,R` -- Felmlee 2010 Table I ("Renal reabsorption Michaelis-Menten constant") and eqs. (4), (11)-(13).
- **Example models:** `Felmlee_2010_ghb_rat.R` (`Km,R = 0.46 mg/mL`, estimated; the paper notes it agrees with the in-vitro MCT1 uptake `Km` of 0.48 mg/mL).
- **Notes:** Ratified 2026-09-02 (operator sidecar `oare_PMC2895455` request-001 / response-001, question q1, option A). A transport `Km` for reabsorption is generally far away from a metabolic `km` for the same drug -- in the founding example the reabsorption `km_reab` of 0.46 mg/mL sits against a metabolic `km` of 0.054 mg/mL in the same model, an eightfold separation -- so never fold them into one name. The reabsorption value is the transporter-affinity one: Felmlee 2010 Discussion places it inside the in-vitro MCT-isoform range of 0.48-1.9 mg/mL and notes that SMCT1, which also carries GHB in the proximal tubule, has a higher affinity still (0.17 and 0.07 mg/mL), close enough that the two transport processes are not separable in vivo.

### logitfdepot (**canonical logit-transformed depot fraction**)
- **Type:** log-transformed-pk
- **Role:** Logit-scale encoding of the fraction of a dose that reaches the depot compartment -- the same quantity as `lfdepot`, held on the logit scale so that it is bounded in (0, 1) regardless of the covariate and eta combination. Inside `model` the bare form is `f = exp(phi) / (1 + exp(phi))` where `phi` collects the fixed effect, covariate shifts, and IIV on the logit scale; IIV is therefore additive on the logit scale, not multiplicative.
- **Source aliases:** `F1` (estimated as `THETA` with `F1 = EXP(PHI)/(1+EXP(PHI))`) -- Kurup 2024 Table 1 and supplement eqs. (8)-(9), (13).
- **Example models:** `Kurup_2024_DZIF10c.R` (inhaled lung deposition, human `F1h = 29.6%`; the species-specific sibling is `logitfdepot_macaque` at `F1m = 1.00%`).
- **Notes:** Follows the `logit`-transform-prefix family (`logitffo`, `logitemax`, `logitfm`) applied to the existing `fdepot` root. Prefer `logitfdepot` over `lfdepot` when the source paper's control stream estimates the fraction on the logit scale, for the reason given in the `logitfm` entry: a log-scale encoding of a bounded quantity can leak above 1 under moderate eta or covariate values. Species-specific and stratum-specific values take the usual suffix (`logitfdepot_macaque`, per the `lfdepot_macaque` precedent in `Nagy_2017_obiltoxaximab.R`), with the reference species keeping the bare name.

### pcmilk (**canonical milk:plasma partition coefficient**)
- **Type:** paper-named-param
- **Role:** Milk:plasma partition coefficient -- the fraction of the drug in the central compartment that is freely distributed into breast milk. Unitless, normally in (0, 1]. Enters the central-to-milk micro-rate constant as `k = (q_milk / vc) * pcmilk`, so at pseudo-equilibrium the milk:plasma concentration ratio (and hence the milk:plasma AUC ratio) equals `pcmilk` exactly. That identity is the natural falsifier for a lactation extraction: a published milk:plasma AUC ratio that does not equal the published partition coefficient means the model that produced it had an additional milk sink. Metabolite forms take the standard suffix (`pcmilk_<metab>`).
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine.R`, `Wattanakul_2024_primaquine_motherinfant.R` (`PC_PQ = 0.376`, `PC_CPQ = 0.00889` per Wattanakul 2024 Table 2; Table 3 reports the simulated milk:plasma AUC ratio as 0.376 (0.375-0.377), equal to `pcmilk`).
- **Notes:** Distinct from a "milk-to-plasma ratio" reported as a raw observed concentration ratio: `pcmilk` is a model parameter that generates such a ratio, not the ratio itself.

### cfcap (**canonical capillary:venous conversion factor**)
- **Type:** paper-named-param
- **Role:** Proportional conversion factor between a venous-plasma model prediction and the corresponding capillary (finger-prick) plasma prediction, estimated when a study assays both matrices and fits them simultaneously. Unitless; the capillary observable is `Ccap = Cc * cfcap`. Metabolite forms take the standard suffix (`cfcap_<metab>`).
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine.R`, `Wattanakul_2024_primaquine_motherinfant.R` (`CF_PQ = 0.898`, `CF_CPQ = 1.06` per Wattanakul 2024 Table 2; the associated observation variables are `Ccap` and `Ccap_cpq`).
- **Notes:** Capillary sampling is standard in field malaria and paediatric pharmacokinetics, so this is expected to recur. Do NOT reuse `cfcap` for a plasma:whole-blood or plasma:serum conversion; those are different matrices and should get their own canonical.

### kmilkinf (**canonical breast-milk-to-infant transfer rate constant**)
- **Type:** paper-named-param
- **Role:** First-order rate constant (1/time) transferring drug from a mother's `milk` compartment into the breastfed infant's `infant_depot` compartment during a feeding window in a mother-to-infant dyad model. Normally fixed to a large value so that essentially all of the milk-compartment amount transfers within the window rather than estimated, because the feed is fast relative to every other process in the model.
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine_motherinfant.R` (`MTINF_PQ = MTINF_CPQ = 100 /h`, fixed high so that more than 95% of the amount in the breast-milk compartment is transferred during the 24-minute feeding window per Wattanakul 2024 Results, 'Predicting infant concentrations').

### feed_n, feed_window, feed_first, milk_intake (**canonical breastfeeding-pattern constants**)
- **Type:** paper-named-param
- **Role:** The four design constants that define a breastfeeding schedule in a lactation-transfer model. `feed_n` is the number of feeds per day (giving a feeding cycle of `24 / feed_n` hours); `feed_window` is the duration of each feed (time); `feed_first` is the time of the first feed after time zero (time); `milk_intake` is the average daily breast-milk intake per kg of infant body weight (volume / mass / time), from which the milk-compartment volume is derived as `milk_intake * WT_INFANT / feed_n`. All four are normally `fixed` -- they describe the simulated feeding scenario, not an estimated property of the drug -- and are the natural handles for a feeding-pattern sensitivity analysis.
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine.R`, `Wattanakul_2024_primaquine_motherinfant.R` (`feed_n = 10` feeds/day, `feed_window = 0.4` h, `feed_first = 1` h, `milk_intake = 0.15` L/kg/day, giving a 2.4 h cycle with a 24-minute feeding window; Wattanakul 2024 Eq. 1 and Results, 'Predicting infant concentrations').
- **Notes:** Kept as separate scalars rather than folded into one schedule object so that each can be varied independently in a sensitivity analysis, which is exactly what the founding paper does (Supplementary Figs. 2-4 vary feeding frequency and volume per feed).

### lmtt_infant (**canonical log-transformed breastfed-infant mean transit time**)
- **Type:** log-transformed-pk
- **Role:** Mean transit absorption time of the breastfed infant in a mother-to-infant dyad model, when the infant's absorption is not scaled from the mother's but taken from a separate paediatric source. Uses the `_infant` suffix in the same sense as the `infant_<canonical>` compartment namespace and the `WT_INFANT` / `AGE_INFANT` covariates: it marks the quantity as belonging to the dyad PARTNER rather than to the modelled subject, whose own mean transit time is the plain `lmtt`.
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine_motherinfant.R` (`MTT_infant = 0.706 h`, fixed from a paediatric literature model to account for more rapid absorption in infants, with 2 transit compartments so `ktr_infant = 3 / MTT_infant`; Wattanakul 2024 Methods, 'Predicting infant concentrations').
- **Notes:** Any other canonical parameter that a dyad model must carry separately for the infant partner takes the same `_infant` suffix; parameters that are merely SCALED from the mother's (allometric clearances and volumes) do not get their own `ini()` name at all, because they are derived inside `model()`.

### kin (**canonical indirect-response production rate**)
- **Type:** paper-named-param
- **Role:** Zero-order production rate of an indirect-response / turnover pool (Dayneka 1993; Jusko traditions).
- **Source aliases:** none.
- **Example models:** indirect-response PD models.
- **Notes:** Codified 2026-05-28 per the naming audit. The compartment-suffixed forms `kin_<compartment>` / `kout_<compartment>` (with `lkin_<compartment>` / `lkout_<compartment>` as the log-transformed primaries) are the canonical family for **tissue-exchange rate constants** in semi-physiologic / dosimetry models that parameterise transport as first-order rate constants rather than as clearances: `kin_<tissue>` is central-to-tissue uptake and `kout_<tissue>` is tissue-to-central return, both 1/time. `<compartment>` must be a canonical compartment from `compartment-names.md`. Prefer these role-based names over the source paper's numeric micro-constants (`k12`, `k21`, `k13`, ...) whenever the compartments are anatomically named -- the canonical `k12` / `k21` / `k13` / `k31` entries mean central-to-`peripheral1` / `peripheral2` exchange specifically, so reusing them for organ compartments would silently mislead. Founding examples: `Lindauer_2017_pembrolizumab.R` and `Yamazaki_2008_crizotinib_mouse.R` (`kin_tumor` / `kout_tumor`); `Siebinga_2023_lu177psma617.R` (`kin_salivary_gland`, `kin_kidney`, `kin_liver`, `kin_tumor`, `kin_other` and their `kout_` partners, replacing Siebinga's `k12`/`k21`/`k13`/`k31`/`k14`/`k41`/`k15`/`k51`/`k16`/`k61`); `Nguyen_2026_linezolid.R` (`kin_saliva` / `kout_saliva`, replacing Nguyen's `Kabs` / `Kreabs`, i.e. `K23` / `K32`). The family name carries no commitment to mass balance: Nguyen 2026's saliva state is a DRIVEN (non-depleting) hypothetical effect compartment, so `kin_saliva * central` appears only in the saliva equation and `kout_saliva * saliva` appears in neither the central equation nor anywhere else -- the paired name still records which direction each rate constant describes, which is the point.

### kout (**canonical indirect-response elimination rate**)
- **Type:** paper-named-param
- **Role:** First-order elimination rate of an indirect-response turnover pool (1 / time). Also the tissue-to-central return leg of the `kin_<compartment>` / `kout_<compartment>` tissue-exchange family -- see the `kin` entry above.
- **Source aliases:** none.
- **Example models:** indirect-response PD models; `Lindauer_2017_pembrolizumab.R`, `Siebinga_2023_lu177psma617.R` (tissue-exchange form), `Nguyen_2026_linezolid.R` (`kin_saliva` / `kout_saliva`, the central-to-saliva secretion and saliva-to-central reabsorption legs of a saliva-TDM model).

### lkel_saliva, kel_saliva (**canonical irreversible salivary elimination rate constant**)
- **Type:** paper-named-param
- **Role:** First-order rate constant (1 / time) for IRREVERSIBLE loss of drug out of the `saliva` compartment -- physically swallowing plus salivary flow past the sampling site. It is the third leg of a saliva-TDM compartment, alongside the `kin_saliva` / `kout_saliva` pair from the `kin_<compartment>` / `kout_<compartment>` tissue-exchange family: `d/dt(saliva) = kin_saliva * central - (kout_saliva + kel_saliva) * saliva`. `lkel_saliva` is the log-transformed `ini()` primary and `kel_saliva` the bare form used inside `model()`.
- **Source aliases:**
  - `Kel` -- Nguyen 2026 Table 2 and Figure 1 ("Kel salivary elimination rate"); `K30` in that paper's supplementary Table S7 control stream.
- **Example models:** `Nguyen_2026_linezolid.R` (founding example; `kel_saliva` = 2.13 1/h).
- **Notes:** The compartment suffix is what makes this name necessary rather than a convenience. The source paper's own symbol for the term is "Kel", which is the same word as the canonical bare central elimination rate constant `kel`; writing it unsuffixed in a model that also carries `kel <- cl / vc` would be a silent collision between two different physical processes. Registered as a compartment-suffixed member of the existing `kel_<suffix>` pattern (`kel_renal`, `kel_other`, `kel_gx`, `kel_dihydroart`, drug-suffixed `kel_<drug>`), and it generalises to any second sampling matrix with a one-way loss path -- a future `kel_<compartment>` for another matrix takes the same form. Note that `kout_saliva` and `kel_saliva` are NOT interchangeable even though both appear as first-order losses from the same state: `kout_saliva` is return towards plasma and `kel_saliva` leaves the system, and it is their sum that sets the saliva disposition rate while the RATIO to `kin_saliva` sets the steady-state saliva:plasma exposure ratio. Ratified 2026-09-02 (sidecar request-001 q1) with the Nguyen 2026 extraction; see `saliva` in `compartment-names.md`.

### lsmax, smax (**canonical maximum zero-order production rate of a drug-STIMULATED turnover pool**)
- **Type:** paper-named-param
- **Role:** Saturating maximum of the zero-order PRODUCTION term of an indirect-response pool whose synthesis is stimulated by drug, entering as `d/dt(<pool>) <- smax * Cc / (ec50 + Cc) - kout * <pool>` (or the sigmoid form with `lhill`). Unlike the unitless `Imax`-style fractional multipliers, `smax` carries the pool's own rate units -- `[pool concentration] / time` (IU/mL/h in the founding example) -- because it IS the production rate rather than a scale factor on one. Inside `model()` the bare name is `smax`; the log-transformed `lsmax` form is the `ini()` primary, since a production rate is strictly positive. Pairs with `lec50` (the half-maximal driving concentration) and `lkout` (the pool's first-order elimination rate).
- **Source aliases:**
  - `Smax` -- Benson 2010 Table 2 and equation 1 ("Smax is the maximum rate of synthesis of IFN-alpha (in international units per milliliter per hour)").
  - `SC50` -- the partner half-maximal concentration in the same `Smax` / `SC50` notational convention; maps to the registered `lec50` / `ec50`, NOT to a new name.
- **Example models:** `Benson_2010_bhma_mouse.R` (founding example -- `lsmax <- log(294)` IU/mL/h, the maximum IFN-alpha synthesis rate induced in mice by the TLR-7 agonist BHMA, with `lec50 <- log(135)` ng/mL and `lkout <- log(0.958)` 1/h; the paper's own consistency check is that the maximum attainable response is `smax / kout` = 306 IU/mL).
- **Notes:** Ratified 2026-09-02 (operator sidecar `oare_PMC2825998` request-001 q2, option B). Deliberately NOT `kin`: the registered `kin` is the DRUG-FREE baseline production rate of a turnover pool, and reusing it for a quantity that is zero in the absence of drug would invert its meaning -- in the founding model the baseline production is exactly zero (unstimulated IFN-alpha was below the assay LOQ), so there is no `kin` to speak of. Deliberately NOT `emax` either: the operator's ratification records that `emax` / `lemax` has no entry in this register at all, and the `emax`-shaped names that do appear in model files are unitless fractional multipliers, so overloading that stem with rate units would mislead. The paper-named `smax_<suffix>` forms already in the registry (`Cao_2025_ferricCarboxymaltose_rat.R` `lsmax_iron`, `Garonzik_2016_daptomycin.R` `smax_s`, `Guiastrennec_2016_gastric_emptying.R` `smax_cckf`) are the multi-arm extension of this bare canonical: a model with a single stimulated pool uses `lsmax` / `smax`, and only a model that must carry two or more stimulated production maxima jointly takes the `_<suffix>` form. Where a paper's `Smax` carries `1/time` units because it scales a first-order rather than a zero-order term, it is a rate CONSTANT and belongs under a `kmax`-style name instead -- see `Mody_2023_doxorubicin_dexrazoxane_mdamb468.R`, which maps "Smax,DOX (1/h)" onto `lkmax_dox` for exactly this reason.

### lkinf_rbc, kinf_rbc (**canonical first-order influx rate constant into the red-cell analyte pool**)
- **Type:** paper-named-param
- **Role:** First-order rate constant (1 / time) for LINEAR influx of drug or active metabolite from plasma into the intracellular red-cell pool (`rbc_<analyte>`, see `compartment-names.md`). Enters as `d/dt(rbc_<analyte>) = kinf_rbc * Cc - keff_rbc * rbc_<analyte>`, so the plasma concentration `Cc` must be expressed in the same concentration units as the red-cell state for the term to be dimensionally consistent.
- **Source aliases:**
  - `Kin` -- used in `Gebhard_2023_methotrexate.R` (paper symbol `K_in^MTX`).
- **Example models:** `Gebhard_2023_methotrexate.R` (`K_in^MTX = 0.031 1/day`), `Gebhard_2023_mercaptopurine_anc.R`, `Yu_2026_tenofovir.R` (`K_TFV-TFVdp-DBS = 0.0068 1/h`, driven by plasma tenofovir in a file whose dosed parent is tenofovir alafenamide -- see the Notes).
- **Notes:** The `kinf` stem is deliberately NOT `kin`: the canonical `kin` is a ZERO-ORDER production rate into an indirect-response turnover pool, which is a different quantity with different units, so reusing that stem would actively mislead. Distinct also from `clin` / `clef` (`Campagne_2019_cyclophosphamide_mouse.R`), which are plasma-to-tissue influx / efflux CLEARANCES (volume / time) paired with an `ecf` compartment; `kinf_rbc` / `keff_rbc` are first-order RATE CONSTANTS. The `_rbc` suffix marks the destination pool; no analyte suffix is carried because each model file holds a single drug arm. A future model fitting two red-cell INFLUX ARMS jointly would extend to `lkinf_rbc_<analyte>` -- but only then. Operator sidecar `oare_PMC12783228` (request-001 / response-001, question q2, answered B) settled the adjacent case: a model file that carries two prodrug arms but only ONE red-cell influx keeps the registered bare `lkinf_rbc` / `lkeff_rbc`, with the driving analyte named in the `label()` and the source-trace comment, because forking to a suffixed name away from an entry that already exists is the opposite of reuse-before-register. `Yu_2026_tenofovir.R` is that case: tenofovir alafenamide is the dosed parent and holds the bare compartment names, yet the single red-cell influx is driven by plasma tenofovir and is still `lkinf_rbc`. When a second influx genuinely does appear, add the suffixed sibling ALONGSIDE the bare entry rather than replacing it -- which is exactly how the [[lkinf_pbmc]] / [[lkinf_pbmc_<analyte>]] pair is structured.

### lkeff_rbc, keff_rbc (**canonical first-order efflux / elimination rate constant out of the red-cell analyte pool**)
- **Type:** paper-named-param
- **Role:** First-order rate constant (1 / time) for loss of drug or active metabolite from the intracellular red-cell pool (`rbc_<analyte>`), combining efflux back to plasma, intracellular catabolism, and red-cell turnover into a single lumped rate. Enters as the `- keff_rbc * rbc_<analyte>` term of the red-cell ODE.
- **Source aliases:**
  - `Keff` -- used in `Gebhard_2023_methotrexate.R` and `Gebhard_2023_mercaptopurine.R` (paper symbols `K_eff^MTX`, `K_eff^6MP`).
- **Example models:** `Gebhard_2023_methotrexate.R` (`K_eff^MTX = 0.018 1/day`), `Gebhard_2023_mercaptopurine.R` (`K_eff^6MP = 0.041 1/day`), `Gebhard_2023_mercaptopurine_anc.R` (`K_eff^6MP = 0.050 1/day`).
- **Notes:** Because the rate is lumped, its reciprocal half-life is a red-cell residence property rather than a pure membrane-transport property -- Gebhard 2023's Discussion validates `K_eff^MTX = 0.018 1/day` against literature red-cell methotrexate half-lives of 30-40 days (0.017-0.023 1/day).

### lkinf_pbmc, kinf_pbmc (**canonical first-order influx rate constant into the PBMC analyte pool**)
- **Type:** paper-named-param
- **Role:** First-order rate constant (1 / time) for LINEAR influx of drug or active metabolite from plasma into the intracellular peripheral-blood-mononuclear-cell pool (`pbmc_<analyte>`, see `compartment-names.md`). The PBMC parallel of `lkinf_rbc`: it enters as `d/dt(pbmc_<analyte>) = kinf_pbmc * Cc - keff_pbmc * pbmc_<analyte>`, so the plasma concentration driving it must be expressed in the same concentration units as the PBMC state. The BARE form names the influx driven by the dosed PARENT, per the standing parent-wins naming rule; when a single model file carries two prodrug or two parent arms feeding the same PBMC pool, the second influx takes the source-analyte suffix `lkinf_pbmc_<analyte>` (see the next entry).
- **Source aliases:**
  - `K_TAF-TFVdp-PBMC` -- used in `Yu_2026_tenofovir.R` (Yu 2026 Results 3.2.1 / Table 2, "intracellular uptake and activation rate from plasma TAF to PBMC TFV-dp").
- **Example models:** `Yu_2026_tenofovir.R` (`K_TAF-TFVdp-PBMC = 1.59 1/h`, the cathepsin-A-mediated route that makes tenofovir alafenamide load PBMCs far more efficiently than tenofovir disoproxil fumarate does).
- **Notes:** The `kinf` stem is deliberately NOT `kin`: the canonical `kin` is a ZERO-ORDER production rate into an indirect-response turnover pool, a different quantity with different units. Distinct from `clin` / `clef` (`Campagne_2019_cyclophosphamide_mouse.R`), which are plasma-to-tissue influx / efflux CLEARANCES (volume / time). The `_pbmc` suffix marks the destination pool and must match the canonical compartment stem `pbmc`. Intracellular PBMC nucleotide-analogue pools are a recurring shape in the HIV treatment and PrEP literature (Baheti 2011, Chen 2016, Yu 2026), so this family is expected to be reused rather than re-invented per paper.

### lkinf_pbmc_<analyte>, kinf_pbmc_<analyte> (**canonical source-suffixed PBMC influx rate constant**)
- **Type:** paper-named-param
- **Role:** The `lkinf_pbmc` influx rate constant, disambiguated by the plasma SOURCE species that drives it, for models in which more than one circulating moiety feeds the same PBMC pool. `<analyte>` is a registered metabolite suffix from `compartment-names.md`. The bare `lkinf_pbmc` is reserved for the parent-driven influx in the same file, so a file that carries both an `lkinf_pbmc` and an `lkinf_pbmc_<analyte>` reads unambiguously: bare means parent-driven, suffixed means driven by the named analyte.
- **Source aliases:**
  - `K_TFV-TFVdp-PBMC` -- used in `Yu_2026_tenofovir.R` as `lkinf_pbmc_tfv` (Yu 2026 Results 3.2.1 / Table 2, "intracellular uptake and activation rate from plasma TFV to PBMC TFV-dp").
- **Example models:** `Yu_2026_tenofovir.R` (`lkinf_pbmc_tfv`, `K_TFV-TFVdp-PBMC = 0.0116 1/h`, two orders of magnitude below the parent-driven `lkinf_pbmc`; the two rate constants together give a PBMC pool that is about 98% tenofovir-alafenamide-driven after a tenofovir alafenamide dose).
- **Notes:** This is the PBMC realisation of the extension axis the `lkinf_rbc` entry anticipated ("a future model fitting two red-cell arms jointly would extend to `lkinf_rbc_<analyte>`"). Add such a suffixed sibling ALONGSIDE the bare entry, never as a replacement for it -- q2 of the same sidecar (answered B) ruled that a red-cell pool with a single influx keeps the registered bare `lkinf_rbc` even when the file's dosed parent is a different species, with the driving analyte stated in the `label()` and source-trace comment instead.

### lkeff_pbmc, keff_pbmc (**canonical first-order efflux / elimination rate constant out of the PBMC analyte pool**)
- **Type:** paper-named-param
- **Role:** First-order rate constant (1 / time) for loss of drug or active metabolite from the intracellular PBMC pool (`pbmc_<analyte>`), lumping efflux back to plasma, intracellular dephosphorylation or catabolism, and cell turnover into a single rate. Enters as the `- keff_pbmc * pbmc_<analyte>` term of the PBMC ODE. A pool has only one loss term, so this name stays bare even when the pool has two named influxes.
- **Source aliases:**
  - `Kel-PBMC` -- used in `Yu_2026_tenofovir.R` (Yu 2026 Results 3.2.1 / Table 2, "elimination rate from PBMC").
- **Example models:** `Yu_2026_tenofovir.R` (`Kel-PBMC = 0.0127 1/h`).
- **Notes:** The PBMC parallel of `lkeff_rbc`, and lumped in the same sense: `log(2) / keff_pbmc` is an intracellular residence half-life rather than a pure membrane-transport property. Yu 2026 quotes a PBMC tenofovir-diphosphate half-life of 60 h against `log(2) / 0.0127 = 54.6 h`, because the published figure is an effective half-life read from the simulated concentration-time profile rather than `log(2) / k`.

### lkeff_milk, keff_milk (**canonical first-order elimination rate constant out of the breast-milk compartment**)
- **Type:** paper-named-param
- **Role:** First-order rate constant (1 / time) for loss of drug out of the canonical `milk` compartment to a sink, lumping removal of milk from the breast by feeding or expression together with milk turnover into a single rate. Enters as the `- keff_milk * milk` term of the milk ODE, exactly as its `lkeff_rbc` / `lkeff_pbmc` siblings enter theirs. Use it only for a UNIDIRECTIONAL milk model in which nothing returns from milk to plasma; when the milk compartment exchanges back with `central`, the transfer is a `k_<from>_<to>` pair (`k_central_milk` / `k_milk_central`) or a `q_milk` / `pcmilk` clearance parameterisation instead.
- **Source aliases:**
  - `kmilk_e` -- used in `Baklouti_2026_amoxicillin.R` (paper symbol `k_milk,e`, "the elimination rate constant of amoxicillin from milk").
- **Example models:** `Baklouti_2026_amoxicillin.R` (`kmilk_e = 0.33 1/h`).
- **Notes:** Registered 2026-09-02 (operator sidecar `oare_PMC13206287` request-001 / response-001, question q1, option A). Extends the `keff_<pool>` efflux family from intracellular analyte pools (`rbc_<analyte>`, `pbmc_<analyte>`) to a secretory matrix compartment; the family's Role sentence already described the quantity verbatim. Distinct from `kmilkinf` (`Wattanakul_2024_primaquine_motherinfant.R`), which is the dyad milk-to-infant-depot transfer normally fixed high, and from `k_milk_central`, which returns drug to plasma. Lumped in the same sense as its siblings, and in this founding case deliberately so: Baklouti 2026's Discussion notes that "in our model, the drug is continuously eliminated from the milk compartment, whereas in reality, this elimination occurs only at the time of feedings", which makes the resulting relative-infant-dose estimate conservative. The ratio `k_central_milk / keff_milk` is the model's steady-state milk:plasma AUC ratio whenever `vmilk == vc`, which is the natural closed-form falsifier for an extraction using this name.

### lvmax_rbc, vmax_rbc (**canonical saturable-influx maximum rate into the red-cell analyte pool**)
- **Type:** paper-named-param
- **Role:** Maximum rate (concentration / time, e.g. umol/L/day) of SATURABLE Michaelis-Menten influx from plasma into the intracellular red-cell pool: `d/dt(rbc_<analyte>) = vmax_rbc * Cc / (km_rbc + Cc) - keff_rbc * rbc_<analyte>`.
- **Source aliases:**
  - `Vmm` -- used in `Gebhard_2023_mercaptopurine.R` (paper symbol `V_mm^6MP`).
- **Example models:** `Gebhard_2023_mercaptopurine.R` (`V_mm^6MP = 0.096 umol/L/day`), `Gebhard_2023_mercaptopurine_anc.R` (`V_mm^6MP = 0.21 umol/L/day`), `Morse_2012_ghb_rbc_invitro.R` (`Vmax = 20.9 nmol/mg protein/min at pH 7.4`, `5.3` at pH 6.5; carried as the experimental-condition-suffixed pair `lvmax_rbc_ph74` / `lvmax_rbc_ph65`).
- **Notes:** Extends the blessed `vmax_<suffix>` / `km_<suffix>` disambiguation pattern to a saturable INFLUX. Distinct from the canonical bare `vmax` / log `lvmax`, which are registered for saturable ELIMINATION and carry amount/time units; `vmax_rbc` is an influx into a concentration state and therefore carries concentration/time units -- except where the destination `rbc_<analyte>` state is itself normalised per mg of red-cell protein, as in an in vitro uptake assay, in which case it carries amount / mass protein / time (see `rbc_ghb` in `compartment-names.md`). **Experimental-condition suffixes.** When a source fits the same influx structure independently under two or more deliberately varied experimental conditions, and shares no parameter between them, append a condition token to the registered stem -- `lvmax_rbc_ph74` / `lvmax_rbc_ph65`, and likewise for `lkm_rbc_*` and `lkinf_rbc_*` -- and switch on the corresponding covariate inside `model()`. This is the same device `HernandezLozano_2025_apramycin_invitro.R` uses for its per-strain-and-pH drug-effect parameters, and it is distinct from the `lkinf_rbc_<analyte>` extension reserved for a model fitting two red-cell influx ARMS jointly: a condition suffix indexes one arm characterised repeatedly, an analyte suffix indexes genuinely different arms. Keep the bare registered name whenever only one condition was studied.

### lkm_rbc, km_rbc (**canonical saturable-influx half-saturation concentration for the red-cell analyte pool**)
- **Type:** paper-named-param
- **Role:** Michaelis constant (concentration, e.g. umol/L) of the saturable influx into the intracellular red-cell pool; the plasma concentration at which influx reaches half of `vmax_rbc`.
- **Source aliases:**
  - `Kmm` -- used in `Gebhard_2023_mercaptopurine.R` (paper symbol `K_mm^6MP`).
- **Example models:** `Gebhard_2023_mercaptopurine.R` (`K_mm^6MP = 0.016 umol/L`), `Gebhard_2023_mercaptopurine_anc.R` (`K_mm^6MP = 0.14 umol/L`).
- **Notes:** Paired with `lvmax_rbc`; both are meaningless alone.

### linieff, inieff (**canonical mid-therapy initialisation fraction for a turnover / maturation chain**)
- **Type:** paper-named-param
- **Role:** Dimensionless fraction of the turnover baseline (`rbase`) at which a turnover or maturation chain is initialised when the observation record starts DURING ongoing therapy rather than at drug-free baseline. The terminal state is initialised at `inieff * rbase` and the upstream chain states at the steady-state values implied by that terminal value, so the chain begins at a treatment steady state carrying an unobserved historical drug effect rather than at the untreated baseline.
- **Source aliases:**
  - `inieff` -- used in `Gebhard_2023_mercaptopurine_anc.R` (Gebhard 2023 main text: "we assume to have reached a treatment steady state with an additional parameter `inieff` describing the drug effect at the initial time point and `X_ma(0) = inieff * base`").
- **Example models:** `Gebhard_2023_mercaptopurine_anc.R` (`inieff = 0.87`, RSE 22%, CV 54%; Friberg chain initialised mid-maintenance-therapy for childhood ALL).
- **Notes:** Estimated with IIV in Gebhard 2023, which is what distinguishes it from a fixed initialisation convention: it absorbs the between-patient spread in how suppressed each patient's ANC already was when their record began. A value of 1 recovers the ordinary drug-free-baseline initialisation. Distinct from `rbase`, which is the untreated baseline itself.

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
- **Notes:** Kept on the bare (untransformed) linear scale rather than log-transformed because the published between-subject variability on it is ADDITIVE (Zhang 2023 Methods: "the placebo effect... is modeled to have an additive variability (i.e., PLAC_i = PLAC_TV + eta_i), with the assumption that the placebo effect in any individual may be higher or lower than the typical value... by an equal probability"), so an `etaiplac` enters as `iplac + etaiplac`. Distinct from `imax`, which is the MAXIMUM of a driver-dependent (sigmoid) inhibition; `iplac` is a driver-free constant. Distinct also from `lplac` as used in `Schoemaker_2018_levetiracetam.R`, which is a log-scale placebo effect on a seizure RATE rather than a fractional inhibition of an IDR input function. Later placebo-controlled IDR extractions in the IBD / CDAI family (certolizumab, ustekinumab, mirikizumab,...) should reuse this canonical.

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

### fe (**canonical fraction of drug excreted unchanged in urine**)
- **Type:** paper-named-param
- **Role:** Fraction of systemically available drug excreted unchanged in the urine; unitless, in [0, 1]. The renal share of total clearance, and therefore the natural coefficient for splitting a total clearance into a renal and a non-renal (usually hepatic-metabolic) arm: `CLi = CL * (fe * RENALFUNC_REL + (1 - fe) * HEPFUNC_REL)`. Normally held `fixed`, because it is a mass-balance property of the compound (a urinary-recovery or radiolabel ADME measurement, or -- where no human value exists -- an animal renal-clearance-to-total-clearance ratio) rather than a quantity the model estimates.
- **Source aliases:**
  - `NFE` / `nfe` -- "normal fraction excreted unchanged (in healthy volunteers)", the Edsim++ / MwPharm `ORG` organ-function-object symbol. Same orientation and scale, no transformation. The "normal / healthy volunteer" qualifier is a statement about the *study population* the value was measured in, not a different quantity: record it in the model's `description` and in the parameter's source-trace comment, following the register's convention of carrying population context in `covariateData` / `description` rather than in the canonical name.
  - `fe,u`, `Fe`, `fu,r` -- assorted paper spellings for urinary recovery of unchanged drug.
- **Example models:** `Visscher_2025_parathyroidHormone.R` (founding example; `fe <- fixed(0.3)`, the healthy-subject renal share of rhPTH(1-84) total clearance, taken from Hruska et al's canine estimate that renal clearance is on average 30% of total clearance, and used as the reference against which `RENALFUNC_REL` scales the renal arm of Visscher 2025 Eq. 3).
- **Notes:** Distinct from `fu` (fraction unbound in *plasma*, a protein-binding property) and from `fm` (fraction metabolised through a named pathway): `fe` is a whole-body excretion mass-balance fraction. Distinct also from the covariates `RENALFUNC_REL` / `HEPFUNC_REL`, which are per-subject organ-function multipliers -- `fe` sets how much of the clearance each of those multipliers governs, while the covariates set how far that arm departs from normal. Pairs naturally with both.

### fu (**canonical fraction unbound in plasma**)
- **Type:** paper-named-param
- **Role:** Fraction unbound in plasma; fixed unitless multiplier (in [0, 1]) used in cerebral-microdialysis-style CNS-distribution models where only free drug crosses the BBB.
- **Source aliases:** none.
- **Example models:** `Campagne_2019_cyclophosphamide_mouse.R`.
- **Notes:** Usually held fixed at the in-vitro equilibrium-dialysis-derived value. The BBB transfer term is `CLin * fu * Cp`. Distinct from `fe` (fraction excreted unchanged in urine, a whole-body mass-balance fraction). When the source reports fu as a *function of concentration* rather than a single number, use the saturable family `fumin` / `fumax` / `cup50` below instead of a scalar `fu`; the derived `model`-block quantity is still named `fu`.

### fumin, fumax, cup50 (**canonical saturable concentration-dependent unbound fraction**)
- **Type:** paper-named-param
- **Role:** Three-parameter description of an unbound fraction that *rises with concentration* because the binding protein saturates, as happens with drugs bound principally to alpha-1-acid glycoprotein (a low-capacity, high-affinity binder). The derived quantity inside `model` keeps the canonical name `fu`, computed on the total (not unbound) concentration in the units the source tabulates the binding data on:

  ```r
  fu <- fumin + (fumax - fumin) * cconc / (cup50 + cconc)
  ```

  `fumin` is the lower-asymptote unbound fraction approached at low concentration and `fumax` the upper asymptote approached at saturating concentration (both unitless, in [0, 1]); `cup50` is the total plasma concentration giving the half-maximal rise from `fumin` to `fumax`, and carries concentration units. All three are normally `fixed()`, because they are fitted to a printed in-vitro binding table rather than estimated from the PK data.
- **Source aliases:** none as distinct symbols -- sources almost always print a single column headed `fu` (or `fu,p`) tabulated against a concentration column, and leave the functional form implicit. Record the source's own concentration units in the `label()` for `cup50`, since papers move freely between mg/L and uM.
- **Example models:** `Jaiswal_2025_dordaviprone.R` (`fumin` 0.0283, `fumax` 0.0628, `cup50` 5.34 mg/L, fitted to the five printed fu-vs-concentration points in Jaiswal 2025 Table 1 and reproducing each within 4%).
- **Notes:** Load-bearing rather than cosmetic wherever apparent clearance is proportional to fu: in the founding example `CL/F = (CLint/F) * fu` exactly, so this one relationship reproduces the paper's own dose non-proportionality (its simulated CL/F rises about 16% from 125 to 625 mg, which the model reproduces to within one percentage point) and dropping it over-predicts the 625 mg AUC by about 26%. Distinct from `lbfu`, which is the slope of a *linearly time-varying* unbound fraction -- a different mechanism (fu drifting over time) and a different functional form; and distinct from the scalar `fu` above, which is the time- and concentration-invariant case. The names follow the Emax-family shape (`e0` / `emax` / `ec50`) but deliberately do not reuse it: `fumin` rather than `fu0` because the lower asymptote is a fitted limit rather than a value anchored at zero concentration, and `cup50` rather than `fuc50` because the quantity is a concentration, not a fraction.

### mic (**canonical minimum inhibitory concentration**)
- **Type:** paper-named-param
- **Role:** Minimum inhibitory concentration of the drug against the target pathogen, in the model's concentration units (note that mg/L and ug/mL are numerically identical, which is why antimicrobial papers move between them freely). Used as the denominator of the PK/PD indices that drive antimicrobial PK/PD-integration, probability-of-target-attainment and PK/PD-cutoff models: `AUC/MIC`, `fAUC/MIC`, `Cmax/MIC`, and as the threshold for `T>MIC` / `fT>MIC`. Normally held `fixed()`, because it is a measured property of the isolate (CLSI / EUCAST broth microdilution) rather than an estimated model parameter; papers that sweep a range of MICs expose it so the user can override it per simulation.
  Because the same isolate has different MICs by different reference methods, a model whose exposure-response was fitted against one susceptibility methodology must state which one in the `label()`; when a paper fits parallel exposure-response relationships against two methodologies, those are separate model files (see `Beredaki_2023_micafungin_clsi.R` / `..._eucast.R`).
- **Source aliases:**
  - `MIC` -- the near-universal paper symbol.
  - `MIC50`, `MIC90` -- population-distribution percentiles; use `mic` only for a single organism's own MIC, and record which percentile the value is in the `label` and the source-trace comment.
  - `MEC` -- minimum effective concentration, the echinocandin/mould analogue of the MIC, reported for filamentous fungi.
  - `IC90`, `IC90,plasma`, `IC50` -- the in vitro inhibitory-concentration analogue reported for intracellular protozoa and other organisms for which broth microdilution is not applicable, and which are assayed by image-based percent-inhibition readouts instead. Functionally identical in the model: a measured in vitro potency held `fixed()` and used as the threshold of a `T>threshold` index. Record the percent-inhibition level, the strain, the host-cell line, and any protein-binding correction in the `label()` and the source-trace comment, because an in vitro IC90 read in assay medium is not the same number as the plasma-equivalent threshold. Used in `Assmus_2025_benznidazole_mouse.R` (6.427 ug/mL: in vitro IC90 against *Trypanosoma cruzi* amastigotes, Tulahuen strain, in 3T3 host cells, corrected for protein binding to assay medium and mouse plasma).
- **Example models:** `Chen_2023_tilmicosin.R` (fixed at 0.25 ug/mL, the CLSI 2013 broth-microdilution MIC of *Pasteurella multocida* strain C44-15, serovar D:7; forms the `AUC24h/MIC` index that drives a sigmoidal Emax kill model), `Lallemand_2023_benzylpenicillin_horse.R` (fixed at 0.25 mg/L, the paper's concluded PK/PD cutoff, and overridden via `params = c(mic = ...)` to sweep 0.0625-2 mg/L; used both as the `fAUC/MIC` denominator and as the threshold integrated by the `t_above_mic` state), `Beredaki_2023_micafungin_clsi.R` (0.008 mg/L, CLSI M27 median MIC of *Candida albicans* CA 580) and `Beredaki_2023_micafungin_eucast.R` (0.016 mg/L, EUCAST E.Def 7.3 median MIC of the same isolate), both driving an `fAUC0-24/MIC` sigmoidal Emax.
- **Notes:** Keep the bare lower-case name -- do NOT log-transform (an `lmic` would imply the MIC is being estimated) and do NOT rename to a generic `thres` / `ec50`, which would lose the microbiological meaning and collide with the sigmoidal-PD canonicals. Distinct from `ec50` / `lec50`: `ec50` is the *fitted* index value producing half-maximal effect (and in an `AUC/MIC`-indexed model is unitless, or has units of time), whereas `mic` is a *measured* concentration used to form the index in the first place. Pair with `fu` when the index is defined on free rather than total drug: apply the unbound fraction to the exposure numerator and leave `mic` on the total-drug scale the susceptibility method actually reports, since the free-drug index `fu * Cc` against `mic` is equivalent to total `Cc` against `mic / fu`. For models that carry more than one pathogen or isolate, suffix per analyte in the usual way (`mic_<organism>`).

### f<matrix> (**canonical secondary-matrix:plasma equilibrium concentration ratio**)
- **Type:** paper-named-param
- **Role:** Unitless multiplicative scale factor converting a plasma / central-compartment concentration into the concentration observed in a second sampling matrix that is **not** modelled as a kinetically distinct compartment: `C<matrix> <- f<matrix> * Cc` (or `f<matrix> * <state> / v<state>` when the matrix has its own volume). The log-transform prefix is `lf<matrix>`; IIV, when the source paper estimates it, is `etalf<matrix>`. Named for the matrix, not the drug: `fsaliva`, `fcortical`, `fcancellous`. Use this whenever a paper concludes that a secondary matrix equilibrates rapidly enough with plasma that a scale factor outperforms a separate compartment -- the pairing with the `C<matrix>` output name and a `propSd_C<matrix>` residual is what makes the encoding self-describing.
- **Source aliases:**
  - `scale-factor` / `Scale-factor` -- Xu 2023 Table 3 and Figure 1 (saliva).
  - `bone:serum ratio` -- Landersdorfer 2009 Table 1 (cortical and cancellous bone).
- **Example models:** `Landersdorfer_2009_moxifloxacin.R` (founding example: `lfcortical` = log(0.803) and `lfcancellous` = log(0.775), the cortical- and cancellous-bone:serum equilibrium concentration ratios, each with estimated IIV), `Xu_2023_busulfan.R` (`lfsaliva` = log(0.88), the saliva:plasma ratio; no IIV, the authors' estimate shrank 100% and was dropped).
- **Notes:** Distinct from `sf<n>` (which corrects a bottom-up *predicted* PBPK partition coefficient and is cluster-indexed, not matrix-named) and from `kpu<n>` / `kp_<tissue>` (unbound tissue:plasma partition coefficients inside a mechanistic distribution model, which carry real tissue kinetics). Also distinct from `cal_slope_<assay>` / `cal_int_<assay>`: those recalibrate two *assays* measuring the same physical state onto a common scale and describe the measurement process, whereas `f<matrix>` is a biological distribution ratio between two genuinely different matrices and takes no intercept. `f<matrix>` deliberately implies **no** kinetics: the ratio is time-invariant, so the matrix profile is the plasma profile rescaled and the two are parallel on a log scale. Do not use it when the paper retains a separate compartment for the matrix with its own transfer and elimination rate constants -- that is an ordinary compartment plus a penetration ratio (see `csf`, `lesion`, `ecf` in `compartment-names.md`). Documented alongside the Xu 2023 busulfan saliva extraction; the pattern itself dates to the Landersdorfer 2009 moxifloxacin extraction.

### kpu<n> (**canonical clustered unbound tissue:plasma partition coefficient**)
- **Type:** paper-named-param
- **Role:** Estimated unbound tissue-to-plasma partition coefficient (Kpu) shared by a *cluster* of PBPK tissues, numbered `kpu1`, `kpu2`, `kpu3`, `kpu4`,... (log-transform prefix `lkpu<n>`). Unitless. Used by "steady-state commonality" PBPK simplifications that keep the full whole-body kinetic structure but reduce the number of unknown distribution parameters by assigning one common Kpu to every tissue in a composition-derived cluster. The numeric suffix indexes the cluster, not a compartment; the cluster-to-tissue mapping must be written out explicitly in `model()` (e.g. `kpu_bone <- kpu1`) and cited to the source table footnote. Convert to the tissue:blood coefficient with `kb = kpu * fu_p / BP`.
- **Source aliases:**
  - `KPU1` .. `KPU4` -- NONMEM `$PK` names in Yau 2023 Appendix S1.
  - `Kpu1` .. `Kpu4` -- Yau 2023 Table 2 / Table S6.
- **Example models:** `Yau_2023_diazepam_pbpk_kpu_human.R`, `Yau_2023_diazepam_pbpk_kpu_rat.R`, `Yau_2023_diazepam_pbpk_lumped_human.R`, `Yau_2023_diazepam_pbpk_lumped_rat.R`.
- **Notes:** Distinct from the per-tissue `kp_<tissue>` family (a separate Kp value named for one specific organ, as in the Gaohua 2012 pregnancy PBPK models) -- `kpu<n>` is deliberately cluster-indexed because the whole point of the parameterisation is that several anatomically distinct tissues share one estimate. Also distinct from `sf<n>` below, which scales a *predicted* Kpu rather than replacing it. In kinetically lumped models the same `kpu<n>` names index the lumped compartments (`kpu1` = central lump, `kpu2` = peripheral1 lump,...). Introduced 2026-07-26 with the Yau 2023 diazepam PBPK extraction.

### sf<n> (**canonical clustered Kpu scaling factor**)
- **Type:** paper-named-param
- **Role:** Estimated multiplicative correction factor applied to a *bottom-up predicted* unbound tissue:plasma partition coefficient, shared by a cluster of PBPK tissues and numbered `sf1`, `sf2`, `sf3`, `sf4`,... (log-transform prefix `lsf<n>`). Unitless. Implements `Kpu_i = Kpu_pred,i * SF_cluster(i)`, where `Kpu_pred,i` comes from a mechanistic tissue-composition model (Rodgers and Rowland, Poulin and Theil, Berezhkovskiy,...) evaluated inside `model()`. The scalar therefore *quantifies the systematic bias* of the bottom-up prediction instead of discarding it, which is what makes this parameterisation attractive for interspecies translation (the bias is assumed to be conserved across species while the composition inputs change).
- **Source aliases:**
  - `SF1`.. `SF4` -- Yau 2023 Table 2 / Table S6 and Eq 10.
- **Example models:** `Yau_2023_diazepam_pbpk_scalar_human.R`, `Yau_2023_diazepam_pbpk_scalar_rat.R`.
- **Notes:** A value of 1 means the bottom-up prediction needs no correction, so `sf<n>` estimates are directly interpretable as fold-bias. Do NOT reuse `sf<n>` for a generic "scaling factor" on a non-partition-coefficient parameter -- allometric exponents use `e_wt_<param>`, and covariate multipliers use the wider `e_<cov>_<param>` family that `e_wt_<param>` belongs to. Paired with the `kpurr_<tissue>` derived `model` quantities in the founding examples. Introduced 2026-07-26 with the Yau 2023 diazepam PBPK extraction.

### eh (**canonical hepatic extraction ratio**)
- **Type:** paper-named-param
- **Role:** Hepatic extraction ratio in the well-stirred liver model (`CL_H = FQ * eh`, `F_HEP = 1 - eh`). Unitless, bounded in [0, 1]. Use whenever a paper parameterises hepatic clearance through the extraction-ratio physiology rather than through `lcl` / `lcl_nonren` directly. The companion log-transform prefix is `logiteh` (logit-scale `ini()` typical value) because the linear-scale `eh` is bounded; logit-additive eta on `logiteh` keeps every individual `eh_i` inside the (0, 1) box.
- **Source aliases:** `EH` -- used in `Chan_2008_maraviroc.R` and Brussee 2018 mAb PBPK (the latter as a derived `model()`-block quantity rather than an `ini()` parameter, so the canonical applies to the `ini()` use case introduced by Chan 2008).
- **Example models:** `Chan_2008_maraviroc.R` (estimated `logiteh = logit(0.662)` with logit-additive IIV per Chan 2008 Eq 9; downstream `eh` enters `clh = fq * eh` and `fhep = 1 - eh`).
- **Notes:** Paired with the canonical compartment / pseudo-parameter `fq` (hepatic plasma flow, fixed at a literature value) and the canonical `clr` (renal clearance, often fixed) when the paper decomposes total CL into renal + hepatic with hepatic-extraction physiology. Distinct from `lcl_nonren` (additive renal + non-renal CL decomposition without an explicit extraction-ratio bound): use `eh` only when the source paper writes hepatic clearance as `CL_H = FQ * E_H` and constrains E_H in [0, 1] (e.g., via a logit-form IIV transformation). The gut-wall sibling is `egut` below; an oral model that resolves first-pass loss into its two anatomical sites carries both.

### egut (**canonical gut-wall extraction ratio**)
- **Type:** paper-named-param
- **Role:** Gut-wall (enterocyte) extraction ratio, the fraction of an absorbed oral dose lost to intestinal metabolism on first pass. Unitless, bounded in [0, 1]; gut availability is `F_GUT = 1 - egut`. The exact anatomical sibling of `eh`: together with the fraction absorbed they factor oral bioavailability as `F = f_a * (1 - egut) * (1 - eh)`. Use whenever a paper resolves first-pass loss into its gut and hepatic components rather than reporting a single lumped `F`, which is the norm for CYP3A4 substrates because the gut wall carries a large share of CYP3A4 activity.
- **Source aliases:** `FG` / `F_gut` -- the usual paper symbols for the *availability* `1 - egut` rather than the extraction ratio itself; convert before use. `EG` -- the extraction ratio directly.
- **Example models:** `Jaiswal_2025_dordaviprone.R` (`egut` 0.149 fixed, with `eh` 0.220, jointly pinned by the printed CL/F and the calibrated baseline bioavailability; both scale with a coadministered CYP3A4 modulator's relative activity through the well-stirred relation). `vandenBerg_2021_uprifosbuvir_pbpk.R` uses the identical name as a bare derived `model()` quantity computed from `egut <- (fug * cl_int_g) / (qgut + fug * cl_int_g)`.
- **Notes:** Registering this is largely a backfill -- the name was already shipping bare in `vandenBerg_2021_uprifosbuvir_pbpk.R` before it had a register entry. Keep `egut` for the extraction *ratio* and never for the availability, so that the `1 - e` convention reads the same way for gut and liver; a paper that prints only `FG` should be transcribed as `egut = 1 - FG` with the arithmetic shown in the source-trace comment. A DDI model scales the underlying intrinsic clearance and re-derives the ratio through the well-stirred relation (multiplying intrinsic clearance by A takes E to `E*A / (1 - E + E*A)`), rather than scaling the ratio directly, which would let availability leave [0, 1].

### clint (**canonical unbound intrinsic clearance**)
- **Type:** paper-named-param
- **Role:** Intrinsic clearance defined on *unbound* drug, i.e. the fu-independent term that must be multiplied by an unbound fraction to become a clearance: `CL = clint * fu`. Volume/time units. The log-transformed form `lclint` is the usual `ini` encoding. Use whenever the model needs clearance to respond to a varying unbound fraction -- for the founding example a *concentration-dependent* fu (see `fumin` / `fumax` / `cup50`), but equally a covariate-driven or disease-driven one. In an oral model written in apparent parameters the entry is `clint/F` and the label should say so.
- **Source aliases:** `CLint`, `CL_int`, `CLu,int` -- standard notation across PBPK and hepatic-clearance papers.
- **Example models:** `Jaiswal_2025_dordaviprone.R` (`lclint` = log(968.18) L/h apparent, set so the model returns the printed CL/F of 29 L/h at 125 mg), `Pei_2023_tacrolimus_pbpk.R` (`lclint` = log(11535)), `Mehta_2023_bedaquiline_mpbpk.R` (`lclint` = log(60.3) with allometric weight scaling and IIV `etalclint`).
- **Notes:** Distinct from `lcl` (total clearance already inclusive of binding) and from the `lcl_<component>` arms (`lcl_renal`, `lcl_nonren`, `lcl_2b6`,...), which are clearances in flow units that sum to a total: `clint` does not become a clearance until multiplied by fu, so the two are not interconvertible without knowing the binding. Pairs naturally with `eh` / `egut`, since the well-stirred model expresses each extraction ratio in terms of `fu * clint` against the relevant flow. **Name collision to be aware of:** `Rymut_2023_anti_tryptase.R` also carries an `lclint`, but there it denotes the internalisation clearance of a drug-target complex (L/day) and is *not* an unbound intrinsic clearance; that model predates this entry and is not retrofitted here.

### nowsmax (**canonical bare NAS-natural-history maximum-baseline term**)
- **Type:** paper-named-param
- **Role:** Multiplicative baseline of the natural neonatal-abstinence-syndrome (NAS) severity-decay-with-postnatal-age term `NOWST = nowsmax * exp(-nowsm * PNA_days)` used inside the Eudy-Byrne 2021 buprenorphine indirect-response PD model of MOTHER NAS scores. `nowsmax` is unitless (a multiplier on `kin/kout`) and enters both the ODE production term `kin * (1 + nowst)` and the drug-free quasi-steady-state initial condition `nows(0) = kin * (1 + nowst) / kout`. Log-transformed form is `lnowsmax`.
- **Source aliases:**
  - `NOWSMAX` -- Eudy-Byrne 2021 Table S3 notation (all-caps).
- **Example models:** `EudyByrne_2021_buprenorphine.R` (typical `nowsmax = 1.92`, unitless; 95% CI 1.76-2.08; log-normal IIV omega^2 = 1.14).
- **Notes:** Paired with `nowsm` (natural decay rate) and the canonical PD state / observation `nows`. Distinct from `rbase` (generic IDR baseline) because `nowsmax` multiplies an exponential decay with postnatal age rather than being a fixed steady-state baseline.

### lnowsmax (**canonical log-transformed NAS-natural-history maximum-baseline term**)
- **Type:** paper-named-param
- **Role:** Log-transformed form of `nowsmax` for `ini` with log-normal IIV. Inside `model` the bare name is `nowsmax = exp(lnowsmax + etalnowsmax)`.
- **Source aliases:**
  - `log NOWSMAX`, `lNOWSMAX` -- equivalent paper notation.
- **Example models:** `EudyByrne_2021_buprenorphine.R`.

### nowsm (**canonical bare NAS-natural-history decay rate with postnatal age**)
- **Type:** paper-named-param
- **Role:** First-order decay rate constant (1/day) of the natural NAS-severity term `NOWST = nowsmax * exp(-nowsm * PNA_days)` with chronological postnatal age. As PNA grows, NOWST decays toward zero and the drug-free quasi-steady-state NAS score `nows0 = kin * (1 + nowst) / kout` approaches its long-term floor `kin / kout`. Log-transformed form is `lnowsm`.
- **Source aliases:**
  - `NOWSM` -- Eudy-Byrne 2021 Table S3 notation (all-caps).
- **Example models:** `EudyByrne_2021_buprenorphine.R` (typical `nowsm = 0.107 1/day`; 95% CI 0.102-0.112; log-normal IIV omega^2 = 1.42; correlated with `nowsmax` at corr = 0.778).
- **Notes:** Units are 1/day. When paired with the canonical `PNA` covariate (which is in MONTHS), inside `model` compute `pna_days = PNA * 30.4375` before use so the exponent stays dimensionless (same pattern as Zhao 2018).

### lnowsm (**canonical log-transformed NAS-natural-history decay rate**)
- **Type:** paper-named-param
- **Role:** Log-transformed form of `nowsm` for `ini` with log-normal IIV. Inside `model` the bare name is `nowsm = exp(lnowsm + etalnowsm)`.
- **Source aliases:**
  - `log NOWSM`, `lNOWSM` -- equivalent paper notation.
- **Example models:** `EudyByrne_2021_buprenorphine.R`.

### cal_slope_<assay> (**canonical assay cross-calibration slope**)
- **Type:** paper-named-param
- **Role:** Slope of a linear recalibration mapping one measurement modality's predictions onto another modality's scale, used when a single structural state is observed by two assays with different bias / gain and the paper estimates the mapping as part of the model: `pred_reference = cal_slope_<assay> * pred_<assay> + cal_int_<assay>`. Unitless. `<assay>` is a short lower-case token naming the *non-reference* modality (`spect`, `pet`, `dbs`, `saliva`,...).
- **Source aliases:**
  - `beta` -- Siebinga 2023 Equation 1 (`Cpred = Cpred_SPECT * beta + alpha`).
- **Example models:** `Siebinga_2023_lu177psma617.R` (`cal_slope_spect` = 0.828, fixed; recalibrates SPECT/CT-derived blood activity onto the blood-sample scale).
- **Notes:** Distinct from a covariate effect (`e_<cov>_<param>`): the calibration parameters describe the *measurement* process, not a biological source of variability, and belong with the residual-error block rather than with the structural PK. Typically fixed after estimation in a data-source-only sub-model.

### cal_int_<assay> (**canonical assay cross-calibration intercept**)
- **Type:** paper-named-param
- **Role:** Intercept of the linear recalibration described under `cal_slope_<assay>`, in the units of the observation.
- **Source aliases:**
  - `alpha` -- Siebinga 2023 Equation 1.
- **Example models:** `Siebinga_2023_lu177psma617.R` (`cal_int_spect` = 6.27 MBq/L, fixed).

### cal_bias_<matrix> (**canonical structural measurement bias**)
- **Type:** paper-named-param
- **Role:** Additive structural measurement offset applied to the model prediction for a given sampling matrix before the residual-error model, in the units of the observation. Encodes a known / assumed systematic difference between the measured and the true quantity (assay calibration bias, matrix interference, background signal) that the paper estimates as a structural parameter rather than absorbing into the residual error. `<matrix>` names the sampling matrix (`blood`, `plasma`, `urine`,...).
- **Source aliases:**
  - `gamma` -- Siebinga 2023 Equation 2 (`Cobs = (Cpred + gamma) * (1 + eps_p) + eps_add`).
- **Example models:** `Siebinga_2023_lu177psma617.R` (`cal_bias_blood` = 0.273 MBq/L, fixed; the paper attributes it to calibration uncertainty from extreme calibration ranges for blood samples).
- **Notes:** A positive `cal_bias_<matrix>` raises the prediction, so it forces predictions above the drug-free baseline; Siebinga 2023 notes this is the source of the apparent under-prediction of low blood observations in its CWRES plots.

### kaff (**canonical equilibrium affinity constant of a saturable binding site**)
- **Type:** paper-named-param
- **Role:** Equilibrium *affinity* constant of a saturable (Langmuir) plasma- or blood-protein binding site, in reciprocal-concentration units (1 / concentration), so that `kaff * C` is dimensionless. It is the reciprocal of a dissociation constant (`kaff = 1 / kd`) and appears in the affinity-parameterised binding form `Ctot = Cu * (1 + bmax / (1 + kaff * Cu) + kns)`, where the saturable term falls from `bmax` at `Cu << 1 / kaff` toward zero at `Cu >> 1 / kaff`. Pairs with `bmax` / `lbmax` (capacity) and `kns` (linear non-saturable term).
- **Source aliases:**
  - `K1` -- Bouazza 2025 Table 2 and Petersen et al.'s prednisolone corticosteroid-binding-globulin (CBG) binding model; the same symbol appears in the wider corticosteroid-binding literature.
  - `Ka`, `K_A` -- general Langmuir / Scatchard notation for an association-equilibrium constant.
- **Example models:** `Bouazza_2025_prednisolone.R` (`kaff = fixed(0.0095)` L/nmol, CBG affinity fixed from Petersen et al.; used to recover total prednisolone from the unbound concentration the disposition is parameterised on).
- **Notes:** **Do NOT reuse `k1` for this role.** The registered `k1` is a second-order *association rate* constant with units 1 / (concentration * time), paired with `k2` and `kd = k2 / k1`; an equilibrium affinity constant has units 1 / concentration and no time dimension, so reusing `k1` would be a genuine role collision of exactly the kind that entry's Notes warn about. `kaff` and `kd` are algebraically interchangeable (`kd = 1 / kaff`, with the capacity rescaling `bmax_conc = bmax / kaff` when moving between the dimensionless-`bmax` affinity form and the concentration-`bmax` dissociation form), so a paper that prints the dissociation constant should use the already-registered `kd` / `lkd` and only a paper that prints the affinity constant should use `kaff`. Preferring the form the paper prints keeps every `ini()` value readable straight off the source table instead of requiring the reader to redo the algebra -- which is the reason this canonical was ratified rather than reparameterising Bouazza 2025 onto `bmax` + `kd`. Held `fixed()` whenever the constant is carried in from an external in-vitro or literature binding study rather than estimated from the popPK data.

### kns (**canonical linear non-saturable protein-binding constant**)
- **Type:** paper-named-param
- **Role:** Dimensionless constant of the *linear*, non-saturable arm of a saturable-plus-linear protein-binding model -- typically albumin binding, which does not saturate over the clinical concentration range, sitting alongside a saturable high-affinity site (CBG, alpha-1-acid glycoprotein, an erythrocyte binder). In the affinity form `Ctot = Cu * (1 + bmax / (1 + kaff * Cu) + kns)` it is the concentration-independent floor of the total:unbound ratio: as the saturable site saturates, `Ctot / Cu` approaches `1 + kns`. In the dissociation form `Cbound = bmax * C / (kd + C) + kns * C` it is the slope of the non-saturable bound-vs-free line.
- **Source aliases:**
  - `Kns` -- Bouazza 2025 Table 2 (prednisolone albumin binding, from Petersen et al.); also `TerHeine_2018_everolimus.R`, which hardcodes it as a literal inside `model` because there it is a pre-fit data transform rather than a parameter of the fitted model.
  - `Nplasma` -- `Sikma_2020_tacrolimus_unbound_plasma.R` used this name (as `lnplasma`) for the same linear non-saturable role. Note the offset differs by parameterisation: Sikma's `Nplasma` is the full total:unbound ratio of the non-saturable arm (`Ctpc = nplasma * Cupc`), whereas `kns` is that ratio minus the unbound contribution already carried by the leading `1` in the Bouazza / Petersen form. Convert accordingly rather than transcribing the number across parameterisations. Existing files are not migrated.
- **Example models:** `Bouazza_2025_prednisolone.R` (`kns = fixed(0.8)`, albumin binding fixed from Petersen et al.).
- **Notes:** Unitless by construction, so it takes no `l` prefix requirement and is written bare in `ini()`. Usually `fixed()` from an external binding study. Before proposing any *further* canonical for a "scale factor" in a protein-binding equation, check for the molecular-weight-ratio signature `MW_drug / MW_protein * 1e6`: `Said_2025_imatinib.R`'s `L = 11700` looks like a new binding parameter but is a unit conversion (493.6 / 42000 * 1e6), correctly hardcoded in `model()` with the arithmetic in a comment.

### kaff2 (**canonical equilibrium affinity constant of a second binding-site class**)
- **Type:** paper-named-param
- **Role:** Equilibrium affinity constant of the *second* class of binding sites in a two-class (biphasic) protein-binding model, in the same reciprocal-concentration units as `kaff`, so that `kaff2 * Cu` is dimensionless. Use it when a source paper prints both `Ka1` and `Ka2` for one protein -- the Rosenthal / Scatchard plot is biphasic -- and both constants are needed to reconstruct the binding. Follows the established numbering idiom for a second instance of one role (`vp` / `vp2`, `q` / `q2`): the first class keeps the bare `kaff` and the second takes the `2` suffix, so the two read straight off the paper's own `Ka1` / `Ka2` labels.
- **Source aliases:**
  - `Ka2` -- Li 2017 Table 1 (naproxen binding to rat plasma albumin, second class of sites).
- **Example models:** `Li_2017_naproxen_rat.R` (`kaff2_cf` = 0.0041, `kaff2_hf` = 0.0043, `kaff2_cm` = 0.0056, `kaff2_hm` = 0.0054 L/umol, all `fixed()` from Table 1 and selected in `model()` by `SEXF` and `DIS_CIA`; paired with `kaff_*` for `Ka1`).
- **Notes:** Ratified 2026-09-02 (task `oare_PMC5399645` sidecar request-001 q2, answer A). Naming a second affinity constant is what lets the *capacity* of each class stay proportional to a live protein covariate: in Li 2017 the linear arm is `kns = n2 * Pt * Ka2`, built inside `model()` from the canonical `ALB` column, so freezing the arm as a per-group `kns` constant would have made the paper's own mechanism -- albumin drives the nonlinearity -- unexercisable. Check the `kns` entry's warning before minting any further binding-equation scalar: a "constant" of the form `MW_drug / MW_protein * 1e6` is a unit conversion and belongs hardcoded in `model()`, not in `ini()`. A second class that the source paper's own approximation collapses to a *linear* arm (Li 2017's stated `Ka2 * Cup << 1`) still keeps `kaff2` rather than being folded into a bare `kns`, because the printed table value is the affinity constant. Do not extend past `kaff2` without a sidecar: a third class of sites is rare enough that the naming should be decided against the paper that needs it.

### f_alb_isf (**canonical interstitial-fluid-to-plasma albumin concentration ratio**)
- **Type:** paper-named-param
- **Role:** Unitless ratio of the albumin (binding-protein) concentration in tissue interstitial fluid to that in plasma, used to scale a plasma-derived binding capacity to the interstitial space so that one set of binding constants serves both compartments. Sits in the `f_` family of unitless fractional scalars (`f_pah`, `f_ntcp`, `f_venous`, `f_exposed_bone`) and is always below 1 in the published values. Distinct from the `kp_<tissue>` family, which is an equilibrium tissue:plasma ratio of the *drug*; `f_alb_isf` is a ratio of the *binding protein*. Sibling ratios for other proteins or spaces take the same shape, `f_<protein>_<space>`.
- **Source aliases:**
  - `E/P` -- Li 2017 Methods (after eq. 4) and Discussion; also the notation of Fleishaker and McNamara 1985 and Day et al. 1995, from which the human values are quoted.
- **Example models:** `Li_2017_naproxen_rat.R` (`f_alb_isf_cia` = 0.9, `f_alb_isf_healthy` = 0.5, both `fixed()` from the literature and selected by `DIS_CIA`; multiplies the albumin concentration in both the saturable capacity `n1 * Pt` and the linear arm `n2 * Pt * Ka2` when eq. 4 is re-solved for the tissue compartment).
- **Notes:** Ratified 2026-09-02 (task `oare_PMC5399645` sidecar request-001 q1, answer A). Typically `fixed()` from physiologic literature rather than estimated, since interstitial albumin is not observed in a PK study; a paper that does estimate it should still use this name, on the log scale (`lf_alb_isf`) only if it reports an exponential typical-value form. The ratio is disease-sensitive -- inflammation raises microvascular permeability and therefore interstitial protein -- so a model spanning inflamed and healthy cohorts is expected to carry two values selected by a disease indicator, as Li 2017 does (0.9 vs 0.5 in rats; the paper quotes 0.73 vs 0.32 for human RA patients and healthy subjects). Rejected alternative recorded so it is not re-opened: `kp_alb_isf` would overload `kp_<tissue>`, which the register defines for drug tissue:plasma ratios, and the register additionally bounds that family away from bare `kp_*` names that qualify rather than name a tissue.

### kurine (**canonical urinary-excretion rate constant**)
- **Type:** paper-named-param
- **Role:** First-order rate constant (1 / time) moving parent drug or a metabolite out of a plasma compartment and into the matching canonical `urine` / `urine_<metab>` cumulative-excretion compartment, used when the source paper carries urinary excretion as a primary estimated rate constant rather than deriving it from a renal clearance. Per-metabolite members take the registered metabolite suffix: `kurine_<metab>`.
- **Source aliases:**
  - `K14`, `K35` -- Thoueille 2026 Data S1 `$PK` / `$DES` (parent and metabolite urinary excretion respectively); Table 2 prints them as `k14` and `k35`.
- **Example models:** `Thoueille_2026_salmeterol.R` (`kurine = 0.000943 1/h` for salmeterol and `kurine_ohsal = 0.0147 1/h` for alpha-hydroxysalmeterol, each feeding its own cumulative urine compartment; the parent constant carries an athlete covariate effect and its own IIV; founding example). Also has a log-transformed primary form `lkurine` / `lkurine_<metab>`.
- **Notes:** Named for the destination compartment it feeds, mirroring how the registered biliary canonical `kbm` names the biliary route. **Do NOT reparameterise onto the registered `lcl_renal`.** The two are algebraically related (`cl_renal = kurine * vc`) but not interchangeable when the paper estimates the rate constant directly: a derived `kurine = cl_renal / vc` inherits the volume's random effect, so in the founding example the IIV on the excretion constant would inflate from variance 0.081 to 0.081 + 0.0262 and the model would no longer reproduce the published parameter table. Use `lcl_renal` when the paper reports a renal *clearance* and `kurine` when it reports a *rate constant*. The alternative spellings `kue` (opaque) and `kren` (invites confusion with `cl_renal`, a different quantity) were both considered and rejected.

### urprod (**canonical urine-production rate**)
- **Type:** paper-named-param
- **Role:** Zero-order rate of urine production (volume / time) driving the canonical `urine_vol` compartment, so that `d/dt(urine_vol) <- urprod` accumulates the volume voided since the last micturition and a urinary concentration can be formed as `urine / urine_vol`. Use when the paper *estimates* the production rate from the data; the log-transformed primary form is `lurprod`.
- **Source aliases:**
  - `UR_PROD` -- Thoueille 2026 Data S1 `$PK` / `$DES` and Table 2.
- **Example models:** `Thoueille_2026_salmeterol.R` (`urprod = 0.079 L/h`, i.e. 79 mL/h, estimated jointly with the PK to approximate physiologic micturition for the studies that did not record urine volumes; IIV of 73% CV retained only for records whose concentration was not corrected for urine specific gravity; founding example).
- **Notes:** Registered in the unpunctuated bare form `urprod` rather than the underscored `ur_prod`, consistent with the register's bare-PK style (`fdepot`, `bmax`, `vmax`) and with the `fumin` / `fumax` ruling. `Heuberger_2018_salbutamol.R` writes a `ur_prod_h` local inside `model()` for the same physical quantity, but there it is *derived* from cardiac-output physiology rather than estimated, so it is a `model()` intermediate and not an `ini()` parameter; that file is not migrated. **Do NOT reuse the registered `kpro`**: that canonical is a first-order biological-synthesis rate constant for a paper-mechanistic production term, and overloading it with a zero-order urine-volume rate would give one name two meanings and two dimensionalities -- the pattern rejected for `enzyme_liver`, `ppc` and `thres`.

### sdep_gsh (**canonical glutathione-pool depletion scaling factor**)
- **Type:** paper-named-param
- **Role:** Founding member of an `sdep_<pool>` family -- "Scaling factor coupling metabolic flux to DEPletion of a named co-substrate pool". It converts the rate at which a drug is metabolised into the rate of change of the paired baseline-normalised `<pool>_pool` compartment, in the form `d/dt(gsh_pool) <- (sdep_gsh / vc) * gsh_pool * kel * central`. Units are volume per amount (L/mg in the founding example), since the pool itself is dimensionless and the flux `kel * central` has units of amount per time.
- **Source aliases:**
  - `S_GSH` / `SGSH` -- Cao 2025 Equation 2, Table 2, and Supplementary Text S2 `$PK` (`TVS_GSH = THETA(9)`).
- **Example models:** `Cao_2025_busulfan.R` (`sdep_gsh = fixed(0.00259)` L/mg, carried over from Langenhorst 2020, with the covariate effect `e_gst_sdep_gsh = 0.28` on baseline GST activity).
- **Notes:** Written bare (no `l` prefix) because it is a mechanistic constant reported as a point value rather than an estimated positive-constrained parameter, per the "Endogenous / mechanistic parameters" guidance. Two rejected alternatives are recorded here so the next reader does not re-open them: `s_gsh` (verbatim from the paper) founds an `s_` prefix with no precedent in this register AND collides visually with NONMEM `$PK` compartment *scale* factors, which appear in this very control stream as `S1 = V1/1000`; and `kdep_gsh` would sit visually with the registered `kdeg` / `kdes` / `kint` / `krel` family, but every member of that family is a first-order RATE constant with units 1/time, whereas this parameter is a scaling factor with units L/mg -- naming it `k*` would be a units lie. Pairs 1:1 with the `gsh_pool` compartment (see `compartment-names.md`); a new co-substrate takes a matched `sdep_<pool>` / `<pool>_pool` pair. The covariate effect follows the ordinary `e_<cov>_<param>` grammar with a shortened covariate token, `e_gst_sdep_gsh`, matching the `e_dpp4_bmaxc` precedent in `Retlich_2015_linagliptin.R` rather than spelling out the full `GST_BL_NMOL_MIN_ML` column name.

### kon_placebo, koff_placebo (**canonical placebo kinetic-system onset and offset rate constants**)
- **Type:** paper-named-param
- **Role:** The two first-order rate constants (1 / time) of a hypothetical placebo "kinetic" system: `kon_placebo` moves a dimensionless dummy dose out of `depot_placebo` into `central_placebo`, and `koff_placebo` removes it from `central_placebo` (see the `depot_placebo` / `central_placebo` pair in `compartment-names.md`). Together they set the time-course -- onset, peak and washout -- of the placebo term of the effect model, whose magnitude is set separately by `slope_placebo`. Strictly positive, so `ini()` normally carries the log-transformed `lkon_placebo` / `lkoff_placebo` form; use it whenever the source reports the IIV as a percent CV.
- **Source aliases:**
  - `k_on`, `KON` -- "rate constant for onset of placebo response" (Ravva 2015 Supplementary Table S1).
  - `k_off`, `KOFF` -- "rate constant for offset of placebo response" (Ravva 2015 Supplementary Table S1).
- **Example models:** `Ravva_2015_varenicline.R` (`lkon_placebo = fixed(log(0.0112))`, `lkoff_placebo = fixed(log(0.130))`, both fixed from a placebo-data-only fit before the combined drug fit).
- **Notes:** Ratified 2026-09-03 alongside the Ravva 2015 varenicline craving extraction (sidecar `oare_PMC4832970` request-001 q1, operator answer B), spelled with the full `placebo` token to match `fpg_placebo` / `hba1c_placebo` / `skit_pla`; no `_pbo` suffix exists in these registers. The `_placebo` qualifier is what keeps these clear of the register's deliberate refusal to canonicalise a bare `kon` (see "Why `kon` is NOT canonicalised" below) -- these two carry exactly one dimensionality, 1 / time, because the driven quantity is a dimensionless dummy amount. Distinct from `k1` / `k2`, which are the association and dissociation rate constants of a real reversible binding reaction and whose association member is second-order. Distinct from `lkplacebo` in `Johnston_2019_empagliflozin.R` and `kpbo` in `Nguyen_2025_valbenazine_tmc.R`, which are single onset-rate constants of an *algebraic* asymptotic-exponential placebo term with no kinetic states behind them; use `kon_placebo` / `koff_placebo` only when the paper actually builds the two-state dummy-dose system, because only then does the rise-and-fall shape exist.

### slope_placebo, slope_drug (**canonical linear placebo and drug effect slopes**)
- **Type:** paper-named-param
- **Role:** The two coefficients of an additive linear effect model in which the response decomposes as baseline plus placebo plus drug, `E = e0 + slope_placebo * <placebo driver> + slope_drug * <drug concentration>`. `slope_drug` multiplies a real drug concentration, so its units are effect per concentration unit; `slope_placebo` multiplies the placebo kinetic-system driver `central_placebo`, so its units are effect per unit dummy amount and are only interpretable together with the unit-dose convention of that system. Both are kept on the bare, untransformed scale: published values are typically negative (an improvement), which forbids a log transform, and the accompanying IIV is additive.
- **Source aliases:**
  - `PSLP`, "placebo response parameter" -- Ravva 2015 Supplementary Table S1.
  - `DSLP`, "drug response parameter" -- Ravva 2015 Supplementary Table S1.
- **Example models:** `Ravva_2015_varenicline.R` (`slope_placebo = -3.32` and `slope_drug = -0.192`, both on the logit scale of a 0-100 craving score, each with an additive eta).
- **Notes:** Ratified 2026-09-03 alongside the Ravva 2015 varenicline craving extraction (sidecar `oare_PMC4832970` request-001 q2, operator answer A: `slope_drug`, with its partner spelled `slope_placebo` to match the full-token placebo convention). The two are registered as a pair because the naming only makes sense as a pair -- a model with a single linear drug slope and no placebo arm should use the plain canonical `slope`. An additive eta is load-bearing rather than incidental when the typical slope is negative but small relative to its IIV: it lets a subject's individual slope take either sign, which is what per-subject IPRED panels of such papers show. When a paper fits several parallel linear models of different endpoints, suffix by endpoint (`slope_drug_<endpoint>`) the same way multi-output residual SDs are suffixed.
### sd_ratio_cl_m5, sd_ratio_kgrow (**canonical ratio of two random-effect standard deviations sharing one eta**)
- **Type:** paper-named-param
- **Role:** Multiplier applied to a shared eta so that a second parameter carries a *different* IIV magnitude while remaining perfectly correlated with the first. Used when a source paper cannot estimate a second variance reliably (zero gradients, a near-singular OMEGA block, or a correlation that estimates to 1) and therefore reduces the model to one random effect plus a scale factor. The family shape is `sd_ratio_<param>`, where `<param>` is the canonical name of the parameter whose IIV is being *rescaled*; the eta that appears in the expression belongs to the other, reference parameter. Carried on the LINEAR scale, and encoded as `<param> <- exp(l<param> + sd_ratio_<param> * eta<reference>)`. Because the ratio is of standard deviations, the implied variance is `sd_ratio^2` times the reference variance -- the identity papers usually print as a footnote.
- **Source aliases:**
  - `Factor IIV CL_M-5` -- Keunecke 2020 Table 3 row label.
- **Example models:** `Courlet_2023_cabamiquine.R` (`sd_ratio_kgrow = fixed(0.9171881)`, the liver-to-blood growth-rate IIV SD ratio, derived as `sqrt(0.0120274 / 0.0142973)` from Table 1 with the correlation fixed at 1; founding example), `Keunecke_2020_regorafenib_phase3.R` (`sd_ratio_cl_m5 = 2.21`, RSE 3.41%, so that `cl_m5 <- exp(lcl_m5 + sd_ratio_cl_m5 * etalcl_m2) * sex_cl_m5`; the Table 3 footnote's `0.385 * 2.21^2 = 1.88` is the implied M-5 variance; doi:10.1111/bcp.14334).
- **Notes:** Heading added 2026-09-02 alongside the Keunecke 2020 regorafenib extraction, registering a shape that `Courlet_2023_cabamiquine.R` had already shipped unregistered -- per the register's "extend, never duplicate" rule, rather than founding a second name for the same construct. Distinct from a correlated omega block (`etaA + etaB ~ c(...)`), which estimates the covariance: `sd_ratio_<param>` is the *reduced* form used when that block is not identifiable, and it forces the correlation to exactly 1, so the second parameter's eta is deterministic given the first. Prefer the full omega block whenever the paper reports one; reach for this name only when the paper itself reports the reduction. Do not use `f_` (reserved for fractions and bioavailability) or a `k`-prefix (reserved for rate constants with 1/time units) -- the quantity is a unitless ratio.

## Neural-ODE (NODE) network parameters

Registered 2026-09-02 with `Bram_2026_warfarin_node.R`, the first neural-ODE
model in this library (operator sidecar `oare_PMC13291805` request-001 q1 = A).

A low-dimensional neural ODE puts a small feed-forward neural network on the
right-hand side of a differential equation in place of a mechanistic function.
The network's fitted weights and biases are ordinary population parameters and
need canonical names. The grammar is `nn_<slot>_<role>_<j>`:

| Token | Meaning |
|---|---|
| `nn_w1_<role>_<j>` | input-to-hidden weight for hidden neuron `j` |
| `nn_b1_<role>_<j>` | bias of hidden neuron `j` |
| `nn_w2_<role>_<j>` | hidden-to-output weight for hidden neuron `j` |
| `nn_beta_<role>` | softplus sharpness of the activation function |

`<role>` is the network's **role label** -- what the network computes -- not the
paper's symbol, matching how the rest of this register names roles rather than
notation. `<j>` indexes the hidden neuron, starting at 1. Extending the family
to a new network needs no fresh sidecar so long as the role token is
descriptive and the slot tokens above are reused verbatim; a network with a
different architecture (more than one hidden layer, a non-softplus activation)
does need one, because the slot vocabulary would have to grow.

These are `paper-named-param`, not `bare-pk`: weights and biases are **signed**,
take negative values, and must never be log-transformed or carry multiplicative
IIV. Their inter-individual variability is additive on the natural scale
(Monolix `distribution=normal`), so the eta is `etann_w1_<role>_<j>` and the
`model()` line is `w1_<role>_<j> <- nn_w1_<role>_<j> + etann_w1_<role>_<j>`.

Do **not** confuse the `nn_` prefix with `lnn` / `nn_fix`, which is a
sigmoidicity exponent in specific BDE / morphine-like models, or with `ntr`,
the transit-chain length. The trailing `_<role>_<j>` slot structure is what
distinguishes a network parameter from either.

### nn_w1_rc_1, nn_w1_rc_2, nn_w1_rc_3, nn_w1_rc_4, nn_w1_rc_5 (**canonical NODE input-to-hidden weights, response-compartment production network**)
- **Type:** paper-named-param
- **Role:** Input-to-hidden-layer weights of the five-neuron network that replaces the production term of an indirect-response (`rc` = response-compartment) model. Unitless. In the pmxNODE encoding the weight enters as `-(w^2)`, i.e. squared and negated to enforce a monotone-decreasing response, so only the magnitude is identifiable.
- **Source aliases:**
  - `Wrc_11` .. `Wrc_15` -- pmxNODE / Monolix names in the Braem 2026 Data S2 deposit.
- **Example models:** `Bram_2026_warfarin_node.R`.

### nn_b1_rc_1, nn_b1_rc_2, nn_b1_rc_3, nn_b1_rc_4, nn_b1_rc_5 (**canonical NODE hidden-layer biases, response-compartment production network**)
- **Type:** paper-named-param
- **Role:** Hidden-layer biases of the same five-neuron network. Unitless.
- **Source aliases:**
  - `brc_11` .. `brc_15` -- pmxNODE / Monolix names in the Braem 2026 Data S2 deposit.
- **Example models:** `Bram_2026_warfarin_node.R`.

### nn_w2_rc_1, nn_w2_rc_2, nn_w2_rc_3, nn_w2_rc_4, nn_w2_rc_5 (**canonical NODE hidden-to-output weights, response-compartment production network**)
- **Type:** paper-named-param
- **Role:** Hidden-to-output-layer weights of the same five-neuron network. These carry the units of the quantity the network returns -- here a production rate, so response units per time.
- **Source aliases:**
  - `Wrc_21` .. `Wrc_25` -- pmxNODE / Monolix names in the Braem 2026 Data S2 deposit.
- **Example models:** `Bram_2026_warfarin_node.R`.

### nn_beta_rc (**canonical NODE softplus sharpness, response-compartment production network**)
- **Type:** paper-named-param
- **Role:** Sharpness parameter of the softplus activation `softplus(x) = log(1 + exp(beta * x)) / beta`, which approaches a rectified linear unit as beta grows. Unitless, and normally fixed rather than estimated.
- **Source aliases:**
  - `beta` -- Braem 2026 Supporting Information Equation 2.
- **Example models:** `Bram_2026_warfarin_node.R` (`nn_beta_rc = fixed(20)`, the pmxNODE default).

### lkadair1, kadair1, lkadair2, kadair2 (**canonical constants of the Adair biphasic effect function**)
- **Type:** paper-named-param
- **Role:** The two denominator constants of an Adair biphasic (bell-shaped) drug-effect function, `Sd = Smax * C / (kadair1 + C + kadair2 * C^2)`. `kadair1` is the LINEAR denominator constant and carries concentration units; `kadair2` is the QUADRATIC self-inhibition constant and carries RECIPROCAL concentration units, and it is the term that creates the descending limb, so the effect peaks at `C = sqrt(kadair1 / kadair2)`. Both are strictly positive and are log-transformed in `ini()` as `lkadair1` / `lkadair2`, bare inside `model()`. `Smax` uses the existing `lsmax` / `smax` naming. Use this family for any biphasic concentration-response of Adair or substrate-inhibition shape, not only for the glucose-insulin case.
- **Source aliases:**
  - `k1` / `k2` -- Gao 2012 eq. 7 and Table 4 ("First receptor binding constant", "Second receptor binding constant").
- **Example models:** `Gao_2012_exenatide_glucose_insulin_rat.R` (`k1 = 0.826 nM`, `k2 = 0.0153 1/nM`; exendin-4's insulinotropic effect rises to a maximum near 7.4 nM and falls at higher concentrations, matching the bell-shaped in vitro and in vivo dose-response the paper cites; founding example).
- **Notes:** Ratified 2026-09-02 (task `oare_PMC3336795` sidecar request-001 q1, answer A). The `adair` infix is load-bearing and deliberately NOT the paper's bare `k1` / `k2`: this register already carries `k1` and `k2` as the ASSOCIATION and DISSOCIATION rate constants of reversible binding, which are 1/(conc*time) and 1/time, so reusing them here would give one name two dimensionalities -- and the numbered pair would additionally be read as the `k12` / `k21` micro-constants. The two constants are numbered rather than role-named (`kadair` / `kadair_inh` was the rejected option B) to keep the direct visual mapping onto the paper's printed `k1` / `k2`; the role is recorded in the `label()`. Reusing `km` for `kadair1` (rejected option C) was declined because `km` denotes a true Michaelis constant of a saturable process, whereas `kadair1` is a denominator coefficient of a stimulation factor. Note that the paper's own mechanistic reading of the two constants (Discussion: a 1:1 complex that may dissociate or bind a second drug molecule) is an interpretation of an empirical fit, not a separately identified binding model -- the TMDD `kon` / `koff` of the same paper are the real binding constants.

### lsstim_&lt;driver&gt;_&lt;target&gt;, sstim_&lt;driver&gt;_&lt;target&gt; (**canonical linear inter-pool stimulation factor**)
- **Type:** paper-named-param
- **Role:** Linear stimulation slope of one turnover pool on another pool's turnover, entering as `(1 + sstim_<driver>_<target> * (<driver> - <driver_baseline>))`. Units are the reciprocal of the DRIVER's concentration units, so two members of the same feedback loop routinely carry different scales. `<driver>` and `<target>` are the canonical compartment / analyte stems from `compartment-names.md`, named in that order, so the direction of the effect is readable without consulting the source. Log-transformed in `ini()` as `lsstim_<driver>_<target>`, bare inside `model()`. Pairs naturally with the `lkin_<analyte>` / `lkout_<analyte>` turnover family.
- **Source aliases:**
  - `SGlu` / `SIns` -- Gao 2012 eqs. 5-6 and Table 4 ("Stimulation factor of glucose on insulin secretion", "Stimulation factor of insulin on glucose disposal").
- **Example models:** `Gao_2012_exenatide_glucose_insulin_rat.R` (`sstim_glucose_insulin = 0.0684 1/mM`, `sstim_insulin_glucose = 0.157 1/nM`; founding example).
- **Notes:** Ratified 2026-09-02 (task `oare_PMC3336795` sidecar request-001 q2, answer A). The paper's own short names `sglu` / `sins` (rejected option B) were declined because they do not say which pool is driver and which is target -- fatal in a feedback loop where both pools are simultaneously both. `lslope_<target>` (rejected option C) was declined for the same ambiguity and because it would overload `slope`, which in this register denotes a concentration-response LINE (`lslope0` / `lslope_inf` in `Liesenfeld_2006_dabigatran_ECT.R`) or an assay cross-calibration slope (`cal_slope_<assay>`), neither of which is an indirect-response coefficient. Do NOT use this family for a drug's effect on a pool: a drug-driven stimulation of production or loss is `lsmax` / `lemax` with its own potency term. This family is specifically for a BIOMARKER driving another biomarker's turnover. When the driver's elevation is rectified at zero (the effect switching off below baseline), that gating belongs in `model()` as a `(delta > 0)` factor, not in a second parameter.

## Nested (multi-level) random effects

model in this library to carry a second hierarchical level of random effects.
This section documents a naming pattern, not a canonical parameter, so it has
no H3 entries and contributes nothing to the runtime name vectors.

When a paper fits a **second hierarchical level** of random effects on top of
the per-subject etas -- a between-study, between-site or between-centre random
effect whose draw is shared by every subject in that group -- name the nested
eta `eta<lparam>_<level>`, where `<lparam>` is the same transformed parameter
name the subject-level eta uses and `<level>` is a short lowercase token naming
the grouping level. Declare it in `ini` with rxode2's native nesting syntax,
piping the grouping column after the variance:

```r
etalcl ~ 0.112896 # subject-level IIV on CL/F
etalcl_study ~ 0.066049 | SIDN # study-level (nested) random effect on CL/F
```

and add it to the same log-scale sum as the subject-level eta in `model`:

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
with only one level `rxSolve` fails with `The following parameter(s) are
required for solving: THETA[2], THETA[1]`, an error naming thetas the model does
not have and giving no hint that the nesting is the cause. And `omega` must be
passed explicitly (`omega = mod$omega`), because a nested omega is a *list* of
matrices keyed by level rather than a single matrix.

### boxcox_<param> (**canonical Box-Cox shape parameter for an interindividual random-effect distribution**)
- **Type:** paper-named-param
- **Role:** Estimated shape parameter (lambda) of a Box-Cox transformation applied to the *interindividual* random effect of `<param>`, used when a paper reports that the eta distribution deviates visibly from log-normal. Unitless and signed. Because rxode2's `boxCox` attaches to the residual-error model and cannot be applied to an eta, the transform is written out explicitly in `model` in the Petersson (2009) form `phi_<param> <- (exp(eta<param>)^lambda - 1) / lambda`, followed by `<param> <- exp(l<param> + phi_<param>)`. A Box-Cox shape is a structural THETA, so it belongs in `ini` (wrapped in `fixed` only when the paper fixed it).
- **Source aliases:**
  - `Box-cox shape parameter` -- Siebinga 2024 Table 3.
  - `LAMBDA`, `SHAPE` -- common NONMEM `$THETA` labels for the same quantity.
- **Example models:** `Siebinga_2024_lu177psmaIT.R` (`boxcox_lkgrow` = -0.822, RSE 28.3%, the shape of the PSA-growth-rate eta; the authors added it because the kG random effect deviated strongly from log-normal, dOFV -28.7).
- **Notes:** Name the entry after the parameter whose eta is transformed, spelled as that parameter appears in `ini()`, so `boxcox_lkgrow` transforms `etalkgrow`. **Distinct from `lambda_<output>`** -- e.g. `lambda_ca125` in `Wilbaux_2014_ovarianCancer_ca125.R`, which is a Box-Cox transform of the *residual error*, declared through `rxode2::boxCox()` inside the error model and carried in the `paper_specific_residual_sds` metadata field. That mechanism cannot express an eta transform, which is why this separate canonical exists; do not conflate the two.

### lkd_direct, kd_direct, lkd_delay, kd_delay (**canonical linear drug-effect coefficients on a turnover loss rate**)
- **Type:** paper-named-param
- **Role:** Proportionality coefficients of a *linear* drug effect that multiply a driving exposure to give a first-order loss rate out of a PD turnover compartment (`E = kd_<route> * C`, entering the response ODE as `- E * R`). Units are reciprocal exposure per unit time and therefore depend on the driver's units, so the `label` must state them explicitly. `kd_direct` is driven by the instantaneous concentration in the source compartment; `kd_delay` is driven by an effect-compartment concentration equilibrating at `ke0`, and is the form to use when a paper reports response continuing after the drug has washed out. Both are strictly positive, so the log-transformed `lkd_<route>` form carries the `ini` typical value and any exponential IIV.
- **Source aliases:**
  - `kD,direct` -- Siebinga 2024 Table 3 and Eq. 9.
  - `kD,delay`, `kD,delayed` -- Siebinga 2024 Table 3 and Eq. 12.
- **Example models:** `Siebinga_2024_lu177psmaIT.R` (`lkd_direct` = log(0.00335) L/day/GBq driven by the tumor radioactivity concentration, and `lkd_delay` = log(0.0000328) L/day/MBq driven by an effect compartment at `ke0` = 0.00128 1/h; both act on the PSA compartment, per Eq. 12).
- **Notes:** **Distinct from `kd`**, the mechanistic *dissociation* rate (1 / time) used in TMDD-type models: these are exposure-scaled slopes with different units, so do not alias them onto `kd`. The `_direct` / `_delay` suffixes name the driver and follow the register's established `<stem>_<suffix>` pattern (`lcl_renal` / `lcl_nonren`, `kge_ctdna` / `kse_ctdna`). Use these names only for the *linear* drug-effect form; when a paper instead fits an Emax or sigmoid-Emax effect, the potency and shape parameters belong to the `ec50` / `lec50` / `hill` family.

## Unit spellings

The machine-readable `units` block wrote the same time
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

Generic structural models (`PK_1cmt`, `PK_2cmt`,...) are dimensionless by
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

Vocabulary for permeability-limited whole-body PBPK models, in which each tissue is resolved into residual blood cells, residual plasma, extracellular water and intracellular water (compartment suffixes `_bc` / `_plasma` / `_ew` / `_iw`; see `compartment-names.md`). Two conventions define the family. First, system physiology is carried as *fractions* -- `fvol_<organ>` of body weight and `fq_<organ>` of cardiac output -- so that a published table of percentages transcribes directly and both columns can be checked to sum to 100%. Second, every permeability-surface-area product and every clearance is expressed as a **fold-multiple of the compartment's own plasma flow** (`fold_*`), which is how these papers parameterise their what-if scenarios: one multiplier per mechanism, swept around the tissue blood flow.

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

---

## Extracorporeal circuit parameters (CRRT / ECMO / apheresis)

Prescription, geometry and transport parameters of an extracorporeal blood circuit carried as *structure* (see the `circuit_*` / `reservoir_*` compartment namespace in `inst/references/compartment-names.md`). Like that namespace these names are modality-neutral, so an ECMO oxygenator or an apheresis column reuses them. The five `q*` flows are the settings a clinician dials into the machine and are exposed as fixed `ini()` entries so a user can re-point a circuit model at a different prescription with `rxSolve(params = ...)` rather than editing the model. Registered 2026-09-02 alongside the compartment family (task `oare_PMC13274749` sidecar question q1, answer A).

### qbfr (**canonical extracorporeal circuit blood flow rate**)
- **Type:** paper-named-param
- **Role:** Blood flow rate drawn from the patient (or the ex vivo reservoir) through the circuit. Splits on the systemic haematocrit into the red-cell and plasma delivery fluxes.
- **Source aliases:** `Q_BFR`, `BFR`, `blood flow rate`
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)
- **Notes:** Papers print this in mL/min while every other circuit flow is printed in mL/h; convert both to the model's `units$time` and say so in the parameter comment.

### qdia (**canonical dialysate flow rate**)
- **Type:** paper-named-param
- **Role:** Countercurrent dialysate flow rate through a hemofilter or dialyser. Enters the effluent flow, `qeff = qdia + qfil`.
- **Source aliases:** `Q_DIA`, `DIA`, `dialysate flow rate`
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)
- **Notes:** Zero for a pure hemofiltration (CVVH) prescription; non-zero for CVVHD and CVVHDF.

### qpbp (**canonical pre-blood pump flow rate**)
- **Type:** paper-named-param
- **Role:** Pre-blood pump replacement fluid flow rate, infused into the blood line upstream of the filter. Dilutes the blood entering the filter, which is why the pre-filter haematocrit is `hct * qbfr / (qbfr + qpbp)`.
- **Source aliases:** `Q_PBP`, `PBP`, `pre-blood pump`
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)

### qrep (**canonical replacement fluid flow rate**)
- **Type:** paper-named-param
- **Role:** Total post-filter replacement fluid flow rate. Contributes to the blood-to-dialysate filtration flow `qfil = qpbp + qrep + qpfr` and is subtracted from the return-line flow.
- **Source aliases:** `Q_REP`, `REP`, `replacement fluids`
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)

### qpfr (**canonical patient fluid removal rate**)
- **Type:** paper-named-param
- **Role:** Net ultrafiltration rate removed from the patient (fluid removal prescription). Contributes to `qfil` and is subtracted from the return-line flow. Zero in a bench circuit that is not net-removing fluid.
- **Source aliases:** `Q_PFR`, `UFR`, `patient fluid removal`
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)
- **Notes:** Papers sometimes tabulate the same quantity as `UFR` (ultrafiltration rate); check whether the tabulated number is the *net* patient removal or the *gross* filtrate flow before mapping it here, because the gross flow is `qfil`, not `qpfr`.

### vreservoir (**canonical ex vivo blood reservoir volume**)
- **Type:** paper-named-param
- **Role:** Volume of the closed-loop blood reservoir of an ex vivo bench circuit; splits on the haematocrit into the `reservoir_plasma` and `reservoir_rbc` volumes.
- **Source aliases:** `reservoir volume`
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)

### vfilter (**canonical circuit blood-side volume**)
- **Type:** paper-named-param
- **Role:** Blood-side volume of the circuit's drug-contacting body -- the hemofilter (or oxygenator) together with the tubing prime volume; splits on the mean filter haematocrit into the `circuit_plasma` and `circuit_rbc` volumes.
- **Source aliases:** `filter volume`, `Filter Blood Volume`, `Circuit Blood Volume`
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)
- **Notes:** Device product guides tabulate the filter blood volume and the circuit blood volume separately; a model that lumps the tubing into the filter state should carry the *circuit* volume here and say which it used.

### vdialysate (**canonical hemofilter dialysate-side volume**)
- **Type:** paper-named-param
- **Role:** Dialysate-side volume of a hemofilter or dialyser; the volume of the `circuit_dialysate` state.
- **Source aliases:** `dialysate volume`, `filter prime volume`
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)

### safilter (**canonical hemofilter membrane surface area**)
- **Type:** paper-named-param
- **Role:** Semipermeable membrane surface area of a hemofilter, dialyser or oxygenator. Enters the diffusive membrane transfer term `dfilter * safilter / thfilter`.
- **Source aliases:** `SA`, `Membrane Surface Area`
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)
- **Notes:** Printed in m2 by device product guides; carry it in whatever length unit makes `dfilter * safilter / thfilter` come out in the model's flow units, and record the conversion in the parameter comment.

### thfilter (**canonical hemofilter membrane thickness**)
- **Type:** paper-named-param
- **Role:** Semipermeable membrane (fiber wall) thickness of a hemofilter, dialyser or oxygenator; the diffusion path length in `dfilter * safilter / thfilter`.
- **Source aliases:** `h`, `membrane thickness`, `Fiber Wall Thickness`
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)

### dfilter (**canonical membrane diffusion coefficient**)
- **Type:** paper-named-param
- **Role:** Diffusion coefficient of the drug across the hemofilter membrane (area / time). Together with `safilter` and `thfilter` it sets the diffusive permeability-surface product that drives drug from circuit plasma into the dialysate down the unbound concentration gradient.
- **Source aliases:** `D`, `MembraneDiffusionCoefficient`
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)
- **Notes:** Usually an *estimated* quantity (optimized against ex vivo circuit data) rather than a device constant, so it is drug-specific as well as filter-specific.

### krbc (**canonical red-cell to plasma concentration ratio**)
- **Type:** paper-named-param
- **Role:** Equilibrium ratio of the drug concentration in red blood cells to that in plasma; the target of the permeability-limited plasma / red-cell exchange, `flux = PS * (C_plasma - C_rbc / krbc)`. Related to the blood-to-plasma ratio by `B:P = krbc * hct + 1 - hct`.
- **Source aliases:** `K_rbc`, `Partition coefficient (blood cells/plasma)`
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)
- **Notes:** PK-Sim and similar platforms compute this internally from lipophilicity and the red-cell composition rather than exposing it as an input; when extracting from such a paper, evaluate the platform's own stored formula rather than substituting a literature blood-to-plasma ratio, and check the implied B:P against the published value.

### permrbc (**canonical plasma to red-cell membrane permeability**)
- **Type:** paper-named-param
- **Role:** Membrane permeability (length / time) governing passive diffusion of drug between plasma and red blood cells. Multiplied by the exchange surface area and the fraction unbound to give the permeability-surface product.
- **Source aliases:** `P (plasma->blood cells)`, `P (blood cells->plasma)`
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)
- **Notes:** For small molecules the resulting permeability-surface product is typically orders of magnitude larger than the blood flow, so the two phases stay at equilibrium; carrying the exchange explicitly rather than assuming instant equilibrium keeps the assumption visible and testable.

### satovrbc (**canonical red-cell exchange area per unit red-cell volume**)
- **Type:** paper-named-param
- **Role:** Effective plasma / red-cell exchange surface area per unit red-cell volume (inverse length). Multiplied by the haematocrit and the compartment volume it gives the exchange surface area for `permrbc`.
- **Source aliases:** `A_to_V_bc`, `Surface/Volume ratio (blood cells)`
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)
- **Notes:** Platform implementations fold a fixed accessibility factor into this constant (Open Systems Pharmacology multiplies the geometric surface-to-volume ratio by 0.6); record whether the registered value is the geometric ratio or the effective one.

### vmaxcrrt (**canonical maximum circuit adsorption rate**)
- **Type:** paper-named-param
- **Role:** Maximum rate (amount / time) of the saturable drug-loss process in an extracorporeal circuit -- adsorption to the membrane and tubing plus any filter clearance and degradation not modelled separately -- in `vmaxcrrt * C / (C + kmcrrt)`.
- **Source aliases:** `V_max` (McKnite 2026 Results 3.3 and MoBi reaction `MidazolamDegradation`)
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)
- **Notes:** Kept distinct from the registered `vmax` (metabolic Michaelis-Menten maximum rate) because it is a *device* property, drug- and circuit-specific rather than enzyme-specific, and because a circuit model may need both at once when the module is grafted onto a body model that metabolises the drug. Pairs with [[kmcrrt]].

### kmcrrt (**canonical circuit adsorption half-saturation concentration**)
- **Type:** paper-named-param
- **Role:** Concentration at which the saturable circuit drug-loss process runs at half its maximum rate, in `vmaxcrrt * C / (C + kmcrrt)`. The low-concentration limit `vmaxcrrt / kmcrrt` is the effective first-order circuit clearance and is the quantity to compare against the flow-limited filter clearance `fu * qeff`.
- **Source aliases:** `K_m` (McKnite 2026 Results 3.3 and MoBi reaction `MidazolamDegradation`)
- **Example models:** `McKnite_2026_midazolam_crrt.R` (founding example)
- **Notes:** Papers report this in molar units (pmol/mL, umol/L) even when the model is built in mass units; convert with the molecular weight and record the arithmetic in the parameter comment. Pairs with [[vmaxcrrt]].
