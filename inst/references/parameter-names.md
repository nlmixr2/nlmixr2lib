# Canonical parameter names

This file is the authoritative register of structural and paper-mechanistic parameter names used in nlmixr2lib models. Every parameter referenced inside a model's `ini()` and `model()` blocks is expected to match one of the canonical names below (or, for covariate effects, the documented `e_<cov>_<param>` shape — see `R/conventions.R::covEffectPattern`). The register is seeded from the 2026-05-28 naming audit and extended whenever a new paper introduces a parameter that isn't yet registered.

## How to use this register

1. **Before adding a parameter to a new model**, search this file (by canonical name and by source alias) for the concept you need.
2. **If the canonical name exists**, use it exactly. Document any source-paper rename in a code comment or in the model's `description` block.
3. **If the source paper uses an alias listed under an existing canonical name**, prefer the canonical name. Aliases are documented for cross-reference, not as a free pass to introduce the deprecated form in new models.
4. **If the parameter is not in this register at all**, propose a new entry with a canonical name, type, role, source aliases, and example models. Verify with the user before committing. The addition is part of the model's PR.
5. **Do not modify existing model files when you discover a missing entry**; simply register the canonical here. Retrofitting existing models is a separate effort.

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

- `log-transformed-pk` → `pkParams` (used by `.isPkParam` and the `l<base>` convention check)
- `bare-pk` → `pkBareParams` (used by `.isPkBareParam`, the bare-counterpart check, and covariate-effect shared-exponent detection)
- `paper-named-param` → `paperNamedParams` (paper-mechanistic parameters that fall outside the standard `ka`/`cl`/`vc` shape but recur across published models)

## Regex constants (kept in R, not migrated)

The following pattern constants remain hard-coded in `R/conventions.R::.nlmixr2libConventionsStatic` because they are structural regular expressions rather than name lists. They are documented here for cross-reference:

- `covEffectPattern = "^e_[A-Za-z0-9]+(_[A-Za-z0-9]+){1,5}$"` — canonical shape of a covariate-effect parameter name: a leading `e_`, a covariate token, and one to five additional tokens (final token is the affected PK parameter, optional middle tokens carry metabolite / shared-exponent / CL-component context).
- `clComponents = c("ss", "time", "renal", "nonren")` — suffix tokens permitted on multi-component CL parameters (e.g., `cl_renal + cl_nonren`).
- `transformPrefixes = c("l", "logit", "probit")` — accepted parameter-name transform prefixes.
- `residualError = c("propSd", "addSd", "expSd")` — canonical residual-error SD names. `propSd` / `addSd` are the proportional and additive SDs for the standard `~ prop(...)`, `~ add(...)`, and combined `~ prop(...) + add(...)` forms. `expSd` is the log-scale residual SD used with `~ lnorm(...)`.
- `deprecatedResidualError`, `deprecatedIivPrefixes`, `deprecatedVolumeNames`, `deprecatedVmaxNames`, `deprecatedParentSuffix` — deprecation lists also remain in R.

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

### lq (**canonical log-transformed first inter-compartmental clearance**)
- **Type:** log-transformed-pk
- **Role:** Inter-compartmental clearance between central and first peripheral compartment (volume / time).
- **Source aliases:** none.
- **Example models:** universal in 2-compartment+ popPK extractions.

### lq2 (**canonical log-transformed second inter-compartmental clearance**)
- **Type:** log-transformed-pk
- **Role:** Inter-compartmental clearance between central and second peripheral compartment (volume / time).
- **Source aliases:** none.
- **Example models:** 3-compartment popPK extractions.

### lq3 (**canonical log-transformed third inter-compartmental clearance**)
- **Type:** log-transformed-pk
- **Role:** Inter-compartmental clearance between central and third peripheral compartment (volume / time) in 4-compartment models.
- **Source aliases:** none.
- **Example models:** `Schmitt_2018_vinflunine.R`, `Li_2017_brentuximab.R`, `Weatherley_2009_maraviroc_iv.R`.
- **Notes:** Pairs with `peripheral3` and `lvp3`.

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

### ltmax_abs (**canonical log-transformed saturable absorption Vmax**)
- **Type:** log-transformed-pk
- **Role:** Log-scale maximum rate of saturable Michaelis-Menten ABSORPTION from the depot / absorption compartment into central (amount / time). Distinct from `lvmax`, which is MM ELIMINATION from central.
- **Source aliases:**
  - `Tmax` -- used in `Jansson_2008_eflornithine_rat.R` for the saturable-absorption Vmax. The paper symbol Tmax is not a time-to-Cmax descriptor but the maximum absorption rate.
- **Example models:** `Jansson_2008_eflornithine_rat.R`.
- **Notes:** Paired with `lkt_abs` (half-saturation amount). The `_abs` suffix disambiguates from `lvmax` (MM elimination).

### lkt_abs (**canonical log-transformed saturable absorption half-saturation amount**)
- **Type:** log-transformed-pk
- **Role:** Log-scale amount of drug in the depot / absorption compartment at which the saturable absorption rate equals half its maximum (amount). Functionally a Km but expressed in amount-in-compartment units rather than concentration.
- **Source aliases:**
  - `Kt` -- used in `Jansson_2008_eflornithine_rat.R`.
- **Example models:** `Jansson_2008_eflornithine_rat.R`.
- **Notes:** Paired with `ltmax_abs`. The amount-in-compartment formulation differs from a concentration Km because the absorption compartment may not have a well-defined volume.

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

### lkel (**canonical log-transformed elimination rate constant (K-PD)**)
- **Type:** log-transformed-pk
- **Role:** First-order elimination rate constant used when no explicit `vc` is estimated (K-PD or single-rate-constant elimination form).
- **Source aliases:**
  - `lke` -- legacy name; replaced 2026-05-28 by the naming audit.
  - `lkde` -- paper-named (Mazzocco 2015 / Shoji 2017 KDE) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
  - `lkp` -- paper-named (van Hasselt 2015 KP) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
- **Example models:** `Mazzocco_2015_temozolomide.R`, `Shoji_2017_fosdagrocorat_oc.R`, `Shoji_2017_fosdagrocorat_p1np.R`, `vanHasselt_2015_eribulin.R`.
- **Notes:** Canonical `lkel` adopted 2026-05-28 per the naming audit.

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
- **Role:** Inter-compartmental clearance between central and `peripheral2` (volume / time).
- **Source aliases:** none.
- **Example models:** universal in 3-compartment popPK extractions.

### q3 (**canonical bare third inter-compartmental clearance**)
- **Type:** bare-pk
- **Role:** Inter-compartmental clearance between central and `peripheral3` (volume / time) in 4-compartment popPK models.
- **Source aliases:** none.
- **Example models:** `Schmitt_2018_vinflunine.R`, `Li_2017_brentuximab.R`, `Weatherley_2009_maraviroc_iv.R`.

### kel (**canonical bare elimination rate constant (K-PD)**)
- **Type:** bare-pk
- **Role:** First-order elimination rate constant in K-PD / single-rate-constant elimination models with no explicit `vc`.
- **Source aliases:**
  - `ke` -- legacy.
  - `kde` -- paper-named (Mazzocco 2015 / Shoji 2017 / Xia 2024 KDE) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
  - `kp` -- paper-named (van Hasselt 2015 KP) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
  - `ps_elim`, `pc_elim` -- paper-named (Wilson 2015 p_S / p_C) bare drug-specific K-PD elim rates; replaced 2026-05-30 by `kel_sunitinib` / `kel_irinotecan`.
- **Example models:** `Mazzocco_2015_temozolomide.R`, `Shoji_2017_fosdagrocorat_oc.R`, `Shoji_2017_fosdagrocorat_p1np.R`, `vanHasselt_2015_eribulin.R`, `Wilson_2015_sunitinib_irinotecan_mouse.R` (bare drug-suffixed `kel_<drug>`), `Xia_2024_warfarin.R`.

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

---

## Paper-named mechanistic parameters

Parameters that don't fit the standard `ka` / `cl` / `vc` shape but recur across published models. Each entry is treated as a canonical bare name; the log-transformed form (`l<name>`) is acceptable wherever the parameter is strictly positive and the source paper reports an exponential typical-value form. Add to this list rather than introducing a new ad-hoc pattern.

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
- **Example models:** `Kleijn_2011_sugammadex_rocuronium.R` (`lke0 = log(0.134) = -2.01`, giving `ke0 = 0.134 1/min` for the rocuronium effect-compartment equilibration; allometric scaling `(BW/70)^-0.25`), `Savic_2010_warfarin.R` (founding registry example -- 1-compartment PK with effect-compartment-driven proportional-odds PD), `Schindler_2016_sunitinib.R` (effect-compartment equilibration with daily AUC at ln(2)/50 1/h).
- **Notes:** Distinct from `lke_kpd` (which is K-PD elimination rate and was deprecated in favour of the canonical `lkel`). `lke0` is specifically the effect-compartment equilibration parameter for hysteresis PK-PD modelling.

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
- **Source aliases:** none.
- **Example models:** `Kleijn_2011_sugammadex_rocuronium.R`.

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
- **Role:** Proliferation / growth-rate constant (1 / time). Paper-mechanistic meaning depends on the model.
- **Source aliases:** none.
- **Example models:** TGI / cell-population PD models.

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
- **Notes:** First-order transit rate constant `ktr = n_transit / mtt`.

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

### kin (**canonical indirect-response production rate**)
- **Type:** paper-named-param
- **Role:** Zero-order production rate of an indirect-response / turnover pool (Dayneka 1993; Jusko traditions).
- **Source aliases:** none.
- **Example models:** indirect-response PD models.
- **Notes:** Codified 2026-05-28 per the naming audit.

### kout (**canonical indirect-response elimination rate**)
- **Type:** paper-named-param
- **Role:** First-order elimination rate of an indirect-response turnover pool (1 / time).
- **Source aliases:** none.
- **Example models:** indirect-response PD models.

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
- **Notes:** Usually held fixed at the in-vitro equilibrium-dialysis-derived value. The BBB transfer term is `CLin * fu * Cp`.

### eh (**canonical hepatic extraction ratio**)
- **Type:** paper-named-param
- **Role:** Hepatic extraction ratio in the well-stirred liver model (`CL_H = FQ * eh`, `F_HEP = 1 - eh`). Unitless, bounded in [0, 1]. Use whenever a paper parameterises hepatic clearance through the extraction-ratio physiology rather than through `lcl` / `lcl_nonren` directly. The companion log-transform prefix is `logiteh` (logit-scale `ini()` typical value) because the linear-scale `eh` is bounded; logit-additive eta on `logiteh` keeps every individual `eh_i` inside the (0, 1) box.
- **Source aliases:** `EH` -- used in `Chan_2008_maraviroc.R` and Brussee 2018 mAb PBPK (the latter as a derived `model()`-block quantity rather than an `ini()` parameter, so the canonical applies to the `ini()` use case introduced by Chan 2008).
- **Example models:** `Chan_2008_maraviroc.R` (estimated `logiteh = logit(0.662)` with logit-additive IIV per Chan 2008 Eq 9; downstream `eh` enters `clh = fq * eh` and `fhep = 1 - eh`).
- **Notes:** Paired with the canonical compartment / pseudo-parameter `fq` (hepatic plasma flow, fixed at a literature value) and the canonical `clr` (renal clearance, often fixed) when the paper decomposes total CL into renal + hepatic with hepatic-extraction physiology. Distinct from `lcl_nonren` (additive renal + non-renal CL decomposition without an explicit extraction-ratio bound): use `eh` only when the source paper writes hepatic clearance as `CL_H = FQ * E_H` and constrains E_H in [0, 1] (e.g., via a logit-form IIV transformation).
