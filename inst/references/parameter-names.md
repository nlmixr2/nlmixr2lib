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
- **Example models:** `Jonckheere_2019_cefepime.R` and similar renally cleared small-molecule popPK extractions.

### lcl_nonren (**canonical log-transformed non-renal clearance arm**)
- **Type:** log-transformed-pk
- **Role:** Non-renal (hepatic / metabolic / extra-renal) component of an additive renal + non-renal clearance decomposition.
- **Source aliases:** none.
- **Example models:** `Jonckheere_2019_cefepime.R` and similar.

### lkel (**canonical log-transformed elimination rate constant (K-PD)**)
- **Type:** log-transformed-pk
- **Role:** First-order elimination rate constant used when no explicit `vc` is estimated (K-PD or single-rate-constant elimination form).
- **Source aliases:**
  - `lke` -- legacy name; replaced 2026-05-28 by the naming audit.
  - `lkde` -- paper-named (Mazzocco 2015 / Shoji 2017 KDE) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
  - `lkp` -- paper-named (van Hasselt 2015 KP) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
- **Example models:** `Mazzocco_2015_temozolomide.R`, `Shoji_2017_fosdagrocorat_oc.R`, `Shoji_2017_fosdagrocorat_p1np.R`, `vanHasselt_2015_eribulin.R`, K-PD templates. The drug-suffixed `lkel_<drug>` (Wilson 2015: `lkel_sunitinib`, `lkel_irinotecan`) carries the per-drug K-PD elimination rate in combination K-PD models.
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

### kel (**canonical bare elimination rate constant (K-PD)**)
- **Type:** bare-pk
- **Role:** First-order elimination rate constant in K-PD / single-rate-constant elimination models with no explicit `vc`.
- **Source aliases:**
  - `ke` -- legacy.
  - `kde` -- paper-named (Mazzocco 2015 / Shoji 2017 / Xia 2024 KDE) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
  - `kp` -- paper-named (van Hasselt 2015 KP) form; replaced 2026-05-30 by the K-PD canonical-name retrofit.
  - `ps_elim`, `pc_elim` -- paper-named (Wilson 2015 p_S / p_C) bare drug-specific K-PD elim rates; replaced 2026-05-30 by `kel_sunitinib` / `kel_irinotecan`.
- **Example models:** K-PD templates plus `Mazzocco_2015_temozolomide.R`, `Shoji_2017_fosdagrocorat_oc.R`, `Shoji_2017_fosdagrocorat_p1np.R`, `vanHasselt_2015_eribulin.R`, `Wilson_2015_sunitinib_irinotecan_mouse.R` (bare drug-suffixed `kel_<drug>`), `Xia_2024_warfarin.R`.

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
- **Example models:** `Svensson_2016_bedaquiline.R`, transit-absorption popPK extractions.

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
