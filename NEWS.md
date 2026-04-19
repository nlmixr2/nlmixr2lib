# nlmixr2lib

# development version

* `checkModelConventions()` — new function that reports deviations from the package's parameter-naming, covariate, compartment, and metadata conventions for a single model or the entire `modeldb`. Called automatically during `buildModelDb()` so convention drift surfaces at package-build time (existing grandfathered deviations continue to build). Canonical standards (PK parameter prefixes, `eta`-prefixed IIV, `propSd`/`addSd` residual error, canonical covariate column register with aliases, canonical compartment vocabulary) are codified in an internal `.nlmixr2libConventions` list that mirrors the `extract-literature-model` skill references. Addresses issue #39.
* Add Thakre 2022 risankizumab two-compartment population PK model with first-order SC absorption in patients with active psoriatic arthritis, pooled across one phase 1, one phase 2, and two pivotal phase 3 studies (KEEPsAKE 1 / KEEPsAKE 2). Final covariates are body weight, age, serum albumin, serum creatinine, and hsCRP on clearance, and body weight on central volume. Built to `extract-literature-model` skill conventions (`etalcl + etalvc` block correlation, `etalka` diagonal, proportional residual error, canonical `WT` / `AGE` / `ALB` / `CREAT` / `hsCRP` covariate columns) with a companion validation vignette that replicates the Table 2 per-dosing-interval exposures and runs PKNCA on the third dosing interval.
* `addEta()`, `addResErr()`, `addDepot()`, `removeDepot()`, `addTransit()`, and `removeTransit()` now accept `model` as a deprecated alias for `ui` (issue #84). Passing `model = ...` emits a deprecation warning; passing both `ui` and `model` is an error.
* `addDepot()` and `addTransit()` now work correctly when `d/dt(central)` or `d/dt(depot)` appears at the beginning or end of the model block, or when transit-compartment ODEs and residual-error (`~`) specs are interleaved with assignment lines. The newly introduced helper and ODE lines are inserted immediately adjacent to the modified ODE so that the relative order of every pre-existing model line is preserved (#77, #78).
* Markov modeling creation functions including `createMarkovModel()` were added
* Add Fasanmade 2009 infliximab two-compartment population PK model in ulcerative colitis (ACT 1 / ACT 2 pooled), with albumin, ADA (ATI), and sex effects on CL and weight and sex effects on Vc. Built to the `extract-literature-model` skill conventions: structured `covariateData`, `population` metadata, canonical `etalcl` / `etalvc`, `propSd` / `addSd`, `peripheral1`, and `SEXF` / `ADA_POS` covariate columns. Includes a validation vignette with source-trace table and PKNCA check stratified by ADA status.
* Add Hu 2026 clesrovimab two-compartment population PK model for preterm and full-term infants with allometric weight scaling, postnatal age maturation function, and race effects on clearance
* Add Clegg 2024 nirsevimab two-compartment population PK model for preterm and term infants with allometric weight scaling, postmenstrual age maturation, race, season, and ADA effects
* Verified all published-literature specific-drug and mAb-consensus models against their source papers and fixed several parameter-encoding bugs that had been latent in the package since their original addition:
  - **CarlssonPetri 2021 liraglutide**: fixed categorical covariate encoding that was zeroing individual clearance for subjects not in the indexed group. `(1 - SEXF)^e_sex_cl` → `e_sex_cl^(1 - SEXF)` (previously evaluated `0^1.12 = 0` for females); `CHILD^e_age_child_cl * ADOLESCENT^e_age_adolescent_cl` → `e_age_child_cl^CHILD * e_age_adolescent_cl^ADOLESCENT` (previously evaluated `0^1.11 * 0^1.06 = 0` for adults). IIV rewritten as `omega^2 = log(1 + CV^2)` per Table 3's explicit `%CV = sqrt(exp(omega^2) - 1) * 100` footnote.
  - **Zhu 2017 lebrikizumab**: fixed IIV variance-covariance block that was storing `sqrt(variance)` (SDs) instead of variances/covariances from Table 3.
  - **Soehoel 2022 tralokinumab**: fixed IIV block (SDs → variance-covariance matrix from Table 2 footnote `IIV = sqrt(exp(omega^2) - 1)`, correlation 0.61 applied on variances not SDs); corrected body-weight exponent on V2/V3 from `0.793` to Table 2's `0.783`.
  - **Kovalenko 2020 dupilumab**: squared the five IIV values so they store variances. Paper Methods explicitly defines `omega` as "the standard deviation [SD] of between-subject variability", and nlmixr2 `etaX ~ value` stores the variance (omega^2), so SDs needed squaring.
  - **Davda 2014 mAb consensus (PK_2cmt_mAb_Davda_2014)**: all parameters verified; no changes required.
* Filled in previously-TODO `population` metadata blocks for the five verified models above with demographics from each paper's Table 1.
* Added validation vignettes for the five verified models above (`CarlssonPetri_2021_liraglutide.Rmd`, `Zhu_2017_lebrikizumab.Rmd`, `Soehoel_2022_tralokinumab.Rmd`, `Kovalenko_2020_dupilumab.Rmd`, `PK_2cmt_mAb_Davda_2014.Rmd`). Each vignette follows the `extract-literature-model` skill conventions with a population description, source trace, virtual cohort, simulation, figure replication, and PKNCA-based NCA validation.
* Retrofit Cirincione 2017 exenatide model to the `extract-literature-model` skill conventions and fix parameter encoding bugs: `ka_max` corrected from `0.0813` to paper value `12.8` /hr, `Km` rescaled to ng/mL so units are consistent with `Cc = central / vc`, and IIV variances rewritten as `log(1 + CV^2)` rather than the `log(1 + CV)` shortcut. Replaces the character-valued `DVID` covariate with `STUDY1` / `STUDY5` binary indicators and adds a companion validation vignette with PKNCA checks against the paper's Figure 5 typical values.
* Fix endogenous-model bugs surfaced during a parameter audit:
  - **Charbonneau 2021 phenylalanine**: fixed `vd_phe` typo on the `f_gut_plasma` line (the undefined symbol made `rxSolve()` fail) — should read `vd`; removed a stray `vd *` factor from the `daily_phe_intake` augmentation output so the reported value is in mg/day rather than (L/kg)·mg/day. All Table 4 parameter values verified against the authors' Zenodo companion notebook.
  - **Kim 2006 IgG**: corrected `V1` units label from `(mg/kg)` to `(mL/kg)` (paper Table 1); closed a missing `)` in the `ljmax` label; removed a redundant `igg_0 <- css` line. Parameter values verified against Table 1.
* Extended the `extract-literature-model` skill to cover endogenous, mechanistic, and turnover models: added naming conventions for Vmax/Km/baseline/rate parameters, a new `references/endogenous-validation.md` reference covering steady-state / perturbation-recovery / mass-balance / dimensional-analysis patterns, and an explicit dimensional-analysis item in `references/verification-checklist.md` (driven by the Kim units-label and Charbonneau `daily_phe_intake` bugs, which were only catchable by dimensional analysis).
* Retrofit the remaining specific-drug models and the Davda 2014 mAb 2-compartment model to the `extract-literature-model` skill conventions (no parameter-value changes). Structured `covariateData` with `description` / `units` / `type` / `reference_category` / `notes` / `source_name`, canonical `units` list (time/dosing/concentration), and `population` metadata blocks added (with TODO placeholders where not yet sourced). IIV etas renamed to `eta` + transformed-parameter name (e.g., `etacl` → `etalcl`, `iiv_lka` → `etalka`, `bsv_fpla_*` → `etalfpla_*`). Race columns renamed to the `RACE_` prefix (`BLACK` → `RACE_BLACK`, `ASIAN` → `RACE_ASIAN`, `MULTIRACIAL` → `RACE_MULTI`, `BLACK_OTH` → `RACE_BLACK_OTH`, `ASIAN_AMIND_MULTI` → `RACE_ASIAN_AMIND_MULTI`); `ADA` → `ADA_POS`; `SEXM` → `SEXF` (value inverted to keep the original effect magnitude). Touched vignettes updated to use the canonical column names. Models touched: `CarlssonPetri_2021_liraglutide`, `Clegg_2024_nirsevimab`, `Grimm_2023_gantenerumab`, `Grimm_2023_trontinemab`, `Hu_2026_clesrovimab`, `Kovalenko_2020_dupilumab`, `Kyhl_2016_nalmefene`, `PK_2cmt_mAb_Davda_2014`, `Soehoel_2022_tralokinumab`, `Xie_2019_agomelatine`, `Zhu_2017_lebrikizumab`.

# Version 0.3.2

* Add Kim 2006 model for IgG metabolism
* Add Xie 2019 agomelatine PK model
* Drop `qs` since it will be archived
* Update tumor growth inhibition models
* `addResErr()` now works with multiple-endpoint models
* Additional testing

# Version 0.3.1

* Bug fix for replacement of multiplicative expressions in `nlmixr2lib`
* phenylalanine_charbonneau_2021 had its net protein breakdown parameter corrected
* Kyhl_2016_nalmefene model was added

# Version 0.3.0

* Added ability to choose style type when modifying models.  Currently
  supported styles are: "camel" for `variablesLikeThis`, "snake" for
  `variables_like_this`, "dot" for `variables.like.this` and "blank"
  for `variableslikethis`.  This can be selected with
  `setCombineType()`.

* With the new combination style, you can change how `eta` variables
  are constructed with the `option(nlmixr2lib.etaCombineType="camel")`
  or whatever you wish it to the variable style to be.

* Added new model building framework for building models

  - **PK model building functions**

     - `addTransit()`/`removeTransit()` which were present before, but now modified and
       made a bit more robust, more closely matching literature method
       of transit compartments.

     - `addDepot()`/`removeDepot()` which were present before, but
       modified to be a bit more robust.

     - `addWeibullAbs()` which adds a Weibull absorption to a PK model

     - `convertMM()` converts linear elimination to Michaelis-Menten elimination

     - `transPK()` converts the `cl` style parameter transformations
       to various other PK transformations like `k`, `aob`, `alpha`,
       `k12`

  - **PD model building functions**

   - `addIndirectLin()` -- this adds an indirect effect model to a PK
     model that has a concentration `Cc` in the model.  This purposely
     uses a simple linear effect of `Cc*Ek` or `Cc*Ik` so it will be
     easy to parse and turn into other functional forms (like `Emax`
     or `Hill`).  If the PK model is not present it will use `Cc` as a
     covariate in a purely PD models.

   - `addIndirect()` -- this builds on `addIndirectLin()` and adds
     `Emax` or `Hill` models to a PK model. You can also set `imax=1`
     or `emax=1` to drop these parameters from being estimated in the
     model.  Additionally `hill=TRUE` will add a Hill coefficient to
     the sigmoid model.

   - `addEffectCmtLin()` -- this adds an effect compartment based on
     the `Cc` in the model.  The linear effect can be modified into
     other function forms.

   - `addDirectLin()` -- this adds a direct effect model based on the
     `Cc` in the model.

   - **Changing functional forms of Effect models**

     - `convertEmax()` changes linear effect models to Emax models

     - `convertEmaxHill()` changes linear effect models to Hill models

     - `convertQuad()` changes linear effect models to quadratic models

     - `convertLogLin()` changes linear effect models to log-linear models

   - **Changing functional forms of Baselines in non-indirect response models**

     - `addBaselineConst()` changes the zero baseline to a estimated
       constant

     - `addBaselineLin()` changes the zero baseline to a estimated
       constant and a linear constant with respect to `time`.

     - `addBaselineExp()` changes the zero baseline to a exponential
       decay with respect to time

     - `addBaseline1exp()` -- the baseline effect is changed from zero
       to to an exponential approaching to a constant (with respect to
       time).

   - **Changing model properties** (all use `addCmtProp()`)

      - `addBioavailability()` adds bioavailability property to a
        compartment

      - `addRate()` adds a modeled rate to a compartment

      - `addDur()` adds modeled duration to a compartment

      - `addIni()` adds an initial value to a compartment

      - `addLag()` adds a lag time to the a compartment

* Add Carlsson Petri (2021) liraglutide PK model
* Add Cirincione (2017) exenatide immediate-release PK model
* Add a variety of indirect response models
* Add a variety of tumor growth inhibition models and move all oncology models
  into a new model database directory
* Add a variety of double-absorption PK models
* `cp` and related `cpddSd` and `cppropSd` were renamed to `Cc`, `CcAddSd` and
  `CcPropSd` (fix #70).
* Multiple-endpoint models will have the `DV` column in the modeldb separated by
  commas.

# Version 0.2.0

* Work with the new `rxode2` version 2.0.12 `model()` and `ini()` assignment
  methods.
* Therapeutic-area specific models have begun being added.
* Models can now give the user some additional information load via the
  `message` meta-data.
* Models can now be in different directories.  The change is for ease of
  maintaining the library, it is not a change that affects users.
* A regression where `addEta()` did not change the parameter, related to a
  change in `rxode2`, was fixed.
* `addEta()` detects where to add etas more robustly when covariates are on the
  parameter.

## Models added

* Add Davda (2014) mAb consensus model
* Add Liu (2017) time-dependent clearance model based on nivolumab
* Add Kovalenko (2020) dupilumab PK model
* Add Soehoel (2022) tralokinumab PK model
* Add Zhu (2017) lebrikizumab PK model

# Version 0.1.0

* Initial version
