# Canonical covariate columns

This file is the authoritative register of covariate column names used in nlmixr2lib models. Every covariate referenced inside a model's `model()` block must have an entry here. The register is seeded from a full audit of `inst/modeldb/` and extended whenever a new paper introduces a covariate that isn't yet registered.

## How to use this register

1. **Before adding a covariate to a new model**, search this file (by canonical name and by source alias) for the concept you need.
2. **If the canonical name exists**, use it exactly. Document any source-paper rename in the model's `covariateData[[name]]$source_name` field.
3. **If the source paper uses an alias listed under an existing canonical name**, use the canonical name and note in `covariateData[[name]]$notes` whether the values must be transformed (e.g., `SEXM -> SEXF` inverts values; the effect coefficient sign / reference category must be inverted as well).
4. **If the covariate is not in this register at all**, propose a new entry with a canonical name, description, units, type, scope, and source aliases. Verify with the user before committing. The addition is part of the model's PR.
5. **Do not modify existing model files when you discover an alias**; simply document the mapping in the register. Retrofitting existing models is a separate effort.

## Scope: general vs specific

Each entry has a `Scope:` field declaring whether it is **general** (any model may use it without warning) or **specific** (only the models listed under `Example models` may use it; other usage is flagged by `checkModelConventions()`). This prevents accidentally reusing a covariate name whose meaning is tied to a particular paper.

- **Scope: general** -- the covariate has a stable, paper-independent meaning and any future model may use it. Examples: `WT`, `AGE`, `SEXF`, `CREAT`, `ADA_POS`, `CRP`, `CRCL`.
- **Scope: specific** -- the covariate's semantics depend on a particular study's design (a specific study indicator, a drug-product variant, a composite race grouping, a tumor-type decomposition, etc.). If a new paper needs the same concept with different semantics, register a new canonical name; if the concept matches, extend the `Example models` list (and consider promoting to general).

When adding or updating an entry, choose the most conservative scope: if in doubt, start with `specific` and promote when a second model legitimately ratifies the name.

## Case convention

Covariate column names should be ALL CAPS. Current non-all-caps canonical names are `dilution` and `nonECZTRA` (both scope: specific), preserved from their source files with "future rename" notes. New entries should default to all caps.

## Entry schema

```yaml
- name: <CANONICAL_NAME>
  description: <one-sentence definition>
  units: <unit string, or "(binary)" / "(categorical)">
  type: continuous | binary | categorical | count
  scope: general | specific
  reference_category: <the 0 group for binary/categorical, or NULL>
  source_aliases:
    - <ALIAS_NAME> (<transformation if any>) -- used in <model.R>
  example_models:
    - <model.R>
  notes: <free text>
```

---

## Demographics

### WT (**canonical for body weight (baseline or time-varying)**)
- **Description:** Body weight (baseline or time-varying).
- **Units:** kg
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with allometric scaling `(WT / ref_wt)^exponent`. Reference weights observed: 70 kg (adults), 75 kg, 84.8 kg, 5 kg (infants).
- **Source aliases:** none known.
- **Example models:** `Clegg_2024_nirsevimab.R`, `Hu_2026_clesrovimab.R`, `Zhu_2017_lebrikizumab.R`, `Kovalenko_2020_dupilumab.R`, `CarlssonPetri_2021_liraglutide.R`, `Cirincione_2017_exenatide.R`, `Grimm_2023_gantenerumab.R`, `Grimm_2023_trontinemab.R`, `Kyhl_2016_nalmefene.R`, `Soehoel_2022_tralokinumab.R`, `Xie_2019_agomelatine.R`, `PK_2cmt_mAb_Davda_2014.R`, `phenylalanine_charbonneau_2021.R`, `Chua_2025_mirikizumab.R`, `Jackson_2022_ixekizumab.R`, `Kotani_2022_astegolimab.R`, `Ma_2020_sarilumab_anc.R`, `Ma_2020_sarilumab_das28crp.R`, `Moein_2022_etrolizumab.R`, `Tiraboschi_2025_amlitelimab.R`, `Robbie_2012_palivizumab.R`, `Bajaj_2017_nivolumab.R`, `Quartino_2019_trastuzumab.R`, `Wang_2020_ontamalimab.R`, `Fau_2020_isatuximab.R`, `Okada_2025_rocatinlimab.R`, `Kunisawa_2014_olprinone.R`, `Xu_2020_daratumumab.R` (reference 78.6 kg; power exponents 0.451 on linear CL and 0.375 on V1).
- **Notes:** Universal. Verify time-varying vs. baseline-only against the source paper.

### AGE (**canonical for subject age**)
- **Description:** Subject age in years.
- **Units:** years
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a
- **Source aliases:** none.
- **Example models:** `Archary_2019_lamivudine.R`, `Budha_2023_tislelizumab.R`, `Chakraborty_2012_canakinumab.R`, `Chen_2020_luspatercept.R`, `Conrado_2014_alzheimer.R`, `Diepstraten_2013_propofol.R`, `Gandhi_2021_abatacept.R`, `Goel_2016_Sonidegib.R`, `Hennig_2013_tobra.R`, `Hong_2025_datopotamab.R`, `Ide_2020_elotuzumab.R`, `Koopman_2023_factorix.R`, `Kuchimanchi_2024_dostarlimab.R`, `Kunarajah_2017_doxorubicin.R`, `Kyhl_2016_nalmefene.R`, `Lahu_2010_roflumilast.R`, `Li_2006_meropenem.R`, `Li_2017_cediranib.R`, `Li_2019_abatacept.R`, `Lin_2024_casirivimab.R`, `Martinez_2019_alirocumab.R`, `Melhem_2022_dostarlimab.R`, `Mulyukov_2018_ranibizumab.R`, `NA_NA_tte_gompertz.R`, `NA_NA_tte_lognormal.R`, `Retlich_2015_linagliptin.R`, `Rosario_2015_vedolizumab.R`, `Svensson_2016_bedaquiline.R`, `Thakre_2022_risankizumab.R`, `Wu_2024_inotuzumab.R`, `Yassen_2025_asundexian.R`, `Yu_2022_ofatumumab.R`, `Zhong_2026_abatacept.R`, `Zhou_2021_belimumab.R`, `Zhu_2017_lebrikizumab.R`.
- **Notes:** Zhu 2017 normalizes as `AGE/40`.

### LBM (**canonical for lean body mass**)
- **Description:** Lean body mass.
- **Units:** kg
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a
- **Source aliases:**
  - `LBW` (lean body weight) -- synonym; same biological quantity (total body weight minus body fat). Hemophilia popPK literature typically uses `LBW` (Hume or James formula) where mAb / general literature uses `LBM`. Used in `Garmann_2017_BAY81_8973.R` (reference 51.1 kg).
- **Example models:** `Kyhl_2016_nalmefene.R` (reference 56.28 kg, exponent 0.626 on CL), `Garmann_2017_BAY81_8973.R` (alias `LBW`; reference 51.1 kg, exponents 0.610 on CL and 0.950 on Vc).

### FFM (**canonical for fat-free mass**)
- **Description:** Fat-free mass derived from body weight, height, and sex via the Janmahasatian et al. formula (Clin Pharmacokinet 2005;44:1051-1065).
- **Units:** kg
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(FFM / ref)^exponent`. Reference values observed: 40.69 kg (Zhou 2021 belimumab pooled adult+pediatric SLE), 45 kg (Aguiar 2021, Crohn's disease cohort median).
- **Source aliases:** none; `FFM` is the universal abbreviation.
- **Example models:** `Zhou_2021_belimumab.R` (reference 40.69 kg; exponents 0.673 on CL and 0.891 on V1), `Aguiar_2021_ustekinumab.R` (reference 45 kg; power exponents 0.598 on CL, 0.590 on Vc, 0.586 on Vp).
- **Notes:** Distinct from `LBM` (lean body mass) which is sometimes computed by the Boer or Hume formulae. When the source paper reports the body-composition formula it used (e.g., Janmahasatian for FFM), record it in `covariateData[[FFM]]$notes`. FFM is preferred over total body weight when scaling monoclonal-antibody PK because mAb distribution is largely confined to extracellular fluid; muscle / lean tissue tracks extracellular volume better than total weight in heavier patients.

### BSA (**canonical for body surface area**)
- **Description:** Body surface area (typically computed by DuBois, Mosteller, or Haycock from height and weight).
- **Units:** m^2
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(BSA / ref)^exponent`.
- **Source aliases:** none.
- **Example models:** `Yamada_2025_zolbetuximab.R` (reference 1.70 m^2; exponents 1.06 on clearances and 0.968 on volumes).
- **Notes:** Oncology mAbs dosed by BSA (mg/m^2) often use BSA in place of body weight for allometric-style scaling. Document the BSA computation formula (DuBois / Mosteller / Haycock) the source paper used; if unstated, record "unspecified."

### BMI (**canonical for body mass index**)
- **Description:** Body mass index at baseline.
- **Units:** kg/m^2
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with a linear-deviation form (`1 + e * (BMI - ref)`) or a power form (`(BMI / ref)^e`). Document the reference value in `covariateData[[BMI]]$notes`.
- **Source aliases:** none known.
- **Example models:** `Chua_2025_mirikizumab.R` (reference 24.75 kg/m^2; linear-deviation effect on logit of bioavailability), `NA_NA_lidocaine.R` (DDMODEL00000281; binary stratification at threshold 27.93 kg/m^2 adding +0.939 to the GX rate constant K30 in the BMI > 27.93 cohort).
- **Notes:** Universal clinical-trial demographic. Derived as `WT / (height_m)^2`; assume time-fixed at baseline unless the source paper states otherwise.

### BMIZ (**canonical for body-mass-index z-score (age- and sex-standardised)**)
- **Description:** Age- and sex-standardised body-mass-index z-score (number of standard deviations above or below the reference-population mean BMI for the subject's age and sex). Distinct from raw `BMI` (kg/m^2): `BMIZ` is unitless and centred at 0 in the reference population, so the reference value used in linear-deviation effects is 0 (not a population BMI in kg/m^2). Time-varying when the source paper carries a per-visit z-score; document baseline-vs-time-varying status in `covariateData[[BMIZ]]$notes`.
- **Units:** unitless (z-score; standard-deviation units)
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with a linear-deviation form `(1 + e * (BMIZ - 0))` so the reference is 0 (population mean for the subject's age/sex). Effect coefficients are interpreted as fractional change per 1 z-score-unit deviation from the reference.
- **Source aliases:**
  - `BMI` -- when a paper uses the column name `BMI` for what is actually a z-score (e.g., the Harun 2019 NMTRAN control stream column `BMI` is documented in the dataset header as "body-mass index z-score"). The canonical column is `BMIZ`; the source-paper column name is recorded in `covariateData[[BMIZ]]$source_name`.
- **Example models:** `Harun_2019_cysticFibrosis.R` (time-varying per-visit BMI z-score; linear-deviation effect on baseline FEV1% predicted with reference 0 and coefficient +0.0382 per z-score unit).
- **Notes:** Distinct from `BMI` (raw kg/m^2 used in adult populations). Paediatric and adolescent studies routinely report BMI as a z-score relative to a growth reference (WHO 2007 Growth Reference for school-aged children, CDC 2000, etc.); document the reference standard the source paper used in `covariateData[[BMIZ]]$notes`. Specific scope until a second paediatric model ratifies the name; at that point promote to `general`.

### SEXF (**canonical for sex**)
- **Description:** Biological sex indicator, 1 = female, 0 = male.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (male).
- **Source aliases:**
  - `SEXM` (values inverted: `SEXF = 1 - SEXM`; effect coefficient sign and reference category both invert) -- used in `CarlssonPetri_2021_liraglutide.R`.
  - `SEX` with `"M"`/`"F"` strings -- derive `SEXF = as.integer(SEX == "F")`.
  - `SEX` with `1`=male / `2`=female numeric coding -- derive `SEXF = as.integer(SEX == 2)`. Used in `Netterberg_2017_docetaxel.R` and `NA_NA_miridesap.R` (DDMODEL00000262 source bundle; Sahota 2015 NONMEM convention).
  - `FEM` (1 = female, 0 = male; same orientation as canonical, no transformation) -- used in `Guiastrennec_2016_gastric_emptying.R`.
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical), `CarlssonPetri_2021_liraglutide.R` (alias `SEXM`), `Bajaj_2017_nivolumab.R` (male-indicator source; effect applied as `exp(coef * (1 - SEXF))` to preserve the paper's female-reference CL_REF / VC_REF), `Fau_2020_isatuximab.R` (exponential effect on Vc; reference category 0 = male), `Netterberg_2017_docetaxel.R` (multiplicative effect on baseline ANC: `BACOV *= (1 + theta * SEXF)`; source column `SEX` with 1 = male, 2 = female encoding, decomposed via `SEXF = as.integer(SEX == 2)`), `NA_NA_miridesap.R` (DDMODEL00000262 / Sahota 2015; multiplicative effect on baseline SAP via `SAP_BASE_ref * (1 + e_sexf_sap0 * SEXF)` with `e_sexf_sap0 = -0.30`; female baseline is ~30% lower than male), `Xu_2020_daratumumab.R` (additive shift on V1 `(1 + e_sexf_vc * SEXF)` with `e_sexf_vc = -0.205`: female V1 is 20.5% lower than male, reference category 0 = male), `Guiastrennec_2016_gastric_emptying.R` (multiplicative +40.7% strengthening of the caloric-feedback slope SLPCAL on gastric emptying in females; `SLPCAL_eff = SLPCAL * (1 + 0.407 * SEXF)`).
- **Notes:** When translating a model that used `SEXM`, flag the sign/reference-category inversion to the user.

### PREG (**canonical for pregnancy status indicator**)
- **Description:** 1 = pregnant, 0 = non-pregnant. Time-fixed per subject in trial cohorts that enrol pregnant and non-pregnant women in parallel; not a time-varying flag.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (non-pregnant).
- **Source aliases:** none known; source NONMEM control streams typically use `PREG` directly.
- **Example models:** `Birgersson_2019_artesunate.R` (multiplicative effect on dihydroartemisinin clearance; the published structural CLM = 190 L/h is reported with the source-paper reference category PREG = 1, so the model file applies the effect via `(1 + e_preg_cl_dha * (1 - PREG))` with `e_preg_cl_dha = -0.214` to preserve verbatim source values; non-pregnant women have ~21% lower CLM relative to pregnant women).
- **Notes:** Use this canonical for adult clinical-trial models that test a pregnancy-vs-non-pregnancy contrast (typical settings: malaria-in-pregnancy PK, antiviral-in-pregnancy PK). Trimester or gestational-age stratification within the pregnant cohort should use a separate canonical (e.g., gestational-age weeks via `GA` or a trimester indicator, ratified separately when needed). The canonical convention is reference category 0 (non-pregnant) following the broader pharmacology default; source papers that use the pregnant cohort as the reference (Birgersson 2019) preserve their published structural values via a `(1 - PREG)` form on the effect coefficient. Ratified canonically on 2026-05-07.

### CHILD (**canonical for child age-cohort indicator**)
- **Description:** 1 = subject is a child, 0 = not a child.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (not child, i.e., adult baseline).
- **Source aliases:**
  - `PED` -- used in the Schoemaker 2018 LEV / BRV pediatric extrapolation (DDMODEL00000239) as the pediatric-vs-adult indicator that gates the Markov-amplitude term, the overdispersion IIV, and the four pediatric offsets on log baseline rate / mixture / placebo / Emax / EC50.
- **Example models:** `CarlssonPetri_2021_liraglutide.R`, `Schoemaker_2018_levetiracetam.R` (DDMODEL00000239).
- **Notes:** Age-group indicator used alongside `ADOLESCENT`; paper's age cutoffs must be captured in `covariateData[[CHILD]]$notes`.

### ADOLESCENT (**canonical for adolescent age-cohort indicator**)
- **Description:** 1 = subject is an adolescent, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0.
- **Example models:** `CarlssonPetri_2021_liraglutide.R`.
- **Notes:** Paired with `CHILD`. Document age cutoffs.

## Pediatric / maturation

### PAGE (**canonical for postmenstrual age**)
- **Description:** Postmenstrual age in months (`GA_weeks / 4.35 + postnatal_months`). Time-varying.
- **Units:** months
- **Type:** continuous
- **Scope:** general
- **Example models:** `Clegg_2024_nirsevimab.R`, `Robbie_2012_palivizumab.R`.

### PNA (**canonical for postnatal age**)
- **Description:** Postnatal age (chronological since birth). Time-varying.
- **Units:** months
- **Type:** continuous
- **Scope:** general
- **Example models:** `Hu_2026_clesrovimab.R`.

### GA (**canonical for gestational age at birth**)
- **Description:** Gestational age at birth. Time-fixed per subject.
- **Units:** weeks
- **Type:** continuous
- **Scope:** general
- **Example models:** `Hu_2026_clesrovimab.R`, `Clegg_2024_nirsevimab.R` (folded into PAGE).

### WT_BIRTH (**canonical for birth weight**)
- **Description:** Body weight at birth. Time-fixed per subject. Distinct from current body weight (`WT`), which is time-varying.
- **Units:** kg
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with linear-deviation `(1 + e * (WT_BIRTH - ref))` or power scaling. Reference value observed: 2.59 kg (Voller 2017 newborn cohort).
- **Source aliases:**
  - `BWEIGHT` -- used in `Voller_2017_phenobarbital.R` (Voller 2017 source data column for birth weight in kg).
- **Example models:** `Voller_2017_phenobarbital.R` (linear-deviation effect on CL: `clbw = 1 + 0.369 * (WT_BIRTH - 2.59)`).
- **Notes:** Time-fixed at birth; characterises pre-/term-newborn cohorts. Pairs with `GA` (gestational age at birth) when both are reported. The conventional clinical-PK abbreviation `BWT` is intentionally NOT used as the canonical name because it is already used across the codebase (Gandhi 2021, Li 2019, Chen 2022, Wojciechowski 2022, Lu 2019) as a source-name alias for body weight (`WT`). The `WT_BIRTH` form keeps the `WT` root consistent with the existing body-weight canonical and avoids the `BWT` ambiguity.

### AGE_DPF (**canonical for zebrafish-larval age in days post-fertilization**)
- **Description:** Age of a zebrafish (Danio rerio) larva in days post-fertilization (dpf). Time-fixed per subject in destructive-sampling designs (each larva is harvested at exactly one observation, so its dpf is fixed at the value assigned at exposure-start).
- **Units:** days post-fertilization (dpf)
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used in van Wijk 2019 with a step-form effect on the absorption rate (`IF (AGE_DPF > 3) k12 = k12_3 * (1 + e_age_dpf_k12)`) and a per-day power-form effect on the elimination rate (`k25 = k25_3 * (1 + e_age_dpf_k25)^(AGE_DPF - 3)`). Reference age is 3 dpf (the youngest cohort in the study).
- **Source aliases:**
  - `AGE` -- van Wijk 2019 NONMEM column. Renamed to canonical `AGE_DPF` because the human-PK canonical `AGE` denotes subject age in years; the zebrafish dpf semantic is incompatible and would silently corrupt any future model that mixed them. Same orientation, no value transformation.
- **Example models:** `vanWijk_2019_paracetamol.R`.
- **Notes:** Specific scope because the meaning is bounded to zebrafish-larval-development PK studies. Distinct from canonical `AGE` (human age in years), `PNA` (postnatal age in months), `PAGE` (postmenstrual age in months), and `GA` (gestational age in weeks) -- none of those are appropriate for a non-mammalian organism whose developmental clock is anchored at fertilization rather than birth. Integer values 3, 4, 5 in the van Wijk 2019 dataset, but treated as a continuous covariate in the elimination-rate power form. Future zebrafish-PK or other non-mammalian-developmental-age models should reuse this canonical only when the covariate is indeed dpf-anchored; other developmental-time conventions (e.g., somite-stage, hpf, dph) would warrant separate canonicals. Ratified canonically on 2026-05-07.

## Nutritional status

### MAL_NOURISH (**canonical for malnutrition status indicator**)
- **Description:** 1 = subject is malnourished at study entry per a paper-defined anthropometric criterion (e.g., WHO height-for-age and weight-for-age Z-scores both < -2.0 in Tikiso 2021); 0 = not malnourished. Subject-level baseline indicator; the time-decaying recovery during nutritional supplementation is carried by the paired `T_NUT_SUPP` column rather than as a time-varying value of `MAL_NOURISH`.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (not malnourished).
- **Source aliases:**
  - `MAL` -- used in `Tikiso_2021_abacavir.R` (the dataset's paper-defined indicator, 1 = malnourished, 0 = not malnourished).
- **Example models:** `Tikiso_2021_abacavir.R` (gates the time-decaying malnutrition effect: `mal_decay = MAL_NOURISH * exp(-T_NUT_SUPP * log(2) / 12.2)`, which then drives multiplicative shifts of `+115%` on F and `-64%` on CL at the start of nutritional supplementation, decaying with a 12.2-day half-life).
- **Notes:** Specific scope because the malnutrition definition (WHO Z-score thresholds, mid-upper arm circumference, weight-for-height vs height-for-age, etc.) is paper-defined; per-model `covariateData[[MAL_NOURISH]]$notes` must document the criterion used. Pairs with `T_NUT_SUPP` (days on nutritional supplementation) when the model uses a time-decaying recovery function; otherwise `MAL_NOURISH` alone serves as a static baseline indicator. Distinct from generic body-weight Z-scores (which are continuous anthropometric metrics rather than a binarised malnutrition indicator).

### T_NUT_SUPP (**canonical for time on nutritional supplementation**)
- **Description:** Time elapsed since the start of nutritional / refeeding supplementation, in days. 0 at the start of supplementation; increases with time. Used by population PK models that describe the time-varying recovery of PK parameters during nutritional rehabilitation in malnourished cohorts.
- **Units:** days
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- enters as the argument of an exponential decay function `exp(-T_NUT_SUPP * log(2) / T_half)` whose half-life `T_half` is an estimated parameter (12.2 days in Tikiso 2021).
- **Source aliases:**
  - `TNUTRI` -- used in `Tikiso_2021_abacavir.R` (the dataset's paper-defined column for days since start of nutritional supplementation; same orientation as the canonical, 0 = start of supplementation, increasing with time on supplementation).
- **Example models:** `Tikiso_2021_abacavir.R` (paired with `MAL_NOURISH`; drives the recovery decay `exp(-T_NUT_SUPP * log(2) / 12.2)` of the malnutrition effect on F and CL).
- **Notes:** Specific scope because the underlying recovery dynamics are tied to the nutritional-rehabilitation protocol of the source study. For non-malnourished subjects (`MAL_NOURISH == 0`) the value is irrelevant because the malnutrition effect is gated by `MAL_NOURISH`; supply 0 as a default. For fully-recovered malnourished subjects, supply a large value (e.g., `>= 100` days, well beyond the 12.2-day Tikiso 2021 half-life) so the decay function reaches near zero and the effect vanishes.

### DRINK_OGTT (**canonical for oral-glucose-tolerance-test-only drink indicator**)
- **Description:** Binary indicator that the postprandial test drink is glucose-only (oral glucose tolerance test, OGTT). 1 = the drink is glucose-only with no fat content (e.g., 25 / 75 / 125 g OGTT); 0 = otherwise (water, or any drink containing fat). Per-occasion (per-test-drink-administration) covariate -- a single subject in a crossover challenge protocol receives different drink types across occasions, so the indicator varies per dose-record rather than per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (water or fat-containing drink). Mutually exclusive with `DRINK_FAT` -- both indicators cannot be 1 simultaneously.
- **Source aliases:** paper narrative "glucose solution" / "OGTT" cohort labels driving the gastric-emptying-onset T50OGTT selection in Guiastrennec 2016.
- **Example models:** `Guiastrennec_2016_gastric_emptying.R` (selects the OGTT-specific half-onset time T50OGTT = 15.7 min for the gastric-emptying delay Hill function; water is recovered when DRINK_OGTT = DRINK_FAT = 0 with the onset factor pinned to 1).
- **Notes:** Specific scope because the OGTT / fat-containing partition is tied to the Guiastrennec 2016 postprandial-challenge design (Studies B and C OGTT arms). Set to 1 for Study B's 25 / 75 / 125 g OGTT drinks and the Study C 75 g OGTT arm; set to 0 for Study A water and for all fat-containing drinks. Pairs with `DRINK_FAT`: the two indicators jointly select the appropriate gastric-emptying-delay T50 parameter (T50OGTT vs T50Fat) for the Hill onset function. Ratified canonically alongside the Guiastrennec 2016 gastric-emptying / CCK / GBE extraction.

### DRINK_FAT (**canonical for fat-containing test-drink indicator**)
- **Description:** Binary indicator that the postprandial test drink contains fat (any nonzero fat content). 1 = the drink contains fat (e.g., the Study C low / medium / high-fat isocaloric drinks, or the Study D medium-high-fat drink); 0 = otherwise (water or glucose-only OGTT drinks). Per-occasion covariate -- a single subject in a crossover challenge protocol receives different drink types across occasions.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (water or glucose-only drink). Mutually exclusive with `DRINK_OGTT` -- both indicators cannot be 1 simultaneously.
- **Source aliases:** paper narrative "low-fat" / "medium-fat" / "high-fat" / "medium-high-fat" cohort labels driving the gastric-emptying-onset T50Fat selection in Guiastrennec 2016.
- **Example models:** `Guiastrennec_2016_gastric_emptying.R` (selects the fat-specific half-onset time T50Fat = 23.1 min for the gastric-emptying delay Hill function).
- **Notes:** Specific scope because the fat-containing partition is tied to the Guiastrennec 2016 postprandial-challenge design. Pairs with `DRINK_OGTT`: the two indicators jointly select the appropriate gastric-emptying-delay T50 parameter (T50OGTT vs T50Fat) for the Hill onset function. Ratified canonically alongside the Guiastrennec 2016 gastric-emptying / CCK / GBE extraction.

## Pregnancy / hormonal status

### TERM_BIRTH (**canonical for term-vs-preterm birth indicator**)
- **Description:** Binary indicator of term-vs-preterm birth status; `1` = term birth (>= 37 weeks gestation), `0` = preterm birth (< 37 weeks gestation). Time-fixed per subject. In Allegaert 2015 (paracetamol PK in young women) the indicator is used to select between two typical-value clearances for the sulphate-formation pathway.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (preterm). Effect form is the symmetric `CL = TERM_BIRTH * theta_term + (1 - TERM_BIRTH) * theta_preterm` selection -- so neither category is "the multiplicative reference", and both per-stratum clearances are estimated parameters.
- **Source aliases:**
  - `TERM` (Allegaert 2015 NONMEM column; same orientation, no transformation) -- used in `Allegaert_2015_paracetamol.R`.
- **Example models:** `Allegaert_2015_paracetamol.R`.
- **Notes:** Distinct from `GA` (continuous gestational age in weeks): `TERM_BIRTH` is the binarized version with the conventional 37-week cutoff. Use `GA` when the source paper carries gestational age as a continuous covariate; use `TERM_BIRTH` only when the paper itself dichotomizes. Do not derive `TERM_BIRTH` from `GA` programmatically inside `model()` -- the term-cutoff convention belongs in data assembly, not the model file.

### CONMED_BIRTHCONTROL (**canonical for oral-contraceptive use indicator**)
- **Description:** Binary indicator of oral hormonal contraceptive use; `1` = currently taking an oral contraceptive (estrogen-progestin or progestin-only pill), `0` = not on hormonal contraception. Time-varying as women cycle on/off contraception across study occasions.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no oral contraceptive). Effect form in Allegaert 2015 is multiplicative: `CL_glucuronide *= theta_CONMED_BIRTHCONTROL` when `CONMED_BIRTHCONTROL == 1`, with `theta_CONMED_BIRTHCONTROL = 1.46` (estrogen-driven UGT2B7 induction).
- **Source aliases:**
  - `BC` (Allegaert 2015 NONMEM column; same orientation, no transformation) -- used in `Allegaert_2015_paracetamol.R`.
- **Example models:** `Allegaert_2015_paracetamol.R`.
- **Notes:** Specific scope because the canonical encoding pools all oral contraceptive types (combined / progestin-only) into a single binary; future models that need to distinguish formulations should register a finer-grained canonical (e.g., `CONMED_BIRTHCONTROL_COMBINED`, `CONMED_BIRTHCONTROL_PROGESTIN`). The full-word canonical name was chosen over a shorter `BC_USE` form for clarity in source traces.

### DIS_EOPE (**canonical for early-onset pre-eclampsia indicator**)
- **Description:** Binary indicator of early-onset pre-eclampsia (eoPE); `1` = eoPE diagnosed before 34 weeks gestation, `0` = not eoPE. Time-fixed per subject within the gestational PK study window. Used by population PK models that compare drug disposition in pregnant women with vs without early-onset pre-eclampsia.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no eoPE). Effect form in Schoenmakers 2025 is multiplicative on CL: `CL_eoPE = CL * ThetaPE^DIS_EOPE` with `ThetaPE = 0.617` (38% reduction in betamethasone CL when eoPE is present); encoded in nlmixr2 as the log-additive shift `cl <- exp(lcl + etalcl + e_eope_cl * DIS_EOPE)` with `e_eope_cl = log(0.617)`.
- **Source aliases:**
  - `PE` -- common abbreviation in obstetric pharmacology papers when the cohort restriction is to early-onset PE only; used in `Schoenmakers_2025_betamethasone.R` (paper notation: `eoPE` / `ThetaPE`).
- **Example models:** `Schoenmakers_2025_betamethasone.R` (multiplicative effect on CL, encoded via the log-additive form on `lcl`; reduces apparent betamethasone clearance from 15.6 L/h to 9.6 L/h).
- **Notes:** Distinct from a broader `PREECL` indicator that would pool early-onset, late-onset and postpartum pre-eclampsia. The "early-onset" specifier corresponds to diagnosis before 34 weeks gestation, the conventional clinical cutoff (Phipps 2019 Nat Rev Nephrol). Future papers that enrol mixed early-/late-onset cohorts or that report PE status without the 34-week stratification should register a separate canonical (e.g., `PREECL` for any-onset PE, or `LOPE` for late-onset PE) rather than reusing `DIS_EOPE` with relaxed semantics. Distinct from `PREG` (pregnancy status indicator): `DIS_EOPE` is a complication-of-pregnancy stratifier within a pregnant cohort, whereas `PREG` discriminates pregnant-vs-non-pregnant subjects. Ratified canonically on 2026-05-11 alongside the Schoenmakers 2025 betamethasone extraction.

## Vital signs

### BODYTEMP (**canonical for body temperature**)
- **Description:** Subject body temperature (typically axillary or oral) at the relevant clinical observation. Captured at study admission in acute-infection PK studies (fever as a marker of acute illness severity); may be time-varying when serial temperature measurements are recorded across visits. Document baseline-vs-time-varying status in `covariateData[[BODYTEMP]]$notes` per model.
- **Units:** degC
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with linear-deviation forms `(1 + e * (BODYTEMP - ref))`. Reference values observed: 36.9 degC (Kloprogge 2013 lumefantrine; pooled-cohort median in Ugandan pregnant + non-pregnant women with uncomplicated P. falciparum malaria).
- **Source aliases:**
  - `TEMP` -- common short form in malaria / infectious-disease NONMEM control streams; used in `Kloprogge_2013_lumefantrine.R` (same orientation as the canonical, no value transformation).
- **Example models:** `Kloprogge_2013_lumefantrine.R` (linear-deviation effect on mean absorption transit time MTT: `MTT_indiv = TVMTT * (1 + e_bodytemp_mtt * (BODYTEMP - 36.9))` with `e_bodytemp_mtt = 0.165` per degC; mean transit time increases ~16.5% per degC over 36.0-39.8 degC, plausibly reflecting reduced gut motility / prolonged absorption in feverish malaria patients).
- **Notes:** General scope because body temperature is a universally applicable vital sign; the Kloprogge 2013 reference value 36.9 degC is cohort-specific (Ugandan malaria cohort median) and future models should document their own reference in `covariateData[[BODYTEMP]]$notes`. Units are degrees Celsius; convert from Fahrenheit (degF) at data-assembly time, not inside `model()`. Distinct from `BODYTEMP_FEBRILE` (not yet registered) which would be a binary fever indicator if a future paper dichotomises at the conventional 37.5 / 38.0 degC threshold. Ratified canonically on 2026-05-16 alongside the Kloprogge 2013 lumefantrine extraction.

## Renal / hepatic function

### URINE_FLOW (**canonical for instantaneous urine flow rate**)
- **Description:** Instantaneous urine flow rate (mL/h) measured over the urine collection interval that includes the current observation. Time-varying.
- **Units:** mL/h
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with a centered-linear-effect form `CL_renal = base + theta_URINE_FLOW * (URINE_FLOW - URINE_FLOW_ref)` with `URINE_FLOW_ref = 100 mL/h` in Allegaert 2015. A value of `0` is a sentinel for "no urine collected during the interval" (i.e., the urine pathway contribution is dropped); the linear-effect term is gated by `URINE_FLOW > 0` and not extrapolated below the centering reference.
- **Source aliases:**
  - `UF` (Allegaert 2015 NONMEM column; same orientation, no transformation) -- used in `Allegaert_2015_paracetamol.R`.
- **Example models:** `Allegaert_2015_paracetamol.R`.
- **Notes:** Specific scope because the centered-linear effect form with the `URINE_FLOW == 0` sentinel-zero rule reflects an Allegaert-specific convention rather than a universally-agreed-upon parameterization. A second model that uses a different effect form (e.g., direct `URINE_FLOW / URINE_FLOW_ref` proportional scaling, no zero-sentinel) should register its own canonical (e.g., `URINE_FLOW_PROP`) rather than reusing `URINE_FLOW` with conflicting semantics. The full-word canonical name was chosen over the bare `UF` source-data abbreviation for clarity in source traces.

### CRCL (**canonical for creatinine-based renal function, BSA-normalized**)
- **Description:** Creatinine-based renal function expressed in mL/min/1.73 m^2. Accepts either an MDRD-/CKD-EPI-estimated glomerular filtration rate or a measured creatinine clearance that has been BSA-normalized as `1.73 x CrCl / BSA`. The per-model `covariateData[[CRCL]]$description` and `notes` must state which method the source paper used.
- **Units:** mL/min/1.73 m^2
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(CRCL / ref)^exponent`. Reference values observed: 80 mL/min/1.73 m^2 (Cirincione 2017, MDRD eGFR), 90 mL/min/1.73 m^2 (Li 2019, calculated GFR), 100 mL/min/1.73 m^2 (Xu 2019, measured-CrCl BSA-normalized).
- **Reference category:** n/a -- used with power scaling `(CRCL / ref)^exponent`. Reference values observed: 80 mL/min/1.73 m^2 (Cirincione 2017, MDRD eGFR), 100 mL/min/1.73 m^2 (Xu 2019, measured-CrCl BSA-normalized), 90 mL/min/1.73 m^2 (Bajaj 2017, CKD-EPI eGFR).
- **Source aliases:**
  - `eGFR` -- MDRD-estimated glomerular filtration rate; used in `Cirincione_2017_exenatide.R` and `Kotani_2022_astegolimab.R`. `Bajaj_2017_nivolumab.R` uses the CKD-EPI variant.
  - `EGFR` -- all-caps variant.
  - `CRCL_BSA` -- BSA-normalized creatinine clearance (measured CrCl / BSA x 1.73); used in `Xu_2019_sarilumab.R`.
  - `1.73*CrCl/BSA` -- the formula form appearing in Xu 2019 Eq. for Vm.
  - `cGFR` -- calculated/estimated GFR, BSA-normalized; used in `Li_2019_abatacept.R`.
  - `CLCR` -- raw Cockcroft-Gault creatinine clearance in mL/min (NOT BSA-normalized); used in `Delattre_2010_amikacin.R` (median 55.5 mL/min in critically ill septic adults).
- **Example models:** `Cirincione_2017_exenatide.R` (MDRD eGFR), `Xu_2019_sarilumab.R` (measured CrCl BSA-normalized), `Kotani_2022_astegolimab.R` (MDRD eGFR), `Li_2019_abatacept.R` (cGFR), `Bajaj_2017_nivolumab.R` (CKD-EPI eGFR, reference 90 mL/min/1.73 m^2), `NA_NA_lidocaine.R` (DDMODEL00000281; binary stratification at threshold 52.7 mL/min adding -0.319 to the GX rate constant K30 in the CRCL <= 52.7 cohort; the source `.ctl` does not state the BSA-normalisation method), `Delattre_2010_amikacin.R` (raw Cockcroft-Gault mL/min, NOT BSA-normalized; reference 55.5 mL/min population median; additive linear effect 1.42 L/h per (CRCL/55.5) on CL).
- **Notes:** The two estimation methods (MDRD/CKD-EPI vs measured CrCl) produce values in the same units and are operationally interchangeable as a covariate on clearance. Document the method explicitly in each model's `covariateData[[CRCL]]$description` so future reviewers can trace the source assay.

### CREAT (**canonical for serum creatinine**)
- **Description:** Serum creatinine concentration (baseline or time-varying).
- **Units:** umol/L or mg/dL -- document the unit used in each model via `covariateData[[CREAT]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(CREAT / ref)^exponent`.
- **Source aliases:**
  - `CRE` (umol/L, reference 70.73) -- used in `Thakre_2022_risankizumab.R`.
  - `SCR` -- common clinical-PK abbreviation; also Llanos-Paez 2020 source column for the patient's individual serum creatinine.
- **Example models:** `Thakre_2022_risankizumab.R`, `Hennig_2013_tobra.R` (umol/L; paired with `CREAT_REF` for the SCR_mean / SCR ratio used in the Hennig 2013 renal-function factor), `Llanos_2017_gentamicin.R` (umol/L; standardized per-patient against `CREAT_REF` rather than a fixed cohort reference), `Llanos-Paez_2020_gentamicin.R` (umol/L; used as the patient's `SCR_i` in the renal-function ratio `(CREAT_REF / CREAT)^0.58` on CL).
- **Notes:** `CREAT` chosen over the shorter `CRE`/`SCR` as the NONMEM/clinical-PK convention that is unambiguous. Per-model reference values must be documented in `covariateData[[CREAT]]$notes`.

### CREAT_REF (**canonical for sex/age/size-expected normal-mean serum creatinine**)
- **Description:** Externally-computed reference serum creatinine for the individual (the expected normal SCR for a healthy person of the same sex, age and body size). Used as the numerator of a ratio against the patient's measured `CREAT` to define a renal-function factor on clearance.
- **Units:** umol/L or mg/dL -- must match the unit of the paired `CREAT` column. Document via `covariateData[[CREAT_REF]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used as `(CREAT_REF / CREAT)^exponent` so that a patient with measured SCR equal to the population-expected normal SCR has factor 1.
- **Source aliases:**
  - `SCR_mean` -- used in `Hennig 2013` (Eq. 5: `f_SCR = (SCR_mean / SCR)^theta_SCR`); also Llanos-Paez 2020 paper notation for the Ceriotti 2008 age/sex-matched physiological mean SCR.
  - `Scrmean` -- Llanos-Paez 2017 paper notation; computed from Ceriotti et al. 2008 age- and sex-stratified medians (Clin Chem 54:559-566, doi:10.1373/clinchem.2007.099648).
- **Example models:** `Hennig_2013_tobra.R`, `Llanos_2017_gentamicin.R` (umol/L; computed externally per Ceriotti et al. 2008), `Llanos-Paez_2020_gentamicin.R` (umol/L; ratio `(CREAT_REF / CREAT)^0.58` multiplies the maturation-scaled CL).
- **Notes:** Specific scope because the formula used to derive the reference value is paper-defined (Hennig 2013 cites a combination of Ceriotti 2008, Junge 2004 and Johansson 2011 reference-interval relationships; Llanos-Paez 2017 and 2020 both use Ceriotti 2008); a future paper that uses a different reference-SCR derivation (e.g., a CKD-EPI-style adult-only reference, or a Schwartz-derived paediatric-only reference) should pin its formula in `covariateData[[CREAT_REF]]$notes` so that a user assembling a virtual cohort can reproduce it. When no covariate data are available to compute `CREAT_REF`, set `CREAT_REF = CREAT` so the renal-function factor evaluates to 1 (matching the Hennig 2013 'covariate set to 1 for missing data' rule). Ratified canonically on 2026-05-08 alongside the Hennig 2013 tobramycin extraction.

### HEMODIAL (**canonical for intermittent-hemodialysis treatment-status indicator**)
- **Description:** 1 = the subject was undergoing intermittent hemodialysis during the modeled period; 0 = no intermittent hemodialysis. Treatment-status flag rather than a measured renal-function value; used as a multiplicative covariate on PK parameters that change with chronic dialysis (typically CL and Vc -- intermittent hemodialysis decreases vancomycin-class CL and reduces interstitial volume overload, lowering Vc). Per-subject indicator in the source data; in Goti 2018 it is treated as time-fixed at the subject level (the cohort either was or was not receiving intermittent hemodialysis during the admission), and Goti 2018 explicitly notes that actual hemodialysis-session timing was not used because of documentation limitations.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no intermittent hemodialysis).
- **Source aliases:**
  - `DIAL` -- used in `Goti_2018_vancomycin.R` (binary indicator on CL and Vc in a 2-compartment vancomycin popPK model). Goti 2018 Methods notes the indicator was created for the routine-TDM cohort (n = 336 hemodialysis subjects of 1812 total) and that all hemodialysis procedures were intermittent and used high-flux membranes.
- **Example models:** `Goti_2018_vancomycin.R` (multiplicative factors on CL and Vc: `0.7^HEMODIAL` on CL and `0.5^HEMODIAL` on Vc, so dialysis subjects have 30% lower CL and 50% lower central volume than non-dialysis subjects).
- **Notes:** Specific to intermittent hemodialysis (IHD). Distinct from peritoneal dialysis (PD) and from continuous renal replacement therapy (CRRT), each of which has different drug-extraction kinetics and would warrant its own canonical (`PERIT_DIAL`, `CRRT_STATUS`) if a future paper retains them as covariates. Goti 2018 treats `HEMODIAL` as time-fixed per subject because session-level dialysis timing was not reliably documented in the source EHR data; a future paper that resolves drug clearance during versus between dialysis sessions would use a time-varying form (or a separate per-session covariate) and the per-model `covariateData[[HEMODIAL]]$notes` would document the time resolution. When pairing `HEMODIAL` with `CRCL`, note that the Cockcroft-Gault CRCL of an anuric hemodialysis patient is by convention very low or set per institution to a small floor value (Goti 2018 truncated CRCL > 150 mL/min to 150 mL/min and corrected SCr < 1 mg/dL in elderly subjects); residual renal function in hemodialysis subjects is highly variable and the dialysis indicator captures the bulk PK shift on top of the CRCL covariate. Ratified canonically on 2026-05-16 alongside the Goti 2018 vancomycin extraction.

### ALB (**canonical for serum albumin**)
- **Description:** Serum albumin concentration.
- **Units:** g/dL or g/L -- document the unit used in each model via `covariateData[[ALB]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(ALB / ref)^exponent`.
- **Source aliases:**
  - `BALB` (baseline albumin) -- used in `Zhou_2021_belimumab.R`. Maps directly to `ALB`; baseline-vs-time-varying status documented in per-model notes.
- **Example models:** `Fasanmade_2009_infliximab.R` (g/dL, reference 4.1), `Thakre_2022_risankizumab.R` (g/L, reference 45), `Chua_2025_mirikizumab.R`, `Moein_2022_etrolizumab.R`, `Tiraboschi_2025_amlitelimab.R`, `Yamada_2025_zolbetuximab.R`, `Li_2019_abatacept.R` (g/dL, reference 4.0; the Li 2019 Methods states 'mg/dL' which is a publication typo -- see the model's `covariateData[[ALB]]$notes`), `Quartino_2019_trastuzumab.R` (g/dL, reference 4; source column `ALBU`; negative exponent -0.998 on linear CL), `Wang_2020_ontamalimab.R` (g/L, reference 39), `Zhou_2021_belimumab.R` (g/L, reference 40; baseline-only, source column `BALB`), `Okada_2025_rocatinlimab.R` (g/L, reference 44; source column `ALBU`; power exponent -1.30 on linear CL), `Xu_2020_daratumumab.R` (g/L, reference 37.0; power exponent -1.149 on linear CL).
- **Notes:** Ratified canonically on 2026-04-19 after cross-model review. Unit varies by paper (g/dL in US-convention papers, g/L in SI-convention papers); the per-model `covariateData[[ALB]]$units` field is load-bearing. Effect-coefficient magnitude is meaningless without the unit.

### TPRO (**canonical for total serum protein**)
- **Description:** Total serum protein concentration (sum of albumin + globulins; baseline or time-varying).
- **Units:** g/L or g/dL -- document the unit used in each model via `covariateData[[TPRO]]$units` (1 g/dL = 10 g/L).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(TPRO / ref)^exponent`. Reference value observed: 74 g/L (Frey 2010 pooled-cohort median).
- **Source aliases:**
  - `PROT` -- Frey 2010 abbreviation in the final-model equation.
  - `TP` -- common clinical-chemistry abbreviation.
- **Example models:** `Frey_2010_tocilizumab.R` (g/L, reference 74; exponent -1.1 on V1).
- **Notes:** Distinct from `ALB` (serum albumin, the largest single component of total protein). Frey 2010 retains both `TPRO` and `ALB` on V1 as separate covariates with opposite signs (TPRO negative, ALB positive) and notes there is no clear mechanistic explanation; the joint effect may reflect serum-volume modifications. `TPRO` ratified canonically on 2026-04-28 alongside the Frey 2010 extraction.

### IGG (**canonical for serum immunoglobulin G**)
- **Description:** Serum total immunoglobulin G concentration (baseline or time-varying). Used in mAb PK analyses as a competition-for-FcRn-recycling covariate on therapeutic-mAb clearance -- high endogenous IgG is hypothesized to displace the therapeutic mAb from FcRn salvage and increase its catabolic clearance.
- **Units:** g/L (typical in SI-convention papers); also reported as mg/dL in US-convention papers -- document the unit used in each model via `covariateData[[IGG]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(IGG / ref)^exponent`. Reference values observed: 14.8 g/L (Zhou 2021), 9.65 g/L (Yang 2021).
- **Source aliases:**
  - `BIGG` (baseline IgG) -- used in `Zhou_2021_belimumab.R`.
  - `IGGBL` (baseline IgG) -- used in `Yang_2021_cemiplimab.R`.
- **Example models:** `Zhou_2021_belimumab.R` (g/L, reference 14.8; baseline-only; exponent 0.293 on CL), `Yang_2021_cemiplimab.R` (g/L, reference 9.65; small positive exponent 0.184 on shared CL/Q).
- **Notes:** Mechanistically meaningful for monoclonal-antibody PK because endogenous IgG competes with the therapeutic mAb for FcRn-mediated recycling. The per-model `covariateData[[IGG]]$units` field is load-bearing (1 g/L ~= 100 mg/dL). Baseline-vs-time-varying status documented in `covariateData[[IGG]]$notes`. Distinct from `lIgG0` / IgG-as-a-state in mechanistic FcRn-competition TMDD models (e.g., `Valenzuela_2025_nipocalimab.R`), where IgG is a dynamic state, not a baseline covariate; use `IGG` only when the source paper treats IgG as a static (baseline) covariate column.

### IGM (**canonical for serum immunoglobulin M**)
- **Description:** Serum total immunoglobulin M (IgM) concentration (baseline). Used in IgRT population-PK analyses as a proxy for B-cell antibody-producing capacity / humoral function -- IgM is the first antibody produced after B-cell activation, so circulating IgM reflects ongoing B-cell activity prior to class-switching to IgG.
- **Units:** g/L (typical SI-convention reporting); also reported as mg/dL in US-convention papers -- document the unit used in each model via `covariateData[[IGM]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with power scaling `(IGM / ref)^exponent`. Reference values observed: 0.21 g/L (Cheng 2026, pooled PID + SAD pediatric cohort median).
- **Source aliases:** none known.
- **Example models:** `Cheng_2026_immunoglobulin.R` (g/L, reference 0.21; baseline-only; power exponent 0.11 on baseline IgG (CBAS) -- IgM enters as a humoral-capacity proxy that informs the endogenous-IgG baseline rather than directly modifying clearance).
- **Notes:** IgM is the immune-globulin class produced by activated B cells before class-switching, so it remains detectable in patients with hypogammaglobulinaemia who still have residual B-cell function. Scope: specific because the relevance of IgM as a covariate depends on the paper's mechanistic interpretation (in Cheng 2026 it acts on the endogenous-IgG baseline; future use cases may differ). Promote to general if a second paper retains IgM with consistent semantics. Ratified canonically on 2026-04-28.

### TBILI (**canonical for total bilirubin**)
- **Description:** Total serum bilirubin concentration.
- **Units:** mg/dL or umol/L -- document the unit used in each model via `covariateData[[TBILI]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(TBILI / ref)^exponent`.
- **Source aliases:**
  - `BIL` (legacy NONMEM short label for total bilirubin) -- used in `NA_NA_lidocaine.R` (DDMODEL00000281; binarised at threshold 0.53 mg/dL with `BIL_HIGH = as.integer(BIL > 0.53)`).
  - `BILT` (Urien 2005 capecitabine paper's NONMEM short label for "total bilirubin") -- used in `Urien_2005_capecitabine.R` (umol/L, reference 8.8; power scaling on the capecitabine non-transformation CL10 and on the 5'-DFUR -> 5-FU rate constant K34).
- **Example models:** `Yamada_2025_zolbetuximab.R` (mg/dL, reference 0.38; small positive exponent 0.0347 on V1), `NA_NA_lidocaine.R` (mg/dL, source column `BIL`; binary effect at threshold 0.53 mg/dL on the GX elimination rate constant K30), `Urien_2005_capecitabine.R` (umol/L, reference 8.8; source column `BILT`; positive exponent +0.32 on capecitabine non-transformation CL10 and negative exponent -0.36 on the 5'-DFUR -> 5-FU rate constant K34).
- **Notes:** Hepatic-function marker. Unit varies by paper (US convention mg/dL, SI convention umol/L; 1 mg/dL ~= 17.1 umol/L). The per-model `covariateData[[TBILI]]$units` field is load-bearing.

### DBIL (**canonical for direct (conjugated) bilirubin**)
- **Description:** Direct (conjugated) serum bilirubin concentration. Distinct from `TBILI`: direct bilirubin is the water-soluble glucuronide-conjugated fraction processed by hepatocytes and excreted in bile, so a rise in DBIL specifically flags impaired biliary excretion / cholestasis or intrahepatic shunting, whereas total bilirubin also captures unconjugated (indirect) hyperbilirubinaemia from haemolysis or Gilbert-type conjugation defects.
- **Units:** mg/dL or umol/L -- document the unit used in each model via `covariateData[[DBIL]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with power scaling `(DBIL / ref)^exponent`. Reference values observed: 2.6 umol/L (Chen 2015 voriconazole Chinese ICU cohort population median).
- **Source aliases:** none known.
- **Example models:** `Chen_2015_voriconazole.R` (umol/L, reference 2.6; negative exponent -0.40 on CL: `CL = TVCL * (DBIL / 2.6)^-0.40`).
- **Notes:** Hepatic-function / cholestasis-specific marker. Unit varies by paper (US convention mg/dL, SI convention umol/L; 1 mg/dL ~= 17.1 umol/L). Distinct entry from `TBILI` because direct vs total are not interchangeable: total = direct + indirect, and the two fractions track different pathophysiologic processes. Scope kept `specific` pending a second model that ratifies DBIL with consistent semantics; promote to `general` once corroborated.

### AST (**canonical for aspartate aminotransferase**)
- **Description:** Serum aspartate aminotransferase activity (baseline or time-varying).
- **Units:** U/L (IU/L; the two labels are used interchangeably in the clinical-PK literature). Document per-model via `covariateData[[AST]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(AST / ref)^exponent`.
- **Source aliases:**
  - `SGOT` (serum glutamic-oxaloacetic transaminase; the legacy clinical-chemistry name for AST) -- used in `Quartino_2019_trastuzumab.R`.
- **Example models:** `Lu_2014_trastuzumabemtansine.R` (U/L, reference 27; small positive exponent 0.071 on CL), `Quartino_2019_trastuzumab.R` (IU/L, reference 24; source column `SGOT`; positive exponent 0.205 on linear CL).
- **Notes:** Hepatic-function marker. Commonly reported alongside `ALT` and `TBILI`; register a separate `ALT` canonical if a future paper requires it. `SGOT` is the older lab-reporting name; values and units are identical to `AST`.

### ALT (**canonical for alanine aminotransferase**)
- **Description:** Serum alanine aminotransferase activity (baseline or time-varying).
- **Units:** U/L (IU/L; the two labels are used interchangeably in the clinical-PK literature). Document per-model via `covariateData[[ALT]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(ALT / ref)^exponent`.
- **Source aliases:**
  - `SGPT` (serum glutamic-pyruvic transaminase; the legacy clinical-chemistry name for ALT, paralleling `SGOT` -> `AST`) -- used in `NA_NA_lidocaine.R` (DDMODEL00000281; binarised at threshold 11 with `SGPT_HIGH = as.integer(SGPT > 11)`).
- **Example models:** `Nikanjam_2019_siltuximab.R` (U/L, reference 19; small negative exponent -0.096 on CL), `Melhem_2022_dostarlimab.R` (U/L, reference 18; small negative exponent -0.0585 on CL, time-varying), `NA_NA_lidocaine.R` (source column `SGPT`; binary effects at threshold 11 on the GX rate constant K30 and on the 2,6-xylidide rate constant K40).
- **Notes:** Hepatic-function marker. Commonly reported alongside `AST` and `TBILI`. Ratified canonically on 2026-04-24. `SGPT` is the older lab-reporting name; values and units are identical to `ALT`.

### ALP (**canonical for alkaline phosphatase**)
- **Description:** Serum alkaline phosphatase activity (baseline or time-varying). Liver-function / cholestasis marker; often used in popPK covariate models either as a continuous concentration with power scaling or as a binary above/below upper-limit-of-normal (ULN) indicator. When binarized inline, document the ULN threshold used (typically ~120 U/L for adults; varies by lab, age, and sex).
- **Units:** U/L (IU/L; the two labels are used interchangeably in the clinical-PK literature). Document per-model via `covariateData[[ALP]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(ALP / ref)^exponent`, with a linear-deviation form, or binarized inline as `alp_high <- (ALP > uln)` for a binary >ULN indicator.
- **Source aliases:** none known.
- **Example models:** `Gupta_2016_lenvatinib.R` (binarized inline as `alp_high <- (ALP > 120)`; the source paper enters `ALP` as a 0/1 NONMEM indicator with `ALP = 1` when the ratio ALP/ULN > 1; multiplicative effect on CL/F: `0.883^alp_high`).
- **Notes:** Liver-function / cholestasis marker; routine clinical-chemistry covariate. Commonly tested alongside `ALT` / `AST` / `GGT` / `TBILI`. Ratified canonically alongside the Gupta 2016 lenvatinib extraction.

### GGT (**canonical for gamma-glutamyltransferase**)
- **Description:** Serum gamma-glutamyltransferase activity (baseline or time-varying); hepatic / cholestatic biliary-enzyme marker.
- **Units:** U/L (IU/L; the two labels are used interchangeably in the clinical-PK literature). Document per-model via `covariateData[[GGT]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with linear-deviation form `1 + theta * (GGT - ref)` or power scaling `(GGT / ref)^exponent`. Reference values observed: 33 U/L (Retlich 2015 popPK linagliptin median), 32.3 U/L (Retlich 2015 popPK/PD linagliptin median).
- **Source aliases:** none known.
- **Example models:** `Retlich_2015_linagliptin.R` (U/L, reference 33; linear-deviation effect on linagliptin CL with coefficient -0.0339 % per U/L deviation. The PK/PD layer uses GGT (reference 32.3 U/L) as a piecewise covariate on baseline DPP-4 activity BSL with a linear-deviation effect below GGT = 175 U/L and a constant +21.3% effect above the threshold).
- **Notes:** Liver-function / cholestasis marker; routine clinical-chemistry covariate. Commonly tested alongside `ALT` / `AST` / `ALP` / `TBILI`. The piecewise above/below-threshold form in Retlich 2015 reflects empirical saturation of the GGT-vs-DPP-4-activity relationship at extreme values. Ratified canonically alongside the Retlich 2015 linagliptin extraction.

### LDH (**canonical for serum lactate dehydrogenase**)
- **Description:** Serum lactate dehydrogenase activity (baseline or time-varying). General-purpose marker of tissue / cellular turnover; in oncology PK analyses it is interpreted as a disease-burden / cell-turnover proxy.
- **Units:** U/L (IU/L; interchangeable). Document per-model via `covariateData[[LDH]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(LDH / ref)^exponent` or with an additive linear-on-log form `exp(coef * (log(LDH) - log(ref)))` (algebraically equivalent to `(log(LDH) / log(ref))^coef`). Reference values observed: 217 U/L (Sanghavi 2020).
- **Source aliases:**
  - `BLDH` (baseline LDH) -- used in `Sanghavi_2020_ipilimumab.R`.
- **Example models:** `Sanghavi_2020_ipilimumab.R` (linear-on-log form on CL with reference 217 U/L; coefficient 0.703), `NA_NA_lidocaine.R` (DDMODEL00000281; binary stratification at threshold 195 U/L switching the typical-value baseline of the 2,6-xylidide rate constant K40).
- **Notes:** Universal lab marker. Sanghavi 2020 log-transforms LDH because the distribution is heavily right-skewed (range 74-6,245 U/L over a median of 217); other papers may use a simple `(LDH/ref)^exponent` form. Document the functional form in `covariateData[[LDH]]$notes`.

### HEPIMP_MILD (**canonical for mild hepatic impairment indicator**)
- **Description:** 1 = mild hepatic impairment per the National Cancer Institute Organ Dysfunction Working Group (NCI ODWG) criteria, 0 = normal hepatic function or non-mild category. NCI ODWG mild = total bilirubin <= ULN with AST > ULN, OR total bilirubin > 1.0xULN to <= 1.5xULN with any AST.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (normal hepatic function; the moderate / severe categories are typically pooled into the reference for population PK analyses where mild impairment is the only category with non-trivial sample size).
- **Source aliases:**
  - `HEPIMP` (with values `1 = mild / 0 = others`) -- used in `Lin_2024_casirivimab.R`.
- **Example models:** `Lin_2024_casirivimab.R` (multiplicative fractional change on CL), `Lu_2022_patritumab.R` (paired with `HEPIMP_MOD_MISSING`; multiplicative fractional effect 0.706 on CLDXd for mild impairment vs the normal-hepatic-function reference).
- **Notes:** Use this column when a model dichotomizes hepatic-impairment status as "mild vs. others" (i.e., normal + the rare moderate/severe cases pooled into the reference). For models that test moderate or severe as separate categories, register additional canonicals `HEPIMP_MOD` / `HEPIMP_SEV` rather than overloading this entry.

### HEPIMP_MOD_MISSING (**canonical for composite moderate-or-data-missing hepatic impairment indicator**)
- **Description:** 1 = moderate hepatic impairment per the NCI ODWG criteria OR baseline hepatic-function data missing/unknown; 0 = normal hepatic function or any other (non-moderate, non-missing) category. Composite indicator used by source papers that pool the moderate-impairment subgroup with patients whose hepatic-function data are missing because both subgroups are individually too small to estimate as separate effects.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (normal hepatic function; mild impairment is typically captured by a separate `HEPIMP_MILD` indicator paired with this column, so all-zero corresponds to NCI ODWG group 1 = normal).
- **Source aliases:**
  - `HEPATIC` / `HEPATIC_MOD_MISSING` -- informal NONMEM names for the composite group.
- **Example models:** `Lu_2022_patritumab.R` (paired with `HEPIMP_MILD`; multiplicative fractional effect 0.532 on CLDXd; the composite group pools n = 6 moderate-impairment patients with n = 6 missing/unknown patients per Lu 2022 Table S5).
- **Notes:** Specific scope because the composition of the "moderate or missing" group is paper-defined and the missing/unknown subgroup may have a different distribution of true hepatic-function status across studies. Use only when the source paper explicitly pools the moderate-impairment cases with missing-data cases under a single coefficient; for models that estimate moderate impairment separately (without pooling missing data), register a `HEPIMP_MOD` canonical instead. Ratified canonically on 2026-04-28.

### B2M (**canonical for serum beta-2-microglobulin**)
- **Description:** Serum beta-2-microglobulin concentration. Low-molecular-weight (~12 kDa) protein freely filtered at the glomerulus and reabsorbed in the proximal tubule; serum levels rise with renal impairment, with increased plasma-cell turnover in multiple myeloma, and with broader lymphoid-cell turnover. Used in oncology PK analyses both as a renal-function proxy and as a tumor-burden / disease-severity covariate.
- **Units:** mg/L
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(B2M / ref)^exponent`. Reference values observed: 3.90 mg/L (Fau 2020 multiple-myeloma cohort median).
- **Source aliases:** none; `B2M` is the universal abbreviation.
- **Example models:** `Fau_2020_isatuximab.R` (mg/L, reference 3.90; exponent 0.343 on the steady-state linear clearance CLinf).
- **Notes:** In multiple myeloma B2M is part of the International Staging System (ISS); in routine PK analyses it is interpreted simultaneously as a renal-function and disease-burden marker. Document the interpretation per-model via `covariateData[[B2M]]$notes`.

### HEPIMP (**canonical for hepatic-impairment indicator (NCI ODWG classification)**)
- **Description:** Baseline hepatic-impairment indicator per the National Cancer Institute Organ Dysfunction Working Group (NCI ODWG) classification: 1 = mild or worse hepatic impairment (group >= 2 = mild, moderate, or severe), 0 = normal hepatic function (group 1).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (normal hepatic function, NCI ODWG group 1).
- **Source aliases:**
  - `BHPTGRPN` (categorical: 1 = normal, 2 = mild, 3 = moderate, 4 = severe; 9999 = missing) -- used in `Lu_2019_polatuzumab.R`. Decompose: `HEPIMP = as.integer(BHPTGRPN > 1.5 & BHPTGRPN != 9999)`.
  - `HEP_IMP` -- retired canonical name; replaced by `HEPIMP` for consistency with the `HEPIMP_MILD` family.
- **Example models:** `Lu_2019_polatuzumab.R` (multiplicative effect on FRAC_NS = 1.19, applied as `1.19^HEPIMP`).
- **Notes:** NCI ODWG classification (Ramalingam SS et al., J Clin Oncol 2010;28:4507) groups subjects by total bilirubin and AST: group 1 = normal, group 2 = mild (TBILI <= ULN and AST > ULN, or TBILI > 1-1.5 x ULN), group 3 = moderate (TBILI > 1.5-3 x ULN), group 4 = severe (TBILI > 3 x ULN). Source papers typically pool groups 2-4 versus group 1 for a binary indicator because the impaired-liver subgroups are individually small. If a future model needs finer resolution (separate effects for mild vs moderate-or-worse), add a parallel `HEPIMP_MOD` canonical rather than overloading this one.

### HEPIMP_SEV (**canonical for severe hepatic impairment indicator**)
- **Description:** 1 = severe hepatic impairment, 0 = normal hepatic function or less-than-severe category. The classification scheme that defines "severe" is paper-specific and must be documented in per-model `covariateData[[HEPIMP_SEV]]$notes`. Two schemes are commonly encountered:
  - **NCI ODWG group 4**: total bilirubin > 3 x ULN with any AST (Ramalingam SS et al., J Clin Oncol 2010;28:4507).
  - **Child-Pugh Class C**: composite score 10-15 across bilirubin, albumin, INR, ascites, and encephalopathy.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (any non-severe category: normal, mild, or moderate; the model typically uses HEPIMP_SEV alongside other severity-specific indicators that partition the non-severe pool further).
- **Source aliases:**
  - `Child-Pugh Class C` -- used in `vanderWalt_2013_dapagliflozin.R` (covariate effects on CLP_M15 and V2M; the paper dichotomizes severe hepatic impairment per the Child-Pugh classification).
- **Example models:** `vanderWalt_2013_dapagliflozin.R` (Child-Pugh Class C; multiplicative fractional effects -0.422 on the dapagliflozin -> D3OG metabolic clearance and +1.33 on the D3OG central volume of distribution; paper text "With severe HI (Child-Pugh Class C), CLP M15 decreased by 41% and V2M increased by 134%").
- **Notes:** Use this column when a model dichotomizes severe hepatic impairment as a separate indicator from milder categories. The classification scheme (NCI ODWG vs Child-Pugh vs other) is paper-specific and must be documented per-model. For composite "moderate-or-severe" pooled indicators, use the parallel `HEPIMP_MODSEV` canonical rather than overloading this entry. Companion to `HEPIMP_MILD` (mild only) and `HEPIMP_MODSEV` (moderate + severe pooled); the SKILL.md anticipates each severity level as its own canonical when the source paper tests them as separate covariates.

### HEPIMP_MODSEV (**canonical for composite moderate-or-severe hepatic impairment indicator**)
- **Description:** 1 = moderate or severe hepatic impairment, 0 = normal hepatic function or mild impairment. Composite indicator used by source papers that pool the moderate and severe subgroups because the severe subgroup alone is too small to support a separate covariate-effect estimate. Distinct from `HEPIMP_MOD_MISSING` (which pools moderate cases with missing-data cases, not with severe cases). The classification scheme that defines the cut points is paper-specific and must be documented in per-model `covariateData[[HEPIMP_MODSEV]]$notes`. Two schemes are commonly encountered:
  - **NCI ODWG groups 3-4 pooled**: total bilirubin > 1.5 x ULN with any AST.
  - **Child-Pugh Class B or C pooled**: composite score >= 7.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (normal hepatic function or mild impairment; the indicator is mutually exclusive with HEPIMP_MILD, so all-zero on both indicators corresponds to the normal-function reference and HEPIMP_MILD=1 with HEPIMP_MODSEV=0 corresponds to mild-only).
- **Source aliases:**
  - `Child-Pugh Class B,C` -- used in `vanderWalt_2013_dapagliflozin.R` (covariate effects on V3P and CLM; the paper dichotomizes moderate-or-severe hepatic impairment per the Child-Pugh classification).
- **Example models:** `vanderWalt_2013_dapagliflozin.R` (Child-Pugh Class B or C; multiplicative fractional effects -0.600 on the dapagliflozin peripheral volume of distribution V3P and -0.293 on the D3OG renal clearance CLM; paper text "Moderate or severe HI (Child-Pugh Class B or C) decreased CLM and the peripheral volume of distribution of dapagliflozin (V3P) by 29 and 60%, respectively").
- **Notes:** Use this column when a model pools moderate-and-severe hepatic impairment under a single coefficient (typically because the severe subgroup alone is too small to estimate as its own effect). The classification scheme (NCI ODWG vs Child-Pugh vs other) is paper-specific and must be documented per-model. Companion to `HEPIMP_MILD` (mild only) and `HEPIMP_SEV` (severe only). Distinct from `HEPIMP_MOD_MISSING` (which pools moderate cases with subjects whose hepatic-function data are missing/unknown, not with severe cases). The composite mod-or-sev pooling is a different load-bearing convention than the mod-or-missing pooling, so the two canonicals must remain separate.

### CPK (**canonical for serum creatine phosphokinase / creatine kinase**)
- **Description:** Serum creatine phosphokinase (also called creatine kinase, CK) activity (baseline or time-varying). Skeletal-muscle / cardiac-muscle injury and turnover marker; in macrophage-targeted PK/PD analyses (axatilimab, anti-CSF-1R) it is interpreted as a Kupffer-cell / tissue-macrophage clearance surrogate because Kupffer cells participate in the elimination of circulating muscle-derived enzymes.
- **Units:** U/L (IU/L; interchangeable). Document per-model via `covariateData[[CPK]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(CPK / ref)^exponent`. Reference values observed: 63 U/L (Yang 2024 axatilimab; pooled-cohort median).
- **Source aliases:**
  - `BLCPK` (baseline CPK) -- informal usage in Yang 2024.
- **Example models:** `Yang_2024_axatilimab.R` (baseline-only covariate on baseline NCMC concentration `BL_NCMC` with power exponent 0.376; reference 63 U/L).
- **Notes:** Muscle-origin enzyme distinct from `AST` / `ALT` (hepatic) and `LDH` (general tissue turnover). Yang 2024 uses CPK alongside `AST` and `LDH` as tracked safety biomarkers. Per-model `covariateData[[CPK]]$notes` should document baseline-vs-time-varying status and the clinical interpretation in the source population (skeletal-muscle injury, macrophage-clearance surrogate, or both). Distinct from any model state variable representing CPK time-course dynamics -- covariate column is the pre-dose laboratory observation.

### DIAL (**canonical for hemodialysis-active indicator (time-varying)**)
- **Description:** Within-subject time-varying indicator for whether an extracorporeal renal-replacement-therapy session (intermittent hemodialysis, hemofiltration, or hemodiafiltration) is currently running. 1 during the session; 0 in the interdialytic interval and in non-dialysed subjects.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no dialysis running). Models that compose an additional `CL_dialysis` term should gate it by `DIAL = 1` so that the interdialytic clearance reduces to the intrinsic body clearance.
- **Source aliases:** none known; the source paper for the ratification model uses `DIAL` as the data column name directly (Liesenfeld 2013 Methods).
- **Example models:** `Liesenfeld_2013_dabigatran.R` (Michaels-equation gate; `cl_total <- cl + DIAL * Michaels(BFR, DFR, KoA)`).
- **Notes:** Distinct from a renal-impairment indicator (subject-level baseline class) and from a renal-replacement-therapy modality indicator (categorical: IHD vs CRRT vs SLED) -- `DIAL` is the per-time-point gate that turns the dialysis-clearance term on and off. Pair with `BFR` and `DFR` when the dialysis clearance depends on flow rates; pair with a filter-specific mass-transfer coefficient (estimated `lkoa` in the model, not a covariate) when the Michaels parameterisation is used. Ratified canonically on 2026-05-16 alongside the Liesenfeld 2013 dabigatran extraction.

### BFR (**canonical for blood flow rate through the extracorporeal circuit during dialysis**)
- **Description:** Instantaneous blood flow rate through the extracorporeal circuit during an active dialysis session. Time-varying within subject; meaningful only when `DIAL = 1` -- in the interdialytic period the value is sentinel and the Michaels-equation term is gated off by `DIAL`.
- **Units:** mL/min
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- enters the Michaels equation together with `DFR` and a hemodialyzer mass-transfer-area coefficient. Values investigated in the ratification source were 200, 300, and 400 mL/min (Liesenfeld 2013 Methods, Study Design; Table 1).
- **Source aliases:** none known.
- **Example models:** `Liesenfeld_2013_dabigatran.R`.
- **Notes:** Pairs with `DIAL` (binary on/off gate) and `DFR` (dialysate flow rate). Ratified canonically on 2026-05-16 alongside the Liesenfeld 2013 dabigatran extraction.

### DFR (**canonical for dialysate flow rate through the extracorporeal circuit during dialysis**)
- **Description:** Instantaneous dialysate flow rate through the extracorporeal circuit during an active dialysis session. Time-varying within subject; meaningful only when `DIAL = 1`.
- **Units:** mL/min
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- enters the Michaels equation together with `BFR`. The ratification source fixed DFR at 700 mL/min throughout (Liesenfeld 2013 Methods, Study Design) and additionally simulated 500 mL/min (Methods, Simulations).
- **Source aliases:** none known.
- **Example models:** `Liesenfeld_2013_dabigatran.R`.
- **Notes:** Pairs with `DIAL` and `BFR`. Ratified canonically on 2026-05-16 alongside the Liesenfeld 2013 dabigatran extraction.

## Hematology

### HGB (**canonical for hemoglobin**)
- **Description:** Blood hemoglobin concentration.
- **Units:** g/L or g/dL -- document the unit used in each model via `covariateData[[HGB]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(HGB / ref)^exponent`.
- **Source aliases:** none; `HGB` is the common NONMEM / clinical-PK abbreviation.
- **Example models:** `Yamada_2025_zolbetuximab.R` (g/L, reference 118; exponent -0.374 on V1).
- **Notes:** Unit varies by paper (SI g/L, US g/dL; 1 g/dL = 10 g/L). The per-model `covariateData[[HGB]]$units` field is load-bearing.

### WBC (**canonical for white blood cell count**)
- **Description:** Total white blood cell count (baseline or time-varying). In chronic lymphocytic leukaemia (CLL) populations the value is elevated because circulating leukaemic B-cells make up the majority of the count, so WBC can serve as a biomarker of target-cell burden rather than general hematology.
- **Units:** 10^9 cells/L (equivalent to 10^3 cells/uL). Document per-model via `covariateData[[WBC]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(WBC / ref)^exponent`. Reference values observed: 10 x 10^9/L (Mould 2007, typical CLL Vmax normalization).
- **Source aliases:** none known (`WBC` is the universal clinical-PK abbreviation).
- **Example models:** `Mould_2007_alemtuzumab.R` (reference 10 x 10^9/L; exponent 0.194 on Vmax).
- **Notes:** Time-varying in treatment studies where the drug depletes the leukaemic clone (e.g., alemtuzumab in CLL): WBC must be supplied at every observation time in the event dataset. In diseases where WBC is not therapeutically targeted the column can be treated as a baseline-only covariate; record the per-model convention in `covariateData[[WBC]]$notes`.

### NLR (**canonical for neutrophil-to-lymphocyte ratio**)
- **Description:** Ratio of absolute neutrophil count to absolute lymphocyte count from a complete blood count with differential. Used as a peripheral inflammation marker. May be reported as baseline only or as a time-varying covariate.
- **Units:** ratio (unitless)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(NLR / ref)^exponent` or exponential effects. Reference values observed: 2.11 (Lin 2024, median in pooled COVID-19 + non-infected cohort).
- **Source aliases:** none; `NLR` is the universal abbreviation in clinical-PK and inflammation-biomarker literature.
- **Example models:** `Lin_2024_casirivimab.R` (time-varying; reference 2.11; small positive exponent +0.029 on CL).
- **Notes:** Document baseline-vs-time-varying status in `covariateData[[NLR]]$notes`. Although it derives from `WBC` differential counts, register it as its own canonical because the ratio (not the absolute counts) is what the model uses.

### HCT (**canonical for hematocrit**)
- **Description:** Hematocrit -- packed red blood cell volume fraction (baseline or time-varying).
- **Units:** % (volume fraction times 100). Document per-model via `covariateData[[HCT]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(HCT / ref)^exponent`. Reference values observed: 45 % (Nestorov 2014, study-population median for severe hemophilia A adults).
- **Source aliases:** none; `HCT` is the universal NONMEM / clinical-PK abbreviation.
- **Example models:** `Nestorov_2014_factorviii.R` (reference 45 %, exponent -0.419 on V1).
- **Notes:** Higher HCT (more red-cell volume) leaves a smaller plasma fraction within total body volume; for plasma-restricted distribution (e.g., factor VIII activity, which circulates in plasma) the central volume of distribution decreases as HCT rises, so the exponent is negative. Document baseline-vs-time-varying status in `covariateData[[HCT]]$notes`. Distinct from `HGB` (mass concentration of hemoglobin); the two correlate but enter different mechanistic relationships.

### THB_MASS (**canonical for total hemoglobin mass**)
- **Description:** Subject-level total hemoglobin mass in grams. Plasma-volume-independent quantity measured by the optimised CO-rebreathing method (Schmidt 2005, Pottgiesser 2008); distinct from `HGB` (mass concentration in plasma) and `HCT` (volume fraction). Used in erythropoiesis / RBC-regeneration models as the steady-state set point Base that drives the negative-feedback term and seeds the steady-state initial conditions of the precursor compartments.
- **Units:** g
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- subject-level baseline supplied as a covariate column. Reference values observed: 885.42 g (Tetschke 2018 Table 1 example subject); Pottgiesser 2008 cohort mean ~870 g across 28 estimable adult-male volunteers.
- **Source aliases:** `Base` (Tetschke 2018 paper symbol).
- **Example models:** `Tetschke_2018_erythropoiesis.R` (reference 885.42 g; Pottgiesser 2008 dataset of 29 healthy adult male volunteers).
- **Notes:** Specific scope because total hemoglobin mass requires the optimised CO-rebreathing method to obtain (Schmidt 2005), which is a specialised technique not present in routine clinical labs; promote to `general` if a second model registers this quantity. Distinct from `HGB` (g/L or g/dL plasma concentration) and `HCT` (RBC volume fraction): `THB_MASS` is the absolute body-pool mass and is not perturbed by short-term plasma-volume fluctuations (Pottgiesser 2008 Section 3.2 explicitly motivates the choice of mass over concentration). Sex-dimorphic: typical value in adult males is meaningfully higher than in adult females; document the sex composition of the population in `covariateData[[THB_MASS]]$notes`.

### NEUT (**canonical for absolute neutrophil count**)
- **Description:** Absolute neutrophil count, typically as a baseline covariate (entered via centred-deviation `(NEUT - ref)` or power scaling `(NEUT / ref)^exponent`) or, in semi-mechanistic myelosuppression models, as a per-subject initial-condition value for the proliferation, transit, and circulating compartments.
- **Units:** cells/mm^3 (equivalent to cells/uL; i.e., the same value reported in 10^9/L x 1000). Document per-model via `covariateData[[NEUT]]$units` if the source paper uses a different unit (e.g., `10^9 cells/L` for `Ozawa_2007_docetaxel.R` per the paper's Table-3 reporting unit).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used in centred-deviation form `exp(coef * (NEUT - ref))`, in power scaling `(NEUT / ref)^exponent`, or as a direct per-subject initial-condition assignment in semi-mechanistic Friberg-family models. Reference values observed: 4133 cells/mm^3 (BAST PTTE 2017 simulated cohort median; `NA_NA_tte_gompertz.R` Event 1 base hazard model); 5 x 10^9/L (Ozawa 2007 typical Japanese cancer cohort, used as the initial condition for the proliferation, transit, and circulation compartments).
- **Source aliases:**
  - `BASE` -- per-subject baseline ANC supplied as a NONMEM data column (used in `Ozawa_2007_docetaxel.R`; Appendix I $INPUT).
- **Example models:** `NA_NA_tte_gompertz.R` (BAST PTTE 2017 / DDMODEL00000243 Event 1 hazard model; centred at NEUT = 4133/mm^3; coefficient -1.56e-4 on the NONMEM rescaled scale, equivalent to `exp(-1.56e-4 * (NEUT - 4133))` on the hazard), `Ozawa_2007_docetaxel.R` (Friberg-extension myelosuppression PD; per-subject baseline ANC supplied via the `NEUT` column, used as the initial condition for the proliferating, three transit, and circulating compartments per the Methods text 'Circ (t = 0) was fixed at its observed value').
- **Notes:** General scope because absolute neutrophil count is a routine clinical-laboratory measurement that recurs across cytotoxic-chemotherapy myelosuppression models (centred-deviation hazard models, Friberg-family per-subject baseline initial conditions, and time-varying ANC outputs). Promoted to general scope on 2026-05-10 to support `Ozawa_2007_docetaxel.R`. The NEUT canonical units are cells/mm^3, but the reporting unit `10^9 cells/L` (numerically NEUT_per_mm3 / 1000) is also common in oncology papers; per-model `covariateData[[NEUT]]$units` documents the per-paper unit. Distinct from `WBC` (total white blood cell count, of which neutrophils are the largest fraction in healthy adults) -- `NEUT` is a specific differential-count subfraction. Also distinct from `NLR` (neutrophil-to-lymphocyte ratio), which is a derived ratio.

## Coagulation / hemostasis biomarkers

### INR_BASE (**canonical for baseline international normalized ratio**)
- **Description:** Pre-medication baseline INR (international normalized ratio of prothrombin time). Time-fixed per subject (measured once, before the first warfarin dose). Used directly in the warfarin K-PD INR equation as an additive constant (`INR = INR_BASE + inrmax * (1 - (coag_s3 + coag_l3)/2)` per Xia 2024 supplement Section 1.1) so the simulated INR returns to the subject-specific baseline when the drug is removed.
- **Units:** (unitless ratio; INR has no units)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- subject-specific baseline. Default simulation value documented per-model in `covariateData[[INR_BASE]]$notes`; the Xia 2024 simulation uses the total-cohort mean of 1.13 (Table 1).
- **Source aliases:**
  - `INR_BASE`, `BL_INR`, `INRBASE` -- pre-medication INR column in NONMEM data sets; document the source-column name per-model in `covariateData[[INR_BASE]]$source_name`.
- **Example models:** `Xia_2024_warfarin.R` (additive baseline in the INR observation equation; cohort mean 1.13, SD 0.59 per Xia 2024 Table 1).
- **Notes:** Distinct from a time-varying INR observation (the model's observed `INR` variable). Healthy subjects with no anticoagulation typically have INR around 1.0; the Hamberg / Xia 2024 model treats deviations from 1.0 as a subject-specific covariate rather than an estimated parameter so the model returns to the observed baseline when warfarin is withdrawn. Ratified canonically on 2026-05-16 alongside the Xia 2024 warfarin extraction.

### VWF (**canonical for von Willebrand factor concentration**)
- **Description:** Plasma concentration (or activity) of von Willebrand factor (VWF) -- the multimeric carrier protein that binds and protects circulating factor VIII (FVIII) from proteolytic degradation and rapid clearance. Used as a covariate on FVIII (and FVIII-Fc) clearance because the vast majority (>95%) of circulating FVIII is in complex with VWF.
- **Units:** IU/dL (equivalent to % of pooled normal plasma); document per-model via `covariateData[[VWF]]$units`. Some sources report `VWF:Ag` (antigen) versus `VWF:RCo` (ristocetin cofactor activity); record which assay was used in `covariateData[[VWF]]$notes`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(VWF / ref)^exponent`. Reference values observed: 118 IU/dL (Nestorov 2014, study-population median).
- **Source aliases:** none; `VWF` is the universal abbreviation. Source papers may write `vWF` (lowercase v) or specify the assay (`VWF:Ag`).
- **Example models:** `Nestorov_2014_factorviii.R` (reference 118 IU/dL, exponent -0.343 on CL; VWF antigen).
- **Notes:** Higher VWF protects FVIII from clearance, so the exponent on CL is negative. VWF is time-varying within an individual (acute-phase response, age, blood group, etc.), but most published population PK models use baseline-only VWF when the within-subject dynamics are not characterized; document the per-model convention in `covariateData[[VWF]]$notes`.

### DDIMER (**canonical for plasma D-dimer concentration**)
- **Description:** Plasma D-dimer protein concentration, the fibrin-degradation peptide produced by plasmin-mediated cleavage of cross-linked fibrin. Used in vascular / coagulation-pathology models as a circulating biomarker of fibrin turnover, intra-aneurysmal thrombus burden, or systemic fibrinolytic activity.
- **Units:** ng/mL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used either with log10-transformed proportional scaling `log10(DDIMER) / median(log10(DDIMER))` (Sherer 2012) or with categorical strata (Sherer 2012 sensitivity analysis groups: <=150, 151-300, 301-900, >900 ng/mL). Reference values observed: 326 ng/mL (Sherer 2012 cohort median; log10 approx 2.513).
- **Source aliases:**
  - `C^(D-dimer)` -- used in `Sherer_2012_AAA.R` (the symbol in Sherer 2012 Methods equation page 2).
- **Example models:** `Sherer_2012_AAA.R` (proportional log10-transformed covariate on the baseline AAA growth rate `beta1` (`e_ddimer_b1 = 0.90 mm/year`) and on the first derivative of growth rate with size `beta2` (`e_ddimer_b2 = 0.37/year`)).
- **Notes:** Specific scope until a second model registers the canonical. Time-fixed (baseline-only) in Sherer 2012 because the HIMS cohort had a single follow-up D-dimer measurement; the source paper flags this as a limitation. Cohort interquartile range 142-785 ng/mL; extrapolation outside this range is not validated by the source. The log10 transformation reflects Sherer 2012's finding that "differences in AAA growth were predominantly driven by patients with the highest plasma D-dimer concentrations." Distinct from the time-varying biomarker columns in indirect-response / TMDD models -- DDIMER enters Sherer 2012 as a baseline regression covariate, not as a dynamic exposure / response variable. Ratified canonically on 2026-05-16 alongside the Sherer 2012 extraction.

### AAA_DIAM (**canonical for baseline abdominal aortic aneurysm diameter**)
- **Description:** Abdominal aortic aneurysm (AAA) maximum infrarenal diameter, ascertained by ultrasound at study entry. Used in vascular disease-progression models as the per-subject baseline severity covariate that anchors the typical-value regression for individual-level growth parameters.
- **Units:** mm
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used in proportional form `AAA_DIAM / median(AAA_DIAM)` so the effect coefficients represent the contribution at the cohort median. Reference value observed: 32.7 mm (Sherer 2012 cohort median; q1 30.8, q3 36.0).
- **Source aliases:**
  - `Y(0)` -- used in `Sherer_2012_AAA.R` (the symbol in Sherer 2012 Methods equation page 2; the baseline screening ultrasound diameter).
- **Example models:** `Sherer_2012_AAA.R` (proportional covariate on all three individual-level parameters: `e_aaadiam_b0 = 32.6 mm` on baseline size beta0, `e_aaadiam_b1 = 2.03 mm/year` on baseline growth rate beta1, and `e_aaadiam_b2 = 0.59/year` on the first derivative of growth rate with size beta2).
- **Notes:** Specific scope until a second model registers the canonical. Time-fixed (baseline-only) per subject -- the value is the single screening ultrasound diameter; the time-evolving AAA diameter during follow-up is the model's observation, not the covariate. Sherer 2012 inclusion criterion (HIMS cohort): 30-49 mm small AAA, so extrapolation outside this range to <30 mm (non-aneurysmal aorta) or >=50 mm (surgical-referral threshold) is not validated by the source. Ratified canonically on 2026-05-16 alongside the Sherer 2012 extraction.

## Disease severity scores

### EASI (**canonical for Eczema Area and Severity Index**)
- **Description:** Eczema Area and Severity Index score (atopic-dermatitis severity composite; bounded continuous, scale 0-72 with higher values = more severe disease).
- **Units:** (score)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- healthy volunteers have EASI = 0. Effect enters as an additive term in models that pool AD patients with HV (e.g., `Tiraboschi_2025_amlitelimab.R`).
- **Source aliases:**
  - `BEASI` (baseline EASI) -- used in `Tiraboschi_2025_amlitelimab.R`.
- **Example models:** `Tiraboschi_2025_amlitelimab.R`.
- **Notes:** When used as a time-invariant baseline covariate (`BEASI`), document in `covariateData[[EASI]]$notes`. Canonical name is `EASI` without the `B` prefix to match the `AGE` / `WT` / `ALB` pattern where baseline vs time-varying status is recorded in notes rather than the column name.

### MGADL (**canonical for Myasthenia Gravis Activities of Daily Living score**)
- **Description:** Myasthenia Gravis Activities of Daily Living score -- eight-item patient-reported outcome measure (each item 0-3), total 0-24, higher values = greater symptom severity and functional limitation.
- **Units:** (score)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- healthy participants (no gMG) have `MGADL = 0` by definition. Effect enters as a baseline covariate on MG-ADL response parameters in gMG cohorts.
- **Source aliases:** none known.
- **Example models:** `Valenzuela_2025_nipocalimab.R` (reference 7 points; power-form effect on `IDecplacebo` and on the slope between MG-ADL change and IgG reduction).
- **Notes:** Baseline-only in Valenzuela 2025 (the observation is the absolute change from baseline MG-ADL). When used time-varying (e.g., in pure PD models driven by disease-progression dynamics), document in `covariateData[[MGADL]]$notes`. Canonical name is `MGADL` without a `BL` prefix to match the `EASI` / `AGE` / `WT` / `ALB` pattern.

### BCVA (**canonical for best-corrected visual acuity**)
- **Description:** Best-corrected visual acuity score measured on the Early Treatment Diabetic Retinopathy Study (ETDRS) chart, expressed as the number of letters read correctly (0-100; higher values = better vision). Used as a baseline severity covariate in ophthalmology PK/PD models of anti-VEGF treatment.
- **Units:** ETDRS letters (0-100)
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used as a baseline input to set the initial condition of an indirect-response BCVA state or as a power-form effect on response parameters. Reference value observed: 55 letters (Mulyukov 2018 narrative: mean study-population baseline BCVA).
- **Source aliases:**
  - `BVA` (baseline visual acuity) -- used in `Mulyukov_2018_ranibizumab.R`.
- **Example models:** `Mulyukov_2018_ranibizumab.R` (baseline BCVA used as the center for the initial-condition draw `g0 = BCVA + eta_g0`).
- **Notes:** Ophthalmology-specific. Baseline-only in Mulyukov 2018 (carried once per subject and used only as the starting BCVA for the indirect-response model). Canonical name drops the `B` prefix to match the `EASI` / `AGE` / `WT` / `ALB` pattern (baseline-vs-time-varying status recorded in `covariateData[[BCVA]]$notes`). Scope is `specific` until a second ophthalmology model ratifies the name; at that point promote to `general`.

### PREV_AE_SCORE (**canonical for previous-time-step ordinal adverse-event score**)
- **Description:** Ordinal adverse-event score recorded at the immediately preceding observation time, used as a Markov-state covariate that conditions the current-time logit / probability calculation on the previous outcome. Integer 0..N where N is the maximum AE grade in the source-paper grading scheme; convention is that the value is 0 at the first observation (no prior AE).
- **Units:** (ordinal score; document per-model the scale in `covariateData[[PREV_AE_SCORE]]$notes`)
- **Type:** count
- **Scope:** specific
- **Reference category:** n/a -- typically used as a categorical conditioner (`IF (PREV_AE_SCORE == 0) ...`) or via piecewise-FPS indicator decomposition; the natural reference is `PREV_AE_SCORE = 0` (no prior AE).
- **Source aliases:**
  - `PREVSCOR` -- used in `Girard_2012_pimasertib.R` (CTCAE 0..3 ocular-AE score).
- **Example models:** `Girard_2012_pimasertib.R` (Markov-state covariate that selects per-previous-score logit thresholds `b01`/`b11`/`b21` and `b02`/`b12`/`b22` and per-previous-score `emax` levels; reset to 0 at TIME = 0 per source `IF (TIME.EQ.0) PREVSCOR=0`).
- **Notes:** Specific scope because the ordinal scale is paper-specific (different AE-grading schemes, different number of categories, different grouping rules -- Girard 2012 collapses CTCAE grades 1+2 into a single "1-2" category and treats grades >=3 as a third stratum). When assembling the simulation event table, set `PREV_AE_SCORE = 0` at the first observation of every subject and update each subsequent observation to the previous observation's sampled score -- matching the NONMEM `IF (TIME.EQ.0) PREVSCOR=0` / `PREVSCOR = DV` carry-forward idiom. Distinct from `PAIN` (continuous baseline pain score) and from `MGADL` / `EASI` (continuous severity scores not modelled as Markov states).

### CDR_SOB (**canonical for Clinical Dementia Rating - Sum of Boxes score**)
- **Description:** Clinical Dementia Rating scale - Sum of Boxes (CDR-SOB) score, the unweighted sum of the six CDR-box scores (memory, orientation, judgement and problem solving, community affairs, home and hobbies, personal care; each scored 0 / 0.5 / 1 / 2 / 3). Total ranges 0-18; higher values indicate more severe cognitive and functional impairment. Widely used as a primary efficacy endpoint in Alzheimer's-disease (AD) clinical trials and as a disease-progression biomarker in mild-cognitive-impairment (MCI) / prodromal-AD populations.
- **Units:** (CDR-SOB units, 0-18 score)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- the source paper centres covariate effects on dataset medians (e.g., Delor 2013 uses CDR_SOB / 2 in the DOT power-form and CDR_SOB - 1 in the mixture-logit additive form).
- **Source aliases:**
  - `CDR_bsl` -- used in `Delor_2013_alzheimer.R` (baseline CDR-SOB at study entry).
  - `CDR` -- alternative bare-name often seen in ADNI / CAMD-style NONMEM datasets.
- **Example models:** `Delor_2013_alzheimer.R` (time-fixed baseline covariate; enters both the per-subject DOT power form (`(CDR_SOB / 2)^e_cdr_sob_dot` with `e_cdr_sob_dot = -0.072`) and the per-subject slow-progression mixture-logit additive form (`+ e_cdr_sob_slow * (CDR_SOB - 1)` with `e_cdr_sob_slow = -1.27`)).
- **Notes:** Canonical name is `CDR_SOB` without a `_BL` suffix to match the `EASI` / `MGADL` / `BCVA` pattern (baseline-vs-time-varying status recorded in `covariateData[[CDR_SOB]]$notes`). The CDR sum-of-boxes form is distinct from the global CDR rating (`CDR_GLOBAL`, a 0 / 0.5 / 1 / 2 / 3 ordinal); the sum-of-boxes is preferred in disease-progression modelling for its finer granularity. Ratified canonically on 2026-05-16 alongside the Delor 2013 extraction.

### ADAS_COG (**canonical for ADAS-cog total cognitive subscale score**)
- **Description:** Alzheimer's Disease Assessment Scale - cognitive subscale (ADAS-cog) total score. The total-11 form ranges 0-70; the modernised total-13 form ranges 0-85. Higher values = more cognitive impairment. The ADAS-cog is the most widely used cognitive endpoint in AD clinical trials and disease-progression modelling.
- **Units:** (ADAS-cog units; document the form (total-11 / total-13) per-model in `covariateData[[ADAS_COG]]$units`)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- covariate effects typically centred on a dataset median (e.g., Delor 2013 centres on ADAS_COG = 12.67).
- **Source aliases:**
  - `ADAS_bsl` -- used in `Delor_2013_alzheimer.R` (baseline ADAS-cog total-11 at study entry).
  - `ADAS`, `ADAS_COG_11`, `ADAS_COG_13` -- alternative bare-name forms seen across ADNI / CAMD datasets.
- **Example models:** `Delor_2013_alzheimer.R` (time-fixed baseline covariate; ADAS-cog total-11 form; enters the per-subject DOT power form `(ADAS_COG / 12.67)^e_adas_cog_dot` with `e_adas_cog_dot = -0.0439`).
- **Notes:** Canonical name is `ADAS_COG` (no `_BL` suffix; baseline-vs-time-varying recorded in notes). The total-11 vs total-13 form must be documented per-model in `covariateData[[ADAS_COG]]$units` because the same numeric ADAS_COG value has different clinical interpretation across the two forms. Conrado 2014 uses ADAS-cog as the modelled observation (response variable) rather than as a baseline covariate; that model file therefore does not list ADAS_COG in its `covariateData`. Ratified canonically on 2026-05-16 alongside the Delor 2013 extraction.

### MMSE (**canonical for Mini Mental State Examination score**)
- **Description:** Mini Mental State Examination total score (0-30; higher values = better cognitive function). A widely used cognitive-screening instrument in AD and MCI populations.
- **Units:** (MMSE units, 0-30 score)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- covariate effects typically centred on a dataset median (e.g., Delor 2013 centres on MMSE = 26).
- **Source aliases:**
  - `MMSE_bsl` -- used in `Delor_2013_alzheimer.R` (baseline MMSE at study entry).
- **Example models:** `Delor_2013_alzheimer.R` (time-fixed baseline covariate; modifies the per-subject disease-progression acceleration parameter alpha via a power form `(MMSE / 26)^e_mmse_alpha` with `e_mmse_alpha = -2.01`).
- **Notes:** Canonical name is `MMSE` (no `_BL` suffix; baseline-vs-time-varying recorded in notes). MMSE is the inverse-direction counterpart of CDR_SOB / ADAS_COG (MMSE high = healthy; CDR_SOB / ADAS_COG high = impaired); covariate-effect coefficient signs are therefore typically opposite to those for CDR_SOB / ADAS_COG. Ratified canonically on 2026-05-16 alongside the Delor 2013 extraction.

### FAQ (**canonical for Functional Assessment Questionnaire score**)
- **Description:** Functional Assessment Questionnaire (Pfeffer FAQ) total score: sum of ten functional-activities items (each scored 0 = normal to 3 = dependent), total 0-30; higher values = greater functional impairment. Used as a functional-status covariate alongside cognitive scores in AD and MCI populations.
- **Units:** (FAQ units, 0-30 score)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- covariate effects typically centred on a dataset median (e.g., Delor 2013 centres on FAQ = 1).
- **Source aliases:**
  - `FAQ_bsl` -- used in `Delor_2013_alzheimer.R` (baseline FAQ at study entry).
- **Example models:** `Delor_2013_alzheimer.R` (time-fixed baseline covariate; enters the per-subject slow-progression mixture-logit additive form `+ e_faq_slow * (FAQ - 1)` with `e_faq_slow = -0.341`).
- **Notes:** Canonical name is `FAQ` (no `_BL` suffix; baseline-vs-time-varying recorded in notes). Distinct from the cognitive scores (CDR_SOB / ADAS_COG / MMSE): FAQ measures instrumental activities of daily living rather than cognitive performance, and adds incremental information about disease-stage severity in MCI cohorts. Ratified canonically on 2026-05-16 alongside the Delor 2013 extraction.

### RHPNM (**canonical for normalized hippocampal volume (head-size and age-adjusted)**)
- **Description:** Normalised hippocampal volume: the subject's average left+right hippocampal volume divided by the value expected for a healthy subject of the same age and estimated intracranial volume (head size). 1.0 corresponds to the healthy reference; values below 1.0 indicate hippocampal atrophy. Derived from MRI volumetry and a healthy-subject regression on age and intracranial volume.
- **Units:** (unitless ratio; 1.0 = healthy reference)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- covariate effects typically centred on RHPNM = 1 (the healthy reference, e.g., Delor 2013).
- **Source aliases:**
  - `RHPNM` -- used in `Delor_2013_alzheimer.R` (baseline normalized hippocampal volume; Delor 2013 derivation: `RHPNMbsl_i = HIPVbsl_i / HPNMbsl_i` where `HPNMbsl_i = Age_i * (-26.6268 + EICVbsl_i * 0.0016 + 3340.4395)`).
- **Example models:** `Delor_2013_alzheimer.R` (time-fixed baseline covariate; enters the per-subject slow-progression mixture-logit additive form `+ e_rhpnm_slow * (RHPNM - 1)` with `e_rhpnm_slow = 7.5`, a strongly positive effect indicating that less atrophic hippocampi (RHPNM closer to 1) are associated with a higher probability of being in the slow-progressing subpopulation).
- **Notes:** Distinct from the raw hippocampal volume (which would be a `HIPV` canonical not yet registered; raw HIPV is confounded with head size and age, hence the need for the normalisation). The Delor 2013 paper notes that the same effect is only marginally significant with unnormalised HIPV (P = 0.02) but strongly significant with the normalised form (RHPNM). The exact age / EICV regression coefficients are paper-specific and any future model adopting this canonical should re-derive the normalisation for its own population or document why the Delor 2013 regression is reused. Ratified canonically on 2026-05-16 alongside the Delor 2013 extraction.

### ACUTE_MED_DAYS (**canonical for baseline number of days/month of acute migraine medication use**)
- **Description:** Baseline number of days per month on which acute migraine medication (triptans or ergot compounds) was used during the 28-day run-in period prior to first dose. Enters as a piecewise-linear shift on baseline migraine or moderate-to-severe headache days in migraine exposure-response models.
- **Units:** days/month
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- piecewise-linear shift with breakpoint at 5 d/mo: contributes 0 below 5 and `slope * (ACUTE_MED_DAYS - 5)` above 5 (Fiedler-Kelly 2020). The 5-day breakpoint reflects the clinical guideline for medication-overuse headache.
- **Source aliases:** "Baseline days/month of acute medications" -- used in `FiedlerKelly_2020_fremanezumab_em.R` and `FiedlerKelly_2020_fremanezumab_cm.R`.
- **Example models:** `FiedlerKelly_2020_fremanezumab_em.R` (slope 0.438 d/d, episodic migraine), `FiedlerKelly_2020_fremanezumab_cm.R` (slope 0.460 d/d, chronic migraine).
- **Notes:** Specific scope because the variable is migraine-domain-bound. Time-fixed per subject (baseline-only). When future migraine E-R models register additional aliases or alternative breakpoints, document them per-model and consider promoting to `general`.

### AMLOAD (**canonical for systemic-amyloidosis whole-body amyloid load grade**)
- **Description:** Ordinal whole-body amyloid load score in patients with systemic amyloidosis. Integer 0-3: 0 = no amyloid (healthy volunteers), 1 = small amyloid load, 2 = moderate amyloid load, 3 = large amyloid load. The grading combines organ-by-organ amyloid presence (liver, spleen, bone, adrenals, gut, heart) and SAP-scintigraphy uptake into a single per-patient severity grade at baseline. Time-fixed per subject.
- **Units:** (categorical 0-3)
- **Type:** categorical
- **Scope:** specific
- **Reference category:** 0 (no amyloid). Sahota 2015 Eq. 2 also treats category 1 as part of the reference for the V4 effect (categories 0 and 1 share a V4 multiplier of 1); only categories 2 and 3 carry a non-zero effect via the cumulative parameters `e_amload2_vp_sap` and `e_amload2_vp_sap + e_amload3_vp_sap` respectively.
- **Source aliases:** none beyond the source-paper `AMLOAD` column.
- **Example models:** `NA_NA_miridesap.R` (DDMODEL00000262; Sahota 2015 Eq. 2 multiplicative effect on SAP peripheral volume V4: V4 = V4_ref * (1 + e_amload2_vp_sap * I(AMLOAD>=2) + e_amload3_vp_sap * I(AMLOAD>=3)); reported effects e_amload2_vp_sap = 6.39 / e_amload3_vp_sap = 26.39 yielding ~7.4x V4 at moderate load and ~33.8x at large load), `Sahota_2015_miridesap.R` (paper-only extraction of the same Sahota 2015 final model with identical Eq. 2 effect on V4; values 6.39 / 26.39 taken from Table 2).
- **Notes:** Scope: specific because the grading scheme is amyloidosis-specific (Sahota 2015 Methods: "The whole body amyloid load covariate, AMLOAD, was a categorical score: 0 for no amyloid in healthy volunteers, 1 for small, 2 for moderate, and 3 for large"). The cumulative monotonic parameterisation in Sahota 2015 Eq. 2 encodes a positive-only step at each grade increment. Co-used with `AMLIVER` for the binary hepatic-involvement modifier. Ratified canonically on 2026-05-15 alongside the DDMODEL00000262 / Sahota 2015 extraction.

### AMLIVER (**canonical for hepatic amyloid involvement indicator**)
- **Description:** Binary indicator for the presence of amyloid in the liver as a separate organ involvement from the overall whole-body amyloid load. 1 = liver amyloid present at baseline; 0 = no liver amyloid. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no liver amyloid).
- **Source aliases:** none beyond the source-paper `AMLIVER` column.
- **Example models:** `NA_NA_miridesap.R` (DDMODEL00000262; Sahota 2015 Eq. 2 multiplicative effect on SAP intercompartmental clearance Q4: Q4 = Q4_ref * (1 + e_amliver_q4 * AMLIVER); reported effect 4.01, yielding ~5x Q4 in patients with hepatic amyloid), `Sahota_2015_miridesap.R` (paper-only extraction of the same Sahota 2015 final model with identical Eq. 2 effect on Q4; value 4.01 from Table 2).
- **Notes:** Scope: specific because the covariate is amyloidosis-specific (Sahota 2015 Methods names AMLIVER alongside AMSPLEEN / AMHEART as organ-specific amyloid-involvement binary indicators; only AMLIVER was retained as a covariate in the final model). Used in combination with the global `AMLOAD` grade so the model can express both general amyloid burden and the specific hepatic-clearance modifier (the SAP-CPHPC complex is cleared by the liver, motivating the hepatic-amyloid-specific Q4 effect). Ratified canonically on 2026-05-15 alongside the DDMODEL00000262 / Sahota 2015 extraction.

## Critical-illness severity

### ORG_FAIL_COUNT (**canonical for number of organs failing in critically ill patients**)
- **Description:** Integer count of failing organs in a critically ill patient at a given observation day, ascertained per-day and reported as the worst-of-day count. The count is decomposed into mutually exclusive strata 0 / 1 / 2 / 3 / >=4 (Vet 2016 used the strata 0, 1, 2, 3, and 4-or-5) that select per-stratum typical clearance values; the strata are not collapsed onto a single linear or power covariate effect because the underlying organ-failure mechanisms (cardiovascular, pulmonary, renal, hepatic, neurologic, hematologic) impair drug elimination heterogeneously.
- **Units:** (count)
- **Type:** categorical
- **Scope:** specific
- **Reference category:** 0 (no organs failing). Per-stratum typical CL values are estimated for ORG_FAIL_COUNT = 1, 2, 3, and >=4; ORG_FAIL_COUNT = 0 sets the baseline typical CL (frequently FIXED at the baseline value, as in Vet 2016).
- **Source aliases:**
  - `ORGF` -- used in `Vet_2016_midazolam.R` (DDMODEL00000249 NMTRAN `$INPUT` column; values 0..>=4). Renamed to canonical `ORG_FAIL_COUNT` when assembling input data for the packaged model.
- **Example models:** `Vet_2016_midazolam.R` (per-stratum typical CL values: ORG_FAIL_COUNT=0 fixed at 1.6 L/h for a 5 kg child with CRP=32 mg/L; ORG_FAIL_COUNT=1 -> 1.29 L/h; ORG_FAIL_COUNT=2 -> 0.957 L/h; ORG_FAIL_COUNT=3 -> 0.842 L/h; ORG_FAIL_COUNT>=4 -> 0.678 L/h).
- **Notes:** Specific scope because the variable is critical-care-population-bound (PICU / ICU). Time-varying within subject (re-evaluated each ICU day). The Vet 2016 organ-failure ascertainment follows the Wilkinson 1987 paediatric multiple organ system failure (MOSF) criteria -- operator-confirmed per-paper criteria should be documented in each model's `covariateData[["ORG_FAIL_COUNT"]]$notes`. Decompose inside `model()` into binary indicators (`orgf1 <- (ORG_FAIL_COUNT == 1)`, `orgf2 <- (ORG_FAIL_COUNT == 2)`, `orgf3 <- (ORG_FAIL_COUNT == 3)`, `orgf_ge4 <- (ORG_FAIL_COUNT >= 4)`) and select per-stratum CL with mutually-exclusive multiplicative-flag arithmetic. Ratified canonically on 2026-05-06.

## Interferon / biomarker panels

### BGENE21 (**canonical for 21-gene type I interferon signature score**)
- **Description:** Baseline 21-gene type I interferon signature score -- a composite transcriptomic score summarising the expression of 21 interferon-regulated genes in whole blood relative to a healthy-donor reference, used as a biomarker of type I IFN pathway activation in SLE and related autoimmune conditions.
- **Units:** unitless fold-change score (relative to healthy-donor reference).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with power scaling `(BGENE21 / ref)^exponent`. Reference values observed: 32 in Narwal 2013 (study-population median was 33), 12.04 in Zheng 2016 (median of the SLE phase IIb cohort, range 0.32-38.59).
- **Source aliases:** none.
- **Example models:** `Narwal_2013_sifalimumab.R` (reference 32, exponent 0.0558 on CL), `Zheng_2016_sifalimumab.R` (reference 12.04, power effect on CL with exponent 0.09).
- **Notes:** Specific to drugs whose mechanism targets the type I IFN pathway (e.g., anti-IFN-alpha antibodies like sifalimumab, anifrolumab). Higher BGENE21 indicates stronger target engagement / disease activity and is associated with increased drug clearance via target-mediated mechanisms. The 21-gene panel composition is tied to the MedImmune/AstraZeneca SLE development programme; a different IFN gene signature (e.g., a 4-gene or 5-gene panel) should be registered under its own canonical name (`BGENE4`, `IFN_SIG`, ...) to avoid conflating panel definitions.

### BGENE21_HIGH (**canonical for binary high-vs-low IFN-21-gene indicator**)
- **Description:** 1 = subject's baseline 21-gene type I IFN signature score is at or above the paper-specified high/low cut-off, 0 = below cut-off.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (low).
- **Source aliases:** none.
- **Example models:** `Almquist_2022_anifrolumab.R` (binary high-IFN indicator on CL).
- **Notes:** Pair with continuous `BGENE21` when the paper reports both. The high/low cut-off is paper-specific (commonly the population median) and must be documented in `covariateData[[BGENE21_HIGH]]$notes` for every model that uses this covariate. Operator decision (2026-04-28): use `BGENE21_HIGH` (not `IFNGS_HIGH`) so the link to the existing `BGENE21` register entry is explicit while the binary nature stays visible in the column name.

## Inflammation markers

### EOS (**canonical for blood eosinophil count**)
- **Description:** Blood eosinophil count (baseline or time-varying).
- **Units:** cells/uL
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(EOS / ref)^exponent`.
- **Source aliases:**
  - `BEOS` (baseline EOS) -- used in `Kotani_2022_astegolimab.R`.
- **Example models:** `Kotani_2022_astegolimab.R` (reference 180 cells/uL, baseline).
- **Notes:** Used as a surrogate of inflammatory burden that correlates with protein turnover and therefore mAb clearance. Canonical name drops the `B` prefix to match the `EASI` / `AGE` / `WT` / `ALB` pattern; baseline-vs-time-varying status is documented in `covariateData[[EOS]]$notes`.

### BLBCELL (**canonical for baseline CD19+ B cell count**)
- **Description:** Baseline CD19+ B cell count (cells/uL) measured by fluorescence-activated cell sorting (FACS) prior to first dose. Used as a covariate / scaling biomarker for B-cell-targeted antibody PK-PD models (e.g., anti-CD20 mAbs in multiple sclerosis or B cell malignancies).
- **Units:** cells/uL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with power scaling `(BLBCELL / ref)^exponent`. Reference value observed: 200 cells/uL (Yu 2022, median of the pooled five-study cohort).
- **Source aliases:**
  - `Bcell0` -- used in `Yu_2022_ofatumumab.R`.
  - `BBCC` (NHL Phase I/Ib/II convention; values in 10^6 cells/L = cells/uL) -- used in `Lu_2019_polatuzumab.R`.
- **Example models:** `Yu_2022_ofatumumab.R` (power effect on the maximum B-cell-lysis stimulatory effect Emax, exponent 0.275, reference 200 cells/uL), `Lu_2019_polatuzumab.R` (two distinct effects: power on CL_INF with input floored at 1 cell/uL, and a thresholded power on CL_T with the BLBCELL/121-cells/uL ratio floored at 1).
- **Notes:** Distinct from a *time-varying* B cell count, which is the PD response variable rather than a covariate. Scope: specific because the clinically relevant baseline depends on the surface marker (CD19, CD20, CD22) and whether the panel reports total B cells or memory/naive subsets -- register a new canonical name if a future paper uses a different marker. Both Yu 2022 (anti-CD20 ofatumumab) and Lu 2019 (anti-CD79b polatuzumab vedotin) use CD19+ counts, so the canonical is reused; subtype-specific differences are documented in each model's `covariateData[[BLBCELL]]$notes`.

### BL_PARP_PBL (**canonical for baseline poly(ADP-ribose) polymerase activity in peripheral blood lymphocytes**)
- **Description:** Subject-specific baseline (pre-dose) poly(ADP-ribose) polymerase (PARP) activity in peripheral blood lymphocytes (PBL), measured as picomoles of PAR polymer formed per 10^6 PBL by an enzyme-activity assay. Used as a covariate / scaling biomarker on the maximal-inhibition residual activity (Emin) parameter in PARP-inhibitor PK/PD models.
- **Units:** pmol/10^6 PBL (document per-model via `covariateData[[BL_PARP_PBL]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with power scaling `(BL_PARP_PBL / ref)^exponent`. Reference value observed: 90.8 pmol/10^6 PBL (Wang 2015 rucaparib; population typical baseline E0 used as a stand-in for the unreported study-cohort median because the paper reports only the typical E0 and the exponent value, not the numeric BLB_median).
- **Source aliases:**
  - `BLB` (Wang 2015's notation for "baseline level in blood / PBL") -- used in `Wang_2015_rucaparib.R`.
- **Example models:** `Wang_2015_rucaparib.R` (power effect on residual maximum-inhibition parameter Emin; exponent 0.620 with the form `Emin = TV(Emin) * (BL_PARP_PBL / BLB_median)^alpha`; PBL paired with a separate tumor-tissue PARP activity covariate that is not yet a registered canonical because the units differ -- pmol/mg protein for tumor vs pmol/10^6 PBL for blood).
- **Notes:** Specific scope because the column is meaningful only for PARP-inhibitor PK/PD models (rucaparib, olaparib, niraparib, talazoparib, veliparib, etc.) and the units are tied to the PBL-specific assay format. A future tumor-tissue PARP activity covariate would need a separate canonical because the units differ (pmol/mg protein) and the biology of cellular PARP activity per mg of protein is not numerically interchangeable with PBL-normalized PARP activity. The Wang 2015 model uses BL_PARP_PBL only on Emin (residual maximum inhibition) and not on E0 or IC50; per-paper effects must be documented in each model's `covariateData[[BL_PARP_PBL]]$notes`. The paper does not publish the numeric study-cohort median of BLB used to center the covariate; the model file uses 90.8 pmol/10^6 PBL (the population typical baseline E0 reported in Wang 2015 Table 2) as a defensible default reference and documents the assumption in the vignette's Assumptions and deviations section.

### CSF1 (**canonical for colony-stimulating factor 1 / macrophage-colony-stimulating factor concentration**)
- **Description:** Plasma colony-stimulating factor 1 (CSF-1, also known as macrophage colony-stimulating factor, M-CSF) concentration (baseline or time-varying). The hematopoietic cytokine that signals through CSF-1R to drive monocyte / macrophage differentiation and survival; used as both a target-engagement biomarker (anti-CSF-1R mAbs increase circulating free CSF-1) and a baseline covariate on PK / PD parameters in CSF-1R-pathway PopPK/PD models.
- **Units:** pg/mL (= ng/L; the two labels are numerically equivalent). Document per-model via `covariateData[[CSF1]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with power scaling `(CSF1 / ref)^exponent`. Reference values observed: 549 pg/mL (Yang 2024 axatilimab; pooled-cohort median).
- **Source aliases:** none formally; informal aliases include `BLCSF1` (baseline CSF-1) and `BL_CSF1` (model-parameter notation in Monolix / NONMEM control streams).
- **Example models:** `Yang_2024_axatilimab.R` (baseline-only covariate on linear clearance `CL` with power exponent 0.912 and on the model parameter `BL_CSF1` with power exponent 0.656; reference 549 pg/mL).
- **Notes:** Specific scope because the column is meaningful only for CSF-1R-pathway-targeting drugs (axatilimab and future anti-CSF-1R molecules). Distinct from any CSF-1 model state representing time-course dynamics -- covariate column is the pre-dose laboratory observation, typically measured by an ELISA assay (Yang 2024 used the R&D Systems Quantikine ELISA). Per-model `covariateData[[CSF1]]$notes` should document the assay used and any LOQ-related imputation for samples below the assay's limit of detection.

### CRP (**canonical for C-reactive protein**)
- **Description:** C-reactive protein concentration. Covers both standard and high-sensitivity (hs-CRP) assays and both baseline and time-varying usages. Each model's `covariateData[[CRP]]$description` and `notes` must state the assay type (standard vs hs-CRP) and whether the column carries a baseline-only or time-varying value, including the paper-specific reference value used for power scaling.
- **Units:** mg/L (document per-model via `covariateData[[CRP]]$units`).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(CRP / ref)^exponent` or exponential effects `exp(coef * (CRP - ref))`. Reference values observed: 4.23 mg/L (Moein 2022, IBD standard assay), 4.31 mg/L (Moein 2022 Table 3 median), 5.21 mg/L (Thakre 2022, baseline hs-CRP), 7.41 mg/L (Chua 2025, baseline standard assay), 14.2 mg/L (Xu 2019, baseline standard assay), 15.7 mg/L (Ma 2020, baseline standard assay), 0.837 mg/dL = 8.37 mg/L (Wang 2020, IBD standard assay; the model carries the source unit mg/dL).
- **Source aliases:**
  - `hsCRP` -- high-sensitivity CRP (mixed-case preserved from earlier register drafts).
  - `HSCRP` -- all-caps variant.
  - `CRPHS` -- used in `Thakre_2022_risankizumab.R` (baseline, high-sensitivity assay).
  - `BLCRP` -- baseline CRP; used in `Xu_2019_sarilumab.R` and `Ma_2020_sarilumab_das28crp.R`.
- **Example models:** `Thakre_2022_risankizumab.R`, `Xu_2019_sarilumab.R`, `Chua_2025_mirikizumab.R`, `Moein_2022_etrolizumab.R`, `Ma_2020_sarilumab_das28crp.R`, `Wang_2020_ontamalimab.R` (mg/dL, reference 0.837).
- **Notes:** The prior separate `hsCRP`, `BLCRP`, and standard-assay `CRP` canonicals were merged on 2026-04-20 to a single general-scope `CRP` canonical. Assay type (standard vs hs-CRP), baseline-vs-time-varying status, and the paper-specific reference value all live in each model's `covariateData[[CRP]]$description` / `notes`. Only aggregate values from hs-validated assays as CRP when the downstream analysis relies on low-range sensitivity; for most inflammatory-disease cohorts (IBD, RA/PsA), baseline CRP is well above the hs-sensitivity range and the distinction is moot.

### AAG (**canonical for alpha-1 acid glycoprotein concentration**)
- **Description:** Serum alpha-1 acid glycoprotein (AAG; orosomucoid; ORM1) concentration, an acute-phase plasma glycoprotein that binds basic and lipophilic drugs (including taxanes such as docetaxel and paclitaxel). Elevated in cancer, inflammation, and infection; influences free-drug fraction and downstream PD effects in cytotoxic-chemotherapy myelosuppression models.
- **Units:** g/L (= mg/mL; 1 g/L is the conventional clinical-PK reporting unit). Document per-model via `covariateData[[AAG]]$units` if a different unit is used.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used in piecewise-linear, power, or exponential effect forms (e.g., the Kloft 2006 cytotoxic-chemotherapy myelosuppression family fits a piecewise-linear effect with breakpoint at the cohort median 1.34 g/L: separate slopes apply for `AAG <= 1.34` and `AAG > 1.34`). Reference values observed: 1.34 g/L (Kloft 2006 / Netterberg 2017, cohort median in mixed adult-cancer cohort).
- **Source aliases:**
  - `AAG` -- used in `Netterberg_2017_docetaxel.R` (per the bundle's NM-TRAN $INPUT block; matching Kloft 2006).
  - `AGP1` -- used in `Ozawa_2007_docetaxel.R` (Appendix I $INPUT block); reported in mg/dL with conversion to canonical g/L via `AAG_g_per_L = AGP_mg_per_dL / 100`.
- **Example models:** `Netterberg_2017_docetaxel.R` (piecewise-linear effects on baseline ANC with separate low-AAG and high-AAG slopes around median 1.34 g/L; linear effect on the drug-effect slope SL via `(1 + theta * (AAG - 1.34))`), `Ozawa_2007_docetaxel.R` (multiplicative power-form effect on the linear drug-effect slope: `SLOPE = theta_SLOPE * (AAG / 0.94)^e_aag_slope` with `e_aag_slope = -1.38`; reference value 0.94 g/L from the published NONMEM control stream `AGPm = 94` mg/dL).
- **Notes:** General scope because serum AAG is a routine clinical-laboratory measurement that recurs across cytotoxic-chemotherapy population-PK / PD analyses (Bruno 1996/1998 docetaxel popPK uses AAG as a CL covariate; Kloft 2006 and downstream Friberg-family myelosuppression models use it on baseline ANC and drug-effect slope; Ozawa 2007 uses a power-form effect on the drug-effect slope only). Time-fixed at baseline unless the source paper states otherwise. The Kloft 2006 piecewise-linear breakpoint at 1.34 g/L corresponds to the population median in their pooled cancer cohort; the Ozawa 2007 normalisation reference 0.94 g/L is the cohort median used inside that paper's NONMEM control stream. Future papers may use different breakpoints / reference values, so document the per-model reference in `covariateData[[AAG]]$notes`. Distinct from `CRP` (a different acute-phase reactant with different binding properties).

### IL6 (**canonical for serum interleukin-6 concentration**)
- **Description:** Serum (or plasma) interleukin-6 (IL-6) concentration. Pro-inflammatory cytokine; elevated in rheumatoid arthritis, Castleman's disease, sepsis, COVID-19, and other inflammatory conditions. Both baseline and time-varying usages are covered; document per-model in `covariateData[[IL6]]$notes` whether the column is baseline-only or time-varying.
- **Units:** pg/mL (= ng/L; the two labels are numerically equivalent).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(IL6 / ref)^exponent` or with the log-transformed form `(log(IL6 * 1000) / log(ref * 1000))^exponent` that some legacy NONMEM analyses adopt. Reference values observed: 20 pg/mL (Frey 2013 baseline; the formula `(log(IL-6 * 1000)/9.9)^exponent` is the algebraic equivalent of `(log(IL-6) / log(20))^exponent` after the constant-factor rescaling that ties the reference log to 9.9 = log(20000)).
- **Source aliases:**
  - `BLIL6`, `bIL6`, `IL6_BASE` -- baseline IL-6 (used in some NONMEM control streams; canonical drops the `BL` prefix per the `EOS` / `EASI` / `AGE` convention with baseline-vs-time-varying status documented in per-model notes).
  - `IL-6`, `IL_6` -- punctuation variants seen in publication tables and figures.
- **Example models:** `Frey_2013_tocilizumab.R` (baseline IL-6, reference 20 pg/mL; log-transformed power effects on EC50, BASE, and the DMARD background-effect parameter).
- **Notes:** Anti-IL-6 / anti-IL-6R drugs (tocilizumab, sarilumab, siltuximab) often see large transient increases in measured IL-6 after dosing -- when the column is time-varying for these drugs, the assay typically detects total (free + drug-bound) IL-6 because the antibody complex slows clearance of IL-6. Document the assay type (free vs. total IL-6) and whether the value is pre-dose (baseline-only) or post-dose / time-varying in each model's `covariateData[[IL6]]$notes`. Frey 2013's formula relies on the relative log-IL-6 ratio rather than the linear concentration, so the column units must be pg/mL exactly (not ng/mL) for the published exponents to apply unchanged.

## Cardiometabolic / target biomarkers

### HDLC (**canonical for high-density lipoprotein cholesterol**)
- **Description:** Serum high-density lipoprotein cholesterol concentration (baseline or time-varying).
- **Units:** mg/dL or mmol/L -- document the unit used in each model via `covariateData[[HDLC]]$units` (1 mmol/L ~= 38.67 mg/dL for cholesterol).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(HDLC / ref)^exponent`. Reference value observed: 54 mg/dL (Frey 2010 pooled-cohort median).
- **Source aliases:**
  - `HDL-C` -- Frey 2010 spelling with hyphen.
  - `HDL_C` -- common alternative spelling.
- **Example models:** `Frey_2010_tocilizumab.R` (mg/dL, reference 54; small negative exponent -0.2 on linear CL; the paper interprets the effect as a body-size surrogate rather than a mechanism).
- **Notes:** Cardiometabolic lipid-panel covariate. In Frey 2010 it correlates with body size (HDL-C is lower in larger patients) and the small CL effect (-14% to +15% across the observed 23-135 mg/dL range) was retained but not interpreted as mechanistic.

### TCHOL (**canonical for total serum cholesterol**)
- **Description:** Serum (or plasma) total cholesterol concentration (baseline or time-varying). Sum of HDL, LDL, VLDL, and other lipoprotein-associated cholesterol fractions; distinct from `HDLC` which captures only the HDL fraction.
- **Units:** mmol/L or mg/dL -- document the unit used in each model via `covariateData[[TCHOL]]$units` (1 mmol/L ~= 38.67 mg/dL for cholesterol).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with linear-deviation forms `1 + theta * (TCHOL - ref)` or power scaling `(TCHOL / ref)^exponent`. Reference value observed: 3 mmol/L (Archary 2018, severely malnourished pediatric LPV cohort, baseline mean 2.7-2.9 mmol/L).
- **Source aliases:**
  - `CHOL` -- Archary 2018 NONMEM column abbreviation; the universal short form.
  - `TC` -- alternative abbreviation common in lipid-panel literature.
- **Example models:** `Archary_2018_lopinavir.R` (mmol/L, reference 3; linear effect on apparent CL/F: `1 + 0.207 * (TCHOL - 3)`; serves as a surrogate for nutritional / hepatic-function recovery rather than a mechanistic effect).
- **Notes:** In severe acute malnutrition, total cholesterol tracks lipid-pool repletion and hepatic recovery and may serve as a surrogate covariate when the actual driver of bioavailability or clearance variability is not directly measured. Lopinavir is highly protein-bound (to albumin and AAG) and lipophilic, so circulating cholesterol can correlate with binding capacity for lipophilic drugs. Distinct from `HDLC` (HDL fraction only) and from `CRP` / `ALB` / `AAG` (separate canonicals for related malnutrition / inflammation markers).

### TRIG (**canonical for serum triglyceride concentration**)
- **Description:** Serum triglyceride concentration (baseline or time-varying).
- **Units:** mmol/L or mg/dL -- document the unit used in each model via `covariateData[[TRIG]]$units` (1 mmol/L ~= 88.5 mg/dL for triglyceride).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with linear-deviation form `(1 + theta * (TRIG - ref))` or power scaling `(TRIG / ref)^exponent`. Reference value observed: 5.3 mmol/L (Archary 2019 lamivudine equation centring -- the equation centring in the source is reported as the cohort average; cohort median in Archary 2019 Table 1 is 2.2-2.3 mmol/L, see model-file Errata).
- **Source aliases:** none known.
- **Example models:** `Archary_2019_lamivudine.R` (mmol/L, reference 5.3; linear-deviation effect on Vc/F with coefficient -0.13 per mmol/L deviation from the reference; lower triglyceride implies higher apparent central volume).
- **Notes:** Cardiometabolic lipid-panel covariate. In Archary 2019 the inverse triglyceride--Vc relationship is interpreted as a nutritional-status / hydrophilic-drug-distribution surrogate: triglycerides rise during nutritional rehabilitation, and lamivudine (a hydrophilic drug) shows decreasing apparent volume as nutritional status improves. The linear-deviation form preserves the source's published parameterization; per-model `covariateData[[TRIG]]$units` is load-bearing because the centring reference and slope are unit-specific.

### LDLC (**canonical for low-density lipoprotein cholesterol**)
- **Description:** Serum low-density lipoprotein cholesterol concentration (baseline or time-varying). The pharmacologically meaningful endpoint for lipid-lowering therapies (statins, PCSK9 inhibitors, ANGPTL3 inhibitors); also serves as a baseline covariate or as the time-varying PD response in indirect-response exposure-response models.
- **Units:** mg/dL or mmol/L -- document the unit used in each model via `covariateData[[LDLC]]$units` (1 mmol/L ~= 38.67 mg/dL for cholesterol).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(LDLC / ref)^exponent` for the baseline-LDLC covariate role, or with no reference (used directly as the PD response state) when it is the modelled output. Reference values observed: 211 mg/dL (Pu 2021 HoFH typical-patient definition).
- **Source aliases:**
  - `LDL-C` -- common spelling with hyphen.
  - `LDL_C` -- common alternative spelling.
  - `LDLBL` (baseline LDL-C) -- used in `Pu_2021_evinacumab.R` (Pu 2021 NM-TRAN $INPUT column for centred baseline LDL-C as a covariate on IC50).
- **Example models:** `Pu_2021_evinacumab.R` (mg/dL, baseline reference 211 mg/dL; power exponent -1.17 on IC50, where higher baseline LDL-C predicts a smaller IC50 and therefore greater sensitivity to evinacumab; LDL-C is also the PD output state initialised at the baseline value).
- **Notes:** Cardiometabolic lipid-panel covariate. Distinct from `HDLC` (high-density lipoprotein cholesterol) and from any total-cholesterol or non-HDL-C derivation. When LDL-C is both the response variable AND a covariate (as in Pu 2021, where baseline LDLC drives IC50 and the time-varying state is the modelled PD), document the dual role in `covariateData[[LDLC]]$notes`.

### ANGPTL3 (**canonical for angiopoietin-like protein 3 concentration**)
- **Description:** Total serum angiopoietin-like protein 3 (ANGPTL3) concentration. ANGPTL3 is the pharmacological target for anti-ANGPTL3 monoclonal antibodies (evinacumab) and antisense oligonucleotides (vupanorsen). Baseline ANGPTL3 acts as a soluble-target biomarker that contributes to target-mediated drug disposition; higher baseline target predicts a higher saturable Vmax.
- **Units:** mg/L (equivalent to ug/mL). Document per-model via `covariateData[[ANGPTL3]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with power scaling `(ANGPTL3 / ref)^exponent`. Reference value observed: 0.08 mg/L (Pu 2021 typical-patient median).
- **Source aliases:**
  - `ANGBL` (baseline ANGPTL3) -- used in `Pu_2021_evinacumab.R` (Pu 2021 NM-TRAN $INPUT column for centred baseline ANGPTL3).
  - `ANGBL0` -- alternative Pu 2021 raw column.
- **Example models:** `Pu_2021_evinacumab.R` (mg/L, baseline reference 0.08 mg/L; power exponent +0.405 on Vmax, where higher baseline target predicts a faster saturable elimination -- biologically consistent with evinacumab being co-cleared along with bound ANGPTL3).
- **Notes:** Specific scope because the column is meaningful only for drugs whose mechanism involves ANGPTL3 (anti-ANGPTL3 mAbs / ASOs). Reusing the name for another anti-ANGPTL3 agent is acceptable (extend the example-models list). The assay in Pu 2021 detects both free and target-bound ANGPTL3 after acid pretreatment of serum; document the assay type (free vs. total) per model in `covariateData[[ANGPTL3]]$notes`.

### FPCSK9 (**canonical for free (unbound) proprotein convertase subtilisin/kexin type 9 concentration**)
- **Description:** Free (unbound, non-drug-bound) serum proprotein convertase subtilisin/kexin type 9 (PCSK9) concentration. For anti-PCSK9 monoclonal antibodies (alirocumab, evolocumab, bococizumab) the free-PCSK9 pool is the pharmacologically active target fraction; drug-target binding reduces FPCSK9 relative to total PCSK9.
- **Units:** ng/mL (document per-model if a paper reports a different unit).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with linear-deviation forms `TVPARAM + theta * (FPCSK9 / ref)` or power-form `(FPCSK9 / ref)^theta`. Reference values observed: 72.9 ng/mL (Martinez 2019 time-varying median).
- **Source aliases:** none known.
- **Example models:** `Martinez_2019_alirocumab.R` (time-varying; additive-linear effect on `Km` with slope -0.541 per (FPCSK9/72.9), reference 72.9 ng/mL).
- **Notes:** Specific scope because the column is meaningful only for drugs whose mechanism is PCSK9 inhibition; reusing the name for a different anti-PCSK9 agent is acceptable (add to Example models). For non-PCSK9 drugs that use a similar target-concentration biomarker, register a new canonical (e.g., `FIL6R`, `FTNF`) rather than overloading `FPCSK9`. Per-model `covariateData[[FPCSK9]]$notes` should state whether the value is baseline-only or time-varying and how missing values were imputed (Martinez 2019 used LOCF).

### SBCMA (**canonical for soluble B-cell maturation antigen concentration**)
- **Description:** Baseline serum (or plasma) concentration of soluble B-cell maturation antigen (sBCMA), the shed extracellular domain of the BCMA receptor (TNFRSF17). Serves as a soluble-target biomarker for BCMA-directed therapeutics (anti-BCMA antibody-drug conjugates such as belantamab mafodotin, BCMA-targeted bispecifics, and BCMA CAR-T) -- sBCMA is elevated in multiple-myeloma and reflects tumour burden, and contributes to target-mediated drug disposition by sequestering circulating drug.
- **Units:** ng/mL (equivalent to ug/L; 1 ng/mL = 1 ug/L). Document per-model via `covariateData[[SBCMA]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with power scaling `(SBCMA / ref)^exponent`. Reference value observed: 50 ng/mL (Papathanasiou 2025 typical-patient definition).
- **Source aliases:**
  - `SBCMABL` (baseline soluble BCMA) -- used in `Papathanasiou_2025_belantamab.R`.
- **Example models:** `Papathanasiou_2025_belantamab.R` (ng/mL, reference 50; power exponents on initial CL +0.113, on ADC Vc +0.0401, on Imax +0.160).
- **Notes:** Specific scope because the column is meaningful only for drugs whose mechanism involves the BCMA receptor (and thus a circulating soluble-target pool). Reusing the name for another anti-BCMA agent is acceptable (extend the example-models list). For other oncology TMDD targets register a new canonical (e.g., `HER2_ECD` already exists for HER2; an analogous `SCD20`, `SCD38` would follow the same pattern). Multiple myeloma populations show sBCMA spanning roughly 2 to 2,000 ng/mL, so the (SBCMA/50)^exponent form should be evaluated with care over the full clinical range.
## Drug exposure metrics

### CAV (**canonical for average drug plasma concentration over a dosing interval**)
- **Description:** Average plasma concentration of the modelled drug over a dosing interval (Cav = AUC_tau / tau). Used as the time-varying or per-period exposure metric in exposure-response models that feed individual empirical-Bayes PK predictions from a previously published population PK model into a downstream PD model.
- **Units:** ug/mL (document per-model via `covariateData[[CAV]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used in Emax/EC50 (e.g., `Emax * CAV / (EC50 + CAV)`) or power (e.g., `(CAV / CavMedian)^exponent`) drug-effect terms. Set to 0 for placebo periods.
- **Source aliases:** `CAV`, `Cav`, `CAVG`.
- **Example models:** `FiedlerKelly_2020_fremanezumab_em.R`, `FiedlerKelly_2020_fremanezumab_cm.R`, `Schoemaker_2018_levetiracetam.R` (DDMODEL00000239; LEV plasma concentration in mg/L).
- **Notes:** Specific scope because the value is intrinsically tied to the modelled drug -- there is no shared meaning across drugs or studies. Each model's `covariateData[[CAV]]$notes` should state how the Cav values are derived (e.g., empirical-Bayes from a referenced population PK model) and that the column is set to 0 for placebo periods.

### DOSE_PHT_MGKGD (**canonical for daily phenytoin dose per kg body weight**)
- **Description:** Patient's own total daily dose of phenytoin (mg) divided by current body weight (kg), expressed as mg/kg/d. Per-dose-record covariate; constant within an inter-dose interval and updated when the prescriber alters the daily dose.
- **Units:** mg/kg/d
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- enters as a self-dose-rate regressor in the dose-dependent powder bioavailability formula `F_powder = 1 - exp(-theta / DOSE_PHT_MGKGD)`. Has no effect when paired with `FORM_POWDER = 0` (tablet); a non-NA non-zero placeholder must still be supplied.
- **Source aliases:**
  - `Dij` -- used in `Yukawa_1990_phenytoin.R` (paper's per-record daily-dose-per-weight regressor, mg/kg/d, in the powder bioavailability equation 4 of Yukawa 1990).
- **Example models:** `Yukawa_1990_phenytoin.R` (powder bioavailability `F_powder = 1 - exp(-9.92 / DOSE_PHT_MGKGD)`; F approaches 1 below ~2 mg/kg/d and decreases monotonically as the daily dose increases, reflecting the lower wettability of the Aleviatin brand phenytoin powder formulation).
- **Notes:** Specific scope because the value is intrinsically tied to phenytoin (PHT) and the Yukawa 1990 powder-vs-tablet bioavailability contrast. Drug-self-dose covariates for other drugs should register sibling canonicals (e.g., `DOSE_<DRUG>_MGKGD`) rather than reuse this name -- the absolute coefficient (theta_BA2 = 9.92 in Yukawa 1990) is not transferable across drugs. Computed as the total daily dose summed across the 2-3 daily phenytoin doses (mg/d) divided by the patient's body weight at the dose record (kg). Ratified canonically on 2026-05-10 alongside the Yukawa 1990 phenytoin extraction.

### PRED_DOSE (**canonical for concomitant oral prednisolone daily dose**)
- **Description:** Concomitant oral prednisolone (or prednisolone-equivalent glucocorticoid) daily dose. Time-varying across the dosing period as the post-transplant conmed_steroid taper progresses.
- **Units:** mg/day. Document per-model via `covariateData[[PRED_DOSE]]$units` when a paper uses a different unit (mg/kg/day) or a different glucocorticoid (methylprednisolone, dexamethasone, hydrocortisone) -- in the latter case convert to prednisolone-equivalent mg/day before populating the column and record the conversion factor in `covariateData[[PRED_DOSE]]$notes`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- continuous, with 0 mg/day (no prednisolone) the natural reference value. In Storset 2014 the effect on tacrolimus oral bioavailability is modeled as a sigmoid-Emax fractional reduction `(1 - Pred_max * PRED_DOSE / (Pred_50 + PRED_DOSE))` with `Pred_max = 0.67` and `Pred_50 = 35 mg/day` (Hill = 1); document the per-model functional form in `covariateData[[PRED_DOSE]]$notes`.
- **Source aliases:**
  - `Prednisolone dose` -- used in `Storset_2014_tacrolimus.R` (mg/day).
- **Example models:** `Storset_2014_tacrolimus.R` (Emax-style fractional reduction in tacrolimus oral bioavailability via prednisolone-driven induction of intestinal CYP3A / P-glycoprotein; Storset 2014 Methods Equations 4 + 6 with Hill = 1).
- **Notes:** Distinct from `PRED_CMAX_FREE` (free prednisolone Cmax co-medication exposure) -- `PRED_CMAX_FREE` is the modelled-from-data peak free concentration, whereas `PRED_DOSE` is the administered daily-dose level supplied directly from the dosing record. Both can coexist in a future model that simultaneously tests dose-driven and exposure-driven effects of prednisolone on a tacrolimus PK parameter. Distinct from `CONMED_STEROID` (binary baseline / concomitant corticosteroid use indicator) and `PRICORT` (binary prior corticosteroid use indicator) -- `PRED_DOSE` carries the daily dose value, not just an on / off flag. Time-varying because tacrolimus PK depends on the conmed_steroid dose at the time of each tacrolimus observation; the conmed_steroid taper schedule must be supplied as a per-time-row covariate column. The corresponding methylprednisolone single-dose induction-bolus indicator (Storset 2014 binary covariate, not retained in the final model) would warrant a separate canonical (e.g. `MPRED_BOLUS`) if a future model retains it. Ratified canonically on 2026-05-08.

### PRED_CMAX_FREE (**canonical for free prednisolone Cmax co-medication exposure**)
- **Description:** Maximum free (unbound, ultrafiltrable) plasma prednisolone concentration over a co-medication dosing interval, used as a co-medication-exposure covariate on a primary modelled drug (e.g., tacrolimus) when concomitant prednisolone is suspected to alter the primary drug's PK via membrane-permeability or fluid-balance effects. Time-fixed per subject in the source paper (one Cmax value per subject, derived from limited-sampling concentrations).
- **Units:** nmol/L
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used in centred-deviation form `(1 + e * (PRED_CMAX_FREE - ref))`. Reference value observed: 155.5 nmol/L (Bergmann 2014 study median).
- **Source aliases:**
  - `PredCmax,free` / `PREDCFR` -- used in `Bergmann_2014_tacrolimus.R` (Bergmann 2014 Table 2 footnote; 162 nmol/L median per Table 1, 155.5 nmol/L centring value per Table 2 equation).
- **Example models:** `Bergmann_2014_tacrolimus.R` (linear deviation effect on tacrolimus apparent central volume V1/F: every 1 nmol/L increase from 155.5 nmol/L decreases V1/F by 0.28%).
- **Notes:** Specific scope because the value is intrinsically tied to a particular co-medication (oral prednisolone) and a particular reference population (kidney transplant recipients on a tapering conmed_steroid regimen). Future popPK models that test free prednisolone Cmax as a covariate on a different primary drug should reuse this canonical; if a future model uses a different prednisolone exposure metric (e.g., `PRED_AUC0_12_FREE` for free AUC 0-12 h, or `PRED_CMAX_TOTAL` for total Cmax), register parallel canonicals rather than overload `PRED_CMAX_FREE`. The source-paper Cmax is derived from limited-sampling concentrations at 1 / 2 / 4 hours postdose per Bergmann 2014 Methods (validated against earlier full-profile data from the same cohort). Distinct from `CAV` (average concentration of the modelled drug) and `CP_MGL` (instantaneous concentration of the modelled drug as a time-varying PD driver) -- `PRED_CMAX_FREE` is the maximum concentration of a co-medication, used as a per-subject covariate. Ratified canonically on 2026-05-08 alongside the Bergmann 2014 extraction.

### AUC_CARBO (**canonical for per-cycle average AUC of carboplatin**)
- **Description:** Per-cycle average AUC of carboplatin used as the time-varying drug-exposure covariate driving cytotoxic tumour-death rates in tumour-size-dynamics models of platinum-based chemotherapy. The value is held step-wise constant over each chemotherapy cycle and resets at the start of the next cycle.
- **Units:** carboplatin AUC units (typically `mg*min/mL`); document per-model via `covariateData[[AUC_CARBO]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- set to 0 in cycles where carboplatin is not administered (e.g., post-discontinuation or non-platinum arms).
- **Source aliases:**
  - `CB` (NONMEM `$INPUT` column in DDMODEL00000217 / DDMODEL00000218) -- used in `Zecchin_2016_tumorovarian.R` and `Zecchin_2016_survival.R`. The DDMORE bundles ship the simulated datasets with the column re-labelled `AUC0`; downstream consumers should map `AUC0` -> `AUC_CARBO`.
- **Example models:** `Zecchin_2016_tumorovarian.R` (Zecchin 2016 SLD model for advanced ovarian cancer, DDMODEL00000217), `Zecchin_2016_survival.R` (Zecchin 2016 OS model, DDMODEL00000218; the OS model integrates the same SLD ODE inline, with the prior IPP-fit subject-level KG/KD0/KD1/IBASE supplied via the dataset).
- **Notes:** Specific scope because the column meaning is tied to a particular cytotoxic agent (carboplatin) and a particular per-cycle averaging convention. Reusing the same column for another platinum analogue (cisplatin, oxaliplatin) is not appropriate -- register a sibling canonical (`AUC_CISPLATIN`, `AUC_OXALIPLATIN`) when needed. The Zecchin 2016 SLD and OS models use the value directly in the death-rate term `kd0 * AUC_CARBO * tumorSize`, with an internal `/1000` numerical scaling carried verbatim from the source `$DES` block.

### AUC_BAST_FW (**canonical for first-week AUC in the BAST PTTE 2017 teaching dataset**)
- **Description:** Area under the plasma concentration-time curve over the first week of treatment for the unspecified hypothetical drug used in the BAST PTTE 2017 teaching guiding-document (DDMODEL00000243). Per-subject, time-fixed (a single early-exposure summary value carried as a baseline covariate). Used to drive the Event 2 Gompertz hazard via centred-deviation form `exp(coef * (AUC_BAST_FW - 3065.5))`.
- **Units:** ug*h/L (per BAST PTTE guiding document Section  2.2.2)
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used in centred-deviation form. Reference value observed: 3065.5 ug*h/L (BAST PTTE 2017 simulated cohort median; runEV2_105 base hazard).
- **Source aliases:**
  - `AUC` -- verbatim NM-TRAN `$INPUT` column name in DDMODEL00000243's executable .mod files and in the bundle's `Simulated_event_data.csv`. Renamed to `AUC_BAST_FW` in the canonical register so that future models using a generic `AUC` column with different semantics will not silently collide.
- **Example models:** `NA_NA_tte_gompertz_ev2.R` (BAST PTTE 2017 / DDMODEL00000243 Event 2 hazard model; centred at 3065.5 ug*h/L, coefficient 3.09e-4 on the NONMEM rescaled scale, equivalent to `exp(3.09e-4 * (AUC_BAST_FW - 3065.5))` on the hazard).
- **Notes:** Specific scope because the value is intrinsically tied to an unspecified hypothetical drug and a fictional simulated population in the BAST PTTE 2017 guiding document -- the column has no shared meaning across drugs or studies. Future models should register a sibling canonical (e.g., `AUC_<DRUG>`) with explicit drug semantics rather than overload this name. The BAST guiding document Section  2.2.2 defines this as "AUC of drug treatment given within the first week (ug*h/L)."

### AUC_GEM (**canonical for per-cycle average AUC of gemcitabine**)
- **Description:** Per-cycle average AUC of gemcitabine (sum of parent and active metabolite exposure, per Zecchin 2016 Methods) used as the time-varying drug-exposure covariate driving cytotoxic tumour-death rates in tumour-size-dynamics models of gemcitabine-containing chemotherapy.
- **Units:** gemcitabine AUC units (typically `mg*h/L` or the paper's `mol*day / 10^6 cells` scaling for the parent-plus-active-metabolite composite); document per-model via `covariateData[[AUC_GEM]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- set to 0 in cycles where gemcitabine is not administered (e.g., carboplatin-monotherapy arms).
- **Source aliases:**
  - `G` (NONMEM `$INPUT` column in DDMODEL00000217 / DDMODEL00000218) -- used in `Zecchin_2016_tumorovarian.R` and `Zecchin_2016_survival.R`. The DDMORE bundles ship the simulated datasets with the column re-labelled `AUC1`; downstream consumers should map `AUC1` -> `AUC_GEM`.
- **Example models:** `Zecchin_2016_tumorovarian.R` (Zecchin 2016 SLD model for advanced ovarian cancer, DDMODEL00000217), `Zecchin_2016_survival.R` (Zecchin 2016 OS model, DDMODEL00000218).
- **Notes:** Specific scope. The Zecchin 2016 SLD and OS models use the value directly in the death-rate term `kd1 * AUC_GEM * tumorSize`, with an internal `/100` numerical scaling carried verbatim from the source `$DES` block.

### AUC_GCV (**canonical for per-q12h-interval AUC of ganciclovir**)
- **Description:** Time-varying ganciclovir AUC over a q12h dosing interval (AUC_0-12), used as the drug-exposure input to indirect-response viral-turnover PK/PD models of cytomegalovirus (CMV) viral load decline under (val)ganciclovir treatment. The Koloskoff 2025 source computes individual AUC_0-12 from an upstream popPK model (Franck 2021) and feeds it to the PD model as a Monolix "varying input"; the PD model itself does not integrate a PK ODE, so AUC_GCV is supplied to nlmixr2 as a time-varying data column.
- **Units:** `mg*h/L` (document per-model via `covariateData[[AUC_GCV]]$units` if a different exposure unit is reported).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- set to 0 in pre-treatment / off-treatment records so the drug-stimulation term `Emax * AUC_GCV / (EC50 + AUC_GCV)` vanishes and the viral load returns to the `kin / kout` steady-state baseline.
- **Source aliases:**
  - `AUC_0-12` -- the printed variable name in Koloskoff 2025 (Methods Section 2.3, Eq. 1, and Table 1). Q24h dosing intervals are entered as `AUC_0-24 / 2` so all data live in the q12h framework (Koloskoff 2025 Methods Section 2.1).
- **Example models:** `Koloskoff_2025_ganciclovir.R` (Koloskoff 2025 indirect viral turnover model for CMV viral load in pediatric SOT / HSCT recipients; AUC_GCV enters the ODE via `kout * (1 + Emax * AUC_GCV / (EC50 + AUC_GCV)) * viralLoad`).
- **Notes:** Specific scope -- the column meaning is tied to ganciclovir as the drug and to a q12h interval-averaging convention. Sibling drug-specific AUC canonicals (`AUC_CARBO`, `AUC_GEM`, `AUC_BAST_FW`, `AUC_PAZO`) follow the same `AUC_<DRUG>` naming pattern; a future PK/PD model that uses a different exposure metric for ganciclovir (e.g., trough concentration, instantaneous concentration) should register a parallel canonical rather than overload `AUC_GCV`. Koloskoff 2025 Monte Carlo simulations are reported under AUC_0-24 (Tables 3 and 4) assuming AUC_0-24 = 2 x AUC_0-12 at steady state; nlmixr2 simulations should set AUC_GCV to the q12h-interval value (i.e., AUC_0-24 / 2).

### AUC_PAZO (**canonical for per-period mean AUC of pazopanib**)
- **Description:** Per-period (per-dose-group in preclinical xenograft studies; per-subject mean dose-adjusted in clinical studies) mean AUC of pazopanib used as the drug-exposure covariate driving the antiangiogenic and cytotoxic effect rates in semi-mechanistic tumour-growth / angiogenesis-inhibition (TGI) models of pazopanib in renal-cell carcinoma. Time-varying step-wise (held constant within a treatment period and resetting when dose level changes or treatment ends).
- **Units:** `ug*h/mL` (`= mg*h/L`). Document per-model via `covariateData[[AUC_PAZO]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- enters via power form `a = a0 * AUC_PAZO^e_auc_pazo_a` and `c = c0 * AUC_PAZO^e_auc_pazo_c`. Set to 0 in periods where pazopanib is not administered; the model rate terms reduce to baseline drug-off (a, c -> 0 when AUC_PAZO -> 0) under that convention.
- **Source aliases:**
  - `AUC` -- used in `Ouerdani_2015_pazopanib_mouse.R` and `Ouerdani_2015_pazopanib.R` (Ouerdani 2015 Methods Equations 2-3; preclinical values 220.2, 656.8, 1140.8 ug*h/mL for the 10, 30, 100 mg/kg mouse dose groups; clinical mean 771.6 ug*h/mL for 800 mg QD pazopanib in RCC patients, with per-subject values 629.4-802.4 ug*h/mL derived from an Emax fit to mean AUCs at the patient's dose history).
- **Example models:** `Ouerdani_2015_pazopanib_mouse.R` (preclinical TGI in CAKI-2 xenograft mice; AUC enters as `c = c0 * AUC_PAZO^0.332` only -- the cytotoxic exponent `e_auc_pazo_a` was fixed to 0 because the in-mouse cytotoxic effect did not vary with exposure across the 10-100 mg/kg range), `Ouerdani_2015_pazopanib.R` (clinical TGI in RCC patients; AUC enters as both `a = a0 * AUC_PAZO^0.125` and `c = c0 * AUC_PAZO^0.142`).
- **Notes:** Specific scope because the column meaning is tied to a particular drug (pazopanib) and the power-form parameterisation that the Ouerdani 2015 model uses. Reusing the name for a different tyrosine-kinase inhibitor is not appropriate -- register a sibling canonical (`AUC_SORAF` for sorafenib, `AUC_SUNI` for sunitinib, etc.) when needed. The Ouerdani 2015 paper reports the preclinical AUCs in `ug*h/mL` from a separate preclinical PK study (cited as the FDA Pharmacology Review for pazopanib NDA 022465); the clinical AUCs come from an Emax fit (Equation derived from Methods) to pooled mean AUCs at varying daily doses (5 mg to 2000 mg) across five prior pazopanib trials. Ratified canonically on 2026-05-12 alongside the Ouerdani 2015 pazopanib mouse and clinical extractions.

### CLI (**canonical for individual posthoc clearance from an upstream popPK fit**)
- **Description:** Subject-specific empirical-Bayes (posthoc) total plasma clearance from a separately published population PK model that the current PD model treats as a fixed input. Used as a per-subject (time-fixed) covariate in PD-only models that derive a per-cycle exposure metric (e.g., AUC = DOSE / CLI) without instantiating a PK ODE.
- **Units:** L/h (document per-model via `covariateData[[CLI]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- appears directly in derived-exposure expressions; not a covariate effect coefficient.
- **Source aliases:**
  - `CL` -- used in the Hansson 2013 sunitinib biomarker / TGI / fatigue PD-model family (DDMODEL00000197 and siblings, including DDMODEL00000222) and in `Schindler_2016_sunitinib.R` (DDMODEL00000221) as the posthoc CL column from the paper's upstream 2-compartment popPK fit.
- **Example models:** `Hansson_2013a_sunitinib.R` (DDMODEL00000197; typical-value reference 32.819 L/h, drawn from the bundle's simulated dataset for subject 1 -- broadly consistent with Houk et al. 2010 typical sunitinib CL), `Hansson_2013b_sunitinib.R` (DDMODEL00000198; tumor growth inhibition variant; same per-subject `CL` column as Hansson 2013a/c), `Hansson_2013c_sunitinib.R` (DDMODEL00000222; uses a per-record `CL` column with subject-specific values 30-43 L/h in the bundle's three-subject simulated dataset), `Schindler_2016_sunitinib.R` (DDMODEL00000221; per-subject post-hoc CL fed in as the `CL` column, vignette uses 50 L/h literature-typical sunitinib CL/F per Houk 2010).
- **Notes:** Specific scope because the values are intrinsically tied to a specific upstream popPK fit (sunitinib in this case; another model adopting CLI would carry its own upstream-PK lineage). Renamed from the source's `CL` column because `cl` is the canonical nlmixr2 PK parameter name (a parameter, not a data column). Each model's `covariateData[[CLI]]$notes` should cite the upstream popPK source (paper or DDMORE ID) and explain how to populate the column for new simulations (typically: simulate the upstream popPK first to obtain individual CL, or set every subject to the typical-value CL for typical-trajectory simulations). Distinct from `DOSE` (current administered dose level) -- the two columns jointly carry a per-cycle drug-exposure summary (`AUC = DOSE / CLI`) for PD-only models that consume posthoc PK from an upstream popPK fit instead of instantiating their own PK ODE.

### CL_INDIV (**canonical for per-subject empirical-Bayes drug clearance**)
- **Description:** Individual point estimate of the modelled drug's clearance, supplied per-subject as a fixed data column. Used in sequential PK->PD models where the PK structure has been fixed from a previously-published population PK analysis and the per-subject empirical-Bayes (POSTHOC) clearances are passed through to drive the PD layer rather than being re-estimated alongside the PD parameters.
- **Units:** L/h (document per-model via `covariateData[[CL_INDIV]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used directly inside `model()` as `CL_INDIV` in place of an estimated `cl <- exp(lcl + etalcl)`.
- **Source aliases:** `CLI` -- used in `Friberg_2002_paclitaxel.R` (NM-TRAN data column for per-subject paclitaxel CL).
- **Example models:** `Friberg_2002_paclitaxel.R`.
- **Notes:** Specific scope because the value is intrinsically tied to the modelled drug -- there is no shared meaning across drugs. Each model's `covariateData[[CL_INDIV]]$notes` should state which upstream popPK source the EBE values come from (e.g., Henningsson 2001 paclitaxel popPK, fixed in the DDMORE encoding) and whether placebo periods are present. Companion volumes are registered as `VC_INDIV` and `VP_INDIV`.

### CMAX_M1 (**canonical for maximum drug plasma concentration during month 1**)
- **Description:** Empirical-Bayes maximum plasma concentration (Cmax) reached during the first month (or the first cycle) of dosing for the modelled drug, used as a per-subject early-exposure covariate in PD / safety models that condition later-time outcomes on early exposure. Continuous, time-invariant per subject (set once from the month-1 / cycle-1 PK simulation).
- **Units:** ng/mL (document per-model via `covariateData[[CMAX_M1]]$units` if a different concentration unit is reported).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used as an additive logit shift `(CMAX_M1 - CMAX_M1_REF) * theta` (Girard 2012) or as a power scaling depending on the source. Reference value observed: 0 ng/mL (Girard 2012 sets `MED17 = 0` as the centering reference).
- **Source aliases:**
  - `CMAXM1` -- used in `Girard_2012_pimasertib.R`.
- **Example models:** `Girard_2012_pimasertib.R` (additive logit shift on the cumulative-logit AE-score model: `theta_cmaxm1 * CMAXM1`).
- **Notes:** Specific scope because the value depends on the upstream drug-specific population-PK model used to derive the empirical-Bayes Cmax. Sibling of the existing `CAV` (average dosing-interval concentration); both are derived exposure metrics fed into downstream PD / safety models. Document the upstream PK model in `covariateData[[CMAX_M1]]$notes` for any future user.

### VC_INDIV (**canonical for per-subject empirical-Bayes central volume of distribution**)
- **Description:** Individual point estimate of the modelled drug's central volume of distribution, supplied per-subject as a fixed data column. Companion to `CL_INDIV` in sequential PK->PD encodings.
- **Units:** L (document per-model via `covariateData[[VC_INDIV]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used directly inside `model()` as `VC_INDIV` in place of an estimated `vc <- exp(lvc + etalvc)`.
- **Source aliases:** `V1I` -- used in `Friberg_2002_paclitaxel.R` (NM-TRAN data column for per-subject paclitaxel V1).
- **Example models:** `Friberg_2002_paclitaxel.R`.
- **Notes:** See `CL_INDIV` notes for the broader convention.

### VP_INDIV (**canonical for per-subject empirical-Bayes peripheral volume of distribution**)
- **Description:** Individual point estimate of the modelled drug's first peripheral volume of distribution, supplied per-subject as a fixed data column. Companion to `CL_INDIV` and `VC_INDIV` in sequential PK->PD encodings.
- **Units:** L (document per-model via `covariateData[[VP_INDIV]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used directly inside `model()` as `VP_INDIV` in place of an estimated `vp <- exp(lvp + etalvp)`.
- **Source aliases:** `V2I` -- used in `Friberg_2002_paclitaxel.R` (NM-TRAN data column for per-subject paclitaxel V2).
- **Example models:** `Friberg_2002_paclitaxel.R`.
- **Notes:** See `CL_INDIV` notes for the broader convention. For models requiring a second peripheral compartment, register `VP2_INDIV` (and add a follow-on entry to this register) when a second model legitimately needs it.

### BAS_SVEGFR3 (**canonical for individual posthoc baseline soluble VEGFR-3 concentration from an upstream PD fit**)
- **Description:** Subject-specific empirical-Bayes (posthoc) baseline plasma sVEGFR-3 (soluble vascular endothelial growth factor receptor 3) concentration from a separately published indirect-response biomarker model that the current downstream PD model treats as a fixed input. Used as a per-subject (time-fixed) covariate in fatigue / adverse-event PD models that consume the upstream biomarker dynamics as data covariates without instantiating the biomarker ODE.
- **Units:** pg/mL (document per-model via `covariateData[[BAS_SVEGFR3]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- appears directly inside `model()` as the initial condition `svegfr3(0) <- BAS_SVEGFR3` and inside the relative-change driver `bm = (svegfr3 - BAS_SVEGFR3) / BAS_SVEGFR3`.
- **Source aliases:**
  - `BAS3` -- used in `Hansson_2013b_sunitinib.R` (DDMODEL00000198) and `Hansson_2013c_sunitinib.R` (DDMODEL00000222) as the posthoc sVEGFR-3 baseline column from the paper's upstream Hansson 2013a biomarker indirect-response fit (DDMODEL00000197).
- **Example models:** `Hansson_2013b_sunitinib.R` (DDMODEL00000198; tumor growth inhibition with sVEGFR-3 driven shrinkage), `Hansson_2013c_sunitinib.R` (DDMODEL00000222; bundle's three-subject simulated dataset reports BAS_SVEGFR3 values 42554-57365 pg/mL).
- **Notes:** Specific scope because the value is intrinsically tied to a specific upstream biomarker model (sVEGFR-3 indirect response under sunitinib in this case). The downstream fatigue model only consumes individual posthoc baseline / MRT / EC50 of the upstream biomarker; it does not re-fit them. Each model's `covariateData[[BAS_SVEGFR3]]$notes` should cite the upstream biomarker-PD source (paper or DDMORE ID) and explain how to populate the column for new simulations (typically: simulate from the upstream biomarker model to obtain individual posthoc baselines, or set every subject to the typical-value baseline for typical-trajectory simulations).

### MRT_SVEGFR3 (**canonical for individual posthoc mean residence time of soluble VEGFR-3 from an upstream PD fit**)
- **Description:** Subject-specific empirical-Bayes (posthoc) mean residence time of plasma sVEGFR-3 from a separately published indirect-response biomarker model that the current downstream PD model treats as a fixed input. Used as a per-subject (time-fixed) covariate in fatigue / adverse-event PD models that consume the upstream biomarker dynamics as data covariates without instantiating the biomarker ODE; appears as `kout3 = 1 / MRT_SVEGFR3` inside `model()`.
- **Units:** h (hours) -- document per-model via `covariateData[[MRT_SVEGFR3]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- appears directly in derived rate-constant expressions; not a covariate effect coefficient.
- **Source aliases:**
  - `MRT3` -- used in `Hansson_2013b_sunitinib.R` (DDMODEL00000198) and `Hansson_2013c_sunitinib.R` (DDMODEL00000222) as the posthoc sVEGFR-3 MRT column from the paper's upstream Hansson 2013a biomarker indirect-response fit (DDMODEL00000197).
- **Example models:** `Hansson_2013b_sunitinib.R` (DDMODEL00000198), `Hansson_2013c_sunitinib.R` (DDMODEL00000222; bundle's three-subject simulated dataset reports MRT_SVEGFR3 values 313-408 h, broadly consistent with the Hansson 2013a typical sVEGFR-3 MRT of 401 h).
- **Notes:** Specific scope; same upstream-biomarker dependency rationale as `BAS_SVEGFR3`. The downstream fatigue model consumes the upstream MRT directly without re-fitting it.

### EC50_SVEGFR3 (**canonical for individual posthoc drug-effect EC50 on soluble VEGFR-3 from an upstream PD fit**)
- **Description:** Subject-specific empirical-Bayes (posthoc) half-maximum-effect concentration of the modelled drug on the production of sVEGFR-3 (or, depending on the upstream model parameterization, on the analogous Imax inhibition pathway) from a separately published indirect-response biomarker model that the current downstream PD model treats as a fixed input. Used as a per-subject (time-fixed) covariate in fatigue / adverse-event PD models that consume the upstream biomarker dynamics as data covariates without instantiating the biomarker ODE.
- **Units:** mg*h/L (when the per-cycle drug-exposure summary is `auc = DOSE / CLI` in mg*h/L, EC50_SVEGFR3 carries the same units; document per-model via `covariateData[[EC50_SVEGFR3]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- appears directly inside `model()` in the simple-Imax drug-effect term `eff3 = auc / (EC50_SVEGFR3 + auc)`.
- **Source aliases:**
  - `EC53` -- used in `Hansson_2013b_sunitinib.R` (DDMODEL00000198) and `Hansson_2013c_sunitinib.R` (DDMODEL00000222) as the posthoc sVEGFR-3 EC50 column from the paper's upstream Hansson 2013a biomarker indirect-response fit (DDMODEL00000197).
- **Example models:** `Hansson_2013b_sunitinib.R` (DDMODEL00000198), `Hansson_2013c_sunitinib.R` (DDMODEL00000222; bundle's three-subject simulated dataset reports EC50_SVEGFR3 values 1.0-2.8 mg*h/L, consistent with the Hansson 2013a typical sVEGFR-3 IC50 typical value of 1.0 mg*h/L).
- **Notes:** Specific scope; same upstream-biomarker dependency rationale as `BAS_SVEGFR3`. The downstream fatigue model consumes the upstream EC50 directly without re-fitting it.

### BAS_SKIT (**canonical for individual posthoc baseline soluble KIT concentration from an upstream PD fit**)
- **Description:** Subject-specific empirical-Bayes (posthoc) baseline plasma sKIT (soluble stem cell factor receptor) concentration from a separately published indirect-response biomarker model that the current downstream PD model treats as a fixed input. Used as a per-subject (time-fixed) covariate in tumor-growth-inhibition / fatigue / adverse-event PD models that consume the upstream sKIT dynamics as data covariates without instantiating the biomarker ODE for sKIT in isolation.
- **Units:** pg/mL (document per-model via `covariateData[[BAS_SKIT]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- appears directly inside `model()` as the initial condition for treated and placebo sKIT compartments and inside the relative-change driver `(skit_pla - skit_drug) / skit_pla` (or the analogous BAS_SKIT-denominated form when only one sKIT compartment is simulated).
- **Source aliases:**
  - `SBAS` -- used in `Hansson_2013b_sunitinib.R` (DDMODEL00000198) as the posthoc sKIT baseline column from the paper's upstream Hansson 2013a biomarker indirect-response fit (DDMODEL00000197).
- **Example models:** `Hansson_2013b_sunitinib.R` (DDMODEL00000198; the Hansson 2013 e84 paper Table 2 reports a typical sKIT baseline of 39200 pg/mL with ~50% CV, matching `Hansson_2013a_sunitinib`'s typical value).
- **Notes:** Specific scope because the value is intrinsically tied to a specific upstream biomarker model (sKIT indirect response under sunitinib in this case). The downstream tumor-growth-inhibition model only consumes individual posthoc baseline / MRT / EC50 / DP-slope of the upstream biomarker; it does not re-fit them. Each model's `covariateData[[BAS_SKIT]]$notes` should cite the upstream biomarker-PD source (paper or DDMORE ID) and explain how to populate the column for new simulations (typically: simulate from the upstream biomarker model to obtain individual posthoc baselines, or set every subject to the typical-value baseline for typical-trajectory simulations). Sister covariates: `MRT_SKIT`, `EC50_SKIT`, `SLOPE_SKIT` (companions for the same upstream biomarker fit).

### MRT_SKIT (**canonical for individual posthoc mean residence time of soluble KIT from an upstream PD fit**)
- **Description:** Subject-specific empirical-Bayes (posthoc) mean residence time of plasma sKIT from a separately published indirect-response biomarker model that the current downstream PD model treats as a fixed input. Used as a per-subject (time-fixed) covariate in tumor-growth-inhibition / fatigue / adverse-event PD models that consume the upstream sKIT dynamics as data covariates without instantiating the biomarker ODE for sKIT in isolation; appears as `kout_skit = 1 / MRT_SKIT` inside `model()`.
- **Units:** h (hours) -- document per-model via `covariateData[[MRT_SKIT]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- appears directly in derived rate-constant expressions; not a covariate effect coefficient.
- **Source aliases:**
  - `SMRT` -- used in `Hansson_2013b_sunitinib.R` (DDMODEL00000198) as the posthoc sKIT MRT column from the paper's upstream Hansson 2013a biomarker indirect-response fit (DDMODEL00000197).
- **Example models:** `Hansson_2013b_sunitinib.R` (DDMODEL00000198; the Hansson 2013 e84 paper Table 2 reports a typical sKIT MRT of 101 days = 2424 h, matching `Hansson_2013a_sunitinib`'s typical value of 2430 h).
- **Notes:** Specific scope; same upstream-biomarker dependency rationale as `BAS_SKIT`. The downstream tumor-growth-inhibition model consumes the upstream MRT directly without re-fitting it.

### EC50_SKIT (**canonical for individual posthoc drug-effect EC50 on soluble KIT from an upstream PD fit**)
- **Description:** Subject-specific empirical-Bayes (posthoc) half-maximum-effect concentration of the modelled drug on the production of sKIT (or, depending on the upstream model parameterization, on the analogous Imax inhibition pathway) from a separately published indirect-response biomarker model that the current downstream PD model treats as a fixed input. Used as a per-subject (time-fixed) covariate in tumor-growth-inhibition / fatigue / adverse-event PD models that consume the upstream sKIT dynamics as data covariates without instantiating the biomarker ODE for sKIT in isolation.
- **Units:** mg*h/L (when the per-cycle drug-exposure summary is `auc = DOSE / CLI` in mg*h/L, EC50_SKIT carries the same units; document per-model via `covariateData[[EC50_SKIT]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- appears directly inside `model()` in the simple-Imax drug-effect term `eff_skit = auc / (EC50_SKIT + auc)`.
- **Source aliases:**
  - `SEC5` -- used in `Hansson_2013b_sunitinib.R` (DDMODEL00000198) as the posthoc sKIT EC50 column from the paper's upstream Hansson 2013a biomarker indirect-response fit (DDMODEL00000197).
- **Example models:** `Hansson_2013b_sunitinib.R` (DDMODEL00000198; the Hansson 2013 e84 paper Table 2 reports a typical (common across all four biomarkers) IC50 of 1.0 mg*h/L, matching `Hansson_2013a_sunitinib`'s shared typical value).
- **Notes:** Specific scope; same upstream-biomarker dependency rationale as `BAS_SKIT`. The downstream tumor-growth-inhibition model consumes the upstream EC50 directly without re-fitting it.

### SLOPE_SKIT (**canonical for individual posthoc linear disease-progression slope on soluble KIT from an upstream PD fit**)
- **Description:** Subject-specific empirical-Bayes (posthoc) linear disease-progression slope on the placebo / untreated sKIT compartment from a separately published indirect-response biomarker model that the current downstream PD model treats as a fixed input. Used as a per-subject (time-fixed) covariate in tumor-growth-inhibition models that need to simulate the placebo-arm sKIT trajectory as a comparator for the drug-arm sKIT compartment; appears as `dps_skit = BAS_SKIT * (1 + SLOPE_SKIT * t)` and `kin_skit = dps_skit * kout_skit` inside `model()`.
- **Units:** 1/h (document per-model via `covariateData[[SLOPE_SKIT]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- appears directly in the placebo-arm Kin expression for sKIT; not a covariate effect coefficient on a structural rate.
- **Source aliases:**
  - `SLO` -- used in `Hansson_2013b_sunitinib.R` (DDMODEL00000198) as the posthoc sKIT linear disease-progression slope column from the paper's upstream Hansson 2013a biomarker indirect-response fit (DDMODEL00000197).
- **Example models:** `Hansson_2013b_sunitinib.R` (DDMODEL00000198; the Hansson 2013 e84 paper Table 2 reports a typical disease-progression slope of 0.0261/month shared between VEGF and sKIT, which equals approximately 3.5e-5/h, matching `Hansson_2013a_sunitinib`'s typical value).
- **Notes:** Specific scope; same upstream-biomarker dependency rationale as `BAS_SKIT`. Sign convention follows the upstream biomarker fit -- positive slope means the placebo / natural-history sKIT trajectory drifts upward over time (capturing disease progression).

### CSS_RBV (**canonical for individual posthoc ribavirin steady-state plasma concentration from an upstream PK fit**)
- **Description:** Subject-specific empirical-Bayes (posthoc) modal estimate of the ribavirin (RBV) steady-state plasma trough concentration from a separately fitted population PK model (e.g., the Laouenan 2015 upstream RBV popPK fit) that the current downstream PD model treats as a fixed input. Used together with `K_RBV` inside `model()` to reconstruct the individual ribavirin concentration time-course analytically as `riba(t) = CSS_RBV * (1 - exp(-K_RBV * t))`.
- **Units:** ng/mL (document per-model via `covariateData[[CSS_RBV]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- appears directly inside `model()` in the analytical RBV-concentration expression `riba = CSS_RBV * (1 - exp(-K_RBV * t))` that drives the inhibition term `riba / (riba + ec50)`.
- **Source aliases:**
  - `css_mode` (modal posterior estimate of individual ribavirin Css from the upstream popPK fit, in ng/mL) -- used in `Laouenan_2015_ribavirin.R` (DDMODEL00000285). Rename `css_mode` -> `CSS_RBV` before passing the dataset to `rxSolve`.
- **Example models:** `Laouenan_2015_ribavirin.R` (DDMODEL00000285; the bundle's `Simulated_Laouenant_2015_CPTPSP_hb_RBV.txt` carries `css_mode` values 2,400-4,000 ng/mL in the 15-subject ANRS-CO20-CUPIC cohort).
- **Notes:** Specific scope because the value is intrinsically tied to a specific upstream popPK fit (ribavirin in HCV-cirrhotic patients on triple therapy in this case) and to the analytical `Css*(1-exp(-k*t))` parameterization; another drug or another popPK parameterization would carry its own canonical. Each model's `covariateData[[CSS_RBV]]$notes` should cite the upstream popPK source (paper or DDMORE ID) and explain how to populate the column for new simulations (typically: simulate the upstream popPK first to obtain individual Css and approach-rate, or set every subject to the typical Css for typical-trajectory simulations). Companion column: `K_RBV`. Distinct from `CAV` (dosing-interval-averaged concentration used in Emax / EC50 PD models with a single per-period exposure number) -- `CSS_RBV` carries the asymptotic steady-state value paired with the approach-rate constant, supporting the full time-course reconstruction.

### K_RBV (**canonical for individual posthoc ribavirin approach-to-steady-state rate constant from an upstream PK fit**)
- **Description:** Subject-specific empirical-Bayes (posthoc) modal estimate of the first-order rate constant governing the exponential approach of ribavirin plasma trough concentrations to steady state, from the same upstream popPK fit that supplies `CSS_RBV`. Used together with `CSS_RBV` inside `model()` to reconstruct the individual ribavirin concentration time-course analytically as `riba(t) = CSS_RBV * (1 - exp(-K_RBV * t))`.
- **Units:** 1/day (document per-model via `covariateData[[K_RBV]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- appears directly inside the analytical RBV-concentration expression alongside `CSS_RBV`.
- **Source aliases:**
  - `k_mode` (modal posterior estimate of individual ribavirin approach-to-Css rate constant from the upstream popPK fit, in 1/day) -- used in `Laouenan_2015_ribavirin.R` (DDMODEL00000285). Rename `k_mode` -> `K_RBV` before passing the dataset to `rxSolve`.
- **Example models:** `Laouenan_2015_ribavirin.R` (DDMODEL00000285; the bundle's `Simulated_Laouenant_2015_CPTPSP_hb_RBV.txt` carries `k_mode` values 0.013-0.47 day^-1 across the 15-subject cohort, corresponding to approach-to-Css half-lives of 1.5-55 days).
- **Notes:** Specific scope; same upstream-PK-dependency rationale as `CSS_RBV`. Distinct from any structural `kel` PK parameter -- `K_RBV` is an *apparent* approach-to-Css rate from a lumped exponential parameterization of the trough time-course, not the elimination-rate constant of a one-compartment IV model (the lumped form absorbs absorption, distribution, and elimination into a single first-order rate). Companion column: `CSS_RBV`.

### CP_MGL (**canonical for instantaneous drug plasma concentration as a time-varying PD driver**)
- **Description:** Instantaneous (per-event-record) plasma concentration of the modeled drug, supplied directly as a time-varying covariate column rather than computed from a coupled PK model. Used in PD-only myelosuppression / toxicity / response models that consume an upstream PK trajectory as input -- typically because the source analysis was done as a sequential PK-then-PD fit, with the PK model fixed from a previously published popPK analysis (e.g., a Bruno-style docetaxel popPK feeding a Friberg / Kloft myelosuppression PD model).
- **Units:** mg/L (= ug/mL; the two labels are numerically equivalent and the canonical reporting convention for cytotoxic-chemotherapy concentrations in Kloft 2006 / Friberg 2002 family analyses). Document per-model via `covariateData[[CP_MGL]]$units` if a different unit is used.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- typically enters as a linear term in the drug-effect expression (e.g., `drug = SL * CP_MGL`, where `SL` has units of `1/(mg/L)`). Set to 0 for placebo periods or any time outside the drug-exposure window. Reference values observed: docetaxel typical Cmax ~3 mg/L after a 100 mg/m^2 1-hour IV infusion (Kloft 2006 / Netterberg 2017 simulated dataset).
- **Source aliases:**
  - `CP` (Kloft 2006 / Netterberg 2017 NM-TRAN $INPUT convention for "predicted drug concentration"; values in mg/L) -- used in `Netterberg_2017_docetaxel.R`.
- **Example models:** `Netterberg_2017_docetaxel.R` (linear drug effect on the proliferation rate of the Friberg myelosuppression chain: `(1 - SL * CP_MGL)` with `SL = 19.27 (mg/L)^-1` after Kloft 2006's `THETA(3)/808*1000` MW-808 conversion; CP_MGL supplied per event row from an upstream docetaxel popPK simulation).
- **Notes:** Specific scope because the covariate's mechanistic meaning is bound to the modeled drug and the source paper's chosen PK input (e.g., a Bruno 1996/1998 docetaxel popPK trajectory for Netterberg 2017). Distinct from `CAV` (dosing-interval-averaged exposure used in Emax / EC50 PD models) -- `CP_MGL` is the instantaneous concentration, sampled at every PD event time. When a future paper requires the same time-varying-PK-as-PD-input pattern for a different drug, register a drug-specific canonical (e.g., `CP_PACL_MGL` for paclitaxel) rather than overloading this name; `CP_MGL` retains the implicit "drug = the modeled drug under the PD analysis" semantics. When a paper supplies the time-varying PK as separate per-subject empirical-Bayes PK parameters (e.g., `CL_INDIV`, `VC_INDIV`, `VP_INDIV`), use those columns in a coupled PK-PD ODE model (see `Friberg_2002_paclitaxel.R`) rather than reducing to `CP_MGL`. The choice between PK-as-covariate (this canonical) and PK-as-EBE-parameters depends on whether the source paper's NM-TRAN dataset shipped Cp directly or shipped the upstream individual PK parameters.

### CP_OXY_NGML (**canonical for instantaneous oxypurinol plasma concentration as a time-varying PD driver**)
- **Description:** Instantaneous plasma concentration of oxypurinol (the active metabolite of allopurinol; xanthine oxidase inhibitor) supplied directly as a time-varying covariate column rather than computed from a coupled PK model. Used in semi-mechanistic uric-acid disposition models that take XOI exposure as input to the production-inhibition equation.
- **Units:** ng/mL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- enters as a Hill term in the production-inhibition expression `1 - Rmax * CP_OXY_NGML / (CP_OXY_NGML + p50)`. Set to 0 outside the drug-exposure window or for non-XOI scenarios. Reference values observed: mean daily concentration on 300 mg/day allopurinol is approximately 10,000 ng/mL (Aksenov 2018, Eq. 13).
- **Source aliases:**
  - `[P]_PIN` (oxypurinol) -- the symbol used in Aksenov 2018 Eq. 9 for the production-inhibitor concentration when the inhibitor is oxypurinol.
- **Example models:** `Aksenov_2018_uricAcid.R` (Hill-type production inhibition with `rmax_oxy = 0.84` and `p50_oxy = 14000 ng/mL` per Aksenov 2018 Table 1).
- **Notes:** Specific scope because the canonical name is bound to oxypurinol; allopurinol's PK is conventionally summarized via the active metabolite oxypurinol (Day et al. 2007). Distinct from `CP_FBX_NGML` (febuxostat) and `CP_LSN_NGML` (lesinurad). When the source paper supplies an upstream popPK for oxypurinol (e.g., Wright et al. 2013, Anzai & Endou 2012), the user simulates that PK to populate this column; otherwise a steady-state value can be used. Ratified canonically on 2026-05-08 alongside the Aksenov 2018 extraction.

### CP_FBX_NGML (**canonical for instantaneous febuxostat plasma concentration as a time-varying PD driver**)
- **Description:** Instantaneous plasma concentration of febuxostat (xanthine oxidase inhibitor) supplied directly as a time-varying covariate column rather than computed from a coupled PK model. Used in semi-mechanistic uric-acid disposition models that take XOI exposure as input to the production-inhibition equation.
- **Units:** ng/mL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- enters as a Hill term in the production-inhibition expression `1 - Rmax * CP_FBX_NGML / (CP_FBX_NGML + p50)`. Set to 0 outside the drug-exposure window or for non-XOI scenarios. Reference values observed: mean daily concentration on 40 mg/day febuxostat is approximately 1000-2000 ng/mL (Aksenov 2018, Bhattaram & Gobburu 2017 regulatory review).
- **Source aliases:**
  - `[P]_PIN` (febuxostat) -- the symbol used in Aksenov 2018 Eq. 9 for the production-inhibitor concentration when the inhibitor is febuxostat.
- **Example models:** `Aksenov_2018_uricAcid.R` (Hill-type production inhibition with `rmax_fbx = 1` (fixed) and `p50_fbx = 120 ng/mL` for hyperuricemic subjects (or 87 ng/mL for normouricemic subjects) per Aksenov 2018 Table 1).
- **Notes:** Specific scope; febuxostat-specific. The `p50` parameter differs between hyperuricemic and normouricemic populations in Aksenov 2018; `Rmax` is fixed at 1 per Bhattaram & Gobburu 2017. Distinct from `CP_OXY_NGML` (oxypurinol) and `CP_LSN_NGML` (lesinurad). Ratified canonically on 2026-05-08 alongside the Aksenov 2018 extraction.

### CP_LSN_NGML (**canonical for instantaneous lesinurad plasma concentration as a time-varying PD driver**)
- **Description:** Instantaneous plasma concentration of lesinurad (uricosuric URAT1 inhibitor) supplied directly as a time-varying covariate column rather than computed from a coupled PK model. Used in semi-mechanistic uric-acid disposition models that take uricosuric exposure as input to the fractional-excretion-increase equation.
- **Units:** ng/mL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- enters as a Hill term in the fractional-excretion expression `FE = FE0 + Fmax * CP_LSN_NGML / (CP_LSN_NGML + p50)`. Set to 0 outside the drug-exposure window. Reference values observed: peak plasma concentration after single dose 200 mg lesinurad is approximately 6,000-9,000 ng/mL (Fleischmann et al. 2014; Shen et al. 2015).
- **Source aliases:**
  - `[P]_RIN` (lesinurad) -- the symbol used in Aksenov 2018 Eq. 10 for the reabsorption-inhibitor concentration when the inhibitor is lesinurad.
- **Example models:** `Aksenov_2018_uricAcid.R` (Hill-type increase in fractional excretion with `fmax_lsn = 0.56` (fixed) and `p50_lsn = 23000 ng/mL` for hyperuricemic subjects (or 11000 ng/mL for normouricemic subjects) per Aksenov 2018 Table 1).
- **Notes:** Specific scope; lesinurad-specific. The `p50` parameter differs between hyperuricemic and normouricemic populations in Aksenov 2018; `Fmax` was fixed during estimation. Distinct from `CP_OXY_NGML` (oxypurinol) and `CP_FBX_NGML` (febuxostat). Ratified canonically on 2026-05-08 alongside the Aksenov 2018 extraction.

### CP_MORPH_NGML (**canonical for instantaneous morphine plasma concentration as a time-varying PD driver**)
- **Description:** Instantaneous plasma concentration of morphine supplied directly as a time-varying covariate column rather than computed from a coupled PK model. Used in PD-only IRT / latent-pain models that take morphine exposure as an external input to a concentration-effect equation.
- **Units:** ng/mL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- enters linearly into the latent-pain equation `pain = pain_state - e_morph_pain * CP_MORPH_NGML + e_time_pain * time` (Valitalo 2017). Reference values observed: most individual predicted concentrations in Valitalo 2017 were within 0-60 ng/mL (Figure 2a); the IRT linear morphine-effect slope is 0.0091 (ng/mL)^-1, so a 20 ng/mL morphine exposure reduces the latent pain by ~0.18 latent-variable units.
- **Source aliases:**
  - `CP` (Valitalo 2017 NM-TRAN $INPUT convention for "morphine plasma concentration"; values in ng/mL) -- used in `Valitalo_2017_morphine.R`.
- **Example models:** `Valitalo_2017_morphine.R` (linear morphine concentration-effect on the IRT latent pain variable; CP_MORPH_NGML supplied per event row from an upstream morphine popPK simulation, typically `Knibbe_2009_morphine.R`).
- **Notes:** Specific scope; morphine-specific. The drug-specific naming follows the existing `CP_OXY_NGML` / `CP_FBX_NGML` / `CP_LSN_NGML` precedent established with Aksenov 2018. Distinct from the broader `CP_MGL` (mg/L PD-driver convention used in Netterberg 2017 docetaxel myelosuppression and similar) because the IRT PD models in this family use ng/mL natively. When a future morphine PD analysis uses mg/L, the conversion is `CP_MORPH_NGML = CP_MORPH_MGL * 1000`. Ratified canonically alongside the Valitalo 2017 morphine extraction (DDMODEL00000247).

### FPG (**canonical for baseline fasting plasma glucose**)
- **Description:** Fasting plasma glucose concentration at baseline (or time-varying baseline-style observation; document per-model). Distinct from `GLU` (time-varying plasma glucose regressor input for mechanistic glucose-kinetics models).
- **Units:** mmol/L (or mg/dL -- 1 mmol/L glucose is approximately 18.02 mg/dL). Document per-model via `covariateData[[FPG]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with linear-deviation form `1 + theta * (FPG - ref)`. Reference values observed: 8.90 mmol/L (Retlich 2015 popPK/PD linagliptin median fasting glucose at baseline).
- **Source aliases:** none known.
- **Example models:** `Retlich_2015_linagliptin.R` (mmol/L, reference 8.90; linear-deviation effect on baseline DPP-4 activity BSL with coefficient 1.46 % per mmol/L deviation).
- **Notes:** Glycemic-control covariate (baseline FPG); routinely reported alongside HbA1c in T2DM populations. Distinct from `GLU` which is a time-varying within-subject glucose regressor for mechanistic glucose-kinetics models (Bizzotto 2016). Ratified canonically alongside the Retlich 2015 linagliptin extraction.

### DPP4_BL_RFU (**canonical for baseline plasma dipeptidyl peptidase-4 activity in relative fluorescence units**)
- **Description:** Baseline plasma dipeptidyl peptidase-4 (DPP-4) enzymatic activity, measured by relative fluorescence units (RFU). DPP-4 is the pharmacological target of the gliptin (DPP-4-inhibitor) drug class; baseline activity correlates with circulating DPP-4 protein concentration in the central compartment and serves as a covariate on the central-compartment binding-site concentration in TMDD models for gliptins.
- **Units:** RFU (assay-specific; document the assay in `covariateData[[DPP4_BL_RFU]]$notes` since RFU values are not directly comparable across assays).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with linear-deviation form `1 + theta * (DPP4_BL_RFU - ref)`. Reference values observed: 12,497 RFU (Retlich 2015 popPK linagliptin median baseline), 11,600 RFU (Retlich 2015 popPK/PD linagliptin median baseline, applied to the individual-predicted BSL_i parameter on EC50).
- **Source aliases:** none known.
- **Example models:** `Retlich_2015_linagliptin.R` (RFU, reference 12,497; linear-deviation effect on the central-compartment binding-site concentration Bmax,C with coefficient 0.00332 % per RFU deviation -- captures the inter-individual correlation between baseline DPP-4 protein concentration and the apparent saturable-binding amplitude).
- **Notes:** Specific scope because DPP-4-activity values in RFU are assay-specific and not directly transferable between studies / instruments. Future gliptin extractions reporting DPP-4 activity in the same assay can reuse this canonical; extractions in absolute enzymatic-rate units (pmol AMC per minute) or normalised units should register a sibling canonical (e.g., `DPP4_BL_PMOL_MIN`). Ratified canonically alongside the Retlich 2015 linagliptin extraction.

### GLU (**canonical for plasma glucose time-course regressor**)
- **Description:** Plasma glucose concentration as a time-varying *regressor* input that drives a mechanistic glucose-kinetics model. Not a covariate that modifies a parameter; the model integrates `GLU` directly through a smoothing filter into a site-of-action glucose variable.
- **Units:** mmol/L
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used as a time-varying regressor. The model declares `linear(GLU)` so rxode2 linearly interpolates `GLU` between dataset rows.
- **Source aliases:**
  - `iglu` (glucose at the current row time) -- used in the DDMORE bundle's `Simulated_glucoseKinetics.csv` for `DDMODEL00000227`. Rename `iglu` -> `GLU` before passing to `rxSolve`.
- **Example models:** `Bizzotto_2016_glucose.R` (driving regressor for the glucose-at-site-of-action delay).
- **Notes:** Specific scope because `GLU` is meaningful only for glucose-kinetics or glucose-PD models that take plasma glucose as an exogenous regressor. The DDMORE bundle's hand-rolled piecewise-linear interpolation (`GL = (t-T1)/(TOBS-T1)*(GLU-GLU1)+GLU1` with bracketing columns `iglu / glun / td / tn`) is replaced in nlmixr2 by `linear(GLU)` declared in `model()`; the bracketing columns are not required.

### INS (**canonical for plasma insulin time-course regressor**)
- **Description:** Plasma insulin concentration as a time-varying *regressor* input that drives a mechanistic glucose-kinetics model. Not a covariate that modifies a parameter; the model integrates `INS` directly through a smoothing filter into a site-of-action insulin variable.
- **Units:** pmol/L
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used as a time-varying regressor. The model declares `linear(INS)` so rxode2 linearly interpolates `INS` between dataset rows.
- **Source aliases:**
  - `iins` (insulin at the current row time) -- used in the DDMORE bundle's `Simulated_glucoseKinetics.csv` for `DDMODEL00000227`. Rename `iins` -> `INS` before passing to `rxSolve`.
  - `INSU` -- used in the DDMORE bundle's `Simulated_ddmoremockdata2.txt` for `DDMODEL00000228`. Rename `INSU` -> `INS` before passing to `rxSolve`.
- **Example models:** `Bizzotto_2016_glucose.R` (driving regressor for the insulin-at-site-of-action delay), `NA_NA_paracetamol.R` (DDMODEL00000228 OGTT model: drives the insulin-on-glucose-elimination first-order effect compartment via `kie * (INS / 6.945 - effect_ins)`).
- **Notes:** Specific scope because `INS` is meaningful only for glucose-kinetics or insulin-PD models that take plasma insulin as an exogenous regressor. For drugs that *modify* circulating insulin as a downstream effect, use a different mechanism-specific name. The DDMORE bundle's hand-rolled piecewise-linear interpolation (`I = (t-T1)/(TOBS-T1)*(INS-INS1)+INS1` with bracketing columns `iins / insn / td / tn`) is replaced in nlmixr2 by `linear(INS)` declared in `model()`; the bracketing columns are not required.

### INS_BL (**canonical for baseline (fasting) plasma insulin concentration**)
- **Description:** Baseline (fasting) plasma insulin concentration, time-fixed per subject. Used as the per-subject anchor for steady-state insulin-driven processes (e.g., initial condition of an insulin-on-elimination effect compartment) and as the per-subject baseline-state insulin input to a baseline-glucose-production rate calculation.
- **Units:** pmol/L (or uU/mL; document per-model via `covariateData[[INS_BL]]$units`). The example model rescales via `INS_BL / 6.945` to convert pmol/L to uU/mL.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a.
- **Source aliases:**
  - `BASI` (baseline insulin) -- used in the DDMORE bundle's `Simulated_ddmoremockdata2.txt` for `DDMODEL00000228`. Rename `BASI` -> `INS_BL` before passing to `rxSolve`.
- **Example models:** `NA_NA_paracetamol.R` (DDMODEL00000228 OGTT model: initialises the insulin-on-elimination effect compartment `effect_ins(0) = INS_BL / 6.945` and feeds the steady-state baseline-glucose-production rate `gpro = gss * (kg + kgi * INS_BL / 6.945) * vg * 180 / 1000`).
- **Notes:** Distinct from `INS` (time-varying regressor); `INS_BL` is a per-subject baseline-state anchor used in initial conditions and steady-state derived quantities, not the dynamic regressor itself. Specific scope because the conversion factor (1/6.945) and the rescaled-units interpretation are paper-specific; future extractions that report baseline insulin in mIU/L or pmol/L directly without rescaling can ratify the same canonical and document the per-model units / conversion in `covariateData[[INS_BL]]$units` / `notes`. Companion concept to `FPG` (baseline fasting plasma glucose).

### IGE (**canonical for serum total immunoglobulin E concentration**)
- **Description:** Baseline serum total immunoglobulin E concentration (free IgE plus, in patients on anti-IgE therapy, omalizumab-IgE complex). For anti-IgE monoclonal antibodies (omalizumab, ligelizumab) IgE is the pharmacologic target; baseline IgE sets the magnitude of the target sink and modifies free-IgE clearance and the rate of IgE production in mechanism-based binding/turnover models.
- **Units:** ng/mL (typical clinical-PK convention). Pretreatment values reported in `IU/mL` are converted via `1 IU/mL = 2.42 ng/mL` (Hayashi 2007 Methods). Document per-model via `covariateData[[IGE]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(IGE / ref)^exponent`. Reference value observed: 482.4 ng/mL (Hayashi 2007 Japanese atopic-asthma cohort).
- **Source aliases:**
  - `IgE0` (baseline IgE concentration) -- used in `Hayashi_2007_omalizumab.R`.
- **Example models:** `Hayashi_2007_omalizumab.R` (ng/mL, reference 482.4; power exponents -0.281 on apparent CL of free IgE and +0.657 on apparent IgE production rate; also used as the initial value for the total-IgE state at t = 0).
- **Notes:** General scope because baseline serum total IgE is a routine clinical-laboratory measurement, not a target tied to one drug. In mechanism-based anti-IgE binding/turnover models the in-model IgE state is a separate dynamic variable (`X_TE`, in nmol or nmol/L) -- `IGE` is the per-subject baseline column used for covariate scaling and (when applicable) state initialization, not the dynamic state itself. For models that use the alternative reporting unit `IU/mL`, multiply by 2.42 before applying the canonical-units (ng/mL) reference value, or document the per-model unit choice in `covariateData[[IGE]]$units` so downstream tooling can interpret the values correctly.

## Count / Markov-feedback PD covariates

These columns are specific to count / Markov / time-to-event PD models that
fit per-record event counts (e.g., daily or monthly seizure counts) with
optional dependence on the previous-period count. Register names retain
their source-paper conventions where those names are unambiguous and
readable.

### TRT_PHASE (**canonical for double-blind active-treatment-phase indicator**)
- **Description:** 1 = the record falls within the active double-blind treatment phase (drug + placebo effects are switched on); 0 = baseline / run-in / off-treatment (drug + placebo effects are zeroed).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (baseline / off-treatment).
- **Source aliases:**
  - `Q2` -- used in the Schoemaker 2018 LEV / BRV pediatric extrapolation (DDMODEL00000239) as the treatment-phase gating multiplier on the combined placebo + drug-effect log-rate term (`LE = LS0 + Q2 * LTRTE`).
- **Example models:** `Schoemaker_2018_levetiracetam.R` (DDMODEL00000239).
- **Notes:** New canonical because the source name `Q2` collides with the canonical PK parameter `q2` (inter-compartmental clearance to peripheral2) -- `q2` could not be used as a covariate column without confusing source-trace lookups. Useful for any phase-gated PD model where placebo and drug effects are constrained to a specific study period; analogous to a "treatment-on" indicator.

### PDV (**canonical for previous-period observed seizure count**)
- **Description:** Previous-period observed seizure count, supplied per-record as a covariate input to capture Markov dependence on the previous count. For daily-count datasets PDV is the observation on the immediately preceding day; for monthly (or other aggregate) records the previous-interval count would apply if used.
- **Units:** (count, integer)
- **Type:** count
- **Scope:** specific
- **Reference category:** n/a -- typically used as `<param> * PDV / (ES50 + PDV)` (Markov-Hill) or analogous form.
- **Source aliases:** `PDV` -- used in the Schoemaker 2018 LEV/BRV pediatric extrapolation (DDMODEL00000239) and in the Ahn 2010 pregabalin Markov seizure-count model that the Schoemaker 2018 publication cites as the precursor structure.
- **Example models:** `Schoemaker_2018_levetiracetam.R` (DDMODEL00000239).
- **Notes:** Specific scope because the column's interpretation (a Markov feedback state) is intrinsic to count-likelihood Markov-feedback PD models. nlmixr2 / rxode2 cannot natively express a Markov dependence of an observation on the immediately preceding observation as a model state -- supplying PDV as a per-record data column is the operator-approved (sidecar response-001 Q2) way to preserve the published structure in this batch. For records where the Markov term should not contribute (e.g., monthly counts in a mixed daily / monthly cohort) the bundle convention is to set PDV = -99 as a sentinel; in `Schoemaker_2018_levetiracetam.R` the model gates the Markov contribution on `CHILD = 1` so the sentinel is multiplied by 0 and is harmless.

### NDAYS (**canonical for number of days in the count-record interval**)
- **Description:** Number of days in the interval over which the observed seizure count was tabulated. Multiplies the per-day rate to give the expected count for the record (e.g., `LAMB = exp(LE) * NDAYS`).
- **Units:** days
- **Type:** count
- **Scope:** general
- **Reference category:** n/a -- appears as a multiplier on a per-day rate.
- **Source aliases:** `NDAYS` -- used in the Schoemaker 2018 LEV/BRV pediatric extrapolation (DDMODEL00000239).
- **Example models:** `Schoemaker_2018_levetiracetam.R` (DDMODEL00000239).
- **Notes:** General scope because the count-interval-length concept is shared across any count or rate-based PD model that mixes record granularity (e.g., daily and monthly counts in the same dataset). For pure-daily or pure-monthly cohorts the column is constant; including it as a covariate keeps the model usable in mixed-granularity simulations.

### MOMENT (**canonical for endotracheal suctioning procedural state**)
- **Description:** Procedural state at the time of the pain / distress assessment in invasive-ventilation studies, used to select between separate pre-procedure / intra-procedure / post-procedure baseline parameters. Coded 1 = before suctioning, 2 = during suctioning, 3 = after suctioning.
- **Units:** (categorical, 3 levels)
- **Type:** categorical
- **Scope:** specific
- **Reference category:** 1 (before suctioning).
- **Source aliases:** `MOMENT` -- used in the Valitalo 2017 IRT morphine PD model (DDMODEL00000247).
- **Example models:** `Valitalo_2017_morphine.R` (DDMODEL00000247; selects between the three baseline pain typical values `presuct` / `suct` / `aftsuct` and the matching 3x3 correlated etas).
- **Notes:** Specific scope because the column's meaning is tied to a particular study procedure (endotracheal suctioning during mechanical ventilation in neonates). The Simons 2003 cohort that Valitalo 2017 re-analysed scheduled pain assessments around suctioning events, so MOMENT changes within-subject at each scheduled assessment. Ratified canonically alongside the Valitalo 2017 morphine extraction.

### ITEM (**canonical for pain-assessment item identifier in IRT graded-response models**)
- **Description:** Identifier of the specific pain-assessment item being scored at each observation row, used to dispatch between the IRT discrimination / difficulty parameter sets in a graded-response model. Valitalo 2017 coding: 1 = COMFORT-B alertness; 2 = COMFORT-B calmness/agitation; 3 = COMFORT-B respiratory response; 5 = COMFORT-B body movement; 7 = COMFORT-B facial tension; 12 = VAS (cm, range 0-10); 25 = PIPP brow bulge; 26 = PIPP eye squeeze; 27 = PIPP nasolabial furrow; 28 = NIPS total.
- **Units:** (categorical, 9-10 levels depending on cohort)
- **Type:** categorical
- **Scope:** specific
- **Reference category:** n/a -- selects per-item parameter sets rather than acting as a reference contrast.
- **Source aliases:** `ITEM` -- used in the Valitalo 2017 IRT morphine PD model (DDMODEL00000247).
- **Example models:** `Valitalo_2017_morphine.R` (DDMODEL00000247; switches the IRT graded-response discrimination / difficulty parameters per row).
- **Notes:** Specific scope because the integer-to-item mapping is tied to the Valitalo 2017 NM-TRAN dataset's coding. Other IRT graded-response models in the library (when they are added) may use different integer codings and should register their own canonical (e.g., `ITEM_<study>`) when the codings collide. The COMFORT-B "muscle tension" item is omitted from the Valitalo 2017 coding because it could not be assessed from video recordings.

### OBSTYPE (**canonical for VAS observer type**)
- **Description:** Observer type for visual-analogue-scale pain assessments: 1 = investigator (video-based), 2 = bedside nurse. Used to select between observer-specific VAS difficulty / discrimination THETAs and between observer-specific residual-error SDs.
- **Units:** (categorical, 2 levels)
- **Type:** categorical
- **Scope:** specific
- **Reference category:** 1 (video investigator).
- **Source aliases:** `OBSTYPE` -- used in the Valitalo 2017 IRT morphine PD model (DDMODEL00000247).
- **Example models:** `Valitalo_2017_morphine.R` (DDMODEL00000247; selects between `diff_vas_video` vs `diff_vas_bedside`, `discr_vas_video` vs `discr_vas_bedside`, and `addSd_vas_video` vs `addSd_vas_bedside`).
- **Notes:** Specific scope because the binary observer-coding is tied to the Valitalo 2017 NM-TRAN dataset. Future IRT models with a different observer split (e.g., parent / nurse / physician) should register their own canonical rather than overloading this 2-level coding.

## Race / ethnicity

**Canonical pattern: `RACE_<GROUP>`.** Use one indicator per race/ethnicity group the source models. Reference category is the implicit 0 = all other groups; document explicitly which groups are in the reference. When the source uses composite groups (e.g., "Black or Other"), name them accordingly (`RACE_BLACK_OTHER`) and list the components in `notes`. The base `RACE_<GROUP>` indicators are scope: general; composite groupings are scope: specific because the grouping is tied to the study's analysis plan.

### RACE_BLACK (**canonical for Black / African American race indicator**)
- **Description:** 1 = Black / African American, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (document the actual reference groups used).
- **Source aliases:**
  - `BLACK` -- used in `Hu_2026_clesrovimab.R`, `Robbie_2012_palivizumab.R`.
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical form), `Robbie_2012_palivizumab.R`.

### RACE_WHITE (**canonical for White race indicator**)
- **Description:** 1 = White, 0 = non-White.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (non-White; complement composition depends on the source paper, typically pooling Black/African American, Asian, American Indian/Alaska Native, Native Hawaiian/Pacific Islander, Other, Not reported, Unknown). Some papers (e.g., Hu 2014) instead use the Caucasian (RACE_WHITE = 1) subgroup as the typical-value reference; the column encoding is unchanged but the model implements the effect on `(1 - RACE_WHITE)`.
- **Source aliases:**
  - `RACE` (with values `1 = White / 0 = non-White`) -- used in `Lin_2024_casirivimab.R`. Source column name `RACE` is generic; the canonical name is intentionally explicit because some other models use `RACE` for a different dichotomy.
  - `RACE` (Caucasian-vs-non-Caucasian dichotomy as named in Hu 2014 Table 2) -- used in `Hu_2014_bapineuzumab.R`. Same canonical column name and 1 = White / 0 = non-White encoding; the typical-value reference is the Caucasian subgroup, so the model implements the 15% non-Caucasian effect on `(1 - RACE_WHITE)`.
- **Example models:** `Lin_2024_casirivimab.R` (multiplicative fractional change on CL relative to non-White reference), `Hu_2014_bapineuzumab.R` (multiplicative 15% increase in CL for non-Caucasian relative to Caucasian reference).
- **Notes:** Used by papers that dichotomize race as White vs. non-White rather than decomposing into separate group indicators. Sign and reference-category interpretation are inverted relative to `RACE_BLACK` / `RACE_ASIAN` / etc.; do NOT combine `RACE_WHITE` with the decomposed indicators in the same model. The model's typical-value reference category (which subgroup gets the unmodified `lcl` / `lvc`) varies between papers -- Lin 2024 uses non-White as the reference, Hu 2014 uses Caucasian (White) as the reference; both share the same canonical column encoding.

### RACE_BLACK_OTH (**canonical for composite Black/Other group**)
- **Description:** 1 = Black/African American or Other race, 0 = other groups.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = White or Native Hawaiian/Pacific Islander (Clegg 2024 grouping).
- **Source aliases:** `BLACK_OTH` -- used in `Clegg_2024_nirsevimab.R`.
- **Example models:** `Clegg_2024_nirsevimab.R`.
- **Notes:** Kept distinct from `RACE_BLACK` because the composite is not interchangeable.

### RACE_ASIAN (**canonical for Asian race indicator**)
- **Description:** 1 = Asian, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Source aliases:** `ASIAN` -- used in `Hu_2026_clesrovimab.R`, `Robbie_2012_palivizumab.R`, `Fau_2020_isatuximab.R`. `RAAS` (race-Asian-vs-other indicator as named in Bajaj 2017 Table 1) -- used in `Bajaj_2017_nivolumab.R`. `RACEN` (race-numeric indicator as named in Lu 2019 / Shi 2020 NONMEM control stream; ASIAN = 1 if RACEN == 1) -- used in `Lu_2019_polatuzumab.R`.
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical form), `Robbie_2012_palivizumab.R`, `Bajaj_2017_nivolumab.R`, `Fau_2020_isatuximab.R`, `Lu_2019_polatuzumab.R` (multiplicative factor `e_asian_vc = 0.929` on acMMAE Vc, i.e., 7.1% lower V1 in Asian patients; verbatim Shi 2020 (PMID 32770353) ethnicity-sensitivity re-quote of the Lu 2019 popPK Asian-race covariate).

### RACE_ASIAN_AMIND_OTH (**canonical for composite Asian / American Indian / Other group**)
- **Description:** 1 = Asian, American Indian / Alaska Native, or Other race; 0 = White or Black. Composite indicator that pools the smaller-N race groups in a population dominated by White and Black subjects, with White + Black serving as the reference category. Distinct from `RACE_ASIAN_AMIND_MULTI` (Clegg 2024 grouping that includes Multiracial; pooled against a different reference), from `RACE_ASIAN_OTH` (within-Asian-population sub-indicator), and from `RACE_BLACK_OTH` (different composite).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = White or Black (the larger-N pooled group used as the reference in the source paper).
- **Source aliases:** none formally; Frey 2013's NONMEM control stream uses the inline race classification rather than a separate named column.
- **Example models:** `Frey_2013_tocilizumab.R` (multiplicative fractional effect on the DAS28 first-order loss rate Kout: `Kout * (1 - 0.25 * RACE_ASIAN_AMIND_OTH)` -- Kout is 25% lower in the Asian/AmInd/Other composite group relative to the White+Black reference).
- **Notes:** Specific scope because the composite grouping is defined by the source paper's analysis plan rather than by a uniform external standard. Do not combine with the decomposed `RACE_ASIAN`, `RACE_OTHER`, etc. indicators in the same model; the composite indicator is mutually exclusive with the decomposition. Ratified canonically on 2026-04-29 in support of the Frey 2013 tocilizumab DAS28 PKPD model. Frey 2013 uses TWO distinct race covariates: this `RACE_ASIAN_AMIND_OTH` indicator on Kout (DAS28-PD-side; pools Asian + AmInd + Other vs White+Black) AND the within-Asian `RACE_ASIAN_OTH` indicator on CL (PK-side; isolates the "Other Asian" subgroup within the Asian-only cohort).

### RACE_ASIAN_AMIND_MULTI (**canonical for composite group**)
- **Description:** 1 = Asian, American Indian / Alaskan Native, or Multiple races, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = White, Black / African American, Native Hawaiian / Pacific Islander, or Other (Clegg 2024 grouping).
- **Source aliases:** `ASIAN_AMIND_MULTI` -- used in `Clegg_2024_nirsevimab.R`.
- **Example models:** `Clegg_2024_nirsevimab.R`.
- **Notes:** Clegg 2024 applies this covariate to both CL and V2 with different coefficients.

### RACE_ASIAN_OTH (**canonical for Asian-other composite race indicator**)
- **Description:** 1 = subject self-identifies as Asian-other (Asian heritage outside the locally-dominant Asian subgroup, e.g. non-Chinese in a Chinese-dominant cohort, or "Other Asian" as a study-defined catch-all category). 0 = otherwise. Reference category is the cohort's dominant race grouping (typically Chinese or White, depending on the study).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0. Document the dominant subgroup explicitly in `covariateData[[RACE_ASIAN_OTH]]$notes` for every model that uses this covariate.
- **Source aliases:** none.
- **Example models:** `Brown_2017_osimertinib.R` (paper's "Asian (not Japanese or Chinese)" composite indicator with linear effect on apparent clearance of the AZ5104 metabolite; reference category Caucasian).
- **Notes:** Distinct from `RACE_ASIAN_AMIND_MULTI` (a 4-way composite of Asian + American Indian + Multiple Races), `RACE_ASIAN_AMIND_OTH` (a 3-way Asian + AmInd + Other composite against a White+Black reference, used in Frey 2013), and `RACE_BLACK_OTH` (different composite). `RACE_ASIAN_OTH` is a within-Asian-population sub-indicator, not a multi-race composite. Operator decision (2026-04-28): kept separate from `RACE_ASIAN` because the paper's "Other Asian" category is its own grouping, not an alias of "Asian (any)". Brown 2017 uses Caucasian (not Chinese) as the dominant reference.

### RACE_NEAS (**canonical for North East Asian composite race indicator**)
- **Description:** 1 = North East Asian heritage (worldwide Chinese, Japanese, or Korean), 0 = non-North East Asian. Composite indicator analogous to `RACE_ASIAN` but specifically restricted to the East Asian subgroup most-relevant to ICH E5 ethnic-sensitivity / Asian-region bridging analyses.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (any non-North East Asian race, including South / Southeast Asian, White, Black, etc.).
- **Source aliases:**
  - `RAC4` -- used in `Zhou_2021_belimumab.R` (Zhou 2021 Table 2 footnote d).
- **Example models:** `Zhou_2021_belimumab.R` (multiplicative factor 1.07 on V1).
- **Notes:** Distinct from the broader `RACE_ASIAN` (which can include South / Southeast Asian populations) because Zhou 2021 specifically tested whether Chinese/Japanese/Korean patients had different PK from the rest of the dataset; the analysis explicitly compared `RAC4` (North East Asian) against alternative race definitions and chose `RAC4` by AIC.

### RACE_MULTI (**canonical for multiracial indicator**)
- **Description:** 1 = multiracial, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Source aliases:** `MULTIRACIAL` -- used in `Hu_2026_clesrovimab.R`.
- **Example models:** `Hu_2026_clesrovimab.R`.

### RACE_OTHER (**canonical for race-category 'Other' indicator**)
- **Description:** 1 = race category "Other," 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Source aliases:**
  - `OTHER` -- used in `Robbie_2012_palivizumab.R`.
- **Example models:** `Zhu_2017_lebrikizumab.R`, `Robbie_2012_palivizumab.R`.

### RACE_HISPANIC (**canonical for Hispanic / Latino ethnicity indicator**)
- **Description:** 1 = Hispanic / Latino, 0 = non-Hispanic. Used by papers that report Hispanic as a separate category alongside Black, Asian, and Other rather than as a distinct ethnicity dimension.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (non-Hispanic; document the paper-specific reference race composition per-model).
- **Source aliases:**
  - `HISPANIC` -- used in `Robbie_2012_palivizumab.R`.
- **Example models:** `Robbie_2012_palivizumab.R` (fractional effect on CL; additional effect on Vc).
- **Notes:** In the US Office-of-Management-and-Budget (OMB) classification Hispanic is an ethnicity rather than a race, but clinical PK analyses frequently treat it as one of the race indicators. When a paper treats Hispanic as a race, use this column; otherwise encode ethnicity separately. Register-wise, this follows the `RACE_<GROUP>` indicator-decomposition pattern.

### RACE_JAPANESE (**canonical for Japanese-heritage race indicator**)
- **Description:** 1 = Japanese heritage, 0 = non-Japanese. Used when Japanese subjects form a distinct subgroup in the study design (e.g., ICH E5 bridging analyses or studies with a dedicated Japanese healthy-volunteer cohort).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (non-Japanese).
- **Source aliases:**
  - `JAPANESE_HV` -- used in `Wang_2017_benralizumab.R` (Japanese healthy-volunteer cohort indicator; the healthy-volunteer vs. asthma-patient distinction is captured separately, not in this covariate).
- **Example models:** `Wade_2015_certolizumab.R` (multiplicative fractional effect on V/F; Wade 2015 breaks Japanese [RACE.EQ.8] out separately from RACE_ASIAN), `Wang_2017_benralizumab.R` (multiplicative factor 1.34 on Vc).
- **Notes:** Distinct from `RACE_NEAS` (North East Asian composite, includes Chinese, Japanese, and Korean) and from `RACE_ASIAN`. Use `RACE_JAPANESE` only when the source paper breaks out Japanese heritage as its own indicator; do not aggregate with other Asian groups when the paper keeps them separate. Ratified canonically on 2026-04-26.

### RACE_CHINESE (**canonical for Chinese-heritage race indicator**)
- **Description:** 1 = Chinese heritage, 0 = non-Chinese. Used when Chinese subjects form a distinct subgroup alongside Japanese, Asian-other, and other race categories (e.g., multiregional oncology trials enrolling Chinese, Japanese, and other Asian cohorts as separate strata).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (non-Chinese; document the paper-specific reference race composition per-model).
- **Source aliases:** none formally; source NONMEM control streams typically use ad-hoc names (e.g., `CHINESE`, `RACE.EQ.X`).
- **Example models:** `Brown_2017_osimertinib.R` (linear additive effect `(1 + 0.17 * RACE_CHINESE)` on apparent clearance of the AZ5104 metabolite; reference category Caucasian).
- **Notes:** Distinct from `RACE_NEAS` (North East Asian composite, includes Chinese, Japanese, and Korean) and from `RACE_ASIAN`. Use `RACE_CHINESE` only when the source paper breaks out Chinese heritage as its own indicator alongside `RACE_JAPANESE` and `RACE_ASIAN_OTH`; do not aggregate with other Asian groups when the paper keeps them separate. Parallels the established `RACE_JAPANESE` entry. Ratified canonically on 2026-05-09.

## Geographic / enrollment-country indicators

Geographical study-site region indicators. Distinct from race / ethnicity (`RACE_*`), which describe subject ancestry; these describe the geographical location of the clinical trial site that enrolled the subject. Used in multi-regional studies (typically those including bridging analyses for Japan or East Asia) to capture region-specific clinical-practice or unmeasured-environment effects on PK that remain after accounting for body weight, race, and laboratory covariates. Encoded as a set of mutually exclusive binary indicators with US as the implicit reference category (all indicators = 0). When a paper groups some non-US regions with US (e.g., Hong 2025 groups US and Japan as the DXd CL reference), the model code uses only the indicators that distinguish the non-reference groups; the data column for the grouped region (e.g., `REGION_JAPAN`) is still recorded so the same dataset can serve other parameters that do separate that group.

### REGION_JAPAN (**canonical for Japan study-site / enrollment-country indicator**)
- **Description:** 1 = study site in Japan (or patient enrolled in Japan, depending on source paper's reporting), 0 = study site / enrollment country outside Japan. Geographical Japan indicator.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-Japan sites / enrollment; specific reference set varies per model -- e.g., Hong 2025 Dato-DXd CL uses "any non-Japan region", whereas Hong 2025 DXd CL groups Japan with US into the reference).
- **Source aliases:**
  - `COUNTRY_JPN` -- retired canonical; used in `Yin_2021_trastuzumabDeruxtecan.R` as an enrollment-country (not study-site-region) indicator. Some papers report country-of-enrollment rather than site region; both map to `REGION_JAPAN` when the binary contrast is Japan vs. non-Japan.
- **Example models:** `Hong_2025_datopotamab.R` (multiplicative effect 1 + (-0.219) = 0.781 on Dato-DXd linear clearance), `Yin_2021_trastuzumabDeruxtecan.R` (multiplicative effect 0.903 on CL_intact and 0.738 on V2_intact when REGION_JAPAN = 1; Yin 2021 retained Japan enrollment-country over Japanese race because the two were highly confounded, correlation -0.81).
- **Notes:** Distinct from `RACE_JAPANESE` (subject ancestry). A subject of Japanese ancestry enrolled at a US site has `RACE_JAPANESE = 1` but `REGION_JAPAN = 0`. Some papers (e.g., Yin 2021) report enrollment country rather than study-site region; both are encoded as `REGION_JAPAN` when the binary contrast is Japan vs. non-Japan. Paired with `REGION_EUROPE` and `REGION_ROW` in multi-regional studies (e.g., Hong 2025); `REGION_JAPAN = 0` for a US-only cohort.

### REGION_EUROPE (**canonical for Europe study-site indicator**)
- **Description:** 1 = study site in Europe, 0 = study site outside Europe. Geographical study-site region.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-Europe study sites; specific reference set varies per model).
- **Source aliases:** none yet; canonical name preferred.
- **Example models:** `Hong_2025_datopotamab.R` (multiplicative effect 1 + 0.240 = 1.240 on DXd clearance versus US/Japan reference).
- **Notes:** Pair with `REGION_JAPAN` and `REGION_ROW` to encode multi-regional study membership; subjects with all three indicators = 0 are in the "US" reference group.

### REGION_ROW (**canonical for Rest-of-World study-site indicator**)
- **Description:** 1 = study site in Rest of World (i.e., not US, Japan, or Europe), 0 = study site in US / Japan / Europe. Residual-region indicator for multi-regional studies.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (US / Japan / Europe study sites; specific reference set varies per model).
- **Source aliases:** none yet; canonical name preferred.
- **Example models:** `Hong_2025_datopotamab.R` (multiplicative effect 1 + 0.196 = 1.196 on DXd clearance versus US/Japan reference).
- **Notes:** "Rest of the World" composition is paper-specific (e.g., Hong 2025 = study sites outside US, Japan, and Europe). Document the subject set in `covariateData[[REGION_ROW]]$notes` per model.

## Pediatric comorbidities

### CLD_PREM (**canonical for chronic lung disease of prematurity**)
- **Description:** 1 = chronic lung disease of prematurity (bronchopulmonary dysplasia, BPD), 0 = no CLD. Time-fixed per subject (diagnosis at study entry).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no CLD of prematurity).
- **Source aliases:**
  - `CLD` -- used in `Robbie_2012_palivizumab.R`.
  - `BPD` -- bronchopulmonary-dysplasia shorthand.
- **Example models:** `Robbie_2012_palivizumab.R` (fractional +20% effect on CL).
- **Notes:** Standard pediatric / neonatology comorbidity flag; ties to palivizumab's label population (high-risk preterm infants) and may re-appear in future pediatric mAb PK analyses (RSV, parenteral nutrition, etc.).

## Comorbidities

### DIAB (**canonical for diabetes-mellitus comorbidity indicator**)
- **Description:** 1 = patient has diabetes mellitus comorbidity (Type 1 or Type 2 not distinguished), 0 = no diabetes comorbidity. Time-fixed at study entry per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no diabetes comorbidity).
- **Source aliases:**
  - `DIAB` -- used in `Chen_2022_guselkumab.R`.
  - `Diabetes` -- used in `Sherer_2012_AAA.R` (Sherer 2012 Methods page 2 symbol "Diabetes").
- **Example models:**
  - `Chen_2022_guselkumab.R` (multiplicative effect on CL/F: 1.15^DIAB, +15% in patients with diabetes).
  - `Sherer_2012_AAA.R` (additive shift on the first derivative of AAA growth rate with size beta2: `e_diab_b2 = -0.32/year` for diabetics; cohort prevalence 14%).
- **Notes:** Captures pre-existing diabetes mellitus as a comorbidity in non-diabetes-primary indications (e.g., psoriatic arthritis, psoriasis, vascular disease). Distinct from a primary disease-state indicator like `DIS_UC`. Type 1 vs Type 2 mellitus is not separated unless the source paper distinguishes them; in pooled-population PK analyses, the covariate is typically a single binary flag derived from medical history. Diabetic patients tend to have higher inflammation and altered IgG turnover, which can manifest as modest changes in monoclonal-antibody clearance. In vascular populations (Sherer 2012) diabetes is associated with slower AAA growth, possibly via aberrant monocyte-matrix interactions (Golledge 2008 mechanism cited in Sherer 2012 Discussion).

### T2DM (**canonical for type-2-diabetes-mellitus-specific indicator**)
- **Description:** 1 = patient has Type-2 diabetes mellitus specifically; 0 = normal-glucose-tolerance control (or other reference cohort pooled in the source analysis). Time-fixed at study entry per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (normal-glucose-tolerance control).
- **Source aliases:**
  - `T2DM` -- used in `NA_NA_paracetamol.R` (DDMODEL00000228).
  - `T2D` -- used in `Guiastrennec_2016_gastric_emptying.R` (matched-cohort flag, 1 = T2D patient vs 0 = matched nondiabetic control).
- **Example models:**
  - `NA_NA_paracetamol.R` (DDMODEL00000228).
  - `Guiastrennec_2016_gastric_emptying.R` (multiplicative -81.1% depression of POTcarbC, the carbohydrate potency on CCK release; all other parameters are common across cohorts).
- **Notes:** Distinct from the existing `DIAB` canonical (which deliberately does not distinguish Type 1 vs Type 2). Specific scope because the reference cohort is study-specific and the mechanism in the example model is a Type-2-versus-healthy stratification of OGTT response; a future T2DM-specific study (e.g., a popPK/PD analysis stratifying by HbA1c level) can ratify the same canonical and document the reference cohort in `covariateData[[T2DM]]$notes`.

### HYPERT (**canonical for hypertension comorbidity / medical-history indicator**)
- **Description:** 1 = patient has a history of (or current) hypertension as a comorbidity; 0 = no hypertension. Time-fixed at study entry per subject (medical-history flag rather than time-varying blood-pressure measurement).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no hypertension comorbidity).
- **Source aliases:**
  - `MHHY` (medical history of hypertension) -- used in `Girard_2012_pimasertib.R`.
- **Example models:** `Girard_2012_pimasertib.R` (additive shift on the cumulative-logit AE-score model: `theta_mhhy * HYPERT`; +0.539 logit units in patients with prior hypertension).
- **Notes:** Companion to `DIAB` (diabetes-mellitus comorbidity); both are baseline binary medical-history flags collected from clinical-history forms. Captures any prior or current hypertension diagnosis, regardless of treatment status; if a future model needs to separate treated vs untreated hypertension, register a refinement (`HYPERT_TREATED`).

## Surgical history / disease state

### POD (**canonical for post-operative day**)
- **Description:** Days elapsed since the qualifying surgical event (e.g., solid-organ transplantation, major resection). Time-varying within subject; integer- or fractional-day valued; rises monotonically from 0 at the day of surgery. Captures time-since-surgery effects on PK that are not explained by other covariates -- e.g., post-transplant clearance of immunosuppressants typically declines toward a steady value over the first weeks-to-months as graft function, fluid status, hematocrit, and corticosteroid taper stabilise.
- **Units:** days
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- typically used in centred-deviation form `(1 + e_pod_param * (POD - ref_pod))`, sometimes with an upper cap (e.g., values > 180 fixed to 180 days when the residual time-varying effect plateaus).
- **Source aliases:**
  - `POD` -- used in `Bergmann_2014_tacrolimus.R` (capped at 180 days; centred at 22.7 days, the dataset mean).
- **Example models:** `Bergmann_2014_tacrolimus.R` (linear deviation from POD = 22.7 days on tacrolimus CL/F; coefficient -0.0021 per day implies a 0.21% per-day decrease in apparent oral clearance with a 180-day plateau).
- **Notes:** Time-varying within subject; the per-row value is the integer day count from the date of surgery to the observation date. For solid-organ-transplant cohorts, `POD` is the conventional NONMEM `$INPUT` column name. When the source paper reports a different name (`DPT` for "days post-transplant", `TX_DAY`, `T_POSTOP`), record the alias here. Distinct from `TIME` (rxode2 time clock) and from `OCC` (integer-valued occasion / period indicator for IOV). When a paper uses `POD` jointly with an IOV occasion column, both can coexist in the dataset: `POD` enters the typical-value covariate equation, `OCC` multiplexes the IOV etas. The 180-day cap in Bergmann 2014 is data-driven (most observations are within the first 90 days post-transplant, so the linear effect is identifiable only over that window) -- document any per-model cap in `covariateData[[POD]]$notes`. Ratified canonically on 2026-05-08 alongside the Bergmann 2014 extraction.

### TTD (**canonical for time to death**)
- **Description:** Days remaining from the observation time to the patient's recorded time of death (TTD >= 0; falls to 0 on the day of death). Time-varying per subject. Available only retrospectively in observational palliative-care / end-of-life cohorts; not usable as a prospective predictor.
- **Units:** days
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used in a first-order exponential decay form on a structural PK parameter. The Franken 2015 morphine model parameterises `CL(TTD) = CL_pop - theta_D * exp(-theta_rate * TTD)`, so the decay term vanishes as TTD -> infinity (asymptotic CL far from death) and reaches its peak drop `theta_D` as TTD -> 0 (day of death).
- **Source aliases:**
  - `TTD` -- used in `Franken_2015_morphine.R` (Franken 2015 NONMEM column for time-to-death in days; paper Eq. 3 and Table 2).
- **Example models:** `Franken_2015_morphine.R` (Franken 2015 Clin Pharmacokinet; first-order exponential decay term on morphine clearance with theta_D = 17.6 L/h and theta_rate = 0.13 /day; CL drops from 47.5 to 29.9 L/h as TTD goes from infinity to 0).
- **Notes:** Specific scope because the covariate's semantics depend on a palliative-care / observational-cohort study design where time of death is known by retrospective abstraction. The Franken 2015 model interprets the TTD-dependent CL change as a composite of physiological end-of-life processes (e.g., reduced hepatic blood flow, cachexia) that are not captured by standard blood-chemistry covariates. For prospective simulation of a virtual cohort without a known time of death, set TTD to a large value (e.g., > 50 days) to recover the asymptotic CL_pop. Distinct from the canonical `POD` (post-operative day, monotonically increasing from a surgery date) -- TTD counts down to death rather than up from a procedure. Ratified canonically on 2026-05-16 alongside the Franken 2015 morphine extraction.

### POSTTX_DAY1 (**canonical for first-24-hours-post-transplant indicator**)
- **Description:** Binary indicator for the first 24 hours post-transplant. 1 = the observation falls within the first 24 hours (day 1) post-transplant; 0 = otherwise. Time-varying per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (any observation outside the first 24 hours post-transplant).
- **Source aliases:**
  - `first day post-transplant` -- used in `Storset_2014_tacrolimus.R`.
- **Example models:** `Storset_2014_tacrolimus.R` (multiplicative ~2.68-fold increase in oral bioavailability on day 1: `fdepot *= 2.68^POSTTX_DAY1`; Storset 2014 Table 2 final theory-based model retains the day-1 factor with subject-level IIV of 57% CV on the day-1 multiplier).
- **Notes:** Distinct from the continuous `POD` (post-operative day) canonical above -- `POSTTX_DAY1` is a derived binary indicator (operationally `POSTTX_DAY1 = as.integer(POD < 1)` when both columns are present in the dataset). The two coexist in the same model when the source paper uses POD-based continuous effects on some parameters and a separate binary day-1 effect on others (Storset 2014 retains the binary day-1 factor explicitly because the day-1 oral bioavailability is ~2.68-fold higher than the rest of the post-transplant period and is not well captured by a continuous POD effect). Storset 2014 Discussion attributes the day-1 oral-bioavailability spike to candidate mechanisms including methylprednisolone-bolus inhibition of intestinal CYP3A / P-glycoprotein, surgery-related inflammation, anaesthesia / opioid effects on gut motility, and reduced food intake -- but no single mechanism was identifiable in the data. The 2.68-fold factor was retained because it produced a 209-point OFV decrease and was crucial for predicting concentrations measured on the first post-transplant day. In Storset 2014 the day-1 effect carries its own subject-level eta (BSV 57% CV on the day-1 factor); only subjects with day-1 observations contribute to that eta. Ratified canonically on 2026-05-08.

### PRIOR_GAST (**canonical for prior gastrectomy**)
- **Description:** Prior (partial or total) gastrectomy indicator, 1 = prior gastrectomy, 0 = no prior gastrectomy. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no prior gastrectomy).
- **Source aliases:**
  - `GAST` -- used in `Yamada_2025_zolbetuximab.R`.
- **Example models:** `Yamada_2025_zolbetuximab.R` (fractional effects on CLss, CLT, V1).
- **Notes:** Renamed from `GAST` on 2026-04-20 to follow the `PRIOR_TNF` / `PRICORT` naming pattern for prior-treatment and surgical-history indicators. Applicable to any PK model where gastrointestinal anatomy affects absorption, first-pass, or protein turnover; not inherently oncology-specific. No distinction between partial vs total gastrectomy unless the source paper separates them.

## Disease state (cross-population indicators)

### DIS_UC (**canonical for ulcerative colitis disease-state indicator**)
- **Description:** 1 = ulcerative colitis patient, 0 = non-UC (e.g., healthy volunteer or non-IBD indication). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-UC subject; the complement group is defined per-model -- typically healthy volunteers and/or patients with another indication such as asthma).
- **Source aliases:**
  - `UC` -- used in `Hua_2015_anrukinzumab.R`.
- **Example models:** `Hua_2015_anrukinzumab.R` (multiplicative fractional increase in CL, +72.8%, on top of weight and albumin effects).
- **Notes:** Used when a population PK model pools UC patients with a non-UC reference population (e.g., Hua 2015: healthy volunteers + asthma patients + UC patients) and UC disease status is tested as a PK covariate. Distinct from `DISEXT_EP` / `DISEXT_OTHER`, which operate *within* a UC-only cohort (disease extension). Start as scope: specific; promote to general if a second paper pools UC with a non-UC reference.

### DIS_SASTHMA (**canonical for moderate-to-severe asthma disease-state indicator**)
- **Description:** 1 = moderate-to-severe asthma patient, 0 = not (e.g., healthy volunteer, mild-to-moderate asthma, or other indication). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-moderate-to-severe-asthma subject; the complement group is defined per-model).
- **Source aliases:**
  - `sAsthma` -- used in `Hua_2015_anrukinzumab.R`.
- **Example models:** `Hua_2015_anrukinzumab.R` (multiplicative fractional change in SC bioavailability, -30.9%).
- **Notes:** The moderate-to-severe vs. mild-to-moderate asthma cutoff is protocol-defined; Hua 2015 uses FEV1 55-80% and ACQ-5 >= 2 for "moderate to severe" (study 4) versus FEV1 > 70% and ACQ-5 <= 1 for "mild to moderate" (study 1). Scope: specific because the severity threshold is tied to a particular analysis plan; future asthma-severity indicators with different thresholds should register as separate canonicals.

### DIS_PJIA (**canonical for polyarticular juvenile idiopathic arthritis disease-state indicator**)
- **Description:** 1 = polyarticular juvenile idiopathic arthritis (pJIA) patient, 0 = non-pJIA (e.g., adult rheumatoid arthritis or other indication). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-pJIA subject; the complement group is defined per-model -- typically adult RA in pooled abatacept analyses).
- **Source aliases:**
  - `JIA` -- used in `Gandhi_2021_abatacept.R` and `Zhong_2026_abatacept.R`.
- **Example models:** `Gandhi_2021_abatacept.R` (additive coefficient on logit-F: pJIA patients have markedly higher SC bioavailability than RA reference), `Zhong_2026_abatacept.R` (additive coefficient +3.08 on logit-F transferred verbatim from a previous internal JIA PPK model that matches Gandhi 2021's published value).
- **Notes:** Used when a population PK model pools pJIA patients with a non-pJIA reference population (e.g., Gandhi 2021: pooled adult RA + pediatric pJIA; Zhong 2026: pooled adult RA + pediatric pJIA + adult/pediatric HM) and pJIA disease/age status is tested as a PK covariate (typically on bioavailability rather than CL). Distinct from `CHILD` and `ADOLESCENT`, which are pure age-band indicators independent of indication. Scope: specific; promote to general if a third paper pools pJIA with a non-pJIA reference and the reference category remains adult RA.

### DIS_CANCER (**canonical for advanced-solid-tumor / oncology cohort indicator**)
- **Description:** 1 = patient with an advanced or metastatic solid tumor (the oncology cohort in a pooled multi-indication PK/PD analysis), 0 = non-oncology subject (healthy volunteer or non-oncology disease cohort pooled in the source analysis). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-oncology subject; the complement group is paper-defined -- typically the union of healthy volunteers and a non-oncology disease cohort such as cGVHD pooled in the source analysis).
- **Source aliases:** none; source NONMEM / Monolix control streams typically derive the indicator from a `POP` or `STUDY` categorical alongside `DIS_HEALTHY`.
- **Example models:** `Yang_2024_axatilimab.R` (multiplicative effect on baseline NCMC: `BL_NCMC x exp(1.22 x DIS_CANCER + 0.618 x DIS_HEALTHY)`; reference category cGVHD when both indicators are 0).
- **Notes:** Used together with `DIS_HEALTHY` to decompose a three-level "participant population" categorical (cGVHD reference, advanced solid tumor, healthy volunteer) into two orthogonal binary indicators. Scope: specific because the disease-pooling reference category is paper-defined (Yang 2024 reference is patients with cGVHD). Ratified canonically on 2026-04-28.

### DIS_CANCER_PED (**canonical for pediatric oncology cohort indicator**)
- **Description:** 1 = pediatric patient receiving cancer-directed therapy (any malignancy, including hematologic cancers such as leukemia and lymphoma as well as solid tumors / blastomas), 0 = pediatric patient admitted for a non-oncology indication (e.g. infection, surgery, transplant). Time-fixed per subject. Distinct from `DIS_CANCER`, which is restricted to advanced/metastatic solid tumors in adults; `DIS_CANCER_PED` is the pediatric variant in the `DIS_CANCER*` family and explicitly covers leukemia-dominant pediatric cohorts.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-oncology pediatric patient; the complement group is paper-defined -- in Llanos-Paez 2020 the complement is pediatric patients admitted for various non-oncology indications, with appendicitis and kidney disease / urinary tract infection the most common).
- **Source aliases:**
  - `ONCOLOGY` -- Llanos-Paez 2020 NONMEM column with the same orientation (1 = oncology, 0 = nononcology); maps directly to `DIS_CANCER_PED`.
- **Example models:** `Llanos-Paez_2020_gentamicin.R` (multiplicative cohort shifts on V1 (-0.154) and Q (-0.321) relative to the nononcology baseline; CL has no oncology effect).
- **Notes:** Use `DIS_CANCER_PED` rather than `DIS_CANCER` whenever the source paper's "oncology" cohort includes hematologic malignancies (leukemia / lymphoma) or pediatric blastomas, because `DIS_CANCER` is canonically restricted to advanced/metastatic solid tumors. Reference-category complement is paper-defined (Llanos-Paez 2020 complement is the pooled pediatric non-oncology admissions cohort). Scope: specific because the complement is paper-defined. Covariate-effect parameters drop the `DIS_` prefix per the existing `DIS_CANCER` -> `e_cancer_*` convention (Yang 2024); use `e_cancer_ped_<param>`.

### DIS_HEALTHY (**canonical for healthy-participant cohort indicator**)
- **Description:** 1 = healthy participant (no diagnosis), 0 = patient (any diagnosis represented in the pooled cohort). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (patient subject; the complement group is the union of disease cohorts pooled in the source analysis).
- **Source aliases:** none known; healthy-participant indicators in source NONMEM control streams typically use ad-hoc names (e.g., `HV`, `HEALTHY`, `DIS_HV`).
- **Example models:** `Nikanjam_2019_siltuximab.R` (multiplicative effects: 0.77 on CL, 0.83 on Vss; reference category is the pooled non-healthy oncology cohort), `Okada_2025_rocatinlimab.R` (multiplicative shift `1 - 0.532` on Vmax when 1; reference complement is the pooled atopic-dermatitis + ulcerative-colitis + plaque-psoriasis patient cohort), `Yang_2024_axatilimab.R` (multiplicative effect on baseline NCMC: `BL_NCMC x exp(1.22 x DIS_CANCER + 0.618 x DIS_HEALTHY)`; reference category cGVHD), `Goel_2016_Sonidegib.R` (multiplicative power-form effect on CL/F: `2.96^DIS_HEALTHY`; reference category is the pooled cancer-patient cohort across X2101 / X1101 / A2201), `Brown_2017_osimertinib.R` (linear factor `(1 + 0.44 x DIS_HEALTHY)` on apparent osimertinib clearance and `(1 + 1.25 x DIS_HEALTHY)` on apparent AZ5104 clearance; reference category is the pooled NSCLC cohort across AURA / AURA2), `Lu_2015_vismodegib.R` (additive-on-log-scale shift on ka via `exp(0.671 * DIS_HEALTHY)` and on F via `exp(0.881 * DIS_HEALTHY)` gated on the Phase I formulation indicator; reference category is the pooled cancer-patient cohort across SHH3925g / SHH4610g / SHH4476g), `Gupta_2016_lenvatinib.R` (multiplicative power-form effect on CL/F: `1.15^DIS_HEALTHY`; reference category is the pooled solid-tumor / thyroid-cancer patient cohort across 15 phase 1-3 studies; healthy subjects show +15 percent CL/F vs cancer patients).
- **Notes:** Used when a population PK model pools healthy participants with patients across heterogeneous indications and the healthy-vs-patient contrast is retained as a covariate. Scope: specific because the complement reference category is paper-defined (Nikanjam 2019 reference is "all non-healthy, non-Castleman, non-SMM tumor types"; Okada 2025 reference is the pooled AD+UC+psoriasis patient cohort; Yang 2024 reference is patients with cGVHD; Goel 2016 reference is the pooled cancer-patient cohort with advanced solid tumors or BCC; Brown 2017 reference is the pooled advanced NSCLC cohort; Lu 2015 reference is the pooled cancer-patient cohort with advanced solid tumors / metastatic or locally-advanced BCC). The retired canonical name `DIS_HV` (healthy-volunteer) was renamed on 2026-05-11 because "volunteer" terminology is discouraged for clinical-trial participants. Ratified canonically on 2026-04-24.

### DIS_CASTLEMAN (**canonical for Castleman's disease indicator**)
- **Description:** 1 = Castleman's disease (multicentric or unicentric), 0 = not Castleman's disease. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-Castleman subject; the complement group is the union of healthy volunteers and other indications pooled in the source analysis).
- **Source aliases:** none known; source NONMEM control streams typically use ad-hoc names (e.g., `CD`, `CASTLEMAN`).
- **Example models:** `Nikanjam_2019_siltuximab.R` (multiplicative +24% effect on CL; no Vss effect).
- **Notes:** Castleman's disease is a lymphoproliferative disorder strongly associated with elevated IL-6 levels; it is the only FDA-approved indication for siltuximab. Scope: specific because the disease-pooling reference category is paper-defined. Ratified canonically on 2026-04-24.

### DIS_HOFH (**canonical for homozygous familial hypercholesterolemia patient indicator**)
- **Description:** 1 = patient with homozygous familial hypercholesterolemia (HoFH), 0 = non-HoFH subject (typically healthy volunteer or another reference cohort pooled in the source analysis). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-HoFH subject; the complement group is paper-defined -- for Pu 2021 the reference is the pooled healthy-volunteer cohort).
- **Source aliases:**
  - `DISTYPN` -- used in `Pu_2021_evinacumab.R` (Pu 2021 NM-TRAN $INPUT column for HoFH-vs-HV disease type, 1 = HoFH).
- **Example models:** `Pu_2021_evinacumab.R` (multiplicative `exp(theta * DIS_HOFH)` factor on Vmax with theta = -0.289, i.e. HoFH patients show ~25% lower target-mediated Vmax than the HV reference; biologically consistent with the LDLR-pathway disruption in HoFH altering ANGPTL3 catabolic kinetics).
- **Notes:** Used when a population PK model pools HoFH patients with healthy volunteers (or another non-HoFH cohort) and HoFH disease status is retained as a covariate. Distinct from a heterozygous-FH (HeFH) indicator because HoFH patients have markedly higher baseline LDL-C (untreated levels often > 500 mg/dL) and a more pronounced response to LDLR-independent therapies. Scope: specific because the reference category is paper-defined.

### DIS_HAE (**canonical for hereditary angioedema patient indicator**)
- **Description:** 1 = patient with hereditary angioedema (HAE-C1INH-Type1, HAE-C1INH-Type2, or HAE-nC1INH), 0 = healthy volunteer (or other non-HAE reference cohort pooled in the source analysis). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-HAE subject; the complement group is paper-defined -- for Diep 2026 the reference is the pooled healthy-volunteer cohort from NCT03263507 and ISIS 721744-CS9).
- **Source aliases:** paper narrative "patient with HAE" / "healthy volunteer" subgroup labels driving the Diep 2026 disease-status covariate effects on Vc/F, Q/F, baseline PKK, and IC50.
- **Example models:** `Diep_2026_donidalorsen.R` (linear `(1 + theta * DIS_HAE)` multiplicative effects on apparent central volume Vc/F (theta = +0.426, +42.6%), apparent intercompartmental clearance Q/F (theta = -0.261, -26.1%), baseline plasma prekallikrein BL (theta = -0.132, -13.2%), and donidalorsen IC50 on PKK production (theta = +0.770, +77.0%) for patients with HAE vs healthy volunteers).
- **Notes:** Used when a population PK/PD model pools HAE patients with healthy volunteers and HAE disease status is retained as a covariate. The three molecular HAE subtypes (HAE-C1INH-Type1, HAE-C1INH-Type2, HAE-nC1INH) are pooled in this indicator following the Diep 2026 analysis; if a future paper resolves subtype-specific covariate effects, separate canonical indicators (e.g., `DIS_HAE_C1INH_T1`) can be added without conflicting with this pooled indicator. Scope: specific because the complement reference category is paper-defined.

### DIS_DMD (**canonical for Duchenne muscular dystrophy patient indicator**)
- **Description:** 1 = patient with Duchenne muscular dystrophy (DMD), 0 = non-DMD subject (healthy volunteer or other reference cohort). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-DMD subject; the complement group is the reference cohort the source analysis pools alongside the DMD population -- typically healthy adult volunteers).
- **Source aliases:** `SPOP` (Wojciechowski 2022 study-population indicator with the same orientation: 1 = DMD pediatric patient, 0 = healthy adult volunteer).
- **Example models:** `Wojciechowski_2022_domagrozumab.R` (additive `1 + theta` shift on baseline myostatin and on the joint kdeg/kint axis; theta_BASE = -0.641, theta_kdegkint = -0.900).
- **Notes:** Used when a population PK/PD model pools DMD patients with a non-DMD reference population and DMD disease status is retained as a covariate. Scope: specific because the reference category is paper-defined. Ratified canonically on 2026-04-26.

### DIS_SMM (**canonical for smoldering multiple myeloma indicator**)
- **Description:** 1 = smoldering (asymptomatic) multiple myeloma, 0 = not smoldering MM. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-SMM subject; the complement group is the union of healthy volunteers and other indications pooled in the source analysis).
- **Source aliases:** none known; source NONMEM control streams typically use ad-hoc names (e.g., `SMM`, `SMOLDMM`).
- **Example models:** `Nikanjam_2019_siltuximab.R` (multiplicative -23% effect on Vss; no CL effect).
- **Notes:** Smoldering multiple myeloma is an asymptomatic plasma-cell disorder distinct from active multiple myeloma; pooled with the Nikanjam 2019 cohort that also included MGUS, multiple myeloma, RCC, ovarian, and other tumor types. Scope: specific because the disease-pooling reference category is paper-defined. Ratified canonically on 2026-04-24.

### MM (**canonical for active multiple myeloma disease indicator**)
- **Description:** 1 = active (non-smoldering) multiple myeloma, 0 = other hematologic malignancy or reference group. Time-fixed per subject.

### DIS_PNH (**canonical for paroxysmal nocturnal hemoglobinuria indicator**)
- **Description:** 1 = paroxysmal nocturnal hemoglobinuria (PNH) patient, 0 = non-PNH subject (healthy volunteer or another indication pooled in the source analysis). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-PNH subject; the complement group is paper-defined -- for Lin 2024 it pools healthy adult volunteers and CHAPLE disease patients).
- **Source aliases:** none known; source NONMEM control streams typically use ad-hoc names (e.g., `PNH`, `DPNH`).
- **Example models:** `Lin_2024_pozelimab.R` (additive-fractional +34.07% effect on Vc; no CL or Vp effect; reference category pools healthy volunteers and CHAPLE patients).
- **Notes:** Paroxysmal nocturnal hemoglobinuria is a rare hematological disease characterized by uncontrolled complement activation on red blood cells; treated with C5-targeted complement inhibitors (eculizumab, ravulizumab, pozelimab). Scope: specific because the disease-pooling reference category is paper-defined. Ratified canonically on 2026-04-27.

### MDSAML (**canonical for MDS or AML disease-type indicator**)
- **Description:** 1 = patient with myelodysplastic syndrome (MDS) or acute myeloid leukemia (AML), 0 = other hematologic malignancy or reference group. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-MDS/AML subjects; the complement group is defined per-model -- typically multiple myeloma and non-Hodgkin lymphoma in Ogasawara 2020).
- **Source aliases:** none; `MDSAML` is the combined indicator used directly in source analyses.
- **Example models:** `Ogasawara_2020_durvalumab.R` (multiplicative factor 1.26 on CL; reference group is the union of MM and NHL subjects).
- **Notes:** Use `MDSAML` as a combined MDS+AML indicator when the source paper collapses the two diagnoses into one covariate. If a future paper separates MDS and AML as distinct indicators, register `DIS_MDS` and `DIS_AML` separately. Scope: specific because the reference category is paper-defined. Ratified canonically on 2026-04-26.

### DIS_AML (**canonical for acute myeloid leukemia disease-state indicator**)
- **Description:** 1 = patient with acute myeloid leukemia (AML), 0 = non-AML subject (the complement group in a pooled multi-indication PK analysis). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-AML subject; the complement group is paper-defined -- for Xu 2023 the reference is patients with advanced solid tumors, alongside the parallel `DIS_MDS` and `DIS_CMML` indicators that decompose the hematologic-malignancy cohort).
- **Source aliases:**
  - `DISEASE_abb == "AML"` -- used in `Xu_2023_MBG453.R` (the Monolix supplement Appendix S2 encodes disease as the categorical column `DISEASE_abb` with categories `{AML, CMML, MDS, Solid_Tumor}` and reference `Solid_Tumor`; the canonical column carries the binary `as.integer(DISEASE_abb == "AML")`).
- **Example models:** `Xu_2023_MBG453.R` (exponential effect on CL: `exp(-0.0146 * DIS_AML)`; not statistically significant in the full covariate model but retained because Xu 2023 used the full-covariate-model approach).
- **Notes:** Use `DIS_AML` (rather than the combined `MDSAML`) when the source paper separates AML from MDS as distinct indicators. Scope: specific because the disease-pooling reference category is paper-defined.

### DIS_MDS (**canonical for myelodysplastic syndrome disease-state indicator**)
- **Description:** 1 = patient with myelodysplastic syndrome (MDS), 0 = non-MDS subject (the complement group in a pooled multi-indication PK analysis). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-MDS subject; the complement group is paper-defined -- for Xu 2023 the reference is patients with advanced solid tumors, alongside the parallel `DIS_AML` and `DIS_CMML` indicators).
- **Source aliases:**
  - `DISEASE_abb == "MDS"` -- used in `Xu_2023_MBG453.R` (Monolix supplement Appendix S2 categorical column; reference category `Solid_Tumor`).
- **Example models:** `Xu_2023_MBG453.R` (exponential effect on CL: `exp(-0.149 * DIS_MDS)`; statistically significant, p = 0.021 -- patients with MDS have ~14% lower CL than the solid-tumor reference).
- **Notes:** Use `DIS_MDS` (rather than the combined `MDSAML`) when the source paper separates MDS from AML as distinct indicators. Scope: specific because the disease-pooling reference category is paper-defined.

### DIS_CMML (**canonical for chronic myelomonocytic leukemia disease-state indicator**)
- **Description:** 1 = patient with chronic myelomonocytic leukemia (CMML), 0 = non-CMML subject (the complement group in a pooled multi-indication PK analysis). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-CMML subject; the complement group is paper-defined -- for Xu 2023 the reference is patients with advanced solid tumors, alongside the parallel `DIS_AML` and `DIS_MDS` indicators that decompose the hematologic-malignancy cohort).
- **Source aliases:**
  - `DISEASE_abb == "CMML"` -- used in `Xu_2023_MBG453.R` (Monolix supplement Appendix S2 categorical column; reference category `Solid_Tumor`).
- **Example models:** `Xu_2023_MBG453.R` (exponential effect on CL: `exp(-0.0411 * DIS_CMML)`; not statistically significant in the full covariate model but retained because Xu 2023 used the full-covariate-model approach).
- **Notes:** CMML is a clonal myeloid malignancy with overlapping features of MDS and myeloproliferative neoplasms. Scope: specific because the disease-pooling reference category is paper-defined.

### DIS_BCPALL (**canonical for B-cell precursor acute lymphoblastic leukemia disease-state indicator**)
- **Description:** 1 = B-cell precursor acute lymphoblastic leukemia (BCP-ALL), 0 = B-cell non-Hodgkin's lymphoma (NHL) or other non-BCP-ALL indication pooled in the source analysis. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-BCP-ALL; in the Wu 2024 cohort this is the adult B-cell NHL stratum).
- **Source aliases:**
  - `ALL` -- used in `Wu_2024_inotuzumab.R` (Wu 2024 calls it the "ALL effect" and notes it bundles disease type with the corresponding bioanalytical assay difference).
- **Example models:** `Wu_2024_inotuzumab.R` (additive fractional-change effects on CL1 (-0.767) and CL2 (-0.362), and gates the BLSTABL and AGE effects on kdes; for kdes itself a -0.924 fractional change for BCP-ALL).
- **Notes:** Used when a population PK model pools BCP-ALL patients with a non-BCP-ALL reference (e.g., Wu 2024: pooled adult B-cell NHL + adult BCP-ALL + pediatric BCP-ALL). Scope: specific because the complement reference category is paper-defined (Wu 2024 reference is pooled adult B-cell NHL). The "ALL effect" theta in Wu 2024 conflates two physiologically distinct sources of variation -- B-cell tumor type (NHL vs ALL surface CD22 burden) and bioanalytical method (ELISA for adult NHL vs HPLC-MS for ALL) -- and cannot be split with the available data; document this confounding when comparing across populations. Ratified canonically on 2026-04-26.

### DIS_SAD (**canonical for secondary antibody deficiency indicator**)
- **Description:** 1 = secondary antibody deficiency (SAD) patient (hypogammaglobulinaemia from external causes such as B-cell-depleting therapy, haematological malignancy, or other immunosuppression), 0 = primary immunodeficiency (PID) patient. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (PID patient; the complement category is the genetic / inborn-error-of-immunity primary immunodeficiency cohort pooled with SAD in the source analysis).
- **Source aliases:** none known; source NONMEM control streams typically use ad-hoc names (e.g., `SAD`, `IMD`, `DIS`).
- **Example models:** `Cheng_2026_immunoglobulin.R` (multiplicative `theta^DIS_SAD` factors on CL (0.542) and on baseline IgG (CBAS, 0.541); reference category PID).
- **Notes:** Used when a population PK model pools PID and SAD pediatric or adult patients receiving immunoglobulin replacement therapy (IgRT) and tests SAD-vs-PID as a covariate. Distinct from the disease-state indicators that pool oncology / autoimmune indications: `DIS_SAD` specifically partitions hypogammaglobulinaemia by its underlying mechanism (genetic vs. acquired). Scope: specific because the SAD cohort composition is paper-defined (in Cheng 2026, 75% post-rituximab and 25% post-CAR-T cell therapy). Ratified canonically on 2026-04-28.

### DIS_AD (**canonical for Alzheimer's disease patient indicator**)
- **Description:** 1 = participant with Alzheimer's disease (clinical AD diagnosis), 0 = non-AD subject (typically healthy volunteer pooled in the source analysis). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-AD subject; the complement group is paper-defined -- for Perez-Ruixo 2025 the reference is the pooled healthy-volunteer cohort).
- **Source aliases:** none known; source NONMEM control streams typically use ad-hoc names (e.g., `AD`, `STATUS`, `DISGRP`).
- **Example models:** `PerezRuixo_2025_posdinemab.R` (acts on baseline free p217+tau in CSF, R0; healthy R0 = 0.793 pmol/L vs AD R0 = 5.995 pmol/L, a 656% relative increase, no PK-parameter effects).
- **Notes:** Used when a population PK/PD model pools healthy volunteers with Alzheimer's disease patients and the AD-vs-HV contrast is retained as a covariate on a target-related parameter (e.g., baseline p-tau, baseline p217+tau). Scope: specific because the complement reference category is paper-defined. Ratified canonically on 2026-04-28.

### DIS_COPD (**canonical for chronic obstructive pulmonary disease patient indicator**)
- **Description:** 1 = patient with chronic obstructive pulmonary disease (clinical COPD diagnosis, typically moderate-to-severe per GOLD criteria), 0 = non-COPD subject (typically healthy volunteer pooled in the source analysis). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-COPD subject; the complement group is paper-defined -- for Lahu 2010 the reference is the pooled phase I healthy-volunteer cohort).
- **Source aliases:**
  - `COPD` -- used in `Lahu_2010_roflumilast.R` (paper text covariate symbol in equation 6 and 7).
- **Example models:** `Lahu_2010_roflumilast.R` (linear additive effects on roflumilast parent CL (-39.4%) and V1 (+184%) and on roflumilast N-oxide CL (-7.9%) and Vd (-21.4%); reference category 0 = pooled phase I healthy volunteers, 1 = pooled phase II/III moderate-to-severe COPD patient).
- **Notes:** Used when a population PK/PD model pools healthy volunteers with COPD patients and the COPD-vs-HV contrast is retained as a covariate on PK parameters. Scope: specific because the complement reference category and the COPD-severity inclusion criteria are paper-defined.

### DIS_OBESE_MORBID (**canonical for morbidly obese cohort indicator**)
- **Description:** 1 = morbidly obese patient (BMI > 40 kg/m^2 in the canonical definition; typical pooled-analysis enrollment criterion is bariatric-surgery patients), 0 = non-obese subject (typically healthy volunteer pooled in the source analysis). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-obese subject; the complement group is paper-defined -- for de Hoogd 2017 the reference is the pooled healthy-volunteer cohort from Sarton 2000 and Romberg 2004).
- **Source aliases:** none known; source NONMEM control streams typically use ad-hoc names (e.g., `OBESE`, `MO`, `COHORT`).
- **Example models:** `deHoogd_2017_morphine.R` (selects per-cohort proportional residual error magnitudes for each of three observed species -- morphine, M3G, M6G -- after a pooled-cohort fit of 20 morbidly obese surgical patients and 20 healthy volunteers).
- **Notes:** Used when a population PK or PK/PD model pools morbidly obese patients with a non-obese reference population (typically healthy volunteers) and the cohort indicator selects per-cohort parameter values (residual error magnitudes, study-specific bioavailability, or similar). Distinct from `BMI` (which is the continuous body-mass-index covariate used for parameter scaling) -- `DIS_OBESE_MORBID` is the binary cohort-membership flag and does not encode a specific BMI threshold for general use; the threshold is paper-defined. Scope: specific because the complement reference category is paper-defined. Ratified canonically on 2026-05-11.

### HSCT_URD_7OF8 (**canonical for hematopoietic stem cell transplant from a 7-of-8 HLA-matched unrelated donor**)
- **Description:** 1 = patient received an allogeneic hematopoietic stem cell transplant (HSCT) from an unrelated donor (URD) HLA-matched at 7 of 8 alleles (single-allele mismatch), 0 = otherwise (the union of patients not in this transplant cohort, including non-HSCT patients and HSCT recipients matched at all 8 alleles). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (not in the 7-of-8-matched HSCT cohort; the complement group is paper-defined -- for Zhong 2026 it pools the RA + pJIA studies and the 8-of-8-matched HSCT cohort, with the latter encoded by the parallel `HSCT_URD_8OF8` indicator).
- **Source aliases:**
  - `COHORT7` -- used in `Zhong_2026_abatacept.R` (Zhong 2026 NM-TRAN indicator for ABA2 study Cohort 7/8).
- **Example models:** `Zhong_2026_abatacept.R` (exponential coefficient -0.326 on CL; the single-allele-mismatch HSCT cohort exhibits ~28% lower abatacept clearance than the reference complement).
- **Notes:** Used together with `HSCT_URD_8OF8` to decompose a three-level "transplant cohort" categorical (non-HSCT-cohort / 7-of-8 / 8-of-8) into two orthogonal binary indicators. The 7-of-8 cohort represents a higher GvHD-risk population because of the single-allele HLA mismatch. Scope: specific because the reference complement (the union of non-transplant disease cohorts pooled in the source analysis) is paper-defined. Ratified canonically on 2026-04-29.

### HSCT_URD_8OF8 (**canonical for hematopoietic stem cell transplant from an 8-of-8 HLA-matched unrelated donor**)
- **Description:** 1 = patient received an allogeneic hematopoietic stem cell transplant (HSCT) from an unrelated donor (URD) HLA-matched at all 8 alleles (full match), 0 = otherwise (the union of patients not in this transplant cohort, including non-HSCT patients and HSCT recipients matched at 7 of 8 alleles). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (not in the 8-of-8-matched HSCT cohort; the complement group is paper-defined -- for Zhong 2026 it pools the RA + pJIA studies and the 7-of-8-matched HSCT cohort, with the latter encoded by the parallel `HSCT_URD_7OF8` indicator).
- **Source aliases:**
  - `COHORT8` -- used in `Zhong_2026_abatacept.R` (Zhong 2026 NM-TRAN indicator for ABA2 study Cohort 8/8).
- **Example models:** `Zhong_2026_abatacept.R` (exponential coefficient -0.0934 on CL and +0.257 on VC; the fully-HLA-matched HSCT cohort exhibits a small CL decrease and a larger VC increase relative to the reference complement).
- **Notes:** Used together with `HSCT_URD_7OF8` to decompose a three-level "transplant cohort" categorical (non-HSCT-cohort / 7-of-8 / 8-of-8) into two orthogonal binary indicators. The 8-of-8 cohort is the lower-risk HLA-matching configuration. Scope: specific because the reference complement (the union of non-transplant disease cohorts pooled in the source analysis) is paper-defined. Ratified canonically on 2026-04-29.

### DIS_PSORIASIS (**canonical for plaque psoriasis disease-state indicator**)
- **Description:** 1 = plaque psoriasis patient, 0 = non-psoriasis subject (e.g., atopic dermatitis, ulcerative colitis, or healthy volunteer). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-psoriasis subject; the complement group is paper-defined -- the union of other disease cohorts pooled in the source analysis).
- **Source aliases:** none known; source NONMEM control streams typically use a categorical `DIS` indicator (e.g., Okada 2025: `DIS=1` for psoriasis, `DIS=0` for healthy, `DIS=2` for UC, `DIS=3` for AD), decomposed into a binary `DIS_PSORIASIS` indicator at ingestion.
- **Example models:** `Okada_2025_rocatinlimab.R` (multiplicative shift `1 - 0.372` on linear CL when 1; reference complement is the pooled atopic dermatitis + ulcerative colitis + healthy-volunteer cohort).
- **Notes:** Used when a population PK model pools plaque-psoriasis patients with a non-psoriasis reference population and psoriasis disease status is retained as a covariate. Scope: specific because the disease-pooling reference category is paper-defined. Ratified canonically on 2026-04-27.

### CARRAGEENAN (**canonical for intraplantar-carrageenan inflammatory-challenge indicator**)
- **Description:** Binary indicator for intraplantar injection of carrageenan suspension as an experimental inflammatory / hyperalgesic challenge. 1 = subject received an intraplantar carrageenan injection at the start of the experiment (the carrageenan-induced peripheral inflammation / thermal-hyperalgesia paradigm); 0 = subject received an intraplantar saline injection (sham control). Time-fixed per subject within an experiment.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (saline-injected sham animal; no induced inflammation).
- **Source aliases:**
  - `CARRAGEENAN` -- used in `VasquezBahena_2009_lumiracoxib_rat.R` (1 = groups II-IX, 100 uL of 1% carrageenan suspension into the right hind paw; 0 = group I, 100 uL of 0.9% saline solution).
- **Example models:** `VasquezBahena_2009_lumiracoxib_rat.R` (switches the COX-2 synthesis-rate model: CARRAGEENAN = 0 selects the constant saline synthesis rate `ks_cox2_saline` and CARRAGEENAN = 1 selects the time-variant gamma function `ks_cox2(t) = A * t^alpha * exp(-beta * t)` driving the carrageenan-induced inflammation profile).
- **Notes:** The carrageenan-induced peripheral inflammation / hyperalgesia model (Winter 1962; Hargreaves 1988 thermal-hyperalgesia variant) is one of the most widely used preclinical assays for screening anti-inflammatory and analgesic drugs in rodents, so the canonical name is reusable for future preclinical extractions. Scope: specific until a second model ratifies the binary-switch semantics on a different PD framework. The Hargreaves-test thermal-hyperalgesia readout (`LT` paw withdrawal latency, seconds) is the typical observable when this indicator is in use, but other readouts (paw oedema, mechanical-allodynia von Frey threshold) are equally valid. Distinct from disease-state indicators (`DIS_*`) because the inflammatory state is experimentally induced at a defined time, not a chronic patient condition.

## Pulmonary / lung-disease biomarkers

### FEV1 (**canonical for forced expiratory volume in 1 second**)
- **Description:** Baseline forced expiratory volume in 1 second (FEV1) reported as an absolute volume in litres. Pulmonary-function spirometry endpoint; reflects large-airway airflow obstruction and is a standard covariate in chronic obstructive pulmonary disease, alpha-1 antitrypsin deficiency, asthma, and cystic-fibrosis disease-progression analyses. Distinct from FEV1 percent-predicted (the % predicted value standardises absolute FEV1 by reference equations from sex / age / height / ethnicity); use this canonical only for the absolute-litre value supplied as a covariate column.
- **Units:** L (absolute volume; document the exhalation-effort standard the source paper cites in `covariateData[[FEV1]]$notes` if non-default -- ATS / ERS post-bronchodilator is the typical convention).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with linear-deviation forms `(theta + e_fev1_param * (FEV1 - ref))` or power forms `(FEV1 / ref)^exponent`. Reference values observed: 1.6 L (Tortorici 2017, population median across the RAPID-RCT/RAPID-OLE A1-PI augmentation cohort).
- **Source aliases:** none yet; canonical name preferred. Source papers typically use the abbreviation `FEV1` directly.
- **Example models:** `Tortorici_2017_a1pi.R` (linear-deviation effect on the lung-density decline rate: `theta5 * (FEV1 - 1.6)` with `theta5 = +0.56 (g/L/year per L FEV1)`; lower-FEV1 patients have steeper natural decline rates independent of A1-PI exposure).
- **Notes:** Distinct from FEV1 percent-predicted (which is a derived ratio with the reference-equation denominator built in; FEV1% is the outcome variable in `Harun_2019_cysticFibrosis.R` rather than a covariate column). Use `FEV1` only when the source paper supplies the absolute-volume value; if the source supplies a percent-predicted value as a covariate, the canonical for that surface is `FEV1_PCTPRED`. Scope: specific until a second model ratifies the absolute-litre semantics; promote to general at that point. Ratified canonically on 2026-05-09 alongside the Tortorici 2017 extraction.

### FEV1_PCTPRED (**canonical for baseline FEV1 as percent of the predicted value**)
- **Description:** Baseline forced expiratory volume in 1 second expressed as a percent of the sex / age / height / ethnicity reference-equation predicted value (FEV1% predicted). The percent-predicted scaling normalises the absolute-litre FEV1 measurement against a healthy-reference standard so a single covariate value is comparable across patient ages and body sizes. Standard pulmonary-function covariate in cystic-fibrosis disease-severity and inhaled-therapy popPK analyses. Time-fixed at study entry / baseline.
- **Units:** % predicted (numeric percentage, e.g. 62.1 for the Ting 2014 population median; not a fraction 0.621).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with power-form effects `(FEV1_PCTPRED / ref)^exponent` or linear-deviation forms `(1 + e * (FEV1_PCTPRED - ref))`. Reference values observed: 62.1 % (Ting 2014, population median across the combined three-study cystic-fibrosis cohort).
- **Source aliases:**
  - `FEV1% predicted` -- used in `Ting_2014_tobramycin_inhaled.R` (paper Table 1 / Table 2). Free-text label, not a typical NMTRAN column header; downstream data sets should use the canonical column name `FEV1_PCTPRED`.
- **Example models:** `Ting_2014_tobramycin_inhaled.R` (power-form effect on apparent central volume of distribution: `(FEV1_PCTPRED / 62.1)^-0.303`; lower lung-function patients have larger apparent Vd/F, consistent with the paper's hypothesis that worsening lung disease increases central-airway aerosol deposition).
- **Notes:** Distinct from absolute-litre `FEV1` (canonical for the unstandardised volume) and from `FEV1` as an outcome variable in disease-progression models (e.g. `Harun_2019_cysticFibrosis.R`, where percent-predicted FEV1 is the dependent variable rather than a covariate). The reference equation used to derive the percent-predicted value is paper-specific; document any non-default reference standard (Hankinson 1999 NHANES III, GLI 2012, Wang 1993, etc.) in `covariateData[[FEV1_PCTPRED]]$notes` per model. Scope: specific until a second model ratifies the percent-predicted semantics; promote to general at that point. Ratified canonically on 2026-05-12 alongside the Ting 2014 extraction.

### A1PI (**canonical for serum alpha-1 proteinase inhibitor concentration**)
- **Description:** Baseline serum alpha-1 proteinase inhibitor (A1-PI; also known as alpha-1 antitrypsin, AAT) concentration. Used in alpha-1 antitrypsin deficiency (AATD) augmentation-therapy modelling as a per-subject pre-treatment exposure covariate (the subject's endogenous A1-PI level at study entry, before any augmentation infusions). Time-fixed per subject. Distinct from a time-course of A1-PI used as a state variable -- when the source paper carries A1-PI as the dynamic dependent variable (the augmentation model's PD output), use `Cc` rather than `A1PI`.
- **Units:** umol/L (typical SI-convention reporting; also reported as mg/dL in US-convention papers -- document the unit used in each model via `covariateData[[A1PI]]$units`). Conversion: 1 umol/L A1-PI ~= 5.2 mg/dL (using MW ~52 kDa). The 11 umol/L "putative protective threshold" used clinically corresponds to ~57 mg/dL.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with power-form effects `(A1PI / ref)^exponent` on the post-treatment exposure intercept and slope. Reference values observed: 5.5 umol/L (Tortorici 2017, approximate median pre-treatment A1-PI among RAPID-RCT placebo-randomised patients, used as the normalisation denominator for power-form covariate effects).
- **Source aliases:** `Cbase` -- used in Tortorici 2017's published equation 6 to denote the baseline pre-treatment A1-PI value; the column name in the modelled dataset would be `A1PI`.
- **Example models:** `Tortorici_2017_a1pi.R` (two power-form effects: `(A1PI/5.5)^theta5` with `theta5 = +0.73` on the placebo-arm post-treatment exposure intercept, and `(A1PI/5.5)^theta4` with `theta4 = -0.12` on the dose-rate slope; together they encode the modest dose-exposure dependence on each subject's endogenous A1-PI level).
- **Notes:** AATD enrolment criteria typically restrict A1PI to <= 11 umol/L (severe deficiency); reference / heterozygous PI*MZ phenotypes have higher levels. Document the source paper's AATD-genotype enrolment criteria in `covariateData[[A1PI]]$notes` per model. Specific scope because the canonical is tied to AATD augmentation-therapy modelling; future PK / PD analyses of A1-PI (or AAT) in non-AATD contexts (acute-phase response, smoking-induced inflammation) should ratify general scope at that time. Ratified canonically on 2026-05-09 alongside the Tortorici 2017 extraction.

## Cystic fibrosis lung-disease indicators

### AIR_TRAP_5Y (**canonical for severe air trapping on chest HRCT scan at age 5 years**)
- **Description:** 1 = subject had a non-zero "air trapping" component score on the validated Brody-II chest high-resolution computed tomography (HRCT) scan performed at age 5 years; 0 = air trapping component score of 0 (absent). Time-fixed per subject (the indicator captures the single end-of-study HRCT performed at age 5 in the Australasian Cystic Fibrosis Bronchoalveolar Lavage (ACFBAL) study). The Brody-II scoring system reports air trapping as a percentage of maximum possible HRCT score; the binary "present vs absent" dichotomy follows the same convention as Rosenow 2015 (Am J Respir Crit Care Med 191:1158-65).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no air trapping detected on the age-5 HRCT scan).
- **Source aliases:**
  - `ATS5C` -- used in `Harun_2019_cysticFibrosis.R` (Harun 2019 NMTRAN `$INPUT` column, "presence of air trapping at age 5"; values 0 = absent, 1 = present).
- **Example models:** `Harun_2019_cysticFibrosis.R` (linear-deviation effect on baseline FEV1% predicted at age 5: `(1 + e_at_baseline * AIR_TRAP_5Y)` with coefficient -0.0417, i.e., subjects with severe air trapping at age 5 have a baseline FEV1% predicted approximately 4.17% lower than those without).
- **Notes:** Specific scope because the indicator is tied to the ACFBAL study's standardised HRCT-at-age-5 protocol and the Brody-II scoring rubric; future paediatric CF lung-disease studies that score HRCT at a different age or use a different scoring system should register a separate canonical (`AIR_TRAP_8Y`, `AIR_TRAP_PRAGMA`, etc.) rather than overload this name. Ratified canonically on 2026-05-08 alongside the Harun 2019 extraction.

### HOSPRA (**canonical for hospitalisation due to a pulmonary exacerbation**)
- **Description:** Time-varying binary indicator of inpatient hospitalisation for management of a pulmonary exacerbation at the time of the FEV1% predicted measurement: 1 = subject is hospitalised because of a pulmonary exacerbation when the spirometry value is recorded, 0 = not hospitalised at that visit. Used as a per-visit covariate in disease-progression models of FEV1% decline in cystic fibrosis where pulmonary exacerbations are tracked from the Australian Cystic Fibrosis Data Registry (ACFDR) inpatient records.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (not hospitalised for a pulmonary exacerbation at the time of the FEV1% measurement).
- **Source aliases:**
  - `HOSPRA` -- used in `Harun_2019_cysticFibrosis.R` (Harun 2019 NMTRAN `$INPUT` column, "hospitalisation at the time of FEV1% predicted measurement; 0=no, 1=yes").
- **Example models:** `Harun_2019_cysticFibrosis.R` (linear-deviation effects on the disease-progression maximum drop and on the half-effect age: `(1 + e_hpe_dmax * HOSPRA)` with coefficient -0.22 on the maximum FEV1% drop and `(1 + e_hpe_t50max * HOSPRA)` with coefficient -0.235 on the age at which 50% of the maximum drop occurs; hospitalised visits accelerate both the magnitude and the onset of FEV1% decline).
- **Notes:** Specific scope because the canonical encoding pools all pulmonary exacerbation-driven hospitalisations into a single binary regardless of severity, duration, or treatment intensity; future CF / chronic-respiratory-disease studies that need to distinguish exacerbation severity (e.g., requiring intravenous antibiotics vs oral) should register a finer-grained canonical. Ratified canonically on 2026-05-08 alongside the Harun 2019 extraction.

## Infectious disease

### LNPC (**canonical for log-transformed admission Plasmodium parasitaemia**)
- **Description:** Natural logarithm of the asexual Plasmodium parasite count (parasites per microlitre of blood) at study admission. Time-fixed per subject (one value per subject, captured at enrolment).
- **Units:** log(parasites/uL)
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with linear-deviation forms `(1 + e * (LNPC - ref))`. Reference values observed: 5.88 log(parasites/uL) (Birgersson 2019, population median in the pooled pregnant + non-pregnant Burkina Faso cohort).
- **Source aliases:** none formally; companion column `PARA` (raw asexual parasite count per microlitre) is provided alongside `LNPC` in the Birgersson 2019 NONMEM dataset but the model uses `LNPC = log(PARA)` as the active covariate.
- **Example models:** `Birgersson_2019_artesunate.R` (linear-deviation effect on relative bioavailability `F1`: `F1LNPC = 1 + e_lnpc_f * (LNPC - 5.88)`; positive coefficient `e_lnpc_f = +0.138` per unit increase in log-parasite-count, reflecting increased oral artesunate bioavailability with higher parasite burden).
- **Notes:** Disease-severity covariate specific to malaria PK models. Higher parasitaemia is a marker of more severe acute malaria infection and has been associated in the source publication with altered oral bioavailability of artesunate (presumably via gut-mucosal / first-pass effects of the febrile parasitised state). Scope: specific because the canonical reference value (5.88) is the Birgersson 2019 cohort median; future malaria-in-pregnancy or malaria-in-children PK models may legitimately reuse `LNPC` but should document their own cohort-specific reference value in `covariateData[[LNPC]]$notes`. Ratified canonically on 2026-05-07.

### SARS_VLOAD (**canonical for SARS-CoV-2 baseline viral load**)
- **Description:** Baseline (pre-treatment) SARS-CoV-2 viral load measured from nasopharyngeal swab by RT-qPCR, reported as log10 RNA copies/mL.
- **Units:** log10 copies/mL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with power scaling `(SARS_VLOAD / ref)^exponent`. Reference values observed: 6.4 log10 copies/mL (Lin 2024, median in pooled COVID-19 cohort).
- **Source aliases:**
  - `VIRAL` -- used in `Lin_2024_casirivimab.R`.
- **Example models:** `Lin_2024_casirivimab.R` (small negative exponent -0.0075 on CL).
- **Notes:** SARS-CoV-2-specific. For non-infected subjects, the value is encoded as 0 in the source dataset (below assay detection); the population-PK exponent is small enough that this 0 is absorbed by the reference shift. Register a parallel canonical for any future paper that uses a different infection (e.g., RSV, influenza).

### SARS_SEROPOS (**canonical for SARS-CoV-2 baseline serostatus positive**)
- **Description:** 1 = SARS-CoV-2 spike or nucleocapsid antibody positive at baseline (prior infection or prior vaccination), 0 = seronegative or other / unknown.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (seronegative; "Other" / unknown serostatus is typically pooled into the reference per the source paper's analysis plan).
- **Source aliases:**
  - `SERPOS` -- used in `Lin_2024_casirivimab.R`.
- **Example models:** `Lin_2024_casirivimab.R` (multiplicative fractional change on CL).
- **Notes:** SARS-CoV-2-specific. The exact assay (anti-spike vs anti-nucleocapsid; vendor) varies by study; document per-model in `covariateData[[SARS_SEROPOS]]$notes`.

### OXYSUP_LOW (**canonical for low-flow supplemental oxygen indicator**)
- **Description:** 1 = subject is receiving low-flow supplemental oxygen at baseline (e.g., nasal cannula, simple face mask), 0 = no supplemental oxygen at baseline.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no supplemental oxygen at baseline; the high-flow / mechanical-ventilation categories are encoded by the parallel `OXYSUP_HIGH` indicator).
- **Source aliases:**
  - `OXYSTAT1` -- used in `Lin_2024_casirivimab.R`.
- **Example models:** `Lin_2024_casirivimab.R` (multiplicative fractional change on CL; +10.6%).
- **Notes:** Decomposed indicator from a 4-level ordered categorical (no oxygen / low-flow / high-flow / mechanical ventilation). Use with the parallel `OXYSUP_HIGH` indicator. Register a separate `OXYSUP_VENT` canonical if a future analysis splits mechanical ventilation from high-flow oxygen.

### OXYSUP_HIGH (**canonical for high-flow supplemental oxygen indicator**)
- **Description:** 1 = subject is receiving high-flow supplemental oxygen at baseline (high-flow nasal cannula, non-rebreather mask, non-invasive positive-pressure ventilation, OR mechanical ventilation pooled into the high-flow category), 0 = otherwise (no supplemental oxygen or low-flow).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no supplemental oxygen at baseline; document whether the source paper pooled mechanical ventilation into this indicator or treated it separately).
- **Source aliases:**
  - `OXYSTAT2` -- used in `Lin_2024_casirivimab.R`.
- **Example models:** `Lin_2024_casirivimab.R` (multiplicative fractional change on CL; +38.0%).
- **Notes:** Companion indicator to `OXYSUP_LOW`. In Lin 2024 the rare mechanical-ventilation cases were pooled into the high-flow indicator (n = 24 across the 7598-subject dataset).

### HIV_POS (**canonical for HIV-positive comorbidity indicator**)
- **Description:** 1 = HIV-1 antibody positive at study entry, 0 = HIV-negative. Time-fixed per subject. Used as a binary comorbidity indicator on PK parameters (typically bioavailability or clearance) when a study population pools HIV-positive and HIV-negative subjects on a non-HIV primary indication (tuberculosis treatment, hepatitis treatment, etc.).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (HIV-negative).
- **Source aliases:**
  - `HIV` -- used in `Jonsson_2011_ethambutol.R` (DDMODEL00000220 NMTRAN `$INPUT` column with values 0 = HIV negative, 1 = HIV positive; same orientation as the canonical).
- **Example models:** `Jonsson_2011_ethambutol.R` (multiplicative `1 + e_hiv_pos_f * HIV_POS` shift on bioavailability; HIV-positive patients exhibit a 15.5% reduction in ethambutol bioavailability versus HIV-negative reference).
- **Notes:** Parallels the `_POS` suffix convention used by `ADA_POS`, `SARS_SEROPOS`, and other serostatus / antibody-positivity indicators. Distinct from a primary disease-state indicator like `DIS_HIV` (not yet registered) -- `HIV_POS` is a comorbidity flag in non-HIV-primary indications where HIV-vs-non-HIV is tested as a PK covariate. Ratified canonically on 2026-05-06.

### EARLY_ART (**canonical for early-vs-delayed antiretroviral-treatment-initiation arm indicator**)
- **Description:** Trial randomization-arm indicator: 1 = subject was randomized to initiate antiretroviral treatment (ART) early (within the first 14 days of admission, before nutritional recovery), 0 = subject was randomized to delayed ART initiation (after nutritional recovery, > 14 days from admission). Time-fixed per subject within the trial. The indicator captures the early-vs-delayed-ART contrast tested in the Archary 2019 / MATCH (Malnutrition and ART Timing in Children with HIV) trial in severely malnourished HIV-infected children; the early-ART arm exhibits ~31% higher abacavir bioavailability than the delayed-ART arm.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (delayed ART initiation; F = 1 typical-value reference).
- **Source aliases:**
  - `EARLY` -- Archary 2019 source description ("randomized to early ART"); same orientation as the canonical (1 = early arm).
- **Example models:** `Archary_2019_abacavir.R` (multiplicative additive shift on bioavailability F: `f(depot) <- (1 + e_earlyart_f * EARLY_ART) * exp(etalfdepot)` with `e_earlyart_f = 0.31`; the early-arm IIV on F is 21.4%, the delayed-arm IIV on F is fixed at 0 in the source -- see model-file vignette Errata).
- **Notes:** Specific scope because the early/delayed cutoff (14 days from admission) and the underlying nutritional-rehabilitation context are tied to the MATCH-trial design; future studies that test a similar early-vs-delayed contrast with a different time cutoff or non-nutritional context should register a new canonical (e.g., `EARLY_ART_28D`). Distinct from `TRT_PHASE` (which gates active-vs-baseline study-phase contributions on a per-record basis) and from `DAY14` (which is a within-subject time-since-treatment-initiation landmark, not an enrolment arm). Ratified canonically on 2026-05-08.

## Oncology

### RCFB1MAX (**canonical for week-1 maximum across-lesion relative change from baseline in SUVmax**)
- **Description:** Per-subject scalar predictor entering the overall-survival hazard in the Schindler 2016 sunitinib joint SUVmax / SLD / OS-TTE / dropout model. Defined as the maximum (across the up-to-five tracked target lesions) of `(SUVmax(t = 168 h) - SUVmax(0)) / SUVmax(0)`, i.e., the most negative (largest reduction) relative change in `[18F]FDG-PET` standardized uptake value at one week of sunitinib therapy. The OS Weibull hazard is `lambh * alphh * t^(alphh - 1) * exp(theta_pred * RCFB1MAX)` with `alphh = 1` (degenerates to constant baseline hazard).
- **Units:** unitless (relative change; typically negative for responders).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- entered directly into the OS hazard exponent. Sign convention: more negative `RCFB1MAX` (greater week-1 SUVmax suppression) reduces the OS hazard (Schindler 2016 reports a positive theta_pred = 5.36, so `exp(5.36 * RCFB1MAX)` is < 1 for negative RCFB1MAX).
- **Source aliases:**
  - `RCFB1MAX` (Schindler 2016 NONMEM `$ERROR` block intermediate; "max relative change in SUVmax from baseline at week 1 across lesions"). Computed inline in the source `.mod` from the on-the-fly SUVmax states at TIME = 168 h and reused in subsequent records.
- **Example models:** `Schindler_2016_sunitinib.R` (DDMODEL00000221).
- **Notes:** Specific scope because the metric is tied to the Schindler 2016 GIST-on-sunitinib joint biomarker / OS analysis. In the source NONMEM `.mod` `RCFB1MAX` is a record-loop state, not a true subject-level covariate -- it is captured at FLAG = 1 / TIME = 168 h from the running SUVmax compartment values and reused on every subsequent record. nlmixr2 / rxode2 do not have an idiomatic equivalent of NONMEM's record-loop persistent state, so the model file consumes `RCFB1MAX` as a per-subject input covariate; reproducing the source's behavior requires a two-stage simulation (run the SUVmax + SLD ODEs first, compute `RCFB1MAX` per subject from the t = 168 h SUVmax values, then run the OS / dropout TTE arms with `RCFB1MAX` bound). The vignette virtual cohort follows this pattern.

### TUMSZ (**canonical for baseline tumor size**)
- **Description:** Baseline tumor size. For solid tumors, the sum of diameters of target lesions per RECIST; for classical Hodgkin lymphoma and lymphoma generally, the sum of products of perpendicular diameters (SPPD) or the sum of linear diameters of target lesions, depending on the source paper.
- **Units:** mm (for linear-diameter constructs); mm^2 (for SPPD constructs; record the per-model convention in `covariateData[[TUMSZ]]$units` and `notes`).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(TUMSZ / ref)^exponent` for continuous effects, or with a paper-specific threshold for categorical-stratum indicators (e.g., Gibiansky 2014 splits BSIZ at 1750 mm^2 into low- vs high-burden strata). Reference values observed: 41 mm (Zhou 2025); 63 mm (Budha 2023); 90 mm (Lu 2014, source reference 9 cm converted to mm); 1750 mm^2 threshold (Gibiansky 2014, SPPD).
- **Source aliases:**
  - `LDIAM` (Zhou 2025; pediatric lymphoma "linear diameter" of target lesions in mm).
  - `TMBD` (originally in cm; `TUMSZ_mm = TMBD_cm * 10`) -- used in `Lu_2014_trastuzumabemtansine.R`.
  - `BSIZ` (Gibiansky 2014; baseline tumor size as the sum of products of perpendicular diameters of target lesions, mm^2; used as the categorical indicator `(BSIZ <= 1750)` in the obinutuzumab popPK model rather than as a continuous power covariate).
- **Example models:** `Budha_2023_tislelizumab.R` (reference 63 mm), `Lu_2014_trastuzumabemtansine.R` (reference 90 mm; source column TMBD in cm, values converted to mm on ingestion), `Zhou_2025_brentuximab.R` (reference 41 mm; source column LDIAM is the sum of linear diameters of target lesions; effect on ADC clearance only), `Gibiansky_2014_obinutuzumab.R` (SPPD in mm^2; used as a categorical indicator `(TUMSZ <= 1750)` on the time-dependent clearance decay rate kdes, not as a continuous power covariate), `Hansson_2013b_sunitinib.R` (DDMODEL00000198; observed baseline tumor SLD used as the per-subject IC of the tumor-size ODE via the IPP-style proportional baseline-residual construction `tumor(0) = TUMSZ * (1 + etaibase * propSd)`; the source .mod reads OBASE from DV at TIME=0/FLAG=4, but nlmixr2 / rxode2 cannot replicate the in-record assignment idiom and consumes the observed baseline as a covariate instead).
- **Notes:** Promoted to scope: general on 2026-04-20 as a conventional oncology baseline-tumor-size measure (RECIST for solid tumors, SPPD or sum-of-linear-diameters for lymphomas). The SPPD vs sum-of-diameters vs sum-of-linear-diameters convention is pooled onto a single column; document the per-model mixture where relevant. When the source paper reports tumor size in cm, convert to mm (the canonical unit) on data ingestion and scale the per-model reference accordingly so `(TUMSZ / ref)^exp` is numerically invariant. For SPPD constructs the natural unit is mm^2 (a product of two perpendicular diameters in mm); record that in the per-model `covariateData[[TUMSZ]]$units` field and do NOT cross-mix mm and mm^2 within a single ingest. When a source paper specifically reports the RECIST 1.1 "sum of longest diameters" of target lesions, use the more specific `TUM_SLD` canonical instead -- `TUMSZ` remains the pooled-tumor-burden register.

### TUM_SLD (**canonical for sum of longest diameters of target lesions**)
- **Description:** Baseline sum of longest diameters of target lesions per RECIST 1.1. More specific than the pooled `TUMSZ` canonical; use `TUM_SLD` when the source paper explicitly reports "sum of longest diameters" (or "sum of lesions") as the tumor-burden metric, distinct from the pooled "sum of diameters / SPPD / sum of linear diameters" mixture covered by `TUMSZ`.
- **Units:** mm
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(TUM_SLD / ref)^exponent`. Reference values observed: 70.0 mm (de Vries Schultink 2020 zenocutuzumab population median).
- **Source aliases:**
  - `SoL` / "sum of lesions" (de Vries Schultink 2020 zenocutuzumab) -- same construct, mm.
- **Example models:** `deVriesSchultink_2020_zenocutuzumab.R` (reference 70.0 mm; power exponent 0.447 on Vmax of the parallel non-linear / Michaelis-Menten clearance).
- **Notes:** Distinct from `TUMSZ` (pooled tumor-size canonical covering RECIST sum-of-diameters / SPPD / sum-of-linear-diameters); `TUM_SLD` is the precise RECIST 1.1 sum-of-longest-diameters metric. Ratified canonically on 2026-04-29 alongside the pilot bispecific extraction (de Vries Schultink 2020 zenocutuzumab). When the source paper reports tumor burden in cm, convert to mm (the canonical unit) on data ingestion and scale the per-model reference accordingly so `(TUM_SLD / ref)^exp` is numerically invariant. Also used as the per-subject initial-condition input for tumour-growth / angiogenesis-inhibition (TGI) ODE models where the source paper sets the SLD state at time zero from the observed baseline SLD rather than estimating a typical-value baseline (e.g., `Ouerdani_2015_pazopanib.R` uses `tumorSize(0) <- TUM_SLD`).

### TUM_VOL (**canonical for caliper-measured baseline tumor volume in xenograft / preclinical studies**)
- **Description:** Baseline tumor volume measured by handheld caliper in subcutaneous-xenograft preclinical studies and similar small-animal models, computed from longest length and orthogonal width via `volume = (length * width^2) / 2` (the standard ellipsoid-approximation formula used in xenograft pharmacology). Per-subject (per-animal) baseline measurement at randomisation into a dosing group. Used both as a covariate stratifier (size at randomisation) and as the per-subject initial-condition input for TGI ODE models where the source paper sets the tumor-volume state at time zero from the observed measurement rather than estimating a typical-value baseline.
- **Units:** `mm^3`
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- typically used directly as the ODE initial state (`tumorSize(0) <- TUM_VOL`); not a covariate effect coefficient.
- **Source aliases:**
  - `P0` (Ouerdani 2015 mouse model symbol for the observed initial tumour volume; mm^3 in the CAKI-2 xenograft cohort; range 100-250 mm^3 at randomisation per the paper's preclinical Methods).
- **Example models:** `Ouerdani_2015_pazopanib_mouse.R` (Ouerdani 2015 preclinical TGI in CAKI-2 xenograft mice; per-subject `TUM_VOL` initialises the `tumorSize` state and is held constant per individual across the 24-day dosing window).
- **Notes:** General scope because xenograft tumor-volume baselines have a shared meaning across preclinical TGI papers regardless of the drug or cell line. Use `TUM_VOL` whenever a preclinical paper supplies per-animal caliper tumor volumes as the per-subject initial state of a tumor-volume ODE; use `TUM_SLD` for clinical RECIST sum-of-longest-diameters (a length, not a volume) and `TUMSZ` for the pooled "baseline tumor burden as a covariate on PK" use case in clinical models. Ratified canonically on 2026-05-12 alongside the Ouerdani 2015 pazopanib mouse extraction.

### TUMTP_CHL (**canonical for classical Hodgkin lymphoma tumor-type indicator**)
- **Description:** 1 = classical Hodgkin lymphoma (cHL) or Hodgkin lymphoma generally, 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (e.g., NSCLC, EC, HCC, UC, GC, CRC, NPC, OC, "Other" solid tumors in the Budha 2023 cohort; systemic anaplastic large-cell lymphoma in the Zhou 2025 pediatric cohort).
- **Source aliases:**
  - `TUMTP` (categorical column with levels like `cHL`, `GC`, ...) -- decompose into `TUMTP_CHL = as.integer(TUMTP == "cHL")`.
  - `DIS` (Zhou 2025; integer code with `DIS == 1` flagging HL) -- decompose into `TUMTP_CHL = as.integer(DIS == 1)`. Zhou 2025 calls the complement "non-HL"; in the Zhou 2025 cohort the non-HL group is exclusively sALCL.
- **Example models:** `Budha_2023_tislelizumab.R`, `Zhou_2025_brentuximab.R` (effects on ADC Q2, MMAE central volume VM, and the ADC->MMAE conversion-decay rate ALFM; the Zhou 2025 paper anchors typical-value parameters to HL patients so the model uses `(1 - TUMTP_CHL)` as the on-effect indicator with reference category 1 = HL).
- **Notes:** Paired with `TUMTP_GC` in Budha 2023; a patient can have at most one of the indicators set to 1 (the remaining tumor types collapse into the reference 0 group). The reference category is the off-encoded value (0) by definition; when a source paper anchors typical-value parameters to the HL group rather than the non-HL group (as Zhou 2025 does), encode the effect as `coef^(1 - TUMTP_CHL)` so the canonical column meaning (1 = cHL/HL) is preserved while the paper's reference (HL) still receives multiplier 1.

### HER2_ECD (**canonical for HER2 shed extracellular domain concentration**)
- **Description:** Baseline serum concentration of the shed extracellular domain of human epidermal growth factor receptor 2 (HER2). Serves as a soluble-antigen biomarker of HER2-mediated target-mediated drug disposition for HER2-directed mAbs / ADCs.
- **Units:** ng/mL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with power scaling `(HER2_ECD / ref)^exponent`. Reference 25 ng/mL used in Lu 2014; reference 8.23 ng/mL (population median) used in Bruno 2005, with a 200 ng/mL plateau cap.
- **Source aliases:**
  - `ECD` -- used in `Lu_2014_trastuzumabemtansine.R` and `Bruno_2005_trastuzumab.R`.
- **Example models:** `Lu_2014_trastuzumabemtansine.R` (reference 25 ng/mL; exponent 0.035 on CL), `Bruno_2005_trastuzumab.R` (reference 8.23 ng/mL; power exponents 0.041 on CL and 0.105 on V; HER2_ECD capped at 200 ng/mL inside model() to reflect the Bruno 2005 plateau observation).
- **Notes:** Scoped specific because the covariate is meaningful only for HER2-targeted agents; if a non-HER2 paper uses a shed-antigen analog for a different target, register a target-specific canonical (e.g., `EGFR_ECD`) rather than reusing this one. Disambiguated from the covariate-columns register by the explicit `HER2_` prefix.

### TRAST_BL (**canonical for baseline trastuzumab concentration from prior therapy**)
- **Description:** Baseline serum concentration of residual unconjugated trastuzumab remaining from prior trastuzumab-containing therapy, measured at the start of a subsequent anti-HER2 treatment (e.g., trastuzumab emtansine). Encodes the magnitude of residual HER2-site competition from a previous trastuzumab exposure.
- **Units:** ug/mL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used as a linear covariate on log(CL) via `exp(coef * TRAST_BL)`. TRAST_BL = 0 corresponds to no detectable residual trastuzumab (reference condition used in Lu 2014).
- **Source aliases:**
  - `TBL` -- used in `Lu_2014_trastuzumabemtansine.R`.
- **Example models:** `Lu_2014_trastuzumabemtansine.R` (linear-on-log coefficient -0.002 per ug/mL on CL).
- **Notes:** Scoped specific because the covariate is meaningful only for drugs that compete with trastuzumab at the HER2 binding site. Expect values clustered at 0 for trastuzumab-naive patients; Lu 2014 observed 0 at the 5th percentile and 54 ug/mL at the 95th percentile.

### TUMTP_GC (**canonical for gastric-cancer tumor-type indicator**)
- **Description:** 1 = gastric cancer (GC) or adenocarcinoma of the gastroesophageal junction (GEJ), 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (same reference group as `TUMTP_CHL`).
- **Source aliases:**
  - `TUMTP` (categorical column) -- decompose into `TUMTP_GC = as.integer(TUMTP == "GC")`.
  - `TTYPE` (Quartino 2019; categorical column with levels `MBC`, `EBC`, `HV`, `AGC`, `Others`) -- decompose into `TUMTP_GC = as.integer(TTYPE == "AGC")`.
  - `TTYPE4` (Wang 2024; level 4 of a five-level tumor-type factor labelled "GCGEJ" in the source) -- decompose into `TUMTP_GC = as.integer(TTYPE4 == 1)`.
- **Example models:** `Budha_2023_tislelizumab.R`, `Quartino_2019_trastuzumab.R` (advanced gastric cancer; per-group typical-value switch on linear CL and Vc rather than an exponential multiplier), `Wang_2024_sugemalimab.R` (gastric + GEJ adenocarcinoma pooled; exponential coefficient log(1.13) on CL and log(1.14) on Vc).
- **Notes:** Follows the `RACE_<GROUP>` indicator-decomposition pattern. New oncology tumor types should be added as additional `TUMTP_<GROUP>` entries so the reference set stays explicit. "Advanced gastric cancer" (AGC), "gastric cancer" (GC), and "GC or adenocarcinoma of the gastroesophageal junction" (GCGEJ) are pooled onto a single `TUMTP_GC` indicator; document the per-paper stage-of-disease and GEJ-inclusion detail in `covariateData[[TUMTP_GC]]$notes`. ESCC (squamous histology) is captured by the separate `TUMTP_ESCC` indicator and is not pooled here.

### TUMTP_OTH (**canonical for 'other tumor types' residual indicator**)
- **Description:** 1 = heterogeneous "other" tumor-type pool (typically NSCLC plus miscellaneous solid tumors such as prostate, ovarian, and colorectal), 0 = one of the named tumor-type groups in the same analysis.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = named tumor-type groups defined per-paper (e.g., MBC, EBC, HV, AGC in Quartino 2019). The complement of all `TUMTP_<GROUP>` indicators defined in the same model.
- **Source aliases:**
  - `TTYPE` (Quartino 2019) -- decompose into `TUMTP_OTH = as.integer(TTYPE == "Others")`.
  - `PAT2` (Sathe 2024) -- integer-coded tumor type column with levels 1 (mTNBC), 2 (mUC or HR+/HER2- mBC), 4 (Other epithelial); the source NONMEM control stream collapses PAT2 = 1 and PAT2 = 2 into the reference (no effect) and applies the deviation only when PAT2 = 4. Decompose into `TUMTP_OTH = as.integer(PAT2 == 4)`.
- **Example models:** `Quartino_2019_trastuzumab.R` (per-group typical-value switch on linear CL; NSCLC plus a small residual group of prostate, ovarian, and other histologies), `Sathe_2024_sacituzumab.R` (multiplicative effect on tAB CL: -13.4% when TUMTP_OTH = 1; "Other" pool = small-cell and non-small-cell lung cancer, colorectal cancer, esophageal cancer, pancreatic ductal adenocarcinoma, etc., n = 184; reference = pooled mTNBC + mUC + HR+/HER2- mBC, n = 345).
- **Notes:** Scope: specific because the set of histologies collapsed into "Others" is defined by the analysis plan of the source paper; two papers' `TUMTP_OTH` columns are not interchangeable. Document the exact per-paper composition (e.g., "NSCLC + prostate + ovarian + other, n = 107 in Quartino 2019"; "small-cell + non-small-cell lung + CRC + esophageal + pancreatic ductal adenocarcinoma, n = 184 in Sathe 2024") in `covariateData[[TUMTP_OTH]]$notes`. A given subject can have at most one of the `TUMTP_<GROUP>` indicators (including `TUMTP_OTH`) set to 1; all-zero means the reference group.

### SPDL1 (**canonical for soluble PD-L1 concentration**)
- **Description:** Baseline (or time-varying) serum concentration of soluble programmed death-ligand 1 (sPD-L1). Serves as a circulating biomarker of target burden and immune activation for anti-PD-1/PD-L1 antibodies.
- **Units:** pg/mL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with power scaling `(SPDL1 / ref)^exponent`. Reference value observed: 173.8 pg/mL (study-population median in Ogasawara 2020).
- **Source aliases:** none; `SPDL1` is the standard abbreviation used directly in source analyses.
- **Example models:** `Ogasawara_2020_durvalumab.R` (power effect on CL, exponent 0.0617, reference 173.8 pg/mL; time-varying; values below LOD imputed as LOD/2 = 33.55 pg/mL), `Quartino_2019_trastuzumab.R` (per-group typical-value switch on linear CL; NSCLC plus a small residual group of prostate, ovarian, and other histologies), `Wang_2024_sugemalimab.R` (heterogeneous solid-tumor residual group of n = 174; exponential coefficient log(0.885) on CL and log(0.926) on Vc; NSCLC is the reference group, not part of `TUMTP_OTH`).
- **Notes:** Scope: specific because sPD-L1 is meaningful only for drugs targeting the PD-1/PD-L1 pathway. For other checkpoint biomarkers (e.g., soluble CTLA-4, soluble LAG-3) register new dedicated canonicals rather than reusing this one. Ratified canonically on 2026-04-26.
  - `TTYPE3` (Wang 2024; level 3 of a five-level tumor-type factor labelled "Other" in the source) -- decompose into `TUMTP_OTH = as.integer(TTYPE3 == 1)`.
- **Example models:** `Quartino_2019_trastuzumab.R` (per-group typical-value switch on linear CL; NSCLC plus a small residual group of prostate, ovarian, and other histologies), `Wang_2024_sugemalimab.R` (heterogeneous solid-tumor residual group of n = 174; exponential coefficient log(0.885) on CL and log(0.926) on Vc; NSCLC is the reference group, not part of `TUMTP_OTH`), `Ogasawara_2020_durvalumab.R` (power effect on CL, exponent 0.0617, reference 173.8 pg/mL; time-varying; values below LOD imputed as LOD/2 = 33.55 pg/mL).
- **Notes:** Scope: specific because the set of histologies collapsed into "Others" is defined by the analysis plan of the source paper; two papers' `TUMTP_OTH` columns are not interchangeable. Document the exact per-paper composition (e.g., "NSCLC + prostate + ovarian + other, n = 107 in Quartino 2019"; "miscellaneous solid tumors excluding NSCLC, lymphoma, GCGEJ, and ESCC, n = 174 in Wang 2024") in `covariateData[[TUMTP_OTH]]$notes`. A given subject can have at most one of the `TUMTP_<GROUP>` indicators (including `TUMTP_OTH`) set to 1; all-zero means the reference group.

### MCPROT (**canonical for serum monoclonal (M) protein concentration**)
- **Description:** Serum monoclonal (M) protein concentration. Multiple-myeloma plasma-cell-burden marker secreted by the tumor clone; elevated MCPROT reflects higher tumor burden and (for IgG-secreting MM) competes with therapeutic IgG mAbs for FcRn-mediated salvage and target-mediated elimination. Typically time-varying -- measured at multiple visits over the treatment course and supplied at every PK observation time via linear interpolation between measurements.
- **Units:** g/dL (US-convention; equivalent to 10 g/L SI). Document the unit used in each model via `covariateData[[MCPROT]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used as a continuous, log-linear effect on the Vmax of target-mediated elimination via `exp(theta * MCPROT)` (i.e., MCPROT enters un-log-transformed). Reference values observed: 0 g/dL (Ide 2020 Vmax,REF reference) and 2.0 g/dL (Ide 2020 figure-1 reference patient).
- **Source aliases:**
  - `TMCPROT` (time-varying serum M-protein concentration) -- used in `Ide_2020_elotuzumab.R`. NONMEM column with imputation sentinel `-99` for missing observations, replaced by population median 2.1 g/dL via `IF(TMCPROT.EQ.-99) TMCPROT = 2.1`.
- **Example models:** `Ide_2020_elotuzumab.R` (g/dL, time-varying; entered un-log-transformed as `exp(0.277 * MCPROT)` on Vmax of the Michaelis-Menten target-mediated elimination from the central compartment).
- **Notes:** Specific scope because the column is mechanistically meaningful only for plasma-cell-targeting therapies in multiple myeloma (e.g., elotuzumab anti-SLAMF7, daratumumab anti-CD38, isatuximab anti-CD38, belantamab anti-BCMA, and BCMA-bispecifics / CAR-T). MCPROT decreases with treatment response; the time-varying form is the only correct way to capture the diminishing target-mediated-elimination component as the tumor regresses. In NONMEM datasets MCPROT is supplied at each event-row time, with linear interpolation between observations and last-observation-carried-forward beyond the last sample (Ide 2020 Methods). Distinct from `MM_NIGG` (which is the immunoglobulin subtype, an MM-disease stratifier that is time-fixed), `SBCMA` (soluble BCMA, a different MM tumor-burden biomarker for BCMA-targeting drugs), and `B2M` (beta-2-microglobulin, a renal-function-and-MM-disease-burden marker). The 1 g/dL = 10 g/L conversion lets future SI-convention papers register the same canonical with their own unit string.

### LMET (**canonical for baseline presence of liver metastases**)
- **Description:** Binary indicator of radiologically documented liver metastases at baseline, 1 = liver metastases present, 0 = no liver metastases.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no liver metastases at baseline).
- **Source aliases:** none; `LMET` is the common NONMEM / clinical-PK abbreviation used directly by the source papers.
- **Example models:** `Quartino_2019_trastuzumab.R` (exponential effect on linear CL; +16.4% typical CL when LMET = 1, per Quartino 2019 Table 1 theta12 = 0.152).
- **Notes:** Liver metastases are associated with hepatic protein-synthesis impairment and altered IgG catabolism, making `LMET` a commonly tested covariate in oncology mAb population PK analyses. Scope: general so future oncology papers can reuse the canonical column. Time-fixed baseline indicator; if a source paper treats it as time-varying (progression during treatment), document in `covariateData[[LMET]]$notes`.

### MET_GE4 (**canonical for baseline number of metastatic sites >= 4 indicator**)
- **Description:** Binary indicator dichotomising the count of baseline metastatic sites at 4, 1 = patient has four or more documented metastatic sites at baseline, 0 = patient has zero to three metastatic sites at baseline. Time-fixed per subject. Treated as a surrogate for tumor burden in oncology mAb popPK analyses.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (fewer than four metastatic sites at baseline).
- **Source aliases:**
  - `MET` -- used in `Bruno_2005_trastuzumab.R` (Bruno 2005 Methods: "MET=1 if number of metastatic sites=4 or greater; otherwise MET=0"). When a source paper supplies the raw integer count column (`NMET`, `N_METS`, etc.) rather than the pre-binarised indicator, derive `MET_GE4 = as.integer(NMET >= 4)`.
- **Example models:** `Bruno_2005_trastuzumab.R` (multiplicative effect on linear CL: typical CL multiplied by `(1 + 0.221 * MET_GE4)`, i.e. +22.1% CL in MET_GE4 = 1 patients; reference MET_GE4 = 0 per Bruno 2005 Table 3).
- **Notes:** The >= 4 split is the dichotomisation used by Bruno 2005 to capture "patients with four or more metastatic sites" as a high-tumor-burden subgroup. Scope: general so future oncology popPK papers using the same dichotomisation can reuse this canonical column. If a future paper uses a different split (e.g. >= 2 or >= 5), register a separate `MET_GEN` canonical rather than overloading this entry. A `MET_GE4 = 1` patient may also have `LMET = 1`; the two columns are not mutually exclusive (liver-only metastases would have `LMET = 1`, `MET_GE4 = 0`).

### ECOG_GE1 (**canonical for Eastern Cooperative Oncology Group performance-status indicator, >= 1**)
- **Description:** 1 if baseline Eastern Cooperative Oncology Group (ECOG) performance status is greater than or equal to 1, 0 if ECOG = 0. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (ECOG performance status = 0, i.e., fully active / asymptomatic).
- **Source aliases:**
  - `PS` / `BPS` -- used in `Bajaj_2017_nivolumab.R` (BPS = "baseline performance status"; the one study using Karnofsky Performance Status was mapped to ECOG via Oken 1982 before binarization) and `Zhang_2019_nivolumab.R` (paper's binary collapse PS=0 vs. PS>0).
  - `ECOG_1` -- alternative explicit form; equivalent to `ECOG_GE1` when ECOG only takes values 0, 1, 2 in the analysis dataset (the typical oncology case).
  - `ECOG_PS_GT0` -- retired name used in earlier register drafts; semantically identical (`>= 1` equals `> 0` for integer ECOG scores).
  - `ECOG101` (categorical 0/1/2 score with thresholding `IF(ECOG101.GT.0.5)`) -- used in `Ide_2020_elotuzumab.R`. Decompose: `ECOG_GE1 = as.integer(ECOG101 >= 1)`.
- **Example models:** `Bajaj_2017_nivolumab.R` (exponential effect on CL with coefficient 0.172), `Zhang_2019_nivolumab.R` (exponential effect exp(0.181) on baseline CL; additive effect -0.138 on the time-varying-CL Emax parameter), `Ide_2020_elotuzumab.R` (multiplicative effect on CL = 1.03; paired with `ECOG_GE2` for separate ECOG=1 vs ECOG>=2 effects), `Netterberg_2017_docetaxel.R` (multiplicative effect on baseline ANC of the Friberg myelosuppression chain: `BACOV *= (1 + theta * ECOG_GE1)` with theta = 0.130; source column `PERF` with ordinal ECOG 0/1/2 values, binarized via `ECOG_GE1 = as.integer(PERF >= 1)` per Kloft 2006).
- **Notes:** Oncology papers conventionally report ECOG as an integer (0-5) but binarize at >= 1 because ECOG >= 2 is rare in trial cohorts. When a source paper provides the ordinal ECOG score separately, derive `ECOG_GE1 = as.integer(ECOG >= 1)`. Zhang 2019 uses `ECOG_GE1` on both baseline CL and the time-varying Emax parameter (unlike Bajaj 2017, which uses it on CL only); document the structural role in each model's `covariateData[[ECOG_GE1]]$notes`. When a paper retains separate effects for ECOG = 1 vs ECOG >= 2 (Ide 2020), pair this column with `ECOG_GE2` and supply both indicators in the event dataset.

### ECOG_GE2 (**canonical for Eastern Cooperative Oncology Group performance-status indicator, >= 2**)
- **Description:** 1 if baseline Eastern Cooperative Oncology Group (ECOG) performance status is greater than or equal to 2, 0 if ECOG <= 1. Time-fixed per subject. Used in models that retain separate effects for ECOG = 1 vs ECOG >= 2 by pairing this column with `ECOG_GE1`.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (ECOG performance status <= 1; in models that pair `ECOG_GE1` and `ECOG_GE2`, both indicators = 0 corresponds to ECOG = 0 and (`ECOG_GE1` = 1, `ECOG_GE2` = 0) corresponds to ECOG = 1).
- **Source aliases:**
  - `ECOG101` (categorical 0/1/2 score with thresholding `IF(ECOG101.GT.1.5)`) -- used in `Ide_2020_elotuzumab.R`. Decompose: `ECOG_GE2 = as.integer(ECOG101 >= 2)`.
- **Example models:** `Ide_2020_elotuzumab.R` (multiplicative effect on CL = 1.15; paired with `ECOG_GE1` to retain separate ECOG = 1 vs ECOG >= 2 effects).
- **Notes:** Parallels `ECOG_GE1`. Use only when the source paper reports a separate effect for ECOG >= 2 in addition to ECOG_GE1; otherwise `ECOG_GE1` alone is sufficient. The paired (`ECOG_GE1`, `ECOG_GE2`) decomposition reproduces a three-level (`ECOG = 0`, `ECOG = 1`, `ECOG >= 2`) ordinal effect with two binaries.

### TUMTP_SCLC (**canonical for small-cell-lung-cancer tumor-type indicator**)
- **Description:** 1 = small cell lung cancer (SCLC), 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (e.g., melanoma, NSCLC, RCC, HCC, CRC in the Sanghavi 2020 cohort; reference category is melanoma).
- **Source aliases:**
  - `TUMTP` (categorical column with levels including `melanoma`, `NSCLC`, `SCLC`, `CRC`, `HCC`, `RCC`) -- decompose into `TUMTP_SCLC = as.integer(TUMTP == "SCLC")`.
- **Example models:** `Sanghavi_2020_ipilimumab.R` (exponential coefficient -0.124 on CL).
- **Notes:** Follows the `TUMTP_CHL` / `TUMTP_GC` decomposition pattern. SCLC is the only retained tumor-type indicator in the Sanghavi 2020 final model after backward elimination; the other tumor types collapse into the reference (melanoma) group.

### TUMTP_LYMPH (**canonical for lymphoma (pooled) tumor-type indicator**)
- **Description:** 1 = lymphoma (heterogeneous lymphoma pool spanning multiple lymphoma histologies -- e.g., classical Hodgkin lymphoma combined with extranodal NK/T-cell lymphoma; or any-histology lymphoma pooled with solid-tumor and leukemia cohorts), 0 = solid tumor or other tumor type.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 = non-lymphoma tumor type (per source paper; e.g., NSCLC, GC/GEJ, ESCC, "Other" solid tumors in the Wang 2024 cohort, with NSCLC as the implicit reference when paired with the other Wang 2024 `TUMTP_*` indicators; or leukemia as the implicit reference in the Akbar 2025 cohort).
- **Source aliases:**
  - `TTYPE1` (Wang 2024) -- decompose into `TUMTP_LYMPH = as.integer(TTYPE1 == 1)`. The Wang 2024 source paper uses a multi-level `TTYPE` factor with levels 1 = lymphoma, 2 = lung cancer (reference), 3 = other, 4 = GCGEJ, 5 = ESCC.
  - Categorical column "type of cancer" with level "Lymphoma" (Akbar 2025) -- decompose into `TUMTP_LYMPH = as.integer(cancer_type == "Lymphoma")`.
- **Example models:** `Wang_2024_sugemalimab.R` (exponential coefficient log(0.877) on baseline CL and log(0.879) on Vc), `Akbar_2025_voriconazole.R` (additive-fractional +1.91% effect on CL relative to leukemia reference; 95% CI spans zero).
- **Notes:** Distinct from `TUMTP_CHL` (which is specifically classical Hodgkin lymphoma). Wang 2024 pools two lymphoma histologies (extranodal NK/T-cell lymphoma from CS1001-201 / NCT03595657 and classical Hodgkin lymphoma from CS1001-202 / NCT03505996) into a single lymphoma indicator; the indicator therefore captures a generic "hematologic-vs-solid-tumor" contrast rather than a histology-specific effect. Akbar 2025 uses a single "Lymphoma" category alongside leukemia, sarcoma, breast cancer, myeloma, and glioma in a heterogeneous-cancer TDM cohort. When a future paper studies a single lymphoma histology distinct from cHL, register a more specific canonical (e.g., `TUMTP_ENKTL`, `TUMTP_NHL`) rather than overloading this one. Document the per-paper histology composition in `covariateData[[TUMTP_LYMPH]]$notes`. Promoted from `Scope: specific` to `Scope: general` on 2026-05-09 alongside the Akbar 2025 voriconazole extraction so that any heterogeneous-cancer-cohort PK analysis can use this canonical name without scope-violation.

### TUMTP_BC (**canonical for breast-cancer tumor-type indicator**)
- **Description:** 1 = breast cancer (any histology / receptor status), 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 = all other tumor types (per source paper; in Lu 2022 the implicit reference is NSCLC, with colorectal cancer pooled into the reference because its CL effect was found insignificant relative to NSCLC).
- **Source aliases:**
  - `TUMTP` (categorical column with levels including `BC`, `NSCLC`, `CRC`) -- decompose into `TUMTP_BC = as.integer(TUMTP == "BC")`.
- **Example models:** `Lu_2022_patritumab.R` (multiplicative fractional effect 0.811 on CLlin of DXd-conjugated antibody for breast-cancer patients vs the NSCLC reference; CRC effect was tested and found insignificant so CRC is pooled into the reference).
- **Notes:** Follows the `TUMTP_CHL` / `TUMTP_GC` / `TUMTP_SCLC` decomposition pattern. Registers the breast-cancer arm of an oncology-cohort tumor-type contrast; pair with sister `TUMTP_<GROUP>` indicators (e.g., `TUMTP_NSCLC`, `TUMTP_CRC`) when a future paper retains separate effects for additional tumor types beyond the implicit reference. Ratified canonically on 2026-04-28.

### TUMTP_ESCC (**canonical for oesophageal-squamous-cell-carcinoma tumor-type indicator**)
- **Description:** 1 = oesophageal squamous cell carcinoma (ESCC), 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (per source paper; in Wang 2024, the implicit reference is NSCLC when all the other `TUMTP_*` indicators are also 0).
- **Source aliases:**
  - `TTYPE5` (Wang 2024) -- decompose into `TUMTP_ESCC = as.integer(TTYPE5 == 1)`.
- **Example models:** `Wang_2024_sugemalimab.R` (exponential coefficient log(0.99) on baseline CL and log(1.08) on Vc).
- **Notes:** Follows the `TUMTP_CHL` / `TUMTP_GC` / `TUMTP_SCLC` decomposition pattern. Distinct from gastroesophageal-junction adenocarcinoma (which is captured by the broader `TUMTP_GC` indicator that pools GC and GEJ adenocarcinomas) -- ESCC is a squamous-cell histology, not adenocarcinoma. Document the per-paper histology composition in `covariateData[[TUMTP_ESCC]]$notes`.

### TUMTP_PCALCL (**canonical for primary cutaneous anaplastic large-cell lymphoma indicator**)
- **Description:** 1 = primary cutaneous anaplastic large-cell lymphoma (pcALCL), 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (in Suri 2018, the non-pcALCL reference comprises Hodgkin lymphoma, systemic ALCL, mycosis fungoides, and other CD30+ hematologic malignancies pooled together).
- **Source aliases:**
  - `PCALCL` -- used in `Suri_2018_brentuximab.R`. Suri 2018 reports the effect as a power-form multiplier `cl_adc *= 0.728^TUMTP_PCALCL` (pcALCL CL ~27% lower than non-pcALCL).
- **Example models:** `Suri_2018_brentuximab.R` (effect on ADC clearance only).
- **Notes:** Follows the `TUMTP_CHL` / `TUMTP_GC` / `TUMTP_SCLC` decomposition pattern. Distinct from `TUMTP_LYMPH` (heterogeneous lymphoma pool) and `TUMTP_CHL` (classical Hodgkin lymphoma). pcALCL is one of two histologies pooled into the broader CTCL category in Suri 2018 (alongside mycosis fungoides); the model singles out pcALCL because Suri 2018 backward elimination retained pcALCL as a separate effect on ADC clearance after exploring the broader CTCL contrast. Ratified canonically on 2026-04-28.

### TUMTP_SARC (**canonical for sarcoma tumor-type indicator**)
- **Description:** 1 = sarcoma (any histology -- soft-tissue or bone sarcoma pooled), 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (per source paper; in Akbar 2025 the implicit reference is leukemia when paired with the other Akbar `TUMTP_*` indicators all = 0).
- **Source aliases:**
  - Categorical column "type of cancer" with level "Sarcoma" -- decompose into `TUMTP_SARC = as.integer(cancer_type == "Sarcoma")`. Used in `Akbar_2025_voriconazole.R`.
- **Example models:** `Akbar_2025_voriconazole.R` (additive-fractional +18.5% effect on CL relative to leukemia reference).
- **Notes:** Follows the `TUMTP_CHL` / `TUMTP_GC` / `TUMTP_SCLC` decomposition pattern. Akbar 2025 pools soft-tissue and bone sarcoma histologies into a single sarcoma category. Scope: specific because the reference category (leukemia in Akbar 2025) is paper-defined. Ratified canonically on 2026-05-09.

### TUMTP_MYELO (**canonical for multiple myeloma tumor-type indicator (used in heterogeneous-cancer pooled cohorts)**)
- **Description:** 1 = multiple myeloma, 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (per source paper; in Akbar 2025 the implicit reference is leukemia when paired with the other Akbar `TUMTP_*` indicators all = 0).
- **Source aliases:**
  - Categorical column "type of cancer" with level "Myeloma" -- decompose into `TUMTP_MYELO = as.integer(cancer_type == "Myeloma")`. Used in `Akbar_2025_voriconazole.R`.
- **Example models:** `Akbar_2025_voriconazole.R` (additive-fractional -2.33% effect on CL relative to leukemia reference; the 95% CI spans zero).
- **Notes:** Distinct from the stub `MM` entry (which is reserved for multiple-myeloma-as-primary-indication PK studies; the `MM` definition lacks a complete schema and predates the TUMTP_* convention) and from `DIS_SMM` (smoldering multiple myeloma, an asymptomatic plasma-cell disorder). Use `TUMTP_MYELO` when the source paper pools multiple myeloma alongside other tumor types in a heterogeneous oncology cohort and treats `cancer type` as a many-level categorical covariate. Scope: specific because the reference category (leukemia in Akbar 2025) is paper-defined. Ratified canonically on 2026-05-09.

### TUMTP_GLIO (**canonical for glioma tumor-type indicator**)
- **Description:** 1 = glioma (any grade / histology -- e.g., glioblastoma, anaplastic astrocytoma, oligodendroglioma pooled), 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (per source paper; in Akbar 2025 the implicit reference is leukemia when paired with the other Akbar `TUMTP_*` indicators all = 0).
- **Source aliases:**
  - Categorical column "type of cancer" with level "Glioma" -- decompose into `TUMTP_GLIO = as.integer(cancer_type == "Glioma")`. Used in `Akbar_2025_voriconazole.R`.
- **Example models:** `Akbar_2025_voriconazole.R` (additive-fractional +8.81% effect on CL relative to leukemia reference; the 95% CI spans zero).
- **Notes:** Follows the `TUMTP_CHL` / `TUMTP_GC` / `TUMTP_SCLC` decomposition pattern. Akbar 2025 reports a single "glioma" category without further subdivision by histology or grade. Scope: specific because the reference category (leukemia in Akbar 2025) is paper-defined. Ratified canonically on 2026-05-09.

### TUMTP_LEUK (**canonical for leukemia tumor-type indicator**)
- **Description:** 1 = leukemia (any subtype -- AML / ALL / CLL / CML pooled), 0 = other tumor types. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (per source paper). In Akbar 2025 leukemia is the implicit reference (column not used directly), but the canonical name is registered so that future papers that retain leukemia as a non-reference contrast can use it.
- **Source aliases:**
  - Categorical column "type of cancer" with level "Leukemia" -- decompose into `TUMTP_LEUK = as.integer(cancer_type == "Leukemia")`. Implicit reference category in `Akbar_2025_voriconazole.R` (so the model file does not consume this column directly; it is registered for future heterogeneous-cancer-cohort analyses).
- **Example models:** `Akbar_2025_voriconazole.R`.
- **Notes:** Distinct from the more specific `DIS_AML`, `DIS_BCPALL`, `DIS_CMML`, `MDSAML` entries -- those are for leukemia-only or leukemia-vs-leukemia contrasts; `TUMTP_LEUK` is for heterogeneous-cancer pooled cohorts where leukemia is one of several tumor types and the analysis treats `cancer type` as a many-level categorical. Akbar 2025 had leukemia as 56.8% of the cohort and used it as the reference category. Scope: specific because the reference category in any source paper is paper-defined. Ratified canonically on 2026-05-09.

### TUMTP_BCL (**canonical for B-cell lymphoma (pooled residual) tumor-type indicator**)
- **Description:** 1 = B-cell lymphoma (BCL), 0 = other tumor types. Time-fixed per subject. In Gibiansky 2014 the BCL category is a pooled residual indolent-B-cell-lymphoma group that includes follicular lymphoma (FL was the primary indication in GAUDI; the four-level DIS column in the NONMEM control stream splits B-cell histologies into CLL = 1, BCL = 2 (residual indolent B-cell-lymphoma pool including FL), DLBCL = 3, MCL = 4).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (in Gibiansky 2014 the implicit reference is CLL when paired with `TUMTP_DLBCL` and `TUMTP_MCL` all = 0; the residual indolent-B-cell-lymphoma pool in this paper is dominated by follicular lymphoma).
- **Source aliases:**
  - `DIS` (Gibiansky 2014; integer code with `DIS == 2` flagging BCL) -- decompose into `TUMTP_BCL = as.integer(DIS == 2)`.
- **Example models:** `Gibiansky_2014_obinutuzumab.R` (effect on time-dependent clearance decay rate kdes via the composite `(TUMTP_BCL + TUMTP_DLBCL + TUMTP_MCL)` (any-NHL effect; ratio 2.08) and on both time-dependent CL_T and steady-state CL_inf via the composite `(TUMTP_BCL + TUMTP_DLBCL)` (shared BCL/DLBCL effect; ratio 0.834 in the reverse direction, i.e., 16.6% lower CL than CLL)).
- **Notes:** Follows the `TUMTP_CHL` / `TUMTP_GC` / `TUMTP_SCLC` decomposition pattern. Distinct from `TUMTP_LYMPH` (a broader heterogeneous lymphoma pool that lumps cHL with NHL histologies) -- `TUMTP_BCL` is specifically B-cell lymphoma and pairs with sibling `TUMTP_DLBCL` and `TUMTP_MCL` for histology-specific contrasts within the NHL family. When a future paper studies follicular lymphoma in isolation (rather than pooled into BCL), register a more specific canonical (e.g., `TUMTP_FL`) rather than overloading this one. Ratified canonically on 2026-05-11.

### TUMTP_DLBCL (**canonical for diffuse large B-cell lymphoma indicator**)
- **Description:** 1 = diffuse large B-cell lymphoma (DLBCL), 0 = other tumor types. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (in Gibiansky 2014 the implicit reference is CLL when paired with `TUMTP_BCL` and `TUMTP_MCL` all = 0).
- **Source aliases:**
  - `DIS` (Gibiansky 2014; integer code with `DIS == 3` flagging DLBCL) -- decompose into `TUMTP_DLBCL = as.integer(DIS == 3)`.
- **Example models:** `Gibiansky_2014_obinutuzumab.R` (effect on kdes via the any-NHL composite indicator; effect on CL_T and CL_inf via the shared BCL/DLBCL composite indicator).
- **Notes:** Follows the `TUMTP_CHL` / `TUMTP_GC` / `TUMTP_SCLC` decomposition pattern. Distinct from `TUMTP_LYMPH` (broader lymphoma pool) and `TUMTP_PCALCL` (primary cutaneous anaplastic large-cell lymphoma; a CD30+ T-cell-lineage entity unrelated to DLBCL). DLBCL is the most common high-grade B-cell-NHL subtype; the Gibiansky 2014 cohort had only 30 DLBCL patients (4.4%), so a single estimated effect on CL is shared with BCL (a much larger pooled group; see `TUMTP_BCL` notes). Ratified canonically on 2026-05-11.

### TUMTP_MCL (**canonical for mantle cell lymphoma indicator**)
- **Description:** 1 = mantle cell lymphoma (MCL), 0 = other tumor types. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (in Gibiansky 2014 the implicit reference is CLL when paired with `TUMTP_BCL` and `TUMTP_DLBCL` all = 0).
- **Source aliases:**
  - `DIS` (Gibiansky 2014; integer code with `DIS == 4` flagging MCL) -- decompose into `TUMTP_MCL = as.integer(DIS == 4)`.
- **Example models:** `Gibiansky_2014_obinutuzumab.R` (effect on kdes via the any-NHL composite indicator; separate effect on CL_T and CL_inf via the standalone MCL indicator (ratio 1.75, i.e., 75% higher CL than CLL)).
- **Notes:** Follows the `TUMTP_CHL` / `TUMTP_GC` / `TUMTP_SCLC` decomposition pattern. Distinct from `TUMTP_LYMPH` (broader lymphoma pool). The Gibiansky 2014 cohort had only 20 MCL patients (2.9%); the paper reports the highest obinutuzumab CL among the four B-cell-malignancy histologies for MCL, consistent with the highest CD20 expression density on MCL B-cells relative to the other histologies. Ratified canonically on 2026-05-11.

### LINE_1L (**canonical for first-line-therapy indicator**)
- **Description:** 1 = first-line therapy (1L) / treatment-naive, 0 = second-line or greater (2L+) / relapsed-or-refractory.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (2L+, second-line or greater / relapsed-refractory).
- **Source aliases:**
  - `LINE` (categorical column with levels `1L`, `2L`, `3L+`, ...) -- decompose into `LINE_1L = as.integer(LINE == "1L")`.
  - `RRFN` (relapsed/refractory flag; treatment-naive corresponds to RRFN == 0) -- used in `Lu_2019_polatuzumab.R`. Decompose: `LINE_1L = as.integer(RRFN == 0)`.
- **Example models:** `Sanghavi_2020_ipilimumab.R` (exponential coefficient -0.0949 on CL), `Lu_2019_polatuzumab.R` (multiplicative effects on V1 = 1.20, kdes = 3.38, CL_T = 3.53, FRAC_NS = 0.756; the same pooled-trial NHL cohort mixes 415 R/R and 45 first-line patients).
- **Notes:** Promoted to scope: general on 2026-04-26 after Lu 2019 polatuzumab vedotin ratified the same 1L vs 2L+ binarization that Sanghavi 2020 ipilimumab introduced. The two papers use different indicator semantics (Sanghavi reports the effect as `exp(-0.0949 * LINE_1L)` and Lu reports `theta^LINE_1L` with theta < or > 1 depending on the parameter); both reduce to the same canonical 0/1 column. If a future paper requires finer resolution (separate effects for 2L vs 3L+), add a parallel `LINE_2L` canonical rather than overloading this one.

### NIVO_1Q3W (**canonical for nivolumab 1 mg/kg every 3 weeks co-administration indicator**)
- **Description:** 1 = ipilimumab co-administered with nivolumab 1 mg/kg every 3 weeks; 0 = otherwise (monotherapy or any other nivolumab regimen).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no nivolumab or any non-1Q3W nivolumab regimen).
- **Source aliases:**
  - `NIVO_REGIMEN` (categorical column with levels `none`, `0.3 mg/kg Q3W`, `1 mg/kg Q2W`, `1 mg/kg Q3W`, `3 mg/kg Q2W`, `3 mg/kg Q3W`) -- decompose into `NIVO_1Q3W = as.integer(NIVO_REGIMEN == "1 mg/kg Q3W")`.
- **Example models:** `Sanghavi_2020_ipilimumab.R` (exponential coefficient 0.0950 on ipilimumab CL).
- **Notes:** Paired with `NIVO_3Q2W` in the Sanghavi 2020 final model; both decomposed indicators are 0 for ipilimumab monotherapy. Other nivolumab regimens (0.3 mg/kg Q3W, 1 mg/kg Q2W, 3 mg/kg Q3W) were tested but not retained in the final model and collapse into the reference 0 group.

### NIVO_3Q2W (**canonical for nivolumab 3 mg/kg every 2 weeks co-administration indicator**)
- **Description:** 1 = ipilimumab co-administered with nivolumab 3 mg/kg every 2 weeks; 0 = otherwise (monotherapy or any other nivolumab regimen).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no nivolumab or any non-3Q2W nivolumab regimen).
- **Source aliases:**
  - `NIVO_REGIMEN` (categorical column) -- decompose into `NIVO_3Q2W = as.integer(NIVO_REGIMEN == "3 mg/kg Q2W")`.
- **Example models:** `Sanghavi_2020_ipilimumab.R` (exponential coefficient 0.191 on ipilimumab CL).
- **Notes:** Paired with `NIVO_1Q3W`; same reference grouping convention.

### COMBO_NIVO (**canonical for any-regimen nivolumab combination-therapy indicator**)
- **Description:** 1 = ipilimumab co-administered with any nivolumab regimen, 0 = ipilimumab monotherapy.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (ipilimumab monotherapy).
- **Source aliases:**
  - `COMBO` -- used in `Sanghavi_2020_ipilimumab.R`. Equivalently derivable from `NIVO_REGIMEN` as `COMBO_NIVO = as.integer(NIVO_REGIMEN != "none")`.
- **Example models:** `Sanghavi_2020_ipilimumab.R` (additive effect -0.202 on the Emax parameter of the time-varying CL function).
- **Notes:** Distinct from the per-regimen `NIVO_1Q3W` / `NIVO_3Q2W` indicators on baseline CL: `COMBO_NIVO` aggregates across all nivolumab regimens and acts on the time-varying-CL Emax parameter, whereas the per-regimen indicators act on baseline (time-zero) CL.

### BLSTABL (**canonical for baseline absolute blast counts in peripheral blood**)
- **Description:** Baseline absolute count of blasts (immature lymphoid/myeloid precursor cells) circulating in peripheral blood. Time-fixed baseline value.
- **Units:** 10^9 counts/L (equivalently 10^9 counts; reported as "x 10^9 counts" in Wu 2024 Table 2).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with power scaling `(BLSTABL / <ref>)^exponent`. Reference value observed: 0.352 x 10^9 counts (Wu 2024 Table 3, BCP-ALL median).
- **Source aliases:**
  - `BLSTABL` -- used in `Wu_2024_inotuzumab.R`.
- **Example models:** `Wu_2024_inotuzumab.R` (power exponent -0.0484 on kdes for BCP-ALL patients only; the effect is gated off for B-cell NHL by multiplying the exponent by `DIS_BCPALL`).
- **Notes:** Distinct from blasts in bone marrow (different specimen) and from `BLSTPB` (percentage of blasts in peripheral blood, used by the predecessor Garrett 2019 adult model). Not applicable for B-cell NHL patients in pooled BCP-ALL + NHL analyses (Wu 2024 retains the effect only in BCP-ALL patients via the DIS_BCPALL gate). When supplying BLSTABL for an NHL subject, set the value to the BCP-ALL reference (0.352) so the gated power term evaluates to 1 numerically. Scope: specific because the covariate is most meaningful in B-cell-leukemia population PK analyses; promote to general if a second paper retains it.

### COMBO_RG (**canonical for anti-CD20 (rituximab or obinutuzumab) combination-therapy indicator**)
- **Description:** 1 = polatuzumab vedotin co-administered with rituximab OR obinutuzumab, 0 = single-agent polatuzumab vedotin (or any other regimen lacking an anti-CD20 partner). Time-fixed per subject in the source paper's analysis cohort.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (single-agent polatuzumab vedotin or no anti-CD20 partner).
- **Source aliases:**
  - `COMBO` (categorical: 0 = single agent, 1 = + rituximab, 2 = + obinutuzumab) -- used in `Lu_2019_polatuzumab.R`. Decompose: `COMBO_RG = as.integer(COMBO == 1 | COMBO == 2)`. The Lu 2019 NONMEM separately defines `RTX = as.integer(COMBO == 1)` and `GA101 = as.integer(COMBO == 2)` and applies effects as `theta^(RTX + GA101)`; because RTX and GA101 are mutually exclusive, RTX + GA101 takes values {0, 1} and the effect collapses to `theta^COMBO_RG`.
- **Example models:** `Lu_2019_polatuzumab.R` (multiplicative effects on CL_INF = 0.844, kdes = 0.932, FRAC_NS = 0.709).
- **Notes:** Rituximab and obinutuzumab both bind CD20 on B cells (rituximab is a Type I anti-CD20 mAb, obinutuzumab a glycoengineered Type II), so co-administration is hypothesized to alter polatuzumab vedotin disposition through depletion of CD79b+ target B cells. The Lu 2019 final model fits a single combined effect rather than separate rituximab- and obinutuzumab-specific effects. Scope: specific because the relevant combination partners (CD20-directed mAbs) are tied to NHL pathway; if a future paper distinguishes rituximab from obinutuzumab combinations, register `COMBO_R` and `COMBO_G` separately rather than overloading this canonical.

### COMBO_DURVA (**canonical for durvalumab combination-therapy indicator**)
- **Description:** 1 = the analyzed therapeutic mAb is co-administered with durvalumab (anti-PD-L1 IgG1), 0 = monotherapy.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (monotherapy).
- **Source aliases:**
  - `COMB` -- used in `Hwang_2022_tremelimumab.R` ($INPUT NM-TRAN data item; control-stream switch `IF(COMB.EQ.0)` selects monotherapy parameters and `IF(COMB.EQ.1)` selects combination-therapy parameters).
- **Example models:** `Hwang_2022_tremelimumab.R` (selects between monotherapy and combination-with-durvalumab values of the time-varying-CL Tmax and lambda parameters).
- **Notes:** Parallels `COMBO_NIVO` but for durvalumab rather than nivolumab co-administration. Acts on the time-varying-CL component (Tmax and lambda); baseline CL is shared between monotherapy and combination groups in Hwang 2022.

### COMBO_LEN_DEX (**canonical for lenalidomide plus dexamethasone combination-therapy indicator**)
- **Description:** 1 = the analyzed therapeutic mAb (or other agent under PK study) is co-administered with the lenalidomide + low-dose-dexamethasone (Ld) backbone, 0 = monotherapy or any non-Ld regimen. Lenalidomide is an immunomodulatory imide (IMiD) that activates natural killer cells; dexamethasone is an immunosuppressant glucocorticoid. The Ld backbone is a standard combination partner in multiple-myeloma and other hematologic-malignancy regimens.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (monotherapy or any non-Ld regimen). When a paper reports the reference patient as "with Ld coadministration" (as Ide 2020 does), the model still stores the canonical 0/1 column and applies the effect as `exp(theta * (COMBO_LEN_DEX - 1))` so that COMBO_LEN_DEX = 1 yields factor 1 (paper's reference) and COMBO_LEN_DEX = 0 activates the effect.
- **Source aliases:**
  - `LENDEX` (1 = with Ld, 0 = without Ld; in Ide 2020 derived from `STUDY != 204011` because study 204011 was the Ld-free elotuzumab-monotherapy cohort) -- used in `Ide_2020_elotuzumab.R`.
  - `COMBO_LD` (**retired** canonical name; renamed to `COMBO_LEN_DEX` on 2026-04-27 for clarity).
- **Example models:** `Ide_2020_elotuzumab.R` (multiplicative effects: CLLd = 0.74 on nonspecific CL, encoded as `exp(log(0.74) * (COMBO_LEN_DEX - 1))`; KINTLd = 10.1 on the second-order target-mediated elimination rate from the peripheral compartment, encoded as `exp(log(10.1) * (COMBO_LEN_DEX - 1))`).
- **Notes:** Specific scope because the canonical's mechanistic relevance is hematologic-malignancy-domain-bound (multiple myeloma and related plasma-cell or B-cell disorders). Distinct from `COMBO_BELAMAF` (which pools Ld with bortezomib-dex and pomalidomide-dex into a single broader "any-combination" belantamab indicator); `COMBO_LEN_DEX` is the per-backbone Ld-only flag. If a future paper distinguishes "Ld-only" from a broader "any-IMiD-plus-dex" backbone with separate effects, register a parallel canonical (e.g., `COMBO_PD` for pomalidomide-dex, `COMBO_VD` for bortezomib-dex). Sign of the exponential coefficient is paper-dependent: Ide 2020 reports `CLLd = 0.74` so the Ld-coadministration arm has 26% lower nonspecific CL than the Ld-free arm, but `KINTLd = 10.1` so the Ld arm has 10x higher second-order target-mediated elimination -- both are mechanistically interpretable (dexamethasone suppresses non-specific catabolic clearance; lenalidomide-activated NK cells increase target-cell-binding-mediated elimination).

### COMBO_BELAMAF (**canonical for any-combination belantamab mafodotin therapy indicator**)
- **Description:** 1 = belantamab mafodotin administered as part of a combination regimen (with bortezomib + dexamethasone, lenalidomide + dexamethasone, or pomalidomide + dexamethasone) in the relapsed/refractory multiple myeloma setting; 0 = belantamab mafodotin monotherapy.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (belantamab mafodotin monotherapy).
- **Source aliases:**
  - `COMBO` (when the source dataset uses a generic combination flag for belantamab mafodotin pooled regimens) -- used in `Papathanasiou_2025_belantamab.R`.
- **Example models:** `Papathanasiou_2025_belantamab.R` (multiplicative factor theta = 1.44 on the Imax parameter of the time-varying CL function -- combination therapy increases the steady-state CL reduction from 33.2 % to 44.0 %).
- **Notes:** Pools the three combination backbones tested in DREAMM-6 / DREAMM-7 / DREAMM-8 (Bor-Dex, Len-Dex, Pom-Dex) into a single binary because Papathanasiou 2025 reports no meaningful per-backbone difference in cycle-1 ADC exposure. If a future paper tests per-backbone combination effects, register dedicated indicators (`COMBO_BELAMAF_BORDEX`, etc.) rather than overloading this aggregate.

### KG (**canonical for subject-specific tumour-growth rate constant from a prior IPP fit**)
- **Description:** Empirical-Bayes posterior estimate of the subject-specific tumour-size first-order growth rate constant carried over from an upstream tumour-size-dynamics population model. Supplied per subject in the dataset and used directly inside a downstream model (e.g., overall-survival hazard) that integrates the tumour-size ODE inline conditional on each subject's growth/death rates.
- **Units:** internal scaled rate; the Zecchin 2016 OS model carries the source convention `KG / 1000 * tumorSize` in the SLD ODE so the column units are `(1/day) * 1000` as published in the source NONMEM run. Document per-model via `covariateData[[KG]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used directly inside the SLD ODE (see `Zecchin_2016_survival.R`).
- **Source aliases:**
  - `KG` (NONMEM `$INPUT` column in DDMODEL00000218; identical column shipped in the bundle's Simulated_OS.csv).
- **Example models:** `Zecchin_2016_survival.R` (Zecchin 2016 OS model, DDMODEL00000218).
- **Notes:** Specific scope because the column is the empirical-Bayes output of a particular upstream fit (`modellib('Zecchin_2016_tumorovarian')` / DDMODEL00000217). When this OS model is used standalone, the user must supply `KG` per subject -- typically by first fitting the SLD model and extracting the per-subject empirical-Bayes posterior. The internal `/1000` scaling is preserved verbatim from the source `$DES` block to maintain numerical equivalence with the published estimates.

### KD0 (**canonical for subject-specific carboplatin-related tumour-death rate constant from a prior IPP fit**)
- **Description:** Empirical-Bayes posterior estimate of the subject-specific carboplatin-driven tumour-size death rate constant carried over from an upstream tumour-size-dynamics population model. Pairs with the time-varying `AUC_CARBO` covariate inside the SLD ODE term `KD0 * AUC_CARBO * tumorSize`.
- **Units:** internal scaled rate; the Zecchin 2016 OS model carries the source convention `KD0 / 1000 * AUC_CARBO * tumorSize` so the column units are `(1/day per AUC unit) * 1000` as published in the source NONMEM run. Document per-model via `covariateData[[KD0]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used directly inside the SLD ODE (see `Zecchin_2016_survival.R`).
- **Source aliases:**
  - `KD0` (NONMEM `$INPUT` column in DDMODEL00000218).
- **Example models:** `Zecchin_2016_survival.R` (Zecchin 2016 OS model, DDMODEL00000218).
- **Notes:** Specific scope because the column is tied to a specific drug (carboplatin) and a specific upstream IPP fit. The internal `/1000` scaling is preserved verbatim from the source `$DES` block.

### KD1 (**canonical for subject-specific gemcitabine-related tumour-death rate constant from a prior IPP fit**)
- **Description:** Empirical-Bayes posterior estimate of the subject-specific gemcitabine-driven tumour-size death rate constant carried over from an upstream tumour-size-dynamics population model. Pairs with the time-varying `AUC_GEM` covariate inside the SLD ODE term `KD1 * AUC_GEM * tumorSize`.
- **Units:** internal scaled rate; the Zecchin 2016 OS model carries the source convention `KD1 / 100 * AUC_GEM * tumorSize` so the column units are `(1/day per AUC unit) * 100` as published in the source NONMEM run. Document per-model via `covariateData[[KD1]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used directly inside the SLD ODE (see `Zecchin_2016_survival.R`).
- **Source aliases:**
  - `KD1` (NONMEM `$INPUT` column in DDMODEL00000218).
- **Example models:** `Zecchin_2016_survival.R` (Zecchin 2016 OS model, DDMODEL00000218).
- **Notes:** Specific scope because the column is tied to a specific drug (gemcitabine) and a specific upstream IPP fit. The internal `/100` scaling is preserved verbatim from the source `$DES` block.

### IBASE (**canonical for subject-specific baseline tumour-size estimate from a prior IPP fit**)
- **Description:** Empirical-Bayes posterior estimate of the subject-specific baseline sum-of-longest-diameters (SLD) tumour size carried over from an upstream tumour-size-dynamics population model. Used both to set the SLD ODE initial state (`tumorSize(0) <- IBASE * 1000` in the Zecchin 2016 model: source convention multiplies by 1000 to convert the internal value to mm) and to scale the time-varying tumour-size ratio (`mmbas <- IBASE * 1000`) inside the OS hazard.
- **Units:** internal scaled length; the Zecchin 2016 OS model carries the source convention `IBASE * 1000 = mm` so the column itself is in metres (1 m = 1000 mm). Document per-model via `covariateData[[IBASE]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used as the SLD ODE initial condition and as the denominator in `(tumorSize - mmbas) / mmbas`.
- **Source aliases:**
  - `IBASE` (NONMEM `$INPUT` column in DDMODEL00000218; identical column shipped in the bundle's Simulated_OS.csv with values typically in the 0.04-0.50 m range).
- **Example models:** `Zecchin_2016_survival.R` (Zecchin 2016 OS model, DDMODEL00000218).
- **Notes:** Distinct from the canonical `TUM_SLD` column. `TUM_SLD` carries the *measured* baseline tumour size in mm (used to compute the time-fixed `NSLD0 = TUM_SLD / 70` covariate term in the Zecchin 2016 OS hazard), whereas `IBASE` carries the empirical-Bayes *fitted* baseline from the upstream SLD model (used to initialise the integrated SLD trajectory and to define the time-varying TSR(t) reference). The two are correlated but not equal because the upstream IPP fit smooths measurement noise away from the observed SLD0. Specific scope because the column is the empirical-Bayes output of a specific upstream model fit and the internal `*1000` unit-conversion is tied to the source NONMEM coding convention.

### NWLS (**canonical for time-varying new-lesion appearance indicator**)
- **Description:** Time-varying binary indicator of whether a new (non-target) RECIST lesion has appeared since enrolment. 1 = new lesion present at the current observation time; 0 = no new lesion as of the current observation time. Once `NWLS` flips to 1 it stays 1 for subsequent observation times in that subject (a step-function flag, not a transient pulse).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no new lesion appeared as of the current time).
- **Source aliases:**
  - `NWLS` (NONMEM `$INPUT` column in DDMODEL00000218); the bundle's Simulated_OS.csv re-labels the same column `NWLSCOV`. Downstream consumers should map `NWLSCOV` -> `NWLS`.
- **Example models:** `Zecchin_2016_survival.R` (Zecchin 2016 OS model, DDMODEL00000218; multiplicative effect on the Weibull hazard via `exp(e_nwls_haz * NWLS)` with `e_nwls_haz = 1.23` per Output_real_OS.lst FINAL TH5 / Table 2 of Zecchin 2016).
- **Notes:** Specific scope because the column encodes a per-paper RECIST-style binary that is supplied by the dataset; it is not a generic "any new lesion" indicator the user can populate from routine clinical data without an explicit lesion-appearance imaging schedule. The Zecchin 2016 OS model uses `NWLS` directly (no time-gating), which is faithful to the simulated dataset shipped in the DDMORE bundle. The source `Output_real_OS.lst` (the listing on the original real dataset) gates the indicator with an additional `TNWLS` (lesion-appearance-time) column not shipped in the bundle's simulated dataset; the two encodings are functionally equivalent when the dataset's `NWLS` column is constructed as a 0/1 step that flips at the lesion-appearance time. The bundle's simulated dataset uses the simpler step-function form, and that is the form the nlmixr2lib model expects.

### MM_NIGG (**canonical for non-IgG multiple myeloma immunoglobulin-type indicator**)
- **Description:** 1 = patient with non-IgG-secreting multiple myeloma (e.g., IgA, IgD, IgE, IgM, light-chain-only / Bence Jones, or non-secretory MM), 0 = patient with IgG-secreting multiple myeloma.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (IgG MM).
- **Source aliases:**
  - `Ig_type` -- used in `Fau_2020_isatuximab.R`. Values 0 / 1 with the same orientation as the canonical (1 = non-IgG MM).
- **Example models:** `Fau_2020_isatuximab.R` (exponential effect on the steady-state linear CL CLinf with coefficient -0.751, and on the time-varying-CL half-time KCL with coefficient -0.931), `Xu_2020_daratumumab.R` (additive shift `(1 + 0.806 * (1 - MM_NIGG))` on linear CL — Xu 2020 parameterises with non-IgG MM as reference, so an IgG MM patient receives an 80.6% higher linear CL than a non-IgG MM patient; canonical column semantics 1 = non-IgG / 0 = IgG are preserved).
- **Notes:** Within-disease (multiple-myeloma) immunoglobulin-subtype stratifier. The mechanistic rationale (Fau 2020) is that endogenous IgG monoclonal protein in IgG-MM patients competes with the therapeutic IgG mAb for FcRn-mediated salvage, raising the therapeutic mAb's catabolic clearance; non-IgG-MM patients lack that competition and exhibit lower therapeutic-mAb clearance. Distinct from the disease-state indicators (`DIS_SMM` = smoldering MM); applies only after a multiple-myeloma diagnosis is established. Scope: specific because the comparison is a within-MM stratifier rather than a cross-population indicator. Reference category at the model level (which value of MM_NIGG corresponds to TVCL = base) varies between papers: Fau 2020 anchors to 0 (IgG MM) and Xu 2020 anchors to 1 (non-IgG MM); the canonical column orientation (1 = non-IgG) is fixed across papers and the per-model `covariateData[[MM_NIGG]]$reference_category` field records which anchor each model uses.

### TUM_TP53_MUT (**canonical for tumour TP53 / p53 mutation indicator**)
- **Description:** Binary indicator of tumour-cell TP53 mutational status as assessed in the source paper. 1 = TP53 mutant tumour (typically a missense mutation detected by sequencing, or p53 protein overexpression by immunohistochemistry used as a surrogate for TP53 missense mutation per Gillet et al. J Neurooncol 2014); 0 = TP53 wild-type tumour. Time-fixed per subject (a somatic tumour-genotype call made at diagnosis, not a germline genotype).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (TP53 wild-type).
- **Source aliases:**
  - "p53 mutation" / "p53 mutated" (Mazzocco 2015, paper text; p53 overexpression by IHC used as a surrogate for TP53 missense mutation).
- **Example models:** `Mazzocco_2015_temozolomide.R` (exponential effect on the temozolomide tumour-cell-death rate constant `gamma`: `gamma = gamma0 * exp(beta_p53 * TUM_TP53_MUT)` with `beta_p53 = log(0.143 / 0.254) = -0.574`; TP53-mutant LGG tumours are 44% less sensitive to TMZ than TP53-wild-type tumours).
- **Notes:** Specific scope because the column encodes a somatic tumour-genotype call (mechanism = altered DNA-damage response in tumour cells) rather than a germline pharmacogenomic variant. Distinct from the `SNP_<GENE>_<RSID>` family, which encodes inherited host germline genotypes affecting drug PK; `TUM_TP53_MUT` encodes a tumour-cell mutation and only makes mechanistic sense for drugs whose effect depends on a functional p53 pathway in the tumour. The Mazzocco 2015 cohort assayed p53 status by IHC overexpression; future extractions that use direct TP53 sequencing should still record their values under this canonical and document the assay method in `covariateData[[TUM_TP53_MUT]]$notes`. Ratified canonically on 2026-05-17 alongside the Mazzocco 2015 temozolomide extraction.

### TUM_1P19Q_CODEL (**canonical for tumour 1p/19q chromosomal codeletion indicator**)
- **Description:** Binary indicator of tumour-cell 1p/19q chromosomal codeletion status. 1 = tumour carries the combined loss of the short arm of chromosome 1 (1p) and the long arm of chromosome 19 (19q); 0 = non-codeleted tumour (intact 1p and/or 19q). Time-fixed per subject (a somatic tumour-genotype call made at diagnosis, typically by fluorescence in-situ hybridisation, comparative genomic hybridisation, or loss-of-heterozygosity assay).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-codeleted; intact 1p and/or 19q).
- **Source aliases:**
  - "1p/19q codeletion" / "1p/19q codeleted" (Mazzocco 2015, paper text).
- **Example models:** `Mazzocco_2015_temozolomide.R` (exponential effect on the damaged-quiescent-to-proliferative repair rate constant `kQpP`: `kQpP = kQpP0 * exp(beta_1p19q * TUM_1P19Q_CODEL)` with `beta_1p19q = log(0.00807 / 0.00947) = -0.160`; 1p/19q-codeleted LGG tumours have 15% lower kQpP than non-codeleted tumours, consistent with longer reported duration of response in codeleted patients).
- **Notes:** Specific scope because the column encodes a brain-tumour-specific somatic chromosomal alteration with mechanistic relevance to DNA-repair capacity in glioma cells; it is not a generic "any tumour-chromosomal-alteration" indicator. The 1p/19q codeletion is a defining molecular feature of oligodendrogliomas (per the 2016 WHO classification of CNS tumours) and is mutually exclusive with TP53 missense mutation in the Mazzocco 2015 cohort (Ricard 2007, ref 12 of Mazzocco 2015). Ratified canonically on 2026-05-17 alongside the Mazzocco 2015 temozolomide extraction.

## Laboratory / disease-activity

### ALBR (**canonical for serum albumin normalized to the laboratory's upper limit of normal**)
- **Description:** Serum albumin normalized to the laboratory's upper limit of normal (`albumin_observed / ULN_albumin`).
- **Units:** (unitless ratio)
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used as a power term `(ALBR / <ref>)^exponent`. Reference 0.78 used in Xu 2019 (corresponds to a median serum albumin of 38 g/L at a typical ULN of ~48.7 g/L).
- **Source aliases:** none.
- **Example models:** `Xu_2019_sarilumab.R`.
- **Notes:** Xu 2019 normalizes to each site's ULN so that values across multiple labs with different reference ranges can be pooled. Scoped specific because the ULN-normalization convention is tied to the Xu 2019 analysis plan; future papers using the same ratio should either add themselves to the example_models list or promote this entry to general.

## Inflammatory-bowel-disease disease-activity covariates

### CALPRO (**canonical for fecal calprotectin**)
- **Description:** Fecal calprotectin, a gut-inflammation biomarker (baseline or time-fixed per subject unless a paper explicitly uses a time-varying value).
- **Units:** mg/kg stool (equivalent to ug/g). Document per-model via `covariateData[[CALPRO]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(CALPRO / ref)^exponent`. Reference 700 mg/kg used in Rosario 2015 (overall population median).
- **Source aliases:** none.
- **Example models:** `Rosario_2015_vedolizumab.R` (reference 700 mg/kg; exponent +0.0310 on linear clearance CLL).
- **Notes:** Common IBD severity biomarker (inflammation of the gut epithelium). Assays typically report in ug/g stool; 1 ug/g = 1 mg/kg. Document per-model whether baseline-only or time-varying in `covariateData[[CALPRO]]$notes`.

### CDAI (**canonical for Crohn's Disease Activity Index**)
- **Description:** Crohn's Disease Activity Index composite score. Higher values indicate more active disease; <150 remission, 150-219 mild, 220-450 moderate, >450 severe. Defined only for patients with a CD diagnosis; set to the reference value (or gate via the indicator) for UC patients.
- **Units:** (score, 0-600)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(CDAI / ref)^exponent`. Reference 300 used in Rosario 2015 (typical moderate-CD score).
- **Source aliases:** none.
- **Example models:** `Rosario_2015_vedolizumab.R` (reference 300; exponent -0.0515 on CLL gated by `IBD_CD` so the effect applies only to CD patients).
- **Notes:** Mutually exclusive with `PMAYO` in pooled UC+CD populations: each patient has exactly one disease-activity score (CDAI for CD, partial Mayo for UC). Gate via the `IBD_CD` indicator when pooling.

### PMAYO (**canonical for partial Mayo score**)
- **Description:** Partial Mayo score for ulcerative colitis (sum of stool-frequency, rectal-bleeding, and physician-global subscores, range 0-9). Higher values indicate more active disease. Defined only for patients with a UC diagnosis; gate via the `IBD_CD` indicator when pooling UC+CD.
- **Units:** (score, 0-9)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(PMAYO / ref)^exponent`. Reference 6 used in Rosario 2015 (typical moderate-UC score).
- **Source aliases:** none.
- **Example models:** `Rosario_2015_vedolizumab.R` (reference 6; exponent +0.0408 on CLL gated by `(1 - IBD_CD)` so the effect applies only to UC patients).
- **Notes:** The partial Mayo score excludes the endoscopy subscore (the full Mayo score is 0-12). Mutually exclusive with `CDAI` in pooled UC+CD populations.

### MAYO_E (**canonical for baseline Mayo endoscopic subscore**)
- **Description:** Mayo endoscopic subscore at baseline for ulcerative colitis, integer 0-3 (0 = normal / inactive disease, 1 = mild, 2 = moderate, 3 = severe). The endoscopic subscore is one of the four components of the full Mayo score (0-12); the partial Mayo score (`PMAYO`) excludes it. Time-fixed per subject.
- **Units:** (score, 0-3)
- **Type:** categorical
- **Scope:** general
- **Reference category:** depends on per-model encoding -- papers that include the full 0-3 range typically reference Mayo 0 or 1 (mild), while papers restricted to moderate-to-severe UC (the typical biologic-induction-therapy population) reference Mayo 2 or 3. Document the per-model reference category in `covariateData[[MAYO_E]]$reference_category`.
- **Source aliases:**
  - `MPRE` -- used in `Faelens_2021_infliximab.R` (NONMEM column for "Mayo endoscopic score pre-induction"). The Faelens 2021 dataset additionally codes a "missing" sentinel `MPRE = -99`; treat this as out-of-domain when applying the model and document per-model.
- **Example models:** `Faelens_2021_infliximab.R` (categorical effect on KE: separate typical KE for Mayo 1 / Mayo 2 / Mayo 3; reference category Mayo 2).
- **Notes:** Distinct from the full Mayo score (0-12) and the partial Mayo score (`PMAYO`, 0-9). The endoscopic subscore alone is the core inclusion criterion in many UC induction-therapy popPK datasets (typically Mayo 2 or 3 = moderate-to-severe disease). Mutually compatible with `PMAYO` in datasets that report both.

### ENDO_ULCER (**canonical for endoscopically active luminal disease at baseline**)
- **Description:** 1 = mucosal ulcerations confirmed at baseline ileocolonoscopy / endoscopy (endoscopically active luminal disease), 0 = no mucosal ulcerations at baseline. Time-fixed (assessed at study entry).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no mucosal ulcerations at baseline).
- **Source aliases:** none known.
- **Example models:** `Aguiar_2021_ustekinumab.R` (Aguiar 2021 Table 3; covariate on baseline fecal calprotectin FC0: 213 mg/kg with ulcers vs 102 mg/kg without).
- **Notes:** IBD-specific structural-disease-activity indicator, distinct from the symptom-driven CDAI / PMAYO and the biomarker-driven CALPRO / CRP. Useful as a covariate on baseline biomarkers (FC, CRP) and as a stratifier for simulations of biochemical remission.

## Inflammatory-bowel-disease diagnosis

### IBD_CD (**canonical for Crohn's disease indicator**)
- **Description:** 1 = Crohn's disease (CD) diagnosis, 0 = ulcerative colitis (UC) diagnosis. Used to gate pooled UC+CD population models where some parameters or covariate effects differ by diagnosis.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (ulcerative colitis).
- **Source aliases:**
  - `DIAGNOSIS`, `DX` (categorical `"UC"`/`"CD"`) -- derive `IBD_CD = as.integer(DX == "CD")`.
- **Example models:** `Rosario_2015_vedolizumab.R` (two typical-CLL switch between UC and CD; gates `PMAYO` and `CDAI` effects; multiplicative +1% effect on Vc).
- **Notes:** Rosario 2015 models separate typical CLL for UC vs CD and gates the partial-Mayo (UC-only) and CDAI (CD-only) disease-activity covariates via this indicator.

## Concomitant / prior medication

### CONMED_AD (**canonical for concomitant Alzheimer's-symptomatic medication**)
- **Description:** 1 = subject is on concomitant Alzheimer's-symptomatic medication (typically a cholinesterase inhibitor and/or memantine) at the observation, 0 = not on such medication. Time-fixed in the Conrado 2014 source dataset (one indicator per subject) and treated as such here; future models that need a time-varying form may register a successor canonical or document time-varying use in `covariateData[[CONMED_AD]]$notes`.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 1 (on concomitant Alzheimer's-symptomatic medication) in the Conrado 2014 model -- the source paper labels CONMED_AD = 1 as the "most common" category and centres the slope-effect coefficient on this group, so the multiplicative factor on disease-progression slope is `1` for `CONMED_AD = 1` and `1 + e_sl_conmed_ad_off` for `CONMED_AD = 0`. The non-standard "most-common-as-reference" convention is preserved from the source for traceability; future Alzheimer's models that adopt the more common "off-treatment as reference" convention should register a successor canonical (`CONMED_AD_TREATED`) rather than overloading this one.
- **Source aliases:**
  - `COMED2` -- used in `Conrado_2014_alzheimer.R` (DDMORE Foundation Model Repository entry DDMODEL00000290). The `2` suffix in the source distinguishes this binary from upstream `COMED` and `PRIMCOMED` columns (free-text concomitant-medication entries) in the same NONMEM input dataset; the binary `COMED2` flag is what enters the model.
- **Example models:** `Conrado_2014_alzheimer.R` (multiplicative factor on the typical-value disease-progression slope: `slope_factor = 1 + (1 - CONMED_AD) * e_sl_conmed_ad_off` with `e_sl_conmed_ad_off = -0.302`, i.e. ~30% smaller progression slope for the off-treatment reference cohort).
- **Notes:** The source `.mod` does not specify the symptomatic-medication class beyond the binary flag; the Conrado 2014 publication context (CAMD ADAS-Cog disease-progression dataset, 2014) makes cholinesterase-inhibitor / memantine the dominant interpretation. Treating `CONMED_AD = 1` as the reference category is unusual relative to the rest of the `CONMED_*` family (`CONMED_PARA`, `CONMED_NSAID`, `CONMED_AZA`, etc.) which all use 0 = not-on as reference; the inversion is preserved here only because the source paper's coefficient was estimated with the "most common = 1" convention. Ratified canonically on 2026-05-06 alongside the Conrado 2014 DDMORE extraction.

### CONMED_AED (**canonical for any concomitant antiepileptic drug coadministration**)
- **Description:** 1 = subject is on at least one concomitant antiepileptic drug (AED) other than the modelled AED at the PK observation, 0 = on the modelled AED as monotherapy. Time-varying when concurrent AEDs cycle on / off across study occasions, time-fixed when the source paper analyses chronic-maintenance cohorts whose concurrent therapy is unchanged across the analysis window.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (modelled-AED monotherapy).
- **Source aliases:**
  - `CO_AED` -- used in `Yukawa_1990_phenytoin.R` (paper's `CO` indicator inverted: source `CO = 1` if PHT alone, theta_co otherwise; canonical `CONMED_AED = 1 - CO_indicator` so 0 is the PHT-monotherapy reference).
- **Example models:** `Yukawa_1990_phenytoin.R` (multiplicative `^CONMED_AED` factor on Vmax (1.08) and Km (1.32) for chronic phenytoin patients on at least one of phenobarbital, carbamazepine, valproate, primidone, clonazepam, sultiame, ethotoin, ethosuximide, acetazolamide, or diazepam).
- **Notes:** Generic concomitant-AED indicator covering the heterogeneous mix of older AEDs (PB, CBZ, VPA, primidone, clonazepam, etc.) studied alongside the modelled drug; the per-paper list of qualifying AEDs must be documented in `covariateData[[CONMED_AED]]$notes`. Distinct from drug-specific concomitant-AED indicators (e.g., a future `CONMED_PB` for concomitant phenobarbital alone) which would warrant separate canonicals when a paper distinguishes effects by AED class. Follows the `CONMED_*` family pattern (`CONMED_AZA`, `CONMED_NSAID`, etc.). Ratified canonically on 2026-05-10 alongside the Yukawa 1990 phenytoin extraction.

### CONMED_AMIO (**canonical for concomitant amiodarone coadministration indicator**)
- **Description:** 1 = subject is coadministered amiodarone (Class III antiarrhythmic; CYP3A4 / CYP2C9 / P-gp inhibitor) during the study, 0 = no concomitant amiodarone. Time-varying when amiodarone start / stop events are captured; the Xia 2024 source treats amiodarone as time-fixed at the analysis baseline.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant amiodarone).
- **Source aliases:**
  - `CM1` -- used in `Xia_2024_warfarin.R` (Xia 2024 Figure 1 / Table 2 covariate-screening variable labelling: `CM1 = 1` indicates combined amiodarone, `0` no combination).
- **Example models:** `Xia_2024_warfarin.R` (piecewise multiplicative effect on warfarin EC50: `ec50 *= (1 + e_amio_ec50 * CONMED_AMIO)` with `e_amio_ec50 = -0.602`, i.e. amiodarone reduces EC50 by ~60% in the Han Chinese cohort).
- **Notes:** Amiodarone is the canonical CYP2C9 / CYP3A4 inhibitor that potentiates warfarin's anticoagulant effect in clinical practice; the Xia 2024 cohort prevalence was 25.7% (Table 1). The per-paper definition (any amiodarone use vs current loading-dose use vs steady-state use) should be documented in `covariateData[[CONMED_AMIO]]$notes`. Ratified canonically on 2026-05-16 alongside the Xia 2024 warfarin extraction.

### CONMED_AMINO (**canonical for concomitant aminosalicylate therapy**)
- **Description:** 1 = on concomitant aminosalicylate (5-aminosalicylic acid / mesalamine / mesalazine / olsalazine / sulfasalazine etc.) therapy at the PK observation, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant aminosalicylate).
- **Source aliases:** `AMINO` -- used in `Rosario_2015_vedolizumab.R`.
- **Example models:** `Rosario_2015_vedolizumab.R` (power-form on CLL: `CLL * 1.02^CONMED_AMINO`).
- **Notes:** Covers the full aminosalicylate class (5-ASA is the single active moiety shared by most agents); use `CONMED_AMINO` rather than `CONMED_5ASA` unless the source paper explicitly restricts the indicator to 5-ASA monotherapy.

### CONMED_AVD (**canonical for brentuximab vedotin + AVD (adriamycin/doxorubicin, vinblastine, dacarbazine) combination indicator**)
- **Description:** 1 = subject is receiving brentuximab vedotin in combination with the AVD chemotherapy backbone (adriamycin a.k.a. doxorubicin, vinblastine, dacarbazine) for newly diagnosed advanced-stage Hodgkin lymphoma; 0 = otherwise (single-agent brentuximab vedotin). Encodes the A+AVD frontline regimen as a study-design covariate on ADC clearance.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (single-agent brentuximab vedotin -- no AVD coadministration).
- **Source aliases:**
  - `DOX` -- used in `Zhou_2025_brentuximab.R` (the NONMEM dataset uses the doxorubicin-administration flag as the AVD-coadministration indicator because doxorubicin is given on the same days as the other AVD agents in this regimen).
- **Example models:** `Zhou_2025_brentuximab.R` (power-form effect on ADC clearance: `CL * 2.12^CONMED_AVD` -- ADC clearance is ~2.1-fold higher under A+AVD vs single-agent BV).
- **Notes:** Distinct from `CONMED_CHEMO` (which is nivolumab + platinum-based chemotherapy). The A+AVD regimen is the standard chemotherapy backbone for frontline classical Hodgkin lymphoma; promote to general scope if a second BV paper reports the same A+AVD covariate effect with a comparable encoding.

### CONMED_AZA (**canonical for concomitant azathioprine**)
- **Description:** 1 = on concomitant azathioprine at the PK observation, 0 = not on azathioprine.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant azathioprine).
- **Source aliases:** `AZA` -- used in `Rosario_2015_vedolizumab.R`.
- **Example models:** `Rosario_2015_vedolizumab.R` (power-form on CLL: `CLL * 0.998^CONMED_AZA`; effect ~= null).
- **Notes:** Thiopurine immunomodulator used as maintenance therapy in IBD. Standard convention is baseline-use-only, but time-varying use is permitted; document per-model.

### CONMED_AZOLE (**canonical for concomitant azole antifungal therapy (CYP3A4/P-gp inhibitor)**)
- **Description:** 1 = patient coadministered an azole antifungal (itraconazole, voriconazole, fluconazole, ketoconazole, posaconazole, isavuconazole, or another systemic azole) during the observation interval, 0 = no concomitant azole antifungal. Time-varying per subject because azole exposure starts and stops during the observation period.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant azole antifungal).
- **Source aliases:**
  - `AZOLE` -- used in `Kirubakaran_2022_tacrolimus.R`.
- **Example models:** `Kirubakaran_2022_tacrolimus.R` (state-dependent typical CL/F: 21.1 L/h without azole and 4.2 L/h with azole, an 80% reduction; also a state-dependent BSV magnitude on CL/F: 61% CV without azole vs 89.5% CV with azole).
- **Notes:** Azoles are mechanism-based CYP3A4 and P-glycoprotein inhibitors with different inhibitor potencies (itraconazole > voriconazole > fluconazole). The per-model `covariateData[[CONMED_AZOLE]]$notes` must document (1) which azoles are pooled into the indicator, (2) any post-cessation lag (Kirubakaran 2022 carries `CONMED_AZOLE = 1` for one week after azole discontinuation to allow tacrolimus apparent clearance to stabilize given itraconazole's long half-life), and (3) whether the indicator is a baseline-only proxy or a true time-varying flag.

### CONMED_CHEMO (**canonical for anti-PD-(L)1 mAb + chemotherapy combination indicator**)
- **Description:** 1 = subject is receiving an anti-PD-(L)1 monoclonal antibody in combination with platinum-based chemotherapy (gemcitabine + cisplatin, pemetrexed + cisplatin, paclitaxel + carboplatin, or platinum-doublet); 0 = otherwise. Encodes chemotherapy coadministration as a study-design covariate on the antibody's CL.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no chemotherapy coadministration -- monotherapy or, where applicable, a non-chemotherapy combination such as anti-PD-1 + anti-CTLA-4).
- **Source aliases:**
  - `CHEMO` -- used in `Zhang_2019_nivolumab.R`.
  - `MONOTR` -- used in `Kuchimanchi_2024_dostarlimab.R` (the paper's structural-equation indicator for monotherapy; the canonical column carries the inverse value, i.e. `CONMED_CHEMO = 1 - MONOTR`, so the canonical column is 1 when the patient is on combo-chemotherapy).
- **Example models:** `Zhang_2019_nivolumab.R` (exponential effect on baseline CL: `exp(-0.104)` ~= 0.90 fold, i.e. ~9.7% lower CL relative to monotherapy), `Kuchimanchi_2024_dostarlimab.R` (multiplicative effect on baseline CL: `1 - 0.0779` = 0.922, i.e. 7.79% lower CL on dostarlimab + carboplatin/paclitaxel relative to dostarlimab monotherapy).
- **Notes:** Promoted from specific to general scope on 2026-04-27 after the Kuchimanchi 2024 dostarlimab + carboplatin/paclitaxel analysis ratified the same pooling convention (any chemotherapy backbone collapsed into a single binary indicator). The two papers use different functional forms for the effect on CL -- Zhang 2019 uses `exp(theta * CONMED_CHEMO)` (exponential) and Kuchimanchi 2024 uses `(1 + theta * CONMED_CHEMO)` (multiplicative); these are different parameterisations of the same underlying study-design indicator and the canonical column meaning is unchanged. Document the per-model functional form in `covariateData[[CONMED_CHEMO]]$notes`.

### CONMED_EFV (**canonical for concomitant efavirenz indicator**)
- **Description:** 1 = subject is receiving efavirenz (EFV)-based antiretroviral therapy as the third agent in a combination ART regimen; 0 = subject is on the comparator regimen specified by the source paper (typically a protease-inhibitor-based regimen such as standard lopinavir/ritonavir 4:1). Efavirenz is a CYP3A and UGT inducer, so the indicator is used to flag PXR-mediated induction of metabolic clearance for co-administered antiretrovirals.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-EFV reference regimen, paper-defined; e.g., LPV/r 4:1 in Tikiso 2021).
- **Source aliases:**
  - `EFV` -- used in `Tikiso_2021_abacavir.R` (the dataset's paper-defined indicator, 1 = on EFV-based ART, 0 = on standard LPV/r 4:1).
- **Example models:** `Tikiso_2021_abacavir.R` (multiplicative effect on apparent oral clearance: `cl *= (1 + 0.120 * CONMED_EFV)`; +12.0% relative to the LPV/r 4:1 reference).
- **Notes:** Specific scope because the comparator regimen (LPV/r 4:1 in Tikiso 2021) is paper-defined; future ART population-PK models that test EFV-vs-other contrasts should extend the example list when the comparator matches, or register a finer-grained sibling indicator otherwise.

### CONMED_EIAED (**canonical for concomitant enzyme-inducing antiepileptic drug indicator**)
- **Description:** 1 = subject is taking at least one enzyme-inducing antiepileptic drug (EIAED) such as carbamazepine, phenobarbital, or phenytoin during the study; 0 = no EIAED coadministration. EIAEDs induce hepatic metabolism (CYP3A4/UGT-mediated pathways) and increase the clearance of co-administered antiepileptic drugs and their active metabolites.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no EIAED coadministration).
- **Source aliases:**
  - `MED` (with the value convention inverted: source codes 1 = absence of EIAEDs and 0 = presence; canonical inverts this so 1 = presence and 0 = absence) -- used in `Rodrigues_2017_oxcarbazepine.R`.
- **Example models:** `Rodrigues_2017_oxcarbazepine.R` (exponential effect on MHD apparent clearance: `cl_mhd *= exp(e_eiaed_cl_mhd * (1 - CONMED_EIAED))`; CL_MHD is 29.3% higher with EIAEDs vs without, encoded as +0.257 on the absence-indicator in the source paper).
- **Notes:** Per-model `covariateData[[CONMED_EIAED]]$notes` should list the specific EIAEDs counted as "EIAED = 1" since inclusion criteria vary across antiepileptic-drug studies. Rodrigues 2017 counts carbamazepine, phenobarbital, and phenytoin as EIAEDs; other AEDs in the dataset (vigabatrin, clobazam, valproic acid, clonazepam, lamotrigine, diazepam, ethosuccimide, progabide) are not. The source paper uses an "absence of EIAED" indicator (`MED = 1` if no EIAED, `MED = 0` if EIAED present); the canonical column inverts this so that the 1 group is "on EIAED", matching the convention for other CONMED_* indicators.

### CONMED_EZE (**canonical for concomitant ezetimibe coadministration indicator**)
- **Description:** 1 = patient is taking ezetimibe (with or without other lipid-lowering comedication), 0 = not on ezetimibe.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (not on ezetimibe).
- **Source aliases:**
  - Derived from an ezetimibe-identifier column in the source.
- **Example models:** `Kuchimanchi_2018_evolocumab.R` (multiplicative effect 1.20 on Vmax: `Vmax * 1.20^CONMED_EZE`; labeled "Statin + ezetimibe exponent" in Kuchimanchi 2018 Table 3 because ~99% of ezetimibe users in the dataset were also on a conmed_statin, so the effect effectively captures combination therapy).
- **Notes:** Scope: specific because Kuchimanchi 2018 interprets the ezetimibe indicator as a combination-therapy marker rather than a pure ezetimibe effect. Future popPK/PD models with cleaner ezetimibe separation should add themselves here or register a more specific canonical.

### CONMED_H2RA (**canonical for concomitant H2-receptor-antagonist use**)
- **Description:** 1 = patient on concomitant histamine H2-receptor-antagonist therapy (e.g., ranitidine, famotidine), 0 = no CONMED_H2RA use. Captures another class of gastric-pH-modifying co-medication that may reduce bioavailability of pH-sensitive orally administered drugs.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no CONMED_H2RA use).
- **Source aliases:**
  - `H2` -- used in `Goel_2016_Sonidegib.R` (Goel 2016 dataset; defined as "significant" CONMED_H2RA use, i.e. duration of CONMED_H2RA use >= 80% of the PK assessment phase).
- **Example models:** `Goel_2016_Sonidegib.R` (multiplicative effect on F: `0.996^CONMED_H2RA`, no clinically meaningful effect; reported alongside `CONMED_PPI` for completeness).
- **Notes:** Per-model `covariateData[[CONMED_H2RA]]$notes` must document the operational definition (Goel 2016: >= 80% of PK assessment phase). Distinct from `CONMED_PPI`.

### CONMED_IFNB1A (**canonical for concomitant interferon beta-1a coadministration indicator**)
- **Description:** 1 = patient coadministered subcutaneous interferon beta-1a (Rebif or equivalent recombinant IFN beta-1a product) during the observation interval, 0 = no concomitant IFN beta-1a. Time-varying per subject because the source studies enrol both monotherapy and IFN beta-1a combination periods.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant IFN beta-1a).
- **Source aliases:**
  - `IFNB1A` -- used in `Savic_2017_cladribine.R`.
- **Example models:** `Savic_2017_cladribine.R` (multiplicative effect on cladribine non-renal clearance: `cl_nonrenal *= (1 + e_ifn_clnr * CONMED_IFNB1A)` with `e_ifn_clnr = 0.21`, i.e. a 21% increase in non-renal CL when coadministered with IFN beta-1a).
- **Notes:** Captured at the dose-record level in Savic 2017 (multiple-dose study 26486 alternated between cladribine alone and cladribine + IFN beta-1a periods). The interaction mechanism is not definitively established in the source paper; Savic 2017 discusses that the observed effect could alternatively be modelled on bioavailability or interpreted as a period-effect / interoccasion variability artefact. Future cladribine + immunomodulator extractions should reuse this canonical when IFN beta-1a is the specific concomitant agent; register a sibling canonical (e.g. `CONMED_IFNB1B`, `CONMED_IFNALPHA`) if a different interferon species is intended.

### CONMED_IPI_1Q6W (**canonical for nivolumab + ipilimumab 1 mg/kg q6w combination indicator**)
- **Description:** 1 = subject is receiving nivolumab in combination with ipilimumab 1 mg/kg every 6 weeks (continuous maintenance); 0 = otherwise.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (any non-1Q6W regimen -- monotherapy, chemotherapy combination, or another ipilimumab schedule).
- **Source aliases:**
  - `IPI1Q6W` -- used in `Zhang_2019_nivolumab.R`.
- **Example models:** `Zhang_2019_nivolumab.R` (exponential effect on baseline CL: `exp(0.159)` ~= 1.17 fold increase relative to monotherapy).
- **Notes:** Paired with `CONMED_IPI_3Q3W`. See the CONMED_IPI_3Q3W note for how the other ipilimumab schedules collapse into the reference group.

### CONMED_IPI_3Q3W (**canonical for nivolumab + ipilimumab 3 mg/kg q3w combination indicator**)
- **Description:** 1 = subject is receiving nivolumab in combination with ipilimumab 3 mg/kg every 3 weeks (4-dose induction); 0 = otherwise. Encodes the high-intensity ipilimumab combination regimen as a study-design covariate on nivolumab CL.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (any non-3Q3W regimen -- monotherapy, chemotherapy combination, or another ipilimumab schedule).
- **Source aliases:**
  - `IPI3Q3W` -- used in `Zhang_2019_nivolumab.R`.
- **Example models:** `Zhang_2019_nivolumab.R` (exponential effect on baseline CL: `exp(0.227)` ~= 1.25 fold increase relative to monotherapy).
- **Notes:** Paired with `CONMED_IPI_1Q6W`; both indicators can coexist in one population, but a single subject has at most one set to 1 in the Zhang 2019 cohort. The remaining ipilimumab schedules (1 mg/kg q3w x 4 induction, 1 mg/kg q12w) had no statistically significant effect on nivolumab CL and were therefore folded into the reference (0) group along with monotherapy, leaving only IPI3Q3W and IPI1Q6W as named non-reference indicators.

### CONMED_IPI_ANY (**canonical for any-ipilimumab-coadministration indicator**)
- **Description:** 1 = subject is receiving nivolumab in combination with any ipilimumab regimen (regardless of dose or schedule); 0 = nivolumab monotherapy or nivolumab + chemotherapy. Encodes the "is there ipilimumab in the regimen" question as a single binary covariate, distinct from the regimen-specific CONMED_IPI_3Q3W and CONMED_IPI_1Q6W indicators above.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no ipilimumab coadministration).
- **Source aliases:**
  - `IPICO` -- used in `Zhang_2019_nivolumab.R`.
- **Example models:** `Zhang_2019_nivolumab.R` (additive effect on the time-varying-CL Emax parameter: Emax += -0.0668 when CONMED_IPI_ANY = 1).
- **Notes:** Logically the union of the regimen-specific indicators (CONMED_IPI_3Q3W, CONMED_IPI_1Q6W, plus the unmodeled 1 mg/kg q3wx4 and 1 mg/kg q12w schedules). Zhang 2019 uses it on the *time-varying* Emax (a different structural parameter from baseline CL), which is why it coexists with the regimen-specific indicators on baseline CL rather than substituting for them.

### CONMED_METFORMIN (**canonical for concomitant metformin co-administration indicator**)
- **Description:** 1 = subject is on concomitant metformin during the modelled treatment period, 0 = not. Time-fixed in source datasets where metformin is a study-design "add-on" arm (e.g., Retlich 2015 Study 4 add-on-to-metformin design); permits a time-varying form for cohorts with on/off metformin transitions, document per-model via `covariateData[[CONMED_METFORMIN]]$notes`.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant metformin).
- **Source aliases:** none known.
- **Example models:** `Retlich_2015_linagliptin.R` (1 = study 4 add-on-to-metformin cohort; multiplicative effect on linagliptin relative bioavailability F: +69% F for metformin co-administration vs the monotherapy reference; the effect is attributed to a metformin -- linagliptin drug-drug interaction consistent with a separately published DDI study, Graefe-Mody 2009).
- **Notes:** Follows the `CONMED_*` concomitant-medication pattern (AZA / MP / MTX / AMINO / NSAID / PARA / AD / RITUX / AED / CHEMO / EIAED / EFV / AZOLE). Metformin is a widely-co-prescribed first-line T2DM oral antidiabetic; future T2DM-popPK / -DDI extractions should reuse this canonical. Ratified canonically alongside the Retlich 2015 linagliptin extraction.

### CONMED_MP (**canonical for concomitant 6-mercaptopurine**)
- **Description:** 1 = on concomitant 6-mercaptopurine (6-MP), 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant 6-MP).
- **Source aliases:** `MP` -- used in `Rosario_2015_vedolizumab.R`.
- **Example models:** `Rosario_2015_vedolizumab.R` (power-form on CLL: `CLL * 1.04^CONMED_MP`).
- **Notes:** Second thiopurine immunomodulator used in IBD maintenance; typically mutually exclusive with `CONMED_AZA` for a given subject.

### CONMED_MTX (**canonical for concomitant methotrexate**)
- **Description:** 1 = on concomitant methotrexate, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant methotrexate).
- **Source aliases:** `MTX` -- used in `Rosario_2015_vedolizumab.R`.
- **Example models:** `Rosario_2015_vedolizumab.R` (power-form on CLL: `CLL * 0.983^CONMED_MTX`).
- **Notes:** Immunomodulator used especially in CD maintenance. Generic concomitant-MTX indicator that may also appear in non-IBD models; start as scope: general.

### CONMED_NSAID (**canonical for concomitant NSAID use**)
- **Description:** 1 = on concomitant non-steroidal anti-inflammatory drug (NSAID) therapy at baseline, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant NSAID use; typical patient).
- **Source aliases:**
  - `NSAID` -- used in `Li_2019_abatacept.R`.
- **Example models:** `Li_2019_abatacept.R` (exponential effect on CL: `CL * exp(0.0640 * CONMED_NSAID)`; ~6.6% higher CL, not clinically relevant per Li 2019).
- **Notes:** Baseline-use-only in Li 2019; time-varying use is permitted, document per-model. Follows the `CONMED_*` concomitant-medication pattern established for IBD models (AZA / MP / MTX / AMINO).

### CONMED_PARA (**canonical for concomitant paracetamol (acetaminophen) use**)
- **Description:** 1 = on concomitant paracetamol (acetaminophen) at the observation, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant paracetamol).
- **Source aliases:**
  - `PCM` -- used in `Plan_2012_pain.R` (DDMORE Foundation Model Repository entry DDMODEL00000194).
- **Example models:** `Plan_2012_pain.R` (additive shift on the placebo-effect logit of the typical pain-score lambda: `phl = logit(TVLAM) + 0.364 * CONMED_PARA`; mean pain score ~9% higher on the 0-10 Likert scale during paracetamol use).
- **Notes:** Distinct from `CONMED_NSAID` -- paracetamol is not classed as an NSAID (no anti-inflammatory mechanism, distinct AE profile). Time-varying use is permitted; the daily Likert measurements in Plan 2012 carry PCM as a per-observation flag. Follows the `CONMED_*` concomitant-medication pattern established for IBD models (AZA / MP / MTX / AMINO) and `CONMED_NSAID` (Li 2019). Ratified canonically alongside the Plan 2012 DDMORE extraction.

### CONMED_PPI (**canonical for concomitant proton-pump inhibitor use**)
- **Description:** 1 = patient on concomitant proton-pump inhibitor (CONMED_PPI) therapy, 0 = no CONMED_PPI use. Captures gastric-pH-elevating co-medication that can reduce the bioavailability of solubility-limited (typically weakly basic) orally administered drugs.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no CONMED_PPI use).
- **Source aliases:**
  - `CONMED_PPI` -- used in `Goel_2016_Sonidegib.R` (Goel 2016 dataset; defined as "significant" CONMED_PPI use, i.e. duration of CONMED_PPI use >= 80% of the PK assessment phase).
- **Example models:** `Goel_2016_Sonidegib.R` (multiplicative effect on F: `0.696^CONMED_PPI`, ~30% lower bioavailability under CONMED_PPI coadministration).
- **Notes:** Per-model `covariateData[[CONMED_PPI]]$notes` must document the operational definition (e.g., Goel 2016 requires CONMED_PPI use covering >= 80% of the PK assessment window; other studies may use a simpler ever-vs-never indicator or a per-record time-varying flag). Distinct from `CONMED_H2RA` (H2-receptor antagonist) which acts on gastric pH via a different mechanism.

### CONMED_PROBENECID (**canonical for concomitant probenecid co-administration indicator**)
- **Description:** 1 = subject coadministered probenecid (organic-anion transport inhibitor; multidrug resistance protein (MRP) family modulator with documented activity at the blood-brain barrier, blood-CSF barrier, and renal tubular secretion sites), 0 = no concomitant probenecid. Time-varying per subject because probenecid exposure can start and stop within the observation period.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no probenecid coadministration).
- **Source aliases:**
  - Paper-specific day-2 treatment-arm indicator -- used in `Xie_2000_m3g_rat.R` (the paper distinguishes a probenecid-treatment arm from a control arm by experimental day; the canonical column is 1 during the day-2 probenecid co-infusion and 0 otherwise, including on day 1 of the probenecid arm before the probenecid loading dose).
- **Example models:** `Xie_2000_m3g_rat.R` (multiplicative exponential effect on the unbound BBB influx clearance CL_u,in: `cluin *= exp(e_conmed_probenecid_cluin * CONMED_PROBENECID)` with `e_conmed_probenecid_cluin = log(0.17 / 0.11) = 0.4353`, i.e. a 1.55-fold increase in CL_u,in into rat brain ECF under probenecid co-administration; CL_u,out and the intercompartmental brain clearance Q_br are not statistically affected by probenecid in the paper's model selection and so are not paired with this canonical).
- **Notes:** Scope: specific because the only on-disk source is a preclinical microdialysis rat BBB-transport paper (Xie 2000) and the column meaning is intrinsically tied to the day-2 probenecid co-infusion design rather than a general clinical-coadministration indicator. Per-model `covariateData[[CONMED_PROBENECID]]$notes` must document the dose / regimen of probenecid used (Xie 2000: 70 umol/kg IV loading dose plus 70 umol/kg/h constant infusion in male Sprague-Dawley rats), the parameter the indicator modifies, and any per-subject vs per-record time-varying convention.

### CONMED_RIF_LPVR4 (**canonical for concomitant rifampicin-based antitubercular treatment with super-boosted lopinavir/ritonavir 4:4 indicator**)
- **Description:** 1 = subject is receiving rifampicin-based antitubercular treatment together with super-boosted lopinavir/ritonavir 4:4 (extra ritonavir added to standard 4:1 LPV/r to counter rifampicin-driven LPV induction); 0 = subject is on the comparator regimen specified by the source paper (typically standard LPV/r 4:1 without rifampicin in Tikiso 2021). Used to flag the combined induction effect of rifampicin (PXR-mediated UGT induction) plus extra ritonavir on co-administered antiretroviral PK.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-RIF + standard LPV/r 4:1 reference).
- **Source aliases:**
  - `RIF` -- used in `Tikiso_2021_abacavir.R` (the dataset's paper-defined indicator, 1 = on rifampicin-based TB treatment with super-boosted LPV/r 4:4, 0 = standard LPV/r 4:1 or EFV).
- **Example models:** `Tikiso_2021_abacavir.R` (multiplicative effect on bioavailability: `f_depot *= (1 + (-0.294) * CONMED_RIF_LPVR4)`; -29.4% relative to the LPV/r 4:1 reference).
- **Notes:** Specific scope because the joint rifampicin + super-boosted-LPV/r contrast is paper-defined; the underlying induction is plausibly rifampicin-driven (the paper's discussion notes that LPV concentrations were similar with vs without super-boosting, weakening a separate ritonavir contribution), but the canonical column captures the joint indicator the paper modeled. Future TB/HIV co-treatment models that test the same regimen contrast should extend the example list rather than register a new canonical.

### CONMED_RITUX (**canonical for concomitant rituximab combination therapy**)
- **Description:** 1 = on concomitant rituximab combination therapy (with or without backbone chemotherapy), 0 = not on rituximab.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant rituximab; single-agent or non-rituximab combination).
- **Source aliases:**
  - `RITUX` -- used in `Wu_2024_inotuzumab.R`.
- **Example models:** `Wu_2024_inotuzumab.R` (additive fractional change on CL1: `CL1 * (1 + (-0.132) * CONMED_RITUX)` ~= 13% lower CL1 with concomitant rituximab).
- **Notes:** Wu 2024 Table 3 footnote b explicitly flips the reference category vs. the predecessor Garrett 2019 adult model: in Garrett 2019 the reference was "with rituximab" (RITUX = 0 meant on-rituximab), whereas in Wu 2024 the reference is "without rituximab" (RITUX = 0 means no concomitant rituximab). Future models that pool an analogous rituximab-combination cohort with a single-agent reference should use this canonical with the Wu 2024 sign convention; if a paper retains the Garrett 2019 reverse-coded convention, document the value transformation in `covariateData[[CONMED_RITUX]]$notes` (`CONMED_RITUX = 1 - source$RITUX`) rather than registering a second canonical. Ratified canonically on 2026-04-26.

### CONMED_SPART (**canonical for spartalizumab (PDR001, anti-PD-1) coadministration indicator**)
- **Description:** 1 = the analyzed therapeutic mAb is coadministered with spartalizumab (PDR001, anti-PD-1 IgG4), 0 = no spartalizumab coadministration. Time-fixed per subject in source analyses to date.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no spartalizumab coadministration -- monotherapy or combination with non-spartalizumab agents such as a hypomethylating agent).
- **Source aliases:**
  - `HASPDR` -- used in `Xu_2023_MBG453.R` (Monolix supplement Appendix S2; the source describes the column as "this patient HAS received PDR001 [spartalizumab, anti PD-1 mAb]").
- **Example models:** `Xu_2023_MBG453.R` (exponential effect on CL: `exp(0.0194 * CONMED_SPART)`; not statistically significant in the full covariate model but retained because Xu 2023 used the full-covariate-model approach).
- **Notes:** Parallels `COMBO_NIVO` (ipilimumab + nivolumab) and `COMBO_DURVA` (durvalumab combinations) but for spartalizumab. Promote to general scope if a second paper reports a spartalizumab-coadministration covariate with comparable encoding.

### CONMED_STATIN (**canonical for concomitant conmed_statin (HMG-CoA reductase inhibitor) therapy**)
- **Description:** 1 = patient coadministered a conmed_statin (HMG-CoA reductase inhibitor) during the study, 0 = no conmed_statin coadministration.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no conmed_statin coadministration).
- **Source aliases:** none known.
- **Example models:** `Martinez_2019_alirocumab.R` (additive effect on linear clearance CLL: `CLL = TVCLL + COV1*(WT-82.9) + COV2*CONMED_STATIN`; +0.00644 L/h when conmed_statin is coadministered).
- **Notes:** Per-model `covariateData[[CONMED_STATIN]]$notes` must document which statins and dose thresholds are included in the "CONMED_STATIN = 1" category, since inclusion criteria vary by study. Martinez 2019 codes CONMED_STATIN = 1 for coadministration of rosuvastatin (< 20 mg/day), atorvastatin (< 40 mg/day), or simvastatin (any dose); other conmed_statin regimens are coded as 0.

### CONMED_STATIN_MONO (**canonical for concomitant statin-monotherapy indicator**)
- **Description:** 1 = patient is on a conmed_statin and no other lipid-lowering comedication (conmed_statin monotherapy), 0 = not on conmed_statin monotherapy (either on no lipid-lowering therapy or on a multi-drug lipid-lowering combination such as conmed_statin + ezetimibe).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (not on conmed_statin monotherapy).
- **Source aliases:**
  - Derived from a conmed_statin-identifier column in the source (any of atorvastatin, rosuvastatin, simvastatin, lovastatin, pravastatin, pitavastatin, fluvastatin) with an AND over "no other lipid-lowering comedication".
- **Example models:** `Kuchimanchi_2018_evolocumab.R` (multiplicative effect 1.13 on Vmax: `Vmax * 1.13^CONMED_STATIN_MONO`).
- **Notes:** Scope: specific because Kuchimanchi 2018 narrowly defines the conmed_statin covariate as monotherapy only ("patients on a conmed_statin only and no other comedication"). Mutually compatible with `CONMED_EZE`: a subject on conmed_statin+ezetimibe has `CONMED_STATIN_MONO = 0` and `CONMED_EZE = 1`; a subject on conmed_statin alone has `CONMED_STATIN_MONO = 1` and `CONMED_EZE = 0`; a subject on no lipid-lowering therapy has both 0. Future popPK/PD models that adopt a broader "any conmed_statin" definition should register a separate `CONMED_STATIN` or `CONMED_STATIN` canonical rather than reusing this name.

### CONMED_STEROID (**canonical for baseline/concomitant systemic corticosteroid use**)
- **Description:** 1 = patient on systemic corticosteroid therapy at baseline (typically continued as concomitant medication during the study), 0 = no baseline corticosteroid use. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no baseline corticosteroid use).
- **Source aliases:**
  - `BSTEROID` -- used in `Narwal_2013_sifalimumab.R` and `Zheng_2016_sifalimumab.R`.
- **Example models:** `Narwal_2013_sifalimumab.R` (multiplicative on CL: `CL * (1 + 0.195 * CONMED_STEROID)`), `Zheng_2016_sifalimumab.R` (multiplicative on CL `(1 + 0.11 * CONMED_STEROID)` and on V1 `(1 - 0.09 * CONMED_STEROID)` in the SLE phase IIb cohort, which was ~85% conmed_steroid-treated at baseline).
- **Notes:** Distinct from `PRICORT`, which is strictly a prior (pre-study) indicator. `CONMED_STEROID` captures concurrent corticosteroid use at / from study baseline in diseases where background conmed_steroid use is standard of care (SLE, severe asthma, etc.). When a future paper needs the two jointly, both can coexist on the same subject. The name `STEROID_BL` was used as an alias in earlier register drafts and is retired; use `CONMED_STEROID` for all future models.

### PRICORT (**canonical for prior corticosteroid use indicator**)
- **Description:** 1 = patient received systemic corticosteroid treatment prior to study entry, 0 = no prior corticosteroid use. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no prior corticosteroid use).
- **Source aliases:** none.
- **Example models:** `Ma_2020_sarilumab_das28crp.R` (multiplicative on DAS28-CRP Kout: `Kout * theta^PRICORT`), `Ma_2020_sarilumab_anc.R` (power-form on Emax: `Emax * 0.819^PRICORT`).
- **Notes:** Ma 2020 applies it as a multiplicative effect of the form `param * theta^PRICORT` in both DAS28-CRP and ANC PD models. Generally applicable clinical-history indicator.

### PRIOR_ANTHRACYCLINE_DOSE (**canonical for prior cumulative anthracycline dose**)
- **Description:** Cumulative dose of anthracycline chemotherapy received by the subject prior to the first dose analysed in the current popPK / popPK-PD model, expressed in doxorubicin-equivalent body-surface-area-normalised mg/m^2. Time-fixed per subject for the analysis window (the running cumulative anthracycline dose at the first observed dose).
- **Units:** mg/m^2 (doxorubicin-equivalent)
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- typically used as a linear shift (`(1 + theta * (PRIOR_ANTHRACYCLINE_DOSE - ref))`) on a baseline parameter (e.g., baseline cardiac troponin I before the next anthracycline cycle). Reference values observed: 90 mg/m^2 (Kunarajah 2017, cohort median).
- **Source aliases:**
  - `PCAMT` -- Kunarajah 2017 NM-TRAN convention ("Prior Cumulative Anthracyclines aMounT"; doxorubicin-equivalent mg/m^2).
- **Example models:** `Kunarajah_2017_doxorubicin.R` (linear shift on baseline cardiac troponin I: `bl_cTnI * (1 + 0.00308 * (PRIOR_ANTHRACYCLINE_DOSE - 90))` -- ~0.31% increase in baseline cTnI per 1 mg/m^2 of prior cumulative anthracycline exposure).
- **Notes:** Distinct from `PRIOR_ANTICANCER` (a binary modality indicator, 1 = any prior anticancer therapy) -- `PRIOR_ANTHRACYCLINE_DOSE` carries the actual cumulative dose, restricted to the anthracycline drug class (doxorubicin, daunorubicin, epirubicin, idarubicin), and is the column needed when the source paper's effect is dose-response in the prior-exposure regime rather than presence / absence. When a paper records anthracycline exposure as anthracycline-class-by-class doses and the model effect aggregates them, sum to a single doxorubicin-equivalent value before populating this column (use the published bone-marrow / cardiotoxicity isoeffective conversion factors). When a paper distinguishes the type of anthracycline (e.g., doxorubicin vs daunorubicin separately), register parallel canonicals (`PRIOR_DOXORUBICIN_DOSE`, `PRIOR_DAUNORUBICIN_DOSE`) rather than overloading this name. Scope: specific because the column meaning is intrinsically tied to anthracycline-class chemotherapy exposure; promote to general if a second paper ratifies the same definition.

### PRIOR_ANTICANCER (**canonical for prior anticancer therapy of any modality**)
- **Description:** 1 = subject received any prior anticancer therapy (cytotoxic chemotherapy, radiotherapy, hormonal therapy, targeted therapy, immunotherapy, or surgical debulking with adjuvant intent) before the start of the analyzed treatment, 0 = treatment-naive. Broader than `LINE_1L`, which is specifically a treatment-line indicator restricted to systemic drug therapy lines. `PRIOR_ANTICANCER` captures the full clinical concept of prior cancer treatment exposure as used in cytotoxic-chemotherapy myelosuppression analyses, where any prior anticancer modality may have depleted the bone-marrow proliferating pool and therefore affects baseline ANC.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (treatment-naive; no prior anticancer therapy of any modality).
- **Source aliases:**
  - `PC` (Kloft 2006 / Netterberg 2017 NM-TRAN convention for "previous anticancer therapy"; values 0 = naive, 1 = had prior anticancer therapy) -- used in `Netterberg_2017_docetaxel.R`.
- **Example models:** `Netterberg_2017_docetaxel.R` (multiplicative effect on baseline ANC: `BACOV *= (1 + theta * PRIOR_ANTICANCER)` with theta = -0.147; prior-anticancer patients have ~14.7% lower baseline ANC than treatment-naive patients).
- **Notes:** Distinct from `LINE_1L` (which is the inverse semantics for systemic-drug therapy lines only) and `PRIOR_TNF` / `PRIOR_BIO` (which are modality-specific to anti-TNF / biologic exposure in inflammatory-disease cohorts). Use `PRIOR_ANTICANCER` when the source paper's covariate counts any anticancer modality (including radiotherapy and surgery) as prior exposure. When a future paper restricts the indicator to cytotoxic chemotherapy alone, use `LINE_1L` (with values inverted: paper's `PRIOR_CHEMO = 1 - LINE_1L`). When a paper distinguishes prior chemotherapy from prior radiotherapy, register a parallel `PRIOR_RADIATION` canonical.

### PRIOR_BIO (**canonical for prior biologic exposure**)
- **Description:** 1 = subject previously treated with any biologic (broader than `PRIOR_TNF`: includes anti-TNF agents plus anti-integrin, anti-IL-12/23, anti-IL-17, anti-IL-23, anti-IL-6, etc.), 0 = biologic-naive.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (biologic-naive).
- **Source aliases:**
  - `bio-naive` (Aguiar 2021 source-paper variable, with the indicator inverted: paper's `bio-naive = 1 - PRIOR_BIO`). Effect coefficients in the source paper apply to `bio-naive`; in `model()` derive `bio_naive <- 1 - PRIOR_BIO` to preserve the paper's reported coefficient.
- **Example models:** `Aguiar_2021_ustekinumab.R` (Aguiar 2021 Table 2 footnote a; multiplicative effect on CL: factor `(1 - 0.227 * (1 - PRIOR_BIO))`, so bio-naive patients have ~23% lower CL than previously-exposed patients).
- **Notes:** Distinct from `PRIOR_TNF` (a strict subset). Use `PRIOR_BIO` when the source paper's covariate counts any biologic as prior exposure (anti-TNF, anti-integrin, anti-IL-12/23, anti-IL-17, anti-IL-23, anti-IL-6, etc.); use `PRIOR_TNF` when the source paper specifically tested anti-TNF exposure. When the source paper uses the inverted "bio-naive" indicator (1 = naive), document the inversion in `covariateData[[PRIOR_BIO]]$notes` and apply `1 - PRIOR_BIO` in `model()` so the canonical column stores 1 = previously exposed.

### PRIOR_TNF (**canonical for prior anti-TNF biologic exposure indicator**)
- **Description:** 1 = subject previously treated with an anti-TNF (tumor necrosis factor) inhibitor, 0 = TNF-naive.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (TNF-naive).
- **Source aliases:**
  - `PRIORTNF` (all caps, no underscore) -- acceptable alternative spelling.
- **Example models:** `Moein_2022_etrolizumab.R` (multiplicative fractional effect on CL, +4.9%).
- **Notes:** Use when the source paper reports a binary "prior anti-TNF inhibitor" covariate on any PK parameter. Generally applicable across RA/PsA/IBD/axSpA biologic PK models.

## Rheumatoid-arthritis disease-activity covariates

### RHEUMATOID_FACTOR (**canonical for serum rheumatoid factor concentration**)
- **Description:** Serum rheumatoid factor (an autoantibody, predominantly IgM, directed against the Fc portion of IgG) concentration. Baseline value typical; document time-varying use in per-model `notes`.
- **Units:** U/mL or IU/mL (interchangeable in the clinical-PK literature). Document per-model via `covariateData[[RHEUMATOID_FACTOR]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling on a log-transformed value: `(log(RHEUMATOID_FACTOR) / log(ref))^exponent` (or, equivalently, the source-paper form `(LRF / log(ref))^exponent` where `LRF = log(RHEUMATOID_FACTOR)`). Reference value observed: 110 U/mL (Frey 2010, corresponding to LRF = 4.7 in the paper's final-model equation).
- **Source aliases:**
  - `LRF` -- log-transformed RF (natural log of the value in U/mL); Frey 2010 fits the covariate on the log scale and reports the reference as `LRF = 4.7` (i.e., `log(110) ~= 4.7`). The canonical column carries the raw RF concentration in U/mL; the log transform is applied inside `model()`.
  - `RF` -- universal NONMEM/clinical-PK abbreviation; rejected as the canonical name on 2026-04-28 because the bare two-letter abbreviation is uncommon in published popPK papers and could be confused with other shortenings.
- **Example models:** `Frey_2010_tocilizumab.R` (U/mL, reference 110 U/mL == LRF = 4.7; small positive exponent +0.1 on linear CL applied to `log(RHEUMATOID_FACTOR)`).
- **Notes:** RF concentrations span several orders of magnitude across the rheumatoid-arthritis population (Frey 2010 observed range 15-11,800 U/mL across the four phase-III studies; reference paper: Frey 2010 Table I), motivating the log transform before power scaling. The mechanistic rationale (Frey 2010 Discussion, p764) is that RF -- being an anti-IgG autoantibody -- could in principle bind the Fc region of the therapeutic IgG monoclonal antibody and accelerate clearance, but the observed CL effect was small (-4.9% to +6.5% across the observed RF range) and the paper acknowledges that high RF concentrations may also reduce the assay's ability to detect the drug, leading to an apparent CL increase. Ratified canonically on 2026-04-28 alongside the Frey 2010 extraction.

### BLPHYVAS (**canonical for baseline physician's global assessment VAS**)
- **Description:** Baseline Physician's Global Assessment of Disease Activity, 100-mm visual analogue scale (0 = no disease activity, 100 = maximum). Time-fixed per subject.
- **Units:** mm (0-100 VAS)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used as a power term `(BLPHYVAS / <ref>)^exponent`. Reference 66 used in Ma 2020.
- **Source aliases:** none.
- **Example models:** `Ma_2020_sarilumab_das28crp.R`.
- **Notes:** One of the components of the DAS28 composite score; in Ma 2020 it appears as a baseline covariate on the DAS28-CRP disease-activity BASE rather than on the score itself. Applicable to any rheumatology model where baseline physician-assessed disease activity is used as a PK/PD covariate.

### BLHAQ (**canonical for baseline HAQ-DI score**)
- **Description:** Baseline Health Assessment Questionnaire Disability Index (HAQ-DI; 0 = no disability, 3 = maximum disability). Time-fixed per subject.
- **Units:** unitless (0-3 composite score)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used as a power term `(BLHAQ / <ref>)^exponent`. Reference 1.75 used in Ma 2020.
- **Source aliases:** none.
- **Example models:** `Ma_2020_sarilumab_das28crp.R`.
- **Notes:** Patient-reported disability score frequently used as a baseline covariate in rheumatoid-arthritis PK/PD analyses.

### PAIN (**canonical for patient-reported global pain visual analogue score**)
- **Description:** Patient-reported global pain on a 100-mm visual analogue scale (PAIN; 0 = no pain, 100 = worst imaginable pain). Distinct from `BLPHYVAS` (the *physician*'s global assessment of disease activity). Both baseline and time-varying usages are covered; document per-model in `covariateData[[PAIN]]$notes` whether the column is baseline-only.
- **Units:** mm (0-100 VAS).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(PAIN / ref)^exponent`. Reference value observed: 60 in Frey 2013 (approximate dataset median across OPTION/TOWARD).
- **Source aliases:** none known (the PAIN column name is used directly in the Frey 2013 NONMEM control stream).
- **Example models:** `Frey_2013_tocilizumab.R` (baseline; power effect on the indirect-response BASE parameter with reference 60 and exponent 0.062).
- **Notes:** PAIN can take the value 0 in real cohorts (a patient who reports no pain at the assessment), which makes the bare power form `(PAIN/60)^exp` return 0. Frey 2013 documents an explicit 0.010 floor on PAIN in its Table 2 covariate range; the model file applies the same floor inside `model()` so simulation under PAIN = 0 returns a finite BASE rather than collapsing to zero. Canonical name follows the `EOS` / `EASI` convention (no `BL` prefix); baseline-vs-time-varying status is per-model.

### SWOL_28JOINT (**canonical for 28-joint swollen joint count**)
- **Description:** Swollen joint count on the 28-joint (DAS28) scale (integer 0-28; component of the DAS28 composite). Baseline value is typical; document time-varying use in per-model `notes`.
- **Units:** count (0-28)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used as a shifted power term `((SWOL_28JOINT + 1)/(<ref> + 1))^exponent` to avoid the zero-count edge case. Reference value observed: 16 in Li 2019 (approximate dataset median of the popPK cohort).
- **Source aliases:**
  - `SWOL` -- used in `Li_2019_abatacept.R` (Li 2019 Methods abbreviation).
- **Example models:** `Li_2019_abatacept.R` (power effect on CL with exponent 0.0965; not clinically relevant per Li 2019).
- **Notes:** The `_28JOINT` suffix distinguishes this from the 66/68-joint swollen count used in some earlier RA scales -- register a separate canonical (`SWOL_66JOINT` or similar) if a future paper uses a different joint-count scale. Canonical name drops the `BL` prefix to match the `EASI` / `AGE` / `WT` / `ALB` convention where baseline-vs-time-varying status is documented in `covariateData` notes rather than the column name.

## Pharmacogenetics

### CYP2D6 (**canonical for CYP2D6 individual metabolic-activity score**)
- **Description:** Continuous individual-level CYP2D6 metabolic-activity score. The intent is a single canonical column for any CYP2D6 phenotype proxy that the source paper reports as a continuous number (probe-substrate model-based individual clearance, copy-number-corrected expression score, activity-score sum from `*allele` genotypes, etc.); the per-model `covariateData[[CYP2D6]]$units`, `description`, and `notes` document which proxy is in force and the population-median reference value used inside the model. Time-invariant in all known examples (germline genotype or one-time probe-substrate measurement).
- **Units:** Paper-specific -- document per-model (e.g., `ng/L` in Ter Heine 2014 where the value is the dextromethorphan-probe model-based individual CYP2D6 clearance).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a (continuous). Models center on a population median (e.g., 1560 ng/L in Ter Heine 2014); document the reference value per-model in `covariateData[[CYP2D6]]$notes`.
- **Source aliases:**
  - `CYP2D6` -- used directly in `TerHeine_2014_tamoxifen.R`.
- **Example models:** `TerHeine_2014_tamoxifen.R` (power-law effect on the tamoxifen -> endoxifen formation clearance: `(CYP2D6 / 1560)^e_CYP2D6_cl_endx`).
- **Notes:** TODO -- consolidate the various CYP2D6 phenotype encodings (continuous probe-derived activity, copy-number-corrected score, categorical PM/IM/EM/UM phenotype, activity-score sum) into a single canonical `CYP2D6` column with an enumerated `notes`-documented `encoding` field, OR introduce companion canonicals (`CYP2D6_PHENO_GROUP` for the categorical PM/IM/EM/UM grouping, `CYP2D6_ACTSCORE` for the AS sum) when a future model needs the categorical form. The current general-scope continuous canonical is sufficient for the Ter Heine 2014 use case but will likely need refinement as more CYP2D6-aware models are added.

### CYP3A4 (**canonical for CYP3A4 / CYP3A4-and-CYP3A5 individual metabolic-activity score**)
- **Description:** Continuous individual-level CYP3A4 (or combined CYP3A4 + CYP3A5) metabolic-activity score. Same intent and documentation policy as `CYP2D6` above. Some sources measure CYP3A4 alone via probe substrate; others (Ter Heine 2014) report a combined CYP3A4/5 activity because the chosen probe (dextromethorphan N-demethylation) cannot distinguish CYP3A4 from CYP3A5. The per-model `notes` field documents whether the value is CYP3A4-only or the CYP3A4 + CYP3A5 combined score.
- **Units:** Paper-specific -- document per-model (e.g., `ng/L` in Ter Heine 2014 where the value is the dextromethorphan-probe model-based individual CYP3A4/5 clearance).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a (continuous). Models center on a population median (e.g., 44.7 ng/L in Ter Heine 2014); document the reference value per-model.
- **Source aliases:**
  - `CYP3A4` -- used directly in `TerHeine_2014_tamoxifen.R` (the column carries combined CYP3A4 + CYP3A5 activity per the source paper).
  - `CYP3A4/5` -- long form sometimes used in source manuscripts; standardize the column name to `CYP3A4` and document the combined-isoform semantics in per-model `notes`.
- **Example models:** `TerHeine_2014_tamoxifen.R` (power-law effect on the tamoxifen -> endoxifen formation clearance: `(CYP3A4 / 44.7)^e_CYP3A4_cl_endx`).
- **Notes:** TODO -- register the rest of the canonical drug-metabolizing-CYP set prospectively (`CYP1A2`, `CYP2A6`, `CYP2B6`, `CYP2C8`, `CYP2C9`, `CYP2C19`, `CYP2E1`, `CYP3A5`) using the same continuous-individual-activity-score pattern, so future popPK models that report CYP-probe-derived covariates can drop straight into the existing convention rather than re-deliberating the encoding each time. Coordinate with the consolidation TODO on `CYP2D6` so categorical-vs-continuous encoding is handled uniformly across all CYPs.

### CYP3A5_EXPR (**canonical for CYP3A5 expresser status**)
- **Description:** 1 = subject carries at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3, equivalent to one or two A alleles at rs776746); 0 = homozygous CYP3A5*3/*3 (G/G at rs776746) -- i.e., a nonexpresser. Time-fixed per subject (germline genotype).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (CYP3A5*3/*3 nonexpresser).
- **Source aliases:**
  - `X` -- used in `Bergmann_2014_tacrolimus.R` (Bergmann 2014 Table 2 footnote: `X = 1` for *1/*1 and *1/*3; `X = 0` for *3/*3).
  - `CYP3A5 expresser` -- used in `Storset_2014_tacrolimus.R` (Storset 2014 Table 2 final theory-based model; `*1/*1` and `*1/*3` pooled as expressers because Storset 2014 had only n = 3 *1/*1 subjects).
- **Example models:** `Bergmann_2014_tacrolimus.R` (multiplicative effect on tacrolimus CL/F: `theta_CYP3A5 ^ CYP3A5_EXPR`, with `theta_CYP3A5 = 1.60`; expressers have 60% higher apparent oral clearance than nonexpressers), `Storset_2014_tacrolimus.R` (multiplicative effects on apparent plasma clearance: `cl *= 1.30^CYP3A5_EXPR`; and on oral bioavailability: `fdepot *= 0.82^CYP3A5_EXPR`; Storset 2014 Table 2 final theory-based model).
- **Notes:** Distinct from the SNP-pattern canonical `SNP_<GENE>_<RSID>` (which encodes "mutant allele presence" -- 1 = at least one variant allele). For CYP3A5 the *3 allele (rs776746 G) is the variant that abolishes function, so a literal "mutant-allele-presence" indicator (1 = any G allele) would group *1/*3 heterozygotes with the *3/*3 nonexpressers, which is the **opposite** of the clinically meaningful expresser-vs-nonexpresser dichotomy used by every CYP3A5-aware popPK model. The `CYP3A5_EXPR` canonical preserves the expresser-equals-1 orientation directly. Future CYP3A5 papers using a *3/*3 indicator (rather than *1 carrier) should still record their values under `CYP3A5_EXPR` and document the value inversion in `notes` (`CYP3A5_EXPR = 1 - source_indicator`); registering a parallel `CYP3A5_NONEXPR` is discouraged. The canonical name follows the `<gene>_<phenotype>` rather than the `<gene>_<rsid>` pattern because the column captures derived metabolic phenotype rather than raw genotype. Distinct from `CYP3A4` (continuous individual-activity score for CYP3A4 / CYP3A4 + CYP3A5 combined): the binary `CYP3A5_EXPR` is the right fit for source papers that report only the rs776746 genotype, while the continuous `CYP3A4` is for sources that report a probe-substrate-derived activity number. Ratified canonically on 2026-05-08 alongside the Bergmann 2014 extraction.

### CYP3A4_INH (**canonical for concomitant CYP3A4 inhibitor coadministration indicator**)
- **Description:** 1 = subject coadministered any CYP3A4 inhibitor during the study, 0 = no concomitant CYP3A4 inhibitor. Distinct from the `CYP3A4` continuous-activity-score canonical above: `CYP3A4_INH` captures concomitant-medication exposure (a drug-drug-interaction indicator), not intrinsic enzyme activity. Use this canonical when the source paper enters CYP3A4-inhibitor coadministration into the popPK model as a binary indicator, regardless of which inhibitor strengths (strong / moderate / weak) the paper pools into the `1` category.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no CYP3A4 inhibitor coadministration).
- **Source aliases:** none standardized; source datasets typically encode the column as `CYP3AI`, `CYP3A4I`, `CYP3AINH`, or a free-text concomitant-medication indicator. Document the source-column name per-model in `covariateData[[CYP3A4_INH]]$source_name`.
- **Example models:** `Yassen_2025_asundexian.R` (proportional-shift effect on CL/F: `(1 + e_cyp3a4_inh_cl * CYP3A4_INH)` with `e_cyp3a4_inh_cl = -0.0531`; the asundexian dataset pools weak + moderate CYP3A4 inhibitors into the `CYP3A4_INH = 1` category because strong inhibitors were a Phase II exclusion criterion).
- **Notes:** Per-model `covariateData[[CYP3A4_INH]]$notes` must document which inhibitor strengths (strong / moderate / weak) and which specific drug examples are pooled into the `CYP3A4_INH = 1` category, since inclusion criteria vary by study. Future models that need stratified encoding (separate strong / moderate / weak indicators) should register companion canonicals (e.g. `CYP3A4_INH_STRONG`, `CYP3A4_INH_MOD`, `CYP3A4_INH_WEAK`) rather than overloading `CYP3A4_INH`. The complementary CYP3A4-inducer indicator should follow the same pattern as a separate canonical (`CYP3A4_IND`) when first needed. Ratified canonically on 2026-05-08 alongside the Yassen 2025 asundexian extraction.

### CYP3A4_IND (**canonical for concomitant CYP3A4 inducer coadministration indicator**)
- **Description:** 1 = subject coadministered any CYP3A4 inducer during the study, 0 = no concomitant CYP3A4 inducer. Sibling indicator to `CYP3A4_INH`; both capture concomitant-medication exposure (a drug-drug-interaction indicator), not intrinsic enzyme activity (distinct from the `CYP3A4` continuous-activity-score canonical above). Use this canonical when the source paper enters CYP3A4-inducer coadministration into the popPK model as a binary indicator, regardless of which inducer strengths (strong / moderate / weak) the paper pools into the `1` category.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no CYP3A4 inducer coadministration).
- **Source aliases:** none standardized; source datasets typically encode the column as `CYP3AIND`, `CYP3A4IND`, `INDU`, `INDUCER`, or a free-text concomitant-medication indicator. Document the source-column name per-model in `covariateData[[CYP3A4_IND]]$source_name`.
- **Example models:** `Gupta_2016_lenvatinib.R` (multiplicative power-form effect on CL/F: `1.30^CYP3A4_IND` with `e_cyp3a4_ind_cl = log(1.30) ~ 0.262`; the Gupta dataset pools any concomitant CYP3A4 inducer reported in the per-subject medication log into the `CYP3A4_IND = 1` category, with `n = 19` (2.4%) of the 779-subject pooled cohort flagged positive).
- **Notes:** Per-model `covariateData[[CYP3A4_IND]]$notes` must document which inducer strengths (strong / moderate / weak) and which specific drug examples are pooled into the `CYP3A4_IND = 1` category, since inclusion criteria vary by study. Future models that need stratified encoding (separate strong / moderate / weak indicators) should register companion canonicals (e.g. `CYP3A4_IND_STRONG`, `CYP3A4_IND_MOD`, `CYP3A4_IND_WEAK`) rather than overloading `CYP3A4_IND`. Sibling canonical to `CYP3A4_INH` (anticipated by that entry's notes). Ratified canonically alongside the Gupta 2016 lenvatinib extraction.

### APOE4_COUNT (**canonical for APOE-epsilon4 allele count**)
- **Description:** Continuous individual-level APOE-epsilon4 allele count: 0 = non-carrier, 1 = heterozygous (one epsilon4 allele), 2 = homozygous (two epsilon4 alleles). Time-invariant (germline genotype). Models in the Alzheimer's-disease-progression literature treat the 0 / 1 / 2 count as a continuous effect on baseline cognitive score and / or disease-progression slope, with the population-mean carrier-allele count used as the centring value (e.g., 0.72 in the Conrado 2014 CAMD cohort).
- **Units:** (count, 0 / 1 / 2 alleles per subject; population-mean centring value documented per-model)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a (continuous). Models centre on the dataset / population mean APOE-epsilon4 count; document the centring value per-model in `covariateData[[APOE4_COUNT]]$notes`.
- **Source aliases:**
  - `APOE4C` -- used directly in `Conrado_2014_alzheimer.R`. The "C" suffix in the source distinguishes the cleaned continuous APOE-epsilon4 count column (0 / 1 / 2 with `unknown` recoded to the population mean) from the upstream raw `APOE4` column (0 = non-carrier, 1 = heterozygous, 2 = homozygous, 3 = unknown).
- **Example models:** `Conrado_2014_alzheimer.R` (centring 0.72; multiplicative effect on baseline ADAS-Cog and on disease-progression slope: `factor = 1 + e * (APOE4_COUNT - 0.72)` with `e_blapoe4 = 0.0372` on baseline and `e_slapoe4 = 0.195` on slope).
- **Notes:** APOE-epsilon4 carrier status is the strongest established genetic risk factor for late-onset Alzheimer's disease; the allele-count form (rather than a binary carrier indicator) is preferred when the source paper distinguishes heterozygous from homozygous carriers. Future models that report only a binary carrier indicator (any-epsilon4 vs none) should register a separate canonical (`APOE4_CARRIER`) rather than overloading `APOE4_COUNT`. The `unknown` category (often recorded as `APOE4 = 3` in CDISC datasets) is conventionally recoded by the source paper to the population-mean count to avoid dropping subjects; document the recoding rule used per-model. Ratified canonically on 2026-05-06 alongside the Conrado 2014 DDMORE extraction.

### FCGR3A_VV (**canonical for FCGR3A 158 V/V homozygote indicator**)
- **Description:** 1 = subject is homozygous for valine at amino-acid position 158 of the FcgammaRIIIa receptor (V/V), encoded by the rs396991 polymorphism in the FCGR3A gene; 0 = otherwise (heterozygote V/F or homozygote F/F pooled). The dominant V/V vs (V/F + F/F) grouping is the encoding used in the Aguiar 2021 source paper after testing dominant and recessive groupings during covariate model building.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 = V/F heterozygote or F/F homozygote (combined).
- **Source aliases:**
  - `FCGR3A` (genotype string, e.g., `"V/V"` / `"V/F"` / `"F/F"`): derive `FCGR3A_VV = as.integer(FCGR3A == "V/V")`.
  - `rs396991` (raw allele coding, often `"AA"` / `"AC"` / `"CC"` or `"GG"` / `"GT"` / `"TT"` depending on assay strand): map V allele -> 1, F allele -> 0 with the assay-specific allele convention; derive `FCGR3A_VV = as.integer(genotype is V-homozygous)`.
- **Example models:** `Aguiar_2021_ustekinumab.R` (Aguiar 2021 Table 2 footnote; logit-scale effect on subcutaneous bioavailability F: 88.8% in V/V vs 71.0% in V/F + F/F).
- **Notes:** rs396991 (FCGR3A 158V>F) is a well-studied pharmacogenetic polymorphism affecting FcgammaRIIIa-IgG affinity and has been associated with response to several IgG monoclonal antibodies (rituximab, infliximab, ustekinumab). The V allele is the higher-affinity variant. Document the assay-strand allele convention used in the source paper in `covariateData[[FCGR3A_VV]]$notes`. Future models that use a recessive (F/F vs V/* combined) or codominant (additive 0/1/2) coding should register a separate canonical (e.g., `FCGR3A_FF`, `FCGR3A_VCT` for V-allele count) rather than overloading `FCGR3A_VV`.

## Immunogenicity

### ADA_POS (**canonical for anti-drug antibody positive status indicator**)
- **Description:** 1 = antidrug-antibody-positive, 0 = ADA-negative.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (ADA-negative).
- **Source aliases:**
  - `ADA` (semantically "ever positive") -- used in `Zhu_2017_lebrikizumab.R`. When translating from a paper that uses `ADA` as "ever positive," verify the time-frame matches ADA_POS semantics before renaming.
  - `ADA` (time-varying positivity, primary covariate in Xu 2019) -- used in `Xu_2019_sarilumab.R`.
  - `NAB` (neutralizing antibody positive -- used in `Petrov_2024_romiplostim.R`). Strictly a subset of total ADA-positive (ADA antibodies that neutralize the drug's biological effect). Document per-model when the source assay measured NAB only and the canonical column thus excludes binding-only ADA.
  - `ATAPOSNEW` (ADA-positive in the newer/updated-assay cohort) -- used in `Suri_2018_brentuximab.R` as the modern-assay arm of a two-era ADA decomposition.
  - `ADA_POSNEW` (**retired** intermediate name; renamed to `ADA_POS` on 2026-04-29 for consistency across single- and multi-assay models).
- **Example models:** `Clegg_2024_nirsevimab.R`, `Hu_2026_clesrovimab.R`, `Petrov_2024_romiplostim.R`, `Suri_2018_brentuximab.R` (multi-assay; paired with `ADA_POSOLD` and `ADA_MISSING`; `cl *= (1 + 0.125 * ADA_POS)`), `Xu_2019_sarilumab.R`.
- **Notes:** In single-assay studies this is a straightforward binary. In studies pooling data across assay generations (different sensitivity / drug-tolerance characteristics), `ADA_POS` represents modern/current-assay positivity; companion indicators `ADA_POSOLD` and `ADA_MISSING` capture historical-assay-positive and missing-result sub-groups respectively. All three are mutually exclusive; reference is ADA-negative (all three = 0).

### ADA_POSOLD (**canonical for ADA-positive in older-assay study indicator**)
- **Description:** 1 = subject is anti-drug-antibody-positive in a study that used the older lower-sensitivity, lower-drug-tolerance ADA assay; 0 = otherwise.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (ADA-negative or not in an older-assay study). Mutually exclusive with `ADA_POS` and `ADA_MISSING`.
- **Source aliases:**
  - `ATAPOSOLD` -- used in `Suri_2018_brentuximab.R`.
- **Example models:** `Suri_2018_brentuximab.R` (multiplicative additive effect on ADC clearance: `cl *= (1 + 0.177 * ADA_POSOLD)`).
- **Notes:** Companion to `ADA_POS` (multi-assay form); see that entry's Notes for the decomposition rationale. The "newer" vs "older" assay split is paper-specific (Suri 2018 newer assay: sensitivity 23.573 ng/mL, drug tolerance 25 ug/mL; older assay: sensitivity 4 ng/mL, drug tolerance 3,125 ng/mL). Time-varying once positive. Ratified canonically on 2026-04-28.

### ADA_MISSING (**canonical for ADA-result-missing indicator**)
- **Description:** 1 = ADA value is missing (subject did not have a measured ADA result, distinct from a measured negative); 0 = ADA result is reported (positive or negative).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (ADA result reported). Mutually exclusive with `ADA_POS` and `ADA_POSOLD`.
- **Source aliases:**
  - `ATAMISSING` -- used in `Suri_2018_brentuximab.R`.
- **Example models:** `Suri_2018_brentuximab.R` (multiplicative additive effect on ADC clearance: `cl *= (1 + 0.192 * ADA_MISSING)`).
- **Notes:** Used when a substantial fraction of the pooled cohort has no ADA measurement (Suri 2018: 205 of 380 patients) and the modeler retains ADA-missing as a separate level rather than collapsing missingness onto the ADA-negative reference. The non-zero positive estimate of `e_adam_adc_cl` indicates ADA-missing patients are not exchangeable with the ADA-negative reference -- interpret with caution given the missingness mechanism is not random. Ratified canonically on 2026-04-28.

### ADA_TITER (**canonical for continuous antidrug-antibody titer/titre**)
- **Description:** Continuous antidrug-antibody titer/titre (time-varying; matched in time to the PK sample). Covers both the British-spelling reciprocal-dilution convention (`ADA_TITRE`, with `ADA_TITRE = 1` for ADA-negative so `log_e(1) = 0` cancels a log-linear effect) and the American-spelling linear-titer convention (`ADA_TITER`, with `ADA_TITER = 0` for ADA-negative). The per-model `covariateData[[ADA_TITER]]$description` and `notes` must state which zero-encoding convention is in force so the covariate column cannot be misinterpreted.
- **Units:** Reciprocal dilution (e.g., 10, 20, 40, ..., 2560) OR assay units (log2 or arbitrary) -- document per-model in `covariateData[[ADA_TITER]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- ADA-negative encoded per-model (see zero-encoding note).
- **Source aliases:**
  - `ADA_TITRE` -- British spelling (reciprocal-dilution convention; `1` for negative).
  - `ADA titre` -- British spelling long form.
  - `ADAT` -- used in `Moein_2022_etrolizumab.R` (American linear-titer convention; `0` for negative).
- **Example models:** `Jackson_2022_ixekizumab.R` (reciprocal-dilution reference convention with `ADA_TITER = 1` for negatives and `(1 + coef * log_e(ADA_TITER))` on CL), `Moein_2022_etrolizumab.R` (linear-titer convention with `ADA_TITER = 0` for negatives and `exp(theta * ADA_TITER)` on CL, per-unit-titer theta = 0.0365), `Robbie_2012_palivizumab.R` (reciprocal-dilution values 0/10/20/40/>=80 with category-specific multiplicative effects per titer bin; 0 = ADA negative reference).
- **Notes:** The prior separate `ADA_TITRE` (British, `1` = negative) and `ADA_TITER` (American, `0` = negative) canonicals were merged on 2026-04-20 into a single general-scope `ADA_TITER`. The zero-encoding convention is the load-bearing semantic and must be documented per-model. Distinct from `ADA_POS` (binary presence/absence); when the paper reports both, the final model usually keeps only one. Imputation rules (LOCF / NOCB / baseline-as-negative) should be documented per-model.

## Disease / treatment history

### HCT_COND_RIC (**canonical for reduced-intensity conditioning regimen indicator**)
- **Description:** 1 = subject received reduced-intensity conditioning (RIC) chemotherapy prior to allogeneic hematopoietic cell transplantation, 0 = subject received myeloablative conditioning (MAC). Conditioning intensity is fixed per subject for the analysis window (the conditioning regimen was completed before transplantation, before any of the post-transplant tacrolimus PK observations). RIC regimens use lower-dose chemotherapy / radiotherapy to preserve some host haematopoiesis and rely on graft-versus-tumour effect for cytoreduction; MAC regimens deliver high-dose chemotherapy and / or total-body irradiation that fully ablates host marrow.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (myeloablative conditioning, MAC).
- **Source aliases:**
  - `RIC` (Dunlap 2025 NM-TRAN convention; binary 0 / 1) -- used directly in `Dunlap_2025_tacrolimus.R`.
- **Example models:** `Dunlap_2025_tacrolimus.R` (Dunlap 2025 Table 2 reduced-covariate-model column; exponential effect on apparent oral clearance: `cl *= 0.63 ^ HCT_COND_RIC`, so RIC recipients have ~37% lower apparent oral tacrolimus clearance than MAC recipients).
- **Notes:** Conditioning regimen intensity has been reported to associate with post-transplant tacrolimus apparent clearance, likely via gut / hepatic CYP3A activity, GVHD-related inflammatory response, and post-transplant haematopoietic state. The paper-specific definition of "RIC" follows the source publication's own classification (e.g., Dunlap 2025 follows the institutional protocol at UNCMC, which pools non-myeloablative conditioning regimens into the RIC category when assigning the binary indicator); document the source paper's RIC criteria in `covariateData[[HCT_COND_RIC]]$notes`. When a future paper distinguishes a third intensity tier (non-myeloablative, NMA) as a separate covariate level rather than pooling NMA into RIC, register a parallel canonical (e.g. `HCT_COND_NMA`) instead of overloading `HCT_COND_RIC`. Scope: specific because the column is meaningful only for allo-HCT recipients. Ratified canonically on 2026-05-09 alongside the Dunlap 2025 tacrolimus extraction.

### DISEXT_EP (**canonical for extensive colitis / pancolitis indicator**)
- **Description:** 1 = extensive colitis or pancolitis disease extension, 0 = otherwise (any non-extensive disease extension, e.g. left-sided colitis or proctitis).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0. In papers that decompose the disease-extension categorical into both `DISEXT_EP` and `DISEXT_OTHER`, `DISEXT_EP = 0 AND DISEXT_OTHER = 0` corresponds to the left-sided-colitis reference group; in papers that use a single binary indicator for extensive colitis, the reference is pooled non-extensive (left-sided + any other extension).
- **Source aliases:**
  - `EXTCOL` -- used in `Faelens_2021_infliximab.R` (binary 0/1 for extensive colitis at baseline; no separate "other" category).
  - Derived from a multi-level `DISEXT` column in the source (levels: left-sided colitis, extensive/pancolitis, other): `DISEXT_EP = as.integer(DISEXT == "extensive/pancolitis")`.
- **Example models:** `Moein_2022_etrolizumab.R` (paired with `DISEXT_OTHER`; multiplicative effect on CL, +8.2% vs. left-sided colitis), `Faelens_2021_infliximab.R` (single-binary encoding; multiplicative fold-change on V of 1.25 when DISEXT_EP = 1).
- **Notes:** Optionally paired with `DISEXT_OTHER` when the source paper decomposes a three-level disease-extension categorical (left-sided / extensive-pancolitis / other) into two indicators; the pairing is paper-specific and not required. Promoted from scope: specific to scope: general on 2026-04-27 because the binary "extensive colitis vs not" semantics generalize across UC popPK papers regardless of whether the original dataset additionally distinguished an "other" disease-extension category.

### DISEXT_OTHER (**canonical for 'other disease extension' indicator**)
- **Description:** 1 = disease extension other than left-sided colitis or extensive/pancolitis, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (left-sided colitis, when paired with `DISEXT_EP = 0`).
- **Source aliases:** Derived from a multi-level `DISEXT` column: `DISEXT_OTHER = as.integer(DISEXT == "other")`.
- **Example models:** `Moein_2022_etrolizumab.R` (multiplicative effect on CL, +18% vs. left-sided colitis; large uncertainty due to 2% prevalence).
- **Notes:** Paired with `DISEXT_EP`; together they encode the three-level disease-extension categorical.

### PRIOR_TAXANE (**canonical for binary prior-taxane chemotherapy indicator**)
- **Description:** 1 = subject received any prior taxane regimen (docetaxel, paclitaxel, cabazitaxel, etc.) before study entry, 0 = taxane-naive. Time-invariant within a subject (records treatment history at baseline). Oncology-pretreatment indicator: relevant for cohorts in which prior taxane exposure plausibly alters disease biology (e.g., advanced / castration-resistant prostate cancer, where taxane pretreatment is associated with more advanced disease and selects for taxane-resistant clones).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (taxane-naive).
- **Source aliases:**
  - `PTAX` -- used in `vanHasselt_2015_eribulin.R` (van Hasselt 2015 paper notation; binary 0 / 1 for prior docetaxel pretreatment).
- **Example models:** `vanHasselt_2015_eribulin.R` (multiplicative effect on baseline serum PSA: `psa0 = exp(lpsa0 + etalpsa0) * e_prior_taxane_psa0^PRIOR_TAXANE` with `e_prior_taxane_psa0 = 3.23` -- prior-taxane patients have ~3.2x higher baseline PSA than taxane-naive patients, consistent with more advanced disease at study entry).
- **Notes:** Pairs with `PRIOR_TAXANE_DAYS` when a continuous duration-of-pretreatment is also relevant (e.g., van Hasselt 2015 uses both: PTAX on PSA0 and NTRT on KD). Distinct from `CONMED_*` (which is concomitant medication during the study, not pretreatment history) and from generic `PRIOR_*` chemotherapy indicators (which would warrant a separate canonical when a paper differentiates by drug class rather than collapsing to taxanes). Scope: specific because the population semantics (CRPC) and the "any prior taxane" pooling are tied to van Hasselt 2015; future papers that distinguish per-drug pretreatment (docetaxel vs paclitaxel vs cabazitaxel separately) should register parallel canonicals rather than overloading `PRIOR_TAXANE`.

### PRIOR_TAXANE_DAYS (**canonical for cumulative days of prior taxane treatment**)
- **Description:** Cumulative number of days of prior taxane chemotherapy at study entry. 0 for taxane-naive patients (i.e., for any subject with `PRIOR_TAXANE = 0`). Time-invariant within a subject (records the pretreatment history at baseline).
- **Units:** days
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- normalised by a population reference value (van Hasselt 2015 uses the population median of 720 days among prior-taxane-pretreated patients) and entered as a power covariate `(1 + PRIOR_TAXANE_DAYS / 720) ^ e_prior_taxane_days_<param>`. The `+1` inside the bracket makes the covariate effect collapse to a multiplier of 1 for taxane-naive patients (PRIOR_TAXANE_DAYS = 0) regardless of the estimated exponent.
- **Source aliases:**
  - `NTRT` -- used in `vanHasselt_2015_eribulin.R` (van Hasselt 2015 paper notation; cumulative number of days of prior taxane treatment).
- **Example models:** `vanHasselt_2015_eribulin.R` (van Hasselt 2015 Eq. 4 power covariate on drug PSA inhibition rate KD0: `kd0 = exp(lkd0 + etalkd0) * (1 + PRIOR_TAXANE_DAYS / 720)^e_prior_taxane_days_kd0` with `e_prior_taxane_days_kd0 = -4.00` -- KD0 decreases with longer prior-taxane exposure, encoding cross-resistance between docetaxel and eribulin via the shared microtubule-inhibition mechanism).
- **Notes:** Pairs with `PRIOR_TAXANE` (binary). The paper also considered an alternative continuous parameterisation in cycles of prior taxane (`NCYCL`, median 30 cycles) which was deemed slightly less informative (dOFV = -8 for NCYCL vs -10 for NTRT) and was not retained in the final model. If a future model needs the cycle-count form, register a parallel canonical (e.g. `PRIOR_TAXANE_CYCLES`) rather than overloading this one. Scope: specific because the 720-day normalisation reference is tied to the van Hasselt 2015 study population (post-docetaxel mCRPC patients).

## Hypercholesterolemia biomarkers

### PCSK9 (**canonical for baseline unbound serum PCSK9 concentration**)
- **Description:** Baseline unbound serum proprotein convertase subtilisin/kexin type 9 (PCSK9) concentration.
- **Units:** ng/mL (document per-model in `covariateData[[PCSK9]]$units` if a different unit -- typically nM -- is used in a given model; conversion uses a PCSK9 molecular weight of ~72 kDa, so 1 nM ~= 72 ng/mL).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(PCSK9 / ref)^exponent`. Reference values observed: 425 ng/mL (= 5.9 nM) in `Kuchimanchi_2018_evolocumab.R` (population median).
- **Source aliases:** none known.
- **Example models:** `Kuchimanchi_2018_evolocumab.R` (power exponent 0.194 on Vmax: `Vmax * (PCSK9/425)^0.194`).
- **Notes:** PCSK9 is the pharmacological target of anti-PCSK9 monoclonal antibodies (evolocumab, alirocumab, etc.); baseline PCSK9 drives the magnitude of target-mediated elimination and is a recurring covariate in anti-PCSK9 popPK models. Baseline (time-fixed) covariate; patients with missing baseline PCSK9 are typically excluded from analyses that include PCSK9 as a covariate.

## Pharmacogenomic SNPs

**Canonical pattern: `SNP_<GENE>_<RSID>`.** Use one binary indicator per SNP genotype that the source paper tests as a model covariate. The `SNP_` prefix makes the category unambiguous; the gene symbol disambiguates rsIDs grouped by gene; the rsID provides a globally unique identifier. Encoding follows the most common pharmacogenomic convention (also used by Papachristos 2020): `1` = mutant allele present (heterozygous or homozygous mutant); `0` = homozygous wild-type. When a paper uses a different encoding (e.g., per-allele dosage `0/1/2`, dominant model with mutant homozygotes only, or recessive model), document the encoding explicitly in `covariateData[[<COL>]]$notes` and consider registering a separate canonical name. SNP indicators default to scope: specific because the parameter on which they act and the encoded reference category are tied to the source paper's analysis plan; promote to general when a second paper ratifies identical semantics.

### SNP_ICAM1_RS1799969 (**canonical for ICAM-1 rs1799969 mutant indicator**)
- **Description:** Binary genotype indicator for the *ICAM1* rs1799969 single-nucleotide polymorphism (G > A; Gly241Arg / K469E in some references). 1 = at least one mutant (A) allele present (heterozygous or homozygous mutant); 0 = homozygous wild-type (GG).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (homozygous wild-type).
- **Source aliases:**
  - `cat` -- Papachristos 2020 (the paper writes the indicator as `cat` in the CL covariate equation; no formal column name is given in the published narrative).
- **Example models:** `Papachristos_2020_bevacizumab_pk.R`, `Papachristos_2020_bevacizumab_qss.R`, `Papachristos_2020_bevacizumab_pkpd.R` (multiplicative effect on bevacizumab CL: `CL * exp(-0.423 * SNP_ICAM1_RS1799969)` in the PK and PK/PD models; `CL * exp(-0.33 * SNP_ICAM1_RS1799969)` in the binding QSS model -- mutant carriers have lower CL and higher trough levels).
- **Notes:** Time-fixed per subject. Mutant carrier rate in the Papachristos 2020 mCRC cohort: 20% (`Notes` Table 1 of the paper). The biological mechanism by which the *ICAM1* mutant slows bevacizumab clearance is unknown; the association is empirical and may be specific to mCRC.

### SNP_VEGFA_RS1570360 (**canonical for VEGF-A rs1570360 mutant indicator**)
- **Description:** Binary genotype indicator for the *VEGFA* rs1570360 single-nucleotide polymorphism (-1154 G > A; promoter region). 1 = at least one mutant (A) allele present; 0 = homozygous wild-type (GG).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (homozygous wild-type).
- **Source aliases:**
  - `cat1` -- Papachristos 2020 (used as the first categorical indicator in the inter-compartmental clearance equation of the PK model; no formal column name in the narrative).
- **Example models:** `Papachristos_2020_bevacizumab_pk.R` (multiplicative effect on bevacizumab Q: `Q * exp(0.378 * SNP_VEGFA_RS1570360)` -- mutant carriers have higher inter-compartmental clearance).
- **Notes:** Time-fixed per subject. Mutant carrier rate in the Papachristos 2020 mCRC cohort: 33%. The covariate is significant in the standalone PK model but does not appear in the binding QSS or PK/PD models because in those models the inter-compartmental clearance covariate effects are absorbed into the rs699947 effect on Q (PK/PD) or into the K_ss / BM0 effects (QSS).

### SNP_VEGFA_RS699947 (**canonical for VEGF-A rs699947 mutant indicator**)
- **Description:** Binary genotype indicator for the *VEGFA* rs699947 single-nucleotide polymorphism (-2578 C > A; promoter region). 1 = at least one mutant (A) allele present; 0 = homozygous wild-type (CC).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (homozygous wild-type).
- **Source aliases:**
  - `cat2` -- Papachristos 2020 PK model (second categorical indicator on Q).
  - `cat` -- Papachristos 2020 binding QSS model (effect on K_ss and BM0) and PK/PD model (effect on Q).
- **Example models:** `Papachristos_2020_bevacizumab_pk.R` (effect on Q: -0.429), `Papachristos_2020_bevacizumab_qss.R` (effect on K_ss: +1.22, on BM0: -0.851), `Papachristos_2020_bevacizumab_pkpd.R` (effect on Q: -0.414).
- **Notes:** Time-fixed per subject. Mutant carrier rate in the Papachristos 2020 mCRC cohort: 52%. The mutant allele is associated with lower baseline free VEGF-A levels and a higher in-vivo affinity (higher K_ss), consistent with reports that rs699947 mutants have prolonged overall survival on bevacizumab-based therapy.

### SNP_SLCO1B1_RS11045819 (**canonical for SLCO1B1 rs11045819 mutant indicator**)
- **Description:** Binary genotype indicator for the *SLCO1B1* rs11045819 single-nucleotide polymorphism (C > A; OATP1B1 transporter, exon 5, P155T). 1 = at least one mutant (A) allele present (heterozygous AC or homozygous AA); 0 = homozygous wild-type (CC). The Hennig 2015 cohort (n = 35 successfully genotyped of 44) reported 5 AC heterozygotes, 30 CC homozygotes, and 0 AA homozygotes, so the indicator is effectively heterozygous-vs-CC in that cohort.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (homozygous wild-type CC).
- **Source aliases:**
  - `SLCO1B1 rs11045819 genotype` -- Hennig 2015 (paper text; the source NONMEM control stream is in the unrecovered AAC supplement, so the formal column name is not on disk).
- **Example models:** `Hennig_2015_rifabutin.R` (multiplicative effect on rifabutin bioavailability F: `F * (1 + 0.304 * SNP_SLCO1B1_RS11045819)` -- AC carriers have ~30% higher rifabutin F than CC reference; dOFV = -6.5).
- **Notes:** Time-fixed per subject. Carrier rate in the Hennig 2015 South-African HIV/TB cohort: 14% (5 of 35 genotyped). SLCO1B1 encodes OATP1B1, a hepatic uptake transporter; rs11045819 has been associated with reduced rifampicin and lopinavir concentrations in prior studies but in Hennig 2015 was associated with INCREASED rifabutin bioavailability (note opposite direction of effect across rifamycins).

### SLCO1B1_HAP15_HET (**canonical for SLCO1B1*15 haplotype heterozygote indicator**)
- **Description:** Binary haplotype indicator for the *SLCO1B1* `*15` reduced-function haplotype (the cis combination of the 388A>G / rs2306283 and 521T>C / rs4149056 variants; encodes the OATP1B1 N130D + V174A double mutant). 1 = subject carries exactly one *15 allele (heterozygous: `*1a/*15` or `*1b/*15`), 0 = otherwise (the union of *15-noncarriers and *15-homozygotes; the paired indicator `SLCO1B1_HAP15_HOM` flags the homozygous group). Time-fixed per subject (germline haplotype).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (homozygous *1a/*1a, *1a/*1b, or *1b/*1b -- i.e., no *15 allele). The reference group is the union of the three non-*15 diplotypes; `SLCO1B1_HAP15_HOM` flags the *15/*15 group.
- **Source aliases:**
  - `HT` -- Ide 2009 (paper text Eq. for `Frel = 1 * theta1^HT * theta2^HM` where `HT = 1` for heterozygotes `*1a/*15` and `*1b/*15`).
- **Example models:** `Ide_2009_pravastatin.R` (multiplicative effect on relative bioavailability Frel: `Frel = 1.50^SLCO1B1_HAP15_HET * 1.95^SLCO1B1_HAP15_HOM` -- *15 heterozygotes have 50% higher Frel than *15-noncarriers; dOFV = 32.2 in backward elimination, p < 0.001).
- **Notes:** Paired with `SLCO1B1_HAP15_HOM` to encode a three-level haplotype categorical (noncarrier / heterozygote / homozygote) with `*15`-noncarrier as the implicit reference (both indicators = 0). Distinct from the SNP-level canonical `SNP_SLCO1B1_RS11045819` (which encodes only the C>A variant at a different position; rs11045819 = P155T) and from the *15 component SNPs rs2306283 (388A>G) and rs4149056 (521T>C) individually: future Ide-style extractions that pool *5 (521T>C only) with *15 should still record their values under this canonical and document the pooling rule in `covariateData[[SLCO1B1_HAP15_HET]]$notes`. Distribution in the Ide 2009 cohort of 57 healthy Japanese male volunteers (Table I): 28 noncarriers, 23 heterozygotes, 6 homozygotes. Ratified canonically on 2026-05-12 alongside the Ide 2009 extraction.

### SLCO1B1_HAP15_HOM (**canonical for SLCO1B1*15 haplotype homozygote indicator**)
- **Description:** Binary haplotype indicator for the *SLCO1B1* `*15` reduced-function haplotype. 1 = subject carries two *15 alleles (homozygous: `*15/*15`), 0 = otherwise (the union of *15-noncarriers and *15-heterozygotes; the paired indicator `SLCO1B1_HAP15_HET` flags the heterozygous group). Time-fixed per subject (germline haplotype).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (homozygous *1a/*1a, *1a/*1b, or *1b/*1b -- i.e., no *15 allele). The reference group is the union of the three non-*15 diplotypes; `SLCO1B1_HAP15_HET` flags the *15-heterozygote group.
- **Source aliases:**
  - `HM` -- Ide 2009 (paper text Eq. for `Frel = 1 * theta1^HT * theta2^HM` where `HM = 1` for homozygotes `*15/*15`).
- **Example models:** `Ide_2009_pravastatin.R` (multiplicative effect on relative bioavailability Frel: `Frel = 1.50^SLCO1B1_HAP15_HET * 1.95^SLCO1B1_HAP15_HOM` -- *15 homozygotes have 95% higher Frel than *15-noncarriers; dOFV = 33.7 in backward elimination, p < 0.001).
- **Notes:** Paired with `SLCO1B1_HAP15_HET` to encode a three-level haplotype categorical (noncarrier / heterozygote / homozygote) with `*15`-noncarrier as the implicit reference (both indicators = 0). See `SLCO1B1_HAP15_HET` Notes for the broader context; population distribution in Ide 2009 was 6 of 57 (10.5%) homozygotes. Ratified canonically on 2026-05-12 alongside the Ide 2009 extraction.

### ABCB1_HAP_TTT (**canonical for ABCB1 TTT haplotype carrier indicator**)
- **Description:** Binary haplotype indicator for the *ABCB1* `TTT` haplotype across the rs1128503 (1236C>T, exon 12, synonymous Gly412Gly) / rs2032582 (2677G>T/A, exon 21, Ala893Ser/Thr) / rs1045642 (3435C>T, exon 26, synonymous Ile1145Ile) SNP block. 1 = subject carries at least one `TTT` haplotype (heterozygous or homozygous; pooled because the homozygote frequency was < 0.1 in the de Wit 2016 cohort); 0 = no `TTT` haplotype. Time-fixed per subject (germline haplotype).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no `TTT` haplotype).
- **Source aliases:**
  - `ABCB1 TTT haplotype` -- de Wit 2016 (paper text Methods 'Pharmacogenetic analysis' and Table 2 final-model `theta TTT on F` row; haplotypes phased in gPLINK with certainty > 0.97).
- **Example models:** `deWit_2016_everolimus.R` (multiplicative effect on apparent bioavailability F: `F = 1 * 0.792^ABCB1_HAP_TTT` -- carriers have 20.8% lower F than non-carriers; dOFV = 9.6 in backward elimination, P < 0.01).
- **Notes:** The `TTT` haplotype of *ABCB1* (P-glycoprotein, MDR1 efflux transporter) is associated with enhanced P-gp efflux activity and reduced everolimus bioavailability (de Wit 2016 Discussion paragraph 5); de Wit cites prior evidence that the same `TTT` haplotype also reduces exposure / efficacy of other P-gp substrates [refs 20-22 in the paper], although directionally inconsistent results have been reported [refs 7, 23, 24]. Het and hom carriers were pooled in de Wit 2016 because the homozygote frequency was < 0.1 (Methods 'Pharmacogenetic analysis'); future extractions that estimate separate het / hom effects should register paired `ABCB1_HAP_TTT_HET` and `ABCB1_HAP_TTT_HOM` indicators following the `SLCO1B1_HAP15_HET` / `SLCO1B1_HAP15_HOM` precedent above. Distinct from any individual ABCB1 SNP indicator (rs1128503, rs2032582, rs1045642 alone) because the haplotype is the cis combination tested jointly. Ratified canonically on 2026-05-16 alongside the de Wit 2016 everolimus extraction.

## Lifestyle / medical history

### SMOKE (**canonical for current-smoker binary indicator**)
- **Description:** 1 = current smoker at baseline, 0 = non-smoker.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (non-smoker).
- **Source aliases:**
  - `Smoking` (case-insensitive) -- used in `Ma_2020_sarilumab_anc.R`.
- **Example models:** `Ma_2020_sarilumab_anc.R` (power-form on baseline ANC: `BASE * 1.15^SMOKE`).
- **Notes:** Baseline-only indicator; does not track within-study smoking-cessation changes. Use this two-level (current vs non-smoker) encoding when the source paper does not split former and never smokers. When the source uses a 3-level smoking-status categorical (never / former / current), use the paired `SMOKE_CURRENT` + `SMOKE_NEVER` indicators below instead -- the 3-level encoding cannot be reduced to a single `SMOKE` column without losing information.

### SMOKE_CURRENT (**canonical for current-smoker indicator (paired with SMOKE_NEVER)**)
- **Description:** 1 = current smoker at baseline, 0 otherwise (former or never smoker). Paired with `SMOKE_NEVER` to encode a 3-level smoking-status categorical with former smoker as the implicit reference (both indicators = 0).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (former smoker, when paired with `SMOKE_NEVER` = 0). The pairing follows the `RACE_<GROUP>` convention for paired indicators.
- **Source aliases:**
  - `Smoking status = Current` / `SMOK = 2` (case-insensitive) -- derived from a 3-level smoking-status column.
- **Example models:** `Hwang_2023_monalizumab.R` (proportional-shift effect on V1: `(1 + 0.0484)^SMOKE_CURRENT`; reference category former smoker).
- **Notes:** Baseline-only indicator. See also `SMOKE_NEVER` (paired indicator) and `SMOKE` (binary current-vs-non-smoker encoding when the source paper does not split former vs never).

### SMOKE_NEVER (**canonical for never-smoker indicator (paired with SMOKE_CURRENT)**)
- **Description:** 1 = never smoker at baseline, 0 otherwise (former or current smoker). Paired with `SMOKE_CURRENT` to encode a 3-level smoking-status categorical with former smoker as the implicit reference (both indicators = 0).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (former smoker, when paired with `SMOKE_CURRENT` = 0). The pairing follows the `RACE_<GROUP>` convention for paired indicators.
- **Source aliases:**
  - `Smoking status = Never` / `SMOK = 0` (case-insensitive) -- derived from a 3-level smoking-status column.
- **Example models:** `Hwang_2023_monalizumab.R` (proportional-shift effect on V1: `(1 - 0.141)^SMOKE_NEVER`; reference category former smoker).
- **Notes:** Baseline-only indicator. See also `SMOKE_CURRENT` (paired indicator) and `SMOKE` (binary current-vs-non-smoker encoding when the source paper does not split former vs never).

## Formulation / assay / study

### ROUTE_IV (**canonical for IV-vs-SC administration route indicator**)
- **Description:** 1 = subject received intravenous (IV) administration, 0 = subcutaneous (SC) administration. Per-subject (study-fixed) covariate flagging the dosing route when a population analysis pools cohorts that differ by route, with covariate effects on PK parameters that capture route-specific disposition behaviour.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (SC).
- **Source aliases:**
  - "Admin route = IV" (categorical effect column in Yu 2022 covariate equations).
  - `IV` -- used in `Zierhut_2008_osteoprotegerin.R` (DDMODEL00000233 `dataObj` column flagging IV vs SC cohort, switching the PK observation residual SD between `CcpropSdIV` and `CcpropSdSC`).
  - "IV" / "SC arm" -- used in `Wang_2021_pertuzumab.R` (per-subject FeDeriCa arm indicator P+H IV vs PH FDC SC, switching the proportional residual SD between `CcpropSdIv` and `CcpropSdSc`).
- **Example models:**
  - `Yu_2022_ofatumumab.R` (exponential effect on R0, CL, Q, ksyninf).
  - `Zierhut_2008_osteoprotegerin.R` (per-subject indicator switching the PK observation residual SD between the IV cohort (`CcpropSdIV`) and the SC cohort (`CcpropSdSC`)).
  - `Wang_2021_pertuzumab.R` (per-subject indicator switching the proportional residual SD between the IV (`CcpropSdIv` = 0.175) and SC (`CcpropSdSc` = 0.155) cohorts of the FeDeriCa popPK).
  - `Fiedler-Kelly_2019_fremanezumab.R` (per-subject indicator switching both the central volume of distribution (Vc,IV = 2.98 L FIXED vs Vc,SC = 1.88 L) and the residual-error structure (IV: proportional-only with SD = sqrt(0.0467) = 0.21610; SC: combined additive sqrt(0.204) + proportional sqrt(0.0531)) in the pooled phase 1/2b/3 fremanezumab popPK).
- **Notes:** This is the per-subject covariate-equation indicator, distinct from the dosing-event `cmt` column that names the target compartment. When simulating, set `ROUTE_IV = 1` for IV cohorts and dose into the central compartment; set `ROUTE_IV = 0` for SC cohorts and dose into the depot. Scope: specific because the set of parameters that differ by route is paper-specific (Yu 2022 carries route-specific exponential effects on disposition parameters; Zierhut 2008, Wang 2021, and Fiedler-Kelly 2019 carry route-specific PK observation residual SDs; Fiedler-Kelly 2019 additionally carries a route-specific central volume of distribution).

### DEVICE_AI (**canonical for autoinjector-vs-prefilled-syringe SC device indicator**)
- **Description:** 1 = subject's SC dose delivered via autoinjector (AI), 0 = prefilled syringe (PFS). Per-subject (study-fixed) covariate flagging the SC delivery device when a model carries device-specific PK effects.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (PFS).
- **Source aliases:** "Formulation = AI" (categorical effect column in Yu 2022 covariate equations).
- **Example models:** `Yu_2022_ofatumumab.R` (exponential effect on k_e(P) and R0), `Diep_2026_donidalorsen.R` (Phoenix linear-effect `(1 + e_device_ai_ka * DEVICE_AI)` on the typical SC absorption rate constant with theta = +0.262 -> multiplier 1.262 for autoinjector vs vial-and-syringe reference; characterized in the ISIS 721744-CS9 single-dose bioequivalence cohort).
- **Notes:** Set to 0 (PFS / vial reference) for IV subjects, since the device is undefined for IV; the IV-specific effects are captured by `ROUTE_IV` instead. Scope: specific because the AI / vial / PFS contrast and which parameters it affects depend on the study's device-comparison design.

### INJSITE_ARM (**canonical for SC injection-site = arm indicator**)
- **Description:** 1 = subject's SC dose injected into the arm, 0 = abdomen (the universal SC reference site across the popPK literature). Per-dose-record covariate flagging the SC injection site when a population analysis estimates site-specific absorption parameters.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (abdomen).
- **Source aliases:** paper narrative "arm" / "abdomen" subgroup labels driving site-specific ka in Diep 2022.
- **Example models:** `Diep_2022_eplontersen.R` (additive log-shift `e_injsite_arm_ka = log(ka_arm / ka_ab)` on the typical absorption rate constant: ka_arm = 0.217 1/h vs ka_ab = 0.282 1/h, ~30% higher ka for abdomen; INJSITE_ARM = 1 selects the arm typical value), `Diep_2026_donidalorsen.R` (Phoenix linear-effect `(1 + e_injsite_arm_ka * INJSITE_ARM)` on the typical absorption rate constant with theta = -0.338 -> multiplier 0.662 for arm; the paper's reference category is "abdomen or thigh" rather than "abdomen" alone, but is consistent with the canonical reference because abdomen is the universal SC reference site and the thigh effect is pooled into the reference category by the Diep 2026 model).
- **Notes:** Specific scope because the arm-vs-abdomen contrast is paper-specific. Sister canonical to `INJSITE_THIGH` (thigh-vs-abdomen indicator anticipated for future SC-route models with thigh-specific absorption). Per-administration rather than per-subject -- a subject in a multi-dose simulation can switch SC injection sites between doses; supply the indicator on each dose record. Distinct from `ROUTE_IV` (IV vs SC route, not within-SC site) and from `DEVICE_AI` (autoinjector vs prefilled syringe, device rather than anatomical site).

### STUDY_APLIOS (**canonical for APLIOS bioequivalence study indicator**)
- **Description:** 1 = subject enrolled in the APLIOS bioequivalence study (NCT03560739; phase 2; ofatumumab AI vs PFS in RMS), 0 = other study in the Yu 2022 pooled analysis.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-APLIOS studies: OMS115102, MIRROR, ASCLEPIOS I, ASCLEPIOS II).
- **Source aliases:** "Study = APLIOS" (categorical effect column in Yu 2022 covariate equations).
- **Example models:** `Yu_2022_ofatumumab.R` (exponential effect on Emax of B cell lysis).
- **Notes:** Captures a between-study shift in the maximum B-cell lysis stimulatory effect not explained by the other covariates in the final model.

### STUDY_MIRROR (**canonical for MIRROR dose-finding study indicator**)
- **Description:** 1 = subject enrolled in the MIRROR dose-finding study (NCT01457924; phase 2; SC ofatumumab dose-ranging in RRMS), 0 = other study in the Yu 2022 pooled analysis.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-MIRROR studies: OMS115102, APLIOS, ASCLEPIOS I, ASCLEPIOS II).
- **Source aliases:** "Study = MIRROR" (categorical effect column in Yu 2022 covariate equations).
- **Example models:** `Yu_2022_ofatumumab.R` (exponential effect on B cell elimination rate kout).
- **Notes:** Captures a between-study shift in the B cell elimination rate not explained by the other covariates in the final model.

### REGI_BID (**canonical for twice-daily dosing-regimen indicator**)
- **Description:** 1 = subject's dosing regimen is BID (twice daily), 0 = QD (once daily) or other non-BID regimen. Per-subject (regimen-fixed) categorical indicator for population analyses that pool QD and BID arms and test regimen as a covariate.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-BID regimen -- typically QD).
- **Source aliases:**
  - `BID` -- used in `Girard_2012_pimasertib.R`.
- **Example models:** `Girard_2012_pimasertib.R` (additive shift on the cumulative-logit AE-score model: `theta_bid * REGI_BID`; -0.399 logit units for BID vs QD).
- **Notes:** Specific scope because the QD-vs-BID contrast is study-specific; future regimen-comparison models that contrast different schedules should either extend this entry's example list (when QD is the reference) or register a sibling indicator (`REGI_TID`, `REGI_QW`) following the same pattern. Distinct from `DOSE` (dose level in mg) and from total-daily-dose aggregates: a 60 mg/day cohort can include either a 60 mg QD subgroup or a 30 mg BID subgroup, and both share the same `DOSE = 60` while differing in `REGI_BID`.

### FED (**canonical for fed-vs-fasted dose-record indicator**)
- **Description:** 1 = fed state at dosing, 0 = fasted.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (fasted).
- **Example models:** `Kyhl_2016_nalmefene.R`, `Goel_2016_Sonidegib.R` (the Goel 2016 healthy-fasted F effect is applied via `e_healthy_fast_f ^ (DIS_HEALTHY * (1 - FED))`; cancer-patient records have FED = 1, healthy-fasted-arm records have FED = 0, high-fat-meal arm records have FED = 1 and the additional FED_HIGHFAT = 1 indicator).

### FED_HIGHFAT (**canonical for high-fat-meal-at-dosing indicator**)
- **Description:** 1 = oral dose administered after a high-fat meal, 0 = oral dose administered under any other meal condition (typically fasted or light meal). Refines the more general `FED` indicator for studies that specifically test the high-fat-meal food effect on bioavailability or absorption rate (a common solubility-limited absorption phenotype).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (non-high-fat meal condition; most often "fasted" or "2 h post light meal" as defined per study protocol). Document the operational reference per model.
- **Source aliases:**
  - `Fatmeal` -- used in `Goel_2016_Sonidegib.R` (covariate on F).
  - `FATM` -- used in `Goel_2016_Sonidegib.R` (covariate on Ka; same indicator as `Fatmeal` in that paper).
- **Example models:** `Goel_2016_Sonidegib.R` (multiplicative effect on F: `5.74^FED_HIGHFAT` -- ~5.7-fold higher F under high-fat meal vs 2 h post-light-meal reference; multiplicative effect on Ka: `1.01^FED_HIGHFAT` -- no meaningful effect).
- **Notes:** Distinct from `FED` (binary fed-vs-fasted): `FED_HIGHFAT` carries the specific "high-fat meal" semantic (typically >= 800 kcal, >= 50% calories from fat per FDA guidance). Document the per-protocol meal definition in `covariateData[[FED_HIGHFAT]]$notes`. When a paper reports a high-fat-meal arm and a separate fasted-healthy arm (e.g., Goel 2016), use `FED_HIGHFAT` for the high-fat semantic and the existing `FED` + `DIS_HEALTHY` indicators for the composite healthy-fasted effect (Goel 2016 applies the e_healthy_fast_f effect via `(DIS_HEALTHY * (1 - FED))`); the retired `HV_FAST` composite indicator was deleted on 2026-05-11.

### MULTI_DOSE_PT (**canonical for multiple-dose-phase-in-patients indicator**)
- **Description:** 1 = dose record from the multiple-dose phase of a clinical-pharmacology study in patients, 0 = otherwise (single-dose run-in records, healthy-volunteer records, or first-dose records in patient studies). Captures any systematic shift in apparent bioavailability between the controlled run-in and the longer multiple-dose phase, typically driven by variable food-restriction compliance over many dosing days.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (single-dose / run-in / healthy-volunteer dose records).
- **Source aliases:**
  - `FMDD` -- used in `Goel_2016_Sonidegib.R` (Goel 2016 covariate on F).
- **Example models:** `Goel_2016_Sonidegib.R` (multiplicative effect on F: `1.16^MULTI_DOSE_PT` -- ~16% higher apparent F during the multiple-dose phase relative to first dose, attributed in the paper to occasional non-fasting compliance).
- **Notes:** Specific scope because the indicator's exact definition (dose-record level vs subject level, run-in inclusion, occasion boundary) is paper-specific. In Goel 2016, the dataset distinguishes the run-in single dose from the daily multiple-dose phase; the indicator switches at the start of the multiple-dose phase for cancer patients. Distinct from `FED` and `FED_HIGHFAT` (which are per-record meal-state indicators) and from `REGI_BID` (regimen indicator). Future models that need a generic "occasion boundary" effect should consider the existing `ooc<n>` IOV pattern instead.

### FORM_TABLET (**canonical for tablet vs solution formulation indicator**)
- **Description:** 1 = tablet formulation, 0 = solution.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (solution).
- **Source aliases:**
  - `TABLET` -- earlier name used by `Kyhl_2016_nalmefene.R` and `Tikiso_2021_abacavir.R`; renamed to `FORM_TABLET` for consistency with the `FORM_*` family (`FORM_CAPSULE`, future `FORM_SUSPENSION`, etc.).
- **Example models:** `Kyhl_2016_nalmefene.R` (additive shift on residual error: tablet vs solution), `Tikiso_2021_abacavir.R` (multiplicative effect on the absorption mean transit time MTT: tablet (abacavir + lamivudine fixed-dose-combination tablet) vs liquid solution; `mtt *= (1 + 0.249 * FORM_TABLET)`, i.e. 24.9% slower absorption for the FDC tablet relative to the abacavir liquid reference).
- **Notes:** Scoped specific because the "tablet vs solution" contrast is tied to formulation-comparison study designs (Kyhl 2016 nalmefene tablet-vs-oral-solution; Tikiso 2021 abacavir + lamivudine FDC tablet vs abacavir liquid). Distinct from the `FORM_FDC` canonical (Wilkins 2008 antitubercular fixed-dose-combination of multiple drugs, contrasted against single-drug tablets) because here the comparator is a non-tablet liquid / solution rather than a separate tablet product. Future formulation-comparison models should either extend this entry's example list when the comparator is a liquid / solution, or register a sibling canonical when contrasting two tablet products.

### FORM_CAPSULE (**canonical for capsule formulation indicator**)
- **Description:** 1 = capsule formulation, 0 = the per-paper comparator non-capsule formulation. The complement formulation is paper-defined: solution for Hennig 2006 / 2007 itraconazole, tablet for Gupta 2016 lenvatinib. Document the comparator and the reference-category bioavailability per-model in `covariateData[[FORM_CAPSULE]]$notes`.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (per-paper non-capsule comparator: solution in Hennig 2006 with `fdepot` fixed to 1 for the solution arm; tablet in Gupta 2016 with F fixed to 1 for the tablet arm).
- **Source aliases:**
  - `PREP` -- used in `Hennig_2006_itraconazole.R` (Clin Pharmacokinet 2006;45(11):1099-1114; PREP = 1 = capsule, PREP = 0 = oral solution) and in `Hennig_2007_itraconazole.R` (Br J Clin Pharmacol 2007;63(4):438-450; DOI 10.1111/j.1365-2125.2006.02778.x; same orientation, capsule typical absorption parameters as the published reference).
  - `CAPSULE` -- earlier name used in the `Hennig_2006_itraconazole.R` and `Hennig_2007_itraconazole.R` model files before the `FORM_*` rename.
- **Example models:** `Hennig_2006_itraconazole.R` (Hennig 2006 Table II final estimates: capsule `ka` 0.09 h^-1 vs solution `ka` 0.96 h^-1 and capsule relative bioavailability 0.55 vs solution 1; the `etalfdepot` IIV applies only to the capsule arm); `Hennig_2007_itraconazole.R` (selects between `lka_cap` and `lka_sol` typical-value absorption rate constants and applies `f(depot) <- (1 - FORM_CAPSULE) + FORM_CAPSULE * fdepot` so the relative bioavailability `F_rel = 0.817` is applied only to the capsule arm); `Gupta_2016_lenvatinib.R` (relative bioavailability of capsule vs tablet is 0.896; F1 fixed to 1 for the tablet reference; the `etalfcap` IIV (30.2% CV) applies only to the capsule arm).
- **Notes:** Scoped specific because the complement reference category is paper-defined (solution for Hennig itraconazole, tablet for Gupta lenvatinib). Sibling to `FORM_TABLET` (Kyhl 2016 / Tikiso 2021 tablet vs solution) under the `FORM_*` family. Future formulation-comparison models that need a capsule indicator should reuse this canonical, extending the example list and documenting the comparator in per-model notes.

### FORM_SUSPENSION (**canonical for tablet-suspended-in-water formulation indicator**)
- **Description:** 1 = subject received the modelled drug as solid oral tablets suspended in a small volume of water at bedside immediately before swallowing (i.e. tablets disintegrated into an extemporaneously-prepared liquid suspension rather than a formally manufactured suspension product); 0 = subject received the same solid oral tablets swallowed whole. Per-dose-occasion (not per-subject) indicator because a single participant can receive both formulations across study occasions in a crossover bioequivalence design.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (tablets swallowed whole; the typical-value MAT reference in Svensson 2018 Table 2).
- **Source aliases:**
  - `FORM` -- used in `Svensson_2018_bedaquiline.R` (paper's narrative "whole vs suspended" with the suspended formulation as the 1 level in the analytical control stream).
- **Example models:** `Svensson_2018_bedaquiline.R` (multiplicative effect on the typical mean absorption time MAT: `mat_typ = exp(lmat) * (1 + e_susp_mat * FORM_SUSPENSION)` with `e_susp_mat = +0.23` (suspended-tablet MAT is 23% longer than whole-tablet MAT; Svensson 2018 Table 2, 95% CI 2.1-48%, P = 0.03); relative bioavailability F is held identical between formulations because the paper found no statistically significant difference (95% nonparametric CI 94-108% within the 80-125% bioequivalence criteria)).
- **Notes:** Specific scope because the bedside-suspension manipulation is tied to the paediatric-dosing motivation of Svensson 2018 (tablets crushed or suspended to enable administration to children before a paediatric formulation is available). Distinct from `FORM_TABLET` (Kyhl 2016 / Tikiso 2021 tablet vs liquid solution) because the contrast here is *the same tablet swallowed two different ways* rather than two different drug products; distinct from `FORM_CAPSULE` and `FORM_POWDER` for the same reason. Future tablet-suspended-in-water bioequivalence extractions (a common paediatric-pharmacokinetics manipulation in tuberculosis, HIV, and oncology paediatric trials) should reuse this canonical and extend the example list. Ratified canonically on 2026-05-16 alongside the Svensson 2018 bedaquiline extraction.

### FORM_POWDER (**canonical for oral powder formulation indicator**)
- **Description:** 1 = subject received the oral powder formulation of the modelled drug; 0 = subject received the comparator solid oral formulation (typically a tablet, but document the per-paper comparator in `covariateData[[FORM_POWDER]]$notes`).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (tablet; F = 1 fixed in Yukawa 1990 Model 2).
- **Source aliases:**
  - `FORM_POWDER` -- used in `Yukawa_1990_phenytoin.R` (paper's `BA` indicator inverted: source `BA = 1` if tablet, 0 if powder; canonical `FORM_POWDER = 1 - BA_indicator` so 0 is the tablet reference).
- **Example models:** `Yukawa_1990_phenytoin.R` (Yukawa 1990 Model 2 dose-dependent powder bioavailability `F_powder = 1 - exp(-9.92 / DOSE_PHT_MGKGD)`; tablet F fixed at 1), `Retlich_2015_linagliptin.R` (multiplicative shift on the linagliptin first-order absorption rate constant Ka: powder-in-bottle Ka = 0.933 1/h vs tablet formulation 2 reference Ka = 0.441 1/h; the tablet formulation 1 comparator is captured by the sibling canonical `FORM_LINAG_TAB1`).
- **Notes:** Specific scope because the "powder vs tablet" contrast is tied to a particular drug-product manufacturing comparison (Yukawa 1990 contrasts Aleviatin brand phenytoin powder with Aleviatin tablets, both from Dainippon Pharmaceutical Co.; Retlich 2015 contrasts an early-phase linagliptin powder-in-bottle formulation against the marketed linagliptin tablet). Mirrors the sibling `FORM_TABLET` (Kyhl 2016 / Tikiso 2021 tablet vs solution) and `FORM_CAPSULE` (Hennig 2006 / Hennig 2007 capsule vs solution) under the `FORM_*` family. Future powder-formulation models should reuse this canonical, extending the example list and documenting the per-paper comparator. Ratified canonically on 2026-05-10 alongside the Yukawa 1990 phenytoin extraction.

### FORM_FDC (**canonical for fixed-dose-combination antitubercular formulation indicator**)
- **Description:** 1 = subject received the rifampicin-containing fixed-dose-combination (FDC) antitubercular product (rifampicin co-formulated with isoniazid, pyrazinamide, and optionally ethambutol in a single tablet); 0 = subject received the same drugs as separate single-drug tablets ("SDC", separate-drug-combination). Per-subject (regimen-fixed) categorical covariate flagging the formulation when a population analysis pools FDC and SDC arms and tests formulation as a covariate on absorption / disposition parameters.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 1 (FDC; the most-common formulation in the Wilkins 2008 cohort and the typical-value reference for `lmtt` / `lcl`).
- **Source aliases:**
  - `FDC` -- used in `Wilkins_2008_rifampicin.R` (DDMODEL00000280 NMTRAN `$INPUT` column; values 0 / 1 with the same orientation as the canonical, 1 = FDC).
- **Example models:** `Wilkins_2008_rifampicin.R` (multiplicative `(1 + e_fdc0_mtt * (1 - FORM_FDC))` shift on MTT and `(1 + e_fdc0_cl * (1 - FORM_FDC))` shift on CL; SDC subjects (FDC = 0) had 104% longer MTT and 23.6% higher CL than the FDC = 1 reference per Wilkins 2008 final estimates).
- **Notes:** Specific scope because the FDC vs SDC contrast is tied to the antitubercular fixed-dose-combination context (rifampicin + isoniazid + pyrazinamide +/- ethambutol). Future antitubercular-FDC models should extend the example list rather than register a new canonical. Distinct from `FORM_TABLET` (Kyhl 2016 nalmefene tablet vs solution) and the rest of the `FORM_*` (drug-product-version) family because the FDC-vs-SDC contrast compares two tablet products, not a tablet vs a non-tablet, and the underlying mechanism is co-formulation-driven absorption rather than drug-product manufacturing.

### RIA_ASSAY (**canonical for radioimmunoassay vs LC-MS/MS bioanalytical method indicator**)
- **Description:** 1 = radioimmunoassay; 0 = LC-MS/MS.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Example models:** `Kyhl_2016_nalmefene.R`.
- **Notes:** Switches the additive residual-error magnitude.

### FORM_NS0 (**canonical for NS0 cell-line formulation indicator**)
- **Description:** 1 = NS0 cell-line formulation, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Example models:** `Zhu_2017_lebrikizumab.R`.

### FORM_CHO_PHASE2 (**canonical for CHO Phase 2 formulation indicator**)
- **Description:** 1 = CHO Phase 2 formulation, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Example models:** `Zhu_2017_lebrikizumab.R`.

### CONMED_EOX (**canonical for concomitant EOX (epirubicin + oxaliplatin + capecitabine) chemotherapy backbone indicator**)
- **Description:** 1 = concomitant epirubicin + oxaliplatin + capecitabine (EOX) chemotherapy backbone, 0 = other backbone (e.g., mFOLFOX6, CAPOX, or single-agent).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-EOX backbone).
- **Source aliases:** `COMB` (used by Yamada 2025 Table 1 with the EOX level coded as the non-reference category; renamed to `CONMED_EOX` to preserve the semantic meaning of the 1-level).
- **Example models:** `Yamada_2025_zolbetuximab.R` (fractional effect on V1).
- **Notes:** Disease-backbone indicator. If a future model needs more backbone categories, encode each as its own indicator (`COMB_CAPOX`, `COMB_FOLFOX`, ...) with a single reference group.

### FORM_P2F2 (**canonical for isatuximab P2F2 drug-material indicator**)
- **Description:** 1 = isatuximab P2F2 drug material (intended commercial / phase III material, used in the EFC14335 / ICARIA-MM study), 0 = P1F1 drug material (early-phase material used in TED10893 / TED14154 / TCD14079).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (P1F1).
- **Source aliases:**
  - `Drug_mat` -- used in `Fau_2020_isatuximab.R`. Values 0 / 1 with the same orientation as the canonical (1 = P2F2 / commercial-bound material).
- **Example models:** `Fau_2020_isatuximab.R` (exponential effect on Vc with coefficient -0.137; P2F2 patients had ~13% lower Vc than P1F1).
- **Notes:** Phase III / commercial-bound formulation indicator for isatuximab; the FORM_* family stays scope-specific per nlmixr2lib policy that drug-product-version indicators are kept model-specific unless they generalize across multiple drugs. Set to 1 to simulate the marketed material.

### FORM_LINAG_TAB1 (**canonical for linagliptin tablet formulation 1 indicator**)
- **Description:** 1 = subject received the linagliptin "tablet formulation 1" (used in Retlich 2015 Study 2), 0 = subject received tablet formulation 2 (the marketed linagliptin tablet, used in Studies 3 and 4) OR the powder-in-bottle formulation (used in Study 1). The powder-vs-tablet contrast is captured by the sibling canonical `FORM_POWDER`; this indicator switches between the two tablet formulations.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (tablet formulation 2 = marketed linagliptin tablet, the typical-value Ka reference; or `FORM_POWDER = 1` for the powder).
- **Source aliases:** none known.
- **Example models:** `Retlich_2015_linagliptin.R` (multiplicative shift on the linagliptin first-order absorption rate constant Ka; typical Ka = 0.441 1/h for tablet formulation 2 (reference), 0.795 1/h for tablet formulation 1, 0.933 1/h for the powder formulation).
- **Notes:** Specific scope because the linagliptin "tablet 1 vs tablet 2" distinction is a drug-product-version comparison local to the Retlich 2015 popPK dataset; tablet formulation 1 was a development formulation that is not marketed. Mirrors the `FORM_DP2` (sarilumab) and `FORM_P2F2` (isatuximab) entries under the `FORM_*` family. Set to 0 for routine marketed-formulation simulation. Ratified canonically alongside the Retlich 2015 linagliptin extraction.

### FORM_DP2 (**canonical for sarilumab drug-product version 2 indicator**)
- **Description:** 1 = sarilumab drug product 2 formulation (used in some phase I studies and the dose-ranging phase II study), 0 = other drug product (DP1 or DP3; DP3 is the commercial formulation).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (DP1 or DP3).
- **Source aliases:**
  - `DP2` -- used in `Xu_2019_sarilumab.R`.
- **Example models:** `Xu_2019_sarilumab.R`.
- **Notes:** Affects both CLO/F (1.30x multiplier) and Ka (0.663x multiplier) in Xu 2019. Set to 0 for routine commercial-formulation simulation.

### FORM_VISMO_PHASEI (**canonical for vismodegib Phase I (dry-blend capsule) formulation indicator**)
- **Description:** 1 = subject received the Phase I clinical-development vismodegib formulation (dry-blend capsules used in early-phase studies SHH3925g and the Phase I cohort of SHH4610g); 0 = subject received the Phase II / commercial-bound formulation (wet-granulation capsules used in SHH4476g, SHH4433g, and most of SHH4683g / SHH4871g). Per-subject (regimen-fixed) categorical indicator.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (Phase II / commercial wet-granulation capsule; F = 1 fixed as the reference in Lu 2015 Eq. 4 and the typical-value `ka` reference in Eq. 4).
- **Source aliases:** none known; Lu 2015 reports formulation as a paper-defined "Phase I vs Phase II formulation" categorical without committing to a single column name.
- **Example models:** `Lu_2015_vismodegib.R` (multiplicative shifts on Ka and relative bioavailability F: typical Ka = 9.025 1/day for the Phase II reference in patients, with `exp(-0.602) = 0.55x` for the Phase I formulation and `exp(0.671) = 1.96x` for healthy volunteers; F = 1 for the Phase II reference, F = 0.346 for the Phase I formulation in patients and 0.836 in HV).
- **Notes:** Specific scope because the "Phase I (dry blend) vs Phase II (wet granulation) capsule" contrast is tied to the vismodegib drug-product-version comparison in Lu 2015. Mirrors the `FORM_DP2` (sarilumab), `FORM_P2F2` (isatuximab), and `FORM_LINAG_TAB1` (linagliptin) entries under the `FORM_*` family of drug-product-version indicators. Set to 0 for routine commercial-formulation simulation. Ratified canonically alongside the Lu 2015 vismodegib extraction.

### FORM_TAC_IR (**canonical for tacrolimus immediate-release vs prolonged-release formulation indicator**)
- **Description:** 1 = subject received the twice-daily immediate-release tacrolimus formulation (Prograf, Astellas) administered every 12 hours; 0 = subject received the once-daily prolonged-release tacrolimus formulation (Advagraf in Europe, Astagraf XL in the US; both Astellas) administered every 24 hours. Per-subject (regimen-fixed) categorical indicator used in popPK analyses that pool the two oral tacrolimus formulations and test formulation as a covariate on absorption and disposition parameters.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (Advagraf prolonged-release; the typical-value reference for Ktr and Vc/F in Woillard 2011 Table 4).
- **Source aliases:**
  - `study` -- used in `Woillard_2011_tacrolimus.R` (Woillard 2011 Methods: "study factor (assumed to be similar to drug formulation) ... study = 1 for the Prograf cohort, study = 0 for the Advagraf cohort"). The source paper's `study` and `formulation` factors are aliased; the canonical preserves the source paper's orientation (Prograf = 1 = IR, Advagraf = 0 = PR), which also matches the standard convention of treating the older immediate-release formulation as the non-reference level.
- **Example models:** `Woillard_2011_tacrolimus.R` (multiplicative effects per Woillard 2011 Table 4: `Ktr = theta1 * theta2^FORM_TAC_IR` with `theta1 = 3.34/h` and `theta2 = 1.53` (Prograf ~53 % faster absorption than Advagraf); `Vc/F = theta6 * theta7^FORM_TAC_IR` with `theta6 = 486 L` and `theta7 = 0.29` (Prograf Vc/F is 29 % of the Advagraf reference)).
- **Notes:** Specific scope because the Prograf-vs-Advagraf contrast is tied to oral tacrolimus formulation comparisons; mirrors the existing `FORM_DP2` (sarilumab), `FORM_P2F2` (isatuximab), `FORM_LINAG_TAB1` (linagliptin), and `FORM_VISMO_PHASEI` (vismodegib) entries under the `FORM_*` family of drug-product-version indicators. The Woillard 2011 paper notes the formulation effect partially confounds with time-post-transplant (Prograf cohort sampled within the first 6 months post-transplant, Advagraf cohort > 12 months post-transplant). Future tacrolimus models that include Envarsus XR (modified-release once-daily granules) or LCP-Tacro (life-cycle-pharma melt-extrusion tablets) should register a sibling canonical (e.g., `FORM_TAC_ENVARSUS`) rather than overloading `FORM_TAC_IR`. Ratified canonically alongside the Woillard 2011 tacrolimus extraction.

### dilution (**canonical for diluted-drug-product indicator**)
- **Description:** 1 = drug diluted (Soehoel 2022 study D2213C00001), 0 = not diluted.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Example models:** `Soehoel_2022_tralokinumab.R`.
- **Notes:** Lower-case preserved from source; future models should rename to `DILUTION`. Kept as alias here to match existing file.

### nonECZTRA (**canonical for non-ECZTRA-trial indicator**)
- **Description:** 1 = not the ECZTRA trial; 0 = ECZTRA.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Example models:** `Soehoel_2022_tralokinumab.R`.
- **Notes:** Mixed case preserved from source; future models should rename to `NON_ECZTRA` or `STUDY_NON_ECZTRA`.

### SEASON2 (**canonical for second RSV season at dosing indicator**)
- **Description:** 1 = second RSV season at dosing, 0 = first RSV season.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Example models:** `Clegg_2024_nirsevimab.R`.
- **Notes:** Study-specific but semantically general (second-exposure indicator). Promote to general if a second RSV-season model adopts the same semantics.

### COHDOSE (**canonical for randomized dose cohort (mg/kg)**)
- **Description:** Randomized dose cohort expressed in mg/kg. Subject-level (time-fixed) covariate carrying the per-subject cohort dose in a study where each subject remained on a single escalating-cohort dose for the full dosing period.
- **Units:** mg/kg
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with power scaling `(COHDOSE / ref)^exponent`. Reference value observed: 1 mg/kg in Narwal 2013.
- **Source aliases:**
  - `DOSE` -- used in `Narwal_2013_sifalimumab.R` (the paper's Eq. 3 variable name; renamed to `COHDOSE` here to avoid colliding with the rxode2/nlmixr2 event-column convention where `DOSE` or `AMT` carries the administered dose).
- **Example models:** `Narwal_2013_sifalimumab.R` (reference 1 mg/kg, exponent 0.0542 on CL).
- **Notes:** Scope: specific because the interpretation depends on a study design where each subject stays on a single dose cohort. For fixed-dose simulations, set `COHDOSE = nominal_dose_mg / WT` per subject. When the subject receives a weight-based dose, `COHDOSE` is the mg/kg label (0.3, 1, 3, or 10 mg/kg for the MI-CP152 cohorts).

### DOSE (**canonical for current administered dose level supplied as a data column**)
- **Description:** Continuous covariate carrying the administered dose level (in mg) as a per-record data column. Two complementary use cases share this canonical:
  (a) **Per-subject assigned dose** -- each subject's fixed assigned dose level (in mg) across the study, used as a power-style or stratified covariate when the population PK model detects a dose-dependent shift in a PK parameter (central volume, clearance, etc.).
  (b) **Time-varying current administered dose** -- the current daily dose at the time of the record, used in PD-only models that derive a per-cycle exposure metric (e.g., `AUC = DOSE / CLI`) from a posthoc-CL covariate without instantiating a PK ODE. The column is set to 0 during off-treatment periods (drug holidays, placebo arm) so the derived exposure becomes 0.
- **Units:** mg (document per-model in `covariateData[[DOSE]]$units` if a different dose unit is used).
- **Type:** continuous
- **Reference category:** n/a -- used with power scaling `(DOSE / ref)^exponent` for use case (a), or directly inside derived-exposure expressions for use case (b). Reference values observed: 600 mg in `Zheng_2016_sifalimumab.R` (middle of the 200/600/1200 mg phase IIb dose range).
- **Source aliases:**
  - `Dose` -- used in `Zheng_2016_sifalimumab.R` and `Castro-Surez_2020_nimotuzumab.R`.
  - `DOS` -- used in the Hansson 2013 sunitinib biomarker / TGI / fatigue PD-model family (DDMODEL00000197 and siblings) as a per-record sunitinib dose column.
- **Example models:** `Zheng_2016_sifalimumab.R` (power effect on V1 with exponent 0.06), `Castro-Surez_2020_nimotuzumab.R` (binary-indicator usage `(DOSE == 50)` applying a 53 % decrease in V1 for the 50 mg cohort), `Hansson_2013a_sunitinib.R` (DDMODEL00000197; time-varying record-level dose feeding `AUC = DOSE / CLI`), `Hansson_2013b_sunitinib.R` (DDMODEL00000198; same time-varying `AUC = DOSE / CLI` form for the tumor growth inhibition model), `Schindler_2016_sunitinib.R` (DDMODEL00000221; same `AUC = DOSE / CLI` form, with the daily-dose column toggling between 50 mg/day on-cycle and 0 on dose-holiday records), `Girard_2012_pimasertib.R` (linear coefficient on the dropout-hazard log-rate: `exp(beta * DOSE)` Weibull multiplier; per-subject daily dose, observed range 1-255 mg).
- **Notes:** Distinct from `DOSE_70MG` (binary indicator for a specific dose group in a trinary-dose design) and from the rxode2/nlmixr2 event column `amt` (which carries the administered dose at dose events). For use case (a), the values are typically time-fixed per subject; for use case (b), the values are time-varying with on/off cycling -- for sunitinib 4-weeks-on / 2-weeks-off cycling, set `DOSE = nominal_daily_mg` (e.g., 50) during on-cycles and 0 during off-cycles or for the placebo arm. Per-model `covariateData[[DOSE]]$notes` should state which use case applies.

### CD (**canonical for cumulative cladribine dose (time-varying)**)
- **Description:** Time-varying cumulative cladribine dose (mg total dose, not body-weight-normalized) administered to each subject up to the current observation time. Stays at zero during the placebo arm and during the pre-dose baseline phase, rises stepwise across the dosing schedule (cladribine is given as short oral pulses), and remains constant between dose events.
- **Units:** mg
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used inside `EXPS = CD * 104.5 / CRL` as an exposure surrogate driving an Emax-style symptomatic effect on disease progression, not as a power-form covariate.
- **Source aliases:** `CD` -- column name used in the DDMODEL00000223 input dataset (`Simulated_Novakovic_2016_multiplesclerosis_cladribine_irt.csv`).
- **Example models:** `Novakovic_2017_cladribine.R`.
- **Notes:** Distinct from `DOSE` (per-subject assigned dose level, time-fixed) and `COHDOSE` (mg/kg cohort label, time-fixed). `CD` is the *cumulative* dose accrued at each timepoint, supplied as a time-varying covariate column rather than via dosing events because the Novakovic 2017 model does not carry an explicit cladribine-PK compartment. Scope: specific because the constant 104.5 inside the exposure-surrogate equation is hard-coded for cladribine in the source.

### TRT (**canonical for treatment-cohort indicator**)
- **Description:** Treatment-cohort indicator used in the Novakovic 2017 cladribine IRT model. 0 = placebo, 1 = cladribine 3.5 mg/kg cumulative-dose cohort, 2 = cladribine 5.25 mg/kg cumulative-dose cohort.
- **Units:** (categorical)
- **Type:** categorical
- **Scope:** specific
- **Reference category:** 0 (placebo).
- **Source aliases:** `TRT` -- column name used in the DDMODEL00000223 input dataset (`Simulated_Novakovic_2016_multiplesclerosis_cladribine_irt.csv`).
- **Example models:** `Novakovic_2017_cladribine.R`.
- **Notes:** Gates the symptomatic and protective drug-effect terms via `TRT >= 1 && t > 0`; the categorical level (1 vs 2) is informational because the dose-response is driven by the time-varying `CD` covariate and the per-cohort dosing schedule, not by `TRT` itself. Scope: specific because the cohort labelling (3.5 vs 5.25 mg/kg cumulative dose over 2 years) is tied to the CLARITY-program cladribine dosing schedule. Future models that need a generic on-treatment indicator should register a new canonical name (e.g., `ON_TREATMENT`) rather than reusing `TRT`.

### DRUG_ORMU (**canonical for Ormutivimab vs HRIG drug-product comparator indicator**)
- **Description:** 1 = subject received Ormutivimab (a recombinant human anti-rabies IgG1 monoclonal antibody, also referred to in the source paper as rHRIG); 0 = subject received plasma-derived human rabies immunoglobulin (HRIG) comparator (or placebo + vaccine without passive antibody). Per-subject (time-fixed) categorical indicator carrying the head-to-head drug-arm assignment in a randomized rabies-vaccine pharmacodynamics study where rabies virus neutralizing antibody (RVNA) activity is modelled with a time-dependent Emax response to the vaccine and the passive antibody product modifies the typical-value Emax / ET50.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (HRIG / no passive antibody).
- **Source aliases:**
  - `antibody type` -- Zhang 2022 Results section 3.4 covariate-effect prose ("antibody type was determined as the covariate significantly affecting the model"); the source NMTRAN column name is not separately reported.
- **Example models:** `Zhang_2022_ormutivimab.R` (additive typical-value shifts on the linear-scale Emax and ET50 of the vaccine-induced RVNA Emax model: `Emax_tv = exp(lEmax) + e_drug_ormu_Emax * DRUG_ORMU` with `e_drug_ormu_Emax = +0.143 IU/mL` and `ET50_tv = exp(lET50) + e_drug_ormu_ET50 * DRUG_ORMU` with `e_drug_ormu_ET50 = -3.8 day`, yielding a higher and faster vaccine-induced antibody peak in the Ormutivimab arms relative to the HRIG comparator).
- **Notes:** Distinct from the `FORM_*` family (within-product formulation-version contrasts of a single drug) because the contrast here is between two biologically distinct products: HRIG is a polyclonal plasma-derived immunoglobulin, while Ormutivimab is a recombinant monoclonal antibody (CHO-cell-produced; the first rhRIG approved in China). Specific scope because the head-to-head HRIG-vs-rHRIG comparator design is tied to the Zhang 2022 phase II rabies-vaccine study. Future head-to-head rhRIG-vs-HRIG popPD models (e.g., SII Rabishield or Twinrab against HRIG) should register a sibling canonical (`DRUG_SIIRMAB`, `DRUG_TWINRAB`) rather than overloading `DRUG_ORMU`; cross-product comparisons that need both indicators in the same dataset can carry them as independent binaries with HRIG as the shared reference. Ratified canonically alongside the Zhang 2022 ormutivimab extraction.

### DOSE_70MG (**canonical for 70 mg dose regimen indicator**)
- **Description:** 1 = subject is on the 70 mg SC Q4W dose regimen, 0 = subject is on the 210 or 490 mg SC Q4W regimen.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (210 mg or 490 mg Q4W regimen).
- **Source aliases:** derived per subject from the trial-assigned dose level.
- **Example models:** `Kotani_2022_astegolimab.R`.
- **Notes:** Zenyatta-study categorical covariate flagging the 70 mg group (lowest dose), modeled as a -15.3% relative change on relative bioavailability. Modeled by Kotani 2022 as `70 mg vs {210 mg, 490 mg}` combined reference.

### DOSE_50MG (**canonical for 50 mg dose administration indicator**)
- **Description:** 1 = dose record is a 50 mg SC administration, 0 = all other SC doses (100, 150, 200, 300 mg) and all IV doses.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (100-300 mg SC or any IV dose).
- **Source aliases:** derived per dose record from the administered amount (`AMT`).
- **Example models:** `Othman_2014_daclizumab.R`, `Diao_2016_daclizumab_cd25.R`, `Diao_2016_daclizumab_cd56bright.R`, `Diao_2016_daclizumab_treg.R`.
- **Notes:** Othman 2014 estimated two separate absolute bioavailabilities because of non-linear dose-normalized exposure at the 50 mg SC dose -- F = 0.84 for the therapeutic 100-300 mg SC range and F = 0.57 for the 50 mg SC cohort. Encoded here as a record-level indicator so the covariate effect `e_dose_50mg_f = 0.57/0.84 - 1 = -0.321` scales bioavailability only on 50 mg SC doses. For clinical-range simulation (150 mg SC Q4W Phase III regimen) leave `DOSE_50MG = 0`. The Diao 2016 PK/PD models inherit the Othman 2014 PK backbone verbatim; they carry `DOSE_50MG` even though the Diao 2016 RRMS regimens are 150 / 300 mg SC only.

### STUDY1 (**canonical for Study-1 cohort indicator**)
- **Description:** 1 = subject enrolled in Study 1 of the Cirincione 2017 pooled analysis, 0 = other. Used to switch the residual-error magnitude per study.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (all other studies; combined with `STUDY5 = 0` selects the pooled "other" residual error).
- **Source aliases:**
  - `DVID = "study1"` (character-valued study identifier; `STUDY1 = as.integer(DVID == "study1")`) -- legacy form previously used in `Cirincione_2017_exenatide.R`.
- **Example models:** `Cirincione_2017_exenatide.R`.
- **Notes:** Paired with `STUDY5`. When both are 0, the subject is in the pooled "other studies" residual-error group.

### STUDY5 (**canonical for Study-5 cohort indicator**)
- **Description:** 1 = subject enrolled in Study 5 of the Cirincione 2017 pooled analysis, 0 = other. Used to switch the residual-error magnitude per study.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (all other studies).
- **Source aliases:**
  - `DVID = "study5"` (character-valued study identifier; `STUDY5 = as.integer(DVID == "study5")`) -- legacy form previously used in `Cirincione_2017_exenatide.R`.
- **Example models:** `Cirincione_2017_exenatide.R`.
- **Notes:** Paired with `STUDY1`. When both are 0, the subject is in the pooled "other studies" residual-error group.

### STUDY_PKU015 (**canonical for sapropterin PKU-015 pediatric study cohort indicator**)
- **Description:** 1 = subject enrolled in study PKU-015 (pediatric population pharmacokinetic study of sapropterin in infants and young children, 0-6 years old, of the Qi 2014 pooled analysis); 0 = study PKU-004 (adolescent / adult open-label extension study, >= 9 years old). Used to switch the residual-error magnitude per study under the log-transform-both-sides (LTBS) constant-CV residual model.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (PKU-004 adolescent / adult cohort).
- **Source aliases:** derived per subject from the trial identifier (`PKU-015` -> 1, `PKU-004` -> 0).
- **Example models:** `Qi_2014_sapropterin.R`.
- **Notes:** Qi 2014 Table 3 reports separate residual-error estimates for the two studies -- PKU-004 = 21.1% CV, PKU-015 = 30.2% CV under the LTBS approach. The `STUDY_PKU015` indicator selects between them. Specific scope because the indicator is tied to the BioMarin sapropterin clinical-development program (PKU-004 = phase 3b extension, PKU-015 = phase 3b pediatric).

### PHASE2 (**canonical for Phase II study cohort indicator**)
- **Description:** 1 = subject enrolled in the Phase II study (MORAb-003-002) of the Farrell 2012 pooled analysis; 0 = Phase I study (MORAb-003-001). Used to switch the residual-error magnitude per study.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (Phase I).
- **Source aliases:** derived per subject from the trial identifier (`MORAb-003-001` -> 0, `MORAb-003-002` -> 1).
- **Example models:** `Farrell_2012_farletuzumab.R`.
- **Notes:** Farrell 2012 Table 3 reports separate residual-error estimates for the two studies -- Phase I uses a proportional-only model (sigma = 20.5%); Phase II uses a combined additive + proportional model (sigma_prop = 34.9%, sigma_add = 7.94 ug/mL). The `PHASE2` indicator selects between them.

### PHASE1 (**canonical for Phase I study cohort indicator**)
- **Description:** 1 = subject enrolled in a Phase 1 study of the Valenzuela 2025 pooled analysis (MOM-M281-001, MOM-M281-007, MOM-M281-010, EDI1001, EDI1002 -- healthy participants); 0 = Phase 2 study (MOM-M281-004 / Vivacity-MG -- participants with gMG). Used to switch the proportional PK residual-error magnitude per study phase.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (Phase 2).
- **Source aliases:** derived per subject from the trial identifier (Phase 1 protocols -> 1, `NCT03772587` Vivacity-MG -> 0).
- **Example models:** `Valenzuela_2025_nipocalimab.R`.
- **Notes:** Valenzuela 2025 Table 3 reports proportional PK residual 0.0834 (Phase 1) vs 0.367 (Phase 2). Distinct from Farrell 2012 `PHASE2` -- the reference category is inverted (Valenzuela 2025 picks Phase 1 as the 1-level).

### STUDY_C2201 (**canonical for Bienczak 2025 ligelizumab study C2201 cohort indicator**)
- **Description:** 1 = subject enrolled in study C2201 (NCT02477332; Novartis Phase 2b ligelizumab dose-finding study in adult CSU patients) of the Bienczak 2025 pooled ligelizumab PopPK analysis; 0 = any other study in the pool (A2103, C2101, C2202, C2302, or C2303). Used to switch the typical CL/F magnitude in study C2201.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (any non-C2201 cohort in the Bienczak 2025 pool: pooled adult / adolescent CSU patients from C2202 / C2302 / C2303 and adult healthy volunteers from A2103 / C2101).
- **Source aliases:** derived per subject from the trial identifier (`C2201` -> 1, else -> 0).
- **Example models:** `Bienczak_2025_ligelizumab.R` (Table S6: study C2201 on CL/F = 0.176, log-additive; `cl *= exp(0.176)` for C2201 subjects).
- **Notes:** Specific scope because the contrast is tied to the Novartis ligelizumab CSU development program. Subject-level / time-fixed; set once from the trial identifier on each subject record. The C2201 effect was retained in the final model because the residual unexplained CL/F differed between C2201 and the other studies after accounting for body weight, IgE, ADA, and disease-state covariates.

### STUDY_LBSL (**canonical for early-phase belimumab LBSL01 / LBSL02 study indicator**)
- **Description:** 1 = subject enrolled in study LBSL01 (NCT00657007) or LBSL02 (NCT00071487) -- the two early-phase belimumab studies that used a different ELISA-based bioanalytical assay; 0 = any other belimumab study in the Zhou 2021 pooled analysis. Used to switch CL and V1 magnitudes per study group (effectively an assay / early-development PK adjustment).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (later-phase studies using the electrochemiluminescence assay).
- **Source aliases:**
  - `INDR` -- used in `Zhou_2021_belimumab.R` (Zhou 2021 Table 2 footnote: study indicator).
- **Example models:** `Zhou_2021_belimumab.R` (multiplicative factors 1.63 on CL and 1.26 on V1 when STUDY_LBSL = 1).
- **Notes:** Conceptually similar to `STUDY1` / `PHASE2` / `ELISA` / `PHASE1` (per-study switches) but specific to the belimumab program. Subject-level (time-fixed); set from the trial identifier on each subject record.

### ELISA (**canonical for ELISA-vs-ECLIA bioanalytical assay indicator**)
- **Description:** 1 = serum nipocalimab concentration measured by ELISA assay (LLOQ 0.150 ug/mL; studies MOM-M281-001, MOM-M281-007, MOM-M281-010); 0 = measured by ECLIA assay (LLOQ 0.010 ug/mL; studies EDI1001, EDI1002, MOM-M281-004). Used to switch the additive PK residual-error magnitude per assay.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (ECLIA).
- **Source aliases:** derived per study from the bioanalytical method (Text S2 of Valenzuela 2025).
- **Example models:** `Valenzuela_2025_nipocalimab.R`.
- **Notes:** Valenzuela 2025 Table 3 reports additive PK residual 0.445 nmol/L (ELISA) vs 0.0342 nmol/L (ECLIA). Assay choice is study-fixed; include as a per-observation indicator so the correct assay residual is applied even when observations from different studies are pooled.

### STUDY_M281_004 (**canonical for MOM-M281-004 (Vivacity-MG) study indicator**)
- **Description:** 1 = subject enrolled in study MOM-M281-004 (Vivacity-MG; NCT03772587; phase 2 in generalized myasthenia gravis); 0 = any other study in the Valenzuela 2025 pooled analysis. Used to switch the IgG-baseline scaling factor.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-Vivacity-MG studies).
- **Source aliases:** derived per subject from the trial identifier (`NCT03772587` -> 1, else -> 0).
- **Example models:** `Valenzuela_2025_nipocalimab.R`.
- **Notes:** Valenzuela 2025 estimated a slightly lower IgG baseline in MOM-M281-004 participants (baseline scaled by `FRIgG0_M281_004 = 0.777` vs. 1 in other studies). Distinct from the disease-state indicator implied by `gMG` -- it is specifically the Vivacity-MG study flag because the IgG baseline factor was only estimated for that study.

### STUDY_ABA2_HLA78 (**canonical for ABA2 trial HLA 7/8-matched cohort indicator**)
- **Description:** 1 = subject enrolled in the ABA2 hematopoietic-cell-transplant trial (IM101311; NCT01743131) HLA 7/8 (one-allele-mismatched donor) cohort, 0 = any other study in the Takahashi 2023 pooled abatacept population PK analysis (RA/JIA reference and ABA2 8/8 cohort).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-ABA2-7/8 -- pooled adult RA / pediatric JIA cohort, plus ABA2 HLA 8/8 cohort which is flagged separately by `STUDY_ABA2_HLA88`).
- **Source aliases:** derived per subject from the trial-cohort identifier (`Cohort = ABA2 7/8` -> 1, else -> 0). Takahashi 2023 Supplemental Table 4 names the corresponding theta `thetaCohort_CL` / `thetaCohort_Vc`.
- **Example models:** `Takahashi_2023_abatacept.R` (multiplicative `Ratio` factors on CL = 0.70 and on Vc = 0.99 vs the RA/JIA reference; values from Takahashi 2023 Supplemental Table 4).
- **Notes:** Pairs with `STUDY_ABA2_HLA88` to reproduce the three-level cohort categorical (RA/JIA, ABA2 7/8, ABA2 8/8) the paper reports as the only retained categorical PK covariate. At most one of `STUDY_ABA2_HLA78` and `STUDY_ABA2_HLA88` is 1 per subject; both 0 reproduces the RA/JIA reference. Scope: specific because the contrast is tied to the ABA2-vs-RA/JIA pooling design.

### STUDY_ABA2_HLA88 (**canonical for ABA2 trial HLA 8/8-matched cohort indicator**)
- **Description:** 1 = subject enrolled in the ABA2 hematopoietic-cell-transplant trial (IM101311; NCT01743131) HLA 8/8 (allele-matched donor) cohort, 0 = any other study in the Takahashi 2023 pooled abatacept population PK analysis (RA/JIA reference and ABA2 7/8 cohort).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-ABA2-8/8 -- pooled adult RA / pediatric JIA cohort, plus ABA2 HLA 7/8 cohort which is flagged separately by `STUDY_ABA2_HLA78`).
- **Source aliases:** derived per subject from the trial-cohort identifier (`Cohort = ABA2 8/8` -> 1, else -> 0). Takahashi 2023 Supplemental Table 4 names the corresponding theta `thetaCohort_CL` / `thetaCohort_Vc`.
- **Example models:** `Takahashi_2023_abatacept.R` (multiplicative `Ratio` factors on CL = 0.91 and on Vc = 1.32 vs the RA/JIA reference; values from Takahashi 2023 Supplemental Table 4).
- **Notes:** Pairs with `STUDY_ABA2_HLA78`. At most one of the two indicators is 1 per subject; both 0 reproduces the RA/JIA reference. Scope: specific.

### DLVL (**canonical for source-protocol dose-level integer indicator**)
- **Description:** Integer dose-level / protocol-arm indicator carried per subject in the DDMODEL00000281 lidocaine bundle's simulated dataset (`Simulated_Lid_B04_ddmore.csv`). Values 1-4 (or higher) flag distinct study-protocol dose / regimen tiers; the `NA_NA_lidocaine.R` model binarises as `DLVL > 2` to switch the typical-value baselines for both the GX elimination rate constant K30 and the lidocaine apparent central volume V1 between a "low" (DLVL <= 2) and a "high" (DLVL > 2) regimen.
- **Units:** (integer-coded categorical)
- **Type:** categorical
- **Scope:** specific
- **Reference category:** the binary form `DLVL <= 2` is the reference (THETA(4) for K30 and THETA(14) for V1 in the source `.ctl`); `DLVL > 2` selects the higher-exposure regimen (THETA(5) and THETA(15) respectively).
- **Source aliases:** `DLVL` -- the column header used in the DDMORE bundle's `.ctl` `$INPUT` and the Simulated_Lid_B04_ddmore.csv data file.
- **Example models:** `NA_NA_lidocaine.R` (DDMODEL00000281; binary derivation `DLVL_HIGH = as.integer(DLVL > 2)` on K30 base and V1 base).
- **Notes:** Specific scope because the integer-coded dose / regimen tiers are paper-specific to the lidocaine BAST.dat ("4-cRUN249") study and the linked publication is not on disk for this extraction. The binary threshold `> 2` reproduces the source `.ctl` line `IF(DLVL.GT.2)P1=0`. If a future model needs a different dose-level binarisation or a continuous treatment, register a distinct canonical name rather than overloading `DLVL`.

### S1A2 (**canonical for source-protocol CYP1A2 substrate / co-medication categorical indicator**)
- **Description:** Categorical CYP1A2-modifying co-medication / phenotype indicator carried per subject in the DDMODEL00000281 lidocaine bundle's simulated dataset. Integer code; in the `NA_NA_lidocaine.R` model the value `S1A2 == 3` selects the "CYP1A2 inducer present" sub-cohort (lidocaine N-deethylation to MEGX is CYP1A2-mediated, so the modifier acts on the GX elimination rate constant K30 in the source's parameterisation). Other integer codes (0, 1, 2) are pooled into the reference.
- **Units:** (integer-coded categorical)
- **Type:** categorical
- **Scope:** specific
- **Reference category:** `S1A2 != 3` (i.e., values 0, 1, 2) -- pooled into the reference.
- **Source aliases:** `S1A2` -- the column header used in the DDMORE bundle's `.ctl` `$INPUT` and the simulated dataset, with sibling columns `D1A2` and `H1A2` carried in the data file but dropped via `=DROP` in the source `.ctl`.
- **Example models:** `NA_NA_lidocaine.R` (DDMODEL00000281; binary derivation `S1A2_IND = as.integer(S1A2 == 3)` on the GX elimination rate constant K30).
- **Notes:** Specific scope because the integer codes for `S1A2` are paper-specific to the lidocaine BAST.dat study and the linked publication is not on disk for this extraction. The exact biological meaning of each integer level (0/1/2/3) is not fully reconstructable from the bundle alone -- the natural interpretation, given the column name encodes "CYP1A2" and the model attaches a sizeable positive K30 modifier of +0.853 to the level-3 cohort, is a CYP1A2-induction or smoking / inducer co-medication indicator. Sibling columns `D1A2` (donor / inhibitor?) and `H1A2` (host / inhibitor?) are dropped in the source `.ctl` so only the level-3 indicator is structurally identifiable from the surviving model code. If a future model needs a richer encoding of CYP1A2 modulation, register a separate canonical (e.g., `CYP1A2_IND`) rather than overloading `S1A2`.

## Occasion / period (IOV)

### OCC (**canonical for the integer-valued occasion / period column**)
- **Description:** Integer-valued occasion / period indicator for inter-occasion-variability (IOV) modelling. Values `1`, `2`, ..., `N` identify the occasion to which each observation belongs (typically a dosing visit, study period, or sampling occasion). Time-varying within subject; constant within an occasion.
- **Units:** (count)
- **Type:** categorical
- **Scope:** general
- **Reference category:** n/a -- `OCC` is decomposed inside `model()` into mutually-exclusive binary indicators, e.g., `oc1 <- (OCC == 1)`, `oc2 <- (OCC == 2)`, ..., that are then multiplied against per-occasion `eta*` slots.
- **Source aliases:**
  - `OCC` -- used in `Jonsson_2011_ethambutol.R` (DDMODEL00000220 NMTRAN `$INPUT` column; values 1..4).
- **Example models:** `Jonsson_2011_ethambutol.R` (4-occasion IOV on log-CL; `cl <- exp(lcl + etalcl + oc1 * etalcl_oc1 + oc2 * etalcl_oc2 + oc3 * etalcl_oc3 + oc4 * etalcl_oc4) * (WT/50)^0.75`, where each `etalcl_oc<k>` is a separate `~ fix(0.127)` after the first to encode NONMEM `$OMEGA BLOCK(1) SAME`).
- **Notes:** `OCC` is the recommended canonical for new IOV-using models -- the binary `ooc1..oocN` indicators below remain canonical for legacy / pre-existing models that ship the data already-decomposed. Ratified canonically on 2026-05-06.

### ooc1, ooc2, ooc3, ooc4 (**canonical for mutually-exclusive crossover occasion indicators**)
- **Description:** Mutually exclusive occasion indicators for a crossover / multi-period design. Exactly one is 1 per observation.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Example models:** `Xie_2019_agomelatine.R`.
- **Notes:** Lower case preserved from source file. Pre-existing legacy form; new models should prefer the integer-valued `OCC` canonical above and decompose into binary indicators inside `model()`.

### DAY14 (**canonical for day-14-post-treatment-initiation landmark indicator**)
- **Description:** Binary within-subject landmark indicator: 1 = the observation falls on or after day 14 of treatment (post-nutritional-rehabilitation steady state in the Archary 2019 / MATCH trial of severely malnourished HIV-infected children), 0 = the observation falls before day 14 (acute / pre-rehabilitation baseline; e.g., day 1 of antiretroviral treatment). Time-varying within subject -- gates a step change in typical-value PK parameters that the source paper attributes to nutritional recovery + auto-induction of hepatic metabolism over the first ~2 weeks of treatment.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (day 1 / pre-rehabilitation).
- **Source aliases:** none known; Archary 2019 reports the day-1-vs-day-14 contrast directly in the parameter table (separate typical-value rows for "day 1" and "day 14").
- **Example models:**
  - `Archary_2019_abacavir.R` (multiplicative additive shift on CL/F: `cl <- exp(lcl + etalcl) * (WT/7)^0.75 * (1 + e_day14_cl * DAY14)` with `e_day14_cl = 0.760`, encoding the day-1 typical CL/F = 3.33 -> day-14 CL/F = 5.86 L/h per 7 kg step).
  - `Archary_2019_lamivudine.R` (multiplicative additive shift on ka: `ka <- exp(lka) * (1 + e_day14_ka * DAY14) * exp(etalka)` with `e_day14_ka = 0.133`, encoding the day-1 typical ka = 0.30 -> day-14 ka = 0.34 /h step).
- **Notes:** Specific scope because the day-14 cutoff is tied to the MATCH-trial nutritional-rehabilitation timeline (severely malnourished children, two-week re-feeding window per WHO guidelines); future studies that test a similar within-subject step change with a different cutoff (day 7, day 28, etc.) should register a new canonical (e.g., `DAY28`). Distinct from `OCC` and `ooc<n>` (which decompose multi-occasion sampling for IOV), from `TRT_PHASE` (which gates active-vs-baseline study-phase contributions), and from `EARLY_ART` (which is a between-subject randomization-arm indicator, not a within-subject landmark). Data assemblers can derive `DAY14 = as.integer(time_post_treatment_start_days >= 14)` for a regularly-sampled multi-day study. Ratified canonically on 2026-05-08.

### CYCLE (**canonical for dose-number / treatment-cycle counter**)
- **Description:** Dose-number / treatment-cycle counter (1 = first dose or cycle, 2 = second, ...). Integer count, time-varying across a multi-dose / multi-cycle treatment course, incremented at each new administration.
- **Units:** (count)
- **Type:** count
- **Scope:** general
- **Reference category:** n/a -- used either with a power-covariate form `CYCLE^Fm` (Fm typically negative) to capture cycle-over-cycle decline in a derived quantity such as ADC-to-payload conversion fraction (Li 2017 brentuximab vedotin), or with a piecewise indicator `CYCLE == 1 vs CYCLE >= 2` to capture a step change in PK between the first and subsequent administrations (Hong 2025 datopotamab deruxtecan; Huynh 2026 VRC07-523LS).
- **Source aliases:** `CYCLE` -- used in `Li_2017_brentuximab.R`, `Hong_2025_datopotamab.R`, `Lu_2022_patritumab.R`, and `Huynh_2026_VRC07523LS.R` with the same canonical name.
- **Example models:**
  - `Li_2017_brentuximab.R` (exponent on the fraction of ADC that converts to MMAE by proteolytic degradation, Fm = -0.261, to reflect tumor-burden reduction across successive treatment cycles).
  - `Hong_2025_datopotamab.R` (cycle-1 vs cycle-2+ piecewise scaling Factor1 = 0.696 on the DAR equation that drives DXd formation rate from total Dato-DXd elimination).
  - `Lu_2022_patritumab.R` (cycle-1 vs cycle-2+ piecewise scaling Factor1 = theta = 0.648 on the payload-to-intact-drug ratio PIR that scales DXd release rate from intact ADC).
- **Notes:** Must be >= 1 throughout (`CYCLE^Fm` is undefined at 0; the piecewise form requires `CYCLE` to be a positive integer at every observation row). Distinct from `ooc<n>` binary-occasion indicators: `CYCLE` is an integer count, not a mutually-exclusive set of indicator columns. Data-assembly helper: set `CYCLE = floor((TIME - TIME_FIRST_DOSE) / cycle_length_days) + 1` for a fixed-interval dosing regimen.

### T_ENTRY (**canonical for per-subject study-entry time on the model integration axis**)
- **Description:** Per-subject study-entry time, expressed on the same integration-time axis the model is solved on (years for AD disease-progression models). The dataset's `TIME` column carries the integration-time variable (with `TIME = 0` the integration origin, typically well before any subject's disease onset); `T_ENTRY` records, per subject, the integration-time value at which the subject's first observation falls. Used in models whose pre-study and post-study time semantics differ -- e.g., a disease-progression activation function defined on global disease time plus a placebo-effect term defined on time-since-study-entry need access to both reference clocks within the same `model()` block.
- **Units:** year (or whatever time unit the model's `units$time` field declares)
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- `T_ENTRY` is a per-subject time-offset covariate that anchors the time-since-study-entry clock used by post-entry-only model terms (placebo / learning effects, study-design transient drops).
- **Source aliases:** none standardized; this canonical originates with the Delor 2013 AD disease-progression extraction. The source paper's NONMEM dataset handles the global / study-time split internally (the dataset is staged with TIME = study time and the model uses a fixed offset to align with DOT); when porting to rxode2 / nlmixr2 the per-subject offset is exposed as a covariate so simulation event tables can carry it explicitly.
- **Example models:** `Delor_2013_alzheimer.R` (placebo-term clock: `t_pl_raw <- time - T_ENTRY; post_entry <- t_pl_raw > 0; placebo_term <- pl_indiv * (1 - exp(-kpl_indiv * t_pl_raw * post_entry)) * post_entry`).
- **Notes:** Scope: specific because the variable is paper-domain-bound -- it only makes sense for models that operate on a global disease-time axis distinct from per-subject study-entry timing. For a typical-value reproduction the user supplies `T_ENTRY` per subject in the simulation event table; a reasonable construction is `T_ENTRY = DOT_individual + a few years of established disease` so the patient is observed during disease progression (see the Delor 2013 vignette for the construction used to reproduce Figures 2-4). Ratified canonically on 2026-05-16 alongside the Delor 2013 extraction.
