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

### WT
- **Description:** Body weight (baseline or time-varying).
- **Units:** kg
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with allometric scaling `(WT / ref_wt)^exponent`. Reference weights observed: 70 kg (adults), 75 kg, 84.8 kg, 5 kg (infants).
- **Source aliases:** none known.
- **Example models:** `Clegg_2024_nirsevimab.R`, `Hu_2026_clesrovimab.R`, `Zhu_2017_lebrikizumab.R`, `Kovalenko_2020_dupilumab.R`, `CarlssonPetri_2021_liraglutide.R`, `Cirincione_2017_exenatide.R`, `Grimm_2023_gantenerumab.R`, `Grimm_2023_trontinemab.R`, `Kyhl_2016_nalmefene.R`, `Soehoel_2022_tralokinumab.R`, `Xie_2019_agomelatine.R`, `PK_2cmt_mAb_Davda_2014.R`, `phenylalanine_charbonneau_2021.R`, `Chua_2025_mirikizumab.R`, `Jackson_2022_ixekizumab.R`, `Kotani_2022_astegolimab.R`, `Ma_2020_sarilumab_anc.R`, `Ma_2020_sarilumab_das28crp.R`, `Moein_2022_etrolizumab.R`, `Tiraboschi_2025_amlitelimab.R`, `Robbie_2012_palivizumab.R`, `Bajaj_2017_nivolumab.R`, `Quartino_2019_trastuzumab.R`, `Wang_2020_ontamalimab.R`, `Fau_2020_isatuximab.R`.
- **Example models:** `Clegg_2024_nirsevimab.R`, `Hu_2026_clesrovimab.R`, `Zhu_2017_lebrikizumab.R`, `Kovalenko_2020_dupilumab.R`, `CarlssonPetri_2021_liraglutide.R`, `Cirincione_2017_exenatide.R`, `Grimm_2023_gantenerumab.R`, `Grimm_2023_trontinemab.R`, `Kyhl_2016_nalmefene.R`, `Soehoel_2022_tralokinumab.R`, `Xie_2019_agomelatine.R`, `PK_2cmt_mAb_Davda_2014.R`, `phenylalanine_charbonneau_2021.R`, `Chua_2025_mirikizumab.R`, `Jackson_2022_ixekizumab.R`, `Kotani_2022_astegolimab.R`, `Ma_2020_sarilumab_anc.R`, `Ma_2020_sarilumab_das28crp.R`, `Moein_2022_etrolizumab.R`, `Tiraboschi_2025_amlitelimab.R`, `Robbie_2012_palivizumab.R`, `Bajaj_2017_nivolumab.R`, `Quartino_2019_trastuzumab.R`, `Wang_2020_ontamalimab.R`, `Okada_2025_rocatinlimab.R`.
- **Notes:** Universal. Verify time-varying vs. baseline-only against the source paper.

### AGE
- **Description:** Subject age in years.
- **Units:** years
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a
- **Source aliases:** none.
- **Example models:** `Kyhl_2016_nalmefene.R`, `Zhu_2017_lebrikizumab.R`.
- **Notes:** Zhu 2017 normalizes as `AGE/40`.

### LBM
- **Description:** Lean body mass.
- **Units:** kg
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a
- **Source aliases:** none.
- **Example models:** `Kyhl_2016_nalmefene.R` (reference 56.28 kg, exponent 0.626 on CL).

### FFM (**canonical for fat-free mass**)
- **Description:** Fat-free mass derived from body weight, height, and sex via the Janmahasatian et al. formula (Clin Pharmacokinet 2005;44:1051-1065).
- **Units:** kg
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(FFM / ref)^exponent`. Reference values observed: 40.69 kg (Zhou 2021 belimumab pooled adult+pediatric SLE), 45 kg (Aguiar 2021, Crohn's disease cohort median).
- **Source aliases:** none; `FFM` is the universal abbreviation.
- **Example models:** `Zhou_2021_belimumab.R` (reference 40.69 kg; exponents 0.673 on CL and 0.891 on V1), `Aguiar_2021_ustekinumab.R` (reference 45 kg; power exponents 0.598 on CL, 0.590 on Vc, 0.586 on Vp).
- **Notes:** Distinct from `LBM` (lean body mass) which is sometimes computed by the Boer or Hume formulae. When the source paper reports the body-composition formula it used (e.g., Janmahasatian for FFM), record it in `covariateData[[FFM]]$notes`. FFM is preferred over total body weight when scaling monoclonal-antibody PK because mAb distribution is largely confined to extracellular fluid; muscle / lean tissue tracks extracellular volume better than total weight in heavier patients.

### BSA
- **Description:** Body surface area (typically computed by DuBois, Mosteller, or Haycock from height and weight).
- **Units:** m^2
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(BSA / ref)^exponent`.
- **Source aliases:** none.
- **Example models:** `Yamada_2025_zolbetuximab.R` (reference 1.70 m^2; exponents 1.06 on clearances and 0.968 on volumes).
- **Notes:** Oncology mAbs dosed by BSA (mg/m^2) often use BSA in place of body weight for allometric-style scaling. Document the BSA computation formula (DuBois / Mosteller / Haycock) the source paper used; if unstated, record "unspecified."

### BMI
- **Description:** Body mass index at baseline.
- **Units:** kg/m^2
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with a linear-deviation form (`1 + e * (BMI - ref)`) or a power form (`(BMI / ref)^e`). Document the reference value in `covariateData[[BMI]]$notes`.
- **Source aliases:** none known.
- **Example models:** `Chua_2025_mirikizumab.R` (reference 24.75 kg/m^2; linear-deviation effect on logit of bioavailability), `NA_NA_lidocaine.R` (DDMODEL00000281; binary stratification at threshold 27.93 kg/m^2 adding +0.939 to the GX rate constant K30 in the BMI > 27.93 cohort).
- **Notes:** Universal clinical-trial demographic. Derived as `WT / (height_m)^2`; assume time-fixed at baseline unless the source paper states otherwise.

### SEXF (**canonical for sex**)
- **Description:** Biological sex indicator, 1 = female, 0 = male.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (male).
- **Source aliases:**
  - `SEXM` (values inverted: `SEXF = 1 - SEXM`; effect coefficient sign and reference category both invert) -- used in `CarlssonPetri_2021_liraglutide.R`.
  - `SEX` with `"M"`/`"F"` strings -- derive `SEXF = as.integer(SEX == "F")`.
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical), `CarlssonPetri_2021_liraglutide.R` (alias `SEXM`), `Bajaj_2017_nivolumab.R` (male-indicator source; effect applied as `exp(coef * (1 - SEXF))` to preserve the paper's female-reference CL_REF / VC_REF), `Fau_2020_isatuximab.R` (exponential effect on Vc; reference category 0 = male), `Netterberg_2017_docetaxel.R` (multiplicative effect on baseline ANC: `BACOV *= (1 + theta * SEXF)`; source column `SEX` with 1 = male, 2 = female encoding, decomposed via `SEXF = as.integer(SEX == 2)`).
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

### CHILD
- **Description:** 1 = subject is a child, 0 = not a child.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (not child, i.e., adult baseline).
- **Source aliases:**
  - `PED` -- used in the Schoemaker 2018 LEV / BRV pediatric extrapolation (DDMODEL00000239) as the pediatric-vs-adult indicator that gates the Markov-amplitude term, the overdispersion IIV, and the four pediatric offsets on log baseline rate / mixture / placebo / Emax / EC50.
- **Example models:** `CarlssonPetri_2021_liraglutide.R`, `Schoemaker_2018_levetiracetam.R` (DDMODEL00000239).
- **Notes:** Age-group indicator used alongside `ADOLESCENT`; paper's age cutoffs must be captured in `covariateData[[CHILD]]$notes`.

### ADOLESCENT
- **Description:** 1 = subject is an adolescent, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0.
- **Example models:** `CarlssonPetri_2021_liraglutide.R`.
- **Notes:** Paired with `CHILD`. Document age cutoffs.

## Pediatric / maturation

### PAGE
- **Description:** Postmenstrual age in months (`GA_weeks / 4.35 + postnatal_months`). Time-varying.
- **Units:** months
- **Type:** continuous
- **Scope:** general
- **Example models:** `Clegg_2024_nirsevimab.R`, `Robbie_2012_palivizumab.R`.

### PNA
- **Description:** Postnatal age (chronological since birth). Time-varying.
- **Units:** months
- **Type:** continuous
- **Scope:** general
- **Example models:** `Hu_2026_clesrovimab.R`.

### GA
- **Description:** Gestational age at birth. Time-fixed per subject.
- **Units:** weeks
- **Type:** continuous
- **Scope:** general
- **Example models:** `Hu_2026_clesrovimab.R`, used implicitly in `Clegg_2024_nirsevimab.R` (folded into PAGE).

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

### BC_USE (**canonical for oral-contraceptive use indicator**)
- **Description:** Binary indicator of oral hormonal contraceptive use; `1` = currently taking an oral contraceptive (estrogen-progestin or progestin-only pill), `0` = not on hormonal contraception. Time-varying as women cycle on/off contraception across study occasions.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no oral contraceptive). Effect form in Allegaert 2015 is multiplicative: `CL_glucuronide *= theta_BC_USE` when `BC_USE == 1`, with `theta_BC_USE = 1.46` (estrogen-driven UGT2B7 induction).
- **Source aliases:**
  - `BC` (Allegaert 2015 NONMEM column; same orientation, no transformation) -- used in `Allegaert_2015_paracetamol.R`.
- **Example models:** `Allegaert_2015_paracetamol.R`.
- **Notes:** Specific scope because the canonical encoding pools all oral contraceptive types (combined / progestin-only) into a single binary; future models that need to distinguish formulations should register a finer-grained canonical (e.g., `BC_COMBINED`, `BC_PROGESTIN`).

## Renal / hepatic function

### UF (**canonical for instantaneous urine flow rate**)
- **Description:** Instantaneous urine flow rate (mL/h) measured over the urine collection interval that includes the current observation. Time-varying.
- **Units:** mL/h
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with a centered-linear-effect form `CL_renal = base + theta_UF * (UF - UF_ref)` with `UF_ref = 100 mL/h` in Allegaert 2015. A value of `0` is a sentinel for "no urine collected during the interval" (i.e., the urine pathway contribution is dropped); the linear-effect term is gated by `UF > 0` and not extrapolated below the centering reference.
- **Source aliases:**
  - `UF` (Allegaert 2015 NONMEM column; same orientation, no transformation) -- used in `Allegaert_2015_paracetamol.R`.
- **Example models:** `Allegaert_2015_paracetamol.R`.
- **Notes:** Specific scope because the centered-linear effect form with the `UF == 0` sentinel-zero rule reflects an Allegaert-specific convention rather than a universally-agreed-upon parameterization. A second model that uses a different effect form (e.g., direct `UF / UF_ref` proportional scaling, no zero-sentinel) should register its own canonical (e.g., `UF_PROP`) rather than reusing `UF` with conflicting semantics.

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
- **Example models:** `Cirincione_2017_exenatide.R` (MDRD eGFR), `Xu_2019_sarilumab.R` (measured CrCl BSA-normalized), `Kotani_2022_astegolimab.R` (MDRD eGFR), `Li_2019_abatacept.R` (cGFR).
- **Example models:** `Cirincione_2017_exenatide.R` (MDRD eGFR), `Xu_2019_sarilumab.R` (measured CrCl BSA-normalized), `Kotani_2022_astegolimab.R` (MDRD eGFR), `Bajaj_2017_nivolumab.R` (CKD-EPI eGFR, reference 90 mL/min/1.73 m^2), `NA_NA_lidocaine.R` (DDMODEL00000281; binary stratification at threshold 52.7 mL/min adding -0.319 to the GX rate constant K30 in the CRCL <= 52.7 cohort; the source `.ctl` does not state the BSA-normalisation method).
- **Notes:** The two estimation methods (MDRD/CKD-EPI vs measured CrCl) produce values in the same units and are operationally interchangeable as a covariate on clearance. Document the method explicitly in each model's `covariateData[[CRCL]]$description` so future reviewers can trace the source assay.

### CREAT (**canonical for serum creatinine**)
- **Description:** Serum creatinine concentration (baseline or time-varying).
- **Units:** umol/L or mg/dL -- document the unit used in each model via `covariateData[[CREAT]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(CREAT / ref)^exponent`.
- **Source aliases:**
  - `CRE` (umol/L, reference 70.73) -- used in `Thakre_2022_risankizumab.R`.
  - `SCR` -- common clinical-PK abbreviation.
- **Example models:** `Thakre_2022_risankizumab.R`.
- **Notes:** `CREAT` chosen over the shorter `CRE`/`SCR` as the NONMEM/clinical-PK convention that is unambiguous. Per-model reference values must be documented in `covariateData[[CREAT]]$notes`.

### ALB (**canonical for serum albumin**)
- **Description:** Serum albumin concentration.
- **Units:** g/dL or g/L -- document the unit used in each model via `covariateData[[ALB]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(ALB / ref)^exponent`.
- **Source aliases:**
  - `BALB` (baseline albumin) -- used in `Zhou_2021_belimumab.R`. Maps directly to `ALB`; baseline-vs-time-varying status documented in per-model notes.
- **Example models:** `Fasanmade_2009_infliximab.R` (g/dL, reference 4.1), `Thakre_2022_risankizumab.R` (g/L, reference 45), `Chua_2025_mirikizumab.R`, `Moein_2022_etrolizumab.R`, `Tiraboschi_2025_amlitelimab.R`, `Yamada_2025_zolbetuximab.R`, `Li_2019_abatacept.R` (g/dL, reference 4.0; the Li 2019 Methods states 'mg/dL' which is a publication typo -- see the model's `covariateData[[ALB]]$notes`), `Quartino_2019_trastuzumab.R` (g/dL, reference 4; source column `ALBU`; negative exponent -0.998 on linear CL), `Wang_2020_ontamalimab.R` (g/L, reference 39), `Zhou_2021_belimumab.R` (g/L, reference 40; baseline-only, source column `BALB`), `Okada_2025_rocatinlimab.R` (g/L, reference 44; source column `ALBU`; power exponent -1.30 on linear CL).
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
- **Example models:** `Yamada_2025_zolbetuximab.R` (mg/dL, reference 0.38; small positive exponent 0.0347 on V1), `NA_NA_lidocaine.R` (mg/dL, source column `BIL`; binary effect at threshold 0.53 mg/dL on the GX elimination rate constant K30).
- **Notes:** Hepatic-function marker. Unit varies by paper (US convention mg/dL, SI convention umol/L; 1 mg/dL ~= 17.1 umol/L). The per-model `covariateData[[TBILI]]$units` field is load-bearing.

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

### NEUT (**canonical for absolute neutrophil count**)
- **Description:** Absolute neutrophil count, typically as a baseline covariate (entered via centred-deviation `(NEUT - ref)` or power scaling `(NEUT / ref)^exponent`). Distinct from `NLR` (neutrophil-to-lymphocyte ratio): `NEUT` is the absolute count itself.
- **Units:** cells/mm^3 (equivalent to cells/uL; i.e., the same value reported in 10^9/L x 1000). Document per-model via `covariateData[[NEUT]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used in centred-deviation form `exp(coef * (NEUT - ref))`. Reference values observed: 4133 cells/mm^3 (BAST PTTE 2017 simulated cohort median; `NA_NA_tte_gompertz.R` Event 1 base hazard model).
- **Source aliases:** none.
- **Example models:** `NA_NA_tte_gompertz.R` (BAST PTTE 2017 / DDMODEL00000243 Event 1 hazard model; centred at NEUT = 4133/mm^3; coefficient -1.56e-4 on the NONMEM rescaled scale, equivalent to `exp(-1.56e-4 * (NEUT - 4133))` on the hazard).
- **Notes:** Specific scope because the BAST 2017 PTTE bundle is a teaching guiding-document with a fictional simulated population (N = 200 hypothetical patients, no real drug or indication); the centred-deviation form and reference value 4133 are tied to that specific cohort and should not be reused without re-deriving the reference from a new population's median. When a future model uses `NEUT` with general clinical-laboratory semantics (real-population hematology covariate), promote the scope to `general`. Distinct from `WBC` (total white blood cell count, of which neutrophils are the largest fraction in healthy adults) -- `NEUT` is a specific differential-count subfraction. Also distinct from `NLR` (neutrophil-to-lymphocyte ratio), which is a derived ratio. The BAST PTTE bundle uses `NEUT` as a baseline (time-fixed) covariate.

## Coagulation / hemostasis biomarkers

### VWF (**canonical for von Willebrand factor concentration**)
- **Description:** Plasma concentration (or activity) of von Willebrand factor (VWF) -- the multimeric carrier protein that binds and protects circulating factor VIII (FVIII) from proteolytic degradation and rapid clearance. Used as a covariate on FVIII (and FVIII-Fc) clearance because the vast majority (>95%) of circulating FVIII is in complex with VWF.
- **Units:** IU/dL (equivalent to % of pooled normal plasma); document per-model via `covariateData[[VWF]]$units`. Some sources report `VWF:Ag` (antigen) versus `VWF:RCo` (ristocetin cofactor activity); record which assay was used in `covariateData[[VWF]]$notes`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(VWF / ref)^exponent`. Reference values observed: 118 IU/dL (Nestorov 2014, study-population median).
- **Source aliases:** none; `VWF` is the universal abbreviation. Source papers may write `vWF` (lowercase v) or specify the assay (`VWF:Ag`).
- **Example models:** `Nestorov_2014_factorviii.R` (reference 118 IU/dL, exponent -0.343 on CL; VWF antigen).
- **Notes:** Higher VWF protects FVIII from clearance, so the exponent on CL is negative. VWF is time-varying within an individual (acute-phase response, age, blood group, etc.), but most published population PK models use baseline-only VWF when the within-subject dynamics are not characterized; document the per-model convention in `covariateData[[VWF]]$notes`.

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

### ACUTE_MED_DAYS (**canonical for baseline number of days/month of acute migraine medication use**)
- **Description:** Baseline number of days per month on which acute migraine medication (triptans or ergot compounds) was used during the 28-day run-in period prior to first dose. Enters as a piecewise-linear shift on baseline migraine or moderate-to-severe headache days in migraine exposure-response models.
- **Units:** days/month
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- piecewise-linear shift with breakpoint at 5 d/mo: contributes 0 below 5 and `slope * (ACUTE_MED_DAYS - 5)` above 5 (Fiedler-Kelly 2020). The 5-day breakpoint reflects the clinical guideline for medication-overuse headache.
- **Source aliases:** "Baseline days/month of acute medications" -- used in `FiedlerKelly_2020_fremanezumab_em.R` and `FiedlerKelly_2020_fremanezumab_cm.R`.
- **Example models:** `FiedlerKelly_2020_fremanezumab_em.R` (slope 0.438 d/d, episodic migraine), `FiedlerKelly_2020_fremanezumab_cm.R` (slope 0.460 d/d, chronic migraine).
- **Notes:** Specific scope because the variable is migraine-domain-bound. Time-fixed per subject (baseline-only). When future migraine E-R models register additional aliases or alternative breakpoints, document them per-model and consider promoting to `general`.

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

### IL6 (**canonical for serum IL-6 concentration**)
- **Description:** Baseline serum interleukin-6 concentration.
- **Units:** pg/mL
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(IL6 / ref)^exponent`.
- **Source aliases:**
  - `bIL6`, `IL6_BASE` (baseline-IL6 column variants).
- **Example models:** `Frey_2013_tocilizumab.R` (power effect on baseline target receptor).
- **Notes:** Time-fixed at baseline unless the source paper states otherwise. Time-varying IL-6 in PD models should be the model's predicted state (a `d/dt(IL6)` trajectory) rather than a covariate; use this register entry only for the baseline / observed-covariate role.

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
- **Example models:** `Netterberg_2017_docetaxel.R` (piecewise-linear effects on baseline ANC with separate low-AAG and high-AAG slopes around median 1.34 g/L; linear effect on the drug-effect slope SL via `(1 + theta * (AAG - 1.34))`).
- **Notes:** General scope because serum AAG is a routine clinical-laboratory measurement that recurs across cytotoxic-chemotherapy population-PK / PD analyses (Bruno 1996/1998 docetaxel popPK uses AAG as a CL covariate; Kloft 2006 and downstream Friberg-family myelosuppression models use it on baseline ANC and drug-effect slope). Time-fixed at baseline unless the source paper states otherwise. The Kloft 2006 piecewise-linear breakpoint at 1.34 g/L corresponds to the population median in their pooled cancer cohort; future papers may use different breakpoints, so document the per-model breakpoint in `covariateData[[AAG]]$notes`. Distinct from `CRP` (a different acute-phase reactant with different binding properties).

### IL6 (**canonical for serum interleukin-6 concentration**)
- **Description:** Serum (or plasma) interleukin-6 (IL-6) concentration. Pro-inflammatory cytokine; elevated in rheumatoid arthritis, Castleman's disease, sepsis, COVID-19, and other inflammatory conditions. Both baseline and time-varying usages are covered; document per-model in `covariateData[[IL6]]$notes` whether the column is baseline-only or time-varying.
- **Units:** pg/mL (= ng/L; the two labels are numerically equivalent).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(IL6 / ref)^exponent` or with the log-transformed form `(log(IL6 * 1000) / log(ref * 1000))^exponent` that some legacy NONMEM analyses adopt. Reference values observed: 20 pg/mL (Frey 2013 baseline; the formula `(log(IL-6 * 1000)/9.9)^exponent` is the algebraic equivalent of `(log(IL-6) / log(20))^exponent` after the constant-factor rescaling that ties the reference log to 9.9 = log(20000)).
- **Source aliases:**
  - `BLIL6` -- baseline IL-6 (used in some NONMEM control streams; canonical drops the `BL` prefix per the `EOS` / `EASI` / `AGE` convention with baseline-vs-time-varying status documented in per-model notes).
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

### CLI (**canonical for individual posthoc clearance from an upstream popPK fit**)
- **Description:** Subject-specific empirical-Bayes (posthoc) total plasma clearance from a separately published population PK model that the current PD model treats as a fixed input. Used as a per-subject (time-fixed) covariate in PD-only models that derive a per-cycle exposure metric (e.g., AUC = DOSE / CLI) without instantiating a PK ODE.
- **Units:** L/h (document per-model via `covariateData[[CLI]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- appears directly in derived-exposure expressions; not a covariate effect coefficient.
- **Source aliases:**
  - `CL` -- used in the Hansson 2013 sunitinib biomarker / TGI / fatigue PD-model family (DDMODEL00000197 and siblings, including DDMODEL00000222) and in `Schindler_2016_sunitinib.R` (DDMODEL00000221) as the posthoc CL column from the paper's upstream 2-compartment popPK fit.
- **Example models:** `Hansson_2013a_sunitinib.R` (DDMODEL00000197; typical-value reference 32.819 L/h, drawn from the bundle's simulated dataset for subject 1 -- broadly consistent with Houk et al. 2010 typical sunitinib CL); `Hansson_2013c_sunitinib.R` (DDMODEL00000222; uses a per-record `CL` column with subject-specific values 30-43 L/h in the bundle's three-subject simulated dataset); `Schindler_2016_sunitinib.R` (DDMODEL00000221; per-subject post-hoc CL fed in as the `CL` column, vignette uses 50 L/h literature-typical sunitinib CL/F per Houk 2010).
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
  - `BAS3` -- used in `Hansson_2013c_sunitinib.R` (DDMODEL00000222) as the posthoc sVEGFR-3 baseline column from the paper's upstream Hansson 2013a biomarker indirect-response fit (DDMODEL00000197).
- **Example models:** `Hansson_2013c_sunitinib.R` (DDMODEL00000222; bundle's three-subject simulated dataset reports BAS_SVEGFR3 values 42554-57365 pg/mL).
- **Notes:** Specific scope because the value is intrinsically tied to a specific upstream biomarker model (sVEGFR-3 indirect response under sunitinib in this case). The downstream fatigue model only consumes individual posthoc baseline / MRT / EC50 of the upstream biomarker; it does not re-fit them. Each model's `covariateData[[BAS_SVEGFR3]]$notes` should cite the upstream biomarker-PD source (paper or DDMORE ID) and explain how to populate the column for new simulations (typically: simulate from the upstream biomarker model to obtain individual posthoc baselines, or set every subject to the typical-value baseline for typical-trajectory simulations).

### MRT_SVEGFR3 (**canonical for individual posthoc mean residence time of soluble VEGFR-3 from an upstream PD fit**)
- **Description:** Subject-specific empirical-Bayes (posthoc) mean residence time of plasma sVEGFR-3 from a separately published indirect-response biomarker model that the current downstream PD model treats as a fixed input. Used as a per-subject (time-fixed) covariate in fatigue / adverse-event PD models that consume the upstream biomarker dynamics as data covariates without instantiating the biomarker ODE; appears as `kout3 = 1 / MRT_SVEGFR3` inside `model()`.
- **Units:** h (hours) -- document per-model via `covariateData[[MRT_SVEGFR3]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- appears directly in derived rate-constant expressions; not a covariate effect coefficient.
- **Source aliases:**
  - `MRT3` -- used in `Hansson_2013c_sunitinib.R` (DDMODEL00000222) as the posthoc sVEGFR-3 MRT column from the paper's upstream Hansson 2013a biomarker indirect-response fit (DDMODEL00000197).
- **Example models:** `Hansson_2013c_sunitinib.R` (DDMODEL00000222; bundle's three-subject simulated dataset reports MRT_SVEGFR3 values 313-408 h, broadly consistent with the Hansson 2013a typical sVEGFR-3 MRT of 401 h).
- **Notes:** Specific scope; same upstream-biomarker dependency rationale as `BAS_SVEGFR3`. The downstream fatigue model consumes the upstream MRT directly without re-fitting it.

### EC50_SVEGFR3 (**canonical for individual posthoc drug-effect EC50 on soluble VEGFR-3 from an upstream PD fit**)
- **Description:** Subject-specific empirical-Bayes (posthoc) half-maximum-effect concentration of the modelled drug on the production of sVEGFR-3 (or, depending on the upstream model parameterization, on the analogous Imax inhibition pathway) from a separately published indirect-response biomarker model that the current downstream PD model treats as a fixed input. Used as a per-subject (time-fixed) covariate in fatigue / adverse-event PD models that consume the upstream biomarker dynamics as data covariates without instantiating the biomarker ODE.
- **Units:** mg*h/L (when the per-cycle drug-exposure summary is `auc = DOSE / CLI` in mg*h/L, EC50_SVEGFR3 carries the same units; document per-model via `covariateData[[EC50_SVEGFR3]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- appears directly inside `model()` in the simple-Imax drug-effect term `eff3 = auc / (EC50_SVEGFR3 + auc)`.
- **Source aliases:**
  - `EC53` -- used in `Hansson_2013c_sunitinib.R` (DDMODEL00000222) as the posthoc sVEGFR-3 EC50 column from the paper's upstream Hansson 2013a biomarker indirect-response fit (DDMODEL00000197).
- **Example models:** `Hansson_2013c_sunitinib.R` (DDMODEL00000222; bundle's three-subject simulated dataset reports EC50_SVEGFR3 values 1.0-2.8 mg*h/L, consistent with the Hansson 2013a typical sVEGFR-3 IC50 typical value of 1.0 mg*h/L).
- **Notes:** Specific scope; same upstream-biomarker dependency rationale as `BAS_SVEGFR3`. The downstream fatigue model consumes the upstream EC50 directly without re-fitting it.

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
- **Example models:** `Bizzotto_2016_glucose.R` (driving regressor for the insulin-at-site-of-action delay).
- **Notes:** Specific scope because `INS` is meaningful only for glucose-kinetics or insulin-PD models that take plasma insulin as an exogenous regressor. For drugs that *modify* circulating insulin as a downstream effect, use a different mechanism-specific name. The DDMORE bundle's hand-rolled piecewise-linear interpolation (`I = (t-T1)/(TOBS-T1)*(INS-INS1)+INS1` with bracketing columns `iins / insn / td / tn`) is replaced in nlmixr2 by `linear(INS)` declared in `model()`; the bracketing columns are not required.

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

## Race / ethnicity

**Canonical pattern: `RACE_<GROUP>`.** Use one indicator per race/ethnicity group the source models. Reference category is the implicit 0 = all other groups; document explicitly which groups are in the reference. When the source uses composite groups (e.g., "Black or Other"), name them accordingly (`RACE_BLACK_OTHER`) and list the components in `notes`. The base `RACE_<GROUP>` indicators are scope: general; composite groupings are scope: specific because the grouping is tied to the study's analysis plan.

### RACE_BLACK (**canonical**)
- **Description:** 1 = Black / African American, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (document the actual reference groups used).
- **Source aliases:**
  - `BLACK` -- used in `Hu_2026_clesrovimab.R`, `Robbie_2012_palivizumab.R`.
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical form), `Robbie_2012_palivizumab.R`.

### RACE_WHITE (**canonical**)
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

### RACE_ASIAN (**canonical**)
- **Description:** 1 = Asian, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Source aliases:** `ASIAN` -- used in `Hu_2026_clesrovimab.R`, `Robbie_2012_palivizumab.R`, `Fau_2020_isatuximab.R`. `RAAS` (race-Asian-vs-other indicator as named in Bajaj 2017 Table 1) -- used in `Bajaj_2017_nivolumab.R`.
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical form), `Robbie_2012_palivizumab.R`, `Bajaj_2017_nivolumab.R`, `Fau_2020_isatuximab.R`.

### RACE_ASIAN_OTH (**canonical for composite Asian / American Indian / Other group**)
- **Description:** 1 = Asian, American Indian / Alaska Native, or Other race; 0 = White or Black. Composite indicator that pools the smaller-N race groups in a population dominated by White and Black subjects, with White + Black serving as the reference category. Distinct from `RACE_ASIAN_AMIND_MULTI` (Clegg 2024 grouping that includes Multiracial; pooled against a different reference) and from `RACE_BLACK_OTH` (different composite).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = White or Black (the larger-N pooled group used as the reference in the source paper).
- **Source aliases:** none formally; Frey 2013's NONMEM control stream uses the inline race classification rather than a separate named column.
- **Example models:** `Frey_2013_tocilizumab.R` (multiplicative fractional effect on the DAS28 first-order loss rate Kout: `Kout * (1 - 0.25 * RACE_ASIAN_OTH)` -- Kout is 25% lower in the Asian/AmInd/Other composite group relative to the White+Black reference).
- **Notes:** Specific scope because the composite grouping is defined by the source paper's analysis plan rather than by a uniform external standard. Do not combine with the decomposed `RACE_ASIAN`, `RACE_OTHER`, etc. indicators in the same model; the composite indicator is mutually exclusive with the decomposition. Ratified canonically on 2026-04-29 in support of the Frey 2013 tocilizumab DAS28 PKPD model.

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
- **Example models:** `Frey_2013_tocilizumab.R` (paper's "Other Asian" composite indicator on CL).
- **Notes:** Distinct from `RACE_ASIAN_AMIND_MULTI` (a 4-way composite of Asian + American Indian + Multiple Races) because the underlying paper's grouping rule is different -- `RACE_ASIAN_OTH` is a within-Asian-population sub-indicator, not a multi-race composite. Operator decision (2026-04-28): kept separate from `RACE_ASIAN` because the paper's "Other Asian" category is its own grouping, not an alias of "Asian (any)".

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

### RACE_MULTI
- **Description:** 1 = multiracial, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Source aliases:** `MULTIRACIAL` -- used in `Hu_2026_clesrovimab.R`.
- **Example models:** `Hu_2026_clesrovimab.R`.

### RACE_OTHER
- **Description:** 1 = race category "Other," 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Source aliases:**
  - `OTHER` -- used in `Robbie_2012_palivizumab.R`.
- **Example models:** `Zhu_2017_lebrikizumab.R`, `Robbie_2012_palivizumab.R`.

### RACE_HISPANIC
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
## Geographic / enrollment-country indicators

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
- **Example models:** `Chen_2022_guselkumab.R` (multiplicative effect on CL/F: 1.15^DIAB, +15% in patients with diabetes).
- **Notes:** Captures pre-existing diabetes mellitus as a comorbidity in non-diabetes-primary indications (e.g., psoriatic arthritis, psoriasis). Distinct from a primary disease-state indicator like `DIS_UC`. Type 1 vs Type 2 mellitus is not separated unless the source paper distinguishes them; in pooled-population PK analyses, the covariate is typically a single binary flag derived from medical history. Diabetic patients tend to have higher inflammation and altered IgG turnover, which can manifest as modest changes in monoclonal-antibody clearance.

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
- **Example models:** `Gandhi_2021_abatacept.R` (additive coefficient on logit-F: pJIA patients have markedly higher SC bioavailability than RA reference); `Zhong_2026_abatacept.R` (additive coefficient +3.08 on logit-F transferred verbatim from a previous internal JIA PPK model that matches Gandhi 2021's published value).
- **Notes:** Used when a population PK model pools pJIA patients with a non-pJIA reference population (e.g., Gandhi 2021: pooled adult RA + pediatric pJIA; Zhong 2026: pooled adult RA + pediatric pJIA + adult/pediatric HM) and pJIA disease/age status is tested as a PK covariate (typically on bioavailability rather than CL). Distinct from `CHILD` and `ADOLESCENT`, which are pure age-band indicators independent of indication. Scope: specific; promote to general if a third paper pools pJIA with a non-pJIA reference and the reference category remains adult RA.

### DIS_CANCER (**canonical for advanced-solid-tumor / oncology cohort indicator**)
- **Description:** 1 = patient with an advanced or metastatic solid tumor (the oncology cohort in a pooled multi-indication PK/PD analysis), 0 = non-oncology subject (healthy volunteer or non-oncology disease cohort pooled in the source analysis). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-oncology subject; the complement group is paper-defined -- typically the union of healthy volunteers and a non-oncology disease cohort such as cGVHD pooled in the source analysis).
- **Source aliases:** none; source NONMEM / Monolix control streams typically derive the indicator from a `POP` or `STUDY` categorical alongside `DIS_HV`.
- **Example models:** `Yang_2024_axatilimab.R` (multiplicative effect on baseline NCMC: `BL_NCMC x exp(1.22 x DIS_CANCER + 0.618 x DIS_HV)`; reference category cGVHD when both indicators are 0).
- **Notes:** Used together with `DIS_HV` to decompose a three-level "participant population" categorical (cGVHD reference, advanced solid tumor, healthy volunteer) into two orthogonal binary indicators. Scope: specific because the disease-pooling reference category is paper-defined (Yang 2024 reference is patients with cGVHD). Ratified canonically on 2026-04-28.

### DIS_HV (**canonical for healthy-volunteer cohort indicator**)
- **Description:** 1 = healthy volunteer (no diagnosis), 0 = patient (any diagnosis represented in the pooled cohort). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (patient subject; the complement group is the union of disease cohorts pooled in the source analysis).
- **Source aliases:** none known; healthy-volunteer indicators in source NONMEM control streams typically use ad-hoc names (e.g., `HV`, `HEALTHY`).
- **Example models:** `Nikanjam_2019_siltuximab.R` (multiplicative effects: 0.77 on CL, 0.83 on Vss; reference category is the pooled non-HV oncology cohort), `Okada_2025_rocatinlimab.R` (multiplicative shift `1 - 0.532` on Vmax when 1; reference complement is the pooled atopic-dermatitis + ulcerative-colitis + plaque-psoriasis patient cohort).
- **Notes:** Used when a population PK model pools healthy volunteers with patients across heterogeneous indications and the HV-vs-patient contrast is retained as a covariate. Scope: specific because the complement reference category is paper-defined (Nikanjam 2019 reference is "all non-HV, non-Castleman, non-SMM tumor types"; Okada 2025 reference is the pooled AD+UC+psoriasis patient cohort). Ratified canonically on 2026-04-24.
- **Example models:** `Nikanjam_2019_siltuximab.R` (multiplicative effects: 0.77 on CL, 0.83 on Vss; reference category is the pooled non-HV oncology cohort), `Yang_2024_axatilimab.R` (multiplicative effect on baseline NCMC: `BL_NCMC x exp(1.22 x DIS_CANCER + 0.618 x DIS_HV)`; reference category cGVHD).
- **Notes:** Used when a population PK model pools healthy volunteers with patients across heterogeneous indications and the HV-vs-patient contrast is retained as a covariate. Scope: specific because the complement reference category is paper-defined (Nikanjam 2019 reference is "all non-HV, non-Castleman, non-SMM tumor types"; Yang 2024 reference is patients with cGVHD). Ratified canonically on 2026-04-24.

### DIS_CASTLEMAN (**canonical for Castleman's disease indicator**)
- **Description:** 1 = Castleman's disease (multicentric or unicentric), 0 = not Castleman's disease. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-Castleman subject; the complement group is the union of healthy volunteers and other indications pooled in the source analysis).
- **Source aliases:** none known; source NONMEM control streams typically use ad-hoc names (e.g., `CD`, `CASTLEMAN`).
- **Example models:** `Nikanjam_2019_siltuximab.R` (multiplicative +24% effect on CL; no Vss effect).
- **Notes:** Castleman's disease is a lymphoproliferative disorder strongly associated with elevated IL-6 levels; it is the only FDA-approved indication for siltuximab. Scope: specific because the disease-pooling reference category is paper-defined. Ratified canonically on 2026-04-24.

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

## Infectious disease (Plasmodium / malaria)

### LNPC (**canonical for log-transformed admission Plasmodium parasitaemia**)
- **Description:** Natural logarithm of the asexual Plasmodium parasite count (parasites per microlitre of blood) at study admission. Time-fixed per subject (one value per subject, captured at enrolment).
- **Units:** log(parasites/uL)
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used with linear-deviation forms `(1 + e * (LNPC - ref))`. Reference values observed: 5.88 log(parasites/uL) (Birgersson 2019, population median in the pooled pregnant + non-pregnant Burkina Faso cohort).
- **Source aliases:** none formally; companion column `PARA` (raw asexual parasite count per microlitre) is provided alongside `LNPC` in the Birgersson 2019 NONMEM dataset but the model uses `LNPC = log(PARA)` as the active covariate.
- **Example models:** `Birgersson_2019_artesunate.R` (linear-deviation effect on relative bioavailability `F1`: `F1LNPC = 1 + e_lnpc_f * (LNPC - 5.88)`; positive coefficient `e_lnpc_f = +0.138` per unit increase in log-parasite-count, reflecting increased oral artesunate bioavailability with higher parasite burden).
- **Notes:** Disease-severity covariate specific to malaria PK models. Higher parasitaemia is a marker of more severe acute malaria infection and has been associated in the source publication with altered oral bioavailability of artesunate (presumably via gut-mucosal / first-pass effects of the febrile parasitised state). Scope: specific because the canonical reference value (5.88) is the Birgersson 2019 cohort median; future malaria-in-pregnancy or malaria-in-children PK models may legitimately reuse `LNPC` but should document their own cohort-specific reference value in `covariateData[[LNPC]]$notes`. Ratified canonically on 2026-05-07.

## Infectious disease (SARS-CoV-2 / COVID-19)

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

## Infectious disease (HIV)

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

### MM_NIGG (**canonical for non-IgG multiple myeloma immunoglobulin-type indicator**)
- **Description:** 1 = patient with non-IgG-secreting multiple myeloma (e.g., IgA, IgD, IgE, IgM, light-chain-only / Bence Jones, or non-secretory MM), 0 = patient with IgG-secreting multiple myeloma.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (IgG MM).
- **Source aliases:**
  - `Ig_type` -- used in `Fau_2020_isatuximab.R`. Values 0 / 1 with the same orientation as the canonical (1 = non-IgG MM).
- **Example models:** `Fau_2020_isatuximab.R` (exponential effect on the steady-state linear CL CLinf with coefficient -0.751, and on the time-varying-CL half-time KCL with coefficient -0.931).
- **Notes:** Within-disease (multiple-myeloma) immunoglobulin-subtype stratifier. The mechanistic rationale (Fau 2020) is that endogenous IgG monoclonal protein in IgG-MM patients competes with the therapeutic IgG mAb for FcRn-mediated salvage, raising the therapeutic mAb's catabolic clearance; non-IgG-MM patients lack that competition and exhibit lower therapeutic-mAb clearance. Distinct from the disease-state indicators (`DIS_SMM` = smoldering MM); applies only after a multiple-myeloma diagnosis is established. Scope: specific because the comparison is a within-MM stratifier rather than a cross-population indicator.
### DIS_PSORIASIS (**canonical for plaque psoriasis disease-state indicator**)
- **Description:** 1 = plaque psoriasis patient, 0 = non-psoriasis subject (e.g., atopic dermatitis, ulcerative colitis, or healthy volunteer). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-psoriasis subject; the complement group is paper-defined -- the union of other disease cohorts pooled in the source analysis).
- **Source aliases:** none known; source NONMEM control streams typically use a categorical `DIS` indicator (e.g., Okada 2025: `DIS=1` for psoriasis, `DIS=0` for healthy, `DIS=2` for UC, `DIS=3` for AD), decomposed into a binary `DIS_PSORIASIS` indicator at ingestion.
- **Example models:** `Okada_2025_rocatinlimab.R` (multiplicative shift `1 - 0.372` on linear CL when 1; reference complement is the pooled atopic dermatitis + ulcerative colitis + healthy-volunteer cohort).
- **Notes:** Used when a population PK model pools plaque-psoriasis patients with a non-psoriasis reference population and psoriasis disease status is retained as a covariate. Scope: specific because the disease-pooling reference category is paper-defined. Ratified canonically on 2026-04-27.

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
- **Units:** mm
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(TUMSZ / ref)^exponent`. Reference values observed: 41 mm (Zhou 2025); 63 mm (Budha 2023); 90 mm (Lu 2014, source reference 9 cm converted to mm).
- **Source aliases:**
  - `LDIAM` (Zhou 2025; pediatric lymphoma "linear diameter" of target lesions in mm).
  - `TMBD` (originally in cm; `TUMSZ_mm = TMBD_cm * 10`) -- used in `Lu_2014_trastuzumabemtansine.R`.
- **Example models:** `Budha_2023_tislelizumab.R` (reference 63 mm), `Lu_2014_trastuzumabemtansine.R` (reference 90 mm; source column TMBD in cm, values converted to mm on ingestion), `Zhou_2025_brentuximab.R` (reference 41 mm; source column LDIAM is the sum of linear diameters of target lesions; effect on ADC clearance only).
- **Notes:** Promoted to scope: general on 2026-04-20 as a conventional oncology baseline-tumor-size measure (RECIST for solid tumors, SPPD or sum-of-linear-diameters for lymphomas). The SPPD vs sum-of-diameters vs sum-of-linear-diameters convention is pooled onto a single column; document the per-model mixture where relevant. When the source paper reports tumor size in cm, convert to mm (the canonical unit) on data ingestion and scale the per-model reference accordingly so `(TUMSZ / ref)^exp` is numerically invariant. When a source paper specifically reports the RECIST 1.1 "sum of longest diameters" of target lesions, use the more specific `TUM_SLD` canonical instead -- `TUMSZ` remains the pooled-tumor-burden register.

### TUM_SLD (**canonical for sum of longest diameters of target lesions**)
- **Description:** Baseline sum of longest diameters of target lesions per RECIST 1.1. More specific than the pooled `TUMSZ` canonical; use `TUM_SLD` when the source paper explicitly reports "sum of longest diameters" (or "sum of lesions") as the tumor-burden metric, distinct from the pooled "sum of diameters / SPPD / sum of linear diameters" mixture covered by `TUMSZ`.
- **Units:** mm
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used with power scaling `(TUM_SLD / ref)^exponent`. Reference values observed: 70.0 mm (de Vries Schultink 2020 zenocutuzumab population median).
- **Source aliases:**
  - `SoL` / "sum of lesions" (de Vries Schultink 2020 zenocutuzumab) -- same construct, mm.
- **Example models:** `deVriesSchultink_2020_zenocutuzumab.R` (reference 70.0 mm; power exponent 0.447 on Vmax of the parallel non-linear / Michaelis-Menten clearance).
- **Notes:** Distinct from `TUMSZ` (pooled tumor-size canonical covering RECIST sum-of-diameters / SPPD / sum-of-linear-diameters); `TUM_SLD` is the precise RECIST 1.1 sum-of-longest-diameters metric. Ratified canonically on 2026-04-29 alongside the pilot bispecific extraction (de Vries Schultink 2020 zenocutuzumab). When the source paper reports tumor burden in cm, convert to mm (the canonical unit) on data ingestion and scale the per-model reference accordingly so `(TUM_SLD / ref)^exp` is numerically invariant.

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
- **Reference category:** n/a -- used with power scaling `(HER2_ECD / ref)^exponent`. Reference 25 ng/mL used in Lu 2014.
- **Source aliases:**
  - `ECD` -- used in `Lu_2014_trastuzumabemtansine.R`.
- **Example models:** `Lu_2014_trastuzumabemtansine.R` (reference 25 ng/mL; exponent 0.035 on CL).
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
- **Example models:** `Ogasawara_2020_durvalumab.R` (power effect on CL, exponent 0.0617, reference 173.8 pg/mL; time-varying; values below LOD imputed as LOD/2 = 33.55 pg/mL).
- **Notes:** Scope: specific because sPD-L1 is meaningful only for drugs targeting the PD-1/PD-L1 pathway. For other checkpoint biomarkers (e.g., soluble CTLA-4, soluble LAG-3) register new dedicated canonicals rather than reusing this one. Ratified canonically on 2026-04-26.
  - `TTYPE3` (Wang 2024; level 3 of a five-level tumor-type factor labelled "Other" in the source) -- decompose into `TUMTP_OTH = as.integer(TTYPE3 == 1)`.
- **Example models:** `Quartino_2019_trastuzumab.R` (per-group typical-value switch on linear CL; NSCLC plus a small residual group of prostate, ovarian, and other histologies), `Wang_2024_sugemalimab.R` (heterogeneous solid-tumor residual group of n = 174; exponential coefficient log(0.885) on CL and log(0.926) on Vc; NSCLC is the reference group, not part of `TUMTP_OTH`).
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
- **Description:** 1 = lymphoma (heterogeneous lymphoma pool spanning multiple lymphoma histologies -- e.g., classical Hodgkin lymphoma combined with extranodal NK/T-cell lymphoma), 0 = solid tumor or other tumor type.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = non-lymphoma tumor type (per source paper; e.g., NSCLC, GC/GEJ, ESCC, "Other" solid tumors in the Wang 2024 cohort, with NSCLC as the implicit reference when paired with the other Wang 2024 `TUMTP_*` indicators).
- **Source aliases:**
  - `TTYPE1` (Wang 2024) -- decompose into `TUMTP_LYMPH = as.integer(TTYPE1 == 1)`. The Wang 2024 source paper uses a multi-level `TTYPE` factor with levels 1 = lymphoma, 2 = lung cancer (reference), 3 = other, 4 = GCGEJ, 5 = ESCC.
- **Example models:** `Wang_2024_sugemalimab.R` (exponential coefficient log(0.877) on baseline CL and log(0.879) on Vc).
- **Notes:** Distinct from `TUMTP_CHL` (which is specifically classical Hodgkin lymphoma). Wang 2024 pools two lymphoma histologies (extranodal NK/T-cell lymphoma from CS1001-201 / NCT03595657 and classical Hodgkin lymphoma from CS1001-202 / NCT03505996) into a single lymphoma indicator; the indicator therefore captures a generic "hematologic-vs-solid-tumor" contrast rather than a histology-specific effect. When a future paper studies a single lymphoma histology distinct from cHL, register a more specific canonical (e.g., `TUMTP_ENKTL`, `TUMTP_NHL`) rather than overloading this one. Document the per-paper histology composition in `covariateData[[TUMTP_LYMPH]]$notes`.

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

## Laboratory / disease-activity

### ALBR
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

## Concomitant IBD medications

### CONMED_AZA (**canonical for concomitant azathioprine**)
- **Description:** 1 = on concomitant azathioprine at the PK observation, 0 = not on azathioprine.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant azathioprine).
- **Source aliases:** `AZA` -- used in `Rosario_2015_vedolizumab.R`.
- **Example models:** `Rosario_2015_vedolizumab.R` (power-form on CLL: `CLL * 0.998^CONMED_AZA`; effect ~= null).
- **Notes:** Thiopurine immunomodulator used as maintenance therapy in IBD. Standard convention is baseline-use-only, but time-varying use is permitted; document per-model.

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

### CONMED_AMINO (**canonical for concomitant aminosalicylate therapy**)
- **Description:** 1 = on concomitant aminosalicylate (5-aminosalicylic acid / mesalamine / mesalazine / olsalazine / sulfasalazine etc.) therapy at the PK observation, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant aminosalicylate).
- **Source aliases:** `AMINO` -- used in `Rosario_2015_vedolizumab.R`.
- **Example models:** `Rosario_2015_vedolizumab.R` (power-form on CLL: `CLL * 1.02^CONMED_AMINO`).
- **Notes:** Covers the full aminosalicylate class (5-ASA is the single active moiety shared by most agents); use `CONMED_AMINO` rather than `CONMED_5ASA` unless the source paper explicitly restricts the indicator to 5-ASA monotherapy.

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

### BLPHYVAS
- **Description:** Baseline Physician's Global Assessment of Disease Activity, 100-mm visual analogue scale (0 = no disease activity, 100 = maximum). Time-fixed per subject.
- **Units:** mm (0-100 VAS)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a -- used as a power term `(BLPHYVAS / <ref>)^exponent`. Reference 66 used in Ma 2020.
- **Source aliases:** none.
- **Example models:** `Ma_2020_sarilumab_das28crp.R`.
- **Notes:** One of the components of the DAS28 composite score; in Ma 2020 it appears as a baseline covariate on the DAS28-CRP disease-activity BASE rather than on the score itself. Applicable to any rheumatology model where baseline physician-assessed disease activity is used as a PK/PD covariate.

### BLHAQ
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

## Concomitant / prior medication

### PRICORT
- **Description:** 1 = patient received systemic corticosteroid treatment prior to study entry, 0 = no prior corticosteroid use. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no prior corticosteroid use).
- **Source aliases:** none.
- **Example models:** `Ma_2020_sarilumab_das28crp.R` (multiplicative on DAS28-CRP Kout: `Kout * theta^PRICORT`); `Ma_2020_sarilumab_anc.R` (power-form on Emax: `Emax * 0.819^PRICORT`).
- **Notes:** Ma 2020 applies it as a multiplicative effect of the form `param * theta^PRICORT` in both DAS28-CRP and ANC PD models. Generally applicable clinical-history indicator.

### STEROID (**canonical for baseline/concomitant systemic corticosteroid use**)
- **Description:** 1 = patient on systemic corticosteroid therapy at baseline (typically continued as concomitant medication during the study), 0 = no baseline corticosteroid use. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no baseline corticosteroid use).
- **Source aliases:**
  - `BSTEROID` -- used in `Narwal_2013_sifalimumab.R` and `Zheng_2016_sifalimumab.R`.
- **Example models:** `Narwal_2013_sifalimumab.R` (multiplicative on CL: `CL * (1 + 0.195 * STEROID)`), `Zheng_2016_sifalimumab.R` (multiplicative on CL `(1 + 0.11 * STEROID)` and on V1 `(1 - 0.09 * STEROID)` in the SLE phase IIb cohort, which was ~85% steroid-treated at baseline).
- **Notes:** Distinct from `PRICORT`, which is strictly a prior (pre-study) indicator. `STEROID` captures concurrent corticosteroid use at / from study baseline in diseases where background steroid use is standard of care (SLE, severe asthma, etc.). When a future paper needs the two jointly, both can coexist on the same subject. The name `STEROID_BL` was used as an alias in earlier register drafts and is retired; use `STEROID` for all future models.

### COADMIN_IPI_3Q3W (**canonical for nivolumab + ipilimumab 3 mg/kg q3w combination indicator**)
- **Description:** 1 = subject is receiving nivolumab in combination with ipilimumab 3 mg/kg every 3 weeks (4-dose induction); 0 = otherwise. Encodes the high-intensity ipilimumab combination regimen as a study-design covariate on nivolumab CL.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (any non-3Q3W regimen -- monotherapy, chemotherapy combination, or another ipilimumab schedule).
- **Source aliases:**
  - `IPI3Q3W` -- used in `Zhang_2019_nivolumab.R`.
- **Example models:** `Zhang_2019_nivolumab.R` (exponential effect on baseline CL: `exp(0.227)` ~= 1.25 fold increase relative to monotherapy).
- **Notes:** Paired with `COADMIN_IPI_1Q6W`; both indicators can coexist in one population, but a single subject has at most one set to 1 in the Zhang 2019 cohort. The remaining ipilimumab schedules (1 mg/kg q3w x 4 induction, 1 mg/kg q12w) had no statistically significant effect on nivolumab CL and were therefore folded into the reference (0) group along with monotherapy, leaving only IPI3Q3W and IPI1Q6W as named non-reference indicators.

### COADMIN_IPI_1Q6W (**canonical for nivolumab + ipilimumab 1 mg/kg q6w combination indicator**)
- **Description:** 1 = subject is receiving nivolumab in combination with ipilimumab 1 mg/kg every 6 weeks (continuous maintenance); 0 = otherwise.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (any non-1Q6W regimen -- monotherapy, chemotherapy combination, or another ipilimumab schedule).
- **Source aliases:**
  - `IPI1Q6W` -- used in `Zhang_2019_nivolumab.R`.
- **Example models:** `Zhang_2019_nivolumab.R` (exponential effect on baseline CL: `exp(0.159)` ~= 1.17 fold increase relative to monotherapy).
- **Notes:** Paired with `COADMIN_IPI_3Q3W`. See the COADMIN_IPI_3Q3W note for how the other ipilimumab schedules collapse into the reference group.

### COADMIN_CHEMO (**canonical for anti-PD-(L)1 mAb + chemotherapy combination indicator**)
- **Description:** 1 = subject is receiving an anti-PD-(L)1 monoclonal antibody in combination with platinum-based chemotherapy (gemcitabine + cisplatin, pemetrexed + cisplatin, paclitaxel + carboplatin, or platinum-doublet); 0 = otherwise. Encodes chemotherapy coadministration as a study-design covariate on the antibody's CL.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no chemotherapy coadministration -- monotherapy or, where applicable, a non-chemotherapy combination such as anti-PD-1 + anti-CTLA-4).
- **Source aliases:**
  - `CHEMO` -- used in `Zhang_2019_nivolumab.R`.
  - `MONOTR` -- used in `Kuchimanchi_2024_dostarlimab.R` (the paper's structural-equation indicator for monotherapy; the canonical column carries the inverse value, i.e. `COADMIN_CHEMO = 1 - MONOTR`, so the canonical column is 1 when the patient is on combo-chemotherapy).
- **Example models:** `Zhang_2019_nivolumab.R` (exponential effect on baseline CL: `exp(-0.104)` ~= 0.90 fold, i.e. ~9.7% lower CL relative to monotherapy), `Kuchimanchi_2024_dostarlimab.R` (multiplicative effect on baseline CL: `1 - 0.0779` = 0.922, i.e. 7.79% lower CL on dostarlimab + carboplatin/paclitaxel relative to dostarlimab monotherapy).
- **Notes:** Promoted from specific to general scope on 2026-04-27 after the Kuchimanchi 2024 dostarlimab + carboplatin/paclitaxel analysis ratified the same pooling convention (any chemotherapy backbone collapsed into a single binary indicator). The two papers use different functional forms for the effect on CL -- Zhang 2019 uses `exp(theta * COADMIN_CHEMO)` (exponential) and Kuchimanchi 2024 uses `(1 + theta * COADMIN_CHEMO)` (multiplicative); these are different parameterisations of the same underlying study-design indicator and the canonical column meaning is unchanged. Document the per-model functional form in `covariateData[[COADMIN_CHEMO]]$notes`.

### COADMIN_IPI_ANY (**canonical for any-ipilimumab-coadministration indicator**)
- **Description:** 1 = subject is receiving nivolumab in combination with any ipilimumab regimen (regardless of dose or schedule); 0 = nivolumab monotherapy or nivolumab + chemotherapy. Encodes the "is there ipilimumab in the regimen" question as a single binary covariate, distinct from the regimen-specific COADMIN_IPI_3Q3W and COADMIN_IPI_1Q6W indicators above.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no ipilimumab coadministration).
- **Source aliases:**
  - `IPICO` -- used in `Zhang_2019_nivolumab.R`.
- **Example models:** `Zhang_2019_nivolumab.R` (additive effect on the time-varying-CL Emax parameter: Emax += -0.0668 when COADMIN_IPI_ANY = 1).
- **Notes:** Logically the union of the regimen-specific indicators (COADMIN_IPI_3Q3W, COADMIN_IPI_1Q6W, plus the unmodeled 1 mg/kg q3wx4 and 1 mg/kg q12w schedules). Zhang 2019 uses it on the *time-varying* Emax (a different structural parameter from baseline CL), which is why it coexists with the regimen-specific indicators on baseline CL rather than substituting for them.

### COADMIN_AVD (**canonical for brentuximab vedotin + AVD (adriamycin/doxorubicin, vinblastine, dacarbazine) combination indicator**)
- **Description:** 1 = subject is receiving brentuximab vedotin in combination with the AVD chemotherapy backbone (adriamycin a.k.a. doxorubicin, vinblastine, dacarbazine) for newly diagnosed advanced-stage Hodgkin lymphoma; 0 = otherwise (single-agent brentuximab vedotin). Encodes the A+AVD frontline regimen as a study-design covariate on ADC clearance.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (single-agent brentuximab vedotin -- no AVD coadministration).
- **Source aliases:**
  - `DOX` -- used in `Zhou_2025_brentuximab.R` (the NONMEM dataset uses the doxorubicin-administration flag as the AVD-coadministration indicator because doxorubicin is given on the same days as the other AVD agents in this regimen).
- **Example models:** `Zhou_2025_brentuximab.R` (power-form effect on ADC clearance: `CL * 2.12^COADMIN_AVD` -- ADC clearance is ~2.1-fold higher under A+AVD vs single-agent BV).
- **Notes:** Distinct from `COADMIN_CHEMO` (which is nivolumab + platinum-based chemotherapy). The A+AVD regimen is the standard chemotherapy backbone for frontline classical Hodgkin lymphoma; promote to general scope if a second BV paper reports the same A+AVD covariate effect with a comparable encoding.

### COADMIN_SPART (**canonical for spartalizumab (PDR001, anti-PD-1) coadministration indicator**)
- **Description:** 1 = the analyzed therapeutic mAb is coadministered with spartalizumab (PDR001, anti-PD-1 IgG4), 0 = no spartalizumab coadministration. Time-fixed per subject in source analyses to date.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no spartalizumab coadministration -- monotherapy or combination with non-spartalizumab agents such as a hypomethylating agent).
- **Source aliases:**
  - `HASPDR` -- used in `Xu_2023_MBG453.R` (Monolix supplement Appendix S2; the source describes the column as "this patient HAS received PDR001 [spartalizumab, anti PD-1 mAb]").
- **Example models:** `Xu_2023_MBG453.R` (exponential effect on CL: `exp(0.0194 * COADMIN_SPART)`; not statistically significant in the full covariate model but retained because Xu 2023 used the full-covariate-model approach).
- **Notes:** Parallels `COMBO_NIVO` (ipilimumab + nivolumab) and `COMBO_DURVA` (durvalumab combinations) but for spartalizumab. Promote to general scope if a second paper reports a spartalizumab-coadministration covariate with comparable encoding.

### STATIN (**canonical for concomitant statin (HMG-CoA reductase inhibitor) therapy**)
- **Description:** 1 = patient coadministered a statin (HMG-CoA reductase inhibitor) during the study, 0 = no statin coadministration.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no statin coadministration).
- **Source aliases:** none known.
- **Example models:** `Martinez_2019_alirocumab.R` (additive effect on linear clearance CLL: `CLL = TVCLL + COV1*(WT-82.9) + COV2*STATIN`; +0.00644 L/h when statin is coadministered).
- **Notes:** Per-model `covariateData[[STATIN]]$notes` must document which statins and dose thresholds are included in the "STATIN = 1" category, since inclusion criteria vary by study. Martinez 2019 codes STATIN = 1 for coadministration of rosuvastatin (< 20 mg/day), atorvastatin (< 40 mg/day), or simvastatin (any dose); other statin regimens are coded as 0.

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

### ADA_POS (**canonical**)
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

### PRIOR_TNF (**canonical**)
- **Description:** 1 = subject previously treated with an anti-TNF (tumor necrosis factor) inhibitor, 0 = TNF-naive.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (TNF-naive).
- **Source aliases:**
  - `PRIORTNF` (all caps, no underscore) -- acceptable alternative spelling.
- **Example models:** `Moein_2022_etrolizumab.R` (multiplicative fractional effect on CL, +4.9%).
- **Notes:** Use when the source paper reports a binary "prior anti-TNF inhibitor" covariate on any PK parameter. Generally applicable across RA/PsA/IBD/axSpA biologic PK models.

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

### DISEXT_EP (**canonical for extensive colitis / pancolitis indicator**)
- **Description:** 1 = extensive colitis or pancolitis disease extension, 0 = otherwise (any non-extensive disease extension, e.g. left-sided colitis or proctitis).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0. In papers that decompose the disease-extension categorical into both `DISEXT_EP` and `DISEXT_OTHER`, `DISEXT_EP = 0 AND DISEXT_OTHER = 0` corresponds to the left-sided-colitis reference group; in papers that use a single binary indicator for extensive colitis, the reference is pooled non-extensive (left-sided + any other extension).
- **Source aliases:**
  - `EXTCOL` -- used in `Faelens_2021_infliximab.R` (binary 0/1 for extensive colitis at baseline; no separate "other" category).
  - Derived from a multi-level `DISEXT` column in the source (levels: left-sided colitis, extensive/pancolitis, other): `DISEXT_EP = as.integer(DISEXT == "extensive/pancolitis")`.
- **Example models:** `Moein_2022_etrolizumab.R` (paired with `DISEXT_OTHER`; multiplicative effect on CL, +8.2% vs. left-sided colitis); `Faelens_2021_infliximab.R` (single-binary encoding; multiplicative fold-change on V of 1.25 when DISEXT_EP = 1).
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

## Concomitant lipid-lowering medication

### STATIN_MONO
- **Description:** 1 = patient is on a statin and no other lipid-lowering comedication (statin monotherapy), 0 = not on statin monotherapy (either on no lipid-lowering therapy or on a multi-drug lipid-lowering combination such as statin + ezetimibe).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (not on statin monotherapy).
- **Source aliases:**
  - Derived from a statin-identifier column in the source (any of atorvastatin, rosuvastatin, simvastatin, lovastatin, pravastatin, pitavastatin, fluvastatin) with an AND over "no other lipid-lowering comedication".
- **Example models:** `Kuchimanchi_2018_evolocumab.R` (multiplicative effect 1.13 on Vmax: `Vmax * 1.13^STATIN_MONO`).
- **Notes:** Scope: specific because Kuchimanchi 2018 narrowly defines the statin covariate as monotherapy only ("patients on a statin only and no other comedication"). Mutually compatible with `EZE`: a subject on statin+ezetimibe has `STATIN_MONO = 0` and `EZE = 1`; a subject on statin alone has `STATIN_MONO = 1` and `EZE = 0`; a subject on no lipid-lowering therapy has both 0. Future popPK/PD models that adopt a broader "any statin" definition should register a separate `STATIN` or `CONMED_STATIN` canonical rather than reusing this name.

### EZE
- **Description:** 1 = patient is taking ezetimibe (with or without other lipid-lowering comedication), 0 = not on ezetimibe.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (not on ezetimibe).
- **Source aliases:**
  - Derived from an ezetimibe-identifier column in the source.
- **Example models:** `Kuchimanchi_2018_evolocumab.R` (multiplicative effect 1.20 on Vmax: `Vmax * 1.20^EZE`; labeled "Statin + ezetimibe exponent" in Kuchimanchi 2018 Table 3 because ~99% of ezetimibe users in the dataset were also on a statin, so the effect effectively captures combination therapy).
- **Notes:** Scope: specific because Kuchimanchi 2018 interprets the ezetimibe indicator as a combination-therapy marker rather than a pure ezetimibe effect. Future popPK/PD models with cleaner ezetimibe separation should add themselves here or register a more specific canonical.

## Hypercholesterolemia biomarkers

### PCSK9
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

## Lifestyle / medical history

### SMOKE
- **Description:** 1 = current smoker at baseline, 0 = non-smoker.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (non-smoker).
- **Source aliases:**
  - `Smoking` (case-insensitive) -- used in `Ma_2020_sarilumab_anc.R`.
- **Example models:** `Ma_2020_sarilumab_anc.R` (power-form on baseline ANC: `BASE * 1.15^SMOKE`).
- **Notes:** Baseline-only indicator; does not track within-study smoking-cessation changes. Use this two-level (current vs non-smoker) encoding when the source paper does not split former and never smokers. When the source uses a 3-level smoking-status categorical (never / former / current), use the paired `SMOKE_CURRENT` + `SMOKE_NEVER` indicators below instead -- the 3-level encoding cannot be reduced to a single `SMOKE` column without losing information.

### SMOKE_CURRENT
- **Description:** 1 = current smoker at baseline, 0 otherwise (former or never smoker). Paired with `SMOKE_NEVER` to encode a 3-level smoking-status categorical with former smoker as the implicit reference (both indicators = 0).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (former smoker, when paired with `SMOKE_NEVER` = 0). The pairing follows the `RACE_<GROUP>` convention for paired indicators.
- **Source aliases:**
  - `Smoking status = Current` / `SMOK = 2` (case-insensitive) -- derived from a 3-level smoking-status column.
- **Example models:** `Hwang_2023_monalizumab.R` (proportional-shift effect on V1: `(1 + 0.0484)^SMOKE_CURRENT`; reference category former smoker).
- **Notes:** Baseline-only indicator. See also `SMOKE_NEVER` (paired indicator) and `SMOKE` (binary current-vs-non-smoker encoding when the source paper does not split former vs never).

### SMOKE_NEVER
- **Description:** 1 = never smoker at baseline, 0 otherwise (former or current smoker). Paired with `SMOKE_CURRENT` to encode a 3-level smoking-status categorical with former smoker as the implicit reference (both indicators = 0).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (former smoker, when paired with `SMOKE_CURRENT` = 0). The pairing follows the `RACE_<GROUP>` convention for paired indicators.
- **Source aliases:**
  - `Smoking status = Never` / `SMOK = 0` (case-insensitive) -- derived from a 3-level smoking-status column.
- **Example models:** `Hwang_2023_monalizumab.R` (proportional-shift effect on V1: `(1 - 0.141)^SMOKE_NEVER`; reference category former smoker).
- **Notes:** Baseline-only indicator. See also `SMOKE_CURRENT` (paired indicator) and `SMOKE` (binary current-vs-non-smoker encoding when the source paper does not split former vs never).

### SMOKE_CURRENT
- **Description:** 1 = subject is a current smoker at baseline, 0 = not. Paired with `SMOKE_NEVER`; both 0 = former smoker. Reference category for the paired pair is "former smoker" so the two indicators capture a 3-level smoking-status covariate without an intercept clash.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (not current smoker).
- **Source aliases:** none known; paper-specific raw forms include `SMOKE = 'CURRENT'` strings or one-hot pivoted columns.
- **Example models:** (forthcoming).
- **Notes:** Use this paired-indicator form when the source paper reports a 3-level smoking covariate (current / former / never). When the paper reports only a 2-level current-vs-not covariate, use the simpler binary `SMOKE` instead (already registered above).

### SMOKE_NEVER
- **Description:** 1 = subject is a never-smoker at baseline, 0 = not. Paired with `SMOKE_CURRENT`; both 0 = former smoker.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (not never-smoker).
- **Source aliases:** none known.
- **Example models:** (forthcoming).
- **Notes:** Reference category for the paired pair is "former smoker." See `SMOKE_CURRENT` for the 2-level alternative.

### HSCT_URD_7OF8 (**canonical for unrelated-donor 7/8 HLA-match indicator**)
- **Description:** 1 = subject's hematopoietic stem-cell transplant donor was unrelated and HLA-matched at 7 of 8 loci (HLA-A, -B, -C, -DRB1), 0 = otherwise. Paired with `HSCT_URD_8OF8`; reference category (both 0) = matched-related donor or autologous transplant.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (matched-related or auto-transplant).
- **Source aliases:** none known.
- **Example models:** `Zhong_2026_abatacept.R` (pediatric aGvHD prophylaxis; donor-match covariate on PK).
- **Notes:** Specific scope because the donor-typing schema is transplant-protocol-dependent (8/8 loci, 10/10 loci, KIR matching, etc.). When a paper uses a different match-resolution definition (e.g. 9/10 vs 10/10), register a new canonical name rather than reusing `HSCT_URD_7OF8`.

### HSCT_URD_8OF8 (**canonical for unrelated-donor 8/8 HLA-match indicator**)
- **Description:** 1 = subject's hematopoietic stem-cell transplant donor was unrelated and HLA-matched at 8 of 8 loci, 0 = otherwise. Paired with `HSCT_URD_7OF8`.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (matched-related or auto-transplant).
- **Source aliases:** none known.
- **Example models:** `Zhong_2026_abatacept.R`.
- **Notes:** See `HSCT_URD_7OF8` for the paired-pair convention and the rationale for specific-scope.

### PAIN (**canonical for baseline pain score**)
- **Description:** Baseline patient-reported pain score on a paper-specified scale (NRS 0-10, VAS 0-100, BPI, or paper-specific scale). Document the scale used in `covariateData[[PAIN]]$notes` for every model that uses this covariate.
- **Units:** (paper-specified -- typically dimensionless score)
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- usually used with a linear-deviation form `1 + e * (PAIN - ref)` or a power form `(PAIN / ref)^e`. Document the reference value in `covariateData[[PAIN]]$notes`.
- **Source aliases:** none known.
- **Example models:** (forthcoming).
- **Notes:** Specific scope because the underlying scale varies across papers; if a future paper uses the same scale, extend the Example models list rather than register a new canonical name. If two papers use different scales for the same conceptual covariate, prefer registering a scale-suffixed name (`PAIN_NRS`, `PAIN_VAS`) so a downstream consumer can convert between scales.

## Formulation / assay / study

### ROUTE_IV
- **Description:** 1 = subject received intravenous (IV) administration, 0 = subcutaneous (SC) administration. Per-subject (study-fixed) covariate flagging the dosing route when a population analysis pools cohorts that differ by route, with covariate effects on PK parameters that capture route-specific disposition behaviour.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (SC).
- **Source aliases:**
  - "Admin route = IV" (categorical effect column in Yu 2022 covariate equations).
  - `IV` -- used in `Zierhut_2008_osteoprotegerin.R` (DDMODEL00000233 `dataObj` column flagging IV vs SC cohort, switching the PK observation residual SD between `CcpropSdIV` and `CcpropSdSC`).
- **Example models:**
  - `Yu_2022_ofatumumab.R` (exponential effect on R0, CL, Q, ksyninf).
  - `Zierhut_2008_osteoprotegerin.R` (per-subject indicator switching the PK observation residual SD between the IV cohort (`CcpropSdIV`) and the SC cohort (`CcpropSdSC`)).
- **Notes:** This is the per-subject covariate-equation indicator, distinct from the dosing-event `cmt` column that names the target compartment. When simulating, set `ROUTE_IV = 1` for IV cohorts and dose into the central compartment; set `ROUTE_IV = 0` for SC cohorts and dose into the depot. Scope: specific because the set of parameters that differ by route is paper-specific (Yu 2022 carries route-specific exponential effects on disposition parameters; Zierhut 2008 carries a route-specific PK observation residual SD).

### DEVICE_AI
- **Description:** 1 = subject's SC dose delivered via autoinjector (AI), 0 = prefilled syringe (PFS). Per-subject (study-fixed) covariate flagging the SC delivery device when a model carries device-specific PK effects.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (PFS).
- **Source aliases:** "Formulation = AI" (categorical effect column in Yu 2022 covariate equations).
- **Example models:** `Yu_2022_ofatumumab.R` (exponential effect on k_e(P) and R0).
- **Notes:** Set to 0 (PFS reference) for IV subjects, since the device is undefined for IV; the IV-specific effects are captured by `ROUTE_IV` instead. Scope: specific because the AI/PFS contrast and which parameters it affects depend on the study's device-comparison design.

### STUDY_APLIOS
- **Description:** 1 = subject enrolled in the APLIOS bioequivalence study (NCT03560739; phase 2; ofatumumab AI vs PFS in RMS), 0 = other study in the Yu 2022 pooled analysis.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-APLIOS studies: OMS115102, MIRROR, ASCLEPIOS I, ASCLEPIOS II).
- **Source aliases:** "Study = APLIOS" (categorical effect column in Yu 2022 covariate equations).
- **Example models:** `Yu_2022_ofatumumab.R` (exponential effect on Emax of B cell lysis).
- **Notes:** Captures a between-study shift in the maximum B-cell lysis stimulatory effect not explained by the other covariates in the final model.

### STUDY_MIRROR
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

### FED
- **Description:** 1 = fed state at dosing, 0 = fasted.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (fasted).
- **Example models:** `Kyhl_2016_nalmefene.R`.

### TABLET
- **Description:** 1 = tablet formulation, 0 = solution.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (solution).
- **Example models:** `Kyhl_2016_nalmefene.R`.
- **Notes:** Scoped specific because the "tablet vs solution" contrast is tied to Kyhl 2016's formulation-comparison design. Future formulation-comparison models should either add themselves here or register a more general `FORM_TABLET` entry.

### FDC (**canonical for fixed-dose-combination antitubercular formulation indicator**)
- **Description:** 1 = subject received the rifampicin-containing fixed-dose-combination (FDC) antitubercular product (rifampicin co-formulated with isoniazid, pyrazinamide, and optionally ethambutol in a single tablet); 0 = subject received the same drugs as separate single-drug tablets ("SDC", separate-drug-combination). Per-subject (regimen-fixed) categorical covariate flagging the formulation when a population analysis pools FDC and SDC arms and tests formulation as a covariate on absorption / disposition parameters.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 1 (FDC; the most-common formulation in the Wilkins 2008 cohort and the typical-value reference for `lmtt` / `lcl`).
- **Source aliases:**
  - `FDC` -- used in `Wilkins_2008_rifampicin.R` (DDMODEL00000280 NMTRAN `$INPUT` column; values 0 / 1 with the same orientation as the canonical, 1 = FDC).
- **Example models:** `Wilkins_2008_rifampicin.R` (multiplicative `(1 + e_fdc0_mtt * (1 - FDC))` shift on MTT and `(1 + e_fdc0_cl * (1 - FDC))` shift on CL; SDC subjects (FDC = 0) had 104% longer MTT and 23.6% higher CL than the FDC = 1 reference per Wilkins 2008 final estimates).
- **Notes:** Specific scope because the FDC vs SDC contrast is tied to the antitubercular fixed-dose-combination context (rifampicin + isoniazid + pyrazinamide +/- ethambutol). Future antitubercular-FDC models should extend the example list rather than register a new canonical. Distinct from the generic `TABLET` (Kyhl 2016 nalmefene tablet vs solution) and `FORM_*` (drug-product-version) indicators because FDC vs SDC compares two tablet products, not a tablet vs a non-tablet, and the underlying mechanism is co-formulation-driven absorption rather than drug-product manufacturing.

### RIA_ASSAY
- **Description:** 1 = radioimmunoassay; 0 = LC-MS/MS.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Example models:** `Kyhl_2016_nalmefene.R`.
- **Notes:** Switches the additive residual-error magnitude.

### FORM_NS0
- **Description:** 1 = NS0 cell-line formulation, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Example models:** `Zhu_2017_lebrikizumab.R`.

### FORM_CHO_PHASE2
- **Description:** 1 = CHO Phase 2 formulation, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Example models:** `Zhu_2017_lebrikizumab.R`.

### COMB_EOX
- **Description:** 1 = concomitant epirubicin + oxaliplatin + capecitabine (EOX) chemotherapy backbone, 0 = other backbone (e.g., mFOLFOX6, CAPOX, or single-agent).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-EOX backbone).
- **Source aliases:** `COMB` (used by Yamada 2025 Table 1 with the EOX level coded as the non-reference category; renamed to `COMB_EOX` to preserve the semantic meaning of the 1-level).
- **Example models:** `Yamada_2025_zolbetuximab.R` (fractional effect on V1).
- **Notes:** Disease-backbone indicator. If a future model needs more backbone categories, encode each as its own indicator (`COMB_CAPOX`, `COMB_FOLFOX`, ...) with a single reference group.

### FORM_P2F2
- **Description:** 1 = isatuximab P2F2 drug material (intended commercial / phase III material, used in the EFC14335 / ICARIA-MM study), 0 = P1F1 drug material (early-phase material used in TED10893 / TED14154 / TCD14079).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (P1F1).
- **Source aliases:**
  - `Drug_mat` -- used in `Fau_2020_isatuximab.R`. Values 0 / 1 with the same orientation as the canonical (1 = P2F2 / commercial-bound material).
- **Example models:** `Fau_2020_isatuximab.R` (exponential effect on Vc with coefficient -0.137; P2F2 patients had ~13% lower Vc than P1F1).
- **Notes:** Phase III / commercial-bound formulation indicator for isatuximab; the FORM_* family stays scope-specific per nlmixr2lib policy that drug-product-version indicators are kept model-specific unless they generalize across multiple drugs. Set to 1 to simulate the marketed material.

### FORM_DP2
- **Description:** 1 = sarilumab drug product 2 formulation (used in some phase I studies and the dose-ranging phase II study), 0 = other drug product (DP1 or DP3; DP3 is the commercial formulation).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (DP1 or DP3).
- **Source aliases:**
  - `DP2` -- used in `Xu_2019_sarilumab.R`.
- **Example models:** `Xu_2019_sarilumab.R`.
- **Notes:** Affects both CLO/F (1.30x multiplier) and Ka (0.663x multiplier) in Xu 2019. Set to 0 for routine commercial-formulation simulation.

### dilution
- **Description:** 1 = drug diluted (Soehoel 2022 study D2213C00001), 0 = not diluted.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Example models:** `Soehoel_2022_tralokinumab.R`.
- **Notes:** Lower-case preserved from source; future models should rename to `DILUTION`. Kept as alias here to match existing file.

### nonECZTRA
- **Description:** 1 = not the ECZTRA trial; 0 = ECZTRA.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Example models:** `Soehoel_2022_tralokinumab.R`.
- **Notes:** Mixed case preserved from source; future models should rename to `NON_ECZTRA` or `STUDY_NON_ECZTRA`.

### SEASON2
- **Description:** 1 = second RSV season at dosing, 0 = first RSV season.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Example models:** `Clegg_2024_nirsevimab.R`.
- **Notes:** Study-specific but semantically general (second-exposure indicator). Promote to general if a second RSV-season model adopts the same semantics.

### COHDOSE
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
- **Example models:** `Zheng_2016_sifalimumab.R` (power effect on V1 with exponent 0.06), `Castro-Surez_2020_nimotuzumab.R` (binary-indicator usage `(DOSE == 50)` applying a 53 % decrease in V1 for the 50 mg cohort), `Hansson_2013a_sunitinib.R` (DDMODEL00000197; time-varying record-level dose feeding `AUC = DOSE / CLI`), `Schindler_2016_sunitinib.R` (DDMODEL00000221; same `AUC = DOSE / CLI` form, with the daily-dose column toggling between 50 mg/day on-cycle and 0 on dose-holiday records), `Girard_2012_pimasertib.R` (linear coefficient on the dropout-hazard log-rate: `exp(beta * DOSE)` Weibull multiplier; per-subject daily dose, observed range 1-255 mg).
- **Notes:** Distinct from `DOSE_70MG` (binary indicator for a specific dose group in a trinary-dose design) and from the rxode2/nlmixr2 event column `amt` (which carries the administered dose at dose events). For use case (a), the values are typically time-fixed per subject; for use case (b), the values are time-varying with on/off cycling -- for sunitinib 4-weeks-on / 2-weeks-off cycling, set `DOSE = nominal_daily_mg` (e.g., 50) during on-cycles and 0 during off-cycles or for the placebo arm. Per-model `covariateData[[DOSE]]$notes` should state which use case applies.

### CD
- **Description:** Time-varying cumulative cladribine dose (mg total dose, not body-weight-normalized) administered to each subject up to the current observation time. Stays at zero during the placebo arm and during the pre-dose baseline phase, rises stepwise across the dosing schedule (cladribine is given as short oral pulses), and remains constant between dose events.
- **Units:** mg
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a -- used inside `EXPS = CD * 104.5 / CRL` as an exposure surrogate driving an Emax-style symptomatic effect on disease progression, not as a power-form covariate.
- **Source aliases:** `CD` -- column name used in the DDMODEL00000223 input dataset (`Simulated_Novakovic_2016_multiplesclerosis_cladribine_irt.csv`).
- **Example models:** `Novakovic_2017_cladribine.R`.
- **Notes:** Distinct from `DOSE` (per-subject assigned dose level, time-fixed) and `COHDOSE` (mg/kg cohort label, time-fixed). `CD` is the *cumulative* dose accrued at each timepoint, supplied as a time-varying covariate column rather than via dosing events because the Novakovic 2017 model does not carry an explicit cladribine-PK compartment. Scope: specific because the constant 104.5 inside the exposure-surrogate equation is hard-coded for cladribine in the source.

### TRT
- **Description:** Treatment-cohort indicator used in the Novakovic 2017 cladribine IRT model. 0 = placebo, 1 = cladribine 3.5 mg/kg cumulative-dose cohort, 2 = cladribine 5.25 mg/kg cumulative-dose cohort.
- **Units:** (categorical)
- **Type:** categorical
- **Scope:** specific
- **Reference category:** 0 (placebo).
- **Source aliases:** `TRT` -- column name used in the DDMODEL00000223 input dataset (`Simulated_Novakovic_2016_multiplesclerosis_cladribine_irt.csv`).
- **Example models:** `Novakovic_2017_cladribine.R`.
- **Notes:** Gates the symptomatic and protective drug-effect terms via `TRT >= 1 && t > 0`; the categorical level (1 vs 2) is informational because the dose-response is driven by the time-varying `CD` covariate and the per-cohort dosing schedule, not by `TRT` itself. Scope: specific because the cohort labelling (3.5 vs 5.25 mg/kg cumulative dose over 2 years) is tied to the CLARITY-program cladribine dosing schedule. Future models that need a generic on-treatment indicator should register a new canonical name (e.g., `ON_TREATMENT`) rather than reusing `TRT`.

### DOSE_70MG
- **Description:** 1 = subject is on the 70 mg SC Q4W dose regimen, 0 = subject is on the 210 or 490 mg SC Q4W regimen.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (210 mg or 490 mg Q4W regimen).
- **Source aliases:** derived per subject from the trial-assigned dose level.
- **Example models:** `Kotani_2022_astegolimab.R`.
- **Notes:** Zenyatta-study categorical covariate flagging the 70 mg group (lowest dose), modeled as a -15.3% relative change on relative bioavailability. Modeled by Kotani 2022 as `70 mg vs {210 mg, 490 mg}` combined reference.

### DOSE_50MG
- **Description:** 1 = dose record is a 50 mg SC administration, 0 = all other SC doses (100, 150, 200, 300 mg) and all IV doses.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (100-300 mg SC or any IV dose).
- **Source aliases:** derived per dose record from the administered amount (`AMT`).
- **Example models:** `Othman_2014_daclizumab.R`, `Diao_2016_daclizumab_cd25.R`, `Diao_2016_daclizumab_cd56bright.R`, `Diao_2016_daclizumab_treg.R`.
- **Notes:** Othman 2014 estimated two separate absolute bioavailabilities because of non-linear dose-normalized exposure at the 50 mg SC dose -- F = 0.84 for the therapeutic 100-300 mg SC range and F = 0.57 for the 50 mg SC cohort. Encoded here as a record-level indicator so the covariate effect `e_dose_50mg_f = 0.57/0.84 - 1 = -0.321` scales bioavailability only on 50 mg SC doses. For clinical-range simulation (150 mg SC Q4W Phase III regimen) leave `DOSE_50MG = 0`. The Diao 2016 PK/PD models inherit the Othman 2014 PK backbone verbatim; they carry `DOSE_50MG` even though the Diao 2016 RRMS regimens are 150 / 300 mg SC only.

### STUDY1
- **Description:** 1 = subject enrolled in Study 1 of the Cirincione 2017 pooled analysis, 0 = other. Used to switch the residual-error magnitude per study.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (all other studies; combined with `STUDY5 = 0` selects the pooled "other" residual error).
- **Source aliases:**
  - `DVID = "study1"` (character-valued study identifier; `STUDY1 = as.integer(DVID == "study1")`) -- legacy form previously used in `Cirincione_2017_exenatide.R`.
- **Example models:** `Cirincione_2017_exenatide.R`.
- **Notes:** Paired with `STUDY5`. When both are 0, the subject is in the pooled "other studies" residual-error group.

### STUDY5
- **Description:** 1 = subject enrolled in Study 5 of the Cirincione 2017 pooled analysis, 0 = other. Used to switch the residual-error magnitude per study.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (all other studies).
- **Source aliases:**
  - `DVID = "study5"` (character-valued study identifier; `STUDY5 = as.integer(DVID == "study5")`) -- legacy form previously used in `Cirincione_2017_exenatide.R`.
- **Example models:** `Cirincione_2017_exenatide.R`.
- **Notes:** Paired with `STUDY1`. When both are 0, the subject is in the pooled "other studies" residual-error group.

### PHASE2
- **Description:** 1 = subject enrolled in the Phase II study (MORAb-003-002) of the Farrell 2012 pooled analysis; 0 = Phase I study (MORAb-003-001). Used to switch the residual-error magnitude per study.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (Phase I).
- **Source aliases:** derived per subject from the trial identifier (`MORAb-003-001` -> 0, `MORAb-003-002` -> 1).
- **Example models:** `Farrell_2012_farletuzumab.R`.
- **Notes:** Farrell 2012 Table 3 reports separate residual-error estimates for the two studies -- Phase I uses a proportional-only model (sigma = 20.5%); Phase II uses a combined additive + proportional model (sigma_prop = 34.9%, sigma_add = 7.94 ug/mL). The `PHASE2` indicator selects between them.

### PHASE1
- **Description:** 1 = subject enrolled in a Phase 1 study of the Valenzuela 2025 pooled analysis (MOM-M281-001, MOM-M281-007, MOM-M281-010, EDI1001, EDI1002 -- healthy participants); 0 = Phase 2 study (MOM-M281-004 / Vivacity-MG -- participants with gMG). Used to switch the proportional PK residual-error magnitude per study phase.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (Phase 2).
- **Source aliases:** derived per subject from the trial identifier (Phase 1 protocols -> 1, `NCT03772587` Vivacity-MG -> 0).
- **Example models:** `Valenzuela_2025_nipocalimab.R`.
- **Notes:** Valenzuela 2025 Table 3 reports proportional PK residual 0.0834 (Phase 1) vs 0.367 (Phase 2). Distinct from Farrell 2012 `PHASE2` -- the reference category is inverted (Valenzuela 2025 picks Phase 1 as the 1-level).

### STDY_LBSL (**canonical for early-phase belimumab LBSL01 / LBSL02 study indicator**)
- **Description:** 1 = subject enrolled in study LBSL01 (NCT00657007) or LBSL02 (NCT00071487) -- the two early-phase belimumab studies that used a different ELISA-based bioanalytical assay; 0 = any other belimumab study in the Zhou 2021 pooled analysis. Used to switch CL and V1 magnitudes per study group (effectively an assay / early-development PK adjustment).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (later-phase studies using the electrochemiluminescence assay).
- **Source aliases:**
  - `INDR` -- used in `Zhou_2021_belimumab.R` (Zhou 2021 Table 2 footnote: study indicator).
- **Example models:** `Zhou_2021_belimumab.R` (multiplicative factors 1.63 on CL and 1.26 on V1 when STDY_LBSL = 1).
- **Notes:** Conceptually similar to `STUDY1` / `PHASE2` / `ELISA` / `PHASE1` (per-study switches) but specific to the belimumab program. Subject-level (time-fixed); set from the trial identifier on each subject record.

### ELISA
- **Description:** 1 = serum nipocalimab concentration measured by ELISA assay (LLOQ 0.150 ug/mL; studies MOM-M281-001, MOM-M281-007, MOM-M281-010); 0 = measured by ECLIA assay (LLOQ 0.010 ug/mL; studies EDI1001, EDI1002, MOM-M281-004). Used to switch the additive PK residual-error magnitude per assay.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (ECLIA).
- **Source aliases:** derived per study from the bioanalytical method (Text S2 of Valenzuela 2025).
- **Example models:** `Valenzuela_2025_nipocalimab.R`.
- **Notes:** Valenzuela 2025 Table 3 reports additive PK residual 0.445 nmol/L (ELISA) vs 0.0342 nmol/L (ECLIA). Assay choice is study-fixed; include as a per-observation indicator so the correct assay residual is applied even when observations from different studies are pooled.

### STUDY_M281_004
- **Description:** 1 = subject enrolled in study MOM-M281-004 (Vivacity-MG; NCT03772587; phase 2 in generalized myasthenia gravis); 0 = any other study in the Valenzuela 2025 pooled analysis. Used to switch the IgG-baseline scaling factor.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-Vivacity-MG studies).
- **Source aliases:** derived per subject from the trial identifier (`NCT03772587` -> 1, else -> 0).
- **Example models:** `Valenzuela_2025_nipocalimab.R`.
- **Notes:** Valenzuela 2025 estimated a slightly lower IgG baseline in MOM-M281-004 participants (baseline scaled by `FRIgG0_M281_004 = 0.777` vs. 1 in other studies). Distinct from the disease-state indicator implied by `gMG` -- it is specifically the Vivacity-MG study flag because the IgG baseline factor was only estimated for that study.

### STUDY_ABA2_HLA78
- **Description:** 1 = subject enrolled in the ABA2 hematopoietic-cell-transplant trial (IM101311; NCT01743131) HLA 7/8 (one-allele-mismatched donor) cohort, 0 = any other study in the Takahashi 2023 pooled abatacept population PK analysis (RA/JIA reference and ABA2 8/8 cohort).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-ABA2-7/8 -- pooled adult RA / pediatric JIA cohort, plus ABA2 HLA 8/8 cohort which is flagged separately by `STUDY_ABA2_HLA88`).
- **Source aliases:** derived per subject from the trial-cohort identifier (`Cohort = ABA2 7/8` -> 1, else -> 0). Takahashi 2023 Supplemental Table 4 names the corresponding theta `thetaCohort_CL` / `thetaCohort_Vc`.
- **Example models:** `Takahashi_2023_abatacept.R` (multiplicative `Ratio` factors on CL = 0.70 and on Vc = 0.99 vs the RA/JIA reference; values from Takahashi 2023 Supplemental Table 4).
- **Notes:** Pairs with `STUDY_ABA2_HLA88` to reproduce the three-level cohort categorical (RA/JIA, ABA2 7/8, ABA2 8/8) the paper reports as the only retained categorical PK covariate. At most one of `STUDY_ABA2_HLA78` and `STUDY_ABA2_HLA88` is 1 per subject; both 0 reproduces the RA/JIA reference. Scope: specific because the contrast is tied to the ABA2-vs-RA/JIA pooling design.

### STUDY_ABA2_HLA88
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

## Study-site region

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

### ooc1, ooc2, ooc3, ooc4
- **Description:** Mutually exclusive occasion indicators for a crossover / multi-period design. Exactly one is 1 per observation.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Example models:** `Xie_2019_agomelatine.R`.
- **Notes:** Lower case preserved from source file. Pre-existing legacy form; new models should prefer the integer-valued `OCC` canonical above and decompose into binary indicators inside `model()`.

### CYCLE
- **Description:** Treatment cycle number (1 = first dosing cycle, 2 = second, ...). Integer count, time-varying across a multi-cycle treatment course, incremented at each new dosing cycle.
- **Units:** (count)
- **Type:** count
- **Scope:** specific
- **Reference category:** n/a -- used either with a power-covariate form `CYCLE^Fm` (Fm typically negative) to capture cycle-over-cycle decline in a derived quantity such as ADC-to-payload conversion fraction (Li 2017 brentuximab vedotin), or with a piecewise indicator `CYCLE == 1 vs CYCLE >= 2` to capture a step change in DAR scaling between cycle 1 and later cycles (Hong 2025 datopotamab deruxtecan).
- **Source aliases:** `CYCLE` -- used in `Li_2017_brentuximab.R`, `Hong_2025_datopotamab.R`, and `Lu_2022_patritumab.R` with the same canonical name.
- **Example models:**
  - `Li_2017_brentuximab.R` (exponent on the fraction of ADC that converts to MMAE by proteolytic degradation, Fm = -0.261, to reflect tumor-burden reduction across successive treatment cycles).
  - `Hong_2025_datopotamab.R` (cycle-1 vs cycle-2+ piecewise scaling Factor1 = 0.696 on the DAR equation that drives DXd formation rate from total Dato-DXd elimination).
  - `Lu_2022_patritumab.R` (cycle-1 vs cycle-2+ piecewise scaling Factor1 = theta = 0.648 on the payload-to-intact-drug ratio PIR that scales DXd release rate from intact ADC).
- **Notes:** Must be >= 1 throughout (`CYCLE^Fm` is undefined at 0; the piecewise form requires `CYCLE` to be a positive integer at every observation row). Distinct from `ooc<n>` binary-occasion indicators: `CYCLE` is an integer count, not a mutually-exclusive set of indicator columns. Data-assembly helper: set `CYCLE = floor((TIME - TIME_FIRST_DOSE) / cycle_length_days) + 1` for a fixed-interval dosing regimen.

---

## Change log

- **2026-05-06** -- Added `DLVL` (specific scope under `Formulation / assay / study`; integer-coded source-protocol dose-level / regimen indicator) and `S1A2` (specific scope under `Formulation / assay / study`; integer-coded source-protocol CYP1A2 substrate / co-medication indicator) canonical entries while extracting `NA_NA_lidocaine.R` (DDMODEL00000281). Both names match the source-data column names used in the BAST.dat ("4-cRUN249") lidocaine study and survive `=DROP` in the `.ctl` `$INPUT` declaration. The `DLVL` binary derivation `DLVL > 2` switches typical-value baselines for the GX rate constant K30 and the lidocaine apparent central volume V1; the `S1A2` binary derivation `S1A2 == 3` selects a CYP1A2-induction-style modifier on K30. Sibling columns `D1A2` and `H1A2` are dropped in the `.ctl` so only the level-3 indicator is structurally identifiable. Source aliases mapped: `BIL` (legacy total-bilirubin label) -> `TBILI`; `SGPT` (legacy ALT label, paralleling `SGOT` -> `AST`) -> `ALT`. The linked publication is not on disk for this extraction; the absence of an external parameter cross-check is documented in the model file's vignette Errata.
- **2026-05-06** -- Added `FDC` (specific scope under `Formulation / assay / study`; binary fixed-dose-combination antitubercular formulation indicator with FDC = 1 as the typical-value reference) canonical entry while extracting `Wilkins_2008_rifampicin.R` (DDMODEL00000280). The 1-level convention follows the Wilkins 2008 source (FDC = 1 = co-formulated tablet, 0 = single-drug tablets) and the bundle's `Most common` annotation in the `.mod` `IF(FDC.EQ.1) ... = 0` block. Distinct from `TABLET` (Kyhl 2016 tablet-vs-solution) and the `FORM_*` family because the FDC-vs-SDC contrast is a co-formulation contrast, not a drug-product-version one.
- **2026-05-06** -- Added `ORG_FAIL_COUNT` (specific scope under a new `Critical-illness severity` H2 section; integer count of failing organs in critically ill patients, decomposed into per-stratum binary indicators inside `model()` to select per-stratum typical clearance) canonical entry while extracting `Vet_2016_midazolam.R` (DDMODEL00000249). Source alias `ORGF` -> canonical `ORG_FAIL_COUNT` (renamed for descriptive clarity over the source NMTRAN abbreviation). Reference category 0 (no organs failing); per-stratum strata 1 / 2 / 3 / >=4. The new section sits between `Disease severity scores` and `Interferon / biomarker panels`. The Vet 2016 paper PDF is not on disk in this worktree; the absence of paper-on-disk parameter cross-checks is documented in the model file's vignette Errata.
- **2026-05-06** -- Added `TERM_BIRTH` (specific scope under a new `Pregnancy / hormonal status` H2 section; binary term-vs-preterm birth indicator), `BC_USE` (specific scope under the same new `Pregnancy / hormonal status` H2 section; binary oral-contraceptive use indicator), and `UF` (specific scope under `Renal / hepatic function`; continuous instantaneous urine flow rate, mL/h, with centering reference 100 mL/h and a `UF == 0` sentinel-zero rule that gates the linear effect on renal CL) canonical entries while extracting `Allegaert_2015_paracetamol.R` (DDMODEL00000267). Source aliases mapped: `TERM` -> `TERM_BIRTH` (same orientation, no transformation), `BC` -> `BC_USE` (same orientation, no transformation), `UF` -> `UF` (no rename). The Allegaert 2015 BMC Anesthesiol PDF is not on disk in this worktree; covariate semantics, units, and reference values are taken from the bundle's `Executable_OriginalModelCode.mod` `$INPUT` comments and `$PK` block.
- **2026-05-06** -- Added `HIV_POS` (general scope under a new `Infectious disease (HIV)` H2 section; binary HIV-positive comorbidity indicator, parallels the `_POS` suffix convention used by `ADA_POS` / `SARS_SEROPOS`) and `OCC` (general scope under `Occasion / period (IOV)`; integer-valued occasion / period column, decomposed into binary indicators inside `model()` for IOV multiplexing) canonical entries while extracting `Jonsson_2011_ethambutol.R` (DDMODEL00000220). Source aliases mapped: `HIV` -> `HIV_POS` (same orientation, 1 = HIV-positive), `OCC` -> `OCC` (no rename). The new `OCC` integer-valued canonical is the recommended form for new IOV-using models per the existing register-note steer (`ooc1..oocN` binary form retained as the legacy canonical).
- **2026-05-06** -- Added `KG`, `KD0`, `KD1`, `IBASE`, and `NWLS` (all specific scope, under `Oncology`) canonical entries while extracting `Zecchin_2016_survival.R` (Zecchin 2016 OS model, DDMODEL00000218). `KG` / `KD0` / `KD1` / `IBASE` are subject-level empirical-Bayes posterior estimates from the upstream Zecchin 2016 SLD model (DDMODEL00000217 / `Zecchin_2016_tumorovarian.R`) carried into the OS model via the dataset; `NWLS` is the time-varying RECIST-style step-function indicator for new (non-target) lesion appearance. Source aliases recorded: `KG`/`KD0`/`KD1`/`IBASE`/`NWLS` (`NWLSCOV` in the bundle's Simulated_OS.csv) -> corresponding canonicals. The internal `/1000` and `/100` numerical scalings on `KG` / `KD0` / `KD1` and the `*1000` scaling on `IBASE` are preserved verbatim from the source `$DES` block to maintain numerical equivalence with the published estimates (Zecchin 2016 BJCP 82(3):717-727; doi:10.1111/bcp.12994). Extended `AUC_CARBO` and `AUC_GEM` example-models lists to include `Zecchin_2016_survival.R` (the OS model integrates the same SLD ODE inline).
- **2026-05-06** -- Added `AAG` (general-scope serum alpha-1 acid glycoprotein concentration; placed under `Inflammation markers` after `CRP`; reference 1.34 g/L from the Kloft 2006 cancer-cohort median; piecewise-linear effect form with breakpoint at the median documented), `PRIOR_ANTICANCER` (general-scope binary indicator for any prior anticancer therapy of any modality including cytotoxic chemotherapy / radiotherapy / surgery / hormonal therapy / targeted therapy / immunotherapy; placed under `Disease / treatment history` before `PRIOR_BIO`; broader than `LINE_1L` which is restricted to systemic-drug therapy lines), and `CP_MGL` (specific-scope instantaneous drug plasma concentration as a time-varying PD driver, mg/L = ug/mL; placed under `Drug exposure metrics` after `CAV`; distinct from `CAV` which is dosing-interval-averaged) canonical entries while extracting `Netterberg_2017_docetaxel.R` (DDMODEL00000224). Source aliases mapped: `SEX` (1 = male, 2 = female; values inverted via `SEXF = as.integer(SEX == 2)`) -> `SEXF`; `PERF` (ordinal ECOG 0/1/2 with binarization at >= 1) -> `ECOG_GE1`; `PC` (0/1) -> `PRIOR_ANTICANCER`; `AAG` -> `AAG`; `CP` -> `CP_MGL`. The model is the Friberg-style myelosuppression PD model with parameter values fixed from Kloft 2006; Netterberg 2017 used the model unchanged for an ANC-prediction-methodology study (no re-fit) and submitted it to the DDMORE Foundation Model Repository. Neither Netterberg 2017 nor Kloft 2006 is on disk in this worktree; the absence of paper-on-disk parameter cross-checks is documented in the model file's vignette Errata.
- **2026-05-06** -- Added `CYP2D6` and `CYP3A4` canonical entries (both general scope; continuous individual-activity scores) under `Pharmacogenetics` while extracting `TerHeine_2014_tamoxifen.R` (DDMODEL00000212). Both names match the source-data column names used in the Ter Heine 2014 dataset; the dextromethorphan-probe-derived model-based individual clearance values (median 1560 ng/L for CYP2D6; 44.7 ng/L for CYP3A4/5) drive a power-law effect on the tamoxifen -> endoxifen formation clearance. Per-entry `notes` carry an explicit TODO to consolidate the variety of CYP-phenotype encodings into a uniform pattern and to register the rest of the canonical drug-metabolizing-CYP set (`CYP1A2`, `CYP2A6`, `CYP2B6`, `CYP2C8`, `CYP2C9`, `CYP2C19`, `CYP2E1`, `CYP3A5`) prospectively rather than ad hoc as each new model lands.
- **2026-05-06** -- Added `RCFB1MAX` (specific scope, under `Oncology`) canonical entry while extracting `Schindler_2016_sunitinib.R` (DDMODEL00000221). The Schindler 2016 PR (#300) was based on a state predating `DOSE`/`CLI` consolidation and originally proposed drug-specific `CL_SUNITINIB` and `DAILY_DOSE` canonicals; on merge those are folded into the existing general `DOSE` and specific `CLI` canonicals, with the model file and vignette renamed to match. Source aliases on Schindler: `CL` -> `CLI` (existing alias), `DOS` -> `DOSE` (existing alias), `RCFB1MAX` -> `RCFB1MAX` (no rename).
- **2026-05-06** -- Added a new `## Count / Markov-feedback PD covariates` section with `TRT_PHASE` (general scope; double-blind active-treatment-phase indicator), `PDV` (specific scope; previous-period observed seizure count, used to express Markov dependence as a per-record covariate when the underlying inference framework cannot represent observation-to-state feedback as a model state), and `NDAYS` (general scope; number of days in the count-record interval) canonical entries while extracting `Schoemaker_2018_levetiracetam.R` (DDMODEL00000239). `Q2` source column renamed to canonical `TRT_PHASE` because the source name `Q2` collides with the canonical PK parameter `q2` (inter-compartmental clearance to peripheral2). `PDV` retained as the canonical name because the abbreviation is the established convention in the count-likelihood Markov-feedback literature. Extended the `CHILD` source-alias list to record `PED` (Schoemaker 2018), and the `CAV` example-models list to record `Schoemaker_2018_levetiracetam.R` (LEV in mg/L, EC50 about 31.5 mg/L).
- **2026-04-28** -- Extended `RACE_WHITE` (general scope) example-models list and source aliases to record `Hu_2014_bapineuzumab.R` (Caucasian-vs-non-Caucasian dichotomy with the Caucasian subgroup as the typical-value reference). The canonical column encoding (1 = White / 0 = non-White) is unchanged; the model implements the 15% non-Caucasian effect on `(1 - RACE_WHITE)`. The change log notes that the typical-value reference category may legitimately differ between papers using `RACE_WHITE`.
- **2026-05-06** -- Added `CL_INDIV`, `VC_INDIV`, and `VP_INDIV` (all specific scope) under `Drug exposure metrics` while extracting `Friberg_2002_paclitaxel.R` (DDMODEL00000186). The three columns carry per-subject empirical-Bayes paclitaxel PK estimates (clearance, central volume, peripheral volume) that drive the myelosuppression PD model in place of an estimated PK layer; the pattern parallels the existing `CAV` exposure-metric convention but for the structural PK parameters themselves rather than a derived exposure summary. Source aliases mapped: `CLI` -> `CL_INDIV`, `V1I` -> `VC_INDIV`, `V2I` -> `VP_INDIV`. Reserved `VP2_INDIV` for future registration when a second peripheral compartment is needed.
- **2026-04-29** -- Added `IL6` (general-scope serum interleukin-6 cytokine biomarker, pg/mL, under `Inflammation markers`), `PAIN` (general-scope patient-reported global pain VAS 0-100, under `Rheumatoid-arthritis disease-activity covariates`), and `RACE_ASIAN_OTH` (specific-scope composite race indicator pooling Asian / American Indian / Other against a White + Black reference, under `Race / ethnicity`) canonical entries while extracting `Frey_2013_tocilizumab.R` (PMID 23436260). Source aliases mapped: `BLIL6` -> `IL6`. Frey 2013 uses log-transformed `(log(IL-6 * 1000) / 9.9)^exp` with reference IL-6 ~= 20 pg/mL; the canonical column carries the raw IL-6 in pg/mL and the log transformation is applied inside `model()`.
- **2026-04-28** -- Added `TUMTP_PCALCL` (specific scope under `Oncology`; primary cutaneous anaplastic large-cell lymphoma indicator following the `TUMTP_<GROUP>` decomposition pattern) and three ADA-status x assay-era indicators (`ADA_POS` [modern-assay arm; general scope], `ADA_POSOLD`, `ADA_MISSING` [both specific scope] under `Immunogenicity`) while extracting `Suri_2018_brentuximab.R`. Initially named `ADA_POSNEW`; renamed to `ADA_POS` on 2026-04-29 for consistency with the existing general `ADA_POS` canonical. The three indicators are mutually exclusive and decompose Suri 2018's four-level ADA-status x assay-era factor with ADA-negative as the reference; the multiplicative additive form `cl *= (1 + theta * ind)` from supplement 1's ATA-status equation is documented per model.
- **2026-04-27** -- Added `STUDY_ABA2_HLA78` and `STUDY_ABA2_HLA88` (both specific scope) canonical entries while extracting `Takahashi_2023_abatacept.R`. The two binary indicators jointly reproduce the three-level RA/JIA-vs-ABA2-7/8-vs-ABA2-8/8 cohort categorical that Takahashi 2023 Supplemental Table 4 retains as the only categorical PK covariate (multiplicative `Ratio` thetas on CL and on Vc; RA/JIA = both indicators 0 = ratio 1 fixed).
- **2026-04-27** -- Added a new `## Study-site region` section with `REGION_JAPAN`, `REGION_EUROPE`, `REGION_ROW` canonical entries (all scope: specific, binary indicators) while extracting `Hong_2025_datopotamab.R`. US is the implicit reference category (all three indicators = 0). The new entries are distinct from the existing `RACE_*` family because they encode trial-site geography rather than subject ancestry; a Japanese-ancestry subject enrolled at a US site has `RACE_JAPANESE = 1` and `REGION_JAPAN = 0`.
- **2026-04-26** -- Added `B2M` (general-scope serum beta-2-microglobulin under `Renal / hepatic function`; reference 3.90 mg/L from the multiple-myeloma cohort median), `MM_NIGG` (specific-scope non-IgG-MM-vs-IgG-MM within-disease immunoglobulin-type indicator under `Oncology`), and `FORM_P2F2` (specific-scope isatuximab phase III / commercial-bound drug-material indicator, placed alongside the existing `FORM_*` entries) canonical entries while extracting `Fau_2020_isatuximab.R`. Source aliases mapped: `Ig_type`->`MM_NIGG`, `Drug_mat`->`FORM_P2F2`. `Fau_2020_isatuximab.R` added to the `WT`, `SEXF`, and `RACE_ASIAN` example-model lists.
- **2026-04-27** -- Promoted `COADMIN_CHEMO` from specific to general scope after the Kuchimanchi 2024 dostarlimab + carboplatin/paclitaxel population-PK analysis ratified the same pooling convention (any chemotherapy backbone collapsed into a single binary indicator) on a second mAb. The two example models use different functional forms -- Zhang 2019 uses `exp(theta * COADMIN_CHEMO)` and Kuchimanchi 2024 uses `(1 + theta * COADMIN_CHEMO)` -- but the canonical column meaning is unchanged; the per-model functional form is documented in the model files' `notes`. Source-alias list updated: `MONOTR` (Kuchimanchi 2024) added with the inverse-value relation `COADMIN_CHEMO = 1 - MONOTR`. Description and example list extended to include `Kuchimanchi_2024_dostarlimab.R`.
- **2026-04-27** -- Added `ACUTE_MED_DAYS` (specific scope under `Disease severity scores`; baseline number of days/month of acute migraine medication use, piecewise-linear with breakpoint at 5 d/mo per medication-overuse-headache convention) and `CAV` (specific scope under a new `Drug exposure metrics` H2 section; average modelled-drug plasma concentration over a dosing interval, used as the exposure covariate in PD exposure-response models that consume empirical-Bayes Cav from a previously-published popPK model) canonical entries while extracting the two Fiedler-Kelly 2020 fremanezumab exposure-response models (`FiedlerKelly_2020_fremanezumab_em.R`, `FiedlerKelly_2020_fremanezumab_cm.R`).
- **2026-04-26** -- Deduplicated four covariate entries that were registered twice during sequential model-PR merges: (1) merged second `IGG` entry (Yang 2021 cemiplimab) into the first (Zhou 2021 belimumab) -- both papers use the same canonical `IGG` column, so the entry now lists both source aliases (`BIGG`, `IGGBL`) and both reference values (14.8 g/L, 9.65 g/L); (2) merged second `BGENE21` entry (Zheng 2016 sifalimumab, misfiled under Inflammation markers) into the first (Narwal 2013 sifalimumab, under Interferon / biomarker panels); (3) merged `STEROID_BL` (Zheng 2016 sifalimumab) into `STEROID` (Narwal 2013 sifalimumab) -- both encode the same baseline/concomitant corticosteroid indicator, `STEROID` is the preferred shorter name, `STEROID_BL` retired; `Zheng_2016_sifalimumab.R` updated accordingly; (4) merged `ECOG_PS_GT0` (Zhang 2019 nivolumab) into `ECOG_GE1` (Bajaj 2017 nivolumab) -- `>= 1` equals `> 0` for integer ECOG scores, `ECOG_GE1` is the preferred name (consistent with a potential future `ECOG_GE2`), `ECOG_PS_GT0` retired; `Zhang_2019_nivolumab.R` updated accordingly.
- **2026-04-26** -- Added `RACE_JAPANESE` (general scope) canonical entry to formalize the entry that was noted in the changelog but lacked a proper H3 section. Source alias `JAPANESE_HV` from `Wang_2017_benralizumab.R` mapped; `Wade_2015_certolizumab.R` uses `RACE_JAPANESE` directly. Added `MM` (active multiple myeloma, specific scope), `MDSAML` (MDS or AML combined indicator, specific scope), and `SPDL1` (soluble PD-L1, specific scope) canonical entries while fixing convention issues in `Ogasawara_2020_durvalumab.R`.
- **2026-04-26** -- Added `DIS_DMD` (specific scope) canonical entry while extracting `Wojciechowski_2022_domagrozumab.R`. Source alias `SPOP`; orientation matches the source (1 = DMD pediatric patient, 0 = healthy adult volunteer reference). The covariate enters as a `(1 + theta * DIS_DMD)` multiplicative shift (additive on the linear scale, not exponentiated) rather than the typical `theta^DIS_DMD` form, matching Eqs. 7-8 of the paper.
- **2026-04-21** -- Added `RACE_HISPANIC` (general) and `CLD_PREM` (general) canonical entries while extracting `Robbie_2012_palivizumab.R`. Extended `ADA_TITER` example_models with the Robbie 2012 category-by-titer-bin usage, and added `Robbie_2012_palivizumab.R` to the `WT`, `PAGE`, `RACE_BLACK`, `RACE_ASIAN`, and `RACE_OTHER` example lists. `HISPANIC`, `CLD`, and `BPD` recorded as source aliases.
- **Initial seed**: Every covariate observed in `inst/modeldb/` as of the audit. Canonical names established: `SEXF`, `ADA_POS`, `RACE_<GROUP>` prefix. Aliases documented but existing model files not modified.
- **2026-04-19** -- Added `CREAT`, `hsCRP`, `ALB` canonical entries after the
  Phase 6 pilot extracted Fasanmade 2009 infliximab and Thakre 2022
  risankizumab. `ALB` had been used informally in two models; now ratified
  with per-model unit documentation (g/dL vs g/L). `CREAT` chosen over
  `CRE`/`SCR`; `hsCRP` preserves lowercase `hs` prefix per the `eGFR`
  precedent. See `tracking/decision_log.md` in the mab_human_consensus
  project for the deliberation.
- **2026-04-19** -- Added `BSA`, `HGB`, `TBILI`, `GAST`, and `COMB_EOX`
  canonical entries in support of Yamada 2025 zolbetuximab. New sections
  created for `Hematology` and `Surgical history / disease state`.
- **2026-04-19** -- Added `BEOS` (baseline blood eosinophil count, cells/uL)
  and `DOSE_70MG` (Zenyatta 70 mg dose-group indicator) canonical entries
  after extracting the Kotani 2022 astegolimab population PK model.
- **2026-04-19** -- Added `ADA_TITER`, `PRIOR_TNF`, `DISEXT_EP`,
  `DISEXT_OTHER` canonical entries while extracting Moein 2022 etrolizumab.
- **2026-04-20** -- Added `CRP` and `BMI` canonical entries for Chua 2025
  mirikizumab.
- **2026-04-20** -- Added `SMOKE` canonical entry from the Ma 2020 sarilumab
  ANC PopPK/PD extraction.
- **2026-04-20** -- Added `EASI` canonical entry for the Eczema Area and
  Severity Index, introduced by `Tiraboschi_2025_amlitelimab.R`.
- **2026-04-20** -- Added `TUMSZ`, `TUMTP_CHL`, `TUMTP_GC` canonical entries
  with the Budha 2023 tislelizumab extraction.
- **2026-04-20** -- Added `ADA_TITRE` canonical entry while extracting the
  Jackson 2022 ixekizumab paediatric psoriasis PopPK model.
- **Xu 2019 sarilumab**: Added canonical entries `ALBR` (albumin / ULN ratio), `CRCL_BSA` (BSA-normalized creatinine clearance), `BLCRP` (baseline C-reactive protein), and `FORM_DP2` (sarilumab drug product 2 indicator). Extended the `ADA_POS` alias list to include the time-varying `ADA` column used in Xu 2019.
- **Ma 2020 sarilumab DAS28-CRP**: Added canonical entries `BLPHYVAS`, `BLHAQ`, and `PRICORT`. Extended the `BLCRP` entry to record Ma 2020 as a second example model.
- **2026-04-21** -- Added `PHASE2` canonical entry (specific scope) for the Farrell 2012 farletuzumab per-study residual-error switch, analogous to `STUDY1`/`STUDY5` in Cirincione 2017.
- **2026-04-20 (scope + mergers)** -- Introduced the `Scope:` field (general vs specific) on every entry; `checkModelConventions()` now warns when a scope-specific covariate appears in a model that is not in its `Example models` list. Merged `hsCRP` + `BLCRP` + the previously separate standard-assay `CRP` into a single general-scope `CRP` canonical (assay and baseline-vs-time-varying nuances now live in per-model `covariateData[[CRP]]$description` / `notes`). Merged `eGFR` + `CRCL_BSA` into a single general-scope `CRCL` canonical (MDRD-estimated vs measured-CrCl-BSA-normalized nuances now live per-model). Merged `ADA_TITRE` + `ADA_TITER` into a single general-scope `ADA_TITER` canonical (British-reciprocal-dilution vs American-linear-titer zero-encoding conventions now live per-model). Renamed `BEOS` -> `EOS` and `GAST` -> `PRIOR_GAST` to follow the `EASI` / `AGE` / `ALB` (baseline status in notes, not column name) and `PRIOR_TNF` / `PRICORT` (prior-treatment / surgical-history) naming patterns. Promoted `TUMSZ` from Budha-specific to general oncology.
- **2026-04-21** -- Added `WBC` canonical entry (total white blood cell
  count) under `Hematology` while extracting `Mould_2007_alemtuzumab.R`.
  Scope: general; reference 10 x 10^9/L used for the Vmax power covariate
  effect in CLL.
- **2026-04-21** -- Added `MGADL` (Myasthenia Gravis Activities of Daily
  Living score; general scope; baseline-only in Valenzuela 2025),
  `ELISA` / `PHASE1` / `STUDY_M281_004` (all specific scope; assay and
  study flags for Valenzuela 2025 nipocalimab) while extracting
  `Valenzuela_2025_nipocalimab.R`.
- **2026-04-21** -- Added `RACE_JAPANESE` canonical entry (scope: general)
  while extracting `Wade_2015_certolizumab.R`, where Japanese is broken
  out as a distinct category from the broader `RACE_ASIAN` indicator
  (both non-overlapping sub-indicators point to different per-population
  reference groups).
- **2026-04-21** -- Added `STEROID` (general-scope; baseline systemic corticosteroid use, distinct from the existing `PRICORT` which captures prior corticosteroid use), `BGENE21` (specific-scope; baseline 21-gene type I interferon signature under a new `Interferon / biomarker panels` section), and `COHDOSE` (specific-scope; randomized dose cohort in mg/kg) canonical entries while extracting Narwal 2013 sifalimumab. Source aliases mapped: `BSTEROID`->`STEROID`, `DOSE`->`COHDOSE`.
- **2026-04-21** -- Added `AST` (general, hepatic enzyme), `HER2_ECD`
  (specific, HER2 shed ECD biomarker), and `TRAST_BL` (specific,
  baseline trastuzumab concentration from prior therapy) canonical
  entries while extracting `Lu_2014_trastuzumabemtansine.R`.
- **2026-04-21** -- Added `DIS_UC` (ulcerative colitis disease-state
  indicator) and `DIS_SASTHMA` (moderate-to-severe asthma disease-state
  indicator) canonical entries under a new `Disease state
  (cross-population indicators)` section while extracting
  `Hua_2015_anrukinzumab.R`. Both scope: specific, decomposed from a
  four-level disease-state categorical (healthy volunteer, mild-to-moderate
  asthma, moderate-to-severe asthma, UC).
- **2026-04-21** -- Added IBD canonical entries for Rosario 2015 vedolizumab
  extraction: `CALPRO` (fecal calprotectin, general), `CDAI` (Crohn's
  Disease Activity Index, general), `PMAYO` (partial Mayo score, general),
  `IBD_CD` (Crohn's-disease-vs-UC indicator, general), and four
  concomitant-medication indicators `CONMED_AZA`, `CONMED_MP`,
  `CONMED_MTX`, `CONMED_AMINO` (all general scope). New H2 sections
  `Inflammatory-bowel-disease disease-activity covariates`,
  `Inflammatory-bowel-disease diagnosis`, and `Concomitant IBD
  medications`.
- **2026-04-21** -- Added `STEROID_BL` (general-scope baseline/concomitant
  steroid-use indicator; distinct from prior-only `PRICORT`), `BGENE21`
  (specific-scope 21-gene type I interferon signature score), and `DOSE`
  (specific-scope per-subject assigned dose level in mg) canonical entries
  while extracting `Zheng_2016_sifalimumab.R`.
- **2026-04-24** -- Added `CONMED_NSAID` (general-scope concomitant NSAID
  indicator; extends the `CONMED_*` IBD pattern) and `SWOL_28JOINT`
  (general-scope 28-joint swollen joint count; RA disease-activity
  section) canonical entries while extracting `Li_2019_abatacept.R`.
  The abatacept-specific SC phase-2 formulation indicator
  (`FORM_ABA_PHASE2`) was kept model-specific per the nlmixr2lib
  global policy that `FORM_*` covariates are not promoted to canonical
  unless they clearly generalize across multiple drugs.
- **2026-04-24** -- Added `BCVA` canonical entry (best-corrected visual acuity
  score in ETDRS letters; scope: specific; ophthalmology-specific baseline
  input to indirect-response BCVA PD models) while extracting
  `Mulyukov_2018_ranibizumab.R`. Source alias `BVA` mapped. Reference value
  55 letters (study-population mean baseline BCVA).
- **2026-04-24** -- Added `ECOG_GE1` (general-scope Eastern Cooperative Oncology Group performance-status indicator; 1 if ECOG >= 1, reference 0) canonical entry under `Oncology` while extracting `Bajaj_2017_nivolumab.R`. Source alias `PS` / `BPS` mapped to `ECOG_GE1`; Bajaj 2017's ECOG derivation (KPS-to-ECOG crosswalk for one study per Oken 1982) documented in the model's `covariateData` notes, not in the register. Added `Bajaj_2017_nivolumab.R` to the `WT`, `SEXF`, `RACE_ASIAN`, and `CRCL` example-model lists.
- **2026-04-24** -- Added `STATIN_MONO` and `EZE` (both scope: specific; new "Concomitant lipid-lowering medication" section) and `PCSK9` (scope: general; new "Hypercholesterolemia biomarkers" section) canonical entries while extracting `Kuchimanchi_2018_evolocumab.R`. `STATIN_MONO` carries Kuchimanchi-specific "statin monotherapy only" semantics; `EZE` captures ezetimibe use (in the Kuchimanchi 2018 dataset overwhelmingly statin+ezetimibe combination therapy).
- **2026-04-24** -- Added `STATIN` (general-scope concomitant statin indicator) under `Concomitant / prior medication` and `FPCSK9` (specific-scope free PCSK9 concentration) under a new `Cardiometabolic / target biomarkers` section while extracting `Martinez_2019_alirocumab.R`.
- **2026-04-24** -- Added `CYCLE` canonical entry (scope: specific; integer count; power-covariate effect `CYCLE^Fm` on the ADC-to-MMAE proteolytic conversion fraction) while extracting `Li_2017_brentuximab.R`. New entry placed under `Occasion / period (IOV)` since it is conceptually an IOV-like index, but distinct from the binary `ooc<n>` pattern.
- **2026-04-27** -- Added `COMBO_DURVA` (specific scope) canonical entry while extracting `Hwang_2022_tremelimumab.R`. Parallels `COMBO_NIVO` (Sanghavi 2020) but for durvalumab co-administration; selects between monotherapy and combination-with-durvalumab values of the time-varying-CL Tmax and lambda parameters. Source alias `COMB` (NM-TRAN $INPUT data item) mapped.
- **2026-04-24** -- Added `LMET` (general-scope binary indicator for baseline liver metastases) and `TUMTP_OTH` (specific-scope residual "other tumor types" indicator, complement of the named `TUMTP_<GROUP>` indicators in the same model) canonical entries while extracting `Quartino_2019_trastuzumab.R`. Extended `TUMTP_GC` with the Quartino 2019 source alias `TTYPE == "AGC"` and `AST` with the legacy clinical-chemistry alias `SGOT`.
- **2026-04-24** -- Added `LDH` (general-scope serum lactate dehydrogenase),
  `TUMTP_SCLC` (specific-scope small-cell-lung-cancer tumor-type
  indicator extending the `TUMTP_*` decomposition pattern), `LINE_1L`
  (specific-scope first-line vs second-line-or-greater therapy
  indicator), `NIVO_1Q3W` and `NIVO_3Q2W` (specific-scope per-regimen
  nivolumab co-administration indicators on baseline CL), and
  `COMBO_NIVO` (specific-scope any-regimen nivolumab combination-therapy
  indicator on the time-varying CL Emax) canonical entries while
  extracting `Sanghavi_2020_ipilimumab.R`.
- **2026-04-26** -- Added `BLBCELL` (specific-scope baseline CD19+ B cell count under `Inflammation markers`), `ROUTE_IV`, `DEVICE_AI`, `STUDY_APLIOS`, and `STUDY_MIRROR` (all specific-scope under `Formulation / assay / study`) canonical entries while extracting `Yu_2022_ofatumumab.R`. Source aliases mapped: `Bcell0`->`BLBCELL`, "Admin route = IV"->`ROUTE_IV`, "Formulation = AI"->`DEVICE_AI`, "Study = APLIOS"->`STUDY_APLIOS`, "Study = MIRROR"->`STUDY_MIRROR`.
- **2026-04-25** -- Added `ECOG_PS_GT0` (general-scope binary indicator for ECOG performance status > 0; Oncology section), `COADMIN_IPI_3Q3W`, `COADMIN_IPI_1Q6W`, `COADMIN_CHEMO`, and `COADMIN_IPI_ANY` (specific-scope coadministration-regimen indicators; Concomitant / prior medication section) canonical entries while extracting `Zhang_2019_nivolumab.R`. Source aliases mapped: `PS`->`ECOG_PS_GT0`, `IPI3Q3W`->`COADMIN_IPI_3Q3W`, `IPI1Q6W`->`COADMIN_IPI_1Q6W`, `CHEMO`->`COADMIN_CHEMO`, `IPICO`->`COADMIN_IPI_ANY`.
- **2026-04-25** -- Added `FFM` (general-scope, fat-free mass, Janmahasatian formula), `IGG` (general-scope, serum immunoglobulin G), `RACE_NEAS` (specific-scope, North East Asian composite race indicator), and `STDY_LBSL` (specific-scope, early-phase belimumab LBSL01/02 study indicator) canonical entries while extracting `Zhou_2021_belimumab.R`. Source aliases mapped: `BALB`->`ALB`, `BIGG`->`IGG`, `RAC4`->`RACE_NEAS`, `INDR`->`STDY_LBSL`.
- **2026-04-25** -- Added `IGG` (general-scope endogenous serum immunoglobulin G concentration; placed under `Renal / hepatic function` near `LDH` since the cemiplimab paper uses it as a baseline lab covariate alongside ALB, ALT, BMI, and WT) canonical entry while extracting `Yang_2021_cemiplimab.R`. Source alias `IGGBL`->`IGG` mapped. Reference value 9.65 g/L.
- **2026-04-25** -- Added `DIS_PJIA` (polyarticular juvenile idiopathic arthritis disease-state indicator; scope: specific) under the existing `Disease state (cross-population indicators)` section while extracting `Gandhi_2021_abatacept.R`, where Gandhi 2021 pools adult RA with pJIA patients and tests pJIA-vs-RA on bioavailability (additive on logit-F: +3.08). Source alias `JIA` mapped. Reused the existing `SWOL_28JOINT` canonical for Gandhi 2021's swollen-joint-count covariate per operator decision: Gandhi 2021's reported reference SJC = 15 is consistent with the 28-joint scale, the same author group used the 28-joint count in Li 2019 (RA-only), and the paper text does not explicitly identify the joint-count scale.
- **2026-04-25** -- Added new top-level section `Pharmacogenomic SNPs` introducing the canonical pattern `SNP_<GENE>_<RSID>` for binary mutant-allele-presence genotype indicators. Three new entries: `SNP_ICAM1_RS1799969`, `SNP_VEGFA_RS1570360`, `SNP_VEGFA_RS699947` (all scope: specific) -- the first pharmacogenomic-genotype covariates in the register, introduced while extracting the three Papachristos 2020 bevacizumab models (PK / binding QSS / PK/PD). Encoding: 1 = at least one mutant allele present; 0 = homozygous wild-type. Source-paper indicator names `cat`, `cat1`, `cat2` (which are positional within each model's covariate equation rather than formal column names) are recorded as aliases.
- **2026-04-25** -- Added `DIAB` (general-scope binary diabetes-mellitus comorbidity indicator) canonical entry under a new `Comorbidities` H2 section while extracting `Chen_2022_guselkumab.R`. Distinct from a primary-disease indicator (`DIS_*`); used in non-diabetes-primary indications where diabetes is tested as a covariate. Source alias `DIAB` mapped.
- **2026-04-26** -- Added `RACE_WHITE` (general; White vs non-White dichotomy), `HEPIMP_MILD` (general; mild hepatic impairment per NCI ODWG criteria), `NLR` (general; neutrophil-to-lymphocyte ratio, hematology), `SARS_VLOAD` (specific; SARS-CoV-2 baseline viral load), `SARS_SEROPOS` (specific; SARS-CoV-2 baseline serostatus), `OXYSUP_LOW` and `OXYSUP_HIGH` (both specific; decomposed from a 4-level supplemental-oxygen categorical) canonical entries while extracting `Lin_2024_casirivimab.R`. New H2 section `Infectious disease (SARS-CoV-2 / COVID-19)` introduced for the SARS_* and OXYSUP_* entries. Source aliases mapped: `RACE`->`RACE_WHITE` (Lin 2024 binary form), `HEPIMP`->`HEPIMP_MILD`, `VIRAL`->`SARS_VLOAD`, `SERPOS`->`SARS_SEROPOS`, `OXYSTAT1`->`OXYSUP_LOW`, `OXYSTAT2`->`OXYSUP_HIGH`.
- **2026-04-27** -- Added `MAYO_E` (general-scope baseline Mayo endoscopic subscore, integer 0-3) canonical entry under `Inflammatory-bowel-disease disease-activity covariates` while extracting `Faelens_2021_infliximab.R`. Source alias `MPRE` mapped (Faelens 2021 NONMEM column for "Mayo endoscopic score pre-induction"; the dataset's `MPRE = -99` missing sentinel is documented as out-of-domain in the per-model `notes`). Distinct from the full Mayo score (0-12) and partial Mayo `PMAYO` (0-9). Also promoted `DISEXT_EP` from scope: specific to scope: general (the binary "extensive colitis vs not" semantics generalize across UC popPK papers regardless of whether the source dataset also distinguishes an "other" disease-extension category) and added Faelens 2021 source alias `EXTCOL` and example-model entry.
- **2026-04-26** -- Added `HEP_IMP` (general-scope binary indicator for NCI ODWG hepatic impairment, mild or worse vs. normal) under `Renal / hepatic function` and `COMBO_RG` (specific-scope binary indicator for anti-CD20 (rituximab or obinutuzumab) combination therapy) under `Oncology` while extracting `Lu_2019_polatuzumab.R`. Source aliases mapped: `BHPTGRPN` (categorical NCI ODWG group with 9999 missing-value sentinel)->`HEP_IMP`; `COMBO` (Lu 2019 categorical 0/1/2)->`COMBO_RG`.
- **2026-04-27** -- Added `SBCMA` (specific-scope, baseline soluble B-cell maturation antigen, in `Cardiometabolic / target biomarkers`; reference 50 ng/mL) and `COMBO_BELAMAF` (specific-scope, any-combination belantamab mafodotin therapy indicator on Imax of the time-varying CL function, in `Concomitant / prior medication`) canonical entries while extracting `Papathanasiou_2025_belantamab.R`. Source aliases mapped: `SBCMABL`->`SBCMA`, `COMBO`->`COMBO_BELAMAF`.
- **2026-04-28** -- Added `CSF1` (specific-scope plasma colony-stimulating factor 1 / M-CSF concentration; placed under `Inflammation markers`), `CPK` (general-scope serum creatine phosphokinase / creatine kinase; placed under `Renal / hepatic function` next to AST/ALT/LDH although mechanistically a muscle-origin enzyme), and `DIS_CANCER` (specific-scope advanced-solid-tumor cohort indicator; placed under `Disease state (cross-population indicators)`) canonical entries while extracting `Yang_2024_axatilimab.R`. Extended `DIS_HV` example_models with Yang 2024 as the second user (cGVHD-reference complement to Nikanjam 2019's non-HV-oncology reference). Source aliases mapped: `BLCSF1` / `BL_CSF1` (Monolix model-parameter name) -> `CSF1`; `BLCPK` -> `CPK`. Reference values: 549 pg/mL (CSF1), 63 U/L (CPK), pooled-cohort medians from Yang 2024 Table S3.
- **2026-04-28** -- Added `DIS_SAD` (specific-scope binary indicator partitioning hypogammaglobulinaemia patients by mechanism: 1 = secondary antibody deficiency, 0 = primary immunodeficiency; placed under `Disease state (cross-population indicators)`) and `IGM` (specific-scope serum immunoglobulin M concentration as a B-cell humoral-capacity proxy; placed under `Renal / hepatic function` next to `IGG`) canonical entries while extracting `Cheng_2026_immunoglobulin.R`. Source aliases mapped: none (the source paper uses bare prose names "type of immunodeficiency" and "IgM level"). Reference values: 0.21 g/L (IGM), Cheng 2026 pooled-cohort median.
- **2026-04-28** -- Added `HDLC` (general-scope serum HDL cholesterol; placed under `Cardiometabolic / target biomarkers`), `TPRO` (general-scope total serum protein; placed under `Renal / hepatic function` next to `ALB`), and `RHEUMATOID_FACTOR` (general-scope serum rheumatoid factor concentration with log-transform-then-power scaling; placed under `Rheumatoid-arthritis disease-activity covariates`) canonical entries while extracting `Frey_2010_tocilizumab.R`. Source aliases mapped: `HDL-C` / `HDL_C` -> `HDLC`; `PROT` / `TP` -> `TPRO`; `LRF` (Frey 2010 final-equation log-RF variable) and `RF` (universal NONMEM/clinical-PK abbreviation, but explicitly rejected on 2026-04-28 as the canonical name in favour of the unambiguous `RHEUMATOID_FACTOR`) -> `RHEUMATOID_FACTOR`. Reference values: 54 mg/dL (HDLC), 74 g/L (TPRO), 110 U/mL == LRF = 4.7 (RHEUMATOID_FACTOR), all from the Frey 2010 pooled-cohort medians.
- **2026-04-28** -- Added `IGE` (general-scope, baseline serum total immunoglobulin E concentration, in `Cardiometabolic / target biomarkers`; reference 482.4 ng/mL with optional `IU/mL` reporting via `1 IU/mL = 2.42 ng/mL` documented per-model) canonical entry while extracting `Hayashi_2007_omalizumab.R`. Source alias mapped: `IgE0`->`IGE`. Distinguished in the H3 notes from the in-model dynamic `X_TE` IgE state used by mechanism-based binding/turnover anti-IgE models.
- **2026-04-28** -- Added `ECOG_GE2` (general-scope ECOG performance-status >= 2 indicator; pairs with `ECOG_GE1` to retain three-level ECOG = 0 / 1 / >=2 ordinal effects in models that test separate >=1 and >=2 thresholds; placed under `Oncology` directly after `ECOG_GE1`), `MCPROT` (specific-scope serum monoclonal protein concentration; tumor-burden marker for plasma-cell-targeting therapies in multiple myeloma; placed under `Oncology` near other tumor-burden markers; reference values 0 g/dL paper-Vmax-reference / 2.0 g/dL paper-figure-reference; time-varying), and `COMBO_LEN_DEX` (specific-scope lenalidomide + dexamethasone combination-therapy indicator on hematologic-malignancy backbone; placed under `Concomitant / prior medication`; distinct from `COMBO_BELAMAF` which pools Ld with bortezomib-dex and pomalidomide-dex) canonical entries while extracting `Ide_2020_elotuzumab.R`. Source aliases mapped: `ECOG101` (with thresholding `IF(.GT.0.5)` for `ECOG_GE1` and `IF(.GT.1.5)` for `ECOG_GE2`)->both `ECOG_GE1` and `ECOG_GE2`; `TMCPROT`->`MCPROT`; `LENDEX`->`COMBO_LEN_DEX`. Extended `ECOG_GE1` example-models list with `Ide_2020_elotuzumab.R`. The Ide 2020 model also extends the example-model lists for `WT`, `AGE`, `SEXF`, `RACE_ASIAN`, `CRCL` (eGFR), `LDH`, `ALB`, `B2M`, `HEPIMP`, and `LINE_1L`.
- **2026-04-28** -- Added `DIS_AD` (specific-scope, Alzheimer's disease patient indicator) under `Disease state (cross-population indicators)` while extracting `PerezRuixo_2025_posdinemab.R`. Source acts on baseline free p217+tau (R0) in CSF only; PK parameters were unaffected by AD status. Reference category is the pooled healthy-volunteer cohort.
- **2026-04-27** -- Renamed `COMBO_LD` -> `COMBO_LEN_DEX` for clarity (lenalidomide+dexamethasone spelled out to avoid ambiguity with other Ld-like abbreviations). `COMBO_LD` recorded as a retired source alias. All references in `Ide_2020_elotuzumab.R` and its validation vignette updated.
- **2026-04-28** -- Added `SMOKE_CURRENT` and `SMOKE_NEVER` (both general-scope, paired binary indicators for a 3-level smoking-status categorical with former smoker as the implicit reference) canonical entries under `Lifestyle / medical history`, alongside the existing 2-level `SMOKE`. Pattern follows the `RACE_<GROUP>` convention for paired indicators. Introduced while extracting `Hwang_2023_monalizumab.R`, where Table 2 reports separate proportional-shift coefficients on V1 for current smoker (+0.0484) and never smoker (-0.141) with former smoker (n=319/507) as the most-common reference category. The existing `SMOKE` entry remains in use for 2-level (current vs non-smoker) encodings; the per-model documentation cross-references the alternative encoding.
- **2026-04-28** -- Added `TUMTP_BC` (general-scope breast-cancer tumor-type indicator) and `HEPIMP_MOD_MISSING` (specific-scope composite moderate-or-data-missing hepatic impairment indicator) canonical entries while extracting `Lu_2022_patritumab.R`. Source aliases mapped: `TUMTP` (BC level)->`TUMTP_BC`; `HEPATIC`->`HEPIMP_MOD_MISSING` (composite of NCI ODWG group 3 = moderate plus n = 6 missing-data patients pooled together by the Lu 2022 covariate analysis). Extended `HEPIMP_MILD` example_models with Lu 2022.
- **2026-04-29** -- Added paired `HSCT_URD_7OF8` and `HSCT_URD_8OF8` (both specific-scope) canonical entries under `Disease state (cross-population indicators)` while extracting `Zhong_2026_abatacept.R`. The two indicators jointly decompose the 3-level Study IM101311 (ABA2) "cohort" categorical (non-HSCT-cohort / 7-of-8 / 8-of-8) into two orthogonal binary indicators in the same style as `DIS_CANCER` + `DIS_HV`. Source aliases mapped: `COHORT7` -> `HSCT_URD_7OF8`, `COHORT8` -> `HSCT_URD_8OF8`. Reused the existing `DIS_PJIA` canonical (source alias `JIA`) and the existing `AST` and `CRCL` canonicals (source aliases `AST` and `cGFR`). Names match the abbreviations used in the Zhong 2026 Figure 1 caption and describe the actual HLA-matching transplant treatment received.
- **2026-05-06** -- Extended canonical `DOSE` to general scope and folded in a second use case (time-varying current administered dose level supplied as a per-record data column, used in PD-only models that derive a per-cycle exposure metric from a posthoc-CL covariate without instantiating a PK ODE) while extracting the Hansson 2013 sunitinib biomarker PD model (DDMODEL00000197, `Hansson_2013a_sunitinib.R`). Added new specific-scope canonical `CLI` (individual posthoc clearance from an upstream popPK fit) under `Drug exposure metrics`. The two columns jointly carry the per-cycle drug-exposure summary (`AUC = DOSE / CLI`) for a PD-only model that consumes posthoc PK from a separately published popPK fit. Source aliases mapped: `DOS` -> `DOSE`, `CL` -> `CLI`. `CLI` is renamed from the source's `CL` column because `cl` is the canonical nlmixr2 PK parameter name. `CLI` is designed to be reusable by the queue siblings `Hansson_2013b_sunitinib` (DDMODEL00000198, TGI) and any future PD-only model that consumes upstream-PK posthoc estimates as data covariates.
- **2026-05-06** -- Added `CMAX_M1` (specific-scope, empirical-Bayes month-1 / cycle-1 Cmax under `Drug exposure metrics` next to `CAV`), `HYPERT` (general-scope hypertension comorbidity / medical-history indicator under `Comorbidities` next to `DIAB`), `REGI_BID` (specific-scope twice-daily dosing-regimen indicator under `Formulation / assay / study`), and `PREV_AE_SCORE` (specific-scope previous-time-step ordinal AE-score Markov-state covariate under `Disease severity scores`) canonical entries while extracting `Girard_2012_pimasertib.R` (DDMODEL00000215; PAGE 21 (2012) Abstr 2458). Source aliases mapped: `CMAXM1` -> `CMAX_M1`, `MHHY` -> `HYPERT`, `BID` -> `REGI_BID`, `PREVSCOR` -> `PREV_AE_SCORE`. The CTCAE 0-3 ocular-AE grading scheme used in Girard 2012 collapses grades 1+2 into a single 1-2 category in the cumulative-logit conditioning; documented in the per-model notes rather than the register so future Markov-AE models can adopt their own grading scheme.
- **2026-05-06** -- Added `AUC_CARBO` and `AUC_GEM` (both specific-scope per-cycle drug-exposure metrics) canonical entries under `Drug exposure metrics` while extracting `Zecchin_2016_tumorovarian.R` (DDMODEL00000217). Source aliases mapped: `CB` -> `AUC_CARBO` (with the DDMORE bundle's simulated dataset using `AUC0`), `G` -> `AUC_GEM` (with the bundle's simulated dataset using `AUC1`). The new entries follow the `CAV` template but document per-cycle AUC semantics for two specific cytotoxics; sibling entries (`AUC_CISPLATIN`, etc.) should be registered alongside, not by overloading these.
- **2026-05-06** -- Added `BAS_SVEGFR3`, `MRT_SVEGFR3`, `EC50_SVEGFR3` (all three specific-scope, under `Drug exposure metrics` because they are upstream-PD posthoc inputs that drive the per-cycle drug-effect summary in the downstream fatigue model) canonical entries while extracting the Hansson 2013c fatigue / adverse-event Markov + proportional-odds model (DDMODEL00000222, `Hansson_2013c_sunitinib.R`). The three sVEGFR-3 columns jointly carry the upstream Hansson 2013a biomarker indirect-response posthoc estimates (initial condition + turnover MRT + drug-effect EC50) for a downstream PD model that consumes the biomarker trajectory as data covariates instead of instantiating its own biomarker ODE. Source aliases mapped: `BAS3` -> `BAS_SVEGFR3`, `MRT3` -> `MRT_SVEGFR3`, `EC53` -> `EC50_SVEGFR3`. Hansson_2013c also added to existing `DOSE` and `CLI` example-models lists.
- **2026-05-06** -- Added `INS` (specific-scope, plasma insulin time-course regressor in pmol/L) and `GLU` (specific-scope, plasma glucose time-course regressor in mmol/L) canonical entries under `Cardiometabolic / target biomarkers` while extracting the DDMORE bundle `DDMODEL00000227` as `Bizzotto_2016_glucose.R`. Both are time-varying regressor inputs (not covariates that modify a parameter); the model declares `linear(INS, GLU)` so rxode2 linearly interpolates them between dataset rows. Source aliases mapped: `iins` -> `INS`, `iglu` -> `GLU`. The bundle's bracketing columns `insn / glun / td / tn` (used in the bundle's hand-rolled piecewise-linear interpolation) are intentionally not registered -- they are not consumed by the nlmixr2 translation, which uses rxode2's native `linear()` declaration instead.
- **2026-05-06** -- Added `CSS_RBV` (specific-scope, individual posthoc ribavirin steady-state trough concentration from an upstream popPK fit; ng/mL) and `K_RBV` (specific-scope, individual posthoc ribavirin approach-to-Css rate constant from the same upstream popPK fit; 1/day) canonical entries under `Drug exposure metrics` while extracting the DDMORE bundle `DDMODEL00000285` as `Laouenan_2015_ribavirin.R`. The two columns jointly carry per-subject empirical-Bayes (modal) estimates from a Laouenan 2015 upstream RBV popPK fit; the downstream PD model uses them inside an analytical `riba = CSS_RBV * (1 - exp(-K_RBV * t))` expression that drives an Imax-style inhibition of hemoglobin synthesis (kin) in a turnover model. Source aliases mapped: `css_mode` -> `CSS_RBV`, `k_mode` -> `K_RBV`. The lumped exponential parameterization replaces a full one- or two-compartment PK ODE; `K_RBV` is *not* a structural elimination-rate constant.
- **2026-05-06** -- Added `NEUT` (specific-scope, baseline absolute neutrophil count in cells/mm^3; placed under `Hematology` next to `NLR`) and `AUC_BAST_FW` (specific-scope, first-week drug AUC in ug*h/L; placed under `Drug exposure metrics` next to `AUC_CARBO`) canonical entries while extracting the BAST PTTE 2017 teaching guiding-document bundle (DDMODEL00000243) as four separate parametric time-to-event models (`NA_NA_tte_gompertz.R` Event 1 exponential / NEUT + AGE; `NA_NA_tte_gompertz_ev2.R` Event 2 Gompertz / AUC; `NA_NA_tte_lognormal.R` Competing Event 1 log-normal / AGE; `NA_NA_tte_loglogistic.R` Competing Event 2 log-logistic / no covariate). Source alias mapped: `AUC` (NM-TRAN $INPUT column name in DDMODEL00000243) -> `AUC_BAST_FW` (renamed in the register so a generic `AUC` column for a different drug in a future model does not collide). The four models share a common 200-subject hypothetical simulated cohort with covariates AGE (24-84 years), NEUT (1030-14,888 cells/mm^3, ref 4133), PRE_TRE (1-6 prior treatments), MAX_LEG (3-12 mm largest-lesion diameter), AUC_BAST_FW (859-7673 ug*h/L, ref 3065.5), and CMAX (100-813 ug/L); only AGE, NEUT, and AUC_BAST_FW survived covariate testing in the four final models. PRE_TRE / MAX_LEG / CMAX appear in the bundle dataset but in no final model and are therefore not registered.
- **2026-05-06** -- Added `WT_BIRTH` (general-scope, time-fixed birth weight in kg under `Pediatric / maturation` after `GA`) canonical entry while extracting `Voller_2017_phenobarbital.R` (DDMODEL00000256). Source alias `BWEIGHT` (Voller 2017) mapped. The conventional `BWT` abbreviation was rejected as the canonical because it is already used across five existing models (Gandhi 2021, Li 2019, Chen 2022, Wojciechowski 2022, Lu 2019) as a source-name alias for body weight `WT`; reusing `BWT` for birth weight would silently break those mappings. The chosen `WT_BIRTH` form preserves the `WT` root and follows the existing `<concept>_<modifier>` register convention. Reference value 2.59 kg observed in Voller 2017's preterm/term-newborn cohort. Operator-confirmed naming choice via runner sidecar (response-001, 2026-05-06).
- Subsequent additions: append new canonical entries as new papers are processed. When adding, bump the audit-completed count in the summary below.

## Summary

- Files audited: 69 R files under `inst/modeldb/` (20 of which reference covariates).
- Canonical H3 entries: 59 (62 parsed entries when counting each `ooc<n>` individually; was ~57 before the 2026-04-20 mergers of `hsCRP`+`BLCRP`+standard-`CRP` -> `CRP`, `eGFR`+`CRCL_BSA` -> `CRCL`, and `ADA_TITRE`+`ADA_TITER` -> `ADA_TITER`; +1 `PHASE2` and +1 `WBC` on 2026-04-21 from Farrell 2012 / Mould 2007; +3 on 2026-04-21 for `STEROID`, `BGENE21`, `COHDOSE` from Narwal 2013; +1 `WT_BIRTH` on 2026-05-06 from Voller 2017).
- Scope: general: 36. Scope: specific: 26 (counting each `ooc<n>` individually, or 23 counting the `### ooc1, ooc2, ooc3, ooc4` heading as one entry).
- Aliases mapped (selected): SEXM->SEXF, ADA->ADA_POS, ADA_TITRE/ADAT->ADA_TITER, BLACK->RACE_BLACK, ASIAN->RACE_ASIAN, MULTIRACIAL->RACE_MULTI, BLACK_OTH->RACE_BLACK_OTH, ASIAN_AMIND_MULTI->RACE_ASIAN_AMIND_MULTI, DVID->STUDY1/STUDY5, CRE->CREAT, hsCRP/HSCRP/CRPHS/BLCRP->CRP, eGFR/EGFR/CRCL_BSA/1.73*CrCl/BSA->CRCL, DP2->FORM_DP2, DISEXT->DISEXT_EP/DISEXT_OTHER, BEASI->EASI, BEOS->EOS, GAST->PRIOR_GAST, COMB->COMB_EOX, TUMTP->TUMTP_CHL/TUMTP_GC, BSTEROID->STEROID, DOSE->COHDOSE.
- Canonical H3 entries: 53 (56 parsed entries when counting each `ooc<n>` individually; was ~57 before the 2026-04-20 mergers of `hsCRP`+`BLCRP`+standard-`CRP` -> `CRP`, `eGFR`+`CRCL_BSA` -> `CRCL`, and `ADA_TITRE`+`ADA_TITER` -> `ADA_TITER`).
- Scope: general: 33. Scope: specific: 23 (counting each `ooc<n>` individually, or 20 counting the `### ooc1, ooc2, ooc3, ooc4` heading as one entry).
- Aliases mapped (selected): SEXM->SEXF, ADA->ADA_POS, ADA_TITRE/ADAT->ADA_TITER, BLACK->RACE_BLACK, ASIAN->RACE_ASIAN, MULTIRACIAL->RACE_MULTI, BLACK_OTH->RACE_BLACK_OTH, ASIAN_AMIND_MULTI->RACE_ASIAN_AMIND_MULTI, DVID->STUDY1/STUDY5, CRE->CREAT, hsCRP/HSCRP/CRPHS/BLCRP->CRP, eGFR/EGFR/CRCL_BSA/1.73*CrCl/BSA->CRCL, DP2->FORM_DP2, DISEXT->DISEXT_EP/DISEXT_OTHER, BEASI->EASI, BEOS->EOS, GAST->PRIOR_GAST, COMB->COMB_EOX, TUMTP->TUMTP_CHL/TUMTP_GC, DX->IBD_CD, AZA->CONMED_AZA, MP->CONMED_MP, MTX->CONMED_MTX, AMINO->CONMED_AMINO.
