# Canonical covariate columns

This file is the authoritative register of covariate column names used in nlmixr2lib models. Every covariate referenced inside a model's `model()` block must have an entry here. The register is seeded from a full audit of `inst/modeldb/` and extended whenever a new paper introduces a covariate that isn't yet registered.

## How to use this register

1. **Before adding a covariate to a new model**, search this file (by canonical name and by source alias) for the concept you need.
2. **If the canonical name exists**, use it exactly. Document any source-paper rename in the model's `covariateData[[name]]$source_name` field.
3. **If the source paper uses an alias listed under an existing canonical name**, use the canonical name and note in `covariateData[[name]]$notes` whether the values must be transformed (e.g., `SEXM → SEXF` inverts values; the effect coefficient sign / reference category must be inverted as well).
4. **If the covariate is not in this register at all**, propose a new entry with a canonical name, description, units, type, and source aliases. Verify with the user before committing. The addition is part of the model's PR.
5. **Do not modify existing model files when you discover an alias**; simply document the mapping in the register. Retrofitting existing models is a separate effort.

## Case convention

Covariate column names should be ALL CAPS unless the source paper uses a specific case that is well-known and unambiguous (e.g., `eGFR` is widely recognized in that case). Current registered non-all-caps names are noted below; new entries should default to all caps.

## Entry schema

```yaml
- name: <CANONICAL_NAME>
  description: <one-sentence definition>
  units: <unit string, or "(binary)" / "(categorical)">
  type: continuous | binary | categorical | count
  reference_category: <the 0 group for binary/categorical, or NULL>
  source_aliases:
    - <ALIAS_NAME> (<transformation if any>) — used in <model.R>
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
- **Reference category:** n/a — used with allometric scaling `(WT / ref_wt)^exponent`. Reference weights observed: 70 kg (adults), 75 kg, 84.8 kg, 5 kg (infants).
- **Source aliases:** none known.
- **Example models:** `Clegg_2024_nirsevimab.R`, `Hu_2026_clesrovimab.R`, `Zhu_2017_lebrikizumab.R`, `Kovalenko_2020_dupilumab.R`, `CarlssonPetri_2021_liraglutide.R`, `Cirincione_2017_exenatide.R`, `Grimm_2023_gantenerumab.R`, `Grimm_2023_trontinemab.R`, `Kyhl_2016_nalmefene.R`, `Soehoel_2022_tralokinumab.R`, `Xie_2019_agomelatine.R`, `PK_2cmt_mAb_Davda_2014.R`, `phenylalanine_charbonneau_2021.R`.
- **Notes:** Universal. Verify time-varying vs. baseline-only against the source paper.

### AGE
- **Description:** Subject age in years.
- **Units:** years
- **Type:** continuous
- **Reference category:** n/a
- **Source aliases:** none.
- **Example models:** `Kyhl_2016_nalmefene.R`, `Zhu_2017_lebrikizumab.R`.
- **Notes:** Zhu 2017 normalizes as `AGE/40`.

### LBM
- **Description:** Lean body mass.
- **Units:** kg
- **Type:** continuous
- **Reference category:** n/a
- **Source aliases:** none.
- **Example models:** `Kyhl_2016_nalmefene.R` (reference 56.28 kg, exponent 0.626 on CL).

### SEXF (**canonical for sex**)
- **Description:** Biological sex indicator, 1 = female, 0 = male.
- **Units:** (binary)
- **Type:** binary
- **Reference category:** 0 (male).
- **Source aliases:**
  - `SEXM` (values inverted: `SEXF = 1 - SEXM`; effect coefficient sign and reference category both invert) — used in `CarlssonPetri_2021_liraglutide.R`.
  - `SEX` with `"M"`/`"F"` strings — derive `SEXF = as.integer(SEX == "F")`.
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical), `CarlssonPetri_2021_liraglutide.R` (alias `SEXM`).
- **Notes:** When translating a model that used `SEXM`, flag the sign/reference-category inversion to the user.

### CHILD
- **Description:** 1 = subject is a child, 0 = not a child.
- **Units:** (binary)
- **Type:** binary
- **Reference category:** 0 (not child, i.e., adult baseline).
- **Example models:** `CarlssonPetri_2021_liraglutide.R`.
- **Notes:** Age-group indicator used alongside `ADOLESCENT`; paper's age cutoffs must be captured in `covariateData[[CHILD]]$notes`.

### ADOLESCENT
- **Description:** 1 = subject is an adolescent, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Reference category:** 0.
- **Example models:** `CarlssonPetri_2021_liraglutide.R`.
- **Notes:** Paired with `CHILD`. Document age cutoffs.

## Pediatric / maturation

### PAGE
- **Description:** Postmenstrual age in months (`GA_weeks / 4.35 + postnatal_months`). Time-varying.
- **Units:** months
- **Type:** continuous
- **Example models:** `Clegg_2024_nirsevimab.R`.

### PNA
- **Description:** Postnatal age (chronological since birth). Time-varying.
- **Units:** months
- **Type:** continuous
- **Example models:** `Hu_2026_clesrovimab.R`.

### GA
- **Description:** Gestational age at birth. Time-fixed per subject.
- **Units:** weeks
- **Type:** continuous
- **Example models:** `Hu_2026_clesrovimab.R`, used implicitly in `Clegg_2024_nirsevimab.R` (folded into PAGE).

## Renal / hepatic function

### eGFR
- **Description:** Estimated glomerular filtration rate (MDRD equation).
- **Units:** mL/min/1.73 m²
- **Type:** continuous
- **Example models:** `Cirincione_2017_exenatide.R` (reference 80 mL/min/1.73 m², exponent 0.838).
- **Notes:** Case preserved (lowercase `e`) because `eGFR` is the widely-used clinical notation.

### CREAT (**canonical for serum creatinine**)
- **Description:** Serum creatinine concentration (baseline or time-varying).
- **Units:** µmol/L or mg/dL — document the unit used in each model via `covariateData[[CREAT]]$units`.
- **Type:** continuous
- **Reference category:** n/a — used with power scaling `(CREAT / ref)^exponent`.
- **Source aliases:**
  - `CRE` (µmol/L, reference 70.73) — used in `Thakre_2022_risankizumab.R`.
  - `SCR` — common clinical-PK abbreviation.
- **Example models:** `Thakre_2022_risankizumab.R`.
- **Notes:** `CREAT` chosen over the shorter `CRE`/`SCR` as the NONMEM/clinical-PK convention that is unambiguous. Per-model reference values must be documented in `covariateData[[CREAT]]$notes`.

### ALB (**canonical for serum albumin**)
- **Description:** Serum albumin concentration.
- **Units:** g/dL or g/L — document the unit used in each model via `covariateData[[ALB]]$units`.
- **Type:** continuous
- **Reference category:** n/a — used with power scaling `(ALB / ref)^exponent`.
- **Source aliases:** none; `ALB` is the universal abbreviation.
- **Example models:** `Fasanmade_2009_infliximab.R` (g/dL, reference 4.1), `Thakre_2022_risankizumab.R` (g/L, reference 45).
- **Notes:** Ratified canonically on 2026-04-19 after cross-model review. Unit varies by paper (g/dL in US-convention papers, g/L in SI-convention papers); the per-model `covariateData[[ALB]]$units` field is load-bearing. Effect-coefficient magnitude is meaningless without the unit.

## Inflammation markers

### hsCRP (**canonical for high-sensitivity C-reactive protein**)
- **Description:** High-sensitivity C-reactive protein concentration (baseline or time-varying).
- **Units:** mg/L (document per-model via `covariateData[[hsCRP]]$units`).
- **Type:** continuous
- **Reference category:** n/a — used with power scaling `(hsCRP / ref)^exponent`.
- **Source aliases:**
  - `CRPHS` (mg/L, reference 5.21) — used in `Thakre_2022_risankizumab.R`.
  - `HSCRP` (all caps).
- **Example models:** `Thakre_2022_risankizumab.R`.
- **Notes:** Case preserved (`hsCRP`) because that is the widely recognized clinical notation, mirroring the `eGFR` precedent. Distinct from non-hs CRP — only use for assays validated for the low-range sensitivity.

## Race / ethnicity

**Canonical pattern: `RACE_<GROUP>`.** Use one indicator per race/ethnicity group the source models. Reference category is the implicit 0 = all other groups; document explicitly which groups are in the reference. When the source uses composite groups (e.g., "Black or Other"), name them accordingly (`RACE_BLACK_OTHER`) and list the components in `notes`.

### RACE_BLACK (**canonical**)
- **Description:** 1 = Black / African American, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Reference category:** 0 (document the actual reference groups used).
- **Source aliases:**
  - `BLACK` — used in `Hu_2026_clesrovimab.R`.
  - `BLACK_OTH` (Black or Other composite; reference: White or Native Hawaiian/Pacific Islander) — used in `Clegg_2024_nirsevimab.R`. Register as its own canonical `RACE_BLACK_OTH` since the grouping is not equivalent.
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical form).

### RACE_BLACK_OTH (**canonical for composite Black/Other group**)
- **Description:** 1 = Black/African American or Other race, 0 = other groups.
- **Reference category:** 0 = White or Native Hawaiian/Pacific Islander (Clegg 2024 grouping).
- **Source aliases:** `BLACK_OTH` — used in `Clegg_2024_nirsevimab.R`.
- **Notes:** Kept distinct from `RACE_BLACK` because the composite is not interchangeable.

### RACE_ASIAN (**canonical**)
- **Description:** 1 = Asian, 0 = other.
- **Source aliases:** `ASIAN` — used in `Hu_2026_clesrovimab.R`.
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical form).

### RACE_ASIAN_AMIND_MULTI (**canonical for composite group**)
- **Description:** 1 = Asian, American Indian / Alaskan Native, or Multiple races, 0 = other.
- **Reference category:** 0 = White, Black / African American, Native Hawaiian / Pacific Islander, or Other (Clegg 2024 grouping).
- **Source aliases:** `ASIAN_AMIND_MULTI` — used in `Clegg_2024_nirsevimab.R`.
- **Notes:** Clegg 2024 applies this covariate to both CL and V2 with different coefficients.

### RACE_MULTI
- **Description:** 1 = multiracial, 0 = other.
- **Source aliases:** `MULTIRACIAL` — used in `Hu_2026_clesrovimab.R`.

### RACE_OTHER
- **Description:** 1 = race category "Other," 0 = not.
- **Example models:** `Zhu_2017_lebrikizumab.R`.

## Laboratory / disease-activity

### ALBR
- **Description:** Serum albumin normalized to the laboratory's upper limit of normal (`albumin_observed / ULN_albumin`).
- **Units:** (unitless ratio)
- **Type:** continuous
- **Reference category:** n/a — used as a power term `(ALBR / <ref>)^exponent`. Reference 0.78 used in Xu 2019 (corresponds to a median serum albumin of 38 g/L at a typical ULN of ~48.7 g/L).
- **Source aliases:** none.
- **Example models:** `Xu_2019_sarilumab.R`.
- **Notes:** Xu 2019 normalizes to each site's ULN so that values across multiple labs with different reference ranges can be pooled.

### CRCL_BSA
- **Description:** Body-surface-area-normalized creatinine clearance, computed as `1.73 * CrCl / BSA` where `CrCl` is in mL/min and `BSA` is in m^2.
- **Units:** mL/min/1.73 m^2
- **Type:** continuous
- **Reference category:** n/a — used as a power term `(CRCL_BSA / <ref>)^exponent`. Reference 100 mL/min/1.73 m^2 in Xu 2019.
- **Source aliases:**
  - `1.73*CrCl/BSA` (the formula appearing in Xu 2019 Eq. for Vm) — used in `Xu_2019_sarilumab.R`.
- **Example models:** `Xu_2019_sarilumab.R`.
- **Notes:** Distinct from `eGFR`: CRCL_BSA is derived from a measured creatinine clearance rather than an MDRD / CKD-EPI estimate.

### BLCRP
- **Description:** Baseline (pre-treatment) C-reactive protein concentration.
- **Units:** mg/L
- **Type:** continuous
- **Reference category:** n/a — used as a power term `(BLCRP / <ref>)^exponent`. Reference 14.2 mg/L in Xu 2019; 15.7 mg/L in Ma 2020.
- **Source aliases:** none.
- **Example models:** `Xu_2019_sarilumab.R`, `Ma_2020_sarilumab_das28crp.R`.
- **Notes:** Time-fixed per subject (baseline value only). Xu 2019 reports it as a significant covariate on Vm with a small exponent (0.0299). Ma 2020 uses it as a power covariate on DAS28-CRP BASE and as an additive log-linear effect on logit(Emax).

## Rheumatoid-arthritis disease-activity covariates

### BLPHYVAS
- **Description:** Baseline Physician's Global Assessment of Disease Activity, 100-mm visual analogue scale (0 = no disease activity, 100 = maximum). Time-fixed per subject.
- **Units:** mm (0-100 VAS)
- **Type:** continuous
- **Reference category:** n/a — used as a power term `(BLPHYVAS / <ref>)^exponent`. Reference 66 used in Ma 2020.
- **Source aliases:** none.
- **Example models:** `Ma_2020_sarilumab_das28crp.R`.
- **Notes:** One of the components of the DAS28 composite score; in Ma 2020 it appears as a baseline covariate on the DAS28-CRP disease-activity BASE rather than on the score itself.

### BLHAQ
- **Description:** Baseline Health Assessment Questionnaire Disability Index (HAQ-DI; 0 = no disability, 3 = maximum disability). Time-fixed per subject.
- **Units:** unitless (0-3 composite score)
- **Type:** continuous
- **Reference category:** n/a — used as a power term `(BLHAQ / <ref>)^exponent`. Reference 1.75 used in Ma 2020.
- **Source aliases:** none.
- **Example models:** `Ma_2020_sarilumab_das28crp.R`.
- **Notes:** Patient-reported disability score frequently used as a baseline covariate in rheumatoid-arthritis PK/PD analyses.

## Concomitant / prior medication

### PRICORT
- **Description:** 1 = patient received systemic corticosteroid treatment prior to study entry, 0 = no prior corticosteroid use. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Reference category:** 0 (no prior corticosteroid use).
- **Source aliases:** none.
- **Example models:** `Ma_2020_sarilumab_das28crp.R` (multiplicative 1.26 on DAS28-CRP Kout).
- **Notes:** Ma 2020 applies it as a multiplicative effect of the form `Kout_i = Kout * theta^PRICORT`.

## Immunogenicity

### ADA_POS (**canonical**)
- **Description:** 1 = antidrug-antibody-positive, 0 = ADA-negative.
- **Units:** (binary)
- **Type:** binary
- **Reference category:** 0 (ADA-negative).
- **Source aliases:**
  - `ADA` (semantically "ever positive") — used in `Zhu_2017_lebrikizumab.R`. When translating from a paper that uses `ADA` as "ever positive," verify the time-frame matches ADA_POS semantics before renaming.
  - `ADA` (time-varying positivity, primary covariate in Xu 2019) — used in `Xu_2019_sarilumab.R`.
- **Example models:** `Clegg_2024_nirsevimab.R`, `Hu_2026_clesrovimab.R`, `Xu_2019_sarilumab.R`.

## Formulation / assay / study

### FED
- **Description:** 1 = fed state at dosing, 0 = fasted.
- **Type:** binary
- **Reference category:** 0 (fasted).
- **Example models:** `Kyhl_2016_nalmefene.R`.

### TABLET
- **Description:** 1 = tablet formulation, 0 = solution.
- **Type:** binary
- **Reference category:** 0 (solution).
- **Example models:** `Kyhl_2016_nalmefene.R`.

### RIA_ASSAY
- **Description:** 1 = radioimmunoassay; 0 = LC-MS/MS.
- **Type:** binary
- **Example models:** `Kyhl_2016_nalmefene.R`.
- **Notes:** Switches the additive residual-error magnitude.

### FORM_NS0
- **Description:** 1 = NS0 cell-line formulation, 0 = other.
- **Type:** binary
- **Example models:** `Zhu_2017_lebrikizumab.R`.

### FORM_CHO_PHASE2
- **Description:** 1 = CHO Phase 2 formulation, 0 = other.
- **Type:** binary
- **Example models:** `Zhu_2017_lebrikizumab.R`.

### FORM_DP2
- **Description:** 1 = sarilumab drug product 2 formulation (used in some phase I studies and the dose-ranging phase II study), 0 = other drug product (DP1 or DP3; DP3 is the commercial formulation).
- **Units:** (binary)
- **Type:** binary
- **Reference category:** 0 (DP1 or DP3).
- **Source aliases:**
  - `DP2` — used in `Xu_2019_sarilumab.R`.
- **Example models:** `Xu_2019_sarilumab.R`.
- **Notes:** Affects both CLO/F (1.30x multiplier) and Ka (0.663x multiplier) in Xu 2019. Set to 0 for routine commercial-formulation simulation.

### dilution
- **Description:** 1 = drug diluted (Soehoel 2022 study D2213C00001), 0 = not diluted.
- **Type:** binary
- **Example models:** `Soehoel_2022_tralokinumab.R`.
- **Notes:** Lower-case preserved from source; future models should rename to `DILUTION`. Kept as alias here to match existing file.

### nonECZTRA
- **Description:** 1 = not the ECZTRA trial; 0 = ECZTRA.
- **Type:** binary
- **Example models:** `Soehoel_2022_tralokinumab.R`.
- **Notes:** Mixed case preserved from source; future models should rename to `NON_ECZTRA` or `STUDY_NON_ECZTRA`.

### SEASON2
- **Description:** 1 = second RSV season at dosing, 0 = first RSV season.
- **Type:** binary
- **Example models:** `Clegg_2024_nirsevimab.R`.
- **Notes:** Study-specific but semantically general (second-exposure indicator).

### STUDY1
- **Description:** 1 = subject enrolled in Study 1 of the Cirincione 2017 pooled analysis, 0 = other. Used to switch the residual-error magnitude per study.
- **Type:** binary
- **Reference category:** 0 (all other studies; combined with `STUDY5 = 0` selects the pooled "other" residual error).
- **Source aliases:**
  - `DVID = "study1"` (character-valued study identifier; `STUDY1 = as.integer(DVID == "study1")`) — legacy form previously used in `Cirincione_2017_exenatide.R`.
- **Example models:** `Cirincione_2017_exenatide.R`.
- **Notes:** Paired with `STUDY5`. When both are 0, the subject is in the pooled "other studies" residual-error group.

### STUDY5
- **Description:** 1 = subject enrolled in Study 5 of the Cirincione 2017 pooled analysis, 0 = other. Used to switch the residual-error magnitude per study.
- **Type:** binary
- **Reference category:** 0 (all other studies).
- **Source aliases:**
  - `DVID = "study5"` (character-valued study identifier; `STUDY5 = as.integer(DVID == "study5")`) — legacy form previously used in `Cirincione_2017_exenatide.R`.
- **Example models:** `Cirincione_2017_exenatide.R`.
- **Notes:** Paired with `STUDY1`. When both are 0, the subject is in the pooled "other studies" residual-error group.

## Occasion / period (IOV)

### ooc1, ooc2, ooc3, ooc4
- **Description:** Mutually exclusive occasion indicators for a crossover / multi-period design. Exactly one is 1 per observation.
- **Type:** binary
- **Example models:** `Xie_2019_agomelatine.R`.
- **Notes:** Lower case preserved from source file. Future models should standardize on `OCC` as a categorical column with integer values (`OCC = 1`, `2`, …) and use `IOV` on the appropriate parameters.

---

## Change log

- **Initial seed** (this PR): Every covariate observed in `inst/modeldb/` as of the audit. Canonical names established: `SEXF`, `ADA_POS`, `RACE_<GROUP>` prefix. Aliases documented but existing model files not modified.
- **2026-04-19** — Added `CREAT`, `hsCRP`, `ALB` canonical entries after the
  Phase 6 pilot extracted Fasanmade 2009 infliximab and Thakre 2022
  risankizumab. `ALB` had been used informally in two models; now ratified
  with per-model unit documentation (g/dL vs g/L). `CREAT` chosen over
  `CRE`/`SCR`; `hsCRP` preserves lowercase `hs` prefix per the `eGFR`
  precedent. See `tracking/decision_log.md` in the mab_human_consensus
  project for the deliberation.
- Subsequent additions: append new canonical entries as new papers are processed. When adding, bump the audit-completed count in the summary below.
- **Xu 2019 sarilumab**: Added canonical entries `ALBR` (albumin / ULN ratio), `CRCL_BSA` (BSA-normalized creatinine clearance), `BLCRP` (baseline C-reactive protein), and `FORM_DP2` (sarilumab drug product 2 indicator). Extended the `ADA_POS` alias list to include the time-varying `ADA` column used in Xu 2019.
- **Ma 2020 sarilumab DAS28-CRP**: Added canonical entries `BLPHYVAS` (baseline Physician's Global Assessment of Disease Activity, 100-mm VAS), `BLHAQ` (baseline HAQ-DI), and `PRICORT` (prior corticosteroid treatment). Extended the `BLCRP` entry to record Ma 2020 as a second example model (reference 15.7 mg/L, covariate on DAS28-CRP BASE and log(Emax)).

## Summary

- Files audited: 61 R files under `inst/modeldb/` (12 of which reference covariates).
- Canonical entries: 33.
- Aliases mapped: 12 (including SEXM→SEXF, ADA→ADA_POS, BLACK→RACE_BLACK, ASIAN→RACE_ASIAN, MULTIRACIAL→RACE_MULTI, BLACK_OTH→RACE_BLACK_OTH, ASIAN_AMIND_MULTI→RACE_ASIAN_AMIND_MULTI, DVID→STUDY1/STUDY5, CRE→CREAT, CRPHS→hsCRP, 1.73*CrCl/BSA→CRCL_BSA, DP2→FORM_DP2).
