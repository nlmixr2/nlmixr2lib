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

### BMI
- **Description:** Body mass index at baseline.
- **Units:** kg/m²
- **Type:** continuous
- **Reference category:** n/a — used with a linear-deviation form (`1 + e * (BMI - ref)`) or a power form (`(BMI / ref)^e`). Document the reference value in `covariateData[[BMI]]$notes`.
- **Source aliases:** none known.
- **Example models:** `Chua_2025_mirikizumab.R` (reference 24.75 kg/m²; linear-deviation effect on logit of bioavailability).
- **Notes:** Universal clinical-trial demographic. Derived as `WT / (height_m)^2`; assume time-fixed at baseline unless the source paper states otherwise.

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

## Disease severity scores

### EASI (**canonical for Eczema Area and Severity Index**)
- **Description:** Eczema Area and Severity Index score (atopic-dermatitis severity composite; bounded continuous, scale 0-72 with higher values = more severe disease).
- **Units:** (score)
- **Type:** continuous
- **Reference category:** n/a — healthy volunteers have EASI = 0. Effect enters as an additive term in models that pool AD patients with HV (e.g., `Tiraboschi_2025_amlitelimab.R`).
- **Source aliases:**
  - `BEASI` (baseline EASI) — used in `Tiraboschi_2025_amlitelimab.R`.
- **Example models:** `Tiraboschi_2025_amlitelimab.R`.
- **Notes:** When used as a time-invariant baseline covariate (`BEASI`), document in `covariateData[[EASI]]$notes`. Canonical name is `EASI` without the `B` prefix to match the `AGE` / `WT` / `ALB` pattern where baseline vs time-varying status is recorded in notes rather than the column name.

## Inflammation markers

### CRP (**canonical for C-reactive protein**)
- **Description:** C-reactive protein concentration (baseline or time-varying) from a standard (not low-range / high-sensitivity) assay.
- **Units:** mg/L (document per-model via `covariateData[[CRP]]$units`).
- **Type:** continuous
- **Reference category:** n/a — document per-model reference value in `covariateData[[CRP]]$notes`.
- **Source aliases:** none known.
- **Example models:** `Chua_2025_mirikizumab.R` (mg/L, reference 7.41); `Moein_2022_etrolizumab.R` (reference 4.23 mg/L, exponential effect on CL).
- **Notes:** Use when the source paper reports "CRP" without specifying a high-sensitivity assay. In IBD and other chronic-inflammation populations, baseline CRP is typically well above the hs-CRP sensitivity range, so a standard assay is adequate. Distinct from `hsCRP`; do not conflate.

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

## Oncology

### TUMSZ (**canonical for baseline tumor size**)
- **Description:** Baseline tumor size. For solid tumors, the sum of diameters of target lesions per RECIST; for classical Hodgkin lymphoma, the sum of products of perpendicular diameters (SPPD).
- **Units:** mm
- **Type:** continuous
- **Reference category:** n/a — used with power scaling `(TUMSZ / ref)^exponent`. Reference values observed: 63 mm (Budha 2023).
- **Source aliases:** none.
- **Example models:** `Budha_2023_tislelizumab.R` (reference 63 mm).
- **Notes:** The SPPD convention for cHL and the sum-of-diameters convention for solid tumors are pooled onto a single column by convention; document the per-model mixture where relevant.

### TUMTP_CHL (**canonical for classical Hodgkin lymphoma tumor-type indicator**)
- **Description:** 1 = classical Hodgkin lymphoma (cHL), 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Reference category:** 0 = all other tumor types (e.g., NSCLC, EC, HCC, UC, GC, CRC, NPC, OC, "Other" solid tumors in the Budha 2023 cohort).
- **Source aliases:**
  - `TUMTP` (categorical column with levels like `cHL`, `GC`, ...) — decompose into `TUMTP_CHL = as.integer(TUMTP == "cHL")`.
- **Example models:** `Budha_2023_tislelizumab.R`.
- **Notes:** Paired with `TUMTP_GC` in Budha 2023; a patient can have at most one of the indicators set to 1 (the remaining tumor types collapse into the reference 0 group).

### TUMTP_GC (**canonical for gastric-cancer tumor-type indicator**)
- **Description:** 1 = gastric cancer (GC), 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Reference category:** 0 = all other tumor types (same reference group as `TUMTP_CHL`).
- **Source aliases:**
  - `TUMTP` (categorical column) — decompose into `TUMTP_GC = as.integer(TUMTP == "GC")`.
- **Example models:** `Budha_2023_tislelizumab.R`.
- **Notes:** Follows the `RACE_<GROUP>` indicator-decomposition pattern. New oncology tumor types should be added as additional `TUMTP_<GROUP>` entries so the reference set stays explicit.

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
- **Example models:** `Ma_2020_sarilumab_das28crp.R` (multiplicative on DAS28-CRP Kout: `Kout * theta^PRICORT`); `Ma_2020_sarilumab_anc.R` (power-form on Emax: `Emax * 0.819^PRICORT`).
- **Notes:** Ma 2020 applies it as a multiplicative effect of the form `param * theta^PRICORT` in both DAS28-CRP and ANC PD models.

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

### ADA_TITER (**canonical**)
- **Description:** Continuous antidrug antibody titer (time-varying). Zero for ADA-negative observations.
- **Units:** (titer units; typically log2 or arbitrary assay units — document per-model in `covariateData[[ADA_TITER]]$units`).
- **Type:** continuous
- **Reference category:** n/a — ADA-negative corresponds to titer = 0.
- **Source aliases:**
  - `ADAT` — used in `Moein_2022_etrolizumab.R`.
- **Example models:** `Moein_2022_etrolizumab.R` (exponential effect on CL, per-unit-titer theta = 0.0365).
- **Notes:** Paired conceptually with `ADA_POS` (binary). When the paper reports both, the final model usually keeps only one. Imputation rules (LOCF / NOCB / baseline-as-negative) should be documented per-model.

## Disease / treatment history

### PRIOR_TNF (**canonical**)
- **Description:** 1 = subject previously treated with an anti-TNF (tumor necrosis factor) inhibitor, 0 = TNF-naive.
- **Units:** (binary)
- **Type:** binary
- **Reference category:** 0 (TNF-naive).
- **Source aliases:**
  - `PRIORTNF` (all caps, no underscore) — acceptable alternative spelling.
- **Example models:** `Moein_2022_etrolizumab.R` (multiplicative fractional effect on CL, +4.9%).
- **Notes:** Use when the source paper reports a binary "prior anti-TNF inhibitor" covariate on any PK parameter.

### DISEXT_EP (**canonical for extensive colitis / pancolitis indicator**)
- **Description:** 1 = extensive colitis or pancolitis disease extension, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Reference category:** 0 (left-sided colitis, when paired with `DISEXT_OTHER = 0`).
- **Source aliases:**
  - Derived from a multi-level `DISEXT` column in the source (levels: left-sided colitis, extensive/pancolitis, other): `DISEXT_EP = as.integer(DISEXT == "extensive/pancolitis")`.
- **Example models:** `Moein_2022_etrolizumab.R` (multiplicative effect on CL, +8.2% vs. left-sided colitis).
- **Notes:** Paired with `DISEXT_OTHER`. `DISEXT_EP = DISEXT_OTHER = 0` corresponds to the left-sided colitis reference group. UC-specific covariate; analogous indicators could be used for other inflammatory bowel disease models.

### DISEXT_OTHER (**canonical for 'other disease extension' indicator**)
- **Description:** 1 = disease extension other than left-sided colitis or extensive/pancolitis, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Reference category:** 0 (left-sided colitis, when paired with `DISEXT_EP = 0`).
- **Source aliases:** Derived from a multi-level `DISEXT` column: `DISEXT_OTHER = as.integer(DISEXT == "other")`.
- **Example models:** `Moein_2022_etrolizumab.R` (multiplicative effect on CL, +18% vs. left-sided colitis; large uncertainty due to 2% prevalence).
- **Notes:** Paired with `DISEXT_EP`; together they encode the three-level disease-extension categorical.

## Lifestyle / medical history

### SMOKE
- **Description:** 1 = current smoker at baseline, 0 = non-smoker.
- **Units:** (binary)
- **Type:** binary
- **Reference category:** 0 (non-smoker).
- **Source aliases:**
  - `Smoking` (case-insensitive) — used in `Ma_2020_sarilumab_anc.R`.
- **Example models:** `Ma_2020_sarilumab_anc.R` (power-form on baseline ANC: `BASE * 1.15^SMOKE`).
- **Notes:** Baseline-only indicator; does not track within-study smoking-cessation changes.

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
- **2026-04-19** — Added `ADA_TITER`, `PRIOR_TNF`, `DISEXT_EP`,
  `DISEXT_OTHER` canonical entries while extracting Moein 2022 etrolizumab.
  `ADA_TITER` is the continuous-titer companion to binary `ADA_POS`;
  `PRIOR_TNF` captures prior anti-TNF therapy (common covariate in IBD
  biologic PK models); `DISEXT_EP` / `DISEXT_OTHER` are paired indicators
  that encode the three-level ulcerative-colitis disease-extension
  categorical (reference = left-sided colitis). Extended the `CRP` entry
  to add Moein 2022 as a second example model.
- **2026-04-20** — Added `CRP` and `BMI` canonical entries for Chua 2025
  mirikizumab (VIVID-1 popPK). `CRP` is for papers that report a standard
  CRP assay where the baseline values are well above the hs-CRP sensitivity
  range, as is typical in moderate-to-severe IBD populations.
- **2026-04-20** — Added `SMOKE` canonical entry from the Ma 2020 sarilumab
  ANC PopPK/PD extraction. Binary baseline-only indicator used as a
  power-form covariate (`BASE * 1.15^SMOKE` on baseline ANC). Extended the
  `PRICORT` entry to record the ANC model as a second example.
- **2026-04-20** — Added `EASI` canonical entry for the Eczema Area and
  Severity Index, introduced by `Tiraboschi_2025_amlitelimab.R` (source
  alias `BEASI` for baseline EASI). Canonical name omits the `B` prefix to
  match the `AGE`/`WT`/`ALB` pattern where baseline vs time-varying status
  is noted in `covariateData[[...]]$notes` rather than the column name.
- **2026-04-20** — Added `TUMSZ`, `TUMTP_CHL`, `TUMTP_GC` canonical entries
  with the Budha 2023 tislelizumab extraction. `TUMSZ` centralizes the
  baseline-tumor-size continuous covariate; `TUMTP_<GROUP>` mirrors the
  `RACE_<GROUP>` decomposition so categorical tumor-type effects are
  stored as indicator columns with an explicit "all other tumor types"
  reference.
- Subsequent additions: append new canonical entries as new papers are processed. When adding, bump the audit-completed count in the summary below.
- **Xu 2019 sarilumab**: Added canonical entries `ALBR` (albumin / ULN ratio), `CRCL_BSA` (BSA-normalized creatinine clearance), `BLCRP` (baseline C-reactive protein), and `FORM_DP2` (sarilumab drug product 2 indicator). Extended the `ADA_POS` alias list to include the time-varying `ADA` column used in Xu 2019.
- **Ma 2020 sarilumab DAS28-CRP**: Added canonical entries `BLPHYVAS` (baseline Physician's Global Assessment of Disease Activity, 100-mm VAS), `BLHAQ` (baseline HAQ-DI), and `PRICORT` (prior corticosteroid treatment). Extended the `BLCRP` entry to record Ma 2020 as a second example model (reference 15.7 mg/L, covariate on DAS28-CRP BASE and log(Emax)).

## Summary

- Files audited: 63 R files under `inst/modeldb/` (14 of which reference covariates).
- Canonical entries: 49.
- Aliases mapped: 16 (including SEXM→SEXF, ADA→ADA_POS, ADAT→ADA_TITER, BLACK→RACE_BLACK, ASIAN→RACE_ASIAN, MULTIRACIAL→RACE_MULTI, BLACK_OTH→RACE_BLACK_OTH, ASIAN_AMIND_MULTI→RACE_ASIAN_AMIND_MULTI, DVID→STUDY1/STUDY5, CRE→CREAT, CRPHS→hsCRP, 1.73*CrCl/BSA→CRCL_BSA, DP2→FORM_DP2, DISEXT→DISEXT_EP/DISEXT_OTHER, BEASI→EASI, TUMTP→TUMTP_CHL/TUMTP_GC).
