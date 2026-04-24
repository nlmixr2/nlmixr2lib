# Canonical covariate columns

This file is the authoritative register of covariate column names used in nlmixr2lib models. Every covariate referenced inside a model's `model()` block must have an entry here. The register is seeded from a full audit of `inst/modeldb/` and extended whenever a new paper introduces a covariate that isn't yet registered.

## How to use this register

1. **Before adding a covariate to a new model**, search this file (by canonical name and by source alias) for the concept you need.
2. **If the canonical name exists**, use it exactly. Document any source-paper rename in the model's `covariateData[[name]]$source_name` field.
3. **If the source paper uses an alias listed under an existing canonical name**, use the canonical name and note in `covariateData[[name]]$notes` whether the values must be transformed (e.g., `SEXM → SEXF` inverts values; the effect coefficient sign / reference category must be inverted as well).
4. **If the covariate is not in this register at all**, propose a new entry with a canonical name, description, units, type, scope, and source aliases. Verify with the user before committing. The addition is part of the model's PR.
5. **Do not modify existing model files when you discover an alias**; simply document the mapping in the register. Retrofitting existing models is a separate effort.

## Scope: general vs specific

Each entry has a `Scope:` field declaring whether it is **general** (any model may use it without warning) or **specific** (only the models listed under `Example models` may use it; other usage is flagged by `checkModelConventions()`). This prevents accidentally reusing a covariate name whose meaning is tied to a particular paper.

- **Scope: general** — the covariate has a stable, paper-independent meaning and any future model may use it. Examples: `WT`, `AGE`, `SEXF`, `CREAT`, `ADA_POS`, `CRP`, `CRCL`.
- **Scope: specific** — the covariate's semantics depend on a particular study's design (a specific study indicator, a drug-product variant, a composite race grouping, a tumor-type decomposition, etc.). If a new paper needs the same concept with different semantics, register a new canonical name; if the concept matches, extend the `Example models` list (and consider promoting to general).

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
- **Scope:** general
- **Reference category:** n/a — used with allometric scaling `(WT / ref_wt)^exponent`. Reference weights observed: 70 kg (adults), 75 kg, 84.8 kg, 5 kg (infants).
- **Source aliases:** none known.
- **Example models:** `Clegg_2024_nirsevimab.R`, `Hu_2026_clesrovimab.R`, `Zhu_2017_lebrikizumab.R`, `Kovalenko_2020_dupilumab.R`, `CarlssonPetri_2021_liraglutide.R`, `Cirincione_2017_exenatide.R`, `Grimm_2023_gantenerumab.R`, `Grimm_2023_trontinemab.R`, `Kyhl_2016_nalmefene.R`, `Soehoel_2022_tralokinumab.R`, `Xie_2019_agomelatine.R`, `PK_2cmt_mAb_Davda_2014.R`, `phenylalanine_charbonneau_2021.R`, `Chua_2025_mirikizumab.R`, `Jackson_2022_ixekizumab.R`, `Kotani_2022_astegolimab.R`, `Ma_2020_sarilumab_anc.R`, `Ma_2020_sarilumab_das28crp.R`, `Moein_2022_etrolizumab.R`, `Tiraboschi_2025_amlitelimab.R`, `Robbie_2012_palivizumab.R`, `Bajaj_2017_nivolumab.R`.
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

### BSA
- **Description:** Body surface area (typically computed by DuBois, Mosteller, or Haycock from height and weight).
- **Units:** m^2
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(BSA / ref)^exponent`.
- **Source aliases:** none.
- **Example models:** `Yamada_2025_zolbetuximab.R` (reference 1.70 m^2; exponents 1.06 on clearances and 0.968 on volumes).
- **Notes:** Oncology mAbs dosed by BSA (mg/m^2) often use BSA in place of body weight for allometric-style scaling. Document the BSA computation formula (DuBois / Mosteller / Haycock) the source paper used; if unstated, record "unspecified."

### BMI
- **Description:** Body mass index at baseline.
- **Units:** kg/m²
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with a linear-deviation form (`1 + e * (BMI - ref)`) or a power form (`(BMI / ref)^e`). Document the reference value in `covariateData[[BMI]]$notes`.
- **Source aliases:** none known.
- **Example models:** `Chua_2025_mirikizumab.R` (reference 24.75 kg/m²; linear-deviation effect on logit of bioavailability).
- **Notes:** Universal clinical-trial demographic. Derived as `WT / (height_m)^2`; assume time-fixed at baseline unless the source paper states otherwise.

### SEXF (**canonical for sex**)
- **Description:** Biological sex indicator, 1 = female, 0 = male.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (male).
- **Source aliases:**
  - `SEXM` (values inverted: `SEXF = 1 - SEXM`; effect coefficient sign and reference category both invert) — used in `CarlssonPetri_2021_liraglutide.R`.
  - `SEX` with `"M"`/`"F"` strings — derive `SEXF = as.integer(SEX == "F")`.
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical), `CarlssonPetri_2021_liraglutide.R` (alias `SEXM`), `Bajaj_2017_nivolumab.R` (male-indicator source; effect applied as `exp(coef * (1 - SEXF))` to preserve the paper's female-reference CL_REF / VC_REF).
- **Notes:** When translating a model that used `SEXM`, flag the sign/reference-category inversion to the user.

### CHILD
- **Description:** 1 = subject is a child, 0 = not a child.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (not child, i.e., adult baseline).
- **Example models:** `CarlssonPetri_2021_liraglutide.R`.
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

## Renal / hepatic function

### CRCL (**canonical for creatinine-based renal function, BSA-normalized**)
- **Description:** Creatinine-based renal function expressed in mL/min/1.73 m². Accepts either an MDRD-/CKD-EPI-estimated glomerular filtration rate or a measured creatinine clearance that has been BSA-normalized as `1.73 × CrCl / BSA`. The per-model `covariateData[[CRCL]]$description` and `notes` must state which method the source paper used.
- **Units:** mL/min/1.73 m²
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(CRCL / ref)^exponent`. Reference values observed: 80 mL/min/1.73 m² (Cirincione 2017, MDRD eGFR), 90 mL/min/1.73 m² (Li 2019, calculated GFR), 100 mL/min/1.73 m² (Xu 2019, measured-CrCl BSA-normalized).
- **Reference category:** n/a — used with power scaling `(CRCL / ref)^exponent`. Reference values observed: 80 mL/min/1.73 m² (Cirincione 2017, MDRD eGFR), 100 mL/min/1.73 m² (Xu 2019, measured-CrCl BSA-normalized), 90 mL/min/1.73 m² (Bajaj 2017, CKD-EPI eGFR).
- **Source aliases:**
  - `eGFR` — MDRD-estimated glomerular filtration rate; used in `Cirincione_2017_exenatide.R` and `Kotani_2022_astegolimab.R`. `Bajaj_2017_nivolumab.R` uses the CKD-EPI variant.
  - `EGFR` — all-caps variant.
  - `CRCL_BSA` — BSA-normalized creatinine clearance (measured CrCl ÷ BSA × 1.73); used in `Xu_2019_sarilumab.R`.
  - `1.73*CrCl/BSA` — the formula form appearing in Xu 2019 Eq. for Vm.
  - `cGFR` — calculated/estimated GFR, BSA-normalized; used in `Li_2019_abatacept.R`.
- **Example models:** `Cirincione_2017_exenatide.R` (MDRD eGFR), `Xu_2019_sarilumab.R` (measured CrCl BSA-normalized), `Kotani_2022_astegolimab.R` (MDRD eGFR), `Li_2019_abatacept.R` (cGFR).
- **Example models:** `Cirincione_2017_exenatide.R` (MDRD eGFR), `Xu_2019_sarilumab.R` (measured CrCl BSA-normalized), `Kotani_2022_astegolimab.R` (MDRD eGFR), `Bajaj_2017_nivolumab.R` (CKD-EPI eGFR, reference 90 mL/min/1.73 m²).
- **Notes:** The two estimation methods (MDRD/CKD-EPI vs measured CrCl) produce values in the same units and are operationally interchangeable as a covariate on clearance. Document the method explicitly in each model's `covariateData[[CRCL]]$description` so future reviewers can trace the source assay.

### CREAT (**canonical for serum creatinine**)
- **Description:** Serum creatinine concentration (baseline or time-varying).
- **Units:** µmol/L or mg/dL — document the unit used in each model via `covariateData[[CREAT]]$units`.
- **Type:** continuous
- **Scope:** general
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
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(ALB / ref)^exponent`.
- **Source aliases:** none; `ALB` is the universal abbreviation.
- **Example models:** `Fasanmade_2009_infliximab.R` (g/dL, reference 4.1), `Thakre_2022_risankizumab.R` (g/L, reference 45), `Chua_2025_mirikizumab.R`, `Moein_2022_etrolizumab.R`, `Tiraboschi_2025_amlitelimab.R`, `Yamada_2025_zolbetuximab.R`, `Li_2019_abatacept.R` (g/dL, reference 4.0; the Li 2019 Methods states 'mg/dL' which is a publication typo — see the model's `covariateData[[ALB]]$notes`).
- **Notes:** Ratified canonically on 2026-04-19 after cross-model review. Unit varies by paper (g/dL in US-convention papers, g/L in SI-convention papers); the per-model `covariateData[[ALB]]$units` field is load-bearing. Effect-coefficient magnitude is meaningless without the unit.

### TBILI (**canonical for total bilirubin**)
- **Description:** Total serum bilirubin concentration.
- **Units:** mg/dL or umol/L — document the unit used in each model via `covariateData[[TBILI]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(TBILI / ref)^exponent`.
- **Source aliases:** none; `TBILI` is the common NONMEM / clinical-PK abbreviation.
- **Example models:** `Yamada_2025_zolbetuximab.R` (mg/dL, reference 0.38; small positive exponent 0.0347 on V1).
- **Notes:** Hepatic-function marker. Unit varies by paper (US convention mg/dL, SI convention umol/L; 1 mg/dL ~= 17.1 umol/L). The per-model `covariateData[[TBILI]]$units` field is load-bearing.

### AST (**canonical for aspartate aminotransferase**)
- **Description:** Serum aspartate aminotransferase activity (baseline or time-varying).
- **Units:** U/L (IU/L; the two labels are used interchangeably in the clinical-PK literature). Document per-model via `covariateData[[AST]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(AST / ref)^exponent`.
- **Source aliases:** none; `AST` is the universal clinical-PK abbreviation.
- **Example models:** `Lu_2014_trastuzumabemtansine.R` (U/L, reference 27; small positive exponent 0.071 on CL).
- **Notes:** Hepatic-function marker. Commonly reported alongside `ALT` and `TBILI`; register a separate `ALT` canonical if a future paper requires it.

## Hematology

### HGB (**canonical for hemoglobin**)
- **Description:** Blood hemoglobin concentration.
- **Units:** g/L or g/dL — document the unit used in each model via `covariateData[[HGB]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(HGB / ref)^exponent`.
- **Source aliases:** none; `HGB` is the common NONMEM / clinical-PK abbreviation.
- **Example models:** `Yamada_2025_zolbetuximab.R` (g/L, reference 118; exponent -0.374 on V1).
- **Notes:** Unit varies by paper (SI g/L, US g/dL; 1 g/dL = 10 g/L). The per-model `covariateData[[HGB]]$units` field is load-bearing.

### WBC (**canonical for white blood cell count**)
- **Description:** Total white blood cell count (baseline or time-varying). In chronic lymphocytic leukaemia (CLL) populations the value is elevated because circulating leukaemic B-cells make up the majority of the count, so WBC can serve as a biomarker of target-cell burden rather than general hematology.
- **Units:** 10^9 cells/L (equivalent to 10^3 cells/µL). Document per-model via `covariateData[[WBC]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(WBC / ref)^exponent`. Reference values observed: 10 × 10^9/L (Mould 2007, typical CLL Vmax normalization).
- **Source aliases:** none known (`WBC` is the universal clinical-PK abbreviation).
- **Example models:** `Mould_2007_alemtuzumab.R` (reference 10 × 10^9/L; exponent 0.194 on Vmax).
- **Notes:** Time-varying in treatment studies where the drug depletes the leukaemic clone (e.g., alemtuzumab in CLL): WBC must be supplied at every observation time in the event dataset. In diseases where WBC is not therapeutically targeted the column can be treated as a baseline-only covariate; record the per-model convention in `covariateData[[WBC]]$notes`.

## Disease severity scores

### EASI (**canonical for Eczema Area and Severity Index**)
- **Description:** Eczema Area and Severity Index score (atopic-dermatitis severity composite; bounded continuous, scale 0-72 with higher values = more severe disease).
- **Units:** (score)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — healthy volunteers have EASI = 0. Effect enters as an additive term in models that pool AD patients with HV (e.g., `Tiraboschi_2025_amlitelimab.R`).
- **Source aliases:**
  - `BEASI` (baseline EASI) — used in `Tiraboschi_2025_amlitelimab.R`.
- **Example models:** `Tiraboschi_2025_amlitelimab.R`.
- **Notes:** When used as a time-invariant baseline covariate (`BEASI`), document in `covariateData[[EASI]]$notes`. Canonical name is `EASI` without the `B` prefix to match the `AGE` / `WT` / `ALB` pattern where baseline vs time-varying status is recorded in notes rather than the column name.

### MGADL (**canonical for Myasthenia Gravis Activities of Daily Living score**)
- **Description:** Myasthenia Gravis Activities of Daily Living score — eight-item patient-reported outcome measure (each item 0-3), total 0-24, higher values = greater symptom severity and functional limitation.
- **Units:** (score)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — healthy participants (no gMG) have `MGADL = 0` by definition. Effect enters as a baseline covariate on MG-ADL response parameters in gMG cohorts.
- **Source aliases:** none known.
- **Example models:** `Valenzuela_2025_nipocalimab.R` (reference 7 points; power-form effect on `IDecplacebo` and on the slope between MG-ADL change and IgG reduction).
- **Notes:** Baseline-only in Valenzuela 2025 (the observation is the absolute change from baseline MG-ADL). When used time-varying (e.g., in pure PD models driven by disease-progression dynamics), document in `covariateData[[MGADL]]$notes`. Canonical name is `MGADL` without a `BL` prefix to match the `EASI` / `AGE` / `WT` / `ALB` pattern.

### BCVA (**canonical for best-corrected visual acuity**)
- **Description:** Best-corrected visual acuity score measured on the Early Treatment Diabetic Retinopathy Study (ETDRS) chart, expressed as the number of letters read correctly (0-100; higher values = better vision). Used as a baseline severity covariate in ophthalmology PK/PD models of anti-VEGF treatment.
- **Units:** ETDRS letters (0-100)
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used as a baseline input to set the initial condition of an indirect-response BCVA state or as a power-form effect on response parameters. Reference value observed: 55 letters (Mulyukov 2018 narrative: mean study-population baseline BCVA).
- **Source aliases:**
  - `BVA` (baseline visual acuity) — used in `Mulyukov_2018_ranibizumab.R`.
- **Example models:** `Mulyukov_2018_ranibizumab.R` (baseline BCVA used as the center for the initial-condition draw `g0 = BCVA + eta_g0`).
- **Notes:** Ophthalmology-specific. Baseline-only in Mulyukov 2018 (carried once per subject and used only as the starting BCVA for the indirect-response model). Canonical name drops the `B` prefix to match the `EASI` / `AGE` / `WT` / `ALB` pattern (baseline-vs-time-varying status recorded in `covariateData[[BCVA]]$notes`). Scope is `specific` until a second ophthalmology model ratifies the name; at that point promote to `general`.

## Interferon / biomarker panels

### BGENE21
- **Description:** Baseline type I interferon gene signature computed from the expression of 21 IFN-alpha / type-I-IFN-inducible transcripts in peripheral blood mononuclear cells. Used as a continuous biomarker of type I IFN pathway activation in SLE (and mechanistically-related indications).
- **Units:** (gene-signature score, unitless)
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used with power scaling `(BGENE21 / ref)^exponent`. Reference value observed: 32 in Narwal 2013 (study-population median was 33).
- **Source aliases:** none.
- **Example models:** `Narwal_2013_sifalimumab.R` (reference 32, exponent 0.0558 on CL).
- **Notes:** Specific to studies that report a 21-gene IFN signature (subset / expansion of the "BGENE4" 4-gene signature). If another paper uses the 4-gene signature or a differently-constructed IFN score, register a distinct canonical (`BGENE4`, `IFN_SIG`, ...).

## Inflammation markers

### EOS (**canonical for blood eosinophil count**)
- **Description:** Blood eosinophil count (baseline or time-varying).
- **Units:** cells/µL
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(EOS / ref)^exponent`.
- **Source aliases:**
  - `BEOS` (baseline EOS) — used in `Kotani_2022_astegolimab.R`.
- **Example models:** `Kotani_2022_astegolimab.R` (reference 180 cells/µL, baseline).
- **Notes:** Used as a surrogate of inflammatory burden that correlates with protein turnover and therefore mAb clearance. Canonical name drops the `B` prefix to match the `EASI` / `AGE` / `WT` / `ALB` pattern; baseline-vs-time-varying status is documented in `covariateData[[EOS]]$notes`.

### BGENE21 (**canonical for 21-gene type I interferon signature score**)
- **Description:** Baseline 21-gene type I interferon signature score — a composite transcriptomic score summarising the expression of 21 interferon-regulated genes in whole blood relative to a healthy-donor reference, used as a biomarker of type I IFN pathway activation in SLE and related autoimmune conditions.
- **Units:** unitless fold-change score (relative to healthy-donor reference).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used with power scaling `(BGENE21 / ref)^exponent`. Reference value observed: 12.04 in `Zheng_2016_sifalimumab.R` (median of the SLE phase IIb cohort, range 0.32-38.59).
- **Source aliases:** none.
- **Example models:** `Zheng_2016_sifalimumab.R` (power effect on CL with exponent 0.09).
- **Notes:** Specific to drugs whose mechanism targets the type I IFN pathway (e.g., anti-IFN-alpha antibodies like sifalimumab, anifrolumab). Higher BGENE21 indicates stronger target engagement / disease activity and is associated with increased drug clearance via target-mediated mechanisms. Scope: specific because the 21-gene panel composition is tied to the MedImmune/AstraZeneca SLE development programme; a different IFN gene signature (e.g., a 4-gene or 5-gene panel) should be registered under its own canonical name to avoid conflating panel definitions.

### CRP (**canonical for C-reactive protein**)
- **Description:** C-reactive protein concentration. Covers both standard and high-sensitivity (hs-CRP) assays and both baseline and time-varying usages. Each model's `covariateData[[CRP]]$description` and `notes` must state the assay type (standard vs hs-CRP) and whether the column carries a baseline-only or time-varying value, including the paper-specific reference value used for power scaling.
- **Units:** mg/L (document per-model via `covariateData[[CRP]]$units`).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(CRP / ref)^exponent` or exponential effects `exp(coef * (CRP - ref))`. Reference values observed: 4.23 mg/L (Moein 2022, IBD standard assay), 4.31 mg/L (Moein 2022 Table 3 median), 5.21 mg/L (Thakre 2022, baseline hs-CRP), 7.41 mg/L (Chua 2025, baseline standard assay), 14.2 mg/L (Xu 2019, baseline standard assay), 15.7 mg/L (Ma 2020, baseline standard assay).
- **Source aliases:**
  - `hsCRP` — high-sensitivity CRP (mixed-case preserved from earlier register drafts).
  - `HSCRP` — all-caps variant.
  - `CRPHS` — used in `Thakre_2022_risankizumab.R` (baseline, high-sensitivity assay).
  - `BLCRP` — baseline CRP; used in `Xu_2019_sarilumab.R` and `Ma_2020_sarilumab_das28crp.R`.
- **Example models:** `Thakre_2022_risankizumab.R`, `Xu_2019_sarilumab.R`, `Chua_2025_mirikizumab.R`, `Moein_2022_etrolizumab.R`, `Ma_2020_sarilumab_das28crp.R`.
- **Notes:** The prior separate `hsCRP`, `BLCRP`, and standard-assay `CRP` canonicals were merged on 2026-04-20 to a single general-scope `CRP` canonical. Assay type (standard vs hs-CRP), baseline-vs-time-varying status, and the paper-specific reference value all live in each model's `covariateData[[CRP]]$description` / `notes`. Only aggregate values from hs-validated assays as CRP when the downstream analysis relies on low-range sensitivity; for most inflammatory-disease cohorts (IBD, RA/PsA), baseline CRP is well above the hs-sensitivity range and the distinction is moot.

## Race / ethnicity

**Canonical pattern: `RACE_<GROUP>`.** Use one indicator per race/ethnicity group the source models. Reference category is the implicit 0 = all other groups; document explicitly which groups are in the reference. When the source uses composite groups (e.g., "Black or Other"), name them accordingly (`RACE_BLACK_OTHER`) and list the components in `notes`. The base `RACE_<GROUP>` indicators are scope: general; composite groupings are scope: specific because the grouping is tied to the study's analysis plan.

### RACE_BLACK (**canonical**)
- **Description:** 1 = Black / African American, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (document the actual reference groups used).
- **Source aliases:**
  - `BLACK` — used in `Hu_2026_clesrovimab.R`, `Robbie_2012_palivizumab.R`.
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical form), `Robbie_2012_palivizumab.R`.

### RACE_BLACK_OTH (**canonical for composite Black/Other group**)
- **Description:** 1 = Black/African American or Other race, 0 = other groups.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = White or Native Hawaiian/Pacific Islander (Clegg 2024 grouping).
- **Source aliases:** `BLACK_OTH` — used in `Clegg_2024_nirsevimab.R`.
- **Example models:** `Clegg_2024_nirsevimab.R`.
- **Notes:** Kept distinct from `RACE_BLACK` because the composite is not interchangeable.

### RACE_ASIAN (**canonical**)
- **Description:** 1 = Asian, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Source aliases:** `ASIAN` — used in `Hu_2026_clesrovimab.R`, `Robbie_2012_palivizumab.R`. `RAAS` (race-Asian-vs-other indicator as named in Bajaj 2017 Table 1) — used in `Bajaj_2017_nivolumab.R`.
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical form), `Robbie_2012_palivizumab.R`, `Bajaj_2017_nivolumab.R`.

### RACE_ASIAN_AMIND_MULTI (**canonical for composite group**)
- **Description:** 1 = Asian, American Indian / Alaskan Native, or Multiple races, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = White, Black / African American, Native Hawaiian / Pacific Islander, or Other (Clegg 2024 grouping).
- **Source aliases:** `ASIAN_AMIND_MULTI` — used in `Clegg_2024_nirsevimab.R`.
- **Example models:** `Clegg_2024_nirsevimab.R`.
- **Notes:** Clegg 2024 applies this covariate to both CL and V2 with different coefficients.

### RACE_MULTI
- **Description:** 1 = multiracial, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Source aliases:** `MULTIRACIAL` — used in `Hu_2026_clesrovimab.R`.
- **Example models:** `Hu_2026_clesrovimab.R`.

### RACE_OTHER
- **Description:** 1 = race category "Other," 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Source aliases:**
  - `OTHER` — used in `Robbie_2012_palivizumab.R`.
- **Example models:** `Zhu_2017_lebrikizumab.R`, `Robbie_2012_palivizumab.R`.

### RACE_HISPANIC
- **Description:** 1 = Hispanic / Latino, 0 = non-Hispanic. Used by papers that report Hispanic as a separate category alongside Black, Asian, and Other rather than as a distinct ethnicity dimension.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (non-Hispanic; document the paper-specific reference race composition per-model).
- **Source aliases:**
  - `HISPANIC` — used in `Robbie_2012_palivizumab.R`.
- **Example models:** `Robbie_2012_palivizumab.R` (fractional effect on CL; additional effect on Vc).
- **Notes:** In the US Office-of-Management-and-Budget (OMB) classification Hispanic is an ethnicity rather than a race, but clinical PK analyses frequently treat it as one of the race indicators. When a paper treats Hispanic as a race, use this column; otherwise encode ethnicity separately. Register-wise, this follows the `RACE_<GROUP>` indicator-decomposition pattern.

## Pediatric comorbidities

### CLD_PREM (**canonical for chronic lung disease of prematurity**)
- **Description:** 1 = chronic lung disease of prematurity (bronchopulmonary dysplasia, BPD), 0 = no CLD. Time-fixed per subject (diagnosis at study entry).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no CLD of prematurity).
- **Source aliases:**
  - `CLD` — used in `Robbie_2012_palivizumab.R`.
  - `BPD` — bronchopulmonary-dysplasia shorthand.
- **Example models:** `Robbie_2012_palivizumab.R` (fractional +20% effect on CL).
- **Notes:** Standard pediatric / neonatology comorbidity flag; ties to palivizumab's label population (high-risk preterm infants) and may re-appear in future pediatric mAb PK analyses (RSV, parenteral nutrition, etc.).

## Surgical history / disease state

### PRIOR_GAST (**canonical for prior gastrectomy**)
- **Description:** Prior (partial or total) gastrectomy indicator, 1 = prior gastrectomy, 0 = no prior gastrectomy. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no prior gastrectomy).
- **Source aliases:**
  - `GAST` — used in `Yamada_2025_zolbetuximab.R`.
- **Example models:** `Yamada_2025_zolbetuximab.R` (fractional effects on CLss, CLT, V1).
- **Notes:** Renamed from `GAST` on 2026-04-20 to follow the `PRIOR_TNF` / `PRICORT` naming pattern for prior-treatment and surgical-history indicators. Applicable to any PK model where gastrointestinal anatomy affects absorption, first-pass, or protein turnover; not inherently oncology-specific. No distinction between partial vs total gastrectomy unless the source paper separates them.

## Disease state (cross-population indicators)

### DIS_UC (**canonical for ulcerative colitis disease-state indicator**)
- **Description:** 1 = ulcerative colitis patient, 0 = non-UC (e.g., healthy volunteer or non-IBD indication). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-UC subject; the complement group is defined per-model — typically healthy volunteers and/or patients with another indication such as asthma).
- **Source aliases:**
  - `UC` — used in `Hua_2015_anrukinzumab.R`.
- **Example models:** `Hua_2015_anrukinzumab.R` (multiplicative fractional increase in CL, +72.8%, on top of weight and albumin effects).
- **Notes:** Used when a population PK model pools UC patients with a non-UC reference population (e.g., Hua 2015: healthy volunteers + asthma patients + UC patients) and UC disease status is tested as a PK covariate. Distinct from `DISEXT_EP` / `DISEXT_OTHER`, which operate *within* a UC-only cohort (disease extension). Start as scope: specific; promote to general if a second paper pools UC with a non-UC reference.

### DIS_SASTHMA (**canonical for moderate-to-severe asthma disease-state indicator**)
- **Description:** 1 = moderate-to-severe asthma patient, 0 = not (e.g., healthy volunteer, mild-to-moderate asthma, or other indication). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-moderate-to-severe-asthma subject; the complement group is defined per-model).
- **Source aliases:**
  - `sAsthma` — used in `Hua_2015_anrukinzumab.R`.
- **Example models:** `Hua_2015_anrukinzumab.R` (multiplicative fractional change in SC bioavailability, -30.9%).
- **Notes:** The moderate-to-severe vs. mild-to-moderate asthma cutoff is protocol-defined; Hua 2015 uses FEV1 55-80% and ACQ-5 >= 2 for "moderate to severe" (study 4) versus FEV1 > 70% and ACQ-5 <= 1 for "mild to moderate" (study 1). Scope: specific because the severity threshold is tied to a particular analysis plan; future asthma-severity indicators with different thresholds should register as separate canonicals.

## Oncology

### TUMSZ (**canonical for baseline tumor size**)
- **Description:** Baseline tumor size. For solid tumors, the sum of diameters of target lesions per RECIST; for classical Hodgkin lymphoma, the sum of products of perpendicular diameters (SPPD).
- **Units:** mm
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(TUMSZ / ref)^exponent`. Reference values observed: 63 mm (Budha 2023); 90 mm (Lu 2014, source reference 9 cm converted to mm).
- **Source aliases:**
  - `TMBD` (originally in cm; `TUMSZ_mm = TMBD_cm * 10`) — used in `Lu_2014_trastuzumabemtansine.R`.
- **Example models:** `Budha_2023_tislelizumab.R` (reference 63 mm), `Lu_2014_trastuzumabemtansine.R` (reference 90 mm; source column TMBD in cm, values converted to mm on ingestion).
- **Notes:** Promoted to scope: general on 2026-04-20 as a conventional oncology baseline-tumor-size measure (RECIST for solid tumors, SPPD for cHL). The SPPD vs sum-of-diameters convention is pooled onto a single column; document the per-model mixture where relevant. When the source paper reports tumor size in cm, convert to mm (the canonical unit) on data ingestion and scale the per-model reference accordingly so `(TUMSZ / ref)^exp` is numerically invariant.

### TUMTP_CHL (**canonical for classical Hodgkin lymphoma tumor-type indicator**)
- **Description:** 1 = classical Hodgkin lymphoma (cHL), 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (e.g., NSCLC, EC, HCC, UC, GC, CRC, NPC, OC, "Other" solid tumors in the Budha 2023 cohort).
- **Source aliases:**
  - `TUMTP` (categorical column with levels like `cHL`, `GC`, ...) — decompose into `TUMTP_CHL = as.integer(TUMTP == "cHL")`.
- **Example models:** `Budha_2023_tislelizumab.R`.
- **Notes:** Paired with `TUMTP_GC` in Budha 2023; a patient can have at most one of the indicators set to 1 (the remaining tumor types collapse into the reference 0 group).

### HER2_ECD (**canonical for HER2 shed extracellular domain concentration**)
- **Description:** Baseline serum concentration of the shed extracellular domain of human epidermal growth factor receptor 2 (HER2). Serves as a soluble-antigen biomarker of HER2-mediated target-mediated drug disposition for HER2-directed mAbs / ADCs.
- **Units:** ng/mL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used with power scaling `(HER2_ECD / ref)^exponent`. Reference 25 ng/mL used in Lu 2014.
- **Source aliases:**
  - `ECD` — used in `Lu_2014_trastuzumabemtansine.R`.
- **Example models:** `Lu_2014_trastuzumabemtansine.R` (reference 25 ng/mL; exponent 0.035 on CL).
- **Notes:** Scoped specific because the covariate is meaningful only for HER2-targeted agents; if a non-HER2 paper uses a shed-antigen analog for a different target, register a target-specific canonical (e.g., `EGFR_ECD`) rather than reusing this one. Disambiguated from the covariate-columns register by the explicit `HER2_` prefix.

### TRAST_BL (**canonical for baseline trastuzumab concentration from prior therapy**)
- **Description:** Baseline serum concentration of residual unconjugated trastuzumab remaining from prior trastuzumab-containing therapy, measured at the start of a subsequent anti-HER2 treatment (e.g., trastuzumab emtansine). Encodes the magnitude of residual HER2-site competition from a previous trastuzumab exposure.
- **Units:** ug/mL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used as a linear covariate on log(CL) via `exp(coef * TRAST_BL)`. TRAST_BL = 0 corresponds to no detectable residual trastuzumab (reference condition used in Lu 2014).
- **Source aliases:**
  - `TBL` — used in `Lu_2014_trastuzumabemtansine.R`.
- **Example models:** `Lu_2014_trastuzumabemtansine.R` (linear-on-log coefficient -0.002 per ug/mL on CL).
- **Notes:** Scoped specific because the covariate is meaningful only for drugs that compete with trastuzumab at the HER2 binding site. Expect values clustered at 0 for trastuzumab-naive patients; Lu 2014 observed 0 at the 5th percentile and 54 ug/mL at the 95th percentile.

### TUMTP_GC (**canonical for gastric-cancer tumor-type indicator**)
- **Description:** 1 = gastric cancer (GC), 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (same reference group as `TUMTP_CHL`).
- **Source aliases:**
  - `TUMTP` (categorical column) — decompose into `TUMTP_GC = as.integer(TUMTP == "GC")`.
- **Example models:** `Budha_2023_tislelizumab.R`.
- **Notes:** Follows the `RACE_<GROUP>` indicator-decomposition pattern. New oncology tumor types should be added as additional `TUMTP_<GROUP>` entries so the reference set stays explicit.

### ECOG_GE1 (**canonical for Eastern Cooperative Oncology Group performance-status indicator, >= 1**)
- **Description:** 1 if baseline Eastern Cooperative Oncology Group (ECOG) performance status is greater than or equal to 1, 0 if ECOG = 0.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (ECOG performance status = 0, i.e., fully active / asymptomatic).
- **Source aliases:**
  - `PS` / `BPS` — the Bajaj 2017 nivolumab analysis reports "baseline performance status" (BPS) as a binary ECOG-derived indicator. In Bajaj 2017 the one study using Karnofsky Performance Status (KPS) values was mapped to the ECOG scale per Oken 1982 before thresholding.
  - `ECOG_1` — alternative explicit form; equivalent to `ECOG_GE1` when ECOG only takes values 0, 1, 2 in the analysis dataset (which is the typical oncology case).
- **Example models:** `Bajaj_2017_nivolumab.R` (exponential effect on CL with coefficient 0.172, reference 0).
- **Notes:** Oncology papers conventionally report ECOG as an integer (0-5) but binarize at >= 1 because ECOG >= 2 is rare in trial cohorts. When a source paper provides the ordinal ECOG score separately, derive `ECOG_GE1 = as.integer(ECOG >= 1)`. If a future paper needs finer resolution (e.g., separate effects for ECOG 1 vs ECOG 2), add a parallel `ECOG_GE2` canonical rather than overloading this one.

## Laboratory / disease-activity

### ALBR
- **Description:** Serum albumin normalized to the laboratory's upper limit of normal (`albumin_observed / ULN_albumin`).
- **Units:** (unitless ratio)
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used as a power term `(ALBR / <ref>)^exponent`. Reference 0.78 used in Xu 2019 (corresponds to a median serum albumin of 38 g/L at a typical ULN of ~48.7 g/L).
- **Source aliases:** none.
- **Example models:** `Xu_2019_sarilumab.R`.
- **Notes:** Xu 2019 normalizes to each site's ULN so that values across multiple labs with different reference ranges can be pooled. Scoped specific because the ULN-normalization convention is tied to the Xu 2019 analysis plan; future papers using the same ratio should either add themselves to the example_models list or promote this entry to general.

## Inflammatory-bowel-disease disease-activity covariates

### CALPRO (**canonical for fecal calprotectin**)
- **Description:** Fecal calprotectin, a gut-inflammation biomarker (baseline or time-fixed per subject unless a paper explicitly uses a time-varying value).
- **Units:** mg/kg stool (equivalent to µg/g). Document per-model via `covariateData[[CALPRO]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(CALPRO / ref)^exponent`. Reference 700 mg/kg used in Rosario 2015 (overall population median).
- **Source aliases:** none.
- **Example models:** `Rosario_2015_vedolizumab.R` (reference 700 mg/kg; exponent +0.0310 on linear clearance CLL).
- **Notes:** Common IBD severity biomarker (inflammation of the gut epithelium). Assays typically report in µg/g stool; 1 µg/g = 1 mg/kg. Document per-model whether baseline-only or time-varying in `covariateData[[CALPRO]]$notes`.

### CDAI (**canonical for Crohn's Disease Activity Index**)
- **Description:** Crohn's Disease Activity Index composite score. Higher values indicate more active disease; <150 remission, 150–219 mild, 220–450 moderate, >450 severe. Defined only for patients with a CD diagnosis; set to the reference value (or gate via the indicator) for UC patients.
- **Units:** (score, 0–600)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(CDAI / ref)^exponent`. Reference 300 used in Rosario 2015 (typical moderate-CD score).
- **Source aliases:** none.
- **Example models:** `Rosario_2015_vedolizumab.R` (reference 300; exponent −0.0515 on CLL gated by `IBD_CD` so the effect applies only to CD patients).
- **Notes:** Mutually exclusive with `PMAYO` in pooled UC+CD populations: each patient has exactly one disease-activity score (CDAI for CD, partial Mayo for UC). Gate via the `IBD_CD` indicator when pooling.

### PMAYO (**canonical for partial Mayo score**)
- **Description:** Partial Mayo score for ulcerative colitis (sum of stool-frequency, rectal-bleeding, and physician-global subscores, range 0–9). Higher values indicate more active disease. Defined only for patients with a UC diagnosis; gate via the `IBD_CD` indicator when pooling UC+CD.
- **Units:** (score, 0–9)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(PMAYO / ref)^exponent`. Reference 6 used in Rosario 2015 (typical moderate-UC score).
- **Source aliases:** none.
- **Example models:** `Rosario_2015_vedolizumab.R` (reference 6; exponent +0.0408 on CLL gated by `(1 - IBD_CD)` so the effect applies only to UC patients).
- **Notes:** The partial Mayo score excludes the endoscopy subscore (the full Mayo score is 0–12). Mutually exclusive with `CDAI` in pooled UC+CD populations.

## Inflammatory-bowel-disease diagnosis

### IBD_CD (**canonical for Crohn's disease indicator**)
- **Description:** 1 = Crohn's disease (CD) diagnosis, 0 = ulcerative colitis (UC) diagnosis. Used to gate pooled UC+CD population models where some parameters or covariate effects differ by diagnosis.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (ulcerative colitis).
- **Source aliases:**
  - `DIAGNOSIS`, `DX` (categorical `"UC"`/`"CD"`) — derive `IBD_CD = as.integer(DX == "CD")`.
- **Example models:** `Rosario_2015_vedolizumab.R` (two typical-CLL switch between UC and CD; gates `PMAYO` and `CDAI` effects; multiplicative +1% effect on Vc).
- **Notes:** Rosario 2015 models separate typical CLL for UC vs CD and gates the partial-Mayo (UC-only) and CDAI (CD-only) disease-activity covariates via this indicator.

## Concomitant IBD medications

### CONMED_AZA (**canonical for concomitant azathioprine**)
- **Description:** 1 = on concomitant azathioprine at the PK observation, 0 = not on azathioprine.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant azathioprine).
- **Source aliases:** `AZA` — used in `Rosario_2015_vedolizumab.R`.
- **Example models:** `Rosario_2015_vedolizumab.R` (power-form on CLL: `CLL * 0.998^CONMED_AZA`; effect ≈ null).
- **Notes:** Thiopurine immunomodulator used as maintenance therapy in IBD. Standard convention is baseline-use-only, but time-varying use is permitted; document per-model.

### CONMED_MP (**canonical for concomitant 6-mercaptopurine**)
- **Description:** 1 = on concomitant 6-mercaptopurine (6-MP), 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant 6-MP).
- **Source aliases:** `MP` — used in `Rosario_2015_vedolizumab.R`.
- **Example models:** `Rosario_2015_vedolizumab.R` (power-form on CLL: `CLL * 1.04^CONMED_MP`).
- **Notes:** Second thiopurine immunomodulator used in IBD maintenance; typically mutually exclusive with `CONMED_AZA` for a given subject.

### CONMED_MTX (**canonical for concomitant methotrexate**)
- **Description:** 1 = on concomitant methotrexate, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant methotrexate).
- **Source aliases:** `MTX` — used in `Rosario_2015_vedolizumab.R`.
- **Example models:** `Rosario_2015_vedolizumab.R` (power-form on CLL: `CLL * 0.983^CONMED_MTX`).
- **Notes:** Immunomodulator used especially in CD maintenance. Generic concomitant-MTX indicator that may also appear in non-IBD models; start as scope: general.

### CONMED_AMINO (**canonical for concomitant aminosalicylate therapy**)
- **Description:** 1 = on concomitant aminosalicylate (5-aminosalicylic acid / mesalamine / mesalazine / olsalazine / sulfasalazine etc.) therapy at the PK observation, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant aminosalicylate).
- **Source aliases:** `AMINO` — used in `Rosario_2015_vedolizumab.R`.
- **Example models:** `Rosario_2015_vedolizumab.R` (power-form on CLL: `CLL * 1.02^CONMED_AMINO`).
- **Notes:** Covers the full aminosalicylate class (5-ASA is the single active moiety shared by most agents); use `CONMED_AMINO` rather than `CONMED_5ASA` unless the source paper explicitly restricts the indicator to 5-ASA monotherapy.

### CONMED_NSAID (**canonical for concomitant NSAID use**)
- **Description:** 1 = on concomitant non-steroidal anti-inflammatory drug (NSAID) therapy at baseline, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant NSAID use; typical patient).
- **Source aliases:**
  - `NSAID` — used in `Li_2019_abatacept.R`.
- **Example models:** `Li_2019_abatacept.R` (exponential effect on CL: `CL * exp(0.0640 * CONMED_NSAID)`; ~6.6% higher CL, not clinically relevant per Li 2019).
- **Notes:** Baseline-use-only in Li 2019; time-varying use is permitted, document per-model. Follows the `CONMED_*` concomitant-medication pattern established for IBD models (AZA / MP / MTX / AMINO).

## Rheumatoid-arthritis disease-activity covariates

### BLPHYVAS
- **Description:** Baseline Physician's Global Assessment of Disease Activity, 100-mm visual analogue scale (0 = no disease activity, 100 = maximum). Time-fixed per subject.
- **Units:** mm (0-100 VAS)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used as a power term `(BLPHYVAS / <ref>)^exponent`. Reference 66 used in Ma 2020.
- **Source aliases:** none.
- **Example models:** `Ma_2020_sarilumab_das28crp.R`.
- **Notes:** One of the components of the DAS28 composite score; in Ma 2020 it appears as a baseline covariate on the DAS28-CRP disease-activity BASE rather than on the score itself. Applicable to any rheumatology model where baseline physician-assessed disease activity is used as a PK/PD covariate.

### BLHAQ
- **Description:** Baseline Health Assessment Questionnaire Disability Index (HAQ-DI; 0 = no disability, 3 = maximum disability). Time-fixed per subject.
- **Units:** unitless (0-3 composite score)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used as a power term `(BLHAQ / <ref>)^exponent`. Reference 1.75 used in Ma 2020.
- **Source aliases:** none.
- **Example models:** `Ma_2020_sarilumab_das28crp.R`.
- **Notes:** Patient-reported disability score frequently used as a baseline covariate in rheumatoid-arthritis PK/PD analyses.

### SWOL_28JOINT (**canonical for 28-joint swollen joint count**)
- **Description:** Swollen joint count on the 28-joint (DAS28) scale (integer 0-28; component of the DAS28 composite). Baseline value is typical; document time-varying use in per-model `notes`.
- **Units:** count (0-28)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used as a shifted power term `((SWOL_28JOINT + 1)/(<ref> + 1))^exponent` to avoid the zero-count edge case. Reference value observed: 16 in Li 2019 (approximate dataset median of the popPK cohort).
- **Source aliases:**
  - `SWOL` — used in `Li_2019_abatacept.R` (Li 2019 Methods abbreviation).
- **Example models:** `Li_2019_abatacept.R` (power effect on CL with exponent 0.0965; not clinically relevant per Li 2019).
- **Notes:** The `_28JOINT` suffix distinguishes this from the 66/68-joint swollen count used in some earlier RA scales — register a separate canonical (`SWOL_66JOINT` or similar) if a future paper uses a different joint-count scale. Canonical name drops the `BL` prefix to match the `EASI` / `AGE` / `WT` / `ALB` convention where baseline-vs-time-varying status is documented in `covariateData` notes rather than the column name.

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

### STEROID
- **Description:** 1 = patient on systemic corticosteroid therapy at baseline (typically continued as concomitant medication during the study), 0 = no baseline corticosteroid use.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no baseline corticosteroid use).
- **Source aliases:**
  - `BSTEROID` — used in `Narwal_2013_sifalimumab.R`.
- **Example models:** `Narwal_2013_sifalimumab.R` (multiplicative on CL: `CL * (1 + 0.195 * STEROID)`).
- **Notes:** Distinct from `PRICORT`, which is strictly a prior (pre-study) indicator. `STEROID` captures concurrent corticosteroid use at / from study baseline. When a future paper needs the two jointly, both can coexist on the same subject.
### STEROID_BL (**canonical for baseline/concomitant systemic steroid use**)
- **Description:** 1 = patient is on systemic corticosteroid treatment at study entry (baseline concomitant use), 0 = not on steroids at baseline. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (not on steroids at baseline).
- **Source aliases:**
  - `BSTEROID` — used in `Zheng_2016_sifalimumab.R`.
- **Example models:** `Zheng_2016_sifalimumab.R` (multiplicative effects on CL `(1 + 0.11 * STEROID_BL)` and on V1 `(1 - 0.09 * STEROID_BL)` in the SLE phase IIb cohort, which was ~85% steroid-treated at baseline).
- **Notes:** Distinct from `PRICORT`, which captures systemic corticosteroid use *prior* to study entry as a one-time clinical-history indicator. `STEROID_BL` captures ongoing/concomitant steroid therapy at baseline in diseases where background steroid use is standard of care (SLE, severe asthma, etc.). In some studies the prior-vs-baseline distinction is immaterial (patients on steroids at baseline have also taken them prior) but the semantic difference is preserved here so that future models can use the appropriate indicator.

## Immunogenicity

### ADA_POS (**canonical**)
- **Description:** 1 = antidrug-antibody-positive, 0 = ADA-negative.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (ADA-negative).
- **Source aliases:**
  - `ADA` (semantically "ever positive") — used in `Zhu_2017_lebrikizumab.R`. When translating from a paper that uses `ADA` as "ever positive," verify the time-frame matches ADA_POS semantics before renaming.
  - `ADA` (time-varying positivity, primary covariate in Xu 2019) — used in `Xu_2019_sarilumab.R`.
- **Example models:** `Clegg_2024_nirsevimab.R`, `Hu_2026_clesrovimab.R`, `Xu_2019_sarilumab.R`.

### ADA_TITER (**canonical for continuous antidrug-antibody titer/titre**)
- **Description:** Continuous antidrug-antibody titer/titre (time-varying; matched in time to the PK sample). Covers both the British-spelling reciprocal-dilution convention (`ADA_TITRE`, with `ADA_TITRE = 1` for ADA-negative so `log_e(1) = 0` cancels a log-linear effect) and the American-spelling linear-titer convention (`ADA_TITER`, with `ADA_TITER = 0` for ADA-negative). The per-model `covariateData[[ADA_TITER]]$description` and `notes` must state which zero-encoding convention is in force so the covariate column cannot be misinterpreted.
- **Units:** Reciprocal dilution (e.g., 10, 20, 40, …, 2560) OR assay units (log2 or arbitrary) — document per-model in `covariateData[[ADA_TITER]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — ADA-negative encoded per-model (see zero-encoding note).
- **Source aliases:**
  - `ADA_TITRE` — British spelling (reciprocal-dilution convention; `1` for negative).
  - `ADA titre` — British spelling long form.
  - `ADAT` — used in `Moein_2022_etrolizumab.R` (American linear-titer convention; `0` for negative).
- **Example models:** `Jackson_2022_ixekizumab.R` (reciprocal-dilution reference convention with `ADA_TITER = 1` for negatives and `(1 + coef * log_e(ADA_TITER))` on CL), `Moein_2022_etrolizumab.R` (linear-titer convention with `ADA_TITER = 0` for negatives and `exp(theta * ADA_TITER)` on CL, per-unit-titer theta = 0.0365), `Robbie_2012_palivizumab.R` (reciprocal-dilution values 0/10/20/40/≥80 with category-specific multiplicative effects per titer bin; 0 = ADA negative reference).
- **Notes:** The prior separate `ADA_TITRE` (British, `1` = negative) and `ADA_TITER` (American, `0` = negative) canonicals were merged on 2026-04-20 into a single general-scope `ADA_TITER`. The zero-encoding convention is the load-bearing semantic and must be documented per-model. Distinct from `ADA_POS` (binary presence/absence); when the paper reports both, the final model usually keeps only one. Imputation rules (LOCF / NOCB / baseline-as-negative) should be documented per-model.

## Disease / treatment history

### PRIOR_TNF (**canonical**)
- **Description:** 1 = subject previously treated with an anti-TNF (tumor necrosis factor) inhibitor, 0 = TNF-naive.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (TNF-naive).
- **Source aliases:**
  - `PRIORTNF` (all caps, no underscore) — acceptable alternative spelling.
- **Example models:** `Moein_2022_etrolizumab.R` (multiplicative fractional effect on CL, +4.9%).
- **Notes:** Use when the source paper reports a binary "prior anti-TNF inhibitor" covariate on any PK parameter. Generally applicable across RA/PsA/IBD/axSpA biologic PK models.

### DISEXT_EP (**canonical for extensive colitis / pancolitis indicator**)
- **Description:** 1 = extensive colitis or pancolitis disease extension, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (left-sided colitis, when paired with `DISEXT_OTHER = 0`).
- **Source aliases:**
  - Derived from a multi-level `DISEXT` column in the source (levels: left-sided colitis, extensive/pancolitis, other): `DISEXT_EP = as.integer(DISEXT == "extensive/pancolitis")`.
- **Example models:** `Moein_2022_etrolizumab.R` (multiplicative effect on CL, +8.2% vs. left-sided colitis).
- **Notes:** Paired with `DISEXT_OTHER`. `DISEXT_EP = DISEXT_OTHER = 0` corresponds to the left-sided colitis reference group. UC-specific covariate; analogous indicators could be registered for other inflammatory-bowel-disease models by promoting this entry or adding new `DISEXT_<CATEGORY>` canonicals.

### DISEXT_OTHER (**canonical for 'other disease extension' indicator**)
- **Description:** 1 = disease extension other than left-sided colitis or extensive/pancolitis, 0 = not.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (left-sided colitis, when paired with `DISEXT_EP = 0`).
- **Source aliases:** Derived from a multi-level `DISEXT` column: `DISEXT_OTHER = as.integer(DISEXT == "other")`.
- **Example models:** `Moein_2022_etrolizumab.R` (multiplicative effect on CL, +18% vs. left-sided colitis; large uncertainty due to 2% prevalence).
- **Notes:** Paired with `DISEXT_EP`; together they encode the three-level disease-extension categorical.

## Lifestyle / medical history

### SMOKE
- **Description:** 1 = current smoker at baseline, 0 = non-smoker.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (non-smoker).
- **Source aliases:**
  - `Smoking` (case-insensitive) — used in `Ma_2020_sarilumab_anc.R`.
- **Example models:** `Ma_2020_sarilumab_anc.R` (power-form on baseline ANC: `BASE * 1.15^SMOKE`).
- **Notes:** Baseline-only indicator; does not track within-study smoking-cessation changes.

## Formulation / assay / study

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
- **Notes:** Disease-backbone indicator. If a future model needs more backbone categories, encode each as its own indicator (`COMB_CAPOX`, `COMB_FOLFOX`, …) with a single reference group.

### FORM_DP2
- **Description:** 1 = sarilumab drug product 2 formulation (used in some phase I studies and the dose-ranging phase II study), 0 = other drug product (DP1 or DP3; DP3 is the commercial formulation).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (DP1 or DP3).
- **Source aliases:**
  - `DP2` — used in `Xu_2019_sarilumab.R`.
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
- **Reference category:** n/a — used with power scaling `(COHDOSE / ref)^exponent`. Reference value observed: 1 mg/kg in Narwal 2013.
- **Source aliases:**
  - `DOSE` — used in `Narwal_2013_sifalimumab.R` (the paper's Eq. 3 variable name; renamed to `COHDOSE` here to avoid colliding with the rxode2/nlmixr2 event-column convention where `DOSE` or `AMT` carries the administered dose).
- **Example models:** `Narwal_2013_sifalimumab.R` (reference 1 mg/kg, exponent 0.0542 on CL).
- **Notes:** Scope: specific because the interpretation depends on a study design where each subject stays on a single dose cohort. For fixed-dose simulations, set `COHDOSE = nominal_dose_mg / WT` per subject. When the subject receives a weight-based dose, `COHDOSE` is the mg/kg label (0.3, 1, 3, or 10 mg/kg for the MI-CP152 cohorts).
### DOSE (**canonical for per-subject assigned dose level**)
- **Description:** Continuous covariate carrying each subject's assigned fixed dose level (in mg) across the study. Used as a power-style covariate when the population PK model detects a dose-dependent shift in a PK parameter (central volume, clearance, etc.).
- **Units:** mg (document per-model in `covariateData[[DOSE]]$units` if a different dose unit is used).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used with power scaling `(DOSE / ref)^exponent`. Reference values observed: 600 mg in `Zheng_2016_sifalimumab.R` (middle of the 200/600/1200 mg phase IIb dose range).
- **Source aliases:**
  - `Dose` — used in `Zheng_2016_sifalimumab.R`.
- **Example models:** `Zheng_2016_sifalimumab.R` (power effect on V1 with exponent 0.06).
- **Notes:** Distinct from `DOSE_70MG` (binary indicator for a specific dose group in a trinary-dose design). Use `DOSE` when the source paper treats dose as a continuous covariate acting on a PK parameter. Scope: specific because the semantics (units, reference dose, whether the covariate applies to CL or V) are study-specific; future models using continuous dose as a covariate should either extend this entry's example-models list or register a more specific variant.

### DOSE_70MG
- **Description:** 1 = subject is on the 70 mg SC Q4W dose regimen, 0 = subject is on the 210 or 490 mg SC Q4W regimen.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (210 mg or 490 mg Q4W regimen).
- **Source aliases:** derived per subject from the trial-assigned dose level.
- **Example models:** `Kotani_2022_astegolimab.R`.
- **Notes:** Zenyatta-study categorical covariate flagging the 70 mg group (lowest dose), modeled as a −15.3% relative change on relative bioavailability. Modeled by Kotani 2022 as `70 mg vs {210 mg, 490 mg}` combined reference.

### STUDY1
- **Description:** 1 = subject enrolled in Study 1 of the Cirincione 2017 pooled analysis, 0 = other. Used to switch the residual-error magnitude per study.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (all other studies; combined with `STUDY5 = 0` selects the pooled "other" residual error).
- **Source aliases:**
  - `DVID = "study1"` (character-valued study identifier; `STUDY1 = as.integer(DVID == "study1")`) — legacy form previously used in `Cirincione_2017_exenatide.R`.
- **Example models:** `Cirincione_2017_exenatide.R`.
- **Notes:** Paired with `STUDY5`. When both are 0, the subject is in the pooled "other studies" residual-error group.

### STUDY5
- **Description:** 1 = subject enrolled in Study 5 of the Cirincione 2017 pooled analysis, 0 = other. Used to switch the residual-error magnitude per study.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (all other studies).
- **Source aliases:**
  - `DVID = "study5"` (character-valued study identifier; `STUDY5 = as.integer(DVID == "study5")`) — legacy form previously used in `Cirincione_2017_exenatide.R`.
- **Example models:** `Cirincione_2017_exenatide.R`.
- **Notes:** Paired with `STUDY1`. When both are 0, the subject is in the pooled "other studies" residual-error group.

### PHASE2
- **Description:** 1 = subject enrolled in the Phase II study (MORAb-003-002) of the Farrell 2012 pooled analysis; 0 = Phase I study (MORAb-003-001). Used to switch the residual-error magnitude per study.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (Phase I).
- **Source aliases:** derived per subject from the trial identifier (`MORAb-003-001` → 0, `MORAb-003-002` → 1).
- **Example models:** `Farrell_2012_farletuzumab.R`.
- **Notes:** Farrell 2012 Table 3 reports separate residual-error estimates for the two studies — Phase I uses a proportional-only model (σ = 20.5%); Phase II uses a combined additive + proportional model (σ_prop = 34.9%, σ_add = 7.94 µg/mL). The `PHASE2` indicator selects between them.

### PHASE1
- **Description:** 1 = subject enrolled in a Phase 1 study of the Valenzuela 2025 pooled analysis (MOM-M281-001, MOM-M281-007, MOM-M281-010, EDI1001, EDI1002 — healthy participants); 0 = Phase 2 study (MOM-M281-004 / Vivacity-MG — participants with gMG). Used to switch the proportional PK residual-error magnitude per study phase.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (Phase 2).
- **Source aliases:** derived per subject from the trial identifier (Phase 1 protocols → 1, `NCT03772587` Vivacity-MG → 0).
- **Example models:** `Valenzuela_2025_nipocalimab.R`.
- **Notes:** Valenzuela 2025 Table 3 reports proportional PK residual 0.0834 (Phase 1) vs 0.367 (Phase 2). Distinct from Farrell 2012 `PHASE2` — the reference category is inverted (Valenzuela 2025 picks Phase 1 as the 1-level).

### ELISA
- **Description:** 1 = serum nipocalimab concentration measured by ELISA assay (LLOQ 0.150 µg/mL; studies MOM-M281-001, MOM-M281-007, MOM-M281-010); 0 = measured by ECLIA assay (LLOQ 0.010 µg/mL; studies EDI1001, EDI1002, MOM-M281-004). Used to switch the additive PK residual-error magnitude per assay.
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
- **Source aliases:** derived per subject from the trial identifier (`NCT03772587` → 1, else → 0).
- **Example models:** `Valenzuela_2025_nipocalimab.R`.
- **Notes:** Valenzuela 2025 estimated a slightly lower IgG baseline in MOM-M281-004 participants (baseline scaled by `FRIgG0_M281_004 = 0.777` vs. 1 in other studies). Distinct from the disease-state indicator implied by `gMG` — it is specifically the Vivacity-MG study flag because the IgG baseline factor was only estimated for that study.

## Occasion / period (IOV)

### ooc1, ooc2, ooc3, ooc4
- **Description:** Mutually exclusive occasion indicators for a crossover / multi-period design. Exactly one is 1 per observation.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Example models:** `Xie_2019_agomelatine.R`.
- **Notes:** Lower case preserved from source file. Future models should standardize on `OCC` as a categorical column with integer values (`OCC = 1`, `2`, …) and use `IOV` on the appropriate parameters.

---

## Change log

- **2026-04-21** — Added `RACE_HISPANIC` (general) and `CLD_PREM` (general) canonical entries while extracting `Robbie_2012_palivizumab.R`. Extended `ADA_TITER` example_models with the Robbie 2012 category-by-titer-bin usage, and added `Robbie_2012_palivizumab.R` to the `WT`, `PAGE`, `RACE_BLACK`, `RACE_ASIAN`, and `RACE_OTHER` example lists. `HISPANIC`, `CLD`, and `BPD` recorded as source aliases.
- **Initial seed**: Every covariate observed in `inst/modeldb/` as of the audit. Canonical names established: `SEXF`, `ADA_POS`, `RACE_<GROUP>` prefix. Aliases documented but existing model files not modified.
- **2026-04-19** — Added `CREAT`, `hsCRP`, `ALB` canonical entries after the
  Phase 6 pilot extracted Fasanmade 2009 infliximab and Thakre 2022
  risankizumab. `ALB` had been used informally in two models; now ratified
  with per-model unit documentation (g/dL vs g/L). `CREAT` chosen over
  `CRE`/`SCR`; `hsCRP` preserves lowercase `hs` prefix per the `eGFR`
  precedent. See `tracking/decision_log.md` in the mab_human_consensus
  project for the deliberation.
- **2026-04-19** — Added `BSA`, `HGB`, `TBILI`, `GAST`, and `COMB_EOX`
  canonical entries in support of Yamada 2025 zolbetuximab. New sections
  created for `Hematology` and `Surgical history / disease state`.
- **2026-04-19** — Added `BEOS` (baseline blood eosinophil count, cells/µL)
  and `DOSE_70MG` (Zenyatta 70 mg dose-group indicator) canonical entries
  after extracting the Kotani 2022 astegolimab population PK model.
- **2026-04-19** — Added `ADA_TITER`, `PRIOR_TNF`, `DISEXT_EP`,
  `DISEXT_OTHER` canonical entries while extracting Moein 2022 etrolizumab.
- **2026-04-20** — Added `CRP` and `BMI` canonical entries for Chua 2025
  mirikizumab.
- **2026-04-20** — Added `SMOKE` canonical entry from the Ma 2020 sarilumab
  ANC PopPK/PD extraction.
- **2026-04-20** — Added `EASI` canonical entry for the Eczema Area and
  Severity Index, introduced by `Tiraboschi_2025_amlitelimab.R`.
- **2026-04-20** — Added `TUMSZ`, `TUMTP_CHL`, `TUMTP_GC` canonical entries
  with the Budha 2023 tislelizumab extraction.
- **2026-04-20** — Added `ADA_TITRE` canonical entry while extracting the
  Jackson 2022 ixekizumab paediatric psoriasis PopPK model.
- **Xu 2019 sarilumab**: Added canonical entries `ALBR` (albumin / ULN ratio), `CRCL_BSA` (BSA-normalized creatinine clearance), `BLCRP` (baseline C-reactive protein), and `FORM_DP2` (sarilumab drug product 2 indicator). Extended the `ADA_POS` alias list to include the time-varying `ADA` column used in Xu 2019.
- **Ma 2020 sarilumab DAS28-CRP**: Added canonical entries `BLPHYVAS`, `BLHAQ`, and `PRICORT`. Extended the `BLCRP` entry to record Ma 2020 as a second example model.
- **2026-04-21** — Added `PHASE2` canonical entry (specific scope) for the Farrell 2012 farletuzumab per-study residual-error switch, analogous to `STUDY1`/`STUDY5` in Cirincione 2017.
- **2026-04-20 (scope + mergers)** — Introduced the `Scope:` field (general vs specific) on every entry; `checkModelConventions()` now warns when a scope-specific covariate appears in a model that is not in its `Example models` list. Merged `hsCRP` + `BLCRP` + the previously separate standard-assay `CRP` into a single general-scope `CRP` canonical (assay and baseline-vs-time-varying nuances now live in per-model `covariateData[[CRP]]$description` / `notes`). Merged `eGFR` + `CRCL_BSA` into a single general-scope `CRCL` canonical (MDRD-estimated vs measured-CrCl-BSA-normalized nuances now live per-model). Merged `ADA_TITRE` + `ADA_TITER` into a single general-scope `ADA_TITER` canonical (British-reciprocal-dilution vs American-linear-titer zero-encoding conventions now live per-model). Renamed `BEOS` → `EOS` and `GAST` → `PRIOR_GAST` to follow the `EASI` / `AGE` / `ALB` (baseline status in notes, not column name) and `PRIOR_TNF` / `PRICORT` (prior-treatment / surgical-history) naming patterns. Promoted `TUMSZ` from Budha-specific to general oncology.
- **2026-04-21** — Added `WBC` canonical entry (total white blood cell
  count) under `Hematology` while extracting `Mould_2007_alemtuzumab.R`.
  Scope: general; reference 10 × 10⁹/L used for the Vmax power covariate
  effect in CLL.
- **2026-04-21** — Added `MGADL` (Myasthenia Gravis Activities of Daily
  Living score; general scope; baseline-only in Valenzuela 2025),
  `ELISA` / `PHASE1` / `STUDY_M281_004` (all specific scope; assay and
  study flags for Valenzuela 2025 nipocalimab) while extracting
  `Valenzuela_2025_nipocalimab.R`.
- **2026-04-21** — Added `RACE_JAPANESE` canonical entry (scope: general)
  while extracting `Wade_2015_certolizumab.R`, where Japanese is broken
  out as a distinct category from the broader `RACE_ASIAN` indicator
  (both non-overlapping sub-indicators point to different per-population
  reference groups).
- **2026-04-21** — Added `STEROID` (general-scope; baseline systemic corticosteroid use, distinct from the existing `PRICORT` which captures prior corticosteroid use), `BGENE21` (specific-scope; baseline 21-gene type I interferon signature under a new `Interferon / biomarker panels` section), and `COHDOSE` (specific-scope; randomized dose cohort in mg/kg) canonical entries while extracting Narwal 2013 sifalimumab. Source aliases mapped: `BSTEROID`→`STEROID`, `DOSE`→`COHDOSE`.
- **2026-04-21** — Added `AST` (general, hepatic enzyme), `HER2_ECD`
  (specific, HER2 shed ECD biomarker), and `TRAST_BL` (specific,
  baseline trastuzumab concentration from prior therapy) canonical
  entries while extracting `Lu_2014_trastuzumabemtansine.R`.
- **2026-04-21** — Added `DIS_UC` (ulcerative colitis disease-state
  indicator) and `DIS_SASTHMA` (moderate-to-severe asthma disease-state
  indicator) canonical entries under a new `Disease state
  (cross-population indicators)` section while extracting
  `Hua_2015_anrukinzumab.R`. Both scope: specific, decomposed from a
  four-level disease-state categorical (healthy volunteer, mild-to-moderate
  asthma, moderate-to-severe asthma, UC).
- **2026-04-21** — Added IBD canonical entries for Rosario 2015 vedolizumab
  extraction: `CALPRO` (fecal calprotectin, general), `CDAI` (Crohn's
  Disease Activity Index, general), `PMAYO` (partial Mayo score, general),
  `IBD_CD` (Crohn's-disease-vs-UC indicator, general), and four
  concomitant-medication indicators `CONMED_AZA`, `CONMED_MP`,
  `CONMED_MTX`, `CONMED_AMINO` (all general scope). New H2 sections
  `Inflammatory-bowel-disease disease-activity covariates`,
  `Inflammatory-bowel-disease diagnosis`, and `Concomitant IBD
  medications`.
- **2026-04-21** — Added `STEROID_BL` (general-scope baseline/concomitant
  steroid-use indicator; distinct from prior-only `PRICORT`), `BGENE21`
  (specific-scope 21-gene type I interferon signature score), and `DOSE`
  (specific-scope per-subject assigned dose level in mg) canonical entries
  while extracting `Zheng_2016_sifalimumab.R`.
- **2026-04-24** — Added `CONMED_NSAID` (general-scope concomitant NSAID
  indicator; extends the `CONMED_*` IBD pattern) and `SWOL_28JOINT`
  (general-scope 28-joint swollen joint count; RA disease-activity
  section) canonical entries while extracting `Li_2019_abatacept.R`.
  The abatacept-specific SC phase-2 formulation indicator
  (`FORM_ABA_PHASE2`) was kept model-specific per the nlmixr2lib
  global policy that `FORM_*` covariates are not promoted to canonical
  unless they clearly generalize across multiple drugs.
- **2026-04-24** — Added `BCVA` canonical entry (best-corrected visual acuity
  score in ETDRS letters; scope: specific; ophthalmology-specific baseline
  input to indirect-response BCVA PD models) while extracting
  `Mulyukov_2018_ranibizumab.R`. Source alias `BVA` mapped. Reference value
  55 letters (study-population mean baseline BCVA).
- **2026-04-24** — Added `ECOG_GE1` (general-scope Eastern Cooperative Oncology Group performance-status indicator; 1 if ECOG >= 1, reference 0) canonical entry under `Oncology` while extracting `Bajaj_2017_nivolumab.R`. Source alias `PS` / `BPS` mapped to `ECOG_GE1`; Bajaj 2017's ECOG derivation (KPS-to-ECOG crosswalk for one study per Oken 1982) documented in the model's `covariateData` notes, not in the register. Added `Bajaj_2017_nivolumab.R` to the `WT`, `SEXF`, `RACE_ASIAN`, and `CRCL` example-model lists.
- Subsequent additions: append new canonical entries as new papers are processed. When adding, bump the audit-completed count in the summary below.

## Summary

- Files audited: 69 R files under `inst/modeldb/` (20 of which reference covariates).
- Canonical H3 entries: 58 (61 parsed entries when counting each `ooc<n>` individually; was ~57 before the 2026-04-20 mergers of `hsCRP`+`BLCRP`+standard-`CRP` → `CRP`, `eGFR`+`CRCL_BSA` → `CRCL`, and `ADA_TITRE`+`ADA_TITER` → `ADA_TITER`; +1 `PHASE2` and +1 `WBC` on 2026-04-21 from Farrell 2012 / Mould 2007; +3 on 2026-04-21 for `STEROID`, `BGENE21`, `COHDOSE` from Narwal 2013).
- Scope: general: 35. Scope: specific: 26 (counting each `ooc<n>` individually, or 23 counting the `### ooc1, ooc2, ooc3, ooc4` heading as one entry).
- Aliases mapped (selected): SEXM→SEXF, ADA→ADA_POS, ADA_TITRE/ADAT→ADA_TITER, BLACK→RACE_BLACK, ASIAN→RACE_ASIAN, MULTIRACIAL→RACE_MULTI, BLACK_OTH→RACE_BLACK_OTH, ASIAN_AMIND_MULTI→RACE_ASIAN_AMIND_MULTI, DVID→STUDY1/STUDY5, CRE→CREAT, hsCRP/HSCRP/CRPHS/BLCRP→CRP, eGFR/EGFR/CRCL_BSA/1.73*CrCl/BSA→CRCL, DP2→FORM_DP2, DISEXT→DISEXT_EP/DISEXT_OTHER, BEASI→EASI, BEOS→EOS, GAST→PRIOR_GAST, COMB→COMB_EOX, TUMTP→TUMTP_CHL/TUMTP_GC, BSTEROID→STEROID, DOSE→COHDOSE.
- Canonical H3 entries: 53 (56 parsed entries when counting each `ooc<n>` individually; was ~57 before the 2026-04-20 mergers of `hsCRP`+`BLCRP`+standard-`CRP` → `CRP`, `eGFR`+`CRCL_BSA` → `CRCL`, and `ADA_TITRE`+`ADA_TITER` → `ADA_TITER`).
- Scope: general: 33. Scope: specific: 23 (counting each `ooc<n>` individually, or 20 counting the `### ooc1, ooc2, ooc3, ooc4` heading as one entry).
- Aliases mapped (selected): SEXM→SEXF, ADA→ADA_POS, ADA_TITRE/ADAT→ADA_TITER, BLACK→RACE_BLACK, ASIAN→RACE_ASIAN, MULTIRACIAL→RACE_MULTI, BLACK_OTH→RACE_BLACK_OTH, ASIAN_AMIND_MULTI→RACE_ASIAN_AMIND_MULTI, DVID→STUDY1/STUDY5, CRE→CREAT, hsCRP/HSCRP/CRPHS/BLCRP→CRP, eGFR/EGFR/CRCL_BSA/1.73*CrCl/BSA→CRCL, DP2→FORM_DP2, DISEXT→DISEXT_EP/DISEXT_OTHER, BEASI→EASI, BEOS→EOS, GAST→PRIOR_GAST, COMB→COMB_EOX, TUMTP→TUMTP_CHL/TUMTP_GC, DX→IBD_CD, AZA→CONMED_AZA, MP→CONMED_MP, MTX→CONMED_MTX, AMINO→CONMED_AMINO.
