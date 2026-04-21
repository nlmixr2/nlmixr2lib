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
- **Example models:** `Clegg_2024_nirsevimab.R`, `Hu_2026_clesrovimab.R`, `Zhu_2017_lebrikizumab.R`, `Kovalenko_2020_dupilumab.R`, `CarlssonPetri_2021_liraglutide.R`, `Cirincione_2017_exenatide.R`, `Grimm_2023_gantenerumab.R`, `Grimm_2023_trontinemab.R`, `Kyhl_2016_nalmefene.R`, `Soehoel_2022_tralokinumab.R`, `Xie_2019_agomelatine.R`, `PK_2cmt_mAb_Davda_2014.R`, `phenylalanine_charbonneau_2021.R`, `Chua_2025_mirikizumab.R`, `Jackson_2022_ixekizumab.R`, `Kotani_2022_astegolimab.R`, `Ma_2020_sarilumab_anc.R`, `Ma_2020_sarilumab_das28crp.R`, `Moein_2022_etrolizumab.R`, `Tiraboschi_2025_amlitelimab.R`.
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
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical), `CarlssonPetri_2021_liraglutide.R` (alias `SEXM`).
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
- **Example models:** `Clegg_2024_nirsevimab.R`.

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
- **Reference category:** n/a — used with power scaling `(CRCL / ref)^exponent`. Reference values observed: 80 mL/min/1.73 m² (Cirincione 2017, MDRD eGFR), 100 mL/min/1.73 m² (Xu 2019, measured-CrCl BSA-normalized).
- **Source aliases:**
  - `eGFR` — MDRD-estimated glomerular filtration rate; used in `Cirincione_2017_exenatide.R` and `Kotani_2022_astegolimab.R`.
  - `EGFR` — all-caps variant.
  - `CRCL_BSA` — BSA-normalized creatinine clearance (measured CrCl ÷ BSA × 1.73); used in `Xu_2019_sarilumab.R`.
  - `1.73*CrCl/BSA` — the formula form appearing in Xu 2019 Eq. for Vm.
- **Example models:** `Cirincione_2017_exenatide.R` (MDRD eGFR), `Xu_2019_sarilumab.R` (measured CrCl BSA-normalized), `Kotani_2022_astegolimab.R` (MDRD eGFR).
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
- **Example models:** `Fasanmade_2009_infliximab.R` (g/dL, reference 4.1), `Thakre_2022_risankizumab.R` (g/L, reference 45), `Chua_2025_mirikizumab.R`, `Moein_2022_etrolizumab.R`, `Tiraboschi_2025_amlitelimab.R`, `Yamada_2025_zolbetuximab.R`.
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
  - `BLACK` — used in `Hu_2026_clesrovimab.R`.
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical form).

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
- **Source aliases:** `ASIAN` — used in `Hu_2026_clesrovimab.R`.
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical form).

### RACE_ASIAN_AMIND_MULTI (**canonical for composite group**)
- **Description:** 1 = Asian, American Indian / Alaskan Native, or Multiple races, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = White, Black / African American, Native Hawaiian / Pacific Islander, or Other (Clegg 2024 grouping).
- **Source aliases:** `ASIAN_AMIND_MULTI` — used in `Clegg_2024_nirsevimab.R`.
- **Example models:** `Clegg_2024_nirsevimab.R`.
- **Notes:** Clegg 2024 applies this covariate to both CL and V2 with different coefficients.

### RACE_JAPANESE (**canonical**)
- **Description:** 1 = Japanese ethnicity, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (document the actual reference groups used).
- **Source aliases:**
  - `RACE` level `8` (integer-coded, `IF(RACE.EQ.8) VRACE = ...`) — used in `Wade_2015_certolizumab.R` (derive `RACE_JAPANESE = as.integer(RACE == 8)`).
- **Example models:** `Wade_2015_certolizumab.R` (fractional effect on V, -25% vs. the pooled reference of White, Indian Asian, Hispanic, Other).
- **Notes:** Follows the `RACE_BLACK` / `RACE_ASIAN` indicator-decomposition pattern. Useful when a paper breaks out Japanese as a distinct category from a broader "Asian" group (common in multinational trials with dedicated Japanese cohorts). When paired with `RACE_ASIAN`, the `RACE_ASIAN` column must encode "Asian excluding Japanese" — document the exclusion in the per-model `covariateData[[RACE_ASIAN]]$notes`.

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
- **Example models:** `Zhu_2017_lebrikizumab.R`.

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

## Oncology

### TUMSZ (**canonical for baseline tumor size**)
- **Description:** Baseline tumor size. For solid tumors, the sum of diameters of target lesions per RECIST; for classical Hodgkin lymphoma, the sum of products of perpendicular diameters (SPPD).
- **Units:** mm
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(TUMSZ / ref)^exponent`. Reference values observed: 63 mm (Budha 2023).
- **Source aliases:** none.
- **Example models:** `Budha_2023_tislelizumab.R` (reference 63 mm).
- **Notes:** Promoted to scope: general on 2026-04-20 as a conventional oncology baseline-tumor-size measure (RECIST for solid tumors, SPPD for cHL). The SPPD vs sum-of-diameters convention is pooled onto a single column; document the per-model mixture where relevant.

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
- **Example models:** `Jackson_2022_ixekizumab.R` (reciprocal-dilution reference convention with `ADA_TITER = 1` for negatives and `(1 + coef * log_e(ADA_TITER))` on CL), `Moein_2022_etrolizumab.R` (linear-titer convention with `ADA_TITER = 0` for negatives and `exp(theta * ADA_TITER)` on CL, per-unit-titer theta = 0.0365).
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
- **2026-04-21** — Added `RACE_JAPANESE` canonical entry (scope: general)
  while extracting `Wade_2015_certolizumab.R`, where Japanese is broken
  out as a distinct category from the broader `RACE_ASIAN` indicator
  (both non-overlapping sub-indicators point to different per-population
  reference groups).
- Subsequent additions: append new canonical entries as new papers are processed. When adding, bump the audit-completed count in the summary below.

## Summary

- Files audited: 69 R files under `inst/modeldb/` (20 of which reference covariates).
- Canonical H3 entries: 53 (56 parsed entries when counting each `ooc<n>` individually; was ~57 before the 2026-04-20 mergers of `hsCRP`+`BLCRP`+standard-`CRP` → `CRP`, `eGFR`+`CRCL_BSA` → `CRCL`, and `ADA_TITRE`+`ADA_TITER` → `ADA_TITER`).
- Scope: general: 33. Scope: specific: 23 (counting each `ooc<n>` individually, or 20 counting the `### ooc1, ooc2, ooc3, ooc4` heading as one entry).
- Aliases mapped (selected): SEXM→SEXF, ADA→ADA_POS, ADA_TITRE/ADAT→ADA_TITER, BLACK→RACE_BLACK, ASIAN→RACE_ASIAN, MULTIRACIAL→RACE_MULTI, BLACK_OTH→RACE_BLACK_OTH, ASIAN_AMIND_MULTI→RACE_ASIAN_AMIND_MULTI, DVID→STUDY1/STUDY5, CRE→CREAT, hsCRP/HSCRP/CRPHS/BLCRP→CRP, eGFR/EGFR/CRCL_BSA/1.73*CrCl/BSA→CRCL, DP2→FORM_DP2, DISEXT→DISEXT_EP/DISEXT_OTHER, BEASI→EASI, BEOS→EOS, GAST→PRIOR_GAST, COMB→COMB_EOX, TUMTP→TUMTP_CHL/TUMTP_GC.
