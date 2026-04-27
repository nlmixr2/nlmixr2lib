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
- **Example models:** `Clegg_2024_nirsevimab.R`, `Hu_2026_clesrovimab.R`, `Zhu_2017_lebrikizumab.R`, `Kovalenko_2020_dupilumab.R`, `CarlssonPetri_2021_liraglutide.R`, `Cirincione_2017_exenatide.R`, `Grimm_2023_gantenerumab.R`, `Grimm_2023_trontinemab.R`, `Kyhl_2016_nalmefene.R`, `Soehoel_2022_tralokinumab.R`, `Xie_2019_agomelatine.R`, `PK_2cmt_mAb_Davda_2014.R`, `phenylalanine_charbonneau_2021.R`, `Chua_2025_mirikizumab.R`, `Jackson_2022_ixekizumab.R`, `Kotani_2022_astegolimab.R`, `Ma_2020_sarilumab_anc.R`, `Ma_2020_sarilumab_das28crp.R`, `Moein_2022_etrolizumab.R`, `Tiraboschi_2025_amlitelimab.R`, `Robbie_2012_palivizumab.R`, `Bajaj_2017_nivolumab.R`, `Quartino_2019_trastuzumab.R`, `Wang_2020_ontamalimab.R`.
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
- **Reference category:** n/a — used with power scaling `(FFM / ref)^exponent`. Reference values observed: 40.69 kg (Zhou 2021 belimumab pooled adult+pediatric SLE), 45 kg (Aguiar 2021, Crohn's disease cohort median).
- **Source aliases:** none; `FFM` is the universal abbreviation.
- **Example models:** `Zhou_2021_belimumab.R` (reference 40.69 kg; exponents 0.673 on CL and 0.891 on V1), `Aguiar_2021_ustekinumab.R` (reference 45 kg; power exponents 0.598 on CL, 0.590 on Vc, 0.586 on Vp).
- **Notes:** Distinct from `LBM` (lean body mass) which is sometimes computed by the Boer or Hume formulae. When the source paper reports the body-composition formula it used (e.g., Janmahasatian for FFM), record it in `covariateData[[FFM]]$notes`. FFM is preferred over total body weight when scaling monoclonal-antibody PK because mAb distribution is largely confined to extracellular fluid; muscle / lean tissue tracks extracellular volume better than total weight in heavier patients.

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
- **Source aliases:**
  - `BALB` (baseline albumin) — used in `Zhou_2021_belimumab.R`. Maps directly to `ALB`; baseline-vs-time-varying status documented in per-model notes.
- **Example models:** `Fasanmade_2009_infliximab.R` (g/dL, reference 4.1), `Thakre_2022_risankizumab.R` (g/L, reference 45), `Chua_2025_mirikizumab.R`, `Moein_2022_etrolizumab.R`, `Tiraboschi_2025_amlitelimab.R`, `Yamada_2025_zolbetuximab.R`, `Li_2019_abatacept.R` (g/dL, reference 4.0; the Li 2019 Methods states 'mg/dL' which is a publication typo — see the model's `covariateData[[ALB]]$notes`), `Quartino_2019_trastuzumab.R` (g/dL, reference 4; source column `ALBU`; negative exponent -0.998 on linear CL), `Wang_2020_ontamalimab.R` (g/L, reference 39), `Zhou_2021_belimumab.R` (g/L, reference 40; baseline-only, source column `BALB`).
- **Notes:** Ratified canonically on 2026-04-19 after cross-model review. Unit varies by paper (g/dL in US-convention papers, g/L in SI-convention papers); the per-model `covariateData[[ALB]]$units` field is load-bearing. Effect-coefficient magnitude is meaningless without the unit.

### IGG (**canonical for serum immunoglobulin G**)
- **Description:** Serum total immunoglobulin G concentration (baseline or time-varying). Used in mAb PK analyses as a competition-for-FcRn-recycling covariate on therapeutic-mAb clearance — high endogenous IgG is hypothesized to displace the therapeutic mAb from FcRn salvage and increase its catabolic clearance.
- **Units:** g/L (typical in SI-convention papers); also reported as mg/dL in US-convention papers — document the unit used in each model via `covariateData[[IGG]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(IGG / ref)^exponent`. Reference values observed: 14.8 g/L (Zhou 2021), 9.65 g/L (Yang 2021).
- **Source aliases:**
  - `BIGG` (baseline IgG) — used in `Zhou_2021_belimumab.R`.
  - `IGGBL` (baseline IgG) — used in `Yang_2021_cemiplimab.R`.
- **Example models:** `Zhou_2021_belimumab.R` (g/L, reference 14.8; baseline-only; exponent 0.293 on CL), `Yang_2021_cemiplimab.R` (g/L, reference 9.65; small positive exponent 0.184 on shared CL/Q).
- **Notes:** Mechanistically meaningful for monoclonal-antibody PK because endogenous IgG competes with the therapeutic mAb for FcRn-mediated recycling. The per-model `covariateData[[IGG]]$units` field is load-bearing (1 g/L ≈ 100 mg/dL). Baseline-vs-time-varying status documented in `covariateData[[IGG]]$notes`. Distinct from `lIgG0` / IgG-as-a-state in mechanistic FcRn-competition TMDD models (e.g., `Valenzuela_2025_nipocalimab.R`), where IgG is a dynamic state, not a baseline covariate; use `IGG` only when the source paper treats IgG as a static (baseline) covariate column.

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
- **Source aliases:**
  - `SGOT` (serum glutamic-oxaloacetic transaminase; the legacy clinical-chemistry name for AST) — used in `Quartino_2019_trastuzumab.R`.
- **Example models:** `Lu_2014_trastuzumabemtansine.R` (U/L, reference 27; small positive exponent 0.071 on CL), `Quartino_2019_trastuzumab.R` (IU/L, reference 24; source column `SGOT`; positive exponent 0.205 on linear CL).
- **Notes:** Hepatic-function marker. Commonly reported alongside `ALT` and `TBILI`; register a separate `ALT` canonical if a future paper requires it. `SGOT` is the older lab-reporting name; values and units are identical to `AST`.

### ALT (**canonical for alanine aminotransferase**)
- **Description:** Serum alanine aminotransferase activity (baseline or time-varying).
- **Units:** U/L (IU/L; the two labels are used interchangeably in the clinical-PK literature). Document per-model via `covariateData[[ALT]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(ALT / ref)^exponent`.
- **Source aliases:** none; `ALT` is the universal clinical-PK abbreviation.
- **Example models:** `Nikanjam_2019_siltuximab.R` (U/L, reference 19; small negative exponent -0.096 on CL), `Melhem_2022_dostarlimab.R` (U/L, reference 18; small negative exponent -0.0585 on CL, time-varying).
- **Notes:** Hepatic-function marker. Commonly reported alongside `AST` and `TBILI`. Ratified canonically on 2026-04-24.

### LDH (**canonical for serum lactate dehydrogenase**)
- **Description:** Serum lactate dehydrogenase activity (baseline or time-varying). General-purpose marker of tissue / cellular turnover; in oncology PK analyses it is interpreted as a disease-burden / cell-turnover proxy.
- **Units:** U/L (IU/L; interchangeable). Document per-model via `covariateData[[LDH]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(LDH / ref)^exponent` or with an additive linear-on-log form `exp(coef * (log(LDH) - log(ref)))` (algebraically equivalent to `(log(LDH) / log(ref))^coef`). Reference values observed: 217 U/L (Sanghavi 2020).
- **Source aliases:**
  - `BLDH` (baseline LDH) — used in `Sanghavi_2020_ipilimumab.R`.
- **Example models:** `Sanghavi_2020_ipilimumab.R` (linear-on-log form on CL with reference 217 U/L; coefficient 0.703).
- **Notes:** Universal lab marker. Sanghavi 2020 log-transforms LDH because the distribution is heavily right-skewed (range 74-6,245 U/L over a median of 217); other papers may use a simple `(LDH/ref)^exponent` form. Document the functional form in `covariateData[[LDH]]$notes`.

### HEPIMP_MILD (**canonical for mild hepatic impairment indicator**)
- **Description:** 1 = mild hepatic impairment per the National Cancer Institute Organ Dysfunction Working Group (NCI ODWG) criteria, 0 = normal hepatic function or non-mild category. NCI ODWG mild = total bilirubin ≤ ULN with AST > ULN, OR total bilirubin > 1.0×ULN to ≤ 1.5×ULN with any AST.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (normal hepatic function; the moderate / severe categories are typically pooled into the reference for population PK analyses where mild impairment is the only category with non-trivial sample size).
- **Source aliases:**
  - `HEPIMP` (with values `1 = mild / 0 = others`) — used in `Lin_2024_casirivimab.R`.
- **Example models:** `Lin_2024_casirivimab.R` (multiplicative fractional change on CL).
- **Notes:** Use this column when a model dichotomizes hepatic-impairment status as "mild vs. others" (i.e., normal + the rare moderate/severe cases pooled into the reference). For models that test moderate or severe as separate categories, register additional canonicals `HEPIMP_MOD` / `HEPIMP_SEV` rather than overloading this entry.

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

### NLR (**canonical for neutrophil-to-lymphocyte ratio**)
- **Description:** Ratio of absolute neutrophil count to absolute lymphocyte count from a complete blood count with differential. Used as a peripheral inflammation marker. May be reported as baseline only or as a time-varying covariate.
- **Units:** ratio (unitless)
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(NLR / ref)^exponent` or exponential effects. Reference values observed: 2.11 (Lin 2024, median in pooled COVID-19 + non-infected cohort).
- **Source aliases:** none; `NLR` is the universal abbreviation in clinical-PK and inflammation-biomarker literature.
- **Example models:** `Lin_2024_casirivimab.R` (time-varying; reference 2.11; small positive exponent +0.029 on CL).
- **Notes:** Document baseline-vs-time-varying status in `covariateData[[NLR]]$notes`. Although it derives from `WBC` differential counts, register it as its own canonical because the ratio (not the absolute counts) is what the model uses.

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

### BGENE21 (**canonical for 21-gene type I interferon signature score**)
- **Description:** Baseline 21-gene type I interferon signature score — a composite transcriptomic score summarising the expression of 21 interferon-regulated genes in whole blood relative to a healthy-donor reference, used as a biomarker of type I IFN pathway activation in SLE and related autoimmune conditions.
- **Units:** unitless fold-change score (relative to healthy-donor reference).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used with power scaling `(BGENE21 / ref)^exponent`. Reference values observed: 32 in Narwal 2013 (study-population median was 33), 12.04 in Zheng 2016 (median of the SLE phase IIb cohort, range 0.32-38.59).
- **Source aliases:** none.
- **Example models:** `Narwal_2013_sifalimumab.R` (reference 32, exponent 0.0558 on CL), `Zheng_2016_sifalimumab.R` (reference 12.04, power effect on CL with exponent 0.09).
- **Notes:** Specific to drugs whose mechanism targets the type I IFN pathway (e.g., anti-IFN-alpha antibodies like sifalimumab, anifrolumab). Higher BGENE21 indicates stronger target engagement / disease activity and is associated with increased drug clearance via target-mediated mechanisms. The 21-gene panel composition is tied to the MedImmune/AstraZeneca SLE development programme; a different IFN gene signature (e.g., a 4-gene or 5-gene panel) should be registered under its own canonical name (`BGENE4`, `IFN_SIG`, ...) to avoid conflating panel definitions.

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

### BLBCELL (**canonical for baseline CD19+ B cell count**)
- **Description:** Baseline CD19+ B cell count (cells/µL) measured by fluorescence-activated cell sorting (FACS) prior to first dose. Used as a covariate / scaling biomarker for B-cell-targeted antibody PK-PD models (e.g., anti-CD20 mAbs in multiple sclerosis or B cell malignancies).
- **Units:** cells/µL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used with power scaling `(BLBCELL / ref)^exponent`. Reference value observed: 200 cells/µL (Yu 2022, median of the pooled five-study cohort).
- **Source aliases:** `Bcell0` — used in `Yu_2022_ofatumumab.R`.
- **Example models:** `Yu_2022_ofatumumab.R` (power effect on the maximum B-cell-lysis stimulatory effect Emax, exponent 0.275, reference 200 cells/µL).
- **Notes:** Distinct from a *time-varying* B cell count, which is the PD response variable rather than a covariate. Scope: specific because the clinically relevant baseline depends on the surface marker (CD19, CD20, CD22) and whether the panel reports total B cells or memory/naive subsets — register a new canonical name if a future paper uses a different marker.

### CRP (**canonical for C-reactive protein**)
- **Description:** C-reactive protein concentration. Covers both standard and high-sensitivity (hs-CRP) assays and both baseline and time-varying usages. Each model's `covariateData[[CRP]]$description` and `notes` must state the assay type (standard vs hs-CRP) and whether the column carries a baseline-only or time-varying value, including the paper-specific reference value used for power scaling.
- **Units:** mg/L (document per-model via `covariateData[[CRP]]$units`).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(CRP / ref)^exponent` or exponential effects `exp(coef * (CRP - ref))`. Reference values observed: 4.23 mg/L (Moein 2022, IBD standard assay), 4.31 mg/L (Moein 2022 Table 3 median), 5.21 mg/L (Thakre 2022, baseline hs-CRP), 7.41 mg/L (Chua 2025, baseline standard assay), 14.2 mg/L (Xu 2019, baseline standard assay), 15.7 mg/L (Ma 2020, baseline standard assay), 0.837 mg/dL = 8.37 mg/L (Wang 2020, IBD standard assay; the model carries the source unit mg/dL).
- **Source aliases:**
  - `hsCRP` — high-sensitivity CRP (mixed-case preserved from earlier register drafts).
  - `HSCRP` — all-caps variant.
  - `CRPHS` — used in `Thakre_2022_risankizumab.R` (baseline, high-sensitivity assay).
  - `BLCRP` — baseline CRP; used in `Xu_2019_sarilumab.R` and `Ma_2020_sarilumab_das28crp.R`.
- **Example models:** `Thakre_2022_risankizumab.R`, `Xu_2019_sarilumab.R`, `Chua_2025_mirikizumab.R`, `Moein_2022_etrolizumab.R`, `Ma_2020_sarilumab_das28crp.R`, `Wang_2020_ontamalimab.R` (mg/dL, reference 0.837).
- **Notes:** The prior separate `hsCRP`, `BLCRP`, and standard-assay `CRP` canonicals were merged on 2026-04-20 to a single general-scope `CRP` canonical. Assay type (standard vs hs-CRP), baseline-vs-time-varying status, and the paper-specific reference value all live in each model's `covariateData[[CRP]]$description` / `notes`. Only aggregate values from hs-validated assays as CRP when the downstream analysis relies on low-range sensitivity; for most inflammatory-disease cohorts (IBD, RA/PsA), baseline CRP is well above the hs-sensitivity range and the distinction is moot.

## Cardiometabolic / target biomarkers

### FPCSK9 (**canonical for free (unbound) proprotein convertase subtilisin/kexin type 9 concentration**)
- **Description:** Free (unbound, non-drug-bound) serum proprotein convertase subtilisin/kexin type 9 (PCSK9) concentration. For anti-PCSK9 monoclonal antibodies (alirocumab, evolocumab, bococizumab) the free-PCSK9 pool is the pharmacologically active target fraction; drug–target binding reduces FPCSK9 relative to total PCSK9.
- **Units:** ng/mL (document per-model if a paper reports a different unit).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used with linear-deviation forms `TVPARAM + theta * (FPCSK9 / ref)` or power-form `(FPCSK9 / ref)^theta`. Reference values observed: 72.9 ng/mL (Martinez 2019 time-varying median).
- **Source aliases:** none known.
- **Example models:** `Martinez_2019_alirocumab.R` (time-varying; additive-linear effect on `Km` with slope −0.541 per (FPCSK9/72.9), reference 72.9 ng/mL).
- **Notes:** Specific scope because the column is meaningful only for drugs whose mechanism is PCSK9 inhibition; reusing the name for a different anti-PCSK9 agent is acceptable (add to Example models). For non-PCSK9 drugs that use a similar target-concentration biomarker, register a new canonical (e.g., `FIL6R`, `FTNF`) rather than overloading `FPCSK9`. Per-model `covariateData[[FPCSK9]]$notes` should state whether the value is baseline-only or time-varying and how missing values were imputed (Martinez 2019 used LOCF).

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

### RACE_WHITE (**canonical**)
- **Description:** 1 = White, 0 = non-White.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (non-White; complement composition depends on the source paper, typically pooling Black/African American, Asian, American Indian/Alaska Native, Native Hawaiian/Pacific Islander, Other, Not reported, Unknown).
- **Source aliases:**
  - `RACE` (with values `1 = White / 0 = non-White`) — used in `Lin_2024_casirivimab.R`. Source column name `RACE` is generic; the canonical name is intentionally explicit because some other models use `RACE` for a different dichotomy.
- **Example models:** `Lin_2024_casirivimab.R` (multiplicative fractional change on CL relative to non-White reference).
- **Notes:** Used by papers that dichotomize race as White vs. non-White rather than decomposing into separate group indicators. Sign and reference-category interpretation are inverted relative to `RACE_BLACK` / `RACE_ASIAN` / etc.; do NOT combine `RACE_WHITE` with the decomposed indicators in the same model.

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

### RACE_NEAS (**canonical for North East Asian composite race indicator**)
- **Description:** 1 = North East Asian heritage (worldwide Chinese, Japanese, or Korean), 0 = non-North East Asian. Composite indicator analogous to `RACE_ASIAN` but specifically restricted to the East Asian subgroup most-relevant to ICH E5 ethnic-sensitivity / Asian-region bridging analyses.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (any non-North East Asian race, including South / Southeast Asian, White, Black, etc.).
- **Source aliases:**
  - `RAC4` — used in `Zhou_2021_belimumab.R` (Zhou 2021 Table 2 footnote d).
- **Example models:** `Zhou_2021_belimumab.R` (multiplicative factor 1.07 on V1).
- **Notes:** Distinct from the broader `RACE_ASIAN` (which can include South / Southeast Asian populations) because Zhou 2021 specifically tested whether Chinese/Japanese/Korean patients had different PK from the rest of the dataset; the analysis explicitly compared `RAC4` (North East Asian) against alternative race definitions and chose `RAC4` by AIC.

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

### RACE_JAPANESE (**canonical for Japanese-heritage race indicator**)
- **Description:** 1 = Japanese heritage, 0 = non-Japanese. Used when Japanese subjects form a distinct subgroup in the study design (e.g., ICH E5 bridging analyses or studies with a dedicated Japanese healthy-volunteer cohort).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (non-Japanese).
- **Source aliases:**
  - `JAPANESE_HV` — used in `Wang_2017_benralizumab.R` (Japanese healthy-volunteer cohort indicator; the healthy-volunteer vs. asthma-patient distinction is captured separately, not in this covariate).
- **Example models:** `Wade_2015_certolizumab.R` (multiplicative fractional effect on V/F; Wade 2015 breaks Japanese [RACE.EQ.8] out separately from RACE_ASIAN), `Wang_2017_benralizumab.R` (multiplicative factor 1.34 on Vc).
- **Notes:** Distinct from `RACE_NEAS` (North East Asian composite, includes Chinese, Japanese, and Korean) and from `RACE_ASIAN`. Use `RACE_JAPANESE` only when the source paper breaks out Japanese heritage as its own indicator; do not aggregate with other Asian groups when the paper keeps them separate. Ratified canonically on 2026-04-26.

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

## Comorbidities

### DIAB (**canonical for diabetes-mellitus comorbidity indicator**)
- **Description:** 1 = patient has diabetes mellitus comorbidity (Type 1 or Type 2 not distinguished), 0 = no diabetes comorbidity. Time-fixed at study entry per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no diabetes comorbidity).
- **Source aliases:**
  - `DIAB` — used in `Chen_2022_guselkumab.R`.
- **Example models:** `Chen_2022_guselkumab.R` (multiplicative effect on CL/F: 1.15^DIAB, +15% in patients with diabetes).
- **Notes:** Captures pre-existing diabetes mellitus as a comorbidity in non-diabetes-primary indications (e.g., psoriatic arthritis, psoriasis). Distinct from a primary disease-state indicator like `DIS_UC`. Type 1 vs Type 2 mellitus is not separated unless the source paper distinguishes them; in pooled-population PK analyses, the covariate is typically a single binary flag derived from medical history. Diabetic patients tend to have higher inflammation and altered IgG turnover, which can manifest as modest changes in monoclonal-antibody clearance.

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

### DIS_PJIA (**canonical for polyarticular juvenile idiopathic arthritis disease-state indicator**)
- **Description:** 1 = polyarticular juvenile idiopathic arthritis (pJIA) patient, 0 = non-pJIA (e.g., adult rheumatoid arthritis or other indication). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-pJIA subject; the complement group is defined per-model — typically adult RA in pooled abatacept analyses).
- **Source aliases:**
  - `JIA` — used in `Gandhi_2021_abatacept.R`.
- **Example models:** `Gandhi_2021_abatacept.R` (additive coefficient on logit-F: pJIA patients have markedly higher SC bioavailability than RA reference).
- **Notes:** Used when a population PK model pools pJIA patients with a non-pJIA reference population (e.g., Gandhi 2021: pooled adult RA + pediatric pJIA) and pJIA disease/age status is tested as a PK covariate (in Gandhi 2021, on bioavailability rather than CL — disease-vs-CL was not clinically relevant). Distinct from `CHILD` and `ADOLESCENT`, which are pure age-band indicators independent of indication. Scope: specific; promote to general if a second paper pools pJIA with a non-pJIA reference.

### DIS_HV (**canonical for healthy-volunteer cohort indicator**)
- **Description:** 1 = healthy volunteer (no diagnosis), 0 = patient (any diagnosis represented in the pooled cohort). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (patient subject; the complement group is the union of disease cohorts pooled in the source analysis).
- **Source aliases:** none known; healthy-volunteer indicators in source NONMEM control streams typically use ad-hoc names (e.g., `HV`, `HEALTHY`).
- **Example models:** `Nikanjam_2019_siltuximab.R` (multiplicative effects: 0.77 on CL, 0.83 on Vss; reference category is the pooled non-HV oncology cohort).
- **Notes:** Used when a population PK model pools healthy volunteers with patients across heterogeneous indications and the HV-vs-patient contrast is retained as a covariate. Scope: specific because the complement reference category is paper-defined (Nikanjam 2019 reference is "all non-HV, non-Castleman, non-SMM tumor types"). Ratified canonically on 2026-04-24.

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
- **Reference category:** 0 (non-DMD subject; the complement group is the reference cohort the source analysis pools alongside the DMD population — typically healthy adult volunteers).
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
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-MM subjects; the complement group is defined per-model — see per-model `covariateData[[MM]]$notes`).
- **Source aliases:** none; `MM` is the common abbreviation used directly in source analyses.
- **Example models:** `Ogasawara_2020_durvalumab.R` (multiplicative factor 0.820 on Vc; reference group is the union of MDS, AML, and NHL subjects).
- **Notes:** Distinct from `DIS_SMM` (smoldering / asymptomatic MM). Use `MM` when the source paper pools active MM patients with other hematologic malignancies and retains an MM indicator in the final model. Scope: specific because the reference category is paper-defined. Ratified canonically on 2026-04-26.

### MDSAML (**canonical for MDS or AML disease-type indicator**)
- **Description:** 1 = patient with myelodysplastic syndrome (MDS) or acute myeloid leukemia (AML), 0 = other hematologic malignancy or reference group. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-MDS/AML subjects; the complement group is defined per-model — typically multiple myeloma and non-Hodgkin lymphoma in Ogasawara 2020).
- **Source aliases:** none; `MDSAML` is the combined indicator used directly in source analyses.
- **Example models:** `Ogasawara_2020_durvalumab.R` (multiplicative factor 1.26 on CL; reference group is the union of MM and NHL subjects).
- **Notes:** Use `MDSAML` as a combined MDS+AML indicator when the source paper collapses the two diagnoses into one covariate. If a future paper separates MDS and AML as distinct indicators, register `DIS_MDS` and `DIS_AML` separately. Scope: specific because the reference category is paper-defined. Ratified canonically on 2026-04-26.
### DIS_BCPALL (**canonical for B-cell precursor acute lymphoblastic leukemia disease-state indicator**)
- **Description:** 1 = B-cell precursor acute lymphoblastic leukemia (BCP-ALL), 0 = B-cell non-Hodgkin's lymphoma (NHL) or other non-BCP-ALL indication pooled in the source analysis. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-BCP-ALL; in the Wu 2024 cohort this is the adult B-cell NHL stratum).
- **Source aliases:**
  - `ALL` — used in `Wu_2024_inotuzumab.R` (Wu 2024 calls it the "ALL effect" and notes it bundles disease type with the corresponding bioanalytical assay difference).
- **Example models:** `Wu_2024_inotuzumab.R` (additive fractional-change effects on CL1 (-0.767) and CL2 (-0.362), and gates the BLSTABL and AGE effects on kdes; for kdes itself a -0.924 fractional change for BCP-ALL).
- **Notes:** Used when a population PK model pools BCP-ALL patients with a non-BCP-ALL reference (e.g., Wu 2024: pooled adult B-cell NHL + adult BCP-ALL + pediatric BCP-ALL). Scope: specific because the complement reference category is paper-defined (Wu 2024 reference is pooled adult B-cell NHL). The "ALL effect" theta in Wu 2024 conflates two physiologically distinct sources of variation — B-cell tumor type (NHL vs ALL surface CD22 burden) and bioanalytical method (ELISA for adult NHL vs HPLC-MS for ALL) — and cannot be split with the available data; document this confounding when comparing across populations. Ratified canonically on 2026-04-26.
## Infectious disease (SARS-CoV-2 / COVID-19)

### SARS_VLOAD (**canonical for SARS-CoV-2 baseline viral load**)
- **Description:** Baseline (pre-treatment) SARS-CoV-2 viral load measured from nasopharyngeal swab by RT-qPCR, reported as log10 RNA copies/mL.
- **Units:** log10 copies/mL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used with power scaling `(SARS_VLOAD / ref)^exponent`. Reference values observed: 6.4 log10 copies/mL (Lin 2024, median in pooled COVID-19 cohort).
- **Source aliases:**
  - `VIRAL` — used in `Lin_2024_casirivimab.R`.
- **Example models:** `Lin_2024_casirivimab.R` (small negative exponent -0.0075 on CL).
- **Notes:** SARS-CoV-2-specific. For non-infected subjects, the value is encoded as 0 in the source dataset (below assay detection); the population-PK exponent is small enough that this 0 is absorbed by the reference shift. Register a parallel canonical for any future paper that uses a different infection (e.g., RSV, influenza).

### SARS_SEROPOS (**canonical for SARS-CoV-2 baseline serostatus positive**)
- **Description:** 1 = SARS-CoV-2 spike or nucleocapsid antibody positive at baseline (prior infection or prior vaccination), 0 = seronegative or other / unknown.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (seronegative; "Other" / unknown serostatus is typically pooled into the reference per the source paper's analysis plan).
- **Source aliases:**
  - `SERPOS` — used in `Lin_2024_casirivimab.R`.
- **Example models:** `Lin_2024_casirivimab.R` (multiplicative fractional change on CL).
- **Notes:** SARS-CoV-2-specific. The exact assay (anti-spike vs anti-nucleocapsid; vendor) varies by study; document per-model in `covariateData[[SARS_SEROPOS]]$notes`.

### OXYSUP_LOW (**canonical for low-flow supplemental oxygen indicator**)
- **Description:** 1 = subject is receiving low-flow supplemental oxygen at baseline (e.g., nasal cannula, simple face mask), 0 = no supplemental oxygen at baseline.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no supplemental oxygen at baseline; the high-flow / mechanical-ventilation categories are encoded by the parallel `OXYSUP_HIGH` indicator).
- **Source aliases:**
  - `OXYSTAT1` — used in `Lin_2024_casirivimab.R`.
- **Example models:** `Lin_2024_casirivimab.R` (multiplicative fractional change on CL; +10.6%).
- **Notes:** Decomposed indicator from a 4-level ordered categorical (no oxygen / low-flow / high-flow / mechanical ventilation). Use with the parallel `OXYSUP_HIGH` indicator. Register a separate `OXYSUP_VENT` canonical if a future analysis splits mechanical ventilation from high-flow oxygen.

### OXYSUP_HIGH (**canonical for high-flow supplemental oxygen indicator**)
- **Description:** 1 = subject is receiving high-flow supplemental oxygen at baseline (high-flow nasal cannula, non-rebreather mask, non-invasive positive-pressure ventilation, OR mechanical ventilation pooled into the high-flow category), 0 = otherwise (no supplemental oxygen or low-flow).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no supplemental oxygen at baseline; document whether the source paper pooled mechanical ventilation into this indicator or treated it separately).
- **Source aliases:**
  - `OXYSTAT2` — used in `Lin_2024_casirivimab.R`.
- **Example models:** `Lin_2024_casirivimab.R` (multiplicative fractional change on CL; +38.0%).
- **Notes:** Companion indicator to `OXYSUP_LOW`. In Lin 2024 the rare mechanical-ventilation cases were pooled into the high-flow indicator (n = 24 across the 7598-subject dataset).

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
  - `TTYPE` (Quartino 2019; categorical column with levels `MBC`, `EBC`, `HV`, `AGC`, `Others`) — decompose into `TUMTP_GC = as.integer(TTYPE == "AGC")`.
- **Example models:** `Budha_2023_tislelizumab.R`, `Quartino_2019_trastuzumab.R` (advanced gastric cancer; per-group typical-value switch on linear CL and Vc rather than an exponential multiplier).
- **Notes:** Follows the `RACE_<GROUP>` indicator-decomposition pattern. New oncology tumor types should be added as additional `TUMTP_<GROUP>` entries so the reference set stays explicit. "Advanced gastric cancer" (AGC) and "gastric cancer" (GC) are pooled onto a single `TUMTP_GC` indicator; document the per-paper stage-of-disease detail in `covariateData[[TUMTP_GC]]$notes`.

### TUMTP_OTH (**canonical for 'other tumor types' residual indicator**)
- **Description:** 1 = heterogeneous "other" tumor-type pool (typically NSCLC plus miscellaneous solid tumors such as prostate, ovarian, and colorectal), 0 = one of the named tumor-type groups in the same analysis.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = named tumor-type groups defined per-paper (e.g., MBC, EBC, HV, AGC in Quartino 2019). The complement of all `TUMTP_<GROUP>` indicators defined in the same model.
- **Source aliases:**
  - `TTYPE` (Quartino 2019) — decompose into `TUMTP_OTH = as.integer(TTYPE == "Others")`.
- **Example models:** `Quartino_2019_trastuzumab.R` (per-group typical-value switch on linear CL; NSCLC plus a small residual group of prostate, ovarian, and other histologies).
- **Notes:** Scope: specific because the set of histologies collapsed into "Others" is defined by the analysis plan of the source paper; two papers' `TUMTP_OTH` columns are not interchangeable. Document the exact per-paper composition (e.g., "NSCLC + prostate + ovarian + other, n = 107 in Quartino 2019") in `covariateData[[TUMTP_OTH]]$notes`. A given subject can have at most one of the `TUMTP_<GROUP>` indicators (including `TUMTP_OTH`) set to 1; all-zero means the reference group.

### SPDL1 (**canonical for soluble PD-L1 concentration**)
- **Description:** Baseline (or time-varying) serum concentration of soluble programmed death-ligand 1 (sPD-L1). Serves as a circulating biomarker of target burden and immune activation for anti-PD-1/PD-L1 antibodies.
- **Units:** pg/mL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used with power scaling `(SPDL1 / ref)^exponent`. Reference value observed: 173.8 pg/mL (study-population median in Ogasawara 2020).
- **Source aliases:** none; `SPDL1` is the standard abbreviation used directly in source analyses.
- **Example models:** `Ogasawara_2020_durvalumab.R` (power effect on CL, exponent 0.0617, reference 173.8 pg/mL; time-varying; values below LOD imputed as LOD/2 = 33.55 pg/mL).
- **Notes:** Scope: specific because sPD-L1 is meaningful only for drugs targeting the PD-1/PD-L1 pathway. For other checkpoint biomarkers (e.g., soluble CTLA-4, soluble LAG-3) register new dedicated canonicals rather than reusing this one. Ratified canonically on 2026-04-26.

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
  - `PS` / `BPS` — used in `Bajaj_2017_nivolumab.R` (BPS = "baseline performance status"; the one study using Karnofsky Performance Status was mapped to ECOG via Oken 1982 before binarization) and `Zhang_2019_nivolumab.R` (paper's binary collapse PS=0 vs. PS>0).
  - `ECOG_1` — alternative explicit form; equivalent to `ECOG_GE1` when ECOG only takes values 0, 1, 2 in the analysis dataset (the typical oncology case).
  - `ECOG_PS_GT0` — retired name used in earlier register drafts; semantically identical (`>= 1` equals `> 0` for integer ECOG scores).
- **Example models:** `Bajaj_2017_nivolumab.R` (exponential effect on CL with coefficient 0.172), `Zhang_2019_nivolumab.R` (exponential effect exp(0.181) on baseline CL; additive effect -0.138 on the time-varying-CL Emax parameter).
- **Notes:** Oncology papers conventionally report ECOG as an integer (0-5) but binarize at >= 1 because ECOG >= 2 is rare in trial cohorts. When a source paper provides the ordinal ECOG score separately, derive `ECOG_GE1 = as.integer(ECOG >= 1)`. Zhang 2019 uses `ECOG_GE1` on both baseline CL and the time-varying Emax parameter (unlike Bajaj 2017, which uses it on CL only); document the structural role in each model's `covariateData[[ECOG_GE1]]$notes`. If a future paper needs finer resolution (e.g., separate effects for ECOG 1 vs ECOG 2), add a parallel `ECOG_GE2` canonical rather than overloading this one.

### TUMTP_SCLC (**canonical for small-cell-lung-cancer tumor-type indicator**)
- **Description:** 1 = small cell lung cancer (SCLC), 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (e.g., melanoma, NSCLC, RCC, HCC, CRC in the Sanghavi 2020 cohort; reference category is melanoma).
- **Source aliases:**
  - `TUMTP` (categorical column with levels including `melanoma`, `NSCLC`, `SCLC`, `CRC`, `HCC`, `RCC`) — decompose into `TUMTP_SCLC = as.integer(TUMTP == "SCLC")`.
- **Example models:** `Sanghavi_2020_ipilimumab.R` (exponential coefficient -0.124 on CL).
- **Notes:** Follows the `TUMTP_CHL` / `TUMTP_GC` decomposition pattern. SCLC is the only retained tumor-type indicator in the Sanghavi 2020 final model after backward elimination; the other tumor types collapse into the reference (melanoma) group.

### LINE_1L (**canonical for first-line-therapy indicator**)
- **Description:** 1 = first-line therapy (1L), 0 = second-line or greater (2L+).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (2L+, second-line or greater).
- **Source aliases:**
  - `LINE` (categorical column with levels `1L`, `2L`, `3L+`, ...) — decompose into `LINE_1L = as.integer(LINE == "1L")`.
- **Example models:** `Sanghavi_2020_ipilimumab.R` (exponential coefficient -0.0949 on CL).
- **Notes:** Common oncology PK covariate; scope: specific because the 2L+ reference grouping (versus separate 2L vs 3L+ indicators) is tied to a particular study's design. Promote to general if a second model legitimately ratifies the same 1L vs 2L+ binarization.

### NIVO_1Q3W (**canonical for nivolumab 1 mg/kg every 3 weeks co-administration indicator**)
- **Description:** 1 = ipilimumab co-administered with nivolumab 1 mg/kg every 3 weeks; 0 = otherwise (monotherapy or any other nivolumab regimen).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no nivolumab or any non-1Q3W nivolumab regimen).
- **Source aliases:**
  - `NIVO_REGIMEN` (categorical column with levels `none`, `0.3 mg/kg Q3W`, `1 mg/kg Q2W`, `1 mg/kg Q3W`, `3 mg/kg Q2W`, `3 mg/kg Q3W`) — decompose into `NIVO_1Q3W = as.integer(NIVO_REGIMEN == "1 mg/kg Q3W")`.
- **Example models:** `Sanghavi_2020_ipilimumab.R` (exponential coefficient 0.0950 on ipilimumab CL).
- **Notes:** Paired with `NIVO_3Q2W` in the Sanghavi 2020 final model; both decomposed indicators are 0 for ipilimumab monotherapy. Other nivolumab regimens (0.3 mg/kg Q3W, 1 mg/kg Q2W, 3 mg/kg Q3W) were tested but not retained in the final model and collapse into the reference 0 group.

### NIVO_3Q2W (**canonical for nivolumab 3 mg/kg every 2 weeks co-administration indicator**)
- **Description:** 1 = ipilimumab co-administered with nivolumab 3 mg/kg every 2 weeks; 0 = otherwise (monotherapy or any other nivolumab regimen).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no nivolumab or any non-3Q2W nivolumab regimen).
- **Source aliases:**
  - `NIVO_REGIMEN` (categorical column) — decompose into `NIVO_3Q2W = as.integer(NIVO_REGIMEN == "3 mg/kg Q2W")`.
- **Example models:** `Sanghavi_2020_ipilimumab.R` (exponential coefficient 0.191 on ipilimumab CL).
- **Notes:** Paired with `NIVO_1Q3W`; same reference grouping convention.

### COMBO_NIVO (**canonical for any-regimen nivolumab combination-therapy indicator**)
- **Description:** 1 = ipilimumab co-administered with any nivolumab regimen, 0 = ipilimumab monotherapy.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (ipilimumab monotherapy).
- **Source aliases:**
  - `COMBO` — used in `Sanghavi_2020_ipilimumab.R`. Equivalently derivable from `NIVO_REGIMEN` as `COMBO_NIVO = as.integer(NIVO_REGIMEN != "none")`.
- **Example models:** `Sanghavi_2020_ipilimumab.R` (additive effect -0.202 on the Emax parameter of the time-varying CL function).
- **Notes:** Distinct from the per-regimen `NIVO_1Q3W` / `NIVO_3Q2W` indicators on baseline CL: `COMBO_NIVO` aggregates across all nivolumab regimens and acts on the time-varying-CL Emax parameter, whereas the per-regimen indicators act on baseline (time-zero) CL.

### ECOG_PS_GT0 (**canonical for ECOG performance status > 0 indicator**)
- **Description:** 1 = ECOG (Eastern Cooperative Oncology Group) performance status greater than 0 at baseline (i.e., PS 1, 2, 3, or 4); 0 = ECOG PS 0 (fully active). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (ECOG PS 0; fully active).
- **Source aliases:**
  - `PS` — used in `Zhang_2019_nivolumab.R` (paper's binary collapse of the four-level ECOG scale into PS=0 vs. PS>0).
- **Example models:** `Zhang_2019_nivolumab.R` (multiplicative effect on baseline CL, additive effect on Emax of time-varying CL).
- **Notes:** Many oncology PopPK papers report a binary-collapsed performance-status indicator rather than the full ECOG scale. The "GT0" suffix preserves the dichotomization explicitly. If a future paper uses a different cut-point (e.g., PS<=1 vs. PS>=2), register it as `ECOG_PS_GT1` rather than reusing this column.

### BLSTABL (**canonical for baseline absolute blast counts in peripheral blood**)
- **Description:** Baseline absolute count of blasts (immature lymphoid/myeloid precursor cells) circulating in peripheral blood. Time-fixed baseline value.
- **Units:** 10^9 counts/L (equivalently 10^9 counts; reported as "x 10^9 counts" in Wu 2024 Table 2).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used with power scaling `(BLSTABL / <ref>)^exponent`. Reference value observed: 0.352 x 10^9 counts (Wu 2024 Table 3, BCP-ALL median).
- **Source aliases:**
  - `BLSTABL` — used in `Wu_2024_inotuzumab.R`.
- **Example models:** `Wu_2024_inotuzumab.R` (power exponent -0.0484 on kdes for BCP-ALL patients only; the effect is gated off for B-cell NHL by multiplying the exponent by `DIS_BCPALL`).
- **Notes:** Distinct from blasts in bone marrow (different specimen) and from `BLSTPB` (percentage of blasts in peripheral blood, used by the predecessor Garrett 2019 adult model). Not applicable for B-cell NHL patients in pooled BCP-ALL + NHL analyses (Wu 2024 retains the effect only in BCP-ALL patients via the DIS_BCPALL gate). When supplying BLSTABL for an NHL subject, set the value to the BCP-ALL reference (0.352) so the gated power term evaluates to 1 numerically. Scope: specific because the covariate is most meaningful in B-cell-leukemia population PK analyses; promote to general if a second paper retains it.

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

### CONMED_RITUX (**canonical for concomitant rituximab combination therapy**)
- **Description:** 1 = on concomitant rituximab combination therapy (with or without backbone chemotherapy), 0 = not on rituximab.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no concomitant rituximab; single-agent or non-rituximab combination).
- **Source aliases:**
  - `RITUX` — used in `Wu_2024_inotuzumab.R`.
- **Example models:** `Wu_2024_inotuzumab.R` (additive fractional change on CL1: `CL1 * (1 + (-0.132) * CONMED_RITUX)` ≈ 13% lower CL1 with concomitant rituximab).
- **Notes:** Wu 2024 Table 3 footnote b explicitly flips the reference category vs. the predecessor Garrett 2019 adult model: in Garrett 2019 the reference was "with rituximab" (RITUX = 0 meant on-rituximab), whereas in Wu 2024 the reference is "without rituximab" (RITUX = 0 means no concomitant rituximab). Future models that pool an analogous rituximab-combination cohort with a single-agent reference should use this canonical with the Wu 2024 sign convention; if a paper retains the Garrett 2019 reverse-coded convention, document the value transformation in `covariateData[[CONMED_RITUX]]$notes` (`CONMED_RITUX = 1 - source$RITUX`) rather than registering a second canonical. Ratified canonically on 2026-04-26.

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

### STEROID (**canonical for baseline/concomitant systemic corticosteroid use**)
- **Description:** 1 = patient on systemic corticosteroid therapy at baseline (typically continued as concomitant medication during the study), 0 = no baseline corticosteroid use. Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no baseline corticosteroid use).
- **Source aliases:**
  - `BSTEROID` — used in `Narwal_2013_sifalimumab.R` and `Zheng_2016_sifalimumab.R`.
- **Example models:** `Narwal_2013_sifalimumab.R` (multiplicative on CL: `CL * (1 + 0.195 * STEROID)`), `Zheng_2016_sifalimumab.R` (multiplicative on CL `(1 + 0.11 * STEROID)` and on V1 `(1 - 0.09 * STEROID)` in the SLE phase IIb cohort, which was ~85% steroid-treated at baseline).
- **Notes:** Distinct from `PRICORT`, which is strictly a prior (pre-study) indicator. `STEROID` captures concurrent corticosteroid use at / from study baseline in diseases where background steroid use is standard of care (SLE, severe asthma, etc.). When a future paper needs the two jointly, both can coexist on the same subject. The name `STEROID_BL` was used as an alias in earlier register drafts and is retired; use `STEROID` for all future models.

### COADMIN_IPI_3Q3W (**canonical for nivolumab + ipilimumab 3 mg/kg q3w combination indicator**)
- **Description:** 1 = subject is receiving nivolumab in combination with ipilimumab 3 mg/kg every 3 weeks (4-dose induction); 0 = otherwise. Encodes the high-intensity ipilimumab combination regimen as a study-design covariate on nivolumab CL.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (any non-3Q3W regimen — monotherapy, chemotherapy combination, or another ipilimumab schedule).
- **Source aliases:**
  - `IPI3Q3W` — used in `Zhang_2019_nivolumab.R`.
- **Example models:** `Zhang_2019_nivolumab.R` (exponential effect on baseline CL: `exp(0.227)` ≈ 1.25 fold increase relative to monotherapy).
- **Notes:** Paired with `COADMIN_IPI_1Q6W`; both indicators can coexist in one population, but a single subject has at most one set to 1 in the Zhang 2019 cohort. The remaining ipilimumab schedules (1 mg/kg q3w x 4 induction, 1 mg/kg q12w) had no statistically significant effect on nivolumab CL and were therefore folded into the reference (0) group along with monotherapy, leaving only IPI3Q3W and IPI1Q6W as named non-reference indicators.

### COADMIN_IPI_1Q6W (**canonical for nivolumab + ipilimumab 1 mg/kg q6w combination indicator**)
- **Description:** 1 = subject is receiving nivolumab in combination with ipilimumab 1 mg/kg every 6 weeks (continuous maintenance); 0 = otherwise.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (any non-1Q6W regimen — monotherapy, chemotherapy combination, or another ipilimumab schedule).
- **Source aliases:**
  - `IPI1Q6W` — used in `Zhang_2019_nivolumab.R`.
- **Example models:** `Zhang_2019_nivolumab.R` (exponential effect on baseline CL: `exp(0.159)` ≈ 1.17 fold increase relative to monotherapy).
- **Notes:** Paired with `COADMIN_IPI_3Q3W`. See the COADMIN_IPI_3Q3W note for how the other ipilimumab schedules collapse into the reference group.

### COADMIN_CHEMO (**canonical for nivolumab + chemotherapy combination indicator**)
- **Description:** 1 = subject is receiving nivolumab in combination with platinum-based chemotherapy (gemcitabine + cisplatin, pemetrexed + cisplatin, paclitaxel + carboplatin, or platinum-doublet); 0 = otherwise. Encodes any chemotherapy coadministration as a study-design covariate on nivolumab CL.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no chemotherapy coadministration — monotherapy or any ipilimumab combination).
- **Source aliases:**
  - `CHEMO` — used in `Zhang_2019_nivolumab.R`.
- **Example models:** `Zhang_2019_nivolumab.R` (exponential effect on baseline CL: `exp(-0.104)` ≈ 0.90 fold, i.e. ~9.7% lower CL relative to monotherapy).
- **Notes:** Pooled across the four chemotherapy backbones contributing to the Zhang 2019 cohort. Promote to general scope only after a second paper reports a chemotherapy-coadministration effect with a comparable pooling convention.

### COADMIN_IPI_ANY (**canonical for any-ipilimumab-coadministration indicator**)
- **Description:** 1 = subject is receiving nivolumab in combination with any ipilimumab regimen (regardless of dose or schedule); 0 = nivolumab monotherapy or nivolumab + chemotherapy. Encodes the "is there ipilimumab in the regimen" question as a single binary covariate, distinct from the regimen-specific COADMIN_IPI_3Q3W and COADMIN_IPI_1Q6W indicators above.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (no ipilimumab coadministration).
- **Source aliases:**
  - `IPICO` — used in `Zhang_2019_nivolumab.R`.
- **Example models:** `Zhang_2019_nivolumab.R` (additive effect on the time-varying-CL Emax parameter: Emax += -0.0668 when COADMIN_IPI_ANY = 1).
- **Notes:** Logically the union of the regimen-specific indicators (COADMIN_IPI_3Q3W, COADMIN_IPI_1Q6W, plus the unmodeled 1 mg/kg q3wx4 and 1 mg/kg q12w schedules). Zhang 2019 uses it on the *time-varying* Emax (a different structural parameter from baseline CL), which is why it coexists with the regimen-specific indicators on baseline CL rather than substituting for them.

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

### FCGR3A_VV (**canonical for FCGR3A 158 V/V homozygote indicator**)
- **Description:** 1 = subject is homozygous for valine at amino-acid position 158 of the FcγRIIIa receptor (V/V), encoded by the rs396991 polymorphism in the FCGR3A gene; 0 = otherwise (heterozygote V/F or homozygote F/F pooled). The dominant V/V vs (V/F + F/F) grouping is the encoding used in the Aguiar 2021 source paper after testing dominant and recessive groupings during covariate model building.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 = V/F heterozygote or F/F homozygote (combined).
- **Source aliases:**
  - `FCGR3A` (genotype string, e.g., `"V/V"` / `"V/F"` / `"F/F"`): derive `FCGR3A_VV = as.integer(FCGR3A == "V/V")`.
  - `rs396991` (raw allele coding, often `"AA"` / `"AC"` / `"CC"` or `"GG"` / `"GT"` / `"TT"` depending on assay strand): map V allele -> 1, F allele -> 0 with the assay-specific allele convention; derive `FCGR3A_VV = as.integer(genotype is V-homozygous)`.
- **Example models:** `Aguiar_2021_ustekinumab.R` (Aguiar 2021 Table 2 footnote; logit-scale effect on subcutaneous bioavailability F: 88.8% in V/V vs 71.0% in V/F + F/F).
- **Notes:** rs396991 (FCGR3A 158V>F) is a well-studied pharmacogenetic polymorphism affecting FcγRIIIa-IgG affinity and has been associated with response to several IgG monoclonal antibodies (rituximab, infliximab, ustekinumab). The V allele is the higher-affinity variant. Document the assay-strand allele convention used in the source paper in `covariateData[[FCGR3A_VV]]$notes`. Future models that use a recessive (F/F vs V/* combined) or codominant (additive 0/1/2) coding should register a separate canonical (e.g., `FCGR3A_FF`, `FCGR3A_VCT` for V-allele count) rather than overloading `FCGR3A_VV`.

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
- **Units:** ng/mL (document per-model in `covariateData[[PCSK9]]$units` if a different unit — typically nM — is used in a given model; conversion uses a PCSK9 molecular weight of ~72 kDa, so 1 nM ≈ 72 ng/mL).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(PCSK9 / ref)^exponent`. Reference values observed: 425 ng/mL (= 5.9 nM) in `Kuchimanchi_2018_evolocumab.R` (population median).
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
  - `cat` — Papachristos 2020 (the paper writes the indicator as `cat` in the CL covariate equation; no formal column name is given in the published narrative).
- **Example models:** `Papachristos_2020_bevacizumab_pk.R`, `Papachristos_2020_bevacizumab_qss.R`, `Papachristos_2020_bevacizumab_pkpd.R` (multiplicative effect on bevacizumab CL: `CL * exp(-0.423 * SNP_ICAM1_RS1799969)` in the PK and PK/PD models; `CL * exp(-0.33 * SNP_ICAM1_RS1799969)` in the binding QSS model — mutant carriers have lower CL and higher trough levels).
- **Notes:** Time-fixed per subject. Mutant carrier rate in the Papachristos 2020 mCRC cohort: 20% (`Notes` Table 1 of the paper). The biological mechanism by which the *ICAM1* mutant slows bevacizumab clearance is unknown; the association is empirical and may be specific to mCRC.

### SNP_VEGFA_RS1570360 (**canonical for VEGF-A rs1570360 mutant indicator**)
- **Description:** Binary genotype indicator for the *VEGFA* rs1570360 single-nucleotide polymorphism (-1154 G > A; promoter region). 1 = at least one mutant (A) allele present; 0 = homozygous wild-type (GG).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (homozygous wild-type).
- **Source aliases:**
  - `cat1` — Papachristos 2020 (used as the first categorical indicator in the inter-compartmental clearance equation of the PK model; no formal column name in the narrative).
- **Example models:** `Papachristos_2020_bevacizumab_pk.R` (multiplicative effect on bevacizumab Q: `Q * exp(0.378 * SNP_VEGFA_RS1570360)` — mutant carriers have higher inter-compartmental clearance).
- **Notes:** Time-fixed per subject. Mutant carrier rate in the Papachristos 2020 mCRC cohort: 33%. The covariate is significant in the standalone PK model but does not appear in the binding QSS or PK/PD models because in those models the inter-compartmental clearance covariate effects are absorbed into the rs699947 effect on Q (PK/PD) or into the K_ss / BM0 effects (QSS).

### SNP_VEGFA_RS699947 (**canonical for VEGF-A rs699947 mutant indicator**)
- **Description:** Binary genotype indicator for the *VEGFA* rs699947 single-nucleotide polymorphism (-2578 C > A; promoter region). 1 = at least one mutant (A) allele present; 0 = homozygous wild-type (CC).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (homozygous wild-type).
- **Source aliases:**
  - `cat2` — Papachristos 2020 PK model (second categorical indicator on Q).
  - `cat` — Papachristos 2020 binding QSS model (effect on K_ss and BM0) and PK/PD model (effect on Q).
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
  - `Smoking` (case-insensitive) — used in `Ma_2020_sarilumab_anc.R`.
- **Example models:** `Ma_2020_sarilumab_anc.R` (power-form on baseline ANC: `BASE * 1.15^SMOKE`).
- **Notes:** Baseline-only indicator; does not track within-study smoking-cessation changes.

## Formulation / assay / study

### ROUTE_IV
- **Description:** 1 = subject received intravenous (IV) administration, 0 = subcutaneous (SC) administration. Per-subject (study-fixed) covariate flagging the dosing route when a population analysis pools cohorts that differ by route, with covariate effects on PK parameters that capture route-specific disposition behaviour.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (SC).
- **Source aliases:** "Admin route = IV" (categorical effect column in Yu 2022 covariate equations).
- **Example models:** `Yu_2022_ofatumumab.R` (exponential effect on R0, CL, Q, ksyn∞).
- **Notes:** This is the per-subject covariate-equation indicator, distinct from the dosing-event `cmt` column that names the target compartment. When simulating, set `ROUTE_IV = 1` for IV cohorts and dose into the central compartment; set `ROUTE_IV = 0` for SC cohorts and dose into the depot. Scope: specific because the set of parameters that differ by route is paper-specific.

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
  - `Dose` — used in `Zheng_2016_sifalimumab.R` and `Castro-Surez_2020_nimotuzumab.R`.
- **Example models:** `Zheng_2016_sifalimumab.R` (power effect on V1 with exponent 0.06), `Castro-Surez_2020_nimotuzumab.R` (binary-indicator usage `(DOSE == 50)` applying a 53 % decrease in V1 for the 50 mg cohort).
- **Notes:** Distinct from `DOSE_70MG` (binary indicator for a specific dose group in a trinary-dose design). Use `DOSE` when the source paper treats dose as a continuous covariate acting on a PK parameter, or when the paper reports a single-cohort dose-stratified shift in a PK parameter (encoded as `(DOSE == <level>)`). Scope: specific because the semantics (units, reference dose, whether the covariate applies to CL or V) are study-specific; future models using continuous dose as a covariate should either extend this entry's example-models list or register a more specific variant.

### DOSE_70MG
- **Description:** 1 = subject is on the 70 mg SC Q4W dose regimen, 0 = subject is on the 210 or 490 mg SC Q4W regimen.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (210 mg or 490 mg Q4W regimen).
- **Source aliases:** derived per subject from the trial-assigned dose level.
- **Example models:** `Kotani_2022_astegolimab.R`.
- **Notes:** Zenyatta-study categorical covariate flagging the 70 mg group (lowest dose), modeled as a −15.3% relative change on relative bioavailability. Modeled by Kotani 2022 as `70 mg vs {210 mg, 490 mg}` combined reference.

### DOSE_50MG
- **Description:** 1 = dose record is a 50 mg SC administration, 0 = all other SC doses (100, 150, 200, 300 mg) and all IV doses.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (100-300 mg SC or any IV dose).
- **Source aliases:** derived per dose record from the administered amount (`AMT`).
- **Example models:** `Othman_2014_daclizumab.R`.
- **Notes:** Othman 2014 estimated two separate absolute bioavailabilities because of non-linear dose-normalized exposure at the 50 mg SC dose — F = 0.84 for the therapeutic 100-300 mg SC range and F = 0.57 for the 50 mg SC cohort. Encoded here as a record-level indicator so the covariate effect `e_dose_50mg_f = 0.57/0.84 - 1 = -0.321` scales bioavailability only on 50 mg SC doses. For clinical-range simulation (150 mg SC Q4W Phase III regimen) leave `DOSE_50MG = 0`.

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

### STDY_LBSL (**canonical for early-phase belimumab LBSL01 / LBSL02 study indicator**)
- **Description:** 1 = subject enrolled in study LBSL01 (NCT00657007) or LBSL02 (NCT00071487) — the two early-phase belimumab studies that used a different ELISA-based bioanalytical assay; 0 = any other belimumab study in the Zhou 2021 pooled analysis. Used to switch CL and V1 magnitudes per study group (effectively an assay / early-development PK adjustment).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (later-phase studies using the electrochemiluminescence assay).
- **Source aliases:**
  - `INDR` — used in `Zhou_2021_belimumab.R` (Zhou 2021 Table 2 footnote: study indicator).
- **Example models:** `Zhou_2021_belimumab.R` (multiplicative factors 1.63 on CL and 1.26 on V1 when STDY_LBSL = 1).
- **Notes:** Conceptually similar to `STUDY1` / `PHASE2` / `ELISA` / `PHASE1` (per-study switches) but specific to the belimumab program. Subject-level (time-fixed); set from the trial identifier on each subject record.

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

### CYCLE
- **Description:** Treatment cycle number (1 = first dosing cycle, 2 = second, ...). Integer count, time-varying across a multi-cycle treatment course, incremented at each new dosing cycle.
- **Units:** (count)
- **Type:** count
- **Scope:** specific
- **Reference category:** n/a — used with a power-covariate form `CYCLE^Fm` (Fm typically negative) to capture cycle-over-cycle decline in a derived quantity such as ADC-to-payload conversion fraction.
- **Source aliases:** `CYCLE` — used in `Li_2017_brentuximab.R` with the same canonical name.
- **Example models:** `Li_2017_brentuximab.R` (exponent on the fraction of ADC that converts to MMAE by proteolytic degradation, Fm = -0.261, to reflect tumor-burden reduction across successive treatment cycles).
- **Notes:** Must be >= 1 throughout (`CYCLE^Fm` is undefined at 0). Distinct from `ooc<n>` binary-occasion indicators: `CYCLE` is an integer count, not a mutually-exclusive set of indicator columns. Data-assembly helper: set `CYCLE = floor((TIME - TIME_FIRST_DOSE) / cycle_length_days) + 1` for a fixed-interval dosing regimen.

---

## Change log

- **2026-04-26** — Deduplicated four covariate entries that were registered twice during sequential model-PR merges: (1) merged second `IGG` entry (Yang 2021 cemiplimab) into the first (Zhou 2021 belimumab) — both papers use the same canonical `IGG` column, so the entry now lists both source aliases (`BIGG`, `IGGBL`) and both reference values (14.8 g/L, 9.65 g/L); (2) merged second `BGENE21` entry (Zheng 2016 sifalimumab, misfiled under Inflammation markers) into the first (Narwal 2013 sifalimumab, under Interferon / biomarker panels); (3) merged `STEROID_BL` (Zheng 2016 sifalimumab) into `STEROID` (Narwal 2013 sifalimumab) — both encode the same baseline/concomitant corticosteroid indicator, `STEROID` is the preferred shorter name, `STEROID_BL` retired; `Zheng_2016_sifalimumab.R` updated accordingly; (4) merged `ECOG_PS_GT0` (Zhang 2019 nivolumab) into `ECOG_GE1` (Bajaj 2017 nivolumab) — `>= 1` equals `> 0` for integer ECOG scores, `ECOG_GE1` is the preferred name (consistent with a potential future `ECOG_GE2`), `ECOG_PS_GT0` retired; `Zhang_2019_nivolumab.R` updated accordingly.
- **2026-04-26** — Added `RACE_JAPANESE` (general scope) canonical entry to formalize the entry that was noted in the changelog but lacked a proper H3 section. Source alias `JAPANESE_HV` from `Wang_2017_benralizumab.R` mapped; `Wade_2015_certolizumab.R` uses `RACE_JAPANESE` directly. Added `MM` (active multiple myeloma, specific scope), `MDSAML` (MDS or AML combined indicator, specific scope), and `SPDL1` (soluble PD-L1, specific scope) canonical entries while fixing convention issues in `Ogasawara_2020_durvalumab.R`.
- **2026-04-26** — Added `DIS_DMD` (specific scope) canonical entry while extracting `Wojciechowski_2022_domagrozumab.R`. Source alias `SPOP`; orientation matches the source (1 = DMD pediatric patient, 0 = healthy adult volunteer reference). The covariate enters as a `(1 + theta * DIS_DMD)` multiplicative shift (additive on the linear scale, not exponentiated) rather than the typical `theta^DIS_DMD` form, matching Eqs. 7-8 of the paper.
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
- **2026-04-24** — Added `STATIN_MONO` and `EZE` (both scope: specific; new "Concomitant lipid-lowering medication" section) and `PCSK9` (scope: general; new "Hypercholesterolemia biomarkers" section) canonical entries while extracting `Kuchimanchi_2018_evolocumab.R`. `STATIN_MONO` carries Kuchimanchi-specific "statin monotherapy only" semantics; `EZE` captures ezetimibe use (in the Kuchimanchi 2018 dataset overwhelmingly statin+ezetimibe combination therapy).
- **2026-04-24** — Added `STATIN` (general-scope concomitant statin indicator) under `Concomitant / prior medication` and `FPCSK9` (specific-scope free PCSK9 concentration) under a new `Cardiometabolic / target biomarkers` section while extracting `Martinez_2019_alirocumab.R`.
- **2026-04-24** — Added `CYCLE` canonical entry (scope: specific; integer count; power-covariate effect `CYCLE^Fm` on the ADC-to-MMAE proteolytic conversion fraction) while extracting `Li_2017_brentuximab.R`. New entry placed under `Occasion / period (IOV)` since it is conceptually an IOV-like index, but distinct from the binary `ooc<n>` pattern.
- **2026-04-24** — Added `LMET` (general-scope binary indicator for baseline liver metastases) and `TUMTP_OTH` (specific-scope residual "other tumor types" indicator, complement of the named `TUMTP_<GROUP>` indicators in the same model) canonical entries while extracting `Quartino_2019_trastuzumab.R`. Extended `TUMTP_GC` with the Quartino 2019 source alias `TTYPE == "AGC"` and `AST` with the legacy clinical-chemistry alias `SGOT`.
- **2026-04-24** — Added `LDH` (general-scope serum lactate dehydrogenase),
  `TUMTP_SCLC` (specific-scope small-cell-lung-cancer tumor-type
  indicator extending the `TUMTP_*` decomposition pattern), `LINE_1L`
  (specific-scope first-line vs second-line-or-greater therapy
  indicator), `NIVO_1Q3W` and `NIVO_3Q2W` (specific-scope per-regimen
  nivolumab co-administration indicators on baseline CL), and
  `COMBO_NIVO` (specific-scope any-regimen nivolumab combination-therapy
  indicator on the time-varying CL Emax) canonical entries while
  extracting `Sanghavi_2020_ipilimumab.R`.
- **2026-04-26** — Added `BLBCELL` (specific-scope baseline CD19+ B cell count under `Inflammation markers`), `ROUTE_IV`, `DEVICE_AI`, `STUDY_APLIOS`, and `STUDY_MIRROR` (all specific-scope under `Formulation / assay / study`) canonical entries while extracting `Yu_2022_ofatumumab.R`. Source aliases mapped: `Bcell0`→`BLBCELL`, "Admin route = IV"→`ROUTE_IV`, "Formulation = AI"→`DEVICE_AI`, "Study = APLIOS"→`STUDY_APLIOS`, "Study = MIRROR"→`STUDY_MIRROR`.
- **2026-04-25** — Added `ECOG_PS_GT0` (general-scope binary indicator for ECOG performance status > 0; Oncology section), `COADMIN_IPI_3Q3W`, `COADMIN_IPI_1Q6W`, `COADMIN_CHEMO`, and `COADMIN_IPI_ANY` (specific-scope coadministration-regimen indicators; Concomitant / prior medication section) canonical entries while extracting `Zhang_2019_nivolumab.R`. Source aliases mapped: `PS`→`ECOG_PS_GT0`, `IPI3Q3W`→`COADMIN_IPI_3Q3W`, `IPI1Q6W`→`COADMIN_IPI_1Q6W`, `CHEMO`→`COADMIN_CHEMO`, `IPICO`→`COADMIN_IPI_ANY`.
- **2026-04-25** — Added `FFM` (general-scope, fat-free mass, Janmahasatian formula), `IGG` (general-scope, serum immunoglobulin G), `RACE_NEAS` (specific-scope, North East Asian composite race indicator), and `STDY_LBSL` (specific-scope, early-phase belimumab LBSL01/02 study indicator) canonical entries while extracting `Zhou_2021_belimumab.R`. Source aliases mapped: `BALB`→`ALB`, `BIGG`→`IGG`, `RAC4`→`RACE_NEAS`, `INDR`→`STDY_LBSL`.
- **2026-04-25** — Added `IGG` (general-scope endogenous serum immunoglobulin G concentration; placed under `Renal / hepatic function` near `LDH` since the cemiplimab paper uses it as a baseline lab covariate alongside ALB, ALT, BMI, and WT) canonical entry while extracting `Yang_2021_cemiplimab.R`. Source alias `IGGBL`→`IGG` mapped. Reference value 9.65 g/L.
- **2026-04-25** — Added `DIS_PJIA` (polyarticular juvenile idiopathic arthritis disease-state indicator; scope: specific) under the existing `Disease state (cross-population indicators)` section while extracting `Gandhi_2021_abatacept.R`, where Gandhi 2021 pools adult RA with pJIA patients and tests pJIA-vs-RA on bioavailability (additive on logit-F: +3.08). Source alias `JIA` mapped. Reused the existing `SWOL_28JOINT` canonical for Gandhi 2021's swollen-joint-count covariate per operator decision: Gandhi 2021's reported reference SJC = 15 is consistent with the 28-joint scale, the same author group used the 28-joint count in Li 2019 (RA-only), and the paper text does not explicitly identify the joint-count scale.
- **2026-04-25** — Added new top-level section `Pharmacogenomic SNPs` introducing the canonical pattern `SNP_<GENE>_<RSID>` for binary mutant-allele-presence genotype indicators. Three new entries: `SNP_ICAM1_RS1799969`, `SNP_VEGFA_RS1570360`, `SNP_VEGFA_RS699947` (all scope: specific) — the first pharmacogenomic-genotype covariates in the register, introduced while extracting the three Papachristos 2020 bevacizumab models (PK / binding QSS / PK/PD). Encoding: 1 = at least one mutant allele present; 0 = homozygous wild-type. Source-paper indicator names `cat`, `cat1`, `cat2` (which are positional within each model's covariate equation rather than formal column names) are recorded as aliases.
- **2026-04-25** — Added `DIAB` (general-scope binary diabetes-mellitus comorbidity indicator) canonical entry under a new `Comorbidities` H2 section while extracting `Chen_2022_guselkumab.R`. Distinct from a primary-disease indicator (`DIS_*`); used in non-diabetes-primary indications where diabetes is tested as a covariate. Source alias `DIAB` mapped.
- **2026-04-26** — Added `RACE_WHITE` (general; White vs non-White dichotomy), `HEPIMP_MILD` (general; mild hepatic impairment per NCI ODWG criteria), `NLR` (general; neutrophil-to-lymphocyte ratio, hematology), `SARS_VLOAD` (specific; SARS-CoV-2 baseline viral load), `SARS_SEROPOS` (specific; SARS-CoV-2 baseline serostatus), `OXYSUP_LOW` and `OXYSUP_HIGH` (both specific; decomposed from a 4-level supplemental-oxygen categorical) canonical entries while extracting `Lin_2024_casirivimab.R`. New H2 section `Infectious disease (SARS-CoV-2 / COVID-19)` introduced for the SARS_* and OXYSUP_* entries. Source aliases mapped: `RACE`→`RACE_WHITE` (Lin 2024 binary form), `HEPIMP`→`HEPIMP_MILD`, `VIRAL`→`SARS_VLOAD`, `SERPOS`→`SARS_SEROPOS`, `OXYSTAT1`→`OXYSUP_LOW`, `OXYSTAT2`→`OXYSUP_HIGH`.
- Subsequent additions: append new canonical entries as new papers are processed. When adding, bump the audit-completed count in the summary below.

## Summary

- Files audited: 69 R files under `inst/modeldb/` (20 of which reference covariates).
- Canonical H3 entries: 58 (61 parsed entries when counting each `ooc<n>` individually; was ~57 before the 2026-04-20 mergers of `hsCRP`+`BLCRP`+standard-`CRP` → `CRP`, `eGFR`+`CRCL_BSA` → `CRCL`, and `ADA_TITRE`+`ADA_TITER` → `ADA_TITER`; +1 `PHASE2` and +1 `WBC` on 2026-04-21 from Farrell 2012 / Mould 2007; +3 on 2026-04-21 for `STEROID`, `BGENE21`, `COHDOSE` from Narwal 2013).
- Scope: general: 35. Scope: specific: 26 (counting each `ooc<n>` individually, or 23 counting the `### ooc1, ooc2, ooc3, ooc4` heading as one entry).
- Aliases mapped (selected): SEXM→SEXF, ADA→ADA_POS, ADA_TITRE/ADAT→ADA_TITER, BLACK→RACE_BLACK, ASIAN→RACE_ASIAN, MULTIRACIAL→RACE_MULTI, BLACK_OTH→RACE_BLACK_OTH, ASIAN_AMIND_MULTI→RACE_ASIAN_AMIND_MULTI, DVID→STUDY1/STUDY5, CRE→CREAT, hsCRP/HSCRP/CRPHS/BLCRP→CRP, eGFR/EGFR/CRCL_BSA/1.73*CrCl/BSA→CRCL, DP2→FORM_DP2, DISEXT→DISEXT_EP/DISEXT_OTHER, BEASI→EASI, BEOS→EOS, GAST→PRIOR_GAST, COMB→COMB_EOX, TUMTP→TUMTP_CHL/TUMTP_GC, BSTEROID→STEROID, DOSE→COHDOSE.
- Canonical H3 entries: 53 (56 parsed entries when counting each `ooc<n>` individually; was ~57 before the 2026-04-20 mergers of `hsCRP`+`BLCRP`+standard-`CRP` → `CRP`, `eGFR`+`CRCL_BSA` → `CRCL`, and `ADA_TITRE`+`ADA_TITER` → `ADA_TITER`).
- Scope: general: 33. Scope: specific: 23 (counting each `ooc<n>` individually, or 20 counting the `### ooc1, ooc2, ooc3, ooc4` heading as one entry).
- Aliases mapped (selected): SEXM→SEXF, ADA→ADA_POS, ADA_TITRE/ADAT→ADA_TITER, BLACK→RACE_BLACK, ASIAN→RACE_ASIAN, MULTIRACIAL→RACE_MULTI, BLACK_OTH→RACE_BLACK_OTH, ASIAN_AMIND_MULTI→RACE_ASIAN_AMIND_MULTI, DVID→STUDY1/STUDY5, CRE→CREAT, hsCRP/HSCRP/CRPHS/BLCRP→CRP, eGFR/EGFR/CRCL_BSA/1.73*CrCl/BSA→CRCL, DP2→FORM_DP2, DISEXT→DISEXT_EP/DISEXT_OTHER, BEASI→EASI, BEOS→EOS, GAST→PRIOR_GAST, COMB→COMB_EOX, TUMTP→TUMTP_CHL/TUMTP_GC, DX→IBD_CD, AZA→CONMED_AZA, MP→CONMED_MP, MTX→CONMED_MTX, AMINO→CONMED_AMINO.
