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
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical), `CarlssonPetri_2021_liraglutide.R` (alias `SEXM`), `Bajaj_2017_nivolumab.R` (male-indicator source; effect applied as `exp(coef * (1 - SEXF))` to preserve the paper's female-reference CL_REF / VC_REF), `Fau_2020_isatuximab.R` (exponential effect on Vc; reference category 0 = male).
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
- **Example models:** `Fasanmade_2009_infliximab.R` (g/dL, reference 4.1), `Thakre_2022_risankizumab.R` (g/L, reference 45), `Chua_2025_mirikizumab.R`, `Moein_2022_etrolizumab.R`, `Tiraboschi_2025_amlitelimab.R`, `Yamada_2025_zolbetuximab.R`, `Li_2019_abatacept.R` (g/dL, reference 4.0; the Li 2019 Methods states 'mg/dL' which is a publication typo — see the model's `covariateData[[ALB]]$notes`), `Quartino_2019_trastuzumab.R` (g/dL, reference 4; source column `ALBU`; negative exponent -0.998 on linear CL), `Wang_2020_ontamalimab.R` (g/L, reference 39), `Zhou_2021_belimumab.R` (g/L, reference 40; baseline-only, source column `BALB`), `Okada_2025_rocatinlimab.R` (g/L, reference 44; source column `ALBU`; power exponent -1.30 on linear CL).
- **Notes:** Ratified canonically on 2026-04-19 after cross-model review. Unit varies by paper (g/dL in US-convention papers, g/L in SI-convention papers); the per-model `covariateData[[ALB]]$units` field is load-bearing. Effect-coefficient magnitude is meaningless without the unit.

### TPRO (**canonical for total serum protein**)
- **Description:** Total serum protein concentration (sum of albumin + globulins; baseline or time-varying).
- **Units:** g/L or g/dL — document the unit used in each model via `covariateData[[TPRO]]$units` (1 g/dL = 10 g/L).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(TPRO / ref)^exponent`. Reference value observed: 74 g/L (Frey 2010 pooled-cohort median).
- **Source aliases:**
  - `PROT` — Frey 2010 abbreviation in the final-model equation.
  - `TP` — common clinical-chemistry abbreviation.
- **Example models:** `Frey_2010_tocilizumab.R` (g/L, reference 74; exponent -1.1 on V1).
- **Notes:** Distinct from `ALB` (serum albumin, the largest single component of total protein). Frey 2010 retains both `TPRO` and `ALB` on V1 as separate covariates with opposite signs (TPRO negative, ALB positive) and notes there is no clear mechanistic explanation; the joint effect may reflect serum-volume modifications. `TPRO` ratified canonically on 2026-04-28 alongside the Frey 2010 extraction.

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

### IGM (**canonical for serum immunoglobulin M**)
- **Description:** Serum total immunoglobulin M (IgM) concentration (baseline). Used in IgRT population-PK analyses as a proxy for B-cell antibody-producing capacity / humoral function — IgM is the first antibody produced after B-cell activation, so circulating IgM reflects ongoing B-cell activity prior to class-switching to IgG.
- **Units:** g/L (typical SI-convention reporting); also reported as mg/dL in US-convention papers — document the unit used in each model via `covariateData[[IGM]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used with power scaling `(IGM / ref)^exponent`. Reference values observed: 0.21 g/L (Cheng 2026, pooled PID + SAD pediatric cohort median).
- **Source aliases:** none known.
- **Example models:** `Cheng_2026_immunoglobulin.R` (g/L, reference 0.21; baseline-only; power exponent 0.11 on baseline IgG (CBAS) — IgM enters as a humoral-capacity proxy that informs the endogenous-IgG baseline rather than directly modifying clearance).
- **Notes:** IgM is the immune-globulin class produced by activated B cells before class-switching, so it remains detectable in patients with hypogammaglobulinaemia who still have residual B-cell function. Scope: specific because the relevance of IgM as a covariate depends on the paper's mechanistic interpretation (in Cheng 2026 it acts on the endogenous-IgG baseline; future use cases may differ). Promote to general if a second paper retains IgM with consistent semantics. Ratified canonically on 2026-04-28.

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
### B2M (**canonical for serum beta-2-microglobulin**)
- **Description:** Serum beta-2-microglobulin concentration. Low-molecular-weight (~12 kDa) protein freely filtered at the glomerulus and reabsorbed in the proximal tubule; serum levels rise with renal impairment, with increased plasma-cell turnover in multiple myeloma, and with broader lymphoid-cell turnover. Used in oncology PK analyses both as a renal-function proxy and as a tumor-burden / disease-severity covariate.
- **Units:** mg/L
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(B2M / ref)^exponent`. Reference values observed: 3.90 mg/L (Fau 2020 multiple-myeloma cohort median).
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
  - `BHPTGRPN` (categorical: 1 = normal, 2 = mild, 3 = moderate, 4 = severe; 9999 = missing) — used in `Lu_2019_polatuzumab.R`. Decompose: `HEPIMP = as.integer(BHPTGRPN > 1.5 & BHPTGRPN != 9999)`.
  - `HEP_IMP` — retired canonical name; replaced by `HEPIMP` for consistency with the `HEPIMP_MILD` family.
- **Example models:** `Lu_2019_polatuzumab.R` (multiplicative effect on FRAC_NS = 1.19, applied as `1.19^HEPIMP`).
- **Notes:** NCI ODWG classification (Ramalingam SS et al., J Clin Oncol 2010;28:4507) groups subjects by total bilirubin and AST: group 1 = normal, group 2 = mild (TBILI <= ULN and AST > ULN, or TBILI > 1-1.5 x ULN), group 3 = moderate (TBILI > 1.5-3 x ULN), group 4 = severe (TBILI > 3 x ULN). Source papers typically pool groups 2-4 versus group 1 for a binary indicator because the impaired-liver subgroups are individually small. If a future model needs finer resolution (separate effects for mild vs moderate-or-worse), add a parallel `HEPIMP_MOD` canonical rather than overloading this one.
### CPK (**canonical for serum creatine phosphokinase / creatine kinase**)
- **Description:** Serum creatine phosphokinase (also called creatine kinase, CK) activity (baseline or time-varying). Skeletal-muscle / cardiac-muscle injury and turnover marker; in macrophage-targeted PK/PD analyses (axatilimab, anti-CSF-1R) it is interpreted as a Kupffer-cell / tissue-macrophage clearance surrogate because Kupffer cells participate in the elimination of circulating muscle-derived enzymes.
- **Units:** U/L (IU/L; interchangeable). Document per-model via `covariateData[[CPK]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(CPK / ref)^exponent`. Reference values observed: 63 U/L (Yang 2024 axatilimab; pooled-cohort median).
- **Source aliases:**
  - `BLCPK` (baseline CPK) — informal usage in Yang 2024.
- **Example models:** `Yang_2024_axatilimab.R` (baseline-only covariate on baseline NCMC concentration `BL_NCMC` with power exponent 0.376; reference 63 U/L).
- **Notes:** Muscle-origin enzyme distinct from `AST` / `ALT` (hepatic) and `LDH` (general tissue turnover). Yang 2024 uses CPK alongside `AST` and `LDH` as tracked safety biomarkers. Per-model `covariateData[[CPK]]$notes` should document baseline-vs-time-varying status and the clinical interpretation in the source population (skeletal-muscle injury, macrophage-clearance surrogate, or both). Distinct from any model state variable representing CPK time-course dynamics — covariate column is the pre-dose laboratory observation.

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

### HCT (**canonical for hematocrit**)
- **Description:** Hematocrit — packed red blood cell volume fraction (baseline or time-varying).
- **Units:** % (volume fraction times 100). Document per-model via `covariateData[[HCT]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(HCT / ref)^exponent`. Reference values observed: 45 % (Nestorov 2014, study-population median for severe hemophilia A adults).
- **Source aliases:** none; `HCT` is the universal NONMEM / clinical-PK abbreviation.
- **Example models:** `Nestorov_2014_factorviii.R` (reference 45 %, exponent -0.419 on V1).
- **Notes:** Higher HCT (more red-cell volume) leaves a smaller plasma fraction within total body volume; for plasma-restricted distribution (e.g., factor VIII activity, which circulates in plasma) the central volume of distribution decreases as HCT rises, so the exponent is negative. Document baseline-vs-time-varying status in `covariateData[[HCT]]$notes`. Distinct from `HGB` (mass concentration of hemoglobin); the two correlate but enter different mechanistic relationships.

## Coagulation / hemostasis biomarkers

### VWF (**canonical for von Willebrand factor concentration**)
- **Description:** Plasma concentration (or activity) of von Willebrand factor (VWF) — the multimeric carrier protein that binds and protects circulating factor VIII (FVIII) from proteolytic degradation and rapid clearance. Used as a covariate on FVIII (and FVIII-Fc) clearance because the vast majority (>95%) of circulating FVIII is in complex with VWF.
- **Units:** IU/dL (equivalent to % of pooled normal plasma); document per-model via `covariateData[[VWF]]$units`. Some sources report `VWF:Ag` (antigen) versus `VWF:RCo` (ristocetin cofactor activity); record which assay was used in `covariateData[[VWF]]$notes`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(VWF / ref)^exponent`. Reference values observed: 118 IU/dL (Nestorov 2014, study-population median).
- **Source aliases:** none; `VWF` is the universal abbreviation. Source papers may write `vWF` (lowercase v) or specify the assay (`VWF:Ag`).
- **Example models:** `Nestorov_2014_factorviii.R` (reference 118 IU/dL, exponent -0.343 on CL; VWF antigen).
- **Notes:** Higher VWF protects FVIII from clearance, so the exponent on CL is negative. VWF is time-varying within an individual (acute-phase response, age, blood group, etc.), but most published population PK models use baseline-only VWF when the within-subject dynamics are not characterized; document the per-model convention in `covariateData[[VWF]]$notes`.

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

### ACUTE_MED_DAYS (**canonical for baseline number of days/month of acute migraine medication use**)
- **Description:** Baseline number of days per month on which acute migraine medication (triptans or ergot compounds) was used during the 28-day run-in period prior to first dose. Enters as a piecewise-linear shift on baseline migraine or moderate-to-severe headache days in migraine exposure-response models.
- **Units:** days/month
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — piecewise-linear shift with breakpoint at 5 d/mo: contributes 0 below 5 and `slope * (ACUTE_MED_DAYS - 5)` above 5 (Fiedler-Kelly 2020). The 5-day breakpoint reflects the clinical guideline for medication-overuse headache.
- **Source aliases:** "Baseline days/month of acute medications" — used in `FiedlerKelly_2020_fremanezumab_em.R` and `FiedlerKelly_2020_fremanezumab_cm.R`.
- **Example models:** `FiedlerKelly_2020_fremanezumab_em.R` (slope 0.438 d/d, episodic migraine), `FiedlerKelly_2020_fremanezumab_cm.R` (slope 0.460 d/d, chronic migraine).
- **Notes:** Specific scope because the variable is migraine-domain-bound. Time-fixed per subject (baseline-only). When future migraine E-R models register additional aliases or alternative breakpoints, document them per-model and consider promoting to `general`.

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
- **Reference category:** n/a — used with power scaling `(IL6 / ref)^exponent`.
- **Source aliases:**
  - `bIL6`, `IL6_BASE` (baseline-IL6 column variants).
- **Example models:** `Frey_2013_tocilizumab.R` (power effect on baseline target receptor).
- **Notes:** Time-fixed at baseline unless the source paper states otherwise. Time-varying IL-6 in PD models should be the model's predicted state (a `d/dt(IL6)` trajectory) rather than a covariate; use this register entry only for the baseline / observed-covariate role.

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
- **Source aliases:**
  - `Bcell0` — used in `Yu_2022_ofatumumab.R`.
  - `BBCC` (NHL Phase I/Ib/II convention; values in 10^6 cells/L = cells/uL) — used in `Lu_2019_polatuzumab.R`.
- **Example models:** `Yu_2022_ofatumumab.R` (power effect on the maximum B-cell-lysis stimulatory effect Emax, exponent 0.275, reference 200 cells/µL), `Lu_2019_polatuzumab.R` (two distinct effects: power on CL_INF with input floored at 1 cell/uL, and a thresholded power on CL_T with the BLBCELL/121-cells/uL ratio floored at 1).
- **Notes:** Distinct from a *time-varying* B cell count, which is the PD response variable rather than a covariate. Scope: specific because the clinically relevant baseline depends on the surface marker (CD19, CD20, CD22) and whether the panel reports total B cells or memory/naive subsets — register a new canonical name if a future paper uses a different marker. Both Yu 2022 (anti-CD20 ofatumumab) and Lu 2019 (anti-CD79b polatuzumab vedotin) use CD19+ counts, so the canonical is reused; subtype-specific differences are documented in each model's `covariateData[[BLBCELL]]$notes`.

### CSF1 (**canonical for colony-stimulating factor 1 / macrophage-colony-stimulating factor concentration**)
- **Description:** Plasma colony-stimulating factor 1 (CSF-1, also known as macrophage colony-stimulating factor, M-CSF) concentration (baseline or time-varying). The hematopoietic cytokine that signals through CSF-1R to drive monocyte / macrophage differentiation and survival; used as both a target-engagement biomarker (anti-CSF-1R mAbs increase circulating free CSF-1) and a baseline covariate on PK / PD parameters in CSF-1R-pathway PopPK/PD models.
- **Units:** pg/mL (= ng/L; the two labels are numerically equivalent). Document per-model via `covariateData[[CSF1]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used with power scaling `(CSF1 / ref)^exponent`. Reference values observed: 549 pg/mL (Yang 2024 axatilimab; pooled-cohort median).
- **Source aliases:** none formally; informal aliases include `BLCSF1` (baseline CSF-1) and `BL_CSF1` (model-parameter notation in Monolix / NONMEM control streams).
- **Example models:** `Yang_2024_axatilimab.R` (baseline-only covariate on linear clearance `CL` with power exponent 0.912 and on the model parameter `BL_CSF1` with power exponent 0.656; reference 549 pg/mL).
- **Notes:** Specific scope because the column is meaningful only for CSF-1R-pathway-targeting drugs (axatilimab and future anti-CSF-1R molecules). Distinct from any CSF-1 model state representing time-course dynamics — covariate column is the pre-dose laboratory observation, typically measured by an ELISA assay (Yang 2024 used the R&D Systems Quantikine ELISA). Per-model `covariateData[[CSF1]]$notes` should document the assay used and any LOQ-related imputation for samples below the assay's limit of detection.

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

### HDLC (**canonical for high-density lipoprotein cholesterol**)
- **Description:** Serum high-density lipoprotein cholesterol concentration (baseline or time-varying).
- **Units:** mg/dL or mmol/L — document the unit used in each model via `covariateData[[HDLC]]$units` (1 mmol/L ≈ 38.67 mg/dL for cholesterol).
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(HDLC / ref)^exponent`. Reference value observed: 54 mg/dL (Frey 2010 pooled-cohort median).
- **Source aliases:**
  - `HDL-C` — Frey 2010 spelling with hyphen.
  - `HDL_C` — common alternative spelling.
- **Example models:** `Frey_2010_tocilizumab.R` (mg/dL, reference 54; small negative exponent -0.2 on linear CL; the paper interprets the effect as a body-size surrogate rather than a mechanism).
- **Notes:** Cardiometabolic lipid-panel covariate. In Frey 2010 it correlates with body size (HDL-C is lower in larger patients) and the small CL effect (-14% to +15% across the observed 23-135 mg/dL range) was retained but not interpreted as mechanistic.

### FPCSK9 (**canonical for free (unbound) proprotein convertase subtilisin/kexin type 9 concentration**)
- **Description:** Free (unbound, non-drug-bound) serum proprotein convertase subtilisin/kexin type 9 (PCSK9) concentration. For anti-PCSK9 monoclonal antibodies (alirocumab, evolocumab, bococizumab) the free-PCSK9 pool is the pharmacologically active target fraction; drug–target binding reduces FPCSK9 relative to total PCSK9.
- **Units:** ng/mL (document per-model if a paper reports a different unit).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used with linear-deviation forms `TVPARAM + theta * (FPCSK9 / ref)` or power-form `(FPCSK9 / ref)^theta`. Reference values observed: 72.9 ng/mL (Martinez 2019 time-varying median).
- **Source aliases:** none known.
- **Example models:** `Martinez_2019_alirocumab.R` (time-varying; additive-linear effect on `Km` with slope −0.541 per (FPCSK9/72.9), reference 72.9 ng/mL).
- **Notes:** Specific scope because the column is meaningful only for drugs whose mechanism is PCSK9 inhibition; reusing the name for a different anti-PCSK9 agent is acceptable (add to Example models). For non-PCSK9 drugs that use a similar target-concentration biomarker, register a new canonical (e.g., `FIL6R`, `FTNF`) rather than overloading `FPCSK9`. Per-model `covariateData[[FPCSK9]]$notes` should state whether the value is baseline-only or time-varying and how missing values were imputed (Martinez 2019 used LOCF).

### SBCMA (**canonical for soluble B-cell maturation antigen concentration**)
- **Description:** Baseline serum (or plasma) concentration of soluble B-cell maturation antigen (sBCMA), the shed extracellular domain of the BCMA receptor (TNFRSF17). Serves as a soluble-target biomarker for BCMA-directed therapeutics (anti-BCMA antibody-drug conjugates such as belantamab mafodotin, BCMA-targeted bispecifics, and BCMA CAR-T) — sBCMA is elevated in multiple-myeloma and reflects tumour burden, and contributes to target-mediated drug disposition by sequestering circulating drug.
- **Units:** ng/mL (equivalent to μg/L; 1 ng/mL = 1 μg/L). Document per-model via `covariateData[[SBCMA]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used with power scaling `(SBCMA / ref)^exponent`. Reference value observed: 50 ng/mL (Papathanasiou 2025 typical-patient definition).
- **Source aliases:**
  - `SBCMABL` (baseline soluble BCMA) — used in `Papathanasiou_2025_belantamab.R`.
- **Example models:** `Papathanasiou_2025_belantamab.R` (ng/mL, reference 50; power exponents on initial CL +0.113, on ADC Vc +0.0401, on Imax +0.160).
- **Notes:** Specific scope because the column is meaningful only for drugs whose mechanism involves the BCMA receptor (and thus a circulating soluble-target pool). Reusing the name for another anti-BCMA agent is acceptable (extend the example-models list). For other oncology TMDD targets register a new canonical (e.g., `HER2_ECD` already exists for HER2; an analogous `SCD20`, `SCD38` would follow the same pattern). Multiple myeloma populations show sBCMA spanning roughly 2 to 2,000 ng/mL, so the (SBCMA/50)^exponent form should be evaluated with care over the full clinical range.
## Drug exposure metrics

### CAV (**canonical for average drug plasma concentration over a dosing interval**)
- **Description:** Average plasma concentration of the modelled drug over a dosing interval (Cav = AUC_tau / tau). Used as the time-varying or per-period exposure metric in exposure-response models that feed individual empirical-Bayes PK predictions from a previously published population PK model into a downstream PD model.
- **Units:** ug/mL (document per-model via `covariateData[[CAV]]$units`).
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used in Emax/EC50 (e.g., `Emax * CAV / (EC50 + CAV)`) or power (e.g., `(CAV / CavMedian)^exponent`) drug-effect terms. Set to 0 for placebo periods.
- **Source aliases:** `CAV`, `Cav`, `CAVG`.
- **Example models:** `FiedlerKelly_2020_fremanezumab_em.R`, `FiedlerKelly_2020_fremanezumab_cm.R`.
- **Notes:** Specific scope because the value is intrinsically tied to the modelled drug — there is no shared meaning across drugs or studies. Each model's `covariateData[[CAV]]$notes` should state how the Cav values are derived (e.g., empirical-Bayes from a referenced population PK model) and that the column is set to 0 for placebo periods.

### IGE (**canonical for serum total immunoglobulin E concentration**)
- **Description:** Baseline serum total immunoglobulin E concentration (free IgE plus, in patients on anti-IgE therapy, omalizumab–IgE complex). For anti-IgE monoclonal antibodies (omalizumab, ligelizumab) IgE is the pharmacologic target; baseline IgE sets the magnitude of the target sink and modifies free-IgE clearance and the rate of IgE production in mechanism-based binding/turnover models.
- **Units:** ng/mL (typical clinical-PK convention). Pretreatment values reported in `IU/mL` are converted via `1 IU/mL = 2.42 ng/mL` (Hayashi 2007 Methods). Document per-model via `covariateData[[IGE]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(IGE / ref)^exponent`. Reference value observed: 482.4 ng/mL (Hayashi 2007 Japanese atopic-asthma cohort).
- **Source aliases:**
  - `IgE0` (baseline IgE concentration) — used in `Hayashi_2007_omalizumab.R`.
- **Example models:** `Hayashi_2007_omalizumab.R` (ng/mL, reference 482.4; power exponents −0.281 on apparent CL of free IgE and +0.657 on apparent IgE production rate; also used as the initial value for the total-IgE state at t = 0).
- **Notes:** General scope because baseline serum total IgE is a routine clinical-laboratory measurement, not a target tied to one drug. In mechanism-based anti-IgE binding/turnover models the in-model IgE state is a separate dynamic variable (`X_TE`, in nmol or nmol/L) — `IGE` is the per-subject baseline column used for covariate scaling and (when applicable) state initialization, not the dynamic state itself. For models that use the alternative reporting unit `IU/mL`, multiply by 2.42 before applying the canonical-units (ng/mL) reference value, or document the per-model unit choice in `covariateData[[IGE]]$units` so downstream tooling can interpret the values correctly.

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
- **Reference category:** 0 (non-White; complement composition depends on the source paper, typically pooling Black/African American, Asian, American Indian/Alaska Native, Native Hawaiian/Pacific Islander, Other, Not reported, Unknown). Some papers (e.g., Hu 2014) instead use the Caucasian (RACE_WHITE = 1) subgroup as the typical-value reference; the column encoding is unchanged but the model implements the effect on `(1 - RACE_WHITE)`.
- **Source aliases:**
  - `RACE` (with values `1 = White / 0 = non-White`) — used in `Lin_2024_casirivimab.R`. Source column name `RACE` is generic; the canonical name is intentionally explicit because some other models use `RACE` for a different dichotomy.
  - `RACE` (Caucasian-vs-non-Caucasian dichotomy as named in Hu 2014 Table 2) — used in `Hu_2014_bapineuzumab.R`. Same canonical column name and 1 = White / 0 = non-White encoding; the typical-value reference is the Caucasian subgroup, so the model implements the 15% non-Caucasian effect on `(1 - RACE_WHITE)`.
- **Example models:** `Lin_2024_casirivimab.R` (multiplicative fractional change on CL relative to non-White reference), `Hu_2014_bapineuzumab.R` (multiplicative 15% increase in CL for non-Caucasian relative to Caucasian reference).
- **Notes:** Used by papers that dichotomize race as White vs. non-White rather than decomposing into separate group indicators. Sign and reference-category interpretation are inverted relative to `RACE_BLACK` / `RACE_ASIAN` / etc.; do NOT combine `RACE_WHITE` with the decomposed indicators in the same model. The model's typical-value reference category (which subgroup gets the unmodified `lcl` / `lvc`) varies between papers — Lin 2024 uses non-White as the reference, Hu 2014 uses Caucasian (White) as the reference; both share the same canonical column encoding.

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
- **Source aliases:** `ASIAN` — used in `Hu_2026_clesrovimab.R`, `Robbie_2012_palivizumab.R`, `Fau_2020_isatuximab.R`. `RAAS` (race-Asian-vs-other indicator as named in Bajaj 2017 Table 1) — used in `Bajaj_2017_nivolumab.R`.
- **Example models:** `Zhu_2017_lebrikizumab.R` (canonical form), `Robbie_2012_palivizumab.R`, `Bajaj_2017_nivolumab.R`, `Fau_2020_isatuximab.R`.

### RACE_ASIAN_AMIND_MULTI (**canonical for composite group**)
- **Description:** 1 = Asian, American Indian / Alaskan Native, or Multiple races, 0 = other.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = White, Black / African American, Native Hawaiian / Pacific Islander, or Other (Clegg 2024 grouping).
- **Source aliases:** `ASIAN_AMIND_MULTI` — used in `Clegg_2024_nirsevimab.R`.
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
- **Notes:** Distinct from `RACE_ASIAN_AMIND_MULTI` (a 4-way composite of Asian + American Indian + Multiple Races) because the underlying paper's grouping rule is different — `RACE_ASIAN_OTH` is a within-Asian-population sub-indicator, not a multi-race composite. Operator decision (2026-04-28): kept separate from `RACE_ASIAN` because the paper's "Other Asian" category is its own grouping, not an alias of "Asian (any)".

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
## Geographic / enrollment-country indicators

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

### DIS_CANCER (**canonical for advanced-solid-tumor / oncology cohort indicator**)
- **Description:** 1 = patient with an advanced or metastatic solid tumor (the oncology cohort in a pooled multi-indication PK/PD analysis), 0 = non-oncology subject (healthy volunteer or non-oncology disease cohort pooled in the source analysis). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-oncology subject; the complement group is paper-defined — typically the union of healthy volunteers and a non-oncology disease cohort such as cGVHD pooled in the source analysis).
- **Source aliases:** none; source NONMEM / Monolix control streams typically derive the indicator from a `POP` or `STUDY` categorical alongside `DIS_HV`.
- **Example models:** `Yang_2024_axatilimab.R` (multiplicative effect on baseline NCMC: `BL_NCMC × exp(1.22 × DIS_CANCER + 0.618 × DIS_HV)`; reference category cGVHD when both indicators are 0).
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
- **Example models:** `Nikanjam_2019_siltuximab.R` (multiplicative effects: 0.77 on CL, 0.83 on Vss; reference category is the pooled non-HV oncology cohort), `Yang_2024_axatilimab.R` (multiplicative effect on baseline NCMC: `BL_NCMC × exp(1.22 × DIS_CANCER + 0.618 × DIS_HV)`; reference category cGVHD).
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
### DIS_PNH (**canonical for paroxysmal nocturnal hemoglobinuria indicator**)
- **Description:** 1 = paroxysmal nocturnal hemoglobinuria (PNH) patient, 0 = non-PNH subject (healthy volunteer or another indication pooled in the source analysis). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-PNH subject; the complement group is paper-defined — for Lin 2024 it pools healthy adult volunteers and CHAPLE disease patients).
- **Source aliases:** none known; source NONMEM control streams typically use ad-hoc names (e.g., `PNH`, `DPNH`).
- **Example models:** `Lin_2024_pozelimab.R` (additive-fractional +34.07% effect on Vc; no CL or Vp effect; reference category pools healthy volunteers and CHAPLE patients).
- **Notes:** Paroxysmal nocturnal hemoglobinuria is a rare hematological disease characterized by uncontrolled complement activation on red blood cells; treated with C5-targeted complement inhibitors (eculizumab, ravulizumab, pozelimab). Scope: specific because the disease-pooling reference category is paper-defined. Ratified canonically on 2026-04-27.

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
- **Reference category:** 0 (non-AD subject; the complement group is paper-defined — for Pérez-Ruixo 2025 the reference is the pooled healthy-volunteer cohort).
- **Source aliases:** none known; source NONMEM control streams typically use ad-hoc names (e.g., `AD`, `STATUS`, `DISGRP`).
- **Example models:** `PerezRuixo_2025_posdinemab.R` (acts on baseline free p217+tau in CSF, R0; healthy R0 = 0.793 pmol/L vs AD R0 = 5.995 pmol/L, a 656% relative increase, no PK-parameter effects).
- **Notes:** Used when a population PK/PD model pools healthy volunteers with Alzheimer's disease patients and the AD-vs-HV contrast is retained as a covariate on a target-related parameter (e.g., baseline p-tau, baseline p217+tau). Scope: specific because the complement reference category is paper-defined. Ratified canonically on 2026-04-28.

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
### MM_NIGG (**canonical for non-IgG multiple myeloma immunoglobulin-type indicator**)
- **Description:** 1 = patient with non-IgG-secreting multiple myeloma (e.g., IgA, IgD, IgE, IgM, light-chain-only / Bence Jones, or non-secretory MM), 0 = patient with IgG-secreting multiple myeloma.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (IgG MM).
- **Source aliases:**
  - `Ig_type` — used in `Fau_2020_isatuximab.R`. Values 0 / 1 with the same orientation as the canonical (1 = non-IgG MM).
- **Example models:** `Fau_2020_isatuximab.R` (exponential effect on the steady-state linear CL CLinf with coefficient -0.751, and on the time-varying-CL half-time KCL with coefficient -0.931).
- **Notes:** Within-disease (multiple-myeloma) immunoglobulin-subtype stratifier. The mechanistic rationale (Fau 2020) is that endogenous IgG monoclonal protein in IgG-MM patients competes with the therapeutic IgG mAb for FcRn-mediated salvage, raising the therapeutic mAb's catabolic clearance; non-IgG-MM patients lack that competition and exhibit lower therapeutic-mAb clearance. Distinct from the disease-state indicators (`DIS_SMM` = smoldering MM); applies only after a multiple-myeloma diagnosis is established. Scope: specific because the comparison is a within-MM stratifier rather than a cross-population indicator.
### DIS_PSORIASIS (**canonical for plaque psoriasis disease-state indicator**)
- **Description:** 1 = plaque psoriasis patient, 0 = non-psoriasis subject (e.g., atopic dermatitis, ulcerative colitis, or healthy volunteer). Time-fixed per subject.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-psoriasis subject; the complement group is paper-defined — the union of other disease cohorts pooled in the source analysis).
- **Source aliases:** none known; source NONMEM control streams typically use a categorical `DIS` indicator (e.g., Okada 2025: `DIS=1` for psoriasis, `DIS=0` for healthy, `DIS=2` for UC, `DIS=3` for AD), decomposed into a binary `DIS_PSORIASIS` indicator at ingestion.
- **Example models:** `Okada_2025_rocatinlimab.R` (multiplicative shift `1 - 0.372` on linear CL when 1; reference complement is the pooled atopic dermatitis + ulcerative colitis + healthy-volunteer cohort).
- **Notes:** Used when a population PK model pools plaque-psoriasis patients with a non-psoriasis reference population and psoriasis disease status is retained as a covariate. Scope: specific because the disease-pooling reference category is paper-defined. Ratified canonically on 2026-04-27.

## Oncology

### TUMSZ (**canonical for baseline tumor size**)
- **Description:** Baseline tumor size. For solid tumors, the sum of diameters of target lesions per RECIST; for classical Hodgkin lymphoma and lymphoma generally, the sum of products of perpendicular diameters (SPPD) or the sum of linear diameters of target lesions, depending on the source paper.
- **Units:** mm
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(TUMSZ / ref)^exponent`. Reference values observed: 41 mm (Zhou 2025); 63 mm (Budha 2023); 90 mm (Lu 2014, source reference 9 cm converted to mm).
- **Source aliases:**
  - `LDIAM` (Zhou 2025; pediatric lymphoma "linear diameter" of target lesions in mm).
  - `TMBD` (originally in cm; `TUMSZ_mm = TMBD_cm * 10`) — used in `Lu_2014_trastuzumabemtansine.R`.
- **Example models:** `Budha_2023_tislelizumab.R` (reference 63 mm), `Lu_2014_trastuzumabemtansine.R` (reference 90 mm; source column TMBD in cm, values converted to mm on ingestion), `Zhou_2025_brentuximab.R` (reference 41 mm; source column LDIAM is the sum of linear diameters of target lesions; effect on ADC clearance only).
- **Notes:** Promoted to scope: general on 2026-04-20 as a conventional oncology baseline-tumor-size measure (RECIST for solid tumors, SPPD or sum-of-linear-diameters for lymphomas). The SPPD vs sum-of-diameters vs sum-of-linear-diameters convention is pooled onto a single column; document the per-model mixture where relevant. When the source paper reports tumor size in cm, convert to mm (the canonical unit) on data ingestion and scale the per-model reference accordingly so `(TUMSZ / ref)^exp` is numerically invariant. When a source paper specifically reports the RECIST 1.1 "sum of longest diameters" of target lesions, use the more specific `TUM_SLD` canonical instead — `TUMSZ` remains the pooled-tumor-burden register.

### TUM_SLD (**canonical for sum of longest diameters of target lesions**)
- **Description:** Baseline sum of longest diameters of target lesions per RECIST 1.1. More specific than the pooled `TUMSZ` canonical; use `TUM_SLD` when the source paper explicitly reports "sum of longest diameters" (or "sum of lesions") as the tumor-burden metric, distinct from the pooled "sum of diameters / SPPD / sum of linear diameters" mixture covered by `TUMSZ`.
- **Units:** mm
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling `(TUM_SLD / ref)^exponent`. Reference values observed: 70.0 mm (de Vries Schultink 2020 zenocutuzumab population median).
- **Source aliases:**
  - `SoL` / "sum of lesions" (de Vries Schultink 2020 zenocutuzumab) — same construct, mm.
- **Example models:** `deVriesSchultink_2020_zenocutuzumab.R` (reference 70.0 mm; power exponent 0.447 on Vmax of the parallel non-linear / Michaelis-Menten clearance).
- **Notes:** Distinct from `TUMSZ` (pooled tumor-size canonical covering RECIST sum-of-diameters / SPPD / sum-of-linear-diameters); `TUM_SLD` is the precise RECIST 1.1 sum-of-longest-diameters metric. Ratified canonically on 2026-04-29 alongside the pilot bispecific extraction (de Vries Schultink 2020 zenocutuzumab). When the source paper reports tumor burden in cm, convert to mm (the canonical unit) on data ingestion and scale the per-model reference accordingly so `(TUM_SLD / ref)^exp` is numerically invariant.

### TUMTP_CHL (**canonical for classical Hodgkin lymphoma tumor-type indicator**)
- **Description:** 1 = classical Hodgkin lymphoma (cHL) or Hodgkin lymphoma generally, 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (e.g., NSCLC, EC, HCC, UC, GC, CRC, NPC, OC, "Other" solid tumors in the Budha 2023 cohort; systemic anaplastic large-cell lymphoma in the Zhou 2025 pediatric cohort).
- **Source aliases:**
  - `TUMTP` (categorical column with levels like `cHL`, `GC`, ...) — decompose into `TUMTP_CHL = as.integer(TUMTP == "cHL")`.
  - `DIS` (Zhou 2025; integer code with `DIS == 1` flagging HL) — decompose into `TUMTP_CHL = as.integer(DIS == 1)`. Zhou 2025 calls the complement "non-HL"; in the Zhou 2025 cohort the non-HL group is exclusively sALCL.
- **Example models:** `Budha_2023_tislelizumab.R`, `Zhou_2025_brentuximab.R` (effects on ADC Q2, MMAE central volume VM, and the ADC->MMAE conversion-decay rate ALFM; the Zhou 2025 paper anchors typical-value parameters to HL patients so the model uses `(1 - TUMTP_CHL)` as the on-effect indicator with reference category 1 = HL).
- **Notes:** Paired with `TUMTP_GC` in Budha 2023; a patient can have at most one of the indicators set to 1 (the remaining tumor types collapse into the reference 0 group). The reference category is the off-encoded value (0) by definition; when a source paper anchors typical-value parameters to the HL group rather than the non-HL group (as Zhou 2025 does), encode the effect as `coef^(1 - TUMTP_CHL)` so the canonical column meaning (1 = cHL/HL) is preserved while the paper's reference (HL) still receives multiplier 1.

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
- **Description:** 1 = gastric cancer (GC) or adenocarcinoma of the gastroesophageal junction (GEJ), 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (same reference group as `TUMTP_CHL`).
- **Source aliases:**
  - `TUMTP` (categorical column) — decompose into `TUMTP_GC = as.integer(TUMTP == "GC")`.
  - `TTYPE` (Quartino 2019; categorical column with levels `MBC`, `EBC`, `HV`, `AGC`, `Others`) — decompose into `TUMTP_GC = as.integer(TTYPE == "AGC")`.
  - `TTYPE4` (Wang 2024; level 4 of a five-level tumor-type factor labelled "GCGEJ" in the source) — decompose into `TUMTP_GC = as.integer(TTYPE4 == 1)`.
- **Example models:** `Budha_2023_tislelizumab.R`, `Quartino_2019_trastuzumab.R` (advanced gastric cancer; per-group typical-value switch on linear CL and Vc rather than an exponential multiplier), `Wang_2024_sugemalimab.R` (gastric + GEJ adenocarcinoma pooled; exponential coefficient log(1.13) on CL and log(1.14) on Vc).
- **Notes:** Follows the `RACE_<GROUP>` indicator-decomposition pattern. New oncology tumor types should be added as additional `TUMTP_<GROUP>` entries so the reference set stays explicit. "Advanced gastric cancer" (AGC), "gastric cancer" (GC), and "GC or adenocarcinoma of the gastroesophageal junction" (GCGEJ) are pooled onto a single `TUMTP_GC` indicator; document the per-paper stage-of-disease and GEJ-inclusion detail in `covariateData[[TUMTP_GC]]$notes`. ESCC (squamous histology) is captured by the separate `TUMTP_ESCC` indicator and is not pooled here.

### TUMTP_OTH (**canonical for 'other tumor types' residual indicator**)
- **Description:** 1 = heterogeneous "other" tumor-type pool (typically NSCLC plus miscellaneous solid tumors such as prostate, ovarian, and colorectal), 0 = one of the named tumor-type groups in the same analysis.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = named tumor-type groups defined per-paper (e.g., MBC, EBC, HV, AGC in Quartino 2019). The complement of all `TUMTP_<GROUP>` indicators defined in the same model.
- **Source aliases:**
  - `TTYPE` (Quartino 2019) — decompose into `TUMTP_OTH = as.integer(TTYPE == "Others")`.
  - `PAT2` (Sathe 2024) — integer-coded tumor type column with levels 1 (mTNBC), 2 (mUC or HR+/HER2- mBC), 4 (Other epithelial); the source NONMEM control stream collapses PAT2 = 1 and PAT2 = 2 into the reference (no effect) and applies the deviation only when PAT2 = 4. Decompose into `TUMTP_OTH = as.integer(PAT2 == 4)`.
- **Example models:** `Quartino_2019_trastuzumab.R` (per-group typical-value switch on linear CL; NSCLC plus a small residual group of prostate, ovarian, and other histologies), `Sathe_2024_sacituzumab.R` (multiplicative effect on tAB CL: -13.4% when TUMTP_OTH = 1; "Other" pool = small-cell and non-small-cell lung cancer, colorectal cancer, esophageal cancer, pancreatic ductal adenocarcinoma, etc., n = 184; reference = pooled mTNBC + mUC + HR+/HER2- mBC, n = 345).
- **Notes:** Scope: specific because the set of histologies collapsed into "Others" is defined by the analysis plan of the source paper; two papers' `TUMTP_OTH` columns are not interchangeable. Document the exact per-paper composition (e.g., "NSCLC + prostate + ovarian + other, n = 107 in Quartino 2019"; "small-cell + non-small-cell lung + CRC + esophageal + pancreatic ductal adenocarcinoma, n = 184 in Sathe 2024") in `covariateData[[TUMTP_OTH]]$notes`. A given subject can have at most one of the `TUMTP_<GROUP>` indicators (including `TUMTP_OTH`) set to 1; all-zero means the reference group.

### SPDL1 (**canonical for soluble PD-L1 concentration**)
- **Description:** Baseline (or time-varying) serum concentration of soluble programmed death-ligand 1 (sPD-L1). Serves as a circulating biomarker of target burden and immune activation for anti-PD-1/PD-L1 antibodies.
- **Units:** pg/mL
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used with power scaling `(SPDL1 / ref)^exponent`. Reference value observed: 173.8 pg/mL (study-population median in Ogasawara 2020).
- **Source aliases:** none; `SPDL1` is the standard abbreviation used directly in source analyses.
- **Example models:** `Ogasawara_2020_durvalumab.R` (power effect on CL, exponent 0.0617, reference 173.8 pg/mL; time-varying; values below LOD imputed as LOD/2 = 33.55 pg/mL).
- **Notes:** Scope: specific because sPD-L1 is meaningful only for drugs targeting the PD-1/PD-L1 pathway. For other checkpoint biomarkers (e.g., soluble CTLA-4, soluble LAG-3) register new dedicated canonicals rather than reusing this one. Ratified canonically on 2026-04-26.
  - `TTYPE3` (Wang 2024; level 3 of a five-level tumor-type factor labelled "Other" in the source) — decompose into `TUMTP_OTH = as.integer(TTYPE3 == 1)`.
- **Example models:** `Quartino_2019_trastuzumab.R` (per-group typical-value switch on linear CL; NSCLC plus a small residual group of prostate, ovarian, and other histologies), `Wang_2024_sugemalimab.R` (heterogeneous solid-tumor residual group of n = 174; exponential coefficient log(0.885) on CL and log(0.926) on Vc; NSCLC is the reference group, not part of `TUMTP_OTH`).
- **Notes:** Scope: specific because the set of histologies collapsed into "Others" is defined by the analysis plan of the source paper; two papers' `TUMTP_OTH` columns are not interchangeable. Document the exact per-paper composition (e.g., "NSCLC + prostate + ovarian + other, n = 107 in Quartino 2019"; "miscellaneous solid tumors excluding NSCLC, lymphoma, GCGEJ, and ESCC, n = 174 in Wang 2024") in `covariateData[[TUMTP_OTH]]$notes`. A given subject can have at most one of the `TUMTP_<GROUP>` indicators (including `TUMTP_OTH`) set to 1; all-zero means the reference group.

### MCPROT (**canonical for serum monoclonal (M) protein concentration**)
- **Description:** Serum monoclonal (M) protein concentration. Multiple-myeloma plasma-cell-burden marker secreted by the tumor clone; elevated MCPROT reflects higher tumor burden and (for IgG-secreting MM) competes with therapeutic IgG mAbs for FcRn-mediated salvage and target-mediated elimination. Typically time-varying — measured at multiple visits over the treatment course and supplied at every PK observation time via linear interpolation between measurements.
- **Units:** g/dL (US-convention; equivalent to 10 g/L SI). Document the unit used in each model via `covariateData[[MCPROT]]$units`.
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — used as a continuous, log-linear effect on the Vmax of target-mediated elimination via `exp(theta * MCPROT)` (i.e., MCPROT enters un-log-transformed). Reference values observed: 0 g/dL (Ide 2020 Vmax,REF reference) and 2.0 g/dL (Ide 2020 figure-1 reference patient).
- **Source aliases:**
  - `TMCPROT` (time-varying serum M-protein concentration) — used in `Ide_2020_elotuzumab.R`. NONMEM column with imputation sentinel `-99` for missing observations, replaced by population median 2.1 g/dL via `IF(TMCPROT.EQ.-99) TMCPROT = 2.1`.
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
  - `PS` / `BPS` — used in `Bajaj_2017_nivolumab.R` (BPS = "baseline performance status"; the one study using Karnofsky Performance Status was mapped to ECOG via Oken 1982 before binarization) and `Zhang_2019_nivolumab.R` (paper's binary collapse PS=0 vs. PS>0).
  - `ECOG_1` — alternative explicit form; equivalent to `ECOG_GE1` when ECOG only takes values 0, 1, 2 in the analysis dataset (the typical oncology case).
  - `ECOG_PS_GT0` — retired name used in earlier register drafts; semantically identical (`>= 1` equals `> 0` for integer ECOG scores).
  - `ECOG101` (categorical 0/1/2 score with thresholding `IF(ECOG101.GT.0.5)`) — used in `Ide_2020_elotuzumab.R`. Decompose: `ECOG_GE1 = as.integer(ECOG101 >= 1)`.
- **Example models:** `Bajaj_2017_nivolumab.R` (exponential effect on CL with coefficient 0.172), `Zhang_2019_nivolumab.R` (exponential effect exp(0.181) on baseline CL; additive effect -0.138 on the time-varying-CL Emax parameter), `Ide_2020_elotuzumab.R` (multiplicative effect on CL = 1.03; paired with `ECOG_GE2` for separate ECOG=1 vs ECOG>=2 effects).
- **Notes:** Oncology papers conventionally report ECOG as an integer (0-5) but binarize at >= 1 because ECOG >= 2 is rare in trial cohorts. When a source paper provides the ordinal ECOG score separately, derive `ECOG_GE1 = as.integer(ECOG >= 1)`. Zhang 2019 uses `ECOG_GE1` on both baseline CL and the time-varying Emax parameter (unlike Bajaj 2017, which uses it on CL only); document the structural role in each model's `covariateData[[ECOG_GE1]]$notes`. When a paper retains separate effects for ECOG = 1 vs ECOG >= 2 (Ide 2020), pair this column with `ECOG_GE2` and supply both indicators in the event dataset.

### ECOG_GE2 (**canonical for Eastern Cooperative Oncology Group performance-status indicator, >= 2**)
- **Description:** 1 if baseline Eastern Cooperative Oncology Group (ECOG) performance status is greater than or equal to 2, 0 if ECOG <= 1. Time-fixed per subject. Used in models that retain separate effects for ECOG = 1 vs ECOG >= 2 by pairing this column with `ECOG_GE1`.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (ECOG performance status <= 1; in models that pair `ECOG_GE1` and `ECOG_GE2`, both indicators = 0 corresponds to ECOG = 0 and (`ECOG_GE1` = 1, `ECOG_GE2` = 0) corresponds to ECOG = 1).
- **Source aliases:**
  - `ECOG101` (categorical 0/1/2 score with thresholding `IF(ECOG101.GT.1.5)`) — used in `Ide_2020_elotuzumab.R`. Decompose: `ECOG_GE2 = as.integer(ECOG101 >= 2)`.
- **Example models:** `Ide_2020_elotuzumab.R` (multiplicative effect on CL = 1.15; paired with `ECOG_GE1` to retain separate ECOG = 1 vs ECOG >= 2 effects).
- **Notes:** Parallels `ECOG_GE1`. Use only when the source paper reports a separate effect for ECOG >= 2 in addition to ECOG_GE1; otherwise `ECOG_GE1` alone is sufficient. The paired (`ECOG_GE1`, `ECOG_GE2`) decomposition reproduces a three-level (`ECOG = 0`, `ECOG = 1`, `ECOG >= 2`) ordinal effect with two binaries.

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

### TUMTP_LYMPH (**canonical for lymphoma (pooled) tumor-type indicator**)
- **Description:** 1 = lymphoma (heterogeneous lymphoma pool spanning multiple lymphoma histologies — e.g., classical Hodgkin lymphoma combined with extranodal NK/T-cell lymphoma), 0 = solid tumor or other tumor type.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = non-lymphoma tumor type (per source paper; e.g., NSCLC, GC/GEJ, ESCC, "Other" solid tumors in the Wang 2024 cohort, with NSCLC as the implicit reference when paired with the other Wang 2024 `TUMTP_*` indicators).
- **Source aliases:**
  - `TTYPE1` (Wang 2024) — decompose into `TUMTP_LYMPH = as.integer(TTYPE1 == 1)`. The Wang 2024 source paper uses a multi-level `TTYPE` factor with levels 1 = lymphoma, 2 = lung cancer (reference), 3 = other, 4 = GCGEJ, 5 = ESCC.
- **Example models:** `Wang_2024_sugemalimab.R` (exponential coefficient log(0.877) on baseline CL and log(0.879) on Vc).
- **Notes:** Distinct from `TUMTP_CHL` (which is specifically classical Hodgkin lymphoma). Wang 2024 pools two lymphoma histologies (extranodal NK/T-cell lymphoma from CS1001-201 / NCT03595657 and classical Hodgkin lymphoma from CS1001-202 / NCT03505996) into a single lymphoma indicator; the indicator therefore captures a generic "hematologic-vs-solid-tumor" contrast rather than a histology-specific effect. When a future paper studies a single lymphoma histology distinct from cHL, register a more specific canonical (e.g., `TUMTP_ENKTL`, `TUMTP_NHL`) rather than overloading this one. Document the per-paper histology composition in `covariateData[[TUMTP_LYMPH]]$notes`.

### TUMTP_ESCC (**canonical for oesophageal-squamous-cell-carcinoma tumor-type indicator**)
- **Description:** 1 = oesophageal squamous cell carcinoma (ESCC), 0 = other tumor types.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 = all other tumor types (per source paper; in Wang 2024, the implicit reference is NSCLC when all the other `TUMTP_*` indicators are also 0).
- **Source aliases:**
  - `TTYPE5` (Wang 2024) — decompose into `TUMTP_ESCC = as.integer(TTYPE5 == 1)`.
- **Example models:** `Wang_2024_sugemalimab.R` (exponential coefficient log(0.99) on baseline CL and log(1.08) on Vc).
- **Notes:** Follows the `TUMTP_CHL` / `TUMTP_GC` / `TUMTP_SCLC` decomposition pattern. Distinct from gastroesophageal-junction adenocarcinoma (which is captured by the broader `TUMTP_GC` indicator that pools GC and GEJ adenocarcinomas) — ESCC is a squamous-cell histology, not adenocarcinoma. Document the per-paper histology composition in `covariateData[[TUMTP_ESCC]]$notes`.

### LINE_1L (**canonical for first-line-therapy indicator**)
- **Description:** 1 = first-line therapy (1L) / treatment-naive, 0 = second-line or greater (2L+) / relapsed-or-refractory.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (2L+, second-line or greater / relapsed-refractory).
- **Source aliases:**
  - `LINE` (categorical column with levels `1L`, `2L`, `3L+`, ...) — decompose into `LINE_1L = as.integer(LINE == "1L")`.
  - `RRFN` (relapsed/refractory flag; treatment-naive corresponds to RRFN == 0) — used in `Lu_2019_polatuzumab.R`. Decompose: `LINE_1L = as.integer(RRFN == 0)`.
- **Example models:** `Sanghavi_2020_ipilimumab.R` (exponential coefficient -0.0949 on CL), `Lu_2019_polatuzumab.R` (multiplicative effects on V1 = 1.20, kdes = 3.38, CL_T = 3.53, FRAC_NS = 0.756; the same pooled-trial NHL cohort mixes 415 R/R and 45 first-line patients).
- **Notes:** Promoted to scope: general on 2026-04-26 after Lu 2019 polatuzumab vedotin ratified the same 1L vs 2L+ binarization that Sanghavi 2020 ipilimumab introduced. The two papers use different indicator semantics (Sanghavi reports the effect as `exp(-0.0949 * LINE_1L)` and Lu reports `theta^LINE_1L` with theta < or > 1 depending on the parameter); both reduce to the same canonical 0/1 column. If a future paper requires finer resolution (separate effects for 2L vs 3L+), add a parallel `LINE_2L` canonical rather than overloading this one.

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
### COMBO_RG (**canonical for anti-CD20 (rituximab or obinutuzumab) combination-therapy indicator**)
- **Description:** 1 = polatuzumab vedotin co-administered with rituximab OR obinutuzumab, 0 = single-agent polatuzumab vedotin (or any other regimen lacking an anti-CD20 partner). Time-fixed per subject in the source paper's analysis cohort.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (single-agent polatuzumab vedotin or no anti-CD20 partner).
- **Source aliases:**
  - `COMBO` (categorical: 0 = single agent, 1 = + rituximab, 2 = + obinutuzumab) — used in `Lu_2019_polatuzumab.R`. Decompose: `COMBO_RG = as.integer(COMBO == 1 | COMBO == 2)`. The Lu 2019 NONMEM separately defines `RTX = as.integer(COMBO == 1)` and `GA101 = as.integer(COMBO == 2)` and applies effects as `theta^(RTX + GA101)`; because RTX and GA101 are mutually exclusive, RTX + GA101 takes values {0, 1} and the effect collapses to `theta^COMBO_RG`.
- **Example models:** `Lu_2019_polatuzumab.R` (multiplicative effects on CL_INF = 0.844, kdes = 0.932, FRAC_NS = 0.709).
- **Notes:** Rituximab and obinutuzumab both bind CD20 on B cells (rituximab is a Type I anti-CD20 mAb, obinutuzumab a glycoengineered Type II), so co-administration is hypothesized to alter polatuzumab vedotin disposition through depletion of CD79b+ target B cells. The Lu 2019 final model fits a single combined effect rather than separate rituximab- and obinutuzumab-specific effects. Scope: specific because the relevant combination partners (CD20-directed mAbs) are tied to NHL pathway; if a future paper distinguishes rituximab from obinutuzumab combinations, register `COMBO_R` and `COMBO_G` separately rather than overloading this canonical.
### COMBO_DURVA (**canonical for durvalumab combination-therapy indicator**)
- **Description:** 1 = the analyzed therapeutic mAb is co-administered with durvalumab (anti-PD-L1 IgG1), 0 = monotherapy.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (monotherapy).
- **Source aliases:**
  - `COMB` — used in `Hwang_2022_tremelimumab.R` ($INPUT NM-TRAN data item; control-stream switch `IF(COMB.EQ.0)` selects monotherapy parameters and `IF(COMB.EQ.1)` selects combination-therapy parameters).
- **Example models:** `Hwang_2022_tremelimumab.R` (selects between monotherapy and combination-with-durvalumab values of the time-varying-CL Tmax and lambda parameters).
- **Notes:** Parallels `COMBO_NIVO` but for durvalumab rather than nivolumab co-administration. Acts on the time-varying-CL component (Tmax and lambda); baseline CL is shared between monotherapy and combination groups in Hwang 2022.
### COMBO_LEN_DEX (**canonical for lenalidomide plus dexamethasone combination-therapy indicator**)
- **Description:** 1 = the analyzed therapeutic mAb (or other agent under PK study) is co-administered with the lenalidomide + low-dose-dexamethasone (Ld) backbone, 0 = monotherapy or any non-Ld regimen. Lenalidomide is an immunomodulatory imide (IMiD) that activates natural killer cells; dexamethasone is an immunosuppressant glucocorticoid. The Ld backbone is a standard combination partner in multiple-myeloma and other hematologic-malignancy regimens.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (monotherapy or any non-Ld regimen). When a paper reports the reference patient as "with Ld coadministration" (as Ide 2020 does), the model still stores the canonical 0/1 column and applies the effect as `exp(theta * (COMBO_LEN_DEX - 1))` so that COMBO_LEN_DEX = 1 yields factor 1 (paper's reference) and COMBO_LEN_DEX = 0 activates the effect.
- **Source aliases:**
  - `LENDEX` (1 = with Ld, 0 = without Ld; in Ide 2020 derived from `STUDY != 204011` because study 204011 was the Ld-free elotuzumab-monotherapy cohort) — used in `Ide_2020_elotuzumab.R`.
  - `COMBO_LD` (**retired** canonical name; renamed to `COMBO_LEN_DEX` on 2026-04-27 for clarity).
- **Example models:** `Ide_2020_elotuzumab.R` (multiplicative effects: CLLd = 0.74 on nonspecific CL, encoded as `exp(log(0.74) * (COMBO_LEN_DEX - 1))`; KINTLd = 10.1 on the second-order target-mediated elimination rate from the peripheral compartment, encoded as `exp(log(10.1) * (COMBO_LEN_DEX - 1))`).
- **Notes:** Specific scope because the canonical's mechanistic relevance is hematologic-malignancy-domain-bound (multiple myeloma and related plasma-cell or B-cell disorders). Distinct from `COMBO_BELAMAF` (which pools Ld with bortezomib-dex and pomalidomide-dex into a single broader "any-combination" belantamab indicator); `COMBO_LEN_DEX` is the per-backbone Ld-only flag. If a future paper distinguishes "Ld-only" from a broader "any-IMiD-plus-dex" backbone with separate effects, register a parallel canonical (e.g., `COMBO_PD` for pomalidomide-dex, `COMBO_VD` for bortezomib-dex). Sign of the exponential coefficient is paper-dependent: Ide 2020 reports `CLLd = 0.74` so the Ld-coadministration arm has 26% lower nonspecific CL than the Ld-free arm, but `KINTLd = 10.1` so the Ld arm has 10x higher second-order target-mediated elimination — both are mechanistically interpretable (dexamethasone suppresses non-specific catabolic clearance; lenalidomide-activated NK cells increase target-cell-binding-mediated elimination).

### COMBO_BELAMAF (**canonical for any-combination belantamab mafodotin therapy indicator**)
- **Description:** 1 = belantamab mafodotin administered as part of a combination regimen (with bortezomib + dexamethasone, lenalidomide + dexamethasone, or pomalidomide + dexamethasone) in the relapsed/refractory multiple myeloma setting; 0 = belantamab mafodotin monotherapy.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (belantamab mafodotin monotherapy).
- **Source aliases:**
  - `COMBO` (when the source dataset uses a generic combination flag for belantamab mafodotin pooled regimens) — used in `Papathanasiou_2025_belantamab.R`.
- **Example models:** `Papathanasiou_2025_belantamab.R` (multiplicative factor θ = 1.44 on the Imax parameter of the time-varying CL function — combination therapy increases the steady-state CL reduction from 33.2 % to 44.0 %).
- **Notes:** Pools the three combination backbones tested in DREAMM-6 / DREAMM-7 / DREAMM-8 (Bor-Dex, Len-Dex, Pom-Dex) into a single binary because Papathanasiou 2025 reports no meaningful per-backbone difference in cycle-1 ADC exposure. If a future paper tests per-backbone combination effects, register dedicated indicators (`COMBO_BELAMAF_BORDEX`, etc.) rather than overloading this aggregate.

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

### MAYO_E (**canonical for baseline Mayo endoscopic subscore**)
- **Description:** Mayo endoscopic subscore at baseline for ulcerative colitis, integer 0-3 (0 = normal / inactive disease, 1 = mild, 2 = moderate, 3 = severe). The endoscopic subscore is one of the four components of the full Mayo score (0-12); the partial Mayo score (`PMAYO`) excludes it. Time-fixed per subject.
- **Units:** (score, 0-3)
- **Type:** categorical
- **Scope:** general
- **Reference category:** depends on per-model encoding — papers that include the full 0-3 range typically reference Mayo 0 or 1 (mild), while papers restricted to moderate-to-severe UC (the typical biologic-induction-therapy population) reference Mayo 2 or 3. Document the per-model reference category in `covariateData[[MAYO_E]]$reference_category`.
- **Source aliases:**
  - `MPRE` — used in `Faelens_2021_infliximab.R` (NONMEM column for "Mayo endoscopic score pre-induction"). The Faelens 2021 dataset additionally codes a "missing" sentinel `MPRE = -99`; treat this as out-of-domain when applying the model and document per-model.
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

### RHEUMATOID_FACTOR (**canonical for serum rheumatoid factor concentration**)
- **Description:** Serum rheumatoid factor (an autoantibody, predominantly IgM, directed against the Fc portion of IgG) concentration. Baseline value typical; document time-varying use in per-model `notes`.
- **Units:** U/mL or IU/mL (interchangeable in the clinical-PK literature). Document per-model via `covariateData[[RHEUMATOID_FACTOR]]$units`.
- **Type:** continuous
- **Scope:** general
- **Reference category:** n/a — used with power scaling on a log-transformed value: `(log(RHEUMATOID_FACTOR) / log(ref))^exponent` (or, equivalently, the source-paper form `(LRF / log(ref))^exponent` where `LRF = log(RHEUMATOID_FACTOR)`). Reference value observed: 110 U/mL (Frey 2010, corresponding to LRF = 4.7 in the paper's final-model equation).
- **Source aliases:**
  - `LRF` — log-transformed RF (natural log of the value in U/mL); Frey 2010 fits the covariate on the log scale and reports the reference as `LRF = 4.7` (i.e., `log(110) ≈ 4.7`). The canonical column carries the raw RF concentration in U/mL; the log transform is applied inside `model()`.
  - `RF` — universal NONMEM/clinical-PK abbreviation; rejected as the canonical name on 2026-04-28 because the bare two-letter abbreviation is uncommon in published popPK papers and could be confused with other shortenings.
- **Example models:** `Frey_2010_tocilizumab.R` (U/mL, reference 110 U/mL ≡ LRF = 4.7; small positive exponent +0.1 on linear CL applied to `log(RHEUMATOID_FACTOR)`).
- **Notes:** RF concentrations span several orders of magnitude across the rheumatoid-arthritis population (Frey 2010 observed range 15-11,800 U/mL across the four phase-III studies; reference paper: Frey 2010 Table I), motivating the log transform before power scaling. The mechanistic rationale (Frey 2010 Discussion, p764) is that RF — being an anti-IgG autoantibody — could in principle bind the Fc region of the therapeutic IgG monoclonal antibody and accelerate clearance, but the observed CL effect was small (-4.9% to +6.5% across the observed RF range) and the paper acknowledges that high RF concentrations may also reduce the assay's ability to detect the drug, leading to an apparent CL increase. Ratified canonically on 2026-04-28 alongside the Frey 2010 extraction.

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

### COADMIN_CHEMO (**canonical for anti-PD-(L)1 mAb + chemotherapy combination indicator**)
- **Description:** 1 = subject is receiving an anti-PD-(L)1 monoclonal antibody in combination with platinum-based chemotherapy (gemcitabine + cisplatin, pemetrexed + cisplatin, paclitaxel + carboplatin, or platinum-doublet); 0 = otherwise. Encodes chemotherapy coadministration as a study-design covariate on the antibody's CL.
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0 (no chemotherapy coadministration — monotherapy or, where applicable, a non-chemotherapy combination such as anti-PD-1 + anti-CTLA-4).
- **Source aliases:**
  - `CHEMO` — used in `Zhang_2019_nivolumab.R`.
  - `MONOTR` — used in `Kuchimanchi_2024_dostarlimab.R` (the paper's structural-equation indicator for monotherapy; the canonical column carries the inverse value, i.e. `COADMIN_CHEMO = 1 - MONOTR`, so the canonical column is 1 when the patient is on combo-chemotherapy).
- **Example models:** `Zhang_2019_nivolumab.R` (exponential effect on baseline CL: `exp(-0.104)` ≈ 0.90 fold, i.e. ~9.7% lower CL relative to monotherapy), `Kuchimanchi_2024_dostarlimab.R` (multiplicative effect on baseline CL: `1 - 0.0779` = 0.922, i.e. 7.79% lower CL on dostarlimab + carboplatin/paclitaxel relative to dostarlimab monotherapy).
- **Notes:** Promoted from specific to general scope on 2026-04-27 after the Kuchimanchi 2024 dostarlimab + carboplatin/paclitaxel analysis ratified the same pooling convention (any chemotherapy backbone collapsed into a single binary indicator). The two papers use different functional forms for the effect on CL — Zhang 2019 uses `exp(theta * COADMIN_CHEMO)` (exponential) and Kuchimanchi 2024 uses `(1 + theta * COADMIN_CHEMO)` (multiplicative); these are different parameterisations of the same underlying study-design indicator and the canonical column meaning is unchanged. Document the per-model functional form in `covariateData[[COADMIN_CHEMO]]$notes`.

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

### COADMIN_AVD (**canonical for brentuximab vedotin + AVD (adriamycin/doxorubicin, vinblastine, dacarbazine) combination indicator**)
- **Description:** 1 = subject is receiving brentuximab vedotin in combination with the AVD chemotherapy backbone (adriamycin a.k.a. doxorubicin, vinblastine, dacarbazine) for newly diagnosed advanced-stage Hodgkin lymphoma; 0 = otherwise (single-agent brentuximab vedotin). Encodes the A+AVD frontline regimen as a study-design covariate on ADC clearance.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (single-agent brentuximab vedotin — no AVD coadministration).
- **Source aliases:**
  - `DOX` — used in `Zhou_2025_brentuximab.R` (the NONMEM dataset uses the doxorubicin-administration flag as the AVD-coadministration indicator because doxorubicin is given on the same days as the other AVD agents in this regimen).
- **Example models:** `Zhou_2025_brentuximab.R` (power-form effect on ADC clearance: `CL * 2.12^COADMIN_AVD` — ADC clearance is ~2.1-fold higher under A+AVD vs single-agent BV).
- **Notes:** Distinct from `COADMIN_CHEMO` (which is nivolumab + platinum-based chemotherapy). The A+AVD regimen is the standard chemotherapy backbone for frontline classical Hodgkin lymphoma; promote to general scope if a second BV paper reports the same A+AVD covariate effect with a comparable encoding.

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
  - `NAB` (neutralizing antibody positive — used in `Petrov_2024_romiplostim.R`). Strictly a subset of total ADA-positive (ADA antibodies that neutralize the drug's biological effect). Document per-model when the source assay measured NAB only and the canonical column thus excludes binding-only ADA.
- **Example models:** `Clegg_2024_nirsevimab.R`, `Hu_2026_clesrovimab.R`, `Petrov_2024_romiplostim.R`, `Xu_2019_sarilumab.R`.

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
- **Description:** 1 = extensive colitis or pancolitis disease extension, 0 = otherwise (any non-extensive disease extension, e.g. left-sided colitis or proctitis).
- **Units:** (binary)
- **Type:** binary
- **Scope:** general
- **Reference category:** 0. In papers that decompose the disease-extension categorical into both `DISEXT_EP` and `DISEXT_OTHER`, `DISEXT_EP = 0 AND DISEXT_OTHER = 0` corresponds to the left-sided-colitis reference group; in papers that use a single binary indicator for extensive colitis, the reference is pooled non-extensive (left-sided + any other extension).
- **Source aliases:**
  - `EXTCOL` — used in `Faelens_2021_infliximab.R` (binary 0/1 for extensive colitis at baseline; no separate "other" category).
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
- **Units:** (paper-specified — typically dimensionless score)
- **Type:** continuous
- **Scope:** specific
- **Reference category:** n/a — usually used with a linear-deviation form `1 + e * (PAIN - ref)` or a power form `(PAIN / ref)^e`. Document the reference value in `covariateData[[PAIN]]$notes`.
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

### FORM_P2F2
- **Description:** 1 = isatuximab P2F2 drug material (intended commercial / phase III material, used in the EFC14335 / ICARIA-MM study), 0 = P1F1 drug material (early-phase material used in TED10893 / TED14154 / TCD14079).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (P1F1).
- **Source aliases:**
  - `Drug_mat` — used in `Fau_2020_isatuximab.R`. Values 0 / 1 with the same orientation as the canonical (1 = P2F2 / commercial-bound material).
- **Example models:** `Fau_2020_isatuximab.R` (exponential effect on Vc with coefficient -0.137; P2F2 patients had ~13% lower Vc than P1F1).
- **Notes:** Phase III / commercial-bound formulation indicator for isatuximab; the FORM_* family stays scope-specific per nlmixr2lib policy that drug-product-version indicators are kept model-specific unless they generalize across multiple drugs. Set to 1 to simulate the marketed material.

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
- **Example models:** `Othman_2014_daclizumab.R`, `Diao_2016_daclizumab_cd25.R`, `Diao_2016_daclizumab_cd56bright.R`, `Diao_2016_daclizumab_treg.R`.
- **Notes:** Othman 2014 estimated two separate absolute bioavailabilities because of non-linear dose-normalized exposure at the 50 mg SC dose — F = 0.84 for the therapeutic 100-300 mg SC range and F = 0.57 for the 50 mg SC cohort. Encoded here as a record-level indicator so the covariate effect `e_dose_50mg_f = 0.57/0.84 - 1 = -0.321` scales bioavailability only on 50 mg SC doses. For clinical-range simulation (150 mg SC Q4W Phase III regimen) leave `DOSE_50MG = 0`. The Diao 2016 PK/PD models inherit the Othman 2014 PK backbone verbatim; they carry `DOSE_50MG` even though the Diao 2016 RRMS regimens are 150 / 300 mg SC only.

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

### STUDY_ABA2_HLA78
- **Description:** 1 = subject enrolled in the ABA2 hematopoietic-cell-transplant trial (IM101311; NCT01743131) HLA 7/8 (one-allele-mismatched donor) cohort, 0 = any other study in the Takahashi 2023 pooled abatacept population PK analysis (RA/JIA reference and ABA2 8/8 cohort).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-ABA2-7/8 — pooled adult RA / pediatric JIA cohort, plus ABA2 HLA 8/8 cohort which is flagged separately by `STUDY_ABA2_HLA88`).
- **Source aliases:** derived per subject from the trial-cohort identifier (`Cohort = ABA2 7/8` → 1, else → 0). Takahashi 2023 Supplemental Table 4 names the corresponding theta `θCohort_CL` / `θCohort_Vc`.
- **Example models:** `Takahashi_2023_abatacept.R` (multiplicative `Ratio` factors on CL = 0.70 and on Vc = 0.99 vs the RA/JIA reference; values from Takahashi 2023 Supplemental Table 4).
- **Notes:** Pairs with `STUDY_ABA2_HLA88` to reproduce the three-level cohort categorical (RA/JIA, ABA2 7/8, ABA2 8/8) the paper reports as the only retained categorical PK covariate. At most one of `STUDY_ABA2_HLA78` and `STUDY_ABA2_HLA88` is 1 per subject; both 0 reproduces the RA/JIA reference. Scope: specific because the contrast is tied to the ABA2-vs-RA/JIA pooling design.

### STUDY_ABA2_HLA88
- **Description:** 1 = subject enrolled in the ABA2 hematopoietic-cell-transplant trial (IM101311; NCT01743131) HLA 8/8 (allele-matched donor) cohort, 0 = any other study in the Takahashi 2023 pooled abatacept population PK analysis (RA/JIA reference and ABA2 7/8 cohort).
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-ABA2-8/8 — pooled adult RA / pediatric JIA cohort, plus ABA2 HLA 7/8 cohort which is flagged separately by `STUDY_ABA2_HLA78`).
- **Source aliases:** derived per subject from the trial-cohort identifier (`Cohort = ABA2 8/8` → 1, else → 0). Takahashi 2023 Supplemental Table 4 names the corresponding theta `θCohort_CL` / `θCohort_Vc`.
- **Example models:** `Takahashi_2023_abatacept.R` (multiplicative `Ratio` factors on CL = 0.91 and on Vc = 1.32 vs the RA/JIA reference; values from Takahashi 2023 Supplemental Table 4).
- **Notes:** Pairs with `STUDY_ABA2_HLA78`. At most one of the two indicators is 1 per subject; both 0 reproduces the RA/JIA reference. Scope: specific.
## Study-site region

Geographical study-site region indicators. Distinct from race / ethnicity (`RACE_*`), which describe subject ancestry; these describe the geographical location of the clinical trial site that enrolled the subject. Used in multi-regional studies (typically those including bridging analyses for Japan or East Asia) to capture region-specific clinical-practice or unmeasured-environment effects on PK that remain after accounting for body weight, race, and laboratory covariates. Encoded as a set of mutually exclusive binary indicators with US as the implicit reference category (all indicators = 0). When a paper groups some non-US regions with US (e.g., Hong 2025 groups US and Japan as the DXd CL reference), the model code uses only the indicators that distinguish the non-reference groups; the data column for the grouped region (e.g., `REGION_JAPAN`) is still recorded so the same dataset can serve other parameters that do separate that group.

### REGION_JAPAN (**canonical for Japan study-site / enrollment-country indicator**)
- **Description:** 1 = study site in Japan (or patient enrolled in Japan, depending on source paper's reporting), 0 = study site / enrollment country outside Japan. Geographical Japan indicator.
- **Units:** (binary)
- **Type:** binary
- **Scope:** specific
- **Reference category:** 0 (non-Japan sites / enrollment; specific reference set varies per model — e.g., Hong 2025 Dato-DXd CL uses "any non-Japan region", whereas Hong 2025 DXd CL groups Japan with US into the reference).
- **Source aliases:**
  - `COUNTRY_JPN` — retired canonical; used in `Yin_2021_trastuzumabDeruxtecan.R` as an enrollment-country (not study-site-region) indicator. Some papers report country-of-enrollment rather than site region; both map to `REGION_JAPAN` when the binary contrast is Japan vs. non-Japan.
- **Example models:** `Hong_2025_datopotamab.R` (multiplicative effect 1 + (−0.219) = 0.781 on Dato-DXd linear clearance), `Yin_2021_trastuzumabDeruxtecan.R` (multiplicative effect 0.903 on CL_intact and 0.738 on V2_intact when REGION_JAPAN = 1; Yin 2021 retained Japan enrollment-country over Japanese race because the two were highly confounded, correlation −0.81).
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
- **Reference category:** n/a — used either with a power-covariate form `CYCLE^Fm` (Fm typically negative) to capture cycle-over-cycle decline in a derived quantity such as ADC-to-payload conversion fraction (Li 2017 brentuximab vedotin), or with a piecewise indicator `CYCLE == 1 vs CYCLE >= 2` to capture a step change in DAR scaling between cycle 1 and later cycles (Hong 2025 datopotamab deruxtecan).
- **Source aliases:** `CYCLE` — used in `Li_2017_brentuximab.R` and `Hong_2025_datopotamab.R` with the same canonical name.
- **Example models:**
  - `Li_2017_brentuximab.R` (exponent on the fraction of ADC that converts to MMAE by proteolytic degradation, Fm = -0.261, to reflect tumor-burden reduction across successive treatment cycles).
  - `Hong_2025_datopotamab.R` (cycle-1 vs cycle-2+ piecewise scaling Factor1 = 0.696 on the DAR equation that drives DXd formation rate from total Dato-DXd elimination).
- **Notes:** Must be >= 1 throughout (`CYCLE^Fm` is undefined at 0; the piecewise form requires `CYCLE` to be a positive integer at every observation row). Distinct from `ooc<n>` binary-occasion indicators: `CYCLE` is an integer count, not a mutually-exclusive set of indicator columns. Data-assembly helper: set `CYCLE = floor((TIME - TIME_FIRST_DOSE) / cycle_length_days) + 1` for a fixed-interval dosing regimen.

---

## Change log

- **2026-04-28** — Extended `RACE_WHITE` (general scope) example-models list and source aliases to record `Hu_2014_bapineuzumab.R` (Caucasian-vs-non-Caucasian dichotomy with the Caucasian subgroup as the typical-value reference). The canonical column encoding (1 = White / 0 = non-White) is unchanged; the model implements the 15% non-Caucasian effect on `(1 - RACE_WHITE)`. The change log notes that the typical-value reference category may legitimately differ between papers using `RACE_WHITE`.
- **2026-04-27** — Added `STUDY_ABA2_HLA78` and `STUDY_ABA2_HLA88` (both specific scope) canonical entries while extracting `Takahashi_2023_abatacept.R`. The two binary indicators jointly reproduce the three-level RA/JIA-vs-ABA2-7/8-vs-ABA2-8/8 cohort categorical that Takahashi 2023 Supplemental Table 4 retains as the only categorical PK covariate (multiplicative `Ratio` thetas on CL and on Vc; RA/JIA = both indicators 0 = ratio 1 fixed).
- **2026-04-27** — Added a new `## Study-site region` section with `REGION_JAPAN`, `REGION_EUROPE`, `REGION_ROW` canonical entries (all scope: specific, binary indicators) while extracting `Hong_2025_datopotamab.R`. US is the implicit reference category (all three indicators = 0). The new entries are distinct from the existing `RACE_*` family because they encode trial-site geography rather than subject ancestry; a Japanese-ancestry subject enrolled at a US site has `RACE_JAPANESE = 1` and `REGION_JAPAN = 0`.
- **2026-04-26** — Added `B2M` (general-scope serum beta-2-microglobulin under `Renal / hepatic function`; reference 3.90 mg/L from the multiple-myeloma cohort median), `MM_NIGG` (specific-scope non-IgG-MM-vs-IgG-MM within-disease immunoglobulin-type indicator under `Oncology`), and `FORM_P2F2` (specific-scope isatuximab phase III / commercial-bound drug-material indicator, placed alongside the existing `FORM_*` entries) canonical entries while extracting `Fau_2020_isatuximab.R`. Source aliases mapped: `Ig_type`→`MM_NIGG`, `Drug_mat`→`FORM_P2F2`. `Fau_2020_isatuximab.R` added to the `WT`, `SEXF`, and `RACE_ASIAN` example-model lists.
- **2026-04-27** — Promoted `COADMIN_CHEMO` from specific to general scope after the Kuchimanchi 2024 dostarlimab + carboplatin/paclitaxel population-PK analysis ratified the same pooling convention (any chemotherapy backbone collapsed into a single binary indicator) on a second mAb. The two example models use different functional forms — Zhang 2019 uses `exp(theta * COADMIN_CHEMO)` and Kuchimanchi 2024 uses `(1 + theta * COADMIN_CHEMO)` — but the canonical column meaning is unchanged; the per-model functional form is documented in the model files' `notes`. Source-alias list updated: `MONOTR` (Kuchimanchi 2024) added with the inverse-value relation `COADMIN_CHEMO = 1 - MONOTR`. Description and example list extended to include `Kuchimanchi_2024_dostarlimab.R`.
- **2026-04-27** — Added `ACUTE_MED_DAYS` (specific scope under `Disease severity scores`; baseline number of days/month of acute migraine medication use, piecewise-linear with breakpoint at 5 d/mo per medication-overuse-headache convention) and `CAV` (specific scope under a new `Drug exposure metrics` H2 section; average modelled-drug plasma concentration over a dosing interval, used as the exposure covariate in PD exposure-response models that consume empirical-Bayes Cav from a previously-published popPK model) canonical entries while extracting the two Fiedler-Kelly 2020 fremanezumab exposure-response models (`FiedlerKelly_2020_fremanezumab_em.R`, `FiedlerKelly_2020_fremanezumab_cm.R`).
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
- **2026-04-27** — Added `COMBO_DURVA` (specific scope) canonical entry while extracting `Hwang_2022_tremelimumab.R`. Parallels `COMBO_NIVO` (Sanghavi 2020) but for durvalumab co-administration; selects between monotherapy and combination-with-durvalumab values of the time-varying-CL Tmax and lambda parameters. Source alias `COMB` (NM-TRAN $INPUT data item) mapped.
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
- **2026-04-27** — Added `MAYO_E` (general-scope baseline Mayo endoscopic subscore, integer 0-3) canonical entry under `Inflammatory-bowel-disease disease-activity covariates` while extracting `Faelens_2021_infliximab.R`. Source alias `MPRE` mapped (Faelens 2021 NONMEM column for "Mayo endoscopic score pre-induction"; the dataset's `MPRE = -99` missing sentinel is documented as out-of-domain in the per-model `notes`). Distinct from the full Mayo score (0-12) and partial Mayo `PMAYO` (0-9). Also promoted `DISEXT_EP` from scope: specific to scope: general (the binary "extensive colitis vs not" semantics generalize across UC popPK papers regardless of whether the source dataset also distinguishes an "other" disease-extension category) and added Faelens 2021 source alias `EXTCOL` and example-model entry.
- **2026-04-26** — Added `HEP_IMP` (general-scope binary indicator for NCI ODWG hepatic impairment, mild or worse vs. normal) under `Renal / hepatic function` and `COMBO_RG` (specific-scope binary indicator for anti-CD20 (rituximab or obinutuzumab) combination therapy) under `Oncology` while extracting `Lu_2019_polatuzumab.R`. Source aliases mapped: `BHPTGRPN` (categorical NCI ODWG group with 9999 missing-value sentinel)→`HEP_IMP`; `COMBO` (Lu 2019 categorical 0/1/2)→`COMBO_RG`.
- **2026-04-27** — Added `SBCMA` (specific-scope, baseline soluble B-cell maturation antigen, in `Cardiometabolic / target biomarkers`; reference 50 ng/mL) and `COMBO_BELAMAF` (specific-scope, any-combination belantamab mafodotin therapy indicator on Imax of the time-varying CL function, in `Concomitant / prior medication`) canonical entries while extracting `Papathanasiou_2025_belantamab.R`. Source aliases mapped: `SBCMABL`→`SBCMA`, `COMBO`→`COMBO_BELAMAF`.
- **2026-04-28** — Added `CSF1` (specific-scope plasma colony-stimulating factor 1 / M-CSF concentration; placed under `Inflammation markers`), `CPK` (general-scope serum creatine phosphokinase / creatine kinase; placed under `Renal / hepatic function` next to AST/ALT/LDH although mechanistically a muscle-origin enzyme), and `DIS_CANCER` (specific-scope advanced-solid-tumor cohort indicator; placed under `Disease state (cross-population indicators)`) canonical entries while extracting `Yang_2024_axatilimab.R`. Extended `DIS_HV` example_models with Yang 2024 as the second user (cGVHD-reference complement to Nikanjam 2019's non-HV-oncology reference). Source aliases mapped: `BLCSF1` / `BL_CSF1` (Monolix model-parameter name) → `CSF1`; `BLCPK` → `CPK`. Reference values: 549 pg/mL (CSF1), 63 U/L (CPK), pooled-cohort medians from Yang 2024 Table S3.
- **2026-04-28** — Added `DIS_SAD` (specific-scope binary indicator partitioning hypogammaglobulinaemia patients by mechanism: 1 = secondary antibody deficiency, 0 = primary immunodeficiency; placed under `Disease state (cross-population indicators)`) and `IGM` (specific-scope serum immunoglobulin M concentration as a B-cell humoral-capacity proxy; placed under `Renal / hepatic function` next to `IGG`) canonical entries while extracting `Cheng_2026_immunoglobulin.R`. Source aliases mapped: none (the source paper uses bare prose names "type of immunodeficiency" and "IgM level"). Reference values: 0.21 g/L (IGM), Cheng 2026 pooled-cohort median.
- **2026-04-28** — Added `HDLC` (general-scope serum HDL cholesterol; placed under `Cardiometabolic / target biomarkers`), `TPRO` (general-scope total serum protein; placed under `Renal / hepatic function` next to `ALB`), and `RHEUMATOID_FACTOR` (general-scope serum rheumatoid factor concentration with log-transform-then-power scaling; placed under `Rheumatoid-arthritis disease-activity covariates`) canonical entries while extracting `Frey_2010_tocilizumab.R`. Source aliases mapped: `HDL-C` / `HDL_C` → `HDLC`; `PROT` / `TP` → `TPRO`; `LRF` (Frey 2010 final-equation log-RF variable) and `RF` (universal NONMEM/clinical-PK abbreviation, but explicitly rejected on 2026-04-28 as the canonical name in favour of the unambiguous `RHEUMATOID_FACTOR`) → `RHEUMATOID_FACTOR`. Reference values: 54 mg/dL (HDLC), 74 g/L (TPRO), 110 U/mL ≡ LRF = 4.7 (RHEUMATOID_FACTOR), all from the Frey 2010 pooled-cohort medians.
- **2026-04-28** — Added `IGE` (general-scope, baseline serum total immunoglobulin E concentration, in `Cardiometabolic / target biomarkers`; reference 482.4 ng/mL with optional `IU/mL` reporting via `1 IU/mL = 2.42 ng/mL` documented per-model) canonical entry while extracting `Hayashi_2007_omalizumab.R`. Source alias mapped: `IgE0`→`IGE`. Distinguished in the H3 notes from the in-model dynamic `X_TE` IgE state used by mechanism-based binding/turnover anti-IgE models.
- **2026-04-28** — Added `ECOG_GE2` (general-scope ECOG performance-status >= 2 indicator; pairs with `ECOG_GE1` to retain three-level ECOG = 0 / 1 / >=2 ordinal effects in models that test separate >=1 and >=2 thresholds; placed under `Oncology` directly after `ECOG_GE1`), `MCPROT` (specific-scope serum monoclonal protein concentration; tumor-burden marker for plasma-cell-targeting therapies in multiple myeloma; placed under `Oncology` near other tumor-burden markers; reference values 0 g/dL paper-Vmax-reference / 2.0 g/dL paper-figure-reference; time-varying), and `COMBO_LEN_DEX` (specific-scope lenalidomide + dexamethasone combination-therapy indicator on hematologic-malignancy backbone; placed under `Concomitant / prior medication`; distinct from `COMBO_BELAMAF` which pools Ld with bortezomib-dex and pomalidomide-dex) canonical entries while extracting `Ide_2020_elotuzumab.R`. Source aliases mapped: `ECOG101` (with thresholding `IF(.GT.0.5)` for `ECOG_GE1` and `IF(.GT.1.5)` for `ECOG_GE2`)→both `ECOG_GE1` and `ECOG_GE2`; `TMCPROT`→`MCPROT`; `LENDEX`→`COMBO_LEN_DEX`. Extended `ECOG_GE1` example-models list with `Ide_2020_elotuzumab.R`. The Ide 2020 model also extends the example-model lists for `WT`, `AGE`, `SEXF`, `RACE_ASIAN`, `CRCL` (eGFR), `LDH`, `ALB`, `B2M`, `HEPIMP`, and `LINE_1L`.
- **2026-04-28** — Added `DIS_AD` (specific-scope, Alzheimer's disease patient indicator) under `Disease state (cross-population indicators)` while extracting `PerezRuixo_2025_posdinemab.R`. Source acts on baseline free p217+tau (R0) in CSF only; PK parameters were unaffected by AD status. Reference category is the pooled healthy-volunteer cohort.
- **2026-04-27** — Renamed `COMBO_LD` → `COMBO_LEN_DEX` for clarity (lenalidomide+dexamethasone spelled out to avoid ambiguity with other Ld-like abbreviations). `COMBO_LD` recorded as a retired source alias. All references in `Ide_2020_elotuzumab.R` and its validation vignette updated.
- Subsequent additions: append new canonical entries as new papers are processed. When adding, bump the audit-completed count in the summary below.

## Summary

- Files audited: 69 R files under `inst/modeldb/` (20 of which reference covariates).
- Canonical H3 entries: 58 (61 parsed entries when counting each `ooc<n>` individually; was ~57 before the 2026-04-20 mergers of `hsCRP`+`BLCRP`+standard-`CRP` → `CRP`, `eGFR`+`CRCL_BSA` → `CRCL`, and `ADA_TITRE`+`ADA_TITER` → `ADA_TITER`; +1 `PHASE2` and +1 `WBC` on 2026-04-21 from Farrell 2012 / Mould 2007; +3 on 2026-04-21 for `STEROID`, `BGENE21`, `COHDOSE` from Narwal 2013).
- Scope: general: 35. Scope: specific: 26 (counting each `ooc<n>` individually, or 23 counting the `### ooc1, ooc2, ooc3, ooc4` heading as one entry).
- Aliases mapped (selected): SEXM→SEXF, ADA→ADA_POS, ADA_TITRE/ADAT→ADA_TITER, BLACK→RACE_BLACK, ASIAN→RACE_ASIAN, MULTIRACIAL→RACE_MULTI, BLACK_OTH→RACE_BLACK_OTH, ASIAN_AMIND_MULTI→RACE_ASIAN_AMIND_MULTI, DVID→STUDY1/STUDY5, CRE→CREAT, hsCRP/HSCRP/CRPHS/BLCRP→CRP, eGFR/EGFR/CRCL_BSA/1.73*CrCl/BSA→CRCL, DP2→FORM_DP2, DISEXT→DISEXT_EP/DISEXT_OTHER, BEASI→EASI, BEOS→EOS, GAST→PRIOR_GAST, COMB→COMB_EOX, TUMTP→TUMTP_CHL/TUMTP_GC, BSTEROID→STEROID, DOSE→COHDOSE.
- Canonical H3 entries: 53 (56 parsed entries when counting each `ooc<n>` individually; was ~57 before the 2026-04-20 mergers of `hsCRP`+`BLCRP`+standard-`CRP` → `CRP`, `eGFR`+`CRCL_BSA` → `CRCL`, and `ADA_TITRE`+`ADA_TITER` → `ADA_TITER`).
- Scope: general: 33. Scope: specific: 23 (counting each `ooc<n>` individually, or 20 counting the `### ooc1, ooc2, ooc3, ooc4` heading as one entry).
- Aliases mapped (selected): SEXM→SEXF, ADA→ADA_POS, ADA_TITRE/ADAT→ADA_TITER, BLACK→RACE_BLACK, ASIAN→RACE_ASIAN, MULTIRACIAL→RACE_MULTI, BLACK_OTH→RACE_BLACK_OTH, ASIAN_AMIND_MULTI→RACE_ASIAN_AMIND_MULTI, DVID→STUDY1/STUDY5, CRE→CREAT, hsCRP/HSCRP/CRPHS/BLCRP→CRP, eGFR/EGFR/CRCL_BSA/1.73*CrCl/BSA→CRCL, DP2→FORM_DP2, DISEXT→DISEXT_EP/DISEXT_OTHER, BEASI→EASI, BEOS→EOS, GAST→PRIOR_GAST, COMB→COMB_EOX, TUMTP→TUMTP_CHL/TUMTP_GC, DX→IBD_CD, AZA→CONMED_AZA, MP→CONMED_MP, MTX→CONMED_MTX, AMINO→CONMED_AMINO.
