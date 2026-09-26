Andrews_2019_tacrolimus <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for twice-daily oral immediate-release tacrolimus (Prograft) in adult kidney transplant recipients during the first 3 months post-transplantation (Andrews 2019 FINAL model). Bioavailability could not be estimated and was fixed to 1, so all disposition parameters are apparent (CL/F, V1/F, Q/F, V2/F). Seven covariate effects act on apparent oral clearance CL/F: CYP3A5 expresser status (1.631 multiplier for *1/*1 or *1/*3 carriers vs the *3/*3 non-expresser reference), CYP3A4*22 carriage (0.8 multiplier vs CYP3A4*1 or unknown), age (power exponent -0.43 centred at 55.72 years), serum albumin (power exponent 0.43 centred at 42 g/L), body surface area (power exponent 0.88 centred at 1.93 m^2), serum creatinine (power exponent -0.14 centred at 134.98 umol/L) and haematocrit (power exponent -0.76 centred at 34 percent). A single covariate acts on apparent central volume V1/F: lean body mass (power exponent 1.52 centred at 58.94 kg). Together the CL/F covariates explained 30 percent of the inter-individual variability in CL/F. Inter-individual variability is diagonal on CL/F, V1/F, V2/F and Q/F. Residual variability is a combined additive plus proportional error for immunoassay-measured samples and a proportional-only error for LC-MS/MS-measured samples, selected per sample by the IMMUNOASSAY indicator. Inter-occasion variability on CL/F (13.6 percent) reported in Table 2 is NOT encoded structurally here, per the nlmixr2lib convention for models with no operational occasion column; see the validation vignette. The companion reduced-covariate model intended for dose selection before transplantation is modellib('Andrews_2019_tacrolimus_startingdose')."
  reference <- "Andrews LM, Hesselink DA, van Schaik RHN, van Gelder T, de Fijter JW, Lloberas N, Elens L, Moes DJAR, de Winter BCM. A population pharmacokinetic model to predict the individual starting dose of tacrolimus in adult renal transplant recipients. Br J Clin Pharmacol. 2019;85(3):601-615. doi:10.1111/bcp.13838"
  vignette <- "Andrews_2019_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description = "Recipient age at transplantation",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per subject over the 3-month analysis window. Enters CL/F as the power term (AGE/55.72)^-0.43. The centring value 55.72 years is the model-building cohort median as hardcoded in the Supporting Information Data S1 NONMEM control stream ('CLAGE = ((AGE/55.72)**THETA(9))'); the published Equation (1) prints it rounded to 56 years. Younger patients have higher apparent oral clearance: Andrews 2019 Results reports that an increase in age from 25 to 65 years gives a 34 percent lower CL/F, which reproduces exactly as (65/25)^-0.43 = 0.663. Per-cohort medians (Table 1): Rotterdam 58.5 years (range 19.4-79.4), Leiden 54.0 years (range 15.0-77.0).",
      source_name = "AGE"
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying over the 3-month analysis window. Reported by Andrews 2019 in g/L, which is the canonical register scale, so no unit conversion is applied. Enters CL/F as the power term (ALB/42)^0.43; the centring value 42 g/L appears both in the published Equation (1) and in the Supporting Information Data S1 control stream ('CLALB = ((ALB/42)**THETA(10))'). The positive exponent means higher albumin gives HIGHER apparent oral clearance, which is the opposite of the protein-binding expectation and is flagged as such by the authors (Discussion): they found no corresponding effect on V1/F, and hypothesise that hypoalbuminaemia is a marker of an underlying inflammatory response that suppresses CYP3A activity. Albumin was measured only in the Rotterdam cohort (Table 1 median 42 g/L, range 12-57); it was unavailable in the Leiden and external-validation cohorts, where the authors fixed it to the population median.",
      source_name = "ALB"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying over the 3-month analysis window. Enters CL/F as the power term (BSA/1.93)^0.88; the centring value 1.93 m^2 appears both in the published Equation (1) and in the Supporting Information Data S1 control stream ('CLBSA = ((BSA/1.93)**THETA(16))'). Andrews 2019 Discussion states this is the first tacrolimus popPK model to carry BSA as a covariate, on the rationale that BSA is a better indicator of metabolic mass than total body weight because it is less affected by abnormal adipose mass. Andrews 2019 Results reports that a change in BSA between 2.25 and 1.5 m^2 changes CL/F by 43 percent, which reproduces as (2.25/1.5)^0.88 = 1.429. Per-cohort medians (Table 1): Rotterdam 2.03 m^2 (range 1.24-2.66), Leiden 1.90 m^2 (range 1.33-2.48). The paper does not name the BSA formula used.",
      source_name = "BSA"
    ),
    CREAT = list(
      description = "Serum creatinine concentration",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying over the 3-month analysis window; changes substantially in the days after transplantation as graft function establishes. Reported by Andrews 2019 in umol/L. Enters CL/F as the power term (CREAT/134.98)^-0.14. The centring value 134.98 umol/L is the model-building cohort median as hardcoded in the Supporting Information Data S1 control stream ('CLCREAT = ((CREAT/134.98)**THETA(11))'); the published Equation (1) prints it rounded to 135. The negative exponent means lower creatinine gives higher apparent oral clearance. Tacrolimus undergoes almost no renal elimination, so the authors state in the Discussion that the mechanism is unclear, and cite intrarenal CYP3A5-dependent tacrolimus metabolism and indirect effects of uraemic toxins on hepatic metabolism as candidate explanations. Per-cohort medians (Table 1): Rotterdam 137 umol/L (range 38-1885), Leiden 124 umol/L (range 62-920). To use a dataset recording creatinine in mg/dL, multiply by 88.4 before passing it to this model.",
      source_name = "CREAT"
    ),
    HCT = list(
      description = "Haematocrit, packed red blood cell volume fraction",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying over the 3-month analysis window. IMPORTANT SCALE NOTE: Andrews 2019 reports haematocrit as a volume FRACTION in L/L (cohort median 0.34, range 0.15-0.80 in Table 1) and centres Equation (1) at 0.34 L/L, but the canonical register column HCT is on the PERCENT scale. Both the per-subject column values and the centring value are therefore multiplied by 100 for this implementation: the reference is recorded as 34 percent and the column must be supplied in percent. This follows the Zhou_2025_tacrolimus.R precedent. Because the effect enters as a ratio (HCT/34)^-0.76, mixing the percent and fraction scales between column and reference would misscale CL/F by 100^-0.76, a factor of about 30. Enters CL/F as the power term (HCT/34)^-0.76. The negative exponent reflects that roughly 70-80 percent of whole-blood tacrolimus is bound to erythrocytes (Andrews 2019 Discussion), so a higher haematocrit sequesters more of the measured whole-blood drug and lowers apparent whole-blood clearance; this is a consequence of the partitioning and of the whole-blood assay matrix, not a change in intrinsic clearance.",
      source_name = "HCT"
    ),
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reported by Andrews 2019 under the label 'lean bodyweight (LBW)', a registered synonym of the canonical LBM column. Time-varying over the 3-month analysis window. This is the only covariate acting on a volume: it enters V1/F as the power term (LBM/58.94)^1.52. Neither the exponent's centring value nor the V1/F covariate equation itself appears in the main paper -- the published Equation (1) covers CL/F only, and Table 2 reports the exponent 1.52 (RSE 20 percent) with no reference value. The centring value 58.94 kg is recovered from the Supporting Information Data S1 NONMEM control stream, which hardcodes it as 'V2LBW = ((LBW/58.94)**THETA(15))' (NONMEM ADVAN4 TRANS4 names the central compartment V2, corresponding to the paper's V1/F). Andrews 2019 does not name the body-composition formula used to compute LBW, so a downstream user should check that whichever formula they apply (James, Boer or Hume) reproduces the source distribution: per-cohort medians (Table 1) are Rotterdam 64.0 kg (range 33.1-85.3) and Leiden 55.9 kg (range 33.6-81.7), and the control-stream centring value 58.94 kg sits between them, consistent with the combined-cohort median.",
      source_name = "LBW"
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser indicator: 1 if the patient carries at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3), 0 if a non-expresser (*3/*3 or *3/*6).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5 non-expresser: *3/*3 or *3/*6)",
      notes = "Time-fixed germline genotype (CYP3A5*3, rs776746). Andrews 2019 Results states that graphical analysis showed no difference in the effect on CL/F between the two expresser genotypes, so *1/*1 and *1/*3 are pooled into a single expresser stratum; the published Equation (1) writes this pooling explicitly as '1.631 if CYP3A5*1/*3 or CYP3A5*1/*1'. The Supporting Information Data S1 control stream parameterises it as CLCYP3A5 = 1 + THETA(13), i.e. the non-expresser reference multiplier is structurally exactly 1.0 and only the expresser increment is estimated. Genotype distribution (Table 1), Rotterdam / Leiden: *1/*1 9 (3.8%) / 4 (4.0%); *1/*3 56 (23.6%) / 17 (17.0%); *3/*3 172 (72.6%) / 76 (76.0%); *3/*6 0 (0%) / 3 (3.0%). The *3/*6 genotype is a non-expresser and is pooled with *3/*3. There was no deviation from Hardy-Weinberg equilibrium.",
      source_name = "CYP3A5"
    ),
    SNP_CYP3A4_RS35599367 = list(
      description = "CYP3A4*22 carrier indicator (rs35599367): 1 if the patient carries at least one *22 allele, 0 if CYP3A4*1 homozygous wild-type or if the genotype was not determined.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A4*1 wild-type, or genotype unknown -- Andrews 2019 Equation (1) pools these two groups)",
      notes = "Time-fixed germline genotype. The published Equation (1) pools unknown-genotype subjects with the CYP3A4*1 wild-type reference, writing the reference stratum as '1.0 if CYP3A4*1 or unknown'; this implementation follows that grouping, so CYP3A4_RS35599367 = 0 means wild-type OR unknown. The Supporting Information Data S1 control stream confirms the pooling with an explicit missing-data branch ('IF(CYP3A4.EQ.-99) CYP3A4 = 1 ; Missing data') and parameterises the effect as CLCYP3A4 = 1 + THETA(14), so the wild-type reference multiplier is structurally exactly 1.0. Genotype distribution (Table 1), Rotterdam / Leiden: *1 205 (86.5%) / 91 (91.0%); *22 22 (9.3%) / 9 (9.0%); Unknown 10 (4.2%) / 0 (0%). Andrews 2019 Discussion notes that carriers required about 20 percent less tacrolimus independent of CYP3A5 genotype status.",
      source_name = "CYP3A4"
    ),
    IMMUNOASSAY = list(
      description = "Per-sample bioanalytical assay indicator: 1 if the tacrolimus whole-blood concentration was measured by immunoassay (antibody-conjugated magnetic immunoassay ACMIA or enzyme multiplied immunoassay technique EMIT); 0 if measured by LC-MS/MS.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (LC-MS/MS reference method)",
      notes = "Per-sample (per-row) indicator, not a subject-level covariate. Andrews 2019 Methods reports that the Rotterdam cohort was assayed by two immunoassays (ACMIA, LLOQ 1.5 ng/mL; EMIT, LLOQ 2.0 ng/mL; shared ULOQ 30.0 ng/mL) while the Leiden cohort was assayed by validated LC-MS/MS (LLOQ 1.0 ng/mL, ULOQ 50.0 ng/mL). The authors state that building the different analytical techniques into the residual error model improved the base model, and the Supporting Information Data S1 control stream implements it as a DVID branch: DVID 1 (immunoassay) gets a combined proportional plus additive error, DVID 2 (LC-MS/MS) gets a proportional-only error. Andrews 2019 pooled the two immunoassay platforms (ACMIA and EMIT) under a single residual-error magnitude rather than separating them; a dataset that needs the ACMIA-vs-EMIT contrast instead would use the ASSAY_CMIA-style two-immunoassay canonical. In a pure-LC-MS/MS prospective dataset, set IMMUNOASSAY = 0 for every row, at which point the immunoassay residual parameters become non-identifiable.",
      source_name = "FLAG / DVID"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 337L,
    n_studies = 2L,
    age_range = "15.0-79.4 years (Rotterdam 19.4-79.4; Leiden 15.0-77.0)",
    age_median = "58.5 years (Rotterdam) and 54.0 years (Leiden); combined-cohort median 55.72 years per the Data S1 control-stream centring value",
    weight_range = "37.6-132.0 kg (Rotterdam) and 40.0-114.0 kg (Leiden)",
    weight_median = "79.4 kg (Rotterdam) and 74.0 kg (Leiden)",
    sex_female_pct = 39.5,
    race_ethnicity = c(Caucasian = 78.3, Asian = 9.2, `African descent` = 7.1, Other = 5.3),
    disease_state = "Adult renal transplant recipients during the first 3 months following kidney transplantation. All patients received oral twice-daily immediate-release tacrolimus (Prograft, Astellas) in combination with mycophenolic acid, with doses tailored by therapeutic drug monitoring.",
    dose_range = "Oral twice-daily immediate-release tacrolimus titrated by therapeutic drug monitoring. Rotterdam target pre-dose concentration (C0) 10.0-15.0 ng/mL in weeks 1-2, 8.0-12.0 ng/mL in weeks 3-4, and 5.0-10.0 ng/mL after week 4. Leiden targeted an AUC of 210 ng*h/mL (corresponding C0 10.0-15.0 ng/mL) for the first 6 weeks, then an AUC of 125 ng*h/mL (corresponding C0 4.0-9.0 ng/mL). The simulation trial in the paper compares a standard 0.2 mg/kg/day bodyweight-based dose against the starting-dose model.",
    regions = "The Netherlands (Erasmus MC, Rotterdam -- model building group 1, n = 237; Leiden University Medical Center -- model building group 2, n = 100).",
    n_observations = 4527L,
    sampling_window = "Tacrolimus whole-blood concentrations collected from 24 h before transplantation until 3 months post-transplantation; 4527 samples over 337 patients (concentration range 1.6-96.0 ng/mL). The Rotterdam cohort contributed 3661 samples, predominantly pre-dose concentrations (C0) starting on day 3 post-transplant; the Leiden cohort contributed 866 samples including full pharmacokinetic curves (pre-dose and 1, 2, 3, 4, 5, 6 h post-ingestion) drawn at steady state a median of 2 weeks after transplantation. Roughly 81 percent of the model-building observations are pre-dose concentrations, which the authors name as a study limitation. Three samples (0.07%) fell below the immunoassay LLOQ and were discarded; 40 samples (0.88%) exceeded the ULOQ and were estimated by NONMEM.",
    assay = "Whole-blood tacrolimus. Rotterdam cohort: antibody-conjugated magnetic immunoassay (ACMIA, LLOQ 1.5 ng/mL) and enzyme multiplied immunoassay technique (EMIT, LLOQ 2.0 ng/mL), shared ULOQ 30.0 ng/mL. Leiden cohort and external validation cohort: validated LC-MS/MS (LLOQ 1.0 ng/mL, ULOQ 50.0 ng/mL). The residual-error model carries separate magnitudes for the pooled immunoassay arm and the LC-MS/MS arm, switched per sample by the IMMUNOASSAY indicator.",
    iov_structure = "Inter-occasion variability was estimated on CL/F (13.6 percent CV, Table 2 final model), with an occasion defined as the measurement of a pre-dose concentration (C0); the Data S1 control stream implements it as three SAME-blocked OMEGA elements with a comment noting it was extended to 41 occasions in the production run. This model file does NOT encode IOV structurally, following the nlmixr2lib convention for source models with no operational occasion column (Andrews_2017_tacrolimus.R and Bukkems_2021_raltegravir.R precedents). A downstream user who wants to simulate IOV can add an occasion indicator and a per-occasion eta on CL/F with variance log(1 + 0.136^2) = 0.018386; see the validation vignette.",
    external_validation = "An independent cohort of 304 adult renal transplant recipients, entirely Caucasian and entirely LC-MS/MS-assayed (1334 samples), none of whom were in the model-building cohort. 297 contributed pre-dose concentrations only; 7 participated in the Symphony-Elite study and contributed full curves (pre-dose and 0.3, 0.7, 1.3, 2, 3, 4, 6, 8, 10, 12 h post-ingestion). Albumin was not measured in this cohort, so the authors fixed albumin to the population median for the external validation and note that the albumin effect could not be externally validated.",
    variability_explained = "The seven CL/F covariates together explained 30 percent of the inter-individual variability in CL/F (Andrews 2019 Results and Conclusions).",
    notes = "Pooled analysis of two Dutch cohorts: 237 patients from a randomized controlled trial at Erasmus MC Rotterdam (with additional PK data retrieved retrospectively from medical records) and 100 patients receiving routine clinical care at LUMC Leiden. Baseline demographics are in Table 1, which reports each cohort separately and does not publish a combined-cohort row; the combined-cohort covariate medians used as the model's centring values are recovered from the Supporting Information Data S1 control stream. Albumin, ASAT, bilirubin, total protein, CRP, ABCB1 and POR*28 genotype were measured only in the Rotterdam cohort. The model is intended for IMMEDIATE-RELEASE twice-daily oral tacrolimus in ADULTS during the FIRST 3 MONTHS post-transplant; it is not validated for once-daily extended-release formulations, for children, or for the stable maintenance phase."
  )

  ini({
    # Final-model fixed-effect estimates from Andrews 2019 Table 2, column
    # 'Final model (RSE %) [shrinkage]'. Bioavailability could not be
    # estimated and was fixed to 1 (Methods, Base model development), so every
    # clearance and volume below is an APPARENT value (CL/F, V1/F, Q/F, V2/F).
    # The reference subject for the typical values is a CYP3A5 non-expresser
    # (*3/*3), CYP3A4*1-or-unknown patient at the combined-cohort covariate
    # medians: age 55.72 years, albumin 42 g/L, BSA 1.93 m^2, creatinine
    # 134.98 umol/L, haematocrit 34 percent, lean body mass 58.94 kg.
    ltlag <- log(0.38); label("Absorption lag time tlag (h)") # Andrews 2019 Table 2 final model tlag = 0.38 h (RSE 49%)
    lka <- log(3.58); label("Absorption rate constant ka (1/h)") # Andrews 2019 Table 2 final model ka = 3.58 (RSE 40%); the Table 2 row label reads 'l h-1' which is a typo -- a first-order absorption rate constant has units 1/h
    lcl <- log(23.0); label("Apparent oral clearance CL/F at the reference subject (L/h)") # Andrews 2019 Table 2 final model CL/F = 23.0 L/h (RSE 3%); also the leading coefficient of Equation (1)
    lvc <- log(692); label("Apparent central volume V1/F at the reference subject (L)") # Andrews 2019 Table 2 final model V1/F = 692 L (RSE 8%)
    lq <- log(11.6); label("Apparent inter-compartmental clearance Q/F (L/h)") # Andrews 2019 Table 2 final model Q/F = 11.6 L/h (RSE 10%)
    lvp <- log(5340); label("Apparent peripheral volume V2/F (L)") # Andrews 2019 Table 2 final model V2/F = 5340 L (RSE 22%)

    # Covariate effects on CL/F. Andrews 2019 Equation (1):
    #   CL/F = 23 * [(1.0 if CYP3A5*3/*3) or (1.631 if CYP3A5*1/*3 or *1/*1)]
    #             * [(1.0 if CYP3A4*1 or unknown) or (0.8 if CYP3A4*22)]
    #             * (Age/56)^-0.43 * (Albumin/42)^0.43 * (BSA/1.93)^0.88
    #             * (Creatinine/135)^-0.14 * (Hematocrit/0.34)^-0.76
    # The five exponent SIGNS are stated in Table 2 and independently confirmed
    # by the Results sentence 'Higher body surface area (BSA), lower serum
    # creatinine, younger age, higher albumin and lower haematocrit levels were
    # identified as covariates enhancing tacrolimus clearance'.
    # Centring values: Equation (1) prints age and creatinine rounded (56, 135);
    # this file uses the unrounded control-stream constants 55.72 and 134.98
    # (Data S1), which are the values the model was actually fitted with. The
    # haematocrit centring value is converted from the paper's 0.34 L/L to the
    # canonical percent scale (34%).
    e_cyp3a5_expr_cl <- 1.631; label("CYP3A5 expresser (*1/*1 or *1/*3) multiplier on CL/F") # Andrews 2019 Equation (1) = 1.631; Table 2 final model CYP3A5*1 = 1.63 (RSE 15%)
    e_cyp3a4_22_cl <- 0.8; label("CYP3A4*22 carrier multiplier on CL/F") # Andrews 2019 Equation (1) = 0.8; Table 2 final model CYP3A4*22 = 0.80 (RSE 32%)
    e_age_cl <- -0.43; label("Age power exponent on CL/F, centred at 55.72 years (unitless)") # Andrews 2019 Table 2 final model Age = -0.43 (RSE 19%); sign confirmed by the Results statement that younger age enhances clearance
    e_alb_cl <- 0.43; label("Serum albumin power exponent on CL/F, centred at 42 g/L (unitless)") # Andrews 2019 Table 2 final model Albumin = 0.43 (RSE 30%); sign confirmed by the Results statement that higher albumin enhances clearance
    e_bsa_cl <- 0.88; label("Body surface area power exponent on CL/F, centred at 1.93 m^2 (unitless)") # Andrews 2019 Table 2 final model BSA = 0.88 (RSE 24%); sign confirmed by the Results statement that higher BSA enhances clearance
    e_creat_cl <- -0.14; label("Serum creatinine power exponent on CL/F, centred at 134.98 umol/L (unitless)") # Andrews 2019 Table 2 final model Creatinine = -0.14 (RSE 26%); sign confirmed by the Results statement that lower creatinine enhances clearance
    e_hct_cl <- -0.76; label("Haematocrit power exponent on CL/F, centred at 34 percent (unitless)") # Andrews 2019 Table 2 final model Haematocrit = -0.76 (RSE 11%); sign confirmed by the Results statement that lower haematocrit enhances clearance

    # Covariate effect on V1/F. Neither this equation nor its centring value is
    # printed in the main paper: Equation (1) covers CL/F only and Table 2
    # reports just the exponent. Both are recovered from the Supporting
    # Information Data S1 NONMEM control stream, which hardcodes
    # 'V2LBW = ((LBW/58.94)**THETA(15))' -- NONMEM ADVAN4 TRANS4 calls the
    # central compartment V2, which is the paper's V1/F.
    e_lbm_vc <- 1.52; label("Lean body mass power exponent on V1/F, centred at 58.94 kg (unitless)") # Andrews 2019 Table 2 final model Lean bodyweight = 1.52 (RSE 20%); centring value 58.94 kg from Supporting Information Data S1 control stream

    # Diagonal inter-individual variability. Andrews 2019 Table 2 reports IIV
    # as %CV and documents no inter-eta correlations (the Data S1 control
    # stream uses a diagonal $OMEGA for the four disposition etas). Variances
    # on the internal log scale are omega^2 = log(1 + CV^2):
    #   CL/F CV 38.6% -> log(1 + 0.386^2) = 0.138889
    #   V1/F CV 49.2% -> log(1 + 0.492^2) = 0.216775
    #   V2/F CV 53.0% -> log(1 + 0.530^2) = 0.247563
    #   Q/F  CV 78.7% -> log(1 + 0.787^2) = 0.482037
    etalcl ~ 0.138889 # Andrews 2019 Table 2 final model IIV CL/F 38.6% CV [8% shrinkage]
    etalvc ~ 0.216775 # Andrews 2019 Table 2 final model IIV V1/F 49.2% CV [25% shrinkage]
    etalvp ~ 0.247563 # Andrews 2019 Table 2 final model IIV V2/F 53.0% CV [39% shrinkage]
    etalq ~ 0.482037 # Andrews 2019 Table 2 final model IIV Q/F 78.7% CV [28% shrinkage]

    # Residual error, switched per sample by the IMMUNOASSAY indicator. The
    # Data S1 control stream branches on DVID: DVID 1 (immunoassay) gets
    # Y = F + F*EPS*THETA(2) + EPS*THETA(17), a combined proportional plus
    # additive error; DVID 2 (LC-MS/MS) gets Y = F + F*EPS*THETA(1), a
    # proportional-only error. $SIGMA is 1 FIX throughout, so the THETAs are
    # the residual standard deviations directly. The LC-MS/MS additive SD is
    # therefore structurally absent, encoded here as fixed(0).
    propSd_immuno <- 0.177; label("Proportional residual SD for immunoassay samples (fraction)") # Andrews 2019 Table 2 final model Proportional Immunoassay = 17.7% (RSE 7%) [22% shrinkage]
    addSd_immuno <- 0.88; label("Additive residual SD for immunoassay samples (ng/mL)") # Andrews 2019 Table 2 final model Additive Immunoassay = 0.88 ug/L (RSE 13%) [22% shrinkage]; 1 ug/L = 1 ng/mL
    propSd_lcms <- 0.245; label("Proportional residual SD for LC-MS/MS samples (fraction)") # Andrews 2019 Table 2 final model Proportional LC-MS/MS = 24.5% (RSE 5%) [12% shrinkage]
    addSd_lcms <- fixed(0); label("Additive residual SD for LC-MS/MS samples (ng/mL)") # Structurally absent: Data S1 control stream gives DVID 2 (LC-MS/MS) a proportional-only error, and Table 2 reports an additive term only for the immunoassay arm
  })

  model({
    # CYP3A5 expresser multiplier on CL/F. The Data S1 control stream
    # parameterises this as CLCYP3A5 = 1 + THETA(13), so the non-expresser
    # reference multiplier is structurally exactly 1.0 and is not estimated.
    # *1/*1 and *1/*3 are pooled into the single expresser stratum, and *3/*6
    # is pooled with *3/*3 as a non-expresser.
    f_cyp3a5 <- 1 + (e_cyp3a5_expr_cl - 1) * CYP3A5_EXPR

    # CYP3A4*22 carrier multiplier on CL/F, likewise parameterised as
    # CLCYP3A4 = 1 + THETA(14) in the control stream. Unknown-genotype
    # subjects are pooled with the CYP3A4*1 wild-type reference, per
    # Equation (1) and the control stream's explicit missing-data branch.
    f_cyp3a4 <- 1 + (e_cyp3a4_22_cl - 1) * SNP_CYP3A4_RS35599367

    # Continuous covariate power terms on CL/F. Centring values are the
    # combined-cohort medians hardcoded in the Data S1 control stream; HCT is
    # on the canonical percent scale, so its reference is 34 (the paper's
    # 0.34 L/L times 100).
    f_age <- (AGE / 55.72)^e_age_cl
    f_alb <- (ALB / 42)^e_alb_cl
    f_bsa <- (BSA / 1.93)^e_bsa_cl
    f_creat <- (CREAT / 134.98)^e_creat_cl
    f_hct <- (HCT / 34)^e_hct_cl

    # The single covariate on a volume: lean body mass on apparent central
    # volume, centred at the control stream's 58.94 kg.
    f_lbm <- (LBM / 58.94)^e_lbm_vc

    # Individual PK parameters. Andrews 2019 Equation (1) multiplies the
    # typical CL/F by all seven covariate factors.
    tlag <- exp(ltlag)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * f_cyp3a5 * f_cyp3a4 * f_age * f_alb * f_bsa * f_creat * f_hct
    vc <- exp(lvc + etalvc) * f_lbm
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # Two-compartment oral disposition (NONMEM ADVAN4 TRANS4). The dose lands
    # in `depot`; bioavailability was fixed to 1, and is implicit in the
    # apparent CL/F and V/F parameterisation.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Tacrolimus whole-blood concentrations are reported in ng/mL. The dose is
    # in mg and vc is in L, so central/vc is mg/L = ug/mL; multiply by 1000 to
    # get ng/mL. This reproduces the Data S1 control stream's scaling
    # 'S2 = V2 / 1000'.
    Cc <- central / vc * 1000

    # Per-sample assay-conditional residual error (IMMUNOASSAY: 1 =
    # immunoassay, 0 = LC-MS/MS reference), following the
    # Andrews_2017_tacrolimus.R precedent for the same analytical-method
    # switch. addSd_lcms is fixed at 0, so LC-MS/MS rows reduce to the
    # proportional-only error the paper specifies.
    addSd <- addSd_immuno * IMMUNOASSAY + addSd_lcms * (1 - IMMUNOASSAY)
    propSd <- propSd_immuno * IMMUNOASSAY + propSd_lcms * (1 - IMMUNOASSAY)
    Cc ~ add(addSd) + prop(propSd)
  })
}
