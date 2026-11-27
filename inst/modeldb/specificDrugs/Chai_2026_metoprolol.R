Chai_2026_metoprolol <- function() {
  description <- "One-compartment population PK model for oral immediate-release metoprolol tartrate with first-order absorption (rate fixed from Eugene 2016) and first-order elimination in older Chinese adults with cardiovascular-kidney-metabolic (CKM) syndrome stage 3-4; CYP2D6 rs1065852 T/T genotype and a CKM2S2-BAG score above 11 each lower apparent clearance (Chai 2026)."
  reference <- paste(
    "Chai H, Sia JEV, Hu D, Jia Y, Hu S, Wu X, Liu D, Lai X, Cui C. (2026).",
    "Metoprolol Population Pharmacokinetics in Older Chinese Patients With",
    "CKM Syndrome: Joint Effects of rs1065852 and CKM2S2-BAG Score on Clearance.",
    "Clin Transl Sci 19(8):e70677.",
    "doi:10.1111/cts.70677.",
    sep = " "
  )
  vignette <- "Chai_2026_metoprolol"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "metoprolol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "metoprolol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CYP2D6_STAR10_HOM = list(
      description = "CYP2D6 rs1065852 (c.100C>T, *10 marker) homozygous T/T genotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (rs1065852 C/C or C/T, pooled)",
      notes = "Chai 2026 Results: 'CYP2D6 rs1065852 genotype (0 = C/C or C/T; 1 = T/T) was incorporated a priori as a fixed covariate on apparent clearance'. Recessive encoding: only the T/T stratum carries an effect, so CYP2D6_STAR10_HET is not used. The paper genotyped rs1065852 alone and treats it as a marker related to CYP2D6*10, not as a full star-allele or phenotype call (Methods 'Study Population and Data Collection'); a subject carrying *4 (which also carries the rs1065852 T allele) would be T-positive here. Cohort: C/C or C/T 30 (71.4%), T/T 12 (28.6%) (Table 2).",
      source_name = "rs1065852 (GENOTYPE in the Figure 3 legend)"
    ),
    SCORE_CKM2S2_BAG = list(
      description = "CKM2S2-BAG cardiovascular-risk score (Yang 2025, iScience 28:112780), integer points",
      units = "points",
      type = "continuous",
      reference_category = NULL,
      notes = "Composite score of Chai 2026 Table 1 (after Yang 2025): 1 point per 10 years of age above 40-49, 2 points for male sex, 2 points for current smoking, and 1 point per 10-unit step of SBP above 120 mmHg, FBG above 100 mg/dL, non-HDL-C above 160 mg/dL, and of eGFR below 90 mL/min/1.73 m^2. The model uses only the cohort-median split, score > 11 versus <= 11 (Chai 2026 Methods 'PopPK Modeling' and Table S3 footnote a); the dichotomy is computed inside model(). The original risk strata (<= 9 vs > 9) were screened and not retained (Table S3 row 11, dOFV -3.487). Cohort median 11, range 3-22; 18 of 42 subjects (42.9%) scored above 11 (Table 2).",
      source_name = "CKM2S2-BAG score (SCORE in the Figure 3 legend, already dichotomised)"
    )
  )

  covariatesDataExcluded <- list(
    DIS_CHF = list(
      description = "Chronic heart failure indicator",
      units = "(binary)",
      type = "binary",
      notes = "Added in the second SCM forward step (Table S3 row 32, dOFV -4.346) but removed at backward elimination (row 50, dOFV +4.346 < 6.63); not in the final model."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (eGFR)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Significant univariately on CL/F (Table S3 row 3, dOFV -5.359) but not retained once the CKM2S2-BAG split was in the model (row 20, dOFV -2.419)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 42L,
    n_studies = 1L,
    n_observations = 83L,
    age_range = "60-93 years (mean 75.5, SD 8.30); 21 aged 60-74 and 21 aged 76-93",
    weight_range = "44.0-102 kg (mean 68.3, SD 12.5)",
    sex_female_pct = 40.5,
    race_ethnicity = "Chinese (100%)",
    disease_state = "Cardiovascular-kidney-metabolic (CKM) syndrome stage 3 (33.3%), 4a (59.5%) or 4b (7.1%); 38.1% chronic heart failure, 35.7% diabetes, 19.0% CKD, 23.8% acute myocardial infarction",
    dose_range = "Immediate-release metoprolol tartrate 6.25-50 mg per administration, mostly 12.5-25 mg twice daily (one subject once daily)",
    regions = "China (Peking University Third Hospital, Beijing; NCT05893849)",
    renal_function = "eGFR median 72 (range 24-104) mL/min/1.73 m^2",
    co_medication = "CYP2D6 inhibitors and CYP3A4 inducers excluded (Appendix S1)",
    sampling = "Sparse steady-state troughs, mostly 14-16 h after dose (0.632-50.1 ng/mL); no absorption-phase samples",
    notes = "Real-world observational study; demographics from Chai 2026 Table 2 and Table S2. Covariates frequencies: rs1065852 T/T 28.6%; CKM2S2-BAG score > 11 42.9%."
  )

  ini({
    # Chai 2026 Table 3 ('Final' column) and the CL/F equation in Results
    # 'PopPK Model Analysis Identifying Significant Factors'.
    lka <- fixed(log(0.235)); label("Absorption rate constant (1/h)") # Table 3: Ka = 0.235 h-1 (fixed), taken from Eugene 2016 (ref 31) because no absorption-phase samples were taken
    lcl <- log(126); label("Apparent clearance CL/F, rs1065852 C/C or C/T and CKM2S2-BAG <= 11 (L/h)") # Table 3 Final: CL/F = 126 L/h (RSE 35%)
    lvc <- log(570); label("Apparent volume of distribution V/F (L)") # Table 3 Final: V/F = 570 L (RSE 78%)

    e_cyp2d6_star10_hom_cl <- -0.318; label("Fractional change in CL/F for rs1065852 T/T (unitless)") # Table 3 Final: Genotype on CL/F = -0.318 (RSE 50%); CL/F equation factor (1 - 0.318)^rs1065852
    e_ckm2s2_bag_gt11_cl <- -0.302; label("Fractional change in CL/F for CKM2S2-BAG score > 11 (unitless)") # Table 3 Final: CKM2S2-BAG score on CL/F = -0.302 (RSE 34%); CL/F equation factor (1 - 0.302)^CKM2S2-BAG

    # Table 3 Final: IIV CL/F 39.9% (RSE 23%, shrinkage 7%). The CL/F equation
    # prints the exponential term as e^0.159, and 0.399^2 = 0.159, so the
    # reported % is sqrt(omega^2) and omega^2 = 0.159.
    etalcl ~ 0.159

    propSd <- 0.347; label("Proportional residual error (fraction)") # Table 3 Final: RUV proportional CV 34.7% (RSE 13%)
  })

  model({
    # Cohort-median split of the CKM2S2-BAG score (Chai 2026 Methods): 1 if > 11.
    ckm2s2_bag_gt11 <- 0
    if (SCORE_CKM2S2_BAG > 11) ckm2s2_bag_gt11 <- 1

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (1 + e_cyp2d6_star10_hom_cl)^CYP2D6_STAR10_HOM *
      (1 + e_ckm2s2_bag_gt11_cl)^ckm2s2_bag_gt11
    vc <- exp(lvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg of metoprolol tartrate (the paper's dosing unit; CL/F and V/F
    # are apparent values on that basis), volume in L: mg/L x 1000 = ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
