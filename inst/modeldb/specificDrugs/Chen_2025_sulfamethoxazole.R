Chen_2025_sulfamethoxazole <- function() {
  description <- paste(
    "One-compartment population PK model for INTRAVENOUS sulfamethoxazole in",
    "adults treated with co-trimoxazole for Pneumocystis jirovecii pneumonia",
    "(Chen 2025). Both the typical volume and the typical clearance are",
    "reported per kilogram of body weight (L/kg and L/kg/h), so body weight",
    "enters both parameters linearly; this is how the 'WT on CL' step",
    "retained in the covariate search is carried, and it is why the printed",
    "model equations show no explicit weight term. Creatinine clearance",
    "scales clearance as a power function normalized to the cohort median of",
    "75.7 mL/min, and continuous renal replacement therapy multiplies",
    "clearance by exp(0.59) = 1.80 because sulfamethoxazole is removed by",
    "ultrafiltration and is not reabsorbed in the ultrafiltrate. Doses are",
    "the SULFAMETHOXAZOLE component of the combination product: a",
    "co-trimoxazole dose of X mg/kg/day delivers 5X/6 mg/kg/day of",
    "sulfamethoxazole and X/6 mg/kg/day of trimethoprim. The trimethoprim",
    "model of the same paper is a separate file; see",
    "modellib('Chen_2025_trimethoprim').",
    sep = " "
  )
  reference <- paste(
    "Chen B, Chen Y, Chen M, Mao Y, Huang Y, Zhou L, Wu W, Li X, Wu X,",
    "Cheng Y, Qiu H. Population pharmacokinetics and Monte Carlo-based",
    "dosing optimization of trimethoprim-sulfamethoxazole.",
    "Antimicrob Agents Chemother. 2025;69(11):e00519-25.",
    "doi:10.1128/aac.00519-25.",
    "Structural equations from Eq. 1 and Eq. 2 of the Results section",
    "('Population pharmacokinetic analysis'); fixed effects, interindividual",
    "variability and residual error from Table 5. The equations are rendered",
    "as images in the publisher PDF and were recovered with",
    "'pdftotext -layout'.",
    sep = " "
  )
  vignette <- "Chen_2025_cotrimoxazole"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters BOTH clearance and volume LINEARLY (exponent exactly 1),",
        "because Table 5 reports the typical values per kilogram: tv V in",
        "L/kg and tv CL in L/kg/h. Weight is therefore not visible as a",
        "separate term in the printed equations (Eq. 1 and Eq. 2), but it is",
        "the 'WT on CL' covariate that the forward/backward search retained",
        "(Table 3: adding WT on CL dropped the objective function value by",
        "7.52, and removing it from the full model raised it by 6.57,",
        "P < 0.05). Because the paper's dosing recommendations are all in",
        "mg/kg and both parameters are per-kg, simulated concentrations are",
        "independent of the weight supplied when the dose is scaled by the",
        "same weight.",
        "Cohort median 60.0 kg (Table 2)."
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = paste(
        "Creatinine clearance calculated with the Cockcroft-Gault equation.",
        "RAW mL/min, NOT normalized to 1.73 m^2 body surface area."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power term (CRCL / 75.7)^0.17 on clearance (Eq. 2). The text after",
        "the equations states explicitly that '75.7 represents the median",
        "CrCL value', matching the cohort median in Table 2 (75.7 mL/min,",
        "51.0-93.7). Methods: 'CrCL ... calculated using the Cockcroft-Gault",
        "equation', which returns raw mL/min; do NOT supply a BSA-normalized",
        "value here. The exponent is shallow (0.17) because sulfamethoxazole",
        "is largely cleared by hepatic N-acetylation and CYP2C9 hydroxylation",
        "rather than by glomerular filtration - it is the renally excreted",
        "metabolite N-acetyl sulfamethoxazole that accumulates in renal",
        "impairment (Discussion).",
        "The Monte Carlo simulations stratify on CrCL bands of <15, 15-29,",
        "30-49, 50-79 and 80-120 mL/min (Table 7).",
        sep = " "
      ),
      source_name        = "CrCL"
    ),
    RRT_CRRT_STATUS = list(
      description        = "Continuous renal replacement therapy during co-trimoxazole treatment",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CRRT)",
      notes              = paste(
        "1 = receiving CRRT, 0 = not. Enters clearance as",
        "exp(dCLdCRRT * (CRRT == 1)) (Eq. 2), i.e. a multiplicative factor of",
        "exp(0.59) = 1.80 rather than a directly reported ratio; the",
        "coefficient 0.59 in Table 5 is on the LOG scale. 19 of 79 patients",
        "(24.1%) received CRRT (Table 2). The Discussion attributes the",
        "increase to loss of tubular reabsorption plus direct removal by",
        "ultrafiltration, and cites Curkovic et al. reporting that",
        "sulfamethoxazole clearance on CRRT exceeds normal renal clearance",
        "while trimethoprim clearance does not - which is why the same",
        "covariate is absent from the trimethoprim model of this paper.",
        "The Introduction quotes an independent 3.5-fold CRRT effect",
        "(83.9 vs 24.4 mL/min) from the literature, larger than the 1.80-fold",
        "estimated here. A CRCL value must still be supplied for CRRT",
        "subjects because the power term is evaluated unconditionally.",
        sep = " "
      ),
      source_name        = "CRRT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      notes              = "Screened but not retained in the sulfamethoxazole model (Methods covariate list; absent from Table 3, Eq. 2 and Table 5). Cohort median 64 years (Table 2). Retained in the TRIMETHOPRIM search chain but likewise absent from that final model - see modellib('Chen_2025_trimethoprim')."
    ),
    HT = list(
      description        = "Height",
      units              = "cm",
      type               = "continuous",
      notes              = "Screened but not retained (Methods covariate list). Cohort median 170.0 cm (Table 2)."
    ),
    SEXF = list(
      description        = "Female sex",
      units              = "(binary)",
      type               = "binary",
      notes              = "Screened as 'gender' but not retained (Methods covariate list). 18 of 79 patients female (Table 2)."
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      notes              = "Screened but not retained (Methods covariate list). Cohort median 34.4 g/L (Table 2). The Discussion notes a literature correlation between albumin and sulfamethoxazole clearance that this cohort did not reproduce."
    ),
    TBILI = list(
      description        = "Total bilirubin",
      units              = "umol/L",
      type               = "continuous",
      notes              = "Screened but not retained (Methods covariate list). Cohort median 21.7 umol/L (Table 2)."
    ),
    DBIL = list(
      description        = "Direct bilirubin",
      units              = "umol/L",
      type               = "continuous",
      notes              = "Screened but not retained (Methods covariate list). Cohort median 10.9 umol/L (Table 2)."
    ),
    HEPIMP_SEV = list(
      description        = "Hepatic impairment severity (Child-Pugh class)",
      units              = "(categorical)",
      type               = "categorical",
      notes              = "Screened as the Child-Pugh classification (A/B/C) but not retained (Methods covariate list; Results 'Liver function ... did not exhibit statistically significant effects'). Child-Pugh A 54, B 18, C 7 (Table 2). Total protein, ALT, AST, GGT and ALP were screened alongside it and are likewise absent from the final model; NAT2 acetylator phenotype (rapid/intermediate/slow) and CYP2C9 metabolizer phenotype (normal/intermediate/poor) were genotyped and screened with no significant effect (Results, Discussion). None of these carries a reported point estimate, so none can be encoded."
    )
  )

  compartmentData <- list(
    central = list(analyte = "sulfamethoxazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 79,
    n_studies      = 1,
    age_median     = "64 years (54-73)",
    height_median  = "170.0 cm (165-175)",
    weight_median  = "60.0 kg (55-70)",
    sex_female_pct = 22.8,
    race_ethnicity = "Not reported; single-center Chinese cohort.",
    disease_state  = paste(
      "Adults (>= 18 years) with a confirmed diagnosis of Pneumocystis",
      "jirovecii pneumonia receiving intravenous co-trimoxazole, in intensive",
      "care units and general inpatient wards. Patients taking hepatic enzyme",
      "inducers (carbamazepine, rifampin) were excluded.",
      sep = " "
    ),
    renal_function = paste(
      "Cockcroft-Gault CrCL median 75.7 mL/min (51.0-93.7); serum creatinine",
      "median 77.8 umol/L (59.0-135.0). 47 of 79 patients (59.5%) had",
      "CrCL < 80 mL/min and 32 (40.5%) had CrCL >= 80 mL/min. 19 patients",
      "(24.1%) received continuous renal replacement therapy.",
      sep = " "
    ),
    hepatic_function = paste(
      "Child-Pugh A 54 (68.3%), B 18 (22.8%), C 7 (8.9%). Total bilirubin",
      "median 21.7 umol/L, albumin 34.4 g/L, ALT 34.5 U/L, AST 32.0 U/L.",
      sep = " "
    ),
    pharmacogenomics = paste(
      "NAT2 acetylator phenotype from rs1799929, rs1799930 and rs1799931:",
      "rapid 30 (37.9%), intermediate 35 (44.4%), slow 14 (17.7%). CYP2C9",
      "phenotype from rs1057910 (*3): normal 68 (86.1%), intermediate 7",
      "(8.8%), poor 4 (5.1%). Neither had a statistically significant effect.",
      sep = " "
    ),
    dose_range     = paste(
      "Routine care rather than protocol-assigned; intravenous infusions of",
      "approximately 1 h given every 6, 8 or 12 h. The Monte Carlo",
      "simulations explore co-trimoxazole 50, 55, 65, 70 and 90 mg/kg/day",
      "given twice, three times or four times daily.",
      sep = " "
    ),
    regions        = "China (Fujian Medical University Union Hospital, Fuzhou), March 2023 to October 2024.",
    n_observations = paste(
      "232 post-dose plasma concentrations from 79 patients: two to three",
      "samples per patient (an end-of-infusion peak, a pre-dose trough and/or",
      "an intermediate sample; Table 1). Sulfamethoxazole was assayed by",
      "LC-MS/MS over a calibrated range of 3.12-400.0 mg/L.",
      sep = " "
    ),
    notes          = paste(
      "Prospective single-center study; demographics from Table 2. Estimation",
      "was by first-order conditional estimation with extended least squares",
      "in Phoenix NLME 8.0, and the reported parameter precision comes from",
      "1,000-sample nonparametric bootstrapping. The efficacy/toxicity window",
      "the paper simulates against is a steady-state peak of 100-200 mg/L for",
      "sulfamethoxazole. Several arithmetic inconsistencies in the published",
      "tables are catalogued in the vignette Errata; none affects the",
      "parameters encoded here.",
      sep = " "
    )
  )

  ini({
    # Structural parameters, Table 5 ("Population pharmacokinetic model
    # estimates and bootstrap results for SMX"), "Final model / Estimate"
    # column. Both are per kilogram of body weight; model() multiplies by WT.
    lvc <- log(0.32); label("Typical weight-normalized central volume of distribution (L/kg)")                        # Table 5 "tv V (L/kg)": 0.32, RSE 3.95%, bootstrap median 0.33 (95% CI 0.30-0.35). Discussion: "the typical apparent volumes of V were 0.32 L/kg for SMX ..., consistent with prior adult pharmacokinetic studies (SMX: 0.17-0.34 L/kg)".
    lcl <- log(0.02); label("Typical weight-normalized clearance at CrCL 75.7 mL/min without CRRT (L/kg/h)")          # Table 5 "tv CL (L/kg/h)": 0.02, RSE 5.0%, bootstrap median 0.02. Discussion: "typical population CL values of 0.02 L/kg/h for SMX ..., aligning with previously reported ranges (SMX: 0.013-0.024 L/kg/h)". The Table 5 bootstrap 95% CI is printed as "0.18-0.22", which is off by a factor of ten from its own median of 0.02; see the vignette Errata.

    # Covariate effects on clearance, Eq. 2:
    #   CL = tvCL * (CrCL/75.7)^dCLdCrCL * exp[dCLdCRRT * (CRRT == 1)] * exp(etaCL)
    e_crcl_cl            <- 0.17; label("Power exponent on (CRCL/75.7) for CL (unitless)")                            # Table 5 "dCLdCrCL": 0.17, RSE 27.75%, bootstrap median 0.18 (95% CI 0.05-0.32)
    e_rrt_crrt_status_cl <- 0.59; label("Log-scale coefficient for CRRT on CL; multiplies CL by exp(0.59) = 1.80 (unitless)") # Table 5 "dCLdCRRT": 0.59, RSE 15.81%, bootstrap median 0.61 (95% CI 0.41-0.78). Enters INSIDE exp() per Eq. 2, so 0.59 is not itself the fold-change.

    # Interindividual variability. Table 5 labels these rows "omega^2" and the
    # footnote confirms "omega^2, variance of interindividual variability", so
    # the printed numbers are variances on the log scale and go into ini()
    # unchanged. Methods: "Interindividual variability was characterised using
    # an exponential residual model", matching the exp(eta) form of Eq. 1-2.
    etalvc ~ 0.08                                                                                                     # Table 5 "omega^2 V": 0.08, RSE 26.83%, bootstrap median 0.13 (95% CI 0.06-0.22); omega = 28.3% on the log scale
    etalcl ~ 0.16                                                                                                     # Table 5 "omega^2 CL": 0.16, RSE 3.91%, bootstrap median 0.17 (95% CI 0.08-0.25); omega = 40.0% on the log scale

    # Residual error. Results: "residual variability was appropriately
    # characterised using a proportional error model" for both analytes;
    # Table 3 confirms the proportional model was selected on OFV (2,091.64 vs
    # 2,091.68 additive and 2,094.38 combined).
    propSd <- 0.09; label("Proportional residual error SD (fraction)")                                                # Table 5 "Proportional error": 0.09, RSE 16.46%, bootstrap median 0.09 (95% CI 0.06-0.12)
  })

  model({
    # Cohort median creatinine clearance used as the normalizing constant.
    # Stated in the sentence following Eq. 1-4: "75.7 represents the median
    # CrCL value"; matches the Table 2 median of 75.7 mL/min.
    ref_crcl <- 75.7

    # Eq. 1: V = tvV * exp(etaV), with tvV per kg, so V (L) = tvV * WT.
    vc <- exp(lvc + etalvc) * WT

    # Eq. 2: CL = tvCL * (CrCL/75.7)^dCLdCrCL * exp[dCLdCRRT * (CRRT == 1)] *
    # exp(etaCL), with tvCL per kg, so CL (L/h) = tvCL * WT * ... The CRRT
    # indicator is used directly rather than as an equality test because
    # RRT_CRRT_STATUS is already coded 0/1.
    cl <- exp(lcl + etalcl) * WT *
      (CRCL / ref_crcl)^e_crcl_cl *
      exp(e_rrt_crrt_status_cl * RRT_CRRT_STATUS)

    kel <- cl / vc

    # One compartment with first-order elimination (Results: "The PopPK
    # characteristics of both SMX and TMP were best described by a
    # one-compartment model with first-order elimination kinetics"; Table 3
    # model 1 OFV 2,091.68 versus 2,123.61 for two compartments). Intravenous
    # only - the analysis excludes oral administration - so there is no depot
    # and no bioavailability term; dose directly into central.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
