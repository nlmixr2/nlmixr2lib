Li_2023_amisulpride <- function() {
  description <- "Population PK model for oral amisulpride in Chinese adults with schizophrenia (Li 2023). One-compartment disposition with first-order absorption and linear elimination, parameterised on the apparent (oral) scale as CL/F and V/F. Renal function, as the race-free CKD-EPI estimated glomerular filtration rate (mL/min/1.73 m^2), was included a priori in the structural model and enters CL/F as a power function centred on 113.84 mL/min/1.73 m^2; no other covariate was retained. Exponential inter-individual variability on CL/F and V/F (none on Ka) and a combined additive plus proportional residual error. Developed on 776 therapeutic-drug-monitoring patients and externally validated on a further 145."
  reference <- "Li A, Mak WY, Ruan T, Dong F, Zheng N, Gu M, Guo W, Zhang J, Cheng H, Ruan C, Shi Y, Zang Y, Zhu X, He Q, Xiang X, Wang G, Zhu X. Population pharmacokinetics of Amisulpride in Chinese patients with schizophrenia with external validation: the impact of renal function. Front Pharmacol. 2023;14:1215065. doi:10.3389/fphar.2023.1215065."
  vignette <- "Li_2023_amisulpride"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "amisulpride", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "amisulpride", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Estimated glomerular filtration rate from the 2021 creatinine CKD-EPI equation WITHOUT the race component, BSA-normalised (Li 2023 Table 1 footnote: 'eGFR calculated using the CKD-EPI equation without the race component'; Supplementary Table S1 CKD-EPI row '142 x (Scr/A)^B x 0.9938^AGE x 1.012 (if female)').",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline value per subject; the source could not model time-varying renal function (Li 2023",
        "Methods 2.1). Included a priori in the structural model (dOFV = -103.7, Results 3.2); CKD-EPI was",
        "selected over Cockcroft-Gault, lean-body-weight Cockcroft-Gault, MDRD and Chinese-adjusted MDRD by",
        "AIC (Supplementary Table S2, model 9 'Ami final', AIC 27190.8). Applied as",
        "CL/F = 60.5 * (eGFR / 113.84)^0.817 (Li 2023 Table 2 equation; NONMEM code in Supplementary",
        "Figure S6 'TVCL = THETA(1)*(GFR_CKD/113.84)**THETA(6)'). Development-cohort median 114",
        "(IQR 101-130) mL/min/1.73 m^2; 89.0 percent CKD stage G1 and only 5 subjects in G3a (Table 1), so",
        "predictions below about 60 mL/min/1.73 m^2 are extrapolations.",
        sep = " "
      ),
      source_name = "eGFR (GFR_CKD in the NONMEM code)"
    )
  )

  covariatesDataExcluded <- list(
    SMOKE = list(
      description = "Current-smoker indicator",
      units = "(binary)",
      type = "binary",
      notes = "The only covariate that exploratory analysis flagged for further testing (Li 2023 Results 3.2; Supplementary Figure S5); not retained after stepwise covariate modelling (forward p < 0.01, backward p < 0.001). 699 of 776 development-cohort patients (90.1 percent) were recorded as smokers (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 776L,
    n_studies = 1L,
    age_range = "Median 33 years (IQR 23-47) in the development cohort (Li 2023 Table 1)",
    age_median = "33 years",
    weight_range = "Median 66 kg (IQR 57-79) in the development cohort (Li 2023 Table 1)",
    weight_median = "66 kg",
    sex_female_pct = 61.2,
    race_ethnicity = "Han Chinese 96.4 percent, other 3.6 percent (development cohort, Li 2023 Table 1)",
    disease_state = "Adults with schizophrenia receiving oral amisulpride under routine therapeutic drug monitoring at Beijing Anding Hospital.",
    dose_range = "Oral amisulpride, daily dose 50-1220 mg across both cohorts (development cohort mean 566 +/- 264, median 575 mg/day; Li 2023 Table 1 and Results 3.1).",
    regions = "China (Beijing), December 2018 to March 2022",
    renal_function = "CKD-EPI eGFR (no race component) median 114 (IQR 101-130) mL/min/1.73 m^2 in the development cohort; 25.4 percent above 130 mL/min/1.73 m^2 (hyperfiltration); CKD stage G1 89.0 percent, G2 10.3 percent, G3a 0.6 percent (Li 2023 Table 1).",
    co_medication = "Lithium co-administered in 96 of 776 development-cohort patients (12.4 percent; Li 2023 Table 1); not retained as a covariate.",
    notes = paste(
      "Retrospective TDM analysis; about 3 random (non-protocol) samples per patient (median 3, IQR 2-4).",
      "Total plasma amisulpride by UPLC-MS/MS, LLOQ 15.625 ng/mL (Li 2023 Methods 2.2).",
      "Development cohort n = 776 (admitted after 15 November 2019); external validation cohort n = 145",
      "(admitted before that date), evaluated with MAXEVAL = 0 (RMSE 53.7 percent, MRPE 10.3 percent).",
      "NONMEM 7.5, FOCE-I.",
      sep = " "
    )
  )

  ini({
    # Li 2023 Table 2 (final estimates). The NONMEM code (Supplementary
    # Figure S6) is ADVAN2 TRANS2 with the dose in mg and S2 = V/1000, i.e.
    # concentrations in ng/mL. Oral data only, so F is absorbed into the
    # apparent parameters.
    lka <- log(0.106); label("First-order absorption rate constant Ka (1/h)")  # Li 2023 Table 2: ka = 0.106 1/h (RSE 12%); bootstrap 0.106 [0.0916-0.1210]
    lcl <- log(60.5); label("Apparent clearance CL/F at eGFR 113.84 mL/min/1.73 m^2 (L/h)")  # Li 2023 Table 2: CL/F = 60.5 L/h (RSE 2%); bootstrap 60.4 [58.62-62.37]
    lvc <- log(645); label("Apparent volume of distribution V/F (L)")  # Li 2023 Table 2: V/F = 645 L (RSE 18%); bootstrap 650 [521.3-769.6]

    # Covariate exponent; centring value 113.84 from the Table 2 equation and
    # Supplementary Figure S6 (TVCL = THETA(1)*(GFR_CKD/113.84)**THETA(6)).
    e_crcl_cl <- 0.817; label("Power exponent of CKD-EPI eGFR on CL/F (unitless)")  # Li 2023 Table 2: 'eGFR on CL/F RCL' = 0.817 (RSE 10%); bootstrap 0.821 [0.7030-0.9319]

    # IIV is exponential (Methods Eq. 1; Supplementary Figure S6
    # CL = TVCL*EXP(ETA(1)), V = THETA(2)*EXP(ETA(2))). Table 2 prints %CV only;
    # converted with omega^2 = log(1 + CV^2). The control-stream $OMEGA
    # starting values (0.187, 0.5) are not near-converged so they cannot pick
    # between that and omega^2 = CV^2 (0.1289 / 1.7135); see the vignette.
    # ETA(3) on KA is '0 FIX' in the code, i.e. no IIV on Ka, so it is omitted.
    etalcl ~ log(1 + 0.359^2)  # Li 2023 Table 2: IIV CL/F = 35.9% CV (RSE 7%, shrinkage 20%)
    etalvc ~ log(1 + 1.309^2)  # Li 2023 Table 2: IIV V/F = 130.9% CV (RSE 25%, shrinkage 46%)

    # Residual error: Supplementary Figure S6 '$ERROR W = SQRT(THETA(4)**2*IPRED**2
    # + THETA(5)**2)', Y = IPRED + W*EPS(1) with SIGMA 1 FIX, so THETA(4) and
    # THETA(5) are SDs combined on the variance scale (combined2).
    propSd <- 0.346; label("Proportional residual error (fraction)")  # Li 2023 Table 2: Prop.RE = 34.6% (RSE 3%); bootstrap 34.35% [32.69-36.46%]
    addSd <- 13.7; label("Additive residual error (ng/mL)")  # Li 2023 Table 2: Add.RE = 13.7 (RSE 36%); bootstrap 13.80 [5.048-22.397]
  })

  model({
    # Li 2023 Table 2 equation: CL/F = 60.5 * (eGFR/113.84)^0.817; V/F = 645; ka = 0.106
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (CRCL / 113.84)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and V/F in L give mg/L; x 1000 gives ng/mL (S2 = V/1000).
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
