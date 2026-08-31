Taylor_2026_methotrexate <- function() {
  description <- "Two-compartment population PK model for intravenous high-dose methotrexate (1.5-8 g/m^2 given over a 4 h infusion) in adult patients with lymphoma; clearance carries a power effect of the CKD-EPI creatinine-cystatin C 2021 eGFR equation (reference 76 mL/min/1.73 m^2, modelled longitudinally) and of baseline serum albumin (reference 4 g/dL). The cystatin C-inclusive eGFR was the paper's primary finding: it outperformed serum creatinine and the creatinine-only and cystatin C-only CKD-EPI equations for predicting methotrexate clearance (Taylor 2026)"
  reference <- paste(
    "Taylor ZL, Barreto EF, Cole KC, Rule AD, Kashani KB, Leung N,",
    "Thompson CA, Witzig TE, Ramsey LB, Barreto JN. (2026). Use of a",
    "Cystatin C-Based GFR Equation in a Population Pharmacokinetic Model of",
    "Methotrexate Clearance in Adult Patients with Lymphoma.",
    "Clin Pharmacokinet 65:465-477.",
    "doi:10.1007/s40262-026-01618-4.",
    sep = " "
  )
  vignette <- "Taylor_2026_methotrexate"
  units    <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Serum methotrexate concentrations were quantified by
  # LC-MS (lower limit of detection 10 ng/mL, ~0.02 umol/L) and reported in
  # umol/L throughout Taylor 2026 (Sect. 2.1, Tables 3 and 4), so the amount
  # unit that pairs with the published V1 (L) and CL (L/h) is umol.
  compartmentData <- list(
    central     = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate from the CKD-EPI creatinine-cystatin C 2021 equation (time-varying)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "BSA-normalized eGFR computed with the CKD-EPI creatinine-cystatin C 2021 equation from paired serum",
        "creatinine and serum cystatin C (Taylor 2026 Sect. 2.2.2). Modelled as a LONGITUDINAL (baseline and",
        "time-dependent) covariate: creatinine and cystatin C were drawn alongside each methotrexate sample",
        "(4, 12, 24, 48, 72 and 96 h after the start of the infusion) and re-converted to eGFR at each timepoint.",
        "Enters CL as a power function with reference 76 mL/min/1.73 m^2 and exponent 0.80 (Taylor 2026 Eq. 1).",
        "The 76 is double-sourced: it is the denominator printed in Eq. 1 and the ESM Table S2 model 8 row",
        "records the tested form as '(eGFR/76)^theta'.",
        "The reference is the MEDIAN of the LONGITUDINAL eGFR record, not the baseline median: Taylor 2026",
        "Sect. 3.2.2 states the reference value was set to the median GFR value, and baseline medians are higher",
        "(82.5 and 80.2 mL/min/1.73 m^2 in the two dose strata, Table 3) because eGFR falls over the post-infusion",
        "monitoring window. Taylor 2026 Table 1 heads the renal block 'Estimated kidney function (mL/min)', which",
        "lumps the BSA-normalized CKD-EPI equations in with the raw Cockcroft-Gault estimate; Sect. 2.2.2 is",
        "explicit that the CKD-EPI values are mL/min/1.73 m^2 and that is the unit used here.",
        "The cystatin C-inclusive equation was selected over CKD-EPI creatinine alone and CKD-EPI cystatin C alone",
        "because it gave the largest OFV drop (dOFV = 117.1, p < 0.0001) and explained ~11.5% of the total",
        "inter-individual variability in CL (Taylor 2026 Sect. 3.2.2).",
        sep = " "
      ),
      source_name        = "eGFR (CKD-EPI Cr-CysC)"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Taylor 2026 reports albumin in US-convention g/dL (cohort median 3.8 g/dL, IQR 3.5-4.2; Table 1) and",
        "centers Eq. 1 on 4 g/dL. The canonical register unit is SI g/L, so model() applies the inline",
        "conversion 'alb_gdL <- ALB * 0.1' before the power term; supply this column as 40 g/L to sit at the",
        "paper's reference. BASELINE value only (not time-varying), unlike CRCL. Albumin was the only variable",
        "with missingness: values for 17 of 80 patients (21%) were imputed from the population median or from",
        "regression-based estimates, with sensitivity analyses confirming the parameter estimates were stable",
        "across imputation methods (Taylor 2026 Sect. 2.2.2 and 3.1). Adding albumin on top of the eGFR model",
        "gave dOFV = 21.1 (p < 0.0001) and explained a further 2.9% of the inter-individual variability in CL.",
        "The exponent is POSITIVE (0.69): higher albumin predicts faster methotrexate clearance.",
        sep = " "
      ),
      source_name        = "Alb"
    )
  )

  # Screened during the exploratory / stepwise covariate analysis but NOT
  # retained in the final model (Taylor 2026 Sect. 2.2.2 and 3.2.2). Recorded
  # for provenance only; none of these appear in model().
  covariatesDataExcluded <- list(
    CREAT = list(
      description = "Serum creatinine (time-varying)",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Cohort baseline 0.8 +/- 0.4 mg/dL (Taylor 2026 Table 1). Statistically associated with CL in the",
        "univariate screen and carried into stepwise covariate modelling, but rejected: as a raw laboratory",
        "value it reached OFV = 13.5 versus -25.5 for serum cystatin C (ESM Table S2 models 5 and 4; the",
        "corresponding dOFV against the base model are -42.9 and -81.5), and explained only 3.7% of the",
        "inter-individual variability in CL versus 10.4% for cystatin C. It is subsumed into the retained CRCL",
        "column, which is computed from creatinine AND cystatin C together. Rejecting creatinine-alone renal",
        "markers is the paper's central claim.",
        sep = " "
      ),
      source_name = "SCr"
    ),
    CYSC = list(
      description = "Serum cystatin C (time-varying)",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Cohort baseline 1.1 +/- 0.5 mg/dL (Taylor 2026 Table 1; note the paper reports cystatin C in mg/dL,",
        "whereas mg/L is the more common clinical convention -- the numeric values, ~1 in a cohort with near-normal",
        "renal function, are those of mg/L). Outperformed serum creatinine as a raw laboratory covariate on CL",
        "(OFV = -25.5, i.e. dOFV = -81.5 against the base model, 10.4% of the inter-individual variability",
        "explained; ESM Table S2 model 4, tested as '(CysC/1.15)^theta') but was itself superseded by the",
        "combined CKD-EPI creatinine-cystatin C eGFR, which gave a much larger dOFV of 117.1.",
        sep = " "
      ),
      source_name = "CysC"
    ),
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Median 68.6 years (IQR 59.3-75.9); 64% of the cohort was 65 or older (Taylor 2026 Table 1).",
        "Statistically associated with CL in the univariate screen and tested as '(Age/65)^theta'",
        "(dOFV = -21.7; ESM Table S2 model 2), but once the CKD-EPI Cr-CysC eGFR was in the model age added",
        "only dOFV = -3.0 (ESM Table S2 model 9) and was not retained in the final model.",
        sep = " "
      ),
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight at baseline",
      units       = "kg",
      type        = "continuous",
      notes       = "Median 80.5 kg (IQR 69.7-92.3; Taylor 2026 Table 1). Screened as a demographic covariate and NOT advanced to stepwise covariate modelling. Note that this model therefore carries no body-size term at all -- neither allometric nor BSA-normalized -- on any structural parameter.",
      source_name = "WT"
    ),
    BSA = list(
      description = "Body surface area at baseline (Du Bois method)",
      units       = "m^2",
      type        = "continuous",
      notes       = "Median 1.97 m^2 (IQR 1.80-2.14; 2.04 in men, 1.77 in women; Taylor 2026 Table 1 and Sect. 3.1). Used only to convert the protocol dose (g/m^2) into a delivered dose (g); screened as a covariate on CL and V1 and NOT advanced to stepwise covariate modelling.",
      source_name = "BSA"
    ),
    BMI = list(
      description = "Body mass index at baseline",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened both as a continuous variable and as a categorical variable using the CDC adult underweight / healthy / overweight / obese bands (Taylor 2026 Sect. 2.2.2); not advanced to stepwise covariate modelling.",
      source_name = "BMI"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "categorical",
      notes       = "54 of 80 (68%) of the cohort were male, i.e. 32% female (Taylor 2026 Table 1). Screened as a demographic covariate and not advanced to stepwise covariate modelling.",
      source_name = "Sex"
    )
  )

  population <- list(
    species           = "human",
    n_subjects        = 80L,
    n_studies         = 1L,
    n_administrations = 80L,
    n_observations    = 427L,
    age_range         = "adults aged >= 18 years; median 68.6 years (IQR 59.3-75.9), 64% aged 65 or older",
    age_median        = "68.6 years",
    weight_range      = "median 80.5 kg (IQR 69.7-92.3)",
    weight_median     = "80.5 kg",
    bsa_median        = "1.97 m^2 (IQR 1.80-2.14), Du Bois method",
    sex_female_pct    = 32,
    race_ethnicity    = "74 of 80 (93%) Caucasian; race was not evaluated as a covariate",
    disease_state     = "Histologically confirmed lymphoma: diffuse large B-cell lymphoma 57%, primary DLBCL of the CNS 30%, EBV-positive DLBCL of the elderly 8%, other 9%. Disease burden systemic 46%, CNS 38%, both 16%. Bone marrow involvement in 12.5%.",
    dose_range        = "Protocol dose 1.5 g/m^2 (11%), 3.5 g/m^2 (43%) or 8 g/m^2 (46%) intravenously over a 4 h infusion. Delivered dose median 7.55 g (IQR 4.83-11.25); 5.1 g (IQR 3.4-6.9) in the <= 3.5 g/m^2 stratum and 11.5 g (IQR 9.1-13.75) in the 8 g/m^2 stratum. Regimens: HDMTX 8 g/m^2 with rituximab and temozolomide (46%), HDMTX 3.5 g/m^2 with R-CHOP (36%), HDMTX monotherapy (18%).",
    renal_function    = "No patient had AKI of any stage at baseline and none were on renal replacement therapy (exclusion criteria). Five patients (6%) had a history of chronic kidney disease. Baseline serum creatinine 0.8 +/- 0.4 mg/dL, cystatin C 1.1 +/- 0.5 mg/dL; baseline CKD-EPI creatinine-cystatin C eGFR 88 +/- 24, CKD-EPI creatinine 93 +/- 27, CKD-EPI cystatin C 83 +/- 26 mL/min/1.73 m^2, Cockcroft-Gault eCrCl 99 +/- 46 mL/min (Taylor 2026 Table 1).",
    regions           = "Single center: Mayo Clinic, Rochester, MN, USA",
    notes             = paste(
      "Prospective single-center study of consecutive adults admitted for intravenous HDMTX between January 2018",
      "and December 2019 (Taylor 2026 Table 1 and Sect. 2.1). Patients were excluded if the infusion extended",
      "beyond 4 h, if they had AKI of any stage at baseline, or if they were on renal replacement therapy.",
      "One sampled administration per patient (80 patients, 80 administrations, 427 serum methotrexate",
      "concentrations), although 48% of patients were enrolled at their first methotrexate dose, 32% at their",
      "second and 20% beyond their second. Samples were drawn at predetermined times (4, 12, 24, 48, 72 and 96 h",
      "after the start of the infusion), with serum creatinine and cystatin C drawn alongside; over 90% of",
      "requested samples were obtained. All patients received protocolized hyperhydration and urine alkalinization,",
      "so urine volume and urine pH were collected but deliberately excluded from the covariate analysis.",
      "Concurrent medications were not evaluated because the treatment protocols hold known perpetrators before",
      "and during HDMTX therapy. Seventeen patients (21%) developed any-stage AKI, the earliest 9 h after the",
      "start of the infusion; the exposure-response analysis was deliberately restricted to PK data through 4 h to",
      "break the feedback loop between AKI and slowed methotrexate clearance. The model was not externally",
      "validated. NONMEM 7.5, FOCE with interaction, on log-transformed concentrations.",
      sep = " "
    )
  )

  ini({
    # Structural parameters. Typical values apply at the covariate reference
    # point eGFR = 76 mL/min/1.73 m^2 and serum albumin = 4 g/dL (Taylor 2026
    # Eq. 1). Two-compartment model, NONMEM ADVAN3 TRANS4, parameterized as
    # CL / V1 / Q / V2 (Taylor 2026 Sect. 2.2.1). RSE and bootstrap columns are
    # quoted from Table 2 so a reviewer can see the precision alongside the
    # point estimate.
    lcl <- log(11.8) ; label("Clearance (L/h)")                                  # Taylor 2026 Table 2: CL 11.8 L/h, RSE 5.1%, bootstrap median 11.8 (95% CI 10.8-12.8)
    lvc <- log(42.2) ; label("Central volume of distribution, V1 (L)")           # Taylor 2026 Table 2: V1 42.2 L, RSE 6.4%, bootstrap median 41.9 (95% CI 37.8-46.5)
    lq  <- log(0.35) ; label("Intercompartmental clearance, Q (L/h)")            # Taylor 2026 Table 2: Q 0.35 L/h, RSE 11.6%, bootstrap median 0.35 (95% CI 0.29-0.43)
    lvp <- log(6.67) ; label("Peripheral volume of distribution, V2 (L)")        # Taylor 2026 Table 2: V2 6.67 L, RSE 10.4%, bootstrap median 6.63 (95% CI 5.56-7.88)

    # Covariate effects on CL only (Taylor 2026 Eq. 1, page 470):
    #   CL_i = CL_p * (eGFR_i / 76)^theta_eGFR * (Alb_i / 4)^theta_Albumin
    # Both are POSITIVE exponents: better renal function and higher albumin
    # both predict faster methotrexate clearance. Both centering constants are
    # double-sourced: they are the denominators printed in Eq. 1, and ESM
    # Table S2 records the tested forms as '(eGFR/76)^theta' (model 8) and
    # '(Alb/4)^theta' (model 3).
    e_crcl_cl <- 0.80 ; label("Power exponent of CKD-EPI Cr-CysC eGFR (CRCL / 76) on CL (unitless)")  # Taylor 2026 Table 2: theta_eGFR 0.80, RSE 11.5%, bootstrap median 0.79 (95% CI 0.64-0.95)
    e_alb_cl  <- 0.69 ; label("Power exponent of baseline serum albumin (ALB / 4 g/dL) on CL (unitless)")  # Taylor 2026 Table 2: theta_Albumin 0.69, RSE 21.3%, bootstrap median 0.67 (95% CI 0.42-0.97)

    # Inter-individual variability. The Table 2 footnote states "IIV and
    # residual error are reported as the variance of the selected parameter",
    # so these values are omega^2 on the log scale and are used directly.
    #
    # Internal consistency check on the 0.02: Taylor 2026 Sect. 3.2.2 reports
    # the IIV on CL falling "from 28.5 to 17%" when the eGFR term was added and
    # a further "2.9%" when albumin was added, i.e. 17 - 2.9 = 14.1% in the
    # final model. sqrt(0.02) = 0.1414 = 14.1%, so the percentages in the text
    # are omega (the SD on the log scale) and Table 2's 0.02 is its square.
    #
    # Only CL and V2 carry IIV. Taylor 2026 Sect. 3.2.1: the V1 and Q omegas
    # "approached 0", a CL-V1 two-parameter block gave 82% RSE and 37%
    # shrinkage, and a full OMEGA block was not supported, so sensitivity
    # testing settled on a parsimonious CL + V2 parameterization.
    etalcl ~ 0.02  # Taylor 2026 Table 2: IIV CL variance 0.02, RSE 24.3%, shrinkage 13.3%, bootstrap median 0.02 (95% CI 0.01-0.03)
    etalvp ~ 0.03  # Taylor 2026 Table 2: IIV V2 variance 0.03, RSE 23.9%, shrinkage 31%, bootstrap median 0.03 (95% CI 0.02-0.05)

    # Residual unexplained variability. Taylor 2026 Sect. 2.2.1 specifies "a
    # proportional error structure implemented via the log-transformed both
    # sides approach". Under log-transformation of both sides the NONMEM
    # additive-on-log-scale sigma maps onto a proportional error in nlmixr2's
    # linear space, and the Table 2 footnote again reports the VARIANCE, so
    # propSd = sqrt(0.20) = 0.447.
    propSd <- sqrt(0.20); label("Proportional residual error (fraction)")  # Taylor 2026 Table 2: residual error variance 0.20, RSE 11.6%, shrinkage 12.2%, bootstrap median 0.20 (95% CI 0.16-0.23)
  })

  model({
    # Serum albumin is supplied in the canonical SI unit g/L; Taylor 2026
    # calibrated the exponent against US-convention g/dL with a reference of
    # 4 g/dL, so convert before forming the power term.
    alb_gdL <- ALB * 0.1

    # Individual parameters. Taylor 2026 Eq. 1 puts both covariates on CL only;
    # V1, Q and V2 carry no covariate. CRCL is time-varying (see
    # covariateData), so cl is re-evaluated at every record.
    cl <- exp(lcl + etalcl) * (CRCL / 76)^e_crcl_cl * (alb_gdL / 4)^e_alb_cl
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

    # Two-compartment IV disposition with first-order elimination from the
    # central compartment (NONMEM ADVAN3 mass balance). The dose is delivered
    # into `central` over the 4 h infusion via the event table's rate / dur
    # column; there is no absorption compartment.
    d/dt(central)     <- q / vp * peripheral1 - (cl + q) / vc * central
    d/dt(peripheral1) <- q / vc * central     - q / vp * peripheral1

    # Serum methotrexate concentration. The state holds umol and vc is in L,
    # so Cc is in umol/L, matching every concentration reported in Taylor 2026.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
