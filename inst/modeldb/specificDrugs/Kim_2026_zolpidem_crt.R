Kim_2026_zolpidem_crt <- function() {
  description <- paste(
    "Sequential population PK-PD model linking oral zolpidem exposure to",
    "choice reaction time (CRT) in 30 healthy Korean volunteers given a",
    "single 10 mg dose (Kim 2026). The PK layer is the companion",
    "`Kim_2026_zolpidem` model held entirely fixed (typical values, IIV and",
    "residual error); only the CRT parameters are estimated. CRT is a",
    "direct (no effect compartment) sigmoid Emax response on plasma",
    "concentration expressed as a fold-change from baseline:",
    "CRT = BASE * (1 + MAXdrug * Cp^HILL / (EC50^HILL + Cp^HILL)).",
    "MAXdrug is fixed at 5.6 because the observed CRT range spanned a",
    "6.6-fold ratio. Unlike DSST, CRT showed no meaningful learning effect",
    "and none is included. This is the most steeply concentration-dependent",
    "of the three endpoints (HILL 6) and the least potent (EC50 282 ug/L).",
    "Higher body weight, lower albumin and slower baseline CRT all increase",
    "potency; lower albumin and slower baseline CRT also steepen the",
    "response; older age raises the CRT baseline. The `effect` compartment",
    "declared in the deposited control stream is vestigial: it carries no",
    "differential equation and is not referenced, so it is not reproduced",
    "here."
  )
  reference <- paste(
    "Kim HC, Sunwoo J, Yoon S, Jang I-J, Chung J-Y. Population",
    "Pharmacokinetic-Pharmacodynamic Analysis of Zolpidem in Healthy",
    "Volunteers: Covariate Contributions to Variability in CNS Responses.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70208.",
    "doi:10.1002/psp4.70208.",
    "PK layer fixed from the companion final PK model; see",
    "modellib('Kim_2026_zolpidem')."
  )
  vignette <- "Kim_2026_zolpidem"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot   = list(analyte = "zolpidem", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "zolpidem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters EC50 as the LINEAR term (1 + -0.0193 * (WT - 60.05)),",
        "where 60.05 kg is the cohort median (Kim 2026 Table 1). Note the",
        "contrast with the companion DSST model, where weight enters IC50",
        "as a POWER term -- Table 3 footnote c prints the linear form here",
        "and Data S3 confirms it, EC50WEI = (1 + THETA(18)*(WEI - 60.05)).",
        "Reproduces Results 3.4: at the 90th weight percentile (78.0 kg)",
        "the point estimate is 0.65 [0.58-0.72], and",
        "1 - 0.0193 * (78.0 - 60.05) = 0.654. Time-fixed per subject.",
        "Source column WEI."
      ),
      source_name = "WEI"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L (canonical register unit); Kim 2026 reports g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters BOTH EC50 and HILL as power terms on (ALB_gdL / 4.4),",
        "with exponents 3.89 and -5.26 respectively; 4.4 g/dL is the",
        "cohort median (Kim 2026 Table 1). The canonical register unit for",
        "ALB is g/L, so model() applies the inline conversion",
        "alb_gdL <- ALB * 0.1 required by the register entry, keeping both",
        "exponents aligned with their original g/dL calibration; a",
        "model-data column of 44 g/L is the reference subject. Because both",
        "terms are ratios, the exponents are unchanged by the conversion.",
        "Reproduces Results 3.4 for HILL: (4.0/4.4)^-5.26 = 1.65 vs the",
        "reported 1.6 [1.3-2.1] at the 10th percentile, and",
        "(4.8/4.4)^-5.26 = 0.633 vs the reported 0.65 [0.53-0.77] at the",
        "90th. Time-fixed per subject."
      ),
      source_name = "ALB"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the CRT baseline as the linear term",
        "(1 + 0.0112 * (AGE - 29.5)), where 29.5 years is the cohort",
        "median (Kim 2026 Table 1). Older subjects have a slower (higher)",
        "baseline CRT -- note the sign is opposite to the DSST baseline,",
        "as expected for a latency rather than a score. Time-fixed per",
        "subject."
      ),
      source_name = "AGE"
    ),
    CRT_BL = list(
      description = "Baseline choice reaction time measured at 1 h on Day -1",
      units = "msec",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters BOTH EC50 (power term (CRT_BL / 438.5)^-1.28) and HILL",
        "(linear term 1 + 0.00637 * (CRT_BL - 438.5)); 438.5 msec is the",
        "cohort median (Kim 2026 Table 1). Data S3 $INPUT names the column",
        "CRTB1, 'CRT baseline 1h'. This is a per-subject observed baseline",
        "readout carried in as data, distinct from the estimated BASE",
        "parameter of this same model. Reproduces Results 3.4 for HILL: at",
        "the 90th percentile (504 msec) the reported point estimate is",
        "1.4 [1.3-1.6] and 1 + 0.00637 * (504 - 438.5) = 1.417. Slower",
        "baseline responders are both more sensitive and more steeply",
        "responsive, which the Discussion frames as reduced",
        "neurophysiological reserve. Time-fixed per subject."
      ),
      source_name = "CRTB1"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Screened in the stepwise covariate analysis (Kim 2026 Methods",
        "2.4, Tables S3 and S4) but not retained on any CRT parameter. No",
        "point estimate reported. Source column SEX coded 1 = M, 2 = F."
      ),
      source_name = "SEX"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened (Kim 2026 Table S1) but not retained; no point estimate reported.",
      source_name = "ALT"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened (Kim 2026 Table S1) but not retained; no point estimate reported.",
      source_name = "AST"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL as reported by Kim 2026 (canonical register unit is umol/L)",
      type = "continuous",
      notes = "Screened (Kim 2026 Table S1) but not retained; no point estimate reported.",
      source_name = "TBIL"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened (Kim 2026 Table S1) but not retained; no point estimate reported.",
      source_name = "CREA"
    ),
    DSST_BL = list(
      description = "Baseline digit symbol substitution test score at 1 h on Day -1",
      units = "score",
      type = "continuous",
      notes = "Screened (Kim 2026 Table S1) but not retained on any CRT parameter; no point estimate reported.",
      source_name = "DSSTB1"
    ),
    VAS_SEDATION_BL = list(
      description = "Baseline sedation visual analog scale at 1 h on Day -1",
      units = "mm",
      type = "continuous",
      notes = "Screened (Kim 2026 Table S1) but not retained; no point estimate reported.",
      source_name = "VASB1"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 30,
    n_studies = 1,
    age_range = "20-44 years",
    age_median = "29.5 years",
    weight_range = "50.0-83.0 kg",
    weight_median = "60.05 kg",
    sex_female_pct = 50,
    race_ethnicity = c(Asian = 100),
    disease_state = "healthy volunteers",
    dose_range = "single oral 10 mg immediate-release zolpidem",
    regions = "Korea",
    n_observations = "325 plasma zolpidem concentrations and 388 CRT observations",
    laboratory = paste(
      "Median (range) at baseline (Kim 2026 Table 1): albumin 4.4",
      "(3.9-5.0) g/dL; CRT at 1 h on Day -1 438.5 (372-709) msec."
    ),
    notes = paste(
      "CRT was measured at -23, -22, -20 and -12 h on Day -1 (baseline)",
      "and at 0.5, 1, 1.5, 2, 3, 4, 6, 8 and 12 h on Day 1 (Kim 2026",
      "Methods 2.1). Results 3.2.3: the observed CRT range was 337 to",
      "2219 msec, a 6.6-fold ratio, which is what MAXdrug = 5.6 was fixed",
      "from."
    )
  )

  ini({
    # ====================================================================
    # PHARMACOKINETIC LAYER -- FIXED FROM THE FINAL PK MODEL
    # ====================================================================
    # Kim 2026 Results 3.2: each PD model was linked sequentially to the
    # final PK model. Data S3 $THETA(1)-$THETA(8) and $OMEGA(1,1)-
    # $OMEGA(6,6) all carry the FIX flag at the Data S1 / Table 2 final
    # values. See modellib("Kim_2026_zolpidem") for the estimating run.

    lka <- fixed(log(11.7))
    label("Absorption rate constant, depot -> central (Ka, 1/h)")                    # Kim 2026 Table 2: Ka = 11.7 1/h. Data S3 $THETA(3) = "(11.7) FIX".

    lcl <- fixed(log(18.0))
    label("Apparent clearance (CL/F, L/h)")                                          # Kim 2026 Table 2: CL/F = 18.0 L/h. Data S3 $THETA(1) = "(18) FIX".

    lvc <- fixed(log(64.0))
    label("Apparent central volume of distribution (Vc/F, L)")                       # Kim 2026 Table 2: Vc/F = 64.0 L. Data S3 $THETA(2) = "(64) FIX".

    lmtt <- fixed(log(0.25))
    label("Mean transit time of the absorption chain (MTT, h)")                      # Kim 2026 Table 2: MTT = 0.25 h. Data S3 $THETA(4) = "(0.25) FIX".

    lnn <- fixed(log(19.4))
    label("Number of transit compartments (NN, unitless)")                           # Kim 2026 Table 2: NN = 19.4. Data S3 $THETA(5) = "(19.4) FIX".

    lfdepot <- fixed(log(1))
    label("Relative bioavailability (BIO, fraction)")                                # Kim 2026 Table 2: BIO = "1 Fixed". Data S3 $THETA(6) = "(1) FIX".

    etalka ~ fixed(2.48)
    # Kim 2026 Table 2: IIV Ka = 330.8 %CV. Data S3 $OMEGA(3,3) = "2.48 FIX".

    etalcl ~ fixed(0.0565)
    # Kim 2026 Table 2: IIV CL/F = 24.1 %CV. Data S3 $OMEGA(1,1) = "0.0565 FIX".

    etalmtt ~ fixed(0.133)
    # Kim 2026 Table 2: IIV MTT = 37.7 %CV. Data S3 $OMEGA(4,4) = "0.133 FIX".

    etalfdepot ~ fixed(0.0647)
    # Kim 2026 Table 2: IIV BIO = 25.9 %CV. Data S3 $OMEGA(6,6) = "0.0647 FIX".

    propSd <- fixed(0.189)
    label("Proportional residual error on plasma zolpidem (fraction)")               # Kim 2026 Table 2: proportional residual error = 18.9%. Data S3 $THETA(7) = "(0.189) FIX".

    addSd <- fixed(0.00001)
    label("Additive residual error on plasma zolpidem (ug/L)")                       # Data S3 $THETA(8) = "(0.00001) FIX"; a numerical guard, not reported in Kim 2026 Table 2.

    # ====================================================================
    # CRT RESPONSE PARAMETERS
    # ====================================================================
    # Typical values are the Final model column of Kim 2026 Table 3,
    # Choice reaction time (CRT) block, reproduced in the Data S3 $THETA
    # block.

    lec50 <- log(282)
    label("Half-maximal effective concentration for CRT (EC50, ug/L)")               # Kim 2026 Table 3: EC50 = 282 ug/L (RSE 8.5%; SIR median 286, 95% CI 252-328). Data S3 $THETA(9) "C_EC50". Base model without covariates was 209.

    lhill <- log(6)
    label("Hill coefficient of the CRT concentration-response (HILL, unitless)")     # Kim 2026 Table 3: HILL = 6 (RSE 8.8%; SIR median 5.97, 95% CI 5.17-6.98). Data S3 $THETA(12) "C_HILL". Base model without covariates was 7.73. Results 3.2.3 notes the CRT relationship is the steepest of the three endpoints.

    lrbase <- log(437)
    label("Baseline choice reaction time (BASE, msec)")                              # Kim 2026 Table 3: BASE = 437 msec (RSE 1.8%; SIR median 437, 95% CI 423-452). Data S3 $THETA(10) "C_CRT0".

    lemax <- fixed(log(5.6))
    label("Maximum fractional drug effect on CRT (MAXdrug, unitless)")               # Kim 2026 Table 3: MAXdrug = "5.6 Fixed". Data S3 $THETA(11) = "(5.6) FIX". Results 3.2.3: "The MAXdrug parameter was fixed at 5.6, based on the observation that the maximum CRT value (2219 ms) was 6.6-fold higher than the minimum CRT value (337 ms)" -- i.e. a 6.6-fold maximum corresponds to a fractional increase of 5.6 over baseline.

    # ====================================================================
    # COVARIATE EFFECTS ON THE CRT PARAMETERS
    # ====================================================================
    # Data S3 $PK composes them multiplicatively:
    #   EC50COV = EC50ALB * EC50CRTB1 * EC50WEI
    #   HILLCOV = HILLALB * HILLCRTB1
    #   CRT0COV = CRT0AGE
    # Kim 2026 Table 3 footnotes c, d and e print exactly these forms.

    e_wt_ec50 <- -0.0193
    label("Linear coefficient on (WT - 60.05 kg) for CRT EC50 (per kg)")             # Kim 2026 Table 3 "Effect of weight on EC50" = -0.0193 (RSE 14.1%; SIR median -0.0196, 95% CI -0.0242 to -0.0145). Data S3 $THETA(18) "EC50WEI1": EC50WEI = (1 + THETA(18)*(WEI - 60.05)). Footnote c prints the same linear form.

    e_alb_ec50 <- 3.89
    label("Power exponent on (ALB / 4.4 g/dL) for CRT EC50 (unitless)")              # Kim 2026 Table 3 "Effect of albumin on EC50" = 3.89 (RSE 29.8%; SIR median 3.81, 95% CI 1.64-5.98). Data S3 $THETA(16) "EC50ALB1": EC50ALB = ((ALB/4.4)**THETA(16)). Footnote c.

    e_crt_bl_ec50 <- -1.28
    label("Power exponent on (CRT_BL / 438.5 msec) for CRT EC50 (unitless)")         # Kim 2026 Table 3 "Effect of CRTB1 on EC50" = -1.28 (RSE 26.1%; SIR median -1.33, 95% CI -2.05 to -0.64). Data S3 $THETA(17) "EC50CRTB11": EC50CRTB1 = ((CRTB1/438.5)**THETA(17)). Footnote c.

    e_alb_hill <- -5.26
    label("Power exponent on (ALB / 4.4 g/dL) for CRT HILL (unitless)")              # Kim 2026 Table 3 "Effect of albumin on HILL" = -5.26 (RSE 21.7%; SIR median -5.13, 95% CI -7.23 to -2.96). Data S3 $THETA(19) "HILLALB1": HILLALB = ((ALB/4.4)**THETA(19)). Footnote d.

    e_crt_bl_hill <- 0.00637
    label("Linear coefficient on (CRT_BL - 438.5 msec) for CRT HILL (per msec)")     # Kim 2026 Table 3 "Effect of CRTB1 on HILL" = 0.00637 (RSE 15.7%; SIR median 0.00644, 95% CI 0.00420-0.00840). Data S3 $THETA(20) "HILLCRTB11": HILLCRTB1 = (1 + THETA(20)*(CRTB1 - 438.5)). Footnote d.

    e_age_rbase <- 0.0112
    label("Linear coefficient on (AGE - 29.5 years) for CRT BASE (per year)")        # Kim 2026 Table 3 "Effect of age on BASE" = 0.0112 (RSE 22.2%; SIR median 0.0112, 95% CI 0.0059-0.0163). Data S3 $THETA(15) "CRT0AGE1": CRT0AGE = (1 + THETA(15)*(AGE - 29.5)). Footnote e: BASE_CRT (msec) = 437 x (1 + 0.0112 x (Age - 29.5)).

    # ====================================================================
    # CRT INTER-INDIVIDUAL VARIABILITY AND RESIDUAL ERROR
    # ====================================================================
    # Data S3 $OMEGA(9,9) (MAXdrug) and $OMEGA(10,10) (HILL) are both
    # "0 FIX", so no eta appears for those two parameters. Results 3.2.3:
    # albumin and CRTB1 together "contributed to a 49.9% reduction in the
    # IIV of HILL, which corresponds to the IIV estimated in the base
    # PK-CRT model before covariate inclusion" -- i.e. the covariates
    # explained the HILL IIV entirely, which is why Kim 2026 Table 3 shows
    # 49.9 %CV for the base model and "-" for the final model.

    etalec50 ~ 0.0436
    # Kim 2026 Table 3: IIV EC50 = 21.1 %CV (RSE 35.1%, shrinkage 20.2%). Data S3 $OMEGA(7,7) = 0.0436; 100*sqrt(exp(0.0436)-1) = 21.1%. Results 3.2.3: the three EC50 covariates reduced this IIV by 7.1% relative to the base model's 28.2 %CV.

    etalrbase ~ 0.00871
    # Kim 2026 Table 3: IIV BASE = 9.4 %CV (RSE 27.6%, shrinkage 4.4%). Data S3 $OMEGA(8,8) = 0.00871; 100*sqrt(exp(0.00871)-1) = 9.35%. Results 3.2.3: age reduced this IIV by 1.7% relative to the base model's 11.1 %CV.

    propSd_CRT <- 0.0942
    label("Proportional residual error on CRT (fraction)")                           # Kim 2026 Table 3: proportional residual error = 9.42% (RSE 7.2%; SIR median 9.50, 95% CI 8.88-10.23). Data S3 $THETA(13) "C_Prop.RE".

    addSd_CRT <- fixed(0.00001)
    label("Additive residual error on CRT (msec)")                                   # Data S3 $THETA(14) "C_Add.RE" = "(0, 0.00001) FIX". Results 3.2.3 describes the CRT residual error as proportional; the fixed 1e-5 additive term is a numerical guard and is not reported in Kim 2026 Table 3.
  })

  model({
    # --- 1. Derived covariate terms ---------------------------------------
    # The canonical ALB column is in SI g/L; Kim 2026 calibrated both
    # albumin exponents and the 4.4 reference against US-convention g/dL,
    # so the register-mandated inline conversion is applied here.
    alb_gdL <- ALB * 0.1

    cov_ec50  <- (alb_gdL / 4.4)^e_alb_ec50 *
      (CRT_BL / 438.5)^e_crt_bl_ec50 *
      (1 + e_wt_ec50 * (WT - 60.05))
    cov_hill  <- (alb_gdL / 4.4)^e_alb_hill *
      (1 + e_crt_bl_hill * (CRT_BL - 438.5))
    cov_rbase <- 1 + e_age_rbase * (AGE - 29.5)

    # --- 2. Individual parameters -----------------------------------------
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc)
    mtt    <- exp(lmtt + etalmtt)
    nn     <- exp(lnn)
    fdepot <- exp(lfdepot + etalfdepot)

    ec50  <- exp(lec50 + etalec50) * cov_ec50
    hill  <- exp(lhill) * cov_hill
    rbase <- exp(lrbase + etalrbase) * cov_rbase
    emax  <- exp(lemax)

    # --- 3. Micro-constants ------------------------------------------------
    kel <- cl / vc

    # --- 4. ODE system ------------------------------------------------------
    # Identical to the companion PK model; see modellib("Kim_2026_zolpidem").
    d/dt(depot)   <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - kel * central

    # --- 5. Bioavailability -------------------------------------------------
    f(depot) <- 0

    # --- 6. Observations and residual error ---------------------------------
    Cc <- 1000 * central / vc

    # Data S3 $ERROR:
    #   EDRUG = CRT0 * (1 + CRTM * (CP**HILL) / (EC50**HILL + CP**HILL))
    #   CRT   = EDRUG
    # No learning term: Results 3.2.3 reports that adding one improved the
    # OFV by only 3.887, below the retention threshold.
    CRT <- rbase * (1 + emax * Cc^hill / (ec50^hill + Cc^hill))

    Cc  ~ prop(propSd) + add(addSd)
    CRT ~ prop(propSd_CRT) + add(addSd_CRT)
  })
}
