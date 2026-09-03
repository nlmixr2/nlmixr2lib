Kim_2026_zolpidem_dsst <- function() {
  description <- paste(
    "Sequential population PK-PD model linking oral zolpidem exposure to",
    "the digit symbol substitution test (DSST) in 30 healthy Korean",
    "volunteers given a single 10 mg dose (Kim 2026). The PK layer is the",
    "companion `Kim_2026_zolpidem` model held entirely fixed (typical",
    "values, IIV and residual error); only the DSST parameters are",
    "estimated. DSST is a direct (no effect compartment) sigmoid Imax",
    "response on plasma concentration, multiplied onto a baseline that",
    "itself rises over the study with a two-parameter hyperbolic learning",
    "(practice) effect: DSST = (BASE + MAXlearn * TIME / (TIME + LRPAR)) *",
    "(1 - Cp^HILL / (IC50^HILL + Cp^HILL)). TIME is absolute study time,",
    "not time after dose -- the authors relabelled the first Day -1",
    "assessment as t = 1 h so that dosing falls at t = 24 h, at which point",
    "about 85% of the maximal learning effect has accrued. Higher body",
    "weight and lower albumin both increase potency (lower IC50), and older",
    "age lowers the DSST baseline. The `effect` compartment declared in the",
    "deposited control stream is vestigial: it carries no differential",
    "equation and is not referenced, so it is not reproduced here."
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
        "Enters IC50 as the power term (WT / 60.05)^-1.57, where 60.05 kg",
        "is the cohort median (Kim 2026 Table 1). Time-fixed per subject.",
        "IMPORTANT -- Kim 2026 Table 3 footnote a prints this term as",
        "'205 x (Weight/60.05) x (1 + 0.541 x (Albumin - 4.4))' with the",
        "exponent lost in typesetting; the orphaned '-1.57' that trails the",
        "Table 3 abbreviations block is that dropped superscript. Data S2",
        "gives the true form, IC50WEI = ((WEI/60.05)**THETA(18)) with",
        "THETA(18) = -1.57, and Table 3's own 'Effect of weight on IC50'",
        "row reports -1.57. The power form is what reproduces the paper's",
        "forest plot: Results 3.4 gives 0.66 [0.59-0.74] at the 90th",
        "weight percentile (78.0 kg), and (78.0/60.05)^-1.57 = 0.663,",
        "whereas the linear form printed in the footnote would give 1.30.",
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
        "Enters IC50 as the linear term (1 + 0.541 * (ALB_gdL - 4.4)),",
        "where 4.4 g/dL is the cohort median (Kim 2026 Table 1). The",
        "canonical register unit for ALB is g/L, so model() applies the",
        "inline conversion alb_gdL <- ALB * 0.1 required by the register",
        "entry, keeping the 0.541 coefficient aligned with its original",
        "g/dL calibration. A model-data column of 44 g/L is therefore the",
        "reference subject. Lower albumin lowers IC50 (greater potency),",
        "which the Discussion attributes to the free-drug hypothesis given",
        "zolpidem's 92.5% plasma protein binding. Time-fixed per subject."
      ),
      source_name = "ALB"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the DSST baseline as the linear term",
        "(1 + -0.0143 * (AGE - 29.5)), where 29.5 years is the cohort",
        "median (Kim 2026 Table 1). Older subjects have a lower baseline",
        "DSST score. Time-fixed per subject."
      ),
      source_name = "AGE"
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
        "2.4, Tables S3 and S4) but not retained on any DSST parameter.",
        "The paper's premise is that the apparent sex differences are",
        "mediated by weight, albumin and age rather than by sex itself.",
        "No point estimate reported. Source column SEX coded 1 = M, 2 = F."
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
      notes = paste(
        "Screened as a candidate covariate on the DSST parameters (Kim",
        "2026 Methods 2.4, Table S1; control-stream $INPUT column DSSTB1,",
        "'DSST baseline 1h') but not retained -- the DSST baseline is",
        "estimated as a model parameter (BASE) with age as its covariate",
        "rather than carried in as an observed baseline. Cohort median 72.5",
        "(range 46-104). No point estimate reported."
      ),
      source_name = "DSSTB1"
    ),
    CRT_BL = list(
      description = "Baseline choice reaction time at 1 h on Day -1",
      units = "msec",
      type = "continuous",
      notes = paste(
        "Screened (Kim 2026 Table S1) but not retained on any DSST",
        "parameter; it IS retained on EC50 and HILL in the companion",
        "PK-CRT model. No point estimate reported."
      ),
      source_name = "CRTB1"
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
    n_observations = "325 plasma zolpidem concentrations and 388 DSST observations",
    laboratory = paste(
      "Median (range) at baseline (Kim 2026 Table 1): albumin 4.4",
      "(3.9-5.0) g/dL; DSST at 1 h on Day -1 72.5 (46-104) score."
    ),
    notes = paste(
      "DSST was measured at -23, -22, -20 and -12 h on Day -1 (baseline)",
      "and at 0.5, 1, 1.5, 2, 3, 4, 6, 8 and 12 h on Day 1 (Kim 2026",
      "Methods 2.1). For the analysis the -23 h Day -1 time point was",
      "redefined as t = 1 h, placing the dose at t = 24 h; the learning",
      "term is a function of that absolute study time."
    )
  )

  ini({
    # ====================================================================
    # PHARMACOKINETIC LAYER -- FIXED FROM THE FINAL PK MODEL
    # ====================================================================
    # Kim 2026 Results 3.2: "The population PK-PD model was developed
    # sequentially by linking each PD model for the DSST, CRT, and VAS to
    # the final PK model." Data S2 $THETA(1)-$THETA(8) and $OMEGA(1,1)-
    # $OMEGA(6,6) all carry the FIX flag at the Data S1 / Table 2 final
    # values, so every PK quantity below is fixed rather than estimated.
    # See modellib("Kim_2026_zolpidem") for the estimating run.

    lka <- fixed(log(11.7))
    label("Absorption rate constant, depot -> central (Ka, 1/h)")                    # Kim 2026 Table 2: Ka = 11.7 1/h. Data S2 $THETA(3) = "(11.7) FIX".

    lcl <- fixed(log(18.0))
    label("Apparent clearance (CL/F, L/h)")                                          # Kim 2026 Table 2: CL/F = 18.0 L/h. Data S2 $THETA(1) = "(18) FIX".

    lvc <- fixed(log(64.0))
    label("Apparent central volume of distribution (Vc/F, L)")                       # Kim 2026 Table 2: Vc/F = 64.0 L. Data S2 $THETA(2) = "(64) FIX".

    lmtt <- fixed(log(0.25))
    label("Mean transit time of the absorption chain (MTT, h)")                      # Kim 2026 Table 2: MTT = 0.25 h. Data S2 $THETA(4) = "(0.25) FIX".

    lnn <- fixed(log(19.4))
    label("Number of transit compartments (NN, unitless)")                           # Kim 2026 Table 2: NN = 19.4. Data S2 $THETA(5) = "(19.4) FIX".

    lfdepot <- fixed(log(1))
    label("Relative bioavailability (BIO, fraction)")                                # Kim 2026 Table 2: BIO = "1 Fixed". Data S2 $THETA(6) = "(1) FIX".

    etalka ~ fixed(2.48)
    # Kim 2026 Table 2: IIV Ka = 330.8 %CV. Data S2 $OMEGA(3,3) = "2.48 FIX".

    etalcl ~ fixed(0.0565)
    # Kim 2026 Table 2: IIV CL/F = 24.1 %CV. Data S2 $OMEGA(1,1) = "0.0565 FIX".

    etalmtt ~ fixed(0.133)
    # Kim 2026 Table 2: IIV MTT = 37.7 %CV. Data S2 $OMEGA(4,4) = "0.133 FIX".

    etalfdepot ~ fixed(0.0647)
    # Kim 2026 Table 2: IIV BIO = 25.9 %CV. Data S2 $OMEGA(6,6) = "0.0647 FIX".

    propSd <- fixed(0.189)
    label("Proportional residual error on plasma zolpidem (fraction)")               # Kim 2026 Table 2: proportional residual error = 18.9%. Data S2 $THETA(7) = "(0.189) FIX".

    addSd <- fixed(0.00001)
    label("Additive residual error on plasma zolpidem (ug/L)")                       # Data S2 $THETA(8) = "(0.00001) FIX"; a numerical guard rather than an estimated component, and not reported in Kim 2026 Table 2.

    # ====================================================================
    # DSST RESPONSE PARAMETERS
    # ====================================================================
    # Typical values are the Final model column of Kim 2026 Table 3
    # ("Population pharmacodynamic parameter estimates and sampling
    # importance resampling results"), Digit symbol substitution test
    # (DSST) block, reproduced in the Data S2 $THETA block.

    lic50 <- log(205)
    label("Half-maximal inhibitory concentration for DSST (IC50, ug/L)")             # Kim 2026 Table 3: IC50 = 205 ug/L (RSE 7.5%; SIR median 204, 95% CI 183-232). Data S2 $THETA(9) "D_IC50". Base model without covariates was 188.

    lhill <- log(2.54)
    label("Hill coefficient of the DSST concentration-response (HILL, unitless)")    # Kim 2026 Table 3: HILL = 2.54 (RSE 11.1%; SIR median 2.55, 95% CI 2.17-2.99). Data S2 $THETA(10) "D_HILL".

    lrbase <- log(68.2)
    label("Baseline DSST score before any learning effect (BASE, score)")            # Kim 2026 Table 3: BASE = 68.2 score (RSE 3.3%; SIR median 67.9, 95% CI 63.6-72.8). Data S2 $THETA(11) "D_BASE".

    lemax_learn <- log(14.4)
    label("Maximum learning (practice) effect on DSST (MAXlearn, score)")            # Kim 2026 Table 3: MAXlearn = 14.4 score (RSE 10.6%; SIR median 14.8, 95% CI 11.6-17.5). Data S2 $THETA(12) "D_MAX".

    lt50_learn <- log(4.3)
    label("Time to half-maximal learning effect on DSST (LRPAR, h)")                 # Kim 2026 Table 3: LRPAR = 4.3 h (RSE 38.6%; SIR median 4.29, 95% CI 1.62-9.22). Data S2 $THETA(13) "D_LRPAR". Results 3.2.2 / Discussion: 4.3 h implies "approximately 85% of the maximal learning effect would be achieved within 24 h" -- 24/(24 + 4.3) = 0.848, which confirms that the learning term runs on absolute study time.

    # ====================================================================
    # COVARIATE EFFECTS ON THE DSST PARAMETERS
    # ====================================================================

    e_wt_ic50 <- -1.57
    label("Power exponent on (WT / 60.05) for DSST IC50 (unitless)")                 # Kim 2026 Table 3 "Effect of weight on IC50" = -1.57 (RSE 24.3%; SIR median -1.57, 95% CI -2.12 to -1.01). Data S2 $THETA(18) "IC50WEI1", entering as IC50WEI = ((WEI/60.05)**THETA(18)). Table 3 footnote a drops this exponent in typesetting; see covariateData[[WT]]$notes.

    e_alb_ic50 <- 0.541
    label("Linear coefficient on (ALB - 4.4 g/dL) for DSST IC50 (per g/dL)")         # Kim 2026 Table 3 "Effect of albumin on IC50" = 0.541 (RSE 41.6%; SIR median 0.520, 95% CI 0.183-0.893). Data S2 $THETA(17) "IC50ALB1", entering as IC50ALB = (1 + THETA(17)*(ALB - 4.4)). Also printed in Table 3 footnote a.

    e_age_rbase <- -0.0143
    label("Linear coefficient on (AGE - 29.5 years) for DSST BASE (per year)")       # Kim 2026 Table 3 "Effect of age on BASE" = -0.0143 (RSE 23.1%; SIR median -0.0143, 95% CI -0.0211 to -0.00697). Data S2 $THETA(16) "DSSTBAGE1", entering as DSSTBAGE = (1 + THETA(16)*(AGE - 29.5)). Table 3 footnote b: BASE_DSST (score) = 68.2 x (1 - 0.0143 x (Age - 29.5)).

    # ====================================================================
    # DSST INTER-INDIVIDUAL VARIABILITY AND RESIDUAL ERROR
    # ====================================================================
    # Data S2 $OMEGA(8,8) (HILL), $OMEGA(10,10) (MAXlearn) and
    # $OMEGA(11,11) (LRPAR) are all "0 FIX", so no eta appears for those
    # three parameters.

    etalic50 ~ 0.0484
    # Kim 2026 Table 3: IIV IC50 = 22.3 %CV (RSE 28.9%, shrinkage 15.2%). Data S2 $OMEGA(7,7) = 0.0484; 100*sqrt(exp(0.0484)-1) = 22.3%. Results 3.2.2: weight and albumin reduced this IIV by 11.5% relative to the base model's 33.8 %CV.

    etalrbase ~ 0.0171
    # Kim 2026 Table 3: IIV BASE = 13.1 %CV (RSE 28.9%, shrinkage 3%). Data S2 $OMEGA(9,9) = 0.0171; 100*sqrt(exp(0.0171)-1) = 13.1%. Results 3.2.2: age reduced this IIV by 2.7% relative to the base model's 15.8 %CV.

    propSd_DSST <- fixed(0.00001)
    label("Proportional residual error on DSST (fraction)")                          # Data S2 $THETA(14) "D_Prop.RE" = "(0, 0.00001) FIX". Results 3.2.2 describes the DSST residual error as additive; the fixed 1e-5 proportional term is a numerical guard and is not reported in Kim 2026 Table 3.

    addSd_DSST <- 6
    label("Additive residual error on DSST (score)")                                 # Kim 2026 Table 3: additive residual error = 6 score (RSE 5.8%; SIR median 6.04, 95% CI 5.58-6.53). Data S2 $THETA(15) "D_Add.RE".
  })

  model({
    # --- 1. Derived covariate terms ---------------------------------------
    # The canonical ALB column is in SI g/L; Kim 2026 calibrated the 0.541
    # coefficient and the 4.4 reference against US-convention g/dL, so the
    # register-mandated inline conversion is applied here.
    alb_gdL <- ALB * 0.1

    # Data S2 $PK:
    #   IC50WEI = ((WEI/60.05)**THETA(18))
    #   IC50ALB = (1 + THETA(17)*(ALB - 4.4))
    #   IC50COV = IC50ALB * IC50WEI
    #   DSSTBAGE = (1 + THETA(16)*(AGE - 29.5))
    cov_ic50  <- (WT / 60.05)^e_wt_ic50 * (1 + e_alb_ic50 * (alb_gdL - 4.4))
    cov_rbase <- 1 + e_age_rbase * (AGE - 29.5)

    # --- 2. Individual parameters -----------------------------------------
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc)
    mtt    <- exp(lmtt + etalmtt)
    nn     <- exp(lnn)
    fdepot <- exp(lfdepot + etalfdepot)

    ic50       <- exp(lic50 + etalic50) * cov_ic50
    hill       <- exp(lhill)
    rbase      <- exp(lrbase + etalrbase) * cov_rbase
    emax_learn <- exp(lemax_learn)
    t50_learn  <- exp(lt50_learn)

    # --- 3. Micro-constants ------------------------------------------------
    kel <- cl / vc

    # --- 4. ODE system ------------------------------------------------------
    # Identical to the companion PK model; see modellib("Kim_2026_zolpidem")
    # for the derivation of the Savic transit input from the Data S1 $DES
    # block.
    d/dt(depot)   <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - kel * central

    # --- 5. Bioavailability -------------------------------------------------
    f(depot) <- 0

    # --- 6. Observations and residual error ---------------------------------
    Cc <- 1000 * central / vc

    # Data S2 $ERROR:
    #   EDRUG  = (CP**HILL) / (IC50**HILL + CP**HILL)
    #   ELEARN = DSSTB + DSSTM * (TIME / (TIME + LRPAR))
    #   DSST   = ELEARN * (1 - EDRUG)
    # TIME is absolute study time. Kim 2026 Methods 2.1 relabels the first
    # Day -1 assessment (-23 h) as t = 1 h so that dosing occurs at
    # t = 24 h; event tables driving this model must place the dose there
    # for the learning term to be on the scale the authors estimated.
    edrug  <- Cc^hill / (ic50^hill + Cc^hill)
    elearn <- rbase + emax_learn * time / (time + t50_learn)
    DSST   <- elearn * (1 - edrug)

    Cc   ~ prop(propSd) + add(addSd)
    DSST ~ prop(propSd_DSST) + add(addSd_DSST)
  })
}
