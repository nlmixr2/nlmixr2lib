Zou_2025_edoxaban <- function() {
  description <- paste(
    "Joint population PK / PD model for oral edoxaban in pediatric patients",
    "(0 to <18 years) with venous thromboembolism or cardiac disease at risk",
    "of thromboembolic events (Zou 2025; pooled phase 1 PK/PD study,",
    "Hokusai-VTE PEDIATRICS and ENNOBLE-ATE). Disposition is two-compartment",
    "with linear elimination; absorption is a chain of 15 transit",
    "compartments emptying at a common rate Ktr followed by a first-order",
    "step Ka into the central compartment. Apparent clearance and",
    "inter-compartmental clearance are allometrically scaled to body weight",
    "(exponent fixed to 0.75), apparent central and peripheral volumes to",
    "body weight (exponent fixed to 1). Clearance additionally carries a",
    "power term on bedside-Schwartz eGFR (reference 110 mL/min/1.73 m^2) and",
    "a Rhodin-style postmenstrual-age renal maturation Hill function with",
    "TM50 and Hill both fixed (47.7 weeks, 3.40). Three direct-response",
    "pharmacodynamic endpoints are driven by the plasma concentration: an",
    "Emax model for anti-factor Xa activity (baseline fixed to 0.1 IU/mL)",
    "and linear models for activated partial thromboplastin time and",
    "prothrombin time, each with variability on its own baseline or maximum",
    "effect. The paper's PD layers were fitted sequentially against OBSERVED",
    "edoxaban concentrations; they are wired to the model-predicted",
    "concentration here so the published PK/PD relationships can be",
    "simulated as one system.",
    sep = " "
  )
  reference <- paste(
    "Zou P, Atluri A, Chang P, Goedecke M, Leil TA.",
    "Population pharmacokinetics and pharmacodynamics of edoxaban in",
    "pediatric patients. CPT Pharmacometrics Syst Pharmacol.",
    "2025;14(1):118-129. doi:10.1002/psp4.13248",
    sep = " "
  )
  vignette <- "Zou_2025_edoxaban"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    transit4    = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    transit5    = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    transit6    = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    transit7    = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    transit8    = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    transit9    = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    transit10   = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    transit11   = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    transit12   = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    transit13   = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    transit14   = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    transit15   = list(analyte = "edoxaban", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "edoxaban", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "edoxaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight at the corresponding visit",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling with reference 70 kg. Supplement NONMEM $PK block:",
        "TVBW = 70; ALLMCL_BW = (BW/TVBW)**THETA(9) with THETA(9) = 0.75 FIX,",
        "applied to CL and Q; ALLMV_BW = (BW/TVBW)**THETA(10) with",
        "THETA(10) = 1 FIX, applied to Vc and Vp. Zou 2025 Results:",
        "'The scaling exponents of Vc/F and Vp/F with body weight were",
        "estimated to be not significantly different from 1; therefore, they",
        "were fixed to 1. The estimated scaling exponent of CL/F and Q/F did",
        "not substantially improve model fit and the exponent was thus fixed",
        "to 0.75.' Cohort median 21.1 kg (range 2.60-157; Table S6).",
        sep = " "
      ),
      source_name        = "BW"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age at birth plus postnatal age)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "WEEKS, not the register-default months: the Rhodin 2009 maturation",
        "constants this model uses (TM50 = 47.7, Hill = 3.40) are only",
        "meaningful on the week scale, matching the Germovsek_2018_meropenem",
        "and Riccobene_2017_ceftaroline precedents. Drives the maturation",
        "factor FPMA = PMA^3.4 / (47.7^3.4 + PMA^3.4) on CL/F only (supplement",
        "NONMEM $PK: TVCL = THETA(3)*CLEGFR*ALLMCL_BW*FPMA). Note this is the",
        "BARE Hill fraction, NOT normalised to a reference PMA, so CL/F_TYP",
        "is the fully matured value approached as PMA grows large.",
        "Cohort median 353 weeks (range 38.6-970; Table S6).",
        sep = " "
      ),
      source_name        = "PMA"
    ),
    CRCL = list(
      description        = "Bedside-Schwartz estimated glomerular filtration rate, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Bedside Schwartz formula, given in Table S3 footnote a as",
        "eGFR (mL/min/1.73 m^2) = 0.413 x height (cm) / serum creatinine (mg/dL);",
        "Zou 2025 Methods 'Renal function measure' cites it as Equation 1 and",
        "states eGFR was the renal-function covariate for both the PopPK and",
        "PopPK/PD analyses. Enters CL/F as the power term (CRCL/110)^0.268",
        "(supplement NONMEM $PK: TVEGFR = 110;",
        "CLEGFR = (EGFR/TVEGFR)**THETA(11)). Cohort median",
        "111 mL/min/1.73 m^2 (range 29.5-774; Table S6). The Schwartz",
        "creatinine clearance of Methods Equation 2 is a DIFFERENT quantity,",
        "used only to compare pediatric with adult renal function in the",
        "Discussion, and is not a model covariate.",
        sep = " "
      ),
      source_name        = "EGFR"
    )
  )

  # Screened in the covariate search and not retained in the final model.
  # Documentation only -- none of these is referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at the corresponding visit",
      units       = "years",
      type        = "continuous",
      notes       = "Tested on volume and clearance (Table S3) and on PD baseline / slope / maximum effect (Table S4); not retained. Body size and postmenstrual age carry the age signal in the final model."
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Paper covariate 'Gender', tested on volume, clearance, bioavailability (Table S3) and on the PD parameters (Table S4); not retained. Cohort 39.9% female (Table S5)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "ukat/L",
      type        = "continuous",
      notes       = "Tested on clearance (Table S3); not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "ukat/L",
      type        = "continuous",
      notes       = "Tested on clearance (Table S3); not retained."
    ),
    BILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Tested on clearance (Table S3); not retained."
    ),
    HGB = list(
      description = "Hemoglobin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Tested on volume (Table S3); not retained."
    ),
    HCT = list(
      description = "Hematocrit",
      units       = "(fraction)",
      type        = "continuous",
      notes       = "Tested on volume (Table S3); not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 208L,
    n_studies      = 3L,
    n_observations = 589L,
    age_range      = "0.011-17.9 years (postmenstrual age 38.6-970 weeks)",
    age_median     = "6.06 years (postmenstrual age 353 weeks)",
    weight_range   = "2.60-157 kg",
    weight_median  = "21.1 kg",
    sex_female_pct = 39.9,
    race_ethnicity = c(White = 60.6, Asian = 13.5, Black = 9.6, Other = 10.6, Unknown = 5.8),
    disease_state  = paste(
      "Pediatric patients 0 to <18 years who required or were on anticoagulant",
      "therapy: confirmed venous thromboembolism (Hokusai-VTE PEDIATRICS,",
      "N = 69) or cardiac disease at risk of thromboembolic events",
      "(ENNOBLE-ATE, N = 73), plus a single-dose phase 1 PK/PD study in",
      "children requiring anticoagulation (N = 66).",
      sep = " "
    ),
    dose_range     = paste(
      "Oral once daily, age- weight- and renal-function-banded (Table S2).",
      "Tablet for 12 to <18 years: 60 mg (>=60 kg), 45 mg (30 to <60 kg),",
      "30 mg (<30 kg). Oral suspension below 12 years: 60 mg (>=60 kg) or",
      "1.2 mg/kg capped at 45 mg (6 to <12 years), 1.4 mg/kg capped at 45 mg",
      "(2 to <6 years), 1.5 mg/kg capped at 45 mg (6 months to <2 years),",
      "0.8 mg/kg capped at 12 mg (0 to <6 months). Subjects with eGFR 30-50%",
      "of normal for age, or on a P-gp inhibitor other than amiodarone,",
      "received a 25-50% reduction. The phase 1 study gave single 30 mg- or",
      "60 mg-equivalent doses.",
      sep = " "
    ),
    renal_function = "eGFR (bedside Schwartz) median 111 mL/min/1.73 m^2, range 29.5-774 (Table S6). Renal function skews supranormal relative to the adult reference cohort.",
    regions        = "Multinational (NCT02303431 phase 1 PK/PD, NCT02798471 Hokusai-VTE PEDIATRICS, NCT03395639 ENNOBLE-ATE).",
    notes          = paste(
      "Demographics from Zou 2025 Tables S5 and S6. Of 605 plasma",
      "concentrations from 208 subjects, 7 (1.1%) were excluded as",
      "unrealistically high or low and 9 (1.5%) as below the 0.764 ng/mL",
      "LLOQ, leaving 589 observations. The PD datasets are subsets with",
      "time-matched concentrations: 233 anti-FXa observations from 122",
      "subjects, 431 aPTT observations from 197 subjects, and 432 PT",
      "observations from 198 subjects. 73.6% of subjects received the oral",
      "suspension and 26.4% tablets; 6.7% dosed fed. Only 7 subjects (3.4%)",
      "used a P-gp inhibitor and 1 (0.5%) a P-gp inducer, so those effects",
      "were not tested. Age bands: 0 to <0.5 y (N = 21), 0.5 to <2 y",
      "(N = 33), 2 to <6 y (N = 49), 6 to <12 y (N = 51), 12 to <18 y",
      "(N = 54).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # PK structural parameters. Zou 2025 Table 1, reported for a typical
    # pediatric subject with body weight 70 kg and eGFR 110 mL/min/1.73 m^2
    # (and, because the maturation factor is the bare Hill fraction, at
    # full renal maturity).
    # =====================================================================
    lcl  <- log(42.87); label("Apparent clearance CL/F for a 70 kg, eGFR 110 mL/min/1.73 m^2 subject (L/h)")                 # Table 1: CL/F = 42.87 L/h (RSE 3%)
    lvc  <- log(261);   label("Apparent central volume of distribution Vc/F for a 70 kg subject (L)")                        # Table 1: Vc/F = 261 L (RSE 0.9%)
    lq   <- log(8.59);  label("Apparent inter-compartmental clearance Q/F for a 70 kg subject (L/h)")                        # Table 1: Q/F = 8.59 L/h (RSE 2.3%)
    lvp  <- log(343.5); label("Apparent peripheral volume of distribution Vp/F for a 70 kg subject (L)")                     # Table 1: Vp/F = 343.5 L (RSE 10.7%)
    lka  <- log(3.71);  label("First-order absorption rate constant from the last transit compartment, Ka (1/h)")            # Table 1: Ka = 3.71 1/h (RSE 0.6%)
    lktr <- log(47.5);  label("Transit rate constant through the absorption chain, Ktr (1/h)")                               # Table 1: Ktr = 47.5 1/h (RSE 1%)

    # =====================================================================
    # Fixed allometric exponents. Supplement NONMEM $THETA block:
    #   0.75 FIX ; 9 ALLMCL BW      (applied to CL and Q)
    #   1 FIX    ; 10 ALLMV BW      (applied to Vc and Vp)
    # =====================================================================
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on CL/F and Q/F (unitless)")     # Table 1 'Fixed exponents for body weight-based scaling': 0.75 (CL/F and Q/F)
    e_wt_vc_vp <- fixed(1.00); label("Allometric exponent on Vc/F and Vp/F (unitless)")    # Table 1 'Fixed exponents for body weight-based scaling': 1.0

    # =====================================================================
    # Renal-function and maturation covariate effects on CL/F. Table 1
    # footnote a gives the complete covariate model:
    #   CL/F = CL/F_TYP * (WT/70)^0.75 * (eGFR/110)^0.268
    #                   * [PMA^3.4 / (47.7^3.4 + PMA^3.4)]
    # TM50 and the Hill coefficient are fixed to the Rhodin 2009 GFR
    # maturation values (Zou 2025 Methods 'PopPK model development':
    # "The values of time to half maturation (TM50) and Hill coefficient
    # were fixed as 47.7 weeks and 3.40, respectively.").
    # =====================================================================
    e_crcl_cl <- 0.268;      label("Exponent on (CRCL / 110) for CL/F (unitless)")                # Table 1: eGFR effect on clearance = 0.268 (RSE 16%)
    tmat50    <- fixed(47.7); label("Postmenstrual age at 50% renal maturation, TM50 (weeks; Rhodin 2009)")  # Table 1: TM50 = 47.7 (Fixed); Rhodin 2009
    hill_mat  <- fixed(3.40); label("Hill coefficient of the renal maturation function (unitless; Rhodin 2009)")  # Table 1: HILL = 3.40 (Fixed); Rhodin 2009

    # =====================================================================
    # PK inter-individual variability. Diagonal OMEGA (supplement NONMEM
    # $OMEGA has five diagonal entries and no block); all five parameters
    # enter log-normally (CL = TVCL*EXP(ETA(1)) etc. in $PK).
    #
    # Table 1 reports each IIV as a percent. That percent is the LOG-NORMAL
    # coefficient of variation, %CV = sqrt(exp(omega^2) - 1) * 100, not
    # sqrt(omega^2) * 100. The supplement's $OMEGA initial estimates -- which
    # are the previous run's near-final values -- settle this: the Vc entry
    # 0.119597 gives sqrt(exp(0.119597) - 1) = 35.64%, matching the reported
    # 35.6% exactly, where sqrt(0.119597) = 34.6% does not. So omega^2 is
    # recovered as log(1 + CV^2), written out below.
    # =====================================================================
    etalcl  ~ log(1 + 0.318^2)   # Table 1: IIV Clearance 31.8% (RSE 10%, shrinkage 15%)
    etalvc  ~ log(1 + 0.356^2)   # Table 1: IIV Central compartment volume 35.6% (RSE 10%, shrinkage 35%)
    etalq   ~ log(1 + 0.767^2)   # Table 1: IIV Inter-compartmental clearance 76.7% (RSE 7%, shrinkage 37%)
    etalktr ~ log(1 + 0.796^2)   # Table 1: IIV Transit rate constant 79.6% (RSE 9%, shrinkage 60%)
    etalka  ~ log(1 + 1.449^2)   # Table 1: IIV Absorption rate 144.9% (RSE 2%, shrinkage 49%)

    # =====================================================================
    # PK residual error. Supplement $ERROR:
    #   W = SQRT(THETA(1)**2 * IPRED**2 + THETA(2)**2); Y = IPRED + W*EPS(1)
    # i.e. combined proportional + additive standard deviations, with
    # $SIGMA 1 FIX. The additive term was fixed (Table 1 '0.71 FIX').
    # =====================================================================
    propSd <- 0.228;        label("Proportional residual error on edoxaban concentration (fraction)")  # Table 1: Proportional error 22.8% (RSE 8%)
    addSd  <- fixed(0.71);  label("Additive residual error on edoxaban concentration (ng/mL)")  # Table 1: Additive error 0.71 FIX

    # =====================================================================
    # PD -- anti-factor Xa activity (IU/mL). Emax model with IIV on Emax.
    # Zou 2025 Results 'PopPK/PD analyses': "Anti-FXa data were best fit
    # with an Emax model with IIV on Emax. Baseline anti-FXa (E0) was fixed
    # to 0.1 IU/mL because of the small number of pre-treatment anti-FXa
    # measurements (N = 9). An additive error model was used."
    # =====================================================================
    lrbase_antiFXa <- fixed(log(0.10)); label("Baseline anti-factor Xa activity, E0 (IU/mL)")           # Table 1 PD (Anti-FXa): Baseline = 0.10 (fixed)
    lemax_antiFXa  <- log(8.65);        label("Maximum edoxaban effect on anti-factor Xa activity, Emax (IU/mL)")  # Table 1 PD (Anti-FXa): Maximum effect = 8.65 (RSE 42.3%)
    lec50_antiFXa  <- log(631);         label("Edoxaban concentration at half-maximal anti-factor Xa effect, EC50 (ng/mL)")  # Table 1 PD (Anti-FXa): concentration at half max = 631 (RSE 50.0%)

    # Table 1 prints the PD IIV entries without a percent sign (14.8, 30.7,
    # 14.5) while the PK IIV entries carry one. They are percent CVs on the
    # same log-normal convention: read as variances they would imply
    # CVs of order 10^5, and read as log-scale SDs they would imply
    # omega = 14.8, both physically impossible for these parameters.
    etalemax_antiFXa ~ log(1 + 0.148^2)   # Table 1 PD (Anti-FXa): IIV on maximum effect 14.8% (RSE 28.9%, shrinkage 42.5%)

    addSd_antiFXa <- 0.247; label("Additive residual error on anti-factor Xa activity (IU/mL)")  # Table 1 PD (Anti-FXa): Additive error = 0.247 (RSE 15.9%)

    # =====================================================================
    # PD -- activated partial thromboplastin time (s). Linear model with
    # IIV on the intercept and proportional residual error (Zou 2025
    # Results 'PopPK/PD analyses').
    # =====================================================================
    lrbase_aPTT <- log(35.5);   label("Baseline activated partial thromboplastin time, E0 (s)")   # Table 1 PD (aPTT): Baseline = 35.5 s (RSE 2.30%)
    lslope_aPTT <- log(0.0467); label("Slope of the edoxaban effect on aPTT (s per ng/mL)")       # Table 1 PD (aPTT): slope = 0.0467 s mL/ng (RSE 10.3%)
    etalrbase_aPTT ~ log(1 + 0.307^2)   # Table 1 PD (aPTT): IIV on baseline aPTT 30.7% (RSE 13.7%, shrinkage 13.9%)
    propSd_aPTT <- 0.197; label("Proportional residual error on aPTT (fraction)")                 # Table 1 PD (aPTT): Proportional error 19.7% CV (RSE 15.3%)

    # =====================================================================
    # PD -- prothrombin time (s). Linear model with IIV on the intercept
    # and proportional residual error (Zou 2025 Results 'PopPK/PD
    # analyses').
    # =====================================================================
    lrbase_PT <- log(14.9);   label("Baseline prothrombin time, E0 (s)")                          # Table 1 PD (PT): Baseline = 14.9 s (RSE 1.47%)
    lslope_PT <- log(0.0415); label("Slope of the edoxaban effect on prothrombin time (s per ng/mL)")  # Table 1 PD (PT): slope = 0.0415 s mL/ng (RSE 3.92%)
    etalrbase_PT ~ log(1 + 0.145^2)     # Table 1 PD (PT): IIV on baseline PT 14.5% (RSE 20.3%, shrinkage 29.8%)
    propSd_PT <- 0.159; label("Proportional residual error on prothrombin time (fraction)")       # Table 1 PD (PT): Proportional error 15.9% CV (RSE 28.1%)
  })

  model({
    # --- 1. Derived covariate terms ---------------------------------------
    # Rhodin-style renal maturation on CL/F. Table 1 footnote a writes this
    # as the BARE Hill fraction PMA^3.4 / (47.7^3.4 + PMA^3.4), with no
    # normalisation to a reference postmenstrual age, and the supplement
    # NONMEM $PK block agrees:
    #   HILL = 3.4 ; TM50 = 47.7 ; FPMA = PMA**HILL/(TM50**HILL+PMA**HILL)
    # The factor therefore approaches 1 at full maturity and equals 0.5 at
    # PMA = 47.7 weeks, so lcl is the fully matured typical clearance.
    fmat <- PAGE^hill_mat / (tmat50^hill_mat + PAGE^hill_mat)

    # Renal-function power term, reference 110 mL/min/1.73 m^2.
    fegfr <- (CRCL / 110)^e_crcl_cl

    # --- 2. Individual PK parameters --------------------------------------
    # Allometric size on the 70 kg reference, per the supplement $PK block.
    cl  <- exp(lcl  + etalcl)  * (WT / 70)^e_wt_cl_q  * fegfr * fmat
    vc  <- exp(lvc  + etalvc)  * (WT / 70)^e_wt_vc_vp
    q   <- exp(lq   + etalq)   * (WT / 70)^e_wt_cl_q
    vp  <- exp(lvp)            * (WT / 70)^e_wt_vc_vp
    ka  <- exp(lka  + etalka)
    ktr <- exp(lktr + etalktr)

    # --- 3. Micro-constants ------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # --- 4. ODE system -----------------------------------------------------
    # Absorption chain exactly as the supplement $MODEL / $PK blocks build
    # it: COMP=(GUT) COMP=(CONC) COMP=(PERIPHERAL) COMP=(transit1) ...
    # COMP=(transit15), with K14 = K45 = ... = K17T18 = KTR carrying the
    # dose from GUT through the fifteen transit compartments, and
    # K18T2 = KA moving drug from the last transit compartment into the
    # central compartment. So the chain has fifteen Ktr-rate steps and one
    # terminal Ka step -- Ka is a separate estimated parameter, not Ktr.
    #
    # The paper's reported mean transit time of 0.337 h is the Savic
    # (NN + 1) / Ktr = 16 / 47.5 convention on NN = 15; the fifteen
    # Ktr-rate steps encoded here sum to 15 / 47.5 = 0.316 h, followed by
    # a mean 1 / Ka = 0.270 h absorption step. See the vignette Errata.
    d/dt(depot)     <- -ktr * depot
    d/dt(transit1)  <-  ktr * depot     - ktr * transit1
    d/dt(transit2)  <-  ktr * transit1  - ktr * transit2
    d/dt(transit3)  <-  ktr * transit2  - ktr * transit3
    d/dt(transit4)  <-  ktr * transit3  - ktr * transit4
    d/dt(transit5)  <-  ktr * transit4  - ktr * transit5
    d/dt(transit6)  <-  ktr * transit5  - ktr * transit6
    d/dt(transit7)  <-  ktr * transit6  - ktr * transit7
    d/dt(transit8)  <-  ktr * transit7  - ktr * transit8
    d/dt(transit9)  <-  ktr * transit8  - ktr * transit9
    d/dt(transit10) <-  ktr * transit9  - ktr * transit10
    d/dt(transit11) <-  ktr * transit10 - ktr * transit11
    d/dt(transit12) <-  ktr * transit11 - ktr * transit12
    d/dt(transit13) <-  ktr * transit12 - ktr * transit13
    d/dt(transit14) <-  ktr * transit13 - ktr * transit14
    d/dt(transit15) <-  ktr * transit14 - ka  * transit15

    d/dt(central)     <- ka * transit15 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                  k12 * central - k21 * peripheral1

    # --- 5. Observation and residual error ---------------------------------
    # Dose in mg and vc in L give mg/L; 1 mg/L = 1000 ng/mL, which matches
    # the supplement's scaling statement S2 = V/1000 (CONC = A2/S2).
    Cc <- 1000 * central / vc

    # Direct-response PD. Each layer is an algebraic function of the plasma
    # concentration -- the paper fitted them sequentially against observed
    # concentrations, so there is no effect compartment and no delay.
    rbase_antiFXa <- exp(lrbase_antiFXa)
    emax_antiFXa  <- exp(lemax_antiFXa + etalemax_antiFXa)
    ec50_antiFXa  <- exp(lec50_antiFXa)
    antiFXa <- rbase_antiFXa + emax_antiFXa * Cc / (ec50_antiFXa + Cc)

    rbase_aPTT <- exp(lrbase_aPTT + etalrbase_aPTT)
    slope_aPTT <- exp(lslope_aPTT)
    aPTT <- rbase_aPTT + slope_aPTT * Cc

    rbase_PT <- exp(lrbase_PT + etalrbase_PT)
    slope_PT <- exp(lslope_PT)
    PT <- rbase_PT + slope_PT * Cc

    Cc      ~ add(addSd) + prop(propSd)
    antiFXa ~ add(addSd_antiFXa)
    aPTT    ~ prop(propSd_aPTT)
    PT      ~ prop(propSd_PT)
  })
}
