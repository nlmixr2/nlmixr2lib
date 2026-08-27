Song_2013_olmesartan <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model with first-order absorption and an absorption lag time for",
    "olmesartan (the active moiety of olmesartan medoxomil) in hypertensive adults and healthy volunteers, pooled",
    "across the CS-8663 (olmesartan + amlodipine) and CS-8635 (olmesartan + amlodipine + hydrochlorothiazide)",
    "development programs. Apparent clearance rises with creatinine clearance; apparent central and peripheral",
    "volumes rise with body weight. The model exists to supply per-subject steady-state exposure",
    "AUCss = Dose / (CL/F) to the companion CS-8635 blood-pressure exposure-response models, so its clearance",
    "term is the load-bearing part. Sister model files from the same paper: modellib('Song_2013_amlodipine'),",
    "modellib('Song_2013_hydrochlorothiazide'),",
    "modellib('Song_2013_olmesartan_amlodipine_hydrochlorothiazide_dbp'),",
    "modellib('Song_2013_olmesartan_amlodipine_hydrochlorothiazide_sbp')."
  )
  reference <- paste(
    "Song S, Carrothers TJ, Moore H, Green M, Miller R, Rohatagi S, Lee J, Wang A, Melino M, Patel M,",
    "Heyrman R, Salazar DE.",
    "Model-supported development of CS-8635: a fixed-dose combination of olmesartan, amlodipine, and",
    "hydrochlorothiazide.",
    "Clin Pharmacol Drug Dev. 2013;2(2):103-112.",
    "doi:10.1002/cpdd.17.",
    "Structural parameters and inter-subject variability from Supplemental Table S3; the clearance covariate",
    "model from main-text Equation 8.",
    sep = " "
  )
  vignette <- "Song_2013_olmesartan_amlodipine_hydrochlorothiazide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "olmesartan", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "olmesartan", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "olmesartan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance calculated by the Cockcroft-Gault method.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "RAW Cockcroft-Gault creatinine clearance in mL/min -- NOT BSA-normalized. Song 2013 states the method",
        "explicitly in the sentence following Equation 10 ('CLCR represent ... creatinine clearance calculated by",
        "Cockroft and Gault method'). Time-fixed per subject. Enters apparent clearance in power form",
        "(CL/F) = 6.32 * (CRCL / 111)^0.425 (Song 2013 Equation 8). The centering value 111 mL/min is printed in",
        "the equation itself and is the median of the olmesartan population PK dataset; for reference, the",
        "N-weighted MEAN CLCR of that dataset computed from Supplemental Table S1 is 115 mL/min, so 111 is",
        "consistent with a median of a right-skewed distribution."
      ),
      source_name        = "CLCR"
    ),
    WT = list(
      description        = "Total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters both apparent volumes in power form, (WT / 91.5)^0.681 on Vc/F and",
        "(WT / 91.5)^0.405 on Vp/F (Song 2013 Supplemental Table S3 rows 'Vc,WTKG' and 'Vp,WTKG', with the",
        "median-normalized power form specified by main-text Equation 6). IMPORTANT: unlike the creatinine-",
        "clearance centering value, Song 2013 never prints the weight centering value. The 91.5 kg used here is",
        "the N-weighted MEAN body weight of the olmesartan population PK dataset computed from Supplemental",
        "Table S1 (415 phase I subjects at 78.0 kg plus 1512 phase III subjects at 95.2 kg), used as a proxy for",
        "the unreported median. The proxy is calibrated: applying the same mean-for-median substitution to the",
        "covariates whose centering values ARE printed recovers them to within 1-3% (age 50.2 vs 50.9 printed for",
        "amlodipine, age 49.3 vs 49.5 for hydrochlorothiazide). This affects the shape of the concentration-time",
        "profile only -- it does not touch CL/F and therefore does not touch the AUCss that drives the companion",
        "exposure-response models."
      ),
      source_name        = "WTKG"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1927L,
    n_studies      = 12L,
    age_range      = "phase I mean 30.9 (SD 7.9) years; phase III mean 55.4 (SD 11) years",
    weight_range   = "phase I mean 77 (SD 13) kg; phase III mean 95.2 (SD 22) kg",
    sex_female_pct = 37.0,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Healthy volunteers (phase I clinical pharmacology studies) plus adults with hypertension",
      "(phase III sparse-sampling subsets of CS8663-A-U301 (COACH) and CS8635-A-U301 (TRINITY))"
    ),
    dose_range     = "olmesartan medoxomil 20-40 mg once daily",
    regions        = "United States and Europe",
    notes          = paste(
      "The olmesartan population PK dataset is the union of the CS-8663 and CS-8635 program studies listed in",
      "Song 2013 Supplemental Table S1 (phase I: CS8663-A-U101/U110/U111/U112 and CS8635-A-U101/U102/U103/U104/",
      "A-E105/U106, n = 415; phase III PK subsets: CS8663-A-U301 n = 556 and CS8635-A-U301 n = 956), giving",
      "n = 1927. Mean creatinine clearance was 127 mL/min in phase I and 111 mL/min in phase III;",
      "13.7-16.3% of phase III subjects were diabetic. The sex split shown here is the N-weighted female",
      "fraction across those same studies. Estimation used FOCE in NONMEM VI level 1."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters: Song 2013 Supplemental Table S3, 'Population
    # Mean / Estimate' column. All are apparent (oral) parameters -- the
    # paper writes CL/F and does not report absolute bioavailability.
    # ------------------------------------------------------------------
    lcl   <- log(6.32)  ; label("Apparent clearance CL/F at CRCL 111 mL/min (L/h)")            # Table S3 CLTYP = 6.32 (SE 1.1% CV); also main-text Eq. 8 intercept
    lvc   <- log(36.8)  ; label("Apparent central volume Vc/F at 91.5 kg (L)")                 # Table S3 Vc = 36.8 (SE 1.4% CV)
    lvp   <- log(29.0)  ; label("Apparent peripheral volume Vp/F at 91.5 kg (L)")              # Table S3 Vp = 29.0 (SE 3.0% CV)
    lq    <- log(1.64)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")            # Table S3 Q = 1.64 (SE 2.3% CV)
    lka   <- log(1.25)  ; label("First-order absorption rate constant (1/h)")                  # Table S3 Ka = 1.25 (SE 3.3% CV)
    ltlag <- log(0.406) ; label("Absorption lag time (h)")                                     # Table S3 ALAG1 = 0.406 (SE 0.5% CV)

    # ------------------------------------------------------------------
    # Covariate effects. Continuous covariates enter as power functions of
    # the median-normalized value (Song 2013 Equation 6).
    # ------------------------------------------------------------------
    e_crcl_cl <- 0.425 ; label("Power exponent for creatinine clearance on CL/F (unitless)")   # Table S3 CLCLCR = 0.425 (SE 7.4% CV); main-text Eq. 8 exponent
    e_wt_vc   <- 0.681 ; label("Power exponent for body weight on Vc/F (unitless)")            # Table S3 Vc,WTKG = 0.681 (SE 8.3% CV)
    e_wt_vp   <- 0.405 ; label("Power exponent for body weight on Vp/F (unitless)")            # Table S3 Vp,WTKG = 0.405 (SE 34% CV)

    # ------------------------------------------------------------------
    # IIV. Table S3 reports inter-subject variability as an "approximate
    # percent coefficient of variation" (footnote b), so the log-normal
    # convention omega^2 = log(1 + CV^2) is used. Table S3 reports no IIV
    # on the absorption lag time.
    # ------------------------------------------------------------------
    etalcl ~ 0.1415864  # Table S3 CLTYP IIV 39% CV -> log(1 + 0.39^2)
    etalvc ~ 0.0654131  # Table S3 Vc    IIV 26% CV -> log(1 + 0.26^2)
    etalvp ~ 0.2311911  # Table S3 Vp    IIV 51% CV -> log(1 + 0.51^2)
    etalq  ~ 0.1415864  # Table S3 Q     IIV 39% CV -> log(1 + 0.39^2)
    etalka ~ 0.7830336  # Table S3 Ka    IIV 109% CV -> log(1 + 1.09^2)

    # ------------------------------------------------------------------
    # Residual error. Table S3 reports three sigmas as VARIANCES: sigma^2_1
    # multiplicative = 0.0839, sigma^2_2 additive = 0.0627 (ng/mL)^2, and a
    # study-specific sigma^2_3 additive for CS8635-A-U301 = 298 (ng/mL)^2.
    # nlmixr2's error-model syntax cannot switch an additive SD on a
    # covariate, so the two generally-applicable terms are encoded and the
    # CS8635-A-U301-specific additive term (17.3 ng/mL) is documented in the
    # vignette deviations section instead.
    # ------------------------------------------------------------------
    propSd <- 0.289655 ; label("Proportional residual error (fraction)")   # Table S3 sigma^2_1 = 0.0839 -> sqrt = 0.2897
    addSd  <- 0.250400 ; label("Additive residual error (ng/mL)")          # Table S3 sigma^2_2 = 0.0627 (ng/mL)^2 -> sqrt = 0.2504
  })

  model({
    # ---- Covariate centering values -----------------------------------
    # ref_crcl is printed in Song 2013 Equation 8. ref_wt is NOT printed
    # anywhere in the paper; it is the N-weighted mean body weight of the
    # olmesartan PK dataset derived from Supplemental Table S1 and used as
    # a proxy for the unreported median (see covariateData[["WT"]]$notes).
    ref_crcl <- 111
    ref_wt   <- 91.5

    # ---- Individual parameters ----------------------------------------
    # Song 2013 Eq. 8: (CL/F)_i = 6.32 * (CLCR_i / 111)^0.425.
    cl <- exp(lcl + etalcl) * (CRCL / ref_crcl)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / ref_wt)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / ref_wt)^e_wt_vp
    q  <- exp(lq + etalq)
    ka <- exp(lka + etalka)

    # ---- Micro-constants ----------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Two-compartment ODE system with first-order absorption --------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # NM-TRAN ALAG1: absorption does not begin until ltlag hours after the dose.
    alag(depot) <- exp(ltlag)

    # ---- Observation ---------------------------------------------------
    # Doses are in mg and volumes in L, so central / vc is in mg/L; the
    # factor 1000 converts to the ng/mL scale of the bioanalytical assay
    # (validated 1-1000 ng/mL, LLOQ 1 ng/mL).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
