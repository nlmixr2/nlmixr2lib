Song_2013_amlodipine <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model with first-order absorption and an absorption lag time for",
    "amlodipine in hypertensive adults and healthy volunteers, pooled across the CS-8663 (olmesartan + amlodipine)",
    "and CS-8635 (olmesartan + amlodipine + hydrochlorothiazide) development programs. Apparent clearance falls",
    "with age; the apparent central volume rises with body weight. The model exists to supply per-subject",
    "steady-state exposure AUCss = Dose / (CL/F) to the companion CS-8635 blood-pressure exposure-response models,",
    "so its clearance term is the load-bearing part. Sister model files from the same paper:",
    "modellib('Song_2013_olmesartan'), modellib('Song_2013_hydrochlorothiazide'),",
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
    "Structural parameters and inter-subject variability from Supplemental Table S4; the clearance covariate",
    "model from main-text Equation 9.",
    sep = " "
  )
  vignette <- "Song_2013_olmesartan_amlodipine_hydrochlorothiazide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "amlodipine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "amlodipine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "amlodipine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description        = "Subject age.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters apparent clearance in power form",
        "(CL/F) = 23.4 * (AGE / 50.9)^-0.349 (Song 2013 Equation 9), i.e. amlodipine clearance falls with age.",
        "The centering value 50.9 years is printed in the equation itself and is the median of the amlodipine",
        "population PK dataset; the N-weighted MEAN age of that dataset computed from Supplemental Table S1 is",
        "50.2 years, which is the calibration used to justify the mean-for-median substitution applied to body",
        "weight (see covariateData[['WT']]$notes). The dataset is strongly bimodal in age -- young healthy",
        "phase I volunteers (mean 31 years) pooled with older hypertensive phase III patients (mean 55 years)."
      ),
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters the apparent central volume in power form (WT / 91.5)^0.285 (Song 2013",
        "Supplemental Table S4 row 'Vc,WTKG', with the median-normalized power form specified by main-text",
        "Equation 6). Table S4 reports NO weight effect on the peripheral volume, in contrast to the olmesartan",
        "and hydrochlorothiazide models. Song 2013 never prints the weight centering value; 91.5 kg is the",
        "N-weighted MEAN body weight of the amlodipine population PK dataset computed from Supplemental Table S1",
        "(the same dataset as olmesartan), used as a proxy for the unreported median."
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
    dose_range     = "amlodipine 5-10 mg once daily",
    regions        = "United States and Europe",
    notes          = paste(
      "The amlodipine population PK dataset is the same union of CS-8663 and CS-8635 program studies as the",
      "olmesartan model (Song 2013 Supplemental Table S1), n = 1927: 415 phase I subjects plus the phase III PK",
      "subsets CS8663-A-U301 (n = 556) and CS8635-A-U301 (n = 956). Estimation used FOCE in NONMEM VI level 1."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters: Song 2013 Supplemental Table S4, 'Population
    # Mean / Estimate' column. All are apparent (oral) parameters.
    # ------------------------------------------------------------------
    lcl   <- log(23.4)  ; label("Apparent clearance CL/F at 50.9 years (L/h)")                 # Table S4 CLTYP = 23.4 (SE 1.1% CV); also main-text Eq. 9 intercept
    lvc   <- log(1060)  ; label("Apparent central volume Vc/F at 91.5 kg (L)")                 # Table S4 Vc = 1060 (SE 1.7% CV)
    lvp   <- log(465)   ; label("Apparent peripheral volume Vp/F (L)")                         # Table S4 Vp = 465 (SE 2.8% CV)
    lq    <- log(26.6)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")            # Table S4 Q = 26.6 (SE 4.1% CV)
    lka   <- log(0.215) ; label("First-order absorption rate constant (1/h)")                  # Table S4 Ka = 0.215 (SE 2.6% CV)
    ltlag <- log(0.315) ; label("Absorption lag time (h)")                                     # Table S4 ALAG1 = 0.315 (SE 1.7% CV)

    # ------------------------------------------------------------------
    # Covariate effects, power functions of the median-normalized value
    # (Song 2013 Equation 6).
    # ------------------------------------------------------------------
    e_age_cl <- -0.349 ; label("Power exponent for age on CL/F (unitless)")                    # Table S4 CLAGE = -0.349 (SE 7.2% CV); main-text Eq. 9 exponent
    e_wt_vc  <-  0.285 ; label("Power exponent for body weight on Vc/F (unitless)")            # Table S4 Vc,WTKG = 0.285 (SE 18% CV)

    # ------------------------------------------------------------------
    # IIV. Table S4 reports inter-subject variability as an "approximate
    # percent coefficient of variation" (footnote b), so the log-normal
    # convention omega^2 = log(1 + CV^2) is used. Table S4 leaves the IIV
    # cells for Q and ALAG1 blank -- no random effect on either.
    # ------------------------------------------------------------------
    etalcl ~ 0.1415864  # Table S4 CLTYP IIV 39% CV -> log(1 + 0.39^2)
    etalvc ~ 0.0861777  # Table S4 Vc    IIV 30% CV -> log(1 + 0.30^2)
    etalvp ~ 0.0284903  # Table S4 Vp    IIV 17% CV -> log(1 + 0.17^2)
    etalka ~ 0.1348805  # Table S4 Ka    IIV 38% CV -> log(1 + 0.38^2)

    # ------------------------------------------------------------------
    # Residual error. Table S4 reports only two sigmas as VARIANCES:
    # sigma^2_1 multiplicative = 0.045 and a study-specific sigma^2_3
    # additive for CS8635-A-U301 = 1.6 (ng/mL)^2. There is no generally
    # applicable additive term, so only the proportional error is encoded;
    # the CS8635-A-U301-specific additive term (1.26 ng/mL) is documented
    # in the vignette deviations section.
    # ------------------------------------------------------------------
    propSd <- 0.212132 ; label("Proportional residual error (fraction)")   # Table S4 sigma^2_1 = 0.045 -> sqrt = 0.2121
  })

  model({
    # ---- Covariate centering values -----------------------------------
    # ref_age is printed in Song 2013 Equation 9. ref_wt is NOT printed
    # anywhere in the paper; it is the N-weighted mean body weight of the
    # amlodipine PK dataset derived from Supplemental Table S1 and used as
    # a proxy for the unreported median (see covariateData[["WT"]]$notes).
    ref_age <- 50.9
    ref_wt  <- 91.5

    # ---- Individual parameters ----------------------------------------
    # Song 2013 Eq. 9: (CL/F)_i = 23.4 * (Age_i / 50.9)^-0.349.
    cl <- exp(lcl + etalcl) * (AGE / ref_age)^e_age_cl
    vc <- exp(lvc + etalvc) * (WT / ref_wt)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)
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
    # (validated 0.05-50 ng/mL, LLOQ 0.05 ng/mL).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
