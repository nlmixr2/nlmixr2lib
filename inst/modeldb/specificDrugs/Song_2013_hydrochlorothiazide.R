Song_2013_hydrochlorothiazide <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model with first-order absorption and an absorption lag time for",
    "hydrochlorothiazide in hypertensive adults and healthy volunteers, pooled across the CS-866 monotherapy and",
    "CS-8635 (olmesartan + amlodipine + hydrochlorothiazide) development programs. Apparent clearance rises with",
    "creatinine clearance and falls with age and in female subjects; both apparent volumes rise with body weight.",
    "The model exists to supply per-subject steady-state exposure AUCss = Dose / (CL/F) to the companion CS-8635",
    "blood-pressure exposure-response models, so its clearance term is the load-bearing part. Sister model files",
    "from the same paper: modellib('Song_2013_olmesartan'), modellib('Song_2013_amlodipine'),",
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
    "Structural parameters and inter-subject variability from Supplemental Table S5; the clearance covariate",
    "model from main-text Equation 10.",
    sep = " "
  )
  vignette <- "Song_2013_olmesartan_amlodipine_hydrochlorothiazide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "hydrochlorothiazide", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "hydrochlorothiazide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "hydrochlorothiazide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance calculated by the Cockcroft-Gault method.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "RAW Cockcroft-Gault creatinine clearance in mL/min -- NOT BSA-normalized (Song 2013, sentence following",
        "Equation 10). Time-fixed per subject. Enters apparent clearance in power form (CRCL / 117.5)^0.499.",
        "The centering value 117.5 mL/min is printed in Equation 10 and is the median of the hydrochlorothiazide",
        "population PK dataset; the N-weighted MEAN of that dataset computed from Supplemental Table S1 is",
        "120 mL/min. Note the centering value differs from the 111 mL/min used in the olmesartan model because",
        "the two drugs were fitted on different (overlapping) study sets."
      ),
      source_name        = "CLCR"
    ),
    AGE = list(
      description        = "Subject age.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters apparent clearance in power form (AGE / 49.5)^-0.214 (Song 2013",
        "Equation 10). The centering value 49.5 years is printed in the equation; the N-weighted MEAN age of the",
        "hydrochlorothiazide PK dataset computed from Supplemental Table S1 is 49.3 years. Song 2013 states this",
        "age effect is retained 'independently of any relationships between gender and age with renal function',",
        "i.e. it is not a proxy for the creatinine-clearance term that appears in the same equation."
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Time-fixed per subject. Enters apparent clearance EXPONENTIALLY as exp(-0.219 * SEXF), i.e. female",
        "subjects have 20% lower hydrochlorothiazide CL/F than males of the same age and renal function.",
        "The exponential form is taken from the printed Equation 10 ('20.3 x e^(-0.219 x Female)'), which is the",
        "authoritative form: Song 2013's generic categorical-covariate Equation 7 gives the fractional form",
        "(1 + coefficient x Covariate), and where a paper's generic template and its printed compound-specific",
        "equation disagree, the printed equation wins. The two forms differ by only 2% at this coefficient",
        "(exp(-0.219) = 0.803 vs 1 - 0.219 = 0.781). The source column is named 'Female' in Equation 10 and",
        "'CLSEX' in Supplemental Table S5; both carry the canonical SEXF orientation (1 = female), so no value",
        "transformation is required."
      ),
      source_name        = "Female"
    ),
    WT = list(
      description        = "Total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters both apparent volumes in power form, (WT / 90.7)^1.92 on Vc/F and",
        "(WT / 90.7)^0.846 on Vp/F (Song 2013 Supplemental Table S5 rows 'Vc,WTKG' and 'Vp,WTKG', with the",
        "median-normalized power form specified by main-text Equation 6). The exponent 1.92 on the central volume",
        "is far above the allometric 1.0 and is reported without an IIV estimate. Song 2013 never prints the",
        "weight centering value; 90.7 kg is the N-weighted MEAN body weight of the hydrochlorothiazide population",
        "PK dataset computed from Supplemental Table S1 (322 phase I subjects at 76.6 kg plus 956 CS8635-A-U301",
        "phase III subjects at 95.5 kg), used as a proxy for the unreported median. Because the exponent is",
        "large, the profile shape is more sensitive to this proxy than in the sibling models -- but CL/F, and",
        "hence the AUCss driving the exposure-response models, is untouched by it."
      ),
      source_name        = "WTKG"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1278L,
    n_studies      = 10L,
    age_range      = "phase I mean 29.9 years; CS8635-A-U301 phase III mean 55.8 (SD 10) years",
    weight_range   = "phase I mean 76.6 kg; CS8635-A-U301 phase III mean 95.5 (SD 22) kg",
    sex_female_pct = 40.3,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Healthy volunteers (phase I clinical pharmacology studies) plus adults with hypertension",
      "(phase III sparse-sampling subset of CS8635-A-U301 (TRINITY))"
    ),
    dose_range     = "hydrochlorothiazide 12.5-25 mg once daily",
    regions        = "United States and Europe",
    notes          = paste(
      "The hydrochlorothiazide population PK dataset differs from its olmesartan and amlodipine siblings: it is",
      "the union of the three CS-866 hydrochlorothiazide phase I studies (CS866-126 n = 30, CS866-127 n = 18,",
      "CS866-134 n = 29), the six CS-8635 phase I studies (n = 245), and the CS8635-A-U301 phase III PK subset",
      "(n = 956), giving n = 1278. The CS8663-A-U301 (COACH) study contributed no hydrochlorothiazide data",
      "because that program studied olmesartan + amlodipine only. Estimation used FOCE in NONMEM VI level 1."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters: Song 2013 Supplemental Table S5, 'Population
    # Mean / Estimate' column. Table S5 labels the volumes with the NONMEM
    # names V2 and V3, which map to the canonical Vc and Vp respectively
    # (confirmed by the same table's 'Vc,WTKG' / 'Vp,WTKG' covariate rows).
    # ------------------------------------------------------------------
    lcl   <- log(20.3)  ; label("Apparent clearance CL/F for a male at CRCL 117.5 mL/min and 49.5 years (L/h)")  # Table S5 CLTYP = 20.3 (SE 1.4% CV); also main-text Eq. 10 intercept
    lvc   <- log(27.7)  ; label("Apparent central volume Vc/F at 90.7 kg (L)")                                   # Table S5 V2 = 27.7 (SE 5.5% CV)
    lvp   <- log(174)   ; label("Apparent peripheral volume Vp/F at 90.7 kg (L)")                                # Table S5 V3 = 174 (SE 1.7% CV)
    lq    <- log(18.3)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")                              # Table S5 Q = 18.3 (SE 2.0% CV)
    lka   <- log(0.364) ; label("First-order absorption rate constant (1/h)")                                    # Table S5 Ka = 0.364 (SE 2.0% CV)
    ltlag <- log(0.419) ; label("Absorption lag time (h)")                                                       # Table S5 ALAG1 = 0.419 (SE 1.6% CV)

    # ------------------------------------------------------------------
    # Covariate effects. Continuous covariates enter as power functions of
    # the median-normalized value (Song 2013 Equation 6); the sex effect is
    # exponential per the printed Equation 10.
    # ------------------------------------------------------------------
    e_crcl_cl <-  0.499 ; label("Power exponent for creatinine clearance on CL/F (unitless)")                  # Table S5 CLCLCR = 0.499 (SE 9.5% CV); main-text Eq. 10 exponent
    e_age_cl  <- -0.214 ; label("Power exponent for age on CL/F (unitless)")                                   # Table S5 CLAGE = -0.214 (SE 13% CV); main-text Eq. 10 exponent
    e_sexf_cl <- -0.219 ; label("Log-scale female effect on CL/F, applied as exp(e_sexf_cl * SEXF) (unitless)") # Table S5 CLSEX = -0.219 (SE 11% CV); main-text Eq. 10 e^(-0.219 x Female)
    e_wt_vc   <-  1.92  ; label("Power exponent for body weight on Vc/F (unitless)")                           # Table S5 Vc,WTKG = 1.92 (SE 8.8% CV)
    e_wt_vp   <-  0.846 ; label("Power exponent for body weight on Vp/F (unitless)")                           # Table S5 Vp,WTKG = 0.846 (SE 10% CV)

    # ------------------------------------------------------------------
    # IIV. Table S5 reports inter-subject variability as an "approximate
    # percent coefficient of variation" (footnote b), so the log-normal
    # convention omega^2 = log(1 + CV^2) is used. Table S5 reports no IIV
    # on the absorption lag time.
    # ------------------------------------------------------------------
    etalcl ~ 0.0861777  # Table S5 CLTYP IIV 30% CV -> log(1 + 0.30^2)
    etalvc ~ 0.5734623  # Table S5 V2    IIV 88% CV -> log(1 + 0.88^2)
    etalvp ~ 0.0560022  # Table S5 V3    IIV 24% CV -> log(1 + 0.24^2)
    etalq  ~ 0.0472652  # Table S5 Q     IIV 22% CV -> log(1 + 0.22^2)
    etalka ~ 0.0318862  # Table S5 Ka    IIV 18% CV -> log(1 + 0.18^2)

    # ------------------------------------------------------------------
    # Residual error. Table S5 reports two sigmas as VARIANCES, both
    # MULTIPLICATIVE and split by study phase: sigma^2_1 (phase I) = 0.0595
    # and sigma^2_2 (phase III) = 0.0819. There is no additive term at all.
    # nlmixr2's error-model syntax cannot switch a residual SD on a
    # covariate, so the phase I (full-profile) value is encoded and the
    # larger phase III sparse-sampling value (28.6%) is documented in the
    # vignette deviations section.
    # ------------------------------------------------------------------
    propSd <- 0.243926 ; label("Proportional residual error, phase I full-profile sampling (fraction)")  # Table S5 sigma^2_1 = 0.0595 -> sqrt = 0.2439
  })

  model({
    # ---- Covariate centering values -----------------------------------
    # ref_crcl and ref_age are printed in Song 2013 Equation 10. ref_wt is
    # NOT printed anywhere in the paper; it is the N-weighted mean body
    # weight of the hydrochlorothiazide PK dataset derived from Supplemental
    # Table S1, used as a proxy for the unreported median (see
    # covariateData[["WT"]]$notes).
    ref_crcl <- 117.5
    ref_age  <- 49.5
    ref_wt   <- 90.7

    # ---- Individual parameters ----------------------------------------
    # Song 2013 Eq. 10:
    #   (CL/F)_i = 20.3 * e^(-0.219 * Female) * (CLCR_i / 117.5)^0.499
    #                                          * (Age_i / 49.5)^-0.214
    cl <- exp(lcl + etalcl + e_sexf_cl * SEXF) *
      (CRCL / ref_crcl)^e_crcl_cl *
      (AGE / ref_age)^e_age_cl
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
    Cc ~ prop(propSd)
  })
}
