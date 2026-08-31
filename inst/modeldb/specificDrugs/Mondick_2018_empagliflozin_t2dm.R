Mondick_2018_empagliflozin_t2dm <- function() {
  description <- paste(
    "Joint population PK-PD model of empagliflozin in adults with type 2",
    "diabetes mellitus, in the form used for the type 1 vs type 2 comparison of",
    "Mondick 2018. The mechanistic renal glucose-reabsorption PD parameters are",
    "the type 2 diabetes estimates of Mondick 2016 as reproduced in Mondick",
    "2018 Supplementary Table 2; the two-compartment oral PK layer is the type 1",
    "diabetes PK model of Mondick 2018, which that paper deliberately reuses for",
    "both populations so the comparison isolates PD differences.",
    sep = " "
  )
  reference <- paste(
    "PD parameters: Mondick J, Riggs M, Sasaki T, Sarashina A, Broedl UC,",
    "Retlich S. Mixed-effects modelling to quantify the effect of empagliflozin",
    "on renal glucose reabsorption in patients with type 2 diabetes.",
    "Diabetes Obes Metab. 2016;18:241-248, as tabulated in Mondick J, Riggs M,",
    "Kaspers S, Soleymanlou N, Marquard J, Nock V. J Clin Pharmacol.",
    "2018;58(5):640-649 Supplementary Table 2, doi:10.1002/jcph.1051.",
    "PK parameters: Mondick 2018 Supplementary Table 1; see",
    "modellib('Mondick_2018_empagliflozin_t1dm').",
    sep = " "
  )
  vignette <- "Mondick_2018_empagliflozin"
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")

  # `glu_urine` follows the naming already established for cumulative excreted
  # urinary glucose in Lu_2014_sglt_qsp.R (same analyte, same role).
  paper_specific_compartments <- c("glu_urine")

  covariateData <- list(
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-scaled against a reference of 25 kg/m^2 on CL/F, Vc/F, Q/F and",
        "Vp/F. The PK layer, including this covariate, is the type 1 diabetes",
        "PK model of Mondick 2018 Supplementary Table 1; the cohort statistic",
        "quoted there (mean 25.6 kg/m^2) is the T1DM cohort's, not the T2DM",
        "cohort's."
      ),
      source_name        = "BMI"
    ),
    CRCL = list(
      description        = paste(
        "BSA-normalized estimated glomerular filtration rate (eGFR). The",
        "estimating equation is not named in Mondick 2018."
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Mechanistic driver rather than a covariate on a PK parameter: the",
        "model de-normalizes it with the individual BSA and multiplies it by",
        "the plasma-minus-reabsorbed glucose concentration difference to give",
        "the urinary glucose excretion rate. Mondick 2018 does not report the",
        "T2DM cohort's eGFR distribution."
      ),
      source_name        = "eGFR"
    ),
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used solely to convert the BSA-normalized CRCL back to an absolute",
        "filtration rate, as CRCL * BSA / 1.73."
      ),
      source_name        = "BSA"
    ),
    GLU = list(
      description        = "Plasma glucose concentration, time-varying regressor",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Exogenous driving regressor, linearly interpolated between observed",
        "time points via linear(GLU). Converted internally from mmol/L to the",
        "paper's mg/dL (factor 18.016). For its type 1 vs type 2 comparison",
        "Mondick 2018 fed BOTH populations the median plasma-glucose profile",
        "observed in the T1DM study, so that the comparison could not be driven",
        "by differences in glucose sampling."
      ),
      source_name        = "PG"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-scaled against the T2DM model's 58-year reference patient on",
        "both Gmax and Imax (Mondick 2018 Supplementary Table 2 footnote:",
        "'Reference age was 42 years for T1DM model and 58 years for T2DM",
        "model'). Note the reference differs from the T1DM sibling model."
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "categorical",
      reference_category = "0 = male (the T2DM reference patient is a 58-year-old male)",
      notes              = paste(
        "Enters as a multiplicative factor raised to the SEXF indicator on both",
        "Gmax and Imax (Mondick 2018 Supplementary Table 2, T2DM column)."
      ),
      source_name        = "SEX"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "empagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "empagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "empagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    glu_urine   = list(analyte = "glucose", units = "mg", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 4,
    age_median     = "reference patient 58 years",
    disease_state  = "type 2 diabetes mellitus",
    dose_range     = "empagliflozin 1 mg to 100 mg",
    notes          = paste(
      "The PD parameter set was estimated in Mondick 2016 (Diabetes Obes Metab",
      "18:241-248) from 4 phase 1 and phase 2 empagliflozin trials spanning",
      "1-100 mg, with 1 further study used for model evaluation. Mondick 2018",
      "reproduces the complete parameter set in its Supplementary Table 2 and",
      "reuses it for the type 1 vs type 2 comparison in its Figures 4 and 5,",
      "but does not restate the T2DM cohort size or its baseline demographics;",
      "n_subjects is therefore NA here. Mondick 2016 is not on disk for this",
      "extraction, so every value in this file is traced to Mondick 2018",
      "Supplementary Table 2 and should be re-verified against Mondick 2016 if",
      "that paper is acquired."
    )
  )

  ini({
    # ---- Population PK -----------------------------------------------------
    # Identical to Mondick_2018_empagliflozin_t1dm: Mondick 2018 states that
    # for the T1DM-vs-T2DM comparison "the PK parameter estimates from the
    # analysis presented here were used to simulate PK profiles for both
    # populations", so this is the PK layer the authors themselves paired with
    # the T2DM PD parameters. All values from Suppl Table 1.
    lcl   <- log(12.3)  ; label("Apparent oral clearance CL/F at BMI 25 kg/m^2 (L/h)")                # Suppl Table 1, 12.3 L/h (%RSE 1.14)
    lvc   <- log(3.47)  ; label("Apparent central volume Vc/F at BMI 25 kg/m^2 (L)")                  # Suppl Table 1, 3.47 L (%RSE 7.10)
    lq    <- log(7.95)  ; label("Apparent intercompartmental clearance Q/F at BMI 25 kg/m^2 (L/h)")   # Suppl Table 1, 7.95 L/h (%RSE 3.23)
    lvp   <- log(88.0)  ; label("Apparent peripheral volume Vp/F at BMI 25 kg/m^2 (L)")               # Suppl Table 1, 88.0 L (%RSE 0.952)
    lka   <- log(0.275) ; label("First-order absorption rate constant ka (1/h)")                      # Suppl Table 1, 0.275 1/h (%RSE 3.17)
    ld1   <- log(0.607) ; label("Duration of the zero-order input into the depot D1 (h)")             # Suppl Table 1, 0.607 h (%RSE 15.1)
    ltlag <- log(0.189) ; label("Absorption lag time ALAG1 (h)")                                      # Suppl Table 1, 0.189 h (%RSE 5.59)

    e_bmi_cl <-  0.554 ; label("Power exponent on (BMI/25) for CL/F (unitless)")  # Suppl Table 1, BMI effect on CL/F 0.554 (%RSE 32.6)
    e_bmi_vc <- -0.241 ; label("Power exponent on (BMI/25) for Vc/F (unitless)")  # Suppl Table 1, BMI effect on Vc/F -0.241 (%RSE 223)
    e_bmi_q  <-  1.77  ; label("Power exponent on (BMI/25) for Q/F (unitless)")   # Suppl Table 1, BMI effect on Q/F 1.77 (%RSE 16.4)
    e_bmi_vp <-  0.775 ; label("Power exponent on (BMI/25) for Vp/F (unitless)")  # Suppl Table 1, BMI effect on Vp/F 0.775 (%RSE 49.9)

    # ---- Population PD: Suppl Table 2, "Patients with type 2 diabetes" ------
    lgmax <- log(330)  ; label("Maximum reabsorbed glucose concentration Gmax (mg/dL)")                            # Suppl Table 2 T2DM, 330 mg/dL (95% CI 311, 350)
    lkm   <- log(105)  ; label("Plasma glucose concentration at half-maximal reabsorption KM (mg/dL)")             # Suppl Table 2 T2DM, 105 mg/dL (95% CI 82.1, 123)
    lic50 <- log(4.61) ; label("Empagliflozin concentration at half-maximal reabsorption inhibition IC50 (nmol/L)")  # Suppl Table 2 T2DM, 4.61 nmol/L (95% CI 2.81, 7.20)

    # Logit-scale, as in the T1DM sibling model; see that file for the
    # (1 - theta) vs (1 + theta) denominator note.
    logitimax <- log(0.573 / (1 - 0.573)) ; label("Logit of maximum inhibition of renal glucose reabsorption Imax (unitless)")  # Suppl Table 2 T2DM, Imax 0.573 (95% CI 0.554, 0.607)
    logitfrac <- log(0.999 / (1 - 0.999)) ; label("Logit of the fraction of filtered glucose reabsorbed below RTG, FRAC (unitless)")  # Suppl Table 2 T2DM, FRAC 0.999 (95% CI 0.998, 0.999)

    # Full covariate model on Gmax and Imax, referenced to a 58-year-old male.
    e_sexf_gmax <-  1.12    ; label("Multiplicative factor on Gmax for female sex (unitless)")  # Suppl Table 2 T2DM, 1.12 (95% CI 1.06, 1.20)
    e_age_gmax  <- -0.0554  ; label("Power exponent on (AGE/58) for Gmax (unitless)")           # Suppl Table 2 T2DM, -0.0554 (95% CI -0.125, 0.0699)
    e_sexf_imax <-  0.984   ; label("Multiplicative factor on Imax for female sex (unitless)")  # Suppl Table 2 T2DM, 0.984 (95% CI 0.898, 1.05)
    e_age_imax  <- -0.170   ; label("Power exponent on (AGE/58) for Imax (unitless)")           # Suppl Table 2 T2DM, -0.170 (95% CI -0.284, 0.00908)

    # ---- Interindividual variability ---------------------------------------
    # PK IIV from Suppl Table 1 (the T1DM PK layer, reused here per the paper).
    etalcl   ~ 0.0387  # Suppl Table 1, IIV CL/F   0.0387 (%CV 19.9)
    etalq    ~ 0.0422  # Suppl Table 1, IIV Q/F    0.0422 (%CV 20.8)
    etalvp   ~ 0.0772  # Suppl Table 1, IIV Vp/F   0.0772 (%CV 28.3)
    etalka   ~ 0.0286  # Suppl Table 1, IIV Ka     0.0286 (%CV 17.0)
    etald1   ~ 0.174   # Suppl Table 1, IIV D1     0.174  (%CV 43.6)
    etaltlag ~ 0.309   # Suppl Table 1, IIV ALAG1  0.309  (%CV 60.2)

    # PD IIV block, T2DM column. Correlation check:
    # -0.0131 / sqrt(0.00895 * 1.35) = -0.119, the tabulated CORR.
    etalgmax + etalogitfrac ~ c(0.00895,
                                -0.0131, 1.35)  # Suppl Table 2 T2DM, IIV Gmax 0.00895; cov -0.0131; IIV FRAC 1.35

    # ---- Residual error ----------------------------------------------------
    # sqrt(0.0584) = 0.2417 (24.2% CV); sqrt(0.207) = 0.4550 (45.5% CV).
    propSd     <- 0.24166 ; label("Proportional residual error on empagliflozin concentration (fraction)")  # Suppl Table 1, err_prop 0.0584 (%CV 24.2)
    propSd_UGE <- 0.45497 ; label("Proportional residual error on urinary glucose excretion (fraction)")    # Suppl Table 2 T2DM, err_prop 0.207 (%CV 45.5)
  })

  model({
    linear(GLU)

    # Physical constants, not fitted quantities. mwEmpa (empagliflozin molar
    # mass, from the formula C23H27ClO7) is NOT reported in the paper and is
    # needed only to dose in mg while reporting concentrations in nmol/L;
    # glucoseConv converts the register's canonical mmol/L to the paper's
    # mg/dL.
    mwEmpa      <- 450.91   # g/mol
    glucoseConv <- 18.016   # (mg/dL) per (mmol/L)

    cl   <- exp(lcl + etalcl)     * (BMI / 25)^e_bmi_cl
    vc   <- exp(lvc)              * (BMI / 25)^e_bmi_vc
    q    <- exp(lq + etalq)       * (BMI / 25)^e_bmi_q
    vp   <- exp(lvp + etalvp)     * (BMI / 25)^e_bmi_vp
    ka   <- exp(lka + etalka)
    d1   <- exp(ld1 + etald1)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    dur(depot)  <- d1
    alag(depot) <- tlag

    Cc <- central / vc * 1e6 / mwEmpa

    # Reference age is 58 years for the T2DM model (Suppl Table 2 footnote).
    gmax <- exp(lgmax + etalgmax) * e_sexf_gmax^SEXF * (AGE / 58)^e_age_gmax
    km   <- exp(lkm)
    ic50 <- exp(lic50)
    imax <- expit(logitimax) * e_sexf_imax^SEXF * (AGE / 58)^e_age_imax
    frac <- expit(logitfrac + etalogitfrac)

    iPG  <- GLU * glucoseConv
    egfr <- CRCL * BSA / 1.73 * 0.6

    inhib  <- imax * Cc / (ic50 + Cc)
    reabsG <- gmax * iPG / (km + iPG) * (1 - inhib)

    if (iPG < reabsG) {
      ugeRate <- egfr * iPG * (1 - frac)
    } else {
      ugeRate <- egfr * (iPG - reabsG)
    }

    d/dt(glu_urine) <- ugeRate

    UGE <- glu_urine
    RTG <- gmax * (1 - inhib) - km

    Cc  ~ prop(propSd)
    UGE ~ prop(propSd_UGE)
  })
}
