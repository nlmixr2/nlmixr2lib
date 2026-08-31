Mondick_2018_empagliflozin_t1dm <- function() {
  description <- paste(
    "Joint population PK-PD model of empagliflozin in adults with type 1",
    "diabetes mellitus. Two-compartment oral PK with sequential zero- and",
    "first-order absorption and an absorption lag time, coupled to a",
    "mechanistic renal glucose-reabsorption model in which urinary glucose",
    "excretion is the difference between the glomerular filtration of plasma",
    "glucose and a saturable (Michaelis-Menten) tubular reabsorption capacity",
    "that empagliflozin inhibits through an Emax function, thereby lowering",
    "the renal threshold for glucose (RTG).",
    sep = " "
  )
  reference <- paste(
    "Mondick J, Riggs M, Kaspers S, Soleymanlou N, Marquard J, Nock V.",
    "Population pharmacokinetic-pharmacodynamic analysis to characterize the",
    "effect of empagliflozin on renal glucose threshold in patients with type 1",
    "diabetes mellitus. J Clin Pharmacol. 2018;58(5):640-649.",
    "doi:10.1002/jcph.1051.",
    "The PD structure follows the type 2 diabetes model of Mondick J, Riggs M,",
    "Sasaki T, Sarashina A, Broedl UC, Retlich S. Diabetes Obes Metab.",
    "2016;18:241-248; see modellib('Mondick_2018_empagliflozin_t2dm') for that",
    "population's parameter set as reproduced in Mondick 2018 Supplementary",
    "Table 2.",
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
        "Vp/F (Mondick 2018 Population PK Modeling equations). Cohort mean (SD)",
        "25.6 (3.69) kg/m^2. Body size enters the model only through BMI; there",
        "is no separate allometric weight term."
      ),
      source_name        = "BMI"
    ),
    CRCL = list(
      description        = paste(
        "BSA-normalized estimated glomerular filtration rate (eGFR). Mondick",
        "2018 does not name the estimating equation; the value is reported as",
        "mL/(min * 1.73 m^2)."
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "NOT a covariate on a PK parameter. It is a mechanistic driver: the",
        "model de-normalizes it to an absolute filtration rate using the",
        "individual BSA (Mondick 2018: 'eGFR values were converted from",
        "mL/(min * 1.73 m2) to mL/min using the individual body surface",
        "values') and multiplies it by the plasma-minus-reabsorbed glucose",
        "concentration difference to give the urinary glucose excretion rate.",
        "Cohort mean (SD) 102 (13.8) mL/min/1.73 m^2."
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
        "filtration rate, as CRCL * BSA / 1.73. Mondick 2018 does not report",
        "the BSA formula used nor the cohort BSA distribution."
      ),
      source_name        = "BSA"
    ),
    GLU = list(
      description        = "Plasma glucose concentration, time-varying regressor",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Exogenous driving regressor, not a covariate on a parameter. Mondick",
        "2018 fitted eight-point daily plasma-glucose profiles and linearly",
        "interpolated between observed time points ('iPG'); the model declares",
        "linear(GLU) so rxode2 reproduces that interpolation. The paper works",
        "in mg/dL throughout; the model converts mmol/L -> mg/dL internally",
        "(factor 18.016), matching the register's canonical unit and the",
        "Bosch_2025_glp1ra_hba1c.R precedent."
      ),
      source_name        = "PG"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-scaled against the 42-year reference patient on both Gmax and",
        "Imax (Mondick 2018 Supplementary Table 2 footnote). Cohort mean (SD)",
        "41.0 (10.9) years; study range 20-60 years. Both age exponents are",
        "estimated with 95% CIs spanning zero, so the paper concludes no",
        "age effect is identifiable; the terms are retained because the paper",
        "reports them as part of a full covariate model."
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "categorical",
      reference_category = "0 = male (the reference patient is a 42-year-old male)",
      notes              = paste(
        "Enters as a multiplicative factor raised to the SEXF indicator on both",
        "Gmax and Imax (Mondick 2018 Supplementary Table 2). Cohort 22/75",
        "female (29.3%)."
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
    n_subjects     = 75,
    n_studies      = 1,
    age_range      = "20-60 years",
    age_median     = "41.0 years (mean; SD 10.9). Reference patient 42 years",
    weight_median  = "79.3 kg (mean; SD 14.3)",
    bmi_median     = "25.6 kg/m^2 (mean; SD 3.69)",
    sex_female_pct = 29.3,
    renal_function = "eGFR mean 102 mL/min/1.73 m^2 (SD 13.8)",
    disease_state  = "type 1 diabetes mellitus, on background insulin",
    dose_range     = "placebo, or empagliflozin 2.5, 10 or 25 mg orally once daily for 28 days, as adjunct to insulin",
    regions        = "Germany, Austria",
    notes          = paste(
      "EASE-1 (NCT01969747), a randomized, placebo-controlled, parallel-group",
      "phase 2 study run at one German and one Austrian centre between",
      "25 November 2013 and 20 April 2014. Randomization 1:1:1:1 to placebo",
      "(n = 19), 2.5 mg (n = 19), 10 mg (n = 19) and 25 mg (n = 18). The",
      "population PK data set held 1814 empagliflozin concentrations and the PD",
      "data set 900 urinary glucose excretion records. Baseline demographics",
      "are in Mondick 2018 Results, 'Analysis Population'. The PK and PD models",
      "were fitted SEQUENTIALLY: PK by FOCE with eta-epsilon interaction, then",
      "the PD model by SAEM with the individual empirical Bayes PK estimates",
      "fixed as input."
    )
  )

  ini({
    # ---- Population PK: Mondick 2018 Supplementary Table 1 ------------------
    # Structural values are apparent (oral) quantities at the BMI = 25 kg/m^2
    # reference. IIV entries are NONMEM OMEGA variances on the log scale; the
    # table's "%CV" column is sqrt(exp(omega^2) - 1), verified for every row
    # (e.g. CL/F: sqrt(exp(0.0387) - 1) = 0.199 = the tabulated 19.9%).
    lcl   <- log(12.3)  ; label("Apparent oral clearance CL/F at BMI 25 kg/m^2 (L/h)")                # Suppl Table 1, 12.3 L/h (%RSE 1.14)
    lvc   <- log(3.47)  ; label("Apparent central volume Vc/F at BMI 25 kg/m^2 (L)")                  # Suppl Table 1, 3.47 L (%RSE 7.10)
    lq    <- log(7.95)  ; label("Apparent intercompartmental clearance Q/F at BMI 25 kg/m^2 (L/h)")   # Suppl Table 1, 7.95 L/h (%RSE 3.23)
    lvp   <- log(88.0)  ; label("Apparent peripheral volume Vp/F at BMI 25 kg/m^2 (L)")               # Suppl Table 1, 88.0 L (%RSE 0.952)
    lka   <- log(0.275) ; label("First-order absorption rate constant ka (1/h)")                      # Suppl Table 1, 0.275 1/h (%RSE 3.17)
    ld1   <- log(0.607) ; label("Duration of the zero-order input into the depot D1 (h)")             # Suppl Table 1, 0.607 h (%RSE 15.1); tabulated as "zero-order absorption rate" but reported in hours, i.e. a duration
    ltlag <- log(0.189) ; label("Absorption lag time ALAG1 (h)")                                      # Suppl Table 1, 0.189 h (%RSE 5.59)

    # BMI power exponents, reference 25 kg/m^2. The parameter <-> value binding
    # is taken from the Supplementary Table 1 row LABELS, which name the
    # parameter and its value on the same line. Mondick 2018's main-text PK
    # equations number these theta_10 -> Vp/F and theta_11 -> Q/F, while
    # Supplementary Table 1 numbers theta_10 -> Q/F and theta_11 -> Vp/F. Only
    # the index labels disagree; both sources agree on which VALUE belongs to
    # which parameter, because the table states it directly. See vignette
    # Errata.
    e_bmi_cl <-  0.554 ; label("Power exponent on (BMI/25) for CL/F (unitless)")  # Suppl Table 1, BMI effect on CL/F 0.554 (%RSE 32.6)
    e_bmi_vc <- -0.241 ; label("Power exponent on (BMI/25) for Vc/F (unitless)")  # Suppl Table 1, BMI effect on Vc/F -0.241 (%RSE 223)
    e_bmi_q  <-  1.77  ; label("Power exponent on (BMI/25) for Q/F (unitless)")   # Suppl Table 1, BMI effect on Q/F 1.77 (%RSE 16.4)
    e_bmi_vp <-  0.775 ; label("Power exponent on (BMI/25) for Vp/F (unitless)")  # Suppl Table 1, BMI effect on Vp/F 0.775 (%RSE 49.9)

    # ---- Population PD: Mondick 2018 Supplementary Table 2, T1DM column -----
    lgmax      <- log(317)  ; label("Maximum reabsorbed glucose concentration Gmax (mg/dL)")                       # Suppl Table 2, 317 mg/dL (95% CI 286, 359)
    lkm        <- log(136)  ; label("Plasma glucose concentration at half-maximal reabsorption KM (mg/dL)")        # Suppl Table 2, 136 mg/dL (95% CI 103, 181)
    lic50      <- log(5.84) ; label("Empagliflozin concentration at half-maximal reabsorption inhibition IC50 (nmol/L)")  # Suppl Table 2, 5.84 nmol/L (95% CI 4.65, 7.71)

    # FRAC and Imax are bounded to (0, 1) by a logistic transform in Mondick
    # 2018 (Methods, the two equations following equation 3), so both are
    # carried on the logit scale here. The paper writes lambda_FRAC = ln(theta/
    # (1 - theta)) + eta and FRAC = expit(lambda_FRAC), so theta IS the typical
    # value on the natural scale; likewise for Imax, whose printed logistic
    # reduces to the identity. NOTE: the Supplementary Table 2 "Model" column
    # renders these denominators as (1 + theta); the main-text equations on
    # p.643 show (1 - theta), which is the form that reproduces the tabulated
    # estimates (expit(logit(0.995)) = 0.995) and is therefore used.
    logitimax  <- log(0.676 / (1 - 0.676)) ; label("Logit of maximum inhibition of renal glucose reabsorption Imax (unitless)")  # Suppl Table 2, Imax 0.676 (95% CI 0.648, 0.707)
    logitfrac  <- log(0.995 / (1 - 0.995)) ; label("Logit of the fraction of filtered glucose reabsorbed below RTG, FRAC (unitless)")  # Suppl Table 2, FRAC 0.995 (95% CI 0.992, 0.997)

    # Full covariate model on Gmax and Imax, referenced to a 42-year-old male.
    # Both AGE exponents have 95% CIs spanning zero; Mondick 2018 concludes no
    # age effect could be identified but reports the estimates, so they are
    # retained here at their published point values.
    e_sexf_gmax <-  1.13    ; label("Multiplicative factor on Gmax for female sex (unitless)")     # Suppl Table 2, sex effect on Gmax 1.13 (95% CI 1.05, 1.23)
    e_age_gmax  <- -0.0205  ; label("Power exponent on (AGE/42) for Gmax (unitless)")              # Suppl Table 2, age effect on Gmax -0.0205 (95% CI -0.145, 0.106)
    e_sexf_imax <-  0.901   ; label("Multiplicative factor on Imax for female sex (unitless)")     # Suppl Table 2, sex effect on Imax 0.901 (95% CI 0.841, 0.975)
    e_age_imax  <-  0.0447  ; label("Power exponent on (AGE/42) for Imax (unitless)")              # Suppl Table 2, age effect on Imax 0.0447 (95% CI -0.0664, 0.129)

    # ---- Interindividual variability ---------------------------------------
    # PK: OMEGA variances, log scale. Vc/F carries no random effect (Mondick
    # 2018: "Random effects parameters were estimated for all structural
    # parameters except Vc/F").
    etalcl   ~ 0.0387  # Suppl Table 1, IIV CL/F   0.0387 (%CV 19.9)
    etalq    ~ 0.0422  # Suppl Table 1, IIV Q/F    0.0422 (%CV 20.8)
    etalvp   ~ 0.0772  # Suppl Table 1, IIV Vp/F   0.0772 (%CV 28.3)
    etalka   ~ 0.0286  # Suppl Table 1, IIV Ka     0.0286 (%CV 17.0)
    etald1   ~ 0.174   # Suppl Table 1, IIV D1     0.174  (%CV 43.6)
    etaltlag ~ 0.309   # Suppl Table 1, IIV ALAG1  0.309  (%CV 60.2)

    # PD: a correlated 2x2 OMEGA block. etalgmax is on the log scale,
    # etalogitfrac on the logit scale, matching how each parameter is
    # transformed. Correlation check: -0.0643 / sqrt(0.0101 * 3.22) = -0.357,
    # the tabulated CORR of -0.356.
    etalgmax + etalogitfrac ~ c(0.0101,
                                -0.0643, 3.22)  # Suppl Table 2, IIV Gmax 0.0101; cov -0.0643; IIV FRAC 3.22

    # ---- Residual error ----------------------------------------------------
    # Both are NONMEM SIGMA variances in the source; propSd is a standard
    # deviation, so each is entered as sqrt(sigma^2). Check: sqrt(0.0584) =
    # 0.2417 = the tabulated 24.2% CV; sqrt(0.169) = 0.4111 = 41.1%.
    propSd     <- 0.24166 ; label("Proportional residual error on empagliflozin concentration (fraction)")  # Suppl Table 1, err_prop 0.0584 (%CV 24.2)
    propSd_UGE <- 0.41110 ; label("Proportional residual error on urinary glucose excretion (fraction)")    # Suppl Table 2, err_prop 0.169 (%CV 41.1)
  })

  model({
    # Plasma glucose is an exogenous, linearly interpolated driving regressor
    # ("iPG" in Mondick 2018), matching the paper's linear interpolation of the
    # eight-point daily profiles between observed time points.
    linear(GLU)

    # ---- Unit constants ----------------------------------------------------
    # Neither constant is a fitted quantity; both are physical constants
    # supplying conversions that Mondick 2018 performs implicitly.
    #   mwEmpa   -- empagliflozin molar mass, 450.91 g/mol. NOT reported in the
    #               paper; taken from the compound's molecular formula
    #               C23H27ClO7. Needed only because doses are administered in
    #               mg while the paper reports concentrations and IC50 in
    #               nmol/L. CL/F and V/F are unit-invariant, so converting at
    #               the observation is exactly equivalent to dosing in nmol.
    #   glucoseConv -- 1 mmol/L glucose = 18.016 mg/dL (glucose molar mass
    #               180.16 g/mol). The paper works in mg/dL; the GLU register
    #               canonical unit is mmol/L.
    mwEmpa      <- 450.91   # g/mol
    glucoseConv <- 18.016   # (mg/dL) per (mmol/L)

    # ---- Individual PK parameters ------------------------------------------
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

    # ---- PK ODEs -----------------------------------------------------------
    # Sequential zero- then first-order absorption: the dose enters `depot` as
    # a zero-order input of duration D1 after a lag of ALAG1, then leaves it
    # first-order at rate ka.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    dur(depot)  <- d1
    alag(depot) <- tlag

    # Plasma empagliflozin in nmol/L: (mg / L) * (1e6 nmol/mmol / (mg/mmol)).
    Cc <- central / vc * 1e6 / mwEmpa

    # ---- Individual PD parameters ------------------------------------------
    # Covariate factors are applied on the natural scale, exactly as printed in
    # Supplementary Table 2 (Gmax * theta^SEXF * (AGE/42)^theta).
    gmax <- exp(lgmax + etalgmax) * e_sexf_gmax^SEXF * (AGE / 42)^e_age_gmax
    km   <- exp(lkm)
    ic50 <- exp(lic50)
    imax <- expit(logitimax) * e_sexf_imax^SEXF * (AGE / 42)^e_age_imax
    frac <- expit(logitfrac + etalogitfrac)

    # ---- Renal glucose handling --------------------------------------------
    # iPG: interpolated plasma glucose in mg/dL.
    # egfr: de-normalize the BSA-indexed eGFR to this subject's absolute
    #   filtration rate and express it in dL/h, so that egfr * (mg/dL) is a
    #   glucose mass rate in mg/h.
    #   (mL/min/1.73 m^2) * (m^2 / 1.73) = mL/min; * 60/100 = dL/h.
    iPG  <- GLU * glucoseConv
    egfr <- CRCL * BSA / 1.73 * 0.6

    # Saturable tubular reabsorption capacity with empagliflozin inhibition
    # (Mondick 2018 equation 3). ReabsG has units of mg/dL.
    inhib  <- imax * Cc / (ic50 + Cc)
    reabsG <- gmax * iPG / (km + iPG) * (1 - inhib)

    # Mondick 2018 equations 1 and 2. Below the reabsorption capacity, only the
    # small non-reabsorbed fraction (1 - FRAC) of the filtered load is
    # excreted; at or above it, reabsorption saturates and the excess is
    # excreted.
    if (iPG < reabsG) {
      ugeRate <- egfr * iPG * (1 - frac)
    } else {
      ugeRate <- egfr * (iPG - reabsG)
    }

    d/dt(glu_urine) <- ugeRate

    # Cumulative urinary glucose excreted, in mg. Mondick 2018 observed 24-hour
    # collections; difference this state across a collection interval to obtain
    # the interval amount.
    UGE <- glu_urine

    # Renal threshold for glucose: the plasma glucose concentration at which
    # d[UGE]/dt = 0, obtained in Mondick 2018 by solving equation 2 for PG.
    # Reported as a derived quantity, not fitted.
    RTG <- gmax * (1 - inhib) - km

    Cc  ~ prop(propSd)
    UGE ~ prop(propSd_UGE)
  })
}
