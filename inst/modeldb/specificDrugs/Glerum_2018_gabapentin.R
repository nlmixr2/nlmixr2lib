Glerum_2018_gabapentin <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption and an absorption lag time for",
    "single 800 mg oral immediate-release gabapentin in healthy adults, fitted nonparametrically",
    "(Pmetrics NPAG) to a four-way crossover comparative bioavailability study of the originator and",
    "three generic formulations. The elimination rate constant scales linearly through the origin",
    "on Cockcroft-Gault creatinine clearance. The source estimated every parameter separately per",
    "formulation ('96 pseudo subjects') but published only one pooled set of weighted medians and",
    "the post-hoc between-occasion CVs; these are encoded as typical values plus 4-occasion",
    "inter-occasion variability with no inter-individual variability. WARNING: the published",
    "medians reproduce neither the paper's own observed NCA (Cmax and AUC about 2-fold too low,",
    "terminal half-life about 3-fold too short) nor its VPC; see the vignette before using it.",
    sep = " "
  )
  reference <- paste(
    "Glerum PJ, Yu Y, Yamada WM, Neely MN, Maliepaard M, Burger DM, Neef C.",
    "Interchangeability of Generic Drugs: A Nonparametric Pharmacokinetic Model of Gabapentin",
    "Generic Drugs. Clin Pharmacol Ther. 2018;104(5):966. doi:10.1002/cpt.1023.",
    "Structure from Equations 2-3 and Methods 'Model building in Pmetrics'; parameter values from",
    "Table 3; error model from Methods and Discussion 'Model uncertainties'.",
    sep = " "
  )
  vignette <- "Glerum_2018_gabapentin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Results Equation 2: Ke = Kr * CrCl, a through-origin proportionality. Cohort mean 120",
        "mL/min (range 71.2-183.3; Methods 'Gabapentin pharmacokinetic data'). Table 3 prints Ke",
        "'recalculated from rate per ml/min creatinine clearance' without stating the CrCl used;",
        "the cohort mean 120 mL/min is taken as that reference, so the encoded form is",
        "kel = 0.32 * (CRCL / 120), i.e. Kr = 0.32 / 120 = 0.00267 per h per mL/min. This lies",
        "inside the Pmetrics search range for Kr (upper bound 0.008 h.ml-1, Methods). The",
        "Cockcroft-Gault weight convention and serum-creatinine unit are not stated.",
        sep = " "
      ),
      source_name = "CrCl"
    ),
    OCC = list(
      description = "Crossover period (occasion 1-4) for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "One occasion per single-dose period of the four-way crossover; each period is one",
        "formulation, so in the source design occasion and formulation are confounded (Discussion).",
        "Values outside 1-4 switch inter-occasion variability off.",
        sep = " "
      ),
      source_name = "pseudo subject (subject x period)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Tested as a multiplicative term on V (Equation 1, V0 = V * WEIGHT) and rejected:",
        "AIC +25 alone and +12 with renal function (Table 1). Cohort mean 71.0 kg, range",
        "52.7-97.0 kg.",
        sep = " "
      )
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened by linear regression against posterior parameters (Methods) and not taken forward; it enters only through the Cockcroft-Gault CRCL."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened by linear regression against posterior parameters (Methods) and not taken forward."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Screened by linear regression against posterior parameters (Methods) and not taken forward."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "not stated",
      type = "continuous",
      notes = "Screened by linear regression (predose screening and follow-up values; Methods); enters the final model only through the Cockcroft-Gault CRCL."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "gabapentin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "gabapentin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "gabapentin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "21-55 years",
    age_mean = "35 years",
    weight_range = "52.7-97.0 kg",
    weight_mean = "71.0 kg",
    sex_female_pct = 58.3,
    race_ethnicity = "Not reported",
    disease_state = "Healthy volunteers",
    dose_range = paste(
      "Single oral 800 mg immediate-release gabapentin in each of four periods (originator and three",
      "generics registered in the Netherlands as bioequivalent)",
      sep = " "
    ),
    regions = "Netherlands (study submitted to the Dutch Medicines Evaluation Board)",
    bmi_range = "19.9-30.0 kg/m2",
    bmi_mean = "23.6 kg/m2",
    renal_function = "Cockcroft-Gault creatinine clearance mean 120 mL/min (range 71.2-183.3)",
    n_concentrations = 1610L,
    notes = paste(
      "Methods 'Gabapentin pharmacokinetic data': 14 females and 10 males completed all four",
      "periods; plasma sampled predose and at 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 8, 12, 24,",
      "36 and 48 h. The four periods were analysed as 96 'pseudo subjects'. An external validation",
      "set (22 subjects, 800 mg single dose, 352 concentrations) was used only for prediction.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural parameters: Table 3, weighted median of the NPAG support-point
    # distribution (95% CI around the median in parentheses). The source fitted
    # these per formulation; only the pooled medians are published. Pmetrics
    # estimates a discrete joint density, encoded here by its marginal medians.
    # ------------------------------------------------------------------------
    lka <- log(0.26); label("Absorption rate constant (1/h)")
    # Table 3, Ka = 0.26 h-1 (0.23-0.31)
    ltlag <- log(0.30); label("Absorption lag time (h)")
    # Table 3, Tlag = 0.30 h (0.21-0.35)
    lvc <- log(90.0); label("Central volume of distribution (L)")
    # Table 3, V = 90.0 L (83.5-106)
    lk12 <- log(0.16); label("Central-to-peripheral rate constant KCP (1/h)")
    # Table 3, KCP = 0.16 h-1 (0.12-0.36)
    lk21 <- log(0.82); label("Peripheral-to-central rate constant KPC (1/h)")
    # Table 3, KPC = 0.82 h-1 (0.30-1.69)
    lkel <- log(0.32); label("Elimination rate constant at CRCL = 120 mL/min (1/h)")
    # Table 3, Ke = 0.32 h-1 (0.29-0.36), 'recalculated from rate per ml/min
    # creatinine clearance'; the CrCl is unstated and taken as the cohort mean
    # 120 mL/min (Methods), Kr = 0.32/120 = 0.00267 h-1 per mL/min.
    lfdepot <- fixed(log(1)); label("Oral bioavailability Fabs (fraction)")
    # Equation 3 'Fabs * dose(T - Tlag)'; no value published. Discussion
    # 'Model uncertainties': F was not estimated ('treat F as a fixed ...
    # population parameter'), so 1 is assumed.

    # ------------------------------------------------------------------------
    # Inter-occasion variability: Table 3, BOV column, 'calculated as the
    # population mean from the estimates of each subjects' four pseudo
    # subjects coefficient of variation'. A descriptive CV of the
    # within-subject estimates, converted to a lognormal variance with
    # omega^2 = log(CV^2 + 1). One draw per crossover period; occasions 2-4
    # share the occasion-1 variance. No inter-individual variability is
    # published, so none is encoded.
    # ------------------------------------------------------------------------
    etaiov_ka_1 ~ 0.1155583 # Table 3, Ka BOV 35%: log(0.35^2 + 1)
    etaiov_ka_2 ~ fixed(0.1155583) # same variance as occasion 1
    etaiov_ka_3 ~ fixed(0.1155583) # same variance as occasion 1
    etaiov_ka_4 ~ fixed(0.1155583) # same variance as occasion 1
    etaiov_tlag_1 ~ 0.4462871 # Table 3, Tlag BOV 75%: log(0.75^2 + 1)
    etaiov_tlag_2 ~ fixed(0.4462871) # same variance as occasion 1
    etaiov_tlag_3 ~ fixed(0.4462871) # same variance as occasion 1
    etaiov_tlag_4 ~ fixed(0.4462871) # same variance as occasion 1
    etaiov_vc_1 ~ 0.1093924 # Table 3, V BOV 34%: log(0.34^2 + 1)
    etaiov_vc_2 ~ fixed(0.1093924) # same variance as occasion 1
    etaiov_vc_3 ~ fixed(0.1093924) # same variance as occasion 1
    etaiov_vc_4 ~ fixed(0.1093924) # same variance as occasion 1
    etaiov_k12_1 ~ 0.7929925 # Table 3, KCP BOV 110%: log(1.10^2 + 1)
    etaiov_k12_2 ~ fixed(0.7929925) # same variance as occasion 1
    etaiov_k12_3 ~ fixed(0.7929925) # same variance as occasion 1
    etaiov_k12_4 ~ fixed(0.7929925) # same variance as occasion 1
    etaiov_k21_1 ~ 0.3800789 # Table 3, KPC BOV 68%: log(0.68^2 + 1)
    etaiov_k21_2 ~ fixed(0.3800789) # same variance as occasion 1
    etaiov_k21_3 ~ fixed(0.3800789) # same variance as occasion 1
    etaiov_k21_4 ~ fixed(0.3800789) # same variance as occasion 1
    etaiov_kel_1 ~ 0.1093924 # Table 3, Ke BOV 34%: log(0.34^2 + 1)
    etaiov_kel_2 ~ fixed(0.1093924) # same variance as occasion 1
    etaiov_kel_3 ~ fixed(0.1093924) # same variance as occasion 1
    etaiov_kel_4 ~ fixed(0.1093924) # same variance as occasion 1

    # ------------------------------------------------------------------------
    # Residual error: Pmetrics lambda model. Methods: 'polynomial equation for
    # the SD = 0.02 + 0.0075*(concentration), and error = (SD2 + Lambda2)0.5'
    # (exponents lost in typesetting). The polynomial was set by the authors
    # from other gabapentin assays (Discussion 'Model uncertainties'), hence
    # fixed(); lambda was fitted by NPAG.
    # ------------------------------------------------------------------------
    addSd <- fixed(0.02); label("Assay-error polynomial intercept C0 (mg/L)")
    # Methods 'Model building in Pmetrics', SD = 0.02 + 0.0075*C
    propSd <- fixed(0.0075); label("Assay-error polynomial slope C1 (fraction)")
    # Methods 'Model building in Pmetrics', SD = 0.02 + 0.0075*C
    lambdaSd <- 0.2990; label("Pmetrics lambda additive process noise (mg/L)")
    # Discussion 'Model uncertainties': 'The Lambda term was fitted by NPAG in
    # the final model to be 0.2990'
  })
  model({
    # 1. Inter-occasion random effects, one draw per crossover period.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 + oc4 * etaiov_ka_4
    iov_tlag <- oc1 * etaiov_tlag_1 + oc2 * etaiov_tlag_2 + oc3 * etaiov_tlag_3 + oc4 * etaiov_tlag_4
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2 + oc3 * etaiov_vc_3 + oc4 * etaiov_vc_4
    iov_k12 <- oc1 * etaiov_k12_1 + oc2 * etaiov_k12_2 + oc3 * etaiov_k12_3 + oc4 * etaiov_k12_4
    iov_k21 <- oc1 * etaiov_k21_1 + oc2 * etaiov_k21_2 + oc3 * etaiov_k21_3 + oc4 * etaiov_k21_4
    iov_kel <- oc1 * etaiov_kel_1 + oc2 * etaiov_kel_2 + oc3 * etaiov_kel_3 + oc4 * etaiov_kel_4

    # 2. Individual parameters. Equation 2: Ke = Kr * CrCl, written as the
    #    Table 3 Ke at the 120 mL/min cohort mean times CRCL / 120.
    ka <- exp(lka + iov_ka)
    tlag <- exp(ltlag + iov_tlag)
    vc <- exp(lvc + iov_vc)
    k12 <- exp(lk12 + iov_k12)
    k21 <- exp(lk21 + iov_k21)
    kel <- exp(lkel + iov_kel) * (CRCL / 120)
    fdepot <- exp(lfdepot)

    # 3. ODE system, Equation 3 with KCP as the central outflow. The printed
    #    central equation reads '- ammnt_blood * (KPC + Ke)' while the
    #    peripheral equation gains 'KCP * ammnt_blood'; that pairing does not
    #    conserve mass and contradicts the Methods definition of KCP as the
    #    central-to-peripheral rate, so KCP is used in both (see vignette).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot + k21 * peripheral1 - (k12 + kel) * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag
    f(depot) <- fdepot

    # 4. Observation and error. concentration_blood = ammnt_blood / vol_blood
    #    (mg / L = mg/L). Pmetrics lambda error:
    #    SD = sqrt((C0 + C1 * C)^2 + lambda^2).
    Cc <- central / vc
    sdCc <- sqrt((addSd + propSd * Cc)^2 + lambdaSd^2)
    Cc ~ add(sdCc)
  })
}
