Zhang_2022_linezolid <- function() {
  description <- paste(
    "Two-compartment population PK model for oral linezolid in Chinese adults",
    "with multidrug-resistant tuberculosis (Zhang 2022), with first-order",
    "absorption after a lag time and linear elimination. Body weight scales CL",
    "(exponent 0.75) and the central volume (exponent 1.0) about an assumed",
    "70 kg reference (the paper prints no reference weight); type 2 diabetes",
    "adds an additive shift to CL and to the central volume. Residual error is",
    "additive. The intercompartmental clearance is encoded as printed",
    "(Q = 0.9 L/h); with it the model under-predicts the paper's own posthoc",
    "trough concentrations about twofold, and a value of 9 L/h would reproduce",
    "them (see the vignette Errata).",
    sep = " "
  )
  reference <- paste(
    "Zhang H, He Y, Davies Forsman L, Paues J, Werngren J, Niward K, Schon T,",
    "Bruchfeld J, Alffenaar J-W, Hu Y (2023). Population pharmacokinetics and",
    "dose evaluations of linezolid in the treatment of multidrug-resistant",
    "tuberculosis. Frontiers in Pharmacology 13:1032674.",
    "doi:10.3389/fphar.2022.1032674.",
    "Parameter estimates are from Table 2 (final model, Phoenix NLME 8.0).",
    "The covariate equation forms and the reference weight are not printed.",
    sep = " "
  )
  vignette <- "Zhang_2022_linezolid"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L"
  )

  compartmentData <- list(
    depot = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "linezolid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "linezolid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power model on CL (exponent 0.75) and the central volume (exponent",
        "1.0), both printed without an RSE in Zhang 2022 Table 2 and so held",
        "constant. The reference weight is not printed; 70 kg is assumed here.",
        "The development-cohort median is 53 kg (IQR 47-66) and the",
        "validation-cohort median 59 kg (IQR 53-64) (Table 1).",
        sep = " "
      ),
      source_name = "Bodyweight"
    ),
    DIS_DIAB = list(
      description = "Type 2 diabetes mellitus indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no type 2 diabetes",
      notes = paste(
        "Type 2 diabetes ('Diabetes type 2', Zhang 2022 Table 1 and Results",
        "3.3). Additive shifts of +0.1 L/h on CL and +5.0 L on the central",
        "volume (Table 2 theta(Diabetes)-CL and -Vd). The equation form is not",
        "printed; a 5.0 coefficient is implausible as an exponential (148-fold),",
        "multiplicative (5-fold) or fractional (6-fold) effect, so both are",
        "encoded as additive shifts in the parameter's own units. Prevalence",
        "19.0% (development) and 17.0% (validation).",
        sep = " "
      ),
      source_name = "Diabetes type 2"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Listed as a candidate effect on CL in Zhang 2022 Table 2 ('theta (Age)-CL') with no estimate, i.e. not retained in the final model.",
      source_name = "Age"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 168,
    n_studies = 1,
    n_observations = 1008,
    age_range = "18 to < 70 years (inclusion); development-cohort median 41 (IQR 33-45)",
    weight_range = "development-cohort median 53 kg (IQR 47-66)",
    sex_female_pct = 25.6,
    race_ethnicity = "Chinese",
    disease_state = "Laboratory-confirmed multidrug-resistant tuberculosis, excluding critically ill, pregnant and HIV-, HBV- or HCV-infected patients",
    dose_range = "Linezolid 600 mg orally once daily",
    regions = "China (Jiangsu, Guizhou and Sichuan provinces)",
    notes = paste(
      "Development cohort of 168 participants from a prospective multicentre",
      "cohort (Zheng 2021), sampled after 2 weeks of inpatient treatment at",
      "predose and 1, 2, 4, 6 and 8 h post-dose. Co-medication was",
      "bedaquiline, a fluoroquinolone and a background regimen. A validation",
      "cohort of 159 participants (Sichuan, Jiangsu, Henan; median weight",
      "59 kg) was sampled predose, 2 and 6 h and analysed by Bayesian",
      "estimation. Type 2 diabetes 19.0% (Zhang 2022 Table 1). Assay range",
      "0.05-30 mg/L.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- Zhang 2022 Table 2 'Typical value of
    # Population parameters' (Frontiers prints values without the leading
    # zero, e.g. '.9').
    lka <- log(2.0)
    label("Absorption rate constant (1/h)")
    # Table 2 'Ka (/h)' = 2.0 (RSD 11.1%)

    lcl <- log(5.6)
    label("Clearance for a 70 kg participant without type 2 diabetes (L/h)")
    # Table 2 'CL (L/h)' = 5.6 (RSD 1.6%)

    lvc <- log(35.8)
    label("Central volume of distribution for a 70 kg participant without type 2 diabetes (L)")
    # Table 2 'V d (L)' = 35.8 (RSD 2.0%)

    lq <- log(0.9)
    label("Intercompartmental clearance (L/h)")
    # Table 2 'Q (L/h)' = '.9' (RSD 13.7%); encoded as printed although the paper's own Table 3 troughs and Figure 6 favour 9 L/h (vignette Errata)

    lvp <- log(58.6)
    label("Peripheral volume of distribution (L)")
    # Table 2 'Vp(L)' = 58.6 (RSD 13.2%)

    ltlag <- log(0.6)
    label("Absorption lag time (h)")
    # Table 2 'Tlag(h)' = '.6' (RSD 8.4%)

    # Covariate effects -- Table 2 'Covariate'
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL (unitless)")
    # Table 2 'theta (Bodyweight)-CL' = .75, printed without an RSD

    e_wt_vc <- fixed(1.0)
    label("Allometric exponent of body weight on the central volume (unitless)")
    # Table 2 'theta (Bodyweight)-Vd' = 1.0, printed without an RSD

    e_dis_diab_cl <- 0.1
    label("Additive shift in CL with type 2 diabetes (L/h)")
    # Table 2 'theta (Diabetes)-CL' = .1 (RSD 10.5%); additive form assumed

    e_dis_diab_vc <- 5.0
    label("Additive shift in the central volume with type 2 diabetes (L)")
    # Table 2 'theta (Diabetes)-Vd' = 5.0 (RSD 11.1%); additive form assumed

    # Between-subject variability -- Table 2 'Variation of parameters between
    # individuals' (%CV), converted as omega^2 = log(CV^2 + 1). No IIV on Tlag.
    etalcl ~ 0.045600 # Table 2 CL IIV 21.6 %CV (RSD 14.6%); log(0.216^2 + 1)
    etalvc ~ 0.057831 # Table 2 Vd IIV 24.4 %CV (RSD 11.1%); log(0.244^2 + 1)

    # Residual error -- Table 2 'Additional residual (mg/L)' (additive).
    addSd <- 0.5
    label("Additive residual error (mg/L)")
    # Table 2 'Additional residual (mg/L)' = .5 (RSD 4.2%)
  })

  model({
    # Covariate model: power scaling of CL and Vc by weight (exponents from
    # Table 2) about an assumed 70 kg reference, plus additive type 2
    # diabetes shifts in L/h and L. The paper prints neither form.
    cl <- (exp(lcl) * (WT / 70)^e_wt_cl + e_dis_diab_cl * DIS_DIAB) * exp(etalcl)
    vc <- (exp(lvc) * (WT / 70)^e_wt_vc + e_dis_diab_vc * DIS_DIAB) * exp(etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka)
    tlag <- exp(ltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    # Dose in mg and volume in L give mg/L.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
