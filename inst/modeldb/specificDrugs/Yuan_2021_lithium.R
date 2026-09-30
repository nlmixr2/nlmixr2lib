Yuan_2021_lithium <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model for lithium after a",
    "single oral dose of lithium carbonate (12 mg/kg) in 52 Chinese children",
    "aged 4-10 years with intellectual disability, from Yuan 2021. Absorption",
    "is a chain of six transit compartments with a common transit rate",
    "constant ktr = (6 + 1) / MTT; relative bioavailability is fixed to unity",
    "with between-subject variability. Body weight enters as a fixed",
    "allometric function (exponent 0.75 on CL/F and Q/F, 1 on Vc/F and",
    "Vp/F) centred on the study-median 20 kg. Doses are in mmol of lithium",
    "ion (1 mg lithium carbonate = 2 / 73.89 mmol Li) and concentrations in",
    "mmol/L; the fit was to baseline-subtracted (pre-dose endogenous lithium",
    "removed) serum concentrations.",
    sep = " "
  )
  reference <- paste(
    "Yuan J, Zhang B, Xu Y, Zhang X, Song J, Zhou W, Hu K, Zhu D, Zhang L,",
    "Shao F, Zhang S, Ding J, Zhu C.",
    "Population Pharmacokinetics of Lithium in Young Pediatric Patients With",
    "Intellectual Disability.",
    "Front Pharmacol. 2021;12:650298. doi:10.3389/fphar.2021.650298.",
    "PMC8082156. Parameter estimates are in Table 1; the structural model is",
    "Figure 1; covariate equations are Methods Equations 5-7.",
    sep = " "
  )
  vignette <- "Yuan_2021_lithium"
  units <- list(time = "h", dosing = "mmol", concentration = "mmol/L")

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fixed allometric function on all clearance parameters (CL/F, Q/F;",
        "exponent 0.75, Methods Equation 5) and all volume parameters (Vc/F,",
        "Vp/F; exponent 1, Equation 6), centred on the study-median body",
        "weight of 20 kg (Methods, 'Covariate Modeling'). Recorded before",
        "the single lithium dose (Methods, 'Medication Dosing'). Set WT = 20",
        "to recover the Table 1 typical values.",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  # Covariates the source screened and did not retain in the final model.
  # Documentation only -- none is referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested as a saturable maturation function on CL,",
        "CL = CL_STD * AGE / (AGE50 + AGE) * (BW / 20)^0.75 (Methods",
        "Equation 7). 'Including age-dependent maturation on clearance did",
        "not improve model fit further' (Results, 'Population PK'); no AGE50",
        "estimate is reported.",
        sep = " "
      )
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Tested as a proportional effect on all parameters with male = 0,",
        "female = 1 (Methods Equation 8); 'gender did not have a significant",
        "impact on lithium PK properties' (Results). No coefficient reported.",
        sep = " "
      )
    )
  )

  compartmentData <- list(
    depot = list(analyte = "lithium", units = "mmol", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "lithium", units = "mmol", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "lithium", units = "mmol", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "lithium", units = "mmol", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "lithium", units = "mmol", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "lithium", units = "mmol", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "lithium", units = "mmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lithium", units = "mmol", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "lithium", units = "mmol", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 52L,
    n_studies = 1L,
    age_range = "48-128 months (4-10 years)",
    age_mean = "84.8 +/- 21.7 months",
    weight_range = "16-44 kg",
    weight_mean = "23.0 +/- 6.2 kg",
    weight_median = "20 kg",
    sex_female_pct = 26.9,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Intellectual disability (DSM-5; IQ < 70) with normal liver, renal",
      "and thyroid function.",
      sep = " "
    ),
    dose_range = "Single oral dose of lithium carbonate 12 mg/kg after a fast of at least 4 h.",
    regions = "Zhengzhou, Henan, China",
    n_observations = 382L,
    notes = paste(
      "Single-centre study at the Third Affiliated Hospital of Zhengzhou",
      "University. 16 children (8 boys, 8 girls) gave an intensive profile",
      "(0.5, 1, 1.5, 2, 4, 8, 12, 24, 36 and 48 h); 36 children (30 boys,",
      "6 girls) gave at least three samples from the same time points.",
      "Serum lithium by ion chromatography (LLOQ 0.00144 mmol/L; no sample",
      "below it). The pre-dose endogenous lithium concentration was",
      "subtracted from every post-dose concentration before fitting.",
      "NONMEM 7.4, FOCE-I, log-transformed data. Demographics are given in",
      "Results, 'Demographics' (38 male, 14 female).",
      sep = " "
    )
  )

  ini({
    # All values: Yuan 2021 Table 1 ('Final population PK parameter
    # estimates of lithium in children'), column 'NONMEM estimates
    # (%RSE)', for a typical child of 20 kg. All disposition parameters
    # are apparent (relative to F, fixed to unity).

    # ---- Absorption ------------------------------------------------
    lmtt <- log(0.52)
    label("Mean transit time through the depot plus six-transit-compartment absorption chain MTT (h)")
    # Table 1: MTT (h) = 0.52 (%RSE 9.9)
    # Table 1: 'Number of transit compartment' = 6 fixed (integer, used as
    # the literal 6 in model()).

    # ---- Disposition -----------------------------------------------
    lcl <- log(0.98)
    label("Apparent elimination clearance CL/F at WT = 20 kg (L/h)")
    # Table 1: CL/F (L/h) = 0.98 (%RSE 4.6)

    lvc <- log(13.1)
    label("Apparent central volume of distribution Vc/F at WT = 20 kg (L)")
    # Table 1: VC/F (L) = 13.1 (%RSE 7.2)

    lq <- log(0.84)
    label("Apparent inter-compartmental clearance Q/F at WT = 20 kg (L/h)")
    # Table 1: Q/F (L/h) = 0.84 (%RSE 9.5)

    lvp <- log(8.2)
    label("Apparent peripheral volume of distribution Vp/F at WT = 20 kg (L)")
    # Table 1: Vp/F (L) = 8.2 (%RSE 17.7)

    # ---- Relative bioavailability ----------------------------------
    lfdepot <- fixed(log(1))
    label("Relative bioavailability F (unitless)")
    # Table 1: F (%) = '100 fixed'; Methods: 'Relative bioavailability (F)
    # was fixed to unity in the population in order to allow investigation
    # of the IIV of absorption'.

    # ---- Allometric scaling (fixed, Methods Equations 5 and 6) -------
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F and Q/F, referenced to 20 kg (unitless)")
    # Methods Equation 5: exponent 0.75 on clearance parameters

    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on Vc/F and Vp/F, referenced to 20 kg (unitless)")
    # Methods Equation 6: exponent 1 on volume parameters

    # ---- Between-subject variability -------------------------------
    # Table 1 column 'CV for IIV (%RSE)'. Exponential IIV (Methods
    # Equation 1), so omega^2 = log(CV^2 + 1). No IIV on CL/F or Q/F:
    # removed from the final model (CL IIV RSE 195%; Q IIV near zero;
    # Results, 'Population PK').
    etalfdepot ~ 0.0878360
    # Table 1: F IIV = 30.3% CV (%RSE 10.9); log(0.303^2 + 1) = 0.0878360
    etalmtt ~ 0.3524159
    # Table 1: MTT IIV = 65.0% CV (%RSE 13.1); log(0.650^2 + 1) = 0.3524159
    etalvc ~ 0.0678689
    # Table 1: VC/F IIV = 26.5% CV (%RSE 16.1); log(0.265^2 + 1) = 0.0678689
    etalvp ~ 0.8791989
    # Table 1: Vp/F IIV = 118.7% CV (%RSE 15.4); log(1.187^2 + 1) = 0.8791989

    # ---- Residual error --------------------------------------------
    # Additive error on log-transformed concentrations (Methods Equation
    # 2, 'equivalent to an exponential residual error on an arithmetic
    # scale'), i.e. a log-normal residual. Table 1 prints sigma (the
    # Methods define the variance as sigma^2), so 0.091 is the log-scale SD.
    expSd <- 0.091
    label("Log-scale residual standard deviation for serum lithium (unitless)")
    # Table 1: sigma = 0.091 (%RSE 19.0), 'additive residue error on a log scale'
  })

  model({
    # ---- Allometric scaling on the study-median 20 kg ---------------
    allom_cl <- (WT / 20)^e_wt_cl
    allom_v <- (WT / 20)^e_wt_vc

    # ---- Individual parameters --------------------------------------
    mtt <- exp(lmtt + etalmtt)
    cl <- exp(lcl) * allom_cl
    vc <- exp(lvc + etalvc) * allom_v
    q <- exp(lq) * allom_cl
    vp <- exp(lvp + etalvp) * allom_v
    fdepot <- exp(lfdepot + etalfdepot)

    # Figure 1: the dose enters through ktr into a1, then a1 -> ... -> a6
    # -> central, all at the same ktr (no separately estimated ka). The
    # dosing compartment plus six transit compartments give seven
    # sequential transfers, so ktr = (n + 1) / MTT with n = 6 (Savic 2007
    # transit-compartment convention, cited in the Methods).
    ktr <- (6 + 1) / mtt

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- ODE system --------------------------------------------------
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ktr * transit6
    d/dt(central) <- ktr * transit6 - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # ---- Observation -------------------------------------------------
    # Doses in mmol lithium and volumes in L, so central / vc is mmol/L.
    Cc <- central / vc

    Cc ~ lnorm(expSd)
  })
}
