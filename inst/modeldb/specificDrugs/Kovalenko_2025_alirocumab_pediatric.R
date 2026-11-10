Kovalenko_2025_alirocumab_pediatric <- function() {
  description <- "Full target-mediated drug disposition (TMDD) population PK model for total alirocumab and free PCSK9 in children and adolescents (8-17 years) with heterozygous familial hypercholesterolaemia (Kovalenko 2025, primary pediatric covariate model). Two-compartment disposition of free alirocumab, free PCSK9 and both alirocumab:PCSK9 complexes (1:1 and 1:2) sharing the same intercompartmental rate constants; SC absorption through one depot-transit compartment; two-step binding of alirocumab to one and then a second PCSK9 molecule in the central compartment; zero-order PCSK9 synthesis. Weight on Vc (power, 50 kg reference) and high-dose statin on PCSK9 synthesis rate. Absorption, distribution and bioavailability fixed to the adult estimates; kon and koff fixed to Biacore values."
  reference <- "Kovalenko P, Harnisch L, Mendell J, Wang Y, Davis JD, DiCioccio AT. Target-Mediated Modeling of Alirocumab in Adolescents and Children >=8 to <12 Years of Age Using Phase II and III Data. Clin Pharmacol Drug Dev. 2025;14(5):347. doi:10.1002/cpdd.1523"
  vignette <- "Kovalenko_2025_alirocumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # complex_db is the 1:2 alirocumab:PCSK9:PCSK9 complex (the paper's L:R:R)
  # formed when the second arm of the bivalent antibody binds a second PCSK9;
  # complex_db_peripheral1 is its peripheral-compartment counterpart. The 1:1
  # complex is the canonical `complex`.
  paper_specific_compartments <- c("complex_db", "complex_db_peripheral1")

  covariateData <- list(
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline (time-fixed) weight; the paper repeated covariate selection with time-dependent weight and kept the baseline form because time-dependent covariates gave unstable OFVs (Methods, 'Population PK Analysis'). Power effect on Vc centred at 50 kg: the dataset median 50.1 kg was replaced by 50 kg because treatment groups were assigned at the 50 kg cut-off (Methods). Analysis-set weight mean 52.0 kg (SD 18.8).",
      source_name = "weight"
    ),
    CONMED_STATIN_HI = list(
      description = "Concomitant high-dose statin therapy",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no statin or low-dose statin)",
      notes = "The paper's HDS covariate (Table 2 'ksyn ~ HDS'). The paper does not state the dose threshold that defines a high-dose statin; it reports that 19 participants received HDS and that pediatric statin doses were mostly <=20 mg once daily (Discussion). Multiplicative exponential effect on PCSK9 synthesis rate, ksyn = 5.76 * exp(0.576 * HDS), i.e. 10.3 nM/day on HDS (Table 2 derived row; Results 'ksyn was 79% higher').",
      source_name = "HDS"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in the stepwise covariate search (Methods) and not retained in the final model.",
      source_name = "age"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in the stepwise covariate search (Methods) and not retained in the final model.",
      source_name = "albumin"
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in the stepwise covariate search (Methods) and not retained in the final model.",
      source_name = "creatinine clearance"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in the stepwise covariate search (Methods) and not retained in the final model.",
      source_name = "BMI"
    ),
    ADA_POS = list(
      description = "Anti-drug antibody positive at any time",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ADA negative)",
      notes = "Tested in the stepwise covariate search; not statistically significant or clinically meaningful (Results).",
      source_name = "ADA at any time"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Tested in the stepwise covariate search (Methods) and not retained in the final model.",
      source_name = "sex"
    ),
    CONMED_EZE = list(
      description = "Concomitant ezetimibe",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no ezetimibe)",
      notes = "Tested in the stepwise covariate search (Methods) and not retained in the final model.",
      source_name = "ezetimibe"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "alirocumab", units = "nmol", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "alirocumab", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "alirocumab (free)", units = "nmol", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "alirocumab (free)", units = "nmol", specimen = "not applicable", verified = TRUE),
    target = list(analyte = "PCSK9 (free)", units = "nmol", specimen = "serum", verified = TRUE),
    target_peripheral1 = list(analyte = "PCSK9 (free)", units = "nmol", specimen = "not applicable", verified = TRUE),
    complex = list(analyte = "alirocumab:PCSK9 (1:1 complex)", units = "nmol", specimen = "serum", verified = TRUE),
    complex_peripheral1 = list(
      analyte = "alirocumab:PCSK9 (1:1 complex)",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    complex_db = list(
      analyte = "alirocumab:PCSK9:PCSK9 (1:2 complex)",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    complex_db_peripheral1 = list(
      analyte = "alirocumab:PCSK9:PCSK9 (1:2 complex)",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 119L,
    n_studies = 2L,
    n_observations = "390 total alirocumab and 531 free PCSK9 concentrations",
    age_range = "8-17 years",
    age_mean = "13.1 years (SD 2.86)",
    weight_mean = "52.0 kg (SD 18.8)",
    sex_female_pct = 57.5,
    disease_state = "Children and adolescents with heterozygous familial hypercholesterolaemia (HeFH)",
    dose_range = "Weight-banded SC alirocumab: <50 kg 30 or 40 mg Q2W, 75 or 150 mg Q4W (up-titration to 75 mg Q2W); >=50 kg 50 or 75 mg Q2W, 150 or 300 mg Q4W (up-titration to 150 mg Q2W) (Table 1)",
    co_medication = "More than 90% received statins; 19 participants received a high-dose statin (Discussion, Results)",
    studies = "Phase II NCT02890992 (DFI14223, 8-week dose finding, 42 treated) and phase III NCT03510884 (EFC14643, 99 treated, 20 also in phase II)",
    notes = "121 patients received active treatment across the two studies; one had no PK observations and one had only one post-dose observation, leaving 119 in the analysis (Methods). 42.5% male (Methods). Absorption, distribution and bioavailability parameters were fixed to the adult model fitted to two dense-sampling phase I studies in healthy adults (see Kovalenko_2025_alirocumab_adult)."
  )

  ini({
    # ---- Absorption and distribution: fixed to the adult estimates (Table 2) ----
    lktr <- fixed(log(22.5)); label("Transit rate constant depot -> transit1, ktr (1/day)") # Table 2 'ktr' pediatric 22.5 (fixed to adult)
    lka <- fixed(log(0.333)); label("Absorption rate constant transit1 -> central, ka (1/day)") # Table 2 'ka' pediatric 0.333 (fixed to adult)
    lfdepot <- fixed(log(0.655)); label("SC bioavailability F (fraction)") # Table 2 'F' pediatric 0.655 (fixed to adult)
    lk12 <- fixed(log(0.194)); label("Central-to-peripheral rate constant kcp (1/day)") # Table 2 'kcp' pediatric 0.194 (fixed to adult)
    lk21 <- fixed(log(0.280)); label("Peripheral-to-central rate constant kpc (1/day)") # Table 2 'kpc' pediatric 0.280 (fixed to adult)

    # ---- Estimated linear PK and PCSK9 turnover (Table 2, pediatric) ----
    lvc <- log(1.64); label("Central volume Vc at 50 kg (L)") # Table 2 'Vc (L)' pediatric 1.64
    lkel <- log(0.0800); label("Elimination rate constant of free alirocumab keL (1/day)") # Table 2 'keL (1/day)' pediatric 0.0800
    lkint <- log(0.372); label("Elimination rate constant of alirocumab:PCSK9 complexes keLR (1/day)") # Table 2 'keLR (1/day)' pediatric 0.372
    lksyn <- log(5.76); label("PCSK9 synthesis rate ksyn without high-dose statin (nM/day)") # Table 2 'ksyn (nM/day), HDS = No' pediatric 5.76
    lr0 <- log(2.43); label("Baseline free PCSK9 concentration R0 (nM)") # Table 2 'R0 (nM)' pediatric 2.43

    # ---- Binding: fixed to Biacore in-vivo values (Methods; Table S11) ----
    lkon <- fixed(log(57.9)); label("Association rate constant kon (1/nM/day)") # Table S11 'kon (1/day/nM)' 57.9 (fixed, Biacore)
    lkoff <- fixed(log(106)); label("Dissociation rate constant koff (1/day)") # Table S11 'koff (1/day)' 106 (fixed, Biacore)

    # ---- Covariate effects (Table 2) ----
    e_wt_vc <- 0.724; label("Power exponent of body weight (WT / 50 kg) on Vc (unitless)") # Table 2 'Vc ~ weight' pediatric 0.724
    e_conmed_statin_hi_ksyn <- 0.576; label("Exponential effect of high-dose statin on ksyn, exp(theta * HDS) (unitless)") # Table 2 'ksyn ~ HDS' pediatric 0.576

    # ---- Between-subject variability (Table 2 SD of log-parameter, squared) ----
    etalvc ~ 0.011664 # Table 2 'SD(ln[Vc])' pediatric 0.108 -> variance 0.108^2
    etalkint ~ 0.126736 # Table 2 'SD(ln[keLR])' pediatric 0.356 -> variance 0.356^2
    etalksyn ~ 0.379456 # Table 2 'SD(ln[ksyn])' pediatric 0.616 -> variance 0.616^2
    etalr0 ~ 0.076729 # Table 2 'SD(ln[R0])' pediatric 0.277 -> variance 0.277^2
    etalkel ~ 0.026896 # Table 2 'SD(ln[keL])' pediatric 0.164 -> variance 0.164^2

    # ---- Residual error: proportional per analyte (Table 2; Table S1 'b2', 'b5') ----
    propSd <- 0.247; label("Proportional residual error, total alirocumab (fraction)") # Table 2 'SD(alirocumab)' pediatric 0.247
    propSd_Ctarget <- 0.414; label("Proportional residual error, free PCSK9 (fraction)") # Table 2 'SD(PCSK9)' pediatric 0.414
  })

  model({
    # Alirocumab molecular weight. Kovalenko 2025 states parameters in nM and
    # concentrations in mg/L without giving the conversion; approximately
    # 146 kDa is the value on the Praluent (alirocumab) label.
    mw_alirocumab <- 146000 # g/mol
    nmol_per_mg <- 1e6 / mw_alirocumab

    # ---- Individual parameters ----
    ktr <- exp(lktr)
    ka <- exp(lka)
    fdepot <- exp(lfdepot)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    vc <- exp(lvc + etalvc) * (WT / 50)^e_wt_vc
    kel <- exp(lkel + etalkel)
    kint <- exp(lkint + etalkint)
    ksyn <- exp(lksyn + etalksyn) * exp(e_conmed_statin_hi_ksyn * CONMED_STATIN_HI)
    r0 <- exp(lr0 + etalr0)
    kon <- exp(lkon)
    koff <- exp(lkoff)
    # Methods: keR = ksyn / R0, so free PCSK9 starts at its steady state
    kdeg <- ksyn / r0

    # ---- Binding fluxes (nmol/day) in the central compartment (Figure 1) ----
    # L + R <-> L:R and L:R + R <-> L:R:R, both with kon and koff.
    bind1 <- kon * central * target / vc - koff * complex
    bind2 <- kon * complex * target / vc - koff * complex_db

    # ---- ODEs: amounts in nmol; every species shares kcp / kpc (Figure 1) ----
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ka * transit1
    d/dt(central) <- ka * transit1 - kel * central - k12 * central + k21 * peripheral1 - bind1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(target) <- ksyn * vc - kdeg * target - k12 * target + k21 * target_peripheral1 - bind1 - bind2
    d/dt(target_peripheral1) <- k12 * target - k21 * target_peripheral1
    d/dt(complex) <- bind1 - bind2 - kint * complex - k12 * complex + k21 * complex_peripheral1
    d/dt(complex_peripheral1) <- k12 * complex - k21 * complex_peripheral1
    d/dt(complex_db) <- bind2 - kint * complex_db - k12 * complex_db + k21 * complex_db_peripheral1
    d/dt(complex_db_peripheral1) <- k12 * complex_db - k21 * complex_db_peripheral1

    # Pre-dose steady state: R0 in the central compartment and the matching
    # distribution equilibrium in the peripheral compartment.
    target(0) <- r0 * vc
    target_peripheral1(0) <- r0 * vc * k12 / k21

    # Doses are in mg; SC doses enter the depot (bioavailability F), IV doses
    # the central compartment. Both are converted to nmol here, so give IV
    # infusions with `dur` rather than `rate`.
    f(depot) <- fdepot * nmol_per_mg
    f(central) <- nmol_per_mg

    # ---- Observations ----
    # Total alirocumab assay measures free plus 1- and 2-PCSK9-bound drug
    # (Analytical Methods); reported in mg/L. Free PCSK9 is reported in nM.
    Cc <- (central + complex + complex_db) / vc / nmol_per_mg
    Ctarget <- target / vc
    Cc ~ prop(propSd)
    Ctarget ~ prop(propSd_Ctarget)
  })
}
