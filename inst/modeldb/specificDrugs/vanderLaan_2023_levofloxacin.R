vanderLaan_2023_levofloxacin <- function() {
  description <- paste(
    "Two-compartment population PK model for oral levofloxacin in young",
    "South African children (aged 0.16 to 6 years) receiving levofloxacin",
    "preventive therapy after household exposure to rifampicin-resistant",
    "tuberculosis (van der Laan 2023; n = 25; randomised crossover of the",
    "100 mg dispersible tablet against the crushed 250 mg non-dispersible",
    "adult tablet). Savic transit-compartment absorption (N fixed to 10,",
    "MTT fixed to 0.167 h, written as an explicit compartment chain)",
    "feeds first-order absorption into central.",
    "CL and Q scale allometrically with weight (exponent 0.75) and Vc and",
    "Vp with exponent 1, all referenced to the cohort-median 12.2 kg, and",
    "CL matures with postmenstrual age through a Hill function (TM50 = 12.7",
    "months, Hill = 3.95). The crushed non-dispersible tablet has 21.5%",
    "lower bioavailability than the dispersible tablet (F fixed to 1).",
    "Between-subject variability is on CL only; between-occasion",
    "variability is on ka, MTT and F, with the BOV on F inflated 3.24-fold",
    "on the occasion of the routine (non-study) dose that precedes the",
    "first pre-dose sample.",
    sep = " "
  )
  reference <- paste(
    "van der Laan LE, Hesseling AC, Schaaf HS, Palmer M, Draper HR,",
    "Wiesner L, Denti P, Garcia-Prats AJ.",
    "Pharmacokinetics and optimized dosing of dispersible and",
    "non-dispersible levofloxacin formulations in young children.",
    "J Antimicrob Chemother. 2023;78(10):2481.",
    "doi:10.1093/jac/dkad257.",
    sep = " "
  )
  vignette <- "vanderLaan_2023_levofloxacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    transit7 = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    transit8 = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    transit9 = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    transit10 = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    transit11 = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "levofloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "levofloxacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of CL and Q (exponent 0.75) and of Vc and Vp",
        "(exponent 1), referenced to the cohort-median 12.2 kg (Table 2",
        "footnote 'Estimates reported for a child with the median weight",
        "of 12.2 kg'). The paper states that allometric scaling was",
        "applied to all clearance and volume parameters (Methods, citing",
        "Anderson and Holford) without printing the exponents; the",
        "standard 0.75 / 1 values are those of the parent model",
        "(Denti 2018) that this analysis started from."
      ),
      source_name = "WT"
    ),
    AGE = list(
      description = "Postnatal age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the postmenstrual-age maturation of CL. Postmenstrual age",
        "in months is computed inside model() as AGE * 12 + 9, i.e.",
        "assuming term gestation of 9 months; gestational age was not",
        "reported. The 9-month offset is the one the parent model",
        "(Denti 2018) states, and it agrees with the Discussion's",
        "'half its mature value at ~4 months (~13 months post-menstrual",
        "age)' for TM50 = 12.7 months."
      ),
      source_name = "AGE"
    ),
    FORM_LEVOFLOXACIN_NONDISPERSIBLE = list(
      description = paste(
        "Levofloxacin formulation indicator: 1 = crushed 250 mg",
        "non-dispersible adult tablet, 0 = 100 mg dispersible paediatric",
        "tablet"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (100 mg dispersible tablet, F fixed to 1)",
      notes = paste(
        "Per-dose record. Bioavailability is multiplied by",
        "(1 + e_form_levofloxacin_nondispersible_fdepot *",
        "FORM_LEVOFLOXACIN_NONDISPERSIBLE) with the Table 2 effect",
        "-21.5%. Both formulations (Macleods) were mixed in 2.5-5 mL of",
        "water and followed by the same volume of rinse; the",
        "non-dispersible tablet was crushed first (Methods). No",
        "formulation effect on ka, MTT or on the variability was",
        "retained (Results)."
      ),
      source_name = "formulation"
    ),
    OCC = list(
      description = paste(
        "Dosing occasion: 1 = routine (non-study) levofloxacin dose taken",
        "the day before a pharmacokinetic day, observed through its",
        "pre-dose sample; 2 = study dose on pharmacokinetic day 1;",
        "3 = study dose on pharmacokinetic day 2"
      ),
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Selects the between-occasion etas on ka, MTT and F. Methods: 'The",
        "pharmacokinetic samples collected pre-dose were treated as a",
        "separate occasion in the model'. Occasion 1 carries the BOV on F",
        "scaled by 3.24 (Table 2 'Scaling of between-occasion variability",
        "in bioavailability for pre-doses'). OCC = 0 (or any value outside",
        "1-3) zeroes every between-occasion eta, which is how the",
        "typical-value and steady-state simulations in the vignette are run.",
        "The three-occasion numbering is a construct of this encoding; the",
        "paper defines the occasions but does not print their indices."
      ),
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 25L,
    n_studies = 1L,
    age_range = "0.16 to 6.01 years",
    age_median = "2.56 years (IQR 1.58-4.03)",
    weight_range = "4.02 to 20.3 kg",
    weight_median = "12.2 kg (IQR 10.7-15.0)",
    sex_female_pct = 40,
    race_ethnicity = c(MixedRace = 68, Black = 32),
    disease_state = paste(
      "HIV-negative children under 6 years receiving levofloxacin",
      "preventive therapy (for at least 2 weeks) after household exposure",
      "to an adult with rifampicin-resistant tuberculosis; all had normal",
      "renal function. Routine co-medication: high-dose isoniazid",
      "(19 / 25) and ethambutol (8 / 25)."
    ),
    dose_range = paste(
      "Single oral doses of about 15-20 mg/kg on each of two",
      "pharmacokinetic days, one with the 100 mg dispersible tablet and one",
      "with the crushed 250 mg non-dispersible tablet (Macleods), in",
      "randomised order; routine levofloxacin (Austell, 23 / 25) on other",
      "days."
    ),
    sampling = paste(
      "Pre-dose and 1, 2, 4, 6 and 8 h post-dose on both pharmacokinetic",
      "days; 286 samples, LLOQ 0.0781 mg/L, BLQ handled by Beal M3."
    ),
    regions = "South Africa (Cape Town).",
    notes = paste(
      "Demographics from van der Laan 2023 Table 1 (median creatinine",
      "24 umol/L, median fat-free mass 9.89 kg). Age strata enrolled:",
      "0 to <1, 1 to <2 and 2 to <6 years."
    )
  )

  ini({
    # Structural parameters, reported for a child of the median weight
    # 12.2 kg (Table 2 footnote). CL is the fully mature value: the
    # maturation factor below is not renormalised to a reference age.
    lcl <- log(4.17); label("Clearance CL/F at 12.2 kg, fully mature (L/h)") # Table 2 'Clearance (L/h)' = 4.17 (3.69, 4.64)
    lvc <- log(19.7); label("Central volume Vc/F at 12.2 kg (L)") # Table 2 'Volume of distribution (L)' = 19.7 (17.7, 21.9)
    lq <- log(0.638); label("Intercompartmental clearance Q/F at 12.2 kg (L/h)") # Table 2 'Intercompartmental clearance (L/h)' = 0.638 (0.365, 0.951)
    lvp <- log(5.53); label("Peripheral volume Vp/F at 12.2 kg (L)") # Table 2 'Peripheral volume of distribution (L)' = 5.53 (3.55, 8.38)
    lka <- log(2.52); label("Absorption rate constant ka (1/h)") # Table 2 'Absorption rate constant' = 2.52 (2.07, 2.90)
    lmtt <- fixed(log(0.167)); label("Mean transit time of the Savic transit-absorption chain (h)") # Table 2 'Mean transit time (h)' = 0.167, fixed
    lfdepot <- fixed(log(1)); label("Bioavailability of the dispersible tablet (fraction)") # Table 2 'Bioavailability for dispersible tablet (%)' = 100, fixed (reference)

    # Allometric exponents (Methods; values from the parent Denti 2018 model)
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)") # Methods: allometric scaling on all clearances (Anderson and Holford)
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on Vc and Vp (unitless)") # Methods: allometric scaling on all volumes (Anderson and Holford)

    # Postmenstrual-age maturation of CL (Hill function)
    ltm50_cl <- log(12.7); label("Postmenstrual age at 50% maturation of CL (months)") # Table 2 'Postmenstrual age at which 50% of maturation is reached' = 12.7 (11.2, 14.1)
    lhill_cl <- log(3.95); label("Hill coefficient of the CL maturation function (unitless)") # Table 2 'Shape factor of the curve' = 3.95 (2.74, 4.79)

    # Formulation effect on bioavailability
    e_form_levofloxacin_nondispersible_fdepot <- -0.215; label("Fractional change in F for the crushed non-dispersible tablet (unitless)") # Table 2 'Effect of crushed non-dispersible tablet on bioavailability (%)' = -21.5 (-31.4 to -11.7)

    # Ratio of the occasion-1 (routine dose) to the study-occasion SD of the
    # bioavailability BOV; multiplies the eta, so the variance scales by 3.24^2.
    sd_ratio_fdepot_occ1 <- 3.24; label("Scaling of the BOV SD in F for the routine-dose occasion (unitless)") # Table 2 'Scaling of between-occasion variability in bioavailability for pre-doses' = 3.24-fold (1.98, 5.40)

    # Variability. Table 2 prints each variability as a percentage equal to
    # sqrt(omega^2) * 100 (the reporting convention of this group's
    # analyses; see Sawe 2025 and Tsirizani 2025, whose deposited control
    # streams confirm it), so omega^2 = (pct / 100)^2.
    etalcl ~ 0.008649 # Table 2 'Between-subject variability (%) for clearance' 9.30 -> 0.093^2

    # Between-occasion variability: one shared magnitude per parameter
    # across the three occasions ($OMEGA BLOCK(1) SAME).
    etaiov_ka_1 ~ 0.293764 # Table 2 BOV 'Absorption rate constant' 54.2% -> 0.542^2
    etaiov_ka_2 ~ fixed(0.293764) # occasion 2; same magnitude
    etaiov_ka_3 ~ fixed(0.293764) # occasion 3; same magnitude
    etaiov_mtt_1 ~ 0.964324 # Table 2 BOV 'Mean transit time' 98.2% -> 0.982^2
    etaiov_mtt_2 ~ fixed(0.964324) # occasion 2; same magnitude
    etaiov_mtt_3 ~ fixed(0.964324) # occasion 3; same magnitude
    etaiov_fdepot_1 ~ 0.051076 # Table 2 BOV 'Bioavailability' 22.6% -> 0.226^2 (scaled by sd_ratio_fdepot_occ1 in model())
    etaiov_fdepot_2 ~ fixed(0.051076) # occasion 2; same magnitude
    etaiov_fdepot_3 ~ fixed(0.051076) # occasion 3; same magnitude

    # Residual error: combined additive and proportional (Methods)
    addSd <- 0.112; label("Additive residual error (mg/L)") # Table 2 'Additive error (mg/L)' = 0.112 (0.0850, 0.161)
    propSd <- 0.0533; label("Proportional residual error (fraction)") # Table 2 'Proportional error (%)' = 5.33 (4.08, 6.70)
  })

  model({
    # 1. Occasion indicators and between-occasion etas. OCC = 0 zeroes them.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3
    iov_fdepot <- oc1 * sd_ratio_fdepot_occ1 * etaiov_fdepot_1 +
      oc2 * etaiov_fdepot_2 + oc3 * etaiov_fdepot_3

    # 2. Maturation of CL with postmenstrual age (months), PMA = postnatal
    # age + 9 months of assumed term gestation.
    pma <- AGE * 12 + 9
    tm50_cl <- exp(ltm50_cl)
    hill_cl <- exp(lhill_cl)
    mat_cl <- pma^hill_cl / (pma^hill_cl + tm50_cl^hill_cl)

    # 3. Individual parameters, allometric reference 12.2 kg
    cl <- exp(lcl + etalcl) * (WT / 12.2)^e_wt_cl_q * mat_cl
    vc <- exp(lvc) * (WT / 12.2)^e_wt_vc_vp
    q <- exp(lq) * (WT / 12.2)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 12.2)^e_wt_vc_vp
    ka <- exp(lka + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    fbio <- exp(lfdepot + iov_fdepot) *
      (1 + e_form_levofloxacin_nondispersible_fdepot * FORM_LEVOFLOXACIN_NONDISPERSIBLE)

    # N = 10 transit compartments, fixed (Table 2 'No. of transit compartments'
    # = 10, fixed); the chain below has exactly that many.
    ktr <- (10 + 1) / mtt
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODEs. Savic transit absorption with N = 10 written as an explicit
    # chain: the dose compartment (depot) and the N transit compartments all
    # empty at ktr = (N + 1) / MTT, so the input to the absorption
    # compartment (transit11, emptied by ka) is exactly Savic's
    # F * Dose * ktr * (ktr * t)^N * exp(-ktr * t) / N!. The explicit chain
    # superposes doses correctly and stays accurate when the occasion eta
    # makes MTT very short, where the analytical transit() input is a spike
    # the solver can step over.
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ktr * transit6
    d/dt(transit7) <- ktr * transit6 - ktr * transit7
    d/dt(transit8) <- ktr * transit7 - ktr * transit8
    d/dt(transit9) <- ktr * transit8 - ktr * transit9
    d/dt(transit10) <- ktr * transit9 - ktr * transit10
    d/dt(transit11) <- ktr * transit10 - ka * transit11
    d/dt(central) <- ka * transit11 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fbio

    # 5. Observation: dose in mg and V in L give mg/L
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
