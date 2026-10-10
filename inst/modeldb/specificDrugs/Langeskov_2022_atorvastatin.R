Langeskov_2022_atorvastatin <- function() {
  description <- "Two-compartment population PK model with one transit compartment followed by first-order absorption for a single 40 mg oral dose of atorvastatin in healthy adults with and without steady-state once-weekly subcutaneous semaglutide 1.0 mg (Langeskov 2022). Semaglutide co-administration reduces ka and the transit rate constant ktr by fractions of 0.829 and 0.791; fixed allometric body-weight scaling (exponent 0.75 on clearances, 1 on volumes, 75 kg reference) and proportional residual error."
  reference <- paste(
    "Langeskov EK, Kristensen K. Population pharmacokinetic of paracetamol",
    "and atorvastatin with co-administration of semaglutide.",
    "Pharmacol Res Perspect. 2022;10(4):e00962. doi:10.1002/prp2.962"
  )
  vignette <- "Langeskov_2022_semaglutide_ddi"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot = list(analyte = "atorvastatin", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "atorvastatin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "atorvastatin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "atorvastatin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of CL/F, Cl2/F (exponent 0.75) and V1/F, V2/F (exponent 1) on (WT / MedianBW), Langeskov 2022 Section 3.2 PML equations. The Table 3 footnote states the disposition parameters are 'For an individual of 75 kg', so the reference weight is 75 kg (cohort mean 75.4 kg, Table 1).",
      source_name = "BW"
    ),
    CONMED_SEMAGLUTIDE = list(
      description = "Semaglutide co-administration indicator: 1 = atorvastatin dosed at semaglutide steady state (once-weekly SC semaglutide escalated 0.25 -> 0.5 -> 1.0 mg, dosed in week 13), 0 = atorvastatin dosed before semaglutide was started.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (atorvastatin alone)",
      notes = "Per-occasion indicator in a one-sequence crossover. Coded like the paracetamol analysis (data-set column `placebo` = 1 for semaglutide co-administration, Section 3.2), so CONMED_SEMAGLUTIDE equals the source column. Enters ka and ktr as tv * (1 - theta * CONMED_SEMAGLUTIDE).",
      source_name = "placebo"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 1L,
    n_observations = 713L,
    age_range = "25-55 years",
    age_mean = "45 years",
    weight_range = "53.6-102 kg",
    weight_mean = "75.4 kg",
    bmi_range = "20.4-29.8 kg/m^2 (mean 25.2)",
    sex_female_pct = 51.6,
    race_ethnicity = "Not reported in source paper.",
    disease_state = "Healthy adults (BMI 20-30 kg/m^2).",
    dose_range = "Single oral 40 mg atorvastatin before semaglutide dose escalation and again at the end of week 13 at semaglutide steady state (once-weekly SC 0.25 mg x 4 weeks, 0.5 mg x 4 weeks, 1.0 mg thereafter).",
    regions = "Single centre (trial NCT02243098, Hausner 2017).",
    notes = "Open-label one-sequence crossover; 31 subjects received atorvastatin alone and 26 with semaglutide (five withdrew). Samples at 0.5, 1, 2, 3, 4, 6, 8, 10, 12, 18, 24, 36 and 48 h post-dose; pre-dose and 72 h samples were all below LLOQ and excluded. Demographics from Langeskov 2022 Table 1 (15 male / 16 female). Fitted in Phoenix NLME 8.1 with FOCE-ELS. Only the parent acid was modelled; the lactone metabolites were not."
  )

  ini({
    # Final-model estimates, Langeskov 2022 Table 3 (atorvastatin).
    lka <- log(5.72); label("Absorption rate constant without semaglutide (1/h)") # Table 3 'ka (h-1) (Placebo)' = 5.72
    e_conmed_semaglutide_ka <- 0.829; label("Fractional reduction of ka with semaglutide co-administration (unitless)") # Table 3 'Theta kacovariate' = 0.829; Section 3.2 dKadplacebo
    lktr <- log(7.28); label("Transit rate constant without semaglutide (1/h)") # Table 3 'ktr (h-1) (Placebo)' = 7.28
    e_conmed_semaglutide_ktr <- 0.791; label("Fractional reduction of ktr with semaglutide co-administration (unitless)") # Table 3 'Theta ktrcovariate' = 0.791; Section 3.2 dKtrdplacebo
    lvc <- log(1843); label("Apparent central volume V1/F for a 75 kg individual (L)") # Table 3 'V1/F (L)' = 1843 (footnote a: 75 kg)
    lcl <- log(620); label("Apparent oral clearance CL/F for a 75 kg individual (L/h)") # Table 3 'Cl/F (L/h)' = 620 (footnote a: 75 kg)
    lvp <- log(4184); label("Apparent peripheral volume V2/F for a 75 kg individual (L)") # Table 3 'V2/F (L)' = 4184 (footnote a: 75 kg)
    lq <- log(873); label("Apparent intercompartmental clearance Cl2/F for a 75 kg individual (L/h)") # Table 3 'Cl2/F (L/h)' = 873 (footnote a: 75 kg)
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on CL/F and Cl2/F (unitless)") # Section 3.2 PML equation for CL/F, exponent 0.75; applied to Cl2/F per the same section
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on V1/F and V2/F (unitless)") # Section 3.2 PML equation for V1/F, exponent 1; applied to V2/F per the same section

    # Between-subject variability: Table 3 reports omega^2 (variances of
    # exponential etas, Section 2.1). No BSV on Cl2/F (Section 3.1).
    etalka ~ 2.03 # Table 3 'omega2 Ka' = 2.03
    etalktr ~ 1.82 # Table 3 'omega2 Ktr' = 1.82
    etalvc ~ 0.625 # Table 3 'omega2 V1/F' = 0.625
    etalcl ~ 0.151 # Table 3 'omega2 Cl/f' = 0.151
    etalvp ~ 0.111 # Table 3 'omega2 V2/F' = 0.111

    # Proportional residual error. Table 3 prints 'Ceps' = 0.322; read as
    # the Phoenix CEps standard deviation, the same scale established for
    # the paracetamol model from its Figure 3A scatter.
    propSd <- 0.322; label("Proportional residual error (fraction)") # Table 3 'Residual unexplained variability (Ceps)' = 0.322
  })

  model({
    # Individual parameters
    ka <- exp(lka + etalka) * (1 - e_conmed_semaglutide_ka * CONMED_SEMAGLUTIDE)
    ktr <- exp(lktr + etalktr) * (1 - e_conmed_semaglutide_ktr * CONMED_SEMAGLUTIDE)
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl_q
    q <- exp(lq) * (WT / 75)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc_vp
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Figure 2: dose -> transit compartment -(ktr)-> absorption compartment
    # -(ka)-> central. The dose enters `depot` (the paper's transit
    # compartment); `transit1` is the paper's absorption compartment.
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ka * transit1
    d/dt(central) <- ka * transit1 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L; x 1000 gives ug/L (Figures 4 and 6).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
