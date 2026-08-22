Mo_2024_shikimicAcid_pig_oral <- function() {
  description <- paste(
    "Preclinical (pig). One-compartment first-order-absorption pharmacokinetic",
    "model for shikimic acid (SA) in growing Landrace x Large White pigs after",
    "a single 50 mg/kg intragastric (gavage) dose (Mo 2024). Mo 2024 selected",
    "this structure over a two-compartment and over lag-time alternatives by",
    "AIC (Section 2.6, supplementary Table S2) and describes the profile as C =",
    "M*(exp(-ke*t) - exp(-ka*t)); the absorption and elimination rate constants",
    "carried here are back-calculated from the reported absorption and",
    "elimination half-lives, and the volume from the reported apparent",
    "clearance. Because only the extravascular route is observed, clearance and",
    "volume are apparent (Cl/F and V/F); the absolute bioavailability of SA in",
    "these pigs was 21.68 percent. All states are expressed per kg body weight",
    "(volume in mL/kg, clearance in mL/h/kg, amounts in ng/kg). Two of the six",
    "pigs showed a double peak that this single-depot model does not reproduce.",
    "Mo 2024 fitted each pig individually in Phoenix WinNonlin and reported only",
    "the mean and SD of the individual estimates, so no between-subject",
    "variability or residual-error model is available; every parameter is fixed",
    "at the published mean and the residual SD is fixed at zero. No",
    "pharmacodynamic model accompanies this route - Mo 2024 fitted the",
    "sigmoid-Emax immune-response models to the intravenous group only; see",
    "modellib('Mo_2024_shikimicAcid_pig_iv')."
  )
  reference <- paste(
    "Mo K, Shen Y, Su D, Lv L, Du J, Ding H, Huang X (2024).",
    "Pharmacokinetic-Pharmacodynamic Modeling of the Immune-Enhancing Effect of",
    "Shikimic Acid in Growing Pigs. J Agric Food Chem 72:26224-26235.",
    "doi:10.1021/acs.jafc.4c09250"
  )
  vignette <- "Mo_2024_shikimicAcid"
  units <- list(time = "h", dosing = "ng/kg", concentration = "ng/mL")

  compartmentData <- list(
    depot   = list(analyte = "shikimic acid", units = "ng", specimen = "administration site", verified = TRUE),
    central = list(analyte = "shikimic acid", units = "ng", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species       = "pig (Landrace x Large White)",
    n_subjects    = 6,
    n_studies     = 1,
    age_median    = "approximately 90 days",
    sex           = "male and female",
    disease_state = "healthy growing pigs",
    dose_range    = "shikimic acid 50 mg/kg as a single intragastric (gavage) dose",
    regions       = "China",
    notes         = paste(
      "Mo 2024 Section 2.2 and Table 1: six growing pigs of comparable health",
      "status and genetic background, approximately 90 days old, allocated to a",
      "two-period two-sequence crossover (sequence A n = 3, sequence B n = 3)",
      "with a 7-day washout; every pig therefore contributes both an",
      "intragastric and an intravenous profile. The same six pigs support the",
      "companion model Mo_2024_shikimicAcid_pig_iv. The gavage dose was given as",
      "a 50 mg/mL aqueous solution. Body weights are given in supplementary",
      "Table S1, which is not part of the open-access record; the paper states",
      "only that weight did not differ between periods (P = 0.320). Plasma was",
      "sampled predose and at 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12 and",
      "24 h (Section 2.3). Both sexes were studied and the authors note sex may",
      "affect SA disposition, but no sex covariate was fitted. No race/ethnicity",
      "data apply."
    )
  )

  ini({
    # One-compartment first-order absorption. Mo 2024 Table 9 reports, as the
    # mean of the six individual intragastric fits, an absorption half-life
    # t1/2ka = 0.85 h, an elimination half-life t1/2ke = 1.81 h and an apparent
    # clearance Cl_F = 1086.52 mL/h/kg for the 50 mg/kg (5e7 ng/kg) dose. The
    # rate constants follow as ka = ln(2)/0.85 = 0.8155 /h and ke = ln(2)/1.81 =
    # 0.3830 /h, and the apparent volume as V/F = Cl_F/ke = 2837.21 mL/kg.
    #
    # These reproduce the observations Mo 2024 tabulates independently: the
    # analytic Tmax = ln(ka/ke)/(ka - ke) is 1.75 h against the reported
    # 1.78 +/- 0.66 h (Table 9), and Dose/Cl_F is 46,018 h*ng/mL against the
    # reported AUC of 50,084 +/- 16,217 h*ng/mL. The Tmax agreement also settles
    # the sign in the printed absorption equation: Section 3.2 prints
    # C = M*(exp(-ke*t) + exp(-ka*t)), a sum that decays monotonically from t = 0
    # and admits no Tmax at all; the Bateman difference form
    # C = M*(exp(-ke*t) - exp(-ka*t)) is the model actually fitted, and it is
    # what the parameterisation below encodes.
    lka <- fixed(log(0.815467)); label("First-order absorption rate constant (1/h)")  # ln(2)/0.85 h; Table 9 t1/2ka
    lcl <- fixed(log(1086.52)); label("Apparent clearance Cl/F (mL/h/kg)")  # Table 9, Cl_F
    lvc <- fixed(log(2837.21)); label("Apparent central volume of distribution V/F (mL/kg)")  # Cl_F / (ln(2)/1.81 h); Table 9 Cl_F and t1/2ke

    addSd <- fixed(0); label("Additive residual SD on shikimic acid plasma concentration (ng/mL; not reported)")  # Mo 2024 reports no residual-error model
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc

    Cc ~ add(addSd)
  })
}
