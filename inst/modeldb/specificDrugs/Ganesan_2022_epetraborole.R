Ganesan_2022_epetraborole <- function() {
  description <- paste(
    "Three-compartment population PK model for epetraborole after intravenous",
    "and oral dosing in healthy adults and adults with complicated urinary",
    "tract or intra-abdominal infection, pooled from five phase 1 and two",
    "phase 2 studies (138 subjects, 2637 plasma samples; Ganesan 2022).",
    "Linear elimination. Body weight scales all clearances and volumes",
    "allometrically against a 70 kg reference. Oral doses pass through a",
    "depot and two transit compartments at a shared rate equal to the",
    "absorption rate constant, which takes a separate value in the fed and",
    "fasted states; absolute bioavailability is logit-normal. An epithelial",
    "lining fluid (ELF) compartment receives drug from the central compartment",
    "(K14) and loses it to a sink (K40); ELF concentration is the ELF amount",
    "divided by the central volume. IIV on CL, Vc, Vp, Vp2, ka, F and K14;",
    "two-occasion IOV on ka; proportional residual error for plasma and ELF."
  )
  reference <- paste(
    "Ganesan H, Safir MC, Bhavnani SM, Krause KM, Rubino CM. 593. Population",
    "Pharmacokinetic Model Development for Epetraborole and Mycobacterium",
    "avium Complex (MAC) Lung Disease Patients Using Data from Phase 1 and 2",
    "Studies. Open Forum Infect Dis. 2022;9(Suppl 2):S325-S326.",
    "doi:10.1093/ofid/ofac492.645 (IDWeek 2022 poster 593).",
    "Parameter estimates: abstract Table 2 (identical in the poster).",
    sep = " "
  )
  vignette <- "Ganesan_2022_epetraborole"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Doses in mg and volumes in L, so concentrations are mg/L (= ug/mL); the
  # source figures plot ng/mL. Oral doses go to "depot" and intravenous doses
  # to "central".
  dosing <- c("depot", "central")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 2 reports CL, CLd and CLd2 in 'L/h/70 kg' and Vc, Vp and Vp2",
        "in 'L/70 kg'; the Results state that body weight acts through 'an",
        "allometric scaling approach'. The exponents are not printed; the",
        "standard fixed values 0.75 (clearances) and 1 (volumes) are used."
      ),
      source_name = "WT"
    ),
    FED = list(
      description = "Fed state at the oral dose (1 = fed, 0 = fasted)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Selects between the two independently estimated absorption rate",
        "constants of Table 2 ('Ka, fast' = 2.89 1/h, 'Ka, fed' = 1.34 1/h).",
        "Bioavailability does not depend on food in the model. Irrelevant for",
        "intravenous doses."
      ),
      source_name = "fed / fasted administration"
    ),
    OCC = list(
      description = "Occasion index for the inter-occasion variability on ka",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Table 2 reports one IOV variance for ka (0.0343) with a separate",
        "shrinkage for occasion 1 and occasion 2, so two occasions share one",
        "variance (NONMEM BLOCK(1) SAME). The source does not define an",
        "occasion; OCC = 1 for the first oral dosing occasion and 2 for a",
        "later one is the natural reading. OCC values other than 1 or 2 give",
        "no IOV contribution."
      ),
      source_name = "occasion"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "epetraborole", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "epetraborole", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "epetraborole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "epetraborole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "epetraborole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "epetraborole", units = "mg", specimen = "plasma", verified = TRUE),
    elf = list(analyte = "epetraborole", units = "mg", specimen = "epithelial lining fluid", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 138,
    n_studies = 7,
    n_observations = "2637 epetraborole PK samples (plasma, plus epithelial lining fluid from the phase 1 ELF study)",
    disease_state = paste(
      "Healthy adults (five phase 1 studies) and adults with complicated",
      "urinary tract infection / acute pyelonephritis or complicated",
      "intra-abdominal infection (two phase 2 studies). No patients with",
      "Mycobacterium avium complex lung disease contributed data; the model",
      "was built to simulate that population."
    ),
    dose_range = paste(
      "IV 200-3000 mg single dose and 500-1600 mg BID (first-in-human),",
      "1500 mg single dose and BID (ELF study), 1500 mg (mass balance),",
      "750 or 1500 mg BID (phase 2); oral 500-4000 mg single dose and",
      "2000-3000 mg BID for 10 days (food-effect study) and 250-750 mg QD or",
      "500 mg QOD tablets for 28 days (phase 1b, NCT04892641)."
    ),
    regions = "Not reported (phase 1b conducted in Australia)",
    notes = paste(
      "Abstract Table 1 lists the seven studies. Demographics (age, weight,",
      "sex, race) are not reported. Only the first four cohorts of the",
      "phase 1b study were included."
    )
  )

  ini({
    # Structural disposition parameters, typical values for a 70 kg adult.
    lcl <- log(15.3); label("Systemic clearance CL for 70 kg (L/h)") # Table 2 'CL (L/h/70 kg)' = 15.3 (%SEM 4.91)
    lq <- log(23); label("Distributional clearance to peripheral compartment 1, CLd, for 70 kg (L/h)") # Table 2 'Cld (L/h/70 kg)' = 23 (%SEM 7.18)
    lq2 <- log(43.3); label("Distributional clearance to peripheral compartment 2, CLd2, for 70 kg (L/h)") # Table 2 'Cld2 (L/h/70 kg)' = 43.3 (%SEM 6.68)
    lvc <- log(15.6); label("Central volume Vc for 70 kg (L)") # Table 2 'Vc (L/70 kg)' = 15.6 (%SEM 6)
    lvp <- log(140); label("Peripheral volume Vp for 70 kg (L)") # Table 2 'Vp (L/70 kg)' = 140 (%SEM 3.02)
    lvp2 <- log(33.2); label("Second peripheral volume Vp2 for 70 kg (L)") # Table 2 'Vp2 (L/70 kg)' = 33.2 (%SEM 9.62)

    # Allometric exponents: not printed; standard fixed values (see covariateData$WT).
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL, CLd and CLd2 (unitless)") # Results 'allometric scaling approach'; exponent not reported
    e_wt_vc <- fixed(1); label("Allometric exponent on Vc, Vp and Vp2 (unitless)") # Results 'allometric scaling approach'; exponent not reported

    # Oral absorption: depot -> transit1 -> transit2 -> central, every step at ka.
    lka_fasted <- log(2.89); label("Absorption rate constant ka in the fasted state (1/h, FED = 0)") # Table 2 'Ka, fast (1/h)' = 2.89 (%SEM 5.33)
    lka_fed <- log(1.34); label("Absorption rate constant ka in the fed state (1/h, FED = 1)") # Table 2 'Ka, fed (1/h)' = 1.34 (%SEM 9.18)
    logitfdepot <- log(0.565 / (1 - 0.565)); label("Logit of absolute oral bioavailability F (unitless)") # Table 2 'F' = 0.565 (%SEM 3.07)

    # Epithelial lining fluid compartment (compartment 4 in the source numbering).
    lk_central_elf <- log(0.266); label("Rate constant from central to ELF compartment, K14 (1/h)") # Table 2 'K14 (1/h)' = 0.266 (%SEM 17.2)
    lkeff_elf <- log(0.661); label("Rate constant out of the ELF compartment, K40 (1/h)") # Table 2 'K40 (1/h)' = 0.661 (%SEM 14.7)

    # IIV. Table 2 'Magnitude of IIV' is sqrt(omega^2) x 100 (the IOV row
    # pairs 0.0343 with 18.5), so each variance is (magnitude / 100)^2.
    etalcl ~ 0.006241 # Table 2 CL IIV 7.9 -> 0.079^2
    etalvc ~ 0.138384 # Table 2 Vc IIV 37.2 -> 0.372^2
    etalvp ~ 0.00753424 # Table 2 Vp IIV 8.68 -> 0.0868^2
    etalvp2 ~ 0.035721 # Table 2 Vp2 IIV 18.9 -> 0.189^2
    etalka ~ 0.087025 # Table 2 'Ka, fast' IIV 29.5 -> 0.295^2; shared by the fed ka
    etalogitfdepot ~ 0.106929 # Table 2 F IIV 32.7 -> 0.327^2, logit scale
    etalk_central_elf ~ 0.334084 # Table 2 K14 IIV 57.80 -> 0.578^2

    # Inter-occasion variability on ka, two occasions sharing one variance.
    etaiov_ka_1 ~ 0.0343 # Table 2 'IOV in Ka' = 0.0343 (%SEM 26.1; 18.5)
    etaiov_ka_2 ~ fix(0.0343) # Table 2 'IOV in Ka' shrinkage reported for occasion 2; same variance

    # Residual error. Table 2 prints the variance; SD = sqrt(variance).
    propSd <- 0.242281; label("Plasma proportional residual error (fraction)") # Table 2 'RVplasma' = 0.0587 (24.2%) -> sqrt(0.0587)
    propSd_Celf <- 0.256125; label("ELF proportional residual error (fraction)") # Table 2 'RVELF' = 0.0656 (25.6%) -> sqrt(0.0656)
  })

  model({
    # Individual disposition parameters with allometric scaling to 70 kg.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    q <- exp(lq) * (WT / 70)^e_wt_cl
    q2 <- exp(lq2) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc
    vp2 <- exp(lvp2 + etalvp2) * (WT / 70)^e_wt_vc

    # Absorption: the fed or fasted ka, with IIV and two-occasion IOV.
    iov_ka <- etaiov_ka_1 * (OCC == 1) + etaiov_ka_2 * (OCC == 2)
    ka <- exp((1 - FED) * lka_fasted + FED * lka_fed + etalka + iov_ka)
    fdepot <- expit(logitfdepot + etalogitfdepot)

    k_central_elf <- exp(lk_central_elf + etalk_central_elf)
    keff_elf <- exp(lkeff_elf)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(transit2) <- ka * transit1 - ka * transit2
    d/dt(central) <- ka * transit2 - (kel + k12 + k13 + k_central_elf) * central +
      k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    # K14 moves drug out of central and K40 removes it from the ELF to a sink
    # (NONMEM general-linear microconstant notation; compartment 0 is outside).
    d/dt(elf) <- k_central_elf * central - keff_elf * elf

    f(depot) <- fdepot

    Cc <- central / vc
    Celf <- elf / vc

    Cc ~ prop(propSd)
    Celf ~ prop(propSd_Celf)
  })
}
