Goulooze_2020_DRL17822 <- function() {
  description <- "Two-compartment population PK model with a six-transit-compartment absorption chain for the oral cholesteryl ester transfer protein (CETP) inhibitor DRL-17822 in healthy male volunteers, with the food-by-formulation interaction on relative bioavailability binned into four levels (nanocrystal vs amorphous solid dispersion; fasted, low-fat, high-fat or continental breakfast), a slower absorption rate for the amorphous solid dispersion after a high-fat breakfast, BMI on central volume, and correlated inter-individual and inter-occasion variability on relative bioavailability and absorption rate"
  reference <- "Goulooze SC, Kruithof AC, Alikunju S, Gautam A, Burggraaf J, Kamerling IMC, Stevens J. The effect of food and formulation on the population pharmacokinetics of cholesteryl ester transferase protein inhibitor DRL-17822 in healthy male volunteers. Br J Clin Pharmacol. 2020;86(10):2095-2101. doi:10.1111/bcp.14297"
  vignette <- "Goulooze_2020_DRL17822"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on central volume centred at 23.3 kg/m^2 (supplement control stream: V2 = THETA(4) * EXP(ETA(1)) * (1 + (BMI - 23.3) * THETA(11))). The analysis covered BMI 18.8-29.9 kg/m^2 (Table S1).",
      source_name = "BMI"
    ),
    FED = list(
      description = "Fed-versus-fasted state at the dose record",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (overnight fast)",
      notes = "1 for any breakfast before the dose (low-fat, high-fat or continental). Together with FED_HIGHFAT, FED_LOWFAT and FORM_DRL17822_ASD it selects the relative-bioavailability bin. A continental breakfast is FED = 1 with FED_HIGHFAT = FED_LOWFAT = 0. Time-varying per dose record in the crossover studies.",
      source_name = "FEED (FEED = 0 fasted; any non-zero value fed)"
    ),
    FED_HIGHFAT = list(
      description = "High-fat breakfast before the dose",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted, low-fat or continental breakfast)",
      notes = "Control-stream FEED = 1. For the nanocrystal it gives the reference bioavailability (same bin as continental); for the amorphous solid dispersion it gives Fmedium and also scales ka by 0.599. Mutually exclusive with FED_LOWFAT; implies FED = 1.",
      source_name = "FEED == 1"
    ),
    FED_LOWFAT = list(
      description = "Low-fat breakfast before the dose",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted, high-fat or continental breakfast)",
      notes = "Control-stream FEED = 3. Gives Fmedium for both formulations (the control stream sets FM = FLF whatever the formulation). Mutually exclusive with FED_HIGHFAT; implies FED = 1.",
      source_name = "FEED == 3"
    ),
    FORM_DRL17822_ASD = list(
      description = "Amorphous solid dispersion formulation of DRL-17822",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (nanocrystal formulation)",
      notes = "Only study 4 (parts A and B) gave the amorphous solid dispersion, fasted, after a low-fat breakfast or after a high-fat breakfast (Table S2). The formulation was never given after a continental breakfast; the control stream then falls through to the reference bin (F = 1), which is reproduced here but is an extrapolation.",
      source_name = "FORM"
    ),
    OCC = list(
      description = "Dosing occasion for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = "n/a",
      notes = "Values 1-4 (crossover periods). The control stream applies IOV only to subjects with ID > 200, i.e. not to study 1, which had one occasion per subject; set OCC = 0 for study-1-like single-occasion records so every occasion indicator is zero and no IOV is applied.",
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Methods 2.2) but not retained."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Methods 2.2) but not retained; BMI was retained on Vc instead."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Methods 2.2) but not retained. The source data set carried height in m (HGT)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "DRL-17822", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "DRL-17822", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "DRL-17822", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "DRL-17822", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "DRL-17822", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "DRL-17822", units = "mg", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "DRL-17822", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "DRL-17822", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "DRL-17822", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 95L,
    n_studies = 4L,
    age_range = "18-53 years",
    weight_range = "59.4-116.7 kg",
    bmi_range = "18.8-29.9 kg/m^2",
    sex_female_pct = 0,
    disease_state = "healthy male volunteers",
    dose_range = "5-1000 mg single oral dose (study 1); 50-450 mg once daily for 2 weeks (study 3); 150 mg single dose (studies 2 and 4)",
    regions = "The Netherlands",
    notes = "Four phase I studies (study 4 in parts A and B): single ascending dose, food interaction, multiple ascending dose and food-by-formulation crossover (Table S2); 2816 plasma concentrations from 95 drug-treated subjects (Results 3.1); demographics by study in Table S1. Nanocrystal (NC) and amorphous solid dispersion (ASD) formulations, dosed fasted or after a low-fat, high-fat or continental breakfast."
  )

  ini({
    lka <- log(1.74)
    label("Absorption and transit rate constant (1/h)") # Table 1 'ka' = 1.74 1/h [RSE 2.5%]; control stream K12, shared by the depot and all six transit compartments
    lvc <- log(86.7)
    label("Apparent central volume at BMI 23.3 kg/m^2 (L)") # Table 1 'Vc' = 86.7 L [6.7%]
    lkel <- log(0.100)
    label("Elimination rate constant from central (1/h)") # Table 1 'kel' = 0.100 1/h [3.1%]
    lvp <- log(599)
    label("Apparent peripheral volume (L)") # Table 1 'Vp' = 599 L [6.7%]
    lq <- log(5.11)
    label("Apparent intercompartmental clearance (L/h)") # Table 1 'Qc/p' = 5.11 L/h [6.3%]

    lfdepot <- fixed(log(1))
    label("Relative bioavailability, reference bin: nanocrystal after a high-fat or continental breakfast (unitless)") # Table 1 'F reference' = 1.0 [fixed]
    e_fed_medium_fdepot <- 0.532
    label("Relative bioavailability, Fmedium bin: ASD after a low- or high-fat breakfast, or nanocrystal after a low-fat breakfast (unitless)") # Table 1 'F medium' = 0.532 [3.2%]
    e_fasted_asd_fdepot <- 0.151
    label("Relative bioavailability, Fmedium-low bin: ASD fasted (unitless)") # Table 1 'F medium-low' = 0.151 [10.7%]
    e_fasted_nc_fdepot <- 0.056
    label("Relative bioavailability, Flow bin: nanocrystal fasted (unitless)") # Table 1 'F low' = 0.056 [8.7%]

    e_bmi_vc <- 0.041
    label("Fractional change in Vc per kg/m^2 of BMI above 23.3 (1/(kg/m^2))") # Table 1 'COV BMI, Vc' = 0.041 [16.3%]; linear form and 23.3 centre from the supplement control stream
    e_asd_highfat_ka <- 0.599
    label("Fraction of ka for the ASD after a high-fat breakfast (unitless)") # Table 1 'k a,asd,HF' = 0.599 [11.2%]

    # IIV: Table 1 reports CV% = sqrt(exp(omega^2) - 1) (footnote a), so
    # omega^2 = log(1 + CV^2); covariance = correlation * sqrt(var1 * var2).
    # Supplement control stream $OMEGA BLOCK(3) over ETA(1) V2, ETA(2) ka,
    # ETA(3) F1 in that order.
    etalvc + etalka + etalfdepot ~ c(
      0.0688626,
      0.0463897, 0.0625202,
      0, -0.0437934, 0.391274
    ) # Table 1 IIV 'Vc' 26.7%, 'ka' 25.4%, 'F reference' 69.2%; 'Cor. IIV Vc-ka' 0.707; 'Cor. IIV F-ka' -0.280; no Vc-F correlation reported (0)
    etalkel ~ 0.0550978 # Table 1 IIV 'kel' 23.8%; control stream $OMEGA BLOCK(1) ETA(12)

    # IOV: control stream $OMEGA BLOCK(2) on F1 and ka followed by three
    # BLOCK(2) SAME blocks (occasions 2-4 share occasion 1's matrix).
    etaiov_lfdepot_1 + etaiov_lka_1 ~ c(
      0.182908,
      -0.00787035, 0.0427539
    ) # Table 1 IOV 'F reference' 44.8%, 'ka' 20.9%; 'Cor. IOV F-ka' -0.089
    etaiov_lfdepot_2 + etaiov_lka_2 ~ c(
      fixed(0.182908),
      fixed(-0.00787035), fixed(0.0427539)
    ) # control stream '$OMEGA BLOCK (2) SAME' occasion 2
    etaiov_lfdepot_3 + etaiov_lka_3 ~ c(
      fixed(0.182908),
      fixed(-0.00787035), fixed(0.0427539)
    ) # control stream '$OMEGA BLOCK (2) SAME' occasion 3
    etaiov_lfdepot_4 + etaiov_lka_4 ~ c(
      fixed(0.182908),
      fixed(-0.00787035), fixed(0.0427539)
    ) # control stream '$OMEGA BLOCK (2) SAME' occasion 4

    propSd <- 0.337639
    label("Proportional residual error (fraction)") # Table 1 'Proportional error (sigma^2)' = 0.114 [4.8%]; SD = sqrt(0.114); the additive $SIGMA term is 0 FIX in the control stream
  })

  model({
    # Occasion indicators (OCC = 0 -> no IOV, as for study 1 in the source)
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iovlfdepot <- oc1 * etaiov_lfdepot_1 + oc2 * etaiov_lfdepot_2 +
      oc3 * etaiov_lfdepot_3 + oc4 * etaiov_lfdepot_4
    iovlka <- oc1 * etaiov_lka_1 + oc2 * etaiov_lka_2 +
      oc3 * etaiov_lka_3 + oc4 * etaiov_lka_4

    # Food-by-formulation bins for relative bioavailability (control stream FM)
    bin_fasted_nc <- (1 - FED) * (1 - FORM_DRL17822_ASD)
    bin_fasted_asd <- (1 - FED) * FORM_DRL17822_ASD
    bin_fed_medium <- FED_LOWFAT + FED_HIGHFAT * FORM_DRL17822_ASD
    bin_reference <- 1 - bin_fasted_nc - bin_fasted_asd - bin_fed_medium
    frel <- exp(lfdepot) * bin_reference +
      e_fed_medium_fdepot * bin_fed_medium +
      e_fasted_asd_fdepot * bin_fasted_asd +
      e_fasted_nc_fdepot * bin_fasted_nc

    ka <- exp(lka + etalka + iovlka) *
      e_asd_highfat_ka^(FED_HIGHFAT * FORM_DRL17822_ASD)
    vc <- exp(lvc + etalvc) * (1 + (BMI - 23.3) * e_bmi_vc)
    kel <- exp(lkel + etalkel)
    vp <- exp(lvp)
    q <- exp(lq)

    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(transit2) <- ka * transit1 - ka * transit2
    d/dt(transit3) <- ka * transit2 - ka * transit3
    d/dt(transit4) <- ka * transit3 - ka * transit4
    d/dt(transit5) <- ka * transit4 - ka * transit5
    d/dt(transit6) <- ka * transit5 - ka * transit6
    d/dt(central) <- ka * transit6 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- frel * exp(etalfdepot + iovlfdepot)

    # mg / L * 1000 = ng/mL (control stream S2 = V2 / 1000)
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
