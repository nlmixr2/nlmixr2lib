Lakshminarayana_2015_cipargamin_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, NMRI, female). Naive-pooled two-compartment population PK",
    "model for oral cipargamin (KAE609, NITD609), the spiroindolone PfATP4 inhibitor,",
    "in Plasmodium berghei-infected and uninfected NMRI mice after a single 5.3 mg/kg",
    "oral dose: first-order absorption with bioavailability assumed complete, and",
    "saturable Michaelis-Menten elimination CL = Vmax / (Km + Cp). Fitted without",
    "between-animal variability. The authors used this model to simulate the",
    "Cmax/TRE, AUC0-48/TRE and %T>TRE PK-PD indices of their dose-fractionation",
    "study (TRE = 2 x the P. berghei IC99 = 48 ng/mL); see the sibling index-driven",
    "parasitemia models Lakshminarayana_2015_cipargamin_cmaxtre,",
    "Lakshminarayana_2015_cipargamin_auctre and Lakshminarayana_2015_cipargamin_ttre."
  )
  reference <- paste(
    "Lakshminarayana SB, Freymond C, Fischli C, Yu J, Weber S, Goh A, Yeung BKS,",
    "Ho PC, Dartois V, Diagana TT, Rottmann M, Blasco F. (2015).",
    "Pharmacokinetic-pharmacodynamic analysis of spiroindolone analogs and KAE609",
    "in a murine malaria model. Antimicrob Agents Chemother 59(2):1200-1210.",
    "doi:10.1128/AAC.03274-14. PMCID: PMC4335872."
  )
  vignette <- "Lakshminarayana_2015_spiroindolones"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "cipargamin", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cipargamin", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cipargamin", units = "mg/kg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse (NMRI, female, specific-pathogen-free)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = paste(
      "Pooled Plasmodium berghei GFP ANKA (PbGFPCON)-infected mice (dosed 24 h after",
      "intravenous inoculation with 2 x 10^7 parasitized erythrocytes) and uninfected",
      "mice; Table 4 shows no difference in PK between the two groups"
    ),
    weight_range = "20-22 g (efficacy-study NMRI mice; PK-study weights not stated)",
    sex_female_pct = 100,
    dose_range = "5.3 mg/kg single oral dose in 10% ethanol, 30% PEG 400, 60% (10% vitamin E TPGS)",
    sampling = "0.25, 0.5, 1, 2, 4, 8, 24 and 48 h post-dose; sparse destructive sampling, about three mice per time point per infection group (Figure 2)",
    regions = "Switzerland (Swiss Tropical and Public Health Institute, Basel); bioanalysis Singapore (Novartis Institute for Tropical Diseases)",
    notes = paste(
      "Naive pooling of the individual concentrations (NONMEM VI, no",
      "intersubject variability). The number of animals is not stated; Figure 2",
      "shows about three per time point per group. NMRI exposure was about 4-fold",
      "higher than in CD-1 mice at a comparable dose, so the NMRI data (the strain",
      "used in the efficacy studies) were used for the model. Oral bioavailability",
      "was assumed to be 100% because no intravenous NMRI data were collected."
    )
  )

  ini({
    lka <- log(1.78)
    label("Absorption rate constant ka (1/h)") # Table 4, Ka = 1.78 1/h
    lvc <- log(4.4)
    label("Apparent central volume of distribution Vc (L/kg)") # Table 4, Vc = 4.4 liters/kg
    lvp <- log(5.1)
    label("Apparent peripheral volume of distribution Vp (L/kg)") # Table 4, Vp = 5.1 liters/kg
    lq <- log(0.61)
    label("Apparent intercompartmental clearance Q (L/h/kg)") # Table 4, Q = 0.61 liters/h/kg
    lvmax <- log(143)
    label("Maximum elimination rate Vmax (ug/h/kg)") # Table 4, Vmax = 143 ug/h/kg
    lkm <- log(34.1)
    label("Michaelis-Menten constant Km (ug/L)") # Table 4, Km = 34.1 ug/liter

    # Naive-pooled fit: the paper reports no residual error magnitude.
    propSd <- fixed(0)
    label("Proportional residual error (fraction; not reported, held at 0)")
  })

  model({
    ka <- exp(lka)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    vmax <- exp(lvmax)
    km <- exp(lkm)

    # Amounts are mg/kg and volumes L/kg, so central / vc is in mg/L (= ug/mL).
    Cc <- central / vc

    # Equation 2: CL = Vmax / (Km + Cp), with Vmax in ug/h/kg and Km and Cp in
    # ug/L. Cp in ug/L is 1000 * Cc, and the resulting CL is in L/h/kg.
    cl <- vmax / (km + 1000 * Cc)
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc ~ prop(propSd)
  })
}
