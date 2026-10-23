Cosson_2018_basimglurant <- function() {
  description <- paste(
    "Two-compartment population PK model for oral modified-release",
    "basimglurant (RO4917523, an mGlu5 negative allosteric modulator) in",
    "healthy adults and adults with major depressive disorder (Cosson 2018,",
    "288 participants pooled from four phase I trials and one phase II",
    "trial). Absorption is a chain of three first-order transit",
    "compartments (depot plus two transit compartments, all at rate Ktr)",
    "after a lag time, with separate Ktr and lag-time estimates in the fed",
    "and fasted states; elimination is first order. Covariates: Asian race",
    "on fasted Ktr and on relative bioavailability, smoking and male sex on",
    "CL/F, body weight and smoking on Vc/F, male sex on Q/F, and body",
    "weight and male sex on Vp/F. The typical subject is a non-Asian female",
    "non-smoker weighing 71 kg. The companion Cmax-driven dizziness model",
    "is Cosson_2018_basimglurant_dizziness."
  )
  reference <- paste(
    "Cosson V, Schaedeli-Stark F, Arab-Alameddine M, Chavanne C, Guerini E,",
    "Derks M, Jaeschke G, Lindemann L, Umbricht D, Santarelli L.",
    "Population Pharmacokinetic and Exposure-dizziness Modeling for a",
    "Metabotropic Glutamate Receptor Subtype 5 Negative Allosteric Modulator",
    "in Major Depressive Disorder Patients.",
    "Clin Transl Sci. 2018;11(5):523-531. doi:10.1111/cts.12566.",
    sep = " "
  )
  vignette <- "Cosson_2018_basimglurant_exposure_dizziness"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline body weight. Power function normalised to 71 kg on Vc/F",
        "(exponent 0.879) and Vp/F (exponent 1.71); no weight effect on",
        "CL/F or Q/F. 71 kg is the reference weight printed in the",
        "covariate equations and the Figure 2 caption (the analysis",
        "population median was 71.9 kg)."
      ),
      source_name = "WT"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female) -- the paper's typical subject",
      notes = paste(
        "Cosson 2018 codes SEX as a male indicator (1 = male, 0 = female)",
        "with female as the typical category (Figure 2 caption and",
        "Figure S3). To keep the paper's female-reference typical values",
        "and coefficients unchanged, model() applies each sex effect to",
        "(1 - SEXF): CL/F * (1 + 0.401 * (1 - SEXF)),",
        "Q/F * (1 - 0.213 * (1 - SEXF)), Vp/F * (1 - 0.434 * (1 - SEXF))."
      ),
      source_name = "SEX"
    ),
    SMOKE = list(
      description = "Current-smoker indicator, 1 = smoker, 0 = non-smoker.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-smoker)",
      notes = paste(
        "Linear-proportional effects on CL/F (1 + 1.08 * SMOKE) and Vc/F",
        "(1 + 0.397 * SMOKE). Cosson 2018 attributes the clearance effect",
        "to CYP1A2 induction by smoking."
      ),
      source_name = "SMOK"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator, 1 = Asian, 0 = other.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = paste(
        "Cosson 2018's ETHN indicator: an effect first seen in Japanese",
        "participants and extended to all Asian participants (Methods).",
        "Slows the FASTED transit rate (Ktr_fasted * (1 - 0.588 * ETHN)) and",
        "raises relative bioavailability by a factor of 1.26. The",
        "bioavailability equation is printed as F = 1 + 1.26 * ETHN, which",
        "would give 2.26; the paper's Figure 2 (Asian -21% on every",
        "apparent parameter, i.e. 1/1.26) and Discussion ('Asians have",
        "slightly higher bioavailability (26%)') both give a factor of 1.26,",
        "which is the encoding used here. See the vignette Errata."
      ),
      source_name = "ETHN"
    ),
    FED = list(
      description = "Fed-state indicator for the dose record, 1 = fed, 0 = fasted.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Switches between two separately estimated absorption parameter",
        "sets: Ktr 4.52 vs 0.999 1/h and lag time 0.212 vs 0.120 h (fasted",
        "vs fed), each with its own between-subject variability on Ktr. The",
        "Asian effect on Ktr applies to the fasted state only. Most of the",
        "analysis population (and all phase II patients) took basimglurant",
        "with food."
      ),
      source_name = "FOOD"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "basimglurant", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "basimglurant", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "basimglurant", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "basimglurant", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "basimglurant", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 288L,
    n_studies = 5L,
    n_observations = 3533L,
    age_median = "43.5 years",
    weight_median = "71.9 kg",
    sex_female_pct = 58,
    race_ethnicity = c(Asian = 17),
    disease_state = paste(
      "88 healthy adults (four phase I trials, including healthy Japanese",
      "and Caucasian participants) and 200 adults with major depressive",
      "disorder with inadequate response to ongoing antidepressant",
      "treatment (phase II trial NCT01437657)"
    ),
    dose_range = paste(
      "0.2, 0.5, 0.7, 1 or 1.5 mg basimglurant modified-release tablet",
      "orally once daily, single and multiple doses"
    ),
    regions = "Phase I sites in France and the United States; phase II sites not stated",
    co_medication = paste(
      "Two phase I trials co-administered fluvoxamine or carbamazepine",
      "(drug-interaction studies); phase II patients continued their",
      "ongoing antidepressant"
    ),
    notes = paste(
      "Cosson 2018 Methods 'Patient population and study design' and",
      "Results 'PK results'. 273 participants received basimglurant with",
      "food, 35 fasted, and 20 under both conditions on two occasions.",
      "Rich PK sampling in healthy participants; sparse sampling (predose,",
      "4 and 6 h on Days 1 and 42; predose and 4 h on Days 14 and 28; any",
      "time on Day 63) in patients. The paper reports no age or weight",
      "ranges; the 5th and 95th weight percentiles were 52 and 100 kg",
      "(Figure 2 caption)."
    )
  )

  ini({
    # Structural parameters (Cosson 2018 Table 1 and the covariate equations
    # in Results 'PK results'); typical subject = female, non-Asian,
    # non-smoker, 71 kg. All clearances and volumes are apparent (/F).
    lcl <- log(8.87); label("Apparent clearance CL/F (L/h)") # Table 1 'Clearance (CL) [L/h]' 8.87, RSE 6.58%
    lvc <- log(231); label("Apparent central volume Vc/F (L)") # Table 1 'Volume of distribution of the central compartment (Vc/V4) [L]' 231, RSE 4.13%
    lq <- log(31.0); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 1 'Intercompartmental clearance (Q) [L/h]' 31.0, RSE 6.23%
    lvp <- log(1140); label("Apparent peripheral volume Vp/F (L)") # Table 1 'Volume of distribution of the peripheral compartment (Vp/V5) [L]' 1140, RSE 5.08%
    lktr_fasted <- log(4.52); label("Transit rate constant Ktr in the fasted state (1/h)") # Table 1 'Ktr under fasted state' 4.52 1/h, RSE 12.8%
    lktr_fed <- log(0.999); label("Transit rate constant Ktr in the fed state (1/h)") # Table 1 'Ktr under fed state' 0.999 1/h, RSE 5.46%
    ltlag_fasted <- log(0.212); label("Absorption lag time in the fasted state (h)") # Table 1 'Lag time under fasted state' 0.212 h, RSE 10.6%
    ltlag_fed <- log(0.120); label("Absorption lag time in the fed state (h)") # Table 1 'Lag time under fed state' 0.120 h, RSE 9.25%

    # Covariate effects. Linear-proportional (1 + e * COV) forms unless noted.
    e_race_asian_ktr_fasted <- -0.588; label("Fractional change in fasted Ktr for Asian participants (unitless)") # Table 1 'Ethnicity on Ktr_fasted' 0.588, RSE 15.2%; sign from the Ktr equation (1 - 0.588 * ETHN)
    e_sex_cl <- 0.401; label("Fractional change in CL/F for male vs female sex (unitless; applied to 1 - SEXF)") # Table 1 'Gender on CL' 0.401, RSE 33.5%; CL/F equation (1 + 0.401 * SEX)
    e_smoke_cl <- 1.08; label("Fractional change in CL/F for current smokers (unitless)") # Table 1 'Smoking status on CL' 1.08, RSE 20.0%; CL/F equation (1 + 1.08 * SMOK)
    e_wt_vc <- 0.879; label("Power exponent of (WT/71) on Vc/F (unitless)") # Table 1 'Body weight on Vc/V4' 0.879, RSE 17.5%; Vc/F equation (WT/71)^0.879
    e_smoke_vc <- 0.397; label("Fractional change in Vc/F for current smokers (unitless)") # Table 1 'Smoking status on Vc/V4' 0.397, RSE 34.5%; Vc/F equation (1 + 0.397 * SMOK)
    e_sex_q <- -0.213; label("Fractional change in Q/F for male vs female sex (unitless; applied to 1 - SEXF)") # Table 1 'Gender on Q' -0.213, RSE 33.5%; Q/F equation (1 - 0.213 * SEX)
    e_wt_vp <- 1.71; label("Power exponent of (WT/71) on Vp/F (unitless)") # Table 1 'Body weight on Vp/V5' 1.71, RSE 9.29%; Vp/F equation (WT/71)^1.71
    e_sex_vp <- -0.434; label("Fractional change in Vp/F for male vs female sex (unitless; applied to 1 - SEXF)") # Table 1 'Gender on Vp/V5' -0.434, RSE 9.38%; Vp/F equation (1 - 0.434 * SEX)
    e_race_asian_fdepot <- 1.26; label("Relative bioavailability of Asian vs non-Asian participants (ratio; applied as e^RACE_ASIAN)") # Table 1 'Ethnicity on relative bioavailability (F)' 1.26, RSE 6.75%; Figure 2 (-21% on every apparent parameter = 1/1.26) and Discussion '26%' fix the ratio reading over the printed F = 1 + 1.26 * ETHN

    # Between-subject variability: Table 1 reports %CV for log-normal etas;
    # omega^2 = log(CV^2 + 1). Covariances = r * sqrt(omega_i^2 * omega_j^2)
    # from the six Table 1 correlations. No BSV on lag time (Results).
    # Order: CL, Vc, Q, Vp. CVs 71.3, 25.9, 36.3, 37.4%; correlations
    # CL-Vc 0.651, CL-Q 0.518, CL-Vp 0.248, Vc-Q 0.721, Vc-Vp 0.350, Q-Vp 0.550.
    etalcl + etalvc + etalq + etalvp ~ c(
      0.41103,
      0.10635, 0.064927,
      0.11684, 0.064636, 0.12378,
      0.057529, 0.032269, 0.070015, 0.13092
    ) # Table 1 'BSV on CL/Vc/Q/Vp [% CV]' and the six 'Correlation between BSV' rows
    etalktr_fasted ~ 0.16175 # Table 1 'BSV on Ktr_fasted [% CV]' 41.9 -> log(0.419^2 + 1)
    etalktr_fed ~ 0.25421 # Table 1 'BSV on Ktr_fed [% CV]' 53.8 -> log(0.538^2 + 1)

    # Residual error: proportional only (Results 'PK results').
    propSd <- 0.262; label("Proportional residual error (fraction)") # Table 1 'Proportional residual error [% CV]' printed 0.262 -- a fraction (26.2%), see vignette Errata
  })

  model({
    # Covariate multipliers (female, non-Asian, non-smoker, 71 kg reference)
    male <- 1 - SEXF

    # Individual parameters
    ktr_fasted <- exp(lktr_fasted + etalktr_fasted) * (1 + e_race_asian_ktr_fasted * RACE_ASIAN)
    ktr_fed <- exp(lktr_fed + etalktr_fed)
    ktr <- ktr_fasted * (1 - FED) + ktr_fed * FED
    tlag <- exp(ltlag_fasted) * (1 - FED) + exp(ltlag_fed) * FED

    cl <- exp(lcl + etalcl) * (1 + e_smoke_cl * SMOKE) * (1 + e_sex_cl * male)
    vc <- exp(lvc + etalvc) * (WT / 71)^e_wt_vc * (1 + e_smoke_vc * SMOKE)
    q <- exp(lq + etalq) * (1 + e_sex_q * male)
    vp <- exp(lvp + etalvp) * (1 + e_sex_vp * male) * (WT / 71)^e_wt_vp
    fdepot <- e_race_asian_fdepot^RACE_ASIAN

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Depot plus two transit compartments, all emptying at Ktr
    # ('a chain of three compartments, including the depot compartment')
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(central) <- ktr * transit2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # mg / L -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
