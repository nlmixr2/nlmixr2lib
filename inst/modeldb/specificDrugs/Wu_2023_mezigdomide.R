Wu_2023_mezigdomide <- function() {
  description <- paste(
    "Two-compartment population PK model for oral mezigdomide (CC-92480),",
    "a cereblon E3 ligase modulator, in healthy subjects pooled from two",
    "phase 1 single-dose studies (Wu 2023). First-order absorption with a",
    "lag time and first-order elimination. Relative bioavailability carries",
    "high-fat-meal, proton-pump-inhibitor (rabeprazole) and formulation B",
    "effects; the absorption rate constant carries high-fat-meal and",
    "formulation-dependent PPI effects; the lag time carries a high-fat-meal",
    "effect. Inter-individual variability on CL, V1 (correlated), Q and V2;",
    "inter-occasion variability (four occasions) on ka and Tlag; combined",
    "additive plus proportional residual error."
  )
  reference <- paste(
    "Wu F, Liu L, Gaudy A, Wang X, Carayannopoulos L, Pourdehnad M, Lamba M.",
    "Model based assessment of food and acid reducing agent effects on oral",
    "absorption of mezigdomide (CC-92480), a novel cereblon E3 ligase",
    "modulator. CPT Pharmacometrics Syst Pharmacol. 2023;12(10):1473-1484.",
    "doi:10.1002/psp4.13024"
  )
  vignette <- "Wu_2023_mezigdomide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    FED_HIGHFAT = list(
      description = "High-fat meal at dosing: 1 = FDA-standard high-fat, high-calorie meal, 0 = fasted",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Source column FOOD (Methods 'Model development and methodology':",
        "'FOOD = 1 if a high-fat meal is given or otherwise 0'). Only the",
        "Study I part 2 food-effect crossover (0.8 mg, N = 16) dosed fed.",
        "Text S2 codes the raw column as fasting = 1 / fed = 2; supply the",
        "0/1 form here. Acts on F, ka and Tlag."
      ),
      source_name = "FOOD"
    ),
    CONMED_PPI = list(
      description = "Co-administration with the proton-pump inhibitor rabeprazole: 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no PPI)",
      notes = paste(
        "Source column PPI. Study II periods 3 and 4 dosed mezigdomide after",
        "at least 7 days of rabeprazole (Methods 'Study designs'), so the",
        "indicator is set per dosing occasion. Acts on F directly and, split",
        "by formulation (the paper's PPIA = PPI * (1 - FORMU) and",
        "PPIB = PPI * FORMU), on ka."
      ),
      source_name = "PPI"
    ),
    FORM_MEZIGDOMIDE_B = list(
      description = "Mezigdomide formulation B: 1 = formulation B, 0 = formulation A",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (formulation A, the first-in-human formulation)",
      notes = paste(
        "Source column FORMU ('FORMU = 0 if mezigdomide is dosed in",
        "formulation A and FORMU = 1 if in formulation B'). Text S2 codes",
        "the raw column as A = 1 / B = 2; supply the 0/1 form here. Only",
        "Study II dosed formulation B. The paper labels formulation A the",
        "reference (Methods) and the test formulation (Figure 1 caption);",
        "the coding above does not depend on which label is right. Acts on",
        "F, and on ka through the PPI-by-formulation split."
      ),
      source_name = "FORMU"
    ),
    OCC = list(
      description = "Occasion (dosing period) index, 1-4, for the inter-occasion variability on ka and Tlag",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Each single dose is its own occasion: Study I part 1 has one",
        "occasion, Study I part 2 two crossover periods, Study II four",
        "periods (Figure 1). Records with OCC outside 1-4 carry no IOV."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL, V1 and V2 (Table S2, runs PK17-PK19) and not retained (BICc did not improve).",
      source_name = "WEIGHT"
    ),
    SEXF = list(
      description = "Sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Tested on CL and V1 (Table S2, runs PK15-PK16) and not retained. Text S2 codes the raw column SEXN as male = 1, female = 2.",
      source_name = "SEXN"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL (Table S2, run PK20; range 51.3-112 umol/L) and not retained.",
      source_name = "CREAT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on V2 (Table S2, run PK21) and not retained.",
      source_name = "AGE"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "mezigdomide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mezigdomide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mezigdomide", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 64L,
    n_studies = 2L,
    age_range = "18-55 years (mean 38.5)",
    weight_range = "54.8-103 kg (mean 80.1)",
    sex_female_pct = 9.4,
    race_ethnicity = c(White = 56.3, Black = 40.6, Asian = 1.6, `American Indian or Alaska Native` = 1.6),
    ethnicity = "Hispanic or Latino 40.6%",
    disease_state = "Healthy subjects (men, and women of non-childbearing potential)",
    dose_range = paste(
      "Single oral doses: 0.4, 0.8, 1.6 or 3.2 mg fasted (Study I part 1);",
      "0.8 mg fasted and fed crossover (Study I part 2); 1.6 mg of",
      "formulation A or B, each with and without rabeprazole (Study II)"
    ),
    regions = "United States",
    n_observations = "2038 plasma concentrations, 460 (22.5%) below the 0.05 ng/mL LLOQ and handled as left-censored",
    notes = paste(
      "Studies NCT03803644 (N = 40) and NCT04211545 (N = 24); demographics",
      "from Wu 2023 Table 1. Estimation by SAEM in Monolix 2020R1."
    )
  )

  ini({
    # Structural parameters: Wu 2023 Figure 3a, 'Fixed effects' (final model).
    lka <- log(1.18); label("Absorption rate constant ka (1/h), fasted, no PPI")  # Figure 3a, 'ka pop (1/h)' = 1.18 (RSE 8.99%)
    ltlag <- log(0.423); label("Absorption lag time Tlag (h), fasted")  # Figure 3a, 'Tlag pop (h)' = 0.423 (RSE 3.32%)
    lcl <- log(35.1); label("Apparent clearance CL/F (L/h)")  # Figure 3a, 'CL pop (L/h)' = 35.10 (RSE 5.10%)
    lvc <- log(440); label("Apparent central volume V1/F (L)")  # Figure 3a, 'V1 pop (L)' = 440 (RSE 8.52%)
    lq <- log(36.8); label("Apparent intercompartmental clearance Q/F (L/h)")  # Figure 3a, 'Q pop (L/h)' = 36.80 (RSE 8.45%)
    lvp <- log(243); label("Apparent peripheral volume V2/F (L)")  # Figure 3a, 'V2 pop (L)' = 243 (RSE 6.28%)
    lfdepot <- fixed(log(1)); label("Relative bioavailability F, formulation A fasted without PPI (unitless)")  # Figure 3a, 'F pop' = 1 (fixed)

    # Covariate effects on the log scale (Methods equations for log F, log ka
    # and log Tlag; Figure 3a beta values).
    e_fed_highfat_fdepot <- 0.234; label("High-fat-meal effect on log F (unitless)")  # Figure 3a, 'beta F,FOOD' = 0.234 (RSE 12.9%)
    e_form_mezigdomide_b_fdepot <- -0.236; label("Formulation B effect on log F (unitless)")  # Figure 3a, 'beta F,FORMU' = -0.236 (RSE 7.88%)
    e_conmed_ppi_fdepot <- -1.06; label("PPI effect on log F (unitless)")  # Figure 3a, 'beta F,PPI' = -1.06 (RSE 2.18%)
    e_fed_highfat_tlag <- 1.11; label("High-fat-meal effect on log Tlag (unitless)")  # Figure 3a, 'beta Tlag,FOOD' = 1.11 (RSE 8.09%)
    e_fed_highfat_ka <- -1.15; label("High-fat-meal effect on log ka (unitless)")  # Figure 3a, 'beta ka,FOOD' = -1.15 (RSE 17.3%)
    e_conmed_ppi_forma_ka <- -0.592; label("PPI effect on log ka with formulation A (unitless)")  # Figure 3a, 'beta ka,PPIA' = -0.592 (RSE 32.2%)
    e_conmed_ppi_formb_ka <- -1.99; label("PPI effect on log ka with formulation B (unitless)")  # Figure 3a, 'beta ka,PPIB' = -1.99 (RSE 9.48%)

    # Inter-individual variability. Figure 3a reports standard deviations
    # (omega); variances are omega^2 and the CL-V1 covariance is
    # corr * omega_CL * omega_V1 = 0.862 * 0.399 * 0.634.
    etalcl + etalvc ~ c(0.159201, 0.218057, 0.401956) # Figure 3a, omega CL = 0.399, omega V1 = 0.634, corr V1,CL = 0.862
    etalq ~ 0.068644 # Figure 3a, omega Q = 0.262
    etalvp ~ 0.065025 # Figure 3a, omega V2 = 0.255

    # Inter-occasion variability on ka and Tlag (Text S1; Figure 3a gamma
    # values are standard deviations). The paper estimates no separate IIV on
    # ka or Tlag, so these occasion-level etas carry both (Text S1). One slot
    # per occasion; occasions 2-4 repeat occasion 1's variance.
    etaiov_ka_1 ~ 0.4225 # Figure 3a, gamma ka = 0.650
    etaiov_ka_2 ~ fixed(0.4225) # same variance as occasion 1
    etaiov_ka_3 ~ fixed(0.4225) # same variance as occasion 1
    etaiov_ka_4 ~ fixed(0.4225) # same variance as occasion 1
    etaiov_tlag_1 ~ 0.0961 # Figure 3a, gamma Tlag = 0.310
    etaiov_tlag_2 ~ fixed(0.0961) # same variance as occasion 1
    etaiov_tlag_3 ~ fixed(0.0961) # same variance as occasion 1
    etaiov_tlag_4 ~ fixed(0.0961) # same variance as occasion 1

    # Residual error: Monolix combined2, y = f + e * sqrt(a^2 + b^2 * f^2) (Text S1).
    addSd <- 0.0257; label("Additive residual SD (ng/mL)")  # Figure 3a, error model 'a' = 0.0257 (RSE 6.46%)
    propSd <- 0.210; label("Proportional residual SD (fraction)")  # Figure 3a, error model 'b' = 0.210 (RSE 2.14%)
  })

  model({
    # Occasion indicators and inter-occasion etas.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 + oc4 * etaiov_ka_4
    iov_tlag <- oc1 * etaiov_tlag_1 + oc2 * etaiov_tlag_2 + oc3 * etaiov_tlag_3 + oc4 * etaiov_tlag_4

    # Formulation-specific PPI indicators (the paper's PPIA and PPIB).
    ppi_forma <- CONMED_PPI * (1 - FORM_MEZIGDOMIDE_B)
    ppi_formb <- CONMED_PPI * FORM_MEZIGDOMIDE_B

    # Individual parameters (Methods, final-model covariate equations).
    ka <- exp(lka + e_fed_highfat_ka * FED_HIGHFAT + e_conmed_ppi_forma_ka * ppi_forma +
      e_conmed_ppi_formb_ka * ppi_formb + iov_ka)
    tlag <- exp(ltlag + e_fed_highfat_tlag * FED_HIGHFAT + iov_tlag)
    fdepot <- exp(lfdepot + e_fed_highfat_fdepot * FED_HIGHFAT + e_conmed_ppi_fdepot * CONMED_PPI +
      e_form_mezigdomide_b_fdepot * FORM_MEZIGDOMIDE_B)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Text S1 structural model (A0 = depot, A1 = central, A2 = peripheral1).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # Dose in mg and volume in L give mg/L; x1000 gives ng/mL (Figure 2 axis).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
