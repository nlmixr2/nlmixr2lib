Yao_2025_simnotrelvir <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption and an",
    "absorption lag for oral simnotrelvir co-administered with ritonavir",
    "100 mg in Chinese healthy adults and adults with mild-to-moderate",
    "COVID-19 (Yao 2025 'Model 2', the final model pooled over the phase Ia,",
    "Ib and II/III trials; 182 subjects, 1634 concentrations). Creatinine",
    "clearance power effect on CL (reference 113 mL/min); a power effect of",
    "the simnotrelvir dose (mg) and a fed-state factor on Ka; fed-state and",
    "ritonavir-dosed-12-h-earlier factors on relative bioavailability;",
    "exponential IIV on CL, Vc and Ka plus inter-occasion variability on Vc",
    "between the first and second intensive sampling occasions; combined",
    "proportional and additive residual error."
  )
  reference <- paste(
    "Yao BF, Yang Y, Xu SS, Tang BH, Chen J, Guo ZJ, Hu HL, Zhang W, Fu SM,",
    "Zhang XF, Hao GX, Yang XM, Song LL, Ye PP, Liu L, Zhu SW, Zheng Y, Zhao W.",
    "Model-informed drug development in public health emergency of",
    "international concern: accelerating marketing authorization of",
    "simnotrelvir. Antimicrob Agents Chemother. 2025;69(11):e00614-25.",
    "doi:10.1128/aac.00614-25"
  )
  vignette <- "Yao_2025_simnotrelvir"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL as the power term (CRCL / 113)^0.39 (Table S5 footnote a:",
        "'CL (L/h) = 22.6 x (CRCL/113)^0.39'); 113 mL/min is the pooled",
        "cohort median (Table 1, range 64.8-230). The paper reports the",
        "column as 'Creatinine clearance rate, mL/min' without naming the",
        "estimating equation; it is used here exactly as reported, in",
        "absolute mL/min."
      ),
      source_name = "CRCL"
    ),
    DOSE = list(
      description = "Simnotrelvir dose amount per administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters Ka as the un-normalized power term DOSE^-0.229 with DOSE in",
        "mg (Table S5 footnote b: 'Ka (h-1) = 0.755 x DOSE^-0.229 x",
        "FOOD-Ka'). Model 2 doses: 150-1200 mg (Table S1). Time-fixed per",
        "subject in the source trials (one dose level per subject).",
        "Place the event columns (id, time, amt, evid, cmt) before DOSE in",
        "a data frame passed to rxSolve()."
      ),
      source_name = "DOSE"
    ),
    FED = list(
      description = "Fed-vs-fasted state at the simnotrelvir dose",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Multiplies Ka by 1.51 and F by 1.36 when fed (Table S5,",
        "'FOOD-Ka' and 'FOOD-F1', fed vs fasted). Only the phase Ia",
        "food-effect part dosed in the fed state (12 of 182 subjects,",
        "Table 1); the patient trials were fasted."
      ),
      source_name = "Food status"
    ),
    CONMED_RTV_PRE12H = list(
      description = paste(
        "Ritonavir had been dosed 12 h before this simnotrelvir dose",
        "(1) versus given simultaneously with it (0)"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ritonavir given together with the simnotrelvir dose)",
      notes = paste(
        "Multiplies F by 1.15 (Table S5 'Ritonavir-F1': 'Effect of relative",
        "administration time of ritonavir on F1 (Administration time of",
        "ritonavir was before simnotrelvir, compared with simultaneous",
        "administration)'). Supply it on dose records. With twice-daily",
        "co-dosing the first simnotrelvir dose is the only one without",
        "ritonavir on board from 12 h earlier, so set 0 for the first dose",
        "and 1 for every later dose; this reading reproduces the first-dose",
        "and steady-state exposures of Table S6 (see the vignette)."
      ),
      source_name = "Ritonavir-F1"
    ),
    OCC = list(
      description = "Integer-valued occasion indicator for inter-occasion variability on Vc",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Two occasions: 'inter-occasion variability was estimated for V2 to",
        "describe the variability of PK characteristic between first and",
        "second intensive plasma sampling' (Supplement eResults, Model 2",
        "covariate analysis). OCC = 1 is the first-dose sampling occasion,",
        "OCC = 2 the steady-state occasion (Table S1 sampling schedules).",
        "Decomposed inside model() into the indicators oc1 / oc2 that",
        "multiplex the IOV etas, following the Jonsson_2011_ethambutol.R",
        "precedent. Pass OCC = 0 to switch IOV off."
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "simnotrelvir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "simnotrelvir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "simnotrelvir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 182L,
    n_studies = 3L,
    n_observations = 1634L,
    age_range = "18-67 years",
    age_median = "35 years",
    weight_range = "39.0-115 kg",
    weight_median = "65.0 kg",
    sex_female_pct = 35.7,
    race_ethnicity = "Chinese (all three trials were conducted in China)",
    disease_state = paste(
      "Healthy volunteers (phase Ia, n = 60) and adults with mild-to-moderate",
      "COVID-19 (phase Ib, n = 24; phase II/III, n = 98)"
    ),
    dose_range = paste(
      "Simnotrelvir 150-1200 mg with ritonavir 100 mg: single doses of",
      "250/750/1200 mg, twice-daily 150/300/750 mg (phase Ia), 300 or 750 mg",
      "twice daily (phase Ib) and 750 mg twice daily (phase II/III)"
    ),
    regions = "China",
    renal_function = "Creatinine clearance median 113 mL/min (range 64.8-230)",
    notes = paste(
      "Yao 2025 Table 1 ('Model 2' columns). Simnotrelvir monotherapy data",
      "from phase Ia were removed from Model 2 because the drug is marketed",
      "co-packaged with ritonavir. 170 subjects fasted, 12 fed. 87 of 1634",
      "concentrations (5.3%) were below the 5 ng/mL quantification limit",
      "and were replaced by 2.5 ng/mL (Supplement eResults)."
    )
  )

  ini({
    lka <- log(0.755); label("Absorption rate constant coefficient for DOSE in mg, fasted (1/h)") # Table S5: Ka = 0.755 h-1 (RSE 23.7%)
    lcl <- log(22.6); label("Apparent clearance at CRCL = 113 mL/min (L/h)") # Table S5: CL = 22.6 L/h (RSE 4.3%)
    lvc <- log(9.48); label("Apparent central volume of distribution (L)") # Table S5: V2 = 9.48 L (RSE 21.0%)
    lvp <- log(3.92); label("Apparent peripheral volume of distribution (L)") # Table S5: V3 = 3.92 L (RSE 14.7%)
    lq <- log(17.4); label("Apparent intercompartmental clearance (L/h)") # Table S5: Q = 17.4 L/h (RSE 14.4%)
    lfdepot <- fixed(log(1)); label("Relative bioavailability, fasted, ritonavir given with the dose (unitless)") # Table S5: F1 = 1, RSE 0 (reference)
    ltlag <- log(0.447); label("Absorption lag time (h)") # Table S5: ALAG = 0.447 h (RSE 3.4%)

    e_crcl_cl <- 0.39; label("Power exponent of (CRCL/113) on CL (unitless)") # Table S5: CRCL-CL = 0.39 (RSE 26.7%)
    e_dose_ka <- -0.229; label("Power exponent of DOSE (mg) on Ka (unitless)") # Table S5: DOSE-Ka = -0.229 (RSE 16.2%)
    e_fed_ka <- 1.51; label("Multiplicative factor on Ka, fed vs fasted (unitless)") # Table S5: FOOD-Ka = 1.51 (RSE 3.6%)
    e_fed_fdepot <- 1.36; label("Multiplicative factor on F, fed vs fasted (unitless)") # Table S5: FOOD-F1 = 1.36 (RSE 3.0%)
    e_conmed_rtv_pre12h_fdepot <- 1.15; label("Multiplicative factor on F, ritonavir dosed 12 h before vs with the dose (unitless)") # Table S5: Ritonavir-F1 = 1.15 (RSE 3.9%)

    etalcl ~ 0.07037 # Table S5: IIV CL 27.0% CV; omega^2 = log(1 + 0.270^2)
    etalvc ~ 0.51034 # Table S5: IIV V2 81.6% CV; omega^2 = log(1 + 0.816^2)
    etalka ~ 0.02528 # Table S5: IIV Ka 16.0% CV; omega^2 = log(1 + 0.160^2)

    # Inter-occasion variability on Vc, one eta per sampling occasion with a
    # shared variance (NONMEM $OMEGA BLOCK SAME).
    etaiov_vc_1 ~ 0.86736 # Table S5: IOV V2 117.5% CV; omega^2 = log(1 + 1.175^2)
    etaiov_vc_2 ~ fixed(0.86736) # Table S5: IOV V2 117.5% CV; shared variance

    propSd <- 0.359; label("Proportional residual error (fraction)") # Table S5: proportional 35.9% (RSE 3.6%)
    addSd <- 1.463; label("Additive residual error (ng/mL)") # Table S5: additive 1.463 ng/mL (RSE 9.7%)
  })
  model({
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2

    ka <- exp(lka + etalka) * DOSE^e_dose_ka * e_fed_ka^FED
    cl <- exp(lcl + etalcl) * (CRCL / 113)^e_crcl_cl
    vc <- exp(lvc + etalvc + iov_vc)
    vp <- exp(lvp)
    q <- exp(lq)
    tlag <- exp(ltlag)
    fdepot <- exp(lfdepot) * e_fed_fdepot^FED * e_conmed_rtv_pre12h_fdepot^CONMED_RTV_PRE12H

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # Dose in mg and volume in L give mg/L; x 1000 converts to ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
