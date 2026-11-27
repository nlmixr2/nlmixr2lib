Yao_2025_simnotrelvir_healthy <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption and an",
    "absorption lag for oral simnotrelvir given alone or with ritonavir",
    "100 mg in Chinese healthy adults (Yao 2025 'Model 1', fitted to the",
    "first-in-human phase Ia trial; 72 subjects, 1742 concentrations) and",
    "used to select the 750/100 mg twice-daily patient dose. Ritonavir",
    "co-administration reduces CL to 24.1% and raises relative",
    "bioavailability (by 36% when ritonavir is given with the simnotrelvir",
    "dose, 61% when it was given 12 h before); Ka follows a power function",
    "of the simnotrelvir dose (mg) with a separate exponent with and without",
    "ritonavir and a fed-state factor; F has a fed-state factor; exponential",
    "IIV on CL, Vc and Ka plus inter-occasion variability on Vc; combined",
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
    CONMED_RTV = list(
      description = "Ritonavir 100 mg co-administration (1) versus simnotrelvir monotherapy (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (simnotrelvir monotherapy)",
      notes = paste(
        "Multiplies CL by 0.241 (Table 2 footnote b: 'CL (L/h) = 174 x",
        "Ritonavir-CL'), switches the dose exponent on Ka from -0.127 to",
        "-0.138 (Table 2 footnotes d/e) and switches on the ritonavir F",
        "factors (Table 2 footnotes g/h). Phase Ia: 24 monotherapy and 60",
        "ritonavir subjects (12 subjects contributed to both, Table 1)."
      ),
      source_name = "Ritonavir"
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
        "Read only when CONMED_RTV = 1. Selects the ritonavir F factor:",
        "1.36 when ritonavir is given with the simnotrelvir dose (Table 2",
        "footnote g, 'administered with the first simnotrelvir dose') and",
        "1.61 when it was given 12 h before (footnote h, 'administered 12 h",
        "before simnotrelvir'), both relative to monotherapy. Supply it on",
        "dose records: 1 when a ritonavir dose was given 12 h before this",
        "simnotrelvir dose. The authors' twice-daily dose-selection",
        "simulation (Table S3) is reproduced with 0 on every dose; coding",
        "later doses 1 raises steady-state exposure by about 20-25% (see the",
        "vignette)."
      ),
      source_name = "Ritonavir-F1"
    ),
    DOSE = list(
      description = "Simnotrelvir dose amount per administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters Ka as the un-normalized power term DOSE^theta with DOSE in",
        "mg (Table 2 footnote c: 'Ka (h-1) = 0.425 x DOSE^-0.127 or -0.138",
        "x FOOD-Ka'). Model 1 doses: 150-3000 mg (Table S1). Time-fixed per",
        "subject in the source trial. Place the event columns (id, time,",
        "amt, evid, cmt) before DOSE in a data frame passed to rxSolve()."
      ),
      source_name = "DOSE"
    ),
    FED = list(
      description = "Fed-vs-fasted state at the simnotrelvir dose",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Multiplies Ka by 1.49 and F by 1.39 when fed (Table 2 'FOOD-Ka'",
        "and 'FOOD-F1', fed vs fasted). Only the phase Ia food-effect part",
        "(750/100 mg single dose) dosed in the fed state."
      ),
      source_name = "Food status"
    ),
    OCC = list(
      description = "Integer-valued occasion indicator for inter-occasion variability on Vc",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Two occasions: 'inter-occasion variability was estimated for V2 to",
        "describe the variability of PK characteristic between first and",
        "second intensive plasma sampling' (Supplement eResults, Model 1",
        "covariate analysis). OCC = 1 is the first-dose sampling occasion,",
        "OCC = 2 the last-dose (steady-state) occasion of the multiple-dose",
        "parts (Table S1). Decomposed inside model() into the indicators",
        "oc1 / oc2 that multiplex the IOV etas, following the",
        "Jonsson_2011_ethambutol.R precedent. Pass OCC = 0 to switch IOV off."
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
    n_subjects = 72L,
    n_studies = 1L,
    n_observations = 1742L,
    age_range = "20-50 years",
    age_median = "30 years",
    weight_range = "46.1-82.9 kg",
    weight_median = "64.4 kg",
    sex_female_pct = 26.4,
    race_ethnicity = "Chinese (phase Ia trial conducted in China)",
    disease_state = "Healthy volunteers",
    dose_range = paste(
      "Simnotrelvir alone 150-3000 mg single dose or 750-1500 mg twice daily;",
      "with ritonavir 100 mg, 250-1200 mg single dose, 150-750 mg twice",
      "daily, and 750 mg single dose fasted and fed (food-effect part)"
    ),
    regions = "China",
    renal_function = "Creatinine clearance median 117 mL/min (range 72.9-178)",
    notes = paste(
      "Yao 2025 Table 1 ('Model 1' column) and Table S1 (NCT05339646). The",
      "main text quotes 60 volunteers for Model 1 while Table 1 and the",
      "Supplement eResults give 72; 60 is the ritonavir-treated subset. 95 of",
      "1742 concentrations (5.5%) were below the 5 ng/mL quantification",
      "limit and were replaced by 2.5 ng/mL (Supplement eResults)."
    )
  )

  ini({
    lka <- log(0.425); label("Absorption rate constant coefficient for DOSE in mg, fasted (1/h)") # Table 2: Ka = 0.425 h-1 (RSE 22.5%)
    lcl <- log(174); label("Apparent clearance, simnotrelvir monotherapy (L/h)") # Table 2: CL = 174 L/h (RSE 7.0%)
    lvc <- log(33.7); label("Apparent central volume of distribution (L)") # Table 2: V2 = 33.7 L (RSE 26.3%)
    lvp <- log(3.84); label("Apparent peripheral volume of distribution (L)") # Table 2: V3 = 3.84 L (RSE 17.5%)
    lq <- log(58.6); label("Apparent intercompartmental clearance (L/h)") # Table 2: Q = 58.6 L/h (RSE 9.9%)
    lfdepot <- fixed(log(1)); label("Relative bioavailability, fasted monotherapy (unitless)") # Table 2: F1 = 1, RSE 0 (reference)
    ltlag <- log(0.223); label("Absorption lag time (h)") # Table 2: ALAG = 0.223 h (RSE 2.0%)

    e_conmed_rtv_cl <- 0.241; label("Multiplicative factor on CL with ritonavir (unitless)") # Table 2: Ritonavir-CL = 0.241 (RSE 20.6%)
    e_dose_ka_mono <- -0.127; label("Power exponent of DOSE (mg) on Ka, simnotrelvir alone (unitless)") # Table 2: DOSE-Ka (footnote d, simnotrelvir alone) = -0.127 (RSE 25.6%)
    e_dose_ka_rtv <- -0.138; label("Power exponent of DOSE (mg) on Ka, with ritonavir (unitless)") # Table 2: DOSE-Ka (footnote e, with ritonavir) = -0.138 (RSE 26.4%)
    e_fed_ka <- 1.49; label("Multiplicative factor on Ka, fed vs fasted (unitless)") # Table 2: FOOD-Ka = 1.49 (RSE 3.7%)
    e_fed_fdepot <- 1.39; label("Multiplicative factor on F, fed vs fasted (unitless)") # Table 2: FOOD-F1 = 1.39 (RSE 3.0%)
    e_conmed_rtv_fdepot <- 1.36; label("Multiplicative factor on F, ritonavir given with the dose vs monotherapy (unitless)") # Table 2: Ritonavir-F1 (footnote g) = 1.36 (RSE 24.0%)
    e_conmed_rtv_pre12h_fdepot <- 1.61; label("Multiplicative factor on F, ritonavir dosed 12 h before vs monotherapy (unitless)") # Table 2: Ritonavir-F1 (footnote h) = 1.61 (RSE 23.2%)

    etalcl ~ 0.07087 # Table 2: IIV CL 27.1% CV; omega^2 = log(1 + 0.271^2)
    etalvc ~ 0.52408 # Table 2: IIV V2 83.0% CV; omega^2 = log(1 + 0.830^2)
    etalka ~ 0.03695 # Table 2: IIV Ka 19.4% CV; omega^2 = log(1 + 0.194^2)

    # Inter-occasion variability on Vc, one eta per sampling occasion with a
    # shared variance (NONMEM $OMEGA BLOCK SAME).
    etaiov_vc_1 ~ 0.39034 # Table 2: IOV V2 69.1% CV; omega^2 = log(1 + 0.691^2)
    etaiov_vc_2 ~ fixed(0.39034) # Table 2: IOV V2 69.1% CV; shared variance

    propSd <- 0.46; label("Proportional residual error (fraction)") # Table 2: proportional 46.0% (RSE 3.7%)
    addSd <- 1.556; label("Additive residual error (ng/mL)") # Table 2: additive 1.556 ng/mL (RSE 11.1%)
  })
  model({
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2

    e_dose_ka <- e_dose_ka_mono * (1 - CONMED_RTV) + e_dose_ka_rtv * CONMED_RTV
    ka <- exp(lka + etalka) * DOSE^e_dose_ka * e_fed_ka^FED
    cl <- exp(lcl + etalcl) * e_conmed_rtv_cl^CONMED_RTV
    vc <- exp(lvc + etalvc + iov_vc)
    vp <- exp(lvp)
    q <- exp(lq)
    tlag <- exp(ltlag)
    rtv_with <- CONMED_RTV * (1 - CONMED_RTV_PRE12H)
    rtv_pre <- CONMED_RTV * CONMED_RTV_PRE12H
    fdepot <- exp(lfdepot) * e_fed_fdepot^FED * e_conmed_rtv_fdepot^rtv_with *
      e_conmed_rtv_pre12h_fdepot^rtv_pre

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
