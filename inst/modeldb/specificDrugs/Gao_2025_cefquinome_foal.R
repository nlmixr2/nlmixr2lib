Gao_2025_cefquinome_foal <- function() {
  description <- paste0(
    "Preclinical (horse, Ili foal). One-compartment PK model with first-order intramuscular ",
    "absorption for cefquinome, a fourth-generation veterinary cephalosporin, in 7-month to ",
    "1-year-old Ili foals given a single 1 mg/kg dose intravenously or intramuscularly in a ",
    "two-phase crossover. Gao 2025 analysed the serum concentrations non-compartmentally in ",
    "WinNonlin 5.2.1 (Section 2.4) and published NO structural compartmental model, so every ",
    "parameter here is either a directly reported Table 1 statistic or is derived from Table 1 by ",
    "a standard one-compartment identity. CL = 0.09 L/h/kg is the reported intravenous clearance. ",
    "The volume is derived as vc = CL / kel with kel = ln(2) / T1/2beta = ln(2) / 2.35 h, giving ",
    "0.3051 L/kg; this choice is corroborated by Figure 2, whose last quantifiable intravenous ",
    "point sits at about 0.09 ug/mL at 12 h against a model prediction of 0.095 ug/mL, whereas the ",
    "alternative volume CL x MRT0-last = 0.240 L/kg predicts 0.046 ug/mL and is refuted. The ",
    "absorption rate constant ka = 0.6851 1/h is derived by solving the one-compartment Tmax ",
    "identity ln(ka/kel) / (ka - kel) = Tmax at the reported Tmax of 2.16 h; it is corroborated by ",
    "the reported Cmax, since the identity Cmax = F x Dose / vc x exp(-kel x Tmax) returns ",
    "0.760 ug/mL against the published 0.89 +/- 0.14 ug/mL (inside one SD). F = 43.86% is the ",
    "reported absolute bioavailability, and equals the published AUC ratio 5.41 / 12.33 = 43.88%. ",
    "Dose the `central` compartment for the intravenous route and the `depot` compartment for the ",
    "intramuscular route. The model deliberately does NOT reproduce the intramuscular T1/2beta of ",
    "4.16 h reported in Table 1: a flip-flop parameterisation that matches it (ka = ln(2)/4.16 = ",
    "0.1666 1/h) predicts Tmax = 4.45 h against the published 2.16 h and Cmax = 0.387 ug/mL ",
    "against the published 0.89 ug/mL, a 2.3-fold error, and is contradicted by Figure 2, in which ",
    "the intravenous and intramuscular terminal slopes are visibly parallel. See the vignette ",
    "Assumptions and deviations section. Gao 2025 fitted no population model and reported no ",
    "between-subject variability or residual error magnitude -- the +/- values in Table 1 are ",
    "standard deviations of the per-foal non-compartmental estimates -- so there are no eta ",
    "parameters and addSd is FIXED at 0; the model is intended for typical-value simulation."
  )
  reference <- paste(
    "Gao T, Liu X, Qiu D, Li Y, Qiu Z, Qi J, Li S, Guo X, Zhang Y, Wang Z, Gao X, Ma Y, Ma T. (2025).",
    "Ex vivo pharmacokinetic/pharmacodynamic integration model of cefquinome",
    "against Escherichia coli in foals.",
    "Veterinary Sciences 12(4):294.",
    "doi:10.3390/vetsci12040294.",
    sep = " "
  )
  vignette <- "Gao_2025_cefquinome"
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "ug/mL"
  )

  compartmentData <- list(
    depot   = list(analyte = "cefquinome", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cefquinome", units = "mg/kg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species        = "horse (Ili foal)",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "7 months to 1 year",
    weight_range   = "191 +/- 21.7 kg (mean +/- SD)",
    disease_state  = "Healthy. Foals were acclimatised for 2 weeks, passed a fitness assessment (physical examination, heart rate, rectal temperature) and had received no medical treatment for at least 2 weeks before enrolment",
    dose_range     = "Single 1 mg/kg dose of cefquinome sulphate injection, the label-recommended foal dose",
    design         = "Two-phase crossover. Ten foals randomised to two groups of five; in phase 1 one group received 1 mg/kg intramuscularly in the cervical region and the other 1 mg/kg intravenously via the jugular vein, and after a 2-week washout the routes were reversed",
    sampling       = "Left jugular vein blood at 0.083, 0.167, 0.25, 0.5, 0.75, 1, 2, 3, 6, 9, 12 and 24 h after dosing; serum assayed by HPLC with diode-array detection at 268 nm (LLOQ 0.01 ug/mL, LOD 0.005 ug/mL, linear 0.02-5 ug/mL, recovery >80%, CV <10%)",
    regions        = "China (Zhaosu County, Yili, Xinjiang; animals supplied by Zhaosu County Xiyu Horse Industry Co., Ltd.)",
    notes          = "Ethics approval XY20230322 (Animal Ethics Committee of Zhaosu County Xiyu Horse Industry Co., Ltd.); the study reports compliance with the ARRIVE guidelines. Gao 2025 Table 1 reports the non-compartmental results as mean +/- SD over the 10 foals. Intravenous: T1/2beta 2.35 +/- 0.38 h, AUC0-last 12.33 +/- 0.69 h*ug/mL, MRT0-last 2.67 +/- 0.13 h, CL 0.09 +/- 0.01 L/h/kg. Intramuscular: Cmax 0.89 +/- 0.14 ug/mL, Tmax 2.16 +/- 0.75 h, T1/2beta 4.16 +/- 0.21 h, AUC0-last 5.41 +/- 0.81 h*ug/mL, MRT0-last 4.92 +/- 0.15 h, CL 0.15 +/- 0.02 L/h/kg, F 43.86 +/- 5.62%. Although blood was drawn through 24 h, Figure 2 stops at 12 h, so the last quantifiable sample -- and hence the '0-last' window -- is 12 h; with a 2.35 h terminal half-life the 12-24 h contribution is below the 0.01 ug/mL limit of quantification and AUC0-last is numerically interchangeable with AUC0-24h. The intramuscular CL of 0.15 L/h/kg is an apparent CL/F, and is not mutually consistent with the intravenous CL of 0.09 L/h/kg and F of 43.86% (which imply CL/F = 0.205 L/h/kg); this is the usual mean-of-ratios versus ratio-of-means gap across separately averaged per-foal estimates and it propagates into the paper's dose equation, as discussed in the vignette Errata."
  )

  ini({
    # =================================================================
    # Disposition -- Gao 2025 Table 1, IV column
    # =================================================================
    # Clearance is reported directly and is the only disposition
    # parameter taken verbatim from the paper.
    lcl <- log(0.09)
    label("Clearance (L/h/kg)")  # Gao 2025 Table 1, IV column, CL = 0.09 +/- 0.01 L/h/kg

    # DERIVED (not printed in Gao 2025): vc = CL / kel with
    # kel = ln(2) / T1/2beta = ln(2) / 2.35 h = 0.29496 1/h, giving
    # 0.09 / 0.29496 = 0.30513 L/kg. Both inputs are Table 1, IV column
    # values. Gao 2025 reports no volume of distribution of any kind.
    #
    # The competing derivation vc = CL x MRT0-last = 0.09 x 2.67 =
    # 0.2403 L/kg is REFUTED by Figure 2: at 12 h, the last quantifiable
    # sample, the intravenous curve reads about 0.09 ug/mL, against
    # 0.095 ug/mL predicted here and 0.046 ug/mL predicted by the MRT
    # derivation. The two derivations disagree because the intravenous
    # profile in Figure 2 is visibly biexponential over the first 2-3 h,
    # so MRT0-last (2.67 h) is shorter than 1/kel (3.39 h). A two-
    # compartment model needs four disposition parameters and Gao 2025
    # publishes only three independent intravenous summaries, so the
    # distribution phase is not identifiable and is not modelled.
    lvc <- log(0.3051)
    label("Volume of distribution (L/kg)")  # DERIVED from Gao 2025 Table 1, IV column: vc = CL / (log(2) / T1/2beta) = 0.09 / (log(2) / 2.35)

    # =================================================================
    # Intramuscular absorption -- Gao 2025 Table 1, IM column
    # =================================================================
    # DERIVED (not printed in Gao 2025): ka is the unique root of the
    # one-compartment Tmax identity
    #   ln(ka / kel) / (ka - kel) = Tmax
    # at kel = 0.29496 1/h and the Table 1 Tmax of 2.16 h, giving
    # ka = 0.68514 1/h. Both inputs are Table 1 values.
    #
    # Corroborated by the independently reported Cmax through the exact
    # one-compartment identity Cmax = (F x Dose / vc) x exp(-kel x Tmax),
    # which returns 0.760 ug/mL against the published 0.89 +/- 0.14
    # ug/mL -- inside one standard deviation, and a value Gao 2025 was
    # not used to derive ka.
    #
    # The alternative flip-flop root, ka = ln(2) / 4.16 = 0.16662 1/h
    # (which would instead reproduce the Table 1 intramuscular T1/2beta
    # of 4.16 h), is REFUTED twice over: it puts Tmax at 4.45 h against
    # the published 2.16 h, and Cmax at 0.387 ug/mL against the
    # published 0.89 ug/mL, a 2.3-fold underprediction. Figure 2 agrees
    # with the root chosen here, showing parallel intravenous and
    # intramuscular terminal slopes rather than the shallower
    # intramuscular slope flip-flop absorption would produce.
    lka <- log(0.6851)
    label("Intramuscular first-order absorption rate constant (1/h)")  # DERIVED from Gao 2025 Table 1, IM column: root of log(ka/kel)/(ka-kel) = Tmax = 2.16 h at kel = log(2)/2.35

    # Absolute bioavailability, reported directly. Equals the published
    # AUC0-last ratio 5.41 / 12.33 = 43.88%, confirming how it was
    # obtained.
    lfdepot <- log(0.4386)
    label("Absolute bioavailability of the intramuscular dose (fraction)")  # Gao 2025 Table 1, IM column, F = 43.86 +/- 5.62%

    # =================================================================
    # Residual error
    # =================================================================
    # Gao 2025 fitted no population model and reported no residual error
    # magnitude, only the mean +/- SD of the per-foal non-compartmental
    # estimates. The residual SD is therefore held at zero for
    # deterministic typical-value simulation rather than inventing a
    # variance. The Table 1 SDs are reproduced in the vignette so a user
    # can add variability deliberately.
    addSd <- fixed(0)
    label("Additive residual SD on serum concentration (ug/mL; not reported in Gao 2025)")  # Gao 2025 reports no residual error; only NCA mean +/- SD
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    ka <- exp(lka)

    kel <- cl / vc

    # One-compartment model with first-order absorption. Dose `central`
    # for the intravenous route and `depot` for the intramuscular route;
    # the bioavailability below applies only to the depot input, so the
    # intravenous route is unaffected by it.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    f(depot) <- exp(lfdepot)

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
