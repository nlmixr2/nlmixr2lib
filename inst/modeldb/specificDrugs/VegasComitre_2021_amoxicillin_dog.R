VegasComitre_2021_amoxicillin_dog <- function() {
  description <- paste(
    "Veterinary (dog). Three-compartment population PK model for",
    "intravenous amoxicillin (given as amoxicillin-clavulanic acid) in",
    "12 healthy laboratory beagles (single 20 mg/kg IV bolus) and 12",
    "critically ill client-owned dogs in intensive care (20 mg/kg as a",
    "0.5 h IV infusion every 8 h for at least 48 h), fitted jointly in",
    "Phoenix NLME. All disposition parameters are body-weight-normalised",
    "(L/kg, L/h/kg), so the dose supplied to this model is mg of",
    "amoxicillin per kg. Critical illness (DIS_CRITILL) lowers clearance",
    "to 43.6% of the healthy value and lowers the clearance to the",
    "superficial peripheral compartment by exp(-9.946), which removes",
    "that compartment from the disposition of sick dogs. Between-occasion",
    "variability on clearance applies to sick dogs only, with a larger",
    "variance on occasion 3. Residual error is combined additive plus",
    "proportional with separate magnitudes for healthy and sick dogs.",
    "Protein binding was negligible, so Cc is also the free concentration."
  )
  reference <- paste(
    "Vegas Comitre MD, Cortellini S, Cherlet M, Devreese M, Roques BB,",
    "Bousquet-Melou A, Toutain P-L, Pelligand L. (2021).",
    "Population Pharmacokinetics of Intravenous Amoxicillin Combined With",
    "Clavulanic Acid in Healthy and Critically Ill Dogs.",
    "Frontiers in Veterinary Science 8:770202.",
    "doi:10.3389/fvets.2021.770202.",
    sep = " "
  )
  vignette <- "VegasComitre_2021_amoxicillin_dog"

  # Amounts are mg of amoxicillin per kg of body weight and volumes are L/kg,
  # so Cc is mg/L. The paper reports concentrations in ng/mL; the two
  # additive residual SDs below are converted from ng/mL to mg/L (/1000).
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "mg/L"
  )

  compartmentData <- list(
    central = list(analyte = "amoxicillin", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "amoxicillin", units = "mg/kg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "amoxicillin", units = "mg/kg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    DIS_CRITILL = list(
      description = "Critical illness (1 = critically ill dog hospitalised in the intensive care unit; 0 = healthy dog)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy)",
      notes = paste(
        "Vegas Comitre 2021 calls this covariate 'Health'. Table 4 footnote:",
        "theta(Health) on Cl = -0.829 means that the clearance in sick dogs",
        "is the healthy value times exp(-0.829), i.e. 0.147 L/kg/h or 43.6%",
        "of control, so the paper's coding is 1 = sick and no sign change is",
        "needed. The contrast is fully confounded with study and breed: all",
        "healthy dogs were intact female laboratory beagles (Toulouse) and",
        "all sick dogs were mixed-breed client-owned ICU patients (Royal",
        "Veterinary College). 11 of 12 sick dogs met SIRS criteria and 10 of",
        "12 had a confirmed infectious focus. Also gates the between-occasion",
        "variability and selects the residual-error magnitudes."
      ),
      source_name = "Health"
    ),
    OCC = list(
      description = "Occasion index for between-occasion variability on clearance in sick dogs (1-4)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Only read when DIS_CRITILL = 1. Estimation occasions (Methods,",
        "Pharmacokinetic Analysis): 1 = doses before and up to the first PK",
        "sample; 2 = doses including the PK period with its 8 h trough and",
        "the next dose; 3 = doses including the 24 h trough and the next",
        "dose; 4 = doses including the 48 h trough and onwards. The paper's",
        "Monte Carlo simulations instead used 24 h blocks from the first",
        "dose (1 = 0-24 h, 2 = 24-48 h, 3 = 48-72 h, 4 = 72-96 h). Occasion",
        "3 carries its own, larger variance; occasions 1, 2 and 4 share one.",
        "Any other value (e.g. 0) switches the occasion term off."
      ),
      source_name = "Occasion"
    )
  )

  # Only Health and Occasion entered the stepwise covariate search. The
  # clinical covariates below were tested afterwards by Spearman correlation
  # against the model-derived AUC0-64h of the sick dogs and none was
  # significant (the APPLEfast score, vasopressor use, crystalloid rate and
  # AKI grade were tested the same way; they have no register entry).
  covariatesDataExcluded <- list(
    ALB = list(
      description = "Serum albumin at admission",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Spearman correlation with AUC0-64h not significant; not tested as",
        "a model covariate. Sick-dog median 21.6 g/L (IQR 16.2-23.3), Table 3."
      )
    ),
    CREAT = list(
      description = "Serum creatinine at admission",
      units = "umol/L",
      type = "continuous",
      notes = paste(
        "Spearman correlation with AUC0-64h not significant; not tested as",
        "a model covariate. Sick-dog median 62.5 umol/L (IQR 39.0-132.5),",
        "Table 3."
      )
    )
  )

  population <- list(
    species = "dog (healthy laboratory beagles + critically ill client-owned mixed breeds)",
    n_subjects = 24L,
    n_studies = 2L,
    n_observations = 218L,
    age_range = "healthy: median 24 months; sick: average 40.2 months (range 4.8-114)",
    weight_range = "healthy: average 11.5 kg (9.9-13.2); sick: average 20.6 kg (11.4-42.0)",
    sex_female_pct = 58.3,
    disease_state = paste(
      "12 healthy intact female beagles, and 12 critically ill dogs in the",
      "intensive care unit: septic peritonitis (7), pyothorax (3), acute",
      "haemorrhagic diarrhoea syndrome (1), burns (1); acute kidney injury",
      "in 2 and mechanical ventilation in 1 (Table 2)."
    ),
    dose_range = paste(
      "Amoxicillin-clavulanic acid 20 mg/kg (16.667 mg/kg amoxicillin).",
      "Healthy: single IV bolus (amoxicillin 16.95 +/- 0.37 mg/kg).",
      "Sick: 0.5 h IV infusion every 8 h for at least 48 h."
    ),
    regions = "United Kingdom (Royal Veterinary College, sick dogs); France (Toulouse Veterinary School, healthy dogs)",
    notes = paste(
      "Table 1 (study designs) and Table 2 (sick-dog demographics). Sick",
      "dogs: 4 male, 6 male neutered, 1 female, 1 female spayed. Healthy",
      "dogs sampled at 0, 0.03, 0.17, 0.42, 0.67, 1, 2, 4, 6, 8, 10 and",
      "12 h; sick dogs at end of infusion, 1, 2 and 4 h plus troughs at",
      "8, 24 and 48 h. LLOQ 50 ng/mL, handled by the M3 method. Plasma",
      "protein binding measured by ultrafiltration was negligible."
    )
  )

  ini({
    # --- Structural disposition of healthy (reference) dogs, per kg of body
    # weight (Vegas Comitre 2021 Table 4, 'Estimate' column). V2/Cl2 are the
    # superficial and V3/Cl3 the deep peripheral compartment (Results).
    lcl <- log(0.336)
    label("Clearance, healthy dog (L/h/kg)") # Table 4 Cl = 0.336 L/kg/h [bootstrap median 0.332 (0.293-0.384)]
    lvc <- log(0.173)
    label("Central volume of distribution (L/kg)") # Table 4 V1 = 0.173 L/kg [bootstrap 0.182 (0.119-0.212)]
    lq <- log(0.861)
    label("Intercompartmental clearance to superficial peripheral compartment, healthy dog (L/h/kg)") # Table 4 Cl2 = 0.861 L/kg/h [bootstrap 0.796 (0.717-1.007)]
    lvp <- log(0.174)
    label("Superficial peripheral volume of distribution (L/kg)") # Table 4 V2 = 0.174 L/kg [bootstrap 0.172 (0.166-0.193)]
    lq2 <- log(0.0373)
    label("Intercompartmental clearance to deep peripheral compartment (L/h/kg)") # Table 4 Cl3 = 0.0373 L/kg/h [bootstrap 0.0368 (0.0291-0.0515)]
    lvp2 <- log(0.0776)
    label("Deep peripheral volume of distribution (L/kg)") # Table 4 V3 = 0.0776 L/kg [bootstrap 0.0754 (0.0627-0.0915)]

    # --- Critical-illness effects, exponential on the log scale (Table 4
    # footnote: sick-dog Cl = healthy Cl * exp(-0.829) = 0.147 L/kg/h).
    e_dis_critill_cl <- -0.829
    label("Log-scale effect of critical illness on clearance (unitless)") # Table 4 theta(Health) on Cl = -0.829 [bootstrap -0.795 (-0.975, -0.570)]
    e_dis_critill_q <- -9.946
    label("Log-scale effect of critical illness on superficial intercompartmental clearance (unitless)") # Table 4 theta(Health) on Cl2 = -9.946 [bootstrap -8.665 (-11.989, -6.910)]

    # --- Between-subject variability. Table 4 reports CV%; Eq. (2) defines
    # CV = 100 * sqrt(exp(omega^2) - 1), so omega^2 = log(1 + CV^2).
    # Off-diagonal terms were not reported, so the matrix is diagonal.
    etalcl ~ 0.0399935 # Table 4 BSV Cl 20.2% CV; log(1 + 0.202^2)
    etalvc ~ 0.309251 # Table 4 BSV V1 60.2% CV; log(1 + 0.602^2)
    etalq ~ 0.164606 # Table 4 BSV Cl2 42.3% CV; log(1 + 0.423^2)
    etalvp ~ 0.00842838 # Table 4 BSV V2 9.2% CV; log(1 + 0.092^2)

    # --- Between-occasion variability on clearance, sick dogs only. Occasions
    # 1, 2 and 4 share one variance; occasion 3 has its own (Table 4).
    etaiov_cl_1 ~ 0.0321307 # Table 4 BOV on Cl (all but Occasion 3) 18.07% CV; log(1 + 0.1807^2)
    etaiov_cl_2 ~ fixed(0.0321307) # shared with occasion 1
    etaiov_cl_3 ~ 0.16346 # Table 4 BOV on Cl (Occasion 3) 42.14% CV; log(1 + 0.4214^2)
    etaiov_cl_4 ~ fixed(0.0321307) # shared with occasion 1

    # --- Residual error, combined additive + proportional, estimated
    # separately for healthy and sick dogs (Table 4). Additive SDs converted
    # from ng/mL to mg/L.
    propSd_healthy <- 0.042
    label("Proportional residual error, healthy dogs (fraction)") # Table 4 proportional error healthy 4.2% [bootstrap 4.1% (3.4-5.1%)]
    addSd_healthy <- 0.0306
    label("Additive residual error, healthy dogs (mg/L)") # Table 4 additive error healthy 30.6 ng/mL [bootstrap 24.9 (8.6-43.1)]
    propSd_critill <- 0.359
    label("Proportional residual error, critically ill dogs (fraction)") # Table 4 proportional error sick 35.9% [bootstrap 34.6% (25.2-48.7%)]
    addSd_critill <- 0.00358
    label("Additive residual error, critically ill dogs (mg/L)") # Table 4 additive error sick 3.58 ng/mL [bootstrap 3.60 (2.53-7.02)]
  })

  model({
    # 1. Occasion indicators for the sick-dog between-occasion variability.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_cl <- DIS_CRITILL *
      (oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4)

    # 2. Individual parameters (Eq. 1, P_i = tvP * exp(eta_Pi)); critical
    #    illness enters exponentially on Cl and Cl2.
    cl <- exp(lcl + etalcl + iov_cl + e_dis_critill_cl * DIS_CRITILL)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq + e_dis_critill_q * DIS_CRITILL)
    vp <- exp(lvp + etalvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 4. Three-compartment disposition; IV bolus or infusion into central.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # 5. Observation. Protein binding was negligible, so Cc is total and free.
    Cc <- central / vc

    propSd_i <- propSd_healthy * (1 - DIS_CRITILL) + propSd_critill * DIS_CRITILL
    addSd_i <- addSd_healthy * (1 - DIS_CRITILL) + addSd_critill * DIS_CRITILL
    Cc ~ add(addSd_i) + prop(propSd_i)
  })
}
