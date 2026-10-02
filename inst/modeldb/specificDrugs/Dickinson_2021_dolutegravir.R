Dickinson_2021_dolutegravir <- function() {
  description <- "Maternal population PK model for oral dolutegravir in women living with HIV who started treatment late in pregnancy (DolPHIN-1), fitted simultaneously to maternal plasma (third trimester, delivery and postpartum), umbilical cord and breast-milk concentrations. Two-compartment disposition with first-order absorption; the maternal central compartment is linked by first-order rate constants to a fetal (umbilical cord) compartment of negligible volume that does not deplete the mother, and to a breast-milk compartment of fixed 0.125 L volume. No covariates were retained. The breastfed-infant extension is modellib('Dickinson_2021_dolutegravir_motherinfant')."
  reference <- paste(
    "Dickinson L, Walimbwa S, Singh Y, Kaboggoza J, Kintu K, Sihlangu M,",
    "Coombs JA, Malaba TR, Byamugisha J, Pertinez H, Amara A, Gini J, Else L,",
    "Heiberg C, Hodel EM, Reynolds H, Myer L, Waitt C, Khoo S, Lamorde M,",
    "Orrell C; DolPHIN-1 Study Group.",
    "Infant exposure to dolutegravir through placental and breast milk",
    "transfer: a population pharmacokinetic analysis of DolPHIN-1.",
    "Clin Infect Dis. 2021;73(5):e1200-e1207. doi:10.1093/cid/ciaa1861.",
    "Structure from the Results text and Figure 3; parameter values from",
    "Table 2 and its footnote b; estimation and BLQ methods from the",
    "Supplementary Material.",
    sep = " "
  )
  vignette <- "Dickinson_2021_dolutegravir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The fetal state is not in the canonical compartment register; it follows
  # the Fauchet_2015_lopinavir_placental.R precedent for a cord-blood
  # compartment that carries a concentration and does not deplete the mother.
  paper_specific_compartments <- c("fetal")

  covariateData <- list(
    OCC = list(
      description = "Integer-valued occasion indicator for between-occasion variability on maternal CL/F",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Dickinson 2021 Table 2 footnote b reports a single interoccasion variability of 21.0% on CL/F without defining the occasions. The DolPHIN-1 sampling visits were the intensive third-trimester profile (day 14 of treatment), the delivery sample, and the intensive postpartum profile (Methods, 'Study Design and Pharmacokinetic Sampling'), so three occasion slots are encoded; occasions 2 and 3 fix their variance to the occasion-1 value (the NONMEM $OMEGA BLOCK(1) SAME idiom). For a single-occasion record pass OCC = 1.",
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    PREG = list(
      description = "Pregnancy status (third trimester vs postpartum)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (postpartum)",
      notes = "Screened on the maternal model (Supplementary Material, 'Covariate analysis') but not retained: there was no significant difference in CL/F between the third trimester and postpartum (Discussion).",
      source_name = "pregnancy (third trimester vs postpartum)"
    ),
    WT = list(
      description = "Maternal body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on the maternal model but not retained (Results, 'Population Pharmacokinetic Modeling': no covariates had a significant association with Vc/F).",
      source_name = "bodyweight"
    ),
    AGE = list(
      description = "Maternal age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on the maternal model but not retained.",
      source_name = "age"
    ),
    TPP = list(
      description = "Time postpartum",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "The paper tested time postpartum in DAYS (the canonical TPP is in weeks). Days postpartum on CL/F was the only univariable-significant covariate (dOFV = -7.34) but did not remain after backward elimination (Results, 'Population Pharmacokinetic Modeling'). Mode of delivery (vaginal vs cesarean) was also screened and not retained.",
      source_name = "days postpartum"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "dolutegravir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dolutegravir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dolutegravir", units = "mg", specimen = "plasma", verified = TRUE),
    fetal = list(analyte = "dolutegravir", units = "mg/L", specimen = "whole blood", verified = FALSE),
    milk = list(analyte = "dolutegravir", units = "mg", specimen = "milk", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 28,
    n_studies = 1,
    age_range = "19-42 years",
    age_median = "27 years",
    weight_range = "44-160 kg",
    weight_median = "67 kg",
    sex_female_pct = 100,
    disease_state = "Treatment-naive pregnant women living with HIV-1, diagnosed late in pregnancy (28-36 weeks gestation), randomized to dolutegravir-based therapy; followed through delivery and up to 18 days postpartum while breastfeeding",
    dose_range = "Dolutegravir 50 mg orally once daily with a tenofovir disoproxil fumarate plus emtricitabine or lamivudine backbone, from the third trimester until a switch to efavirenz-based therapy 2-18 days postpartum",
    regions = "Uganda (Kampala; 14 women) and South Africa (Cape Town; 14 women)",
    n_observations = "533 maternal plasma (250 third trimester, 18 delivery, 265 postpartum), 16 umbilical cord and 80 breast-milk concentrations",
    notes = "DolPHIN-1 (NCT02245022). Baseline demographics from Table 1. Mode of delivery: 23 vaginal (82%), 3 cesarean (11%), 2 missing. Postpartum sampling interval median 7 days (range 2-18). LLQ 0.01 mg/L in all matrices; 39% of breast-milk samples were BLQ and were handled with the M3 method."
  )

  ini({
    # Maternal plasma (Table 2). The analytical-solution maternal model was
    # estimated with the Laplacian method in NONMEM 7.4 (Supplementary
    # Material). CL/F, Vc/F, Q/F and Vp/F are apparent (oral) values.
    lcl <- log(1.50)
    label("Apparent maternal clearance CL/F (L/h)") # Table 2 'CL/F' = 1.50 L/h (RSE 1.9%)
    lvc <- log(24.6)
    label("Apparent maternal central volume Vc/F (L)") # Table 2 'Vc/F' = 24.6 L (RSE 4.7%)
    lq <- log(0.0138)
    label("Apparent maternal intercompartmental clearance Q/F (L/h)") # Table 2 'Q/F' = 0.0138 L/h (RSE 3.7%)
    lvp <- log(2.01)
    label("Apparent maternal peripheral volume Vp/F (L)") # Table 2 'Vp/F' = 2.01 L (RSE 16.7%)
    lka <- log(0.75)
    label("Absorption rate constant ka (1/h)") # Table 2 'ka' = 0.75 1/h (RSE 4.9%)

    # Umbilical cord (Table 2). No IIV could be estimated on the transfer
    # rate constants (Results).
    lk_central_fetal <- log(2.81)
    label("Mother-to-fetus transfer rate constant kM-F (1/h)") # Table 2 'kM-F' = 2.81 1/h (RSE 3.6%)
    lk_fetal_central <- log(2.20)
    label("Fetus-to-mother transfer rate constant kF-M (1/h)") # Table 2 'kF-M' = 2.20 1/h (RSE 40.9%)

    # Breast milk (Table 2).
    lk_central_milk <- log(0.0027)
    label("Mother-to-breast-milk transfer rate constant kM-BM (1/h)") # Table 2 'kM-BM' = 0.0027 1/h (RSE 32.3%)
    lk_milk_central <- log(16.3)
    label("Breast-milk-to-mother transfer rate constant kBM-M (1/h)") # Table 2 'kBM-M' = 16.3 1/h (RSE 4.6%)
    lvmilk <- fixed(log(0.125))
    label("Breast-milk compartment volume VBM (L)") # Table 2 'VBM' = 0.125 L (Fixed); Results: fixed from refs 16 and 18 for identifiability

    # IIV (Table 2 footnote b), reported as percentages and converted with
    # omega^2 = log(1 + CV^2). The CL/F-Vc/F correlation of 0.88 gives the
    # covariance 0.88 * sqrt(0.020240 * 0.041965) = 0.025647.
    etalcl + etalvc ~ c(0.020240, 0.025647, 0.041965) # Table 2 footnote b: IIV CL/F 14.3% (RSE 6.0%), Vc/F 20.7% (RSE 23.7%); Table 2 'Correlation between CL/F and Vc/F' = 0.88

    # Interoccasion variability on CL/F (Table 2 footnote b): 21.0% (RSE 0.6%);
    # omega^2 = log(1 + 0.210^2) = 0.043155. One variance shared by every
    # occasion; occasions 2 and 3 are fixed to the occasion-1 value.
    etaiov_cl_1 ~ 0.043155
    etaiov_cl_2 ~ fix(0.043155) # SAME-equivalent: equal to the occasion-1 variance
    etaiov_cl_3 ~ fix(0.043155) # SAME-equivalent: equal to the occasion-1 variance

    # Residual error (Table 2 'Residual error, %'): proportional in every
    # matrix (Results, 'Proportional error models were used throughout').
    propSd <- 0.365
    label("Proportional residual error, maternal plasma (fraction)") # Table 2 'Maternal plasma' = 36.5% (RSE 0.3%)
    propSd_Cfetal <- 0.330
    label("Proportional residual error, umbilical cord (fraction)") # Table 2 'Umbilical cord' = 33.0% (RSE 0.4%)
    propSd_Cmilk <- 0.584
    label("Proportional residual error, breast milk (fraction)") # Table 2 'Breast milk' = 58.4% (RSE 0.5%)
  })

  model({
    # 1. Occasion multiplexing of the between-occasion eta on CL/F.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3

    # 2. Individual parameters. No covariate effects were retained.
    cl <- exp(lcl + etalcl + iov_cl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka)
    k_central_fetal <- exp(lk_central_fetal)
    k_fetal_central <- exp(lk_fetal_central)
    k_central_milk <- exp(lk_central_milk)
    k_milk_central <- exp(lk_milk_central)
    vmilk <- exp(lvmilk)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    Cc <- central / vc

    # 4. ODEs. The breast-milk compartment exchanges mass with the maternal
    #    central compartment (Figure 3); at the fitted values its steady-state
    #    content is kM-BM / kBM-M = 1.7e-4 of the central amount, so the
    #    exchange does not perceptibly alter maternal plasma. The fetal
    #    compartment has 'negligible volume (which did not alter the maternal
    #    compartment)' (Results), so it is written as a non-depleting link on
    #    the concentration scale: at steady state Cfetal / Cc = kM-F / kF-M =
    #    1.277, the published cord:maternal ratio of 1.279.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central +
      k21 * peripheral1 - k_central_milk * central + k_milk_central * milk
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(fetal) <- k_central_fetal * Cc - k_fetal_central * fetal
    d/dt(milk) <- k_central_milk * central - k_milk_central * milk

    # 5. Observations. Dose in mg and volumes in L give mg/L.
    Cfetal <- fetal
    Cmilk <- milk / vmilk

    Cc ~ prop(propSd)
    Cfetal ~ prop(propSd_Cfetal)
    Cmilk ~ prop(propSd_Cmilk)
  })
}
