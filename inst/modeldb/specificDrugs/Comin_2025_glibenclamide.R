Comin_2025_glibenclamide <- function() {
  description <- paste(
    "Seven-compartment semi-physiological population PK model for",
    "intravenous microdose [11C]glibenclamide (glyburide) in healthy adult",
    "volunteers, fitted to whole-body dynamic PET concentrations in arterial",
    "blood, liver, kidneys, spleen and pancreas over the first 40 min after",
    "injection. Blood (central) exchanges first-order with a liver",
    "compartment and a kidney compartment (each with a catenary deep",
    "sub-compartment, the two sub-compartments sharing one organ volume) and",
    "with one-compartment spleen and pancreas; first-order elimination from",
    "blood; the dose enters blood as a zero-order input of estimated",
    "duration. A single intravenous rifampicin dose (OATP inhibition)",
    "multiplies the liver-to-blood rate constant by exp(0.96) = 2.6.",
    "Inter-occasion variability on that rate constant is fixed (SD 0.01)."
  )
  reference <- paste(
    "Comin L, Marie S, Ursino M, Zohar S, Tournier N, Comets E. Modeling",
    "Whole-Body Dynamic PET Microdosing Data to Predict the Whole-Body",
    "Pharmacokinetics of Glyburide in Humans. Clin Pharmacokinet.",
    "2025;64(11):1709-1722. doi:10.1007/s40262-025-01562-9 (PMC12618397)."
  )
  vignette <- "Comin_2025_glibenclamide"

  # Every state is an amount in ug (the injected radioactivity was converted
  # to a glibenclamide mass with the specific radioactivity and the molar mass
  # 494 ug/umol, Methods eq. 1). Volumes are in L, so concentrations are ug/L,
  # the unit of the additive residual errors in Table 2.
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ug/L"
  )

  covariateData <- list(
    CONMED_RIFAMPICIN_SD = list(
      description = paste(
        "1 = [11C]glibenclamide injected immediately after a single",
        "intravenous rifampicin infusion (9 mg/kg over 30 min; second PET",
        "occasion), 0 = [11C]glibenclamide alone (baseline occasion)"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Acute OATP inhibition phase. Enters as",
        "k_liver_central = k_liver_central,pop * exp(0.96 * CONMED_RIFAMPICIN_SD)",
        "(Table 2 legend equation), a 2.6-fold faster return from the liver",
        "to blood. The indicator applies to the studied design only",
        "(single 9 mg/kg IV rifampicin given immediately before the",
        "radiotracer)."
      ),
      source_name = "Rif"
    ),
    OCC = list(
      description = "PET acquisition occasion (1 = baseline, 2 = after rifampicin)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Selects the occasion-level random effect on the liver-to-blood rate",
        "constant (kappa_io in Methods eq. 2). In the source study occasion 2",
        "always coincides with CONMED_RIFAMPICIN_SD = 1. Records with OCC",
        "outside 1-2 carry no inter-occasion effect."
      ),
      source_name = "occasion (o)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "glibenclamide", units = "ug", specimen = "whole blood", verified = TRUE),
    liver = list(analyte = "glibenclamide", units = "ug", specimen = "tissue", verified = TRUE),
    liver_deep = list(analyte = "glibenclamide", units = "ug", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "glibenclamide", units = "ug", specimen = "tissue", verified = TRUE),
    kidney_deep = list(analyte = "glibenclamide", units = "ug", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "glibenclamide", units = "ug", specimen = "tissue", verified = TRUE),
    pancreas = list(analyte = "glibenclamide", units = "ug", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 15,
    n_studies = 1,
    age_range = "20-65 years",
    weight_range = "53-88 kg",
    sex_female_pct = 26.7,
    disease_state = "healthy volunteers",
    dose_range = paste(
      "Intravenous bolus microdose of [11C]glibenclamide (12.5 +/- 6.2 ug),",
      "alone and, in 10 of 16 enrolled subjects, again at least 3 h later",
      "immediately after a 30-min IV infusion of rifampicin 9 mg/kg"
    ),
    regions = "France",
    notes = paste(
      "16 healthy volunteers enrolled in three groups (Table 1): men > 50",
      "years (n = 4, median age 59, weight 81 kg), women > 50 years (n = 5,",
      "60 years, 66 kg) and men < 30 years (n = 7, 25 years, 73 kg). One",
      "woman homozygous for SLCO1B1 c.521T>C (CC) was excluded from model",
      "building, so the final model was fitted to 15 subjects (4 women);",
      "Supplementary Table S1 refits it with that subject included. One man",
      "had been splenectomised and contributed no spleen data. Whole-body",
      "dynamic PET-MR (GE Signa) over 40 min; arterial blood from a left",
      "ventricle + aorta image-derived input. Total radioactivity was",
      "attributed to the parent drug (radiometabolites < 10% by the end of",
      "acquisition). EudraCT 2017-001703-69."
    )
  )

  ini({
    # Blood (central) and elimination - Table 2
    lvc <- log(3.6); label("Blood (central) volume VB (L)") # Table 2 'VB,pop' 3.6 L (RSE 8.3%)
    lkel <- log(0.79); label("Elimination rate constant from blood ke (1/h)") # Table 2 'ke,pop' 0.79 1/h (RSE 15%)
    ld1 <- log(0.015); label("Duration of the zero-order input into blood Tk0 (h)") # Table 2 'Tk0,pop' 0.015 h (RSE 13%)

    # Liver - Table 2
    lv_liver <- log(1.1); label("Liver volume VL, shared by both liver sub-compartments (L)") # Table 2 'VL,pop' 1.1 L (RSE 15%)
    lk_liver_central <- log(21); label("Liver (first sub-compartment) to blood rate constant kL1B (1/h)") # Table 2 'kL1B,pop' 21 1/h (RSE 8.7%)
    lk_central_liver <- log(9.0); label("Blood to liver (first sub-compartment) rate constant kBL1 (1/h)") # Table 2 'kBL1,pop' 9.0 1/h (RSE 14%)
    lk_liver_liver_deep <- log(8.2); label("Liver first to second sub-compartment rate constant kL1L2 (1/h)") # Table 2 'kL1L2,pop' 8.2 1/h (RSE 9.6%)
    lk_liver_deep_liver <- log(4.5); label("Liver second to first sub-compartment rate constant kL2L1 (1/h)") # Table 2 'kL2L1,pop' 4.5 1/h (RSE 8.9%)
    e_conmed_rifampicin_sd_k_liver_central <- 0.96; label("Rifampicin effect on log kL1B (unitless)") # Table 2 'beta L1B,Rif' 0.96 (RSE 1.6%)

    # Kidneys - Table 2
    lv_kidney <- log(1.6); label("Kidney volume VK, shared by both kidney sub-compartments (L)") # Table 2 'VK,pop' 1.6 L (RSE 6.1%)
    lk_kidney_central <- log(179); label("Kidney (first sub-compartment) to blood rate constant kK1B (1/h)") # Table 2 'kK1B,pop' 179 1/h (RSE 10%)
    lk_central_kidney <- log(62); label("Blood to kidney (first sub-compartment) rate constant kBK1 (1/h)") # Table 2 'kBK1,pop' 62 1/h (RSE 13%)
    lk_kidney_kidney_deep <- log(19); label("Kidney first to second sub-compartment rate constant kK1K2 (1/h)") # Table 2 'kK1K2,pop' 19 1/h (RSE 9.8%)
    lk_kidney_deep_kidney <- log(8.2); label("Kidney second to first sub-compartment rate constant kK2K1 (1/h)") # Table 2 'kK2K1,pop' 8.2 1/h (RSE 6.3%)

    # Spleen - Table 2 (VS and kBS poorly identified, RSE 2.7e5% and 3.0e5%)
    lv_spleen <- log(0.34); label("Spleen volume VS (L)") # Table 2 'VS,pop' 0.34 L (RSE 2.7e+5%)
    lk_central_spleen <- log(5.9); label("Blood to spleen rate constant kBS (1/h)") # Table 2 'kBS,pop' 5.9 1/h (RSE 3.0e+5%)
    lk_spleen_central <- log(131); label("Spleen to blood rate constant kSB (1/h)") # Table 2 'kSB,pop' 131 1/h (RSE 6.0%)

    # Pancreas - Table 2 (VP and kBP poorly identified, RSE 3.7e5% and 4.0e5%)
    lv_pancreas <- log(0.23); label("Pancreas volume VP (L)") # Table 2 'VP,pop' 0.23 L (RSE 3.7e+5%)
    lk_central_pancreas <- log(5.9); label("Blood to pancreas rate constant kBP (1/h)") # Table 2 'kBP,pop' 5.9 1/h (RSE 4.0e+5%)
    lk_pancreas_central <- log(138); label("Pancreas to blood rate constant kPB (1/h)") # Table 2 'kPB,pop' 138 1/h (RSE 11%)

    # IIV - Table 2 reports omega as the SD of the log-scale random effect;
    # variances below are omega^2. No correlations were retained.
    etalvc ~ 0.0324 # Table 2 omega VB = 0.18 (CV 18%)
    etalkel ~ 0.1225 # Table 2 omega ke = 0.35 (CV 36%)
    etald1 ~ 0.2116 # Table 2 omega Tk0 = 0.46 (CV 49%)
    etalv_liver ~ 0.2916 # Table 2 omega VL = 0.54 (CV 58%)
    etalk_liver_central ~ 0.0484 # Table 2 omega L1B = 0.22 (CV 22%)
    etalk_central_liver ~ 0.16 # Table 2 omega BL1 = 0.40 (CV 42%)
    etalk_liver_deep_liver ~ 0.0484 # Table 2 omega L2L1 = 0.22 (CV 22%)
    etalv_kidney ~ 0.008281 # Table 2 omega VK = 0.091 (CV 9.1%)
    etalk_central_kidney ~ 0.0676 # Table 2 omega BK1 = 0.26 (CV 27%)
    etalk_kidney_kidney_deep ~ 0.0225 # Table 2 omega K1K2 = 0.15 (CV 15%)
    etalk_kidney_deep_kidney ~ 0.0144 # Table 2 omega K2K1 = 0.12 (CV 12%)
    etalk_central_spleen ~ 0.04 # Table 2 omega BS = 0.20 (CV 20%)
    etalk_central_pancreas ~ 0.1444 # Table 2 omega BP = 0.38 (CV 39%)
    etalk_pancreas_central ~ 0.09 # Table 2 omega PB = 0.30 (CV 31%)

    # IOV on kL1B only, held at SD 0.01 (variance 1e-4) - Results and Table 2 legend
    etaiov_lk_liver_central_1 ~ fixed(0.0001) # Table 2 legend: IOV gamma kL1B = 0.01 (SD), variance 1e-4
    etaiov_lk_liver_central_2 ~ fixed(0.0001) # same IOV variance on occasion 2

    # Residual error, Monolix combined1 g = a + b*f (Methods); liver and
    # kidney constant error only - Table 2
    addSd <- 0.062; label("Blood additive residual error aB (ug/L)") # Table 2 'aB' 0.062 ug/L (RSE 6.2%)
    propSd <- 0.26; label("Blood proportional residual error bB (fraction)") # Table 2 'bB' 0.26 (RSE 0.72%)
    addSd_Cliver <- 0.12; label("Liver additive residual error aL (ug/L)") # Table 2 'aL' 0.12 ug/L (RSE 1.7%)
    addSd_Ckidney <- 0.23; label("Kidney additive residual error aK (ug/L)") # Table 2 'aK' 0.23 ug/L (RSE 0.41%)
    addSd_Cspleen <- 0.038; label("Spleen additive residual error aS (ug/L)") # Table 2 'aS' 0.038 ug/L (RSE 9.0%)
    propSd_Cspleen <- 0.24; label("Spleen proportional residual error bS (fraction)") # Table 2 'bS' 0.24 (RSE 1.0%)
    addSd_Cpancreas <- 0.026; label("Pancreas additive residual error aP (ug/L)") # Table 2 'aP' 0.026 ug/L (RSE 11%)
    propSd_Cpancreas <- 0.19; label("Pancreas proportional residual error bP (fraction)") # Table 2 'bP' 0.19 (RSE 4.8%)
  })

  model({
    # Occasion indicators for the inter-occasion random effect
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    # Individual parameters (log-normal; Methods eq. 2 and Table 2 legend)
    vc <- exp(lvc + etalvc)
    kel <- exp(lkel + etalkel)
    d1 <- exp(ld1 + etald1)

    v_liver <- exp(lv_liver + etalv_liver)
    k_liver_central <- exp(
      lk_liver_central + e_conmed_rifampicin_sd_k_liver_central * CONMED_RIFAMPICIN_SD +
        etalk_liver_central + oc1 * etaiov_lk_liver_central_1 + oc2 * etaiov_lk_liver_central_2
    )
    k_central_liver <- exp(lk_central_liver + etalk_central_liver)
    k_liver_liver_deep <- exp(lk_liver_liver_deep)
    k_liver_deep_liver <- exp(lk_liver_deep_liver + etalk_liver_deep_liver)

    v_kidney <- exp(lv_kidney + etalv_kidney)
    k_kidney_central <- exp(lk_kidney_central)
    k_central_kidney <- exp(lk_central_kidney + etalk_central_kidney)
    k_kidney_kidney_deep <- exp(lk_kidney_kidney_deep + etalk_kidney_kidney_deep)
    k_kidney_deep_kidney <- exp(lk_kidney_deep_kidney + etalk_kidney_deep_kidney)

    v_spleen <- exp(lv_spleen)
    k_central_spleen <- exp(lk_central_spleen + etalk_central_spleen)
    k_spleen_central <- exp(lk_spleen_central)

    v_pancreas <- exp(lv_pancreas)
    k_central_pancreas <- exp(lk_central_pancreas + etalk_central_pancreas)
    k_pancreas_central <- exp(lk_pancreas_central + etalk_pancreas_central)

    # Amount ODEs (Results eq. 5); the dose is a zero-order input of
    # duration Tk0 into blood (dose records need rate = -2)
    d/dt(central) <- k_liver_central * liver + k_kidney_central * kidney +
      k_spleen_central * spleen + k_pancreas_central * pancreas -
      (k_central_liver + k_central_kidney + k_central_spleen + k_central_pancreas + kel) * central
    d/dt(liver) <- k_central_liver * central + k_liver_deep_liver * liver_deep -
      (k_liver_central + k_liver_liver_deep) * liver
    d/dt(liver_deep) <- k_liver_liver_deep * liver - k_liver_deep_liver * liver_deep
    d/dt(kidney) <- k_central_kidney * central + k_kidney_deep_kidney * kidney_deep -
      (k_kidney_central + k_kidney_kidney_deep) * kidney
    d/dt(kidney_deep) <- k_kidney_kidney_deep * kidney - k_kidney_deep_kidney * kidney_deep
    d/dt(spleen) <- k_central_spleen * central - k_spleen_central * spleen
    d/dt(pancreas) <- k_central_pancreas * central - k_pancreas_central * pancreas

    dur(central) <- d1

    # Organ concentrations (Results eq. 5): each two-sub-compartment organ
    # is observed as the sum of its sub-compartments over one organ volume
    Cc <- central / vc
    Cliver <- (liver + liver_deep) / v_liver
    Ckidney <- (kidney + kidney_deep) / v_kidney
    Cspleen <- spleen / v_spleen
    Cpancreas <- pancreas / v_pancreas

    Cc ~ add(addSd) + prop(propSd) + combined1()
    Cliver ~ add(addSd_Cliver)
    Ckidney ~ add(addSd_Ckidney)
    Cspleen ~ add(addSd_Cspleen) + prop(propSd_Cspleen) + combined1()
    Cpancreas ~ add(addSd_Cpancreas) + prop(propSd_Cpancreas) + combined1()
  })
}
