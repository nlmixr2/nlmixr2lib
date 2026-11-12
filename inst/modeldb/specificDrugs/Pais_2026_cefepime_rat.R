Pais_2026_cefepime_rat <- function() {
  description <- paste(
    "Preclinical (rat). Three-compartment population PK model of intravenous cefepime",
    "in male Sprague-Dawley rats with folic-acid-induced acute kidney injury: a plasma",
    "(central) compartment exchanging first-order with two sampled brain-tissue",
    "compartments, the cerebral cortex and the hippocampus, with first-order elimination",
    "from plasma. Body weight scales the elimination rate constant (exponent 0.75) and",
    "the central volume (exponent 1), both centred on the cohort weight. Fit in Monolix",
    "2024R1 (SAEM). The paper reports no between-animal variances and no residual-error",
    "magnitudes, so the model is typical-value only. Companion logistic exposure-seizure",
    "models from the same paper: Pais_2026_cefepime_seizure_*_rat.",
    sep = " "
  )
  reference <- paste(
    "Pais GM, Lesnicki E, Marianski S, Valdez K, Gibson Z, Christopher J, Lepard K,",
    "Gilchrist A, Mizuno T, Scheetz MH.",
    "Elevated hippocampal and cortical cefepime concentrations correlate with seizures",
    "in an acute kidney injury rat model.",
    "Antimicrob Agents Chemother. 2026;70(4):e01005-25. doi:10.1128/aac.01005-25.",
    sep = " "
  )
  vignette <- "Pais_2026_cefepime_rat"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters as log(WT / 0.282) on log(kel) (coefficient 0.75) and log(vc) (coefficient 1).",
        "The paper centres on 'the median rat weight of the data set' without printing it;",
        "0.282 kg is the reported cohort mean (282 +/- 17.0 g) and is used as the reference.",
        "Supply in kg (a 282 g rat is WT = 0.282)."
      ),
      source_name = "WT (Table 2 'beta_kel_logtWT', 'beta_Vcentral_logtWT')"
    )
  )

  compartmentData <- list(
    central = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE),
    brain_cortex = list(analyte = "cefepime", units = "mg", specimen = "tissue", verified = TRUE),
    brain_hippocampus = list(analyte = "cefepime", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley; male)",
    n_subjects = 18L,
    n_studies = 1L,
    n_observations = "127 plasma (of 129; two >3 SD outliers excluded), 18 cerebral cortex and 18 hippocampus concentrations (one brain harvest per rat at euthanasia)",
    age_range = "approximately 9-10 weeks",
    weight_range = "mean 282 g (SD 17.0 g)",
    sex_female_pct = 0,
    disease_state = paste(
      "Folic-acid-induced acute kidney injury (250 mg/kg intraperitoneally on day 1, then",
      "100 mg/kg daily in the multiple-dose arm), which prolongs the plasma cefepime",
      "half-life about 10-fold. Seizures were scored on a modified Racine scale (0-7)."
    ),
    dose_range = paste(
      "Cefepime 1,593 mg/kg once (single folic-acid dose, serial sacrifice on day 1) or",
      "1,250 mg/kg once daily for up to 5 days, intravenously via a jugular catheter 30 min",
      "after folic acid."
    ),
    regions = "United States (Midwestern University, Downers Grove, Illinois)",
    notes = paste(
      "Plasma sampled at 0, 15, 30 and 120 min after each daily dose; cortex and hippocampus",
      "harvested once at euthanasia (30 or 120 min on day 1, or 120 min on day 3 or 5) after",
      "saline perfusion. Brain concentrations are total tissue concentrations in ug/g, treated",
      "as mg/L by the source model. Concentrations are total (not unbound)."
    )
  )

  ini({
    # Pais 2026 Table 2 'Median parameter values from the final model'. Monolix
    # parameterises each typical value directly, so every structural value is the
    # printed median. The ODEs are Figure 3 / Methods 'Ordinary differential
    # equations'; see model() for the two typographical errors in the printed
    # equations and how the paper's own exposures settle them.
    lkel <- log(0.14); label("Elimination rate constant from plasma, kel (1/h)") # Table 2 kel = 0.14 1/h (RSE 38.7%)
    lvc <- log(0.12); label("Central (plasma) volume, Vcentral (L)") # Table 2 Vcentral = 0.12 L (RSE 6.85%)
    lv_brain_cortex <- log(0.0039); label("Cerebral cortex volume, Vcortex (L)") # Table 2 Vcortex = 0.0039 L (RSE 18.6%)
    lv_brain_hippocampus <- log(0.00039); label("Hippocampus volume, Vhippocampus (L)") # Table 2 Vhippocampus = 0.00039 L (RSE 14.8%)

    # k12 has no RSE in Table 2; Discussion: 'the transfer rate constant from the
    # central compartment to the cortex was fixed to reduce RSE%'. Table 1 row
    # 'Three-compartment (fixed k12 and corr Vhippocampus kel)' is the selected model.
    lk_central_brain_cortex <- fixed(log(0.00013)); label("Plasma-to-cortex rate constant, k12 (1/h)") # Table 2 k12 = 0.00013 1/h, no RSE (fixed)
    lk_brain_cortex_central <- log(0.17); label("Cortex-to-plasma rate constant, k21 (1/h)") # Table 2 k21 = 0.17 1/h (RSE 18.5%)
    lk_central_brain_hippocampus <- log(0.000015); label("Plasma-to-hippocampus rate constant, k13 (1/h)") # Table 2 k13 = 0.000015 1/h (RSE 24.9%)
    lk_brain_hippocampus_central <- log(0.15); label("Hippocampus-to-plasma rate constant, k31 (1/h)") # Table 2 k31 = 0.15 1/h (RSE 31.2%)

    # Weight effects: Monolix log-transformed covariate logtWT = log(WT / WTref).
    # Neither coefficient carries an RSE, so both were fixed.
    e_wt_kel <- fixed(0.75); label("Exponent of WT/0.282 kg on kel (unitless)") # Table 2 beta_kel_logtWT = 0.75, no RSE (fixed)
    e_wt_vc <- fixed(1); label("Exponent of WT/0.282 kg on Vcentral (unitless)") # Table 2 beta_Vcentral_logtWT = 1, no RSE (fixed)

    # Between-animal variability: Table 2 reports only the correlation
    # corr_kel_Vhippocampus = 0.92 (RSE 20.8%) between the kel and Vhippocampus
    # random effects; no standard deviation of any random effect is printed, so
    # no eta is encoded (a correlation without variances cannot be carried).

    # Residual error: neither the error model nor its magnitudes are reported.
    # Both a proportional and an additive term are carried per output at zero.
    propSd <- fixed(0); label("Proportional residual SD, plasma (fraction; not reported)") # not reported in the source
    addSd <- fixed(0); label("Additive residual SD, plasma (mg/L; not reported)") # not reported in the source
    propSd_Cbrain_cortex <- fixed(0); label("Proportional residual SD, cortex (fraction; not reported)") # not reported in the source
    addSd_Cbrain_cortex <- fixed(0); label("Additive residual SD, cortex (mg/L; not reported)") # not reported in the source
    propSd_Cbrain_hippocampus <- fixed(0); label("Proportional residual SD, hippocampus (fraction; not reported)") # not reported in the source
    addSd_Cbrain_hippocampus <- fixed(0); label("Additive residual SD, hippocampus (mg/L; not reported)") # not reported in the source
  })

  model({
    # Weight scaling centred on the cohort weight (0.282 kg; see covariateData).
    lwtn <- log(WT / 0.282)
    kel <- exp(lkel + e_wt_kel * lwtn)
    vc <- exp(lvc + e_wt_vc * lwtn)
    v_brain_cortex <- exp(lv_brain_cortex)
    v_brain_hippocampus <- exp(lv_brain_hippocampus)
    k_central_brain_cortex <- exp(lk_central_brain_cortex)
    k_brain_cortex_central <- exp(lk_brain_cortex_central)
    k_central_brain_hippocampus <- exp(lk_central_brain_hippocampus)
    k_brain_hippocampus_central <- exp(lk_brain_hippocampus_central)

    # Figure 3 schematic: central <-> cortex (k12 / k21), central <-> hippocampus
    # (k13 / k31), kel out of central. The printed cortex equation reads
    # 'ddt_Acortex = k12*Acentral - (Acortex*k12)' and the printed output reads
    # 'Ccortex = Acentral/Vcortex'. Both are typographical: k21 is estimated
    # (Table 2) yet would appear nowhere, and Acentral/Vcortex would put cortex
    # concentrations about 30-fold above plasma, whereas Results report a cortex
    # AUC0-24 of 565 against a plasma AUC0-24 of 28,314 mg*h/L (ratio 0.020).
    # The schematic's reading gives the ratio k12*Vc/(k21*Vcortex) = 0.024, in
    # line with the paper; the hippocampus equations are printed consistently.
    d/dt(central) <- -(kel + k_central_brain_cortex + k_central_brain_hippocampus) * central +
      k_brain_cortex_central * brain_cortex +
      k_brain_hippocampus_central * brain_hippocampus
    d/dt(brain_cortex) <- k_central_brain_cortex * central - k_brain_cortex_central * brain_cortex
    d/dt(brain_hippocampus) <- k_central_brain_hippocampus * central - k_brain_hippocampus_central * brain_hippocampus

    Cc <- central / vc
    Cbrain_cortex <- brain_cortex / v_brain_cortex
    Cbrain_hippocampus <- brain_hippocampus / v_brain_hippocampus

    Cc ~ add(addSd) + prop(propSd)
    Cbrain_cortex ~ add(addSd_Cbrain_cortex) + prop(propSd_Cbrain_cortex)
    Cbrain_hippocampus ~ add(addSd_Cbrain_hippocampus) + prop(propSd_Cbrain_hippocampus)
  })
}
