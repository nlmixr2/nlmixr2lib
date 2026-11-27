Friedrich_2014_evacetrapib_hdl <- function() {
  description <- paste0(
    "Exposure-response model for the percent change from baseline in ",
    "high-density lipoprotein cholesterol (HDL-C) during once-daily oral ",
    "evacetrapib, a cholesteryl ester transfer protein (CETP) inhibitor, ",
    "given alone or with a statin to adults with hypercholesterolemia or ",
    "low HDL-C (Friedrich 2014, Eq. 1). Algebraic, with no ODE states: an ",
    "Emax function of the individual steady-state evacetrapib AUC, plus a ",
    "placebo term, approaches its plateau as 1 - exp(-K * time) from the ",
    "first dose. The AUC giving half the maximal effect rises with baseline ",
    "HDL-C, so patients with lower baseline HDL-C respond more at the same ",
    "exposure. The population placebo and statin effects and the Hill ",
    "coefficient were held at 0, 0 and 1 in the final model. Additive IIV on ",
    "Emax and the placebo term, exponential IIV on K, additive residual ",
    "error. Supply AUC_EVACETRAPIB from the companion PK model ",
    "Friedrich_2014_evacetrapib.R as 1000 x DOSE / CL/F."
  )
  reference <- paste0(
    "Friedrich S, Kastelein JJP, James D, Waterhouse T, Nissen SE, ",
    "Nicholls SJ, Krueger KA. The pharmacokinetics and ",
    "pharmacokinetic/pharmacodynamic relationships of evacetrapib ",
    "administered as monotherapy or in combination with statins. CPT ",
    "Pharmacometrics Syst Pharmacol. 2014;3(1):e94. doi:10.1038/psp.2013.70. ",
    "Final estimates: Table 2 and Supplementary Data ",
    "'Friedrich_HDL_model_code_and_output' (NONMEM 7.2 $PRED and final ",
    "estimates). PK layer: modellib('Friedrich_2014_evacetrapib')."
  )
  vignette <- "Friedrich_2014_evacetrapib"
  units <- list(
    time = "h",
    dosing = "n/a (no dosing events; evacetrapib exposure is supplied as the AUC_EVACETRAPIB covariate in ng*h/mL)",
    concentration = "percent (hdlcfb is the percentage change from baseline in HDL-C)"
  )

  covariateData <- list(
    AUC_EVACETRAPIB = list(
      description = "Individual steady-state evacetrapib AUC over the 24 h once-daily dosing interval",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Time-fixed per subject. The authors fixed each patient's post hoc ",
        "AUC from the final PK model in the analysis dataset (Methods). ",
        "Compute it from modellib('Friedrich_2014_evacetrapib') as ",
        "1000 * DOSE / cl, with DOSE in mg and the individual cl in L/h. ",
        "Set to 0 for placebo and statin-only patients. Geometric mean ",
        "values in the study: 2,300, 5,900 and 19,700 ng*h/mL at 30, 100 ",
        "and 500 mg monotherapy."
      ),
      source_name = "AUC"
    ),
    HDLC = list(
      description = "Baseline HDL-C",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Time-fixed baseline value. Enters the log EAUC50 multiplicatively, ",
        "log(EAUC50) = lauc50 * exp(e_hdlc_lauc50 * (HDLC - 52.2)), ",
        "52.2 mg/dL being the median baseline HDL-C (Table 2 footnote f; ",
        "NONMEM $PRED). Lower baseline HDL-C gives a lower EAUC50 and so a ",
        "larger percent rise at the same AUC."
      ),
      source_name = "bHDL (G43b in the NONMEM $PRED)"
    ),
    CONMED_STATIN = list(
      description = "Statin coadministered, 1 = atorvastatin 20 mg, simvastatin 40 mg or rosuvastatin 10 mg once daily",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no statin)",
      notes = paste0(
        "1 for every patient in a statin arm, with or without evacetrapib ",
        "(the NONMEM $PRED sets I1 = 1 for treatment groups above 4). It ",
        "multiplies the statin effect e_conmed_statin_hdl, which the final ",
        "model holds at 0, so the column has no effect on predictions; it is ",
        "kept so that the model reproduces the authors' code."
      ),
      source_name = "I1 (from GRP)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on the HDL-C parameters (Methods) and not retained."
    ),
    SEXF = list(
      description = "Sex, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Screened as gender (Methods) and not retained."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Methods) and not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened (Methods) and not retained."
    ),
    TRIG = list(
      description = "Baseline triglycerides",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened (Methods) and not retained in the HDL-C model."
    ),
    LDLC = list(
      description = "Baseline LDL-C",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened (Methods) and not retained in the HDL-C model."
    ),
    APOA1 = list(
      description = "Baseline apolipoprotein A-I",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened (Methods) and not retained in the HDL-C model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 391L,
    n_observations = 1882L,
    n_studies = 1L,
    age_range = "28-83 years (whole randomised study)",
    age_mean = "58.3 years (whole randomised study)",
    weight_range = "44-164 kg (whole randomised study)",
    weight_mean = "83.6 kg (whole randomised study)",
    sex_female_pct = 56.0,
    race_ethnicity = c(White = 92.9),
    disease_state = paste0(
      "Adults with hypercholesterolemia or low HDL-C; mean baseline HDL-C ",
      "55.1 mg/dL (median 52.2 mg/dL)."
    ),
    dose_range = paste0(
      "Placebo; evacetrapib 30, 100 or 500 mg once daily; atorvastatin ",
      "20 mg, simvastatin 40 mg or rosuvastatin 10 mg once daily alone or ",
      "with evacetrapib 100 mg once daily; 12 weeks."
    ),
    notes = paste0(
      "HDL-C sampled at weeks 2, 4, 8 and 12; 1,882 observations from 391 ",
      "patients. The endpoint is the percent change from each patient's own ",
      "baseline HDL-C."
    )
  )

  ini({
    # Supplementary Data 'Friedrich_HDL_model_code_and_output' (NONMEM 7.2
    # $PRED, FOCE-I) carries the final estimates; Table 2 prints the same
    # values on the natural scale.
    emax <- 177; label("Maximal percent change from baseline in HDL-C at steady state (%)") # Table 2 'Emax' = 177 (%SEE 6.55); NONMEM THETA(1) 177

    # The NONMEM code estimates log(EAUC50): LY = EMAX*AUC/(EXP(EAUC50) +
    # AUC) with EAUC50 = THETA(2)*EXP(THETA(7)*(bHDL - 52.2)), so the
    # baseline-HDL effect multiplies the LOG of EAUC50. Table 2 footnote f
    # writes it on the natural scale instead (5,380 * exp(0.00317 * (bHDL -
    # 52.2))); only the code's form reproduces Figure 2 and the Discussion's
    # 110% -> 127% example (see the vignette).
    lauc50 <- 8.59; label("Log of the steady-state AUC giving half the maximal effect at baseline HDL-C 52.2 mg/dL (log ng*h/mL)") # NONMEM THETA(2) 8.59 (exp = 5,379); Table 2 'EAUC50' = 5,380 h*ng/ml (%SEE 1.54)
    e_hdlc_lauc50 <- 0.00317; label("Exponential effect of baseline HDL-C on log EAUC50 (1/(mg/dL))") # Table 2 'Effect of baseline HDL on EAUC50' = 0.00317 (%SEE 9.97); NONMEM THETA(7)
    lhill <- fixed(log(1)); label("Hill coefficient of the AUC Emax function (unitless)") # NONMEM THETA(3) GAMMA = 1 (held constant); Results: including GAM did not improve the fit
    plac <- fixed(0); label("Typical placebo percent change from baseline in HDL-C (%)") # Table 2 'Placebo effect' = 0 (Fixed); NONMEM THETA(4) 0
    e_conmed_statin_hdl <- fixed(0); label("Statin percent change from baseline in HDL-C (%)") # NONMEM THETA(6) STAT = 0 (held constant); Results: a preliminary estimate of 2.28 (%SEE 83.3) did not improve the fit
    lkhdl <- log(0.00684); label("First-order rate constant of the approach to the HDL-C plateau (1/h)") # Table 2 'Time course constant (K)' = 0.00684 hour-1 (%SEE 7.43); NONMEM THETA(5)

    etaemax ~ 3680 # NONMEM OMEGA(1) 3680, additive on Emax; Table 2 IIV 34.3% = 100*sqrt(3680)/177
    etalkhdl ~ 0.281 # NONMEM OMEGA(2) 0.281, exponential on K; Table 2 IIV 57% = 100*sqrt(exp(0.281)-1)
    etaplac ~ 59.7 # NONMEM OMEGA(3) 59.7, additive on the placebo term; Table 2 IIV 7.7 = sqrt(59.7) percentage points

    addSd <- 14.18; label("Additive residual error (percentage points)") # NONMEM SIGMA 201 (variance) -> sqrt(201) = 14.18; Table 2 'Residual error (additive, %)' = 14.2 (%SEE 11.2)
  })
  model({
    emax_i <- emax + etaemax
    plac_i <- plac + etaplac
    khdl <- exp(lkhdl + etalkhdl)
    hill <- exp(lhill)

    lauc50_i <- lauc50 * exp(e_hdlc_lauc50 * (HDLC - 52.2))
    auc50 <- exp(lauc50_i)

    drug <- emax_i * AUC_EVACETRAPIB^hill / (auc50^hill + AUC_EVACETRAPIB^hill)
    hdlss <- plac_i + e_conmed_statin_hdl * CONMED_STATIN + drug

    # time is measured from the first dose of study treatment.
    hdlcfb <- hdlss * (1 - exp(-khdl * time))
    hdlcfb ~ add(addSd)
  })
}
