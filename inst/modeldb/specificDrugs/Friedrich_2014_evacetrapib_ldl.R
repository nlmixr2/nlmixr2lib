Friedrich_2014_evacetrapib_ldl <- function() {
  description <- paste0(
    "Exposure-response model for the percent change from baseline in ",
    "low-density lipoprotein cholesterol (LDL-C) during once-daily oral ",
    "evacetrapib, a cholesteryl ester transfer protein (CETP) inhibitor, ",
    "given alone or with atorvastatin, simvastatin or rosuvastatin to adults ",
    "with hypercholesterolemia or low HDL-C (Friedrich 2014, Eqs. 2 and 3). ",
    "Algebraic and time-independent (the response is taken to be at steady ",
    "state at every on-treatment visit): an Emax function of the individual ",
    "steady-state evacetrapib AUC whose maximum scales with baseline ",
    "apolipoprotein A-I, a placebo term that depends on baseline LDL-C and ",
    "triglycerides, a statin effect, and an interaction term for the ",
    "combination arms. The interaction coefficient of -0.997 makes the ",
    "combination effect the product of the two independent fractional ",
    "reductions, i.e. evacetrapib and statins act independently. Additive ",
    "IIV on the placebo term and additive residual error. Supply ",
    "AUC_EVACETRAPIB from the companion PK model Friedrich_2014_evacetrapib.R ",
    "as 1000 x DOSE / CL/F."
  )
  reference <- paste0(
    "Friedrich S, Kastelein JJP, James D, Waterhouse T, Nissen SE, ",
    "Nicholls SJ, Krueger KA. The pharmacokinetics and ",
    "pharmacokinetic/pharmacodynamic relationships of evacetrapib ",
    "administered as monotherapy or in combination with statins. CPT ",
    "Pharmacometrics Syst Pharmacol. 2014;3(1):e94. doi:10.1038/psp.2013.70. ",
    "Final estimates: Table 3 and Supplementary Data ",
    "'Friedrich_LDL_model_code_and_output' (NONMEM 7.2 $PRED and final ",
    "estimates). PK layer: modellib('Friedrich_2014_evacetrapib')."
  )
  vignette <- "Friedrich_2014_evacetrapib"
  units <- list(
    time = "h",
    dosing = "n/a (no dosing events; evacetrapib exposure is supplied as the AUC_EVACETRAPIB covariate in ng*h/mL)",
    concentration = "percent (ldlcfb is the percentage change from baseline in LDL-C)"
  )

  covariateData <- list(
    AUC_EVACETRAPIB = list(
      description = "Individual steady-state evacetrapib AUC over the 24 h once-daily dosing interval",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Time-fixed per subject; the authors fixed each patient's post hoc ",
        "AUC from the final PK model in the analysis dataset (Methods). ",
        "Compute it from modellib('Friedrich_2014_evacetrapib') as ",
        "1000 * DOSE / cl, with DOSE in mg and the individual cl in L/h. ",
        "Set to 0 for placebo and statin-only patients. Besides driving the ",
        "Emax term, AUC_EVACETRAPIB > 0 identifies the evacetrapib-treated ",
        "patients when the combination indicator is derived from ",
        "CONMED_STATIN (see that entry)."
      ),
      source_name = "AUC"
    ),
    APOA1 = list(
      description = "Baseline serum apolipoprotein A-I",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Time-fixed baseline value. Power effect on Emax, ",
        "Emax = -44.1 * (APOA1 / 153)^-0.989, 153 mg/dL being the median ",
        "baseline ApoA1 (Table 3 footnote d; NONMEM $PRED). Lower baseline ",
        "ApoA1 gives a more negative Emax (larger LDL-C reduction)."
      ),
      source_name = "bApoA1 (H20b in the NONMEM $PRED)"
    ),
    LDLC = list(
      description = "Baseline LDL-C",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Time-fixed baseline value. Linear effect on the placebo term, ",
        "-0.119 * (LDLC - 146.56); the NONMEM $PRED centres on 146.56 mg/dL, ",
        "printed as the rounded median 147 in Table 3 footnote e. Higher ",
        "baseline LDL-C gives a lower (more negative) placebo response. ",
        "LDLC here is the covariate; the modelled output is ldlcfb."
      ),
      source_name = "bLDL (G39b in the NONMEM $PRED)"
    ),
    TRIG = list(
      description = "Baseline serum triglycerides",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Time-fixed baseline value. Power effect on the placebo term, ",
        "(TRIG / 120.01)^1.87 - 1; the NONMEM $PRED centres on 120.01 mg/dL, ",
        "printed as the rounded median 120 in Table 3 footnote e. Higher ",
        "baseline triglycerides give a higher placebo response."
      ),
      source_name = "bTRIG (F95b in the NONMEM $PRED)"
    ),
    CONMED_STATIN = list(
      description = "Statin coadministered, 1 = atorvastatin 20 mg, simvastatin 40 mg or rosuvastatin 10 mg once daily",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no statin)",
      notes = paste0(
        "1 for every patient in a statin arm, with or without evacetrapib. ",
        "The paper's two arm indicators are derived from it: I_COMB = ",
        "CONMED_STATIN * (AUC_EVACETRAPIB > 0) for the evacetrapib + statin ",
        "arms and I_ST = CONMED_STATIN - I_COMB for the statin-only arms ",
        "(NONMEM $PRED: I1 for groups 5, 7, 9 and I2 for groups 6, 8, 10). ",
        "The three statins were not distinguished: the model found no ",
        "difference in LDL-C response between them when tested ",
        "individually (Results)."
      ),
      source_name = "I_ST and I_COMB (I1 and I2 in the NONMEM $PRED, from GRP)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on the LDL-C parameters (Methods) and not retained."
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
    HDLC = list(
      description = "Baseline HDL-C",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened (Methods) and not retained in the LDL-C model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 388L,
    n_observations = 1469L,
    n_studies = 1L,
    age_range = "28-83 years (whole randomised study)",
    age_mean = "58.3 years (whole randomised study)",
    weight_range = "44-164 kg (whole randomised study)",
    weight_mean = "83.6 kg (whole randomised study)",
    sex_female_pct = 56.0,
    race_ethnicity = c(White = 92.9),
    disease_state = paste0(
      "Adults with hypercholesterolemia or low HDL-C; mean baseline LDL-C ",
      "144.3 mg/dL (median 146.56), mean triglycerides 139.8 mg/dL (median ",
      "120.01), median ApoA1 153 mg/dL."
    ),
    dose_range = paste0(
      "Placebo; evacetrapib 30, 100 or 500 mg once daily; atorvastatin ",
      "20 mg, simvastatin 40 mg or rosuvastatin 10 mg once daily alone or ",
      "with evacetrapib 100 mg once daily; 12 weeks."
    ),
    co_medication = "Atorvastatin 20 mg, simvastatin 40 mg or rosuvastatin 10 mg once daily in the statin arms.",
    notes = paste0(
      "LDL-C sampled at weeks 2, 4, 8 and 12; 1,469 observations from 388 ",
      "patients. The endpoint is the percent change from each patient's own ",
      "baseline LDL-C; no consistent change between week 2 and week 12 was ",
      "found, so the model has no time course."
    )
  )

  ini({
    # Supplementary Data 'Friedrich_LDL_model_code_and_output' (NONMEM 7.2
    # $PRED, FOCE-I) carries the final estimates; Table 3 prints the same
    # values on the natural scale.
    emax <- -44.1; label("Maximal percent change from baseline in LDL-C from evacetrapib at baseline ApoA1 153 mg/dL (%)") # Table 3 'Emax' = -44.1 (%SEE 8.93); NONMEM THETA(1)
    lauc50 <- 8.35; label("Log of the steady-state AUC giving half the maximal effect (log ng*h/mL)") # NONMEM THETA(2) 8.35, used as EXP(EAUC50) (exp = 4,230); Table 3 'EAUC50' = 4,230 h*ng/ml (%SEE 3.23)
    e_apoa1_emax <- -0.989; label("Power exponent of baseline ApoA1 on Emax (unitless)") # Table 3 'Effect of baseline ApoA1 on Emax' = -0.989 (%SEE 17.1); NONMEM THETA(5); footnote d
    plac <- fixed(0); label("Typical placebo percent change from baseline in LDL-C at median covariates (%)") # Table 3 'PLAC' = 0 (Fixed); NONMEM TVPLAC = 0. Results: a preliminary estimate of 3.88 (%SEE 45.6) did not improve the fit
    e_ldlc_plac <- -0.119; label("Linear effect of baseline LDL-C on the placebo term (%/(mg/dL))") # Table 3 'Effect of baseline LDL on PLAC' = -0.119 (%SEE 19.6); NONMEM THETA(6); footnote e
    e_trig_plac <- 1.87; label("Power exponent of baseline triglycerides in the placebo term (unitless)") # Table 3 'Effect of baseline triglyceride on PLAC' = 1.87 (%SEE 10.5); NONMEM THETA(7); footnote e
    e_conmed_statin_ldl <- -38.7; label("Statin percent change from baseline in LDL-C (%)") # Table 3 'Statin effect' = -38.7 (%SEE 2.92); NONMEM THETA(3)
    e_inter <- -0.997; label("Evacetrapib-statin pharmacodynamic interaction coefficient (unitless)") # Table 3 'Interaction effect' = -0.997 (%SEE 3.42); NONMEM THETA(4). -1 is exact pharmacologic independence

    etaplac ~ 140 # NONMEM OMEGA(1) 140, additive on the placebo term; Table 3 IIV 11.8 = sqrt(140) percentage points

    addSd <- 11.31; label("Additive residual error (percentage points)") # NONMEM SIGMA 128 (variance) -> sqrt(128) = 11.31; Table 3 'Residual error (additive, %)' = 11.3 (%SEE 12.1)
  })
  model({
    emax_i <- emax * (APOA1 / 153)^e_apoa1_emax
    auc50 <- exp(lauc50)
    plac_i <- plac + e_ldlc_plac * (LDLC - 146.56) + (TRIG / 120.01)^e_trig_plac - 1 + etaplac

    drug <- emax_i * AUC_EVACETRAPIB / (auc50 + AUC_EVACETRAPIB)

    # Arm indicators of Eq. 2: statin with evacetrapib (comb) and statin
    # alone (stat_only).
    comb <- CONMED_STATIN * (AUC_EVACETRAPIB > 0)
    stat_only <- CONMED_STATIN - comb

    # Eq. 2. With e_inter = -1 the last term is the percent change of the
    # product of the two fractional responses, (1 + STAT/100) * (1 + LY/100).
    ldlcfb <- plac_i + e_conmed_statin_ldl * stat_only + drug * (1 - comb) +
      e_inter * comb * 100 * (1 - (1 + e_conmed_statin_ldl / 100) * (1 + drug / 100))
    ldlcfb ~ add(addSd)
  })
}
