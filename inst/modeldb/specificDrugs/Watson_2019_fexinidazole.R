Watson_2019_fexinidazole <- function() {
  description <- paste(
    "Population PK model of fexinidazole sulfone (M2), the major active",
    "metabolite of oral fexinidazole, in healthy volunteers, g-HAT patients",
    "and chronic Chagas disease patients (one-compartment disposition with",
    "M2 formation through a depot and two transit compartments, allometric",
    "body weight, population-specific relative bioavailability), with the",
    "cumulative area under the log M2 concentration-time curve (AUlogC) and",
    "the twelve four-parameter sigmoid exposure-response curves that relate",
    "AUlogC to nadir haematology and peak liver transaminases in chronic",
    "Chagas disease patients"
  )
  reference <- paste(
    "Watson JA, Strub-Wourgaft N, Tarral A, Ribeiro I, Tarning J, White NJ.",
    "Pharmacokinetic-pharmacodynamic assessment of the hepatic and bone",
    "marrow toxicities of the new trypanoside fexinidazole. Antimicrob Agents",
    "Chemother. 2019;63(4):e02515-18. doi:10.1128/AAC.02515-18"
  )
  vignette <- "Watson_2019_fexinidazole"
  units <- list(
    time = "h",
    dosing = "umol",
    concentration = "umol/L"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of CL/F (exponent 0.75, fixed) and V/F",
        "(exponent 1, fixed) about the reference weight BWMED = 50 kg of the",
        "supplement NONMEM code."
      ),
      source_name = "BW"
    ),
    DIS_HAT_GAMBIENSE = list(
      description = paste(
        "Patient with human African trypanosomiasis caused by Trypanosoma",
        "brucei gambiense (g-HAT) indicator"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 = not a g-HAT patient (healthy volunteer or Chagas patient)",
      notes = paste(
        "1 = subject from the g-HAT field trials FEX004 (adults, late-stage)",
        "or FEX006 (children 6-14 years > 20 kg and adults, stage 1/2).",
        "Scales the relative bioavailability as F = 1 + e_dis_hat_gambiense_fdepot",
        "(supplement code: F1 = TVF1 * (1 + TVHAT * HAT + TVCHAGAS * CHAGAS)).",
        "DIS_HAT_GAMBIENSE and DIS_CHAGAS_CHRONIC are mutually exclusive; both 0 = healthy",
        "volunteer (phase 1 FEX001-FEX003)."
      ),
      source_name = "HAT"
    ),
    DIS_CHAGAS_CHRONIC = list(
      description = "Patient with chronic indeterminate Chagas disease indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = not a Chagas patient (healthy volunteer or g-HAT patient)",
      notes = paste(
        "1 = subject from the Bolivian dose-finding trial CH-FEX001",
        "(NCT02498782). Scales the relative bioavailability as",
        "F = 1 + e_dis_chagas_chronic_fdepot."
      ),
      source_name = "CHAGAS"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "fexinidazole", units = "umol", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "fexinidazole", units = "umol", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "fexinidazole", units = "umol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "fexinidazole sulfone (M2)", units = "umol", specimen = "plasma", verified = TRUE),
    auc_log_central = list(
      analyte = "fexinidazole sulfone (M2)",
      units = "ln(nmol/L)*day",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 462,
    n_studies = 6,
    age_range = paste(
      "adult healthy volunteers; g-HAT adults (> 15 years) and children",
      "6-14 years weighing > 20 kg; Chagas adults"
    ),
    weight_range = "about 20-100 kg (Figure 2, right panels)",
    weight_median = "50 kg (reference weight BWMED of the NONMEM code)",
    sex_female_pct = 33.8,
    race_ethnicity = paste(
      "phase 1 volunteers of sub-Saharan African origin (France); g-HAT",
      "patients from the Democratic Republic of Congo and the Central African",
      "Republic; Chagas patients from Bolivia"
    ),
    disease_state = paste(
      "healthy volunteers (105), g-HAT sleeping sickness (317) and chronic",
      "indeterminate Chagas disease (40)"
    ),
    dose_range = paste(
      "oral fexinidazole 100-3600 mg single doses and 1200-3600 mg daily",
      "for 14 days (phase 1); 1800 mg x 4 days then 1200 mg x 6 days",
      "(g-HAT adults; 1200 then 600 mg for children 20-35 kg); 1200 or",
      "1800 mg daily for 2, 4 or 8 weeks (Chagas)"
    ),
    regions = "France, Democratic Republic of Congo, Central African Republic, Bolivia",
    notes = paste(
      "Watson 2019 Table 1: FEX001 71, FEX002 12, FEX003 22 (all male",
      "healthy volunteers), FEX004 203 (125 M / 78 F), FEX006 114 (64 M / 50 F),",
      "CH-FEX001 40 (12 M / 28 F); 156 of 462 female. PK fitted to M2",
      "concentrations only (> 4,500 samples; plasma in phase 1, dried blood",
      "spots in the field trials, capillary-blood:plasma ratio 0.97). The",
      "exposure-response curves were fitted to the 40 fexinidazole-treated",
      "Chagas patients only."
    )
  )

  ini({
    # Supplement Table 1 (final NONMEM estimates) and supplement section 2.1
    # (final NONMEM control stream). Amounts are molar: the code takes
    # IPRED = A(2)/V1 and compares log(IPRED) with the natural log of the
    # molar M2 concentration, with no scaling factor (S2 = V1).
    lmtt <- log(5.65); label("Mean transit time of M2 formation, MTT (h)") # Supplement Table 1 TH1 TVMT = 5.65 h (RSE 4.1%)
    lcl <- log(6.37); label("Apparent clearance of M2 at 50 kg, CL/F (L/h)") # Supplement Table 1 TH2 TVCL = 6.37 L/h/F (RSE 5.4%)
    lvc <- log(239); label("Apparent volume of M2 at 50 kg, V/F (L)") # Supplement Table 1 TH3 TVV1 = 239 L/F (RSE 6%)
    lfdepot <- fixed(log(1)); label("Relative bioavailability in healthy volunteers (unitless)") # Supplement Table 1 TH4 TVF1 = 1 FIX
    e_dis_hat_gambiense_fdepot <- 0.956; label("Fractional increase in relative bioavailability in g-HAT patients (unitless)") # Supplement Table 1 TH5 TVHAT = 0.956 (RSE 11.7%)
    e_dis_chagas_chronic_fdepot <- 1.55; label("Fractional increase in relative bioavailability in Chagas patients (unitless)") # Supplement Table 1 TH6 TVCHAGAS = 1.55 (RSE 14.6%)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)") # Supplement code: CL = TVCL * EXP(ETA(2)) * (BW/BWMED)**(0.75)
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V/F (unitless)") # Supplement code: V1 = TVV1 * EXP(ETA(3)) * (BW/BWMED)**1

    etalmtt ~ 0.129 # Supplement Table 1 OM1 IIV MT = 0.129 (RSE 13.6%)
    etalcl ~ 0.0143 # Supplement Table 1 OM2 IIV CL = 0.0143 (RSE 48.3%)
    etalvc ~ 0.0419 # Supplement Table 1 OM3 IIV V1 = 0.0419 (RSE 22%)
    etalfdepot ~ 0.152 # Supplement Table 1 OM4 IIV F1 = 0.152 (RSE 10.6%)

    # Supplement code $ERROR: W = SQRT(SIGMA(1,1)); Y = IPRED + W*ERR(1), with
    # ERR(1) ~ N(0, SIGMA(1,1)). The residual SD on the log scale is therefore
    # W * sqrt(SIGMA(1,1)) = SIGMA(1,1) itself.
    expSd <- 0.391; label("Residual SD of log M2 concentration (log scale)") # Supplement Table 1 SI1 = 0.391 (RSE 4.2%), applied as SD per the $ERROR block

    # Exposure-response curves (Methods Eq. 1 and supplement Stan code):
    # response = e0 + emax / (1 + exp(-exp(lslope) * (AUlogC - auc50))).
    # e0 is the Stan max_effect (low-exposure asymptote), emax the Stan
    # min_effect - max_effect (signed increment to the high-exposure
    # asymptote), auc50 the Stan ed50 and lslope the Stan log_slope. The
    # posterior summaries were not printed except for the 80% credible
    # intervals of ed50 (Table 2); every value below is a least-squares fit of
    # the sigmoid to the posterior-mean curve digitised from Figure 4 or 7.
    # Peak ALT / AST absolute models: log10(peak, U/L)
    e0_log10_peak_alt <- 1.515; label("Peak ALT, low-exposure asymptote (log10 U/L)") # Fig 4 top-left mean curve, digitised (268 points)
    emax_log10_peak_alt <- 0.8042; label("Maximal change in peak ALT (log10 U/L)") # Fig 4 top-left, digitised
    lauc50_log10_peak_alt <- log(261.4); label("AUlogC at half-maximal peak-ALT change (ln(nM)*day)") # Fig 4 top-left, digitised; Table 2 ALT absolute EC50 CI80 164-362
    lslope_log10_peak_alt <- -4.284; label("Log logistic slope of the peak-ALT curve (log of 1/(ln(nM)*day))") # Fig 4 top-left, digitised
    e0_log10_peak_ast <- 1.540; label("Peak AST, low-exposure asymptote (log10 U/L)") # Fig 4 bottom-left mean curve, digitised (278 points)
    emax_log10_peak_ast <- 0.6182; label("Maximal change in peak AST (log10 U/L)") # Fig 4 bottom-left, digitised
    lauc50_log10_peak_ast <- log(307.4); label("AUlogC at half-maximal peak-AST change (ln(nM)*day)") # Fig 4 bottom-left, digitised; Table 2 AST absolute EC50 CI80 244-371
    lslope_log10_peak_ast <- -3.873; label("Log logistic slope of the peak-AST curve (log of 1/(ln(nM)*day))") # Fig 4 bottom-left, digitised
    # Peak ALT / AST relative models: log2(peak / baseline)
    e0_log2_peak_alt_fold <- 0.01828; label("Peak ALT fold change, low-exposure asymptote (log2)") # Fig 4 top-right mean curve, digitised (246 points)
    emax_log2_peak_alt_fold <- 4.068; label("Maximal change in peak ALT fold change (log2)") # Fig 4 top-right, digitised
    lauc50_log2_peak_alt_fold <- log(275.9); label("AUlogC at half-maximal ALT fold change (ln(nM)*day)") # Fig 4 top-right, digitised; Table 2 ALT relative EC50 CI80 119-371
    lslope_log2_peak_alt_fold <- -4.942; label("Log logistic slope of the ALT fold-change curve (log of 1/(ln(nM)*day))") # Fig 4 top-right, digitised
    e0_log2_peak_ast_fold <- 0.7527; label("Peak AST fold change, low-exposure asymptote (log2)") # Fig 4 bottom-right mean curve, digitised (251 points)
    emax_log2_peak_ast_fold <- 2.378; label("Maximal change in peak AST fold change (log2)") # Fig 4 bottom-right, digitised
    lauc50_log2_peak_ast_fold <- log(351.5); label("AUlogC at half-maximal AST fold change (ln(nM)*day)") # Fig 4 bottom-right, digitised; Table 2 AST relative EC50 CI80 256-396
    lslope_log2_peak_ast_fold <- -3.848; label("Log logistic slope of the AST fold-change curve (log of 1/(ln(nM)*day))") # Fig 4 bottom-right, digitised
    # Nadir haematology absolute models
    e0_nadir_neutrophil_count <- 2.293; label("Nadir neutrophil count, low-exposure asymptote (10^3/uL)") # Fig 7 row 1 left mean curve, digitised (237 points)
    emax_nadir_neutrophil_count <- -1.815; label("Maximal change in nadir neutrophil count (10^3/uL)") # Fig 7 row 1 left, digitised
    lauc50_nadir_neutrophil_count <- log(299.2); label("AUlogC at half-maximal nadir-neutrophil change (ln(nM)*day)") # Fig 7 row 1 left, digitised; Table 2 neutrophils absolute EC50 CI80 262-342
    lslope_nadir_neutrophil_count <- -3.676; label("Log logistic slope of the nadir-neutrophil curve (log of 1/(ln(nM)*day))") # Fig 7 row 1 left, digitised
    e0_nadir_platelet_count <- 246.0; label("Nadir platelet count, low-exposure asymptote (10^3/uL)") # Fig 7 row 2 left mean curve, digitised (254 points)
    emax_nadir_platelet_count <- -115.0; label("Maximal change in nadir platelet count (10^3/uL)") # Fig 7 row 2 left, digitised
    lauc50_nadir_platelet_count <- log(310.7); label("AUlogC at half-maximal nadir-platelet change (ln(nM)*day)") # Fig 7 row 2 left, digitised; Table 2 platelets absolute EC50 CI80 146-240 (see vignette)
    lslope_nadir_platelet_count <- -4.611; label("Log logistic slope of the nadir-platelet curve (log of 1/(ln(nM)*day))") # Fig 7 row 2 left, digitised
    e0_nadir_lymphocyte_count <- 1.491; label("Nadir lymphocyte count, low-exposure asymptote (10^3/uL)") # Fig 7 row 3 left mean curve, digitised (261 points)
    emax_nadir_lymphocyte_count <- -0.2232; label("Maximal change in nadir lymphocyte count (10^3/uL)") # Fig 7 row 3 left, digitised
    lauc50_nadir_lymphocyte_count <- log(829.6); label("AUlogC at half-maximal nadir-lymphocyte change (ln(nM)*day)") # Fig 7 row 3 left, digitised; Table 2 lymphocytes absolute EC50 CI80 74-1,304
    lslope_nadir_lymphocyte_count <- -6.096; label("Log logistic slope of the nadir-lymphocyte curve (log of 1/(ln(nM)*day))") # Fig 7 row 3 left, digitised
    e0_nadir_hb <- 14.03; label("Nadir hemoglobin, exposure-independent (g/dL)") # Fig 7 row 4 left mean curve, digitised (235 points), flat to within one pixel
    # Nadir haematology relative models: 100 * (nadir - baseline) / baseline
    e0_nadir_neutrophil_pct_change <- 7.396; label("Nadir neutrophil change from baseline, low-exposure asymptote (%)") # Fig 7 row 1 right mean curve, digitised (237 points)
    emax_nadir_neutrophil_pct_change <- -45.53; label("Maximal change in nadir neutrophil percent change (%)") # Fig 7 row 1 right, digitised
    lauc50_nadir_neutrophil_pct_change <- log(273.7); label("AUlogC at half-maximal neutrophil percent change (ln(nM)*day)") # Fig 7 row 1 right, digitised; Table 2 neutrophils relative EC50 CI80 208-352
    lslope_nadir_neutrophil_pct_change <- -3.521; label("Log logistic slope of the neutrophil percent-change curve (log of 1/(ln(nM)*day))") # Fig 7 row 1 right, digitised
    e0_nadir_platelet_pct_change <- -11.67; label("Nadir platelet change from baseline, low-exposure asymptote (%)") # Fig 7 row 2 right mean curve, digitised (248 points)
    emax_nadir_platelet_pct_change <- -31.03; label("Maximal change in nadir platelet percent change (%)") # Fig 7 row 2 right, digitised
    lauc50_nadir_platelet_pct_change <- log(219.0); label("AUlogC at half-maximal platelet percent change (ln(nM)*day)") # Fig 7 row 2 right, digitised; Table 2 platelets relative EC50 CI80 199-240
    lslope_nadir_platelet_pct_change <- -3.101; label("Log logistic slope of the platelet percent-change curve (log of 1/(ln(nM)*day))") # Fig 7 row 2 right, digitised
    e0_nadir_lymphocyte_pct_change <- -20.50; label("Nadir lymphocyte change from baseline, low-exposure asymptote (%)") # Fig 7 row 3 right mean curve, digitised (227 points)
    emax_nadir_lymphocyte_pct_change <- -16.37; label("Maximal change in nadir lymphocyte percent change (%)") # Fig 7 row 3 right, digitised
    lauc50_nadir_lymphocyte_pct_change <- log(211.8); label("AUlogC at half-maximal lymphocyte percent change (ln(nM)*day)") # Fig 7 row 3 right, digitised; Table 2 lymphocytes relative EC50 CI80 163-396
    lslope_nadir_lymphocyte_pct_change <- -2.993; label("Log logistic slope of the lymphocyte percent-change curve (log of 1/(ln(nM)*day))") # Fig 7 row 3 right, digitised
    e0_nadir_hb_pct_change <- -8.865; label("Nadir hemoglobin change from baseline, low-exposure asymptote (%)") # Fig 7 row 4 right mean curve, digitised (237 points)
    emax_nadir_hb_pct_change <- -5.093; label("Maximal change in nadir hemoglobin percent change (%)") # Fig 7 row 4 right, digitised
    lauc50_nadir_hb_pct_change <- log(303.8); label("AUlogC at half-maximal hemoglobin percent change (ln(nM)*day)") # Fig 7 row 4 right, digitised; Table 2 hemoglobin relative EC50 CI80 137-842
    lslope_nadir_hb_pct_change <- -3.752; label("Log logistic slope of the hemoglobin percent-change curve (log of 1/(ln(nM)*day))") # Fig 7 row 4 right, digitised
  })

  model({
    # Supplement code: BWMED = 50; NN = 2; KTR = (NN+1)/MT
    mtt <- exp(lmtt + etalmtt)
    ktr <- 3 / mtt
    cl <- exp(lcl + etalcl) * (WT / 50)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 50)^e_wt_vc
    kel <- cl / vc

    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(central) <- ktr * transit2 - kel * central

    # F1 = TVF1 * (1 + TVHAT * HAT + TVCHAGAS * CHAGAS) * EXP(ETA(4))
    f(depot) <- exp(lfdepot + etalfdepot) * (1 + e_dis_hat_gambiense_fdepot * DIS_HAT_GAMBIENSE + e_dis_chagas_chronic_fdepot * DIS_CHAGAS_CHRONIC)

    Cc <- central / vc

    # AUlogC: area under the natural-log M2 concentration (nmol/L) curve in
    # days, integrating only while the concentration exceeds 1 nmol/L
    # (ln C > 0); see the vignette for how this definition was checked
    # against Figure 2.
    cc_nmol <- Cc * 1000
    lncc_pos <- log(max(cc_nmol, 1))
    d/dt(auc_log_central) <- lncc_pos / 24

    # Exposure-response curves (Methods Eq. 1); each output is the typical
    # (posterior-mean) summary outcome for the cumulative AUlogC reached so far.
    log10_peak_alt <- e0_log10_peak_alt + emax_log10_peak_alt / (1 + exp(-exp(lslope_log10_peak_alt) * (auc_log_central - exp(lauc50_log10_peak_alt))))
    log10_peak_ast <- e0_log10_peak_ast + emax_log10_peak_ast / (1 + exp(-exp(lslope_log10_peak_ast) * (auc_log_central - exp(lauc50_log10_peak_ast))))
    log2_peak_alt_fold <- e0_log2_peak_alt_fold + emax_log2_peak_alt_fold / (1 + exp(-exp(lslope_log2_peak_alt_fold) * (auc_log_central - exp(lauc50_log2_peak_alt_fold))))
    log2_peak_ast_fold <- e0_log2_peak_ast_fold + emax_log2_peak_ast_fold / (1 + exp(-exp(lslope_log2_peak_ast_fold) * (auc_log_central - exp(lauc50_log2_peak_ast_fold))))
    nadir_neutrophil_count <- e0_nadir_neutrophil_count + emax_nadir_neutrophil_count / (1 + exp(-exp(lslope_nadir_neutrophil_count) * (auc_log_central - exp(lauc50_nadir_neutrophil_count))))
    nadir_platelet_count <- e0_nadir_platelet_count + emax_nadir_platelet_count / (1 + exp(-exp(lslope_nadir_platelet_count) * (auc_log_central - exp(lauc50_nadir_platelet_count))))
    nadir_lymphocyte_count <- e0_nadir_lymphocyte_count + emax_nadir_lymphocyte_count / (1 + exp(-exp(lslope_nadir_lymphocyte_count) * (auc_log_central - exp(lauc50_nadir_lymphocyte_count))))
    nadir_hb <- e0_nadir_hb
    nadir_neutrophil_pct_change <- e0_nadir_neutrophil_pct_change + emax_nadir_neutrophil_pct_change / (1 + exp(-exp(lslope_nadir_neutrophil_pct_change) * (auc_log_central - exp(lauc50_nadir_neutrophil_pct_change))))
    nadir_platelet_pct_change <- e0_nadir_platelet_pct_change + emax_nadir_platelet_pct_change / (1 + exp(-exp(lslope_nadir_platelet_pct_change) * (auc_log_central - exp(lauc50_nadir_platelet_pct_change))))
    nadir_lymphocyte_pct_change <- e0_nadir_lymphocyte_pct_change + emax_nadir_lymphocyte_pct_change / (1 + exp(-exp(lslope_nadir_lymphocyte_pct_change) * (auc_log_central - exp(lauc50_nadir_lymphocyte_pct_change))))
    nadir_hb_pct_change <- e0_nadir_hb_pct_change + emax_nadir_hb_pct_change / (1 + exp(-exp(lslope_nadir_hb_pct_change) * (auc_log_central - exp(lauc50_nadir_hb_pct_change))))

    Cc ~ lnorm(expSd)
  })
}
