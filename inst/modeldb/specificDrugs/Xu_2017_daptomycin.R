Xu_2017_daptomycin <- function() {
  description <- paste(
    "Two-compartment intravenous population pharmacokinetic model for",
    "daptomycin in adults spanning normal-to-impaired renal function and four",
    "dialysis modalities, updated from the Chaves 2014 renal-impairment model",
    "with patients on continuous veno-venous haemodialysis (CVVHD, n = 9) and",
    "continuous veno-venous haemodiafiltration (CVVHDF, n = 8) (Xu 2017).",
    "Drug input is a zero-order infusion of estimated duration D1 = 0.41 h.",
    "Clearance follows two separate equations: not-on-dialysis patients (Eq 6)",
    "carry a power effect of baseline creatinine clearance referenced to",
    "80 mL/min, while dialysis patients (Eq 1) instead take a modality-specific",
    "clearance (haemodialysis, CAPD, CVVHD or CVVHDF) multiplied by a",
    "dialysis-membrane flux factor. Both equations share a body-temperature",
    "power term referenced to 37 degC, a female multiplier, and five",
    "adjudicated-diagnosis multipliers. CVVHD and CVVHDF additionally take",
    "their own central volume, peripheral volume and inter-compartmental",
    "clearance; body weight scales Q2 and Vp allometrically at estimated",
    "exponents referenced to 70 kg, and a confirmed Gram-positive infection",
    "raises Vp 1.75-fold. Inter-individual variability is a 2x2 block on CL and",
    "Vc plus diagonal terms on Q2 and Vp; the additive residual error switches",
    "between the LC-MS/MS and the non-LC-MS/MS bioanalytical assay."
  )
  reference <- paste(
    "Xu X, Khadzhynov D, Peters H, Chaves RL, Hamed K, Levi M, Corti N (2017).",
    "Population pharmacokinetics of daptomycin in adult patients undergoing",
    "continuous renal replacement therapy. Br J Clin Pharmacol 83(3):498-509.",
    "doi:10.1111/bcp.13131"
  )
  vignette <- "Xu_2017_daptomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight at baseline; allometric power term on Q2 (Eq 3) and Vp (Eq 4), both referenced to 70 kg",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject. Xu 2017 Table 1 gives a pooled median of 75 kg",
        "(range 42.0-152.8); the CVVHD and CVVHDF subgroups had medians of 74 kg",
        "(42-100) and 82 kg (63-120). The exponents are ESTIMATED, not fixed at",
        "an allometric 0.75: 0.797 on Q2 (theta 10) and 0.764 on Vp (theta 11).",
        "Body weight does not act on CL or Vc in this model."
      ),
      source_name = "WT"
    ),
    CRCL = list(
      description = "Baseline creatinine clearance; power effect on clearance in NOT-on-dialysis patients only (Eq 6), referenced to 80 mL/min",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Raw mL/min, NOT BSA-normalized; the source writes it CLC0 ('creatinine",
        "clearance at baseline') and does not name the estimating equation.",
        "Reference 80 mL/min. The term appears ONLY in the not-on-dialysis",
        "clearance equation (Xu 2017 Eq 6); the dialysis equation (Eq 1) has no",
        "renal-function term, consistent with creatinine clearance not being",
        "defined for dialysis-dependent patients. Of the 459 pooled subjects,",
        "385 had CrCl at or above 30 mL/min."
      ),
      source_name = "CLC0"
    ),
    SEXF = list(
      description = "Female sex indicator; multiplicative factor 0.867 on clearance in both the dialysis (Eq 1) and the not-on-dialysis (Eq 6) equation",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Xu 2017 Table 1: 272 male (59%) and 187 female (41%) of 459 pooled",
        "subjects. The source parameterises the effect as theta 8 raised to the",
        "power Sex[Female], so females have 13.3% lower clearance than males."
      ),
      source_name = "Sex[Female]"
    ),
    BODYTEMP = list(
      description = "Body temperature; power term (BODYTEMP / 37)^2.28 on clearance in both clearance equations",
      units = "degC",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline value, time-fixed per subject in this analysis. Reference",
        "37 degC. Xu 2017 Table 1 gives a pooled median of 37.1 degC (range",
        "35.1-40.1); CVVHD median 37.2 (36.5-38.3) and CVVHDF median 36.8",
        "(35.8-37.9). The exponent 2.28 is the least precisely estimated fixed",
        "effect in the model (60% RSE, Table S2), so the term moves clearance by",
        "only about +/-8% across the observed temperature range."
      ),
      source_name = "TEMP"
    ),
    RRT_HEMODIAL_STATUS = list(
      description = "Intermittent-haemodialysis indicator; selects the haemodialysis clearance 0.219 L/h in place of the not-on-dialysis clearance",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not receiving intermittent haemodialysis)",
      notes = paste(
        "REPLACEMENT rule, not an additive arm: when the indicator is 1 the",
        "whole not-on-dialysis clearance term (including its creatinine-clearance",
        "power effect) is switched out and the modality-specific clearance is",
        "used instead (Xu 2017 Eq 1 vs Eq 6). Time-fixed per subject; 40 of the",
        "459 pooled subjects were on haemodialysis Q48h or thrice weekly.",
        "Mutually exclusive with PERIT_DIAL, RRT_CVVHD_STATUS and",
        "RRT_CVVHDF_STATUS; all four zero selects the not-on-dialysis equation."
      ),
      source_name = "HD"
    ),
    PERIT_DIAL = list(
      description = "Continuous ambulatory peritoneal dialysis indicator; selects the CAPD clearance 0.237 L/h in place of the not-on-dialysis clearance",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not on peritoneal dialysis)",
      notes = paste(
        "REPLACEMENT rule, as for RRT_HEMODIAL_STATUS. Time-fixed per subject;",
        "14 of the 459 pooled subjects were on CAPD. Unlike the Takama 2007",
        "founding example, the effect here is on clearance rather than on central",
        "volume: CAPD patients retain the same Vc, Q2 and Vp as the",
        "not-on-dialysis and haemodialysis strata."
      ),
      source_name = "CAPD"
    ),
    RRT_CVVHD_STATUS = list(
      description = "Continuous veno-venous haemodialysis indicator; selects the CVVHD-specific CL, Vc, Q2 and Vp",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not on CVVHD)",
      notes = paste(
        "REPLACEMENT rule on all four disposition parameters. Time-fixed per",
        "subject; 9 of the 459 pooled subjects. CVVHD clearance (0.936 L/h) is",
        "1.25-fold the not-on-dialysis typical value and roughly 4-fold the",
        "haemodialysis and CAPD values, consistent with continuous rather than",
        "intermittent solute removal. Xu 2017 Methods describe the procedures as",
        "a high-flux 1.4-1.8 m^2 filter with blood flow 100-200 mL/min and a",
        "target dialysis flow of 30-40 mL/kg/h; the dose recommendation is stated",
        "to apply only to comparable procedures."
      ),
      source_name = "CVVHD"
    ),
    RRT_CVVHDF_STATUS = list(
      description = "Continuous veno-venous haemodiafiltration indicator; selects the CVVHDF-specific CL, Vc, Q2 and Vp",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not on CVVHDF)",
      notes = paste(
        "REPLACEMENT rule on all four disposition parameters. Time-fixed per",
        "subject; 8 of the 459 pooled subjects, all of them on a high-flux",
        "membrane (Xu 2017 Table 1). The CVVHDF clearance (0.528 L/h) is 29%",
        "LOWER than the not-on-dialysis typical value, which the authors",
        "attribute to prefilter fluid substitution and to lower flow rates than",
        "in the CVVHD cohort. Because every CVVHDF subject carried the high-flux",
        "membrane factor 1.36, the CVVHDF clearance actually realised in the",
        "source cohort is about 0.72 L/h, not 0.528."
      ),
      source_name = "CVVHDF"
    ),
    FILT_FLUX_LOW = list(
      description = "Low-flux dialysis-membrane indicator; multiplicative factor 0.96 on clearance in dialysis patients (Eq 1)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (high-flux membrane or membrane not recorded)",
      notes = paste(
        "Paired with FILT_FLUX_HIGH; the two are mutually exclusive and both are",
        "0 when the membrane was not recorded, which is the reference level and",
        "covers 424 of 459 subjects (Xu 2017 Table 1: low flux 7, high flux 28,",
        "not available 424). The term appears only in the dialysis clearance",
        "equation (Eq 1), not in the not-on-dialysis equation (Eq 6)."
      ),
      source_name = "DIAM[Low flux]"
    ),
    FILT_FLUX_HIGH = list(
      description = "High-flux dialysis-membrane indicator; multiplicative factor 1.36 on clearance in dialysis patients (Eq 1)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (low-flux membrane or membrane not recorded)",
      notes = paste(
        "See FILT_FLUX_LOW. All 8 CVVHDF subjects and 2 of the 9 CVVHD subjects",
        "carried a high-flux membrane, so this 36% clearance increase is",
        "load-bearing when reproducing the CVVHDF exposures of Xu 2017 Table 3."
      ),
      source_name = "DIAM[High flux]"
    ),
    DIS_IEAC_LEFT = list(
      description = "Left-sided infective endocarditis indicator (IEAC category 1); multiplicative factor 1.16 on clearance",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (diagnosis not adjudicated / not available, the all-zero reference of the five IEAC indicators)",
      notes = paste(
        "One of five mutually exclusive indicators decomposing the independent",
        "external adjudication committee (IEAC) final diagnosis. Xu 2017 Table 1",
        "gives 9 (2%) of 459 subjects; 337 (73%) had no adjudicated diagnosis and",
        "form the reference level. The factor enters both clearance equations."
      ),
      source_name = "IEAC[IEAC 1]"
    ),
    DIS_IEAC_RIGHT_COMP = list(
      description = "Complicated right-sided infective endocarditis indicator (IEAC category 2); multiplicative factor 1.3 on clearance",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (diagnosis not adjudicated / not available)",
      notes = "Xu 2017 Table 1: 13 (3%) of 459 subjects. See DIS_IEAC_LEFT for the shared reference level.",
      source_name = "IEAC[IEAC 2]"
    ),
    DIS_IEAC_RIGHT_UNCOMP = list(
      description = "Uncomplicated right-sided infective endocarditis indicator (IEAC category 3); multiplicative factor 1.3 on clearance",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (diagnosis not adjudicated / not available)",
      notes = paste(
        "Xu 2017 Table 1: 5 (1%) of 459 subjects. The point estimate coincides",
        "with the complicated-RIE factor at 1.3 but the two were estimated",
        "separately (theta 16 and theta 17, Table S2)."
      ),
      source_name = "IEAC[IEAC 3]"
    ),
    DIS_BACTEREMIA_COMP = list(
      description = "Complicated Staphylococcus aureus bacteraemia indicator (IEAC category 4); multiplicative factor 1.10 on clearance",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (diagnosis not adjudicated / not available)",
      notes = "Xu 2017 Table 1: 58 (13%) of 459 subjects, the most frequent adjudicated diagnosis.",
      source_name = "IEAC[IEAC 4]"
    ),
    DIS_BACTEREMIA_UNCOMP = list(
      description = "Uncomplicated Staphylococcus aureus bacteraemia indicator (IEAC category 5); multiplicative factor 1.13 on clearance",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (diagnosis not adjudicated / not available)",
      notes = "Xu 2017 Table 1: 37 (8%) of 459 subjects.",
      source_name = "IEAC[IEAC 5]"
    ),
    DIS_INFECT_ACTIVE = list(
      description = "Confirmed Gram-positive infection indicator; multiplicative factor 1.75 on peripheral volume (Eq 4)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no Gram-positive infection; healthy volunteers and subjects only suspected of infection)",
      notes = paste(
        "Time-fixed per subject in this analysis. Xu 2017 writes it INFN,",
        "'presence of Gram-positive infection'; Table 1 gives 273 (59%) yes and",
        "186 (41%) no, the 'no' stratum being footnoted as healthy volunteers and",
        "some subjects only SUSPECTED of infection. All 17 CRRT subjects were",
        "infected. This is the largest single covariate effect in the model: a",
        "75% increase in peripheral volume, consistent with the interstitial",
        "fluid expansion of acute infection."
      ),
      source_name = "INFN"
    ),
    ASSAY_LCMSMS = list(
      description = "LC-MS/MS bioanalytical-method indicator; selects the smaller of the two additive residual-error magnitudes",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-LC-MS/MS chromatographic assay, i.e. the HPLC methods of the older studies)",
      notes = paste(
        "Per-sample indicator. Xu 2017 supplementary Table S2 footnotes the two",
        "residual variances as 'the study that used LC/MS/MS assay' (ASSY1,",
        "sigma^2 = 2.07 mg^2/L^2, SD 1.44 mg/L) and 'studies that did not use the",
        "LC/MS/MS assay' (ASSY0, sigma^2 = 5.31 mg^2/L^2, SD 2.30 mg/L). The",
        "non-LC-MS/MS comparator is HPLC, which the Table S2 abbreviation list",
        "spells out; it is NOT an immunoassay, so the sibling IMMUNOASSAY and",
        "RIA_ASSAY canonicals do not apply. In a pure LC-MS/MS prospective",
        "dataset set ASSAY_LCMSMS = 1 for every row."
      ),
      source_name = "ASSY"
    )
  )

  compartmentData <- list(
    central = list(analyte = "daptomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "daptomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 459,
    n_studies = 3,
    age_range = "adults; the 17 CRRT subjects were 56-90 years (supplementary Table S1)",
    weight_range = "42.0-152.8 kg (median 75)",
    sex_female_pct = 41,
    disease_state = paste(
      "Gram-positive infection (complicated and uncomplicated Staphylococcus",
      "aureus bacteraemia, right- and left-sided infective endocarditis) plus",
      "healthy volunteers and subjects suspected of infection; renal function",
      "spanning CrCl at or above 30 mL/min (n = 385), intermittent",
      "haemodialysis (n = 40), CAPD (n = 14) and CRRT (CVVHD n = 9, CVVHDF n = 8)"
    ),
    dose_range = "intravenous daptomycin; the simulations cover 4, 6, 8, 10 and 12 mg/kg Q24h and Q48h",
    regions = "pooled multinational clinical-trial database (Chaves 2014) plus two European CRRT studies (Zurich, Berlin)",
    assay = paste(
      "Plasma daptomycin by LC-MS/MS in one study and by non-LC-MS/MS (HPLC)",
      "methods in the others; the assay indicator was retained in the final",
      "residual-error model with a separate additive magnitude for each."
    ),
    notes = paste(
      "Xu 2017 Methods and Results: PK data from patients on CVVHD (n = 9) and",
      "CVVHDF (n = 8), reported by Corti 2013 and Khadzhynov 2011, were pooled",
      "with the database of the Chaves 2014 renal-impairment daptomycin model,",
      "itself an update of the Dvorchik 2004 model. Three subjects whose dialysis",
      "status was unknown were EXCLUDED from the present analysis, leaving 456",
      "analysed of the 459 tabulated. Estimation used NONMEM 7.2.0 with FOCEI.",
      "Eta shrinkage was 10.0% (CL), 7.52% (Vc), 28.5% (Q2) and 40.9% (Vp)."
    )
  )

  ini({
    # ---- Clearance, not-on-dialysis stratum (Xu 2017 Eq 6) ----
    lcl <- log(0.751); label("Clearance, not on dialysis (L/h)")  # Table S2, CL NOT-ON-DIALYSIS = theta 6 = 0.751 (3% RSE); Table 2 prints 0.75
    e_crcl_cl <- 0.540; label("Exponent on (CRCL / 80 mL/min) for clearance, not-on-dialysis only (unitless)")  # Table S2, *(CLC0/80)^theta 7 = 0.540 (7% RSE)

    # ---- Clearance, dialysis strata (Xu 2017 Eq 1; replacement rule) ----
    lcl_hemodialysis <- log(0.219); label("Clearance, intermittent haemodialysis (L/h)")  # Table S2, CL HD = theta 20 = 0.219 (6% RSE); Table 2 prints 0.22
    lcl_capd <- log(0.237); label("Clearance, continuous ambulatory peritoneal dialysis (L/h)")  # Table S2, CL CAPD = theta 21 = 0.237 (8% RSE); Table 2 prints 0.24
    lcl_cvvhd <- log(0.936); label("Clearance, CVVHD (L/h)")  # Table S2, CL CVVHD = theta 22 = 0.936 (6% RSE); Table 2 prints 0.94
    lcl_cvvhdf <- log(0.528); label("Clearance, CVVHDF (L/h)")  # Table S2, CL CVVHDF = theta 23 = 0.528 (14% RSE); Table 2 prints 0.53
    e_filt_flux_low_cl <- 0.96; label("Low-flux dialysis membrane factor on clearance (unitless)")  # Table S2, *theta 13 DIAM [Low Flux] = 0.96 (14% RSE)
    e_filt_flux_high_cl <- 1.36; label("High-flux dialysis membrane factor on clearance (unitless)")  # Table S2, *theta 14 DIAM [High Flux] = 1.36 (8% RSE)

    # ---- Clearance covariates shared by Eq 1 and Eq 6 ----
    e_bodytemp_cl <- 2.28; label("Exponent on (BODYTEMP / 37 degC) for clearance (unitless)")  # Table S2, *(TEMP/37)^theta 9 = 2.28 (60% RSE)
    e_sexf_cl <- 0.867; label("Female factor on clearance (unitless)")  # Table S2, *theta 8 SEX [Female] = 0.867 (3% RSE)
    e_ieac_left_cl <- 1.16; label("Left-sided infective endocarditis factor on clearance (unitless)")  # Table S2, *theta 15 IEAC [LIE] = 1.16 (12% RSE)
    e_ieac_right_comp_cl <- 1.3; label("Complicated right-sided infective endocarditis factor on clearance (unitless)")  # Table S2, *theta 16 IEAC [Complicated RIE] = 1.3 (9% RSE)
    e_ieac_right_uncomp_cl <- 1.3; label("Uncomplicated right-sided infective endocarditis factor on clearance (unitless)")  # Table S2, *theta 17 IEAC [Uncomplicated RIE] = 1.3 (10% RSE)
    e_bacteremia_comp_cl <- 1.10; label("Complicated bacteraemia factor on clearance (unitless)")  # Table S2, *theta 18 IEAC [Complicated Bacteraemia] = 1.10 (4% RSE)
    e_bacteremia_uncomp_cl <- 1.13; label("Uncomplicated bacteraemia factor on clearance (unitless)")  # Table S2, *theta 19 IEAC [Uncomplicated Bacteraemia] = 1.13 (5% RSE)

    # ---- Central volume (Xu 2017 Eq 2) ----
    lvc <- log(4.86); label("Central volume, all strata except CVVHD and CVVHDF (L)")  # Table S2, V1 other = theta 2 = 4.86 (4% RSE)
    lvc_cvvhd <- log(5.736); label("Central volume, CVVHD (L)")  # Table S2, V1 CVVHD = theta 24 = 5.736 (8% RSE); Table 2 prints 5.74
    lvc_cvvhdf <- log(6.526); label("Central volume, CVVHDF (L)")  # Table S2, V1 CVVHDF = theta 25 = 6.526 (6% RSE); Table 2 prints 6.53

    # ---- Inter-compartmental clearance (Xu 2017 Eq 3) ----
    lq <- log(3.69); label("Inter-compartmental clearance, all strata except CVVHD and CVVHDF (L/h)")  # Table S2, Q other = theta 3 = 3.69 (6% RSE)
    lq_cvvhd <- log(7.111); label("Inter-compartmental clearance, CVVHD (L/h)")  # Table S2, Q CVVHD = theta 26 = 7.111 (15% RSE); Table 2 prints 7.11
    lq_cvvhdf <- log(2.879); label("Inter-compartmental clearance, CVVHDF (L/h)")  # Table S2, Q CVVHDF = theta 27 = 2.879 (35% RSE); Table 2 prints 2.88
    e_wt_q <- 0.797; label("Exponent on (WT / 70 kg) for inter-compartmental clearance (unitless)")  # Table S2, *(WT/70)^theta 10 = 0.797 (28% RSE)

    # ---- Peripheral volume (Xu 2017 Eq 4) ----
    lvp <- log(3.2); label("Peripheral volume, all strata except CVVHD and CVVHDF (L)")  # Table S2, V2 other = theta 4 = 3.2 (3% RSE)
    lvp_cvvhd <- log(4.891); label("Peripheral volume, CVVHD (L)")  # Table S2, V2 CVVHD = theta 28 = 4.891 (7% RSE); Table 2 prints 4.89
    lvp_cvvhdf <- log(3.846); label("Peripheral volume, CVVHDF (L)")  # Table S2, V2 CVVHDF = theta 29 = 3.846 (16% RSE); Table 2 prints 3.85
    e_wt_vp <- 0.764; label("Exponent on (WT / 70 kg) for peripheral volume (unitless)")  # Table S2, *(WT/70)^theta 11 = 0.764 (14% RSE)
    e_infect_active_vp <- 1.75; label("Confirmed Gram-positive infection factor on peripheral volume (unitless)")  # Table S2, *theta 12 INFN [INFN1] = 1.75 (6% RSE)

    # ---- Zero-order input duration (Xu 2017 Eq 5) ----
    ld1 <- log(0.41); label("Zero-order intravenous input duration D1 (h)")  # Table S2, D1 (h) = theta 5 = 0.41 (0.2% RSE)

    # ---- Inter-individual variability (Xu 2017 Table S2, 'Inter-individual variance') ----
    # 2x2 block on CL and Vc. The covariance is not printed; it is recovered
    # from the printed correlation r = 0.52 as 0.52 * sqrt(0.296 * 0.654).
    etalcl + etalvc ~ c(
      0.296,
      0.2288, 0.654
    )
    etalq ~ 0.669  # Table S2, omega^2 Q = 0.669
    etalvp ~ 0.267  # Table S2, omega^2 V2 = 0.267

    # ---- Residual error (Xu 2017 Table S2, 'Residual variance'; additive only) ----
    addSd_lcmsms <- 1.4387; label("Additive residual SD, LC-MS/MS assay (mg/L)")  # Table S2, sigma^2 add[ASSY1] = 2.07 -> SD = sqrt(2.07)
    addSd_other <- 2.3043; label("Additive residual SD, non-LC-MS/MS assay (mg/L)")  # Table S2, sigma^2 add[ASSY0] = 5.31 -> SD = sqrt(5.31)
  })

  model({
    # Dialysis-modality selectors. The four indicators are mutually exclusive;
    # all four zero selects the not-on-dialysis clearance equation (Eq 6).
    onCvvhd <- RRT_CVVHD_STATUS
    onCvvhdf <- RRT_CVVHDF_STATUS
    onDialysis <- RRT_HEMODIAL_STATUS + PERIT_DIAL + onCvvhd + onCvvhdf
    otherStratum <- 1 - onCvvhd - onCvvhdf

    # Covariate multipliers shared by Eq 1 and Eq 6.
    covTempCl <- (BODYTEMP / 37)^e_bodytemp_cl
    covSexCl <- e_sexf_cl^SEXF
    covIeacCl <- e_ieac_left_cl^DIS_IEAC_LEFT *
      e_ieac_right_comp_cl^DIS_IEAC_RIGHT_COMP *
      e_ieac_right_uncomp_cl^DIS_IEAC_RIGHT_UNCOMP *
      e_bacteremia_comp_cl^DIS_BACTEREMIA_COMP *
      e_bacteremia_uncomp_cl^DIS_BACTEREMIA_UNCOMP

    # Eq 1: dialysis clearance, modality-specific and scaled by membrane flux.
    clDialysis <- (RRT_HEMODIAL_STATUS * exp(lcl_hemodialysis) +
      PERIT_DIAL * exp(lcl_capd) +
      onCvvhd * exp(lcl_cvvhd) +
      onCvvhdf * exp(lcl_cvvhdf)) *
      e_filt_flux_low_cl^FILT_FLUX_LOW *
      e_filt_flux_high_cl^FILT_FLUX_HIGH

    # Eq 6: not-on-dialysis clearance, driven by baseline creatinine clearance.
    clNotOnDialysis <- exp(lcl) * (CRCL / 80)^e_crcl_cl

    cl <- (onDialysis * clDialysis + (1 - onDialysis) * clNotOnDialysis) *
      covTempCl * covSexCl * covIeacCl * exp(etalcl)

    # Eq 2: central volume; CVVHD and CVVHDF have their own values.
    vc <- (otherStratum * exp(lvc) +
      onCvvhd * exp(lvc_cvvhd) +
      onCvvhdf * exp(lvc_cvvhdf)) * exp(etalvc)

    # Eq 3: inter-compartmental clearance, allometric on body weight.
    q <- (otherStratum * exp(lq) +
      onCvvhd * exp(lq_cvvhd) +
      onCvvhdf * exp(lq_cvvhdf)) *
      (WT / 70)^e_wt_q * exp(etalq)

    # Eq 4: peripheral volume, allometric on body weight and raised by infection.
    vp <- (otherStratum * exp(lvp) +
      onCvvhd * exp(lvp_cvvhd) +
      onCvvhdf * exp(lvp_cvvhdf)) *
      (WT / 70)^e_wt_vp *
      e_infect_active_vp^DIS_INFECT_ACTIVE * exp(etalvp)

    # Eq 5: estimated zero-order input duration. Dose records must carry
    # rate = -2 for rxode2 to honour dur(); a plain bolus ignores it.
    d1 <- exp(ld1)
    dur(central) <- d1

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Eq 7: plasma concentration in the central compartment.
    Cc <- central / vc

    # Additive residual error switched per sample by the bioanalytical method.
    addSd <- ASSAY_LCMSMS * addSd_lcmsms + (1 - ASSAY_LCMSMS) * addSd_other
    Cc ~ add(addSd)
  })
}
