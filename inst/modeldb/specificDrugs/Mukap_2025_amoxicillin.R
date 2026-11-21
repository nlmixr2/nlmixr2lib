Mukap_2025_amoxicillin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous amoxicillin in",
    "neonates with confirmed or suspected sepsis in Papua New Guinea",
    "(Mukap 2025), fitted to dried-blood-spot (DBS) concentrations with an",
    "estimated plasma-to-DBS ratio for venous plasma samples. Allometric",
    "scaling on body weight (reference 2.8 kg; exponents fixed at 0.75 on",
    "CL and Q and 1 on Vc and Vp), exponential maturation of CL with",
    "postmenstrual age (centred on 39 weeks) and a 57% step increase in CL",
    "for postnatal age above 2 days. Interindividual variability on CL, Vc",
    "and the plasma-to-DBS ratio, interoccasion variability on CL and Vc",
    "with each dose an occasion, and additive residual error on the log",
    "scale for DBS and plasma."
  )
  reference <- paste(
    "Mukap M, Sprod C, Yoo O, Tefuarani N, Vince J, Laman M, Page-Sharp M,",
    "Moore BR, Batty KT, Davis TME, Salman S, Manning L. A simplified",
    "amoxicillin regimen with dose frequency based on post-natal age in",
    "neonates with confirmed or suspected infection. Antimicrob Agents",
    "Chemother. 2025;69(4):e01491-24. doi:10.1128/aac.01491-24.",
    sep = " "
  )
  vignette <- "Mukap_2025_amoxicillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "amoxicillin",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "amoxicillin",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling (WT / 2.8)^0.75 on CL and Q and (WT / 2.8)^1 on",
        "Vc and Vp, exponents fixed (Supplement section 2: 'Allometric",
        "scaling was pre-established ... with set exponents of 3/4 for",
        "clearance and 1 for volume'). 2.8 kg is the cohort median (Table",
        "1; Table 2 parameter units 'per 2.8 kg'). Cohort range 0.9-5.5 kg."
      ),
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Carried in WEEKS for this model (the register default is months):",
        "the maturation term FMAT = exp(0.0816 * (PMA - 39)) of Table 2 is",
        "centred on 39 weeks. Table 1 reports PMA in days (median 273,",
        "range 196-328 days, i.e. 28-47 weeks); missing GA was replaced by",
        "the cohort median GA when forming PMA (Table 1 footnote). The",
        "Discussion states the exponential form is only considered plausible",
        "over the studied PMA range."
      ),
      source_name = "PMA"
    ),
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Supplied in months, as the register requires. The paper's CL",
        "covariate is dichotomous on postnatal age in DAYS: CLCOV = 1 for",
        "PNA <= 2 days and 1 + 0.57 for PNA > 2 days (Table 2). model()",
        "compares PNA against 2 / 30.4375 months. Results prose says 'babies",
        "> 3 days old'; the Table 2 equation (> 2 days, i.e. from day 3 for",
        "integer ages) is encoded. Cohort PNA median 0 days (IQR 0-2, range",
        "0-48)."
      ),
      source_name = "PNA"
    ),
    OCC = list(
      description = "Integer-valued occasion indicator for interoccasion variability on CL and Vc (one occasion per dose)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Supplement section 2: 'Each dose was considered as an individual",
        "event, thus incorporating inter-occasion variability into the",
        "clearance and volume of distribution.' OCC = k labels the dosing",
        "interval that starts at the k-th dose (records before the first",
        "dose take OCC = 1). The paper reports one IOV magnitude per",
        "parameter and no occasion count; 21 occasions are encoded here,",
        "covering 7 days of 8-hourly dosing. Occasions 2-21 are fixed equal",
        "to occasion 1. Records with OCC outside 1..21 carry no IOV."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Not retained: Results 'No other significant covariate relationships were discovered.'"
    ),
    HCT = list(
      description = "Haematocrit",
      units = "fraction",
      type = "continuous",
      notes = paste(
        "Measured (median 0.49, range 0.06-0.74; Results) and discussed as",
        "the likely explanation of the DBS/plasma ratio, but not retained as",
        "a covariate in the final model (Results: no other significant",
        "covariate relationships)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 177L,
    n_studies = 1L,
    age_range = "Postnatal age 0-48 days (median 0, IQR 0-2 days)",
    weight_range = "0.9-5.5 kg (median 2.8, IQR 2.1-3.3 kg)",
    ga_range = "28-42 weeks (median 38, IQR 37-40)",
    pma_range = "196-328 days, i.e. 28-47 weeks (median 273 days)",
    sex_female_pct = 43,
    race_ethnicity = "Papua New Guinean neonates (not further reported).",
    disease_state = paste(
      "Neonates (<= 2 months) with bacteriologically confirmed (9.0%) or",
      "clinically suspected sepsis; 61.0% clinically confirmed and 30.0%",
      "with at least one sepsis indicator. Twelve (6.8%) died. More than",
      "three-quarters were term babies with early-onset sepsis."
    ),
    dose_range = paste(
      "30 mg/kg intravenous amoxicillin as a 3-min slow injection (median",
      "80 mg, range 15-96 mg; 9.4-59 mg/kg), twice daily for the first",
      "7 days of life and thrice daily thereafter (PMGH practice)."
    ),
    regions = "Papua New Guinea (Port Moresby General Hospital neonatal unit).",
    notes = paste(
      "Prospective observational PK study. DBS heel-prick samples at",
      "baseline, 0.5, 1, 4, 12, 24 and 48 h after the first dose (1,145",
      "retained) and one venous plasma sample at about 24 h (134 retained).",
      "Sex: 57% male among those with known sex (12 unknown); the 43%",
      "female figure assumes the same denominator. Fitted in NONMEM 7.5.1",
      "with FOCE-I (Supplement section 2)."
    )
  )

  ini({
    # Structural parameters: Mukap 2025 Table 2 final model, for a 2.8 kg
    # neonate at PMA 39 weeks and PNA <= 2 days.
    lcl <- log(0.856)
    label("Clearance for a 2.8 kg neonate at PMA 39 weeks, PNA <= 2 days (L/h)")
    # Table 2 row 'theta1' 0.856 (RSE 4.8%), CL per 2.8 kg
    lvc <- log(4.11)
    label("Central volume of distribution for a 2.8 kg neonate (L)")
    # Table 2 row 'VC, l 2.8 kg-1' 4.11 (RSE 4.2%)
    lq <- log(0.406)
    label("Intercompartmental clearance for a 2.8 kg neonate (L/h)")
    # Table 2 row 'Q, lh-1 2.8 kg-1' 0.406 (RSE 12.3%)
    lvp <- log(1.35)
    label("Peripheral volume of distribution for a 2.8 kg neonate (L)")
    # Table 2 row 'VP, l 2.8 kg-1' 1.35 (RSE 8.8%)
    lratio_plasma_central <- log(0.783)
    label("Plasma-to-DBS concentration ratio (unitless)")
    # Table 2 row 'Ratio in plasma to in DBS' 0.783 (RSE 7%)

    # Allometric exponents, fixed (Supplement section 2).
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL and Q (unitless)")
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on Vc and Vp (unitless)")

    # Covariate effects on CL (Table 2).
    e_page_cl <- 0.0816
    label("Exponential maturation coefficient of CL per week of PMA (1/week)")
    # Table 2 row 'theta exp, PMA' 0.0816 (RSE 17%)
    e_pna_cl <- 0.57
    label("Fractional increase in CL for postnatal age above 2 days (unitless)")
    # Table 2 row 'theta Dichotomous, PNA' 0.57 (RSE 27.7%)

    # Interindividual variability. Table 2 footnote: variability is reported
    # as 100 x sqrt(variance), so variance = (value / 100)^2.
    etalcl ~ 0.245025
    # Table 2 IIV row 'CL' 49.5; 0.495^2
    etalvc ~ 0.178084
    # Table 2 IIV row 'VC' 42.2; 0.422^2
    etalratio_plasma_central ~ 0.227529
    # Table 2 IIV row 'Ratio' 47.7; 0.477^2

    # Interoccasion variability on CL, one occasion per dose. One shared
    # variance; occasions 2-21 repeat it.
    etaiov_cl_1 ~ 0.101124
    # Table 2 IOV row 'CL' 31.8; 0.318^2
    etaiov_cl_2 ~ fixed(0.101124)
    etaiov_cl_3 ~ fixed(0.101124)
    etaiov_cl_4 ~ fixed(0.101124)
    etaiov_cl_5 ~ fixed(0.101124)
    etaiov_cl_6 ~ fixed(0.101124)
    etaiov_cl_7 ~ fixed(0.101124)
    etaiov_cl_8 ~ fixed(0.101124)
    etaiov_cl_9 ~ fixed(0.101124)
    etaiov_cl_10 ~ fixed(0.101124)
    etaiov_cl_11 ~ fixed(0.101124)
    etaiov_cl_12 ~ fixed(0.101124)
    etaiov_cl_13 ~ fixed(0.101124)
    etaiov_cl_14 ~ fixed(0.101124)
    etaiov_cl_15 ~ fixed(0.101124)
    etaiov_cl_16 ~ fixed(0.101124)
    etaiov_cl_17 ~ fixed(0.101124)
    etaiov_cl_18 ~ fixed(0.101124)
    etaiov_cl_19 ~ fixed(0.101124)
    etaiov_cl_20 ~ fixed(0.101124)
    etaiov_cl_21 ~ fixed(0.101124)

    # Interoccasion variability on Vc, one occasion per dose.
    etaiov_vc_1 ~ 0.119716
    # Table 2 IOV row 'VC' 34.6; 0.346^2
    etaiov_vc_2 ~ fixed(0.119716)
    etaiov_vc_3 ~ fixed(0.119716)
    etaiov_vc_4 ~ fixed(0.119716)
    etaiov_vc_5 ~ fixed(0.119716)
    etaiov_vc_6 ~ fixed(0.119716)
    etaiov_vc_7 ~ fixed(0.119716)
    etaiov_vc_8 ~ fixed(0.119716)
    etaiov_vc_9 ~ fixed(0.119716)
    etaiov_vc_10 ~ fixed(0.119716)
    etaiov_vc_11 ~ fixed(0.119716)
    etaiov_vc_12 ~ fixed(0.119716)
    etaiov_vc_13 ~ fixed(0.119716)
    etaiov_vc_14 ~ fixed(0.119716)
    etaiov_vc_15 ~ fixed(0.119716)
    etaiov_vc_16 ~ fixed(0.119716)
    etaiov_vc_17 ~ fixed(0.119716)
    etaiov_vc_18 ~ fixed(0.119716)
    etaiov_vc_19 ~ fixed(0.119716)
    etaiov_vc_20 ~ fixed(0.119716)
    etaiov_vc_21 ~ fixed(0.119716)

    # Residual error: additive on the log scale for both matrices
    # (Supplement section 2; Table 2 footnote defines AddLnDBS and
    # AddLnPlasma).
    expSd <- 0.182
    label("Additive residual error on log DBS concentration (log scale SD)")
    # Table 2 row 'AddLnDBS' 18.2
    expSd_plasma <- 0.584
    label("Additive residual error on log plasma concentration (log scale SD)")
    # Table 2 row 'PropPlasma' 58.4, defined as AddLnPlasma in the footnote
  })

  model({
    # Occasion indicators; a record with OCC outside 1..21 carries no IOV.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    oc8 <- (OCC == 8)
    oc9 <- (OCC == 9)
    oc10 <- (OCC == 10)
    oc11 <- (OCC == 11)
    oc12 <- (OCC == 12)
    oc13 <- (OCC == 13)
    oc14 <- (OCC == 14)
    oc15 <- (OCC == 15)
    oc16 <- (OCC == 16)
    oc17 <- (OCC == 17)
    oc18 <- (OCC == 18)
    oc19 <- (OCC == 19)
    oc20 <- (OCC == 20)
    oc21 <- (OCC == 21)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6 +
      oc7 * etaiov_cl_7 + oc8 * etaiov_cl_8 + oc9 * etaiov_cl_9 +
      oc10 * etaiov_cl_10 + oc11 * etaiov_cl_11 + oc12 * etaiov_cl_12 +
      oc13 * etaiov_cl_13 + oc14 * etaiov_cl_14 + oc15 * etaiov_cl_15 +
      oc16 * etaiov_cl_16 + oc17 * etaiov_cl_17 + oc18 * etaiov_cl_18 +
      oc19 * etaiov_cl_19 + oc20 * etaiov_cl_20 + oc21 * etaiov_cl_21
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2 + oc3 * etaiov_vc_3 +
      oc4 * etaiov_vc_4 + oc5 * etaiov_vc_5 + oc6 * etaiov_vc_6 +
      oc7 * etaiov_vc_7 + oc8 * etaiov_vc_8 + oc9 * etaiov_vc_9 +
      oc10 * etaiov_vc_10 + oc11 * etaiov_vc_11 + oc12 * etaiov_vc_12 +
      oc13 * etaiov_vc_13 + oc14 * etaiov_vc_14 + oc15 * etaiov_vc_15 +
      oc16 * etaiov_vc_16 + oc17 * etaiov_vc_17 + oc18 * etaiov_vc_18 +
      oc19 * etaiov_vc_19 + oc20 * etaiov_vc_20 + oc21 * etaiov_vc_21

    # Table 2: FMAT = exp(theta_exp,PMA * (PMA - 39)), PMA in weeks.
    fmat <- exp(e_page_cl * (PAGE - 39))
    # Table 2: CLCOV = 1 for PNA <= 2 days, 1 + theta for PNA > 2 days.
    # PNA is in months; 2 days = 2 / 30.4375 months.
    clcov <- 1 + e_pna_cl * (PNA > 2 / 30.4375)

    # Table 2: CL = theta1 * FMAT * CLCOV per 2.8 kg, allometric exponents
    # fixed at 0.75 (CL, Q) and 1 (Vc, Vp).
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 2.8)^e_wt_cl * fmat * clcov
    vc <- exp(lvc + etalvc + iov_vc) * (WT / 2.8)^e_wt_vc
    q <- exp(lq) * (WT / 2.8)^e_wt_cl
    vp <- exp(lvp) * (WT / 2.8)^e_wt_vc
    ratio_plasma_central <- exp(lratio_plasma_central + etalratio_plasma_central)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3-min intravenous injections enter central (Methods).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # DBS concentration: mg / L = mg/L. The structural model is fitted to
    # the DBS (capillary whole-blood) concentrations.
    Cc <- central / vc
    # Plasma concentration: Table 2 'Ratio in plasma to in DBS' and Methods
    # (Simulations): 'the simulated DBS concentrations were then adjusted to
    # free plasma concentrations by multiplication'. Encoded as printed; the
    # paper's paired samples (DBS/plasma 0.48) point the other way -- see
    # the vignette.
    Cc_plasma <- Cc * ratio_plasma_central

    Cc ~ lnorm(expSd)
    Cc_plasma ~ lnorm(expSd_plasma)
  })
}
