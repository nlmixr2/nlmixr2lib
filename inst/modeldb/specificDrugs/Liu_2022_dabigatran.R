Liu_2022_dabigatran <- function() {
  description <- paste(
    "Integrated population PK/PD model for oral dabigatran etexilate in",
    "healthy Chinese adults and Chinese adults with non-valvular atrial",
    "fibrillation (Liu 2022, 'Integrated model' columns of Table 2). PK:",
    "two-compartment model with first-order absorption and an absorption lag",
    "for total dabigatran, with age and HDL cholesterol power effects on",
    "apparent clearance; the PK was fitted to the healthy-volunteer data and",
    "held fixed in the integrated fit. PD: direct linear effects of the",
    "dabigatran plasma concentration on activated partial thromboplastin time",
    "(aPTT), prothrombin time (PT) and anti-factor IIa activity, with age and",
    "creatinine clearance on the aPTT baseline and creatinine clearance on the",
    "PT and anti-factor IIa slopes. The sibling Liu_2022_dabigatran_healthy",
    "carries the healthy-volunteer-only PD estimates."
  )
  reference <- paste(
    "Liu Y-O, Xie Q-F, Liu Z-Y, Wang Z, Mu G-Y, Zhang Y-T, Zhao Z-N, Yuan D-D,",
    "Guo L-P, Wang N, Xiang J, Song H-T, Jiang J, Xiang Q, Cui Y-M. Population",
    "pharmacokinetic analysis for dabigatran etexilate in Chinese patients with",
    "non-valvular atrial fibrillation. Front Cardiovasc Med. 2022;9:998751.",
    "doi:10.3389/fcvm.2022.998751"
  )
  vignette <- "Liu_2022_dabigatran"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on apparent clearance, (22 / AGE)^e_age_cl, and on the",
        "aPTT baseline, (22 / AGE)^e_age_rbase_aPTT. The reference 22 years is",
        "the healthy-volunteer median (Liu 2022 Table 1). The printed",
        "Equation 3 form theta_cov^(COV / COVpop) is not used; see the",
        "vignette 'Assumptions and deviations' for the reasoning."
      ),
      source_name = "AGE"
    ),
    HDLC = list(
      description = "Serum high-density lipoprotein cholesterol",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on apparent clearance, (1.42 / HDLC)^e_hdlc_cl, with",
        "reference 1.42 mmol/L, the healthy-volunteer median (Liu 2022",
        "Table 1). Supply in mmol/L (divide mg/dL by 38.67)."
      ),
      source_name = "HDL-C"
    ),
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Cockcroft-Gault creatinine clearance (Liu 2022 Discussion, 'calculated",
        "using the CG formula'), raw mL/min. Power effects (120.18 / CRCL)^e",
        "on the aPTT baseline and on the PT and anti-factor IIa slopes; the",
        "reference 120.18 mL/min is the healthy-volunteer median (Liu 2022",
        "Table 1), which the paper's CrCL 30 and 50 mL/min simulations",
        "reproduce (vignette)."
      ),
      source_name = "CrCL"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "dabigatran etexilate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(
      analyte = "dabigatran (total, after alkaline hydrolysis of the acyl glucuronide)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(analyte = "dabigatran (total)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 285L,
    n_studies = 2L,
    age_range = "18-41 years (healthy volunteers); 39-91 years (patients)",
    age_median = "22 years (healthy volunteers); 73 years (patients)",
    weight_range = "50.3-73.8 kg (healthy volunteers); 44-115 kg (patients)",
    weight_median = "60.3 kg (healthy volunteers); 70 kg (patients)",
    sex_female_pct = 33.7,
    race_ethnicity = "Chinese (all subjects)",
    disease_state = "118 healthy volunteers and 167 adults with non-valvular atrial fibrillation (NVAF)",
    dose_range = paste(
      "Healthy volunteers: single oral 150 mg dabigatran etexilate (Pradaxa),",
      "fasted (n = 59) or after a high-fat meal (n = 59). Patients: 110 mg",
      "twice daily at steady state (Liu 2022 Table 1)."
    ),
    regions = "China (Peking University First Hospital and four further Chinese hospitals)",
    renal_function = "CrCL median 120.18 (86.96-171.02) mL/min in healthy volunteers and 64.48 (23.00-179.55) mL/min in patients (Liu 2022 Table 1)",
    notes = paste(
      "Baseline characteristics from Liu 2022 Table 1. PK data: 1,849",
      "concentrations from the healthy volunteers (rich sampling, 17-18",
      "time points); 23 patients contributed sparse PK used only for external",
      "validation. PD data: 472 aPTT, 468 PT and 447 anti-factor IIa",
      "observations from healthy volunteers (0, 2, 8, 12 h) and 329 aPTT,",
      "330 PT and 331 anti-factor IIa observations from patients (trough",
      "> 10 h post-dose and peak 2 h post-dose at steady state). Female",
      "percentage pooled from Table 1 (33 + 63 of 285). Study registered as",
      "NCT03161496."
    )
  )

  ini({
    # ---- PK (Liu 2022 Table 2, 'Integrated model' column; all held FIX
    # at the healthy-volunteer estimates) ---------------------------------
    lka <- fixed(log(0.261))
    label("Absorption rate constant (1/h)") # Table 2 Ka = 0.261, FIX (healthy model 0.261, RSE 3.3%)
    ltlag <- fixed(log(0.634))
    label("Absorption lag time (h)") # Table 2 ALAG1 = 0.634, FIX (taken from the literature, Results 'PopPK model')
    lfdepot <- fixed(log(1))
    label("Relative bioavailability (unitless)") # Table 2 F1 = 1, FIX
    lcl <- fixed(log(166))
    label("Apparent clearance CL/F at AGE 22 y and HDLC 1.42 mmol/L (L/h)") # Table 2 CL/F = 166, FIX (healthy model 166, RSE 3.7%)
    lvc <- fixed(log(235))
    label("Apparent central volume V2/F (L)") # Table 2 V2/F = 235, FIX (healthy model 235, RSE 11.4%)
    lq <- fixed(log(35.5))
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 Q/F = 35.5, FIX
    lvp <- fixed(log(345))
    label("Apparent peripheral volume V3/F (L)") # Table 2 V3/F = 345, FIX

    e_age_cl <- fixed(0.506)
    label("Power exponent of (22 / AGE) on CL/F (unitless)") # Table 2 'Impact of AGE on CL/F' = 0.506, FIX
    e_hdlc_cl <- fixed(0.694)
    label("Power exponent of (1.42 / HDLC) on CL/F (unitless)") # Table 2 'Impact of HDL on CL/F' = 0.694, FIX

    # ---- aPTT (s) ---------------------------------------------------------
    lrbase_aPTT <- log(32)
    label("Baseline aPTT E0 at AGE 22 y and CRCL 120.18 mL/min (s)") # Table 2 APTT E0 = 32 (RSE 1.7%)
    lslope_aPTT <- log(0.0643)
    label("Linear aPTT slope K1 (s per ng/mL)") # Table 2 K1 = 0.0643 (RSE 17%)
    e_age_rbase_aPTT <- 0.0381
    label("Power exponent of (22 / AGE) on the aPTT baseline (unitless)") # Table 2 'Impact of AGE on E0' = 0.0381 (RSE 27%)
    e_crcl_rbase_aPTT <- 0.252
    label("Power exponent of (120.18 / CRCL) on the aPTT baseline (unitless)") # Table 2 'Impact of CRCL on E0' = 0.252 (RSE 12.2%)

    # ---- PT (s) -----------------------------------------------------------
    lrbase_PT <- log(11.8)
    label("Baseline PT E1 (s)") # Table 2 PT E1 = 11.8 (RSE 0.6%)
    lslope_PT <- log(0.006)
    label("Linear PT slope K2 at CRCL 120.18 mL/min (s per ng/mL)") # Table 2 K2 = 0.006 (RSE 16.2%)
    e_crcl_slope_PT <- 2.02
    label("Power exponent of (120.18 / CRCL) on the PT slope (unitless)") # Table 2 'Impact of CRCL on K2' = 2.02 (RSE 8.8%)

    # ---- Anti-factor IIa activity (ng/mL) -----------------------------------
    lrbase_antiFIIa <- log(20.5)
    label("Baseline anti-factor IIa activity E2 (ng/mL)") # Table 2 Anti-IIa E2 = 20.5 (RSE 9.3%)
    lslope_antiFIIa <- log(0.562)
    label("Linear anti-factor IIa slope K3 at CRCL 120.18 mL/min (ng/mL per ng/mL)") # Table 2 K3 = 0.562 (RSE 11.7%)
    e_crcl_slope_antiFIIa <- 1.35
    label("Power exponent of (120.18 / CRCL) on the anti-factor IIa slope (unitless)") # Table 2 'Impact of CRCL on K3' = 1.35 (RSE 13.8%)

    # ---- IIV (Table 2 'IIV' columns are NONMEM omega^2 variances) ----------
    etalcl ~ fixed(0.142) # Table 2 healthy model CL/F IIV 0.142 (RSE 15.4%); integrated column prints a dash; carried forward with the held-constant PK
    etalvc ~ fixed(1.14) # Table 2 integrated V2/F IIV 1.14, held constant
    etalrbase_aPTT ~ 0.0258 # Table 2 integrated APTT E0 IIV 0.0258 (RSE 17%)
    etalslope_PT ~ 1.25 # Table 2 integrated K2 IIV 1.25 (RSE 20.4%)
    etalslope_antiFIIa ~ 1.1 # Table 2 integrated K3 IIV 1.1 (RSE 17%)

    # ---- Residual error (Table 2 values are NONMEM sigma variances; SD = sqrt) ----
    addSd <- fixed(sqrt(221))
    label("Additive residual SD for dabigatran (ng/mL)") # Table 2 additive residual error 221 (variance), FIX -> SD 14.87
    addSd_aPTT <- sqrt(30.1)
    label("Additive residual SD for aPTT (s)") # Table 2 healthy-model aPTT residual 30.1 (variance) -> SD 5.49; integrated cell misprinted as 0.0258
    addSd_PT <- sqrt(0.89)
    label("Additive residual SD for PT (s)") # Table 2 integrated PT residual 0.89 (variance, RSE 9.9%) -> SD 0.943
    addSd_antiFIIa <- sqrt(1190)
    label("Additive residual SD for anti-factor IIa activity (ng/mL)") # Table 2 integrated anti-IIa residual 1,190 (variance, RSE 10.3%) -> SD 34.5
  })
  model({
    # PK: Liu 2022 Table 2. Covariates enter as (reference / covariate)^exponent
    # with the healthy-volunteer medians as references (see vignette).
    ka <- exp(lka)
    tlag <- exp(ltlag)
    fdepot <- exp(lfdepot)
    cl <- exp(lcl + etalcl) * (22 / AGE)^e_age_cl * (1.42 / HDLC)^e_hdlc_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag
    f(depot) <- fdepot

    # mg/L -> ng/mL
    Cc <- 1000 * central / vc

    # PD: Liu 2022 Equation 5, PD = E0 + K * X1, with X1 the dabigatran
    # plasma concentration (ng/mL).
    rbase_aPTT <- exp(lrbase_aPTT + etalrbase_aPTT) *
      (22 / AGE)^e_age_rbase_aPTT * (120.18 / CRCL)^e_crcl_rbase_aPTT
    slope_aPTT <- exp(lslope_aPTT)
    rbase_PT <- exp(lrbase_PT)
    slope_PT <- exp(lslope_PT + etalslope_PT) * (120.18 / CRCL)^e_crcl_slope_PT
    rbase_antiFIIa <- exp(lrbase_antiFIIa)
    slope_antiFIIa <- exp(lslope_antiFIIa + etalslope_antiFIIa) *
      (120.18 / CRCL)^e_crcl_slope_antiFIIa

    aPTT <- rbase_aPTT + slope_aPTT * Cc
    PT <- rbase_PT + slope_PT * Cc
    antiFIIa <- rbase_antiFIIa + slope_antiFIIa * Cc

    Cc ~ add(addSd)
    aPTT ~ add(addSd_aPTT)
    PT ~ add(addSd_PT)
    antiFIIa ~ add(addSd_antiFIIa)
  })
}
