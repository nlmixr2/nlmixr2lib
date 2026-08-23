Gracia_2025_cr51edta_tc99mdtpa <- function() {
  description <- paste(
    "One-compartment IV-bolus population PK model fitted simultaneously to",
    "51Cr-EDTA and 99mTc-DTPA plasma concentrations in 59 oncopediatric",
    "children (40 51Cr-EDTA / 19 99mTc-DTPA) receiving cisplatin and/or",
    "ifosfamide (Gracia 2025). Clearance is the measured glomerular",
    "filtration rate and follows the CYSPED power equation on serum",
    "creatinine, plasma cystatin C and body weight; volume is proportional",
    "to body weight. A binary per-record tracer indicator",
    "(TRACER_TC99M_DTPA; the paper's MES covariate) enters both CL and V",
    "multiplicatively as (1 - theta * TRACER_TC99M_DTPA) and carries the",
    "paper's estimand, the inter-tracer bias of -0.9% (95% CI +/- 11.4%)",
    "on clearance. Inter-occasion variability is carried on CL, and the",
    "proportional residual error is tracer-specific."
  )
  reference <- paste(
    "Gracia M, Ankaoua V, Alonso M, Pasquet M, Chatelut E (2025).",
    "A population pharmacokinetic approach to compare 51Cr-EDTA and",
    "99mTc-DTPA clearances in measuring renal glomerular filtration rate",
    "in oncopediatrics. Pediatric Nephrology 40(10):3163-3168.",
    "doi:10.1007/s00467-025-06828-9.",
    sep = " "
  )
  vignette <- "Gracia_2025_cr51edta_tc99mdtpa"

  # Gracia 2025 reports CL in mL/min and V in mL, with sampling at 90, 110,
  # 130 and 150 min after an IV bolus expressed in MBq. The paper never
  # states the plasma-concentration unit; MBq/mL follows from dose in MBq
  # divided by a volume in mL. Any consistent radioactivity unit gives the
  # same CL, because CL = dose / AUC is unit-invariant in the tracer amount.
  units <- list(time = "min", dosing = "MBq", concentration = "MBq/mL")

  compartmentData <- list(
    central = list(
      analyte  = "51Cr-EDTA or 99mTc-DTPA (the tracer selected by TRACER_TC99M_DTPA)",
      units    = "MBq",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source symbol Scr. Enters CL as the power term (Scr/42.25)^theta3",
        "with theta3 = -0.568 (Gracia 2025 Table 2). The 42.25 umol/L",
        "centering constant is inherited unchanged from the CYSPED equation",
        "quoted in Methods 'Pharmacokinetic analysis'. Measured by an",
        "enzymatic creatininase method on a Cobas 8000 analyzer (Methods",
        "'Patients and data'). Cohort mean 42 umol/L (range 11-102) in the",
        "51Cr-EDTA group and 42 umol/L (range 8-76) in the 99mTc-DTPA group",
        "(Table 1). Time-varying: collected at baseline and at the end of",
        "treatment, i.e. once per GFR measurement occasion."
      ),
      source_name        = "Scr"
    ),
    CYSC = list(
      description        = "Plasma cystatin C",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source symbol PcysC. Enters CL as the power term",
        "(PcysC/0.7738)^theta4 with theta4 = -0.295 (Gracia 2025 Table 2).",
        "Measured by automated particle-enhanced nephelometric immunoassay",
        "(PENIA) (Methods 'Patients and data'). Cohort mean 0.76 mg/L",
        "(range 0.26-1.36) in the 51Cr-EDTA group and 0.82 mg/L (range",
        "0.56-1.14) in the 99mTc-DTPA group (Table 1); the pooled 59-patient",
        "mean is (40*0.76 + 19*0.82)/59 = 0.779 mg/L, consistent with the",
        "0.7738 mg/L centering constant printed in Table 2. Note that the",
        "Discussion restates the same equation with the older 40-patient",
        "CYSPED constant 0.76 -- Table 2 is the final-model parameter table",
        "and is used here; see the vignette Errata.",
        "Time-varying: one value per GFR measurement occasion."
      ),
      source_name        = "PcysC"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source symbol BW. Enters CL as the power term (BW/35.70)^theta5",
        "with theta5 = +1.06, and enters V linearly as V = theta7 * BW with",
        "theta7 = 256 mL/kg (Gracia 2025 Table 2). The volume term is not",
        "centered -- theta7 is a per-kilogram coefficient, so V is",
        "proportional to body weight with no reference weight. Cohort mean",
        "34.3 kg (range 9.07-72.0) in the 51Cr-EDTA group and 42.9 kg",
        "(range 10.0-91.8) in the 99mTc-DTPA group (Table 1). The",
        "Discussion restates the CL equation with the older 40-patient",
        "CYSPED constant 34.3 kg; Table 2's 35.70 kg is used here (see the",
        "vignette Errata). Time-varying across occasions.",
        "Body weight also sets the 99mTc-DTPA bolus dose (6, 9 or 12 MBq",
        "for < 10, 10-20 and > 20 kg)."
      ),
      source_name        = "BW"
    ),
    TRACER_TC99M_DTPA = list(
      description        = "Exogenous GFR tracer identity: 1 = 99mTc-DTPA, 0 = 51Cr-EDTA",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (51Cr-EDTA, the historical reference tracer and the group supplying 40 of the 59 children)",
      notes              = paste(
        "Source column MES ('GFR measurement'), defined in Gracia 2025",
        "Methods 'Pharmacokinetic analysis' and in the Table 2 footnote as",
        "'MES = 0 or = 1, respectively, for 51Cr-EDTA or 99mTc-DTPA data'.",
        "Enters both CL and V in the paper's printed multiplicative form",
        "(1 - theta * MES), with theta6 = -0.009 on CL and theta8 = +0.017",
        "on V, and additionally switches the proportional residual error",
        "between 4.99% CV (51Cr-EDTA) and 7.18% CV (99mTc-DTPA) (Table 2).",
        "Because theta6 is negative, a 99mTc-DTPA record has CL multiplied",
        "by 1.009 -- the paper's headline result, a non-significant -0.9%",
        "(95% CI +/- 11.4%) bias with 99mTc-DTPA clearance slightly",
        "overestimated. Per-record covariate. In the source data the tracer",
        "is perfectly confounded with the source study (Cysped NCT02822404",
        "supplied all 51Cr-EDTA data; CysPedVal RnIPH 2020-101 all",
        "99mTc-DTPA data; 'No patient successively received both",
        "radioisotopic tracers'), but the estimand is the tracer effect, so",
        "the tracer-identity canonical is used rather than a STUDY_*",
        "indicator."
      ),
      source_name        = "MES"
    ),
    OCC = list(
      description        = "GFR measurement occasion",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "One occasion per GFR measurement. Gracia 2025 Methods 'Patients",
        "and data': 'Between one and five (median three) GFR measurements",
        "were taken per patient depending on the number of cycles of",
        "chemotherapy administered', so OCC takes values 1..5. Decomposed",
        "inside model() into binary indicators oc1..oc5 that multiplex the",
        "five inter-occasion-variability etas on log-CL (Methods:",
        "'Interoccasion variability of CL was included')."
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 59,
    n_studies      = 2,
    age_range      = "1.4-17.8 years (51Cr-EDTA group); 2.0-16.3 years (99mTc-DTPA group)",
    age_mean       = "10.6 years (51Cr-EDTA group); 10.1 years (99mTc-DTPA group)",
    weight_range   = "9.07-72.0 kg (51Cr-EDTA group); 10.0-91.8 kg (99mTc-DTPA group)",
    weight_mean    = "34.3 kg (51Cr-EDTA group); 42.9 kg (99mTc-DTPA group)",
    sex_female_pct = 55.9,
    disease_state  = paste(
      "Pediatric malignancy (osteosarcoma 33, neuroblastoma 9, germinal",
      "tumor 7, rhabdomyosarcoma 3, other 7) treated with cisplatin and/or",
      "ifosfamide and requiring GFR monitoring for chemotherapy",
      "nephrotoxicity"
    ),
    dose_range     = paste(
      "2.072 MBq 51Cr-EDTA IV bolus; 6, 9 or 12 MBq 99mTc-DTPA IV bolus for",
      "body weight < 10 kg, 10-20 kg and > 20 kg respectively"
    ),
    regions        = "France (Toulouse University Hospital)",
    renal_function = "Not restricted at inclusion; GFR monitored for nephrotoxicity of cisplatin and/or ifosfamide",
    notes          = paste(
      "Baseline demographics from Gracia 2025 Table 1 (characteristics",
      "before the first GFR measurement). 40 children from the Cysped",
      "clinical trial (NCT02822404, 51Cr-EDTA, February 2012 to September",
      "2015) and 19 from the CysPedVal analysis (RnIPH 2020-101,",
      "99mTc-DTPA, December 2020 to June 2022). Inclusion criteria: age",
      "0-18 years and body weight over 4 kg. Sex 26 male / 33 female",
      "overall (18/22 and 8/11 by group). Each tracer was given as an IV",
      "bolus twice per patient (before chemotherapy and at end of",
      "treatment) with four blood samples drawn from the contralateral arm",
      "at 90, 110, 130 and 150 min. Between one and five (median three)",
      "GFR measurements per patient. No patient received both tracers."
    )
  )

  ini({
    # Structural typical values at the reference covariate set
    # (Scr = 42.25 umol/L, PcysC = 0.7738 mg/L, BW = 35.70 kg, 51Cr-EDTA).
    lcl <- log(80.8)
    label("Typical clearance, i.e. measured glomerular filtration rate, at the reference covariate set (mL/min)")  # Gracia 2025 Table 2: theta1 = 80.8 +/- 3.1
    lvc <- log(256)
    label("Typical volume of distribution per kilogram of body weight (mL/kg)")  # Gracia 2025 Table 2: theta7 = 256 +/- 20

    # Covariate effects on CL -- power form, Table 2 CL equation, inherited
    # from the CYSPED equation and re-estimated on the pooled 59 children.
    e_creat_cl <- -0.568
    label("Power exponent on serum creatinine (CREAT/42.25) for CL (unitless)")  # Gracia 2025 Table 2: theta3 = -0.568 +/- 0.156
    e_cysc_cl <- -0.295
    label("Power exponent on plasma cystatin C (CYSC/0.7738) for CL (unitless)")  # Gracia 2025 Table 2: theta4 = -0.295 +/- 0.168
    e_wt_cl <- 1.06
    label("Power exponent on body weight (WT/35.70) for CL (unitless)")  # Gracia 2025 Table 2: theta5 = +1.06 +/- 0.145

    # Tracer effect. Entered exactly in the paper's printed multiplicative
    # form (1 - theta * TRACER_TC99M_DTPA), so the SIGN convention is the
    # paper's: a NEGATIVE theta raises the 99mTc-DTPA value.
    e_tracer_tc99m_dtpa_cl <- -0.009
    label("99mTc-DTPA-vs-51Cr-EDTA effect on CL, applied as (1 - theta * TRACER_TC99M_DTPA) (unitless)")  # Gracia 2025 Table 2: theta6 = -0.009 +/- 0.114
    e_tracer_tc99m_dtpa_vc <- 0.017
    label("99mTc-DTPA-vs-51Cr-EDTA effect on V, applied as (1 - theta * TRACER_TC99M_DTPA) (unitless)")  # Gracia 2025 Table 2: theta8 = +0.017 +/- 0.143

    # Inter-individual variability. Gracia 2025 Methods reports CV% for the
    # random effects; the log-scale variance that reproduces the published
    # CV is omega^2 = log(CV^2 + 1).
    etalcl ~ 0.0120274  # Gracia 2025 Table 2: interindividual variability on CL = 11% CV -> log(0.11^2 + 1)
    etalvc ~ 0.0427539  # Gracia 2025 Table 2: interindividual variability on V = 20.9% CV -> log(0.209^2 + 1)

    # Inter-occasion variability on CL (Gracia 2025 Methods: "Interoccasion
    # variability of CL was included"; Table 2: 14.5% CV). NONMEM shares one
    # variance across occasions via $OMEGA BLOCK(1) SAME; nlmixr2 has no SAME
    # shortcut, so each occasion gets its own eta with the variance fixed to
    # the occasion-1 value after the first (the convention used in
    # Jonsson_2011_ethambutol.R, Chen_2023_nemonoxacin.R and
    # Waterhouse_2024_vedolizumab.R).
    etaiov_cl_1 ~ 0.0208070       # Gracia 2025 Table 2: IOV on CL = 14.5% CV -> log(0.145^2 + 1)
    etaiov_cl_2 ~ fix(0.0208070)  # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_3 ~ fix(0.0208070)  # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_4 ~ fix(0.0208070)  # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_5 ~ fix(0.0208070)  # SAME-equivalent: equal to the occasion-1 IOV variance

    # Tracer-specific proportional residual error. Gracia 2025 Methods:
    # "Considering that some observed plasma concentrations corresponded to
    # 51Cr-EDTA data and others to 99mTc-DTPA data, specific residual
    # variability was allocated to 51Cr-EDTA and 99mTc-DTPA data."
    propSd_cr51edta <- 0.0499
    label("Proportional residual SD for 51Cr-EDTA plasma concentrations (fraction)")  # Gracia 2025 Table 2: residual error, 51Cr-EDTA = 4.99% CV
    propSd_tc99mdtpa <- 0.0718
    label("Proportional residual SD for 99mTc-DTPA plasma concentrations (fraction)")  # Gracia 2025 Table 2: residual error, 99mTc-DTPA = 7.18% CV
  })

  model({
    # Decompose the occasion column into binary indicators to multiplex the
    # five inter-occasion-variability etas onto log-CL.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5

    # Individual parameters. CL is the measured glomerular filtration rate;
    # the covariate model is the CYSPED power equation re-estimated on the
    # pooled 59-patient dataset (Gracia 2025 Table 2):
    #   TVCL = theta1 * (Scr/42.25)^theta3 * (PcysC/0.7738)^theta4 *
    #          (BW/35.70)^theta5 * (1 - theta6 * MES)
    #   TVV  = theta7 * BW * (1 - theta8 * MES)
    cl <- exp(lcl + etalcl + iov_cl) *
      (CREAT / 42.25)^e_creat_cl *
      (CYSC / 0.7738)^e_cysc_cl *
      (WT / 35.70)^e_wt_cl *
      (1 - e_tracer_tc99m_dtpa_cl * TRACER_TC99M_DTPA)
    vc <- exp(lvc + etalvc) * WT *
      (1 - e_tracer_tc99m_dtpa_vc * TRACER_TC99M_DTPA)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Observation. The proportional residual SD is switched by the tracer
    # indicator; TRACER_TC99M_DTPA is binary, so exactly one term survives.
    propSd_eff <- propSd_cr51edta * (1 - TRACER_TC99M_DTPA) +
      propSd_tc99mdtpa * TRACER_TC99M_DTPA

    Cc <- central / vc
    Cc ~ prop(propSd_eff)
  })
}
