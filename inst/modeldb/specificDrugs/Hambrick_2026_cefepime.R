Hambrick_2026_cefepime <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order elimination for",
    "intravenously infused cefepime in critically ill children and young",
    "adults at high risk of acute kidney injury (renal angina index >= 8)",
    "in a pediatric intensive care unit (Hambrick 2026; 49 patients, 296",
    "scavenged plasma concentrations). It re-fits the Morales Junior 2025",
    "PICU cefepime model with every structural, IIV and residual-error",
    "parameter fixed to that model, re-estimating only the power exponent",
    "on BSA-normalized eGFR (0.40) and adding a urine NGAL effect:",
    "clearance is multiplied by exp(-0.42), about 34% lower, when the",
    "urine NGAL concentration at PICU admission is >= 500 ng/mL.",
    "Clearance is 6.38 L/h at 70 kg, eGFR 147.6 mL/min/1.73 m^2 and urine",
    "NGAL < 500 ng/mL; the central volume is 15 L at 70 kg and expands",
    "exponentially with cumulative percentage fluid balance. Allometric",
    "exponents are fixed at 0.75 (CL, Q) and 1 (V1, V2) at 70 kg.",
    "Inter-occasion variability on CL (30.3% CV) and V1 (14.2% CV)",
    "distinguishes the first 48 h of PICU admission (occasion 1) from the",
    "time thereafter (occasion 2). Proportional residual error only.",
    sep = " "
  )
  reference <- paste(
    "Hambrick HR, Morales Junior R, Curry C, Collins M, Johnson L, Mizuno T,",
    "Krallman KA, Goldstein SL, Tang Girdwood S. Urine NGAL adds to serum",
    "creatinine in predicting cefepime clearance in critically ill children",
    "at high risk of acute kidney injury. Int J Antimicrob Agents.",
    "2026;67(5):107741. doi:10.1016/j.ijantimicag.2026.107741. PMID",
    "41690641; PMC13325201. Final parameter estimates and the printed model",
    "equations are from Table 2; the parameters Table 2 marks as fixed are",
    "the final estimates of Morales Junior R et al. Clin Pharmacokinet.",
    "2025;64(4):553-564 (doi:10.1007/s40262-025-01485-5), packaged as",
    "MoralesJunior_2025_cefepime.",
    sep = " "
  )
  vignette <- "Hambrick_2026_cefepime"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # TOTAL plasma cefepime was quantified by HPLC (Hambrick 2026 Methods,
  # 'Blood sampling and cefepime concentration assessment'); the unbound
  # fraction of 0.80 used for the fT>MIC targets is applied outside the PK
  # model (population$protein_binding).
  compartmentData <- list(
    central = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling with exponents FIXED at 0.75 for CL and Q and 1",
        "for V1 and V2, standardised to 70 kg (Hambrick 2026 Table 2",
        "equations, inherited from Morales Junior 2025). Table 1 reports",
        "the weight on study day 1: median 29.4 kg (IQR 20-50). The paper",
        "does not say whether WT was carried time-varying in the fit.",
        sep = " "
      ),
      source_name = "WT"
    ),
    CRCL = list(
      description = "BSA-normalized estimated glomerular filtration rate (bedside Schwartz equation for age < 18 years, race-neutral CKD-EPI for age >= 18 years)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Serum-creatinine eGFR, NOT a Cockcroft-Gault clearance: the",
        "bedside Schwartz equation below 18 years and the race-neutral",
        "SCr-based CKD-EPI equation at 18 years and older (Hambrick 2026",
        "Methods, 'Clinical data collection'). Enters clearance as",
        "(CRCL / 147.6)^0.40 (Table 2); 147.6 is the normalising constant",
        "inherited from Morales Junior 2025. Missing values were filled by",
        "last observation carried backwards, so the column is naturally",
        "time-varying. Median eGFR at PICU admission 62.2 mL/min/1.73 m^2",
        "(IQR 41-122; Table 1). Must be strictly positive (power term).",
        "Patients on any dialysis or extracorporeal therapy at screening",
        "were excluded, and data after such therapy started were dropped.",
        sep = " "
      ),
      source_name = "eGFR"
    ),
    URINE_NGAL = list(
      description = "Urine neutrophil gelatinase-associated lipocalin concentration at PICU admission, not normalised to urine creatinine",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Absolute spot-urine NGAL measured once at PICU admission (reflex",
        "order for renal angina index >= 8), so time-fixed within subject.",
        "Enters clearance through the binary threshold",
        "exp(-0.42 * [URINE_NGAL >= 500]) (Hambrick 2026 Table 2 'e^(beta_NGAL",
        "x NGAL(1 or 0))' and Results: 'Modelling uNGAL as a binary variable",
        "using a threshold of 500 ng/mL yielded the greatest reduction in",
        "the OFV'). The threshold is applied inside model(), so the column",
        "holds the measured concentration. Values below the 50 ng/mL limit",
        "of quantification were set to 50 ng/mL. Table 1: median 154.5",
        "ng/mL (IQR 118-2065); 15 of 50 patients (30%) were >= 500 ng/mL.",
        "Log-transformed and inverse continuous forms and cut-offs of 150",
        "and 1000 ng/mL were also screened (Methods) and not retained. The",
        "early-vs-late (48 h) PICU period was not significant as a fixed",
        "effect on CL and enters only through the IOV (Results). Not",
        "normalised to urine creatinine (unavailable for most patients,",
        "Discussion); use UNGALCR for the creatinine ratio.",
        sep = " "
      ),
      source_name = "uNGAL"
    ),
    CUM_FLUID_BAL_PCT = list(
      description = "Cumulative fluid balance expressed as a percentage of study-day-1 body weight",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Sum of all previous days' and the present day's net fluid balance",
        "(intake minus output, midnight to midnight) divided by the body",
        "weight on study day 1, x 100 (Hambrick 2026 Methods, 'Clinical",
        "data collection'). Naturally time-varying and can be negative.",
        "Enters V1 as exp(0.026 * CUM_FLUID_BAL_PCT), the coefficient fixed",
        "to Morales Junior 2025. The Monte Carlo simulations set it to 0%",
        "for every simulated patient.",
        sep = " "
      ),
      source_name = "Cum%FB"
    ),
    OCC = list(
      description = "Occasion index for inter-occasion variability (1 = first 48 h of PICU admission, 2 = after 48 h)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Hambrick 2026 Methods: 'defining the data within the first 48 h of",
        "PICU admission as one occasion and the data thereafter as a second",
        "occasion, estimating potential inter-occasion variability (IOV) on",
        "both CL and V1'. Time is measured from PICU ADMISSION, not from",
        "the first cefepime dose (cefepime could start up to 48 h before",
        "or after admission). Decomposed inside model() into occ1 / occ2,",
        "which select the per-occasion IOV etas; records with any other",
        "value receive no IOV. Of the 296 concentrations, 175 (59%) were in",
        "occasion 1 (Results).",
        sep = " "
      ),
      source_name = "occasion"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 49L,
    n_enrolled = 50L,
    n_studies = 1L,
    n_sites = 1L,
    n_concentrations = 296L,
    age_range = "3 months to 31 years",
    age_median = "11.5 years (IQR 5-17; Table 1)",
    weight_median = "29.4 kg on study day 1 (IQR 20-50; Table 1)",
    height_median = "127 cm (IQR 101-151; Table 1)",
    sex_female_pct = 50,
    race_ethnicity = "Not reported.",
    disease_state = paste(
      "Critically ill children and young adults admitted to the pediatric",
      "intensive care unit at Cincinnati Children's Hospital Medical Center",
      "with a renal angina index >= 8 (elevated risk of severe acute kidney",
      "injury), cefepime started within 48 h of PICU admission and a urine",
      "NGAL measured as standard of care. 68% had KDIGO SCr-defined AKI of",
      "any stage on PICU admission (66% stage 2-3); 74% were mechanically",
      "ventilated and 66% received vasopressors; median 2 organ failures",
      "(Table 1). Only 10 (20%) had a positive culture (Table S1); most",
      "received cefepime empirically. Patients on dialysis or any",
      "extracorporeal therapy at screening were excluded.",
      sep = " "
    ),
    renal_function = "SCr-eGFR on PICU admission median 62.2 mL/min/1.73 m^2 (IQR 41-122). Urine NGAL on PICU admission median 154.5 ng/mL (IQR 118-2065); 25 (50%) >= 150 ng/mL and 15 (30%) >= 500 ng/mL (Table 1, Results).",
    dose_range = "Clinician-directed. Institutional standard 50 mg/kg per dose (maximum 2000 mg) as a 30-min infusion every 8 h for eGFR > 60, every 12 h for eGFR 30-60 and every 24 h for eGFR < 30 mL/min/1.73 m^2. Monte Carlo simulations covered 50 mg/kg (max 2000 mg) every 8 or 12 h over 30 min or 3 h, and 100 or 150 mg/kg/day (max 4000 / 6000 mg/day) as 24-h continuous infusions.",
    protein_binding = "Not fitted. TOTAL plasma cefepime was measured; the fT>MIC targets assume 20% protein binding, i.e. an unbound fraction of 0.80 (Methods, 'Monte Carlo simulations'). Multiply Cc by 0.80 for free drug.",
    sampling = "Scavenged residual plasma samples from routine care, 2-11 per patient (median 6); 175 of 296 (59%) within 48 h of PICU admission. HPLC assay linear 0.5-200 mg/L, CV < 15%.",
    exclusions = "One patient with paraplegia (supra-normal eGFR of about 200 mL/min/1.73 m^2 from low muscle mass) was excluded from the PK analysis, leaving 49 of 50 enrolled.",
    regions = "United States (single centre, Cincinnati, Ohio)",
    notes = "Parameter estimation by SAEM in Monolix 2024R1; Monte Carlo simulations in Simulx 2024R1. Not externally validated (Discussion)."
  )

  ini({
    # -----------------------------------------------------------------
    # Structural parameters. Hambrick 2026 Table 2 marks CLpop, V1pop,
    # beta_Cum%FB, Qpop, V2pop, IIV CL, IIV V1 and b with an asterisk:
    # 'Fixed to the value from the cefepime model previously developed and
    # validated by Morales Junior et al.' Methods: 'In all model
    # iterations, PK parameters other than tested covariate effects were
    # fixed to the values estimated in the original cefepime popPK model
    # ... The only exception was the beta coefficient corresponding to
    # eGFR's effect on clearance, which was re-estimated'.
    # -----------------------------------------------------------------

    lcl <- fixed(log(6.38))
    label("Clearance at WT = 70 kg, eGFR = 147.6 mL/min/1.73 m^2 and urine NGAL < 500 ng/mL (L/h)")
    # Hambrick 2026 Table 2, 'CLpop (L/h/70kg0.75)' = 6.38*, fixed to
    # Morales Junior 2025.

    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on (WT / 70) shared by CL and Q (unitless)")
    # Hambrick 2026 Table 2 CL and Q equations: '(WT/70)0.75' and
    # 'Qpop x (WT/70)0.75'. A structural constant, not estimated.

    e_crcl_cl <- 0.40
    label("Power exponent on (CRCL / 147.6) for CL (unitless)")
    # Hambrick 2026 Table 2, 'beta eGFR' = 0.40 (RSE 25.5%; bootstrap
    # median 0.37, 95% CI 0.22-0.54). Re-estimated (Morales Junior 2025:
    # 0.66) to absorb collinearity with urine NGAL.

    e_urine_ngal_cl <- -0.42
    label("Exponential coefficient on CL for urine NGAL >= 500 ng/mL at PICU admission (unitless)")
    # Hambrick 2026 Table 2, 'beta NGAL' = -0.42 (RSE 36.9%; bootstrap
    # median -0.47, 95% CI -0.76 to -0.13), in the term
    # 'e^(beta NGAL x NGAL (1 or 0))'. exp(-0.42) = 0.657, the '~35% lower
    # CL' of the Abstract and Discussion.

    lvc <- fixed(log(15))
    label("Central volume of distribution at WT = 70 kg and zero cumulative fluid balance (L)")
    # Hambrick 2026 Table 2, 'V1pop (L/70 kg)' = 15*, fixed to Morales
    # Junior 2025.

    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on (WT / 70) shared by V1 and V2 (unitless)")
    # Hambrick 2026 Table 2 V1 and V2 equations print a bare '(WT/70)',
    # i.e. exponent 1.

    e_cum_fluid_bal_pct_vc <- fixed(0.026)
    label("Exponential coefficient of cumulative percentage of fluid balance on V1 (per %)")
    # Hambrick 2026 Table 2, 'beta Cum%FB' = 0.026*, fixed to Morales
    # Junior 2025. Exponential form per the Morales Junior 2025 Table 2
    # equation (the Hambrick Table 2 V1 row is typeset without the 'e^').

    lq <- fixed(log(3.65))
    label("Intercompartmental clearance at WT = 70 kg (L/h)")
    # Hambrick 2026 Table 2, 'Qpop (L/h/70kg0.75)' = 3.65*, fixed to
    # Morales Junior 2025.

    lvp <- fixed(log(8.91))
    label("Peripheral volume of distribution at WT = 70 kg (L)")
    # Hambrick 2026 Table 2, 'V2pop (L/70 kg)' = 8.91*, fixed to Morales
    # Junior 2025.

    # -----------------------------------------------------------------
    # Variability. Table 2 footnote: IIV and IOV are CVs computed as
    # sqrt(e^(omega^2) - 1) x 100, so omega^2 = log(1 + CV^2):
    #   IIV CL 38.1%: log(1 + 0.381^2) = 0.1355452
    #   IIV V1 14.9%: log(1 + 0.149^2) = 0.0219581
    #   IOV CL 30.3%: log(1 + 0.303^2) = 0.0878360
    #   IOV V1 14.2%: log(1 + 0.142^2) = 0.0199634
    # No random effects on Q or V2 (Methods: sparse sampling precluded
    # them). No correlations reported.
    # -----------------------------------------------------------------

    etalcl ~ fixed(0.1355452)
    # Hambrick 2026 Table 2, 'IIV CL' = 38.1%*, the Morales Junior 2025 value.

    etalvc ~ fixed(0.0219581)
    # Hambrick 2026 Table 2, 'IIV V1' = 14.9%*, the Morales Junior 2025 value.

    # Inter-occasion variability (occasion 1 = first 48 h of PICU
    # admission, occasion 2 = thereafter). Occasion 2 repeats the
    # occasion-1 variance (Monolix estimates one IOV variance per
    # parameter).
    etaiov_cl_1 ~ 0.0878360
    # Hambrick 2026 Table 2, 'IOV CL' = 30.3% (RSE 16.0%; bootstrap median 27%).
    etaiov_cl_2 ~ fixed(0.0878360)
    # Same variance as occasion 1.
    etaiov_vc_1 ~ 0.0199634
    # Hambrick 2026 Table 2, 'IOV V1' = 14.2% (RSE 73.3%; bootstrap median 21%).
    etaiov_vc_2 ~ fixed(0.0199634)
    # Same variance as occasion 1.

    propSd <- fixed(0.319)
    label("Proportional residual error (fraction)")
    # Hambrick 2026 Table 2, 'Error model parameter, proportional only',
    # row 'B' = 31.9%*, fixed to Morales Junior 2025.
  })

  model({
    # Normalising constants as printed in the Table 2 equations.
    wt_ref <- 70 # kg
    crcl_ref <- 147.6 # mL/min/1.73 m^2
    urine_ngal_cut <- 500 # ng/mL, Results: binary threshold retained

    # Binary urine NGAL indicator (1 = admission uNGAL >= 500 ng/mL).
    urine_ngal_high <- (URINE_NGAL >= urine_ngal_cut)

    # Inter-occasion variability selected by the occasion index.
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    iov_cl <- occ1 * etaiov_cl_1 + occ2 * etaiov_cl_2
    iov_vc <- occ1 * etaiov_vc_1 + occ2 * etaiov_vc_2

    # Table 2: CL = CLpop x (WT/70)^0.75 x (eGFR/147.6)^beta
    #               x e^(beta_NGAL x NGAL(1 or 0)) x e^etaCl
    cl <- exp(lcl + etalcl + iov_cl + e_urine_ngal_cl * urine_ngal_high) *
      (WT / wt_ref)^e_wt_cl_q * (CRCL / crcl_ref)^e_crcl_cl

    # Table 2: V1 = V1pop x (WT/70) x e^(beta x Cum%FB) x e^etaV1
    vc <- exp(lvc + etalvc + iov_vc) *
      (WT / wt_ref)^e_wt_vc_vp *
      exp(e_cum_fluid_bal_pct_vc * CUM_FLUID_BAL_PCT)

    # Table 2: Q = Qpop x (WT/70)^0.75 and V2 = V2pop x (WT/70).
    q <- exp(lq) * (WT / wt_ref)^e_wt_cl_q
    vp <- exp(lvp) * (WT / wt_ref)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # TOTAL plasma cefepime (what the assay measured); free = 0.80 x Cc.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
