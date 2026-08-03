Djerada_2014_nefopam <- function() {
  description <- paste(
    "Joint parent + metabolite population PK model for nefopam",
    "(2-compartment; IV zero-order infusion into central) and its",
    "N-demethyl metabolite desmethyl-nefopam (1-compartment; fed by an",
    "apparent metabolic clearance K13 from the nefopam central",
    "compartment) in 48 elderly patients (65-99 years) with normal,",
    "moderate, or severe renal impairment after a single 20 mg 30 min IV",
    "nefopam infusion for postoperative analgesia (hip-fracture repair).",
    "Fit in Monolix v4.1 with SAEM. No covariates were retained in the",
    "final model despite an extensive screen (age, sex, weight, height,",
    "BMI, LBW, FFM, IBW, BSA, iohexol-clearance GFR, MDRD GFR, and",
    "Cockroft-Gault GFR); severe renal impairment (GFR <= 30 mL/min)",
    "correlated post-hoc with reduced nefopam clearance and higher",
    "AUC0-infinity by Spearman test but did not meet the model-building",
    "inclusion criteria."
  )
  reference <- paste(
    "Djerada Z, Fournet-Fayard A, Gozalo C, Lelarge C, Lamiable D,",
    "Millart H, Malinovsky J-M. Population pharmacokinetics of nefopam",
    "in elderly, with or without renal impairment, and its link to",
    "treatment response. Br J Clin Pharmacol. 2014 Jun;77(6):1027-1038.",
    "doi:10.1111/bcp.12291.",
    sep = " "
  )
  vignette <- "Djerada_2014_nefopam"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Table 1: global 83 +/- 8 years (range 65-99).",
        "Screened per Methods (stepwise covariate modelling, forward",
        "LRT P = 0.05, backward LRT P = 0.01, RSE < 20%);",
        "did not meet the inclusion criteria and was not in the final",
        "model."
      )
    ),
    SEXF = list(
      description       = "Female sex indicator (1 = female, 0 = male)",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (male)",
      notes             = paste(
        "Source coded sex as 1 = male, 2 = female (Methods).",
        "Global sex ratio (M:F) = 0.45 => 15 males, 33 females",
        "(Table 1). Canonical SEXF = 1 - (source_sex - 1)",
        "= 2 - source_sex when source is on the 1/2 scale.",
        "Screened; not retained in the final model."
      )
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Table 1: global 63 +/- 13 kg (range 40-100).",
        "Screened; not retained in the final model."
      )
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Table 1: global 161 +/- 7 cm (range 144-183).",
        "Screened; not retained in the final model."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Table 1: global 23.9 +/- 4.7 kg/m^2 (range 15.6-39.1).",
        "Screened; not retained in the final model."
      )
    ),
    LBW = list(
      description = "Lean body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened per Methods (author-cited derivation) alongside FFM,",
        "IBW, and BSA. None retained in the final model."
      )
    ),
    FFM = list(
      description = "Fat-free mass",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened; not retained in the final model."
    ),
    IBW = list(
      description = "Ideal body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened; not retained in the final model."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened; not retained in the final model."
    ),
    GFR_IO = list(
      description = "Glomerular filtration rate estimated from iohexol plasma clearance",
      units       = "mL/min per 1.73 m^2",
      type        = "continuous",
      notes       = paste(
        "Iohexol was administered as a 1 min IV infusion of 5 mL of",
        "Omnipaque 180 starting 29 min after nefopam infusion began;",
        "iohexol plasma clearance is a gold-standard GFR marker in the",
        "elderly (Methods, reference [19]).",
        "Table 1: global GFR_IO 48.36 +/- 24.19 mL/min",
        "(range 10.77-131.50) with three renal-function subgroups",
        "(normal > 60; moderate < 60; severe <= 30 mL/min per 1.73 m^2",
        "using the Cockroft-Gault definition).",
        "Screened as a continuous covariate on CL_nef; did not meet the",
        "inclusion criteria in the population model. A separate",
        "post-hoc Spearman test on individual EBE-derived CL_nef showed",
        "GFR_IO <= 30 mL/min correlated with reduced CL_nef",
        "(Results, Table 2)."
      )
    ),
    GFR_MDRD = list(
      description = "Glomerular filtration rate estimated with the 4-variable Modification of Diet in Renal Disease equation",
      units       = "mL/min per 1.73 m^2",
      type        = "continuous",
      notes       = paste(
        "Screened as an alternative renal-function descriptor.",
        "Not retained in the final model."
      )
    ),
    GFR_CG = list(
      description = "Glomerular filtration rate estimated with the Cockroft-Gault formula",
      units       = "mL/min per 1.73 m^2",
      type        = "continuous",
      notes       = paste(
        "Used to classify patients into normal / moderate / severe",
        "renal-function subgroups a priori (Methods).",
        "Screened as a continuous covariate on PK parameters; not",
        "retained in the final model."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 48,
    n_studies      = 1,
    age_range      = "65-99 years",
    age_median     = "83 years (mean +/- SD = 83 +/- 8)",
    weight_range   = "40-100 kg",
    weight_median  = "63 kg (mean +/- SD = 63 +/- 13)",
    height_range   = "144-183 cm (mean +/- SD = 161 +/- 7)",
    bmi_range      = "15.6-39.1 kg/m^2 (mean +/- SD = 23.9 +/- 4.7)",
    sex_female_pct = 68.75,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Elderly patients scheduled for postoperative repair of a",
      "fractured hip; enrolled between June 2006 and August 2011.",
      "Renal function stratified as normal (Cockroft-Gault GFR",
      "> 60 mL/min per 1.73 m^2; n = 10), moderate impairment",
      "(GFR < 60; n = 28), or severe impairment (GFR <= 30; n = 10).",
      "Iohexol-clearance GFR (GFR_IO) global mean 48.36 +/- 24.19",
      "mL/min (range 10.77-131.50)."
    ),
    dose_range    = "20 mg nefopam hydrochloride diluted in 0.9% saline as a single 30 min IV infusion; iohexol 5 mL of Omnipaque 180 as a 1 min IV infusion starting 29 min after nefopam infusion for GFR estimation.",
    regions       = "France (single-centre: Reims University Hospital).",
    notes         = paste(
      "Table 1 baseline demographics: global (n=48) age 83 +/- 8",
      "years, weight 63 +/- 13 kg, height 161 +/- 7 cm,",
      "BMI 23.9 +/- 4.7 kg/m^2, sex ratio (male:female) = 0.45",
      "(15 males, 33 females = 68.75% female).",
      "Exclusion criteria: nefopam or iodinated-contrast",
      "hypersensitivity, hepatic insufficiency, epilepsy, glaucoma,",
      "prostate adenoma, congestive heart failure, unstable angina,",
      "cardiac-infarction sequelae, uncontrolled arrhythmia,",
      "haemodialysis.",
      "Sampling: 0, 5, 10, 15, 30, 60, 120, 180, 240, 480, 1440 min",
      "after nefopam start; nefopam, desmethyl-nefopam, and iohexol",
      "plasma concentrations assayed per prior published methods.",
      "Fitting: Monolix v4.1.3 with SAEM + MCMC (2 chains), all runs",
      "carried out > 6 times to ensure stability."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Nefopam structural PK parameters (2-compartment; total CL_nef,
    # inter-compartmental Q_nef, central V1_nef, peripheral V2_nef).
    # Values are the "estimated mean +/- SEM" typical values reported
    # in Results (Nefopam pharmacokinetic model, p 1030-1031).
    # BSVs are on the CV-percent scale; omega^2 = log(CV^2 + 1) with
    # CV expressed as a fraction (CV_percent / 100).
    # -----------------------------------------------------------------

    lcl <- log(17.3)
    label("Nefopam total clearance CL_nef (L/h)")
    # Djerada 2014 Results p 1031: CL_nef = 17.3 +/- 1.9 L/h.

    lvc <- log(114)
    label("Nefopam central volume V1_nef (L)")
    # Djerada 2014 Results p 1031: V1_nef = 114 +/- 21 L.

    lq  <- log(80.7)
    label("Nefopam inter-compartmental clearance Q_nef (L/h)")
    # Djerada 2014 Results p 1031: Q_nef = 80.7 +/- 15.0 L/h.

    lvp <- log(208)
    label("Nefopam peripheral volume V2_nef (L)")
    # Djerada 2014 Results p 1031: V2_nef = 208 +/- 31 L.

    # -----------------------------------------------------------------
    # Desmethyl-nefopam structural PK parameters (1-compartment; total
    # CL_dnef, distribution volume V_dnef, and apparent metabolic
    # clearance K13 of nefopam-to-desmethyl-nefopam). Values from
    # Results (Desmethyl-nefopam pharmacokinetic model, p 1031).
    # -----------------------------------------------------------------

    lcl_dnef <- log(17.81)
    label("Desmethyl-nefopam total clearance CL_dnef (L/h)")
    # Djerada 2014 Results p 1031: CL_dnef = 17.81 +/- 0.6 L/h.

    lvc_dnef <- log(12.5)
    label("Desmethyl-nefopam distribution volume V_dnef (L)")
    # Djerada 2014 Results p 1031: V_dnef = 12.5 +/- 0.01 L.

    lcl_form_dnef <- log(0.35)
    label("Apparent nefopam-to-desmethyl-nefopam formation clearance K13 (L/h)")
    # Djerada 2014 Figure 1 legend + Results p 1031: K13 =
    # 0.35 +/- 0.022 L/h. Named K13 in the source; encoded here as the
    # canonical lcl_form_<metab> pattern (see references/parameter-names.md).
    # The apparent value absorbs the metabolite-formation fraction and
    # the metabolite bioavailability into a single positive coefficient
    # (typical "apparent clearance of metabolisation" parameterisation).

    # -----------------------------------------------------------------
    # Between-subject variability. Exponential model: theta_i = theta_TV
    # * exp(eta_i) with diagonal Omega. Correlations were tested but not
    # retained (Methods p 1030). BSVs reported as CV percent; convert to
    # variance via omega^2 = log(CV^2 + 1).
    # -----------------------------------------------------------------

    etalcl           ~ 0.24762
    # Djerada 2014 Results p 1031: BSV(CL_nef) = 53 +/- 9 %CV;
    # omega^2 = log(0.53^2 + 1) = 0.24762.

    etalvc           ~ 0.90184
    # Djerada 2014 Results p 1031: BSV(V1_nef) = 121 +/- 21 %CV;
    # omega^2 = log(1.21^2 + 1) = 0.90184.

    etalq            ~ 0.48489
    # Djerada 2014 Results p 1031: BSV(Q_nef) = 79 +/- 17 %CV;
    # omega^2 = log(0.79^2 + 1) = 0.48489.

    etalvp           ~ 0.34333
    # Djerada 2014 Results p 1031: BSV(V2_nef) = 64 +/- 15 %CV;
    # omega^2 = log(0.64^2 + 1) = 0.34333.

    etalcl_dnef      ~ 1.12250
    # Djerada 2014 Results p 1031: BSV(CL_dnef) = 144 +/- 19 %CV;
    # omega^2 = log(1.44^2 + 1) = 1.12250.

    etalvc_dnef      ~ 0.00010
    # Djerada 2014 Results p 1031: BSV(V_dnef) = 1 +/- 0.9 %CV;
    # omega^2 = log(0.01^2 + 1) = 9.9995e-05, rounded to 1e-04. The
    # near-zero point estimate probably reflects near-unidentifiability
    # of V_dnef separately from K13 in the apparent-clearance
    # parameterisation; the paper reports the value at face and it is
    # preserved here.

    etalcl_form_dnef ~ 0.89200
    # Djerada 2014 Results p 1031: BSV(K13) = 120 +/- 67 %CV;
    # omega^2 = log(1.20^2 + 1) = 0.89200.

    # -----------------------------------------------------------------
    # Residual error: combined additive + proportional on each output.
    # -----------------------------------------------------------------

    addSd <- 1.33
    label("Nefopam additive residual SD (ug/L)")
    # Djerada 2014 Results p 1031: additive a = 1.33 ug/L.

    propSd <- 0.228
    label("Nefopam proportional residual SD (fraction)")
    # Djerada 2014 Results p 1031: proportional coefficient = 0.228.

    addSd_dnef <- 0.01
    label("Desmethyl-nefopam additive residual SD (ug/L)")
    # Djerada 2014 Results p 1031: additive a = 0.01 ug/L.

    propSd_dnef <- 0.68
    label("Desmethyl-nefopam proportional residual SD (fraction)")
    # Djerada 2014 Results p 1031: proportional coefficient = 0.68.
  })

  model({
    # Individual PK parameters. Exponential form with diagonal Omega
    # per Djerada 2014 Methods (Population pharmacokinetics).
    cl           <- exp(lcl + etalcl)
    vc           <- exp(lvc + etalvc)
    q            <- exp(lq  + etalq)
    vp           <- exp(lvp + etalvp)
    cl_dnef      <- exp(lcl_dnef      + etalcl_dnef)
    vc_dnef      <- exp(lvc_dnef      + etalvc_dnef)
    cl_form_dnef <- exp(lcl_form_dnef + etalcl_form_dnef)

    # ODE system. Nefopam: 2-compartment with total elimination CL_nef
    # already accounting for all pathways (Djerada 2014 parameterises
    # CL_nef and K13 as independent apparent parameters; K13 is NOT
    # subtracted from CL_nef in the parent central ODE). Doses arrive on
    # the central compartment as a zero-order IV infusion via the event
    # table's rate/duration columns (see the validation vignette).
    d/dt(central) <-
      q * peripheral1 / vp - (cl + q) * central / vc
    d/dt(peripheral1) <-
      q * central / vc - q * peripheral1 / vp

    # Desmethyl-nefopam: 1-compartment. Apparent formation flux
    # cl_form_dnef * C_nef delivers a mass rate into the metabolite
    # compartment; cl_dnef is total elimination of the metabolite.
    d/dt(central_dnef) <-
      cl_form_dnef * central / vc - cl_dnef * central_dnef / vc_dnef

    # Plasma concentrations. Dose in mg on a volume in L gives
    # concentration in mg/L; the assay reports ug/L (equivalent to
    # ng/mL) so a factor of 1000 converts mg/L -> ug/L.
    Cc      <- 1000 * central      / vc
    Cc_dnef <- 1000 * central_dnef / vc_dnef

    Cc      ~ add(addSd)      + prop(propSd)
    Cc_dnef ~ add(addSd_dnef) + prop(propSd_dnef)
  })
}
