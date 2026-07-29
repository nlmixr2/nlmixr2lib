Nath_2007_melphalan <- function() {
  description <- "Two-compartment IV population PK model for melphalan in paediatric blood or marrow transplant recipients (Nath 2007). Structural CL is a linear additive function of body weight, prior-carboplatin therapy, and 99mTc-DTPA-tracer-measured GFR; central volume Vc is a linear additive function of body weight; intercompartmental rate constants k12 and k21 are estimated directly (not as Q/Vc and Q/Vp)."
  reference   <- paste(
    "Nath CE, Shaw PJ, Montgomery K, Earl JW.",
    "Population pharmacokinetics of melphalan in paediatric blood or marrow",
    "transplant recipients.",
    "Br J Clin Pharmacol. 2007;64(2):151-164.",
    "doi:10.1111/j.1365-2125.2007.02862.x.",
    sep = " "
  )
  vignette <- "Nath_2007_melphalan"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (paediatric; baseline at the time of melphalan dose).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear (not allometric power) scaling on both CL and Vc per Table 4:",
        "CL = theta5 * WT + theta6 * CPT + theta7 * GFR;",
        "Vc = theta2 + theta8 * WT.",
        "Development cohort median 18.8 kg (IQR 13.5-28.1); overall study weight",
        "range 7.7-104 kg (Discussion). No reference weight is used because the",
        "WT scaling is linear (not divisive)."
      ),
      source_name        = "WT"
    ),
    PRIOR_CARBOPLATIN = list(
      description        = paste(
        "Binary indicator for prior carboplatin chemotherapy administered as",
        "part of the BMT conditioning block on each of the 5 days preceding",
        "the melphalan dose. Carboplatin dose was determined by the Calvert",
        "formula targeting AUC = 4 mg/mL/min using each child's GFR and BSA."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no prior carboplatin)",
      notes              = paste(
        "Additive linear effect on CL with coefficient theta6 = -3.17 L/h",
        "(Table 4). Prior-carboplatin recipients have melphalan CL reduced by",
        "3.17 L/h compared with carboplatin-naive recipients of the same body",
        "weight and GFR; the mechanism is residual renal-tubular injury from",
        "the 5-day carboplatin pre-treatment block. Development cohort 17/39",
        "received prior carboplatin."
      ),
      source_name        = "CPT"
    ),
    CRCL = list(
      description        = paste(
        "Glomerular filtration rate measured directly by 99mTc-DTPA",
        "(99mTc-diethylenetriaminepentacetic acid) plasma-clearance tracer,",
        "BSA-normalised to mL/min/1.73 m^2. Gold-standard direct measurement",
        "of true GFR; chosen by the source paper because creatinine-based",
        "estimates are unreliable in paediatric oncology / BMT cohorts where",
        "muscle mass and serum-creatinine generation are highly variable."
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Additive linear effect on CL with coefficient theta7 = 0.0377 L/h per",
        "(mL/min/1.73 m^2) (Table 4); no centering / no reference value (the",
        "raw GFR multiplies the slope directly). Development cohort median 115",
        "mL/min/1.73 m^2 (IQR 94-139). Pooled with creatinine-based renal",
        "function under the CRCL canonical per the 2026-06-17 register policy",
        "(BSA-normalised renal function (mL/min/1.73 m^2): creatinine-based",
        "estimate OR tracer-measured GFR); the assay method (99mTc-DTPA) is",
        "documented here so future reviewers can trace the source."
      ),
      source_name        = "GFR (99mTc-DTPA tracer plasma clearance)"
    )
  )

  covariatesDataExcluded <- list(
    HT = list(
      description = "Body height. Screened during covariate model building; not retained in the final model.",
      units       = "cm",
      type        = "continuous",
      notes       = "Univariately reduced OFV by > 6.63 on CL, V, and k12 but was not retained after backward elimination once WT was in the model (Results 'Development of a covariate population pharmacokinetic model')."
    ),
    BSA = list(
      description = "Body surface area. Screened during covariate model building; not retained.",
      units       = "m^2",
      type        = "continuous",
      notes       = "Univariately reduced OFV by > 6.63 on CL but was not retained after backward elimination once WT was in the model."
    ),
    WT_ALLO = list(
      description = "Body weight raised to the 0.75 power (allometric size). Screened; not retained.",
      units       = "kg^0.75",
      type        = "continuous",
      notes       = "WT^0.75 was tested as an alternative size descriptor on CL and V but the linear WT form gave the better fit and was retained."
    ),
    AGE = list(
      description = "Age in years. Screened during covariate model building; not retained.",
      units       = "years",
      type        = "continuous",
      notes       = "Univariately reduced OFV by > 6.63 on CL but was not retained after backward elimination once WT was in the model. Defines population scope (0.3-18 yrs) rather than entering as an effect."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male). Screened during covariate model building; not retained.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained. Development cohort 28/11 male/female (28% female)."
    ),
    PRIOR_TBI = list(
      description = "Binary indicator for prior total-body irradiation as part of BMT conditioning. Screened; not retained.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Univariately reduced OFV by > 6.63 on CL but was not retained after",
        "backward elimination (Results 'Development of a covariate population",
        "pharmacokinetic model'). The Discussion notes that the previous",
        "two-stage analysis [9] found TBI significant on CL but the population",
        "analysis did not. Development cohort 16/39 received prior TBI."
      )
    ),
    PRIOR_BUSULFAN = list(
      description = "Binary indicator for prior busulphan as part of BMT conditioning. Screened; not retained.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Not significant in the preliminary screening phase (did not reduce OFV by 6.63 on any PK parameter). Development cohort 7/39 received prior busulphan."
    ),
    DOSE_GROUP = list(
      description = "Mass-normalised dose group (mg/m^2). Screened; not retained.",
      units       = "mg/m^2",
      type        = "continuous",
      notes       = "Tested as a covariate to assess dose-linearity. Not significant in the preliminary screening phase, supporting linear pharmacokinetics of melphalan over the 30-180 mg/m^2 range."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 59L,
    n_studies      = 1L,
    n_subjects_development = 39L,
    n_subjects_validation  = 20L,
    n_observations         = 848L,
    n_observations_development = 571L,
    n_observations_validation  = 277L,
    age_range      = "0.3-18 years",
    age_median     = "5.4 months (development cohort median; values reported in months in Table 1)",
    weight_range   = "7.7-104 kg (overall study; Discussion)",
    weight_median  = "18.8 kg (development cohort, IQR 13.5-28.1)",
    height_median  = "110 cm (development cohort, IQR 92-137)",
    bsa_median     = "0.76 m^2 (development cohort, IQR 0.60-1.0)",
    crcl_median    = "115 mL/min/1.73 m^2 (development cohort, IQR 94-139; 99mTc-DTPA tracer GFR)",
    sex_female_pct = 28L,
    race_ethnicity = "Not reported (single-centre Australian cohort).",
    disease_state  = paste(
      "Paediatric autologous or allogeneic blood or marrow transplant",
      "recipients with malignant diseases. Diagnoses (development cohort):",
      "neuroblastoma (13), Acute lymphoblastic leukaemia (6), Acute myeloid",
      "leukaemia (5), rhabdomyosarcoma (5), non-Hodgkin's lymphoma (4), soft",
      "tissue sarcoma (2), Ewing's sarcoma (1), hepatoblastoma (1),",
      "retinoblastoma (1), mediastinal large-cell lymphoma (1), and others",
      "(Table 1)."
    ),
    dose_range     = paste(
      "Single high dose 140 or 180 mg/m^2, or divided dose schedule (3 days",
      "of 70 mg/m^2 or 4 days of 30 mg/m^2). All doses administered as a",
      "15-min intravenous infusion with double maintenance fluids."
    ),
    bmt_conditioning_comedications = paste(
      "Carboplatin (CPT; 17/39 development cohort, on each of 5 days preceding",
      "melphalan via the Calvert formula targeting AUC = 4 mg/mL/min);",
      "total body irradiation (TBI; 16/39); busulphan (7/39).",
      "Carboplatin shifts CL by an additive -3.17 L/h (Table 4); TBI and",
      "busulphan were screened and dropped."
    ),
    regions        = "Australia (The Children's Hospital at Westmead, Sydney; 1994-2003).",
    notes          = paste(
      "Prospective single-centre paediatric BMT cohort. Sampling schedule:",
      "pre-infusion, then 0, 5, 10, 15, 20, 30, 40, 50 min and 1, 2, 3, 4,",
      "6, 12, 24 h after end of infusion; median 15 samples per subject in",
      "the development cohort. Melphalan was assayed by HPLC with LOQ 0.5",
      "mg/L. Random allocation to development (n=39) vs validation (n=20)",
      "subcohorts with no significant baseline differences (Table 1).",
      "NONMEM v5.1.1, FOCE with eta-epsilon interaction."
    )
  )

  ini({
    # Structural population PK parameters (Nath 2007 Table 4 final covariate
    # population PK model). The CL model is a linear additive function of
    # body weight, prior-carboplatin status, and tracer-measured GFR; the V
    # model is a linear additive function of body weight; k12 and k21 are
    # estimated directly (rate constants in 1/h, not as Q / Vc / Vp).
    #
    # Structural model:
    #   CL = theta5 * WT + theta6 * CPT + theta7 * GFR    (theta1 fixed to 0)
    #   Vc = theta2 + theta8 * WT
    #   k12 = theta3
    #   k21 = theta4

    # Body-weight coefficient on CL (Table 4 theta5).
    lcl <- log(0.34);   label("Body-weight coefficient on CL (theta5, L/h per kg)")  # Table 4 theta5 = 0.34 (95% CI 0.26, 0.42)

    # Intercept (constant term) on Vc (Table 4 theta2).
    lvc <- log(1.12);   label("Constant term in Vc (theta2, L)")                      # Table 4 theta2 = 1.12 (95% CI 0.11, 2.13)

    # Body-weight coefficient on Vc (Table 4 theta8).
    e_wt_vc <- 0.178;   label("Body-weight coefficient on Vc (theta8, L per kg)")     # Table 4 theta8 = 0.178 (95% CI 0.120, 0.236)

    # Prior-carboplatin additive effect on CL (Table 4 theta6).
    e_prior_carboplatin_cl <- -3.17; label("Prior-carboplatin effect on CL (theta6, L/h)")  # Table 4 theta6 = -3.17 (95% CI -4.62, -1.72)

    # Tracer-measured GFR additive linear effect on CL (Table 4 theta7).
    e_crcl_cl <- 0.0377;             label("GFR effect on CL (theta7, L/h per (mL/min/1.73 m^2))")  # Table 4 theta7 = 0.0377 (95% CI 0.021, 0.054)

    # Intercompartmental rate constants (Nath 2007 Table 4 theta3, theta4).
    # The paper estimates k12 and k21 directly as rate constants (1/h) rather
    # than via Q / Vc and Q / Vp. theta1 (constant in CL model) is reported in
    # Table 4 as "Fixed to zero" and is not encoded as a separate parameter.
    lk12 <- log(1.70);  label("Intercompartmental rate constant k12 (theta3, 1/h)")  # Table 4 theta3 = 1.70 (95% CI 1.16, 2.24)
    lk21 <- log(1.84);  label("Intercompartmental rate constant k21 (theta4, 1/h)")  # Table 4 theta4 = 1.84 (95% CI 1.39, 2.29)

    # Inter-individual variability (Nath 2007 Table 4 IIV column). The paper
    # explicitly defines %CV as sqrt(eta variance) * 100 -- i.e. omega is
    # reported directly as CV (the simple-approximation form), not as the
    # log-normal-true-CV form. Therefore omega^2 = (CV/100)^2.
    #   IIV CL  : 27.3% CV -> omega^2 = 0.273^2 = 0.074529
    #   IIV V   : 33.8% CV -> omega^2 = 0.338^2 = 0.114244
    #   IIV k12 : 52.2% CV -> omega^2 = 0.522^2 = 0.272484
    #   IIV k21 : 61.7% CV -> omega^2 = 0.617^2 = 0.380689
    # Exponential (log-normal) random-effects model: theta_i = theta * exp(eta_i).
    # The paper does not report off-diagonal correlations for the final model,
    # so the IIV block is diagonal.
    etalcl  ~ 0.074529   # Table 4 IIV CL  27.3% CV (95% CI 21, 32)
    etalvc  ~ 0.114244   # Table 4 IIV V   33.8% CV (95% CI 20, 43)
    etalk12 ~ 0.272484   # Table 4 IIV k12 52.2% CV (95% CI 30, 68)
    etalk21 ~ 0.380689   # Table 4 IIV k21 61.7% CV (95% CI 28, 83)

    # Residual variability (Nath 2007 Table 4 residual variability rows).
    # Combined additive + proportional error: Y = Y_hat * (1 + e1) + e2.
    propSd <- 0.093;     label("Proportional residual error (fraction)")   # Table 4 proportional 9.3% CV (95% CI 7, 11)
    addSd  <- 0.00731;   label("Additive residual error (mg/L)")            # Table 4 additive 0.00731 mg/L (95% CI -0.004, 0.019)
  })

  model({
    # Individual PK parameters. The structural model from Table 4 is preserved
    # exactly: CL has WT-linear, prior-carboplatin-additive, and GFR-additive
    # components; Vc has an intercept plus WT-linear term; k12 and k21 are
    # rate constants (1/h) with their own log-normal IIV.
    cl  <- (exp(lcl) * WT + e_prior_carboplatin_cl * PRIOR_CARBOPLATIN + e_crcl_cl * CRCL) * exp(etalcl)
    vc  <- (exp(lvc) + e_wt_vc * WT) * exp(etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)

    kel <- cl / vc

    # Two-compartment IV PK. Melphalan is administered as a 15-min IV infusion
    # in the source paper; the library model does not hard-code the infusion
    # duration, so users specify rate / dur per dose in their event table.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central                 - k21 * peripheral1

    # Plasma concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
