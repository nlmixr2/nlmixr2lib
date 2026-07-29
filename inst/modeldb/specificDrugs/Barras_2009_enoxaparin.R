Barras_2009_enoxaparin <- function() {
  description <- paste(
    "Two-compartment first-order absorption population PK model",
    "of anti-factor Xa activity in 118 adults (PE / DVT / ACS /",
    "atrial fibrillation) receiving subcutaneous enoxaparin",
    "treatment doses (1 mg/kg BID by total or lean body weight,",
    "1.5 mg/kg BID for LBW-based obese dosing) under conventional",
    "vs lean-body-weight-and-renal-function individualised dosing",
    "(Barras 2009 randomised controlled trial). CL is a composite",
    "renal + non-renal model with LBW substituted into the",
    "Cockcroft-Gault CrCl equation; central volume scales linearly",
    "with LBW. The paper additionally reports a three-category",
    "proportional-odds bleeding / bruising adverse-event PD model",
    "with logit(P[S<=1]) = 2.83 - 2.75*(Age/61) - 0.536*(cAUC/23)",
    "and logit(P[S<=2]) = logit(P[S<=1]) + 2.05, driven by patient",
    "Age and cumulative AUC (cAUC) of anti-Xa activity from first",
    "dose to event. The proportional-odds PD layer is NOT encoded",
    "in this model file -- it requires canonical parameter names",
    "for cumulative-logit / proportional-odds PD models that are",
    "not yet registered in references/parameter-names.md. The PD",
    "equation is reproduced in the validation vignette, where it",
    "is applied deterministically to cAUC values derived from the",
    "simulated PK profile (see vignette Source trace and",
    "Assumptions and deviations sections)."
  )
  reference <- paste(
    "Barras MA, Duffull SB, Atherton JJ, Green B. Modelling the",
    "occurrence and severity of enoxaparin-induced bleeding and",
    "bruising events. Br J Clin Pharmacol. 2009;68(5):700-711.",
    "doi:10.1111/j.1365-2125.2009.03518.x"
  )
  vignette <- "Barras_2009_enoxaparin"
  units    <- list(time = "hour", dosing = "IU", concentration = "IU/mL")

  covariateData <- list(
    LBM = list(
      description        = "Lean body weight calculated by the Janmahasatian (2005) formula from total body weight, height, and sex",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Stored under the canonical LBM column per",
        "inst/references/covariate-columns.md; LBW is a registered",
        "source alias of LBM. Reference value 55 kg (Barras 2009",
        "population median LBW per Table 3, range 30-86 kg in the",
        "n=118 PK cohort, range 34-80 kg in the n=103 PD cohort).",
        "LBW enters the model two ways: (a) as a linear scalar on",
        "the central volume Vc (Eq. 3: Vc = 3.43 * (LBW / 55)); and",
        "(b) as the body-size descriptor substituted for total",
        "weight in the Cockcroft-Gault CrCl equation that feeds the",
        "renal arm of CL (Methods, Covariate model: 'CLCR was",
        "calculated according to the C-G equation [Cockcroft and",
        "Gault 1976] and where Wt, IBW and LBW were included",
        "separately into the equation', with the LBW-substituted",
        "form retained in the final model). LBW is also a linear",
        "scalar on the non-renal arm of CL (Eq. 2: 0.42 *",
        "(LBW / 55))."
      ),
      source_name        = "LBW"
    ),
    CRCL = list(
      description        = paste(
        "Creatinine clearance estimated by the Cockcroft-Gault",
        "(1976) equation with LBW (Janmahasatian 2005 formula)",
        "substituted for total body weight as the body-size",
        "descriptor. Raw mL/min, NOT BSA-normalised."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference value 70 mL/min (Barras 2009 population median",
        "for the LBW-substituted C-G CrCl per Table 3, range",
        "10-170 mL/min in the n=118 PK cohort, range 17-166 mL/min",
        "in the n=103 PD cohort). Stored under the canonical CRCL",
        "column per inst/references/covariate-columns.md; this is",
        "the same operational form as Delattre 2010 amikacin (raw",
        "C-G mL/min, not BSA-normalised) and is distinct from the",
        "MDRD / CKD-EPI eGFR forms used in many recent mAb models.",
        "Used as a linear scalar on the renal arm of CL only",
        "(Eq. 2: 0.3 * (CRCL / 70)); the non-renal arm scales with",
        "LBW (not CRCL)."
      ),
      source_name        = "CLCR (LBW)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 118L,
    n_studies       = 1L,
    age_range       = "23-91 years",
    age_median      = "61 years",
    weight_range    = "43-120 kg",
    weight_median   = "77 kg",
    lbw_range       = "30-86 kg",
    lbw_median      = "55 kg",
    ibw_range       = "39-83 kg",
    ibw_median      = "64 kg",
    sex_female_pct  = 38,
    race_ethnicity  = "Not reported (Australian university-hospital cohort)",
    disease_state   = paste(
      "Hospitalised adults treated for pulmonary embolism, deep",
      "vein thrombosis, acute coronary syndrome, or atrial",
      "fibrillation; subjects with severe renal impairment",
      "(CLCR < 10 mL/min by Cockcroft-Gault), pregnancy, recent",
      "warfarin / heparin (within 7 days), liver enzymes > 2x",
      "normal, INR > 1.2, aPTT > 60 s, or haemodialysis were",
      "excluded. Co-administration of warfarin in 31% and",
      "antiplatelet drugs (aspirin / clopidogrel) in 42% per Table",
      "3."
    ),
    dose_range      = paste(
      "Subcutaneous enoxaparin twice daily under conventional",
      "(product label) or individualised dosing. Individualised",
      "arm: 1 mg/kg total body weight BID for subjects < 100 kg,",
      "1.5 mg/kg LBW BID for subjects >= 100 kg (considered obese,",
      "LBW-based dose); subjects with CLCR < 50 mL/min dose-",
      "reduced after 48 h based on renal function. Conventional",
      "arm: product-label total-body-weight-based dosing as",
      "prescribed. Mean (SD) treatment duration 3.5 +/- 2.3 days.",
      "The mg dose was converted to anti-factor Xa international",
      "units at data-assembly time using the standard enoxaparin",
      "specific activity (approximately 100 IU anti-Xa per mg) so",
      "the model parameters describe anti-Xa activity in IU/mL,",
      "matching the chromogenic assay readout."
    ),
    regions         = "Australia (Royal Brisbane and Women's Hospital and Mater Health Services)",
    renal_function  = paste(
      "LBW-substituted C-G CrCl median 70 mL/min (range 10-170)",
      "in the n=118 PK cohort; 73 mL/min (range 17-166) in the",
      "n=103 PD cohort. Total-weight C-G CrCl median 85 mL/min",
      "(range 15-244). Subjects with severe renal impairment",
      "(CLCR < 30 mL/min) represented only 7% of the population",
      "and contributed 5% of the anti-Xa observations per the",
      "Discussion."
    ),
    co_medication   = "Concomitant warfarin in 31% of subjects; concomitant antiplatelet drugs (aspirin and / or clopidogrel) in 42% (per Table 3).",
    notes           = paste(
      "Prospective randomised controlled trial: 118 subjects",
      "contributed 349 anti-Xa activity observations (mean 3 per",
      "subject, range 1-4), with 93% of samples in the first 48 h",
      "of therapy. Samples scheduled by D-optimal design",
      "(predose; 15-30 min; 60-120 min; 180-300 min post dose).",
      "Anti-Xa activity measured by the chromogenic STA-Rotachrom",
      "assay (Diagnostica Stago) at the Queensland Health",
      "Pathology Service (LOQ 0.1 IU/mL, inter-assay precision 5.9",
      "%CV). NONMEM V.1.1 with G77 FORTRAN, FOCE-I estimation.",
      "Six samples (1.7%) below LOQ excluded in the initial model-",
      "building phase; addition at half-LOQ did not affect the",
      "final parameter estimates (per Methods 'Data below the",
      "level of quantification'). Bootstrap 873 / 1000 successful",
      "minimisations supported the asymptotic 95% CIs (Table 4).",
      "A separate n=103 subset (those with bleeding / bruising",
      "assessment beyond baseline) was used for the proportional-",
      "odds PD model."
    )
  )

  ini({
    # ------------------------------------------------------------
    # Structural PK parameters at the population median
    # (LBM = 55 kg, CRCL = 70 mL/min). Barras 2009 Table 4
    # 'Covariate model (95% CI)' column; final model values match
    # the equations on page 705 (Eq. 2 for CL, Eq. 3 for Vc).
    # ------------------------------------------------------------
    lka <- log(0.26); label("Absorption rate constant (Ka, 1/h)")              # Barras 2009 Table 4: Ka = 0.26 (95% CI 0.17, 0.35)
    lvc <- log(3.43); label("Apparent central volume at LBM = 55 kg (Vc/F, L)") # Barras 2009 Table 4 / Eq. 3: Vc = 3.43 (95% CI 2.20, 4.66)
    lvp <- log(5.77); label("Apparent peripheral volume (Vp/F, L)")            # Barras 2009 Table 4: Vp = 5.77 (95% CI 0.75, 10.8)
    lq  <- log(0.31); label("Apparent inter-compartmental clearance (Q/F, L/h)") # Barras 2009 Table 4: Q = 0.31 (95% CI 0.15, 0.47)

    # ------------------------------------------------------------
    # Composite renal + non-renal CL (Barras 2009 Eq. 2). The
    # paper parameterises total CL/F at the population median as
    # CL = 0.3 * (CRCL / 70) + 0.42 * (LBM / 55), i.e. a renal arm
    # scaled by LBW-substituted C-G CrCl and a non-renal arm
    # scaled by LBW. The renal arm intercept 0.3 L/h corresponds
    # to a typical subject with CRCL = 70 mL/min; the non-renal
    # arm intercept 0.42 L/h corresponds to LBM = 55 kg. Log-
    # transformed separately so the two arms can be added inside
    # model() with a single shared IIV multiplier on the total.
    # ------------------------------------------------------------
    lcl_renal  <- log(0.30); label("Renal arm of CL at CRCL = 70 mL/min (L/h)")  # Barras 2009 Table 4: CL_renal = 0.3 (95% CI 0.14, 0.45)
    lcl_nonren <- log(0.42); label("Non-renal arm of CL at LBM = 55 kg (L/h)")   # Barras 2009 Table 4: CL_nonren = 0.42 (95% CI 0.24, 0.61)

    # ------------------------------------------------------------
    # Inter-individual variability. Barras 2009 Table 4 reports
    # IIV as %CV from a log-normal (exponential) parameterisation
    # on each PK parameter; omega^2 = log(CV^2 + 1).
    #   CL: 37.8% CV -> log(0.378^2 + 1) = 0.13342
    #   Vc: 35.6% CV -> log(0.356^2 + 1) = 0.11907
    #   Ka: 30.3% CV -> log(0.303^2 + 1) = 0.08763
    # The Results 'PK analysis -- Model building' paragraph states
    # 'Vc and CL allowed to co-vary' but the covariance value is
    # not reported in Table 4 (only the diagonal CV%s). This
    # implementation uses diagonal IIVs; see the vignette
    # Assumptions and deviations section.
    # ------------------------------------------------------------
    etalcl ~ 0.13342 # log(0.378^2 + 1); 37.8% CV on CL/F (95% CI 22.3, 48.0), Barras 2009 Table 4
    etalvc ~ 0.11907 # log(0.356^2 + 1); 35.6% CV on Vc/F (95% CI 7, 51.0),    Barras 2009 Table 4
    etalka ~ 0.08763 # log(0.303^2 + 1); 30.3% CV on Ka   (95% CI 5.3, 42.5),  Barras 2009 Table 4

    # ------------------------------------------------------------
    # Residual error (Barras 2009 Table 4: additive error model,
    # epsilon = 0.09 IU/mL (95% CI 0.086, 0.1)). The Methods
    # 'Base heterogeneity and residual error model' paragraph and
    # the Results 'PK analysis -- Model building' paragraph both
    # state that an additive (not proportional or combined) error
    # model best described residual variability on anti-Xa
    # activity.
    # ------------------------------------------------------------
    addSd <- 0.09; label("Additive residual error on anti-Xa activity (IU/mL)") # Barras 2009 Table 4: epsilon = 0.09 (95% CI 0.086, 0.1)
  })
  model({
    # ------------------------------------------------------------
    # Composite renal + non-renal CL (Barras 2009 Eq. 2 with IIV
    # added). The two arms are added on the linear scale (paper
    # parameterisation); IIV enters as a multiplicative log-normal
    # factor on the total. Renal arm scales linearly with CRCL
    # (LBW-substituted C-G); non-renal arm scales linearly with
    # LBM. Vc scales linearly with LBM (Eq. 3).
    # ------------------------------------------------------------
    cl <- (exp(lcl_renal) * (CRCL / 70) + exp(lcl_nonren) * (LBM / 55)) * exp(etalcl)
    vc <- exp(lvc + etalvc) * (LBM / 55)
    vp <- exp(lvp)
    q  <- exp(lq)
    ka <- exp(lka + etalka)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Dose in IU, volumes in L. central/vc has units IU/L; divide
    # by 1000 to convert to IU/mL (the paper's reporting unit for
    # anti-Xa activity).
    Cc <- central / (vc * 1000)
    Cc ~ add(addSd)
  })
}
