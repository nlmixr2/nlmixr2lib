Hanke_2020_cimetidine <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous and fasted oral",
    "cimetidine in adults (Hanke 2020, Electronic Supplementary Material",
    "section 3), fitted in NONMEM to mean concentration-time profiles from 25",
    "published studies (100-800 mg). The fasted-state double peak is described",
    "by splitting each oral dose into two portions that share one first-order",
    "absorption rate constant: the first portion (fraction 1 - VF2, typical",
    "71.2 percent) is absorbed from depot without delay and the second",
    "(fraction VF2, typical 28.8 percent) from depot2 after a lag time of",
    "1.54 h. Total oral bioavailability is 90.2 percent and elimination is",
    "first order from the central compartment. Random effects are",
    "between-study variability on VF2 (logit scale), the second-portion lag",
    "time, CL and Vc; there are no covariates. This is the population-PK",
    "analysis the paper used to derive the split-dose input for its PK-Sim",
    "whole-body PBPK model; the PBPK layer itself is not reproduced here."
  )
  reference <- paste(
    "Hanke N, Turk D, Selzer D, Ishiguro N, Ebner T, Wiebe S, Muller F,",
    "Stopfer P, Nock V, Lehr T. A Comprehensive Whole-Body Physiologically",
    "Based Pharmacokinetic Drug-Drug-Gene Interaction Model of Metformin and",
    "Cimetidine in Healthy Adults and Renally Impaired Individuals.",
    "Clin Pharmacokinet. 2020;59(11):1419-1431. doi:10.1007/s40262-020-00896-w.",
    "Population-PK parameters from Electronic Supplementary Material Table",
    "S3.4.1; structure and residual-error model from the NONMEM control",
    "stream in ESM section 3.5.",
    sep = " "
  )
  vignette <- "Hanke_2020_cimetidine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # An oral dose is given as TWO dose records of the full total cimetidine
  # dose (TCD), one to depot and one to depot2; f(depot) and f(depot2) then
  # route (1 - VF2) * FTOT and VF2 * FTOT of it, as in the NONMEM $PK block
  # (F1 = VF1 * FTOT, F2 = VF2 * FTOT). Intravenous doses go to central.
  dosing <- c("depot", "depot2", "central")

  compartmentData <- list(
    depot = list(analyte = "cimetidine", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "cimetidine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cimetidine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cimetidine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_studies = 25,
    n_subjects = 215,
    age_range = "19-80 years",
    disease_state = "healthy volunteers and peptic ulcer patients",
    dose_range = "100-800 mg; intravenous (bolus or 2-30 min infusion) and oral in the fasted state (solution, capsule, tablet)",
    regions = "not reported (published literature studies, 1975-1998)",
    unit_of_analysis = "study-mean concentration-time profiles, plus 3 published representative individual profiles",
    notes = paste(
      "Hanke 2020 ESM Table S3.3.1 lists 25 studies (9 intravenous, 16 fasted",
      "oral) with 215 subject-arm entries; subjects overlap between the",
      "intravenous and oral arms of crossover studies, so the number of",
      "distinct individuals is smaller. ESM section 3.3.1: 'Average cimetidine",
      "concentration-time profiles were used for model building, together",
      "with 3 individual profiles that were published as representative",
      "examples of larger study populations.' Section 3.3.2 therefore",
      "interprets the random effects as interstudy variability (ISV) rather",
      "than between-subject variability. ESM Table S4.2.1 marks the studies",
      "as healthy-volunteer or peptic-ulcer-patient cohorts; Jonsson 1982",
      "(400 mg IV, n = 19, 43-80 years) appears only in Table S3.3.1. Some",
      "studies measured whole blood rather than plasma (Table S4.5.2); the",
      "population-PK analysis pooled them without a matrix term.",
      sep = " "
    )
  )

  ini({
    # ========================================================================
    # Absorption and bioavailability (ESM Table S3.4.1; control stream
    # section 3.5 $PK). KA1 is shared by both oral portions; ALAG2 delays
    # only the second portion; FTOT is the absolute bioavailability applied
    # to the whole oral dose.
    # ========================================================================
    lka <- log(0.753); label("First-order absorption rate constant, both oral portions (1/h)") # Table S3.4.1 'KA' = 0.753 1/h (RSE 5.1%)
    lfdepot <- log(0.902); label("Absolute oral bioavailability FTOT (fraction)") # Table S3.4.1 'FTOT' = 90.2% (RSE 6.3%)
    logitfrac <- log(0.288 / (1 - 0.288)); label("Logit of the fraction of the total oral dose attributed to the second portion VF2 (unitless)") # Table S3.4.1 'VF2' = 0.288 (RSE 19.4%); $PK PHI_2 = LOG(VF2/(1-VF2))
    ltlag2 <- log(1.54); label("Lag time of the second oral portion ALAG2 (h)") # Table S3.4.1 'ALAG' = 1.54 h (RSE printed as 0%; $THETA(2) carries no FIX flag)

    # ========================================================================
    # Disposition: two-compartment model with first-order elimination from
    # central (ESM section 3.4; $DES K30 = CL/V3, K34 = Q/V3, K43 = Q/V4).
    # ========================================================================
    lcl <- log(41.2); label("Clearance (L/h)") # Table S3.4.1 'CL' = 41.2 L/h (RSE 6%)
    lvc <- log(32.6); label("Central volume of distribution (L)") # Table S3.4.1 'V3' = 32.6 L (RSE 11.8%)
    lq <- log(45.4); label("Intercompartmental clearance (L/h)") # Table S3.4.1 'Q' = 45.4 L/h (RSE 7.3%)
    lvp <- log(46); label("Peripheral volume of distribution (L)") # Table S3.4.1 'V4' = 46 L (RSE 4.8%)

    # ========================================================================
    # Interstudy variability. ESM section 3.3.2: 'ISVs were modeled
    # exponentially'; the control stream puts ETA(2)-ETA(4) on CL, V3 and
    # ALAG2 as THETA * EXP(ETA) and ETA(1) on the logit of VF2 (PHI_2 + ETA(1)).
    # Table S3.4.1 reports each ISV as a %CV; the ESM uses the log-normal
    # relation CV = sqrt(exp(omega^2) - 1) elsewhere (Table S9.0.1 footnote:
    # '35 % CV ... (= 1.40 GSD)'), so omega^2 = log(1 + CV^2). The same column
    # and conversion are applied to the logit-scale VF2 row. The $OMEGA block
    # is diagonal and its printed numbers are initial estimates, not finals.
    # ========================================================================
    etalogitfrac ~ log(1 + 0.848^2) # Table S3.4.1 'ISV VF2' = 84.8 %CV (RSE 28.8%)
    etaltlag2 ~ log(1 + 0.205^2) # Table S3.4.1 'ISV ALAG' = 20.5 %CV (RSE 18%)
    etalcl ~ log(1 + 0.216^2) # Table S3.4.1 'ISV CL' = 21.6 %CV (RSE 20.9%)
    etalvc ~ log(1 + 0.387^2) # Table S3.4.1 'ISV V3' = 38.7 %CV (RSE 15.7%)

    # ========================================================================
    # Residual error: $ERROR Y = IPRED + IPRED * EPS(1) + EPS(2), two
    # independent epsilons (combined proportional plus additive). The
    # additive $SIGMA 0.000323 carries the FIX flag, so its control-stream
    # value is the final value; SD = sqrt(0.000323) = 0.017972 mg/L. Table
    # S3.4.1 prints it as 'Add RE +- 1.8', i.e. scaled by 100 like the
    # percent-valued proportional row.
    # ========================================================================
    propSd <- 0.122; label("Proportional residual SD (fraction)") # Table S3.4.1 'Prop RE' = 12.2% (RSE 7.8%)
    addSd <- fixed(0.017972); label("Additive residual SD (mg/L)") # $SIGMA EPS(2) = 0.000323 with FIX flag; sqrt(0.000323) = 0.017972; Table S3.4.1 'Add RE' = 1.8
  })

  model({
    # ---- Individual parameters (control stream $PK) ------------------------
    ka <- exp(lka)
    tlag2 <- exp(ltlag2 + etaltlag2)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    fdepot <- exp(lfdepot)

    # VF2_2 = EXP(PHI_2 + ETA(1)) / (1 + EXP(PHI_2 + ETA(1)))
    frac <- expit(logitfrac + etalogitfrac)

    # ---- Micro-constants ($DES) --------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- ODE system ($DES DADT(1)-DADT(4)) ---------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka * depot2
    d/dt(central) <- ka * depot + ka * depot2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- Bioavailability and lag ($PK F1, F2, ALAG2) -----------------------
    f(depot) <- (1 - frac) * fdepot
    f(depot2) <- frac * fdepot
    alag(depot2) <- tlag2

    # ---- Observation and residual error ($ERROR; S3 = V3) ------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
