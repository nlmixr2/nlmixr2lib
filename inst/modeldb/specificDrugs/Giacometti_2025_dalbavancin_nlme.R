Giacometti_2025_dalbavancin_nlme <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for dalbavancin, estimated by nonlinear",
    "mixed-effects (Monolix, maximum likelihood) in 218 adults undergoing therapeutic drug",
    "monitoring during long-acting dalbavancin therapy for skin and soft tissue, bone and joint,",
    "vascular prosthetic and endocarditis staphylococcal infections. This is the NLME arm of a",
    "three-way methodological comparison (naive-pooled least squares vs NLME vs a Neural ODE); the",
    "companion least-squares fit to the SAME cohort is modellib('Giacometti_2025_dalbavancin_nlls').",
    "NO COVARIATES are carried: age, height, weight, sex and serum creatinine were all collected",
    "and screened, but Appendix B reports that adding them did not improve the fit (KS p = 0.45),",
    "so the authors deliberately retained the covariate-free model as the more robust and scalable",
    "one. INTER-INDIVIDUAL VARIABILITY AND RESIDUAL ERROR ARE NOT REPORTED: Table 4 prints only",
    "the four structural fixed effects, and its plus-or-minus figures are the standard deviation",
    "ACROSS the six cross-validation folds, not omegas and not RSEs. Simulations from this file are",
    "therefore typical-value only. Predictions are compared in the source against an efficacy",
    "threshold of 8.04 mg/L total plasma dalbavancin.",
    sep = " "
  )
  reference <- paste(
    "Giacometti T, Rocchi E, Cojutti PG, Magnani F, Remondini D, Pea F, Castellani G.",
    "Leveraging Neural ODEs for Population Pharmacokinetics of Dalbavancin in Sparse Clinical",
    "Data. Entropy. 2025;27(6):602. doi:10.3390/e27060602",
    sep = " "
  )
  vignette <- "Giacometti_2025_dalbavancin"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Giacometti 2025 Figure 2 is the schematic: drug is
  # administered into the central compartment, which is the only compartment
  # that eliminates.
  compartmentData <- list(
    central     = list(analyte = "dalbavancin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dalbavancin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # No covariate enters model(). See covariatesDataExcluded below.
  covariateData <- NULL

  # Every covariate the study collected was screened on the NLME parameters and
  # NONE was retained. Giacometti 2025 Methods 2.1: "All the dataset's
  # covariates collected were: age, height, weight, sex and serum creatine
  # concentration." Methods 2.2: "The covariates were not included in the NLME
  # and two-compartment models because, when tested, their inclusion did not
  # improve the results; see Appendix B." Appendix B: "The KS test for the two
  # distributions yields a p-value of 0.45, indicating no significant difference
  # between them. Consequently, while we acknowledge that the exclusion of
  # clinically relevant covariates may limit the physiological interpretability
  # of the NLME model, we preferred to select the less complex model since it
  # represents a robust and scalable version of the model and does not rely on
  # the specific choices of the covariates to include for each parameter, which
  # vary a lot depending on the specific iterations of the cross-validation."
  # Documentation only -- none of these is referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on the NLME parameters and not retained (Methods 2.2, Appendix B). Table 1: range 18-92 y, mean 64 +/- 16 y."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened on the NLME parameters and not retained (Methods 2.2, Appendix B). Table 1: range 145-190 cm, mean 171 +/- 9 cm."
    ),
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened on the NLME parameters and not retained (Methods 2.2, Appendix B), so this model",
        "applies NO allometric scaling -- unlike the sibling dalbavancin model",
        "modellib('Baiardi_2025_dalbavancin'), which fixes allometric exponents of 0.75 on the",
        "clearances and 1 on the volumes. Table 1: range 40-140 kg, mean 77 +/- 16 kg."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on the NLME parameters and not retained (Methods 2.2, Appendix B). Methods 2.1: 145 males and 73 females of 218 patients (33.5% female). Source column 'sex'; the paper does not state its coding direction, which is immaterial here because the covariate is not used."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Screened on the NLME parameters and not retained (Methods 2.2, Appendix B). Reported in SI",
        "units: Table 1 gives 'Creatine conc. (umol/L)' range 19-411, mean 94 +/- 48 (the paper",
        "spells it 'creatine' in Table 1 and 'serum creatine concentration' in Methods 2.1;",
        "creatinINE is meant -- 94 umol/L is 1.06 mg/dL, a normal serum creatinine, whereas the",
        "quoted range is far outside any plasma creatine concentration). No derived renal-function",
        "measure (eGFR, creatinine clearance) is reported, so this model carries no renal covariate",
        "-- contrast modellib('Cojutti_2024_dalbavancin'), from the same clinical unit, which does",
        "retain CKD-EPI eGFR on CL."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 218L,
    n_studies      = 1L,
    n_observations = "669 dalbavancin plasma concentration measurements over 703 recorded administrations (Methods 2.1)",
    age_range      = "18-92 years",
    age_mean       = "64 +/- 16 years (Table 1)",
    height_range   = "145-190 cm (mean 171 +/- 9 cm, Table 1)",
    weight_range   = "40-140 kg",
    weight_mean    = "77 +/- 16 kg (Table 1)",
    sex_female_pct = 33.5,
    race_ethnicity = "Not reported; single-centre Italian cohort.",
    disease_state  = paste(
      "Adults receiving dalbavancin under therapeutic drug monitoring. Methods 2.1: dalbavancin is",
      "'a long-acting antibiotic used for the treatment of skin and soft tissue infections and as a",
      "second-line agent in patients with staphylococcal infections, such as bone and joint",
      "infections, vascular prosthetic joint infections and endocarditis.' The paper does not",
      "tabulate the per-infection-type counts."
    ),
    renal_function = "Reported only as serum creatinine: range 19-411 umol/L, mean 94 +/- 48 umol/L (Table 1). No eGFR or creatinine clearance is derived, and no renal covariate is carried by the model.",
    dose_range     = "350-1500 mg intravenously, 703 administrations across the 218 patients (Methods 2.1). The paper does not report the dosing intervals or the infusion duration; the NLME model was configured in Monolix as 'two-compartment distribution, infusion, no delay and linear elimination' (Methods 2.2).",
    regions        = "Italy (IRCCS Azienda Ospedaliero-Universitaria di Bologna)",
    notes          = paste(
      "Retrospective single-centre therapeutic-drug-monitoring cohort, April 2021 to December 2024",
      "(Ethics Committee 897/2021/Oss/AOUBo). The reported estimates are the MEAN over the six folds",
      "of a 6-fold cross-validation, and the plus-or-minus figures in Table 4 are the standard",
      "deviation ACROSS those folds -- they are neither standard errors nor between-subject",
      "variances. Estimation was maximum likelihood in Monolix. Sampling is explicitly sparse and",
      "irregular: Methods 2.2 notes that 'frequently, patients undergo multiple drug administrations",
      "before any pharmacokinetic measurement is recorded, with some measurements occurring only",
      "after several doses.' The paper reports no goodness-of-fit statistic for the NLME arm beyond",
      "the residual distributions of Figures 4-6 and the Kolmogorov-Smirnov comparisons of Tables",
      "2-3; the R^2 values of Appendix C Table A2 (0.83-0.87) belong to the Neural ODE arm. The",
      "Neural ODE that is the paper's headline model is NOT extractable: it replaces the ODE",
      "right-hand side with a feed-forward neural network whose several hundred trained weights are",
      "not published (only the learned volume V_NODE = 5.7 +/- 0.1 L, Equation 5, and",
      "V_NODE,cov = 5.2 +/- 0.8 L, Equation 6). Analysis code is at",
      "https://github.com/TommyGiak/pharmacoNODE but the trained weights and the patient data are",
      "not distributed. Table 4 also carries a third 'Literature' column transcribed from Cojutti",
      "2024 (Clin Pharmacokinet 63:1271-1282); that column is NOT re-extracted here because it is",
      "already available as modellib('Cojutti_2024_dalbavancin')."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters -- Giacometti 2025 Table 4, "NLME Model"
    # column. The table caption defines exactly what these are: "Mean
    # values of the estimated parameters for the two-compartment model and
    # NLME among the six iterations of the cross-validation. The
    # associated error is the standard deviation." So each value below is
    # the across-fold MEAN, and the quoted spread is the across-fold SD --
    # NOT a standard error, NOT an RSE, and NOT a between-subject
    # variance. Monolix reports the fixed effects on the linear scale;
    # they are log-transformed here per the nlmixr2lib convention.
    #
    # The paper's V1 / V2 map onto the canonical vc / vp; its Cl and Q map
    # onto cl and q. Units are stated in Table 4's row labels: Cl and Q in
    # L/h, V1 and V2 in L.
    # ------------------------------------------------------------------
    lcl <- log(0.0367); label("Clearance (L/h)")                        # Giacometti 2025 Table 4, NLME column: Cl = 0.0367 +/- 0.0006 L/h (across-fold SD)
    lvc <- log(6.32);   label("Central volume of distribution V1 (L)")      # Giacometti 2025 Table 4, NLME column: V1 = 6.32 +/- 0.15 L (across-fold SD)
    lq  <- log(0.028);  label("Intercompartmental clearance (L/h)")      # Giacometti 2025 Table 4, NLME column: Q  = 0.028 +/- 0.003 L/h (across-fold SD)
    lvp <- log(13.9);   label("Peripheral volume of distribution V2 (L)")   # Giacometti 2025 Table 4, NLME column: V2 = 13.9 +/- 0.6 L (across-fold SD)

    # ------------------------------------------------------------------
    # Inter-individual variability: NOT REPORTED.
    #
    # This is a genuine NLME fit, so random effects were certainly
    # estimated -- Methods 2.2 describes the framework ("random effects,
    # which account for inter-individual variability [...] typically
    # assumed to follow a multivariate normal distribution with zero mean
    # and covariance matrix G") and states that Monolix estimated the
    # model by maximum likelihood. But Table 4 is the paper's only
    # parameter table and it prints four rows -- Cl, Q, V1, V2 -- with no
    # omega block, no shrinkage, and no CV%. The Appendix C tables report
    # MSE and R^2 for the Neural ODE arm only.
    #
    # No eta terms are declared, rather than being written as
    # `etalcl ~ fixed(0)`, because a zero-variance diagonal makes OMEGA
    # singular and breaks the Cholesky sampler used by rxSolve. Sibling
    # precedent: inst/modeldb/specificDrugs/Thoueille_2026_salmeterol.R.
    # Simulations from this model are therefore TYPICAL-VALUE ONLY. No
    # variance is invented; see the vignette Errata.
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Residual unexplained variability: NOT REPORTED either.
    #
    # Methods 2.2 gives the objective function as a plain squared loss
    # over observed minus predicted concentrations, and the model was
    # fitted in Monolix by maximum likelihood, so an error model existed;
    # but neither its form (constant / proportional / combined) nor any
    # magnitude appears anywhere in the article. Both terms are therefore
    # encoded as exact zeros rather than invented, following
    # inst/modeldb/specificDrugs/Takada_2025_vancomycin.R. Both a
    # proportional and an additive term are carried so that a downstream
    # user who has a residual-error estimate can supply either without
    # editing model(); the standing convention for an ambiguous residual
    # form is to encode both.
    # ------------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")  # Giacometti 2025: no residual-error estimate is published for the NLME arm
    addSd  <- fixed(0); label("Additive residual SD (mg/L; 0 -- not reported in the source)")          # Giacometti 2025: no residual-error estimate is published for the NLME arm
  })

  model({
    # 1. Individual PK parameters. No covariate and no random effect --
    #    both by the source's own design decision (Appendix B) and by its
    #    reporting gap (no omegas printed), respectively.
    cl <- exp(lcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # 2. Two-compartment IV disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 3. ODE system, Giacometti 2025 Equation 2a. The printed equation is
    #
    #      dCc/dt = -Cc*Cl/V1 - Cc*Q/V1 + Cp*Q/V2
    #      dCp/dt =  Cc*Q/V1  - Cp*Q/V2
    #
    #    which is EXACTLY the canonical mass-conserving two-compartment
    #    system written on AMOUNTS, with the two amount states named as if
    #    they were the concentrations: substituting Cc -> A1 and Cp -> A2
    #    gives dA1/dt = -(CL/V1)A1 - (Q/V1)A1 + (Q/V2)A2 and
    #    dA2/dt = (Q/V1)A1 - (Q/V2)A2, i.e. d(A1 + A2)/dt = -(CL/V1)A1.
    #    Read literally on concentrations it would not conserve mass (the
    #    influx and efflux of the peripheral compartment would be divided
    #    by different volumes), so the labelling in the surrounding prose
    #    -- "Cc and Cp the concentrations of the two compartments" -- is
    #    the slip, not the algebra. This is also the model Monolix was
    #    configured with ("two-compartment distribution, infusion, no
    #    delay and linear elimination", Methods 2.2), which is the standard
    #    mass-conserving form, and it matches the Figure 2 schematic in
    #    which the central compartment is "the only compartment that can
    #    eliminate the antibiotic at rate Cl". Encoded below in amounts,
    #    the nlmixr2lib convention.
    #
    #    Dalbavancin is given intravenously, so doses enter `central`
    #    directly and there is no absorption compartment; Methods 2.2
    #    specifies "infusion, no delay", so the dose is delivered as a
    #    zero-order infusion set in the event table (rate/dur), not by any
    #    model-side term.
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-    k12       * central - k21 * peripheral1

    # 4. Total plasma dalbavancin concentration in mg/L (dose in mg,
    #    volumes in L). The source's assay and its 8.04 mg/L efficacy
    #    threshold are both on the total concentration.
    Cc <- central / vc

    # 5. Observation. Both residual magnitudes are zero (unreported), so
    #    this reduces to a deterministic prediction; the structure is kept
    #    so a user with an error estimate can fill it in via ini().
    Cc ~ prop(propSd) + add(addSd)
  })
}
