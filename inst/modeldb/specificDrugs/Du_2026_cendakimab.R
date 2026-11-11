Du_2026_cendakimab <- function() {
  description <- paste(
    "Exposure-response model for cendakimab (anti-IL-13 mAb) and dysphagia",
    "days (DD) in adults and adolescents with eosinophilic esophagitis (Du",
    "2026, phase 3 study CC-93538-EE-001). DD, the number of days with",
    "dysphagia in a diary window of NDAYS evaluable days (nominally 14), follows",
    "a combined uniform-binomial (CUB) distribution: with probability pi the",
    "count is Binomial(NDAYS, xi), otherwise it is uniform on 0..NDAYS. The",
    "daily dysphagia probability xi is the inverse logit of a latent disease",
    "variable, logit(xi) = BASE_i - Pmax * (1 - exp(-kpbo * t)) - Emax * (1 - R),",
    "with R a turnover latent variable (R(0) = 1, kin = kout) whose production is",
    "inhibited by cendakimab, Cc / (EC50 + Cc). The between-subject distribution",
    "of BASE is a logit-transformed normal (two shape parameters); Pmax depends",
    "on the individual BASE and Emax on steroid inadequate response or",
    "intolerance. Baseline counts are conditioned on the enrolment criterion",
    "DD >= 4. Cendakimab PK (first-order SC absorption with lag into a",
    "two-compartment linear model) is driven by per-subject empirical-Bayes PK",
    "parameters supplied as data columns, exactly as in the authors' NONMEM",
    "code; the authors' population PK model is not published."
  )
  reference <- paste(
    "Du S, Wojciechowski J, Zhang P, Aras U, Murthy B, Shen J, Kondic A, Hu C.",
    "Latent Variable Indirect Response Modeling of Cendakimab Exposure-Response",
    "for Longitudinal Dysphagia Days Using a Combined Uniform-Binomial",
    "Likelihood Framework. CPT Pharmacometrics Syst Pharmacol. 2026;15(2):e70199.",
    "doi:10.1002/psp4.70199"
  )
  vignette <- "Du_2026_cendakimab"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ug/mL (Cc); dysphagia_days is a count of days"
  )

  covariateData <- list(
    CL_INDIV = list(
      description = "Per-subject empirical-Bayes cendakimab clearance",
      units = "L/h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "NONMEM data column 'CLi' (Data S1 control stream, '$PK CL = CLi').",
        "Individual estimate from the authors' cendakimab population PK model,",
        "which Du 2026 describes as 'not shown' and which has no publication.",
        "For simulation, supply each subject's value; placebo subjects receive",
        "no doses, so any positive placeholder works (the authors used 1). The",
        "vignette uses an illustrative typical PK built by the maintainers from",
        "the simulated concentration profiles of Du 2026 Figure 2."
      ),
      source_name = "CLi"
    ),
    VC_INDIV = list(
      description = "Per-subject empirical-Bayes cendakimab central volume",
      units = "L",
      type = "continuous",
      reference_category = NULL,
      notes = "NONMEM data column 'VCi'. Same upstream source and simulation guidance as CL_INDIV.",
      source_name = "VCi"
    ),
    Q_INDIV = list(
      description = "Per-subject empirical-Bayes cendakimab intercompartmental clearance",
      units = "L/h",
      type = "continuous",
      reference_category = NULL,
      notes = "NONMEM data column 'Qi'. Same upstream source and simulation guidance as CL_INDIV.",
      source_name = "Qi"
    ),
    VP_INDIV = list(
      description = "Per-subject empirical-Bayes cendakimab peripheral volume",
      units = "L",
      type = "continuous",
      reference_category = NULL,
      notes = "NONMEM data column 'VPi'. Same upstream source and simulation guidance as CL_INDIV.",
      source_name = "VPi"
    ),
    KA_INDIV = list(
      description = "Per-subject empirical-Bayes cendakimab first-order SC absorption rate constant",
      units = "1/h",
      type = "continuous",
      reference_category = NULL,
      notes = "NONMEM data column 'KAi'. Same upstream source and simulation guidance as CL_INDIV.",
      source_name = "KAi"
    ),
    TLAG_INDIV = list(
      description = "Per-subject empirical-Bayes cendakimab SC absorption lag time",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = "NONMEM data column 'ALAG1i' (applied as ALAG1 on the depot). Same upstream source and simulation guidance as CL_INDIV.",
      source_name = "ALAG1i"
    ),
    FDEPOT_INDIV = list(
      description = "Per-subject empirical-Bayes cendakimab SC bioavailability",
      units = "(fraction)",
      type = "continuous",
      reference_category = NULL,
      notes = "NONMEM data column 'F1i' (applied as F1 on the depot). Same upstream source and simulation guidance as CL_INDIV.",
      source_name = "F1i"
    ),
    STEROID_IR = list(
      description = "Steroid inadequate responder or intolerant (1) versus steroid responder or naive (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (steroid responder or steroid-naive, Steroid R/N)",
      notes = paste(
        "NONMEM data column 'SIRIN' ('IF (SIRIN.EQ.1) COVSIRINEMAX = SIRINEMAX1',",
        "reference 'no prior history/intolerance'). Pre-specified randomisation",
        "stratification factor of study CC-93538-EE-001, retained in the final",
        "model on clinical grounds (RSE 88.5%). 283 of 427 subjects (66.3%) were",
        "Steroid IR/I (Du 2026 Table 1). Multiplies Emax by exp(0.2377)."
      ),
      source_name = "SIRIN"
    ),
    NDAYS = list(
      description = "Number of evaluable diary days in the dysphagia-days window",
      units = "days",
      type = "count",
      reference_category = NULL,
      notes = paste(
        "NONMEM data column 'NDAYS' ('Period which Dysphagia Days was",
        "evaluated'). The binomial size of the CUB likelihood; nominally 14 but",
        "smaller when diary entries are missing (Du 2026 Methods 'Analysis",
        "Data'). Set to 14 for simulation."
      ),
      source_name = "NDAYS"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "cendakimab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cendakimab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "cendakimab", units = "mg", specimen = "tissue", verified = TRUE),
    latent_dysphagia_days = list(
      analyte = "latent drug-effect variable DE(t)",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 427L,
    n_studies = 1L,
    age_range = "13-72 years",
    age_median = "35 years",
    weight_range = "41.8-176.9 kg",
    weight_median = "81.0 kg",
    sex_female_pct = 30.4,
    race_ethnicity = "White 91.3%, Asian 5.9% (Japanese 5.2%), Black or African American 1.2%, Other 0.7%, missing 0.9%",
    disease_state = "Eosinophilic esophagitis (EoE) with baseline DD >= 4 in a 14-day window; 66.3% steroid inadequate responders or intolerant; 33.7% with a prior history of esophageal dilation.",
    dose_range = "Placebo (n = 143), cendakimab 360 mg SC once weekly for 48 weeks (n = 143), or 360 mg SC once weekly for 24 weeks then every 2 weeks for 24 weeks (n = 141).",
    regions = "Multinational (phase 3 study CC-93538-EE-001).",
    notes = "Du 2026 Table 1 and Table S1: 427 of 430 randomised subjects who received at least one dose and had at least one post-baseline DD assessment; 3 untreated subjects excluded. DD was recorded through Week 48. Adolescents: 34 subjects (9 placebo, 9 QW, 16 QW/Q2W)."
  )

  ini({
    # Typical baseline on the logit scale of the daily dysphagia probability xi.
    logitrbase <- 4.521; label("Typical baseline latent disease variable BASE, logit of the daily dysphagia probability (logit units)") # Table 2: logit-transformed baseline (BASE) = 4.521 (RSE 14.80%); THETA(1)

    # Shape of the logit-transformed IIV on BASE (Du 2026 Equations 9-10).
    logitskew_logitrbase <- -2.419; label("Logit-transformation skewness parameter LGT1 of the IIV on BASE (logit units)") # Table 2: LGT1 = -2.419 (RSE 7.310%); THETA(7)
    logitwidth_logitrbase <- 61.24; label("Logit-transformation width parameter LGT2 of the IIV on BASE (logit units)") # Table 2: LGT2 = 61.24 (RSE 8.550%); THETA(8)

    # Placebo effect: estimated as a half-life in weeks; NONMEM
    # POPPBOKON = LOG(2)/(168*THETA(2)) per hour.
    lkpbo <- log(log(2) / (168 * 27.66)); label("Placebo-onset rate constant kPBOKON (1/h); half-life 27.66 weeks") # Table 2: half-life of PBO onset = 27.66 weeks (RSE 31.60%); THETA(2)
    lpmax <- log(8.092); label("Maximum placebo effect Pmax at the typical baseline (logit units)") # Table 2: Pmax = 8.092 (RSE 7.320%); THETA(3)
    e_logitrbase_pmax <- 0.9115; label("Change in Pmax per logit unit of individual BASE above the typical BASE (unitless)") # Table 2: BASEPmax = 0.9115 (RSE 0.9780%); THETA(10)

    # Drug effect: turnover latent variable with kin = kout = DEKON.
    lkout <- log(log(2) / (168 * 1.408)); label("Drug-effect onset and offset rate constant kDEKON (1/h); half-life 1.408 weeks") # Table 2: half-life of drug onset = 1.408 weeks (RSE 22.60%); THETA(4)
    lemax <- log(1.865); label("Maximum drug effect Emax in steroid responders or steroid-naive subjects (logit units)") # Table 2: Emax = 1.865 (RSE 16.80%); THETA(5)
    e_steroid_ir_emax <- 0.2377; label("Log-scale shift in Emax for steroid inadequate responders or intolerant subjects (unitless)") # Table 2: SIRINEmax = 0.2377 (RSE 88.50%); THETA(11)
    lec50 <- 4.337; label("Cendakimab serum concentration at half-maximal drug effect, EC50 (log ug/mL; 76.5 ug/mL)") # Table 2: log-transformed EC50 = 4.337 (RSE 9.000%); THETA(6)

    # CUB mixture weight on the logit scale; PIP = EXP(PIPT)/(1+EXP(PIPT)).
    logitpinform <- 2.568; label("Logit of the probability pi that a DD record is informative (binomial) rather than uniform (logit units)") # Table 2: uncertainty parameter pi = 2.568 (RSE 3.380%); THETA(9)

    etalogitrbase + etalkpbo ~ c(
      4.309,
      1.178, 4.041
    ) # Table 2: omega2 BASE 4.309, covariance BASE-kPBOKON 1.178, omega2 kPBOKON 4.041
  })

  model({
    # Cendakimab PK from per-subject empirical-Bayes parameters (Data S1 $PK
    # and $DES): first-order SC absorption with lag into a two-compartment
    # linear model.
    cl <- CL_INDIV
    vc <- VC_INDIV
    q <- Q_INDIV
    vp <- VP_INDIV
    ka <- KA_INDIV
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- TLAG_INDIV
    f(depot) <- FDEPOT_INDIV

    # dose mg / volume L = mg/L = ug/mL
    Cc <- central / vc

    # Individual baseline with logit-transformed IIV (Du 2026 Equations 9-10):
    # BASE_i = BASE + (expit(LGT1 + eta) - expit(LGT1)) * LGT2.
    rbase <- logitrbase +
      (expit(logitskew_logitrbase + etalogitrbase) - expit(logitskew_logitrbase)) *
        logitwidth_logitrbase

    kpbo <- exp(lkpbo + etalkpbo)
    pmax <- exp(lpmax) + (rbase - logitrbase) * e_logitrbase_pmax
    kout <- exp(lkout)
    emax <- exp(lemax + e_steroid_ir_emax * STEROID_IR)
    ec50 <- exp(lec50)
    pinform <- expit(logitpinform)

    # Drug-effect latent variable (Du 2026 Equation 5; Data S1
    # DADT(5) = DEKON*(1-CENTC/(EC50+CENTC))-DEKON*A(5), A_0(5) = 1).
    d/dt(latent_dysphagia_days) <- kout * (1 - Cc / (ec50 + Cc)) - kout * latent_dysphagia_days
    latent_dysphagia_days(0) <- 1

    # Placebo effect (Du 2026 Equations 3-4): PBO(t) = exp(-kpbo * t) solves
    # dPBO/dt = -kpbo * PBO with PBO(0) = 1; t is time after the first dose.
    placebo_eff <- pmax * (1 - exp(-kpbo * t))
    drug_eff <- emax * (1 - latent_dysphagia_days)

    # Daily dysphagia probability (Du 2026 Equation 2; Data S1 $ERROR IPRET)
    # and the expected count of the binomial component (IPRE).
    logit_xi <- rbase - placebo_eff - drug_eff
    xi <- expit(logit_xi)
    dd_binom_mean <- NDAYS * xi

    # CUB likelihood of the observed count DV out of NDAYS days (Du 2026
    # Equation 1). Baseline records (t <= 0) are conditioned on the enrolment
    # criterion DD >= 4 (Equations 7-8; Data S1
    # CDFLB = PIP*(1-BINOMIALCDF(3))+((NDAYS-4+1)*(1-PIP))/(NDAYS+1)).
    pdf_dd <- pinform * exp(llikBinom(DV, NDAYS, xi)) + (1 - pinform) / (NDAYS + 1)
    cdf3_dd <- exp(llikBinom(0, NDAYS, xi)) + exp(llikBinom(1, NDAYS, xi)) +
      exp(llikBinom(2, NDAYS, xi)) + exp(llikBinom(3, NDAYS, xi))
    if (t <= 0) {
      lik_dd <- pdf_dd / (pinform * (1 - cdf3_dd) + (NDAYS - 3) * (1 - pinform) / (NDAYS + 1))
    } else {
      lik_dd <- pdf_dd
    }
    ll(dysphagia_days) ~ log(lik_dd)
  })
}
