Collins_2023_belantamab_mprotein <- function() {
  description <- paste(
    "Concentration-driven tumor growth inhibition (TGI) model for serum",
    "M-protein in patients with relapsed/refractory multiple myeloma treated",
    "with the antibody-drug conjugate belantamab mafodotin, with logistic",
    "growth plus a modified weak Allee term, a saturating effect-compartment",
    "driven kill term, an exponential resistance decay, and covariate effects",
    "of baseline beta-2-microglobulin on the growth rate, baseline M-protein",
    "on the kill rate, and extramedullary disease plus baseline soluble BCMA",
    "on the effect-compartment rate constant (Collins 2023). The embedded",
    "two-compartment time-varying-clearance PK layer is a typical-value ADC",
    "driver taken from Collins 2023 Table S1; for the full covariate",
    "population PK model of belantamab mafodotin see",
    "modellib('Papathanasiou_2025_belantamab'). The paper's ocular-safety",
    "discrete time Markov model is not included (see vignette).",
    sep = " "
  )
  reference <- paste(
    "Collins J, van Noort M, Rathi C, Post TM, Struemper H, Jewell RC,",
    "Ferron-Brady G. Longitudinal efficacy and safety modeling and simulation",
    "framework to aid dose selection of belantamab mafodotin for patients with",
    "multiple myeloma. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(10):1411-1424. doi:10.1002/psp4.13016.",
    "The two-compartment PK structure originates with Struemper et al. (2019)",
    "and is re-estimated in Collins 2023 Table S1; the full covariate",
    "population PK model for the same drug is available as",
    "modellib('Papathanasiou_2025_belantamab').",
    sep = " "
  )
  vignette <- "Collins_2023_belantamab_mprotein"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    MCPROT = list(
      description        = "Baseline serum monoclonal (M) protein concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline value, reported in g/L throughout Collins 2023",
        "(register default is the US-convention g/dL; 1 g/dL = 10 g/L).",
        "Serves two roles in this model: (1) it seeds the initial condition of",
        "the `tumor` state, which holds the serum M-protein concentration",
        "itself, and (2) it is binarized inline as `mcprot_low <- (MCPROT < 20)`",
        "to carry the multiplicative factor 1.41 on the kill rate constant KD",
        "(Collins 2023 Table 1 row 'Effect of bMPROT < 20 on KD' and the",
        "printed equation KDpop = TVKD * theta_bMPROT<20). Patients required a",
        "baseline M-protein >= 5 g/L to enter the M-protein analysis dataset.",
        "Note that this model does NOT use MCPROT as a time-varying covariate",
        "(contrast Ide_2020_elotuzumab, where it modulates target-mediated",
        "elimination) -- here the M-protein time course is the modelled state.",
        sep = " "
      ),
      source_name        = "bMPROT"
    ),
    B2M = list(
      description        = "Baseline serum beta-2-microglobulin concentration",
      units              = "nmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline value, reported in nM (nmol/L) by Collins 2023",
        "(register default is mg/L; ~1 mg/L = 85 nmol/L for the 11.8 kDa",
        "monomer). Power exponent 0.219 on the growth rate constant KGR with",
        "reference 350 nmol/L per the printed equation",
        "KGRpop = TVKGR * (bB2M / 350)^theta_KGR_BB2M (Collins 2023 Table 1).",
        "Interpreted here as a disease-burden marker: Collins 2023 reports high",
        "baseline beta-2-microglobulin as the strongest covariate trend on a",
        "larger growth rate parameter. Simulation-population median 364 nmol/L",
        "(5th-95th percentile 170-1034; Collins 2023 Table S2).",
        sep = " "
      ),
      source_name        = "bB2M"
    ),
    SBCMA = list(
      description        = "Baseline serum soluble B-cell maturation antigen (sBCMA) concentration",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline value in ng/mL, matching the register default",
        "(1 ng/mL = 1 ug/L). Power exponent -0.414 on the effect-compartment",
        "rate constant KEO with reference 100 ng/mL per the printed equation",
        "KEOpop = TVKEO * theta_Extramedullary * (baseline sBCMA / 100)^",
        "theta_KEO_baseline_sBCMA (Collins 2023 Table 1). The 100 ng/mL",
        "reference for the M-protein model is distinct from the 88.2 ng/mL",
        "reference that Collins 2023 Table 1 footnote b uses for the",
        "sBCMA-on-Emax effect in the ocular-events Markov model (out of scope",
        "for this file). Simulation-population median 88.2 ng/mL (5th-95th",
        "percentile 6.28-587; Collins 2023 Table S2).",
        sep = " "
      ),
      source_name        = "baseline sBCMA"
    ),
    DIS_EMD = list(
      description        = "Extramedullary disease at screening",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (medullary / bone-marrow-confined multiple myeloma)",
      notes              = paste(
        "Time-fixed screening flag. Multiplicative factor 0.108 on the",
        "effect-compartment rate constant KEO when DIS_EMD = 1, i.e. a roughly",
        "nine-fold slower plasma-to-effect-site equilibration in patients with",
        "extramedullary disease (Collins 2023 Table 1 row 'Effect of MEDFL on",
        "KEO', 95% CI 0.0617-0.187, excludes 1). Encoded as",
        "e_emd_ke0^DIS_EMD so that DIS_EMD = 0 contributes a factor of 1.",
        "Source column MEDFL maps to canonical DIS_EMD.",
        sep = " "
      ),
      source_name        = "MEDFL"
    )
  )

  covariatesDataExcluded <- list(
    LOT = list(
      description = "Number of prior lines of anti-myeloma therapy",
      units       = "(count)",
      type        = "count",
      notes       = paste(
        "Screened graphically against individual M-protein model parameter",
        "estimates but not retained in the final model; no point estimate is",
        "reported (Collins 2023 Methods, 'Other relevant covariates ... were",
        "assessed with graphical relationships between the individual estimates",
        "of random effects and the covariate values').",
        sep = " "
      )
    ),
    ECOG = list(
      description = "Baseline Eastern Cooperative Oncology Group performance status",
      units       = "(ordinal 0-5)",
      type        = "categorical",
      notes       = paste(
        "Screened graphically but not retained in the final M-protein model;",
        "no point estimate is reported (Collins 2023 Methods).",
        sep = " "
      )
    ),
    MM_NIGG = list(
      description = "Non-IgG-secreting multiple myeloma indicator (immunoglobulin type)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Immunoglobulin type was screened graphically but not retained in the",
        "final M-protein model; no point estimate is reported (Collins 2023",
        "Methods).",
        sep = " "
      )
    ),
    CRP = list(
      description = "Baseline C-reactive protein concentration",
      units       = "mg/L",
      type        = "continuous",
      notes       = paste(
        "Screened graphically but not retained in the final M-protein model;",
        "no point estimate is reported (Collins 2023 Methods).",
        sep = " "
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 169L,
    n_studies      = 2L,
    disease_state  = paste(
      "Relapsed/refractory multiple myeloma (RRMM). Only patients who were",
      "followed for response by serum M-protein and had a baseline serum",
      "M-protein concentration >= 5 g/L entered the M-protein model dataset.",
      sep = " "
    ),
    dose_range     = paste(
      "Belantamab mafodotin by intravenous infusion every 3 weeks:",
      "0.03-4.6 mg/kg in DREAMM-1 (phase I dose escalation) and 2.5 or",
      "3.4 mg/kg in DREAMM-2 (phase II). The approved reference regimen is",
      "2.5 mg/kg every 3 weeks.",
      sep = " "
    ),
    notes          = paste(
      "M-protein model estimation dataset: 169 patients -- DREAMM-1",
      "(NCT02064387; n = 42) and DREAMM-2 (NCT03525678; n = 127) -- per",
      "Collins 2023 Results. The 218-patient DREAMM-2 population summarized",
      "below (and in Collins 2023 Table S2) is the simulation population, not",
      "the estimation population; per-study baseline demographics for the",
      "estimation population are reported only in the primary DREAMM-1 and",
      "DREAMM-2 publications (Trudel 2018; Lonial 2020) and not in Collins",
      "2023 itself. The ocular-events discrete time Markov model used all 218",
      "DREAMM-2 patients and 2273 ocular assessments; it is not part of this",
      "model file.",
      sep = " "
    ),
    simulation_population = paste(
      "All 218 DREAMM-2 patients, retaining their baseline covariate values",
      "from the M-protein / ocular-events model-building set. Mean (SD),",
      "median (5th, 95th percentile) per Collins 2023 Table S2: weight 76.2",
      "(17.3), 73.3 (53.0, 108) kg; albumin 38.3 (4.76), 39 (30, 45) g/L;",
      "beta-2-microglobulin 470 (431), 364 (170, 1034) nmol/L; sBCMA 173 (244),",
      "88.2 (6.28, 587) ng/mL; M-protein 18.5 (15.6), 13.0 (3.0, 50.3) g/L.",
      "Missing baseline serum M-protein values were imputed from the",
      "relationship between serum M-protein and baseline albumin, sBCMA, and",
      "IgG. Each simulation covered 315 days from first dose with 500",
      "replications per regimen.",
      sep = " "
    ),
    biomarker_summary = paste(
      "Baseline (DREAMM-2 simulation population, Collins 2023 Table S2):",
      "beta-2-microglobulin median 364 nmol/L, sBCMA median 88.2 ng/mL,",
      "M-protein median 13.0 g/L.",
      sep = " "
    ),
    dropout        = paste(
      "Patient dropout attributed to progressive disease was included in the",
      "paper's simulations: patients were removed at the time of disease",
      "progression, defined by IMWG criteria for serum M-protein (an increase",
      "of > 25% over the lowest on-treatment M-protein concentration with an",
      "absolute increase of > 0.5 g/dL). Disease progression was the only",
      "discontinuation reason modelled. Dropout is a simulation-harness rule",
      "rather than a model equation and is therefore not encoded in this file.",
      sep = " "
    ),
    software       = "NONMEM 7.3 (M-protein model and simulations); R 3.2.5 or higher."
  )

  ini({
    # ------------------------------------------------------------------
    # Embedded ADC PK driver -- two-compartment model with sigmoidal
    # time-varying clearance, typical values only (Collins 2023 Table S1,
    # "Updated population PK parameters used in the two-compartment
    # model", re-estimated from data available as of December 2019).
    #
    # Covariate effects: Collins 2023 Table S1 reports 11 covariate
    # exponents / factors (BWT, bALB, sex, sBCMA, bIgG, dose < 1 mg/kg,
    # and DREAMM-1 study on CL / V1 / V2), but it does NOT report the
    # reference (centering) value for any of the power terms, so those
    # terms cannot be encoded faithfully. The full covariate population
    # PK model for the same drug -- with published reference values -- is
    # modellib('Papathanasiou_2025_belantamab'). See vignette Errata.
    #
    # Inter-individual variability: Table S1 reports omega^2 CL 0.111,
    # V1 0.0191, Q 0.0916, V2 0.21, IMAX 0.134, TI50 0.168, but its
    # footnote a states that "for the DTMM and M-protein models, post hoc
    # empirical Bayes estimates based on PK data were used to simulate
    # the individual trajectories; random effects variances were not
    # used." The PK layer is therefore encoded as a typical-value driver
    # with no etas, matching how the paper used it. See vignette Errata.
    # ------------------------------------------------------------------
    lcl   <- log(0.915); label("Typical ADC clearance CL at time zero (L/day)")                        # Collins 2023 Table S1: TVCL = 0.915 L/day
    lvc   <- log(4.46);  label("Typical ADC central volume of distribution V1 (L)")                    # Collins 2023 Table S1: TVV1 = 4.46 L
    lq    <- log(0.740); label("Typical ADC intercompartmental clearance Q (L/day)")                   # Collins 2023 Table S1: TVQ = 0.740 L/day
    lvp   <- log(6.49);  label("Typical ADC peripheral volume of distribution V2 (L)")                 # Collins 2023 Table S1: TVV2 = 6.49 L
    cl_hill_max  <- -0.309;     label("Maximal log-fold change in ADC clearance over time IMAX (unitless)")   # Collins 2023 Table S1: IMAX = -0.309
    lcl_hill_t50 <- log(56.9);  label("Log of time at which half the change in ADC clearance has occurred (log days)")  # Collins 2023 Table S1: TI50 = 56.9 days
    cl_hill_gamma <- 3.81;       label("Sigmoidicity (Hill) exponent of time on ADC clearance (unitless)")     # Collins 2023 Table S1: Gamma = 3.81

    # ------------------------------------------------------------------
    # M-protein TGI structural parameters (Collins 2023 Table 1,
    # "M-protein model fixed-effects parameter"; equations printed in
    # Figure 1a).
    # ------------------------------------------------------------------
    lkg      <- log(0.0174);  label("Typical M-protein growth rate constant KGR (1/day)")              # Collins 2023 Table 1: TVKGR = 0.0174 /day (95% CI 0.0146, 0.0208)
    lkd      <- log(0.0207);  label("Typical M-protein kill (death) rate constant KD (1/day)")         # Collins 2023 Table 1: TVKD = 0.0207 /day (95% CI 0.0159, 0.0268)
    llambda  <- log(0.0024);  label("Typical resistance rate constant LAMBDA (1/day)")                 # Collins 2023 Table 1: TVLAMBDA = 0.0024 /day (95% CI 0.00123, 0.00480)
    lslope   <- log(0.619);   label("Slope of the saturating drug-effect term on the kill rate (mL/ug)")  # Collins 2023 Table 1: Slope = 0.619 mL/ug (95% CI 0.554, 0.691)
    lke0     <- log(0.0284);  label("Typical effect-compartment equilibration rate constant KEO (1/day)")  # Collins 2023 Table 1: TVKEO = 0.0284 /day (95% CI 0.0248, 0.0326)

    # A50 is the M-protein concentration at which the growth rate is
    # halved by the modified weak Allee term. Collins 2023 Table 1 labels
    # the row "A50, ug/ml", but the Methods text states "a modified weak
    # Allee effect (50% growth reduction at M-protein < 1.16 g/L)" and
    # M-protein is reported in g/L everywhere else in the paper (baseline
    # mean 18.5 g/L, carrying capacity 150 g/L). The Table 1 unit is a
    # three-orders-of-magnitude typographical error; g/L is used here so
    # that A50 is dimensionally consistent with the `tumor` state. See
    # vignette Errata.
    la50     <- log(1.16);    label("M-protein concentration halving the growth rate A50 (g/L)")       # Collins 2023 Methods text ("50% growth reduction at M-protein < 1.16 g/L"); Table 1 row "A50" value 1.16 (95% CI 0.839, 1.481)

    # Logistic carrying capacity, fixed from the literature rather than
    # estimated: "with the carrying capacity fixed at 150 g/L based on
    # existing literature" (Collins 2023 Methods, citing ref 18).
    ltsmax   <- fixed(log(150)); label("Logistic carrying capacity for serum M-protein Kmax (g/L)")    # Collins 2023 Methods: carrying capacity fixed at 150 g/L; printed in the Figure 1a KGR'(t) equation

    # ------------------------------------------------------------------
    # Covariate effects on the M-protein model (Collins 2023 Table 1 and
    # the three population equations printed beneath it:
    #   KGRpop = TVKGR * (bB2M / 350)^theta_KGR_BB2M
    #   KDpop  = TVKD  * theta_bMPROT<20
    #   KEOpop = TVKEO * theta_Extramedullary
    #            * (baseline sBCMA / 100)^theta_KEO_baseline_sBCMA
    # The two indicator-driven factors (theta_bMPROT<20 and
    # theta_Extramedullary) are printed without an exponent; they are
    # encoded here as factor^indicator so that a patient in the reference
    # category contributes a factor of exactly 1.
    # ------------------------------------------------------------------
    e_b2m_kg     <-  0.219;  label("Power exponent of baseline beta-2-microglobulin on KGR")           # Collins 2023 Table 1: Effect of bB2M on KGR = 0.219 (95% CI 0.011, 0.427)
    e_mcprot_kd  <-  1.41;   label("Multiplicative factor on KD when baseline M-protein < 20 g/L")     # Collins 2023 Table 1: Effect of bMPROT < 20 on KD = 1.41 (95% CI 1.08, 1.85)
    e_emd_ke0    <-  0.108;  label("Multiplicative factor on KEO for extramedullary disease")          # Collins 2023 Table 1: Effect of MEDFL on KEO = 0.108 (95% CI 0.0617, 0.187)
    e_sbcma_ke0  <- -0.414;  label("Power exponent of baseline soluble BCMA on KEO")                   # Collins 2023 Table 1: Effect of baseline sBCMA on KEO = -0.414 (95% CI -0.561, -0.267)

    # ------------------------------------------------------------------
    # Inter-individual variability (Collins 2023 Table 1, "M-protein
    # model random-effects parameter"). KGR, KD, and LAMBDA form a full
    # 3x3 block; the reported variances and covariances are supplied
    # directly as the lower triangle in row-major order:
    #
    #             KGR     KD    LAMBDA
    #   KGR      0.677  0.546   0.217
    #   KD       0.546  1.070   1.120
    #   LAMBDA   0.217  1.120   3.040
    #
    # (correlations 0.642 KGR-KD, 0.151 KGR-LAMBDA, 0.621 KD-LAMBDA; the
    # matrix is positive definite). The KEO, Slope, and A50 variances
    # were held at 0.015 rather than estimated ("FIXED" / "fixed" in the
    # Table 1 rows) and are wrapped in fixed() accordingly.
    # ------------------------------------------------------------------
    etalkg + etalkd + etallambda ~ c(0.677,
                                     0.546, 1.07,
                                     0.217, 1.12, 3.04)  # Collins 2023 Table 1: omega^2 KGR 0.677 (%RSE 15.5, shrinkage 15.8), Covariance(KGR~KD) 0.546 (20.0), omega^2 KD 1.07 (21.1) [19.7], Covariance(KGR~LAMBDA) 0.217 (78.3), Covariance(KD~LAMBDA) 1.12 (46.4), omega^2 LAMBDA 3.04 (32.9) [35.6]

    etalke0   ~ fixed(0.015)  # Collins 2023 Table 1: omega^2 KEO FIXED = 0.015
    etalslope ~ fixed(0.015)  # Collins 2023 Table 1: omega^2 Slope fixed = 0.015
    etala50   ~ fixed(0.015)  # Collins 2023 Table 1: omega^2 A50 fixed = 0.015

    # ------------------------------------------------------------------
    # Residual error.
    #
    # M-protein: Collins 2023 Table 1 reports a combined error model with
    # an additive SD of 0.832 (in g/L, the M-protein unit) and a
    # proportional term of 0.0633 (a fraction, i.e. 6.33%).
    #
    # ADC concentration: Collins 2023 Table S1 reports "RES ERR, additive
    # sigma on log scale" = 0.026. Additive on the natural-log scale is
    # equivalent to proportional in linear space for small errors, so it
    # is encoded as a proportional SD of 0.026. (Note that Table S1
    # reports sigma directly rather than sigma^2; the companion
    # Papathanasiou 2025 model reports the variance 0.0633 for the same
    # log-additive form, i.e. an SD of 0.252, so the two publications
    # differ by roughly an order of magnitude on the ADC residual error.
    # The Collins 2023 value is used here because it is the value on disk
    # for this paper. See vignette Errata.)
    # ------------------------------------------------------------------
    propSd       <- 0.026;  label("ADC proportional residual error (fraction; additive on the log scale)")  # Collins 2023 Table S1: RES ERR, additive sigma on log scale = 0.026
    propSd_tumor <- 0.0633; label("M-protein proportional residual error (fraction)")                  # Collins 2023 Table 1: %RES ERR, proportional = 0.0633 (%RSE 6.2)
    addSd_tumor  <- 0.832;  label("M-protein additive residual error (g/L)")                           # Collins 2023 Table 1: %RES ERR, additive SD = 0.832 (%RSE 3.8)
  })

  model({
    # ------------------------------------------------------------------
    # 1. ADC PK driver -- typical-value two-compartment model with
    #    sigmoidal time-varying clearance (Collins 2023 Table S1; the
    #    functional form CL(t) = CL_0 * exp(IMAX * t^Gamma /
    #    (TI50^Gamma + t^Gamma)) is shared with the companion
    #    Papathanasiou 2025 population PK model of the same drug). The
    #    multiplier is 1 at t = 0 and approaches exp(IMAX) = 0.734 as
    #    t -> infinity, i.e. a 26.6% decrease in clearance over the
    #    treatment course. rxode2's `t` is time since the start of the
    #    simulation, which matches the paper's "Time" since first dose.
    # ------------------------------------------------------------------
    cl_hill_t50 <- exp(lcl_hill_t50)
    cl   <- exp(lcl) * exp(cl_hill_max * t^cl_hill_gamma / (cl_hill_t50^cl_hill_gamma + t^cl_hill_gamma))
    vc   <- exp(lvc)
    q    <- exp(lq)
    vp   <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ------------------------------------------------------------------
    # 2. Individual M-protein model parameters. Covariate factors follow
    #    the three population equations printed beneath Collins 2023
    #    Table 1. The two indicator-driven factors are written as
    #    factor^indicator so the reference category contributes 1.
    # ------------------------------------------------------------------
    mcprot_low <- (MCPROT < 20)   # binary indicator "bMPROT < 20" (g/L) per Collins 2023 Table 1

    kg     <- exp(lkg + etalkg) * (B2M / 350)^e_b2m_kg
    kd     <- exp(lkd + etalkd) * e_mcprot_kd^mcprot_low
    lambda <- exp(llambda + etallambda)
    ke0    <- exp(lke0 + etalke0) * e_emd_ke0^DIS_EMD * (SBCMA / 100)^e_sbcma_ke0
    slope  <- exp(lslope + etalslope)
    a50    <- exp(la50 + etala50)
    tsmax  <- exp(ltsmax)

    # ------------------------------------------------------------------
    # 3. ODE system.
    #
    #    PK: linear two-compartment ADC kinetics. Dose in mg and volumes
    #    in L give central / vc in mg/L = ug/mL, which is the unit of the
    #    Slope parameter's reciprocal (mL/ug).
    #
    #    Effect compartment (Collins 2023 Figure 1a; Sheiner 1979 form):
    #    the effect-site ADC concentration C lags plasma with rate
    #    constant KEO. Collins 2023 abbreviation list defines "C, effect
    #    concentration", and Figure 1a drives the kill term with C.
    #
    #    M-protein TGI (Collins 2023 Figure 1a, verbatim):
    #      dM(t)/dt = KGR'(t) * M(t) - KD'(t) * M(t)
    #      KGR'(t)  = KGR * M(t) / (M(t) + A50) * (1 - M(t) / 150)
    #      KD'(t)   = KD * (1 - exp(-Slope * C)) * exp(-lambda * t)
    #    The M(t) / (M(t) + A50) factor is the modified weak Allee term
    #    (growth slows at very low M-protein) and (1 - M(t) / 150) is the
    #    logistic term enforcing the 150 g/L carrying capacity.
    # ------------------------------------------------------------------
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    Cc <- central / vc

    d/dt(effect) <- ke0 * (Cc - effect)

    kgr_t <- kg * (tumor / (tumor + a50)) * (1 - tumor / tsmax)
    kd_t  <- kd * (1 - exp(-slope * effect)) * exp(-lambda * t)

    d/dt(tumor) <- (kgr_t - kd_t) * tumor

    # ------------------------------------------------------------------
    # 4. Initial conditions. The `tumor` state holds the serum M-protein
    #    concentration in g/L and starts at the subject's baseline value.
    #    The effect compartment starts empty (no drug before the first
    #    dose).
    # ------------------------------------------------------------------
    tumor(0) <- MCPROT

    # ------------------------------------------------------------------
    # 5. Observations and residual error.
    # ------------------------------------------------------------------
    Cc    ~ prop(propSd)
    tumor ~ prop(propSd_tumor) + add(addSd_tumor)
  })
}
