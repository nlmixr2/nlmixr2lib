Othman_2013_ABT_102 <- function() {
  description <- "Population PK/PD model of body-temperature effects of ABT-102, a TRPV1 antagonist, in 108 healthy adult volunteers across three phase 1 trials (Othman 2013). PK is a one-compartment model with one transit absorption compartment, first-order elimination, and formulation-dependent absorption lag (0.3 h solution, 0.6 h solid dispersion) and relative bioavailability (40% solution vs solid-dispersion reference); PK parameter values are taken from the upstream popPK analysis (Othman 2012, J Clin Pharmacol). The PD layer models body temperature as the additive sum of (a) a measurement-type-dependent baseline (oral thermometer vs core ingestible capsule), (b) a 24-h circadian rhythm (cosine in time with measurement-type-dependent amplitude and a shared 7.6-h phase shift), and (c) an Emax drug effect on plasma concentration with time-driven exponential tolerance (Emax decays with half-life T50 = 28 h). Two parallel outputs (BT_oral, BT_core) are produced with measurement-type-dependent additive residual error; for a given subject only one output is realised (oral thermometer subjects use BT_oral; core ingestible-capsule subjects use BT_core)."

  reference <- paste(
    "Othman AA, Nothaft W, Awni WM, Dutta S.",
    "Effects of the TRPV1 antagonist ABT-102 on body temperature in healthy volunteers:",
    "pharmacokinetic/pharmacodynamic analysis of three phase 1 trials.",
    "Br J Clin Pharmacol. 2013 Apr;75(4):1029-1040.",
    "doi:10.1111/j.1365-2125.2012.04405.x.",
    "PK structure and parameter values adapted from Othman AA, Nothaft W, Awni WM, Dutta S,",
    "Pharmacokinetics of the TRPV1 antagonist ABT-102 in healthy human volunteers:",
    "population analysis of data from 3 phase 1 trials,",
    "J Clin Pharmacol. 2012; 52: 1028-1041",
    "(the upstream popPK reference [25] of the PD paper; PK is fixed from individual",
    "empirical-Bayes estimates of that model in the Othman 2013 PD analysis).",
    sep = " "
  )
  vignette <- "Othman_2013_ABT_102"

  # Etas pair additively with paper-symbol baselines, circadian amplitudes and phase
  # shift rather than with log-transformed parents. The paper (Othman 2013, last
  # paragraph before Table 2) states explicitly: "Inter-subject variability in the
  # baseline and amplitude and phase shift of the circadian rhythm were estimated
  # using additive error models [for parameter x, Px = TVPx + eta_x]". Only EC50 and
  # T50 use exponential / log-normal IIV; those parents are log-transformed
  # (lec50 / lt50) and their etas (etalec50 / etalt50) follow the standard
  # nlmixr2lib convention.
  paper_specific_etas <- c("etabl", "etaamp", "etaps")

  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "ng/mL",
    dosing_notes  = "Oral ABT-102; FORM_SOLUTION = 1 -> solution formulation (Studies 1 and 2), = 0 -> solid-dispersion formulation (Study 3).",
    concentration_notes = "ABT-102 plasma Cc; central is internally in mg with vc in L, Cc = 1000 * central / vc to express in ng/mL. Body temperature outputs BT_oral and BT_core are in degC."
  )

  covariateData <- list(
    FORM_SOLUTION = list(
      description        = "Formulation indicator for ABT-102: 1 = oral solution (Studies 1 and 2), 0 = solid-dispersion (Study 3, the bioavailability and lag-time anchor in Othman 2012 / 2013). The solid-dispersion formulation is the F = 1 reference; the oral solution has a relative bioavailability of 40% and a shorter absorption lag (0.3 h vs 0.6 h) per Othman 2013 PK/PD-model Results paragraph 1 (citing the upstream Othman 2012 popPK fit).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (solid-dispersion formulation; F = 1 anchor and lag = 0.6 h)",
      notes              = "Per-subject categorical covariate fixed by study assignment. In Othman 2013 Studies 1 and 2 the oral-solution formulation was used (single-dose escalation 2-40 mg; multiple-dose 2/4/8 mg twice daily for 7 days) and Study 3 used the solid-dispersion formulation (multiple-dose 1/2/4 mg twice daily for 7 days). The two formulations differ only in absorption-lag time and relative bioavailability; no other formulation-driven differences in PK or PD were retained in the final model (Othman 2013 Discussion paragraph 2: 'No differences between the two formulations in the body temperature effect were distinguishable once the differences in exposure were accounted for'). Encoded as multiplicative log-scale effects on tlag and fdepot in model() below.",
      source_name        = "Formulation (solid-dispersion vs oral-solution)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 108L,
    n_studies       = 3L,
    age_range       = "Adult healthy volunteers (specific range not retabulated in the PD paper)",
    weight_range    = "Adult healthy volunteers (specific range not retabulated in the PD paper)",
    sex_female_pct  = NULL,
    race_ethnicity  = NULL,
    disease_state   = "Healthy adult volunteers (no diagnosed condition; randomized 2:1 to ABT-102:placebo within each dose group in all three studies).",
    dose_range      = paste(
      "Study 1: single dose escalation 2, 6, 18, 30, 40 mg ABT-102 oral solution (45 subjects total, 9 per dose group);",
      "Study 2: multiple twice-daily dose 2, 4, 8 mg ABT-102 oral solution for 7 days (27 subjects, 9 per dose group);",
      "Study 3: multiple twice-daily dose 1, 2, 4 mg ABT-102 solid-dispersion formulation for 7 days (36 subjects, 12 per dose group).",
      sep = " "
    ),
    regions         = "Not specified",
    notes           = paste(
      "108 subjects total contributing 7493 body-temperature measurements (2696 oral thermometer in Studies 1 and 2; 4797 core ingestible-capsule in Study 3, after exclusion of 51 erroneous core values below 34 degC that coincided with ingestion of cold liquids and were treated as measurement artifacts).",
      "Demographics and disposition of the 108 subjects were not retabulated in the body-temperature PK/PD paper; readers are referred to Othman 2012 (J Clin Pharmacol 52:1028-1041) for the demographic detail.",
      "Oral body temperatures ranged from 34.3 to 38.4 degC (mean 36.5, median 36.4); core body temperatures ranged from 34.0 to 38.7 degC (median 37.2).",
      sep = " "
    ),
    nonmem_method   = "FOCE with interaction (NONMEM VI; Icon Development Solutions, Ellicott City, MD); ADVAN6 user-defined subroutine.",
    pd_observations_oral = 2696L,
    pd_observations_core = 4797L
  )

  ini({
    # =====================================================================
    # PK component -- fixed-structure values reproduced inline from the
    # upstream popPK analysis (Othman 2012, J Clin Pharmacol 52:1028-1041,
    # the reference [25] of the PD paper). The Othman 2013 PD analysis used
    # individual empirical-Bayes PK estimates from that fit as fixed inputs;
    # the parameter values below are the population estimates re-reported in
    # the PD paper's PK/PD-model Results paragraph 1 (page 1033 in the
    # journal printing; bootstrap 95% CIs in parentheses are from the same
    # paragraph). One-compartment disposition with a single transit
    # compartment for absorption, first-order elimination, and formulation-
    # dependent lag and relative bioavailability.
    # =====================================================================

    lcl  <- log(16)        ; label("Apparent clearance CL/F (L/h)")                       # Othman 2013 PK/PD-model Results paragraph 1 (CL/F 16 L/h, 95% CI 14-18); upstream Othman 2012 popPK
    lvc  <- log(215)       ; label("Apparent central volume V/F (L)")                      # Othman 2013 PK/PD-model Results paragraph 1 (V/F 215 L, 95% CI 192-237); upstream Othman 2012 popPK
    lktr <- log(1.4)       ; label("Transit-absorption rate constant ktr (1/h)")           # Othman 2013 PK/PD-model Results paragraph 1 (ktr 1.4 1/h, 95% CI 1.3-1.6); upstream Othman 2012 popPK
    ltlag <- log(0.6)      ; label("Absorption lag time for solid-dispersion (h)")         # Othman 2013 PK/PD-model Results paragraph 1 (solid-dispersion lag 0.6 h, 95% CI 0.5-0.8); upstream Othman 2012 popPK

    # Formulation effect on absorption lag (log fold-change for oral solution vs solid-dispersion reference).
    e_form_solution_tlag <- log(0.3 / 0.6) ; label("Log fold-change in absorption lag time for solution vs solid-dispersion")  # Othman 2013 PK/PD-model Results paragraph 1 (solution lag 0.3 h vs solid-dispersion lag 0.6 h); upstream Othman 2012 popPK

    # Bioavailability anchor: solid-dispersion = 1 (reference); oral-solution = 0.40.
    lfdepot <- fixed(log(1)) ; label("Bioavailability anchor for solid-dispersion reference (log scale)")  # Othman 2013 PK/PD-model Results paragraph 1 (solid-dispersion is the bioavailability reference); upstream Othman 2012 popPK
    e_form_solution_fdepot <- log(0.40)    ; label("Log relative bioavailability of solution vs solid-dispersion") # Othman 2013 PK/PD-model Results paragraph 1 (solution Frel 40%, 95% CI 35-45%); upstream Othman 2012 popPK

    # PK IIV (additive on log-scale; CV% converted via omega^2 = log(CV^2 + 1)).
    etalcl   ~ 0.1033      # CV = 33% (Othman 2013 PK/PD-model Results paragraph 1); omega^2 = log(0.33^2 + 1) = 0.1033. Upstream Othman 2012 popPK
    etalvc   ~ 0.0560      # CV = 24% (Othman 2013 PK/PD-model Results paragraph 1); omega^2 = log(0.24^2 + 1) = 0.0560. Upstream Othman 2012 popPK
    etalktr  ~ 0.1769      # CV = 44% (Othman 2013 PK/PD-model Results paragraph 1); omega^2 = log(0.44^2 + 1) = 0.1769. Upstream Othman 2012 popPK
    etaltlag ~ 0.2152      # CV = 49% (Othman 2013 PK/PD-model Results paragraph 1); omega^2 = log(0.49^2 + 1) = 0.2152. Upstream Othman 2012 popPK

    # =====================================================================
    # PD component -- body-temperature parameters (Othman 2013 Table 2 final
    # estimates; 95% bootstrap CIs from the same table, computed from 1000
    # successfully converging bootstrap replicates).
    # =====================================================================

    # Baselines (additive IIV on linear scale; Othman 2013 reports a common
    # inter-subject variance shared across oral and core; see model-building
    # step Model 4 in Othman 2013 Table 1 "different baselines for core vs.
    # oral body temperature; same ISV variance").
    bl_oral <- 36.3        ; label("Baseline body temperature, oral thermometer (degC)")              # Othman 2013 Table 2 (baseline oral 36.3 degC, RSE 0.1%, bootstrap 95% CI 36.3-36.4)
    bl_core <- 37.0        ; label("Baseline body temperature, core ingestible capsule (degC)")        # Othman 2013 Results paragraph 4 + Abstract (baseline core 37.0 degC, 95% CI 37.0-37.1); estimated jointly with bl_oral via measurement-type-dependent fixed effect (Model 4 in Table 1)

    # Circadian-rhythm amplitudes (additive IIV; Othman 2013 reports
    # measurement-type-dependent amplitudes with a common ISV variance --
    # Model-building step Final model in Othman 2013 Table 1 "different
    # circadian rhythm amplitudes for oral vs. core temperature measurements;
    # same ISV variance").
    amp_oral <- 0.25       ; label("Circadian-rhythm amplitude, oral (degC)")                          # Othman 2013 Table 2 (amplitude oral 0.25 degC, RSE 6%, bootstrap 95% CI 0.22-0.28)
    amp_core <- 0.31       ; label("Circadian-rhythm amplitude, core (degC)")                          # Othman 2013 Table 2 (amplitude core 0.31 degC, RSE 5%, bootstrap 95% CI 0.28-0.34)

    # Phase shift (shared across measurement types; additive IIV).
    ps <- 7.6              ; label("Circadian-rhythm phase shift in the cosine argument (h)")           # Othman 2013 Table 2 (phase-shift 7.6 h, RSE 2%, bootstrap 95% CI 7.3-7.9). Dosing was performed at approximately 08:00; the cosine peak at t = ps = 7.6 h corresponds to ~15:36 (3:36 PM), matching the typical late-afternoon body-temperature maximum (Othman 2013 Discussion paragraph 4).

    # ABT-102 drug-effect parameters (Emax with time-driven exponential
    # tolerance; subject-independent Emax per Othman 2013 model equation).
    emax  <- 2.2           ; label("ABT-102 maximum body-temperature increase at time = 0 (degC)")     # Othman 2013 Table 2 (Emax 2.2 degC, RSE 8%, bootstrap 95% CI 1.9-2.7)
    lec50 <- log(20)       ; label("ABT-102 plasma concentration at half-maximum effect, EC50 (ng/mL)") # Othman 2013 Table 2 (EC50 20 ng/mL, RSE 15%, bootstrap 95% CI 15-28)
    lt50  <- log(28)       ; label("Tolerance half-time T50 for Emax exponential decay (h)")            # Othman 2013 Table 2 (T50 28 h, RSE 17%, bootstrap 95% CI 20-43)

    # PD IIV.
    # Baseline IIV is shared between oral and core ("same ISV variance" -- Model 4 in Table 1).
    etabl   ~ 0.05          # Othman 2013 Table 2 row omega^2_Baseline (0.05 degC^2, RSE 15%, bootstrap 95% CI 0.03-0.06); additive eta on bl shared by bl_oral and bl_core. SD ~ 0.22 degC matches Discussion paragraph 4 (standard deviation for baseline approximately 0.2 degC).
    # Amplitude IIV is shared between oral and core (same ISV variance -- Final model in Table 1).
    etaamp  ~ 0.008         # Othman 2013 Table 2 row omega^2_Amplitude (0.008 degC^2, RSE 21%, bootstrap 95% CI 0.005-0.011); additive eta on amplitude. SD ~ 0.089 degC matches Discussion paragraph 4 (approximately 0.1 degC).
    etaps   ~ 0.69          # Othman 2013 Table 2 row omega^2_Phase-shift (0.69 h^2, RSE 31%, bootstrap 95% CI 0.32-1.35); additive eta on phase shift. SD ~ 0.83 h.
    etalec50 ~ 1.57         # Othman 2013 Table 2 row omega^2_EC50 (1.57 log-scale variance, RSE 18%, bootstrap 95% CI 1.05-2.15); exponential IIV log-normal. CV approximately sqrt(exp(1.57) - 1) * 100 = 200%.
    etalt50  ~ 0.34         # Othman 2013 Table 2 row omega^2_Tolerance-T50 (0.34 log-scale variance, RSE 28%, bootstrap 95% CI 0.18-0.53); exponential IIV. CV approximately sqrt(exp(0.34) - 1) * 100 = 64%.

    # Residual error (per measurement type).
    # Othman 2013 used a combined additive plus inverse-proportional residual error model
    # of the form Res_kij = (39.5 - BT_ij) * kappa_k_2 * eps_2 + kappa_k_1 * eps_1
    # where the additive component (eps_1 with variance sigma^2_1) is constant and the
    # inverse-proportional component (eps_2 with variance sigma^2_2) is multiplied by
    # (39.5 - BT) so that residual variability is larger at low body temperatures, a
    # device-artifact pattern motivated by erroneous low readings from incomplete
    # equilibration (oral thermometer) or ingested cold liquids (core capsule)
    # (Othman 2013 Discussion paragraph 5 + Results equations). The Othman 2013 Table 2
    # estimates are:
    #   Oral: sigma^2_1 = 0.22 (additive variance, degC^2),   sigma^2_2 = 0.04 (inverse-proportional variance)
    #   Core: sigma^2_1 = 0.33 (additive variance, degC^2),   sigma^2_2 = 0.02 (inverse-proportional variance)
    # nlmixr2's residual-error syntax does not directly express the (39.5 - BT) state-
    # dependent coefficient on a residual variance. The encoding below uses a single
    # additive effective SD per measurement type, evaluated at the cohort-median body
    # temperature of 37 degC:
    #   addSd_BT_oral = sqrt(sigma^2_1 + (39.5 - 37)^2 * sigma^2_2) = sqrt(0.22 + 6.25 * 0.04) = sqrt(0.47) ~= 0.686 degC
    #   addSd_BT_core = sqrt(sigma^2_1 + (39.5 - 37)^2 * sigma^2_2) = sqrt(0.33 + 6.25 * 0.02) = sqrt(0.455) ~= 0.675 degC
    # See vignette Errata for the rationale and the implication that residuals at very
    # low BT (around 34 degC) will be under-dispersed and residuals at very high BT
    # (around 39 degC) will be over-dispersed relative to the paper's published error
    # model. For most simulation use cases (BT near the cohort median), the difference
    # is negligible (< 0.1 degC RMS).
    addSd_BT_oral <- 0.686    ; label("Effective additive residual SD for oral body-temperature observations (degC)")  # Othman 2013 Table 2 oral sigma^2_1 = 0.22, sigma^2_2 = 0.04; effective at BT = 37 degC
    addSd_BT_core <- 0.675    ; label("Effective additive residual SD for core body-temperature observations (degC)")  # Othman 2013 Table 2 core sigma^2_1 = 0.33, sigma^2_2 = 0.02; effective at BT = 37 degC
  })

  model({
    # ---- PK individual parameters (log-normal IIV on the log-scale parents).
    cl     <- exp(lcl   + etalcl)
    vc     <- exp(lvc   + etalvc)
    ktr    <- exp(lktr  + etalktr)
    tlag   <- exp(ltlag + etaltlag + e_form_solution_tlag * FORM_SOLUTION)
    fdepot <- exp(lfdepot           + e_form_solution_fdepot * FORM_SOLUTION)
    kel    <- cl / vc

    # ---- PK ODEs: one-compartment disposition with a single transit
    # compartment between depot and central, both governed by the transit
    # rate constant ktr (Savic chain with n = 1 transit).
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot - ktr * transit1
    d/dt(central)  <-  ktr * transit1 - kel * central

    # ---- Bioavailability and absorption lag (per FORM_SOLUTION).
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # ---- Plasma ABT-102 concentration. central is in mg and vc in L so
    # central / vc is in mg/L; multiply by 1000 to express Cc in ng/mL (the
    # paper's reporting unit and the unit of EC50 below).
    Cc <- 1000 * central / vc

    # ---- PD individual parameters.
    bl_oral_i  <- bl_oral  + etabl     # additive IIV; shared variance across measurement types
    bl_core_i  <- bl_core  + etabl
    amp_oral_i <- amp_oral + etaamp    # additive IIV; shared variance across measurement types
    amp_core_i <- amp_core + etaamp
    ps_i       <- ps       + etaps
    ec50_i     <- exp(lec50 + etalec50)
    t50_i      <- exp(lt50  + etalt50)

    # ---- Circadian-rhythm contribution. Cosine with a 24-h period; t carries
    # the rxode2-default time-since-study-start convention (t = 0 corresponds
    # to the first dose at approximately 08:00 in the source studies).
    cr_oral <- amp_oral_i * cos(2 * pi / 24 * (t - ps_i))
    cr_core <- amp_core_i * cos(2 * pi / 24 * (t - ps_i))

    # ---- ABT-102 drug-effect contribution: Emax of plasma concentration with
    # time-driven exponential tolerance (Emax decays with half-life T50 from
    # an initial value at t = 0). Subject-independent Emax (no eta).
    de <- (emax * exp(-log(2) * t / t50_i) * Cc) / (ec50_i + Cc)

    # ---- Body-temperature outputs (two parallel outputs; for a given
    # subject only one is realised, selected by the measurement device).
    BT_oral <- bl_oral_i + cr_oral + de
    BT_core <- bl_core_i + cr_core + de

    # ---- Residual error.
    BT_oral ~ add(addSd_BT_oral)
    BT_core ~ add(addSd_BT_core)
  })
}
