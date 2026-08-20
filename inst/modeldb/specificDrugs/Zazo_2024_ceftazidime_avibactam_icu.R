Zazo_2024_ceftazidime_avibactam_icu <- function() {
  description <- "Semi-mechanistic PK/PD model of ceftazidime-avibactam against Pseudomonas aeruginosa (strain 2154) in the CRITICALLY ILL (intensive care unit) population. One-compartment linear PK for each drug, with volume proportional to body weight and clearance a linear function of creatinine clearance, driving a two-state bacterial growth model (an actively growing population carried on the log10 scale and a resting population carried on the linear CFU/mL scale) under a shared logistic capacity limit. Ceftazidime kill follows a sigmoidal Emax function whose EC50 is lowered biexponentially by avibactam; exponential delay functions retard growth and initial kill; ceftazidime is degraded by a bacterial-density-dependent Hill function that avibactam inhibits. Sibling model for the non-ICU control population: Zazo_2024_ceftazidime_avibactam_control."
  reference <- "Zazo H, Aguazul Y, Lanao JM. Dosing Evaluation of Ceftazidime-Avibactam in Intensive Care Unit Patients Based on Pharmacokinetic/Pharmacodynamic (PK/PD) Modeling and Simulation. Antibiotics (Basel). 2024;13(9):861. doi:10.3390/antibiotics13090861. Structural equations: Section 4.2, Eqs 3-9. PK parameters: Table 5, 'ICU Patient Population' columns (derived by the authors from refs 17-19: Stein 2019, Merdjan 2017, Welage 1984). Bacterial growth/kill parameters: Table 6 (adapted from Sy SKB et al. J Antimicrob Chemother 2018;73:1295-1304)."
  vignette <- "Zazo_2024_ceftazidime_avibactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Ceftazidime and avibactam are co-administered drugs; neither is a
  # metabolite of the other, so neither maps onto the registered
  # `<canonical>_<metab>` pattern. The two bacterial states are
  # paper-mechanistic. Same declaration style as the sibling
  # ceftazidime/avibactam extraction Kroemer_2024_..._hfim (conc_caz /
  # conc_avi / bact_*) and Mohamed_2012_colistin (central_col / bact_s).
  paper_specific_compartments <- c(
    "central_caz", "central_avi", "bact_active", "bact_resting"
  )

  # NOTE ON THE STATE SCALES. Eq 5 is written as d(log10(P1))/dt, so
  # `bact_active` holds log10(P1) in log10 CFU/mL. Eq 7 is written as
  # d(P2)/dt with linear (P2) and (P1), so `bact_resting` holds P2 in
  # CFU/mL. That asymmetry is what the paper prints, and it is the
  # reading that reproduces the Table 1 bacterial densities (see the
  # vignette "Assumptions and deviations" section for the arithmetic).
  compartmentData <- list(
    central_caz  = list(analyte = "ceftazidime", units = "mg", specimen = "plasma", verified = TRUE),
    central_avi  = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE),
    bact_active  = list(analyte = "Pseudomonas aeruginosa, actively growing population P1 (state is log10 CFU/mL)", units = "log10 CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_resting = list(analyte = "Pseudomonas aeruginosa, resting population P2 (state is linear CFU/mL)", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight. Enters Eq 3 as a direct linear multiplier of the distribution-volume coefficient Vd (L/kg): V (L) = Vd * Weight (kg). No allometric exponent and no reference weight are used.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline (time-fixed). The Monte Carlo simulations sampled body weight from a log-normal distribution with mean 75 kg and CV 13% (Section 4.1). Vd * 75 kg reproduces the Table 2 volumes exactly: 0.46 * 75 = 34.5 L vs 34.5 L (ceftazidime) and 0.68 * 75 = 51.0 L vs 51.0 L (avibactam).",
      source_name = "Weight"
    ),
    CRCL = list(
      description = "Creatinine clearance. Enters Eq 4 as the linear predictor of total serum clearance for each drug: Cl (L/h) = Cli + CLs * Clcr (mL/min).",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "RAW creatinine clearance in mL/min, NOT normalised to 1.73 m^2 body surface area -- the Eq 4 slope CLs is reported per mL/min and the six simulated subpopulations are defined by raw Clcr (100, 60, 40, 20, 10 and 3 mL/min). Same raw Cockcroft-Gault-style usage as Delattre_2010_amikacin. The paper does not state which creatinine-clearance formula was used. Sampled log-normally about each subpopulation mean with CV 20% (Section 4.1).",
      source_name = "Clcr"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1000L,
    n_studies = 0L,
    disease_state = "Simulated critically ill (intensive care unit) adults infected with Pseudomonas aeruginosa (strain 2154; ATCC 15442 is the ICU-acquired reference strain named in the Introduction). Six renal-function subpopulations defined by creatinine clearance: 100, 60, 40, 20, 10 and 3 mL/min.",
    weight_median = "75 kg (log-normal, CV 13%)",
    renal_function = "Clcr 100, 60, 40, 20, 10 and 3 mL/min (log-normal about each mean, CV 20%)",
    dose_range = "Summary-of-product-characteristics regimens, 2 h IV infusion: (2.0 g CAZ / 0.5 g AVI) q8h for Clcr 100 and 60; (1.0/0.25) q8h for Clcr 40; (0.75/0.19) q12h for Clcr 20; (0.75/0.19) q24h for Clcr 10; (0.75/0.19) q48h for Clcr 3 (Table 4). Shortened-interval regimens proposed for Clcr <= 10 mL/min: (0.75/0.19) q12h and q24h.",
    notes = paste(
      "This is a simulation study, not a fit to individual patient data: 1000 Monte Carlo runs per subpopulation over the first week (168 h) of treatment, executed in GoldSim Pro v10.5.",
      "The ICU PK parameters (Table 5) were derived by the authors from published studies in critically ill subjects (refs 17-19); the sibling model Zazo_2024_ceftazidime_avibactam_control carries the non-ICU parameter set. Relative to the control population the ICU volumes are roughly doubled -- the paper attributes this to the systemic inflammatory response syndrome (Discussion) -- while clearance depends less steeply on creatinine clearance.",
      "All PK/PD parameters were sampled from log-normal distributions using the tabulated geometric mean and geometric standard deviation (Section 4.2, 'Monte Carlo simulations').",
      "Model validation was a numerical predictive check against observed data from Falcone 2020 (ref 7), Fresan 2022 (ref 11) and Sy 2018 (ref 12), summarised as average fold errors of 1.01-1.40 in Table 1.",
      "PK/PD targets: T > MIC > 90% and Cmin/MIC > 1.3, with MIC50 for meropenem-non-susceptible P. aeruginosa taken as 4 mg/L (ceftazidime) and 1 mg/L (avibactam)."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Pharmacokinetics -- Table 5, "Control Population" columns.
    # Eq 3: V (L)    = Vd (L/kg) * Weight (kg)
    # Eq 4: Cl (L/h) = Cli (L/h) + CLs * Clcr (mL/min)
    # Log-transformed here because the Monte Carlo sampled every
    # parameter from a log-normal distribution (Section 4.2).
    # ---------------------------------------------------------------
    lvd_caz  <- log(0.46); label("Ceftazidime distribution-volume coefficient Vd (L/kg)")           # Table 5, ICU CAZ: 0.46 (30.4)
    lvd_avi  <- log(0.68); label("Avibactam distribution-volume coefficient Vd (L/kg)")             # Table 5, ICU AVI: 0.68 (27.9)
    lcli_caz <- log(1.15); label("Ceftazidime clearance intercept Cli (L/h)")                       # Table 5, ICU CAZ: 1.15 (54.8)
    lcli_avi <- log(0.89); label("Avibactam clearance intercept Cli (L/h)")                         # Table 5, ICU AVI: 0.89 (65.2)
    lcls_caz <- log(0.04); label("Ceftazidime clearance slope CLs ((L/h) per mL/min of Clcr)")      # Table 5, ICU CAZ: 0.04 (37.2)
    lcls_avi <- log(0.10); label("Avibactam clearance slope CLs ((L/h) per mL/min of Clcr)")        # Table 5, ICU AVI: 0.10 (30.0)

    # ---------------------------------------------------------------
    # Bacterial growth and kill -- Table 6 (adapted from Sy et al. 2018).
    #
    # EXPONENT-SIGN ERRATUM. Table 6 of the published PDF prints
    # "9.67 E2", "4.23 E2", "5.0 E3" and "7.71 E2" with a superscript 2
    # or 3 and NO minus sign. All four exponents must be NEGATIVE:
    #   * with the positive exponents the model predicts a bacterial
    #     density of ~0 log10 CFU/mL for all three static time-kill
    #     conditions against the 9.23 / 8.78 / 6.57 printed in Table 1
    #     (sum of squared errors 205 vs 0.04 for the negative reading);
    #   * k1-2 is added to kdeath,1 in Eq 5 and so must carry units of
    #     1/h -- 5000/h would eradicate the inoculum instantly;
    #   * Degmax = 771/h would give ceftazidime a 3-second half-life.
    # See the vignette "Assumptions and deviations" section.
    # ---------------------------------------------------------------
    llog10nmax <- log(9.89);    label("Maximum bacterial load capacity, log10 Nmax (log10 CFU/mL)")  # Table 6 Nmax: 9.89 [2.06]; see units note below
    lkgrowth1  <- log(0.346);   label("First-order growth rate constant of the active population P1 (1/h)")  # Table 6 kgrowth,1: 0.346 [20.9]
    lemax      <- log(0.240);   label("Maximum ceftazidime kill-rate constant (1/h)")                 # Table 6 Emax: 0.240 [16.0]
    lec50_a    <- log(52.3);    label("First coefficient A of the biexponential ceftazidime EC50 function (mg/L)")   # Table 6 A: 52.3 [17.2]
    lec50_b    <- log(12.6);    label("Second coefficient B of the biexponential ceftazidime EC50 function (mg/L)")  # Table 6 B: 12.6 [26.0]
    lalpha     <- log(2.38);    label("Exponential constant on A relating avibactam concentration to ceftazidime potency (L/mg)")  # Table 6 alpha: 2.38 [119]
    lbeta      <- log(9.67e-2); label("Exponential constant on B relating avibactam concentration to ceftazidime potency (L/mg)")  # Table 6 beta: printed "9.67 E2"; negative exponent -- see erratum note above
    lhill      <- log(2.60);    label("Hill coefficient of the sigmoidal ceftazidime kill function (unitless)")       # Table 6 gamma: 2.60 [34.23]

    # Delay functions, Eq 8. delta1 is reported without a CV.
    ldelta1 <- fixed(log(4.23e-2)); label("Exponential delay constant retarding growth of the active population (1/h)")  # Table 6 delta1: printed "4.23 E2 (fixed)"; negative exponent -- see erratum note above
    ldelta2 <- log(0.213);          label("Exponential delay constant slowing the initial kill of the active population (1/h)")  # Table 6 delta2: 0.213 [17.46]

    # Active -> resting conversion, Eqs 5 and 7. Reported without a CV.
    lkconv <- fixed(log(5.0e-3)); label("Rate constant for conversion of bacteria from the active to the resting state (1/h)")  # Table 6 k1-2: printed "5.0 E3 (fixed)"; negative exponent -- see erratum note above

    # Bacterial degradation of ceftazidime, Eq 9.
    ldegmax   <- log(7.71e-2); label("Maximum bacterial degradation rate constant of ceftazidime (1/h)")  # Table 6 Degmax: printed "7.71 E2 [51.9]"; negative exponent -- see erratum note above
    lkm_deg   <- fixed(log(8.5)); label("Bacterial density giving 50% of the maximum ceftazidime degradation rate (log10 CFU/mL)")  # Table 6 Km: 8.5 (fixed)
    lhill_deg <- log(1.46);    label("Hill coefficient of the ceftazidime-degradation function (unitless)")  # Table 6 phi: 1.46 [82.2]
    lic50_deg <- log(1.96);    label("Avibactam concentration halving the ceftazidime degradation rate (mg/L)")  # Table 6 IC50: 1.96 [58.2]

    # ---------------------------------------------------------------
    # Variability. The paper reports a CV(%) beside every PK and PD
    # parameter and states that "the parameters were sampled from a
    # log-normal distribution" for each of the 1000 virtual patients
    # (Section 4.2). Each CV is therefore carried here as an
    # exponential inter-individual variability, with
    #   omega^2 = log(CV^2 + 1).
    # Note the provenance difference the paper draws in Section 4.2:
    # the Table 5 CVs on Cli and CLs are ESTIMATION ERRORS (parameter
    # uncertainty) while the CV on Vd is inter-individual variability.
    # The Monte Carlo drew both the same way, so both are encoded as
    # etas; see the vignette "Assumptions and deviations" section.
    # ---------------------------------------------------------------
    etalvd_caz  ~ 0.08839176  # Table 5, ICU CAZ Vd  CV 30.4%; log(0.304^2 + 1)
    etalvd_avi  ~ 0.07495997  # Table 5, ICU AVI Vd  CV 27.9%; log(0.279^2 + 1)
    etalcli_caz ~ 0.26259808  # Table 5, ICU CAZ Cli CV 54.8%; log(0.548^2 + 1)
    etalcli_avi ~ 0.35424479  # Table 5, ICU AVI Cli CV 65.2%; log(0.652^2 + 1)
    etalcls_caz ~ 0.12960971  # Table 5, ICU CAZ CLs CV 37.2%; log(0.372^2 + 1)
    etalcls_avi ~ 0.08617770  # Table 5, ICU AVI CLs CV 30.0%; log(0.300^2 + 1)

    etallog10nmax ~ 0.00042427  # Table 6 Nmax     CV  2.06%; log(0.0206^2 + 1)
    etalkgrowth1  ~ 0.04275389  # Table 6 kgrowth1 CV 20.9%;  log(0.209^2 + 1)
    etalemax      ~ 0.02527781  # Table 6 Emax     CV 16.0%;  log(0.160^2 + 1)
    etalec50_a    ~ 0.02915484  # Table 6 A        CV 17.2%;  log(0.172^2 + 1)
    etalec50_b    ~ 0.06541314  # Table 6 B        CV 26.0%;  log(0.260^2 + 1)
    etalalpha     ~ 0.88215467  # Table 6 alpha    CV 119%;   log(1.19^2 + 1)
    etalbeta      ~ 0.00490198  # Table 6 beta     CV  7.01%; log(0.0701^2 + 1)
    etalhill      ~ 0.11079807  # Table 6 gamma    CV 34.23%; log(0.3423^2 + 1)
    etaldelta2    ~ 0.03002972  # Table 6 delta2   CV 17.46%; log(0.1746^2 + 1)
    etaldegmax    ~ 0.23851362  # Table 6 Degmax   CV 51.9%;  log(0.519^2 + 1)
    etalhill_deg  ~ 0.51622144  # Table 6 phi      CV 82.2%;  log(0.822^2 + 1)
    etalic50_deg  ~ 0.29171692  # Table 6 IC50     CV 58.2%;  log(0.582^2 + 1)

    # The paper reports no residual-error model: it is a forward Monte
    # Carlo simulation in GoldSim, not an estimation run, so there is
    # no $SIGMA to transcribe. Fixed at zero rather than invented.
    addSd_Cc_caz    <- fixed(0); label("Additive residual SD on the ceftazidime concentration (mg/L; not reported by the source)")
    addSd_Cc_avi    <- fixed(0); label("Additive residual SD on the avibactam concentration (mg/L; not reported by the source)")
    addSd_bactTotal <- fixed(0); label("Additive residual SD on the total bacterial density (log10 CFU/mL; not reported by the source)")
  })

  model({
    # ---- 1. Individual PK parameters (Eqs 3 and 4) ----
    vc_caz <- exp(lvd_caz + etalvd_caz) * WT
    vc_avi <- exp(lvd_avi + etalvd_avi) * WT
    # Eq 4's intercept and slope each get their own simple line so that both
    # thetas sit in a mu-referenced position; combining them inline as
    # `exp(a + eta_a) + exp(b + eta_b) * CRCL` puts the slope theta outside a
    # mu-ref position and rxode2 warns that the etas defaulted to non-mu.
    cli_caz <- exp(lcli_caz + etalcli_caz)
    cls_caz <- exp(lcls_caz + etalcls_caz)
    cli_avi <- exp(lcli_avi + etalcli_avi)
    cls_avi <- exp(lcls_avi + etalcls_avi)
    cl_caz <- cli_caz + cls_caz * CRCL
    cl_avi <- cli_avi + cls_avi * CRCL

    # ---- 2. Individual PD parameters (Table 6) ----
    # log10nmax is log10(Nmax), i.e. 9.89 log10 CFU/mL. Table 6 labels
    # the Nmax units "CFU/mL", but Eq 5 divides by log10(Nmax) and only
    # the log10 reading reproduces Table 1: taking Nmax = 9.89 CFU/mL
    # literally gives log10(Nmax) = 0.995 and a plateau near
    # 0.7 log10 CFU/mL instead of the published 6.57-9.23.
    log10nmax <- exp(llog10nmax + etallog10nmax)
    kgrowth1  <- exp(lkgrowth1 + etalkgrowth1)
    # Table 6: kgrowth,2 = (1/10^7) * kgrowth,1 -- derived, not estimated.
    kgrowth2  <- kgrowth1 / 1e7
    emax      <- exp(lemax + etalemax)
    ec50_a    <- exp(lec50_a + etalec50_a)
    ec50_b    <- exp(lec50_b + etalec50_b)
    alpha     <- exp(lalpha + etalalpha)
    beta      <- exp(lbeta + etalbeta)
    hill      <- exp(lhill + etalhill)
    delta1    <- exp(ldelta1)
    delta2    <- exp(ldelta2 + etaldelta2)
    kconv     <- exp(lkconv)
    degmax    <- exp(ldegmax + etaldegmax)
    km_deg    <- exp(lkm_deg)
    hill_deg  <- exp(lhill_deg + etalhill_deg)
    ic50_deg  <- exp(lic50_deg + etalic50_deg)

    # ---- 3. Plasma concentrations ----
    Cc_caz <- central_caz / vc_caz
    Cc_avi <- central_avi / vc_avi

    # ---- 4. Ceftazidime kill of the active population (Eq 6) ----
    # Avibactam lowers the ceftazidime EC50 biexponentially.
    ec50_caz <- ec50_a * exp(-alpha * Cc_avi) + ec50_b * exp(-beta * Cc_avi)
    kdeath1  <- emax * Cc_caz^hill / (ec50_caz^hill + Cc_caz^hill)

    # ---- 5. Delay functions (Eq 8) ----
    # I = 1 in the absence of ceftazidime (the drug-free growth control);
    # otherwise I ramps up from 0 as 1 - exp(-delta_i * t).
    i1 <- ifelse(Cc_caz <= 0, 1, 1 - exp(-delta1 * t))
    i2 <- ifelse(Cc_caz <= 0, 1, 1 - exp(-delta2 * t))

    # ---- 6. Bacterial populations (Eqs 5 and 7) ----
    # bact_active holds log10(P1); bact_resting holds P2 on the linear
    # CFU/mL scale (see the state-scale note above the compartmentData
    # block). The logistic bracket is shared by both equations and is
    # written on the total density log10(P1 + P2) in both.
    p1  <- 10^bact_active
    p2  <- max(bact_resting, 0)
    cap <- 1 - log10(p1 + p2) / log10nmax

    d/dt(bact_active)  <- i1 * kgrowth1 * cap * bact_active -
      (i2 * kdeath1 + kconv) * bact_active
    d/dt(bact_resting) <- kgrowth2 * cap * p2 + kconv * p1

    # ---- 7. Drug disposition, with bacterial degradation of CAZ (Eq 9) ----
    # Eq 9 is written on the ceftazidime CONCENTRATION; because vc_caz is
    # constant the same first-order rate applies to the amount. Km is the
    # log10-transformed CFU density (text under Eq 9), so the Hill term is
    # driven by bact_active rather than by the linear count p1.
    degrate <- degmax * bact_active^hill_deg /
      (km_deg^hill_deg + bact_active^hill_deg) *
      (1 - Cc_avi / (ic50_deg + Cc_avi))

    d/dt(central_caz) <- -(cl_caz / vc_caz) * central_caz - degrate * central_caz
    d/dt(central_avi) <- -(cl_avi / vc_avi) * central_avi

    # ---- 8. Initial conditions (Section 4.2, "Pharmacodynamic model") ----
    # P1: initial inoculum 10^6 CFU/mL -> 6 on the log10 state scale.
    # P2: initial inoculum 1/10^7 CFU/mL on the linear state scale.
    bact_active(0)  <- 6
    bact_resting(0) <- 1e-7

    # ---- 9. Outputs ----
    # Total bacterial density, the quantity reported in Tables 1 and 3
    # and plotted in Figures 1-4 panel 3.
    bactTotal <- log10(p1 + p2)

    Cc_caz ~ add(addSd_Cc_caz)
    Cc_avi ~ add(addSd_Cc_avi)
    bactTotal ~ add(addSd_bactTotal)
  })
}
