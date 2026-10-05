Srimani_2022_ixazomib_mprotein <- function() {
  description <- "Joint exposure-efficacy model for oral ixazomib plus lenalidomide-dexamethasone (LenDex) in relapsed/refractory multiple myeloma from the phase III TOURMALINE-MM1 trial (Srimani 2022). Ixazomib plasma concentrations come from the three-compartment population PK model of Gupta 2017 (fixed). Serum M-protein is the sum of a drug-sensitive population, relative to the observed baseline, that follows a type-1 indirect response towards a concentration-dependent steady-state nadir, and a resistant population that grows exponentially from the individual time of nadir. Three hazard sub-models are carried as cumulative-hazard states: an empirical dropout hazard driven by the relative steady-state M-protein and by M-protein rising above 116% of baseline, a log-logistic accelerated-failure-time relapse hazard whose scale depends on ixazomib concentration and on the M-protein response rate constant, and a progression-free-survival hazard driven by the M-protein response parameters, baseline M-protein and the weekly ixazomib AUC. Covariates are high-risk cytogenetics (on the response rate constant) and prior immunomodulatory-drug therapy (on the steady-state nadir)."
  reference <- paste(
    "Srimani JK, Diderichsen PM, Hanley MJ, Venkatakrishnan K, Labotka R, Gupta N.",
    "Population pharmacokinetic/pharmacodynamic joint modeling of ixazomib efficacy",
    "and safety using data from the pivotal phase III TOURMALINE-MM1 study in",
    "multiple myeloma patients. CPT Pharmacometrics Syst Pharmacol.",
    "2022;11(8):1085-1099. doi:10.1002/psp4.12815.",
    "Ixazomib PK layer: Gupta N, Diderichsen PM, Hanley MJ, et al. Clin Pharmacokinet.",
    "2017;56(11):1355-1368. doi:10.1007/s40262-017-0526-4",
    "(also available as modellib('Gupta_2017_ixazomib'))."
  )
  vignette <- "Srimani_2022_ixazomib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL", mprotein = "g/L")

  paper_specific_compartments <- c("mprotein_rel", "cumhaz_relapse")

  compartmentData <- list(
    depot = list(analyte = "ixazomib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ixazomib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ixazomib", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "ixazomib", units = "mg", specimen = "tissue", verified = TRUE),
    auc_central = list(analyte = "ixazomib", units = "ng*h/mL", specimen = "not applicable", verified = TRUE),
    mprotein_rel = list(
      analyte = "drug-sensitive serum M-protein relative to baseline",
      units = "(fraction of baseline)",
      specimen = "serum",
      verified = TRUE
    ),
    cumhaz_drop = list(
      analyte = "dropout cumulative hazard",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    ),
    cumhaz_relapse = list(
      analyte = "relapse cumulative hazard",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    ),
    cumhaz_pfs = list(
      analyte = "progression-free-survival cumulative hazard",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area at baseline.",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only by the Gupta 2017 ixazomib PK layer: power covariate on the second peripheral volume, reference 1.87 m^2, exponent 2.06.",
      source_name = "BSA"
    ),
    MCPROT = list(
      description = "Observed baseline serum M-protein concentration (Y_BL).",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline, g/L (not the register default g/dL). Scales the relative M-protein states to g/L (Srimani 2022 Equation 2), sets the M-protein threshold of the dropout hazard (alpha_BL * Y_BL, Equation 8) and enters the PFS hazard as a power covariate centred on the efficacy-dataset median of 23 g/L (Equation 13; Table 1). Patients entered the efficacy dataset only with a baseline M-protein of at least 10 g/L.",
      source_name = "Y_BL"
    ),
    TUM_CYTOGENETIC_HIGH_RISK = list(
      description = "High-risk cytogenetics indicator (1 = high risk, 0 = standard risk or risk not available).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (standard cytogenetic risk or risk unavailable; Srimani 2022 Equation 6).",
      notes = "Linear multiplicative effect on the M-protein response rate constant kR: kR_i = theta_kR * (1 + 0.590 * TUM_CYTOGENETIC_HIGH_RISK) (Equation 6; Table 2 row 'kR (BCYABCAT)'). 20.8% of the efficacy dataset was high risk (Table 1). TOURMALINE-MM1 defined high risk as del(17p), t(4;14) or t(14;16).",
      source_name = "BCYABCAT"
    ),
    PRIOR_IMID = list(
      description = "Prior immunomodulatory-drug (IMiD: thalidomide, lenalidomide, pomalidomide) therapy indicator (1 = exposed, 0 = IMiD-naive).",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (IMiD-exposed); the paper's effect is on the naive group.",
      notes = "Linear multiplicative effect on the relative steady-state M-protein nadir: Yss_i = theta_Yss * (1 - 0.427 * (1 - PRIOR_IMID)), i.e. 42.7% lower in IMiD-naive patients (Equation 7; Table 2 row 'Yss (PIMID)'). 55.9% of the efficacy dataset was IMiD-exposed (Table 1).",
      source_name = "PIMID"
    ),
    T_NADIR = list(
      description = "Individual time of the serum M-protein nadir after the first dose; the resistant M-protein population grows from this time.",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = "Not a model parameter: Srimani 2022 derived t_nadir for each patient before the PK/PD analysis by fitting exponential functions to the observed M-protein profile and taking the time of the minimum fitted value plus 0.5 g/L (Methods, M-protein PK/PD model). For simulation it can be drawn from the relapse hazard carried by this model (the relapse sub-model is the time-to-event model that stands in for the nadir time; see the vignette). Set it beyond the end of the simulation for a patient whose M-protein never regrows.",
      source_name = "t_nadir"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 467L,
    n_studies = 1L,
    age_range = "40-91 years (median 66)",
    age_median = "66 years",
    sex_female_pct = 43.7,
    race_ethnicity = c(White = 87.6, Black = 1.5, Asian = 7.3, `Not reported` = 2.4, Other = 1.3),
    disease_state = "Relapsed and/or refractory multiple myeloma after 1-3 prior lines of therapy, with measurable serum M-protein (baseline at least 10 g/L and at least 3 M-protein observations for this efficacy dataset).",
    dose_range = "Ixazomib 4 mg or matching placebo orally on days 1, 8 and 15 of 28-day cycles, both arms with lenalidomide 25 mg (10 mg for reduced creatinine clearance) on days 1-21 and dexamethasone 40 mg on days 1, 8, 15 and 22, until progression or unacceptable toxicity.",
    regions = "Global (TOURMALINE-MM1, C16010).",
    baseline_labs = "Baseline M-protein 23.0 g/L (10.0-102); creatinine clearance 80.8 mL/min (22.9-231); hemoglobin 113 g/L (68.0-167). Medians [range] from Srimani 2022 Table 1.",
    biomarkers = "High-risk cytogenetics 20.8%, standard 55.0%, not available 24.2%; ISS stage III 13.3%; prior IMiD 55.9%; prior proteasome inhibitor 67.9%; 1 prior line 59.7% (Table 1).",
    notes = "Exposure-efficacy dataset: 467 of the 720 treated patients (240 ixazomib, 227 placebo). Demographics from Srimani 2022 Table 1."
  )

  ini({
    # Ixazomib PK -- Gupta 2017 Table 3, fixed: Srimani 2022 used Bayesian
    # (post-hoc) individual PK parameters from this model rather than
    # re-estimating it (Methods, Pharmacokinetic/pharmacodynamic modeling).
    lka <- fixed(log(0.34)); label("Ixazomib first-order absorption rate constant (1/h)") # Gupta 2017 Table 3: Ka = 0.34/h
    lcl <- fixed(log(1.86)); label("Ixazomib clearance (L/h)") # Gupta 2017 Table 3: CL = 1.86 L/h
    lvc <- fixed(log(13.7)); label("Ixazomib central volume (L)") # Gupta 2017 Table 3: V2 = 13.7 L
    lfdepot <- fixed(log(0.58)); label("Ixazomib oral bioavailability (fraction)") # Gupta 2017 Table 3: F = 58%
    lq <- fixed(log(5.18)); label("Ixazomib intercompartmental clearance to peripheral1 (L/h)") # Gupta 2017 Table 3: Q3 = 5.18 L/h
    lvp <- fixed(log(309)); label("Ixazomib first peripheral volume (L)") # Gupta 2017 Table 3: V3 = 309 L
    lq2 <- fixed(log(26.1)); label("Ixazomib intercompartmental clearance to peripheral2 (L/h)") # Gupta 2017 Table 3: Q4 = 26.1 L/h
    lvp2 <- fixed(log(205)); label("Ixazomib second peripheral volume (L)") # Gupta 2017 Table 3: V4 = 205 L
    ltlag <- fixed(log(13 / 60)); label("Ixazomib absorption lag time (h)") # Gupta 2017 Table 3: TLAG = 13 min
    e_bsa_vp2 <- fixed(2.06); label("Power exponent of BSA on the second peripheral volume (unitless)") # Gupta 2017 Table 3: V4[BSA] = 2.06
    # Gupta 2017 Table 3 IIV as %CV, omega^2 = log(1 + CV^2); CL-F correlation 0.82.
    etalcl + etalfdepot ~ fixed(c(0.17697, 0.22550, 0.42726)) # Gupta 2017 Table 3: CL 44%CV, F 73%CV, rho 0.82
    etalvp2 ~ fixed(0.48492) # Gupta 2017 Table 3: V4 79%CV

    # M-protein PK/PD model -- Srimani 2022 Table 2. NONMEM time unit is the
    # hour: exp(-6.70) * 168 = 0.207/week matches the printed 0.206/week.
    lkr <- -6.70; label("Log M-protein response rate constant kR of the drug-sensitive population (log 1/h)") # Table 2 kR = -6.70 (0.206/week)
    lyss <- -1.95; label("Log relative steady-state M-protein nadir Yss without ixazomib (log fraction of baseline)") # Table 2 Yss = -1.95 (14.3%)
    lkl <- -9.78; label("Log growth rate constant kL of the resistant population after the nadir (log 1/h)") # Table 2 KL = -9.78 (0.00951/week)
    imax <- 0.758; label("Maximum fractional ixazomib inhibition of the steady-state nadir (fraction)") # Table 2 Imax = 0.758
    lic50 <- 1.19; label("Log ixazomib concentration at half-maximal inhibition (log ng/mL)") # Table 2 IC50 = 1.19 (3.29 ng/mL)
    e_cytohr_kr <- 0.590; label("Fractional change in kR for high-risk cytogenetics (unitless)") # Table 2 kR (BCYABCAT) = 0.590
    e_imidnaive_yss <- -0.427; label("Fractional change in Yss for IMiD-naive patients (unitless)") # Table 2 Yss (PIMID) = -0.427

    # Table 2 prints the variance; the %CV column is its square root
    # (sqrt(0.655) = 0.809 -> 81.0%CV).
    etalkr ~ 0.655 # Table 2 IIV kR = 0.655 (81.0%CV)
    etalyss ~ 2.39 # Table 2 IIV Yss = 2.39 (155%CV)
    etalkl ~ 1.34 # Table 2 IIV KL = 1.34 (116%CV)

    propSd_mprotein <- 0.218; label("Proportional residual error on M-protein (fraction)") # Table 2 Prop. Error = 0.218
    addSd_mprotein <- fixed(0.5); label("Additive residual error on M-protein (g/L)") # Table 2 Add. Error = 0.500 (fix)

    # Dropout hazard -- Equations 8-9, Table 2 'Dropout model'.
    llambda0_drop <- -10.1; label("Log baseline dropout hazard (log 1/h)") # Table 2 lambda0 = -10.1 (0.00697/week)
    lambda_rss_drop <- 0.0000125; label("Dropout hazard per unit relative steady-state M-protein Rss (1/h)") # Table 2 lambda_RSS = 0.0000125
    let50_drop <- 8.14; label("Log time to half-maximal dropout onset after t0 (log h)") # Table 2 ET50 = 8.14 (20.4 week)
    lt0_drop <- 6.86; label("Log dropout hazard lag time t0 (log h)") # Table 2 t0 = 6.86 (5.68 week)
    llambda_mprot_drop <- -6.60; label("Log dropout hazard per unit M-protein relative to baseline above the threshold (log 1/h)") # Table 2 lambda_M-protein = -6.60 (0.228/week)
    lalpha_bl_drop <- 0.148; label("Log M-protein threshold relative to baseline that switches on the M-protein dropout hazard (log fraction)") # Table 2 alpha_BL = 0.148 (116%)

    # Relapse hazard -- Equations 10-12, Table 2 'Relapse model'.
    llambda0_relapse <- -1.23; label("Log multiplier on the log-logistic relapse hazard (log unitless)") # Table 2 lambda0 = -1.23
    lalpha_relapse <- 10.9; label("Log baseline log-logistic scale alpha0 of time to relapse (log h)") # Table 2 alpha = 10.9 (321 week)
    lbeta_relapse <- 0.938; label("Log log-logistic shape beta of time to relapse (log unitless)") # Table 2 beta = 0.938 (2.56)
    probit_tnadir0 <- -1.67; label("Probit of the fraction of patients with no estimable nadir time (probit)") # Table 2 P_tnadir,0 = -1.67 (4.71%)
    alpha_ixa_relapse <- 0.0509; label("Linear ixazomib-concentration effect on the relapse scale alpha (mL/ng)") # Table 2 alpha_IXA = 0.0509
    alpha_kr_relapse <- -2960; label("Exponential kR coefficient on the relapse scale alpha (h)") # Table 2 alpha_kR = -2960 (-17.6 week)
    alpha_kr0_relapse <- 0.0509; label("Intercept added to the exponential kR term of the relapse scale alpha (unitless)") # Table 2 alpha_kR,0 = 0.0509

    # PFS hazard -- Equation 13, Table 2 'PFS model'.
    llambda0_pfs <- -9.21; label("Log baseline PFS hazard for the median patient (log 1/h)") # Table 2 lambda0 = -9.21 (0.0168/week)
    lambda_ixa_pfs <- -0.00108; label("Log-hazard coefficient of weekly ixazomib AUC on PFS (mL/(ng*h))") # Table 2 lambda_IXA = -0.00108
    let50_pfs <- 7.82; label("Log time to half-maximal PFS hazard onset (log h)") # Table 2 ET50 = 7.82 (14.8 week)
    e_kr_pfs <- 0.798; label("Power exponent of kR on the PFS hazard (unitless)") # Table 2 lambda0 (kR) = 0.798
    e_yss_pfs <- 0.680; label("Power exponent of Yss on the PFS hazard (unitless)") # Table 2 lambda0 (Yss) = 0.680
    e_mcprot_pfs <- 0.784; label("Power exponent of baseline M-protein on the PFS hazard (unitless)") # Table 2 lambda0 (M-protein) = 0.784
  })

  model({
    # Ixazomib PK (Gupta 2017)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2 + etalvp2) * (BSA / 1.87)^e_bsa_vp2
    fdepot <- exp(lfdepot + etalfdepot)
    tlag <- exp(ltlag)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (cl / vc) * central -
      (q / vc) * central + (q / vp) * peripheral1 -
      (q2 / vc) * central + (q2 / vp2) * peripheral2
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1
    d/dt(peripheral2) <- (q2 / vc) * central - (q2 / vp2) * peripheral2
    f(depot) <- fdepot
    alag(depot) <- tlag

    # mg / L * 1000 = ng/mL
    Cc <- 1000 * central / vc

    # Weekly AUC (ng*h/mL) over the trailing 168 h, the exposure metric of
    # the PFS hazard. delay() returns the zero initial history before 168 h.
    # It forces the non-stiff dense dop853 solver, which can fail for an
    # occasional subject at the default tolerances; atol = rtol = 1e-6 in
    # rxSolve() avoids that.
    d/dt(auc_central) <- Cc
    aucwk <- auc_central - delay(auc_central, 168)

    # M-protein (Equations 1-7)
    kr <- exp(lkr + etalkr) * (1 + e_cytohr_kr * TUM_CYTOGENETIC_HIGH_RISK)
    yss <- exp(lyss + etalyss) * (1 + e_imidnaive_yss * (1 - PRIOR_IMID))
    kl <- exp(lkl + etalkl)
    ic50 <- exp(lic50)
    rss <- yss * (1 - imax * Cc / (ic50 + Cc))

    mprotein_rel(0) <- 1
    d/dt(mprotein_rel) <- kr * (rss - mprotein_rel)

    tpost <- (t - T_NADIR) * (t > T_NADIR)
    rplus <- exp(kl * tpost) - 1
    mprotein <- MCPROT * (mprotein_rel + rplus)

    # Dropout hazard (Equations 8-9): the whole sum is modulated by the onset
    # function gamma(t), which is zero until t0.
    t0_drop <- exp(lt0_drop)
    et50_drop <- exp(let50_drop)
    tstar_drop <- (t - t0_drop) * (t > t0_drop)
    gamma_drop <- tstar_drop / (et50_drop + tstar_drop)
    alpha_bl_drop <- exp(lalpha_bl_drop)
    hazard_mprot_drop <- exp(llambda_mprot_drop) * (mprotein / MCPROT) * (mprotein > alpha_bl_drop * MCPROT)
    hazard_drop <- (exp(llambda0_drop) + lambda_rss_drop * rss + hazard_mprot_drop) * gamma_drop
    d/dt(cumhaz_drop) <- hazard_drop
    surv_drop <- exp(-cumhaz_drop)

    # Relapse hazard (Equations 10-12): log-logistic hazard whose scale alpha
    # is accelerated by the instantaneous ixazomib concentration and by kR.
    alpha_relapse <- exp(lalpha_relapse) * (1 + alpha_ixa_relapse * Cc) *
      (exp(alpha_kr_relapse * kr) + alpha_kr0_relapse)
    beta_relapse <- exp(lbeta_relapse)
    tsc_relapse <- t / alpha_relapse
    hazard_relapse <- exp(llambda0_relapse) * (beta_relapse / alpha_relapse) *
      tsc_relapse^(beta_relapse - 1) / (1 + tsc_relapse^beta_relapse)
    d/dt(cumhaz_relapse) <- hazard_relapse
    surv_relapse <- exp(-cumhaz_relapse)
    p_tnadir0 <- phi(probit_tnadir0)

    # PFS hazard (Equation 13): covariates centred on the efficacy-dataset
    # medians kR = 0.001488/h, Yss = 0.15 and baseline M-protein = 23 g/L.
    et50_pfs <- exp(let50_pfs)
    hazard_pfs <- exp(llambda0_pfs) *
      (kr / 0.001488)^e_kr_pfs * (yss / 0.15)^e_yss_pfs * (MCPROT / 23)^e_mcprot_pfs *
      (t / (t + et50_pfs)) * exp(lambda_ixa_pfs * aucwk)
    d/dt(cumhaz_pfs) <- hazard_pfs
    surv_pfs <- exp(-cumhaz_pfs)

    mprotein ~ add(addSd_mprotein) + prop(propSd_mprotein)
  })
}
