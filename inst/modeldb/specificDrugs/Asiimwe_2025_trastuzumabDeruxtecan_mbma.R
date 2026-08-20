Asiimwe_2025_trastuzumabDeruxtecan_mbma <- function() {
  description <- "MBMA. Two-compartment linear population PK model of trastuzumab deruxtecan (T-DXd, DS-8201, Enhertu, anti-HER2 antibody-drug conjugate) fitted by model-based meta-analysis to summary-level concentration-time data digitised from 4 published clinical trials in patients with HER2-positive breast, non-small-cell lung, gastric/GEJ, colorectal and other solid tumors. Between-study variability (BSV) is a study-level random effect on CL and Vc (block correlation 0.915) representing differences in inclusion criteria across trials and dose-escalation arms treated as separate studies. Residual error is proportional + additive on serum concentration (ng/mL) and was weighted by the square root of each trial's sample size during fitting; the tabulated parameter estimates are the unweighted values and simulation of a study of N subjects should scale the residual SD by 1/sqrt(N). Suitable simulation scope is study-arm-mean cycle-1 concentration-time profiles, NOT individual-patient concentrations. Parameter values are the T-DXd column of Asiimwe 2025 Table 1 (Monolix 2024R1)."

  reference <- paste(
    "Asiimwe IG, Chtiba N, Mouksassi S, Pillai G(C), Peter RM, Yuen E, Pilla Reddy V.",
    "Postmarketing Assessment of Antibody-Drug Conjugates:",
    "Proof-of-Concept Using Model-Based Meta-Analysis and a Clinical Utility Index Approach.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(11):1957-1970.",
    "doi:10.1002/psp4.70013",
    sep = " "
  )
  vignette <- "Asiimwe_2025_trastuzumab_ADCs_mbma"
  units <- list(
    time          = "day",
    dosing        = "mg (IV infusion; approved regimens 5.4 mg/kg q3w for HER2-positive breast cancer and NSCLC, 6.4 mg/kg q3w for gastric/GEJ cancer; infused over 90 min in Cycle 1 and 30 min in subsequent cycles)",
    concentration = "ng/mL (serum concentration of the intact antibody-drug conjugate T-DXd; residual error terms and reported constant SD are in ng/mL)"
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 4L,
    weight_median  = "59.0 kg (median across the 4 T-DXd Pop-PK studies; used when computing the mg-per-patient dose from a mg/kg regimen)",
    disease_state  = "HER2-positive advanced or metastatic solid tumors: breast cancer (dominant), non-small-cell lung, gastric / gastro-oesophageal junction, colorectal, and other HER2-expressing tumors (see Asiimwe 2025 Table S3 for the 4 Pop-PK studies)",
    dose_range     = "0.8 to 8.0 mg/kg IV; three-weekly (q3w) and weekly dosing were included for the Pop-PK fit (Cycle 1 data only; Cycle 1 approximates steady state given the ~6 day elimination half-life and 21-day cycles)",
    regions        = "global (studies included cohorts from North America, Europe, and Japan; see Asiimwe 2025 Table S3)",
    design         = "MBMA of summary-level data (mean concentrations per timepoint per study arm), weighted by the square root of each trial's sample size. Different doses within the same study were treated as separate 'studies' (Between-Treatment-Arm Variability absorbed into BSV) to match Phase I dose-escalation designs.",
    reference_subject = "median body weight 59.0 kg, no covariate effects modelled (weight enters only through dose = mg/kg x WT); no other covariates were tested in the primary Pop-PK model (Asiimwe 2025 Methods 2.2).",
    notes          = "Between-study variability (BSV) here is a study-level random effect, NOT individual between-subject variability. The T-DXd BSV magnitudes are small (omega_CL = 0.084, omega_Vc = 0.038) with high RSEs (48.3% / 38.8%) because only 4 studies were available; the sensitivity model with BTAV separated from BSV did not converge (Asiimwe 2025 Results 3.2). The paper also attempted a nonlinear-clearance / target-mediated model and a payload (DXd) sub-model but neither converged with acceptable precision; only the linear conjugate model is reproduced here."
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "trastuzumabDeruxtecan", units = NA_character_, specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "trastuzumabDeruxtecan", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  ini({
    # ============================================================
    # Two-compartment linear MBMA Pop-PK model (Asiimwe 2025 Methods 2.2,
    # Table 1 T-DXd column). Structural parameters in L / day; residual
    # error and additive constant error in ng/mL (paper's stated units).
    # Between-study random effects (log-normal) on CL and Vc are
    # correlated (rho = 0.915). Fit in Monolix 2024R1 via the
    # lixoftConnectors R interface.
    # ============================================================

    # ----- Structural PK parameters (Asiimwe 2025 Table 1, T-DXd) -----
    lcl <- log(0.585); label("Clearance CL at reference (L/day)")                       # Asiimwe 2025 Table 1 T-DXd: CL = 0.585 L/day (RSE 3.7%)
    lvc <- log(2.785); label("Central volume of distribution Vc at reference (L)")      # Asiimwe 2025 Table 1 T-DXd: Vc = 2.785 L (RSE 1.5%)
    lvp <- log(1.243); label("Peripheral volume of distribution Vp (L)")                # Asiimwe 2025 Table 1 T-DXd: Vp = 1.243 L (RSE 9.2%)
    lq  <- log(0.652); label("Intercompartmental clearance Q (L/day)")                  # Asiimwe 2025 Table 1 T-DXd: Q = 0.652 L/day (RSE 9.7%)

    # ============================================================
    # Between-study variability (BSV) as a correlated log-normal
    # random effect on CL and Vc (Asiimwe 2025 Methods 2.2 and Table 1).
    # Reported values are omega (standard deviation of the underlying
    # normal deviate); ini() takes variances on the diagonal and the
    # covariance off-diagonal, so:
    #   Var(eta_lcl) = 0.084^2 = 0.007056
    #   Var(eta_lvc) = 0.038^2 = 0.001444
    #   Cov         = rho * omega_cl * omega_vc
    #              = 0.915 * 0.084 * 0.038 = 0.002921
    # These are STUDY-level etas (MBMA), NOT individual between-subject
    # variability -- see pbpk-qsp-mbma.md. The BSV magnitudes are small
    # here because only 4 studies contributed and the reported RSEs
    # (~40-50%) reflect that limited data.
    # ============================================================
    eta_study_lcl + eta_study_lvc ~ c(0.007056,
                                      0.002921, 0.001444)  # Asiimwe 2025 Table 1 T-DXd: omega_CL = 0.084 (RSE 48.3%), omega_Vc = 0.038 (RSE 38.8%), corr(CL,Vc) = 0.915 (RSE 36.8%); variances and covariance computed as omega^2 and rho*omega_cl*omega_vc

    # ============================================================
    # Residual error (Asiimwe 2025 Methods 2.2 and Table 1 T-DXd).
    # Combined proportional + additive on Cc (ng/mL). Reported
    # 'weighted proportional error' 0.338 (33.8%) and 'constant error'
    # 3043.004 ng/mL are the point estimates BEFORE the sqrt(N)
    # per-trial weighting. Simulation for a study of size N should
    # scale these SDs by 1 / sqrt(N).
    # ============================================================
    propSd <- 0.338;     label("Proportional residual error SD (fraction; unweighted study-level)")  # Asiimwe 2025 Table 1 T-DXd: weighted proportional error = 0.338 (RSE 13.4%)
    addSd  <- 3043.004;  label("Additive residual error SD on Cc (ng/mL; unweighted study-level)")   # Asiimwe 2025 Table 1 T-DXd: constant error = 3043.004 ng/mL (RSE 8.2%)
  })

  model({
    # Individual (study-level) PK parameters. eta_study_* are the
    # MBMA study-level random effects on log-CL and log-Vc; there
    # is no covariate model in the primary Asiimwe 2025 T-DXd fit.
    cl <- exp(lcl + eta_study_lcl)
    vc <- exp(lvc + eta_study_lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                    k12 * central - k21 * peripheral1

    # Central compartment amount in mg; volume in L -> mg/L = ug/mL.
    # Reported concentrations and additive residual are in ng/mL, so
    # multiply by 1000 to keep Cc in ng/mL (matches Asiimwe 2025 Table 1).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
