Asiimwe_2025_trastuzumabEmtansine_mbma <- function() {
  description <- "MBMA. Two-compartment linear population PK model of trastuzumab emtansine (T-DM1, anti-HER2 antibody-drug conjugate) fitted by model-based meta-analysis to summary-level concentration-time data digitised from 14 published clinical trials in patients with HER2-positive metastatic breast cancer and other HER2-positive solid tumors. Between-study variability (BSV) is a study-level random effect on CL and Vc (block correlation 0.825) representing differences in inclusion criteria and the pooling of Phase I dose-escalation arms as separate 'studies'. Residual error is proportional + additive on serum concentration (ng/mL) and was weighted by the square root of each trial's sample size during fitting; the tabulated parameter estimates are the unweighted values and simulation of a study of N subjects should scale the residual SD by 1/sqrt(N). Suitable simulation scope is study-arm-mean cycle-1 concentration-time profiles, NOT individual-patient concentrations. Parameter values are the T-DM1 column of Asiimwe 2025 Table 1 (Monolix 2024R1)."

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
    dosing        = "mg (IV infusion; typical 3.6 mg/kg q3w for HER2-positive breast cancer, infused over 90 min in Cycle 1 and 30 min in subsequent cycles)",
    concentration = "ng/mL (serum concentration of the intact antibody-drug conjugate T-DM1; residual error terms and reported constant SD are in ng/mL)"
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 14L,
    weight_median  = "69.4 kg (median across the 14 T-DM1 Pop-PK studies; used when computing the mg-per-patient dose from a mg/kg regimen)",
    disease_state  = "predominantly HER2-positive metastatic breast cancer; a minority of arms include HER2-positive urothelial/bladder, pancreatic/cholangiocarcinoma, and other HER2-positive solid tumors (see Asiimwe 2025 Table S1 for the 14 Pop-PK studies)",
    dose_range     = "0.3 to 4.8 mg/kg IV; three-weekly (q3w) and weekly dosing were included for the Pop-PK fit (Cycle 1 data only; Cycle 1 approximates steady state given the ~4 day elimination half-life and 21-day cycles)",
    regions        = "global (studies included cohorts from North America, Europe, and Japan/China; see Asiimwe 2025 Table S1)",
    design         = "MBMA of summary-level data (mean concentrations per timepoint per study arm), weighted by the square root of each trial's sample size. Different doses within the same study were treated as separate 'studies' (Between-Treatment-Arm Variability absorbed into BSV) to match Phase I dose-escalation designs where the same patients are not randomised across arms.",
    reference_subject = "median body weight 69.4 kg, no covariate effects modelled (weight enters only through dose = mg/kg x WT); no other covariates were tested in the primary Pop-PK model (Asiimwe 2025 Methods 2.2).",
    notes          = "Between-study variability (BSV) here is a study-level random effect, NOT individual between-subject variability. The paper also fitted a two-compartment nonlinear-clearance / target-mediated model and a payload (DM1) sub-model but neither converged with acceptable precision; only the linear conjugate model is reproduced here. Parameter values in the sensitivity analyses that separated BTAV from BSV (Table S5) or restricted to dose-finding studies (Table S6) were similar to the primary values, with wider precision."
  )

  ini({
    # ============================================================
    # Two-compartment linear MBMA Pop-PK model (Asiimwe 2025 Methods 2.2,
    # Table 1 T-DM1 column). Structural parameters in L / day; residual
    # error and additive constant error in ng/mL (paper's stated units).
    # Between-study random effects (log-normal) on CL and Vc are
    # correlated (rho = 0.825). Fit in Monolix 2024R1 via the
    # lixoftConnectors R interface.
    # ============================================================

    # ----- Structural PK parameters (Asiimwe 2025 Table 1, T-DM1) -----
    lcl <- log(0.809); label("Clearance CL at reference (L/day)")                       # Asiimwe 2025 Table 1 T-DM1: CL = 0.809 L/day (RSE 7.3%)
    lvc <- log(3.283); label("Central volume of distribution Vc at reference (L)")      # Asiimwe 2025 Table 1 T-DM1: Vc = 3.283 L (RSE 4.8%)
    lvp <- log(0.748); label("Peripheral volume of distribution Vp (L)")                # Asiimwe 2025 Table 1 T-DM1: Vp = 0.748 L (RSE 8.3%)
    lq  <- log(1.120); label("Intercompartmental clearance Q (L/day)")                  # Asiimwe 2025 Table 1 T-DM1: Q = 1.120 L/day (RSE 26.3%)

    # ============================================================
    # Between-study variability (BSV) as a correlated log-normal
    # random effect on CL and Vc (Asiimwe 2025 Methods 2.2 and Table 1).
    # Reported values are omega (standard deviation of the underlying
    # normal deviate); ini() takes variances on the diagonal and the
    # covariance off-diagonal, so:
    #   Var(eta_lcl) = 0.334^2 = 0.111556
    #   Var(eta_lvc) = 0.221^2 = 0.048841
    #   Cov         = rho * omega_cl * omega_vc
    #              = 0.825 * 0.334 * 0.221 = 0.0609
    # These are STUDY-level etas (MBMA), NOT individual between-subject
    # variability -- see pbpk-qsp-mbma.md.
    # ============================================================
    eta_study_lcl + eta_study_lvc ~ c(0.111556,
                                      0.0609, 0.048841)  # Asiimwe 2025 Table 1 T-DM1: omega_CL = 0.334 (RSE 15.5%), omega_Vc = 0.221 (RSE 14.8%), corr(CL,Vc) = 0.825 (RSE 10.1%); variances and covariance computed as omega^2 and rho*omega_cl*omega_vc

    # ============================================================
    # Residual error (Asiimwe 2025 Methods 2.2 and Table 1 T-DM1).
    # Combined proportional + additive on Cc (ng/mL). Reported
    # 'weighted proportional error' 0.430 (43%) and 'constant error'
    # 2633.766 ng/mL are the point estimates BEFORE the sqrt(N)
    # per-trial weighting. Simulation for a study of size N should
    # scale these SDs by 1 / sqrt(N).
    # ============================================================
    propSd <- 0.430;     label("Proportional residual error SD (fraction; unweighted study-level)")  # Asiimwe 2025 Table 1 T-DM1: weighted proportional error = 0.430 (RSE 9.1%)
    addSd  <- 2633.766;  label("Additive residual error SD on Cc (ng/mL; unweighted study-level)")   # Asiimwe 2025 Table 1 T-DM1: constant error = 2633.766 ng/mL (RSE 6.1%)
  })

  model({
    # Individual (study-level) PK parameters. eta_study_* are the
    # MBMA study-level random effects on log-CL and log-Vc; there
    # is no covariate model in the primary Asiimwe 2025 T-DM1 fit.
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
