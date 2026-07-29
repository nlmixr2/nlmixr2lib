Almquist_2022_anifrolumab <- function() {
  description <- "Two-compartment QSS-TMDD population PK model for anifrolumab (anti-IFNAR1 IgG1-kappa) in healthy volunteers and adults with systemic lupus erythematosus (Almquist 2022): linear plus quasi-steady-state target-mediated elimination via a dynamic IFNAR1 receptor pool, time-varying linear clearance (Emax-on-time), and IFNGS-high/low and body-weight covariate effects."
  reference <- "Almquist J, Kuruvilla D, Mai T, et al. Nonlinear Population Pharmacokinetics of Anifrolumab in Healthy Volunteers and Patients With Systemic Lupus Erythematosus. J Clin Pharmacol. 2022;62(9):1106-1120. doi:10.1002/jcph.2055"
  vignette <- "Almquist_2022_anifrolumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline). Used for allometric scaling of CL and Vc with reference weight 69.1 kg (median of pooled population, Almquist 2022 Table S2).",
      source_name        = "BLWGHT"
    ),
    BGENE21_HIGH = list(
      description        = "Baseline 21-gene type I interferon signature status (IFNGS-high vs IFNGS-low)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (IFNGS-high)",
      notes              = "1 = baseline 21-gene IFN signature score >= 4-gene-test cut-off (IFNGS-high; reference category, factor on CL = 1); 0 = IFNGS-low (factor on CL = e_ifngslow_cl < 1, so CL is reduced for IFNGS-low subjects). The 4-gene IFNGS test (IFI27, IFI44, IFI44L, RSAD2) classifies high vs low using an analytically validated cut-off; >=75% of patients across the SLE studies were IFNGS-high (Almquist 2022 Methods, 21-Gene IFNGS PD Assay section).",
      source_name        = "BGENEFLG"
    )
  )

  population <- list(
    n_subjects     = 670L,
    n_sle          = 664L,
    n_healthy      = 6L,
    n_studies      = 5L,
    studies        = "Phase 1 healthy-volunteer SAD (NCT02601625), Phase 2 study in Japanese patients (NCT01559090), Phase 2b MUSE (NCT01438489), Phase 3 TULIP-1 (NCT02446912), Phase 3 TULIP-2 (NCT02446899).",
    age_range      = "18-69 years (pooled)",
    age_median     = "41 years (pooled, Table S2)",
    weight_range   = "42.0-144.8 kg (pooled)",
    weight_median  = "68.6 kg (pooled, Table S2)",
    sex_female_pct = 92.5,
    race_ethnicity = c(White = 57.2, Black = 13.6, Asian = 11.0, Other = 17.0),
    disease_state  = "Moderate-to-severe systemic lupus erythematosus (664 SLE patients) plus 6 healthy adult volunteers (single-dose phase 1).",
    dose_range     = "100-1000 mg IV every 4 weeks across SLE studies; single-dose 30 mg SC, 300 mg IV/SC, 600 mg SC in the phase 1 study (subcutaneous arms excluded from popPK).",
    regions        = "Global (North America, Europe, Asia including Japan).",
    ifngs_high_pct = "75.5-88.2% across SLE studies",
    ada_status     = "ADA prevalence (any visit) 6.9% (46/670) overall; no significant impact on PK in pooled TULIP analysis.",
    notes          = "Demographics and dosing summarised in Almquist 2022 Tables S1-S2. The popPK analysis dataset comprises 6049 anifrolumab serum concentrations; PK from SC arms in phase 1 was excluded."
  )

  ini({
    # Structural PK -- typical values for the IFNGS-high reference subject at 69.1 kg (Almquist 2022 Table 1)
    lcl    <- log(0.193);  label("Linear systemic clearance for IFNGS-high reference subject at 69.1 kg (CL, L/day)")  # Table 1: CL = 0.193 L/day
    lvc    <- log(2.93);   label("Central volume of distribution at 69.1 kg (Vc, L)")                                  # Table 1: Vc = 2.93 L
    lq     <- log(0.937);  label("Inter-compartmental clearance (Q, L/day)")                                           # Table 1: Q = 0.937 L/day
    lvp    <- log(3.30);   label("Peripheral volume of distribution (Vp, L)")                                          # Table 1: Vp = 3.30 L

    # QSS-TMDD binding parameters (target = type I IFN receptor IFNAR1, in nmol/L)
    lkss   <- log(0.712);  label("Quasi-steady-state binding constant Kss = (koff + kint)/kon (nmol/L)")               # Table 1: Kss = 0.712 nmol/L
    lrbase    <- log(0.0999); label("Baseline IFNAR1 concentration R0 (nmol/L)")                                          # Table 1: R0 = 0.0999 nmol/L

    # Receptor turnover: kdeg fixed from confocal-imaging studies (Wang 2013, supplement ref 1).
    # Supplement appendix (line 140) sets the internalisation rate kint equal to kdeg = 77.4 d^-1
    # so the drug-target complex is removed at the same rate as the unbound receptor.
    # supplement-derived (JCPH appendix; not in main paper's results table 1 except as the kint=kdeg constraint)
    kdeg   <- fixed(77.4); label("Receptor degradation rate kdeg (1/day, fixed; equals internalisation rate kint per supplement)")  # Supplement appendix: kdeg = kint = 77.4 d^-1 fixed

    # Covariate effects on linear CL and Vc
    e_ifngslow_cl <- 0.793;  label("Multiplicative factor on CL for IFNGS-low subjects (1 = high reference, <1 = low)")  # Table 1: F_IFNGS-low = 0.793
    e_wt_cl       <- 0.601;  label("Allometric power exponent of (WT/69.1) on CL (unitless)")                            # Table 1: BW on CL = 0.601
    e_wt_vc       <- 0.764;  label("Allometric power exponent of (WT/69.1) on Vc (unitless)")                            # Table 1: BW on Vc = 0.764

    # Time-varying linear CL: empirical Emax-on-time multiplicative factor
    # F_emp(t) = exp((tmax + eta_tmax) * t / (tc50 + t))   per Almquist 2022 main-text Eq. for F_EMPIR (paper text reduces the supplement's Hill form to Hill power = 1)
    tmax   <- -0.155;  label("Maximal change in log(CL) over time (asymptotic; unitless)")                              # Table 1: Tmax = -0.155
    tc50   <- 380;     label("Time to reach half the maximal change in log(CL) (days)")                                 # Table 1: TC50 = 380 days

    # IIV (variances; Tmax has additive eta on the linear scale, all others are exp-additive on log-scale)
    etalcl   ~ 0.109   # Table 1: omega^2 = 0.109 (CV ~33.0%) on lcl
    etalvc   ~ 0.0723  # Table 1: omega^2 = 0.0723 (CV ~26.9%) on lvc
    etalrbase   ~ 0.0882  # Table 1: omega^2 = 0.0882 (CV ~29.7%) on lrbase
    etatmax  ~ 0.146   # Table 1: omega^2 = 0.146 (CV ~38.2%) on tmax (additive eta, NONMEM ETA(4))

    # Residual error (paper Table 1 reports SDs in ng/mL; this implementation works in ug/mL so addSd is /1000)
    propSd <- 0.305;     label("Proportional residual error (fraction)")                                                # Table 1: proportional SD = 0.305
    addSd  <- 0.0201;    label("Additive residual error (ug/mL); paper reports 20.1 ng/mL = 0.0201 ug/mL")              # Table 1: additive SD = 20.1 (ng/mL)
  })

  model({
    # Anifrolumab molecular weight (148 kDa = 0.148 mg/nmol). Used to convert between
    # drug mass (mg, the dosing/CL/Vc unit) and drug amount (nmol, the binding-kinetics unit).
    mwdrug <- 0.148  # mg/nmol; anifrolumab IgG1-kappa, ~148 kDa

    # Reference body weight (median of pooled popPK dataset, Table S2)
    wtref  <- 69.1   # kg

    # Covariate-corrected individual PK parameters
    f_ifngs <- e_ifngslow_cl^(1 - BGENE21_HIGH)               # 1 if BGENE21_HIGH=1 (high), else e_ifngslow_cl
    f_emp   <- exp((tmax + etatmax) * t / (tc50 + t))         # Emax-on-time multiplicative factor for linear CL

    cl    <- exp(lcl + etalcl) * f_ifngs * (WT / wtref)^e_wt_cl * f_emp
    vc    <- exp(lvc + etalvc) * (WT / wtref)^e_wt_vc
    q     <- exp(lq)
    vp    <- exp(lvp)
    kss   <- exp(lkss)                                        # nmol/L
    rbase    <- exp(lrbase + etalrbase)                                # nmol/L
    kint  <- kdeg                                             # supplement constraint: complex internalisation rate equals receptor degradation rate

    # Working concentrations (with explicit unit handling)
    cdrug_ugmL <- central     / vc                            # ug/mL = mg/L (free drug in central)
    cdrug_nM   <- cdrug_ugmL  / mwdrug                        # nmol/L (free drug in central)
    cperi_ugmL <- peripheral1 / vp                            # ug/mL (free drug in peripheral)
    crtot_nM   <- total_target / vc                           # nmol/L (total IFNAR1 = free + bound)

    # Quasi-steady-state drug:target complex concentration (nmol/L)
    complex_nM <- crtot_nM * cdrug_nM / (kss + cdrug_nM)

    # QSS partitioning denominator: 1 + d[Ab*R]/dAb (Mager & Krzyzanski 2005; supplement appendix Eq. for DADT(2))
    qss_denom  <- 1 + crtot_nM * kss / (kss + cdrug_nM)^2

    # Aggregate free-drug mass-loss rate from the central compartment (mg/day):
    #   linear CL out of Vc + complex internalisation + intercompartmental flux to peripheral
    # Internalisation: kint (1/day) * complex_nM (nmol/L) * vc (L) * mwdrug (mg/nmol) = mg/day
    ab_total_loss <- cl * cdrug_ugmL +
                     kint * complex_nM * vc * mwdrug +
                     q * (cdrug_ugmL - cperi_ugmL)

    # ODE system: free drug central + peripheral; total receptor pool with linear turnover.
    # Receptor pool: dRtot/dt = kdeg * (R0*Vc - Rtot)  -- supplement appendix DADT(4) under kint = kdeg.
    d/dt(central)      <- -ab_total_loss / qss_denom
    d/dt(peripheral1)  <-  q * (cdrug_ugmL - cperi_ugmL)
    d/dt(total_target) <-  kdeg * (rbase * vc - total_target)
    total_target(0)    <-  rbase * vc                            # nmol; supplement appendix A0(4) = R0*V2

    # Observation: free anifrolumab serum concentration in ug/mL
    Cc <- cdrug_ugmL
    Cc ~ add(addSd) + prop(propSd)
  })
}
