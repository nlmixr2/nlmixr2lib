Talke_2018_dexmedetomidine <- function() {
  description <- "Three-compartment IV population PK plus effect-compartment sigmoid Emax PD model for dexmedetomidine-induced peripheral vasoconstriction (ADC units from finger photoplethysmography) in healthy adult volunteers, with a priori allometric body-weight scaling on CL, Q2, Q3 (exponent 0.75) and V1, V2, V3 (exponent 1) at a 70 kg reference weight (Talke and Anderson 2018, Tables 3 and 4)"
  reference <- "Talke P, Anderson BJ. Pharmacokinetics and pharmacodynamics of dexmedetomidine-induced vasoconstriction in healthy volunteers. Br J Clin Pharmacol. 2018;84(7):1364-1372. doi:10.1111/bcp.13571"
  vignette <- "Talke_2018_dexmedetomidine"
  units <- list(time = "min", dosing = "ug", concentration = "ug/L") # Methods + Tables: dose ug, plasma in ug/L (= ng/mL); time in minutes (CL in L/min, t1/2 keo in min)

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Used for a priori allometric scaling per Talke 2018 Methods: CL, Q2, Q3 scale as (WT/70)^0.75 and V1, V2, V3 scale as (WT/70)^1; reference weight 70 kg.",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 10L,                        # Talke 2018 Abstract / Methods / Table 1
    n_studies      = 1L,                         # Single-centre study at UCSF
    age_range      = "21-36 years",              # Talke 2018 Table 1 (mean 29, SD 4)
    age_median     = "29 years (mean)",          # Talke 2018 Table 1
    weight_range   = "52-89 kg",                 # Talke 2018 Table 1 (mean 72, SD 13)
    weight_median  = "72 kg (mean)",             # Talke 2018 Table 1
    sex_female_pct = 50,                         # Talke 2018 Table 1: 5 male / 5 female
    race_ethnicity = NULL,                       # Not reported in Talke 2018 Table 1
    disease_state  = "Healthy adult volunteers; exclusions: history of cardiac, pulmonary, hepatic or renal disease; alcohol or drug abuse; prescription medications; age >45 years; weight >130% of normal.",
    dose_range     = "Single computer-controlled (STANPUMP) IV infusion of dexmedetomidine over 15 min targeting plasma concentration 0.3 ng/mL; observed cumulative dose 0.28 ug/kg.",
    regions        = "USA (single centre, University of California San Francisco).",
    n_observations = "120 plasma concentrations for PK; 4560 ADC observations for PD.",
    notes          = "BMI 20-27 kg/m^2 (mean 23, SD 2), height 160-183 cm (mean 174, SD 9). Sympathetic fibres of the left arm were blocked with axillary perivascular brachial plexus block (30 mL 1% mepivacaine) prior to dexmedetomidine infusion to isolate the direct peripheral vasoconstrictive effect (Talke 2018 Methods)."
  )

  ini({
    # Structural PK parameters - reference weight 70 kg (Talke 2018 Methods + Table 3A)
    lcl  <- log(0.751); label("Clearance CL at 70 kg (L/min)")                          # Talke 2018 Table 3A: CL = 0.751 L/min/70kg (SE 4.6%; bootstrap 0.712, 95% CI 0.107-0.872)
    lvc  <- log(12);    label("Central volume of distribution V1 at 70 kg (L)")         # Talke 2018 Table 3A: V1 = 12 L/70kg (SE 6.3%; bootstrap 12.8, 95% CI 9.57-16.9)
    lq   <- log(1.01);  label("Inter-compartmental clearance Q2 at 70 kg (L/min)")      # Talke 2018 Table 3A: Q2 = 1.01 L/min/70kg (SE 11.2%; bootstrap 0.96, 95% CI 0.531-1.816)
    lvp  <- log(11.7);  label("Peripheral volume of distribution V2 at 70 kg (L)")      # Talke 2018 Table 3A: V2 = 11.7 L/70kg (SE 13.1%; bootstrap 13.1, 95% CI 6.369-35.61)
    lq2  <- log(1.03);  label("Inter-compartmental clearance Q3 at 70 kg (L/min)")      # Talke 2018 Table 3A: Q3 = 1.03 L/min/70kg (SE 3.7%; bootstrap 1.06, 95% CI 0.674-1.39); no IIV reported
    lvp2 <- log(47.2);  label("Second peripheral volume of distribution V3 at 70 kg (L)") # Talke 2018 Table 3A: V3 = 47.2 L/70kg (SE 11.2%; bootstrap 48.4, 95% CI 35.295-217.025); no IIV reported

    # PD parameters - sigmoid Emax effect-compartment model (Talke 2018 Methods + Table 4A)
    le0      <- log(5860);  label("Baseline ADC effect E0 (ADC units)")                     # Talke 2018 Table 4A: E0 = 5860 ADC units (bootstrap 5666, 95% CI 4549-7527)
    lemax    <- log(3090);  label("Maximum vasoconstriction effect EMAX (ADC units)")       # Talke 2018 Table 4A: EMAX = 3090 ADC units (bootstrap 3160, 95% CI 1540-7777)
    lc50     <- log(0.233); label("Effect-compartment concentration giving 50% EMAX, C50 (ug/L)") # Talke 2018 Table 4A: C50 = 0.233 ug/L (bootstrap 0.211, 95% CI 0.098-0.390)
    lhill    <- log(2.83);  label("Hill coefficient N (unitless)")                          # Talke 2018 Table 4A: N = 2.83 (bootstrap 3.28, 95% CI 1.70-4.35); no IIV reported
    lt12keo  <- log(2.16);  label("Effect-compartment equilibration half-life t1/2 keo (min)") # Talke 2018 Table 4A: t1/2 keo = 2.16 min (bootstrap 2.01, 95% CI 1.37-2.90); keo = log(2)/t1/2 keo

    # Allometric scaling exponents - fixed a priori per Talke 2018 Methods (theoretical allometric)
    e_wt_cl <- fixed(0.75); label("Body-weight allometric exponent shared by CL, Q2, Q3")  # Talke 2018 Methods: PWR = 0.75 for clearance terms; fixed (no SE reported)
    e_wt_vc <- fixed(1);    label("Body-weight allometric exponent shared by V1, V2, V3")  # Talke 2018 Methods: PWR = 1 for distribution volumes; fixed (no SE reported)

    # PK IIV - full BLOCK(4) on log(CL), log(V1), log(Q2), log(V2).
    # PPV reported in Talke 2018 Table 3A as sqrt(NONMEM Omega) on the log scale -> variance = PPV^2.
    # Off-diagonal correlations from Talke 2018 Table 3B; covariance = rho * sigma_i * sigma_j.
    # Lower-triangular row order: var(CL), cov(CL,V1), var(V1), cov(CL,Q2), cov(V1,Q2), var(Q2),
    # cov(CL,V2), cov(V1,V2), cov(Q2,V2), var(V2). PPV CL=0.158, V1=0.339, Q2=0.590, V2=0.868.
    # Correlations: rho(CL,V1)=0.420, rho(CL,Q2)=0.685, rho(V1,Q2)=0.367,
    # rho(CL,V2)=0.971, rho(V1,V2)=0.210, rho(Q2,V2)=0.833.
    # NOTE: the Table 3B correlations encoded with the rounded Table 3A PPVs
    # are not jointly positive-definite (min eigenvalue -2.1e-3, a rounding
    # artifact), making the IIV covariance invalid for simulation (chol()
    # fails). The block is adjusted to the nearest positive-definite matrix
    # (Higham 2002; Matrix::nearPD). Changes are confined to the 3rd-4th
    # significant figure except var(CL) (0.024964 -> 0.026892), which rises
    # slightly to admit the strong CL correlations.
    etalcl + etalvc + etalq + etalvp ~ c(0.026892177,
                                         0.022190997, 0.114969390,
                                         0.064167228, 0.073373298, 0.348148950,
                                         0.132676120, 0.061870386, 0.426517130, 0.753547660)

    # PD IIV - full BLOCK(4) on log(E0), log(EMAX), log(C50), log(t1/2 keo).
    # PPV from Talke 2018 Table 4A; correlations from Table 4B. Methods explicitly notes
    # "The covariance of EMAX and C50 variability was incorporated into the model"; remaining
    # covariances reported in Table 4B are encoded as a full block (see vignette Assumptions).
    # Lower-triangular row order: var(E0), cov(E0,EMAX), var(EMAX), cov(E0,C50),
    # cov(EMAX,C50), var(C50), cov(E0,t12), cov(EMAX,t12), cov(C50,t12), var(t12).
    # PPV E0=0.462, EMAX=1.349, C50=1.628, t12=0.841.
    # Correlations: rho(E0,EMAX)=0.116, rho(E0,C50)=0.374, rho(EMAX,C50)=0.897,
    # rho(E0,t12)=0.458, rho(EMAX,t12)=0.716, rho(C50,t12)=0.913.
    etale0 + etalemax + etalc50 + etalt12keo ~ c(0.213444,
                                                 0.072296, 1.819801,
                                                 0.281299, 1.969641, 2.650384,
                                                 0.177952, 0.812308, 1.250032, 0.707281)

    # PK residual error - combined additive + proportional on plasma concentration (ug/L)
    addSd  <- 0.011; label("Additive residual error on Cc (ug/L)")                       # Talke 2018 Table 3A ERRADD = 0.011 ug/L (bootstrap 0.010, 95% CI 0.001-0.132); see vignette for note on eta-on-RUV
    propSd <- 0.104; label("Proportional residual error on Cc (fraction)")               # Talke 2018 Table 3A ERRPROP = 10.4% (bootstrap 11.8%, 95% CI 7.02-19.36)

    # PD residual error - combined additive + proportional on ADC units
    addSd_vasocon  <- 113;   label("Additive residual error on ADC vasoconstriction (ADC units)") # Talke 2018 Table 4A ERRADD = 113 ADC units (bootstrap 103, 95% CI 72-138); see vignette for note on eta-on-RUV
    propSd_vasocon <- 0.618; label("Proportional residual error on ADC vasoconstriction (fraction)") # Talke 2018 Table 4A ERRPROP = 61.8% (bootstrap 61.8%, 95% CI 36.4-326)
  })
  model({
    # Individual PK parameters with a priori allometric weight scaling (reference 70 kg)
    cl  <- exp(lcl  + etalcl) * (WT / 70)^e_wt_cl                  # Talke 2018 Methods: CL  scales as (WT/70)^0.75
    vc  <- exp(lvc  + etalvc) * (WT / 70)^e_wt_vc                  # Talke 2018 Methods: V1  scales as (WT/70)^1
    q   <- exp(lq   + etalq)  * (WT / 70)^e_wt_cl                  # Talke 2018 Methods: Q2  scales as (WT/70)^0.75
    vp  <- exp(lvp  + etalvp) * (WT / 70)^e_wt_vc                  # Talke 2018 Methods: V2  scales as (WT/70)^1
    q2  <- exp(lq2)           * (WT / 70)^e_wt_cl                  # Talke 2018 Methods: Q3  scales as (WT/70)^0.75 (no IIV per Table 3A)
    vp2 <- exp(lvp2)          * (WT / 70)^e_wt_vc                  # Talke 2018 Methods: V3  scales as (WT/70)^1    (no IIV per Table 3A)

    # Individual PD parameters
    e0      <- exp(le0     + etale0)                               # baseline ADC effect (ADC units)
    emax    <- exp(lemax   + etalemax)                             # maximum vasoconstriction effect (ADC units)
    c50     <- exp(lc50    + etalc50)                              # half-maximal effect-site concentration (ug/L)
    hill    <- exp(lhill)                                          # Hill coefficient (no IIV per Table 4A)
    t12keo  <- exp(lt12keo + etalt12keo)                           # t1/2 keo (min)
    keo     <- log(2) / t12keo                                     # Talke 2018 Methods: keo derived from t1/2 keo

    # Micro-rate constants for explicit ODE form
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # 3-compartment IV disposition (dosing into central as IV infusion)
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-                                                      k13 * central - k31 * peripheral2

    # Effect compartment driven by plasma concentration; Ce equilibrates to Cc with rate keo
    Cc <- central / vc                                             # plasma concentration (ug/L)
    d/dt(effect) <- keo * (Cc - effect)                            # Ce in the effect compartment (ug/L)

    # Sigmoid Emax PD: ADC vasoconstriction response (Talke 2018 Methods Eq. for Effect)
    vasocon <- e0 + emax * effect^hill / (c50^hill + effect^hill)

    # Observation models
    Cc      ~ add(addSd)         + prop(propSd)                    # plasma concentration residual error
    vasocon ~ add(addSd_vasocon) + prop(propSd_vasocon)            # ADC vasoconstriction residual error
  })
}
