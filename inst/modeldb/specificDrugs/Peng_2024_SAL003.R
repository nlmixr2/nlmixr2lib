Peng_2024_SAL003 <- function() {
  description <- "Two-compartment population PK model for SAL003, a novel anti-PCSK9 IgG4 monoclonal antibody, with first-order SC absorption (with lag time), saturable Michaelis-Menten elimination from the central compartment, and a body-weight effect on central volume, in Chinese healthy volunteers and patients with hyperlipidemia (Peng 2024)"
  reference <- "Peng J, Huang J, Tan H, Kuang Y, Yang G, Huang Z. Model-Informed Dose Selection for a Novel Human Immunoglobulin G4 Derived Monoclonal Antibody Targeting Proprotein Convertase Kwashiorkor Type 9: Insights from Population Pharmacokinetics-Pharmacodynamics and Systems Pharmacology. ACS Pharmacol Transl Sci. 2024;7(2):406-420. doi:10.1021/acsptsci.3c00256"
  vignette <- "Peng_2024_SAL003"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate (WT/70)^0.77 on central volume of distribution Vc (Peng 2024 Table 2: dVdWeight = 0.77, RSE 20.74%). The paper does not specify the reference body weight associated with the tvV estimate of 5.07 L; 70 kg is used here as the conventional allometric reference and the choice is documented in the vignette's Assumptions and deviations section. The studied population's median body weight is approximately 65 kg.",
      source_name        = "WT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 40L,
    n_studies       = 2L,
    age_range       = "19-65 years",
    weight_range    = "45.7-88.8 kg",
    weight_median   = "approximately 65 kg (61.3 kg in SAD, 67.65 kg in MAD)",
    sex_female_pct  = 38,
    race_ethnicity  = c(HanChinese = 95, OtherChinese = 5),
    disease_state   = "Pooled single-ascending-dose (SAD) phase 1 in healthy Chinese volunteers (n = 20 on SAL003) and multiple-ascending-dose (MAD) phase 1 in Chinese patients with primary hypercholesterolemia or mixed hyperlipidemia (n = 20 on SAL003).",
    dose_range      = "SAD: single subcutaneous 70, 140, or 420 mg (the 35 and 280 mg SAD cohorts were not measured for PK and were excluded from the popPK model). MAD: subcutaneous 140 mg every 4 weeks (Q4W) or 420 mg every 8 weeks (Q8W) over 16 weeks, on a stable atorvastatin background.",
    regions         = "China (Third Xiangya Hospital, Central South University).",
    notes           = "Pooled n = 40 SAL003-treated subjects across SAD and MAD phase 1 studies registered at chinadrugtrials.org.cn (SAD: CTR20200225; MAD: CTR20212013). Baseline demographics from Peng 2024 Table 1; SAL003 PK characterised via ELISA on plasma samples."
  )

  ini({
    # ---- Structural PK parameters (Peng 2024 Table 2, final popPK model) ----
    # Paper's native parameterization uses micro-rate constants K12 and K21 in 1/h;
    # we re-parameterize to Q (L/day) and Vp (L) to follow the nlmixr2lib convention
    # of macro-rate constants (lvc / lq / lvp) for structural log-parameters. The
    # numerical conversions are:
    #   Q  = K12 * Vc = 0.0018 /h * 5.07 L = 0.009126 L/h = 0.21902 L/day
    #   Vp = Q / K21  = 0.009126 / 0.0055 = 1.659 L
    # which match the supplement (Table S3) reported Vper = 1.68 L (rounding).
    # Rates reported in the paper are converted to day-units (x 24) for nlmixr2lib.
    lka     <- log(0.02 * 24)
    label("First-order SC absorption rate Ka (1/day)")                                                       # Peng 2024 Table 2: Ka = 0.02 1/h (RSE 6.50%)

    lvc     <- log(5.07)
    label("Central volume of distribution Vc at 70 kg reference body weight (L)")                            # Peng 2024 Table 2: V = 5.07 L (RSE 5.15%)

    lq      <- log(0.0018 * 5.07 * 24)
    label("Inter-compartmental clearance Q (L/day; derived as K12 * Vc * 24)")                               # Peng 2024 Table 2: K12 = 0.0018 1/h (RSE 22.67%); derived Q = 0.21902 L/day

    lvp     <- log(0.0018 * 5.07 / 0.0055)
    label("Peripheral volume of distribution Vp (L; derived as Q/K21)")                                      # Peng 2024 Table 2: K21 = 0.0055 1/h (RSE 13.33%); derived Vp = 1.659 L (paper supplement Table S3: 1.68 L)

    lvmax   <- log(0.50 * 24)
    label("Maximum Michaelis-Menten elimination rate Vmax (mg/day)")                                         # Peng 2024 Table 2: Vmax = 0.50 mg/h (RSE 5.72%)

    lkm     <- log(9883.96 / 1000)
    label("Michaelis-Menten constant Km (mg/L; reported in paper as 9883.96 ng/mL)")                         # Peng 2024 Table 2: Km = 9883.96 ng/mL (RSE 10.16%)

    ltlag   <- log(1.76 / 24)
    label("SC absorption lag time Tlag (day)")                                                                # Peng 2024 Table 2: Tlag = 1.76 h (RSE 13.19%)

    lfdepot <- fixed(log(0.783))
    label("Subcutaneous bioavailability F (fraction; FIXED at preclinical Macaca fascicularis value 78.3%)") # Peng 2024 Table 2 footnote a: F fixed at preclinical 78.3% from Macaca fascicularis; supplement Table S3: 0.783

    # ---- Covariate effect on Vc ----
    # Peng 2024 Table 2: dVdWeight = 0.77 (RSE 20.74%). The paper does not state the
    # reference body weight; 70 kg is used as the conventional allometric reference.
    # The form is power-law: Vc_i = tvVc * (WT / 70)^e_wt_vc * exp(etalvc).
    e_wt_vc <- 0.77
    label("Power exponent of (WT/70) on central volume of distribution Vc (unitless)")                       # Peng 2024 Table 2: dVdWeight = 0.77 (RSE 20.74%)

    # ---- Inter-individual variability (Peng 2024 Table 2; interpreted as omega^2) ----
    # The paper reports OMEGA values directly for V, Ka, K12, K21, Vmax, Km, and Tlag.
    # Re-parameterizing from (Vc, K12, K21) to (Vc, Q, Vp) makes the original etas
    # linear combinations of the new ones:
    #   eta_Q  = eta_K12 + eta_Vc                  (since log Q  = log K12 + log Vc)
    #   eta_Vp = eta_K12 + eta_Vc - eta_K21        (since log Vp = log Q  - log K21)
    # so the implied joint distribution of (etalvc, etalq, etalvp) is a 3x3 block
    # with off-diagonal covariances equal to the shared eta variances:
    #   var(etalvc) = 0.14
    #   var(etalq)  = 1.57 + 0.14            = 1.71
    #   var(etalvp) = 1.57 + 0.14 + 0.70     = 2.41
    #   cov(etalvc, etalq)  = var(eta_Vc)              = 0.14
    #   cov(etalvc, etalvp) = var(eta_Vc)              = 0.14
    #   cov(etalq,  etalvp) = var(eta_K12) + var(eta_Vc) = 1.71
    # This 3x3 block exactly preserves the joint distribution of individual model
    # predictions under the paper's K12/K21 parameterisation.
    etalvc + etalq + etalvp ~ c(
      0.14,
      0.14, 1.71,
      0.14, 1.71, 2.41
    )
    etalka   ~ 0.26                                                                                          # Peng 2024 Table 2: omega^2(Ka)   = 0.26
    etalvmax ~ 0.25                                                                                          # Peng 2024 Table 2: omega^2(Vmax) = 0.25
    etalkm   ~ 0.65                                                                                          # Peng 2024 Table 2: omega^2(Km)   = 0.65
    etaltlag ~ 0.83                                                                                          # Peng 2024 Table 2: omega^2(Tlag) = 0.83 (shrinkage 28.84%)

    # ---- Residual error (Peng 2024 Table 2) ----
    # 'sigma_proportional PK = 0.12' (RSE 1.49%, shrinkage 5.89%); bootstrap (Table S2)
    # mean = 0.12, 95% CI [0.10, 0.14]. Encoded as a proportional SD on the linear scale.
    propSd <- 0.12
    label("Proportional residual error (SD, fraction)")
  })

  model({
    # ---- Individual PK parameters ----
    ka     <- exp(lka + etalka)
    vc     <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q      <- exp(lq  + etalq)
    vp     <- exp(lvp + etalvp)
    vmax   <- exp(lvmax + etalvmax)
    km     <- exp(lkm + etalkm)
    fdepot <- exp(lfdepot)
    tlag   <- exp(ltlag + etaltlag)

    # ---- Micro-rate constants (derived from macro parameters) ----
    k12 <- q / vc
    k21 <- q / vp

    # ---- Observation ----
    # Dose is in mg, vc in L, so Cc = central / vc has units mg/L (= 1000 ng/mL).
    Cc <- central / vc

    # ---- 2-compartment ODEs with Michaelis-Menten elimination from central ----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - vmax * Cc / (km + Cc) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ---- Bioavailability and absorption lag on SC depot ----
    f(depot)    <- fdepot
    alag(depot) <- tlag

    Cc ~ prop(propSd)
  })
}
