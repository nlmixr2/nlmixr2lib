VasquezBahena_2009_lumiracoxib_rat <- function() {
  description <- "Preclinical (rat). Two-compartment population PK plus indirect-response PK/PD model for the antinociceptive effect of oral lumiracoxib in carrageenan-induced thermal hyperalgesia in female Wistar rats (Vasquez-Bahena 2009). PK: first-order absorption with lag time and dose-dependent relative bioavailability. PD: time-variant (gamma function) carrageenan-induced COX-2 synthesis with first-order COX-2 degradation; lumiracoxib reversibly inactivates COX-2 via a competitive binding model (COX-2_act = KD * COX-2 / (KD + Cp)). The level of inflammatory mediators (MED) equals the active COX-2 amount and drives the paw withdrawal latency response LT = LT0 / (1 + MED)."
  reference <- "Vasquez-Bahena DA, Salazar-Morales UE, Ortiz MI, Castaneda-Hernandez G, Troconiz IF. Pharmacokinetic-pharmacodynamic modelling of the analgesic effects of lumiracoxib, a selective inhibitor of cyclooxygenase-2, in rats. Br J Pharmacol. 2010;159(1):176-187. doi:10.1111/j.1476-5381.2009.00508.x"
  vignette <- "VasquezBahena_2009_lumiracoxib_rat"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    DOSE = list(
      description        = "Per-subject assigned oral lumiracoxib dose level (mg/kg). Enters the dose-dependent relative-bioavailability formula Frel = 1 - IMAX * DOSE / (D50 + DOSE).",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reported per-subject mg/kg dose level (Methods, Study design: 1, 3, 10 or 30 mg/kg in experiment I; 10 or 30 mg/kg in experiment II). The mg amount delivered into the depot at the dose record (`amt`) equals DOSE * body weight in kg; the mg/kg covariate is required separately because the D50 = 4.3 estimate in Table 1 is on the mg/kg scale.",
      source_name        = "DOSE"
    ),
    CARRAGEENAN = list(
      description        = "Binary indicator for intraplantar carrageenan injection (1) vs saline injection (0). Switches the COX-2 synthesis-rate model from a constant rate (saline) to a gamma-function time-variant rate (carrageenan).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (saline injection / no inflammatory stimulus)",
      notes              = "Methods, Study design: animals in groups II-IX received a single intraplantar injection of 1% carrageenan suspension (100 uL) into the right hind paw at the experiment start; group I (and group VII in experiment II for the pre-vehicle window) received saline. CARRAGEENAN = 0 selects ks_cox2 = ks_cox2_saline (Table 3 footnote a); CARRAGEENAN = 1 selects the time-variant gamma function ks_cox2(t) = A * t^alpha * exp(-beta * t) (paper Eq. 3, Table 3 'Experiment I and II').",
      source_name        = "CARRAGEENAN"
    )
  )

  population <- list(
    species         = "rat (Wistar, female, fasted)",
    n_subjects      = 80L,
    n_studies       = 2L,
    age_range       = "Not reported",
    weight_range    = "0.180-0.200 kg",
    sex_female_pct  = 100,
    disease_state   = "Carrageenan-induced thermal hyperalgesia in the right hind paw (1% carrageenan suspension, 100 uL intraplantar) assessed by the Hargreaves plantar test",
    dose_range      = "0, 1, 3, 10 or 30 mg/kg oral lumiracoxib (single dose; 0.5% carboxymethylcellulose / Tween 80 suspension, 4 mL/kg)",
    regions         = "Mexico (Centro de Investigacion y de Estudios Avanzados del Instituto Politecnico Nacional, Mexico City)",
    notes           = "Methods, Animals + Study design. Experiment I: 60 female Wistar rats randomly allocated to six groups (I = saline; II = carrageenan only; III-VI = carrageenan + lumiracoxib 1/3/10/30 mg/kg co-administered with carrageenan). Experiment II: 20 female rats randomly allocated to three groups (VII = carrageenan + vehicle; VIII = carrageenan + 10 mg/kg lumiracoxib at 4 h post-carrageenan; IX = carrageenan + 30 mg/kg lumiracoxib at 4 h post-carrageenan). PK sampling was performed only in experiment I (13 nominal time points between 0.083 and 10 h after lumiracoxib). Latency-time PD measurements were collected in both experiments (9 time points in experiment I covering 1-10 h, 16 time points in experiment II covering 1-24 h after carrageenan)."
  )

  ini({
    # ------------------------------------------------------------------
    # PK structural parameters -- final population estimates from Table 1
    # (p. 181, "Parameter estimates from the final population pharmacokinetic
    # model"; fit to experiment I plasma data).
    # ------------------------------------------------------------------
    lka    <- log(15.2)  ; label("First-order absorption rate constant (1/h)")            # Table 1: ka = 15.2 (RSE 23%)
    lcl    <- log(0.16)  ; label("Plasma clearance (L/h)")                                # Table 1: CL = 0.16 (RSE 14%)
    lvc    <- log(0.49)  ; label("Central volume of distribution (L)")                    # Table 1: Vc = 0.49 (RSE 9%)
    lq     <- log(0.42)  ; label("Inter-compartmental clearance (L/h)")                   # Table 1: Q = 0.42 (RSE 12%)
    lvp    <- log(0.99)  ; label("Peripheral volume of distribution (L)")                 # Table 1: VT = 0.99 (RSE 19%)
    ltlag  <- log(0.06)  ; label("Absorption lag time (h)")                               # Table 1: Tlag = 0.06 (RSE 10%)

    # Dose-dependent relative bioavailability Frel = exp(lfrel + etalfrel) * (1 - IMAX * DOSE / (D50 + DOSE)).
    # Final estimates from Table 1. lfrel is fixed at log(1) and anchors the
    # exponential IIV on Frel (Methods, Data analysis: "Inter-individual
    # variability (IIV) was modelled exponentially and expressed as coefficient
    # of variation"); the dose-dependent reduction is applied as a derived
    # multiplicative factor (Ide 2009 pravastatin pattern in this package).
    lfrel      <- fixed(log(1)) ; label("Relative-bioavailability anchor F0 (fixed at 1; the dose-dependent reduction is applied separately)")  # paired with etalfrel; Table 1 IIV row
    limax_frel <- log(0.67)     ; label("Maximum fractional reduction in relative bioavailability IMAX (dimensionless)") # Table 1: IMAX = 0.67 (RSE 6%)
    ld50_frel  <- log(4.3)      ; label("Dose eliciting half-maximal reduction in Frel D50 (mg/kg)")                     # Table 1: D50 = 4.3 mg/kg (RSE 43%)

    # ------------------------------------------------------------------
    # PD structural parameters -- final population estimates from Table 3
    # (p. 183, "Experiment I and II" columns, fit to pooled latency-time data).
    # ------------------------------------------------------------------
    # Constant COX-2 synthesis rate for the saline-only group (Table 3 footnote a).
    lks_cox2_saline <- log(0.21) ; label("Constant COX-2 synthesis rate in saline-injected rats (COX-2/h)")  # Table 3 footnote a: ks_COX-2 (saline) = 0.21 (RSE 38%)

    # Gamma-function carrageenan-induced COX-2 synthesis: ks(t) = A * t^alpha * exp(-beta * t).
    la_cox2     <- log(1.42) ; label("Gamma-function amplitude A for carrageenan-induced COX-2 synthesis (COX-2/h^(1+alpha))") # Table 3: A = 1.42 (RSE 15%)
    lalpha_cox2 <- log(1.3)  ; label("Gamma-function shape exponent alpha (dimensionless)")                                    # Table 3: alpha = 1.3 (RSE 21%)
    lbeta_cox2  <- log(0.33) ; label("Gamma-function decay-rate constant beta (1/h)")                                          # Table 3: beta = 0.33 (RSE 17%)

    lkd_cox2 <- log(0.89) ; label("COX-2 degradation rate constant kD_COX-2 (1/h)")                                            # Table 3: kD_COX-2 = 0.89 (RSE 39%)

    # Drug-receptor binding and baseline latency.
    lkd_drug <- log(0.24) ; label("Equilibrium dissociation constant of the COX-2-lumiracoxib complex KD (ug/mL)")             # Table 3: KD = 0.24 (RSE 21%)
    llt0     <- log(14.3) ; label("Paw withdrawal latency time at baseline LT0 (s)")                                           # Table 3: LT0 = 14.3 (RSE 3%)

    # ------------------------------------------------------------------
    # Inter-individual variability -- exponential IIV expressed as CV%
    # (Methods, Data analysis: "Inter-individual variability (IIV) was
    # modelled exponentially and expressed as coefficient of variation
    # [CV (%)]"). Internal log-scale variance: omega^2 = log(CV^2 + 1).
    # ------------------------------------------------------------------
    etalcl   ~ log(0.36^2 + 1)  # Table 1: IIV(CL) = 36%
    etalvc   ~ log(0.39^2 + 1)  # Table 1: IIV(Vc) = 39%
    etalka   ~ log(0.92^2 + 1)  # Table 1: IIV(kA) = 92%
    etalfrel ~ log(0.18^2 + 1)  # Table 1: IIV on Frel (whole formula row) = 18%

    etalkd_drug ~ log(0.62^2 + 1)  # Table 3: IIV(KD) experiment I+II = 62%
    etallt0     ~ log(0.13^2 + 1)  # Table 3: IIV(LT0) experiment I+II = 13%

    # ------------------------------------------------------------------
    # Residual error -- Table 1 (PK, combined) and Table 3 (PD, additive).
    # ------------------------------------------------------------------
    addSd  <- 0.037 ; label("Additive residual error on lumiracoxib plasma concentration (ug/mL)") # Table 1: Additive = 0.037 (RSE 20%)
    propSd <- 0.22  ; label("Proportional residual error on lumiracoxib plasma concentration (fraction)") # Table 1: Proportional = 22% (RSE 10%)

    addSd_LT <- 2.4 ; label("Additive residual error on paw withdrawal latency response (s)")     # Table 3: Residual error = 2.4 s (RSE 10%, experiment I+II)
  })

  model({
    # ------------------------------------------------------------------
    # Individual PK parameters
    # ------------------------------------------------------------------
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)
    tlag <- exp(ltlag)

    # Dose-dependent relative bioavailability with exponential IIV anchored
    # by lfrel (fixed at log(1)) and a derived dose-dependent reduction.
    imax_frel <- exp(limax_frel)
    d50_frel  <- exp(ld50_frel)
    frel      <- exp(lfrel + etalfrel) * (1 - imax_frel * DOSE / (d50_frel + DOSE))

    # ------------------------------------------------------------------
    # Individual PD parameters
    # ------------------------------------------------------------------
    ks_cox2_saline <- exp(lks_cox2_saline)
    a_cox2         <- exp(la_cox2)
    alpha_cox2     <- exp(lalpha_cox2)
    beta_cox2      <- exp(lbeta_cox2)
    kd_cox2        <- exp(lkd_cox2)
    kd_drug        <- exp(lkd_drug + etalkd_drug)
    lt0            <- exp(llt0 + etallt0)

    # Time-variant COX-2 synthesis rate (Eq. 3). t is the global simulation
    # clock measured from the carrageenan injection at t = 0. The gamma
    # term evaluates to 0 at t = 0 (since alpha_cox2 > 0), matching the
    # paper's assumption that COX-2 expression is below detection prior
    # to the inflammatory stimulus (Methods, COX-2 model, paragraph 2).
    ks_cox2_gamma <- a_cox2 * t^alpha_cox2 * exp(-beta_cox2 * t)
    ks_cox2       <- (1 - CARRAGEENAN) * ks_cox2_saline + CARRAGEENAN * ks_cox2_gamma

    # ------------------------------------------------------------------
    # ODE system
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # COX-2 turnover (Eq. 2). Initial condition COX-2(0) = 0 (Methods,
    # COX-2 model, paragraph 2).
    d/dt(cox2) <- ks_cox2 - kd_cox2 * cox2

    f(depot)    <- frel
    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # Observations
    # ------------------------------------------------------------------
    # Lumiracoxib plasma concentration (central amt in mg / vc in L =
    # mg/L = ug/mL).
    Cc <- central / vc

    # Active (inhibitable) COX-2 in the absence vs presence of lumiracoxib
    # via competitive binding (Eq. 4). MED equals COX-2_act (Methods,
    # Pharmacokinetic/pharmacodynamic modelling, final paragraph).
    cox2_act <- (kd_drug * cox2) / (kd_drug + Cc)
    med      <- cox2_act

    # Paw withdrawal latency response (Eq. 1).
    LT <- lt0 / (1 + med)

    Cc ~ add(addSd) + prop(propSd)
    LT ~ add(addSd_LT)
  })
}
