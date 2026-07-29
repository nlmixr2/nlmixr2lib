Canini_2015_danoprevir <- function() {
  description <- paste0(
    "Combined pharmacokinetic / viral-kinetic (PK/VK) model for the HCV NS3/4A protease inhibitor danoprevir in ",
    "adults with chronic hepatitis C (Canini 2015). PK is a two-compartment disposition (central = plasma, ",
    "peripheral1 = tissue including liver) with zero-order absorption of duration Tk0 into central after a ",
    "post-dose lag Tlag, followed by first-order elimination. The viral-kinetic layer follows the standard ",
    "Neumann 1998 model with target cells assumed to remain constant (only the productively-infected cell ",
    "pool `infected` and the free-virion pool `virus` are carried as ODE states). Drug effectiveness ",
    "e(t) = C1 / (EC50 + C1) is a Hill-Emax function with Emax fixed at 1 and Hill coefficient fixed at 1 ",
    "(the Canini 2015 model-selection winner over h = 2 and h = 3 by AICc/BIC); e(t) blocks virion ",
    "production so d/dt(virus) = (1 - e) * p * infected - c * virus. Because Canini 2015 does not estimate ",
    "the virion production rate p or the infection rate beta separately, the packaged model fixes p = 1 as ",
    "a mathematical scaling (the identifiability constraint p * beta_T = d * c at the pre-treatment steady ",
    "state fully determines beta_T = d * c and I(0) = c * V0). All 10 estimated typical values from ",
    "Canini 2015 Table 2 (Tlag, Tk0, V1, ke, k12, k21, V0, c, d, EC50) and all 10 IIVs are carried; ",
    "cohort-specific V0 offsets (Table 2 cohorts 2-5) and the screened HCV genotype covariate are documented ",
    "in vignette Errata but not applied as covariates in the packaged model (see notes in population and in ",
    "the vignette Assumptions section). PK is in per-hour units to match the paper's Table 2 reporting; ",
    "the viral-kinetic /day-scale rate constants c and d are converted inline to /h so the integrated ",
    "system runs on a single hour-scale time axis."
  )

  reference <- paste(
    "Canini L, Chatterjee A, Guedj J, Lemenuel-Diot A, Brennan B, Smith PF, Perelson AS.",
    "(2015). A pharmacokinetic/viral kinetic model to evaluate the treatment effectiveness",
    "of danoprevir against chronic HCV. Antiviral Therapy 20(5):469-477.",
    "doi:10.3851/IMP2879.",
    sep = " "
  )

  vignette <- "Canini_2015_danoprevir"

  units <- list(
    time          = "hr",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    HCV_GT1B = list(
      description        = "HCV genotype-1 subtype indicator. 1 = patient infected with HCV genotype 1B; 0 = patient infected with HCV genotype 1A. Screened in Canini 2015 as a candidate covariate on PK/VK parameters (Methods, Parameter estimation and statistical methods) but not retained in the final PK/VK model reported in Table 2.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HCV GT1A; 30 percent of the Canini 2015 cohort, 12 of 40 patients; 55 percent GT1B, 15 percent GT1 with subtype undetermined).",
      notes              = "Documented in covariatesDataExcluded because the paper explicitly reports that HCV genotype was tested as a covariate ('HCV genotype was tested as a covariate in the model to study its effect on the PK/VK parameters', Methods) but no genotype effect appears in Canini 2015 Table 2. Not referenced in model()."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32L,
    n_studies      = 1L,
    n_cohorts      = 5L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "Adults with chronic HCV infection (30 percent genotype-1a, 55 percent genotype-1b, 15 percent genotype-1 with undetermined subtype; 32 treatment-naive across cohorts 1-4 and 8 prior peginterferon-alpha / ribavirin non-responders in cohort 5, of whom 5-6 non-responders received active drug versus placebo per the 8:2 randomisation).",
    dose_range     = "Danoprevir oral 100 mg BID (cohort 1), 100 mg TID (cohort 2), 200 mg BID (cohort 3), 200 mg TID (cohort 4), and 300 mg BID (cohort 5); each cohort dosed for 14 days of monotherapy.",
    dose_regimens  = "100 mg twice daily (cohort 1), 100 mg three times daily (cohort 2), 200 mg twice daily (cohort 3), 200 mg three times daily (cohort 4), 300 mg twice daily (cohort 5). All oral. Cohorts 1-4 were treatment-naive; cohort 5 comprised prior peginterferon-alpha / ribavirin non-responders.",
    regions        = NA_character_,
    notes          = paste0(
      "Data are the 32 active-drug patients from a randomised (8:2 active:placebo per cohort) Phase 1 single-ascending-dose study; placebo patients were excluded from the modelling dataset. Post day-14 the same patients were allowed to start peginterferon-alpha / ribavirin; only days 0-14 of danoprevir monotherapy contribute to the fit. Cohort-specific typical baseline viral loads V0 reported in Canini 2015 Table 2 (5.4e5, 1.9e6, 1.8e6, 2.6e6, 2.7e6 IU/mL for cohorts 1-5 respectively) are subsumed by the paper's very large 118 percent IIV on V0 in this packaged model; the cohort-1 reference value 5.4e5 IU/mL is used as the pooled typical value. See the validation vignette Errata section for the abstract-vs-Methods contradiction over whether the effectiveness driver concentration is plasma C1 or tissue C2 (the packaged model follows the explicit Methods statement C1)."
    )
  )

  ini({
    # =========================================================================
    # PK typical values (Canini 2015 Table 2, cohort-1 reference). All PK rate
    # constants are in /h and volumes in L; these units match Table 2 as
    # reported. IIVs are reported as %CV in Canini 2015 Table 2 and are
    # translated to the internal log-normal variance omega^2 = log(1 + CV^2).
    # =========================================================================
    ltlag <- log(0.50)
    label("Post-dose absorption lag Tlag (h)")                 # Canini 2015 Table 2: Tlag = 0.50 h (IIV 63.5% CV, se 8.4)

    ld1   <- log(0.58)
    label("Zero-order absorption duration Tk0 (h)")            # Canini 2015 Table 2: Tk0  = 0.58 h (IIV 122% CV, se 17)

    lvc   <- log(2100)
    label("Effective central-compartment volume V1/F1 (L)")    # Canini 2015 Table 2: V1  = 2100 L (IIV 90% CV, se 11); V1 redefined as V1/F1 per Methods (F1 not identifiable)

    lke   <- log(1.02)
    label("First-order central-compartment elimination rate ke (1/h)")  # Canini 2015 Table 2: ke  = 1.02 /h (IIV 30.3% CV, se 6.8)

    lk12  <- log(0.25)
    label("Central-to-peripheral distribution rate k12 (1/h)")  # Canini 2015 Table 2: k12 = 0.25 /h (IIV 79.4% CV, se 25)

    lk21  <- log(0.81)
    label("Peripheral-to-central distribution rate k21 (1/h)")  # Canini 2015 Table 2: k21 = 0.81 /h (IIV 89.2% CV, se 16)

    # =========================================================================
    # Viral-kinetic typical values (Canini 2015 Table 2). c and d are reported
    # in /day; they are converted to /h inside model() by dividing by 24 so
    # the whole system runs on the hour-scale time axis used by PK.
    # =========================================================================
    lc    <- log(5.28)
    label("Virion clearance rate constant c (1/day; converted to 1/h in model())")  # Canini 2015 Table 2: c = 5.28 /day (IIV 38.4% CV, se 5.8)

    ld    <- log(0.18)
    label("Infected-cell loss rate constant d (1/day; converted to 1/h in model())")  # Canini 2015 Table 2: d = 0.18 /day (IIV 55.9% CV, se 13)

    lrbase <- log(5.4e5)
    label("Baseline (pre-treatment) viral load V0 (IU/mL)")     # Canini 2015 Table 2: V0 = 5.4e5 IU/mL (cohort-1 reference; IIV 118% CV, se 14). Cohort-specific V0 offsets for cohorts 2-5 are subsumed by the 118% IIV in this packaged model (see population notes and vignette Errata).

    lec50 <- log(0.0082)
    label("Danoprevir EC50 for virion-production inhibition (ng/mL)")  # Canini 2015 Table 2: EC50 = 0.0082 ng/mL (IIV 152% CV, se 22); the Hill coefficient h was fixed to 1 after AICc/BIC comparison of h = 1 (best), 2, 3.

    # =========================================================================
    # Fixed structural constants (not estimated in Canini 2015).
    #   h (Hill exponent): fixed to 1 after model selection over 1, 2, 3.
    #   p (virion production rate) and beta * T (infection rate scaled by
    #     target-cell concentration): NOT estimated in the paper. At the pre-
    #     treatment steady state with V > 0, the identifiability constraint
    #     p * beta_T = d * c fixes the product but leaves one degree of freedom.
    #     Following Neumann-family convention the packaged model fixes p = 1
    #     as a mathematical scaling; then I(0) = c * V0 / p = c * V0 and
    #     beta_T = d * c / p = d * c. p = 1 is not a biological rate; users
    #     interested in the biological I on hepatocyte-count scale would need
    #     to rescale I and beta_T by their choice of p.
    # =========================================================================
    lhill <- fixed(log(1))
    label("Hill exponent for the Emax effectiveness function (unitless; FIXED at 1 in Canini 2015 after AICc/BIC selection)")  # Canini 2015 Results, PK/VK model: 'Fitting the PK/VK model for h=1 (AICc=6097), h=2 (AICc=6123) and h=3 (AICc=6154), showed that h=1 provided the best model'.

    lp    <- fixed(log(1))
    label("Virion production rate p (virions per infected cell per unit time; FIXED at 1 as a Neumann-family scaling because Canini 2015 does not identify p separately)")

    # =========================================================================
    # IIV (Canini 2015 Table 2). All %CV values reported in the 'IIV, % +- se'
    # row are converted to log-normal variance via omega^2 = log(1 + CV^2).
    # =========================================================================
    etaltlag ~ log(1 + 0.635^2)  # Canini 2015 Table 2: Tlag IIV = 63.5%  CV -> log(1 + 0.635^2) = 0.339
    etald1   ~ log(1 + 1.22^2)   # Canini 2015 Table 2: Tk0  IIV = 122%   CV -> log(1 + 1.22^2)  = 0.912
    etalvc   ~ log(1 + 0.90^2)   # Canini 2015 Table 2: V1   IIV = 90%    CV -> log(1 + 0.90^2)  = 0.593
    etalke   ~ log(1 + 0.303^2)  # Canini 2015 Table 2: ke   IIV = 30.3%  CV -> log(1 + 0.303^2) = 0.0878
    etalk12  ~ log(1 + 0.794^2)  # Canini 2015 Table 2: k12  IIV = 79.4%  CV -> log(1 + 0.794^2) = 0.489
    etalk21  ~ log(1 + 0.892^2)  # Canini 2015 Table 2: k21  IIV = 89.2%  CV -> log(1 + 0.892^2) = 0.585
    etalc    ~ log(1 + 0.384^2)  # Canini 2015 Table 2: c    IIV = 38.4%  CV -> log(1 + 0.384^2) = 0.138
    etald    ~ log(1 + 0.559^2)  # Canini 2015 Table 2: d    IIV = 55.9%  CV -> log(1 + 0.559^2) = 0.272
    etalrbase ~ log(1 + 1.18^2)  # Canini 2015 Table 2: V0   IIV = 118%   CV -> log(1 + 1.18^2)  = 0.872
    etalec50 ~ log(1 + 1.52^2)   # Canini 2015 Table 2: EC50 IIV = 152%   CV -> log(1 + 1.52^2)  = 1.197

    # =========================================================================
    # Residual error (Canini 2015 Results, PK/VK model). PK residual is a
    # combined additive + proportional error model on the linear-scale
    # danoprevir plasma concentration (ng/mL). VK residual is additive on
    # log10(viral load).
    # =========================================================================
    addSd    <- 0.25
    label("PK additive residual SD on plasma danoprevir concentration (ng/mL)")  # Canini 2015 Results, PK/VK model: 'additive error term (a = 0.25 +- 0.012 ng/ml)'

    propSd   <- 0.61
    label("PK proportional residual SD on plasma danoprevir concentration (fraction)")  # Canini 2015 Results, PK/VK model: 'proportional error term (b = 0.61 +- 0.02)'

    addSd_Vlog10 <- 0.29
    label("Viral-load additive residual SD on log10(HCV RNA IU/mL)")  # Canini 2015 Results, PK/VK model: 'An additive error model was found to best describe the residual error for the VK data (a = 0.29 +- 0.0094 log10 IU/ml)'.
  })

  model({
    # =========================================================================
    # 1. Individual PK parameters (all in /h units matching the paper's
    # Table 2 reporting; see units$time = 'hr').
    # =========================================================================
    tlag <- exp(ltlag + etaltlag)
    d1   <- exp(ld1   + etald1)
    vc   <- exp(lvc   + etalvc)
    ke   <- exp(lke   + etalke)
    k12i <- exp(lk12  + etalk12)
    k21i <- exp(lk21  + etalk21)

    # =========================================================================
    # 2. Individual viral-kinetic rate constants. c and d are stored in /day
    # (Table 2 reporting); divide by 24 to convert to /h to match the PK
    # time axis. V0, EC50, Hill, and p are already dimensionless or scale-
    # neutral with respect to the time axis.
    # =========================================================================
    c_indiv    <- exp(lc + etalc) / 24
    d_indiv    <- exp(ld + etald) / 24
    v0_indiv   <- exp(lrbase + etalrbase)
    ec50_indiv <- exp(lec50 + etalec50)
    hill_indiv <- exp(lhill)
    p_indiv    <- exp(lp)

    # =========================================================================
    # 3. PK ODE system: two-compartment disposition with zero-order absorption
    # of duration d1 (= Tk0) into central after a post-dose lag tlag; first-
    # order elimination at rate ke from central and reversible distribution
    # to peripheral1 with forward rate k12i and backward rate k21i.
    #
    # Event-table consumers should administer each dose as:
    #   evid = 1, cmt = 'central', amt = <mg>, rate = -2
    # so rxode2 uses dur(central) = d1 to compute the constant zero-order
    # infusion rate D / d1 over duration d1 after the lag.
    # =========================================================================
    d/dt(central)     <- -ke * central - k12i * central + k21i * peripheral1
    d/dt(peripheral1) <-  k12i * central - k21i * peripheral1

    alag(central) <- tlag
    dur(central)  <- d1

    # =========================================================================
    # 4. Plasma danoprevir concentration in ng/mL. central is in mg (matches
    # the dose unit) and vc is in L, so central / vc has units mg/L = ug/mL;
    # multiplying by 1000 converts to ng/mL to match Table 2 EC50 and the
    # paper's Cmax / Cmin reporting scale.
    # =========================================================================
    cc_ngml <- (central / vc) * 1000

    # =========================================================================
    # 5. Time-varying effectiveness (Canini 2015 pharmacodynamics section:
    # 'e(t) varies as a function of C1, according to the Emax model'). Emax
    # is fixed at 1. Hill exponent is fixed at 1 by AICc/BIC selection
    # (h = 1 vs h = 2 vs h = 3). The 1e-30 numerical floor on the
    # denominator prevents division-by-zero when both cc_ngml and ec50 are
    # at the lower numerical limit early in the simulation (before any
    # dosing, cc_ngml = 0 and ec50 > 0 so the ratio is well-defined; the
    # floor is defensive).
    # =========================================================================
    e_t <- cc_ngml^hill_indiv /
           (ec50_indiv^hill_indiv + cc_ngml^hill_indiv + 1e-30)

    # =========================================================================
    # 6. Viral-kinetic ODE system (Neumann 1998 standard model with target
    # cells assumed constant per Canini 2015 Methods). Because Canini 2015
    # does not identify p and beta separately, p is fixed at 1 and beta * T
    # is set to d * c to satisfy the pre-treatment steady-state constraint
    # p * beta_T = d * c. With p = 1 the derived quantities are:
    #   beta_T = d * c   (infection-rate * target-cell concentration)
    #   I(0)   = c * V0  (infected-cell steady-state initial condition)
    #   V(0)   = V0      (pre-treatment viral load)
    # d/dt(infected) = beta_T * virus - d * infected
    #                = d * c * virus - d * infected
    # d/dt(virus)    = (1 - e_t) * p * infected - c * virus
    #                = (1 - e_t)     * infected - c * virus  (p = 1)
    # =========================================================================
    beta_t <- d_indiv * c_indiv / p_indiv

    d/dt(infected) <- beta_t * virus - d_indiv * infected
    d/dt(virus)    <- (1 - e_t) * p_indiv * infected - c_indiv * virus

    # Initial conditions at pre-treatment steady state (V > 0):
    infected(0) <- c_indiv * v0_indiv / p_indiv
    virus(0)    <- v0_indiv

    # =========================================================================
    # 7. Observation outputs and residual error. Cc is the plasma danoprevir
    # concentration in ng/mL (paper's reporting unit); Vlog10 is log10 of
    # the free-virion concentration (log10(IU/mL)). The 1e-12 floor inside
    # log10() prevents -Inf when virus approaches zero after sustained
    # effective treatment.
    # =========================================================================
    Cc     <- cc_ngml
    Vlog10 <- log10(virus + 1e-12)

    Cc     ~ add(addSd) + prop(propSd)
    Vlog10 ~ add(addSd_Vlog10)
  })
}
