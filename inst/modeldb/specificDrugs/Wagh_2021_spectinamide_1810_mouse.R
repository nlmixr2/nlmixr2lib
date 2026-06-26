Wagh_2021_spectinamide_1810_mouse <- function() {
  description <- "Preclinical (BALB/c mouse, Mycobacterium tuberculosis infection). Population PK + PK/PD model for subcutaneous spectinamide 1810 in a murine TB efficacy / dose-fractionation study. PK is a two-compartment first-order absorption model with all volumes and clearances expressed per kg body weight (mg/kg dosing, mg/L plasma); IIV is carried on CL/F only (13.2% CV in infected animals). PK/PD couples the plasma central concentration Cc to a hypothetical PAE (post-antibiotic effect) compartment that tracks Cc whenever Cc exceeds the PAE concentration and otherwise decays first-order at rate K_PAE; the PAE concentration drives bacterial killing through a sigmoidal Emax (K_kill_max, EC50, Hill g) on a one-population logistic-growth Mycobacterium tuberculosis model (K_gs net growth rate, log10 N_max carrying capacity, log10 baseline CFU at aerosol infection time). A binary STUDY_WAGH_2 covariate switches K_kill_max from the study 1 typical value to the study 2 value via a 1.15 multiplicative factor; study-specific log10 CFU residual SDs are exposed as parameters."
  reference <- paste(
    "Wagh S, Rathi C, Lukka PB, Parmar K, Temrikar Z, Liu J,",
    "Scherman MS, Lee RE, Robertson GT, Lenaerts AJ, Meibohm B. (2021).",
    "Model-Based Exposure-Response Assessment for Spectinamide 1810 in",
    "a Mouse Model of Tuberculosis.",
    "Antimicrob Agents Chemother 65(11):e01744-20.",
    "doi:10.1128/AAC.01744-20. PMID 34424046.",
    sep = " "
  )
  vignette <- "Wagh_2021_spectinamide_1810_mouse"
  units <- list(
    time          = "hour",
    dosing        = "mg/kg (per-kg dosing; PK volumes and clearances are also per-kg)",
    concentration = "mg/L for plasma Cc and PAE compartment; log10(CFU/lung) for the bacterial PD observation"
  )

  covariateData <- list(
    STUDY_WAGH_2 = list(
      description        = "Binary indicator that the mouse cohort is study 2 of Wagh 2021 (versus the reference study 1 cohort).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (study 1: 7.08 log CFU initial-treatment bacterial load, RUV1 = 0.284 log10 CFU)",
      notes              = "Drives the multiplicative study-2 effect on K_kill_max (1.15-fold higher killing). Wagh 2021 Table 3: study-specific coefficient for K_kill_max = 1.15 (RSE 4.0%). Subject-level (time-fixed); set once from the trial identifier on each subject record. See vignette Assumptions and deviations for the per-study residual-error treatment.",
      source_name        = "STUDY (Wagh 2021 Tables 4-5: study 1 vs study 2)"
    )
  )

  population <- list(
    species        = "mouse (BALB/c, female)",
    n_subjects     = 376L,
    n_studies      = 2L,
    age_range      = "8 weeks at study start (acclimatized 72 h prior to dosing)",
    weight_range   = "18-22 g body weight (Wagh 2021 Methods)",
    sex_female_pct = 100,
    race_ethnicity = NA_character_,
    disease_state  = "Females BALB/c mice chronically infected with Mycobacterium tuberculosis (ATCC 35801 / TMC 107 / Erdman; MIC of spectinamide 1810 = 1.6 mg/L) via low-dose aerosol exposure (~100 CFU/lung deposition); treatment initiated at day 34 post-infection. Healthy uninfected controls (n=84) provided dense-sampling PK that seeded the Bayesian priors in the infected-animal sparse-sampling fit. The current model is the final integrated PK/PD model in the infected animals.",
    dose_range     = "Spectinamide 1810 subcutaneous, 20-800 mg/kg per dose, weekly totals 20-4000 mg/kg, frequencies QW / BIW / TIW / QD / BID. 29 dose-group combinations across studies 1 and 2 (Wagh 2021 Tables 4 and 5).",
    regions        = "USA (Univ. of Tennessee Health Science Center for healthy-animal PK; Colorado State University BSL-3 for infected efficacy / PK)",
    notes          = paste(
      "Total infected mice n = 196 (147 with evaluable PK + log CFU in study 1 +",
      "49 from study 2). 84 healthy mice provided PK only. The integrated final",
      "popPK model used Bayesian priors from the healthy-animal fit (via the",
      "PRIOR NWPRI subroutine in NONMEM) to stabilise estimation in the",
      "sparse-sampling infected-animal data (two PK samples per animal at 0.25",
      "and 8 h post-dose). Bacterial counts (log CFU per lung) from each mouse",
      "were obtained via lung harvest after 4 weeks of treatment + 2-day washout,",
      "with one PD sample per animal (destructive sampling). PK parameters in",
      "the final model are the infected-animal estimates (Wagh 2021 Table 1,",
      "Infected column); healthy-animal estimates are not encoded (the paper",
      "concluded there was no relevant PK difference between healthy and",
      "infected animals -- see vignette Errata)."
    )
  )

  ini({
    # =========================================================================
    # Plasma PK (Wagh 2021 Table 1, Infected column). Time = hours; volumes
    # and clearances are per kg of body weight. Spectinamide 1810 doses are
    # subcutaneous and reach the central compartment via first-order
    # absorption with rate ka. CL/F and V/F are apparent (F-relative); the
    # bioavailability anchor lfdepot is fixed at log(1) because the absolute
    # subcutaneous F is not estimated in the paper.
    # =========================================================================
    lka     <- log(15.4)
    label("Absorption rate constant K_a (1/h)")
    # Wagh 2021 Table 1: K_a (infected) = 15.4 1/h (RSE 8%).

    lvc     <- log(0.217)
    label("Apparent central volume V_c/F (L/kg)")
    # Wagh 2021 Table 1: V_c/F (infected) = 0.217 L/kg (RSE 3%).

    lvp     <- log(0.160)
    label("Apparent peripheral volume V_p/F (L/kg)")
    # Wagh 2021 Table 1: V_p/F (infected) = 0.160 L/kg (RSE 12%).

    lcl     <- log(0.697)
    label("Apparent clearance CL/F (L/h/kg)")
    # Wagh 2021 Table 1: CL/F (infected) = 0.697 L/h/kg (RSE 2%).

    lq      <- log(0.0089)
    label("Apparent intercompartmental clearance Q/F (L/h/kg)")
    # Wagh 2021 Table 1: Q/F (infected) = 0.0089 L/h/kg (RSE 6%).

    lfdepot <- fixed(log(1))
    label("Bioavailability anchor F = 1 (CL/F and V/F are apparent F-relative)")
    # Subcutaneous F is not estimated; PK parameters are apparent (F-relative).

    # =========================================================================
    # Bacterial natural growth (Wagh 2021 Table 2). All FIXED during the
    # PK/PD analysis (paper Methods 'PK/PD modeling with PAE estimation':
    # "The parameters related to the bacterial natural growth model were
    # fixed during the PK/PD analysis to the previously established values").
    # The paper reports N_max and initial CFU on the log10 scale; the model
    # stores their natural-log forms so the bare bmax and rbase used in
    # model() are the linear-scale carrying capacity and initial CFU.
    # =========================================================================
    lkgs    <- fixed(log(0.0357))
    label("Net bacterial growth rate K_gs (1/h; FIXED from natural-growth fit)")
    # Wagh 2021 Table 2: K_gs = 0.0357 1/h (RSE 4%), corresponding to a
    # bacterial doubling time of log(2)/K_gs = 19.4 h (the paper reports a
    # mean generation time of 28.0 h via 1/K_gs * ln(2) when expressed on
    # logistic-growth terms; both forms are equivalent under the logistic
    # parameterisation used here).

    lbmax   <- fixed(log(10) * 6.11)
    label("Log of carrying capacity (natural-log of N_max linear; CFU per lung; FIXED)")
    # Wagh 2021 Table 2: Log N_max = 6.11 (RSE 1%) on log10 scale ->
    # N_max = 10^6.11 = 1.288e6 CFU per lung. Stored as natural log of the
    # linear value so bmax = exp(lbmax) returns the linear carrying capacity.

    lrbase  <- fixed(log(10) * 1.77)
    label("Log of initial bacterial count at aerosol infection time (natural-log of CFU per lung; FIXED)")
    # Wagh 2021 Table 2: Log CFU initial = 1.77 (RSE 3%) on log10 scale ->
    # 10^1.77 = 58.9 CFU per lung. Used as the bacteria(0) initial condition.

    # =========================================================================
    # PK/PD with PAE (Wagh 2021 Table 3). The PAE compartment tracks the
    # central plasma concentration Cc during the rising portion of each
    # dosing interval and decays first-order at rate K_PAE after each peak;
    # its concentration drives the sigmoidal Emax bacterial-kill term.
    # K_kill_max corresponds to the study 1 typical value; the
    # study-specific scalar (e_study_killmax) multiplies K_kill_max for
    # mice in study 2 (STUDY_WAGH_2 = 1).
    # =========================================================================
    lkillmax <- log(0.0374)
    label("Max bacterial kill rate K_kill_max in study 1 (1/h)")
    # Wagh 2021 Table 3: K_kill_max = 0.0374 1/h (RSE 5.3%).

    lec50   <- log(79.6)
    label("EC50 of bacterial kill on PAE concentration (mg/L)")
    # Wagh 2021 Table 3: EC50 = 79.6 mg/L (RSE 60.2%).

    lhill   <- log(1.58)
    label("Hill coefficient g of bacterial kill function (unitless)")
    # Wagh 2021 Table 3: g = 1.58 (RSE 23.1%).

    lkpae   <- log(0.0142)
    label("First-order PAE decay rate K_PAE (1/h)")
    # Wagh 2021 Table 3: K_PAE = 0.0142 1/h (RSE 74.6%); corresponds to a
    # PAE half-life of log(2) / K_PAE = 48.8 h (Wagh 2021 Discussion).

    e_study_killmax <- fixed(1.15)
    label("Study-2 multiplicative effect on K_kill_max relative to study 1 (FIXED)")
    # Wagh 2021 Table 3: study-specific coefficient for K_kill_max = 1.15
    # (RSE 4.0%); applied multiplicatively as K_kill_max_study2 = K_kill_max *
    # (1 + STUDY_WAGH_2 * (e_study_killmax - 1)) = K_kill_max * 1.15 when
    # STUDY_WAGH_2 = 1.

    # =========================================================================
    # Between-animal variability. Wagh 2021 Methods: "Estimation of
    # between-animal variability was limited to clearance." The PD model
    # could not separately estimate between-animal variability because only
    # one CFU measurement was available per animal (paper Methods 'PK/PD
    # modeling with PAE estimation').
    # omega^2 = log(1 + CV^2) for log-normal IIV:
    # omega^2 = log(1 + 0.132^2) = 0.01727.
    # =========================================================================
    etalcl  ~ 0.01727
    # Wagh 2021 Table 1: BAV(CL/F) (infected) = 13.2% CV (RSE 24%, shrinkage 59%).

    # =========================================================================
    # Residual error. The PK residual is reported as 60% CV on the LTBS
    # (log-transform-both-sides) scale (Wagh 2021 Table 1, RUV infected
    # column), which corresponds to a proportional residual on the linear
    # concentration scale in rxode2 / nlmixr2. The PD residuals are
    # additive on the log10 CFU scale (Wagh 2021 Methods 'PK/PD modeling
    # with PAE estimation': "the CFU data were log transformed by taking
    # the decadic logarithm to CFU... RUV was characterized using an
    # additive error"). The model carries the study 1 residual SD as the
    # bare addSd_log_cfu and the study 2 residual SD as a paper-specific
    # alternate parameter; users can swap between the two for stochastic
    # VPC simulation of either cohort.
    # =========================================================================
    propSd  <- 0.600
    label("Proportional residual SD on plasma Cc (LTBS scale, infected animals)")
    # Wagh 2021 Table 1: RUV (infected) = 60.0% CV (RSE 8%, shrinkage 4%).

    addSd_log_cfu <- 0.284
    label("Additive residual SD on log10(CFU/lung) for study 1 mice (study-1 RUV1)")
    # Wagh 2021 Table 3: RUV1 (study 1) = 0.284 (RSE 12.1%).
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Individual PK parameters with IIV on CL only (paper restriction).
    # -----------------------------------------------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    # Fixed bacterial-growth and drug-effect parameters back-transformed.
    kgs      <- exp(lkgs)
    bmax     <- exp(lbmax)
    rbase    <- exp(lrbase)
    killmax_base <- exp(lkillmax)
    ec50     <- exp(lec50)
    hill     <- exp(lhill)
    kpae     <- exp(lkpae)

    # -----------------------------------------------------------------------
    # 2. Study-specific effective K_kill_max. STUDY_WAGH_2 = 0 (study 1
    #    reference) -> killmax = killmax_base; STUDY_WAGH_2 = 1 (study 2) ->
    #    killmax = killmax_base * e_study_killmax = killmax_base * 1.15.
    # -----------------------------------------------------------------------
    killmax <- killmax_base * (1 + STUDY_WAGH_2 * (e_study_killmax - 1))

    # -----------------------------------------------------------------------
    # 3. Plasma central concentration (mg/L).
    # -----------------------------------------------------------------------
    Cc <- central / vc

    # -----------------------------------------------------------------------
    # 4. Two-compartment first-order absorption PK ODE system. Central
    #    derivative is computed explicitly so the PAE compartment can
    #    track d/dt(Cc) during the rising portion of each dosing interval
    #    (see equation block 5 below).
    # -----------------------------------------------------------------------
    dcentral_dt   <- ka * depot - (cl + q) / vc * central + q / vp * peripheral1
    d/dt(depot)      <- -ka * depot
    d/dt(central)    <-  dcentral_dt
    d/dt(peripheral1) <-  q / vc * central - q / vp * peripheral1

    # -----------------------------------------------------------------------
    # 5. PAE compartment (Wagh 2021 Methods 'PK/PD modeling with PAE
    #    estimation' + Figure 5). The PAE concentration C_PAE behaves as
    #    a high-water-mark filter of the plasma Cc with first-order
    #    decay: conceptually C_PAE tracks Cc while Cc is rising and
    #    decays first-order at K_PAE after each in-dose peak (Wagh 2021
    #    Discussion: PAE half-life ~ 48.8 h). Encoded as a fast-load /
    #    slow-decay effect compartment:
    #      d C_PAE / dt = k_load_pae * weight_up * (Cc - C_PAE) - K_PAE * C_PAE
    #    where weight_up is a sigmoidal switch that is approximately 1
    #    when Cc > C_PAE and approximately 0 when Cc < C_PAE
    #    (sigmoidal-step approximation of the paper's piecewise
    #    high-water-mark to avoid the discontinuity that would otherwise
    #    break the ODE solver). k_load_pae is fixed at 100/h, well above
    #    K_a = 15.4/h, so effect tracks rising Cc tightly during
    #    absorption; below 1% lag at the in-dose peak under any
    #    Wagh 2021 dosing regimen. effect(0) = 0 because there is no
    #    drug at the aerosol infection time.
    # -----------------------------------------------------------------------
    k_load_pae <- 100
    weight_up  <- 1 / (1 + exp(-50 * (Cc - effect)))
    d/dt(effect) <- k_load_pae * weight_up * (Cc - effect) - kpae * effect
    Cpae <- effect

    # -----------------------------------------------------------------------
    # 6. Bacterial killing by spectinamide 1810. Sigmoidal Emax-Hill on
    #    the PAE concentration; the kill rate is per-bacterium (1/h).
    # -----------------------------------------------------------------------
    kill_rate <- killmax * Cpae ^ hill / (ec50 ^ hill + Cpae ^ hill)

    # -----------------------------------------------------------------------
    # 7. One-population logistic-growth bacterial dynamics (Wagh 2021
    #    Eq. 4):
    #      d N / dt = K_gs * (1 - N / N_max) * N - kill_rate * N
    #    N is the bacterial count per lung (linear scale). The initial
    #    condition matches the aerosol-infection-time bacterial load
    #    fitted by the natural-growth model (Wagh 2021 Table 2,
    #    log_CFU_initial = 1.77 on log10 scale = 58.9 CFU per lung).
    # -----------------------------------------------------------------------
    d/dt(bacteria) <- kgs * (1 - bacteria / bmax) * bacteria - kill_rate * bacteria

    bacteria(0) <- rbase

    # -----------------------------------------------------------------------
    # 8. Bioavailability anchor (F = 1; CL/F and V/F are apparent).
    # -----------------------------------------------------------------------
    f(depot) <- exp(lfdepot)

    # -----------------------------------------------------------------------
    # 9. Observations.
    #    - Cc: plasma spectinamide 1810 concentration (mg/L), proportional
    #          residual on the LTBS scale (Wagh 2021 Table 1, infected
    #          animals).
    #    - log_cfu: log10 of bacterial count per lung (the paper's PD
    #          observation, "the CFU data were log transformed by taking
    #          the decadic logarithm"). Additive residual on the log10
    #          scale; the bare addSd_log_cfu corresponds to study 1
    #          (RUV1). Study 2 (RUV2 = 0.173) is documented in vignette
    #          Errata; users can override addSd_log_cfu for study-2
    #          stochastic VPCs.
    # -----------------------------------------------------------------------
    log_cfu <- log10(bacteria)
    Cc       ~ prop(propSd)
    log_cfu  ~ add(addSd_log_cfu)
  })
}
