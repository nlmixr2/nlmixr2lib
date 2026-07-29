Cao_2017_dha_ring_early <- function() {
  description <- paste(
    "In vitro (P. falciparum 3D7 laboratory strain, early-ring-stage parasites).",
    "Dynamic stress PD model from Cao 2017 capturing the delayed dihydroartemisinin (DHA)",
    "killing effect observed in tightly age-synchronized parasite cultures; one of four",
    "stage-specific NLME fits in Table 1 (early-ring stage corresponds to 2 h post-infection",
    "per the Klonis 2013 experimental design that supplied the viability data). The killing",
    "rate k = kmax(S) * C^hill / (Kc(S)^hill + C^hill) is modulated by a dynamic stress",
    "variable S(t) that accumulates while drug concentration C exceeds C* = 0.1 nM",
    "(dS/dt = lambda*(1 - S)) and resets toward zero otherwise. Stress-dependent",
    "kmax(S) = alpha*S and Kc(S) = beta1*(1 - S) + beta2 (paper eq 7 and 8). DHA",
    "concentration evolves in the central compartment with first-order decay (default",
    "kdrug = log(2)/8 /h for in vitro experiments; override for in vivo simulations).",
    "Parasite count N(t) is normalized to N(0) = 1, so the deterministic `parasites`",
    "state IS the surviving viability fraction (paper eq 17). See modellib('Cao_2017_dha_ring_mid'),",
    "modellib('Cao_2017_dha_troph_early'), modellib('Cao_2017_dha_troph_late') for the",
    "other three stage-specific fits combined in the in vivo PK-PD simulation of Fig 6."
  )
  reference <- paste(
    "Cao P, Klonis N, Zaloumis S, Dogovski C, Xie SC, Saralamba S, White LJ, Fowkes FJI,",
    "Tilley L, Simpson JA, McCaw JM. (2017). A dynamic stress model explains the delayed",
    "drug effect in artemisinin treatment of Plasmodium falciparum.",
    "Antimicrob Agents Chemother 61(12):e00618-17.",
    "doi:10.1128/AAC.00618-17."
  )
  vignette <- "Cao_2017_dha_artemisinin_stress_model"
  paper_specific_compartments <- c("stress", "parasites")
  units <- list(
    time          = "h",
    dosing        = "nM (initial DHA concentration deposited into central)",
    concentration = "nM (central DHA); unitless fraction (parasites = viability)"
  )

  covariateData <- list()

  population <- list(
    species        = "in vitro (P. falciparum 3D7 laboratory strain)",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    disease_state  = paste(
      "Tightly age-synchronized 3D7 parasites (>80% within a 1-h age window) at the",
      "early ring stage (2 h post-infection per Klonis et al 2013, Cao 2017 reference 7).",
      "Cultures grown in human red blood cells; viability assayed in the trophozoite stage",
      "of the following life cycle, 48 h after pulse start, with M3 censoring at the",
      "0.005 limit of detection (Materials and methods)."
    ),
    dose_range     = paste(
      "Initial DHA pulse concentrations approximately 39 nM and 300 nM in Cao 2017 Fig 1,",
      "with additional concentrations in the supplemental data of Klonis 2013 (Cao 2017",
      "reference 7); pulse durations 1, 2, 4, and 6 h."
    ),
    notes          = paste(
      "The model was fit by NLME (NONMEM 7.3.0 + Perl-speaks-NONMEM 3.7.6) separately to",
      "viability data for each parasite life-cycle stage. The Hill coefficient `hill`",
      "(paper symbol gamma) was fixed across stages at 1.7892, the pooled mean of the",
      "stage-specific estimates in Table S1 (Table 1 footnote a). Residual error was",
      "partitioned into between- and within-duplicate normal components (variances",
      "sigma_b^2 and sigma_w^2) but the numeric variances are NOT reported in Cao 2017",
      "Table 1. Per nlmixr2lib policy on unreported variance components for mechanistic",
      "models, the file is encoded as a typical-value deterministic mechanism without",
      "residual error or IIV; see vignette Errata."
    )
  )

  ini({
    # =====================================================================
    # Stage-specific NLME estimates -- EARLY RING stage (Cao 2017 Table 1).
    # All wrapped in fixed() because they are published point estimates
    # reproduced for typical-value simulation; nlmixr2lib does not refit.
    # =====================================================================

    lambda <- fixed(6.2504)
    label("Stress accumulation rate lambda (1/h)")
    # Cao 2017 Table 1, Early ring stage: lambda = 6.2504 /h
    # (SE 0.5745; model-based 95% CI 5.1243-7.3765). The half-life of the
    # unstressed state ln(2)/lambda = 0.11 h is the parenthetical column 2
    # entry in Table 1.

    alpha <- fixed(1.6915)
    label("Maximum achievable killing rate alpha (1/h)")
    # Cao 2017 Table 1, Early ring stage: alpha = 1.6915 /h
    # (SE 0.1378; model-based 95% CI 1.4215-1.9616). alpha is the maximum
    # killing rate achieved at S = 1 and saturating drug concentration
    # (paper eq 7: kmax(S) = alpha * S).

    beta1 <- fixed(990.84)
    label("Stress-dependent Kc slope beta1 (nM)")
    # Cao 2017 Table 1, Early ring stage: beta1 = 990.84 nM
    # (SE 373.49; model-based 95% CI 258.81-1722.9). beta1 governs the
    # stress-modulation of Kc per paper eq 8: Kc(S) = beta1*(1-S) + beta2,
    # so Kc(S=0) = beta1 + beta2 and Kc(S=1) = beta2. The parametric
    # bootstrap CI for early ring (-20377, 1876.5) includes negative values
    # but the asymptotic point estimate is positive and is used here.

    beta2 <- fixed(12.519)
    label("Stationary half-maximal killing concentration beta2 (nM)")
    # Cao 2017 Table 1, Early ring stage: beta2 = 12.519 nM
    # (SE 1.0631; model-based 95% CI 10.435-14.602). beta2 = Kc(S=1)
    # (paper eq 8), the half-maximal killing concentration once stress
    # has saturated.

    hill <- fixed(1.7892)
    label("Hill coefficient (paper symbol gamma; unitless; fixed across stages)")
    # Cao 2017 Table 1 footnote a: "gamma is fixed to be 1.7892 based on
    # estimates in Table S1 in the supplemental material" -- the stage-
    # pooled mean of the four individually-estimated gamma values.

    Cstar <- fixed(0.1)
    label("Drug-concentration threshold for stress accumulation C* (nM)")
    # Cao 2017 Fig 6 legend: "The modulatory variable S is assumed to
    # follow equation 6 only when DHA concentration, C, is >= 0.1 nM
    # (i.e., C* = 0.1 nM), and S is immediately reset to zero when the
    # DHA concentration drops below 0.1 nM."

    kdrug_pulse <- fixed(log(2) / 8)
    label("DHA elimination rate during a drug pulse (1/h); default in vitro t1/2 = 8 h")
    # Cao 2017 page 5: "the in vitro half-life of DHA (t1/2) was measured
    # to be about 8 h (21)"; kdrug_pulse = log(2)/8 ~ 0.08664 /h drives
    # the exponential decay of central DHA during a pulse (paper eq 13:
    # C(t) = C0 * 0.5^(t / t1/2)). For in vivo simulation override to
    # log(2)/0.9 ~ 0.770 /h via rxSolve(params = c(kdrug_pulse=log(2)/0.9))
    # to match the 0.9 h plasma DHA half-life cited from references 2
    # and 21 of Cao 2017 (eq 18, in vivo PK).

    kreset <- fixed(100)
    label("Fast first-order rate constant for wash-out and stress reset (1/h)")
    # Encoded as a fast first-order decay (rate 100 /h, effective half-
    # life ~25 s) to approximate (a) the paper's instantaneous in vitro
    # drug-wash at the end of a pulse and (b) the instantaneous stress
    # reset specification: "S is immediately reset to zero when the DHA
    # concentration drops below 0.1 nM" (Fig 6 legend). The ODE-solver-
    # friendly continuous approximation is mathematically near-equivalent
    # on the timescales of parasite killing (1-6 h pulses; 24 h dosing
    # intervals in vivo).

    tend_pulse <- fixed(6)
    label("End time of in vitro DHA pulse (h); user-overridable for non-default pulse durations or in vivo simulation")
    # Klonis 2013 (Cao 2017 reference 7) tested pulse durations of 1, 2,
    # 4, and 6 h; the default 6 h is the longest. Override via
    # rxSolve(params = c(tend_pulse = <value>)) per simulation. For in
    # vivo simulation where there is no in vitro wash (the drug profile
    # is determined by the repeated dosing events and the elimination
    # rate kdrug_pulse), set tend_pulse to a value larger than the entire
    # simulation horizon (e.g., 1e6) so the wash-out branch is never
    # triggered and central decays purely at kdrug_pulse between doses.
  })

  model({
    # =====================================================================
    # 1. Central-compartment DHA PK (paper eq 13 in vitro form).
    #    During the pulse (time < tend_pulse) DHA decays at kdrug_pulse
    #    (default in vitro log(2)/8 /h, t1/2 = 8 h). For time >=
    #    tend_pulse the elimination rate jumps to kreset (~100 /h, t1/2 ~
    #    25 s) to approximate the paper's instantaneous in vitro wash-out.
    #    For in vivo simulation set tend_pulse to a value larger than the
    #    simulation horizon and override kdrug_pulse to log(2)/0.9 /h so
    #    central decays purely at the in vivo half-life between bolus
    #    doses. Initial DHA concentration C0 is set via an amt event into
    #    the central compartment at time = 0.
    # =====================================================================
    in_pulse <- (time < tend_pulse)
    # Active decay = either during the pulse (slow kdrug_pulse) OR post-
    # pulse with central still measurable (fast kreset). The (central >
    # 1e-9) floor stops the post-pulse decay once central is numerically
    # indistinguishable from zero, preventing subnormal-arithmetic NaN
    # propagation during long post-wash observation horizons.
    kdrug <- in_pulse * kdrug_pulse +
             (1 - in_pulse) * kreset * (central > 1e-9)
    d/dt(central) <- -kdrug * central

    # =====================================================================
    # 2. Dynamic stress (paper eq 6 + Fig 6 legend reset specification).
    #    Stress accumulates with first-order kinetics toward 1 while drug
    #    concentration exceeds C*, and decays back to 0 with rate kreset
    #    when drug falls below C*. The (1 - above_thresh) factor pivots
    #    the branch via a 0/1 indicator evaluated at solver time.
    # =====================================================================
    above_thresh <- (central > Cstar)
    # Stress reset only kicks in while drug is below C* AND stress is
    # itself still > 1e-9 (the same numerical-noise floor used on central).
    # Stops the decay once stress is essentially zero so subnormal
    # arithmetic does not contaminate downstream Hill / parasite terms.
    d/dt(stress) <- above_thresh * (lambda * (1 - stress)) +
                    (1 - above_thresh) * (-kreset * stress) *
                      (stress > 1e-9)
    stress(0) <- 0

    # =====================================================================
    # 3. Stress-modulated kmax and Kc (paper eq 7 and 8).
    # =====================================================================
    kmax <- alpha * stress
    Kc   <- beta1 * (1 - stress) + beta2

    # =====================================================================
    # 4. Killing rate (paper eq 5). Small epsilon terms keep the Hill
    #    expression numerically well-behaved near central = 0 (during and
    #    after washout) and near Kc = 0 (mid-ring stage where beta2 is
    #    near zero and stress is high).
    # =====================================================================
    # Clip central to non-negative before raising to hill power; numerical
    # noise in the ODE solver can briefly drive central slightly negative
    # post-wash, and (negative)^hill produces NaN. The clip uses the
    # 0/1 indicator (central > 0) to zero-out the contribution exactly
    # when the solver is in the post-wash near-zero regime.
    Ceps   <- (central > 0) * central + 1e-30
    Kceps  <- Kc + 1e-30
    k_kill <- kmax * Ceps^hill / (Kceps^hill + Ceps^hill)

    # =====================================================================
    # 5. Parasite count (paper eq 4 for tightly age-synchronized parasites:
    #    dN/dt = -k(C, S) * N). N(0) is normalized to 1, so the surviving
    #    `parasites` state IS the viability fraction (paper eq 17).
    # =====================================================================
    d/dt(parasites) <- -k_kill * parasites
    parasites(0) <- 1
  })
}
