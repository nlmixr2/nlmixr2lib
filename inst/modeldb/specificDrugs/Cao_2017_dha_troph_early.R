Cao_2017_dha_troph_early <- function() {
  description <- paste(
    "In vitro (P. falciparum 3D7 laboratory strain, early-trophozoite-stage parasites).",
    "Dynamic stress PD model from Cao 2017 capturing the delayed dihydroartemisinin (DHA)",
    "killing effect; one of four stage-specific NLME fits in Table 1 (early-trophozoite",
    "stage corresponds to 24 h post-infection). The killing rate",
    "k = kmax(S) * C^hill / (Kc(S)^hill + C^hill) is modulated by a dynamic stress",
    "variable S(t) that accumulates while drug concentration C exceeds C* = 0.1 nM",
    "(dS/dt = lambda*(1 - S)) and resets toward zero otherwise. Stress-dependent",
    "kmax(S) = alpha*S and Kc(S) = beta1*(1 - S) + beta2 (paper eq 7 and 8). Trophozoite",
    "stages have substantially higher maximum killing rate alpha than ring stages",
    "(Fig 4B), reflecting greater drug susceptibility once stress has accumulated. DHA",
    "concentration evolves in the central compartment with first-order decay (default",
    "kdrug = log(2)/8 /h for in vitro experiments; override for in vivo simulations).",
    "Parasite count N(t) is normalized to N(0) = 1, so the deterministic `parasites`",
    "state IS the surviving viability fraction (paper eq 17). See modellib('Cao_2017_dha_ring_early'),",
    "modellib('Cao_2017_dha_ring_mid'), modellib('Cao_2017_dha_troph_late') for the",
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
      "early trophozoite stage (24 h post-infection per Klonis et al 2013, Cao 2017",
      "reference 7). Cultures grown in human red blood cells; viability assayed in the",
      "trophozoite stage of the following life cycle, 48 h after pulse start, with M3",
      "censoring at the 0.005 limit of detection."
    ),
    dose_range     = paste(
      "Initial DHA pulse concentrations approximately 39 nM and 300 nM in Cao 2017 Fig 1,",
      "with additional concentrations in the supplemental data of Klonis 2013 (Cao 2017",
      "reference 7); pulse durations 1, 2, 4, and 6 h."
    ),
    notes          = paste(
      "Trophozoite stages exhibit higher alpha and faster lambda than ring stages",
      "(Fig 4); early trophozoite has the highest alpha across all four stages",
      "(5.7434 /h vs 1.6915, 1.1224, 2.8626 for early-ring, mid-ring, late-trophozoite",
      "respectively). The stage-pooled Hill coefficient gamma was fixed at 1.7892",
      "(Table 1 footnote a). Residual error variances were estimated by NONMEM but",
      "are not numerically reported in Table 1; per nlmixr2lib policy the file is encoded",
      "as a typical-value deterministic mechanism without residual error or IIV; see",
      "vignette Errata."
    )
  )

  ini({
    # =====================================================================
    # Stage-specific NLME estimates -- EARLY TROPHOZOITE stage (Cao 2017 Table 1).
    # =====================================================================

    lambda <- fixed(1.2290)
    label("Stress accumulation rate lambda (1/h)")
    # Cao 2017 Table 1, Early trophozoite stage: lambda = 1.2290 /h
    # (SE 0.2249; model-based 95% CI 0.7882-1.6698). Half-life of the
    # unstressed state ln(2)/lambda = 0.56 h (Table 1 parenthetical).

    alpha <- fixed(5.7434)
    label("Maximum achievable killing rate alpha (1/h)")
    # Cao 2017 Table 1, Early trophozoite stage: alpha = 5.7434 /h
    # (SE 0.7460; model-based 95% CI 4.2813-7.2054). Largest alpha across
    # all four stages, indicating that early-trophozoite parasites are the
    # most susceptible to DHA killing at saturating stress (Fig 4B).

    beta1 <- fixed(317.64)
    label("Stress-dependent Kc slope beta1 (nM)")
    # Cao 2017 Table 1, Early trophozoite stage: beta1 = 317.64 nM
    # (SE 86.143; model-based 95% CI 148.80-486.48).

    beta2 <- fixed(39.570)
    label("Stationary half-maximal killing concentration beta2 (nM)")
    # Cao 2017 Table 1, Early trophozoite stage: beta2 = 39.570 nM
    # (SE 4.6038; model-based 95% CI 30.546-48.593). beta2 = Kc(S = 1)
    # per eq 8; the half-maximal killing concentration once stress has
    # saturated.

    hill <- fixed(1.7892)
    label("Hill coefficient (paper symbol gamma; unitless; fixed across stages)")
    # Cao 2017 Table 1 footnote a: fixed at 1.7892, stage-pooled mean of
    # the four individually-estimated gamma values from Table S1.

    Cstar <- fixed(0.1)
    label("Drug-concentration threshold for stress accumulation C* (nM)")
    # Cao 2017 Fig 6 legend: C* = 0.1 nM; S follows eq 6 only when DHA
    # concentration C >= 0.1 nM; S resets to zero when C < 0.1 nM.

    kdrug_pulse <- fixed(log(2) / 8)
    label("DHA elimination rate during a drug pulse (1/h); default in vitro t1/2 = 8 h")
    # Cao 2017 page 5: in vitro DHA half-life ~8 h (reference 21).
    # Override to log(2)/0.9 ~ 0.770 /h for in vivo simulation (eq 18).

    kreset <- fixed(100)
    label("Fast first-order rate for wash-out and stress reset (1/h)")
    # Fast first-order approximation of the paper's instantaneous reset
    # specification (Fig 6 legend).

    tend_pulse <- fixed(6)
    label("End time of in vitro DHA pulse (h); user-overridable")
    # Klonis 2013 (Cao 2017 reference 7) tested pulses of 1, 2, 4, 6 h.
    # For in vivo simulation set tend_pulse large.
  })

  model({
    # =====================================================================
    # 1. Central-compartment DHA PK (paper eq 13 in vitro form). During
    #    the pulse (time < tend_pulse) DHA decays at kdrug_pulse; for
    #    time >= tend_pulse the elimination rate jumps to kreset to
    #    approximate instantaneous wash-out.
    # =====================================================================
    in_pulse <- (time < tend_pulse)
    # The (central > 1e-9) floor stops the post-pulse decay once central
    # is numerically indistinguishable from zero, preventing subnormal-
    # arithmetic NaN propagation during long post-wash observations.
    kdrug <- in_pulse * kdrug_pulse +
             (1 - in_pulse) * kreset * (central > 1e-9)
    d/dt(central) <- -kdrug * central

    # =====================================================================
    # 2. Dynamic stress (paper eq 6 + Fig 6 legend reset specification).
    # =====================================================================
    above_thresh <- (central > Cstar)
    # Stress reset only kicks in while drug is below C* AND stress is
    # itself still > 1e-9 (same numerical-noise floor used on central).
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
    # 4. Killing rate (paper eq 5).
    # =====================================================================
    # Clip central to non-negative before raising to hill power; numerical
    # noise in the ODE solver can briefly drive central slightly negative
    # post-wash, and (negative)^hill produces NaN.
    Ceps   <- (central > 0) * central + 1e-30
    Kceps  <- Kc + 1e-30
    k_kill <- kmax * Ceps^hill / (Kceps^hill + Ceps^hill)

    # =====================================================================
    # 5. Parasite count (paper eq 4; N(0) = 1 normalization -> viability).
    # =====================================================================
    d/dt(parasites) <- -k_kill * parasites
    parasites(0) <- 1
  })
}
