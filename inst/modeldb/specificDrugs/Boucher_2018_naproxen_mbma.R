Boucher_2018_naproxen_mbma <- function() {
  description <- "MBMA. Model-based meta-analysis longitudinal time-course Emax model for the Western Ontario and McMaster Universities (WOMAC) pain score (0-10 scale) in adults with osteoarthritis, fitted to study-arm-mean data from 18 randomized double-blind placebo-controlled trials of naproxen vs placebo (12 flare designs, 6 non-flare). The WOMAC pain response over time follows a three-parameter Emax model in time: pain = E0 + Emax * time / (ET50 + time), where ET50 is the time to half-maximal effect. Flare design shifts both baseline E0 and Emax; naproxen treatment shifts Emax and shortens ET50 (faster onset: ET50 0.21 week vs placebo 0.69 week). Between-study variability is carried as study-arm-level random effects on E0 (SD 0.62) and Emax (SD 0.74); the residual describes study-arm-mean variability weighted by each arm's observed standard error (sigma fixed to 1). Suitable simulation scope is study-arm-mean WOMAC pain time-course, NOT individual-patient pain scores. Parameter values are the NONMEM column of Table 2 (the same model was fit in NONMEM, BUGS, and R with closely agreeing estimates)."

  reference <- paste(
    "Boucher M, Bennetts M.",
    "The Many Flavors of Model-Based Meta-Analysis:",
    "Part II: Modeling Summary Level Longitudinal Responses.",
    "CPT Pharmacometrics Syst Pharmacol. 2018 May;7(5):288-297.",
    "doi:10.1002/psp4.12299.",
    sep = " "
  )
  vignette <- "Boucher_2018_naproxen_mbma"
  units <- list(
    time          = "week (time since start of treatment; the Emax is a function of time, not of dose or concentration)",
    dosing        = "not applicable (MBMA summary-level time-course model; treatment is encoded by the NAPROXEN study-arm indicator covariate, not by rxode2 dose events)",
    concentration = "WOMAC pain units / arm (0-10 normalised Western Ontario and McMaster Universities osteoarthritis index pain subscale; output Cc is the study-arm mean WOMAC pain score, NOT a drug concentration; the slash is only to satisfy checkModelConventions unit parsing)"
  )

  covariateData <- list(
    FLARE = list(
      description        = "Study-arm flare-design indicator: 1 if the trial used a flare design (subjects washed out of pain medication and required a predefined pain flare-up to be randomized), 0 for a non-flare design.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-flare design)",
      notes              = "MBMA study-arm-level covariate (a property of the trial design, not of an individual patient). Documented inline rather than in the individual-level pop-PK canonical register inst/references/covariate-columns.md, following the Vargo_2014_statins_ezetimibe_mbma and Sadouki_2025 in-file-documentation precedent for MBMA / multi-drug study-arm covariates. Encoded as the indicator 'If' in Boucher 2018 Eqs 2 and 3: flare shifts both baseline E0 (e_flare_e0) and Emax (e_flare_emax). Of the 18 trials, 12 were flare designs and 6 were non-flare.",
      source_name        = "If (Boucher 2018 Eqs 2-3)"
    ),
    NAPROXEN = list(
      description        = "Study-arm treatment indicator: 1 if the arm received naproxen, 0 if the arm received placebo.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (placebo arm)",
      notes              = "MBMA study-arm-level covariate (a property of the trial arm). Documented inline per the Vargo_2014_statins_ezetimibe_mbma / Sadouki_2025 precedent rather than in the individual-level register. Encoded as the indicator 'In' in Boucher 2018 Eqs 3 and 4: naproxen shifts Emax (e_naproxen_emax) and shortens ET50 (e_naproxen_et50, additive on the log scale). All 18 included trials had both a naproxen and a placebo arm.",
      source_name        = "In (Boucher 2018 Eqs 3-4)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 18L,
    disease_state  = "adults with osteoarthritis (knee/hip) and OA pain; endpoint is the WOMAC pain subscale on a 0-10 scale",
    design         = "18 randomized double-blind placebo-controlled parallel-group trials, each with both a naproxen and a placebo arm; 12 trials used a flare design (washout plus required pain flare-up for eligibility) and 6 did not. Source: internal clinical study reports plus publicly available literature.",
    treatments     = "naproxen vs placebo (per-arm naproxen dose was not modeled; the model characterizes the time-course of response pooled across the naproxen doses studied)",
    timepoints     = "WOMAC pain reported over a time-course up to about 13 weeks; the number of trials reporting at weeks 2, 6, and 12 was 13, 9, and 7 respectively (Boucher 2018 Results)",
    notes          = "MBMA at the study-arm level: each modeled data point is the mean WOMAC pain in one trial arm at one timepoint, weighted by its observed standard error (variance SD^2/n). The model is intended for simulating study-arm-mean WOMAC pain time-courses and is NOT suitable for individual-subject simulation. Total patient count is reported only in Supplementary Table S2, which is not on disk; n_subjects is left NA. See Boucher 2018 'Example dataset' section and Figure 1 for the design and time-course."
  )

  ini({
    # ============================================================
    # Three-parameter Emax-in-time model (Boucher 2018 Eq 1):
    #   WOMAC_ijk = E0 + eta1 + (Emax + eta2) * t / (ET50 + t) + eps
    # with E0, Emax, ET50 built from design (flare) and treatment
    # (naproxen) covariates via Eqs 2-4. All values are the NONMEM
    # column of Table 2 (the paper produced diagnostics from the
    # NONMEM fit; BUGS and R(NLME) agreed closely).
    # ============================================================

    # ----- Baseline WOMAC pain E0 (Eq 2): E0 = E0nf + If * dE0f -----
    e0 <- 5.20
    label("Baseline WOMAC pain, non-flare reference (E0 nonflare; 0-10 units)")  # Boucher 2018 Table 2 NONMEM, E0 (nonflare)

    e_flare_e0 <- 0.96
    label("Flare-design additive shift on baseline WOMAC pain (delta E0 flare; 0-10 units)")  # Boucher 2018 Table 2 NONMEM, dE0 (flare)

    # ----- Maximal change Emax (Eq 3, final additive model):
    #   Emax = Emax_p,nf + If * dEmax_pf + In * dEmax_n
    # The flare-by-treatment interaction on Emax (Eq 5) was tested and
    # not significant, so Table 2/Table 3 estimates exclude it.
    emax <- -1.16
    label("Maximal change from baseline, placebo non-flare reference (Emax_p nonflare; 0-10 units, negative = pain reduction)")  # Boucher 2018 Table 2 NONMEM, Emax_p (nonflare)

    e_flare_emax <- -0.82
    label("Flare-design additive shift on Emax (delta Emax_p flare; 0-10 units)")  # Boucher 2018 Table 2 NONMEM, dEmax_p (flare)

    e_naproxen_emax <- -0.79
    label("Naproxen additive shift on Emax (delta Emax naproxen; 0-10 units)")  # Boucher 2018 Table 2 NONMEM, dEmax_n

    # ----- ET50 time to half-maximal effect (Eq 4, log scale):
    #   ln(ET50) = ln(ET50p) + In * ln(dET50n)
    # ET50 fitted in log space to keep it positive.
    let50 <- -0.37
    label("Log time to half-maximal effect, placebo reference (Ln ET50p; ET50 = exp value = 0.69 week)")  # Boucher 2018 Table 2 NONMEM, Ln(ET50p); exp(-0.37) = 0.69 week

    e_naproxen_et50 <- -1.17
    label("Naproxen additive shift on log ET50 (Ln dET50n; naproxen ET50 = exp(-0.37-1.17) = 0.21 week)")  # Boucher 2018 Table 2 NONMEM, Ln(dET50n); naproxen ET50 = exp(let50 + this) = 0.21 week (faster onset than placebo)

    # ============================================================
    # Between-study (study-arm-level) random effects on E0 and Emax
    # (Boucher 2018 Eq 1: eta1 on E0, eta2 on Emax; both ~ N(0, tau^2)).
    # Encoded as MBMA study-level etas (NOT individual between-subject
    # variability) per the SKILL Phase-1 Step-3a MBMA guidance. The
    # ini() value is the VARIANCE = tau^2; Table 2 reports tau (the SD).
    # ============================================================
    eta_study_e0   ~ 0.3844   # Boucher 2018 Table 2 NONMEM, s1 = 0.62 (SD of eta1 on E0); variance = 0.62^2 = 0.3844
    eta_study_emax ~ 0.5476   # Boucher 2018 Table 2 NONMEM, s2 = 0.74 (SD of eta2 on Emax); variance = 0.74^2 = 0.5476

    # ============================================================
    # Residual error (Boucher 2018 Eq 1 / Methods): residuals ~ N(0,
    # SD_ijk^2 / n_ijk). Because the weights were the observed standard
    # errors of the study-arm means, sigma was FIXED to 1. The operative
    # per-observation residual SD is therefore sigma * SE_ijk = SE_ijk,
    # the observed SD_ijk/sqrt(n_ijk) supplied with the data. This file
    # exposes the fixed sigma and leaves per-arm SE reweighting to
    # downstream simulation code (same pattern as Vargo_2014 MBMA).
    # ============================================================
    addSd <- fixed(1)
    label("Residual error multiplier sigma (FIXED to 1); per study-arm-mean observation the residual SD is sigma * SE_ijk = SE_ijk (observed SD/sqrt(n))")  # Boucher 2018 Methods/Eq 1: "residuals ... variance SD_ijk^2/n_ijk ... sigma was fixed to 1"
  })

  model({
    # Study-arm covariates (binary): FLARE (1 = flare design), NAPROXEN
    # (1 = naproxen arm, 0 = placebo).

    # Baseline WOMAC pain E0 (Eq 2): non-flare reference + flare shift.
    # eta1 (between-study) enters additively on E0 per Eq 1.
    baseline <- e0 + e_flare_e0 * FLARE + eta_study_e0

    # Maximal change Emax (Eq 3, final additive model): placebo-non-flare
    # reference + flare shift + naproxen shift. eta2 (between-study)
    # enters additively on Emax per Eq 1.
    emax_arm <- emax + e_flare_emax * FLARE + e_naproxen_emax * NAPROXEN + eta_study_emax

    # ET50 (Eq 4, log scale): placebo reference shifted by naproxen.
    et50 <- exp(let50 + e_naproxen_et50 * NAPROXEN)

    # WOMAC pain time-course (Eq 1 three-parameter Emax in time).
    # At time 0 the fraction is 0 so Cc = baseline; as time grows large
    # Cc approaches baseline + emax_arm. Cc is the study-arm-mean WOMAC
    # pain score (0-10), NOT a drug concentration.
    Cc <- baseline + emax_arm * time / (et50 + time)
    Cc ~ add(addSd)
  })
}
