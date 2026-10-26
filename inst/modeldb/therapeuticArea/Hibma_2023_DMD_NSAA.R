Hibma_2023_DMD_NSAA <- function() {
  description <- paste(
    "Joint disease-progression model for the North Star Ambulatory",
    "Assessment (NSAA, total score 0-34) and loss of ambulation (LOA) in",
    "boys with Duchenne muscular dystrophy (DMD), fit to the placebo and",
    "domagrozumab arms of the phase II trial B5161002 pooled with the CINRG",
    "natural-history study B5161011 (Hibma 2023). The NSAA is a modified",
    "indirect-response state integrated over age from birth (NSAA = 0 at",
    "age 0): a zero-order production rate kin and a dissipation rate",
    "kout * age that grows linearly with age, so the score rises to a",
    "peak near age 5.5 years and then declines. LOA is a parametric",
    "time-to-event model on age with a log-normal baseline hazard; the",
    "individual NSAA dissipation rate constant enters the hazard as a",
    "power covariate (kout / median kout)^theta. Correlated exponential",
    "IIV on kin and kout; additive residual error on the NSAA total",
    "score. No covariates were retained and the model carries no drug",
    "effect."
  )
  reference <- paste(
    "Hibma JE, Jayachandran P, Neelakantan S, Harnisch LO.",
    "Disease progression modeling of the North Star Ambulatory Assessment",
    "for Duchenne Muscular Dystrophy.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(3):375-386.",
    "doi:10.1002/psp4.12921"
  )
  vignette <- "Hibma_2023_DMD_NSAA"
  units <- list(
    time = "year (subject age)",
    dosing = "n/a (disease-progression model with no drug input)",
    concentration = "NSAA total score (points, observation nsaa)"
  )

  compartmentData <- list(
    nsaa = list(
      analyte = "North Star Ambulatory Assessment total score",
      units = "points",
      specimen = "not applicable",
      verified = TRUE
    ),
    cumhaz = list(
      analyte = "Cumulative hazard of loss of ambulation",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    RACE_ASIAN = list(
      description = "Asian race indicator: 1 = Asian, 0 = non-Asian.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Tested by stepwise covariate modeling on the joint NSAA-LOA model and not retained (Hibma 2023 Results, 'Covariate analysis'). Race was split Asian vs non-Asian because Asian was the only non-White race contributing more than 10% of patients."
    ),
    CONMED_STEROID = list(
      description = "Corticosteroid use indicator: 1 = on corticosteroids, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no corticosteroid)",
      notes = "Steroid-naive at baseline, 3 months of stable steroid use at baseline, and steroid use at the time of observation were each screened (Appendix S1) and not retained. Fewer than 10% of participants were steroid-naive (Hibma 2023 Discussion)."
    ),
    STUDY_B5161002 = list(
      description = "Study indicator: 1 = domagrozumab phase II trial B5161002, 0 = CINRG natural-history study B5161011.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (B5161011, CINRG natural-history study)",
      notes = "Study type (clinical trial vs natural history) and the B5161002 treatment sequence (domagrozumab yes/no and by sequence) were screened and not retained. Dystrophin mutation status (Exon 44/45/51/53 skip-amenable, deletion 3-7) was screened as an exploratory covariate and not retained; mutation status was unknown for over 70% of participants."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 276L,
    n_studies = 2L,
    n_observations = "2001 NSAA total scores (1404 from 120 participants in B5161002; 597 from 156 participants in B5161011); 29 participants lost ambulation.",
    age_range = "Baseline age 6.0-16.0 years (B5161002) and 4.0-20.3 years (B5161011)",
    age_median = "Baseline age 9.1 years (B5161002) and 7.4 years (B5161011)",
    weight_range = "13.4-108.5 kg (both studies pooled; Table S3)",
    weight_median = "28.9 kg (B5161002) and 23.3 kg (B5161011)",
    sex_female_pct = 0,
    race_ethnicity = "Predominantly White; Asian 26 and 15 participants in the two Table S1 columns (Table S1 column headers appear transposed for some rows -- see vignette).",
    disease_state = "Duchenne muscular dystrophy; ambulatory at baseline (one B5161011 participant non-ambulatory at baseline). Baseline NSAA total score median 21.5 (B5161002) and 23.0 (B5161011).",
    dose_range = "Not applicable -- natural-history disease-progression model. B5161002 participants received domagrozumab (5, 20 or 40 mg/kg) or placebo in a sequence design; no treatment effect was retained.",
    regions = "International (Pfizer phase II trial and the Cooperative International Neuromuscular Research Group natural-history network).",
    notes = "Baseline demographics in Hibma 2023 Tables S1-S4. NSAA values of zero (29 below the LLOQ of 1) were excluded from the NSAA fit and treated as LOA events; NSAA observations above the ULOQ of 33 were handled with the M3 method. B5161002 participants were assumed to be on stable corticosteroids from 6 months before study start."
  )

  ini({
    # Hibma 2023 Table 1 (final joint NSAA-LOA model). NONMEM 7.4.3.
    lkin <- log(8.757)
    label("Log zero-order NSAA production rate constant kin (NSAA points/year)") # Table 1: theta kin = 8.757 NSAA/year (RSE 0.2956%)
    lkout <- log(0.05767)
    label("Log NSAA dissipation rate constant kout, multiplied by age in the ODE (1/year^2)") # Table 1: theta kout = 0.05767 years^-1 (RSE 0.3526%); footnote a: (theta kout x age) is the dissipation rate constant

    # Log-normal LOA baseline hazard on age (Hibma 2023 Equation 8).
    ltmed_loa <- 3.216
    label("Log-normal LOA hazard location mu: log of the median age at LOA for a subject at the median kout (log-years)") # Table 1: theta mu location = 3.216 (RSE 4.271%)
    sdlogt_loa <- 0.3551
    label("Log-normal LOA hazard scale sigma: SD of log age at LOA (unitless)") # Table 1: theta sigma scale = 0.3551 (RSE 13.04%)
    e_kout_haz <- 2.476
    label("Power exponent of (kout / median kout) on the LOA hazard (unitless)") # Table 1: theta koutHAZ = 2.476 (RSE 19.04%); Equations 9-10

    # Table 1 prints IIV as 100*sqrt(omega^2): its 95% CIs are the square
    # roots of a Wald CI on omega^2 (kin: sqrt(0.15801 -/+ 1.96*0.02034*0.15801)
    # = 0.3895, 0.4053, printed (38.95, 40.53)).
    # Covariance = 0.9697 * 0.3975 * 0.9282 = 0.35778.
    etalkin + etalkout ~ c(0.15801, 0.35778, 0.86156) # Table 1: IIV kin 39.75%, IIV kout 92.82%, correlation kin-kout 96.97%

    addSd <- 2.782
    label("Additive residual error on the NSAA total score (points)") # Table 1: theta residual error = 2.782 total score (RSE 0.3585%)
  })

  model({
    kin <- exp(lkin + etalkin)
    kout <- exp(lkout + etalkout)

    # Equation 1: dNSAA/dt = kin - kout * Age * NSAA. The integration
    # variable t is age in years and NSAA starts at zero at birth
    # (Discussion: the initial condition is set at age zero).
    d/dt(nsaa) <- kin - kout * t * nsaa

    # Equations 9-10: h(t) = h0(t) * (kout / median(kout))^theta. With
    # log-normal IIV the population median of kout is its typical value
    # exp(lkout), so the ratio reduces to exp(etalkout).
    cov_kout_haz <- (kout / exp(lkout))^e_kout_haz

    # Equation 8: log-normal hazard on age. h0 is 0 at age 0 (its limit);
    # the branch keeps log(t) out of the t = 0 evaluation.
    if (t > 0) {
      z_loa <- (log(t) - ltmed_loa) / sdlogt_loa
      h0_loa <- exp(-0.5 * z_loa^2) / (t * sdlogt_loa * sqrt(2 * pi)) / (1 - phi(z_loa))
    } else {
      h0_loa <- 0
    }
    hazard <- h0_loa * cov_kout_haz

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0

    # Probability of remaining ambulatory at age t.
    sur <- exp(-cumhaz)

    nsaa ~ add(addSd)
  })
}
