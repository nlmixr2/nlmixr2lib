`Catalan-Latorre_2018_taurine_rat` <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Population PK model for taurine (2-aminoethylsulphonic acid) in male",
    "Wistar rats after IV bolus or oral gavage administration (1, 10, or",
    "100 mg per animal). Two-compartment disposition (central and",
    "peripheral1) with zero-order endogenous formation Q0, first-order",
    "passive oral absorption ka, first-order inter-compartmental",
    "distribution (K12, K21), and non-linear renal elimination described",
    "as two parallel Michaelis-Menten processes: saturable tubular",
    "secretion (Vms, Kms) and saturable tubular reabsorption (Vmr, Kmr),",
    "with net elimination = secretion - reabsorption. Oral bioavailability",
    "was modelled as 100% (passive diffusion; not altered by nutritional",
    "status). Protein-energy undernutrition (MAL_NOURISH = 1) reduces the",
    "secretion Vmax by 9.4% relative to well-nourished animals; no other",
    "PK parameter depends on nutritional status. Initial conditions in the",
    "central and peripheral compartments are set from the analytic",
    "positive root of the no-dose steady-state quadratic so that the",
    "endogenous taurine concentration is reproduced at t = 0."
  )
  reference <- paste(
    "Catalan-Latorre A, Nacher A, Merino V, Diez O, Merino Sanjuan M.",
    "A preclinical study to model taurine pharmacokinetics in the",
    "undernourished rat.",
    "British Journal of Nutrition. 2018;119(7):732-741.",
    "doi:10.1017/S0007114518000156.",
    sep = " "
  )
  vignette <- "Catalan-Latorre_2018_taurine_rat"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    MAL_NOURISH = list(
      description        = "Protein-energy undernutrition status indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = protein-energy undernourished (UN) per the paper's",
        "dual-criterion definition: end-of-adaptation body weight below",
        "80% of the well-nourished mean AND serum albumin below 23 g/L",
        "after 23-25 days on the TD 99168 5%-protein restricted diet",
        "(10 g per d, 159 kJ per d); 0 = well-nourished (WN, fed",
        "standard 14%-protein chow at 20 g per d, 251.9 kJ per d).",
        "Subject-level baseline indicator; the time-decay during",
        "supplementation is not modelled. Reduces the secretion Vmax",
        "(`vms`) by 9.4% in UN animals (Catalan-Latorre 2018 Table 3:",
        "FVms_DN = 0.906)."
      ),
      source_name        = "NUT"
    )
  )

  population <- list(
    species        = "rat (male Wistar)",
    n_subjects     = 64L,
    n_studies      = 1L,
    age_range      = "8-9 weeks at study start",
    weight_range   = "mean 235.5 (SD 7.9) g overall; well-nourished mean ~ 284.9 g; undernourished mean ~ 212.4 g at end of 23-25 d adaptation",
    sex_female_pct = 0,
    disease_state  = paste(
      "Healthy male Wistar rats randomised to one of two diets for a",
      "23-25 d adaptation period: well-nourished (n = 32; standard 14%",
      "protein chow, 20 g per d, 251.9 kJ per d) or undernourished",
      "(n = 32; TD 99168 5% protein restricted chow, 10 g per d, 159 kJ",
      "per d). 81% of the undernourished cohort (n = 26) developed",
      "moderate-to-severe undernutrition per the dual weight + serum",
      "albumin criterion."
    ),
    dose_range     = paste(
      "Single dose IV bolus or oral gavage of taurine (Sigma-Aldrich,",
      "Welwyn Garden City) dissolved in saline. Twelve groups of",
      "4-6 animals each crossing two routes (IV / oral) x three doses",
      "(1, 10, 100 mg per animal) x two nutritional groups (WN / UN)."
    ),
    regions        = "Spain (University of Valencia)",
    notes          = paste(
      "Baseline taurine plasma concentration was not statistically",
      "different between groups: WN 50.97 (SD 14.54) mg/L vs UN 48.86",
      "(SD 22.11) mg/L (P > 0.05). Plasma taurine quantified by HPLC",
      "with fluorimetric detection after OPA-MPA derivatisation,",
      "calibration range 50-750 uM (~6.3-94 mg/L), LOQ 57.11 uM.",
      "Modelled in NONMEM VI (FO method, ADVAN 3 / 11 / 9 explored).",
      "See Catalan-Latorre 2018 Methods (Protocol, animals and",
      "experimental procedures; Pharmacokinetic modelling and",
      "statistical analysis) and Tables 1-3 for the experimental",
      "design and fitted model output."
    )
  )

  ini({
    # ----------------------------------------------------------
    # Structural parameters -- Catalan-Latorre 2018 Table 3
    # (original-dataset Estimate column). Values are per animal,
    # not per kg; the mean rat body weight was 235.5 g.
    # CV% values from Table 3 are reproduced as trailing comments.
    #
    # Q0, Vms, Vmr are mass-flux rates (mg/h) entering the
    # Michaelis-Menten elimination form Vmax * Cc / (Km + Cc),
    # which has units of mg/h when Vmax is mg/h, Cc and Km in
    # mg/L. The paper's Table 3 header "(mg/l per h)" for Vms /
    # Vmr is inconsistent with the parameter-value sanity check
    # (a mass-rate interpretation reproduces the observed
    # ~50 mg/L baseline; a concentration-rate interpretation
    # would give a baseline two orders of magnitude too low).
    # See vignette Assumptions and deviations.
    # ----------------------------------------------------------
    lq0  <- log(13.7)    ; label("Endogenous taurine formation rate Q0 (mg/h per animal)")  # Table 3 Q0 = 13.7 (CV 36.0)
    lvc  <- log(0.0416)  ; label("Central volume of distribution Vc (L per animal)")        # Table 3 Vc = 0.0416 L (CV 17.5)
    lk12 <- log(2.61)    ; label("Inter-compartmental rate constant K12 (1/h)")             # Table 3 K12 = 2.61 (CV 44.8)
    lk21 <- log(0.73)    ; label("Inter-compartmental rate constant K21 (1/h)")             # Table 3 K21 = 0.73 (CV 48.1)
    lvms <- log(192.0)   ; label("Saturable tubular secretion Vmax Vms (mg/h per animal)")  # Table 3 Vms = 192.0 (CV 56.3)
    lkms <- log(399)     ; label("Saturable tubular secretion Km Kms (mg/L)")               # Table 3 Kms = 399 (CV 110)
    lvmr <- log(16.9)    ; label("Saturable tubular reabsorption Vmax Vmr (mg/h per animal)")  # Table 3 Vmr = 16.9 (CV 457)
    lkmr <- log(96.1)    ; label("Saturable tubular reabsorption Km Kmr (mg/L)")            # Table 3 Kmr = 96.1 (CV 190)
    lka  <- log(1.19)    ; label("First-order oral absorption rate ka (1/h)")               # Table 3 Ka  = 1.19 (CV 18.2)

    # ----------------------------------------------------------
    # Covariate effect on Vms.
    # Paper parameterisation (Eq. 3, Fig. 2 caption): Vms is
    # multiplied by FVms_DN^NUT, with FVms_DN = 0.906 and NUT = 1
    # for undernourished. nlmixr2lib canonical encoding uses the
    # additive offset form vms * (1 + e_<cov>_<param> * <COV>) so
    # that the reference subgroup (MAL_NOURISH = 0) recovers
    # exp(lvms); the offset coefficient is FVms_DN - 1 = -0.094,
    # i.e. a 9.4% reduction in Vms for UN animals.
    # ----------------------------------------------------------
    e_mal_nourish_vms <- -0.094 ; label("Fractional shift in Vms for malnutrition (UN) status (unitless)")  # Table 3 FVms_DN = 0.906 (CV 15.1), encoded as 0.906 - 1

    # ----------------------------------------------------------
    # IIV (BSV in paper) -- omega^2 = log(CV^2 + 1) conversion.
    # Table 3 reports BSV as CV% on the natural-scale parameter.
    #   Q0:  CV  25.7% -> log(1 + 0.257^2)  = 0.06397
    #   Vc:  CV  50.6% -> log(1 + 0.506^2)  = 0.22810
    #   K12: CV 120.4% -> log(1 + 1.204^2)  = 0.89571
    #   K21: CV  93.4% -> log(1 + 0.934^2)  = 0.62727
    #   Vms: CV  17.1% -> log(1 + 0.171^2)  = 0.02882
    #   ka:  CV  69.6% -> log(1 + 0.696^2)  = 0.39503
    # No IIV was estimated for Kms, Kmr, Vmr, or FVms_DN.
    # ----------------------------------------------------------
    etalq0  ~ 0.06397   # Table 3 BSV Q0  = 25.7%  CV (CV 104)
    etalvc  ~ 0.22810   # Table 3 BSV Vc  = 50.6%  CV (CV 52.7)
    etalk12 ~ 0.89571   # Table 3 BSV K12 = 120.4% CV (CV 38.7)
    etalk21 ~ 0.62727   # Table 3 BSV K21 = 93.4%  CV (CV 293)
    etalvms ~ 0.02882   # Table 3 BSV Vms = 17.1%  CV (CV 285)
    etalka  ~ 0.39503   # Table 3 BSV ka  = 69.6%  CV (CV 42.3)

    # ----------------------------------------------------------
    # Residual error -- proportional (CV 22.4%).
    # Methods describes a "slope-intercept" combined model with
    # a CV (proportional) and an SD (additive) term, but Table 3
    # reports only the single proportional sigma = 22.4%.
    # Encoded as proportional only. See vignette Assumptions.
    # ----------------------------------------------------------
    propSd <- 0.224 ; label("Proportional residual error on plasma taurine (fraction)")  # Table 3 sigma = 22.4% (CV 10.0)
  })

  model({
    # ----------------------------------------------------------
    # 1. Individual parameters
    # ----------------------------------------------------------
    q0  <- exp(lq0  + etalq0)
    vc  <- exp(lvc  + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    vms <- exp(lvms + etalvms)
    kms <- exp(lkms)
    vmr <- exp(lvmr)
    kmr <- exp(lkmr)
    ka  <- exp(lka  + etalka)

    # ----------------------------------------------------------
    # 2. Covariate effect: malnutrition reduces secretion Vmax.
    # ----------------------------------------------------------
    vms_eff <- vms * (1 + e_mal_nourish_vms * MAL_NOURISH)

    # ----------------------------------------------------------
    # 3. Endogenous-baseline initialisation.
    # At no-dose steady state (dCentral/dt = 0, depot = 0), the
    # central-compartment balance reduces -- after multiplying
    # out the two Michaelis-Menten denominators -- to a quadratic
    # in the steady-state plasma concentration Cp0:
    #
    #   A * Cp0^2 + B * Cp0 + C = 0
    #
    # where
    #   A = q0 - vms_eff + vmr
    #   B = q0 * (kms + kmr) - vms_eff * kmr + vmr * kms
    #   C = q0 * kms * kmr
    #
    # Under the typical-value parameter set A is negative
    # (vms_eff dominates q0 + vmr) and the discriminant is
    # positive, so the physiologically meaningful (positive)
    # root is (-B - sqrt(B^2 - 4AC)) / (2A). For WN
    # (MAL_NOURISH = 0) this gives Cp0 ~ 43.8 mg/L and for UN
    # (MAL_NOURISH = 1) Cp0 ~ 50.4 mg/L -- both within 1 SD of
    # the observed baselines (WN 50.97 +/- 14.54, UN 48.86 +/-
    # 22.11 mg/L; P > 0.05 between groups per Results).
    # Peripheral A2 is set so that the inter-compartmental flux
    # is balanced (K12 * A1 = K21 * A2) at t = 0.
    # ----------------------------------------------------------
    Acoef <- q0 - vms_eff + vmr
    Bcoef <- q0 * (kms + kmr) - vms_eff * kmr + vmr * kms
    Ccoef <- q0 * kms * kmr
    Cp0   <- (-Bcoef - sqrt(Bcoef * Bcoef - 4 * Acoef * Ccoef)) / (2 * Acoef)

    central(0)     <- Cp0 * vc
    peripheral1(0) <- k12 * Cp0 * vc / k21

    # ----------------------------------------------------------
    # 4. ODE system (Catalan-Latorre 2018 Eq. 3 + Fig. 2).
    # Parallel saturable tubular secretion (out) and reabsorption
    # (in) act on plasma; net elimination = secretion -
    # reabsorption. Oral bioavailability = 1.
    # ----------------------------------------------------------
    Cc <- central / vc

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot + q0 -
                          vms_eff * Cc / (kms + Cc) +
                          vmr     * Cc / (kmr + Cc) -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----------------------------------------------------------
    # 5. Observation and error
    # ----------------------------------------------------------
    Cc ~ prop(propSd)
  })
}
