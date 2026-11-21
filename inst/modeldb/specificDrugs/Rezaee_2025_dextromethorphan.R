Rezaee_2025_dextromethorphan <- function() {
  description <- paste(
    "Joint parent-metabolite population PK/PD model for oral dextromethorphan",
    "and its active O-demethylated metabolite dextrorphan, with an",
    "effect-compartment model of citric-acid-induced cough suppression in",
    "healthy CYP2D6 normal-metaboliser volunteers (Rezaee 2025). Dextromethorphan",
    "dose: one depot draining to the one-compartment parent at ka and directly",
    "to the one-compartment metabolite (first-pass formation, no dose",
    "apportionment) at a separate rate constant, with systemic",
    "parent-to-metabolite conversion and equal apparent volumes for both.",
    "Dextrorphan dosed on its own: first-order absorption into a",
    "two-compartment disposition sharing the metabolite elimination rate",
    "constant. Cough count: baseline minus a placebo time course",
    "scale * k * t * exp(-k * t), multiplied by a sigmoid Imax term driven by",
    "the potency-weighted sum of the dextromethorphan and dextrorphan",
    "effect-site concentrations (relative potency of dextrorphan 0.26)."
  )
  reference <- paste(
    "Rezaee S, Wright CE, Morice AH, Rostami-Hodjegan A.",
    "Dextromethorphan Versus Dextrorphan: A Quantitative Comparison of",
    "Antitussive Potency Following Separate Administration of Metabolite.",
    "J Clin Pharmacol. 2025;65(11):1616-1625.",
    "doi:10.1002/jcph.70049."
  )
  vignette <- "Rezaee_2025_dextromethorphan"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    TRT_DEXTROMETHORPHAN = list(
      description = "Treatment-period indicator: 1 = the period in which the subject received oral dextromethorphan hydrobromide, 0 = otherwise (placebo or dextrorphan period)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo period when TRT_DEXTRORPHAN is also 0)",
      notes = paste(
        "Selects the additive residual SD of the cough count for the",
        "dextromethorphan arm (Table 3, 'Additive(DEX)' = 2.0). The drug",
        "effect itself is concentration-driven and needs no indicator.",
        "Mutually exclusive with TRT_DEXTRORPHAN within a period of the",
        "three-period crossover."
      ),
      source_name = "treatment arm (DEX)"
    ),
    TRT_DEXTRORPHAN = list(
      description = "Treatment-period indicator: 1 = the period in which the subject received oral dextrorphan D-tartrate on its own, 0 = otherwise (placebo or dextromethorphan period)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo or dextromethorphan period)",
      notes = paste(
        "Switches the dextrorphan disposition between the two sub-models the",
        "authors fitted: 0 = formed dextrorphan in a one-compartment",
        "disposition with the volume shared with the parent (Table 2 'V/F'),",
        "1 = directly administered dextrorphan in a two-compartment",
        "disposition (Table 2 'V1/F(DOR)', 'k12(DOR)', 'k21(DOR)'). Also",
        "selects the proportional residual error for dextrorphan",
        "(0.3 formed vs 0.2 directly administered) and the additive residual",
        "SD of the cough count for the dextrorphan arm (2.9)."
      ),
      source_name = "treatment arm (DOR)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "dextromethorphan", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dextromethorphan", units = "mg", specimen = "plasma", verified = TRUE),
    depot_dxor = list(analyte = "dextrorphan", units = "mg", specimen = "administration site", verified = TRUE),
    central_dxor = list(analyte = "dextrorphan", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_dxor = list(analyte = "dextrorphan", units = "mg", specimen = "plasma", verified = TRUE),
    effect = list(analyte = "dextromethorphan", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    effect_dxor = list(analyte = "dextrorphan", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 23L,
    n_studies = 1L,
    age_range = "19-51 years",
    sex_female_pct = 47.8,
    disease_state = paste(
      "Healthy non-smoking volunteers, genotyped CYP2D6 normal",
      "metabolisers, FEV1 > 80% predicted, screened to cough 7-14 times to",
      "five inhalations of 10% citric acid with <= 20% intra- and inter-day",
      "variability."
    ),
    dose_range = paste(
      "Single oral doses in a double-blind, randomised, placebo-controlled",
      "three-period crossover (>= 7 days between periods, Latin-square",
      "order): placebo; 40.5 mg dextromethorphan hydrobromide; 44.6 mg",
      "dextrorphan D-tartrate."
    ),
    regions = "United Kingdom (Sheffield)",
    notes = paste(
      "12 of 23 subjects male (Methods, 'Subjects and study design').",
      "Plasma dextromethorphan and unconjugated dextrorphan by HPLC with",
      "fluorometric detection (LLOQ 0.5 ng/mL dextromethorphan, 0.3 ng/mL",
      "dextrorphan) at 0, 1, 2, 3, 4, 5, 6, 8, 10, 12 and 24 h; cough",
      "counts at baseline (t = -1 h) and at the same times. BLQ",
      "concentrations left-censored. Sequential PK-then-PD fit in Monolix",
      "2024R1. Body weight and race are not reported."
    )
  )

  ini({
    # ---- Pharmacokinetics: Table 2 -------------------------------------
    # Dextromethorphan arm (joint parent-metabolite model, Figure 1).
    lka <- fixed(log(2.6))
    label("Absorption rate constant of dextromethorphan from the depot, ka(DEX) (1/h)") # Table 2: ka(DEX) = 2.6 h-1 (fixed); fixed to Moghadamnia 2003 per Results
    lk_dxor_form <- log(9.1)
    label("First-pass formation rate constant of dextrorphan from the dextromethorphan depot, ka(DORDEX) (1/h)") # Table 2: ka(DORDEX) = 9.1 h-1 (RSE 28%)
    lvc <- log(4874)
    label("Apparent volume shared by dextromethorphan and formed dextrorphan, V/F (L)") # Table 2: V/F = 4874 L (RSE 10%)
    lkel <- log(0.1)
    label("Elimination rate constant of dextromethorphan by routes other than dextrorphan formation, k(DEX) (1/h)") # Table 2: k(DEX) = 0.1 h-1 (RSE 36%)
    lkmet <- log(0.04)
    label("Systemic dextromethorphan-to-dextrorphan transformation rate constant, kpm (1/h)") # Table 2: kpm = 0.04 h-1 (RSE 44%)
    lkel_dxor <- log(0.2)
    label("Elimination rate constant of dextrorphan, k(DOR) (1/h)") # Table 2: k(DOR) = 0.2 h-1 (RSE 8%)

    # Dextrorphan administered on its own (two-compartment sub-model).
    lka_dxor <- fixed(log(2.6))
    label("Absorption rate constant of directly administered dextrorphan, ka(DOR) (1/h)") # Table 2: ka(DOR) = 2.6 h-1 (fixed)
    lvc_dxor <- log(6984)
    label("Apparent central volume of directly administered dextrorphan, V1/F(DOR) (L)") # Table 2: V1/F(DOR) = 6984 L (RSE 10%)
    lk12_dxor <- log(0.57)
    label("Central-to-peripheral rate constant of directly administered dextrorphan, k12(DOR) (1/h)") # Table 2: k12(DOR) = 0.57 h-1 (RSE 39%)
    lk21_dxor <- log(0.60)
    label("Peripheral-to-central rate constant of directly administered dextrorphan, k21(DOR) (1/h)") # Table 2: k21(DOR) = 0.60 h-1 (RSE 44%)

    # ---- Placebo cough model: Table 3 and Methods equation 1 ------------
    lrbase <- log(10.7)
    label("Baseline cough count before any treatment, Baseline_cough (coughs)") # Table 3: Baseline_cough = 10.7 (RSE 3%)
    lscale_pbo <- log(4.5)
    label("Magnitude of the placebo cough reduction, scale (coughs)") # Table 3: scale = 4.5 (RSE 32%)
    lkpbo <- log(0.2)
    label("Rate constant of placebo cough suppression and return to baseline, K (1/h)") # Table 3: K = 0.2 h-1 (RSE 29%)

    # ---- PK-PD link: Table 3 and Methods equation 2 ---------------------
    lke0 <- log(0.08)
    label("Effect-site equilibration rate constant of dextromethorphan, ke0(DEX) (1/h)") # Table 3: ke0(DEX) = 0.08 h-1 (RSE 110%)
    lke0_dxor <- log(0.02)
    label("Effect-site equilibration rate constant of dextrorphan, ke0(DOR) (1/h)") # Table 3: ke0(DOR) = 0.02 h-1 (RSE 344%)
    lic50 <- log(0.3)
    label("Effect-site dextromethorphan-equivalent concentration giving half of Imax, IC50 (ng/mL)") # Table 3: IC50 = 0.3 ng/mL (RSE 327%)
    limax <- log(0.23)
    label("Maximum fractional cough inhibition by drug, Imax (fraction of the placebo-model count)") # Table 3: Imax = 0.23 (RSE 27%)
    lhill <- log(2.8)
    label("Hill coefficient of the inhibition, n (unitless)") # Table 3: N = 2.8 (RSE 98%)
    lrelpot_dxor <- log(0.26)
    label("Potency of dextrorphan relative to dextromethorphan, Pot (unitless)") # Table 3: Relative potency = 0.26 (RSE 96%)

    # ---- Between-subject variability: log-normal, omega^2 = log(1 + CV^2)
    etalk_dxor_form ~ log(1 + 1.98^2) # Table 2: BSV ka(DORDEX) 198% CV
    etalvc ~ log(1 + 0.48^2) # Table 2: BSV V/F 48% CV
    etalkel ~ log(1 + 1.08^2) # Table 2: BSV k(DEX) 108% CV
    etalkel_dxor ~ log(1 + 0.32^2) # Table 2: BSV k(DOR) 32% CV
    etalkmet ~ log(1 + 0.43^2) # Table 2: BSV kpm 43% CV
    etalvc_dxor ~ log(1 + 0.32^2) # Table 2: BSV V1/F(DOR) 32% CV
    etalk12_dxor ~ log(1 + 1.96^2) # Table 2: BSV k12(DOR) 196% CV
    etalk21_dxor ~ log(1 + 2.99^2) # Table 2: BSV k21(DOR) 299% CV
    etalrbase ~ log(1 + 0.14^2) # Table 3: BSV Baseline_cough 14% CV
    etalscale_pbo ~ log(1 + 1.15^2) # Table 3: BSV Scale 115% CV
    etalkpbo ~ log(1 + 0.85^2) # Table 3: BSV K 85% CV
    etalke0 ~ log(1 + 1.76^2) # Table 3: BSV ke0(DEX) 176% CV
    etalke0_dxor ~ log(1 + 1.71^2) # Table 3: BSV ke0(DOR) 171% CV
    etalic50 ~ log(1 + 1074.25^2) # Table 3: BSV IC50 107,425% CV as printed
    etalimax ~ log(1 + 0.68^2) # Table 3: BSV Imax 68% CV
    etalhill ~ log(1 + 2.05^2) # Table 3: BSV N 205% CV
    etalrelpot_dxor ~ log(1 + 0.27^2) # Table 3: BSV Relative potency 27% CV

    # ---- Residual error ------------------------------------------------
    propSd <- 0.3
    label("Proportional residual error, dextromethorphan (fraction)") # Table 2: Proportional (DEX) = 0.3 (RSE 8%)
    propSd_dxor <- 0.3
    label("Proportional residual error, dextrorphan formed after dextromethorphan (fraction)") # Table 2: Proportional (DORDEX) = 0.3 (RSE 6%)
    propSd_dxor_direct <- 0.2
    label("Proportional residual error, dextrorphan after its own administration (fraction)") # Table 2: Proportional (DOR) = 0.2 (RSE 7%)
    addSd_cough_placebo <- 1.4
    label("Additive residual error of the cough count, placebo period (coughs)") # Table 3: Additive(Placebo) = 1.4 (RSE 5%)
    addSd_cough_dextromethorphan <- 2.0
    label("Additive residual error of the cough count, dextromethorphan period (coughs)") # Table 3: Additive(DEX) = 2.0 (RSE 5%)
    addSd_cough_dextrorphan <- 2.9
    label("Additive residual error of the cough count, dextrorphan period (coughs)") # Table 3: Additive(DOR) = 2.9 (RSE 5%)
  })

  model({
    # Individual PK parameters
    ka <- exp(lka)
    k_dxor_form <- exp(lk_dxor_form + etalk_dxor_form)
    vc <- exp(lvc + etalvc)
    kel <- exp(lkel + etalkel)
    kmet <- exp(lkmet + etalkmet)
    kel_dxor <- exp(lkel_dxor + etalkel_dxor)
    ka_dxor <- exp(lka_dxor)
    vc_dxor <- exp(lvc_dxor + etalvc_dxor)
    k12_dxor <- exp(lk12_dxor + etalk12_dxor)
    k21_dxor <- exp(lk21_dxor + etalk21_dxor)

    # Dextrorphan disposition differs by period: formed dextrorphan shares the
    # parent's apparent volume and is one-compartment; dextrorphan given on its
    # own has its own central volume and a peripheral compartment. k(DOR) is
    # common to both.
    v_dxor <- vc * (1 - TRT_DEXTRORPHAN) + vc_dxor * TRT_DEXTRORPHAN
    k12_dxor_on <- k12_dxor * TRT_DEXTRORPHAN

    # Individual PD parameters
    rbase <- exp(lrbase + etalrbase)
    scale_pbo <- exp(lscale_pbo + etalscale_pbo)
    kpbo <- exp(lkpbo + etalkpbo)
    ke0 <- exp(lke0 + etalke0)
    ke0_dxor <- exp(lke0_dxor + etalke0_dxor)
    ic50 <- exp(lic50 + etalic50)
    imax <- exp(limax + etalimax)
    hill <- exp(lhill + etalhill)
    relpot_dxor <- exp(lrelpot_dxor + etalrelpot_dxor)

    # Dextromethorphan dose: a single depot drains in parallel to the parent
    # (ka) and, by first-pass conversion, to the metabolite (k_dxor_form);
    # no dose apportionment (Figure 1).
    d/dt(depot) <- -(ka + k_dxor_form) * depot
    d/dt(central) <- ka * depot - (kel + kmet) * central
    d/dt(depot_dxor) <- -ka_dxor * depot_dxor
    d/dt(central_dxor) <- k_dxor_form * depot + kmet * central + ka_dxor * depot_dxor -
      (kel_dxor + k12_dxor_on) * central_dxor + k21_dxor * peripheral1_dxor
    d/dt(peripheral1_dxor) <- k12_dxor_on * central_dxor - k21_dxor * peripheral1_dxor

    # Plasma concentrations: amounts in mg, volumes in L -> x 1000 for ng/mL
    Cc <- 1000 * central / vc
    Cc_dxor <- 1000 * central_dxor / v_dxor

    # Effect-site concentrations (ng/mL)
    d/dt(effect) <- ke0 * (Cc - effect)
    d/dt(effect_dxor) <- ke0_dxor * (Cc_dxor - effect_dxor)

    # Placebo time course: the cough count falls by scale * k * t * exp(-k * t)
    # (maximum scale / e at t = 1 / k) and returns to baseline. Time is from the
    # dose, so the model must be simulated with the dose at t = 0.
    cough_placebo <- rbase - scale_pbo * kpbo * t * exp(-kpbo * t)

    # Drug inhibition by the potency-weighted sum of the two effect-site
    # concentrations acting on a common receptor (Methods equation 2)
    drive <- effect / ic50 + effect_dxor * relpot_dxor / ic50
    # imax * D^n / (1 + D^n) written as imax / (1 + D^-n): identical for D > 0,
    # and finite when the very wide IC50 variability (omega^2 = 13.96) makes D^n
    # overflow to Inf, where the printed form gives Inf / Inf = NaN. At D = 0
    # it evaluates to imax / Inf = 0.
    inh <- imax / (1 + drive^(-hill))
    cough <- cough_placebo * (1 - inh)

    # Residual errors are period-specific
    sd_dxor <- propSd_dxor * (1 - TRT_DEXTRORPHAN) + propSd_dxor_direct * TRT_DEXTRORPHAN
    sd_cough <- addSd_cough_placebo * (1 - TRT_DEXTROMETHORPHAN - TRT_DEXTRORPHAN) +
      addSd_cough_dextromethorphan * TRT_DEXTROMETHORPHAN +
      addSd_cough_dextrorphan * TRT_DEXTRORPHAN

    Cc ~ prop(propSd)
    Cc_dxor ~ prop(sd_dxor)
    cough ~ add(sd_cough)
  })
}
