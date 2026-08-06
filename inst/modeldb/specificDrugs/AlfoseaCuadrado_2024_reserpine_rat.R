AlfoseaCuadrado_2024_reserpine_rat <- function() {
  description <- paste(
    "Preclinical (rat). One-compartment population PK for reserpine with parallel",
    "first-order and zero-order absorption, coupled to a Sharma precursor-pool PK-PD",
    "model with a parallel three-transit chain for monoamine depletion in the amygdala,",
    "prefrontal cortex and spinal cord (reserpine-induced myalgia model of fibromyalgia)"
  )
  reference <- paste(
    "Alfosea-Cuadrado GM, Zarzoso-Foj J, Adell A, Valverde-Navarro AA,",
    "Gonzalez-Soler EM, Mangas-Sanjuan V, Blasco-Serra A.",
    "Population Pharmacokinetic-Pharmacodynamic Analysis of a Reserpine-Induced",
    "Myalgia Model in Rats. Pharmaceutics. 2024;16(8):1101.",
    "doi:10.3390/pharmaceutics16081101"
  )
  vignette <- "AlfoseaCuadrado_2024_reserpine_rat"
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. PK states verified against Figure 1 (schematic) and the
  # Figure 2 pc-VPC axis; PD states against Section 3.2.2 and Figure 1.
  compartmentData <- list(
    depot = list(
      analyte = "reserpine", units = "mg/kg", specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "reserpine", units = "mg/kg", specimen = "plasma",
      verified = TRUE
    ),
    precursor1 = list(
      analyte = "monoamine precursor pool", units = "mg/L", specimen = "not applicable",
      verified = TRUE
    ),
    effect = list(
      analyte = "monoamines (serotonin, dopamine, noradrenaline)", units = "mg/L",
      specimen = "tissue", verified = TRUE
    ),
    transit1 = list(
      analyte = "monoamine-degradation signal", units = "(unitless)",
      specimen = "not applicable", verified = TRUE
    ),
    transit2 = list(
      analyte = "monoamine-degradation signal", units = "(unitless)",
      specimen = "not applicable", verified = TRUE
    ),
    transit3 = list(
      analyte = "monoamine-degradation signal", units = "(unitless)",
      specimen = "not applicable", verified = TRUE
    )
  )

  covariateData <- list(
    CNSREG_PFC = list(
      description        = "Medial prefrontal cortex sampling-region indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (amygdala, the reference sampling region)",
      notes              = paste(
        "Per-observation anatomical sampling region. Selects the prefrontal-cortex",
        "typical value of the precursor production rate kpin. Mutually exclusive with",
        "CNSREG_SC; both 0 selects the amygdala reference. The only covariate retained",
        "in the final model (Section 3.2.2); body weight, age, breed and sex were all",
        "screened and rejected (Section 4, study limitations)."
      ),
      source_name        = "BRAIN AREA (level 'PFC')"
    ),
    CNSREG_SC = list(
      description        = "Lumbar spinal cord sampling-region indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (amygdala, the reference sampling region)",
      notes              = paste(
        "Per-observation anatomical sampling region. Selects the spinal-cord typical",
        "value of the precursor production rate kpin. Mutually exclusive with",
        "CNSREG_PFC; both 0 selects the amygdala reference."
      ),
      source_name        = "BRAIN AREA (level 'SC')"
    )
  )

  # Screened but not retained in the final model (Section 4): the paper reports no
  # usable point estimates for these, so they are documentation only.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "g",
      type        = "continuous",
      notes       = paste(
        "Rats weighed 300-450 g (Section 2.1). Screened as a covariate on the PK-PD",
        "parameters but not statistically significant, attributed by the authors to",
        "the low variation across animals (Section 4)."
      )
    ),
    AGE = list(
      description = "Age",
      units       = "weeks",
      type        = "continuous",
      notes       = "Screened and rejected (Section 4); no point estimate reported."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened and rejected (Section 4). All 120 animals were in fact male",
        "Sprague-Dawley rats (Section 2.1), so the screen carried no information."
      )
    )
  )

  population <- list(
    species       = "rat (male Sprague-Dawley)",
    n_subjects    = 120,
    n_studies     = 1,
    weight_range  = "300-450 g",
    sex_female_pct = 0,
    disease_state = "reserpine-induced myalgia (RIM) model of fibromyalgia syndrome",
    dose_range    = "0.1, 0.5 and 1 mg/kg once daily for three consecutive days",
    regions       = "Spain (University of Valencia)",
    n_observations = "120 PK and 828 PD observations, collected 48-96 h after the first dose",
    sampling      = paste(
      "Destructive sampling: one terminal sample per animal, at pre-dose or 30 min, 2, 4,",
      "24 or 48 h after the third dose (Section 2.1, Figure S1). Longitudinal profiles",
      "were assembled by pooling samples from different animals (Section 3.1), which the",
      "authors identify as the reason for the large inter-animal variability (Section 4)."
    ),
    endpoints     = paste(
      "Plasma reserpine by LC-MS (LLOQ 0.1 ug/mL); serotonin, dopamine and noradrenaline",
      "by HPLC with electrochemical detection in amygdala, medial prefrontal cortex and",
      "lumbar spinal cord homogenates."
    ),
    notes         = paste(
      "Sample counts per dose group, neurotransmitter and tissue are in Table 1.",
      "The analysis dataset is openly deposited at doi:10.5281/zenodo.11206173."
    )
  )

  ini({
    # -- Pharmacokinetics (Table 2, 'Population PKPD Model Estimates') ---------
    # Absorption is split between a first-order route (fraction F1) and a
    # zero-order route (fraction 1 - F1) delivered over a duration Tk0_2.
    lka <- fixed(log(19.14))
    label("First-order absorption rate constant (ka1, 1/h)")  # Table 2: ka1 = 19.14, FIX

    # Table 2 heads this row 'ka2 (mg/h/kg)', but both the Table 2 footnote and the
    # Figure 1 legend define it as the DURATION of zero-order absorption (Tk0_2).
    # The duration reading is the one that reproduces the Figure 2 pc-VPC tail; see
    # the vignette Errata.
    ld1 <- log(45.43)
    label("Duration of zero-order absorption (Tk0_2, h)")  # Table 2: ka2 = 45.43 (RSE 12%)

    # Logit rather than log scale so the two absorption fractions F1 and (1 - F1)
    # always split a non-negative dose; see the vignette Errata.
    logitfdepot <- log(0.95 / (1 - 0.95))
    label("Fraction of dose absorbed by the first-order route (F1, unitless)")  # Table 2: F1 = 0.95 (RSE 3%)

    lvc <- log(1.3e-3)
    label("Apparent central volume of distribution (V, L/kg)")  # Table 2: V = 1.3 mL/kg

    lcl <- log(4.5e-4)
    label("Apparent elimination clearance (CL, L/h/kg)")  # Table 2: CL = 4.5e-1 mL/h/kg

    # -- Precursor-pool pharmacodynamics (Table 2) ----------------------------
    # Region-specific typical values of the precursor production rate. Amygdala is
    # the reference level; the three rows share one inter-animal variability term.
    lkpin <- log(6.97)
    label("Precursor production rate, amygdala (kin, mg/L/h)")  # Table 2: kin AMY = 6.97 (RSE 18%)

    lkpin_pfc <- log(2.10)
    label("Precursor production rate, prefrontal cortex (kin, mg/L/h)")  # Table 2: kin PFC = 2.10 (RSE 18%)

    lkpin_sc <- log(1.78)
    label("Precursor production rate, spinal cord (kin, mg/L/h)")  # Table 2: kin SC = 1.78 (RSE 19%)

    lkin <- log(8.6e-4)
    label("Response production rate constant, precursor to response (kp, 1/h)")  # Table 2: kp = 8.6e-4 (RSE 14%)

    lkout <- log(2.7e-2)
    label("Response degradation rate constant (kout, 1/h)")  # Table 2: kout = 2.7e-2 (RSE 11%)

    lktr <- log(1.9e-1)
    label("Transit-chain rate constant (k0 / ktr, 1/h)")  # Table 2: k0 = 1.9e-1 (RSE 6%)

    lslp1 <- log(1.1e-1)
    label("Linear slope of reserpine on kp stimulation (SLP1, L/mg)")  # Table 2: SLP1 = 1.1e-1 (RSE 47%)

    lslp2 <- log(1.25)
    label("Linear slope of reserpine on transit-chain production (SLP2, L/mg)")  # Table 2: SLP2 = 1.25 (RSE 20%)

    # -- Inter-animal variability (Table 2, 'Inter-Animal Variability') -------
    # Reported as CV%; converted with omega^2 = log(1 + CV^2).
    etalka          ~ 1.80965   # Table 2: ka1 IAV 226% -> log(1 + 2.26^2)
    etald1          ~ 0.09750   # Table 2: ka2 IAV 32%  -> log(1 + 0.32^2)
    etalogitfdepot  ~ 1.43617   # Table 2: F1 IAV 179%  -> log(1 + 1.79^2)
    etalvc          ~ 0.29866   # Table 2: V IAV 59%    -> log(1 + 0.59^2)
    etalcl          ~ 0.12830   # Table 2: CL IAV 37%   -> log(1 + 0.37^2)
    etalkpin        ~ 0.66330   # Table 2: kin IAV 97%  -> log(1 + 0.97^2); shared by all three regions
    etalkin         ~ 0.08075   # Table 2: kp IAV 29%   -> log(1 + 0.29^2)
    etalkout        ~ 0.04727   # Table 2: kout IAV 22% -> log(1 + 0.22^2)
    etalktr         ~ 0.00807   # Table 2: k0 IAV 9%    -> log(1 + 0.09^2)
    etalslp1        ~ 2.62589   # Table 2: SLP1 IAV 358% -> log(1 + 3.58^2)
    etalslp2        ~ 0.43671   # Table 2: SLP2 IAV 74% -> log(1 + 0.74^2)

    # -- Residual unexplained variability (Table 2) ---------------------------
    propSd <- 0.54
    label("Proportional residual error, reserpine plasma concentration (fraction)")  # Table 2: PK 54% (RSE 9%)

    propSd_effect <- 0.71
    label("Proportional residual error, monoamine levels (fraction)")  # Table 2: PD 71% (RSE 9%)
  })

  model({
    # 1. Individual parameters
    ka     <- exp(lka + etalka)
    d1     <- exp(ld1 + etald1)
    fdepot <- exp(logitfdepot + etalogitfdepot) / (1 + exp(logitfdepot + etalogitfdepot))
    vc     <- exp(lvc + etalvc)
    cl     <- exp(lcl + etalcl)

    # Amygdala is the reference region (both indicators 0); the shared
    # inter-animal variability term applies to whichever region is selected.
    kpin <- (exp(lkpin) * (1 - CNSREG_PFC - CNSREG_SC) +
               exp(lkpin_pfc) * CNSREG_PFC +
               exp(lkpin_sc) * CNSREG_SC) * exp(etalkpin)

    kin  <- exp(lkin + etalkin)
    kout <- exp(lkout + etalkout)
    ktr  <- exp(lktr + etalktr)
    slp1 <- exp(lslp1 + etalslp1)
    slp2 <- exp(lslp2 + etalslp2)

    # 2. Micro-constants and the driving concentration
    kel <- cl / vc
    Cc  <- central / vc

    # 3. Absorption: the first-order fraction enters the depot, the complementary
    #    zero-order fraction is infused directly into central over Tk0_2
    #    (Figure 1, 'Depot 1 (F1)' / 'Depot 2 (1-F1)').
    f(depot)     <- fdepot
    f(central)   <- 1 - fdepot
    # NOTE (rxode2 5.1.3): this model declares two endpoints (Cc and effect), and
    # in that configuration the modelled duration is not evaluated for a dose
    # record flagged rate = -2 -- the solver reports "Duration is zero/negative".
    # Simulations must therefore supply the duration on the dose record itself
    # (`et(..., cmt = "central", dur = <Tk0_2>)`), which is numerically identical;
    # pass a per-subject value to retain the inter-animal variability on Tk0_2.
    # The declaration is kept because it is the structural statement of the model.
    dur(central) <- d1

    # 4. ODE system. Precursor-pool structure of Sharma 1998 (reference 40):
    #    reserpine stimulates the precursor-to-response transfer (SLP1) and, via a
    #    three-transit chain, the degradation of the response (SLP2). Figure 1.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    d/dt(precursor1) <- kpin - kin * (1 + slp1 * Cc) * precursor1
    d/dt(effect)     <- kin * (1 + slp1 * Cc) * precursor1 - kout * transit3 * effect

    d/dt(transit1) <- ktr * (1 + slp2 * Cc) - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3

    # 5. Baselines. The pre-dose system is at steady state, so precursor1(0) and
    #    effect(0) follow from setting the derivatives to zero. The paper states the
    #    transit-chain initial condition M10 = M20 = M30 = 1 (Section 3.2.2), which
    #    also makes the chain's own production and loss balance at ktr.
    precursor1(0) <- kpin / kin
    effect(0)     <- kpin / kout
    transit1(0)   <- 1
    transit2(0)   <- 1
    transit3(0)   <- 1

    # 6. Observations
    Cc ~ prop(propSd)
    effect ~ prop(propSd_effect)
  })
}
