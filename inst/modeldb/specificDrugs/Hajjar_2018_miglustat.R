Hajjar_2018_miglustat <- function() {
  description <- paste(
    "Two-compartment population PK model for oral miglustat (AT2221;",
    "N-butyl-1-deoxynojirimycin) administered as a pharmacological",
    "chaperone for cipaglucosidase alfa in adult patients with Pompe",
    "disease (Hajjar 2018 ACCP poster, phase 1/2 study ATB200-02 /",
    "NCT02675465). Absorption is described by a sequential",
    "zero-order release into the depot (duration D1 = 0.459 h) followed",
    "by first-order absorption (Ka = 0.485 /h) into the central",
    "compartment. Apparent disposition parameters (CL/F = 8.55 L/h,",
    "Vc/F = 36.3 L, Q/F = 3.16 L/h, Vp/F = 45.6 L) are reported at the",
    "70 kg reference body weight and allometrically scaled with",
    "exponents 0.75 fixed on clearances and 1 fixed on volumes.",
    "Bioavailability F is not estimable from oral-only data and is",
    "anchored at 1 (apparent CL/F and V/F parameterisation). Residual",
    "error is proportional (variance 0.0408, SD 0.202 on the",
    "linear-scale concentration)."
  )
  reference <- paste(
    "Hajjar JL, Mondick JT, Gastonguay MR, Mulberg AE, Johnson FK.",
    "Population pharmacokinetic modeling of enzyme replacement therapy",
    "ATB200 and pharmacological chaperone AT2221 in adult patients with",
    "Pompe disease and simulation to predict adolescent exposures. ACCP",
    "Annual Meeting Poster; September 2018; Seattle, WA.",
    "https://metrumrg.com/wp-content/uploads/Pubs/2018-ACCP-Population-PK-of-ATB200-AT221-in-Pompe-Patients_2018-09-18-Poster_L1e.pdf"
  )
  vignette <- "Hajjar_2018_pompe_disease"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight (kg). Drives the allometric scaling of all",
        "apparent clearances (exponent 0.75 fixed) and all apparent",
        "volumes (exponent 1.0 fixed) with reference 70 kg per Hajjar",
        "2018 Methods 'Modeling' bullet 'Clearance and volume parameters",
        "were allometrically scaled by individual weights normalized to",
        "70 kg body weight, with exponents fixed to values of 0.75 and",
        "1, respectively'. Adult cohort baseline mean age 49.4 years;",
        "sex 10 M / 5 F across the 15 adults included in the modeling",
        "(Hajjar 2018 Table 1)."
      ),
      source_name        = "WT"
    )
  )

  covariatesDataExcluded <- list(
    BACT = list(
      description        = "Prior enzyme-replacement-therapy experience indicator (1 = ERT-experienced with alglucosidase alfa; 0 = ERT-naive)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ERT-naive)",
      notes              = paste(
        "Screened as a categorical covariate on apparent CL/F via the",
        "model CL * theta_ERT (Hajjar 2018 Methods 'Modeling' last two",
        "bullets). NOT retained in the final model: the mean ERT",
        "covariate effect was 0.82 with 95% CI [0.41, 1.23] for AT2221,",
        "and the small ERT-naive sample size (n = 5) was deemed",
        "insufficient to estimate an effect of prior ERT experience",
        "(Results 'Population Pharmacokinetic Models' bullets 4-6).",
        "Documented here for covariate-screen provenance; not referenced",
        "inside model()."
      ),
      source_name        = "ERT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 15L,
    n_studies      = 1L,
    age_range      = "24-66 years",
    age_mean       = "49.4 years",
    weight_range   = "not reported (paper states only that allometric scaling normalised to 70 kg)",
    sex_female_pct = 100 * 6 / 15,
    disease_state  = paste(
      "Adults with Pompe disease (genetic deficiency of acid",
      "alpha-glucosidase / GAA). 10 ERT-experienced adults previously",
      "treated with alglucosidase alfa (mean 4.8 years on prior ERT,",
      "SD 1.42) and 5 ERT-naive adults. One additional ERT-experienced",
      "subject was excluded from the modeling for missing PK values",
      "(Table 1 footnote a). Baseline 6-Minute Walk Test 392 m",
      "(SD 93) ERT-experienced vs 400 m (SD 84) ERT-naive; baseline",
      "forced vital capacity (upright) 52% predicted (SD 13) vs 53%",
      "predicted (SD 20)."
    ),
    dose_range     = paste(
      "ERT-experienced (n = 10): two AT2221 dose-occasion levels of",
      "130 mg and 260 mg, each co-administered with 20 mg/kg ATB200.",
      "ERT-naive (n = 5): a single 260 mg AT2221 dose co-administered",
      "with 20 mg/kg ATB200. Plasma samples collected over 24 h periods",
      "for both ATB200 and AT2221 concentration assays."
    ),
    regions        = "United States; ATB200-02 (NCT02675465), a phase 1/2 study sponsored by Amicus Therapeutics.",
    notes          = paste(
      "Adult cohort only; the same poster reports Monte-Carlo",
      "simulations forward-projecting to adolescents (12 to <18 years",
      "old). The structural model is fit to adult data only.",
      "Parameter estimation used NONMEM 7.3 with First Order Conditional",
      "Estimation (FOCE) and Maximum Likelihood Estimation. 500-trial",
      "Monte-Carlo predictive checks supported the structural model.",
      "Note: miglustat is the INN for AT2221 (N-butyl-1-deoxynojirimycin);",
      "it is approved for type-1 Gaucher disease at higher doses (Zavesca",
      "100 mg t.i.d.). At the 130 / 260 mg single-dose levels used in",
      "ATB200-02 it is co-administered as a pharmacological chaperone",
      "to stabilise the recombinant GAA (cipaglucosidase alfa). The",
      "2018 ACCP poster source carries no peer-reviewed DOI; the on-disk",
      "metadata DOI 10.1038/mt.2009.53 is a vendor-index defect and is",
      "not used as the canonical reference."
    )
  )

  ini({
    # =========================================================================
    # AT2221 (miglustat, oral) structural parameters at the reference subject
    # WT = 70 kg. Estimated point values from Hajjar 2018 Table 2 'AT2221 model
    # parameter estimates' column. All values entered as point estimates (no
    # fixed() wrapper) because the poster reports a non-zero RSE for each.
    # =========================================================================
    lcl  <- log(8.55); label("Apparent clearance CL/F (L/h, at WT 70 kg)")                                        # Hajjar 2018 Table 2: CL/F = 8.55 L/h (RSE 11.6%)
    lvc  <- log(36.3); label("Apparent central volume of distribution Vc/F (L, at WT 70 kg)")                     # Hajjar 2018 Table 2: Vc/F = 36.3 L (RSE 21.3%)
    lq   <- log(3.16); label("Apparent inter-compartmental clearance Q/F (L/h, at WT 70 kg)")                     # Hajjar 2018 Table 2: Q/F = 3.16 L/h (RSE 26.3%)
    lvp  <- log(45.6); label("Apparent peripheral volume of distribution Vp/F (L, at WT 70 kg)")                  # Hajjar 2018 Table 2: Vp/F = 45.6 L (RSE 29.1%)

    # =========================================================================
    # Absorption: sequential zero-order release into depot (duration D1) +
    # first-order absorption from depot to central at rate Ka. Hajjar 2018
    # Methods 'Modeling' bullet 'AT2221 disposition was described by a
    # 2-compartment model with sequential zero and first-order absorption'.
    # =========================================================================
    lka <- log(0.485); label("First-order absorption rate constant Ka (1/h)")                                     # Hajjar 2018 Table 2: Ka = 0.485 /h (RSE 52.0%)
    ld1 <- log(0.459); label("Zero-order absorption duration D1 (h) into the depot")                              # Hajjar 2018 Table 2: D1 = 0.459 h (RSE 27.9%)

    # Bioavailability anchor. Hajjar 2018 reports apparent CL/F and V/F (oral
    # study only, no IV reference arm), so absolute F is not identifiable.
    # Anchor F = 1 so the population disposition parameters retain their
    # 'apparent' meaning.
    lfdepot <- fixed(log(1)); label("Oral bioavailability F (anchored at 1; apparent CL/F and V/F parameterisation)")  # Apparent-parameter convention; F not estimable from oral-only data

    # =========================================================================
    # Allometric exponents (shared with the ATB200 model; Hajjar 2018 Methods
    # 'Modeling' bullet: exponents fixed to 0.75 and 1).
    # =========================================================================
    allo_cl <- fixed(0.75);  label("Allometric exponent on clearances (unitless; fixed at theoretical 0.75)")     # Hajjar 2018 Methods 'Modeling'
    allo_v  <- fixed(1.00);  label("Allometric exponent on volumes (unitless; fixed at theoretical 1.00)")        # Hajjar 2018 Methods 'Modeling'

    # =========================================================================
    # Between-subject variability (Hajjar 2018 Table 2 'BSV' column). Lognormal
    # variances computed as omega^2 = log(1 + CV^2):
    #   CL/F 17.50% -> log(1 + 0.175^2) = 0.03016541 (RSE of BSV 90.1%)
    #   Q/F  43.40% -> log(1 + 0.434^2) = 0.17257080 (RSE of BSV 108%)
    #   Ka   44.30% -> log(1 + 0.443^2) = 0.17919080 (RSE of BSV 70.8%)
    # Vc/F, Vp/F, and D1 do not have a BSV reported (single-value cells), so
    # no eta is included for them. Single etas; no correlation block reported.
    # =========================================================================
    etalcl ~ 0.03016541  # Hajjar 2018 Table 2: BSV CL/F 17.50% (RSE 90.1%)
    etalq  ~ 0.17257080  # Hajjar 2018 Table 2: BSV Q/F 43.40% (RSE 108%)
    etalka ~ 0.17919080  # Hajjar 2018 Table 2: BSV Ka 44.30% (RSE 70.8%)

    # =========================================================================
    # Residual error. Table 2 reports 'Variance of residual error = 0.0408
    # (RSE 4.88%)'. AT2221 concentrations are in the 1-25 ug/mL range
    # (Figure 2 visual predictive check axes), so the small variance is
    # consistent with a proportional / log-additive residual rather than a
    # linear-scale additive variance. propSd = sqrt(0.0408) = 0.202.
    # =========================================================================
    propSd <- 0.2019901; label("Proportional residual error (fraction CV) on linear-scale concentration")          # Hajjar 2018 Table 2: variance = 0.0408; SD = sqrt(0.0408)
  })

  model({
    # Reference body weight for allometric scaling.
    ref_wt <- 70

    # ---- Individual PK parameters ----------------------------------------
    cl <- exp(lcl + etalcl) * (WT / ref_wt) ^ allo_cl
    vc <- exp(lvc)          * (WT / ref_wt) ^ allo_v
    q  <- exp(lq  + etalq)  * (WT / ref_wt) ^ allo_cl
    vp <- exp(lvp)          * (WT / ref_wt) ^ allo_v
    ka <- exp(lka + etalka)
    d1 <- exp(ld1)
    fdepot <- exp(lfdepot)

    # Micro-rate constants for the explicit two-compartment ODE form.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- ODE system ------------------------------------------------------
    # Sequential zero-order release into depot (duration D1, set via
    # `dur(depot) <- d1`) followed by first-order absorption Ka into the
    # central compartment.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    dur(depot) <- d1
    f(depot)   <- fdepot

    # Concentration in the central compartment (dose in mg, vc in L => mg/L,
    # which equals ug/mL for matching the paper's plasma-concentration axes).
    Cc <- central / vc

    # Residual error: proportional on the linear-scale concentration.
    Cc ~ prop(propSd)
  })
}
