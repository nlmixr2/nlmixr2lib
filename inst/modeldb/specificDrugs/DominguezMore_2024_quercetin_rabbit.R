DominguezMore_2024_quercetin_rabbit <- function() {
  description <- paste(
    "Preclinical (rabbit). Double first-order absorption, two-compartment oral",
    "population PK model for quercetin in male New Zealand White rabbits",
    "(Dominguez More 2024). The measured analyte is total quercetin released by",
    "enzymatic deconjugation of the rutin conjugates quercetin-3-O-glucuronide",
    "and quercetin-3-O-sulfate, expressed as rutin equivalents; free quercetin",
    "was never detected before deconjugation. Two parallel absorption sites",
    "reproduce the observed double concentration peak: a fast site (small",
    "intestine, ka_fast) taking the fraction frel of the dose, and a slow site",
    "(large intestine, after luminal efflux and microbial hydrolysis, ka_slow)",
    "taking the remainder and switched on only after a delay tlag. The",
    "P. peruviana calyx extract matrix acts as a categorical covariate on both",
    "absorption rate constants, on V and on both distribution micro-constants.",
    "NOTE: driven with the administered rutin dose this model over-predicts the",
    "published concentrations by roughly three orders of magnitude, because the",
    "paper reports no absolute bioavailability term for the metabolite; see",
    "population$notes and the validation vignette."
  )
  reference <- paste(
    "Dominguez More GP, Rey DP, Valderrama IH, Ospina LF, Aragon DM.",
    "Rutin and Physalis peruviana extract: population pharmacokinetics in",
    "New Zealand rabbits. Pharmaceutics. 2024;16(10):1241.",
    "doi:10.3390/pharmaceutics16101241",
    sep = " "
  )
  vignette <- "DominguezMore_2024_rutin_physalis"

  units <- list(time = "h", dosing = "ug/kg", concentration = "ng/mL")

  # Both absorption sites receive a dose record carrying the full oral dose;
  # f(depot) / f(depot2) split it. Declared explicitly because the automatic
  # detection only recognises `depot` and `central`.
  dosing <- c("depot", "depot2")

  compartmentData <- list(
    depot       = list(analyte = "quercetin conjugates", units = "ug/kg", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "quercetin conjugates", units = "ug/kg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "quercetin",            units = "ug/kg", specimen = "plasma",              verified = TRUE),
    peripheral1 = list(analyte = "quercetin",            units = "ug/kg", specimen = "plasma",              verified = TRUE)
  )

  covariateData <- list(
    FORM_RUTIN_EXTRACT = list(
      description        = "Source of the administered rutin: within the Physalis peruviana calyx extract matrix versus the isolated pure compound",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pure rutin)",
      notes              = paste(
        "Time-fixed per animal; each rabbit received only one treatment. The",
        "oral extract dose (500 mg/kg of extract, containing 14.80 ug rutin",
        "per mg) delivers 7.4 mg/kg of rutin against the 100 mg/kg pure-rutin",
        "dose, so the covariate is confounded with dose level in this study",
        "design; the authors did not fit a separate dose effect. Dominguez",
        "More 2024 Section 3.3 reports the extract affects ka1, ka2, V, k12",
        "and k21, but not k, F1 or Tlag2."
      ),
      source_name        = "source of rutin"
    )
  )

  population <- list(
    species        = "rabbit (New Zealand White, male)",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "9-10 weeks",
    weight_range   = "1.8-2.2 kg",
    sex_female_pct = 0,
    race_ethnicity = "not applicable",
    disease_state  = "healthy",
    dose_range     = paste(
      "Single oral gavage dose (0.5 mL/kg): pure rutin 100 mg/kg, or",
      "P. peruviana calyx extract 500 mg/kg equivalent to 7.4 mg/kg of rutin"
    ),
    regions        = "Colombia (Universidad Nacional de Colombia, Bogota)",
    notes          = paste(
      "The two oral arms (n = 5 each) of the four-arm study; the two",
      "intravenous arms fed the companion parent model",
      "DominguezMore_2024_rutin_rabbit. Plasma was sampled at 0, 0.083, 0.25,",
      "0.30, 0.75, 1, 2, 3, 4, 6, 8, 12, 24 and 48 h. Each sample was assayed",
      "twice: directly, for free rutin and free quercetin, and again after a",
      "30 min beta-glucuronidase/arylsulfatase deconjugation at pH 5.5 and",
      "37 degrees C. Free quercetin was never detected before deconjugation,",
      "so the modelled observation is the deconjugation product and represents",
      "the sum of quercetin-3-O-glucuronide and quercetin-3-O-sulfate, on a",
      "rutin-equivalent basis. Fitted in Monolix 2024R1; AIC 1244; all RSEs",
      "below 33%; validated by VPC and by a 1000-run bootstrap (100 runs did",
      "not converge). See Dominguez More 2024 Sections 2.2.4 and 3.3, and",
      "Table 4.",
      "",
      "ABSOLUTE-SCALE CAVEAT. Table 4 reports V as a primary (not apparent)",
      "volume and the model carries no absolute bioavailability parameter --",
      "frel and 1 - frel split the whole administered dose between the two",
      "absorption sites. Driven with the administered rutin dose, the",
      "published parameters therefore predict AUC(0-inf) = Dose / (k * V) =",
      "12,570 mg*h/L for the pure-rutin arm against the 9.28 mg*h/L reported",
      "in Table 2 (a factor of about 1360), and 6,280 mg*h/L against",
      "8.27 mg*h/L for the extract arm (a factor of about 760). Because the",
      "two factors differ by about 1.8-fold they cannot be a single units",
      "error: an unreported, arm-dependent bioavailability is folded into the",
      "published V. The model is shipped exactly as published; users who need",
      "absolute concentrations must supply their own bioavailability term,",
      "and the validation vignette therefore checks profile SHAPE (peak",
      "times, the double peak, terminal slope, MRT) rather than absolute",
      "exposure."
    )
  )

  ini({
    # ----------------------------------------------------------------------
    # Absorption -- Table 4, "Model Estimations / Population Value" column,
    # with the structure of Figure 3b: two parallel first-order absorption
    # sites, the second gated by a delay Tlag2 through the indicator
    # delta(t - Tlag2) (0 before Tlag2, 1 after). Section 4 identifies the
    # fast site as the small intestine (conjugates absorbed directly) and the
    # slow site as the large intestine, reached after luminal efflux of the
    # conjugates and microbial hydrolysis to the aglycone.
    #
    # ka1 / ka2 map onto the registered parallel-route canonicals lka_fast /
    # lka_slow rather than onto a positional lka1 / lka2 (see
    # inst/references/parameter-names.md).
    # ----------------------------------------------------------------------
    lka_fast <- log(11.146) ; label("Fast (small-intestinal) absorption rate constant ka1 (1/h)")        # Table 4: ka1 = 11.146 1/h (RSE 6.3%)
    lka_slow <- log(0.094)  ; label("Slow (large-intestinal) absorption rate constant ka2 (1/h)")        # Table 4: ka2 = 0.094 1/h (RSE 21.0%)
    ltlag    <- log(2.971)  ; label("Delay before the second absorption process starts, Tlag2 (h)")      # Table 4: Tlag2 = 2.971 h (RSE 2.8%)

    # F1 is the fraction absorbed at the FIRST site. Eq. 13 puts it on the
    # logit scale -- logit(F1) = logit(F1pop) + eta_F1 -- which is exactly the
    # registered logitfrel encoding, so the tabulated 0.270 is converted here.
    logitfrel <- log(0.270 / (1 - 0.270)) ; label("Logit of the fast-site absorbed fraction F1 (logit units)")  # Table 4: F1 = 0.270 (RSE 9.4%) -> logit = -0.9948

    # ----------------------------------------------------------------------
    # Disposition -- micro-constant parameterisation, as for the parent model.
    # ----------------------------------------------------------------------
    lvc  <- log(0.036) ; label("Central volume of distribution V (L/kg)")                                # Table 4: V = 0.036 L/kg (RSE 14.8%)
    lkel <- log(0.221) ; label("Elimination rate constant k (1/h)")                                      # Table 4: k = 0.221 1/h (RSE 20.6%)
    lk12 <- log(0.251) ; label("Central-to-peripheral distribution rate constant k12 (1/h)")             # Table 4: k12 = 0.251 1/h (RSE 28.3%)
    lk21 <- log(0.040) ; label("Peripheral-to-central distribution rate constant k21 (1/h)")             # Table 4: k21 = 0.040 1/h (RSE 32.4%)

    # ----------------------------------------------------------------------
    # Extract-matrix covariate effects (beta in Eq. 4). Section 3.3: "the
    # extract affected the two absorption rate constants (ka1 and ka2), the V,
    # and the two micro-constants of distribution (k12, k21)". Unlike the
    # intravenous parent model, EVERY parameter of this model that carries a
    # transformation is lognormal (Section 3.3: "a logit-normal transformation
    # applied to F1 and a lognormal transformation for the other parameters"),
    # so all five betas below are additive on the LOG scale -- Eqs. 11, 12,
    # 15, 17 and 18. No beta was estimated on k, F1 or Tlag2.
    # ----------------------------------------------------------------------
    e_form_rutin_extract_ka_fast <- -0.949 ; label("Extract-matrix shift on log(ka1) (log units)")       # Table 4: beta_ka1 = -0.949 (RSE 18.1%)
    e_form_rutin_extract_ka_slow <-  3.528 ; label("Extract-matrix shift on log(ka2) (log units)")       # Table 4: beta_ka2 = 3.528 (RSE 6.7%)
    e_form_rutin_extract_vc      <- -1.910 ; label("Extract-matrix shift on log(V) (log units)")         # Table 4: beta_V = -1.910 (RSE 8.0%)
    e_form_rutin_extract_k12     <-  1.076 ; label("Extract-matrix shift on log(k12) (log units)")       # Table 4: beta_k12 = 1.076 (RSE 33.5%)
    e_form_rutin_extract_k21     <-  2.067 ; label("Extract-matrix shift on log(k21) (log units)")       # Table 4: beta_k21 = 2.067 (RSE 16.0%)

    # ----------------------------------------------------------------------
    # Inter-individual variability. Table 4's sub-heading is "Standard
    # deviation of the Random Effects", so the tabulated omegas are SDs and
    # are squared here because nlmixr2's `~` takes a VARIANCE. Only F1, Tlag2
    # and k carry a random effect (Eqs. 13, 14 and 16); ka1, ka2, V, k12 and
    # k21 do not, and Table 4 lists no Omega for them. The F1 eta is additive
    # on the LOGIT scale, the other two on the log scale.
    # ----------------------------------------------------------------------
    etalogitfrel ~ 0.025600  # Table 4: Omega_F1    = 0.160 (SD, RSE 23.9%) -> 0.160^2
    etaltlag     ~ 0.005929  # Table 4: Omega_Tlag2 = 0.077 (SD, RSE 33.3%) -> 0.077^2
    etalkel      ~ 0.151321  # Table 4: Omega_k     = 0.389 (SD, RSE 23.9%) -> 0.389^2

    # ----------------------------------------------------------------------
    # Residual error -- combined (Section 3.3: "the residual error was better
    # described by a combined model"). Monolix's combined1 model is
    # y = f + (a + b * f) * e, so `a` is the additive SD in the concentration
    # unit and `b` the proportional SD.
    # ----------------------------------------------------------------------
    addSd  <- 20.661 ; label("Additive residual error (ng/mL)")                                          # Table 4: a = 20.661 (RSE 15.4%)
    propSd <- 0.044  ; label("Proportional residual error (fraction)")                                   # Table 4: b = 0.044 (RSE 32.4%)
  })

  model({
    # ----------------------------------------------------------------------
    # 1. Individual parameters.
    # ----------------------------------------------------------------------
    ka_fast <- exp(lka_fast + e_form_rutin_extract_ka_fast * FORM_RUTIN_EXTRACT)
    ka_slow <- exp(lka_slow + e_form_rutin_extract_ka_slow * FORM_RUTIN_EXTRACT)
    frel    <- expit(logitfrel + etalogitfrel)
    tlag    <- exp(ltlag + etaltlag)

    vc  <- exp(lvc  + e_form_rutin_extract_vc  * FORM_RUTIN_EXTRACT)
    kel <- exp(lkel + etalkel)
    k12 <- exp(lk12 + e_form_rutin_extract_k12 * FORM_RUTIN_EXTRACT)
    k21 <- exp(lk21 + e_form_rutin_extract_k21 * FORM_RUTIN_EXTRACT)

    # ----------------------------------------------------------------------
    # 2. ODE system (Figure 3b). depot is the first (fast, small-intestinal)
    #    absorption site and depot2 the second (slow, large-intestinal) one.
    #
    #    The paper writes the second site as dXa2/dt = -ka2 * delta(t-Tlag2) *
    #    Xa2, i.e. the amount sits inert in the site until Tlag2 and only then
    #    begins to empty. alag(depot2) is the identical process: no mass
    #    leaves depot2 before tlag, and the full (1 - frel) share is released
    #    from tlag onward at ka_slow.
    # ----------------------------------------------------------------------
    d/dt(depot)       <- -ka_fast * depot
    d/dt(depot2)      <- -ka_slow * depot2
    d/dt(central)     <-  ka_fast * depot + ka_slow * depot2 -
      (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)     <- frel
    f(depot2)    <- 1 - frel
    alag(depot2) <- tlag

    # ----------------------------------------------------------------------
    # 3. Observation. Doses are in ug/kg and vc is in L/kg, so central / vc is
    #    in ug/L = ng/mL. See the absolute-scale caveat in population$notes.
    # ----------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
