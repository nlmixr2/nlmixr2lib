Bulitta_2009_cefuroxime_axetil <- function() {
  description <- paste(
    "Semiphysiological population PK model for oral cefuroxime axetil",
    "(acetoxyethyl-ester prodrug of cefuroxime) in healthy adult male",
    "volunteers after a standardized high-fat breakfast. Three drug",
    "compartments (stomach -> intestine -> central): a saturable,",
    "time-dependent Michaelis-Menten release from the stomach to the",
    "intestine followed by first-order absorption from the intestine to",
    "the central compartment, with one-compartment linear disposition.",
    "The maximum gastric-release rate Vmax is modulated over time-past-",
    "meal by a sigmoidal (Hill) function whose maximum fractional",
    "change Emax is logit-transformed to range over [-1, 9]. Vmax and",
    "Km are estimated as fractions of dose (Bulitta 2009 Results,",
    "'estimated and are reported as fractions of the cefuroxime dose'),",
    "so the absolute Vmax (mg/h) and Km (mg) scale linearly with the",
    "stomach-compartment dose amount. Parameter values reproduced here",
    "are from the S-ADAPT importance-sampling Monte Carlo EM fit",
    "(Bulitta 2009 Table 2, 'S-ADAPT Population mean' column), which",
    "the authors recommend as the best parametric fit; NONMEM and NPAG",
    "estimates are reported alongside in the paper for comparison."
  )
  reference <- paste(
    "Bulitta JB, Landersdorfer CB, Kinzig M, Holzgrabe U, Sorgel F.",
    "New semiphysiological absorption model to assess the",
    "pharmacodynamic profile of cefuroxime axetil using nonparametric",
    "and parametric population pharmacokinetics. Antimicrobial Agents",
    "and Chemotherapy. 2009;53(8):3462-3471.",
    "doi:10.1128/AAC.00054-09"
  )
  vignette <- "Bulitta_2009_cefuroxime_axetil"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species         = "human",
    n_subjects      = 24L,
    n_studies       = 1L,
    age_range       = "18-31 years",
    age_median      = "24.5 years (mean)",
    weight_range    = "58.2-93.6 kg",
    weight_median   = "73.8 kg (mean)",
    height_range    = "166-193 cm (mean 179, SD 8.0)",
    sex_female_pct  = 0,
    race_ethnicity  = "Caucasian (24/24)",
    disease_state   = "Healthy volunteers",
    dose_range      = paste(
      "Single oral dose of 300.72 mg cefuroxime axetil suspension",
      "(equivalent to 250 mg cefuroxime) with 240 mL water immediately",
      "after a standardized high-fat breakfast (Bulitta 2009 'Study",
      "design and drug administration')."
    ),
    regions         = "Germany (Institute for Biomedical and Pharmaceutical Research, Nurnberg)",
    notes           = paste(
      "Single-dose, single-centre healthy-volunteer bioavailability",
      "study (Bulitta 2009 'Subjects'). Blood samples drawn at predose",
      "and 0.5, 1, 1.5, 2, 2.33, 2.67, 3, 3.33, 3.67, 4, 4.5, 5, 6, 8,",
      "10, and 12 h post-dose. Cefuroxime quantified by LC-MS/MS",
      "(LLOQ 0.00900 mg/L). Parameter estimates reported here are the",
      "S-ADAPT importance-sampling parametric Monte Carlo EM fit",
      "(pmethod=8) of the final semiphysiological absorption model."
    )
  )

  ini({
    # =============================================================
    # Structural parameters. Population means (typical values) from
    # the S-ADAPT final fit (Bulitta 2009 Table 2, 'S-ADAPT
    # Population mean' column). Relative standard errors in % from
    # the same column are annotated in source-trace comments.
    # =============================================================
    lcl     <- log(21.7);  label("Apparent total clearance CL/F (L/h)")                                                                    # Table 2 S-ADAPT mean: 21.7 (4.1% SE)
    lvc     <- log(38.7);  label("Apparent central volume V/F (L)")                                                                        # Table 2 S-ADAPT mean: 38.7 (4.1% SE)
    lvmax   <- log(0.505); label("Maximum gastric-release rate Vmax0 expressed as a fraction of dose (1/h)")                               # Table 2 S-ADAPT mean: 0.505 (22% SE); multiplied by podo(stomach) inside model() so absolute Vmax (mg/h) is dose-proportional
    lkm     <- log(0.426); label("Michaelis constant Km for stomach release expressed as a fraction of dose (unitless)")                   # Table 2 S-ADAPT mean: 42.6% = 0.426 (39% SE); multiplied by podo(stomach) inside model() so absolute Km (mg) is dose-proportional
    ltabs   <- log(9.34);  label("Absorption half-life T_abs from intestine to central (min)")                                              # Table 2 S-ADAPT mean: 9.34 min (15% SE); footnote f gives ka = ln(2) / (T_abs / 60) in 1/h
    ltc50   <- log(1.61);  label("Time past meal at which Vmax changed by 50% TC50 (h)")                                                   # Table 2 S-ADAPT mean: 1.61 h (20% SE)
    lgtemax <- -0.762;     label("Lg_Emax: logit-transformed maximum fractional change in Vmax (transformed scale, unitless)")              # Table 2 S-ADAPT mean: -0.762 (31% SE); back-transform Emax = 10*exp(lgtemax)/(1+exp(lgtemax)) - 1 maps to (-1, 9) per 'Parameter variability and observation model'
    lhill   <- fixed(log(10)); label("Hill coefficient gamma for the time-dependent Vmax modulator (unitless)")                             # Table 2 footnote g: initial estimates 10-15; held fixed at 10 in the final fit to improve stability

    # =============================================================
    # Inter-individual variability: full 7x7 variance-covariance
    # matrix on the transformed (log / logit-Emax) scale, copied
    # directly from Bulitta 2009 Table 3, 'Variance-covariance
    # matrix for the final population PK model in S-ADAPT'.
    # The diagonals match (BSV from Table 2)^2 to within rounding;
    # Lg_Emax variance is on the logit-transformed scale per
    # Table 2 footnote d.
    # =============================================================
    # Row 1: Var(CL/F);
    # Row 2: Cov(V/F, CL/F), Var(V/F);
    # Row 3: Cov(Vmax0/dose, .), Var(Vmax0/dose);
    # Row 4: Cov(Km/dose,   .), Var(Km/dose);
    # Row 5: Cov(T_abs,     .), Var(T_abs);
    # Row 6: Cov(TC50,      .), Var(TC50);
    # Row 7: Cov(Lg_Emax,   .), Var(Lg_Emax).
    etalcl + etalvc + etalvmax + etalkm + etaltabs + etaltc50 + etalgtemax ~ c(
      0.0393,
      0.0326,  0.0333,
      0.1143,  0.0877,  0.8100,
      0.2693,  0.2066,  1.5925,  3.2560,
      0.0083,  0.0142, -0.1036, -0.1640,  0.3604,
      0.0931,  0.0389,  0.5324,  1.0200,  0.0459,  0.8945,
      0.1530,  0.1131,  0.5702,  1.1824,  0.0414,  0.8197, 0.9793
    )

    # =============================================================
    # Residual error: combined additive + proportional model
    # (Bulitta 2009 'Parameter variability and observation model':
    # 'a combined additive and proportional error model').
    # =============================================================
    propSd <- 0.0851; label("Proportional residual error (fraction)")    # Table 2 S-ADAPT: 0.0851 (5.5% SE on variance)
    addSd  <- 0.0029; label("Additive residual error (mg/L)")            # Table 2 S-ADAPT: 0.0029 mg/L (48% SE on variance)
  })

  model({
    # -----------------------------------------------------------
    # Individual fixed-effect parameters on the natural / dose-
    # fraction scale.
    # -----------------------------------------------------------
    cl         <- exp(lcl     + etalcl)
    vc         <- exp(lvc     + etalvc)
    vmax0_frac <- exp(lvmax   + etalvmax)         # Vmax0 / dose (1/h)
    km_frac    <- exp(lkm     + etalkm)           # Km   / dose (unitless)
    tabs       <- exp(ltabs   + etaltabs)         # absorption half-life from intestine (min)
    tc50       <- exp(ltc50   + etaltc50)         # h
    lgtemax_i  <- lgtemax     + etalgtemax        # individual logit-transformed Emax
    hill       <- exp(lhill)                       # fixed at 10

    # Logit-style back-transform: Emax in (-1, 9) per Bulitta 2009
    # 'Parameter variability and observation model'.
    emax <- 10.0 * exp(lgtemax_i) / (1.0 + exp(lgtemax_i)) - 1.0

    # Absolute Vmax0 (mg/h) and Km (mg): the paper estimates Vmax0/
    # dose and Km/dose because the rate of stomach release is
    # primarily determined by the meal, not the cefuroxime dose
    # (Bulitta 2009 Results, 'Vmax0 and Km were estimated and are
    # reported as fractions of the cefuroxime dose'). The dose
    # amount comes from podo(stomach), so re-scaling to different
    # dose sizes preserves the dose-independent Tmax.
    dose_stomach <- podo(stomach)
    vmax0 <- vmax0_frac * dose_stomach
    km    <- km_frac    * dose_stomach

    # First-order absorption rate constant ka (1/h) derived from
    # the absorption half-life T_abs (min): footnote f.
    ka  <- log(2) * 60.0 / tabs

    # Elimination rate (1/h).
    kel <- cl / vc

    # -----------------------------------------------------------
    # Time-dependent Vmax: time-past-meal TPM coincides with time
    # since dose because the cefuroxime axetil suspension is given
    # immediately after the standardized breakfast. Sigmoidal
    # modulator scales Vmax0 by (1 + Emax * TPM^hill / (TC50^hill
    # + TPM^hill)) per Bulitta 2009 'Structural model' equation.
    # -----------------------------------------------------------
    tpm    <- tad(stomach)
    vmax_t <- vmax0 * (1.0 + emax * (tpm^hill) / (tc50^hill + tpm^hill))

    # -----------------------------------------------------------
    # ODE system. Bulitta 2009 'Structural model':
    #   dA1/dt = -Vmax(TPM) * A1 / (Km + A1)
    #   dA2/dt =  Vmax(TPM) * A1 / (Km + A1) - ka * A2
    #   dA3/dt =  ka * A2 - (CL/F) / (V/F) * A3
    # All initial conditions are zero; the stomach receives the
    # bolus dose.
    # -----------------------------------------------------------
    rel_rate        <- vmax_t * stomach / (km + stomach)
    d/dt(stomach)   <- -rel_rate
    d/dt(intestine) <-  rel_rate - ka * intestine
    d/dt(central)   <-  ka * intestine - kel * central

    # -----------------------------------------------------------
    # Observation and error model.
    # Cc has units mg/L (central in mg, vc in L).
    # -----------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
