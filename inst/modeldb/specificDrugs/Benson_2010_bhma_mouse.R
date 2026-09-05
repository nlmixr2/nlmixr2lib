Benson_2010_bhma_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, male CD-1 outbred).",
    "Sequential PK-PD model for the Toll-like receptor 7 (TLR-7) agonist",
    "BHMA [9-benzyl-8-hydroxy-2-(2-methoxyethoxy) adenine] with plasma",
    "alpha interferon (IFN-alpha) as the pharmacodynamic biomarker",
    "(Benson 2010). Stage 1: two-compartment apparent oral PK with",
    "first-order absorption (NONMEM ADVAN4), all disposition parameters",
    "reported per kilogram body weight; the absorption rate constant was",
    "fixed at 35 1/h because absorption was complete before the earliest",
    "PK sample at 0.1 h, and bioavailability could not be estimated",
    "without intravenous data, so CL, Vc, Q and Vp are apparent (/F)",
    "values. Stage 2: indirect-response (turnover) model for plasma",
    "IFN-alpha in which the drug STIMULATES synthesis through an ordinary",
    "Emax (non-sigmoid) function of the predicted plasma BHMA",
    "concentration and IFN-alpha is removed by a first-order rate",
    "constant, d(IFN)/dt = Smax * Cc / (SC50 + Cc) - kout * IFN. The",
    "IFN-alpha baseline is fixed at zero because unstimulated levels were",
    "below the 0.3 IU/mL assay LOQ in every animal. The Emax form beat a",
    "linear drug-effect alternative by dAIC = -22.9. Parameters from",
    "Benson 2010 Table 1 (PK) and Table 2 (PK-PD); no covariates were",
    "retained in either stage."
  )

  reference <- paste(
    "Benson N, de Jongh J, Duckworth JD, Jones HM, Pertinez HE, Rawal JK,",
    "van Steeg TJ, Van der Graaf PH. Pharmacokinetic-pharmacodynamic",
    "modeling of alpha interferon response induced by a Toll-like 7",
    "receptor agonist in mice.",
    "Antimicrob Agents Chemother. 2010;54(3):1179-1185.",
    "doi:10.1128/AAC.00551-09"
  )

  vignette <- "Benson_2010_bhma_mouse"

  units <- list(
    time          = "h",
    dosing        = "mg/kg",
    concentration = "ng/mL (BHMA in plasma); IU/mL (IFN-alpha in plasma)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. All four entries were checked against Benson 2010
  # (Materials and Methods, 'In vivo experiments with BHMA' and 'PK and
  # PK-PD modeling'; Fig. 1 schematic).
  compartmentData <- list(
    depot       = list(analyte = "BHMA",      units = "mg/kg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "BHMA",      units = "mg/kg", specimen = "plasma",              verified = TRUE),
    peripheral1 = list(analyte = "BHMA",      units = "mg/kg", specimen = "plasma",              verified = TRUE),
    ifna        = list(analyte = "IFN-alpha", units = "IU/mL", specimen = "plasma",              verified = TRUE)
  )

  # No covariates were tested or retained in either stage of the Benson
  # 2010 analysis; the mice were a single male CD-1 cohort dosed on a
  # mg/kg basis, so body weight is absorbed into the per-kg
  # parameterisation rather than entering as a covariate.
  covariateData <- list()

  population <- list(
    species        = "mouse (male CD-1 outbred)",
    n_subjects     = 71L,
    n_studies      = 1L,
    disease_state  = paste(
      "Healthy, uninfected male CD-1 outbred mice (Charles River",
      "Laboratories UK, Margate, Kent). No disease model was used;",
      "IFN-alpha induction was measured in naive animals."
    ),
    dose_range     = paste(
      "BHMA administered orally in 5% dimethyl sulfoxide / 50%",
      "polyethylene glycol 200 / 45% water. PD groups (n per group in",
      "parentheses): 0.1 mg/kg (3), 0.3 mg/kg (2), 0.5 mg/kg (2),",
      "1 mg/kg (10), 2.5 mg/kg (4), 5 mg/kg (20), 10 mg/kg (5) -- 46",
      "animals total. PK-only groups: 0.1, 1 and 5 mg/kg, 25 animals",
      "total. Single doses only."
    ),
    regions        = "United Kingdom (Pfizer Global Research and Development, Sandwich).",
    notes          = paste(
      "Benson 2010 Materials and Methods, 'In vivo experiments with",
      "BHMA'. A matrix (composite) sampling design was used: two 50 uL",
      "saphenous-vein samples plus a terminal 1 mL sample per animal,",
      "so PK and biomarker profiles are composite rather than fully",
      "serial. PK-only sampling times were 0.1, 0.25, 1, 2, 6, 8, 10 and",
      "24 h at 5 mg/kg and 1, 2, 6, 8 and 10 h at 1 and 0.1 mg/kg.",
      "BHMA was quantified by LC-MS/MS (LOQ 0.5 ng/mL, assay CV 20%) and",
      "IFN-alpha by a modified PBL Biomedical ELISA (LOQ 0.3 IU/mL,",
      "assay CV 20%; 1 IU/mL was taken to equal 15 pg/mL). The plasma",
      "free fraction of BHMA in mouse was 0.33 (SE 0.006) by equilibrium",
      "dialysis, which the paper uses to convert the total SC50 of",
      "135 ng/mL to an unbound potency of ca. 125 nM. Estimation was in",
      "NONMEM VI with FOCE-interaction; ADVAN4 for the PK stage and",
      "ADVAN6 for the PK-PD stage, fitted sequentially with the PK",
      "parameters fixed at their typical population values."
    )
  )

  ini({
    # ==================================================================
    # Stage 1 -- apparent two-compartment oral PK (Benson 2010 Table 1).
    # Every disposition parameter is reported per kilogram body weight;
    # the Discussion states that without intravenous data "the
    # bioavailability (F) cannot be estimated, and hence, the volume of
    # distribution (V) of ca. 40 liters/kg must be regarded as the V/F
    # ratio", so all four are apparent (/F) values.
    # ==================================================================
    lka <- fixed(log(35))
    label("Apparent first-order oral absorption rate constant (ka, 1/h)")
    # Benson 2010 Table 1 row 3: Ka = 35 h^-1, marked "Fix". Results,
    # "PK of BHMA in mice": "The absorption rate constant was fixed at
    # 35 h-1, reflecting complete absorption prior to the collection of
    # the earliest PK samples (at 0.1 h)."

    lcl <- log(3.45)
    label("Apparent clearance per kg body weight (CL/F, L/h/kg)")
    # Benson 2010 Table 1 row 1: CL = 3.45 liters/h/kg, CV 16%.

    lvc <- log(4.31)
    label("Apparent central volume of distribution per kg body weight (V2/F, L/kg)")
    # Benson 2010 Table 1 row 2: V2 = 4.31 liters/kg, CV 24%.

    lq <- log(2.78)
    label("Apparent intercompartmental clearance per kg body weight (Q/F, L/h/kg)")
    # Benson 2010 Table 1 row 4: Q = 2.78 liters/kg/h, CV 22%.

    lvp <- log(39.7)
    label("Apparent peripheral volume of distribution per kg body weight (V3/F, L/kg)")
    # Benson 2010 Table 1 row 5: V3 = 39.7 liters/kg, CV 36%. The
    # Discussion refers to this as the "ca. 40 liters/kg" V/F ratio.

    # ==================================================================
    # Stage 2 -- IFN-alpha indirect response (Benson 2010 Table 2,
    # equation 1). Zero baseline: Results, "PK and PD of IFN induction
    # by BHMA in mice" -- "in the absence of a TLR-7 agonist, the basal
    # level of IFN-alpha was, in all cases, below the LOQ of our assay
    # (0.3 IU/ml), and hence it was fixed at zero for the purpose of
    # PK-PD modeling."
    # ==================================================================
    lsmax <- log(294)
    label("Maximum zero-order IFN-alpha synthesis rate (Smax, IU/mL/h)")
    # Benson 2010 Table 2 row 3: Smax = 294 IU/ml/h, CV 8%. The
    # Discussion cross-checks this against kout: "the maximum response
    # is given by the Smax/kout ratio of 306 IU/ml"
    # (294 / 0.958 = 306.9, consistent).

    lec50 <- log(135)
    label("Total plasma BHMA concentration giving half-maximal IFN-alpha synthesis (SC50, ng/mL)")
    # Benson 2010 Table 2 row 2: SC50 = 135 ng/ml, CV 24%. This is a
    # TOTAL plasma concentration; the Results note that "corrected for
    # the free fraction" (fu = 0.33) it "equates to 125 nM" unbound.

    lkout <- log(0.958)
    label("First-order IFN-alpha elimination rate constant (kout, 1/h)")
    # Benson 2010 Table 2 row 1: kout = 0.958 h^-1, CV 0.1%. The
    # Discussion converts this to a half-life: "The apparent degradation
    # rate constant for IFN-alpha was 0.96 h-1, giving a half-life
    # (t1/2) of approximately 0.7 h" (log(2)/0.958 = 0.723 h).

    # ==================================================================
    # Interindividual variability. Benson 2010 equation 3 is the
    # multiplicative exponential model Pi = theta * exp(eta_i), and the
    # Methods state "The derived IIV is expressed as the percent CV", so
    # the tabulated percentages are lognormal CVs and the variances are
    # recovered as omega^2 = log(CV^2 + 1).
    # ==================================================================
    etalcl ~ 0.134880
    # Benson 2010 Table 1 row 6: IIV of CL = 38% (CV of the estimate
    # 38%). omega^2 = log(1 + 0.38^2) = 0.134880.

    etalsmax ~ 0.398776
    # Benson 2010 Table 2 row 4: IIV of Smax = 70% (CV of the estimate
    # 26%). omega^2 = log(1 + 0.70^2) = 0.398776. Results: "IIV of the
    # Smax was high at 70%"; the Discussion attributes it to a genuine
    # interindividual PD difference rather than assay variability (the
    # IFN-alpha assay CV was 20%). No IIV was retained on kout or SC50.
    # ==================================================================
    # Residual error. Benson 2010 gives both an additive (equation 4)
    # and a proportional (equation 5) residual model; each stage uses
    # one of them, identified by the units in which the table reports
    # the estimate.
    # ==================================================================
    propSd <- 0.46
    label("Proportional residual SD on plasma BHMA concentration (fraction)")
    # Benson 2010 Table 1 row 7: "Residual error (%) 46", CV 32%. The
    # percentage unit identifies the proportional model of equation 5.

    addSd_ifna <- 65
    label("Additive residual SD on plasma IFN-alpha concentration (IU/mL)")
    # Benson 2010 Table 2 row 5: "Residual error (IU/ml) 65", CV 19%.
    # The IU/mL unit identifies the additive model of equation 4;
    # Results: "the residual error was 65 IU/ml".
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters. All are per kilogram body weight, so
    #    doses must be supplied in mg/kg. Only CL carries IIV (Table 1).
    # ------------------------------------------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # ------------------------------------------------------------------
    # 2. Individual PD parameters. Only Smax carries IIV (Table 2).
    # ------------------------------------------------------------------
    smax <- exp(lsmax + etalsmax)
    ec50 <- exp(lec50)
    kout <- exp(lkout)

    # ------------------------------------------------------------------
    # 3. Micro-constants for the two-compartment disposition (ADVAN4).
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ------------------------------------------------------------------
    # 4. BHMA PK. Amounts are in mg/kg because volumes are in L/kg.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1

    # ------------------------------------------------------------------
    # 5. Predicted plasma BHMA concentration. central/vc is mg/kg per
    #    L/kg = mg/L = ug/mL; the factor 1000 converts to ng/mL, the
    #    units in which Benson 2010 reports concentrations (LOQ
    #    0.5 ng/mL) and in which SC50 = 135 ng/mL was estimated.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc

    # ------------------------------------------------------------------
    # 6. IFN-alpha indirect response, Benson 2010 equation 1:
    #
    #        dR      Smax * Cp
    #        --  =  -----------  -  kout * R
    #        dt      SC50 + Cp
    #
    #    An ordinary (non-sigmoid) Emax stimulation of synthesis with no
    #    Hill exponent -- Table 2 lists none, and the Results report that
    #    "comparing the ordinary Emax model with a linear model gave a
    #    dAIC of -22.9 and that, therefore, the ordinary Emax model had
    #    the superior fit". The baseline is fixed at zero (see ini()),
    #    which is why the production term is a bare Smax * Cc / (SC50 +
    #    Cc) rather than the usual kin * (1 + stimulation).
    # ------------------------------------------------------------------
    d/dt(ifna) <- smax * Cc / (ec50 + Cc) - kout * ifna
    ifna(0)    <- 0

    # ------------------------------------------------------------------
    # 7. Two observations: plasma BHMA (proportional error, equation 5)
    #    and plasma IFN-alpha (additive error, equation 4).
    # ------------------------------------------------------------------
    Cc   ~ prop(propSd)
    ifna ~ add(addSd_ifna)
  })
}
