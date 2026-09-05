Chotsiri_2024_primaquine <- function() {
  description <- "One-compartment population PK model for primaquine after a single low oral dose (0.40-0.50 mg base/kg) in 28 adult African males in Mali, roughly half of them G6PD-deficient. Absorption is a five-transit chain (ktr = 6/MTT) into a one-compartment disposition model. Fixed allometric body-weight scaling on clearance (0.75) and volume (1) referenced at WT = 62.5 kg, the cohort median. Venous plasma and capillary plasma were fitted simultaneously: the capillary prediction is the venous prediction multiplied by an estimated conversion factor CF, and each matrix carries its own exponential (additive-on-log-scale) residual error. Amounts and concentrations are molar, matching the source model, in which the same molar primaquine dose drives the parent and both metabolite models under an assumed 1:1 transformation."
  reference <- paste(
    "Chotsiri P, Mahamar A, Diawara H, Fasinu PS, Diarra K, Sanogo K,",
    "Bousema T, Walker LA, Brown JM, Dicko A, Gosling R, Chen I, Tarning J.",
    "Population pharmacokinetics of primaquine and its metabolites in African males.",
    "Malar J. 2024;23:159. doi:10.1186/s12936-024-04979-y.",
    "Structural model from Results, 'Pharmacokinetic properties of PQ and its",
    "metabolites'; parameter values from Table 2 ('Primaquine' block);",
    "allometric function and log-normal parameter model from Methods,",
    "'Population pharmacokinetic model'.",
    sep = " "
  )
  vignette <- "Chotsiri_2024_primaquine_metabolites"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline. Reference WT = 62.5 kg, 'according to median body weight in the population' (Methods, 'Population pharmacokinetic model'). Allometric exponents were fixed a priori, not estimated: 'Individual body-weight (BW_i) was introduced into the pharmacokinetic model as a fixed allometric function on all clearance (n = 0.75) and volume (n = 1.00) parameters'. Cohort weights ranged 44.0-83.0 kg across the three dose groups (Table 1). Source column BW.",
      source_name        = "BW"
    )
  )

  # Covariates the authors screened by stepwise addition (p < 0.05) and
  # elimination (p < 0.001) but did NOT retain: "All other covariates of
  # biological relevance (i.e. age, malaria status, G6PD genotypes, G6PD
  # phenotypes, and CYP2D6 phenotype) were evaluated" (Methods) and "None of
  # the other covariates were statistically significant. The G6PD genotype
  # (A-variant and wild-type) and G6PD phenotype (determined by a
  # semiquantitative test) were not statistically significant on any
  # pharmacokinetic parameters of PQ or its metabolites" (Results). No point
  # estimates are published for any of them, so none can be encoded; they are
  # documented here so the paper's covariate screen is not lost.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened, not retained. Cohort median 20-39 years by dose group, range 18-50 years (Table 1)."
    ),
    PARA = list(
      description = "Plasmodium falciparum infection detected by PCR",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened as 'malaria status', not retained. 14 of 28 participants were PCR-positive without microscopically detected parasitaemia (Table 1, 'Positive malaria PCR')."
    ),
    G6PD_A_MINUS = list(
      description = "Glucose-6-phosphate dehydrogenase A- genotype (SNPs 202A and 376G) versus wild-type",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened, not retained. 16 A- and 12 wild-type (Table 1). Documentation only: this name is not in the covariate register because no shipped model uses it."
    ),
    G6PD_DEF = list(
      description = "G6PD-deficient phenotype by semi-quantitative test (0-6.5 U/g haemoglobin)",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened, not retained. 17 deficient and 10 normal, one not reported (Table 1 and its footnote a). Documentation only: this name is not in the covariate register because no shipped model uses it."
    ),
    CYP2D6_PM = list(
      description = "CYP2D6 poor-metabolizer phenotype",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened as 'CYP2D6 phenotype', not retained. Cohort: 1 poor, 13 intermediate, 13 extensive, 1 extensive/ultra-rapid metabolizer (Table 1). The paper reports the screen only at the level of the whole phenotype covariate, so the individual indicator levels have no separate published estimates."
    )
  )

  compartmentData <- list(
    depot    = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    central  = list(analyte = "primaquine", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 28,
    n_studies      = 1,
    age_range      = "18-50 years",
    age_median     = "20 years (group 1), 32 years (group 2), 39 years (group 3)",
    weight_range   = "44.0-83.0 kg",
    weight_median  = "62.5 kg",
    sex_female_pct = 0,
    disease_state  = "Healthy adult males without microscopically detected malaria parasite infection; 17/28 G6PD-deficient and 10/28 G6PD-normal by semi-quantitative test, 16/28 carrying the G6PD A- genotype (SNPs 202A and 376G).",
    dose_range     = "Single oral dose of primaquine 0.40 mg/kg (n = 7), 0.45 mg/kg (n = 7) or 0.50 mg/kg (n = 14), i.e. 17.6-41.5 mg, given after a fatty snack as a crushed 15-mg tablet suspended in 15 mL of water.",
    regions        = "Mali (Ouelessebougou area); adult arm (Part I) of an open-label, non-randomised dose-adjustment safety trial, ClinicalTrials.gov NCT02535767.",
    notes          = "Baseline demographics from Table 1. Only the 28 adult males of the parent trial's Part I contributed pharmacokinetic samples. Sampling: venous plasma (4 mL) pre-dose and at 1, 4, 8 and 24 h; capillary plasma (0.5 mL) at 2, 4 and 6 h. All 196 samples were above the LLOQ (5 ng/mL for primaquine). Haemoglobin and methaemoglobin were followed for 28 days but were not modelled: 'None of these correlations are significantly different from the zero-slope' (Figs. 4 and 5)."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters, Table 2 'Primaquine' block. Footnote a:
    # "Computed population mean parameter estimates from NONMEM were
    # calculated for a typical individual at a body weight of 62.5 kg."
    # Apparent (per-bioavailability) parameters throughout: F is fixed
    # to 1, so CL and Vc are CL/F and Vc/F.
    # ------------------------------------------------------------------
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability F (unitless)")                                   # Table 2 'F  1 fixed'
    lmtt    <- log(0.563);    label("Mean absorption transit time MTT (h)")                                         # Table 2 'MTT (h) 0.563 (18.7% RSE)'
    lcl     <- log(15.4);     label("Apparent elimination clearance CL/F at WT = 62.5 kg (L/h)")                    # Table 2 'CL/F (L/h) 15.4 (9.57% RSE)'
    lvc     <- log(163);      label("Apparent central volume of distribution Vc/F at WT = 62.5 kg (L)")             # Table 2 'V_C/F (L) 163 (10.3% RSE)'

    # Capillary:venous conversion factor. Table 2 prints it in the '%'
    # column as 32.9 with a 28.7-37.5 confidence interval; the sibling
    # carboxy-primaquine and glucuronide rows print the same quantity as
    # a point estimate in percent with the interval already expressed as
    # a fraction (69.1 with 0.643-0.746; 40.1 with 0.353-0.465), which
    # settles that the underlying parameter is the fraction 0.329, not
    # a percentage used as-is.
    cfcap <- 0.329; label("Capillary:venous plasma conversion factor CF (unitless)")                                # Table 2 'CF (%) 32.9 (6.89% RSE)'

    # Allometric exponents, fixed a priori (not estimated).
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F (unitless)")   # Methods: "a fixed allometric function on all clearance (n = 0.75) ... parameters, centralized to 62.5 kg"
    e_wt_vc <- fixed(1);    label("Allometric exponent on Vc/F (unitless)")   # Methods: "... and volume (n = 1.00) parameters, centralized to 62.5 kg"

    # ------------------------------------------------------------------
    # Inter-individual variability. Table 2 reports IIV as %CV; the
    # footnote gives the transform used, 100 x sqrt(exp(omega^2) - 1),
    # so omega^2 = log(1 + CV^2). Each variance below is that inversion
    # of the published %CV.
    #
    # No IIV is carried on Vc/F or CF: Table 2 prints '-' in the IIV
    # column for both, and Methods states "Estimated inter-individual
    # variability below 10% was fixed to zero." They are omitted rather
    # than written as fixed(0), which would make OMEGA singular and break
    # the Cholesky sampler used by rxSolve.
    # ------------------------------------------------------------------
    etalfdepot ~ 0.2467359   # Table 2 'F  IIV 52.9% (15.4% RSE)';   log(1 + 0.529^2)
    etalmtt    ~ 0.3369643   # Table 2 'MTT IIV 63.3% (23.5% RSE)';  log(1 + 0.633^2)
    etalcl     ~ 0.0142973   # Table 2 'CL/F IIV 12% (12.4% RSE)';   log(1 + 0.120^2)

    # ------------------------------------------------------------------
    # Residual error. Methods: "Unexplained residual errors were modelled
    # separately for capillary and venous plasma concentrations and
    # implemented as an additive error model on the log-transformed
    # concentrations, equivalent to an exponential error on arithmetic
    # scale" -- that is nlmixr2's `~ lnorm(expSd)`.
    #
    # The Table 2 footnote defines both symbols as VARIANCES ("sigma_CP
    # is the variance of an exponential residual error of the capillary
    # samples, and sigma_VP is the variance ... of the venous samples"),
    # so the SD nlmixr2 wants is the square root. The published number is
    # kept visible by writing it inside sqrt().
    # ------------------------------------------------------------------
    expSd      <- sqrt(0.173); label("Residual SD on the natural-log scale, venous plasma (log units)")     # Table 2 'sigma_VP 0.173 (9.37% RSE)'
    expSd_Ccap <- sqrt(0.226); label("Residual SD on the natural-log scale, capillary plasma (log units)")  # Table 2 'sigma_CP 0.226 (9.33% RSE)'
  })

  model({
    # ---- Individual parameters ----------------------------------------
    # Methods Eq. 2: theta_i = theta_TV * exp(eta_i,theta) * (BW_i/62.5)^n
    fdepot <- exp(lfdepot + etalfdepot)
    mtt    <- exp(lmtt    + etalmtt)
    cl     <- exp(lcl     + etalcl) * (WT / 62.5)^e_wt_cl
    vc     <- exp(lvc)              * (WT / 62.5)^e_wt_vc

    # Transit-chain rate constant. Results: "Absorption of PQ and PQCG
    # were best described by 5 transit compartment models". The chain is
    # depot -> transit1 -> ... -> transit5 -> central, i.e. six equal
    # first-order transfers whose mean total time is MTT, so
    # ktr = (5 + 1)/MTT. The paper's own carboxy-primaquine model, whose
    # prose says "2 transit compartment model", pins this reading: with
    # ktr = 3/MTT it reproduces that compound's published T_MAX and
    # C_MAX, and with ktr = 2/MTT it does not.
    ktr <- 6 / mtt

    kel <- cl / vc

    # ---- ODE system (amounts in nmol, volumes in L) --------------------
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(transit3) <-  ktr * transit2 - ktr * transit3
    d/dt(transit4) <-  ktr * transit3 - ktr * transit4
    d/dt(transit5) <-  ktr * transit4 - ktr * transit5
    d/dt(central)  <-  ktr * transit5 - kel * central

    f(depot) <- fdepot

    # ---- Observations --------------------------------------------------
    # Two simultaneously fitted matrices. Results: "A linear association
    # between capillary and venous plasma concentrations of all compounds
    # was assumed and modelled using an estimated conversion factor at
    # the population level."
    Cc   <- central / vc
    Ccap <- Cc * cfcap

    Cc   ~ lnorm(expSd)
    Ccap ~ lnorm(expSd_Ccap)
  })
}
