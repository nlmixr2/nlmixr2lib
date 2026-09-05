Chotsiri_2024_primaquineCarbamoylGlucuronide <- function() {
  description <- "One-compartment population PK model for primaquine carbamoyl-glucuronide, a phase-II metabolite of primaquine whose pharmacokinetics are reported here for the first time, after a single low oral dose of primaquine (0.40-0.50 mg base/kg) in 28 adult African males in Mali. Fitted independently of the parent-drug model: a five-transit absorption chain (ktr = 6/MTT) feeds the metabolite disposition compartment, dosed with the molar primaquine amount under an assumed 1:1 transformation factor, so no fraction metabolised is identifiable and CL and Vc are apparent values relative to that molar dose. Fixed allometric body-weight scaling on clearance (0.75) and volume (1) referenced at WT = 62.5 kg. Venous and capillary plasma are fitted simultaneously through an estimated capillary:venous conversion factor, each matrix carrying its own exponential residual error."
  reference <- paste(
    "Chotsiri P, Mahamar A, Diawara H, Fasinu PS, Diarra K, Sanogo K,",
    "Bousema T, Walker LA, Brown JM, Dicko A, Gosling R, Chen I, Tarning J.",
    "Population pharmacokinetics of primaquine and its metabolites in African males.",
    "Malar J. 2024;23:159. doi:10.1186/s12936-024-04979-y.",
    "Structural model from Results, 'Pharmacokinetic properties of PQ and its",
    "metabolites'; parameter values from Table 2 ('Primaquine",
    "carbamoyl-glucuronide' block); allometric function and log-normal",
    "parameter model from Methods, 'Population pharmacokinetic model'.",
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

  # Screened by stepwise addition (p < 0.05) / elimination (p < 0.001) and
  # not retained for any compound: "None of the other covariates were
  # statistically significant" (Results). No point estimates are published,
  # so none can be encoded. See Chotsiri_2024_primaquine.R for the same list.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened, not retained. Range 18-50 years (Table 1)."
    ),
    PARA = list(
      description = "Plasmodium falciparum infection detected by PCR",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened as 'malaria status', not retained. 14 of 28 participants PCR-positive (Table 1)."
    ),
    G6PD_A_MINUS = list(
      description = "Glucose-6-phosphate dehydrogenase A- genotype (SNPs 202A and 376G) versus wild-type",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened, not retained. 16 A- and 12 wild-type (Table 1). Documentation only: not in the covariate register because no shipped model uses it."
    ),
    G6PD_DEF = list(
      description = "G6PD-deficient phenotype by semi-quantitative test (0-6.5 U/g haemoglobin)",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened, not retained. 17 deficient and 10 normal (Table 1). Documentation only: not in the covariate register because no shipped model uses it."
    ),
    CYP2D6_PM = list(
      description = "CYP2D6 poor-metabolizer phenotype",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened as 'CYP2D6 phenotype', not retained. 1 poor, 13 intermediate, 13 extensive, 1 extensive/ultra-rapid metabolizer (Table 1)."
    )
  )

  compartmentData <- list(
    depot    = list(analyte = "primaquine carbamoyl-glucuronide", units = "nmol", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "primaquine carbamoyl-glucuronide", units = "nmol", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "primaquine carbamoyl-glucuronide", units = "nmol", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "primaquine carbamoyl-glucuronide", units = "nmol", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "primaquine carbamoyl-glucuronide", units = "nmol", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "primaquine carbamoyl-glucuronide", units = "nmol", specimen = "administration site", verified = TRUE),
    central  = list(analyte = "primaquine carbamoyl-glucuronide", units = "nmol", specimen = "plasma", verified = TRUE)
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
    dose_range     = "Single oral dose of primaquine 0.40 mg/kg (n = 7), 0.45 mg/kg (n = 7) or 0.50 mg/kg (n = 14), i.e. 17.6-41.5 mg. The metabolite model is dosed with the same amount expressed in moles.",
    regions        = "Mali (Ouelessebougou area); adult arm (Part I) of an open-label, non-randomised dose-adjustment safety trial, ClinicalTrials.gov NCT02535767.",
    notes          = "Baseline demographics from Table 1. Sampling: venous plasma pre-dose and at 1, 4, 8 and 24 h; capillary plasma at 2, 4 and 6 h; LLOQ 5 ng/mL for primaquine carbamoyl-glucuronide and all 196 samples were above it. 'Pharmacokinetic properties of the carbamoyl-glucuronide metabolite are reported here for the first time' (Discussion)."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters, Table 2 'Primaquine carbamoyl-glucuronide'
    # block, for a typical individual of 62.5 kg (footnote a).
    #
    # CL and Vc are apparent values relative to the molar primaquine
    # dose: the fraction of primaquine converted to this metabolite is
    # not identifiable from metabolite data alone -- "Pharmacokinetic
    # models of PQ and its metabolites were fitted separately due to the
    # unknown fraction of total primaquine elimination resulting in each
    # specific metabolic pathway" (Discussion).
    # ------------------------------------------------------------------
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability F (unitless)")                                   # Table 2 'F  1 fixed'
    lmtt    <- log(1.13);     label("Mean absorption transit time MTT (h)")                                         # Table 2 'MTT (h) 1.13 (8.84% RSE)'
    lcl     <- log(2.83);     label("Apparent elimination clearance CL/F at WT = 62.5 kg (L/h)")                    # Table 2 'CL/F (L/h) 2.83 (15.8% RSE)'
    lvc     <- log(55.4);     label("Apparent central volume of distribution Vc/F at WT = 62.5 kg (L)")             # Table 2 'V_C/F (L) 55.4 (12.4% RSE)'

    # Capillary:venous conversion factor. Table 2's '%' column carries
    # 40.1 while its confidence interval is printed as the fraction
    # 0.353-0.465, which settles the parameter as the fraction 0.401.
    cfcap <- 0.401; label("Capillary:venous plasma conversion factor CF (unitless)")                                # Table 2 'CF (%) 40.1 (6.97% RSE)'

    # Allometric exponents, fixed a priori (not estimated).
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F (unitless)")   # Methods: "a fixed allometric function on all clearance (n = 0.75) ... parameters, centralized to 62.5 kg"
    e_wt_vc <- fixed(1);    label("Allometric exponent on Vc/F (unitless)")   # Methods: "... and volume (n = 1.00) parameters, centralized to 62.5 kg"

    # ------------------------------------------------------------------
    # Inter-individual variability, inverted from Table 2's %CV column
    # through the published transform 100 x sqrt(exp(omega^2) - 1), i.e.
    # omega^2 = log(1 + CV^2). Table 2 prints '-' for Vc/F and CF, whose
    # IIV was estimated below 10% and fixed to zero; those etas are
    # omitted rather than written as fixed(0) to keep OMEGA non-singular.
    # ------------------------------------------------------------------
    etalfdepot ~ 0.3387732   # Table 2 'F  IIV 63.5% (12.1% RSE)';   log(1 + 0.635^2)
    etalmtt    ~ 0.1118419   # Table 2 'MTT IIV 34.4% (25.1% RSE)';  log(1 + 0.344^2)
    etalcl     ~ 0.2882449   # Table 2 'CL/F IIV 57.8% (12.1% RSE)'; log(1 + 0.578^2)

    # ------------------------------------------------------------------
    # Residual error: additive on the log scale, i.e. `~ lnorm()`. The
    # Table 2 footnote defines sigma_VP and sigma_CP as VARIANCES, so the
    # SD is the square root; the published number stays visible inside
    # sqrt().
    # ------------------------------------------------------------------
    expSd      <- sqrt(0.108); label("Residual SD on the natural-log scale, venous plasma (log units)")     # Table 2 'sigma_VP 0.108 (12.3% RSE)'
    expSd_Ccap <- sqrt(0.242); label("Residual SD on the natural-log scale, capillary plasma (log units)")  # Table 2 'sigma_CP 0.242 (8.00% RSE)'
  })

  model({
    # Methods Eq. 2: theta_i = theta_TV * exp(eta_i,theta) * (BW_i/62.5)^n
    fdepot <- exp(lfdepot + etalfdepot)
    mtt    <- exp(lmtt    + etalmtt)
    cl     <- exp(lcl     + etalcl) * (WT / 62.5)^e_wt_cl
    vc     <- exp(lvc)              * (WT / 62.5)^e_wt_vc

    # Results: "Absorption of PQ and PQCG were best described by 5 transit
    # compartment models". Six equal first-order transfers of mean total
    # time MTT, so ktr = (5 + 1)/MTT -- the same reading validated on the
    # paper's carboxy-primaquine model, whose prose transit count and
    # published secondary estimates agree only under ktr = (n + 1)/MTT.
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
    Cc   <- central / vc
    Ccap <- Cc * cfcap

    Cc   ~ lnorm(expSd)
    Ccap ~ lnorm(expSd_Ccap)
  })
}
