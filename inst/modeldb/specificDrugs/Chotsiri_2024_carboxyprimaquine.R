Chotsiri_2024_carboxyprimaquine <- function() {
  description <- "One-compartment population PK model for carboxy-primaquine, the major monoamine-oxidase-A metabolite of primaquine, after a single low oral dose of primaquine (0.40-0.50 mg base/kg) in 28 adult African males in Mali. Fitted independently of the parent-drug model: the metabolite disposition compartment is driven by a two-transit absorption chain (ktr = 3/MTT) dosed with the molar primaquine amount under an assumed 1:1 transformation factor, so no fraction metabolised is identifiable and CL and Vc are apparent values relative to that molar dose. Fixed allometric body-weight scaling on clearance (0.75) and volume (1) referenced at WT = 62.5 kg. Venous and capillary plasma are fitted simultaneously through an estimated capillary:venous conversion factor, each matrix carrying its own exponential residual error."
  reference <- paste(
    "Chotsiri P, Mahamar A, Diawara H, Fasinu PS, Diarra K, Sanogo K,",
    "Bousema T, Walker LA, Brown JM, Dicko A, Gosling R, Chen I, Tarning J.",
    "Population pharmacokinetics of primaquine and its metabolites in African males.",
    "Malar J. 2024;23:159. doi:10.1186/s12936-024-04979-y.",
    "Structural model from Results, 'Pharmacokinetic properties of PQ and its",
    "metabolites'; parameter values from Table 2 ('Carboxy-primaquine' block);",
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
    depot    = list(analyte = "carboxyprimaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "carboxyprimaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "carboxyprimaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    central  = list(analyte = "carboxyprimaquine", units = "nmol", specimen = "plasma", verified = TRUE)
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
    notes          = "Baseline demographics from Table 1. Sampling: venous plasma pre-dose and at 1, 4, 8 and 24 h; capillary plasma at 2, 4 and 6 h; LLOQ 1 ng/mL for carboxy-primaquine and all 196 samples were above it. The 24-hour sampling window does not span the metabolite's terminal phase, which the authors flag explicitly as a limitation: 'the sampling schedule did not capture fully the elimination phase of CPQ' (Discussion). The resulting terminal half-life, median 469 h (325-1870), is therefore poorly determined and far longer than the 15.6 h reported in an earlier healthy-volunteer study."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters, Table 2 'Carboxy-primaquine' block, for a
    # typical individual of 62.5 kg (footnote a).
    #
    # There is no F row in this block. Relative bioavailability was fixed
    # to 1 with no inter-individual variability (Methods: "Estimated
    # inter-individual variability below 10% was fixed to zero"), so no
    # f(depot) statement appears in model() below. CL and Vc remain
    # apparent values: the dose driving this model is the molar
    # primaquine dose, not a measured carboxy-primaquine dose, because
    # the fraction of primaquine converted is not identifiable from
    # metabolite data alone -- "Pharmacokinetic models of PQ and its
    # metabolites were fitted separately due to the unknown fraction of
    # total primaquine elimination resulting in each specific metabolic
    # pathway" (Discussion).
    # ------------------------------------------------------------------
    lmtt <- log(1.24);  label("Mean absorption transit time MTT (h)")                                   # Table 2 'MTT (h) 1.24 (11.7% RSE)'
    lcl  <- log(0.129); label("Apparent elimination clearance CL/F at WT = 62.5 kg (L/h)")              # Table 2 'CL/F (L/h) 0.129 (28.5% RSE)'
    lvc  <- log(93.3);  label("Apparent central volume of distribution Vc/F at WT = 62.5 kg (L)")       # Table 2 'V_C/F (L) 93.3 (7.00% RSE)'

    # Capillary:venous conversion factor. Table 2's '%' column carries
    # 69.1 while its confidence interval is printed as the fraction
    # 0.643-0.746, which settles the parameter as the fraction 0.691.
    cfcap <- 0.691; label("Capillary:venous plasma conversion factor CF (unitless)")                    # Table 2 'CF (%) 69.1 (3.83% RSE)'

    # Allometric exponents, fixed a priori (not estimated).
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F (unitless)")   # Methods: "a fixed allometric function on all clearance (n = 0.75) ... parameters, centralized to 62.5 kg"
    e_wt_vc <- fixed(1);    label("Allometric exponent on Vc/F (unitless)")   # Methods: "... and volume (n = 1.00) parameters, centralized to 62.5 kg"

    # ------------------------------------------------------------------
    # Inter-individual variability, inverted from Table 2's %CV column
    # through the published transform 100 x sqrt(exp(omega^2) - 1), i.e.
    # omega^2 = log(1 + CV^2). Table 2 prints '-' for CL/F and CF, whose
    # IIV was estimated below 10% and fixed to zero; those etas are
    # omitted rather than written as fixed(0) to keep OMEGA non-singular.
    # ------------------------------------------------------------------
    etalmtt ~ 0.3551601   # Table 2 'MTT IIV 65.3% (13.2% RSE)';   log(1 + 0.653^2)
    etalvc  ~ 0.1309195   # Table 2 'V_C/F IIV 37.4% (15.7% RSE)'; log(1 + 0.374^2)

    # ------------------------------------------------------------------
    # Residual error: additive on the log scale, i.e. `~ lnorm()`. The
    # Table 2 footnote defines sigma_VP and sigma_CP as VARIANCES, so the
    # SD is the square root; the published number stays visible inside
    # sqrt().
    # ------------------------------------------------------------------
    expSd      <- sqrt(0.0328); label("Residual SD on the natural-log scale, venous plasma (log units)")     # Table 2 'sigma_VP 0.0328 (8.96% RSE)'
    expSd_Ccap <- sqrt(0.101);  label("Residual SD on the natural-log scale, capillary plasma (log units)")  # Table 2 'sigma_CP 0.101 (8.84% RSE)'
  })

  model({
    # Methods Eq. 2: theta_i = theta_TV * exp(eta_i,theta) * (BW_i/62.5)^n
    mtt <- exp(lmtt + etalmtt)
    cl  <- exp(lcl)             * (WT / 62.5)^e_wt_cl
    vc  <- exp(lvc  + etalvc)   * (WT / 62.5)^e_wt_vc

    # Results: "the absorption of CPQ was best described by a 2 transit
    # compartment model". The chain depot -> transit1 -> transit2 ->
    # central is three equal first-order transfers of mean total time
    # MTT, so ktr = (2 + 1)/MTT. This reading reproduces the published
    # secondary estimates for this compound: T_MAX 4.8 h against a
    # reported median 4.57 h and C_MAX 1233 nmol/L against a reported
    # median 338 ng/mL = 1232 nmol/L at MW 274.32 g/mol.
    ktr <- 3 / mtt

    kel <- cl / vc

    # ---- ODE system (amounts in nmol, volumes in L) --------------------
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(central)  <-  ktr * transit2 - kel * central

    # ---- Observations --------------------------------------------------
    Cc   <- central / vc
    Ccap <- Cc * cfcap

    Cc   ~ lnorm(expSd)
    Ccap ~ lnorm(expSd_Ccap)
  })
}
