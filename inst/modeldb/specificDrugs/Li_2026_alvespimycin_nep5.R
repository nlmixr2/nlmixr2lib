Li_2026_alvespimycin_nep5 <- function() {
  description <- "One-compartment population PK model for the heat shock protein 90 inhibitor 17-DMAG (alvespimycin) given as an IV infusion to adult patients with advanced solid tumors, recovered as the most parsimonious non-dominated solution (5 estimated parameters, OFV 9813.4) on the NSGA-II multi-objective Pareto front of Li 2026, with first-order elimination, log-normal IIV on CL and Vc, and an exponential residual error model. This is a machine-search model structure from a multi-objective model-selection study, not an expert-developed final model."
  reference <- paste(
    "Li X, Sale M, Craig J, Nieforth K, Mazur A, Bies RR. (2026).",
    "Multi-objective optimization in population pharmacokinetic model",
    "selection and optimization: application of NSGA-II in pyDarwin.",
    "J Pharmacokinet Pharmacodyn 53(1):26.",
    "doi:10.1007/s10928-026-10036-9.",
    "Parameter values from Supplementary Table S2 (Model a, NEP = 5).",
    "Structural parameterization (occasion structure, residual-error form,",
    "weight centering) is read from the pyDarwin template and token files",
    "for this same DMAG search space, deposited as Supplementary Material 1,",
    "2 and 4 of Li X, Sale M, Nieforth K, Craig J, Wang F, Solit D, Feng K,",
    "Hu M, Bies RR, Zhao L. (2024). pyDarwin machine learning algorithms",
    "application and comparison in nonlinear mixed-effect model selection",
    "and optimization. J Pharmacokinet Pharmacodyn 51(6):785-796.",
    "doi:10.1007/s10928-024-09932-9.",
    sep = " "
  )
  vignette <- "Li_2026_alvespimycin"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "alvespimycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # Model a retained no covariate effects (Li 2026 Supplementary Table S3, DMAG
  # NSGA-II rows, "Covariate" column is `-` at NEP = 5).
  covariateData <- list()

  # Covariates that the NSGA-II search space screened for the DMAG data set but
  # that this particular non-dominated solution did not retain (Li 2026 Table 3
  # "Searched covariates" for DMAG; the pyDarwin token file of Li 2024
  # Supplementary Material 2 enumerates the same set as `V~WT`, `CL~WT`,
  # `Q~WT`, `Q2~WT`, `V2~WT`, `V3~WT`, `V~AGE`, `CL~AGE`, `V~SEX`, `CL~SEX`,
  # `CL~SCR`). Documentation only; not referenced in `model()`.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = "In the search space as a power function of weight centred on 81 kg on Vc, CL, Q2, Q3, V2 and V3; not retained at NEP = 5."
    ),
    AGE = list(
      description = "Age.",
      units       = "years",
      type        = "continuous",
      notes       = "In the search space as a power function of age centred on 60 years on Vc and CL; not retained at NEP = 5."
    ),
    SEXF = list(
      description = "Sex indicator.",
      units       = "(binary)",
      type        = "categorical",
      notes       = "In the search space as an exponential effect of the NONMEM SEX column on Vc and CL; not retained at NEP = 5. The Li 2024 control stream codes the column as SEX with a 0/1 encoding whose reference level is not stated in either publication."
    ),
    CREAT = list(
      description = "Serum creatinine.",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "In the search space as a power function of serum creatinine centred on 0.9 mg/dL on CL (Li 2024 control stream `CSCR = SCR/0.9`); not retained at NEP = 5."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 66L,
    n_observations = 951L,
    n_studies      = 1L,
    age_range      = "28-82 years (median 63)",
    age_median     = "63 years",
    weight_range   = "48.2-136.5 kg (median 80.3)",
    weight_median  = "80.3 kg",
    sex_female_pct = 38,
    disease_state  = "Adult patients with advanced solid tumors.",
    dose_range     = "IV infusion of 17-DMAG; median 33, range 2.2-413 per the dose column of Li 2026 Table 1 (see notes on the unit of that column).",
    notes          = "Demographics from Li 2026 Supplementary Table S1 (DMAG panel): age 63 (28-82) years, weight 80.3 (48.2-136.5) kg, blood urea nitrogen 15 (5-70) mg/dL, creatinine 1 (0.6-1.8) mg/dL, 41 (62%) male and 25 (38%) female. Sampling was relatively rich, averaging more than 15 samples per participant over up to 102 h. This is the same 17-DMAG data set analysed by the stepwise model of Aregbe 2012 (registered here as `Aregbe_2012_alvespimycin`), which reports the per-dose range as 2.2-413 mg/m^2; Li 2026 Table 1 prints the same two numbers labelled 'mg'."
  )

  ini({
    # Structural PK parameters from Li 2026 Supplementary Table S2, Model a
    # (NEP = 5, OFV = 9813.4). NONMEM ADVAN1 / TRANS2 (Li 2026 Table S3, DMAG
    # rows, "COM" = 1 at NEP = 5); the ADVAN1 branch of the Li 2024 token file
    # parameterizes the single compartment as CL and V.
    lcl <- log(8.41) ; label("Clearance (L/h)")                      # Li 2026 Table S2 Model a: CL = 8.41 L/hr (RSE 5%)
    lvc <- log(98.4) ; label("Central volume of distribution (L)")   # Li 2026 Table S2 Model a: V  = 98.4 L    (RSE 4%)

    # Between-subject variability. Li 2026 Table S2 reports BSV as a CV%; the
    # internal log-normal variance is omega^2 = log(CV^2 + 1). BSV was retained
    # on V and CL only (Li 2026 Table S3, "BSV" column reads "V, CL").
    etalcl ~ 0.120592  # Li 2026 Table S2 Model a: BSV CL 35.8% (RSE 10), variance = log(0.358^2 + 1)
    etalvc ~ 0.093462  # Li 2026 Table S2 Model a: BSV V  31.3% (RSE 10), variance = log(0.313^2 + 1)

    # Residual error. Li 2026 Table S3 labels the NEP = 5 residual model
    # "prop", which is the RESERR option 0 of the deposited pyDarwin token file
    # (Li 2024 Supplementary Material 2): `IOBS = F*EXP(EPS(1))`. That is a
    # log-normal (exponential) residual, so it maps to `lnorm(expSd)` rather
    # than `prop(propSd)`. The tabulated 0.214 is the NONMEM $SIGMA variance
    # (see the vignette Errata for the cross-check that fixes the scale), so
    # expSd = sqrt(0.214).
    expSd <- 0.462601 ; label("Exponential (log-normal) residual error SD")  # Li 2026 Table S2 Model a: proportional error 0.214 (RSE 6); SD = sqrt(0.214)
  })

  model({
    # Individual PK parameters; no covariate effects at NEP = 5.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # One-compartment IV PK with first-order elimination (NONMEM ADVAN1 /
    # TRANS2). Dose lands in `central` via the data set's cmt column; the
    # infusion duration is encoded on the dose record (rate or dur).
    d/dt(central) <- -cl / vc * central

    # Plasma concentration. `central` is in mg and `vc` in L, so central / vc
    # is mg/L; the factor of 1000 converts to the ng/mL scale on which the
    # 17-DMAG assay data were modelled (the Li 2024 control stream sets
    # S1 = V1 with no scaling, and the prediction-corrected VPC panels of
    # Li 2026 Figure 2 run to roughly 2500 on the concentration axis).
    Cc <- 1000 * central / vc

    # Log-normal residual error on the concentration scale.
    Cc ~ lnorm(expSd)
  })
}
