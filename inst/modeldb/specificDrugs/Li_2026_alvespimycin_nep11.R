Li_2026_alvespimycin_nep11 <- function() {
  description <- "Three-compartment population PK model for the heat shock protein 90 inhibitor 17-DMAG (alvespimycin) given as an IV infusion to adult patients with advanced solid tumors, recovered as the 11-estimated-parameter non-dominated solution (OFV 8166.5) on the NSGA-II multi-objective Pareto front of Li 2026, with first-order elimination, log-normal IIV on CL, Vc and Vp, between-occasion variability on clearance, and an exponential residual error model. This is a machine-search model structure from a multi-objective model-selection study, not an expert-developed final model."
  reference <- paste(
    "Li X, Sale M, Craig J, Nieforth K, Mazur A, Bies RR. (2026).",
    "Multi-objective optimization in population pharmacokinetic model",
    "selection and optimization: application of NSGA-II in pyDarwin.",
    "J Pharmacokinet Pharmacodyn 53(1):26.",
    "doi:10.1007/s10928-026-10036-9.",
    "Parameter values from Supplementary Table S2 (Model c, NEP = 11).",
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

  compartmentData <- list(
    central     = list(analyte = "alvespimycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "alvespimycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "alvespimycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    OCC = list(
      description        = "Integer-valued occasion indicator for between-occasion variability multiplexing.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Values 1, 2, 3 identify the dosing occasion within subject. The deposited pyDarwin control stream for this data set (Li 2024 Supplementary Material 4, `$PK`) writes the BOV blocks as `IF(OCC.EQ.1) ... IF(OCC.EQ.2) ... IF(OCC.EQ.3)`, so three occasions are carried. Decomposed inside `model()` into binary indicators `oc1`, `oc2`, `oc3` that multiplex the BOV etas.",
      source_name        = "OCC"
    )
  )

  # Screened but not retained at NEP = 11 (Li 2026 Table 3 "Searched
  # covariates" for DMAG; Li 2026 Table S3 shows an empty "Covariate" cell at
  # NEP = 11).
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = "In the search space as a power function of weight centred on 81 kg on Vc, CL, Q2, Q3, V2, V3; not retained at NEP = 11."
    ),
    AGE = list(
      description = "Age.",
      units       = "years",
      type        = "continuous",
      notes       = "In the search space as a power function of age centred on 60 years on Vc and CL; not retained at NEP = 11."
    ),
    SEXF = list(
      description = "Sex indicator.",
      units       = "(binary)",
      type        = "categorical",
      notes       = "In the search space as an exponential effect of the NONMEM SEX column on Vc and CL; not retained at NEP = 11. The Li 2024 control stream codes the column as SEX with a 0/1 encoding whose reference level is not stated in either publication."
    ),
    CREAT = list(
      description = "Serum creatinine.",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "In the search space as a power function of serum creatinine centred on 0.9 mg/dL on CL (Li 2024 control stream `CSCR = SCR/0.9`); not retained at NEP = 11."
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
    # Structural PK parameters from Li 2026 Supplementary Table S2, Model c
    # (NEP = 11, OFV = 8166.5). NONMEM ADVAN11 / TRANS4 (Li 2026 Table S3,
    # DMAG rows, "COM" = 3 at NEP = 11). NONMEM names map to nlmixr2
    # conventions as CL -> lcl, V -> lvc, Q2 -> lq, V2 -> lvp, Q3 -> lq2,
    # V3 -> lvp2.
    lcl  <- log(9.26) ; label("Clearance (L/h)")                                  # Li 2026 Table S2 Model c: CL = 9.26 L/hr (RSE 6%)
    lvc  <- log(28.1) ; label("Central volume of distribution (L)")               # Li 2026 Table S2 Model c: V  = 28.1 L    (RSE 10%)
    lq   <- log(81.8) ; label("Inter-compartmental clearance Q2 (L/h)")           # Li 2026 Table S2 Model c: Q2 = 81.8 L/hr (RSE 4%)
    lvp  <- log(75.6) ; label("Peripheral volume of distribution V2 (L)")         # Li 2026 Table S2 Model c: V2 = 75.6 L    (RSE 9%)
    lq2  <- log(9.47) ; label("Inter-compartmental clearance Q3 (L/h)")           # Li 2026 Table S2 Model c: Q3 = 9.47 L/hr (RSE 8%)
    lvp2 <- log(85.3) ; label("Peripheral volume of distribution V3 (L)")         # Li 2026 Table S2 Model c: V3 = 85.3 L    (RSE 7%)

    # Between-subject variability. Li 2026 Table S2 reports BSV as a CV%; the
    # internal log-normal variance is omega^2 = log(CV^2 + 1). BSV was retained
    # on V, CL and V2 (Li 2026 Table S3, "BSV" column reads "V, CL, V2").
    etalcl ~ 0.173302  # Li 2026 Table S2 Model c: BSV CL 43.5% (RSE 11), variance = log(0.435^2 + 1)
    etalvc ~ 0.400656  # Li 2026 Table S2 Model c: BSV V  70.2% (RSE 11), variance = log(0.702^2 + 1)
    etalvp ~ 0.336964  # Li 2026 Table S2 Model c: BSV V2 63.3% (RSE 10), variance = log(0.633^2 + 1)

    # Between-occasion variability on CL (Li 2026 Table S3, "BOV" column reads
    # "CL" at NEP = 11; Table S2 puts the 27.7% in the between-occasion column
    # of the CL row). Implemented as a NONMEM `$OMEGA BLOCK(1)` reused across
    # occasions via `BLOCK SAME`; nlmixr2 has no `SAME` shortcut, so the first
    # eta carries the estimated variance and the remainder are fixed to it.
    etaiov_cl_1 ~ 0.073928       # Li 2026 Table S2 Model c: BOV CL 27.7% (RSE 12), variance = log(0.277^2 + 1)
    etaiov_cl_2 ~ fix(0.073928)  # NONMEM $OMEGA BLOCK SAME (occasion 2 shares the occasion-1 variance)
    etaiov_cl_3 ~ fix(0.073928)  # NONMEM $OMEGA BLOCK SAME (occasion 3 shares the occasion-1 variance)

    # Residual error. Li 2026 Table S3 labels the NEP = 11 residual model
    # "prop", which is the RESERR option 0 of the deposited pyDarwin token file
    # (Li 2024 Supplementary Material 2): `IOBS = F*EXP(EPS(1))`, a log-normal
    # (exponential) residual that maps to `lnorm(expSd)`. The tabulated 0.029
    # is the NONMEM $SIGMA variance, so expSd = sqrt(0.029) = 0.170, which is
    # the same residual magnitude the stepwise three-compartment model of
    # Aregbe 2012 reports for this data set (16.1%).
    expSd <- 0.170294 ; label("Exponential (log-normal) residual error SD")  # Li 2026 Table S2 Model c: proportional error 0.029 (RSE 6); SD = sqrt(0.029)
  })

  model({
    # Decompose the integer-valued occasion column into binary indicators for
    # BOV multiplexing on log-CL (three occasions per the deposited control
    # stream). Unused indicator rows simply zero out.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)

    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3

    # Individual PK parameters; no covariate effects at NEP = 11. CL carries
    # both BSV and BOV; Q2, Q3 and V3 carry neither.
    cl  <- exp(lcl  + etalcl + iov_cl)
    vc  <- exp(lvc  + etalvc)
    q   <- exp(lq)
    vp  <- exp(lvp  + etalvp)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2)

    # Three-compartment IV PK with first-order elimination from the central
    # compartment (NONMEM ADVAN11 / TRANS4 mass-balance ODEs). Dose lands in
    # `central` via the data set's cmt column; the infusion duration is
    # encoded on the dose record (rate or dur).
    d/dt(central)     <-  q  / vp  * peripheral1 + q2 / vp2 * peripheral2 -
                          (cl + q + q2) / vc * central
    d/dt(peripheral1) <-  q  / vc  * central     - q  / vp  * peripheral1
    d/dt(peripheral2) <-  q2 / vc  * central     - q2 / vp2 * peripheral2

    # Plasma concentration. `central` is in mg and `vc` in L, so central / vc
    # is mg/L; the factor of 1000 converts to the ng/mL scale on which the
    # 17-DMAG assay data were modelled.
    Cc <- 1000 * central / vc

    # Log-normal residual error on the concentration scale.
    Cc ~ lnorm(expSd)
  })
}
