Li_2026_alvespimycin_nep16 <- function() {
  description <- "Three-compartment population PK model for the heat shock protein 90 inhibitor 17-DMAG (alvespimycin) given as an IV infusion to adult patients with advanced solid tumors, recovered as the 16-estimated-parameter non-dominated solution (OFV 8041.3) on the NSGA-II multi-objective Pareto front of Li 2026 and identical to the optimum found by the single-objective hybrid genetic algorithm, with first-order elimination, a power effect of body weight on central volume, log-normal IIV on CL, Vc, Vp and Q3, between-occasion variability on CL, Vc and Q2, and a combined residual error model. This is a machine-search model structure from a multi-objective model-selection study, not an expert-developed final model."
  reference <- paste(
    "Li X, Sale M, Craig J, Nieforth K, Mazur A, Bies RR. (2026).",
    "Multi-objective optimization in population pharmacokinetic model",
    "selection and optimization: application of NSGA-II in pyDarwin.",
    "J Pharmacokinet Pharmacodyn 53(1):26.",
    "doi:10.1007/s10928-026-10036-9.",
    "Parameter values from Supplementary Table S2 (Model d, NEP = 16).",
    "Structural parameterization (weight centering, occasion structure,",
    "residual-error form) and the between-occasion / between-subject",
    "assignment of the Q2 and Q3 variance terms are read from the deposited",
    "pyDarwin template, token and final control files for this same model,",
    "Supplementary Material 1, 2 and 4 of Li X, Sale M, Nieforth K, Craig J,",
    "Wang F, Solit D, Feng K, Hu M, Bies RR, Zhao L. (2024). pyDarwin machine",
    "learning algorithms application and comparison in nonlinear mixed-effect",
    "model selection and optimization.",
    "J Pharmacokinet Pharmacodyn 51(6):785-796.",
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
    WT = list(
      description        = "Total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters as a power function of weight centred on 81 kg acting on the central volume of distribution. The centring constant is not printed in Li 2026; it is read from the deposited pyDarwin template and final control stream for this same DMAG search space (Li 2024 Supplementary Material 1 and 4, `CWTKGONE = WT/81 ;; WEIGHT CENTERED ON ONE`, applied as `TVV = THETA(1)*CWTKGONE**THETA(7)`). Note that 81 kg is close to, but not identical to, the cohort median of 80.3 kg reported in Li 2026 Supplementary Table S1.",
      source_name        = "WT"
    ),
    OCC = list(
      description        = "Integer-valued occasion indicator for between-occasion variability multiplexing.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Values 1, 2, 3 identify the dosing occasion within subject. The deposited pyDarwin control stream for this model (Li 2024 Supplementary Material 4, `$PK`) writes the three BOV blocks as `IF(OCC.EQ.1) ... IF(OCC.EQ.2) ... IF(OCC.EQ.3)`. Decomposed inside `model()` into binary indicators `oc1`, `oc2`, `oc3` that multiplex the BOV etas on log-CL, log-Vc and log-Q2.",
      source_name        = "OCC"
    )
  )

  # Screened but not retained at NEP = 16 (Li 2026 Table 3 "Searched
  # covariates" for DMAG; Li 2026 Table S3 shows `V~WT` as the only retained
  # relationship at NEP = 16, which is confirmed by the deposited control
  # stream's phenotype line).
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units       = "years",
      type        = "continuous",
      notes       = "In the search space as a power function of age centred on 60 years on Vc and CL; not retained at NEP = 16."
    ),
    SEXF = list(
      description = "Sex indicator.",
      units       = "(binary)",
      type        = "categorical",
      notes       = "In the search space as an exponential effect of the NONMEM SEX column on Vc and CL; not retained at NEP = 16. The Li 2024 control stream codes the column as SEX with a 0/1 encoding whose reference level is not stated in either publication."
    ),
    CREAT = list(
      description = "Serum creatinine.",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "In the search space as a power function of serum creatinine centred on 0.9 mg/dL on CL (Li 2024 control stream `CSCR = SCR/0.9`); not retained at NEP = 16, though it does appear in the more complex non-dominated solutions further along the Pareto front."
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
    # Structural PK parameters from Li 2026 Supplementary Table S2, Model d
    # (NEP = 16, OFV = 8041.3). NONMEM ADVAN11 / TRANS4 (Li 2026 Table S3,
    # DMAG rows, "COM" = 3 at NEP = 16). NONMEM names map to nlmixr2
    # conventions as CL -> lcl, V1 -> lvc, Q2 -> lq, V2 -> lvp, Q3 -> lq2,
    # V3 -> lvp2. `lvc` is the typical central volume at the 81 kg centring
    # weight, not at the cohort median.
    lcl  <- log(8.64) ; label("Clearance (L/h)")                                          # Li 2026 Table S2 Model d: CL = 8.64 L/hr (RSE 6%)
    lvc  <- log(29.1) ; label("Central volume of distribution at WT = 81 kg (L)")         # Li 2026 Table S2 Model d: V  = 29.1 L    (RSE 7%)
    lq   <- log(79.3) ; label("Inter-compartmental clearance Q2 (L/h)")                   # Li 2026 Table S2 Model d: Q2 = 79.3 L/hr (RSE 4%)
    lvp  <- log(71.9) ; label("Peripheral volume of distribution V2 (L)")                 # Li 2026 Table S2 Model d: V2 = 71.9 L    (RSE 7%)
    lq2  <- log(9.74) ; label("Inter-compartmental clearance Q3 (L/h)")                   # Li 2026 Table S2 Model d: Q3 = 9.74 L/hr (RSE 7%)
    lvp2 <- log(136)  ; label("Peripheral volume of distribution V3 (L)")                 # Li 2026 Table S2 Model d: V3 = 136 L     (RSE 7%)

    # Covariate effect: power function of weight centred on 81 kg acting on the
    # central volume, i.e. `TVV = THETA(1) * (WT/81)^e_wt_vc` in the deposited
    # control stream (Li 2024 Supplementary Material 4, `$PK`).
    e_wt_vc <- 1.31 ; label("Power exponent of (WT/81 kg) on central volume (unitless)")  # Li 2026 Table S2 Model d: power relationship between weight and volume = 1.31 (RSE 22%)

    # Between-subject variability. Li 2026 Table S2 reports BSV as a CV%; the
    # internal log-normal variance is omega^2 = log(CV^2 + 1). BSV was retained
    # on V, CL, V2 and Q3 (Li 2026 Table S3, "BSV" column reads "V, CL, V2,
    # Q3"; the deposited control stream writes `Q3 = THETA(4)*EXP(ETA(12))`).
    etalcl ~ 0.225548  # Li 2026 Table S2 Model d: BSV CL 50.3% (RSE 11), variance = log(0.503^2 + 1)
    etalvc ~ 0.168933  # Li 2026 Table S2 Model d: BSV V  42.9% (RSE 24), variance = log(0.429^2 + 1)
    etalvp ~ 0.296074  # Li 2026 Table S2 Model d: BSV V2 58.7% (RSE 11), variance = log(0.587^2 + 1)
    etalq2 ~ 0.453978  # Li 2026 Table S2 Model d: BSV Q3 75.8% (RSE 10), variance = log(0.758^2 + 1)

    # Between-occasion variability on CL, Vc and Q2 (Li 2026 Table S3, "BOV"
    # column reads "V, CL, Q2" at NEP = 16; the deposited control stream
    # defines `IOVCL`, `IOVV` and `IOVQ2`). NOTE: Li 2026 Supplementary Table
    # S2 prints the 28.6% BOV term on the Q3 row rather than the Q2 row; the
    # Li 2024 table of the identical model, Li 2026's own Table S3, and the
    # deposited control stream all place it on Q2. See the vignette Errata.
    # Implemented as NONMEM `$OMEGA BLOCK(1)` reused across occasions via
    # `BLOCK SAME`; nlmixr2 has no `SAME` shortcut, so the first eta carries
    # the estimated variance and the remainder are fixed equal to it.
    etaiov_cl_1 ~ 0.087836       # Li 2026 Table S2 Model d: BOV CL 30.3% (RSE 13), variance = log(0.303^2 + 1)
    etaiov_cl_2 ~ fix(0.087836)  # NONMEM $OMEGA BLOCK SAME (occasion 2 shares the occasion-1 variance)
    etaiov_cl_3 ~ fix(0.087836)  # NONMEM $OMEGA BLOCK SAME (occasion 3 shares the occasion-1 variance)

    etaiov_vc_1 ~ 0.092325       # Li 2026 Table S2 Model d: BOV V 31.1% (RSE 34), variance = log(0.311^2 + 1)
    etaiov_vc_2 ~ fix(0.092325)  # NONMEM $OMEGA BLOCK SAME (occasion 2 shares the occasion-1 variance)
    etaiov_vc_3 ~ fix(0.092325)  # NONMEM $OMEGA BLOCK SAME (occasion 3 shares the occasion-1 variance)

    etaiov_q_1 ~ 0.078623        # Li 2026 Table S2 Model d: BOV 28.6% (RSE 13), on Q2 per Li 2024 and the control stream, variance = log(0.286^2 + 1)
    etaiov_q_2 ~ fix(0.078623)   # NONMEM $OMEGA BLOCK SAME (occasion 2 shares the occasion-1 variance)
    etaiov_q_3 ~ fix(0.078623)   # NONMEM $OMEGA BLOCK SAME (occasion 3 shares the occasion-1 variance)

    # Residual error. Li 2026 Table S3 labels the NEP = 16 residual model
    # "comb", which is RESERR option 1 of the deposited pyDarwin token file
    # (Li 2024 Supplementary Material 2): `IOBS = F*EXP(EPS(1)) + EPS(2)`. The
    # multiplicative half is therefore exponential rather than strictly
    # proportional, but rxode2 cannot solve an `lnorm() + add()` endpoint, so
    # it is encoded as `prop() + add()`; the two agree to second order at this
    # residual magnitude (exp(0.13) - 1 = 0.139). Both tabulated values are
    # NONMEM $SIGMA variances, so the SDs are sqrt(0.017) and sqrt(14.7), the
    # latter on the ng/mL concentration scale. See the vignette Errata.
    propSd <- 0.130384 ; label("Proportional residual error (fraction)")     # Li 2026 Table S2 Model d: proportional error 0.017 (RSE 7); SD = sqrt(0.017) = 0.130, matching the 13.2% that Li 2024 reports for this same model
    addSd  <- 3.834058 ; label("Additive residual error SD (ng/mL)")         # Li 2026 Table S2 Model d: additive error 14.7 (RSE 36); SD = sqrt(14.7)
  })

  model({
    # Decompose the integer-valued occasion column into binary indicators for
    # BOV multiplexing on log-CL, log-Vc and log-Q2 (three occasions per the
    # deposited control stream). Unused indicator rows simply zero out.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)

    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2 + oc3 * etaiov_vc_3
    iov_q  <- oc1 * etaiov_q_1  + oc2 * etaiov_q_2  + oc3 * etaiov_q_3

    # Weight is centred on 81 kg, so a subject at the centring weight has the
    # typical volume `exp(lvc)` (Li 2024 control stream: `CWTKGONE = WT/81`,
    # `TVV = THETA(1)*CWTKGONE**THETA(7)`).
    ref_wt <- 81

    # Individual PK parameters. Vc carries the weight effect, BSV and BOV; CL
    # carries BSV and BOV; Q2 carries BOV only; Q3 and V2 carry BSV only; V3
    # carries neither.
    cl  <- exp(lcl  + etalcl + iov_cl)
    vc  <- exp(lvc  + etalvc + iov_vc) * (WT / ref_wt)^e_wt_vc
    q   <- exp(lq            + iov_q)
    vp  <- exp(lvp  + etalvp)
    q2  <- exp(lq2  + etalq2)
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
    # 17-DMAG assay data were modelled and on which `addSd` is expressed.
    Cc <- 1000 * central / vc

    # Combined residual error on the concentration scale.
    Cc ~ prop(propSd) + add(addSd)
  })
}
