Bizzotto_2016_glucose <- function() {
  description <- "Mechanistic model of glucose tracer kinetics in humans driven by time-varying plasma insulin and glucose regressors (Bizzotto 2016). Glucose uptake is a Michaelis-Menten function of glucose at the site of action whose maximum rate Vmax is itself a Hill (sigmoidal) function of insulin at the site of action; this captures the observation that hyperglycemia suppresses the glucose-clearance response to hyperinsulinemia. Two-compartment delays smooth plasma insulin and glucose into their site-of-action analogues, and a heart-lung block plus a three-channel periphery block give the tracer disposition. Distributed in the DDMORE Foundation Model Repository (DDMODEL00000227) as a simulation-only implementation; the linked publication fits the same equations to real data from 123 subjects spanning normal-tolerant, impaired-glucose-tolerance, and type 2 diabetic adults."
  reference <- paste(
    "Bizzotto R, Natali A, Gastaldelli A, Muscelli E, Brehm A, Roden M, Ferrannini E, Mari A. (2016).",
    "Glucose uptake saturation explains glucose kinetics profiles measured by different tests.",
    "Am J Physiol Endocrinol Metab 311(2):E346-E357.",
    "doi:10.1152/ajpendo.00045.2016.",
    "DDMORE Foundation Model Repository: DDMODEL00000227.",
    sep = " "
  )
  vignette <- "Bizzotto_2016_glucose"
  units <- list(time = "min", dosing = "umol/m^2", concentration = "mmol/L")
  ddmore_id    <- "DDMODEL00000227"
  replicate_of <- NULL

  covariateData <- list(
    INS = list(
      description        = "Plasma insulin concentration time-course (regressor input)",
      units              = "pmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying regressor; supplied at every observation/event row in the dataset. Linearly interpolated between rows via the `linear(INS)` declaration in `model()`. Drives the two-compartment delay (Z1 -> Z) representing insulin at the site of action. The DDMORE bundle's Simulated_glucoseKinetics.csv carries this column as `iins` (insulin at the current row time); rename to `INS` before passing the dataset to rxSolve.",
      source_name        = "iins"
    ),
    GLU = list(
      description        = "Plasma glucose concentration time-course (regressor input)",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying regressor; supplied at every observation/event row in the dataset. Linearly interpolated between rows via the `linear(GLU)` declaration in `model()`. Drives the two-compartment delay (X1 -> X) representing glucose at the site of action. The DDMORE bundle's Simulated_glucoseKinetics.csv carries this column as `iglu` (glucose at the current row time); rename to `GLU` before passing the dataset to rxSolve.",
      source_name        = "iglu"
    )
  )

  population <- list(
    n_subjects     = 123L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "Adults spanning the glucose-tolerance spectrum (normal-tolerant, impaired-glucose-tolerance, and type 2 diabetic). 123 subjects across five experimental tests: three-step hyperglycemic-hyperinsulinemic clamp (HGclamp, n=8), two-step isoglycemic-hyperinsulinemic clamp (ISOclamp, n=8), paired oral glucose tolerance test plus euglycemic clamp (OGTT/clamp, n=8), mixed-meal test (MTT, n=91), and paired mixed-meal test plus hyperglycemic clamp (MTT/clamp, n=8). Specific demographic detail (age range, weight range, sex distribution) is described in the Bizzotto 2016 publication but not reproduced in the DDMORE bundle and the publication PDF is not on disk in this worktree.",
    dose_range     = "Glucose tracer (e.g. [6,6-d2]glucose) administered as an intravenous bolus and continuous infusion; tracer doses normalised to body surface area (umol/m^2 for bolus, umol/min/m^2 for infusion). The DDMORE bundle's Simulated_glucoseKinetics.csv ships test-typical bolus + infusion regimens (e.g. ~1000 umol/m^2 bolus + ~8 umol/min/m^2 infusion in the clamp arms).",
    regions        = NA_character_,
    notes          = "Population n and test-by-test n derived from the Long_technical_model_description_glucoseKinetics.txt file in the DDMORE bundle for DDMODEL00000227; the bundle states the test mix and total subject count (123) match the related publication. Demographic detail (age, weight, sex, region) is in the Bizzotto 2016 publication but the publication PDF is not on disk in this worktree, so finer-grained population descriptors are recorded as NA. See the validation vignette's Errata section for the full list of bundle-versus-publication caveats."
  )

  ini({
    # Structural parameters - DDMODEL00000227 glucoseKinetics.mdl `gu_v1_par` `STRUCTURAL`
    # block (lines 47-61), labelled "final parameter estimates from related publication"
    # (Bizzotto 2016). The bundle is shipped as a simulation-only model: the
    # Output_real_glucoseKinetics.txt listing is from re-fitting the simulated
    # dataset (which the bundle itself flags as meaningless for estimation), so
    # the .mdl `parObj` STRUCTURAL block - not the Output_real listing - holds
    # the publication-derived point estimates that drive the simulation.
    lkmg   <- log(3.88)  ; label("Glucose Michaelis-Menten Km for peripheral uptake (mmol/L)")               # typ_KmG
    lvmax0 <- log(338)   ; label("Basal Vmax of glucose uptake (umol/min/m^2)")                              # typ_Vmax0
    lemax  <- log(4812)  ; label("Insulin-driven asymptotic increase in Vmax (umol/min/m^2)")                # typ_Emax
    lhill <- log(1.62)  ; label("Hill exponent on insulin-driven Vmax saturation (unitless)")              # typ_gamma
    lkmi   <- log(784)   ; label("Insulin half-saturation for Vmax modulation (pmol/L)")                     # typ_KmI
    lt12i  <- log(15.9)  ; label("Insulin-at-action two-compartment delay half-life (min)")                  # typ_t12I
    lt12g  <- fixed(log(0.7))  ; label("Glucose-at-action two-compartment delay half-life (min) - FIXED")    # typ_t12G fix=true
    lvtot  <- log(12648) ; label("Total apparent volume of glucose distribution (mL)")                       # typ_V
    lflambda3 <- log(0.0582 / (1 - 0.0582)) ; label("Logit of flambda3 (ratio of channel-3 to channel-2 periphery rate constants)")  # typ_flambda3
    lflambda2 <- log(0.154  / (1 - 0.154))  ; label("Logit of flambda2 (ratio of channel-2 to channel-1 periphery rate constants)")  # typ_flambda2
    lw1       <- log(0.609  / (1 - 0.609))  ; label("Logit of w1 (fraction of periphery flow into channel 1)")                       # typ_w1
    lfw2      <- log(0.901  / (1 - 0.901))  ; label("Logit of fw2 (fraction of remaining periphery flow into channel 2 after w1)")   # typ_fw2
    lpflow <- fixed(log(2688)) ; label("Peripheral flow rate F (mL/min/m^2) - FIXED at 3200*0.84")           # typ_F fix=true (3200*0.84)

    # Inter-individual variability - DDMODEL00000227 glucoseKinetics.mdl `gu_v1_par`
    # `VARIABILITY` block (lines 63-80). All etas are Normal(0, var) on the log-scale
    # (`type is var`). Components with `var = 0, fix=true` in the .mdl have no eta in
    # nlmixr2 (no IIV is the same as a fixed-zero variance):
    #   var_Vmax0, var_t12G, var_flambda2, var_fw2, var_F all = 0 FIXED.
    # The bundle declares one correlation: corr_gamma_KmI = -0.44 between eta_gamma
    # and eta_KmI; the corresponding covariance for the nlmixr2 block matrix is
    #   cov = corr * sqrt(var_gamma * var_KmI) = -0.44 * sqrt(0.111 * 0.263) = -0.0752.
    etalkmg                    ~ 0.219                              # var_KmG
    etalemax                   ~ 0.112                              # var_Emax
    etalhill + etalkmi        ~ c(0.111, -0.0752, 0.263)           # var_gamma, cov(hill,KmI), var_KmI
    etalt12i                   ~ 0.151                              # var_t12I
    etalvtot                   ~ 0.0557                             # var_V
    etalflambda3               ~ 0.179                              # var_flambda3
    etalw1                     ~ 0.773                              # var_w1

    # Residual error - DDMODEL00000227 glucoseKinetics.mdl `gu_v1_par` `VARIABILITY`
    # alpha = 0.014 mmol/L, additive on the tracer-concentration observation Y=G with
    # epsilon ~ Normal(0, sigma) and sigma = 1 fixed (Monolix-standard form). That is a
    # plain additive residual SD of 0.014 mmol/L on the tracer concentration G in mmol/L.
    addSd <- 0.014 ; label("Additive residual error on tracer concentration (mmol/L)")  # alpha
  })
  model({
    # Mark the regressor inputs as linearly interpolated between dataset rows; the model
    # mirrors the .mdl regressors INS (plasma insulin, pmol/L) and GLU (plasma glucose,
    # mmol/L). The .mdl uses a hand-rolled piecewise-linear interpolation via T1 / TOBS /
    # GLU1 / INS1 columns to enforce linear behaviour irrespective of the simulator's
    # default; rxode2's `linear(...)` declaration achieves the same thing natively, so
    # only the two driving regressors INS and GLU are required in the input dataset.
    linear(INS, GLU)

    # Mechanistic constants - DDMODEL00000227 glucoseKinetics.mdl `MODEL_PREDICTION`
    # block (lines 172-174). VHL is the heart-lung volume per body surface area, deltaHL
    # the heart-lung block elimination rate, and delta the periphery-block sink rate.
    VHL     <- 700  # mL/m^2
    deltaHL <- 15   # 1/min
    delta   <- 10   # 1/min

    # Individual parameters - log-normal etas around the typical values, except for the
    # logit-normal fractions flambda2/flambda3/w1/fw2 which use logit(typ) + eta on the
    # unconstrained scale and back-transform with the inverse-logit.
    kmg      <- exp(lkmg   + etalkmg)
    vmax0    <- exp(lvmax0)                                    # var_Vmax0  fixed at 0  -> no eta
    emax     <- exp(lemax  + etalemax)
    hill    <- exp(lhill + etalhill)
    kmi      <- exp(lkmi   + etalkmi)
    t12i     <- exp(lt12i  + etalt12i)
    t12g     <- exp(lt12g)                                     # var_t12G   fixed at 0  -> no eta
    vtot     <- exp(lvtot  + etalvtot)
    # Logit-scale individual raw values (mu-referenced) before the inverse-logit
    # transform. Splitting these onto their own lines keeps nlmixr2's mu-reference
    # detection happy for the etas on flambda3 and w1.
    lflambda3i <- lflambda3 + etalflambda3
    lw1i       <- lw1       + etalw1
    flambda3 <- 1 / (1 + exp(-lflambda3i))
    flambda2 <- 1 / (1 + exp(-lflambda2))                      # var_flambda2 fixed at 0 -> no eta
    w1       <- 1 / (1 + exp(-lw1i))
    fw2      <- 1 / (1 + exp(-lfw2))                           # var_fw2    fixed at 0  -> no eta
    pflow    <- exp(lpflow)                                    # var_F      fixed at 0  -> no eta

    # Derived periphery-block rate constants - .mdl MODEL_PREDICTION lines 177-183.
    # w2 = (1 - w1) * fw2 and w3 = 1 - w1 - w2 give a simplex (w1 + w2 + w3 = 1) over the
    # three periphery channels. lambda1/2/3 are the channel-specific rate constants
    # derived from the "fast / medium / slow" ratios flambda2 = lambda2/lambda1 and
    # flambda3 = lambda3/lambda2; the explicit form keeps lambda3 < lambda2 < lambda1 by
    # construction (flambda{2,3} < 1). c1/c2 are heart-lung output coefficients.
    w2       <- (1 - w1) * fw2
    w3       <- 1 - w1 - w2
    lambda1  <- (w1 * flambda2 * flambda3 + w2 * flambda3 + w3) / (flambda2 * flambda3) *
                delta * pflow / (delta * (vtot - VHL) - pflow)
    lambda2  <- lambda1 * flambda2
    lambda3  <- lambda2 * flambda3
    c1       <- deltaHL * pflow / (deltaHL * VHL - 2 * pflow)
    c2       <- -deltaHL * pflow / (deltaHL * VHL - pflow)

    # Site-of-action insulin (Z) and glucose (X), assembled from the dynamic states.
    # Vmax_eff is the Hill-modulated maximum glucose uptake; cl is the saturable
    # peripheral clearance; E is the per-pass fractional extraction (clearance over
    # peripheral flow). Defined here so the periphery-block ODEs below see the
    # current-time values.
    Vmax_eff <- vmax0 + emax * Z^hill / (kmi^hill + Z^hill)
    cl       <- Vmax_eff / (kmg + X)
    E        <- cl / pflow

    # Tracer concentration in the heart-lung output (G, mmol/L) and the total flux back
    # into the heart-lung block from the periphery (Gv); the heart-lung ODEs and the
    # periphery ODEs reference both.
    G  <- c1 * xHL1 - c1 * xHL2
    Gv <- delta * xPER4

    # Bioavailability scaling on the tracer-dose targets - .mdl `COMPARTMENT` block
    # `phi1: {finput=1/F}` and `phi2: {finput=1/F}`. A bolus or infusion of tracer
    # delivered to xHL1 / xHL2 enters the compartment scaled by 1 / pflow, which converts
    # the tracer-mass dose (umol/m^2 or umol/min/m^2) into the (mmol/L)*min units carried
    # by the heart-lung states.
    f(xHL1) <- 1 / pflow
    f(xHL2) <- 1 / pflow

    # Two-compartment delays smoothing GL and INS into the site-of-action quantities
    # X (glucose) and Z (insulin). .mdl `DEQ` block, lines 204-207.
    d/dt(X1)  <- (GLU - X1) * log(2) / t12g
    d/dt(X)   <- (X1  - X)  * log(2) / t12g
    d/dt(Z1)  <- (INS - Z1) * log(2) / t12i
    d/dt(Z)   <- (Z1  - Z)  * log(2) / t12i

    # Heart-lung block - .mdl DEQ lines 218-219. Tracer enters here via dosing on
    # xHL1 (bolus or infusion) and / or xHL2 (typically an offsetting infusion in the
    # bundle's regimens), and recirculates from the periphery via Gv.
    d/dt(xHL1) <- c2       * xHL1 + Gv
    d/dt(xHL2) <- -deltaHL * xHL2 + Gv

    # Periphery block - three parallel channels with weights w1 / w2 / w3 and channel-
    # specific drains lambda1 / lambda2 / lambda3, summed into a single return reservoir
    # xPER4 that drains at rate delta back into the heart-lung block as Gv. .mdl DEQ
    # lines 223-228.
    d/dt(xPER1) <- -lambda1 * xPER1 + w1 * (1 - E) * G
    d/dt(xPER2) <- -lambda2 * xPER2 + w2 * (1 - E) * G
    d/dt(xPER3) <- -lambda3 * xPER3 + w3 * (1 - E) * G
    d/dt(xPER4) <- lambda1  * xPER1 + lambda2 * xPER2 + lambda3 * xPER3 - delta * xPER4

    # Observation: tracer concentration G in mmol/L with additive residual error.
    # Y = G + alpha * epsilon, epsilon ~ N(0, sigma=1) in the .mdl OBSERVATION block
    # (line 239), which is a plain `add(addSd)` model in nlmixr2 with addSd = alpha.
    G ~ add(addSd)
  })
}
