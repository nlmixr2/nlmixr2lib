Vorontsova_2022_misoprostol <- function() {
  description <- "One-compartment population pharmacokinetic model of misoprostol acid (MPA, the active metabolite) following buccal or vaginal misoprostol tablet administration for full-term labor induction, with mixed linear plus Michaelis-Menten clearance. The absorption rate constant is estimated separately by route (buccal vs vaginal) and dose level (25 vs 50 microgram); a relative bioavailability multiplier captures the higher vaginal exposure, with inter-occasion variability across the first two dose events. IIV is estimated on apparent clearance, apparent central volume of distribution, and the absorption rate constant. Development cohort: 47 women at term gestation (>=37 weeks) undergoing labor induction in the IMPROVE trial (NCT02408315)."
  reference <- paste(
    "Vorontsova Y, Haas DM, Flannery K, Masters AR, Silva LL, Pierson RC, Yeley B,",
    "Hogg G, Guise D, Heathman M, Quinney SK. (2022). Pharmacokinetics of vaginal",
    "versus buccal misoprostol for labor induction at term.",
    "Clin Transl Sci 15(8):1937-1945. doi:10.1111/cts.13306.",
    sep = " "
  )
  vignette <- "Vorontsova_2022_misoprostol"
  units <- list(time = "h", dosing = "ng", concentration = "pg/mL", volume = "L", clearance = "L/h")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot   = list(analyte = "misoprostol", units = "ng", specimen = "administration site", verified = FALSE),
    central = list(analyte = "misoprostol", units = "ng", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    ROUTE_VAGINAL = list(
      description        = "Administration route (vaginal vs buccal): 1 = vaginal misoprostol, 0 = buccal misoprostol.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (buccal)",
      notes              = "Per-subject: each participant was randomized to receive misoprostol via either the vaginal or buccal route for the entire study; placebo was administered via the opposite route so participants and providers were blinded. `ROUTE_VAGINAL = 1` selects the vaginal absorption rate constants and applies the relative bioavailability multiplier `F_v/b` (typical value 2.3) on the depot dose; `ROUTE_VAGINAL = 0` (buccal reference) applies F = 1 and the buccal ka values. New canonical registered alongside this extraction; see `inst/references/covariate-columns.md` H3 entry `ROUTE_VAGINAL`.",
      source_name        = "ROUTE"
    ),
    DOSE = list(
      description        = "Nominal per-dose amount of misoprostol tablet.",
      units              = "microgram",
      type               = "continuous",
      reference_category = NULL,
      notes              = "In the source trial the first (loading) dose was 25 microgram and subsequent doses (up to 6 additional doses) were 50 microgram every 4 h. Together with `ROUTE_VAGINAL` this covariate selects between the four estimated absorption rate constants (buccal 25 microgram, buccal 50 microgram, vaginal 25 microgram, vaginal 50 microgram). Only 25 and 50 microgram values were studied; simulations with other DOSE values collapse both `(DOSE == 25)` and `(DOSE == 50)` indicators to zero and yield ka = exp(0 + etalka), which does not reflect the paper -- restrict simulated dose values to {25, 50}. DOSE is a per-dose-record covariate that should be carried forward (last-observation-carried-forward) between successive dose events so ka correctly reflects the current dose.",
      source_name        = "DOSE"
    ),
    OCC = list(
      description        = "Integer-valued dosing-occasion indicator for inter-occasion-variability multiplexing on relative bioavailability.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Values 1 and 2 identify the first (25 microgram) and second (50 microgram) dose occasions respectively, per the source protocol which collected PK samples for the first two doses only. Decomposed inside `model()` into binary indicators `oc1` and `oc2` that multiplex the two IOV etas on log `F_v/b`. Occasion 2 IOV variance is fixed equal to the occasion-1 estimate (NONMEM `$OMEGA BLOCK(1) SAME` convention) because nlmixr2 has no direct SAME shortcut. For records outside the first two dose intervals pass OCC = 1 or OCC = 2 as appropriate; values outside {1, 2} collapse both indicators to zero and remove the IOV shift.",
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index (screened, not retained).",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Assessed as a covariate on CL/Fb, V/Fb, and ka but no relationship identified with empirical Bayes estimates (Results paragraph 'No correlations were identified among BMI or age and empirical Bayes estimates for V/Fb, CL/Fb, or ka (Figure S1).'). Not retained. Cohort BMI median 35 (buccal) / 36 (vaginal) kg/m^2 per Table 1."
    ),
    AGE = list(
      description = "Maternal age (screened, not retained).",
      units       = "years",
      type        = "continuous",
      notes       = "Assessed as a covariate but no relationship identified (Figure S1). Cohort age median 26 (buccal) / 24 (vaginal) years per Table 1."
    ),
    GA = list(
      description = "Gestational age (screened, not retained).",
      units       = "weeks",
      type        = "continuous",
      notes       = "Assessed as a covariate but not retained. Cohort gestational age median 39.0 (buccal) / 39.8 (vaginal) weeks per Table 1."
    ),
    RACE_BLACK = list(
      description = "African American / Black race indicator (screened, not retained).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Race was reported as African American/Black (30-33% across arms), White (46-52%), Other (17-21%) per Table 1. Race was assessed as a covariate but not retained ('body mass index, race, gestational age, and maternal age were assessed as covariates' -- none retained in the final model)."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 47L,
    n_studies        = 1L,
    n_observations   = 469L,
    age_range        = "median 24 (95% CI 21.6-26.4) years vaginal arm; median 26 (95% CI 23.2-28.7) years buccal arm; eligibility >=14 years",
    weight_range     = "not reported; BMI median 35-36 kg/m^2 (95% CI 31.6-39.8 across arms)",
    gestational_age_range = "median 39.0 (95% CI 38.5-39.6) weeks buccal arm; 39.8 (95% CI 39.3-40.3) weeks vaginal arm; eligibility >=37 0/7 weeks",
    sex_female_pct   = 100,
    race_ethnicity   = c(African_American_or_Black = 32, White = 49, Other = 19),
    disease_state    = "Term pregnancy (>=37 weeks gestation) undergoing labor induction. Women were eligible if >=14 years old with a singleton cephalic-presentation pregnancy, modified Bishop Score <=6, and either a medically indicated induction beyond 37 0/7 weeks or an elective induction after 39 0/7 weeks. Excluded: known misoprostol allergy, prior uterine scarring, untreated cervical infection, prior induction/cervical ripening this pregnancy, contraindications to labor induction or misoprostol, planned cesarean for maternal/fetal condition, known major fetal congenital anomaly, or category 2/3 fetal tracing before initiating induction.",
    dose_range       = "Misoprostol 25 microgram initial dose then 50 microgram every 4 h if clinically indicated, up to 24 h (maximum 7 doses total). Tablets split from a 100 microgram Novel Laboratories tablet by investigational pharmacy. PK sampling: pre-dose and ~0.25, 0.5, 1, 2, 3, and 4 h after each of the first two doses.",
    regions          = "USA (Eskenazi Health and Indiana University Health Methodist Hospitals, Indianapolis, IN).",
    study_ids        = "NCT02408315 (IMPROVE trial); FDA IND 122727.",
    notes            = "PK substudy nested within the IMPROVE randomized triple-masked placebo-controlled noninferiority trial of buccal vs vaginal misoprostol for labor induction at term. Fifty women initially enrolled; three excluded (two vaginal-arm subjects expulsed the tablet shortly after placement; one subject was withdrawn shortly after the first dose at provider request); 47 analyzed (24 buccal, 23 vaginal). MPA quantified by LC-MS/MS at the IU Simon Cancer Center CPAC Laboratory (LOQ = 0.3 pg/mL; three BLQ observations after the initial dose excluded); 469 plasma MPA observations retained. NONMEM 7.5 FOCEI, PsN, Pirana 3.0.0; R 3.5.3 for data management and diagnostics. Bootstrap resampling n = 200 for parameter precision; VPC n = 200 stratified by route."
  )

  ini({
    # Structural PK parameters -- Vorontsova 2022 Table 2 "Population parameter
    # estimates" NONMEM final Estimate column. Table caption unit conventions
    # (with in-file interpretations noted for the two typographical anomalies):
    #   CL/Fb   in L/h
    #   V/Fb    in L
    #   F_v/b   unitless
    #   ka      in 1/h (four separate rows for buccal/vaginal x 25/50 microgram)
    #   Vmax/Fb in pg/mL as printed; interpreted as pg/(mL*h) (concentration-rate
    #           form of the Michaelis-Menten equation, dCc/dt|_MM = -Vmax*Cc/(Km+Cc)).
    #           The "/h" is missing from the table caption but is required for
    #           dimensional consistency with dCc/dt. Documented in vignette Errata.
    #   Km      in pg as printed; interpreted as pg/mL (concentration form
    #           matching Cc). The "/mL" is missing from the table caption but
    #           is required for dimensional consistency. Documented in vignette Errata.
    #
    # Buccal is the structural reference for F (F_b = 1 by construction; paper
    # Results 'Observed MPA concentration versus time plots demonstrated that
    # absorption and bioavailability differed by route of administration
    # (Figure 1). Thus, we included a relative bioavailability term (F_v/b) for
    # the vaginal route relative to buccal (F_b, assumed to be one)').
    lcl              <- log(730)      ; label("Apparent linear clearance CL/Fb (L/h)")                    # Table 2 row 'CL/Fb, L/h' Estimate = 730 (RSE 22.5%)
    lvc              <- log(610)      ; label("Apparent central volume of distribution V/Fb (L)")         # Table 2 row 'V/Fb, L' Estimate = 610 (RSE 33.4%)
    lvmax            <- log(5.45)     ; label("Apparent Michaelis-Menten Vmax/Fb (pg/mL/h)")              # Table 2 row 'Vmax/Fb, pg/ml' Estimate = 5.45 (RSE 12.8%); units interpretation pg/(mL*h) documented in vignette Errata
    lkm              <- log(2.5)      ; label("Michaelis-Menten constant Km (pg/mL)")                     # Table 2 row 'Km, pg' Estimate = 2.5 (RSE 41.2%); units interpretation pg/mL documented in vignette Errata
    lfvb             <- log(2.3)      ; label("Relative bioavailability vaginal vs buccal, F_v/b (unitless)")        # Table 2 row 'F_v/b' Estimate = 2.3 (RSE 20.8%)
    lka_buccal_25    <- log(0.709)    ; label("ka buccal 25 microgram dose (1/h)")                        # Table 2 row 'ka, 1/h (buccal, 25 microgram)' Estimate = 0.709 (RSE 15.7%)
    lka_buccal_50    <- log(0.537)    ; label("ka buccal 50 microgram dose (1/h)")                        # Table 2 row 'ka, 1/h (buccal, 50 microgram)' Estimate = 0.537 (RSE 16.2%)
    lka_vaginal_25   <- log(0.464)    ; label("ka vaginal 25 microgram dose (1/h)")                       # Table 2 row 'ka, 1/h (vaginal, 25 microgram)' Estimate = 0.464 (RSE 36%)
    lka_vaginal_50   <- log(0.240)    ; label("ka vaginal 50 microgram dose (1/h)")                       # Table 2 row 'ka, 1/h (vaginal, 50 microgram)' Estimate = 0.240 (RSE 28.9%)

    # Inter-individual variability (IIV). Paper Methods 'Pharmacokinetic
    # analysis' paragraph: 'Interindividual and interoccasion variabilities
    # were assumed to be log-normally distributed.' IIV was included on CL, V,
    # and ka only (Results 'Interindividual variability (IIV) was estimated for
    # apparent clearance (CL/Fb), apparent volume of distribution (V/Fb), and
    # ka'). Table 2 Omega block Estimate column values are interpreted as
    # omega^2 (log-scale variances).
    etalcl  ~ 0.602      # Table 2 Omega row 'CL, L/h' Estimate = 0.602 (RSE 33.7%, shrinkage 15.6%)
    etalvc  ~ 1.55       # Table 2 Omega row 'V, L' Estimate = 1.55 (RSE 23.5%, shrinkage 16%)
    etalka  ~ 0.0599     # Table 2 Omega row 'ka, 1/h' Estimate = 0.0599 (RSE 106.2%, shrinkage 38.3%)

    # Inter-occasion variability (IOV) on log F_v/b across the two dose
    # occasions (paper Results 'incorporation of IOV on F_v/b significantly
    # decreased OFV'). nlmixr2 has no NONMEM `$OMEGA BLOCK(1) SAME` shortcut,
    # so occasion 2 gets its own eta with variance fixed equal to the
    # occasion-1 estimate (matching the SAME convention used in
    # Jonsson_2011_ethambutol.R and Aregbe_2012_alvespimycin.R).
    etaiov_fvb_1 ~ 0.0848         # Table 2 Omega row 'IOV' Estimate = 0.0848 (RSE 60.3%, shrinkage 39.4%) -- estimated for occasion 1
    etaiov_fvb_2 ~ fix(0.0848)    # SAME-equivalent: equal to occasion-1 IOV variance

    # Residual variability -- proportional error only (paper Methods:
    # 'Residual variability was modeled by a proportional error model.').
    # Table 2 Sigma row 'PROP' Estimate = 0.262 is the NONMEM sigma^2 on the
    # linear (untransformed) scale; nlmixr2 uses SD form so propSd =
    # sqrt(sigma^2).
    propSd <- sqrt(0.262) ; label("Proportional residual error (fraction)")   # Table 2 Sigma row 'PROP' Estimate = 0.262 (RSE 10.3%, shrinkage 9%); sigma^2 -> SD via sqrt
  })

  model({
    # Decompose the integer-valued occasion covariate into binary indicators
    # for IOV multiplexing on log F_v/b (paper: two dose occasions sampled).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fvb <- oc1 * etaiov_fvb_1 + oc2 * etaiov_fvb_2

    # Dose- and route-dependent absorption rate constant. Four separate ka
    # values indexed by (ROUTE_VAGINAL, DOSE) share a single IIV eta on the
    # log scale (paper Results 'Addition of separate ka for each dose (25 and
    # 50 microgram) and route of administration ... significantly decreased
    # OFV').
    lka <- (1 - ROUTE_VAGINAL) * ((DOSE == 25) * lka_buccal_25 + (DOSE == 50) * lka_buccal_50) +
                ROUTE_VAGINAL  * ((DOSE == 25) * lka_vaginal_25 + (DOSE == 50) * lka_vaginal_50)
    ka <- exp(lka + etalka)

    # Individual PK parameters
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    vmax <- exp(lvmax)
    km   <- exp(lkm)
    fvb  <- exp(lfvb + iov_fvb)

    # Route-based bioavailability applied to the depot dose. Buccal
    # (ROUTE_VAGINAL = 0) uses F = 1; vaginal (ROUTE_VAGINAL = 1) applies the
    # relative bioavailability multiplier F_v/b (with occasion-specific IOV).
    fdepot <- 1 - ROUTE_VAGINAL + ROUTE_VAGINAL * fvb

    # One-compartment PK with first-order absorption. Elimination combines a
    # linear clearance term (CL / Vc) and a saturable Michaelis-Menten term
    # (dCc/dt|_MM = -Vmax * Cc / (Km + Cc); scaled by Vc to give amount
    # rate). Amounts in ng and Vc in L imply Cc = ng/L = pg/mL (numerically
    # the same units the paper reports).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot -
                      (cl / vc) * central -
                      vc * vmax * (central / vc) / (km + central / vc)

    # Route-based bioavailability on the depot dose
    f(depot) <- fdepot

    # Concentration in plasma
    Cc <- central / vc

    # Proportional residual error on the linear scale
    Cc ~ prop(propSd)
  })
}
