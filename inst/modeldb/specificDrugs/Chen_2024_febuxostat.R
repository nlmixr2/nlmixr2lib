Chen_2024_febuxostat <- function() {
  description <- "Two-compartment population PK model for febuxostat in 128 healthy Chinese adult volunteers (Chen 2024), with a dual-depot parallel first-order absorption model built to reproduce the multiple-peak plasma profiles of this BCS class II compound. A fraction F1 of the dose is absorbed from a first depot (rate ka1, lag time Tlag1) and the remaining 1 - F1 from a second depot (rate ka2, lag time Tlag2). A standardised high-fat high-calorie meal slows both absorption rate constants and lengthens Tlag1; body weight scales the apparent central and peripheral volumes."
  reference <- paste(
    "Chen W, Jiang B, Ruan Z, Yang D, Hu Y, Lou H.",
    "Population pharmacokinetic analysis of febuxostat with high focus on",
    "absorption kinetics and food effect.",
    "BMC Pharmacol Toxicol. 2024;25(1):57.",
    "doi:10.1186/s40360-024-00783-1.",
    sep = " "
  )
  vignette <- "Chen_2024_febuxostat"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "febuxostat", units = "mg", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "febuxostat", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "febuxostat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "febuxostat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters Vc/F and Vp/F through the paper's Eq. (2)",
        "log-linear (power) form log(P) = log(P_TV) + theta_COV * log(COV_i / COV_m),",
        "i.e. P = P_TV * (WT / WT_ref)^theta. The paper states that COV_m is the median",
        "of the covariate but never prints the pooled median weight; Table 1 gives",
        "per-study medians of 61.7, 63.4, 59.5 and 60.4 kg (n-weighted mean 61.3 kg) and",
        "the authors' own typical-subject simulation in Fig. 5 uses 60 kg. WT_ref is",
        "therefore assumed to be 60 kg here; see the vignette Errata.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    FED_HIGHFAT = list(
      description        = "Prandial state at dosing: 1 = single oral dose taken with a standardised high-fat high-calorie meal, 0 = fasted",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted; at least 10 h fast before dosing, first meal permitted 4 h post-dose)",
      notes              = paste(
        "The paper's covariate is named 'Food' (Methods: 'Food is 0 for the fasted state,",
        "and food equals 1 for the fed state'). The fed arm is the FDA-style high-fat",
        "high-calorie breakfast (800-1000 kcal; fried eggs in vegetable oil, 50 g fried",
        "bacon, 75 g toast, 15 g French fries, 240 mg whole milk per Supplementary",
        "Information section 1), consumed within 30 min of dosing, so the canonical",
        "FED_HIGHFAT indicator is used in preference to the generic FED. Enters ka1, ka2",
        "and Tlag1 through the paper's Eq. (1) exponential form",
        "log(P) = log(P_TV) + theta_COV * Food. Per-dose-record covariate.",
        sep = " "
      ),
      source_name        = "Food"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL/F and volume of distribution by stepwise covariate modelling (Methods, 'Population pharmacokinetic analysis') but not retained in the final model; no coefficient is reported."
    ),
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on CL/F and volume of distribution by stepwise covariate modelling but not retained in the final model; no coefficient is reported."
    ),
    BMI = list(
      description = "Body mass index at baseline",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened on CL/F and volume of distribution by stepwise covariate modelling but not retained in the final model; no coefficient is reported."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 128L,
    n_studies      = 2L,
    age_range      = "18-44 years",
    age_median     = "26.0-29.0 years across the four study arms (Table 1)",
    weight_range   = "45.4-79.0 kg",
    weight_median  = "59.5-63.4 kg across the four study arms (Table 1)",
    sex_female_pct = 36.7,
    race_ethnicity = "100% Chinese (Han-majority single-centre cohort; race was not modelled as a covariate)",
    disease_state  = "Healthy adult volunteers (no hyperuricaemia or gout)",
    dose_range     = "Single oral dose of 20 mg or 80 mg febuxostat tablet with 240 mL water, fasted or fed",
    regions        = "China (Center of Clinical Pharmacology, Second Affiliated Hospital of Zhejiang University School of Medicine, Hangzhou)",
    notes          = paste(
      "2455 plasma concentration records from 128 volunteers (81 male, 47 female) pooled",
      "from two open-label, single-dose, randomised, two-period crossover bioequivalence",
      "studies; only the reference-product periods were used in the popPK analysis.",
      "Arm sizes (Table 1): 20 mg fasting n = 36, 20 mg fed n = 31, 80 mg fasting n = 31,",
      "80 mg fed n = 30. Concentrations below the lower limit of quantification were",
      "omitted. Sampling ran from 30 min pre-dose to 48 h (20 mg study) or 36 h (80 mg",
      "study) post-dose. Estimation by SAEM in Monolix Suite 2023R1.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Structural parameters - Chen 2024 Table 3, "Final model estimates
    # (RSE%)" column. All values are typical-value (population) estimates
    # reported on the linear scale; ini() carries log(value) per the
    # nlmixr2lib log-transformed-PK convention.
    #
    # Structure (Results, "PopPK analysis" + Fig. 2): the absorbed dose is
    # split between two depots. Depot 1 receives fraction F1, is absorbed
    # first-order at ka1 after lag Tlag1; depot 2 receives the remaining
    # (1 - F1) and is absorbed first-order at ka2 after lag Tlag2. Both
    # feed a two-compartment disposition model with first-order
    # elimination.
    # =====================================================================

    # Absorption
    lka     <- log(4.02); label("First-depot absorption rate constant ka1 (1/h)")                  # Table 3: ka1 = 4.02 1/h (RSE 11.4%)
    lka2    <- log(4.34); label("Second-depot absorption rate constant ka2 (1/h)")                 # Table 3: ka2 = 4.34 1/h (RSE 15.6%)
    lfdepot <- log(0.51); label("Fraction of the dose absorbed via the first depot F1 (unitless)") # Table 3: F1 = 0.51 (RSE 3.7%); the remaining 1 - F1 enters the second depot
    ltlag   <- log(0.22); label("First-depot absorption lag time Tlag1 (h)")                       # Table 3: Tlag1 = 0.22 h (RSE 8.93%)
    ltlag2  <- log(0.92); label("Second-depot absorption lag time Tlag2 (h)")                      # Table 3: Tlag2 = 0.92 h (RSE 7.63%)

    # Disposition (apparent, i.e. /F - the model has no absolute
    # bioavailability parameter; F1 only splits the dose between depots)
    lcl <- log(6.14);  label("Apparent clearance CL/F (L/h)")                                      # Table 3: CL/F = 6.14 L/h (RSE 2.26%)
    lvc <- log(10.38); label("Apparent central volume Vc/F (L, 60 kg reference)")                  # Table 3: Vc/F = 10.38 L (RSE 2.07%)
    lq  <- log(1.93);  label("Apparent intercompartmental clearance Q/F (L/h)")                    # Table 3: Q/F = 1.93 L/h (RSE 4.25%)
    lvp <- log(10.22); label("Apparent peripheral volume Vp/F (L, 60 kg reference)")               # Table 3: Vp/F = 10.22 L (RSE 3.21%)

    # =====================================================================
    # Covariate effects.
    #
    # Categorical (Food), Eq. (1): log(P) = log(P_TV) + theta_COV * Food,
    #   i.e. P = P_TV * exp(theta * FED_HIGHFAT), Food = 0 fasted / 1 fed.
    #   Cross-check on the reported effect sizes (Discussion): "food reduced
    #   ka1 and ka2 by 87% and 66%, respectively" -- exp(-2.04) = 0.130
    #   (-87.0%) and exp(-1.08) = 0.340 (-66.0%). Both reproduce exactly,
    #   confirming the exponential form and the sign convention.
    #
    # Continuous (WT), Eq. (2):
    #   log(P) = log(P_TV) + theta_COV * log(COV_i / COV_m),
    #   i.e. P = P_TV * (WT / WT_ref)^theta -- a power (allometric-style)
    #   model, not an additive-linear one, despite the Methods wording
    #   "Linear functions were used for continuous covariates" (linear on
    #   the log-log scale). WT_ref = 60 kg is assumed; see covariateData.
    # =====================================================================
    e_fed_highfat_ka   <- -2.04; label("Effect of a high-fat meal on ka1 (log-scale shift)")       # Table 3: Food_ka1 = -2.04 (RSE 8%)
    e_fed_highfat_ka2  <- -1.08; label("Effect of a high-fat meal on ka2 (log-scale shift)")       # Table 3: Food_ka2 = -1.08 (RSE 19.9%)
    e_fed_highfat_tlag <- 0.73;  label("Effect of a high-fat meal on Tlag1 (log-scale shift)")     # Table 3: Food_Tlag1 = 0.73 (RSE 18.2%)
    e_wt_vc            <- 0.58;  label("Power exponent of body weight on Vc/F (unitless)")         # Table 3: WT_Vc/F = 0.58 (RSE 27.9%)
    e_wt_vp            <- 1.15;  label("Power exponent of body weight on Vp/F (unitless)")         # Table 3: WT_Vp/F = 1.15 (RSE 22.5%)

    # =====================================================================
    # Inter-individual variability - Chen 2024 Table 3, "IIV (CV%) (RSE%)"
    # column. Methods: "Lognormal distributions were assumed for all model
    # parameters, with inter-individual variability (IIV) modeled
    # exponentially", i.e. P_i = P_TV * exp(eta_i).
    #
    # The table reports IIV as a coefficient of variation, so the internal
    # nlmixr2 variance is omega^2 = log(CV^2 + 1) (Monolix's lognormal
    # CV convention). Worked values below. The alternative reading of the
    # column as 100 * omega is ruled out by the observed Cmax variability
    # in Table 2 (~28-33% CV across the four arms): omega^2 = log(CV^2 + 1)
    # reproduces it, while treating the column as 100 * omega inflates
    # simulated Cmax CV several-fold. See the vignette Errata.
    #
    # IIV is uncorrelated - the paper reports no off-diagonal omega terms
    # and no correlation matrix.
    #
    # CAVEAT on etalfdepot: F1 is a fraction bounded on (0, 1), but the
    # paper states a lognormal distribution for every parameter, which is
    # unbounded above. Encoded literally (as here), P(F1 > 1) = 1.8%, and
    # those subjects receive a negative dose fraction into depot2. Monolix
    # would more usually declare a fraction as logit-normal, but the paper
    # neither says so nor reports a logit-scale omega, so no such value is
    # invented here. Typical-value simulations (zeroRe() / omega = NA) are
    # unaffected. See the vignette Errata.
    # =====================================================================
    etalka     ~ 0.704247  # CV 101.11% (RSE 7.62%)  -> log(1.0111^2 + 1) = 0.704247 (omega 0.8392)
    etalka2    ~ 1.187690  # CV 150.98% (RSE 7.72%)  -> log(1.5098^2 + 1) = 1.187690 (omega 1.0898)
    etalfdepot ~ 0.102596  # CV  32.87% (RSE 7.24%)  -> log(0.3287^2 + 1) = 0.102596 (omega 0.3203)
    etaltlag   ~ 0.505737  # CV  81.13% (RSE 6.91%)  -> log(0.8113^2 + 1) = 0.505737 (omega 0.7112)
    etaltlag2  ~ 0.723643  # CV 103.05% (RSE 6.64%)  -> log(1.0305^2 + 1) = 0.723643 (omega 0.8507)
    etalcl     ~ 0.064008  # CV  25.71% (RSE 6.44%)  -> log(0.2571^2 + 1) = 0.064008 (omega 0.2530)
    etalvc     ~ 0.031331  # CV  17.84% (RSE 10.6%)  -> log(0.1784^2 + 1) = 0.031331 (omega 0.1770)
    etalq      ~ 0.133621  # CV  37.81% (RSE 11.1%)  -> log(0.3781^2 + 1) = 0.133621 (omega 0.3655)
    etalvp     ~ 0.058848  # CV  24.62% (RSE 10.9%)  -> log(0.2462^2 + 1) = 0.058848 (omega 0.2426)

    # =====================================================================
    # Residual error - Chen 2024 Table 3, "Combined error model
    # parameters". Methods: "Residual variability was assessed using
    # proportional, additional, and combined proportional and additional
    # error models"; the final model uses the combined form. The paper does
    # not state whether Monolix's combined1 (a + b*f) or combined2
    # (sqrt(a^2 + (b*f)^2)) variance form was used; the rxode2 default
    # combined2 form is used here. The additive term is negligible against
    # the observed concentration range (hundreds to thousands of ng/mL), so
    # the choice has no practical consequence. See the vignette Errata.
    # =====================================================================
    addSd  <- 2.14; label("Additive residual SD (ng/mL)")             # Table 3: Add. = 2.14 ng/mL (RSE 5.82%)
    propSd <- 0.12; label("Proportional residual SD (fraction)")      # Table 3: Prop. = 0.12 (RSE 2.27%)
  })

  model({
    # ---- Individual parameters ------------------------------------
    # Food (FED_HIGHFAT) acts multiplicatively on ka1, ka2 and Tlag1
    # via Eq. (1); body weight acts as a power function on the two
    # apparent volumes via Eq. (2), referenced to 60 kg.
    ka     <- exp(lka  + etalka)  * exp(e_fed_highfat_ka  * FED_HIGHFAT)
    ka2    <- exp(lka2 + etalka2) * exp(e_fed_highfat_ka2 * FED_HIGHFAT)
    tlag   <- exp(ltlag + etaltlag) * exp(e_fed_highfat_tlag * FED_HIGHFAT)
    tlag2  <- exp(ltlag2 + etaltlag2)
    fdepot <- exp(lfdepot + etalfdepot)

    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (WT / 60)^e_wt_vc
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp) * (WT / 60)^e_wt_vp

    # ---- Micro-constants ------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- ODE system (Fig. 2) --------------------------------------
    d/dt(depot)       <- -ka  * depot
    d/dt(depot2)      <- -ka2 * depot2
    d/dt(central)     <-  ka * depot + ka2 * depot2 -
                          kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ---- Dose split and lag times ---------------------------------
    # The same dose record is given to both depots; F1 and 1 - F1
    # partition it. Total bioavailability is 1 by construction, so
    # CL/F, Vc/F, Q/F and Vp/F remain apparent parameters.
    f(depot)     <- fdepot
    f(depot2)    <- 1 - fdepot
    alag(depot)  <- tlag
    alag(depot2) <- tlag2

    # ---- Observation ----------------------------------------------
    # central is in mg and vc in L, so central/vc is mg/L; the factor
    # 1000 converts to the reported ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
