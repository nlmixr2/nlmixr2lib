Hu_2026_utreloxastat <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model with an eight-transit-",
    "compartment first-order oral absorption chain and time-varying clearance",
    "for utreloxastat (PTC857), a 15-lipoxygenase inhibitor in development for",
    "amyotrophic lateral sclerosis, in healthy adult volunteers (Hu 2026, the",
    "first-in-human study: single ascending doses 100-1000 mg, multiple",
    "ascending doses 150-500 mg BID for 14 days, and a three-period food-effect",
    "crossover at 500 mg). Apparent clearance fell and the terminal half-life",
    "lengthened over repeated dosing (23 h on Day 1 to over 33 h on Day 14), so",
    "the authors screened twelve candidate time-varying clearance functions and",
    "selected the exponential-in-log-time form (their Model 8),",
    "CL(t) = CLTV * (1 + exp(-KTMCL * ln(t))), which is the canonical",
    "cl_exp_inf + cl_exp_component * exp(-kdes * .) decomposition with the",
    "decaying component's amplitude equal to CLTV and the decay running against",
    "natural-log time after first dose. Absorption is more than dose",
    "proportional: relative bioavailability carries a positive power effect of",
    "dose level, plus additive increases of 36% for low-fat and 57% for",
    "high-fat meals. Age, sex, body weight, BMI, race, ethnicity and hepatic /",
    "renal laboratory markers were screened and none was retained."
  )
  reference <- paste(
    "Hu Y, Gao L, Kong R. Population Pharmacokinetic Modeling Practice of",
    "Time-Varying Clearance: Insights from a First-in-Human Study Case.",
    "J Clin Pharmacol. 2026;66(3):e70171. doi:10.1002/jcph.70171"
  )
  vignette <- "Hu_2026_utreloxastat"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Absorption is an explicit Savic-style chain: the dose lands in `depot` and
  # passes through eight transit compartments before entering `central`, every
  # step at the same rate constant `ktr`. That is exactly the structure implied
  # by the Hu 2026 Table 3 footnote, "Mean Transit Time = (number of transit
  # compartment + 1) / transit constant = (8 + 1) / 4.10 = 2.20 h": nine
  # first-order transfers in series (depot plus eight transits) at 4.10 1/h.
  compartmentData <- list(
    depot       = list(analyte = "utreloxastat", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "utreloxastat", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "utreloxastat", units = "mg", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "utreloxastat", units = "mg", specimen = "administration site", verified = TRUE),
    transit4    = list(analyte = "utreloxastat", units = "mg", specimen = "administration site", verified = TRUE),
    transit5    = list(analyte = "utreloxastat", units = "mg", specimen = "administration site", verified = TRUE),
    transit6    = list(analyte = "utreloxastat", units = "mg", specimen = "administration site", verified = TRUE),
    transit7    = list(analyte = "utreloxastat", units = "mg", specimen = "administration site", verified = TRUE),
    transit8    = list(analyte = "utreloxastat", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "utreloxastat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "utreloxastat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DOSE_UTRELOXASTAT_MG = list(
      description        = "Administered utreloxastat dose per administration",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Amount of a single administration, not the daily total: a 250 mg BID subject carries 250.",
        "Enters relative bioavailability as the power term (DOSE_UTRELOXASTAT_MG / 250)^e_dose_fdepot.",
        "Hu 2026 states the power form (Methods, 'A power model was used for continuous covariates')",
        "and reports the exponent (Table 3, theta9 = 0.21) but never prints the normalizing dose.",
        "250 mg is adopted here because it is the reference dose of the paper's own covariate forest",
        "plot (Figure 4 / Table 5, 'Reference: 250 mg/BID/fasted') and the median dose of both the SAD",
        "and MAD cohorts (Table 1). The constant cancels out of every exposure ratio the paper reports,",
        "so no published number can falsify it, but it does set the absolute concentration scale --",
        "change it deliberately if a different normalization is preferred. Study range 100-1000 mg."
      ),
      source_name        = "Dose levels"
    ),
    FED_LOWFAT = list(
      description        = "Low-fat meal at the time of dosing",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste(
        "1 = dose taken with a low-fat meal, 0 = fasted (mutually exclusive with FED_HIGHFAT).",
        "Hu 2026 pooled three meal conditions into this one level: the protocol low-fat low-calorie",
        "meal of the food-effect part; the site-standardized meal given throughout the SAD and MAD",
        "parts (fat providing at least 25% of calories, which the authors describe as 'very close to",
        "the standard low-fat meal as defined by FDA guidance'); and a customized medium-fat meal",
        "whose effect was 'not statistically distinguishable' from low-fat (Introduction and",
        "Discussion). FED_LOWFAT = 1 is therefore the usual state in this study, not a special arm."
      ),
      source_name        = "Low-fat meals"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat high-calorie meal at the time of dosing",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste(
        "1 = dose taken with a high-fat high-calorie meal, 0 = fasted (mutually exclusive with",
        "FED_LOWFAT). Only the 500 mg three-period food-effect crossover contributed this level",
        "(Methods, 'PK Data Collection'); the paper defines the high-fat end of the meal scale as",
        "55-65 g fat."
      ),
      source_name        = "High-fat meals"
    )
  )

  # Screened during the stepwise covariate search and NOT retained: none reduced
  # the objective function by the forward-addition criterion (Results, 'Covariate
  # Impact on PK Exposures'; Model Structure section). No point estimates are
  # published for any of them, so they are documentation only.
  covariatesDataExcluded <- list(
    AGE        = list(description = "Age", units = "years", type = "continuous",
                      notes = "Screened, not retained. Study range 18.0-55.0 years (Table 1)."),
    SEXF       = list(description = "Female sex", units = "(binary)", type = "binary",
                      notes = "Screened, not retained. 47.1% female (Table 1)."),
    WT         = list(description = "Body weight", units = "kg", type = "continuous",
                      notes = "Screened, not retained. Study range 52.4-99.3 kg (Table 1). No allometric scaling is present in the final model."),
    BMI        = list(description = "Body mass index", units = "kg/m^2", type = "continuous",
                      notes = "Screened, not retained. Distribution not tabulated in Hu 2026."),
    RACE_WHITE = list(description = "White race indicator", units = "(binary)", type = "binary",
                      notes = "Screened, not retained. Hu 2026 names race and ethnicity among the screened factors but does not tabulate their distribution."),
    ALB        = list(description = "Serum albumin", units = "g/dL", type = "continuous",
                      notes = "Screened, not retained; hepatic-function marker."),
    ALP        = list(description = "Alkaline phosphatase", units = "U/L", type = "continuous",
                      notes = "Screened, not retained; hepatic-function marker."),
    ALT        = list(description = "Alanine aminotransferase", units = "U/L", type = "continuous",
                      notes = "Screened, not retained; hepatic-function marker."),
    AST        = list(description = "Aspartate aminotransferase", units = "U/L", type = "continuous",
                      notes = "Screened, not retained; hepatic-function marker."),
    TBILI      = list(description = "Total bilirubin", units = "mg/dL", type = "continuous",
                      notes = "Screened, not retained; hepatic-function marker."),
    CREAT      = list(description = "Serum creatinine", units = "mg/dL", type = "continuous",
                      notes = "Screened, not retained; renal-function marker. Elimination of utreloxastat is primarily metabolic with negligible renal clearance (Introduction).")
  )

  population <- list(
    species        = "human",
    n_subjects     = 68L,
    n_studies      = 1L,
    age_range      = "18.0-55.0 years",
    age_median     = "29.5 years",
    weight_range   = "52.4-99.3 kg",
    weight_median  = "75.2 kg",
    sex_female_pct = 47.1,
    race_ethnicity = "Not tabulated in Hu 2026; race and ethnicity were screened as covariates and not retained.",
    disease_state  = "Healthy adult volunteers.",
    dose_range     = paste(
      "Oral utreloxastat in three parts of one first-in-human study: single ascending doses of 100,",
      "250, 500 and 1000 mg with a low-fat breakfast (SAD, n = 32, 8 per dose level); multiple",
      "ascending doses of 150, 250 and 500 mg twice daily for 14 days with low-fat meals (MAD,",
      "n = 24); and a three-period single-dose 500 mg food-effect crossover under fasted, high-fat",
      "high-calorie and low-fat low-calorie conditions (FE, n = 12)."
    ),
    regions        = "Not specified.",
    notes          = paste(
      "Demographics from Hu 2026 Table 1 (overall column). 1463 quantifiable PK observations;",
      "concentrations below the quantitation limit were discarded following the M1 method.",
      "Sampling was intensive to 72 h post dose in the SAD and FE parts, to 24 h after the Day 1",
      "dose and to 72 h after the Day 14 morning dose in the MAD part, plus a single Day 7 predose",
      "sample (Table S3). The near-absence of data between Day 1 and Day 14 is the paper's stated",
      "reason that the curvature of the clearance-versus-time function is only weakly identified."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters. All values are Hu 2026 Table 3, "PK Parameter
    # Estimates of the Final PopPK Model (Model 8)", and are reproduced
    # unchanged in the supplement's Table S1 / S2 reference column.
    # Every disposition parameter is apparent (per unit bioavailability):
    # the paper never estimated absolute F.
    # ---------------------------------------------------------------------

    # Table 3 names theta1 "Ka", but the Table 3 footnote uses that same value
    # as the chain's "transit constant" in MTT = (8 + 1) / 4.10 = 2.20 h, so it
    # is the canonical ktr: one rate constant governs the depot, all eight
    # transit transfers, and entry into central.
    lktr <- log(4.10);     label("Transit-chain rate constant, also the terminal absorption rate into central (ktr, 1/h)")  # Table 3 theta1 (Ka = 4.10 1/h, ASE 0.21, %RSE 5.10, 90% CI 3.75-4.44); footnote: MTT = (8 + 1) / 4.10 = 2.20 h
    lvc  <- log(362.00);   label("Apparent central volume of distribution (V/F, L)")                                        # Table 3 theta2 (V/F = 362.00 L, ASE 21.27, %RSE 5.88, 90% CI 326.89-397.10)
    lvp  <- log(1827.15);  label("Apparent peripheral volume of distribution (V2/F, L)")                                    # Table 3 theta4 (V2/F = 1827.15 L, ASE 140.19, %RSE 7.67, 90% CI 1595.84-2058.46)
    lq   <- log(37.39);    label("Apparent intercompartmental clearance (Q/F, L/h)")                                        # Table 3 theta5 (Q/F = 37.39 L/h, ASE 2.35, %RSE 6.31, 90% CI 33.50-41.28)

    # ---------------------------------------------------------------------
    # Time-varying clearance, Hu 2026 Table 2 Model 8:
    #
    #     CL_TV,t = CL_TV * (1 + exp(-KTMCL * ln(t)))
    #
    # written here in the nlmixr2lib canonical decomposition
    # CL(t) = cl_exp_inf + cl_exp_component * exp(-cl_exp_kdes * .), which for
    # Model 8 has cl_exp_component == cl_exp_inf (Model 8 is the FTMCL-free
    # member of the Table 2 family -- compare Model 6, which carries an
    # estimated FTMCL amplitude; Table S1 lists KTMCL alone for Model 8) and
    # runs the decay against ln(t) rather than t.
    #
    # Arithmetic check against the paper's own statement that CL/F at the
    # plateau is "approximately 71 to 74 L/h" for every one of the twelve
    # candidate models (Results, "Model Structure and Pharmacokinetic
    # Parameters", referring to the blue dashed line in Table 2):
    #     t = 336 h (Day 14):  46.85 * (1 + 336^-0.096) = 73.6 L/h
    #     t = 672 h (4 weeks): 46.85 * (1 + 672^-0.096) = 71.9 L/h
    # which is where the paper says the curve reaches its plateau.
    # ---------------------------------------------------------------------
    lcl_exp_inf  <- log(46.85);  label("Asymptotic apparent clearance, the paper's CL_TV (CL/F, L/h)")                      # Table 3 theta3 (CL/F = 46.85 L/h, ASE 2.98, %RSE 6.34, 90% CI 41.93-51.76)
    lcl_exp_kdes <- log(0.096);  label("Decay coefficient of the clearance component against ln(time), paper KTMCL (unitless)")  # Table 3 theta10 (KTMCL = 0.096, ASE 0.022, %RSE 23.71, 90% CI 0.059-0.13); Table 2 Model 8. Table 3 labels KTMCL "h-1", but it multiplies the dimensionless ln(t) in Model 8 and so is unitless.

    # ---------------------------------------------------------------------
    # Relative bioavailability (Frel). Hu 2026 Methods: "A power model was
    # used for continuous covariates, and a binominal model was used for
    # categorical covariates." The additive (binomial) form of the meal terms
    # is confirmed by the paper's own forest-plot output: Table 5 reports
    # steady-state AUC ratios of 1.37 for low-fat and 1.56 for high-fat versus
    # fasted, matching 1 + 0.36 = 1.36 and 1 + 0.57 = 1.57 to Monte-Carlo
    # noise. The power form of the dose term is confirmed the same way:
    # exposure scales as Dose^(1 + 0.21), giving (150/250)^1.21 = 0.54 and
    # (500/250)^1.21 = 2.31 against the paper's 0.53 and 2.28.
    # ---------------------------------------------------------------------
    lfdepot              <- fixed(log(1));  label("Relative bioavailability at the reference condition, 250 mg fasted (Frel, unitless)")  # Structural anchor: Hu 2026 estimated only relative bioavailability, so absolute F is not identifiable
    e_fed_lowfat_fdepot  <- 0.36;           label("Additive effect of a low-fat meal on Frel (unitless)")                   # Table 3 theta7 (0.36, ASE 0.071, %RSE 19.48, 90% CI 0.25-0.48); Results: "food increases drug absorption by up to 36% (90%CI: 25% to 48%) for low-fat meals"

    # Hu 2026 states this effect inconsistently. Table 3 theta8 and the Results
    # section both give 57% (90% CI 42% to 72%); the Discussion's closing summary
    # instead reads "72% (90% CI: 42% to 72%) for high-fat meals", repeating the
    # confidence interval's upper bound in place of the point estimate. 0.57 is
    # adopted, and the paper's own forest-plot output adjudicates it
    # independently: Table 5 reports a steady-state AUC ratio of 1.56 for
    # high-fat versus fasted, which matches 1 + 0.57 = 1.57 and rules out
    # 1 + 0.72 = 1.72. See the vignette's Assumptions and deviations section.
    e_fed_highfat_fdepot <- 0.57;           label("Additive effect of a high-fat meal on Frel (unitless)")                  # Table 3 theta8 (0.57, ASE 0.092, %RSE 16.19, 90% CI 0.42-0.72); Results: "57% (90%CI: 42% to 72%) for high-fat meals"
    e_dose_fdepot        <- 0.21;           label("Power exponent on (dose / 250 mg) for Frel (unitless)")                  # Table 3 theta9 (0.21, ASE 0.036, %RSE 17.06, 90% CI 0.15-0.27); Results: "the positive power value for dose levels on relative absorption fraction ... indicated that the absorption fraction increased at higher dose levels"

    # ---------------------------------------------------------------------
    # Inter-individual variability. Table 3 reports the two retained etas as
    # variances on the log scale: its "(%CV)" column is
    # sqrt(exp(omega^2) - 1) * 100, which reproduces 16.53% from 0.027
    # (sqrt(exp(0.027) - 1) = 0.1654) and 42.26% from the unrounded IIV-Ka
    # value that Table 3 prints as 0.16. Reading them as SDs instead would
    # make the Ka %CV about 16%, contradicting the printed 42.26%.
    # ---------------------------------------------------------------------
    etalktr        ~ 0.16    # Table 3 ETA/IIV row 1 (IIV-Ka = 0.16, ASE 0.030, %CV 42.26, shrinkage 1.21%)
    etalcl_exp_inf ~ 0.027   # Table 3 ETA/IIV row 2 (IIV-CL/F = 0.027, ASE 0.0059, %CV 16.53, shrinkage 11.24%)

    # ---------------------------------------------------------------------
    # Residual error. Methods give the observation equation as
    #     Ln(Y) = Ln(f(theta, eta, t)) + epsilon
    # i.e. additive on the natural-log scale, which is nlmixr2's lnorm().
    #
    # Table 3 row 6 ("Additive residual") is 0.46 with ASE 0.0091, and it is a
    # STANDARD DEVIATION, not a variance. Proof from the reported precision
    # alone: with 1463 observations the Cramer-Rao lower bound on the standard
    # error of a variance estimate is sigma^2 * sqrt(2/n) = 0.46 * 0.037 =
    # 0.017, so an ASE of 0.0091 is impossible for a variance. The same bound
    # for a standard-deviation estimate is sigma / sqrt(2n) = 0.46 / 54.1 =
    # 0.0085, which 0.0091 clears. The two IIV rows pass the mirror-image test
    # and are variances (see the IIV block above).
    # ---------------------------------------------------------------------
    expSd <- 0.46;  label("Log-scale residual standard deviation (unitless; approx 46% CV)")  # Table 3 theta6 (Additive residual = 0.46, ASE 0.0091, %RSE 1.97, 90% CI 0.45-0.48); Methods equation Ln(Y) = Ln(f) + epsilon
  })

  model({
    # -------------------------------------------------------------------
    # Numerical guard on the ln(t) singularity of Model 8. Hu 2026 flags it
    # in the Discussion: "applying a natural logarithm transformation of
    # predose time (time = 0) results in undefined values (approaching
    # infinity)". The paper does not say how the fit handled it, so the
    # time argument is floored here at 0.01 h (36 s). The floor is
    # numerically inert rather than a modelling choice: after 0.01 h the
    # dose has barely begun to leave `depot` and has to cross eight further
    # transit compartments, so `central` is ~1e-15 of the dose and the value
    # of `cl` over that interval cannot move the solution. The vignette
    # demonstrates the insensitivity across floors of 1e-4, 0.01 and 0.5 h.
    #
    # `time` is rxode2 model time, which with dosing started at t = 0 is the
    # paper's "longitudinal Time" (time after first dose) -- note this is
    # time after FIRST dose, not time after the most recent dose: Table 2
    # spells TAD separately, and uses it only in Model 4.
    # -------------------------------------------------------------------
    tcl <- max(time, 0.01)   # h

    # 1. Individual parameters
    ktr <- exp(lktr + etalktr)
    vc  <- exp(lvc)
    vp  <- exp(lvp)
    q   <- exp(lq)

    # 2. Time-varying clearance, Table 2 Model 8. cl_exp_component equals
    #    cl_exp_inf because Model 8 carries no separate FTMCL amplitude, so
    #    total clearance runs from 2 * CL_TV at t = 1 h down towards CL_TV.
    cl_exp_inf       <- exp(lcl_exp_inf + etalcl_exp_inf)
    cl_exp_component <- cl_exp_inf
    cl_exp_kdes      <- exp(lcl_exp_kdes)
    cl <- cl_exp_inf + cl_exp_component * exp(-cl_exp_kdes * log(tcl))

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. ODE system: depot plus eight transit compartments in series into
    #    central, all at ktr, then two-compartment disposition.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3 - ktr * transit4
    d/dt(transit5)    <-  ktr * transit4 - ktr * transit5
    d/dt(transit6)    <-  ktr * transit5 - ktr * transit6
    d/dt(transit7)    <-  ktr * transit6 - ktr * transit7
    d/dt(transit8)    <-  ktr * transit7 - ktr * transit8
    d/dt(central)     <-  ktr * transit8 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                   k12 * central - k21 * peripheral1

    # 5. Relative bioavailability: additive meal terms times a power term on
    #    dose level, normalized at the paper's 250 mg reference condition.
    f(depot) <- exp(lfdepot) *
      (1 + e_fed_lowfat_fdepot * FED_LOWFAT + e_fed_highfat_fdepot * FED_HIGHFAT) *
      (DOSE_UTRELOXASTAT_MG / 250)^e_dose_fdepot

    # 6. Observation. Doses are in mg and vc is in L, so central/vc is mg/L =
    #    ug/mL; the factor 1000 converts to the ng/mL that Hu 2026 uses
    #    throughout (Table 4 column headers "Cmax (ng mL-1)").
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
