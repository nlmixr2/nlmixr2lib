Koenig_2025_cefiderocol <- function() {
  description <- "Two-compartment population PK model for cefiderocol in hospitalized adults with cystic fibrosis during an acute pulmonary exacerbation, fitted nonparametrically with the Pmetrics NPAG algorithm and parameterised by clearance, central volume and the intercompartmental micro-rate constants k12 and k21"
  reference <- paste(
    "Koenig C, Monogue ML, Shields RK, Sakon CM, Fratoni AJ, Roenfanz HF,",
    "Finklea JD, Pope JS, Nicolau DP, Kuti JL.",
    "Cefiderocol pharmacokinetics during acute pulmonary exacerbations in",
    "hospitalized adult persons with cystic fibrosis.",
    "Antimicrob Agents Chemother. 2025;69(1):e01539-24.",
    "doi:10.1128/aac.01539-24",
    sep = " "
  )
  vignette <- "Koenig_2025_cefiderocol"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The FINAL model carries no covariate effect on any disposition parameter
  # (see covariatesDataExcluded for the screen that rejected them). FU is the
  # one covariate the model consumes, and it acts only on the free-concentration
  # observable that drives the paper's target-attainment analysis.
  covariateData <- list(
    FU = list(
      description = paste(
        "Fraction of cefiderocol unbound in plasma. Measured per subject by",
        "ultrafiltration (Centrifree) of an ex-vivo plasma sample drawn 3 h",
        "after the start of the final infusion, assayed in triplicate.",
        "Enters ONLY the free-concentration observable Cfree, which is what",
        "the paper's %fT > MIC target-attainment analysis is computed on; it",
        "does not scale clearance, volume or any transfer constant."
      ),
      units = "fraction (unitless; 0.52 = 52% unbound)",
      type = "continuous",
      source_name = "Protein binding (%)",
      notes = paste(
        "Koenig 2025 reports PROTEIN BINDING, not fraction unbound. Table 3",
        "gives a per-subject protein-binding percentage (mean 48%, range",
        "35-57%), so convert on ingestion with FU = 1 - PB/100: the cohort",
        "mean 48% bound is FU = 0.52 and the observed 35-57% range is",
        "FU 0.43-0.65. The Methods formula is printed as",
        "'protein binding [%] = 1 - CPFF/CPlasma * 100', which is missing a",
        "bracket; it means (1 - CPFF/CPlasma) * 100, as confirmed by the",
        "Table 3 values lying between 35 and 57.",
        "Per-subject values (Table 3, by subject ID): 1 = 54%, 2 = 45%,",
        "3 = 55%, 4 = 38%, 5 = 47%, 6 = 57%, 7 = 56%, 8 = 49%, 9 = 35%."
      )
    )
  )

  # Covariates that Koenig 2025 screened and REJECTED. Both were tested by
  # linear regression against the individual Bayesian parameter estimates and,
  # for eGFR, by three nested covariate models; none reduced the AIC by more
  # than 2 relative to the base two-compartment model, so the base model was
  # selected as final (Results 'Population pharmacokinetic analyses' and
  # supplemental 'Model development process' table).
  covariatesDataExcluded <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance by the Cockcroft-Gault equation, RAW mL/min and",
        "NOT BSA-normalized, computed on IDEAL body weight (adjusted body",
        "weight where total body weight exceeded ideal by more than 20%, and",
        "total body weight where it fell below ideal) per the supplemental",
        "'Calculation of CrCL using Cockcroft-Gault' section. The paper calls",
        "this column eGFR throughout."
      ),
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Screened and NOT retained, despite being the one statistically",
        "significant correlate found: supplemental Figure S1 regresses CrCL on",
        "total body clearance and reports slope 0.0399 L/h per mL/min",
        "(SE 0.014, t = 2.847, p = 0.0248, 95% CI 0.007-0.073) on an",
        "intercept of 0.9974 L/h. The Discussion quotes this slope as 0.039.",
        "Three nested covariate models were then fitted and all failed to",
        "improve the fit (supplemental 'Model development process' table):",
        "CL = CLi + CLs*CrCL gave AIC 509, CL = CL0*(CrCL/117) gave AIC 501",
        "and CL = CL0*((CrCL/117)^0.75) gave AIC 501, against AIC 501 for the",
        "base model. The authors attribute the null result to the narrow",
        "observed range (71-164 mL/min), noting that other cefiderocol",
        "cohorts spanning wider renal function DO retain eGFR on CL.",
        "Cohort values: mean 117 +/- 24 mL/min, range 71-164 (Table 1);",
        "eGFR < 60 mL/min was an exclusion criterion, so this model carries",
        "no information whatever about renal impairment."
      )
    ),
    WT = list(
      description = "Total body weight at baseline",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened and NOT retained. Supplemental Figure S2 regresses body",
        "weight on central volume and reports slope 0.0427 L/kg (SE 0.137,",
        "t = 0.312, p = 0.764, 95% CI -0.281 to 0.366) on an intercept of",
        "3.1682 L -- i.e. indistinguishable from zero. The Results quote",
        "p = 0.76. No allometric term was imposed either, so this model is",
        "NOT weight-scaled; the Discussion attributes the null result to the",
        "narrow observed range (45-78 kg) and cautions explicitly against",
        "extrapolating to pwCF of higher body weight.",
        "Cohort values: mean 62 +/- 10 kg, range 45-78 (Table 1)."
      )
    )
  )

  compartmentData <- list(
    central = list(analyte = "cefiderocol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefiderocol", units = "mg", specimen = "plasma", verified = FALSE)
  )

  population <- list(
    species = "human",
    n_subjects = 9,
    n_studies = 1,
    age_range = "22-58 years",
    age_median = "33 years (mean, SD 11)",
    weight_range = "45-78 kg",
    weight_median = "62 kg (mean, SD 10)",
    height_range = "157-180 cm (mean 170, SD 7)",
    sex_female_pct = 11.1,
    disease_state = "cystic fibrosis, hospitalized with an acute pulmonary exacerbation",
    dose_range = paste(
      "2 g cefiderocol as a 3 h prolonged intravenous infusion; q8h in 6",
      "subjects and q6h in 3 subjects, the frequency set by the approved",
      "label according to Cockcroft-Gault eGFR"
    ),
    renal_function = paste(
      "normal to augmented: Cockcroft-Gault eGFR 71-164 mL/min (mean 117,",
      "SD 24). eGFR < 60 mL/min, any renal replacement therapy and",
      "haemodialysis were exclusion criteria, and 3 of 9 subjects had",
      "augmented renal function (eGFR > 120 mL/min)"
    ),
    co_medication = paste(
      "all subjects received standard-of-care antibiotics at the discretion of",
      "the attending provider, with cefiderocol added on top; 6 of 9 subjects",
      "were on CFTR modulator therapy during the study"
    ),
    regions = "United States (prospective study across 4 sites)",
    notes = paste(
      "Demographics from Koenig 2025 Table 1, which publishes every subject",
      "individually. Ten subjects were enrolled; one was excluded for a",
      "positive pregnancy test before receiving study drug, leaving 9 who",
      "completed. A total of 80 plasma samples were drawn at 0 (pre-dose),",
      "1.5, 3, 3.25, 3.5, 4, 5, 6 and 8 h after the start of the final dose,",
      "which followed at least 3 prior doses so that sampling is at steady",
      "state. Fitted in Pmetrics 2.1.1 (NPAG) under R 4.3.2. The",
      "two-compartment model was selected on AIC 501 against 548 for a",
      "one-compartment model; goodness-of-fit R^2 was 0.96 with a fitted",
      "gamma of 1.8. The lower limit of quantification was 0.1 mg/L in",
      "plasma and 0.05 mg/L in protein-free filtrate.",
      "This analysis was presented in preliminary form as poster P-1227 at",
      "IDWeek 2024 (abstract ofae631.1409, Open Forum Infect Dis",
      "2025;12(Suppl 1):S784). The poster's estimates DIFFER from the final",
      "paper's and are superseded by them -- see the vignette Errata."
    )
  )

  ini({
    # Structural parameters: Koenig 2025 Table 2, 'Population estimate
    # (mean +/- SD)' column. These are the arithmetic MEAN of the nonparametric
    # joint density, and each reproduces the mean of the 9 published individual
    # MAP Bayesian estimates in the same table to 3 significant figures
    # (CL 5.653, Vc 5.810, k12 4.284, k21 2.249).
    #
    # Encoded as the MEDIAN of a lognormal, following the convention used for
    # every other nonparametric (NPAG) model in this library, so that a
    # typical-value simulation reproduces the published number exactly. The
    # paper reports a mean rather than a median, so a stochastic cohort's
    # arithmetic mean sits above these values by exp(omega^2 / 2). See the
    # vignette Assumptions section, which quantifies the offset per parameter.
    lcl  <- log(5.66); label("Total body clearance (L/h)")                            # Table 2: CL 5.66 +/- 1.28
    lvc  <- log(5.81); label("Central volume of distribution (L)")                     # Table 2: Vc 5.81 +/- 3.52
    lk12 <- log(4.29); label("Central-to-peripheral transfer rate constant (1/h)")     # Table 2: K12 4.29 +/- 3.46
    lk21 <- log(2.25); label("Peripheral-to-central transfer rate constant (1/h)")     # Table 2: K21 2.25 +/- 2.76

    # Interindividual variability. NPAG places every model parameter in the
    # joint density, so Table 2's SD column IS the between-subject SD on the
    # natural scale. Converted to lognormal variances with the standard
    # identity omega^2 = log(CV^2 + 1), which preserves the reported central
    # value as the distribution median.
    #
    # The variances below are computed from the DIAGONAL OF THE SUPPLEMENT
    # COVARIANCE MATRIX rather than from the Table 2 SD column, because the
    # matrix carries the unrounded variances while Table 2 rounds their square
    # roots to three significant figures. The two agree to that rounding
    # (sqrt of 1.647, 12.409, 11.966, 7.592 is 1.283, 3.523, 3.459, 2.755
    # against the printed 1.28, 3.52, 3.46, 2.76), so this is a precision
    # choice, not a different source:
    #   CL   omega^2 = log(1 +  1.647 / 5.66^2) = 0.05013  (CV 22.7%)
    #   Vc   omega^2 = log(1 + 12.409 / 5.81^2) = 0.31306  (CV 60.6%)
    #   k12  omega^2 = log(1 + 11.966 / 4.29^2) = 0.50089  (CV 80.6%)
    #   k21  omega^2 = log(1 +  7.592 / 2.25^2) = 0.91615  (CV 122.5%)
    #
    # These etas are INDEPENDENT even though the paper does publish the full
    # 4x4 covariance matrix of the joint density (supplemental 'Covariance
    # Matrix of final pharmacokinetic model'; its sqrt(diag) reproduces all
    # four reported SDs to 3 significant figures, confirming it is the matrix
    # behind Table 2). That matrix CANNOT be carried over, because the
    # multivariate lognormal it would require does not exist: the reported
    # Vc-k12 correlation of -0.724 lies below the lognormal-feasible lower
    # bound of -0.669 at these CVs, so the moment-matched log-scale matrix is
    # indefinite (smallest eigenvalue -0.085) and chol() fails on it. The
    # other five pairs are individually feasible, but a block omega must be
    # jointly positive definite. The published matrix and this arithmetic are
    # reproduced in the vignette so a user can audit and, if they wish,
    # project it themselves.
    etalcl  ~ 0.05013; label("IIV on clearance (log-scale variance)")
    etalvc  ~ 0.31306; label("IIV on central volume (log-scale variance)")
    etalk12 ~ 0.50089; label("IIV on k12 (log-scale variance)")
    etalk21 ~ 0.91615; label("IIV on k21 (log-scale variance)")

    # Residual error. Pmetrics weights each observation by its assay SD,
    # modelled as the polynomial SD = C0 + C1*[obs] + C2*[obs]^2 + C3*[obs]^3
    # and then scaled by a fitted gamma multiplier. Koenig 2025 Methods
    # 'Pharmacokinetic analyses' gives C0 = 0.0068, C1 = 0.0585, C2 = 0 and
    # C3 = 0, and the Results report the fitted gamma as 1.8, so the total
    # residual SD is 1.8 * (0.0068 + 0.0585 * Cc) = 0.01224 + 0.10530 * Cc.
    # Because C2 = C3 = 0 the polynomial is LINEAR in the observation, which is
    # nlmixr2's combined1() error structure -- Pmetrics sums the additive and
    # proportional parts directly, whereas nlmixr2's default combines them in
    # quadrature.
    #
    # Deliberately NOT fixed(): although C0 and C1 are stated assay constants
    # derived from the inter-run CVs of the LC/MS-MS method, the gamma that
    # multiplies both of them was estimated, so neither product is a
    # held-constant value.
    addSd  <- 0.01224; label("Additive residual SD (mg/L); gamma 1.8 x C0 0.0068")        # Methods: C0 = 0.0068, gamma = 1.8
    propSd <- 0.10530; label("Proportional residual SD (fraction); gamma 1.8 x C1 0.0585") # Methods: C1 = 0.0585, gamma = 1.8
  })
  model({
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)

    # Koenig 2025 parameterises the model as CL, Vc, k12 and k21 (Table 2), so
    # the elimination micro-constant and the peripheral compartment's macro
    # parameters are derived rather than estimated. Q and Vp are exact
    # algebraic identities of the published micro-constants, NOT extra
    # information: Q = k12 * Vc and Vp = Q / k21. At the population estimates
    # they are Q = 24.9 L/h and Vp = 11.1 L.
    #
    # Defining q and vp here is LOAD-BEARING, not cosmetic. rxode2 5.1.7
    # inspects an rxUi model for a recognisable linear-compartment
    # parameterisation and, when it finds one, solves it with its analytic
    # linCmt kernel and IGNORES the d/dt() right-hand sides entirely. A model
    # that defines cl and vc but NOT q and vp is recognised as a
    # ONE-compartment system, so the peripheral compartment below is silently
    # discarded: peripheral1 disappears from the solve output, central
    # plateaus at rate/kel during an infusion and the concentration decays
    # mono-exponentially at kel with no distribution phase. Supplying q and vp
    # makes the kernel solve the correct two-compartment system, which was
    # verified to reproduce the explicit-ODE solution exactly. Note that
    # AUC over a dosing interval is dose/CL either way, so an AUC-based check
    # CANNOT detect this -- see the vignette, which gates on the trough.
    q   <- k12 * vc
    vp  <- q / k21
    kel <- cl / vc

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Dose in mg and vc in L, so central / vc is mg/L (== ug/mL).
    Cc <- central / vc

    # Free (unbound) concentration. This is the quantity the paper's
    # %fT > MIC target-attainment analysis is evaluated against; only total
    # plasma cefiderocol was assayed and fitted, so the residual error below
    # attaches to Cc and not to Cfree.
    Cfree <- FU * Cc

    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
