Bustinduy_2016_praziquantel <- function() {
  description <- paste(
    "Two-compartment oral population PK model with a first-order absorption compartment and an",
    "absorption lag for TOTAL praziquantel (R-PZQ plus S-PZQ racemate) in 59 Ugandan children aged",
    "3-8 years with egg-patent Schistosoma mansoni intestinal schistosomiasis, given a single oral",
    "dose of 40 or 60 mg/kg after a standardised breakfast. Fitted non-parametrically with the NPAG",
    "algorithm in Pmetrics 1.2.6. The authors describe the structure as three compartments (gut,",
    "central, peripheral); the gut is the absorption compartment, so in nlmixr2 terms this is a",
    "two-compartment disposition model with a depot. Clearance, central volume and the two",
    "intercompartmental rate constants Kcp (k12) and Kpc (k21) are primary parameters, each carrying",
    "its own inter-individual variability, as is the absorption lag. Oral bioavailability was NOT",
    "estimated, so clearance and volume are apparent (CL/F 608 L/h, V/F 474 L) and lfdepot is a",
    "fixed unit anchor. No covariate was retained: weight, age and sex were screened against the",
    "Bayesian posterior estimates and none showed a relationship, so the base model is the final",
    "model. Residual unexplained variability is carried as fixed(0) because the Pmetrics assay-error",
    "polynomial was never published. The authors also fitted a separate enantiomer-specific (R and S)",
    "structural model, but reported its parameters as 'data not shown', so only the total-PZQ model",
    "is extractable. The companion logistic exposure-response model for parasitological cure is NOT",
    "encoded here because its intercept is unreported; see the vignette Errata.",
    sep = " "
  )
  reference <- paste(
    "Bustinduy AL, Waterhouse D, de Sousa-Figueiredo JC, Roberts SA, Atuhaire A, Van Dam GJ,",
    "Corstjens PLAM, Scott JT, Stanton MC, Kabatereine NB, Ward S, Hope WW, Stothard JR.",
    "Population pharmacokinetics and pharmacodynamics of praziquantel in Ugandan children with",
    "intestinal schistosomiasis: higher dosages are required for maximal efficacy.",
    "mBio. 2016;7(4):e00227-16. doi:10.1128/mBio.00227-16. PMCID: PMC4992966.",
    sep = " "
  )
  vignette <- "Bustinduy_2016_praziquantel"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot       = list(analyte = "praziquantel", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "praziquantel", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "praziquantel", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened but NOT retained. Methods (Pharmacokinetic population analyses): 'Potential",
        "relationships between each model parameter and covariates (e.g., weight, age, and gender)",
        "were explored by plotting the Bayesian posterior estimate for the parameter against the",
        "covariate'; Results (PZQ pharmacokinetics): 'There was no relationship between any of the",
        "Bayesian estimates of the parameter values and any of the available covariates. Thus, the",
        "standard base model was used.' Body weight nonetheless sets the administered dose, which is",
        "prescribed in mg/kg. The Discussion notes that an allometric 0.75 exponent on clearance",
        "'may well be the most appropriate scaling function but would require PK data from children",
        "with a wider range of weights', i.e. it was not estimable in the 15-34 kg cohort studied.",
        sep = " "
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened against the Bayesian posterior parameter estimates but not retained (Methods;",
        "Results). Cohort range 3-8 years, split by the authors into preschool (3-5 years, n = 17)",
        "and school-aged (6-8 years, n = 43) strata in Table 1.",
        sep = " "
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = paste(
        "Screened as 'gender' against the Bayesian posterior parameter estimates but not retained in",
        "the PK model (Methods; Results). 38 of 60 enrolled children (63.3%) were female. Female sex",
        "does appear in the separate logistic cure-rate analysis (Table 5, unadjusted OR 3.18,",
        "95% CI 0.89-11.3, P = 0.068), which is not encoded in this model file.",
        sep = " "
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 59L,
    n_studies        = 1L,
    age_range        = "3-8 years",
    age_median       = "Means by arm: 6.4 years (40 mg/kg) and 6.3 years (60 mg/kg); median not reported",
    weight_range     = "15-34 kg",
    weight_median    = "Means by arm: 21.8 kg (40 mg/kg) and 23.0 kg (60 mg/kg); median not reported",
    sex_female_pct   = 63.3,
    race_ethnicity   = "Not reported (Ugandan children from Bugoigo and Walukuba villages on the shore of Lake Albert)",
    disease_state    = paste(
      "Egg-patent Schistosoma mansoni intestinal schistosomiasis in a hyperendemic setting.",
      "Baseline intensity by Kato-Katz: heavy (>400 epg) 23 children, medium (100-399 epg) 17,",
      "light (1-100 epg) 19, and 1 child egg-negative but CCA-positive. Arithmetic mean baseline",
      "egg count 950.4 epg (40 mg/kg arm) and 491.0 epg (60 mg/kg arm). 37 of 59 were anaemic",
      "(haemoglobin < 11.5 g/dL) and 45 of 59 had a positive malaria rapid diagnostic test; no child",
      "tested positive for HIV. All children had received 40 mg/kg praziquantel at least one year",
      "before the study and had tolerated it.",
      sep = " "
    ),
    dose_range       = paste(
      "Single oral dose of praziquantel, randomised 1:1 to 40 mg/kg (n = 30) or 60 mg/kg (n = 30),",
      "administered after a breakfast of local foods because food increases PZQ bioavailability.",
      sep = " "
    ),
    regions          = "Uganda (Bugoigo and Walukuba villages, Lake Albert)",
    notes            = paste(
      "Baseline demographics from Bustinduy 2016 Table 1 (60 enrolled). 59 children completed the",
      "protocol and entered the population PK analysis (one child was withdrawn 6 h into sampling",
      "with falciparum malaria); 58 were available at the 24-day pharmacodynamic follow-up. Venous",
      "samples were drawn at 0, 1, 2, 4, 6, 12 and 24 h. R-PZQ and S-PZQ were quantified separately",
      "by LC-MS/MS, linear from 5 to 1500 ng/mL for each isomer; this model was fitted to their sum",
      "(total PZQ). Concentrations are in mg/L on the Fig. 3 observed-versus-predicted axes, while",
      "the Fig. 2 individual profiles and the Table 2 Cmax column are in ng/mL -- the Table 2 header",
      "prints 'ug/ml' for Cmax, which is a unit typo; see the vignette Errata.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # All structural values are the MEAN of the NPAG non-parametric parameter
    # distribution reported in Bustinduy 2016 Table 3 (59 children, both dose
    # arms pooled). Table 3 also prints a median for each parameter; the median
    # is recorded in the trailing comments. The mean is used here for
    # consistency with the sibling Pmetrics/NPAG extraction
    # Setiawan_2023_sulbactam, and because for CL and V the reported mean lies
    # BELOW the reported median, so the non-parametric distribution is
    # left-skewed and no log-normal can reproduce both -- see vignette Errata.
    #
    # Oral bioavailability was not estimated, so CL and V are apparent
    # (CL/F, V/F) and every value below is conditional on F = 1.
    # ------------------------------------------------------------------------
    lka <- log(14.89);  label("Absorption rate constant from the gut to the central compartment (1/h)")
    # Table 3, 'K a (h-1)' mean = 14.89 (median 9.88, SD 13.31, CV 89.36%)
    lcl <- log(608.02); label("Apparent clearance from the central compartment, CL/F (L/h)")
    # Table 3, 'SCL/F (liter/h)' mean = 608.02 (median 677.59, SD 320.10, CV 52.65%)
    lvc <- log(473.97); label("Apparent central volume of distribution, V/F (L)")
    # Table 3, 'V/F (liter)' mean = 473.97 (median 503.90, SD 244.80, CV 51.65%)

    lk12 <- log(25.90); label("Transfer rate constant central -> peripheral1 (1/h)")
    # Table 3, 'K cp (h-1)' mean = 25.90 (median 22.28, SD 18.89, CV 72.92%)
    lk21 <- log(33.30); label("Transfer rate constant peripheral1 -> central (1/h)")
    # Table 3, 'K pc (h-1)' mean = 33.30 (median 25.71, SD 26.08, CV 78.32%)

    ltlag <- log(1.94); label("Absorption lag time (h)")
    # Table 3, 'T lag (h)' mean = 1.94 (median 1.67, SD 1.13, CV 58.34%)

    lfdepot <- fixed(log(1)); label("Oral bioavailability of the depot (unitless)")
    # Methods (Pharmacokinetic population analyses): 'Oral bioavailability was not estimated.'
    # Table 3 footnote a: 'F, oral bioavailability of PZQ, which was not estimated in this study.'
    # Anchored at 1 so that CL and V carry the /F, exactly as the paper reports them.

    # ------------------------------------------------------------------------
    # Inter-individual variability. NPAG estimates a discrete non-parametric
    # distribution rather than a parametric omega. Table 3 reports a mean, a
    # median, an SD and a CV% per parameter, with CV% = SD / mean on the LINEAR
    # scale; the CV% column is carried here as a LOG-NORMAL approximation using
    # omega^2 = log(CV^2 + 1). This is an approximation imposed on a
    # non-parametric distribution -- see vignette Assumptions and deviations.
    # Off-diagonal covariances are not reported and are therefore absent.
    # ------------------------------------------------------------------------
    # Table 3 Ka CV% = 89.36 -> log(0.8936^2 + 1)
    etalka   ~ 0.586965
    # Table 3 SCL/F CV% = 52.65 -> log(0.5265^2 + 1)
    etalcl   ~ 0.244672
    # Table 3 V/F CV% = 51.65 -> log(0.5165^2 + 1)
    etalvc   ~ 0.236472
    # Table 3 Kcp CV% = 72.92 -> log(0.7292^2 + 1)
    etalk12  ~ 0.426400
    # Table 3 Kpc CV% = 78.32 -> log(0.7832^2 + 1)
    etalk21  ~ 0.478345
    # Table 3 Tlag CV% = 58.34 -> log(0.5834^2 + 1)
    etaltlag ~ 0.292935

    # ------------------------------------------------------------------------
    # Residual unexplained variability is NOT reported. Methods state only that
    # the fit was assessed 'by mean weighted error ... and by mean weighted
    # squared error'; neither the Pmetrics assay-error polynomial coefficients
    # (C0-C3) nor a lambda/gamma term appears anywhere in the paper, and the sole
    # supplemental item is Figure S1 (a side-effect frequency plot). Carried as
    # fixed(0) rather than invented -- see vignette Errata.
    # ------------------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")
    addSd  <- fixed(0); label("Additive residual SD (mg/L; 0 -- not reported in the source)")
  })

  model({
    # 1. Individual parameters. No covariate was retained (see
    #    covariatesDataExcluded): the base model is the final model.
    ka   <- exp(lka + etalka)
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    tlag <- exp(ltlag + etaltlag)

    kcp <- exp(lk12 + etalk12)
    kpc <- exp(lk21 + etalk21)

    # 2. Intercompartmental clearance and peripheral volume, derived from the
    #    primary rate constants. The ODEs below MUST be driven by q / vp rather
    #    than by kcp / kpc directly: rxSolve() defaults to useLinCmt = TRUE and,
    #    when the peripheral transfer is written straight from stored
    #    micro-constants, that rewrite can silently drop peripheral1 and solve a
    #    one-compartment model. Routing through q and vp keeps the closed-form
    #    and ODE solvers in agreement (the vignette asserts this identity).
    q  <- kcp * vc
    vp <- q / kpc

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. Two-compartment disposition with a first-order absorption (gut)
    #    compartment. Methods: 'The three compartments were an absorptive
    #    compartment (i.e., gut), a central compartment (bloodstream), and a
    #    peripheral compartment.'
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # 5. Absorption lag and the unit bioavailability anchor.
    alag(depot) <- tlag
    f(depot)    <- exp(lfdepot)

    # 6. Observation. The assay measured total plasma praziquantel as the sum of
    #    the separately quantified R and S enantiomers (Results; Methods).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
