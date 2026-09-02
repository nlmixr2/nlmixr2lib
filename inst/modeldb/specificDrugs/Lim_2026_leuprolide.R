Lim_2026_leuprolide <- function() {
  description <- "One-compartment population PK model for the leuprolide acetate 3-month depot formulation in children with central precocious puberty, with parallel immediate and transit-delayed first-order absorption"
  reference <- paste(
    "Lim CN, Al Yacoub ON, Mostafa NM, Salem AH.",
    "Fixed Dosing of Leuprolide Acetate, a GnRH Agonist, in Children with",
    "Central Precocious Puberty: A Population Pharmacokinetic Justification.",
    "Pediatric Drugs. 2026;28:295-305. doi:10.1007/s40272-025-00733-2",
    sep = " "
  )
  vignette <- "Lim_2026_leuprolide"

  # Amounts are carried in micrograms so that `central / vc` (vc in L) is
  # directly in ng/mL, the unit the paper's assay and LLOQ (0.025 ng/mL) use.
  # An 11.25 mg dose is therefore given as amt = 11250 and a 30 mg dose as
  # amt = 30000.
  units <- list(time = "day", dosing = "ug", concentration = "ng/mL")

  # The final model carries no covariates: "In the covariates assessment, none
  # of the covariates, including body size indices and age, was found to be
  # significant" (Results, paragraph following the structural equations), so
  # `covariateData` is empty and everything the paper screened is documented in
  # `covariatesDataExcluded` below.
  covariateData <- list()

  # Screened in the forward-inclusion / backward-elimination covariate search
  # (Methods 2.5 "Development of the Covariate Model") and in the exploratory
  # empirical-Bayes-estimate-versus-covariate plots (Fig. 2), but NOT retained
  # in the final model. Documented here so the paper's covariate screen is not
  # lost; none of these names is referenced in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as a power model scaled by the population median",
        "(Methods 2.5). Not significant at the 0.01 forward-inclusion level;",
        "Fig. 2 shows the empirical-Bayes CL and V estimates against body",
        "weight with no trend. The absence of a weight effect is the paper's",
        "central conclusion (fixed rather than weight-based dosing)."
      ),
      source_name = "Body weight"
    ),
    AGE = list(
      description = "Age",
      units = "year",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as a power model scaled by the population median",
        "(Methods 2.5); not retained. Fig. 2 shows no trend in the",
        "empirical-Bayes CL or V estimates against age over 1-10 years."
      ),
      source_name = "Age"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = paste(
        "Screened (Methods 2.5); not retained. The cohort was 91.7% female",
        "(Table 1), so the male stratum carried very little information."
      ),
      source_name = "Sex"
    ),
    BSA = list(
      description = "Body surface area, Mosteller formula",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods 2.5, Mosteller formula); not retained.",
      source_name = "Body surface area"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods 2.5); not retained.",
      source_name = "Body mass index"
    ),
    CRCL = list(
      description = "Creatinine clearance, Cockcroft-Gault, raw and NOT BSA-normalized",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened (Methods 2.5, Cockcroft-Gault); not retained. Table 1",
        "reports mean (SD) [range] of 117 (31.6) [39.5-193] mL/min in the",
        "11.25 mg arm and 107 (29.7) [64.1-162] mL/min in the 30 mg arm.",
        "Raw mL/min, not normalized to 1.73 m^2."
      ),
      source_name = "Creatinine clearance"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "not reported",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as one of the liver-function markers (Methods 2.5);",
        "not retained. The paper does not tabulate the values or state the",
        "unit convention, so no unit is asserted here."
      ),
      source_name = "bilirubin"
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "not reported",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as one of the liver-function markers (Methods 2.5);",
        "not retained. Values and unit convention are not reported."
      ),
      source_name = "blood urea nitrogen"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "not reported",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as one of the liver-function markers (Methods 2.5);",
        "not retained. Values and unit convention are not reported."
      ),
      source_name = "aspartate aminotransferase"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "not reported",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as one of the liver-function markers (Methods 2.5);",
        "not retained. Values and unit convention are not reported."
      ),
      source_name = "alanine transaminase"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "leuprolide", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    depot2 = list(
      analyte = "leuprolide", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "leuprolide", units = "ug",
      specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 48,
    n_studies = 1,
    n_observations = 293,
    age_range = "1-10 years",
    age_median = "mean 7.88 years (11.25 mg arm) and 7.75 years (30 mg arm); medians not reported",
    weight_range = "14.1-62.6 kg",
    weight_median = "mean 38.2 kg (11.25 mg arm) and 36.4 kg (30 mg arm); medians not reported",
    sex_female_pct = 91.7,
    race_ethnicity = c(White = 56.3, `Black/African American` = 20.8, Other = 22.9),
    disease_state = "central precocious puberty",
    dose_range = "11.25 mg or 30 mg leuprolide acetate 3-month depot intramuscularly, two injections 3 months apart (day 1 and day 84)",
    regions = "USA, including Puerto Rico (22 sites)",
    renal_function = "creatinine clearance (Cockcroft-Gault) 39.5-193 mL/min",
    notes = paste(
      "Randomized, open-label phase III study NCT00635817 (Methods 2.1).",
      "80 children were enrolled; the PK analysis set is the sampled subset of",
      "48 participants (24 per dose arm) contributing 293 plasma leuprolide",
      "concentrations (Methods 2.2, Results paragraph 1). Baseline",
      "characteristics are Table 1; the two arms are summarized separately",
      "there and the values above pool them where the pooled figure is",
      "unambiguous. Samples were drawn pre-dose and at 0.5 and 1 hour and at",
      "2, 4, 8, 12 and 24 weeks after the first injection. LLOQ 0.025 ng/mL;",
      "about 3.5% of observations were below it and were imputed at LLOQ/2",
      "(0.0125 ng/mL) by the M5 method."
    )
  )

  ini({
    # ---------------- Structural parameters (Table 2, "Original dataset") ----
    # Apparent (CL/F, Vd/F) because only intramuscular depot data were
    # available; F is not separately identifiable. Time unit is the day.
    lcl <- log(181)
    label("Apparent clearance CL/F (L/day)")
    # Table 2: CL/F = 181 L/day (%RSE 34.9; bootstrap median 184, 5th-95th 152-215)

    lvc <- log(7.11)
    label("Apparent central volume of distribution Vd/F (L)")
    # Table 2: Vd/F = 7.11 L (%RSE 35.2; bootstrap median 6.97, 5th-95th 5.37-8.73)

    lka <- log(0.441)
    label("Immediate first-order absorption rate constant Ka1 out of depot (1/day)")
    # Table 2: Ka1 = 0.441 1/day (%RSE 7.81; bootstrap median 0.438, 5th-95th 0.414-0.464)

    lka2 <- log(0.00879)
    label("Delayed first-order absorption rate constant Ka2 out of depot2 (1/day)")
    # Table 2: Ka2 = 0.00879 1/day (%RSE 29.6; bootstrap median 0.00892, 5th-95th 0.00424-0.0196)

    logitffo <- log(0.920 / (1 - 0.920))
    label("Fraction of dose absorbed via the immediate first-order process FRAC (unitless, logit scale)")
    # Table 2: FRAC = 0.920 (%RSE 3.77; bootstrap median 0.919, 5th-95th 0.866-0.940).
    # Fig. 1 names this F1; the delayed fraction is F2 = 1 - FRAC = 0.080.
    # The paper reports FRAC on the natural scale with no IIV; the logit
    # encoding reproduces 0.920 exactly and cannot leak above 1 downstream.

    lmtt <- log(34.1)
    label("Mean transit time MTT of the delayed-absorption transit chain (day)")
    # Table 2: MTT = 34.1 day (%RSE 8.89; bootstrap median 33.9, 5th-95th 27.8-43.7)

    lntr <- fixed(log(3))
    label("Number of transit compartments N in the delayed-absorption chain (unitless)")
    # Table 2: N = 3 (FIXED). Results: "The number of transit compartments was
    # fixed due to insufficient data to estimate the parameter with good
    # precision. The number of compartments was fixed at three based on a
    # successful model run that adequately described the absorption phase".

    # ---------------- Inter-individual variability (Table 2, lower block) ----
    # The Table 2 lower block is headed "Estimate (%CV) [shrinkage, %]", but the
    # parenthetical is the relative standard ERROR of the estimate, not a %CV of
    # the random effect: back-computing sqrt(exp(omega^2) - 1) reproduces the
    # printed number only for CL (42.1%) and misses every other row by a factor
    # of 2 or more, while (upper - lower) / 2 / estimate / 1.645 from the
    # bootstrap 5th-95th percentiles in the same row reproduces all six within
    # about 3 percentage points. The "Estimate" column therefore holds the
    # OMEGA VARIANCE on the log scale, which is also what the covariance row
    # requires to be self-consistent: 0.127 / sqrt(0.163 * 0.323) = 0.554, a
    # valid correlation. See the vignette Errata for the full argument.
    etalcl + etalvc ~ c(
      0.163,
      0.127, 0.323
    )
    # Table 2: IIV CL/F = 0.163 (%RSE 42.1) [shrinkage 15.8%];
    #          Covariance of CL/F and Vd/F = 0.127 (%RSE 55.1);
    #          IIV Vd/F = 0.323 (%RSE 49.1) [shrinkage 14.4%]

    etalka ~ 0.0217
    # Table 2: IIV Ka1 = 0.0217 (%RSE 39.0) [shrinkage 21.2%]

    etalka2 ~ 0.432
    # Table 2: IIV Ka2 = 0.432 (%RSE 43.7) [shrinkage 19.3%]

    etalmtt ~ 0.0763
    # Table 2: IIV MTT = 0.0763 (%RSE 56.6) [shrinkage 34.4%]

    # ---------------- Residual unexplained variability ----------------------
    # Results equations: Cobs_ij = Cpred_ij * (1 + eps1_ij), eps1_ij ~ N(0,
    # sigma^2_1,1) -- a pure proportional error on the linear scale. Table 2
    # prints 0.132 in the same block as the OMEGA variances, so it is the
    # NONMEM $SIGMA VARIANCE; nlmixr2's prop() takes the standard deviation,
    # hence sqrt(0.132) = 0.3633 (36.3% CV).
    propSd <- sqrt(0.132)
    label("Proportional residual error on plasma leuprolide concentration (fraction)")
    # Table 2: Proportional residual error = 0.132 (%RSE 14.8; bootstrap median
    # 0.132, 5th-95th 0.102-0.166)
  })

  model({
    # ---------------- Individual parameters ---------------------------------
    # Results: Ka1_i = theta1 * exp(eta1_i), Ka2_i = theta2 * exp(eta2_i),
    # CL_i = theta3 * exp(eta3_i), V_i = theta4 * exp(eta4_i), N_i = theta5,
    # MTT_i = theta6 * exp(eta5_i). Inter-individual variability was
    # "assumed to be lognormally distributed and modeled using an exponential
    # error structure" (Methods 2.3).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalka)
    ka2 <- exp(lka2 + etalka2)
    mtt <- exp(lmtt + etalmtt)
    ntr <- exp(lntr)

    # FRAC on the natural scale; F2 = 1 - FRAC is the delayed fraction.
    ffo <- 1 / (1 + exp(-logitffo))

    # ---------------- Micro-constants ---------------------------------------
    kel <- cl / vc # Results: K30_i = CL_i / V_i
    ktr <- (ntr + 1) / mtt # Results: Ktr_i = (N_i + 1) / MTT_i

    # ---------------- Delayed-absorption (transit) input --------------------
    # First term of the paper's differential equation (2), reproduced exactly:
    #
    #   Dose * F2 * Ktr * (Ktr * tad)^N * exp(-Ktr * tad)
    #                     / ( sqrt(2 * pi) * N^(N + 0.5) * exp(-N) )
    #
    # The denominator is Stirling's approximation of N!; the paper states the
    # term "was transformed to a logarithmic form to prevent numerical
    # difficulties for a large N", which is the Savic (2007) log-domain
    # transit implementation. It is kept here rather than substituting
    # rxode2's transit(), whose lgamma(N + 1) uses the EXACT 3! = 6: with
    # N = 3 the Stirling denominator is 5.8355, so the published
    # implementation delivers 6 / 5.8355 = 1.0282 times F2 * Dose over the
    # whole profile. Reproducing the authors' 2.8% excess is what makes the
    # published parameter estimates self-consistent.
    #
    # tad(depot) is the paper's t_ad (time after dose) and podo(depot) is the
    # administered dose amount BEFORE bioavailability, matching NONMEM's PODO.
    tad_dose <- tad(depot)
    lfact_ntr <- 0.5 * log(2 * pi) + (ntr + 0.5) * log(ntr) - ntr
    transit_in <- podo(depot) * (1 - ffo) * ktr *
      exp(ntr * log(ktr * tad_dose) - ktr * tad_dose - lfact_ntr)

    # ---------------- ODE system --------------------------------------------
    # Differential equations (1), (2) and (3) of the Results section.
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- transit_in - ka2 * depot2
    d/dt(central) <- ka * depot + ka2 * depot2 - kel * central

    # Fig. 1: F1 = fraction of dose absorbed via the immediate first-order
    # process = FRAC. The delayed fraction F2 enters through transit_in above.
    f(depot) <- ffo

    # ---------------- Observation and residual error ------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
