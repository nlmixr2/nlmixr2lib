Kapralos_2021_octreotide <- function() {
  description <- paste(
    "One-compartment population PK model for octreotide after a single 30 mg",
    "intramuscular injection of the long-acting repeatable (LAR, Sandostatin LAR",
    "Depot) PLGA-microsphere formulation in healthy adult male volunteers. Release",
    "from the depot is a weighted sum of four processes: an initial burst dosed",
    "directly into a first-order absorption compartment, and three parallel",
    "delayed releases each described by the Savic (2007) closed-form transit",
    "density (Stirling approximation of n!) with its own mean transit time and",
    "transit count. The four fractions follow a multivariate logistic-normal",
    "distribution so they stay in (0, 1) and sum to one, the first mean transit",
    "time is logit-constrained below 300 h, and the second and third mean transit",
    "times are built as positive increments on the preceding one so the three",
    "delayed releases stay in sequential order. A binary subpopulation indicator",
    "from a pre-fit shape-respecting k-means clustering of the individual profiles",
    "(13% of subjects, early extended release) shifts the apparent clearance and",
    "the logits of the second and third release fractions. Single-dose model: the",
    "delayed-release inputs are driven by the most recent dose only."
  )
  reference <- paste(
    "Kapralos I, Dokoumetzidis A. Population Pharmacokinetic Modelling of the",
    "Complex Release Kinetics of Octreotide LAR: Defining Sub-Populations by",
    "Cluster Analysis. Pharmaceutics. 2021;13(10):1578.",
    "doi:10.3390/pharmaceutics13101578. Final-model estimates from Table 2;",
    "structural and variability equations from Equations 1-7 and Figure 1.",
    "The transit-compartment closed form follows Savic RM, Jonker DM, Kerbusch T,",
    "Karlsson MO. J Pharmacokinet Pharmacodyn. 2007;34(5):711-726."
  )
  vignette <- "Kapralos_2021_octreotide"
  units <- list(time = "h", dosing = "mg", concentration = "pg/mL")
  # The single IM dose is given to `depot` (the absorption compartment of
  # Figure 1). f(depot) takes only the burst fraction; the three delayed
  # fractions reach `depot` through the closed-form transit inputs, which read
  # the dose amount through podo(depot).
  dosing <- c("depot")

  compartmentData <- list(
    depot = list(analyte = "octreotide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "octreotide", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    MIX_EARLY_REL = list(
      description = "Subpopulation indicator from the pre-fit shape-respecting k-means clustering of individual PK profiles: 1 = early extended-release profile (paper 'cluster 2'), 0 = typical multi-phase release profile (paper 'cluster 1')",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (cluster 1, typical multi-phase release; 103 of 118 subjects)",
      notes = paste(
        "Kapralos 2021 Equation 7 codes the cluster as Parameter_pop = theta1 +",
        "(cluster - 1) * theta2 with cluster = 1 or 2, so MIX_EARLY_REL = cluster - 1.",
        "The cluster was assigned before the NONMEM fit by kmlShape clustering of",
        "each subject's concentration profile normalised by its mean concentration",
        "(Section 2.2.1, lambda = 0.001 per Section 3.1); it is not a measured",
        "clinical covariate and not a NONMEM $MIX class. 15 of 118 subjects",
        "(12.7%, Table 1) were in cluster 2; for population simulation draw",
        "MIX_EARLY_REL ~ Bernoulli(15/118), for a typical subject set it to 0."
      ),
      source_name = "cluster"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on the disposition parameters and not retained (Section 2.2.4 and Discussion: 'Covariates of size and age failed to explain the population variability'). Median 75 kg (Q1-Q3 66-86, Table 1).",
      source_name = "body weight"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on the disposition parameters and not retained (Discussion). Median 28 years (Q1-Q3 23-37, Table 1).",
      source_name = "age"
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Available demographic (Section 2.1); not retained in the final model. Median 175 cm (Q1-Q3 170-178, Table 1).",
      source_name = "height"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Available demographic (Section 2.1); not retained in the final model. Median 24.75 kg/m^2 (Q1-Q3 22.4-27.7, Table 1).",
      source_name = "BMI"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 118,
    n_studies = 1,
    n_observations = 3936,
    age_range = "Q1-Q3 23-37 years",
    age_median = "28 years",
    weight_range = "Q1-Q3 66-86 kg",
    weight_median = "75 kg",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state = "healthy volunteers",
    dose_range = "30 mg octreotide LAR (Sandostatin LAR Depot) single deep intramuscular injection, fasting",
    regions = "Jordan (Triumpharma CRO, Amman)",
    notes = paste(
      "Reference arm of a phase 1 single-dose bioequivalence study (Section 2.1).",
      "Dense sampling: pre-dose plus 36 samples from 0.5 to 2088 h. The cohort",
      "consisted solely of Caucasian males (Section 3). Demographics are medians",
      "(Q1-Q3) from Table 1. Cluster 1: 103 subjects; cluster 2: 15 subjects",
      "(Table 1). Serum octreotide by LC-MS/MS, calibration range",
      "8.835-4010.010 pg/mL (Section 2.1)."
    )
  )

  ini({
    # Disposition (apparent values; F is not identifiable, Section 2.2.4)
    lka <- log(0.27); label("First-order absorption rate constant from the absorption compartment ka (1/h)") # Table 2 ka = 0.27 (RSE 2.2%)
    lcl <- log(32.7); label("Apparent clearance CL/F, cluster 1 (L/h)") # Table 2 CL = 32.7 (RSE 5.8%)
    lvc <- log(15.3); label("Apparent volume of distribution V/F (L)") # Table 2 V = 15.3 (RSE 7.7%)
    e_mix_early_rel_cl <- -8.61; label("Additive shift in CL/F for the early-release cluster (L/h)") # Table 2 CL 'Cluster effect' = -8.61 (RSE 34%); Equation 7

    # Release fractions: additive log-ratio (multivariate logistic-normal)
    # coordinates u1..u3 of Equation 5, with the fourth (slowest) delayed
    # release as the reference category.
    logitfburst <- -5.18; label("Log-ratio u1 of the burst fraction to the third-delayed-release fraction (unitless)") # Table 2 YF1 = -5.18 (RSE 1.8%)
    logitfdel1 <- -3.36; label("Log-ratio u2 of the first-delayed-release fraction to the third-delayed-release fraction, cluster 1 (unitless)") # Table 2 YF2 = -3.36 (RSE 7.9%)
    logitfdel2 <- -1.54; label("Log-ratio u3 of the second-delayed-release fraction to the third-delayed-release fraction, cluster 1 (unitless)") # Table 2 YF3 = -1.54 (RSE 2.8%)
    e_mix_early_rel_logitfdel1 <- 3.06; label("Additive shift in u2 for the early-release cluster (unitless)") # Table 2 YF2 'Cluster effect' = 3.06 (RSE 33%); Equation 7
    e_mix_early_rel_logitfdel2 <- -0.523; label("Additive shift in u3 for the early-release cluster (unitless)") # Table 2 YF3 'Cluster effect' = -0.523 (RSE 26.8%); Equation 7

    # Mean transit times (Equations 3 and 4)
    mtt1max <- fixed(300); label("Upper bound of the first delayed-release mean transit time MTT1 (h)") # Section 2.2.4 and Equation 4: 'MTT1 was constrained to 300 h'
    logitmtt1 <- -0.421; label("Logit of MTT1 / 300 h, YMTT1 (unitless)") # Table 2 YMTT1 = -0.421 (RSE 21.8%); Equation 4
    ldmtt2 <- log(181); label("Increment of MTT2 over MTT1, theta_2 of Equation 3 (h)") # Table 2 'MTT2' = 181 (RSE 3.3%); Equation 3 MTT_j = MTT_(j-1) + theta_j * exp(eta_j)
    ldmtt3 <- log(506); label("Increment of MTT3 over MTT2, theta_3 of Equation 3 (h)") # Table 2 'MTT3' = 506 (RSE 3.8%); Equation 3

    # Transit counts of the three delayed releases (Equation 2)
    lnn1 <- log(3.42); label("Number of transit compartments N1 of the first delayed release (unitless)") # Table 2 N1 = 3.42 (RSE 15%)
    lnn2 <- log(17.9); label("Number of transit compartments N2 of the second delayed release (unitless)") # Table 2 N2 = 17.9 (RSE 6%)
    lnn3 <- log(5.08); label("Number of transit compartments N3 of the third delayed release (unitless)") # Table 2 N3 = 5.08 (RSE 5%)

    # IIV. Table 2 reports IIV as CV% = sqrt(omega^2) * 100 (Section 2.2.4),
    # i.e. the SD of eta times 100, so omega^2 = (CV/100)^2. No IIV on ka
    # (Section 3.2). Table 2 prints only the diagonal; the fraction/MTT
    # covariances retained in the fit (Section 3.2) are not reported.
    etalvc ~ 0.155236 # Table 2 IIV V = 39.4% -> 0.394^2 [shrinkage 16.3%]
    etalcl ~ 0.079524 # Table 2 IIV CL = 28.2% -> 0.282^2 [shrinkage 1%]
    etalogitfburst ~ 0.083521 # Table 2 IIV YF1 = 28.9% -> 0.289^2 [shrinkage 3.4%]
    etalogitfdel1 ~ 1.658944 # Table 2 IIV YF2 = 128.8% -> 1.288^2 [shrinkage 4%]
    etalogitfdel2 ~ 0.042025 # Table 2 IIV YF3 = 20.5% -> 0.205^2 [shrinkage 30.3%]
    etalogitmtt1 ~ 0.361201 # Table 2 IIV YMTT1 = 60.1% -> 0.601^2 [shrinkage 12%]
    etaldmtt2 ~ 0.029929 # Table 2 IIV MTT2 = 17.3% -> 0.173^2 [shrinkage 17.6%]
    etaldmtt3 ~ 0.040804 # Table 2 IIV MTT3 = 20.2% -> 0.202^2 [shrinkage 1.7%]
    etalnn1 ~ 0.506944 # Table 2 IIV N1 = 71.2% -> 0.712^2 [shrinkage 22%]
    etalnn2 ~ 0.068644 # Table 2 IIV N2 = 26.2% -> 0.262^2 [shrinkage 31.2%]
    etalnn3 ~ 0.098596 # Table 2 IIV N3 = 31.4% -> 0.314^2 [shrinkage 9.3%]

    # Residual error (combined proportional + additive)
    propSd <- 0.143; label("Proportional residual error (fraction)") # Table 2 'Proportional Residual Error' = 0.143 (RSE 1.3%)
    addSd <- 28.4; label("Additive residual error (pg/mL)") # Table 2 'Additive Residual Error' = 28.4 (RSE 3.7%)
  })
  model({
    # 1. Disposition. Equation 7 puts the cluster effect on the linear scale of
    #    the typical value; Equation 6 then applies the log-normal eta.
    ka <- exp(lka)
    cl <- (exp(lcl) + e_mix_early_rel_cl * MIX_EARLY_REL) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    # 2. Release fractions, Equation 5 (multivariate logistic-normal,
    #    Tsamandouras 2015). fburst = burst, fdel1..fdel3 = first..third delayed release
    #    (Figure 1); fdel3 is the reference category, so all four lie in (0, 1)
    #    and sum to one.
    u1 <- logitfburst + etalogitfburst
    u2 <- logitfdel1 + e_mix_early_rel_logitfdel1 * MIX_EARLY_REL + etalogitfdel1
    u3 <- logitfdel2 + e_mix_early_rel_logitfdel2 * MIX_EARLY_REL + etalogitfdel2
    fden <- exp(u1) + exp(u2) + exp(u3) + 1
    fburst <- exp(u1) / fden
    fdel1 <- exp(u2) / fden
    fdel2 <- exp(u3) / fden
    fdel3 <- 1 / fden

    # 3. Mean transit times. Equation 4: MTT1 = 300 * e^y / (e^y + 1), y ~ N.
    #    Equation 3: each later MTT is the previous one plus a positive
    #    log-normal increment, which keeps the three releases in order.
    ymtt1 <- logitmtt1 + etalogitmtt1
    mtt1 <- mtt1max * exp(ymtt1) / (exp(ymtt1) + 1)
    dmtt2 <- exp(ldmtt2 + etaldmtt2)
    dmtt3 <- exp(ldmtt3 + etaldmtt3)
    mtt2 <- mtt1 + dmtt2
    mtt3 <- mtt2 + dmtt3

    # 4. Transit counts and rate constants. Figure 1: MTT_i = (n_i + 1) / ktr_i.
    nn1 <- exp(lnn1 + etalnn1)
    nn2 <- exp(lnn2 + etalnn2)
    nn3 <- exp(lnn3 + etalnn3)
    ktr1 <- (nn1 + 1) / mtt1
    ktr2 <- (nn2 + 1) / mtt2
    ktr3 <- (nn3 + 1) / mtt3

    # 5. Equation 2, the Savic closed-form transit density with n! replaced by
    #    Stirling's approximation sqrt(2 * pi) * n^(n + 0.5) * exp(-n), kept in
    #    the as-published form (the density therefore integrates to
    #    Gamma(n + 1) / Stirling(n), about 1 + 1/(12 n), rather than exactly 1).
    #    Evaluated on the log scale for numerical stability. The time argument
    #    is time after the most recent dose; the paper is single-dose.
    tdose <- tad(depot)
    doseamt <- podo(depot)
    lstir1 <- 0.5 * log(2 * pi) + (nn1 + 0.5) * log(nn1) - nn1
    lstir2 <- 0.5 * log(2 * pi) + (nn2 + 0.5) * log(nn2) - nn2
    lstir3 <- 0.5 * log(2 * pi) + (nn3 + 0.5) * log(nn3) - nn3
    rin1 <- exp(log(ktr1) + nn1 * log(ktr1 * tdose) - ktr1 * tdose - lstir1)
    rin2 <- exp(log(ktr2) + nn2 * log(ktr2 * tdose) - ktr2 * tdose - lstir2)
    rin3 <- exp(log(ktr3) + nn3 * log(ktr3 * tdose) - ktr3 * tdose - lstir3)

    # 6. ODEs. Equation 1: the absorption compartment receives the weighted
    #    transit inputs and is drained at ka; the burst fraction is the dose
    #    bolus itself (Figure 1, F1 straight into the absorption compartment).
    d/dt(depot) <- doseamt * (fdel1 * rin1 + fdel2 * rin2 + fdel3 * rin3) - ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- fburst

    # 7. Observation: mg/L -> pg/mL (x 1e6).
    Cc <- central / vc * 1e6
    Cc ~ add(addSd) + prop(propSd)
  })
}
