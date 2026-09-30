Chan_2020_fenebrutinib_das28 <- function() {
  description <- paste(
    "Longitudinal exposure-response model for DAS28 (CRP) in rheumatoid",
    "arthritis patients treated with fenebrutinib (GDC-0853) or placebo in",
    "the phase 2 ANDES trial (Chan 2020). DAS28 is a log-normally",
    "distributed subject baseline plus a hyperbolic (Emax in study day)",
    "time course whose maximum is the placebo maximum plus an Emax function",
    "of the individual steady-state daily fenebrutinib AUC; the onset time",
    "differs between fenebrutinib-treated and placebo patients, and the",
    "total maximum carries a log-normal random effect. Additive residual",
    "error. The AUC is a covariate column produced by",
    "Chan_2020_fenebrutinib."
  )
  reference <- paste(
    "Chan P, Yu J, Chinn L, Prohn M, Huisman J, Matzuka B, Hanley W,",
    "Tuckwell K, Quartino A. Population Pharmacokinetics, Efficacy",
    "Exposure-response Analysis, and Model-based Meta-analysis of",
    "Fenebrutinib in Subjects with Rheumatoid Arthritis. Pharm Res.",
    "2020;37(2):25. doi:10.1007/s11095-019-2752-y"
  )
  vignette <- "Chan_2020_fenebrutinib"
  units <- list(time = "day", dosing = "mg", concentration = "das28 (DAS28-CRP score, unitless)")

  covariateData <- list(
    AUC_FENEBRUTINIB = list(
      description = "Individual steady-state daily (0-24 h) AUC of fenebrutinib",
      units = "h*ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Empirical Bayes (post hoc) estimate from the Chan 2020 popPK model",
        "(Chan_2020_fenebrutinib) at the nominal dose. 0 for placebo. It",
        "also switches the onset time: Model S3 uses the drug onset time",
        "when DOSE > 0 and the placebo onset time when DOSE = 0, and in this",
        "dataset (placebo and fenebrutinib arms only) DOSE > 0 is",
        "equivalent to AUC_FENEBRUTINIB > 0."
      ),
      source_name = "AUC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 467L,
    n_studies = 1L,
    age_range = "19-75 years (median 52)",
    weight_range = "38-153 kg (median 71)",
    sex_female_pct = 80.5,
    disease_state = "Adults with moderate to severe active seropositive rheumatoid arthritis on stable methotrexate, with inadequate response to methotrexate (cohort 1) or to anti-TNF therapy (cohort 2)",
    dose_range = "Placebo or fenebrutinib 50 mg QD, 150 mg QD or 200 mg BID tablets for 12 weeks (adalimumab arm not included)",
    regions = "Eastern Europe 62%, Latin America 33%, US 5%",
    n_observations = "2676 DAS28 (CRP) observations on days 0, 7, 14, 28, 56 and 84 (LOCF imputation, withdrawals excluded)",
    notes = "Chan 2020 Methods 'DAS28 (CRP)', Table S3 (demographics) and Table S5 (estimates); FOCE-I in NONMEM. No covariate analysis was performed."
  )

  ini({
    # Model S3 $PRED. THETA(3), THETA(4) and THETA(5) enter through EXP(),
    # so Table S5 prints them back-transformed (no RSE, CI only).
    lrbase_das28 <- log(5.46); label("Log typical baseline DAS28 (CRP) score")  # Table S5, theta1 'Baseline DAS28 score' = 5.46 (RSE 0.9%); log-normal via C = TVC*EXP(ETA(1))
    emax_pbo_das28 <- -1.36; label("Maximum placebo effect over time on DAS28 (score units)")  # Table S5, theta2 'Maximum placebo effect over time (DAS28)' = -1.36 (RSE 7.9%)
    lt50_pbo_das28 <- log(36.7); label("Log time of 50% of the maximum effect in placebo patients (day)")  # Table S5, theta3 'Time of 50% placebo effect (d)' = 36.7 (95% CI 25.7-52.4)
    lt50_das28 <- log(47); label("Log time of 50% of the maximum effect in fenebrutinib-treated patients (day)")  # Table S5, theta4 'Time of 50% drug effect (d)' = 47 (95% CI 39.2-56.3)
    lec50_das28 <- log(293); label("Log daily steady-state AUC giving 50% of the maximum drug effect (h*ng/mL)")  # Table S5, theta5 'Exposure at which 50% drug effect (AUC ng.hr/mL)' = 293 (95% CI 14.5-5920)
    emax_das28 <- -0.964; label("Maximum fenebrutinib effect over time on DAS28 (score units)")  # Table S5, theta6 'Maximum drug effect over time (DAS28)' = -0.964 (RSE 19.4%)

    etalrbase_das28 ~ 0.0206 # Table S5, omega1.1 'omega2 Baseline' = 0.0206 (RSE 8.6%; shrinkage 8.63%)
    # Model S3 puts ETA(2) on the TOTAL maximum, TMAX = (THETA(2) +
    # THETA(6)*AUC/(EAUC50+AUC)) * EXP(ETA(2)), so it scales the placebo part
    # as well as the drug part.
    etaemax_das28 ~ 0.283 # Table S5, omega2.2 'omega2 Max drug effect' = 0.283 (RSE 13.6%; shrinkage 25.06%)

    # Model S3 'Y = IPRED + EPS(1)'. Table S5 prints sigma = 0.3 with RSE
    # 2.2%, below the asymptotic floor sqrt(2/2676) = 2.7% for an estimated
    # VARIANCE but above the floor for an SD, so 0.3 is read as the SD.
    addSd_das28 <- 0.3; label("Additive residual SD on DAS28 (score units)")  # Table S5, sigma 'Additive residual error in patients' = 0.3 (RSE 2.2%)
  })

  model({
    trt <- (AUC_FENEBRUTINIB > 0)
    base <- exp(lrbase_das28 + etalrbase_das28)
    emaxTot <- (emax_pbo_das28 + emax_das28 * AUC_FENEBRUTINIB /
      (exp(lec50_das28) + AUC_FENEBRUTINIB)) * exp(etaemax_das28)
    # Model S3 fixes THILL = 1; `time` is study day.
    t50 <- trt * exp(lt50_das28) + (1 - trt) * exp(lt50_pbo_das28)
    das28 <- base + emaxTot * time / (time + t50)
    das28 ~ add(addSd_das28)
  })
}
