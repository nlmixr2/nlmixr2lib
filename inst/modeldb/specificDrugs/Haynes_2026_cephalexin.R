Haynes_2026_cephalexin <- function() {
  description <- "One-compartment oral population PK model for cephalexin pooled from infancy to adulthood (Haynes 2026; 71 subjects aged 9 days to 29 years from four studies in the United States, Canada and New Zealand). First-order absorption with an absorption lag time and linear elimination. Apparent volume (Vd/F) and apparent clearance (CL/F) scale with fat-free mass by estimated power exponents (reference 56 kg FFM); CL/F additionally matures with postmenstrual age through a sigmoid Hill function (TM50 42.96 weeks, Hill coefficient 5.06), and the absorption rate constant is reduced by a fixed scalar exp(-0.81) in infants with postnatal age <= 21 days. Log-normal inter-individual variability on lag time, ka and CL/F (none on Vd/F); combined additive plus proportional residual error."
  reference <- "Haynes AS, Gonzalez D, Messacar K, Autmizguine J, Gwee A, Everts R, Anderson P. Cephalexin Model-Informed Dosing Recommendations for Severe Infections: A Pooled Pharmacokinetic Analysis From Infancy to Adulthood. Open Forum Infect Dis. 2026;13(10):ofag611. doi:10.1093/ofid/ofag611"
  vignette <- "Haynes_2026_cephalexin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "cephalexin",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "cephalexin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling on Vd/F (estimated exponent 0.92) and CL/F (estimated exponent 0.82), both normalized to a 56 kg reference FFM (Haynes 2026 Table 2 and Supplementary Appendix C 'Final structural and covariate-adjusted equations'). Supplementary Appendix A describes two FFM approaches: (1) Al-Sallami 2015 for children and adolescents with total body weight used below 3 years and the Janmahasatian 2005 equation from 18 years, and (2) the O'Hanlon 2023 neonate-to-adult method (Janmahasatian-type adult-referenced FFM scaled by an age- and sex-dependent maturation fraction for subjects < 20 years), with final FFM capped at 95% of body weight. The paper does not say explicitly which approach produced the final FFM column; the Table 1 neonatal medians (FFM 2.99 kg at 3.35 kg body weight) are below body weight, which is consistent with approach (2) rather than approach (1). Cohort FFM 1.99-77.18 kg (median 8.06 kg; Table 1).",
      source_name = "FFM"
    ),
    PAGE = list(
      description = "Postmenstrual age (gestational age at birth plus postnatal age)",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying, carried in WEEKS rather than the register-default months because the source maturation function and its 42.96-week TM50 are written in weeks (permitted by the PAGE register entry). Drives a sigmoid Hill maturation of CL/F: PMA^5.06 / (42.96^5.06 + PMA^5.06) (Haynes 2026 Supplementary Appendix C). Gestational age was recorded only for the neonatal study; a gestational age of 40 weeks was assigned to all older subjects (Supplementary Appendix A), so for children and adults PAGE = 40 + postnatal age in weeks and the maturation term is effectively 1.",
      source_name = "PMA"
    ),
    PNA = list(
      description = "Postnatal age (chronological age since birth)",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Enters only through a categorical cut: subjects with PNA <= 21 days have ka multiplied by exp(beta_PNA) = exp(-0.81) (Haynes 2026 Results 'Final Model', Table 2 row 'beta age < 21d_Ka' and Supplementary Appendix C, which prints the cut as 'PNA <= 21 days'). The canonical PNA column is carried in months, so model() recovers days as PNA * 30.4375 before applying the 21-day cut.",
      source_name = "PNA"
    )
  )

  # Screened during the stepwise covariate search (forward addition p < 0.05,
  # backward elimination p < 0.01; Haynes 2026 Methods 'Model Development' and
  # Supplementary Appendix B) but not retained. Also screened without a
  # canonical covariate column: lean body mass (Peters 2011), height,
  # formulation (capsule vs suspension) and feeding status.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as the body-size descriptor on Vd/F and CL/F normalized to 70 kg; FFM was preferred by 4.7 OFV points after PMA was added on CL/F (Haynes 2026 Supplementary Appendix C 'Body size effects'). Cohort 2.2-112 kg (median 10.95 kg; Table 1)."
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained as its own covariate; it enters only as a component of postmenstrual age (PAGE). Available for the neonatal study only; 40 weeks assigned to older subjects (Supplementary Appendix A)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (Haynes 2026 Methods 'Model Development'). Cohort 0.15-1.28 mg/dL (median 0.34; Table 1), with missing values imputed by regression on age (Supplementary Appendix A)."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (BSA-normalized)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Bedside Schwartz (k = 0.31 for infants <= 28 days, 0.357 male / 0.331 female for older infants, 0.413 for children < 18 years) or CKD-EPI 2021 race-free equation from 18 years (Supplementary Appendix A). On CL/F eGFR gave a 34.5-point OFV drop, less than PMA's 57 points, and added only 1.2 points on top of PMA, so it was not retained (Supplementary Appendix C)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 71L,
    n_studies = 4L,
    n_observations = 360L,
    age_range = "9 days to 29.2 years",
    age_median = "1.46 years",
    weight_range = "2.2-112 kg",
    weight_median = "10.95 kg",
    ffm_range = "1.99-77.18 kg",
    ffm_median = "8.06 kg",
    ga_range = "29 3/7 to 40 6/7 weeks gestational age at birth (neonatal study only; 40 weeks assigned to older subjects)",
    sex_female_pct = 39,
    disease_state = "Pooled: (1) 33 hospitalized young infants 9-56 days old receiving antibiotics at Children's Hospital Colorado (NCT04916951), (2) 15 children 9 months-17 years hospitalized with musculoskeletal infections at Children's Hospital Colorado given a single cephalexin dose (NCT03802552), (3) 12 children 1-16 years with musculoskeletal infections at CHU Sainte-Justine, Montreal, sampled after a steady-state dose, and (4) 11 healthy young adults 18-29 years in New Zealand given a single 1000 mg dose after a standardized meal.",
    dose_range = "8.93-51.18 mg/kg (median 24.6 mg/kg); 28.5-1500 mg per dose (median 550 mg). Oral suspension (47 subjects), capsule (20) or unknown (4).",
    renal_function = "eGFR median 101.34 mL/min/1.73 m^2 (range 37.80-187.91)",
    regions = "United States (Colorado), Canada (Montreal), New Zealand",
    notes = "Demographics in Haynes 2026 Table 1. 22 of 360 samples (6%) were below the limit of quantification (7 in Study 3 at <= 0.625 mg/L, 15 in Study 4 at <= 0.2 mg/L) and were handled with the Monolix censored likelihood (M4). The neonatal Study 1 dataset is the same one modelled on its own in Haynes_2025_cephalexin."
  )

  ini({
    # Structural parameters (Haynes 2026 Table 2). Vd/F and CL/F are the values
    # for a subject with 56 kg fat-free mass; CL/F is the fully mature value
    # (the Hill maturation term tends to 1 at high postmenstrual age).
    ltlag <- log(0.53); label("Absorption lag time (h)") # Table 2 'Tlag (h)' = 0.53 (RSE 13.16%)
    lka <- log(1.59); label("Absorption rate constant for PNA > 21 days (1/h)") # Table 2 'Ka (1/h)' = 1.59 (RSE 18.2%)
    lvc <- log(22.04); label("Apparent volume of distribution Vd/F for 56 kg FFM (L)") # Table 2 'Vd/F (L/56 kg)' = 22.04 (RSE 7.57%)
    lcl <- log(16.6); label("Mature apparent clearance CL/F for 56 kg FFM (L/h)") # Table 2 'CL/F (L/h/56 kg)' = 16.6 (RSE 4.25%)

    # Covariate effects (Haynes 2026 Table 2; all estimated)
    e_pnale21d_ka <- -0.81; label("Log-scale shift in ka for postnatal age <= 21 days (ka multiplied by exp(e_pnale21d_ka); unitless)") # Table 2 'beta age < 21d_Ka' = -0.81 (RSE 34.23%); Appendix C Ka = Ka_pop * exp(beta_PNA)
    e_ffm_vc <- 0.92; label("Power exponent of FFM/56 kg on Vd/F (unitless)") # Table 2 'beta FFM_V' = 0.92 (RSE 2.7%)
    e_ffm_cl <- 0.82; label("Power exponent of FFM/56 kg on CL/F (unitless)") # Table 2 'beta FFM_CL' = 0.82 (RSE 2.82%)
    lhill_mat <- log(5.06); label("Hill coefficient of the postmenstrual-age maturation of CL/F (unitless)") # Table 2 'Hill coefficient' = 5.06 (RSE 15.51%)
    ltmat50 <- log(42.96); label("Postmenstrual age at 50% of mature CL/F, TM50 (weeks)") # Table 2 'TM50 (weeks PMA)' = 42.96 (RSE 2.08%)

    # Inter-individual variability. Table 2 reports CV% of log-normal random
    # effects; Monolix CV = sqrt(exp(omega^2) - 1), so omega^2 = log(1 + CV^2).
    # The bootstrap medians in the same rows (0.77, 0.45, 0.16) are on the
    # omega (SD) scale and agree with sqrt(omega^2) = 0.75, 0.49, 0.17.
    # No IIV on Vd/F: fixed to zero for 53.8% shrinkage (Appendix C).
    etaltlag ~ 0.56484 # Table 2 'IIV for Tlag (CV%)' = 87.13%; log(1 + 0.8713^2)
    etalka ~ 0.23827 # Table 2 'IIV for Ka (CV%)' = 51.87%; log(1 + 0.5187^2)
    etalcl ~ 0.029322 # Table 2 'IIV for CL/F (CV%)' = 17.25%; log(1 + 0.1725^2)

    # Residual error: combined (Results 'Final Model'; Appendix C)
    addSd <- 0.31; label("Additive residual error (mg/L)") # Table 2 'Constant (a, mg/L)' = 0.31 (RSE 18.01%)
    propSd <- 0.21; label("Proportional residual error (fraction)") # Table 2 'Proportional (b)' = 0.21 (RSE 9.39%)
  })

  model({
    # Individual parameters, Haynes 2026 Supplementary Appendix C:
    #   Tlag = Tlag_pop * exp(eta_Tlag)
    #   Ka   = Ka_pop * exp(beta_PNA) * exp(eta_Ka)   if PNA <= 21 days
    #   Ka   = Ka_pop * exp(eta_Ka)                   if PNA >  21 days
    #   Vd/F = Vd_pop/F * (FFM / 56)^0.92
    #   CL/F = CL_pop/F * (FFM / 56)^0.82 * PMA^5.06 / (42.96^5.06 + PMA^5.06) * exp(eta_CL)
    tlag <- exp(ltlag + etaltlag)

    # PNA is carried in months; the source cut is 21 days (1 month = 30.4375 days).
    pnaDays <- PNA * 30.4375
    pnaLe21d <- 0
    if (pnaDays <= 21) pnaLe21d <- 1
    ka <- exp(lka + etalka + e_pnale21d_ka * pnaLe21d)

    hill_mat <- exp(lhill_mat)
    tmat50 <- exp(ltmat50)
    fmat <- PAGE^hill_mat / (tmat50^hill_mat + PAGE^hill_mat)

    vc <- exp(lvc) * (FFM / 56)^e_ffm_vc
    cl <- exp(lcl + etalcl) * (FFM / 56)^e_ffm_cl * fmat

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # CL/F and Vd/F are apparent oral parameters, so no bioavailability term.
    # Monolix 'combined1' residual error: SD = a + b * f.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
