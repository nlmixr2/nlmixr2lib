Baek_2014_cyclosporine <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for oral micro-emulsified cyclosporine in adult Korean living-donor renal transplant recipients during the first postoperative week (Baek 2014). The elimination rate constant k is estimated in place of clearance (CL/F = k * V2/F). Absorption changes with postoperative day (POD): relative to POD 2, ka rises 4.56-fold and the lag time 1.15-fold on POD 3, and the relative bioavailability rises 1.09-fold on POD 3 and then falls by a further age-dependent factor on POD 7 (0.807 at 41.2 years). Central volume scales with body weight and age. Between-subject variability is estimated on V2/F, V3/F, Q/F, ka and each POD change in ka and F; residual error is proportional on whole-blood cyclosporine."
  reference <- "Baek H, Han S, Yim DS, Kim SJ, Lee SY, Jang HR, Lee JE, Kim DJ, Kim YG, Oh HY, Huh W. Early non-steady-state population pharmacokinetics of oral cyclosporine in renal transplant recipients. Drug Des Devel Ther. 2014;8:2241. doi:10.2147/DDDT.S70595"
  vignette <- "Baek_2014_cyclosporine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "cyclosporine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cyclosporine", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "cyclosporine", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (WT / 60)^0.689 on V2/F (Table 3). The reference weight is 60 kg (the rounded cohort mean, 58.7 kg; Table 1).",
      source_name = "weight"
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects centred at the cohort mean age of 41.2 years (Table 1): (AGE / 41.2)^-0.227 on V2/F and (AGE / 41.2)^-0.293 on the POD 7 change in relative bioavailability (D7F1) (Table 3).",
      source_name = "age"
    ),
    POD = list(
      description = "Post-operative day (days since renal transplantation)",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = "The source treats POD as a three-level category (PK study days 2, 3 and 7; Table 3). The model derives two indicators from the day count: POD >= 3 switches on the POD 3 changes in ka, lag time and F, and POD >= 7 additionally switches on the POD 7 changes in ka and F, so POD 7 carries the POD 3 factors as well (the cumulative products of Table 3). Oral dosing starts on POD 2, so the model is defined for POD >= 2; days 4-6 were not sampled and take the POD 3 values, and days after 7 take the POD 7 values (extrapolation beyond the source). Time-varying: supply the day on every dose and observation row.",
      source_name = "POD"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "31 of 69 recipients female (Table 1). Tested as a discrete covariate (Methods, Covariate analysis); not retained."
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser status (CYP3A5*1 carrier)",
      units = "(binary)",
      type = "binary",
      notes = "20 expressers (*1/*1 or *1/*3) and 49 nonexpressers (*3/*3) (Table 2). The source codes nonexpressers as 1, the reverse of the canonical orientation (CYP3A5_EXPR = 1 - source code). Tested; not retained."
    ),
    SNP_ABCB1_RS1045642 = list(
      description = "ABCB1 (MDR1) c.3435C>T mutant allele carrier indicator",
      units = "(binary)",
      type = "binary",
      notes = "24 CC wild type (code 0) and 45 CT/TT carriers (code 1) (Table 2). A post hoc t-test found the individual POD 7 change in ka (D7KA) higher in carriers (P = 0.002), but adding the genotype as a covariate inflated the standard errors of the related estimates, so it was not retained (Results). CYP3A4 genotype was not tested because 67 of 69 recipients were *1/*1."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 69,
    n_studies = 1,
    n_observations = 1640,
    age_range = "16-63 years",
    age_mean = "41.2 years (SD 10.6)",
    weight_range = "38.9-97.6 kg",
    weight_mean = "58.7 kg (SD 10.2)",
    sex_female_pct = 44.9,
    race_ethnicity = "Korean (single centre, Seoul).",
    disease_state = "Adults (> 16 years) after a first living-donor renal transplantation (43 related, 26 unrelated donors). Underlying disease: glomerulonephritis 25, diabetes 9, ADPKD 3, renal tuberculosis 1, unknown 31. Pre-transplant renal replacement: pre-emptive 10, haemodialysis 48, CAPD 11; 2.1 (SD 2.4) years on dialysis.",
    dose_range = "Cyclosporine 4 mg/kg IV on POD 1, then micro-emulsified oral cyclosporine (Cipol-N soft capsule) 8 mg/kg/day in two divided doses from POD 2, adjusted to a trough target of 300-350 ng/mL. Co-medication: mycophenolate mofetil 1,500 mg/day and tapering methylprednisolone; no induction therapy.",
    regions = "Republic of Korea (Samsung Medical Center, Seoul).",
    co_medication = "Mycophenolate mofetil; methylprednisolone 500 mg peri-operatively and on POD 1, tapered 50% daily to POD 4, then 60 mg/day IV on POD 5-7.",
    genotype = "ABCB1 3435 CC 24 / CT+TT 45; CYP3A5 expressers 20 / nonexpressers 49; CYP3A4 *1/*1 67 / *1/*1B 2 (Table 2).",
    notes = "Full 12-h PK profiles after the morning dose on POD 2, 3 and 7 (pre-dose and 1, 2, 3, 4, 6, 8, 12 h); two patients with acute rejection had no POD 7 profile. Whole-blood cyclosporine by monoclonal EMIT (LLOQ 25 ng/mL). Demographics in Table 1, genotypes in Table 2."
  )

  ini({
    lkel <- log(0.254); label("Elimination rate constant k (1/h)") # Table 3: k = 0.254 1/h (RSE 3.11%)
    lvc <- log(142); label("Apparent central volume V2/F at 60 kg and 41.2 years, POD 2 bioavailability (L)") # Table 3: COV2 = 142 L (RSE 3.49%)
    e_wt_vc <- 0.689; label("Power exponent of (WT / 60) on V2/F (unitless)") # Table 3: EWTV = 0.689 (RSE 16.5%)
    e_age_vc <- -0.227; label("Power exponent of (AGE / 41.2) on V2/F (unitless)") # Table 3: EAGV = -0.227 (RSE 30.8%)
    lvp <- log(236); label("Apparent peripheral volume V3/F (L)") # Table 3: V3 = 236 L (RSE 11.4%)
    lq <- log(24.9); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3: Q = 24.9 L/h (RSE 6.73%)

    lka <- log(1.09); label("Absorption rate constant ka on POD 2 (1/h)") # Table 3: INKA = 1.09 1/h (RSE 16.5%)
    e_pod3_ka <- 4.56; label("Multiplicative change in ka on POD 3 relative to POD 2 (ratio)") # Table 3: D3KA = 4.56 (RSE 28.1%)
    e_pod7_ka <- fixed(1); label("Multiplicative change in ka on POD 7 relative to POD 3 (ratio)") # Table 3: D7KA = 1 (fixed), kept to carry its BSV

    ltlag <- log(0.811); label("Absorption lag time on POD 2 (h)") # Table 3: INLA = 0.811 h (RSE 3.52%)
    e_pod3_tlag <- 1.15; label("Multiplicative change in lag time on POD 3 relative to POD 2 (ratio)") # Table 3: D3LA = 1.15 (RSE 3.78%)

    lfdepot <- fixed(log(1)); label("Relative oral bioavailability on POD 2 (reference, fraction)") # Table 3: F1 = 1 on POD 2
    e_pod3_fdepot <- 1.09; label("Relative bioavailability on POD 3 versus POD 2 (ratio)") # Table 3: D3F1 = 1.09 (RSE 2.57%)
    e_pod7_fdepot <- 0.807; label("Relative bioavailability on POD 7 versus POD 3 at 41.2 years (ratio)") # Table 3: T7F1 = 0.807 (RSE 3.45%)
    e_age_pod7_fdepot <- -0.293; label("Power exponent of (AGE / 41.2) on the POD 7 bioavailability change (unitless)") # Table 3: EAGF = -0.293 (RSE 40.1%)

    # Table 3 reports between-subject variability as CV%; no definition is
    # given, so the exact log-normal conversion omega^2 = log(1 + CV^2) is used.
    # BSV on k, the lag time and D3LA was fixed to 0 (not estimable) and is
    # omitted.
    etalvc + etalvp ~ c(
      log(1 + 0.134^2),
      0.523 * sqrt(log(1 + 0.134^2) * log(1 + 0.671^2)),
      log(1 + 0.671^2)
    ) # Table 3: BSV V2 13.4%, V3 67.1%, correlation V2~V3 0.523
    etalq ~ log(1 + 0.354^2) # Table 3: BSV Q 35.4%
    etalka + etae_pod3_ka ~ c(
      log(1 + 1.36^2),
      -0.571 * sqrt(log(1 + 1.36^2) * log(1 + 1.74^2)),
      log(1 + 1.74^2)
    ) # Table 3: BSV INKA 136%, D3KA 174%, correlation INKA~D3KA -0.571
    etae_pod7_ka ~ log(1 + 2.68^2) # Table 3: BSV D7KA 268%
    etae_pod3_fdepot ~ log(1 + 0.128^2) # Table 3: BSV D3F1 12.8%
    etae_pod7_fdepot ~ log(1 + 0.212^2) # Table 3: BSV D7F1 21.2%

    propSd <- 0.23; label("Proportional residual error (fraction)") # Table 3: sigma1 (proportional) = 0.23 (RSE 2.18%)
  })
  model({
    # POD indicators. POD 7 carries the POD 3 factors as well: Table 3 gives
    # ka(POD 7) = INKA * D3KA * D7KA and F1(POD 7) = D3F1 * D7F1.
    pod3 <- POD >= 3
    pod7 <- POD >= 7

    kel <- exp(lkel)
    vc <- exp(lvc + etalvc) * (WT / 60)^e_wt_vc * (AGE / 41.2)^e_age_vc
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)

    ka <- exp(lka + etalka) *
      (e_pod3_ka * exp(etae_pod3_ka))^pod3 *
      (e_pod7_ka * exp(etae_pod7_ka))^pod7

    # Table 3 defines the lag time on POD 2 and POD 3 only; POD 7 keeps the
    # POD 3 value (no D7LA term).
    tlag <- exp(ltlag) * e_pod3_tlag^pod3

    fdepot <- exp(lfdepot) *
      (e_pod3_fdepot * exp(etae_pod3_fdepot))^pod3 *
      (e_pod7_fdepot * (AGE / 41.2)^e_age_pod7_fdepot * exp(etae_pod7_fdepot))^pod7

    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag
    f(depot) <- fdepot

    # Dose in mg, volume in L -> mg/L; x 1000 -> ng/mL (whole blood)
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
