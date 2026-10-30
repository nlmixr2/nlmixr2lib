Wang_2024_leritrelvir <- function() {
  description <- paste(
    "Two-compartment population PK model for oral leritrelvir (RAY1216), a",
    "SARS-CoV-2 3CL protease inhibitor, given alone or with 100 mg ritonavir",
    "in healthy Chinese adults. Sequential zero-order input into the depot",
    "after an absorption lag, then first-order absorption; every structural",
    "parameter has a separate typical value with and without ritonavir.",
    "Female sex lowers CL/F, food raises F1, ka and the zero-order duration,",
    "and in repeated-dose regimens CL/F steps down 12 h after the first dose",
    "(auto-inhibition) and back up 96 h after it (auto-induction), with",
    "occasion-gated random effects on CL/F and Vc/F."
  )
  reference <- paste(
    "Wang K, Li H, Li Y, Xu F, Sun Z, Yang Y, Huang J, Chen X.",
    "Model informed dose regimen optimizing in development of leritrelvir",
    "for the treatment of mild or moderate COVID-19.",
    "Front Pharmacol. 2024;15:1449583. doi:10.3389/fphar.2024.1449583"
  )
  vignette <- "Wang_2024_leritrelvir"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CONMED_RTV = list(
      description = "Concomitant ritonavir 100 mg (pharmacokinetic booster) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (leritrelvir alone)",
      notes = paste(
        "Selects between the two independently estimated sets of structural",
        "typical values in Wang 2024 Table 2 ('Mono' vs 'Comb'); the",
        "supplementary NONMEM control stream branches on IF (RTV.EQ.0). The",
        "random effects, covariate effects and residual error are shared",
        "between the two sets. Ritonavir was given as 100 mg with each",
        "leritrelvir dose (single dose in the SAD/food-effect arms, BID in",
        "the MAD arms)."
      ),
      source_name = "RTV"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Enters CL/F as exp(-0.224 * SEX) in the control stream. The paper",
        "does not print the SEX coding, but it reports that female subjects",
        "had 12-51% HIGHER exposure than males (Results 3.3.2; Table S3),",
        "which with a negative coefficient on CL/F requires SEX = 1 for",
        "female. The source column therefore maps one-to-one onto SEXF."
      ),
      source_name = "SEX"
    ),
    FED = list(
      description = "Fed-state dosing indicator (1 = after a high-fat, high-calorie standard meal, 0 = after an overnight fast of at least 10 h)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Only the crossover food-effect arms (FE-1: 800 mg alone; FE-2:",
        "500 mg with ritonavir) contributed fed doses; every other arm was",
        "dosed fasted. Enters F1, ka and the zero-order duration D1 as",
        "exp(theta * FED), and switches the residual error to the 70.2%",
        "fed-absorption proportional SD within 2 h after a fed dose."
      ),
      source_name = "FED"
    ),
    MULTI_DOSE = list(
      description = "Repeated-dose regimen indicator (1 = subject on a multiple-dose regimen, 0 = single dose)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (single dose)",
      notes = paste(
        "The control stream gates both occasion indicators on DOSESCH.GT.1",
        "(multiple-dose regimen): OC1 = 1 when DOSESCH > 1 and TIME >= 12 h,",
        "OC2 = 1 when DOSESCH > 1 and TIME >= 96 h, where TIME is time since",
        "the first dose of the regimen. Here oc1 = MULTI_DOSE * (t >= 12) and",
        "oc2 = MULTI_DOSE * (t >= 96), so the event table must start at the",
        "first dose (t = 0). Single-dose subjects keep MULTI_DOSE = 0 at all",
        "times, which stops the occasion terms switching on during their",
        "post-dose sampling."
      ),
      source_name = "DOSESCH"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "leritrelvir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "leritrelvir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "leritrelvir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 72,
    n_studies = 1,
    age_range = "19-48 years",
    age_median = "35 years",
    weight_range = "48.2-90.8 kg",
    weight_median = "60.6 kg",
    sex_female_pct = 50,
    race_ethnicity = "Chinese (single-country phase I study in China; race not tabulated)",
    disease_state = "Healthy adult volunteers",
    dose_range = paste(
      "Leritrelvir alone: 400, 800 and 1600 mg single dose; 200 mg single dose",
      "(DDI period 1); 800 mg single dose fasted and fed; 400 mg TID for 4",
      "days plus a morning dose on day 5. With 100 mg ritonavir: 200 mg single",
      "dose; 500 mg single dose fasted and fed; 300, 500 and 700 mg BID for 4",
      "days plus a morning dose on day 5."
    ),
    regions = "China",
    n_observations = 1571,
    notes = paste(
      "Phase I SAD / MAD / food-effect / DDI study NCT05829551 (Wang 2024",
      "Figure 1, Table 1). 1571 leritrelvir plasma concentrations from 72",
      "subjects; 150 of 1721 (8.72%) samples below the LLOQ were excluded.",
      "The 30 crossover food-effect subjects were coded as two NONMEM IDs",
      "each (102 IDs): 78 fasted / 24 fed ID-periods and 54 mono / 48",
      "ritonavir ID-periods. Median albumin 43.9 g/L, CrCL 115 mL/min."
    )
  )

  ini({
    # Structural typical values, leritrelvir ALONE (CONMED_RTV = 0)
    # Wang 2024 Table 2 'Mono' column; NONMEM THETA(1), (3), (5), (7), (9), (11), (13)
    lcl <- log(44.4); label("Apparent clearance CL/F, leritrelvir alone (L/h)") # Table 2: CL Mono = 44.4 L/h (RSE 6.68%)
    lvc <- log(18.9); label("Apparent central volume Vc/F, leritrelvir alone (L)") # Table 2: Vc Mono = 18.9 L (RSE 18%); row mislabelled 'peripheral compartment, Vc'
    lka <- log(0.293); label("First-order absorption rate constant, leritrelvir alone (1/h)") # Table 2: Ka Mono = 0.293 1/h (RSE 6.95%)
    lq <- log(5.57); label("Apparent intercompartmental clearance Q/F, leritrelvir alone (L/h)") # Table 2: Q Mono = 5.57 L/h (RSE 6.47%)
    lvp <- log(58.2); label("Apparent peripheral volume Vp/F, leritrelvir alone (L)") # Table 2: Vp Mono = 58.2 L (RSE 7.49%)
    ld1 <- log(0.0722); label("Duration of zero-order input into the depot D1, leritrelvir alone (h)") # Table 2: D1 Mono = 0.0722 h (RSE 38%)
    ltlag <- log(0.189); label("Absorption lag time ALAG1, leritrelvir alone (h)") # Table 2: ALAG1 Mono = 0.189 h (RSE 2.33%)

    # Structural typical values, leritrelvir WITH 100 mg ritonavir (CONMED_RTV = 1)
    # Wang 2024 Table 2 'Comb' column; NONMEM THETA(2), (4), (6), (8), (10), (12), (14)
    lcl_boosted <- log(12.2); label("Apparent clearance CL/F, with ritonavir (L/h)") # Table 2: CL Comb = 12.2 L/h (RSE 8.02%)
    lvc_boosted <- log(9.78); label("Apparent central volume Vc/F, with ritonavir (L)") # Table 2: Vc Comb = 9.78 L (RSE 27.1%)
    lka_boosted <- log(0.169); label("First-order absorption rate constant, with ritonavir (1/h)") # Table 2: Ka Comb = 0.169 1/h (RSE 7.83%)
    lq_boosted <- log(0.611); label("Apparent intercompartmental clearance Q/F, with ritonavir (L/h)") # Table 2: Q Comb = 0.611 L/h (RSE 15.8%)
    lvp_boosted <- log(9.30); label("Apparent peripheral volume Vp/F, with ritonavir (L)") # Table 2: Vp Comb = 9.30 L (RSE 11.5%)
    ld1_boosted <- log(0.0684); label("Duration of zero-order input into the depot D1, with ritonavir (h)") # Table 2: D1 Comb = 0.0684 h (RSE 49.8%)
    ltlag_boosted <- log(0.166); label("Absorption lag time ALAG1, with ritonavir (h)") # Table 2: ALAG1 Comb = 0.166 h (RSE 3.22%)

    lfdepot <- fixed(log(1)); label("Relative bioavailability F1 in the fasted state (unitless)") # Table 2: 'Relative bioavailability, F1' = 1, FIX

    # Covariate effects, all on the log scale: parameter * exp(theta * covariate).
    # Values are the exponents printed in the Results 3.2.1 equations; Table 2
    # lists exp(theta).
    e_fed_fdepot <- 0.409; label("Food effect on F1, log-scale (unitless)") # Results 3.2.1 equation F1 = exp(0.409 * FED); Table 2 F1.FED = 1.51 (RSE 5.9%)
    e_fed_ka <- 0.291; label("Food effect on ka, log-scale (unitless)") # Results 3.2.1 equation Ka = 0.293 * exp(0.291 * FED); control stream KA uses THETA(17); Table 2 row labelled D1.FED = 1.34 (RSE 12.6%)
    e_fed_d1 <- 3.14; label("Food effect on the zero-order duration D1, log-scale (unitless)") # Results 3.2.1 equation D1 = 0.0722 * exp(3.14 * FED); control stream D1 uses THETA(16); Table 2 row labelled Ka.FED = 23.2 (RSE 113%)
    e_sexf_cl <- -0.224; label("Female-sex effect on CL/F, log-scale (unitless)") # Results 3.2.1 equation exp(-0.224 * SEX); Table 2 CL.sex = 0.800 (RSE 7.62%)
    e_occ1_cl <- -0.469; label("Auto-inhibition step on CL/F from 12 h after the first dose of a repeated-dose regimen, log-scale (unitless)") # Results 3.2.1 equation exp(-0.469 * OC1); Table 2 CL.OC1.INH = 0.626 (RSE 7.84%)
    e_occ2_cl <- 0.412; label("Auto-induction step on CL/F from 96 h after the first dose of a repeated-dose regimen, log-scale (unitless)") # Results 3.2.1 equation exp(0.412 * OC2); Table 2 CL.OC2.IND = 1.51 (RSE 6.18%)

    # Inter-individual variability. Table 2 prints IIV (%) = 100 * sqrt(omega^2);
    # the Wald 95% CIs are symmetric on the omega^2 scale, e.g. CL 30.6 [23.6-36.2]:
    # (23.6^2 + 36.2^2) / 2 = 933.7 = 30.6^2 to rounding, whereas log(1 + CV^2)
    # misses the midpoint. Hence omega^2 = (IIV% / 100)^2.
    etalcl ~ 0.093636 # Table 2: IIV CL = 30.6% -> 0.306^2
    etalvc ~ 0.652864 # Table 2: IIV Vc = 80.8% -> 0.808^2
    etalka ~ 0.068121 # Table 2: IIV Ka = 26.1% -> 0.261^2
    etalvp ~ 0.199809 # Table 2: IIV Vp = 44.7% -> 0.447^2
    etald1 ~ 2.4025 # Table 2: IIV D1 = 155% -> 1.55^2

    # Occasion-gated random effects (control stream ETA(8)-ETA(11)). Each is a
    # per-subject eta switched on by oc1 / oc2; after 96 h both apply. Table 2
    # prints only the OC1 rows; the OC2 rows are blank, read as NONMEM $OMEGA
    # SAME, so the OC2 variances repeat the OC1 values.
    etaiov_cl_1 ~ 0.037249 # Table 2: OC1.CL = 19.3% -> 0.193^2
    etaiov_cl_2 ~ 0.037249 # Table 2: OC2.CL row blank; SAME as OC1.CL
    etaiov_vc_1 ~ 0.295936 # Table 2: OC1.Vc = 54.4% -> 0.544^2
    etaiov_vc_2 ~ 0.295936 # Table 2: OC2.Vc row blank; SAME as OC1.Vc

    # Residual error, proportional (supplement equation 2)
    propSd <- 0.244; label("Proportional residual error (fraction)") # Table 2: proportional residual = 24.4% (RSE 2.06%)
    propSd_fedearly <- 0.702; label("Proportional residual error within 2 h after a fed dose (fraction)") # Table 2: proportional error of fed subjects within 2 h after dose = 70.2% (RSE 16.8%)
  })

  model({
    # Typical values switch on ritonavir co-administration (control stream
    # IF (RTV.EQ.0) ... ELSE ... ENDIF)
    tvlcl <- lcl * (1 - CONMED_RTV) + lcl_boosted * CONMED_RTV
    tvlvc <- lvc * (1 - CONMED_RTV) + lvc_boosted * CONMED_RTV
    tvlka <- lka * (1 - CONMED_RTV) + lka_boosted * CONMED_RTV
    tvlq <- lq * (1 - CONMED_RTV) + lq_boosted * CONMED_RTV
    tvlvp <- lvp * (1 - CONMED_RTV) + lvp_boosted * CONMED_RTV
    tvld1 <- ld1 * (1 - CONMED_RTV) + ld1_boosted * CONMED_RTV
    tvltlag <- ltlag * (1 - CONMED_RTV) + ltlag_boosted * CONMED_RTV

    # Occasion indicators: t is time since the regimen's first dose
    oc1 <- MULTI_DOSE * (t >= 12)
    oc2 <- MULTI_DOSE * (t >= 96)

    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2

    cl <- exp(tvlcl + etalcl + iov_cl) *
      exp(e_occ1_cl * oc1 + e_occ2_cl * oc2 + e_sexf_cl * SEXF)
    vc <- exp(tvlvc + etalvc + iov_vc)
    ka <- exp(tvlka + e_fed_ka * FED + etalka)
    q <- exp(tvlq)
    vp <- exp(tvlvp + etalvp)
    d1 <- exp(tvld1 + e_fed_d1 * FED + etald1)
    tlag <- exp(tvltlag)
    fdepot <- exp(lfdepot + e_fed_fdepot * FED)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Zero-order input into the depot over d1 after the lag: dose records
    # must carry rate = -2 for rxode2 to use the modelled duration.
    alag(depot) <- tlag
    dur(depot) <- d1
    f(depot) <- fdepot

    # mg / L * 1000 = ng/mL (control stream S2 = V2 / 1000)
    Cc <- 1000 * central / vc

    # Control stream $ERROR: IF (TFDS.LE.2.AND.FED.EQ.1) the second EPS applies,
    # TFDS being time after the previous dose.
    early_flag <- tad() <= 2
    fedearly <- FED * early_flag
    propSdEff <- propSd_fedearly * fedearly + propSd * (1 - fedearly)
    Cc ~ prop(propSdEff)
  })
}
