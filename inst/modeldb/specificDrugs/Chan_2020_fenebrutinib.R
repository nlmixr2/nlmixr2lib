Chan_2020_fenebrutinib <- function() {
  description <- paste(
    "Three-compartment population PK model for oral fenebrutinib",
    "(GDC-0853), a noncovalent Bruton's tyrosine kinase inhibitor, pooled",
    "from two phase 1 studies in healthy volunteers and the phase 2 ANDES",
    "trial in rheumatoid arthritis patients (Chan 2020). Absorption is a",
    "flexible Savic transit-compartment absorption model (TCAM) that feeds",
    "the central compartment directly; the mean transit time carries fed,",
    "proton-pump-inhibitor (PPI) and PPI-by-fed effects, the number of",
    "transit compartments carries fed and tablet-formulation effects, and",
    "the relative bioavailability carries PPI and PPI-by-fed effects plus",
    "inter-individual and five-occasion inter-occasion variability.",
    "Apparent clearance carries PPI, age and healthy-volunteer effects.",
    "The proportional residual error decays exponentially with time after",
    "dose in healthy volunteers and is constant in patients. This is the",
    "PK model that produced the individual steady-state daily AUC driving",
    "the companion exposure-response models Chan_2020_fenebrutinib_acr and",
    "Chan_2020_fenebrutinib_das28."
  )
  reference <- paste(
    "Chan P, Yu J, Chinn L, Prohn M, Huisman J, Matzuka B, Hanley W,",
    "Tuckwell K, Quartino A. Population Pharmacokinetics, Efficacy",
    "Exposure-response Analysis, and Model-based Meta-analysis of",
    "Fenebrutinib in Subjects with Rheumatoid Arthritis. Pharm Res.",
    "2020;37(2):25. doi:10.1007/s11095-019-2752-y. (A correction notice",
    "revised only the article title; no parameter value is affected.)"
  )
  vignette <- "Chan_2020_fenebrutinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "fenebrutinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fenebrutinib", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "fenebrutinib", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL/F as exp((AGE / 48)^e_age_cl), exactly as printed in",
        "Model S1 ('CLAGE = ((AGE/48)**THETA(22))' then",
        "'CLCOV = EXP(CLAGE)*EXP(CLPAT)'). 48 years is the pooled",
        "analysis-set median (Table S2, 'All' column, median 48, range",
        "18-75). Because of the outer exp(), the factor at AGE = 48 is",
        "e = 2.718, not 1; see the in-file note on lcl and the vignette",
        "Assumptions and deviations section."
      ),
      source_name = "AGE"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator: 1 = healthy volunteer, 0 = rheumatoid arthritis patient",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (rheumatoid arthritis patient, the most common category)",
      notes = paste(
        "Model S1 derives HVPAT = 1 when PAT < 1 ('patients most common').",
        "Healthy volunteers are the 78 subjects of phase 1 studies GA29347",
        "and GP29832; patients are the 307 RA patients of phase 2 GA29350",
        "(Table S2). Multiplies CL/F by 1.52 and switches the residual",
        "error to the time-after-dose-dependent healthy-volunteer form."
      ),
      source_name = "PAT"
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor use at the dose: 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no PPI)",
      notes = paste(
        "Time-varying per Chan 2020 Methods ('Concomitant proton pump",
        "inhibitor (PPI) and food were categorized as time-varying",
        "covariates'). In phase 1 study GP29832 the PPI was rabeprazole",
        "given for 3 days before the fenebrutinib dose (Table S1); in the",
        "phase 2 trial 122 of 307 patients (39.7%) were PPI-comedicated",
        "(Table S2). Acts on CL/F, MTT and F1."
      ),
      source_name = "PPI"
    ),
    FED = list(
      description = "Fed state at the dose: 1 = fed, 0 = fasted",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Time-varying. Food co-administration in the phase 2 trial was",
        "self-reported (Table S1); 305 of 307 patients (99.3%) were",
        "classified as fed (Table S2). Model S1 applies FED (not the",
        "derived FOOD = FED * (PAT == 0) indicator, which it computes but",
        "never uses) to MTT, NTR and, jointly with PPI, to MTT and F1."
      ),
      source_name = "FED"
    ),
    FORM_FENEBRUTINIB_TABLET = list(
      description = "Fenebrutinib formulation: 1 = tablet, 0 = powder-in-capsule",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (powder-in-capsule)",
      notes = paste(
        "Model S1 sets TAB = 1 when FORM = 2. The powder-in-capsule",
        "formulation was given in the phase 1 multiple-ascending-dose study",
        "GA29347 and in one arm of GP29832 Part 1; every other arm, and the",
        "whole phase 2 trial, used the tablet (Table S1). Acts on the",
        "number of transit compartments only."
      ),
      source_name = "FORM"
    ),
    OCC = list(
      description = "Occasion index (1-5) for the inter-occasion variability on F1",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Model S1 defines five occasions, each carrying its own IOV eta on",
        "F1 (ETA(10) to ETA(14), a single reported variance). In the phase",
        "1 studies they map from the source OCC column (1, 2, 3, 11, 13);",
        "in phase 2 they map from the planned time after first dose,",
        "PTAFD = 144 h (occasion 1, day 7), 648 h (occasion 2, day 28",
        "trough), 650 / 653 / 657 h (occasion 3, day 28 post-dose), 1320 h",
        "(occasion 4, day 56) and 1992 h (occasion 5, day 84). Records",
        "outside occasions 1-5 carry no IOV. Supply 1-5 here; the mapping",
        "to the source codes is the user's responsibility."
      ),
      source_name = "OCC / PTAFD"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Part of the predefined covariate set screened in the popPK",
        "analysis (Table S2; pooled median 73 kg, range 38-153) but not",
        "retained in the final model. No estimate is reported."
      ),
      source_name = "WT"
    ),
    SEXF = list(
      description = "Sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened (Table S2; 255 of 385 female) but not retained in the final popPK model.",
      source_name = "SEX"
    ),
    CRP = list(
      description = "Baseline C-reactive protein",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Table S2) but not retained; not available for GA29347.",
      source_name = "CRP"
    ),
    CRCL = list(
      description = "Baseline creatinine clearance",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Table S2; median 111 mL/min) but not retained.",
      source_name = "CRCL"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 385L,
    n_studies = 3L,
    age_range = "18-75 years (pooled; median 48)",
    weight_range = "38-153 kg (pooled; median 73)",
    sex_female_pct = 66.2,
    race_ethnicity = "White 73.5%, Black 4.9%, American Indian or Alaska Native 6.8%, Asian 0.3%, multiple 1.3%, missing 13.2% (Table S2)",
    disease_state = "78 healthy volunteers (phase 1 GA29347 and GP29832) and 307 adults with moderate to severe active rheumatoid arthritis with inadequate response to methotrexate or anti-TNF therapy (phase 2 GA29350 / ANDES)",
    dose_range = "Phase 1: 20-250 mg BID or 500 mg QD powder-in-capsule for 14 days; 200 mg single dose (capsule or tablet, fasted or fed, with or without rabeprazole) and 200 mg BID tablet. Phase 2: 50 mg QD, 150 mg QD or 200 mg BID tablet for 12 weeks.",
    regions = "Phase 2: US, Eastern Europe and Latin America",
    n_observations = "4059 quantifiable concentrations (4565 records, 506 BLQ excluded)",
    notes = "Chan 2020 Results 'Population Pharmacokinetic (popPK) Modeling' and supplementary Tables S1-S2. Estimation FOCE-I in NONMEM 7.3."
  )

  ini({
    # Structural parameters. Chan 2020 Table I prints every log-transformed
    # THETA back-transformed (footnote: 'Log-transformed parameters have been
    # back-transformed'), so each value below is log(printed value).
    #
    # NOTE ON THE CLEARANCE SCALE. Model S1 multiplies CL by
    # EXP((AGE/48)**THETA(22)), which equals e = 2.718 at the reference age
    # rather than 1, so the 19.5 L/h of Table I is exp(THETA(1)) and NOT the
    # typical CL/F. The typical CL/F of a 48-year-old patient without PPI is
    # 19.5 * e = 53.0 L/h. This is encoded exactly as printed because two
    # independent published numbers only agree with the e-scaled reading:
    # the phase 2 AUC medians of Figure S3 and the phase 1 steady-state
    # half-life range of 4.2-9.9 h (see the vignette).
    lcl <- log(19.5); label("exp(THETA(1)) for apparent clearance CL/F (L/h); multiplied by e at the reference age, see comment")  # Table I, theta1 'CL/F: Apparent systemic clearance' = 19.5 L/h (95% CI 18.2-20.9)
    lvc <- log(381); label("Apparent central volume V1/F (L)")  # Table I, theta2 'V1/F' = 381 L (95% CI 333-437)
    lvp <- log(284); label("Apparent first peripheral volume V2/F (L)")  # Table I, theta3 'V2/F' = 284 L (95% CI 254-318)
    lq <- log(52.8); label("Apparent first intercompartmental clearance Q1/F (L/h)")  # Table I, theta4 'Q1/F' = 52.8 L/h (95% CI 44.1-63.2)
    lvp2 <- log(273); label("Apparent second peripheral volume V3/F (L)")  # Table I, theta5 'V3/F' = 273 L (95% CI 222-336)
    lq2 <- log(4.47); label("Apparent second intercompartmental clearance Q2/F (L/h)")  # Table I, theta6 'Q2/F' = 4.47 L/h (95% CI 3.73-5.36)
    lnn <- log(14.9); label("Number of transit compartments NTR (reference: capsule, fasted)")  # Table I, theta7 'NTR' = 14.9 (95% CI 13.0-17.1)
    lmtt <- log(0.849); label("Mean transit time MTT (h) (reference: fasted, no PPI)")  # Table I, theta8 'MTT' = 0.849 h (95% CI 0.755-0.954)
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1 (reference; unitless)")  # Model S1 'TVF1 = 1' (no THETA)

    # Covariate effects, all multiplicative as exp(theta * indicator) in
    # Model S1; values are log(back-transformed Table I value).
    e_fed_mtt <- log(1.43); label("Log fed effect on MTT")  # Table I, theta11 'Food effect on MTT' = 1.43 (95% CI 1.22-1.68)
    e_conmed_ppi_mtt <- log(0.835); label("Log PPI effect on MTT")  # Table I, theta12 'PPI effect on MTT' = 0.835 (95% CI 0.692-1.01)
    e_ppi_fed_mtt <- log(2.26); label("Log PPI-by-fed interaction effect on MTT")  # Table I, theta13 'PPI and food effect on MTT' = 2.26 (95% CI 1.78-2.86)
    e_conmed_ppi_f <- log(0.657); label("Log PPI effect on F1")  # Table I, theta14 'PPI effect on F1' = 0.657 (95% CI 0.568-0.759)
    e_ppi_fed_f <- log(0.693); label("Log PPI-by-fed interaction effect on F1")  # Table I, theta15 'PPI and food effect on F1' = 0.693 (95% CI 0.611-0.785)
    e_fed_nn <- log(0.864); label("Log fed effect on NTR")  # Table I, theta16 'Food effect on NTR' = 0.864 (95% CI 0.45-1.66)
    e_tablet_nn <- log(0.049); label("Log tablet-formulation effect on NTR")  # Table I, theta17 'Tablet effect on NTR' = 0.049 (95% CI 0.0287-0.0838)
    e_conmed_ppi_cl <- log(0.663); label("Log PPI effect on CL/F")  # Table I, theta20 'PPI effect on CL/F' = 0.663 (95% CI 0.650-0.675); Results '33.7%' decrease
    e_age_cl <- -0.161; label("Age exponent inside exp((AGE/48)^e_age_cl) on CL/F (unitless)")  # Table I, theta22 'Age effect on CL/F' = -0.161 (RSE 41.4%; 95% CI -0.291 to -0.0304); not log-transformed
    e_healthy_cl <- log(1.52); label("Log healthy-volunteer effect on CL/F")  # Table I, theta23 'Healthy volunteer effect on CL/F' = 1.52 (95% CI 1.27-1.82); Results '52% higher'

    # Inter-individual variability (Table I omega^2 values). Model S1 also
    # carries ETA(3)-ETA(7) on V2/F, Q1/F, V3/F, Q2/F and NTR, but Table I
    # reports no variance for them; they are omitted (equivalent to fixing
    # them at 0). See the vignette Assumptions and deviations section.
    etalcl ~ 0.0732 # Table I, omega1.1 'omega2 CL/F' = 0.0732 (RSE 14.9%; shrinkage 27.3%)
    etalvc ~ 0.100 # Table I, omega2.2 'omega2 V1/F' = 0.100 (RSE 22.6%; shrinkage 54.4%)
    etalmtt ~ 0.0861 # Table I, omega8.8 'omega2 MTT' = 0.0861 (RSE 29.7%; shrinkage 58.1%)
    etalfdepot ~ 0.131 # Table I, omega9.9 'omega2 F1' = 0.131 (RSE 17.1%; shrinkage 31.9%)

    # Inter-occasion variability on F1, five occasions (ETA(10)-ETA(14) in
    # Model S1) with one reported variance, encoded as the occasion-indicator
    # expansion because rxode2 cannot simulate the `eta ~ var | OCC` form.
    # Occasions 2-5 repeat occasion 1's variance ($OMEGA BLOCK(1) SAME).
    etaiov_f_1 ~ 0.299 # Table I, omega10.1 'omega2 IOV on F1' = 0.299 (RSE 5.30%; shrinkage 26.0%)
    etaiov_f_2 ~ fixed(0.299) # same variance as occasion 1 (SAME)
    etaiov_f_3 ~ fixed(0.299) # same variance as occasion 1 (SAME)
    etaiov_f_4 ~ fixed(0.299) # same variance as occasion 1 (SAME)
    etaiov_f_5 ~ fixed(0.299) # same variance as occasion 1 (SAME)

    # Residual error. Model S1 $ERROR: W = SQRT(PROPERR^2 * IPRED^2 +
    # THETA(10)^2). THETA(10) (additive) is absent from Table I and the
    # Results describe 'a proportional error model', so it is omitted
    # (additive SD 0). Patients: PROPERR = EXP(THETA(21)). Healthy volunteers:
    # PROPERR = EXP(THETA(9)) * (1 - ERRMAX * (1 - EXP(-EXP(THETA(18)) * TAD)))
    # with ERRMAX = EXP(THETA(19)) / (1 + EXP(THETA(19))).
    propSd <- 0.390; label("Proportional residual SD in patients (fraction)")  # Table I, theta21 'Proportional residual error in patients' = 0.390 (95% CI 0.372-0.408)
    propSd_hv <- 1.94; label("Proportional residual SD in healthy volunteers at time after dose 0 (fraction)")  # Table I, theta9 'Proportional residual error' = 1.94 (95% CI 1.29-2.92), exp(THETA(9))
    lkruv_hv <- log(1.94); label("Log rate of decline of the healthy-volunteer proportional residual SD with time after dose (1/h)")  # Table I, theta18 'Residual error rate in healthy volunteers' = 1.94 (95% CI 1.53-2.46)
    logitfruv_hv <- log(6.66); label("Logit of the maximum fractional decline of the healthy-volunteer proportional residual SD (unitless)")  # Table I, theta19 'Maximum residual error in healthy volunteers' = 6.66 (95% CI 4.21-10.6) = exp(THETA(19))
  })

  model({
    # Occasion indicators and IOV on F1 (Model S1 IOVF1).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_f <- oc1 * etaiov_f_1 + oc2 * etaiov_f_2 + oc3 * etaiov_f_3 +
      oc4 * etaiov_f_4 + oc5 * etaiov_f_5

    # Individual disposition parameters (Model S1 $PK).
    cl <- exp(lcl + etalcl) * exp(e_conmed_ppi_cl * CONMED_PPI) *
      exp((AGE / 48)^e_age_cl) * exp(e_healthy_cl * DIS_HEALTHY)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)
    vp2 <- exp(lvp2)
    q2 <- exp(lq2)

    # Transit-compartment absorption (Model S1 NTRCOV / MTTCOV / F1COV).
    nn <- exp(lnn) * exp(e_fed_nn * FED) * exp(e_tablet_nn * FORM_FENEBRUTINIB_TABLET)
    mtt <- exp(lmtt + etalmtt) * exp(e_fed_mtt * FED) *
      exp(e_conmed_ppi_mtt * CONMED_PPI) * exp(e_ppi_fed_mtt * CONMED_PPI * FED)
    fdepot <- exp(lfdepot + etalfdepot + iov_f) * exp(e_conmed_ppi_f * CONMED_PPI) *
      exp(e_ppi_fed_f * CONMED_PPI * FED)
    ktr <- (nn + 1) / mtt

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Savic transit input straight into the central compartment (Model S1
    # $DES INP1; there is no separate absorption depot and no ka). The
    # normaliser is the authors' Stirling approximation to log(NTR!)
    # (Model S1 LOGF1) rather than lgamma(NTR + 1). The dose is recorded on
    # `central` and f(central) <- 0 suppresses the bolus, so the density is
    # the only input; podo()/tad() carry the compartment argument because
    # the bare forms evaluate to 0 through the rxUi path. Like Model S1
    # (COM(1) / COM(3) reset at each dose), only the most recent dose feeds
    # the input.
    tstar <- tad(central)
    if (tstar < 0.001) tstar <- 0.001
    lognfac <- 0.5 * log(2 * 3.1415926 * nn) + nn * log(nn) - nn +
      log(1 + 1 / (12 * nn) + 1 / (288 * nn * nn))
    trInput <- exp(log(podo(central) * fdepot) + nn * log(ktr * tstar) + log(ktr) -
      ktr * tstar - lognfac)

    d/dt(central) <- trInput - (kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    f(central) <- 0

    # Concentration: dose in mg, volume in L -> mg/L; x1000 -> ng/mL.
    Cc <- 1000 * central / vc

    # Residual error (Model S1 $ERROR).
    errmax <- expit(logitfruv_hv)
    propSdHv <- propSd_hv * (1 - errmax * (1 - exp(-exp(lkruv_hv) * tstar)))
    propSdInd <- DIS_HEALTHY * propSdHv + (1 - DIS_HEALTHY) * propSd
    Cc ~ prop(propSdInd)
  })
}
