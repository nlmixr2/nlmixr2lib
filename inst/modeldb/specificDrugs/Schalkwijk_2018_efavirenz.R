Schalkwijk_2018_efavirenz <- function() {
  description <- "Mechanism-based two-compartment population PK model for oral efavirenz in pregnant and non-pregnant women living with HIV (Schalkwijk 2018): three-transit-compartment absorption, a well-stirred liver giving systemic hepatic clearance and first-pass extraction from intrinsic clearance, CYP2B6 metabolizer phenotype (poor / intermediate / extensive) on intrinsic clearance, gestational-age polynomials for serum albumin (driving the unbound fraction) and haematocrit (driving hepatic plasma flow), pregnancy effects on mean absorption time and relative bioavailability, and allometric scaling to 70 kg."
  reference <- paste(
    "Schalkwijk S, ter Heine R, Colbers AC, Huitema ADR, Denti P, Dooley KE,",
    "Capparelli E, Best BM, Cressey TR, Greupink R, Russel FGM, Mirochnick M,",
    "Burger DM. A Mechanism-Based Population Pharmacokinetic Analysis Assessing",
    "the Feasibility of Efavirenz Dose Reduction to 400 mg in Pregnant Women.",
    "Clin Pharmacokinet. 2018;57(11):1421-1433. doi:10.1007/s40262-018-0642-9",
    sep = " "
  )
  vignette <- "Schalkwijk_2018_efavirenz"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of every flow (CLint/F, Q/F and the literature hepatic blood flow Qhep; exponent 0.75) and every volume (Vc/F, Vp/F; exponent 1) to a 70 kg reference (Methods 2.3: 'all flow parameters and volumes were scaled to a total non-pregnant body weight of 70 kg according to allometric theory. The allometric exponents were fixed to 3/4 for flow parameters and 1 for volumes of distribution'). The paper's Table 3 simulations drew weight from a log-normal distribution with geometric mean 62 kg (geometric SD 1.3) for both the pregnant and the non-pregnant cohort.",
      source_name = "WT"
    ),
    EGA = list(
      description = "Maternal estimated gestational age (weeks); 0 for non-pregnant women",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Argument of the gestational polynomials for serum albumin (Eq 5) and haematocrit (Eq 6), which return their non-pregnant values at EGA = 0. Set EGA = 0 on every non-pregnant (including postpartum) record. Table 1 labels the column 'Median gestational age at sampling times, years', but the values (21-39) and Methods 2.7 ('34 +/- 2.3 weeks') show the unit is weeks. The Table 3 third-trimester simulations drew EGA from a normal distribution with mean 34 and SD 2.3 weeks.",
      source_name = "GA"
    ),
    PREG = list(
      description = "Pregnancy status indicator (1 = pregnant, 0 = non-pregnant)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-pregnant, including postpartum)",
      notes = "Dichotomous pregnancy covariate retained on mean absorption time (MAT 1.67 h pregnant vs 2.12 h non-pregnant) and relative bioavailability (116% of non-pregnant) after forward inclusion / backward elimination (Methods 2.4; Results paragraph 4; Table 2). PREG = 1 must coincide with EGA > 0: the albumin and haematocrit effects are carried by EGA, the absorption effects by PREG. Time-varying in the pooled crossover studies (pregnant visit vs postpartum visit).",
      source_name = "pregnancy"
    ),
    CYP2B6_PM = list(
      description = "CYP2B6 poor-metabolizer indicator (1 = two loss-of-function alleles at 516G>T / 983T>C, pooling the slow and very slow genotype classes)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extensive metabolizer is the reference when CYP2B6_PM = CYP2B6_IM = 0)",
      notes = "Selects the poor-metabolizer typical CLint/F of 1380 L/h (Table 2). The paper fitted a three-subpopulation mixture on CLint/F with population frequencies fixed to 14% PM / 36% IM / 50% EM (Results paragraph 2); genotype was available for only 18% of the cohort (Methods 2.5, Study 7 only), and subjects without genotype were assigned to the subpopulation with the highest individual probability. Here the subpopulation is supplied as a covariate. The paper's PM level pools the register's strict CYP2B6_SM and CYP2B6_USM levels (Methods 2.2), so do not supply a strict-SM column to this model.",
      source_name = "PM (mixture subpopulation)"
    ),
    CYP2B6_IM = list(
      description = "CYP2B6 intermediate-metabolizer indicator (1 = a single variant allele at 516 or 983)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extensive metabolizer is the reference when CYP2B6_PM = CYP2B6_IM = 0)",
      notes = "Selects the intermediate-metabolizer typical CLint/F of 3340 L/h (Table 2). Same genotype definition as the register's CYP2B6_IM (Methods 2.2: 'IMs (single variant allele at position 516 or 983)'). Mixture-model frequency fixed to 36% (Results paragraph 2).",
      source_name = "IM (mixture subpopulation)"
    ),
    OCC = list(
      description = "Pharmacokinetic assessment visit (occasion) index, 1-6",
      units = "(count)",
      type = "categorical",
      reference_category = "n/a -- decomposed inside model() into six binary indicators multiplied against the per-occasion etaiov_lfdepot_<k> slots",
      notes = "Inter-occasion variability is on relative bioavailability F only (Results paragraph 1). Methods 2.3: 'Each occasion was defined as a pharmacokinetic assessment visit, ranging from 1 to 6 assessments per individual', so six slots are encoded; slots 2-6 share the occasion-1 variance (NONMEM $OMEGA BLOCK(1) SAME idiom) and are fixed() to it. OCC = 0 on every record switches inter-occasion variability off.",
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "efavirenz", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "efavirenz", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "efavirenz", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "efavirenz", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "efavirenz", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "efavirenz", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 247L,
    n_studies = 6L,
    age_range = "adult women (age not tabulated)",
    weight_range = "37-132 kg (Table 1 per-study ranges)",
    weight_median = "53-83 kg across studies and pregnancy states (Table 1)",
    sex_female_pct = 100,
    race_ethnicity = "Pooled: Black (Studies 3 and 7, 100%), Caucasian (Study 4, 100%), Thai/Caucasian (Study 5), Hispanic/non-Hispanic (Study 6), mixed international (Study 2); Table 1",
    disease_state = "HIV-1 infection on efavirenz-based antiretroviral therapy; women in the second or third trimester of pregnancy and non-pregnant (including postpartum) women",
    dose_range = "Efavirenz 600 mg once daily orally (all but five women; 300, 400 or 800 mg in Studies 4 and 5)",
    regions = "Pooled international (South Africa, Thailand, South America, USA, Europe)",
    gestational_age_range = "21-39 weeks at sampling (Table 1)",
    n_samples = "1968 plasma samples in total, 774 during pregnancy (including the 11-woman external-evaluation Study 1)",
    notes = "Model development used studies 2-7 (247 women); Study 1 (11 women, 219 samples) was held out for external evaluation, giving the 258 women and 1968 samples quoted in the Abstract. 116 women were sampled only when not pregnant; 142 contributed both pregnant and postpartum samples. CYP2B6 genotype was available for 18% of women (Study 7: 10 PM, 25 IM, 26 EM). Women on interacting co-medication (e.g. rifampicin, isoniazid) were excluded. NONMEM 7.3.0, FOCE-I. Table 1."
  )

  ini({
    # Absorption. Three transit compartments with n + 1 = 4 identical rate
    # constants ktr = (n + 1) / MAT (Eq 7; Fig 1 caption).
    lmtt <- log(2.12); label("Mean absorption time MAT in non-pregnant women (h)") # Table 2 row 'MAT (h)' = 2.12
    lmtt_preg <- log(1.67); label("Mean absorption time MAT in pregnant women (h)") # Table 2 row 'MAT (h) in pregnant women' = 1.67

    # Relative bioavailability (typical non-pregnant F fixed to 1, Methods 2.3)
    lfdepot <- fixed(log(1)); label("Relative bioavailability F in non-pregnant women (unitless)") # Methods 2.3: 'the typical value of bioavailability was fixed to 1'
    e_preg_fdepot <- 0.16; label("Fractional increase in relative bioavailability during pregnancy (unitless)") # Table 2 row 'F (%) relative to non-pregnant women' = 116, i.e. 1.16 - 1

    # Intrinsic hepatic clearance by CYP2B6 metabolizer subpopulation, 70 kg
    lclint <- log(4580); label("Apparent intrinsic hepatic clearance CLint/F, extensive metabolizer, 70 kg (L/h)") # Table 2 row 'CLint/F (L/h) Extensive' = 4580
    lclint_im <- log(3340); label("Apparent intrinsic hepatic clearance CLint/F, intermediate metabolizer, 70 kg (L/h)") # Table 2 row 'CLint/F (L/h) Intermediate' = 3340
    lclint_pm <- log(1380); label("Apparent intrinsic hepatic clearance CLint/F, poor metabolizer, 70 kg (L/h)") # Table 2 row 'CLint/F (L/h) Poor' = 1380

    # Distribution, 70 kg
    lvc <- log(133); label("Apparent central volume Vc/F, 70 kg (L)") # Table 2 row 'Vc/F (L)' = 133
    lvp <- log(390); label("Apparent peripheral volume Vp/F, 70 kg (L)") # Table 2 row 'Vp/F (L)' = 390
    lq <- log(35); label("Apparent intercompartmental clearance Q/F, 70 kg (L/h)") # Table 2 row 'Q/F (L/h)' = 35

    # Fixed physiology of the well-stirred liver (Eqs 1-4)
    lqh <- fixed(log(109)); label("Hepatic blood flow Qhep, non-pregnant, 70 kg (L/h)") # Methods 2.2: Qhep 'fixed to the literature values (109 L/h) for non-pregnant women'
    lkd <- fixed(log(2.05)); label("Efavirenz-albumin dissociation constant kD (uM)") # Methods 2.2: kD 'fixed to the in vitro literature value, 2.05 uM'

    # Allometric exponents (Methods 2.3)
    e_wt_cl <- fixed(0.75); label("Allometric exponent on flow parameters CLint/F, Q/F and Qhep (unitless)") # Methods 2.3: 'fixed to 3/4 for flow parameters'
    e_wt_vc <- fixed(1); label("Allometric exponent on volumes Vc/F and Vp/F (unitless)") # Methods 2.3: 'and 1 for volumes of distribution'

    # IIV, reported as %CV; omega^2 = log(1 + CV^2). The CLint/MAT correlation
    # (6%) was not included (Results paragraph 1).
    etalclint ~ 0.09752 # Table 2 row 'IIV CLint/F (%)' = 32; log(1 + 0.32^2) = 0.09752
    etalmtt ~ 0.17700 # Table 2 row 'IIV MAT (%)' = 44; log(1 + 0.44^2) = 0.17700

    # IOV on F, one slot per assessment visit (Methods 2.3: 1 to 6 visits)
    etaiov_lfdepot_1 ~ 0.05600 # Table 2 row 'IOV F (%)' = 24; log(1 + 0.24^2) = 0.05600
    etaiov_lfdepot_2 ~ fixed(0.05600) # same magnitude, occasion 2
    etaiov_lfdepot_3 ~ fixed(0.05600) # same magnitude, occasion 3
    etaiov_lfdepot_4 ~ fixed(0.05600) # same magnitude, occasion 4
    etaiov_lfdepot_5 ~ fixed(0.05600) # same magnitude, occasion 5
    etaiov_lfdepot_6 ~ fixed(0.05600) # same magnitude, occasion 6

    propSd <- 0.18; label("Proportional residual error on total plasma efavirenz (fraction)") # Table 2 row 'Proportional residual error (%)' = 18
  })

  model({
    # 1. Gestational physiology (EGA in weeks; EGA = 0 is non-pregnant)
    alb <- (45.8 - 0.1775 * EGA - 0.0033 * EGA^2) / 0.07 # Eq 5, serum albumin [P] (uM); 0.07 g/L per uM
    hct <- 39.1 - 0.0544 * EGA - 0.0021 * EGA^2 # Eq 6, haematocrit (v/v %)

    # Inter-occasion term on F
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    iov_fdepot <- oc1 * etaiov_lfdepot_1 + oc2 * etaiov_lfdepot_2 +
      oc3 * etaiov_lfdepot_3 + oc4 * etaiov_lfdepot_4 +
      oc5 * etaiov_lfdepot_5 + oc6 * etaiov_lfdepot_6

    # CYP2B6 subpopulation; extensive metabolizer when both indicators are 0
    cyp2b6_em <- 1 - CYP2B6_PM - CYP2B6_IM
    lclint_pheno <- cyp2b6_em * lclint + CYP2B6_IM * lclint_im + CYP2B6_PM * lclint_pm

    # 2. Individual parameters
    clint <- exp(lclint_pheno + etalclint) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    qh <- exp(lqh) * (WT / 70)^e_wt_cl
    kd <- exp(lkd)
    mtt <- exp(lmtt * (1 - PREG) + lmtt_preg * PREG + etalmtt)
    fdepot <- exp(lfdepot + iov_fdepot) * (1 + e_preg_fdepot * PREG)

    # Well-stirred liver (Eqs 1-4). The liver is algebraic: systemic hepatic
    # clearance qhp * eh acts on central, and the absorbed dose loses the
    # fraction eh on its first pass. Vh is not reported (see vignette).
    fu <- kd / (kd + alb) # Eq 4
    qhp <- (1 - hct / 100) * qh # Eq 3, hepatic plasma flow (L/h)
    eh <- clint * fu / (qhp + clint * fu) # Eq 2
    cl <- qhp * eh # Eq 1

    # 3. Micro-constants
    ktr <- (3 + 1) / mtt # Eq 7 with n = 3 transit compartments
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(central) <- (1 - eh) * ktr * transit3 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Bioavailability
    f(depot) <- fdepot

    # 6. Observation and error. Cu is the predicted unbound concentration
    # (Eq 9), not an observed output.
    Cc <- central / vc
    Cu <- fu * Cc
    Cc ~ prop(propSd)
  })
}
