Huang_2026_busulfan <- function() {
  description <- "One-compartment IV popPK model for busulfan in Chinese pediatric allogeneic hematopoietic stem cell transplantation recipients (5 centres), with normal-fat-mass allometric scaling (separately estimated fat fractions for CL and V, exponents 3/4 and 1), a postmenstrual-age sigmoid maturation function on CL, a GSTA1 slow-metabolizer factor on CL, and between-occasion variability on CL (Huang 2026)."
  reference <- "Huang C, Du X, Xue L, Song L, Liu X, Chen Z, Wang X, Qin M, Sun X, Wang D, Yang J, Li D, Jiang H, Ding W, Qu Y, Zhu Z, Zhao L, Yuan J, Dong J, Wu D, Hu S, Miao L. Association of Busulfan Exposure with Event-Free Survival in a Chinese Pediatric Cohort: A Multicenter, Prospective Cohort Study. Drug Des Devel Ther. 2026;20:593835. doi:10.2147/DDDT.S593835"
  vignette <- "Huang_2026_busulfan"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "busulfan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Actual body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed at baseline. Combined with FFM to build the parameter-specific normal fat mass, NFM = FFM + Ffat * (WT - FFM) (Supplementary Equation S2). The standard subject is a 70 kg, 176 cm adult male (Table 1 footnote).",
      source_name = "WT"
    ),
    FFM = list(
      description = "Fat-free mass (Janmahasatian 2005 semi-mechanistic equation, derived from body weight, height and sex)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed at baseline. Supplementary Equation S1: FFM = WHSmax * HT^2 * WT / (WHS50 * HT^2 + WT), HT in metres, with WHSmax = 42.92 and WHS50 = 30.93 kg/m^2 for males and 37.99 and 35.98 kg/m^2 for females. The adult equation is applied to children without a paediatric correction. FFMSTD = 56.1 kg is the male equation at 70 kg and 1.76 m (Table 1 footnote).",
      source_name = "FFM"
    ),
    AGE = list(
      description = "Postnatal age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Combined with GA inside model() to form postmenstrual age, PMA (weeks) = AGE * 52 + GA. The paper names PMA as the maturation covariate (Results) but does not print its derivation; the form used by the same group's predecessor model (Du 2022, Supplementary Equation S6) is adopted.",
      source_name = "Age"
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Second component of postmenstrual age. Not reported for this cohort; supply 40 weeks (term birth) when unknown."
    ),
    GSTA1_RM = list(
      description = "GSTA1 rapid-metabolizer diplotype group indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (rapid metabolizer, the paper's reference group for CL_STD)",
      notes = "The paper dichotomises GSTA1 diplotypes into rapid versus slow metabolizers, where slow includes intermediate and poor metabolizers, following the Ben Hassine 2021 activity grouping (Methods; Table 1 footnote). Slow = 1 - GSTA1_RM, and CL is multiplied by FGSTA1_CL = 0.865 when GSTA1_RM = 0. The paper's reference stratum (rapid) is therefore GSTA1_RM = 1, the opposite of the register's default 0 reference; no coefficient is inverted.",
      source_name = "GSTA1 (rapid vs slow)"
    ),
    OCC = list(
      description = "Dosing occasion",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Integer occasion 1-4 multiplexing the between-occasion variability on CL. The paper reports BOV on CL (11.2%) without defining the occasion; sampling was around the first dose and around the fifth dose (Methods), so occasion = conditioning day (doses 1-4 = 1, doses 5-8 = 2, and so on, up to 4 for the 4-day regimen) is used here. Records with OCC outside 1-4 carry no BOV.",
      source_name = "OCC"
    )
  )

  # Screened in the stepwise covariate search (Methods; Supplementary Table S1)
  # and not retained: age, height and sex as separate covariates, the other GST
  # genes (GSTM1, GSTM5, GSTP1, GSTT1), and concomitant fludarabine, phenytoin
  # and metronidazole. The paper reports no estimates for them.

  population <- list(
    species = "human",
    n_subjects = 457,
    n_studies = 1,
    n_centers = 5,
    n_observations = 2103,
    age_range = "0.5-17.0 years (median 7.5; outcome cohort of 455)",
    weight_range = "5.0-96.5 kg (median 24.5; outcome cohort of 455)",
    sex_female_pct = 36.7,
    race_ethnicity = "Chinese",
    disease_state = "Pediatric patients (< 18 years) undergoing first allogeneic HSCT with IV busulfan-based myeloablative conditioning (BuCy, BuCyFlu or BuFlu); 73.9% malignant (AML, ALL, MDS, JMML, mixed-phenotype acute leukemia) and 26.1% nonmalignant disease (primary immunodeficiency, inherited metabolic disease, thalassemia, bone marrow failure)",
    dose_range = "0.8-1.2 mg/kg per dose (actual or adjusted body weight), 2 h IV infusion four times daily for 3 days (Bu3, n = 326) or 4 days (Bu4, n = 129), no dose adjustment during treatment",
    regions = "China (five pediatric HSCT centres; lead centre First Affiliated Hospital / Children's Hospital of Soochow University)",
    notes = "Prospective multicentre cohort (NCT04786002), January 2018 to January 2023; includes patients from two earlier single-centre studies (Du 2022; Yuan 2021). Samples at 2, 4 and 6 h after the start of the first infusion, with extra samples at 0.5, 1 and 2.5 h in some patients, and repeat sampling 30 min before the fifth dose or around the fifth dose in some patients. LLOQ 10, 80 or 100 ng/mL by centre; 3 of 2103 concentrations were below the LLOQ and were kept at their measured values. Demographics are from Table 2 (n = 455 outcome cohort; 2 of the 457 PK patients were excluded for follow-up < 3 months)."
  )

  ini({
    # Structural parameters standardised to a 70 kg, 176 cm adult male with
    # rapid GSTA1 metabolism (Table 1 and footnote).
    lcl <- log(10.8); label("Typical clearance for a 70 kg, 176 cm adult male (L/h)") # Table 1 CL_STD = 10.8 L/h
    lvc <- log(46.4); label("Typical volume of distribution for a 70 kg, 176 cm adult male (L)") # Table 1 V_STD = 46.4 L

    # Fractions of fat mass in the parameter-specific normal fat mass
    # (Supplementary Equation S2).
    ffat_cl <- 1.08; label("Fraction of fat mass contributing to NFM on CL (unitless)") # Table 1 Ffat_CL = 1.080
    ffat_v <- 0.653; label("Fraction of fat mass contributing to NFM on V (unitless)") # Table 1 Ffat_V = 0.653

    # Theory-based allometric exponents (Supplementary Equations S3 and S4).
    e_nfm_cl <- fixed(0.75); label("Allometric exponent of NFM on CL (unitless)") # Supplementary Equation S3 exponent 0.75
    e_nfm_vc <- fixed(1); label("Allometric exponent of NFM on V (unitless)") # Supplementary Equation S4 exponent 1

    # Sigmoid maturation of CL on postmenstrual age (Results).
    tm50_mat <- 41.6; label("Postmenstrual age at which CL reaches half the adult value (weeks)") # Table 1 TM50 = 41.6 weeks
    hill_mat <- 3.42; label("Hill coefficient of the CL maturation function (unitless)") # Table 1 HILL = 3.42

    # GSTA1 slow metabolizers (intermediate + poor) vs rapid (Equation S6).
    e_gsta1_slow_cl <- 0.865; label("Multiplicative CL factor for GSTA1 slow vs rapid metabolizers (unitless)") # Table 1 FGSTA1_CL = 0.865

    # Between-subject variability (exponential). Table 1 BSV_CL 0.336 and
    # BSV_V 0.252 (Results: 33.6% and 25.2%) are read as the SD of eta.
    etalcl ~ 0.112896 # Table 1 BSV_CL 0.336 -> 0.336^2
    etalvc ~ 0.063504 # Table 1 BSV_V 0.252 -> 0.252^2

    # Between-occasion variability on CL, one eta per occasion sharing a
    # variance (NONMEM SAME). Table 1 BOV_CL 0.112 (Results: 11.2%) read as SD.
    etaiov_cl_1 ~ 0.012544 # Table 1 BOV_CL 0.112 -> 0.112^2
    etaiov_cl_2 ~ fix(0.012544) # same variance as occasion 1
    etaiov_cl_3 ~ fix(0.012544) # same variance as occasion 1
    etaiov_cl_4 ~ fix(0.012544) # same variance as occasion 1

    # Combined residual error (Results: proportional 17.8% and additive
    # 0.0312 mg/L).
    propSd <- 0.178; label("Proportional residual error (fraction)") # Table 1 RUV_PROP = 0.178
    addSd <- 0.0312; label("Additive residual error (mg/L)") # Table 1 RUV_ADD = 0.0312 mg/L
  })

  model({
    # Parameter-specific normal fat mass (Supplementary Equation S2). The
    # standard NFM is the same formula for the 70 kg standard male with
    # FFMSTD = 56.1 kg (Table 1 footnote; Supplementary Equation S1).
    nfm_cl <- FFM + ffat_cl * (WT - FFM)
    nfm_v <- FFM + ffat_v * (WT - FFM)
    nfm_std_cl <- 56.1 + ffat_cl * (70 - 56.1)
    nfm_std_v <- 56.1 + ffat_v * (70 - 56.1)

    # Postmenstrual age (weeks) and sigmoid maturation of CL.
    pma <- AGE * 52 + GA
    fmat <- 1 / (1 + (pma / tm50_mat)^(-hill_mat))

    # GSTA1: slow metabolizers (GSTA1_RM = 0) clear 13.5% less.
    fgsta1 <- 1 + (e_gsta1_slow_cl - 1) * (1 - GSTA1_RM)

    # Between-occasion variability on CL.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4

    # Supplementary Equations S3 and S4.
    cl <- exp(lcl + etalcl + iov_cl) * (nfm_cl / nfm_std_cl)^e_nfm_cl * fmat * fgsta1
    vc <- exp(lvc + etalvc) * (nfm_v / nfm_std_v)^e_nfm_vc

    kel <- cl / vc
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
