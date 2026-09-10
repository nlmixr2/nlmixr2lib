Kang_2025_radotinib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral radotinib (a second-",
    "generation BCR-ABL1 tyrosine kinase inhibitor) in Asian healthy",
    "volunteers and patients with chronic-phase chronic myeloid leukemia",
    "(Kang 2025). Savic transit-compartment absorption (analytical form,",
    "N = 6.58, MTT = 1.88 h) feeds a first-order absorption step (ka) into",
    "a two-compartment disposition model with first-order elimination.",
    "Apparent oral clearance carries a 24 h cosinor circadian rhythm with a",
    "fractional amplitude of 0.683 peaking 7 h after the reference morning",
    "dose (about 16:00 clock time), and is 64.6% faster in healthy",
    "volunteers than in CML patients (37.9 versus 23.0 L/h). Apparent",
    "central volume declines linearly with age at -1.29% of the age-31",
    "typical value per year. Interindividual variability is carried on CL/F,",
    "MTT and N; interoccasion variability is shared between the central and",
    "peripheral volumes across the day-1 and day-14 sampling occasions.",
    sep = " "
  )
  reference <- paste(
    "Kang M, Kim J, Lee Y, Shin JS, Park MS, Jiang Q, Chung EK, Lee JI.",
    "Population Pharmacokinetics of Radotinib in Healthy Volunteers and",
    "Patients with Chronic Myeloid Leukemia. Pharmaceuticals (Basel).",
    "2025 Nov 10;18(11):1705. doi:10.3390/ph18111705. PMCID: PMC12655659.",
    sep = " "
  )
  vignette <- "Kang_2025_radotinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(
      analyte = "radotinib", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "radotinib", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "radotinib", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    DIS_CML = list(
      description        = "Chronic-phase chronic myeloid leukemia disease-state indicator: 1 = patient with CML-CP, 0 = healthy volunteer",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (patient with CML-CP) -- see notes; the published typical CL/F of 23.0 L/h is the CML-CP value, and the estimated coefficient is the healthy-volunteer increment",
      notes              = paste(
        "Kang 2025 equation 1 and Table 2 footnote a: CL/F (L/h) = 23.0 *",
        "[1 + 0.646 * (1 - CML)] * [1 + 0.683 * (circadian effect)], with",
        "'disease status: 0 = healthy volunteers, 1 = patients with chronic",
        "myeloid leukemia'. The coefficient therefore multiplies the",
        "COMPLEMENT (1 - DIS_CML): CL/F is 23.0 L/h in a CML-CP patient and",
        "23.0 * 1.646 = 37.9 L/h in a healthy volunteer, i.e. 39.2% slower",
        "in patients (Abstract, Results section 2.3). The same",
        "apply-via-the-complement pattern is recorded for DIS_AML under the",
        "Vaddady 2024 quizartinib alias. Kang 2025 Discussion (study",
        "limitations, fourth point) cautions that the disease effect is",
        "confounded with ethnicity and sex in this pooled dataset: the",
        "healthy volunteers were all Korean males (n = 23) and all CML-CP",
        "patients were Chinese of both sexes (n = 24), so the estimate may",
        "partly reflect interethnic or sex-related variability. Time-fixed",
        "per subject.",
        sep = " "
      ),
      source_name        = "CML"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Kang 2025 equation 2 and Table 2 footnote b: Vc/F (L) = 383 *",
        "[1 - 0.0129 * (age - 31)]. Centred at 31 years, the pooled cohort",
        "median (Table 1: median 31, range 20-72 years). Results section",
        "2.3 reports the realised Vc/F spread over that range as 180 to",
        "437 L, which the centred linear form reproduces exactly at ages 72",
        "and 20. NOTE the linear form is only valid over roughly the",
        "observed age range: it reaches zero at age 108.5 years and turns",
        "negative beyond it, so do not extrapolate. Kang 2025 Discussion",
        "notes the effect is statistically significant but of limited",
        "clinical magnitude within 20-72 years because steady-state",
        "exposure is governed by clearance, not by Vc/F.",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    OCC = list(
      description        = "Sampling-occasion indicator: 1 = day-1 (non-steady-state) occasion, 2 = day-14 (steady-state) occasion",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Kang 2025 Methods section 4.3: 'Interoccasion variability (IOV)",
        "was incorporated in the model assuming a log-normal distribution",
        "..., where each occasion was defined as a distinct pharmacokinetic",
        "state (i.e., Day 1 as non-steady state, Day 14 as steady state)'.",
        "The CML-CP study (NCT03722420) sampled on both days; the",
        "single-dose healthy-volunteer study (NCT06461078) contributes only",
        "the day-1 occasion, so healthy-volunteer records take OCC = 1.",
        "Decomposed inside model() into binary indicators oc1 and oc2 that",
        "multiplex the two per-occasion IOV etas, whose variances are equal",
        "across occasions per the source's single-variance IOV reporting",
        "(Table 2 reports one 'omega IOV' row).",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 47L,
    n_studies      = 2L,
    n_observations = 640L,
    age_range      = "20-72 years (Table 1: pooled median 31; healthy volunteers median 29, range 20-51; CML-CP patients median 32, range 21-72)",
    age_median     = "31 years",
    weight_range   = "47.0-96.0 kg (Table 1: pooled median 65.7 kg)",
    weight_median  = "65.7 kg",
    height_range   = "155-186 cm (Table 1: pooled median 171 cm)",
    bmi_range      = "18.4-30.3 kg/m^2 (Table 1: pooled median 22.8 kg/m^2)",
    sex_female_pct = 23.4,
    race_ethnicity = c(Korean = 49, Chinese = 51),
    hepatic_function = "ALT 7-80 IU/L (pooled median 20); AST 11-76 IU/L (pooled median 20) (Table 1)",
    renal_function = "Cockcroft-Gault creatinine clearance 72-200 mL/min (pooled median 114 mL/min) (Table 1)",
    disease_state  = paste(
      "23 healthy adult male volunteers (Korean) and 24 patients with",
      "Philadelphia-chromosome-positive chronic-phase chronic myeloid",
      "leukemia (Chinese) diagnosed within the preceding 6 months, with",
      "typical BCR-ABL1 transcript type and normal organ function;",
      "T315I-mutation carriers were excluded."
    ),
    dose_range     = paste(
      "Healthy volunteers: single oral 400 mg radotinib under fasting",
      "conditions. CML-CP patients: 300 mg orally twice daily under fasting",
      "conditions, sampled on days 1 and 14."
    ),
    regions        = "Republic of Korea (healthy volunteers) and China (CML-CP patients); all 47 participants of Asian descent",
    notes          = paste(
      "Pooled analysis of two clinical PK studies, ClinicalTrials.gov",
      "NCT06461078 (healthy volunteers, 306 concentrations) and",
      "NCT03722420 (CML-CP patients, 334 concentrations); 640",
      "concentrations in total (Results section 2.1, Table 1). NONMEM",
      "7.5.1, FOCE with interaction; final OFV 6422. Bioanalytical range",
      "5-2000 ng/mL by HPLC-MS/MS in both studies (Methods section 4.2).",
      "Kang 2025 Discussion flags the small, exclusively Asian cohort as a",
      "limitation and notes that concomitant CYP3A4 modulators were not",
      "tested as independent covariates but are implicitly absorbed into",
      "the disease-status term."
    )
  )

  ini({
    # =====================================================================
    # Structural PK parameters. All values are FINAL population estimates
    # from Kang 2025 Table 2 ('Population Estimates' column) and equations
    # 1-7 of Results section 2.3. Bootstrap medians and 95% CIs from the
    # 1000-replicate nonparametric bootstrap are quoted alongside for
    # reference. The typical CL/F of 23.0 L/h is the value for a CML-CP
    # patient at the circadian mesor (cos term = 0); the healthy-volunteer
    # value 23.0 * 1.646 = 37.9 L/h is recovered through e_dis_cml_cl.
    # =====================================================================
    lcl <- log(23.0)
    label("Typical apparent oral clearance CL/F in a CML-CP patient at the circadian mesor (L/h)")
    # Table 2 row 'CL/F (L/h)' = 23.0 (RSE 7.3%); bootstrap 23.2 (20.0-26.5)

    lvc <- log(383)
    label("Typical apparent central volume of distribution Vc/F at the median age of 31 years (L)")
    # Table 2 row 'Vc/F (L)' = 383 (RSE 13.8%); bootstrap 395 (281-526)

    lq <- log(132)
    label("Apparent intercompartmental clearance Q/F (L/h)")
    # Table 2 row 'Q/F (L/h)' = 132 (RSE 11.6%); bootstrap 126 (88.9-185); equation 3 carries no eta

    lvp <- log(519)
    label("Typical apparent peripheral volume of distribution Vp/F (L)")
    # Table 2 row 'Vp/F (L)' = 519 (RSE 15.6%); bootstrap 529 (398-733)

    lka <- log(1.59)
    label("First-order absorption rate constant ka from the depot to the central compartment (1/h)")
    # Table 2 row 'ka (h-1)' = 1.59 (RSE 21.6%); bootstrap 1.56 (1.07-3.14); equation 5 carries no eta

    lmtt <- log(1.88)
    label("Mean transit time MTT of the Savic transit-absorption chain (h)")
    # Table 2 row 'MTT (h)' = 1.88 (RSE 10.4%); bootstrap 1.86 (1.43-2.27)

    lnn <- log(6.58)
    label("Number of transit compartments N of the Savic chain (non-integer allowed) (count)")
    # Table 2 row 'N' = 6.58 (RSE 17.0%); bootstrap 6.34 (4.46-10.6)

    # =====================================================================
    # Circadian rhythm on CL/F. Kang 2025 equation 1 writes the modulation
    # as the bracket [1 + 0.683 * cos{2*pi*(TIME - 7)/24}], so the amplitude
    # is FRACTIONAL (a multiplier of the mesor), unlike the absolute-units
    # cosinor amplitude of MohammedAli_2025_tacrolimus. CL/F therefore
    # ranges over 23.0 * (1 +/- 0.683) = 7.29 to 38.7 L/h across the day in
    # a CML-CP patient; the bracket never changes sign because 0.683 < 1.
    #
    # The phase shift of 7 h appears ONLY inside equation 1 -- Table 2 does
    # not tabulate it and reports no RSE for it -- so it is encoded as a
    # structural constant. TIME is the NONMEM time variable, i.e. hours
    # since the first dose, NOT clock time: Kang 2025 Discussion states the
    # model puts peak CL/F activity at 16:00, which is 7 h after a ~09:00
    # morning dose. Users simulating a different dosing clock time should
    # shift acrophase_cl accordingly.
    # =====================================================================
    amp_cl <- 0.683
    label("Fractional amplitude of the 24 h cosinor circadian rhythm on CL/F (fraction of the mesor)")
    # Table 2 row 'CL/F Circadian effect' = 0.683 (RSE 14.8%); bootstrap 0.686 (0.502-0.871)

    acrophase_cl <- fixed(7)
    label("Acrophase of the CL/F circadian rhythm, in hours after the first dose (h)")
    # Kang 2025 equation 1: cos{2*pi*(TIME - 7)/24}; not tabulated in Table 2, no RSE reported

    # =====================================================================
    # Covariate effects.
    # =====================================================================
    e_dis_cml_cl <- 0.646
    label("Fractional increase in CL/F for a healthy volunteer (DIS_CML = 0) relative to a CML-CP patient (unitless)")
    # Table 2 row 'CL/F Disease status' = 0.646 (RSE 30.0%); bootstrap 0.625 (0.249-1.02); enters as (1 - DIS_CML) per equation 1

    e_age_vc <- -0.0129
    label("Fractional change in Vc/F per year of age above the median of 31 years (1/year)")
    # Table 2 row 'Vc/F Age' = -0.0129 (RSE 37.6%); bootstrap -0.0127 (-0.0258 to -0.003)

    # =====================================================================
    # Interindividual variability. Kang 2025 Table 2 reports the IIV block
    # as 'omega <parameter>' values that are STANDARD DEVIATIONS, not
    # variances. The paper settles this itself: Results section 2.3 states
    # 'the variances for Vc/F, Q, Vp/F, ka, MTT and N were fixed at 0, 0, 0,
    # 0, 0.1 and 0.2, respectively', and Table 2 prints omega MTT = 0.316 =
    # sqrt(0.1) and omega N = 0.447 = sqrt(0.2). Every entry below is
    # therefore squared to reach the log-scale variance rxode2 expects.
    #
    # IIV on Vc/F, Q/F, Vp/F and ka was fixed at exactly zero, so those
    # parameters carry no eta at all -- a zero-variance eta would make the
    # OMEGA matrix singular and break rxSolve's Cholesky sampler.
    # =====================================================================
    etalcl ~ 0.151321
    # Table 2 row 'omega CL/F' = 0.389 (RSE 20.1%), a SD -> variance 0.389^2 = 0.151321; bootstrap 0.379 (0.301-0.459)
    etalmtt ~ fixed(0.1)
    # Results section 2.3: MTT IIV variance FIXED at 0.1; Table 2 row 'omega MTT' = 0.316 = sqrt(0.1), RSE 'NA'
    etalnn ~ fixed(0.2)
    # Results section 2.3: N IIV variance FIXED at 0.2; Table 2 row 'omega N' = 0.447 = sqrt(0.2), RSE 'NA'

    # =====================================================================
    # Interoccasion variability. Results section 2.3: 'IOV was estimated for
    # volumes of distribution (i.e., Vc/F and Vp/F)'. Equations 2 and 4 both
    # multiply by the SAME symbol e^(eta_IOV) with no distinguishing
    # subscript, and Table 2 reports a single 'omega IOV' row, so one
    # occasion-level random effect is shared by both volumes. With the two
    # occasions of the study design (day 1, day 14) this is encoded as two
    # per-occasion etas of equal variance, the second fixed to the first,
    # matching the NONMEM '$OMEGA BLOCK(1) SAME' idiom used by the sibling
    # OCC models in this library.
    # =====================================================================
    etaiov_v_1 ~ 0.487204
    # Table 2 row 'omega IOV' = 0.698 (RSE 25.1%), a SD -> variance 0.698^2 = 0.487204; bootstrap 0.669 (0.494-0.829) (occasion 1)
    etaiov_v_2 ~ fixed(0.487204)
    # IOV on both log-volumes, occasion 2; variance equal to occasion 1 per the source's single-variance IOV reporting

    # =====================================================================
    # Residual variability. Results section 2.3: 'Residual variability was
    # most appropriately described using a proportional model.'
    # =====================================================================
    propSd <- 0.200
    label("Proportional residual error (fraction)")
    # Table 2 row 'sigma proportional (%)' = 20.0, RSE 'NA' -> 0.200 as a fraction
  })

  model({
    # 1. Occasion indicators (binary decomposition of the OCC column).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    # 2. Occasion-level IOV eta shared by Vc/F and Vp/F (equations 2 and 4).
    iov_v <- oc1 * etaiov_v_1 + oc2 * etaiov_v_2

    # 3. Circadian modulation of CL/F (equation 1). `time` is the rxode2
    #    solve time, which plays the role of NONMEM's TIME: hours since the
    #    first dose. See ini() for why the acrophase is 7 h. The cosine term
    #    is named after its amplitude parameter, following the circadian
    #    precedent in MohammedAli_2025_tacrolimus: the `cl_time_` / `cl_exp_`
    #    canonical stems are reserved for SECULAR (monotone) time-dependent
    #    clearance, and parameter-names.md states that periodic (diurnal /
    #    circadian) variation keeps its own names.
    amp_cl_t <- amp_cl * cos(2 * pi * (time - acrophase_cl) / 24)

    # 4. Individual structural PK parameters.
    #    CL/F: equation 1. The disease coefficient multiplies (1 - DIS_CML)
    #    so that the typical value 23.0 L/h applies to a CML-CP patient and
    #    a healthy volunteer gets 23.0 * 1.646 = 37.9 L/h.
    cl  <- exp(lcl + etalcl) * (1 + e_dis_cml_cl * (1 - DIS_CML)) * (1 + amp_cl_t)
    #    Vc/F: equation 2, linear in age centred at the cohort median 31 y.
    vc  <- exp(lvc + iov_v) * (1 + e_age_vc * (AGE - 31))
    #    Q/F: equation 3, no covariate and no random effect.
    q   <- exp(lq)
    #    Vp/F: equation 4, IOV only.
    vp  <- exp(lvp + iov_v)
    #    ka: equation 5, no covariate and no random effect.
    ka  <- exp(lka)
    #    MTT and N: equations 6 and 7, each with exponential IIV.
    mtt <- exp(lmtt + etalmtt)
    nn  <- exp(lnn + etalnn)
    #    Savic 2007 transit-chain rate: MTT = (N + 1)/ktr, the "+1"
    #    accounting for the depot compartment that the chain feeds.
    ktr <- (nn + 1) / mtt

    # 5. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 6. ODE system. The transit chain is implemented in the analytical
    #    Savic form so that the non-integer N = 6.58 is handled smoothly:
    #      input(t) = D * ktr * (ktr*t)^N * exp(-ktr*t) / Gamma(N + 1)
    #    evaluated in log space via lgammafn(). `tad(depot)` is the time
    #    since the most recent dose and is called once and stored, because
    #    rxode2 fails to parse two time-function calls in one expression;
    #    it is floored away from zero so that log(0) cannot arise. The dose
    #    lands on depot but its bolus is suppressed with f(depot) <- 0, so
    #    the whole dose enters through the transit input rate instead of
    #    being delivered twice.
    tdos  <- tad(depot)
    tdose <- max(tdos, 1e-8)
    input <- exp(
      log(podo(depot)) + log(ktr) + nn * log(ktr * tdose) -
        ktr * tdose - lgammafn(nn + 1)
    )
    d/dt(depot)       <- input - ka * depot
    d/dt(central)     <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                              k12 * central - k21 * peripheral1
    f(depot) <- 0

    # 7. Observation and error. Doses are in mg and volumes in L, so
    #    central/vc is mg/L; radotinib concentrations are reported in ng/mL
    #    throughout Kang 2025 (Tables 2 and 3, Figures 2 and 3), and
    #    1 mg/L = 1000 ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
