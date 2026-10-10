Bell_2026_pirtobrutinib <- function() {
  description <- "Two-compartment population PK model with linear elimination and a four-compartment transit absorption chain for oral pirtobrutinib, a non-covalent Bruton tyrosine kinase inhibitor, in 595 adults with relapsed or refractory B-cell malignancies (mantle cell lymphoma, CLL/SLL and other non-Hodgkin lymphoma) from the phase 1/2 BRUIN study given 25-300 mg once daily (Bell 2026). Apparent clearance CL/F = 2.02 L/h, central volume Vc/F = 32.8 L, intercompartmental clearance Q/F = 8.38 L/h and peripheral volume Vp/F = 19.5 L for a 70 kg patient with eGFR 74.96 mL/min/1.73 m^2 and serum albumin 41.6 g/L; mean transit time MTT = 1.08 h. Body weight scales CL/F and Q/F (shared estimated exponent 0.524) and Vc/F and Vp/F (shared estimated exponent 0.785); CL/F rises exponentially with eGFR and falls with serum albumin (power -0.677), which also lowers Vc/F (power -0.513). Interindividual variability on CL/F (37.9% CV) and MTT (25.0% CV), inter-occasion variability on MTT (45.9% CV) and a proportional residual error (20.5%). None of the covariate effects was judged clinically meaningful and no dose adjustment is recommended."
  reference <- "Bell R, O'Brien LM, Yuen E, Liu D, Chapman SC. Population pharmacokinetic analysis of pirtobrutinib, a non-covalent BTK inhibitor, in patients with hematological malignancies from the Phase 1/2 BRUIN study. Cancer Chemother Pharmacol. 2026;96(1):111. doi:10.1007/s00280-026-04951-4. PMID: 42786239. PMCID: PMC13612640."
  vignette <- "Bell_2026_pirtobrutinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight at study entry",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power (allometric-form) scaling with ESTIMATED exponents, referenced to 70 kg: CL/F and Q/F share (WT/70)^0.524 and Vc/F and Vp/F share (WT/70)^0.785 (Bell 2026 Table 3 rows 'Allometry on CL and Q' and 'Allometry on Vc and Vp', footnotes a and b). Table 3 footnote defines WT as 'body weight at entry', i.e. baseline. Cohort median 77 kg (range 36-153; Table 2); 5th / 50th / 95th percentiles 51.8 / 76.6 / 113 kg (Results, 'Impact of covariates on pirtobrutinib PK'). Missing baseline weight (2 patients) was imputed with the population median (Methods, 'Handling of missing data').",
      source_name = "WT"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate, MDRD-6 (serum creatinine, age, sex, race, serum urea nitrogen and serum albumin), BSA-normalized",
      units = "mL/min/1.73m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Exponential effect on CL/F centred at 74.96 mL/min/1.73 m^2: CL/F = CL * exp(0.00329 * (eGFR - 74.96)) (Bell 2026 Table 3 row 'eGFR (mL/min/1.73 m2; Theta12)' and footnote c). The estimating equation is the six-variable MDRD Study equation, 170 * SCr^-0.999 * age^-0.176 * BUN^-0.17 * ALB(g/dL)^0.318 * 0.762 (female) * 1.18 (Black) (Table 1 footnote a and Table 2 footnote d). The 74.96 centring value is not the Table 2 cohort median (72; range 22-132) and its origin is not stated; it is used as printed. Supplementary Figure S2 uses 72.7 mL/min/1.73 m^2 as the median reference patient and 39.6 / 101 as the 5th / 95th percentiles.",
      source_name = "eGFR"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects referenced to 41.6 g/L: CL/F = CL * (ALB/41.6)^-0.677 and Vc/F = Vc * (ALB/41.6)^-0.513 (Bell 2026 Table 3 rows 'Albumin (g/L; Theta11)' and 'Albumin (g/L; Theta13)', footnotes d and e). Reported in SI g/L, the register's canonical unit, so no conversion. Cohort median 41 g/L (range 19-57; Table 2); 5th / 50th / 95th percentiles 31.0 / 41.0 / 47.8 g/L (Results, 'Impact of covariates on pirtobrutinib PK').",
      source_name = "ALB"
    ),
    OCC = list(
      description = "Occasion index for inter-occasion variability on the absorption mean transit time",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Values 1-4, decomposed inside model() into binary indicators oc1..oc4 that multiplex four IOV etas on log-MTT; occasions 2-4 share the occasion-1 variance (NONMEM $OMEGA BLOCK(1) SAME). Bell 2026 reports one MTT inter-occasion variability (Table 3, 45.9% CV) but does not define an occasion or state how many there were. Four occasions are therefore an extraction-side construct, matching the four BRUIN intensive PK sampling days (Cycle 1 Day 1, Cycle 1 Day 8, Cycle 2 Day 1, Cycle 4 Day 1; Methods, 'Pharmacokinetic sampling'). Assign each dosing occasion its own OCC; records outside 1-4 zero all indicators and drop IOV. Precedents for an unstated occasion count: Sasaki_2022_delamanid.R, Ding_2026_vancomycin.R.",
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on MTT, CL/F and V/F (Table 1); not retained. Results: 'Age, sex, race, ethnicity, cancer type, mild hepatic impairment, and formulation did not significantly affect pirtobrutinib disposition (Figure S1)'. Median 68 years (range 27-95; Table 2). No coefficient reported."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested on MTT, CL/F and V/F (Table 1); not retained. Supplementary Figure S1A shows empirical-Bayes CL/F box plots by sex only. 201 of 595 (34%) female."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Race tested on CL/F (Table 1); not retained. 17 of 595 (3%); Discussion notes the subgroup was too small to detect modest effects. No coefficient reported."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Race tested on CL/F (Table 1); not retained. 39 of 595 (7%); Discussion states the Asian subgroup was large enough to assess a clinically meaningful effect on CL/F and none was observed. No coefficient reported."
    ),
    RACE_HISPANIC = list(
      description = "Hispanic ethnicity indicator",
      units = "(binary)",
      type = "binary",
      notes = "Ethnicity tested on CL/F (Table 1); not retained. 23 of 595 (4%) Hispanic. No coefficient reported."
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic impairment indicator (NCI ODWG)",
      units = "(binary)",
      type = "binary",
      notes = "Hepatic function by NCI ODWG criteria tested on CL/F (Table 1); mild impairment not retained (106 of 595, 18%). Only 13 moderate and 1 severe patient were included, too few to draw conclusions (Results). No coefficient reported."
    ),
    TUMTP_MCL = list(
      description = "Mantle cell lymphoma indicator",
      units = "(binary)",
      type = "binary",
      notes = "Cancer type (MCL 23%, CLL/SLL 44%, other NHL 32%; Table 2) tested on MTT, CL/F and V/F (Table 1); not retained. No coefficient reported."
    ),
    TUMTP_CLL = list(
      description = "Chronic lymphocytic leukemia / small lymphocytic lymphoma indicator",
      units = "(binary)",
      type = "binary",
      notes = "See TUMTP_MCL. Not retained; no coefficient reported."
    ),
    FORM_PIRTOBRUTINIB_T2 = list(
      description = "Commercial tablet formulation (T2) indicator, reference = initial tablet T1",
      units = "(binary)",
      type = "binary",
      notes = "Formulation tested on F and MTT (Table 1); not retained. 210 patients (35%) started on the initial tablet T1 and 385 (65%) on the commercial tablet T2 (Methods, 'Study design and patients'). No coefficient reported."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "pirtobrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "pirtobrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "pirtobrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "pirtobrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "pirtobrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pirtobrutinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pirtobrutinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 595L,
    n_studies = 1L,
    n_observations = 4487L,
    age_range = "27-95 years",
    age_median = "68 years",
    weight_range = "36-153 kg",
    weight_median = "77 kg",
    bmi_median = "26 kg/m^2 (range 14-47)",
    sex_female_pct = 34,
    race_ethnicity = c(White = 86, BlackOrAfricanAmerican = 3, Asian = 7, Other = 5, NotReported = 0.2),
    ethnicity = c(NonHispanic = 92, Hispanic = 4, NotReported = 4),
    disease_state = "Relapsed or refractory B-cell malignancies after failure of or intolerance to standard therapy: mantle cell lymphoma 23%, CLL/SLL 44%, other non-Hodgkin lymphoma 32%",
    dose_range = "Pirtobrutinib monotherapy 25, 50, 100, 150, 200, 250 or 300 mg orally once daily (phase 1 dose escalation); 200 mg once daily (recommended phase 2 dose) in phase 1 expansion and phase 2; 28-day cycles",
    renal_function = "eGFR (MDRD-6) median 72 mL/min/1.73 m^2 (range 22-132); normal 20%, mild impairment 51%, moderate 28%, severe 1%",
    hepatic_function = "NCI ODWG: normal 80%, mild impairment 18%, moderate 2%, severe < 1%",
    albumin = "Serum albumin median 41 g/L (range 19-57)",
    formulation = "Initial tablet T1 35%, commercial tablet T2 65% (formulation of first dose)",
    notes = "Baseline characteristics from Bell 2026 Table 2. BRUIN (NCT03740529), open-label multicentre phase 1/2 trial; interim data cutoff 31 January 2022. Intensive sampling (predose, 1, 2, 4 and 8 h) on Cycle 1 Day 1, Cycle 1 Day 8, Cycle 2 Day 1 and Cycle 4 Day 1 in phase 1 dose escalation; single predose samples on Cycle 1 Day 8 and Cycle 4 Day 1 in phase 1 expansion and phase 2. 4487 of 4867 concentrations from 595 of 611 patients were evaluable. Strong and moderate CYP3A4 inhibitors and inducers were prohibited. NONMEM 7.4.2, FOCE with interaction."
  )

  ini({
    # ---- Absorption (Bell 2026 Table 3) ----------------------------------
    lfdepot <- fixed(log(1))
    label("Oral bioavailability F (unitless)")
    # Table 3 row 'Bioavailability (F, fraction, Theta1)' = 1 fixed

    lmtt <- log(1.08)
    label("Mean absorption transit time MTT (h)")
    # Table 3 row 'MTT (h, Theta2)' = 1.08 (%SEE 2.97)

    # ---- Disposition (Bell 2026 Table 3) ---------------------------------
    # Typical values for WT = 70 kg, eGFR = 74.96 mL/min/1.73 m^2 and
    # ALB = 41.6 g/L (Table 3 footnotes a-e).
    lcl <- log(2.02)
    label("Apparent clearance CL/F (L/h)")
    # Table 3 row 'Clearance (CL, L/h, Theta3)' = 2.02 (%SEE 1.66)

    lvc <- log(32.8)
    label("Apparent central volume Vc/F (L)")
    # Table 3 row 'Central volume of distribution (Vc, L, Theta4)' = 32.8 (%SEE 3.60)

    lq <- log(8.38)
    label("Apparent intercompartmental clearance Q/F (L/h)")
    # Table 3 row 'Intercompartmental clearance (Q, L/h, Theta5)' = 8.38 (%SEE 10.5)

    lvp <- log(19.5)
    label("Apparent peripheral volume Vp/F (L)")
    # Table 3 row 'Peripheral volume of distribution (Vp, L, Theta6)' = 19.5 (%SEE 5.49)

    # ---- Covariate effects (Bell 2026 Table 3) ---------------------------
    e_wt_cl_q <- 0.524
    label("Shared body-weight exponent on CL/F and Q/F (unitless)")
    # Table 3 'Allometry on CL and Q / Body Weight (kg; Theta9)' = 0.524 (%SEE 11.5); footnote a

    e_wt_vc_vp <- 0.785
    label("Shared body-weight exponent on Vc/F and Vp/F (unitless)")
    # Table 3 'Allometry on Vc and Vp / Body Weight (kg; Theta10)' = 0.785 (%SEE 6.31); footnote b

    e_crcl_cl <- 0.00329
    label("Exponential eGFR coefficient on CL/F (per mL/min/1.73 m^2)")
    # Table 3 'Covariate effects on CL / eGFR (mL/min/1.73 m2; Theta12)' = 0.00329 (%SEE 31.0); footnote c

    e_alb_cl <- -0.677
    label("Power exponent of serum albumin on CL/F (unitless)")
    # Table 3 'Covariate effects on CL / Albumin (g/L; Theta11)' = -0.677 (%SEE 16.7); footnote d

    e_alb_vc <- -0.513
    label("Power exponent of serum albumin on Vc/F (unitless)")
    # Table 3 'Covariate effect on Vc / Albumin (g/L; Theta13)' = -0.513 (%SEE 22.0); footnote e

    # ---- Interindividual variability (Bell 2026 Table 3) -----------------
    # Table 3 reports IIV and IOV as CV%; each log-normal variance is
    # recovered as omega^2 = log(1 + CV^2).
    etalmtt ~ 0.060625
    # Table 3 IIV 'MTT (Omega2)' = 25.0% CV (%SEE 30.7); log(1 + 0.250^2) = 0.060625
    etalcl ~ 0.134217
    # Table 3 IIV 'CL (Omega3)' = 37.9% CV (%SEE 7.61); log(1 + 0.379^2) = 0.134217

    # ---- Inter-occasion variability on MTT (Bell 2026 Table 3) -----------
    # One magnitude is reported; occasions 2-4 share the occasion-1 variance
    # (NONMEM $OMEGA BLOCK(1) SAME). The occasion count is not stated by the
    # paper; see covariateData[['OCC']].
    etaiov_mtt_1 ~ 0.191183
    # Table 3 IOV 'MTT' = 45.9% CV (%SEE 13.4); log(1 + 0.459^2) = 0.191183
    etaiov_mtt_2 ~ fixed(0.191183)
    etaiov_mtt_3 ~ fixed(0.191183)
    etaiov_mtt_4 ~ fixed(0.191183)

    # ---- Residual variability (Bell 2026 Table 3) ------------------------
    propSd <- 0.205
    label("Proportional residual error (fraction)")
    # Table 3 'Residual variability / Proportional' = 0.205 (%SEE 2.37); read as an SD
  })

  model({
    # --- 1. Occasion indicators and absorption ----------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 + oc4 * etaiov_mtt_4

    mtt <- exp(lmtt + etalmtt + iov_mtt)
    # Savic parameterisation of the four-transit-compartment chain: the dose
    # enters depot and passes depot -> transit1 -> ... -> transit4 -> central,
    # five first-order steps at the common rate ktr, so MTT = (4 + 1) / ktr is
    # the mean time from dosing to arrival in central.
    ktr <- (4 + 1) / mtt

    # --- 2. Individual disposition parameters (Table 3 footnotes a-e) -----
    cl <- exp(lcl + etalcl) *
      (WT / 70)^e_wt_cl_q *
      exp(e_crcl_cl * (CRCL - 74.96)) *
      (ALB / 41.6)^e_alb_cl
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc) * (WT / 70)^e_wt_vc_vp * (ALB / 41.6)^e_alb_vc
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp

    # --- 3. Micro-constants (Methods half-life equation) -------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # --- 4. ODE system -----------------------------------------------------
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(central) <- ktr * transit4 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # --- 5. Bioavailability -------------------------------------------------
    f(depot) <- exp(lfdepot)

    # --- 6. Observation and error -----------------------------------------
    # Dose in mg and volume in L give mg/L = ug/mL; x 1000 for ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
