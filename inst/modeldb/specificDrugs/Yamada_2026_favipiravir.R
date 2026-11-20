Yamada_2026_favipiravir <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral favipiravir, complete conversion to its hydroxylated metabolite M1 (one-compartment, disposition fixed from Wang 2020), and an indirect-response (turnover) model for serum uric acid in which the unbound-potency-scaled favipiravir and M1 concentrations additively inhibit uric acid elimination through a shared Imax (= 1) function, in Japanese adults treated for COVID-19 or severe fever with thrombocytopenia syndrome (Yamada 2026)."
  reference <- paste(
    "Yamada T, Kawasuji H, Ogami C, Hasegawa C, Kaneda M, Yamaguchi D, Iida S, Aoyama T, Yamamoto Y, Tsuji Y.",
    "Population Pharmacokinetic Analysis and Modelling of Serum Uric Acid Dynamics in Patients Treated with Favipiravir.",
    "Pharmaceuticals (Basel). 2026;19(7):1008. doi:10.3390/ph19071008.",
    "M1 disposition, the M1 renal-function covariate, and the favipiravir-to-M1 molecular-weight scaling are inherited from",
    "Wang Y, Zhong W, Salam A, Tarning J, et al. Phase 2a, open-label, dose-escalating, multi-center pharmacokinetic study of",
    "favipiravir (T-705) in combination with oseltamivir in patients with severe influenza. EBioMedicine. 2020;62:103125.",
    "doi:10.1016/j.ebiom.2020.103125 (Supplementary Appendix 1 NONMEM code)."
  )
  vignette <- "Yamada_2026_favipiravir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "favipiravir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "favipiravir", units = "mg", specimen = "serum", verified = TRUE),
    central_m1 = list(analyte = "favipiravir hydroxide (M1)", units = "mg", specimen = "serum", verified = TRUE),
    urate = list(analyte = "urate", units = "mg/dL", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric size factor FSIZE = (TBW/70)^PWR with PWR fixed to 0.75 for CL_FPV/F and CL_M1 and 1 for Vd_FPV and Vd_M1 (Yamada 2026 Equation 4 and Section 4.6). Also enters the Cockcroft-Gault creatinine clearance (Equation 5). Cohort median 68.0 kg (Table 1).",
      source_name = "TBW"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters only the Cockcroft-Gault creatinine clearance (Equation 5) that drives the M1 renal-function factor. Age was not retained as a covariate on CL_FPV/F (Section 2.2, dOFV = 0.004).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Enters only the Cockcroft-Gault creatinine clearance (x 0.85 if female; Equation 5).",
      source_name = "SEX"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Serum creatinine in mg/dL for the Cockcroft-Gault creatinine clearance CLcr = (140 - AGE) x TBW / (72 x SCr) (x 0.85 if female) (Equation 5), centred on CLcr_STD = 100 mL/min in the M1 renal-function factor RF = 1 + 0.0095 x (CLcr - 100) (Equation 7). Cohort median 0.9 mg/dL (Table 1).",
      source_name = "Serum creatinine"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 23L,
    n_studies = 1L,
    n_observations = "92 favipiravir serum concentrations (12 further BLQ samples excluded); 184 serum uric acid levels",
    age_range = "52-88 years (2.5th-97.5th percentile)",
    age_median = "69 years",
    weight_range = "40.2-87.9 kg (2.5th-97.5th percentile)",
    weight_median = "68.0 kg",
    sex_female_pct = 47.8,
    race_ethnicity = "Not reported (single Japanese centre)",
    disease_state = "Hospitalised patients treated with favipiravir for COVID-19 or severe fever with thrombocytopenia syndrome (SFTS)",
    dose_range = "Oral favipiravir 1800 mg twice daily on day 1, then 800 mg twice daily; median treatment 13 days (range 2-14)",
    regions = "Japan (University of Toyama Hospital, April 2020 - July 2025)",
    baseline_labs = "Median serum creatinine 0.9 mg/dL (0.5-1.5), albumin 3.7 g/dL, AST 46 IU/L, ALT 32 IU/L, CRP 3.12 mg/dL; baseline serum uric acid available in 18 of 23 patients (median 4.43 mg/dL)",
    notes = "Retrospective observational study using residual serum from routine care (Yamada 2026 Section 4.2, Table 1). M1 concentrations were not measured; M1 is predicted from the fixed Wang 2020 disposition."
  )

  ini({
    # Favipiravir PK (Yamada 2026 Equation 1 and Table 2, 70 kg reference)
    lka <- log(log(2) / 1.15); label("First-order absorption rate constant ka = ln(2)/Tabs (1/h)") # Table 2 'Tabs' = 1.15 h; Section 4.5 ka = ln(2)/Tabs
    lcl <- log(1.31); label("Apparent favipiravir clearance CL_FPV/F at 70 kg (L/h)") # Table 2 'CLFPV/F' = 1.31 L/h; Equation 1
    lvc <- log(10.4); label("Favipiravir volume of distribution Vd_FPV at 70 kg (L)") # Table 2 'VdFPV' = 10.4 L; Equation 1
    lfdepot <- fixed(log(1)); label("Oral bioavailability F (fraction)") # Table 2 'F' = 1 FIXED; Section 4.5

    # M1 (favipiravir hydroxide) disposition, fixed from Wang 2020
    lcl_m1 <- fixed(log(16.3)); label("M1 clearance CL_M1 at 70 kg and CLcr 100 mL/min (L/h)") # Table 2 'CLM1' = 16.30 FIXED; Wang 2020 Table 3
    lvc_m1 <- fixed(log(6.44)); label("M1 volume of distribution Vd_M1 at 70 kg (L)") # Table 2 'VdM1' = 6.44 FIXED; Wang 2020 Table 3
    e_crcl_cl_m1 <- fixed(0.0095); label("Linear creatinine-clearance coefficient TRF on CL_M1 (per mL/min)") # Table 2 'TRF' = 0.0095 FIXED; Equation 7

    # Allometric size factor FSIZE = (WT/70)^PWR
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL_FPV/F and CL_M1 (unitless)") # Section 4.6: PWR fixed to 0.75 for CLFPV/F and CLM1
    e_wt_vc <- fixed(1); label("Allometric exponent on Vd_FPV and Vd_M1 (unitless)") # Section 4.6: PWR fixed to 1 for VdFPV and VdM1

    # Uric acid turnover (Equation 8)
    lrbase <- log(4.43); label("Baseline serum uric acid UAbase (mg/dL)") # Table 2 'UAbase' = 4.43 mg/dL
    lkout <- log(log(2) / 17.6); label("Uric acid elimination rate constant kout = ln(2)/T1/2UA (1/h)") # Table 2 'T1/2UA' = 17.6 h; Section 4.7 kout = ln(2)/T1/2UA

    # Imax inhibition of uric acid elimination (Equations 9-10; additive, alpha = 0)
    limax <- fixed(log(1)); label("Maximum fractional inhibition of uric acid elimination Imax (unitless)") # Table 2 'Imax' = 1 FIXED
    lic50 <- fixed(log(126)); label("Unbound favipiravir concentration producing 50% of Imax (mg/L)") # Table 2 'C50FPV' = 126 mg/L FIXED (unbound)
    lic50_m1 <- fixed(log(52)); label("Unbound M1 concentration producing 50% of Imax (mg/L)") # Table 2 'C50M1' = 52 mg/L FIXED (unbound)
    fu <- fixed(0.5); label("Unbound fraction of favipiravir in serum (fraction)") # Table 2 'FUBFPV' = 0.5 FIXED; Equation 10
    fu_m1 <- fixed(0.3); label("Unbound fraction of M1 in serum (fraction)") # Table 2 'FUBM1' = 0.3 FIXED; Equation 10
    lhill <- log(0.41); label("Hill coefficient of the favipiravir inhibition term (unitless)") # Table 2 'GAMMAFPV' = 0.41
    lhill_m1 <- fixed(log(1)); label("Hill coefficient of the M1 inhibition term (unitless)") # Table 2 'GAMMAM1' = 1 FIXED

    # Between-subject variability. Table 2 footnote: BSV = sqrt(OMEGA), so each
    # printed value is an SD on the log scale and is squared here.
    etalcl ~ 0.233289 # Table 2 BSV 'CLFPV/F' = 0.483 -> 0.483^2
    etalcl_m1 ~ fixed(0.044944) # Table 2 BSV 'CLM1' = 0.212 (held constant) -> 0.212^2; Wang 2020 omega2 0.0448
    etalrbase ~ 0.069169 # Table 2 BSV 'UAbase' = 0.263 -> 0.263^2
    etalhill ~ 0.8464 # Table 2 BSV 'GAMMAFPV' = 0.920 -> 0.920^2

    # Residual error (Table 2 footnote: RUV estimated as THETA, i.e. SDs)
    propSd <- 0.445; label("Favipiravir proportional residual SD (fraction)") # Table 2 'RUVPROP_FPV' = 0.445
    addSd <- 0.135; label("Favipiravir additive residual SD (mg/L)") # Table 2 'RUVADD_FPV' = 0.135 mg/L
    addSd_m1 <- fixed(0.241); label("M1 additive residual SD (mg/L)") # Table 2 'RUVADD_M1' = 0.241 mg/L FIXED (M1 not measured)
    addSd_sUA <- 1.3; label("Serum uric acid additive residual SD (mg/dL)") # Table 2 'RUVADD_UA' = 1.300 mg/dL
  })

  model({
    # Cockcroft-Gault creatinine clearance with actual body weight (Equation 5;
    # Wang 2020 Appendix 1 uses the same weight-based form), in mL/min
    crcl <- (140 - AGE) * WT / (72 * CREAT) * (1 - 0.15 * SEXF)

    # M1 renal-function factor centred on CLcr_STD = 100 mL/min (Equations 6-7)
    rf <- 1 + e_crcl_cl_m1 * (crcl - 100)

    # Individual parameters (Equations 1, 2, 4, 6)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    cl_m1 <- exp(lcl_m1 + etalcl_m1) * (WT / 70)^e_wt_cl * rf
    vc_m1 <- exp(lvc_m1) * (WT / 70)^e_wt_vc

    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout)
    kin <- rbase * kout

    imax <- exp(limax)
    hill <- exp(lhill + etalhill)
    hill_m1 <- exp(lhill_m1)
    # Unbound potencies converted to total-concentration potencies (Equation 10)
    ic50_tot <- exp(lic50) / fu
    ic50_tot_m1 <- exp(lic50_m1) / fu_m1

    kel <- cl / vc
    kel_m1 <- cl_m1 / vc_m1

    # Favipiravir is converted completely to M1 (Section 4.5); the mass flux is
    # scaled by the M1/favipiravir molecular-weight ratio as in the Wang 2020
    # control stream (MOLM / MOLP = 173.10 / 157.10).
    mw_ratio_m1 <- 173.10 / 157.10

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(central_m1) <- mw_ratio_m1 * kel * central - kel_m1 * central_m1

    Cc <- central / vc
    Cc_m1 <- central_m1 / vc_m1

    # Additive Imax inhibition by favipiravir and M1 (Equation 9 with the
    # interaction coefficient alpha = 0 in the final model, Section 2.3)
    occ <- (Cc / ic50_tot)^hill
    occ_m1 <- (Cc_m1 / ic50_tot_m1)^hill_m1
    edrug <- imax * (occ + occ_m1) / (1 + occ + occ_m1)

    # Turnover with zero-order production and first-order elimination that the
    # drug inhibits (Equation 8; the printed 'kout x PDI' omits the UA state)
    d/dt(urate) <- kin - kout * (1 - edrug) * urate
    urate(0) <- rbase

    f(depot) <- exp(lfdepot)

    sUA <- urate

    Cc ~ add(addSd) + prop(propSd)
    Cc_m1 ~ add(addSd_m1)
    sUA ~ add(addSd_sUA)
  })
}
