Gibiansky_2022_tisotumab <- function() {
  description <- "Joint population PK model for the antibody-drug conjugate tisotumab vedotin (ADC) and its released payload monomethyl auristatin E (MMAE) in adults with locally advanced or metastatic solid tumors (Gibiansky 2022). The ADC is a two-compartment model with parallel linear and Michaelis-Menten elimination from the central compartment. Linear ADC elimination, multiplied by a drug-to-antibody ratio that decays mono-exponentially with time after dose from 4 to 1, feeds an amount-only MMAE delay compartment in full and the central MMAE compartment by an additional fraction FR1; a fraction FR2 of the Michaelis-Menten flux also feeds the delay compartment. The delay compartment drains first-order into a one-compartment MMAE model with apparent clearance and volume. Covariates: body weight, albumin and sex on ADC CL and Vc; weight on Q and Vp; weight, albumin, tumor type, eGFR, tumor size, ECOG performance status and hepatic impairment on MMAE clearance; weight, ECOG and albumin on MMAE volume; age and weight on the delay rate constant."
  reference <- paste(
    "Gibiansky L, Passey C, Voellinger J, Gunawan R, Hanley WD, Gupta M, Winter H.",
    "Population pharmacokinetic analysis for tisotumab vedotin in patients with locally advanced",
    "and/or metastatic solid tumors. CPT Pharmacometrics Syst Pharmacol. 2022;11(10):1358-1370.",
    "doi:10.1002/psp4.12850. PMID 35932175; PMCID PMC9574719.",
    "Covariate coefficients and variance parameters, which the article does not print, are from",
    "the same final model as tabulated in the FDA BLA 761208 Multi-discipline Review (2021),",
    "Tables 49 and 50."
  )
  vignette <- "Gibiansky_2022_tisotumab"
  paper_specific_etas <- c("etaRUV", "etaRUV_mmae")

  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects normalised to 75 kg (Gibiansky 2022 Table S1) on ADC CL and Q (shared exponent), ADC Vc, ADC Vp, MMAE CL, MMAE V and the delay rate constant ktr. Analysis-population median 70.1 kg, range 33-148 kg (FDA BLA 761208 review Table 45).",
      source_name = "WT"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects normalised to 40 g/L (Gibiansky 2022 Table S1) on ADC CL, ADC Vc, MMAE CL and MMAE V. Analysis-population median 40 g/L, range 27-52 g/L (FDA BLA 761208 review Table 45).",
      source_name = "ALBUM"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = 1,
      notes = "The source codes a male indicator (Table S1 '[theta22 if male]'), so the reference is female. The male multipliers 1.09 (CL) and 1.13 (Vc) are kept as published and applied as e_sexm_*^(1 - SEXF).",
      source_name = "SEX"
    ),
    TUMTP_OTHER = list(
      description = "Non-cervical tumor type indicator (1 = any tumor type other than cervical cancer, 0 = cervical cancer)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Gibiansky 2022 Table S1 '[theta27 if other tumors]' on MMAE clearance; reference is cervical cancer (172 / 399 patients, Table 1). The 'other' pool is the remaining 227 patients of innovaTV 201, 202 and 207 (solid tumors known to express tissue factor).",
      source_name = "other tumors"
    ),
    TUMSZ = list(
      description = "Baseline tumor size as the sum of diameters of target lesions",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Exponential effect exp((TUMSZ / 60 - 1) * theta28) on MMAE clearance (Gibiansky 2022 Table S1, source column SUMDIAM). Analysis-population median 57 mm, range 0-299 mm; 8 (2.0%) values missing in the source data set (FDA BLA 761208 review Table 45).",
      source_name = "SUMDIAM"
    ),
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate by the CKD-EPI equation",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect normalised to 90 mL/min/1.73 m^2 on MMAE clearance (Gibiansky 2022 Table S1, source column CKDEPI). Analysis-population median 91, range 24.1-145 mL/min/1.73 m^2 (FDA BLA 761208 review Table 45).",
      source_name = "CKDEPI"
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status >= 1 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Gibiansky 2022 Table S1 '[theta30 if ECOG > 0]' on MMAE clearance and '[theta32 if ECOG > 0]' on MMAE volume. All analysed patients had ECOG 0 (38.8%) or 1 (61.2%) (Table 1).",
      source_name = "ECOG"
    ),
    HEPIMP = list(
      description = "Hepatic impairment indicator by the NCI Organ Dysfunction Working Group criteria (1 = impaired, 0 = normal)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Gibiansky 2022 Table S1 '[theta31 if hepatic impairment]' on MMAE clearance. Every impaired patient in the analysis had mild impairment (58 / 399, Table 1); the model has no information on moderate or severe impairment.",
      source_name = "hepatic impairment"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect normalised to 60 years on the delay rate constant ktr (Gibiansky 2022 Table S1). Analysis-population mean 56.1 years, range 21-81 years.",
      source_name = "AGE"
    ),
    STUDY_INNOVATV202 = list(
      description = "innovaTV 202 study indicator (NCT02552121; 1 = record from innovaTV 202, 0 = other studies)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Multiplies the ADC residual-error SD by 1.42 (Gibiansky 2022 Table S1 '[theta36 if innovaTV 202 study]'; value from FDA BLA 761208 review Table 49). Affects the residual error only, not the PK. Set to 0 for simulation of a new population.",
      source_name = "innovaTV 202 study"
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "tisotumab vedotin (antibody-drug conjugate)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "tisotumab vedotin (antibody-drug conjugate)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    precursor_mmae = list(
      analyte = "monomethyl auristatin E (MMAE), released and not yet in plasma; amount in ADC-mass-equivalent units scaled by the drug-to-antibody ratio",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    central_mmae = list(
      analyte = "monomethyl auristatin E (MMAE); amount in ADC-mass-equivalent units scaled by the drug-to-antibody ratio",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 399L,
    n_studies = 4L,
    age_range = "21-81 years",
    age_mean = "56.1 years",
    weight_range = "33-148 kg",
    weight_median = "70.1 kg",
    sex_female_pct = 74.2,
    race_ethnicity = c(White = 92.2, Black = 1.5, Asian = 2.5, Other = 3.8),
    disease_state = "Locally advanced and/or metastatic solid tumors known to express tissue factor (cervical cancer 43.1%; other tumors 56.9%).",
    dose_range = "0.3-2.2 mg/kg IV every 3 weeks (innovaTV 201); 1.2 mg/kg on days 1, 8 and 15 of a 28-day cycle or 2.0 mg/kg every 3 weeks (innovaTV 202); 2.0 mg/kg every 3 weeks, at most 200 mg, in innovaTV 204 and 207.",
    regions = "Europe (70.7%), United States (29.3%)",
    renal_function = "Normal 53.9%, mild impairment 35.6%, moderate impairment 10.5%; no severe impairment.",
    hepatic_function = "Normal 85.5%, mild impairment (NCI ODWG) 14.5%; no moderate or severe impairment.",
    notes = "Pooled phase I/II innovaTV 201 (NCT02001623, n = 195), innovaTV 202 (NCT02552121, n = 33), innovaTV 204 (NCT03438396, n = 101, cervical cancer) and innovaTV 207 (NCT03485209, n = 70). 4847 ADC and 5145 MMAE concentrations; BLQ records handled by the M3 method. Antidrug antibodies in 13 patients (3.3%) had no detectable effect. Categorical demographics from Gibiansky 2022 Table 1; continuous covariates from FDA BLA 761208 review Table 45."
  )

  ini({
    # ADC structural parameters (Gibiansky 2022 Table 2; reference patient
    # 75 kg, albumin 40 g/L, female).
    lcl <- log(1.42); label("ADC nonspecific linear clearance CL (L/day)") # Table 2 theta1 = 1.42
    lq <- log(4.01); label("ADC intercompartmental clearance Q (L/day)") # Table 2 theta2 = 4.01
    lvc <- log(3.10); label("ADC central volume Vc (L)") # Table 2 theta3 = 3.10
    lvp <- log(4.47); label("ADC peripheral volume Vp (L)") # Table 2 theta4 = 4.47
    lvmax <- log(3.35); label("ADC maximum Michaelis-Menten elimination rate Vmax (ug/mL/day)") # Table 2 theta5 = 3.35
    lkm <- log(3.44); label("ADC Michaelis-Menten constant KM (ug/mL)") # Table 2 theta6 = 3.44

    # MMAE structural parameters (Gibiansky 2022 Table 2). CL_MMAE and V_MMAE
    # are apparent values: the true ADC-to-MMAE conversion fraction is not
    # identifiable without MMAE dosing.
    lktr <- log(0.271); label("Rate constant of the MMAE delay compartment ktr (1/day)") # Table 2 theta9 = 0.271
    lcl_mmae <- log(42.8); label("Apparent MMAE clearance CL_MMAE (L/day)") # Table 2 theta10 = 42.8
    lvc_mmae <- log(2.09); label("Apparent MMAE central volume V_MMAE (L)") # Table 2 theta11 = 2.09
    lbeta <- log(0.0189); label("Rate constant of drug-to-antibody ratio decay with time after dose beta (1/day)") # Table 2 theta12 = 0.0189
    fr1 <- 0.0205; label("Fraction FR1 of DAR-scaled nonspecific ADC elimination delivered directly to central MMAE (unitless)") # Table 2 theta15 = 0.0205
    fr2 <- 0.0508; label("Fraction FR2 of DAR-scaled Michaelis-Menten ADC elimination delivered to the MMAE delay compartment (unitless)") # Table 2 theta16 = 0.0508
    dar0 <- fixed(3); label("Initial excess of the drug-to-antibody ratio over its floor of 1, DAR = 1 + dar0 * exp(-beta * tad) (unitless)") # Table S1 'DAR0 = 3' (DAR 4 at dosing, 1 at the floor)
    mwratio <- fixed(4.8); label("MMAE output scale factor: MMAE/ADC molecular-weight ratio times 1000 for mg/L -> ng/mL (unitless)") # Table S1 'MW_RATIO = 4.8'; Methods quotes 718/152 = 4.72

    # Covariate effects on ADC parameters (FDA BLA 761208 review Table 49,
    # theta17-theta23; Gibiansky 2022 Table S1 for the functional forms)
    e_wt_cl <- 0.495; label("Power exponent of WT/75 on ADC CL and Q (unitless)") # FDA review Table 49 theta17 = 0.495
    e_wt_vc <- 0.380; label("Power exponent of WT/75 on ADC Vc (unitless)") # FDA review Table 49 theta18 = 0.380
    e_wt_vp <- 0.622; label("Power exponent of WT/75 on ADC Vp (unitless)") # FDA review Table 49 theta19 = 0.622
    e_alb_cl <- -0.396; label("Power exponent of ALB/40 on ADC CL (unitless)") # FDA review Table 49 theta20 = -0.396
    e_alb_vc <- -0.197; label("Power exponent of ALB/40 on ADC Vc (unitless)") # FDA review Table 49 theta21 = -0.197
    e_sexm_cl <- 1.09; label("Multiplicative male effect on ADC CL, applied as e_sexm_cl^(1 - SEXF) (unitless)") # FDA review Table 49 theta22 = 1.09
    e_sexm_vc <- 1.13; label("Multiplicative male effect on ADC Vc, applied as e_sexm_vc^(1 - SEXF) (unitless)") # FDA review Table 49 theta23 = 1.13

    # Covariate effects on MMAE parameters (FDA BLA 761208 review Table 49,
    # theta24-theta35)
    e_wt_cl_mmae <- 0.457; label("Power exponent of WT/75 on CL_MMAE (unitless)") # FDA review Table 49 theta24 = 0.457
    e_wt_vc_mmae <- 0.895; label("Power exponent of WT/75 on V_MMAE (unitless)") # FDA review Table 49 theta25 = 0.895
    e_alb_cl_mmae <- 0.935; label("Power exponent of ALB/40 on CL_MMAE (unitless)") # FDA review Table 49 theta26 = 0.935
    e_tumtp_other_cl_mmae <- 1.22; label("Multiplicative non-cervical tumor effect on CL_MMAE (unitless)") # FDA review Table 49 theta27 = 1.22
    e_tumsz_cl_mmae <- -0.147; label("Exponential coefficient of (TUMSZ/60 - 1) on CL_MMAE (unitless)") # FDA review Table 49 theta28 = -0.147
    e_crcl_cl_mmae <- 0.271; label("Power exponent of CRCL/90 on CL_MMAE (unitless)") # FDA review Table 49 theta29 = 0.271
    e_ecog_cl_mmae <- 0.803; label("Multiplicative ECOG >= 1 effect on CL_MMAE (unitless)") # FDA review Table 49 theta30 = 0.803
    e_hepimp_cl_mmae <- 0.853; label("Multiplicative hepatic-impairment effect on CL_MMAE (unitless)") # FDA review Table 49 theta31 = 0.853
    e_ecog_vc_mmae <- 0.827; label("Multiplicative ECOG >= 1 effect on V_MMAE (unitless)") # FDA review Table 49 theta32 = 0.827
    e_alb_vc_mmae <- 0.575; label("Power exponent of ALB/40 on V_MMAE (unitless)") # FDA review Table 49 theta33 = 0.575
    e_age_ktr <- -0.252; label("Power exponent of AGE/60 on ktr (unitless)") # FDA review Table 49 theta34 = -0.252
    e_wt_ktr <- -0.175; label("Power exponent of WT/75 on ktr (unitless)") # FDA review Table 49 theta35 = -0.175
    e_study_innovatv202_ruv <- 1.42; label("Multiplicative innovaTV 202 study effect on the ADC residual SD (unitless)") # FDA review Table 49 theta36 = 1.42

    # Inter-individual variability (FDA BLA 761208 review Table 50, omega^2;
    # the CV column there is sqrt(omega^2), matching Gibiansky 2022 Results:
    # CL 23.2%, Vc 17.2%, Vp 14.4%, CL_MMAE 54.7%, V_MMAE 46.3%). The CL-Vp
    # covariance is zero (Table S1 omega block).
    etalcl + etalvc + etalvp ~ c(
      0.0538,
      0.0166, 0.0296,
      0, 0.0145, 0.0208
    ) # FDA review Table 50: Omega11 0.0538, Omega12 0.0166 (R 0.415), Omega22 0.0296, Omega23 0.0145 (R 0.586), Omega33 0.0208
    etalktr ~ 0.0212 # FDA review Table 50: Omega55 0.0212 (CV 14.5%)
    etalcl_mmae + etalvc_mmae ~ c(
      0.299,
      0.125, 0.215
    ) # FDA review Table 50: Omega66 0.299, Omega67 0.125 (R 0.495), Omega77 0.215
    etaRUV ~ 0.0561 # FDA review Table 50: Omega44 0.0561, IIV on the ADC residual SD (CV 23.7%)
    etaRUV_mmae ~ 0.0712 # FDA review Table 50: Omega88 0.0712, IIV on the MMAE residual SD (CV 26.7%)

    # Residual error: Y = IPRED + SD * eps, sigma^2 = 1 FIX, with
    # SD = sqrt(IPRED^2 * theta_prop^2 + theta_add^2) * study factor * exp(eta)
    # (Table S1).
    propSd <- 0.129; label("Proportional residual SD on ADC Cc (fraction)") # Table 2 theta7 = 0.129
    addSd <- 0.0173; label("Additive residual SD on ADC Cc (ug/mL)") # Table 2 theta8 = 0.0173
    propSd_mmae <- 0.282; label("Proportional residual SD on MMAE Cc_mmae (fraction)") # Table 2 theta13 = 0.282
    addSd_mmae <- 0.0113; label("Additive residual SD on MMAE Cc_mmae (ng/mL)") # Table 2 theta14 = 0.0113, printed with ug/mL units
  })

  model({
    # 1. Covariate multipliers (Gibiansky 2022 Table S1)
    cov_cl <- (WT / 75)^e_wt_cl * e_sexm_cl^(1 - SEXF) * (ALB / 40)^e_alb_cl
    cov_q <- (WT / 75)^e_wt_cl
    cov_vc <- (WT / 75)^e_wt_vc * e_sexm_vc^(1 - SEXF) * (ALB / 40)^e_alb_vc
    cov_vp <- (WT / 75)^e_wt_vp
    cov_cl_mmae <- (WT / 75)^e_wt_cl_mmae * (ALB / 40)^e_alb_cl_mmae *
      e_tumtp_other_cl_mmae^TUMTP_OTHER * (CRCL / 90)^e_crcl_cl_mmae *
      exp((TUMSZ / 60 - 1) * e_tumsz_cl_mmae) *
      e_ecog_cl_mmae^ECOG_GE1 * e_hepimp_cl_mmae^HEPIMP
    cov_vc_mmae <- (WT / 75)^e_wt_vc_mmae * e_ecog_vc_mmae^ECOG_GE1 *
      (ALB / 40)^e_alb_vc_mmae
    cov_ktr <- (AGE / 60)^e_age_ktr * (WT / 75)^e_wt_ktr

    # 2. Individual parameters
    cl <- exp(lcl + etalcl) * cov_cl
    q <- exp(lq) * cov_q
    vc <- exp(lvc + etalvc) * cov_vc
    vp <- exp(lvp + etalvp) * cov_vp
    vmax <- exp(lvmax)
    km <- exp(lkm)
    ktr <- exp(lktr + etalktr) * cov_ktr
    cl_mmae <- exp(lcl_mmae + etalcl_mmae) * cov_cl_mmae
    vc_mmae <- exp(lvc_mmae + etalvc_mmae) * cov_vc_mmae
    beta <- exp(lbeta)

    # 3. Micro-constants and the drug-to-antibody ratio, which decays from 4
    # at each dose towards 1 with time after the most recent dose
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_mmae <- cl_mmae / vc_mmae
    # tad() is NA before the first dose, when no ADC is present; 0 keeps the
    # ODE right-hand side finite there
    tad_dar <- tad()
    if (is.na(tad_dar)) {
      tad_dar <- 0
    }
    dar <- 1 + dar0 * exp(-beta * tad_dar)

    # Michaelis-Menten elimination flux (mg/day); vmax is a concentration rate,
    # so the amount rate is vmax * central / (km + Cc)
    mm_flux <- vmax * central / (km + central / vc)

    # 4. ODE system (Table S1; the peripheral equation is printed there with
    # -K12 * A1, a sign typo for the +K12 * A1 inflow of Figure 1)
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - mm_flux
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(precursor_mmae) <- dar * kel * central + fr2 * dar * mm_flux - ktr * precursor_mmae
    d/dt(central_mmae) <- fr1 * dar * kel * central + ktr * precursor_mmae - kel_mmae * central_mmae

    # 5. Observations: ADC in ug/mL, MMAE in ng/mL
    Cc <- central / vc
    Cc_mmae <- mwratio * central_mmae / vc_mmae

    ruv_fac <- e_study_innovatv202_ruv^STUDY_INNOVATV202 * exp(etaRUV)
    ruv_fac_mmae <- exp(etaRUV_mmae)
    propSd_i <- propSd * ruv_fac
    addSd_i <- addSd * ruv_fac
    propSd_mmae_i <- propSd_mmae * ruv_fac_mmae
    addSd_mmae_i <- addSd_mmae * ruv_fac_mmae

    Cc ~ add(addSd_i) + prop(propSd_i) + combined2()
    Cc_mmae ~ add(addSd_mmae_i) + prop(propSd_mmae_i) + combined2()
  })
}
