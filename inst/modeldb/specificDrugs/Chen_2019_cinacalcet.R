Chen_2019_cinacalcet <- function() {
  description <- "Two-compartment population PK model with lagged first-order absorption for oral cinacalcet, linked to a semi-mechanistic calcium-sensing-receptor (CaSR) allosteric PK/PD model of intact parathyroid hormone (pseudo-steady-state) and corrected serum calcium (turnover with PTH feedback), in adults and paediatrics with secondary hyperparathyroidism on dialysis and healthy adults"
  reference <- "Chen P, Narayanan A, Wu B, Sohn W, Melhem M. Bridging adults and paediatrics with secondary hyperparathyroidism receiving haemodialysis: a pharmacokinetic-pharmacodynamic analysis of cinacalcet. Br J Clin Pharmacol. 2019;85(6):1312-1325. doi:10.1111/bcp.13900"
  vignette <- "Chen_2019_cinacalcet"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL",
    pth = "pmol/L (intact PTH; 1 pmol/L = 9.43 pg/mL)",
    ca = "mmol/L (corrected serum calcium; 1 mmol/L = 4.0 mg/dL)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling with fixed exponents (0.75 on CL/F, 1 on central V/F). The reference weight is not printed in Chen 2019; 70 kg is assumed (see vignette Assumptions and deviations). Paediatric dosing in the source studies used post-dialysis (dry) weight.",
      source_name = "body weight"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-participant indicator (1 = healthy adult volunteer, study 970241; 0 = patient with CKD and secondary hyperparathyroidism on dialysis)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Selects the healthy-volunteer (HV) or patient set of PD parameters (PTH0, Ca0, Ratio, Kout,Ca, power on PTH, power on RO) and the matching iPTH / cCa residual errors from Chen 2019 Table 4. PK parameters are shared.",
      source_name = "HV / Patient (Table 4 subscripts)"
    ),
    STUDY_PHASE3 = list(
      description = "Phase 3 study indicator (1 = phase 3 study 20000172; 0 = the other seven studies)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Selects the PK residual-error magnitude only (Chen 2019 Table 3: 'sigma Ph3 studies prop' vs 'sigma other studies prop'). Set to 0 for simulation outside the phase 3 sampling design.",
      source_name = "Ph3 studies / other studies (Table 3)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "cinacalcet", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cinacalcet", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cinacalcet", units = "mg", specimen = "plasma", verified = TRUE),
    ca = list(analyte = "corrected calcium", units = "mmol/L", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 648L,
    n_studies = 8L,
    age_range = "28 days to 18 years (paediatric studies); adults (range not reported)",
    weight_range = "not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported",
    disease_state = "Secondary hyperparathyroidism (SHPT) with chronic kidney disease / end-stage renal disease on dialysis (adults and paediatrics); 79 healthy adult volunteers (study 970241)",
    dose_range = "Oral cinacalcet 1-300 mg: single doses, once-daily multiple doses and titrated once-daily regimens",
    regions = "Multinational (Amgen clinical studies)",
    notes = "Chen 2019 Table 1: paediatric studies 20070208 (n = 43, 6-18 y), 20110100 (n = 17, 28 days-6 y), 20030227 (n = 12, 6-17 y), 20090005 (n = 12, 28 days-6 y); adult SHPT studies 20000172 (phase 3, n = 403), 20000187 (n = 22), 980126 (n = 60); healthy-volunteer study 970241 (n = 79). Paediatric subjects contributed <15% of samples (Discussion). No baseline demographics table is printed. 34.7% of PK samples were below the 0.1 ng/mL LLOQ and excluded."
  )

  ini({
    # Pharmacokinetics -- Chen 2019 Table 3 (final population PK model), reference WT 70 kg (assumed)
    lka <- log(0.80); label("First-order absorption rate constant Ka (1/h)") # Table 3 'Ka' = 0.80 (printed unit 'L h-1' is a typo for h-1)
    lcl <- log(287.74); label("Apparent clearance CL/F at 70 kg (L/h)") # Table 3 'CL/F' = 287.74 L/h; Abstract
    lvc <- log(4104.15); label("Apparent central volume V2/F at 70 kg (L)") # Table 3 'V2/F' = 4104.15 L
    lvp <- log(19824.20); label("Apparent peripheral volume V3/F (L)") # Table 3 'V3/F' = 19824.20 L
    lq <- log(296.3); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 'Q/F' = 296.3 L/h
    ltlag <- log(0.36); label("Absorption lag time Tlag (h)") # Table 3 'Tlag' = 0.36 h
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT on CL/F (unitless)") # Results (after Table 2) 'exponent fixed to 0.75 for CL'
    e_wt_vc <- fixed(1); label("Allometric exponent of WT on V2/F (unitless)") # Results (after Table 2) 'and 1 for V1' (central volume per Discussion)

    # Semi-mechanistic CaSR model constants -- Chen 2019 Table 2, Equation 1
    kd_ca <- fixed(1.2); label("Equilibrium dissociation constant of ionized Ca at the CaS receptor KD (mmol/L)") # Table 2 Eq. 1 'KD ... was assumed to be 1.2 mM'
    f_ionized_ca <- fixed(0.45); label("Fraction of corrected Ca that is ionized (unitless)") # Table 2 Eq. 1 'Ionized Ca was calculated empirically as 0.45 of cCa'

    # PD parameters, patient stratum -- Chen 2019 Table 4
    lrbase_pth_patient <- log(54.29); label("Baseline iPTH, patients PTH0 (pmol/L)") # Table 4 'PTH0 Patient' = 54.29 pM
    lrbase_ca_patient <- log(2.45); label("Baseline corrected Ca, patients Ca0 (mmol/L)") # Table 4 'Ca0 Patient' = 2.45 mM
    lratio_patient <- log(0.072); label("Allosteric ratio alpha/Ki, patients (mL/ng; multiplies Cc in ng/mL)") # Table 4 'RatioPatient' = 0.072 (printed 'L mM-1'; see vignette for the concentration unit)
    lkout_ca_patient <- log(0.03); label("First-order Ca elimination rate Kout,Ca, patients (1/h)") # Table 4 'Kout,Ca Patient' = 0.03 (printed 'L h-1'); t1/2 23.1 h (Results)
    lgamma_pth_patient <- log(0.07); label("Power on PTH change for Ca production (gamma), patients (unitless)") # Table 4 'Power on PTHPatient' = 0.07; Table 2 Eq. 3
    lambda_ro_patient <- -3.30; label("Power on receptor-occupancy change for PTH (lambda), patients (unitless)") # Table 4 'Power on ROPatient' = -3.30; Table 2 Eq. 2

    # PD parameters, healthy-volunteer stratum -- Chen 2019 Table 4
    lrbase_pth_hv <- log(2.38); label("Baseline iPTH, healthy volunteers PTH0 (pmol/L)") # Table 4 'PTH0 HV' = 2.38 pM
    lrbase_ca_hv <- log(2.30); label("Baseline corrected Ca, healthy volunteers Ca0 (mmol/L)") # Table 4 'Ca0 HV' = 2.30 mM
    lratio_hv <- log(0.02); label("Allosteric ratio alpha/Ki, healthy volunteers (mL/ng)") # Table 4 'RatioHV' = 0.02
    lkout_ca_hv <- log(0.01); label("First-order Ca elimination rate Kout,Ca, healthy volunteers (1/h)") # Table 4 'Kout,Ca HV' = 0.01
    lgamma_pth_hv <- log(1.80); label("Power on PTH change for Ca production (gamma), healthy volunteers (unitless)") # Table 4 'Power on PTHHV' = 1.80
    lambda_ro_hv <- -10.03; label("Power on receptor-occupancy change for PTH (lambda), healthy volunteers (unitless)") # Table 4 'Power on ROHV' = -10.03

    # IIV -- variances omega^2 = log(1 + CV^2) from the %CV in Tables 3 and 4
    etalka ~ 0.71315 # Table 3 'omega Ka' = 102% CV
    etalcl + etalvc ~ c(0.45590, 0.35814, 0.38940) # Table 3 'omega CL/F' = 76% CV, 'omega V2/F' = 69% CV, correlation 0.85
    etalvp ~ 0.55366 # Table 3 'omega V3/F' = 86% CV
    etalq ~ 0.38008 # Table 3 'omega Q/F' = 68% CV
    etaltlag ~ 0.10939 # Table 3 'omega Tlag' = 34% CV
    etalrbase_pth ~ 0.34331 # Table 4 'omega PTH0' = 64.00% CV
    etalrbase_ca ~ 0.0043465 # Table 4 'omega Ca0' = 6.60% CV
    etalratio ~ 0.49968 # Table 4 'omega Ratio' = 80.51% CV
    etalgamma_pth ~ 0.81795 # Table 4 'omega Power on PTH' = 112.51% CV
    etalambda_ro ~ 0.44696 # Table 4 'omega Power on RO' = 75.07% CV (exponential IIV on a negative typical value)

    # Residual error -- additive on the log scale (lnorm); SD = sqrt(log(1 + CV^2)) from the printed %CV
    expSd_nonph3 <- 0.47238; label("Log-scale residual SD for cinacalcet, non-phase-3 studies") # Table 3 'sigma other studies prop' = 50% CV
    expSd_ph3 <- 0.64626; label("Log-scale residual SD for cinacalcet, phase 3 study 20000172") # Table 3 'sigma Ph3 studies prop' = 72% CV
    expSd_pth_patient <- 0.41190; label("Log-scale residual SD for iPTH, patients") # Table 4 'sigma iPTH Patient' = 43% CV
    expSd_ca_patient <- 0.069914; label("Log-scale residual SD for corrected Ca, patients") # Table 4 'sigma cCa Patient' = 7% CV
    expSd_pth_hv <- 0.22704; label("Log-scale residual SD for iPTH, healthy volunteers") # Table 4 'sigma iPTH HV' = 23% CV
    expSd_ca_hv <- 0.037287; label("Log-scale residual SD for corrected Ca, healthy volunteers") # Table 4 'sigma cCa_HV' = 3.73% CV
  })

  model({
    # Pharmacokinetics (Chen 2019 Results / Figure 1): two compartments,
    # lagged first-order absorption, allometry on CL/F and central V/F only.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # PD parameters: patient or healthy-volunteer stratum (Chen 2019 Table 4)
    hv <- DIS_HEALTHY
    rbase_pth <- exp(
      lrbase_pth_patient * (1 - hv) + lrbase_pth_hv * hv + etalrbase_pth
    )
    rbase_ca <- exp(
      lrbase_ca_patient * (1 - hv) + lrbase_ca_hv * hv + etalrbase_ca
    )
    ratio <- exp(lratio_patient * (1 - hv) + lratio_hv * hv + etalratio)
    kout_ca <- exp(lkout_ca_patient * (1 - hv) + lkout_ca_hv * hv)
    gamma_pth <- exp(
      lgamma_pth_patient * (1 - hv) + lgamma_pth_hv * hv + etalgamma_pth
    )
    lambda_ro <- (lambda_ro_patient * (1 - hv) + lambda_ro_hv * hv) *
      exp(etalambda_ro)

    # Table 2 Eq. 6: Kin,Ca = Ca0 * Kout,Ca
    kin_ca <- rbase_ca * kout_ca

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    # Plasma cinacalcet: mg / L * 1000 = ng/mL
    Cc <- 1000 * central / vc

    # Table 2 Eq. 1: baseline CaSR occupancy by ionized Ca without drug
    ca_rel0 <- f_ionized_ca * rbase_ca / kd_ca
    rho0 <- ca_rel0 / (ca_rel0 + 1)
    # Table 2 Eq. 5: occupancy with cinacalcet, Ratio = alpha / Ki
    ca_rel <- f_ionized_ca * ca / kd_ca * (1 + ratio * Cc)
    rho <- ca_rel / (ca_rel + 1)
    # Table 2 Eq. 2: pseudo-steady-state iPTH
    pth <- rbase_pth * (rho / rho0)^lambda_ro

    # Table 2 Eq. 3: corrected-Ca turnover, production driven by iPTH change
    ca(0) <- rbase_ca
    d/dt(ca) <- kin_ca * (pth / rbase_pth)^gamma_pth - kout_ca * ca

    expSd <- expSd_nonph3 * (1 - STUDY_PHASE3) + expSd_ph3 * STUDY_PHASE3
    expSd_pth <- expSd_pth_patient * (1 - hv) + expSd_pth_hv * hv
    expSd_ca <- expSd_ca_patient * (1 - hv) + expSd_ca_hv * hv

    Cc ~ lnorm(expSd)
    pth ~ lnorm(expSd_pth)
    ca ~ lnorm(expSd_ca)
  })
}
