Cosson_2025_crovalimab <- function() {
  description <- "Two-compartment population PK model of crovalimab (anti-C5 recycling monoclonal antibody) with first-order subcutaneous absorption and logit-normal bioavailability, plus an exponentially decaying additional clearance in patients switching from another C5 inhibitor (transient drug-target-drug immune complexes), in healthy volunteers and patients with paroxysmal nocturnal haemoglobinuria from the COMPOSER and COMMODORE 1, 2 and 3 trials (Cosson 2025)"
  reference <- "Cosson V, Fu R, Kulasekararaj A, Nishimura JI, Panse J, Roth A, Scheinberg P, Tong H, Yoon SS, Beveridge L, Gotanda K, Jaminion F, Henrich A, Lundberg P, Shi D, Sreckovic S, Zhang Y, Zhang Z, Benkali K, Buatois S. Pharmacokinetic characterization and exposure-response relationship of crovalimab in the COMMODORE 1, 2 and 3 and COMPOSER trials of patients with paroxysmal nocturnal haemoglobinuria. Br J Clin Pharmacol. 2025;91(5):1479. doi:10.1111/bcp.16394"
  vignette <- "Cosson_2025_crovalimab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
  dosing <- c("depot", "central")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling centred on 75 kg: (WT/75)^0.75 (fixed) on CL, on the switch-patient clearance CLs and on Q; (WT/75)^0.683824 (estimated) on V2 and V3 (supplement control stream $PK CLWT / VWT and Supplemental Results covariate equations).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the absorption rate constant centred on 36 years: ka = 0.126 * (AGE/36)^-0.459 (supplement control stream KAAGE; Supplemental Results covariate equations; Table S3 'Age on Ka').",
      source_name = "AGE"
    ),
    ADA_TITER = list(
      description = "Anti-crovalimab antibody titre",
      units = "reciprocal dilution titre (ADA-negative = 0)",
      type = "continuous",
      reference_category = NULL,
      notes = "The source model uses a three-group titre categorisation (control stream column ADA3T3G); groups 2 and 3 ('medium or high titer ADA', Supplemental Results) multiply CL by (1 + 0.247352). The Results text identifies these as titres >= 10^4 (8.1% of the analysis population), so the model applies the effect when ADA_TITER >= 10000. ADA-negative subjects carry ADA_TITER = 0. Whether the source dataset carried the titre group as time-varying or per-subject is not stated.",
      source_name = "ADA3T3G"
    ),
    CONMED_ECULIZUMAB_BL = list(
      description = "Receiving eculizumab (standard or higher-than-approved dose) at the start of crovalimab, i.e. switching from eculizumab",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (C5-inhibitor-naive at crovalimab start)",
      notes = "The source control stream switches on the additional clearance CLs when SWITCH != 0, where SWITCH codes switching from standard-dose eculizumab, high-dose eculizumab or ravulizumab. The paper found no difference between the three on CLs (Discussion), so the model uses one switch indicator: switch_c5i = max(CONMED_ECULIZUMAB_BL, CONMED_RAVULIZUMAB_BL).",
      source_name = "SWITCH"
    ),
    CONMED_RAVULIZUMAB_BL = list(
      description = "Receiving ravulizumab at the start of crovalimab, i.e. switching from ravulizumab",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (C5-inhibitor-naive at crovalimab start)",
      notes = "Second level of the source SWITCH column (21 ravulizumab switch patients in COMMODORE 1 arm C, Table S1). Acts identically to CONMED_ECULIZUMAB_BL; at most one of the two should be 1.",
      source_name = "SWITCH"
    ),
    SAMPLE_INTENSIVE = list(
      description = "Rich PK sampling indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (sparse sampling)",
      notes = "Selects the additive residual-error magnitude. RICH = 1 for COMPOSER Part 1 (healthy volunteers, study 39144 PART 1) and study YO42311; all other observations are sparse (supplement control stream $PK RICH). Rich-sampling additive SD 1.01 ug/mL; sparse-sampling additive SD 14.4 ug/mL with between-subject variability on its magnitude.",
      source_name = "RICH"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "crovalimab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "crovalimab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "crovalimab", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 430L,
    n_studies = 4L,
    age_range = "13-85 years",
    age_median = "35.5 years (naive patients), 44.0 years (switch patients), 36.0 years (healthy volunteers)",
    weight_range = "40.6-140 kg",
    weight_median = "65.2 kg (naive patients), 71.0 kg (switch patients), 72.8 kg (healthy volunteers); covariate reference 75 kg",
    disease_state = "Paroxysmal nocturnal haemoglobinuria: 210 patients naive to complement inhibition and 211 patients switching from eculizumab (standard or higher-than-approved dose) or ravulizumab; plus 9 healthy volunteers (COMPOSER Part 1).",
    dose_range = "Phase III regimen (Table 1): 1000 mg IV on day 1 (1500 mg if >= 100 kg), 340 mg SC on days 2, 8, 15 and 22, then 680 mg SC (1020 mg if >= 100 kg) on day 29 and every 4 weeks; COMPOSER used dose-escalation and earlier regimens.",
    regions = "Global (COMMODORE 1 and 2, COMPOSER) and China (COMMODORE 3).",
    notes = "Pooled data from COMPOSER (NCT03157635, phase I/II), COMMODORE 1 (NCT04432584), COMMODORE 2 (NCT04434092) and COMMODORE 3 (NCT04654468) (Table S1). Baseline demographics in Table 2; sex and race distributions are not reported. Total serum crovalimab by ELISA, LLOQ 50.0 ng/mL. ADA titres >= 10^4 in 8.1% of the analysis population."
  )

  ini({
    # Structural parameters (Table S3 final model; reference 75 kg, 36 years)
    lcl_exp_inf <- log(0.0791); label("Clearance CL once the switch-related clearance has decayed (L/day)") # Table S3 'CL' = 0.0791 L/day (RSE 6.70%)
    lcl_exp_component <- log(0.0341); label("Initial additional clearance CLs in switch patients (L/day)") # Table S3 'CLs' = 0.0341 L/day (RSE 17.6%)
    lcl_exp_kdes <- log(0.0144); label("Exponential decay rate constant Kd of CLs (1/day)") # Table S3 'Kd' = 0.0144 1/day (RSE 20.6%)
    lvc <- log(3.23); label("Central volume of distribution V2 (L)") # Table S3 'V2' = 3.23 L (RSE 1.30%)
    lq <- log(0.168); label("Inter-compartmental clearance Q (L/day)") # Table S3 'Q' = 0.168 L/day (RSE 12.4%)
    lvp <- log(2.32); label("Peripheral volume of distribution V3 (L)") # Table S3 'V3' = 2.32 L (RSE 7.17%)
    lka <- log(0.126); label("First-order SC absorption rate constant Ka (1/day)") # Table S3 'Ka' = 0.126 1/day (RSE 13.0%)
    logitfdepot <- log(0.830 / (1 - 0.830)); label("Logit of SC bioavailability F (F = 0.830)") # Table S3 'F' = 0.830 (RSE 7.72%); control stream F1 = expit(logit(TVF1) + ETA(4))

    # Covariate effects (Supplemental Results covariate equations; Table S3)
    e_wt_cl <- fixed(0.75); label("Power exponent of body weight on CL, CLs and Q (unitless; reference 75 kg)") # Table S3 'WT on CL and Q' = 0.750 (FIX); control stream THETA(10) 0.75 FIX, CLS = CLSP*CLWT
    e_wt_vc <- 0.683824; label("Power exponent of body weight on V2 and V3 (unitless; reference 75 kg)") # Supplemental Results equations exponent 0.683824; Table S3 'WT on V2 and V3' = 0.684 (RSE 8.07%)
    e_ada_cl_exp_inf <- 0.247352; label("Fractional increase in CL with ADA titre >= 10^4 (unitless)") # Supplemental Results equation (1 + 0.247352); Table S3 'ADA on CL' = 0.247 (RSE 33.5%)
    e_age_ka <- -0.459; label("Power exponent of age on Ka (unitless; reference 36 years)") # Table S3 'Age on Ka' = -0.459 (RSE 26.5%); Supplemental Results equation (age/36)^-0.459

    # IIV. Table S3 reports IIV on the SD scale (labelled CV; RSE 'on the
    # approximate SD scale'); variance = value^2. The supplement $OMEGA block
    # (0.0418, 0.0496, 0.140, 1.39, 0.494, 0.336, 0.617) confirms the SD reading.
    # IIV on CLs, Kd and the proportional RUV are 0 FIX in the source and omitted.
    etalcl_exp_inf ~ 0.042436 # Table S3 'IIV CL' 0.206; 0.206^2
    etalvc ~ 0.050176 # Table S3 'IIV V2' 0.224; 0.224^2
    etalka ~ 0.146689 # Table S3 'IIV Ka' 0.383; 0.383^2
    etalogitfdepot ~ 1.3456 # Table S3 'IIV F' 1.16 (logit scale); 1.16^2
    etalvp ~ 0.498436 # Table S3 'IIV V3' 0.706; 0.706^2
    etalq ~ 0.341056 # Table S3 'IIV Q' 0.584; 0.584^2
    etaaddSd ~ 0.619369 # Table S3 'IIV additive RUV' 0.787; 0.787^2 (sparse-sampling additive SD only)

    # Residual error: W = sqrt(PROP^2 * IPRED^2 + ADD^2), SIGMA 1 FIX
    propSd <- 0.113; label("Proportional residual error (fraction)") # Table S3 'Proportional RUV' = 0.113 (RSE 3.92%)
    addSdIntensive <- 1.01; label("Additive residual error, rich sampling (ug/mL)") # Table S3 'Additive RUV rich sampling' SD = 1.01 (RSE 23.4%)
    addSdSparse <- 14.4; label("Typical additive residual error, sparse sampling (ug/mL)") # Table S3 'Additive RUV sparse sampling' SD = 14.4 (RSE 10.9%)
  })

  model({
    # 1. Derived covariate terms (supplement control stream $PK)
    switch_c5i <- max(CONMED_ECULIZUMAB_BL, CONMED_RAVULIZUMAB_BL)
    ada_high <- 0
    if (ADA_TITER >= 10000) ada_high <- 1
    wt_cl <- (WT / 75)^e_wt_cl
    wt_v <- (WT / 75)^e_wt_vc

    # 2. Individual parameters
    cl_exp_inf <- exp(lcl_exp_inf + etalcl_exp_inf) * wt_cl * (1 + e_ada_cl_exp_inf * ada_high)
    cl_exp_component <- exp(lcl_exp_component) * wt_cl * switch_c5i
    cl_exp_kdes <- exp(lcl_exp_kdes)
    vc <- exp(lvc + etalvc) * wt_v
    vp <- exp(lvp + etalvp) * wt_v
    q <- exp(lq + etalq) * wt_cl
    ka <- exp(lka + etalka) * (AGE / 36)^e_age_ka
    fdepot <- expit(logitfdepot + etalogitfdepot)

    # Total clearance: switch patients carry an extra clearance from transient
    # drug-target-drug immune complexes that decays exponentially with time
    # since the first crovalimab dose (control stream $DES CLN + CLS*EXP(-KD*T)).
    cl <- cl_exp_inf + cl_exp_component * exp(-cl_exp_kdes * t)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system (IV doses into central, SC doses into depot)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Bioavailability (SC only)
    f(depot) <- fdepot

    # 6. Observation and error (mg/L = ug/mL)
    addSd <- addSdIntensive * SAMPLE_INTENSIVE + addSdSparse * exp(etaaddSd) * (1 - SAMPLE_INTENSIVE)
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
