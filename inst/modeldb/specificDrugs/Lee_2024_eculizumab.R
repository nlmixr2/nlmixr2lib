Lee_2024_eculizumab <- function() {
  description <- "Two-compartment population PK model of eculizumab (SB12 proposed biosimilar and reference Soliris) with a direct inhibitory sigmoid Imax model for terminal complement activity and a direct sigmoid Emax model linking terminal complement activity to serum LDH, in healthy subjects and patients with paroxysmal nocturnal haemoglobinuria (Lee 2024)"
  reference <- "Lee H, Park J, Jang H, Lee SJ, Kim J. Population pharmacokinetic, pharmacodynamic and efficacy modeling of SB12 (proposed eculizumab biosimilar) and reference eculizumab. Eur J Clin Pharmacol. 2024;80(9):1325-1338. doi:10.1007/s00228-024-03703-8"
  vignette <- "Lee_2024_eculizumab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power (allometric) scaling on CL and Vc, centred on the pooled-population median weight of 80.9 kg: CL = theta1 * (WT/80.9)^theta7 and Vc = theta2 * (WT/80.9)^theta8 (Table 3 covariate equations). Exponents are estimated, not fixed: 1.1400 on CL and 0.8630 on Vc. Baseline weight; the paper's covariate screen used a power model centred on the median.",
      source_name        = "WT"
    ),
    DIS_PNH = list(
      description        = "Paroxysmal nocturnal haemoglobinuria patient status (subject group: PNH patient versus healthy subject)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy subject enrolled in the phase I study)",
      notes              = "Subject group was retained as a categorical covariate with separately estimated typical values for each category, on three parameters: Vc (Table 3: 3.47 L healthy versus 5.68 L PNH), baseline terminal complement activity E0 (Table 4: 85.90% healthy versus 101.00% PNH) and maximum inhibition Imax (Table 4: 0.93 healthy versus 0.88 PNH). Encoded here as log-scale shifts relative to the healthy typical value, so e_pnh_vc = log(5.68/3.47), e_pnh_rbase_tca = log(101.00/85.90) and e_pnh_imax_tca = log(0.88/0.93) reproduce the published patient typical values exactly. The efficacy (LDH) sub-model was fitted in PNH patients only.",
      source_name        = "subject group (healthy subjects versus PNH patients)"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "eculizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "eculizumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 289L,
    n_studies      = 2L,
    age_range      = "18-79 years (phase I median 40, range 19-55; phase III median 36, range 18-79)",
    age_median     = "40 years (phase I healthy subjects); 36 years (phase III PNH patients)",
    weight_range   = "43.0-111.0 kg (phase I 70.0-94.3; phase III 43.0-111.0)",
    weight_median  = "82.40 kg (phase I healthy subjects); 63.00 kg (phase III PNH patients); 80.9 kg pooled reference weight used in the covariate model",
    sex_female_pct = 12.8,
    race_ethnicity = "Phase I: 95.8% White, 0.4% Black or African American, 0.4% Asian, 0.4% American Indian or Alaska Native, 2.9% Other. Phase III: 36.7% White, 53.1% Asian, 10.2% Native Hawaiian or Other Pacific Islander, 36.7% Other (Table 2; categories are not mutually exclusive as reported).",
    disease_state  = "240 healthy subjects (phase I) and 49 patients with paroxysmal nocturnal haemoglobinuria (phase III).",
    dose_range     = "Phase I: 300 mg single intravenous infusion over 35 min. Phase III: 600 mg IV every week for 4 weeks (induction), then 900 mg at week 5, then 900 mg every 2 weeks through week 50 (maintenance).",
    regions        = "Multi-national phase III (Germany, India, Republic of Korea, Malaysia, Mexico, Romania, Taiwan, Thailand, Ukraine); phase I region not stated.",
    notes          = "Pooled randomised double-blind phase I three-arm single-dose study in healthy subjects (SB12, EU-sourced Soliris, US-sourced Soliris; 80 per arm) and a randomised double-blind multicentre cross-over phase III study in PNH patients (SB12-to-ECU and ECU-to-SB12 sequences). 4136 quantifiable serum concentrations, 2900 terminal complement activity levels from 289 subjects and 1350 LDH levels from 49 PNH patients. Treatment group (SB12 versus reference eculizumab) was NOT a significant covariate on any PK, PD or efficacy parameter, so one parameter set describes both products. PK LLOQ 0.8 ug/mL (14.1% of serum concentrations were below LLOQ and treated as missing); PD assay range 10-125%."
  )

  ini({
    # -------------------------------------------------------------------------
    # Structural PK parameters, pooled healthy + PNH final model
    # (Lee 2024 Table 3). Reference weight 80.9 kg.
    # -------------------------------------------------------------------------
    lcl <- log(0.0174); label("Clearance CL (L/h)")                                 # Lee 2024 Table 3 (theta1)
    lvc <- log(3.47);   label("Central volume of distribution Vc in healthy subjects (L)")  # Lee 2024 Table 3 (theta2)
    lvp <- log(0.79);   label("Peripheral volume of distribution Vp (L)")           # Lee 2024 Table 3 (Vp)
    lq  <- log(0.0134); label("Inter-compartmental clearance Q (L/h)")              # Lee 2024 Table 3 (Q)

    # -------------------------------------------------------------------------
    # PK covariate effects (Lee 2024 Table 3)
    #   CL = theta1 * (WT/80.9)^theta7
    #   Vc = theta2 * (WT/80.9)^theta8, with theta2 estimated separately for
    #        healthy subjects (3.47 L) and PNH patients (5.68 L)
    # -------------------------------------------------------------------------
    e_wt_cl  <- 1.1400; label("Power exponent of body weight on CL (unitless; reference 80.9 kg)")  # Lee 2024 Table 3 (theta7)
    e_wt_vc  <- 0.8630; label("Power exponent of body weight on Vc (unitless; reference 80.9 kg)")  # Lee 2024 Table 3 (theta8)
    e_pnh_vc <- log(5.68 / 3.47); label("Log-scale shift of Vc in PNH patients relative to healthy subjects (unitless)")  # Lee 2024 Table 3 (theta2_pat 5.68 L vs theta2 3.47 L)

    # -------------------------------------------------------------------------
    # Inter-individual variability on PK parameters (Lee 2024 Table 3, %CV;
    # ETAs applied exponentially per Methods "Base model").
    # omega^2 = log(CV^2 + 1):
    #   15.62 % -> log(0.1562^2 + 1) = 0.0241056  (CL)
    #   12.74 % -> log(0.1274^2 + 1) = 0.0161004  (Vc)
    #   36.80 % -> log(0.3680^2 + 1) = 0.1270061  (Vp)
    # Covariance from the reported correlation rho(CL,Vc) = 0.54:
    #   0.54 * sqrt(0.0241056 * 0.0161004) = 0.0106383
    # -------------------------------------------------------------------------
    etalcl + etalvc ~ c(0.0241056,
                        0.0106383, 0.0161004)  # Lee 2024 Table 3: IIV CL 15.62 %CV, Vc 12.74 %CV, rho 0.54
    etalvp ~ 0.1270061                          # Lee 2024 Table 3: IIV Vp 36.80 %CV

    # -------------------------------------------------------------------------
    # Pharmacodynamics: terminal complement activity (%). Direct response
    # model with an inhibitory sigmoid Imax relationship to serum eculizumab
    # concentration (Lee 2024 Results / Discussion; parameters Table 4).
    #   tca = E0 * (1 - Imax * Cc^H / (IC50^H + Cc^H))
    # E0 and Imax were estimated separately by subject group.
    # -------------------------------------------------------------------------
    lrbase_tca <- log(85.90); label("Baseline terminal complement activity E0 in healthy subjects (%)")  # Lee 2024 Table 4 (E0)
    limax_tca  <- log(0.93);  label("Maximum fractional inhibition of terminal complement activity Imax in healthy subjects (unitless)")  # Lee 2024 Table 4 (Imax)
    lec50_tca  <- log(36.60); label("Serum eculizumab concentration achieving 50% of Imax, IC50 (ug/mL)")  # Lee 2024 Table 4 (IC50)
    lhill_tca  <- log(4.56);  label("Hill coefficient of the concentration-terminal complement activity relationship (unitless)")  # Lee 2024 Table 4 (H)

    e_pnh_rbase_tca <- log(101.00 / 85.90); label("Log-scale shift of baseline terminal complement activity in PNH patients relative to healthy subjects (unitless)")  # Lee 2024 Table 4 (E0_pat 101.00% vs E0 85.90%)
    e_pnh_imax_tca  <- log(0.88 / 0.93);    label("Log-scale shift of Imax in PNH patients relative to healthy subjects (unitless)")  # Lee 2024 Table 4 (Imax_pat 0.88 vs Imax 0.93)

    # omega^2 = log(CV^2 + 1):
    #   15.29 % -> 0.0231093 (E0)
    #    2.65 % -> 0.0007020 (Imax)
    #   22.74 % -> 0.0504181 (IC50)
    etalrbase_tca ~ 0.0231093   # Lee 2024 Table 4: IIV E0   = 15.29 %CV
    etalimax_tca  ~ 0.0007020   # Lee 2024 Table 4: IIV Imax =  2.65 %CV
    etalec50_tca  ~ 0.0504181   # Lee 2024 Table 4: IIV IC50 = 22.74 %CV

    # -------------------------------------------------------------------------
    # Efficacy: serum lactate dehydrogenase (U/L). Direct sigmoid Emax model
    # driven by terminal complement activity (Lee 2024 Results / Discussion;
    # parameters Table 5, fitted in PNH patients only).
    #   ldh = LL0 + LMAX * tca^LGAM / (LC50^LGAM + tca^LGAM)
    # -------------------------------------------------------------------------
    lrbase_ldh <- log(206);  label("Baseline LDH in the terminal complement activity-LDH relationship LL0 (U/L)")  # Lee 2024 Table 5 (LL0)
    lemax_ldh  <- log(1680); label("Maximum LDH increment in the terminal complement activity-LDH relationship LMAX (U/L)")  # Lee 2024 Table 5 (LMAX)
    lec50_ldh  <- log(39);   label("Terminal complement activity achieving 50% of LMAX, LC50 (%)")  # Lee 2024 Table 5 (LC50)
    lhill_ldh  <- log(4.30); label("Hill coefficient of the terminal complement activity-LDH relationship LGAM (unitless)")  # Lee 2024 Table 5 (LGAM)

    # omega^2 = log(CV^2 + 1):
    #   28.97 % -> 0.0805897 (LL0)
    #   81.57 % -> 0.5100452 (LMAX)
    #   55.91 % -> 0.2720044 (LGAM)
    etalrbase_ldh ~ 0.0805897   # Lee 2024 Table 5: IIV LL0  = 28.97 %CV
    etalemax_ldh  ~ 0.5100452   # Lee 2024 Table 5: IIV LMAX = 81.57 %CV
    etalhill_ldh  ~ 0.2720044   # Lee 2024 Table 5: IIV LGAM = 55.91 %CV

    # -------------------------------------------------------------------------
    # Residual error. All three endpoints used a proportional error model
    # (Lee 2024 Tables 3, 4 and 5).
    # -------------------------------------------------------------------------
    propSd     <- 0.1170; label("Proportional residual error for serum eculizumab concentration (fraction)")  # Lee 2024 Table 3 (sigma_prop = 11.70%)
    propSd_tca <- 0.1860; label("Proportional residual error for terminal complement activity (fraction)")    # Lee 2024 Table 4 (sigma_prop = 18.60%)
    propSd_ldh <- 0.3200; label("Proportional residual error for LDH (fraction)")                             # Lee 2024 Table 5 (sigma_prop = 32.00%)
  })

  model({
    # 1. Reference weight for the power covariate model (Lee 2024 Table 3).
    WT_ref <- 80.9

    # 2. Individual PK parameters. Subject group enters Vc as a log-scale
    #    shift so that exp(lvc + e_pnh_vc) reproduces the published PNH
    #    typical value of 5.68 L exactly.
    cl <- exp(lcl + etalcl) * (WT / WT_ref)^e_wt_cl
    vc <- exp(lvc + e_pnh_vc * DIS_PNH + etalvc) * (WT / WT_ref)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system (two-compartment, intravenous, first-order elimination)
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Individual PD (terminal complement activity) parameters
    rbase_tca <- exp(lrbase_tca + e_pnh_rbase_tca * DIS_PNH + etalrbase_tca)
    imax_tca  <- exp(limax_tca + e_pnh_imax_tca * DIS_PNH + etalimax_tca)
    ec50_tca  <- exp(lec50_tca + etalec50_tca)
    hill_tca  <- exp(lhill_tca)

    # 6. Individual efficacy (LDH) parameters
    rbase_ldh <- exp(lrbase_ldh + etalrbase_ldh)
    emax_ldh  <- exp(lemax_ldh + etalemax_ldh)
    ec50_ldh  <- exp(lec50_ldh)
    hill_ldh  <- exp(lhill_ldh + etalhill_ldh)

    # 7. Observations.
    #    Cc  : serum eculizumab concentration (ug/mL = mg/L)
    #    tca : terminal complement activity (%), direct inhibitory sigmoid Imax
    #    ldh : serum lactate dehydrogenase (U/L), direct sigmoid Emax driven by
    #          terminal complement activity
    #    tca_pos floors the complement-activity driver of the LDH relationship
    #    at zero. Imax carries an exponential ETA (Table 4) and is therefore
    #    not bounded above by 1, so a small fraction of simulated subjects can
    #    take imax_tca > 1 and drive tca slightly negative at full complement
    #    blockade; tca^hill_ldh would then be NaN. The floor changes nothing
    #    for any subject with a physically meaningful (non-negative) tca.
    Cc      <- central / vc
    tca     <- rbase_tca * (1 - imax_tca * Cc^hill_tca / (ec50_tca^hill_tca + Cc^hill_tca))
    tca_pos <- max(tca, 0)
    ldh     <- rbase_ldh + emax_ldh * tca_pos^hill_ldh / (ec50_ldh^hill_ldh + tca_pos^hill_ldh)

    Cc  ~ prop(propSd)
    tca ~ prop(propSd_tca)
    ldh ~ prop(propSd_ldh)
  })
}
