Han_2015_sibutramine <- function() {
  description <- "Two-compartment population PK for the active mono-desmethyl metabolite M1 plus a one-compartment PK for the downstream di-desmethyl metabolite M2 of the appetite-suppressant prodrug sibutramine, combined with an asymptotic exposure-response weight-loss PD model in Korean obese adults with metabolic syndrome. Sibutramine is dosed orally and assumed to convert entirely to M1 during absorption; M1 is then metabolised entirely to M2 and M2 is the only elimination pathway. Drug effect inhibits the rate of weight gain via a sigmoid Emax function of the steady-state sum AUC of M1 and M2 (AUC_ss,sum, computed from the current daily dose and the individual M1 and M2 clearances). A constant placebo effect is acknowledged only in female subjects and scales with mean-normalised baseline BMI."
  reference <- paste(
    "Han S., Jeon S., Hong T., Lee J., Bae S. H., Park W.-S., Park G.-J.,",
    "Youn S., Jang D. Y., Kim K.-S., Yim D.-S. (2015).",
    "Exposure-response model for sibutramine and placebo: suggestion for",
    "application to long-term weight-control drug development.",
    "Drug Des Devel Ther 9:5185-5194.",
    "doi:10.2147/DDDT.S85435.",
    sep = " "
  )
  vignette <- "Han_2015_sibutramine"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear deviation around a reference of 35 years applied to M1 clearance: `cl = cl_typ * (1 + e_age_cl * (AGE - 35))` with estimated coefficient 0.0120 1/year (Han 2015 Methods, CL_M1 equation; Table 2 C_AGE = 0.0120 yr^-1). Older patients have higher CL_M1 (i.e., faster M1-to-M2 metabolism), consistent with Hind et al. 2007 finding that the M1/M2 AUC ratio rises with age (Han 2015 Discussion).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male).",
      notes              = "Han 2015 codes SEX as 0 for males and 1 for females, which matches the nlmixr2lib canonical SEXF. SEXF enters the PD model in two places: (1) baseline body weight `base = mw - SEXF * fwc` (Table 2 BASE = MW - SEX * FWC, with MW = 89.1 kg the typical male baseline and FWC = 11.4 kg the female correction), and (2) the placebo effect gate `p_max = (p_fem + etap_fem) * SEXF * (BMI/30.1)^bex`, which makes the placebo effect identically zero for males (Han 2015 Results: 'placebo effect was acknowledged only in female subjects').",
      source_name        = "SEX"
    ),
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on the female placebo response: `p_max = (p_fem + etap_fem) * SEXF * (BMI/30.1)^bex` with reference BMI = 30.1 kg/m^2 (the cohort overall mean per Han 2015 Table 1) and estimated exponent bex = -4.74 (Table 2). Because bex < 0, less-obese females (BMI < 30.1) have a larger placebo effect than more-obese females, matching Han 2015 Results: 'placebo effect was more pronounced in female and relatively less obese patients'. Assumed time-fixed at baseline value.",
      source_name        = "BMI"
    ),
    DOSE = list(
      description        = "Current daily oral sibutramine-base dose (time-varying per subject)",
      units              = "mg (sibutramine base equivalent)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Use case (b) of the canonical DOSE: time-varying current daily dose driving a derived exposure metric inside the PD layer. Han 2015's published PD model drives the body-weight ODE from the steady-state sum AUC of the two active metabolites, AUC_ss,sum = AUC_M1 + AUC_M2 = Dose/CL_M1 + Dose/CL_M2 (paper equation 1, assuming complete sibutramine -> M1 -> M2 conversion). nlmixr2lib reproduces this by deriving AUC_ss,sum inside model() from the user-supplied current daily dose DOSE together with the individual M1 and M2 clearances (in L/hour, derived from the day-based cl and cl_m2 by dividing by 24). Per-subject daily dose in Han 2015: 0 mg/day for placebo; 8.37 mg/day sibutramine base initially in the active arm, escalated to 12.55 mg/day at week 4 if weight loss was less than 2 kg. Update DOSE over time in the user dataset to encode dose escalations. The DOSE covariate is supplied alongside (not in place of) any rxode2 AMT dose events that drive the depot compartment for PK simulation; the two channels are deliberately decoupled so the PD AUC_ss,sum reflects the steady-state assumption underlying Han 2015's exposure-response analysis.",
      source_name        = "DOSE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 120L,
    n_studies      = 2L,
    age_range      = "18-65 years (inclusion criterion); enrolled cohort 38.7 +/- 8.39 years (Table 1)",
    age_median     = "38.7 years (overall mean per Table 1)",
    weight_range   = "82.2 +/- 12.11 kg overall (Table 1)",
    weight_median  = "82.2 kg (overall mean per Table 1)",
    sex_female_pct = 67.5,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Korean abdominally obese adults (waist circumference >= 90 cm in men or >= 85 cm in women per the Korean Society for the Study of Obesity) with metabolic syndrome per the adult treatment panel (ATP) III definition, without clinically significant hypertension, diabetes, or any other underlying disease causing obesity. Subjects taking any drugs that could affect body weight or sibutramine PK/PD were excluded.",
    dose_range     = "Oral sibutramine 8.37 mg base/day (= 11.51 mg sibutramine mesylate; same active content as Reductil 10 mg = 10 mg sibutramine hydrochloride monohydrate) for the first 4 weeks; escalated to 12.55 mg base/day (= 17.26 mg sibutramine mesylate, equivalent to 15 mg sibutramine hydrochloride monohydrate) for weeks 4-24 if weight loss was less than 2 kg at week 4. The placebo arm followed the same titration scheme with matching placebo tablets.",
    regions        = "Republic of Korea (eight hospitals: Seoul St Mary's, Yeouido St Mary's, St Paul's, Uijeongbu St Mary's, St Vincent's, Daejeon St Mary's, Bucheon St Mary's, Incheon St Mary's).",
    notes          = "120 abdominally obese adults with metabolic syndrome were enrolled (60 sibutramine, 57 placebo, 3 excluded from PK; Table 1). To address the sparseness of patient PK samples (sparse 4-point sampling at predose, week 4 single sample, week 8 single sample, and week 24 predose), data from a separate full-PK study in 16 healthy young Korean male subjects (208 observations per metabolite) were merged into the PK dataset, giving 422 patient + 416 healthy = 838 PK observations total. A patient-vs-healthy indicator (ISP) was tested as a covariate and was not significant on any PK or PD parameter; ISP is therefore not included in this model file. Body-weight measurements (the PD endpoint) come from the patient cohort only (every 4 weeks for 24 weeks, n = 120 subjects). All subjects maintained their usual exercise and physical activity; dietary guidelines of 500-600 kcal/day deficit were provided and dietary intake was recorded by 24-hour recall. Cohort details per Han 2015 Table 1 and Study procedures section."
  )

  ini({
    # PK structural parameters. Time unit = day. Paper reports per-hour values
    # for rates and clearances; converted to per-day here by multiplying by 24.
    # Volumes are reported in L and unchanged. Concentrations are in ng/mL.
    lka     <- log(0.348 * 24)   ; label("M1 absorption rate (1/day; sibutramine absorbed and converted to M1)")  # Han 2015 Table 2: k_a = 0.348 1/h (paper Methods Figure 1 caption: ka represents absorption AND metabolism of sibutramine)
    lcl     <- log(158 * 24)     ; label("M1 metabolic clearance to M2 at reference age 35 years (L/day)")  # Han 2015 Table 2: CL_M1,t = 158 L/h
    lvc     <- log(2340)         ; label("M1 central volume of distribution (L)")  # Han 2015 Table 2: V_M1,c = 2,340 L
    lq      <- log(157 * 24)     ; label("M1 inter-compartmental clearance (L/day)")  # Han 2015 Table 2: Q_M1 = 157 L/h
    lvp     <- log(2060)         ; label("M1 peripheral volume of distribution (L)")  # Han 2015 Table 2: V_M1,p = 2,060 L
    lcl_m2  <- log(70.7 * 24)    ; label("M2 elimination clearance (L/day)")  # Han 2015 Table 2: CL_M2 = 70.7 L/h
    lvc_m2  <- log(43.9)         ; label("M2 central volume of distribution (L)")  # Han 2015 Table 2: V_M2 = 43.9 L

    # PK covariate effect (linear deviation around reference age 35 years).
    e_age_cl <- 0.0120           ; label("Linear age effect on M1 clearance (1/year; reference age 35 years)")  # Han 2015 Table 2: C_AGE = 0.0120 yr^-1; Methods: CL_M1 = CL_M1,t * (1 + C_AGE * (AGE - 35))

    # PD structural parameters. Time unit = day. BASE-related parameters in kg;
    # k_out in 1/day; placebo and drug-effect parameters in their published
    # units. The placebo coefficient `p_fem` is kept on linear (un-logged) scale
    # to match the source paper's additive IIV structure on the female placebo
    # coefficient (`p_fem + eta`); bex is left on linear scale because the
    # estimate is negative (Table 2 BEX = -4.74).
    lmw     <- log(89.1)         ; label("Typical baseline body weight in males (kg)")  # Han 2015 Table 2: MW = 89.1 kg (Baseline body weight; BASE = MW - SEX * FWC)
    lfwc    <- log(11.4)         ; label("Female weight correction subtracted from MW for females (kg)")  # Han 2015 Table 2: FWC = 11.4 kg
    lk_out  <- log(0.00947)      ; label("Body-weight loss rate constant (1/day)")  # Han 2015 Table 2: k_out = 0.00947 day^-1 (asymptotic-weight-loss half-life = ln(2)/k_out = 73.2 days; Discussion)
    p_fem   <- 0.0327            ; label("Placebo-effect coefficient in females at reference BMI 30.1 (fraction of k_in inhibited)")  # Han 2015 Table 2: P_fem = 0.0327
    bex     <- -4.74             ; label("Exponent of mean-normalised baseline BMI on placebo effect (unitless)")  # Han 2015 Table 2: BEX = -4.74; negative exponent gives larger placebo effect at lower BMI (Discussion)
    lemax   <- log(0.0735)       ; label("Maximal drug inhibition of weight gain (fraction)")  # Han 2015 Table 2: E_max = 0.0735 (paper Discussion: 'Empirical maximal efficacy was 7.35%')
    lauc50  <- log(106)          ; label("Sum-AUC at half-maximal drug inhibition AUC_50 (h*ng/mL)")  # Han 2015 Table 2: AUC_50 = 106 h*ng/mL

    # Inter-individual variability. Han 2015 Methods: BSV applied exponentially
    # (P_i = P_pop * exp(eta_i)) for the PK parameters and for BASE; additive
    # BSV (P_fem,i = P_fem + eta_i, restricted to females by the SEXF gate)
    # for the female placebo coefficient. The Table 2 column is reported as
    # CV% (between-subject variability); the log-normal variance is
    # omega^2 = log(1 + (CV/100)^2). For the additive eta on p_fem the
    # interpretation is the NONMEM convention CV% = 100 * sigma / theta, so
    # sigma = (CV/100) * theta and omega^2 = sigma^2.
    etalvc    ~ log(1 + 1.275^2)        # Han 2015 Table 2: V_M1,c CV 127.5 %
    etalvp    ~ log(1 + 0.778^2)        # Han 2015 Table 2: V_M1,p CV 77.8 %
    etalcl    ~ log(1 + 0.478^2)        # Han 2015 Table 2: CL_M1,t CV 47.8 %
    etalq     ~ log(1 + 0.820^2)        # Han 2015 Table 2: Q_M1 CV 82.0 %
    etalcl_m2 ~ log(1 + 0.296^2)        # Han 2015 Table 2: CL_M2 CV 29.6 %
    etalka    ~ log(1 + 1.172^2)        # Han 2015 Table 2: k_a CV 117.2 %
    etalmw    ~ log(1 + 0.125^2)        # Han 2015 Table 2: BASE (MW) CV 12.5 %; eta is applied to log MW so the same eta is shared by males (BASE = MW * exp(eta)) and females (BASE = MW * exp(eta) - FWC) per the source equation BASE = MW - SEX*FWC
    etap_fem  ~ (0.0435 * 0.0327)^2     # Han 2015 Table 2: P_fem additive BSV CV 4.35 % (NONMEM CV% = 100 * sigma / theta -> sigma = 0.0435 * 0.0327 = 1.42e-3; omega^2 = sigma^2 = 2.02e-6)

    # Residual error. Han 2015 Results: 'For both metabolites, only a
    # proportional error model was chosen'. Table 2 reports a single
    # proportional residual-variance estimate sigma^2_M1,p = 0.296; the
    # variance for M2 is not given separately and the M1 estimate is reused
    # here for the M2 observation; see vignette Assumptions and deviations.
    # BW residual error is additive (variance = 1.21 kg^2; SD = 1.1 kg).
    propSd    <- sqrt(0.296)     ; label("Proportional residual error on M1 plasma concentration (CV fraction)")  # Han 2015 Table 2: sigma^2_M1,p = 0.296 (so SD on linear-space proportional residual is sqrt(0.296) = 0.544)
    propSd_m2 <- sqrt(0.296)     ; label("Proportional residual error on M2 plasma concentration (CV fraction; paper-shared with M1)")  # Han 2015 Results sentence 'For both metabolites, only a proportional error model was chosen': Table 2 only reports one proportional sigma^2 = 0.296, reused for M2 here; vignette Assumptions documents the deviation
    addSd_BW  <- sqrt(1.21)      ; label("Additive residual error on body weight (kg)")  # Han 2015 Table 2: sigma^2 (PD) = 1.21 kg^2; additive residual on BW
  })

  model({
    # Individual PK parameters. Time axis is days; the age covariate is a
    # linear deviation around 35 years applied to CL_M1.
    age_cl <- 1 + e_age_cl * (AGE - 35)
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl + etalcl) * age_cl
    vc     <- exp(lvc + etalvc)
    q      <- exp(lq  + etalq)
    vp     <- exp(lvp + etalvp)
    cl_m2  <- exp(lcl_m2 + etalcl_m2)
    vc_m2  <- exp(lvc_m2)   # no IIV per Table 2 (ne)

    # PK ODEs. Sibutramine is the orally administered prodrug; the paper
    # assumes it converts entirely to the active mono-desmethyl metabolite M1
    # during absorption ('Instead of sibutramine, the inactive prodrug M1 was
    # assumed to be given orally and converted to M2 thereafter'). M1 is in
    # turn metabolised entirely to the di-desmethyl metabolite M2 (no separate
    # M1 excretion pathway); M2 is then eliminated via cl_m2. The shared
    # absorption-plus-conversion rate constant is `ka`. Dose `depot` should be
    # specified in mg of sibutramine base (matching DOSE_MG below).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          cl * central / vc -
                          q  * central / vc + q * peripheral1 / vp
    d/dt(peripheral1) <-  q  * central / vc - q * peripheral1 / vp
    d/dt(central_m2)  <-  cl * central / vc - cl_m2 * central_m2 / vc_m2

    # Concentration unit-conversion factor: dose is in mg and volumes in L, so
    # central/vc yields mg/L; multiply by 1000 to report Cc in ng/mL matching
    # the paper's reported concentration units (LLOQ 0.05 ng/mL for M1 and
    # 0.1 ng/mL for M2; Han 2015 Methods, Plasma concentration measurements).
    Cc    <- 1000 * central    / vc
    Cc_m2 <- 1000 * central_m2 / vc_m2

    # PD layer: body-weight ODE driven by AUC_ss,sum, the sum of the
    # steady-state daily AUCs of M1 and M2. Han 2015 equation 1:
    #   AUC_ss,sum = Dose_M1/CL_M1 + Dose_M2/CL_M2,
    # with full sibutramine -> M1 -> M2 conversion making Dose_M1 = Dose_M2 =
    # the administered daily sibutramine-base dose. cl and cl_m2 in this model
    # are in L/day; converted to L/h here so the AUC formula yields h*ng/mL
    # (DOSE_MG [mg] * 1e6 ng/mg / CL [L/h] / 1e3 mL/L = DOSE_MG * 1000 / CL).
    cl_h      <- cl    / 24
    cl_m2_h   <- cl_m2 / 24
    AUC_sssum <- 1000 * DOSE / cl_h + 1000 * DOSE / cl_m2_h

    # Individual PD parameters. The single IIV on log MW drives the
    # exponential between-subject spread of the typical male baseline weight;
    # for females the same eta shifts MW in log-space before the additive
    # female-weight correction FWC is subtracted.
    mw    <- exp(lmw + etalmw)
    fwc   <- exp(lfwc)
    base  <- mw - SEXF * fwc
    k_out <- exp(lk_out)
    k_in  <- k_out * base

    # Placebo. Han 2015 Results: 'placebo effect was acknowledged only in
    # female subjects'; SEXF gates the term to zero in males. The additive
    # eta is applied to p_fem so that the additive variability lives only on
    # the female-side coefficient (consistent with the paper's structure).
    p_fem_i <- p_fem + etap_fem
    p_max   <- p_fem_i * SEXF * (BMI / 30.1)^bex

    # Drug effect (sigmoidal Emax on AUC_ss,sum).
    emax   <- exp(lemax)
    auc_50 <- exp(lauc50)
    e_drug <- emax * AUC_sssum / (auc_50 + AUC_sssum)

    # Body-weight ODE (Han 2015 equations 2-5):
    #   dBW/dt = k_in * (1 - E_drug - P_max) - k_out * BW,
    # with k_in = k_out * BASE so that dBW/dt = 0 at baseline (BW = BASE) in
    # the absence of drug and placebo. The asymptotic-weight-loss half-life
    # is ln(2) / k_out ~ 73.2 days (Discussion).
    d/dt(bw) <- k_in * (1 - e_drug - p_max) - k_out * bw
    bw(0)    <- base

    BW <- bw

    # Observations and residual error.
    Cc    ~ prop(propSd)
    Cc_m2 ~ prop(propSd_m2)
    BW    ~ add(addSd_BW)
  })
}
