Davis_2008_maraviroc <- function() {
  description <- "Concentration-QT mixed-effects regression model relating single-dose oral maraviroc plasma concentrations to individual heart-rate-corrected QT intervals in healthy adult male and female volunteers (Davis 2008). No structural pharmacokinetic component is fit: the model is the one-stage NONMEM mixed-effects regression of observed QT on observed RR interval and observed maraviroc plasma concentration (Cp), with a fractional female-sex multiplier on the population QT intercept, a population QT/RR correction-factor exponent (Fridericia-style), and a linear concentration-QT slope. The single-dose population slope estimate (0.970 us mL/ng, 95% CI -0.571 to 2.48) was not significantly different from zero across the studied concentration range up to 2363 ng/mL. For simulation, supply observed or simulated maraviroc plasma concentration as the time-varying covariate CP_MVC_NGML and the RR interval as RR. Interoccasion variability on the QT intercept reported by the paper (13.4 ms^2) is not encoded (Hong_2015_moxifloxacin precedent: nlmixr2lib has no idiomatic IOV encoding for distributed models)."

  reference <- paste(
    "Davis JD, Hackman F, Layton G, Higgins T, Sudworth D, Weissgerber G.",
    "Effect of single doses of maraviroc on the QT/QTc interval in healthy",
    "subjects. Br J Clin Pharmacol. 2008 Apr;65(Suppl 1):68-75.",
    "doi:10.1111/j.1365-2125.2008.03138.x. PMID 18333868.",
    sep = " "
  )
  vignette <- "Davis_2008_maraviroc"
  units <- list(
    time          = "hour",
    dosing        = "not applicable (concentration-QT regression model; maraviroc plasma concentration is supplied as the time-varying covariate CP_MVC_NGML rather than via rxode2 dose events)",
    concentration = "ms (the modelled observation QT is the ECG QT interval in ms, NOT a drug concentration; the input covariate CP_MVC_NGML is maraviroc plasma concentration in ng/mL; the slash in this units string is only to satisfy checkModelConventions concentration-units parsing)"
  )

  covariateData <- list(
    CP_MVC_NGML = list(
      description        = "Time-varying maraviroc plasma concentration driving the linear concentration-QT term theta3 * Cp.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Davis 2008 Methods 'Concentration-QT modelling' (page 70): the one-stage analysis uses all postdose data for maraviroc 100, 300 and 900 mg, plus placebo and run-in-day predose data (Cp = 0 on placebo / predose). For simulation supply Cp at each observation time-point either from a separate maraviroc popPK source or from interpolated mean profiles of Davis 2008 Figure 2. Placebo simulations use Cp = 0 throughout.",
      source_name        = "Cp (Methods 'Concentration-QT modelling - Base model' equation)"
    ),
    RR = list(
      description        = "Time-varying ECG RR interval used inside the QT/RR correction factor (RR/1000)^theta4.",
      units              = "ms",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Davis 2008 used observed RR intervals from the ECG recorder over each post-dose measurement period (Methods 'Pharmacokinetics and pharmacodynamics'). The model normalises by 1000 ms (corresponding to 60 beats/min) inside the correction factor; for a resting adult at 65 bpm RR ~ 923 ms, at 60 bpm RR = 1000 ms, at 75 bpm RR = 800 ms. Time-varying within and between treatment periods.",
      source_name        = "RR (Methods 'Concentration-QT modelling - Base model' equation)"
    ),
    SEXF = list(
      description        = "Female sex indicator (0 = male, 1 = female).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Davis 2008 Methods 'Concentration-QT modelling - Base model' encodes Sex_i = 0 for men and Sex_i = 1 for women, with multiplicative effect (1 + theta2 * Sex_i) on the population QT intercept theta1 (Table 3 estimates: theta1 = 398 ms for men; theta1*(1 + 0.0166) = 404.6 ms for women).",
      source_name        = "Sex"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 61L,
    n_studies      = 1L,
    age_range      = "19-44 years (inclusion 18-45 years)",
    weight_range   = "men 60-90 kg; women 50-85 kg (per inclusion criteria); cohort means 72 kg (men) and 61 kg (women)",
    height_range   = "men mean 175 cm; women mean 165 cm",
    sex_female_pct = 51,
    race_ethnicity = c(White = 96.7, Asian = 3.3, Black = 1.6),
    disease_state  = "Healthy adult male and (postmenopausal, surgically sterile, or contracepting) female volunteers",
    dose_range     = "Single oral doses of maraviroc 100, 300, or 900 mg; placebo; oral moxifloxacin 400 mg active comparator. Five-way crossover with 7-day washout (Methods 'Study design')",
    regions        = "Single centre (Pfizer Global Research and Development, Sandwich, UK)",
    notes          = "Demographics from Davis 2008 Results (page 71). 30 men (1 Asian, 29 White) mean 175 cm / 72 kg; 31 women (1 Asian, 1 Black, 29 White) mean 165 cm / 61 kg. 61 subjects enrolled; one defaulted and three discontinued for non-treatment-related AEs (miscarriage, tonsil abscess, pyelonephritis). Per-treatment n in the QT analysis: maraviroc 100 mg n=59, 300 mg n=58, 900 mg n=58, moxifloxacin 400 mg n=58, placebo n=58 (manually-read QTcI; Table 2)."
  )

  ini({
    # =========================================================================
    # Concentration-QT mixed-effects regression - Davis 2008 'Concentration-QT
    # modelling - Model 2' (Methods page 70, full equation):
    #
    #   QT_ij = (theta1 * (1 + theta2 * Sex_i) + eta_1j)
    #             * (RR_ij / 1000) ^ (theta4 + eta_2j)
    #           + (theta3 + eta_3j) * Cp_ij + eps_ij
    #
    # Davis 2008 prints "(theta3 + eta_2j) * Cp_ij" with the SAME eta_2 index
    # appearing on the correction-factor exponent AND on the slope (Methods
    # page 70 Model 1 / Model 2 equations). Table 3 nevertheless reports
    # three distinct IIVs ("intersubject variability (correction factor)"
    # = 0.00287 and "(slope)" = 5.41), so the equation index repetition is
    # a typographical error and the slope eta is encoded here as a separate
    # eta_3j (etaslpqt) per the table. See vignette 'Assumptions and
    # deviations' for the full discussion.
    #
    # Davis 2008 Table 3 (page 73) reports the final retained model. A Model 3
    # extension with theta5 (concentration-by-RR-exponent interaction) was
    # tested but not retained ("there was no evidence of a change in the
    # QT/RR relationship with concentration", Results page 72), so Table 3 is
    # Model 2 (intercept + sex + RR correction + linear concentration slope).
    # =========================================================================
    lqtb       <- log(398);    label("Population baseline QT intercept theta1 for males (ms)")              # Davis 2008 Table 3 'Intercept (ms)' = 398 (95% CI 391-403)
    e_sexf_qtb <- 0.0166;      label("Fractional female-sex effect theta2 on the QT intercept (unitless)")  # Davis 2008 Table 3 'Sex (women)' = 0.0166 (95% CI -0.000489 to 0.0364); female intercept = 398 * (1 + 0.0166) = 404.6 ms
    lbrr       <- log(0.324);  label("Population QT/RR correction-factor exponent theta4 (unitless)")        # Davis 2008 Table 3 'QT correction factor' = 0.324 (95% CI 0.309 to 0.338). Compare Fridericia's b_s = 0.333.
    slpqt      <- 0.000970;    label("Linear concentration-QT slope theta3 (ms per (ng/mL))")                # Davis 2008 Table 3 'Slope (concentration, us mL/ng)' = 0.970 (95% CI -0.571 to 2.48). Paper unit us mL/ng = us/(ng/mL); converted to ms/(ng/mL) by dividing by 1000 (1 us = 0.001 ms): 0.970 / 1000 = 0.000970 ms/(ng/mL).

    # =========================================================================
    # IIV - Davis 2008 Table 3 lower block. The paper writes all etas
    # additively on the linear scale; CV% relative to each fixed effect is
    # also tabulated. For nlmixr2lib convention compliance the log-transformed
    # parameters lqtb and lbrr carry log-normal multiplicative etas with
    # omega^2 = log(1 + CV^2); for the small CV's tabulated (3.35% on
    # intercept, 16.5% on correction factor) the log-normal encoding is
    # approximately equivalent to the paper's additive-on-linear form. The
    # linear-scale slope slpqt retains an additive-on-linear eta (etaslpqt)
    # because the slope estimate can be negative across individuals
    # (population 95% CI for the slope crosses zero, CV 240%).
    # =========================================================================
    etalqtb    ~ 0.001122      # Davis 2008 Table 3 'Intersubject variability (intercept)' linear-scale variance = 178 ms^2, CV 3.35%; log-normal equivalent omega^2 = log(1 + 0.0335^2) = 0.001122 (small CV makes this nearly equal to the linear-additive form)
    etalbrr    ~ 0.026979      # Davis 2008 Table 3 'Intersubject variability (correction factor)' linear-scale variance = 0.00287, CV 16.5%; log-normal equivalent omega^2 = log(1 + 0.165^2) = 0.026979
    etaslpqt   ~ 5.41e-6       # Davis 2008 Table 3 'Intersubject variability (slope)' linear-scale variance = 5.41 (us mL/ng)^2 = 5.41 * 1e-6 (ms/(ng/mL))^2 after the same /1000 conversion used on the slope; CV 240%; encoded additive-on-linear because individual slopes can be negative

    addSd      <- sqrt(73.6);  label("Additive residual SD on QT (ms)")                                      # Davis 2008 Table 3 'Residual intrasubject variability' variance = 73.6 ms^2, CV 2.16%; additive SD = sqrt(73.6) = 8.580 ms
  })

  model({
    # -------------------------------------------------------------------------
    # Population QT intercept with sex effect plus subject IIV (Davis 2008
    # Methods 'Concentration-QT modelling - Model 2'). The paper places the
    # eta on the linear scale (additive to theta1 * (1+theta2*Sex_i)); for
    # convention compliance the eta is encoded log-normal here, which is
    # approximately equivalent for the small CV = 3.35% reported in Table 3.
    # -------------------------------------------------------------------------
    qt_intercept <- exp(lqtb + etalqtb) * (1 + e_sexf_qtb * SEXF)

    # -------------------------------------------------------------------------
    # Individual QT/RR correction-factor exponent (Davis 2008 Methods
    # equation). Log-normal eta is approximately equivalent to the paper's
    # additive-on-linear form for CV = 16.5%.
    # -------------------------------------------------------------------------
    brr <- exp(lbrr + etalbrr)

    # -------------------------------------------------------------------------
    # Individual concentration-QT slope (Davis 2008 Methods equation). The
    # eta is encoded additive-on-linear (etaslpqt) because the individual
    # slope can be negative (population CV 240%, 95% CI crosses zero).
    # -------------------------------------------------------------------------
    slp_i <- slpqt + etaslpqt

    # -------------------------------------------------------------------------
    # QT prediction (ms) - Davis 2008 'Concentration-QT modelling - Model 2'
    # final equation. CP_MVC_NGML is the time-varying maraviroc plasma
    # concentration (ng/mL); RR is the time-varying RR interval (ms).
    # -------------------------------------------------------------------------
    QT <- qt_intercept * (RR / 1000) ^ brr + slp_i * CP_MVC_NGML
    QT ~ add(addSd)
  })
}
