Abrantes_2020_octocogAlfa <- function() {
  description <- "Integrated population PK and repeated time-to-categorical-event (RTTCE) model for bleeds in children, adolescents and adults with severe hemophilia A on prophylactic octocog alfa (BAY 81-8973, Kovaltry) from the LEOPOLD trials (Abrantes 2020, all bleeds). FVIII activity follows a two-compartment model with lean body weight (power, centred at 51.1 kg) on CL and V1 and a random effect scaling the combined additive + proportional residual error. The bleeding hazard is Gompertz in years since study start with an inhibitory Imax effect of the predicted FVIII activity, h(t) = lam * exp(gam * (t - 1)) * (1 - C / (C + IF50)) * exp(eta), estimated as the hazards at 0.5 and 20 IU/dL one year after study start (2.82 and 1.00 per year; derived lam 2.96 per year, IF50 10.2 IU/dL, gam -0.566 per year; IIV 136% CV). The severity of each bleed (mild / moderate / severe) is a proportional-odds model with one shared logit-scale random effect. The cumulative hazard is an ODE state; the survivor function and the severity probabilities are exposed. This is the prior model of the Abrantes 2019 Bayesian bleed-forecasting study. Sister model: modellib('Abrantes_2020_octocogAlfa_jointBleed') (RTTCE re-estimated on joint spontaneous bleeds only)."
  reference <- paste(
    "Abrantes JA, Solms A, Garmann D, Nielsen EI, Jonsson S, Karlsson MO.",
    "Relationship between factor VIII activity, bleeds and individual",
    "characteristics in severe hemophilia A patients. Haematologica.",
    "2020;105(5):1443-1453. doi:10.3324/haematol.2019.217133. PMID:31371418.",
    "PMCID:PMC7193498. Parameter estimates in Table 2; hazard, lambda/IF50",
    "reparameterisation and proportional-odds severity model in the Online",
    "Supplementary Appendix (equations 1-3 and 6-12).",
    "Applied as the Bayesian prior in: Abrantes JA, Solms A, Garmann D,",
    "Nielsen EI, Jonsson S, Karlsson MO. Bayesian forecasting utilizing",
    "bleeding information to support dose individualization of factor VIII.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8(12):894-903.",
    "doi:10.1002/psp4.12464. PMID:31668021. PMCID:PMC6930854 (PK equations",
    "with the IIV and residual variances in the Supplemental Text; hazard in",
    "Equation 1).",
    sep = " "
  )
  vignette <- "Abrantes_2020_octocogAlfa"
  units <- list(time = "h", dosing = "IU", concentration = "IU/dL")

  # One random effect multiplies the whole residual error (Abrantes 2019
  # Supplemental Text equation 8: y = FVIII + (FVIII * eps1 + eps2) *
  # exp(eta_eps)) and one multiplies the whole bleeding hazard (Abrantes 2020
  # Supplementary Appendix equation 1); neither has an ini() typical value.
  paper_specific_etas <- c("etaresSd", "etalhaz")

  covariateData <- list(
    LBM = list(
      description = "Lean body weight at study start (canonical column LBM; the source uses LBW)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL (exponent 0.646) and V1 (0.958), centred at 51.1 kg (Abrantes 2020 Table 2 footnote a; Abrantes 2019 Supplemental Text equations 1-2), the reference inherited from the Garmann 2017 BAY 81-8973 model (modellib('Garmann_2017_BAY81_8973')). The body-composition formula is not stated in either paper. Cohort median 50.9 kg (range 9.25-79.2; Abrantes 2020 Table 1). Acts on the bleeding hazard only through the predicted FVIII activity.",
      source_name = "LBW"
    )
  )

  # Abrantes 2020 integrated these covariates through a full random effects
  # model (FREM): each covariate was an observation with its own random
  # effect, and parameter-covariate relationships came from the covariance
  # matrix, which is not published. Only correlation coefficients (Figure 2)
  # and graphical effect sizes (Figure 3) are reported, so none is encoded.
  # The IIV in ini() is the FREM marginal (no-covariate) variability.
  covariatesDataExcluded <- list(
    VWF = list(
      description = "Von Willebrand factor level at study start",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "FREM correlation with CL r = -0.54 (reduced unexplained IIV on CL by 4.7%); typical CL 2.30 dL/h at 64% (5th percentile) and 1.59 dL/h at 179% (95th percentile) (Abrantes 2020 Results, Figure 3). Also correlated with bleeding severity (r = 0.28). Missing for 55 patients (all of LEOPOLD kids); median 104% (range 43-242) (Table 1)."
    ),
    ABR_PRESTUDY = list(
      description = "Number of bleeds in the 12-month pre-study period",
      units = "bleeds/year",
      type = "continuous",
      reference_category = NULL,
      notes = "Strongest FREM relationship: r = 0.45 with the bleeding hazard (reduced unexplained IIV by 15%). Against a mean patient with 8.2 pre-study bleeds, 1 bleed (5th percentile) gives a 54% lower and 84 bleeds (95th percentile) a 147% higher hazard, i.e. lambda(0.5 IU/dL) 1.30 and 6.97 per year (Abrantes 2020 Results, Figure 3). Median 11 (range 0-106; Table 1)."
    ),
    RACE_BLACK = list(
      description = "Black race indicator: 1 = Black, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian is the Figure 3 comparator)",
      notes = "FREM correlation with CL r = -0.25; estimated CL 1.41 dL/h vs 1.93 dL/h in Caucasians (27% lower) (Abrantes 2020 Results and Discussion). 8 patients (4.7%; Table 1)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator: 1 = Asian, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian is the Figure 3 comparator)",
      notes = "FREM correlation with CL r = 0.17; estimated CL 2.21 dL/h vs 1.93 dL/h in Caucasians (14.2% higher). 23 patients (13%; Table 1)."
    ),
    AGE = list(
      description = "Age at study start",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "FREM correlation with bleeding severity r = 0.23 (Abrantes 2020 Results). Median 22 years (range 1-61; Table 1)."
    ),
    WT = list(
      description = "Body weight at study start",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the FREM; no relevant relationship reported. Median 60 kg (range 11-124; Table 1)."
    ),
    BMI = list(
      description = "Body mass index at study start",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the FREM; no relevant relationship reported. Median 20.3 (range 13.0-38.3; Table 1)."
    ),
    N_TARGET_JOINTS = list(
      description = "Number of target joints for bleeds at study start",
      units = "(count)",
      type = "count",
      reference_category = NULL,
      notes = "Screened in the FREM, alone and as the ratio of pre-study bleeds to target joints; no relevant relationship reported. Median 1 (range 1-5); missing for 61 patients (Table 1)."
    ),
    PRIOR_PROPHYLAXIS = list(
      description = "Treatment history before the study: 1 = prophylaxis, 0 = on demand",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (on demand)",
      notes = "Screened in the FREM; correlated with pre-study bleeds (r = -0.71) but not reported as relevant for any model parameter. 91 prophylaxis / 81 on demand (Table 1)."
    )
  )

  compartmentData <- list(
    central = list(analyte = "octocog alfa (FVIII activity)", units = "IU", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "octocog alfa (FVIII activity)", units = "IU", specimen = "plasma", verified = TRUE),
    cumhaz = list(
      analyte = "Cumulative hazard of a bleed (any severity)",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 183L,
    n_studies = 3L,
    age_range = "1-61 years",
    age_median = "22 years",
    weight_range = "11-124 kg",
    weight_median = "60 kg",
    sex_female_pct = NULL,
    race_ethnicity = c(White = 77, Black = 4.7, Asian = 13, Hispanic = 4.7, Unknown = 0.6),
    disease_state = "Previously treated patients with severe hemophilia A (endogenous FVIII activity < 1 IU/dL), no history of FVIII inhibitors.",
    dose_range = "Intravenous octocog alfa prophylaxis: 20-50 IU/kg 2-3 times/week (LEOPOLD I), 20-40 IU/kg 2-3 times/week (LEOPOLD II), 25-50 IU/kg at least twice weekly (LEOPOLD kids); single 50 IU/kg PK doses in LEOPOLD I. Median individual dose 37.9 IU/kg (range 4.30-199).",
    regions = "Multinational (LEOPOLD I, LEOPOLD II, LEOPOLD kids).",
    biomarkers = "1535 FVIII activity observations (chromogenic assay) from 183 patients; 633 bleeds in 172 patients with a bleeding observation period (116 with at least one bleed; median 2, range 0-33); 47% mild, 44% moderate, 9.8% severe; 61% spontaneous.",
    notes = "Post-hoc pooled analysis of LEOPOLD I (NCT01029340), LEOPOLD II (NCT01233258) and LEOPOLD kids (NCT01311648) (Abrantes 2020 Table 1, Table S1). Median bleeding observation period 12.0 months (range 3.08-13.1); 6 months in LEOPOLD kids. Eleven patients contributed PK only. PK and RTTCE components estimated simultaneously (IMPMAP, NONMEM 7.4.3); covariates added afterwards by FREM with the fixed effects fixed. The Abrantes 2019 forecasting study used the same 172 patients (51 aged < 12 years, 121 aged >= 12 years)."
  )

  ini({
    # ------------------------------------------------------------------
    # PK (Abrantes 2020 Table 2; Abrantes 2019 Supplemental Text eq. 1-7).
    # ------------------------------------------------------------------
    lcl <- log(1.93); label("Clearance at LBW 51.1 kg (dL/h)") # Abrantes 2020 Table 2 'CL, dL/h' = 1.93 (RSE 3.54%)
    lvc <- log(30.3); label("Central volume of distribution at LBW 51.1 kg (dL)") # Abrantes 2020 Table 2 'V1, dL' = 30.3 (RSE 1.72%)
    lq <- log(1.69); label("Intercompartmental clearance (dL/h)") # Abrantes 2020 Table 2 'Q, dL/h' = 1.69 (RSE 23.0%)
    lvp <- log(6.29); label("Peripheral volume of distribution (dL)") # Abrantes 2020 Table 2 'V2, dL' = 6.29 (RSE 8.79%)
    e_lbm_cl <- 0.646; label("Power exponent of lean body weight on CL, centred at 51.1 kg (unitless)") # Abrantes 2020 Table 2 'Effect of LBW on CL' = 0.646 (RSE 10.0%); 2019 Suppl. eq. 1 prints 0.65
    e_lbm_vc <- 0.958; label("Power exponent of lean body weight on V1, centred at 51.1 kg (unitless)") # Abrantes 2020 Table 2 'Effect of LBW on V1' = 0.958 (RSE 3.13%); 2019 Suppl. eq. 2 prints 0.96

    # IIV variances (Abrantes 2019 Supplemental Text); Abrantes 2020 Table 2
    # reports %CV = 100 * sqrt(omega^2).
    etalcl ~ 0.0913 # Abrantes 2019 Suppl. Text: variance 0.0913; Abrantes 2020 Table 2 IIV CL 30.2 %CV
    etalvc ~ 0.0226 # Abrantes 2019 Suppl. Text: variance 0.0226; Abrantes 2020 Table 2 IIV V1 15.1 %CV
    etaresSd ~ 0.398 # Abrantes 2019 Suppl. Text: variance 0.398; Abrantes 2020 Table 2 IIV residual error 63.1 %CV

    # Residual error, Abrantes 2020 Table 2 (the 2019 Supplemental Text
    # gives the variances 2.31 = 1.52^2 and 0.0411 = 0.203^2 with the
    # proportional / additive labels swapped).
    propSd <- 0.203; label("Proportional residual error (fraction)") # Abrantes 2020 Table 2 'Proportional, %CV' = 20.3 (RSE 5.52%)
    addSd <- 1.52; label("Additive residual error (IU/dL)") # Abrantes 2020 Table 2 'Additive, IU/dL' = 1.52 (RSE 11.0%)

    # ------------------------------------------------------------------
    # RTTCE, all bleeds (Abrantes 2020 Table 2; Supplementary Appendix
    # eq. 1-3 and 6-12). Hazards are per year; time in the Gompertz term is
    # years since study start.
    # ------------------------------------------------------------------
    llam05_haz <- log(2.82); label("Bleeding hazard at FVIII 0.5 IU/dL one year after study start (1/year)") # Abrantes 2020 Table 2 'lambda 0.5 IU/dL, year-1' = 2.82 (RSE 17.6%)
    llam20_haz <- log(1.00); label("Bleeding hazard at FVIII 20 IU/dL one year after study start (1/year)") # Abrantes 2020 Table 2 'lambda 20 IU/dL, year-1' = 1.00 (RSE 17.3%)
    gam_haz <- -0.566; label("Gompertz shape factor of the bleeding hazard (1/year)") # Abrantes 2020 Table 2 'gamma, year-1' = -0.566 (RSE 26.1%); Abrantes 2019 eq. 1

    # The source models the logits of P(moderate or severe) = b1 + eta and
    # P(severe) = b1 + b2 + eta. Stored in the library's cumulative-from-
    # the-bottom form: logit P(mild) = -b1 and logit P(mild or moderate) =
    # -(b1 + b2), so logite0 = -b1 and logitstep1 = -b2; the shared random
    # effect changes sign, which leaves its normal distribution unchanged.
    logite0 <- -0.424; label("Logit of the probability that a bleed is mild (unitless)") # Abrantes 2020 Table 2 'b1' = 0.424 (RSE 68.2%), sign reversed
    logitstep1 <- 3.43; label("Logit increment from P(mild) to P(mild or moderate) (unitless)") # Abrantes 2020 Table 2 'b2' = -3.43 (RSE 6.42%), sign reversed

    # Abrantes 2020 Table 2 IIV lambda 136 %CV; omega^2 = 1.36^2 on the
    # %CV = 100 * sqrt(omega^2) convention of the PK rows. The forecasting
    # paper quotes the same 136% CV (Abrantes 2019 eq. 1).
    etalhaz ~ 1.8496 # Abrantes 2020 Table 2 IIV lambda 136 %CV; 1.36^2 = 1.8496
    etalogite0 ~ 2.7889 # Abrantes 2020 Table 2 IIV 'b1 and b2, SD' = 1.67; 1.67^2 = 2.7889
  })

  model({
    # --- PK ---------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (LBM / 51.1)^e_lbm_cl
    vc <- exp(lvc + etalvc) * (LBM / 51.1)^e_lbm_vc
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    # --- bleeding hazard --------------------------------------------------
    lam05_haz <- exp(llam05_haz)
    lam20_haz <- exp(llam20_haz)
    # Supplementary Appendix eq. 2-3: hazard scale at zero FVIII one year
    # after study start, and the FVIII activity halving the hazard.
    lam_haz <- lam05_haz * lam20_haz * (0.5 - 20) / (lam05_haz * 0.5 - lam20_haz * 20)
    ic50 <- (lam05_haz * 0.5 - lam20_haz * 20) / (lam20_haz - lam05_haz)

    # Model time is hours from study start; the Gompertz term uses years.
    t_year <- t / 8760
    hazard_year <- lam_haz * exp(gam_haz * (t_year - 1)) * (1 - Cc / (Cc + ic50)) * exp(etalhaz)
    hazard <- hazard_year / 8760

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    # Probability of no bleed since the start of the integration window.
    sur <- exp(-cumhaz)

    # --- bleed severity (proportional odds; eq. 6-12) ---------------------
    logit_mild <- logite0 + etalogite0
    prob_bleedsev_mild <- expit(logit_mild)
    prob_bleedsev_le_moderate <- expit(logit_mild + logitstep1)
    prob_bleedsev_moderate <- prob_bleedsev_le_moderate - prob_bleedsev_mild
    prob_bleedsev_severe <- 1 - prob_bleedsev_le_moderate

    # Residual error: eps1 and eps2 independent, so their variances add, and
    # exp(eta) scales the total (2019 Supplemental Text eq. 8).
    resSd <- exp(etaresSd) * sqrt((propSd * Cc)^2 + addSd^2)
    Cc ~ add(resSd)
    # Severity of a bleed: category 1 = mild, 2 = moderate, 3 = severe.
    bleedsev ~ c(prob_bleedsev_mild, prob_bleedsev_moderate)
  })
}
