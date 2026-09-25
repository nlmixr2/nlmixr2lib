Stott_2018_amphotericinB_deoxycholate <- function() {
  description <- "Two-compartment IV-infusion population PK model for amphotericin B deoxycholate in adults with HIV-associated cryptococcal meningitis, with body weight entering clearance and central volume as additive linear intercept-plus-slope terms; disposition is written with explicit k12 / k21 micro-constants and no q / vp pair, so solve it with rxSolve(useLinCmt = FALSE) or the peripheral compartment is silently discarded (Stott 2018)"
  reference <- "Stott KE, Beardsley J, Whalley S, Kibengo FM, Mai NTH, Tung NLN, Cuc NTK, Kolamunnage-Dona R, Hope W, Day J. Population pharmacokinetic model and meta-analysis of outcomes of amphotericin B deoxycholate use in adults with cryptococcal meningitis. Antimicrob Agents Chemother. 2018;62(7):e02526-17. doi:10.1128/AAC.02526-17"
  vignette <- "Stott_2018_amphotericinB_deoxycholate"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Stott 2018 Materials and Methods, "Population
  # pharmacokinetic modeling": X(1) and X(2) are the amounts of amphotericin B
  # in milligrams in the central and peripheral compartments; the assay
  # measured amphotericin B in plasma by HPLC.
  compartmentData <- list(
    central = list(analyte = "amphotericin B", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "amphotericin B", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight at enrollment",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Stott 2018 equations e and f: SCL = Int_c + (Wt * Sl_c) and V = Int_v + (Wt * Sl_v). The weight effect is an ADDITIVE LINEAR term on the natural scale with a separately parameterised intercept and slope, NOT an allometric power term and NOT centred on a reference weight; the intercept is therefore the clearance (or volume) extrapolated to zero body weight and is not a typical-patient value. Time-fixed at the subject level. Cohort median 48 kg, range 32-68 kg (Stott 2018 Table 1, Combined column).",
      source_name = "Wt"
    )
  )

  # Screened during covariate model building but NOT retained in the final
  # model (model 2). Documented here so the paper's covariate screen is
  # preserved without declaring covariates that model() never references.
  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate by the Cockcroft-Gault equation",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Stott 2018 Results, 'Population pharmacokinetic models': univariate linear regression of the model-1 Bayesian posteriors showed a positive association of eGFR with clearance (slope 0.01, 95% CI 0 to 0.02) and volume (slope 0.67, 95% CI 0.36 to 0.98), so model 4 added eGFR to clearance as SCL = Int_c + (Wt * Sl_c) * (eGFR / med_eGFR) (equation g). Model 4 neither raised the log likelihood nor lowered the AIC relative to model 2 (log likelihood -43.1 vs -42.8; AIC 102.7 vs 101.9) and was rejected. Cohort median 76.7 mL/min/1.73 m^2, range 35.4-146.7 (Stott 2018 Table 1).",
      source_name = "eGFR"
    ),
    RACE_ASIAN = list(
      description = "Vietnamese cohort indicator, 1 = recruited in Ho Chi Minh City, Vietnam",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Ugandan cohort)",
      notes = "Stott 2018 Results: Bayesian posterior clearance 2.03 L/h (95% CI 1.69 to 2.38) in Vietnamese vs 2.24 L/h (1.91 to 2.56) in Ugandan patients (P = 0.37, Student t test), and volume 33.55 L (17.96 to 49.13) vs 63.93 L (40.98 to 86.88) (P = 0.09, Mann-Whitney). Ethnicity was not incorporated in the final model. 22 of 42 patients were Vietnamese.",
      source_name = "ethnicity"
    ),
    RACE_BLACK = list(
      description = "Ugandan cohort indicator, 1 = recruited at Masaka General Hospital, Uganda",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Vietnamese cohort)",
      notes = "Complement of RACE_ASIAN above; 20 of 42 patients were Ugandan. Screened and rejected together with RACE_ASIAN (Stott 2018 Results and Materials and Methods, 'Population pharmacokinetic modeling').",
      source_name = "ethnicity"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 42L,
    n_studies = 2L,
    age_range = "20-73 years",
    age_median = "33 years (mean 36)",
    weight_range = "32-68 kg",
    weight_median = "48 kg (mean 48)",
    sex_female_pct = 52.4,
    race_ethnicity = c(Vietnamese = 52.4, Ugandan = 47.6),
    disease_state = "HIV-associated cryptococcal meningitis. Baseline body mass index median 18 kg/m^2 (range 12-25); serum creatinine median 69 umol/L (37-167); Cockcroft-Gault eGFR median 76.7 mL/min/1.73 m^2 (35.4-146.7). Patients with renal failure, pregnancy, gastrointestinal bleeding, more than 7 days of prior anticryptococcal therapy, or requiring corticosteroids were excluded.",
    dose_range = "Amphotericin B deoxycholate 1 mg/kg once daily by intravenous infusion over 5 to 6 h, plus fluconazole 800 mg/day. Two patients also received dexamethasone.",
    regions = "Vietnam (Hospital for Tropical Diseases, Ho Chi Minh City; n = 22) and Uganda (Masaka General Hospital; n = 20)",
    renal_function = "Median Cockcroft-Gault eGFR 76.7 mL/min/1.73 m^2 (range 35.4-146.7); patients in renal failure were excluded",
    notes = "Three patients were recruited from a randomised controlled trial of adjuvant dexamethasone (ISRCTN 59144167) and 39 from a subsequent prospective descriptive study with identical inclusion and exclusion criteria at the same sites (Stott 2018 Materials and Methods, 'Clinical pharmacokinetic studies'). 479 of 553 plasma concentrations were analysable (282 of 312 Vietnamese and 197 of 241 Ugandan; mean 11.4 samples per patient, range 6 to 18); 74 samples were excluded because the sampling time was not recorded. Sampling was predose then 1, 2, 4, 8, 12 and 24 h after the start of infusion (the first five patients additionally at 16 and 20 h), on treatment days 1 or 2 and day 7, with sparse paired samples at lumbar punctures up to 17 days. Amphotericin B was assayed in plasma by HPLC with UV detection, calibration range 0.05 to 8.0 mg/L, limit of quantitation 0.05 mg/L, coefficient of variation below 9.3 percent across the range. Fitting used the nonparametric adaptive grid (NPAG) algorithm of Pmetrics 1.5.0, with observations weighted by the inverse of the estimated assay variance."
  )

  ini({
    # Structural parameters: Stott 2018 Table 2, model 2 (the final model).
    # Pmetrics NPAG reports the MEAN and the standard deviation of each
    # parameter's nonparametric marginal distribution; the means are used as
    # the typical values because Stott 2018 Results states that "the mean
    # parameter values predicted the observed values better than the medians".
    lcl <- log(0.67)
    label("Clearance intercept of the linear CL-versus-weight relation (L/h)") # Stott 2018 Table 2 model 2, 'SCL intercept (liters/h)' mean = 0.67
    e_wt_cl <- 0.03
    label("Body-weight slope of the linear CL-versus-weight relation (L/h per kg)") # Stott 2018 Table 2 model 2, 'SCL slope (liters/h/kg)' mean = 0.03
    lvc <- log(1.76)
    label("Central volume intercept of the linear V-versus-weight relation (L)") # Stott 2018 Table 2 model 2, 'Vc intercept (liters)' mean = 1.76
    e_wt_vc <- 0.82
    label("Body-weight slope of the linear V-versus-weight relation (L per kg)") # Stott 2018 Table 2 model 2, 'Vc slope (liters/kg)' mean = 0.82
    lk12 <- log(5.36)
    label("Central-to-peripheral first-order rate constant (1/h)") # Stott 2018 Table 2 model 2, 'K12 (h-1)' mean = 5.36
    lk21 <- log(9.92)
    label("Peripheral-to-central first-order rate constant (1/h)") # Stott 2018 Table 2 model 2, 'K21 (h-1)' mean = 9.92

    # Inter-individual variability. Pmetrics NPAG estimates a discrete,
    # non-Gaussian joint distribution over the support points and reports only
    # its marginal mean and SD, so the nonparametric shape cannot be
    # reproduced here. Each marginal is approximated by a log-normal with the
    # same coefficient of variation: omega^2 = log(1 + (SD / mean)^2), with the
    # reported mean carried as the median of the log-normal. See the vignette
    # 'Assumptions and deviations' section.
    etalcl ~ 0.000223 # Table 2 model 2, 'SCL intercept' SD 0.01 on mean 0.67 -> log(1 + (0.01/0.67)^2)
    etae_wt_cl ~ 0.105361 # Table 2 model 2, 'SCL slope' SD 0.01 on mean 0.03 -> log(1 + (0.01/0.03)^2)
    etalvc ~ 0.429977 # Table 2 model 2, 'Vc intercept' SD 1.29 on mean 1.76 -> log(1 + (1.29/1.76)^2)
    etae_wt_vc ~ 0.668759 # Table 2 model 2, 'Vc slope' SD 0.80 on mean 0.82 -> log(1 + (0.80/0.82)^2)
    etalk12 ~ 0.951894 # Table 2 model 2, 'K12' SD 6.76 on mean 5.36 -> log(1 + (6.76/5.36)^2); the Abstract prints 6.67 for this SD, a digit transposition
    etalk21 ~ 0.928184 # Table 2 model 2, 'K21' SD 12.27 on mean 9.92 -> log(1 + (12.27/9.92)^2)

    # Residual error. Pmetrics weighted each observation by the inverse of the
    # estimated assay variance and Stott 2018 does not publish the resulting
    # error polynomial or its gamma / lambda scale factor, so no fitted
    # residual magnitude exists to transcribe. The value below is the
    # analytical assay imprecision reported in Materials and Methods,
    # "Measurement of amphotericin B concentrations", used here as a stand-in
    # so that stochastic simulation is possible; it is an assay-precision
    # figure and not a paper-derived residual-error estimate, and it excludes
    # model misspecification.
    propSd <- fixed(0.093)
    label("Proportional residual SD, taken from the reported HPLC assay imprecision (fraction)") # Stott 2018 Materials and Methods: 'The coefficient of variation was < 9.3% over the concentration range of 0.05 to 8 mg/liter'
  })
  model({
    # Individual parameters. Stott 2018 equations e and f parameterise
    # clearance and central volume as additive linear functions of weight with
    # separately estimated intercepts and slopes. In a nonparametric fit every
    # one of those four quantities is a subject-level parameter, so the
    # intercept and the slope each carry their own log-normal eta.
    wt_cl <- e_wt_cl * exp(etae_wt_cl)
    wt_vc <- e_wt_vc * exp(etae_wt_vc)
    cl <- exp(lcl + etalcl) + wt_cl * WT
    vc <- exp(lvc + etalvc) + wt_vc * WT
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    kel <- cl / vc

    # Stott 2018 equations a.2 and b. Dose enters central as a zero-order
    # intravenous infusion R(1), encoded on the event rows.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Stott 2018 equation c: Y(1) = X(1) / V. Amounts in mg and volumes in L
    # give mg/L, the units of the reported concentrations (assay calibration
    # range 0.05 to 8.0 mg/L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
