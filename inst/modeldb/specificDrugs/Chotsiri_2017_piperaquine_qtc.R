Chotsiri_2017_piperaquine_qtc <- function() {
  description <- paste(
    "Linear concentration-QTc model for piperaquine in 16 healthy Thai adult",
    "volunteers, from Chotsiri 2017. The pharmacodynamic endpoint is the",
    "double-delta-corrected prolongation of the individually rate-corrected QT",
    "interval (DeltaDeltaQTcI, ms): each subject's post-dose QTcI minus their",
    "own pre-dose baseline, minus the time-matched same quantity from the",
    "primaquine-alone (placebo) arm, which removes both heart-rate and",
    "circadian effects. The relationship is a direct-response line with no",
    "hysteresis, DeltaDeltaQTcI = 0 + eta + 0.0417 * Cc, so every 100 ng/mL of",
    "plasma piperaquine adds 4.17 ms. The typical baseline was estimated close",
    "to zero and fixed there, but between-subject variability was retained on",
    "it (SD 15.9 ms), which is what lets an individual maximum prolongation",
    "come out negative. Emax, power and effect-compartment forms were tested",
    "and the linear direct-response model was retained; age, sex, potassium,",
    "sodium and primaquine coadministration were screened on the relationship",
    "and none was retained. The three-compartment piperaquine population PK",
    "model of the same paper is embedded so the prediction is driven by",
    "simulated concentration, matching the paper's sequential fit in which",
    "individual pharmacokinetic parameter estimates were imputed directly into",
    "the pharmacodynamic model. Doses are expressed as piperaquine BASE.",
    "Sister model files from the same paper:",
    "modellib('Chotsiri_2017_piperaquine') (population PK) and",
    "modellib('Chotsiri_2017_dihydroartemisinin') (dihydroartemisinin",
    "population PK).",
    sep = " "
  )
  reference <- paste(
    "Chotsiri P, Wattanakul T, Hoglund RM, Hanboonkunupakarn B,",
    "Pukrittayakamee S, Blessborn D, Jittamala P, White NJ, Day NPJ,",
    "Tarning J.",
    "Population pharmacokinetics and electrocardiographic effects of",
    "dihydroartemisinin-piperaquine in healthy volunteers.",
    "Br J Clin Pharmacol. 2017;83(12):2752-2766.",
    "doi:10.1111/bcp.13372. PMC5698590.",
    "Open Access under CC BY-NC 4.0.",
    "The pharmacodynamic parameters are in Table 2 ('Pharmacodynamic",
    "parameters') and the model equation is Equation 8 in Methods,",
    "'Population cardiac electrophysiological pharmacodynamics of",
    "piperaquine'; the embedded pharmacokinetic parameters are in Table 2",
    "('Pharmacokinetic parameters of piperaquine'). Sister model files from",
    "the same paper: modellib('Chotsiri_2017_piperaquine') and",
    "modellib('Chotsiri_2017_dihydroartemisinin').",
    sep = " "
  )
  vignette <- "Chotsiri_2017_dihydroartemisinin_piperaquine"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "(the observation QTcI is the double-delta-corrected prolongation of the individually rate-corrected QT interval, ms; the driving plasma piperaquine concentration Cc is in ng/mL)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Acts on the EMBEDDED pharmacokinetic model only (fixed allometric",
        "exponent 0.75 on all clearances and 1 on all volumes, centred on the",
        "study-median 64 kg), and therefore reaches the QTc prediction only",
        "through piperaquine concentration. Body weight was not among the",
        "covariates screened on the pharmacodynamic parameters. Set WT = 64 to",
        "recover the tabulated typical values.",
        sep = " "
      ),
      source_name = "BW"
    ),
    OCC = list(
      description = paste(
        "Integer-valued occasion index identifying which of the two",
        "dihydroartemisinin-piperaquine dosing occasions a dose belongs to,",
        "for between-occasion variability on the relative bioavailability and",
        "mean transit time of the embedded pharmacokinetic model.",
        sep = " "
      ),
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Acts on the EMBEDDED pharmacokinetic model only. Values 1 and 2 are",
        "multiplexed inside model() onto etaiov_fdepot_1 / etaiov_fdepot_2",
        "and etaiov_mtt_1 / etaiov_mtt_2, each pair sharing one variance.",
        "OCC = 0, or any value outside 1..2, yields the occasion-free",
        "typical value. See",
        "modellib('Chotsiri_2017_piperaquine') for the full annotation.",
        sep = " "
      ),
      source_name = "OCC"
    )
  )

  # Covariates the source screened on the pharmacodynamic parameters and did
  # not retain. Documentation only -- none is referenced in model().
  covariatesDataExcluded <- list(
    CONMED_PRIMAQUINE = list(
      description = "Concomitant single low-dose primaquine administration indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (dihydroartemisinin-piperaquine alone)",
      notes = paste(
        "Screened on the exposure-response relationship and not retained:",
        "'Primaquine did not affect the relationship and no other significant",
        "covariates (age, gender and electrolyte levels) were identified in the",
        "stepwise covariate approach' (Results, 'Electrocardiographic effects",
        "of piperaquine'). The primaquine-alone arm serves as the placebo arm",
        "of the double-delta correction; the paper separately confirmed there",
        "was no concentration-response relationship between primaquine",
        "concentration and DeltaQTc. No coefficient is reported.",
        sep = " "
      )
    ),
    AGE = list(
      description = "Subject age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a linear covariate on the piperaquine-related DeltaDeltaQTc prolongation (Methods) and not retained. No coefficient is reported."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened as a linear covariate on the piperaquine-related DeltaDeltaQTc prolongation (Methods) and not retained. No coefficient is reported. (Sex was screened here even though it was excluded from the pharmacokinetic covariate search.)"
    ),
    POT = list(
      description = "Serum potassium at admission.",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a linear covariate on the piperaquine-related DeltaDeltaQTc prolongation (Methods) and not retained; the Discussion attributes this to the healthy-volunteer setting. Table 1 reports the cohort median as 4.25 mmol/L. No coefficient is reported."
    ),
    SOD = list(
      description = "Serum sodium at admission.",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a linear covariate on the piperaquine-related DeltaDeltaQTc prolongation (Methods) and not retained. Table 1 reports the cohort median as 140 mmol/L. No coefficient is reported."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 16L,
    n_studies = 1L,
    age_range = "22-53 years",
    age_median = "40 years",
    weight_range = "54.0-71.4 kg",
    weight_median = "64.1 kg",
    sex_female_pct = 68.75,
    race_ethnicity = c(Asian = 100),
    disease_state = "Healthy volunteers. Participants with malaria or glucose-6-phosphate dehydrogenase deficiency, and pregnant or lactating women, were excluded.",
    dose_range = "Single oral dose of three co-formulated dihydroartemisinin-piperaquine tablets (40 mg dihydroartemisinin + 320 mg piperaquine phosphate each), i.e. 960 mg piperaquine phosphate = 514 mg piperaquine base, given 30 min after a light meal (~200 kcal, 8 g fat).",
    regions = "Bangkok, Thailand",
    biomarkers = "12-lead ECG (ECG-1250K, Nihon Kohden) recorded at 10 mm/mV and 25 mm/s after at least 20 min of rest: twice before dosing and at 1, 2, 4, 8, 12 and 24 h after each administration. Any automatic readout above 450 ms was adjudicated by an unblinded research physician and a blinded cardiologist. Baseline uncorrected QT was 395 ms (370-446) and baseline QTc 422 ms (386-466) (Table 1).",
    notes = paste(
      "Open-label, randomized, three-way crossover (NCT01525511). Observed QT",
      "was corrected for heart rate three ways -- Bazett (exponent 1/2),",
      "Fridericia (exponent 1/3) and a data-driven individual exponent fitted",
      "by ordinary least squares to each subject's placebo-arm QT-RR pairs",
      "(Equation 5). The individual correction left the fewest subjects with a",
      "QTc-RR regression slope different from zero (5/16 versus 6/16 for both",
      "Bazett and Fridericia) and was carried forward. DeltaQTc is post-dose",
      "QTc minus baseline QTc (Equation 6) and DeltaDeltaQTc subtracts the",
      "time-matched DeltaQTc of the primaquine-alone arm (Equation 7). The",
      "concentration-response analysis found no relationship with DeltaQRS",
      "(P = 0.520), and DeltaJTc and DeltaQTc behaved almost identically.",
      "Estimated in NONMEM 7.3 with the Laplacian method. Eta shrinkage 26.5%,",
      "epsilon shrinkage 2.39%.",
      sep = " "
    )
  )

  ini({
    # ==================================================================
    # EMBEDDED PIPERAQUINE PK -- identical to
    # modellib('Chotsiri_2017_piperaquine'), Chotsiri 2017 Table 2
    # 'Pharmacokinetic parameters of piperaquine'. See that file for the
    # full annotation, including the transit-chain convention
    # (ka = ktr = 3 / MTT), the 64-kg allometric reference, and the
    # piperaquine-phosphate-to-base dose conversion. The residual error
    # of the concentration observation is deliberately absent here: this
    # model observes the electrocardiographic endpoint, and the
    # concentration is a latent driver.
    # ==================================================================
    lmtt <- log(3.13)
    label("Mean transit time through the two-transit-compartment absorption chain MTT (h)")
    # Chotsiri 2017 Table 2: MTT = 3.13 h (%RSE 9.42; bootstrap 95% CI 2.66-3.84)

    lcl <- log(27.4)
    label("Apparent elimination clearance CL/F at WT = 64 kg (L/h)")
    # Chotsiri 2017 Table 2: CL/F = 27.4 L/h (%RSE 5.50; bootstrap 95% CI 24.6-30.4)

    lvc <- log(751)
    label("Apparent central volume of distribution Vc/F at WT = 64 kg (L)")
    # Chotsiri 2017 Table 2: V_C/F = 751 L (%RSE 23.5; bootstrap 95% CI 470-1160)

    lq <- log(206)
    label("Apparent inter-compartmental clearance to peripheral1, Q1/F at WT = 64 kg (L/h)")
    # Chotsiri 2017 Table 2: Q_P1/F = 206 L/h (%RSE 9.56; bootstrap 95% CI 166-242)

    lvp <- log(1900)
    label("Apparent first peripheral volume of distribution Vp1/F at WT = 64 kg (L)")
    # Chotsiri 2017 Table 2: V_P1/F = 1900 L (%RSE 8.23; bootstrap 95% CI 1660-2260)

    lq2 <- log(71.5)
    label("Apparent inter-compartmental clearance to peripheral2, Q2/F at WT = 64 kg (L/h)")
    # Chotsiri 2017 Table 2: Q_P2/F = 71.5 L/h (%RSE 9.01; bootstrap 95% CI 58.5-84.4)

    lvp2 <- log(13500)
    label("Apparent second peripheral volume of distribution Vp2/F at WT = 64 kg (L)")
    # Chotsiri 2017 Table 2: V_P2/F = 13 500 L (%RSE 8.95; bootstrap 95% CI 11 400-16 000)

    lfdepot <- fixed(log(1))
    label("Relative bioavailability F for the typical subject (unitless)")
    # Chotsiri 2017 Table 2: F (%) = '100 Fixed'

    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F, Q1/F and Q2/F, referenced to 64 kg (unitless)")
    # Chotsiri 2017 Methods, Equation 3: clearance exponent 3/4

    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on Vc/F, Vp1/F and Vp2/F, referenced to 64 kg (unitless)")
    # Chotsiri 2017 Methods, Equation 4: volume exponent 1

    etalfdepot ~ 0.0315384
    # Chotsiri 2017 Table 2: F BSV = 17.9% CV (%RSE 34.0); log(0.179^2 + 1) = 0.0315384
    etalcl ~ 0.0118110
    # Chotsiri 2017 Table 2: CL/F BSV = 10.9% CV (%RSE 37.2); log(0.109^2 + 1) = 0.0118110
    etalvc ~ 0.1653246
    # Chotsiri 2017 Table 2: V_C/F BSV = 42.4% CV (%RSE 40.9); log(0.424^2 + 1) = 0.1653246
    etalq2 ~ 0.0564569
    # Chotsiri 2017 Table 2: Q_P2/F BSV = 24.1% CV (%RSE 36.3); log(0.241^2 + 1) = 0.0564569
    # Table 2's F row carries TWO variability lines for piperaquine --
    # a between-subject 17.9% and a between-occasion 19.1% -- and the
    # second is easily lost when the multi-line table cell is collapsed.
    etaiov_fdepot_1 ~ 0.0358313
    # Chotsiri 2017 Table 2: F BOV = 19.1% CV (%RSE 13.3); log(0.191^2 + 1) = 0.0358313
    etaiov_fdepot_2 ~ fixed(0.0358313)
    # Same published BOV variance as occasion 1
    etaiov_mtt_1 ~ 0.0986537
    # Chotsiri 2017 Table 2: MTT BOV = 32.2% CV (%RSE 13.4); log(0.322^2 + 1) = 0.0986537
    etaiov_mtt_2 ~ fixed(0.0986537)
    # Same published BOV variance as occasion 1

    # ==================================================================
    # CONCENTRATION-DeltaDeltaQTcI MODEL -- Chotsiri 2017 Table 2
    # 'Pharmacodynamic parameters' and Equation 8 in Methods:
    #
    #   DeltaDeltaQTc_i = theta_1 + eta_1 + theta_2 * CP + epsilon_i
    #
    # where theta_1 is the typical baseline prolongation, eta_1 its
    # normally distributed between-subject variability, theta_2 the
    # slope of the exposure-response line, and CP the individually
    # predicted piperaquine plasma concentration. A linear direct
    # response was retained over power, Emax and delayed-response
    # (turnover / link) alternatives (Results, 'Electrocardiographic
    # effects of piperaquine').
    # ==================================================================
    e0 <- fixed(0)
    label("Typical baseline double-delta-corrected QTcI prolongation (ms)")
    # Chotsiri 2017 Table 2: BASE (ms) = '0 Fixed'. Results: 'The population
    # baseline DeltaDeltaQTc prolongation was estimated close to zero and
    # therefore fixed to this value but allowed for between-subject
    # variability in the same parameter.' Not log-transformed: the baseline
    # of a double-delta-corrected change can take either sign.

    slope <- 0.0417
    label("Increase in double-delta-corrected QTcI prolongation per ng/mL of plasma piperaquine (ms per ng/mL)")
    # Chotsiri 2017 Table 2: SLOPE [ms (ng/mL)^-1] = 0.0417 (%RSE 12.5;
    # bootstrap 95% CI 0.0313-0.0511), i.e. 4.17 ms per 100 ng/mL, which is
    # the figure quoted in the Abstract and Results. Not log-transformed so
    # that a QTc-shortening drug could reuse this canonical unchanged.

    # Between-subject variability on the baseline. This is the ONE
    # variability term of the pharmacodynamic model: 'No major
    # between-subject variability was observed in other pharmacodynamic
    # parameters in the final model' (Results).
    #
    # Table 2's variability column is headed '%CV of BSV/BOV', but the
    # %CV transform in footnote a cannot apply to a parameter that is
    # additive in ms and fixed at zero -- a coefficient of variation
    # about zero is undefined. The 15.9 is therefore a standard
    # deviation in ms, and the paper's own simulations confirm it: the
    # reported median maximum prolongation after monthly mass drug
    # administration is 18.9 ms with a 95% CI of -6.44 to 49.0 ms, and
    # an individual maximum can only come out NEGATIVE if the baseline
    # eta has an SD of roughly 16 ms (18.9 - 1.96 * 15.9 = -12; the
    # concentration spread narrows this to about -6). A variance
    # reading (SD 3.99 ms) cannot produce a negative maximum at all.
    etae0 ~ 252.81
    # Chotsiri 2017 Table 2: BASE BSV = 15.9 ms (%RSE 33.4; bootstrap 95% CI 0.973-43.11); an SD in ms, so the variance is 15.9^2 = 252.81

    # ==================================================================
    # Residual error. Table 2 footnote: 'sigma_PD, residual additive
    # error variance of DeltaDeltaQTc prolongation'. 146 is therefore a
    # variance in ms^2, giving an SD of 12.1 ms -- the only
    # dimensionally coherent reading, since an additive residual SD of
    # 146 ms on an interval whose baseline is ~420 ms is impossible.
    # ==================================================================
    addSd <- sqrt(146)
    label("Additive residual standard deviation on the double-delta-corrected QTcI prolongation (ms)")
    # Chotsiri 2017 Table 2: sigma_PD (ms) = 146 (%RSE 25.5; bootstrap 95% CI 82.1-220), a variance -> SD = sqrt(146) = 12.1 ms
  })

  model({
    # ---- embedded piperaquine PK -------------------------------------
    iov_fdepot <- (OCC == 1) * etaiov_fdepot_1 + (OCC == 2) * etaiov_fdepot_2
    iov_mtt <- (OCC == 1) * etaiov_mtt_1 + (OCC == 2) * etaiov_mtt_2

    allom_cl <- (WT / 64)^e_wt_cl
    allom_v <- (WT / 64)^e_wt_vc

    mtt <- exp(lmtt + iov_mtt)

    cl <- exp(lcl + etalcl) * allom_cl
    vc <- exp(lvc + etalvc) * allom_v
    q <- exp(lq) * allom_cl
    vp <- exp(lvp) * allom_v
    q2 <- exp(lq2 + etalq2) * allom_cl
    vp2 <- exp(lvp2) * allom_v

    fdepot <- exp(lfdepot + etalfdepot + iov_fdepot)

    ktr <- 3 / mtt

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(central) <- ktr * transit2 - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    f(depot) <- fdepot

    # Plasma piperaquine concentration in ng/mL, the latent driver of
    # the electrocardiographic effect.
    Cc <- 1000 * central / vc

    # ---- linear concentration-DeltaDeltaQTcI relationship -------------
    # Equation 8. No threshold and no saturation: the Emax and power
    # forms were tested and the line was retained, and the paper's own
    # simulations are explicitly 'based on the assumption that a linear
    # concentration-effect relationship continued at piperaquine plasma
    # levels over 500 ng/mL' (Results). The max() is a numerical guard
    # against a solver undershoot only.
    e0_i <- e0 + etae0

    QTcI <- e0_i + slope * max(Cc, 0)

    QTcI ~ add(addSd)
  })
}
