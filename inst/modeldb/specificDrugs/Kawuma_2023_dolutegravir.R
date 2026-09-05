Kawuma_2023_dolutegravir <- function() {
  description <- "Two-compartment population PK model for dolutegravir with lagged first-order absorption in healthy volunteers, quantifying the drug-drug interaction with rifabutin (-33.1% on central volume) alongside the previously reported rifampicin interaction (+143% on clearance)"
  reference <- paste(
    "Kawuma AN, Wasmann RE, Dooley KE, Boffito M, Maartens G, Denti P.",
    "Drug-drug interaction between rifabutin and dolutegravir: A population",
    "pharmacokinetic model. Br J Clin Pharmacol. 2023;89(3):1216-1221.",
    "doi:10.1111/bcp.15604.",
    "Extends the dolutegravir-rifampicin model of Kawuma AN, Wasmann RE,",
    "Dooley KE, Boffito M, Maartens G, Denti P. Population pharmacokinetic",
    "model and alternative dosing regimens for dolutegravir coadministered",
    "with rifampicin. Antimicrob Agents Chemother. 2022;66(6):e00215-22.",
    "doi:10.1128/aac.00215-22, which is the source for the fixed allometric",
    "exponents, the log-normal BSV/BOV random-effect structure, the combined",
    "residual-error form and the study-specific assay LLOQ values; every",
    "parameter VALUE below is from the 2023 paper's own Table 1.",
    sep = " "
  )
  vignette <- "Kawuma_2023_dolutegravir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on all disposition parameters with a reference weight of 70 kg (Kawuma 2023 Table 1 footnote b). Exponents fixed at 0.75 for the clearances and 1 for the volumes per Kawuma 2022 Methods ('The allometric exponents for clearance and volume parameters were fixed at 0.75 and 1, respectively'). Total body weight was preferred over fat-free mass, which did not describe the data any better (Kawuma 2022 Results).",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) in the canonical register; the SOURCE paper's reference category is FEMALE",
      notes              = "The source reports 'Male sex on ka (%) = -38.1', i.e. a male indicator against a female reference, so the printed typical ka of 1.63 /h belongs to women. To keep both the canonical 1 = female orientation and the verbatim printed values, the effect is applied to the complement, ka * (1 + e_sexf_ka * (1 - SEXF)) with e_sexf_ka = -0.381 -- the same construction registered for Bajaj_2017_nivolumab.R and Wada_2023_sparsentan.R. Men absorb dolutegravir 38.1% more slowly; no sex difference in clearance or bioavailability was identified (Kawuma 2022 Results).",
      source_name        = "Male sex"
    ),
    CONMED_RIFAMPICIN = list(
      description        = "Concomitant rifampicin coadministration indicator, 1 = dolutegravir taken with rifampicin 600 mg once daily, 0 = dolutegravir alone",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no rifampicin)",
      notes              = "Rifampicin 600 mg once daily, dosed for at least 7 days before the pharmacokinetic sampling visit in both source studies, so the estimate describes steady-state enzyme/transporter induction rather than the onset. Rifampicin induces UGT1A1 and CYP3A4 (the enzymes that clear dolutegravir) as well as P-glycoprotein and BCRP.",
      source_name        = "Rifampicin co-administration"
    ),
    CONMED_RIFABUTIN = list(
      description        = "Concomitant rifabutin coadministration indicator, 1 = dolutegravir taken with rifabutin 300 mg once daily, 0 = dolutegravir alone",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no rifabutin)",
      notes              = "Rifabutin 300 mg once daily for 14 days (arm B of NCT01231542) before the pharmacokinetic sampling visit, so the estimate describes steady-state induction. This is the interaction the 2023 paper adds: the authors tested rifabutin on clearance, on bioavailability, and on both jointly, and selected the effect on central volume on goodness-of-fit, plausibility and BIC grounds (Kawuma 2023 Results). The mechanism is speculative -- the authors hypothesise induction of P-glycoprotein restricting tissue distribution (Kawuma 2023 Discussion).",
      source_name        = "Rifabutin co-administration"
    ),
    FED = list(
      description        = "Fed-versus-fasted dose-record indicator, 1 = dose taken after a meal, 0 = dose taken fasted",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Kawuma 2023 Table 1 reports two absorption lag times, labelled by study and annotated by prandial state in footnotes c and d: 0.205 h for NCT01231542 ('Dolutegravir taken under fasted conditions', an overnight fast) and 0.986 h for RADIO ('Under fed condition', a standard breakfast). The two studies are perfectly confounded with the two prandial states in the source data, but the authors themselves index the lag time by prandial state rather than by study when simulating (Kawuma 2022 Figure 3 caption: 'For the absorption lag time, we used the value estimated for the NCT01231542 study since this was done under fasted conditions') and attribute the difference mechanistically to food (Kawuma 2022 Discussion, citing Song 2012). Bioavailability did not differ between the two studies, which the authors ascribe to the low fat content of the RADIO breakfast; because the meal is described only as 'a standard breakfast' the generic FED indicator is used rather than FED_LOWFAT.",
      source_name        = "Study (RADIO = fed, NCT01231542 = fasted)"
    ),
    STUDY_RADIO = list(
      description        = "RADIO study indicator, 1 = record from the RADIO healthy-volunteer study, 0 = record from NCT01231542",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (NCT01231542)",
      notes              = "Selects the study-specific additive residual error, which tracks the two studies' bioanalytical assays rather than any physiology: RADIO used a UHPLC-UV method validated over 0.050-10 mg/L and NCT01231542 an LC-MS/MS method validated over 0.020-20 mg/L (Kawuma 2022 'Analytical assay'), and the additive error was constrained to be at least 20% of each study's LLOQ (Kawuma 2022 Methods). The printed additive errors of 0.036 and 0.0485 mg/L are the resulting TOTAL additive standard deviations -- their lower confidence bounds, 0.00435 and 0.0103 mg/L, sit on the 0.2 x LLOQ floors of 0.004 and 0.010 mg/L. Perfectly confounded with FED in the source data; the two indicators are kept separate because they carry different mechanisms and land on different parts of the model.",
      source_name        = "Study (RADIO vs NCT01231542)"
    ),
    OCC = list(
      description        = "Integer-valued occasion indicator for between-occasion-variability multiplexing",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Kawuma 2022 Methods defines an occasion as 'a dosing event leading to at least one observation', with each pharmacokinetic sampling visit contributing two occasions (the predose occasion and the postdose occasion). A single between-occasion variance per absorption parameter is reported and shared by every occasion, so occasion 1 carries the estimated variance and later occasions fix it to the same value (the registered NONMEM $OMEGA BLOCK(1) SAME idiom -- see Blackman_2026_methotrexate.R, Goggin_2004_emfilermin.R, Barnett_2018_coproporphyrin_I.R). Four occasions are encoded, which spans a two-period within-subject crossover with a predose and a postdose occasion in each period; because the variance is common to every occasion, extending the chain is a mechanical copy of the fix() lines. For single-occasion records pass OCC = 1.",
      source_name        = "OCC"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "dolutegravir", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "dolutegravir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dolutegravir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 41L,
    n_studies      = 2L,
    age_median     = "43 years",
    age_range      = "31-50 years (interquartile range)",
    weight_median  = "81.5 kg",
    weight_range   = "69.5-88.6 kg (interquartile range)",
    sex_female_pct = 32,
    disease_state  = "healthy HIV-negative volunteers",
    dose_range     = "dolutegravir 50 mg once daily, 50 mg twice daily and 100 mg once daily, alone or with rifampicin 600 mg once daily or rifabutin 300 mg once daily",
    regions        = "United Kingdom (RADIO, London) and United States (NCT01231542, Baltimore MD)",
    notes          = "Kawuma 2023 Results: 41 volunteers (68% male), 16 in RADIO and 25 in NCT01231542 (12 in arm A, 13 in arm B), contributing 907 dolutegravir plasma concentrations of which 90 were taken during rifabutin coadministration. No sample was below the lower limit of quantification. Fuller demographics are in Kawuma 2022 Table 1. The externally-validating INSPIRING patient cohort of Kawuma 2022 did not contribute to this model: the authors could not reliably co-model healthy-volunteer and patient data and note that the healthy-volunteer model slightly underpredicts trough concentrations in patients, making target-attainment predictions conservative."
  )

  ini({
    # =========================================================================
    # Structural parameters. All are apparent (oral) values for a typical 70 kg
    # individual, from Kawuma 2023 Table 1 'Typical value (95% CI)'; the
    # confidence intervals are non-parametric bootstrap (n = 500 replicates)
    # per Table 1 footnote a.
    # =========================================================================
    lcl <- log(1.03)
    label("Apparent clearance CL/F (L/h)")
    # Kawuma 2023 Table 1 row 'CL/F (L/h)' = 1.03 (0.945, 1.15)

    lvc <- log(13.3)
    label("Apparent central volume of distribution Vc/F (L)")
    # Kawuma 2023 Table 1 row 'Vc/F (L)' = 13.3 (11.7, 14.5)

    lq <- log(0.675)
    label("Apparent inter-compartmental clearance Q/F (L/h)")
    # Kawuma 2023 Table 1 row 'Q/F (L/h)' = 0.675 (0.439, 1.02)

    lvp <- log(3.52)
    label("Apparent peripheral volume of distribution Vp/F (L)")
    # Kawuma 2023 Table 1 row 'Vp/F (L)' = 3.52 (2.69, 4.24)

    lka <- log(1.63)
    label("Absorption rate constant ka in women (1/h)")
    # Kawuma 2023 Table 1 row 'Absorption rate constant (/h)' = 1.63 (1.16, 2.73).
    # The reference category of the sex effect below is female, so this is the
    # typical value for a woman.

    lfdepot <- fixed(log(1))
    label("Relative bioavailability (unitless)")
    # Kawuma 2023 Table 1 row 'Relative bioavailability' = 1 fixed

    ltlag <- log(0.205)
    label("Absorption lag time under fasted dosing (h)")
    # Kawuma 2023 Table 1 row 'Absorption lag time NCT01231542 study (h)' =
    # 0.205 (0.00205, 0.381); Table 1 footnote c: 'Dolutegravir taken under
    # fasted conditions'.

    # =========================================================================
    # Covariate effects.
    # =========================================================================
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on the clearance parameters (unitless)")
    # Kawuma 2023 Table 1 footnote b (allometric scaling with weight, median of
    # 70 kg, on CL, Q, Vc and Vp); exponent value from Kawuma 2022 Methods
    # 'Pharmacokinetic analysis': 'The allometric exponents for clearance and
    # volume parameters were fixed at 0.75 and 1, respectively'.

    e_wt_vc <- fixed(1)
    label("Allometric exponent on the volume parameters (unitless)")
    # Kawuma 2023 Table 1 footnote b; exponent value from Kawuma 2022 Methods
    # (fixed at 1 for volume parameters).

    e_fed_tlag <- log(0.986 / 0.205)
    label("Log-ratio of the fed to the fasted absorption lag time (unitless)")
    # Kawuma 2023 Table 1 rows 'Absorption lag time RADIO study (h)' = 0.986
    # (0.602, 1.43) with footnote d 'Under fed condition' and 'Absorption lag
    # time NCT01231542 study (h)' = 0.205 with footnote c 'fasted'. Encoded as
    # a log-ratio so that exp(ltlag) = 0.205 h fasted and
    # exp(ltlag + e_fed_tlag) = 0.986 h fed reproduce both printed values
    # exactly. log(0.986 / 0.205) = 1.5710.

    e_rifampicin_cl <- 1.43
    label("Fractional increase in CL/F with rifampicin coadministration (unitless)")
    # Kawuma 2023 Table 1 covariate row 'Rifampicin co-administration on CL
    # (%)' = +143 (+126, +155), i.e. clearance is multiplied by 2.43.

    e_rifabutin_vc <- -0.331
    label("Fractional change in Vc/F with rifabutin coadministration (unitless)")
    # Kawuma 2023 Table 1 covariate row 'Rifabutin co-administration on Vc (%)'
    # = -33.1 (-25.1, -42.3), also quoted in the Results as 'a 33.1% decrease
    # (95% CI 25.1-42.3) in dolutegravir's central volume of distribution
    # (dOFV = -56, df = 1, P < .001)'.

    e_sexf_ka <- -0.381
    label("Fractional change in ka in men relative to the female reference (unitless)")
    # Kawuma 2023 Table 1 covariate row 'Male sex on ka (%)' = -38.1 (-15.1,
    # -56.1). The source reference category is female, so this coefficient is
    # applied to (1 - SEXF) in model() to preserve both the canonical
    # 1 = female orientation of SEXF and the printed typical ka of 1.63 /h.

    # =========================================================================
    # Between-subject variability. Kawuma 2023 Table 1 section 'Parameter
    # variability (% CV)', footnote e: 'Calculated by CV% = sqrt(omega^2) * 100'
    # -- so the tabulated percentage is 100 x the log-scale standard deviation
    # omega, and the variance nlmixr2 wants is (CV% / 100)^2. Random effects
    # are log-normal (Kawuma 2022 Table 2 footnote c). BSV was tested on the
    # disposition parameters and retained on clearance only.
    # =========================================================================
    etalcl ~ 0.063001
    # Kawuma 2023 Table 1 row 'BSV clearance' = 25.1 (15.7, 34.4) %CV;
    # omega^2 = 0.251^2 = 0.063001

    # =========================================================================
    # Between-occasion variability on the three absorption parameters (Kawuma
    # 2022 Methods: 'BOV was incorporated on absorption parameters, while BSV
    # was tested on disposition parameters'). One variance per parameter is
    # reported and shared across occasions; nlmixr2 has no $OMEGA BLOCK(1) SAME
    # shortcut, so occasion 1 carries the estimated variance and occasions 2-4
    # fix it to the same value.
    # =========================================================================
    etaiov_fdepot_1 ~ 0.1849
    # Kawuma 2023 Table 1 row 'BOV bioavailability' = 43.0 (37.0, 49.6) %CV;
    # omega^2 = 0.430^2 = 0.1849 (estimated, occasion 1)
    etaiov_fdepot_2 ~ fix(0.1849)  # SAME-equivalent: equal to the occasion-1 variance
    etaiov_fdepot_3 ~ fix(0.1849)  # SAME-equivalent: equal to the occasion-1 variance
    etaiov_fdepot_4 ~ fix(0.1849)  # SAME-equivalent: equal to the occasion-1 variance

    etaiov_tlag_1 ~ 0.219961
    # Kawuma 2023 Table 1 row 'BOV lag time' = 46.9 (16.5, 79.8) %CV;
    # omega^2 = 0.469^2 = 0.219961 (estimated, occasion 1)
    etaiov_tlag_2 ~ fix(0.219961)  # SAME-equivalent: equal to the occasion-1 variance
    etaiov_tlag_3 ~ fix(0.219961)  # SAME-equivalent: equal to the occasion-1 variance
    etaiov_tlag_4 ~ fix(0.219961)  # SAME-equivalent: equal to the occasion-1 variance

    etaiov_ka_1 ~ 0.368449
    # Kawuma 2023 Table 1 row 'BOV absorption rate constant' = 60.7 (41.9,
    # 79.2) %CV; omega^2 = 0.607^2 = 0.368449 (estimated, occasion 1)
    etaiov_ka_2 ~ fix(0.368449)  # SAME-equivalent: equal to the occasion-1 variance
    etaiov_ka_3 ~ fix(0.368449)  # SAME-equivalent: equal to the occasion-1 variance
    etaiov_ka_4 ~ fix(0.368449)  # SAME-equivalent: equal to the occasion-1 variance

    # =========================================================================
    # Residual unexplained variability (Kawuma 2023 Table 1 section 'Residual
    # unexplained variability'). Combined proportional plus additive on the
    # linear mg/L scale; the additive component is study-specific because it
    # was constrained to be at least 20% of each study's assay LLOQ (Kawuma
    # 2022 Methods and 'Analytical assay').
    # =========================================================================
    propSd <- 0.0885
    label("Proportional residual error (fraction)")
    # Kawuma 2023 Table 1 row 'proportional error (%)' = 8.85 (7.32, 9.79)

    addSd <- 0.036
    label("Additive residual error for the NCT01231542 assay (mg/L)")
    # Kawuma 2023 Table 1 row 'Additive error NCT01231542 study (mg/L)' =
    # 0.036 (0.00435, 0.0801)

    addSd_radio <- 0.0485
    label("Additive residual error for the RADIO assay (mg/L)")
    # Kawuma 2023 Table 1 row 'Additive error RADIO study (mg/L)' = 0.0485
    # (0.0103, 0.0861)
  })

  model({
    # 1. Decompose the integer occasion column into binary indicators and
    #    multiplex the between-occasion etas onto the three absorption
    #    parameters. For single-occasion records pass OCC = 1.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)

    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4
    iov_tlag <- oc1 * etaiov_tlag_1 + oc2 * etaiov_tlag_2 +
      oc3 * etaiov_tlag_3 + oc4 * etaiov_tlag_4
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 +
      oc3 * etaiov_ka_3 + oc4 * etaiov_ka_4

    # 2. Individual parameters. Allometric scaling on all four disposition
    #    parameters with a 70 kg reference, exponents fixed at 0.75 (clearances)
    #    and 1 (volumes). Rifampicin raises clearance by 143%; rifabutin lowers
    #    central volume by 33.1%; men absorb 38.1% more slowly than the female
    #    reference; a fed dose has a ~4.8-fold longer absorption lag.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl *
      (1 + e_rifampicin_cl * CONMED_RIFAMPICIN)
    vc <- exp(lvc) * (WT / 70)^e_wt_vc *
      (1 + e_rifabutin_vc * CONMED_RIFABUTIN)
    q <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    ka <- exp(lka + iov_ka) * (1 + e_sexf_ka * (1 - SEXF))
    tlag <- exp(ltlag + e_fed_tlag * FED + iov_tlag)
    fdepot <- exp(lfdepot + iov_fdepot)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. Two-compartment disposition with first-order absorption from an oral
    #    depot and first-order elimination from the central compartment.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Relative bioavailability and absorption lag time on the depot.
    f(depot) <- fdepot
    alag(depot) <- tlag

    # 6. Observation and combined residual error. Dose in mg and vc in L give
    #    Cc in mg/L, the unit the source tabulates concentrations and the
    #    0.064 mg/L PA-IC90 / 0.3 mg/L EC90 targets in.
    Cc <- central / vc

    # The additive component is selected by study; the proportional component
    # is common to both.
    addSd_study <- addSd * (1 - STUDY_RADIO) + addSd_radio * STUDY_RADIO
    Cc ~ add(addSd_study) + prop(propSd)
  })
}
