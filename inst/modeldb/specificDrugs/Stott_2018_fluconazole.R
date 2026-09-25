Stott_2018_fluconazole <- function() {
  description <- paste(
    "Four-compartment population PK model for oral fluconazole in adults with",
    "HIV-associated cryptococcal meningitis (Stott 2018): a gut absorption",
    "compartment feeding a central compartment that exchanges with a CSF/CNS",
    "compartment and with a peripheral compartment through asymmetric",
    "first-order transfer rate constants, with first-order elimination from",
    "central. Non-parametric (Pmetrics NPAG) population estimates; no covariate",
    "was retained in the final model.",
    sep = " "
  )
  reference <- paste(
    "Stott KE, Beardsley J, Kolamunnage-Dona R, Castelazo AS, Kibengo FM,",
    "Mai NTH, Tung NLN, Cuc NTK, Day J, Hope W. Population pharmacokinetics",
    "and cerebrospinal fluid penetration of fluconazole in adults with",
    "cryptococcal meningitis. Antimicrob Agents Chemother. 2018;62(9):e00885-18.",
    "doi:10.1128/AAC.00885-18.",
    sep = " "
  )
  vignette <- "Stott_2018_fluconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Stott 2018 "Population pharmacokinetic covariate screening": patient
  # weight, BMI, sex, ethnicity and baseline eGFR were screened by
  # bidirectional stepwise multivariate linear regression against the Bayesian
  # posterior estimates of volume and clearance. Only weight reached the
  # screening threshold against volume, and then only marginally (slope 0.22,
  # 95% CI -0.06 to 0.51, P = 0.05); a variant model scaling Vc to weight
  # (Table 2, "Model 2") gave comparable log likelihood, AIC, bias and
  # imprecision, so the covariate-free base model ("Model 1") was retained as
  # the final model. These covariates are therefore documented but never
  # referenced in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened; not retained. The only covariate to reach the screening",
        "threshold against volume of distribution (slope 0.22, 95% CI -0.06 to",
        "0.51, P = 0.05), but the weight-scaled variant (Table 2, Model 2) did",
        "not improve the fit and was dropped. Median 48 kg (range 32-68 kg).",
        "The Discussion attributes the null result to the narrow weight range",
        "in this severely wasted cohort (median BMI 18 kg/m^2). Weight still",
        "sets the companion amphotericin B dose but does not scale any",
        "fluconazole model parameter.",
        sep = " "
      ),
      source_name = "weight"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened; not retained. Median 18 kg/m^2 (range 12-25 kg/m^2), n = 35.",
      source_name = "BMI"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = paste(
        "Screened; not retained. Stott 2018 Results: clearance 0.79 L/h in",
        "males versus 0.66 L/h in females (P = 0.09) and volume 18.07 L in",
        "both (P = 0.97). The source reports sex as counts of male and female",
        "patients, so SEXF is the natural canonical encoding. See the",
        "population notes for the Table 1 sex-count discrepancy.",
        sep = " "
      ),
      source_name = "sex"
    ),
    REGION_VIETNAM = list(
      description = "Vietnam enrollment-site indicator (0 = Uganda)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = Uganda (Masaka General Hospital)",
      notes = paste(
        "Screened; not retained. The source screens this two-level",
        "enrollment-site variable under the name 'ethnicity' (Vietnamese",
        "versus Ugandan patients); the two levels are exactly the two study",
        "sites, so the canonical enrollment-country indicator is used here.",
        "Stott 2018 Results: clearance 0.74 L/h (95% CI 0.64-0.83) in",
        "Vietnamese versus 0.71 L/h (0.59-0.82) in Ugandan patients",
        "(P = 0.51); volume 16.88 L (14.33-19.44) versus 19.44 L",
        "(16.88-22.0) (P = 0.16).",
        sep = " "
      ),
      source_name = "ethnicity"
    ),
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate (Cockcroft-Gault)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = paste(
        "Screened; not retained. Reported by the source as eGFR in",
        "mL/min/1.73 m^2 by the Cockcroft-Gault equation; Cockcroft-Gault",
        "natively yields mL/min, so the BSA normalisation stated in Table 1 is",
        "the source's own description and is recorded verbatim rather than",
        "reinterpreted. Median 84.8 mL/min/1.73 m^2 (range 35.4-146.7), n = 33.",
        "The Discussion attributes the absence of a creatinine-clearance",
        "effect on fluconazole clearance to the narrow range of renal function",
        "in this cohort.",
        sep = " "
      ),
      source_name = "eGFR"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "fluconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "fluconazole", units = "mg", specimen = "plasma", verified = TRUE),
    csf = list(analyte = "fluconazole", units = "mg", specimen = "CSF", verified = TRUE),
    peripheral1 = list(analyte = "fluconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 43,
    n_studies = 2,
    age_median = "33 years (range 20-73 years)",
    weight_median = "48 kg (range 32-68 kg)",
    sex_female_pct = 51,
    disease_state = "HIV-associated cryptococcal meningitis; severely wasted (median BMI 18 kg/m^2, range 12-25)",
    dose_range = paste(
      "Fluconazole given orally, or by nasogastric tube where conscious level",
      "precluded swallowing. The majority of patients received 800 mg q24h;",
      "two patients received one-off doses of 400 mg, two received one-off",
      "doses of 600 mg, and one patient's 800 mg q24h regimen was escalated to",
      "1200 mg q24h for 6 days from day 8. All patients also received",
      "amphotericin B deoxycholate 1 mg/kg q24h infused over 5-6 h.",
      sep = " "
    ),
    renal_function = "Baseline Cockcroft-Gault eGFR median 84.8 mL/min/1.73 m^2 (range 35.4-146.7, n = 33); baseline creatinine median 70 umol/L (range 37-167)",
    regions = "Vietnam (Hospital for Tropical Diseases, Ho Chi Minh City) and Uganda (Masaka General Hospital)",
    notes = paste(
      "43 patients (23 Vietnam, 20 Uganda) recruited over 11 months between",
      "January and November 2016: 3 from a multicentre randomised controlled",
      "trial of adjuvant dexamethasone (ISRCTN 59144167) and, after that trial",
      "stopped early, 40 from a prospective descriptive study at the same",
      "sites. 508 plasma observations (312 Vietnam, 196 Uganda) and 167 CSF",
      "observations (52 Vietnam, 115 Uganda), means of 11.8 plasma and 3.9 CSF",
      "samples per patient; one further Ugandan CSF observation was excluded as",
      "unverifiable. Fluconazole was assayed by LC-MS/MS over 1-120 mg/L with a",
      "1 mg/L limit of quantitation. Fitted with the non-parametric adaptive",
      "grid algorithm of Pmetrics 1.5.0 under R 3.1.1.",
      "",
      "SEX COUNTS. Table 1's 'Combined' column reports 23 males and 20 females,",
      "which contradicts both its own per-site columns (13 + 8 = 21 males;",
      "10 + 12 = 22 females) and the Results text, 'Twenty-two patients (52%)",
      "were female'. The per-site sums and the Results text agree, so 22 female",
      "of 43 (51%) is recorded here and the Combined column is treated as a",
      "transcription error.",
      "",
      "ESTIMATED INITIAL CONDITIONS ARE DELIBERATELY NOT ENCODED. Some patients",
      "had taken fluconazole at an undocumented time before enrolment, leaving",
      "detectable drug in the first PK sample. To absorb this, the source",
      "estimated a non-zero initial condition for each of the four compartments,",
      "multiplied by a binary switch set to 1 when fluconazole was detected in",
      "the first PK sample and 0 otherwise (Table 3: ICgut mean 34.67 mg, SD",
      "22.74; ICcentral 35.86 mg, SD 19.67; ICCNS 31.06 mg, SD 23.47;",
      "ICperipheral 34.29 mg, SD 13.21; all four medians sit at 49.96-49.99 mg,",
      "essentially the midpoint of the search grid, so the distributions are",
      "close to unidentified). These are a fitting device for this specific",
      "43-patient data set, not a generalisable patient characteristic: the",
      "source itself states that 'for the simulations, the initial conditions of",
      "all compartments were defaulted to zero', so every simulated result the",
      "paper reports (the AUC distributions of Fig. 3 and the PTA analysis of",
      "Fig. 4) uses zero initial conditions. This forward-simulation model",
      "therefore starts from zero, matching the source's own simulation setup.",
      "Users wishing to reproduce the baseline-positive patients should pre-load",
      "the compartments with the Table 3 amounts. This follows the established",
      "handling of the same artefact in modellib('Debord_2001_cyclosporin')",
      "(residual concentration C0) and modellib('Mosha_2014_lumefantrine')",
      "(residual previous-treatment dose F0).",
      "",
      "Companion model from the same group, same disease and same four-",
      "compartment structure: modellib('Stott_2023_flucytosine').",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- Stott 2018 Table 3, "Mean" column. The source
    # fits a non-parametric (NPAG) support-point distribution and reports its
    # Mean, Median and SD per parameter; the Results state that "the mean
    # parameter estimates better fitted the data than medians and were used to
    # calculate Bayesian estimates of drug exposure for each individual
    # patient", so the mean set is the one carried forward here. All volumes
    # and the clearance are APPARENT (the source's SCL/F, Vc/F, Vcns/F
    # notation) because fluconazole was given orally or by nasogastric tube and
    # no bioavailability term was estimated.
    lka            <- log(8.78);  label("Absorption rate constant Ka, gut to central (1/h)")                # Stott 2018 Table 3: Ka mean 8.78 /h (median 1.73, SD 11.98)
    lcl            <- log(0.72);  label("Apparent clearance SCL/F from central (L/h)")                      # Stott 2018 Table 3: SCL/F mean 0.72 L/h (median 0.65, SD 0.24)
    lvc            <- log(18.07); label("Apparent central volume Vc/F (L)")                                 # Stott 2018 Table 3: Vc/F mean 18.07 L (median 17.41, SD 6.31)
    lk_central_csf <- log(35.43); label("Transfer rate constant central to CSF, paper's Kcs (1/h)")         # Stott 2018 Table 3: Kcs mean 35.43 /h (median 42.55, SD 13.74)
    lk_csf_central <- log(28.63); label("Transfer rate constant CSF to central, paper's Ksc (1/h)")         # Stott 2018 Table 3: Ksc mean 28.63 /h (median 29.04, SD 10.03)
    lk12           <- log(12.20); label("Transfer rate constant central to peripheral1, paper's Kcp (1/h)") # Stott 2018 Table 3: Kcp mean 12.20 /h (median 8.36, SD 11.17)
    lk21           <- log(18.10); label("Transfer rate constant peripheral1 to central, paper's Kpc (1/h)") # Stott 2018 Table 3: Kpc mean 18.10 /h (median 18.34, SD 8.25)
    lvcsf          <- log(32.07); label("Apparent CSF/CNS compartment volume Vcns/F (L)")                   # Stott 2018 Table 3: Vcns/F mean 32.07 L (median 30.49, SD 17.60)

    # IIV. The source is a non-parametric fit, so it publishes a support-point
    # distribution rather than an OMEGA matrix; Table 3 reports only a per-
    # parameter Mean and SD. Each variance below is the log-normal moment match
    # to that pair, omega^2 = log(1 + (SD/Mean)^2), which is the same encoding
    # the companion Stott_2023_flucytosine extraction uses (there the source
    # printed the CV% column directly; here CV is formed as SD/Mean from the
    # same table). The matrix is DIAGONAL because the source publishes no
    # covariances between support-point parameters -- the true NPAG joint
    # density is discrete and correlated, and a diagonal log-normal cannot
    # reproduce it exactly. See the vignette for what this approximation does
    # and does not preserve.
    etalka            ~ 1.05144   # Stott 2018 Table 3: Ka     SD 11.98 / mean  8.78 -> CV 1.36446 -> log(1 + CV^2)
    etalcl            ~ 0.10536   # Stott 2018 Table 3: SCL/F  SD  0.24 / mean  0.72 -> CV 0.33333 -> log(1 + CV^2)
    etalvc            ~ 0.11506   # Stott 2018 Table 3: Vc/F   SD  6.31 / mean 18.07 -> CV 0.34920 -> log(1 + CV^2)
    etalk_central_csf ~ 0.14010   # Stott 2018 Table 3: Kcs    SD 13.74 / mean 35.43 -> CV 0.38781 -> log(1 + CV^2)
    etalk_csf_central ~ 0.11577   # Stott 2018 Table 3: Ksc    SD 10.03 / mean 28.63 -> CV 0.35033 -> log(1 + CV^2)
    etalk12           ~ 0.60883   # Stott 2018 Table 3: Kcp    SD 11.17 / mean 12.20 -> CV 0.91557 -> log(1 + CV^2)
    etalk21           ~ 0.18876   # Stott 2018 Table 3: Kpc    SD  8.25 / mean 18.10 -> CV 0.45580 -> log(1 + CV^2)
    etalvcsf          ~ 0.26327   # Stott 2018 Table 3: Vcns/F SD 17.60 / mean 32.07 -> CV 0.54880 -> log(1 + CV^2)

    # Residual error. Stott 2018 Methods, "Population pharmacokinetic
    # modeling": "Model error was attributed separately to process noise
    # (including errors in sampling times or dosing) and assay variance.
    # Process noise was modeled using lambda, an additive error term. The data
    # were weighted by the inverse of the estimated assay variance." The FORM
    # is therefore stated -- an additive process-noise term plus a
    # concentration-dependent assay-variance weighting -- but no numeric value
    # for lambda or for the assay error polynomial is reported anywhere in the
    # paper, and there is no supplement. Both an additive and a proportional
    # term are carried for each output so the structure is available, and both
    # are fixed to zero rather than invented. The only precision figures the
    # paper gives are LC-MS/MS quality-control CVs (plasma intraday < 3.4%,
    # interday < 6.7%; CSF intraday < 5.2%, interday < 5.3%, over 1-90 mg/L),
    # which characterise the assay and not the model's residual variability.
    propSd      <- fixed(0); label("Proportional residual SD on plasma Cc (fraction; 0 -- not reported in the source)")  # Stott 2018 Methods: additive lambda plus assay-variance weighting declared, no estimate published
    addSd       <- fixed(0); label("Additive residual SD on plasma Cc (mg/L; 0 -- not reported in the source)")          # Stott 2018 Methods: additive lambda plus assay-variance weighting declared, no estimate published
    propSd_Ccsf <- fixed(0); label("Proportional residual SD on CSF Ccsf (fraction; 0 -- not reported in the source)")   # Stott 2018 Methods: additive lambda plus assay-variance weighting declared, no estimate published
    addSd_Ccsf  <- fixed(0); label("Additive residual SD on CSF Ccsf (mg/L; 0 -- not reported in the source)")           # Stott 2018 Methods: additive lambda plus assay-variance weighting declared, no estimate published
  })

  model({
    # Individual parameters. No covariate enters any of them: the final model
    # is the covariate-free base model (Stott 2018 Results, "Covariate
    # investigation", and Table 2 "Model 1").
    ka            <- exp(lka + etalka)
    cl            <- exp(lcl + etalcl)
    vc            <- exp(lvc + etalvc)
    vcsf          <- exp(lvcsf + etalvcsf)
    k_central_csf <- exp(lk_central_csf + etalk_central_csf)
    k_csf_central <- exp(lk_csf_central + etalk_csf_central)
    k12           <- exp(lk12 + etalk12)
    k21           <- exp(lk21 + etalk21)
    kel           <- cl / vc

    # Stott 2018 Equations 1-4. The source numbers its states 1 = gut,
    # 2 = central (c), 3 = CSF (s), 4 = peripheral (p), so its Kcs / Ksc are
    # the central <-> CSF pair and its Kcp / Kpc are the central <-> peripheral
    # pair. Both exchanges are parameterised by an asymmetric pair of
    # first-order rate constants rather than by a single inter-compartmental
    # clearance, and the asymmetry is real: Kcs * Vc = 640 L/h does not equal
    # Ksc * Vcns = 918 L/h, and Kcp * Vc = 220 L/h does not equal Kpc * Vp
    # (Vp is not estimated at all), so neither leg can be re-expressed with a
    # symmetric q.
    d/dt(depot)       <- -ka * depot                                                # Equation 1: dX(1)/dt = -Ka * X(1)
    d/dt(central)     <- ka * depot - (k12 + k_central_csf + kel) * central +
      k_csf_central * csf + k21 * peripheral1                                       # Equation 2: dX(2)/dt = Ka*X(1) - (Kcp + Kcs + SCL/V)*X(2) + Ksc*X(3) + Kpc*X(4)
    d/dt(csf)         <- k_central_csf * central - k_csf_central * csf              # Equation 3: dX(3)/dt = Kcs*X(2) - Ksc*X(3)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1                          # Equation 4: dX(4)/dt = Kcp*X(2) - Kpc*X(4)

    Cc   <- central / vc                                                            # Equation 5: Y(1) = X(2)/V
    Ccsf <- csf / vcsf                                                              # Equation 6: Y(2) = X(3)/Vcns

    Cc ~ add(addSd) + prop(propSd)
    Ccsf ~ add(addSd_Ccsf) + prop(propSd_Ccsf)
  })
}
