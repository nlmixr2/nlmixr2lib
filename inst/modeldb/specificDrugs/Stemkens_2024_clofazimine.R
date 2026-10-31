Stemkens_2024_clofazimine <- function() {
  description <- paste(
    "Two-compartment population PK model for oral clofazimine in adults with",
    "pulmonary or extrapulmonary nontuberculous mycobacterial (NTM) disease,",
    "fitted to the pooled C-LOAD loading-dose study (300 mg once daily for 4",
    "weeks, then 100 mg once daily) and the PERC trial (100 mg once daily,",
    "no loading dose). Absorption runs through a depot and one transit",
    "compartment that share a single rate constant (0.683 1/h). Apparent",
    "clearance (4.17 L/h), central volume (460 L), intercompartmental",
    "clearance (24.6 L/h) and peripheral volume (10100 L) are allometrically",
    "scaled on total body weight (70 kg reference; fixed exponents 0.75 for",
    "clearances and 1 for volumes). Relative bioavailability of each 300 mg",
    "dose is 0.736 times that of a 100 mg dose. Random effects are",
    "between-subject variability on clearance and on the absorption rate",
    "constant and inter-occasion variability on bioavailability of the doses",
    "given on PK sampling days; residual error is combined proportional plus",
    "additive."
  )
  reference <- paste(
    "Stemkens R, Lemson A, Koele SE, Svensson EM, te Brake LHM, van Crevel R,",
    "Boeree MJ, Hoefsloot W, van Ingen J, Aarnoutse RE (2024).",
    "A loading dose of clofazimine to rapidly achieve steady-state-like",
    "concentrations in patients with nontuberculous mycobacterial disease.",
    "J Antimicrob Chemother 79(12):3100-3108. doi:10.1093/jac/dkae309"
  )
  vignette <- "Stemkens_2024_clofazimine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "clofazimine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "clofazimine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "clofazimine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "clofazimine", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of CL and Q (exponent 0.75) and of V and Vp",
        "(exponent 1), reference 70 kg (Table 4 footnote b; Text S3 $PK",
        "'AlloCL = (WT/70)**0.75', 'AlloV = (WT/70)**1'). Total body weight",
        "was preferred over fat-free mass and fat mass as the size descriptor",
        "(Text S2)."
      ),
      source_name = "WT"
    ),
    DOSE_HIGH = list(
      description = "High-dose (300 mg) dose-record indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (100 mg dose)",
      notes = paste(
        "Per dose record: 1 for each 300 mg clofazimine dose (the C-LOAD",
        "loading phase), 0 for each 100 mg dose. Multiplies the",
        "bioavailability of that dose by 0.736 (Table 4 'Effect of",
        "high-dose on bioavailability'; Text S3 $PK 'IF(highdose.EQ.1)",
        "highdoseEFF = THETA(6)'). The source control stream references",
        "'highdose' without listing it in $INPUT; the Results state that the",
        "contrast is 300 mg doses against 100 mg doses. Only 100 mg and 300 mg",
        "were studied, so the indicator is undefined for other dose levels."
      ),
      source_name = "highdose"
    ),
    OCC = list(
      description = "PK-sampling occasion index for the inter-occasion variability on bioavailability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Text S3 data dictionary: 'Occasion (0 for no PK sample, 1/2/3 for",
        "each PK day)'. Each occasion carries its own eta on bioavailability",
        "(ETA(3)-ETA(5), one shared variance via $OMEGA BLOCK(1) SAME); doses",
        "on days without PK sampling (OCC = 0) carry no inter-occasion term.",
        "In C-LOAD the occasions were the last 300 mg loading dose (day 28),",
        "the first 100 mg maintenance dose (day 29) and month 4; PERC had",
        "occasions at month 1 and month 4. Set OCC = 0 on every record to",
        "simulate without inter-occasion variability."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_RIFAMPICIN = list(
      description = "Concomitant rifampicin",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no rifampicin)",
      notes = "Tested in the stepwise covariate search (Text S2) and not retained (Results, PK model development). 7 of 12 C-LOAD participants used rifampicin (Table 1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in the stepwise covariate search (Text S2) and not retained."
    ),
    SEXF = list(
      description = "Female sex",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Tested in the stepwise covariate search (Text S2) and not retained."
    ),
    DISEASE_SITE = list(
      description = "Pulmonary versus extrapulmonary NTM disease",
      units = "(binary)",
      type = "binary",
      reference_category = NULL,
      notes = "Tested in the stepwise covariate search (Text S2) and not retained. Name used here for documentation only; not a registered canonical."
    ),
    STUDY = list(
      description = "Study identifier (C-LOAD = 0, PERC = 1)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (C-LOAD)",
      notes = "Present in the Text S3 $INPUT but unused in the final model; the Results report no significant difference in absorption rate, volume, clearance or bioavailability between the two studies."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 2L,
    age_range = "39-82 years (C-LOAD, median 69)",
    weight_range = "51-105 kg (C-LOAD, median 62)",
    ffm_range = "38-70 kg (C-LOAD, median 50)",
    sex_female_pct = 50,
    race_ethnicity = "C-LOAD: 10 of 12 Caucasian, 2 of 12 Middle Eastern",
    disease_state = paste(
      "Pulmonary (9 of 12) or extrapulmonary (3 of 12) nontuberculous",
      "mycobacterial disease in C-LOAD, mostly M. avium complex; PERC enrolled",
      "patients with M. avium complex pulmonary disease."
    ),
    dose_range = paste(
      "C-LOAD: clofazimine 300 mg once daily for 4 weeks, then 100 mg once",
      "daily to month 4. PERC: 100 mg once daily without loading dose. Taken",
      "with food, combined with other antimycobacterial drugs (azithromycin,",
      "ethambutol, rifampicin, amikacin and others)."
    ),
    regions = "Netherlands (Radboud University Medical Center, Nijmegen)",
    notes = paste(
      "Pooled data: 112 concentrations from 12 C-LOAD participants",
      "(ClinicalTrials.gov NCT05294146) and 67 concentrations from 19 PERC",
      "participants (Results, PK model development). Demographics are",
      "reported for C-LOAD only (Table 1); the PERC cohort is described in",
      "Zweijpfenning 2024 (Chest 165:1082-92)."
    )
  )

  ini({
    # Structural parameters: Table 4 typical values (95% CI from SIR), which
    # are identical to the Text S3 $THETA values. All are apparent (oral)
    # parameters for a 70 kg subject (Table 4 footnote b).
    lcl <- log(4.17)
    label("Apparent clearance CL/F at 70 kg (L/h)") # Table 4 'Clearance (L/h)' 4.17 (2.81-5.79); Text S3 THETA(1)
    lvc <- log(460)
    label("Apparent central volume V/F at 70 kg (L)") # Table 4 'Central volume of distribution (L)' 460 (370-608); Text S3 THETA(2)
    lka <- log(0.683)
    label("First-order rate constant from depot to transit and from transit to central (1/h)") # Table 4 'Absorption-rate constant (1/h)' 0.683 (0.502-0.940); Text S3 THETA(3) TVKTR used for K12 and K23
    lq <- log(24.6)
    label("Apparent intercompartmental clearance Q/F at 70 kg (L/h)") # Table 4 'Intercompartmental clearance (L/h)' 24.6 (18.1-32.9); Text S3 THETA(4)
    lvp <- log(10100)
    label("Apparent peripheral volume Vp/F at 70 kg (L)") # Table 4 'Peripheral volume of distribution (L)' 10100 (7670-13000); Text S3 THETA(5)
    lfdepot <- fixed(log(1))
    label("Bioavailability of a 100 mg dose (fraction)") # Text S3 $PK 'F1 = 1* IOVF *highdoseEFF'; Table 4 footnote b, apparent parameters

    # Covariate effects.
    e_dose_high_fdepot <- 0.736
    label("Bioavailability of a 300 mg dose relative to a 100 mg dose (fraction)") # Table 4 'Effect of high-dose on bioavailability' 73.6 (59.6-91.8); Text S3 THETA(6) 0.736
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL and Q (unitless)") # Table 4 footnote b; Text S3 'AlloCL = (WT/70)**0.75'
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on V and Vp (unitless)") # Table 4 footnote b; Text S3 'AlloV = (WT/70)**1'

    # Between-subject variability. Text S3 $OMEGA carries the log-scale
    # variances; Table 4 footnote c converts them to CV = sqrt(exp(OMEGA) - 1),
    # and the two agree: sqrt(exp(0.55) - 1) = 0.856, sqrt(exp(0.501) - 1) = 0.806.
    etalcl ~ 0.55 # Text S3 $OMEGA '0.55 ; 1 IIV CL'; Table 4 IIV CL 85.6 CV pct (55.9-143)
    etalka ~ 0.501 # Text S3 $OMEGA '0.501 ; 2 IIV KTR'; Table 4 IIV absorption-rate constant 80.6 CV pct (47.8-139)

    # Inter-occasion variability on bioavailability: three occasions sharing
    # one variance ($OMEGA BLOCK(1) 0.488 then two BLOCK(1) SAME). nlmixr2 has
    # no SAME shortcut, so occasions 2 and 3 are fixed to the occasion-1 value.
    etaiov_fdepot_1 ~ 0.488 # Text S3 $OMEGA BLOCK(1) '0.488 ; 3 IOV F'; Table 4 IOV bioavailability 79.3 CV pct (59.7-122)
    etaiov_fdepot_2 ~ fix(0.488) # Text S3 $OMEGA BLOCK(1) SAME (ETA(4))
    etaiov_fdepot_3 ~ fix(0.488) # Text S3 $OMEGA BLOCK(1) SAME (ETA(5))

    # Residual error: Text S3 $ERROR 'Y = IPRED*(1+PROPERR)+ADDERR'.
    propSd <- 0.08
    label("Proportional residual error (fraction)") # Table 4 'Proportional error' 8.00 (5.95-10.5); Text S3 $SIGMA 0.00639 = 0.0799^2
    addSd <- 0.0417
    label("Additive residual error SD (mg/L)") # Text S3 $SIGMA '0.00174 ; Additive error' is a variance, sqrt(0.00174) = 0.0417; Table 4 prints 0.00174 (see vignette Errata)
  })

  model({
    # Inter-occasion bioavailability term: the control stream sets
    # IOVF = EXP(ETA(k)) when OCC == k and IOVF = 1 when OCC == 0.
    iov_fdepot <- (OCC == 1) * etaiov_fdepot_1 +
      (OCC == 2) * etaiov_fdepot_2 +
      (OCC == 3) * etaiov_fdepot_3

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    ka <- exp(lka + etalka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ADVAN5: depot -> transit -> central at the same rate constant
    # (K12 = K23 = KTR); elimination and distribution from central.
    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(central) <- ka * transit1 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # F1 = 1 * IOVF * highdoseEFF, evaluated per dose record.
    f(depot) <- exp(lfdepot + iov_fdepot) * e_dose_high_fdepot^DOSE_HIGH

    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
