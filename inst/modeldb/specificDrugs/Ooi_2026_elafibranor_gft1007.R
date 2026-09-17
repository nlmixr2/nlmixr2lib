Ooi_2026_elafibranor_gft1007 <- function() {
  description <- paste0(
    "Two-compartment population pharmacokinetic model for GFT1007, the ",
    "active metabolite of elafibranor, after oral elafibranor dosing in the ",
    "pooled analysis of 17 clinical trials supporting the primary biliary ",
    "cholangitis (PBC) indication (894 subjects, 10,592 observations; doses ",
    "5-360 mg). The authors fitted the metabolite separately from the parent ",
    "(modellib('Ooi_2026_elafibranor')) because a simultaneous joint fit was ",
    "computationally intractable, so this model has the same structural form ",
    "as the parent: sequential zero-order input into a depot over a duration ",
    "D1 followed by first-order transfer at ka = 1/MAT after an absorption ",
    "lag, then two-compartment disposition with first-order elimination. The ",
    "depot therefore represents formation of GFT1007 rather than absorption ",
    "of the parent, and its parameters are apparent ones conditional on the ",
    "elafibranor dose. Relative bioavailability, D1, MAT and the lag time are ",
    "formulation-specific over six clinical formulations. Apparent clearance ",
    "and apparent inter-compartmental clearance scale allometrically with ",
    "baseline body weight at a fixed exponent of 0.75 and the two apparent ",
    "volumes at 1.0, referenced to 75 kg. Additional covariates are baseline ",
    "alanine aminotransferase on CL/F and Vp/F (power), baseline BMI on CL/F ",
    "and Vp/F, baseline creatinine clearance on CL/F, age on Vc/F, female ",
    "sex and PBC population on Q/F, and fed state on D1 and relative ",
    "bioavailability. Unlike the parent, GFT1007 carries no ",
    "inter-occasion variability, no dose effect on bioavailability and no ",
    "bioanalytical-method effect. The residual error is exponential ",
    "(additive on the log scale) with six magnitudes selected by study phase ",
    "and whether the sample was drawn within 2.5 h of the last dose. ",
    "GFT1007 exposure is about 5 times the parent AUC and 8 times the parent ",
    "Cmax at 80 mg/day; the exposure-response model that consumes both is ",
    "modellib('Ooi_2026_elafibranor_alptb')."
  )

  reference <- paste(
    "Ooi QX, Brendel K, van Beek S, Aguiar Zdovc J, Bardol M, Dehez M.",
    "Population Pharmacokinetics and Pharmacokinetics-Pharmacodynamics",
    "Analyses of Elafibranor to Support Dose Selection in Primary Biliary",
    "Cholangitis. CPT Pharmacometrics Syst Pharmacol. 2026;15(0):e70247.",
    "doi:10.1002/psp4.70247.",
    sep = " "
  )

  vignette <- "Ooi_2026_elafibranor"

  # As for the parent model, the NONMEM data set carried amounts in nmol and
  # volumes in mL; the tabulated L / L/h values are kept here with doses in
  # umol, giving numerically identical concentrations in umol/L. Both
  # analytes share one PK data set, so the dosed amount is the MOLAR
  # ELAFIBRANOR dose (dose_mg / 384.49 * 1000 umol); formation is 1:1 molar.
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Included a priori with allometric scaling referenced to 75 kg ",
        "(Supplementary Datafile S2 $PK, LOG(WTKGBL/75)); exponents fixed at ",
        "0.75 on CL/F and Q/F and 1.00 on Vc/F and Vp/F (Table S3)."
      ),
      source_name = "WTKGBL"
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Power model normalised to 27 U/L on CL/F and Vp/F (Supplementary ",
        "Datafile S2 $PK, CLALTBL = (ALTBL/27)**THETA and V3ALTBL). The ",
        "authors used a power rather than exponential model because the ALT ",
        "distribution is highly positively skewed (Methods 2.2.3). Mean ",
        "38.3 U/L (SD 35.1) in the PK analysis set (Table S2)."
      ),
      source_name = "ALTBL"
    ),
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Exponential effect centred at 44 years on Vc/F (Supplementary Datafile S2 $PK, V2AGEYBL). Mean 43.9 years (SD 16.0) (Table S2).",
      source_name = "AGEYBL"
    ),
    BMI = list(
      description = "Baseline body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Exponential effects centred at 26.44 kg/m^2 on CL/F and Vp/F (Supplementary Datafile S2 $PK, CLBMIBL / V3BMIBL). Mean 27.7 kg/m^2 (SD 5.9) (Table S2).",
      source_name = "BMIBL"
    ),
    CRCL = list(
      description = "Baseline creatinine clearance, body-surface-area normalised",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Exponential effect centred at 97.49 mL/min/1.73 m^2 on CL/F ",
        "(Supplementary Datafile S2 $PK, CLEGFRBL; the source column is ",
        "named EGFRBL but Table S3 and Figure 3 label the covariate CRCL). ",
        "Mean 96.6 mL/min/1.73 m^2 (SD 19.8), missing in 10.1% (Table S2). ",
        "Figure 3B evaluates 15, 30, 60, 90 and 120 mL/min/1.73 m^2."
      ),
      source_name = "EGFRBL"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male, the most common category in the PK analysis set)",
      notes = paste0(
        "Fractional difference on Q/F: 9.23% higher in women (Table S3, ",
        "'Female sex on Q/F (proportional increase) 0.0923'). The source ",
        "column SEXN codes 1 = male (reference) and 2 = female, so ",
        "SEXF = SEXN - 1. Note the sign is OPPOSITE to the parent model, ",
        "where female sex lowers Q/F by 18.3%."
      ),
      source_name = "SEXN"
    ),
    DIS_PBC = list(
      description = "Primary biliary cholangitis patient indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer or non-PBC patient population)",
      notes = "Fractional difference on Q/F: 21.4% higher in patients with PBC (Table S3, 'PBC on Q/F (proportional increase) 0.214'). Corresponds to POPN2 in Supplementary Datafile S2.",
      source_name = "POPN2"
    ),
    FED = list(
      description = "Fed state at the time of dosing",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted, 98.9% of baseline records in Table S2)",
      notes = paste0(
        "Fractional differences on D1 and relative bioavailability: the ",
        "zero-order formation duration is 11.6-fold longer (a 1160% ",
        "proportional increase) and relative bioavailability 25.2% higher ",
        "when the dose is taken fed (Table S3). Time-varying."
      ),
      source_name = "FOODN"
    ),
    FORM_ELA_F1 = list(
      description = "Elafibranor formulation 1 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all six formulation indicators 0 selects formulation 5, the phase II/III formulation)",
      notes = "Table S1 / Table S3 formulation numbering.",
      source_name = "FORMN == 1"
    ),
    FORM_ELA_F2 = list(
      description = "Elafibranor formulation 2 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (formulation 5 reference)",
      notes = "Table S1 / Table S3 formulation numbering.",
      source_name = "FORMN == 2"
    ),
    FORM_ELA_F3 = list(
      description = "Elafibranor formulation 3 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (formulation 5 reference)",
      notes = "Table S1 / Table S3 formulation numbering. GFT1007 relative bioavailability for formulation 3 is also fixed to 1.00, so formulations 3 and 5 share the reference bioavailability but differ in D1, MAT and lag time.",
      source_name = "FORMN == 3"
    ),
    FORM_ELA_F4 = list(
      description = "Elafibranor formulation 4 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (formulation 5 reference)",
      notes = "Table S1 / Table S3 formulation numbering; the phase IIb MASH formulation.",
      source_name = "FORMN == 4"
    ),
    FORM_ELA_F6 = list(
      description = "Elafibranor formulation 6 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (formulation 5 reference)",
      notes = "Table S1 / Table S3 formulation numbering.",
      source_name = "FORMN == 6"
    ),
    STUDY_PHASE2 = list(
      description = "Phase II study-stratum indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase I when STUDY_PHASE3 is also 0)",
      notes = "Selects the residual-error magnitudes only (Table S3 residual-error rows are tabulated per study phase).",
      source_name = "PHASEN == 2"
    ),
    STUDY_PHASE3 = list(
      description = "Phase III study-stratum indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase I when STUDY_PHASE2 is also 0)",
      notes = "Selects the residual-error magnitudes only (Table S3 residual-error rows are tabulated per study phase).",
      source_name = "PHASEN == 3"
    ),
    STUDY_GFT505B_319_1 = list(
      description = "ELATIVE phase III study (GFT505B-319-1, NCT04526665) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any other pooled study)",
      notes = "The absorption lag time was estimated separately for this study and fixed to 0 h (Table S3; Supplementary Datafile S2 $PK, IF(STUDYIDN.EQ.20) ALAG1 = THETA(41) with THETA(41) fixed at 0).",
      source_name = "STUDYIDN == 20"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "GFT1007",
      units = "umol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "GFT1007",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "GFT1007",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 894L,
    n_studies = 17L,
    n_observations = 10592L,
    age_range = "not reported; mean 43.9 years (SD 16.0)",
    age_median = "mean 43.9 years (SD 16.0) (Table S2)",
    weight_range = "not reported; mean 81.0 kg (SD 18.4)",
    weight_median = "mean 81.0 kg (SD 18.4) (Table S2)",
    sex_female_pct = 36.7,
    race_ethnicity = c(
      White = 84.7,
      `Black or African American` = 2.2,
      Asian = 0.8,
      `American Indian or Alaska Native` = 0.2,
      `Native Hawaiian or other Pacific Islander` = 0.2,
      `Multiple or other` = 1.2,
      `Unknown or not reported` = 0.6,
      Missing = 10.1
    ),
    disease_state = paste0(
      "pooled healthy volunteers and patients with renal impairment, ",
      "hepatic impairment, metabolic dysfunction-associated steatohepatitis ",
      "(MASH) or primary biliary cholangitis (PBC)"
    ),
    dose_range = "elafibranor 5-360 mg once daily, single dose to daily dosing for over a year",
    regions = "not reported",
    renal_function = "mean baseline creatinine clearance 96.6 mL/min/1.73 m^2 (SD 19.8); a dedicated renal-impairment study (GFT505-118-13) was included",
    hepatic_function = "NCI hepatic-impairment score 0 in 64.0%, 1 in 27.5%, 2 in 1.6%, 3 in 0.1%; a dedicated hepatic-impairment study (GFT505-118-14) was included",
    iov_structure = "None. Graphical analysis did not identify inter-occasion variability in the GFT1007 absorption parameters, so no IOV was incorporated (Results 3.1.5).",
    notes = paste0(
      "Baseline characteristics from Table S2; the GFT1007 PK data set held ",
      "894 subjects and 10,592 observations (Results section 3.1.1). Study ",
      "list in Table S1."
    )
  )

  # Six residual-error magnitudes, one per (study phase x time-since-last-dose
  # window) stratum in Table S3, combined into the single symbol expSdCc
  # inside model(). Declared here so checkModelConventions() does not read the
  # stratum suffixes as deviant residual-error names.
  paper_specific_residual_sds <- c(
    "expSdP1Late",
    "expSdP1Early",
    "expSdP2Late",
    "expSdP2Early",
    "expSdP3Late",
    "expSdP3Early"
  )

  # Inter-individual variability on the residual error itself (the
  # EPS(k)*EXP(ETA(j)) construction of Supplementary Datafile S2 $ERROR).
  paper_specific_etas <- c(
    "etaruvP1Late",
    "etaruvP1Early",
    "etaruvP2Late",
    "etaruvP2Early",
    "etaruvP3Late"
  )

  ini({
    # ==================================================================
    # DISPOSITION -- Table S3, GFT1007 column
    # ==================================================================
    lcl <- log(11.2)
    label("Apparent clearance CL/F (L/h)")                     # Table S3 'CL/F (L/h)' 11.2 [RSE 1.51%]
    lvc <- log(13.6)
    label("Apparent central volume Vc/F (L)")                  # Table S3 'Vc/F (L)' 13.6 [RSE 4.59%]
    lq <- log(6.83)
    label("Apparent inter-compartmental clearance Q/F (L/h)")  # Table S3 'Q/F (L/h)' 6.83 [RSE 2.02%]
    lvp <- log(77.8)
    label("Apparent peripheral volume Vp/F (L)")               # Table S3 'Vp/F (L)' 77.8 [RSE 2.79%]

    e_wt_cl_q <- fixed(0.750)
    label("Allometric exponent of baseline weight on CL/F and Q/F (unitless)")   # Table S3 'Baseline WT on CL (power)' 0.750 (FIX)
    e_wt_vc_vp <- fixed(1.00)
    label("Allometric exponent of baseline weight on Vc/F and Vp/F (unitless)")  # Table S3 'Baseline WT on V (power)' 1.00 (FIX)

    # ==================================================================
    # FORMULATION-SPECIFIC FORMATION / ABSORPTION -- Table S3
    # Formulation 5 is the reference; GFT1007 relative bioavailability is
    # fixed to 1.00 for formulations 3 AND 5.
    # ==================================================================
    lfdepot_f1 <- log(1.28)
    label("Relative bioavailability, formulation 1 (fraction)")     # Table S3 'Frel formulation 1' 1.28 [RSE 2.53%]
    lfdepot_f2 <- log(1.18)
    label("Relative bioavailability, formulation 2 (fraction)")     # Table S3 'Frel formulation 2' 1.18 [RSE 3.56%]
    lfdepot_f3 <- fixed(log(1.00))
    label("Relative bioavailability, formulation 3 (fraction)")     # Table S3 'Frel formulation 3' 1.00 (FIX)
    lfdepot_f4 <- log(1.14)
    label("Relative bioavailability, formulation 4 (fraction)")     # Table S3 'Frel formulation 4' 1.14 [RSE 3.58%]
    lfdepot_f5 <- fixed(log(1.00))
    label("Relative bioavailability, formulation 5 (fraction, reference)")  # Table S3 'Frel formulation 5' 1.00 (FIX)
    lfdepot_f6 <- log(0.836)
    label("Relative bioavailability, formulation 6 (fraction)")     # Table S3 'Frel formulation 6' 0.836 [RSE 5.48%]

    ld1_f1 <- log(0.216)
    label("Zero-order input duration D1, formulation 1 (h)")        # Table S3 'D1 formulation 1 (h)' 0.216 [RSE 11.6%]
    ld1_f2 <- log(0.511)
    label("Zero-order input duration D1, formulation 2 (h)")        # Table S3 'D1 formulation 2 (h)' 0.511 [RSE 11.4%]
    ld1_f3 <- log(0.779)
    label("Zero-order input duration D1, formulation 3 (h)")        # Table S3 'D1 formulation 3 (h)' 0.779 [RSE 6.97%]
    ld1_f4 <- log(0.716)
    label("Zero-order input duration D1, formulation 4 (h)")        # Table S3 'D1 formulation 4 (h)' 0.716 [RSE 8.39%]
    ld1_f5 <- log(0.294)
    label("Zero-order input duration D1, formulation 5 (h)")        # Table S3 'D1 formulation 5 (h)' 0.294 [RSE 5.78%]
    ld1_f6 <- log(0.420)
    label("Zero-order input duration D1, formulation 6 (h)")        # Table S3 'D1 formulation 6 (h)' 0.420 [RSE 33.0%]

    lmat_f1 <- log(1.17)
    label("Mean absorption time, formulation 1 (h)")                # Table S3 'MAT formulation 1 (h)' 1.17 [RSE 2.80%]
    lmat_f2 <- log(1.15)
    label("Mean absorption time, formulation 2 (h)")                # Table S3 'MAT formulation 2 (h)' 1.15 [RSE 5.29%]
    lmat_f3 <- log(1.16)
    label("Mean absorption time, formulation 3 (h)")                # Table S3 'MAT formulation 3 (h)' 1.16 [RSE 3.15%]
    lmat_f4 <- log(1.15)
    label("Mean absorption time, formulation 4 (h)")                # Table S3 'MAT formulation 4 (h)' 1.15 [RSE 5.74%]
    lmat_f5 <- log(1.25)
    label("Mean absorption time, formulation 5 (h)")                # Table S3 'MAT formulation 5 (h)' 1.25 [RSE 2.85%]
    lmat_f6 <- log(1.46)
    label("Mean absorption time, formulation 6 (h)")                # Table S3 'MAT formulation 6 (h)' 1.46 [RSE 12.3%]

    ltlag_f1 <- fixed(log(0.170))
    label("Absorption lag time, formulation 1 (h)")                 # Table S3 'Lag-time formulation 1 (h)' 0.170 (FIX)
    ltlag_f2 <- fixed(log(0.170))
    label("Absorption lag time, formulation 2 (h)")                 # Table S3 'Lag-time formulation 2 (h)' 0.170 (FIX)
    ltlag_f3 <- log(0.133)
    label("Absorption lag time, formulation 3 (h)")                 # Table S3 'Lag-time formulation 3 (h)' 0.133 [RSE 3.57%]
    ltlag_f4 <- log(0.155)
    label("Absorption lag time, formulation 4 (h)")                 # Table S3 'Lag-time formulation 4 (h)' 0.155 [RSE 7.38%]
    ltlag_f5 <- log(0.137)
    label("Absorption lag time, formulation 5 (h)")                 # Table S3 'Lag-time formulation 5 (h)' 0.137 [RSE 1.03%]
    ltlag_f6 <- log(0.122)
    label("Absorption lag time, formulation 6 (h)")                 # Table S3 'Lag-time formulation 6 (h)' 0.122 [RSE 5.01%]

    # ==================================================================
    # COVARIATE EFFECTS -- Table S3
    # ==================================================================
    e_alt_cl <- -0.0489
    label("Power exponent of baseline ALT on CL/F (unitless)")                # Table S3 'Baseline ALT on CL/F (power)' -0.0489 [RSE 32.8%]
    e_alt_vp <- -0.123
    label("Power exponent of baseline ALT on Vp/F (unitless)")                # Table S3 'Baseline ALT on Vp/F (power)' -0.123 [RSE 35.0%]
    e_fed_d1 <- 11.6
    label("Proportional increase in D1 in the fed state (fraction)")          # Table S3 'Food on D1 (proportional increase)' 11.6 [RSE 17.2%]
    e_fed_fdepot <- 0.252
    label("Proportional increase in relative bioavailability in the fed state (fraction)")  # Table S3 'Food on Frel (proportional increase)' 0.252 [RSE 19.5%]
    e_age_vc <- 0.00487
    label("Exponential coefficient of baseline age on Vc/F (per year)")       # Table S3 'Baseline age on Vc/F (exponential)' 0.00487 [RSE 41.2%]
    e_bmi_cl <- -0.00982
    label("Exponential coefficient of baseline BMI on CL/F (per kg/m^2)")     # Table S3 'Baseline BMI on CL/F (exponential)' -0.00982 [RSE 19.0%]
    e_bmi_vp <- -0.0227
    label("Exponential coefficient of baseline BMI on Vp/F (per kg/m^2)")     # Table S3 'Baseline BMI on Vp/F (exponential)' -0.0227 [RSE 21.8%]
    e_crcl_cl <- 0.00313
    label("Exponential coefficient of baseline CRCL on CL/F (per mL/min/1.73 m^2)")  # Table S3 'Baseline CRCL on CL/F (exponential)' 0.00313 [RSE 11.6%]
    e_sexf_q <- 0.0923
    label("Proportional difference in Q/F for female sex (fraction)")         # Table S3 'Female sex on Q/F (proportional increase)' 0.0923 [RSE 47.5%]
    e_pbc_q <- 0.214
    label("Proportional difference in Q/F in the PBC population (fraction)")  # Table S3 'PBC on Q/F (proportional increase)' 0.214 [RSE 52.9%]

    # ==================================================================
    # INTER-INDIVIDUAL VARIABILITY -- Table S3, GFT1007 column
    # Reported on the approximate SD (CV) scale; variances are its square.
    # ==================================================================
    etalfdepot_f1 ~ 0.017161   # Table S3 'IIV Frel formulation 1 (CV)' 0.131; 0.131^2
    # IIV on Frel is 0 (FIX) for formulations 2 and 4.
    etalfdepot_f3 ~ 0.045796   # Table S3 'IIV Frel formulation 3 (CV)' 0.214; 0.214^2
    etalfdepot_f5 ~ 0.078400   # Table S3 'IIV Frel formulation 5 (CV)' 0.280; 0.280^2
    etalfdepot_f6 ~ 0.025281   # Table S3 'IIV Frel formulation 6 (CV)' 0.159; 0.159^2

    # Formulations 1-4 share a single D1 eta in Supplementary Datafile S2
    # (IF(FORMN.GE.3) ... ETA(21), the slot labelled 'IIV on D1 FORMCAP'),
    # which is why Table S3 prints the same 0.976 for all four.
    etald1_f1to4 ~ 0.952576    # Table S3 'IIV D1 formulation 1-4 (CV)' 0.976; 0.976^2
    etald1_f5 ~ 2.689600       # Table S3 'IIV D1 formulation 5 (CV)' 1.64; 1.64^2
    etald1_f6 ~ 1.587600       # Table S3 'IIV D1 formulation 6 (CV)' 1.26; 1.26^2

    # IIV on MAT was estimated for formulation 6 only; fixed to 0 elsewhere.
    etalmat_f6 ~ 0.261121      # Table S3 'IIV MAT formulation 6 (CV)' 0.511; 0.511^2

    etalcl ~ 0.020164          # Table S3 'IIV CL/F (CV)' 0.142; 0.142^2
    etalvc ~ 0.293764          # Table S3 'IIV Vc/F (CV)' 0.542; 0.542^2
    etalq ~ 0.032041           # Table S3 'IIV Q/F (CV)' 0.179; 0.179^2
    etalvp ~ 0.171396          # Table S3 'IIV Vp/F (CV)' 0.414; 0.414^2

    # ==================================================================
    # RESIDUAL ERROR -- Table S3, GFT1007 column
    # Exponential (additive on the log scale), stratified by study phase and
    # by whether the sample was drawn within 2.5 h of the last dose. Note the
    # 2.5 h split, versus 5 h for the parent.
    # ==================================================================
    expSdP1Late <- 0.251
    label("Residual SD, phase I, >2.5 h since last dose (log scale)")     # Table S3 'RUV Phase I TSLD>2.5h (CV)' 0.251 [RSE 1.94%]
    expSdP1Early <- 0.421
    label("Residual SD, phase I, <=2.5 h since last dose (log scale)")    # Table S3 'RUV Phase I TSLD<=2.5h (CV)' 0.421 [RSE 2.82%]
    expSdP2Late <- 0.429
    label("Residual SD, phase II, >2.5 h since last dose (log scale)")    # Table S3 'RUV Phase II TSLD>2.5h (CV)' 0.429 [RSE 7.58%]
    expSdP2Early <- 0.448
    label("Residual SD, phase II, <=2.5 h since last dose (log scale)")   # Table S3 'RUV Phase II TSLD<=2.5h (CV)' 0.448 [RSE 7.22%]
    expSdP3Late <- 0.496
    label("Residual SD, phase III, >2.5 h since last dose (log scale)")   # Table S3 'RUV Phase III TSLD>2.5h (CV)' 0.496 [RSE 6.56%]
    expSdP3Early <- 0.525
    label("Residual SD, phase III, <=2.5 h since last dose (log scale)")  # Table S3 'RUV Phase III TSLD<=2.5h (CV)' 0.525 [RSE 3.14%]

    # IIV on the residual error. The phase III early stratum is 0 (FIX) and
    # therefore carries no eta.
    etaruvP1Late ~ 0.056169    # Table S3 'IIV RUV Phase I TSLD>2.5h (CV)' 0.237; 0.237^2
    etaruvP1Early ~ 0.160801   # Table S3 'IIV RUV Phase I TSLD<=2.5h (CV)' 0.401; 0.401^2
    etaruvP2Late ~ 0.067081    # Table S3 'IIV RUV Phase II TSLD>2.5h (CV)' 0.259; 0.259^2
    etaruvP2Early ~ 0.250000   # Table S3 'IIV RUV Phase II TSLD<=2.5h (CV)' 0.500; 0.500^2
    etaruvP3Late ~ 0.110889    # Table S3 'IIV RUV Phase III TSLD>2.5h (CV)' 0.333; 0.333^2
  })

  model({
    # ------------------------------------------------------------------
    # 1. Formulation selection (formulation 5 when all indicators are 0)
    # ------------------------------------------------------------------
    isF5 <- 1 - FORM_ELA_F1 - FORM_ELA_F2 - FORM_ELA_F3 - FORM_ELA_F4 - FORM_ELA_F6

    lfdepotSel <-
      lfdepot_f1 * FORM_ELA_F1 + lfdepot_f2 * FORM_ELA_F2 +
      lfdepot_f3 * FORM_ELA_F3 + lfdepot_f4 * FORM_ELA_F4 +
      lfdepot_f5 * isF5 + lfdepot_f6 * FORM_ELA_F6
    # IIV on Frel is 0 (FIX) for formulations 2 and 4.
    etalfdepotSel <-
      etalfdepot_f1 * FORM_ELA_F1 + etalfdepot_f3 * FORM_ELA_F3 +
      etalfdepot_f5 * isF5 + etalfdepot_f6 * FORM_ELA_F6

    ld1Sel <-
      ld1_f1 * FORM_ELA_F1 + ld1_f2 * FORM_ELA_F2 + ld1_f3 * FORM_ELA_F3 +
      ld1_f4 * FORM_ELA_F4 + ld1_f5 * isF5 + ld1_f6 * FORM_ELA_F6
    etald1Sel <-
      etald1_f1to4 * (FORM_ELA_F1 + FORM_ELA_F2 + FORM_ELA_F3 + FORM_ELA_F4) +
      etald1_f5 * isF5 + etald1_f6 * FORM_ELA_F6

    lmatSel <-
      lmat_f1 * FORM_ELA_F1 + lmat_f2 * FORM_ELA_F2 + lmat_f3 * FORM_ELA_F3 +
      lmat_f4 * FORM_ELA_F4 + lmat_f5 * isF5 + lmat_f6 * FORM_ELA_F6
    # IIV on MAT is 0 (FIX) for formulations 1-5.
    etalmatSel <- etalmat_f6 * FORM_ELA_F6

    ltlagSel <-
      ltlag_f1 * FORM_ELA_F1 + ltlag_f2 * FORM_ELA_F2 + ltlag_f3 * FORM_ELA_F3 +
      ltlag_f4 * FORM_ELA_F4 + ltlag_f5 * isF5 + ltlag_f6 * FORM_ELA_F6

    # ------------------------------------------------------------------
    # 2. Individual formation / absorption parameters
    # ------------------------------------------------------------------
    fdepot <- exp(lfdepotSel + etalfdepotSel) * (1 + e_fed_fdepot * FED)
    d1 <- exp(ld1Sel + etald1Sel) * (1 + e_fed_d1 * FED)
    mat <- exp(lmatSel + etalmatSel)
    ka <- 1 / mat
    # The lag time is fixed to 0 in the ELATIVE phase III study.
    tlag <- exp(ltlagSel) * (1 - STUDY_GFT505B_319_1)

    # ------------------------------------------------------------------
    # 3. Individual disposition parameters
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl_q * (ALT / 27)^e_alt_cl *
      exp(e_bmi_cl * (BMI - 26.44)) * exp(e_crcl_cl * (CRCL - 97.49))
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc_vp * exp(e_age_vc * (AGE - 44))
    q <- exp(lq + etalq) * (WT / 75)^e_wt_cl_q *
      (1 + e_sexf_q * SEXF) * (1 + e_pbc_q * DIS_PBC)
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vc_vp * (ALT / 27)^e_alt_vp *
      exp(e_bmi_vp * (BMI - 26.44))

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ------------------------------------------------------------------
    # 4. ODE system -- NONMEM ADVAN4 TRANS4 with a zero-order input into the
    #    depot (Datafile S2 $SUBROUTINE / $PK). The depot represents GFT1007
    #    formation from the administered elafibranor dose.
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central +
      k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    dur(depot) <- d1
    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # 5. Observation and residual error
    # ------------------------------------------------------------------
    Cc <- central / vc

    tsld <- tad()
    isEarly <- (tsld <= 2.5)
    isLate <- 1 - isEarly
    isPh2 <- STUDY_PHASE2
    isPh3 <- STUDY_PHASE3
    isPh1 <- 1 - isPh2 - isPh3

    expSdSel <-
      isPh1 * (isEarly * expSdP1Early + isLate * expSdP1Late) +
      isPh2 * (isEarly * expSdP2Early + isLate * expSdP2Late) +
      isPh3 * (isEarly * expSdP3Early + isLate * expSdP3Late)

    # The phase III early stratum has IIV on the residual error fixed to 0.
    etaruvSel <-
      isPh1 * (isEarly * etaruvP1Early + isLate * etaruvP1Late) +
      isPh2 * (isEarly * etaruvP2Early + isLate * etaruvP2Late) +
      isPh3 * isLate * etaruvP3Late

    expSdCc <- expSdSel * exp(etaruvSel)

    Cc ~ lnorm(expSdCc)
  })
}
