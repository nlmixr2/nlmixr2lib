Ooi_2026_elafibranor <- function() {
  description <- paste0(
    "Two-compartment population pharmacokinetic model for the parent drug ",
    "elafibranor, a PPAR-alpha/delta agonist approved for second-line ",
    "treatment of primary biliary cholangitis (PBC), pooled over 17 clinical ",
    "trials (13 phase I, 2 phase II, 2 phase III; 892 subjects, 12,205 ",
    "observations; doses 5-360 mg; healthy volunteers and patients with ",
    "renal impairment, hepatic impairment, MASH or PBC). Oral absorption is ",
    "sequential zero-order input into the depot over a duration D1, followed ",
    "by first-order transfer at ka = 1/MAT, after an absorption lag time. ",
    "Relative bioavailability, D1, MAT and the lag time are all ",
    "formulation-specific over six clinical formulations (formulation 5, the ",
    "phase II/III formulation, is the reference with Frel fixed to 1). ",
    "Relative bioavailability additionally rises with dose through a ",
    "sigmoidal Emax function of the milligram dose, which reproduces the ",
    "less-than-proportional exposure observed below 50 mg. Apparent ",
    "clearance and apparent inter-compartmental clearance scale ",
    "allometrically with baseline body weight at a fixed exponent of 0.75 ",
    "and the two apparent volumes at 1.0, referenced to 75 kg. Additional ",
    "covariates are baseline albumin on CL/F and Q/F, age on Vp/F, D1 and ",
    "Frel, female sex on Q/F, baseline BMI on Vp/F, PBC population on Vp/F, ",
    "and fed state on MAT. Two bioanalytical-method indicators shift the ",
    "predicted concentration multiplicatively. Inter-individual variability ",
    "is diagonal and, for the absorption parameters, formulation-specific; ",
    "inter-occasion variability is carried on D1 and MAT over ten occasions. ",
    "The residual error is exponential (additive on the log scale) with ",
    "twelve magnitudes selected by bioanalytical method, study phase and ",
    "whether the sample was drawn within 5 h of the last dose. The active ",
    "metabolite is modelled separately in ",
    "modellib('Ooi_2026_elafibranor_gft1007'); the exposure-response model ",
    "that consumes both is modellib('Ooi_2026_elafibranor_alptb')."
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

  # The NONMEM data set carried amounts in nmol and volumes in mL (the
  # LOG(1000) terms in the Supplementary Datafile S1 $PK block convert the
  # tabulated L and L/h values to mL and mL/h), so predicted concentrations
  # are in nmol/mL. This file keeps the tabulated L / L/h values and doses in
  # umol, which gives numerically identical concentrations in umol/L
  # (1 nmol/mL == 1 umol/L). A milligram dose is converted with the
  # elafibranor molecular weight 384.49 g/mol, which the authors use
  # explicitly in Supplementary Datafile S6 ($PK, AUCSSP line): 80 mg is
  # 80 / 384.49 * 1000 = 208.07 umol.
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Baseline (not time-varying) body weight, included a priori as a ",
        "mechanistic covariate with allometric scaling referenced to 75 kg ",
        "(Supplementary Datafile S1 $PK, LOG(WTKGBL/75)). Exponents fixed at ",
        "0.75 on CL/F and Q/F and 1.00 on Vc/F and Vp/F (Table S3). Mean ",
        "81.0 kg (SD 18.4) in the PK analysis set (Table S2)."
      ),
      source_name = "WTKGBL"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Exponential effect centred at 43.30 g/L on both CL/F and Q/F ",
        "(Supplementary Datafile S1 $PK, CLALBBL / QALBBL blocks). Mean ",
        "43.3 g/L (SD 3.5) in the PK analysis set (Table S2). The control ",
        "stream imputes missing albumin as 42 g/L in healthy volunteers and ",
        "46 g/L in the hepatic-impairment population; that imputation is not ",
        "reproduced here because it is a data-handling rule, not part of the ",
        "structural model."
      ),
      source_name = "ALBBL"
    ),
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Exponential effects centred at 44 years on Vp/F, D1 and Frel ",
        "(Supplementary Datafile S1 $PK, V3AGEYBL / D1AGEYBL / FRELAGEYBL). ",
        "Mean 43.9 years (SD 16.0) in the PK analysis set (Table S2)."
      ),
      source_name = "AGEYBL"
    ),
    BMI = list(
      description = "Baseline body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Exponential effect centred at 26.44 kg/m^2 on Vp/F (Supplementary ",
        "Datafile S1 $PK, V3BMIBL). Mean 27.7 kg/m^2 (SD 5.9) in the PK ",
        "analysis set (Table S2)."
      ),
      source_name = "BMIBL"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male, the most common category in the PK analysis set)",
      notes = paste0(
        "Fractional difference on Q/F relative to the most common category: ",
        "Q/F is 18.3% lower in women (Table S3, 'Female sex on Q/F ",
        "(proportional increase) -0.183'). The source column SEXN codes ",
        "1 = male (the reference) and 2 = female, so SEXF = SEXN - 1. ",
        "36.7% of the PK analysis set were female (Table S2)."
      ),
      source_name = "SEXN"
    ),
    DIS_PBC = list(
      description = "Primary biliary cholangitis patient indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer or non-PBC patient population)",
      notes = paste0(
        "Fractional difference on Vp/F: 28.4% lower apparent peripheral ",
        "volume in patients with PBC (Table S3, 'PBC on Vp/F ",
        "(proportional increase) -0.284'). Corresponds to the POPN2 ",
        "indicator of Supplementary Datafile S1."
      ),
      source_name = "POPN2"
    ),
    FED = list(
      description = "Fed state at the time of dosing",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted, 98.9% of baseline records in Table S2)",
      notes = paste0(
        "Fractional difference on MAT: mean absorption time is 155% higher ",
        "(2.55-fold) when the dose is taken fed (Table S3, 'Food on MAT ",
        "(proportional increase) 1.55'). Time-varying: a subject may be ",
        "dosed fed on some occasions and fasted on others."
      ),
      source_name = "FOODN"
    ),
    DOSE_ELA_MG = list(
      description = "Administered elafibranor dose",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Milligram dose of the record, entering the sigmoidal Emax function ",
        "that raises relative bioavailability with dose (Supplementary ",
        "Datafile S1 $PK, DOSE_FREL). Must be the MILLIGRAM dose even though ",
        "amounts in the event table are in umol, because ED50 = 11.0 mg is ",
        "reported on the milligram scale (Table S3). Studied doses were ",
        "5-360 mg (Table S1)."
      ),
      source_name = "DOSEN"
    ),
    FORM_ELA_F1 = list(
      description = "Elafibranor formulation 1 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all six formulation indicators 0 selects formulation 5, the phase II/III formulation)",
      notes = "Table S1 / Table S3 formulation numbering. Used in the phase I studies GFT505-106-1, GFT505-106-2, GFT505-108-4 and GFT505-108-3.",
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
      notes = "Table S1 / Table S3 formulation numbering.",
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
    ASSAY_SEPIP = list(
      description = "Bioanalytical method with separation of the interfering peak only",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (neither separation of the interfering peak nor addition of formic acid)",
      notes = paste0(
        "Per-observation indicator. Shifts the predicted concentration ",
        "multiplicatively by exp(-0.182) and selects its own residual-error ",
        "magnitudes (Table S3). Mutually exclusive with ASSAY_SEPIP_FA. ",
        "24.6% of baseline records (Table S2)."
      ),
      source_name = "BIOANN == 2"
    ),
    ASSAY_SEPIP_FA = list(
      description = "Bioanalytical method with separation of the interfering peak AND addition of formic acid",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (neither separation of the interfering peak nor addition of formic acid)",
      notes = paste0(
        "Per-observation indicator. Shifts the predicted concentration ",
        "multiplicatively by exp(-0.760) and selects its own residual-error ",
        "magnitudes (Table S3). Mutually exclusive with ASSAY_SEPIP. ",
        "3.4% of baseline records carried formic acid (Table S2)."
      ),
      source_name = "BIOANN %in% c(3, 4)"
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
      notes = paste0(
        "The absorption lag time was estimated separately for this study and ",
        "fixed to 0 h (Table S3, 'Lag-time for Study GFT505B-319-1 (h) ",
        "0 (FIX)'; Supplementary Datafile S1 $PK, IF(STUDYIDN.EQ.20) ",
        "ALAG1 = THETA(54) with THETA(54) fixed at 0). Set to 1 to reproduce ",
        "the ELATIVE phase III profiles."
      ),
      source_name = "STUDYIDN == 20"
    ),
    OCC = list(
      description = "Dosing-occasion index for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste0(
        "Integer 1-10. Supplementary Datafile S1 carries ten occasion slots ",
        "for each of the two IOV-bearing absorption parameters ",
        "($ABBREVIATED REPLACE ETA(OCC_MAT)=ETA(14,...,23) and ",
        "ETA(OCC_D1)=ETA(24,...,33)), all held equal through $OMEGA BLOCK(1) ",
        "SAME. Set OCC to the index of the dosing interval; a ",
        "single-occasion simulation may use OCC = 1 throughout, which leaves ",
        "one IOV draw per subject. Occasions beyond 10 are not defined by ",
        "the source."
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "elafibranor",
      units = "umol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "elafibranor",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "elafibranor",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 892L,
    n_studies = 17L,
    n_observations = 12205L,
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
    dose_range = "5-360 mg once daily, single dose to daily dosing for over a year",
    regions = "not reported",
    renal_function = "mean baseline creatinine clearance 96.6 mL/min/1.73 m^2 (SD 19.8); a dedicated renal-impairment study (GFT505-118-13) was included",
    hepatic_function = "NCI hepatic-impairment score 0 in 64.0%, 1 in 27.5%, 2 in 1.6%, 3 in 0.1%; a dedicated hepatic-impairment study (GFT505-118-14) was included",
    iov_structure = "Inter-occasion variability on D1 (93.4% CV) and MAT (41.4% CV) over ten occasions; encoded here via the OCC covariate.",
    notes = paste0(
      "Baseline characteristics from Table S2 (which reports N = 894, the ",
      "union of the elafibranor and GFT1007 analysis sets); the elafibranor ",
      "PK data set itself held 892 subjects and 12,205 observations ",
      "(Results section 3.1.1). Study list in Table S1."
    )
  )

  # Twelve residual-error magnitudes, one per (bioanalytical method x study
  # phase x time-since-last-dose window) stratum reported in Table S3, are
  # combined into the single symbol expSdCc inside model(). Declared here so
  # checkModelConventions() does not read the stratum suffixes as deviant
  # residual-error names (the Rich_2026_momelotinib.R / Friberg_2012_voriconazole.R
  # pattern).
  paper_specific_residual_sds <- c(
    "expSdP1NoneLate",
    "expSdP1NoneEarly",
    "expSdP1SepLate",
    "expSdP1SepEarly",
    "expSdP1SepfaLate",
    "expSdP1SepfaEarly",
    "expSdP2Late",
    "expSdP2Early",
    "expSdP3NoneLate",
    "expSdP3NoneEarly",
    "expSdP3SepLate",
    "expSdP3SepEarly"
  )

  # Inter-individual variability on the residual error itself (the
  # EPS(k)*EXP(ETA(j)) construction of Supplementary Datafile S1 $ERROR) and
  # the per-occasion IOV slots.
  paper_specific_etas <- c(
    "etaruvP1NoneLate",
    "etaruvP1NoneEarly",
    "etaruvP1SepLate",
    "etaruvP1SepEarly",
    "etaruvP1SepfaLate",
    "etaruvP2Late",
    "etaruvP3NoneLate",
    paste0("etaiov_mat_", 1:10),
    paste0("etaiov_d1_", 1:10)
  )

  ini({
    # ==================================================================
    # DISPOSITION -- Table S3, elafibranor column
    # ==================================================================
    lcl <- log(47.1)
    label("Apparent clearance CL/F (L/h)")                         # Table S3 'CL/F (L/h)' 47.1 [RSE 11.9%]
    lvc <- log(68.7)
    label("Apparent central volume Vc/F (L)")                      # Table S3 'Vc/F (L)' 68.7 [RSE 15.5%]
    lq <- log(263)
    label("Apparent inter-compartmental clearance Q/F (L/h)")      # Table S3 'Q/F (L/h)' 263 [RSE 12.0%]
    lvp <- log(4310)
    label("Apparent peripheral volume Vp/F (L)")                   # Table S3 'Vp/F (L)' 4310 [RSE 12.1%]

    e_wt_cl_q <- fixed(0.750)
    label("Allometric exponent of baseline weight on CL/F and Q/F (unitless)")   # Table S3 'Baseline WT on CL (power)' 0.750 (FIX)
    e_wt_vc_vp <- fixed(1.00)
    label("Allometric exponent of baseline weight on Vc/F and Vp/F (unitless)")  # Table S3 'Baseline WT on V (power)' 1.00 (FIX)

    # ==================================================================
    # FORMULATION-SPECIFIC ABSORPTION -- Table S3, elafibranor column
    # Formulation 5 (the phase II/III formulation) is the reference.
    # ==================================================================
    lfdepot_f1 <- log(1.02)
    label("Relative bioavailability, formulation 1 (fraction)")     # Table S3 'Frel formulation 1' 1.02 [RSE 5.73%]
    lfdepot_f2 <- log(1.13)
    label("Relative bioavailability, formulation 2 (fraction)")     # Table S3 'Frel formulation 2' 1.13 [RSE 6.38%]
    lfdepot_f3 <- log(1.02)
    label("Relative bioavailability, formulation 3 (fraction)")     # Table S3 'Frel formulation 3' 1.02 [RSE 5.16%]
    lfdepot_f4 <- log(1.07)
    label("Relative bioavailability, formulation 4 (fraction)")     # Table S3 'Frel formulation 4' 1.07 [RSE 5.21%]
    lfdepot_f5 <- fixed(log(1.00))
    label("Relative bioavailability, formulation 5 (fraction, reference)")  # Table S3 'Frel formulation 5' 1.00 (FIX)
    lfdepot_f6 <- log(0.979)
    label("Relative bioavailability, formulation 6 (fraction)")     # Table S3 'Frel formulation 6' 0.979 [RSE 11.0%]

    ld1_f1 <- log(0.388)
    label("Zero-order input duration D1, formulation 1 (h)")        # Table S3 'D1 formulation 1 (h)' 0.388 [RSE 13.4%]
    ld1_f2 <- log(0.440)
    label("Zero-order input duration D1, formulation 2 (h)")        # Table S3 'D1 formulation 2 (h)' 0.440 [RSE 27.7%]
    ld1_f3 <- log(0.728)
    label("Zero-order input duration D1, formulation 3 (h)")        # Table S3 'D1 formulation 3 (h)' 0.728 [RSE 9.92%]
    ld1_f4 <- log(0.0538)
    label("Zero-order input duration D1, formulation 4 (h)")        # Table S3 'D1 formulation 4 (h)' 0.0538 [RSE 35.2%]
    ld1_f5 <- log(0.187)
    label("Zero-order input duration D1, formulation 5 (h)")        # Table S3 'D1 formulation 5 (h)' 0.187 [RSE 4.31%]
    ld1_f6 <- log(0.131)
    label("Zero-order input duration D1, formulation 6 (h)")        # Table S3 'D1 formulation 6 (h)' 0.131 [RSE 56.1%]

    lmat_f1 <- log(0.629)
    label("Mean absorption time, formulation 1 (h)")                # Table S3 'MAT formulation 1 (h)' 0.629 [RSE 6.02%]
    lmat_f2 <- log(0.579)
    label("Mean absorption time, formulation 2 (h)")                # Table S3 'MAT formulation 2 (h)' 0.579 [RSE 13.5%]
    lmat_f3 <- log(0.646)
    label("Mean absorption time, formulation 3 (h)")                # Table S3 'MAT formulation 3 (h)' 0.646 [RSE 5.36%]
    lmat_f4 <- log(0.687)
    label("Mean absorption time, formulation 4 (h)")                # Table S3 'MAT formulation 4 (h)' 0.687 [RSE 11.6%]
    lmat_f5 <- log(0.911)
    label("Mean absorption time, formulation 5 (h)")                # Table S3 'MAT formulation 5 (h)' 0.911 [RSE 3.44%]
    lmat_f6 <- log(1.15)
    label("Mean absorption time, formulation 6 (h)")                # Table S3 'MAT formulation 6 (h)' 1.15 [RSE 9.89%]

    ltlag_f1 <- fixed(log(0.0450))
    label("Absorption lag time, formulation 1 (h)")                 # Table S3 'Lag-time formulation 1 (h)' 0.0450 (FIX)
    ltlag_f2 <- fixed(log(0.00100))
    label("Absorption lag time, formulation 2 (h)")                 # Table S3 'Lag-time formulation 2 (h)' 0.00100 (FIX)
    ltlag_f3 <- log(0.147)
    label("Absorption lag time, formulation 3 (h)")                 # Table S3 'Lag-time formulation 3 (h)' 0.147 [RSE 1.63%]
    ltlag_f4 <- fixed(log(0.170))
    label("Absorption lag time, formulation 4 (h)")                 # Table S3 'Lag-time formulation 4 (h)' 0.170 (FIX)
    ltlag_f5 <- log(0.137)
    label("Absorption lag time, formulation 5 (h)")                 # Table S3 'Lag-time formulation 5 (h)' 0.137 [RSE 1.33%]
    ltlag_f6 <- log(0.132)
    label("Absorption lag time, formulation 6 (h)")                 # Table S3 'Lag-time formulation 6 (h)' 0.132 [RSE 3.34%]

    # ==================================================================
    # DOSE EFFECT ON RELATIVE BIOAVAILABILITY -- Table S3
    # Frel is multiplied by 1 + Emax * DOSE^gamma / (ED50^gamma + DOSE^gamma),
    # which saturates near 1.8 above ~50 mg and reproduces the
    # less-than-proportional exposure the authors report below 50 mg
    # (Discussion).
    # ==================================================================
    emax_dose_fdepot <- 0.800
    label("Maximum proportional increase in relative bioavailability with dose (fraction)")  # Table S3 'Emax dose on Frel (proportional increase)' 0.800 [RSE 25.3%]
    ed50_dose_fdepot <- 11.0
    label("Dose giving half the maximal increase in relative bioavailability (mg)")          # Table S3 'ED50 dose on Frel (mg)' 11.0 [RSE 12.5%]
    gamma_dose_fdepot <- 2.31
    label("Hill coefficient of the dose effect on relative bioavailability (unitless)")      # Table S3 'Gamma dose on Frel' 2.31 [RSE 31.6%]

    # ==================================================================
    # COVARIATE EFFECTS -- Table S3
    # ==================================================================
    e_fed_mat <- 1.55
    label("Proportional increase in MAT in the fed state (fraction)")            # Table S3 'Food on MAT (proportional increase)' 1.55 [RSE 18.3%]
    e_alb_cl <- -0.0227
    label("Exponential coefficient of baseline albumin on CL/F (per g/L)")       # Table S3 'Baseline ALB on CL/F (exponential)' -0.0227 [RSE 17.5%]
    e_alb_q <- -0.0158
    label("Exponential coefficient of baseline albumin on Q/F (per g/L)")        # Table S3 'Baseline ALB on Q/F (exponential)' -0.0158 [RSE 25.0%]
    e_age_vp <- 0.00589
    label("Exponential coefficient of baseline age on Vp/F (per year)")          # Table S3 'Baseline age on Vp/F (exponential)' 0.00589 [RSE 13.9%]
    e_age_d1 <- -0.0172
    label("Exponential coefficient of baseline age on D1 (per year)")            # Table S3 'Baseline age on D1 (exponential)' -0.0172 [RSE 31.2%]
    e_age_fdepot <- 0.00332
    label("Exponential coefficient of baseline age on relative bioavailability (per year)")  # Table S3 'Baseline age on Frel (exponential)' 0.00332 [RSE 38.0%]
    e_sexf_q <- -0.183
    label("Proportional difference in Q/F for female sex (fraction)")            # Table S3 'Female sex on Q/F (proportional increase)' -0.183 [RSE 13.4%]
    e_bmi_vp <- -0.0151
    label("Exponential coefficient of baseline BMI on Vp/F (per kg/m^2)")        # Table S3 'Baseline BMI on Vp/F (exponential)' -0.0151 [RSE 16.5%]
    e_pbc_vp <- -0.284
    label("Proportional difference in Vp/F in the PBC population (fraction)")    # Table S3 'PBC on Vp/F (proportional increase)' -0.284 [RSE 20.9%]

    e_assay_sepip_cc <- -0.182
    label("Log-scale shift of the predicted concentration for separation of the interfering peak (unitless)")             # Table S3 'Separation of interfering peak (exponential)' -0.182 [RSE 26.7%]
    e_assay_sepip_fa_cc <- -0.760
    label("Log-scale shift of the predicted concentration for separation of the interfering peak plus formic acid (unitless)")  # Table S3 'Separation of interfering peak and addition of formic acid (exponential)' -0.760 [RSE 11.4%]

    # ==================================================================
    # INTER-INDIVIDUAL VARIABILITY -- Table S3
    # Table S3 reports IIV on the approximate SD (CV) scale; variances below
    # are the square of the tabulated value.
    # ==================================================================
    etalfdepot_f1 ~ 0.041616   # Table S3 'IIV Frel formulation 1 (CV)' 0.204; 0.204^2
    etalfdepot_f2 ~ 0.048400   # Table S3 'IIV Frel formulation 2 (CV)' 0.220; 0.220^2
    etalfdepot_f3 ~ 0.084681   # Table S3 'IIV Frel formulation 3 (CV)' 0.291; 0.291^2
    etalfdepot_f4 ~ 0.056169   # Table S3 'IIV Frel formulation 4 (CV)' 0.237; 0.237^2
    etalfdepot_f5 ~ 0.127449   # Table S3 'IIV Frel formulation 5 (CV)' 0.357; 0.357^2
    etalfdepot_f6 ~ 0.116964   # Table S3 'IIV Frel formulation 6 (CV)' 0.342; 0.342^2

    # IIV on D1 was fixed to 0 for formulations 1-3; those etas are therefore
    # absent rather than declared with zero variance.
    etald1_f4 ~ 12.1801        # Table S3 'IIV D1 formulation 4 (CV)' 3.49; 3.49^2
    etald1_f5 ~ 2.7225         # Table S3 'IIV D1 formulation 5 (CV)' 1.65; 1.65^2
    etald1_f6 ~ 4.0000         # Table S3 'IIV D1 formulation 6 (CV)' 2.00; 2.00^2

    # IIV on MAT was estimated for formulation 1 only; fixed to 0 elsewhere.
    etalmat_f1 ~ 0.075076      # Table S3 'IIV MAT formulation 1 (CV)' 0.274; 0.274^2

    etalcl ~ 0.055225          # Table S3 'IIV CL/F (CV)' 0.235; 0.235^2
    etalvc ~ 1.166400          # Table S3 'IIV Vc/F (CV)' 1.08; 1.08^2
    etalq ~ 0.022801           # Table S3 'IIV Q/F (CV)' 0.151; 0.151^2
    etalvp ~ 0.00595984        # Table S3 'IIV Vp/F (CV)' 0.0772; 0.0772^2

    # ==================================================================
    # INTER-OCCASION VARIABILITY -- Table S3 and Supplementary Datafile S1
    # Ten occasion slots per parameter, all constrained equal by
    # $OMEGA BLOCK(1) SAME, so occasions 2-10 are fixed to the occasion-1
    # variance.
    # ==================================================================
    etaiov_mat_1 ~ 0.170958    # Datafile S1 $OMEGA 14 0.170958; sqrt = 0.4135, Table S3 'IOV MAT (CV)' 0.414
    etaiov_mat_2 ~ fix(0.170958)   # $OMEGA BLOCK(1) SAME
    etaiov_mat_3 ~ fix(0.170958)   # $OMEGA BLOCK(1) SAME
    etaiov_mat_4 ~ fix(0.170958)   # $OMEGA BLOCK(1) SAME
    etaiov_mat_5 ~ fix(0.170958)   # $OMEGA BLOCK(1) SAME
    etaiov_mat_6 ~ fix(0.170958)   # $OMEGA BLOCK(1) SAME
    etaiov_mat_7 ~ fix(0.170958)   # $OMEGA BLOCK(1) SAME
    etaiov_mat_8 ~ fix(0.170958)   # $OMEGA BLOCK(1) SAME
    etaiov_mat_9 ~ fix(0.170958)   # $OMEGA BLOCK(1) SAME
    etaiov_mat_10 ~ fix(0.170958)  # $OMEGA BLOCK(1) SAME

    etaiov_d1_1 ~ 0.872344     # Datafile S1 $OMEGA 24 0.872344; sqrt = 0.9340, Table S3 'IOV D1 (CV)' 0.934
    etaiov_d1_2 ~ fix(0.872344)    # $OMEGA BLOCK(1) SAME
    etaiov_d1_3 ~ fix(0.872344)    # $OMEGA BLOCK(1) SAME
    etaiov_d1_4 ~ fix(0.872344)    # $OMEGA BLOCK(1) SAME
    etaiov_d1_5 ~ fix(0.872344)    # $OMEGA BLOCK(1) SAME
    etaiov_d1_6 ~ fix(0.872344)    # $OMEGA BLOCK(1) SAME
    etaiov_d1_7 ~ fix(0.872344)    # $OMEGA BLOCK(1) SAME
    etaiov_d1_8 ~ fix(0.872344)    # $OMEGA BLOCK(1) SAME
    etaiov_d1_9 ~ fix(0.872344)    # $OMEGA BLOCK(1) SAME
    etaiov_d1_10 ~ fix(0.872344)   # $OMEGA BLOCK(1) SAME

    # ==================================================================
    # RESIDUAL ERROR -- Table S3
    # Exponential (implemented by the authors as additive on the log scale),
    # with a separate magnitude per bioanalytical method x study phase x
    # time-since-last-dose window. 'Late' is more than 5 h since the last
    # dose (the disposition phase); 'Early' is 5 h or less (the absorption
    # phase).
    # ==================================================================
    expSdP1NoneLate <- 0.0795
    label("Residual SD, phase I, no peak separation or formic acid, >5 h since last dose (log scale)")   # Table S3 'RUV Phase I TSLD>5h (CV)' 0.0795 [RSE 3.29%]
    expSdP1NoneEarly <- 0.485
    label("Residual SD, phase I, no peak separation or formic acid, <=5 h since last dose (log scale)")  # Table S3 'RUV Phase I TSLD<=5h (CV)' 0.485 [RSE 2.53%]
    expSdP2Late <- 0.260
    label("Residual SD, phase II, >5 h since last dose (log scale)")                                     # Table S3 'RUV Phase II TSLD>5h (CV)' 0.260 [RSE 9.77%]
    expSdP2Early <- 0.380
    label("Residual SD, phase II, <=5 h since last dose (log scale)")                                    # Table S3 'RUV Phase II TSLD<=5h (CV)' 0.380 [RSE 3.33%]
    expSdP3NoneLate <- 0.232
    label("Residual SD, phase III, no peak separation, >5 h since last dose (log scale)")                # Table S3 'RUV Phase III TSLD>5h (CV)' 0.232 [RSE 15.1%]
    expSdP3NoneEarly <- 0.408
    label("Residual SD, phase III, no peak separation, <=5 h since last dose (log scale)")               # Table S3 'RUV Phase III TSLD<=5h (CV)' 0.408 [RSE 4.86%]
    expSdP1SepfaLate <- 0.138
    label("Residual SD, phase I, peak separation plus formic acid, >5 h since last dose (log scale)")    # Table S3 'RUV Phase I TSLD>5h (CV)' 0.138 [RSE 6.30%], 'with separation of interfering peak and addition of formic acid' block
    expSdP1SepfaEarly <- 0.810
    label("Residual SD, phase I, peak separation plus formic acid, <=5 h since last dose (log scale)")   # Table S3 'RUV Phase I TSLD<=5h (CV)' 0.810 [RSE 5.41%], 'with separation of interfering peak and addition of formic acid' block
    expSdP1SepLate <- 0.117
    label("Residual SD, phase I, peak separation only, >5 h since last dose (log scale)")                # Table S3 'RUV Phase I TSLD>5h (CV)' 0.117 [RSE 5.21%], 'with separation of interfering peak only' block
    expSdP1SepEarly <- 0.467
    label("Residual SD, phase I, peak separation only, <=5 h since last dose (log scale)")               # Table S3 'RUV Phase I TSLD<=5h (CV)' 0.467 [RSE 5.25%], 'with separation of interfering peak only' block
    expSdP3SepLate <- 0.235
    label("Residual SD, phase III, peak separation only, >5 h since last dose (log scale)")              # Table S3 'RUV Phase III TSLD>5h (CV)' 0.235 [RSE 2.64%], 'with separation of interfering peak only' block
    expSdP3SepEarly <- 0.355
    label("Residual SD, phase III, peak separation only, <=5 h since last dose (log scale)")             # Table S3 'RUV Phase III TSLD<=5h (CV)' 0.355 [RSE 3.68%], 'with separation of interfering peak only' block

    # IIV on the residual error itself. Strata whose IIV Table S3 reports as
    # 0 (FIX) carry no eta here rather than a zero-variance one.
    etaruvP1NoneLate ~ 0.183184   # Table S3 'IIV RUV Phase I TSLD>5h (CV)' 0.428; 0.428^2
    etaruvP1NoneEarly ~ 0.062500  # Table S3 'IIV RUV Phase I TSLD<=5h (CV)' 0.250; 0.250^2
    etaruvP2Late ~ 0.096721       # Table S3 'IIV RUV Phase II TSLD>5h (CV)' 0.311; 0.311^2
    etaruvP3NoneLate ~ 0.142884   # Table S3 'IIV RUV Phase III TSLD>5h (CV)' 0.378; 0.378^2
    etaruvP1SepfaLate ~ 0.067081  # Table S3 'IIV RUV Phase I TSLD>5h (CV)' 0.259; 0.259^2, formic-acid block
    etaruvP1SepLate ~ 0.160000    # Table S3 'IIV RUV Phase I TSLD>5h (CV)' 0.400; 0.400^2, peak-separation-only block
    etaruvP1SepEarly ~ 0.148996   # Table S3 'IIV RUV Phase I TSLD<=5h (CV)' 0.386; 0.386^2, peak-separation-only block
  })

  model({
    # ------------------------------------------------------------------
    # 1. Formulation selection. All six formulation indicators equal to 0
    #    selects formulation 5, the reference (Table S3).
    # ------------------------------------------------------------------
    isF5 <- 1 - FORM_ELA_F1 - FORM_ELA_F2 - FORM_ELA_F3 - FORM_ELA_F4 - FORM_ELA_F6

    lfdepotSel <-
      lfdepot_f1 * FORM_ELA_F1 + lfdepot_f2 * FORM_ELA_F2 +
      lfdepot_f3 * FORM_ELA_F3 + lfdepot_f4 * FORM_ELA_F4 +
      lfdepot_f5 * isF5 + lfdepot_f6 * FORM_ELA_F6
    etalfdepotSel <-
      etalfdepot_f1 * FORM_ELA_F1 + etalfdepot_f2 * FORM_ELA_F2 +
      etalfdepot_f3 * FORM_ELA_F3 + etalfdepot_f4 * FORM_ELA_F4 +
      etalfdepot_f5 * isF5 + etalfdepot_f6 * FORM_ELA_F6

    ld1Sel <-
      ld1_f1 * FORM_ELA_F1 + ld1_f2 * FORM_ELA_F2 + ld1_f3 * FORM_ELA_F3 +
      ld1_f4 * FORM_ELA_F4 + ld1_f5 * isF5 + ld1_f6 * FORM_ELA_F6
    # IIV on D1 is 0 (FIX) for formulations 1-3.
    etald1Sel <-
      etald1_f4 * FORM_ELA_F4 + etald1_f5 * isF5 + etald1_f6 * FORM_ELA_F6

    lmatSel <-
      lmat_f1 * FORM_ELA_F1 + lmat_f2 * FORM_ELA_F2 + lmat_f3 * FORM_ELA_F3 +
      lmat_f4 * FORM_ELA_F4 + lmat_f5 * isF5 + lmat_f6 * FORM_ELA_F6
    # IIV on MAT is 0 (FIX) for formulations 2-6.
    etalmatSel <- etalmat_f1 * FORM_ELA_F1

    ltlagSel <-
      ltlag_f1 * FORM_ELA_F1 + ltlag_f2 * FORM_ELA_F2 + ltlag_f3 * FORM_ELA_F3 +
      ltlag_f4 * FORM_ELA_F4 + ltlag_f5 * isF5 + ltlag_f6 * FORM_ELA_F6

    # ------------------------------------------------------------------
    # 2. Inter-occasion variability. Ten occasion slots each on MAT and D1,
    #    all sharing one variance (Datafile S1 $OMEGA BLOCK(1) SAME).
    # ------------------------------------------------------------------
    iovMat <-
      etaiov_mat_1 * (OCC == 1) + etaiov_mat_2 * (OCC == 2) +
      etaiov_mat_3 * (OCC == 3) + etaiov_mat_4 * (OCC == 4) +
      etaiov_mat_5 * (OCC == 5) + etaiov_mat_6 * (OCC == 6) +
      etaiov_mat_7 * (OCC == 7) + etaiov_mat_8 * (OCC == 8) +
      etaiov_mat_9 * (OCC == 9) + etaiov_mat_10 * (OCC == 10)
    iovD1 <-
      etaiov_d1_1 * (OCC == 1) + etaiov_d1_2 * (OCC == 2) +
      etaiov_d1_3 * (OCC == 3) + etaiov_d1_4 * (OCC == 4) +
      etaiov_d1_5 * (OCC == 5) + etaiov_d1_6 * (OCC == 6) +
      etaiov_d1_7 * (OCC == 7) + etaiov_d1_8 * (OCC == 8) +
      etaiov_d1_9 * (OCC == 9) + etaiov_d1_10 * (OCC == 10)

    # ------------------------------------------------------------------
    # 3. Individual absorption parameters
    #    Datafile S1 $PK: TVFREL is multiplied by the age effect and by the
    #    sigmoidal dose effect before the eta is applied; TVD1 by the age
    #    effect; TVMAT by the food effect.
    # ------------------------------------------------------------------
    doseFdepot <- 1 + emax_dose_fdepot * DOSE_ELA_MG^gamma_dose_fdepot /
      (ed50_dose_fdepot^gamma_dose_fdepot + DOSE_ELA_MG^gamma_dose_fdepot)

    fdepot <- exp(lfdepotSel + etalfdepotSel) * doseFdepot *
      exp(e_age_fdepot * (AGE - 44))
    d1 <- exp(ld1Sel + etald1Sel + iovD1) * exp(e_age_d1 * (AGE - 44))
    mat <- exp(lmatSel + etalmatSel + iovMat) * (1 + e_fed_mat * FED)
    ka <- 1 / mat
    # The lag time is fixed to 0 in the ELATIVE phase III study.
    tlag <- exp(ltlagSel) * (1 - STUDY_GFT505B_319_1)

    # ------------------------------------------------------------------
    # 4. Individual disposition parameters
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl_q *
      exp(e_alb_cl * (ALB - 43.30))
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc_vp
    q <- exp(lq + etalq) * (WT / 75)^e_wt_cl_q *
      exp(e_alb_q * (ALB - 43.30)) * (1 + e_sexf_q * SEXF)
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vc_vp *
      exp(e_age_vp * (AGE - 44)) * exp(e_bmi_vp * (BMI - 26.44)) *
      (1 + e_pbc_vp * DIS_PBC)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ------------------------------------------------------------------
    # 5. ODE system -- NONMEM ADVAN4 TRANS4 with a zero-order input into the
    #    depot (Datafile S1 $SUBROUTINE / $PK).
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central +
      k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    dur(depot) <- d1
    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # 6. Observation and residual error
    #    Datafile S1 $ERROR shifts the log-scale individual prediction by the
    #    bioanalytical-method coefficient, i.e. multiplies the predicted
    #    concentration by exp(coefficient).
    # ------------------------------------------------------------------
    Cc <- (central / vc) *
      exp(e_assay_sepip_cc * ASSAY_SEPIP + e_assay_sepip_fa_cc * ASSAY_SEPIP_FA)

    # Residual-error stratum selection. Phase II carried only the
    # unmodified bioanalytical method, and phase III only the unmodified
    # method or peak separation, so the source estimates no formic-acid
    # stratum outside phase I.
    tsld <- tad()
    isEarly <- (tsld <= 5)
    isLate <- 1 - isEarly
    isPh2 <- STUDY_PHASE2
    isPh3 <- STUDY_PHASE3
    isPh1 <- 1 - isPh2 - isPh3
    isSep <- ASSAY_SEPIP
    isSepfa <- ASSAY_SEPIP_FA
    isNone <- 1 - isSep - isSepfa

    expSdSel <-
      isPh1 * (
        isNone * (isEarly * expSdP1NoneEarly + isLate * expSdP1NoneLate) +
          isSep * (isEarly * expSdP1SepEarly + isLate * expSdP1SepLate) +
          isSepfa * (isEarly * expSdP1SepfaEarly + isLate * expSdP1SepfaLate)
      ) +
      isPh2 * (isEarly * expSdP2Early + isLate * expSdP2Late) +
      isPh3 * (
        (1 - isSep) * (isEarly * expSdP3NoneEarly + isLate * expSdP3NoneLate) +
          isSep * (isEarly * expSdP3SepEarly + isLate * expSdP3SepLate)
      )

    etaruvSel <-
      isPh1 * (
        isNone * (isEarly * etaruvP1NoneEarly + isLate * etaruvP1NoneLate) +
          isSep * (isEarly * etaruvP1SepEarly + isLate * etaruvP1SepLate) +
          isSepfa * isLate * etaruvP1SepfaLate
      ) +
      isPh2 * isLate * etaruvP2Late +
      isPh3 * (1 - isSep) * isLate * etaruvP3NoneLate

    expSdCc <- expSdSel * exp(etaruvSel)

    Cc ~ lnorm(expSdCc)
  })
}
