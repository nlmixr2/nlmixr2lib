Salinger_2019_pretomanid <- function() {
  description <- "One-compartment population PK model with three-transit-compartment absorption for pretomanid in healthy subjects and subjects with drug-sensitive, multidrug-resistant or extensively drug-resistant pulmonary tuberculosis"
  reference <- "Salinger DH, Subramoney V, Everitt D, Nedelman JR. Population pharmacokinetics of the antituberculosis agent pretomanid. Antimicrob Agents Chemother. 2019;63(10):e00907-19. doi:10.1128/AAC.00907-19"
  vignette <- "Salinger_2019_pretomanid"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  # Oral doses enter the FIRST TRANSIT compartment, not `depot`: `depot` here is
  # the source control stream's ABS compartment, which sits downstream of the
  # three-compartment transit chain and empties into `central` at KA. Declared
  # explicitly because the registry's fallback heuristic would otherwise report
  # `depot` and send users to the wrong compartment.
  dosing <- "transit1"

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling on CL (exponent 0.75) and V2 (exponent 1), both normalised to the 55 kg reference subject of Salinger 2019 Results 'Model application'. Baseline value; Winsorized at median +/- 5 SD as described in Materials and Methods.",
      source_name = "WT"
    ),
    FED = list(
      description = "Fed-versus-fasted state at the dose record",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (fed) is the model reference; F1 = 1 for 200 mg administered fed",
      notes = "The source control stream (Table S3) carries the OPPOSITE polarity column FASTED, so every effect is applied here as e_..._..^(1 - FED). Fed/fasted was not well controlled in studies NC-002, NC-003, NC-005, NC-006 and Nix-TB; per Table S2 footnote a those subjects were assumed FED and separate study indicators absorb the residual food effect.",
      source_name = "FASTED (= 1 - FED)"
    ),
    DOSE_PRETOMANID_MG = list(
      description = "Administered pretomanid oral dose level",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-function effects on F1 (fasted only), KA, MTT (fasted only) and V2, all normalised to the 200 mg reference dose. Dose range 50-1500 mg. A separate indicator (e_dose_1000mg_fdepot) captures the 1000 mg FED records only.",
      source_name = "DOSE"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Apparent clearance in females was 18% lower than in males (Salinger 2019 Results 'Final model').",
      source_name = "FEMALE"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-subject (non-tuberculosis) cohort indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (subject with tuberculosis)",
      notes = "The source control stream derives this as HV = 1 - DS - MDR, i.e. a subject who is neither drug-sensitive-TB nor MDR-or-worse TB. All 162-211 healthy subjects were from North America and none was HIV-positive.",
      source_name = "HV"
    ),
    DIS_TB_MDR = list(
      description = "Rifampicin-and-isoniazid-resistant (MDR or worse) tuberculosis indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (drug-sensitive tuberculosis or healthy subject)",
      notes = "1 for MDR-TB, treatment-intolerant MDR-TB, non-responsive MDR-TB and XDR-TB. Reproduces the source control-stream expression (MDR - NIX + MDRIT + MDRNR + MDRXDR) used on CL, which evaluates to 1 for every drug-resistant subject whether or not they were enrolled in Nix-TB. On V2 the source excludes XDR, which is recovered here as DIS_TB_MDR * (1 - DIS_TB_XDR_STRICT).",
      source_name = "MDR / NIXPT"
    ),
    DIS_TB_XDR_STRICT = list(
      description = "Extensively-drug-resistant (XDR) tuberculosis stratum indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (every other resistance stratum, including MDR and treatment-intolerant / non-responsive MDR)",
      notes = "Derived in the source control stream as MDRXDR = 1 when NIXPT == 3. Carries its own V2 effect (1.75) that replaces, rather than multiplies, the MDR V2 effect. The strict sibling is correct here because Salinger 2019 contrasts XDR against MDR and TI/NR-MDR separately rather than pooling them.",
      source_name = "NIXPT == 3"
    ),
    HIV_POS = list(
      description = "HIV-positive comorbidity indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (HIV-negative)",
      notes = "Acts on BOTH CL (0.842) and F1 (0.789); the net effect on apparent oral clearance CL/F1 is 0.842/0.789 = 1.067, which is the '6% higher' apparent clearance reported in Results 'Final model'. No healthy subject was HIV-positive.",
      source_name = "HIV"
    ),
    CONMED_EFV = list(
      description = "Concomitant efavirenz indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant efavirenz)",
      notes = "Time-varying: the source used time-dependent on/off histories for antiretroviral comedications (Materials and Methods 'Data'). Acts on CL (2.17) and F1 (1.24); together with the HIV_POS effects this reproduces the reported 46% exposure reduction.",
      source_name = "EFV"
    ),
    CONMED_LPV = list(
      description = "Concomitant ritonavir-boosted lopinavir indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant lopinavir/ritonavir)",
      notes = "Time-varying antiretroviral comedication. Acts on CL only (1.14); together with the HIV_POS effects this reproduces the reported 17% exposure reduction.",
      source_name = "LPVR"
    ),
    CONMED_CYP3A4_IND = list(
      description = "Concomitant CYP3A4-inducing antiretroviral indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant CYP3A4-inducing antiretroviral)",
      notes = "Pools the less-prevalent CYP3A4-inducing antiretrovirals other than efavirenz, which the source grouped as INDUC (Materials and Methods 'Pre-Nix model'). Inducing strengths are not enumerated by the source. CYP3A4 inhibitors were also planned as a group but no subject was on an alternative inhibitor, so no companion indicator exists.",
      source_name = "INDUC"
    ),
    CONMED_MOXIFLOXACIN = list(
      description = "Concomitant moxifloxacin indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant moxifloxacin)",
      notes = "Regimen partner. Carries a main effect on CL (0.967) plus moxifloxacin-pyrazinamide and bedaquiline-moxifloxacin-pyrazinamide interaction effects on CL and F1; the interactions are products of the indicators, exactly as the source control stream writes them.",
      source_name = "MOX"
    ),
    CONMED_PYRAZINAMIDE = list(
      description = "Concomitant pyrazinamide indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant pyrazinamide)",
      notes = "Regimen partner. Enters only through the moxifloxacin-pyrazinamide and bedaquiline-moxifloxacin-pyrazinamide interaction products; the source retained no pyrazinamide-alone effect.",
      source_name = "PZA"
    ),
    CONMED_BEDAQUILINE = list(
      description = "Concomitant bedaquiline indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant bedaquiline)",
      notes = "Regimen partner. Enters only through the bedaquiline-moxifloxacin-pyrazinamide (BPaMZ) interaction product on CL and F1; bedaquiline given without moxifloxacin and pyrazinamide (as in the Nix-TB BPaL regimen) therefore has no effect, which is the source's finding that bedaquiline and linezolid together had little impact on pretomanid exposure.",
      source_name = "BDQ"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-function effect on CL normalised to the 35 g/L reference. Baseline (or screening) value; Winsorized at median +/- 5 SD.",
      source_name = "ALB"
    ),
    TBILI = list(
      description = "Baseline total bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-function effect on F1 normalised to the 5 umol/L reference. Reported in SI units by the source (Table S2), so no unit conversion is applied. Baseline (or screening) value; Winsorized at median +/- 5 SD.",
      source_name = "TBIL"
    ),
    STUDY_NC003 = list(
      description = "Study NC-003 cohort indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any other study)",
      notes = "Enters only as a scaler on the magnitude of the MTT random effect, not on any typical value, because fed/fasted conditions in NC-003 were uncertain and the study needed its own absorption variability.",
      source_name = "NC3"
    ),
    STUDY_NC005 = list(
      description = "Study NC-005 cohort indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any other study)",
      notes = "Study-specific adjustment to the absorption model (KA and MTT) because fed/fasted conditions in NC-005 were uncertain. The late median Tmax of 6.5 h in this study (Results 'Tmax') is reproduced by the 0.186-fold KA alone; see the vignette Errata for the MTT boundary estimate.",
      source_name = "NC5"
    ),
    STUDY_NIXTB = list(
      description = "Nix-TB study cohort indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any other study)",
      notes = "Nix-TB was added in the final modelling stage and required its own F1 typical value (1.54) and F1 random-effect magnitude, its own Box-Cox shape parameters for the CL and V2 random effects, and a second step-up in CL from week 6 onward.",
      source_name = "NIX"
    ),
    OOC1 = list(
      description = "Interoccasion-variability occasion-1 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not occasion 1)",
      notes = "Up to three mutually exclusive occasions were defined per trial from trial length and sampling schedule (Table S1); e.g. for Nix-TB occasions 1, 2 and 3 are weeks 2, 8 and 16. Interoccasion random effects act on F1 and CL and are correlated within an occasion.",
      source_name = "OCC1"
    ),
    OOC2 = list(
      description = "Interoccasion-variability occasion-2 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not occasion 2)",
      notes = "See OOC1. Occasion definitions are per-study and listed in Table S1.",
      source_name = "OCC2"
    ),
    OOC3 = list(
      description = "Interoccasion-variability occasion-3 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not occasion 3)",
      notes = "See OOC1. Only trials long enough to support three occasions carry this indicator.",
      source_name = "OCC3"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on CL and F1 (Materials and Methods 'Pre-Nix model') but not retained as significant in the final model (Results 'Final model')."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m2",
      type = "continuous",
      notes = "Tested but not retained; body weight was retained instead."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Tested but not retained. More than 99% of subjects had normal or only mildly impaired renal function by eGFR, and in two mass-balance studies < 1% of parent drug was excreted in urine."
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate",
      units = "mL/min/1.73m2",
      type = "continuous",
      notes = "Tested but not retained; see CRCL."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested on CL and F1 but not retained; total bilirubin and albumin were the hepatic markers retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested on CL and F1 but not retained; see ALT."
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Race was tested but not retained as significant. Source columns CAU and OTR encoded Caucasian and Other against a Black reference; 58.4% of the pooled cohort was Black (Table S2)."
    ),
    CONMED_CLOFAZIMINE = list(
      description = "Concomitant clofazimine indicator",
      units = "(binary)",
      type = "binary",
      notes = "Present in the source analysis data set (control-stream column CFZ) and assessed on clearance and bioavailability, but no clofazimine effect was retained in the final model."
    ),
    STUDY_NC002 = list(
      description = "Study NC-002 cohort indicator",
      units = "(binary)",
      type = "binary",
      notes = "Present in the source analysis data set (control-stream column NC2) and tested as a study-specific food-effect adjustment, but no NC-002 effect was retained in the final model."
    ),
    STUDY_NC006 = list(
      description = "Study NC-006 cohort indicator",
      units = "(binary)",
      type = "binary",
      notes = "Present in the source analysis data set (control-stream column NC6) and tested as a study-specific food-effect adjustment, but no NC-006 effect was retained in the final model."
    )
  )

  compartmentData <- list(
    transit1 = list(analyte = "pretomanid", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "pretomanid", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "pretomanid", units = "mg", specimen = "administration site", verified = TRUE),
    depot = list(analyte = "pretomanid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pretomanid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1054L,
    n_studies = 14L,
    n_observations = 17725L,
    age_range = "18-77 years",
    age_median = "27 years (healthy subjects), 31 years (subjects with tuberculosis)",
    weight_range = "29-121 kg",
    weight_median = "75 kg (healthy subjects), 53 kg (subjects with tuberculosis)",
    sex_female_pct = 35,
    race_ethnicity = c(Black = 58.4, Caucasian = 14.8, Other = 26.8),
    disease_state = "healthy subjects and subjects with drug-sensitive, multidrug-resistant, treatment-intolerant or non-responsive multidrug-resistant, or extensively drug-resistant pulmonary tuberculosis; 24% of subjects with tuberculosis were HIV-positive",
    dose_range = "50-1500 mg single oral dose; 100-1000 mg once daily for 7 days to 6 months",
    regions = "North America (all healthy subjects); sub-Saharan Africa (> 95% of subjects with tuberculosis)",
    co_medication = "pretomanid given as monotherapy or with bedaquiline, clofazimine, moxifloxacin, linezolid and/or pyrazinamide; antiretrovirals including efavirenz and lopinavir/ritonavir in HIV-positive subjects",
    notes = "Six phase 1, six phase 2 and two phase 3 studies (Salinger 2019 Results 'Data'; per-study detail in Table S1, covariate summary in Table S2). Table S2 tabulates 1056 subjects against the 1054 quoted in Results; the Results figure is used here. Renal function: more than 99% of subjects had normal or only mildly impaired renal function by eGFR."
  )

  ini({
    # ---- Absorption ------------------------------------------------------
    # Three transit compartments (rate KTR = 3/MTT) feed a first-order
    # absorption compartment (rate KA); Table S3 $MODEL / Figure S2.
    lka <- log(1.38)
    label("First-order absorption rate constant from the absorption compartment for a 200 mg fed dose (1/h)") # Table 1 'KA, h-1' = 1.38 (95% CI 1.25 to 1.52)
    lmtt <- log(1.25)
    label("Mean transit time through the three transit compartments for a 200 mg fed dose (h)") # Table 1 'MTT, h' = 1.25 (95% CI 1.14 to 1.35)

    # ---- Disposition -----------------------------------------------------
    # Piecewise-constant clearance with two step-ups. Salinger 2019
    # Materials and Methods 'Pre-Nix model' describes a 'two-step CL (with
    # initial and steady-state CL)'; 'Final model' adds a further
    # 'step-up in CL at week 6' for Nix-TB subjects only. Table S3 writes
    # this additively as LOG(THETA(2) + SSD*THETA(11) + SSWK6*THETA(38)),
    # so the two late arms below are sums of reported estimates.
    lcl <- log(3.30)
    label("Apparent oral clearance before the first clearance step, for a 55 kg reference subject (L/h)") # Table 1 'CL, liters/h' = 3.30 (95% CI 3.14 to 3.46)
    lcl_late <- log(3.30 + 0.175)
    label("Apparent oral clearance from study day 5 onward, for a 55 kg reference subject (L/h)") # Table 1 'CL' 3.30 + 'SS CL, liters/h' 0.175 (95% CI 0.135 to 0.214); Table S3 THETA(2) + THETA(11)
    ltclchange <- fixed(log(96))
    label("Time after the first dose at which clearance steps to its day-5 value (h)") # Table S3: SSD = 1 when STDAY >= 5, i.e. 4 days = 96 h after the first dose; a structural breakpoint, not estimated
    lcl_late2 <- log(3.30 + 0.175 + 0.466)
    label("Apparent oral clearance from week 6 onward in Nix-TB subjects, for a 55 kg reference subject (L/h)") # Table 1 'CL' 3.30 + 'SS CL' 0.175 + 'CL NIX WK >= 6, liters/h' 0.466 (95% CI 0.285 to 0.646); Table S3 THETA(2) + THETA(11) + THETA(38)
    ltclchange2 <- fixed(log(1008))
    label("Time after the first dose at which clearance steps again in Nix-TB subjects (h)") # Table S3: SSWK6 = 1 when STUDYNO == 201 and STWEEK >= 6, i.e. 6 weeks = 1008 h; a structural breakpoint, not estimated
    lvc <- log(90.4)
    label("Apparent central volume of distribution for a 55 kg reference subject at 200 mg (L)") # Table 1 'V2, liters' = 90.4 (95% CI 85.8 to 94.9)

    # ---- Relative bioavailability ----------------------------------------
    lfdepot <- fixed(log(1))
    label("Relative bioavailability for 200 mg administered fed (unitless)") # Table 1 'F1 (fixed)' = 1; the reference condition against which every other F1 effect is expressed

    # ---- Allometric exponents (held at the canonical values) -------------
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on (WT/55) for CL (unitless)") # Table 1 'CL ~ WT (fixed)' = 0.75; reported without a confidence interval
    e_wt_vc <- fixed(1)
    label("Allometric exponent on (WT/55) for V2 (unitless)") # Table 1 'V2 ~ WT (fixed)' = 1; reported without a confidence interval

    # ---- Covariate effects on relative bioavailability -------------------
    e_fed_fdepot <- 0.513
    label("Fold change in F1 when dosed fasted rather than fed (unitless)") # Table 1 'F1 ~ FASTED' = 0.513 (95% CI 0.485 to 0.541)
    e_dose_fdepot <- -0.264
    label("Power exponent on (DOSE/200) for F1 in fasted subjects (unitless)") # Table 1 'F1 ~ DOSE&FASTED' = -0.264 (95% CI -0.302 to -0.226)
    e_dose_1000mg_fdepot <- -0.00302
    label("Power exponent on (DOSE/200) for F1 in fed subjects dosed 1000 mg (unitless)") # Table 1 'F1 ~ FED&1,000 mg' = -0.00302 (95% CI -0.148 to 0.142)
    e_moxifloxacin_pyrazinamide_fdepot <- 0.925
    label("Fold change in F1 on concomitant moxifloxacin plus pyrazinamide (unitless)") # Table 1 'F1 ~ MOX * PZA' = 0.925 (95% CI 0.848 to 1.00)
    e_bedaquiline_moxifloxacin_pyrazinamide_fdepot <- 1.25
    label("Additional fold change in F1 on concomitant bedaquiline plus moxifloxacin plus pyrazinamide (unitless)") # Table 1 'F1 ~ BDQ * MOX * PZA' = 1.25 (95% CI 1.04 to 1.46)
    e_efv_fdepot <- 1.24
    label("Fold change in F1 on concomitant efavirenz (unitless)") # Table 1 'F1 ~ EFV' = 1.24 (95% CI 1.03 to 1.45)
    e_hiv_pos_fdepot <- 0.789
    label("Fold change in F1 in HIV-positive subjects (unitless)") # Table 1 'F1 ~ HIV' = 0.789 (95% CI 0.721 to 0.856)
    e_tbili_fdepot <- 0.0880
    label("Power exponent on (TBILI/5) for F1 (unitless)") # Table 1 'F1 ~ TBIL (ref = 5)' = 0.0880 (95% CI 0.0495 to 0.126)
    e_study_nixtb_fdepot <- 1.54
    label("Fold change in F1 in the Nix-TB study (unitless)") # Table 1 'F1 ~ NIX' = 1.54 (95% CI 1.28 to 1.80)

    # ---- Covariate effects on absorption ---------------------------------
    e_fed_ka <- 0.482
    label("Fold change in KA when dosed fasted rather than fed (unitless)") # Table 1 'KA ~ FASTED' = 0.482 (95% CI 0.453 to 0.511)
    e_dose_ka <- -0.128
    label("Power exponent on (DOSE/200) for KA (unitless)") # Table 1 'KA ~ DOSE' = -0.128 (95% CI -0.158 to -0.0990)
    e_study_nc005_ka <- 0.186
    label("Fold change in KA in study NC-005 (unitless)") # Table 1 'KA ~ NC5' = 0.186 (95% CI 0.150 to 0.221)
    e_fed_mtt <- 0.311
    label("Fold change in MTT when dosed fasted rather than fed (unitless)") # Table 1 'MTT ~ FASTED' = 0.311 (95% CI 0.293 to 0.329)
    e_dose_mtt <- -0.155
    label("Power exponent on (DOSE/200) for MTT in fasted subjects (unitless)") # Table 1 'MTT ~ DOSE&FASTED' = -0.155 (95% CI -0.187 to -0.123)
    e_study_nc005_mtt <- 6.95e-07
    label("Fold change in MTT in study NC-005 (unitless)") # Table 1 'MTT ~ NC5' = 6.95E-07; a boundary estimate that collapses MTT to ~1e-06 h in NC-005 (see vignette Errata). Confirmed against Table S4: NC-005 median MTT + 1/KA = 3.8 h equals 1/KA = 3.9 h on its own

    # ---- Covariate effects on clearance ----------------------------------
    e_dis_healthy_cl <- 1.16
    label("Fold change in CL in healthy subjects (unitless)") # Table 1 'CL ~ HS' = 1.16 (95% CI 1.09 to 1.23)
    e_moxifloxacin_cl <- 0.967
    label("Fold change in CL on concomitant moxifloxacin (unitless)") # Table 1 'CL ~ MOX' = 0.967 (95% CI 0.902 to 1.03)
    e_moxifloxacin_pyrazinamide_cl <- 0.733
    label("Additional fold change in CL on concomitant moxifloxacin plus pyrazinamide (unitless)") # Table 1 'CL ~ MOX * PZA' = 0.733 (95% CI 0.663 to 0.804)
    e_dis_tb_mdr_cl <- 1.15
    label("Fold change in CL in MDR, TI/NR MDR or XDR tuberculosis (unitless)") # Table 1 'CL ~ MDR, TI/NR MDR, or XDR' = 1.15 (95% CI 1.04 to 1.27)
    e_bedaquiline_moxifloxacin_pyrazinamide_cl <- 1.32
    label("Additional fold change in CL on concomitant bedaquiline plus moxifloxacin plus pyrazinamide (unitless)") # Table 1 'CL ~ BDQ * MOX * PZA' = 1.32 (95% CI 1.12 to 1.51)
    e_efv_cl <- 2.17
    label("Fold change in CL on concomitant efavirenz (unitless)") # Table 1 'CL ~ EFV' = 2.17 (95% CI 1.89 to 2.45)
    e_lpv_cl <- 1.14
    label("Fold change in CL on concomitant lopinavir/ritonavir (unitless)") # Table 1 'CL ~ LPVR' = 1.14 (95% CI 1.10 to 1.18)
    e_hiv_pos_cl <- 0.842
    label("Fold change in CL in HIV-positive subjects (unitless)") # Table 1 'CL ~ HIV' = 0.842 (95% CI 0.779 to 0.905)
    e_cyp3a4_ind_cl <- 1.35
    label("Fold change in CL on a concomitant CYP3A4-inducing antiretroviral (unitless)") # Table 1 'CL ~ INDUC' = 1.35 (95% CI 1.24 to 1.46)
    e_sexf_cl <- 0.837
    label("Fold change in CL in females (unitless)") # Table 1 'CL ~ FEMALE' = 0.837 (95% CI 0.808 to 0.867)
    e_alb_cl <- 0.200
    label("Power exponent on (ALB/35) for CL (unitless)") # Table 1 'CL ~ ALB (ref = 35)' = 0.200 (95% CI 0.0789 to 0.322)

    # ---- Covariate effects on volume -------------------------------------
    e_dose_vc <- 0.111
    label("Power exponent on (DOSE/200) for V2 (unitless)") # Table 1 'V2 ~ DOSE' = 0.111 (95% CI 0.0845 to 0.137)
    e_dis_tb_mdr_vc <- 1.44
    label("Fold change in V2 in MDR or TI/NR MDR tuberculosis (unitless)") # Table 1 'V2 ~ MDR or TI/NR MDR' = 1.44 (95% CI 1.24 to 1.64); excludes XDR, which carries its own effect
    e_dis_tb_xdr_strict_vc <- 1.75
    label("Fold change in V2 in XDR tuberculosis (unitless)") # Table 1 'V2 ~ XDR' = 1.75 (95% CI 1.39 to 2.11)

    # ---- Random-effect shape and magnitude modifiers ----------------------
    e_study_nixtb_etalfdepot <- 0.919
    label("Log-scale scaler on the magnitude of the F1 random effect in the Nix-TB study (unitless)") # Table 1 'F1 Var ~ NIX' = 0.919 (95% CI 0.558 to 1.28); Table S3 F1V = EXP(THETA(37)*NIX) multiplies ETA(6)
    e_study_nc003_etalmtt <- -0.645
    label("Log-scale scaler on the magnitude of the MTT random effect in study NC-003 (unitless)") # Table 1 'MTT Var ~ NC3' = -0.645 (95% CI -1.06 to -0.226); Table S3 MTTV3 = EXP(THETA(21)*NC3) multiplies ETA(5)
    boxcox_lvc_nonnix <- 9.55
    label("Box-Cox shape parameter for the V2 random effect outside Nix-TB (unitless)") # Table 1 'Box-Cox V2 non-NIX' = 9.55 (95% CI 2.25 to 16.9)
    boxcox_lvc_nixtb <- 26.0
    label("Box-Cox shape parameter for the V2 random effect in Nix-TB (unitless)") # Table 1 'Box-Cox V2 NIX' = 26.0 (95% CI -8.54 to 60.5)
    boxcox_lcl_nonnix <- 1.36
    label("Box-Cox shape parameter for the CL random effect outside Nix-TB (unitless)") # Table 1 'Box-Cox CL non-NIX' = 1.36 (95% CI 0.733 to 1.99)
    boxcox_lcl_nixtb <- 2.78
    label("Box-Cox shape parameter for the CL random effect in Nix-TB (unitless)") # Table 1 'Box-Cox CL NIX' = 2.78 (95% CI 1.36 to 4.20)

    # ---- Interindividual variability -------------------------------------
    # Table 1 'OMEGA matrix terms'. Variances of normally distributed etas;
    # exp(eta) multiplies each parameter, except for CL and V2 where a
    # Box-Cox transform is applied to the normal eta first (Table S3).
    etae_dose_fdepot ~ 0.0274 # Table 1 OMEGA 'F1 ~ dose/fasted' = 0.0274 (95% CI 0.018 to 0.0368); Table S3 ETA(1), which enters the fasted F1 dose slope
    etalcl ~ 0.0373 # Table 1 OMEGA 'CL' = 0.0373 (95% CI 0.029 to 0.0455); Table S3 ETA(2), Box-Cox transformed
    etalvc ~ 0.00892 # Table 1 OMEGA 'V2' = 0.00892 (95% CI 0.0018 to 0.016); Table S3 ETA(3), Box-Cox transformed
    # Correlated block on KA, MTT and F1; Table S3 $OMEGA BLOCK(3) over
    # ETA(4), ETA(5), ETA(6) in that order.
    etalka + etalmtt + etalfdepot ~ c(
      0.309,
      0.0649, 0.593,
      -0.0339, 0.00391, 0.0227
    ) # Table 1 OMEGA rows 'KA' 0.309, 'KA-MTT' 0.0649, 'MTT' 0.593, 'KA-F1' -0.0339, 'MTT-F1' 0.00391, 'F1' 0.0227

    # ---- Interoccasion variability ---------------------------------------
    # Table S3 $OMEGA BLOCK(2) followed by two BLOCK(2) SAME blocks: one
    # correlated F1/CL pair drawn independently on each of up to three
    # occasions, with the same covariance matrix on every occasion. The
    # second and third occasions are fixed to the first occasion's values
    # to reproduce SAME.
    etaiov_lfdepot_1 + etaiov_lcl_1 ~ c(
      0.0412,
      0.0101, 0.0185
    ) # Table 1 OMEGA rows 'IOC F1' 0.0412 (95% CI 0.0362 to 0.0462), 'IOC F1-CL' 0.0101 (95% CI 0.00697 to 0.0132), 'IOC CL' 0.0185 (95% CI 0.0155 to 0.0214)
    etaiov_lfdepot_2 + etaiov_lcl_2 ~ c(
      fixed(0.0412),
      fixed(0.0101), fixed(0.0185)
    ) # Table S3 '$OMEGA BLOCK(2) SAME': occasion 2 shares the occasion-1 covariance matrix
    etaiov_lfdepot_3 + etaiov_lcl_3 ~ c(
      fixed(0.0412),
      fixed(0.0101), fixed(0.0185)
    ) # Table S3 '$OMEGA BLOCK(2) SAME': occasion 3 shares the occasion-1 covariance matrix

    # ---- Residual error --------------------------------------------------
    # Table S3 $ERROR: Y = IPRED + EPS(1)*(IPRED + 1e-04)^THETA(10) + EPS(2)
    # with a diagonal $SIGMA, i.e. independent power and additive terms
    # whose variances add (nlmixr2 combined2). Concentrations are in ng/mL.
    propSd <- 0.740
    label("Power-error coefficient on the predicted concentration (ng/mL^(1-powExp))") # Table 1 'Proportional error (variance)' = 0.548 (95% CI 0.506 to 0.589); SD = sqrt(0.548) = 0.740
    powExp <- 0.795
    label("Exponent of the predicted concentration in the power residual-error term (unitless)") # Table 1 'Error power' = 0.795 (95% CI 0.789 to 0.800); Table S3 THETA(10)
    addSd <- 3.39
    label("Additive residual standard deviation (ng/mL)") # Table 1 'Additive error (variance)' = 11.5 (95% CI 5.84 to 17.2); SD = sqrt(11.5) = 3.39
  })

  model({
    # ---- 1. Derived covariate and time terms -----------------------------
    # The source control stream carries FASTED; FED is the canonical column,
    # so every fasted effect is raised to the (1 - FED) power.
    fasted <- 1 - FED
    doseRatio <- DOSE_PRETOMANID_MG / 200

    # Table S3: DOSEIX = 1 only for the 1000 mg FED records.
    doseix1000fed <- FED * (DOSE_PRETOMANID_MG > 999.5) * (DOSE_PRETOMANID_MG < 1000.5)

    # Table S3: the V2 MDR effect uses (MDR - NIX + MDRIT + MDRNR), which is
    # 1 for MDR and TI/NR MDR subjects but 0 for XDR subjects, who instead
    # carry their own V2 effect.
    mdrNotXdr <- DIS_TB_MDR * (1 - DIS_TB_XDR_STRICT)

    # Piecewise-constant clearance breakpoints, measured from the first dose.
    # tafd() is hoisted onto its own line: rxode2's mu-reference walker
    # indexes the second element of every sub-call, so a zero-argument call
    # sharing a statement with a covariate or parameter aborts the parse.
    tafdNow <- tafd()
    tclchange <- exp(ltclchange)
    tclchange2 <- exp(ltclchange2)
    clStep1 <- (tafdNow >= tclchange)
    clStep2 <- STUDY_NIXTB * (tafdNow >= tclchange2)

    # ---- 2. Random-effect assembly ---------------------------------------
    # Interoccasion random effects: one independent correlated F1/CL draw per
    # occasion, selected by the mutually exclusive occasion indicators.
    iovlfdepot <- OOC1 * etaiov_lfdepot_1 + OOC2 * etaiov_lfdepot_2 + OOC3 * etaiov_lfdepot_3
    iovlcl <- OOC1 * etaiov_lcl_1 + OOC2 * etaiov_lcl_2 + OOC3 * etaiov_lcl_3

    # Box-Cox transformed CL and V2 random effects (Table S3), with separate
    # shape parameters inside and outside Nix-TB.
    boxcoxCl <- boxcox_lcl_nonnix * (1 - STUDY_NIXTB) + boxcox_lcl_nixtb * STUDY_NIXTB
    etaTrCl <- (exp(etalcl)^boxcoxCl - 1) / boxcoxCl
    boxcoxVc <- boxcox_lvc_nonnix * (1 - STUDY_NIXTB) + boxcox_lvc_nixtb * STUDY_NIXTB
    etaTrVc <- (exp(etalvc)^boxcoxVc - 1) / boxcoxVc

    # Study-specific scalers on the magnitude of the F1 and MTT random
    # effects (Table S3 F1V and MTTV3).
    sdScaleFdepot <- exp(e_study_nixtb_etalfdepot * STUDY_NIXTB)
    sdScaleMtt <- exp(e_study_nc003_etalmtt * STUDY_NC003)

    # ---- 3. Individual parameters ----------------------------------------
    fdepot <- exp(lfdepot + sdScaleFdepot * etalfdepot + iovlfdepot) *
      e_fed_fdepot^fasted *
      doseRatio^(e_dose_1000mg_fdepot * doseix1000fed +
        fasted * (e_dose_fdepot + etae_dose_fdepot)) *
      e_moxifloxacin_pyrazinamide_fdepot^(CONMED_MOXIFLOXACIN * CONMED_PYRAZINAMIDE) *
      e_bedaquiline_moxifloxacin_pyrazinamide_fdepot^(CONMED_BEDAQUILINE * CONMED_MOXIFLOXACIN * CONMED_PYRAZINAMIDE) *
      e_efv_fdepot^CONMED_EFV *
      e_hiv_pos_fdepot^HIV_POS *
      (TBILI / 5)^e_tbili_fdepot *
      e_study_nixtb_fdepot^STUDY_NIXTB

    ka <- exp(lka + etalka) *
      e_fed_ka^fasted *
      doseRatio^e_dose_ka *
      e_study_nc005_ka^STUDY_NC005

    mtt <- exp(lmtt + sdScaleMtt * etalmtt) *
      e_fed_mtt^fasted *
      doseRatio^(fasted * e_dose_mtt) *
      e_study_nc005_mtt^STUDY_NC005
    ktr <- 3 / mtt

    clBase <- exp(lcl) * (1 - clStep1) +
      exp(lcl_late) * clStep1 * (1 - clStep2) +
      exp(lcl_late2) * clStep2
    cl <- clBase * exp(etaTrCl + iovlcl) *
      (WT / 55)^e_wt_cl *
      e_dis_healthy_cl^DIS_HEALTHY *
      e_moxifloxacin_cl^CONMED_MOXIFLOXACIN *
      e_moxifloxacin_pyrazinamide_cl^(CONMED_MOXIFLOXACIN * CONMED_PYRAZINAMIDE) *
      e_dis_tb_mdr_cl^DIS_TB_MDR *
      e_bedaquiline_moxifloxacin_pyrazinamide_cl^(CONMED_BEDAQUILINE * CONMED_MOXIFLOXACIN * CONMED_PYRAZINAMIDE) *
      e_efv_cl^CONMED_EFV *
      e_lpv_cl^CONMED_LPV *
      e_hiv_pos_cl^HIV_POS *
      e_cyp3a4_ind_cl^CONMED_CYP3A4_IND *
      e_sexf_cl^SEXF *
      (ALB / 35)^e_alb_cl

    vc <- exp(lvc + etaTrVc) *
      (WT / 55)^e_wt_vc *
      doseRatio^e_dose_vc *
      e_dis_tb_mdr_vc^mdrNotXdr *
      e_dis_tb_xdr_strict_vc^DIS_TB_XDR_STRICT

    # ---- 4. Micro-constants ----------------------------------------------
    kel <- cl / vc

    # ---- 5. ODE system ---------------------------------------------------
    # Table S3 $MODEL: dose enters TRANS1; TRANS1 -> TRANS2 -> TRANS3 all at
    # KTR, TRANS3 -> ABS at KTR, ABS -> CENTRAL at KA, CENTRAL eliminated at
    # KEL. The three KTR-governed compartments are what MTT = 3/KTR refers
    # to, so the total mean absorption time is MTT + 1/KA.
    d/dt(transit1) <- -ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(depot) <- ktr * transit3 - ka * depot
    d/dt(central) <- ka * depot - kel * central

    # ---- 6. Relative bioavailability on the dosing compartment -----------
    f(transit1) <- fdepot

    # ---- 7. Observation and residual error -------------------------------
    # central is in mg and vc in L, so central/vc is mg/L; the source
    # dataset reports concentrations in ng/mL, which the residual-error
    # magnitudes are calibrated to.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + pow(propSd, powExp) + combined2()
  })
}
