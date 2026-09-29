Lu_2016_meropenem <- function() {
  description <- "Three-compartment population PK model for meropenem in the plasma and cerebrospinal fluid of adults with post-neurosurgical bacterial meningitis (Lu 2016): two-compartment plasma disposition with linear elimination from the central compartment, plus a small fixed-volume CSF compartment that exchanges with the central compartment through a very low inter-compartmental clearance scaled by a partition multiplier, and that additionally loses drug through the patient's charted CSF drain."
  reference <- "Lu C, Zhang Y, Chen M, Zhong P, Chen Y, Yu J, Wu X, Wu J, Zhang J. Population pharmacokinetics and dosing regimen optimization of meropenem in cerebrospinal fluid and plasma in patients with meningitis after neurosurgery. Antimicrob Agents Chemother. 2016;60(11):6619-6625. doi:10.1128/AAC.00997-16"
  vignette <- "Lu_2016_meropenem"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L"
  )

  # Issue #482: what molecule each compartment holds, in what units, in what
  # biological matrix. Verified against Lu 2016 Fig. 1 (pharmacokinetic
  # structural model for meropenem) and the Table 2 parameter definitions.
  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    csf = list(analyte = "meropenem", units = "mg", specimen = "CSF", verified = TRUE)
  )

  # The ONE covariate the final model consumes. It is not a screened covariate
  # effect on a structural parameter -- it is a measured prescription quantity
  # that enters model() directly as the CSF compartment's drainage clearance
  # (the same role that RRT_PERIT_DIAL_FILL_VOLUME's register entry sanctions).
  covariateData <- list(
    CSF_DRAIN_VOL_24H = list(
      description = "Volume of cerebrospinal fluid removed through the patient's external ventricular or lumbar cistern drain over 24 hours",
      units = "mL/24h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "REQUIRED model input -- the model will not solve without it. Lu 2016 Results,",
        "'Population PK modeling': 'The intercompartmental rate constants were defined as",
        "K12 = Q1 x PC/V1, K21 = Q1/V2, and k20 = daily CSF drainage volume/24.'",
        "The last clause is dimensionally impossible as printed (a volume per day divided",
        "by 24 has units of volume/time, i.e. a CLEARANCE, not the 1/time of a rate",
        "constant), and Fig. 1 separately names BOTH 'CL2, clearance of CSF compartment'",
        "and 'k20, elimination rate constant for CSF compartment'. The model therefore",
        "reads the clause as CL2 = CSF_DRAIN_VOL_24H / 24 (mL/h, converted to L/h here)",
        "with k20 = CL2 / V2. See the vignette 'Assumptions and deviations' for the",
        "numerical discrimination of this reading against the paper's own published",
        "probability-of-target-attainment claims. Table 1 charts the cohort at",
        "126 +/- 81 mL/24h (range 0-350); the Monte Carlo simulations of the paper set it",
        "to 0, 50, 150 and 250 mL/24h, and the paper's headline recommendation is to keep",
        "it below 150 mL/24h. A patient with no drain in situ takes the value 0, which",
        "correctly removes the drainage term entirely."
      ),
      source_name = "CSF daily drainage vol"
    )
  )

  # Lu 2016 Results, 'Population PK modeling': "No covariate was identified to
  # have significant impact on the model from covariate screening." The final
  # model therefore carries NO covariate effect on any structural parameter.
  # The screened-but-not-retained set (Methods, 'Population PK model
  # development') is documented here so the paper's covariate search is
  # preserved without triggering an unused-covariate warning.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened by stepwise forward inclusion / backward deletion and not retained. Cohort mean 43.4 +/- 13.1 years, range 19-77 (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Screened and not retained. 32 of 82 patients (39.0%) were female (Table 1, 'Sex, male/female 50/32')."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened and not retained; the final model carries no allometric term. Cohort mean 65.2 +/- 11.6 kg, range 41.5-100 (Table 1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened and not retained. Cohort mean 23.1 +/- 3.5 kg/m^2, range 13.7-32.7 (Table 1). Height (mean 167.7 +/- 7.2 cm, range 150-180) is tabulated only as the input to BMI and was not itself a screened covariate."
    ),
    SCR = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened and not retained. No summary value is tabulated; the paper reports it only as the input to the Cockcroft-Gault CLCR. Units are not stated in the paper and are recorded here as the canonical mg/dL."
    ),
    CRCL = list(
      description = "Creatinine clearance, Cockcroft-Gault",
      units = "mL/min",
      type = "continuous",
      notes = "Screened and not retained -- a headline negative result. Methods: 'CLCR was estimated from the Cockcroft-Gault equation (33) using the age, body weight, and serum creatinine level of each subject.' Cohort mean 142.6 +/- 52.8 mL/min, range 57.3-355.7 (Table 1) -- an augmented-renal-clearance cohort. The Discussion attributes the negative result to the cohort: 'All subjects we included had normal or mildly impaired renal function, and that might be the reason that CRCL was not identified as a covariate.' Severe renal dysfunction (CLCR <= 10 mL/min) was an exclusion criterion, so the model carries no information about renal impairment."
    ),
    ALT = list(
      description = "Serum alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened and not retained. No summary value is tabulated. Decompensated liver disease (Child-Pugh B or C) was an exclusion criterion."
    ),
    CSF_WBC = list(
      description = "White blood cell count in cerebrospinal fluid",
      units = "10^6/L",
      type = "continuous",
      notes = "Screened and not retained. Cohort mean 2,139.5 +/- 2,877.7 x 10^6/L, range 1-20,000 (Table 1). Clinically load-bearing outside the model: the dosing regimen a patient received was chosen partly on this count (Methods, 'Dosing regimens' -- 1 g q8h below 1,000 x 10^6/L, 1 g q6h or 2 g q8h above it), and a count above 300 x 10^6/L was part of the probable-meningitis inclusion criterion. Documentation only: no register entry was created, because the column is never referenced in model()."
    ),
    CSF_RBC = list(
      description = "Red blood cell count in cerebrospinal fluid",
      units = "10^6/L",
      type = "continuous",
      notes = "Screened and not retained. No summary value is tabulated. Documentation only; no register entry created."
    ),
    CSF_NEUT = list(
      description = "Absolute neutrophil count in cerebrospinal fluid",
      units = "10^6/L",
      type = "continuous",
      notes = "Screened and not retained. No summary value is tabulated. Documentation only; no register entry created."
    ),
    CSF_GLU = list(
      description = "Glucose concentration in cerebrospinal fluid",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened and not retained. Cohort mean 2.4 +/- 1.6 mmol/L, range 0.3-7.9 (Table 1). Documentation only; no register entry created."
    ),
    CSF_TPRO = list(
      description = "Total protein concentration in cerebrospinal fluid",
      units = "g/L",
      type = "continuous",
      notes = "Screened and not retained. Cohort mean 1.9 +/- 1.5 g/L, range 0.2-8.7 (Table 1) -- markedly raised, as expected with an inflamed blood-CSF barrier. This is the one screened CSF analyte that already carries a canonical register entry."
    ),
    CSF_CHLORIDE = list(
      description = "Chloride concentration in cerebrospinal fluid",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened and not retained. No summary value is tabulated. Documentation only; no register entry created."
    ),
    CONMED_MANNITOL = list(
      description = "Concomitant mannitol indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Screened and not retained. Counts are not tabulated. Documentation only; no register entry created."
    ),
    CONMED_DEXAMETHASONE = list(
      description = "Concomitant dexamethasone indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Screened and not retained. Counts are not tabulated. Documentation only; no register entry created."
    ),
    CONMED_VANCOMYCIN = list(
      description = "Concomitant vancomycin (or norvancomycin) indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Screened and not retained. The most common co-medication: 59 of 82 patients (72.0%) received vancomycin or norvancomycin (Table 1, 'Combination treatment')."
    ),
    CONMED_FOSFOMYCIN = list(
      description = "Concomitant fosfomycin indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Screened and not retained. Counts are not tabulated; Table 1 groups 7 patients under 'Combination treatment, Other'. Documentation only; no register entry created."
    ),
    CONMED_NIMODIPINE = list(
      description = "Concomitant nimodipine indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Screened and not retained. Counts are not tabulated. Documentation only; no register entry created."
    ),
    CONMED_VALPROATE = list(
      description = "Concomitant sodium valproate indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Screened and not retained. Counts are not tabulated. Documentation only; no register entry created. Clinically notable because carbapenems are known to depress valproate concentrations, but the paper screened it only as a covariate on meropenem PK."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 82L,
    n_studies = 1L,
    age_range = "19-77 years",
    age_mean = "43.4 years (SD 13.1)",
    weight_range = "41.5-100 kg",
    weight_mean = "65.2 kg (SD 11.6)",
    sex_female_pct = 39.0,
    disease_state = "Adults with probable or proven bacterial meningitis following neurosurgery (65 tumour, 8 trauma, 9 other underlying disease). Febrile at baseline (body temperature mean 38.9 +/- 0.6 degrees C, range 36.5-40.9); CSF white cell count mean 2,139.5 +/- 2,877.7 x 10^6/L; CSF protein mean 1.9 +/- 1.5 g/L; blood white cell count mean 12.6 +/- 4.4 x 10^9/L. Only 8 of 82 had a positive CSF culture (6 Acinetobacter baumannii, 1 Klebsiella pneumoniae, 1 Stenotrophomonas maltophilia).",
    renal_function = "Cockcroft-Gault creatinine clearance mean 142.6 +/- 52.8 mL/min, range 57.3-355.7 -- normal to augmented. Creatinine clearance <= 10 mL/min and haemodialysis were exclusion criteria, so the model carries no information about renal impairment.",
    dose_range = "Meropenem 1 g every 8 h (42 patients), 1 g every 6 h (19 patients) or 2 g every 8 h (21 patients) by intravenous infusion at a fixed rate of 1 g/h, for 3 to 14 days; total daily dose 3-6 g. The 1 g and 2 g doses were dissolved in 100 mL and 250 mL of 0.9% sodium chloride respectively.",
    regions = "Single centre: Huashan Hospital Affiliated to Fudan University, Shanghai, China. ClinicalTrials.gov NCT02506686.",
    notes = "315 plasma and 297 CSF meropenem concentrations from 82 patients. Blood and CSF were sampled simultaneously after the fourth dose; patients were randomised to one of two sampling schedules (group 1, n = 40: during infusion and 10 min, 2 h and 4 h after the end of infusion; group 2, n = 42: end of infusion, 1 h and 3 h after the end of infusion, and immediately before the next dose). CSF was collected through a lumbar cistern or external ventricular drain, and the charted daily drainage volume was 126 +/- 81 mL (range 0-350). Total (not unbound) meropenem was assayed by HPLC-UV. Estimation used FOCE-I in NONMEM 7.3; the reported confidence intervals are the 2.5th-97.5th percentiles of 1,000 bootstrap replicates. Eta shrinkage ranged from 7.75% to 69.7% and epsilon shrinkage was 12.2%. For the paper's PK/PD analysis the plasma protein binding of meropenem was taken as 2% (unbound fraction 0.98) and all meropenem in CSF was assumed unbound."
  )

  ini({
    # ---- Plasma disposition (Lu 2016 Table 2, "Estimate, Mean" column) ----
    lcl <- log(22.2)
    label("Clearance from the central compartment (L/h)")
    # Table 2: CL1 22.2 L/h; bootstrap 2.5th-97.5th percentile 20.5-24.0
    lvc <- log(17.9)
    label("Central compartment volume (L)")
    # Table 2: V1 17.9 L; bootstrap 16.1-19.5
    lq <- log(1.79)
    label("Inter-compartmental clearance central <-> peripheral (L/h)")
    # Table 2: Q2 1.79 L/h; bootstrap 1.21-2.99
    lvp <- log(3.84)
    label("Peripheral compartment volume (L)")
    # Table 2: V3 3.84 L; bootstrap 3.04-4.95

    # ---- CSF compartment (Lu 2016 Fig. 1 and Table 2) ----
    # V2 is FIXED but nonetheless carries an estimated between-subject
    # variability in Table 2; that is the paper as published, not a
    # transcription error. Results: 'The initial estimation of population CSF
    # volume was set at 0.15 liter (35), with the upper and lower limits from
    # 0.13 to 0.17 liter to reflect potential difference. During the
    # model-building process, 0.13 liter fit the data best, so the CSF
    # compartment volume was fixed to 0.13 liter.'
    lvcsf <- fixed(log(0.13))
    label("CSF compartment volume (L)")
    # Table 2: V2 0.13 L, no bootstrap interval reported (fixed)
    lqcsf <- log(0.010)
    label("Inter-compartmental clearance central <-> CSF (L/h)")
    # Table 2: Q1 0.010; bootstrap 0.010-0.010. Table 2's unit header reads
    # '(1/h)' but the Abstract states it as a clearance -- 'The central,
    # intercentral/peripheral, and intercentral/CSF compartment clearances were
    # 22.2 liters/h, 1.79 liters/h, and 0.01 liter/h, respectively' -- and the
    # paper's own K12 = Q1 x PC/V1 and K21 = Q1/V2 only balance dimensionally
    # if Q1 is L/h. The '(1/h)' header is a typographical error; see the
    # vignette Errata.
    lkp_csf <- log(0.172)
    label("CSF:plasma partition multiplier for transfer into the CSF compartment (unitless)")
    # Table 2: PC 0.172; bootstrap 0.140-0.220. Results: 'PC is defined as
    # transfer multiplier between the central and CSF compartments.' It is
    # provably the equilibrium CSF:plasma partition coefficient: setting
    # K12 x A_central = K21 x A_csf with the paper's own definitions gives
    # PC x Cc = Ccsf at equilibrium.

    # ---- Between-subject variability ----
    # Table 2's second block is headed 'Between-subject variability (%)', so
    # the entries are coefficients of variation. They are converted with the
    # house convention omega^2 = log(CV^2 + 1) (as in Luu_2017_nusinersen.R
    # and Stott_2023_flucytosine.R). See the vignette 'Assumptions and
    # deviations' for the alternative reading (that the printed percentages
    # are omega itself x 100) and why it was not adopted. BSV was estimated on
    # V1, V2, CL1, Q1 and PC only -- V3 and Q2 carry none.
    etalcl ~ 0.048958
    # log(0.224^2 + 1); Table 2 CL1 BSV 22.4%, bootstrap 17.3-26.5
    etalvc ~ 0.017274
    # log(0.132^2 + 1); Table 2 V1 BSV 13.2%, bootstrap 0-20.0
    etalvcsf ~ 0.130919
    # log(0.374^2 + 1); Table 2 V2 BSV 37.4%, bootstrap 0-71.4
    etalqcsf ~ 0.537859
    # log(0.844^2 + 1); Table 2 Q1 BSV 84.4%, bootstrap 50.0-104
    etalkp_csf ~ 0.144305
    # log(0.394^2 + 1); Table 2 PC BSV 39.4%, bootstrap 20.0-48.0

    # ---- Residual error ----
    # Table 2 reports ONE pooled 'Residual error (%)' of 34.9 (bootstrap
    # 31.6-38.7) covering both the plasma and the CSF observations. nlmixr2
    # requires one residual term per endpoint, so the single published value is
    # carried on both. The two identical numbers are the paper as published,
    # not a copy-paste error.
    propSd <- 0.349
    label("Proportional residual SD on plasma meropenem Cc (fraction)")
    # Table 2: Residual error 34.9%; bootstrap 31.6-38.7
    propSd_Ccsf <- 0.349
    label("Proportional residual SD on CSF meropenem Ccsf (fraction)")
    # Table 2: the same single pooled 34.9% residual error
  })

  model({
    # Individual parameters. No covariate enters any structural parameter:
    # 'No covariate was identified to have significant impact on the model from
    # covariate screening' (Results). V3 and Q2 carry no eta.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    vcsf <- exp(lvcsf + etalvcsf)
    qcsf <- exp(lqcsf + etalqcsf)
    kp_csf <- exp(lkp_csf + etalkp_csf)

    # Concentrations (mg/L). States are amounts in mg.
    Cc <- central / vc # Fig. 1: concentration in the central compartment
    Cp <- peripheral1 / vp # Fig. 1: concentration in the peripheral compartment
    Ccsf <- csf / vcsf # Fig. 1: concentration in the CSF compartment

    # Lu 2016 Results: 'The intercompartmental rate constants were defined as
    # K12 = Q1 x PC/V1, K21 = Q1/V2, and k20 = daily CSF drainage volume/24.'
    # K12 moves drug central -> CSF and K21 moves it CSF -> central; the
    # asymmetry introduced by PC is what makes Ccsf/Cc equilibrate at PC
    # rather than at 1.
    k12csf <- qcsf * kp_csf / vc # Results: K12 = Q1 x PC / V1
    k21csf <- qcsf / vcsf # Results: K21 = Q1 / V2

    # Drainage clearance of the CSF compartment, in L/h, from the patient's
    # charted 24-hour drain volume in mL. Fig. 1 names this 'CL2, clearance of
    # CSF compartment' and 'k20, elimination rate constant for CSF
    # compartment', with k20 = CL2 / V2. The Results clause 'k20 = daily CSF
    # drainage volume/24' is dimensionally a clearance, not a rate constant,
    # and is read here as CL2; see covariateData and the vignette Errata.
    clcsf <- CSF_DRAIN_VOL_24H / 1000 / 24

    d/dt(central) <- -cl * Cc - q * Cc + q * Cp - k12csf * central + k21csf * csf
    d/dt(peripheral1) <- q * Cc - q * Cp
    d/dt(csf) <- k12csf * central - k21csf * csf - clcsf * Ccsf

    Cc ~ prop(propSd)
    Ccsf ~ prop(propSd_Ccsf)
  })
}
