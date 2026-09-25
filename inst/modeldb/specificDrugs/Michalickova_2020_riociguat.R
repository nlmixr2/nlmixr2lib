Michalickova_2020_riociguat <- function() {
  description <- "Joint parent (riociguat) + metabolite (M1, desmethylriociguat) population PK model in adults with chronic thromboembolic pulmonary hypertension (CTEPH) treated in routine clinical practice (Michalickova 2020). Riociguat is one-compartment with first-order absorption (ka fixed at 3 1/h from the literature) and two parallel first-order elimination routes from the central compartment: metabolic formation of M1 (CLf,M1/F, power function of total bilirubin centred at 0.69 mg/dL) and all remaining pathways (CLe,r/F, linear in creatinine clearance centred at 70 mL/min). M1 is one-compartment with first-order elimination (CLe,M1/F) and shares the parent's apparent volume of distribution (assumed for identifiability). An absorption lag time of 2.95 h applies only to the six patients whose late post-dose concentrations were unexpectedly high (MIX_LAGGED_ABS = 1). With one sample per patient, inter-individual and residual variability could not be separated, so the model carries no etas and the proportional residual errors absorb both."
  reference <- "Michalickova D, Jansa P, Bursova M, Hlozek T, Cabala R, Hartinger JM, Ambroz D, Aschermann M, Lindner J, Linhart A, Slanar O, Krekels EHJ. Population pharmacokinetics of riociguat and its metabolite in patients with chronic thromboembolic pulmonary hypertension from routine clinical practice. Pulm Circ. 2020;10(1):2045894019898031. doi:10.1177/2045894019898031. The absorption rate constant is fixed to a literature value (Michalickova 2020 ref 17: Saleh S, Becker C, Frey R, et al. Population pharmacokinetics of single-dose riociguat in patients with renal or hepatic impairment. Pulm Circ. 2016;6:S75-S85)."
  vignette <- "Michalickova_2020_riociguat"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance (kidney-function estimate by the CKD-EPI equation)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed (one sample per patient). Michalickova 2020 names the covariate 'creatinine clearance (calculated by CKD-EPI (Chronic Kidney Disease Epidemiology Collaboration) equation)' and reports it in mL/min (Methods 'Covariate analysis'; Table 1 median 70, IQR 59-79 mL/min); whether the CKD-EPI value was de-normalised from mL/min/1.73 m^2 is not stated. Enters the remaining (non-M1) riociguat clearance through the origin as CLe,r/F = 0.66 * (CRCL / 70) (Table 2), i.e. 0.0094 L/h per mL/min -- the Results quote 'increase 0.009 L/h per unit (mL/min) creatinine clearance'.",
      source_name = "CREACL"
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. The canonical TBILI column is in umol/L; Michalickova 2020 modelled total bilirubin in mg/dL, so model() converts with TBILI / 17.1 before applying the published power function CLf,M1/F = 0.665 * (BILTOT / 0.69)^-0.462 (Table 2). The 0.69 mg/dL reference (11.8 umol/L) is the cohort median (Table 1, IQR 0.53-0.98 mg/dL).",
      source_name = "BILTOT (mg/dL)"
    ),
    MIX_LAGGED_ABS = list(
      description = "Delayed-absorption subgroup indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no absorption lag)",
      notes = "Time-fixed per subject. 1 = one of the six patients (ID 6, 10, 11, 19, 25 and 28) whose riociguat concentrations were unexpectedly high at relatively late times after the last dose (Figure 1); for these patients only, an absorption lag time of 2.95 h was estimated (Table 2 row 'Tlag (ID = 6,10,11,19,25,28)'). The class was assigned deterministically by visual inspection of the raw data, not by a NONMEM $MIXTURE block (Methods 'Population PK analysis'). The authors suggest food may be the cause but state their design cannot confirm it. For typical-patient simulation set MIX_LAGGED_ABS = 0; to reproduce the source cohort draw MIX_LAGGED_ABS ~ Bernoulli(6/49 = 0.122).",
      source_name = "ID in (6, 10, 11, 19, 25, 28)"
    )
  )

  # Covariates Michalickova 2020 screened on every structural parameter
  # (CLf,M1/F, CLe,r/F, Vd/F and CLe,M1/F) but did not retain once the two
  # retained relationships were included (Methods 'Covariate analysis';
  # Results 'Population PK analysis'). Documentation only -- none is referenced
  # in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested in linear and exponential equations; not retained. Cohort median 74 (IQR 66-78) years (Table 1)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested in linear and exponential equations; not retained. Cohort median 80 (IQR 67-95) kg (Table 1)."
    ),
    IBW = list(
      description = "Ideal body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested in linear and exponential equations; not retained. Cohort median 61 (IQR 55-74) kg (Table 1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Tested in linear and exponential equations; not retained. Cohort median 27.8 (IQR 23.8-30.8) kg/m^2 (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "(as reported)",
      type = "continuous",
      notes = "Tested in linear and exponential equations; not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "(as reported)",
      type = "continuous",
      notes = "Tested in linear and exponential equations; not retained."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "(as reported)",
      type = "continuous",
      notes = "Tested in linear and exponential equations; not retained."
    ),
    GGT = list(
      description = "Gamma-glutamyl transferase",
      units = "(as reported)",
      type = "continuous",
      notes = "Tested in linear and exponential equations; not retained."
    ),
    NTPROBNP = list(
      description = "N-terminal pro b-type natriuretic peptide",
      units = "ng/L",
      type = "continuous",
      notes = "Tested in linear and exponential equations; not retained. Cohort median 493 (IQR 235-1460) ng/L (Table 1)."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor",
      units = "(binary)",
      type = "binary",
      notes = "Tested as a fractional effect; not retained. 35 of 49 patients (71%) (Table 1)."
    ),
    CONMED_DIURETIC = list(
      description = "Concomitant diuretic",
      units = "(binary)",
      type = "binary",
      notes = "Tested as a fractional effect; not retained. 35 of 49 patients (71%) (Table 1)."
    ),
    CONMED_DIGOXIN = list(
      description = "Concomitant digoxin",
      units = "(binary)",
      type = "binary",
      notes = "Tested as a fractional effect; not retained. 3 of 49 patients (6.1%) (Table 1)."
    ),
    CONMED_ACEI = list(
      description = "Concomitant angiotensin-converting enzyme inhibitor",
      units = "(binary)",
      type = "binary",
      notes = "Tested as a fractional effect; not retained. 11 of 49 patients (22.4%) (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested as a fractional effect; not retained. 24 of 49 patients (49%) female (Table 1)."
    ),
    SMOKE = list(
      description = "Current smoker indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested as a fractional effect; not retained. Only 1 of 49 patients (2%) smoked (Table 1)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "riociguat", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "riociguat", units = "mg", specimen = "serum", verified = TRUE),
    central_m1 = list(analyte = "M1 (desmethylriociguat)", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 49L,
    n_studies = 1L,
    age_range = "median 74 years (IQR 66-78)",
    weight_range = "median 80 kg (IQR 67-95)",
    sex_female_pct = 49,
    race_ethnicity = "(not reported; single-centre Czech cohort)",
    disease_state = "Chronic thromboembolic pulmonary hypertension: inoperable CTEPH (74%) or persistent/recurrent pulmonary hypertension after pulmonary endarterectomy (26%); NYHA class II 42.8%, class III 57.2%.",
    dose_range = "Riociguat 1.5-2.5 mg orally three times daily at a stable dose for at least 3 months (median 7.5 mg/day, IQR 6.75-7.5).",
    regions = "Czech Republic (General University Hospital, Prague)",
    renal_function = "Creatinine clearance (CKD-EPI) median 70 mL/min (IQR 59-79)",
    hepatic_function = "Total bilirubin median 0.69 mg/dL (IQR 0.53-0.98)",
    co_medication = "Warfarin 84%, diuretics 71%, proton-pump inhibitors 71%, ACE inhibitors 22.4%, ARBs 8.2%, digoxin 6.1%, H2 antihistamines 4%",
    n_observations = "One serum riociguat and one M1 concentration per patient (49 of each), sampled 1.25-6.75 h after the last dose at steady state.",
    notes = "Observational routine-practice study (Michalickova 2020 Methods 'Research design'; Table 1). Riociguat concentrations ranged 44-749 ug/L and M1 17-314 ug/L (Results 'Study population and data')."
  )

  ini({
    lka <- fixed(log(3)); label("Absorption rate constant ka (1/h)") # Table 2: 'Ka/F (h-1) 3 FIX'; fixed from the literature (Methods, ref 17)
    lcl_nonmet <- log(0.66); label("Apparent riociguat clearance through pathways other than M1 formation, CLe,r/F, at CRCL = 70 mL/min (L/h)") # Table 2: CLe,rTV = 0.66 L/h (17% RSE)
    lcl_met <- log(0.665); label("Apparent riociguat metabolic formation clearance to M1, CLf,M1/F, at total bilirubin 0.69 mg/dL (L/h)") # Table 2: CLf,M1TV = 0.665 L/h (17% RSE)
    e_tbili_cl_met <- -0.462; label("Power exponent of total bilirubin (mg/dL / 0.69) on CLf,M1/F (unitless)") # Table 2: thetaBILTOT = -0.462 (20% RSE); Results text quotes -0.463
    lvc <- log(3.63); label("Apparent volume of distribution of riociguat, shared by M1 (VP/F = VM/F) (L)") # Table 2: 'VP/F = VM/F (L)' = 3.63 (15% RSE)
    lcl_m1 <- log(1.47); label("Apparent M1 elimination clearance CLe,M1/F (L/h)") # Table 2: CLe,M1TV = 1.47 L/h (19% RSE)
    ltlag <- log(2.95); label("Absorption lag time in the delayed-absorption subgroup (h)") # Table 2: 'Tlag (ID = 6,10,11,19,25,28) (h)' = 2.95 (6% RSE)

    # Table 2 reports the residual-variability VARIANCES (0.152 and 0.268);
    # the SDs below are their square roots. With a single sample per patient
    # these terms represent inter-individual AND residual variability together.
    propSd <- 0.3899; label("Proportional residual error on riociguat concentration (fraction)") # Table 2: riociguat proportional variance 0.152 -> SD sqrt(0.152) = 0.3899
    propSd_m1 <- 0.5177; label("Proportional residual error on M1 concentration (fraction)") # Table 2: M1 proportional variance 0.268 -> SD sqrt(0.268) = 0.5177
  })

  model({
    # Total bilirubin in the paper's mg/dL units (canonical TBILI is umol/L)
    tbili_mgdl <- TBILI / 17.1

    ka <- exp(lka)
    # Remaining riociguat clearance: linear in creatinine clearance through the
    # origin (Table 2: CLe,r/F = CLe,rTV * (CREACL/70))
    cl_nonmet <- exp(lcl_nonmet) * (CRCL / 70)
    # Formation clearance to M1: power function of total bilirubin
    # (Table 2: CLf,M1/F = CLf,M1TV * (BILTOT/0.69)^thetaBILTOT)
    cl_met <- exp(lcl_met) * (tbili_mgdl / 0.69)^e_tbili_cl_met
    cl_m1 <- exp(lcl_m1)
    # Shared apparent volume for parent and metabolite (Methods: 'The same
    # values of volume of distribution (Vd) of riociguat and its metabolite
    # were assumed')
    vc <- exp(lvc)
    vc_m1 <- vc
    # Lag time only in the six late-high-concentration patients
    tlag <- exp(ltlag) * MIX_LAGGED_ABS

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (cl_met + cl_nonmet) / vc * central
    # M1 formed 1:1 on a mass basis (no molecular-weight correction stated)
    d/dt(central_m1) <- cl_met / vc * central - cl_m1 / vc_m1 * central_m1
    alag(depot) <- tlag

    # Dose mg / volume L = mg/L; x 1000 -> ug/L (the paper's reporting unit)
    Cc <- central / vc * 1000
    Cc_m1 <- central_m1 / vc_m1 * 1000

    Cc ~ prop(propSd)
    Cc_m1 ~ prop(propSd_m1)
  })
}
