Revilla_2010_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in critically ill adult medical ICU patients (Revilla 2010). Clearance is the sum of a renal arm proportional to weight-normalised creatinine clearance and a non-renal arm scaling as AGE^-0.24; central volume of distribution is per kg with a >2-fold increase when serum creatinine exceeds 1 mg/dL."
  reference <- "Revilla N, Martin-Suarez A, Paz Perez M, Martin Gonzalez F, Fernandez de Gatta MM. Vancomycin dosing assessment in intensive care unit patients based on a population pharmacokinetic/pharmacodynamic simulation. Br J Clin Pharmacol. 2010;70(2):201-212. doi:10.1111/j.1365-2125.2010.03679.x"
  vignette <- "Revilla_2010_vancomycin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Revilla 2010 Table 2: mean 73.0 kg (SD 13.3, range 45-150). Used to convert the paper's weight-normalised CL (mL/min/kg) and V (L/kg) into total L/h and L. Body weight was added to the model before the evaluation of other covariates per Discussion paragraph 4.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age at start of vancomycin therapy",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Revilla 2010 Table 2: mean 61.1 years (SD 16.3, range 18-85). Enters the non-renal CL arm as AGE^theta2 with theta2 = -0.24 (Table 4 final-model column).",
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = "24-hour measured creatinine clearance (Levey-estimated CrCl substituted for the ~8 percent of records with no 24-hour collection); raw mL/min, not BSA-normalised",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCR. Revilla 2010 Table 2: mean 74.7 mL/min (SD 58.0, range 10-328); the 24-hour measured CrCl was preferred and the Levey-formula estimate was substituted only when the measurement was unavailable (Methods, Patients and study design, ~8 percent of data). Stored under the canonical CRCL covariate per the inst/references/covariate-columns.md precedent that accepts raw mL/min when the source does not BSA-normalise (Delattre 2010 amikacin, MedellinGaribay 2015 gentamicin). Enters the renal CL arm as 0.67 * CRCL/WT in mL/min/kg.",
      source_name        = "CLCR"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CrSe. Revilla 2010 Table 2: mean 1.4 mg/dL (SD 1.0, range 0.6-5.0). A 0.5 mg/dL quantification limit was assumed for values below the assay limit. Enters V via a dichotomous indicator A = (CREAT > 1 mg/dL): V multiplies by theta4 = 2.49 when A = 1 (Table 4 final-model column; threshold direction confirmed by Results paragraph below Table 3 and Discussion paragraph 4).",
      source_name        = "CrSe"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 191L,
    n_studies        = 1L,
    age_range        = "18-85 years",
    age_median       = "61.1 years (mean)",
    weight_range     = "45-150 kg",
    weight_median    = "73.0 kg (mean)",
    sex_female_pct   = 34,
    race_ethnicity   = "Not reported (single-centre Spanish ICU cohort)",
    disease_state    = "Adult medical ICU patients receiving vancomycin (severe trauma n=81, post-surgery n=50, sepsis n=49 of whom n=13 septic shock, respiratory infections / pneumonia n=66, multiorgan failure n=35, respiratory distress syndrome n=30, hypovolaemic shock n=15)",
    dose_range       = "Vancomycin IV; mostly 60-min infusions (n=406 episodes; 42% 1000 mg q12h and 20% 1000 mg q24h, mean 18.4 mg/kg/day); 14 episodes continuous infusion (mean rate 56.9 mg/h); no patient received a loading dose",
    regions          = "Spain (single centre, University Hospital of Salamanca, 1999-2004)",
    apache_ii        = "mean 18.0 (SD 6.9, range 2-41)",
    serum_albumin    = "mean 2.3 g/dL (SD 0.7, range 0.5-6.2)",
    renal_function   = "24-hour measured CrCl mean 74.7 mL/min (SD 58.0, range 10-328); raw mL/min, not BSA-normalised. Patients on renal replacement therapy were excluded.",
    n_concentrations = 569L,
    notes            = "Retrospective non-comparative cohort. 569 concentration-time records from 191 patients (mean 2.98 per patient, range 1-19), of which ~80 percent are C(min) trough samples (0-60 min pre-dose). External validation cohort: 46 patients / 73 concentrations (2007-2008), mean age 58.7 +/- 16.6 years, mean weight 73.1 +/- 19.8 kg. Excluded from model building: neoplasic disorders, prior cardiac surgery, renal replacement therapy, and patients without concentration-time data. Mechanical ventilation 87%; parenteral nutrition 46%; co-medications albumin 21%, furosemide 64%, catecholamines 71%. Demographics in Table 2. Fit in NONMEM v5 level 1.1 (FOCE INTERACTION, ADVAN1 TRANS2). External evaluation: standardised prediction errors 0.14 +/- 0.70 mg/L (95% CI -0.03 to 0.30 included zero); 100% of observed concentrations fell within PRED +/- 2 * SDpop."
  )

  ini({
    # Revilla 2010 Table 4 final-model column (OF = 2420.69). Structural form
    # (Results paragraph below Table 3 / above Table 4):
    #   CL (mL/min/kg) = theta1 * CLcr (mL/min/kg) + AGE(years)^theta2
    #   V  (L/kg)      = theta3 * theta4^A; A = 1 if CrSe > 1 mg/dL else 0
    # The renal-arm coefficient theta1 multiplies the weight-normalised
    # CRCL; the non-renal arm has an implicit unit coefficient of 1 mL/min/kg
    # (i.e. the value of AGE^theta2 evaluated at AGE = 1 year), which we
    # encode explicitly as fixed(log(1)) so the lcl_renal / lcl_nonren
    # additive multi-component-CL pattern parses cleanly under the package
    # conventions. Multiplying by WT and 60/1000 converts mL/min/kg to L/h.
    lcl_renal       <- log(0.67); label("Slope of weight-normalised renal CL on weight-normalised CRCL (mL/min/kg per mL/min/kg)") # Revilla 2010 Table 4 theta1 = 0.67 (RSE 6.6%; 95% CI 0.76, 0.58)
    lcl_nonren      <- fixed(log(1)); label("Non-renal CL anchor (mL/min/kg at AGE = 1 year; structurally fixed to 1 per the paper's equation)") # Revilla 2010 implicit unit coefficient on AGE^theta2 (Results paragraph 2 above Table 4)
    e_age_cl_nonren <- -0.24; label("Age exponent on non-renal CL arm (unitless)") # Revilla 2010 Table 4 theta2 = -0.24 (RSE 6.7%; 95% CI -0.21, -0.27)

    lvc        <- log(0.82); label("Central volume of distribution at CrSe <= 1 mg/dL (L/kg)") # Revilla 2010 Table 4 theta3 = 0.82 L/kg (RSE 2.8%; 95% CI 0.94, 0.70)
    e_creat_vc <- log(2.49); label("Log multiplier on V for CrSe > 1 mg/dL (unitless)")        # Revilla 2010 Table 4 theta4 = 2.49 (RSE 9.9%; 95% CI 2.98, 2.00)

    # Inter-individual variability per Revilla 2010 Methods (Eq. for the
    # exponential IIV model) and Table 4 final-model column.
    # omega^2 = log(1 + CV^2) for log-normal etas.
    etalcl ~ 0.08688 # Revilla 2010 Table 4: omega_CL = 30.13% CV (RSE 16.4%; 95% CI 34.64, 24.82); log(1 + 0.3013^2)
    etalvc ~ 0.05083 # Revilla 2010 Table 4: omega_V  = 22.83% CV (RSE 38.8%; 95% CI 30.28, 11.18); log(1 + 0.2283^2)

    # Additive residual error (Revilla 2010 Methods residual error model and
    # Table 4 final-model column).
    addSd <- 4.23; label("Additive residual error (mg/L)") # Revilla 2010 Table 4: residual SD = 4.23 mg/L (RSE 9.8%; 95% CI 4.62, 3.80)
  })

  model({
    # Renal + non-renal CL arms in mL/min/kg, then scaled to L/h by WT and 60/1000.
    cl_renal_per_kg  <- exp(lcl_renal) * (CRCL / WT)
    cl_nonren_per_kg <- exp(lcl_nonren) * AGE^e_age_cl_nonren
    cl <- (cl_renal_per_kg + cl_nonren_per_kg) * WT * (60 / 1000) * exp(etalcl)

    # V (L/kg) = theta3 * theta4^A with A = (CREAT > 1); multiplied by WT to L.
    crse_hi <- (CREAT > 1)
    vc <- exp(lvc + e_creat_vc * crse_hi + etalvc) * WT

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose mg / vc L gives mg/L (= ug/mL); IV-infusion rate is set on the
    # dosing record (Revilla 2010 used 60-min and continuous infusions).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
