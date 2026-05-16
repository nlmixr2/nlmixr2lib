Goti_2018_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in hospitalized adults with and without intermittent hemodialysis (Goti 2018). Volumes scaled allometrically to body weight (reference 70 kg, fixed linear exponent), CL scaled by Cockcroft-Gault creatinine clearance with a power exponent (reference 120 mL/min), and intermittent hemodialysis acts as a multiplicative factor on CL (0.7) and central volume (0.5)."
  reference <- "Goti V, Chaturvedula A, Fossler MJ, Mok S, Jacob JT. Hospitalized patients with and without hemodialysis have markedly different vancomycin pharmacokinetics: a population pharmacokinetic model-based analysis. Ther Drug Monit. 2018;40(2):212-221. doi:10.1097/FTD.0000000000000490"
  vignette <- "Goti_2018_vancomycin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Goti 2018 Table 1: median 79 kg (range 33-255). Volumes scaled by (WT/70) with fixed linear exponent in the structural model; Goti 2018 Methods state 'Volume parameters were normalized to a 70-kg individual before other covariate relationships were examined.' Table 2 covariate-relationship row explicitly shows TVVc = theta4 * (WT/70) * theta5^DIAL; the same (WT/70) factor is applied to the peripheral volume (Vp) per the Methods text, since the Table 2 row only lists covariates retained on top of the structural WT normalisation.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Goti 2018 Methods: Cockcroft-Gault creatinine clearance with the following preprocessing rules applied upstream of the model (a) values greater than 150 mL/min were truncated to 150 mL/min; (b) for subjects older than 60 years, serum creatinine values below 1 mg/dL were truncated to 1 mg/dL before computing CrCL (correction for falsely depressed SCr in elderly with low muscle mass). Reference value 120 mL/min (Table 2 covariate-relationship row: TVCL = theta1 * (CRCL/120)^theta2 * theta3^DIAL). Population median was 62 mL/min (Table 1, range 4-150 mL/min after truncation). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form).",
      source_name        = "CrCL"
    ),
    HEMODIAL = list(
      description        = "Intermittent-hemodialysis treatment-status indicator (1 = subject on intermittent hemodialysis during the admission, 0 = not)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no intermittent hemodialysis)",
      notes              = "Goti 2018 Methods: all hemodialysis procedures were intermittent and used high-flux membranes. 336 of 1812 subjects (18.5%) were on hemodialysis. Time-fixed per subject in the source data because actual hemodialysis-session timestamps were not used (the paper notes documentation limitations prevented per-session timing). Enters the structural model as multiplicative factors theta3^DIAL on CL (theta3 = 0.7) and theta5^DIAL on Vc (theta5 = 0.5), so dialysis subjects have 30% lower CL and 50% lower Vc than non-dialysis subjects.",
      source_name        = "DIAL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1812L,
    n_studies      = 2L,
    age_range      = "17-101 years",
    age_median     = "57 years",
    weight_range   = "33-255 kg",
    weight_median  = "79 kg",
    sex_female_pct = 46.5,
    race_ethnicity = "Not reported (two 500-bed academic medical centers in Atlanta, Georgia; population with a high prevalence of end-stage renal disease)",
    disease_state  = "Hospitalized adults aged >16 years receiving at least one dose of intravenous vancomycin between January 1, 2011 and December 31, 2012; 64.7% had ICU stay during hospitalization; 18.5% (336/1812) were receiving intermittent hemodialysis",
    dose_range     = "Intravenous vancomycin per local nomogram; therapeutic-drug-monitoring trough concentrations observed in routine care (median observed serum concentration 18 mg/L, range 5.1-120)",
    regions        = "United States (Emory University and Mercer University academic medical centers in Atlanta, Georgia)",
    renal_function = "Cockcroft-Gault creatinine clearance median 62 mL/min (range 4-150 after truncation at 150); 18.5% of subjects on intermittent hemodialysis",
    n_concentrations = 2765L,
    notes          = "Baseline demographics from Goti 2018 Table 1. 1920 patients initially eligible; 108 excluded (concentrations >120 hours post-dose or BQL of 5 mg/L). 1812 retained for model development with 2765 vancomycin concentrations; nearly two thirds of patients had only one concentration. Data extracted from electronic medical records (retrospective routine TDM). Goti 2018 explicitly notes that final-model NONMEM minimization 'terminated due to rounding errors' but the covariance step was successful and all parameter correlations were <0.95; the authors evaluated the covariance output and accepted the final-model results -- see vignette Errata."
  )

  ini({
    # Structural parameters (Goti 2018 Table 2 final-model column, typical-value
    # column "Values (%RSE)"). Reference subject is a 70-kg non-dialysis patient
    # with CRCL = 120 mL/min.
    lcl <- log(4.5);  label("Clearance at CRCL = 120 mL/min, non-dialysis (L/h)") # Goti 2018 Table 2: CL = 4.5 L/h
    lvc <- log(58.4); label("Central volume at WT = 70 kg, non-dialysis (L)")     # Goti 2018 Table 2: Vc = 58.4 L
    lvp <- log(38.4); label("Peripheral volume at WT = 70 kg (L)")                # Goti 2018 Table 2: Vp = 38.4 L
    lq  <- log(6.5);  label("Intercompartmental clearance Q (L/h)")               # Goti 2018 Table 2: Q = 6.5 L/h

    # Covariate effects (Goti 2018 Table 2; covariate-relationship row):
    #   TVCL = theta1 * (CRCL/120)^theta2 * theta3^DIAL
    #   TVVc = theta4 * (WT/70)        * theta5^DIAL
    e_crcl_cl     <- 0.8; label("Power exponent on (CRCL/120) for CL")           # Goti 2018 Table 2: CrCL on CL = 0.8
    e_hemodial_cl <- 0.7; label("Multiplicative factor on CL when HEMODIAL = 1") # Goti 2018 Table 2: DIAL on CL = 0.7
    e_hemodial_vc <- 0.5; label("Multiplicative factor on Vc when HEMODIAL = 1") # Goti 2018 Table 2: DIAL on Vc = 0.5

    # Inter-individual variability (Goti 2018 Table 2 final-model %CV rows);
    # omega^2 = log(CV^2 + 1) for log-normal etas.
    etalcl ~ 0.14705 # log(0.398^2 + 1); 39.8% CV on CL
    etalvc ~ 0.51019 # log(0.816^2 + 1); 81.6% CV on Vc
    etalvp ~ 0.28232 # log(0.571^2 + 1); 57.1% CV on Vp

    # Combined additive + proportional residual error (Goti 2018 Table 2).
    addSd  <- 3.4;   label("Additive residual error (mg/L)")               # Goti 2018 Table 2: additive error SD = 3.4 mg/L
    propSd <- 0.227; label("Proportional residual error (fraction)")       # Goti 2018 Table 2: proportional error CV = 22.7%
  })
  model({
    # Individual PK parameters. CL scales by (CRCL/120)^0.8 and by 0.7^HEMODIAL;
    # Vc scales by (WT/70) and by 0.5^HEMODIAL; Vp scales by (WT/70) per the
    # Methods statement that volume parameters were normalized to a 70-kg
    # individual; Q has no covariates.
    cl <- exp(lcl + etalcl) * (CRCL / 120)^e_crcl_cl * e_hemodial_cl^HEMODIAL
    vc <- exp(lvc + etalvc) * (WT / 70) * e_hemodial_vc^HEMODIAL
    vp <- exp(lvp + etalvp) * (WT / 70)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
