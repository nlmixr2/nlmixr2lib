Roberts_2011_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin administered by continuous infusion in adult septic critically ill ICU patients (Roberts 2011). Volume of distribution scales linearly with total body weight (1.53 L/kg); clearance scales linearly with BSA-normalized 24-hour urinary creatinine clearance referenced to 100 mL/min/1.73 m^2 (4.58 L/h at the reference)."
  reference <- "Roberts JA, Taccone FS, Udy AA, Vincent JL, Jacobs F, Lipman J. Vancomycin dosing in critically ill patients: robust methods for improved continuous-infusion regimens. Antimicrob Agents Chemother. 2011;55(6):2704-2709. doi:10.1128/AAC.01708-10"
  vignette <- "Roberts_2011_vancomycin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Roberts 2011 Table 1: mean 74.8 kg (SD 15.8). Enters the volume-of-distribution structural parameter as a linear factor with no reference scaling (Eq. 3: TVV = theta1 * TBW; theta1 = 1.53 L/kg). Source column TBW.",
      source_name        = "TBW"
    ),
    CRCL = list(
      description        = "Creatinine clearance from 24-hour urinary collection, normalized to body surface area",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Roberts 2011 Methods (Patients and data collection): urinary CrCl collected daily as a routine procedure and normalized to body surface area. Population mean 90.7 mL/min/1.73 m^2 (SD 60.4) per Table 1. Enters CL as a linear factor referenced to 100 mL/min/1.73 m^2 (Eq. 4: TVCL = theta2 * CrCl/100; theta2 = 4.58 L/h). Source column CrCl.",
      source_name        = "CrCl"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 206L,
    n_studies        = 1L,
    age_range        = "adult (>=18 years)",
    age_median       = "58.1 years (mean, SD 14.8)",
    weight_range     = "mean 74.8 kg (SD 15.8)",
    weight_median    = "74.8 kg (mean)",
    sex_female_pct   = 38.4,
    race_ethnicity   = "Not reported (single-center study in Brussels, Belgium)",
    disease_state    = "Adult septic critically ill patients in the intensive care unit, treated empirically or for documented Gram-positive infection (frequently MRSA, methicillin-resistant Staphylococcus epidermidis, or ampicillin-resistant Enterococcus). Exclusion criteria: age <18 years, prior intermittent vancomycin in the 48 h before continuous-infusion onset, renal replacement therapy, continuous-infusion duration <48 h, pregnancy, burns, or cystic fibrosis. APACHE II score 21 (16-27); SOFA score 7.6 (SD 4.2).",
    dose_range       = "Vancomycin by continuous infusion. Loading dose 15 mg/kg over 30 min (rounded to 125 mg) or simplified flat 750 mg (TBW <70 kg) / 1000 mg (TBW >=70 kg) over 30 min; maintenance 30 mg/kg/day or simplified flat 2000 mg/day (TBW <70 kg) / 3000 mg/day (TBW >=70 kg) by 24 h continuous infusion, adapted for renal failure. Per-patient target serum concentration 20-30 mg/L.",
    regions          = "Belgium (single center, Erasme Hospital, Brussels)",
    apache_ii        = "median 21 (IQR 16-27)",
    sofa             = "mean 7.6 (SD 4.2)",
    renal_function   = "BSA-normalized 24-hour urinary creatinine clearance mean 90.7 mL/min/1.73 m^2 (SD 60.4); patients on renal replacement therapy excluded",
    n_concentrations = 579L,
    notes            = "Retrospective single-center cohort, January 2008 to December 2009. Vancomycin assay by fluorescence polarization immunoassay (TDx; Abbott); assay limit 0.6 mg/L. Each patient contributed 2-3 trough samples (daily at 8 a.m., at least 16 h after onset of continuous infusion to capture pseudo-steady state). Model fit in NONMEM 6.1 (FOCEI) via Wings for NONMEM 6.1.3. Bootstrap N = 1000 used to derive 95% confidence intervals (Roberts 2011 Table 2). Covariates screened and not retained: age, sex, SOFA score, BMI."
  )

  ini({
    # Structural parameters (Roberts 2011 Table 2 final-covariate-model column,
    # bootstrap means). Reference subject is a 1-kg-TBW patient (per-kg V slope)
    # with CrCl = 100 mL/min/1.73 m^2.
    lcl <- log(4.58);  label("Clearance at CRCL = 100 mL/min/1.73 m^2 (L/h)")     # Roberts 2011 Eq. 4 + Table 2: CL = 4.58 L/h (95% CI 4.09-5.19)
    lvc <- log(1.53);  label("Volume of distribution per kg total body weight (L/kg)") # Roberts 2011 Eq. 3 + Table 2: V = 1.53 L/kg (95% CI 1.31-1.71)

    # Inter-individual variability (Roberts 2011 Table 2 final-covariate-model
    # %CV rows; exponential variability model per Eq. 2). For log-normal etas,
    # omega^2 = log(CV^2 + 1).
    etalcl ~ 0.140924 # log(0.389^2 + 1); 38.9% CV on CL (Roberts 2011 Table 2: 95% CI 28.3-55.6)
    etalvc ~ 0.130918 # log(0.374^2 + 1); 37.4% CV on V  (Roberts 2011 Table 2: 95% CI 16.6-54.9)

    # Combined proportional + additive residual error (Roberts 2011 Table 2
    # Random error rows).
    propSd <- 0.199; label("Proportional residual error (fraction)") # Roberts 2011 Table 2: residual unexplained variability CV = 19.9% (95% CI 14.5-24.6)
    addSd  <- 2.4;   label("Additive residual error (mg/L)")         # Roberts 2011 Table 2: residual additive SD = 2.4 mg/L (95% CI 1.3-3.0)
  })
  model({
    # Individual PK parameters. CL is linear in BSA-normalized CrCl referenced to
    # 100 mL/min/1.73 m^2; V is linear in TBW with no reference scaling (the
    # paper reports V as 1.53 L/kg directly).
    cl <- exp(lcl + etalcl) * CRCL / 100
    vc <- exp(lvc + etalvc) * WT

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volumes in L -> central/vc has units mg/L. Drug enters the
    # central compartment directly as continuous infusion (rate set on the
    # dosing record); a loading dose is administered as a separate IV
    # infusion to the same compartment.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
