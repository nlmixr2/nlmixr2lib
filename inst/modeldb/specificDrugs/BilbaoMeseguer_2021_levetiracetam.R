BilbaoMeseguer_2021_levetiracetam <- function() {
  description <- "Two-compartment IV population PK model for levetiracetam in critically ill adults with normal or augmented renal clearance, with clearance split into a fixed non-renal arm plus a power function of measured urinary creatinine clearance (Bilbao-Meseguer 2021)"
  reference <- "Bilbao-Meseguer I, Barrasa H, Asin-Prieto E, Alarcia-Lacalle A, Rodriguez-Gascon A, Maynar J, Sanchez-Izquierdo JA, Balziskueta G, Griffith MS-B, Quilez Trasobares N, Solinis MA, Isla A. Population Pharmacokinetics of Levetiracetam and Dosing Evaluation in Critically Ill Patients with Normal or Augmented Renal Function. Pharmaceutics. 2021;13(10):1690. doi:10.3390/pharmaceutics13101690"
  vignette <- "BilbaoMeseguer_2021_levetiracetam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "levetiracetam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "levetiracetam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance MEASURED from a urine collection, raw mL/min (NOT BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Bilbao-Meseguer 2021 Table 1 footnote 2: CrCl = urine creatinine (mg/dL) x urine volume per minute (mL/min) / plasma creatinine (mg/dL). Raw mL/min, not normalized to 1.73 m^2. Cohort median 117 mL/min (range 54-239). Enters as (CRCL / 120)^2.5 L/h added to a 3.5 L/h non-renal arm; the 120 mL/min divisor is the value printed in the final-model equation and the Table 3 header, not the cohort median. Baseline value (covariates were assessed at baseline, Section 2.4).",
      source_name = "CrCl"
    )
  )

  covariatesDataExcluded <- list(
    DIAG_TRAUMA = list(
      description = "Trauma vs non-trauma admission diagnosis",
      units = "(binary)",
      type = "binary",
      notes = "Significant on V1 in forward inclusion (p < 0.05) but dropped in backward elimination (p < 0.01 criterion); not in the final model (Section 3.3)."
    ),
    APACHE_II = list(
      description = "APACHE II score at baseline",
      units = "points",
      type = "continuous",
      notes = "Significant on V1 in forward inclusion but dropped in backward elimination; not in the final model (Section 3.3)."
    ),
    RENAL_ARC = list(
      description = "Augmented renal clearance indicator (CrCl >= 130 mL/min)",
      units = "(binary)",
      type = "binary",
      notes = "Significant on CL as a categorical covariate, but continuous CrCl was preferred because it reduced IIV on CL more (5.6% vs 3.9%) (Section 3.3)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 1L,
    n_observations = 158L,
    age_range = "23-81 years",
    age_median = "60 years",
    weight_range = "58-115 kg",
    weight_median = "80 kg",
    sex_female_pct = 33,
    race_ethnicity = "Not reported (Spanish ICU population)",
    disease_state = "Critically ill adults in the ICU treated with levetiracetam (haemorrhagic stroke 37%, trauma 30%, other neurological diagnoses 33%); APACHE II median 18 (range 5-35)",
    renal_function = "Measured urinary creatinine clearance median 117 mL/min (range 54-239); inclusion required CrCl > 50 mL/min; 10 of 27 (37%) had augmented renal clearance (CrCl > 130 mL/min)",
    dose_range = "500, 1000 or 1500 mg every 12 h as a 30-min IV infusion (18 of 27 patients on 500 mg q12h); sampled at steady state",
    regions = "Spain (Araba University Hospital, Vitoria-Gasteiz; Doce de Octubre Hospital, Madrid)",
    notes = "Baseline demographics per Bilbao-Meseguer 2021 Table 1. Prospective open-label two-centre study, 2019-2020. Median 6 (minimum 5) plasma samples per patient; HPLC-UV assay linear 2-100 mg/L. NONMEM 7.4 FOCE-I."
  )

  ini({
    # Structural parameters: Bilbao-Meseguer 2021 Table 3, 'Final Model Estimate'
    # column, and the final-model equation in Section 3.3.
    lcl_nonren <- log(3.5); label("Non-renal (CrCl-independent) arm of CL (L/h)") # Table 3: theta_nr = 3.5 (RSE 9%)
    e_crcl_cl_renal <- 2.5; label("Exponent of (CRCL/120) giving the renal arm of CL in L/h (unitless)") # Table 3: theta_r = 2.5 (RSE 17%); Section 3.3 equation CL = (3.5 + (CrCl/120)^2.5) x exp(eta1)
    lvc <- log(20.7); label("Central volume of distribution V1 (L)") # Table 3: V1 = 20.7 L (RSE 18%)
    lq <- log(31.9); label("Intercompartmental clearance Q (L/h)") # Table 3: Q = 31.9 L/h (RSE 22%)
    lvp <- log(33.5); label("Peripheral volume of distribution V2 (L)") # Table 3: V2 = 33.5 L (RSE 13%)

    # Inter-individual variability: exponential (log-normal) on CL and V1,
    # diagonal (Section 3.3, 'no correlation was detected'). omega^2 =
    # log(CV^2 + 1) from the Table 3 IIV percentages.
    etalcl ~ 0.10155 # Table 3: IIV_CL = 32.7% (RSE 21%); log(0.327^2 + 1)
    etalvc ~ 0.27293 # Table 3: IIV_V1 = 56.1% (RSE 29%); log(0.561^2 + 1)

    # Proportional residual error: Table 3 'RE_proportional' = 22.3% (RSE 15%).
    propSd <- 0.223; label("Proportional residual error (fraction)") # Table 3: RE_proportional = 22.3%
  })
  model({
    # Section 3.3 final-model equation: CL = (theta_nr + (CrCl/120)^theta_r) x exp(eta1).
    # The renal arm has an implicit coefficient of 1 L/h at CRCL = 120 mL/min.
    cl_nonren <- exp(lcl_nonren)
    cl_renal <- (CRCL / 120)^e_crcl_cl_renal
    cl <- (cl_nonren + cl_renal) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
