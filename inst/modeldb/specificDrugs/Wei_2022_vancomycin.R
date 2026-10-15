Wei_2022_vancomycin <- function() {
  description <- "One-compartment IV-infusion population PK model for vancomycin in Chinese adults after neurosurgery (Wei 2022). Clearance scales with a power function of CKD-EPI estimated glomerular filtration rate (exponent 0.80, reference 115.2 mL/min) and of body weight (exponent 0.30, reference 70 kg) and is multiplied by exp(0.13) under concomitant mannitol; the volume of distribution is fixed at 60.2 L from an earlier Chinese neurosurgical model because the trough-dominated data could not estimate it. Exponential IIV on clearance only and combined additive plus proportional residual error."
  reference <- "Wei S, Zhang D, Zhao Z, Mei S. Population pharmacokinetic model of vancomycin in postoperative neurosurgical patients. Front Pharmacol. 2022;13:1005791. doi:10.3389/fphar.2022.1005791"
  vignette <- "Wei_2022_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Estimated glomerular filtration rate by the creatinine-based CKD-EPI equation, as printed by Wei 2022 (labelled mL/min, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Wei 2022 Eq. (4): eGFR (mL/min) = 144 * (Scr/a)^b * 0.993^age, with Scr in mg/dL; female a = 0.7, b = -0.329 (Scr <= 0.7) or -1.209 (Scr > 0.7); male a = 0.9, b = -0.411 (Scr <= 0.9) or -1.210 (Scr > 0.9). The paper prints the 144 multiplier for both sexes and labels the result mL/min, with no 1.73 m^2 term, so the column is supplied exactly as that equation computes it. Table 1: mean 112.74 (SD 30.91), range 3.52-244.48 mL/min. Reference 115.2 mL/min in Eq. (5) (Methods: continuous covariates centered at their median). Stored under the canonical CRCL, whose register entry lists `eGFR` as a source alias; the Cockcroft-Gault CLcr that the paper also screened (Table 2 model 3) is a different column and is not used.",
      source_name = "eGFR"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Wei 2022 Table 1: mean 69.74 (SD 13.05), range 37.5-130 kg. Reference 70 kg in Eq. (5).",
      source_name = "BW"
    ),
    CONMED_MANNITOL = list(
      description = "Concomitant mannitol indicator (1 = co-medicated with mannitol)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant mannitol)",
      notes = "Wei 2022 Eq. (5): CL multiplied by e^A, A = 0.13 when co-medicated with mannitol, otherwise A = 0. Table 1: 60.32% of patients used mannitol (to relieve cerebral oedema and reduce intracranial pressure). The paper does not say whether the flag was set per record or per patient; Supplementary Appendix S2 splits the VPC by patient-level mannitol use.",
      source_name = "Mannitol"
    )
  )

  # Covariates screened on CL (Table 2, forward addition from the base model)
  # but not retained in the final model. No final-model estimate exists for any.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age. Screened on CL and not retained.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Wei 2022 Table 2 model 6 (dOFV -77.30 against the base model); not carried forward once eGFR, which already contains age, entered. Table 1: mean 52.41 (SD 15.11), range 18-89 years."
    ),
    HT = list(
      description = "Height. Screened on CL and not retained.",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Wei 2022 Table 2 model 9 (dOFV -9.08). Table 1: mean 167.88 (SD 7.98), range 145-192 cm."
    ),
    BMI = list(
      description = "Body mass index. Screened on CL and not retained.",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Wei 2022 Table 2 model 11 (dOFV -4.87, p > 0.01). Table 1: mean 24.64 (SD 3.64), range 15.61-47.75."
    ),
    CREAT = list(
      description = "Serum creatinine. Screened on CL and not retained in favour of the CKD-EPI eGFR.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Wei 2022 Table 2 model 4 (dOFV -396.75); Discussion: the final Scr model (OFV 5536.97) fit worse than the eGFR model (OFV 5464.67). Table 1: mean 64.87 (SD 76.89), range 9.79-957.5 umol/L."
    ),
    SEXF = list(
      description = "Female sex indicator. Screened on CL and not retained.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Wei 2022 Table 2 model 13 (dOFV -0.05). Table 1: 370 male / 190 female."
    ),
    CONMED_MEROPENEM = list(
      description = "Concomitant meropenem indicator. Screened on CL and not retained.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant meropenem)",
      notes = "Wei 2022 Table 2 model 12 (dOFV -0.12). Table 1: 71.32% of patients."
    ),
    CONMED_DIURETIC = list(
      description = "Concomitant diuretic indicator (furosemide or torasemide). Screened on CL and not retained.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant diuretic)",
      notes = "Wei 2022 Methods 2.1 lists 'diuretics such as furosemide or torasemide'; Table 2 model 10 (dOFV -6.44, p > 0.01). Table 1: 15.95% of patients."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 560L,
    n_studies = 1L,
    n_concentrations = 895L,
    age_range = "18-89 years",
    age_mean = "52.41 years (SD 15.11)",
    weight_range = "37.5-130 kg",
    weight_mean = "69.74 kg (SD 13.05)",
    sex_female_pct = 33.9,
    race_ethnicity = "Chinese (single centre, Beijing).",
    disease_state = "Adults treated with intravenous vancomycin with at least one therapeutic-drug-monitoring concentration; 497 of the 560 patients were postoperative neurosurgical patients, most commonly with intracranial space-occupying lesions. Pregnant women and patients with cystic fibrosis were excluded. 68.2% (382/560) had augmented renal clearance during monitoring; 2.72% had severe renal insufficiency, most of them on renal replacement therapy.",
    dose_range = "0.5 or 1 g per dose every 6, 8, 12 or 24 h by intravenous infusion (default 1 h infusion assumed by the authors); dose per administration mean 951.19 mg (SD 152.23), range 50-1500 mg.",
    regions = "China (Beijing Tiantan Hospital), October 2018 to March 2022.",
    renal_function = "CKD-EPI eGFR mean 112.74 mL/min (SD 30.91, range 3.52-244.48); Cockcroft-Gault CLcr mean 152.94 mL/min (SD 74.89); serum creatinine mean 64.87 umol/L (SD 76.89). Renal replacement therapy in 2.72%.",
    co_medication = "Meropenem 71.32%, mannitol 60.32%, diuretics (furosemide or torasemide) 15.95%.",
    notes = "Demographics from Wei 2022 Table 1. Retrospective routine TDM data, mostly steady-state troughs (after more than 4 doses); serum concentrations 14.20 +/- 7.36 mg/L (range 0.91-52.96) by chemiluminescence immunoassay (ADVIA Centaur XP, calibration range 0.67-90 mg/L). Fitted in Phoenix NLME 8.3 with FOCE-ELS; evaluated by GOF plots, a 5000-replicate bootstrap and a 10000-replicate VPC (95.3% of observations inside the 90% prediction interval). No external validation."
  )

  ini({
    # Structural parameters: Wei 2022 Table 3, final-model 'Estimate (%RSE)'
    # column; covariate model Eq. (5)-(6):
    #   CL (L/h) = 7.98 * (eGFR/115.2)^0.8 * (BW/70)^0.3 * e^A, A = 0.13 with mannitol
    #   V (L)    = 60.2
    lcl <- log(7.98); label("Typical clearance at eGFR 115.2 mL/min, WT 70 kg, no mannitol (L/h)") # Table 3: CL = 7.98 L/h (RSE 1.90%, 95% CI 7.68-8.28); Eq. (5)
    # V was not estimable from trough-dominated data and was fixed to the value
    # of the Jing 2020 Chinese neurosurgical model (Methods 2.2.1); Table 3
    # gives it with no RSE or CI in every column.
    lvc <- fixed(log(60.2)); label("Volume of distribution (L)") # Table 3: V = 60.2 L, fixed; Eq. (6); Methods 2.2.1

    e_crcl_cl <- 0.80; label("Power exponent of (eGFR/115.2) on CL (unitless)") # Table 3 'eGFR on CL' = 0.80 (RSE 4.30%, 95% CI 0.74-0.87); Eq. (5)
    e_wt_cl <- 0.30; label("Power exponent of (WT/70) on CL (unitless)") # Table 3 'BW on CL' = 0.30 (RSE 20.19%, 95% CI 0.18-0.42); Eq. (5)
    e_conmed_mannitol_cl <- 0.13; label("Exponential effect of concomitant mannitol on CL (unitless)") # Table 3 'Mannitol on CL' = 0.13 (RSE 17.85%, 95% CI 0.08-0.17); Eq. (5) e^A, A = 0.13

    # IIV: exponential model (Eq. 1) on CL only. Table 3 reports 'IIV CL (CV%)'
    # = 21.45; omega^2 = log(1 + 0.2145^2). See the vignette for the check of
    # this scale against the Table 4 simulated AUC24 intervals.
    etalcl ~ 0.044983 # Table 3: IIV CL 21.45 CV% (final model); log(1 + 0.2145^2)

    # Residual error, Eq. (2): Cobs = Cpred + eps * sqrt(1 + (Cpred * sigma1 / sigma2)^2)
    # with eps ~ N(0, sigma2^2), i.e. SD = sqrt(sigma2^2 + (sigma1 * Cpred)^2).
    propSd <- 0.25; label("Proportional residual error (fraction)") # Table 3: sigma1 (multiplicative, CV) = 0.25 (RSE 6.45%, 95% CI 0.22-0.28)
    addSd <- 1.51; label("Additive residual error (mg/L)") # Table 3: sigma2 (additive) = 1.51 mg/L (RSE 25.83%, 95% CI 0.75-2.28)
  })

  model({
    cl <- exp(lcl + etalcl) * (CRCL / 115.2)^e_crcl_cl * (WT / 70)^e_wt_cl *
      exp(e_conmed_mannitol_cl * CONMED_MANNITOL)
    vc <- exp(lvc)

    kel <- cl / vc

    # One compartment, first-order elimination; doses are 1 h IV infusions
    # into central (Methods 2.1).
    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
