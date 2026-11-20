Briquet_2026_meropenem <- function() {
  description <- paste(
    "Two-compartment population PK model for UNBOUND meropenem in 22",
    "non-critically ill adult inpatients with normal or mildly reduced renal",
    "function (eGFR 64-234 mL/min/1.73 m2) treated with 2 g three times daily",
    "as 20 min intravenous infusions, first every 8 h and then on a daytime-only",
    "q6-6-12 h schedule designed for outpatient parenteral antimicrobial therapy",
    "(OPAT) (Briquet 2026). Concentrations are unbound meropenem measured in",
    "plasma ultrafiltrate, so Cc is the unbound concentration and the volumes",
    "are apparent volumes referenced to it. Linear elimination from the central",
    "compartment; log-normal IIV on clearance and central volume with a",
    "correlation of 0.76; combined additive plus proportional residual error. No",
    "covariate was retained (body weight, height, sex, age, serum creatinine and",
    "eGFR were screened). Fitted with SAEM in Monolix 2024R1."
  )
  reference <- paste(
    "Briquet C, Ngougni Pokem P, Liu X, Wijnant GJ, Cornu O, Yildiz H,",
    "Verroken A, Roberts JA, Yombi JC, Van Bambeke F.",
    "A novel meropenem dosing strategy for Outpatient Parenteral Antimicrobial",
    "Therapy (OPAT): a pilot pharmacokinetic/pharmacodynamic study in",
    "non-critically ill patients.",
    "Antimicrob Agents Chemother. 2026;70(7):e00372-26.",
    "doi:10.1128/aac.00372-26"
  )
  vignette <- "Briquet_2026_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened (Methods, 'Population pharmacokinetic analysis'); no covariate met",
        "the inclusion criteria. Table 1 median 85 kg (range 44-155)."
      )
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened; not retained. Table 1 does not report height."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "unitless",
      type = "categorical",
      notes = "Screened ('gender'); not retained. Table 1: 12 male, 10 female."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened; not retained. Table 1 median 61.5 years (range 26-82)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened; not retained. Table 1 does not report serum creatinine."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (CKD-EPI)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = paste(
        "Screened; not retained. Table 1 median 102 (range 64-234, mean 108.1);",
        "16 patients above 90 and 6 between 60 and 90. Inclusion required an eGFR",
        "of at least 50, so the narrow renal-function spread is the likely reason",
        "no renal covariate was found (Discussion)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 22L,
    n_studies = 1L,
    n_observations = 242L,
    age_range = "26-82 years (median 61.5)",
    age_median = "61.5 years",
    weight_range = "44-155 kg (median 85)",
    weight_median = "85 kg",
    sex_female_pct = 45.5,
    disease_state = paste(
      "Non-critically ill, clinically stable hospitalised adults with",
      "Gram-negative infections (meropenem MIC <= 2 mg/L): bone 15 (68%),",
      "intra-abdominal 5 (23%), skin and soft tissue 2 (9%)."
    ),
    renal_function = paste(
      "CKD-EPI eGFR median 102 mL/min/1.73 m2 (range 64-234); 16 patients",
      "above 90 and 6 with mild impairment (60-90). eGFR >= 50 was an",
      "inclusion criterion."
    ),
    dose_range = paste(
      "Meropenem 2 g as a 20 min intravenous infusion three times daily, given",
      "each patient first every 8 h (q8-8-8 h) and then, once clinically stable,",
      "at 6 h intervals during the day only (q6-6-12 h)."
    ),
    regions = "Belgium (Cliniques universitaires Saint-Luc, Brussels)",
    notes = paste(
      "Single-centre prospective study (EudraCT 2019-003900-10), November 2019",
      "to October 2021. 246 samples were collected (all 22 patients on the",
      "conventional regimen, 20 on the OPAT-adapted regimen); one sample with an",
      "uncertain time and three Tukey 3.5 x IQR outliers were excluded.",
      "Unbound concentrations by HPLC-MS/MS on serum ultrafiltrate. Baseline",
      "characteristics are Table 1."
    )
  )

  ini({
    # Structural parameters: Table 2 'Estimate (%RSE)' column.
    lcl <- log(11.7); label("Clearance of unbound meropenem (L/h)") # Table 2 'CL (L/h)' 11.7 (RSE 11.1%)
    lvc <- log(23.6); label("Central volume of distribution (L)") # Table 2 'V1 (L)' 23.6 (RSE 11.9%)
    lq <- log(1.98); label("Inter-compartmental clearance (L/h)") # Table 2 'Q (L/h)' 1.98 (RSE 39.5%)
    lvp <- log(5.44); label("Peripheral volume of distribution (L)") # Table 2 'V2 (L)' 5.44 (RSE 17.3%)

    # IIV: Table 2 'IIV_CL (%)' 51.7 and 'IIV_V1 (%)' 42.9 are Monolix omega
    # values (SD of the log-scale random effect) x 100, so the variances are
    # 0.517^2 = 0.267289 and 0.429^2 = 0.184041. 'Corr_V1_Cl' 0.76 is a
    # correlation coefficient: cov = 0.76 * 0.517 * 0.429 = 0.16856268.
    etalcl + etalvc ~ c(0.267289, 0.16856268, 0.184041)

    # Residual error: Table 2 'Additive (mg/L)' 0.16 (RSE 48.5%) and
    # 'Proportional' 0.41 (RSE 8.26%), Monolix combined error.
    addSd <- 0.16; label("Additive residual error (mg/L)")
    propSd <- 0.41; label("Proportional residual error (fraction)")
  })

  model({
    # Exponential IIV, theta_j = theta_p * exp(eta_j) (Methods equation).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment model with first-order elimination (Results).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Unbound concentration: dose in mg over volumes in L gives mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
