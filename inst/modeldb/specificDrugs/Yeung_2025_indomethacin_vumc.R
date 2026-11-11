Yeung_2025_indomethacin_vumc <- function() {
  description <- paste(
    "One-compartment population PK model for intravenous indomethacin in",
    "extremely preterm neonates given a single 0.2 mg/kg prophylactic dose",
    "at 12 h of life (CYP2C9-reduced IND popPK model, re-estimated on the",
    "Vanderbilt University Medical Center dried-blood-spot dataset).",
    "Volume, non-renal (hepatic) clearance and its postnatal-age exponent",
    "were re-estimated; the renal-clearance, plasma and urine terms were",
    "held fixed because the cohort had DBS samples only. Same structure as",
    "Yeung_2025_indomethacin_cmh.",
    sep = " "
  )
  reference <- paste(
    "Yeung CHT, Sekulich DC, Scott A, Nolte WM, Gibson K, Su R,",
    "Alrifai MW, Lopata SM, Lewis T, Reese J, Maharaj A (2025).",
    "The Relationship of Indomethacin Exposure With Efficacy and Renal",
    "Toxicity Outcomes for Preterm Infants in the Neonatal Intensive Care",
    "Unit. Clin Transl Sci 18(5):e70251. doi:10.1111/cts.70251.",
    "Structural model from Krzyzanski W, Stockard B, Gaedigk A, et al.",
    "(2023) Developmental pharmacokinetics of indomethacin in preterm",
    "neonates: severely decreased drug clearance in the first week of",
    "life. CPT Pharmacometrics Syst Pharmacol 12:110-121.",
    "doi:10.1002/psp4.12881.",
    sep = " "
  )
  vignette <- "Yeung_2025_indomethacin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    PNA = list(
      description = "Postnatal age, time-varying (updated at every record).",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power terms in DAYS centred on 1.25 days, as in",
        "Yeung_2025_indomethacin_cmh (Krzyzanski 2023 Equation 5); Yeung",
        "2025 states the structural model was retained for the VUMC",
        "re-estimation, so the centring value is carried over. model()",
        "recovers days with PNA * 30.4375. Supply PNA > 0."
      ),
      source_name = "PNA (days)"
    ),
    SAMPLE_DBS = list(
      description = "Per-observation sampling-matrix indicator (1 = concentration from a dried blood spot converted to plasma; 0 = plasma).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (plasma)",
      notes = paste(
        "Every VUMC concentration was a DBS converted to plasma with Yeung",
        "2025 Equation 1, so SAMPLE_DBS = 1 reproduces the fitted residual",
        "error. The plasma residual is carried fixed from the CMH",
        "re-estimate for users who simulate plasma samples."
      ),
      source_name = "DBS"
    )
  )

  compartmentData <- list(
    central = list(analyte = "indomethacin", units = "mg", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "indomethacin", units = "mg", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 30,
    n_studies = 1,
    age_range = "postnatal age at the dose 11.4-20.3 h (median 13 h)",
    weight_range = "370-1040 g at first dose (median 645 g)",
    ga_range = "22.1-27.9 weeks gestational age at birth (median 24.9)",
    sex_female_pct = 53,
    race_ethnicity = c(Black = 36.7, White = 56.7, Other = 3.3, `Not reported` = 3.3),
    disease_state = "Extremely preterm neonates (< 28 weeks GA, < 1000 g birthweight) receiving indomethacin for IVH and PDA prophylaxis",
    dose_range = "0.2 mg/kg IV single dose at 12 h of life",
    regions = "USA (Vanderbilt University Medical Center, Nashville)",
    notes = paste(
      "Prospective cohort recruited April 2021 to March 2023 (Yeung 2025",
      "Methods 'Study Populations', Table 1 VUMC IVH column). On average",
      "8.5 dried-blood-spot concentrations per infant over the first week",
      "of life, median 46 (1.03-209.2) h after the dose."
    )
  )

  ini({
    lvc <- log(0.164); label("Volume of distribution (L)") # Table S2, VUMC 'theta V [L]' = 0.164 (RSE 6.2%)
    lcl_renal <- fixed(log(0.0000172)); label("Renal clearance at PNA = 1.25 days (L/h)") # Table S2, VUMC 'theta CLR' = NA (renal terms held fixed, Table S2 footnote); value from the CMH re-estimate 0.0000172
    lcl_nonren <- log(0.000572); label("Non-renal (hepatic) clearance at PNA = 1.25 days (L/h)") # Table S2, VUMC 'theta CLH [L/h]' = 0.000572 (RSE 22.4%)
    e_pna_cl_renal <- fixed(-0.192); label("Power exponent of PNA on renal clearance (unitless)") # Table S2, VUMC 'theta PNA-CLR' = NA (held fixed); value from the CMH re-estimate -0.192
    e_pna_cl_nonren <- 1.21; label("Power exponent of PNA on non-renal clearance (unitless)") # Table S2, VUMC 'theta PNA-CLH' = 1.21 (RSE 11.5%)

    etalvc ~ log(1 + 0.263^2) # Table S2, VUMC 'IIV V (%)' = 26.3
    etalcl_renal ~ fixed(log(1 + 0.690^2)) # Table S2, VUMC 'IIV CLR (%)' = NA (not re-estimated); CMH re-estimate 69.0
    etalcl_nonren ~ log(1 + 0.698^2) # Table S2, VUMC 'IIV CLH (%)' = 69.8

    expSd <- fixed(0.312730); label("Log-scale residual SD, plasma concentrations (unitless)") # Table S2, VUMC 'epsilon2 plasma' = NA (no plasma data); CMH re-estimate sqrt(0.0978)
    expSdDbs <- 0.425441; label("Log-scale residual SD, dried blood spot concentrations (unitless)") # Table S2, VUMC 'epsilon2 log-normal DBS' = 0.181 (RSE 23.9%); sqrt(0.181)
    expSd_Aurine <- fixed(0.310483); label("Log-scale residual SD, cumulative urine amount (unitless)") # Table S2, VUMC 'epsilon2 urine' = NA (no urine data); CMH re-estimate sqrt(0.0964)
  })

  model({
    pna_days <- PNA * 30.4375

    vc <- exp(lvc + etalvc)
    cl_renal <- exp(lcl_renal + etalcl_renal) * (pna_days / 1.25)^e_pna_cl_renal
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) * (pna_days / 1.25)^e_pna_cl_nonren

    d/dt(central) <- -(cl_renal + cl_nonren) / vc * central
    d/dt(urine) <- cl_renal / vc * central

    Cc <- central / vc
    Aurine <- urine

    expSdCc <- expSd * (1 - SAMPLE_DBS) + expSdDbs * SAMPLE_DBS
    Cc ~ lnorm(expSdCc)
    Aurine ~ lnorm(expSd_Aurine)
  })
}
