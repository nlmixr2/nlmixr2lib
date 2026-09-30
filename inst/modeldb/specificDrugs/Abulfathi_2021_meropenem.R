Abulfathi_2021_meropenem <- function() {
  description <- "Two-compartment population PK model for intravenous meropenem in South African adults with rifampicin-sensitive pulmonary tuberculosis (COMRADE trial; Abulfathi 2021). Allometric scaling on all disposition parameters with total body weight (reference 70 kg; fixed exponents 0.75 on CL and Q, 1 on V1 and V2) and a power effect of weight-standardised Cockcroft-Gault creatinine clearance (CRCL * 70 / WT, reference 115 mL/min) on CL; combined additive + proportional residual error."
  reference <- "Abulfathi AA, de Jager V, van Brakel E, Reuter H, Gupte N, Vanker N, Barnes GL, Nuermberger E, Dorman SE, Diacon AH, Dooley KE, Svensson EM. The Population Pharmacokinetics of Meropenem in Adult Patients With Rifampicin-Sensitive Pulmonary Tuberculosis. Front Pharmacol. 2021;12:637618. doi:10.3389/fphar.2021.637618"
  vignette <- "Abulfathi_2021_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric size descriptor on all disposition parameters, reference",
        "70 kg with fixed theoretical exponents 0.75 (CL, Q) and 1 (V1, V2)",
        "(Abulfathi 2021 Results and Table 2 footnote: TVCL = THETA(1) *",
        "((WTKG/70)**0.75) * ...; TVV1 = THETA(2) * WTKG/70; TVQ = THETA(3) *",
        "((WTKG/70)**0.75); TVV2 = THETA(4) * WTKG/70). WT also standardises",
        "the creatinine clearance to a 70 kg body size (CRCL * 70 / WT)",
        "inside the CL covariate term. Observed range 39.3-76.3 kg, median",
        "52.7 kg (Table 1)."
      ),
      source_name = "WTKG"
    ),
    CRCL = list(
      description = paste(
        "Creatinine clearance by the Cockcroft-Gault equation (raw mL/min,",
        "NOT BSA-normalized)."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Carried as the raw Cockcroft-Gault value in mL/min, not the canonical",
        "BSA-normalized mL/min/1.73 m^2 form (the same deviation documented in",
        "Tsuji_2017_linezolid.R). Inside the model the value is first",
        "size-standardised to a 70 kg body weight, CRCL * 70 / WT, and then",
        "enters CL as a power function normalised to the population median",
        "115 mL/min (Abulfathi 2021 Table 2 footnote: ((CLCR*70/WTKG)/115) **",
        "THETA(7), with THETA(7) = 0.416). Observed range 57.7-203 mL/min,",
        "median 115 mL/min (Table 1)."
      ),
      source_name = "CLCR"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on CL before and during the stepwise covariate search (a previously reported effect) but not retained in the final model (Abulfathi 2021 Results: 'both age and rifampicin had an insignificant impact on meropenem CL')."
    ),
    CONMED_RIFAMPICIN = list(
      description = "Concomitant oral rifampicin (20 mg/kg once daily)",
      units = "(binary)",
      type = "binary",
      notes = "The MACR2X3 arm received rifampicin; its effect on CL was tested and found insignificant (Abulfathi 2021 Results and Discussion), so it is not in the final model."
    ),
    HT = list(
      description = "Height",
      units = "m",
      type = "continuous",
      notes = "Screened on CL in the stepwise covariate search; not retained (Abulfathi 2021 Methods 'Covariate Model' and Results)."
    ),
    HIV_POS = list(
      description = "HIV-positive status",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL in the stepwise covariate search; not retained (Abulfathi 2021 Methods 'Covariate Model' and Results)."
    ),
    SEXF = list(
      description = "Female sex",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL and V1 in the stepwise covariate search; not retained (Abulfathi 2021 Methods 'Covariate Model' and Results)."
    ),
    RACE_BLACK = list(
      description = "Black race (versus mixed Asian ancestry)",
      units = "(binary)",
      type = "binary",
      notes = "Race was screened on CL and V1 in the stepwise covariate search; not retained (Abulfathi 2021 Methods 'Covariate Model' and Results)."
    ),
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      notes = "Tested as the allometric size descriptor in place of total body weight; it gave a smaller OFV drop (20.4 vs 22.2 points) and total body weight was kept (Abulfathi 2021 Discussion)."
    )
  )

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 49L,
    n_studies = 1L,
    age_range = "20.0-62.7 years",
    age_median = "36.0 years",
    weight_range = "39.3-76.3 kg",
    weight_median = "52.7 kg",
    sex_female_pct = 24.5,
    race_ethnicity = c(Black = 32.7, `Mixed Asian ancestry` = 67.3),
    disease_state = "Sputum smear-positive, rifampicin-sensitive pulmonary tuberculosis (22.4% HIV-positive); creatinine clearance 57.7-203 mL/min (median 115)",
    dose_range = paste(
      "Meropenem IV for 14 days in four arms: 2 g over 0.5 h every 8 h plus",
      "oral rifampicin 20 mg/kg once daily (MACR2X3, n = 12); 2 g over 0.5 h",
      "every 8 h (MAC2X3, n = 13); 1 g over 0.5 h every 8 h (MAC1X3, n = 12);",
      "3 g over 1 h once daily (MAC3X1, n = 12). All arms also received oral",
      "amoxicillin/clavulanate with each meropenem dose."
    ),
    regions = "South Africa",
    notes = paste(
      "Demographics from Abulfathi 2021 Table 1 (n = 49 with PK sampling out",
      "of 60 randomised). Phase 2 open-label randomised COMRADE trial.",
      "Intensive sampling on day 14 at pre-dose and 0.5, 1, 1.5, 2, 3, 4, 6",
      "and 8 h post-dose; 404 of 441 concentrations analysed (34 BQL, LLOQ",
      "0.5 mg/L; 3 outliers excluded)."
    )
  )

  ini({
    # Structural parameters for a 70 kg subject with size-standardised CrCl
    # 115 mL/min (Abulfathi 2021 Table 2, 'Population estimate' column).
    lcl <- log(11.8); label("Clearance at 70 kg and CRCL*70/WT = 115 mL/min (L/h)") # Table 2: CL (L/h/70 kg) = 11.8 (RSE 4.9%)
    lvc <- log(14.2); label("Central volume at 70 kg (L)") # Table 2: V1 (L/70 kg) = 14.2 (RSE 3.8%)
    lq <- log(3.26); label("Intercompartmental clearance at 70 kg (L/h)") # Table 2: Q (L/h/70 kg) = 3.26 (RSE 27.5%)
    lvp <- log(3.12); label("Peripheral volume at 70 kg (L)") # Table 2: V2 (L/70 kg) = 3.12 (RSE 10.8%)

    # Fixed theoretical allometric exponents (Results: 'fixed theoretical
    # exponents of 1 for volume of distribution and 0.75 for clearance';
    # Table 2 footnote TVCL/TVQ ** 0.75, TVV1/TVV2 linear in WTKG/70).
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)") # Results / Table 2 footnote
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on V1 and V2 (unitless)") # Results / Table 2 footnote

    # Power exponent of size-standardised CrCl on CL (Table 2 footnote THETA(7)).
    e_crcl_cl <- 0.416; label("Power exponent of (CRCL*70/WT)/115 on CL (unitless)") # Table 2: 'Creatinine clearance on CL' = 0.416 (RSE 30.5%)

    # Inter-individual variability. Table 2 footnote b: %CV = SQRT(EXP(OMEGA)-1)*100,
    # so OMEGA = log(CV^2 + 1). No IIV on Q (Results: 'No significant
    # variability could be detected in Q'); no off-diagonal terms reported.
    etalcl ~ 0.039221 # log(0.20^2 + 1); Table 2 'IIV CL' 20 %CV
    etalvc ~ 0.017016 # log(0.131^2 + 1); Table 2 'IIV V1' 13.1 %CV
    etalvp ~ 0.753120 # log(1.06^2 + 1); Table 2 'IIV V2' 106 %CV

    # Combined additive + proportional residual error (Results; Table 2).
    propSd <- 0.178; label("Proportional residual error (fraction)") # Table 2: 'Proportional residual error (%)' = 0.178 (RSE 14.8%)
    addSd <- 1.16; label("Additive residual error (mg/L)") # Table 2: 'Additive residual error (mg/L)' = 1.16 (RSE 19.6%)
  })
  model({
    # Size-standardised creatinine clearance (mL/min per 70 kg), Table 2 footnote.
    crcl_std <- CRCL * 70 / WT

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q * (crcl_std / 115)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> mg/L (Figure 2 schema: concentration A1/V1).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
