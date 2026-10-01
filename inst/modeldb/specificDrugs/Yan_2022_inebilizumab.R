Yan_2022_inebilizumab <- function() {
  description <- "Two-compartment population PK model for inebilizumab (anti-CD19 afucosylated IgG1k) in adults with neuromyelitis optica spectrum disorder, systemic sclerosis or relapsing multiple sclerosis, with parallel linear and Michaelis-Menten eliminations from the central compartment. The Michaelis-Menten Vmax (CD19-mediated clearance) decays mono-exponentially with time since first dose, reflecting B-cell depletion, and is higher in the systemic-sclerosis study MI-CP200. Body weight scales CL, Vc, Q and Vp by estimated power exponents. A first-order subcutaneous depot with the separately reported absorption half-life and bioavailability is included."
  reference <- "Yan L, Kimko H, Wang B, Cimbora D, Katz E, Rees WA. Population Pharmacokinetic Modeling of Inebilizumab in Subjects with Neuromyelitis Optica Spectrum Disorders, Systemic Sclerosis, or Relapsing Multiple Sclerosis. Clin Pharmacokinet. 2022;61(3):387-400. doi:10.1007/s40262-021-01071-5"
  vignette <- "Yan_2022_inebilizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on CL (exponent 0.57), Vc (0.39), Q (0.84) and Vp (0.40), normalised to the analysis-dataset median of 66.2 kg (Yan 2022 Eq. 3 and Table 2 'Total' median; the FDA clinical pharmacology review of BLA 761142 states the same 66.2 kg reference). Included in the structural model a priori (Yan 2022 Section 3.2).",
      source_name = "WT"
    ),
    STUDY_MICP200 = list(
      description = "Study indicator: 1 = subject enrolled in study MI-CP200 (single ascending IV doses in systemic sclerosis), 0 = subject in study CD-IA-MEDI-551-1102 (relapsing MS) or CD-IA-MEDI-551-1155 (NMOSD)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (studies 1102 and 1155, the most common realisation per Yan 2022 Eq. 4)",
      notes = "Fractional-change effect on Vmax per Yan 2022 Eq. 4: Vmax = Vmax_ref * (1 + 2.10 * STUDY_MICP200), with Table 6 'Study CP200 on Vmax (%)' = 210. The study effect is confounded with disease (all MI-CP200 subjects had systemic sclerosis); the paper labels it a study effect.",
      source_name = "Study CP200"
    )
  )

  covariatesDataExcluded <- list(
    ADA_POS = list(
      description = "Time-varying anti-drug antibody status (titer >= 50)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ADA-negative)",
      notes = "Tested on CL (Yan 2022 Table 4 model 6); the OFV increased by 0.466, so the effect was not retained."
    ),
    CRCL = list(
      description = "Renal function: estimated glomerular filtration rate",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Neither eGFR (mL/min/1.73 m^2) nor unnormalised creatinine clearance (mL/min) correlated with CL (Yan 2022 Section 3.3); not retained."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "No correlation with CL (Yan 2022 Section 3.3); not retained. ALP and AST likewise screened and not retained."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "No significant effect on CL (Yan 2022 Section 3.3); not retained. Sex, race, AQP4-IgG serostatus, EDSS, prior NMOSD attacks, disease duration, baseline CD20+ B-cell count and concomitant paracetamol / diphenhydramine / prednisolone / methylprednisolone were likewise screened and not retained."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "inebilizumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "inebilizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "inebilizumab", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 213L,
    n_studies = 3L,
    age_range = "18-73 years",
    age_median = "44 years",
    weight_range = "38.0-148 kg",
    weight_median = "66.2 kg",
    sex_female_pct = 86.9,
    race_ethnicity = c(White = 58.7, Black = 8.9, Asian = 18.3, American_Indian_or_Alaska_Native = 6.6, Other = 7.5),
    disease_state = "Neuromyelitis optica spectrum disorder (study 1155, n = 174), systemic sclerosis (study MI-CP200, n = 24) and relapsing-remitting multiple sclerosis (study 1102, n = 15).",
    dose_range = "Single IV 0.1-10 mg/kg (MI-CP200); two IV infusions of 30, 100 or 600 mg on days 1 and 15 (1102); two IV infusions of 300 mg on days 1 and 15 (1155). Six MS subjects given single SC 60 or 300 mg were excluded from the IV model and used only to characterise SC absorption.",
    regions = "Multinational",
    notes = "Baseline demographics from Yan 2022 Table 2 (IV analysis set, N = 213; 1617 concentrations analysed). Baseline creatinine clearance median 110 mL/min; ADA-positive 9.9%."
  )

  ini({
    # Structural parameters for the reference subject (WT 66.2 kg, studies
    # 1102/1155). Yan 2022 Table 6 reports CL/Q in mL/day and volumes in mL;
    # converted to L/day and L here so that Cc = central / vc is in mg/L = ug/mL.
    lcl <- log(0.188); label("Linear clearance CL (L/day)") # Yan 2022 Table 6: CL 188 mL/day (RSE 2.2%)
    lvc <- log(2.95); label("Central volume of distribution Vc (L)") # Yan 2022 Table 6: Vc 2950 mL (RSE 1.4%)
    lq <- log(0.363); label("Intercompartmental clearance Q (L/day)") # Yan 2022 Table 6: Q 363 mL/day (RSE 6.0%)
    lvp <- log(2.57); label("Peripheral volume of distribution Vp (L)") # Yan 2022 Table 6: Vp 2570 mL (RSE 2.8%)

    # Time-dependent Michaelis-Menten elimination (Yan 2022 Eqs. 6 and 8):
    # TDVM = VMAX * exp(-Kdec * time). Vmax converted from ug/day to mg/day.
    lvmax <- log(0.832); label("Baseline Michaelis-Menten Vmax, studies 1102/1155 (mg/day)") # Yan 2022 Table 6: Vmax 832 ug/day (RSE 5.3%)
    lkdes <- log(0.00294); label("First-order rate of decrease of Vmax with time, Kdec (1/day)") # Yan 2022 Table 6: Kdec 0.00294 /day (RSE 55.1%)
    lkm <- log(5.89); label("Michaelis-Menten Km (ug/mL)") # Yan 2022 Table 6: Km 5.89 ug/mL (RSE 25.5%)

    # Subcutaneous absorption (Yan 2022 Section 3.4, six MS subjects combined
    # with the IV data after the IV model was final). The paper reports the
    # absorption half-life and absolute bioavailability only.
    lka <- log(log(2) / 4.1); label("First-order SC absorption rate constant ka (1/day)") # Yan 2022 Section 3.4: absorption half-life 4.1 days; ka = ln(2)/4.1 = 0.169 /day
    lfdepot <- log(0.81); label("Absolute SC bioavailability (fraction)") # Yan 2022 Section 3.4: absolute bioavailability 81%

    # Body-weight power exponents (Yan 2022 Eq. 3, weight / 66.2 kg).
    e_wt_cl <- 0.57; label("Power exponent of WT/66.2 on CL (unitless)") # Yan 2022 Table 6: Weight on CL 0.57 (RSE 15.8%)
    e_wt_vc <- 0.39; label("Power exponent of WT/66.2 on Vc (unitless)") # Yan 2022 Table 6: Weight on Vc 0.39 (RSE 22.4%)
    e_wt_q <- 0.84; label("Power exponent of WT/66.2 on Q (unitless)") # Yan 2022 Table 6: Weight on Q 0.84 (RSE 21.1%)
    e_wt_vp <- 0.40; label("Power exponent of WT/66.2 on Vp (unitless)") # Yan 2022 Table 6: Weight on Vp 0.40 (RSE 27.9%)

    # Study MI-CP200 (systemic sclerosis) fractional change on Vmax (Eq. 4).
    e_study_micp200_vmax <- 2.10; label("Fractional change in Vmax for study MI-CP200, applied as 1 + e * STUDY_MICP200 (fraction)") # Yan 2022 Table 6: Study CP200 on Vmax 210% (RSE 19.5%)

    # IIV reported as %CV; omega^2 = log(CV^2 + 1).
    # CL 27% -> 0.07037; Vc 17% -> 0.02849; Vp 16% -> 0.02528; Vmax 30% -> 0.08618.
    # No IIV on Q, Km, Kdec, ka or F. CL-Vmax covariance was not retained
    # (Yan 2022 Table 4 model 7, dOFV -0.049).
    etalcl ~ 0.07037 # Yan 2022 Table 6: IIV CL 27% CV
    etalvc ~ 0.02849 # Yan 2022 Table 6: IIV Vc 17% CV
    etalvp ~ 0.02528 # Yan 2022 Table 6: IIV Vp 16% CV
    etalvmax ~ 0.08618 # Yan 2022 Table 6: IIV Vmax 30% CV

    # Residual error: proportional only in the final model (Yan 2022 Table 6
    # and Section 3.3 'the proportional residual error model in the final
    # pharmacokinetic model was sufficient').
    propSd <- 0.218
    label("Proportional residual error (fraction)") # Yan 2022 Table 6: Proportional error 21.8% CV (RSE 4.6%)
  })
  model({
    # Individual parameters (Yan 2022 Eqs. 1, 3 and 4).
    cl <- exp(lcl + etalcl) * (WT / 66.2)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 66.2)^e_wt_vc
    q <- exp(lq) * (WT / 66.2)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 66.2)^e_wt_vp
    vmax <- exp(lvmax + etalvmax) * (1 + e_study_micp200_vmax * STUDY_MICP200)
    kdes <- exp(lkdes)
    km <- exp(lkm)
    ka <- exp(lka)

    # Time-dependent Vmax (Yan 2022 Eq. 8); t is time since the first dose.
    vmax_t <- vmax * exp(-kdes * t)

    Cc <- central / vc

    # Yan 2022 Eqs. 5-7, with a first-order SC depot (Section 3.4).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot -
      (cl / vc) * central -
      (q / vc) * central +
      (q / vp) * peripheral1 -
      vmax_t * Cc / (km + Cc)
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1

    f(depot) <- exp(lfdepot)

    Cc ~ prop(propSd)
  })
}
