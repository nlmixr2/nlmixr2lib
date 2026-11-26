Gastine_2019_micafungin <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for micafungin in critically ill adult ICU",
    "patients, most with sepsis and two thirds on renal replacement therapy (Gastine 2019).",
    "Clearance decreases linearly with the daily SOFA score (-2.4% per point, centred on the",
    "cohort median of 12) and is 21.1% lower when total bilirubin exceeds 4 mg/dL; central",
    "volume is 30.8% lower when the SOFA score exceeds 10. Renal replacement therapy, body",
    "weight and the other screened covariates were not retained. Exponential IIV on CL, Vc and",
    "Vp; proportional residual error. NONMEM 7.3 ADVAN6 FOCE-I; 36 patients on 100 mg once",
    "daily. Sibling micafungin models: modellib('Martial_2017_micafungin') and",
    "modellib('Maseda_2018_micafungin').",
    sep = " "
  )
  reference <- paste(
    "Gastine S, Lanckohr C, Blessou M, Horn D, Fobker M, Bause D, Hempel G, Ellger B.",
    "Pharmacokinetics of micafungin in critically ill patients.",
    "Sci Rep. 2019;9:17741. doi:10.1038/s41598-019-53093-6.",
    "Parameter estimates from Table 2; covariate effects from Results 'Population",
    "pharmacokinetic modelling'; demographics from Table 1.",
    sep = " "
  )
  vignette <- "Gastine_2019_micafungin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds: total micafungin in plasma, measured by
  # reversed-phase HPLC (Methods 'Sample preparation and analysis of
  # micafungin', LLOQ 0.1 ug/mL); intravenous dosing into the central
  # compartment (Methods 'Study protocol').
  compartmentData <- list(
    central = list(analyte = "micafungin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "micafungin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    SOFA = list(
      description = "Sequential Organ Failure Assessment score, recorded daily (time-varying)",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Methods: SOFA was 'documented each day in the ICU' and 'All covariates were collected as",
        "time-varying covariates, if feasible'. Enters twice. (1) Linear on CL: Results report 'a",
        "linear decrease in clearance with a factor of 0.024 per score point' from the PsN",
        "stepwise covariate modelling; the PsN linear parameterisation (1 + theta * (SOFA -",
        "median)) is used, centred on the Table 1 cohort median of 12 (range 2-22). The centring",
        "is not printed; the uncentred reading (1 - 0.024 * SOFA) would put the typical clearance",
        "of a SOFA-12 patient at 1.11 L/h, against the 1.49 L/h implied by the steady-state",
        "AUC/MIC medians of Figure 7 (see the vignette). (2) Categorical on Vc: 'a 30.8% lower",
        "central volume of distribution in the critically ill group' with SOFA > 10; SOFA <= 10",
        "is the reference.",
        sep = " "
      ),
      source_name = "SOFA"
    ),
    TBILI = list(
      description = "Total serum bilirubin (time-varying)",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters only through a threshold: Results 'a categorical effect for bilirubin on",
        "clearance with a cut-off for bilirubin >4 mg/dl was included, showing a 21.1% lower",
        "clearance in the group with elevated bilirubin levels'. The paper works in mg/dL; the",
        "canonical column is umol/L, converted inline in model() (1 mg/dL = 17.1 umol/L, so the",
        "cut-off is 68.4 umol/L). Bilirubin <= 4 mg/dL is the reference. Table 1 does not report",
        "the bilirubin distribution of the cohort.",
        sep = " "
      ),
      source_name = "bilirubin"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened; not retained. Results: 'An influence of the patients' body weight, especially with obese patients, as reported by Hall et al. was not detected.' Table 1 median 94.5 kg (49.9-162)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Methods); not retained. Table 1 median 65 years (22-84)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened (Methods); not retained. Table 1 median 29.87 kg/m^2 (18.15-93.66)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Methods: 'As categorical covariate, sex was tested'; not retained. 12 of 36 patients were female."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened to depict changes in protein binding (Methods); not retained. Table 1 median 2.5 g/dL (1.5-4.2). Total plasma protein (median 5.3 g/dL, 2.5-7.3) was screened for the same reason and also not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened (Methods); not retained. Table 1 median 53 U/L (8-4615). AST (median 110, 18-12550) and butyrylcholinesterase (median 2214, 509-6357 U/L) were screened as further hepatic markers and not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened (Methods); not retained. Table 1 median 110 U/L (18-12550)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Documented for patients not on renal replacement therapy (Methods); not retained. Table 1 median 1.1 mg/dL (0.3-6.2)."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Documented for patients not on renal replacement therapy, in mL/min (not BSA-normalised); not retained. Table 1 median 97 mL/min (23-286)."
    ),
    APACHE_II = list(
      description = "APACHE II score at ICU admission",
      units = "(score)",
      type = "continuous",
      notes = "Recorded on ICU admission (Methods); not retained. Table 1 median 35 (13-46)."
    ),
    SAPS_II = list(
      description = "Simplified Acute Physiology Score II, recorded daily",
      units = "points",
      type = "continuous",
      notes = "Documented each day in the ICU (Methods); not retained. Table 1 median 55 (14-104)."
    ),
    RRT_CRRT_STATUS = list(
      description = "Renal replacement therapy (CVVHDF, SLEDD or intermittent haemodialysis) indicator",
      units = "(binary)",
      type = "binary",
      notes = "The study's main question. Results: the impact of RRT group was tested on the population parameters and 'revealing only an influence of the SOFA score'; Abstract: 'Renal replacement therapy (RRT) did not influence micafungin clearance or the volumes of distribution.' 23 of 36 patients received RRT (13 CVVHDF, 3 SLEDD, 5 intermittent haemodialysis, 2 CVVHDF and SLEDD)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 36L,
    n_studies = 1L,
    age_range = "22-84 years (median 65)",
    weight_range = "49.9-162 kg (median 94.5)",
    sex_female_pct = 33.3,
    race_ethnicity = "Not reported (single-centre German cohort).",
    disease_state = paste(
      "Critically ill adult ICU patients receiving micafungin as empirical (25) or targeted (11)",
      "therapy of suspected or proven Candida infection. 78% met sepsis criteria, 75% received",
      "vasopressors, 86% were mechanically ventilated, 64% received renal replacement therapy and",
      "16% ECMO support. APACHE II median 35 (13-46), SOFA median 12 (2-22), SAPS II median 55",
      "(14-104).",
      sep = " "
    ),
    dose_range = paste(
      "Micafungin 100 mg intravenously once daily; median treatment length 10.5 days (2-79).",
      "The infusion duration is not stated in the paper.",
      sep = " "
    ),
    regions = "Germany (University Hospital Muenster).",
    indication = "Empirical or targeted treatment of invasive candidiasis in the ICU.",
    notes = paste(
      "Sampling (Methods): day 1 at 1, 4, 8, 12, 16, 20 h and a pre-dose trough at 24 h; from day 2,",
      "a peak one hour after infusion and a trough. Up to 10 days of treatment were sampled.",
      "Plasma micafungin by HPLC (LLOQ 0.1 ug/mL). Table 1 gives the weight range as 49.9-162",
      "kg while the Results text gives 49.4-162 kg; the table value is used here.",
      sep = " "
    )
  )

  ini({
    # Table 2 'Final population pharmacokinetic model estimates'. The typical
    # values apply at SOFA = 12 (CL), SOFA <= 10 (Vc) and bilirubin <= 4 mg/dL.
    lcl <- log(1.56)
    label("Clearance at SOFA = 12 and bilirubin <= 4 mg/dL (L/h)") # Table 2: Clearance 1.56 L/h (RSE 4.8%)
    lvc <- log(16.2)
    label("Central volume at SOFA <= 10 (L)") # Table 2: Central volume of distribution 16.2 L (RSE 12.7%)
    lq <- log(14.4)
    label("Intercompartmental clearance (L/h)") # Table 2: Inter-compartmental clearance 14.4 L/h (RSE 7.7%)
    lvp <- log(13.8)
    label("Peripheral volume (L)") # Table 2: Peripheral volume of distribution 13.8 L (RSE 15.7%)

    e_sofa_cl <- -0.024
    label("Linear fractional change in CL per SOFA point above 12 (1/score)") # Results 'Population pharmacokinetic modelling': linear decrease in clearance with a factor of 0.024 per score point
    e_sofa_vc <- -0.308
    label("Fractional change in Vc when SOFA > 10 (unitless)") # Results: 30.8% lower central volume of distribution with a SOFA > 10
    e_tbili_cl <- -0.211
    label("Fractional change in CL when bilirubin > 4 mg/dL (unitless)") # Results: 21.1% lower clearance with bilirubin > 4 mg/dl

    # IIV: Table 2 reports each as a percentage. The paper states no CV
    # definition, so the exact log-normal conversion omega^2 = log(1 + CV^2)
    # is used. Results: 'a full covariance matrix was estimated', but the
    # off-diagonal elements are not reported, so the etas are independent.
    etalcl ~ log(1 + 0.259^2) # Table 2: Inter-individual Variability Clearance 25.9% (RSE 9.8%)
    etalvc ~ log(1 + 0.489^2) # Table 2: Inter-individual Variability Central volume of distribution 48.9% (RSE 13.4%)
    etalvp ~ log(1 + 0.70^2) # Table 2: Inter-individual Variability Peripheral volume of distribution 70% (RSE 13%)

    propSd <- 0.257
    label("Proportional residual error (fraction)") # Table 2: Proportional residual Error 0.257 (RSE 7.4%), read as a 25.7% SD
  })

  model({
    # Bilirubin threshold is in mg/dL in the source; TBILI is in umol/L.
    tbili_mgdl <- TBILI / 17.1
    bili_high <- 0
    if (tbili_mgdl > 4) bili_high <- 1
    sofa_high <- 0
    if (SOFA > 10) sofa_high <- 1

    cl <- exp(lcl + etalcl) * (1 + e_sofa_cl * (SOFA - 12)) * (1 + e_tbili_cl * bili_high)
    vc <- exp(lvc + etalvc) * (1 + e_sofa_vc * sofa_high)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
