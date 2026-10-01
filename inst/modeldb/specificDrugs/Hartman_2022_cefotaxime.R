Hartman_2022_cefotaxime <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order elimination for",
    "intravenous cefotaxime in critically ill children (0-18 years)",
    "admitted to a paediatric intensive care unit (POPSICLE study,",
    "Radboudumc, the Netherlands). Clearance and central volume carry",
    "estimated power functions of body weight normalised to the cohort",
    "median of 10.95 kg (exponents 1.11 and 1.18); peripheral volume and",
    "inter-compartmental clearance are not weight-scaled. Correlated",
    "log-normal interindividual variability on clearance and central",
    "volume; log-normal (exponential) residual error (Hartman 2022)."
  )
  reference <- paste(
    "Hartman SJF, Upadhyay PJ, Mathot RAA, van der Flier M, Schreuder MF,",
    "Bruggemann RJ, Knibbe CAJ, de Wildt SN. Population pharmacokinetics of",
    "intravenous cefotaxime indicates that higher doses are required for",
    "critically ill children. J Antimicrob Chemother. 2022;77(6):1725-1732.",
    "doi:10.1093/jac/dkac095. The open-access supplementary data",
    "(dkac095_supplementary_data.docx, retrieved from EuropePMC PMC9155601)",
    "provides the final NONMEM control stream used to confirm the model",
    "structure, the parameter scales and the residual-error form, and",
    "Figure S6 used for the dose-evaluation replication in the vignette."
  )
  vignette <- "Hartman_2022_cefotaxime"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds. Confirmed against Hartman 2022: the model is
  # a two-compartment model with intravenous administration only (Methods
  # 'Cefotaxime dosing': bolus intravenous infusions, 'as intravenous
  # push'; supplement control stream $INPUT 'CMT ; Compartment = 1 (doses
  # and measurements in intravenous comparment)'), and the observations are
  # total cefotaxime in plasma (Methods 'Blood sampling, handling and
  # analysis': 'Total cefotaxime concentrations were analysed in plasma
  # samples').
  compartmentData <- list(
    central = list(analyte = "cefotaxime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefotaxime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at the start of the paediatric ICU admission.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column WT ('Bodyweight of the patient at the start of ICU",
        "admission in kg', supplement control stream $INPUT), i.e. a",
        "baseline value held constant within a patient. Enters clearance",
        "and central volume as estimated power functions normalised to",
        "the cohort median weight of 10.95 kg (Table 1 'Weight (kg)",
        "10.95 [5.2-28.5]'; Table 2 equations; control stream",
        "'WT_median = 10.95'). Observed range 2.7-80 kg (Table 1);",
        "simulations outside that range are extrapolation."
      ),
      source_name = "WT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 52,
    n_studies = 1,
    n_observations = 479,
    age_range = "0.03-17.69 years (median 1.61, IQR 0.17-8.63)",
    weight_range = "2.7-80 kg (median 10.95, IQR 5.2-28.5)",
    sex_female_pct = 38.5,
    race_ethnicity = "Predominantly Caucasian (>90%; Discussion 'Limitations').",
    disease_state = paste(
      "Critically ill children in a level-3 paediatric intensive care unit",
      "receiving intravenous cefotaxime; main admission reasons respiratory",
      "failure (63.5%), neurological impairment (17.3%) and circulatory",
      "failure (9.6%). Median PRISM-3 6 (range 0-16); 90.4% mechanically",
      "ventilated and 36.5% on vasopressive co-medication during the study.",
      "ECMO and kidney-replacement therapy were exclusion criteria."
    ),
    dose_range = paste(
      "Prophylactic 100 mg/kg/day (maximum 4 g/day; selective",
      "decontamination of the digestive tract) or therapeutic 150",
      "mg/kg/day (maximum 12 g/day), as three or four intravenous push",
      "doses per day; received dose median 100 (range 50-151.1)",
      "mg/kg/day; 11 patients (21.2%) on the therapeutic dose."
    ),
    renal_function = paste(
      "Baseline serum creatinine median 28 (range 8-87) umol/L; baseline",
      "Schwartz-2012 eGFR median 89 (range 46-398) mL/min/1.73 m2."
    ),
    regions = "Netherlands (Radboudumc, Nijmegen; single centre).",
    notes = paste(
      "POPSICLE study (NCT03248349), June 2017 - May 2019. Rich sampling",
      "(median 10 samples per patient) over up to 14 study days. Total",
      "plasma cefotaxime by LC-MS/MS, LLOQ 0.100 mg/L, M6 handling of",
      "BLQ samples. Demographics from Table 1 of Hartman 2022."
    )
  )

  ini({
    # --- Structural parameters (typical values at WT = 10.95 kg) ---
    lcl <- log(2.8)
    label("Clearance at the median weight of 10.95 kg (L/h)")
    # Hartman 2022 Table 2, 'CLpop (L/h)': 2.8 (RSE 8.5%); bootstrap median
    # 2.724, 95% CI 2.195-3.205. Control stream THETA(1) 'TVCL'.

    lvc <- log(2.62)
    label("Central volume of distribution at the median weight of 10.95 kg (L)")
    # Hartman 2022 Table 2, 'V1pop (L)': 2.62 (RSE 10.5%); bootstrap median
    # 2.60, 95% CI 1.641-3.721. Control stream THETA(2) 'TVV1'.

    lvp <- log(1.55)
    label("Peripheral volume of distribution (L)")
    # Hartman 2022 Table 2, 'V2pop (L)': 1.55 (RSE 16.3%); bootstrap median
    # 1.571, 95% CI 1.093-2.439. Control stream THETA(3) 'TVV2'; not
    # weight-scaled ('V2 = TVV2').

    lq <- log(1.15)
    label("Inter-compartmental clearance (L/h)")
    # Hartman 2022 Table 2, 'Qpop (L/h)': 1.15 (RSE 18.3%); bootstrap median
    # 1.148, 95% CI 0.321-1.865. Control stream THETA(4) 'TVQ'; not
    # weight-scaled ('Q = TVQ').

    # --- Estimated allometric weight exponents ---
    e_wt_cl <- 1.11
    label("Power exponent of body weight on clearance (unitless; reference 10.95 kg)")
    # Hartman 2022 Table 2, 'Theta1' in CLi = CLpop * (WT/10.95)^Theta1:
    # 1.11 (RSE 8.5%); bootstrap median 1.106, 95% CI 0.901-1.362. Control
    # stream THETA(5) 'COV_WT_CL'. Estimated (Results: 'a power function
    # with estimated allometric exponents').

    e_wt_vc <- 1.18
    label("Power exponent of body weight on central volume (unitless; reference 10.95 kg)")
    # Hartman 2022 Table 2, 'Theta2' in V1i = V1pop * (WT/10.95)^Theta2:
    # 1.18 (RSE 9.3%); bootstrap median 1.11, 95% CI 0.683-1.647. Control
    # stream THETA(6) 'COV_WT_V1'.

    # --- Interindividual variability (OMEGA BLOCK(2), variances) ---
    # Hartman 2022 Table 2 'Interindividual variability': CL 0.359 (RSE
    # 21.2%), 'Block matrix' 0.305, V1 0.581 (RSE 20.9%). These are the raw
    # NONMEM omega-squared values: the Results text gives the final-model
    # IIV as 65.7% (CL) and 88.8% (V1), and sqrt(exp(0.359) - 1) = 0.657
    # and sqrt(exp(0.581) - 1) = 0.888 reproduce both exactly. The block
    # entry 0.305 is the CL-V1 covariance (correlation 0.668).
    etalcl + etalvc ~ c(0.359, 0.305, 0.581)

    # --- Residual error (exponential) ---
    expSd <- sqrt(0.307)
    label("Log-normal residual error SD (log scale)")
    # Hartman 2022 Table 2, 'Residual variability / Proportional error':
    # 0.307 (RSE 15.5%); bootstrap median 0.298, 95% CI 0.220-0.390. The
    # control stream codes 'Y = F*EXP(EPS(1))' with '$SIGMA 0.311 ; ERR
    # PROP', and the supplementary results state 'The error model was
    # implemented log normally to ensure predicted concentrations were not
    # below 0 mg/L'. The tabulated 0.307 is the raw $SIGMA, a variance, on
    # the same raw scale as the omega-squared rows of the same table (the
    # initial estimate 0.311 in $SIGMA sits beside it), so the log-scale SD
    # is sqrt(0.307) = 0.554.
  })

  model({
    # 1. Covariate terms (Hartman 2022 Table 2 equations; control stream
    #    'COV_WT_CL = (WT/WT_median)**THETA(5)',
    #    'COV_WT_V1 = (WT/WT_median)**THETA(6)', WT_median = 10.95)
    f_wt_cl <- (WT / 10.95)^e_wt_cl
    f_wt_vc <- (WT / 10.95)^e_wt_vc

    # 2. Individual parameters
    cl <- exp(lcl + etalcl) * f_wt_cl
    vc <- exp(lvc + etalvc) * f_wt_vc
    vp <- exp(lvp)
    q <- exp(lq)

    # 3. Micro-constants (control stream K10 = CL/V1, K12 = Q/V1,
    #    K21 = Q/V2)
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system (control stream $DES, commented in favour of ADVAN5
    #    but identical: DADT(1) = -K10*A(1) - K12*A(1) + K21*A(2),
    #    DADT(2) = K12*A(1) - K21*A(2)). Intravenous dosing only.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation (total plasma cefotaxime, S1 = V1) and log-normal
    #    residual error (Y = F*EXP(EPS(1)))
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
