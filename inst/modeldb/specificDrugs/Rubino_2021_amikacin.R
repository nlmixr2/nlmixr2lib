Rubino_2021_amikacin <- function() {
  description <- paste(
    "One-compartment population PK model for amikacin in serum after once-daily",
    "nebulised amikacin liposome inhalation suspension (ALIS) 590 mg in adults with",
    "treatment-refractory nontuberculous mycobacterial (NTM) lung disease (Rubino 2021).",
    "Zero-order nebulised input into a lung absorption compartment, first-order",
    "absorption from lung to serum, and linear apparent clearance CLt/F from an",
    "apparent central volume Vc/F; a urine compartment accumulates the renally",
    "excreted amount at the true renal clearance CLr (allometrically scaled by body",
    "weight to a 51 kg reference), so the serum concentration and the cumulative",
    "urine amount were fit simultaneously and the systemic bioavailability is implied",
    "by CLr / (CLt/F).",
    sep = " "
  )
  reference <- paste(
    "Rubino CM, Onufrak NJ, van Ingen J, Griffith DE, Bhavnani SM, Yuen DW, Mange KC,",
    "Winthrop KL. Population Pharmacokinetic Evaluation of Amikacin Liposome Inhalation",
    "Suspension in Patients with Treatment-Refractory Nontuberculous Mycobacterial Lung",
    "Disease. Eur J Drug Metab Pharmacokinet. 2021;46(2):277-287.",
    "doi:10.1007/s13318-020-00669-7.",
    "Correction: Eur J Drug Metab Pharmacokinet. 2021. doi:10.1007/s13318-021-00687-z",
    "(corrects the Figure 2 x-axis labels and a figure cross-reference only; no model",
    "value is affected).",
    "The structural model, zero-order lung input and additive residual-error form are",
    "carried from the cystic fibrosis model it was based on: Okusanya OO, Bhavnani SM,",
    "Hammel JP, Forrest A, Bulik CC, Ambrose PG, Gupta R. Antimicrob Agents Chemother.",
    "2014;58(9):5005-5015. doi:10.1128/AAC.02421-13.",
    sep = " "
  )
  vignette <- "Rubino_2021_amikacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline body weight. Scales the renal clearance CLr only, as",
        "CLr = 1.931 * (WT / 51)^0.75 (Rubino 2021 Table 2 footnote; exponent fixed,",
        "no %SEM). No other structural parameter carries a covariate. Cohort median",
        "52.6 kg, range 33.8-80 kg (Table 1).",
        sep = " "
      ),
      source_name = "WTKG"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "amikacin",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "amikacin",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    ),
    urine = list(
      analyte = "amikacin",
      units = "mg",
      specimen = "urine",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 53L,
    n_studies = 2L,
    age_range = "20-84 years",
    age_median = "63 years",
    weight_range = "33.8-80 kg",
    weight_median = "52.6 kg",
    sex_female_pct = 84.9,
    race_ethnicity = c(Japanese = 52.8, White = 47.2),
    disease_state = paste(
      "Treatment-refractory nontuberculous mycobacterial lung disease (sputum",
      "persistently culture-positive for Mycobacterium avium complex, or M. abscessus",
      "in TR02-112, despite >= 6 months of guideline-based multidrug therapy)",
      sep = " "
    ),
    dose_range = paste(
      "Amikacin liposome inhalation suspension 590 mg once daily by eFlow nebuliser",
      "for up to 168 days (TR02-112) or up to 16 months (CONVERT)",
      sep = " "
    ),
    regions = "USA, Canada, Japan",
    renal_function = "eGFR median 88.3 (range 57.4-140) mL/min/1.73 m^2",
    notes = paste(
      "Pooled phase 2 TR02-112 (n = 14, all White, all female; 111 serum and 23",
      "urine samples) and phase 3 CONVERT (n = 39, 28 Japanese and 11 White; 380",
      "serum samples, no urine) pharmacokinetic substudies. 89 of 491 serum",
      "concentrations were below the 0.15 mg/L LLOQ and were fit with the Beal M3",
      "method. Demographics from Rubino 2021 Table 1; sample counts from Results",
      "section 3.2. Fit with MC-PEM in S-ADAPT 1.57.",
      sep = " "
    )
  )

  ini({
    # Structural parameters: Rubino 2021 Table 2, 'Population mean / Final estimate'.
    # The model is written in the apparent (/F) frame of the serum data: the full
    # nebulised dose enters the lung depot and CLt/F, Vc/F carry the unknown
    # systemic bioavailability, while CLr is a true clearance fit to the urine
    # amounts (Figure 1; Methods section 2.3).
    lcl <- log(34.29); label("Apparent total clearance CLt/F (L/h)") # Table 2: CLt/F = 34.29 L/h (16.36 %SEM)
    lvc <- log(272.6); label("Apparent central volume Vc/F (L)") # Table 2: Vc/F = 272.6 L (11.41 %SEM)
    lka <- log(1.866); label("First-order absorption rate constant from lung to serum ka (1/h)") # Table 2: ka = 1.866 1/h (24.86 %SEM)
    lcl_renal <- log(1.931); label("Renal clearance CLr for a 51 kg patient (L/h)") # Table 2: 'CLr coefficient' = 1.931 L/h (14.29 %SEM); footnote: population mean value for CLr in a patient weighing 51 kg
    e_wt_cl_renal <- fixed(0.75); label("Allometric exponent of body weight on CLr (unitless)") # Table 2: 'CLr WTKG power' = 0.75, no %SEM (fixed)

    # Inter-individual variability: Table 2 'Inter-individual variability (%CV)'.
    # Exponential (log-normal) IIV (Methods section 2.3); diagonal, no
    # correlations reported. Converted with omega^2 = log(CV^2 + 1).
    etalcl ~ 0.41595 # Table 2: CLt/F IIV 71.82 %CV -> log(0.7182^2 + 1)
    etalvc ~ 0.35324 # Table 2: Vc/F IIV 65.09 %CV -> log(0.6509^2 + 1)
    etalka ~ 0.15049 # Table 2: ka IIV 40.30 %CV -> log(0.4030^2 + 1)
    etalcl_renal ~ 0.08750 # Table 2: CLr IIV 30.24 %CV -> log(0.3024^2 + 1)

    # Residual error: Table 2. Separate additive error models for serum and urine
    # (Methods section 2.3; the additive form of the base model is stated in
    # Okusanya 2014 Methods 'Pharmacokinetic analysis').
    addSd <- 0.615; label("Additive residual SD on serum amikacin concentration (mg/L)") # Table 2: 'Residual error (serum)' = 0.615 (3.941 %SEM)
    addSd_Aurine <- 14.0; label("Additive residual SD on amikacin amount excreted in urine (mg)") # Table 2: 'Residual error (urine)' = 14.0 (27.30 %SEM)
  })
  model({
    # Individual parameters
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalka)
    cl_renal <- exp(lcl_renal + etalcl_renal) * (WT / 51)^e_wt_cl_renal

    kel <- cl / vc

    # Zero-order nebulised input into the lung (dose record with a duration),
    # first-order lung -> serum absorption, linear elimination at CLt/F. Urine
    # accumulates the renally excreted amount CLr * Cc; it is a parallel
    # bookkeeping state (dotted arrow in Figure 1), not an extra loss from
    # the apparent central compartment, whose total loss is already CLt/F.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(urine) <- cl_renal * central / vc

    Cc <- central / vc
    Aurine <- urine

    Cc ~ add(addSd)
    Aurine ~ add(addSd_Aurine)
  })
}
