Adamiszak_2025_fluconazole <- function() {
  description <- "One-compartment intravenous population PK model for fluconazole in hemato-oncologic pediatric patients receiving once-daily 0.5-1 h infusions for Candida spp. prophylaxis, with body-weight allometric scaling referenced to 70 kg (exponent fixed at 0.75 on CL and at 1.0 on V), log-normal between-subject variability on CL and V, and a proportional residual error. Developed in nlmixr2/FOCEI from 35 plasma concentrations in nine children aged 7 months to 18 years, and used to run probability-of-target-attainment simulations against an fAUC/MIC target."
  reference   <- paste(
    "Adamiszak A, Derwich K, Bartkowska-Sniatkowska A, Pietrzkiewicz K,",
    "Niewiadomska-Wojnalowicz I, Czyrski A, Jusko WJ, Bienert A (2025).",
    "Fluconazole dosing for the prevention of Candida spp. infections in",
    "hemato-oncologic pediatric patients: population pharmacokinetic modeling",
    "and probability of target attainment simulations.",
    "Pharmaceutics 17(4):488.",
    "doi:10.3390/pharmaceutics17040488."
  )
  vignette <- "Adamiszak_2025_fluconazole"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "fluconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight; the only covariate retained in the final model, entering as allometric scaling on both CL and V.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference weight 70 kg (Adamiszak 2025 Equations (1) and (2): CL = CL_T * (WT/70)^0.75,",
        "V = V_T * (WT/70)^1, where CL_T and V_T are typical values for a 70 kg adult).",
        "Observed range in the analysis population 6.0-58.5 kg, median 28.5 kg (Table 1).",
        "The supplementary nlmixr2 code names this data column 'Weight'."
      ),
      source_name        = "Weight"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age in years.",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened by stepwise covariate modeling and by the mlcov package but not retained in the",
        "final model. Reported only as a post hoc regression of body-weight-normalized CL against",
        "age (Adamiszak 2025 Section 3.2 and Figure 3): median CL 0.68 mL/min/kg below the median",
        "age of 9.75 years vs 0.32 mL/min/kg above it. No point estimate of an age coefficient is",
        "published, so the effect cannot be encoded."
      )
    ),
    CRCL = list(
      description = "BSA-normalized renal function as estimated glomerular filtration rate, computed with both the Bedside Schwartz (2009) and the Schwartz 2012 equations.",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = paste(
        "Screened but explicitly not significant: 'The effect of eGFR on fluconazole CL was not",
        "significant' (Adamiszak 2025 Section 3.2). Reported only as a median split at",
        "117.9 mL/min/1.73m2 (CL/BW 0.32 vs 0.56 mL/min/kg, Figure 3). No coefficient is published."
      )
    ),
    SEXF = list(
      description = "Female-sex indicator (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Collected as a baseline characteristic (Adamiszak 2025 Table 1, 6 of 9 female) and screened, but not retained in the final model."
    ),
    CREAT = list(
      description = "Serum creatinine.",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Collected as a baseline characteristic (Adamiszak 2025 Table 1, median 0.31 mg/dL) and",
        "screened. Not retained in the final model, although the Discussion notes that earlier",
        "fluconazole popPK analyses identified serum creatinine as a covariate on CL."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 9L,
    n_studies      = 1L,
    age_range      = "0.50-18.00 years (7 months to 18 years)",
    age_median     = "9.75 years",
    weight_range   = "6.00-58.50 kg",
    weight_median  = "28.50 kg",
    sex_female_pct = 66.7,
    race_ethnicity = NULL,
    disease_state  = "Hemato-oncologic pediatric inpatients receiving intravenous fluconazole for prophylaxis of Candida spp. infections.",
    dose_range     = "3-11 mg/kg intravenous fluconazole once daily as a 0.5 or 1 h infusion (registered SmPC prophylaxis doses).",
    regions        = "Poland (single centre: Karol Jonscher Teaching Hospital, Poznan University of Medical Sciences).",
    renal_function = "Bedside Schwartz eGFR median 151.0 mL/min/1.73m2 (range 115.5-240.3); Schwartz 2012 eGFR median 117.77 mL/min/1.73m2 (range 95.31-170.03) (Table 1).",
    notes          = paste(
      "Adamiszak 2025 Section 3.1 and Table 1. Prospective opportunistic-sampling study",
      "(ClinicalTrials.gov NCT05426499), patients recruited December 2022 to October 2024.",
      "35 plasma concentrations from 9 patients (median 4 samples per patient, range 2-5);",
      "all concentrations were above the 0.5 mg/L lower limit of quantification of the HPLC-UV assay.",
      "All 9 patients contributed to the popPK fit; patient ID 9 (6 kg, 7 months, CL 0.006 L/h/kg)",
      "was excluded only from the post hoc CL-versus-covariate regressions of Figures 2 and 3, where",
      "it is plotted as a red outlier."
    )
  )

  ini({
    # =========================================================================
    # Structural parameters, Adamiszak 2025 Table 2 ("A One-Compartment Model
    # with Allometric Scaling (Final Model)"), referenced to a 70 kg adult per
    # Equations (1) and (2).
    # =========================================================================
    lcl <- log(1.24)   ; label("Clearance CL_T (L/h) for a 70 kg adult")               # Adamiszak 2025 Table 2, Clearance CL = 1.24 L/h (%RSE 23.23); consistent with the reported 0.018 L/h/kg for a 70 kg adult
    lvc <- log(104.07) ; label("Central volume of distribution V_T (L) for a 70 kg adult") # Adamiszak 2025 Table 2, Volume V = 104.07 L (%RSE 21.59); consistent with the reported 1.49 L/kg

    # =========================================================================
    # Allometric exponents, held fixed by the authors (Table 2 rows
    # "Weight (kg), WT on CL" = "fixed 0.75" and "Weight (kg), WT on V" =
    # "fixed 1.00"), matching the exponents written into Equations (1) and (2).
    # =========================================================================
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on CL for WT/70 (unitless)") # Adamiszak 2025 Table 2 "Weight (kg), WT on CL: fixed 0.75"; Equation (1)
    e_wt_vc <- fixed(1.00) ; label("Allometric exponent on V for WT/70 (unitless)")  # Adamiszak 2025 Table 2 "Weight (kg), WT on V: fixed 1.00"; Equation (2)

    # =========================================================================
    # Between-subject variability. Table 2 reports IIV as %CV with the footnote
    # "CV%, coefficient of variation calculated as sqrt(exp(omega^2) - 1) x 100%",
    # so the variance is recovered as omega^2 = log(1 + CV^2):
    #   CL: CV = 88.54% -> omega^2 = log(1 + 0.8854^2) = 0.578821
    #   V : CV = 55.85% -> omega^2 = log(1 + 0.5585^2) = 0.271493
    # The published nlmixr2 code (Supplementary Materials, "nlmixr2 final model
    # code") specifies diagonal (uncorrelated) etas on CL and V.
    # =========================================================================
    etalcl ~ 0.578821 # Adamiszak 2025 Table 2, IIV on CL = 88.54 CV%, back-transformed via the Table 2 footnote
    etalvc ~ 0.271493 # Adamiszak 2025 Table 2, IIV on V = 55.85 CV%, back-transformed via the Table 2 footnote

    # =========================================================================
    # Residual unexplained variability. "The proportional error model resulted
    # in the best data fit" (Section 3.2); Table 2 reports it as 25.19 CV%,
    # i.e. a proportional standard deviation of 0.2519 on the linear scale
    # (Supplementary "nlmixr2 final model code": cp ~ prop(prop.err)).
    # =========================================================================
    propSd <- 0.2519 ; label("Proportional residual error (fraction)") # Adamiszak 2025 Table 2, Proportional residual error = 25.19 CV%
  })

  model({
    # 1. Individual parameters with body-weight allometric scaling to a 70 kg
    #    reference adult (Adamiszak 2025 Equations (1) and (2)).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    # 2. Micro-constant.
    kel <- cl / vc

    # 3. One-compartment disposition with first-order elimination from the
    #    central compartment; fluconazole is given intravenously, so the dose
    #    enters `central` directly (as a zero-order infusion in the source
    #    data, dur = 0.5 or 1 h).
    d/dt(central) <- -kel * central

    # 4. Observation: total plasma fluconazole concentration in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
