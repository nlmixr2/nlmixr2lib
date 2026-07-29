vanRongen_2018_metformin <- function() {
  description <- "One-compartment population PK model for oral metformin in 22 overweight and obese Caucasian adolescents (van Rongen 2018). First-order absorption into a single central compartment with apparent oral clearance (CL/F) and apparent oral volume of distribution (V/F). Total body weight (TBW) enters linearly on CL/F with reference 75.8 kg (study median): CL/F = 1.17 * (1 + 0.0138 * (TBW - 75.8)) L/min. Proportional residual error; IIV on CL/F, V/F, and ka."
  reference <- paste(
    "van Rongen A, van der Aa MP, Matic M, van Schaik RHN, Deneer VHM,",
    "van der Vorst MM, Knibbe CAJ (2018).",
    "Increased Metformin Clearance in Overweight and Obese Adolescents:",
    "A Pharmacokinetic Substudy of a Randomized Controlled Trial.",
    "Pediatric Drugs 20(4):365-374.",
    "doi:10.1007/s40272-018-0293-1.",
    sep = " "
  )
  vignette <- "vanRongen_2018_metformin"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation covariate on CL/F with reference 75.8 kg (study median). Source column 'TBW' renamed to canonical 'WT' on input. Time-fixed at baseline.",
      source_name        = "TBW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 22L,
    n_studies      = 1L,
    age_range      = "11.1-17.5 years",
    age_median     = "14.5 years (mean)",
    weight_range   = "54.7-104.9 kg",
    weight_median  = "75.8 kg (study median, used as reference in covariate model); 79.3 kg (mean)",
    sex_female_pct = 72.7,
    race_ethnicity = c(White = 100),
    disease_state  = "Overweight and obese insulin-resistant adolescents (BMI-SDS > 2.3, HOMA-IR >= 3.4); no T2DM. All of Caucasian (Dutch) descent. 5 overweight, 17 obese.",
    dose_range     = "Oral metformin 500 mg (n=3) or 1000 mg (n=19) twice daily for 37 weeks; PK sub-study sampled an OGTT dose in the fasted state.",
    regions        = "Netherlands (St Antonius Hospital Nieuwegein and Jeroen Bosch Hospital 's-Hertogenbosch)",
    notes          = "Demographics from Table 1 of van Rongen 2018. NCT01487993 / EudraCT 2010-023980-17. The 75.8 kg reference weight used by the final covariate model on CL/F is the study median TBW (paper Section 3.3 / Table 2)."
  )

  ini({
    # Final covariate-model parameter estimates from Table 2 of van Rongen
    # 2018 (page 369), column "Final covariate model (RSE%)". Bootstrap
    # 95% CIs from 995/1000 successful resamples are reported in the
    # right-most column of Table 2 and agree with the point estimates.
    # The paper also presents an "Excess weight covariate model" (Table 2
    # column 2) that performs similarly (OFV -140.8 vs -141.0); the
    # vignette notes the alternative under "Assumptions and deviations".

    # Structural disposition - first-order absorption into single central
    # compartment. ka is reported in 1/min and V/F in L; CL/F in L/min.
    lka  <- log(0.0248) ; label("Absorption rate constant ka (1/min)")            # Table 2 final Ka = 0.0248 1/min (RSE 25%)
    lvc  <- log(485)    ; label("Apparent oral volume of distribution V/F (L)")   # Table 2 final V/F = 485 L (RSE 10%)
    lcl  <- log(1.17)   ; label("Apparent oral clearance CL/F at TBW = 75.8 kg (L/min)") # Table 2 final CL/F_{75.8 kg} = 1.17 L/min (RSE 6%)

    # Covariate effect - linear TBW slope on CL/F. The final covariate
    # model is parameterized as:
    #   CL/F = CL/F_{75.8 kg} * (1 + V_slope * (TBW - 75.8))
    # where V_slope = 0.0138 per kg (Table 2 final V = 0.0138, RSE 44%).
    e_wt_cl <- 0.0138 ; label("Linear TBW slope on CL/F (per kg above 75.8 kg)") # Table 2 final V = 0.0138 (RSE 44%)

    # Inter-individual variability. Table 2 reports IIV as CV% on the
    # log-normal scale; the internal variance is omega^2 = log(1 + CV^2).
    # Shrinkage in square brackets is reproduced from Table 2 for the
    # reader's reference. No IIV is reported on the residual error.
    etalcl ~ 0.06937  # log(1 + 0.268^2) -- Table 2 final IIV CL/F = 26.8% (RSE 31%) [shrinkage 6%]
    etalvc ~ 0.12305  # log(1 + 0.362^2) -- Table 2 final IIV V/F  = 36.2% (RSE 16%) [shrinkage 14%]
    etalka ~ 0.40169  # log(1 + 0.703^2) -- Table 2 final IIV ka   = 70.3% (RSE 22%) [shrinkage 32%]

    # Residual error - proportional only. Table 2 final proportional
    # error = 21.8% CV (RSE 25%) [shrinkage 18%].
    propSd <- 0.218 ; label("Proportional residual SD (fraction)")               # Table 2 final proportional error = 21.8% (RSE 25%)
  })

  model({
    # Individual structural parameters. TBW (WT) enters linearly on CL/F
    # with reference 75.8 kg per the final covariate model (Table 2):
    #   CL/F_i = CL/F_{75.8} * (1 + V * (TBW_i - 75.8)) * exp(etalcl_i)
    # No covariate is retained on V/F or ka in the final model
    # (Section 3.3: "No other covariates were identified as a
    # significant covariate for any of the pharmacokinetic parameters").
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) * (1 + e_wt_cl * (WT - 75.8))

    kel <- cl / vc

    # One-compartment first-order oral absorption. Dose in mg, V/F in L,
    # so Cc = central / V/F is in mg/L, matching the assay calibration
    # range (0.2-5.0 mg/L; Methods Section 2.3).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
