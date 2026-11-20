Shenkoya_2026_enoxaparin <- function() {
  description <- paste(
    "One-compartment first-order absorption population PK model of anti-factor Xa",
    "activity after subcutaneous enoxaparin thromboprophylaxis (every 12 h) in 408",
    "burn patients (386 adults, 22 children). Apparent clearance scales with body",
    "weight (power, 87 kg reference); apparent volume scales with body weight and",
    "with the percentage of total body surface area burned (power, 87 kg and 19%",
    "references). Exponential between-subject variability on CL/F and Vc/F (nearly",
    "perfectly correlated), between-occasion variability on CL/F over four",
    "dose-count occasions, and proportional residual error. No between-subject",
    "variability on ka."
  )
  reference <- paste(
    "Shenkoya B, Yellepeddi V, Prazak AMB, Gopalakrishnan M. A real-world",
    "data-driven approach to optimizing enoxaparin dosing in burn patients.",
    "Br J Clin Pharmacol. 2026;92(4):1126-1136. doi:10.1002/bcp.70347"
  )
  vignette <- "Shenkoya_2026_enoxaparin"
  units <- list(time = "h", dosing = "IU", concentration = "IU/mL")

  compartmentData <- list(
    depot = list(analyte = "enoxaparin", units = "IU", specimen = "administration site", verified = TRUE),
    central = list(analyte = "enoxaparin", units = "IU", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Shenkoya 2026 Table 1: median 86 kg (range 11-208). Power covariate on",
        "CL/F (exponent 0.670) and on Vc/F (exponent 1.07), both normalised to",
        "87 kg as printed in the Table 2 parameter equations (the Methods centre",
        "continuous covariates on their median; 87 kg is the value the final",
        "model uses)."
      ),
      source_name = "Weight"
    ),
    BSA_AFFECTED_PCT = list(
      description = "Total body surface area burned (TBSA), percent",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Shenkoya 2026 Table 1: median 13% (range 0.25-91); <10% in 40%,",
        "10-20% in 32%, >20% in 27% of patients. Power covariate on Vc/F",
        "(exponent 0.473) normalised to 19% as printed in the Table 2 Vc/F",
        "equation. Treated as time-fixed (one value per patient). TBSA on CL/F",
        "was tested and removed at backward elimination (Table S2 model 16)."
      ),
      source_name = "TBSA"
    ),
    OCC = list(
      description = "Occasion index for between-occasion variability on CL/F (1-4)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Shenkoya 2026 Discussion: occasions were grouped into four bins by the",
        "number of enoxaparin doses received, <10, 10-19, 20-39 and >=40 doses.",
        "Encoded as OCC = 1 for doses 1-9, 2 for doses 10-19, 3 for doses 20-39",
        "and 4 for dose 40 onward; an observation takes the occasion of the most",
        "recent dose. The paper does not state the exact bin boundary convention",
        "beyond the four ranges. Records with OCC outside 1-4 receive no",
        "between-occasion variability."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Gender on CL/F was retained at backward elimination (removal +14 OFV)",
        "but dropped from the final model 'for lack of biological rationale,",
        "because weight has already been included' (Results; Table S2 model 20)."
      )
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Age on Vc/F was added (Table S2 model 7) and then eliminated (model 17); age on CL/F was not included (model 6)."
    ),
    CRCL = list(
      description = "Glomerular filtration rate by Cockcroft-Gault (capped at 250 mL/min)",
      units = "mL/min",
      type = "continuous",
      notes = "GFR on CL/F and on Vc/F was tested and not included (Table S2 models 10 and 11)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 408L,
    n_studies = 1L,
    n_observations = 1288L,
    n_doses = 15517L,
    age_range = "1-89 years (median 45); 386 adults and 22 children (<18 years), one child under 6 years",
    weight_range = "11-208 kg (median 86)",
    tbsa_range = "0.25-91% (median 13); <10% 40%, 10-20% 32%, >20% 27%",
    sex_female_pct = 29,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Thermal, electrical or chemical burns admitted to the Burn Trauma ICU;",
      "9% with inhalational injury. Serum creatinine median 0.76 mg/dL (0.21-1.48);",
      "Cockcroft-Gault GFR median 141 mL/min (40-250, capped at 250 in 24",
      "patients), 89% with GFR > 90 mL/min."
    ),
    dose_range = paste(
      "Subcutaneous enoxaparin every 12 h for VTE prophylaxis, 6-120 mg per dose",
      "(15 517 doses). Starting dose by the Lin equation (22.8 + 3.3 x TBSA/10 +",
      "1.89 x WT/10 mg, rounded to 10 mg): 30 mg 20%, 40 mg 51%, 50 mg 21%,",
      "60 mg 6%, other 2%."
    ),
    regions = "United States (University of Utah Burn Center, single centre)",
    notes = paste(
      "Retrospective real-world cohort, July 2014 - January 2023 (Shenkoya 2026",
      "Table 1). Anti-Xa by chromogenic assay, LLOQ 0.1 IU/mL; 84 BLQ samples",
      "(6%) excluded (M1). Samples drawn at median 4.13 h (IQR 3.95-4.42) post",
      "dose. Fitted in Pumas 2.5.1 with FOCE."
    )
  )

  ini({
    lka <- log(1.10); label("Absorption rate constant (1/h)") # Table 2 'Kapop 1/h' = 1.10 (RSE 47.4%); Abstract 1.10 1/h
    lcl <- log(1.56); label("Apparent clearance CL/F at 87 kg (L/h)") # Table 2 'CLpop L/h' = 1.56 (RSE 4.67%)
    lvc <- log(22.1); label("Apparent central volume Vc/F at 87 kg and 19% TBSA (L)") # Table 2 'Vcpop L' = 22.1 (RSE 13.3%)

    e_wt_cl <- 0.670; label("Power exponent of (WT/87) on CL/F (unitless)") # Table 2 'theta WT_CL' = 0.670 (RSE 15.4%)
    e_wt_vc <- 1.07; label("Power exponent of (WT/87) on Vc/F (unitless)") # Table 2 'theta WT_Vc' = 1.07 (RSE 34.1%)
    e_bsa_affected_pct_vc <- 0.473; label("Power exponent of (TBSA/19) on Vc/F (unitless)") # Table 2 'theta TBSA_Vc' = 0.473 (RSE 17.8%)

    # Between-subject variability, exponential. Table 2 'Correlation of random
    # effects' prints rho(etaCL, etaVc) = 1.00, a singular (non-invertible)
    # block; the Results text gives the correlation as .99, which is used here.
    etalcl + etalvc ~ c(log(1 + 0.151^2), 0.99 * sqrt(log(1 + 0.151^2) * log(1 + 0.329^2)), log(1 + 0.329^2)) # Table 2 eta CL CV 15.1%, eta Vc CV 32.9%; Results text correlation 0.99

    # Between-occasion variability on CL/F (Methods Eq. 7, kappa ~ N(0, gamma^2)),
    # one variance shared by the four occasions (Discussion: <10, 10-19, 20-39
    # and >=40 doses). Occasions 2-4 repeat the occasion-1 variance.
    etaiov_cl_1 ~ log(1 + 0.201^2) # Table 2 kappa CL CV 20.1%
    etaiov_cl_2 ~ fixed(log(1 + 0.201^2)) # same variance as occasion 1
    etaiov_cl_3 ~ fixed(log(1 + 0.201^2)) # same variance as occasion 1
    etaiov_cl_4 ~ fixed(log(1 + 0.201^2)) # same variance as occasion 1

    propSd <- 0.287; label("Proportional residual error (fraction)") # Table 2 'Proportional %' = 28.7 (RSE 2.96%)
  })
  model({
    # Occasion indicators for the between-occasion variability on CL/F.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4

    # Table 2: CLi = CLpop * exp(etaCL + kappaCL) * (WT/87)^thetaWT_CL;
    # Vci = Vcpop * exp(etaVc) * (WT/87)^thetaWT_Vc * (TBSA/19)^thetaTBSA_Vc;
    # Kai = Kapop.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 87)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 87)^e_wt_vc * (BSA_AFFECTED_PCT / 19)^e_bsa_affected_pct_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in anti-Xa IU (enoxaparin 1 mg = 100 IU); IU / (L x 1000 mL/L) = IU/mL.
    Cc <- central / (vc * 1000)
    Cc ~ prop(propSd)
  })
}
