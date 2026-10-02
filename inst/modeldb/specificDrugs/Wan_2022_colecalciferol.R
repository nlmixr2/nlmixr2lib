Wan_2022_colecalciferol <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "first-order elimination describing total serum 25-hydroxyvitamin D",
    "(25(OH)D) after oral colecalciferol (vitamin D3) in children aged",
    "1-18 years with chronic kidney disease stages 2-4 (Wan 2022, C3 trial).",
    "The colecalciferol dose (ug) is absorbed directly into an apparent",
    "25(OH)D central compartment whose initial amount is the basal 25(OH)D",
    "concentration times the apparent volume (no endogenous input term), with",
    "a priori allometric weight scaling (exponents 0.75 on CL/F and 1 on V/F,",
    "reference 24 kg), correlated IIV on CL/F and the basal concentration,",
    "and additive residual error on log-transformed concentrations."
  )
  reference <- paste(
    "Wan M, Green B, Iyengar AA, Kamath N, Reddy HV, Sharma J, Singhal J,",
    "Uthup S, Ekambaram S, Selvam S, Rait G, Shroff R, Patel JP.",
    "Population pharmacokinetics and dose optimisation of colecalciferol in",
    "paediatric patients with chronic kidney disease.",
    "Br J Clin Pharmacol. 2022;88(3):1223-1234. doi:10.1111/bcp.15064.",
    sep = " "
  )
  vignette <- "Wan_2022_colecalciferol"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "A priori allometric scaling on CL/F (exponent 0.75) and V/F",
        "(exponent 1), standardised to 24 kg, the median weight of the study",
        "population (Wan 2022 Section 3.1 and Table 2 footnote; supplement",
        "$PK 'TVV = THETA(1)*(WT/24)', 'TVCL = THETA(2)*((WT/24)**0.75)').",
        "Table 1 median 23.9 kg (IQR 16-38)."
      ),
      source_name = "WT"
    )
  )

  # Wan 2022 Section 2.3 tested serum creatinine (scaled by the expected sex-
  # and age-adjusted normal value) and the type of kidney disease on CL/F.
  # Section 3.1: creatinine did not improve the model; glomerular disease
  # improved the OFV by 4.83 in forward inclusion but failed backward
  # elimination (dOFV > 6.635 required). Neither is in the final model and no
  # coefficient is published, so both are documentation only.
  covariatesDataExcluded <- list(
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested on CL/F as the ratio to the sex- and age-adjusted expected",
        "normal creatinine (supplement $INPUT 'ADCR'); did not improve the",
        "model (Section 3.1). Supplement $INPUT gives CR in mg/dL; Table 1",
        "reports a median of 97.3 umol/L (IQR 66.3-153.9)."
      ),
      source_name = "CR"
    ),
    CREAT_REF = list(
      description = "Sex- and age-adjusted expected normal serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Denominator of the scaled creatinine tested on CL/F (Section 2.3,",
        "supplement $INPUT 'ADCR'); not retained."
      ),
      source_name = "ADCR"
    ),
    DIS_GLOMERULAR = list(
      description = "Glomerular (versus non-glomerular) kidney disease indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-glomerular disease)",
      notes = paste(
        "Supplement $INPUT 'NKD ; 0=Non-glomerulus disease 1=Glomerulus",
        "disease'. Improved the OFV by 4.83 on CL/F in forward inclusion but",
        "was removed at backward elimination (Section 3.1); no coefficient is",
        "published. 20 of 83 children (24 percent) had glomerular disease",
        "(Table 1). No canonical register entry exists because no model in",
        "the library retains it; the name is documentation of the paper's",
        "screen only."
      ),
      source_name = "NKD"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "colecalciferol",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "25-hydroxyvitamin D",
      units = "ug",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 83L,
    n_studies = 1L,
    n_observations = 363L,
    age_range = "1-18 years (inclusion); median 9.4 years (IQR 6.2-14)",
    age_median = "9.4 years",
    weight_range = "IQR 16-38 kg",
    weight_median = "23.9 kg",
    sex_female_pct = 30,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Chronic kidney disease stages 2-4 with serum 25(OH)D < 30 ng/mL;",
      "24 percent glomerular, 76 percent non-glomerular disease; median",
      "eGFR 45.2 mL/min/1.73 m^2 (IQR 29-63.6)"
    ),
    dose_range = paste(
      "Oral colecalciferol granules: 3000 IU daily, 25000 IU weekly or",
      "100000 IU monthly for 3-month intensive courses (up to 3 courses),",
      "then 1000 IU daily maintenance for up to 9 months (1 IU = 0.025 ug)"
    ),
    regions = "India (4 centres, 8-18.5 degrees N)",
    baseline_25ohd = "median 18.6 ng/mL (IQR 13.4-23.4)",
    notes = paste(
      "C3 trial (CTRI/2015/11/010180), randomised 1:1:1 (30 daily, 27 weekly,",
      "26 monthly). 363 serum 25(OH)D samples, roughly every 3 months at",
      "assumed steady state, median 4 per child (range 2-5); total 25(OH)D",
      "by ID-LC-MS/MS. Demographics from Wan 2022 Table 1 and Section 3."
    )
  )

  ini({
    # Structural parameters. Final estimates from the supplementary NONMEM
    # control stream ('NONMEM model file with values set to final
    # estimates'); Table 2 prints the same values rounded.
    lka <- fixed(log(0.323)); label("Absorption rate constant Ka (1/h)") # Supplement $THETA(4) '0.323 FIX'; Table 2 Ka 0.323 /h (fixed), from Ocampo-Pelland
    lcl <- log(0.0328); label("Apparent clearance CL/F of 25(OH)D for a 24 kg child (L/h)") # Supplement $THETA(2) 0.0328; Table 2 CL/F 0.033 L/h (RSE 23%); Discussion 0.0328 L/h at 24 kg
    lvc <- log(322); label("Apparent volume of distribution V/F for a 24 kg child (L)") # Supplement $THETA(1) 322; Table 2 V/F 322 L (RSE 31%)
    lrbase <- log(17.2); label("Basal 25(OH)D concentration C0 at time zero (ng/mL)") # Supplement $THETA(3) BC 17.2; Table 2 C0 17.2 ng/mL (RSE 6%)

    # Allometric exponents, added a priori and not estimated.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT/24 on CL/F (unitless)") # Section 3.1 and Table 2 footnote; supplement '((WT/24)**0.75)'
    e_wt_vc <- fixed(1); label("Allometric exponent of WT/24 on V/F (unitless)") # Section 3.1 and Table 2 footnote; supplement 'THETA(1)*(WT/24)'

    # IIV: supplement $OMEGA BLOCK(2) 0.878, -0.273, 0.117 (ETA(1) on CL,
    # ETA(2) on BC). Table 2 prints sqrt(omega) as %CV: sqrt(0.878) = 93.7%,
    # sqrt(0.117) = 34.2%.
    etalcl + etalrbase ~ c(0.878, -0.273, 0.117) # Supplement $OMEGA BLOCK(2); Table 2 'BSV of apparent CL' 93.7% and 'BSV of C0' 34.2%

    # Residual error: NONMEM fitted ln(DV) with Y = log(IPRED) + EPS(1),
    # $SIGMA 0.145 (variance), i.e. log-normal error with SD sqrt(0.145).
    expSd <- 0.381; label("Additive residual error on log 25(OH)D (SD, log scale)") # Supplement $SIGMA 0.145 -> sqrt = 0.381; Table 2 residual 38.1%
  })

  model({
    # Individual parameters (supplement $PK)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 24)^e_wt_cl
    vc <- exp(lvc) * (WT / 24)^e_wt_vc
    rbase <- exp(lrbase + etalrbase)

    kel <- cl / vc

    # One-compartment model with first-order absorption (ADVAN2 TRANS2)
    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central

    # Supplement $PK 'A_0(2)=BC*V': the central compartment starts at the
    # basal concentration. There is no endogenous input, so the basal amount
    # declines with the 25(OH)D elimination rate constant.
    central(0) <- rbase * vc

    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
