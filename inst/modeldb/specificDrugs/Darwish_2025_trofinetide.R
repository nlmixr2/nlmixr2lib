Darwish_2025_trofinetide <- function() {
  description <- "Population PK model for oral trofinetide in Rett syndrome (Darwish 2025): two-compartment with first-order absorption and linear elimination, pooled across healthy volunteers and patients with Rett syndrome, fragile X syndrome, or traumatic brain injury."
  reference <- "Darwish M, Passarell J, Maxwell K, Youakim JM, Bradley H, Bishop KM. Population Pharmacokinetic Modeling to Support Trofinetide Dosing for the Treatment of Rett Syndrome. Advances in Therapy. 2025;42(2):1026-1043. doi:10.1007/s12325-024-03056-9"
  vignette <- "Darwish_2025_trofinetide"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Trofinetide concentrations were quantified in
  # lithium-heparinized WHOLE BLOOD (Methods, "Analysis Dataset"), not plasma.
  compartmentData <- list(
    depot       = list(analyte = "trofinetide", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "trofinetide", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "trofinetide", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on clearance; reference 58 kg, the median body weight of the analysis population (Darwish 2025 Table 2 footnote, 'WTKG/58'). Analysis-population range 13-140 kg.",
      source_name        = "WTKG"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on central volume; reference 22.4 years, the median age of the analysis population (Darwish 2025 Table 2 footnote, 'AGE/22.4'). Analysis-population range 5-64 years.",
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = "Glomerular filtration rate, BSA-normalized (creatinine-based estimate)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on clearance; reference 124 mL/min/1.73 m^2, the median GFR of the analysis population (Darwish 2025 Table 2 footnote, 'GFR/124'). The source column is named GFR and maps to the canonical general-scope CRCL covariate, which covers BSA-normalized renal function from either a creatinine-based estimate or a tracer-measured GFR. Darwish 2025 does not state the estimating equation used.",
      source_name        = "GFR"
    ),
    DIS_RETT = list(
      description        = "Rett syndrome disease-state indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteers and any non-Rett cohort)",
      notes              = "1 = patient with Rett syndrome, 0 = all other subjects (Darwish 2025 Results, definition of Rett_I). Decreases CL by 16.9% and increases Vp by 61.6%. Also selects the disease-cohort residual-error magnitude, which Darwish 2025 pooled across Rett syndrome, fragile X syndrome, and traumatic brain injury. All Rett syndrome participants were female.",
      source_name        = "Rett"
    ),
    DIS_TBI = list(
      description        = "Traumatic brain injury disease-state indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteers and any non-TBI cohort)",
      notes              = "1 = patient with traumatic brain injury, 0 = all other subjects (Darwish 2025 Results, definition of TBI_I). Increases CL by 23.5% and decreases Vp by 75.2%. Also selects the disease-cohort residual-error magnitude. NOTE: distinct from the register's DIS_BURN_RECENT canonical, whose documented source alias 'TBI' denotes recent burn injury rather than traumatic brain injury.",
      source_name        = "TBI"
    ),
    DIS_FXS = list(
      description        = "Fragile X syndrome disease-state indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteers and any non-FXS cohort)",
      notes              = "1 = patient with fragile X syndrome, 0 = all other subjects (Darwish 2025 Results, definition of FXS_I). Increases Vc by 115%. Also selects the disease-cohort residual-error magnitude. All fragile X syndrome participants were male.",
      source_name        = "FXS"
    ),
    FED = list(
      description        = "Fed-vs-fasted state at dosing",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "1 = dose administered in the fed state, 0 = fasted (Darwish 2025 Results, definition of Fed_I). Decreases ka by 9.49% and F1 by 13.3%. Estimated from the dedicated phase 1 food-effect study ACP-2566-006; Darwish 2025 does not state the meal composition, so the general FED canonical applies rather than FED_HIGHFAT.",
      source_name        = "FED"
    ),
    DOSE_18G = list(
      description        = "18 g trofinetide dose-level indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any dose level other than 18 g, including the 6-12 g therapeutic weight-banded doses)",
      notes              = "1 = subject received the 18 g supratherapeutic dose, 0 = all other subjects (Darwish 2025 Results, definition of DoseGrp1_I). Decreases F1 by 13.2%. Arises from the thorough-QTc study ACP-2566-008, which included supratherapeutic 18 g and 24 g doses; the effect captures the less-than-proportional rise in exposure above the therapeutic range. Mutually exclusive with DOSE_24G.",
      source_name        = "DoseGrp1"
    ),
    DOSE_24G = list(
      description        = "24 g trofinetide dose-level indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any dose level other than 24 g, including the 6-12 g therapeutic weight-banded doses)",
      notes              = "1 = subject received the 24 g supratherapeutic dose, 0 = all other subjects (Darwish 2025 Results, definition of DoseGrp2_I). Decreases F1 by 28.4%. Arises from the thorough-QTc study ACP-2566-008. Mutually exclusive with DOSE_18G.",
      source_name        = "DoseGrp2"
    ),
    AE_DIARRHEA = list(
      description        = "Concurrent treatment-emergent diarrhea indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no diarrhea at the time of the dose record)",
      notes              = "1 = subject experiencing diarrhea, 0 = otherwise. Darwish 2025 Results explicitly describes Diar_I as a TIME-VARYING indicator, so the value is carried per dose record rather than per subject. Decreases F1 by 14.8%. Diarrhea is the most common trofinetide adverse event; 52.4% of the Rett syndrome participants in the analysis dataset experienced it at some point during the studies.",
      source_name        = "Diar"
    )
  )

  # Covariates that Darwish 2025 screened in the stepwise covariate analysis
  # (Methods, "Covariates of interest") but did NOT retain in the final model.
  # Documented here so the provenance of the paper's screen survives without
  # declaring covariates that `model()` never references.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as part of the composite 'sex/disease state' covariate (six categories: male and female healthy volunteers, RTT [female only], FXS [male only], male TBI, female TBI). Backward elimination removed every female-healthy-subject shift -- from Vp, then Vc, then CL (Darwish 2025 Table 1, eliminations 3, 4, and 7) -- so sex itself is not in the final model. The retained members of that composite are the disease-state indicators DIS_RETT, DIS_TBI, and DIS_FXS."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Listed among the covariates of interest (Darwish 2025 Methods) but not retained; body weight was the size descriptor carried into the final model."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function marker (Darwish 2025 Methods) but not retained in the final model."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function marker (Darwish 2025 Methods) but not retained in the final model."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function marker (Darwish 2025 Methods) but not retained in the final model. Darwish 2025 does not report the units in which it was screened."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 442,
    n_observations = 5595,
    n_studies      = 13,
    age_range      = "5-64 years (median 21)",
    weight_range   = "13-140 kg (median 62.2)",
    gfr_reference  = "124 mL/min/1.73 m^2 (analysis-population median)",
    disease_state  = "Pooled analysis of 156 healthy adult volunteers, 185 patients with Rett syndrome (female only), 44 patients with fragile X syndrome (male only), and 57 patients with traumatic brain injury.",
    dose_range     = "Oral, gastric-tube, and intravenous (bolus and infusion) trofinetide; oral doses spanned the 6-12 g therapeutic weight-banded range plus supratherapeutic 18 g and 24 g dose levels.",
    regions        = "Not reported",
    notes          = "Darwish 2025 Results, 'Population Pharmacokinetic Model'. Pooled across eight phase 1 studies (five in healthy volunteers), four phase 2 studies (Rett syndrome, fragile X syndrome, traumatic brain injury), and the phase 3 LAVENDER study in girls and women aged 5-20 years with Rett syndrome. Assay: LC-MS/MS in lithium-heparinized whole blood, LLOQ 0.10 ug/mL. The target steady-state exposure range used to confirm the weight-banded regimen was AUC0-12 800-1200 ug*h/mL."
  )

  ini({
    # Structural parameters and covariate effects -- Darwish 2025 Table 2 (final
    # parameter estimates) and the Results typical-value equations for F1, ka,
    # CL, Vc, and Vp. Reference values are the analysis-population medians:
    # WT 58 kg, AGE 22.4 years, GFR 124 mL/min/1.73 m^2.
    lcl           <- log(11.8);   label("Clearance at the reference covariate values (L/h)")                        # Table 2: CL = 11.8 L/h (RSE 2.02%)
    e_wt_cl       <- 0.443;       label("Power exponent on (WT/58) for clearance (unitless)")                       # Table 2: Exponent of (WTKG/58) for CL = 0.443 (RSE 8.97%)
    e_crcl_cl     <- 0.273;       label("Power exponent on (CRCL/124) for clearance (unitless)")                    # Table 2: Exponent of (GFR/124) for CL = 0.273 (RSE 20.2%)
    e_rett_cl <- -0.169;      label("Proportional shift in clearance for Rett syndrome (fraction)")             # Table 2: Proportional shift in CL for Rett = 1 is -0.169 (RSE 25.6%)
    e_tbi_cl  <- 0.235;       label("Proportional shift in clearance for traumatic brain injury (fraction)")    # Table 2: Proportional shift in CL for TBI = 1 is 0.235 (RSE 17.6%)

    lvc           <- log(24.9);   label("Central volume of distribution at the reference age (L)")                  # Table 2: Vc = 24.9 L (RSE 3.89%)
    e_age_vc      <- 0.549;       label("Power exponent on (AGE/22.4) for central volume (unitless)")               # Table 2: Exponent of (AGE/22.4) for Vc = 0.549 (RSE 8.84%)
    e_fxs_vc  <- 1.15;        label("Proportional shift in central volume for fragile X syndrome (fraction)")   # Table 2: Proportional shift in Vc for FXS = 1 is 1.15 (RSE 19.9%)

    lq            <- log(1.44);   label("Intercompartmental clearance (L/h)")                                       # Table 2: Q = 1.44 L/h (RSE 6.17%)

    lvp           <- log(35.3);   label("Peripheral volume of distribution (L)")                                    # Table 2: Vp = 35.3 L (RSE 5.46%)
    e_rett_vp <- 0.616;       label("Proportional shift in peripheral volume for Rett syndrome (fraction)")     # Table 2: Proportional shift in Vp for Rett = 1 is 0.616 (RSE 35.8%)
    e_tbi_vp  <- -0.752;      label("Proportional shift in peripheral volume for traumatic brain injury (fraction)") # Table 2: Proportional shift in Vp for TBI = 1 is -0.752 (RSE 4.22%)

    lka           <- log(0.391);  label("First-order absorption rate constant in the fasted state (1/h)")           # Table 2: ka = 0.391 1/h (RSE 4.09%)
    e_fed_ka      <- -0.0949;     label("Proportional shift in ka for the fed state (fraction)")                    # Table 2: Shift in ka for FED = -0.0949 (RSE 24.1%)

    lfdepot       <- log(0.828);  label("Oral bioavailability in the fasted, therapeutic-dose, diarrhea-free reference state (fraction)") # Table 2: F1 = 0.828 (RSE 3.65%)
    e_fed_f       <- -0.133;      label("Proportional shift in bioavailability for the fed state (fraction)")       # Table 2: Shift in F1 for FED = -0.133 (RSE 7.28%)
    e_dose18g_f  <- -0.132;      label("Proportional shift in bioavailability for the 18 g dose level (fraction)") # Table 2: Shift in F1 for 18-g dose = -0.132 (RSE 19.9%)
    e_dose24g_f  <- -0.284;      label("Proportional shift in bioavailability for the 24 g dose level (fraction)") # Table 2: Shift in F1 for 24-g dose = -0.284 (RSE 5.66%)
    e_diarrhea_f <- -0.148;    label("Proportional shift in bioavailability during diarrhea (fraction)")         # Table 2: Shift in F1 for Diarrhea = -0.148 (RSE 21.0%)

    # IIV -- Darwish 2025 Table 2 reports the magnitude of interindividual
    # variability as %CV. Methods Eq. (1) defines %CV = sqrt(exp(omega^2) - 1) *
    # 100, so omega^2 = log(1 + CV^2). No IIV was estimated on ka ("NE" in
    # Table 2), so ka carries no eta here.
    etalcl     ~ log(1 + 0.136^2)   # Table 2: CL 13.6 %CV (RSE 17.4%; eta-shrinkage 37.1%)
    etalvc     ~ log(1 + 0.308^2)   # Table 2: Vc 30.8 %CV (RSE 12.8%; eta-shrinkage 31.4%)
    etalq      ~ log(1 + 0.644^2)   # Table 2: Q 64.4 %CV (RSE 18.5%; eta-shrinkage 9.5%)
    etalvp     ~ log(1 + 0.282^2)   # Table 2: Vp 28.2 %CV (RSE 14.6%; eta-shrinkage 39.8%)
    etalfdepot ~ log(1 + 0.208^2)   # Table 2: F1 20.8 %CV (RSE 15.7%; eta-shrinkage 55.0%)

    # Residual error -- Darwish 2025 Methods Eq. (2) specifies a log/exponential
    # residual error model, equivalent to nlmixr2's `lnorm()`. Two separate
    # magnitudes were estimated, one for healthy subjects and one for subjects
    # with Rett syndrome, FXS, or TBI. The Table 2 "Population mean" column
    # holds the NONMEM $SIGMA VARIANCE; the "Final estimate" %CV column holds
    # the corresponding log-scale SD (sqrt(0.0789) = 0.281 -> 28.1 %CV;
    # sqrt(0.137) = 0.370 -> 37.0 %CV). `lnorm()` takes the SD, so both are
    # entered as sqrt(sigma^2).
    expSdHealthy <- sqrt(0.0789); label("Log-scale residual SD, healthy subjects (log units)")                        # Table 2: Residual variability healthy subjects sigma^2 = 0.0789 (RSE 1.54%), reported as 28.1 %CV
    expSdDisease <- sqrt(0.137);  label("Log-scale residual SD, subjects with Rett syndrome, TBI, or FXS (log units)") # Table 2: Residual variability subjects with Rett syndrome, TBI, or FXS sigma^2 = 0.137 (RSE 3.50%), reported as 37.0 %CV
  })

  model({
    # Individual PK parameters. Every covariate term below is transcribed from
    # the Darwish 2025 Results typical-value equations; each categorical effect
    # enters as a proportional shift of the form (1 + theta * indicator).
    cl <- exp(lcl + etalcl) * (WT / 58)^e_wt_cl * (CRCL / 124)^e_crcl_cl *
      (1 + e_tbi_cl * DIS_TBI) * (1 + e_rett_cl * DIS_RETT)
    vc <- exp(lvc + etalvc) * (AGE / 22.4)^e_age_vc * (1 + e_fxs_vc * DIS_FXS)
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) * (1 + e_rett_vp * DIS_RETT) * (1 + e_tbi_vp * DIS_TBI)
    ka <- exp(lka) * (1 + e_fed_ka * FED)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1

    # Bioavailability carries the fed-state, supratherapeutic-dose, and
    # diarrhea effects. Applied to the depot so intravenous dosing (directly
    # into `central`) is unaffected, matching the source dataset in which F1
    # was identifiable only because IV data were pooled with oral data.
    fdepot <- exp(lfdepot + etalfdepot) * (1 + e_fed_f * FED) *
      (1 + e_dose18g_f * DOSE_18G) * (1 + e_dose24g_f * DOSE_24G) *
      (1 + e_diarrhea_f * AE_DIARRHEA)
    f(depot) <- fdepot

    # `central` is in mg and `vc` in L, so central/vc is mg/L = ug/mL, the
    # units in which Darwish 2025 reports trofinetide whole-blood concentrations.
    Cc <- central / vc

    # A single residual-error model switched by the disease-state indicators.
    # The three indicators are mutually exclusive, so their sum is 1 for any
    # patient cohort and 0 for healthy volunteers.
    expSd <- expSdDisease * (DIS_RETT + DIS_TBI + DIS_FXS) +
      expSdHealthy * (1 - DIS_RETT - DIS_TBI - DIS_FXS)
    Cc ~ lnorm(expSd)
  })
}
