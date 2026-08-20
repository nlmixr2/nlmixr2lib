Winning_2025_certepetide <- function() {
  description <- paste(
    "Two-compartment population PK model with linear elimination for certepetide",
    "(LSTA1, CEND-1), a cyclic tumor-penetrating RGD peptide, in adults with",
    "metastatic pancreatic ductal adenocarcinoma receiving certepetide with",
    "nab-paclitaxel and gemcitabine. Body weight scales central and peripheral",
    "volume as estimated power functions; baseline Cockcroft-Gault creatinine",
    "clearance splits clearance into an additive non-renal plus renal arm."
  )
  reference <- paste(
    "Winning A, Sietsema WK, Buck KK, Linsmeier A, Wiczling P (2025).",
    "Population Pharmacokinetic Modeling of Certepetide in Human Subjects With",
    "Metastatic Pancreatic Ductal Adenocarcinoma.",
    "Clin Pharmacol Drug Dev 14(3):240-251. doi:10.1002/cpdd.1502."
  )
  vignette <- "Winning_2025_certepetide"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline weight. Enters Vc and Vp as an estimated power",
        "function (WT / 70)^theta, per Winning 2025 Equation (6). The reference",
        "value 70 kg is the standard reference value named in the Equation (6)",
        "definition list ('a standard reference value (ie, 60 years of age,",
        "70 kg body weight)') and is the reference subject weight used for the",
        "Figure 4 forest plots ('a reference 70-kg subject'). The cohort median",
        "was 73.7 kg (range 54.0-121 kg, Winning 2025 Results, Data Disposition).",
        "Certepetide itself is dosed in mg/kg, so WT enters the simulation twice:",
        "once through the dose amount and once through the volume terms."
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = paste(
        "Baseline creatinine clearance computed by the Cockcroft-Gault equation,",
        "in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2)."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source term 'baseline CrCL' / 'BCRCL'. Stored under the canonical CRCL",
        "column per inst/references/covariate-columns.md, which accepts raw",
        "Cockcroft-Gault mL/min when the source paper applies no BSA",
        "normalization, with the assay form documented per model (same",
        "convention as Delattre_2010_amikacin.R). NOT stored as CRCL_BASE:",
        "that canonical is reserved for the Wahlby 2004 BCOV/DCOV decomposition",
        "in which a time-fixed baseline column is paired with a time-varying",
        "CRCL column. Winning 2025 carries no time-varying renal-function",
        "column -- CrCL is simply time-fixed by data design, which the",
        "CRCL_BASE register entry explicitly routes to CRCL. Enters CL as the",
        "additive renal arm CL = CL_NR + CL_R * (CRCL / 90) per Equation (7);",
        "90 mL/min is stated there as the reference value. Cohort median",
        "96.8 mL/min (range 48.2-172 mL/min)."
      ),
      source_name        = "BCRCL"
    )
  )

  # Covariates screened during model development but not retained in the final
  # model. Both were tested and dropped, so they carry no coefficient and are
  # documentation only -- see Winning 2025 Results, Pharmacokinetic Model.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "An age effect on CL was retained in an intermediate model",
        "(dOFV = -17.07) but was replaced in the final model by baseline CrCL,",
        "because CrCL correlates strongly with both age and body weight and the",
        "CrCL parameterization gave a similar OFV with a more parsimonious",
        "structure. No age coefficient is reported for the final model",
        "(Winning 2025 Table 1)."
      )
    ),
    SEXF = list(
      description = "Sex, 1 = female",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Sex was tested on CL and Vc and was not significant",
        "(dOFV < 10.83); removed in backward elimination",
        "(Winning 2025 Results, Pharmacokinetic Model)."
      )
    )
  )

  compartmentData <- list(
    central     = list(analyte = "certepetide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "certepetide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 31,
    n_studies      = 1,
    n_observations = 1142,
    age_range      = "42.6-79.3 years",
    age_median     = "62.1 years",
    weight_range   = "54.0-121 kg",
    weight_median  = "73.7 kg",
    sex_female_pct = 35.5,
    race_ethnicity = c(White = 87.1, Black = 6.5, `Mixed or other` = 6.5),
    disease_state  = paste(
      "Unresectable metastatic exocrine pancreatic ductal adenocarcinoma,",
      "first-line treatment"
    ),
    co_medication  = paste(
      "nab-paclitaxel 125 mg/m^2 IV over 30 min and gemcitabine 1000 mg/m^2",
      "IV over 30 min on days 1, 8 and 15 of each 28-day cycle; certepetide",
      "given immediately after nab-paclitaxel"
    ),
    dose_range     = paste(
      "0.2, 0.8, 1.6 and 3.2 mg/kg certepetide as a slow intravenous push over",
      "1 minute; monotherapy run-in followed by days 1, 8 and 15 of 28-day",
      "chemotherapy cycles"
    ),
    renal_function = paste(
      "61.3% normal (CrCL >= 90 mL/min), 25.8% mild impairment",
      "(CrCL 60-89 mL/min), 12.9% moderate impairment (CrCL 30-59 mL/min);",
      "median baseline CrCL 96.8 mL/min (48.2-172 mL/min)"
    ),
    hepatic_function = paste(
      "67.7% normal (AST and bilirubin <= ULN), 32.3% mild impairment",
      "(AST > ULN or ULN < bilirubin <= 1.5 * ULN)"
    ),
    regions        = "Australia (3 clinical sites)",
    notes          = paste(
      "First-in-human phase 1 rising-dose study NCT03517176 (Dean et al.).",
      "Baseline demographics: Winning 2025 Results, Data Disposition, and",
      "Tables S3-S4. Of 1142 PK observations, 50 (4.4%) were below the 50.0",
      "ng/mL limit of quantification and were excluded by the M1 method; 8",
      "observations with high predose measurements were also excluded. Assay",
      "was LC-MS/MS calibrated over 50.0-2500 ng/mL."
    )
  )

  ini({
    # -------------------------------------------------------------------
    # Structural parameters. Winning 2025 Table 1 ("Pharmacokinetic
    # Parameters for the Certepetide Final Model"). Volumes are the typical
    # values at the 70 kg reference weight; the clearance arms are the
    # typical values of Equation (7) at the 90 mL/min reference CrCL.
    # -------------------------------------------------------------------
    lvc        <- log(5.87); label("Central volume of distribution Vc at WT = 70 kg (L)")             # Table 1 theta_Vc = 5.87 L (95% CI 5.23, 6.51)
    lvp        <- log(10.3); label("Peripheral volume of distribution Vp at WT = 70 kg (L)")          # Table 1 theta_Vp = 10.3 L (95% CI 9.51, 11.1)
    lq         <- log(24.9); label("Intercompartmental clearance Q (L/h)")                            # Table 1 theta_Q = 24.9 L/h (95% CI 21.7, 28.2)
    lcl_nonren <- log(2.56); label("Non-renal clearance CL_NR, independent of CrCL (L/h)")            # Table 1 theta_CLNR = 2.56 L/h (95% CI 1.41, 3.71)
    lcl_renal  <- log(4.32); label("Renal clearance CL_R at CrCL = 90 mL/min (L/h)")                  # Table 1 theta_CLR = 4.32 L/h (95% CI 3.22, 5.42)

    # -------------------------------------------------------------------
    # Covariate effects. Estimated power exponents on (WT / 70), applied as
    # (WT / 70)^e per Equation (6). No weight effect was retained on CL or Q:
    # the final model carries renal function on CL instead (Results,
    # Pharmacokinetic Model).
    # -------------------------------------------------------------------
    e_wt_vc <- 0.933; label("Power exponent on (WT / 70) for Vc (unitless)")                          # Table 1 theta_WT-Vc = 0.933 (95% CI 0.602, 1.26)
    e_wt_vp <- 0.879; label("Power exponent on (WT / 70) for Vp (unitless)")                          # Table 1 theta_WT-Vp = 0.879 (95% CI 0.558, 1.20)

    # -------------------------------------------------------------------
    # Interindividual variability. Log-normal per Equation (1),
    # theta_i = theta_typical * exp(eta_i); a DIAGONAL Omega was estimated,
    # so there are no off-diagonal terms. Table 1 reports omega^2 directly
    # (the bracketed CV% values follow Equation (2),
    # CV% = 100 * sqrt(exp(omega^2) - 1), which reproduces 21.4 / 17.5 /
    # 14.5 exactly from the variances below).
    #
    # NOTE: the paper fixed the variance for Q to 0 ("The variance for Q and
    # the variance of an additive part of the residual error model was fixed
    # to 0, as they tended to 0 during the model-building process"), so Q
    # carries no eta here. Encoding it as `etalq ~ fixed(0)` would make Omega
    # singular and break rxode2's Cholesky sampler; omitting the eta is the
    # faithful and solvable encoding of a variance fixed at zero.
    # -------------------------------------------------------------------
    etalvc ~ 0.0447  # Table 1 IIV-Vc omega^2_Vc = 0.0447 [CV% = 21.4] (95% CI 0.00996, 0.0795); shrinkage 19.8%
    etalcl ~ 0.0301  # Table 1 IIV-CL omega^2_CL = 0.0301 [CV% = 17.5] (95% CI 0.0136, 0.0467); shrinkage 4.99%
    etalvp ~ 0.0209  # Table 1 IIV-Vp omega^2_Vp = 0.0209 [CV% = 14.5] (95% CI 0.00215, 0.0397); shrinkage 25.7%

    # -------------------------------------------------------------------
    # Residual error. Equation (4) is a combined proportional plus additive
    # model, but the additive variance sigma^2_2 was fixed to 0 during model
    # building (Results, Pharmacokinetic Model), leaving a purely
    # proportional error -- as stated in the Abstract, Discussion and
    # Conclusions. Table 1 reports the variance sigma^2_1 = 0.0393 with
    # [CV% = 19.8]; sqrt(0.0393) = 0.198242, which reproduces the reported
    # 19.8% exactly and confirms the tabulated value is a variance on the
    # linear (proportional) scale.
    # -------------------------------------------------------------------
    propSd <- 0.198242; label("Proportional residual error (fraction)")                               # Table 1 sigma^2_1 = 0.0393 [CV% = 19.8] (95% CI 0.0333, 0.0452); sqrt(0.0393) = 0.198242
  })

  model({
    # -------------------------------------------------------------------
    # 1. Individual PK parameters.
    #
    # Volumes: power model on body weight, Equation (6),
    #   theta_P,i = theta_P,typical * (WT / 70)^theta_P,WT
    #
    # Clearance: additive non-renal plus renal decomposition, Equation (7),
    #   theta_CL,i = theta_CL,NR + theta_CL,R * (BCRCL / 90)
    # with the log-normal eta of Equation (1) applied to the resulting total
    # clearance (Table 1 reports a single IIV-CL term, not one per arm).
    # The log transforms keep both clearance arms positive.
    #
    # Q carries no covariate effect and no IIV (Table 1).
    # -------------------------------------------------------------------
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp
    q  <- exp(lq)
    cl <- (exp(lcl_nonren) + exp(lcl_renal) * (CRCL / 90)) * exp(etalcl)

    # -------------------------------------------------------------------
    # 2. Micro-constants.
    # -------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # -------------------------------------------------------------------
    # 3. Two-compartment ODE system with linear elimination from central.
    # Certepetide is given intravenously (a slow push over 1 minute), so
    # doses target `central` directly; there is no depot and no
    # bioavailability term.
    # -------------------------------------------------------------------
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-          k12 * central - k21 * peripheral1

    # -------------------------------------------------------------------
    # 4. Observation and residual error. Concentrations are in ug/mL
    # (= mg/L) given doses in mg and volumes in L; the paper reports
    # certepetide concentrations in ng/mL, which is 1000 * Cc.
    # -------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
