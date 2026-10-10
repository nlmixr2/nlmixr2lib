Zhou_2022_simurosertib <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model for oral simurosertib",
    "(TAK-931, a cell division cycle 7 kinase inhibitor) in 198 adults with",
    "advanced solid tumors, from Zhou 2022. Absorption is a chain of two",
    "transit compartments in which the absorption rate constant equals the",
    "transit rate constant (ktr = 2 / MTT); elimination is first-order",
    "linear. Creatinine clearance (Cockcroft-Gault) and body weight enter",
    "apparent clearance as power functions centred on 90.45 mL/min and",
    "65.95 kg, and body weight also scales the apparent central volume and",
    "intercompartmental clearance. Between-subject variability is carried",
    "on apparent clearance and the mean transit time; residual error is",
    "additive on log-transformed concentrations.",
    sep = " "
  )
  reference <- paste(
    "Zhou X, Ouerdani A, Diderichsen PM, Gupta N.",
    "Population Pharmacokinetics of TAK-931, a Cell Division Cycle 7 Kinase",
    "Inhibitor, in Patients With Advanced Solid Tumors.",
    "J Clin Pharmacol. 2022;62(3):422-433.",
    "doi:10.1002/jcph.1974. PMC9297904.",
    "Open Access under CC BY-NC 4.0.",
    "Parameter estimates are in Table 3 and Equation 7; the model schema is",
    "Supplemental Figure S1 and the simulated day-1 / day-14 exposure",
    "summaries used for validation are Supplemental Table S2.",
    sep = " "
  )
  vignette <- "Zhou_2022_simurosertib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on CL/F (exponent 0.484), Vc/F (0.867) and Q/F",
        "(0.938), centred on 65.95 kg (Equation 7). The Table 2 cohort",
        "median is 65.8 kg (range 29.8-127 kg); the 65.95 kg printed in",
        "Equation 7 and the Figure 5 legend is the analysis-data-set median",
        "used as the centring value, and is the value encoded here.",
        sep = " "
      ),
      source_name = "WGT"
    ),
    CRCL = list(
      description = paste(
        "Baseline creatinine clearance estimated by the Cockcroft-Gault",
        "equation on total body weight, in raw mL/min (NOT BSA-normalized).",
        sep = " "
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on CL/F (exponent 0.325), centred on 90.45 mL/min",
        "(Equation 7; Methods, 'Model-Based Simulations', 'CrCL of 90.45",
        "mL/min (median value)'). Cockcroft-Gault as printed in Equation 3:",
        "(140 - age) * weight / (72 * serum creatinine in mg/dL), times 0.85",
        "for women. Table 2 cohort median 89.9 mL/min (range 35-204 mL/min);",
        "the model was not informed below 35 mL/min.",
        sep = " "
      ),
      source_name = "CrCL"
    )
  )

  # Covariates the source screened and did not retain in the final model.
  # Documentation only -- none is referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the univariate step (Methods, 'Covariate Model Development') and not retained; the Abstract and Conclusions state age (36-88 years) had no impact on CL/F. No coefficient is reported."
    ),
    BSA = list(
      description = "Baseline body surface area.",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a body-size covariate and not retained; body weight was carried instead (Results, 'Covariate Model Development'). No coefficient is reported."
    ),
    BMI = list(
      description = "Baseline body mass index.",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a body-size covariate and not retained; body weight was carried instead. No coefficient is reported."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Not identified as a covariate on CL/F (Results, 'Covariate Model Development'); the apparent trend with sex was attributed to its correlation with body weight. No coefficient is reported."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Race was evaluated as White vs Asian vs other and was not identified as a covariate on CL/F (Results; Figure 1A). No coefficient is reported."
    ),
    TBILI = list(
      description = "Baseline total bilirubin.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the univariate step and not retained. Table 2 cohort median 8.6 umol/L (range 1.7-24.5). No coefficient is reported."
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic impairment indicator (NCI ODWG).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal hepatic function)",
      notes = "Mild hepatic impairment (TB <= ULN and AST > ULN, or ULN < TB <= 1.5 x ULN with any AST; 19.7% of the cohort) was not identified as a covariate on CL/F (Results; Figure 1). No coefficient is reported."
    ),
    WHO_PS = list(
      description = "Baseline ECOG performance status.",
      units = "(integer score)",
      type = "count",
      reference_category = NULL,
      notes = "ECOG 0 (55.6%) or 1 (44.4%); screened and not identified as a covariate on CL/F. No coefficient is reported."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "simurosertib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "simurosertib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "simurosertib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "simurosertib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "simurosertib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 198L,
    n_studies = 3L,
    n_observations = 2678L,
    age_range = "36-88 years",
    age_median = "61 years",
    weight_range = "29.8-127 kg",
    weight_median = "65.8 kg",
    crcl_range = "35-204 mL/min (Cockcroft-Gault)",
    crcl_median = "89.9 mL/min",
    sex_female_pct = 60,
    race_ethnicity = c(Asian = 55.1, White = 33.3, `Black or African American` = 4, Other = 7.6),
    disease_state = paste(
      "Adults with advanced nonhematologic solid tumors, including metastatic",
      "pancreatic and colorectal cancer. Renal function: 50.5% normal, 32.8%",
      "mild and 17.7% moderate impairment; hepatic function: 80.3% normal and",
      "19.7% mild impairment; ECOG 0 (55.6%) or 1 (44.4%).",
      sep = " "
    ),
    dose_range = paste(
      "Oral simurosertib 20-150 mg once daily on several on/off schedules",
      "(14 days on / 7 off, 7 on / 7 off, 21 days continuous, 2 on / 5 off)",
      "as powder-in-capsule, plus 80 mg single doses of capsule and tablet in",
      "a relative-bioavailability crossover. Taken fasted.",
      sep = " "
    ),
    regions = "Japan (first-in-human phase 1 TAK-931-1002) and Western sites (phase 2 TAK-931-2001).",
    notes = paste(
      "Pooled from studies TAK-931-1002 (phase 1, n = 80; NCT02699749),",
      "TAK-931-1003 (phase 1 relative bioavailability, n = 20; NCT03708211)",
      "and TAK-931-2001 (phase 2, n = 98; NCT03261947) (Table 1).",
      "Demographics from Table 2. Capsule and tablet were bioequivalent, so",
      "formulation was not tested as a covariate.",
      sep = " "
    )
  )

  ini({
    # Structural parameters: Zhou 2022 Table 3 ('Final Model Parameter
    # Estimates'), for the reference patient with CrCL 90.45 mL/min and
    # body weight 65.95 kg (Equation 7). All are apparent (relative to F).
    lcl <- log(38.0)
    label("Apparent clearance CL/F at CrCL 90.45 mL/min and 65.95 kg (L/h)") # Table 3 'CL/F, L/h' 38.0 (RSE 0.6%); Equation 7 '38'
    lvc <- log(194)
    label("Apparent central volume Vc/F at 65.95 kg (L)") # Table 3 'Vc/F, L' 194 (RSE 0.5%); Equation 7
    lq <- log(7.71)
    label("Apparent intercompartmental clearance Q/F at 65.95 kg (L/h)") # Table 3 'Q/F, L/h' 7.71 (RSE 2.4%); Equation 7
    lvp <- log(140)
    label("Apparent peripheral volume Vp/F (L)") # Table 3 'Vp/F, L' 140 (RSE 1.6%)
    lmtt <- log(0.756)
    label("Mean transit time through the two transit steps (h)") # Table 3 'MTT, h' 0.756 (RSE 17.3%)

    # Covariate exponents: Table 3 and Equation 7 (power functions centred
    # on the analysis-data-set medians, Equation 4).
    e_crcl_cl <- 0.325
    label("Power exponent of CrCL/90.45 on CL/F (unitless)") # Table 3 'Covariate exponent for creatinine clearance on CL/F' 0.325 (RSE 15%)
    e_wt_cl <- 0.484
    label("Power exponent of WT/65.95 on CL/F (unitless)") # Table 3 'Covariate exponent for body weight on CL/F' 0.484 (RSE 21.8%)
    e_wt_vc <- 0.867
    label("Power exponent of WT/65.95 on Vc/F (unitless)") # Table 3 'Covariate exponent for body weight on Vc/F' 0.867 (RSE 11.4%)
    e_wt_q <- 0.938
    label("Power exponent of WT/65.95 on Q/F (unitless)") # Table 3 'Covariate exponent for body weight on Q/F' 0.938 (RSE 15.7%)

    # IIV: Table 3 reports %CV; the footnote converts with
    # CV = sqrt(exp(omega^2) - 1), so omega^2 = log(CV^2 + 1).
    etalcl ~ 0.050245 # Table 3 'IIV on CL/F, %CV' 22.7 -> log(0.227^2 + 1)
    etalmtt ~ 0.357909 # Table 3 'IIV on MTT, %CV' 65.6 -> log(0.656^2 + 1)

    # Residual error: additive on log-transformed concentrations
    # (Methods Equation 2, 'log-transform both sides').
    expSd <- 0.499
    label("Additive residual error SD on the log scale (unitless)") # Table 3 'Additive residual error (standard deviation) in log scale' 0.499 (RSE 4%)
  })
  model({
    # Equation 7 covariate model (references 90.45 mL/min and 65.95 kg)
    cl <- exp(lcl + etalcl) * (CRCL / 90.45)^e_crcl_cl * (WT / 65.95)^e_wt_cl
    vc <- exp(lvc) * (WT / 65.95)^e_wt_vc
    q <- exp(lq) * (WT / 65.95)^e_wt_q
    vp <- exp(lvp)
    mtt <- exp(lmtt + etalmtt)

    # Supplemental Figure S1: dose -> transit 1 -> transit 2 -> central,
    # every transfer governed by MTT and NTR = 2. The absorption rate
    # constant out of transit 2 was set equal to the transit rate constant
    # (Results, 'Base Model Development'). The paper does not print the
    # ktr-MTT relation. ktr = NTR / MTT = 2 / MTT is used: MTT is the mean
    # time through the NTR transit steps, and the ka step that was later set
    # equal to ktr comes on top. This reading reproduces the paper's own
    # simulated day-1 Cmax summary (Supplemental Table S2: median 171.1,
    # 2.5th percentile 96.5 ng/mL); the alternative ktr = (NTR + 1) / MTT
    # overpredicts the median by ~15% and the 2.5th percentile by ~30%.
    # See the vignette 'Assumptions and deviations'.
    ktr <- 2 / mtt

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(central) <- ktr * transit2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L; x 1000 gives ng/mL
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
