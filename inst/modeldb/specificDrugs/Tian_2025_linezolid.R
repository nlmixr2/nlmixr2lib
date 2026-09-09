Tian_2025_linezolid <- function() {
  description <- paste(
    "One-compartment population PK model with first-order elimination for",
    "intravenous linezolid in Chinese children (0-16 years) with confirmed",
    "or suspected bacterial infections (Tian 2025).",
    "Clearance carries two power covariates referenced to the cohort",
    "medians, CL = 1.91 * (WT/12.5)^0.696 * (eGFR/190.1)^0.291 (Table 3),",
    "where eGFR is the Schwartz-formula estimated glomerular filtration",
    "rate; the cohort median eGFR of 190.1 mL/min/1.73 m^2 places most of",
    "these children in augmented renal clearance. Central volume is a",
    "single typical value with no covariate, V = 10.5 L. Weight and eGFR",
    "were the only covariates retained by stepwise forward selection /",
    "backward elimination; sex, age, height, AST, ALT, albumin, total and",
    "direct bilirubin, serum creatinine, and blood urea nitrogen were all",
    "screened and rejected (see covariatesDataExcluded). Inter-individual",
    "variability is exponential on both CL and V, and residual",
    "variability is exponential (log-normal), so the model is encoded with",
    "lnorm() rather than a proportional residual.",
    sep = " "
  )
  reference <- paste(
    "Tian X, Jiang T, Dong L, Zhang X, Jiao W, Liu G, Li Q, Bi J, You D,",
    "Cao L, Guo W, Jin Z, Zhang Q, Xu Y, Zhao W, Qi H, Zheng Y, Shen A.",
    "Population pharmacokinetics and clinical assessment of linezolid in",
    "pediatric bacterial infections.",
    "Antimicrob Agents Chemother. 2025;69(5):e01299-24.",
    "doi:10.1128/aac.01299-24.",
    "PMCID PMC12057362.",
    sep = " "
  )
  vignette <- "Tian_2025_linezolid"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  compartmentData <- list(
    central = list(analyte = "linezolid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (baseline). Tian 2025 Table 1 reports a",
        "median of 12.5 kg (5th-95th percentile 5.3-45.9 kg), and the",
        "Table 3 footnote names 12.5 kg as the median weight used to",
        "centre the power term on clearance. Weight entered the model in",
        "forward selection with an OFV drop of 77.376 points and was",
        "retained through backward elimination (Results, 'Model",
        "building'). Volume of distribution carries no weight covariate",
        "in the final model.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (Schwartz formula)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (baseline). Computed by the authors with",
        "the Schwartz formula (Methods, 'Study design'), so this is a",
        "creatinine-based, BSA-normalized estimate rather than a measured",
        "clearance. Tian 2025 Table 1 reports a median of 190.1",
        "mL/min/1.73 m^2 (5th-95th percentile 97.0-350.2), and the Table 3",
        "footnote names 190.1 as the median eGFR used to centre the power",
        "term on clearance. That median is far above the 90-130",
        "mL/min/1.73 m^2 the paper defines as normal renal function, so",
        "the reference subject of this model is an augmented-renal-",
        "clearance child; users simulating a normal-renal-function child",
        "should expect CL below the typical value. eGFR entered the model",
        "in forward selection with an OFV drop of 21.119 points and was",
        "retained through backward elimination (Results, 'Model",
        "building').",
        sep = " "
      ),
      source_name        = "eGFR"
    )
  )

  # Covariates screened by the stepwise forward-inclusion / backward-elimination
  # procedure (Tian 2025 Methods, "Population pharmacokinetic modeling of
  # linezolid"; Results, "Model building") but NOT retained in the final model.
  # The paper states that of the 12 screened covariates only WT and eGFR
  # produced a significant OFV drop; no point estimate is reported for any of
  # the rejected covariates. Documented here so the paper's covariate screen is
  # preserved without declaring covariates that model() never references.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained. Table 1: 48/80 male, 32/80 female (40.0% female)."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened; not retained. Table 1 median 3.3 years (5th-95th 0.1-12.6)."
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened; not retained. Table 1 median 99.5 cm (5th-95th 55.1-159.8)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened; not retained. Table 1 median 26.0 U/L (5th-95th 12-117)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened; not retained. Table 1 median 15.3 U/L (5th-95th 5.4-170.0)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened; not retained. Table 1 median 34.0 g/L (5th-95th 25.5-43.5)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened; not retained. Table 1 median 6.5 umol/L (5th-95th 2.4-37.1). Abbreviated TB in the Results text and TBIL in Table 1."
    ),
    DBIL = list(
      description = "Direct (conjugated) bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened; not retained. Table 1 median 1.7 umol/L (5th-95th 0.5-15.0). Abbreviated DB in the Results text and DBIL in Table 1."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Screened; not retained. Table 1 median 24.0 umol/L (5th-95th",
        "9.0-48.0). Serum creatinine is the input to the Schwartz eGFR",
        "that WAS retained, so the two are strongly correlated and only",
        "the derived eGFR survives the covariate screen.",
        sep = " "
      )
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened; not retained. Table 1 median 2.7 mmol/L (5th-95th 1.0-5.9)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 80L,
    n_studies      = 1L,
    age_range      = "0-16 years by inclusion criterion; observed median 3.3 years (5th-95th 0.1-12.6)",
    weight_range   = "median 12.5 kg (5th-95th 5.3-45.9)",
    height_range   = "median 99.5 cm (5th-95th 55.1-159.8)",
    sex_female_pct = 40.0,
    race_ethnicity = c(Asian = 100),
    renal_function = paste(
      "Schwartz eGFR median 190.1 mL/min/1.73 m^2 (5th-95th 97.0-350.2).",
      "Most children were in augmented renal clearance, which the paper",
      "defines as eGFR > 130 mL/min/1.73 m^2; the Discussion reports that",
      "CL and V were both significantly higher in the ARC subgroup than",
      "in children with normal renal function (90 <= eGFR <= 130), while",
      "AUC/MIC did not differ significantly between the two.",
      sep = " "
    ),
    disease_state  = paste(
      "Children with confirmed or suspected bacterial infection treated",
      "with intravenous linezolid for more than 48 h. In the 67-patient",
      "efficacy and safety subset (Table 2) the commonest presentations",
      "were pulmonary infection (32, 47.8%), central nervous system",
      "infection (21, 31.3%), sepsis (21, 31.3%), osteomyelitis (9,",
      "13.4%) and skin / soft tissue infection (8, 11.9%). Cultured",
      "pathogens included Staphylococcus aureus (13), Streptococcus",
      "pneumoniae (5) and methicillin-resistant Staphylococcus aureus",
      "(4). Children receiving levofloxacin were excluded because it is",
      "the internal standard of the linezolid assay.",
      sep = " "
    ),
    dose_range     = paste(
      "Intravenous linezolid (Zyvox, Pfizer) at the label standard dose:",
      "10 mg/kg q8h for children under 12 years (73/80 subjects, median",
      "10.0 mg/kg/dose, 5th-95th 9.4-10.5) and 600 mg q12h for children",
      "aged 12 years and above (7/80 subjects).",
      sep = " "
    ),
    regions        = paste(
      "China (six centres: Beijing Children's Hospital of Capital Medical",
      "University; Baoding Hospital of Beijing Children's Hospital;",
      "Hebei Children's Hospital; Children's Hospital of the Capital",
      "Institute of Pediatrics; Henan Children's Hospital; Tianjin",
      "Children's Hospital)",
      sep = " "
    ),
    notes          = paste(
      "Prospective multi-centre study, March 2021 to June 2022, Chinese",
      "Clinical Trial Registry ChiCTR2200061207. 80 children contributed",
      "157 plasma linezolid concentrations (range 0.25-33.67 ug/mL),",
      "median 2.0 samples per patient, drawn opportunistically once",
      "steady state had been reached (at least 48 h of treatment).",
      "Concentrations measured by HPLC-UV at 254 nm with levofloxacin as",
      "internal standard; lower limit of quantitation 0.25 ug/mL, linear",
      "range 0.25-50 ug/mL. Model fitted in NONMEM 7.4 with FOCE-I and",
      "evaluated by 1000-replicate nonparametric bootstrap, NPDE, and VPC",
      "(NPDE mean 0.0652, variance 1.08, global adjusted P = 0.855).",
      "A separate efficacy and safety analysis in 67 of the 80 children",
      "used empirical Bayes estimates from this model; it reports",
      "descriptive AUC/MIC and Cmin summaries and fits no additional",
      "exposure-response model, so nothing beyond this PK model is",
      "extractable from the paper.",
      sep = " "
    )
  )

  ini({
    # Final-model fixed-effect estimates from Tian 2025 Table 3 ("Final model",
    # Estimate column) and the CL covariate equation printed in the same table.
    # Reference subject: WT = 12.5 kg and eGFR = 190.1 mL/min/1.73 m^2, the
    # cohort medians named in the Table 3 footnote.
    lcl <- log(1.91); label("Clearance at WT = 12.5 kg and eGFR = 190.1 mL/min/1.73 m^2 (L/h)")  # Tian 2025 Table 3 theta1 = 1.91 L/h (RSE 5.1%; bootstrap median 1.88). Consistent with the Discussion's 0.15 +/- 0.06 L/h/kg: 1.91 / 12.5 = 0.153 L/h/kg.
    lvc <- log(10.5); label("Central volume of distribution (L)")                                # Tian 2025 Table 3 theta2 = 10.5 L (RSE 11.8%; bootstrap median 10.5). V = theta2 with no covariate. Consistent with the Discussion's 0.77 +/- 0.50 L/kg: 10.5 / 12.5 = 0.84 L/kg.

    # Power exponents of the two retained clearance covariates.
    # Table 3: CL = theta1 * (WT/12.5)^theta3 * (eGFR/190.1)^theta4.
    e_wt_cl   <- 0.696; label("Power exponent of (WT/12.5) on CL (unitless)")     # Tian 2025 Table 3 theta3 = 0.696 (RSE 10.6%; bootstrap median 0.701, 5th-95th 0.348-1.09). Estimated, not fixed at an allometric 0.75.
    e_crcl_cl <- 0.291; label("Power exponent of (eGFR/190.1) on CL (unitless)")  # Tian 2025 Table 3 theta4 = 0.291 (RSE 36.1%; bootstrap median 0.309, 5th-95th 0.01-0.983)

    # Inter-individual variability. Tian 2025 Methods: "Interindividual
    # variability of the PK parameters was estimated using an exponential
    # model", theta_i = theta_mean * exp(eta_i). Table 3 reports the two IIV
    # terms in a column headed "Inter-individual variability (%)" as 33.02 and
    # 82.28. Those percentages are read here as omega on the log scale (i.e.
    # 100 * sqrt(omega^2)), not as exact log-normal CVs, because the SAME
    # percent column of the same table also carries the residual-variability
    # row, and NONMEM / PsN print that quantity as 100 * sqrt(sigma^2). Reading
    # one row of a column on the SD scale and its neighbours on a
    # CV-transformed scale is not a convention any of that toolchain uses. The
    # competing reading (omega^2 = log(1 + CV^2)) would give 0.1035 and 0.5170;
    # see the vignette Assumptions and deviations section, which states the
    # consequence of the alternative reading.
    etalcl ~ 0.10903  # Tian 2025 Table 3 IIV CL = 33.02% (RSE 37.7%; bootstrap median 33.02, 5th-95th 0.71-168.39) -> omega = 0.3302, omega^2 = 0.3302^2 = 0.10903; lognormal CV sqrt(exp(0.10903) - 1) = 33.9%
    etalvc ~ 0.67700  # Tian 2025 Table 3 IIV V = 82.28% (RSE 33.8%; bootstrap median 81.02, 5th-95th 0.55-132.16) -> omega = 0.8228, omega^2 = 0.8228^2 = 0.67700; lognormal CV sqrt(exp(0.67700) - 1) = 98.4%

    # Residual error. Tian 2025 Methods: "Residual variability was estimated
    # using an exponential model, addition model, and mixed model,
    # respectively. The optimal residual model was determined by analyzing the
    # OFV", and Results, "Model building": "exponential models were used for
    # both interindividual variation and residual variation". An exponential
    # residual, Y = IPRED * exp(eps), is nlmixr2's lnorm() residual, so the
    # 26.23% of Table 3 is the additive SD on the log-transformed
    # concentration scale, 0.2623 -- the same scale reading applied to the two
    # IIV rows above.
    expSd <- 0.2623; label("Additive residual SD on the log-transformed concentration scale (log-normal)")  # Tian 2025 Table 3 residual variability = 26.23% (RSE 69.6%; bootstrap median 28.54, 5th-95th 13.52-106.25)
  })

  model({
    # Covariate centering values, Tian 2025 Table 3 footnote: "In our
    # population, 12.5 kg and 190.1 mL/min/1.73 m^2 are the median WT and
    # eGFR values, respectively."
    wt_norm   <- WT / 12.5
    crcl_norm <- CRCL / 190.1

    # Individual PK parameters, Tian 2025 Table 3.
    # CL = theta1 * (WT/12.5)^theta3 * (eGFR/190.1)^theta4; V = theta2.
    cl <- exp(lcl + etalcl) * wt_norm^e_wt_cl * crcl_norm^e_crcl_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One-compartment disposition with first-order elimination. Linezolid was
    # given by intravenous infusion, so the dose enters `central` directly;
    # there is no absorption compartment and no bioavailability term.
    d/dt(central) <- -kel * central

    # Doses in mg with vc in L give central/vc in mg/L, which equals the
    # ug/mL used throughout Tian 2025 (assay linear range 0.25-50 ug/mL,
    # observed concentrations 0.25-33.67 ug/mL, safety threshold Cmin =
    # 7 ug/mL).
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
