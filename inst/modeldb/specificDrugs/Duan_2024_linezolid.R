Duan_2024_linezolid <- function() {
  description <- paste(
    "One-compartment population PK model with first-order elimination for",
    "intravenous linezolid in Chinese premature neonates undergoing",
    "therapeutic drug monitoring (Duan 2024).",
    "Clearance and central volume both scale with body surface area as a",
    "power function referenced to the cohort mean BSA of 0.127 m^2:",
    "V = 0.783 * (BSA/0.127)^1.066 (Equation 1) and",
    "CL = 0.154 * (BSA/0.127)^1.185 (Equation 2, whose printed exponent",
    "1.186 is the rounded form of the Table 4 estimate 1.185 used here).",
    "BSA was the only covariate retained by stepwise",
    "covariate modelling; gestational age, postnatal age, postmenstrual",
    "age, birth weight, current weight, sex, and the hepatic / renal",
    "laboratory panel were all screened and rejected (see",
    "covariatesDataExcluded). Inter-individual variability was reported",
    "only for clearance (omega^2 = 0.132); Equation 1 also carries an",
    "exp(etaVd) term but Table 4 reports no corresponding variance, so",
    "etalvc is encoded as fixed(0). Residual variability is proportional",
    "with an unusually large magnitude (stdev0 = 1.120, i.e. 112%); see",
    "the vignette Assumptions and deviations section before using this",
    "model for stochastic residual-error simulation.",
    sep = " "
  )
  reference <- paste(
    "Duan LF, Li JJ, Shen LR, Chen XL, Yu YX, Yang ZM, Zhang Q, Cai Y,",
    "Li JH, Wu J, Zhao HZ, Xu JH, Feng ZT, Tang L.",
    "Therapeutic drug monitoring of linezolid in Chinese premature",
    "neonates: a population pharmacokinetic analysis and dosage",
    "optimization.",
    "Antimicrob Agents Chemother. 2024;68(11):e01148-24.",
    "doi:10.1128/aac.01148-24.",
    "PMCID PMC11539233.",
    sep = " "
  )
  vignette <- "Duan_2024_linezolid"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  compartmentData <- list(
    central = list(analyte = "linezolid", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (baseline). Duan 2024 Table 1 reports",
        "cohort BSA = 0.13 +/- 0.03 m^2; Equations 1-2 centre the power",
        "model on 0.127 m^2, which is the unrounded cohort mean and the",
        "BSA of the 1.571 kg typical neonate named in the Table 4",
        "footnote. Computed by the authors with the DuBois and DuBois",
        "formula (Methods, Equation 3). Note that Equation 3 as printed",
        "carries the coefficient 0.07184; the published DuBois",
        "coefficient is 0.007184, so the printed equation is off by a",
        "factor of 10 (a typesetting error). BSA is supplied to this",
        "model directly as a covariate column, so the typo does not",
        "affect the model itself; users deriving BSA themselves should",
        "use 0.007184 * WT^0.425 * HT^0.725 with WT in kg and HT in cm.",
        "Monte Carlo simulations in the paper used BSA = 0.11, 0.13, and",
        "0.15 m^2 (Tables 5 and 6).",
        sep = " "
      ),
      source_name        = "BSA"
    )
  )

  # Covariates screened by the stepwise covariate modelling procedure
  # (Duan 2024 Methods, "Pharmacokinetic analysis and model evaluation")
  # but NOT retained in the final model. Documented here so the paper's
  # covariate screen is preserved without declaring covariates that
  # model() never references. Platelet count and the 1-minute / 5-minute
  # Apgar scores were also screened; neither has a canonical continuous
  # covariate-column name in inst/references/covariate-columns.md and
  # neither was retained, so they are recorded in this comment rather
  # than given register entries.
  covariatesDataExcluded <- list(
    GA = list(
      description = "Gestational age",
      units       = "weeks",
      type        = "continuous",
      notes       = "Screened by SCM; not retained. Table 1: 31.00 +/- 2.74 weeks."
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units       = "weeks",
      type        = "continuous",
      notes       = paste(
        "Screened by SCM; not retained. Table 1: 33.16 +/- 2.77 weeks.",
        "Reported in weeks by the paper; the register documents PAGE in",
        "months, so convert before reuse.",
        sep = " "
      )
    ),
    PNA = list(
      description = "Postnatal age",
      units       = "days",
      type        = "continuous",
      notes       = paste(
        "Screened by SCM; not retained. Table 1 median 13.00 days (IQR",
        "8.00-19.25). The Discussion attributes the absence of a PNA",
        "effect to the narrow 8-19 day range in this cohort. Reported in",
        "days by the paper; the register documents PNA in months, so",
        "convert before reuse.",
        sep = " "
      )
    ),
    WT_BIRTH = list(
      description = "Birth weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened by SCM; not retained. Table 1: 1449.44 +/- 496.73 g."
    ),
    WT = list(
      description = "Current body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened by SCM; not retained (BSA was retained instead).",
        "Table 1: 1571.30 +/- 503.81 g. The Methods note that correlated",
        "covariates were screened and only one of a correlated pair was",
        "carried forward, which is why weight drops out in favour of BSA.",
        sep = " "
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened by SCM; not retained. Table 1: 32/54 male (59.26%), so 40.74% female."
    ),
    HGB = list(
      description = "Haemoglobin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened by SCM; not retained. Table 1: 126.43 +/- 20.90 g/L."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened by SCM; not retained. Table 1 median 9.00 U/L (IQR 5.75-17.00)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened by SCM; not retained. Not tabulated in Table 1."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened by SCM; not retained. Table 1 median 81.36 umol/L (IQR 40.93-136.94)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened by SCM; not retained. Table 1 median 31.95 g/L (IQR 29.33-34.43)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened by SCM; not retained. Table 1: 40.55 +/- 12.96 umol/L (Jaffe method)."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = paste(
        "Screened by SCM; not retained. Table 1 median 29.46",
        "mL/min/1.73 m^2 (IQR 24.79-40.41). Derived by the authors with",
        "the Schwartz formula CLcr = K * L / SCR using K = 0.33 for",
        "premature infants, body length L in cm, and serum creatinine",
        "SCR in mg/dL (Methods, Equation 4).",
        sep = " "
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 54L,
    n_studies      = 1L,
    age_range      = "gestational age 31.00 +/- 2.74 weeks; postmenstrual age 33.16 +/- 2.77 weeks",
    age_median     = "postnatal age 13.00 days (IQR 8.00-19.25)",
    weight_range   = "current weight 1571.30 +/- 503.81 g; birth weight 1449.44 +/- 496.73 g",
    weight_median  = "1.571 kg (typical-subject weight named in the Table 4 footnote)",
    sex_female_pct = 40.74,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Premature neonates in a neonatal intensive care unit treated with",
      "intravenous linezolid for late-onset sepsis. Pneumonia in 53/54",
      "(98.15%), bloodstream infection in 48/54 (88.89%), necrotizing",
      "enterocolitis in 5/54 (9.26%), suppurative meningitis in 1/54",
      "(1.85%). Complications included respiratory failure (77.78%),",
      "respiratory distress syndrome (53.70%), asphyxia (27.78%), and",
      "septic shock (11.11%). A pathogen was identified in 31/54",
      "(57.41%), most commonly coagulase-negative Staphylococcus (27).",
      sep = " "
    ),
    dose_range     = paste(
      "Intravenous linezolid (Zyvox) per the manufacturer's",
      "instructions: 10 mg/kg q12h for preterm neonates with gestational",
      "age <34 weeks and postnatal age <7 days; 10 mg/kg q8h for preterm",
      "neonates with gestational age >=34 weeks and for those with",
      "gestational age <34 weeks and postnatal age >=7 days. Treatment",
      "course 7-14 days (median 10 days, IQR 9-11).",
      sep = " "
    ),
    regions        = "China (single centre: Affiliated Suzhou Hospital of Nanjing Medical University, Suzhou, Jiangsu)",
    notes          = paste(
      "Retrospective therapeutic-drug-monitoring cohort, November 2019",
      "to November 2023. 54 infants contributed 84 serum linezolid",
      "concentrations measured by LC-MS/MS (validated 0.5-50 ug/mL).",
      "One or two samples per infant, at least one of which was a",
      "steady-state trough drawn after the fourth maintenance dose and",
      "30 min before the next dose; the remainder were opportunistic",
      "post-fourth-dose samples. Model fitted in Phoenix NLME 8.3 by",
      "FOCE-ELS. Baseline demographics from Duan 2024 Table 1;",
      "concomitant antibiotics in 23/54 (42.59%), most often meropenem",
      "(17). Because the design is trough-dominated with at most two",
      "samples per infant, the paper notes that population predictions",
      "of the final model were poor (Discussion, limitation iii).",
      sep = " "
    )
  )

  ini({
    # Final-model fixed-effect estimates from Duan 2024 Table 4 (Full model,
    # Estimate column) and Equations 1-2. Reference subject: BSA = 0.127 m^2,
    # the cohort-mean BSA used to centre both power terms, corresponding to
    # the 1.571 kg typical neonate of the Table 4 footnote b.
    lvc <- log(0.783); label("Central volume of distribution at BSA = 0.127 m^2 (L)")  # Duan 2024 Table 4 tvV = 0.783 L (RSE 4.749%, 95% CI 0.709-0.857; bootstrap median 0.787, 95% CI 0.712-0.869); Equation 1
    lcl <- log(0.154); label("Clearance at BSA = 0.127 m^2 (L/h)")                     # Duan 2024 Table 4 tvCL = 0.154 L/h (RSE 5.010%, 95% CI 0.139-0.169; bootstrap median 0.154, 95% CI 0.140-0.172); Equation 2. Text: equivalent to 0.098 L/h/kg for the 1.571 kg typical neonate.

    # BSA power exponents. Equation 1: V = 0.783 * (BSA/0.127)^1.066; Equation
    # 2: CL = 0.154 * (BSA/0.127)^1.186. Table 4 reports the CL exponent as
    # 1.185 with its RSE and confidence intervals; the printed Equation 2
    # rounds it to 1.186. The Table 4 value is used here because it is the
    # estimate column that carries the uncertainty and bootstrap summaries.
    e_bsa_vc <- 1.066; label("Power exponent of (BSA/0.127) on Vd (unitless)")  # Duan 2024 Table 4 theta_BSA-V = 1.066 (RSE 30.750%, 95% CI 0.413-1.720; bootstrap median 1.100, 95% CI 0.311-1.915); Equation 1
    e_bsa_cl <- 1.185; label("Power exponent of (BSA/0.127) on CL (unitless)")  # Duan 2024 Table 4 theta_BSA-CL = 1.185 (RSE 21.774%, 95% CI 0.671-1.700; bootstrap median 1.212, 95% CI 0.555-1.771); Equation 2 prints the rounded 1.186

    # Inter-individual variability. Duan 2024 describes IIV with an
    # exponential random-effects model (Results, "Development of the
    # population pharmacokinetics model"). Table 4 reports a single IIV
    # term, omega^2 CL = 0.132 (shrinkage 5.838%), i.e. an omega of 0.363
    # on the log scale and a lognormal CV of
    # sqrt(exp(0.132) - 1) = 37.6%.
    etalcl ~ 0.132  # Duan 2024 Table 4 omega^2 CL = 0.132 (RSE 18.58%, 95% CI 0.0914-0.172; bootstrap median 0.129, 95% CI 0.0886-0.170); shrinkage 5.838%

    # Equation 1 writes V = 0.783 * (BSA/0.127)^1.066 * exp(etaVd), i.e. it
    # declares an IIV term on volume, but Table 4 reports no omega^2 V and no
    # supplement exists. Encoded as fixed(0) so the published structural
    # declaration is preserved without inventing a variance. Precedent:
    # Azechi_2024_tamibarotene_pediatric.R, Chi_2018_propofol.R.
    etalvc ~ fixed(0)  # Duan 2024 Equation 1 declares exp(etaVd); Table 4 reports no corresponding variance

    # Residual error. Duan 2024 Results: "residual variability was fitted
    # with a proportional residual error model"; Methods state that additive,
    # proportional, and mixed models were tested and the proportional model
    # was selected. Table 4 reports the single residual term as
    # "stdev0 = 1.120" with footnote "stdev0, standard deviation" and, unlike
    # every other row of that table, no unit in the "Parameter (unit)"
    # column -- consistent with the dimensionless SD of the proportional
    # epsilon in Phoenix NLME's Cobs = C * (1 + eps) parameterisation.
    # 1.120 therefore encodes a 112% proportional residual SD. This is very
    # large (it implies ~19% of residual draws are negative), which is
    # consistent with the authors' own statement that population predictions
    # of the final model were poor on this trough-dominated design. See the
    # vignette Assumptions and deviations section: validation in the vignette
    # is performed on the individual-prediction scale, not by adding this
    # residual error.
    propSd <- 1.120; label("Proportional residual SD (fraction)")  # Duan 2024 Table 4 stdev0 = 1.120 (RSE 13.183%, 95% CI 0.826-1.415; bootstrap median 1.052, 95% CI 0.635-1.341)
  })

  model({
    # BSA centering value: 0.127 m^2, the cohort-mean BSA that Duan 2024
    # Equations 1-2 divide by (Table 1 reports the rounded 0.13 +/- 0.03).
    bsa_norm <- BSA / 0.127

    # Individual PK parameters, Duan 2024 Equations 1-2.
    vc <- exp(lvc + etalvc) * bsa_norm^e_bsa_vc
    cl <- exp(lcl + etalcl) * bsa_norm^e_bsa_cl

    kel <- cl / vc

    # One-compartment disposition with first-order elimination. Linezolid was
    # given intravenously, so the dose enters `central` directly; there is no
    # absorption compartment and no bioavailability term.
    d/dt(central) <- -kel * central

    # Doses in mg and vc in L give central/vc in mg/L, which equals the
    # ug/mL used throughout Duan 2024 (target trough 2-8 ug/mL, assay range
    # 0.5-50 ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
