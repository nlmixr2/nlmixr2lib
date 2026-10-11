Tan_2023_luvometinib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral luvometinib (FCN-159), a",
    "MEK1/2 inhibitor, in 45 Chinese adults with advanced NRAS-aberrant",
    "melanoma or neurofibromatosis type 1 (NF1) pooled from two phase I",
    "studies. First-order absorption with an absorption lag time and linear",
    "elimination. Body surface area enters all four disposition parameters",
    "as an estimated power function centered at 1.66 m^2; total serum",
    "protein scales the central volume, female sex enlarges the peripheral",
    "volume, and relative bioavailability is 1.29-fold higher in NF1 than",
    "in melanoma. Log-normal between-subject variability on CL/F, Vc/F, Q/F",
    "and Vp/F with a proportional residual error.",
    sep = " "
  )
  reference <- paste(
    "Tan Y., Cui A., Qian L., Li C., Wu Z., Yang Y., Han P., Huang X.,",
    "Diao L. (2023). Population pharmacokinetics of FCN-159, a MEK1/2",
    "inhibitor, in adult patients with advanced melanoma and",
    "neurofibromatosis type 1 (NF1) and model informed dosing",
    "recommendations for NF1 pediatrics.",
    "Frontiers in Pharmacology 14:1101991.",
    "doi:10.3389/fphar.2023.1101991.",
    sep = " "
  )
  vignette <- "Tan_2023_luvometinib"

  # Tan 2023 Methods 2.2: the LC-MS/MS assay range is '0.2-200 ng/ml'. Doses
  # are in mg and the disposition parameters in L and L/hr (Table 2), so an
  # amount over a volume is mg/L; model() multiplies by 1000 to report Cc in
  # ng/mL (1 mg/L = 1000 ng/mL).
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Every disposition parameter is APPARENT (Table 2 heads them CL/F, VC/F,
  # Q/F, VP/F), so the state amounts are apparent amounts relative to the
  # melanoma bioavailability, which the paper fixes at 1 (Table 2 footnote).
  compartmentData <- list(
    depot = list(
      analyte = "luvometinib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "luvometinib",
      units = "mg (apparent, relative to melanoma bioavailability)",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "luvometinib",
      units = "mg (apparent, relative to melanoma bioavailability)",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area, the allometric size covariate on CL/F, Vc/F, Q/F and Vp/F",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power function (BSA/1.66)^theta on all four disposition parameters",
        "(Tan 2023 Table 1). Methods 2.4.2: 'BSA to the median BSA value of",
        "1.66 m 2'; Table S2 confirms the pooled median 1.66 m^2 (range",
        "1.33-2.16). BSA was included 'by default without statistical test'",
        "as the allometric term, with exponents ESTIMATED rather than fixed",
        "(Table 2 reports an RSE for each). The BSA formula is not stated in",
        "the paper. Baseline value."
      ),
      source_name = "BSA"
    ),
    TPRO = list(
      description = "Total serum (plasma) protein, a power covariate on the apparent central volume Vc/F",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power function (TP/68.95)^-3.1 on Vc/F (Tan 2023 Table 1 and",
        "Table 2 'VC/F - TP'). The centering value 68.95 g/L is the printed",
        "Table 1 equation constant; Table S2 lists a pooled median of 71.75",
        "g/L (range 59-82), so 68.95 is not the Table S2 median. The",
        "equation value is used because the exponent was estimated against",
        "it; Results 3.3 also uses 68.95 g/L as the forest-plot reference",
        "('a patient with total plasma protein of 76.89 g/L would be",
        "expected to have 30% higher Cmax,ss relative to a patient with",
        "total plasma protein of 68.95 g/L'). Baseline value."
      ),
      source_name = "TP"
    ),
    SEXF = list(
      description = "Female sex indicator, a fractional-change covariate on the apparent peripheral volume Vp/F",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Tan 2023 Table 1 footnote: 'sex (SEXi) was categorized into male",
        "(SEXi = 0) and female (SEXi = 1)', and Vp/F carries",
        "(1 + theta13 x SEXi (if female)). The paper's SEX column already has",
        "the canonical SEXF orientation, so no recoding is needed. Table S2:",
        "26 male of 45 (19 female)."
      ),
      source_name = "SEX"
    ),
    TUMTP_NF1 = list(
      description = "Neurofibromatosis type 1 indicator (vs advanced melanoma), a covariate on relative oral bioavailability",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (advanced melanoma; bioavailability fixed at 1)",
      notes = paste(
        "Tan 2023 Table 1 footnote: 'cancer type (TYPEi) was categorized",
        "into melanoma (TYPEi = 0) and NF1 (TYPEi = 1)'. F1 = 1 for melanoma",
        "and theta7 = 1.29 for NF1 (Table 2 footnote: 'F1 is the relative",
        "bioavailability for NF1, where F1 is 1 for melanoma'). Disease and",
        "study are fully confounded: all 33 FCN-159-001 subjects had",
        "melanoma and all 12 FCN-159-002 subjects had NF1 (Table S2), so",
        "the effect is equally a study effect; the Discussion attributes it",
        "to bioavailability because Cmax and AUC were both higher in NF1."
      ),
      source_name = "TYPE"
    )
  )

  # Covariates Tan 2023 screened (Methods 2.4.2) but did not retain in the
  # final model.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Retained on CL/F in the full covariate model with an estimate of",
        "-0.00972 (RSE 23.70%) but removed from the final model as too small",
        "to matter (Results 3.3); its functional form is not stated. Table S2",
        "median 49 years (range 20-71)."
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Retained on CL/F in the full covariate model with an estimate of",
        "-0.0134 (RSE 54%) but removed from the final model (Results 3.3).",
        "Removing age and albumin raised IIV on CL/F from 17.2% to 19% and",
        "OFV from 2995.8 to 3014.724. Table S2 median 45.0 g/L (35.2-53)."
      )
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Tested as the alternative allometric size covariate (normalized to",
        "60 kg); BSA was selected instead (Discussion: the two gave 'quite",
        "similar goodness of fit'). Table S2 median 63 kg (43-106)."
      )
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Screened on CL/F (Methods 2.4.2); not retained. Table S2 median 104 mL/min (49.1-224)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL/F; not retained. Table S2 median 19 U/L (8-49.7)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL/F; not retained. Table S2 median 14.8 U/L (4-35.4)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened on CL/F; not retained. Table S2 median 12.2 umol/L (4.6-29.8)."
    ),
    RBC = list(
      description = "Red blood cell count",
      units = "10^12/L",
      type = "continuous",
      notes = "Screened on the volumes of distribution; not retained. Table S2 median 4.56 x 10^12/L (3.25-6.09)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 45L,
    n_studies = 2L,
    n_observations = "1030 plasma concentrations; 60 of 1090 samples (5.5%) were below the 0.2 ng/mL LLOQ and excluded (Results 3.1)",
    age_range = "20-71 years",
    age_median = "49 years",
    weight_range = "43-106 kg",
    weight_median = "63 kg",
    bsa_range = "1.33-2.16 m^2",
    bsa_median = "1.66 m^2",
    total_protein = "median 71.75 g/L (range 59-82)",
    sex_female_pct = 42.2,
    race_ethnicity = "Chinese (both studies were conducted in China with Chinese patients)",
    disease_state = "advanced melanoma with aberrant or mutant NRAS (FCN-159-001, n = 33) and neurofibromatosis type 1 (FCN-159-002 adult cohort, n = 12)",
    dose_range = "0.2-15 mg orally once daily in the fasting state (melanoma: 0.2-4 mg for 21 days and 6-15 mg for 28 days of 28-day cycles; NF1: 4-12 mg continuous)",
    regions = "China",
    notes = paste(
      "Demographics are Tan 2023 Supplementary Table S2 (overall and by",
      "study); study design and PK sampling are Supplementary Table S1. The",
      "melanoma study opened with a single-dose PK run-in that the",
      "Discussion credits with characterizing the two-compartment",
      "disposition. Table S2 prints the sex and diagnosis percentages",
      "against a denominator of 48 instead of 45 (e.g. 26 male shown as",
      "54.17%); sex_female_pct above is 19/45. Pediatric NF1 exposure was",
      "projected by simulation with the same model (Results 3.5); no",
      "pediatric data entered the fit."
    )
  )

  ini({
    # Final population estimates, Tan 2023 Table 2 ('Final population
    # pharmacokinetic model parameters for FCN-159'), with the covariate
    # equations of Table 1. Typical values are for a male melanoma patient
    # with BSA 1.66 m^2 and total protein 68.95 g/L.
    lka <- log(0.5); label("First-order absorption rate constant (1/h)") # Table 2 'KA (1/hr) = 0.5' (RSE 8%)
    ltlag <- log(0.211); label("Absorption lag time (h)") # Table 2 'ALAG1 (hr) = 0.211' (RSE 2.10%)
    lcl <- log(13.2); label("Apparent clearance CL/F (L/h)") # Table 2 'CL/F (L/hr) = 13.2' (RSE 3.80%)
    lvc <- log(48.7); label("Apparent central volume Vc/F (L)") # Table 2 'VC/F (L) = 48.7' (RSE 13.30%)
    lq <- log(35.1); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 'Q/F (L/hr) = 35.1' (RSE 7.70%)
    lvp <- log(314); label("Apparent peripheral volume Vp/F (L)") # Table 2 'VP/F (L) = 314' (RSE 6.10%)

    # Power exponents on (BSA/1.66), all estimated (Table 1 theta8-theta11).
    e_bsa_cl <- 1.11; label("Power exponent of BSA on CL/F (unitless)") # Table 2 'CL/F - BSA = 1.11' (RSE 26.70%)
    e_bsa_vc <- 2.41; label("Power exponent of BSA on Vc/F (unitless)") # Table 2 'VC/F - BSA = 2.41' (RSE 43.60%)
    e_bsa_q <- 2.20; label("Power exponent of BSA on Q/F (unitless)") # Table 2 'Q/F - BSA = 2.20' (RSE 22.10%)
    e_bsa_vp <- 3.84; label("Power exponent of BSA on Vp/F (unitless)") # Table 2 'VP/F - BSA = 3.84' (RSE 13.30%)

    e_tpro_vc <- -3.1; label("Power exponent of total protein (TP/68.95) on Vc/F (unitless)") # Table 2 'VC/F - TP = -3.1' (RSE 30.60%); Table 1 theta12
    e_sexf_vp <- 0.751; label("Fractional change in Vp/F for females (unitless)") # Table 2 'VP/F - SEX = 0.751' (RSE 24.80%); Table 1 theta13
    e_tumtp_nf1_fdepot <- 1.29; label("Relative bioavailability in NF1 vs melanoma (unitless ratio)") # Table 2 'F1* = 1.29' (RSE 6.50%); footnote: F1 is 1 for melanoma

    # IIV is reported as %CV of a log-normal eta (Methods 2.4.1); variances
    # are omega^2 = log(1 + CV^2). No covariance was reported, so the etas
    # are diagonal.
    etalcl ~ 0.03546 # Table 2 'IIV CL/F = 19%': log(1 + 0.19^2)
    etalvc ~ 0.39315 # Table 2 'IIV VC/F = 69.4%': log(1 + 0.694^2)
    etalq ~ 0.06541 # Table 2 'IIV Q/F = 26%': log(1 + 0.26^2)
    etalvp ~ 0.06787 # Table 2 'IIV VP/F = 26.5%': log(1 + 0.265^2)

    # Results 3.2: 'The residual error was modeled as a proportional residual
    # variability'; Table 2 lists only the proportional term.
    propSd <- 0.258; label("Proportional residual error (fraction)") # Table 2 'Proportional residual = 25.8%' (RSE 5.5%)
  })

  model({
    # Individual parameters (Tan 2023 Table 1).
    ka <- exp(lka)
    tlag <- exp(ltlag)
    cl <- exp(lcl + etalcl) * (BSA / 1.66)^e_bsa_cl
    vc <- exp(lvc + etalvc) * (BSA / 1.66)^e_bsa_vc * (TPRO / 68.95)^e_tpro_vc
    q <- exp(lq + etalq) * (BSA / 1.66)^e_bsa_q
    vp <- exp(lvp + etalvp) * (BSA / 1.66)^e_bsa_vp * (1 + e_sexf_vp * SEXF)

    # F1 = 1 for melanoma and theta7 for NF1.
    fdepot <- 1 * (1 - TUMTP_NF1) + e_tumtp_nf1_fdepot * TUMTP_NF1

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # mg/L -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
