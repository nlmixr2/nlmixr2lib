Mouksassi_2015_thrombomodulinAlfa <- function() {
  description <- paste(
    "One-compartment intravenous population PK model for thrombomodulin",
    "alfa (ART-123, a soluble recombinant human thrombomodulin) pooled",
    "across 24 healthy adults given single 0.02 or 0.06 mg/kg IV doses and",
    "368 adults with sepsis and suspected disseminated intravascular",
    "coagulation given 0.06 mg/kg/day (maximum 6 mg) IV for six consecutive",
    "days (Mouksassi 2015; N = 392). Typical clearance is 0.14 L/h and",
    "typical central volume 5.47 L at the reference 68 kg body weight and",
    "72.6 mL/min creatinine clearance, giving a typical elimination",
    "half-life of 27.08 h. Clearance carries a power effect of body weight",
    "(exponent 0.56) and of baseline Cockcroft-Gault creatinine clearance",
    "truncated at 150 mL/min (exponent 0.27); central volume carries a power",
    "effect of body weight (exponent 0.50). Between-subject variability is",
    "log-normal on clearance (variance 0.11, i.e. 33% on the SD scale) and",
    "on central volume (variance 0.15, i.e. 39%), correlated at 0.57.",
    "Residual variability is proportional (29%). The analysis supports the",
    "conclusion that no dose adjustment is needed in sepsis and DIC patients",
    "with mild to severe renal impairment, because simulated exposures stay",
    "inside the 300-5400 ng/mL therapeutic range.",
    sep = " "
  )
  reference <- paste(
    "Mouksassi MS, Marier JF, Bax L, Osawa Y, Tsuruta K. (2015). Population",
    "Pharmacokinetic Analysis of Thrombomodulin Alfa to Support Dosing",
    "Rationale in Patients with Renal Impairment.",
    "Clin Pharmacol Drug Dev. 4(3):210-217.",
    "doi:10.1002/cpdd.163",
    sep = " "
  )
  vignette <- "Mouksassi_2015_thrombomodulinAlfa"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Confirmed against Mouksassi 2015 Methods "Population PK
  # Modeling" (one-compartment model fitted to intravenously administered
  # ART-123) and "Bioanalytical Assay" (PLASMA concentrations of thrombomodulin
  # alfa by sandwich ELISA, lower determination limit 1 ng/mL).
  compartmentData <- list(
    central = list(analyte = "thrombomodulin alfa", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters both clearance and central volume as a power function",
        "normalized to 68 kg, which is the median body weight of the pooled",
        "392-subject analysis population (Mouksassi 2015 Table 1, Combined",
        "row). Exponents are estimated, not fixed to allometric 0.75 / 1:",
        "Table 2 reports 0.56 (14.15% SE) on CL and 0.50 (21.14% SE) on V.",
        "Treated as a baseline, time-fixed value; the paper screened body",
        "weight both as a continuous and as a categorical (< or >= 100 kg)",
        "covariate and retained the continuous power form.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance by the Cockcroft-Gault equation, raw (not BSA-normalized), truncated at 150 mL/min",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column CRCL. Computed with the Cockcroft-Gault equation in",
        "raw mL/min and NOT BSA-normalized to mL/min/1.73 m^2; stored under",
        "the canonical CRCL column per inst/references/covariate-columns.md,",
        "which admits raw mL/min when the source paper does not normalize",
        "(same convention as Delattre_2010_amikacin.R and",
        "Chen_2023_nemonoxacin.R). Two data-handling rules from Mouksassi",
        "2015 Methods and Results must be reproduced when using this model:",
        "(1) values above 150 mL/min are TRUNCATED to 150 ('The upper range",
        "of CrCL was truncated to 150 mL/min for the covariate analysis in",
        "order to better represent upper physiological range of CrCL'; this",
        "affected 1.6% of Phase I and 9.5% of Phase IIb subjects), and",
        "(2) the BASELINE value is used - the paper screened baseline CrCL",
        "and time-varying CrCL separately and the final model retained",
        "'creatinine clearance at baseline'. Enters clearance as a power",
        "function normalized to 72.6 mL/min. Patients on renal replacement",
        "therapy for end-stage renal disease or acute kidney injury were",
        "excluded from the Phase IIb study, so the model is not qualified",
        "for dialysis or ESRD.",
        sep = " "
      ),
      source_name        = "CRCL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 392L,
    n_studies      = 2L,
    age_range      = "18-93 years",
    age_median     = "59 years",
    weight_range   = "30-150 kg",
    weight_median  = "68 kg",
    sex_female_pct = 37.4,
    race_ethnicity = c(Caucasian = 62.1, Indian = 26.1, Asian = 5.8, Black = 3.4, Other = 2.4, Hispanic = 0.3),
    disease_state  = paste(
      "Pooled: 24 normal healthy volunteers (Phase I single-dose study) and",
      "368 adults with sepsis and suspected disseminated intravascular",
      "coagulation (Phase IIb randomized, double-blind, placebo-controlled",
      "study; the 368 are the ART-123-treated subset with plasma samples out",
      "of 750 randomized). Patients on renal replacement therapy for",
      "end-stage renal disease or acute kidney injury were excluded.",
      sep = " "
    ),
    dose_range     = paste(
      "Phase I: single IV 0.02 or 0.06 mg/kg. Phase IIb: 0.06 mg/kg/day up",
      "to a maximum of 6 mg once daily for six consecutive days, given as an",
      "IV bolus injection or a 15-minute IV infusion. 1 mg of thrombomodulin",
      "alfa is approximately 6400 U.",
      sep = " "
    ),
    regions        = "Not reported by region; the Phase IIb study population includes Caucasian, Indian, Asian, Black, Other and Hispanic subjects",
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance, truncated at 150 mL/min:",
      "combined mean 74.6 mL/min (range 5.0-150). Phase I mean 114.2",
      "(range 74-150); Phase IIb mean 72.0 (range 5-150). The renal-function",
      "strata used for the dosing-rationale simulations are normal",
      "(90-120 mL/min, median 105), mild (60-89, median 75), moderate",
      "(30-59, median 45) and severe (15-30, median 22).",
      sep = " "
    ),
    notes          = paste(
      "Baseline demographics per Mouksassi 2015 Table 1. Sex and race",
      "percentages are reported over 380 subjects (142 female / 238 male;",
      "236 Caucasian, 99 Indian, 22 Asian, 13 Black, 9 Other, 1 Hispanic),",
      "12 fewer than the 392 subjects in the PK analysis - the paper does",
      "not explain the difference. Sampling was rich in Phase I (baseline;",
      "5, 10, 15, 30 min; 1, 2, 4, 6, 8, 12, 24, 36, 48, 72, 96, 120, 144 h)",
      "and sparse in Phase IIb (12 +/- 3 h after the Day 1 dose, < 30 min",
      "before the Day 3 dose, 12 +/- 3 h after the Day 3 dose, Day 7, Day 14",
      "and in some cases Day 28). Endogenous thrombomodulin in the placebo",
      "arm (mean +/- SD 8.6 +/- 19.1 ng/mL) was measured but excluded from",
      "the analysis. Estimation was performed in Phoenix NLME 1.2, not",
      "NONMEM. The therapeutic range of the product is 300-5400 ng/mL: the",
      "lower bound is the effective concentration from non-clinical and",
      "clinical data, the upper bound is the highest concentration with no",
      "bleeding event in non-clinical toxicology.",
      sep = " "
    )
  )

  ini({
    # Structural parameters at the reference subject: WT 68 kg, CRCL 72.6
    # mL/min. Mouksassi 2015 Table 2, "Typical value (% SE)" column, whose CL
    # and V rows print the full covariate equations
    #   CL = 0.14 (2.01) * (WT/68)^0.56 (14.15) * (CRCL/72.6)^0.27 (12.56)
    #   V  = 5.47 (2.71) * (WT/68)^0.50 (21.14)
    # with the % SE of each element in parentheses immediately after it.
    lcl <- log(0.14); label("Clearance at WT 68 kg and CRCL 72.6 mL/min (L/h)")            # Mouksassi 2015 Table 2, CL row: 0.14 L/h (2.01% SE)
    lvc <- log(5.47); label("Central volume of distribution at WT 68 kg (L)")              # Mouksassi 2015 Table 2, V row: 5.47 L (2.71% SE)

    # Covariate effects, all estimated (each carries a % SE in Table 2 and each
    # was added by the stepwise procedure in Supplemental Table S2).
    e_wt_cl   <- 0.56; label("Power exponent on (WT / 68) for clearance (unitless)")       # Mouksassi 2015 Table 2, CL row: 0.56 (14.15% SE); Table S2 step 1
    e_crcl_cl <- 0.27; label("Power exponent on (CRCL / 72.6) for clearance (unitless)")   # Mouksassi 2015 Table 2, CL row: 0.27 (12.56% SE); Table S2 step 2
    e_wt_vc   <- 0.50; label("Power exponent on (WT / 68) for central volume (unitless)")  # Mouksassi 2015 Table 2, V row: 0.50 (21.14% SE); Table S2 step 3

    # Correlated log-normal between-subject variability on CL and Vc.
    #
    # SCALE. Mouksassi 2015 Table 2's BSV column is a VARIANCE - stated by the
    # table footnote ("BSV = between-subject variability (variance)") and
    # confirmed independently by Figure 3, which annotates the same two
    # distributions as "BSV = 33 %" for CL/F and "BSV = 39 %" for V/F:
    # sqrt(0.11) = 0.332 and sqrt(0.15) = 0.387. So no CV-to-variance
    # conversion is applied here; the tabulated numbers ARE the variances.
    #
    # COVARIANCE. The Table 2 footnote reads "Covariance was 0.7 (correlation
    # = 0.57)". A covariance of 0.7 is impossible for these variances - the
    # Cauchy-Schwarz bound is sqrt(0.11 * 0.15) = 0.1285 - so the printed 0.7
    # has lost a leading zero. The printed CORRELATION is admissible and is
    # used to reconstruct the covariance: 0.57 * sqrt(0.11 * 0.15) = 0.0732,
    # which round-trips to the printed 0.57 exactly and to "0.07" under the
    # decimal-slip reading. That a covariance was estimated at all is stated
    # twice more: Methods ("Covariance between ETA terms was estimated in the
    # model where correlations between ETAs were deemed probable") and Results
    # ("A 1-compartment model with shared omega for clearance (CL) and volume
    # of distribution (V) resulted in the lowest objective function").
    etalcl + etalvc ~ c(0.11,
                        0.0732, 0.15)                                                     # Mouksassi 2015 Table 2 BSV column (0.11, 10.59% SE; 0.15, 14.21% SE) and footnote a

    # Residual error. Supplemental Table S1 selected the proportional error
    # model (AIC 22,618) over mixed (22,624) and additive (23,645).
    propSd <- 0.29; label("Proportional residual error (fraction)")                        # Mouksassi 2015 Table 2, "Proportional error" row: 0.29 (2.43% SE)
  })
  model({
    # Individual PK parameters. Both covariate effects are power functions
    # normalized to the reference covariate values printed inside the Table 2
    # equations (68 kg, 72.6 mL/min). CRCL must be supplied as the BASELINE
    # Cockcroft-Gault value truncated at 150 mL/min - see covariateData.
    cl <- exp(lcl + etalcl) * (WT / 68)^e_wt_cl * (CRCL / 72.6)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 68)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Doses are in mg and vc is in L, so central / vc is mg/L. The factor 1000
    # converts to ng/mL, the unit the assay reports (lower determination limit
    # 1 ng/mL) and the unit of the 300-5400 ng/mL therapeutic range and of the
    # Figure 4 y-axis.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
