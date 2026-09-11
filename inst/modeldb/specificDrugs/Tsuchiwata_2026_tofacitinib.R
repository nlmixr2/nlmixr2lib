Tsuchiwata_2026_tofacitinib <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral tofacitinib in adults with active ankylosing spondylitis"
  reference <- paste(
    "Tsuchiwata S, Suzuki A, Wang Q, Kanik K, Fallon L, Menon S.",
    "Population pharmacokinetics of tofacitinib in patients with active",
    "ankylosing spondylitis. Int J Clin Pharmacol Ther. 2026; 64(2): 57-65.",
    "doi:10.5414/CP204781. PMID: 41355395.",
    sep = " "
  )
  vignette <- "Tsuchiwata_2026_tofacitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline (time-fixed). Enters CL/F and V/F as separate power functions",
        "normalized to the reference value of 40 years (the approximate dataset",
        "median; Table 2 footnote b and the reference-patient definition in the",
        "Results). Cohort mean 41.8 years (SD 11.7); 64 years is the 95th",
        "percentile used for the Figure 2 covariate-impact assessment."
      ),
      source_name        = "Age"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male (SEXF = 0)",
      notes              = paste(
        "Table 2 reports the effect as 'Sex: Female (vs. male)', so male is the",
        "reference category and the coefficient is the fractional change in CL/F",
        "for a female patient: CL/F is multiplied by (1 + e_sexf_cl * SEXF).",
        "The 95% CI of this effect contains the null value; it is retained",
        "because the paper reports a full (not reduced) covariate model."
      ),
      source_name        = "Sex"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "non-Asian (RACE_ASIAN = 0)",
      notes              = paste(
        "Table 2 reports the effect as 'Race: Asian (vs. non-Asian)', a fractional",
        "change applied as CL/F * (1 + e_race_asian_cl * RACE_ASIAN). The cohort is",
        "79.9% White, 19.7% Asian and 0.4% not reported (Table 1), so the non-Asian",
        "reference group is predominantly White."
      ),
      source_name        = "Race"
    ),
    CRCL_BASE = list(
      description        = "Baseline creatinine clearance (Cockcroft-Gault)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline, time-fixed, and NOT body-surface-area normalized -- the paper's",
        "BCCL is a raw Cockcroft-Gault creatinine clearance in mL/min, which is why",
        "the canonical column is CRCL_BASE rather than the BSA-normalized CRCL.",
        "Enters CL/F as a power function normalized to 126 mL/min. Cohort mean",
        "129 mL/min (SD 33.7); the lowest value in the analysis dataset was",
        "48.1 mL/min (Figure 2 footnote a)."
      ),
      source_name        = "BCCL"
    ),
    CRP = list(
      description        = "Baseline C-reactive protein",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline, time-fixed, standard (not high-sensitivity) assay. UNITS TRAP:",
        "this paper reports CRP in mg/dL, not the mg/L used by most other entries",
        "in the covariate register -- the reference value 0.851 mg/dL is 8.51 mg/L.",
        "Enters CL/F as a power function normalized to 0.851 mg/dL. Cohort mean",
        "1.41 mg/dL (SD 1.55; Table 1). The 95% CI of this effect contains the",
        "null value; it is retained because the paper reports a full covariate model."
      ),
      source_name        = "BCRP"
    ),
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline, time-fixed. Enters V/F only, as a power function normalized to",
        "78 kg (the dataset median). Body weight was pre-specified as a candidate",
        "predictor of CL/F as well, but was dropped from the final full model",
        "because it correlated with BCCL (r = 0.58) and its inclusion changed the",
        "objective function by only -0.066 (Results, Final full model)."
      ),
      source_name        = "Body weight"
    )
  )

  covariatesDataExcluded <- list(
    RACE_HISPANIC = list(
      description = "Hispanic / Latino ethnicity indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Pre-specified as a candidate predictor of CL/F but not carried into the",
        "final full model because 97.8% of the dataset was non-Hispanic/Latino",
        "(Results, Final full model development). No coefficient is reported."
      )
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "tofacitinib", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "tofacitinib", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 279,
    n_studies      = 2,
    n_observations = 1917,
    age_mean       = "41.8 years (SD 11.7)",
    weight_mean    = "78.1 kg (SD 17.4)",
    weight_median  = "78 kg",
    weight_range   = "54-107 kg (5th-95th percentiles)",
    sex_female_pct = 20.4,
    race_ethnicity = c(White = 79.9, Asian = 19.7, NotAvailable = 0.4),
    disease_state  = "active ankylosing spondylitis (modified New York criteria, BASDAI >= 4 and BASDAI back-pain score >= 4)",
    dose_range     = "2, 5 or 10 mg orally twice daily",
    renal_function = "baseline Cockcroft-Gault creatinine clearance mean 129 mL/min (SD 33.7); lowest observed 48.1 mL/min",
    regions        = "not reported by region; pooled phase 2 (NCT01786668) and phase 3 (NCT03502616) studies",
    notes          = paste(
      "Baseline demographics from Table 1 of Tsuchiwata 2026. Phase 2 (N = 147;",
      "2 mg BID n = 50, 5 mg BID n = 49, 10 mg BID n = 48) and phase 3 (N = 132,",
      "all 5 mg BID). Sampling was sparse: week 4 (pre-dose, 0.5 h and 2 h",
      "post-dose) and week 8 (pre-dose, 0.5, 2 and 3 h post-dose). Assay range",
      "0.1-100 ng/mL; ~1% of observations were below the limit of quantification",
      "and were treated as missing."
    )
  )

  ini({
    # Structural parameters. All are typical values for the paper's reference
    # patient: White, male, 78 kg, 40 years old, BCCL 126 mL/min,
    # BCRP 0.851 mg/dL (Table 2 footnote b; Results, Final full model).
    lka <- log(3.07); label("Absorption rate constant (1/h)")                 # Table 2, 'K a /hour -1' = 3.07 (RSE 10.2%)
    lcl <- log(27.1); label("Apparent oral clearance (L/h)")                  # Table 2, 'CL/F, L/hour' = 27.1 (RSE 2.41%)
    lvc <- log(126); label("Apparent volume of distribution (L)")             # Table 2, 'V/F, L' = 126 (RSE 2.81%)

    # Covariate effects on CL/F. Continuous covariates are power functions
    # normalized to the reference (approximate median) value; categorical
    # covariates are entered as one coefficient giving the fractional change
    # relative to the reference category (Table 2 footnote b).
    e_age_cl <- -0.244; label("Power exponent on (AGE/40) for CL/F (unitless)")                                # Table 2, covariate CL/F ~ Age
    e_sexf_cl <- 0.0237; label("Fractional change in CL/F for female vs male (unitless)")                      # Table 2, covariate CL/F ~ Sex: Female (vs. male)
    e_race_asian_cl <- -0.103; label("Fractional change in CL/F for Asian vs non-Asian (unitless)")            # Table 2, covariate CL/F ~ Race: Asian (vs. non-Asian)
    e_crcl_base_cl <- 0.233; label("Power exponent on (CRCL_BASE/126) for CL/F (unitless)")                    # Table 2, covariate CL/F ~ BCCL
    e_crp_cl <- -0.0185; label("Power exponent on (CRP/0.851) for CL/F (unitless)")                            # Table 2, covariate CL/F ~ BCRP

    # Covariate effects on V/F.
    e_age_vc <- -0.230; label("Power exponent on (AGE/40) for V/F (unitless)")                                 # Table 2, covariate V/F ~ Age
    e_wt_vc <- 0.574; label("Power exponent on (WT/78) for V/F (unitless)")                                    # Table 2, covariate V/F ~ Body weight

    # IIV. The paper fitted an exponential IIV model on CL/F and V/F with an
    # OMEGA BLOCK (Methods, Base structural model). Table 2's 'IIV' column is a
    # percent CV, so the variance is (CV/100)^2 = 0.282^2 and 0.366^2; the
    # off-diagonal is the covariance reported directly on the variance scale.
    # Implied correlation 0.0760 / sqrt(0.282^2 * 0.366^2) = 0.736, and the
    # 2x2 block is positive definite (determinant 0.00488).
    etalcl + etalvc ~ c(0.282^2, 0.0760, 0.366^2)                                                             # Table 2, 'IIV (RSE%)' 28.2 and 36.6 and row 'Covariance, CL/F-V/F' = 0.0760

    # Residual error. Two proportional error models were used, split on time
    # after dose at 9 hours (Methods, Base structural model and random-effects
    # model development).
    propSd_nontrough <- 0.602; label("Proportional residual error SD, time after dose < 9 h (fraction)")       # Table 2, 'Proportional error CV, TAD < 9 hours,%' = 60.2
    propSd_trough <- 0.696; label("Proportional residual error SD, time after dose >= 9 h (fraction)")         # Table 2, 'Proportional error CV, TAD >= 9 hours,%' = 69.6
  })

  model({
    # 1. Covariate models.
    #    Continuous covariates enter as power functions normalized to the
    #    reference (approximate median) value; each categorical covariate enters
    #    as a single coefficient giving the fractional change relative to its
    #    reference category (Table 2 footnote b). Verified against the paper's
    #    own reported covariate impacts (Results, Impact of covariates):
    #      (64/40)^-0.244    = 0.891 -> 10.9% lower CL/F at age 64
    #      1 + 0.0237        = 1.024 ->  2.4% higher CL/F in females
    #      1 - 0.103         = 0.897 -> 10.3% lower CL/F in Asian patients
    #      (50/126)^0.233    = 0.806 -> 19.4% lower CL/F at BCCL 50 mL/min
    #      (8/0.851)^-0.0185 = 0.959 ->  4.1% lower CL/F at BCRP 8 mg/dL
    #      (64/40)^-0.230    = 0.898 -> 10.2% lower V/F at age 64
    #      (54/78)^0.574     = 0.810 -> 19%   lower V/F at 54 kg
    #      (107/78)^0.574    = 1.199 -> 20%   higher V/F at 107 kg
    cov_cl <- (AGE / 40)^e_age_cl *
      (CRCL_BASE / 126)^e_crcl_base_cl *
      (CRP / 0.851)^e_crp_cl *
      (1 + e_sexf_cl * SEXF) *
      (1 + e_race_asian_cl * RACE_ASIAN)
    cov_vc <- (AGE / 40)^e_age_vc * (WT / 78)^e_wt_vc

    # 2. Individual parameters. No IIV was reported on ka (Table 2, 'NA').
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * cov_cl
    vc <- exp(lvc + etalvc) * cov_vc

    # 3. Micro-constants.
    kel <- cl / vc

    # 4. ODE system: one-compartment disposition with first-order absorption.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 5. Observation and error.
    #    Doses are in mg and vc is in L, so central/vc is mg/L; the factor 1000
    #    converts to ng/mL, the assay unit used throughout the paper.
    Cc <- 1000 * central / vc

    #    The proportional residual error switches on time after dose at 9 h, so
    #    the error magnitude is assembled as a model variable and passed to
    #    prop() (same construction as Bukkems_2021_raltegravir). tad() is
    #    evaluated once and reused: calling a time function twice inside one
    #    expression fails to parse in rxode2.
    trough_flag <- tad() >= 9
    propSdTad <- propSd_trough * trough_flag + propSd_nontrough * (1 - trough_flag)
    Cc ~ prop(propSdTad)
  })
}
