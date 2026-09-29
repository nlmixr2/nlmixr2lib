Xie_2019_tofacitinib <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and an",
    "absorption lag time for oral tofacitinib in 650 adults with active",
    "psoriatic arthritis, pooled from the phase 3 OPAL Broaden (NCT01877668)",
    "and OPAL Beyond (NCT01882439) studies (Xie 2019). The model is",
    "parameterized in apparent oral clearance (CL/F) and apparent volume of",
    "distribution (V/F), so bioavailability is folded into both and F is not",
    "separately identifiable. A full (not reduced) covariate model is carried:",
    "CL/F varies with baseline age (power -0.20 on AGE/50), baseline",
    "creatinine clearance (power 0.32 on CRCL_BASE/120), baseline C-reactive",
    "protein (power -0.02 on CRP/0.49), and multiplicative factors for Black",
    "(0.91), Asian (0.95) and Other (0.96) race relative to White, for",
    "non-Hispanic ethnicity (1.12) relative to the Hispanic reference, and for",
    "female sex (1.05) relative to male; the body-weight effect on CL/F was",
    "fixed to zero by the authors because unconstrained runs returned a",
    "biologically implausible negative exponent. V/F varies with baseline age",
    "(power -0.22) and baseline body weight (power 0.68 on WT/83.3).",
    "Inter-individual variability is an exponential eta on CL/F and a large",
    "exponential eta on Ka; V/F carries no eta of its own -- its individual",
    "deviation is constructed as vc_eta_scale * etalcl (the paper's 'scaling",
    "parameter', 0.5), which forces a correlation of exactly 1 between the",
    "CL/F and V/F random effects. Residual error is proportional with a",
    "magnitude that switches on time after dose at 5 hours (22.9% CV at or",
    "before 5 h, 52.6% CV after), and the paper additionally estimated a",
    "65.8% CV inter-individual variability on the residual magnitude itself,",
    "carried here as etaruv. Only baseline creatinine clearance produced a",
    "clinically relevant exposure change; see the validation vignette for the",
    "reproduction of every covariate effect the paper reports.",
    sep = " "
  )
  reference <- paste(
    "Xie R, Deng C, Wang Q, Kanik KS, Nicholas T, Menon S.",
    "Population pharmacokinetics of tofacitinib in patients with psoriatic",
    "arthritis. Int J Clin Pharmacol Ther. 2019; 57(9): 464-473.",
    "doi:10.5414/CP203516. PMID: 31319909. PMCID: PMC6704728.",
    sep = " "
  )
  vignette <- "Xie_2019_tofacitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline (time-fixed). Enters CL/F and V/F as two separately estimated",
        "power functions, both normalized to the reference value of 50 years",
        "(the reference-patient definition in the Table 3 footnote and in",
        "Materials and methods). Cohort mean 49.0 years (SD 12.1), range 18-78",
        "(Table 2). The paper's covariate-impact assessment contrasts an",
        "80-year-old against the 50-year-old reference."
      ),
      source_name = "Age"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline (time-fixed), reported as BWT in the paper. Enters V/F as a",
        "power function normalized to 83.3 kg (the reference patient). It also",
        "appears in the CL/F covariate model, but with its exponent FIXED TO",
        "ZERO, so it has no effect on CL/F: initial runs estimated a",
        "biologically implausible negative exponent, and the authors fixed the",
        "effect to zero after confirming that the fit was essentially unchanged",
        "(Results). The zero-exponent term is retained in model() so that the",
        "authors' full covariate structure is reproduced literally rather than",
        "silently pruned. Cohort mean 84.8 kg (SD 19.0), range 38.1-159.7",
        "(Table 2); 61 and 109 kg are the 10th and 90th percentiles used in the",
        "Figure 3 covariate-impact assessment."
      ),
      source_name = "BWT"
    ),
    CRCL_BASE = list(
      description = "Baseline creatinine clearance (Cockcroft-Gault)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline, time-fixed, and NOT body-surface-area normalized -- the",
        "paper's BCCL is a raw Cockcroft-Gault creatinine clearance in mL/min",
        "(Materials and methods), which is why the canonical column is",
        "CRCL_BASE rather than the BSA-normalized CRCL. Enters CL/F as a power",
        "function normalized to 120 mL/min, described in the Results as the",
        "median value in the analysis dataset. Cohort mean 123.4 mL/min",
        "(SD 37.7), range 49.1-348.5 (Table 2); the paper notes there were no",
        "patients below 49 mL/min. This is the only covariate the authors judged",
        "to produce a clinically relevant change in exposure."
      ),
      source_name = "BCCL"
    ),
    CRP = list(
      description = "Baseline C-reactive protein",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline, time-fixed. UNITS TRAP: this paper reports CRP in mg/dL, not",
        "the mg/L used by most other entries in the covariate register -- the",
        "reference value 0.49 mg/dL is 4.9 mg/L. Enters CL/F as a power function",
        "normalized to 0.49 mg/dL. Cohort mean 1.1 mg/dL (SD 2.0), range",
        "0.0-16.4 (Table 2). The 95% CI of the exponent (-0.04, 0.00) touches",
        "the null value; it is retained because the paper reports a full",
        "covariate model estimated without stepwise selection."
      ),
      source_name = "BCRP"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = paste(
        "Table 3 reports the 'Female' effect as a MULTIPLICATIVE factor on CL/F",
        "(Equations 2 and 3: theta_i = theta_TV * theta_x for a non-reference",
        "category, theta_i = theta_TV for the reference category), so the null",
        "value is 1 and not 0. Male is the reference category. The 95% CI",
        "(0.99, 1.11) contains the null value; it is retained because the paper",
        "reports a full covariate model. Cohort 55.4% female (Table 2)."
      ),
      source_name = "Female"
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "White (all three race indicators = 0)",
      notes = paste(
        "One of three mutually exclusive non-reference race indicators (Black,",
        "Asian, Other) whose common reference category is White. Table 3 reports",
        "a MULTIPLICATIVE factor on CL/F of 0.91, so the null value is 1. Only 3",
        "of 650 patients (0.5%) were Black, which is why the 95% CI",
        "(0.42, 1.78) is very wide; the Discussion argues from a larger RA",
        "analysis that no clinically relevant race difference is expected."
      ),
      source_name = "Black"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "White (all three race indicators = 0)",
      notes = paste(
        "One of three mutually exclusive non-reference race indicators whose",
        "common reference category is White. Table 3 reports a MULTIPLICATIVE",
        "factor on CL/F of 0.95, so the null value is 1. Cohort 3.1% Asian",
        "(20 of 650; Table 2)."
      ),
      source_name = "Asian"
    ),
    RACE_OTHER = list(
      description = "Race-category 'Other' indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "White (all three race indicators = 0)",
      notes = paste(
        "One of three mutually exclusive non-reference race indicators whose",
        "common reference category is White. Table 3 reports a MULTIPLICATIVE",
        "factor on CL/F of 0.96, so the null value is 1. Cohort 2.6% 'Other'",
        "(17 of 650; Table 2)."
      ),
      source_name = "Other race"
    ),
    RACE_HISPANIC = list(
      description = "Hispanic / Latino ethnicity indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "Hispanic (RACE_HISPANIC = 1)",
      notes = paste(
        "POLARITY TRAP: this paper's reference category is the HISPANIC group,",
        "not the non-Hispanic group -- 'For ethnicity, Hispanic patients were",
        "used as the baseline typical patient for simulations' (Materials and",
        "methods), and the reference patient in the Table 3 footnote is",
        "Hispanic. Table 3 therefore reports a single MULTIPLICATIVE factor of",
        "1.12 that applies to NON-Hispanic patients. The canonical column",
        "direction is retained (1 = Hispanic), and the effect is applied through",
        "(1 - RACE_HISPANIC) in model() so that the register's polarity and the",
        "paper's reference category can both be honoured; the same device is",
        "used elsewhere in this package when a source codes the opposite level",
        "of a canonical binary. Only 10.6% of the cohort (69 of 650) was",
        "Hispanic, so the reference group is the small one here."
      ),
      source_name = "Non-Hispanic"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "tofacitinib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "tofacitinib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 650,
    n_studies = 2,
    n_observations = 3252,
    age_mean = "49.0 years (SD 12.1)",
    age_range = "18-78 years",
    weight_mean = "84.8 kg (SD 19.0)",
    weight_range = "38.1-159.7 kg",
    sex_female_pct = 55.4,
    race_ethnicity = c(White = 93.9, Asian = 3.1, Other = 2.6, Black = 0.5),
    ethnicity = c(Hispanic = 10.6, NonHispanic = 89.4),
    disease_state = "active psoriatic arthritis with an inadequate response to a conventional synthetic DMARD (OPAL Broaden, TNF-inhibitor naive) or to a TNF inhibitor (OPAL Beyond); all patients received a background csDMARD",
    dose_range = "5 or 10 mg orally twice daily",
    renal_function = "baseline Cockcroft-Gault creatinine clearance mean 123.4 mL/min (SD 37.7), range 49.1-348.5; no patient below 49 mL/min",
    disease_severity = "baseline C-reactive protein mean 1.1 mg/dL (SD 2.0), range 0.0-16.4",
    regions = "two global phase 3 studies: OPAL Broaden (NCT01877668, 12 months) and OPAL Beyond (NCT01882439, 6 months)",
    notes = paste(
      "Baseline demographics from Table 2 of Xie 2019. The 650 patients",
      "comprise 310 from OPAL Broaden (5 mg BID n = 104, 10 mg BID n = 104,",
      "placebo to 5 mg BID at month 3 n = 52, placebo to 10 mg BID at month 3",
      "n = 50) and 377 from OPAL Beyond (5 mg BID n = 127, 10 mg BID n = 126,",
      "placebo to 5 mg BID n = 64, placebo to 10 mg BID n = 60); the two study",
      "totals sum to 687 treatment-group assignments against 650 patients with",
      "sufficient data for analysis, and the paper does not reconcile the",
      "difference. Sampling was sparse and trough-weighted: month 1 (pre-dose",
      "and 2 h post-dose), month 4 and month 6 (pre-dose and 0.5, 2 and 3 h",
      "post-dose), with pre-dose samples drawn 12 +/- 2 h after the previous",
      "evening dose. Concentrations were measured by a validated LC-MS/MS",
      "assay; the observed range in Figure 1 spans roughly 0.1-200 ng/mL."
    )
  )

  ini({
    # ---- Structural parameters -------------------------------------------
    # All typical values are for the paper's reference patient: White, male,
    # Hispanic, body weight 83.3 kg, age 50 years, BCRP 0.49 mg/dL,
    # BCCL 120 mL/min (Table 3 footnote b; Materials and methods).
    lka <- log(13.8); label("Absorption rate constant (1/h)")                      # Table 3, 'Ka, /h' = 13.8 (RSE 7.9%), bootstrap 95% CI 12.1-16.6
    lcl <- log(20.4); label("Apparent oral clearance CL/F (L/h)")                  # Table 3, 'CL/F, L/h' = 20.4 (RSE 4.7%), bootstrap 95% CI 18.6-21.8
    lvc <- log(110); label("Apparent volume of distribution V/F (L)")              # Table 3, 'V/F, L' = 110 (RSE 1.2%), bootstrap 95% CI 108-113
    ltlag <- log(0.3); label("Absorption lag time (h)")                            # Table 3, 'Lag time, h' = 0.3 (RSE 0.8%), bootstrap 95% CI 0.3-0.3

    # ---- Covariate effects on CL/F ---------------------------------------
    # Continuous covariates are power functions normalized to the reference
    # patient's value (Equation 1). Categorical covariates are MULTIPLICATIVE
    # factors whose null value is 1, not 0 (Equations 2 and 3).
    e_age_cl <- -0.20; label("Power exponent on (AGE/50) for CL/F (unitless)")                                   # Table 3, covariate CL/F ~ Age = -0.20 (95% CI -0.31, -0.08)
    e_wt_cl <- fixed(0); label("Power exponent on (WT/83.3) for CL/F (unitless; zero by the authors' choice)")   # Table 3, covariate CL/F ~ BWT = 0 (FIX); footnote b
    e_crcl_base_cl <- 0.32; label("Power exponent on (CRCL_BASE/120) for CL/F (unitless)")                       # Table 3, covariate CL/F ~ BCCL = 0.32 (95% CI 0.21, 0.42)
    e_crp_cl <- -0.02; label("Power exponent on (CRP/0.49) for CL/F (unitless)")                                 # Table 3, covariate CL/F ~ BCRP = -0.02 (95% CI -0.04, 0.00)
    e_race_black_cl <- 0.91; label("Multiplicative factor on CL/F for Black vs White race (unitless)")           # Table 3, covariate CL/F ~ Black = 0.91 (95% CI 0.42, 1.78)
    e_race_asian_cl <- 0.95; label("Multiplicative factor on CL/F for Asian vs White race (unitless)")           # Table 3, covariate CL/F ~ Asian = 0.95 (95% CI 0.83, 1.08)
    e_race_other_cl <- 0.96; label("Multiplicative factor on CL/F for Other vs White race (unitless)")           # Table 3, covariate CL/F ~ Other race = 0.96 (95% CI 0.86, 1.12)
    e_race_hispanic_cl <- 1.12; label("Multiplicative factor on CL/F for non-Hispanic vs Hispanic (unitless)")   # Table 3, covariate CL/F ~ Non-Hispanic = 1.12 (95% CI 1.04, 1.24)
    e_sexf_cl <- 1.05; label("Multiplicative factor on CL/F for female vs male (unitless)")                      # Table 3, covariate CL/F ~ Female = 1.05 (95% CI 0.99, 1.11)

    # ---- Covariate effects on V/F ----------------------------------------
    e_age_vc <- -0.22; label("Power exponent on (AGE/50) for V/F (unitless)")                                    # Table 3, covariate V/F ~ Age = -0.22 (95% CI -0.32, -0.13)
    e_wt_vc <- 0.68; label("Power exponent on (WT/83.3) for V/F (unitless)")                                     # Table 3, covariate V/F ~ BWT = 0.68 (95% CI 0.55, 0.79)

    # ---- IIV -------------------------------------------------------------
    # The paper estimated an exponential eta on CL/F and on Ka, and reports NO
    # separate eta for V/F (Table 3 shows '-'). Instead 'a scaling parameter
    # was used to describe the impact of the IIV of the V/F on the IIV of
    # CL/F' (Results), i.e. the V/F deviation is a fixed multiple of the CL/F
    # deviation. That is what produces the 'covariance term' the Methods
    # mention, and it forces the CL/F-V/F random-effect correlation to exactly
    # 1. It is implemented in model() as vc_eta_scale * etalcl.
    #
    # SCALE CONVENTION (assumption; see the vignette Assumptions section).
    # Table 3's 'IIV (RSE %)' column is a percent, and the paper does not state
    # whether it is 100*sqrt(omega^2) or the exact log-normal CV
    # 100*sqrt(exp(omega^2)-1). The NONMEM/PsN reporting convention used
    # throughout this paper's sibling analyses is 100*sqrt(omega^2), so the
    # variances below are (IIV%/100)^2. This matters little for CL/F (31.7%)
    # but a great deal for Ka (198%): under the exact-CV reading omega^2 would
    # be log(1+1.98^2) = 1.59 rather than 3.92.
    etalcl ~ 0.317^2                                                                                            # Table 3, 'CL/F, L/h' IIV = 31.7% (RSE 3.7%), bootstrap 95% CI 28.2-34.8
    etalka ~ 1.98^2                                                                                             # Table 3, 'Ka, /h' IIV = 198% (RSE 4.2%), bootstrap 95% CI 180-222
    etaruv ~ 0.658^2                                                                                            # Table 3, IIV on both proportional-error rows = 65.8% (RSE 5.5%), bootstrap 95% CI 57.2-73.1
    vc_eta_scale <- 0.5; label("Scaling factor relating eta_V/F to eta_CL/F (correlation fixed to 1; eta_V/F = vc_eta_scale * etalcl)")  # Table 3, 'Scaling parameter' = 0.5 (RSE 9.0%), bootstrap 95% CI 0.4-0.5

    # ---- Residual error --------------------------------------------------
    # Two proportional-error magnitudes split on time after dose at 5 hours
    # (Results: 'two proportional-error models for concentrations collected
    # > 5 hours postdose and <= 5 hours postdose'). Both share the single
    # 65.8% CV inter-individual variability on the residual magnitude above.
    propSd_early <- 0.229; label("Proportional residual error SD, time after dose <= 5 h (fraction)")            # Table 3, 'Proportional error, TAD <= 5 hours, %' = 22.9 (RSE 3.8%)
    propSd_late <- 0.526; label("Proportional residual error SD, time after dose > 5 h (fraction)")              # Table 3, 'Proportional error, TAD > 5 hours, %' = 52.6 (RSE 4.7%)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Covariate models (Equations 1-3).
    #
    #    Equation 1 (continuous): theta_i = theta_TV * (cov_i/cov_median)^theta_x
    #    Equation 2 (categorical, non-reference level): theta_i = theta_TV * theta_x
    #    Equation 3 (categorical, reference level):     theta_i = theta_TV
    #
    #    The three race indicators are mutually exclusive and share the White
    #    reference level, so the linear (1 + sum of (factor - 1) * indicator)
    #    form below reproduces Equation 2 exactly for whichever indicator is 1
    #    and Equation 3 when all are 0.
    #
    #    The ethnicity term is inverted: the paper's reference level is
    #    Hispanic, so the 1.12 factor applies when RACE_HISPANIC is 0.
    #
    #    The body-weight term on CL/F has an exponent fixed to zero, so it
    #    evaluates to 1 for every subject; it is kept to reproduce the authors'
    #    full covariate structure literally.
    #
    #    Verified against every covariate impact the paper reports in prose:
    #      (80/50)^-0.20      = 0.9103 ->  8.97% lower CL/F at age 80   [paper: 8.97%]
    #      (80/50)^-0.22      = 0.9018 ->  9.82% lower V/F  at age 80   [paper: 9.9%]
    #      (50/120)^0.32      = 0.7557 -> 24.43% lower CL/F at BCCL 50  [paper: 24.3%]
    #      (60/120)^0.32      = 0.8011 -> 19.89% lower CL/F at BCCL 60  [paper: 19.8%]
    #      (3.0/0.49)^-0.02   = 0.9644 ->  3.56% lower CL/F at BCRP 3.0 [paper: 3.6%]
    #      1.05               ->  5.0% higher CL/F in females           [paper: 5.4%]
    #      1.12               -> 12.0% higher CL/F in non-Hispanics     [paper: 12.3%]
    #      (61/83.3)^0.68     = 0.8091 -> 19.1% lower  V/F at 61 kg     [paper: ~19%]
    #      (109/83.3)^0.68    = 1.2006 -> 20.1% higher V/F at 109 kg    [paper: ~20%]
    #      1/(50/120)^0.32    = 1.3233 -> 32.3% higher AUC at BCCL 50   [paper: 32%]
    #    The two categorical values are the only ones that do not reproduce to
    #    the paper's own precision; both are consistent with Table 3's two
    #    decimal places (1.054 rounds to 1.05 and 1.123 rounds to 1.12), i.e.
    #    the prose percentages were computed from unrounded NONMEM estimates.
    # ------------------------------------------------------------------
    cov_cl <- (AGE / 50)^e_age_cl *
      (WT / 83.3)^e_wt_cl *
      (CRCL_BASE / 120)^e_crcl_base_cl *
      (CRP / 0.49)^e_crp_cl *
      (1 + (e_sexf_cl - 1) * SEXF) *
      (1 + (e_race_hispanic_cl - 1) * (1 - RACE_HISPANIC)) *
      (1 + (e_race_black_cl - 1) * RACE_BLACK +
        (e_race_asian_cl - 1) * RACE_ASIAN +
        (e_race_other_cl - 1) * RACE_OTHER)
    cov_vc <- (AGE / 50)^e_age_vc * (WT / 83.3)^e_wt_vc

    # ------------------------------------------------------------------
    # 2. Individual parameters. V/F has no eta of its own: its individual
    #    log-scale deviation is vc_eta_scale * etalcl (perfect correlation
    #    with the CL/F deviation), which is the paper's 'scaling parameter'.
    # ------------------------------------------------------------------
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * cov_cl
    vc <- exp(lvc + vc_eta_scale * etalcl) * cov_vc
    tlag <- exp(ltlag)

    # ------------------------------------------------------------------
    # 3. Micro-constants.
    # ------------------------------------------------------------------
    kel <- cl / vc

    # ------------------------------------------------------------------
    # 4. ODE system: one-compartment disposition with first-order
    #    absorption from a depot preceded by an absorption lag.
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # 5. Observation and error.
    #    Doses are in mg and vc is in L, so central/vc is mg/L; the factor
    #    1000 converts to ng/mL, the assay unit of Figures 1 and 2.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc

    #    The proportional residual magnitude switches on time after dose at
    #    5 h, and carries its own inter-individual variability (etaruv), so
    #    the magnitude is assembled as a model variable and handed to prop()
    #    (the Bukkems_2021_raltegravir TAD-switch precedent combined with the
    #    Ooi_2026_elafibranor eta-on-residual precedent). tad() is evaluated
    #    once and reused: calling a time function twice inside one expression
    #    fails to parse in rxode2.
    early_flag <- tad() <= 5
    propSdTad <- propSd_early * early_flag + propSd_late * (1 - early_flag)
    propSdCc <- propSdTad * exp(etaruv)
    Cc ~ prop(propSdCc)
  })
}
