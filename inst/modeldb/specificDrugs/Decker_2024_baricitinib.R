Decker_2024_baricitinib <- function() {
  description <- "Two-compartment population PK model with zero-order absorption, an absorption lag, and fixed allometric scaling for oral baricitinib in 217 pediatric patients aged 2 to <18 years with polyarticular-course juvenile idiopathic arthritis (JUVE-BASIS, NCT03773978). Apparent total clearance is partitioned semi-mechanistically into an eGFR-dependent apparent renal arm (CLr/F) and an eGFR-independent apparent non-renal arm (CLnr/F). Baseline body weight enters through allometric exponents fixed at 0.75 on all clearance terms (CLr/F, CLnr/F, Q) and at 1 on both volumes (V1/F, V2/F), referenced to 74 kg. The structure was carried unchanged from the adult rheumatoid-arthritis baricitinib population PK model; no further covariate met the stepwise-covariate-modeling inclusion criteria, so the base model is the final model."
  reference <- paste(
    "Decker RL, Steven Ernest C II, Radtke DB, Wang R, Araujo J, Keller SY, Zhang X.",
    "A population pharmacokinetic model using allometric scaling for baricitinib in patients",
    "with juvenile idiopathic arthritis.",
    "CPT Pharmacometrics Syst Pharmacol 2024;13(6):970-981.",
    "doi:10.1002/psp4.13131.",
    sep = " "
  )
  vignette <- "Decker_2024_baricitinib"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight (weight at study entry, WTE in the source paper).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at the study-entry value; Figure S5 legend defines WT as baseline body weight. Allometric reference weight is 74 kg (Table 2 footnotes b-e), which is the adult rheumatoid-arthritis reference weight carried over from the adult model, NOT a pediatric cohort statistic. Cohort mean 50.2 kg, range 11.0-111 kg (Table 1); 87% of patients weighed >=30 kg. Simulations in the paper spanned 10-120 kg.",
      source_name        = "WTE"
    ),
    CRCL_BASE = list(
      description        = "Baseline estimated glomerular filtration rate (bedside Schwartz equation), BSA-normalized and time-fixed per subject.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters CLr/F as the ratio (CRCL_BASE / 93); 93 mL/min/1.73 m^2 is the median eGFR of the previous adult population PK analysis and is the reference at which the reported CLr/F typical value of 6.34 L/h applies (Table 2 footnote f). Estimated with the bedside Schwartz equation (Table 1 abbreviations), a creatinine-based pediatric eGFR already expressed per 1.73 m^2. Cohort mean 119 mL/min/1.73 m^2, range 66.5-201 (Table 1).",
      source_name        = "baseline eGFR"
    ),
    CRCL = list(
      description        = "Time-varying estimated glomerular filtration rate (bedside Schwartz equation), BSA-normalized.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used only through the within-subject deviation (CRCL - CRCL_BASE), which is the 'change in eGFR from baseline' (deltaeGFR) term of Table 2 footnote f -- the Wahlby 2004 baseline/difference (BCOV/DCOV) decomposition of a time-varying covariate. Set CRCL = CRCL_BASE to switch the time-varying arm off; the model then reduces to baseline-eGFR-only scaling. The source paper does not tabulate the distribution of on-treatment eGFR changes.",
      source_name        = "eGFR"
    )
  )

  # Screened in the stepwise covariate model but not retained: no covariate met
  # the forward-inclusion criteria (p <= 0.01, i.e. dOFV >= 6.635, AND a >=5%
  # reduction in the BSV variance of the affected parameter), so backward
  # exclusion was never run and the base model became the final model (Methods,
  # "Pharmacokinetic model development"). Documented here to preserve the
  # provenance of the covariate screen without declaring unused covariates.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at study entry.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a continuous covariate and not retained. The Discussion states explicitly that after the effect of weight was accounted for, age was not a significant covariate on the clearance or volume terms; age and weight were highly correlated (Figure S5). Cohort mean 13.3 years, range 2-17 (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units       = "(binary)",
      type        = "binary",
      reference_category = "male",
      notes       = "Screened as 'gender' in the stepwise covariate model and not retained. Cohort: 150 of 217 female (69%) (Table 1)."
    ),
    RACE_WHITE = list(
      description = "White race indicator.",
      units       = "(binary)",
      type        = "binary",
      reference_category = "non-White",
      notes       = "Screened as part of the categorical 'race' covariate and not retained. 149 of 217 White (69%) (Table 1)."
    ),
    RACE_BLACK = list(
      description = "Black race indicator.",
      units       = "(binary)",
      type        = "binary",
      reference_category = "non-Black",
      notes       = "Screened as part of the categorical 'race' covariate and not retained. 5 of 217 Black (Table 1)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator.",
      units       = "(binary)",
      type        = "binary",
      reference_category = "non-Asian",
      notes       = "Screened as part of the categorical 'race' covariate and not retained. 48 of 217 Asian (Table 1)."
    ),
    RACE_JAPANESE = list(
      description = "Japanese-heritage indicator.",
      units       = "(binary)",
      type        = "binary",
      reference_category = "non-Japanese",
      notes       = "Screened as the 'Japanese versus non-Japanese' categorical covariate and not retained. The source paper does not tabulate the Japanese subgroup size."
    )
  )

  # 'Ethnicity' and 'JIA sub-type' (polyarticular RF+/RF-, extended oligoarticular,
  # enthesitis-related arthritis, juvenile psoriatic arthritis) were also screened
  # and not retained; neither has a canonical covariate-column name in
  # inst/references/covariate-columns.md and neither is used by this model, so they
  # are recorded here in prose rather than as covariatesDataExcluded entries.

  compartmentData <- list(
    central     = list(analyte = "baricitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "baricitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 217L,
    n_studies      = 1L,
    n_observations = 1261L,
    age_range      = "2-17 years (enrolment criterion 2 to <18 years)",
    age_mean       = "13.3 years",
    weight_range   = "11.0-111 kg",
    weight_mean    = "50.2 kg",
    sex_female_pct = 69.1,
    race_ethnicity = c(White = 68.7, Asian = 22.1, `Native American` = 3.2, Black = 2.3, Multiple = 0.9, Missing = 2.8),
    disease_state  = "Polyarticular-course juvenile idiopathic arthritis (polyarticular RF-positive or RF-negative, extended oligoarticular, enthesitis-related arthritis, or juvenile psoriatic arthritis) with inadequate response or intolerance to one or more prior conventional synthetic or biologic DMARDs.",
    renal_function = "Baseline eGFR (bedside Schwartz) mean 119 mL/min/1.73 m^2, range 66.5-201.",
    dose_range     = "Oral baricitinib once daily, age-based during the trial: 4 mg QD for ages 9 to <18 years and 2 mg QD for ages 2 to <9 years. The model-based simulations evaluated 1, 2 and 4 mg QD over 10-120 kg and supported the approved weight-based posology of 2 mg QD for 10 to <30 kg and 4 mg QD for >=30 kg.",
    study          = "JUVE-BASIS (NCT03773978), an international phase 3 randomized, double-blind, placebo-controlled withdrawal efficacy and safety trial. Data pooled from the 2-week safety/PK period and the 12-week open-label lead-in period.",
    regions        = "International (multi-region); Japanese versus non-Japanese was screened as a covariate.",
    notes          = "Demographics from Table 1. Of 1377 plasma concentrations collected, 116 (8%) were excluded (106 below the 0.200 ng/mL LLOQ, 6 collected pre-dose or within the lag-time threshold, 4 giving biologically implausible D1 estimates), leaving 1261 records. Safety/PK-period samples were collected as dried whole blood on a Mitra VAMS microsampling device and converted to plasma equivalents with a study-specific blood/plasma ratio of 1.29 (Figure S2, n = 15 concordance pairs), then pooled with the open-label lead-in plasma samples. Estimation used NONMEM 7.4.2 with PsN 4.8.1; parameter estimates from the previous adult rheumatoid-arthritis analysis were used as priors."
  )

  ini({
    # Structural typical values -- Decker 2024 Table 2, "Population mean (%SEE)" column.
    # The reference individual is 74 kg with a baseline eGFR of 93 mL/min/1.73 m^2.
    ld1        <- log(0.457); label("Zero-order absorption duration D1 (h)")                                     # Table 2: D1 = 0.457 h (%SEE 2.76); bootstrap 0.452 (0.373-0.512)
    boxcox_d1  <- 0.501;      label("Box-Cox shape parameter for the BSV on D1 (unitless)")                      # Table 2: Box-Cox transformation parameter for D1 = 0.501 (%SEE 1.40); bootstrap 0.538 (0.323-0.847)
    ltlag      <- log(0.147); label("Absorption lag time (h)")                                                   # Table 2: LAG = 0.147 h (%SEE 2.64); bootstrap 0.147 (0.144-0.150)
    lcl_nonren <- log(3.36);  label("Apparent non-renal clearance CLnr/F at 74 kg (L/h)")                        # Table 2: CLnr/F = 3.36 L/h (%SEE 4.76); bootstrap 3.32 (3.20-3.46)
    lcl_renal  <- log(6.34);  label("Apparent renal clearance CLr/F at 74 kg, eGFR 93 mL/min/1.73 m^2 (L/h)")    # Table 2: CLr/F = 6.34 L/h (%SEE 2.49); bootstrap 6.41 (6.23-6.52)
    lvc        <- log(92.6);  label("Apparent central volume V1/F at 74 kg (L)")                                 # Table 2: V1/F = 92.6 L (%SEE 1.89); bootstrap 92.6 (90.1-94.5)
    lq         <- log(2.67);  label("Apparent intercompartmental clearance Q at 74 kg (L/h)")                    # Table 2: Q = 2.67 L/h (%SEE 7.64); bootstrap 2.76 (2.62-2.88)
    lvp        <- log(25.6);  label("Apparent peripheral volume V2/F at 74 kg (L)")                              # Table 2: V2/F = 25.6 L (%SEE 12.1); bootstrap 25.1 (21.2-30.1)

    # Allometric exponents, held constant during estimation -- Table 2 rows
    # "Allometric Scaling CL" = 0.75 (FIX) and "Allometric Scaling V" = 1 (FIX).
    # A sensitivity analysis that estimated them instead returned 0.48 and 0.70
    # (Discussion); those values are NOT the final model and are not encoded here.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent of WT/74 on CLr/F, CLnr/F and Q (unitless)")           # Table 2 footnotes b, d: (WTE/74)^0.75
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent of WT/74 on V1/F and V2/F (unitless)")                 # Table 2 footnotes c, e: (WTE/74)^1.00

    # Within-subject change in eGFR acting on the apparent renal clearance arm.
    e_dcrcl_cl_renal <- 0.00586; label("Fractional change in CLr/F per mL/min/1.73 m^2 of (CRCL - CRCL_BASE)")   # Table 2: covariate for change in eGFR on CLr/F = 0.00586 (%SEE 29.0); bootstrap 0.00589 (0.00349-0.00852)

    # Between-subject variability -- Table 2, "BSV (%SEE)" column.
    # Footnote a defines the reported quantity as %CV = (sqrt(exp(OMEGA)-1))*100,
    # so the internal variance is OMEGA = log((%CV/100)^2 + 1).
    # Footnote g states the two covariance entries are already on the omega^2 scale.
    #
    # Lower-triangular block over (CLnr/F, CLr/F, V1/F), row by row:
    #   [1,1] 0.539830  CLnr/F BSV 84.6% CV (%SEE 10.3) -> log(0.846^2 + 1)
    #   [2,1] 0.272     cov(CLnr/F, CLr/F), Table 2 (%SEE 9.41), already omega^2
    #   [2,2] 0.144987  CLr/F  BSV 39.5% CV (%SEE 5.30) -> log(0.395^2 + 1)
    #   [3,1] 0         cov(CLnr/F, V1/F) is not reported in Table 2; fixed to 0
    #   [3,2] -0.00706  cov(CLr/F, V1/F),  Table 2 (%SEE 129), already omega^2
    #   [3,3] 0.108176  V1/F   BSV 33.8% CV (%SEE 23.2) -> log(0.338^2 + 1)
    # The block is positive definite (minimum eigenvalue 0.0059).
    # NOTE: trailing comments must not appear inside this multi-line c(); the
    # rxode2 comment-to-label rewriter turns them into stray semicolons.
    etalcl_nonren + etalcl_renal + etalvc ~ c(
      0.539830,
      0.272,     0.144987,
      fixed(0), -0.00706,  0.108176
    )
    etald1 ~ 1.070026        # D1 BSV 138.4% CV (%SEE 13.1) -> log(1.384^2 + 1) = 1.070026; this is the variance of the pre-Box-Cox eta
    etalq  ~ fixed(0.022545) # Q BSV 15.1% CV, reported as FIX in Table 2 -> log(0.151^2 + 1) = 0.022545
    etalvp ~ 0.277044        # V2/F BSV 56.5% CV (%SEE 25.2) -> log(0.565^2 + 1) = 0.277044

    # Residual error -- Table 2 footnote h identifies the reported value as a
    # standard deviation, on a proportional error model.
    propSd <- 0.361; label("Proportional residual error (fraction)")                                             # Table 2: proportional error = 0.361 (%SEE 2.07); bootstrap 0.360 (0.346-0.377)
  })

  model({
    # Within-subject deviation of renal function from the subject's own baseline --
    # the "change in eGFR from baseline" (deltaeGFR) term of Table 2 footnote f.
    # This is the Wahlby 2004 baseline/difference decomposition of a time-varying
    # covariate; set CRCL = CRCL_BASE to disable the time-varying arm.
    dcrcl <- CRCL - CRCL_BASE

    # Allometric size factors referenced to 74 kg (Table 2 footnotes b-e).
    allocl <- (WT / 74)^e_wt_cl_q
    allov  <- (WT / 74)^e_wt_vc_vp

    # Semi-mechanistic partition of apparent total clearance. Table 2 footnote f:
    #   CLr/F = CLr/F_pop * ((baseline eGFR / 93) + theta_deGFR * deltaeGFR)
    # where 93 mL/min/1.73 m^2 is the median eGFR of the previous adult analysis.
    # Table 2 footnote b applies the 0.75 allometric factor to the CL sum; because
    # both arms carry the same exponent, distributing it over the arms is identical.
    cl_renal  <- exp(lcl_renal + etalcl_renal) *
      (CRCL_BASE / 93 + e_dcrcl_cl_renal * dcrcl) * allocl
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) * allocl
    cl        <- cl_renal + cl_nonren

    q  <- exp(lq  + etalq)  * allocl
    vc <- exp(lvc + etalvc) * allov
    vp <- exp(lvp + etalvp) * allov

    # Box-Cox-transformed BSV on the zero-order absorption duration (Table 2 row
    # "Box-Cox transformation parameter for D1"). NONMEM form:
    #   phi = (exp(eta)^lambda - 1) / lambda ; D1 = TVD1 * exp(phi)
    boxcox_eta_d1 <- ((exp(etald1))^boxcox_d1 - 1) / boxcox_d1
    d1   <- exp(ld1) * exp(boxcox_eta_d1)
    tlag <- exp(ltlag)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Zero-order absorption directly into the central compartment over duration D1,
    # after an absorption lag. Figure S1 of the supplement defines the model as
    # having only a central (C) and a peripheral (Cp) compartment, with D1 the
    # absorption duration -- there is no separate depot state. CL, Q and both
    # volumes are apparent (/F), so bioavailability is folded into them.
    dur(central)  <- d1
    alag(central) <- tlag

    # Doses are in mg and volumes in L, so central/vc is mg/L; x1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
