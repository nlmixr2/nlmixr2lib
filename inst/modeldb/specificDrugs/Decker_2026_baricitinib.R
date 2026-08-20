Decker_2026_baricitinib <- function() {
  description <- "Two-compartment population PK model with zero-order absorption, an absorption lag, and fixed allometric scaling for oral baricitinib in 393 pediatric patients aged 2 to <18 years with moderate-to-severe atopic dermatitis (BREEZE-AD-PEDS, NCT03952559). Apparent total clearance is partitioned semi-mechanistically into an eGFR-dependent apparent renal arm (CLr/F) and an eGFR-independent apparent non-renal arm (CLnr/F). Baseline body weight enters through allometric exponents fixed at 0.75 on all clearance terms (CLr/F, CLnr/F, Q) and at 1 on both volumes (V1/F, V2/F), referenced to 74 kg. The structure was carried unchanged from the adult atopic-dermatitis baricitinib population PK model; no further covariate met the stepwise-covariate-modeling inclusion criteria, so the base model is the final model."
  reference <- paste(
    "Decker RL, Ernest CS II, Radtke DB, Prakash A, Zhang X.",
    "A population pharmacokinetic and exposure-response analysis for baricitinib",
    "in pediatric patients with atopic dermatitis.",
    "Clin Pharmacokinet 2026;65(1):119-131.",
    "doi:10.1007/s40262-025-01563-8.",
    sep = " "
  )
  vignette <- "Decker_2026_baricitinib"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight (weight at study entry, WTE in the source paper).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at the study-entry value (Table S1 footnote a, 'at study entry'). The allometric reference weight is 74 kg (Table 1 footnotes b-e), which is carried over from the adult analysis and is NOT a statistic of this pediatric cohort (whose mean weight is 46.6 kg). Cohort mean 46.6 kg, range 12.0-104 kg (Table S1); 78% of patients weighed >=30 kg. Simulations in the paper spanned 10-120 kg.",
      source_name        = "WTE"
    ),
    CRCL_BASE = list(
      description        = "Baseline estimated glomerular filtration rate (bedside Schwartz equation), BSA-normalized and time-fixed per subject.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters CLr/F as the ratio (CRCL_BASE / 93); 93 mL/min/1.73 m^2 is the median eGFR of the previous adult population PK analysis and is the reference at which the reported CLr/F typical value of 7.9 L/h applies (Table 1 footnote f). Estimated with the bedside Schwartz equation (Methods, 'PK Model Development'), a creatinine-based pediatric eGFR already expressed per 1.73 m^2. Cohort mean 109 mL/min/1.73 m^2, range 54.3-196 (Table S1); by weight group the mean is 125 for <30 kg and 105 for >=30 kg.",
      source_name        = "baseline eGFR"
    ),
    CRCL = list(
      description        = "Time-varying estimated glomerular filtration rate (bedside Schwartz equation), BSA-normalized.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used only through the within-subject deviation (CRCL - CRCL_BASE), which is the 'change in eGFR from baseline' (deltaeGFR) term of Table 1 footnote f -- the Wahlby 2004 baseline/difference (BCOV/DCOV) decomposition of a time-varying covariate. Set CRCL = CRCL_BASE to switch the time-varying arm off; the model then reduces to baseline-eGFR-only scaling. The source paper does not tabulate the distribution of on-treatment eGFR changes.",
      source_name        = "eGFR"
    )
  )

  # Screened in the stepwise covariate model (PsN SCM) but not retained: no
  # covariate met the forward-inclusion criteria (p <= 0.01, i.e. dOFV >= 6.635,
  # AND a >=5% reduction in the BSV variance of the affected parameter), so
  # backward exclusion was never run and the base model became the final model
  # (Methods, "PK Model Development and E-R Relationship"). Documented here to
  # preserve the provenance of the covariate screen without declaring unused
  # covariates.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at study entry.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as the only continuous covariate and not retained. The Discussion states explicitly that age was not a significant covariate after the effect of weight was accounted for; Figure S6 shows the high age-weight correlation. Cohort mean 12.0 years, range 2.00-17.9 (Table S1); 72% of patients were 10 to <18 years."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units       = "(binary)",
      type        = "binary",
      reference_category = "male",
      notes       = "Screened as 'gender' in the stepwise covariate model and not retained. Cohort: 197 of 392 female (50%), 195 male (Table S1)."
    ),
    RACE_WHITE = list(
      description = "White race indicator.",
      units       = "(binary)",
      type        = "binary",
      reference_category = "non-White",
      notes       = "Screened as part of the categorical 'race' covariate and not retained. 287 of 392 White (73%) (Table S1)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator.",
      units       = "(binary)",
      type        = "binary",
      reference_category = "non-Asian",
      notes       = "Screened as part of the categorical 'race' covariate and not retained. 76 of 392 Asian (19%) (Table S1)."
    ),
    RACE_BLACK = list(
      description = "Black race indicator.",
      units       = "(binary)",
      type        = "binary",
      reference_category = "non-Black",
      notes       = "Screened as part of the categorical 'race' covariate and not retained. 11 of 392 Black (3%) (Table S1)."
    )
  )

  # Table S1 also counts 8 Native American, 1 Multiple and 9 Missing among the
  # 392 patients. Those three strata are part of the same screened-and-rejected
  # 'race' covariate; they have no separate canonical covariate-column name in
  # inst/references/covariate-columns.md and are not used by this model, so they
  # are recorded here in prose rather than as covariatesDataExcluded entries.

  compartmentData <- list(
    central     = list(analyte = "baricitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "baricitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 393L,
    n_studies      = 1L,
    n_observations = 2035L,
    age_range      = "2.00-17.9 years (enrolment criterion 2 to <18 years)",
    age_mean       = "12.0 years",
    weight_range   = "12.0-104 kg",
    weight_mean    = "46.6 kg",
    sex_female_pct = 50.3,
    race_ethnicity = c(White = 73.2, Asian = 19.4, Black = 2.8, `Native American` = 2.0, Multiple = 0.3, Missing = 2.3),
    disease_state  = "Moderate-to-severe atopic dermatitis (validated Investigator Global Assessment for Atopic Dermatitis score >=3, Eczema Area and Severity Index >=16, body surface area involvement >=10%) with inadequate response to moderate or higher potency topical corticosteroids and inadequate response or intolerance to topical calcineurin inhibitors, or inadequate response to systemic treatment.",
    renal_function = "Baseline eGFR (bedside Schwartz) mean 109 mL/min/1.73 m^2, range 54.3-196; by weight group, mean 125 for <30 kg and 105 for >=30 kg. Screening eGFR was required to exceed 60 mL/min/1.73 m^2.",
    dose_range     = "Oral baricitinib once daily. Open-label PK lead-in used age-based high doses (4 mg QD for 10 to <18 years, 2 mg QD for 2 to <10 years). The 16-week double-blind period randomised 1:1:1:1 to placebo or low/medium/high dose, which in absolute terms was 1, 2 or 4 mg QD for ages 10 to <18 years and 0.5, 1 or 2 mg QD for ages 2 to <10 years. The model-based simulations evaluated 10-120 kg and supported the approved weight-based posology of 2 mg QD for 10 to <30 kg and 4 mg QD for >=30 kg.",
    study          = "BREEZE-AD-PEDS (NCT03952559), a multicentre phase 3 randomized, double-blind, placebo-controlled trial with an open-label PK lead-in period. Concentration data from both the lead-in and the double-blind treatment periods were pooled for the population PK analysis.",
    regions        = "Multicentre international.",
    notes          = "The population PK dataset held 2035 concentrations from 393 participants; one participant was later excluded from the reported demographics because of an erroneous screening weight, so Table S1 summarises 392. Lead-in samples were collected as dried whole blood on a Mitra VAMS microsampling device and converted to plasma equivalents with a study-specific blood/plasma ratio of 1.32 (Figure S7, 15 concordance pairs from 15 participants; 214 whole-blood samples from 33 participants converted), then pooled with the plasma samples. Assay LLOQ 0.200 ng/mL, ULOQ 200 ng/mL. Estimation used NONMEM 7.5.0 with PsN 4.7.0; parameter estimates from the previous adult atopic-dermatitis analysis were used as priors. A cluster of implausible concentrations attributed to mis-recorded dose times was retained in the final analysis after a sensitivity refit showed minimal parameter impact (Discussion). The exposure-response analysis for vIGA-AD 0/1 at week 16 was a descriptive exposure-quartile summary in 514 patients; it reports no fitted parameters and is therefore not encoded in this model."
  )

  ini({
    # Structural typical values -- Decker 2026 Table 1, "Population mean (%SEE)" column.
    # The reference individual is 74 kg with a baseline eGFR of 93 mL/min/1.73 m^2.
    ld1        <- log(0.263); label("Zero-order absorption duration D1 (h)")                                     # Table 1: D1 = 0.263 h (%SEE 8.86); bootstrap 0.253 (0.232-0.295)
    boxcox_d1  <- 0.311;      label("Box-Cox shape parameter for the BSV on D1 (unitless)")                      # Table 1: Box-Cox transformation parameter for D1 = 0.311 (%SEE 6.85); bootstrap 0.395 (0.167-0.481)
    ltlag      <- log(0.144); label("Absorption lag time (h)")                                                   # Table 1: LAG = 0.144 h (%SEE 1.73); bootstrap 0.146 (0.135-0.151)
    lcl_nonren <- log(2.76);  label("Apparent non-renal clearance CLnr/F at 74 kg (L/h)")                        # Table 1: CLnr/F = 2.76 L/h (%SEE 15.4); bootstrap 2.68 (2.51-3.03)
    lcl_renal  <- log(7.9);   label("Apparent renal clearance CLr/F at 74 kg, eGFR 93 mL/min/1.73 m^2 (L/h)")    # Table 1: CLr/F = 7.9 L/h (%SEE 12.1); bootstrap 8.02 (7.44-8.38)
    lvc        <- log(119);   label("Apparent central volume V1/F at 74 kg (L)")                                 # Table 1: V1/F = 119 L (%SEE 4.12); bootstrap 119 (116-122)
    lq         <- log(2.4);   label("Apparent intercompartmental clearance Q at 74 kg (L/h)")                    # Table 1: Q = 2.4 L/h (%SEE 9.17); bootstrap 2.52 (2.19-2.63)
    lvp        <- log(46.8);  label("Apparent peripheral volume V2/F at 74 kg (L)")                              # Table 1: V2/F = 46.8 L (%SEE 14.7); bootstrap 46.5 (34.3-60.5)

    # Allometric exponents, held constant during estimation -- Table 1 rows
    # "Allometric scaling CL" = 0.75 (FIX) and "Allometric scaling V" = 1 (FIX).
    # A sensitivity analysis that estimated them instead returned 0.59 for
    # clearance and 0.66 for volume (Discussion); those values are NOT the final
    # model and are not encoded here.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent of WT/74 on CLr/F, CLnr/F and Q (unitless)")           # Table 1 footnotes b, d: (WTE/74)^0.75
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent of WT/74 on V1/F and V2/F (unitless)")                 # Table 1 footnotes c, e: (WTE/74)^1.00

    # Within-subject change in eGFR acting on the apparent renal clearance arm.
    e_dcrcl_cl_renal <- 0.00778; label("Fractional change in CLr/F per mL/min/1.73 m^2 of (CRCL - CRCL_BASE)")   # Table 1: covariate for change in eGFR on CLr/F = 0.00778 (%SEE 42.5); bootstrap 0.00738 (0.00462-0.0115)

    # Between-subject variability. The variances below are the "Estimate" column
    # of ESM Table S2 ("Correlation Matrix for the Final Model"), which tabulates
    # the OMEGA variances directly together with their RSE% -- the RSE% values
    # match Table 1's %SEE column for every BSV term (13.9, 4.54, 16.5, 19.9,
    # 44.3, fixed), and the square roots of the variances reproduce Table S2's
    # own diagonal of standard deviations (1.14, 0.559, 0.598, 0.127, 0.921,
    # 0.150) exactly. Table S2's off-diagonal entries are CORRELATIONS; combined
    # with those standard deviations they reproduce Table 1's two reported
    # covariances (0.904*0.5595*0.5983 = 0.3026 ~ 0.303, and
    # -0.349*0.5983*0.1269 = -0.02650 ~ -0.0265). Table S2 also gives
    # corr(CLnr/F, V1/F) = 0.00 explicitly, so that block element is fixed to 0.
    #
    # NOTE: Table 1's "BSV a (%CV)" column is NOT used to derive these variances.
    # Footnote a defines %CV = (sqrt(exp(OMEGA)-1))*100, which round-trips for D1,
    # V1/F, Q and V2/F but NOT for CLnr/F (58.4% implies 0.293, Table S2 says
    # 0.313) or CLr/F (62.3% implies 0.328, Table S2 says 0.358). Reconstructing
    # the block from the Table 1 %CV values while keeping Table 1's own
    # covariances yields a matrix that is not positive definite (minimum
    # eigenvalue -0.0073), i.e. an impossible covariance matrix, whereas the
    # Table S2 variances give a positive definite block. See the vignette
    # "Assumptions and deviations" section.
    #
    # Lower-triangular block over (CLnr/F, CLr/F, V1/F), row by row:
    #   [1,1]  0.313    CLnr/F variance, Table S2 (RSE 4.54)
    #   [2,1]  0.303    cov(CLnr/F, CLr/F), Table 1 (%SEE 10.4), already omega^2
    #   [2,2]  0.358    CLr/F  variance, Table S2 (RSE 16.5)
    #   [3,1]  0        corr(CLnr/F, V1/F) = 0.00 in Table S2; fixed to zero
    #   [3,2] -0.0265   cov(CLr/F, V1/F),  Table 1 (%SEE 15.8), already omega^2
    #   [3,3]  0.0161   V1/F   variance, Table S2 (RSE 19.9)
    # The block is positive definite (minimum eigenvalue 0.0038).
    # NOTE: trailing comments must not appear inside this multi-line c(); the
    # rxode2 comment-to-label rewriter turns them into stray semicolons.
    etalcl_nonren + etalcl_renal + etalvc ~ c(
      0.313,
      0.303,     0.358,
      fixed(0), -0.0265,  0.0161
    )
    etald1 ~ 1.30            # D1 BSV variance, ESM Table S2 (RSE 13.9); this is the variance of the pre-Box-Cox eta
    etalq  ~ fixed(0.0225)   # Q BSV variance, ESM Table S2 (SD diagonal 0.150); Table 1 reports the same term as 15.1 percent CV
    etalvp ~ 0.848           # V2/F BSV variance, ESM Table S2 (RSE 44.3)

    # Residual error -- Table 1 footnote h identifies the reported value as a
    # standard deviation, on a proportional error model.
    propSd <- 0.427; label("Proportional residual error (fraction)")                                             # Table 1: proportional error = 0.427 (%SEE 13.9); bootstrap 0.426 (0.411-0.441)
  })

  model({
    # Within-subject deviation of renal function from the subject's own baseline --
    # the "change in eGFR from baseline" (deltaeGFR) term of Table 1 footnote f.
    # This is the Wahlby 2004 baseline/difference decomposition of a time-varying
    # covariate; set CRCL = CRCL_BASE to disable the time-varying arm.
    dcrcl <- CRCL - CRCL_BASE

    # Allometric size factors referenced to 74 kg (Table 1 footnotes b-e).
    allocl <- (WT / 74)^e_wt_cl_q
    allov  <- (WT / 74)^e_wt_vc_vp

    # Semi-mechanistic partition of apparent total clearance. Table 1 footnote f:
    #   CLr/F = CLr/F_pop * ((baseline eGFR / 93) + theta_deGFR * deltaeGFR)
    # where 93 mL/min/1.73 m^2 is the median eGFR of the previous adult analysis.
    # Table 1 footnote b applies the 0.75 allometric factor to the CL sum; because
    # both arms carry the same exponent, distributing it over the arms is identical.
    cl_renal  <- exp(lcl_renal + etalcl_renal) *
      (CRCL_BASE / 93 + e_dcrcl_cl_renal * dcrcl) * allocl
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) * allocl
    cl        <- cl_renal + cl_nonren

    q  <- exp(lq  + etalq)  * allocl
    vc <- exp(lvc + etalvc) * allov
    vp <- exp(lvp + etalvp) * allov

    # Box-Cox-transformed BSV on the zero-order absorption duration (Table 1 row
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
    # absorption duration and LAG the absorption lag -- there is no separate depot
    # state. CL, Q and both volumes are apparent (/F), so bioavailability is folded
    # into them.
    dur(central)  <- d1
    alag(central) <- tlag

    # Doses are in mg and volumes in L, so central/vc is mg/L; x1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
