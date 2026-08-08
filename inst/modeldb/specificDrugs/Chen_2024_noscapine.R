Chen_2024_noscapine <- function() {
  description <- paste(
    "Semi-mechanistic population PK model for oral noscapine in healthy adults",
    "(Chen 2024) with zero-order release into a four-compartment transit chain,",
    "an explicit liver compartment that carries the profound hepatic first-pass",
    "extraction (apparent clearance CL/F leaves the liver, not the plasma), and",
    "three-compartment systemic disposition. Liver plasma flow Qh (55 L/h) and",
    "liver volume Vh (1.5 L) are held fixed at physiologic values. Apparent",
    "clearance carries a CYP2C9 genotype-predicted-phenotype effect on three",
    "levels keyed to the CPIC activity score (extensive metabolizer reference",
    "958 L/h; intermediate metabolizer with activity score 1.5, 531 L/h; poor",
    "and intermediate metabolizers with activity score 1.0, 343 L/h) plus a",
    "total-body-weight power term (exponent 1.34, 77.3 kg reference);",
    "inter-compartmental clearance to the first peripheral compartment carries",
    "an age power term (exponent 0.348, 29 year reference). Relative",
    "bioavailability of the reformulated test suspension is 82.8% of the",
    "reference suspension. Correlated inter-individual variability on absorption",
    "duration, CL/F, Vc/F and F1; inter-occasion variability on the transit rate",
    "and on Qp1/F across the two crossover periods; exponential residual error.",
    sep = " "
  )
  reference <- paste(
    "Chen Z, Taubert M, Chen C, Boland J, Dong Q, Bilal M, Dokos C, Wachall B,",
    "Wargenau M, Scheidel B, Wiesen MHJ, Schaeffeler E, Tremmel R, Schwab M,",
    "Fuhr U (2024). A Semi-Mechanistic Population Pharmacokinetic Model of",
    "Noscapine in Healthy Subjects Considering Hepatic First-Pass Extraction",
    "and CYP2C9 Genotypes. Drugs in R&D 24:187-199.",
    "doi:10.1007/s40268-024-00466-6.",
    sep = " "
  )
  vignette <- "Chen_2024_noscapine"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  # Declared explicitly: the dose enters the FIRST TRANSIT compartment, not a
  # `depot` or `central`, so buildModelDb()'s depot/central heuristic would
  # otherwise mislabel this model as central-dosed. Dose records must target
  # cmt = "transit1" with rate = -2 (zero-order, dur(transit1) = D1).
  dosing <- "transit1"

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline total body weight (TBW). Drives a power model on apparent",
        "clearance normalised to the cohort MEAN of 77.3 kg (Sect. 3.3.2",
        "equation: CL/F = (TBW/77.3)^1.34 * TVCL_phenotype/F * exp(eta)).",
        "Note the paper centres on the mean 77.3 kg, not the median 77.0 kg",
        "(Table 1). The exponent 1.34 was ESTIMATED (Table 4, RSE 35.1%,",
        "SIR 95% CI 0.485-2.33), not fixed at an allometric 0.75.",
        "Table 1 range 54.3-107 kg."
      ),
      source_name        = "TBW"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline age. Drives a power model on the inter-compartmental",
        "clearance between the central and first peripheral compartment,",
        "normalised to the cohort MEDIAN of 29.0 years (Sect. 3.3.2 equation:",
        "Qp1/F = (Age/29.0)^0.348 * TVQp1/F * exp(eta)). The paper centres on",
        "the median 29 years, not the mean 33 years (Table 1). Exponent 0.348",
        "ESTIMATED (Table 4, RSE 32.3%, SIR 95% CI 0.125-0.568).",
        "Table 1 range 19.0-65.0 years."
      ),
      source_name        = "Age"
    ),
    CYP2C9_IM_AS15 = list(
      description        = paste(
        "CYP2C9 intermediate-metabolizer with CPIC activity score 1.5",
        "(genotype *1/*2) indicator"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (extensive metabolizer, *1/*1 or *1/*9, activity score 2.0) when",
        "CYP2C9_PM_IM_AS10 is also 0"
      ),
      notes              = paste(
        "Mutually exclusive with CYP2C9_PM_IM_AS10; both 0 identifies the",
        "extensive-metabolizer reference stratum whose typical CL/F is the",
        "lcl typical value. Table 2: 10/48 subjects (20.8%) are *1/*2.",
        "Additive log shift on CL/F: log(531/958) = -0.5901 (Table 4",
        "CL/F_IM_1.5 = 531 L/h vs CL/F_EM = 958 L/h)."
      ),
      source_name        = "CYP2C9 genotype-predicted phenotype (IM with AS of 1.5)"
    ),
    CYP2C9_PM_IM_AS10 = list(
      description        = paste(
        "CYP2C9 poor-metabolizer or intermediate-metabolizer with CPIC",
        "activity score 1.0 (genotypes *1/*3, *2/*3, *3/*3) indicator"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (extensive metabolizer, *1/*1 or *1/*9, activity score 2.0) when",
        "CYP2C9_IM_AS15 is also 0"
      ),
      notes              = paste(
        "Mutually exclusive with CYP2C9_IM_AS15. Table 2: 6/48 subjects",
        "(12.6%) carry *1/*3 (n = 3), *2/*3 (n = 2) or *3/*3 (n = 1).",
        "Sect. 3.1: activity-score-1.0 intermediate metabolizers were pooled",
        "with poor metabolizers 'because of the limited sample size of",
        "homozygous carriers'. This is a strict SUBSET of the pooled",
        "CYP2C9_PM_IM canonical, which would additionally capture the",
        "activity-score-1.5 *1/*2 subjects carried here by CYP2C9_IM_AS15.",
        "Additive log shift on CL/F: log(343/958) = -1.0271 (Table 4",
        "CL/F_PM_IM_1.0 = 343 L/h vs CL/F_EM = 958 L/h)."
      ),
      source_name        = "CYP2C9 genotype-predicted phenotype (PM & IM with AS of 1)"
    ),
    FORM_NOSCAPINE_TEST = list(
      description        = paste(
        "Reformulated noscapine oral suspension (InfectoPharm, the bioequi-",
        "valence study test product) vs the Nipaxon 5 mg/mL oral suspension",
        "(McNeil, reference) formulation indicator"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (Nipaxon 5 mg/mL reference oral suspension; the relative-",
        "bioavailability anchor F = 1)"
      ),
      notes              = paste(
        "Per-dose-record indicator: this is a two-period crossover, so every",
        "subject receives the test product on one period and the reference on",
        "the other (Sect. 2.1). Set on each dose row. Drives relative",
        "bioavailability only: F = exp((lf1 + etalf1) * FORM_NOSCAPINE_TEST),",
        "so reference doses get F = 1 exactly and test doses get 82.8% with",
        "34.1% CV inter-individual variability (Table 4 F1). Sect. 3.3.2:",
        "'only apparent bioavailability differed between test and reference",
        "preparations' -- no absorption-rate parameter differs by formulation."
      ),
      source_name        = "formulation (test vs reference)"
    ),
    OCC = list(
      description        = "Integer-valued crossover-period occasion indicator",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "OCC = 1 for the first treatment period and OCC = 2 for the second",
        "(two-period crossover with a 6-14 day washout, Sect. 2.1).",
        "Decomposed inside model() into binary indicators oc1 / oc2 that",
        "multiplex the inter-occasion-variability etas on the log transit rate",
        "and on log Qp1/F (Table 4 IOV rows), following the",
        "Chen_2023_nemonoxacin.R and Jonsson_2011_ethambutol.R precedent.",
        "For single-occasion simulation pass OCC = 1."
      ),
      source_name        = "period"
    )
  )

  covariatesDataExcluded <- list(
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Screened in the forward-addition covariate analysis (Sect. 2.4.2)",
        "but not retained in the final model (Sect. 3.3.2: 'Because of a lack",
        "of significant effects, other clinical characteristics ... were not",
        "incorporated'). Table 1 mean 177 cm, range 159-197 cm."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Screened (Sect. 2.4.2) but not retained (Sect. 3.3.2).",
        "Table 1 mean 24.6, range 18.0-29.5 kg/m^2."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Liver-function marker screened (Sect. 2.4.2) but not retained",
        "(Sect. 3.3.2). Table 1 mean 24.1, range 7.0-96.0 U/L."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Liver-function marker screened (Sect. 2.4.2) but not retained",
        "(Sect. 3.3.2). Table 1 mean 23.5, range 14.0-39.0 U/L."
      )
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Liver-function marker screened (Sect. 2.4.2) but not retained",
        "(Sect. 3.3.2). Table 1 mean 64.8, range 36.0-103 U/L."
      )
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate (CKD-EPI 2009)",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Renal-function marker screened (Sect. 2.4.2) but not retained",
        "(Sect. 3.3.2). Table 1 mean 107, range 76.0-132 mL/min."
      )
    ),
    UREA = list(
      description = "Plasma urea concentration",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Renal-function marker screened (Sect. 2.4.2) but not retained",
        "(Sect. 3.3.2). Table 1 mean 26.7, range 13.0-40.0 mg/dL."
      )
    ),
    CYP2C19_PM_IM = list(
      description = "CYP2C19 intermediate-metabolizer indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Genotyped in all 48 subjects and screened on CL/F (Sect. 2.4.2,",
        "Table 2: 13/48 IM), but no significant effect was found -- Fig. 2",
        "shows no significant Cmax / AUC difference and Sect. 3.3.2 retained",
        "only CYP2C9. Not in the final model."
      )
    ),
    CYP3A4_IM = list(
      description = "CYP3A4 intermediate-metabolizer (*1/*22) indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Genotyped and screened (Table 2: 9/48 are *1/*22) but not retained",
        "(Sect. 3.3.2, Fig. 2)."
      )
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser (*1 carrier) indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Genotyped and screened (Table 2: 12/48 *1/*3 intermediate",
        "metabolizers, 36/48 poor metabolizers, no *1/*1) but not retained",
        "(Sect. 3.3.2, Fig. 2)."
      )
    )
  )

  compartmentData <- list(
    transit1    = list(analyte = "noscapine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "noscapine", units = "mg", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "noscapine", units = "mg", specimen = "administration site", verified = TRUE),
    transit4    = list(analyte = "noscapine", units = "mg", specimen = "administration site", verified = TRUE),
    liver       = list(analyte = "noscapine", units = "mg", specimen = "tissue", verified = TRUE),
    central     = list(analyte = "noscapine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "noscapine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "noscapine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 48L,
    n_studies      = 1L,
    n_observations = 1920L,
    age_range      = "19.0-65.0 years",
    age_median     = "29 years",
    weight_range   = "54.3-107 kg",
    weight_median  = "77.0 kg",
    sex_female_pct = 37.5,
    race_ethnicity = paste(
      "Predominantly European ancestry by genome-wide principal-component",
      "analysis; 2 subjects clustered with South Asian / American populations",
      "and 2 with African populations (Fig. S1 of the ESM)"
    ),
    disease_state  = "healthy volunteers",
    dose_range     = "50 mg single oral dose (5 mg/mL oral suspension) in each of two crossover periods",
    regions        = "Germany (Clinical Pharmacology Unit, University Hospital of Cologne)",
    bmi_range      = "18.0-29.5 kg/m^2",
    genotype       = paste(
      "CYP2C9: 32 extensive metabolizers (*1/*1 n = 31, *1/*9 n = 1),",
      "10 intermediate metabolizers with activity score 1.5 (*1/*2),",
      "6 poor / intermediate metabolizers with activity score 1.0",
      "(*1/*3 n = 3, *2/*3 n = 2, *3/*3 n = 1). CYP2C19, CYP3A4, CYP3A5 and",
      "CYP2E1 were also genotyped (Table 2) but had no significant PK effect."
    ),
    study_design   = paste(
      "Randomized, two-period, two-stage, crossover bioequivalence study of a",
      "reformulated noscapine oral suspension (test) vs Nipaxon 5 mg/mL oral",
      "suspension (reference), 6-14 day washout. DRKS00017760;",
      "EUDRA-CT 2019-002012-12."
    ),
    notes          = paste(
      "Demographics and baseline covariates in Table 1 (n = 48; 30 men, 18",
      "women). Genotype frequencies in Table 2. Rich sampling at baseline and",
      "0.17, 0.33, 0.5, 0.67, 0.83, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 6,",
      "9, 12, 16 and 24 h post-dose (Sect. 2.1). Concentrations below the",
      "0.100 ug/L LLOQ were set to zero before Cmax (96 records) and discarded",
      "after Cmax (17 records)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Absorption (Table 4). Zero-order release of the dose into the first
    # of four sequential transit compartments (Sect. 3.3.1: 'Zero-order
    # absorption described the absorption process best, while using
    # transit compartments instead of simply introducing a lag time
    # better described the delayed gastric emptying ... Manually
    # including four sequential transit compartments proved to be the
    # most robust approach'). Fig. 1 labels the release arm 'Rate
    # (equivalent to dose/absorption duration)' and the transit chain
    # transfer 'TR average transit rate'.
    # ------------------------------------------------------------------
    ld1  <- log(0.0624); label("Log zero-order absorption duration D1 (h)")                                        # Table 4: D1 = 0.0624 h (RSE 28.6%; SIR median 0.0663, 95% CI 0.0353-0.105)
    lktr <- log(13.0)  ; label("Log transit rate constant TR (1/h)")                                              # Table 4: TR = 13.0 1/h (RSE 6.4%; SIR median 13.0, 95% CI 11.4-14.7)

    # ------------------------------------------------------------------
    # Apparent clearance (Table 4). Three CYP2C9 genotype-predicted
    # phenotype strata were estimated as three separate typical values;
    # the extensive metabolizer (activity score 2.0) value is the
    # reference and the other two strata enter as additive log shifts
    # (Kleideiter_2017_cebranopadol.R precedent). Note the CL/F values
    # are NOT proportional to the activity score (531/958 = 0.554 and
    # 343/958 = 0.358 vs activity-score ratios 0.75 and 0.50), so a
    # continuous activity-score covariate cannot reproduce the model --
    # discrete indicators are required.
    # ------------------------------------------------------------------
    lcl <- log(958); label("Log apparent clearance CL/F in CYP2C9 extensive metabolizers (L/h)")                   # Table 4: CL/F_EM = 958 L/h (RSE 10.3%; SIR median 939, 95% CI 747-1144); *1/*1 and *1/*9

    e_cyp2c9im_as15_lcl   <- log(531 / 958); label("Log-shift on CL/F for CYP2C9 intermediate metabolizers, activity score 1.5 (unitless)")        # Table 4: CL/F_IM_1.5 = 531 L/h (RSE 23.0%) vs CL/F_EM 958 L/h -> log(531/958) = -0.5901
    e_cyp2c9pmim_as10_lcl <- log(343 / 958); label("Log-shift on CL/F for CYP2C9 poor and intermediate metabolizers, activity score 1.0 (unitless)") # Table 4: CL/F_PM_IM_1.0 = 343 L/h (RSE 24.8%) vs CL/F_EM 958 L/h -> log(343/958) = -1.0271

    # ------------------------------------------------------------------
    # Systemic disposition (Table 4). Three-compartment: one central and
    # two apparent peripheral compartments (Sect. 3.3.1). All volumes and
    # inter-compartmental clearances are apparent (/F).
    # ------------------------------------------------------------------
    lvc  <- log(49.7); label("Log apparent central volume Vc/F (L)")                                              # Table 4: Vc/F = 49.7 L (RSE 5.9%; SIR median 49.6, 95% CI 43.7-55.2)
    lq   <- log(18.9); label("Log apparent inter-compartmental clearance to peripheral1, Qp1/F (L/h)")             # Table 4: Qp1/F = 18.9 L/h (RSE 4.0%; SIR median 18.8, 95% CI 17.3-20.3)
    lvp  <- log(243) ; label("Log apparent first peripheral volume Vp1/F (L)")                                     # Table 4: Vp1/F = 243 L (RSE 3.6%; SIR median 243, 95% CI 227-261)
    lq2  <- log(25.8); label("Log apparent inter-compartmental clearance to peripheral2, Qp2/F (L/h)")             # Table 4: Qp2/F = 25.8 L/h (RSE 4.3%; SIR median 25.8, 95% CI 23.6-27.8)
    lvp2 <- log(31.2); label("Log apparent second peripheral volume Vp2/F (L)")                                    # Table 4: Vp2/F = 31.2 L (RSE 4.5%; SIR median 31.2, 95% CI 28.4-33.9)

    # ------------------------------------------------------------------
    # Relative bioavailability of the reformulated test suspension
    # (Table 4 F1 = 82.8%). The reference Nipaxon suspension is the
    # anchor at F = 1, so this term applies only when
    # FORM_NOSCAPINE_TEST = 1.
    # ------------------------------------------------------------------
    lf1 <- log(0.828); label("Log relative bioavailability F1 of the test suspension vs reference (fraction)")     # Table 4: F1 = 82.8% (RSE 5.4%; SIR median 83.1, 95% CI 75.2-92.2)

    # ------------------------------------------------------------------
    # Semi-physiologic liver compartment. Both constants were held FIXED
    # at physiologic values, chosen after a sensitivity analysis
    # (Sect. 3.3.3: 'the simpler relationship between liver and central
    # compartments, and predefined values of 1.5 L and 55 L/h for liver
    # volume and Qh, respectively, were applied in the structural model';
    # also Fig. 1 caption 'Qh liver plasma flow (fixed at 55 L/h)' and
    # 'Vh liver volume (fixed at 1.5 L)'). Table S2 of the ESM reports
    # that CL/F, F1, IIV and residual variability were nearly unchanged
    # across alternative liver settings.
    #
    # NOTE these are TRUE physiologic constants while the plasma-side
    # volumes and clearances are APPARENT (/F). Mixing the two scales is
    # the paper's own construction: the liver state therefore holds an
    # apparent amount, and the liver arm sets the SHAPE of the profile
    # while total exposure remains AUC = Dose * F / (CL/F). See the
    # vignette's Assumptions and deviations section.
    # ------------------------------------------------------------------
    lqh     <- fixed(log(55)) ; label("Log liver plasma flow Qh (L/h)")                                            # Fig. 1 caption and Sect. 3.3.3: Qh FIXED at 55 L/h
    lvliver <- fixed(log(1.5)); label("Log liver volume Vh (L)")                                                   # Fig. 1 caption and Sect. 3.3.3: Vh FIXED at 1.5 L

    # ------------------------------------------------------------------
    # Covariate exponents (Table 4 'Covariates'; equations in
    # Sect. 3.3.2). Both were ESTIMATED, not fixed at allometric values.
    #   CL/F  = (TBW/77.3)^1.34  * TVCL_phenotype/F * exp(eta2)
    #   Qp1/F = (Age/29.0)^0.348 * TVQp1/F          * exp(eta1)
    # 77.3 kg is the cohort MEAN total body weight and 29.0 years the
    # cohort MEDIAN age (Table 1).
    # ------------------------------------------------------------------
    e_wt_cl <- 1.34 ; label("Power exponent on CL/F with total body weight, 77.3 kg reference (unitless)")         # Table 4: TBW on CL/F = 1.34 (RSE 35.1%; SIR median 1.35, 95% CI 0.485-2.33)
    e_age_q <- 0.348; label("Power exponent on Qp1/F with age, 29.0 year reference (unitless)")                    # Table 4: Age on Qp1/F = 0.348 (RSE 32.3%; SIR median 0.343, 95% CI 0.125-0.568)

    # ------------------------------------------------------------------
    # Correlated inter-individual variability (Table 4). Footnote a:
    # 'IIV is reported as CV%'; the model is exponential
    # (Sect. 2.4.1: theta_i = theta * exp(eta_i)), so the reported CV%
    # back-transforms as omega^2 = log(1 + CV^2). Footnote b: 'Correlation
    # is shown as covariance between omega^2', so the three reported
    # off-diagonals are covariances on the log scale and enter the block
    # directly. Covariances the paper does not report were not in the
    # NONMEM $OMEGA BLOCK and are structurally zero -> fixed(0).
    #
    # The CV%-as-log-normal-CV reading (rather than CV% = omega * 100)
    # is confirmed against the paper's own Table 3: it predicts a
    # between-subject AUC CV of 77.8% vs the observed 78.1-79.0%,
    # whereas the omega = CV reading overshoots to 82.3%. The gate is
    # reproduced in the vignette.
    # ------------------------------------------------------------------
    # Per-element source trace (comments must stay OUTSIDE the c(...) -- a `#`
    # inside an omega block breaks rxode2's comment-to-label rewriter):
    #   diag row 1  1.809534  Table 4 IIV D1    = 226%  CV (RSE 40.6%, shrinkage 19.1%) -> log(1 + 2.26^2)
    #   diag row 2  0.283922  Table 4 IIV CL/F  = 57.3% CV (RSE 18.8%, shrinkage  0.1%) -> log(1 + 0.573^2)
    #   diag row 3  0.078095  Table 4 IIV Vc/F  = 28.5% CV (RSE 15.5%, shrinkage  5.2%) -> log(1 + 0.285^2)
    #   diag row 4  0.110003  Table 4 IIV F1    = 34.1% CV (RSE 16.4%, shrinkage  1.2%) -> log(1 + 0.341^2)
    #   (2,1)       0.109     Table 4 covariance(D1, CL/F) = 0.109  (RSE 225%)
    #   (3,1)      -0.175     Table 4 covariance(D1, V/F)  = -0.175 (RSE 92.6%)
    #   (4,1)       0.162     Table 4 covariance(D1, F1)   = 0.162  (RSE 131%)
    #   (3,2), (4,2), (4,3) = fixed(0): not reported, i.e. not in the $OMEGA BLOCK
    etald1 + etalcl + etalvc + etalf1 ~
      c(1.809534,
        0.109,     0.283922,
        -0.175,    fixed(0),   0.078095,
        0.162,     fixed(0),   fixed(0),   0.110003)

    # ------------------------------------------------------------------
    # Inter-occasion variability across the two crossover periods
    # (Table 4 'Inter-occasion variability (%)'). Reported on the same
    # CV% scale as the IIVs, so converted identically. nlmixr2 has no
    # NONMEM `$OMEGA BLOCK(1) SAME` shortcut, so occasion 2 gets its own
    # eta with the variance fixed equal to the occasion-1 estimate
    # (Chen_2023_nemonoxacin.R / Jonsson_2011_ethambutol.R precedent).
    # ------------------------------------------------------------------
    etaiov_ktr_1 ~ 0.228795      # Table 4: IOV TR = 50.7% CV (RSE 16.5%, shrinkage 2.7%) -> log(1 + 0.507^2)
    etaiov_ktr_2 ~ fix(0.228795) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_q_1   ~ 0.088404      # Table 4: IOV Qp1/F = 30.4% CV (RSE 13.2%, shrinkage 14.5%) -> log(1 + 0.304^2)
    etaiov_q_2   ~ fix(0.088404) # SAME-equivalent: equal to the occasion-1 IOV variance

    # ------------------------------------------------------------------
    # Residual error. Sect. 2.4.1: 'An additive residual error model was
    # applied to the log-transformed observations versus time data
    # (equivalent to an exponential error model with untransformed
    # data)', and Table 4 footnote d 'Additive error in logarithmic
    # domain' -> Cc ~ lnorm(expSd). The tabulated 0.0420 is the NONMEM
    # $SIGMA VARIANCE on the log scale, so the nlmixr2 SD is its square
    # root: sqrt(0.0420) = 0.2049, i.e. a 20.7% residual CV. Operator-
    # ratified reading (sidecar oare_PMC11315837 request-001 Q2). The
    # variance reading is supported by (a) the NONMEM $SIGMA convention,
    # (b) this being the only variability row left in raw NONMEM units
    # while every other was converted to CV%, (c) the reported 2.4% RSE
    # matching the ~3.2% = sqrt(2/1920) expected for a variance from
    # 1920 observations rather than the ~1.6% an SD would give, and
    # (d) a 4.2% residual CV being implausibly tight for LC-MS/MS plasma
    # data.
    # ------------------------------------------------------------------
    expSd <- sqrt(0.0420); label("Exponential (log-domain additive) residual error SD (fraction)")                 # Table 4: Additive error = 0.0420 (RSE 2.4%, epsilon-shrinkage 9.5%; SIR median 0.0418, 95% CI 0.0397-0.438 -- the upper bound is a probable typo for 0.0438)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Occasion indicators multiplexing the inter-occasion-variability
    #    etas. OCC = 1 is the first crossover period, OCC = 2 the second.
    # ------------------------------------------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_ktr <- oc1 * etaiov_ktr_1 + oc2 * etaiov_ktr_2
    iov_q   <- oc1 * etaiov_q_1   + oc2 * etaiov_q_2

    # ------------------------------------------------------------------
    # 2. Individual parameters. The two CYP2C9 activity-score indicators
    #    are mutually exclusive; both 0 gives the extensive-metabolizer
    #    reference CL/F. Covariate centring values are the paper's own
    #    (77.3 kg mean total body weight, 29.0 year median age).
    # ------------------------------------------------------------------
    d1 <- exp(ld1 + etald1)
    cl <- exp(lcl + etalcl +
                e_cyp2c9im_as15_lcl   * CYP2C9_IM_AS15 +
                e_cyp2c9pmim_as10_lcl * CYP2C9_PM_IM_AS10) * (WT / 77.3)^e_wt_cl
    vc  <- exp(lvc + etalvc)
    q   <- exp(lq + iov_q) * (AGE / 29.0)^e_age_q
    vp  <- exp(lvp)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2)
    ktr <- exp(lktr + iov_ktr)

    # Fixed physiologic liver constants.
    qh     <- exp(lqh)
    vliver <- exp(lvliver)

    # Relative bioavailability: the reference suspension is the anchor at
    # F = 1 exactly, so the F1 term (and its inter-individual
    # variability) applies only on test-formulation dose records.
    frel <- exp((lf1 + etalf1) * FORM_NOSCAPINE_TEST)

    # ------------------------------------------------------------------
    # 3. ODE system (Fig. 1). The dose is released zero-order into the
    #    first of four sequential transit compartments, each transferring
    #    at the average transit rate TR. The chain empties into the LIVER,
    #    not the central compartment: apparent clearance CL/F leaves from
    #    the liver (rate constant CL/F / Vh), which is what generates the
    #    profound first-pass extraction (Sect. 3.3.3, dOFV = -562).
    #    Liver and plasma exchange at the fixed liver plasma flow Qh --
    #    liver to plasma with rate constant Qh/Vh and plasma to liver with
    #    Qh/(Vc/F), exactly as labelled on the Fig. 1 arrows.
    # ------------------------------------------------------------------
    d/dt(transit1)    <- -ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3 - ktr * transit4
    d/dt(liver)       <-  ktr * transit4 + qh * central / vc -
                          qh * liver / vliver - cl * liver / vliver
    d/dt(central)     <-  qh * liver / vliver - qh * central / vc -
                          q  * central / vc + q  * peripheral1 / vp -
                          q2 * central / vc + q2 * peripheral2 / vp2
    d/dt(peripheral1) <-  q  * central / vc - q  * peripheral1 / vp
    d/dt(peripheral2) <-  q2 * central / vc - q2 * peripheral2 / vp2

    # ------------------------------------------------------------------
    # 4. Zero-order release and formulation bioavailability. Dose records
    #    must target cmt = "transit1" with rate = -2 so rxode2 uses the
    #    modelled duration dur(transit1) = D1 (Fig. 1 'Rate ... equivalent
    #    to dose/absorption duration').
    # ------------------------------------------------------------------
    dur(transit1) <- d1
    f(transit1)   <- frel

    # ------------------------------------------------------------------
    # 5. Observation. Doses in mg and volumes in L give mg/L; x 1000
    #    converts to the ug/L of Table 3. (Fig. 3's VPC axis is
    #    ln(concentration) in mg/L, i.e. ln(Cc) - ln(1000).)
    # ------------------------------------------------------------------
    Cc <- central / vc * 1000

    Cc ~ lnorm(expSd)
  })
}
