Crass_2025_pegcetacoplan_ga_progression <- function() {
  description <- paste0(
    "Linear (zero-order) disease-progression model for geographic atrophy ",
    "(GA) lesion area secondary to age-related macular degeneration, fit to ",
    "sham-treated study eyes and untreated fellow eyes pooled across the ",
    "phase II FILLY and phase III OAKS and DERBY pegcetacoplan trials. This ",
    "is the untreated natural-history layer of the Crass 2025 model series -- ",
    "no drug effect is present. Lesion area in each eye is described ",
    "algebraically as an initial area plus a constant growth rate multiplied ",
    "by study time: area = rbase + slope * time. Study eye and fellow eye are ",
    "fit jointly with a 4x4 correlated between-subject random-effect block ",
    "over (initial area, slope) x (study eye, fellow eye); every eta is ",
    "Box-Cox transformed (Petersson 2009 form). Unilateral GA at baseline ",
    "(vs bilateral) is a structural covariate on both study-eye initial area ",
    "and study-eye slope; absence of subfoveal involvement, unifocal (vs ",
    "multifocal) lesion morphology and more than 20 intermediate or large ",
    "drusen groups act on the study-eye slope only. Residual error is ",
    "combined proportional plus additive, estimated separately for the two ",
    "eyes. There are no ODE states and no dosing events -- both outputs are ",
    "algebraic functions of time, so observation records must be tagged with ",
    "dvid rather than cmt."
  )
  reference <- paste(
    "Crass RL, Prem K, Gaudreault F, Lusk E, Ribeiro R, Chapel S, Baumal CR.",
    "Pharmacokinetic/pharmacodynamic analysis of geographic atrophy lesion",
    "area in patients receiving pegcetacoplan treatment or sham.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(2):257-267.",
    "doi:10.1002/psp4.13264.",
    "Parameter values from Supporting Information Data S1, Table S4;",
    "model equations from Data S1 'The following equations describe the",
    "exposure-response model' (drug-effect term omitted for this base model).",
    sep = " "
  )
  vignette <- "Crass_2025_pegcetacoplan_geographic_atrophy"

  # The primary parameters are the linear-scale initial lesion area and the
  # linear-scale progression slope, both carried on the log scale in ini() so
  # the exponential IIV is the canonical eta + l<param> pairing. The four etas
  # are Box-Cox transformed inside model() before entering the exponential.

  units <- list(
    time          = "day",
    dosing        = "(none; this model has no drug input -- it describes untreated GA progression)",
    concentration = "(GA lesion area, mm^2; observations lesionStudy and lesionFellow)"
  )

  covariateData <- list(
    DIS_GA_UNILATERAL = list(
      description        = "1 = unilateral geographic atrophy at baseline (GA present in the study eye only); 0 = bilateral GA (GA present in both eyes at baseline). Time-fixed per subject. In Crass 2025 this is a prespecified structural covariate carried on both the study-eye initial lesion area and the study-eye time slope.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (bilateral GA). 81% of the pooled Crass 2025 population had bilateral disease at baseline (paper Table 1).",
      notes              = "Fellow-eye observations were excluded by the authors for patients with unilateral GA at baseline who converted to bilateral disease during the study, so a subject with DIS_GA_UNILATERAL = 1 contributes study-eye records only in the source dataset. The covariate does not act on either fellow-eye parameter.",
      source_name        = "UNILATGA"
    ),
    DIS_GA_NONSUBFOVEAL = list(
      description        = "1 = the study-eye GA lesion does not involve the fovea (nonsubfoveal); 0 = subfoveal involvement. Time-fixed per subject; assessed at baseline by fundus autofluorescence.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (subfoveal involvement). 63% of study eyes had subfoveal lesions at baseline (paper Table 1).",
      notes              = "Encoded in the source paper's direction (NOFOV = 1 for the nonsubfoveal group) so the published positive coefficient applies directly. Acts on the study-eye time slope only; the corresponding effect on initial lesion area was removed in backward elimination (Table S2 step 3, dOFV 1.343).",
      source_name        = "NOFOV"
    ),
    DIS_GA_UNIFOCAL = list(
      description        = "1 = unifocal study-eye GA lesion (a single atrophic patch); 0 = multifocal lesion (two or more patches). Time-fixed per subject; assessed at baseline by fundus autofluorescence.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (multifocal). 70% of study eyes were multifocal at baseline (paper Table 1).",
      notes              = "Acts on the study-eye time slope only; the corresponding effect on initial lesion area was removed in backward elimination (Table S2 step 1, dOFV 0.584).",
      source_name        = "UNIFOC"
    ),
    DRUSEN_GT20 = list(
      description        = "1 = more than 20 intermediate or large drusen groups (drusen diameter >= 63 micrometres, AREDS simplified severity scale) present in the study eye at baseline; 0 = 20 or fewer. Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (20 or fewer intermediate or large drusen groups). 45% of study eyes had more than 20 at baseline (paper Table 1); two OAKS patients had missing data.",
      notes              = "The source paper supplies the already-dichotomised indicator rather than the raw drusen-group count; the >= 63 micrometre size threshold follows Ferris 2005 (AREDS report No. 18). Acts on the study-eye time slope only; the corresponding effect on initial lesion area was removed in backward elimination (Table S2 step 4, dOFV 3.007).",
      source_name        = "MOREDR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1501L,
    n_studies      = 3L,
    age_range      = "60-100 years (pooled); FILLY enrolled patients aged >= 50 years, OAKS and DERBY >= 60 years",
    age_median     = "79 years",
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = 62,
    race_ethnicity = c(White = 94, Black = 0.4, Asian = 0.5, Other = 0.2, Missing = 5),
    disease_state  = "Geographic atrophy secondary to age-related macular degeneration. Baseline total GA lesion area 2.5-17.5 mm^2 in the study eye with at least one lesion >= 1.25 mm^2 if multifocal; best-corrected visual acuity 24 ETDRS letters or better. 81% bilateral GA, 63% subfoveal, 70% multifocal, 82% pseudodrusen, 45% with more than 20 intermediate or large drusen groups.",
    dose_range     = "n/a for this model -- the estimation dataset for the disease-progression model comprised only sham-treated study eyes and untreated fellow eyes (498 of the 1501 patients received sham; the remaining patients contributed fellow-eye records only). The parent trials randomised patients 2:2:1:1 to intravitreal pegcetacoplan 15 mg monthly, pegcetacoplan 15 mg every other month, sham monthly, or sham every other month.",
    regions        = "(not reported in Crass 2025)",
    notes          = "Pooled from the 12-month phase II FILLY trial (NCT02503332, n = 246) and the 24-month phase III OAKS (NCT03525613, n = 636) and DERBY (NCT03525600, n = 620) trials; paper Table 1 reports n = 1502 including one patient without quantifiable PK samples who was excluded from the PK/PD modelling analyses, leaving 1501. Median baseline lesion area 8 mm^2 (range 2-18) in study eyes. Follow-up: 12 months of treatment plus 6 months off-treatment in FILLY, 24 months of treatment in OAKS and DERBY. Covariates screened by the authors but NOT retained in the final model: low-luminance deficit (removed on initial lesion area, Table S2 step 2, dOFV 1.339; and on time slope, step 5, dOFV 4.663) and, in prespecified post hoc analyses on the final PK/PD model, concomitant anti-vascular endothelial growth factor medication (used by 8% of study eyes and 16% of fellow eyes) and anti-polyethylene-glycol antidrug antibodies (prevalence 67% for the PEG moiety and 5% for the peptide moiety) -- neither was predicted to have a clinically meaningful impact (paper Figures 5 and 6)."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural fixed effects -- Crass 2025 Data S1 Table S4 (disease
    # progression model). Estimates are on the linear scale in the source
    # table and are log-transformed here so the source's exponential IIV
    # (INIT = theta * exp(eta_BC)) becomes exp(l<param> + eta_BC).
    # ------------------------------------------------------------------
    lrbase_study  <- log(7.73)    ; label("Log typical initial GA lesion area, study eye (mm^2)")
    # Table S4 "Study eye initial lesion area, mm2" = 7.73 (95% CI 7.32, 8.14; ASE 0.208; %RSE 2.69)

    lslope_study  <- log(0.00554) ; label("Log typical GA lesion-area progression rate, study eye (mm^2/day)")
    # Table S4 "Study eye time slope, mm2/day" = 0.00554 (95% CI 0.00511, 0.00597; ASE 0.000221; %RSE 3.99); 0.00554 * 365.25 = 2.02 mm^2/year, the paper's headline untreated growth rate

    lrbase_fellow <- log(7.28)    ; label("Log typical initial GA lesion area, fellow eye (mm^2)")
    # Table S4 "Fellow eye initial lesion area, mm2" = 7.28 (95% CI 6.92, 7.64; ASE 0.184; %RSE 2.53)

    lslope_fellow <- log(0.00486) ; label("Log typical GA lesion-area progression rate, fellow eye (mm^2/day)")
    # Table S4 "Fellow eye time slope, mm2/day" = 0.00486 (95% CI 0.00467, 0.00505; ASE 0.0000971; %RSE 2.00)

    # ------------------------------------------------------------------
    # Covariate effects. All are proportional shifts entering as
    # (1 + effect * INDICATOR) per the Data S1 equations for INIT_i,study
    # and SLOPE_i,study, so a value of -0.151 means a 15.1% lower slope.
    # ------------------------------------------------------------------
    e_unilat_rbase_study    <- -0.169 ; label("Proportional shift in study-eye initial lesion area for unilateral vs bilateral GA (fraction)")
    # Table S4 "Unilateral GA on study eye initial lesion area" = -0.169 (95% CI -0.258, -0.079; ASE 0.0458; %RSE 27.1); Data S1 theta_11

    e_unilat_slope_study    <- -0.151 ; label("Proportional shift in study-eye progression rate for unilateral vs bilateral GA (fraction)")
    # Table S4 "Unilateral GA on study eye time slope" = -0.151 (95% CI -0.257, -0.046; ASE 0.0538; %RSE 35.5); Data S1 theta_12

    e_nonsubfov_slope_study <-  0.179 ; label("Proportional shift in study-eye progression rate for nonsubfoveal vs subfoveal GA (fraction)")
    # Table S4 "No subfoveal involvement on time slope" = 0.179 (95% CI 0.083, 0.276; ASE 0.0493; %RSE 27.5); Data S1 theta_19

    e_unifocal_slope_study  <- -0.156 ; label("Proportional shift in study-eye progression rate for unifocal vs multifocal GA (fraction)")
    # Table S4 "Unifocal GA on time slope" = -0.156 (95% CI -0.230, -0.083; ASE 0.0376; %RSE 24.0); Data S1 theta_20

    e_drusen_slope_study    <- -0.147 ; label("Proportional shift in study-eye progression rate for more than 20 vs 20 or fewer intermediate/large drusen groups (fraction)")
    # Table S4 "More (>20) intermediate/large drusen groups on time slope" = -0.147 (95% CI -0.218, -0.075; ASE 0.0365; %RSE 24.9); Data S1 theta_21

    # ------------------------------------------------------------------
    # Box-Cox shape parameters (lambda) for the four between-subject
    # random effects. Data S1: ETIS_i = (exp(eta)^lambda - 1) / lambda,
    # equivalently (exp(lambda * eta) - 1) / lambda (Petersson 2009 form).
    # ------------------------------------------------------------------
    boxcox_rbase_study  <- -0.306 ; label("Box-Cox shape parameter for the IIV on study-eye initial lesion area (unitless)")
    # Table S4 "Box-Cox study eye initial lesion area distribution shape parameter" = -0.306 (95% CI -0.482, -0.130; ASE 0.0899; %RSE 29.4)

    boxcox_slope_study  <- -0.270 ; label("Box-Cox shape parameter for the IIV on study-eye progression rate (unitless)")
    # Table S4 "Box-Cox study eye time slope distribution shape parameter" = -0.270 (95% CI -0.450, -0.090; ASE 0.0918; %RSE 34.0)

    boxcox_rbase_fellow <- -0.540 ; label("Box-Cox shape parameter for the IIV on fellow-eye initial lesion area (unitless)")
    # Table S4 "Box-Cox fellow eye initial lesion area distribution shape parameter" = -0.540 (95% CI -0.595, -0.484; ASE 0.0284; %RSE 5.26)

    boxcox_slope_fellow <- -0.348 ; label("Box-Cox shape parameter for the IIV on fellow-eye progression rate (unitless)")
    # Table S4 "Box-Cox fellow eye time slope distribution shape parameter" = -0.348 (95% CI -0.428, -0.269; ASE 0.0406; %RSE 11.6)

    # ------------------------------------------------------------------
    # 4x4 correlated between-subject variability block, in the source's
    # eta order: 1 = init study, 2 = slope study, 3 = init fellow,
    # 4 = slope fellow. Values are variances (diagonal) and covariances
    # (off-diagonal) from Table S4 "Interindividual Variability".
    # Reported correlations reproduce exactly from these numbers, e.g.
    # 0.236 / sqrt(0.285 * 0.394) = 0.705, the paper's study-vs-fellow
    # slope correlation quoted in the Discussion.
    # eta-shrinkage: 1e-10% (init study), 4.20% (slope study),
    # 0.777% (init fellow), 7.85% (slope fellow).
    # ------------------------------------------------------------------
    etalrbase_study + etalslope_study + etalrbase_fellow + etalslope_fellow ~
      c(0.256,
        0.111, 0.285,
        0.294, 0.161, 0.737,
        0.141, 0.236, 0.220, 0.394)

    # ------------------------------------------------------------------
    # Residual error: combined proportional + additive, estimated
    # separately per eye. The Table S4 "Estimate" column carries the
    # proportional term on the fraction (SD) scale -- the table footnote
    # defines the transformed column as "the product of the estimate and
    # 100% for proportional parameters", and the transformed column gives
    # the extra significant figure used here (2.44% -> 0.0244).
    # ------------------------------------------------------------------
    propSd_lesionStudy  <- 0.0244 ; label("Proportional residual error, study-eye lesion area (fraction)")
    # Table S4 "Study eye proportional, proportion" = 0.024 (95% CI 0.021, 0.028; ASE 0.00163; %RSE 6.68); transformed 2.44% (2.12, 2.76)

    addSd_lesionStudy   <- 0.192  ; label("Additive residual error, study-eye lesion area (mm^2)")
    # Table S4 "Study eye additive, mm2" = 0.192 (95% CI 0.154, 0.231; ASE 0.0198; %RSE 10.3)

    propSd_lesionFellow <- 0.0310 ; label("Proportional residual error, fellow-eye lesion area (fraction)")
    # Table S4 "Fellow eye proportional, proportion" = 0.031 (95% CI 0.028, 0.034; ASE 0.00141; %RSE 4.55); transformed 3.10% (2.82, 3.38)

    addSd_lesionFellow  <- 0.207  ; label("Additive residual error, fellow-eye lesion area (mm^2)")
    # Table S4 "Fellow eye additive, mm2" = 0.207 (95% CI 0.179, 0.234; ASE 0.0140; %RSE 6.77)
  })

  model({
    # Box-Cox transformation of each N(0, omega^2) random effect. The
    # Data S1 equations print this as (exp(eta)^lambda - 1) / lambda;
    # exp(eta)^lambda == exp(lambda * eta), and the latter form is used
    # here because it avoids raising a value to a non-integer power.
    etatr_rbase_study  <- (exp(boxcox_rbase_study  * etalrbase_study)  - 1) / boxcox_rbase_study
    etatr_slope_study  <- (exp(boxcox_slope_study  * etalslope_study)  - 1) / boxcox_slope_study
    etatr_rbase_fellow <- (exp(boxcox_rbase_fellow * etalrbase_fellow) - 1) / boxcox_rbase_fellow
    etatr_slope_fellow <- (exp(boxcox_slope_fellow * etalslope_fellow) - 1) / boxcox_slope_fellow

    # Individual study-eye parameters. Data S1:
    #   INIT_i,study  = theta_init  * (1 + UNILATGA * theta_11) * exp(ETIS_i)
    #   SLOPE_i,study = theta_slope * (1 + UNILATGA * theta_12)
    #                                * (1 + NOFOV    * theta_19)
    #                                * (1 + UNIFOC   * theta_20)
    #                                * (1 + MOREDR   * theta_21) * exp(ETSS_i)
    rbaseStudy <- exp(lrbase_study + etatr_rbase_study) *
      (1 + e_unilat_rbase_study * DIS_GA_UNILATERAL)

    slopeStudy <- exp(lslope_study + etatr_slope_study) *
      (1 + e_unilat_slope_study    * DIS_GA_UNILATERAL) *
      (1 + e_nonsubfov_slope_study * DIS_GA_NONSUBFOVEAL) *
      (1 + e_unifocal_slope_study  * DIS_GA_UNIFOCAL) *
      (1 + e_drusen_slope_study    * DRUSEN_GT20)

    # Individual fellow-eye parameters. The fellow eye is untreated and
    # carries no covariate effects. Data S1:
    #   INIT_i,fellow  = theta_init  * exp(ETIF_i)
    #   SLOPE_i,fellow = theta_slope * exp(ETSF_i)
    # (the Data S1 equations reuse the symbols theta_init / theta_slope
    # for both eyes, but Table S4 reports separate study-eye and
    # fellow-eye fixed effects, which is what is encoded here).
    rbaseFellow <- exp(lrbase_fellow + etatr_rbase_fellow)
    slopeFellow <- exp(lslope_fellow + etatr_slope_fellow)

    # Linear (zero-order) disease progression. There are no ODE states --
    # both outputs are algebraic in time, so simulation event tables must
    # tag observation rows with dvid and leave cmt missing.
    lesionStudy  <- rbaseStudy  + slopeStudy  * time
    lesionFellow <- rbaseFellow + slopeFellow * time

    lesionStudy  ~ add(addSd_lesionStudy)  + prop(propSd_lesionStudy)
    lesionFellow ~ add(addSd_lesionFellow) + prop(propSd_lesionFellow)
  })
}
