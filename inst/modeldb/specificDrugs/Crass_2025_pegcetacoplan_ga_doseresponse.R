Crass_2025_pegcetacoplan_ga_doseresponse <- function() {
  description <- paste0(
    "Dose-response layer of the Crass 2025 geographic atrophy (GA) lesion-area ",
    "model series: the linear (zero-order) disease-progression model with a ",
    "step-function pegcetacoplan treatment effect on the study-eye lesion ",
    "growth rate, estimated separately for each active regimen (intravitreal ",
    "15 mg monthly and 15 mg every other month) against sham. Fit to pooled ",
    "study-eye and fellow-eye data from the phase II FILLY and phase III OAKS ",
    "and DERBY trials. Lesion area in each eye is algebraic in study time: ",
    "area = rbase + slope * time. Study eye and fellow eye are fit jointly ",
    "with a 4x4 correlated between-subject random-effect block over (initial ",
    "area, slope) x (study eye, fellow eye); every eta is Box-Cox transformed ",
    "(Petersson 2009 form). Baseline disease covariates act on the study eye ",
    "only: unilateral GA on both initial area and slope, and nonsubfoveal ",
    "involvement, unifocal morphology and more than 20 intermediate or large ",
    "drusen groups on the slope. The fellow eye is untreated and carries no ",
    "covariate or drug effect. Residual error is combined proportional plus ",
    "additive, estimated separately per eye. There are no ODE states and no ",
    "dosing events -- the regimen enters through binary covariate columns, so ",
    "observation records must be tagged with dvid rather than cmt."
  )
  reference <- paste(
    "Crass RL, Prem K, Gaudreault F, Lusk E, Ribeiro R, Chapel S, Baumal CR.",
    "Pharmacokinetic/pharmacodynamic analysis of geographic atrophy lesion",
    "area in patients receiving pegcetacoplan treatment or sham.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(2):257-267.",
    "doi:10.1002/psp4.13264.",
    "Parameter values from Supporting Information Data S1, Table S5;",
    "model equations from Data S1 'The following equations describe the",
    "exposure-response model', with the exposure-driven DSLOPE term replaced",
    "by the per-regimen step function of Data S1 Table S3 run 1.",
    sep = " "
  )
  vignette <- "Crass_2025_pegcetacoplan_geographic_atrophy"

  units <- list(
    time          = "day",
    dosing        = "(none; the pegcetacoplan regimen enters as the REGI_QM / REGI_Q2M covariate pair, not as a dose record)",
    concentration = "(GA lesion area, mm^2; observations lesionStudy and lesionFellow)"
  )

  covariateData <- list(
    DIS_GA_UNILATERAL = list(
      description        = "1 = unilateral geographic atrophy at baseline (GA present in the study eye only); 0 = bilateral GA. Time-fixed per subject. Prespecified structural covariate on both the study-eye initial lesion area and the study-eye time slope.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (bilateral GA). 81% of the pooled population had bilateral disease at baseline (paper Table 1).",
      notes              = "Fellow-eye observations were excluded by the authors for patients with unilateral GA at baseline who converted to bilateral disease during the study. The covariate does not act on either fellow-eye parameter.",
      source_name        = "UNILATGA"
    ),
    DIS_GA_NONSUBFOVEAL = list(
      description        = "1 = the study-eye GA lesion does not involve the fovea (nonsubfoveal); 0 = subfoveal involvement. Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (subfoveal involvement). 63% of study eyes had subfoveal lesions at baseline (paper Table 1).",
      notes              = "Encoded in the source paper's direction (NOFOV = 1 for the nonsubfoveal group) so the published positive coefficient applies directly. Acts on the study-eye time slope only.",
      source_name        = "NOFOV"
    ),
    DIS_GA_UNIFOCAL = list(
      description        = "1 = unifocal study-eye GA lesion; 0 = multifocal lesion. Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (multifocal). 70% of study eyes were multifocal at baseline (paper Table 1).",
      notes              = "Acts on the study-eye time slope only.",
      source_name        = "UNIFOC"
    ),
    DRUSEN_GT20 = list(
      description        = "1 = more than 20 intermediate or large drusen groups (diameter >= 63 micrometres, AREDS simplified severity scale) in the study eye at baseline; 0 = 20 or fewer.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (20 or fewer intermediate or large drusen groups). 45% of study eyes had more than 20 at baseline (paper Table 1).",
      notes              = "The source paper supplies the already-dichotomised indicator rather than the raw drusen-group count. Acts on the study-eye time slope only.",
      source_name        = "MOREDR"
    ),
    REGI_QM = list(
      description        = "1 = study eye assigned to intravitreal pegcetacoplan 15 mg once monthly; 0 = any other assignment (pegcetacoplan every other month, or sham). Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on the monthly pegcetacoplan regimen). The model's reference state is REGI_QM = REGI_Q2M = 0, i.e. sham treatment. 505 of 1501 patients received pegcetacoplan monthly.",
      notes              = "REGI_QM and REGI_Q2M are mutually exclusive: at most one may be 1 for a given subject, and both are 0 for sham-treated study eyes and for all fellow eyes. Set both to 0 to recover the untreated disease-progression trajectory.",
      source_name        = "(per-arm treatment indicator; the source paper does not name the NONMEM column)"
    ),
    REGI_Q2M = list(
      description        = "1 = study eye assigned to intravitreal pegcetacoplan 15 mg every other month (EOM); 0 = any other assignment (pegcetacoplan monthly, or sham). Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on the every-other-month pegcetacoplan regimen). 498 of 1501 patients received pegcetacoplan EOM.",
      notes              = "Sibling of REGI_QM under the REGI_* regimen-indicator family. Mutually exclusive with REGI_QM; both 0 gives the sham reference.",
      source_name        = "(per-arm treatment indicator; the source paper does not name the NONMEM column)"
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
    disease_state  = "Geographic atrophy secondary to age-related macular degeneration. Baseline total GA lesion area 2.5-17.5 mm^2 in the study eye with at least one lesion >= 1.25 mm^2 if multifocal; best-corrected visual acuity 24 ETDRS letters or better. 81% bilateral GA, 63% subfoveal, 70% multifocal, 45% with more than 20 intermediate or large drusen groups.",
    dose_range     = "Intravitreal pegcetacoplan 15 mg (0.1 mL of 150 mg/mL) monthly (n = 505) or every other month (n = 498), or sham monthly / every other month (n = 498), for up to 24 months.",
    regions        = "(not reported in Crass 2025)",
    notes          = "Pooled from the 12-month phase II FILLY trial (NCT02503332, n = 246) and the 24-month phase III OAKS (NCT03525613, n = 636) and DERBY (NCT03525600, n = 620) trials; paper Table 1 reports n = 1502 including one patient without quantifiable PK samples who was excluded from the PK/PD modelling analyses. The dose-response model was selected over an Emax alternative (Data S1 Table S3 run 3) because EC50 could not be precisely estimated (95% CI 0.000878 to 261000), and over a linear-in-concentration alternative (run 6, dOFV -68.9 vs -81.6 for the dose-response step function)."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural fixed effects -- Crass 2025 Data S1 Table S5
    # (dose-response model). Log-transformed here so the source's
    # exponential IIV becomes exp(l<param> + eta_BC).
    # ------------------------------------------------------------------
    lrbase_study  <- log(7.69)    ; label("Log typical initial GA lesion area, study eye (mm^2)")
    # Table S5 "Study eye initial lesion area, mm2" = 7.69 (95% CI 7.41, 7.97; ASE 0.141; %RSE 1.83)

    lslope_study  <- log(0.00565) ; label("Log typical GA lesion-area progression rate, study eye (mm^2/day)")
    # Table S5 "Study eye time slope, mm2/day" = 0.00565 (95% CI 0.00533, 0.00597; ASE 0.000163; %RSE 2.89)

    lrbase_fellow <- log(7.25)    ; label("Log typical initial GA lesion area, fellow eye (mm^2)")
    # Table S5 "Fellow eye initial lesion area, mm2" = 7.25 (95% CI 6.90, 7.60; ASE 0.179; %RSE 2.47)

    lslope_fellow <- log(0.00491) ; label("Log typical GA lesion-area progression rate, fellow eye (mm^2/day)")
    # Table S5 "Fellow eye time slope, mm2/day" = 0.00491 (95% CI 0.00472, 0.00510; ASE 0.0000968; %RSE 1.97)

    # ------------------------------------------------------------------
    # Baseline-disease covariate effects; proportional shifts entering as
    # (1 + effect * INDICATOR) per the Data S1 equations.
    # ------------------------------------------------------------------
    e_unilat_rbase_study    <- -0.124 ; label("Proportional shift in study-eye initial lesion area for unilateral vs bilateral GA (fraction)")
    # Table S5 "Unilateral GA on study eye initial lesion area" = -0.124 (95% CI -0.183, -0.064; ASE 0.0303; %RSE 24.5); Data S1 theta_11

    e_unilat_slope_study    <- -0.138 ; label("Proportional shift in study-eye progression rate for unilateral vs bilateral GA (fraction)")
    # Table S5 "Unilateral GA on study eye time slope" = -0.138 (95% CI -0.203, -0.074; ASE 0.0331; %RSE 23.9); Data S1 theta_12

    e_nonsubfov_slope_study <-  0.119 ; label("Proportional shift in study-eye progression rate for nonsubfoveal vs subfoveal GA (fraction)")
    # Table S5 "No subfoveal involvement on time slope" = 0.119 (95% CI 0.067, 0.170; ASE 0.0263; %RSE 22.2); Data S1 theta_19

    e_unifocal_slope_study  <- -0.165 ; label("Proportional shift in study-eye progression rate for unifocal vs multifocal GA (fraction)")
    # Table S5 "Unifocal GA on time slope" = -0.165 (95% CI -0.206, -0.124; ASE 0.0209; %RSE 12.6); Data S1 theta_20

    e_drusen_slope_study    <- -0.132 ; label("Proportional shift in study-eye progression rate for more than 20 vs 20 or fewer intermediate/large drusen groups (fraction)")
    # Table S5 "More (>20) intermediate/large drusen groups on time slope" = -0.132 (95% CI -0.172, -0.091; ASE 0.0206; %RSE 15.7); Data S1 theta_21

    # ------------------------------------------------------------------
    # Pegcetacoplan step-function drug effect on the study-eye slope.
    # Both are proportional shifts, so -0.204 gives a 1 - 0.204 = 0.796
    # ("0.80-fold") lower rate of lesion growth vs sham, the paper's
    # headline monthly dose-response result.
    # ------------------------------------------------------------------
    e_regi_qm_slope_study  <- -0.204 ; label("Proportional shift in study-eye progression rate on pegcetacoplan 15 mg intravitreal monthly vs sham (fraction)")
    # Table S5 "Monthly drug effect, proportion" = -0.204 (95% CI -0.247, -0.160; ASE 0.0223; %RSE 10.9); paper Results: 0.80-fold lower rate (95% CI 0.75, 0.84) and a 20.4% reduction in lesion area change from baseline

    e_regi_q2m_slope_study <- -0.172 ; label("Proportional shift in study-eye progression rate on pegcetacoplan 15 mg intravitreal every other month vs sham (fraction)")
    # Table S5 "Every-other-month drug effect, proportion" = -0.172 (95% CI -0.219, -0.126; ASE 0.0238; %RSE 13.8); paper Results: 0.83-fold lower rate (95% CI 0.78, 0.87) and a 17.2% reduction in lesion area change from baseline

    # ------------------------------------------------------------------
    # Box-Cox shape parameters (lambda) for the four random effects.
    # Data S1: ETIS_i = (exp(eta)^lambda - 1) / lambda.
    # ------------------------------------------------------------------
    boxcox_rbase_study  <- -0.312 ; label("Box-Cox shape parameter for the IIV on study-eye initial lesion area (unitless)")
    # Table S5 "Box-Cox study eye initial lesion area distribution shape parameter" = -0.312 (95% CI -0.427, -0.198; ASE 0.0585; %RSE 18.7)

    boxcox_slope_study  <- -0.281 ; label("Box-Cox shape parameter for the IIV on study-eye progression rate (unitless)")
    # Table S5 "Box-Cox study eye time slope distribution shape parameter" = -0.281 (95% CI -0.371, -0.191; ASE 0.0458; %RSE 16.3)

    boxcox_rbase_fellow <- -0.536 ; label("Box-Cox shape parameter for the IIV on fellow-eye initial lesion area (unitless)")
    # Table S5 "Box-Cox fellow eye initial lesion area distribution shape parameter" = -0.536 (95% CI -0.592, -0.481; ASE 0.0282; %RSE 5.26)

    boxcox_slope_fellow <- -0.404 ; label("Box-Cox shape parameter for the IIV on fellow-eye progression rate (unitless)")
    # Table S5 "Box-Cox fellow eye time slope distribution shape parameter" = -0.404 (95% CI -0.478, -0.331; ASE 0.0374; %RSE 9.24)

    # ------------------------------------------------------------------
    # 4x4 correlated between-subject variability block in the source's
    # eta order: 1 = init study, 2 = slope study, 3 = init fellow,
    # 4 = slope fellow. Variances on the diagonal, covariances off it,
    # from Table S5 "Interindividual Variability". The reported
    # correlations reproduce from these numbers, e.g.
    # 0.255 / sqrt(0.306 * 0.399) = 0.729.
    # eta-shrinkage: 0.0815% (init study), 6.36% (slope study),
    # 1.27% (init fellow), 7.51% (slope fellow).
    # ------------------------------------------------------------------
    etalrbase_study + etalslope_study + etalrbase_fellow + etalslope_fellow ~
      c(0.259,
        0.101, 0.306,
        0.295, 0.166, 0.744,
        0.127, 0.255, 0.220, 0.399)

    # ------------------------------------------------------------------
    # Residual error: combined proportional + additive per eye. The
    # proportional term is on the fraction (SD) scale; the extra
    # significant figure comes from the table's transformed column
    # (2.76% -> 0.0276).
    # ------------------------------------------------------------------
    propSd_lesionStudy  <- 0.0276 ; label("Proportional residual error, study-eye lesion area (fraction)")
    # Table S5 "Study eye proportional, proportion" = 0.028 (95% CI 0.026, 0.030; ASE 0.000960; %RSE 3.48); transformed 2.76% (2.57, 2.95)

    addSd_lesionStudy   <- 0.175  ; label("Additive residual error, study-eye lesion area (mm^2)")
    # Table S5 "Study eye additive, mm2" = 0.175 (95% CI 0.153, 0.197; ASE 0.0112; %RSE 6.42)

    propSd_lesionFellow <- 0.0309 ; label("Proportional residual error, fellow-eye lesion area (fraction)")
    # Table S5 "Fellow eye proportional, proportion" = 0.031 (95% CI 0.028, 0.034; ASE 0.00141; %RSE 4.55); transformed 3.09% (2.82, 3.37)

    addSd_lesionFellow  <- 0.207  ; label("Additive residual error, fellow-eye lesion area (mm^2)")
    # Table S5 "Fellow eye additive, mm2" = 0.207 (95% CI 0.180, 0.235; ASE 0.0140; %RSE 6.76)
  })

  model({
    # Box-Cox transformation of each N(0, omega^2) random effect. Data S1
    # prints (exp(eta)^lambda - 1) / lambda; exp(eta)^lambda ==
    # exp(lambda * eta), which is the form used here.
    etatr_rbase_study  <- (exp(boxcox_rbase_study  * etalrbase_study)  - 1) / boxcox_rbase_study
    etatr_slope_study  <- (exp(boxcox_slope_study  * etalslope_study)  - 1) / boxcox_slope_study
    etatr_rbase_fellow <- (exp(boxcox_rbase_fellow * etalrbase_fellow) - 1) / boxcox_rbase_fellow
    etatr_slope_fellow <- (exp(boxcox_slope_fellow * etalslope_fellow) - 1) / boxcox_slope_fellow

    # Individual study-eye parameters (Data S1 INIT_i,study, SLOPE_i,study).
    rbaseStudy <- exp(lrbase_study + etatr_rbase_study) *
      (1 + e_unilat_rbase_study * DIS_GA_UNILATERAL)

    slopeStudy <- exp(lslope_study + etatr_slope_study) *
      (1 + e_unilat_slope_study    * DIS_GA_UNILATERAL) *
      (1 + e_nonsubfov_slope_study * DIS_GA_NONSUBFOVEAL) *
      (1 + e_unifocal_slope_study  * DIS_GA_UNIFOCAL) *
      (1 + e_drusen_slope_study    * DRUSEN_GT20)

    # Treated study-eye slope (Data S1 DSLOPE_i,study) with the drug
    # effect entered as a per-regimen step function. REGI_QM and
    # REGI_Q2M are mutually exclusive, so exactly one term is active for
    # a treated eye and neither is active for sham.
    dslopeStudy <- slopeStudy *
      (1 + e_regi_qm_slope_study  * REGI_QM +
           e_regi_q2m_slope_study * REGI_Q2M)

    # Fellow eye: untreated, no covariate effects (Data S1
    # INIT_i,fellow, SLOPE_i,fellow). The Data S1 equations reuse the
    # symbols theta_init / theta_slope for both eyes, but Table S5
    # reports separate study-eye and fellow-eye fixed effects, which is
    # what is encoded here.
    rbaseFellow <- exp(lrbase_fellow + etatr_rbase_fellow)
    slopeFellow <- exp(lslope_fellow + etatr_slope_fellow)

    # Linear (zero-order) disease progression; no ODE states, so
    # simulation event tables must tag observation rows with dvid and
    # leave cmt missing.
    lesionStudy  <- rbaseStudy  + dslopeStudy * time
    lesionFellow <- rbaseFellow + slopeFellow * time

    lesionStudy  ~ add(addSd_lesionStudy)  + prop(propSd_lesionStudy)
    lesionFellow ~ add(addSd_lesionFellow) + prop(propSd_lesionFellow)
  })
}
