Crass_2025_pegcetacoplan_ga_exposureresponse <- function() {
  description <- paste0(
    "Final PK/PD (exposure-response) layer of the Crass 2025 geographic ",
    "atrophy (GA) lesion-area model series: the linear (zero-order) ",
    "disease-progression model with a linear-in-log-concentration ",
    "pegcetacoplan effect on the study-eye lesion growth rate, driven by the ",
    "instantaneous individual-predicted vitreous humour pegcetacoplan ",
    "concentration. Fit to pooled study-eye and fellow-eye data from the ",
    "phase II FILLY and phase III OAKS and DERBY trials. Lesion area in each ",
    "eye is algebraic in study time: area = rbase + dslope(t) * time, where ",
    "dslope(t) = slope * (1 + effect * log(Cvitreous(t) + 1)). Study eye and ",
    "fellow eye are fit jointly with a 4x4 correlated between-subject ",
    "random-effect block over (initial area, slope) x (study eye, fellow eye); ",
    "every eta is Box-Cox transformed (Petersson 2009 form). Baseline disease ",
    "covariates act on the study eye only. Residual error is combined ",
    "proportional plus additive per eye. There are no ODE states and no ",
    "dosing events: the vitreous concentration is supplied as the CEFFECT ",
    "time-varying covariate column (the source analysis took it from an ",
    "external pegcetacoplan population PK model, so the PK layer is not part ",
    "of this file), and observation records must be tagged with dvid rather ",
    "than cmt."
  )
  reference <- paste(
    "Crass RL, Prem K, Gaudreault F, Lusk E, Ribeiro R, Chapel S, Baumal CR.",
    "Pharmacokinetic/pharmacodynamic analysis of geographic atrophy lesion",
    "area in patients receiving pegcetacoplan treatment or sham.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(2):257-267.",
    "doi:10.1002/psp4.13264.",
    "Parameter values from Supporting Information Data S1, Table S6;",
    "model equations from Data S1 'The following equations describe the",
    "exposure-response model'. Vitreous concentrations were generated in the",
    "source analysis by an external population PK model reported only as a",
    "conference presentation (Crass RL, Smith B, Epling D, Chapel S, Prem K,",
    "Gaudreault F. Population pharmacokinetics of pegcetacoplan in patients",
    "with neovascular age-related macular degeneration and geographic atrophy.",
    "American Conference on Pharmacometrics, November 5-8 2023, National",
    "Harbor MD), which is not reproduced here; supply the concentration",
    "trajectory through the CEFFECT covariate column.",
    sep = " "
  )
  vignette <- "Crass_2025_pegcetacoplan_geographic_atrophy"

  units <- list(
    time          = "day",
    dosing        = "(none; pegcetacoplan exposure enters as the CEFFECT vitreous-concentration covariate, not as a dose record)",
    concentration = "(GA lesion area, mm^2; observations lesionStudy and lesionFellow. The CEFFECT drug-concentration covariate is in ug/mL.)"
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
      notes              = "Encoded in the source paper's direction (NOFOV = 1 for the nonsubfoveal group) so the published positive coefficient applies directly. Acts on the study-eye time slope only. The paper's Figure 2 forest-plot fold-effects are 1 + theta from this table, e.g. 1 + 0.118 = 1.12 for nonsubfoveal involvement.",
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
    CEFFECT = list(
      description        = "Instantaneous individual-predicted pegcetacoplan concentration in the vitreous humour of the study eye, supplied per event record as a time-varying covariate and consumed directly as the pharmacodynamic driver. Set to 0 for sham-treated study eyes, for untreated fellow eyes, and before the first pegcetacoplan injection.",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = "n/a -- enters as log(CEFFECT + 1); CEFFECT = 0 gives log(1) = 0 and recovers the untreated progression rate. Reference values: an intravitreal 15 mg dose into the assumed 4 mL vitreous volume gives an immediate post-dose concentration of about 3750 ug/mL, for which 1 + (-0.026) * log(3751) = 0.786 -- consistent with the 0.80-fold monthly effect estimated by the companion dose-response model.",
      notes              = "The source analysis generated this trajectory from an external pegcetacoplan population PK model (conference presentation, not on disk) using each patient's actual dosing history and the empirical Bayes estimate of the vitreous-to-serum absorption rate constant, with the vitreous volume assumed to be 4 mL (paper Methods, 'Pegcetacoplan exposure' and 'PK/PD model'). That PK model is NOT reproduced in this file; the user must supply CEFFECT. Lesion area was assessed at dosing visits, so the source observations pair each lesion measurement with a near-trough vitreous concentration. Member of the general 'effect-site PD driver' CEFFECT family; here the biophase is the vitreous humour of an intravitreally dosed eye rather than an IV-anaesthetic biophase.",
      source_name        = "Cv(t)"
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
    notes          = "Pooled from the 12-month phase II FILLY trial (NCT02503332, n = 246) and the 24-month phase III OAKS (NCT03525613, n = 636) and DERBY (NCT03525600, n = 620) trials; paper Table 1 reports n = 1502 including one patient without quantifiable PK samples who was excluded from the PK/PD modelling analyses. A linear-in-log-concentration form was selected over a nonlinear Emax form because no Emax parameterisation could be precisely estimated (Data S1 Table S3 run 3: EC50 95% CI 0.000878 to 261000). Prespecified post hoc analyses on this model found no clinically meaningful effect of concomitant anti-vascular endothelial growth factor medication or of anti-polyethylene-glycol antidrug antibodies (paper Figures 5 and 6)."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural fixed effects -- Crass 2025 Data S1 Table S6 (PK/PD
    # model). Log-transformed here so the source's exponential IIV
    # becomes exp(l<param> + eta_BC).
    # ------------------------------------------------------------------
    lrbase_study  <- log(7.69)    ; label("Log typical initial GA lesion area, study eye (mm^2)")
    # Table S6 "Study eye initial lesion area, mm2" = 7.69 (95% CI 7.42, 7.97; ASE 0.141; %RSE 1.83)

    lslope_study  <- log(0.00559) ; label("Log typical GA lesion-area progression rate, study eye (mm^2/day)")
    # Table S6 "Study eye time slope, mm2/day" = 0.00559 (95% CI 0.00527, 0.00591; ASE 0.000163; %RSE 2.92)

    lrbase_fellow <- log(7.25)    ; label("Log typical initial GA lesion area, fellow eye (mm^2)")
    # Table S6 "Fellow eye initial lesion area, mm2" = 7.25 (95% CI 6.90, 7.60; ASE 0.179; %RSE 2.47)

    lslope_fellow <- log(0.00490) ; label("Log typical GA lesion-area progression rate, fellow eye (mm^2/day)")
    # Table S6 "Fellow eye time slope, mm2/day" = 0.00490 (95% CI 0.00471, 0.00509; ASE 0.0000967; %RSE 1.97)

    # ------------------------------------------------------------------
    # Baseline-disease covariate effects; proportional shifts entering as
    # (1 + effect * INDICATOR) per the Data S1 equations. These are the
    # values behind the paper's Figure 2 forest plot (fold-effect =
    # 1 + theta): 1.12 nonsubfoveal, 0.836 unifocal, 0.870 more drusen,
    # 0.861 unilateral on slope.
    # ------------------------------------------------------------------
    e_unilat_rbase_study    <- -0.123 ; label("Proportional shift in study-eye initial lesion area for unilateral vs bilateral GA (fraction)")
    # Table S6 "Unilateral GA on study eye initial lesion area" = -0.123 (95% CI -0.183, -0.064; ASE 0.0303; %RSE 24.5); Data S1 theta_11

    e_unilat_slope_study    <- -0.139 ; label("Proportional shift in study-eye progression rate for unilateral vs bilateral GA (fraction)")
    # Table S6 "Unilateral GA on study eye time slope" = -0.139 (95% CI -0.204, -0.074; ASE 0.0331; %RSE 23.8); Data S1 theta_12; paper Results 0.860-fold (90% CI 0.808, 0.913)

    e_nonsubfov_slope_study <-  0.118 ; label("Proportional shift in study-eye progression rate for nonsubfoveal vs subfoveal GA (fraction)")
    # Table S6 "No subfoveal involvement on time slope" = 0.118 (95% CI 0.066, 0.169; ASE 0.0262; %RSE 22.3); Data S1 theta_19; paper Results 1.12-fold (90% CI 1.07, 1.16)

    e_unifocal_slope_study  <- -0.164 ; label("Proportional shift in study-eye progression rate for unifocal vs multifocal GA (fraction)")
    # Table S6 "Unifocal GA on time slope" = -0.164 (95% CI -0.205, -0.124; ASE 0.0208; %RSE 12.6); Data S1 theta_20; paper Results 0.837-fold (90% CI 0.800, 0.872)

    e_drusen_slope_study    <- -0.130 ; label("Proportional shift in study-eye progression rate for more than 20 vs 20 or fewer intermediate/large drusen groups (fraction)")
    # Table S6 "More (>20) intermediate/large drusen groups on time slope" = -0.130 (95% CI -0.171, -0.090; ASE 0.0206; %RSE 15.8); Data S1 theta_21; paper Results 0.869-fold (90% CI 0.837, 0.904)

    # ------------------------------------------------------------------
    # Pegcetacoplan exposure-response slope. Data S1 theta_23: the
    # proportional shift in the study-eye progression rate per unit of
    # natural-log-transformed individual-predicted vitreous
    # concentration at time t.
    # ------------------------------------------------------------------
    e_ceffect_slope_study <- -0.026 ; label("Proportional shift in study-eye progression rate per unit log vitreous pegcetacoplan concentration (fraction per log(ug/mL))")
    # Table S6 "Drug effect slope, proportion/log(ug/mL)" = -0.026 (95% CI -0.032, -0.020; ASE 0.00304; %RSE 11.7); paper Results: "the rate of GA lesion area growth was estimated to decrease by 2.6% per unit of log-transformed vitreous pegcetacoplan concentration"

    # ------------------------------------------------------------------
    # Box-Cox shape parameters (lambda) for the four random effects.
    # Data S1: ETIS_i = (exp(eta)^lambda - 1) / lambda.
    # ------------------------------------------------------------------
    boxcox_rbase_study  <- -0.313 ; label("Box-Cox shape parameter for the IIV on study-eye initial lesion area (unitless)")
    # Table S6 "Box-Cox study eye initial lesion area distribution shape parameter" = -0.313 (95% CI -0.428, -0.198; ASE 0.0584; %RSE 18.7)

    boxcox_slope_study  <- -0.287 ; label("Box-Cox shape parameter for the IIV on study-eye progression rate (unitless)")
    # Table S6 "Box-Cox study eye time slope distribution shape parameter" = -0.287 (95% CI -0.377, -0.198; ASE 0.0455; %RSE 15.8)

    boxcox_rbase_fellow <- -0.537 ; label("Box-Cox shape parameter for the IIV on fellow-eye initial lesion area (unitless)")
    # Table S6 "Box-Cox fellow eye initial lesion area distribution shape parameter" = -0.537 (95% CI -0.592, -0.481; ASE 0.0282; %RSE 5.26)

    boxcox_slope_fellow <- -0.405 ; label("Box-Cox shape parameter for the IIV on fellow-eye progression rate (unitless)")
    # Table S6 "Box-Cox fellow eye time slope distribution shape parameter" = -0.405 (95% CI -0.478, -0.332; ASE 0.0373; %RSE 9.21)

    # ------------------------------------------------------------------
    # 4x4 correlated between-subject variability block in the source's
    # eta order: 1 = init study, 2 = slope study, 3 = init fellow,
    # 4 = slope fellow. Variances on the diagonal, covariances off it,
    # from Table S6 "Interindividual Variability". The reported
    # correlations reproduce from these numbers, e.g.
    # 0.256 / sqrt(0.307 * 0.399) = 0.731.
    # eta-shrinkage: 0.0803% (init study), 6.35% (slope study),
    # 1.27% (init fellow), 7.54% (slope fellow).
    # ------------------------------------------------------------------
    etalrbase_study + etalslope_study + etalrbase_fellow + etalslope_fellow ~
      c(0.259,
        0.101, 0.307,
        0.295, 0.167, 0.744,
        0.127, 0.256, 0.220, 0.399)

    # ------------------------------------------------------------------
    # Residual error: combined proportional + additive per eye. The
    # proportional term is on the fraction (SD) scale; the extra
    # significant figure comes from the table's transformed column
    # (2.76% -> 0.0276).
    # ------------------------------------------------------------------
    propSd_lesionStudy  <- 0.0276 ; label("Proportional residual error, study-eye lesion area (fraction)")
    # Table S6 "Study eye proportional, proportion" = 0.028 (95% CI 0.026, 0.030; ASE 0.000960; %RSE 3.47); transformed 2.76% (2.58, 2.95)

    addSd_lesionStudy   <- 0.175  ; label("Additive residual error, study-eye lesion area (mm^2)")
    # Table S6 "Study eye additive, mm2" = 0.175 (95% CI 0.153, 0.197; ASE 0.0112; %RSE 6.43)

    propSd_lesionFellow <- 0.0309 ; label("Proportional residual error, fellow-eye lesion area (fraction)")
    # Table S6 "Fellow eye proportional, proportion" = 0.031 (95% CI 0.028, 0.034; ASE 0.00141; %RSE 4.54); transformed 3.09% (2.82, 3.37)

    addSd_lesionFellow  <- 0.207  ; label("Additive residual error, fellow-eye lesion area (mm^2)")
    # Table S6 "Fellow eye additive, mm2" = 0.207 (95% CI 0.180, 0.235; ASE 0.0140; %RSE 6.75)
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

    # Exposure-driven study-eye slope. Data S1:
    #   DSLOPE_i,study(t) = SLOPE_i,study * (1 + theta_23 * log(Cv(t) + 1))
    # The +1 inside the logarithm keeps the untreated state well defined:
    # CEFFECT = 0 gives log(1) = 0 and DSLOPE collapses to SLOPE.
    dslopeStudy <- slopeStudy * (1 + e_ceffect_slope_study * log(CEFFECT + 1))

    # Fellow eye: untreated, no covariate or drug effect (Data S1
    # INIT_i,fellow, SLOPE_i,fellow). The Data S1 equations reuse the
    # symbols theta_init / theta_slope for both eyes, but Table S6
    # reports separate study-eye and fellow-eye fixed effects, which is
    # what is encoded here.
    rbaseFellow <- exp(lrbase_fellow + etatr_rbase_fellow)
    slopeFellow <- exp(lslope_fellow + etatr_slope_fellow)

    # Linear (zero-order) disease progression with a time-varying
    # treated slope. Note that the source model is algebraic, NOT an
    # integrated ODE: Data S1 writes A_i,study(t) = INIT + DSLOPE(t) * t,
    # i.e. the slope evaluated at the current time multiplied by the
    # current time, not the integral of the slope. That form is
    # reproduced verbatim here. There are no ODE states, so simulation
    # event tables must tag observation rows with dvid and leave cmt
    # missing.
    lesionStudy  <- rbaseStudy  + dslopeStudy * time
    lesionFellow <- rbaseFellow + slopeFellow * time

    lesionStudy  ~ add(addSd_lesionStudy)  + prop(propSd_lesionStudy)
    lesionFellow ~ add(addSd_lesionFellow) + prop(propSd_lesionFellow)
  })
}
