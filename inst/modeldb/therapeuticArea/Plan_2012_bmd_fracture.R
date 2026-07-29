Plan_2012_bmd_fracture <- function() {
  description <- paste0(
    "Bayesian joint disease-progression model linking postmenopausal ",
    "bone mineral density (BMD) decline to a repeated time-to-fracture ",
    "hazard, fit by Plan et al. (PAGE 2012 conference poster) to the ",
    "2005-2008 NHANES dataset (1605 postmenopausal women, 204 fracture ",
    "events total; one femoral-neck BMD measure per subject and 0-5 ",
    "fracture events each). The BMD time-course is piecewise linear in ",
    "years since final menstrual period (FMP), with three annualized ",
    "percent-change slopes: menopause transition (t in [-1, 2] yr), ",
    "early postmenopausal (t in [2, 5] yr), and long-term post-FMP (t > ",
    "5 yr). The fracture hazard is exponential (Weibull shape alpha = 1 ",
    "in the final model) with a BMD-centred log-linear modulation ",
    "h(t) = exp(theta_h * (1 + theta_bmd * (BMD(t) - 0.8))). The ",
    "cumulative-hazard ODE d/dt(cumhaz) <- h accumulates from t = 0 ",
    "(FMP), and the survival probability sur = exp(-cumhaz) is exposed ",
    "as a derived output. BMI, NHANES ethnicity, and age at FMP were ",
    "carried as 'centered covariate effects on all parameters' in the ",
    "source Bayesian fit but the per-parameter coefficient magnitudes ",
    "are not reported in the poster; they are recorded in ",
    "covariatesDataExcluded for provenance. See validation vignette ",
    "Errata for the slope-unit interpretation: the printed source ",
    "equation BMD(t) = b + sum(slope * piece) is additive in linear ",
    "g/cm^2 only if the slopes are interpreted as %/yr-of-baseline ",
    "(the Greendale 2012 convention cited as reference [5] in the ",
    "poster); the model() block does the /100 conversion explicitly."
  )
  reference <- paste(
    "Plan EL, Baron KT, Gastonguay MR, French JL, Gillespie WR, Riggs MM.",
    "Bayesian Joint Modeling of Bone Mineral Density and Repeated",
    "Time-To-Fracture Event for Multiscale Bone Systems Model Extension.",
    "PAGE (Population Approach Group in Europe) 2012 conference poster.",
    "URL https://metrumrg.com/wp-content/uploads/2017/10/Riggs_PAGE2012.pdf",
    "(no DOI assigned; PAGE 2012 abstract IV-04).",
    "Piecewise-linear BMD structure cited from",
    "Greendale GA et al. J Bone Miner Res 2012 (poster reference [5]).",
    sep = " "
  )
  vignette <- "Plan_2012_bmd_fracture"

  # `cumhaz` is the canonical hazard-accumulator compartment registered in
  # R/conventions.R; no paper_specific_compartments declaration needed.

  units <- list(
    time          = "year (years since final menstrual period; t = 0 at FMP)",
    dosing        = "n/a (disease-progression model; no drug input)",
    concentration = "g/cm^2 (BMD observation); survival probability for sur (derived output, unitless)"
  )

  # No covariates are used in model() because the source poster does not
  # report per-parameter coefficient magnitudes for the demographic
  # covariates it retained. The three covariates the authors carried in
  # the final fit ('Centered covariate effects added on all parameters')
  # are documented in covariatesDataExcluded so the provenance is
  # preserved without triggering convention warnings.
  covariateData <- list()

  covariatesDataExcluded <- list(
    BMI = list(
      description        = paste0(
        "Body mass index at the time of the NHANES examination. Listed in ",
        "Plan 2012 Methods ('Included covariates: BMI, ethnicity, and ",
        "FMPage') as one of three covariates added as centered effects ",
        "on every BMD-model structural parameter (b, s_tr, s_po, s_fi). ",
        "Documented here in covariatesDataExcluded because the poster ",
        "does not report the per-parameter coefficient magnitudes; ",
        "encoding fixed(0) entries for all 12 unreported covariate ",
        "effects would add clutter without information."
      ),
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Centred form per the source paper. The reference (centring) ",
        "value is not reported in the poster; the typical NHANES ",
        "postmenopausal-women cohort BMI is approximately 28-30 kg/m^2. ",
        "Per-parameter coefficient magnitudes (e_bmi_b, e_bmi_s_tr, ",
        "e_bmi_s_po, e_bmi_s_fi) are not reported in the source PAGE ",
        "2012 poster and would need to be obtained from the authors' ",
        "underlying WinBUGS run for forward simulation."
      ),
      source_name        = "BMI"
    ),
    AGE_FMP = list(
      description        = paste0(
        "Age at the final menstrual period (FMP), i.e. the subject's ",
        "chronological age in years at the time her FMP occurred. ",
        "Distinct from current age: in the Plan 2012 NHANES cohort the ",
        "mean current age is 63 yr (95% IP 27-85) and the mean age at ",
        "FMP is 45 yr (95% IP 26-57). Listed in Plan 2012 Methods as ",
        "one of three retained covariates (centered effect on every BMD-",
        "model structural parameter). The hazard-model investigation ",
        "also includes FMPage as a candidate covariate (Plan 2012 ",
        "Models: 'Investigated covariates: observed BMD, BMD(t), ",
        "FMPage, and time')."
      ),
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Time-fixed per subject. Derived from the NHANES Reproductive ",
        "Health questionnaire (age at last natural menstrual period). ",
        "Per-parameter coefficient magnitudes are not reported in the ",
        "source PAGE 2012 poster. No canonical column exists in ",
        "inst/references/covariate-columns.md for 'age at FMP' (distinct ",
        "from the canonical AGE = subject current age); the paper-",
        "specific name AGE_FMP is used here in covariatesDataExcluded ",
        "for documentation only and is not validated against the ",
        "canonical register."
      ),
      source_name        = "FMPage"
    ),
    ETHNICITY = list(
      description        = paste0(
        "Self-reported ethnicity using the NHANES 2005-2008 categories ",
        "(Non-Hispanic White, Non-Hispanic Black, Mexican American, ",
        "Other Hispanic, Other Race / Multi-Racial). Listed in Plan 2012 ",
        "Methods as one of three retained BMD-model covariates ('BMI, ",
        "ethnicity, and FMPage'). The specific reference category and ",
        "the per-parameter coefficient magnitudes are not reported in ",
        "the source PAGE 2012 poster."
      ),
      units              = "(categorical)",
      type               = "categorical",
      reference_category = "not reported in the source poster (canonical NHANES practice usually treats Non-Hispanic White as the reference)",
      notes              = paste0(
        "Time-fixed per subject. The poster does not specify how the ",
        "five NHANES ethnicity groups were collapsed for modelling ",
        "(e.g., binary White-vs-non-White, full 4-level dummy, or a ",
        "subset). The paper-specific name ETHNICITY is used here in ",
        "covariatesDataExcluded for documentation only; a re-extraction ",
        "with author-supplied coefficient magnitudes would replace this ",
        "with the canonical RACE_BLACK / RACE_HISPANIC / RACE_ASIAN ",
        "indicators registered in inst/references/covariate-columns.md."
      ),
      source_name        = "ethnicity"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 1605L,
    n_studies       = 1L,
    age_range       = "27-85 years (95% interpercentile range; mean 63 yr at NHANES examination)",
    age_median      = "mean 63 years at examination (Plan 2012 Methods/Data)",
    weight_range    = "not transcribed in this extraction (NHANES 2005-2008 anthropometric distributions are public but not summarized in the poster)",
    sex_female_pct  = 100,
    race_ethnicity  = "NHANES 2005-2008 sample weights span Non-Hispanic White, Non-Hispanic Black, Mexican American, Other Hispanic, and Other Race / Multi-Racial. Per-group counts are not reported in the poster.",
    disease_state   = "Postmenopausal women (any age at FMP; women on the natural-menopause trajectory rather than a surgical / induced menopause sub-cohort).",
    dose_range      = "n/a (this model has no drug input; the disease-progression trajectory is driven entirely by years since FMP).",
    regions         = "United States (NHANES 2005-2008 nationally representative sample).",
    notes           = paste0(
      "Plan 2012 PAGE conference poster (IV-04). Data source: 2005-2008 ",
      "National Health and Nutrition Examination Survey (NHANES) ",
      "demographics, dual-energy X-ray absorptiometry, body measures, ",
      "osteoporosis, and reproductive-health datasets. Cohort: 1605 ",
      "postmenopausal women, mean current age 63 yr (95% IP 27-85), ",
      "mean age at final menstrual period 45 yr (95% IP 26-57). Each ",
      "subject contributes one femoral-neck BMD measurement and 0-5 ",
      "self-reported fracture events (204 total fracture events across ",
      "the cohort). Fit using WinBUGS / BlackBox with mrgsim / deSolve ",
      "in R for simulation diagnostics. The Plan 2012 BMD model ",
      "structure is taken from Greendale 2012 (J Bone Miner Res ",
      "reference [5] in the poster); the structural parameter estimates ",
      "are reported as 'close to reported literature values [5]'."
    ),
    n_fractures_total = 204L
  )

  ini({
    # =====================================================================
    # BMD piecewise-linear structural parameters (Plan 2012, Results,
    # 'Structural parameter estimates' row). The source PAGE 2012 poster
    # reports:
    #     b = 0.84, s_tr = -1.66, s_po = -0.85, s_fi = -0.34
    # with no per-parameter uncertainty in the poster Results block. The
    # poster cites Greendale 2012 (J Bone Miner Res; reference [5]) as
    # the source of the piecewise-linear functional form. The Greendale
    # 2012 convention reports BMD slopes as annualized percent change
    # of baseline BMD (%/yr): taking the printed equation
    #     BMD(t) = b + s_tr * t_(-1,2) + s_po * t_(2,5) + s_fi * t_(5,inf)
    # literally additive in linear g/cm^2 gives BMD(t = 60 yr) ~ -25
    # g/cm^2, which is physically impossible; under the %/yr-of-baseline
    # interpretation (which matches the Greendale 2012 SWAN model
    # reporting style) the model() block applies the slopes as
    # b * (1 + sum(slope_k * piece_k) / 100), reproducing the ~ 0.84 -> 0.62
    # g/cm^2 decline shown in poster Figure 3 over t = 0 -> 60 yr post-FMP.
    # =====================================================================

    b_bmd <- 0.84       ; label("Reference BMD at the start of the menopause transition phase t = -1 yr (g/cm^2)")
    # Plan 2012 Results, 'Structural parameter estimates'; b = 0.84.

    s_tr  <- -1.66      ; label("Annualized percent change of baseline BMD during menopause transition phase (t in [-1, 2] yr) (%/yr; see vignette Errata for unit-interpretation rationale)")
    # Plan 2012 Results, 'Structural parameter estimates'; s_tr = -1.66.

    s_po  <- -0.85      ; label("Annualized percent change of baseline BMD during early postmenopausal phase (t in [2, 5] yr) (%/yr)")
    # Plan 2012 Results, 'Structural parameter estimates'; s_po = -0.85.

    s_fi  <- -0.34      ; label("Annualized percent change of baseline BMD during long-term post-FMP phase (t > 5 yr) (%/yr)")
    # Plan 2012 Results, 'Structural parameter estimates'; s_fi = -0.34.

    # =====================================================================
    # Repeated time-to-fracture hazard parameters. Plan 2012 reports the
    # final-model hazard form
    #     h(t) = exp(theta_h * (1 + theta_bmd * (BMD(t) - 0.8)))
    # with the right-hand side caption text 'theta_h = -5.5 (95%CI -5.49
    # to -5.54), theta_bmd = 1.5 (95%CI 1.49-1.52)', BMD reference value
    # 0.8 g/cm^2, and Weibull shape alpha = 1 (i.e., exponential RTTE
    # between events; the poster Models section lists alpha != 1 only as
    # an investigated covariate that did not survive into the final
    # model). The hazard accumulates from t = 0 (FMP); see Plan 2012
    # caption to Figure 5.
    # =====================================================================

    theta_h   <- -5.5   ; label("Log baseline fracture hazard at BMD = bmd_ref (log(1/yr); h(t = 0) = exp(theta_h * 1) if BMD(0) = bmd_ref)")
    # Plan 2012 Results, Figure 5 right-panel caption: 'theta_h = -5.5 (95% CI -5.49 to -5.54)'.

    theta_bmd <- 1.5    ; label("BMD effect on hazard: unitless modulation of theta_h via (1 + theta_bmd * (BMD - bmd_ref))")
    # Plan 2012 Results, Figure 5 right-panel caption: 'theta_BMD = 1.5 (95% CI 1.49-1.52)'.

    bmd_ref   <- fixed(0.8)
    label("Reference BMD for hazard centring (g/cm^2; fixed at the cohort-typical postmenopausal BMD per Plan 2012 Figure 4 caption)")
    # Plan 2012 Results, Figure 4 caption: 'BMD = 0.8 g/cm^2, alpha = 1'.

    # =====================================================================
    # Residual variability. Plan 2012 Results, BMD panel: 'Random effect
    # included as residual variability: sigma = 0.131 (95% credible
    # interval (CI) 0.127-0.136)'. Only one femoral-neck BMD measurement
    # per subject precludes separating IIV from RUV, so the published
    # sigma is the combined random-effect + residual SD on the linear
    # BMD scale (additive in g/cm^2). No separate eta* parameters are
    # added.
    # =====================================================================

    addSd <- 0.131      ; label("Additive residual SD on BMD scale (g/cm^2; combined IIV + RUV per Plan 2012 BMD-panel Results)")
    # Plan 2012 Results, BMD panel: 'sigma = 0.131 (95% CI 0.127-0.136)'.

    # =====================================================================
    # Not encoded as primary ini() parameters:
    #   - Weibull shape alpha is fixed at 1 in the final RTTE model (the
    #     hazard is therefore exponential between events) so the shape
    #     parameter does not appear in h(t). The poster Models section
    #     records 'time (alpha != 1, Weibull distribution)' as an
    #     investigated-but-not-retained covariate.
    #   - BMD-model covariates (BMI, NHANES ethnicity, age at FMP):
    #     retained in the source fit as 'centered covariate effects added
    #     on all parameters' but per-parameter coefficient magnitudes are
    #     not reported in the poster. Documented in
    #     covariatesDataExcluded; see also vignette Errata.
    # =====================================================================
  })

  model({
    # =====================================================================
    # 1. Piecewise-linear pieces. Each piece_<k>(t) is the linear-growth
    #    contribution of phase k to BMD(t), clamped to the phase's
    #    duration so the function is continuous across the breakpoints
    #    and constant outside the phase's [a_k, b_k] support:
    #        piece_k(t) = max(0, min(t - a_k, b_k - a_k))
    #    where (a_k, b_k) = (-1, 2) for the transition phase, (2, 5) for
    #    the early postmenopausal phase, and (5, +inf) for the long-term
    #    phase. ifelse is used because rxode2's parser does not accept
    #    the parallel-vector pmin / pmax functions inside model().
    # =====================================================================

    piece_tr <- ifelse(time < -1, 0, ifelse(time > 2, 3, time + 1))
    piece_po <- ifelse(time <  2, 0, ifelse(time > 5, 3, time - 2))
    piece_fi <- ifelse(time <  5, 0, time - 5)

    # =====================================================================
    # 2. BMD time-course. Slopes are annualized percent change of
    #    baseline BMD (Greendale 2012 / Plan 2012 convention); the /100
    #    conversion is applied here so the ini() values match the
    #    poster's reported numbers one-for-one. The resulting
    #        BMD(t) = b_bmd * (1 + sum(slope_k * piece_k) / 100)
    #    reproduces the ~ 0.84 -> 0.62 g/cm^2 decline over t = 0 -> 60 yr
    #    in poster Figure 3.
    # =====================================================================

    BMD <- b_bmd * (1 + (s_tr * piece_tr + s_po * piece_po + s_fi * piece_fi) / 100)

    # =====================================================================
    # 3. Fracture hazard. Plan 2012 final-model RTTE form, Figure 5
    #    right-panel caption:
    #        h(t) = exp(theta_h * (1 + theta_bmd * (BMD(t) - bmd_ref)))
    #    At BMD = bmd_ref (= 0.8 g/cm^2) the modulating factor is 1 and
    #    h = exp(theta_h) = exp(-5.5) ~ 4.1e-3 events/yr. For BMD <
    #    bmd_ref (low-BMD, high-fracture-risk regime), the factor is < 1
    #    so theta_h * factor is less negative and the hazard rises. For
    #    BMD > bmd_ref the hazard falls.
    # =====================================================================

    hazard <- exp(theta_h * (1 + theta_bmd * (BMD - bmd_ref)))

    # =====================================================================
    # 4. Cumulative-hazard ODE. The poster Results section ('Fracture
    #    risk') states 'h(t) accumulates from t_0 = FMP'. The cumulative
    #    hazard Lambda(t) = integral from 0 to t of h(u) du is the
    #    quantity that determines the survival function S(t) =
    #    exp(-Lambda(t)^alpha) with alpha = 1 in the final model.
    # =====================================================================

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0

    # =====================================================================
    # 5. Derived outputs. sur is the marginal survival probability up to
    #    time t (i.e., the probability of remaining fracture-free from
    #    FMP through time t). For forward simulation of repeated fracture
    #    events, users sample u ~ Uniform(0, 1) and solve
    #    cumhaz(t_event) = -log(u) for the next event time, then reset
    #    cumhaz to 0 and continue.
    # =====================================================================

    sur <- exp(-cumhaz)

    # =====================================================================
    # 6. BMD observation with additive residual error. Plan 2012 fit on
    #    one femoral-neck BMD measurement per subject; the combined
    #    additive SD is 0.131 g/cm^2.
    # =====================================================================

    BMD ~ add(addSd)
  })
}
