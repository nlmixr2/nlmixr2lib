Hamuro_2017_DMD_6MWT <- function() {
  description <- paste0(
    "Natural-history disease-progression model for the six-minute walk ",
    "test (6MWT, meters) in ambulatory boys with Duchenne muscular ",
    "dystrophy (DMD) on stable corticosteroids (Hamuro 2017). The 6MWT ",
    "vs subject age is modelled as the minimum of two simultaneously ",
    "estimated linear lines (Phoenix NLME 'min' function): a ",
    "developmental line with a positive slope (improvement) and a ",
    "disease-induced line with a negative slope (decline). The two ",
    "lines intersect at the population-typical age of maximum 6MWT (10 ",
    "years). Exponential between-subject variability is estimated on ",
    "both slopes (the intercepts have no IIV); residual error is ",
    "additive. Disease-progression model with no drug dosing."
  )
  reference <- paste(
    "Hamuro L, Chan P, Tirucherai G, AbuTarif M.",
    "Developing a Natural History Progression Model for Duchenne",
    "Muscular Dystrophy Using the Six-Minute Walk Test.",
    "CPT Pharmacometrics Syst Pharmacol. 2017 Sep;6(9):596-603.",
    "doi:10.1002/psp4.12220.",
    sep = " "
  )
  vignette <- "Hamuro_2017_DMD_6MWT"
  units <- list(
    time          = "year",
    dosing        = "n/a (disease-progression model with no drug dosing)",
    concentration = "m (six-minute walk test distance, observation walkDist)"
  )

  covariateData <- list(
    AGE = list(
      description        = "Subject age in years at the time of the 6MWT measurement (time-varying, updated at each visit).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Model 4 fits 6MWT directly against the subject's age at the time of each measurement, not against time-since-baseline. For simulation, supply AGE on each observation row; the model uses it without further transformation. Hamuro 2017 Table 1 reports a mean baseline age of about 8-10 years across the contributing studies; the model was developed against a baseline-age range of approximately 4-17 years (Hamuro 2017 Table 1) and the age range covered by observations extends a year or two beyond the baseline.",
      source_name        = "Age"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 88L,
    n_studies      = 2L,
    age_range      = "approximately 4-17 years (combined cohort, observations); baseline-age range 5-15.3 years",
    age_median     = "mean baseline 8.3 years (McDonald 2013) and 9.5 years (Goemans 2013)",
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = 0,
    race_ethnicity = NA_character_,
    disease_state  = "Ambulatory boys with Duchenne muscular dystrophy (DMD) on stable corticosteroids. 100% of subjects from Goemans 2013 (90% deflazacort, 10% prednisone/prednisolone) and 70% of subjects from McDonald 2013 (steroid type and regimen not reported).",
    dose_range     = "Not applicable -- disease-progression model with no drug dosing. Steroid use is implicit in the population.",
    regions        = "International (McDonald 2013 multi-site randomised controlled trial placebo arm) and Belgium (Goemans 2013 Leuven Neuromuscular Reference Center natural-history study)",
    n_observations = "228 6MWT records across 88 subjects (Goemans 2013: 122 records / 35 subjects, 2-6 per subject; McDonald 2013: 106 records / 53 subjects, 2 per subject at baseline and 48 weeks). 13 of 228 records (6%) were below the 50 m quantification limit and handled via M3 in the original fit; the packaged model does not include the M3 censoring component.",
    notes          = "All data were digitised from published figures (plot digitizer) rather than obtained as individual subject records. Hamuro 2017 Methods (Data sources) and Table 1 describe the contributing studies (McDonald 2013 Figure 4, Goemans 2013 Figure 1)."
  )

  ini({
    # ============================================================
    # Final model: 'Linear with simultaneous estimation' (Hamuro
    # 2017 Methods, page 597, Model 4). For each subject and at each
    # age the predicted 6MWT is the minimum of two simultaneously
    # estimated lines, walk_dev (positive slope) and walk_decl
    # (negative slope). Between-subject variability is exponential
    # on both slopes; the intercepts are typical-value-only (Table
    # 2 BSV %CV reported as N/A). Residual error is additive (Table
    # 2 stdev). Parameter estimates from the final selected model
    # are reproduced from Hamuro 2017 Table 2 (column 'Estimate').
    # ============================================================

    # ---------------- Developmental line ----------------
    intercept    <- 270             ; label("Developmental intercept at age 0 (m)")                   # Hamuro 2017 Table 2 (Estimate 270 m, 95% CI 197-344)
    lslope_dev   <- log(19.6)       ; label("Log developmental slope magnitude (m/year)")             # Hamuro 2017 Table 2 (Estimate 19.6 m/year, 95% CI 9.4-29.8)

    # ---------------- Disease-induced line ----------------
    # Stored as the log magnitude (positive). The negative sign is
    # applied in model() so that slope_decl < 0 matches Table 2's
    # signed -84.9 m/year, and exponential IIV multiplies the
    # magnitude rather than introducing a sign change.
    intercept2   <- 1298            ; label("Disease-induced intercept at age 0 (m)")                 # Hamuro 2017 Table 2 (Estimate 1298 m, 95% CI 1158-1437)
    lslope_decl  <- log(84.9)       ; label("Log disease-induced slope magnitude (m/year)")           # Hamuro 2017 Table 2 (Estimate -84.9 m/year, 95% CI -97.6 to -72.2; magnitude stored)

    # ---------------- Between-subject variability ----------------
    # Hamuro 2017 Methods: 'Interpatient variability modeled as
    # exponential' for Model 4. Variance on the natural-log scale
    # is log(1 + CV^2). Table 2 reports BSV %CV of 22 on the
    # developmental slope and 23 on the disease-induced slope.
    etalslope_dev  ~ 0.04727                                                                          # Hamuro 2017 Table 2 BSV %CV 22 -> log(1 + 0.22^2)
    etalslope_decl ~ 0.05153                                                                          # Hamuro 2017 Table 2 BSV %CV 23 -> log(1 + 0.23^2)

    # ---------------- Residual error ----------------
    addSd <- 56.9                   ; label("Additive residual SD on 6MWT (m)")                       # Hamuro 2017 Table 2 stdev (Estimate 56.9 m, 95% CI 52.2-61.6)
  })

  model({
    # ----- Individual structural parameters -----
    # Developmental and disease-induced slopes carry exponential
    # IIV on their magnitudes (positive). The sign on the
    # disease-induced slope is applied at the slope level so that
    # the line trends downward in age.
    slope_dev  <-  exp(lslope_dev  + etalslope_dev)
    slope_decl <- -exp(lslope_decl + etalslope_decl)

    # ----- Two linear lines in age -----
    walk_dev  <- intercept  + slope_dev  * AGE
    walk_decl <- intercept2 + slope_decl * AGE

    # ----- Simultaneous minimum (Hamuro 2017 Methods, Model 4) -----
    # 6MWT = min(walk_dev, walk_decl). The minimum is realised here
    # using a sign-based selector that is equivalent to the
    # `min()` operator at every age except the (measure-zero) point
    # of intersection where the two lines meet. The paper notes
    # that this point of non-differentiability is acceptable for
    # the expectation-maximisation (QRPEM) estimator that the
    # original fit used.
    walkDist <- (walk_dev <= walk_decl) * walk_dev +
                (walk_dev >  walk_decl) * walk_decl

    walkDist ~ add(addSd)
  })
}
