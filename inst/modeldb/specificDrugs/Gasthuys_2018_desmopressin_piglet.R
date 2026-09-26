Gasthuys_2018_desmopressin_piglet <- function() {
  description <- paste(
    "Preclinical (pig). Two-compartment population PK model with a dual,",
    "parallel input function for a 120 ug desmopressin sublingual oral",
    "lyophilisate (Minirin Melt) in growing piglets aged 8 days to 6 months,",
    "used as a juvenile animal model for the pediatric population",
    "(Gasthuys 2018). A fraction Bio of the dose is released zero-order over",
    "a duration D1 into a first depot representing buccal absorption and then",
    "reaches the central compartment first-order at Ka1; the remaining",
    "1 - Bio is assumed to be swallowed and enters a second depot after a",
    "1 h lag, from which it is absorbed first-order at Ka2. The dual input",
    "reproduces the second peak seen in most piglet plasma profiles.",
    "Elimination is linear. Absolute oral bioavailability could not be",
    "estimated because no intravenous data were collected, so clearance and",
    "both volumes are apparent (CL/F, V1/F, V2/F). Body weight is the only",
    "retained covariate, entering CL and V1 as a power function centred on",
    "10 kg."
  )
  reference <- paste(
    "Gasthuys E, Vermeulen A, Croubels S, Millecam J, Schauvliege S,",
    "van Bergen T, De Bruyne P, Vande Walle J, Devreese M.",
    "Population Pharmacokinetic Modeling of a Desmopressin Oral Lyophilisate",
    "in Growing Piglets as a Model for the Pediatric Population.",
    "Front Pharmacol. 2018;9:41. doi:10.3389/fphar.2018.00041.",
    sep = " "
  )
  vignette <- "Gasthuys_2018_desmopressin_piglet"

  # Both depots are dosing targets: the same sublingual dose is split between
  # `depot` (buccal, zero-order release) and `depot2` (swallowed fraction,
  # lagged). Declared explicitly because automatic detection only recognises
  # `depot` and `central`.
  dosing <- c("depot", "depot2")

  # The administered dose is 120 ug and all volumes are reported in L, so the
  # states are ug and central/vc is ug/L. Plasma concentrations are reported
  # throughout the paper in pg/mL (LOQ 4.2 pg/mL, Cmax about 250 pg/mL in the
  # 8-day-old group), so the observation carries an explicit 1000 pg/mL per
  # ug/L conversion; no parameter value is rescaled.
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "pg/mL"
  )

  # Issue #482: what each ODE state holds. All four states are amounts of
  # desmopressin in ug (Figure 3, model scheme).
  compartmentData <- list(
    depot = list(
      analyte = "desmopressin",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    depot2 = list(
      analyte = "desmopressin",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "desmopressin",
      units = "ug",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "desmopressin",
      units = "ug",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Centred on the typical weight of 10 kg (Table 2 caption). Power",
        "function on both CL/F and V1/F, retained after forward inclusion",
        "and backward elimination (dOFV = 14.43 and 39.72 respectively).",
        "Recorded once per animal at the single study occasion, so it is",
        "constant within a subject.",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  # Screened during stepwise covariate modelling on CL/F, V1/F and the
  # absorption parameters but not retained in the final model (Results,
  # "Pharmacokinetic Analysis"). GFR reached significance on forward
  # inclusion but failed backward elimination, so no point estimate is
  # published for any of these and none is referenced in model().
  covariatesDataExcluded <- list(
    BSA = list(
      description = "Body surface area, determined by 3D CT reconstruction",
      units = "m^2",
      type = "continuous",
      notes = "Screened on CL/F, V1/F and absorption parameters; not retained.",
      source_name = "BSA"
    ),
    AGE = list(
      description = "Age of the piglet",
      units = "years",
      type = "continuous",
      notes = paste(
        "Screened as a covariate; not retained. Reported in the paper by",
        "age group (8 days, 4 weeks, 7 weeks, 6 months) rather than as a",
        "continuous column.",
        sep = " "
      ),
      source_name = "age"
    ),
    CRCL = list(
      description = "Glomerular filtration rate, measured as exo-iohexol clearance",
      units = "mL/min/m^2",
      type = "continuous",
      notes = paste(
        "Screened on CL/F. Significant on forward inclusion but did not",
        "reach the backward-elimination criterion, so it was dropped and no",
        "coefficient is published (Supplementary Table 1 run 16). Note the",
        "paper normalises GFR per m^2 of body surface area, not per",
        "1.73 m^2 as the canonical CRCL register entry does.",
        sep = " "
      ),
      source_name = "GFR"
    )
  )

  population <- list(
    species = "pig (Landrace x large white piglet, Seghers Hybrid)",
    n_subjects = 32,
    n_studies = 1,
    age_range = "8 days to 6 months (four groups of n = 8: 8 days, 4 weeks, 7 weeks, 6 months)",
    weight_range = "1.54-124 kg",
    weight_median = "2.21 kg (8 days), 10.0 kg (4 weeks), 15.0 kg (7 weeks), 113.6 kg (6 months)",
    sex_female_pct = 0,
    disease_state = "healthy growing piglets used as a juvenile animal model for the pediatric population",
    dose_range = "120 ug desmopressin oral lyophilisate placed sublingually, single dose",
    regions = "Belgium (Ghent University)",
    notes = paste(
      "Table 1 gives mean +/- SD and median [range] for body weight, body",
      "surface area and GFR by age group. All animals were male, catheterised",
      "in the jugular vein, fasted 1 h before and 1.5 h after dosing, and",
      "sampled richly (0, 5, 15, 30, 60 min and 1.5, 2, 3, 4, 6, 8, 10, 12,",
      "24 h). 8% of concentrations were below the 4.2 pg/mL LOQ and were",
      "excluded. Note that the Methods text reports the 7-week group as",
      "13.9 +/- 2.74 kg while Table 1 reports 15.8 +/- 1.98 kg; Table 1 is",
      "used here.",
      sep = " "
    )
  )

  ini({
    # Structural model, Table 2. Body weight is centred on 10 kg, so these
    # are the typical values for a 10 kg piglet. F could not be estimated
    # (no intravenous data), so CL and both volumes are apparent.
    lcl <- log(395)
    label("Apparent total body clearance CL/F (L/h)") # Table 2: theta1 = 395 L/h (SE 31.8, RSE 8.05%)
    lvc <- log(131)
    label("Apparent central volume of distribution V1/F (L)") # Table 2: theta2 = 131 L (SE 21.1, RSE 16.1%)
    lka <- log(0.275)
    label("First-order absorption rate constant from the buccal depot Ka1 (1/h)") # Table 2: theta3 = 0.275 1/h (SE 0.0272, RSE 9.89%)
    lq <- log(32)
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2: theta4 = 32 L/h (SE 6.8, RSE 21.3%)
    lvp <- log(436)
    label("Apparent peripheral volume of distribution V2/F (L)") # Table 2: theta5 = 436 L (SE 137, RSE 31.4%)
    ld1 <- log(0.16)
    label("Duration of the zero-order release into the buccal depot D1 (h)") # Table 2: theta6 = 0.16 h (SE 0.0473, RSE 29.6%)
    lka2 <- log(0.399)
    label("First-order absorption rate constant from the swallowed depot Ka2 (1/h)") # Table 2: theta7 = 0.399 1/h (SE 0.0677, RSE 17.0%)
    logitfdepot <- logit(0.86)
    label("Fraction of the dose passing through the buccal depot Bio (unitless)") # Table 2: theta8 = 86% (SE 0.0488, RSE 5.67%)
    ltlag <- log(1)
    label("Lag time before the swallowed dose reaches the second depot Tlag (h)") # Table 2: theta9 = 1 h (SE 0.00196, RSE 0.20%)

    # Covariate model, Table 2. Power exponents on (BW/10).
    e_wt_cl <- 1.03
    label("Power exponent on (WT/10) for CL/F (unitless)") # Table 2 'Influence of BW on CL' = 1.03 (SE 0.0627, RSE 6.09%)
    e_wt_vc <- 0.691
    label("Power exponent on (WT/10) for V1/F (unitless)") # Table 2 'Influence of BW on V1' = 0.691 (SE 0.135, RSE 19.54%)

    # Inter-individual variability, Table 2. The printed column is an omega
    # VARIANCE, not a CV: for every IIV row the printed RSE equals
    # (SE / estimate) / 2, which is the RSE on the SD scale of a variance
    # estimate (e.g. CL 0.0462 / 0.175 / 2 = 13.2%, V1 0.264 / 0.641 / 2 =
    # 20.6%, Bio 0.391 / 0.627 / 2 = 31.2%). IIV on Q/F, Ka2 and Tlag was
    # not estimable and was fixed to zero, so those etas are omitted here.
    etalcl ~ 0.175 # Table 2 'IIV CL/F' = 0.175 (shrinkage 1.97%, RSE 13.20%)
    etalvc ~ 0.641 # Table 2 'IIV V1/F' = 0.641 (shrinkage 25.4%, RSE 20.59%)
    etalka ~ 0.0903 # Table 2 'IIV Ka1' = 0.0903 (shrinkage 16.0%, RSE 18.77%)
    etalvp ~ 0.634 # Table 2 'IIV V2/F' = 0.634 (shrinkage 61.8%, RSE 15.93%)
    etald1 ~ 0.485 # Table 2 'IIV D1' = 0.485 (shrinkage 47.4%, RSE 26.29%)
    # Bio is bounded in (0, 1) because the complement 1 - Bio is the dose
    # fraction entering the second depot, so its eta is carried on the logit
    # scale rather than the exponential scale used for the unbounded
    # parameters above. See the model file note below Table 2 discussion and
    # the vignette Assumptions section for the two lines of evidence.
    etalogitfdepot ~ 0.627 # Table 2 'IIV Bio' = 0.627 (shrinkage 42.9%, RSE 31.18%)

    # Residual error. The paper modelled log-transformed concentrations with
    # an additive error ('log-transform both sides'), which the Methods
    # state is a proportional error model in the untransformed domain. The
    # estimate is reported directly on the SD scale: 0.0291 / 0.228 = 12.76%
    # reproduces the printed RSE without a factor of two.
    propSd <- 0.228
    label("Proportional residual error (fraction)") # Table 2 'ADD = theta10' = 22.8% CV (SE 0.0291, RSE 12.76%)
  })

  model({
    # ----------------------------------------------------------------
    # 1. Individual parameters, Table 2 structural-model column. Body
    #    weight is centred on 10 kg (Table 2 caption).
    # ----------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (WT / 10)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 10)^e_wt_vc
    ka <- exp(lka + etalka)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)
    d1 <- exp(ld1 + etald1)
    ka2 <- exp(lka2)
    tlag <- exp(ltlag)
    fdepot <- expit(logitfdepot + etalogitfdepot)

    # ----------------------------------------------------------------
    # 2. Micro-constants.
    # ----------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----------------------------------------------------------------
    # 3. ODE system, Figure 3. Two parallel absorption paths feed one
    #    central compartment:
    #      depot  - buccal; receives Bio * Dose released zero-order over
    #               D1 and empties first-order at Ka1. Dose records must
    #               carry rate = -2 for rxode2 to use dur(depot); a plain
    #               bolus would silently ignore D1.
    #      depot2 - swallowed; receives (1 - Bio) * Dose after a Tlag
    #               delay and empties first-order at Ka2.
    # ----------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka2 * depot2
    d/dt(central) <- ka * depot + ka2 * depot2 -
      kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ----------------------------------------------------------------
    # 4. Dose split, zero-order release duration and lag time.
    # ----------------------------------------------------------------
    f(depot) <- fdepot
    dur(depot) <- d1
    f(depot2) <- 1 - fdepot
    alag(depot2) <- tlag

    # ----------------------------------------------------------------
    # 5. Observation. States are ug and vc is L, so central / vc is ug/L;
    #    the factor of 1000 converts to the pg/mL used throughout the paper.
    # ----------------------------------------------------------------
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
