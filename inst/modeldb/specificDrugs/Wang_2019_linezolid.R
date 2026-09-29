Wang_2019_linezolid <- function() {
  description <- paste0(
    "One-compartment population PK model with intravenous administration ",
    "and first-order elimination for linezolid 600 mg every 12 h in ",
    "Chinese critically ill adults with and without shock (Wang 2019). ",
    "Clearance carries a power-form platelet-count effect referenced to ",
    "200 x 10^9/L with an exponent of 0.261, so clearance rises with ",
    "platelet count. The paper's headline negative result is that shock ",
    "type (no shock, septic, hemorrhagic, neurogenic, cardiogenic) was ",
    "not retained on clearance; of seventeen screened covariates only ",
    "platelet count met the inclusion criterion. Inter-individual ",
    "variability is estimated on both clearance and volume; residual ",
    "error is proportional and very large (about 100% CV), reflecting ",
    "opportunistic therapeutic-drug-monitoring sampling in a 37-patient ",
    "single-centre cohort."
  )
  reference <- paste0(
    "Wang D, Zheng X, Yang Y, Chen X. Population pharmacokinetic analysis ",
    "of linezolid in patients with different types of shock: Effect of ",
    "platelet count. Exp Ther Med. 2019;18(2):1786-1792. ",
    "doi:10.3892/etm.2019.7747. PMCID: PMC6676194."
  )
  vignette <- "Wang_2019_linezolid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    PLT = list(
      description = "Platelet count",
      units = "10^9 cells/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Enters clearance as the power function (PLT / 200)^0.261 (Wang ",
        "2019 Equation e). The centring constant 200 is printed inside ",
        "the equation and is the rounded cohort median: Wang 2019 Table I ",
        "gives the observed median platelet count as 213 x 10^9/L (range ",
        "11-895), and the Methods state that 'Covconstant was fixed at a ",
        "value similar to the population median of the covariate', so 200 ",
        "is a deliberately rounded median rather than a transcription of ",
        "it. The exponent is positive: clearance increases with platelet ",
        "count, spanning a 2.2-fold range across the observed 11-895 ",
        "x 10^9/L window. The paper does not state whether the count is ",
        "the admission value or a time-varying repeat measurement; since ",
        "the concentrations came from routine therapeutic drug ",
        "monitoring, treat it as the value contemporaneous with the ",
        "sampled course. Platelet count was the only one of seventeen ",
        "screened covariates to meet the forward-inclusion criterion ",
        "(dOFV -11.101, Table II), and it survived backward elimination ",
        "(dOFV +11.101, P < 0.01). Its mechanistic interpretation is ",
        "confounded in this cohort: linezolid causes thrombocytopenia, ",
        "and linezolid accumulation from low clearance is the recognised ",
        "driver of that toxicity, so a low platelet count may be a ",
        "consequence of low clearance rather than a cause of it. The ",
        "model is therefore descriptive, not causal."
      ),
      source_name = "PLT"
    )
  )

  # Covariates Wang 2019 screened (Table II) but did NOT retain in the final
  # model. Documentation only: none of these is referenced in model(). Only
  # sex reached P < 0.05 in step 1 (dOFV -5.262) and it failed the step-2
  # test once platelet count was in the model (dOFV -2.181, P > 0.05).
  # Four screened analytes -- globulin, the albumin/globulin ratio, mean
  # corpuscular hemoglobin and mean corpuscular hemoglobin concentration --
  # have no canonical entry in inst/references/covariate-columns.md; because
  # they were not retained, no new canonical is minted for them here. They
  # are recorded with their dOFV values in population$covariate_screen below,
  # together with the four shock strata other than sepsis.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = paste0(
        "Screened on CL. Wang 2019 Table II: dOFV -5.262 (P < 0.05) in ",
        "forward step 1, the only covariate other than platelet count to ",
        "pass; dropped in step 2 because adding it on top of platelet ",
        "count gave only dOFV -2.181 (P > 0.05). The cohort was 27 male / ",
        "10 female (Table I). The paper reports no point estimate for a ",
        "sex effect, so none can be encoded."
      ),
      source_name = "Sex"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL; Wang 2019 Table II dOFV -0.759 (P > 0.05). Cohort median 62 years (range 29-89, Table I).",
      source_name = "Age"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL; Wang 2019 Table II dOFV -2.461 (P > 0.05). Cohort median 32.4 g/L (range 25.4-43.8, Table I).",
      source_name = "ALB"
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL; Wang 2019 Table II dOFV -0.060 (P > 0.05). Cohort median 62.8 g/L (range 45.9-75.1, Table I).",
      source_name = "TP"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL; Wang 2019 Table II dOFV -0.015 (P > 0.05). Cohort median 48 IU/L (range 3-525, Table I).",
      source_name = "ALT"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL; Wang 2019 Table II dOFV -0.019 (P > 0.05). Cohort median 41 IU/L (range 11-181, Table I).",
      source_name = "AST"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = paste0(
        "Screened on CL; Wang 2019 Table II dOFV -0.011 (P > 0.05). ",
        "Cohort median 85 umol/L (range 16-499, Table I). Linezolid is ",
        "predominantly non-renally cleared, so the null result is ",
        "mechanistically unsurprising despite the wide renal-function ",
        "range in the cohort."
      ),
      source_name = "SCR"
    ),
    BUN = list(
      description = "Blood urea",
      units = "mmol/L",
      type = "continuous",
      notes = paste0(
        "Screened on CL; Wang 2019 Table II dOFV -0.230 (P > 0.05). ",
        "Cohort median 11.1 mmol/L (range 2.4-29.9, Table I). Reported by ",
        "the paper as 'urea' in mmol/L, i.e. whole-molecule urea rather ",
        "than the urea-nitrogen basis that the mg/dL BUN convention uses; ",
        "divide by 0.357 to convert mg/dL BUN to mmol/L urea."
      ),
      source_name = "Urea"
    ),
    TBA = list(
      description = "Total serum bile acids",
      units = "umol/L",
      type = "continuous",
      notes = "Screened on CL; Wang 2019 Table II dOFV -3.489 (P > 0.05, i.e. short of the 3.84 threshold). Cohort median 3.8 umol/L (range 1.0-81.5, Table I).",
      source_name = "TBA"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened on CL; Wang 2019 Table II dOFV -3.140 (P > 0.05). Cohort median 11.1 umol/L (range 2.0-414.7, Table I).",
      source_name = "TBIL"
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      notes = "Screened on CL; Wang 2019 Table II dOFV -0.964 (P > 0.05). Cohort median 26.5% (range 19.0-44.8, Table I).",
      source_name = "HCT"
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL; Wang 2019 Table II dOFV -0.497 (P > 0.05). Cohort median 87 g/L (range 63-140, Table I).",
      source_name = "HGB"
    ),
    DIS_SEPSIS = list(
      description = "Septic shock indicator (1 = septic shock)",
      units = "(binary)",
      type = "binary",
      notes = paste0(
        "Screened on CL; Wang 2019 Table II dOFV -2.038 (P > 0.05). ",
        "11 of 37 patients. The canonical DIS_SEPSIS marks active sepsis; ",
        "Wang 2019's stratum is specifically septic SHOCK, a more severe ",
        "subset. The three other shock strata (hemorrhagic n = 1, ",
        "neurogenic n = 4, cardiogenic n = 3) and the no-shock stratum ",
        "(n = 18) have no canonical register entry and, having been ",
        "rejected, none is minted; their dOFV values are recorded in ",
        "population$covariate_screen."
      ),
      source_name = "Septic shock"
    )
  )

  compartmentData <- list(
    central = list(analyte = "linezolid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 37L,
    n_studies = 1L,
    age_range = "29-89 years",
    age_median = "62 years",
    sex_female_pct = 27.0,
    race_ethnicity = c(Asian = 100),
    disease_state = paste0(
      "Critically ill adults treated with linezolid, stratified by shock ",
      "type: no shock (n = 18), septic shock (n = 11), hemorrhagic shock ",
      "(n = 1), neurogenic shock (n = 4), cardiogenic shock (n = 3)"
    ),
    dose_range = "Linezolid 600 mg intravenously every 12 h (a single regimen; no dose-ranging)",
    regions = "China (single centre, Zhongda Hospital affiliated to Southeast University, Nanjing)",
    renal_function = "Wide: serum creatinine median 85 umol/L, range 16-499 (Table I). Renal function was screened on clearance and not retained.",
    hepatic_function = "Wide: ALT median 48 IU/L (range 3-525), total bilirubin median 11.1 umol/L (range 2.0-414.7) (Table I).",
    co_medication = "Patients taking drugs expected to interact with linezolid pharmacokinetics were excluded by design (Methods, 'Patients and data collection').",
    bioanalytical = "Retrospective therapeutic drug monitoring concentrations abstracted from hospital records; assay method and lower limit of quantification are not reported.",
    # Full Wang 2019 Table II forward-inclusion screen, in the paper's own
    # order. Recorded here rather than in covariatesDataExcluded because
    # four of the analytes and four of the shock strata have no canonical
    # covariate-column name and, being rejected, none is minted for them.
    # dOFV is the change in objective function value relative to the base
    # model (OFV 1080.738); the inclusion threshold was dOFV < -3.84.
    covariate_screen = c(
      sex = -5.262,
      age = -0.759,
      no_shock = -0.118,
      septic_shock = -2.038,
      hemorrhagic_shock = -1.605,
      neurogenic_shock = -0.125,
      cardiogenic_shock = -0.345,
      albumin = -2.461,
      globulin = -0.011,
      albumin_globulin_ratio = -0.384,
      alanine_transaminase = -0.015,
      aspartate_transaminase = -0.019,
      serum_creatinine = -0.011,
      urea = -0.230,
      total_protein = -0.060,
      total_bile_acid = -3.489,
      total_bilirubin = -3.140,
      platelets = -11.101,
      hematocrit = -0.964,
      hemoglobin = -0.497,
      mean_corpuscular_hemoglobin = -3.570,
      mean_corpuscular_hemoglobin_concentration = -0.526
    ),
    notes = paste0(
      "Retrospective single-centre analysis of routine therapeutic drug ",
      "monitoring records collected January 2016 to August 2018 (Wang ",
      "2019 Methods; demographics in Table I). Weight was not recorded ",
      "among the screened covariates and no weight range is reported, so ",
      "the model carries no allometric term and its parameters are ",
      "absolute (L/h, L) rather than per-kilogram. The number of ",
      "concentration records is not stated; Figure 3 shows sampling from ",
      "roughly 50 to 820 h after the start of therapy, i.e. at steady ",
      "state under the fixed 600 mg every 12 h regimen. Model evaluation ",
      "used a 1,000-replicate bootstrap (992 successful) and a ",
      "prediction-corrected visual predictive check."
    )
  )

  ini({
    # ----- Structural parameters (Wang 2019 Table III, 'Estimate' column) -----
    # The typical clearance is for a patient with a platelet count of
    # 200 x 10^9/L, at which the covariate factor equals 1.
    lcl <- log(11.8)
    label("Clearance CL at PLT = 200 x 10^9/L (L/h)")  # Wang 2019 Table III: CL 11.8 L/h (bootstrap median 11.2, 95% CI 3.110-16.625)
    lvc <- log(209)
    label("Central volume of distribution V (L)")  # Wang 2019 Table III: V 209 L (bootstrap median 197, 95% CI 53.850-308.000)

    # ----- Covariate effect on CL -----
    e_plt_cl <- 0.261
    label("Power exponent of (PLT / 200) on CL (unitless)")  # Wang 2019 Table III: theta PLT 0.261 (95% CI 0.052-0.425); Equation (e) CL = theta_CL x (PLT/200)^theta_PLT

    # ----- Inter-individual variability -----
    # Wang 2019 Table III reports omega_CL and omega_V with identical point
    # estimates (0.299), identical bootstrap medians (0.287) and identical
    # bias (-4.013%); only the bootstrap confidence limits differ in the
    # third decimal. They are taken as two separate diagonal OMEGA elements
    # that happened to converge to the same value, which is what the table
    # states; no covariance between them is reported, so the block is
    # diagonal.
    #
    # SCALE: the values are read as NONMEM OMEGA variances, not standard
    # deviations, giving 59.0% CV on each parameter. Two independent checks
    # support this. (1) Wang 2019 Figure 2A shows individual predictions
    # spanning about 1.2 to 17.2 mg/L. Under the fixed 600 mg every 12 h
    # regimen at steady state the between-subject part of that spread is
    # driven by 1/CL, and with a within-subject peak-to-trough ratio near 2
    # the between-subject range is about 5.6-fold = exp(1.72); across 37
    # subjects an expected range of about 4.1 omega implies omega near 0.42
    # BEFORE accounting for empirical-Bayes shrinkage, which with this
    # model's very large residual error is substantial and biases the
    # observable spread downwards -- so the true omega is above 0.42, i.e.
    # sqrt(0.299) = 0.547 rather than 0.299. (2) The table reports sigma_1
    # on the same convention, and sigma_1 read as a variance is what
    # reconciles Figure 2A with Figure 2B (see the residual-error block).
    etalcl ~ 0.299  # Wang 2019 Table III: 'omega CL' 0.299 (variance; 95% CI 0.212-0.359) -> sqrt(exp(0.299)-1) = 59.0% CV
    etalvc ~ 0.299  # Wang 2019 Table III: 'omega V' 0.299 (variance; 95% CI 0.211-0.356) -> sqrt(exp(0.299)-1) = 59.0% CV

    # ----- Residual unexplained variability -----
    # Wang 2019 Equation (b) is OB = IP x (1 + eps_1), i.e. purely
    # proportional, and Table III labels sigma_1 'residual variability
    # (proportional error)' with an estimate of 1.020. Read as a NONMEM
    # SIGMA variance this is a proportional standard deviation of
    # sqrt(1.020) = 1.010, i.e. about 101% CV -- extreme, but it is what the
    # published diagnostics show. Verified against Wang 2019 Figure 2:
    #  - Figure 2A (observations vs individual predictions) has observed
    #    concentrations ranging from roughly 0.3x to 4.5x the individual
    #    prediction at any given prediction level, and a lowess curve that
    #    falls far below the line of identity above IPRED 10.
    #  - Figure 2B (|iWRES| vs individual predictions) tops out at 4.1 for
    #    the bulk of the data with one outlier at 9.3. For a proportional
    #    model |iWRES| = |DV/IPRED - 1| / sd, so the 4.5x ratio visible in
    #    panel A implies sd near 3.5/4.1 = 0.85-1.0. A sd of 0.10 would put
    #    those same points at |iWRES| near 35, which panel B rules out.
    #  - The dense row of points at DV = 0 in panel A pairs with the dense
    #    row at |iWRES| = 0 in panel B; those are dosing records, for which
    #    NONMEM writes DV = 0 and IWRES = 0, not below-quantification
    #    observations, and they do not inform sigma.
    # The magnitude is plausible for opportunistic therapeutic-drug-
    # monitoring sampling where the recorded sampling times are unreliable:
    # with a 12.3 h half-life and 12 h dosing interval, a mis-recorded time
    # maps onto a large concentration error.
    propSd <- 1.010
    label("Proportional residual error (fraction)")  # Wang 2019 Table III: sigma 1 = 1.020 (variance; 95% CI 0.864-1.179) -> sd = sqrt(1.020) = 1.010
  })

  model({
    # ----- 1. Individual parameters -----
    # Wang 2019 Equations (e) and (f). The paper writes the exponential
    # inter-individual variability term of Equation (a) as a multiplicative
    # 'omega' factor in Equations (e) and (f); it is exp(eta) here.
    cl <- exp(lcl + etalcl) * (PLT / 200)^e_plt_cl
    vc <- exp(lvc + etalvc)

    # ----- 2. Micro-constants -----
    kel <- cl / vc

    # ----- 3. ODE system -----
    # Intravenous administration directly into the central compartment; the
    # paper reports no absorption phase and no infusion duration.
    d/dt(central) <- -kel * central

    # ----- 4. Observation and error -----
    # central is mg and vc is L, so central/vc is mg/L, matching the
    # 'Observations (mg/l)' axis of Wang 2019 Figure 3.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
