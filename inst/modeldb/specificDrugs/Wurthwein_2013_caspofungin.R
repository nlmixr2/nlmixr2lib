Wurthwein_2013_caspofungin <- function() {
  description <- paste(
    "Linear two-compartment population PK model with proportional residual",
    "error for once-daily 2-hour intravenous caspofungin infusions (70, 100,",
    "150, 200 mg QD) in adults with proven or probable invasive aspergillosis",
    "(Wurthwein 2013). Clearance and central volume share a single linear",
    "body-weight fractional change centred on the cohort median body weight",
    "of 76 kg (CL_i = CL_typ * [1 + 0.0102 * (WT - 76)]; V1_i = V1_typ *",
    "[1 + 0.0102 * (WT - 76)]). Inter-individual variability is modelled",
    "exponentially on CL, V1, and V2 with an estimated CL-V1 covariance",
    "(correlation 0.802). Inter-occasion variability (16% CV) is included",
    "on CL across five sampling occasions (days 1, 4, 7, 14, 28) via the",
    "OCC covariate; downstream users who only need typical-value or IIV-only",
    "simulations can pass OCC = 0 (or any value outside 1..5) so the IOV",
    "terms zero out. Dose-level, gender, age, baseline serum bilirubin and",
    "baseline creatinine clearance were screened but not retained."
  )
  reference <- paste(
    "Wurthwein G, Cornely OA, Trame MN, Vehreschild JJ, Vehreschild MJGT,",
    "Farowski F, Muller C, Boos J, Hempel G, Hallek M, Groll AH.",
    "Population pharmacokinetics of escalating doses of caspofungin in a",
    "phase II study of patients with invasive aspergillosis.",
    "Antimicrob Agents Chemother. 2013;57(4):1664-1671.",
    "doi:10.1128/AAC.01912-12."
  )
  vignette <- "Wurthwein_2013_caspofungin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (recorded on study day 1 per Methods,",
        "'Pharmacokinetic sampling and recording of covariates').",
        "Linear fractional-change covariate on CL and V1 centred on",
        "the cohort median 76 kg (Table 1). Cohort range 43-104 kg",
        "(median 76 kg) per Table 1 'total' column. Same fractional",
        "change retained for CL and V1 (Section 'Population",
        "pharmacokinetics of CAS', third paragraph: 'linear modeling",
        "of different fractional changes for CL and V1 due to the",
        "covariate body weight compared to the model using the same",
        "fractional change for CL and V1 showed no difference')."
      ),
      source_name        = "BW"
    ),
    OCC = list(
      description        = paste(
        "Integer-valued occasion indicator for the five PK sampling",
        "occasions in Wurthwein 2013 (1 = day 1, 2 = day 4, 3 = day 7,",
        "4 = day 14, 5 = day 28)."
      ),
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Time-varying within subject; constant within an occasion.",
        "Decomposed inside model() into binary indicators oc1..oc5 that",
        "multiplex the per-occasion IOV etas on log-CL (16% CV; Table 2",
        "final model). For typical-value or IIV-only simulations pass",
        "OCC = 0 so all five binary indicators evaluate to FALSE and the",
        "IOV terms zero out. Source paper Methods,",
        "'Population pharmacokinetic analysis' paragraph 2: 'Interoccasion",
        "variability (IOV) was tested by exponential models on clearance",
        "(CL) and central volume of distribution (V1) on five occasions",
        "(days 1, 4, 7, 14, and 28)'; only IOV on CL was retained in the",
        "final model (P < 0.01)."
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 46L,
    n_studies        = 1L,
    age_range        = "18-74 years (median 61)",
    age_median       = "61 years",
    weight_range     = "43-104 kg (median 76)",
    weight_median    = "76 kg",
    sex_female_pct   = 54.3,
    disease_state    = paste(
      "Adults with proven or probable invasive aspergillosis as defined",
      "by modified European Organisation for Research and Treatment of",
      "Cancer (EORTC) criteria. 27 of 46 had acute leukemia and 31 of 46",
      "were neutropenic. Excluded if serum bilirubin > 3x upper limit of",
      "normal, AST/ALT > 5x ULN, alkaline phosphatase > 5x ULN, active",
      "veno-occlusive disease, or expected survival < 5 days."
    ),
    dose_range       = paste(
      "Caspofungin 70, 100, 150, or 200 mg administered once daily as a",
      "2-hour intravenous infusion (no loading dose); maximum 28 days of",
      "treatment. Dose group sizes: 70 mg n = 9, 100 mg n = 8,",
      "150 mg n = 9, 200 mg n = 20."
    ),
    regions          = "Germany (three university hospitals; September 2006 - July 2009)",
    notes            = paste(
      "Phase II dose-escalation study (EudraCT 2006-001936-30,",
      "ClinicalTrials.gov NCT00404092). 468 plasma samples collected on",
      "day 1 (pre-dose, 2 h [peak], 3 h, 5-7 h, 24 h [trough]) and at",
      "peak and trough on days 4, 7, 14, and 28; final analysis dataset",
      "462 samples (1 sample excluded for sampling-time uncertainty,",
      "5 implausible peaks/troughs excluded; one subject with extreme",
      "V2 = 43.8 L also excluded from the final-model dataset of 45",
      "subjects). Caspofungin quantified by LC-MS/MS over 84 ug/L to",
      "84,000 ug/L. NONMEM 7.1.0 FOCE with interaction; covariate",
      "selection forward at P = 0.05 and backward at P = 0.01."
    )
  )

  ini({
    # Structural fixed effects -- Wurthwein 2013 Table 2 final-model column.
    # Reference subject: WT = 76 kg (cohort median, used as the centring
    # weight in the linear body-weight covariate on CL and V1).
    lcl <- log(0.411); label("Clearance (L/h)")                                                # Table 2, Final model: CL = 0.411 L/h (5% RSE)
    lvc <- log(5.85);  label("Central volume of distribution V1 (L)")                          # Table 2, Final model: V1 = 5.85 L (4% RSE)
    lq  <- log(0.843); label("Intercompartmental clearance Q (L/h)")                           # Table 2, Final model: Q  = 0.843 L/h (13% RSE)
    lvp <- log(6.53);  label("Peripheral volume of distribution V2 (L)")                       # Table 2, Final model: V2 = 6.53 L (15% RSE)

    # Body-weight covariate effect on CL and V1 (shared fractional change).
    # CL_i = CL_typ * [1 + 0.0102 * (WT - 76)]; V1_i = V1_typ * [1 + 0.0102 * (WT - 76)].
    # Centred on the cohort median 76 kg per Table 1 ('total' column) and
    # the footnote 'a Assuming a median body weight of 76 kg' on Table 2.
    e_wt_cl_vc <- 0.0102; label("Linear fractional change in CL and V1 per kg (1/kg)")         # Table 2, Final model: Factor body wt on CL, V1 = 0.0102 (20% RSE)

    # Inter-individual variability: CL and V1 correlated (correlation 0.802
    # at the omega scale per Table 2 final model), V2 diagonal.
    # omega^2 conversions for exponential / log-normal IIV: omega^2 = log(1 + CV^2).
    #   CL: 28.5% -> log(1 + 0.285^2) = 0.07809
    #   V1: 28.8% -> log(1 + 0.288^2) = 0.07968
    #   V2: 66.8% -> log(1 + 0.668^2) = 0.36896
    # Covariance(CL,V1) = corr * sqrt(omega_CL^2 * omega_V1^2) = 0.802 * sqrt(0.07809 * 0.07968) = 0.06327.
    etalcl + etalvc ~ c(0.07809,
                        0.06327, 0.07968)                                                      # Table 2, Final model: IIV CL 28.5%, IIV V1 28.8%, correlation 0.802
    etalvp          ~ 0.36896                                                                  # Table 2, Final model: IIV V2 66.8%

    # Inter-occasion variability on CL across five sampling occasions
    # (days 1, 4, 7, 14, 28). The Wurthwein paper reports a single IOV CV
    # of 16% (Table 2 final model), implying NONMEM $OMEGA BLOCK(1) SAME
    # across occasions. Encoded here as five separate etas with the
    # variance fixed equal to the estimated occasion-1 variance after the
    # first (matching the Jonsson_2011_ethambutol pattern).
    # omega^2 = log(1 + 0.160^2) = 0.02528.
    etaiov_cl_1 ~ 0.02528                                                                       # Table 2, Final model: IOV CL = 16.0% (13% RSE); occasion-1 estimated
    etaiov_cl_2 ~ fix(0.02528)                                                                  # SAME across the five sampling occasions; no per-occasion CV reported
    etaiov_cl_3 ~ fix(0.02528)                                                                  # SAME across the five sampling occasions; no per-occasion CV reported
    etaiov_cl_4 ~ fix(0.02528)                                                                  # SAME across the five sampling occasions; no per-occasion CV reported
    etaiov_cl_5 ~ fix(0.02528)                                                                  # SAME across the five sampling occasions; no per-occasion CV reported

    # Residual error: proportional only (Table 2 final model row
    # 'Proportional residual error (%)' = 14.3% with no additive component).
    propSd <- 0.143; label("Proportional residual error (fraction)")                            # Table 2, Final model: proportional RUV = 14.3% (10% RSE)
  })

  model({
    # Decompose the integer-valued occasion column into binary indicators
    # for IOV multiplexing on log-CL. Five sampling occasions correspond to
    # days 1, 4, 7, 14, and 28 of the Wurthwein 2013 trial.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
              oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5

    # Linear body-weight effect on CL and V1 (shared fractional change
    # centred on the cohort median 76 kg). Q and V2 carry no covariate.
    wt_effect_cl_vc <- 1 + e_wt_cl_vc * (WT - 76)

    # Individual PK parameters. Exponential IIV (and IOV on CL) on the
    # log-scale with the linear weight effect applied multiplicatively.
    cl <- exp(lcl + etalcl + iov_cl) * wt_effect_cl_vc
    vc <- exp(lvc + etalvc)          * wt_effect_cl_vc
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

    # Linear two-compartment model with intravenous input to the central
    # compartment. Infusion duration (2 h per Methods, 'Study drug
    # treatment') is encoded on the dose record (rate or dur) by the user.
    d/dt(central)     <- -(cl + q) / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <-  q  / vc       * central - q / vp * peripheral1

    # Plasma concentration. Dose mg, vc L -> Cc mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
