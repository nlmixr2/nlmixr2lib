Juul_2019_somapacitan <- function() {
  description <- "One-compartment population PK model with parallel first-order and zero-order (via one transit compartment) subcutaneous absorption and Michaelis-Menten elimination, linked to an indirect-response IGF-I model with an additive Emax stimulation of IGF-I production, for once-weekly somapacitan in healthy adults, adults with growth hormone deficiency and children with growth hormone deficiency (Juul 2019)"
  reference <- "Juul RV, Rasmussen MH, Agerso H, Overgaard RV. Pharmacokinetics and Pharmacodynamics of Once-Weekly Somapacitan in Children and Adults: Supporting Dosing Rationales with a Model-Based Analysis of Three Phase I Trials. Clin Pharmacokinet. 2019;58(1):63-75. doi:10.1007/s40262-018-0662-5"
  vignette <- "Juul_2019_somapacitan"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # What each ODE state holds, verified against Juul 2019 Figure 1: the dose
  # enters the absorption compartment 'abs' (depot), leaves it first-order
  # (Ka) to central and zero-order (K0/F) to a single transit compartment,
  # which empties first-order (Ktr) into central; Michaelis-Menten elimination
  # from central; the central concentration stimulates the production rate of
  # the 'IGF-I' compartment. The IGF-I state is carried as a serum
  # concentration (Kin in ng/mL/h, Kout in 1/h; Table 4).
  compartmentData <- list(
    depot = list(analyte = "somapacitan", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "somapacitan", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "somapacitan", units = "mg", specimen = "serum", verified = TRUE),
    igf1 = list(analyte = "insulin-like growth factor I", units = "ng/mL", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate with reference weight 85 kg (Methods, Covariate",
        "Analysis, equation 1: P_i = P_typ * (BW / 85 kg)^theta * exp(eta)).",
        "One exponent is shared by K0/F, V/F and Vmax/F (Table 3 footnote b),",
        "another by Ka and Ktr (footnote c); Kout and Emax carry their own",
        "exponents (Table 4). Baseline body weight; the trials lasted at",
        "most 4 weeks."
      ),
      source_name = "BW"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator (1 = healthy adult, 0 = subject with growth hormone deficiency)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (adult or child with growth hormone deficiency)",
      notes = paste(
        "Table 4 'theta HV' rows; HV = healthy volunteer (vs. GHD) per the",
        "Table 4 legend. Discrete covariates enter as",
        "P_i = P_typ * exp(theta * Cov) * exp(eta) (Methods, equation 2).",
        "Shifts Kin (theta 0.550) and Kout (theta -0.241), i.e. mainly the",
        "baseline IGF-I level. Equals 1 for every subject of trial",
        "NCT01514500 and 0 for trials NCT01706783 and NCT01973244."
      ),
      source_name = "HV"
    ),
    CHILD = list(
      description = "Child age-cohort indicator (1 = child, 0 = adult)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (adult)",
      notes = paste(
        "Selects the residual error of the log-transformed somapacitan",
        "concentrations: 'Separate residual distributions were fitted to",
        "adults and children' (Methods, Data Handling; Table 3",
        "'Experimental error adults' / 'Experimental error children').",
        "Children were 6-11 years old (Table 2). No structural parameter",
        "depends on it -- the child/adult differences are carried by WT."
      ),
      source_name = "age group"
    ),
    OCC = list(
      description = "Occasion index for the inter-occasion variability of the PK model",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Consumed by the ini() nesting syntax 'etaiov_lvc ~ ... | OCC' (Table 3",
        "IOV column: one IOV shared by K0/F, V1/F and Vmax/F), not by a",
        "model() expression. The paper does not define the occasion; one",
        "occasion per weekly dose (dosing interval) is assumed. An event table",
        "must carry at least two distinct OCC values for rxode2 to draw the",
        "occasion effect."
      ),
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 123,
    n_studies = 3,
    n_observations = "5171 somapacitan concentrations above the LLOQ (0.5 ng/mL) plus serial IGF-I profiles",
    age_range = "6-69 years (children 6-11 years; adults 21-69 years)",
    age_mean = "32.6 years (healthy adults 33.9; adults with GHD 51.4; children 8.3)",
    weight_range = "18.5-120.5 kg (children 18.5-41.2 kg; adults 54.1-120.5 kg)",
    weight_mean = "66.9 kg (healthy adults 74.7; adults with GHD 82.5; children 26.2)",
    sex_female_pct = 13.0,
    race_ethnicity = "Healthy-adult trial: 24 Japanese and 49 non-Japanese (Caucasian); other trials: no Japanese subjects (Table 2)",
    disease_state = "Healthy adults (n = 73), adults with growth hormone deficiency (n = 26) and prepubertal children with growth hormone deficiency (n = 24)",
    dose_range = paste(
      "Once-weekly subcutaneous somapacitan: healthy adults 0.02-0.32 mg/kg",
      "single dose and 0.02-0.24 mg/kg x 4 weekly doses; adults with GHD",
      "0.02-0.12 mg/kg x 4 weekly doses; children with GHD 0.02-0.16 mg/kg",
      "single dose (Table 1). The 0.01 mg/kg healthy-adult cohort was excluded."
    ),
    regions = "Not reported in the article",
    notes = paste(
      "Pooled phase I trials NCT01514500 (healthy adults), NCT01706783",
      "(adults with GHD) and NCT01973244 (children with GHD); Tables 1 and 2.",
      "107 males and 16 females. Baseline IGF-I mean (SD): healthy adults",
      "190.7 (43.2), adults with GHD 93.4 (40.7), children 107.9 (90.1) ng/mL."
    )
  )

  ini({
    # -------------------------------------------------------------------
    # PK (Table 3, final PK model). Typical values for an 85 kg subject
    # with adult GHD. All apparent (/F) quantities; dose in mg.
    # -------------------------------------------------------------------
    lka <- log(0.0253); label("First-order absorption rate constant Ka, 85 kg (1/h)")              # Table 3 'Ka (1/h)' = 0.0253 [0.0223-0.0284]
    lr1 <- log(0.211); label("Zero-order absorption rate K0/F from depot to transit, 85 kg (mg/h)") # Table 3 'K0/F (mg/h)' = 0.211 [0.189-0.233]
    lktr <- log(0.0102); label("Transit rate constant Ktr, 85 kg (1/h)")                             # Table 3 'Ktr (1/h)' = 0.0102 [0.00932-0.0111]
    lvc <- log(12.5); label("Apparent volume of distribution V1/F, 85 kg (L)")                      # Table 3 'V1/F (L)' = 12.5 [11.0-13.9]
    lvmax <- log(0.268); label("Apparent maximum elimination rate Vmax/F, 85 kg (mg/h)")            # Table 3 'Vmax/F (mg/h)' = 0.268 [0.238-0.298]
    lkm <- log(63.0); label("Michaelis constant for saturable elimination Km (ng/mL)")               # Table 3 'Km (ng/mL)' = 63.0 [59.2-66.9]

    e_wt_vc <- 1.28; label("Power exponent of (WT/85) shared by K0/F, V1/F and Vmax/F (unitless)")   # Table 3 'theta BW K0, theta BW V1, theta BW Vmax' = 1.28 [1.13-1.44]
    e_wt_ka <- -0.503; label("Power exponent of (WT/85) shared by Ka and Ktr (unitless)")            # Table 3 'theta BW Ka, theta BW Ktr' = -0.503 [-0.608 to -0.398]

    # -------------------------------------------------------------------
    # PK/PD (Table 4, final PK/PD model): indirect response with an
    # additive Emax stimulation of the IGF-I production rate.
    # -------------------------------------------------------------------
    lkout <- log(0.0252); label("First-order IGF-I turnover rate constant Kout, 85 kg GHD (1/h)")   # Table 4 'Kout (1/h)' = 0.0252 [0.0214-0.0290]
    lkin <- log(2.19); label("Zero-order IGF-I production rate Kin, GHD (ng/mL/h)")                  # Table 4 'Kin (ng/mL/h)' = 2.19 [1.79-2.59]
    lec50 <- log(16.8); label("Somapacitan concentration at half-maximal stimulation EC50 (ng/mL)")  # Table 4 'EC50 (ng/mL)' = 16.8 [14.7-18.9]
    lemax <- log(15.1); label("Maximum increase in IGF-I production rate Emax, 85 kg (ng/mL/h)")    # Table 4 'Emax (ng/mL/h)' = 15.1 [14.2-16.0]

    e_wt_kout <- -0.244; label("Power exponent of (WT/85) on Kout (unitless)")                       # Table 4 'theta BW Kout' = -0.244 [-0.403 to -0.0861]
    e_wt_emax <- 0.460; label("Power exponent of (WT/85) on Emax (unitless)")                        # Table 4 'theta BW Emax' = 0.460 [0.332-0.587]
    e_healthy_kout <- -0.241; label("Healthy-volunteer exponential coefficient on Kout (unitless)")   # Table 4 'theta HV Kout' = -0.241 [-0.377 to -0.106]
    e_healthy_kin <- 0.550; label("Healthy-volunteer exponential coefficient on Kin (unitless)")      # Table 4 'theta HV Kin' = 0.550 [0.363-0.737]

    # -------------------------------------------------------------------
    # IIV and IOV. Tables 3 and 4 report CV%; converted to log-normal
    # variances with omega^2 = log(CV^2 + 1). One IIV and one IOV are shared
    # by K0/F, V1/F and Vmax/F (Table 3 footnote a), i.e. a random effect on
    # 1/F; they carry the V1/F name. EC50 and Km have no IIV.
    # -------------------------------------------------------------------
    etalka ~ log(1 + 0.477^2)          # Table 3 Ka IIV CV 47.7%
    etalvc ~ log(1 + 0.318^2)          # Table 3 K0/F, V1/F, Vmax/F shared IIV CV 31.8%
    etalktr ~ log(1 + 0.328^2)         # Table 3 Ktr IIV CV 32.8%
    etalkout ~ log(1 + 0.221^2)        # Table 4 Kout IIV CV 22.1%
    etalkin ~ log(1 + 0.369^2)         # Table 4 Kin IIV CV 36.9%
    etalemax ~ log(1 + 0.207^2)        # Table 4 Emax IIV CV 20.7%
    # Written as a literal: with a '| OCC' condition, lotri 1.0.5 drops an
    # enclosing log() and would store 1 + 0.216^2 as the variance.
    etaiov_lvc ~ 0.04560 | OCC # Table 3 K0/F, V1/F, Vmax/F shared IOV CV 21.6%: log(1 + 0.216^2) = 0.04560

    # -------------------------------------------------------------------
    # Residual error. Somapacitan: log-transformed data with additive error
    # (Methods, Data Handling), separate for adults and children; the
    # Table 3 'Experimental error' values are read as SDs on the log scale,
    # as Table 4 prints its IGF-I errors as SDs (ng/mL, %).
    # IGF-I: combined additive and proportional.
    # -------------------------------------------------------------------
    expSd_adult <- 0.377; label("Log-scale residual SD of somapacitan, adults")                  # Table 3 'Experimental error adults' = 0.377
    expSd_ped <- 0.425; label("Log-scale residual SD of somapacitan, children")                  # Table 3 'Experimental error children' = 0.425
    addSd_IGF1 <- 7.32; label("Additive residual SD of IGF-I (ng/mL)")                           # Table 4 'Additive error (ng/mL)' = 7.32
    propSd_IGF1 <- 0.143; label("Proportional residual SD of IGF-I (fraction)")                  # Table 4 'Proportional error (%)' = 14.3
  })

  model({
    # Individual PK parameters (Methods, equation 1, reference weight 85 kg).
    ka <- exp(lka + etalka) * (WT / 85)^e_wt_ka
    ktr <- exp(lktr + etalktr) * (WT / 85)^e_wt_ka
    r1 <- exp(lr1 + etalvc + etaiov_lvc) * (WT / 85)^e_wt_vc
    vc <- exp(lvc + etalvc + etaiov_lvc) * (WT / 85)^e_wt_vc
    vmax <- exp(lvmax + etalvc + etaiov_lvc) * (WT / 85)^e_wt_vc
    km <- exp(lkm)

    # Serum concentration: mg / L * 1000 = ng/mL (Km and EC50 are in ng/mL).
    Cc <- 1000 * central / vc

    # Zero-order K0/F transfer from depot to transit1 (Figure 1), which runs
    # at the constant rate r1 while drug remains in the depot and must stop
    # when it is exhausted. The gate depot / (depot + 1e-5 mg) is 1 to within
    # 1e-5 / depot and closes smoothly once the depot holds less than about
    # 10 ng, so the depot cannot be driven negative.
    zo_flux <- r1 * depot / (depot + 1e-5)

    d/dt(depot) <- -ka * depot - zo_flux
    d/dt(transit1) <- zo_flux - ktr * transit1
    d/dt(central) <- ka * depot + ktr * transit1 - vmax * Cc / (km + Cc)

    # IGF-I indirect response: additive Emax stimulation of the production
    # rate by the somapacitan concentration (Figure 1; Emax carries
    # production-rate units, ng/mL/h). Discrete covariates enter as
    # exp(theta * Cov) (Methods, equation 2). The IGF-I state starts at the
    # drug-free steady state Kin / Kout.
    kin <- exp(lkin + etalkin + e_healthy_kin * DIS_HEALTHY)
    kout <- exp(lkout + etalkout + e_healthy_kout * DIS_HEALTHY) * (WT / 85)^e_wt_kout
    emax <- exp(lemax + etalemax) * (WT / 85)^e_wt_emax
    ec50 <- exp(lec50)

    d/dt(igf1) <- kin + emax * Cc / (ec50 + Cc) - kout * igf1
    igf1(0) <- kin / kout

    IGF1 <- igf1

    expSd <- expSd_adult * (1 - CHILD) + expSd_ped * CHILD

    Cc ~ lnorm(expSd)
    IGF1 ~ add(addSd_IGF1) + prop(propSd_IGF1)
  })
}
