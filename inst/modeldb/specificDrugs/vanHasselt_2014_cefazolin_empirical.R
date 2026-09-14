vanHasselt_2014_cefazolin_empirical <- function() {
  description <- "Two-compartment population PK model for free and total cefazolin in pregnant women (empirical gestational covariate model). Clearance is the sum of a constant arm and an arm scaled linearly by gestational age; disposition parameters are referenced to the UNBOUND concentration and total cefazolin is reconstructed algebraically as Cunbound/fu."
  reference <- paste(
    "van Hasselt JGC, Allegaert K, van Calsteren K, Beijnen JH, Schellens JHM, Huitema ADR.",
    "Semiphysiological versus empirical modelling of the population pharmacokinetics of free and total cefazolin during pregnancy.",
    "Biomed Res Int. 2014;2014:897216. doi:10.1155/2014/897216.",
    "Corrigendum: Biomed Res Int. 2015;2015:124035. doi:10.1155/2015/124035",
    "(corrects the Table 2 covariate-equation footnotes, which are switched in the original,",
    "and clarifies that CL, Vc, Vp and Q were fitted on free cefazolin so CLtotal = CLfree * fu).",
    sep = " "
  )
  vignette <- "vanHasselt_2014_cefazolin"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    EGA = list(
      description        = "Maternal estimated gestational age at the time of the observation",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives the empirical gestational effect on clearance via the linear term (1 + EGA/40);",
        "the normalisation factor 40 is the maximum gestational age in the pooled dataset",
        "(paper Methods 2.3: 'For linear GA models this was the maximum GA of 40').",
        "Observed range 17-40 weeks, median 33 weeks (Table 1).",
        "NOTE: unlike the semiphysiological sibling model, this empirical form does NOT collapse",
        "to a non-pregnant clearance at EGA = 0 - the gestational arm contributes cl_renal * 1 there.",
        "It is a within-window linear interpolation and should not be extrapolated below 17 weeks;",
        "the paper makes exactly this point in its Discussion when contrasting the two approaches.",
        "A fixed EGA of 40 weeks was assigned to the term-pregnancy caesarean cohort (Methods 2.1)."
      ),
      source_name        = "GA"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "cefazolin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefazolin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 94L,
    n_studies      = 3L,
    n_observations = 187L,
    age_range      = "20-42 years",
    age_median     = "31 years",
    weight_range   = "54-99 kg",
    weight_median  = "72 kg",
    sex_female_pct = 100,
    race_ethnicity = "not reported in the source paper",
    disease_state  = "pregnant women undergoing in utero surgical intervention, elective caesarean delivery, or fetal intervention; cefazolin given as surgical prophylaxis",
    dose_range     = "1 g or 2 g intravenously; 2 g every 8 h for 2 days in the prospective cohort, single 1 g or 2 g bolus in the two literature cohorts",
    ga_range       = "17-40 weeks (median 33)",
    regions        = "Belgium (University Hospitals Leuven) plus two previously published cohorts",
    renal_function = "serum creatinine median 0.64 mg/dL (range 0.33-0.88); creatinine clearance computed by Cockcroft-Gault using body weight",
    notes          = paste(
      "Pooled from one prospective study and two published studies with individual-level data (Methods 2.1).",
      "Prospective cohort: 41 pregnant women, 153 cefazolin observations, median GA 25 weeks (range 17-34),",
      "2 g every 8 h for 2 days during in utero surgery, free cefazolin available for 84% of observations.",
      "Fiore Mitchell et al.: 24 plasma samples at term caesarean delivery after 1 g i.v. bolus, mean sampling time 1.85 h,",
      "GA fixed at 40 weeks for all patients in that study.",
      "Brown et al.: 10 fetal interventions in 7 women, single 2 g i.v. bolus, mean GA 27 weeks, mean sampling time 0.5 h.",
      "Demographics in Table 1. Substantial missing data imputed by pooled medians (Methods 2.5):",
      "body weight missing for 46% of patients, age for 55%, and some serum creatinine for 58%."
    )
  )

  ini({
    # Structural parameters. Source: van Hasselt 2014 Table 2, column
    # 'Empirical model CL ~ GA'. CL and Q are reported in L/min and the volumes
    # in L, so the model time unit is minutes and the values are carried
    # unconverted.
    #
    # CLEARANCE DECOMPOSITION. The paper writes clearance as the sum of two
    # arms (corrigendum footnote a): CL = theta_CL0 + theta_CLPreg * (1 + GA/40).
    # Methods 2.4.2 names the same theta pair in the sibling semiphysiological
    # model 'non-CrCL-related clearance' and 'GFR-related clearance', which is
    # why the registered renal / non-renal clearance-arm canonicals are used
    # here. Neither arm is the total clearance: the total is their sum.
    lcl_nonren <- log(0.119); label("Non-gestational (non-renal) clearance arm, referenced to unbound cefazolin (L/min)") # Table 2, row 'Clearance' theta_CL0, empirical column = 0.119 L/min (RSE 58%)
    lcl_renal  <- log(0.217); label("Gestation-scaled (renal) clearance arm at GA = 0, referenced to unbound cefazolin (L/min)") # Table 2, row 'Gestation effect on clearance' theta_CLPreg, empirical column = 0.217 (RSE 16%)
    lvc        <- log(33.1);  label("Central volume, referenced to unbound cefazolin (L)")                # Table 2, row 'Central volume' V_C, empirical column = 33.1 L (RSE 17%)
    lvp        <- log(12.8);  label("Peripheral volume, referenced to unbound cefazolin (L)")             # Table 2, row 'Peripheral volume' V_P, empirical column = 12.8 L (RSE 27%)
    lq         <- log(0.326); label("Intercompartmental clearance, referenced to unbound cefazolin (L/min)") # Table 2, row 'Intercompartmental clearance' Q, empirical column = 0.326 L/min (RSE 25%)

    # Protein binding. Results 3.1 eq 7: C_total = C_free / fu, i.e. a constant
    # (linear) binding model; nonlinear binding could not be identified.
    lfu        <- log(0.286); label("Fraction of cefazolin unbound in plasma (unitless)")                 # Table 2, row 'Free fraction' F_U, empirical column = 0.286 (RSE 5%)

    # IIV. Methods eq 1 is exponential: P_i = P * exp(eta_i). Table 2 reports the
    # between-subject variability as CV%, so the internal variance is taken as
    # (CV/100)^2. See the vignette Errata for the alternative exact-lognormal
    # reading log(1 + CV^2), which differs by <1% for these magnitudes.
    etalcl ~ 0.039601  # Table 2, row 'Clearance' omega_CL, empirical column = 19.9 CV% -> 0.199^2
    etalvc ~ 0.226576  # Table 2, row 'Central volume' omega_V1, empirical column = 47.6 CV% -> 0.476^2
    etalvp ~ 0.116964  # Table 2, row 'Peripheral volume' omega_V2, empirical column = 34.2 CV% -> 0.342^2
    etalfu ~ 0.030625  # Table 2, row 'Free fraction' omega_FU, empirical column = 17.5 CV% -> 0.175^2

    # Residual error. Methods eq 2 is a combined proportional + additive model on
    # the linear concentration scale, fitted separately for the free and total
    # cefazolin outputs. Table 2 heads this block 'Residual unexplained
    # variability variances', so each tabulated value is a VARIANCE and the SD
    # entered here is its square root.
    propSd          <- 0.181108; label("Proportional residual error, total cefazolin (fraction)")   # Table 2, row 'Proportional, total concentration' sigma_TP, empirical column = 0.0328 (variance) -> sqrt = 0.181108
    addSd           <- 0.917606; label("Additive residual error, total cefazolin (mg/L)")           # Table 2, row 'Additive, total concentration' sigma_TA, empirical column = 0.842 (variance) -> sqrt = 0.917606
    propSd_Cunbound <- 0.125698; label("Proportional residual error, free cefazolin (fraction)")    # Table 2, row 'Proportional, free concentration' sigma_FP, empirical column = 0.0158 (variance) -> sqrt = 0.125698
    addSd_Cunbound  <- 0.478539; label("Additive residual error, free cefazolin (mg/L)")            # Table 2, row 'Additive, free concentration' sigma_FA, empirical column = 0.229 (variance) -> sqrt = 0.478539
  })

  model({
    # 1. Empirical gestational effect on clearance.
    #    Corrigendum (Biomed Res Int 2015;2015:124035) Table 2 footnote a:
    #      CL = theta_CL0 + theta_CLPreg * (1 + (GA/40))
    #    The ORIGINAL paper prints this equation against the semiphysiological
    #    column; the corrigendum states footnotes a and b were switched, and the
    #    corrected assignment is the one that matches Methods eq 3 (a linear GA
    #    relation normalised by the maximum GA of 40) and reproduces the base
    #    model's CL of 0.49 L/min at the cohort median GA of 33 weeks.
    preg_cl <- 1 + EGA / 40

    # 2. Individual parameters. Methods eq 1 places the exponential IIV on the
    #    covariate-adjusted typical value, so the single reported omega_CL
    #    multiplies the SUM of the two clearance arms rather than either arm.
    cl_nonren <- exp(lcl_nonren)
    cl_renal  <- exp(lcl_renal)
    cl <- (cl_nonren + cl_renal * preg_cl) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)
    fu <- exp(lfu + etalfu)

    # 3. Micro-constants. cl, q, vc and vp are all referenced to the unbound
    #    concentration, so the compartment amounts are TOTAL cefazolin and the
    #    micro-constants are formed in the usual way.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. Two-compartment disposition (Results 3.1: 'A two-compartmental model
    #    best described the data').
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-          k12 * central - k21 * peripheral1

    # 5. Observations. The model was fitted on the FREE cefazolin concentration
    #    (corrigendum clarification), so central/vc is the unbound concentration
    #    and the total concentration follows from the constant binding model of
    #    Results eq 7, C_total = C_free / fu. Equivalently CLtotal = CLfree * fu,
    #    as the corrigendum states.
    Cunbound <- central / vc
    Cc       <- Cunbound / fu

    Cc       ~ add(addSd) + prop(propSd)
    Cunbound ~ add(addSd_Cunbound) + prop(propSd_Cunbound)
  })
}
