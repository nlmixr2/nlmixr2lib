Giacometti_2025_dalbavancin_nlls <- function() {
  description <- paste(
    "Two-compartment intravenous PK model for dalbavancin, estimated by NAIVE-POOLED NONLINEAR",
    "LEAST SQUARES in 218 adults undergoing therapeutic drug monitoring during long-acting",
    "dalbavancin therapy. This is the least-squares arm of a three-way methodological comparison",
    "(naive-pooled least squares vs NLME vs a Neural ODE) fitted to the same cohort as",
    "modellib('Giacometti_2025_dalbavancin_nlme'); it is structurally identical and differs only in",
    "its parameter values. IT IS A NEGATIVE COMPARATOR, NOT A RECOMMENDED MODEL: the source reports",
    "that it 'systematically underestimates long-term drug concentrations' and 'fails to accurately",
    "estimate the inter-compartmental clearance Q' (Q is 15-fold the NLME estimate and 11-fold the",
    "published literature value), giving a terminal half-life near 8 days against roughly 28 days",
    "for the NLME fit. It is packaged so that the paper's central comparison is reproducible;",
    "prefer the NLME sibling, or modellib('Cojutti_2024_dalbavancin'), for any predictive use. This",
    "same fit was also used per fold to generate the synthetic profiles that pre-trained the Neural",
    "ODE. No covariate is carried (Appendix B), and no inter-individual variability exists to carry:",
    "a naive-pooled least-squares fit has no random-effects layer at all. Residual error is not",
    "reported, so simulations are typical-value only.",
    sep = " "
  )
  reference <- paste(
    "Giacometti T, Rocchi E, Cojutti PG, Magnani F, Remondini D, Pea F, Castellani G.",
    "Leveraging Neural ODEs for Population Pharmacokinetics of Dalbavancin in Sparse Clinical",
    "Data. Entropy. 2025;27(6):602. doi:10.3390/e27060602",
    sep = " "
  )
  vignette <- "Giacometti_2025_dalbavancin"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Giacometti 2025 Figure 2 is the schematic.
  compartmentData <- list(
    central     = list(analyte = "dalbavancin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dalbavancin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # No covariate enters model(). See covariatesDataExcluded below.
  covariateData <- NULL

  # Same screen, same outcome, as the NLME sibling. Giacometti 2025 Methods
  # 2.2: "The covariates were not included in the NLME and two-compartment
  # models because, when tested, their inclusion did not improve the results;
  # see Appendix B." Documentation only -- none of these is referenced in
  # model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened and not retained (Methods 2.2, Appendix B). Table 1: range 18-92 y, mean 64 +/- 16 y."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened and not retained (Methods 2.2, Appendix B). Table 1: range 145-190 cm, mean 171 +/- 9 cm."
    ),
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened and not retained (Methods 2.2, Appendix B), so this model applies NO allometric scaling. Table 1: range 40-140 kg, mean 77 +/- 16 kg."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened and not retained (Methods 2.2, Appendix B). Methods 2.1: 145 males and 73 females of 218 patients (33.5% female)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened and not retained (Methods 2.2, Appendix B). Table 1 reports it in SI units as 'Creatine conc. (umol/L)', range 19-411, mean 94 +/- 48; creatinINE is meant (94 umol/L = 1.06 mg/dL). No eGFR or creatinine clearance is derived anywhere in the paper."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 218L,
    n_studies      = 1L,
    n_observations = "669 dalbavancin plasma concentration measurements over 703 recorded administrations (Methods 2.1)",
    age_range      = "18-92 years",
    age_mean       = "64 +/- 16 years (Table 1)",
    height_range   = "145-190 cm (mean 171 +/- 9 cm, Table 1)",
    weight_range   = "40-140 kg",
    weight_mean    = "77 +/- 16 kg (Table 1)",
    sex_female_pct = 33.5,
    race_ethnicity = "Not reported; single-centre Italian cohort.",
    disease_state  = paste(
      "Adults receiving dalbavancin under therapeutic drug monitoring for skin and soft tissue",
      "infections and, as a second-line agent, for staphylococcal bone and joint infections,",
      "vascular prosthetic joint infections and endocarditis (Methods 2.1). The paper does not",
      "tabulate the per-infection-type counts."
    ),
    renal_function = "Reported only as serum creatinine: range 19-411 umol/L, mean 94 +/- 48 umol/L (Table 1). No renal covariate is carried by the model.",
    dose_range     = "350-1500 mg intravenously, 703 administrations across the 218 patients (Methods 2.1). The paper does not report the dosing intervals or the infusion duration.",
    regions        = "Italy (IRCCS Azienda Ospedaliero-Universitaria di Bologna)",
    notes          = paste(
      "Same retrospective single-centre therapeutic-drug-monitoring cohort as the NLME sibling",
      "modellib('Giacometti_2025_dalbavancin_nlme'), April 2021 to December 2024 (Ethics Committee",
      "897/2021/Oss/AOUBo). ESTIMATION IS NAIVE-POOLED NONLINEAR LEAST SQUARES, not a mixed-effects",
      "fit: Methods 2.2 states 'the fitting procedures of the two compartment models have been",
      "performed via nonlinear least squares (NLLS)', so there is no random-effects layer and no",
      "individual-level parameter distribution -- the four values are a single pooled",
      "concentration-time curve through all 669 measurements. As with the NLME arm, the reported",
      "values are the MEAN over six cross-validation folds and the plus-or-minus figures in Table 4",
      "are the across-fold standard deviation. The same per-fold fit was reused as the data",
      "augmentation engine for the Neural ODE: Methods 2.2, 'during each fold of the",
      "cross-validation, a two-compartment pharmacokinetic model was fitted to the available data",
      "using nonlinear least squares estimation. The fitted parameters were then used to simulate",
      "drug concentration-time profiles under various dosing regimens. From each simulated profile,",
      "50 time points were sampled on a logarithmic time scale to generate synthetic training",
      "data.' The paper's verdict on this arm is explicit -- Discussion: 'The two-compartment model",
      "systematically underestimates long-term drug concentrations, as evidenced by the obvious",
      "trend in relative residuals (Figure 4)' and 'the two-compartment model fails to accurately",
      "estimate the inter-compartmental clearance Q' -- and the Kolmogorov-Smirnov tests of Table 2",
      "put its relative-residual distribution significantly apart from both the NLME",
      "(p = 1.4e-19) and the Neural ODE (p = 2.5e-18)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters -- Giacometti 2025 Table 4,
    # "Two-Compartment" column (the naive-pooled NLLS fit). Table 4
    # caption: "Mean values of the estimated parameters for the
    # two-compartment model and NLME among the six iterations of the
    # cross-validation. The associated error is the standard deviation."
    # Each value below is the across-fold MEAN and the quoted spread is
    # the across-fold SD -- not a standard error and not a variance
    # component. Reported on the linear scale and log-transformed here per
    # the nlmixr2lib convention. Units from Table 4's row labels: Cl and Q
    # in L/h, V1 and V2 in L.
    #
    # NOTE how far Q sits from every other estimate of the same quantity:
    # 0.42 L/h here against 0.028 L/h for the NLME arm on the identical
    # data and 0.038 L/h in the literature column (Cojutti 2024). The
    # paper names this as the fit's failure mode rather than a finding.
    # ------------------------------------------------------------------
    lcl <- log(0.054); label("Clearance (L/h)")                        # Giacometti 2025 Table 4, Two-Compartment column: Cl = 0.054 +/- 0.002 L/h (across-fold SD)
    lvc <- log(5.3);   label("Central volume of distribution V1 (L)")      # Giacometti 2025 Table 4, Two-Compartment column: V1 = 5.3 +/- 0.4 L (across-fold SD)
    lq  <- log(0.42);  label("Intercompartmental clearance (L/h)")      # Giacometti 2025 Table 4, Two-Compartment column: Q  = 0.42 +/- 0.14 L/h (across-fold SD); the source calls this estimate inaccurate
    lvp <- log(8.8);   label("Peripheral volume of distribution V2 (L)")   # Giacometti 2025 Table 4, Two-Compartment column: V2 = 8.8 +/- 1.5 L (across-fold SD)

    # ------------------------------------------------------------------
    # Inter-individual variability: STRUCTURALLY ABSENT, not merely
    # unreported. Methods 2.2 fits this arm by naive-pooled nonlinear
    # least squares, which has no random-effects layer -- one curve is
    # fitted through the pooled observations of all 218 patients. There is
    # therefore no omega to transcribe and none is declared. (The NLME
    # sibling is the case where random effects genuinely exist but go
    # unreported.)
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Residual unexplained variability: NOT REPORTED. Methods 2.2 gives
    # the objective as a plain squared loss on observed minus predicted
    # concentration, but no residual SD, %CV or sigma appears anywhere in
    # the article. Encoded as exact zeros rather than invented, following
    # inst/modeldb/specificDrugs/Takada_2025_vancomycin.R; both a
    # proportional and an additive term are carried per the standing
    # convention for an unspecified residual form.
    # ------------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")  # Giacometti 2025: no residual-error estimate is published
    addSd  <- fixed(0); label("Additive residual SD (mg/L; 0 -- not reported in the source)")          # Giacometti 2025: no residual-error estimate is published
  })

  model({
    # 1. PK parameters. No covariate (Appendix B) and no random effect (a
    #    naive-pooled least-squares fit has none).
    cl <- exp(lcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # 2. Two-compartment IV disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 3. ODE system, Giacometti 2025 Equation 2a -- the canonical
    #    mass-conserving two-compartment system, identical to the NLME
    #    sibling. See the extended note in
    #    inst/modeldb/specificDrugs/Giacometti_2025_dalbavancin_nlme.R for
    #    why the printed equation is read on AMOUNTS despite the
    #    surrounding prose naming the states as concentrations. Dalbavancin
    #    is intravenous, so doses enter `central` directly; the infusion is
    #    set in the event table, not by any model-side term.
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-    k12       * central - k21 * peripheral1

    # 4. Total plasma dalbavancin concentration in mg/L (dose in mg,
    #    volumes in L).
    Cc <- central / vc

    # 5. Observation. Both residual magnitudes are zero (unreported), so
    #    this reduces to a deterministic prediction.
    Cc ~ prop(propSd) + add(addSd)
  })
}
