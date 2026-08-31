Mukker_2026_tuvusertib <- function() {
  description <- paste(
    "Two-compartment population PK of the oral ATR inhibitor tuvusertib (M1774)",
    "in 55 patients with advanced solid tumors (DDRiver Solid Tumors 301 Part A1,",
    "5-270 mg QD), with first-order absorption, an absorption lag time, and",
    "concentration-dependent (greater-than-dose-proportional) apparent clearance",
    "implemented as a dimensionless 'clearance compartment' turnover state whose",
    "amount scales CL/F. Residual variability is additive on the log scale."
  )
  reference <- paste(
    "Mukker JK, Diderichsen PM, Hellmann F, Yap TA, Plummer R, Tolcher AW,",
    "de Bono JS, Gounaris I, Szucs Z, Zimmermann A, Kareva I, Bolleddula J,",
    "Seithel-Keuth A, Locatelli G, Enderlin M, Hicking C, Zutshi A, Gao W,",
    "Strotmann R, Benincosa L, Venkatakrishnan K.",
    "Integrated Population Pharmacokinetic, Pharmacodynamic, and Safety Analyses",
    "to Inform Dosage Selection in the Clinical Development of the ATR Inhibitor",
    "Tuvusertib. Clin Pharmacol Ther. 2026;119(3):618-628. doi:10.1002/cpt.70029.",
    sep = " "
  )
  vignette <- "Mukker_2026_tuvusertib"

  # The dimensionless clearance-modulating turnover state that the authors call
  # the "clearance compartment" (A4 in Figure 2a) uses the `clearance_capacity`
  # canonical, ratified with this extraction. It is deliberately NOT mapped onto
  # the `enzyme` / `enz_pool` autoinduction canonicals: those assert an enzyme
  # mechanism, and this paper makes no mechanistic claim about the state at all.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot               = list(analyte = "tuvusertib", units = "mg", specimen = "administration site", verified = TRUE),
    central             = list(analyte = "tuvusertib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1         = list(analyte = "tuvusertib", units = "mg", specimen = "tissue", verified = TRUE),
    clearance_capacity  = list(analyte = "none", units = "unitless (relative to drug-free baseline of 1)", specimen = "not applicable", verified = TRUE)
  )

  # Mukker 2026 Results, "POPPK modeling ...": "The covariate analyses did not
  # identify any discernible effect of the evaluated covariates (age, body
  # weight, sex, race, concomitant medications, laboratory measurements, and
  # organ impairment) on the variability of the POPPK parameters." No covariate
  # is retained in the final model, so none is referenced in model().
  covariateData <- list()
  covariatesDataExcluded <- list(
    AGE = list(description = "Age", units = "y", type = "continuous",
               notes = "Screened in the stepwise covariate search (Supplement S1, 'Development of the population pharmacokinetic model'); not retained in the final model."),
    WT = list(description = "Baseline body weight", units = "kg", type = "continuous",
              notes = "Screened; not retained. Continuous covariates were tested as power functions scaled to the population median or a standard value (e.g. 70 kg)."),
    SEXF = list(description = "Female sex indicator", units = "(binary)", type = "binary",
                reference_category = "male", notes = "Screened as a binary factor; not retained."),
    ECOG = list(description = "Baseline ECOG performance status", units = "(score)", type = "categorical",
                notes = "Screened; not retained."),
    RACE_ASIAN = list(description = "Asian race indicator", units = "(binary)", type = "binary",
                      reference_category = "non-Asian",
                      notes = "Screened; not retained. Ethnic sensitivity was instead assessed post hoc by overlaying observed AUCss for Asian (N = 5) and Black or African American (N = 2) patients on the 90% PI of a dose-AUC power model (Figure 4e).")
  )

  population <- list(
    species        = "human",
    n_subjects     = 55L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = c(Asian = 5, `Black or African American` = 2, `non-Asian` = 48),
    disease_state  = "advanced / metastatic solid tumors",
    dose_range     = "5-270 mg once daily, plus intermittent regimens (180 mg QD 2w on/1w off, 220 mg QD 2w on/1w off, 150 mg BID 4d on/3d off)",
    regions        = "multicenter first-in-human trial (DDRiver Solid Tumors 301, NCT04170153, Part A1)",
    notes          = paste(
      "Race counts are the group sizes reported verbatim in Mukker 2026 'Ethnic",
      "sensitivity assessment in PK and time course of HGB reduction' (Asian",
      "N = 5, Black or African American N = 2, non-Asian N = 48); they are",
      "overlapping / incomplete relative to N = 55 as published and are recorded",
      "here as reported. Age, weight and sex distributions are not tabulated in",
      "this paper -- they belong to the primary Part A1 clinical publication",
      "(reference 10). PK sampling: Cycle 1 Day 1 (pre-dose and 0.5, 1, 2, 3, 4,",
      "6, 8, 12 h), Day 2 pre-dose, Day 8 (same intensive schedule), Day 9",
      "pre-dose, Day 15 pre-dose. Estimation in NONMEM 7.3, FOCE with",
      "interaction; 250-replicate nonparametric bootstrap."
    )
  )

  ini({
    # --- Structural parameters: Mukker 2026 Table 1 (final POPPK model) -------
    # Log-transformed in ini() so re-fits stay on the positive real line; the
    # paper reports linear-scale values.
    lka   <- log(0.441); label("Absorption rate constant (KA, 1/h)")                       # Table 1: KA = 0.441 1/h (bootstrap 95% CI 0.391, 0.505)
    lcl   <- log(55.7);  label("Apparent clearance at the drug-free baseline (CL/F, L/h)") # Table 1: CL/F = 55.7 L/h (48.5, 66.1)
    lvc   <- log(30.0);  label("Apparent central volume of distribution (VC/F, L)")        # Table 1: VC/F = 30.0 L (17.3, 46.6)
    lq    <- log(3.59);  label("Apparent intercompartmental clearance (Q/F, L/h)")         # Table 1: Q/F = 3.59 L/h (2.70, 4.99)
    lvp   <- log(136);   label("Apparent peripheral volume of distribution (VP/F, L)")     # Table 1: VP/F = 136 L (85.3, 191)
    ltlag <- log(0.369); label("Absorption lag time (ALAG1, h)")                           # Table 1: ALAG1 = 0.369 h (0.346, 0.400)

    # Clearance-compartment turnover (Mukker 2026 Figure 2a; Supplement S1).
    # KCL is both the zero-order production rate INTO the clearance compartment
    # and the drug-free first-order loss rate constant out of it, so the state
    # sits at 1 in the absence of drug (see model() for the algebra).
    lkcl  <- log(0.0878); label("Clearance-compartment turnover rate constant (KCL, 1/h)") # Table 1: KCL = 0.0878 1/h (0.0682, 0.115)
    slp   <- 0.303;       label("Tuvusertib effect on clearance-compartment loss (SLP, %/(ng/mL))") # Table 1: SLP = 0.303 %/(ng/mL) (0.216, 0.407)

    # --- Interindividual variability -----------------------------------------
    # Table 1 reports IIV as %CV; omega^2 = log(CV^2 + 1) for a log-normal.
    etalcl ~ 0.05922; label("IIV on CL/F")  # Table 1: IIV CL/F = 24.7 %CV (17.2, 30.4), shrinkage 6.15%; log(0.247^2 + 1) = 0.05922
    etalvc ~ 0.86243; label("IIV on VC/F")  # Table 1: IIV VC/F = 117 %CV (93.6, 149), shrinkage 10.1%; log(1.17^2 + 1) = 0.86243

    # --- Residual unexplained variability ------------------------------------
    # Supplement S1: "the residual variability was described by an additive
    # error model on log scale", i.e. log-normal residual error in nlmixr2.
    expSd <- 0.714; label("Residual SD on the natural-log concentration scale (log ng/mL)") # Table 1: SD of RUV = 0.714 log ng/mL (0.642, 0.763), shrinkage 4.32%
  })

  model({
    # --- Individual parameters -----------------------------------------------
    ka   <- exp(lka)
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    tlag <- exp(ltlag)
    kcl  <- exp(lkcl)

    # --- Observation ----------------------------------------------------------
    # Amounts are in mg and volumes in L, so central/vc is mg/L; x1000 converts
    # to the ng/mL scale that Table 1's SLP and the pCHK1 IC90 (7.9 ng/mL) use.
    Cc <- 1000 * central / vc

    # --- Clearance compartment (Mukker 2026 Figure 2a) -----------------------
    # Figure 2a states the two rate laws verbatim:
    #     KEL = (CL / VC) * A4            (elimination from the central cmt)
    #     production into A4  = KCL
    #     loss constant of A4 = KCL * (1 + SLP * CP)
    # With drug absent the state is at KCL / KCL = 1, so CL/F equals the Table 1
    # value at low doses and A4 is dimensionless. Rising tuvusertib
    # concentrations accelerate loss of A4, A4 falls, CL/F falls, and exposure
    # rises more than dose-proportionally -- exactly the mechanism the
    # Supplement describes ("increasing tuvusertib concentrations led to a
    # reduction of the amount of the clearance compartment, which in turn
    # reduced the tuvusertib clearance").
    # SLP is tabulated in %/(ng/mL), so it enters the dimensionless bracket
    # as slp/100 with Cc in ng/mL.
    d/dt(clearance_capacity) <- kcl - kcl * (1 + slp / 100 * Cc) * clearance_capacity
    clearance_capacity(0)    <- 1

    # --- Disposition ----------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl / vc) * clearance_capacity * central -
                          (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    alag(depot) <- tlag

    Cc ~ lnorm(expSd)
  })
}
