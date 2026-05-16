Thorsted_2016_somatropin_human <- function() {
  description <- paste(
    "Translational (allometrically-scaled rat-to-human) population",
    "PKPD model for recombinant human growth hormone (rhGH /",
    "somatropin) in growth-hormone-deficient adult males. Structural",
    "parameter values are derived from the Thorsted 2016",
    "hypophysectomized-rat PKPD fit by allometric scaling to a 70 kg",
    "reference subject (Table 3 of the source paper): clearance terms",
    "(CL, Q) and Vmax with exponent 0.75; distribution volumes (Vc,",
    "Vp) with exponent 0.9 (the empirically-selected best-fit exponent",
    "for human i.v. data); first-order absorption rate constants (ka1,",
    "ka2) and kout with exponent -0.25; KM unscaled; Emax and EC50",
    "unscaled. The s.c. absorption model is the corrected form",
    "(Table 3 / Figure 5): bioavailability of the ka2 path reduced",
    "from 0.833 (rat) to 0.500, and one transit compartment added to",
    "the ka1 path. The IGF-1 indirect response uses kin = kout * R0",
    "with R0 fixed to 65 ng/mL (human population mean per Laursen",
    "1996) and is driven directly by plasma rhGH (no effect-delay",
    "chain - the rat CPLAG chain is intentionally dropped for the",
    "human prediction). Bodyweight gain is not included in the human",
    "model. Variability is inherited from the rat PKPD fit; residual",
    "error is fixed at the values used for the human-simulation",
    "validation (Methods)."
  )
  reference <- paste(
    "Thorsted A, Thygesen P, Agerso H, Laursen T, Kreilgaard M.",
    "Translational mixed-effects PKPD modelling of recombinant human",
    "growth hormone - from hypophysectomized rat to patients.",
    "Br J Pharmacol. 2016 Jun;173(11):1742-55.",
    "doi:10.1111/bph.13473.",
    "Parameters scaled from the rat fit; see also",
    "modellib('Thorsted_2016_somatropin_rat')."
  )
  vignette <- "Thorsted_2016_somatropin_human"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ng/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling around a reference adult body weight of 70",
        "kg per Thorsted 2016 Table 3 footnote 'a' (70 kg man).",
        "Clearance terms (CL, Q) and the Michaelis-Menten maximum",
        "elimination capacity (Vmax) scale with (WT/70)^0.75;",
        "distribution volumes (Vc, Vp) scale with (WT/70)^0.9 (the",
        "exponent that best fit the human i.v. data per Figure 4);",
        "first-order absorption rate constants (ka1, ka2) and the",
        "IGF-1 degradation rate (kout) scale with (WT/70)^(-0.25). KM",
        "and the Emax / EC50 / F1 / F2 parameters are unscaled. The",
        "exponents are fixed at their canonical allometric values per",
        "Thorsted 2016 Methods (Prediction of human PKPD).",
        sep = " "
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "Adult males (Laursen 1996 validation cohort)",
    weight_range   = "approximately 70 kg (Table 3 reference)",
    sex_female_pct = 0,
    disease_state  = paste(
      "Growth hormone deficient adult males withdrawn from rhGH",
      "treatment for 1 week to washout IGF-1 prior to study (Laursen",
      "1996); model validated against 5 ug/kg i.v. and 33.5 ug/kg s.c.",
      "rhGH single-dose data."
    ),
    dose_range     = "5 ug/kg i.v. bolus and 33.5 ug/kg s.c. single doses (Laursen 1996 validation)",
    notes          = paste(
      "Translational forward-simulation construct rather than a fit to",
      "human data: parameter values are obtained by allometrically",
      "scaling the Thorsted 2016 hypophysectomized-rat fit (Table 2",
      "rat values) to a 70 kg human (Table 3 predicted values) with",
      "the s.c. absorption model empirically corrected (Figure 5).",
      "The Laursen 1996 dataset (10 GH-deficient males) is the",
      "validation reference; predictions of rhGH concentrations and",
      "IGF-1 induction were checked by simulating the dataset 2000",
      "times with the inherited rat IIVs and the fixed residual-error",
      "values described in Methods (Prediction of human PKPD)."
    )
  )

  ini({
    # PK structural parameters - Thorsted 2016 Table 3 predicted human
    # values for a 70 kg man (allometric scaling from rat Table 2:
    # CL, Q, Vmax with exponent 0.75; Vc, Vp with exponent 0.9;
    # ka1, ka2 with exponent -0.25; KM unscaled).
    lcl     <- log(3.88);     label("Linear clearance CL (L/h) at 70 kg")                  # Table 3
    lvmax   <- log(1565);     label("Michaelis-Menten Vmax (ug/h) at 70 kg")               # Table 3
    lkm     <- log(358);      label("Michaelis-Menten KM (ug/L = ng/mL)")                  # Table 3 (text: KM not scaled; Table 3 column lists 403 but Methods state KM was not scaled - see vignette Errata)
    lvc     <- log(2.51);     label("Central volume of distribution Vc (L) at 70 kg")     # Table 3
    lq      <- log(1.37);     label("Inter-compartmental clearance Q (L/h) at 70 kg")     # Table 3
    lvp     <- log(2.94);     label("Peripheral volume of distribution Vp (L) at 70 kg")  # Table 3
    lka1    <- log(0.587);    label("First-order SC absorption rate ka1 (1/h) at 70 kg - direct path with added transit") # Table 3
    lka2    <- log(0.237);    label("First-order SC absorption rate ka2 (1/h) at 70 kg - delayed via transit") # Table 3

    # Bioavailabilities - corrected SC absorption model (paper Figure 5
    # narrative: F2 reduced from 0.833 to 0.50 and one transit
    # compartment added to the ka1 path).
    lfdepot  <- log(0.0316);  label("Bioavailability F1 through ka1 SC path (fraction)")   # Table 3 (inherited from rat)
    lfdepot2 <- log(0.500);   label("Bioavailability F2 through ka2 SC path (fraction) - corrected from rat 0.833 to 0.50") # Results - corrected SC absorption

    # PD structural parameters - Thorsted 2016 Table 3 human predicted.
    # The rat CPLAG effect-delay chain is intentionally dropped for the
    # human prediction (paper: "The effect delay for rats was not
    # included for predictions, as patients had only been deprived of
    # rhGH for 1 week in order to washout IGF-1").
    lkout   <- log(0.0178);   label("First-order IGF-1 degradation rate kout (1/h) at 70 kg") # Table 3
    lr0     <- fixed(log(65)); label("Baseline IGF-1 R0 (ng/mL) - fixed to human population mean (Laursen 1996)") # Table 3 / Methods
    lemax   <- log(9.88);     label("Maximum stimulation of kin relative to baseline (unscaled from rat)") # Table 3
    lec50   <- log(16.3);     label("rhGH concentration for 50% Emax (ug/L = ng/mL; unscaled from rat)") # Table 3

    # Allometric exponents - fixed at canonical values per Methods.
    e_wt_cl   <- fixed(0.75);  label("Allometric exponent on CL, Q, Vmax (unitless)")      # Methods
    e_wt_vc   <- fixed(0.9);   label("Allometric exponent on Vc, Vp (unitless) - human")   # Methods / Figure 4
    e_wt_ka   <- fixed(-0.25); label("Allometric exponent on ka1, ka2, kout (unitless)")   # Methods

    # IIV - inherited from the rat fit (Thorsted 2016 Methods: "with
    # variability as estimated in the PKPD model"). CV% from Table 2
    # converted to omega^2 via log(CV^2 + 1). The CL-Vp block is
    # correlated with rho = -0.568 (Results); off-diagonal is
    # rho * sqrt(var_cl) * sqrt(var_vp). CL 11.6% CV -> 0.01337;
    # Vp 18.4% CV -> 0.03330; cov -> -0.01199.
    etalcl + etalvp ~ c(0.01337,
                        -0.01199, 0.03330)                                                  # Table 2 IIV CL 11.6% CV, Vp 18.4% CV, correlation -0.568
    etalka2  ~ 0.008613                                                                     # Table 2 IIV ka2 9.3% CV (Table 3's 19.3% appears inconsistent - see Errata)
    etalr0   ~ 0.028504                                                                     # Table 2 IIV R0 17.0% CV: log(0.170^2 + 1)

    # Residual error - Thorsted 2016 Methods (Prediction of human
    # PKPD): fixed values used in the simulation-based validation.
    propSd       <- fixed(0.25); label("Proportional residual error on rhGH plasma concentration (fraction)")  # Methods - human validation
    addSd        <- fixed(0.25); label("Additive residual SD on rhGH plasma concentration (ng/mL = ug/L)")      # Methods - human validation
    addSd_IGF1   <- fixed(10);   label("Additive residual SD on IGF-1 plasma concentration (ng/mL)")            # Methods - human validation
  })

  model({
    # 1. Allometric scaling (WT covariate) around reference 70 kg.
    wt_norm <- WT / 70

    # 2. Individual PK parameters
    cl   <- exp(lcl   + etalcl) * wt_norm^e_wt_cl
    vmax <- exp(lvmax)          * wt_norm^e_wt_cl
    km   <- exp(lkm)
    vc   <- exp(lvc)            * wt_norm^e_wt_vc
    q    <- exp(lq)             * wt_norm^e_wt_cl
    vp   <- exp(lvp   + etalvp) * wt_norm^e_wt_vc
    ka1  <- exp(lka1)           * wt_norm^e_wt_ka
    ka2  <- exp(lka2  + etalka2)* wt_norm^e_wt_ka
    fdepot  <- exp(lfdepot)
    fdepot2 <- exp(lfdepot2)

    # Individual PD parameters
    kout <- exp(lkout)          * wt_norm^e_wt_ka
    r0   <- exp(lr0 + etalr0)
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    # 3. Micro-constants for the 2-compartment disposition
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. ODE system
    Cc <- central / vc

    # Parallel SC absorption with one transit compartment on each
    # path (corrected human absorption model per Figure 5 narrative):
    # depot -> transit1 -> central (ka1 path with added transit) and
    # depot2 -> transit2 -> central (ka2 path with original transit).
    d/dt(depot)       <- -ka1 * depot
    d/dt(transit1)    <-  ka1 * depot  - ka1 * transit1
    d/dt(depot2)      <- -ka2 * depot2
    d/dt(transit2)    <-  ka2 * depot2 - ka2 * transit2
    d/dt(central)     <-  ka1 * transit1 + ka2 * transit2 -
                          kel * central - k12 * central + k21 * peripheral1 -
                          (vmax * Cc) / (km + Cc)
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # IGF-1 indirect response (no effect-delay chain for the human
    # prediction): kin = kout * R0 stimulated directly by plasma Cc via
    # Emax / EC50.
    stim <- emax * Cc / (ec50 + Cc)
    d/dt(igf1) <- kout * r0 * (1 + stim) - kout * igf1
    igf1(0)    <- r0

    # 5. Bioavailability: SC doses go to depot (ka1 path) with F1 and
    # to depot2 (ka2 path) with F2. IV doses are placed into central
    # directly with no F applied.
    f(depot)  <- fdepot
    f(depot2) <- fdepot2

    # 6. Observations and residual error
    IGF1 <- igf1

    Cc   ~ add(addSd) + prop(propSd)
    IGF1 ~ add(addSd_IGF1)
  })
}
