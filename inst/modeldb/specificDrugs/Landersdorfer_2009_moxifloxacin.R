Landersdorfer_2009_moxifloxacin <- function() {
  description <- "Population PK model for oral moxifloxacin bone penetration (Landersdorfer 2009): two-compartment plasma disposition with first-order absorption from a gut depot, plus two paper-mechanistic bone matrix compartments (cortical and cancellous bone) connected to the central compartment by fixed transfer rate constants. The bone tissue:serum equilibrium concentration ratio is captured by the multiplicative scale terms fcortical and fcancellous on the cortical and cancellous bone observations. Disposition parameters were MAP-Bayesian estimated against Simon 1997 priors; bone-penetration scale terms used noninformative priors. Single 400 mg oral dose in 24 adults undergoing total hip replacement; serum and femoral bone samples (cortical + cancellous, head + neck) collected 2 to 7 hours post-dose."
  reference <- paste(
    "Landersdorfer CB, Kinzig M, Hennig FF, Bulitta JB,",
    "Holzgrabe U, Drusano GL, Sorgel F, Gusinde J.",
    "Penetration of moxifloxacin into bone evaluated by Monte Carlo",
    "simulation.",
    "Antimicrob Agents Chemother. 2009 May;53(5):2074-81.",
    "doi:10.1128/AAC.01056-08. PMID 19237653.",
    sep = " "
  )
  vignette <- "Landersdorfer_2009_moxifloxacin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  paper_specific_compartments <- c("cortical", "cancellous")

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Population mean 76.8 kg (SD 13.4) per Materials and Methods 'Study participants'; weight was not retained as a covariate on any PK parameter in the final model."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Population mean 63 years (SD 15) per Materials and Methods 'Study participants'; not retained as a covariate in the final model."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "14 of 24 subjects female per Materials and Methods 'Study participants'; not retained as a covariate in the final model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    age_range      = "63 +/- 15 years (mean +/- SD; specific range not tabulated)",
    age_median     = "63 years (mean)",
    weight_range   = "76.8 +/- 13.4 kg (mean +/- SD; specific range not tabulated)",
    weight_median  = "76.8 kg (mean)",
    sex_female_pct = 58,
    disease_state  = "Adults with coxarthrosis undergoing elective total hip replacement, no joint inflammation. Twenty subjects also received perioperative intravenous amoxicillin-clavulanate, three received levofloxacin, and one received clindamycin as standard antibacterial prophylaxis given in parallel to the study dose of moxifloxacin.",
    dose_range     = "Single 400 mg oral dose (Avalox tablet) administered 2 to 7 hours before surgery.",
    regions        = "Germany (University of Erlangen).",
    notes          = "Sparse PK sampling: one serum sample (prior to femoral bone resection) plus one bone sample (femoral head, with or without femoral neck; separated into cortical and cancellous tissue) per subject. Demographic and sampling details from Materials and Methods 'Study participants' and 'Sampling schedule'."
  )

  ini({
    # =========================================================================
    # Disposition parameters - MAP-Bayesian / ADAPT II final estimates.
    # Reported as median (range) across the 24 individual MAP-Bayesian
    # estimates in Table 1 column 'MAP-Bayesian (ADAPT II)'. Apparent (CL/F,
    # V/F) per Materials and Methods 'Structural model' paragraph 3.
    # Informative log-normal priors for the disposition parameters were
    # drawn from Simon et al. 1997 and three other moxifloxacin PK studies;
    # ka was estimated via NONMEM V against these disposition priors.
    # =========================================================================
    lka  <- log(1.6);  label("Absorption rate constant ka (1/h)")                            # Results 'PK analysis' paragraph 1; absorption half-life 26 min from NONMEM V
    lcl  <- log(10.8); label("Apparent total clearance CL/F (L/h)")                          # Table 1 MAP-Bayesian median 10.8 (range 9.85-11.5)
    lvc  <- log(62.0); label("Apparent central volume of distribution V_Central/F (L)")     # Table 1 MAP-Bayesian median 62.0 (range 58.5-65.4)
    lvp  <- log(59.5); label("Apparent peripheral volume of distribution V_Peripheral/F (L)") # Table 1 MAP-Bayesian median 59.5 (range 48.0-71.6)
    lq   <- log(18.9); label("Apparent intercompartmental clearance CL_ic/F (L/h)")          # Table 1 MAP-Bayesian median 18.9 (range 15.3-23.2)

    # =========================================================================
    # Bone tissue:serum equilibrium concentration ratios (paper-mechanistic).
    # F_cortical and F_cancellous are multiplicative scale terms on the
    # observed cortical and cancellous bone concentrations relative to the
    # corresponding compartmental concentrations (Materials and Methods
    # 'Structural model' paragraph 4). 'A value equal to 1 means that
    # concentrations after a continuous infusion at steady state are the
    # same in cortical bone and in serum'. Noninformative uniform priors in
    # the MAP-Bayesian step; CV reported in Table 1 is the log-normal
    # between-subject variability of the individual MAP-Bayesian estimates.
    # =========================================================================
    lfcortical   <- log(0.803); label("Cortical bone:serum equilibrium concentration ratio (unitless)")   # Table 1 MAP-Bayesian median 0.803 (35% CV, range 0.185-1.71)
    lfcancellous <- log(0.775); label("Cancellous bone:serum equilibrium concentration ratio (unitless)") # Table 1 MAP-Bayesian median 0.775 (48% CV, range 0.278-1.56)

    # =========================================================================
    # Inter-individual variability.
    # omega^2 = log(CV^2 + 1) for log-normal IIV.
    # The source paper reports an explicit CV% only for F_cortical (35%) and
    # F_cancellous (48%) in Table 1; the disposition parameters (CL, Vc, Vp,
    # CL_ic, ka) were MAP-Bayesian estimated with informative priors whose
    # between-subject SDs are taken from Simon et al. 1997 and are not
    # reproduced numerically in this paper. IIV on the disposition
    # parameters is therefore not encoded here; see the validation vignette
    # 'Assumptions and deviations' section for the impact on Monte Carlo
    # simulations and the user-facing override pattern.
    # =========================================================================
    etalfcortical   ~ 0.115573  # 35% CV; log(1 + 0.35^2) = 0.115573 - Table 1 MAP-Bayesian
    etalfcancellous ~ 0.207436  # 48% CV; log(1 + 0.48^2) = 0.207436 - Table 1 MAP-Bayesian

    # =========================================================================
    # Residual error - proportional model on each of serum, cortical bone,
    # and cancellous bone concentrations per Materials and Methods
    # 'MAP-Bayesian estimation' paragraph: 'The residual unidentified
    # variability was described by a proportional error model for the serum
    # and bone concentrations.' Numerical residual-error magnitudes are not
    # reported in Table 1 or elsewhere in the source paper. The values
    # below are 0 placeholders so the structural model loads; users running
    # stochastic simulations or fitting to new data should override propSd,
    # propSd_Ccortical, and propSd_Ccancellous. See the validation
    # vignette 'Assumptions and deviations' section.
    # =========================================================================
    propSd             <- fixed(0); label("Proportional residual error on serum Cc (fraction; placeholder)")
    propSd_Ccortical   <- fixed(0); label("Proportional residual error on cortical bone Ccortical (fraction; placeholder)")
    propSd_Ccancellous <- fixed(0); label("Proportional residual error on cancellous bone Ccancellous (fraction; placeholder)")
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Individual PK parameters - no covariates in the final model.
    #    Disposition parameters are typical-value only (no IIV encoded,
    #    see ini() rationale); F_cortical and F_cancellous carry log-normal
    #    IIV from Table 1.
    # -----------------------------------------------------------------------
    ka          <- exp(lka)
    cl          <- exp(lcl)
    vc          <- exp(lvc)
    vp          <- exp(lvp)
    q           <- exp(lq)
    fcortical   <- exp(lfcortical   + etalfcortical)
    fcancellous <- exp(lfcancellous + etalfcancellous)

    # -----------------------------------------------------------------------
    # 2. Paper-fixed bone-compartment structural constants.
    #    Materials and Methods 'PK modeling approach' paragraph 3:
    #    'a volume of distribution of 0.5 liter each for the cortical and
    #    cancellous bone compartments, which is equivalent to fixing k24
    #    and k25 to 0.022 h-1 in our model'.
    #    Equilibration rate constants k42 and k52 (bone -> central) are
    #    fixed to an equilibration half-life of 15 min per Materials and
    #    Methods 'PK modeling approach' paragraph 2.
    # -----------------------------------------------------------------------
    vbone     <- 0.5             # L; fixed - Materials and Methods 'PK modeling approach' paragraph 3
    k_cb_in   <- 0.022           # 1/h; fixed - k24 = k25 (central -> bone) - Materials and Methods 'PK modeling approach' paragraph 3
    k_cb_out  <- log(2) / 0.25   # 1/h; fixed - k42 = k52 (bone -> central), 15 min equilibration half-life - Materials and Methods 'PK modeling approach' paragraph 2

    # -----------------------------------------------------------------------
    # 3. Plasma micro-constants.
    # -----------------------------------------------------------------------
    kel <- cl / vc
    k23 <- q  / vc
    k32 <- q  / vp

    # -----------------------------------------------------------------------
    # 4. ODE system (Materials and Methods 'Structural model' equations).
    #    States carry drug amount (mg). Gut/depot is X1, central X2,
    #    peripheral disposition X3, cortical bone X4, cancellous bone X5.
    # -----------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k23 * central + k32 * peripheral1 - k_cb_in * central + k_cb_out * cortical - k_cb_in * central + k_cb_out * cancellous
    d/dt(peripheral1) <-  k23 * central - k32 * peripheral1
    d/dt(cortical)    <-  k_cb_in * central - k_cb_out * cortical
    d/dt(cancellous)  <-  k_cb_in * central - k_cb_out * cancellous

    # -----------------------------------------------------------------------
    # 5. Observed concentrations (mg/L). Dose in mg and vc in L give
    #    central/vc directly in mg/L. The bone compartments hold amount in
    #    mg; the compartmental bone concentration is cortical/vbone, which
    #    is then scaled by fcortical (fcancellous) to give the observed
    #    bone:serum ratio at equilibrium.
    # -----------------------------------------------------------------------
    Cc          <- central / vc
    Ccortical   <- fcortical   * cortical   / vbone
    Ccancellous <- fcancellous * cancellous / vbone

    Cc          ~ prop(propSd)
    Ccortical   ~ prop(propSd_Ccortical)
    Ccancellous ~ prop(propSd_Ccancellous)
  })
}
