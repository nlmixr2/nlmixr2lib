Clinckers_2008_MHD_rat <- function() {
  description <- paste(
    "Preclinical (rat). Population PK model for 10,11-dihydro-10-hydroxy-",
    "carbamazepine (MHD), the active metabolite of oxcarbazepine, in male",
    "Wistar rat plasma and hippocampal extracellular fluid (Clinckers 2008).",
    "One-compartment central disposition (V2) with combined zero-order",
    "(fraction F1 of dose over duration D2) and lagged first-order (1 - F1,",
    "ka with lag ALAG1) absorption after intraperitoneal bolus, coupled to",
    "a biophase / effect compartment (V3) reached via inter-compartmental",
    "rate constants k23 and k32. Acute focal pilocarpine-induced seizure",
    "activity and local intrahippocampal verapamil (efflux-transporter",
    "blockade) each shrink the biophase volume (V3a -> V3b under seizure;",
    "V3a -> V3c under verapamil); plasma kinetics are unaffected.")
  reference <- paste(
    "Clinckers R, Smolders I, Michotte Y, Ebinger G, Danhof M, Voskuyl RA,",
    "Della Pasqua O. Impact of efflux transporters and of seizures on the",
    "pharmacokinetics of oxcarbazepine metabolite in the rat brain.",
    "Br J Pharmacol. 2008;155(7):1127-1138. doi:10.1038/bjp.2008.366")
  vignette <- "Clinckers_2008_MHD_rat"
  units <- list(
    time          = "min",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    SEIZURE_ACUTE = list(
      description        = "Acute focal pilocarpine-induced seizure activity indicator: 1 = animal is undergoing intrahippocampal pilocarpine-evoked limbic seizures during the modelled observation window; 0 = no seizures.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no acute seizure activity)",
      notes              = "Time-fixed per animal in the source study (each rat is allocated to a single treatment arm; transient seizures are induced 30 min post-dose and subside before the end of sampling). Selects V3b (biophase volume during seizures) over V3a in the model. Mutually exclusive with EFFLUX_INHIB in the original Clinckers 2008 design.",
      source_name        = "A"
    ),
    EFFLUX_INHIB = list(
      description        = "Local intrahippocampal efflux-transporter inhibitor co-perfusion indicator: 1 = the brain microdialysis probe is co-perfused with verapamil (5 mM, P-glycoprotein inhibitor); 0 = no inhibitor.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no local efflux-transporter blockade)",
      notes              = "Time-fixed per animal in the source study. Selects V3c (biophase volume during efflux blockade) over V3a in the model. Mutually exclusive with SEIZURE_ACUTE in the original Clinckers 2008 design.",
      source_name        = "B"
    )
  )

  population <- list(
    species        = "rat (Wistar, male)",
    n_subjects     = 32L,
    n_studies      = 1L,
    age_range      = "Adult; specific age not reported",
    weight_range   = "260-320 g",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = paste(
      "Male Wistar albino rats (Iffa Credo). Subgroups: control",
      "(intrahippocampal vehicle); acute focal pilocarpine-induced limbic",
      "seizures (10 mM pilocarpine perfused 30-70 min post-dose); local",
      "efflux-transporter blockade (5 mM intrahippocampal verapamil)."
    ),
    dose_range     = paste(
      "20, 40, 60, 80, 100, 150 mg/kg single intraperitoneal bolus of",
      "MHD (suspension in propylene glycol/ethanol/saline 6:2:2). Doses",
      "20-60 mg/kg are sub-therapeutic; 80-150 mg/kg are anticonvulsant."
    ),
    regions        = "Brussels, Belgium (single laboratory)",
    notes          = paste(
      "Demographics from Clinckers 2008 Methods (Animals; Study design)",
      "and Table 1 (per-dose-group N). Plasma and hippocampal microdialysate",
      "MHD concentrations were quantified by LC-UV (LOQ 5 ng/mL).",
      "Microdialysis probe relative recovery was determined per sample via",
      "the internal reference technique using mCBZ as the internal standard."
    )
  )

  ini({
    # --------------------------------------------------------------------
    # Structural parameters. Clinckers 2008 Table 2 reports point estimates
    # in terms of rate constants (k, k23, k32) and volumes (V2, V3) rather
    # than CL/Q; the NONMEM ADVAN6 subroutine was used (Methods, Data
    # analysis: 'The ADVAN6 subroutine in NONMEM was used... The PK model
    # was parameterized in terms of rate constants and volumes of
    # distribution.'). Rate-constant parameterization is preserved here
    # so the THETA values map 1:1 onto the published estimates; the
    # asymmetric biophase kinetics (k23 != k32, plus three condition-
    # specific V3 values) cannot be re-expressed as a single Q. The
    # non-canonical names (lkel, lk23, lk32, lv3a/b/c, lalag, lfr, ld2)
    # are listed in the vignette Assumptions and deviations.
    # Volumes from Table 2 are converted from mL to L (V2 = 484 mL = 0.484 L,
    # etc.) so that concentrations expressed as central/vc are in
    # mg/L = ug/mL = 1e3 ng/mL (an explicit 1e3 multiplier converts to
    # the paper's ng/mL units in the observation block).
    # --------------------------------------------------------------------

    lka  <- log(0.22)
    label("First-order absorption rate constant ka (1/min)")                          # Clinckers 2008 Table 2: ka = 0.22 /min (control structural model)

    lkel <- log(0.0114)
    label("First-order elimination rate constant kel (paper k, 1/min)")               # Clinckers 2008 Table 2: k = 0.0114 /min

    lk23 <- log(0.0203)
    label("Inter-compartmental rate central -> biophase k23 (1/min)")                 # Clinckers 2008 Table 2: k23 = 0.0203 /min

    lk32 <- log(0.0306)
    label("Inter-compartmental rate biophase -> central k32 (1/min)")                 # Clinckers 2008 Table 2: k32 = 0.0306 /min

    lvc  <- log(0.484)
    label("Central volume of distribution V2 (L)")                                    # Clinckers 2008 Table 2: V2 = 484 mL = 0.484 L

    lv3a <- log(0.813)
    label("Biophase volume of distribution V3 under control conditions V3a (L)")      # Clinckers 2008 Table 2: V3a = 813 mL = 0.813 L

    lv3b <- log(0.563)
    label("Biophase volume of distribution V3 during acute seizures V3b (L)")         # Clinckers 2008 Table 2: V3b = 563 mL = 0.563 L (P < 0.05 vs V3a)

    lv3c <- log(0.407)
    label("Biophase volume of distribution V3 during efflux-transporter blockade V3c (L)")  # Clinckers 2008 Table 2: V3c = 407 mL = 0.407 L (P < 0.05 vs V3a)

    lalag <- log(3.23)
    label("Lag time for first-order absorption ALAG1 (min)")                          # Clinckers 2008 Table 2: ALAG1 = 3.23 min

    lfr  <- log(0.347)
    label("Log fraction of dose absorbed via zero-order route F1 (unitless)")         # Clinckers 2008 Table 2: F1 = 0.347 (the lognormal IIV form may exceed 1 in tails -- matches NONMEM ADVAN6 parameterisation)

    ld2  <- log(12.6)
    label("Duration of the zero-order absorption infusion D2 (min)")                  # Clinckers 2008 Table 2: D2 = 12.6 min

    # --------------------------------------------------------------------
    # Inter-individual variability. Clinckers 2008 Methods (eq. for P_i =
    # theta * exp(eta_i)) -- log-normal IIV on rate constants and on the
    # absorption fraction / duration. Table 2 reports omega^2 values for
    # k23, k32, V3a, V3c, F1, D2; no omega^2 reported for ka, kel, V2,
    # V3b (V3b applies only during seizures where IIV was not separately
    # estimated). Each rat is allocated to a single experimental condition,
    # so etalv3a applies under control or seizure animals and etalv3c
    # applies under verapamil animals.
    # --------------------------------------------------------------------
    etalk23 ~ 0.14                                                                    # Clinckers 2008 Table 2: omega^2(k23) = 0.14
    etalk32 ~ 0.153                                                                   # Clinckers 2008 Table 2: omega^2(k32) = 0.153
    etalv3a ~ 0.078                                                                   # Clinckers 2008 Table 2: omega^2(V3a) = 0.078
    etalv3c ~ 0.0187                                                                  # Clinckers 2008 Table 2: omega^2(V3c) = 0.0187
    etalfr  ~ 0.0975                                                                  # Clinckers 2008 Table 2: omega^2(F1) = 0.0975
    etald2  ~ 0.0575                                                                  # Clinckers 2008 Table 2: omega^2(D2) = 0.0575

    # --------------------------------------------------------------------
    # Residual error.
    #
    # Plasma (Clinckers 2008 Methods: C_obs = C_pred * (1 + eps1) + eps2,
    # combined proportional + additive). Table 2 reports SIGMA-style
    # variances; convert to SDs via sqrt:
    #   Proportional plasma sigma^2 = 0.0256 -> SD = 0.160 (16% CV)
    #   Additive plasma sigma^2     = 13.0   -> SD = 3.606 ng/mL
    #
    # Brain (Clinckers 2008 Methods: 'proportional error model' with W =
    # sqrt(prop_brain * F^2 + WT_weighing_factor); WT is the residual
    # weighing factor reported separately in Table 2). This W form is
    # algebraically equivalent to a combined proportional + additive
    # error on the brain output:
    #   Proportional brain prop coefficient = 0.605 -> SD = 0.778 (77.8% CV)
    #   Additive brain WT coefficient       = 0.0782 -> SD = 0.280 ng/mL
    # Bootstrap CV% of these two estimates is large (101.3% and 122.3%
    # respectively), reflecting the limited brain microdialysis
    # information content for the residual model (Clinckers 2008 Results
    # 'Population PK model selection': WT was the only parameter not
    # estimated with CV% < 15%, at 59.6%).
    # --------------------------------------------------------------------
    propSd       <- 0.160
    label("Plasma proportional residual SD (fraction)")                               # Clinckers 2008 Table 2: sigma^2 prop plasma = 0.0256 -> SD = sqrt(0.0256) = 0.160

    addSd        <- 3.606
    label("Plasma additive residual SD (ng/mL)")                                      # Clinckers 2008 Table 2: sigma^2 add plasma = 13.0 -> SD = sqrt(13.0) = 3.606

    propSd_Cbrain <- 0.778
    label("Brain (biophase) proportional residual SD (fraction)")                     # Clinckers 2008 Table 2: brain prop coef = 0.605 inside W = sqrt(prop*F^2 + WT) -> SD = sqrt(0.605) = 0.778

    addSd_Cbrain  <- 0.280
    label("Brain (biophase) additive residual SD (ng/mL)")                            # Clinckers 2008 Table 2: WT (residual weighing factor, additive coef inside W) = 0.0782 -> SD = sqrt(0.0782) = 0.280
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters (log-normal IIV per Methods).
    # ------------------------------------------------------------------
    ka  <- exp(lka)
    kel <- exp(lkel)
    k23 <- exp(lk23 + etalk23)
    k32 <- exp(lk32 + etalk32)

    vc  <- exp(lvc)

    # Biophase volume is condition-specific. Each animal sits in exactly
    # one of {control, seizure, verapamil}, so the three V3 components
    # are mutually exclusive switches:
    #   control     (SEIZURE_ACUTE = 0, EFFLUX_INHIB = 0) -> V3a * exp(etalv3a)
    #   seizure     (SEIZURE_ACUTE = 1, EFFLUX_INHIB = 0) -> V3b (no IIV)
    #   verapamil   (SEIZURE_ACUTE = 0, EFFLUX_INHIB = 1) -> V3c * exp(etalv3c)
    # Each per-condition v3 is computed on its own simple line so the
    # mu-referencing parser recognises the eta -> log-parameter mapping.
    v3a <- exp(lv3a + etalv3a)
    v3b <- exp(lv3b)
    v3c <- exp(lv3c + etalv3c)
    v3  <- v3a * (1 - SEIZURE_ACUTE) * (1 - EFFLUX_INHIB) +
           v3b *      SEIZURE_ACUTE  * (1 - EFFLUX_INHIB) +
           v3c * (1 - SEIZURE_ACUTE) *      EFFLUX_INHIB

    # Absorption modifiers
    alag <- exp(lalag)
    fr   <- exp(lfr + etalfr)
    d2   <- exp(ld2 + etald2)

    # ------------------------------------------------------------------
    # 2. ODEs. depot receives the lagged first-order portion (1 - fr) of
    #    the IP bolus and feeds central via ka. central also receives a
    #    direct zero-order infusion (the fr portion delivered over d2
    #    minutes) via rxode2's dur(central). effect is the biophase
    #    compartment (V3) reached via the inter-compartmental k23/k32.
    # ------------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central - k23 * central + k32 * effect
    d/dt(effect)  <-                                 k23 * central - k32 * effect

    # ------------------------------------------------------------------
    # 3. Absorption: split the IP dose between depot (lagged first-order,
    #    fraction 1 - fr) and central (zero-order, fraction fr, duration
    #    d2). The user supplies two simultaneous dose records, one to
    #    each compartment, both with the full IP amt; the bioavailability
    #    fractions handle the actual split.
    # ------------------------------------------------------------------
    f(depot)   <- 1 - fr
    lag(depot) <- alag

    f(central)   <- fr
    dur(central) <- d2

    # ------------------------------------------------------------------
    # 4. Observations. Volumes are in L and amounts in mg, so
    #    amount / volume = mg/L = ug/mL; multiply by 1000 to express
    #    concentrations in ng/mL (the source assay unit, LOQ 5 ng/mL).
    # ------------------------------------------------------------------
    Cc      <- central / vc * 1000
    Cbrain  <- effect  / v3 * 1000

    Cc     ~ prop(propSd)        + add(addSd)
    Cbrain ~ prop(propSd_Cbrain) + add(addSd_Cbrain)
  })
}
