Mohamed_2016_colistin_meropenem <- function() {
  description <- "In vitro time-kill PK/PD model for colistin and meropenem alone and in combination against P. aeruginosa wild-type ATCC 27853 (BACT=2) and meropenem-resistant clinical isolate ARU552 (BACT=1); proliferating/resting bacterial subpopulations with strain-specific drug effects, colistin adaptive resistance, and a meropenem-resistant mutant subpopulation"
  reference <- paste(
    "Mohamed A. F., Kristoffersson A. N., Karvanen M., Nielsen E. I., Cars O., Friberg L. E. (2016).",
    "Dynamic interaction of colistin and meropenem on a WT and a resistant strain of",
    "Pseudomonas aeruginosa as quantified in a PK/PD model.",
    "J Antimicrob Chemother 71(5):1279-1290.",
    "doi:10.1093/jac/dkv488.",
    "DDMORE Foundation Model Repository: DDMODEL00000173.",
    sep = " "
  )
  vignette <- "Mohamed_2016_colistin_meropenem"
  units <- list(time = "h", dosing = "CFU/mL", concentration = "mg/L")

  ddmore_id <- "DDMODEL00000173"
  replicate_of <- NULL

  covariateData <- list(
    BACT = list(
      description = "Bacterial strain indicator: 1 = ARU552 (meropenem-resistant clinical isolate), 2 = ATCC 27853 (wild-type)",
      units = "(categorical)",
      type = "categorical",
      reference_category = "2 (ATCC)",
      notes = "Hard switch on >25 strain-specific structural parameters. Source: Executable_ColistinMeropenem_Interaction.mod $INPUT.",
      source_name = "BACT"
    ),
    TYPE = list(
      description = "Experiment-condition indicator: 1 = meropenem alone, 2 = colistin alone, 3 = combination, 4/5 = control (no drug)",
      units = "(categorical)",
      type = "categorical",
      reference_category = "4 (control)",
      notes = "Drives MERO/COL/CONTR indicator flags inside model(). Source: .mod $INPUT.",
      source_name = "TYPE"
    ),
    Ccol = list(
      description = "Initial colistin concentration in the experimental well at t = 0",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Sets the colistin compartment initial value (A_0(6) in NONMEM). Also drives the discrete KE2 (colistin in vitro elimination) lookup table inside model(). Source: .mod $INPUT line 33 and $PK A_0(6) initialization.",
      source_name = "Ccol"
    ),
    Cmer = list(
      description = "Initial meropenem concentration in the experimental well at t = 0",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Sets the meropenem compartment initial value (A_0(3) in NONMEM). Source: .mod $INPUT line 32 and $PK A_0(3) initialization.",
      source_name = "Cmer"
    ),
    DILcol = list(
      description = "Per-row bath dilution rate applied to the colistin compartment",
      units = "1/h",
      type = "continuous",
      reference_category = NULL,
      notes = "Read directly from the dataset on each integration interval; appears additively with KE2 in d/dt(col). Non-zero only for dynamic (ID2 > 18) experiments.",
      source_name = "DILcol"
    ),
    DILmer = list(
      description = "Per-row bath dilution rate applied to the meropenem compartment",
      units = "1/h",
      type = "continuous",
      reference_category = NULL,
      notes = "Read directly from the dataset; additive with KE1 in d/dt(mero). Non-zero only for dynamic experiments.",
      source_name = "DILmer"
    ),
    # NOTE: REPL (replicate index) and BLOQ (below-LOQ censoring flag) are
    # present in the source dataset but are NOT referenced inside model().
    # REPL drives a per-replicate residual EPS slot in NONMEM that this
    # extraction collapses into a single combined SD (all SIGMAs are FIX
    # equal). BLOQ is handled by nlmixr2 via the standard `cens` event
    # column. They are intentionally omitted from covariateData per the
    # checkModelConventions() rule that covariateData entries must be
    # referenced in model().
    ID2 = list(
      description = "Within-subject experiment index used to flag dynamic-system experiments (ID2 > 18)",
      units = "(categorical)",
      type = "categorical",
      reference_category = NULL,
      notes = "Drives the dynamic-system growth-rate scaling KGDxS = 1 + THETA(35) when ID2 > 18; KGDxS = 1 otherwise. Source: .mod $PK lines 126-127.",
      source_name = "ID2"
    )
  )

  population <- list(
    n_subjects = NA_integer_,
    n_studies = 1,
    age_range = NA_character_,
    weight_range = NA_character_,
    sex_female_pct = NA_real_,
    disease_state = "in vitro bacterial culture (no human subjects)",
    dose_range = "Initial inoculum approximately 2.5e5 to 3.5e5 CFU/mL; static colistin 0.04-24 mg/L (ATCC) or 0.08-24 mg/L (ARU); static meropenem; combination experiments at colistin 0.25-8 mg/L (intended) with paired meropenem; dynamic time-kill experiments use bath-flow dilution rates DILcol / DILmer per the dataset",
    regions = NA_character_,
    notes = "In vitro pharmacodynamic model fit to time-kill experiments against Pseudomonas aeruginosa wild-type ATCC 27853 (BACT=2) and meropenem-resistant clinical isolate ARU552 (BACT=1). Static and dynamic cultures, with replicate-specific residual error. Final estimates from Output_real_ColistinMeropenem_interaction.lst (NONMEM7 LAPLACIAN evaluation, OBJV = 1997.721)."
  )

  ini({
    # Final estimates from Output_real_ColistinMeropenem_interaction.lst
    # FINAL PARAMETER ESTIMATE block (lines 855-857). All values match the
    # .mod $THETA initial values because the production run was a MAXEVAL=0
    # EVALUATION (msfb282-final_tidyup) at the converged estimates. THETA
    # numbers below correspond to .mod $THETA lines 485-519. Parameter naming
    # follows .mod symbology (lowercase) since the model is mechanistic and
    # the values are point estimates; no IIV / no estimation is intended in
    # nlmixr2lib usage.

    # Bacterial growth rates (1/h)
    kga      <- c(0, 1.08, 10);     label("ATCC growth rate KGA (1/h)")          # TH1, .mod L485
    kgp      <- c(0, 0.814, 10);    label("ARU552 growth rate KGP (1/h)")        # TH2, .mod L486
    kk       <- fixed(0.179);       label("Natural bacterial death rate (1/h)")  # TH3 FIX, .mod L487

    # Colistin drug-effect parameters
    emax_col_atcc  <- fixed(282);   label("Emax colistin ATCC (1/h)")            # TH4 FIX, .mod L488
    slop1_col_aru  <- c(0, 17.2, 50); label("Slope colistin ARU (linear, 1/h per (mg/L)^GAMC)") # TH5, .mod L489
    ec51_col_atcc  <- c(0, 1.84, 20); label("EC50 colistin ATCC (mg/L)")        # TH8, .mod L492
    gamc_col       <- fixed(1);     label("Hill coefficient colistin GAMC")      # TH11 FIX, .mod L495

    # Meropenem drug-effect parameters
    # NOTE: THETA(6) (.mod L490, "EMAX ATCC MERO") is FIX 0 in the source and
    #       is never referenced in $DES (the ATCC + meropenem branch uses the
    #       linear-power form SLOPAM*CA2^GAMA, not Emax). The .mod assigns
    #       TVEMAX2 = THETA(6) but EMAX2 is dead code. Omitted from ini().
    emax_mero_aru  <- c(0, 1.74, 100); label("Emax meropenem ARU (1/h)")        # TH7, .mod L491
    slopam_mero_atcc <- c(0, 2.16, 100); label("Slope (power) meropenem ATCC (1/h per (mg/L)^GAMA)") # TH9, .mod L493
    ec53_mero_aru  <- c(0, 17.7, 50); label("EC50 meropenem ARU (mg/L)")         # TH10, .mod L494
    gama_mero_atcc <- c(0, 0.376, 100); label("Hill/power exponent meropenem ATCC GAMA") # TH12, .mod L496
    gamp_mero_aru  <- c(0, 2.79, 10); label("Hill exponent meropenem ARU GAMP")   # TH26, .mod L510

    # Bmax (carrying capacity), parameterized as log10(BMAX)
    bmax_log10     <- c(8, 9.11, 10); label("log10 of maximum bacterial count BMAX (log10 CFU/mL); BMAX = 10^bmax_log10") # TH13, .mod L497

    # Colistin adaptive-resistance parameters
    kon51_col_atcc <- c(0, 2.42, 10); label("KON51 EC50 of colistin-resistance development for ATCC (mg/L)") # TH14, .mod L498
    kon52_col_aru  <- c(0, 18.6, 10000); label("KON52 EC50 of colistin-resistance development for ARU (mg/L)") # TH15, .mod L499
    koffc_col      <- c(0, 0.103, 5); label("KOFFC return-to-susceptibility rate colistin (1/h)") # TH16, .mod L500
    konca_col_atcc <- c(0, 0.779, 100); label("KONCA Emax of colistin-resistance development for ATCC (1/h)") # TH18, .mod L502
    koncp_col_aru  <- c(0, 2.68, 100); label("KONCP Emax of colistin-resistance development for ARU (1/h)") # TH19, .mod L503
    gam3c_col      <- c(0, 1.08, 10); label("GAM3C Hill exponent of colistin-resistance development") # TH22, .mod L506
    gam1_col       <- c(0, 0.0245, 10); label("GAM1 inhibition exponent of colistin drug effect by adaptive resistance") # TH23, .mod L507

    # Meropenem adaptive-resistance parameters (declared but not used in this version of the model)
    koffm_mero     <- fixed(0);     label("KOFFM return-to-susceptibility rate meropenem (FIX 0; meropenem adaptive-resistance branch unused)") # TH17 FIX, .mod L501
    konma_mero_atcc <- fixed(0.113); label("KONMA Emax of meropenem-resistance development for ATCC (FIX, branch unused)") # TH20 FIX, .mod L504
    konmp_mero_aru  <- fixed(0.0821); label("KONMP Emax of meropenem-resistance development for ARU (FIX, branch unused)") # TH21 FIX, .mod L505

    # Combination interaction exponents (can be negative -> can't log-transform)
    inter_atcc <- c(-10, -2.58, 10); label("Interaction exponent ATCC (THETA(27); >0 synergistic, <0 antagonistic)") # TH27, .mod L511
    inter_aru  <- c(-10, 1.30, 10); label("Interaction exponent ARU (THETA(28))") # TH28, .mod L512

    # Mutant subpopulation parameters
    # NOTE: THETA(30) (.mod L514, "-LOG10 MUTF ARU") is FIX 0 and unused in
    #       equations: .mod L148 sets MUTP = MUTA unconditionally, so
    #       THETA(30) does not enter the model. Omitted from ini().
    neg_log10_mutf_atcc <- c(0, 2.83, 10); label("-log10 of mutant fraction MUTF ATCC (THETA(29))") # TH29, .mod L513
    concshift_mut <- c(0, 22.5, 100); label("Mutant-subpopulation concentration-shift fold (THETA(31); EC53M = EC53 * (1 + concshift_mut))") # TH31, .mod L515
    kg_decrease_atcc <- fixed(0); label("Mutant KG fractional decrease ATCC (FIX 0; THETA(33))") # TH33 FIX, .mod L517
    kg_decrease_aru  <- c(0, 0.536, 1); label("Mutant KG fractional decrease ARU (THETA(34))") # TH34, .mod L518

    # Dynamic-system growth-rate scaling
    kg_dyn_increase  <- c(0, 0.332, 10); label("Fractional KG increase for dynamic-system experiments ID2 > 18 (THETA(35))") # TH35, .mod L519

    # NOTE: THETA(24), THETA(25), THETA(32) (.mod L508-509, L516) are FIX
    #       placeholders that the source declares but never references in
    #       any equation. Omitted from ini() for the same reason as
    #       THETA(6) and THETA(30): nlmixr2 requires every ini() entry to
    #       appear in model(). The original .mod $THETA values (1, 1, 0
    #       respectively) are recorded here for traceability.

    # Residual error.
    # NONMEM model: Y = IPRED + LOG(10)*(EPS(1) + AAA*EPS(2) + AAB*EPS(3) + AAC*EPS(4) + AAD*EPS(5))
    # where AAA..AAD are one-hot replicate indicators (REPL=1..4) and IPRED is on natural log scale.
    # SIGMA(1) = 0.422 (common, log10 var); SIGMA(2..5) = 0.0364 each (replicate-specific, log10 var, BLOCK SAME).
    # All SIGMAs are FIXED. Per-row variance is (SIGMA(1) + SIGMA(REPL+1)) on log10 scale; combined SD on
    # natural-log scale = sqrt(SIGMA(1) + SIGMA(2)) * log(10) (since all replicate-specific SIGMAs are equal).
    # Verbatim collapse to a single additive residual on the log-bacteria observation:
    addSd_logBact <- fixed(sqrt(0.422 + 0.0364) * log(10)); label("Combined residual SD on natural-log scale (sqrt(SIGMA(1)+SIGMA(2))*log(10) per .lst SIGMA block lines 877-889)") # SIGMA, .mod L521-526
  })

  model({
    # =====================================================================
    # Verbatim translation of Executable_ColistinMeropenem_Interaction.mod
    # ($PK lines 48-277; $DES lines 280-437). Strain switches on BACT
    # (1 = ARU552, 2 = ATCC 27853); experiment-type switches on TYPE
    # (1 = mero alone, 2 = col alone, 3 = combination, 4/5 = control).
    # All ETA(1) usage is dropped because OMEGA is FIX 0 in the source.
    # =====================================================================

    # Indicator flags from TYPE
    MERO  <- 0
    if (TYPE == 1) MERO <- 1
    if (TYPE == 3) MERO <- 1
    COL   <- 0
    if (TYPE == 2) COL <- 1
    if (TYPE == 3) COL <- 1
    CONTR <- 0
    if (TYPE == 4) CONTR <- 1
    if (TYPE == 5) CONTR <- 1

    # Strain-conditional growth rates and Bmax-driven feedback rate
    KGDxS <- 1
    if (ID2 > 18) KGDxS <- 1 + kg_dyn_increase   # .mod L127
    BMAX <- 10^bmax_log10                          # .mod L97

    KGA  <- 0
    if (BACT == 2) KGA <- KGDxS * kga              # .mod L62, L169 (ETA dropped)
    KGP  <- 0
    if (BACT == 1) KGP <- KGDxS * kgp              # .mod L65, L170

    # Mutant-subpopulation growth rates (KGAM / KGPM)
    KGAFM <- kg_decrease_atcc
    KGPFM <- kg_decrease_aru
    KGAM <- KGDxS * kga * (1 - KGAFM)              # .mod L141
    KGPM <- KGDxS * kgp * (1 - KGPFM)              # .mod L142

    # Mutant fraction (MUTA = MUTP per .mod L148)
    MUTA <- 10^(-neg_log10_mutf_atcc)
    MUTP <- MUTA

    # Drug effect / Hill / interaction parameters strain- and condition-conditional
    EMAX1  <- 0
    if (COL == 1) {
      if (BACT == 2) EMAX1 <- emax_col_atcc        # ATCC + COL, .mod L70
    }
    SLOP1  <- 0
    if (COL == 1) {
      if (BACT == 1) SLOP1 <- slop1_col_aru        # ARU + COL, .mod L73
    }
    EMAX3  <- 0
    if (MERO == 1) {
      if (BACT == 1) EMAX3 <- emax_mero_aru        # ARU + MERO, .mod L79
    }
    EC51   <- 0
    if (COL == 1) {
      if (BACT == 2) EC51 <- ec51_col_atcc         # .mod L82
    }
    SLOPAM <- 0
    if (MERO == 1) {
      if (BACT == 2) SLOPAM <- slopam_mero_atcc    # .mod L85
    }
    EC53   <- 0
    if (MERO == 1) {
      if (BACT == 1) EC53 <- ec53_mero_aru         # .mod L88
    }
    GAMC <- 0
    if (COL == 1) GAMC <- gamc_col                 # .mod L90
    GAMA <- 0
    if (MERO == 1) {
      if (BACT == 2) GAMA <- gama_mero_atcc        # .mod L92
    }
    GAMP <- 0
    if (MERO == 1) {
      if (BACT == 1) GAMP <- gamp_mero_aru         # .mod L94 (NB: .mod uses TVGAMP = THETA(26))
    }

    # Mutant susceptibility shift (.mod L133-134)
    MSLOPAM <- 0
    if (MERO == 1) {
      if (BACT == 2) MSLOPAM <- slopam_mero_atcc / ((1 + concshift_mut)^gama_mero_atcc)
    }
    EC53M <- 0
    if (MERO == 1) {
      if (BACT == 1) EC53M <- ec53_mero_aru * (1 + concshift_mut)
    }

    # Adaptive-resistance KON/KOFF (colistin) and reserved meropenem slots
    KON51 <- 0
    if (COL == 1) {
      if (BACT == 2) KON51 <- kon51_col_atcc       # .mod L99
    }
    KON52 <- 0
    if (COL == 1) {
      if (BACT == 1) KON52 <- kon52_col_aru        # .mod L102
    }
    KOFFC <- 0
    if (COL == 1) KOFFC <- koffc_col               # .mod L105
    KONCA <- 0
    if (COL == 1) {
      if (BACT == 2) KONCA <- konca_col_atcc       # .mod L110
    }
    KONCP <- 0
    if (COL == 1) {
      if (BACT == 1) KONCP <- koncp_col_aru        # .mod L113
    }
    GAM3C <- 0
    if (COL == 1) GAM3C <- gam3c_col               # .mod L121
    GAM1  <- 0
    if (COL == 1) GAM1 <- gam1_col                 # .mod L123

    # Meropenem adaptive-resistance branch (declared but unused; .mod $DES leaves DADT(4)=DADT(5)=0)
    KOFFM <- 0
    if (MERO == 1) KOFFM <- koffm_mero
    KONMA <- 0
    if (MERO == 1) {
      if (BACT == 2) KONMA <- konma_mero_atcc
    }
    KONMP <- 0
    if (MERO == 1) {
      if (BACT == 1) KONMP <- konmp_mero_aru
    }

    # =====================================================================
    # KE2 — colistin in vitro elimination rate (1/h).
    # The .mod ($PK lines 199-276) hard-codes a per-experiment KE2 lookup
    # keyed on the experimental-well initial colistin concentration Ccol
    # with a TIME break at 8 h. The cascade is reproduced verbatim. Any
    # Ccol value not in the cascade leaves KE2 = 0 (matching .mod default).
    # =====================================================================
    KE1 <- 0.02   # Meropenem in vitro elimination rate, .mod L197 (constant, "value obtained from exp")

    KE2 <- 0
    # ATCC + colistin alone (.mod L200-213)
    if (Ccol == 0.04162 && t < 8)  KE2 <- 0.059
    if (Ccol == 0.04162 && t >= 8) KE2 <- 0.003
    if (Ccol == 0.1877  && t < 8)  KE2 <- 0.079
    if (Ccol == 0.1877  && t >= 8) KE2 <- 0.021
    if (Ccol == 0.2960  && t < 8)  KE2 <- 0.042
    if (Ccol == 0.2960  && t >= 8) KE2 <- 0.042
    if (Ccol == 0.7753  && t < 8)  KE2 <- 0.003
    if (Ccol == 0.7753  && t >= 8) KE2 <- 0.017
    if (Ccol == 2.182   && t < 8)  KE2 <- 0.038
    if (Ccol == 2.182   && t >= 8) KE2 <- 0.015
    if (Ccol == 4.305   && t < 8)  KE2 <- 0.0076
    if (Ccol == 4.305   && t >= 8) KE2 <- 0.00061
    if (Ccol == 12.024  && t < 8)  KE2 <- 0.012
    if (Ccol == 12.024  && t >= 8) KE2 <- 0.0033
    # ARU552 + colistin alone (.mod L214-227)
    if (Ccol == 0.08091 && t < 8)  KE2 <- 0.16
    if (Ccol == 0.08091 && t >= 8) KE2 <- 0.039
    if (Ccol == 0.2483  && t < 8)  KE2 <- 0.031
    if (Ccol == 0.2483  && t >= 8) KE2 <- 0.070
    if (Ccol == 0.9582  && t < 8)  KE2 <- 0.049
    if (Ccol == 0.9582  && t >= 8) KE2 <- 0.093
    if (Ccol == 2.672   && t < 8)  KE2 <- 0.037
    if (Ccol == 2.672   && t >= 8) KE2 <- 0.056
    if (Ccol == 5.113   && t < 8)  KE2 <- 0.027
    if (Ccol == 5.113   && t >= 8) KE2 <- 0.022
    if (Ccol == 11.471  && t < 8)  KE2 <- 0.011
    if (Ccol == 11.471  && t >= 8) KE2 <- 0.0058
    if (Ccol == 24.099  && t < 8)  KE2 <- 0.0078
    if (Ccol == 24.099  && t >= 8) KE2 <- 0.0021
    # Combination experiments (.mod L229-276). CA025 / CA026 grouping per .mod L229-238.
    CA025 <- 0
    if (Ccol == 0.205) CA025 <- 1
    if (Ccol == 0.199) CA025 <- 1
    if (Ccol == 0.195) CA025 <- 1
    CA026 <- 0
    if (Ccol == 0.211) CA026 <- 1
    if (Ccol == 0.223) CA026 <- 1
    if (Ccol == 0.187) CA026 <- 1
    if (CA025 == 1 && t < 8)  KE2 <- 0.196
    if (CA025 == 1 && t >= 8) KE2 <- 0.0534
    if (CA026 == 1 && t < 8)  KE2 <- 0.196
    if (CA026 == 1 && t >= 8) KE2 <- 0.0534
    # Intended 0.5 mg/L combined (.mod L239-250)
    if (Ccol == 0.431 && t < 8)  KE2 <- 0.111
    if (Ccol == 0.431 && t >= 8) KE2 <- 0.126
    if (Ccol == 0.411 && t < 8)  KE2 <- 0.111
    if (Ccol == 0.411 && t >= 8) KE2 <- 0.126
    if (Ccol == 0.408 && t < 8)  KE2 <- 0.111
    if (Ccol == 0.408 && t >= 8) KE2 <- 0.126
    if (Ccol == 0.482 && t < 8)  KE2 <- 0.111
    if (Ccol == 0.482 && t >= 8) KE2 <- 0.126
    if (Ccol == 0.479 && t < 8)  KE2 <- 0.111
    if (Ccol == 0.479 && t >= 8) KE2 <- 0.126
    if (Ccol == 0.465 && t < 8)  KE2 <- 0.111
    if (Ccol == 0.465 && t >= 8) KE2 <- 0.126
    # Intended 1 mg/L combined (.mod L251-254)
    if (Ccol == 0.673 && t < 8)  KE2 <- 0.0503
    if (Ccol == 0.673 && t >= 8) KE2 <- 0.0143
    if (Ccol == 0.71  && t < 8)  KE2 <- 0.0503
    if (Ccol == 0.71  && t >= 8) KE2 <- 0.0143
    # Intended 2 mg/L combined (.mod L255-264)
    if (Ccol == 1.467 && t < 8)  KE2 <- 0.0330
    if (Ccol == 1.467 && t >= 8) KE2 <- 0.0367
    if (Ccol == 1.234 && t < 8)  KE2 <- 0.0330
    if (Ccol == 1.234 && t >= 8) KE2 <- 0.0367
    if (Ccol == 1.433 && t < 8)  KE2 <- 0.0330
    if (Ccol == 1.433 && t >= 8) KE2 <- 0.0367
    if (Ccol == 1.657 && t < 8)  KE2 <- 0.0330
    if (Ccol == 1.657 && t >= 8) KE2 <- 0.0367
    if (Ccol == 1.418 && t < 8)  KE2 <- 0.0330
    if (Ccol == 1.418 && t >= 8) KE2 <- 0.0367
    # Intended 4 mg/L combined (.mod L265-272)
    if (Ccol == 3.523 && t < 8)  KE2 <- 0.0351
    if (Ccol == 3.523 && t >= 8) KE2 <- 0.0245
    if (Ccol == 3.447 && t < 8)  KE2 <- 0.0351
    if (Ccol == 3.447 && t >= 8) KE2 <- 0.0245
    if (Ccol == 3.288 && t < 8)  KE2 <- 0.0351
    if (Ccol == 3.288 && t >= 8) KE2 <- 0.0245
    if (Ccol == 3.031 && t < 8)  KE2 <- 0.0351
    if (Ccol == 3.031 && t >= 8) KE2 <- 0.0245
    # Intended 8 mg/L combined (.mod L273-276)
    if (Ccol == 7.669 && t < 8)  KE2 <- 0.00226
    if (Ccol == 7.669 && t >= 8) KE2 <- 0.0001
    if (Ccol == 6.622 && t < 8)  KE2 <- 0.00226
    if (Ccol == 6.622 && t >= 8) KE2 <- 0.0001

    # =====================================================================
    # Initial conditions
    # The dose record (EVID = 1, AMT = inoculum CFU/mL) deposits into the S
    # compartment; the mutant subpopulation is initialized via a separate
    # vignette-side dose to S_mut at AMT = inoculum * MUT (matching .mod
    # A_0(9) = NMUT = AMT * MUTA when MERO == 1). Drug-concentration
    # compartments are initialized from the Ccol / Cmer covariates.
    # =====================================================================
    mero(0)         <- Cmer    # A_0(3) = Cmer, .mod L160
    mero_ce(0)      <- 0       # A_0(4), .mod L161 (not used)
    mero_bindoff(0) <- 1       # A_0(5), .mod L162 (not used)
    col(0)          <- Ccol    # A_0(6) = Ccol, .mod L163
    col_ce(0)       <- 0       # A_0(7), .mod L164
    col_bindoff(0)  <- 1       # A_0(8), .mod L165
    R(0)            <- 0       # A_0(2), .mod L159
    R_mut(0)        <- 0       # A_0(10), .mod L167

    # =====================================================================
    # $DES — drug effects, conversion rates, and ODEs (.mod lines 280-437)
    # Reads compartment values: CA1 = col, CA2 = mero, CE1 = col_ce, CE2 = 0
    # =====================================================================
    AT <- S + R + S_mut + R_mut    # .mod L286 (total bacteria)

    # Single-drug effects on susceptible bacteria
    COLA <- 0
    if (COL == 1) {
      if (BACT == 2) COLA <- EMAX1 * (col)^GAMC / (EC51^GAMC + (col)^GAMC)   # .mod L290
    }
    COLP <- 0
    if (COL == 1) {
      if (BACT == 1) COLP <- SLOP1 * col^GAMC                                # .mod L293
    }
    MEROA <- 0
    if (MERO == 1) {
      if (BACT == 2) MEROA <- SLOPAM * mero^GAMA                             # .mod L296
    }
    MEROP <- 0
    if (MERO == 1) {
      if (BACT == 1) MEROP <- EMAX3 * (mero)^GAMP / (EC53^GAMP + (mero)^GAMP) # .mod L299
    }
    # Mutant variants
    MEROAM <- 0
    if (MERO == 1) {
      if (BACT == 2) MEROAM <- MSLOPAM * mero^GAMA                           # .mod L304
    }
    MEROPM <- 0
    if (MERO == 1) {
      if (BACT == 1) MEROPM <- EMAX3 * (mero)^GAMP / (EC53M^GAMP + (mero)^GAMP) # .mod L307
    }

    # Inhibition by colistin adaptive resistance (CE1 = col_ce); CE2 = 0 in $DES (.mod L282)
    COLA1 <- 0
    if (COL == 1) {
      if (BACT == 2) COLA1 <- COLA * (1 - col_ce^GAM1)                       # .mod L311
    }
    COLP1 <- 0
    if (COL == 1) {
      if (BACT == 1) COLP1 <- COLP * (1 - col_ce^GAM1)                       # .mod L314
    }
    MEROA1 <- 0
    if (MERO == 1) {
      if (BACT == 2) MEROA1 <- MEROA * (1 - 0)                               # CE2=0, .mod L317
    }
    MEROP1 <- 0
    if (MERO == 1) {
      if (BACT == 1) MEROP1 <- MEROP * (1 - 0)                               # CE2=0, .mod L320
    }
    MEROA1M <- 0
    if (MERO == 1) {
      if (BACT == 2) MEROA1M <- MEROAM * (1 - 0)                             # .mod L325
    }
    MEROP1M <- 0
    if (MERO == 1) {
      if (BACT == 1) MEROP1M <- MEROPM * (1 - 0)                             # .mod L328
    }

    # Combination drug effect (TYPE == 3). Guarded against zero-sum denominator
    # by adding a tiny floor; .mod relies on NONMEM's tolerance to silently
    # skip evaluation when MERO==0 / COL==0, but rxode2 always evaluates the
    # expression. With MERO==0 OR COL==0, MEROA1/COLA1 = 0 and the bracketed
    # ratio is 0/0; we use a small denominator floor matching standard
    # numerical practice (effect is gated by the outer if(TYPE==3) so the
    # value is discarded when not in combination mode).
    DRUGA <- 0
    DRUGB <- 0
    DRUGAM <- 0
    DRUGBM <- 0
    if (TYPE == 3) {
      if (BACT == 2) {
        denom_a <- COLA1 + MEROA1 + 1e-30
        DRUGA <- COLA1 * (1 + MEROA1 / denom_a)^inter_atcc +
                 MEROA1 * (1 + COLA1 / denom_a)^inter_atcc                   # .mod L333
        denom_am <- COLA1 + MEROA1M + 1e-30
        DRUGAM <- COLA1 * (1 + MEROA1M / denom_am)^inter_atcc +
                  MEROA1M * (1 + COLA1 / denom_am)^inter_atcc                # .mod L349
      }
      if (BACT == 1) {
        denom_b <- COLP1 + MEROP1 + 1e-30
        DRUGB <- COLP1 * (1 + MEROP1 / denom_b)^inter_aru +
                 MEROP1 * (1 + COLP1 / denom_b)^inter_aru                    # .mod L341
        denom_bm <- COLP1 + MEROP1M + 1e-30
        DRUGBM <- COLP1 * (1 + MEROP1M / denom_bm)^inter_aru +
                  MEROP1M * (1 + COLP1 / denom_bm)^inter_aru                 # .mod L356
      }
    }

    # Conversion rate from S -> R driven by total bacterial mass (.mod L362-373)
    FEEDA <- 0
    if (BACT == 2) FEEDA <- (KGA - kk) / BMAX
    FEEDP <- 0
    if (BACT == 1) FEEDP <- (KGP - kk) / BMAX
    SRA <- 0
    if (BACT == 2) SRA <- FEEDA * AT
    SRP <- 0
    if (BACT == 1) SRP <- FEEDP * AT

    # Colistin-resistance development rates (.mod L376-381)
    KON1CA <- 0
    if (COL == 1) {
      if (BACT == 2) KON1CA <- (KONCA * col^GAM3C) / (col^GAM3C + KON51^GAM3C)
    }
    KON1CP <- 0
    if (COL == 1) {
      if (BACT == 1) KON1CP <- (KONCP * col^GAM3C) / (col^GAM3C + KON52^GAM3C)
    }

    # =====================================================================
    # ODEs (.mod $DES lines 383-437)
    # =====================================================================

    # Susceptible bacteria S (.mod L383-391)
    dSdt <- 0
    if (COL == 1) {
      if (BACT == 2) dSdt <- KGA * S - (kk + COLA1) * S - SRA * S
      if (BACT == 1) dSdt <- KGP * S - (kk + COLP1) * S - SRP * S
    }
    if (MERO == 1) {
      if (BACT == 2) dSdt <- KGA * S - (kk + MEROA1) * S - SRA * S
      if (BACT == 1) dSdt <- KGP * S - (kk + MEROP1) * S - SRP * S
    }
    if (TYPE == 3) {
      if (BACT == 2) dSdt <- KGA * S - (kk + DRUGA) * S - SRA * S
      if (BACT == 1) dSdt <- KGP * S - (kk + DRUGB) * S - SRP * S
    }
    if (CONTR == 1) {
      if (BACT == 2) dSdt <- KGA * S - kk * S - SRA * S
      if (BACT == 1) dSdt <- KGP * S - kk * S - SRP * S
    }
    d/dt(S) <- dSdt

    # Resting bacteria R (.mod L393-395)
    dRdt <- 0
    if (BACT == 2) dRdt <- -kk * R + SRA * S
    if (BACT == 1) dRdt <- -kk * R + SRP * S
    d/dt(R) <- dRdt

    # Meropenem concentration (.mod L397-398)
    dmerodt <- 0
    if (MERO == 1) dmerodt <- -(KE1 + DILmer) * mero
    d/dt(mero) <- dmerodt

    # Meropenem adaptive-resistance compartments — both unused per .mod L401-408
    d/dt(mero_ce)      <- 0
    d/dt(mero_bindoff) <- 0

    # Colistin concentration (.mod L410-411)
    dcoldt <- 0
    if (COL == 1) dcoldt <- -(KE2 + DILcol) * col
    d/dt(col) <- dcoldt

    # Colistin adaptive-resistance ON / OFF (.mod L413-421)
    dcoldcedt <- 0
    if (COL == 1) {
      if (BACT == 2) dcoldcedt <- KON1CA * col_bindoff - KOFFC * col_ce
      if (BACT == 1) dcoldcedt <- KON1CP * col_bindoff - KOFFC * col_ce
    }
    d/dt(col_ce) <- dcoldcedt
    dcolboffdt <- 0
    if (COL == 1) {
      if (BACT == 2) dcolboffdt <- KOFFC * col_ce - KON1CA * col_bindoff
      if (BACT == 1) dcolboffdt <- KOFFC * col_ce - KON1CP * col_bindoff
    }
    d/dt(col_bindoff) <- dcolboffdt

    # Mutant susceptible S_mut (.mod L424-432)
    dSmutdt <- 0
    if (COL == 1) {
      if (BACT == 2) dSmutdt <- KGAM * S_mut - (kk + COLA1) * S_mut - SRA * S_mut
      if (BACT == 1) dSmutdt <- KGPM * S_mut - (kk + COLP1) * S_mut - SRP * S_mut
    }
    if (MERO == 1) {
      if (BACT == 2) dSmutdt <- KGAM * S_mut - (kk + MEROA1M) * S_mut - SRA * S_mut
      if (BACT == 1) dSmutdt <- KGPM * S_mut - (kk + MEROP1M) * S_mut - SRP * S_mut
    }
    if (TYPE == 3) {
      if (BACT == 2) dSmutdt <- KGAM * S_mut - (kk + DRUGAM) * S_mut - SRA * S_mut
      if (BACT == 1) dSmutdt <- KGPM * S_mut - (kk + DRUGBM) * S_mut - SRP * S_mut
    }
    if (CONTR == 1) {
      if (BACT == 2) dSmutdt <- KGAM * S_mut - kk * S_mut - SRA * S_mut
      if (BACT == 1) dSmutdt <- KGPM * S_mut - kk * S_mut - SRP * S_mut
    }
    d/dt(S_mut) <- dSmutdt

    # Mutant resting R_mut (.mod L434-437)
    dRmutdt <- 0
    if (BACT == 2) dRmutdt <- -kk * R_mut + SRA * S_mut
    if (BACT == 1) dRmutdt <- -kk * R_mut + SRP * S_mut
    d/dt(R_mut) <- dRmutdt

    # =====================================================================
    # Observation: log(total bacteria), with a 1e-5 floor to mirror the .mod
    # ATOT clip ($ERROR lines 452-454).
    # =====================================================================
    ATOT <- S + R + S_mut + R_mut
    if (ATOT < 1e-5) ATOT <- 1e-5
    logBact <- log(ATOT)
    logBact ~ add(addSd_logBact)
  })
}
