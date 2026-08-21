HernandezLozano_2025_apramycin_invitro <- function() {
  description <- "In vitro (Escherichia coli EN591 and ATCC 700336). Semi-mechanistic time-kill pharmacodynamic model of the aminoglycoside apramycin against two Escherichia coli urinary isolates at pH 7.4 and pH 6. Two bacterial subpopulations (1 = main, apramycin-susceptible; 2 = a subpopulation with decreased susceptibility) each carry two states: a growing, drug-susceptible state S and a dormant, drug-insusceptible state D that bacteria enter as a response to high population densities at rate ksr = (kg - kd) * Btot / Bmax. Apramycin adds to the death rate of the S states through a power model normalized to the strain-and-pH-specific MIC, kdrug = Slope * (Cu/MIC)^gamma with gamma fixed to 1, and simultaneously drives the transfer of bacteria from subpopulation 1 to subpopulation 2 at rate kada * (Cu/MIC), so no pre-existing resistant fraction is assumed in the inoculum. MIC-normalization makes one parameter set describe ATCC 700336 at both pH levels; EN591 needs separate drug-effect parameters at pH 6 and pH 7.4. Sibling models: HernandezLozano_2025_apramycin_mouse (the same PD structure re-estimated on mouse kidney and bladder cfu data and driven by a mouse subcutaneous PK model) and HernandezLozano_2025_apramycin_human (the in vivo PD component driven by the human population PK model of Zhao 2022)."
  reference <- "Hernandez-Lozano I, Aranzana-Climent V, Cao S, Matias C, Hansen JU, Liepinsh E, Hughes D, Hobbie SN, Vingsbo Lundberg C, Friberg LE. Model-informed drug development for antimicrobials: translational pharmacokinetic-pharmacodynamic modelling of apramycin to facilitate prediction of efficacious dose in complicated urinary tract infections. J Antimicrob Chemother. 2025 Feb 3;80(2):302-311. doi:10.1093/jac/dkae409. PMID: 39545353. PMCID: PMC11695905. In vitro PD structure: Materials and methods, 'PKPD modelling', plus the schematic in Figure 1. Parameter estimates and 95% CIs: Table 1, section 'In vitro PD parameters'. MIC values: Results, 'In vitro time-kill curves and PD modelling'. The natural bacterial death rate kd was fixed to 0.179 per hour from Nielsen EI, Cars O, Friberg LE, Antimicrob Agents Chemother 2011;55:4619-30 (doi:10.1128/AAC.00182-11)."
  vignette <- "HernandezLozano_2025_apramycin"
  units <- list(time = "h", dosing = "mg/L (bath concentration)", concentration = "mg/L")

  # Declared explicitly because the dosing target is neither `depot` nor
  # `central`: this is a static in vitro system, so the "dose" is the bath
  # concentration placed directly into the `apramycin` state at time 0.
  dosing <- c("apramycin")

  # Bacterial states follow the notation of Figure 1: subpopulation 1 is the
  # main, apramycin-susceptible population and subpopulation 2 the
  # less-susceptible ("resistant") one; within each, S is the growing
  # drug-susceptible state and D the dormant drug-insusceptible state. The
  # names do not match the canonical bact_<phenotype><lifecycle-digit> scheme
  # in compartment-names.md because the trailing digit there indexes the
  # Bulitta / Wicha replication life cycle, which is a different biological
  # distinction from this model's growing-versus-dormant states; the sibling
  # Nielsen-lineage model Zhao_2024_ciprofloxacin_colistin_invitro.R declares
  # its bacterial states the same way. `apramycin` holds the static bath
  # concentration in mg/L, following the in-vitro convention established by
  # Nielsen_2007_semimechanistic_antibiotic_pd.R.
  paper_specific_compartments <- c(
    "apramycin",
    "bact_s1", "bact_d1", "bact_s2", "bact_d2"
  )

  compartmentData <- list(
    apramycin = list(analyte = "apramycin", units = "mg/L", specimen = "administration site", verified = TRUE),
    bact_s1 = list(analyte = "Escherichia coli, main apramycin-susceptible subpopulation, growing drug-susceptible (S) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_d1 = list(analyte = "Escherichia coli, main apramycin-susceptible subpopulation, dormant drug-insusceptible (D) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_s2 = list(analyte = "Escherichia coli, subpopulation with decreased apramycin susceptibility, growing drug-susceptible (S) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_d2 = list(analyte = "Escherichia coli, subpopulation with decreased apramycin susceptibility, dormant drug-insusceptible (D) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    BACT = list(
      description = "Bacterial-strain indicator selecting the growth rate constant and, jointly with PH_MEDIUM, the MIC and the three drug-effect parameters (SlopeS, SlopeR, kada) of the in vitro time-kill fit.",
      units = "(categorical)",
      type = "categorical",
      source_name = "STRAIN",
      reference_category = NULL,
      notes = paste(
        "591 = Escherichia coli EN591, an MDR rmtB clinical isolate; MIC 8 mg/L at pH 7.4 and 32 mg/L at pH 6.",
        "700336 = Escherichia coli ATCC 700336 (EN1085), a trimethoprim/sulfamethoxazole-resistant urinary tract infection isolate; MIC 4 mg/L at pH 7.4 and 16 mg/L at pH 6.",
        "There is no reference level: each level selects its own parameter set rather than contributing a contrast against a baseline.",
        "The growth rate constant kg is strain-specific but pH-independent, so a strain shares it across both pH levels. The same codes are used by the sibling in vivo files HernandezLozano_2025_apramycin_mouse.R and HernandezLozano_2025_apramycin_human.R."
      )
    ),
    PH_MEDIUM = list(
      description = "pH of the Mueller-Hinton II growth medium, selecting the MIC for both strains and, for EN591 only, the drug-effect parameter set (SlopeS, SlopeR, kada).",
      units = "pH units",
      type = "continuous",
      source_name = "PH",
      reference_category = NULL,
      notes = paste(
        "The time-kill experiments were run at two buffered pH levels, 7.4 (the conventional neutral condition) and 6 (chosen because urine is acidic); the medium was buffered so pH stayed stable across the 28 h experiment.",
        "Apramycin's MIC is four times higher at pH 6 than at pH 7.4 for both strains (Results, 'In vitro time-kill curves and PD modelling'; Becker 2021, doi:10.1016/j.ebiom.2021.103652).",
        "For ATCC 700336, MIC-normalization of the drug effect let a single parameter set describe both pH levels 'without significantly changing the model fit'; for EN591, 'differences in bacterial regrowth at 4 x MIC between pH 6 and pH 7.4 conditions resulted in different drug effect parameters' (Results). PH_MEDIUM therefore switches SlopeS, SlopeR and kada for EN591 but only the MIC for ATCC 700336.",
        "Only the two studied levels are characterised. The model dispatches on a midpoint threshold of 6.7 rather than interpolating, because no intermediate pH was studied; a value below 6.7 is treated as the pH 6 condition and a value at or above it as the pH 7.4 condition.",
        "The sibling in vivo files omit this covariate: there pH is a property of the organ (kidney assumed 7.4, bladder assumed 6) rather than an experimental factor, and is baked into which per-site MIC applies."
      )
    )
  )

  population <- list(
    species = "in vitro (Escherichia coli EN591 and ATCC 700336)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    organism = "Escherichia coli EN591 (MDR rmtB isolate) and Escherichia coli ATCC 700336 / EN1085 (trimethoprim/sulfamethoxazole-resistant urinary tract infection isolate)",
    system = "Static time-kill experiments in 2 mL propylene tubes of prewarmed Mueller-Hinton II broth, incubated at 37 C, with viable counts at 0, 1, 2, 4, 6, 8, 10, 24 and 28 h",
    medium = "Mueller-Hinton II broth, pH adjusted and buffered to pH 7.4 or pH 6 to keep the pH stable throughout the experiment",
    temperature = "37 C",
    duration = "28 h",
    starting_inoculum = "10^5 to 10^6 CFU aliquots inoculated; the model estimated an initial density of 10^5.53 CFU/mL",
    limit_of_detection = "10 CFU/mL",
    mic_values = c(
      `EN591 pH 7.4` = "8 mg/L",
      `EN591 pH 6` = "32 mg/L",
      `ATCC 700336 pH 7.4` = "4 mg/L",
      `ATCC 700336 pH 6` = "16 mg/L"
    ),
    concentration_range = "0.25 to 8 x MIC",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "The acidic pH arm was included because urine is more acidic than blood and the bactericidal activity of apramycin changes with pH; the MIC is four times higher at pH 6 than at pH 7.4 for both strains (Becker 2021, doi:10.1016/j.ebiom.2021.103652, reference 22 of the source paper).",
      "Model development, VPCs and the sampling-importance resampling that produced the 95% CIs in Table 1 were done in NONMEM 7.5.0 with PsN 3.5.0.",
      "An adaptive-resistance mechanism and an Emax concentration-effect relationship were both tested; the retained model uses the power drug-effect model and the kada susceptible-to-resistant transfer encoded here."
    )
  )

  ini({
    # ---- Simulation design inputs ---------------------------------------
    # The four measured MICs are exposed as fixed ini() entries, matching the
    # sibling in vivo files, so a user can re-point the model at a different
    # isolate with rxSolve(params = c(micEn591ph74 = ...)). Which one applies
    # is selected in model() from BACT and PH_MEDIUM.
    micEn591ph74 <- fixed(8); label("MIC of apramycin for EN591 at pH 7.4 (mg/L)")                            # Results: MIC_EN591 = 8 mg/L at pH 7.4
    micEn591ph6 <- fixed(32); label("MIC of apramycin for EN591 at pH 6 (mg/L)")                              # Results: MIC_EN591 = 32 mg/L at pH 6
    micAtccph74 <- fixed(4); label("MIC of apramycin for ATCC 700336 at pH 7.4 (mg/L)")                       # Results: MIC_ATCC700336 = 4 mg/L at pH 7.4
    micAtccph6 <- fixed(16); label("MIC of apramycin for ATCC 700336 at pH 6 (mg/L)")                         # Results: MIC_ATCC700336 = 16 mg/L at pH 6
    fu <- fixed(1); label("Fraction of apramycin unbound in the growth medium (unitless)")                    # Protein-free Mueller-Hinton II broth: the total bath concentration IS the free concentration. Exposed so the same PD block can be driven by an unbound plasma concentration in the sibling in vivo files

    # ---- Bacterial growth model (Table 1, In vitro PD parameters) --------
    lkgEn591 <- log(1.64); label("Log growth rate constant of the susceptible subpopulation, EN591 (1/h)")    # Table 1: kg = 1.64 (95% CI 1.54-1.77) for EN591
    lkgAtcc <- log(2.27); label("Log growth rate constant of the susceptible subpopulation, ATCC 700336 (1/h)") # Table 1: kg = 2.27 (95% CI 2.05-2.59) for ATCC700336
    redukg <- 0.63; label("Fractional reduction in growth rate for the resistant subpopulation (%)")          # Table 1: Redukg = 0.63 % (95% CI 0.05-5.51), shared by both strains
    lkd <- fixed(log(0.179)); label("Log natural bacterial death rate constant (1/h)")                        # Table 1: kd = 0.179, Fixed. Methods: "set to a fixed value of 0.179 h-1 based on previously published data" (Nielsen 2011, reference 23)
    lbmax <- log(10^9.18); label("Log maximum bacterial density (CFU/mL)")                                    # Table 1: Bmax = 9.18 log10 CFU/mL (95% CI 8.97-9.38); carried on the natural-log scale per the library naming convention
    linoc <- log(10^5.53); label("Log initial bacterial density (CFU/mL)")                                    # Table 1: Inoc = 5.53 log10 CFU/mL (95% CI 5.44-5.64)

    # ---- Drug effect (Table 1, In vitro PD parameters) -------------------
    gam <- fixed(1); label("Power on the MIC-normalized drug effect (unitless)")                              # Table 1: gamma = 1, Fixed
    slopeSEn591ph6 <- 1.31; label("Apramycin effect on the susceptible subpopulation, EN591 at pH 6 (1/h)")    # Table 1: SlopeS = 1.31 (95% CI 1.18-1.45), EN591 pH 6
    slopeSEn591ph74 <- 2.37; label("Apramycin effect on the susceptible subpopulation, EN591 at pH 7.4 (1/h)") # Table 1: SlopeS = 2.37 (95% CI 1.89-2.85), EN591 pH 7
    slopeSAtcc <- 2.99; label("Apramycin effect on the susceptible subpopulation, ATCC 700336 (1/h)")          # Table 1: SlopeS = 2.99 (95% CI 2.38-3.80), ATCC700336 (both pH)
    slopeREn591ph6 <- 0.385; label("Apramycin effect on the resistant subpopulation, EN591 at pH 6 (1/h)")     # Table 1: SlopeR = 0.385 (95% CI 0.353-0.420), EN591 pH 6
    slopeREn591ph74 <- 0.254; label("Apramycin effect on the resistant subpopulation, EN591 at pH 7.4 (1/h)")  # Table 1: SlopeR = 0.254 (95% CI 0.232-0.280), EN591 pH 7
    slopeRAtcc <- 0.776; label("Apramycin effect on the resistant subpopulation, ATCC 700336 (1/h)")           # Table 1: SlopeR = 0.776 (95% CI 0.604-0.967), ATCC700336 (both pH)

    # kada is tabulated as "kada x 1000", so the value carried here is the
    # printed number divided by 1000. That reading is confirmed by the
    # Results text, which states that kada was "notably higher in vivo than
    # in vitro": the in vivo estimates are 0.031 (kidney) and 1.35 (bladder),
    # which exceed 8.0e-5 but would be LOWER than the printed 0.080.
    kadaEn591ph6 <- 0.080 / 1000; label("Drug-driven transfer rate constant from the susceptible to the resistant subpopulation, EN591 at pH 6 (1/h)")    # Table 1: kada x 1000 = 0.080 (95% CI 0.033-0.153), EN591 pH 6
    kadaEn591ph74 <- 0.028 / 1000; label("Drug-driven transfer rate constant from the susceptible to the resistant subpopulation, EN591 at pH 7.4 (1/h)") # Table 1: kada x 1000 = 0.028 (95% CI 0.012-0.063), EN591 pH 7
    kadaAtcc <- 0.210 / 1000; label("Drug-driven transfer rate constant from the susceptible to the resistant subpopulation, ATCC 700336 (1/h)")          # Table 1: kada x 1000 = 0.210 (95% CI 0.013-0.581), ATCC700336 (both pH)

    # ---- Residual error ---------------------------------------------------
    # "For residual unexplained variability, additive error models were used"
    # (Materials and methods, Data analysis and software), but no sigma value
    # is reported anywhere in the paper or the supplement. Held at 0 rather
    # than invented; see the vignette Errata.
    addSd <- fixed(0); label("Additive residual SD on the log10 bacterial count (log10 CFU/mL); magnitude not reported")  # Methods: additive residual error model; no value reported in Table 1 or the supplement
  })

  model({
    # ---- 1. Back-transform and select the strain-and-pH parameter set ----
    kd <- exp(lkd)
    bmax <- exp(lbmax)
    inoc <- exp(linoc)

    # Table 1 indexes the in vitro drug-effect parameters by strain AND by
    # growth-medium pH, so the two covariates are switched on jointly. Only
    # the two studied pH levels (6 and 7.4) are characterised; the midpoint
    # 6.7 dispatches between them rather than interpolating. EN591 takes
    # pH-specific SlopeS / SlopeR / kada; for ATCC 700336 the paper found
    # MIC-normalization sufficient, so only its MIC changes with pH.
    if (BACT == 591) {
      kg1 <- exp(lkgEn591)
      if (PH_MEDIUM < 6.7) {
        mic <- micEn591ph6
        slopeS <- slopeSEn591ph6
        slopeR <- slopeREn591ph6
        kada <- kadaEn591ph6
      } else {
        mic <- micEn591ph74
        slopeS <- slopeSEn591ph74
        slopeR <- slopeREn591ph74
        kada <- kadaEn591ph74
      }
    } else {
      kg1 <- exp(lkgAtcc)
      slopeS <- slopeSAtcc
      slopeR <- slopeRAtcc
      kada <- kadaAtcc
      if (PH_MEDIUM < 6.7) {
        mic <- micAtccph6
      } else {
        mic <- micAtccph74
      }
    }

    # Growth of the resistant subpopulation is reduced by Redukg percent
    # (Table 1 unit column: %).
    kg2 <- kg1 * (1 - redukg / 100)

    # ---- 2. Drug effect ---------------------------------------------------
    # The bath concentration is static in the source experiments; the state
    # holds mg/L directly, so no volume term is needed.
    conc <- fu * apramycin
    cmic <- conc / mic
    kdrug1 <- slopeS * cmic^gam
    kdrug2 <- slopeR * cmic^gam
    # RECONSTRUCTION, NOT A PRINTED EQUATION. The source paper contains no
    # numbered equations; kada appears only in the Methods sentence "A rate
    # constant (k_ada), describing the drug-driven transfer from the
    # susceptible to the resistant population, was included in the model
    # rather than assuming initial pre-existing percentages of each
    # subpopulation" and as a row of Table 1. Figure 1 draws no arrow for it
    # and its legend does not list it. The form below was inferred on
    # dimensional grounds: Table 1 gives kada in 1/h, so the multiplier must
    # be dimensionless, and the MIC-normalized concentration cmic is the only
    # such quantity in the model (a raw-concentration form would leave kada in
    # L/(mg*h), contradicting the table, and a plain first-order form would
    # not be "drug-driven" and would let resistance emerge in the untreated
    # growth-control arms). Four checks the reconstruction was NOT built to
    # satisfy all pass: growth controls plateau at the tabulated Bmax; the
    # in vivo kidney control plateaus at the tabulated Bmax,k; >= 10 mg/kg BID
    # achieves stasis as the Results state; and -- the sharp one -- EN591 at
    # 4 x MIC regrows by 28 h at pH 7.4 but is held down at pH 6, which is
    # precisely the observation the paper cites as its reason for splitting
    # the pH-specific drug-effect parameters. See the vignette Errata.
    kada1to2 <- kada * cmic

    # ---- 3. Density-dependent transfer into the dormant state -------------
    btot <- bact_s1 + bact_d1 + bact_s2 + bact_d2
    ksr1 <- (kg1 - kd) * btot / bmax
    ksr2 <- (kg2 - kd) * btot / bmax

    # ---- 4. Drug in the growth medium -------------------------------------
    # Static exposure: the pH-buffered medium kept apramycin stable over the
    # 28 h experiment, and no degradation rate is reported.
    d/dt(apramycin) <- 0

    # ---- 5. Bacterial system (Figure 1) -----------------------------------
    # Drug acts only on the growing S states; the dormant D states are
    # unsusceptible and only die at kd. There is no back-transfer D -> S.
    d/dt(bact_s1) <- kg1 * bact_s1 - (kd + kdrug1) * bact_s1 -
      ksr1 * bact_s1 - kada1to2 * bact_s1
    d/dt(bact_d1) <- ksr1 * bact_s1 - kd * bact_d1
    d/dt(bact_s2) <- kg2 * bact_s2 - (kd + kdrug2) * bact_s2 -
      ksr2 * bact_s2 + kada1to2 * bact_s1
    d/dt(bact_d2) <- ksr2 * bact_s2 - kd * bact_d2

    # ---- 6. Initial conditions --------------------------------------------
    # All bacteria start in the growing state of the main subpopulation.
    # Methods: kada was included "rather than assuming initial pre-existing
    # percentages of each subpopulation", so subpopulation 2 starts empty.
    bact_s1(0) <- inoc

    # ---- 7. Observation ----------------------------------------------------
    Cc <- log10(btot + 1e-6)
    Cc ~ add(addSd)
  })
}
