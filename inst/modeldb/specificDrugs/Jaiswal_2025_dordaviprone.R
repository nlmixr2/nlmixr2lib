Jaiswal_2025_dordaviprone <- function() {
  description <- paste(
    "Two-compartment oral population pharmacokinetic reduction of the",
    "Simcyp minimal-PBPK-with-single-adjusting-compartment (SAC) model",
    "for the imipridone dordaviprone (ONC201) in healthy adults",
    "(Jaiswal 2025). The source model was built in the Simcyp",
    "Population-based Simulator (V21) with the mechanistic ADAM",
    "absorption model; neither the whole-body mass-balance equations nor",
    "the ADAM gut model are published, so the platform model itself",
    "cannot be encoded here. What IS fully reported is the dordaviprone",
    "compound layer, and it is sufficient to reconstruct the disposition",
    "as an ordinary compartmental model: first-order absorption into a",
    "depot, distribution between a systemic compartment and the SAC (the",
    "paper's k_in / k_out, encoded as the canonical k12 / k21), and",
    "concentration-dependent elimination from the systemic compartment.",
    "Dordaviprone binds saturably to alpha-1-acid glycoprotein, so the",
    "unbound fraction rises with concentration (Table 1) and apparent",
    "oral clearance rises with dose; that relationship is carried",
    "explicitly and reproduces the paper's own dose non-proportionality",
    "(the paper's simulated CL/F rises 15.9% from 125 to 625 mg; this",
    "model gives 16.6%) without any adjustment.",
    "Two parameters could NOT be transcribed and were calibrated against",
    "the paper's own simulated summary statistics, because the source",
    "reports no first-order ka (it used ADAM) and no gut availability:",
    "ka was calibrated to the simulated Tmax and the baseline",
    "bioavailability to the simulated Cmax, both at 625 mg. Both are",
    "corroborated independently -- the calibrated ka of 1.58 1/h is",
    "within 12% of the 1.80 1/h implied by the reported P_eff,man, and",
    "the calibrated bioavailability is the value that makes the reported",
    "itraconazole Cmax ratio physically attainable. Everything else is",
    "either a verbatim Table 1 input or an arithmetic consequence of one.",
    "The reduction reproduces the paper's simulated AUC to within 9.3%,",
    "Cmax to within 5.6% and Tmax to within 2.4% at the two doses that",
    "were NOT used for calibration (125 and 375 mg).",
    "Coadministered CYP3A4 modulators act through a single relative",
    "CYP3A4 activity term that scales gut and hepatic first-pass",
    "extraction and systemic clearance together; each modulator's",
    "activity was back-solved from its published AUC ratio alone, and",
    "the published Cmax ratio is then reproduced within 11% for five of",
    "the six modulators (rifampicin is over-predicted by 23%).",
    "This is a typical-value simulation model: the source reports no",
    "inter-individual variance components and no residual-error model,",
    "so there are no etas and propSd is fixed at zero. The perpetrator",
    "predictions that are part of the paper's contribution (dordaviprone",
    "as an inhibitor of CYP3A4, CYP2C8 and CYP2D6) depend on proprietary",
    "Simcyp substrate compound files and are NOT reproducible from this",
    "model; see the vignette for the full list of deviations.",
    sep = " "
  )
  reference <- paste(
    "Jaiswal S, Patel NK, Jones HM, McFeely S, Faison SL, Tippin T,",
    "Naderer O. (2025). Assessing Cytochrome P450 Drug Interaction Risk",
    "for Dordaviprone Using Physiologically Based Pharmacokinetic",
    "Modeling. CPT Pharmacometrics Syst Pharmacol 14(10):1695-1704.",
    "doi:10.1002/psp4.70093.",
    sep = " "
  )
  vignette <- "Jaiswal_2025_dordaviprone"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Jaiswal 2025 Figure 1 (the PBPK
  # modeling-strategy schematic) and Section 2.3, which describes the
  # minimal distribution model as a systemic compartment exchanging with
  # a single adjusting compartment.
  compartmentData <- list(
    depot = list(
      analyte = "dordaviprone", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "dordaviprone", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "dordaviprone", units = "mg",
      specimen = "tissue", verified = TRUE
    )
  )

  # The six coadministered CYP3A4 modulators the paper simulates. Each is
  # a mutually-exclusive binary flag; all zero = dordaviprone alone. The
  # paper reports one AUC geometric mean ratio per modulator, and each
  # coefficient below is the relative CYP3A4 activity back-solved from
  # that single published ratio (see ini()).
  covariateData <- list(
    CONMED_ITRACONAZOLE = list(
      description        = "Concomitant itraconazole (strong CYP3A4 inhibitor)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no CYP3A4 modulator coadministered",
      notes              = paste(
        "Jaiswal 2025 Table 2 / Table S6: multiple-dose itraconazole",
        "with a single 125 mg dordaviprone dose, clinical study 103.",
        "This is the only modulator arm with observed clinical data;",
        "the other five are model predictions."
      ),
      source_name        = "itraconazole"
    ),
    CONMED_ERYTHROMYCIN = list(
      description        = "Concomitant erythromycin (moderate CYP3A4 inhibitor)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no CYP3A4 modulator coadministered",
      notes              = "Jaiswal 2025 Figure 3: erythromycin 500 mg four times daily, simulated only.",
      source_name        = "erythromycin"
    ),
    CONMED_FLUCONAZOLE = list(
      description        = "Concomitant fluconazole (moderate CYP3A4 inhibitor)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no CYP3A4 modulator coadministered",
      notes              = "Jaiswal 2025 Figure 3: fluconazole 200 mg once daily, simulated only.",
      source_name        = "fluconazole"
    ),
    CONMED_CIMETIDINE = list(
      description        = "Concomitant cimetidine (weak CYP3A4 inhibitor)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no CYP3A4 modulator coadministered",
      notes              = "Jaiswal 2025 Figure 3: cimetidine 400 mg twice daily, simulated only.",
      source_name        = "cimetidine"
    ),
    CONMED_EFV = list(
      description        = "Concomitant efavirenz (moderate CYP3A4 inducer)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no CYP3A4 modulator coadministered",
      notes              = paste(
        "Jaiswal 2025 Figure 3: efavirenz 600 mg once daily, simulated",
        "only. Here the reference category is simply 'no efavirenz'",
        "rather than the alternative-antiretroviral-regimen reference",
        "used by the antiretroviral popPK models that share this column."
      ),
      source_name        = "efavirenz"
    ),
    CONMED_RIFAMPICIN = list(
      description        = "Concomitant rifampicin (strong CYP3A4 inducer)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no CYP3A4 modulator coadministered",
      notes              = "Jaiswal 2025 Figure 3: rifampicin 600 mg once daily, simulated only.",
      source_name        = "rifampicin"
    )
  )

  # Screened by the source model but not carried here.
  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. Jaiswal 2025 Table 1 expresses Vss and Vsac in",
        "L/kg, so the Simcyp model does scale distribution volume with",
        "weight, and the simulations matched each trial's observed",
        "demographics. This reduction fixes the 70 kg reference weight",
        "instead, because the systemic-compartment volume also subtracts",
        "a liver volume whose weight-scaling Simcyp does not publish."
      ),
      units = "kg", type = "continuous",
      notes = "Simcyp North European Caucasian population file; relationship not published."
    ),
    CONMED_PPI = list(
      description = paste(
        "Concomitant proton-pump inhibitor (rabeprazole 20 mg once",
        "daily). Jaiswal 2025 Table 2 / Table S8 simulated and observed",
        "this arm, and both showed no change in dordaviprone exposure",
        "(simulated AUC 18,081 vs 18,086 h*ng/mL without and with",
        "rabeprazole). The effect operates through the pH-solubility",
        "profile inside the ADAM absorption model, which is not",
        "reproducible here; because the published effect is null, the",
        "no-modulator model reproduces both arms."
      ),
      units = "(binary)", type = "binary",
      notes = "Null effect in the source; documented rather than encoded."
    ),
    FOOD = list(
      description = paste(
        "High-fat meal. Jaiswal 2025 Table 2 / Table S5 report a fed",
        "Cmax about 25% below fasted with essentially unchanged AUC",
        "(simulated fed AUC 18,515 vs fasted 18,474 h*ng/mL). The food",
        "effect arises from bile-micelle partitioning and gastric",
        "emptying inside the ADAM model, neither of which survives the",
        "reduction to a single first-order ka, so the fed arm is not",
        "encoded."
      ),
      units = "(binary)", type = "binary",
      notes = "Absorption-model-mediated; not reproducible from a first-order ka."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 15L,
    n_studies      = 4L,
    age_range      = "Healthy adult participants; individual ages not tabulated in the source.",
    weight_median  = "70 kg (reference weight used for the L/kg volume inputs)",
    disease_state  = paste(
      "Healthy participants. The paper argues the same magnitudes of",
      "change would apply to the target H3 K27M-mutant glioma",
      "population, but no patient data were used in model development",
      "or verification."
    ),
    dose_range     = "Single oral doses of 125, 375 and 625 mg; 625 mg is the planned therapeutic dose (once weekly).",
    regions        = "North European Caucasian virtual population (Simcyp default, Howgate 2006) with trial-matched demographics.",
    studies        = paste(
      "ONC201-101 Part A single ascending dose (n = 15 healthy",
      "participants) supplied the 625 mg profile used to fit the",
      "minimal-PBPK distribution parameters and the 125 / 375 mg",
      "profiles used for verification; Part B2 (n = 30) supplied the",
      "high-fat-meal arm. The [14C]-dordaviprone mass-balance study",
      "(n = 6) fixed f_a and the zero renal clearance. Clinical study",
      "103 (n = 18) is the itraconazole drug-drug-interaction study",
      "whose observed CL/F of 29 L/h anchors the retrograde clearance",
      "calculation and whose AUC and Cmax ratios calibrated fm CYP3A4.",
      "Clinical study 107 (n = 16) supplied the rabeprazole arm."
    ),
    notes          = paste(
      "n_subjects records the 15 participants of ONC201-101 Part A, the",
      "study the disposition parameters were fitted against; n_studies",
      "counts the four clinical studies listed above. This is a PBPK",
      "analysis rather than a population-PK fit, so there is no pooled",
      "analysis dataset and no estimated variance components. The",
      "simulated geometric CV% reported in Supporting Information",
      "Tables S1-S8 (about 62% on AUC and 45% on Cmax) is Simcyp",
      "virtual-population output driven by demographic and enzyme-",
      "abundance distributions that are not published, so it is not",
      "encoded as an omega here."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter below is fixed: this is a typical-value simulation
    # model. Values are Jaiswal 2025 Table 1 inputs, arithmetic
    # consequences of them, or (for the two absorption / first-pass
    # parameters that the source does not report at all) calibrated
    # against the paper's own simulated summary statistics. Each case is
    # labelled explicitly below.
    #
    # Shared quantities, all Jaiswal 2025 Table 1:
    #   MW    = 386.49 g/mol
    #   B:P   = 0.778           blood-to-plasma ratio
    #   f_a   = 0.8             fraction absorbed (Section 2.3)
    #   Vss   = 1.11 L/kg       predicted, Rodgers & Rowland "Method 3"
    #   Vsac  = 6.60e-05 L/kg   single adjusting compartment volume
    #   k_in  = 0.11 1/h, k_out = 0.15 1/h   fitted to the 625 mg profile
    #   CL/F  = 29 L/h          observed apparent oral clearance, study 103
    #   fm CYP3A4 = 0.8, fm CYP2D6 = 0.2
    #   P_eff,man = 3.74e-04 cm/s (Simcyp predicted)
    #
    # Un-printed standard physiological constants used in the derivations
    # below: body weight 70 kg, liver volume 1.65 L, hepatic blood flow
    # QH 90 L/h. Only QH enters anything DDI-related, and the
    # gut-versus-liver split it controls is insensitive to it: varying QH
    # by +/-20% moves the maximum attainable bioavailability increase by
    # only 3% (1.400 to 1.442).
    # ------------------------------------------------------------------

    # -- Absorption ----------------------------------------------------
    lka <- fixed(log(1.579834))
    label("First-order absorption rate constant ka (1/h)")
    # CALIBRATED, not transcribed. Jaiswal 2025 used the Simcyp ADAM
    # mechanistic absorption model and reports no first-order ka, so a
    # value must be chosen for the reduction. It was calibrated to
    # reproduce the simulated median Tmax of 1.20 h after 625 mg
    # (Supporting Information Table S1). Independent corroboration: the
    # reported P_eff,man of 3.74e-04 cm/s implies ka = 2 * P_eff / r with
    # the standard human small-intestinal radius r = 1.5 cm, i.e.
    # 2 * 3.74e-04 * 3600 / 1.5 = 1.795 1/h, which is within 12% of the
    # calibrated value. The precedent model Hanley_2024_brigatinib.R
    # carries a ka chosen the same way (its source likewise rejected the
    # permeability-predicted ka in favour of one recovering the observed
    # Tmax).

    # -- Distribution --------------------------------------------------
    # Systemic-compartment volume. Table 1 gives whole-body Vss =
    # 1.11 L/kg and the SAC volume Vsac = 6.60e-05 L/kg. In the Simcyp
    # minimal-PBPK layout the systemic compartment holds what is left
    # after the SAC and the liver:
    #   Vsys = (1.11 - 0.0000660) * 70 - 1.65 = 76.045 L
    # This model is written in APPARENT parameters (everything divided by
    # the baseline oral bioavailability F), so the volume carried here is
    #   vc = Vsys / F = 76.045 / 0.531058 = 143.196 L
    # See `eh` / `egut` below for how F = 0.531058 is obtained.
    lvc <- fixed(log(143.19592))
    label("Apparent central volume Vc/F (L)")

    lk12 <- fixed(log(0.11))
    label("Central-to-SAC transfer rate constant k12 (1/h)")
    # Jaiswal 2025 Table 1 k_in = 0.11 1/h, fitted by the authors to
    # recover the observed plasma profile in the itraconazole study. The
    # rate constants act on drug MASS, so k_in / k_out map exactly onto
    # the canonical k12 / k21 and Vsac itself never enters the plasma
    # prediction.

    lk21 <- fixed(log(0.15))
    label("SAC-to-central transfer rate constant k21 (1/h)")
    # Jaiswal 2025 Table 1 k_out = 0.15 1/h.

    # -- Elimination ---------------------------------------------------
    # Apparent unbound intrinsic clearance. Because dordaviprone binds
    # saturably to AAG, CL/F is proportional to the unbound fraction:
    # writing CL/F = CLint/F * fu and cancelling the hepatic availability
    # that appears in both CL and F gives an exactly linear dependence on
    # fu, with no residual dependence on QH. CLint/F was set so that the
    # model returns the reported CL/F of 29 L/h at a 125 mg dose
    # (Table 1, clinical study 103).
    lclint <- fixed(log(968.1837))
    label("Apparent unbound intrinsic clearance CLint/F (L/h)")

    fm_cyp3a4 <- fixed(0.8)
    label("Fraction of dordaviprone metabolism mediated by CYP3A4")
    # Jaiswal 2025 Table 1 fm CYP3A4 = 0.8 (the remaining 0.2 is CYP2D6).
    # The authors adjusted this in the retrograde calculation to recover
    # the observed itraconazole AUC and Cmax ratios. It sets a ceiling of
    # 1 / 0.2 = 5.0 on the attainable AUC ratio for a complete CYP3A4
    # inhibitor, which the reported strong-inhibitor value of 4.62 sits
    # just under -- an internal consistency check on the reported fm.

    # -- Saturable plasma protein binding ------------------------------
    # Jaiswal 2025 Section 2.3 states that saturable binding to AAG was
    # accounted for in the PBPK model, and Table 1 tabulates fu at five
    # total plasma concentrations:
    #   0.193 mg/L -> 0.0294      3.865 mg/L -> 0.0412
    #   0.386 mg/L -> 0.0303      7.730 mg/L -> 0.0492
    #   1.932 mg/L -> 0.0389
    # A saturable form fitted to those five printed points reproduces
    # every one of them to within 4%.
    # Canonical names ratified by operator sidecar (request-001 q1, answer B):
    # fumin / fumax / cup50, matching the earlier ruling on the Zhang FDA
    # review (oare_PMC11544005) for the same concept. The paper's own symbol
    # for all three is simply "fu" tabulated against concentration.
    fumin <- fixed(0.02827436)
    label("Lower-asymptote unbound fraction (low concentration)")
    fumax <- fixed(0.06279485)
    label("Upper-asymptote unbound fraction (saturating concentration)")
    cup50 <- fixed(5.337066)
    label("Total plasma concentration at half-maximal unbinding (mg/L)")

    # -- Baseline first-pass extraction --------------------------------
    # Two extraction ratios split the baseline oral bioavailability
    # F = f_a * (1 - egut) * (1 - eh) = 0.8 * 0.851 * 0.780 = 0.531058.
    # They are pinned by two simultaneous constraints, with no freedom
    # left over once F is known:
    #   (i)  the reported CL/F = 29 L/h with f_a = 0.8, B:P = 0.778 and
    #        QH = 90 L/h fixes the product that appears in the
    #        well-stirred retrograde calculation;
    #   (ii) F itself is CALIBRATED to the simulated Cmax of 2810 ng/mL
    #        after 625 mg (Supporting Information Table S1), because the
    #        source reports no gut availability and ADAM's gut model is
    #        not published.
    # Corroboration for (ii): the reported itraconazole Cmax ratio of
    # 1.71 can only be attained if bioavailability can rise by at least
    # about 1.49-fold, which requires F <= 0.8 / 1.49 = 0.536. The
    # independently calibrated F of 0.531 satisfies that bound and no
    # larger value would.
    eh <- fixed(0.219944)
    label("Baseline hepatic extraction ratio")
    egut <- fixed(0.148993)
    label("Baseline gut-wall extraction ratio")

    # -- CYP3A4 modulator effects --------------------------------------
    # Each coefficient is log(relative CYP3A4 activity) for that
    # modulator, back-solved from the modulator's published AUC ratio
    # ALONE. The published Cmax ratio was deliberately held out and is
    # reproduced as an out-of-sample check (vignette Table 4):
    #   itraconazole  AUCR 4.62 (Table S6)   Cmax 1.71 -> predicted 1.56
    #   erythromycin  AUCR 2.68 (Figure 3)   Cmax 1.51 -> predicted 1.39
    #   fluconazole   AUCR 2.48 (Figure 3)   Cmax 1.50 -> predicted 1.36
    #   cimetidine    AUCR 1.42 (Figure 3)   Cmax 1.28 -> predicted 1.14
    #   efavirenz     AUCR 0.349 (Figure 3)  Cmax 0.565 -> predicted 0.61
    #   rifampicin    AUCR 0.167 (Figure 3)  Cmax 0.328 -> predicted 0.40
    # The activities ladder monotonically with the modulators' FDA
    # potency classification, which is itself a consistency check.
    e_conmed_itraconazole_cyp3a4 <- fixed(log(0.1378178))
    label("log relative CYP3A4 activity with itraconazole")
    e_conmed_erythromycin_cyp3a4 <- fixed(log(0.3722123))
    label("log relative CYP3A4 activity with erythromycin")
    e_conmed_fluconazole_cyp3a4 <- fixed(log(0.4101628))
    label("log relative CYP3A4 activity with fluconazole")
    e_conmed_cimetidine_cyp3a4 <- fixed(log(0.7386476))
    label("log relative CYP3A4 activity with cimetidine")
    e_conmed_efv_cyp3a4 <- fixed(log(2.0790900))
    label("log relative CYP3A4 activity with efavirenz")
    e_conmed_rifampicin_cyp3a4 <- fixed(log(3.1347450))
    label("log relative CYP3A4 activity with rifampicin")

    propSd <- fixed(0)
    label("Proportional residual error (none reported by the source)")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Typical-value parameters. No random effects.
    # ------------------------------------------------------------------
    ka    <- exp(lka)
    vc    <- exp(lvc)
    k12   <- exp(lk12)
    k21   <- exp(lk21)
    clint <- exp(lclint)

    # ------------------------------------------------------------------
    # 2. Relative CYP3A4 activity contributed by a coadministered
    # modulator. All six indicators zero (dordaviprone alone) gives
    # a3a4 = 1; the terms are additive on the log scale so the flags are
    # multiplicative on activity.
    # ------------------------------------------------------------------
    a3a4 <- exp(
      e_conmed_itraconazole_cyp3a4 * CONMED_ITRACONAZOLE +
        e_conmed_erythromycin_cyp3a4 * CONMED_ERYTHROMYCIN +
        e_conmed_fluconazole_cyp3a4 * CONMED_FLUCONAZOLE +
        e_conmed_cimetidine_cyp3a4 * CONMED_CIMETIDINE +
        e_conmed_efv_cyp3a4 * CONMED_EFV +
        e_conmed_rifampicin_cyp3a4 * CONMED_RIFAMPICIN
    )

    # Relative hepatic intrinsic clearance. Only the CYP3A4 share
    # responds; the CYP2D6 share (1 - fm_cyp3a4) is untouched. The gut
    # wall is treated as entirely CYP3A4, the standard assumption for a
    # CYP3A4 substrate, so gut intrinsic clearance scales with a3a4.
    ah <- 1 - fm_cyp3a4 + fm_cyp3a4 * a3a4

    # Well-stirred scaling of each extraction ratio. Multiplying an
    # intrinsic clearance by a factor A takes an extraction ratio E to
    # E*A / (1 - E + E*A), so the corresponding availability 1 - E is
    # divided by (1 - E + E*A). These two denominators are therefore all
    # that is needed.
    gfac <- 1 - egut + egut * a3a4
    hfac <- 1 - eh + eh * ah

    # ------------------------------------------------------------------
    # 3. Concentration-dependent unbound fraction and apparent clearance.
    # `cconc` is the total plasma concentration in mg/L, which is the
    # scale the Table 1 binding data are tabulated on.
    # ------------------------------------------------------------------
    cconc <- central / vc
    fu    <- fumin + (fumax - fumin) * cconc / (cup50 + cconc)

    # CL/F is exactly proportional to fu (the hepatic availability
    # cancels between numerator and denominator), and a modulator scales
    # it by the hepatic activity times the gut-availability change.
    cl  <- clint * fu * ah * gfac
    kel <- cl / vc

    # ------------------------------------------------------------------
    # 4. ODE system. `peripheral1` is the paper's single adjusting
    # compartment; k12 and k21 act on masses (mg), matching the paper's
    # k_in / k_out definitions.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ------------------------------------------------------------------
    # 5. Bioavailability. Because vc and clint are apparent parameters
    # (already divided by the baseline bioavailability F), this term
    # carries only the CHANGE in bioavailability caused by a modulator,
    # and equals 1 when dordaviprone is given alone.
    # ------------------------------------------------------------------
    fdepot   <- 1 / (gfac * hfac)
    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # 6. Observation. Doses are in mg and vc is in L, so central / vc is
    # in mg/L = ug/mL; multiply by 1000 to report ng/mL, the units used
    # throughout Jaiswal 2025 Table 2 and Supporting Information.
    # ------------------------------------------------------------------
    Cc <- 1000 * cconc
    Cc ~ prop(propSd)
  })
}
