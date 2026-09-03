Saporta_2026_meropenem <- function() {
  description <- "Preclinical (mouse, CD-1 female). Coupled pharmacokinetic-pharmacodynamic model of meropenem against Klebsiella pneumoniae DSM116099 in a standardized murine lung infection model (COMBINE pneumonia model) run at three experimentally-induced immune states: neutropenic, intermediate suppression, and immunocompetent (cyclophosphamide 200+150, 75+50, or 0 mg/kg intraperitoneally at 4 and 1 days before infection). Subcutaneous meropenem PK is a one-compartment plasma model with first-order absorption; lung epithelial lining fluid (ELF) is described by a two-compartment limb (elf plus a second lung compartment) that is DRIVEN by the plasma concentration but does not deplete it, because mass transfer between plasma and lung was deliberately not retained in the final model (Figure 2 draws that link dashed). The only PK parameter that differs by immune state is the apparent central volume (2.19, 2.46, 3.40 L/kg), which is what makes both plasma and ELF profiles differ between immune states. The bacterial system follows the Nielsen 2007 semi-mechanistic lineage extended with an immune-response limb: a growing drug-susceptible state (bact_susceptible), a dormant non-growing drug-insusceptible state (bact_resting) entered at kSD = (S + D) * (kgrowth - kdeath) / Bmax, and a phagocytosed state (bact_phagocytosed) that bacteria enter from both other states at the phagocytosis rate kphag and leave by digestion at kdig (constrained to equal kphag). Phagocytosed bacteria still count toward the observed CFU but not toward Bmax, because phagocytic cells were lysed before plating. kphag is built up additively across immune states as IRneu + IRint + IRcom (0, 0.185 and 0.318 1/h in the neutropenic, intermediate and immunocompetent states), with IRneu fixed to 0 because in neutropenic mice the immune contribution was not differentiable from bacterial growth. Meropenem adds a killing rate on susceptible bacteria only, kdrug = Emax * Cu / (EC50 + Cu), driven by the UNBOUND PLASMA concentration (fu = 0.81) rather than by ELF, which fitted better by 15.7 OFV points. The reduced meropenem contribution in immunocompetent mice arises structurally, not from any change in Emax or EC50: a larger fraction of the bacteria sits in the phagocytosed state where meropenem has no effect. Model time zero is 2 h after infection, the start of treatment, at which each immune state's susceptible state is initialized to the observed median control count."
  reference <- "Saporta R, Tassi N, Biordi V, Ticha O, Ginosyan A, Loryan I, Nielsen EI, Bekeredjian-Ding I, Kerscher B, Friberg LE. Pharmacokinetic-pharmacodynamic modeling to evaluate the relative impact of immune response and meropenem on bacterial killing in vivo. Antimicrob Agents Chemother. 2026 Apr;70(4):e01788-25. doi:10.1128/aac.01788-25. PMCID: PMC13041408. All parameter estimates: Table 1. Bacterial-system equations 1-3 and the kSD definition: Materials and Methods, 'PKPD modeling'. Emax drug-effect equation 5 with gamma = 1: Materials and Methods, 'PKPD modeling', and Results, 'Mouse PD experiments and PKPD modeling' ('The Emax model was selected'). Unbound fraction fu = 0.81: Materials and Methods, 'PKPD modeling' (cited to reference 19). Model structure including the non-depleting plasma-to-ELF link: Figure 2 and Results, 'Mouse PK experiments and PK modeling'. Immune-state phagocytosis rates 0 / 0.185 / 0.318 1/h and the intermediate-state volume interpolation: Results. Median control bacterial counts at 2 h after infection (7.56, 7.80, 7.52 log10 CFU/lung): Results, 'Mouse PD experiments and PKPD modeling'. Klebsiella pneumoniae DSM116099 MIC 0.032 mg/L: Materials and Methods, 'Mouse lung infection model'."
  vignette <- "Saporta_2026_meropenem"
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")

  # Bacterial life-cycle states. The registered
  # ^bact_(susceptible|intermediate|resistant) regex is a DRUG-RESISTANCE
  # phenotype namespace; this model splits one phenotype along a different
  # axis (growing / dormant / phagocytosed), so the states are declared
  # paper-specific exactly as Nielsen_2007_semimechanistic_antibiotic_pd.R
  # does for its growing / resting pair.
  paper_specific_compartments <- c(
    "bact_susceptible", "bact_resting", "bact_phagocytosed"
  )

  compartmentData <- list(
    depot = list(analyte = "meropenem", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "meropenem", units = "mg/kg", specimen = "plasma", verified = TRUE),
    elf = list(analyte = "meropenem", units = "mg/kg", specimen = "epithelial lining fluid", verified = TRUE),
    lung = list(analyte = "meropenem", units = "mg/kg", specimen = "tissue", verified = TRUE),
    bact_susceptible = list(analyte = "Klebsiella pneumoniae DSM116099, growing drug-susceptible (S) state", units = "CFU/lung", specimen = "not applicable", verified = TRUE),
    bact_resting = list(analyte = "Klebsiella pneumoniae DSM116099, dormant non-growing drug-insusceptible (D) state", units = "CFU/lung", specimen = "not applicable", verified = TRUE),
    bact_phagocytosed = list(analyte = "Klebsiella pneumoniae DSM116099, phagocytosed (P) state inside host immune cells", units = "CFU/lung", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    IMMUNE_STATE = list(
      description = "Experimentally-induced immune status of the mouse: 0 = neutropenic, 1 = intermediate suppression, 2 = immunocompetent. Selects the apparent central volume of distribution, the phagocytosis rate constant, and the initial susceptible bacterial count.",
      units = "(categorical)",
      type = "categorical",
      source_name = "immune status",
      reference_category = "0 (neutropenic)",
      notes = paste(
        "Induced by intraperitoneal cyclophosphamide at 4 and 1 days before infection: 200 + 150 mg/kg for the neutropenic state, 75 + 50 mg/kg for the intermediate state, and none for the immunocompetent state (Materials and Methods, 'Mouse lung infection model').",
        "Ordered: the three levels are monotone in immune competence, and the phagocytosis rate is built up additively across them (kphag = IRneu, IRneu + IRint, IRneu + IRint + IRcom).",
        "The paper used the immune status as a categorical covariate rather than scaling by the measured blood granulocyte counts (0.2, 0.4 and 1.9 x 10^3 /mm3 in the three states); the granulocyte-scaled parameterization 'led to a worse model fit' (Results).",
        "Neutropenic is the reference level: it is the standard preclinical infection-model condition, it carries the base central-volume estimate, and its phagocytosis increment IRneu is fixed to 0 because the immune contribution was not differentiable from the bacterial growth rate."
      )
    )
  )

  population <- list(
    species = "mouse (CD-1, female, specific pathogen-free, 6-10 weeks old at intervention start)",
    n_subjects = 240L,
    n_studies = 1L,
    sex_female_pct = 100,
    disease_state = "Klebsiella pneumoniae DSM116099 lung infection (meropenem MIC 0.032 mg/L by EUCAST broth microdilution in triplicate), established by intranasal instillation of 50 uL of bacterial suspension into the nares under isoflurane anesthesia, per the standardized COMBINE murine pneumonia model with the immune status modified",
    dose_range = "Meropenem 40 or 300 mg/kg subcutaneously; a single dose at 2 h post-infection in the PK experiments, and every 4 h starting at 2 h post-infection in the PD experiments",
    regions = "Paul-Ehrlich-Institut, Langen, Germany (in vivo experiments); Uppsala University, Sweden (modeling)",
    immune_states = "Neutropenic (cyclophosphamide 200 + 150 mg/kg intraperitoneally at 4 and 1 days pre-infection), intermediate suppression (75 + 50 mg/kg), and immunocompetent (no cyclophosphamide)",
    notes = paste(
      "60 mice contributed the PK experiments (plasma at 5, 15, 30, 60, 120, 180 and 240 min after administration, 1-2 samples per animal; bronchoalveolar lavage at 5, 15, 30, 180 and 240 min, 1 sample per animal) and 180 mice the PD experiments (blood and whole lung tissue at 2, 4, 8, 12, 20 and 26 h after infection, 1 sample per animal). Table S1.",
      "Meropenem ELF concentrations were not measured directly: they were derived from bronchoalveolar-lavage concentrations by urea correction, C_ELF = C_BAL * Urea_plasma / Urea_BAL. LLOQ 0.05 mg/L in plasma and BAL; data below LLOQ handled by the M3 method.",
      "MODEL TIME ZERO IS 2 h AFTER INFECTION, which is both the start of treatment and the time at which the initial susceptible bacterial counts apply. The paper's '26 h' endpoint (24 h after start of treatment) is therefore model time 24 h, and the PD sampling times 2, 4, 8, 12, 20 and 26 h after infection are model times 0, 2, 6, 10, 18 and 24 h.",
      "No between-subject variability was estimated: the design is destructively sampled (1-2 samples per mouse), so the model carries typical values plus residual error only.",
      "Simulated dose-ranging (20 to 1200 mg/kg every 4 h) and dose-fractionation (equal total doses split every 2, 4, 6 or 8 h) studies used an inoculum of 7.5 log10 CFU/lung and were run without residual variability."
    )
  )

  ini({
    # ---- Plasma PK (Table 1) ---------------------------------------------
    lcl <- log(5.33); label("Log apparent meropenem clearance (L/h/kg)")                                       # Table 1: CL = 5.33 L/(h.kg), RSE 11%
    lka <- log(71.3); label("Log first-order subcutaneous absorption rate constant (1/h)")                     # Table 1: ka = 71.3 1/h, RSE 39%

    # Apparent central volume is the ONLY PK parameter that differs by immune
    # state (Results: "best described by different central volumes of
    # distribution"). The intermediate-state value was not estimated: it was
    # "scaled to 2.46 L/kg by linear interpolation based on median white blood
    # cell counts from the PD experiments" (Results), hence fixed().
    lvcNeu <- log(2.19); label("Log apparent central volume of distribution, neutropenic mice (L/kg)")              # Table 1: Vneu = 2.19 L/kg, RSE 18%
    lvcInt <- fixed(log(2.46)); label("Log apparent central volume of distribution, intermediate suppression (L/kg)") # Table 1: Vint = 2.46 L/kg, Fixed
    lvcCom <- log(3.40); label("Log apparent central volume of distribution, immunocompetent mice (L/kg)")          # Table 1: Vcom = 3.40 L/kg, RSE 18%

    # ---- Lung / ELF limb (Table 1) ---------------------------------------
    lv_elf <- log(0.0554); label("Log apparent epithelial-lining-fluid volume of distribution (L/kg)")         # Table 1: VELF = 0.0554 L/kg, RSE 23%
    lv_lung <- log(9.27); label("Log apparent volume of distribution of the second lung compartment (L/kg)")   # Table 1: VL2 = 9.27 L/kg, RSE 26%
    lq_elf <- log(0.385); label("Log apparent central-to-ELF inter-compartmental clearance (L/h/kg)")          # Table 1: Q = 0.385 L/(h.kg), RSE 19%
    lq_elf_lung <- log(2.48); label("Log apparent ELF-to-lung inter-compartmental clearance (L/h/kg)")         # Table 1: Q2 = 2.48 L/(h.kg), RSE 13%

    # ---- Bacterial system (Table 1) --------------------------------------
    lkgrowth <- log(0.372); label("Log bacterial growth rate constant (1/h)")                                  # Table 1: kgrowth = 0.372 1/h, RSE 11%
    lkdeath <- fixed(log(0.179)); label("Log bacterial natural death rate constant, all states (1/h)")         # Table 1: kdeath = 0.179 1/h, Fixed
    # Bmax is tabulated on the log10 CFU/lung scale but enters kSD on the
    # LINEAR CFU scale, so it is carried here as log(10^9.59).
    lbmax <- fixed(log(10^9.59)); label("Log maximum bacterial count of the system (CFU/lung)")                # Table 1: Bmax = 9.59 log10 CFU/lung, Fixed; fixed "to the median bacterial count in neutropenic controls at 26 hours after infection"

    # ---- Immune response: additive phagocytosis-rate increments ----------
    # Table 1 parameterizes the phagocytosis rate cumulatively across immune
    # states: kphag = IRneu (neutropenic), IRneu + IRint (intermediate),
    # IRneu + IRint + IRcom (immunocompetent). Kept on the natural scale
    # because they are additive increments and IRneu is exactly 0.
    irNeu <- fixed(0); label("Phagocytosis rate constant in the neutropenic state (1/h)")                      # Table 1: IRneu = 0 1/h, Fixed; "not differentiable from kgrowth, and was therefore fixed to 0"
    irInt <- 0.185; label("Increase in phagocytosis rate constant in the intermediate suppression state (1/h)") # Table 1: IRint = 0.185 1/h, RSE 12%
    irCom <- 0.133; label("Increase in phagocytosis rate constant in the immunocompetent state (1/h)")          # Table 1: IRcom = 0.133 1/h, RSE 18%

    # ---- Meropenem effect (Table 1) --------------------------------------
    lemax <- log(0.934); label("Log maximum meropenem killing rate constant (1/h)")                            # Table 1: Emax = 0.934 1/h, RSE 14%
    lec50 <- log(1.62); label("Log unbound meropenem concentration giving half of Emax (mg/L)")                # Table 1: EC50 = 1.62 mg/L, RSE 18%
    fu <- fixed(0.81); label("Fraction of meropenem unbound in mouse plasma (unitless)")                       # Materials and Methods, PKPD modeling: "Unbound meropenem plasma (fu = 0.81) (19)"

    # ---- Initial bacterial counts at 2 h after infection -----------------
    # "For each immune state, the S compartment was initialized to the median
    # bacterial count from control groups at 2 h after infection" (Materials
    # and Methods); the three medians are reported in Results. Observed data
    # summaries, not estimated parameters, hence fixed().
    linocNeu <- fixed(log(10^7.56)); label("Log initial susceptible bacterial count, neutropenic mice (CFU/lung)")        # Results: median observed 7.56 log10 CFU/lung at 2 h after infection
    linocInt <- fixed(log(10^7.80)); label("Log initial susceptible bacterial count, intermediate suppression (CFU/lung)") # Results: median observed 7.80 log10 CFU/lung at 2 h after infection
    linocCom <- fixed(log(10^7.52)); label("Log initial susceptible bacterial count, immunocompetent mice (CFU/lung)")     # Results: median observed 7.52 log10 CFU/lung at 2 h after infection

    # ---- Residual error --------------------------------------------------
    # PK was fitted log-transform-both-sides with additive terms on the
    # log-transformed scale, which is an exponential (log-normal) residual on
    # the linear concentration scale. PD used a log10-transform-both-sides
    # approach, so its residual is additive on the log10 CFU scale.
    #
    # ERRATUM: Table 1's two PK residual rows are internally transposed --
    # the row NAMED RESPlasma is DESCRIBED as "Residual error of ELF PK" and
    # the row named RESELF as "Residual error of plasma PK". The parameter-name
    # column is taken as authoritative here (operator ruling); see the vignette
    # Errata. This assignment is an inference, not a printed fact.
    expSd <- 0.816; label("Exponential residual SD on the plasma meropenem concentration (log scale)")  # Table 1: RESPlasma = 0.816, RSE 17%
    expSd_Celf <- 0.738; label("Exponential residual SD on the ELF meropenem concentration (log scale)") # Table 1: RESELF = 0.738, RSE 22%
    addSd_cfu <- 0.623; label("Additive residual SD on the log10 bacterial count (log10 CFU/lung)")     # Table 1: RESPD = 0.623 log10 CFU/lung, RSE 14%
  })

  model({
    # ---- 1. Immune-state-dependent quantities -----------------------------
    # IMMUNE_STATE: 0 = neutropenic, 1 = intermediate suppression,
    # 2 = immunocompetent. Half-integer cut points so a numeric column that
    # carries the levels as doubles still dispatches correctly.
    if (IMMUNE_STATE < 0.5) {
      vc <- exp(lvcNeu)
      kphag <- irNeu
      inoc <- exp(linocNeu)
    } else if (IMMUNE_STATE < 1.5) {
      vc <- exp(lvcInt)
      kphag <- irNeu + irInt
      inoc <- exp(linocInt)
    } else {
      vc <- exp(lvcCom)
      kphag <- irNeu + irInt + irCom
      inoc <- exp(linocCom)
    }
    # "The digestion rate kdig was not significantly different from kphag"
    # (Results), so the final model constrains the two to be equal. No
    # independent estimate exists in Table 1.
    kdig <- kphag

    # ---- 2. Back-transform the log-parameterized values -------------------
    cl <- exp(lcl)
    ka <- exp(lka)
    v_elf <- exp(lv_elf)
    v_lung <- exp(lv_lung)
    q_elf <- exp(lq_elf)
    q_elf_lung <- exp(lq_elf_lung)
    kgrowth <- exp(lkgrowth)
    kdeath <- exp(lkdeath)
    bmax <- exp(lbmax)
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    kel <- cl / vc

    # ---- 3. Plasma PK: one compartment, first-order SC absorption ---------
    # The Q flux into the lung limb does NOT appear here: "Mass transfer
    # between plasma and lung compartments was not retained in the model to
    # improve stability and provide reliable plasma predictions" (Results),
    # which Figure 2 draws as a dashed C <-> L_ELF link. The lung limb is
    # therefore driven by the plasma concentration without depleting it.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    Cc <- central / vc

    # ---- 4. Lung limb: ELF plus a second lung compartment -----------------
    Celf <- elf / v_elf
    Clung <- lung / v_lung
    d/dt(elf) <- q_elf * (Cc - Celf) - q_elf_lung * (Celf - Clung)
    d/dt(lung) <- q_elf_lung * (Celf - Clung)

    # ---- 5. Meropenem effect (equation 5, gamma = 1) ----------------------
    # Driven by the UNBOUND PLASMA concentration: "The model fit was better
    # (dOFV = -15.7) when the plasma concentration-time profile, rather than
    # the ELF concentration-time profile, was driving the bacterial killing"
    # (Results). Power and Emax models performed similarly (dOFV = -1.8) and
    # the Emax model was selected for biological plausibility.
    cu <- fu * Cc
    kdrug <- emax * cu / (ec50 + cu)

    # ---- 6. Bacterial system (equations 1-3) ------------------------------
    # kSD is linear in the S + D burden with proportionality
    # (kgrowth - kdeath) / Bmax, so an untreated non-phagocytosing system
    # plateaus at Bmax. Phagocytosed bacteria are deliberately excluded from
    # this term: they "contribute to the total CFU counts but not to Bmax, as
    # phagocytic cells were lysed prior to the counting procedure".
    bactSD <- bact_susceptible + bact_resting
    kSD <- (kgrowth - kdeath) * bactSD / bmax

    d/dt(bact_susceptible) <- kgrowth * bact_susceptible -
      (kdeath + kdrug) * bact_susceptible -
      kSD * bact_susceptible -
      kphag * bact_susceptible
    d/dt(bact_resting) <- kSD * bact_susceptible -
      kdeath * bact_resting -
      kphag * bact_resting
    d/dt(bact_phagocytosed) <- kphag * bact_susceptible +
      kphag * bact_resting -
      kdig * bact_phagocytosed

    # ---- 7. Initial conditions --------------------------------------------
    # Model time zero is 2 h after infection (the start of treatment). All
    # bacteria start in the growing susceptible state.
    bact_susceptible(0) <- inoc

    # ---- 8. Observations ---------------------------------------------------
    # Total observed CFU includes the phagocytosed state. The 1e-6 floor keeps
    # log10() finite if the integrator drives the total to numerical zero.
    cfu <- log10(bact_susceptible + bact_resting + bact_phagocytosed + 1e-6)

    Cc ~ lnorm(expSd)
    Celf ~ lnorm(expSd_Celf)
    cfu ~ add(addSd_cfu)
  })
}
