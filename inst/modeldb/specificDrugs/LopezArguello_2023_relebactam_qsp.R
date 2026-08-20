LopezArguello_2023_relebactam_qsp <- function() {
  description <- "QSP. In vitro (Pseudomonas aeruginosa PAO1). Whole-cell penicillin-binding protein (PBP) covalent-binding model for relebactam (diazabicyclooctane beta-lactamase inhibitor): seven coupled ODEs for the rate of net influx of drug across the outer membrane into periplasm and the competitive, mass-balanced acylation of six PBPs (1a, 1b, 2, 3, 4, 5/6) counted as molecules per bacterial cell. The intact parameter switches between the intact whole-cell assay (penetration-limited; drug enters periplasm at Rate_Influx/access) and the lysed isolated-membrane assay (no outer membrane; a vast excess of drug molecules is present at time 0)."
  reference <- "Lopez-Arguello S, Montaner M, Sayed ARM, Oliver A, Bulitta JB, Moya B. Penicillin-Binding Protein 5/6 Acting as a Decoy Target in Pseudomonas aeruginosa Identified by Whole-Cell Receptor Binding and Quantitative Systems Pharmacology. Antimicrob Agents Chemother. 2023 Jun;67(6):e01603-22. doi:10.1128/aac.01603-22. Structural model: main-text Materials and Methods Eqs 1-4, with the SADAPT-TRAN estimation code in supplemental Fig. S8 taken as the authoritative form. Parameter estimates: Table 1 (rate of net influx and PBP access), Table S1 (second-order acylation rate constants), Table S2 (background-noise, initial-condition-multiplier and residual-error parameters)."
  vignette <- "LopezArguello_2023_pbp_binding_pseudomonas"
  units <- list(time = "min", dosing = "molecules/cell", concentration = "molecules/cell")

  # Mechanistic per-cell molecule counts: six penicillin-binding proteins plus the
  # periplasmic drug pool. None map onto a canonical PK compartment name.
  paper_specific_compartments <- c(
    "npbp1a", "npbp1b", "npbp2", "npbp3", "npbp4", "npbp56", "nperi"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    npbp1a = list(analyte = "PBP 1a", units = NA_character_, specimen = "administration site", verified = FALSE),
    npbp1b = list(analyte = "PBP 1b", units = NA_character_, specimen = "administration site", verified = FALSE),
    npbp2  = list(analyte = "PBP 2", units = NA_character_, specimen = "administration site", verified = FALSE),
    npbp3  = list(analyte = "PBP 3", units = NA_character_, specimen = "administration site", verified = FALSE),
    npbp4  = list(analyte = "PBP 4", units = NA_character_, specimen = "administration site", verified = FALSE),
    npbp56 = list(analyte = "PBP 5/6", units = NA_character_, specimen = "administration site", verified = FALSE),
    nperi  = list(analyte = "relebactam", units = NA_character_, specimen = "administration site", verified = FALSE)
  )

  covariateData <- list(
    # Cell-preparation flag. Fig. S8 lines 66-73 branch on it: it gates BOTH the
    # outer-membrane influx term and the periplasmic initial condition.
    CELLS_INTACT = list(
      description = "Cell-preparation flag for the covalent PBP-binding assay: 1 = intact whole cells, in which drug must cross the outer membrane to reach the periplasmic PBPs; 0 = lysed cells (isolated PBP-containing membranes), in which the outer-membrane barrier is absent and drug is present in vast excess.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = lysed cells (isolated membranes)",
      notes = "When CELLS_INTACT = 1 the periplasm starts empty and fills at Rate_Influx/access; when CELLS_INTACT = 0 the influx term is zero and the periplasm is initialised to n_peri_lysed molecules per cell. Both assays were fit simultaneously.",
      source_name = "INTACT"
    ),
    CONC_REL_MGL = list(
      description = "Static extracellular relebactam concentration applied in the PBP-binding assay. Scales the rate of net influx and PBP access (Eq 1: Rate_Influx/access = Rate_Influx/access,scaled x C_drug).",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- studied at 4 mg/L (a fixed concentration within the clinically relevant range; MIC not determined for beta-lactamase inhibitors)",
      notes = "Table 1, column 'Studied extracellular drug conc. in PBP binding assays'. Held static for the 60-min assay (no dosing events).",
      source_name = "C_drug"
    )
  )

  population <- list(
    species = "in vitro (Pseudomonas aeruginosa PAO1 reference strain)",
    n_subjects = 1L,
    n_studies = 1L,
    disease_state = "Wild-type P. aeruginosa PAO1; relebactam MIC not determined (beta-lactamase inhibitor)",
    model_system = "Covalent PBP-binding assay in intact cells and in lysed cells (isolated membranes). Late exponential-phase cultures (7.6 log10 CFU/mL) in cation-adjusted Mueller-Hinton broth at 37 degrees C were sampled at 0, 15, 30 and 60 min; unbound PBPs were then labelled with 25 uM Bocillin FL and quantified by SDS-PAGE band intensity.",
    dose_range = "relebactam at a static extracellular concentration of 4 mg/L (a fixed concentration within the clinically relevant range; MIC not determined for beta-lactamase inhibitors); no dosing events",
    notes = paste(
      "Estimated by population modelling in S-ADAPT v1.57 (importance sampling, pmethod=4);",
      "the intact-cell and lysed-cell data were fit simultaneously.",
      "The 15 drugs were split into five estimation datasets of four drugs each with imipenem as the shared backbone drug; the Noise, Fini and residual-error parameters in this file are those of dataset 2 (IPM, AVI, REL, FOX) (Table S2).",
      "The total number of PBP molecules per cell (1,731) was borrowed from published Escherichia coli data and split across the six PBPs using the P. aeruginosa relative band intensities (Materials and Methods, Mass balance equations; Table S3)."
    )
  )

  ini({
    # --- Rate of net influx and PBP access (Eq 1; intact-cell assay only) ---
    rate_influx_scaled <- 9.85; label("Scaled rate of net influx and PBP access (drug molecules/min per mg/L)") # Table 1, relebactam: 9.85 -- shared estimate with avibactam (Table 1 footnote c: relebactam's rate was difficult to estimate given its limited PBP binding in intact bacteria and was eventually estimated as a shared parameter with avibactam)
    # The studied extracellular concentration enters as the CONC_REL_MGL covariate
    # (Table 1: relebactam studied at 4 mg/L = a fixed concentration within the clinically relevant range; MIC not determined for beta-lactamase inhibitors).

    # --- Second-order acylation rate constants (Table S1; units 1e-3 /min) ---
    # Fig. S8 lines 75-80 divide each tabulated value by 1000 to get 1/min.
    k2_1a <- 2.41; label("Second-order acylation rate constant, PBP1a (1e-3 /min)") # Table S1, relebactam PBP1a: 2.41
    k2_1b <- 1.69; label("Second-order acylation rate constant, PBP1b (1e-3 /min)") # Table S1, relebactam PBP1b: 1.69
    k2_2 <- 15.1;  label("Second-order acylation rate constant, PBP2 (1e-3 /min)") # Table S1, relebactam PBP2: 15.1 -- shared estimate with avibactam (Table S1 footnote b: the apparent acylation rate constant for PBP2 by relebactam was shared with the estimate of avibactam)
    k2_3 <- 3.2;   label("Second-order acylation rate constant, PBP3 (1e-3 /min)") # Table S1, relebactam PBP3: 3.2
    k2_4 <- 5.93;  label("Second-order acylation rate constant, PBP4 (1e-3 /min)") # Table S1, relebactam PBP4: 5.93
    k2_5 <- 4.72;  label("Second-order acylation rate constant, PBP5/6 (1e-3 /min)") # Table S1, relebactam PBP5/6: 4.72

    # --- Assay configuration (Fig. S8 lines 66-73) ---
    n_peri_lysed <- fixed(1e7); label("Periplasmic drug molecules per cell at time 0, lysed-cell assay") # Fig. S8 line 71: N_IN_IC = 10000000 (vast excess vs the 1,731 PBPs)

    # --- Michaelis-Menten constant, common to all six PBPs ---
    km_pbp <- fixed(1000); label("Michaelis-Menten constant for PBP acylation (drug molecules per cell)") # Materials and Methods: "The Km was fixed to 1,000 drug molecules" (only one concentration studied)

    # --- Nominal PBP expression: 1,731 molecules/cell split by relative band intensity ---
    n_pbp1a_0 <- fixed(153);  label("Nominal PBP1a expression (molecules/cell)") # Materials and Methods, Mass balance equations: 153 molecules for PBP1a (8.8% of all PBPs)
    n_pbp1b_0 <- fixed(118);  label("Nominal PBP1b expression (molecules/cell)") # Materials and Methods, Mass balance equations: 118 molecules for PBP1b (6.8% of all PBPs)
    n_pbp2_0 <- fixed(50);    label("Nominal PBP2 expression (molecules/cell)") # Materials and Methods, Mass balance equations: 50 molecules for PBP2 (2.9% of all PBPs)
    n_pbp3_0 <- fixed(79);    label("Nominal PBP3 expression (molecules/cell)") # Materials and Methods, Mass balance equations: 79 molecules for PBP3 (4.6% of all PBPs)
    n_pbp4_0 <- fixed(99);    label("Nominal PBP4 expression (molecules/cell)") # Materials and Methods, Mass balance equations: 99 molecules for PBP4 (5.7% of all PBPs)
    n_pbp56_0 <- fixed(1232); label("Nominal PBP5/6 expression (molecules/cell)") # Materials and Methods, Mass balance equations: 1232 molecules for PBP5/6 (71.2% of all PBPs)

    # --- Non-suppressible background gel-band fraction (Table S2, dataset 2 (IPM, AVI, REL, FOX)) ---
    noise_1a <- 0.0905; label("Background gel-band fraction not suppressible by drug, PBP1a (unitless)") # Table S2 NOISE_1a = 0.0905
    noise_1b <- 0.152;  label("Background gel-band fraction not suppressible by drug, PBP1b (unitless)") # Table S2 NOISE_1b = 0.152
    noise_2 <- 0.287;   label("Background gel-band fraction not suppressible by drug, PBP2 (unitless)") # Table S2 NOISE_2 = 0.287
    noise_3 <- 0.359;   label("Background gel-band fraction not suppressible by drug, PBP3 (unitless)") # Table S2 NOISE_3 = 0.359
    noise_4 <- 0.125;   label("Background gel-band fraction not suppressible by drug, PBP4 (unitless)") # Table S2 NOISE_4 = 0.125
    noise_5 <- 0.042;   label("Background gel-band fraction not suppressible by drug, PBP5/6 (unitless)") # Table S2 NOISE_5 = 0.042

    # --- Random multiplier on the initial number of PBP molecules (Table S2, dataset 2 (IPM, AVI, REL, FOX)) ---
    fini_1a <- 0.977; label("Multiplier on the initial number of PBP1a molecules (unitless)") # Table S2 Fini,1a = 0.977
    fini_1b <- 0.966; label("Multiplier on the initial number of PBP1b molecules (unitless)") # Table S2 Fini,1b = 0.966
    fini_2 <- 0.983;  label("Multiplier on the initial number of PBP2 molecules (unitless)") # Table S2 Fini,2 = 0.983
    fini_3 <- 0.973;  label("Multiplier on the initial number of PBP3 molecules (unitless)") # Table S2 Fini,3 = 0.973
    fini_4 <- 0.951;  label("Multiplier on the initial number of PBP4 molecules (unitless)") # Table S2 Fini,4 = 0.951
    fini_5 <- 0.979;  label("Multiplier on the initial number of PBP5/6 molecules (unitless)") # Table S2 Fini,5 = 0.979

    # --- Residual error: SD = SDin + SDsl * Y (Fig. S8 lines 111-116; Table S2, dataset 2 (IPM, AVI, REL, FOX)) ---
    addSd_PBP1a <- 4.91; label("Additive residual SD, PBP1a (molecules/cell)") # Table S2 SDin1a = 4.91
    addSd_PBP1b <- 9.89; label("Additive residual SD, PBP1b (molecules/cell)") # Table S2 SDin1b = 9.89
    addSd_PBP2 <- 8.52;  label("Additive residual SD, PBP2 (molecules/cell)") # Table S2 SDin2 = 8.52
    addSd_PBP3 <- 6.87;  label("Additive residual SD, PBP3 (molecules/cell)") # Table S2 SDin3 = 6.87
    addSd_PBP4 <- 9.54;  label("Additive residual SD, PBP4 (molecules/cell)") # Table S2 SDin4 = 9.54
    addSd_PBP56 <- 22.7; label("Additive residual SD, PBP5/6 (molecules/cell)") # Table S2 SDin5 = 22.7

    propSd_PBP1a <- 0.16;   label("Proportional residual SD, PBP1a (fraction)") # Table S2 SDsl1a = 0.16
    propSd_PBP1b <- 0.0697; label("Proportional residual SD, PBP1b (fraction)") # Table S2 SDsl1b = 0.0697
    propSd_PBP2 <- 0.0073;  label("Proportional residual SD, PBP2 (fraction)") # Table S2 SDsl2 = 0.0073
    propSd_PBP3 <- 0.135;   label("Proportional residual SD, PBP3 (fraction)") # Table S2 SDsl3 = 0.135
    propSd_PBP4 <- 0.0881;  label("Proportional residual SD, PBP4 (fraction)") # Table S2 SDsl4 = 0.0881
    propSd_PBP56 <- 0.162;  label("Proportional residual SD, PBP5/6 (fraction)") # Table S2 SDsl5 = 0.162
  })

  model({
    # 1. Rate of net influx and PBP access (Eq 1). Fig. S8 lines 66-73: the influx
    #    term is present only in the intact-cell assay; in lysed cells the outer
    #    membrane is absent and the periplasmic pool is instead initialised to a
    #    vast excess of drug molecules.
    rate_influx <- CELLS_INTACT * rate_influx_scaled * CONC_REL_MGL

    # 2. Table S1 tabulates the acylation rate constants in 1e-3 /min; Fig. S8
    #    lines 75-80 rescale them to 1/min.
    kacyl1a <- k2_1a / 1000
    kacyl1b <- k2_1b / 1000
    kacyl2 <- k2_2 / 1000
    kacyl3 <- k2_3 / 1000
    kacyl4 <- k2_4 / 1000
    kacyl5 <- k2_5 / 1000

    # 3. Michaelis-Menten acylation rates (Fig. S8 lines 14-19). Main-text Eq 2
    #    writes MM as Vmax/(Km + Nperi) and carries the extra Nperi factor inside
    #    Eq 3; the Fig. S8 code folds Nperi into MM instead. The two forms give the
    #    same PBP ODEs, but only the code form makes the periplasmic mass balance
    #    (Eq 4 / Fig. S8 line 28) consistent with them, so the code form is used.
    #    With km_pbp = 1,000 molecules the reaction is effectively first-order in
    #    Nperi in intact cells and fully saturated in the 1e7-molecule lysed assay.
    mm1a <- kacyl1a * nperi / (km_pbp + nperi)
    mm1b <- kacyl1b * nperi / (km_pbp + nperi)
    mm2 <- kacyl2 * nperi / (km_pbp + nperi)
    mm3 <- kacyl3 * nperi / (km_pbp + nperi)
    mm4 <- kacyl4 * nperi / (km_pbp + nperi)
    mm5 <- kacyl5 * nperi / (km_pbp + nperi)

    # 4. Unbound PBP molecules per cell (Eq 3; Fig. S8 lines 21-26). Deacylation is
    #    not modelled - no significant dissociation was seen over the 60-min study.
    d/dt(npbp1a) <- -mm1a * npbp1a
    d/dt(npbp1b) <- -mm1b * npbp1b
    d/dt(npbp2) <- -mm2 * npbp2
    d/dt(npbp3) <- -mm3 * npbp3
    d/dt(npbp4) <- -mm4 * npbp4
    d/dt(npbp56) <- -mm5 * npbp56

    # 5. Periplasmic drug molecules per cell (Eq 4; Fig. S8 line 28). Each acylation
    #    event consumes one drug molecule, so the six PBPs compete for the same pool.
    d/dt(nperi) <- rate_influx - mm1a * npbp1a - mm1b * npbp1b - mm2 * npbp2 -
      mm3 * npbp3 - mm4 * npbp4 - mm5 * npbp56

    # 6. Initial conditions (Fig. S8 lines 84-90). Only the drug-suppressible part of
    #    each band is a state; the non-suppressible background is added back in the
    #    output equation, so the observed signal starts at n_pbp*_0 * fini_*.
    npbp1a(0) <- n_pbp1a_0 * (1 - noise_1a) * fini_1a
    npbp1b(0) <- n_pbp1b_0 * (1 - noise_1b) * fini_1b
    npbp2(0) <- n_pbp2_0 * (1 - noise_2) * fini_2
    npbp3(0) <- n_pbp3_0 * (1 - noise_3) * fini_3
    npbp4(0) <- n_pbp4_0 * (1 - noise_4) * fini_4
    npbp56(0) <- n_pbp56_0 * (1 - noise_5) * fini_5
    nperi(0) <- (1 - CELLS_INTACT) * n_peri_lysed

    # 7. Observed number of unbound PBP molecules per cell (Fig. S8 lines 102-107).
    #    max(., 0) reproduces the negative-state guard on Fig. S8 lines 94-100.
    PBP1a <- max(npbp1a, 0) + n_pbp1a_0 * noise_1a * fini_1a
    PBP1b <- max(npbp1b, 0) + n_pbp1b_0 * noise_1b * fini_1b
    PBP2 <- max(npbp2, 0) + n_pbp2_0 * noise_2 * fini_2
    PBP3 <- max(npbp3, 0) + n_pbp3_0 * noise_3 * fini_3
    PBP4 <- max(npbp4, 0) + n_pbp4_0 * noise_4 * fini_4
    PBP56 <- max(npbp56, 0) + n_pbp56_0 * noise_5 * fini_5

    # 8. Residual error: Fig. S8 lines 111-116 give V = (SDin + SDsl * Y)^2, i.e. a
    #    linear sum of the additive and proportional SDs -> combined1() in nlmixr2.
    PBP1a ~ add(addSd_PBP1a) + prop(propSd_PBP1a) + combined1()
    PBP1b ~ add(addSd_PBP1b) + prop(propSd_PBP1b) + combined1()
    PBP2 ~ add(addSd_PBP2) + prop(propSd_PBP2) + combined1()
    PBP3 ~ add(addSd_PBP3) + prop(propSd_PBP3) + combined1()
    PBP4 ~ add(addSd_PBP4) + prop(propSd_PBP4) + combined1()
    PBP56 ~ add(addSd_PBP56) + prop(propSd_PBP56) + combined1()
  })
}
