Hyland_2009_midazolam_hlm <- function() {
  description <- "In vitro (pooled human liver microsomes, HLM101). Two-site substrate-inhibitory cooperative binding model for the direct N-glucuronidation of midazolam to midazolam N-glucuronide (MDZG) by human liver microsomes. Midazolam N-glucuronidation does not follow Michaelis-Menten kinetics: the Eadie-Hofstee plot shows both allosteric activation at low substrate and substrate inhibition at high substrate, so the authors fitted a combination of the Hill equation with the two-site substrate-inhibitory equation of Shou et al. The substrate sits in a static incubation compartment (initial-rate assay, 20 min, 0.5 mg microsomal protein per mL, 0-800 uM midazolam) and MDZG accumulates at the instantaneous rate v; the model observable is that formation rate in pmol/min/mg protein. The rate rises to a maximum of about 229 pmol/min/mg near 100 uM midazolam and then falls back toward Vmax * beta at high substrate. Sibling model: Hyland_2009_midazolam_rugt1a4, the same equation independently fitted to recombinant UGT1A4 Supersomes, which the paper identified as the only UGT of the twelve screened that conjugates midazolam directly."
  reference <- "Hyland R, Osborne T, Payne A, Kempshall S, Logan YR, Ezzeddine K, Jones B. In vitro and in vivo glucuronidation of midazolam in humans. Br J Clin Pharmacol. 2009 Apr;67(4):445-454. doi:10.1111/j.1365-2125.2009.03386.x. PMID: 19371318. PMCID: PMC2679108. Model equation: Equation 2 and its symbol-definition block, Results, 'MDZG enzyme kinetics' (p. 449). Parameter estimates: Table 4, column 'MDZ HLM' (p. 450). The fixed value of beta: Results, 'MDZG enzyme kinetics' -- 'since initial fitting of the data indicated a very low beta value (<0.01), this value was set as a constant of 0.01, to allow for better fitting of other variables'. Assay conditions: Materials and methods, 'Microsomal incubations', subsection 'MDZG kinetics'. Fitted curve reproduced here: Figure 3, open circles (HLM101)."
  vignette <- "Hyland_2009_midazolam_glucuronidation"

  # "min" and the incubation's natural concentration units. The observable is a
  # specific formation RATE rather than a concentration; the `concentration`
  # slot names the observable's units so the source trace is unambiguous.
  units <- list(
    time = "min",
    dosing = "umol/L (midazolam concentration placed in the incubation)",
    concentration = "pmol/min/mg protein (MDZG formation rate)"
  )

  # The dosing target is neither `depot` nor `central`: this is a static in
  # vitro incubation, so the "dose" is the midazolam concentration placed
  # directly into the `midazolam` state at time 0. Same idiom as
  # HernandezLozano_2025_apramycin_invitro.R.
  dosing <- c("midazolam")

  # `midazolam` holds the static incubation substrate concentration in uM,
  # following the in vitro convention of HernandezLozano_2025_apramycin_invitro.R
  # (`apramycin` for the static bath concentration). `mdzg` accumulates the
  # glucuronide formed per mg of microsomal protein. Neither is a canonical
  # role in compartment-names.md: this is an enzyme-kinetic incubation, not a
  # body, so there is no depot / central / peripheral structure to map onto.
  paper_specific_compartments <- c("midazolam", "mdzg")

  compartmentData <- list(
    midazolam = list(analyte = "midazolam", units = "umol/L", specimen = "administration site", verified = TRUE),
    mdzg = list(analyte = "midazolam N-glucuronide", units = "pmol/mg protein", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "in vitro (pooled human liver microsomes)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    system = "Pooled human liver microsomes (HLM101, BD Biosciences) at 0.5 mg protein/mL",
    matrix = "50 mM Tris-HCl pH 7.4 with 5 mM saccharolactone, alamethacin at 50 ug/mg protein and 10 mM MgCl2; the mixture was held on ice for 15 min so alamethacin could form pores in the microsomal membrane before midazolam was added",
    initiation = "Warmed to 37 C and initiated with 5 mM UDP-glucuronic acid (UDPGA)",
    temperature = "37 C",
    incubation_time = "20 min (initial-rate conditions), terminated with 3 volumes of ice-cold acetonitrile",
    substrate_range = "0-800 umol/L midazolam, n = 3 per concentration; 18 substrate concentrations were used because the equation has five unknown variables",
    replicates = "triplicate",
    quantification = "HPLC-MS/MS (Sciex API 3000) against a calibration curve prepared from an authentic MDZG standard isolated from HLM incubations and confirmed by 1H NMR and 1H-13C gHMBC to be conjugated on the alpha-nitrogen of the imidazole ring (90.3% purity)",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "Curve fitting was done in Grafit v5; the two-site substrate-inhibitory cooperative binding model (Equation 2) was compared against equation 7 of Houston and Kenworthy by F-test and gave the better fit (P < 0.001 and lower chi-squared). R-squared for this HLM fit was 0.962 (Table 4).",
      "The same paper's in vivo arm -- six healthy men aged 22-54 years weighing 71-94 kg who received a single 3 mg oral or 1 mg intravenous midazolam dose -- quantified MDZG in urine over 24 h, roughly 1-2% of the administered dose on a molar basis. That arm is descriptive mass-balance accounting, not a fitted pharmacokinetic model, so it is not encoded here; it is reported in the vignette. Note that the intravenous summary row of Table 3 is internally inconsistent: it prints a geometric mean of 26.5 ug, but recomputing from the six tabulated subject values gives 30.6 ug, which is what the Results text (30.7 +/- 5.7 ug) and the Discussion (about 31 ug) both report. The printed 26.5 is subject 10's own intravenous value. The oral column is self-consistent at 48.8 ug. The vignette shows the recomputation.",
      "UGT1A4 was the only enzyme of the twelve recombinant UGTs screened (1A1, 1A3, 1A4, 1A6, 1A7, 1A8, 1A9, 1A10, 2B4, 2B7, 2B15, 2B17, alongside a control) that conjugates midazolam directly, favoured over UGT2B4 and UGT2B7 by more than 100-fold and 400-fold respectively.",
      "Ketoconazole and itraconazole inhibited MDZG formation only weakly, with IC50 values of 150 +/- 6 uM and 308 +/- 14 uM (Figure 4); trifluoperazine activated the pathway to a maximum of 138% at 60 uM before inhibiting it. None of those three effects is encoded in this model -- see the vignette Errata for why."
    )
  )

  ini({
    # ---- Table 4, column 'MDZ HLM' --------------------------------------
    # Vmax, Ks and Ki are positive-constrained kinetic constants and are
    # carried on the natural-log scale per the library convention; alpha and
    # n are dimensionless shape factors and are carried on the natural scale
    # (same treatment as `gam` and `redukg` in
    # HernandezLozano_2025_apramycin_invitro.R). Every value is an estimate
    # from the Grafit fit, reported +/- its standard error, except beta.
    lvmax <- log(445); label("Log maximum catalytic rate at full saturation (pmol/min/mg protein)")  # Table 4, MDZ HLM: Vmax = 445 +/- 1 pmol min-1 mg-1
    lks <- log(46); label("Log dissociation constant for productive substrate binding to the enzyme (umol/L)")  # Table 4, MDZ HLM: Ks = 46 +/- 5 uM. Called "Km (46 uM)" in the Abstract
    lki <- log(58); label("Log dissociation constant for inhibitory substrate binding to the enzyme (umol/L)")  # Table 4, MDZ HLM: Ki = 58 +/- 6 uM
    alpha <- 18.3; label("Ks modifying factor for binding at the second site (unitless)")  # Table 4, MDZ HLM: a = 18.3 +/- 0.0007
    nhill <- 2.3; label("Hill coefficient, the slope of the allosteric effect (unitless)")  # Table 4, MDZ HLM: n = 2.3 +/- 0.05

    # beta is the turnover-rate modifying factor: the catalytic rate of the
    # doubly-occupied enzyme relative to the singly-occupied enzyme, so
    # Vmax * beta is the residual rate the curve decays to at saturating
    # substrate. It is the one parameter the authors did NOT estimate.
    beta <- fixed(0.01); label("Turnover rate modifying factor for the doubly-occupied enzyme (unitless)")  # Results, MDZG enzyme kinetics: "initial fitting of the data indicated a very low b value (<0.01), this value was set as a constant of 0.01"

    # ---- Residual error ---------------------------------------------------
    # The paper reports the fit quality as R-squared = 0.962 (Table 4) and
    # gives no residual standard deviation, and Grafit reports no $SIGMA
    # analogue. Held at 0 rather than invented; see the vignette Errata.
    addSd <- fixed(0); label("Additive residual SD on the MDZG formation rate (pmol/min/mg protein); magnitude not reported")  # Not reported: Table 4 gives R-squared only
  })

  model({
    # ---- 1. Back-transform ------------------------------------------------
    vmax <- exp(lvmax)
    ks <- exp(lks)
    ki <- exp(lki)

    # ---- 2. Equation 2 ----------------------------------------------------
    # DEVIATION FROM THE EQUATION AS TYPESET. Equation 2 is printed as
    #
    #        Vmax * [ 1/Ks + (beta * [S]^n) / (alpha*Ks*Ki) ]
    #   v = --------------------------------------------------
    #        1/[S]^n + 1/Ks + 1/Ki + [S]^n / (alpha*Ks*Ki)
    #
    # i.e. with Ks and Ki carrying no exponent. That form cannot be what was
    # fitted, for three independent reasons:
    #
    #  (a) Dimensional analysis. The denominator adds 1/[S]^n (units uM^-n,
    #      with n = 2.3) to 1/Ks (units uM^-1). Those are only commensurable
    #      when n = 1. Raising Ks and Ki to the n-th power -- the standard
    #      Hill generalisation of a dissociation constant -- is the unique
    #      assignment that makes every term uM^-n.
    #  (b) It does not reproduce the paper's own Figure 3. With the Table 4
    #      values, the equation as typeset peaks at 202 pmol/min/mg at 10.5 uM
    #      midazolam and has fallen to 16 pmol/min/mg by 100 uM. Figure 3
    #      shows the fitted HLM curve peaking at about 230 pmol/min/mg near
    #      100 uM. The Ks^n / Ki^n form gives 229 pmol/min/mg at 98 uM and
    #      tracks the published curve along both limbs (175 at 200 uM, 75 at
    #      400 uM, 22 at 800 uM against roughly 175, 77 and 20 read off
    #      Figure 3).
    #  (c) It contradicts Figure 4, whose uninhibited controls were 150.0-160.2
    #      pmol/min/mg with midazolam "incubated at a concentration
    #      approximating to the Km". At 46 uM the form as typeset gives 58
    #      pmol/min/mg; the Ks^n / Ki^n form gives 170.
    #
    # The source is "a modification of equation 11 from Shou et al." combined
    # with the Hill equation (Materials and methods, Data analysis), and the
    # Shou two-site substrate-inhibition model with [S] replaced by [S]^n is
    # exactly the Ks^n / Ki^n form. The typeset equation has simply lost the
    # exponents on Ks and Ki. The vignette gates all three checks above.
    #
    # Written below multiplied through by [S]^n, which is algebraically
    # identical and well behaved at [S] = 0 (where the printed form evaluates
    # 1/0). This is also the canonical way Shou et al. write the model.
    ksn <- ks^nhill
    kin <- ki^nhill
    sn <- midazolam^nhill

    v <- vmax * (sn / ksn + beta * sn * sn / (alpha * ksn * kin)) /
      (1 + sn / ksn + sn / kin + sn * sn / (alpha * ksn * kin))

    # ---- 3. The incubation ------------------------------------------------
    # Initial-rate assay: the 20 min incubation consumes only a small fraction
    # of the substrate, so the concentration is held constant, exactly as the
    # published rate law assumes. The fraction turned over is
    # v * 20 min * 0.5 mg protein/mL / [S], which reduces to 0.01 * v/[S];
    # over the studied 10-800 umol/L range that peaks at 3.7% near 41 umol/L
    # and is 1.3% at 10 umol/L and 2.3% at the 100 umol/L rate maximum. The
    # vignette computes this explicitly rather than asserting it.
    d/dt(midazolam) <- 0
    d/dt(mdzg) <- v

    # ---- 4. Observation ----------------------------------------------------
    # The paper's observable is the initial rate v (Figure 3 plots v against
    # midazolam concentration), not a concentration. `Cc` is the library's
    # required single-output observation name and carries that rate here;
    # the accumulated glucuronide is available as the `mdzg` state.
    Cc <- v
    Cc ~ add(addSd)
  })
}
