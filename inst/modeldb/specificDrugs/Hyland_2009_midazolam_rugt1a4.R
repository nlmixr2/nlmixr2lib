Hyland_2009_midazolam_rugt1a4 <- function() {
  description <- "In vitro (recombinant human UGT1A4 Supersomes). Two-site substrate-inhibitory cooperative binding model for the direct N-glucuronidation of midazolam to midazolam N-glucuronide (MDZG) by recombinant UGT1A4. UGT1A4 was the only enzyme of the twelve recombinant UGTs screened that conjugates midazolam directly, favoured over UGT2B4 and UGT2B7 by more than 100-fold and 400-fold respectively. As in human liver microsomes the pathway is not Michaelis-Menten: the Eadie-Hofstee plot shows allosteric activation at low substrate and substrate inhibition at high substrate, so the authors fitted a combination of the Hill equation with the two-site substrate-inhibitory equation of Shou et al. The substrate sits in a static incubation compartment (initial-rate assay, 20 min, 0.5 mg protein per mL, 0-800 uM midazolam) and MDZG accumulates at the instantaneous rate v; the model observable is that formation rate in pmol/min/mg protein. The rate rises to a maximum of about 217 pmol/min/mg near 120 uM midazolam and then falls back toward Vmax * beta at high substrate. Sibling model: Hyland_2009_midazolam_hlm, the same equation independently fitted to pooled human liver microsomes."
  reference <- "Hyland R, Osborne T, Payne A, Kempshall S, Logan YR, Ezzeddine K, Jones B. In vitro and in vivo glucuronidation of midazolam in humans. Br J Clin Pharmacol. 2009 Apr;67(4):445-454. doi:10.1111/j.1365-2125.2009.03386.x. PMID: 19371318. PMCID: PMC2679108. Model equation: Equation 2 and its symbol-definition block, Results, 'MDZG enzyme kinetics' (p. 449). Parameter estimates: Table 4, column 'rUGT1A4' (p. 450). The fixed value of beta: Results, 'MDZG enzyme kinetics' -- 'since initial fitting of the data indicated a very low beta value (<0.01), this value was set as a constant of 0.01, to allow for better fitting of other variables'. Assay conditions: Materials and methods, 'Microsomal incubations', subsection 'MDZG kinetics'. Enzyme selectivity screen: Results, 'Initial screening'. Fitted curve reproduced here: Figure 3, closed circles (UGT1A4)."
  vignette <- "Hyland_2009_midazolam_glucuronidation"

  units <- list(
    time = "min",
    dosing = "umol/L (midazolam concentration placed in the incubation)",
    concentration = "pmol/min/mg protein (MDZG formation rate)"
  )

  # Static in vitro incubation: the "dose" is the midazolam concentration
  # placed directly into the `midazolam` state at time 0.
  dosing <- c("midazolam")

  paper_specific_compartments <- c("midazolam", "mdzg")

  compartmentData <- list(
    midazolam = list(analyte = "midazolam", units = "umol/L", specimen = "administration site", verified = TRUE),
    mdzg = list(analyte = "midazolam N-glucuronide", units = "pmol/mg protein", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "in vitro (recombinant human UGT1A4 Supersomes)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    system = "Recombinant human UGT1A4 Supersomes (BD Biosciences) at 0.5 mg protein/mL",
    matrix = "50 mM Tris-HCl pH 7.4 with 5 mM saccharolactone, alamethacin at 50 ug/mg protein and 10 mM MgCl2; the mixture was held on ice for 15 min so alamethacin could form pores in the membrane before midazolam was added",
    initiation = "Warmed to 37 C and initiated with 5 mM UDP-glucuronic acid (UDPGA)",
    temperature = "37 C",
    incubation_time = "20 min (initial-rate conditions), terminated with 3 volumes of ice-cold acetonitrile",
    substrate_range = "0-800 umol/L midazolam, n = 3 per concentration; 18 substrate concentrations were used because the equation has five unknown variables",
    replicates = "triplicate",
    quantification = "HPLC-MS/MS (Sciex API 3000) against a calibration curve prepared from an authentic MDZG standard isolated from HLM incubations and confirmed by 1H NMR and 1H-13C gHMBC to be conjugated on the alpha-nitrogen of the imidazole ring (90.3% purity)",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "Curve fitting was done in Grafit v5; the two-site substrate-inhibitory cooperative binding model (Equation 2) was compared against equation 7 of Houston and Kenworthy by F-test and gave the better fit (P < 0.001 and lower chi-squared). R-squared for this rUGT1A4 fit was 0.986 (Table 4).",
      "The Hill coefficient of 2.6 for the recombinant enzyme (2.3 in microsomes) is discussed as evidence of large conformational change and possible quaternary structure: the paper compares it with haemoglobin, whose coefficient of 2.7 accompanies four subunits.",
      "Substrate inhibition of UGT1A4 was not reported in the earlier midazolam glucuronidation work of Klieber et al.; the authors attribute the difference to studying substrate concentrations up to an order of magnitude higher, and note that substrate inhibition is unlikely to matter in vivo because plasma midazolam does not exceed 0.85 uM even after a 15 mg dose.",
      "Ketoconazole and itraconazole inhibited MDZG formation only weakly, with IC50 values of 150 +/- 6 uM and 308 +/- 14 uM (Figure 4); trifluoperazine activated the pathway to a maximum of 138% at 60 uM before inhibiting it. None of those three effects is encoded in this model -- see the vignette Errata for why."
    )
  )

  ini({
    # ---- Table 4, column 'rUGT1A4' --------------------------------------
    # Vmax, Ks and Ki are carried on the natural-log scale per the library
    # convention; alpha and n are dimensionless shape factors on the natural
    # scale. Every value is an estimate from the Grafit fit, reported +/- its
    # standard error, except beta.
    lvmax <- log(427); label("Log maximum catalytic rate at full saturation (pmol/min/mg protein)")  # Table 4, rUGT1A4: Vmax = 427 +/- 0.5 pmol min-1 mg-1
    lks <- log(64); label("Log dissociation constant for productive substrate binding to the enzyme (umol/L)")  # Table 4, rUGT1A4: Ks = 64 +/- 6.4 uM. Called "Km (64 uM)" in the Abstract
    lki <- log(79); label("Log dissociation constant for inhibitory substrate binding to the enzyme (umol/L)")  # Table 4, rUGT1A4: Ki = 79 +/- 8 uM
    alpha <- 14.9; label("Ks modifying factor for binding at the second site (unitless)")  # Table 4, rUGT1A4: a = 14.9 +/- 0.0008
    nhill <- 2.6; label("Hill coefficient, the slope of the allosteric effect (unitless)")  # Table 4, rUGT1A4: n = 2.6 +/- 0.07

    # The one parameter the authors did NOT estimate; beta was fixed for both
    # enzyme systems.
    beta <- fixed(0.01); label("Turnover rate modifying factor for the doubly-occupied enzyme (unitless)")  # Results, MDZG enzyme kinetics: "initial fitting of the data indicated a very low b value (<0.01), this value was set as a constant of 0.01"

    # ---- Residual error ---------------------------------------------------
    # The paper reports the fit quality as R-squared = 0.986 (Table 4) and
    # gives no residual standard deviation. Held at 0 rather than invented;
    # see the vignette Errata.
    addSd <- fixed(0); label("Additive residual SD on the MDZG formation rate (pmol/min/mg protein); magnitude not reported")  # Not reported: Table 4 gives R-squared only
  })

  model({
    # ---- 1. Back-transform ------------------------------------------------
    vmax <- exp(lvmax)
    ks <- exp(lks)
    ki <- exp(lki)

    # ---- 2. Equation 2 ----------------------------------------------------
    # DEVIATION FROM THE EQUATION AS TYPESET, identical to the one documented
    # at length in the sibling file Hyland_2009_midazolam_hlm.R. Equation 2 is
    # printed with Ks and Ki carrying no exponent, which (a) is dimensionally
    # inconsistent -- the denominator adds 1/[S]^n in uM^-n to 1/Ks in uM^-1,
    # commensurable only when n = 1 -- and (b) does not reproduce the paper's
    # own Figure 3. With the Table 4 rUGT1A4 values the form as typeset peaks
    # at 188 pmol/min/mg at 8.7 uM midazolam and has fallen to 7 pmol/min/mg by
    # 100 uM, whereas Figure 3 shows the fitted rUGT1A4 curve peaking near
    # 217 pmol/min/mg at about 110-120 uM. Raising Ks and Ki to the n-th power
    # -- the standard Hill generalisation of a dissociation constant, and
    # exactly the Shou et al. two-site model with [S] replaced by [S]^n that
    # Materials and methods cites as the source -- is the unique
    # dimensionally-consistent assignment and gives 217 pmol/min/mg at 120 uM.
    # The vignette gates this against the published curve for both enzyme
    # systems.
    #
    # Written below multiplied through by [S]^n, which is algebraically
    # identical and well behaved at [S] = 0.
    ksn <- ks^nhill
    kin <- ki^nhill
    sn <- midazolam^nhill

    v <- vmax * (sn / ksn + beta * sn * sn / (alpha * ksn * kin)) /
      (1 + sn / ksn + sn / kin + sn * sn / (alpha * ksn * kin))

    # ---- 3. The incubation ------------------------------------------------
    # Initial-rate assay: the 20 min incubation consumes only a small fraction
    # of the substrate, so the concentration is held constant, exactly as the
    # published rate law assumes. The fraction turned over is
    # v * 20 min * 0.5 mg protein/mL / [S], which reduces to 0.01 * v/[S] and
    # stays below 4% over the studied 10-800 umol/L range. The vignette
    # computes this explicitly rather than asserting it.
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
