Stites_2018_kras_g12c_qsp <- function() {
  description <- "QSP. In vitro / in silico (human KRAS G12C mutant cancer cell network). Mechanistic 14-state model of the Ras nucleotide-cycling network -- GEF-catalysed exchange, GAP-accelerated and intrinsic GTP hydrolysis, nucleotide association and dissociation, effector binding, and first-order Ras protein production and degradation -- carrying both wild-type Ras (KRAS + NRAS + HRAS) and the KRAS G12C mutant, extended with the two classes of covalent G12C inhibitor. Nucleotide-pocket inhibitors (NPIs, e.g. SML-8-73-1) bind irreversibly to nucleotide-free G12C; switch-II-pocket inhibitors (SIIPIs, e.g. ARS-853 and Ostrem compound 12) bind irreversibly to GDP-bound and nucleotide-free G12C and block GTP loading. Both drug-bound forms remain subject to Ras degradation, which is what makes covalent-inhibitor efficacy turnover-limited. An optional GEF-loading arm lets Ras GEFs load an NPI into the nucleotide pocket in competition with GTP and GDP. Steady-state readouts are the fraction of total Ras in the GTP-bound state and the fraction of effector bound to RasGTP; secondary G12C mutations are simulated by scaling the intrinsic GTPase rate and the intrinsic nucleotide dissociation rates."
  reference <- paste(
    "Stites EC, Shaw AS.",
    "Quantitative Systems Pharmacology Analysis of KRAS G12C Covalent Inhibitors.",
    "CPT Pharmacometrics Syst Pharmacol. 2018;7(5):342-351. doi:10.1002/psp4.12291.",
    "PMCID PMC5980551.",
    "Every ODE, rate constant and derived quantity below is transcribed from the authors'",
    "own deposited MATLAB implementation, distributed as the article's Supplementary",
    "Information archive (PSP4-7-342-s001.zip, folder 'Model Files': G12CDrugModel.m holds",
    "the 14 differential equations, RasG12C_inhibition.m the parameter assignments and",
    "initial conditions, and GenerateFigure2.m through GenerateFigure6.m the published",
    "scenarios). Table 1 of the article tabulates the eleven parameters that are new in",
    "this paper; the remaining biochemical constants come from the authors' earlier Ras",
    "network model, Stites EC, Trampont PC, Ma Z, Ravichandran KS. Network analysis of",
    "oncogenic Ras activation in cancer. Science. 2007;318(5849):463-467.",
    "doi:10.1126/science.1144642, and are reproduced verbatim in the deposited code.",
    sep = " "
  )
  vignette <- "Stites_2018_kras_g12c_covalent_inhibitors"

  # The deposited implementation works entirely in SI base units: molar
  # concentrations and seconds. Rate constants are therefore 1/s (first order)
  # or 1/(M*s) (second order), exactly as Table 1 prints them. There are no
  # dosing events -- inhibitor exposure is a fixed total drug concentration set
  # through npi0 / siipi0, which is how the authors' driver scripts pass it.
  units <- list(time = "s", dosing = "M", concentration = "M")

  covariateData <- list()

  covariatesDataExcluded <- list()

  # Every state is a molecular species of the Ras network rather than a
  # pharmacokinetic compartment, so none of them map onto the canonical
  # depot / central / peripheral vocabulary.
  paper_specific_compartments <- c(
    "ras_gdp",
    "ras_gtp",
    "ras_free",
    "eff",
    "ras_gtp_eff",
    "g12c_gdp",
    "g12c_gtp",
    "g12c_free",
    "g12c_gtp_eff",
    "g12c_gdp_siipi",
    "g12c_free_siipi",
    "g12c_npi",
    "npi",
    "siipi"
  )

  compartmentData <- list(
    ras_gdp = list(
      analyte = "wild-type Ras (KRAS + NRAS + HRAS) bound to GDP",
      units = "M",
      specimen = "not applicable",
      verified = TRUE
    ),
    ras_gtp = list(
      analyte = "wild-type Ras bound to GTP (free, not effector-bound)",
      units = "M",
      specimen = "not applicable",
      verified = TRUE
    ),
    ras_free = list(
      analyte = "wild-type Ras with an empty nucleotide pocket",
      units = "M",
      specimen = "not applicable",
      verified = TRUE
    ),
    eff = list(
      analyte = "free Ras effector protein",
      units = "M",
      specimen = "not applicable",
      verified = TRUE
    ),
    ras_gtp_eff = list(
      analyte = "wild-type RasGTP:effector complex",
      units = "M",
      specimen = "not applicable",
      verified = TRUE
    ),
    g12c_gdp = list(
      analyte = "KRAS G12C mutant bound to GDP",
      units = "M",
      specimen = "not applicable",
      verified = TRUE
    ),
    g12c_gtp = list(
      analyte = "KRAS G12C mutant bound to GTP (free, not effector-bound)",
      units = "M",
      specimen = "not applicable",
      verified = TRUE
    ),
    g12c_free = list(
      analyte = "KRAS G12C mutant with an empty nucleotide pocket",
      units = "M",
      specimen = "not applicable",
      verified = TRUE
    ),
    g12c_gtp_eff = list(
      analyte = "KRAS G12C RasGTP:effector complex",
      units = "M",
      specimen = "not applicable",
      verified = TRUE
    ),
    g12c_gdp_siipi = list(
      analyte = "switch-II-pocket inhibitor covalently bound to GDP-bound KRAS G12C",
      units = "M",
      specimen = "not applicable",
      verified = TRUE
    ),
    g12c_free_siipi = list(
      analyte = "switch-II-pocket inhibitor covalently bound to nucleotide-free KRAS G12C",
      units = "M",
      specimen = "not applicable",
      verified = TRUE
    ),
    g12c_npi = list(
      analyte = "nucleotide-pocket inhibitor covalently bound to KRAS G12C",
      units = "M",
      specimen = "not applicable",
      verified = TRUE
    ),
    npi = list(
      analyte = "free nucleotide-pocket inhibitor",
      units = "M",
      specimen = "not applicable",
      verified = TRUE
    ),
    siipi = list(
      analyte = "free switch-II-pocket inhibitor",
      units = "M",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in silico (human KRAS G12C mutant cancer cell); parameterised from human and in vitro Ras biochemistry",
    n_subjects = 0L,
    n_studies = 0L,
    disease_state = "KRAS G12C mutant cancer (lung adenocarcinoma is the motivating indication)",
    dose_range = "Covalent inhibitor concentrations of 0 to 1e-5 M (0 to 10 micromolar) were simulated; the published dose-response grids use 0, 1e-10, 1e-9, 2e-9, 5e-9, 1e-8, 1.3e-8, 2e-8, 3e-8, 5e-8, 7e-8, 1e-7, 1.3e-7, 2e-7, 3e-7, 5e-7, 7e-7, 1e-6, 1.3e-6, 2e-6, 3e-6, 5e-6, 7e-6 and 1e-5 M.",
    regions = "Not applicable (no clinical study)",
    notes = "No subjects were studied. This is a purely mechanistic model: the paper fits nothing and estimates nothing, and every constant is a biochemical rate taken from the literature or from the authors' earlier Ras network model (Stites 2007 Science). The default cellular composition is the one used throughout the paper: 50 percent of total Ras is KRAS and one of the two KRAS alleles is mutated (Table 1, Percent_KRAS = 50 and Percent_Mut = 50), giving 25 percent of total Ras as KRAS G12C and 75 percent as wild-type Ras. Validation targets are the reported cell-culture potencies of ARS-853: roughly 50 percent inhibition of Ras signal near 1 micromolar and roughly 95 percent inhibition at 10 micromolar (Lito 2016 Science; Patricelli 2016 Cancer Discov reports IC50 of 1 to 2 micromolar)."
  )

  ini({
    # =====================================================================
    # 1. Scenario inputs. These are the knobs the authors' driver scripts
    #    pass to RasG12C_inhibition.m; they are the intended user-varying
    #    inputs and are deliberately NOT wrapped in fixed().
    # =====================================================================
    npi0 <- 0
    label("Total nucleotide-pocket-inhibitor concentration (M)")
    # RasG12C_inhibition.m argument concs(3); swept over the Figure 2 / 5 / 6 dose grid

    siipi0 <- 0
    label("Total switch-II-pocket-inhibitor concentration (M)")
    # RasG12C_inhibition.m argument concs(4); swept over the Figure 2 / 3 / 4 / 6 dose grid

    mut_frac <- 0.25
    label("Fraction of total cellular Ras that is the KRAS G12C mutant (unitless)")
    # Table 1 Percent_KRAS = 50 and Percent_Mut = 50 -> 0.5 * 0.5 = 0.25; concs(1) in every driver script

    wt_frac <- 0.75
    label("Fraction of total cellular Ras that is wild-type (KRAS + NRAS + HRAS) (unitless)")
    # complement of mut_frac; concs(2) = .75 in every driver script

    kon_drug <- 76
    label("Covalent on-rate of whichever inhibitor is present (1/(M*s))")
    # Table 1 'k on,SIIPI,ARS' = 76 (ARS-853, the default). Table 1 also gives 0.12 for
    # Ostrem compound 12 and 2e6 for the NPI; RasG12C_inhibition.m argument 'onrate'
    # feeds k(32) = kaDrug, which is shared by the NPI and SIIPI binding reactions.

    gefload <- 0
    label("GEF loading factor: 0 disables GEF-mediated NPI loading, 1 loads the NPI in proportion to its abundance, >1 favours the NPI over nucleotide (unitless)")
    # RasG12C_inhibition.m argument 'GEFload'; Figure 5b uses 0, 1, 10 and 100, Figure 6 uses 4

    rtkfact <- 1
    label("Multiplier on basal GEF activity, used to model receptor-tyrosine-kinase activation or inhibition (unitless)")
    # RasG12C_inhibition.m argument 'RTKfact'; Figure 3 uses 0.1, 1 and 10, Figure 3cd uses 0.2, 1 and 5, Figure 6c uses 5

    fold_gtpase_impair <- 1
    label("Fold impairment of the G12C intrinsic GTPase rate caused by a secondary mutation (unitless)")
    # RasG12C_inhibition.m argument changes(1); Figure 4a uses 10 and 100, Figure 6b uses 100

    fold_cycling <- 1
    label("Fold increase in the G12C intrinsic nucleotide dissociation rates caused by a secondary mutation (unitless)")
    # RasG12C_inhibition.m argument changes(2); Figure 4b uses 10 and 100, Figure 6a uses 100

    tburn <- 0
    label("Time at which the inhibitor is allowed to engage the network (s); set > 0 to equilibrate drug-free first")
    # RasG12C_inhibition_timecourse.m equilibrates with y(13) = y(14) = 0 before setting the drug

    # =====================================================================
    # 2. Cellular abundances and system constants.
    #    Values are hardcoded at the top of RasG12C_inhibition.m; GTP and
    #    GDP are also printed in Table 1.
    # =====================================================================
    gtot <- fixed(4e-7)
    label("Total cellular Ras concentration, wild-type plus mutant (M)")
    # RasG12C_inhibition.m: GTot = 4e-7

    efftot <- fixed(4e-7)
    label("Total cellular Ras effector protein concentration (M)")
    # RasG12C_inhibition.m: EffTot = 4e-7

    gef <- fixed(2e-10)
    label("Basally active Ras GEF concentration (M)")
    # RasG12C_inhibition.m: GEF = 2e-10, 'fit parameter from original Ras model'

    gap <- fixed(6e-11)
    label("Basally active Ras GAP concentration (M)")
    # RasG12C_inhibition.m: GAP = 6e-11, 'fit parameter from original Ras model'

    gtp <- fixed(1.8e-4)
    label("Cellular GTP concentration (M)")
    # Table 1 GTP = 1.8e-4; RasG12C_inhibition.m: GTP = 180e-6

    gdp <- fixed(1.8e-5)
    label("Cellular GDP concentration (M)")
    # Table 1 GDP = 1.8e-5; RasG12C_inhibition.m: GDP = 18e-6

    kd_eff <- fixed(8e-8)
    label("Dissociation constant of the wild-type RasGTP:effector interaction (M)")
    # RasG12C_inhibition.m: Kd = 80e-9

    volscale <- fixed(250)
    label("Membrane-localisation volume scaling factor applied to every Michaelis constant (unitless)")
    # RasG12C_inhibition.m: volscale = 250, 'value and derivation from work of Kholodenko'

    # =====================================================================
    # 3. Wild-type Ras biochemical rate constants (Stites 2007 Science,
    #    reproduced verbatim in RasG12C_inhibition.m).
    # =====================================================================
    kgtpase <- fixed(3.5e-4)
    label("Intrinsic GTP hydrolysis rate constant of wild-type RasGTP (1/s)")
    # RasG12C_inhibition.m: kint = 3.5e-4

    kdiss_gdp <- fixed(1.1e-4)
    label("Dissociation rate constant of GDP from wild-type RasGDP (1/s)")
    # RasG12C_inhibition.m: kdissD = 1.1e-4

    kdiss_gtp <- fixed(2.5e-4)
    label("Dissociation rate constant of GTP from wild-type RasGTP (1/s)")
    # RasG12C_inhibition.m: kdissT = 2.5e-4

    kass_gdp <- fixed(2.3e6)
    label("Association rate constant of GDP with nucleotide-free Ras (1/(M*s))")
    # RasG12C_inhibition.m: kassD = 2.3e6, multiplied by GDP to give the pseudo-first-order k(10)

    kass_gtp <- fixed(2.2e6)
    label("Association rate constant of GTP with nucleotide-free Ras (1/(M*s))")
    # RasG12C_inhibition.m: kassT = 2.2e6, multiplied by GTP to give the pseudo-first-order k(11)

    kcat_gap <- fixed(5.4)
    label("Catalytic rate constant of the GAP reaction converting RasGTP to RasGDP (1/s)")
    # RasG12C_inhibition.m: kcat = 5.4

    km_gap <- fixed(2.3e-7)
    label("Michaelis constant of the GAP reaction, before volume scaling (M)")
    # RasG12C_inhibition.m: Km = .23e-6/volscale; the pre-scaling value is recorded here

    kcat_gef_gdp <- fixed(3.9)
    label("Catalytic rate constant of the GEF reaction converting RasGDP to RasGTP (1/s)")
    # RasG12C_inhibition.m: kD = 3.9 (and mkD = 3.9 for the mutant)

    km_gef_gdp <- fixed(3.86e-4)
    label("Michaelis constant of the GEF reaction on RasGDP, before volume scaling (M)")
    # RasG12C_inhibition.m: KmD = 3.86e-4/volscale; the pre-scaling value is recorded here

    km_gef_gtp <- fixed(3e-4)
    label("Michaelis constant of the GEF reaction on RasGTP, before volume scaling (M)")
    # RasG12C_inhibition.m: KmT = 3e-4/volscale; the pre-scaling value is recorded here

    kass_eff <- fixed(4.5e7)
    label("Association rate constant of effector protein with RasGTP (1/(M*s))")
    # RasG12C_inhibition.m: kassEff = 4.5e7

    # =====================================================================
    # 4. Protein turnover and the G12C mutant scaling factors.
    # =====================================================================
    kdeg <- fixed(8e-6)
    label("First-order degradation rate constant of Ras in every form (1/s); set to 0 to reproduce the no-turnover model")
    # Table 1 'k deg' = 8e-6 /s, a half-life of about 24 h; RasG12C_inhibition.m: kDeg = 8e-6

    scale_gtpase_g12c <- fixed(49 / 68)
    label("Intrinsic GTPase rate of KRAS G12C relative to wild-type Ras (unitless)")
    # Table 1 'Scalingfactor GTPase,G12C' = 72 percent; RasG12C_inhibition.m: mkint = (49/68)*kint

    scale_kd_eff_g12c <- fixed(67 / 56)
    label("Effector dissociation constant of KRAS G12C relative to wild-type Ras (unitless)")
    # Table 1 'Scalingfactor kdissG12C' = 120 percent; RasG12C_inhibition.m: mKd = (67/56)*Kd
  })

  model({
    # -------------------------------------------------------------------
    # 1. Cellular Ras pools. RasG12C_inhibition.m lines
    #    MutTot = Mutconc*GTot and WTRasTot = WTconc*GTot.
    # -------------------------------------------------------------------
    wt_tot <- wt_frac * gtot
    mut_tot <- mut_frac * gtot

    # -------------------------------------------------------------------
    # 2. Wild-type Ras rate constants. The nucleotide association rates are
    #    collapsed to pseudo-first-order constants by multiplying in the
    #    cellular nucleotide abundance (k(10) and k(11) of the deposited
    #    code), and every Michaelis constant is divided by the membrane
    #    volume scaling factor. The GEF reverse catalytic rate is not an
    #    independent constant: it is fixed by the Haldane relationship for
    #    the intrinsic nucleotide reactions.
    # -------------------------------------------------------------------
    kassd <- kass_gdp * gdp
    kasst <- kass_gtp * gtp
    haldane <- (kassd * kdiss_gtp) / (kdiss_gdp * kasst)

    kmd <- km_gef_gdp / volscale
    kmt <- km_gef_gtp / volscale
    kmg <- km_gap / volscale
    kcat_gef_gtp <- kcat_gef_gdp * kmt * haldane / kmd

    vmaxd <- gef * rtkfact * kcat_gef_gdp
    vmaxt <- gef * rtkfact * kcat_gef_gtp
    vmax_gap <- gap * kcat_gap
    kdiss_eff <- kass_eff * kd_eff

    # -------------------------------------------------------------------
    # 3. KRAS G12C rate constants. Only two biochemical properties differ
    #    from wild type -- the intrinsic GTPase rate and the effector
    #    affinity (Hunter 2015 Mol Cancer Res, via Table 1). A secondary
    #    mutation further divides the GTPase rate by fold_gtpase_impair and
    #    multiplies both intrinsic nucleotide dissociation rates by
    #    fold_cycling. Note that the G12C GAP-catalysed rate is set equal to
    #    the intrinsic GTPase rate (mkcat = mkint in the deposited code),
    #    which is how GAP insensitivity of codon-12 mutants is encoded.
    # -------------------------------------------------------------------
    mkgtpase <- scale_gtpase_g12c * kgtpase / fold_gtpase_impair
    mkdissd <- kdiss_gdp * fold_cycling
    mkdisst <- kdiss_gtp * fold_cycling
    mkassd <- kassd
    mkasst <- kasst

    mkd_eff <- scale_kd_eff_g12c * kd_eff
    mkdiss_eff <- mkd_eff * kass_eff

    mhaldane <- (mkassd * mkdisst) / (mkdissd * mkasst)
    mkmd <- km_gef_gdp / volscale
    mkmt <- km_gef_gtp / volscale
    mkcat_gef_gtp <- kcat_gef_gdp * mkmt * mhaldane / mkmd

    vmaxdv <- gef * rtkfact * kcat_gef_gdp
    vmaxtv <- gef * rtkfact * mkcat_gef_gtp
    vmax_gap_v <- gap * mkgtpase
    kmgv <- km_gap / volscale

    # -------------------------------------------------------------------
    # 4. Ras production. Set so that the drug-free steady-state total of
    #    each species equals its specified cellular abundance (k(30) and
    #    k(31) of the deposited code). Both vanish when kdeg is set to 0,
    #    which is exactly the no-turnover variant the paper compares
    #    against in Figure 2 (RasG12C_inhibition_NT.m sets k(29), k(30) and
    #    k(31) all to 0).
    # -------------------------------------------------------------------
    kpro_wt <- kdeg * wt_tot
    kpro_mut <- kdeg * mut_tot

    # -------------------------------------------------------------------
    # 5. GEF loading of the nucleotide-pocket inhibitor. The deposited code
    #    computes these two maximal velocities once, from the NOMINAL total
    #    NPI concentration passed to the driver (k(33) and k(34) are set
    #    from the scalar DrugConcNPI, not from the time-varying free-drug
    #    state), so they are held constant here as well.
    # -------------------------------------------------------------------
    vmaxd_npi <- vmaxdv * npi0 / gtp * gefload
    vmaxt_npi <- vmaxtv * npi0 / gdp * gefload

    # -------------------------------------------------------------------
    # 6. Initial conditions (RasG12C_inhibition.m y(1) through y(14)). All
    #    Ras starts GDP-bound, all effector free, and the inhibitor pools
    #    start at their total concentrations.
    # -------------------------------------------------------------------
    ras_gdp(0) <- wt_tot
    eff(0) <- efftot
    g12c_gdp(0) <- mut_tot
    npi(0) <- npi0
    siipi(0) <- siipi0

    # -------------------------------------------------------------------
    # 7. Reaction fluxes. Numbering follows G12CDrugModel.m exactly,
    #    including the gaps the authors left (there is no R15 to R17, R19
    #    to R21, R31, R33, R36 or R40). drug_engaged holds every covalent
    #    binding reaction off until tburn so that a time course can start
    #    from the drug-free steady state.
    # -------------------------------------------------------------------
    drug_engaged <- (t >= tburn)

    # Shared GEF competition denominator: wild-type and mutant Ras compete
    # for the same limited pool of GEF.
    gef_denom <- 1 + ras_gdp / kmd + ras_gtp / kmt + g12c_gdp / mkmd + g12c_gtp / mkmt

    r1 <- (vmaxd * ras_gdp / kmd - vmaxt * ras_gtp / kmt) / gef_denom
    r2 <- vmax_gap * ras_gtp / (kmg * (1 + g12c_gtp / kmgv) + ras_gtp)
    r3 <- kgtpase * ras_gtp
    r4 <- kdiss_gdp * ras_gdp - kassd * ras_free
    r5 <- kdiss_gtp * ras_gtp - kasst * ras_free
    r6 <- kass_eff * ras_gtp * eff - kdiss_eff * ras_gtp_eff
    r7 <- kgtpase * ras_gtp_eff

    r8 <- (vmaxdv * g12c_gdp / mkmd - vmaxtv * g12c_gtp / mkmt) / gef_denom
    r9 <- vmax_gap_v * g12c_gtp / (kmgv * (1 + ras_gtp / kmg) + g12c_gtp)
    r10 <- mkgtpase * g12c_gtp
    r11 <- mkdissd * g12c_gdp - mkassd * g12c_free
    r12 <- mkdisst * g12c_gtp - mkasst * g12c_free
    r13 <- kass_eff * g12c_gtp * eff - mkdiss_eff * g12c_gtp_eff
    r14 <- mkgtpase * g12c_gtp_eff

    # Nucleotide exchange on SIIPI-bound G12C: GDP may still leave and
    # rebind (k(27) = k(21) and k(28) = k(23)), but GTP may not load.
    r18 <- mkdissd * g12c_gdp_siipi - mkassd * g12c_free_siipi

    # Degradation of every Ras-containing species
    r22 <- kdeg * ras_gdp
    r23 <- kdeg * ras_gtp
    r24 <- kdeg * ras_free
    r25 <- kdeg * ras_gtp_eff
    r26 <- kdeg * g12c_gdp
    r27 <- kdeg * g12c_gtp
    r28 <- kdeg * g12c_free
    r29 <- kdeg * g12c_gtp_eff
    r30 <- kdeg * g12c_gdp_siipi
    r32 <- kdeg * g12c_free_siipi
    r38 <- kdeg * g12c_npi

    # Zero-order Ras production into the nucleotide-free pools
    r34 <- kpro_wt
    r35 <- kpro_mut

    # Covalent binding. The NPI binds nucleotide-free G12C only; the SIIPI
    # binds both the GDP-bound and the nucleotide-free forms.
    r37 <- kon_drug * g12c_free * npi * drug_engaged
    r39 <- kon_drug * g12c_gdp * siipi * drug_engaged
    r41 <- kon_drug * g12c_free * siipi * drug_engaged

    # GEF-mediated loading of the NPI, in competition with nucleotide
    r42 <- (vmaxd_npi * g12c_gdp / mkmd) / gef_denom * drug_engaged
    r43 <- (vmaxt_npi * g12c_gtp / mkmt) / gef_denom * drug_engaged

    # -------------------------------------------------------------------
    # 8. The 14 differential equations, transcribed from the dydt vector of
    #    G12CDrugModel.m in the same order. Effector released when a
    #    RasGTP:effector complex degrades returns to the free effector pool
    #    (the +r25 and +r29 terms), so total effector is conserved. Drug
    #    released when a drug-bound Ras degrades likewise returns to the
    #    free drug pool (the +r38, +r30 and +r32 terms), so total NPI and
    #    total SIIPI are each conserved.
    # -------------------------------------------------------------------
    d/dt(ras_gdp) <- -r1 + r2 + r3 - r4 + r7 - r22
    d/dt(ras_gtp) <- r1 - r2 - r3 - r5 - r6 - r23
    d/dt(ras_free) <- r4 + r5 - r24 + r34
    d/dt(eff) <- -r6 + r7 - r13 + r14 + r25 + r29
    d/dt(ras_gtp_eff) <- r6 - r7 - r25

    d/dt(g12c_gdp) <- -r8 + r9 + r10 - r11 + r14 - r26 - r39 - r42
    d/dt(g12c_gtp) <- r8 - r9 - r10 - r12 - r13 - r27 - r43
    d/dt(g12c_free) <- r11 + r12 - r28 - r37 + r35 - r41
    d/dt(g12c_gtp_eff) <- r13 - r14 - r29

    d/dt(g12c_gdp_siipi) <- -r18 - r30 + r39
    d/dt(g12c_free_siipi) <- r18 - r32 + r41
    d/dt(g12c_npi) <- r37 - r38 + r42 + r43

    d/dt(npi) <- -r37 + r38 - r42 - r43
    d/dt(siipi) <- -r39 - r41 + r30 + r32

    # -------------------------------------------------------------------
    # 9. Readouts. The two measures of Ras pathway signalling the paper
    #    uses are the fraction of total Ras in the GTP-bound state and the
    #    fraction of total effector bound to RasGTP (RasG12C_inhibition.m
    #    return values fractRasact and fractEffbound). Published dose
    #    responses plot each of these normalised to its own drug-free
    #    value; total_rasgtp is the un-normalised molar quantity plotted in
    #    Figure 3b and Figure 3d.
    # -------------------------------------------------------------------
    total_rasgtp <- ras_gtp + g12c_gtp + ras_gtp_eff + g12c_gtp_eff
    total_ras <- wt_tot + mut_tot
    frac_ras_active <- total_rasgtp / total_ras
    frac_eff_bound <- (ras_gtp_eff + g12c_gtp_eff) / efftot
    frac_eff_bound_wt <- ras_gtp_eff / efftot
    frac_eff_bound_mut <- g12c_gtp_eff / efftot
  })
}
