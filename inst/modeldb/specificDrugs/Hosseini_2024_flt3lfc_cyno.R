Hosseini_2024_flt3lfc_cyno <- function() {
  description <- "QSP. Preclinical (cynomolgus monkey). Minimal PBPK model (Cao/Jusko two-tissue-group mPBPK) for FLT3L-Fc (RO7497987), a half-life-extended effectorless Fc fusion of human FLT3 ligand, with expansion-enhanced target-mediated drug disposition: drug binds one or two monomeric FLT3 receptors in plasma, and the double-bound (homodimer) complex drives sigmoidal expansion of the total FLT3 receptor pool, which in turn amplifies TMDD. Calibrated to single- and repeat-IV-dose PK in cynomolgus monkeys over 0.1-10 mg/kg (Hosseini 2024 Table S2, cyno column)."
  reference   <- "Hosseini I, Fleisher B, Getz J, Decalf J, Kwong M, Ovacik M, Bainbridge TW, Moussion C, Rao GK, Gadkar K, Kamath AV, Ramanujan S. A Minimal PBPK/PD Model with Expansion-Enhanced Target-Mediated Drug Disposition to Support a First-in-Human Clinical Study Design for a FLT3L-Fc Molecule. Pharmaceutics. 2024 May 15;16(5):660. doi:10.3390/pharmaceutics16050660. PMCID PMC11125320. Structural equations from the Supplementary Materials 'Supplemental ODEs and Repeated Assignments' (SimBiology export); parameter values from Supplementary Table S2, 'Cyno Value' column."
  vignette    <- "Hosseini_2024_flt3lfc"

  # FLT3L-Fc is given IV straight into `plasma`. Declared explicitly because
  # the registry's automatic detection only recognises `depot` / `central`.
  dosing <- c("plasma")

  # complex_sb / complex_db are the single-bound (drug:FLT3) and double-bound
  # (FLT3:drug:FLT3 homodimer) receptor complexes of the expansion-enhanced
  # TMDD mechanism; the paper's CRp and CRpR. The canonical `complex` supplies
  # only one slot, so both carry explicit paper-mechanistic names.
  paper_specific_compartments <- c("complex_sb", "complex_db")

  units <- list(time = "h", dosing = "ug", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    plasma     = list(analyte = "FLT3L-Fc", units = "ug", specimen = "plasma", verified = TRUE),
    tight      = list(analyte = "FLT3L-Fc", units = "ug", specimen = "tissue", verified = TRUE),
    leaky      = list(analyte = "FLT3L-Fc", units = "ug", specimen = "tissue", verified = TRUE),
    lymph      = list(analyte = "FLT3L-Fc", units = "ug", specimen = "lymph", verified = TRUE),
    target     = list(analyte = "free FLT3 receptor", units = "nM", specimen = "plasma", verified = TRUE),
    complex_sb = list(analyte = "single-bound FLT3L-Fc:FLT3 complex", units = "nM", specimen = "plasma", verified = TRUE),
    complex_db = list(analyte = "double-bound FLT3:FLT3L-Fc:FLT3 homodimer complex", units = "nM", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species       = "cynomolgus monkey",
    n_subjects    = 9,
    n_studies     = 1,
    weight_range  = "2.6 kg nominal body weight for the species (Hosseini 2024 Table S2)",
    disease_state = "Healthy cynomolgus monkeys (Charles River Laboratories; IACUC protocols 20239581 and 5003546)",
    dose_range    = "Single IV 0.1 mg/kg; repeat IV 1 or 10 mg/kg on days 0 and 21 (preclinical study 1, n = 3 per group)",
    notes         = paste0(
      "Typical-value mechanistic (QSP) simulator: the paper calibrated the model ",
      "by particle-swarm optimisation in SimBiology/gQSPSim and reports point ",
      "estimates with no uncertainty, no inter-individual variability and no ",
      "residual-error model, so none are encoded. Anti-drug antibodies were ",
      "detected in all nine animals by day 14 and depressed late exposure; the ",
      "published model deliberately does not describe ADA, so simulated ",
      "concentrations correspond to the ADA-negative samples only. Preclinical ",
      "study 2 (repeat IV 1 or 3 mg/kg on days 0 and 21, n = 6) was held out as ",
      "an external validation set and was not used for calibration."
    ),
    model_class   = "QSP / minimal PBPK with expansion-enhanced TMDD (7 states)",
    n_states      = 7
  )

  ini({
    # ---- Species physiology (Hosseini 2024 Table S2, Cyno column) ----
    vplasma     <- fixed(0.09);    label("Plasma volume (L)")                                          # Table S2 Vplasma, cyno fit from data
    vlymph      <- fixed(0.086);   label("Lymph volume (L)")                                           # Table S2 Vlymph, cyno fit from data
    visf        <- fixed(0.579);   label("Total interstitial fluid volume (L)")                        # Table S2 ISF, cyno fit from data; table unit label 'mL/kg' is a typo -- the SimBiology project names the parameter ISF_L, i.e. litres
    kp          <- fixed(0.8);     label("Fraction of interstitial space available to IgG1 (unitless)") # Table S2 Kp, Cao 2013
    ltot        <- fixed(0.012);   label("Total lymph flow (L/h)")                                     # Table S2 L, cyno fit from data
    clp         <- fixed(4.2e-4);  label("Nonspecific clearance of FLT3L-Fc from plasma (L/h)")        # Table S2 CLp, cyno fit from data
    sigma_tight <- fixed(0.985);   label("Vascular reflection coefficient, tight tissues (unitless)")  # Table S2 sigma_tight, fit from data
    sigma_leaky <- fixed(0.656);   label("Vascular reflection coefficient, leaky tissues (unitless)")  # Table S2 sigma_leaky, fit from data
    sigma_lymph <- fixed(0.2);     label("Lymphatic reflection coefficient (unitless)")                # Table S2 sigma_lymph, Cao 2013

    # ---- FLT3 target turnover and binding (Table S2, Cyno column) ----
    rbase_target <- fixed(11.3);   label("Baseline free FLT3 receptor concentration in plasma (nM)")   # Table S2 CR,0, fit from data
    kdeg         <- fixed(0.016);  label("Turnover rate constant of FLT3 receptors and of all drug-receptor complexes (1/h)") # Table S2 kdeg, fit from data; supplement repeated assignment kint = kdeg_central
    kon1         <- fixed(0.1);    label("Association rate constant, drug to first FLT3 receptor (1/(nM*h))")  # Table S2 kon1, assumption
    kd1          <- fixed(9);      label("Equilibrium dissociation constant, drug to first FLT3 receptor (nM)") # Table S2 KD1, Verstraete 2011
    kon2         <- fixed(3.6);    label("Association rate constant, single-bound complex to second FLT3 receptor (1/(nM*h))") # Table S2 kon2, in-house Biacore
    kd2          <- fixed(0.2);    label("Equilibrium dissociation constant, single-bound complex to second FLT3 receptor (nM)") # Table S2 KD2, in-house Biacore

    # ---- Expansion of the total FLT3 receptor pool (Table S2, Cyno column) ----
    vmax_prolif <- fixed(1.90e-3); label("Maximum rate of FLT3 receptor pool expansion driven by the double-bound complex (1/h)") # Table S2 vmprolif, fit from data
    km_prolif   <- fixed(0.547);   label("Double-bound receptor occupancy giving half-maximal expansion (fraction of total receptor)") # Table S2 kmprolif, fit from data; the supplement ODE compares km_prolif with RO_DB/100, so it is on the fraction scale despite Table S2 wording it as a percentage
    hill_prolif <- fixed(3.93);    label("Hill coefficient of the receptor-expansion function (unitless)") # Table S2 alpha, fit from data
  })

  model({
    # Molecular weight of FLT3L-Fc: 83000 Da = 83 ug/nmol
    # (Table S2 MWFC, in-house measurement).
    mw_ug_nmol <- 83

    # Off-rates from the reported affinities. The supplement writes the
    # dissociation flux as kon*KD*complex, so koff = kon * KD.
    koff1 <- kon1 * kd1                                       # 0.9  1/h
    koff2 <- kon2 * kd2                                       # 0.72 1/h

    # Complex internalisation equals free-receptor turnover
    # (supplement repeated assignment kint_1h = kdeg_central_1h).
    kint <- kdeg

    # Zero-order receptor synthesis holding the free receptor at baseline in
    # the absence of drug (ksyn = kdeg * CR,0).
    ksyn <- kdeg * rbase_target

    # Cao 2013 mPBPK tissue partition. The tight tissue group holds 65% of the
    # available interstitial volume but receives only 33% of lymph flow; the
    # leaky group is the complement. Table S2's "fleaky = 0.65 / ftight = 0.35"
    # rows carry the volume split under a flow-fraction caption; the
    # supplement's SimBiology project settles it by storing these four
    # assignments verbatim -- "Vtight_L = 0.65*ISF_L*Kp",
    # "Vleaky_L = 0.35*ISF_L*Kp", "L_tight_Lh = 0.33*L_Lh" and
    # "L_leaky_Lh = 0.67*L_Lh" -- which is the Cao 2013 convention already used
    # by Yuan_2019_concizumab.
    vtight <- 0.65 * visf * kp
    vleaky <- 0.35 * visf * kp
    ltight <- 0.33 * ltot
    lleaky <- 0.67 * ltot

    # Concentrations from amounts
    cp_ugL     <- plasma / vplasma
    ctight_ugL <- tight / vtight
    cleaky_ugL <- leaky / vleaky
    clymph_ugL <- lymph / vlymph
    cp_nM      <- cp_ugL / mw_ug_nmol

    # Mass-action binding fluxes (nM/h)
    bind1 <- kon1 * cp_nM * target - koff1 * complex_sb
    bind2 <- kon2 * complex_sb * target - koff2 * complex_db

    # Receptor occupancy. Each double-bound complex sequesters two receptors,
    # so it is counted twice in the total pool (Eq 6).
    target_tot <- target + complex_sb + 2 * complex_db
    ro_sb <- complex_sb / max(target_tot, 1e-18) * 100                    # Eq 7
    ro_db <- 2 * complex_db / max(target_tot, 1e-18) * 100                # Eq 5
    ro    <- (1 - target / max(target_tot, 1e-18)) * 100                  # Eq 8

    # Fractional double-bound occupancy driving expansion. Floored at zero so
    # solver round-off below zero cannot raise a negative base to a fractional
    # power (the supplement applies the same guard as real(max(0, .)) in its
    # PD block).
    rodb_frac <- max(ro_db / 100, 0)
    expansion <- target_tot * vmax_prolif * rodb_frac^hill_prolif /
      (rodb_frac^hill_prolif + km_prolif^hill_prolif)

    # ---- Drug disposition, minimal PBPK (Eq 1 and supplement ODEs) ----
    d/dt(plasma) <- -bind1 * vplasma * mw_ug_nmol -
      (1 - sigma_tight) * ltight * cp_ugL -
      (1 - sigma_leaky) * lleaky * cp_ugL +
      ltot * clymph_ugL -
      clp * cp_ugL

    d/dt(tight) <- (1 - sigma_tight) * ltight * cp_ugL -
      (1 - sigma_lymph) * ltight * ctight_ugL

    d/dt(leaky) <- (1 - sigma_leaky) * lleaky * cp_ugL -
      (1 - sigma_lymph) * lleaky * cleaky_ugL

    d/dt(lymph) <- (1 - sigma_lymph) * ltight * ctight_ugL +
      (1 - sigma_lymph) * lleaky * cleaky_ugL -
      ltot * clymph_ugL

    # ---- Target and complexes in plasma (Eq 2-4) ----
    d/dt(target)     <- ksyn - kdeg * target - bind1 - bind2 + expansion
    d/dt(complex_sb) <- bind1 - bind2 - kint * complex_sb
    d/dt(complex_db) <- bind2 - kint * complex_db

    target(0) <- rbase_target

    # ---- Outputs ----
    Cc      <- cp_ugL / 1000                                  # ug/mL
    CcnM    <- cp_nM                                          # free FLT3L-Fc, nM
    RO      <- ro
    ROsb    <- ro_sb
    ROdb    <- ro_db
    Rtot    <- target_tot
  })
}
