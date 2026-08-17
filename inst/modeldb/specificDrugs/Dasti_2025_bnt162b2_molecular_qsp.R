Dasti_2025_bnt162b2_molecular_qsp <- function() {
  description <- paste(
    "QSP (multiscale mRNA-vaccine immune response, molecular layer). Dasti 2025",
    "single-antigen-presenting-cell model of what happens inside one dendritic",
    "cell after it internalises an LNP-encapsulated BNT162b2 mRNA particle:",
    "endosomal escape of mRNA into the cytosol, ribosomal translation into",
    "spike antigen, proteolytic routing of antigen to peptides, competitive",
    "binding of those peptides to internalised MHC-II, and recycling of free",
    "and peptide-loaded MHC-II between the endosomal pool and the plasma",
    "membrane. 8 ODEs on a minutes timescale. The observable is the number of",
    "peptide-MHC-II complexes displayed on the plasma membrane; the rise and",
    "fall of that exposure curve is what sets the low / medium / high",
    "antigen-exposure transition rates used by the tissue layer",
    "(Dasti_2025_bnt162b2_qsp and its over-60 and mRNA-1273 variants), via the",
    "threshold-crossing construction of Supporting Information Section S3.",
    "There is no dosing event: the initial endosomal mRNA concentration",
    "represents the mRNA delivered by a single internalised LNP. Deterministic",
    "mechanism model with no IIV and no residual error.",
    sep = " "
  )
  reference <- paste(
    "Dasti L, Giampiccolo S, Pettina E, Fiandaca G, Zangani N, Leonardelli L,",
    "De Lima Hedayioglu F, Campanile E, Marchetti L. A Multiscale Quantitative",
    "Systems Pharmacology Model for the Development and Optimization of mRNA",
    "Vaccines. CPT Pharmacometrics Syst Pharmacol. 2025;14(7):1213-1224.",
    "doi:10.1002/psp4.70041.",
    "Molecular-layer state variables from Supporting Information Table S4,",
    "initial values from Table S5, parameters from Table S6, equations from",
    "Supporting Information Section S2.4 (equations 91-98). Cross-checked",
    "against the authors' published MATLAB implementation at",
    "https://github.com/cosbi-research/QSPmRNAVaccines (molecular_layer/).",
    sep = " "
  )
  vignette <- "Dasti_2025_mrna_vaccines"
  units <- list(
    time          = "min",
    dosing        = "not applicable (mRNA enters as the endosomal initial condition, M)",
    concentration = "molecules per cell (peptide-MHC-II complexes on the plasma membrane)"
  )
  # Every state is a paper-mechanistic intracellular species.
  paper_specific_compartments <- c(
    "mRNAe", "mRNAc", "Ag", "P",
    "MHCunPM", "MHCunINT", "MHCbPM", "MHCbINT"
  )
  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Descriptions transcribed from SI Table S4;
  # verified = TRUE means checked against that table. Sub-cellular sites are
  # mapped onto the controlled specimen vocabulary; the analyte text carries
  # the precise location.
  compartmentData <- list(
    mRNAe    = list(analyte = "mRNA in the endosome", units = "M", specimen = "endosome", verified = TRUE),
    mRNAc    = list(analyte = "mRNA in the cytosol of an antigen-presenting cell", units = "M", specimen = "tissue", verified = TRUE),
    Ag       = list(analyte = "translated antigen protein in the cytosol of an antigen-presenting cell", units = "M", specimen = "tissue", verified = TRUE),
    P        = list(analyte = "free antigenic peptide in the cytosol of an antigen-presenting cell", units = "M", specimen = "tissue", verified = TRUE),
    MHCunPM  = list(analyte = "unbound MHC-II on the antigen-presenting cell plasma membrane", units = "molecules", specimen = "tissue", verified = TRUE),
    MHCunINT = list(analyte = "unbound internalized MHC-II", units = "molecules", specimen = "endosome", verified = TRUE),
    MHCbINT  = list(analyte = "internalized MHC-II/peptide complex", units = "molecules", specimen = "endosome", verified = TRUE),
    MHCbPM   = list(analyte = "MHC-II/peptide complex on the antigen-presenting cell plasma membrane", units = "molecules", specimen = "tissue", verified = TRUE)
  )
  covariateData <- list()
  population <- list(
    species       = "human (single antigen-presenting cell; intracellular model)",
    disease_state = "not applicable (subcellular antigen-processing model)",
    vaccine       = "Pfizer-BioNTech BNT162b2",
    n_subjects    = NA_integer_,
    notes         = paste(
      "Not fitted to data. Rate constants were taken from the antigen-",
      "presentation literature (Table S6 references 19-27) except the antigen",
      "production rate kp, which the authors computed for the BNT162b2 mRNA",
      "sequence assuming elongation-limited translation and HEK293 tRNA",
      "abundances following Chu 2014. Recalibrating for another vaccine would",
      "require a sequence-specific kp.",
      sep = " "
    )
  )

  ini({
    kesc     <- fixed(8.3333e-06); label("mRNA endosomal escape rate (L/min)")                              # Table S6 (8.3333e-06 L/min); ref 21
    kdmRNA   <- fixed(log(2) / (10 * 60)); label("mRNA degradation rate in cytosol (1/min; 10-h half-life)") # Table S6 prints 0.0012/min, a 4-dp rounding of log(2)/(10*60) = 0.00115525; the exact form is used here and in the published code (ref 3)
    kp       <- fixed(15.72);      label("Antigen production rate via mRNA translation (1/min)")            # Table S6 (15.72); computed for the BNT162b2 mRNA sequence, SI eq. 92 note
    kdAg     <- fixed(1.00e-02);   label("Native antigen degradation rate (1/min)")                         # Table S6 (1.00E-02/min); ref 25
    kpept    <- fixed(1.20e-02);   label("Rate constant routing antigen to lysosomes (1/min)")              # Table S6 (1.20E-02/min); ref 22
    kaMHC    <- fixed(100);        label("MHC-II/peptide association rate constant (1/(M*min))")     # Table S6 (100 M-1min-1); ref 22
    kdMHC    <- fixed(4.2000e-04); label("MHC-II/peptide dissociation rate constant (1/min)")        # Table S6 (4.2000e-04/min); ref 22
    kin      <- fixed(1.00e-02);   label("MHC-II internalisation rate constant (1/min)")                    # Table S6 (1.00E-02/min); ref 23
    kout     <- fixed(2.00e-02);   label("Vesicle recycling rate to the plasma membrane (1/min)")            # Table S6 (2.00E-02/min); ref 24
    ksyn     <- fixed(1.00e-03);   label("MHC-II synthesis (and matched degradation) rate constant (1/min)") # Table S6 prints 1.10E-03/min but the published code uses 1.0e-03; 1.0e-03 is used here because it reproduces the tissue-layer transition rates of Table S3 more closely and keeps the two layers self-consistent -- see vignette Errata (ref 24)
    Vc       <- fixed(1.00e-15);   label("Cell volume (L)")                                                 # Table S6 (1.00E-15 L); ref 26
    Ve       <- fixed(4.00e-17);   label("Endosomal volume (L)")                                            # Table S6 (4.00E-17 L); ref 27
    NAvo     <- fixed(6.022e23);   label("Avogadro constant (molecules/mole)")                              # SI Table S3 general parameters; ref 7
    mRNAe0   <- fixed(2.0757e-07); label("Initial endosomal mRNA concentration (M)")                        # Table S5 (2.0757e-07 M), calculated as mRNAmolmax/(Ve*NAvo) with mRNAmolmax = 5 molecules per LNP (Table S6, ref 20)
    MHCunPM0 <- fixed(2.00e+05);   label("Initial unbound MHC-II on the plasma membrane (molecules)")        # Table S5 (2.00e+05 molecules); ref 19
    expSd    <- fixed(0);          label("Exponential (log-scale) residual SD -- not applicable")           # deterministic mechanistic layer, not fitted to data; fixed at 0
  })

  model({
    # ==== mRNA release and translation (SI eqs. 91-93) ==================
    # kesc is a volumetric rate (L/min), so it is divided by the endosomal
    # volume for the efflux term and by the cytosolic volume for the influx
    # term, exactly as printed.
    d/dt(mRNAe) <- -kesc / Ve * mRNAe
    d/dt(mRNAc) <- kesc / Vc * mRNAe - kdmRNA * mRNAc
    d/dt(Ag)    <- kp * Vc / Ve * mRNAc - (kdAg + kpept) * Ag

    # ==== Free antigenic peptide in the endosome (SI eq. 94) ============
    # The MHC-II states are molecule counts, so they are converted to
    # endosomal molar concentration via /(Ve*NAvo) where they meet Ag/P.
    d/dt(P) <- kpept * Ag + kdMHC / Ve * MHCbINT / NAvo -
      kaMHC / Ve * MHCunINT / NAvo * P - kdAg * P

    # ==== MHC-II recycling and peptide loading (SI eqs. 95-98) ==========
    d/dt(MHCunPM)  <- kout * MHCunINT - kin * MHCunPM + kdMHC * MHCbPM -
      ksyn * MHCunPM
    d/dt(MHCunINT) <- kin * MHCunPM - kout * MHCunINT + kdMHC * MHCbINT -
      kaMHC * MHCunINT * P + ksyn * (MHCunPM + MHCbPM)
    d/dt(MHCbINT)  <- kaMHC * MHCunINT * P - kdMHC * MHCbINT - kout * MHCbINT
    d/dt(MHCbPM)   <- kout * MHCbINT - ksyn * MHCbPM - kdMHC * MHCbPM

    # ==== Initial conditions (SI Table S5) =============================
    mRNAe(0)   <- mRNAe0
    MHCunPM(0) <- MHCunPM0
    # mRNAc, Ag, P, MHCunINT, MHCbPM and MHCbINT all start at 0 (Table S5)

    # ==== Observation ==================================================
    # Peptide-MHC-II complexes displayed on the plasma membrane: the
    # "MHCII-antigen exposition curve" of SI Figure S1 and main-text
    # Figure 2c. Cc is the canonical single-output observation name; here
    # it is a molecule count per cell, not a drug concentration.
    Cc <- MHCbPM
    Cc ~ lnorm(expSd)
  })
}
