Chen_2014_immunogenicity_qsp <- function() {
  description <- "QSP (multiscale immunogenicity). Chen 2014 theoretical mechanistic model of antidrug-antibody (ADA) formation against therapeutic proteins: subcellular antigen presentation on dendritic cells, cellular T-helper and B-cell activation / proliferation / differentiation, and whole-body antigen and ADA disposition. Encoded with 1 T-epitope (j=1), 1 lumped effective binding MHC-II (representing 3 human MHC-II binders of the paper demo, all sharing the same affinity), 1 lumped ADA affinity class (i=1) and 1-compartment antigen PK -- deterministic simulation with human parameter set. See vignette for the paper's polyclonal 17-clone / 6-MHC-II full form and the reasons for the reduction."
  reference <- paste(
    "Chen X, Hickling TP, Vicini P.",
    "A Mechanistic, Multiscale Mathematical Model of Immunogenicity",
    "for Therapeutic Proteins: Part 1 - Theoretical Model.",
    "CPT Pharmacometrics Syst Pharmacol. 2014;3(9):e133.",
    "doi:10.1038/psp.2014.30.",
    "State variables from Supplementary Table S1.",
    "Parameter values from Supplementary Table S2 (human column).",
    "Equations from Supplementary Materials sections 1-6",
    "(Word document psp4201430-sup-0004.doc, equations 1-29 with sub-equations).",
    sep = " "
  )
  vignette <- "Chen_2014_immunogenicity_qsp"
  units <- list(
    time          = "day",
    dosing        = "pmole (antigenic protein; convert from mg via molecular weight)",
    concentration = "pM (antigenic protein in plasma)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    MS   = list(analyte = "antigenic protein", units = NA_character_, specimen = "administration site", verified = FALSE),
    iDC  = list(analyte = "T-epitope presenting dendritic cells", units = NA_character_, specimen = "administration site", verified = FALSE),
    mDC  = list(analyte = "memory T-epitope presenting dendritic cells", units = NA_character_, specimen = "administration site", verified = FALSE),
    AgE  = list(analyte = "antigen presenting endosomes", units = NA_character_, specimen = "administration site", verified = FALSE),
    pE   = list(analyte = "processed antigen presenting endosomes", units = NA_character_, specimen = "administration site", verified = FALSE),
    cpE  = list(analyte = "cytoplasmic processed antigen presenting endosomes", units = NA_character_, specimen = "administration site", verified = FALSE),
    cptE = list(analyte = "cytoplasmic processed T-epitope presenting endosomes", units = NA_character_, specimen = "administration site", verified = FALSE),
    mhcE = list(analyte = "MHC-II molecules on antigen presenting endosomes", units = NA_character_, specimen = "administration site", verified = FALSE),
    pmE  = list(analyte = "processed antigen presenting macrophages", units = NA_character_, specimen = "administration site", verified = FALSE),
    cpmE = list(analyte = "cytoplasmic processed antigen presenting macrophages", units = NA_character_, specimen = "administration site", verified = FALSE),
    pmM  = list(analyte = "processed MHC-II molecules on macrophages", units = NA_character_, specimen = "administration site", verified = FALSE),
    cpmM = list(analyte = "cytoplasmic processed MHC-II molecules on macrophages", units = NA_character_, specimen = "administration site", verified = FALSE),
    mhcM = list(analyte = "MHC-II molecules on macrophages", units = NA_character_, specimen = "administration site", verified = FALSE),
    NT   = list(analyte = "naive T-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    aTn  = list(analyte = "activated naive T-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    aTm  = list(analyte = "memory T-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    MT   = list(analyte = "memory T-cells", units = NA_character_, specimen = "tumor", verified = FALSE),
    FT   = list(analyte = "free tumor cells", units = NA_character_, specimen = "tumor", verified = FALSE),
    NB   = list(analyte = "naive B-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    aBn  = list(analyte = "activated naive B-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    aBm  = list(analyte = "memory B-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    MB   = list(analyte = "memory B-cells", units = NA_character_, specimen = "tumor", verified = FALSE),
    SP   = list(analyte = "soluble protein", units = NA_character_, specimen = "plasma", verified = FALSE),
    LP   = list(analyte = "lymph node", units = NA_character_, specimen = "not applicable", verified = FALSE),
    Ab   = list(analyte = "antibody", units = NA_character_, specimen = "plasma", verified = FALSE),
    AgIS = list(analyte = "antigenic protein in depot", units = NA_character_, specimen = "administration site", verified = FALSE),
    Ag   = list(analyte = "antigenic protein", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "human (theoretical / hypothetical antigenic protein)",
    n_subjects     = 0,
    n_studies      = 0,
    disease_state  = "None (theoretical simulation, no clinical data fit). Parameters compiled from experimental immunology literature; see Supplementary Table S2 references.",
    dose_range     = "50 mg/kg IV in the paper's Figure 4 hypothetical-Ag simulation (MW = 150 kDa); dose sweep 3.3e1 to 3.3e10 pmole in the Figure 5 sensitivity analysis.",
    regions        = "n/a (no clinical trial)",
    notes          = "Deterministic mechanistic QSP; no IIV, no residual error. See paper Methods 'Simulation of immune responses in human against a hypothetical antigen'. Encoded parameters use the human column of Supplementary Table S2."
  )

  ini({
    # === All parameters fixed (deterministic QSP, no data fit; Table S2 human values) ===

    # 1) Maturation signal (LPS) & dendritic cells (paper section 1)
    betaMS  <- fixed(0.3696);   label("MS (LPS) elimination rate (1/day)")                        # Table S2 (Shi 2008 ref 2)
    betaID  <- fixed(0.0924);   label("Immature DC death & homeostatic production rate (1/day)")  # Table S2 (Merad 2009 ref 4)
    deltaID <- fixed(1.5);      label("Maximum immature DC activation rate (1/day)")              # Table S2 (Lee 2009 ref 5)
    KMS     <- fixed(9852);     label("MS concentration for half-max ID activation (ng/L)")       # Table S2 (Saito 2006 ref 6)
    betaMD  <- fixed(0.2310);   label("Mature DC death rate (1/day)")                             # Table S2 (Diao 2006 ref 7)
    iDC0    <- fixed(5e7);      label("Initial immature DC number (cells)")                       # Table S2 (Castiglione 2005 ref 26)

    # 2) Antigen presentation on DCs (paper section 2)
    alphaAgE  <- fixed(14.4);      label("Antigen internalization rate into DC endosome (1/day)")               # Table S2 (Agrawal 1996 ref 8)
    betaAgE   <- fixed(17.28);     label("Antigen degradation rate in endosome (1/day)")                        # Table S2 (Agrawal 1996 ref 8)
    betaPep   <- fixed(144);       label("Epitope-peptide degradation rate in endosome (1/day)")                # Table S2 (Agrawal 1996 ref 8)
    betaMHC   <- fixed(1.663);     label("MHC-II degradation rate (1/day)")                                     # Table S2 (Cella 1997 ref 9)
    betaPM    <- fixed(0.1663);    label("Peptide-MHC complex degradation rate (1/day)")                        # Table S2 (Cella 1997 ref 9)
    kExt      <- fixed(28.8);      label("Exocytosis rate of peptide-MHC complex (endosome -> membrane, 1/day)")# Table S2 (Agrawal 1996 ref 8)
    kIntCplx  <- fixed(14.4);      label("Internalization rate of peptide-MHC & Ag-BCR (membrane -> endosome, 1/day)")# Table S2 (Agrawal 1996 ref 8)
    Vendo     <- fixed(4e-16);     label("Endosome volume in one dendritic cell (L)")                           # Table S2 (Agrawal 1996 ref 8)
    konMHC    <- fixed(8.64e-3);   label("On-rate constant for peptide-MHC-II binding (1/(pM*day))")            # Table S2 (Foote 1995 ref 1)
    KpMN      <- fixed(400);       label("Number of peptide-MHC on DC for half-max naive T activation (dimensionless)")# Table S2 (Kimachi 1997 ref 10)
    KpMM      <- fixed(40);        label("Number of peptide-MHC on DC for half-max memory T activation (dimensionless)")# Table S2 (Kimachi 1997 ref 10)
    cp0pl     <- fixed(3.025e8);   label("Endogenous competing-protein amount in plasma (pmole)")               # Table S2 (Agrawal 1996 ref 8)
    koffC     <- fixed(34560);     label("Off-rate constant for competing peptide-MHC binding (1/day)")         # Table S2 (Agrawal 1996 ref 8)

    # 3) T helper cell kinetics (paper section 3)
    betaNT   <- fixed(0.0029);     label("Naive T-cell death rate (1/day)")                          # Table S2 (McCune 2000 ref 12)
    deltaNT  <- fixed(1.5);        label("Maximum naive T-cell activation rate (1/day)")             # Table S2 (Lee 2009 ref 5)
    rhoAT    <- fixed(0.5973);     label("Maximum activated T-cell proliferation rate (1/day)")      # Table S2 (Croft 1994 / Sundrud 2004 refs 14-15)
    betaAT   <- fixed(0.18);       label("Activated T-cell death rate (1/day)")                      # Table S2 (Kohler 2007 ref 16)
    deltaMT  <- fixed(1.5);        label("Maximum memory T-cell activation rate (1/day)")            # Table S2 (Lee 2009 ref 5)
    betaMT   <- fixed(2.7397e-4);  label("Memory T-cell death rate (1/day)")                         # Table S2 (Hammarlund 2003 ref 17)
    betaFT   <- fixed(0.18);       label("Functional T-cell death rate (1/day)")                     # Table S2 (Lee 2009 ref 5)
    fracT1   <- fixed(0.5);        label("Fraction of activated T cells differentiating to memory T cells (unitless)")# Table S2 (Bell 1970 ref 18)
    NT0      <- fixed(1445);       label("Initial naive T-cell number (cells, per T-epitope)")       # Table S2 (Castiglione 2005 / Delluc 2011 refs 26, 32)

    # 4) B cell kinetics (paper section 4)
    deltaNB  <- fixed(3);          label("Maximum naive B-cell activation rate (1/day)")             # Table S2 (Lee 2009 ref 5)
    deltaMB  <- fixed(3);          label("Maximum memory B-cell activation rate (1/day)")            # Table S2 (Lee 2009 ref 5)
    betaMB   <- fixed(7.83e-5);    label("Memory B-cell death rate (1/day)")                         # Table S2 (Crotty 2004 ref 25)
    betaSP   <- fixed(0.2310);     label("Short-lived plasma-cell death rate (1/day)")               # Table S2 (Castiglione 2005 ref 26)
    betaLP   <- fixed(0.0050);     label("Long-lived plasma-cell death rate (1/day)")                # Table S2 (Slifka 1998 ref 27)
    rhoABN   <- fixed(0.3333);     label("Maximum activated-B proliferation rate (naive-derived, 1/day)")# Table S2 (Tangye 2003 / Fecteau 2009 refs 22-23)
    rhoABM   <- fixed(0.7273);     label("Maximum activated-B proliferation rate (memory-derived, 1/day)")# Table S2 (Tangye 2003 / Fecteau 2009 refs 22-23)
    betaAB   <- fixed(0.2518);     label("Activated B-cell death rate (1/day)")                      # Table S2 (Lee 2009 / Fecteau 2009 refs 5, 23)
    fracB1   <- fixed(0.5);        label("Fraction of activated B cells differentiating to memory B (unitless)")# Table S2 (Bell 1970 ref 18)
    fracB2   <- fixed(0.4);        label("Fraction of activated B cells differentiating to short-lived plasma (unitless)")# Table S2 (data fitting per Table S2 note)
    CCN      <- fixed(10);         label("Carrying capacity of one FT cell for naive B activation (unitless)")# Table S2 (data fitting per Table S2 note)
    CCM      <- fixed(100);        label("Carrying capacity of one FT cell for memory B activation (unitless)")# Table S2 (Yefenof 1986 ref 21)
    BRN      <- fixed(75000);      label("BCR number per B cell (receptors/cell)")                   # Table S2 (Greer 2004 ref 20)
    KRb      <- fixed(1);          label("Occupied BCR number for half-max naive B activation (unitless)")# Table S2 (Bell 1970 ref 18)
    Kab      <- fixed(1e-6);       label("Ag-BCR/Ab binding association constant (1/pM, lumped middle clone, i=1)")# Table S2 note 1 (middle of 17-clone 2-fold ladder, Bell 1970 ref 18)
    NB0      <- fixed(5200);       label("Initial naive B-cell number (cells, i=1 lumped over 17 sub-clones)")# Table S2 (Castiglione 2005 / Crotty 2003 refs 26, 33)

    # 5) ADA (antibody) disposition (paper section 5)
    alphaAb  <- fixed(8.64e8);     label("Ab secretion rate per plasma cell (molecules/(cell*day))") # Table S2 (Auner / Hibi refs 28-29)
    betaAb   <- fixed(0.0301);     label("Free Ab elimination rate (1/day)")                         # Table S2 (Castiglione 2005 ref 26)
    betaCmp  <- fixed(0.0301);     label("Ag-Ab immune complex elimination rate (1/day) (paper notes as Ag-specific; set equal to betaAb here as the neutral default)")# paper 'Ag-specific'; see vignette Errata

    # 6) Antigen (therapeutic protein) PK (paper section 6, 1-compartment reduction of the extended 2-cpt model)
    kelAg    <- fixed(log(2) / 2); label("Antigen plasma elimination rate (1/day; from paper t1/2 = 2 days)")# paper Methods (hypothetical Ag t1/2 = 2 days)
    kaAg     <- fixed(0);          label("Antigen absorption rate from injection site (1/day; 0 for IV dose)")# paper Methods (IV bolus)
    Vp       <- fixed(2.75);       label("Plasma volume (L)")                                        # Table S2

    # No IIV; residual error nominal (deterministic QSP)
    propSd <- fixed(0.001); label("Residual proportional error (nominal; deterministic model)")
  })

  model({
    # --- 1) Maturation signal & DC dynamics (Eqs 1-3) ---
    actFrac <- MS / (KMS * Vp + MS + 1e-12)
    d/dt(MS)  <- -betaMS * MS
    d/dt(iDC) <- betaID * iDC0 - betaID * iDC - deltaID * actFrac * iDC
    d/dt(mDC) <- deltaID * actFrac * iDC - betaMD * mDC

    # --- 6) Antigen (therapeutic protein) PK; needed early so subcellular can reference AgConc ---
    AgConc <- Ag / Vp

    # --- 2) Antigen presentation per mature DC (Eqs 4-13; k=1 lumped MHC-II, j=1 T-epitope) ---
    # State variables (per-DC intensive):
    #   AgE   = antigenic protein in endosome (pmole)
    #   pE    = T-epitope peptide in endosome (pmole)
    #   mhcE  = free MHC-II in endosome (pmole)
    #   pmE   = T-epitope-MHC complex in endosome (pmole)
    #   pmM   = T-epitope-MHC complex on DC membrane (pmole)
    #   cpE   = endogenous competing protein in endosome (pmole)
    #   cptE  = competing peptide in endosome (pmole)
    #   cpmE  = competing peptide-MHC complex in endosome (pmole)
    #   cpmM  = competing peptide-MHC complex on DC membrane (pmole)
    #   mhcM  = free MHC-II on DC membrane (pmole)
    # Endosomal MHC-II synthesis (order-of-magnitude value; paper 'alphaM' not tabulated).
    # See vignette Errata for the rationale.
    Msynth <- betaMHC * 1e-6

    AgEinflux <- alphaAgE * AgConc * Vendo
    cpEinflux <- alphaAgE * (cp0pl / Vp) * Vendo

    d/dt(AgE)  <- AgEinflux - betaAgE * AgE
    d/dt(pE)   <- betaAgE * AgE - betaPep * pE - konMHC * pE * mhcE + koffC * pmE
    d/dt(cpE)  <- cpEinflux - betaAgE * cpE
    d/dt(cptE) <- betaAgE * cpE - betaPep * cptE - konMHC * cptE * mhcE + koffC * cpmE
    d/dt(mhcE) <- Msynth - betaMHC * mhcE - konMHC * pE * mhcE + koffC * pmE - konMHC * cptE * mhcE + koffC * cpmE + kIntCplx * mhcM
    d/dt(pmE)  <- konMHC * pE * mhcE - koffC * pmE - betaPM * pmE - kExt * pmE
    d/dt(cpmE) <- konMHC * cptE * mhcE - koffC * cpmE - betaPM * cpmE - kExt * cpmE
    d/dt(pmM)  <- kExt * pmE - kIntCplx * pmM
    d/dt(cpmM) <- kExt * cpmE - kIntCplx * cpmM
    d/dt(mhcM) <- -kIntCplx * mhcM

    # Total peptide-MHC on all DC membranes, in molecule number (Eq 14.2 with k=1, j=1):
    #   pM_j = pmM [pmole/DC] * 1e-12 [mol/pmole] * N_A [molec/mol] * mDC [DCs]
    #   1e-12 * 6.02e23 = 6.02e11
    pMtot <- pmM * 6.02e11 * mDC

    # --- 3) T helper cell activation & differentiation (Eqs 14-18) ---
    # D_{N,j} = mDC/(mDC + NT + aTn + aTm + MT) * pM/(pM + KpMN)     [naive activation function]
    # D_{M,j} = mDC/(mDC + NT + aTn + aTm + MT) * pM/(pM + KpMM)     [memory activation function]
    # E_{N,j} = mDC/(mDC + NT + aTn + aTm + MT) * (pM - KpMN)/(pM + KpMN)  [-1..+1: proliferate/differentiate]
    # E_{M,j} = mDC/(mDC + NT + aTn + aTm + MT) * (pM - KpMM)/(pM + KpMM)
    poolT <- mDC + NT + aTn + aTm + MT + 1e-12
    Dn <- (mDC / poolT) * (pMtot / (pMtot + KpMN + 1e-12))
    Dm <- (mDC / poolT) * (pMtot / (pMtot + KpMM + 1e-12))
    En <- (mDC / poolT) * ((pMtot - KpMN) / (pMtot + KpMN + 1e-12))
    Em <- (mDC / poolT) * ((pMtot - KpMM) / (pMtot + KpMM + 1e-12))

    # (E+1)/2 in [0,1] is the proliferation fraction; (1-E)/2 differentiates -> MT or FT
    prolTn <- (En + 1) / 2
    prolTm <- (Em + 1) / 2
    diffTn <- 1 - prolTn
    diffTm <- 1 - prolTm

    d/dt(NT)  <- betaNT * NT0 - betaNT * NT - deltaNT * Dn * NT
    d/dt(aTn) <- deltaNT * Dn * NT + rhoAT * prolTn * aTn - betaAT * aTn
    d/dt(aTm) <- deltaMT * Dm * MT + rhoAT * prolTm * aTm - betaAT * aTm
    d/dt(MT)  <- rhoAT * diffTn * fracT1 * aTn + rhoAT * diffTm * fracT1 * aTm - deltaMT * Dm * MT - betaMT * MT
    d/dt(FT)  <- rhoAT * diffTn * (1 - fracT1) * aTn + rhoAT * diffTm * (1 - fracT1) * aTm - betaFT * FT

    # --- 4) B cell activation & differentiation (Eqs 19-24, i=1 lumped ADA clone) ---
    # Free-Ag / BCR / Ab binding (Eqs 4.1.x), i=1:
    #   ro = Kab*[Ag_f] / (1 + Kab*[Ag_f]); approximate [Ag_f] ~ AgConc.
    ro   <- (Kab * AgConc) / (1 + Kab * AgConc + 1e-12)
    Rbcr <- ro * BRN                                    # occupied BCR per B cell (Eq 19.2)
    Ffun <- Rbcr / (Rbcr + KRb + 1e-12)                 # activation function (Eq 19.1)
    Gfun <- Ffun * (1 - ro)                             # tolerance function (Eq 20.1)
    Hfun <- (Rbcr - KRb) / (Rbcr + KRb + 1e-12)         # proliferation switch (Eq 20.2)

    # P_N, P_M: helper-T support (Eqs 19.3, 21.1); i=1 so sums collapse.
    poolB <- NB + aBn + aBm + MB + 1e-12
    Pn <- (CCN * FT) / (CCN * FT + poolB)
    Pm <- (CCM * FT) / (CCM * FT + poolB)

    prolB  <- (Hfun + 1) / 2
    diffB  <- 1 - prolB

    d/dt(NB)  <- -deltaNB * Ffun * Pn * NB
    d/dt(aBn) <- deltaNB * Gfun * Pn * NB + rhoABN * prolB * Pn * aBn - betaAB * aBn
    d/dt(aBm) <- deltaMB * Gfun * Pm * MB + rhoABM * prolB * Pm * aBm - betaAB * aBm
    d/dt(MB)  <- rhoABN * diffB * Pn * fracB1 * aBn + rhoABM * diffB * Pm * fracB1 * aBm - deltaMB * Ffun * Pm * MB - betaMB * MB
    d/dt(SP)  <- rhoABN * diffB * Pn * fracB2 * aBn + rhoABM * diffB * Pm * fracB2 * aBm - betaSP * SP
    d/dt(LP)  <- rhoABN * diffB * Pn * (1 - fracB1 - fracB2) * aBn + rhoABM * diffB * Pm * (1 - fracB1 - fracB2) * aBm - betaLP * LP

    # --- 5) ADA (Ab) disposition (Eq 25); alphaAb in molecules/(cell.day), so /N_A -> pmole/(cell.day) ---
    AbProd <- (alphaAb / 6.02e11) * (SP + LP)
    freeFrac  <- (1 - ro) * (1 - ro)          # (1 - ro)^2
    boundFrac <- 2 * ro - ro * ro             # 2*ro - ro^2
    d/dt(Ab) <- AbProd - betaAb * freeFrac * Ab - betaCmp * boundFrac * Ab

    # --- 6) Antigen PK with immune-mediated clearance (Eq 27 reduced form) ---
    totBcells <- NB + aBn + aBm + MB
    BCRpmole  <- (BRN / 6.02e11) * totBcells
    AgSink    <- 2 * ro * betaCmp * Ab + kIntCplx * ro * BCRpmole
    d/dt(AgIS) <- -kaAg * AgIS
    d/dt(Ag)   <-  kaAg * AgIS - kelAg * Ag - AgSink

    # === Initial conditions (paper Table S2 human) ===
    iDC(0) <- iDC0
    NT(0)  <- NT0
    NB(0)  <- NB0
    # All others start at 0 (Table S2)

    # Observation: Ag concentration in plasma (pM); paper Figure 4d
    Cc <- AgConc
    Cc ~ prop(propSd)
  })
}
