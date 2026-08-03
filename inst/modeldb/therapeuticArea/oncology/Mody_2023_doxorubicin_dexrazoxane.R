Mody_2023_doxorubicin_dexrazoxane <- function() {
  description <- "In vitro (AC16 immortalized human cardiomyocyte cell line) toxicodynamic (TD) model of doxorubicin (DOX) induced cardiotoxicity and dexrazoxane (DEX) cardioprotection, linked to clinical human pharmacokinetics for the in vitro to in vivo translation of drug-drug interaction and Q3W / Q1W dosing-regimen optimization. Clinical DOX PK is a 3-compartment mammillary model with linear elimination and 15-min IV infusion (parameters from Kontny et al. 2013, ref [30], reproduced inline in Mody 2023 Table 2 for a 1.8 m^2 body surface area typical subject). Clinical DEX PK is a 2-compartment mammillary model with linear elimination and 15-min IV infusion, originally fitted to phase 1 data reported by Earhart et al. 1982 Cancer Res 42(12):5255-61 (Mody 2023 Table 2). The TD model on AC16 has (1) an exponential-growth baseline for the untreated cell viability R, (2) a linear DOX growth-inhibition slope S_DOX on kg, (3) a linear DEX growth-inhibition slope S_DEX on kg, (4) a Hill (capacity-limited) DOX stimulation-of-death signal K_DOX = Kmax * CDOX / (KC50 + CDOX) delayed through three transit compartments K1..K3 with rate 1/tau_DOX, and (5) a Hill (capacity-limited) DEX inhibition of the DOX-death signal K_DEX = Imax * CDEX / (IC50 + CDEX). The DEX inhibition enters the transit-chain source term as a subtraction on Kmax (see Assumptions and deviations in the vignette). No IIV or residual error is fitted (Table 1 reports %RSE on point estimates only); an arbitrary 10% CV IIV is added at simulation time per the paper's Methods. Placeholder additive residual SDs are encoded as fixed() so the model parses cleanly."
  reference <- paste(
    "Mody H, Vaidya TR, Ait-Oudhia S (2023).",
    "In vitro to clinical translational pharmacokinetic/pharmacodynamic",
    "modeling of doxorubicin (DOX) and dexrazoxane (DEX) interactions:",
    "Safety assessment and optimization.",
    "Scientific Reports 13:3100.",
    "doi:10.1038/s41598-023-29964-4.",
    sep = " "
  )
  vignette <- "Mody_2023_doxorubicin_dexrazoxane"

  units <- list(
    time          = "h",
    dosing        = "mg (IV; for a 1.8 m^2 body-surface-area typical subject, mg/m^2 doses are pre-multiplied by 1.8 in the event table)",
    concentration = "mg/L (plasma, both drugs); uM (TD driver, derived); % (cell viability, TD readout)"
  )

  covariateData <- list()

  # AC16 % cell viability is a proliferating-cell pool endpoint that does
  # not map to the Simeoni-family `cycling_cells` (tumour cells) or the
  # Friberg-family `circ` (haematopoietic cells) canonicals; declared as
  # paper-specific per compartment-names.md.
  paper_specific_compartments <- c("viability")

  population <- list(
    species         = "human + in vitro (AC16 immortalized human cardiomyocyte cell line)",
    n_subjects      = 500L,
    n_studies       = 1L,
    disease_state   = "In vitro cardiomyocyte cell viability under DOX (0.5-10 uM), DEX (5-100 uM), or DOX (0.5 uM) + DEX (5-100 uM) exposures over 72 h; clinical PK simulations for a body-surface-area of 1.8 m^2 typical adult receiving DOX and DEX by 15-min IV infusion.",
    dose_range      = "DOX 20-60 mg/m^2 Q3W IV; DEX 500 mg/m^2 (10:1 DEX:DOX at 50 mg/m^2 DOX) up to 2500 mg/m^2 Q3W IV; DOX 16.67 mg/m^2 Q1W dose-fractionation regimen.",
    cell_line       = "AC16, immortalized adult human ventricular cardiomyocyte cell line (Millipore Sigma-EMD), derived from primary human ventricular myocytes; retains key transcriptional factors and myogenic markers.",
    culture         = "DMEM + non-essential amino acids + sodium bicarbonate + penicillin/streptomycin + PBS; incubator 37 C. Cells seeded at 10x10^3 cells/well (100 uL) of a 96-well plate; adhesion 12 h; drug exposure 0-72 h; CCK-8 colorimetric viability assay (OD 450 nm) at 0, 12, 24, 48, 72 h.",
    dosing_frequency = "Q3W (once every 3 weeks) preferred regimen; Q1W (once every 1 week) dose-fractionation alternative; both drugs given as 15-min IV infusions simultaneously; three cycles simulated (9 weeks total).",
    notes           = "Simulation cohort of n=500 subjects with an arbitrary 10% CV IIV on TD parameters (Methods: 'introducing an arbitrary inter-individual variability (IIV) of 10% on the TD parameters'). No demographic covariates modelled; the clinical PK parameters are per 1.8 m^2 typical body-surface-area subject. Table 1 (TD parameters, in vitro AC16 data) is the paper's core deliverable; Table 2 (clinical PK parameters) is plumbing to drive the in vivo TD simulation. Third author (Sihem Ait-Oudhia) affiliation is Merck & Co Quantitative Pharmacology and Pharmacometrics (QP2); first two authors (Hardik Mody and Tanaya R. Vaidya) contributed equally."
  )

  ini({
    # ===================================================================
    # DOXORUBICIN CLINICAL PHARMACOKINETICS (3-compartment IV)
    # ===================================================================
    # From Mody 2023 Table 2 (top). All CL/V values are for a subject
    # with body-surface-area 1.8 m^2 (per-BSA units in the table). The
    # underlying PK model is Kontny et al. 2013 J Clin Oncol
    # doi:10.1200/JCO.2012.44.7466 (paper ref [30]), reproduced inline
    # by Mody 2023 -- values are lifted directly from Table 2, no
    # upstream Kontny paper needed on disk.
    lcl  <- log(53.3)  ; label("Doxorubicin clearance CL (L/h per 1.8 m^2)")     # Mody 2023 Table 2 top
    lvc  <- log(17.7)  ; label("Doxorubicin central volume Vc (L per 1.8 m^2)")  # Mody 2023 Table 2 top
    lq   <- log(58.7)  ; label("Doxorubicin Q2 inter-compartmental clearance (L/h per 1.8 m^2)")  # Mody 2023 Table 2 top
    lvp  <- log(1830)  ; label("Doxorubicin first peripheral volume Vp (L per 1.8 m^2)")          # Mody 2023 Table 2 top
    lq2  <- log(21.8)  ; label("Doxorubicin Q3 inter-compartmental clearance (L/h per 1.8 m^2)")  # Mody 2023 Table 2 top
    lvp2 <- log(71.6)  ; label("Doxorubicin second peripheral volume Vp2 (L per 1.8 m^2)")        # Mody 2023 Table 2 top

    # ===================================================================
    # DEXRAZOXANE CLINICAL PHARMACOKINETICS (2-compartment IV)
    # ===================================================================
    # Table 2 bottom reports micro-constants: kel = 1 (h^-1), k12 = 1
    # (h^-1), k21 = 1 (h^-1), V = 14.6 L (single 4% RSE for kel, 3% for
    # k12, 8% for k21, 5% for V). Under the canonical CL/Vc/Q/Vp macro
    # parameterisation these micro-constants map directly:
    #   CL_dex  = kel * V   = 1  * 14.6 = 14.6 L/h
    #   Vc_dex  = V         =              14.6 L
    #   Q_dex   = k12 * V   = 1  * 14.6 = 14.6 L/h
    #   Vp_dex  = V*k12/k21 = 14.6*1/1  = 14.6 L
    # Original DEX PK was fitted by Mody 2023 (Supp. Fig. 3) to human
    # phase 1 data reported by Earhart RH, Tutsch KD, Koeller JM et al.
    # Cancer Res 1982 42(12):5255-61 (paper ref [29]).
    lcl_dex <- log(14.6) ; label("Dexrazoxane clearance CL_dex (L/h)")               # Mody 2023 Table 2 bottom (= kel * V)
    lvc_dex <- log(14.6) ; label("Dexrazoxane central volume Vc_dex (L)")            # Mody 2023 Table 2 bottom
    lq_dex  <- log(14.6) ; label("Dexrazoxane inter-compartmental clearance Q_dex (L/h)")  # Mody 2023 Table 2 bottom (= k12 * V)
    lvp_dex <- log(14.6) ; label("Dexrazoxane peripheral volume Vp_dex (L)")         # Mody 2023 Table 2 bottom (= V*k12/k21)

    # ===================================================================
    # AC16 CELLULAR-LEVEL TOXICODYNAMIC (TD) PARAMETERS
    # ===================================================================
    # From Mody 2023 Table 1. Simultaneous fit to single-agent DOX,
    # single-agent DEX, and DOX+DEX combination cell-viability data.
    lkg  <- log(0.0115) ; label("AC16 first-order growth rate constant kg (1/h)")                   # Mody 2023 Table 1
    lr0  <- log(101)    ; label("AC16 baseline cell viability R0 (%)")                              # Mody 2023 Table 1
    lsdex <- log(0.00968) ; label("Slope of DEX linear growth-inhibition on kg, S_DEX (1/uM)")      # Mody 2023 Table 1
    lsdox <- log(0.167)   ; label("Slope of DOX linear growth-inhibition on kg, S_DOX (1/uM)")     # Mody 2023 Table 1
    lemax <- log(0.0697)  ; label("Maximal DOX-induced stimulation-of-death rate Kmax,DOX (1/h)")   # Mody 2023 Table 1
    lec50 <- log(0.107)   ; label("DOX concentration for half-maximal stimulation-of-death KC50,DOX (uM)")  # Mody 2023 Table 1
    lktr  <- log(0.126)   ; label("Transit-chain rate constant 1/tau_DOX for DOX stimulation-of-death delay (1/h)")  # Mody 2023 Table 1
    limax_dex <- log(0.0625) ; label("Maximal DEX inhibition on DOX stimulation-of-death Imax,DEXi (1/h)")   # Mody 2023 Table 1
    lic50_dex <- log(39)     ; label("DEX concentration for half-maximal inhibition of DOX-death IC50,DEXi (uM)")  # Mody 2023 Table 1

    # ===================================================================
    # PLASMA-TO-uM UNIT CONVERSION CONSTANTS (FIXED)
    # ===================================================================
    # The TD model was fit with drug concentrations in uM (paper Methods:
    # in vitro exposures 0.5-10 uM DOX and 5-100 uM DEX). Clinical PK is
    # in mg/L (Table 2's CL and V are in L/h and L). Convert the central
    # compartment concentration mg/L -> uM at the PK/TD interface using
    # the molar mass (g/mol) of the free base:
    #   c_uM = Cc (mg/L) * 1000 / MW (g/mol)
    # DOX free base (PubChem CID 31703) MW = 543.52 g/mol -> factor 1.8397 uM per mg/L.
    # DEX free base (PubChem CID  71384) MW = 268.27 g/mol -> factor 3.7275 uM per mg/L.
    # Held FIXED at literature molar mass. Downstream users dosing the
    # HCl salt (DOX HCl MW 579.98) should scale the amt column by
    # 543.52 / 579.98 = 0.9371 at the event-table build time to preserve
    # free-base amounts.
    lconv_dox <- fixed(log(1.8397))  ; label("DOX mg/L -> uM conversion factor (dimensionless; literature)")   # PubChem CID 31703 (DOX free base MW 543.52 g/mol)
    lconv_dex <- fixed(log(3.7275))  ; label("DEX mg/L -> uM conversion factor (dimensionless; literature)")   # PubChem CID 71384 (DEX free base MW 268.27 g/mol)

    # ===================================================================
    # RESIDUAL ERROR (PLACEHOLDERS)
    # ===================================================================
    # Mody 2023 Tables 1 and 2 report %RSE on point estimates only; no
    # residual SD is tabulated for DOX plasma, DEX plasma, or AC16 % cell
    # viability. Small operator-chosen placeholders so the multi-output
    # model parses cleanly; the paper's Methods introduces a 10% CV IIV
    # on TD parameters at simulation time (arbitrary, not fit).
    # operator-chosen placeholder (paper reports no residual SD)
    addSd            <- fixed(1)   ; label("Additive residual SD on DOX plasma Cc (mg/L; placeholder)")
    # operator-chosen placeholder (paper reports no residual SD)
    addSd_dex        <- fixed(1)   ; label("Additive residual SD on DEX plasma Cc_dex (mg/L; placeholder)")
    # operator-chosen placeholder (paper reports no residual SD)
    addSd_viability  <- fixed(5)   ; label("Additive residual SD on AC16 cell viability (%; placeholder)")
  })

  model({
    # ===================================================================
    # 1. Typical-value parameter back-transforms (no etas)
    # ===================================================================
    cl   <- exp(lcl)
    vc   <- exp(lvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    q2   <- exp(lq2)
    vp2  <- exp(lvp2)
    cl_dex <- exp(lcl_dex)
    vc_dex <- exp(lvc_dex)
    q_dex  <- exp(lq_dex)
    vp_dex <- exp(lvp_dex)
    kg   <- exp(lkg)
    r0   <- exp(lr0)
    sdex <- exp(lsdex)
    sdox <- exp(lsdox)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    ktr  <- exp(lktr)
    imax_dex <- exp(limax_dex)
    ic50_dex <- exp(lic50_dex)
    conv_dox <- exp(lconv_dox)
    conv_dex <- exp(lconv_dex)

    # ===================================================================
    # 2. PK micro-constants
    # ===================================================================
    kel     <- cl / vc
    k12     <- q  / vc
    k21     <- q  / vp
    k13     <- q2 / vc
    k31     <- q2 / vp2
    kel_de  <- cl_dex / vc_dex
    k12_de  <- q_dex  / vc_dex
    k21_de  <- q_dex  / vp_dex

    # ===================================================================
    # 3. Clinical PK ODEs
    # ===================================================================
    # DOX 3-compartment mammillary; dose enters central via 15-min IV
    # infusion (event table sets rate).
    d/dt(central)     <- -kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # DEX 2-compartment mammillary; dose enters central_dex via 15-min
    # IV infusion.
    d/dt(central_dex)     <- -kel_de * central_dex -
                              k12_de * central_dex + k21_de * peripheral1_dex
    d/dt(peripheral1_dex) <-  k12_de * central_dex - k21_de * peripheral1_dex

    # ===================================================================
    # 4. PK -> TD interface (mg/L plasma -> uM TD driver)
    # ===================================================================
    Cc     <- central     / vc
    Cc_dex <- central_dex / vc_dex
    c_dox_uM <- Cc     * conv_dox
    c_dex_uM <- Cc_dex * conv_dex

    # ===================================================================
    # 5. Cellular-level toxicodynamic (TD) ODEs
    # ===================================================================
    # Mody 2023 Eq 6a-6f (combination DOX + DEX on AC16):
    #   K_source_stimdeath = (Kmax,DOX - K_DEX) * CDOX / (KC50,DOX + CDOX)
    #   K_DEX              = Imax,DEXi * CDEX / (IC50,DEXi + CDEX)
    # delayed through three transit compartments K1DOX -> K3DOX with
    # first-order transit rate 1/tau_DOX (= ktr in this file), then the
    # terminal transit3 drives -R * transit3 death of AC16 cell
    # viability. Growth uses two multiplicative linear inhibitions
    # (Eqs 4e / 5): (1 - S_DOX * CDOX) * (1 - S_DEX * CDEX). The
    # subtractive-on-Kmax form of the DEX -> DOX-death inhibition is
    # my interpretation of the printed Table 1 units (Kmax and Imax both
    # in 1/h) and text ("capacity limited Hill function K_DEX on DOX
    # cell-death stimulation") given the paper's PDF equations were not
    # OCR-decoded -- see the vignette Assumptions and deviations
    # section for a note on this interpretation.
    k_dex_inhib   <- imax_dex * c_dex_uM / (ic50_dex + c_dex_uM)
    k_dox_stim    <- (emax - k_dex_inhib) * c_dox_uM / (ec50 + c_dox_uM)

    d/dt(transit1) <- ktr * (k_dox_stim - transit1)
    d/dt(transit2) <- ktr * (transit1   - transit2)
    d/dt(transit3) <- ktr * (transit2   - transit3)

    d/dt(viability) <- kg * viability * (1 - sdox * c_dox_uM) * (1 - sdex * c_dex_uM) -
                       transit3 * viability

    # AC16 cell viability starts at baseline R0 = 101%. Transit chain
    # starts at 0 (no delayed death signal at t = 0); PK compartments
    # start at 0 (dose amounts arrive via the event table).
    viability(0) <- r0

    # ===================================================================
    # 6. Observation variables and error models
    # ===================================================================
    Cc         ~ add(addSd)
    Cc_dex     ~ add(addSd_dex)
    viability  ~ add(addSd_viability)
  })
}
