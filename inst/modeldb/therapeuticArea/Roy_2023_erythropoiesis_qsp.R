Roy_2023_erythropoiesis_qsp <- function() {
  description <- paste(
    "QSP. Quantitative systems pharmacology model of human",
    "erythropoiesis and its disruption in anemia due to chronic",
    "kidney disease (CKD). 20 ODE states covering the",
    "Hb-PHD-HIF-EPO negative-feedback loop (Hb-dependent prolyl",
    "hydroxylase activity, PHD-driven HIF-alpha degradation,",
    "HIF-driven renal EPO production), a two-compartment plasma /",
    "peripheral distribution for endogenous EPO and for each",
    "erythropoiesis-stimulating agent, target-mediated disposition",
    "via reversible EPO-EPOR binding on erythroid progenitors, and a",
    "progenitor -> precursor -> reticulocyte -> RBC maturation chain",
    "in which EPO-EPOR complex rescues progenitors from apoptosis",
    "through a 4th-order Hill function. Four therapies are built in:",
    "recombinant human EPO (epoetin alfa, SC), darbepoetin alfa (SC),",
    "and the oral HIF prolyl-hydroxylase inhibitors vadadustat and",
    "daprodustat, the latter with saturable non-linear clearance.",
    "Hemoglobin is an algebraic function of reticulocyte and RBC",
    "counts and their mean corpuscular hemoglobin. The published",
    "model carries five alternate parameterizations (virtual",
    "patients: healthy and CKD stages 1.5, 3, 4, 5) that differ in",
    "five parameters; the values shipped in ini() are the",
    "non-dialysis reference virtual patient (CKD stage 4), and the",
    "full stage table is reproduced in the validation vignette.",
    "All values are hand-calibrated literature / fitted constants",
    "from the authors' SimBiology release; the paper reports no IIV",
    "and no residual-error model (population variability is generated",
    "by log-normal resampling of parameters, not by an omega matrix).",
    sep = " "
  )
  reference <- paste(
    "Roy M, Saroha S, Sarma U, Sarathy H, Kumar R.",
    "Quantitative systems pharmacology model of erythropoiesis to",
    "simulate therapies targeting anemia due to chronic kidney",
    "disease. Front Pharmacol. 2023;14:1274490.",
    "doi:10.3389/fphar.2023.1274490.",
    "Equations 1-28 are from the main text section 3.3; the expanded",
    "PHD and progenitor-apoptosis Hill forms and the EPO molecule",
    "unit conversion are from the Supplementary Appendix",
    "('Detailed Equations'). Every parameter value, unit and",
    "per-CKD-stage alternate parameterization is from the",
    "Supplementary 'Parameter Dashboard' (Table1.XLSX, sheets",
    "'Parameters' and 'Species'). Reaction stoichiometry, rate laws",
    "and assignment rules were read from the authors' own SimBiology",
    "project released as Supplementary DataSheet1.ZIP, file",
    "'Model File and Scripts/ErythropoiesisModel_5Aug2023_Final3.sbproj'.",
    sep = " "
  )
  vignette <- "Roy_2023_erythropoiesis_qsp"

  # Mechanistic states of the published SimBiology model. Names follow the
  # authors' species names (hyphens replaced by underscores so they are
  # valid R symbols) so the source trace is unambiguous.
  paper_specific_compartments <- c(
    "HIFa", "EPO_plasma", "EPO_periphery", "EPO_receptor",
    "EPO_LR_complex", "rhuEPO_LR_complex", "darbe_LR_complex",
    "rhuEPO_SC_Dose", "rhuEPO_plasma", "rhuEPO_periphery",
    "darbe_SC_Dose", "Darbe_plasma", "darbe_periphery",
    "PHI_Dose_Gut", "PHI_Plasma", "PHI_periphery",
    "Progenitors", "Precursors", "Retics_plasma", "RBCM"
  )

  # Declared explicitly: buildModelDb() infers the registry's dosing column
  # only from compartments literally named "depot" or "central", neither of
  # which exists here. rHuEPO and darbepoetin are dosed subcutaneously (ng)
  # and the oral PHIs into the gut (ug).
  dosing <- c("rhuEPO_SC_Dose", "darbe_SC_Dose", "PHI_Dose_Gut")

  units <- list(
    time          = "h",
    dosing        = "ng (rHuEPO / darbepoetin SC) or ug (oral PHI)",
    concentration = "ng/mL"
  )

  # Units are those declared in the Supplementary Parameter Dashboard sheet
  # "Species" (column InitialAmountUnits), so verified = TRUE throughout. The
  # one exception is PHI_periphery, which that sheet declares in ug/mL while
  # PHI_Plasma and every other peripheral species is ng/mL; this
  # implementation reads it as ng/mL (the self-consistent interpretation, and
  # the one that reproduces the paper's vadadustat response), so it is flagged
  # unverified. See the vignette Errata.
  # Progenitors and Precursors physically reside in bone marrow but are
  # modelled in the central compartment, a caveat the Supplement states
  # explicitly; "blood cell" is the nearest specimen in the vocabulary.
  compartmentData <- list(
    HIFa              = list(analyte = "hypoxia-inducible factor 1-alpha", units = "ng/mL", specimen = "tissue", verified = TRUE),
    EPO_plasma        = list(analyte = "endogenous erythropoietin", units = "ng/mL", specimen = "plasma", verified = TRUE),
    EPO_periphery     = list(analyte = "endogenous erythropoietin", units = "ng/mL", specimen = "tissue", verified = TRUE),
    EPO_receptor      = list(analyte = "free erythropoietin receptor", units = "molecule/L", specimen = "blood cell", verified = TRUE),
    EPO_LR_complex    = list(analyte = "endogenous EPO-EPOR complex", units = "molecule/L", specimen = "blood cell", verified = TRUE),
    rhuEPO_LR_complex = list(analyte = "epoetin alfa-EPOR complex", units = "molecule/L", specimen = "blood cell", verified = TRUE),
    darbe_LR_complex  = list(analyte = "darbepoetin alfa-EPOR complex", units = "molecule/L", specimen = "blood cell", verified = TRUE),
    rhuEPO_SC_Dose    = list(analyte = "epoetin alfa", units = "ng", specimen = "administration site", verified = TRUE),
    rhuEPO_plasma     = list(analyte = "epoetin alfa", units = "ng/mL", specimen = "plasma", verified = TRUE),
    rhuEPO_periphery  = list(analyte = "epoetin alfa", units = "ng/mL", specimen = "tissue", verified = TRUE),
    darbe_SC_Dose     = list(analyte = "darbepoetin alfa", units = "ng", specimen = "administration site", verified = TRUE),
    Darbe_plasma      = list(analyte = "darbepoetin alfa", units = "ng/mL", specimen = "plasma", verified = TRUE),
    darbe_periphery   = list(analyte = "darbepoetin alfa", units = "ng/mL", specimen = "tissue", verified = TRUE),
    PHI_Dose_Gut      = list(analyte = "HIF prolyl-hydroxylase inhibitor (vadadustat or daprodustat)", units = "ug", specimen = "administration site", verified = TRUE),
    PHI_Plasma        = list(analyte = "HIF prolyl-hydroxylase inhibitor (vadadustat or daprodustat)", units = "ng/mL", specimen = "plasma", verified = TRUE),
    PHI_periphery     = list(analyte = "HIF prolyl-hydroxylase inhibitor (vadadustat or daprodustat)", units = "ng/mL", specimen = "tissue", verified = FALSE),
    Progenitors       = list(analyte = "erythroid progenitor cells (CFU-E)", units = "molecule/L", specimen = "blood cell", verified = TRUE),
    Precursors        = list(analyte = "erythroid precursor cells", units = "molecule/L", specimen = "blood cell", verified = TRUE),
    Retics_plasma     = list(analyte = "reticulocytes", units = "molecule/L", specimen = "blood cell", verified = TRUE),
    RBCM              = list(analyte = "mature erythrocytes", units = "molecule/L", specimen = "blood cell", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species       = "human",
    n_subjects    = 300,
    n_studies     = 10,
    age_range     = "not reported for the virtual population",
    disease_state = paste(
      "Anemia due to chronic kidney disease. Five reference virtual",
      "patients span healthy and CKD stages 1.5, 3, 4 and 5. Virtual",
      "populations are of two kinds: non-dialysis (ND; CKD stages 3,",
      "4 and 5, erythropoiesis-stimulating-agent naive, built around",
      "a CKD 4 reference patient) and hemodialysis (HD; CKD stage 5",
      "only, previously treated with rHuEPO, built around a CKD 5",
      "reference patient). Iron sufficiency is assumed throughout."
    ),
    dose_range    = paste(
      "rHuEPO (epoetin alfa) SC 10,000 IU QW, 20,000 IU Q2W,",
      "50 IU/kg TIW, 10,500 IU TIW and 40,000 IU QW;",
      "darbepoetin alfa SC 10-100 ug Q2W;",
      "vadadustat PO 150-600 mg QD; daprodustat PO 1-12 mg QD."
    ),
    notes         = paste(
      "This is a virtual population, not an observed cohort. Each",
      "virtual population is 300 patients filtered from 10,000",
      "plausible patients generated by log-normal resampling around a",
      "reference virtual patient (Supplementary, 'Developing Virtual",
      "Population'). Calibration data are the ten clinical trials in",
      "Tables 1 and 2; steady-state Hb, EPO, reticulocyte and RBC",
      "ranges by CKD stage are from Li et al. 2019 and Sheth & Shah",
      "2016 (Figure 3)."
    )
  )

  ini({
    # =====================================================================
    # All values are hand-calibrated literature / fitted constants taken
    # from the Supplementary Parameter Dashboard (Table1.XLSX, sheet
    # "Parameters"), NOT maximum-likelihood estimates. They are exposed as
    # estimable typical values (no fixed() wrapper, following
    # Rao_2023_covid19_qsp / Lu_2023_sglt_qsp) so downstream users can
    # recalibrate. Source tag "PD" = Parameter Dashboard row.
    #
    # Model time unit is hours. Protein concentrations are ng/mL; cellular
    # species are molecule/L; PHI plasma is ng/mL.
    #
    # DEFAULT PARAMETERIZATION = the non-dialysis (ND) reference virtual
    # patient, i.e. CKD stage 4. The five stage-varying parameters are
    # flagged "[VP-varying]"; the full five-stage table is in the vignette.
    # =====================================================================

    # ---- Hb -> PHD -> HIF-alpha limb (Eq. 1-3; Suppl. "Detailed Equations")
    PHD_total <- 52.5    ; label("Total basal PHD concentration (ng/mL)")                                  # PD PHD_total
    k_modulate_HGB_efffect_on_PHD <- 100 ; label("Vmax of the Hb effect on PHD (unitless)")                # PD k_modulate_HGB_efffect_on_PHD (author spelling retained)
    Km_HGB_HIF <- 5      ; label("Km of the Hb effect on PHD (unitless)")                                  # PD Km_HGB_HIF
    n_HGB_PHD  <- 3.3    ; label("Hill coefficient of the Hb effect on PHD (unitless)")                    # PD n_HGB_PHD
    kprod_HIF_basal <- 1.35 ; label("Zero-order HIF-alpha production rate (ng/mL/h)")                      # PD kprod_HIF_basal
    k_modulate_PHD_effect_on_HIF <- 0.34 ; label("Second-order rate of PHD-driven HIF-alpha degradation (mL/ng/h)") # PD k_modulate_PHD_effect_on_HIF

    # ---- HIF-alpha -> EPO production and EPO disposition (Eq. 4-6, 11) ----
    k_production_EPO <- 1.6   ; label("[VP-varying] Vmax of HIF-driven EPO production (ng/mL/h)")          # PD k_production_EPO, CKD 4 column
    Km_prod_EPO      <- 2.23  ; label("Km of HIF-driven EPO production (unitless)")                        # PD Km_prod_EPO
    n2               <- 1     ; label("Hill coefficient of HIF-driven EPO production (unitless)")          # PD n2
    kel_EPO_plasma   <- 0.02  ; label("Non-specific plasma clearance rate of endogenous EPO (1/h)")        # PD kel_EPO_plasma
    k_epo_cp <- 0.1    ; label("Endogenous EPO central-to-peripheral rate constant (1/h)")                 # PD k_epo_cp
    k_epo_pc <- 0.2762 ; label("Endogenous EPO peripheral-to-central rate constant (1/h)")                 # PD k_epo_pc

    # ---- EPO / ESA binding to EPOR (Eq. 7-9, 12-13) -----------------------
    kon_EPO_LR_complex  <- 1e-13 ; label("EPO-EPOR association rate constant (L/molecule/h)")              # PD kon_EPO_LR_complex
    koff_EPO_LR_complex <- 5.5   ; label("rHuEPO-EPOR dissociation rate constant (1/h)")                   # PD koff_EPO_LR_complex
    kd <- 5.5e13 ; label("EPO-EPOR equilibrium dissociation constant (molecule/L)")                        # PD kd; endogenous off-rate is encoded as kon*kd in the SimBiology rate law
    kdeg_EPO_LR <- 1.2 ; label("Internalization rate of the EPO-EPOR and rHuEPO-EPOR complexes (1/h)")      # PD kdeg_EPO_LR
    MW_EPO <- 5.06e-11 ; label("Molecular weight of endogenous EPO and rHuEPO (ng/molecule)")              # PD MW_EPO
    MW_ESA <- 6.16e-11 ; label("Molecular weight of darbepoetin alfa (ng/molecule)")                       # PD MW_ESA

    # ---- rHuEPO (epoetin alfa) PK (section 3.3.5) -------------------------
    ka1_ESA_SCtoPlasma <- 0.0362 ; label("rHuEPO subcutaneous absorption rate constant (1/h)")             # PD ka1_ESA_SCtoPlasma
    kel_ESA_plasma <- 0.02 ; label("rHuEPO non-specific plasma elimination rate constant (1/h)")           # PD kel_ESA_plasma
    k_esa_cp <- 0.038 ; label("rHuEPO central-to-peripheral rate constant (1/h)")                          # PD k_esa_cp
    k_esa_pc <- 0.01  ; label("rHuEPO peripheral-to-central rate constant (1/h)")                          # PD k_esa_pc

    # ---- Darbepoetin alfa PK/binding (section 3.3.6) ----------------------
    ka1_darbe_SCtoPlasma <- 0.03 ; label("Darbepoetin subcutaneous absorption rate constant (1/h)")        # PD ka1_darbe_SCtoPlasma
    kel_darbe_plasma <- 0.018 ; label("Darbepoetin plasma elimination rate constant (1/h)")                # PD kel_darbe_plasma
    k_darbe_cp <- 0.03  ; label("Darbepoetin central-to-peripheral rate constant (1/h)")                   # PD k_darbe_cp
    k_darbe_pc <- 0.006 ; label("Darbepoetin peripheral-to-central rate constant (1/h)")                   # PD k_darbe_pc
    kon_darbe_LR_complex  <- 1e-14 ; label("Darbepoetin-EPOR association rate constant (L/molecule/h)")    # PD kon_darbe_LR_complex
    koff_darbe_LR_complex <- 6.3   ; label("Darbepoetin-EPOR dissociation rate constant (1/h)")            # PD koff_darbe_LR_complex
    kint_ESA_LR <- 0.65 ; label("Internalization rate of the darbepoetin-EPOR complex (1/h)")              # PD kint_ESA_LR

    # ---- HIF prolyl-hydroxylase inhibitor (PHI) PK and effect -------------
    # Defaults are VADADUSTAT. Daprodustat overrides are given in the
    # vignette (Vmax_DRUG_PHD 2.2, ka_gutToPlasma 0.2239, kel_PHI 0.025,
    # CLD_PHI 0.173, Vmax_nonlinearClearance 115, Vc 5.38).
    Vmax_DRUG_PHD <- 1.27 ; label("Vmax of the PHI effect on PHD (unitless; vadadustat)")                  # PD Vmax_DRUG_PHD, Vada row
    Km_PHI_drug   <- 1e-4 ; label("Km of the PHI effect on PHD (unitless, drug conc in ug/mL)")            # PD Km_PHI_drug
    n_PHIeffPHD   <- 0.17 ; label("Hill coefficient of the PHI effect on PHD (unitless)")                  # PD n_PHIeffPHD
    ka_gutToPlasma <- 0.36 ; label("PHI gut-to-plasma absorption rate constant (1/h; vadadustat)")         # PD ka_gutToPlasma, Vada row
    kel_PHI <- 0.1  ; label("PHI linear plasma elimination rate constant (1/h; vadadustat)")               # PD kel_VADA, Vada row
    CLD_PHI <- 900  ; label("PHI central-peripheral distribution clearance (mL/h; vadadustat)")            # PD CLD_VADA, Vada row
    Vmax_nonlinearClearance <- 0 ; label("Vmax of saturable PHI clearance (ng/mL/h; 0 for vadadustat)")     # PD Vmax_non-linearClearanceDapro, Vada row
    Km_nonlinearClearance <- 4.52 ; label("Km of saturable PHI clearance (ng/mL)")                          # PD Km_non-linearClearanceDapro

    # ---- Erythroid maturation chain (Eq. 14-27) ---------------------------
    kprod_prog <- 1.5e9 ; label("[VP-varying] Zero-order progenitor production rate (molecule/L/h)")        # PD kprod_prog, CKD 4 column
    k_deg_prog <- 0.026 ; label("[VP-varying] Baseline progenitor apoptosis rate constant (1/h)")           # PD k_deg_prog, CKD 4 column
    k_progenitorsToPrecursors <- 0.0046 ; label("Progenitor-to-precursor differentiation rate constant (1/h)") # PD k_progenitorsToPrecursors
    progenitorAmplification <- 32 ; label("Precursors generated per differentiating progenitor (unitless)") # sbproj reaction_29 hard-coded factor 32
    Vmax_EPO_EPOR_effective <- 0.99 ; label("Maximal fractional rescue of progenitors by EPO-EPOR (unitless)") # PD Vmax_EPO-EPOR_effective
    Km_EPO_EPOR_Progenitors <- 2.5e12 ; label("[VP-varying] Km of the EPO-EPOR rescue of progenitors (molecule/L)") # PD Km_EPO-EPOR_Progenitors, CKD 4 column
    nEPO_EPOR_to_Progenitors <- 4 ; label("Hill coefficient of the EPO-EPOR rescue of progenitors (unitless)") # PD nEPO-EPOR_to_Progenitors
    kdeg_Precursors <- 0.05 ; label("Precursor loss rate constant (1/h)")                                  # PD kdeg_Precursors
    k_PrecursorsToPretics <- 0.04 ; label("Precursor-to-reticulocyte release rate constant (1/h)")          # PD k_PrecursorsToPretics
    kdeg_Retics_Plasma <- 0.03 ; label("Plasma reticulocyte loss rate constant (1/h)")                     # PD kdeg_Retics_Plasma
    k_ReticsToRBCM <- 0.027 ; label("Reticulocyte-to-erythrocyte maturation rate constant (1/h)")          # PD k_ReticsToRBCM
    kdeg_RBCM <- 0.000473 ; label("[VP-varying] Erythrocyte loss rate constant (1/h)")                     # PD kdeg_RBCM, CKD 4 column

    # ---- Hemoglobin and reporting scalars (Eq. 28; rules 1-2) -------------
    MCH_Reti <- 3e-11   ; label("Mean corpuscular hemoglobin of a reticulocyte (g/molecule)")              # PD MCH_Reti
    MCH_RBCM <- 2.9e-11 ; label("Mean corpuscular hemoglobin of an erythrocyte (g/molecule)")              # PD MCH_RBCM
    scaler7    <- 0.1    ; label("Unit scalar giving Hb in g/dL (L/dL)")                                   # PD scaler7
    scaler_EPO <- 119.05 ; label("EPO mass-to-activity conversion (milliU/ng)")                            # PD scaler_EPO

    # ---- Volumes ----------------------------------------------------------
    Vc <- 4.9573 ; label("Central compartment volume (L)")                                                 # PD Central_compart_0, Default row
    Vp <- 3.067  ; label("Peripheral compartment volume (L)")                                              # PD Peripheral_0
  })

  model({
    # ---------------------------------------------------------------------
    # Unit reconciliation. The published SimBiology project runs with unit
    # conversion enabled, so two conversions that are implicit in the
    # released rate laws must be written out here:
    #   (a) EPO_plasma/MW_EPO evaluates to molecule/mL but is multiplied by
    #       EPO_receptor in molecule/L, so it is scaled by 1000 (mL per L).
    #       The reverse conversion appears on the plasma side of the same
    #       reaction, where the molecule/L/h flux times MW gives ng/L/h and
    #       must be divided by 1000 to give ng/mL/h.
    #   (b) scalar_PHI_DRUG carries units mL/ug, i.e. the PHI Hill function
    #       takes plasma PHI in ug/mL while the PHI_Plasma state is ng/mL.
    # Both are confirmed by the validation in the vignette: with these
    # factors the model reproduces the paper's reported reference-patient
    # Hb responses, and without them the EPO-EPOR limb is inert.
    # ---------------------------------------------------------------------
    mLperL <- 1000
    ngPerUg <- 1000

    # ---- Algebraic outputs (sbproj rules 1, 2, 5, 14) --------------------
    HGB <- (MCH_Reti * Retics_plasma + MCH_RBCM * RBCM) * scaler7          # Eq. 28 / rule_1
    RBC_total <- Retics_plasma + RBCM                                      # rule_5
    Retic_pct <- 100 * Retics_plasma / RBC_total                           # rule_3
    EPO_plasma_IU <- (EPO_plasma + rhuEPO_plasma + Darbe_plasma) * scaler_EPO # rule_2

    # PHI concentration entering the Hill function, in ug/mL. The 1e-24
    # offset is a solver guard only: PHI_Plasma^0.17 has an infinite
    # derivative at zero, which stalls the integrator in the undosed
    # (steady-state) run. It is numerically negligible at any dosed level.
    PHIugmL <- PHI_Plasma / ngPerUg + 1e-24
    PHIeffect <- Vmax_DRUG_PHD * PHIugmL^n_PHIeffPHD /
      (Km_PHI_drug^n_PHIeffPHD + PHIugmL^n_PHIeffPHD)
    HGBeffect <- k_modulate_HGB_efffect_on_PHD / (Km_HGB_HIF + HGB^n_HGB_PHD)
    PHD <- PHD_total / (1 + PHIeffect + HGBeffect)                         # Eq. 1 / rule_14

    # ---- Receptor binding fluxes, molecule/L/h (Eq. 7-9; reactions 8_1, 2_1, 2_2)
    bindEPO <- kon_EPO_LR_complex * EPO_receptor *
      (EPO_plasma / MW_EPO) * mLperL -
      kon_EPO_LR_complex * kd * EPO_LR_complex
    bindRHU <- kon_EPO_LR_complex * EPO_receptor *
      (rhuEPO_plasma / MW_EPO) * mLperL -
      koff_EPO_LR_complex * rhuEPO_LR_complex
    bindDAR <- kon_darbe_LR_complex * EPO_receptor *
      (Darbe_plasma / MW_ESA) * mLperL -
      koff_darbe_LR_complex * darbe_LR_complex

    # ---- Fractional rescue of progenitors from apoptosis ------------------
    # All three ligand-receptor complexes act on the same pool
    # (sbproj reaction_7); Suppl. "Detailed Equations" gives the 4th-order form.
    LRtotal <- EPO_LR_complex + rhuEPO_LR_complex + darbe_LR_complex
    fEPO <- Vmax_EPO_EPOR_effective * LRtotal^nEPO_EPOR_to_Progenitors /
      (Km_EPO_EPOR_Progenitors^nEPO_EPOR_to_Progenitors +
         LRtotal^nEPO_EPOR_to_Progenitors)

    # ---- Hb-PHD-HIF loop (Eq. 2-3) ----------------------------------------
    d/dt(HIFa) <- kprod_HIF_basal -
      HIFa * PHD * k_modulate_PHD_effect_on_HIF

    # ---- Endogenous EPO (Eq. 4-6, 10-11) ----------------------------------
    d/dt(EPO_plasma) <-
      k_production_EPO * HIFa^n2 / (Km_prod_EPO^n2 + HIFa^n2) -
      kel_EPO_plasma * EPO_plasma -
      (k_epo_cp * EPO_plasma - k_epo_pc * EPO_periphery) -
      bindEPO * MW_EPO / mLperL
    d/dt(EPO_periphery) <- k_epo_cp * EPO_plasma - k_epo_pc * EPO_periphery

    # ---- EPOR and ligand-receptor complexes (Eq. 12-13) -------------------
    d/dt(EPO_LR_complex)    <- bindEPO - kdeg_EPO_LR * EPO_LR_complex
    d/dt(rhuEPO_LR_complex) <- bindRHU - kdeg_EPO_LR * rhuEPO_LR_complex
    d/dt(darbe_LR_complex)  <- bindDAR - kint_ESA_LR * darbe_LR_complex
    d/dt(EPO_receptor) <- -bindEPO - bindRHU - bindDAR +
      kdeg_EPO_LR * EPO_LR_complex + kdeg_EPO_LR * rhuEPO_LR_complex +
      kint_ESA_LR * darbe_LR_complex                                       # reaction_28

    # ---- rHuEPO (epoetin alfa) PK -----------------------------------------
    d/dt(rhuEPO_SC_Dose) <- -ka1_ESA_SCtoPlasma * rhuEPO_SC_Dose
    d/dt(rhuEPO_plasma) <-
      ka1_ESA_SCtoPlasma * rhuEPO_SC_Dose / (Vc * mLperL) -
      kel_ESA_plasma * rhuEPO_plasma -
      (k_esa_cp * rhuEPO_plasma - k_esa_pc * rhuEPO_periphery) -
      bindRHU * MW_EPO / mLperL
    d/dt(rhuEPO_periphery) <- k_esa_cp * rhuEPO_plasma -
      k_esa_pc * rhuEPO_periphery

    # ---- Darbepoetin alfa PK ----------------------------------------------
    d/dt(darbe_SC_Dose) <- -ka1_darbe_SCtoPlasma * darbe_SC_Dose
    d/dt(Darbe_plasma) <-
      ka1_darbe_SCtoPlasma * darbe_SC_Dose / (Vc * mLperL) -
      kel_darbe_plasma * Darbe_plasma -
      (k_darbe_cp * Darbe_plasma - k_darbe_pc * darbe_periphery) -
      bindDAR * MW_ESA / mLperL
    d/dt(darbe_periphery) <- k_darbe_cp * Darbe_plasma -
      k_darbe_pc * darbe_periphery

    # ---- Oral PHI PK (reactions 17, 20, 21, 9) ----------------------------
    # PHI_Dose_Gut is an amount in ug; dividing by Vc in L gives ug/L = ng/mL.
    d/dt(PHI_Dose_Gut) <- -ka_gutToPlasma * PHI_Dose_Gut
    d/dt(PHI_Plasma) <- ka_gutToPlasma * PHI_Dose_Gut / Vc -
      kel_PHI * PHI_Plasma -
      CLD_PHI * (PHI_Plasma - PHI_periphery) / (Vc * mLperL) -
      Vmax_nonlinearClearance * PHI_Plasma /
        (Km_nonlinearClearance + PHI_Plasma)
    d/dt(PHI_periphery) <- CLD_PHI * (PHI_Plasma - PHI_periphery) /
      (Vp * mLperL)

    # ---- Erythroid maturation chain (Eq. 14-27) ---------------------------
    d/dt(Progenitors) <- kprod_prog -
      k_progenitorsToPrecursors * Progenitors -
      k_deg_prog * Progenitors * (1 - fEPO)
    d/dt(Precursors) <-
      k_progenitorsToPrecursors * Progenitors * progenitorAmplification -
      kdeg_Precursors * Precursors - k_PrecursorsToPretics * Precursors
    d/dt(Retics_plasma) <- k_PrecursorsToPretics * Precursors -
      k_ReticsToRBCM * Retics_plasma - kdeg_Retics_Plasma * Retics_plasma
    d/dt(RBCM) <- k_ReticsToRBCM * Retics_plasma - kdeg_RBCM * RBCM

    # ---- Initial conditions ------------------------------------------------
    # Supplementary Parameter Dashboard, sheet "Species", column
    # InitialAmount. These are the authors' starting values, not the
    # steady state; the model is run to steady state before dosing
    # (Methods 3.4.3, "simulating steady-state dynamics").
    HIFa(0)              <- 0.4206       # ng/mL
    EPO_receptor(0)      <- 4e13         # molecule/L
    Progenitors(0)       <- 3.6303e9     # molecule/L
    Precursors(0)        <- 5.556988e10  # molecule/L
    Retics_plasma(0)     <- 5.556988e10  # molecule/L
    RBCM(0)              <- 4.8275e12    # molecule/L
  })
}
