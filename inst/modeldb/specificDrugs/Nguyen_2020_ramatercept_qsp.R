Nguyen_2020_ramatercept_qsp <- function() {
  description <- paste(
    "QSP. Myostatin / activin A / ActRIIB systems-pharmacology model (Nguyen",
    "2020) parameterised for ACE-031 (ramatercept), a soluble ActRIIB-Fc",
    "ligand trap, in healthy postmenopausal women.",
    "Structure identical to modellib('Nguyen_2020_fseeefc_qsp'): plasma plus",
    "muscle, anterior-pituitary and other-tissue interstitium linked by lymph",
    "flow; latent and mature myostatin and activin A; drug-ligand binding with",
    "complexes cleared like free drug; ligand binding to ActRIIB in muscle and",
    "pituitary.",
    "Compound PK parameters were fitted to digitised ACE-031 PK data and",
    "the FSH parameters (Vmax, RO50, Hill, kdeg) were estimated on the",
    "ACE-031 serum-FSH data (Figure 3); ligand-binding rate constants were",
    "derived from in vitro Kd (Sako 2010). The muscle-growth parameters are",
    "the adnectin estimates, used by the authors to forward-predict ACE-031",
    "muscle-volume change (Table 3). Deterministic; ligand, receptor and FSH",
    "states start at the drug-free steady state.",
    sep = " "
  )
  reference <- paste(
    "Nguyen HQ, Iskenderian A, Ehmann D, Jasper P, Zhang Z, Rong H, Welty D,",
    "Narayanan R. Leveraging Quantitative Systems Pharmacology Approach into",
    "Development of Human Recombinant Follistatin Fusion Protein for Duchenne",
    "Muscular Dystrophy. CPT Pharmacometrics Syst Pharmacol. 2020;9(6):342-352.",
    "doi:10.1002/psp4.12518. Model equations from the deposited R Markdown",
    "(run_QSPmodel.Rmd, Supplementary Information).",
    sep = " "
  )
  vignette <- "Nguyen_2020_myostatin_activin_dmd_qsp"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Every non-drug state is a paper-specific molecular species in a named
  # physiological space; none has a canonical compartment name. The drug in
  # plasma uses `central` (an amount, mg) and the subcutaneous reservoir uses
  # `depot`; all other states are concentrations in nM, as in the deposited code.
  paper_specific_compartments <- c(
    "drug_muscle",
    "drug_pituitary",
    "drug_other",
    "md_plasma",
    "md_muscle",
    "md_other",
    "ad_plasma",
    "ad_muscle",
    "ad_pituitary",
    "ad_other",
    "myo_plasma",
    "myo_muscle",
    "myo_other",
    "ppmyo_plasma",
    "ppmyo_muscle",
    "ppmyo_other",
    "act_plasma",
    "act_muscle",
    "act_pituitary",
    "act_other",
    "ppact_plasma",
    "ppact_muscle",
    "ppact_pituitary",
    "ppact_other",
    "actriib_muscle",
    "myo_actriib_muscle",
    "act_actriib_muscle",
    "actriib_pituitary",
    "act_actriib_pituitary",
    "muscle_growth",
    "fsh"
  )

  compartmentData <- list(
    depot = list(analyte = "ACE-031 (ActRIIB-Fc)", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ACE-031 (ActRIIB-Fc) (free, not ligand-bound)", units = "mg", specimen = "plasma", verified = TRUE),
    drug_muscle = list(analyte = "ACE-031 (ActRIIB-Fc) (free)", units = "nM", specimen = "tissue", verified = TRUE),
    drug_pituitary = list(analyte = "ACE-031 (ActRIIB-Fc) (free)", units = "nM", specimen = "tissue", verified = TRUE),
    drug_other = list(analyte = "ACE-031 (ActRIIB-Fc) (free)", units = "nM", specimen = "tissue", verified = TRUE),
    md_plasma = list(analyte = "ACE-031 (ActRIIB-Fc)-myostatin complex", units = "nM", specimen = "plasma", verified = TRUE),
    md_muscle = list(analyte = "ACE-031 (ActRIIB-Fc)-myostatin complex", units = "nM", specimen = "tissue", verified = TRUE),
    md_other = list(analyte = "ACE-031 (ActRIIB-Fc)-myostatin complex", units = "nM", specimen = "tissue", verified = TRUE),
    ad_plasma = list(analyte = "ACE-031 (ActRIIB-Fc)-activin A complex", units = "nM", specimen = "plasma", verified = TRUE),
    ad_muscle = list(analyte = "ACE-031 (ActRIIB-Fc)-activin A complex", units = "nM", specimen = "tissue", verified = TRUE),
    ad_pituitary = list(analyte = "ACE-031 (ActRIIB-Fc)-activin A complex", units = "nM", specimen = "tissue", verified = TRUE),
    ad_other = list(analyte = "ACE-031 (ActRIIB-Fc)-activin A complex", units = "nM", specimen = "tissue", verified = TRUE),
    myo_plasma = list(analyte = "mature myostatin (free)", units = "nM", specimen = "plasma", verified = TRUE),
    myo_muscle = list(analyte = "mature myostatin (free)", units = "nM", specimen = "tissue", verified = TRUE),
    myo_other = list(analyte = "mature myostatin (free)", units = "nM", specimen = "tissue", verified = TRUE),
    ppmyo_plasma = list(analyte = "myostatin-propeptide latent complex", units = "nM", specimen = "plasma", verified = TRUE),
    ppmyo_muscle = list(analyte = "myostatin-propeptide latent complex", units = "nM", specimen = "tissue", verified = TRUE),
    ppmyo_other = list(analyte = "myostatin-propeptide latent complex", units = "nM", specimen = "tissue", verified = TRUE),
    act_plasma = list(analyte = "mature activin A (free)", units = "nM", specimen = "plasma", verified = TRUE),
    act_muscle = list(analyte = "mature activin A (free)", units = "nM", specimen = "tissue", verified = TRUE),
    act_pituitary = list(analyte = "mature activin A (free)", units = "nM", specimen = "tissue", verified = TRUE),
    act_other = list(analyte = "mature activin A (free)", units = "nM", specimen = "tissue", verified = TRUE),
    ppact_plasma = list(analyte = "activin A-propeptide complex", units = "nM", specimen = "plasma", verified = TRUE),
    ppact_muscle = list(analyte = "activin A-propeptide complex", units = "nM", specimen = "tissue", verified = TRUE),
    ppact_pituitary = list(analyte = "activin A-propeptide complex", units = "nM", specimen = "tissue", verified = TRUE),
    ppact_other = list(analyte = "activin A-propeptide complex", units = "nM", specimen = "tissue", verified = TRUE),
    actriib_muscle = list(analyte = "ActRIIB receptor (free)", units = "nM", specimen = "tissue", verified = TRUE),
    myo_actriib_muscle = list(analyte = "myostatin-ActRIIB complex", units = "nM", specimen = "tissue", verified = TRUE),
    act_actriib_muscle = list(analyte = "activin A-ActRIIB complex", units = "nM", specimen = "tissue", verified = TRUE),
    actriib_pituitary = list(analyte = "ActRIIB receptor (free)", units = "nM", specimen = "tissue", verified = TRUE),
    act_actriib_pituitary = list(analyte = "activin A-ActRIIB complex", units = "nM", specimen = "tissue", verified = TRUE),
    muscle_growth = list(analyte = "muscle volume change from baseline", units = "%", specimen = "not applicable", verified = TRUE),
    fsh = list(analyte = "follicle-stimulating hormone", units = "ng/mL", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA,
    n_studies = 1,
    disease_state = "healthy postmenopausal women (single ascending-dose ACE-031 study, Attie 2013)",
    sex_female_pct = 100,
    dose_range = "single subcutaneous doses of 0.3 and 3 mg/kg (Figure 3); 1 and 3 mg/kg for muscle volume (Table 3)",
    notes = paste(
      "Parameters were estimated by fitting the QSP model to digitised",
      "aggregate PK and FSH data (Figure 3); subject counts are not reported.",
      "Doses in mg (mg/kg x 71 kg)."
    )
  )

  ini({
    # ---- System physiology (Table 1) --------------------------------------
    vplasma <- fixed(3.126); label("Plasma volume (L)") # Table 1 V_plasma = 3.126 L (Shah 2012)
    vmuscle <- fixed(3.91); label("Interstitial volume of muscle (L)") # Table 1 V_muscle = 3.91 L (Shah 2012)
    vpituitary <- fixed(5.4e-5); label("Interstitial volume of anterior pituitary (L)") # Table 1 V_pituitary = 5.4e-5 L
    lf_muscle <- fixed(0.33469); label("Lymph flow of muscle (L/h)") # Table 1 LF_muscle = 0.335 L/h; deposited code LF_muscle = 0.33469
    lf_pituitary <- fixed(1.34e-6); label("Lymph flow of pituitary (L/h)") # Table 1 LF_pituitary = 1.34e-6 L/h
    lf_other <- fixed(0.29689); label("Lymph flow of other tissues (L/h)") # Table 1 LF_other = 0.297 L/h; deposited code LF_other = 0.29689

    # ---- Myostatin (Table 1) ------------------------------------------------
    mw_myo <- fixed(25000); label("Molecular weight of mature myostatin (g/mol)") # Table 1 MW_Myo = 25,000
    mw_ppmyo <- fixed(80000); label("Molecular weight of myostatin-propeptide complex (g/mol)") # Table 1 MW_ppMyo = 80,000
    ksyn_ppmyo <- fixed(0.809821); label("Synthesis rate of myostatin-propeptide complex in muscle (nmol/h)") # Table 1 k_syn_ppMyo = 0.809 (unit printed 1/hour; enters the muscle balance as an amount rate, see vignette); deposited code 0.809821
    kcleave_ppmyo_muscle <- fixed(0.069768); label("Cleavage rate of myostatin-propeptide complex in muscle (1/h)") # Table 1 k_cleave_ppMyo_m = 0.0698; deposited code 0.069768
    kcleave_ppmyo_plasma <- fixed(0.01); label("Cleavage rate of myostatin-propeptide complex in plasma (1/h)") # Table 1 k_cleave_ppMyo_p = 0.01
    kdeg_ppmyo <- fixed(0.346574); label("Degradation rate of myostatin-propeptide complex, all spaces (1/h)") # Table 1 k_deg_ppMyo_m (or _p, _o) = 0.346; deposited code 0.346574
    kdeg_myo <- fixed(0.346574); label("Degradation rate of mature myostatin, all spaces (1/h)") # Table 1 k_deg_Myo_m (or _p, _o) = 0.346; deposited code 0.346574
    sigma_v_myo <- fixed(0.7); label("Vascular reflection coefficient of myostatin species (unitless)") # Table 1 sigma_m_V_ppMyo (or _myo) = 0.7; footnote a: same in other tissues
    sigma_is_myo <- fixed(0.2); label("Lymphatic reflection coefficient of myostatin species (unitless)") # Table 1 sigma_m_IS_ppMyo (or _myo) = 0.2; footnote a
    kon_myo_actriib <- fixed(1.328); label("Myostatin-ActRIIB association rate constant (1/(nM*h))") # Table 1 kon_Myo_ActRIIB = 1.328 (Sako 2010)
    koff_myo_actriib <- fixed(0.124); label("Myostatin-ActRIIB dissociation rate constant (1/h)") # Table 1 koff_Myo_ActRIIB = 0.124 (Sako 2010)

    # ---- Activin A (Table 1) ------------------------------------------------
    mw_act <- fixed(25000); label("Molecular weight of mature activin A (g/mol)") # Table 1 MW_Act = 25,000
    mw_ppact <- fixed(80000); label("Molecular weight of activin A-propeptide complex (g/mol)") # Table 1 MW_ppAct = 80,000
    ksyn_ppact_pituitary <- fixed(9.1e-6); label("Synthesis rate of activin A-propeptide complex in pituitary (nmol/h)") # Table 1 k_syn_ppAct_Pt = 9.1e-6 (unit printed 1/hour)
    ksyn_ppact_plasma <- fixed(0.120279); label("Synthesis rate of activin A-propeptide complex in plasma (nmol/h)") # Table 1 k_syn_ppAct_p = 0.12; deposited code 0.120279
    ksyn_ppact_muscle <- fixed(0.069841); label("Synthesis rate of activin A-propeptide complex in muscle (nmol/h)") # Table 1 k_syn_ppAct_m = 0.0698; deposited code 0.069841
    kcleave_ppact <- fixed(3113.17); label("Cleavage rate of activin A-propeptide complex, all producing spaces (1/h)") # Table 1 k_cleave_ppAct = k_cleave_ppAct_p = 3,113; deposited code 3113.17
    kdeg_ppact <- fixed(1.38629); label("Degradation rate of activin A-propeptide complex, all spaces (1/h)") # Table 1 k_deg_ppAct_Pt (or _p, _m, _o) = 1.386; deposited code 1.38629
    kdeg_act <- fixed(2.07944); label("Degradation rate of mature activin A, all spaces (1/h)") # Table 1 k_deg_Act_Pt (or _p, _m, _o) = 2.08; deposited code 2.07944
    sigma_v_act <- fixed(0.7); label("Vascular reflection coefficient of activin species (unitless)") # Table 1 sigma_m_V_ppAct (or _Act) = 0.7; footnote b: same in pituitary and other tissues
    sigma_is_act <- fixed(0.2); label("Lymphatic reflection coefficient of activin species (unitless)") # Table 1 sigma_m_IS_ppAct (or _Act) = 0.2; footnote b
    kon_act_actriib <- fixed(14.90); label("Activin A-ActRIIB association rate constant (1/(nM*h))") # Table 1 kon_Act_ActRIIB = 14.90 (Sako 2010)
    koff_act_actriib <- fixed(0.533); label("Activin A-ActRIIB dissociation rate constant (1/h)") # Table 1 koff_Act_ActRIIB = 0.533 (Sako 2010)

    # ---- ActRIIB (Table 1) --------------------------------------------------
    bl_actriib <- fixed(0.138); label("Total ActRIIB concentration in muscle and in pituitary (nM)") # Table 1 ActRIIB concentration = 0.138 nM

    # ---- ACE-031 compound parameters (Table 2, ACE-031 column) ------------
    mw_drug <- fixed(100000); label("Molecular weight of ACE-031 (g/mol)") # Table 2, ACE-031 column MW_Drug = 100,000
    lka <- fixed(log(0.011)); label("Subcutaneous absorption rate constant (1/h)") # Table 2, ACE-031 column k_a = 0.011
    lfdepot <- fixed(log(0.7)); label("Subcutaneous bioavailability (fraction)") # Table 2, ACE-031 column F = 0.7 (footnote b: typical s.c. value for therapeutic proteins)
    vother <- fixed(0.38); label("Volume of other tissues (L)") # Table 2, ACE-031 column V_other = 0.38 L
    kdeg_drug <- fixed(0.0019); label("Degradation rate of free ACE-031, all spaces (1/h)") # Table 2, ACE-031 column k_deg_drug = 0.0019
    kdeg_md <- fixed(0.0019); label("Degradation rate of ACE-031-myostatin complex (1/h)") # Table 2, ACE-031 column k_deg_MD = 0.0019 (assumed ~ k_deg_drug)
    kdeg_ad <- fixed(0.0019); label("Degradation rate of ACE-031-activin A complex (1/h)") # Table 2, ACE-031 column k_deg_AD = 0.0019 (assumed ~ k_deg_drug)
    kon_drug_myo <- fixed(0.6); label("ACE-031-myostatin association rate constant (1/(nM*h))") # Table 2, ACE-031 column kon_drug-Myo = 0.6
    koff_drug_myo <- fixed(0.0561); label("ACE-031-myostatin dissociation rate constant (1/h)") # Table 2, ACE-031 column koff_drug-Myo = 0.0561 (footnote e: from Kd in Sako 2010)
    kon_drug_act <- fixed(0.6); label("ACE-031-activin A association rate constant (1/(nM*h))") # Table 2, ACE-031 column kon_drug-Act = 0.6 (footnote i)
    koff_drug_act <- fixed(0.0214); label("ACE-031-activin A dissociation rate constant (1/h)") # Table 2, ACE-031 column koff_drug-Act = 0.0214 (footnote e: from Kd in Sako 2010)
    sigma_v_drug <- fixed(0.8); label("Vascular reflection coefficient of drug and drug complexes (unitless)") # Table 2, ACE-031 column sigma_m_V_drug = 0.80
    sigma_is_drug <- fixed(0.2); label("Lymphatic reflection coefficient of drug and drug complexes (unitless)") # Table 2, ACE-031 column sigma_m_IS_drug = 0.2 (footnote h: same in pituitary and other tissues)

    # ---- Muscle and FSH pharmacodynamics (Table 2, adnectin / ACE-031 columns)
    vmax_muscle <- fixed(0.005); label("Maximal rate of muscle-volume increase (%/h)") # Table 2 V_max_muscle = 0.005 (adnectin estimate; used to forward-predict ACE-031 muscle volume, Table 3)
    h_muscle <- fixed(6); label("Hill coefficient of the ActRIIB occupancy-muscle link (unitless)") # Table 2 h_muscle = 6
    ro50_muscle <- fixed(23.36); label("Muscle ActRIIB occupancy giving half-maximal effect (%)") # Table 2 RO_50_muscle = 23.36
    kdeg_muscle <- fixed(0.0002); label("First-order loss rate of the muscle-volume change (1/h)") # Table 2 k_deg_muscle = 0.0002
    vmax_fsh <- fixed(25.75); label("Maximal FSH production rate (ng/mL/h)") # Table 2 V_max_FSH = 25.75 (ACE-031 fit)
    h_fsh <- fixed(1.80); label("Hill coefficient of the pituitary occupancy-FSH link (unitless)") # Table 2 h_FSH = 1.80
    ro50_fsh <- fixed(21.78); label("Pituitary ActRIIB occupancy giving half-maximal FSH production (%)") # Table 2 RO_50_FSH = 21.78
    kdeg_fsh <- fixed(1.00); label("First-order FSH elimination rate (1/h)") # Table 2 k_deg_FSH = 1.00

    # Deterministic systems model: no residual-error model is reported.
    propSd <- fixed(0); label("Proportional residual error (fraction)")
  })

  model({
    ka <- exp(lka)
    fdepot <- exp(lfdepot)

    kd_myo_actriib <- koff_myo_actriib / kon_myo_actriib
    kd_act_actriib <- koff_act_actriib / kon_act_actriib

    # ---- Drug-free steady state (initial conditions) -------------------------
    # The deposited code starts every ligand at 1e-6 nM and runs the system
    # drug-free for 1000 h before the first dose. Without drug the receptor
    # binding fluxes vanish at steady state, so every ligand steady state solves
    # a linear plasma-hub system: for each tissue j,
    #   X_j = (S_j + LF_j (1 - sigma_V) X_p) / (k_j V_j + LF_j (1 - sigma_IS)),
    # and the plasma balance then gives X_p. The receptor then partitions by
    # equilibrium binding with the total receptor conserved.
    # Myostatin-propeptide complex (made in muscle; cleaved in muscle and plasma)
    dd_ppmyo_m <- (kcleave_ppmyo_muscle + kdeg_ppmyo) * vmuscle + lf_muscle * (1 - sigma_is_myo)
    dd_ppmyo_o <- kdeg_ppmyo * vother + lf_other * (1 - sigma_is_myo)
    ss_ppmyo_p <- (lf_muscle * (1 - sigma_is_myo) * ksyn_ppmyo / dd_ppmyo_m) /
      ((kcleave_ppmyo_plasma + kdeg_ppmyo) * vplasma +
        lf_muscle * (1 - sigma_v_myo) * (1 - lf_muscle * (1 - sigma_is_myo) / dd_ppmyo_m) +
        lf_other * (1 - sigma_v_myo) * (1 - lf_other * (1 - sigma_is_myo) / dd_ppmyo_o))
    ss_ppmyo_m <- (ksyn_ppmyo + lf_muscle * (1 - sigma_v_myo) * ss_ppmyo_p) / dd_ppmyo_m
    ss_ppmyo_o <- lf_other * (1 - sigma_v_myo) * ss_ppmyo_p / dd_ppmyo_o
    # Mature myostatin
    dd_myo_m <- kdeg_myo * vmuscle + lf_muscle * (1 - sigma_is_myo)
    dd_myo_o <- kdeg_myo * vother + lf_other * (1 - sigma_is_myo)
    ss_myo_p <- (kcleave_ppmyo_plasma * ss_ppmyo_p * vplasma +
      lf_muscle * (1 - sigma_is_myo) * kcleave_ppmyo_muscle * ss_ppmyo_m * vmuscle / dd_myo_m) /
      (kdeg_myo * vplasma +
        lf_muscle * (1 - sigma_v_myo) * (1 - lf_muscle * (1 - sigma_is_myo) / dd_myo_m) +
        lf_other * (1 - sigma_v_myo) * (1 - lf_other * (1 - sigma_is_myo) / dd_myo_o))
    ss_myo_m <- (kcleave_ppmyo_muscle * ss_ppmyo_m * vmuscle + lf_muscle * (1 - sigma_v_myo) * ss_myo_p) / dd_myo_m
    ss_myo_o <- lf_other * (1 - sigma_v_myo) * ss_myo_p / dd_myo_o
    # Activin A-propeptide complex (made in plasma, muscle and pituitary)
    dd_ppact_m <- (kcleave_ppact + kdeg_ppact) * vmuscle + lf_muscle * (1 - sigma_is_act)
    dd_ppact_t <- (kcleave_ppact + kdeg_ppact) * vpituitary + lf_pituitary * (1 - sigma_is_act)
    dd_ppact_o <- kdeg_ppact * vother + lf_other * (1 - sigma_is_act)
    ss_ppact_p <- (ksyn_ppact_plasma +
      lf_muscle * (1 - sigma_is_act) * ksyn_ppact_muscle / dd_ppact_m +
      lf_pituitary * (1 - sigma_is_act) * ksyn_ppact_pituitary / dd_ppact_t) /
      ((kcleave_ppact + kdeg_ppact) * vplasma +
        lf_muscle * (1 - sigma_v_act) * (1 - lf_muscle * (1 - sigma_is_act) / dd_ppact_m) +
        lf_pituitary * (1 - sigma_v_act) * (1 - lf_pituitary * (1 - sigma_is_act) / dd_ppact_t) +
        lf_other * (1 - sigma_v_act) * (1 - lf_other * (1 - sigma_is_act) / dd_ppact_o))
    ss_ppact_m <- (ksyn_ppact_muscle + lf_muscle * (1 - sigma_v_act) * ss_ppact_p) / dd_ppact_m
    ss_ppact_t <- (ksyn_ppact_pituitary + lf_pituitary * (1 - sigma_v_act) * ss_ppact_p) / dd_ppact_t
    ss_ppact_o <- lf_other * (1 - sigma_v_act) * ss_ppact_p / dd_ppact_o
    # Mature activin A
    dd_act_m <- kdeg_act * vmuscle + lf_muscle * (1 - sigma_is_act)
    dd_act_t <- kdeg_act * vpituitary + lf_pituitary * (1 - sigma_is_act)
    dd_act_o <- kdeg_act * vother + lf_other * (1 - sigma_is_act)
    ss_act_p <- (kcleave_ppact * ss_ppact_p * vplasma +
      lf_muscle * (1 - sigma_is_act) * kcleave_ppact * ss_ppact_m * vmuscle / dd_act_m +
      lf_pituitary * (1 - sigma_is_act) * kcleave_ppact * ss_ppact_t * vpituitary / dd_act_t) /
      (kdeg_act * vplasma +
        lf_muscle * (1 - sigma_v_act) * (1 - lf_muscle * (1 - sigma_is_act) / dd_act_m) +
        lf_pituitary * (1 - sigma_v_act) * (1 - lf_pituitary * (1 - sigma_is_act) / dd_act_t) +
        lf_other * (1 - sigma_v_act) * (1 - lf_other * (1 - sigma_is_act) / dd_act_o))
    ss_act_m <- (kcleave_ppact * ss_ppact_m * vmuscle + lf_muscle * (1 - sigma_v_act) * ss_act_p) / dd_act_m
    ss_act_t <- (kcleave_ppact * ss_ppact_t * vpituitary + lf_pituitary * (1 - sigma_v_act) * ss_act_p) / dd_act_t
    ss_act_o <- lf_other * (1 - sigma_v_act) * ss_act_p / dd_act_o
    # Receptors at binding equilibrium, total conserved at bl_actriib
    ss_r_m <- bl_actriib / (1 + ss_myo_m / kd_myo_actriib + ss_act_m / kd_act_actriib)
    ss_r_t <- bl_actriib / (1 + ss_act_t / kd_act_actriib)
    ss_ro_t <- 100 * (bl_actriib - ss_r_t) / bl_actriib

    ppmyo_plasma(0) <- ss_ppmyo_p
    ppmyo_muscle(0) <- ss_ppmyo_m
    ppmyo_other(0) <- ss_ppmyo_o
    myo_plasma(0) <- ss_myo_p
    myo_muscle(0) <- ss_myo_m
    myo_other(0) <- ss_myo_o
    ppact_plasma(0) <- ss_ppact_p
    ppact_muscle(0) <- ss_ppact_m
    ppact_pituitary(0) <- ss_ppact_t
    ppact_other(0) <- ss_ppact_o
    act_plasma(0) <- ss_act_p
    act_muscle(0) <- ss_act_m
    act_pituitary(0) <- ss_act_t
    act_other(0) <- ss_act_o
    actriib_muscle(0) <- ss_r_m
    myo_actriib_muscle(0) <- ss_r_m * ss_myo_m / kd_myo_actriib
    act_actriib_muscle(0) <- ss_r_m * ss_act_m / kd_act_actriib
    actriib_pituitary(0) <- ss_r_t
    act_actriib_pituitary(0) <- ss_r_t * ss_act_t / kd_act_actriib
    fsh(0) <- vmax_fsh * ss_ro_t^h_fsh / (ss_ro_t^h_fsh + ro50_fsh^h_fsh) / kdeg_fsh

    # ---- Free drug concentration in plasma (nM) -------------------------------
    drug_plasma <- central * 1e6 / mw_drug / vplasma

    # ---- Transport fluxes, plasma -> tissue interstitium (nmol/h) ------------
    f_drug_m <- lf_muscle * ((1 - sigma_v_drug) * drug_plasma - (1 - sigma_is_drug) * drug_muscle)
    f_drug_t <- lf_pituitary * ((1 - sigma_v_drug) * drug_plasma - (1 - sigma_is_drug) * drug_pituitary)
    f_drug_o <- lf_other * ((1 - sigma_v_drug) * drug_plasma - (1 - sigma_is_drug) * drug_other)
    f_md_m <- lf_muscle * ((1 - sigma_v_drug) * md_plasma - (1 - sigma_is_drug) * md_muscle)
    f_md_o <- lf_other * ((1 - sigma_v_drug) * md_plasma - (1 - sigma_is_drug) * md_other)
    f_ad_m <- lf_muscle * ((1 - sigma_v_drug) * ad_plasma - (1 - sigma_is_drug) * ad_muscle)
    f_ad_t <- lf_pituitary * ((1 - sigma_v_drug) * ad_plasma - (1 - sigma_is_drug) * ad_pituitary)
    f_ad_o <- lf_other * ((1 - sigma_v_drug) * ad_plasma - (1 - sigma_is_drug) * ad_other)
    f_ppmyo_m <- lf_muscle * ((1 - sigma_v_myo) * ppmyo_plasma - (1 - sigma_is_myo) * ppmyo_muscle)
    f_ppmyo_o <- lf_other * ((1 - sigma_v_myo) * ppmyo_plasma - (1 - sigma_is_myo) * ppmyo_other)
    f_myo_m <- lf_muscle * ((1 - sigma_v_myo) * myo_plasma - (1 - sigma_is_myo) * myo_muscle)
    f_myo_o <- lf_other * ((1 - sigma_v_myo) * myo_plasma - (1 - sigma_is_myo) * myo_other)
    f_ppact_m <- lf_muscle * ((1 - sigma_v_act) * ppact_plasma - (1 - sigma_is_act) * ppact_muscle)
    f_ppact_t <- lf_pituitary * ((1 - sigma_v_act) * ppact_plasma - (1 - sigma_is_act) * ppact_pituitary)
    f_ppact_o <- lf_other * ((1 - sigma_v_act) * ppact_plasma - (1 - sigma_is_act) * ppact_other)
    f_act_m <- lf_muscle * ((1 - sigma_v_act) * act_plasma - (1 - sigma_is_act) * act_muscle)
    f_act_t <- lf_pituitary * ((1 - sigma_v_act) * act_plasma - (1 - sigma_is_act) * act_pituitary)
    f_act_o <- lf_other * ((1 - sigma_v_act) * act_plasma - (1 - sigma_is_act) * act_other)

    # ---- Net binding fluxes (nmol/h) ------------------------------------------
    b_md_p <- (kon_drug_myo * drug_plasma * myo_plasma - koff_drug_myo * md_plasma) * vplasma
    b_md_m <- (kon_drug_myo * drug_muscle * myo_muscle - koff_drug_myo * md_muscle) * vmuscle
    b_md_o <- (kon_drug_myo * drug_other * myo_other - koff_drug_myo * md_other) * vother
    b_ad_p <- (kon_drug_act * drug_plasma * act_plasma - koff_drug_act * ad_plasma) * vplasma
    b_ad_m <- (kon_drug_act * drug_muscle * act_muscle - koff_drug_act * ad_muscle) * vmuscle
    b_ad_t <- (kon_drug_act * drug_pituitary * act_pituitary - koff_drug_act * ad_pituitary) * vpituitary
    b_ad_o <- (kon_drug_act * drug_other * act_other - koff_drug_act * ad_other) * vother
    b_mr_m <- (kon_myo_actriib * actriib_muscle * myo_muscle - koff_myo_actriib * myo_actriib_muscle) * vmuscle
    b_ar_m <- (kon_act_actriib * actriib_muscle * act_muscle - koff_act_actriib * act_actriib_muscle) * vmuscle
    b_ar_t <- (kon_act_actriib * actriib_pituitary * act_pituitary - koff_act_actriib * act_actriib_pituitary) * vpituitary

    # ---- Receptor occupancy (%) -----------------------------------------------
    r_total_m <- actriib_muscle + myo_actriib_muscle + act_actriib_muscle
    ro_muscle_myo <- 100 * myo_actriib_muscle / r_total_m
    ro_muscle_act <- 100 * act_actriib_muscle / r_total_m
    ro_muscle <- ro_muscle_myo + ro_muscle_act
    ro_pituitary <- 100 * act_actriib_pituitary / (act_actriib_pituitary + actriib_pituitary)

    # ---- Drug ---------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    f(depot) <- fdepot
    # Plasma drug balance in nmol/h, converted to mg/h for the amount state
    d/dt(central) <- ka * depot - (kdeg_drug * drug_plasma * vplasma + f_drug_m + f_drug_t + f_drug_o +
      b_md_p + b_ad_p) * mw_drug / 1e6
    d/dt(drug_muscle) <- (f_drug_m - b_md_m - b_ad_m - kdeg_drug * drug_muscle * vmuscle) / vmuscle
    d/dt(drug_pituitary) <- (f_drug_t - b_ad_t - kdeg_drug * drug_pituitary * vpituitary) / vpituitary
    d/dt(drug_other) <- (f_drug_o - b_md_o - b_ad_o - kdeg_drug * drug_other * vother) / vother

    # ---- Drug-ligand complexes --------------------------------------------------
    d/dt(md_plasma) <- (b_md_p - kdeg_md * md_plasma * vplasma - f_md_m - f_md_o) / vplasma
    d/dt(md_muscle) <- (b_md_m - kdeg_md * md_muscle * vmuscle + f_md_m) / vmuscle
    d/dt(md_other) <- (b_md_o - kdeg_md * md_other * vother + f_md_o) / vother
    d/dt(ad_plasma) <- (b_ad_p - kdeg_ad * ad_plasma * vplasma - f_ad_m - f_ad_t - f_ad_o) / vplasma
    d/dt(ad_muscle) <- (b_ad_m - kdeg_ad * ad_muscle * vmuscle + f_ad_m) / vmuscle
    d/dt(ad_pituitary) <- (b_ad_t - kdeg_ad * ad_pituitary * vpituitary + f_ad_t) / vpituitary
    d/dt(ad_other) <- (b_ad_o - kdeg_ad * ad_other * vother + f_ad_o) / vother

    # ---- Myostatin ------------------------------------------------------------
    d/dt(ppmyo_plasma) <- (-(kcleave_ppmyo_plasma + kdeg_ppmyo) * ppmyo_plasma * vplasma - f_ppmyo_m - f_ppmyo_o) / vplasma
    d/dt(ppmyo_muscle) <- (ksyn_ppmyo - (kcleave_ppmyo_muscle + kdeg_ppmyo) * ppmyo_muscle * vmuscle + f_ppmyo_m) / vmuscle
    d/dt(ppmyo_other) <- (-kdeg_ppmyo * ppmyo_other * vother + f_ppmyo_o) / vother
    d/dt(myo_plasma) <- (kcleave_ppmyo_plasma * ppmyo_plasma * vplasma - kdeg_myo * myo_plasma * vplasma -
      f_myo_m - f_myo_o - b_md_p) / vplasma
    d/dt(myo_muscle) <- (kcleave_ppmyo_muscle * ppmyo_muscle * vmuscle - kdeg_myo * myo_muscle * vmuscle +
      f_myo_m - b_md_m - b_mr_m) / vmuscle
    d/dt(myo_other) <- (-kdeg_myo * myo_other * vother + f_myo_o - b_md_o) / vother

    # ---- Activin A ------------------------------------------------------------
    d/dt(ppact_plasma) <- (ksyn_ppact_plasma - (kcleave_ppact + kdeg_ppact) * ppact_plasma * vplasma -
      f_ppact_m - f_ppact_t - f_ppact_o) / vplasma
    d/dt(ppact_muscle) <- (ksyn_ppact_muscle - (kcleave_ppact + kdeg_ppact) * ppact_muscle * vmuscle + f_ppact_m) / vmuscle
    d/dt(ppact_pituitary) <- (ksyn_ppact_pituitary - (kcleave_ppact + kdeg_ppact) * ppact_pituitary * vpituitary +
      f_ppact_t) / vpituitary
    d/dt(ppact_other) <- (-kdeg_ppact * ppact_other * vother + f_ppact_o) / vother
    d/dt(act_plasma) <- (kcleave_ppact * ppact_plasma * vplasma - kdeg_act * act_plasma * vplasma -
      f_act_m - f_act_t - f_act_o - b_ad_p) / vplasma
    d/dt(act_muscle) <- (kcleave_ppact * ppact_muscle * vmuscle - kdeg_act * act_muscle * vmuscle +
      f_act_m - b_ad_m - b_ar_m) / vmuscle
    d/dt(act_pituitary) <- (kcleave_ppact * ppact_pituitary * vpituitary - kdeg_act * act_pituitary * vpituitary +
      f_act_t - b_ad_t - b_ar_t) / vpituitary
    d/dt(act_other) <- (-kdeg_act * act_other * vother + f_act_o - b_ad_o) / vother

    # ---- ActRIIB receptor in muscle and pituitary -------------------------------
    d/dt(actriib_muscle) <- (-b_mr_m - b_ar_m) / vmuscle
    d/dt(myo_actriib_muscle) <- b_mr_m / vmuscle
    d/dt(act_actriib_muscle) <- b_ar_m / vmuscle
    d/dt(actriib_pituitary) <- -b_ar_t / vpituitary
    d/dt(act_actriib_pituitary) <- b_ar_t / vpituitary

    # ---- Pharmacodynamics -------------------------------------------------------
    # Muscle: SI equation for d(%MuscleGrowth)/dt; starts at 0 when dosing
    # starts (deposited code switches it on at the first dose).
    d/dt(muscle_growth) <- vmax_muscle * (1 - ro_muscle^h_muscle / (ro_muscle^h_muscle + ro50_muscle^h_muscle)) -
      kdeg_muscle * muscle_growth
    # FSH: deposited-code form (production rises with pituitary occupancy);
    # the SI prints a (1 - Hill) term here, see vignette.
    d/dt(fsh) <- vmax_fsh * ro_pituitary^h_fsh / (ro_pituitary^h_fsh + ro50_fsh^h_fsh) - kdeg_fsh * fsh

    # ---- Outputs ----------------------------------------------------------------
    # Cc is TOTAL drug in plasma (free + ligand-bound), the quantity the
    # serum assays measured and Figure 4c plots; Cfree is ligand-free drug
    # (the deposited code's Drug_p_ngml). They diverge only when the drug
    # concentration approaches the nM-range ligand concentrations.
    Cc <- (drug_plasma + md_plasma + ad_plasma) * mw_drug / 1000
    Cfree <- drug_plasma * mw_drug / 1000
    myo_free_plasma <- myo_plasma * mw_myo / 1000
    act_free_plasma <- act_plasma * mw_act / 1000
    # Latent + mature (ng/mL), the 'ppMyo+Myo' circulating level of Table S1
    myo_total_plasma <- myo_free_plasma + ppmyo_plasma * mw_ppmyo / 1000
    act_total_plasma <- act_free_plasma + ppact_plasma * mw_ppact / 1000

    Cc ~ prop(propSd)
  })
}
