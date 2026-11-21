Gevertz_2025_teclistamab_qsp <- function() {
  description <- "QSP. Semi-mechanistic PK/PD model of teclistamab, a BCMA x CD3 bispecific T-cell engager, in relapsed/refractory multiple myeloma, built to locate the peak of the bell-shaped trimer-versus-dose curve. Free drug follows two-compartment disposition with first-order subcutaneous absorption and a clearance that decays exponentially from the first dose, and distributes into a bone-marrow (tumour microenvironment) compartment. Drug binds CD3 (on T cells, produced in a thymus compartment) and soluble BCMA in plasma, and CD3, soluble BCMA and membrane-bound BCMA in bone marrow, by mass action, forming dimers and trimers. Only the bone-marrow drug-CD3-membrane BCMA trimer (`trimer`, the pharmacologically active trimer, PAT) is efficacious. The initial conditions are the authors' pre-equilibration starting values, not the drug-free steady state: the authors run the model drug-free for 50 days before the first dose, so dose at day 50 to reproduce the published simulations. Inter-individual variability on the five virtual-patient attributes (Cl1, kshed, baseline membrane BCMA, k1bm and baseline bone-marrow CD3) reproduces the paper's virtual population. Parameter values and equations are the published Table S1 and Eqs (1)-(18). The published figures were produced by the authors' code, which uses a lower membrane BCMA internalisation rate (0.554/day, not 2/day) and CD3 internalisation rate (1.584/day, not 16.6/day) among other differences; with the Table S1 values the median PAT peak of the virtual population moves from about 1.5-1.8 mg/kg to about 5 mg/kg weekly. The validation article rebuilds the code configuration and reproduces the paper's 1.8 mg/kg Q1W and 3 mg/kg Q2W optima."
  reference <- paste(
    "Gevertz JL, Kareva I. Using Virtual Patients to Evaluate Dosing",
    "Strategies for Bispecifics With a Bell-Shaped Efficacy Curve.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(12):2139-2148.",
    "doi:10.1002/psp4.70105.",
    "Model equations (1)-(18) and Table S1 are in Supporting Information Data S1;",
    "the authors' MATLAB implementation is at https://github.com/jgevertz/teclistamab",
    "(commit e6dc5479). The time-varying clearance component and the 74 kg",
    "reference weight come from Miao et al. (2023) Target Oncol 18:667-684",
    "(reference 5 of Data S1), see modellib('Miao_2023_teclistamab')."
  )
  vignette <- "Gevertz_2025_teclistamab"

  # Time is days. Every state except `depot` holds a CONCENTRATION in nM,
  # exactly as Equations (1)-(18) are written; `depot` holds the subcutaneous
  # dose as an AMOUNT in nmol and is absorbed into the central (plasma)
  # concentration via ka * depot / vc. Teclistamab is 146 kDa (main text,
  # Methods), so 1 mg = 1e6 / 146000 nmol = 6.849 nmol, and a mg/kg dose is
  # converted at the 74 kg reference weight the authors' code uses (V1 = 3 L
  # is coded as 3000 / 74 mL/kg; the paper itself states no body weight).
  units <- list(time = "day", dosing = "nmol", concentration = "nM")

  # States with no registered canonical: the thymus and bone-marrow free CD3
  # pools, the free membrane BCMA pool, and the two pharmacologically inactive
  # drug-CD3-soluble BCMA trimers.
  paper_specific_compartments <- c(
    "target_cd3_thymus",
    "target_cd3_bonemarrow",
    "target_bcma_bonemarrow",
    "complex_cd3_sbcma",
    "complex_cd3_sbcma_bonemarrow"
  )

  compartmentData <- list(
    depot = list(
      analyte = "teclistamab (Data S1 Dsc)",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "teclistamab, free (Data S1 Dp)", units = "nM", specimen = "plasma", verified = TRUE),
    peripheral1 = list(
      analyte = "teclistamab, free (Data S1 Dper)",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    target_cd3_thymus = list(analyte = "CD3, free (Data S1 CD3th)", units = "nM", specimen = "tissue", verified = TRUE),
    target_cd3_central = list(analyte = "CD3, free (Data S1 CD3p)", units = "nM", specimen = "plasma", verified = TRUE),
    target = list(analyte = "soluble BCMA, free (Data S1 sBp)", units = "nM", specimen = "plasma", verified = TRUE),
    complex = list(
      analyte = "teclistamab-soluble BCMA dimer (Data S1 DsBp)",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    drug_cd3_central = list(
      analyte = "teclistamab-CD3 dimer (Data S1 D3p)",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    complex_cd3_sbcma = list(
      analyte = "teclistamab-CD3-soluble BCMA inactive trimer (Data S1 DsB3p)",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    target_bcma_bonemarrow = list(
      analyte = "membrane-bound BCMA, free (Data S1 mBbm)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    target_bonemarrow = list(
      analyte = "soluble BCMA, free (Data S1 sBbm)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    target_cd3_bonemarrow = list(
      analyte = "CD3, free (Data S1 CD3bm)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    bonemarrow = list(analyte = "teclistamab, free (Data S1 Dbm)", units = "nM", specimen = "tumor", verified = TRUE),
    drug_cd3_bonemarrow = list(
      analyte = "teclistamab-CD3 dimer (Data S1 D3bm)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    complex_bonemarrow = list(
      analyte = "teclistamab-soluble BCMA dimer (Data S1 DsBbm)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    drug_bcma_bonemarrow = list(
      analyte = "teclistamab-membrane BCMA dimer (Data S1 DmBbm)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    complex_cd3_sbcma_bonemarrow = list(
      analyte = "teclistamab-CD3-soluble BCMA inactive trimer (Data S1 DsB3bm)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    trimer = list(
      analyte = "teclistamab-CD3-membrane BCMA pharmacologically active trimer (Data S1 DmB3bm, PAT)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "relapsed/refractory multiple myeloma (virtual patients)",
    dose_range = "0.1-10 mg/kg SC Q1W or Q2W after step-up doses of 0.06 and 0.3 mg/kg one week apart (Figure 5); 0.72 and 1.5 mg/kg SC Q1W for the PK benchmark (Figure 2)",
    notes = paste(
      "No model was fitted to subject-level data. Structural PK values are taken or derived from the",
      "MajesTEC-1 population PK analysis (Miao 2023) and from Girgis 2022, whose digitised mean",
      "plasma PK at 0.72 and 1.5 mg/kg SC is reproduced in Figure 2A; target-binding and turnover",
      "values are from the literature cited in Table S1. The virtual population is 1200 virtual",
      "patients, each defined by independent lognormal draws of Cl1, kshed, baseline membrane BCMA",
      "in bone marrow, k1bm and baseline bone-marrow CD3 (Figure S2), chosen to mimic the size of a",
      "Phase 3 trial."
    )
  )

  ini({
    # Drug PK -- Table S1 'PK parameters'. Every value is a fixed literature
    # or authors' estimate; nothing was estimated against data in this paper.
    lka <- fixed(log(0.2)); label("SC absorption rate constant k01 (1/day)") # Table S1 k01 = 0.2 1/d
    lvc <- fixed(log(3)); label("Central (plasma) volume V1 (L)") # Table S1 V1 = '3-4' L; k10 = Cl1/V1 = 0.17 and k12 = Cl2/V1 = 0.13 (Table S1) both imply V1 = 3 L
    lvp <- fixed(log(1.34)); label("Peripheral volume V2 (L)") # Table S1 V2 = 1.34 L
    lq <- fixed(log(0.4)); label("Intercompartmental clearance Cl2 (L/day)") # Table S1 Cl2 = 0.4 L/d (k12 = Cl2/V1, k21 = Cl2/V2)
    lcl_exp_inf <- fixed(log(0.5)); label("Time-independent clearance Cl1 (L/day)") # Table S1 Cl1 = 0.5 L/d
    lcl_exp_component <- fixed(log(0.547)); label("Time-dependent clearance component Clvar at the first dose (L/day)") # Table S1 Cl1 footnote names Clvar but does not print it: taken from its ref 5, Miao 2023 Table 2 CL2 = 0.547 L/d (modellib(Miao_2023_teclistamab)); the authors code uses cl_tvar = 0.547 L/d
    lcl_exp_kdes <- fixed(log(0.1)); label("Decay rate constant delta of the time-dependent clearance (1/day)") # Table S1 Cl1 footnote: delta = 0.1
    vbm <- fixed(2.2); label("Bone-marrow volume Vbm (L)") # Table S1 Vbm = 2.2 L
    lktr_ab_c_bm <- fixed(log(0.1)); label("Drug distribution rate constant, plasma to bone marrow, k1bm (1/day)") # Table S1 k1bm = 0.1 1/d
    ktr_ab_bm_c <- fixed(0.001); label("Drug distribution rate constant, bone marrow to plasma, kbm1 (1/day)") # Table S1 kbm1 = 0.001 1/d

    # Binding -- Table S1. The bone-marrow CD3 rate constants kon3,bm / koff3,bm
    # of Eqs (12)-(18) are not tabulated separately; they equal the plasma values.
    kon_cd3 <- fixed(10); label("Drug-CD3 association rate constant kon3 (1/nM/day)") # Table S1 kon3,p = 10 1/nM/d
    koff_cd3 <- fixed(285); label("Drug-CD3 dissociation rate constant koff3 (1/day)") # Table S1 koff3,p = 285 1/d
    kon_sbcma <- fixed(111); label("Drug-soluble BCMA association rate constant kon,sB (1/nM/day)") # Table S1 kon,sB = 111 1/nM/d
    koff_sbcma <- fixed(20); label("Drug-soluble BCMA dissociation rate constant koff,sB (1/day)") # Table S1 koff,sB = 20 1/d
    kon_mbcma <- fixed(111); label("Drug-membrane BCMA association rate constant kon,mB (1/nM/day)") # Table S1 kon,mB = 111 1/nM/d
    koff_mbcma <- fixed(20); label("Drug-membrane BCMA dissociation rate constant koff,mB (1/day)") # Table S1 koff,mB = 20 1/d

    # CD3 -- Table S1 'CD3 parameters'
    ktr_cd3_thymus_c <- fixed(0.5); label("CD3 distribution rate constant, thymus to plasma, kth1 (1/day)") # Table S1 kth1 = 0.5 1/d
    kint_cd3 <- fixed(16.6); label("CD3 internalisation rate constant kint3 (1/day)") # Table S1 kint3 = 16.6 1/d
    ktr_cd3_bm_c <- fixed(0.1); label("Free CD3 distribution rate constant, bone marrow to plasma, k3,bm1 (1/day)") # Table S1 k3,bm1 = 0.1 1/d
    ktr_cd3_c_bm <- fixed(0.01); label("Free CD3 distribution rate constant, plasma to bone marrow, k3,1bm (1/day)") # Table S1 k3,1bm = 0.01 1/d

    # BCMA -- Table S1 'BCMA parameters' and k10,sB / kB,bm1 / kB,1bm / kout
    lkshed <- fixed(log(0.05)); label("Membrane BCMA shedding rate constant kshed (1/day)") # Table S1 kshed = 0.05 1/d
    kint_bcma <- fixed(2); label("Membrane BCMA internalisation rate constant kint,mB (1/day)") # Table S1 kint_mB = 2 1/d
    kdeg_sbcma <- fixed(0.04); label("Systemic clearance rate constant of soluble BCMA k10,sB (1/day)") # Table S1 k10,sB = 0.04 1/d
    ktr_sbcma_bm_c <- fixed(0.1); label("Soluble BCMA distribution rate constant, bone marrow to plasma, kB,bm1 (1/day)") # Table S1 kB,bm1 = 0.1 1/d
    ktr_sbcma_c_bm <- fixed(0.1); label("Soluble BCMA distribution rate constant, plasma to bone marrow, kB,1bm (1/day)") # Table S1 kB,1bm = 0.1 1/d
    kel_complex_bm <- fixed(0.01); label("Clearance rate constant of the bone-marrow inactive trimer kout (1/day)") # Table S1 kout = 0.01 1/d

    # Nonzero initial conditions -- Table S1. These are the values the authors
    # start their 50-day drug-free pre-equilibration from.
    bl_cd3_thymus <- fixed(0.001); label("Initial free CD3 in thymus CD3th(0) (nM)") # Table S1 CD3th(0) = 0.001 nM
    bl_cd3_central <- fixed(0.01); label("Initial free CD3 in plasma CD3p(0) (nM)") # Table S1 CD3p(0) = 0.01 nM
    lbl_cd3_bm <- fixed(log(0.1)); label("Initial free CD3 in bone marrow CD3bm(0) (nM)") # Table S1 CD3bm(0) = 0.1 nM
    bl_sbcma <- fixed(20); label("Initial free soluble BCMA in plasma sBp(0) (nM)") # Table S1 sBp(0) = 20 nM
    lbl_mbcma_bm <- fixed(log(50)); label("Initial membrane-bound BCMA in bone marrow mBbm(0) (nM)") # Table S1 mBbm(0) = 50 nM
    bl_sbcma_bm <- fixed(3); label("Initial soluble BCMA in bone marrow sBbm(0) (nM)") # Table S1 sBbm(0) = 3 nM

    # Virtual-patient variability -- Figure S2 lognormal mean mu and SD sigma,
    # converted to a log-scale variance log(1 + (sigma/mu)^2).
    etalcl_exp_inf ~ fixed(0.518232) # Figure S2A Cl1 mu = 6.0676, sigma = 5 (mL/kg/day): log(1 + (5/6.0676)^2)
    etalkshed ~ fixed(1.085189) # Figure S2B kshed mu = 0.05, sigma = 0.07: log(1 + (0.07/0.05)^2)
    etalbl_mbcma_bm ~ fixed(0.148420) # Figure S2C mBbm0 mu = 50, sigma = 20: log(1 + (20/50)^2)
    etalktr_ab_c_bm ~ fixed(0.039221) # Figure S2D k1bm mu = 0.1, sigma = 0.02: log(1 + (0.02/0.1)^2)
    etalbl_cd3_bm ~ fixed(0.009950) # Figure S2E CD3bm0 mu = 0.1, sigma = 0.01: log(1 + (0.01/0.1)^2)

    # The paper's simulations are deterministic; no residual error is reported.
    propSd <- fixed(0); label("Proportional residual error (fraction; not reported in source)")
  })

  model({
    # Individual parameters
    ka <- exp(lka)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    cl_exp_inf <- exp(lcl_exp_inf + etalcl_exp_inf)
    cl_exp_component <- exp(lcl_exp_component)
    cl_exp_kdes <- exp(lcl_exp_kdes)
    ktr_ab_c_bm <- exp(lktr_ab_c_bm + etalktr_ab_c_bm)
    kshed <- exp(lkshed + etalkshed)
    bl_cd3_bm <- exp(lbl_cd3_bm + etalbl_cd3_bm)
    bl_mbcma_bm <- exp(lbl_mbcma_bm + etalbl_mbcma_bm)

    # Table S1 Cl1: clearance decays from the first dose, Cl1 + Clvar * exp(-delta * t).
    # max(0, tafd()) is 0 before the first dose, when every drug state is zero.
    tfirst <- max(0, tafd())
    cl <- cl_exp_inf + cl_exp_component * exp(-cl_exp_kdes * tfirst)
    kel <- cl / vc # k10
    k12 <- q / vc
    k21 <- q / vp

    # Synthesis rates are set from the baselines so each target pool starts at
    # balance (Data S1 text after Eq 13; Eq 1 for the thymus).
    ksyn_cd3_thymus <- bl_cd3_thymus * ktr_cd3_thymus_c
    ksyn_cd3_bm <- bl_cd3_bm * kint_cd3
    ksyn_mbcma <- bl_mbcma_bm * (kint_bcma + kshed)

    target_cd3_thymus(0) <- bl_cd3_thymus
    target_cd3_central(0) <- bl_cd3_central
    target(0) <- bl_sbcma
    target_bcma_bonemarrow(0) <- bl_mbcma_bm
    target_bonemarrow(0) <- bl_sbcma_bm
    target_cd3_bonemarrow(0) <- bl_cd3_bm

    # Plasma
    d/dt(depot) <- -ka * depot # Eq 4
    d/dt(central) <- ka * depot / vc - kel * central - k12 * central + k21 * vp / vc * peripheral1 -
      kon_cd3 * target_cd3_central * central + koff_cd3 * drug_cd3_central -
      kon_sbcma * target * central + koff_sbcma * complex -
      ktr_ab_c_bm * central + ktr_ab_bm_c * vbm / vc * bonemarrow # Eq 5
    d/dt(peripheral1) <- k12 * vc / vp * central - k21 * peripheral1 - kel * peripheral1 # Eq 6
    d/dt(target_cd3_thymus) <- ksyn_cd3_thymus - ktr_cd3_thymus_c * target_cd3_thymus # Eq 1
    d/dt(target_cd3_central) <- ktr_cd3_thymus_c * target_cd3_thymus - kint_cd3 * target_cd3_central +
      ktr_cd3_bm_c * vbm / vc * target_cd3_bonemarrow - ktr_cd3_c_bm * target_cd3_central -
      kon_cd3 * target_cd3_central * central + koff_cd3 * drug_cd3_central -
      kon_cd3 * target_cd3_central * complex + koff_cd3 * complex_cd3_sbcma # Eq 2
    d/dt(target) <- ktr_sbcma_bm_c * vbm / vc * target_bonemarrow - kdeg_sbcma * target -
      kon_sbcma * target * central + koff_sbcma * complex -
      kon_sbcma * target * drug_cd3_central + koff_sbcma * complex_cd3_sbcma -
      ktr_sbcma_c_bm * target # Eq 3
    d/dt(complex) <- kon_sbcma * target * central - koff_sbcma * complex -
      kon_cd3 * target_cd3_central * complex + koff_cd3 * complex_cd3_sbcma -
      kdeg_sbcma * complex # Eq 7
    d/dt(drug_cd3_central) <- kon_cd3 * target_cd3_central * central - koff_cd3 * drug_cd3_central -
      kon_sbcma * target * drug_cd3_central + koff_sbcma * complex_cd3_sbcma -
      kint_cd3 * drug_cd3_central # Eq 8
    d/dt(complex_cd3_sbcma) <- kon_sbcma * target * drug_cd3_central - koff_sbcma * complex_cd3_sbcma +
      kon_cd3 * target_cd3_central * complex - koff_cd3 * complex_cd3_sbcma -
      kint_cd3 * complex_cd3_sbcma - kdeg_sbcma * complex_cd3_sbcma # Eq 9

    # Bone marrow (tumour microenvironment)
    d/dt(target_bcma_bonemarrow) <- ksyn_mbcma - kshed * target_bcma_bonemarrow - kint_bcma * target_bcma_bonemarrow -
      kon_mbcma * target_bcma_bonemarrow * bonemarrow + koff_mbcma * drug_bcma_bonemarrow -
      kon_mbcma * target_bcma_bonemarrow * drug_cd3_bonemarrow + koff_mbcma * trimer # Eq 10
    d/dt(target_bonemarrow) <- kshed * target_bcma_bonemarrow - ktr_sbcma_bm_c * target_bonemarrow -
      kon_sbcma * target_bonemarrow * bonemarrow + koff_sbcma * complex_bonemarrow -
      kon_sbcma * target_bonemarrow * drug_cd3_bonemarrow + koff_sbcma * complex_cd3_sbcma_bonemarrow +
      ktr_sbcma_c_bm * target # Eq 11
    d/dt(target_cd3_bonemarrow) <- ksyn_cd3_bm - kint_cd3 * target_cd3_bonemarrow -
      ktr_cd3_bm_c * target_cd3_bonemarrow + ktr_cd3_c_bm * vc / vbm * target_cd3_central -
      kon_cd3 * target_cd3_bonemarrow * bonemarrow + koff_cd3 * drug_cd3_bonemarrow -
      kon_cd3 * target_cd3_bonemarrow * complex_bonemarrow + koff_cd3 * complex_cd3_sbcma_bonemarrow -
      kon_cd3 * target_cd3_bonemarrow * drug_bcma_bonemarrow + koff_cd3 * trimer # Eq 12
    d/dt(bonemarrow) <- vc / vbm * ktr_ab_c_bm * central - ktr_ab_bm_c * bonemarrow -
      kon_mbcma * target_bcma_bonemarrow * bonemarrow + koff_mbcma * drug_bcma_bonemarrow -
      kon_sbcma * target_bonemarrow * bonemarrow + koff_sbcma * complex_bonemarrow -
      kon_cd3 * target_cd3_bonemarrow * bonemarrow + koff_cd3 * drug_cd3_bonemarrow # Eq 13
    d/dt(drug_cd3_bonemarrow) <- kon_cd3 * target_cd3_bonemarrow * bonemarrow - koff_cd3 * drug_cd3_bonemarrow -
      kon_sbcma * target_bonemarrow * drug_cd3_bonemarrow + koff_sbcma * complex_cd3_sbcma_bonemarrow -
      kon_mbcma * target_bcma_bonemarrow * drug_cd3_bonemarrow + koff_mbcma * trimer -
      kint_bcma * drug_cd3_bonemarrow # Eq 14
    d/dt(complex_bonemarrow) <- kon_sbcma * target_bonemarrow * bonemarrow - koff_sbcma * complex_bonemarrow -
      kon_cd3 * target_cd3_bonemarrow * complex_bonemarrow + koff_cd3 * complex_cd3_sbcma_bonemarrow # Eq 15
    d/dt(drug_bcma_bonemarrow) <- kon_mbcma * target_bcma_bonemarrow * bonemarrow - koff_mbcma * drug_bcma_bonemarrow -
      kon_cd3 * target_cd3_bonemarrow * drug_bcma_bonemarrow + koff_cd3 * trimer -
      kint_bcma * drug_bcma_bonemarrow # Eq 16
    d/dt(complex_cd3_sbcma_bonemarrow) <- kon_cd3 * target_cd3_bonemarrow * complex_bonemarrow -
      koff_cd3 * complex_cd3_sbcma_bonemarrow +
      kon_sbcma * target_bonemarrow * drug_cd3_bonemarrow - koff_sbcma * complex_cd3_sbcma_bonemarrow -
      kel_complex_bm * complex_cd3_sbcma_bonemarrow # Eq 17
    d/dt(trimer) <- kon_cd3 * target_cd3_bonemarrow * drug_bcma_bonemarrow - koff_cd3 * trimer +
      kon_mbcma * target_bcma_bonemarrow * drug_cd3_bonemarrow - koff_mbcma * trimer -
      kint_bcma * trimer # Eq 18

    # Free drug in plasma (nM); Figure 2A plots this times 146 kDa in ng/mL.
    Cc <- central
    Cc ~ prop(propSd)
  })
}
