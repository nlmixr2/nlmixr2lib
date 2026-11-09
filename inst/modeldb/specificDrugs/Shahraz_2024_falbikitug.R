Shahraz_2024_falbikitug <- function() {
  description <- "QSP. Mechanistic PK/PD model of falbikitug (AZD0171), a humanized IgG1 anti-leukemia inhibitory factor (LIF) antibody, in adults with advanced solid tumors and metastatic pancreatic ductal adenocarcinoma (Shahraz 2024; Phase I NCT03490669 n = 41 and interim Phase II NCT04999969 n = 49). Three compartments: central (serum), a peripheral tissue compartment and a tumor compartment. AZD0171 has linear clearance from the central compartment, distributes to tissue by an intercompartmental clearance and to the tumor by a Krogh-cylinder permeability-diffusion flux. LIF is synthesized only in the tumor at a zero-order rate, eliminated from the central compartment with a 2 h half-life, and exchanges with tissue and tumor. In tissue and tumor LIF binds LIFR and gp130 independently to form the LIF:LIFR:gp130 signaling complex; AZD0171 binds free LIF and LIF:LIFR (in the central compartment free LIF only) and blocks gp130 recruitment. LIFR and gp130 turn over to fitted / fixed steady-state levels, and occupied receptors internalize. All states start receptor-only with no LIF: like the authors (first dose at 240 h in Figure 4), run the drug-free system for at least 240 h before the first dose so that LIF and the receptor complexes reach baseline. Outputs: serum AZD0171 (Cc, ug/mL), serum total LIF (Ctotal_target, pg/mL), and tumor free LIF and signaling complex for the target-engagement predictions."
  reference <- paste(
    "Shahraz A, Penney M, Candido J, Opoku-Ansah G, Neubauer M, Eyles J,",
    "Ojo O, Liu N, Luheshi NM, Phipps A, Vishwanathan K. A mechanistic",
    "PK/PD model of AZD0171 (anti-LIF) to support Phase II dose selection.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(10):1670-1680.",
    "doi:10.1002/psp4.13204.",
    "Parameter values from Table 1 and Supplementary Table S4; equations",
    "from Supplementary Table S5 (fluxes, ODEs and rules) and Appendix S3",
    "(Krogh tumor-penetration flux).",
    sep = " "
  )
  vignette <- "Shahraz_2024_falbikitug"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ug/mL (AZD0171: Cc); pg/mL (total LIF: Ctotal_target)"
  )

  # Tissue and tumor species are per-compartment concentrations (nM), as in
  # the source's SimBiology model; they have no canonical compartment names.
  paper_specific_compartments <- c(
    "drug_tissue",
    "lif_tissue",
    "lifr_tissue",
    "gp130_tissue",
    "lif_lifr_tissue",
    "lif_gp130_tissue",
    "drug_lif_tissue",
    "drug_lif_lifr_tissue",
    "signal_tissue",
    "drug_tumor",
    "lif_tumor",
    "lifr_tumor",
    "gp130_tumor",
    "lif_lifr_tumor",
    "lif_gp130_tumor",
    "drug_lif_tumor",
    "drug_lif_lifr_tumor",
    "signal_tumor"
  )

  compartmentData <- list(
    central = list(analyte = "falbikitug", units = "mg", specimen = "serum", verified = TRUE),
    target = list(analyte = "LIF", units = "nmol", specimen = "serum", verified = TRUE),
    complex = list(analyte = "falbikitug:LIF complex", units = "nmol", specimen = "serum", verified = TRUE),
    drug_tissue = list(analyte = "falbikitug", units = "nM", specimen = "tissue", verified = TRUE),
    lif_tissue = list(analyte = "LIF", units = "nM", specimen = "tissue", verified = TRUE),
    lifr_tissue = list(analyte = "LIFR", units = "nM", specimen = "tissue", verified = TRUE),
    gp130_tissue = list(analyte = "gp130", units = "nM", specimen = "tissue", verified = TRUE),
    lif_lifr_tissue = list(analyte = "LIF:LIFR complex", units = "nM", specimen = "tissue", verified = TRUE),
    lif_gp130_tissue = list(analyte = "LIF:gp130 complex", units = "nM", specimen = "tissue", verified = TRUE),
    drug_lif_tissue = list(analyte = "falbikitug:LIF complex", units = "nM", specimen = "tissue", verified = TRUE),
    drug_lif_lifr_tissue = list(
      analyte = "falbikitug:LIF:LIFR complex",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    signal_tissue = list(
      analyte = "LIF:LIFR:gp130 signaling complex",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    drug_tumor = list(analyte = "falbikitug", units = "nM", specimen = "tumor", verified = TRUE),
    lif_tumor = list(analyte = "LIF", units = "nM", specimen = "tumor", verified = TRUE),
    lifr_tumor = list(analyte = "LIFR", units = "nM", specimen = "tumor", verified = TRUE),
    gp130_tumor = list(analyte = "gp130", units = "nM", specimen = "tumor", verified = TRUE),
    lif_lifr_tumor = list(analyte = "LIF:LIFR complex", units = "nM", specimen = "tumor", verified = TRUE),
    lif_gp130_tumor = list(analyte = "LIF:gp130 complex", units = "nM", specimen = "tumor", verified = TRUE),
    drug_lif_tumor = list(analyte = "falbikitug:LIF complex", units = "nM", specimen = "tumor", verified = TRUE),
    drug_lif_lifr_tumor = list(
      analyte = "falbikitug:LIF:LIFR complex",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    signal_tumor = list(analyte = "LIF:LIFR:gp130 signaling complex", units = "nM", specimen = "tumor", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 90,
    n_studies = 2,
    age_range = "18 years or older",
    disease_state = paste(
      "Advanced, unresectable solid tumors without curative options (Phase I",
      "NCT03490669, n = 41); first-line metastatic pancreatic ductal",
      "adenocarcinoma with tumoral CD8+ T cells, treated in combination with",
      "durvalumab, gemcitabine and nab-paclitaxel (interim Phase II",
      "NCT04999969, n = 49)."
    ),
    dose_range = "75-1500 mg IV every 3 weeks (Phase I dose escalation); 1500 mg IV every 2 weeks (Phase II).",
    regions = "Multinational (Phase II at 40-50 sites globally).",
    notes = paste(
      "The paper reports no demographic table. CL, Vc, ksynth_LIFt, LIFR_ss",
      "and Q were fitted per patient to serum AZD0171 and total LIF in",
      "SimBiology (MATLAB 2020b); the typical values are the medians of the",
      "individual fits (Table 1). Sampling schedules are in Appendix S1",
      "(Supplementary Tables S1-S3)."
    )
  )

  ini({
    # --- Fitted AZD0171 PK and LIF / LIFR parameters (medians of individual fits) ---
    lcl <- log(0.0135); label("AZD0171 clearance from the central compartment (L/h)") # Table 1 'AZD0171 clearance' 0.0135 [0.00877-0.0173]; Table S4 CL_AZD0171
    lvc <- log(3.72); label("Volume of the central compartment (L)") # Table 1 'Volume of the central compartment' 3.72 [3.05-4.69]; Table S4 Vcentral
    lq <- log(0.110); label("Intercompartmental clearance of AZD0171 and AZD0171:LIF, central to tissue (L/h)") # Table 1 'Intercompartmental clearance' 0.110 [0.0455-0.487]; Table S4 Q_drug
    lksyn_lif <- log(1.67); label("Zero-order LIF synthesis rate in the tumor (pmol/h)") # Table 1 'Synthesis Rate of LIF in tumor' 1.67 [0.290-8.22]; Table S4 ksynth_LIFt
    llifr_ss <- log(5.926); label("Steady-state LIFR concentration in tissue and tumor (nM)") # Table 1 'Steady-state LIFR concentration in tumor' 5.926 [1.49-194.5]; Table S4 LIFR_ss

    # --- Fixed system parameters (Table S4) ---
    lvp <- fixed(log(3.72)); label("Volume of the tissue compartment, set equal to Vc (L)") # Table S4 Vtissue 3.72, 'Calculated (Vcentral=Vtissue)'; Table S5 rule 21
    v_tumor <- fixed(0.258); label("Volume of the tumor compartment (L)") # Table S4 Vtumor 0.258, from the median baseline tumor size; Methods
    kon_drug_lif <- fixed(3.6); label("Association rate constant, AZD0171 to LIF (1/(nM*h))") # Table S4 konAZD0171LIF 3.6, 'Fixed value'
    kd_drug_lif <- fixed(0.074); label("Equilibrium dissociation constant, AZD0171 to LIF (nM)") # Table S4 KD_LIFAZD0171 0.074, internal data
    kon_lif_lifr <- fixed(3.6); label("Association rate constant, LIF to LIFR (1/(nM*h))") # Table S4 konLIFLIFR 3.6, 'Fixed value'
    koff_lif_lifr <- fixed(11.16); label("Dissociation rate constant, LIF from LIFR (1/h)") # Table S4 koffLIFLIFR 11.16 (= KD_LIFLIFR x konLIFLIFR)
    kon_lif_gp130 <- fixed(3.6); label("Association rate constant, LIF to gp130 (1/(nM*h))") # Table S4 konLIFgp130 3.6, 'Fixed value'
    koff_lif_gp130 <- fixed(260.568); label("Dissociation rate constant, LIF from gp130 (1/h)") # Table S4 koffLIFgp130 260.568 (= KD_LIFgp130 x konLIFgp130)
    thalf_lif <- fixed(2); label("Half-life of LIF degradation in the central compartment (h)") # Table S4 thalf_LIF 2, Gunawardana 2003
    gp130_ss <- fixed(25); label("Steady-state gp130 concentration in tissue and tumor (nM)") # Table S4 gp130_ss 25, same as mouse models
    kint_gp130 <- fixed(4.158); label("Internalization rate of gp130 and LIF:gp130 (1/h)") # Table S4 kintgp130 4.158
    kint_lifr_slow <- fixed(0.42); label("Internalization rate of LIFR, LIF:LIFR and AZD0171:LIF:LIFR (1/h)") # Table S4 kintLIFRslow 0.42, Hilton and Nicola 1992
    kint_lifr_fast <- fixed(1.89); label("Internalization rate of the LIF:LIFR:gp130 complex (1/h)") # Table S4 kintLIFRfast 1.89, Hilton and Nicola 1992
    r_cap <- fixed(8); label("Radius of tumor blood capillary (um)") # Table S4 Rcap 8
    r_krogh <- fixed(75); label("Average distance between two tumor capillaries (um)") # Table S4 RKrogh 75
    void_frac <- fixed(0.24); label("Tumor void fraction (unitless)") # Table S4 epsilon_tumor 0.24, Shah 2012
    mw_drug <- fixed(150); label("Molecular weight of AZD0171 (kDa)") # Table S4 MW_Drug 150
    mw_lif <- fixed(25); label("Molecular weight of LIF (kDa)") # Table S4 MW_LIF 25
    mw_drug_lif <- fixed(175); label("Molecular weight of the AZD0171:LIF complex (kDa)") # Table S4 MW_DrugLIF 175
    gamma_qlif <- fixed(2.25); label("Coefficient of the Li 2017 size-clearance relation for LIF intercompartmental clearance (L/h)") # Table S4 gamma 2.25, Li 2017
    mu_qlif <- fixed(2.48); label("Exponent of the Li 2017 size-clearance relation for LIF intercompartmental clearance (unitless)") # Table S4 mu 2.48, Li 2017

    # --- Between-patient variability for the virtual-patient simulations ---
    # Results: the five fitted parameters were 'taken from a log-normal
    # distribution with 30% coefficient of variation (CV)', so each variance is
    # log(1 + 0.30^2) = 0.0862. Assumed, not estimated.
    etalcl ~ fixed(0.0862) # Results, virtual patients: 30% CV
    etalvc ~ fixed(0.0862) # Results, virtual patients: 30% CV
    etalq ~ fixed(0.0862) # Results, virtual patients: 30% CV
    etalksyn_lif ~ fixed(0.0862) # Results, virtual patients: 30% CV
    etallifr_ss ~ fixed(0.0862) # Results, virtual patients: 30% CV

    # --- Residual error: not reported (individual SimBiology fits) ---
    propSd <- fixed(0); label("Proportional residual error on serum AZD0171 (fraction; not reported in the source)") # not reported
    propSd_Ctotal_target <- fixed(0); label("Proportional residual error on serum total LIF (fraction; not reported in the source)") # not reported
  })

  model({
    # Amounts in nmol, concentrations in nM, time in h (Table S5 units; the
    # pg/mL and ug/mL conversions are Table S5 rules 10, 11, 20 and 24).
    mg_to_nmol <- 1000 / mw_drug

    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp)
    lifr_ss <- exp(llifr_ss + etallifr_ss)
    # Table S4 gives ksynth_LIFt in pmol/h; the model runs in nmol (Figure 4d
    # legend lists the same rate as 0.0017 nmol/h).
    ksyn_lif <- exp(lksyn_lif + etalksyn_lif) / 1000

    koff_drug_lif <- kon_drug_lif * kd_drug_lif # Table S5 rule 23 (0.2664 1/h)
    kdeg_lif <- log(2) / thalf_lif # Table S5 Reaction_9

    # Li 2019 size relations (Table S5 rules 1-6). Diffusion in cm^2/h,
    # permeability in um/h. Rule 6 uses MW_Drug, not MW_DrugLIF, for the
    # AZD0171:LIF permeability, which Table S4 confirms (4.32 um/h).
    diff_lif <- 114.58 * mw_lif^(-0.918) * 3600 / 1e7
    diff_drug <- 114.58 * mw_drug^(-0.918) * 3600 / 1e7
    diff_drug_lif <- 114.58 * mw_drug_lif^(-0.918) * 3600 / 1e7
    pcap_lif <- (-10.43 * log(mw_lif) + 53.46) * 3600 / 1000
    pcap_drug <- (-10.43 * log(mw_drug) + 53.46) * 3600 / 1000
    pcap_drug_lif <- (-10.43 * log(mw_drug) + 53.46) * 3600 / 1000
    # Li 2017 relation, Table S5 rule 7 (0.271 L/h)
    q_lif <- gamma_qlif * mw_lif^(-mu_qlif) * (70 / 0.028)^0.75
    # Tumor radius from volume, Table S5 rule 8 (L -> cm^3)
    r_tumor <- (3 * v_tumor * 1000 / (4 * pi))^(1 / 3)

    # Krogh-cylinder exchange rate constants (1/h), Table S5 Reactions 13 and
    # 36-38, 44, 47: 2 * P * Rcap / RKrogh^2 + 6 * D / Rtumor^2
    k_krogh_drug <- 2 * pcap_drug * r_cap / r_krogh^2 + 6 * diff_drug / r_tumor^2
    k_krogh_lif <- 2 * pcap_lif * r_cap / r_krogh^2 + 6 * diff_lif / r_tumor^2
    k_krogh_drug_lif <- 2 * pcap_drug_lif * r_cap / r_krogh^2 + 6 * diff_drug_lif / r_tumor^2

    # Central concentrations (nM), Table S5 rules 12-14
    c_drug <- central * mg_to_nmol / vc
    c_lif <- target / vc
    c_drug_lif <- complex / vc

    # --- Central fluxes (nmol/h) ---
    r_bind_c <- kon_drug_lif * c_drug * c_lif * vc - koff_drug_lif * complex # Reaction_6
    r_q_drug <- q * (drug_tissue - c_drug) # Reaction_3
    r_q_drug_lif <- q * (drug_lif_tissue - c_drug_lif) # Reaction_10
    r_q_lif <- q_lif * (lif_tissue - c_lif) # Reaction_43
    r_tu_drug <- k_krogh_drug * (c_drug - drug_tumor / void_frac) # Reaction_37 (x Vtumor = -Reaction_38)
    r_tu_lif <- k_krogh_lif * (c_lif - lif_tumor / void_frac) # Reaction_44 (x Vtumor = -Reaction_36)
    r_tu_drug_lif <- k_krogh_drug_lif * (c_drug_lif - drug_lif_tumor / void_frac) # Reaction_47 (x Vtumor = -Reaction_13)

    # --- Tissue binding fluxes (nM/h) ---
    r11 <- kon_drug_lif * drug_tissue * lif_tissue - koff_drug_lif * drug_lif_tissue
    r2 <- kon_drug_lif * drug_tissue * lif_lifr_tissue - koff_drug_lif * drug_lif_lifr_tissue
    r14 <- kon_lif_lifr * lif_tissue * lifr_tissue - koff_lif_lifr * lif_lifr_tissue
    r15 <- kon_lif_gp130 * lif_tissue * gp130_tissue - koff_lif_gp130 * lif_gp130_tissue
    r16 <- kon_lif_gp130 * lif_lifr_tissue * gp130_tissue - koff_lif_gp130 * signal_tissue
    r17 <- kon_lif_lifr * lifr_tissue * drug_lif_tissue - koff_lif_lifr * drug_lif_lifr_tissue
    r19 <- kon_lif_lifr * lifr_tissue * lif_gp130_tissue - koff_lif_lifr * signal_tissue

    # --- Tumor binding fluxes (nM/h) ---
    r20 <- kon_drug_lif * drug_tumor * lif_tumor - koff_drug_lif * drug_lif_tumor
    r1 <- kon_drug_lif * drug_tumor * lif_lifr_tumor - koff_drug_lif * drug_lif_lifr_tumor
    r21 <- kon_lif_lifr * lif_tumor * lifr_tumor - koff_lif_lifr * lif_lifr_tumor
    r22 <- kon_lif_gp130 * gp130_tumor * lif_tumor - koff_lif_gp130 * lif_gp130_tumor
    r23 <- kon_lif_lifr * lifr_tumor * lif_gp130_tumor - koff_lif_lifr * signal_tumor
    r24 <- kon_lif_gp130 * gp130_tumor * lif_lifr_tumor - koff_lif_gp130 * signal_tumor
    r25 <- kon_lif_lifr * drug_lif_tumor * lifr_tumor - koff_lif_lifr * drug_lif_lifr_tumor

    # --- Central compartment (amounts) ---
    # Table S5 ODE 10 for AZD0171, in mg (divided by mg_to_nmol). The free LIF
    # and AZD0171:LIF ODEs are not printed in Table S5; they are assembled from
    # its flux list (Reactions 6, 8, 9, 10, 13, 36, 43).
    d/dt(central) <- (r_q_drug - cl * c_drug - r_bind_c - r_tu_drug * v_tumor) / mg_to_nmol
    d/dt(target) <- -r_bind_c - kdeg_lif * target + r_q_lif - r_tu_lif * v_tumor
    d/dt(complex) <- r_bind_c - cl * c_drug_lif + r_q_drug_lif - r_tu_drug_lif * v_tumor

    # --- Tissue compartment (concentrations), Table S5 ODEs 1-9 ---
    # LIF synthesis in tissue (Reaction_40, ksynth_LIFt * alpha / Vtissue) is
    # omitted: alpha is not defined in Table S4 and the Methods state that LIF
    # is synthesized only in the tumor.
    d/dt(drug_tissue) <- -r11 - r_q_drug / vp - r2
    d/dt(lif_tissue) <- -r11 - r14 - r15 - r_q_lif / vp
    d/dt(lifr_tissue) <- -r14 - r17 - r19 + kint_lifr_slow * (lifr_ss - lifr_tissue)
    d/dt(gp130_tissue) <- -r15 - r16 + kint_gp130 * (gp130_ss - gp130_tissue)
    d/dt(lif_lifr_tissue) <- r14 - r16 - r2 - kint_lifr_slow * lif_lifr_tissue
    d/dt(lif_gp130_tissue) <- r15 - r19 - kint_gp130 * lif_gp130_tissue
    d/dt(drug_lif_tissue) <- r11 - r17 - r_q_drug_lif / vp
    d/dt(drug_lif_lifr_tissue) <- r17 + r2 - kint_lifr_slow * drug_lif_lifr_tissue
    d/dt(signal_tissue) <- r16 + r19 - kint_lifr_fast * signal_tissue

    # --- Tumor compartment (concentrations) ---
    # Not printed in Table S5; assembled from its flux list (Reactions 1, 20-25,
    # 32-35, 37, 39, 44-47) in the same form as the printed tissue ODEs.
    # LIF synthesis enters per void volume (Reaction_39).
    d/dt(drug_tumor) <- r_tu_drug - r20 - r1
    d/dt(lif_tumor) <- r_tu_lif + ksyn_lif / (void_frac * v_tumor) - r20 - r21 - r22
    d/dt(lifr_tumor) <- kint_lifr_slow * (lifr_ss - lifr_tumor) - r21 - r23 - r25
    d/dt(gp130_tumor) <- kint_gp130 * (gp130_ss - gp130_tumor) - r22 - r24
    d/dt(lif_lifr_tumor) <- r21 - r24 - r1 - kint_lifr_slow * lif_lifr_tumor
    d/dt(lif_gp130_tumor) <- r22 - r23 - kint_gp130 * lif_gp130_tumor
    d/dt(drug_lif_tumor) <- r20 + r_tu_drug_lif - r25
    d/dt(drug_lif_lifr_tumor) <- r25 + r1 - kint_lifr_slow * drug_lif_lifr_tumor
    d/dt(signal_tumor) <- r23 + r24 - kint_lifr_fast * signal_tumor

    lifr_tissue(0) <- lifr_ss
    gp130_tissue(0) <- gp130_ss
    lifr_tumor(0) <- lifr_ss
    gp130_tumor(0) <- gp130_ss

    # --- Derived outputs (Table S5 rules) ---
    lif_free_serum_pgml <- c_lif * mw_lif * 1000 # rule 20
    lif_free_tumor_pgml <- lif_tumor * mw_lif * 1000 # rule 24
    lifr_total_tumor <- lifr_tumor + lif_lifr_tumor + drug_lif_lifr_tumor + signal_tumor # rule 18
    active_lifr_pct_tumor <- signal_tumor * 100 / lifr_total_tumor # rule 19

    # Serum AZD0171 (ug/mL; rule 11 equals central / vc in mg/L)
    Cc <- central / vc
    # Serum total LIF (pg/mL), rule 10 as printed: the complex is weighted by
    # MW_DrugLIF (175 kDa), not by the LIF mass alone.
    Ctotal_target <- (c_lif * mw_lif + c_drug_lif * mw_drug_lif) * 1000

    Cc ~ prop(propSd)
    Ctotal_target ~ prop(propSd_Ctotal_target)
  })
}
