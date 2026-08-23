Patidar_2024_mAb_mouse <- function() {
  description <- "PBPK (minimal, mPBPK; 39 ODEs). Preclinical (mouse). Platform minimal-PBPK model relating monoclonal-antibody physicochemical properties (molecular weight, Stoke's radius, net surface charge, FcRn affinity, antigen affinity) to plasma and tissue disposition. Plasma and lymph compartments plus lumped tight-tissue and leaky-tissue compartments, each split into vascular, endosomal and interstitial sub-spaces; plasma carries a nested endothelial endosome. Includes explicit FcRn binding and recycling at endosomal pH 6, two-pore size-based transcapillary transport, size-based renal clearance, charge-dependent pinocytosis / non-specific membrane binding / interstitial volume, and TMDD against both soluble and membrane-bound antigen. Defaults reproduce Case 1: non-specific FcRn-binding IgG1 (150 kDa, neutral charge, no target) in wild-type mice."
  reference   <- "Patidar K, Pillai N, Dhakal S, Avery LB, Mavroudis PD. A minimal physiologically based pharmacokinetic model to study the combined effect of antibody size, charge, and binding affinity to FcRn/antigen on antibody pharmacokinetics. J Pharmacokinet Pharmacodyn. 2024;51(6):477-492. doi:10.1007/s10928-023-09899-z. PMID 38386198. PMCID PMC11576895. Model equations A1-A48, physiology/kinetic parameter Table 1 (mouse column) and two-pore Table 2 are in the Supplementary Information (10928_2023_9899_MOESM1_ESM.docx). Two-pore permeability-surface-area (PS) and Peclet-number (Pe) derivations are inherited from the cited upstream framework: Li Z, Shah DK. Two-pore physiologically based pharmacokinetic model with de novo derived parameters for predicting plasma PK of different size protein therapeutics. J Pharmacokinet Pharmacodyn. 2019;46(3):305-318. doi:10.1007/s10928-019-09639-2 (Supplementary Material Eqs. 13-30; on disk as LiShah_2019_twopore_10928_2019_9639_MOESM1_ESM.docx)."
  vignette    <- "Patidar_2024_mAb_size_charge_fcrn"

  # All 39 states are paper-mechanistic sub-space / binding species of the
  # Patidar 2024 mPBPK system and do not map onto canonical compartments.
  # Suffixes: _p plasma, _e plasma nested endosome, _v1/_v2 tight/leaky
  # vascular, _e1/_e2 tight/leaky endosome, _is1/_is2 tight/leaky
  # interstitium, _l lymph.
  paper_specific_compartments <- c(
    "a_p", "t_p", "atc_p",
    "a_e", "t_e", "atc_e", "fcrn_e", "fcrna_e", "fcrnatc_e",
    "a_v1", "atc_v1", "t_v1", "tm_v1", "arm_v1", "atm_v1",
    "a_e1", "t_e1", "atc_e1", "fcrn_e1", "fcrna_e1", "fcrnatc_e1",
    "a_is1", "atc_is1",
    "a_v2", "atc_v2", "t_v2", "tm_v2", "arm_v2", "atm_v2",
    "a_e2", "t_e2", "atc_e2", "fcrn_e2", "fcrna_e2", "fcrnatc_e2",
    "a_is2", "atc_is2",
    "a_l", "atc_l"
  )

  units <- list(time = "h", dosing = "nmol", concentration = "nM")

  # Issue #482: every state holds a CONCENTRATION (nM) in its sub-space, as
  # written in Patidar 2024 eqs. A1-A39 (the published system is
  # concentration-based, not amount-based). verified = TRUE means the analyte
  # and specimen were read off the supplement's state-variable table
  # (MOESM1 "State variable / Description" table, p. 5).
  compartmentData <- list(
    a_p        = list(analyte = "free antibody", units = "nM", specimen = "plasma", verified = TRUE),
    t_p        = list(analyte = "free soluble antigen", units = "nM", specimen = "plasma", verified = TRUE),
    atc_p      = list(analyte = "antibody-soluble-antigen complex", units = "nM", specimen = "plasma", verified = TRUE),
    a_e        = list(analyte = "free antibody", units = "nM", specimen = "endosome", verified = TRUE),
    t_e        = list(analyte = "free soluble antigen", units = "nM", specimen = "endosome", verified = TRUE),
    atc_e      = list(analyte = "antibody-soluble-antigen complex", units = "nM", specimen = "endosome", verified = TRUE),
    fcrn_e     = list(analyte = "free FcRn", units = "nM", specimen = "endosome", verified = TRUE),
    fcrna_e    = list(analyte = "FcRn-antibody complex", units = "nM", specimen = "endosome", verified = TRUE),
    fcrnatc_e  = list(analyte = "FcRn-antibody-antigen complex", units = "nM", specimen = "endosome", verified = TRUE),
    a_v1       = list(analyte = "free antibody", units = "nM", specimen = "tissue", verified = TRUE),
    atc_v1     = list(analyte = "antibody-soluble-antigen complex", units = "nM", specimen = "tissue", verified = TRUE),
    t_v1       = list(analyte = "free soluble antigen", units = "nM", specimen = "tissue", verified = TRUE),
    tm_v1      = list(analyte = "free membrane-bound antigen", units = "nM", specimen = "tissue", verified = TRUE),
    arm_v1     = list(analyte = "antibody-membrane-protein complex", units = "nM", specimen = "tissue", verified = TRUE),
    atm_v1     = list(analyte = "antibody-membrane-antigen complex", units = "nM", specimen = "tissue", verified = TRUE),
    a_e1       = list(analyte = "free antibody", units = "nM", specimen = "endosome", verified = TRUE),
    t_e1       = list(analyte = "free soluble antigen", units = "nM", specimen = "endosome", verified = TRUE),
    atc_e1     = list(analyte = "antibody-soluble-antigen complex", units = "nM", specimen = "endosome", verified = TRUE),
    fcrn_e1    = list(analyte = "free FcRn", units = "nM", specimen = "endosome", verified = TRUE),
    fcrna_e1   = list(analyte = "FcRn-antibody complex", units = "nM", specimen = "endosome", verified = TRUE),
    fcrnatc_e1 = list(analyte = "FcRn-antibody-antigen complex", units = "nM", specimen = "endosome", verified = TRUE),
    a_is1      = list(analyte = "free antibody", units = "nM", specimen = "tissue", verified = TRUE),
    atc_is1    = list(analyte = "antibody-soluble-antigen complex", units = "nM", specimen = "tissue", verified = TRUE),
    a_v2       = list(analyte = "free antibody", units = "nM", specimen = "tissue", verified = TRUE),
    atc_v2     = list(analyte = "antibody-soluble-antigen complex", units = "nM", specimen = "tissue", verified = TRUE),
    t_v2       = list(analyte = "free soluble antigen", units = "nM", specimen = "tissue", verified = TRUE),
    tm_v2      = list(analyte = "free membrane-bound antigen", units = "nM", specimen = "tissue", verified = TRUE),
    arm_v2     = list(analyte = "antibody-membrane-protein complex", units = "nM", specimen = "tissue", verified = TRUE),
    atm_v2     = list(analyte = "antibody-membrane-antigen complex", units = "nM", specimen = "tissue", verified = TRUE),
    a_e2       = list(analyte = "free antibody", units = "nM", specimen = "endosome", verified = TRUE),
    t_e2       = list(analyte = "free soluble antigen", units = "nM", specimen = "endosome", verified = TRUE),
    atc_e2     = list(analyte = "antibody-soluble-antigen complex", units = "nM", specimen = "endosome", verified = TRUE),
    fcrn_e2    = list(analyte = "free FcRn", units = "nM", specimen = "endosome", verified = TRUE),
    fcrna_e2   = list(analyte = "FcRn-antibody complex", units = "nM", specimen = "endosome", verified = TRUE),
    fcrnatc_e2 = list(analyte = "FcRn-antibody-antigen complex", units = "nM", specimen = "endosome", verified = TRUE),
    a_is2      = list(analyte = "free antibody", units = "nM", specimen = "tissue", verified = TRUE),
    atc_is2    = list(analyte = "antibody-soluble-antigen complex", units = "nM", specimen = "tissue", verified = TRUE),
    a_l        = list(analyte = "free antibody", units = "nM", specimen = "lymph", verified = TRUE),
    atc_l      = list(analyte = "antibody-soluble-antigen complex", units = "nM", specimen = "lymph", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "mouse (wild-type and FcRn knockout)",
    n_subjects     = NA,
    n_studies      = 6,
    age_range      = NA,
    weight_range   = "28 g reference body weight (Patidar 2024 Supplementary Table 1, mouse column)",
    sex_female_pct = NA,
    race_ethnicity = NA,
    disease_state  = "Healthy wild-type mice, FcRn knockout mice, and human-FcRn transgenic mice; plus tumour-bearing mice for the anti-CEA target-mediated case.",
    dose_range     = "IV bolus. Non-specific IgG1 8 mg/kg (fitting, Patidar Fig. 2) and 3.8 ug/kg (validation, Fig. 3); size variants 5 mg/kg (Fig. 4); charge variants 10 mg/kg (Fig. 5) and 5 mg/kg (Fig. A6); anti-CEA IgG 1, 10 and 25 mg/kg (Fig. 6).",
    regions        = NA,
    notes          = paste0(
      "Typical-value mechanistic (platform mPBPK) simulator: the paper reports no ",
      "between-subject variability and no residual-error model, so this file carries ",
      "neither. Observed data were digitized from published sources with ",
      "WebPlotDigitizer and summarised as min-max intervals. Only the most sensitive ",
      "parameters were estimated (kup, kup_p, sigma_tight, sigma_leaky by local ",
      "sensitivity analysis; Rm_total, KD_NSB, Spino, Kp for charge variants; ",
      "kptm and kint for the anti-CEA case); all other parameters are a priori. ",
      "Defaults here are Case 1 (non-specific FcRn-binding IgG1 in wild-type mice, ",
      "150 kDa, neutral charge, target expression set to zero). Alternative ",
      "published configurations are given in the ini() comments and worked in the ",
      "validation vignette: FcRn knockout (kup_per_h 0.15, or 0.26 for the human-FcRn ",
      "transgenic dataset); size variants (mw 100000 or 50000); charge variants ",
      "(charge -8, 0, +5, with spino2 2.99 for charge +5); anti-CEA TMDD ",
      "(icc_t 0.016 nM, icc_tm 80 nM)."
    ),
    model_class    = "Minimal PBPK / QSP platform model (two-pore size-based transport, explicit FcRn recycling in three endosomal sub-spaces, charge-dependent pinocytosis and non-specific binding, soluble + membrane-bound TMDD)",
    n_states       = 39,
    upstream_model = paste0(
      "Structure follows Cao/Jusko second-generation mPBPK lumping (tight vs leaky ",
      "tissues), the Yuan 2018 nested-plasma-endosome mPBPK, and the Shah & Betts 2012 ",
      "platform mAb PBPK. Two-pore PS and Peclet derivations are inherited from Li & ",
      "Shah 2019 (doi:10.1007/s10928-019-09639-2) Supplementary Material Eqs. 13-30, ",
      "which is on disk in the source directory."
    )
  )

  ini({
    # ---- Compound-specific (physicochemical) inputs ------------------------
    # Patidar 2024 "Model Parameters": the compound-specific model input is
    # molecular weight, Stoke's radius (derived from MW), net surface charge,
    # FcRn affinity and antigen affinity.
    mw     <- fixed(150000); label("Antibody molecular weight (Da)")                                                   # Patidar 2024 Results "PK of size-variants": full-length IgG 150 kDa; 100 kDa one-armed IgG and 50 kDa Fc-less fragment are the other published cases
    charge <- fixed(0);      label("Antibody net surface charge (unitless)")                                           # Patidar 2024 Results "PK of charge-variants": neutral variant has zero net surface charge; published variants are -8, 0 and +5 (validation set -10, -4, +10)

    # ---- Mouse physiology (Patidar 2024 Supplementary Table 1, mouse column)
    bw     <- fixed(0.028);    label("Body weight (kg)")                                                              # Patidar 2024 Suppl. Table 1 (Shah 2012)
    lvp    <- fixed(log(8.5e-4)); label("Plasma volume (L, log-scale)") # Patidar 2024 Suppl. Table 1: 0.85 mL (Shah 2012)
    vlymph <- fixed(1.717e-3); label("Lymph compartment volume (L)")                                                   # Patidar 2024 Suppl. Table 1: 1.717 mL (Shah 2012)
    l1     <- fixed(2.7820e-4); label("Lymph flow rate in tight tissue (L/h)")                                         # Patidar 2024 Suppl. Table 1 (Cao 2013)
    l2     <- fixed(4.6792e-4); label("Lymph flow rate in leaky tissue (L/h)")                                         # Patidar 2024 Suppl. Table 1 (Cao 2013)
    ltot   <- fixed(0.000738);  label("Total lymph flow rate (L/h)")                                                   # Patidar 2024 Suppl. Table 1 (Yuan 2018); the table note defines L = L1 + L2, which gives 7.4612e-4 (1.1% higher than the tabulated value used here)
    ve1    <- fixed(9.404e-5); label("Endosomal volume in tight tissue (L)")                                           # Patidar 2024 Suppl. Table 1: 0.09404 mL (Shah 2012)
    ve2    <- fixed(3.685e-5); label("Endosomal volume in leaky tissue (L)")                                           # Patidar 2024 Suppl. Table 1: 0.03685 mL (Shah 2012)
    vv1    <- fixed(4.695e-4); label("Vascular volume in tight tissue (L)")                                            # Patidar 2024 Suppl. Table 1: 0.4695 mL (Shah 2012)
    vv2    <- fixed(3.477e-4); label("Vascular volume in leaky tissue (L)")                                            # Patidar 2024 Suppl. Table 1: 0.3477 mL (Shah 2012)
    vis1   <- fixed(3.554e-3); label("Interstitial fluid volume in tight tissue (L)")                                   # Patidar 2024 Suppl. Table 1: 3.554 mL (Shah 2012)
    vis2   <- fixed(1.3539e-3); label("Interstitial fluid volume in leaky tissue (L)")                                  # Patidar 2024 Suppl. Table 1: 1.3539 mL (Shah 2012)
    sigma_tight <- fixed(0.9);  label("Vascular reflection coefficient, tight tissue (unitless)")                       # Patidar 2024 Suppl. Table 1, estimated (Results: "The estimated sigma1 and sigma2 values are 0.9 and 0.86")
    sigma_leaky <- fixed(0.86); label("Vascular reflection coefficient, leaky tissue (unitless)")                       # Patidar 2024 Suppl. Table 1, estimated (Results: "The estimated sigma1 and sigma2 values are 0.9 and 0.86")
    fr     <- fixed(0.715);    label("Fraction of endosomal recycling directed to vascular space (unitless)")           # Patidar 2024 Suppl. Table 1 (Shah 2012)
    fcrn_b <- fixed(40000);    label("Total FcRn concentration in endosomes (nM)")                                      # Patidar 2024 Suppl. Table 1: 40 uM (Shah 2012)
    gfr    <- fixed(0.0167);   label("Glomerular filtration rate (L/h)")                                                # Patidar 2024 Suppl. Table 1 (Shah 2012)

    # ---- Kinetic rate constants (Patidar 2024 Suppl. Table 1, mouse column)
    kup_p_per_h <- fixed(0.05);   label("Pinocytosis uptake rate constant into plasma nested endosome (1/h)")           # Patidar 2024 Suppl. Table 1, estimated; Results: "The estimated uptake rate for nested endosomes in plasma (kup,p) ... is 0.05 1/h"
    kup_per_h   <- fixed(0.0276); label("Pinocytosis uptake rate constant into tissue endosomes (1/h)")                 # Patidar 2024 Results: "... and tissue endosomes (kup) is ... 0.0276 1/h" (Suppl. Table 1 rounds this to 0.027); for FcRn knockout mice kup was recalibrated to 0.15 1/h, and to 0.26 1/h for the human-FcRn transgenic dataset (Fig. A8)
    clrec_p     <- fixed(7.27e-6); label("FcRn recycling clearance from the plasma nested endosome (L/h)")              # Patidar 2024 Suppl. Table 1 (Yuan 2018). Table note: "CLrec is calculated using transit time and volume of endosomes"; clrec_p/ve = 5.19 1/h implies an endosomal transit time of 11.5 min, and the same transit time is applied to the tissue endosomes in model()
    kdeg_per_h  <- fixed(42.9);   label("Lysosomal degradation rate constant of free antibody in endosomes (1/h)")      # Patidar 2024 Suppl. Table 1 (Shah 2012)
    clcat       <- fixed(5.3958e-5); label("Catabolic clearance of antigen and antibody-antigen complex in endosomes (L/h)") # Patidar 2024 Suppl. Table 1: 1.295 mL/day = 1.295e-3/24 L/h (Yuan 2018)
    k1on_per_nm_h <- fixed(0.12);  label("Antibody-FcRn association rate constant at endosomal pH 6 (1/(nM*h))")        # Patidar 2024 Suppl. Table 1 (Yuan 2018)
    k1off_per_h   <- fixed(1.8);   label("Antibody-FcRn dissociation rate constant at endosomal pH 6 (1/h)")            # Patidar 2024 Suppl. Table 1 (Yuan 2018)
    kon_per_nm_h  <- fixed(5.9);   label("Antibody-antigen association rate constant at pH 7.4 (1/(nM*h))")             # Patidar 2024 Suppl. Table 1 (Ferl 2005)
    koff_per_h    <- fixed(0.0508); label("Antibody-antigen dissociation rate constant at pH 7.4 (1/h)")                # Patidar 2024 Suppl. Table 1 (Ferl 2005)
    keon_per_nm_h <- fixed(5.904); label("Antibody-antigen association rate constant at endosomal pH 6 (1/(nM*h))")     # Patidar 2024 Suppl. Table 1 (Ferl 2005); the table prints the unit as "1/h" but an association constant is second-order and the human column gives "1/nM.h", so 1/(nM*h) is used
    keoff_per_h   <- fixed(0.0508); label("Antibody-antigen dissociation rate constant at endosomal pH 6 (1/h)")        # Patidar 2024 Suppl. Table 1 (Ferl 2005)
    kpt_per_h     <- fixed(0.3465); label("Degradation rate constant of soluble antigen (1/h)")                         # Patidar 2024 Suppl. Table 1, estimated; corresponds to the assumed 2 h soluble-antigen half-life and reproduces the tabulated ksyn (0.016 nM * 0.3465 1/h = 0.0055 nM/h)
    kptm_per_h    <- fixed(0.0192); label("Degradation rate constant of membrane-bound antigen (1/h)")                   # Patidar 2024 Suppl. Table 1, estimated; Results: membrane-bound antigen half-life estimated at 36 h and degradation rate 0.019 1/h; reproduces the tabulated ksyn,m (80 nM * 0.0192 1/h = 1.536 nM/h)
    kint_per_h    <- fixed(0.015); label("Internalization rate constant of membrane-bound antibody-antigen complex (1/h)") # Patidar 2024 Suppl. Table 1, estimated; Results: "The estimated internalization rate (kint) of membrane-bound mAb-antigen complex is 0.015 1/h"
    rm_total      <- fixed(71.86); label("Total negatively charged membrane protein concentration (nM)")                 # Patidar 2024 Suppl. Table 1, estimated; Results: "Rm,total and KD,NSB were estimated to be 71.86 nM and 8.35 nM"

    # ---- Antigen baseline expression --------------------------------------
    # Defaults are zero: Case 1-3 (non-specific IgG, size variants, charge
    # variants) carry no specific target, so "target-mediated specific
    # clearance was not included in the model" (Patidar 2024 Model fitting
    # and validation).
    icc_t  <- fixed(0); label("Baseline concentration of soluble antigen (nM)")                                          # Patidar 2024 Suppl. Table 1 gives 0.016 nM (Kiyama 1990) for the anti-CEA case; set to 0 here because the default configuration is non-specific IgG with no target
    icc_tm <- fixed(0); label("Baseline concentration of membrane-bound antigen (nM)")                                    # Patidar 2024 Suppl. Table 1 gives 80 nM (Ferl 2005) for the anti-CEA case; Results: "membrane-bound antigens to 80 nM". Set to 0 here for the non-specific IgG default

    # ---- Charge-dependent pinocytosis scaling factors ---------------------
    spino_p <- fixed(1); label("Pinocytosis scaling factor, plasma nested endosome (unitless)")                          # Patidar 2024 Model fitting: "We fixed Spino to 1 in nested endosomes in plasma, as charge did not affect IgG uptake in nested plasma endosomes"
    spino1  <- fixed(1); label("Pinocytosis scaling factor, tight tissue endosome (unitless)")                           # Patidar 2024 Model fitting: Spino,1 = 1 for neutral IgG, fixed to 1 for positively charged IgG, and assumed 1 for negatively charged IgG
    spino2  <- fixed(1); label("Pinocytosis scaling factor, leaky tissue endosome (unitless)")                           # Patidar 2024 Model fitting: Spino,2 = 1 for neutral and negatively charged IgG; estimated as 2.99 for the +5 charge variant (Results "PK of charge-variants")

    # ---- Two-pore transport constants (Patidar 2024 Suppl. Table 2) -------
    r_large <- fixed(22.85); label("Large pore radius (nm)")                                                             # Patidar 2024 Suppl. Table 2 (Li & Shah 2019)
    r_small <- fixed(4.44);  label("Small pore radius (nm)")                                                             # Patidar 2024 Suppl. Table 2 (Li & Shah 2019)
    alpha_l <- fixed(0.042); label("Fraction of hydraulic conductivity through large pores (unitless)")                   # Patidar 2024 Suppl. Table 2 (Li & Shah 2019)
    x_j     <- fixed(0.38);  label("Isogravimetric recirculation flow as a fraction of tissue lymph flow (unitless)")     # Patidar 2024 Suppl. Table 2: Jiso,tight = 0.38*L1 and Jiso,leaky = 0.38*L2 (Li & Shah 2019 Suppl. Eq. 14, where 0.042*0.958*(0.95-0.108)*10/1 = 0.38)
    xp_nm3  <- fixed(13197); label("Two-pore permeability constant Xp (nm^3)")                                            # Li & Shah 2019 Suppl. Eq. 20: Xp = R*T/(6*pi*N) * 8/(Starling force) = 13197 nm^3, with R = 62.363 L*mmHg/K/mol, T = 310 K, N = 6.02e23 1/mol, Starling force 1 mmHg
  })

  model({
    # =====================================================================
    # Derived compound properties
    # =====================================================================
    # Back-transform the log-scale plasma volume
    vp <- exp(lvp)

    # Endosomal volume of the systemic vascular endothelial cells. Patidar
    # 2024 Suppl. Table 1 gives this as a formula rather than a number:
    # Ve = 0.005/100 * BW (Yuan 2018), i.e. 0.005% of body weight.
    ve <- 0.005 / 100 * bw

    # Stoke's radius from molecular weight (Patidar 2024 eq. A44)
    ae <- 0.0483 * mw^0.386

    # Vascular reflection coefficients of large and small pores
    # (Patidar 2024 eqs. A45-A46 == Li & Shah 2019 Suppl. Eqs. 27-28)
    sigma_large <- 0.000035 * mw^0.717
    sigma_small <- 1 - 0.8489 * exp(-0.00004 * mw)

    # Fractional accessible cross-sectional pore area
    # (Li & Shah 2019 Suppl. Eqs. 25-26). The second equation is printed in
    # that supplement with the label "A_S/A_oS" repeated; it is the LARGE-pore
    # expression, identified by its much slower decay with MW (at 150 kDa it
    # gives 0.35, matching the paper's statement that a 150 kDa antibody is
    # "almost completely restricted by small pores" but retains "relatively
    # low permeability through large pores").
    a_small_frac <- 0.2352 * exp(-0.00008295 * mw) + 0.7767 * exp(-0.00053095 * mw)
    a_large_frac <- 0.3429 * exp(-0.00012175 * mw) + 0.6571 * exp(-0.00000421 * mw)

    # Size-independent groupings of the PS expressions
    # (Li & Shah 2019 Suppl. Eqs. 23-24); PS_i = xps * L_i (Eqs. 21-22)
    xps <- xp_nm3 * (1 / ae) * a_small_frac * (1 - alpha_l) / r_small^2
    xpl <- xp_nm3 * (1 / ae) * a_large_frac * alpha_l / r_large^2

    # Peclet numbers (Li & Shah 2019 Suppl. Eqs. 29-30). Tissue lymph flow
    # cancels, so a single Pe applies to both tissue compartments.
    pe_l <- (x_j + alpha_l) * (1 - sigma_large) / xpl
    pe_s <- (-x_j + (1 - alpha_l)) * (1 - sigma_small) / xps

    # Lymph flow through large and small pores (Patidar 2024 eqs. A40-A43).
    # The small-pore flows carry -Jiso, not +Jiso as printed in A41/A43:
    # the isogravimetric flow recirculates out through large pores and back
    # in through small pores, so J_large + J_small must equal the tissue
    # lymph flow. This is the form given by the cited upstream derivation
    # (Li & Shah 2019 Suppl. Eq. 30).
    jl1 <- (x_j + alpha_l) * l1
    js1 <- (-x_j + (1 - alpha_l)) * l1
    jl2 <- (x_j + alpha_l) * l2
    js2 <- (-x_j + (1 - alpha_l)) * l2

    # Permeability-surface area products (Li & Shah 2019 Suppl. Eqs. 21-22)
    ps_l1 <- xpl * l1
    ps_s1 <- xps * l1
    ps_l2 <- xpl * l2
    ps_s2 <- xps * l2

    # Size-based renal clearance (Patidar 2024 eq. 4) with the glomerular
    # sieving coefficient fitted against Stoke's radius in Fig. A2
    theta_sieve <- 0.9259 * exp(-((ae - 1.254) / 1.254)^2)
    cl_renal    <- gfr * theta_sieve

    # Charge-dependent parameters.
    # Kp, the fraction of interstitial volume available (Patidar 2024
    # Fig. A4): reproduces the fitted values 1 (neutral), 0.80 (+5) and
    # 0.62 (-8).
    kp <- -0.0067 * charge^2 - 0.0063 * charge + 1

    # Non-specific binding affinity to negatively charged membrane proteins
    # (Patidar 2024 Fig. A3). The figure plots 1/KD,NSB against net charge
    # and prints y = 1/(0.05906*exp(0.5105*x) + 0.06062); reading the printed
    # equation as KD,NSB reproduces all three fitted values -- 8.35 nM
    # (neutral), 1.21 nM (+5) and 16.22 nM (-8).
    kd_nsb    <- 1 / (0.05906 * exp(0.5105 * charge) + 0.06062)
    # kon,NSB is assumed equal to the antibody-FcRn association rate constant
    # (Patidar 2024 Suppl. Table 3 assumption 10)
    kon_nsb   <- k1on_per_nm_h
    koff_nsb  <- kon_nsb * kd_nsb

    # =====================================================================
    # Derived physiological groupings
    # =====================================================================
    # Endosomal recycling. The tabulated CLrec is the plasma nested-endosome
    # value; the table note states CLrec is computed from the endosomal
    # transit time and the endosomal volume, so the implied first-order
    # recycling rate constant clrec_p/ve applies to all three endosomes and
    # the per-compartment clearances are krec*ve_i.
    krec <- clrec_p / ve

    # Available interstitial volumes. Kp scales the interstitial space only
    # ("The change in volume of distribution for charge-variants is included
    # in the model by multiplying a scaling factor (Kp) to the volume of
    # interstitial space"); eqs. A10/A11/A24/A25 additionally print Kp on the
    # vascular volume, which is not what the text describes.
    vis1_eff <- kp * vis1
    vis2_eff <- kp * vis2

    # Antigen synthesis rates (Patidar 2024 Suppl. Table 1 note:
    # ksyn = ICC_p,T * k_p,T and ksyn,m = ICC_p,Tm * k_p,Tm)
    ksyn  <- icc_t * kpt_per_h
    ksynm <- icc_tm * kptm_per_h

    # =====================================================================
    # Two-pore transcapillary transport (Patidar 2024 eqs. A47-A48).
    # A47/A48 define a clearance CL_TP such that the vascular-to-interstitial
    # flux is CL_TP * C_vascular, carrying the concentration gradient in a
    # (1 - Cis/Cv) factor. Multiplying through by C_vascular gives the
    # equivalent flux form used here, which avoids dividing by a state that
    # is zero before the dose.
    # =====================================================================
    ftp_l1 <- ps_l1 * pe_l / (exp(pe_l) - 1)
    ftp_s1 <- ps_s1 * pe_s / (exp(pe_s) - 1)
    ftp_l2 <- ps_l2 * pe_l / (exp(pe_l) - 1)
    ftp_s2 <- ps_s2 * pe_s / (exp(pe_s) - 1)

    # Diffusive (gradient-driven) plus convective (solvent-drag) transport
    flux_tp1 <- (ftp_l1 + ftp_s1) * (a_v1 - a_is1) +
                (jl1 * (1 - sigma_large) + js1 * (1 - sigma_small)) * a_v1
    flux_tp2 <- (ftp_l2 + ftp_s2) * (a_v2 - a_is2) +
                (jl2 * (1 - sigma_large) + js2 * (1 - sigma_small)) * a_v2

    # =====================================================================
    # Binding fluxes (nM/h in the sub-space where they occur)
    # =====================================================================
    # Plasma: antibody-soluble antigen
    bind_p   <- kon_per_nm_h * a_p * t_p - koff_per_h * atc_p

    # Plasma nested endosome: antibody-antigen at pH 6, antibody-FcRn,
    # and antigen binding onto the FcRn-antibody complex
    bind_e   <- keon_per_nm_h * a_e * t_e - keoff_per_h * atc_e
    fcb_e    <- k1on_per_nm_h * a_e * fcrn_e - k1off_per_h * fcrna_e
    fcb_e_c  <- k1on_per_nm_h * atc_e * fcrn_e - k1off_per_h * fcrnatc_e
    bind_ef  <- keon_per_nm_h * fcrna_e * t_e - keoff_per_h * fcrnatc_e

    # Tight tissue vascular space
    bind_v1  <- kon_per_nm_h * a_v1 * t_v1 - koff_per_h * atc_v1
    bind_m1  <- kon_per_nm_h * a_v1 * tm_v1 - koff_per_h * atm_v1
    bind_r1  <- kon_nsb * a_v1 * rm_total - koff_nsb * arm_v1

    # Tight tissue endosome
    bind_e1  <- keon_per_nm_h * a_e1 * t_e1 - keoff_per_h * atc_e1
    fcb_e1   <- k1on_per_nm_h * a_e1 * fcrn_e1 - k1off_per_h * fcrna_e1
    fcb_e1_c <- k1on_per_nm_h * atc_e1 * fcrn_e1 - k1off_per_h * fcrnatc_e1
    bind_ef1 <- keon_per_nm_h * fcrna_e1 * t_e1 - keoff_per_h * fcrnatc_e1

    # Leaky tissue vascular space
    bind_v2  <- kon_per_nm_h * a_v2 * t_v2 - koff_per_h * atc_v2
    bind_m2  <- kon_per_nm_h * a_v2 * tm_v2 - koff_per_h * atm_v2
    bind_r2  <- kon_nsb * a_v2 * rm_total - koff_nsb * arm_v2

    # Leaky tissue endosome
    bind_e2  <- keon_per_nm_h * a_e2 * t_e2 - keoff_per_h * atc_e2
    fcb_e2   <- k1on_per_nm_h * a_e2 * fcrn_e2 - k1off_per_h * fcrna_e2
    fcb_e2_c <- k1on_per_nm_h * atc_e2 * fcrn_e2 - k1off_per_h * fcrnatc_e2
    bind_ef2 <- keon_per_nm_h * fcrna_e2 * t_e2 - keoff_per_h * fcrnatc_e2

    # =====================================================================
    # Plasma compartment (Patidar 2024 eqs. A1-A3)
    # =====================================================================
    d/dt(a_p) <- krec * ve / vp * fcrna_e -
                 spino_p * kup_p_per_h * a_p +
                 (ltot * a_l - (1 - sigma_tight) * l1 * a_p - (1 - sigma_leaky) * l2 * a_p) / vp -
                 bind_p

    # A2 as printed is -(kpt - Spino*kup_p)*T_p - Spino*kup_p*T_p, which sums
    # to -kpt*T_p: the pinocytosed fraction is carried into the nested
    # endosome by the kup_p term in d/dt(t_e).
    d/dt(t_p) <- ksyn - bind_p - kpt_per_h * t_p

    d/dt(atc_p) <- bind_p +
                   (ltot * atc_l - (1 - sigma_tight) * l1 * atc_p - (1 - sigma_leaky) * l2 * atc_p) / vp +
                   krec * ve / vp * fcrnatc_e -
                   spino_p * kup_p_per_h * atc_p

    # =====================================================================
    # Plasma nested endosome (Patidar 2024 eqs. A4-A9)
    #
    # Pinocytosis fluxes carry an explicit source-volume / endosome-volume
    # ratio, which the printed equations omit. Plasma loses kup_p * C_plasma
    # per unit PLASMA volume (eq. A1), so the endosome must gain
    # kup_p * C_plasma * vp/ve per unit ENDOSOMAL volume for antibody mass to
    # be conserved across the pinocytosis step; the endosome is ~600-fold
    # smaller than plasma in the mouse, so taking eq. A4 literally discards
    # >99% of the pinocytosed antibody. The volume-corrected form is the one
    # that reproduces the paper's observed data: it gives mouse plasma
    # concentrations of 370 / 326 / 298 / 278 nM at 24 / 48 / 96 / 168 h
    # against digitized observed values of 406 / 370 / 302 / 293 nM, and human
    # adalimumab 183 nM at 168 h against ~170 nM in Fig. 7, whereas the
    # literal form gives 21 nM and 0.5 nM respectively. See the vignette
    # "Assumptions and deviations" section.
    #
    # Eq. A4 prints kup where Table 1 and the Results distinguish kup_p
    # (plasma nested endosome, 0.05 1/h) from kup (tissue endosomes,
    # 0.0276 1/h); the plasma loss term in A1 is kup_p, so kup_p is used here.
    # Eqs. A5/A6 name the source species as T_v and ATC_e, where the plasma
    # sources T_p and ATC_p are required.
    # =====================================================================
    d/dt(a_e) <- spino_p * kup_p_per_h * a_p * vp / ve -
                 fcb_e - kdeg_per_h * a_e - bind_e

    d/dt(t_e) <- spino_p * kup_p_per_h * t_p * vp / ve -
                 bind_e - bind_ef - clcat / ve * t_e

    d/dt(atc_e) <- spino_p * kup_p_per_h * atc_p * vp / ve -
                   fcb_e_c - clcat / ve * atc_e + bind_e

    d/dt(fcrn_e) <- -fcb_e - fcb_e_c + krec * (fcrna_e + fcrnatc_e)

    d/dt(fcrna_e) <- fcb_e - bind_ef - krec * fcrna_e

    d/dt(fcrnatc_e) <- fcb_e_c + bind_ef - krec * fcrnatc_e

    # =====================================================================
    # Tight tissue -- vascular space (Patidar 2024 eqs. A10-A15)
    # A13 is printed with k_p,T and +k_off*ATC_v1; Table 1 carries a distinct
    # membrane-antigen degradation constant k_p,Tm, and mass balance with A15
    # requires the membrane complex ATm_v1 here.
    # =====================================================================
    d/dt(a_v1) <- (1 - sigma_tight) * l1 * a_p / vv1 -
                  flux_tp1 / vv1 +
                  fr * krec * ve1 / vv1 * fcrna_e1 -
                  spino1 * kup_per_h * a_v1 -
                  bind_v1 - bind_r1 - bind_m1

    d/dt(atc_v1) <- (1 - sigma_tight) * l1 * atc_p / vv1 +
                    fr * krec * ve1 / vv1 * fcrnatc_e1 -
                    spino1 * kup_per_h * atc_v1 +
                    bind_v1

    d/dt(t_v1) <- ksyn - bind_v1 - kpt_per_h * t_v1

    d/dt(tm_v1) <- ksynm - bind_m1 - kptm_per_h * tm_v1

    d/dt(arm_v1) <- bind_r1

    d/dt(atm_v1) <- bind_m1 - kint_per_h * atm_v1

    # =====================================================================
    # Tight tissue -- endosome (Patidar 2024 eqs. A16-A21)
    # A16 draws uptake from the vascular space only, but A22 removes antibody
    # from the interstitium by the same pinocytosis process, so both sources
    # feed the endosome here; otherwise the interstitial loss is an
    # unaccounted elimination pathway. A17 prints T_v1 inside an endosomal
    # equation where T_e1 is required, and A34 (the leaky-tissue mirror of
    # A20) omits the operator before ke_off*FcRnATC.
    # =====================================================================
    d/dt(a_e1) <- spino1 * kup_per_h * (a_v1 * vv1 + a_is1 * vis1_eff) / ve1 -
                  fcb_e1 - kdeg_per_h * a_e1 - bind_e1

    d/dt(t_e1) <- spino1 * kup_per_h * t_v1 * vv1 / ve1 -
                  bind_e1 - bind_ef1 - clcat / ve1 * t_e1

    d/dt(atc_e1) <- spino1 * kup_per_h * (atc_v1 * vv1 + atc_is1 * vis1_eff) / ve1 -
                    fcb_e1_c - clcat / ve1 * atc_e1 + bind_e1

    d/dt(fcrn_e1) <- -fcb_e1 - fcb_e1_c + krec * (fcrna_e1 + fcrnatc_e1)

    d/dt(fcrna_e1) <- fcb_e1 - bind_ef1 - krec * fcrna_e1

    d/dt(fcrnatc_e1) <- fcb_e1_c + bind_ef1 - krec * fcrnatc_e1

    # =====================================================================
    # Tight tissue -- interstitium (Patidar 2024 eqs. A22-A23)
    # The two-pore influx is driven by the vascular concentration; A22/A36
    # print CL_TP*A_is/(Kp*Vis), i.e. the interstitial state feeding itself.
    # =====================================================================
    d/dt(a_is1) <- -(1 - sigma_large) * l1 * a_is1 / vis1_eff -
                   spino1 * kup_per_h * a_is1 +
                   (1 - fr) * krec * ve1 / vis1_eff * fcrna_e1 +
                   flux_tp1 / vis1_eff

    d/dt(atc_is1) <- -(1 - sigma_large) * l1 * atc_is1 / vis1_eff -
                     spino1 * kup_per_h * atc_is1 +
                     (1 - fr) * krec * ve1 / vis1_eff * fcrnatc_e1

    # =====================================================================
    # Leaky tissue -- vascular space (Patidar 2024 eqs. A24-A29)
    # =====================================================================
    d/dt(a_v2) <- (1 - sigma_leaky) * l2 * a_p / vv2 -
                  (cl_renal * a_v2 + flux_tp2) / vv2 +
                  fr * krec * ve2 / vv2 * fcrna_e2 -
                  spino2 * kup_per_h * a_v2 -
                  bind_v2 - bind_r2 - bind_m2

    d/dt(atc_v2) <- (1 - sigma_leaky) * l2 * atc_p / vv2 +
                    fr * krec * ve2 / vv2 * fcrnatc_e2 -
                    spino2 * kup_per_h * atc_v2 +
                    bind_v2

    d/dt(t_v2) <- ksyn - bind_v2 - kpt_per_h * t_v2

    d/dt(tm_v2) <- ksynm - bind_m2 - kptm_per_h * tm_v2

    d/dt(arm_v2) <- bind_r2

    d/dt(atm_v2) <- bind_m2 - kint_per_h * atm_v2

    # =====================================================================
    # Leaky tissue -- endosome (Patidar 2024 eqs. A30-A35)
    # =====================================================================
    d/dt(a_e2) <- spino2 * kup_per_h * (a_v2 * vv2 + a_is2 * vis2_eff) / ve2 -
                  fcb_e2 - kdeg_per_h * a_e2 - bind_e2

    d/dt(t_e2) <- spino2 * kup_per_h * t_v2 * vv2 / ve2 -
                  bind_e2 - bind_ef2 - clcat / ve2 * t_e2

    d/dt(atc_e2) <- spino2 * kup_per_h * (atc_v2 * vv2 + atc_is2 * vis2_eff) / ve2 -
                    fcb_e2_c - clcat / ve2 * atc_e2 + bind_e2

    d/dt(fcrn_e2) <- -fcb_e2 - fcb_e2_c + krec * (fcrna_e2 + fcrnatc_e2)

    d/dt(fcrna_e2) <- fcb_e2 - bind_ef2 - krec * fcrna_e2

    d/dt(fcrnatc_e2) <- fcb_e2_c + bind_ef2 - krec * fcrnatc_e2

    # =====================================================================
    # Leaky tissue -- interstitium (Patidar 2024 eqs. A36-A37)
    # =====================================================================
    d/dt(a_is2) <- -(1 - sigma_large) * l2 * a_is2 / vis2_eff -
                   spino2 * kup_per_h * a_is2 +
                   (1 - fr) * krec * ve2 / vis2_eff * fcrna_e2 +
                   flux_tp2 / vis2_eff

    d/dt(atc_is2) <- -(1 - sigma_large) * l2 * atc_is2 / vis2_eff -
                     spino2 * kup_per_h * atc_is2 +
                     (1 - fr) * krec * ve2 / vis2_eff * fcrnatc_e2

    # =====================================================================
    # Lymph compartment (Patidar 2024 eqs. A38-A39)
    # =====================================================================
    d/dt(a_l) <- (1 - sigma_large) * l1 * a_is1 / vlymph +
                 (1 - sigma_large) * l2 * a_is2 / vlymph -
                 ltot * a_l / vlymph

    d/dt(atc_l) <- (1 - sigma_large) * l1 * atc_is1 / vlymph +
                   (1 - sigma_large) * l2 * atc_is2 / vlymph -
                   ltot * atc_l / vlymph

    # =====================================================================
    # Dosing and initial conditions
    # =====================================================================
    # The published system is written in concentrations, so an IV dose in
    # nmol enters the plasma state divided by the plasma volume.
    f(a_p) <- 1 / vp

    # FcRn is present at its total concentration in every endosome, and
    # antigen sits at its baseline expression in plasma and in the tissue
    # vascular spaces (Patidar 2024 Suppl. Table 3 assumptions 1, 3 and 4).
    fcrn_e(0)  <- fcrn_b
    fcrn_e1(0) <- fcrn_b
    fcrn_e2(0) <- fcrn_b
    t_p(0)     <- icc_t
    t_v1(0)    <- icc_t
    t_v2(0)    <- icc_t
    tm_v1(0)   <- icc_tm
    tm_v2(0)   <- icc_tm

    # =====================================================================
    # Outputs
    # =====================================================================
    # Total antibody in plasma (free + soluble complex)
    Cc <- a_p + atc_p

    # Total tight- and leaky-tissue antibody concentrations
    # (Patidar 2024 eqs. 1-2), summing every antibody-containing species in
    # each sub-space
    Ctight <- ((a_e1 + atc_e1 + fcrna_e1 + fcrnatc_e1) * ve1 +
               (a_v1 + atc_v1 + arm_v1 + atm_v1) * vv1 +
               (a_is1 + atc_is1) * vis1) / (ve1 + vv1 + vis1)

    Cleaky <- ((a_e2 + atc_e2 + fcrna_e2 + fcrnatc_e2) * ve2 +
               (a_v2 + atc_v2 + arm_v2 + atm_v2) * vv2 +
               (a_is2 + atc_is2) * vis2) / (ve2 + vv2 + vis2)
  })
}
