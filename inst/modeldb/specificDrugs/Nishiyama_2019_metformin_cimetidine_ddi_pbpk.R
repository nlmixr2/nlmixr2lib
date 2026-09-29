Nishiyama_2019_metformin_cimetidine_ddi_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, transporter-mediated drug-drug interaction). The",
    "metformin and cimetidine PBPK models of Nishiyama 2019 solved",
    "simultaneously, with cimetidine acting as a competitive inhibitor of",
    "the metformin carriers OCT1 (hepatic basolateral), OCT2 (renal",
    "basolateral) and MATE1/MATE2-K (renal luminal). This is the model",
    "behind the paper's headline result: because the electrogenic OCT",
    "steps are driven by a CONSTANT rather than a concentration-dependent",
    "membrane potential, the observed interaction is reproduced using",
    "inhibition constants measured in vitro, whereas the earlier",
    "published models needed those constants lowered 8- to 500-fold. The",
    "simulation also identifies MATE inhibition, not OCT2 inhibition, as",
    "the mechanism: the predicted AUC, Cmax and renal-clearance ratios",
    "are insensitive to the OCT2 Ki across its whole reported range and",
    "strongly sensitive to the MATE Ki. Dosing the cimetidine states with",
    "zero amount recovers the control arm exactly, so the same model",
    "object gives both arms of the interaction. The metformin parameters",
    "are the 250 mg fit, which is the dose of the interaction study;",
    "beta_kidney is shipped at 0.1 with the matching RMATE/dif of 183,",
    "and the three other beta_kidney / RMATE/dif pairs the authors",
    "carried are tabulated in the vignette. The Ki values shipped are the",
    "in vitro geometric means; the vignette also runs the fitted in vivo",
    "MATE Ki of Table 3. No between-subject variability or residual-error",
    "model is reported, so propSd is a placeholder and there are no etas."
  )
  reference <- paste(
    "Nishiyama K, Toshimoto K, Lee W, Ishiguro N, Bister B, Sugiyama Y.",
    "Physiologically-Based Pharmacokinetic Modeling Analysis for",
    "Quantitative Prediction of Renal Transporter-Mediated Interactions",
    "Between Metformin and Cimetidine.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8(6):396-406.",
    "doi:10.1002/psp4.12398.",
    "The two ODE systems are transcribed from Supplementary Material S2",
    "(PSP4-8-396-s008.pdf) and the Supplemental Text (PSP4-8-396-s007.pdf);",
    "the competitive-inhibition term is Methods eq. 15. Metformin",
    "parameters are Table S1 (PSP4-8-396-s003.pdf) and Table 1 (250 mg",
    "panel); cimetidine parameters are Table S4 (PSP4-8-396-s006.pdf);",
    "physiology is Tables S2 and S3 (PSP4-8-396-s004.pdf and",
    "PSP4-8-396-s005.pdf); the in vitro Ki values and their ranges are the",
    "footnote to Table 2 and the fitted in vivo MATE Ki values are Table 3.",
    "The clinical interaction data reproduced are Somogyi A et al.,",
    "Br J Clin Pharmacol. 1987;23:545-551.",
    sep = " "
  )
  vignette <- "Nishiyama_2019_metformin_cimetidine"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  covariateData <- list()

  # Compartment vocabulary. Every state carries a `_met` or `_cim` suffix
  # because the two drugs are solved in one system, so none of them is the
  # bare canonical name; the stems are those of the two companion models
  # Nishiyama_2019_metformin_pbpk.R and Nishiyama_2019_cimetidine_pbpk.R,
  # which document them.
  paper_specific_compartments <- c(
    "intestine_cim",
    "a_feces_cim",
    "blood_cim",
    "muscle_cim",
    "skin_cim",
    "adipose_cim",
    "is_liver_cim",
    "int_liver_cim",
    "a_metab_cim",
    "glom_b_cim",
    "glom_u_cim",
    "pt1_b_cim",
    "pt2_b_cim",
    "pt3_b_cim",
    "pt1_c_cim",
    "pt2_c_cim",
    "pt3_c_cim",
    "pt1_u_cim",
    "pt2_u_cim",
    "pt3_u_cim",
    "dt_b_cim",
    "dt_c_cim",
    "dt_u_cim",
    "cd_b_cim",
    "cd_c_cim",
    "cd_u_cim",
    "urine_cim",
    "transit1_met",
    "intestine1_met",
    "intestine2_met",
    "intestine3_met",
    "a_feces_met",
    "plasma_met",
    "rbc_met",
    "muscle_p_met",
    "muscle_e_met",
    "skin_p_met",
    "skin_e_met",
    "adipose_p_met",
    "adipose_e_met",
    "is_liver1_p_met",
    "is_liver2_p_met",
    "is_liver3_p_met",
    "is_liver4_p_met",
    "is_liver5_p_met",
    "is_liver1_e_met",
    "is_liver2_e_met",
    "is_liver3_e_met",
    "is_liver4_e_met",
    "is_liver5_e_met",
    "int_liver1_met",
    "int_liver2_met",
    "int_liver3_met",
    "int_liver4_met",
    "int_liver5_met",
    "a_metab_met",
    "glom_p_met",
    "glom_e_met",
    "glom_u_met",
    "pt1_p_met",
    "pt2_p_met",
    "pt3_p_met",
    "pt1_e_met",
    "pt2_e_met",
    "pt3_e_met",
    "pt1_c_met",
    "pt2_c_met",
    "pt3_c_met",
    "pt1_u_met",
    "pt2_u_met",
    "pt3_u_met",
    "dt_p_met",
    "dt_e_met",
    "dt_c_met",
    "dt_u_met",
    "cd_p_met",
    "cd_e_met",
    "cd_c_met",
    "cd_u_met",
    "urine_met"
  )

  compartmentData <- list(
    intestine_cim = list(analyte = "cimetidine", units = "ug", specimen = "administration site", verified = TRUE),
    a_feces_cim = list(analyte = "cimetidine", units = "ug", specimen = "faeces", verified = TRUE),
    blood_cim = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    muscle_cim = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    skin_cim = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    adipose_cim = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    is_liver_cim = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    int_liver_cim = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    a_metab_cim = list(analyte = "cimetidine", units = "ug", specimen = "not applicable", verified = TRUE),
    glom_b_cim = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    glom_u_cim = list(analyte = "cimetidine", units = "ug", specimen = "urine", verified = TRUE),
    pt1_b_cim = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    pt2_b_cim = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    pt3_b_cim = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    pt1_c_cim = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    pt2_c_cim = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    pt3_c_cim = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    pt1_u_cim = list(analyte = "cimetidine", units = "ug", specimen = "urine", verified = TRUE),
    pt2_u_cim = list(analyte = "cimetidine", units = "ug", specimen = "urine", verified = TRUE),
    pt3_u_cim = list(analyte = "cimetidine", units = "ug", specimen = "urine", verified = TRUE),
    dt_b_cim = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    dt_c_cim = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    dt_u_cim = list(analyte = "cimetidine", units = "ug", specimen = "urine", verified = TRUE),
    cd_b_cim = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    cd_c_cim = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    cd_u_cim = list(analyte = "cimetidine", units = "ug", specimen = "urine", verified = TRUE),
    urine_cim = list(analyte = "cimetidine", units = "ug", specimen = "urine", verified = TRUE),
    transit1_met = list(analyte = "metformin", units = "ug", specimen = "administration site", verified = TRUE),
    intestine1_met = list(analyte = "metformin", units = "ug", specimen = "administration site", verified = TRUE),
    intestine2_met = list(analyte = "metformin", units = "ug", specimen = "administration site", verified = TRUE),
    intestine3_met = list(analyte = "metformin", units = "ug", specimen = "administration site", verified = TRUE),
    a_feces_met = list(analyte = "metformin", units = "ug", specimen = "faeces", verified = TRUE),
    plasma_met = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    rbc_met = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    muscle_p_met = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    muscle_e_met = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    skin_p_met = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    skin_e_met = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    adipose_p_met = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    adipose_e_met = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    is_liver1_p_met = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    is_liver2_p_met = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    is_liver3_p_met = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    is_liver4_p_met = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    is_liver5_p_met = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    is_liver1_e_met = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    is_liver2_e_met = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    is_liver3_e_met = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    is_liver4_e_met = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    is_liver5_e_met = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    int_liver1_met = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver2_met = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver3_met = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver4_met = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver5_met = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    a_metab_met = list(analyte = "metformin", units = "ug", specimen = "not applicable", verified = TRUE),
    glom_p_met = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    glom_e_met = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    glom_u_met = list(analyte = "metformin", units = "ug", specimen = "urine", verified = TRUE),
    pt1_p_met = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    pt2_p_met = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    pt3_p_met = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    pt1_e_met = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    pt2_e_met = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    pt3_e_met = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    pt1_c_met = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    pt2_c_met = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    pt3_c_met = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    pt1_u_met = list(analyte = "metformin", units = "ug", specimen = "urine", verified = TRUE),
    pt2_u_met = list(analyte = "metformin", units = "ug", specimen = "urine", verified = TRUE),
    pt3_u_met = list(analyte = "metformin", units = "ug", specimen = "urine", verified = TRUE),
    dt_p_met = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    dt_e_met = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    dt_c_met = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    dt_u_met = list(analyte = "metformin", units = "ug", specimen = "urine", verified = TRUE),
    cd_p_met = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    cd_e_met = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    cd_c_met = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    cd_u_met = list(analyte = "metformin", units = "ug", specimen = "urine", verified = TRUE),
    urine_met = list(analyte = "metformin", units = "ug", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 6L,
    disease_state = "healthy adults",
    dose_range = "250 mg metformin with and without 400 mg cimetidine, single oral doses",
    notes = paste(
      "The interaction reproduced is the six-subject crossover of",
      "Somogyi et al. 1987 (Br J Clin Pharmacol 23:545-551), in which",
      "cimetidine raised the metformin plasma AUC 1.47-fold and Cmax",
      "1.72-fold and lowered renal clearance to 0.72 of control. Body",
      "physiology is the standard 70 kg adult of Davies & Morris 1993",
      "(Table S2), so the six subjects enter only as the observed",
      "summary statistics the simulation is compared against, not as a",
      "fitted population."
    )
  )

  ini({
    # =========== METFORMIN (victim) -- 250 mg fit ===========

    # ---------------------------------------------------------------
    # Absorption -- fitted (Table 1, 250 mg panel, the dose used in the
    # cimetidine interaction study). ka was estimated from the 1,500 mg
    # data and reused here; ktrans and r_mate_dif are the 250 mg fits.
    # ---------------------------------------------------------------
    lka_met <- log(0.21)
    label("Intestinal absorption rate constant (1/h)") # Table 1, 1,500 mg: ka 0.21 +/- 0.013 /h
    lktrans_met <- log(0.61)
    label("Gastric transit rate constant into the first intestinal segment (1/h)") # Table 1, 250 mg: ktrans 0.61 +/- 0.044 /h
    fafg_met <- fixed(0.84)
    label("Intestinal availability Fa x Fg (fraction)") # Table S1: FaFg 0.84 at 250 mg; back-calculated from bioavailability

    # ---------------------------------------------------------------
    # Distribution -- Table S1. Kp values are the Rodgers-Leahy-Rowland
    # in silico tissue-to-plasma ratios, tabulated by the authors.
    # ---------------------------------------------------------------
    kp_muscle_met <- fixed(2.09)
    label("Muscle-to-plasma concentration ratio (unitless)") # Table S1: Kp,muscle 2.09
    kp_skin_met <- fixed(1.46)
    label("Skin-to-plasma concentration ratio (unitless)") # Table S1: Kp,skin 1.46
    kp_adipose_met <- fixed(0.27)
    label("Adipose-to-plasma concentration ratio (unitless)") # Table S1: Kp,adipose 0.27
    kin_rbc_met <- fixed(0.006)
    label("Plasma-to-erythrocyte partitioning rate constant (1/h)") # Table S1: kin,RBC 0.006 /h
    kout_rbc_met <- fixed(0.02)
    label("Erythrocyte-to-plasma partitioning rate constant (1/h)") # Table S1: kout,RBC 0.02 /h

    # ---------------------------------------------------------------
    # Liver -- Table S1 and the Supplemental Text eqs. (1)-(7).
    # ---------------------------------------------------------------
    clint_all_met <- fixed(10.7)
    label("Overall hepatic intrinsic clearance (L/h)") # Table S1: CLint,all 10.7 L/h, held constant; back-calculated from the 250 mg intravenous data of Tucker 1981
    rdif_met <- fixed(0.186)
    label("Hepatic passive-to-active uptake clearance ratio (unitless)") # Table S1: Rdif 0.186
    beta_liver_met <- fixed(0.5)
    label("Hepatic elimination fraction beta (unitless)") # Table S1: beta,liver 0.5, held constant; fitting was insensitive to 0.2 / 0.5 / 0.8
    r_oct1_met <- fixed(1.32)
    label("OCT1 influx-to-efflux clearance ratio (unitless)") # Table S1: ROCT1,inf/eff 1.32
    km_oct1_um_met <- fixed(1470)
    label("OCT1 Michaelis constant for metformin (umol/L)") # Table S1 in vitro Km: OCT1 1,470 umol/L
    volt_h <- fixed(-40)
    label("Hepatocyte plasma-membrane potential (mV)") # Table S1: membrane potential -40 mV

    # ---------------------------------------------------------------
    # Kidney -- Table S1, Table S3 and Methods eqs. (5)-(14).
    # ---------------------------------------------------------------
    pd_met <- fixed(1.8e-5)
    label("Passive permeability from the PAMPA assay (m/h)") # Table S1: Pd 1.8e-5 m/h
    r_oct2_met <- fixed(1.32)
    label("OCT2 influx-to-efflux clearance ratio (unitless)") # Table S1: ROCT2,inf/eff 1.32
    km_oct2_um_met <- fixed(1178)
    label("OCT2 Michaelis constant for metformin (umol/L)") # Table S1 in vitro Km: OCT2 geometric mean 1,178 umol/L (range 810-1,465)
    km_mate_um_met <- fixed(740)
    label("MATE Michaelis constant for metformin (umol/L)") # Table S1 in vitro Km: MATEs geometric mean 740 umol/L (range 283-1,980)
    beta_kidney_met <- fixed(0.1)
    label("Renal secretion rate-determining fraction beta (unitless)") # Methods: beta,kidney not identifiable, held at 0.1, 0.3, 0.5 or 0.8; 0.1 shipped
    lr_mate_dif_met <- log(183)
    label("Ratio of MATE intrinsic clearance to luminal passive efflux clearance (unitless)") # Table 1, 250 mg, beta,kidney 0.1: RMATE/dif 183 +/- 30.1

    # ---------------------------------------------------------------
    # Body physiology -- Table S2, the standard 70 kg adult of
    # Davies & Morris 1993.
    # ---------------------------------------------------------------
    ht <- fixed(0.44)
    label("Haematocrit (fraction)") # Table S2: Ht 0.44
    vblood <- fixed(5.20)
    label("Blood volume (L)") # Table S2: Vblood 5.20 L
    vhc <- fixed(1.22)
    label("Hepatocyte volume (L)") # Table S2: VHC 1.22 L
    veh <- fixed(0.469)
    label("Hepatic extracellular volume (L)") # Table S2: VEH 0.469 L
    vmuscle <- fixed(35.0)
    label("Muscle volume (L)") # Table S2: Vmuscle 35.0 L
    vskin <- fixed(7.80)
    label("Skin volume (L)") # Table S2: Vskin 7.80 L
    vadipose <- fixed(10.0)
    label("Adipose volume (L)") # Table S2: Vadipose 10.0 L
    vkidney <- fixed(0.28)
    label("Kidney volume (L)") # Table S2: Vr 0.28 L
    qh <- fixed(87.0)
    label("Hepatic blood flow (L/h)") # Table S2: Qh 87.0 L/h
    qmuscle <- fixed(45.0)
    label("Muscle blood flow (L/h)") # Table S2: Qmuscle 45.0 L/h
    qskin <- fixed(18.0)
    label("Skin blood flow (L/h)") # Table S2: Qskin 18.0 L/h
    qadipose <- fixed(15.6)
    label("Adipose blood flow (L/h)") # Table S2: Qadipose 15.6 L/h

    # ---------------------------------------------------------------
    # Kidney physiology -- Table S3.
    # ---------------------------------------------------------------
    qr <- fixed(74.4)
    label("Renal blood flow (L/h)") # Table S3: Qr 74.4 L/h
    qgfr <- fixed(7.50)
    label("Glomerular filtration rate (L/h)") # Table S3: QGFR 7.50 L/h
    qurine <- fixed(0.06)
    label("Urine flow (L/h)") # Table S3: Qurine 0.06 L/h
    sa_r_pt <- fixed(0.81)
    label("Proximal-tubule basolateral surface area (m^2)") # Table S3: proximal tubule vessel 0.81 m^2
    sa_u_pt <- fixed(6.1)
    label("Proximal-tubule luminal surface area (m^2)") # Table S3: proximal tubule lumen 6.1 m^2
    sa_r_dt <- fixed(0.21)
    label("Distal-tubule basolateral surface area (m^2)") # Table S3: distal tubule vessel 0.21 m^2
    sa_u_dt <- fixed(0.21)
    label("Distal-tubule luminal surface area (m^2)") # Table S3: distal tubule lumen 0.21 m^2
    sa_r_cd <- fixed(0.045)
    label("Collecting-duct basolateral surface area (m^2)") # Table S3: collecting duct vessel 0.045 m^2
    sa_u_cd <- fixed(0.045)
    label("Collecting-duct luminal surface area (m^2)") # Table S3: collecting duct lumen 0.045 m^2
    volt_vpt <- fixed(-70)
    label("Proximal-tubule basolateral membrane potential (mV)") # Table S3: proximal tubule vessel -70 mV
    volt_upt <- fixed(-60)
    label("Proximal-tubule luminal membrane potential (mV)") # Table S3: proximal tubule lumen -60 mV
    volt_vdt <- fixed(-60)
    label("Distal-tubule basolateral membrane potential (mV)") # Table S3: distal tubule vessel -60 mV
    volt_udt <- fixed(-50)
    label("Distal-tubule luminal membrane potential (mV)") # Table S3: distal tubule lumen -50 mV
    volt_vcd <- fixed(-70)
    label("Collecting-duct basolateral membrane potential (mV)") # Table S3: collecting duct vessel -70 mV
    volt_ucd <- fixed(-30)
    label("Collecting-duct luminal membrane potential (mV)") # Table S3: collecting duct lumen -30 mV

    # ---------------------------------------------------------------
    # Residual error. The paper fits by weighted least squares in NAPP
    # and reports no residual-error model, so this is a placeholder held
    # constant to keep the model simulable; see the vignette Errata.
    # ---------------------------------------------------------------

    # =========== CIMETIDINE (perpetrator) ===========

    # ---------------------------------------------------------------
    # Physicochemical and absorption -- Table S4. None of the cimetidine
    # parameters was fitted in this paper: Table S4 cites Burt 2016 for
    # the compound layer, so every value is held constant.
    # ---------------------------------------------------------------
    pka_cim <- fixed(6.9)
    label("Acid dissociation constant (unitless)") # Table S4: pKa 6.9
    lam_cim <- fixed(0.1)
    label("Passive diffusion of the ionised relative to the unionised species (unitless)") # Supplemental Text: lambda set to 0.1 after Yoshikado 2017
    lka_cim <- fixed(log(0.70))
    label("Intestinal absorption rate constant (1/h)") # Table S4: ka 0.70 /h
    ltlag_cim <- fixed(log(0.15))
    label("Intestinal absorption lag time (h)") # Table S4: Tlag 0.15 h
    fafg_cim <- fixed(0.92)
    label("Intestinal availability Fa x Fg (fraction)") # Table S4: FaFg 0.92

    # ---------------------------------------------------------------
    # Distribution -- Table S4. The Kp values are referenced to BLOOD,
    # not plasma, because the cimetidine model carries a single blood
    # compartment.
    # ---------------------------------------------------------------
    fu_cim <- fixed(0.80)
    label("Unbound fraction in plasma (fraction)") # Table S4: fu,cim 0.80
    rb_cim <- fixed(0.97)
    label("Blood-to-plasma concentration ratio (unitless)") # Table S4: Rb 0.97
    kp_muscle_cim <- fixed(0.86)
    label("Muscle-to-blood concentration ratio (unitless)") # Table S4: Kp,muscle 0.86
    kp_skin_cim <- fixed(0.72)
    label("Skin-to-blood concentration ratio (unitless)") # Table S4: Kp,skin 0.72
    kp_adipose_cim <- fixed(0.24)
    label("Adipose-to-blood concentration ratio (unitless)") # Table S4: Kp,adipose 0.24

    # ---------------------------------------------------------------
    # Liver -- Table S4. Uptake is linear, so there is no hepatic Km.
    # ---------------------------------------------------------------
    ps_act_cim <- fixed(12.0)
    label("Hepatic active uptake clearance (L/h)") # Table S4: PSact 12.0 L/h
    rdif_cim <- fixed(1.16)
    label("Hepatic passive-to-active uptake clearance ratio (unitless)") # Table S4: Rdif 1.16
    cl_met_h_cim <- fixed(11.3)
    label("Hepatic metabolic intrinsic clearance (L/h)") # Table S4: CLmet 11.3 L/h

    # ---------------------------------------------------------------
    # Kidney -- Table S4. The Vmax values already carry the relative
    # activity factor of Burt 2016, so they are whole-organ maxima.
    # ---------------------------------------------------------------
    pd_cim <- fixed(7.9e-5)
    label("Passive permeability from the PAMPA assay (m/h)") # Table S4: Pd 7.9e-5 m/h
    r_oct2_cim <- fixed(1.32)
    label("OCT2 influx-to-efflux clearance ratio (unitless)") # Table S4: ROCT2,inf/eff 1.32
    km_oct2_um_cim <- fixed(72.6)
    label("OCT2 Michaelis constant for cimetidine (umol/L)") # Table S4: Km,OCT2 72.6 umol/L
    vmax_oct2_umol_cim <- fixed(7265)
    label("OCT2 maximum transport rate (umol/h)") # Table S4: Vmax,OCT2 7,265 umol/h
    km_oat3_um_cim <- fixed(161)
    label("OAT3 Michaelis constant for cimetidine (umol/L)") # Table S4: Km,OAT3 161 umol/L
    vmax_oat3_umol_cim <- fixed(4124)
    label("OAT3 maximum transport rate (umol/h)") # Table S4: Vmax,OAT3 4,124 umol/h
    km_mate_um_cim <- fixed(7.7)
    label("MATE Michaelis constant for cimetidine (umol/L)") # Table S4: Km,MATE 7.7 umol/L
    vmax_mate_umol_cim <- fixed(453)
    label("MATE maximum transport rate (umol/h)") # Table S4: Vmax,MATE 453 umol/h

    # ---------------------------------------------------------------
    # Body physiology -- Table S2 (shared with the metformin model).
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    # Kidney physiology -- Table S3 (shared with the metformin model).
    # ---------------------------------------------------------------
    ph_b_cim <- fixed(7.4)
    label("Vascular pH (unitless)") # Table S3: vessel pH 7.4
    ph_c_cim <- fixed(7.2)
    label("Tubular-cell pH (unitless)") # Table S3: cell pH 7.2
    ph_u_pt_cim <- fixed(7.0)
    label("Proximal-tubule luminal pH (unitless)") # Table S3: proximal tubule lumen pH 7.0
    ph_u_dt_cim <- fixed(6.7)
    label("Distal-tubule luminal pH (unitless)") # Table S3: distal tubule lumen pH 6.7
    ph_u_cd_cim <- fixed(6.4)
    label("Collecting-duct luminal pH (unitless)") # Table S3: collecting duct lumen pH 6.4
    # ---------------------------------------------------------------
    # Inhibition constants of cimetidine for the metformin carriers.
    # Shipped values are the in vitro geometric means of the footnote to
    # Table 2. Ranges are OCT2 72.6-509 umol/L and MATEs 1.22-13.5
    # umol/L; OCT1 has a single reported value.
    # ---------------------------------------------------------------
    ki_oct1_um <- fixed(104)
    label("Cimetidine inhibition constant for OCT1 (umol/L)") # Table 2 footnote: OCT1 Ki 104 umol/L
    ki_oct2_um <- fixed(159)
    label("Cimetidine inhibition constant for OCT2 (umol/L)") # Table 2 footnote: OCT2 Ki geometric mean 159 umol/L (range 72.6-509)
    ki_mate_um <- fixed(3.93)
    label("Cimetidine inhibition constant for MATEs (umol/L)") # Table 2 footnote: MATEs Ki geometric mean 3.93 umol/L (range 1.22-13.5)

    propSd <- fixed(0.1)
    label("Proportional residual error (fraction)") # not reported; placeholder
  })

  model({
    fara <- 96485 # Faraday constant, C/mol
    rgas <- 8.314 # gas constant, J/(mol K)
    tabs <- 310 # body temperature, K
    zval <- 1 # metformin valence at physiological pH
    nh <- zval * volt_h / 1000 * fara / (rgas * tabs)
    enh <- exp(nh)
    nvpt <- zval * volt_vpt / 1000 * fara / (rgas * tabs)
    envpt <- exp(nvpt)
    nupt <- zval * volt_upt / 1000 * fara / (rgas * tabs)
    enupt <- exp(nupt)
    nvdt <- zval * volt_vdt / 1000 * fara / (rgas * tabs)
    envdt <- exp(nvdt)
    nudt <- zval * volt_udt / 1000 * fara / (rgas * tabs)
    enudt <- exp(nudt)
    nvcd <- zval * volt_vcd / 1000 * fara / (rgas * tabs)
    envcd <- exp(nvcd)
    nucd <- zval * volt_ucd / 1000 * fara / (rgas * tabs)
    enucd <- exp(nucd)
    qu1 <- qgfr
    dqu <- (qgfr - qurine) / 5
    qu2 <- qu1 - dqu
    qu3 <- qu2 - dqu
    qu4 <- qu3 - dqu
    qu5 <- qu4 - dqu
    qu6 <- qurine
    qr1 <- qr - qgfr
    qr2 <- qr1 + dqu
    qr3 <- qr2 + dqu
    qr4 <- qr3 + dqu
    qr5 <- qr4 + dqu
    qr6 <- qr - qu6

    # ============ CIMETIDINE (perpetrator) ============

    mw_cim <- 252.34 # cimetidine_cim free_cim base_cim, g_cim/mol_cim; converts_cim Km_cim and_cim Vmax_cim from_cim umol_cim to_cim ug_cim

    ka_cim <- exp(lka_cim)
    tlag_cim <- exp(ltlag_cim)

    fr_ion_cim <- 1 / (1 + 10^(ph_b_cim - pka_cim))
    fr_union_cim <- fr_ion_cim * 10^(ph_b_cim - pka_cim)
    fc_ion_cim <- 1 / (1 + 10^(ph_c_cim - pka_cim))
    fc_union_cim <- fc_ion_cim * 10^(ph_c_cim - pka_cim)
    fu_ion_pt_cim <- 1 / (1 + 10^(ph_u_pt_cim - pka_cim))
    fu_union_pt_cim <- fu_ion_pt_cim * 10^(ph_u_pt_cim - pka_cim)
    fu_ion_dt_cim <- 1 / (1 + 10^(ph_u_dt_cim - pka_cim))
    fu_union_dt_cim <- fu_ion_dt_cim * 10^(ph_u_dt_cim - pka_cim)
    fu_ion_cd_cim <- 1 / (1 + 10^(ph_u_cd_cim - pka_cim))
    fu_union_cd_cim <- fu_ion_cd_cim * 10^(ph_u_cd_cim - pka_cim)


    gamma_h_cim <- (fr_union_cim + lam_cim * fr_ion_cim) / (fc_union_cim + enh * lam_cim * fc_ion_cim)
    ps_h_difinf_cim <- rdif_cim * ps_act_cim
    ps_h_difeff_cim <- ps_h_difinf_cim / gamma_h_cim

    pd_union_cim <- pd_cim / (fr_ion_cim * lam_cim + fr_union_cim)
    pd_ion_cim <- pd_union_cim * lam_cim
    ps_r_pt_difinf_cim <- pd_union_cim * sa_r_pt * 1000 * fr_union_cim +
      pd_ion_cim * sa_r_pt * 1000 * nvpt / (envpt - 1) * fr_ion_cim
    gamma_r_pt_cim <- (lam_cim * fr_ion_cim + fr_union_cim) / (envpt * lam_cim * fc_ion_cim + fc_union_cim)
    ps_r_pt_difeff_cim <- ps_r_pt_difinf_cim / gamma_r_pt_cim
    ps_u_pt_difinf_cim <- pd_union_cim * sa_u_pt * 1000 * fu_union_pt_cim +
      pd_ion_cim * sa_u_pt * 1000 * nupt / (enupt - 1) * fu_ion_pt_cim
    gamma_u_pt_cim <- (lam_cim * fu_ion_pt_cim + fu_union_pt_cim) / (enupt * lam_cim * fc_ion_cim + fc_union_cim)
    ps_u_pt_difeff_cim <- ps_u_pt_difinf_cim / gamma_u_pt_cim
    ps_r_dt_difinf_cim <- pd_union_cim * sa_r_dt * 1000 * fr_union_cim +
      pd_ion_cim * sa_r_dt * 1000 * nvdt / (envdt - 1) * fr_ion_cim
    gamma_r_dt_cim <- (lam_cim * fr_ion_cim + fr_union_cim) / (envdt * lam_cim * fc_ion_cim + fc_union_cim)
    ps_r_dt_difeff_cim <- ps_r_dt_difinf_cim / gamma_r_dt_cim
    ps_u_dt_difinf_cim <- pd_union_cim * sa_u_dt * 1000 * fu_union_dt_cim +
      pd_ion_cim * sa_u_dt * 1000 * nudt / (enudt - 1) * fu_ion_dt_cim
    gamma_u_dt_cim <- (lam_cim * fu_ion_dt_cim + fu_union_dt_cim) / (enudt * lam_cim * fc_ion_cim + fc_union_cim)
    ps_u_dt_difeff_cim <- ps_u_dt_difinf_cim / gamma_u_dt_cim
    ps_r_cd_difinf_cim <- pd_union_cim * sa_r_cd * 1000 * fr_union_cim +
      pd_ion_cim * sa_r_cd * 1000 * nvcd / (envcd - 1) * fr_ion_cim
    gamma_r_cd_cim <- (lam_cim * fr_ion_cim + fr_union_cim) / (envcd * lam_cim * fc_ion_cim + fc_union_cim)
    ps_r_cd_difeff_cim <- ps_r_cd_difinf_cim / gamma_r_cd_cim
    ps_u_cd_difinf_cim <- pd_union_cim * sa_u_cd * 1000 * fu_union_cd_cim +
      pd_ion_cim * sa_u_cd * 1000 * nucd / (enucd - 1) * fu_ion_cd_cim
    gamma_u_cd_cim <- (lam_cim * fu_ion_cd_cim + fu_union_cd_cim) / (enucd * lam_cim * fc_ion_cim + fc_union_cim)
    ps_u_cd_difeff_cim <- ps_u_cd_difinf_cim / gamma_u_cd_cim

    km_oct2_cim <- km_oct2_um_cim * mw_cim
    vmax_oct2_cim <- vmax_oct2_umol_cim * mw_cim
    km_oat3_cim <- km_oat3_um_cim * mw_cim
    vmax_oat3_cim <- vmax_oat3_umol_cim * mw_cim
    km_mate_cim <- km_mate_um_cim * mw_cim
    vmax_mate_cim <- vmax_mate_umol_cim * mw_cim


    vr_b_cim <- 0.30 * vkidney
    vr_cell_cim <- 0.24 * vkidney
    vr_u_cim <- 0.46 * vkidney
    v_glom_b_cim <- 5 / 129 * vr_b_cim
    v_pt_b_cim <- 17 / 129 * vr_b_cim
    v_dt_b_cim <- 11 / 129 * vr_b_cim
    v_cd_b_cim <- 62 / 129 * vr_b_cim
    v_glom_u_cim <- 5 / 129 * vr_u_cim
    v_pt_u_cim <- 17 / 129 * vr_u_cim
    v_dt_u_cim <- 11 / 129 * vr_u_cim
    v_cd_u_cim <- 62 / 129 * vr_u_cim
    v_pt_c_cim <- 17 / 124 * vr_cell_cim
    v_dt_c_cim <- 11 / 124 * vr_cell_cim
    v_cd_c_cim <- 62 / 124 * vr_cell_cim

    c_blood_cim <- blood_cim / vblood
    c_muscle_cim <- muscle_cim / vmuscle
    c_skin_cim <- skin_cim / vskin
    c_adipose_cim <- adipose_cim / vadipose
    c_eh_cim <- is_liver_cim / veh
    c_hc_cim <- int_liver_cim / vhc
    c_glom_b_cim <- glom_b_cim / v_glom_b_cim
    c_glom_u_cim <- glom_u_cim / v_glom_u_cim
    c_pt1_b_cim <- pt1_b_cim / (v_pt_b_cim / 3)
    c_pt2_b_cim <- pt2_b_cim / (v_pt_b_cim / 3)
    c_pt3_b_cim <- pt3_b_cim / (v_pt_b_cim / 3)
    c_pt1_c_cim <- pt1_c_cim / (v_pt_c_cim / 3)
    c_pt2_c_cim <- pt2_c_cim / (v_pt_c_cim / 3)
    c_pt3_c_cim <- pt3_c_cim / (v_pt_c_cim / 3)
    c_pt1_u_cim <- pt1_u_cim / (v_pt_u_cim / 3)
    c_pt2_u_cim <- pt2_u_cim / (v_pt_u_cim / 3)
    c_pt3_u_cim <- pt3_u_cim / (v_pt_u_cim / 3)
    c_dt_b_cim <- dt_b_cim / v_dt_b_cim
    c_dt_c_cim <- dt_c_cim / v_dt_c_cim
    c_dt_u_cim <- dt_u_cim / v_dt_u_cim
    c_cd_b_cim <- cd_b_cim / v_cd_b_cim
    c_cd_c_cim <- cd_c_cim / v_cd_c_cim
    c_cd_u_cim <- cd_u_cim / v_cd_u_cim

    cub1_cim <- fu_cim * fr_ion_cim * c_pt1_b_cim / rb_cim
    cub2_cim <- fu_cim * fr_ion_cim * c_pt2_b_cim / rb_cim
    cub3_cim <- fu_cim * fr_ion_cim * c_pt3_b_cim / rb_cim
    cuc1_cim <- fc_ion_cim * c_pt1_c_cim
    cuc2_cim <- fc_ion_cim * c_pt2_c_cim
    cuc3_cim <- fc_ion_cim * c_pt3_c_cim
    jo2_1_cim <- vmax_oct2_cim / 3 * (cub1_cim / (km_oct2_cim + cub1_cim) - cuc1_cim * envpt / r_oct2_cim / (km_oct2_cim + cuc1_cim))
    jo2_2_cim <- vmax_oct2_cim / 3 * (cub2_cim / (km_oct2_cim + cub2_cim) - cuc2_cim * envpt / r_oct2_cim / (km_oct2_cim + cuc2_cim))
    jo2_3_cim <- vmax_oct2_cim / 3 * (cub3_cim / (km_oct2_cim + cub3_cim) - cuc3_cim * envpt / r_oct2_cim / (km_oct2_cim + cuc3_cim))
    joat_1_cim <- vmax_oat3_cim / 3 * cub1_cim / (km_oat3_cim + cub1_cim)
    joat_2_cim <- vmax_oat3_cim / 3 * cub2_cim / (km_oat3_cim + cub2_cim)
    joat_3_cim <- vmax_oat3_cim / 3 * cub3_cim / (km_oat3_cim + cub3_cim)
    jmate_1_cim <- vmax_mate_cim / 3 * cuc1_cim / (km_mate_cim + cuc1_cim)
    jmate_2_cim <- vmax_mate_cim / 3 * cuc2_cim / (km_mate_cim + cuc2_cim)
    jmate_3_cim <- vmax_mate_cim / 3 * cuc3_cim / (km_mate_cim + cuc3_cim)

    d/dt(intestine_cim) <- -ka_cim / fafg_cim * intestine_cim
    alag(intestine_cim) <- tlag_cim
    d/dt(a_feces_cim) <- ka_cim * (1 - fafg_cim) / fafg_cim * intestine_cim

    d/dt(blood_cim) <- qh * c_eh_cim + qr6 * c_cd_b_cim - qh * c_blood_cim - qr * c_blood_cim -
      qmuscle * (c_blood_cim - c_muscle_cim / kp_muscle_cim) -
      qskin * (c_blood_cim - c_skin_cim / kp_skin_cim) -
      qadipose * (c_blood_cim - c_adipose_cim / kp_adipose_cim)
    d/dt(muscle_cim) <- qmuscle * (c_blood_cim - c_muscle_cim / kp_muscle_cim)
    d/dt(skin_cim) <- qskin * (c_blood_cim - c_skin_cim / kp_skin_cim)
    d/dt(adipose_cim) <- qadipose * (c_blood_cim - c_adipose_cim / kp_adipose_cim)

    d/dt(is_liver_cim) <- ka_cim * intestine_cim + qh * (c_blood_cim - c_eh_cim) + ps_h_difeff_cim * c_hc_cim -
      fu_cim * (ps_h_difinf_cim + ps_act_cim) * c_eh_cim / rb_cim
    d/dt(int_liver_cim) <- fu_cim * (ps_h_difinf_cim + ps_act_cim) * c_eh_cim / rb_cim -
      ps_h_difeff_cim * c_hc_cim - cl_met_h_cim * c_hc_cim
    d/dt(a_metab_cim) <- cl_met_h_cim * c_hc_cim

    d/dt(glom_b_cim) <- qr * c_blood_cim - qr1 * c_glom_b_cim - fu_cim * qu1 * c_glom_b_cim / rb_cim
    d/dt(glom_u_cim) <- qgfr * (fu_cim * c_glom_b_cim / rb_cim - c_glom_u_cim)

    d/dt(pt1_b_cim) <- qr1 * c_glom_b_cim + ps_r_pt_difeff_cim / 3 * c_pt1_c_cim - jo2_1_cim - joat_1_cim -
      ps_r_pt_difinf_cim / 3 * fu_cim * c_pt1_b_cim / rb_cim - qr2 * c_pt1_b_cim
    d/dt(pt2_b_cim) <- qr2 * c_pt1_b_cim + ps_r_pt_difeff_cim / 3 * c_pt2_c_cim - jo2_2_cim - joat_2_cim -
      ps_r_pt_difinf_cim / 3 * fu_cim * c_pt2_b_cim / rb_cim - qr3 * c_pt2_b_cim
    d/dt(pt3_b_cim) <- qr3 * c_pt2_b_cim + ps_r_pt_difeff_cim / 3 * c_pt3_c_cim - jo2_3_cim - joat_3_cim -
      ps_r_pt_difinf_cim / 3 * fu_cim * c_pt3_b_cim / rb_cim - qr4 * c_pt3_b_cim
    d/dt(pt1_c_cim) <- jo2_1_cim + joat_1_cim + ps_r_pt_difinf_cim / 3 * fu_cim * c_pt1_b_cim / rb_cim +
      ps_u_pt_difinf_cim / 3 * c_pt1_u_cim - ps_r_pt_difeff_cim / 3 * c_pt1_c_cim -
      jmate_1_cim - ps_u_pt_difeff_cim / 3 * c_pt1_c_cim
    d/dt(pt2_c_cim) <- jo2_2_cim + joat_2_cim + ps_r_pt_difinf_cim / 3 * fu_cim * c_pt2_b_cim / rb_cim +
      ps_u_pt_difinf_cim / 3 * c_pt2_u_cim - ps_r_pt_difeff_cim / 3 * c_pt2_c_cim -
      jmate_2_cim - ps_u_pt_difeff_cim / 3 * c_pt2_c_cim
    d/dt(pt3_c_cim) <- jo2_3_cim + joat_3_cim + ps_r_pt_difinf_cim / 3 * fu_cim * c_pt3_b_cim / rb_cim +
      ps_u_pt_difinf_cim / 3 * c_pt3_u_cim - ps_r_pt_difeff_cim / 3 * c_pt3_c_cim -
      jmate_3_cim - ps_u_pt_difeff_cim / 3 * c_pt3_c_cim
    d/dt(pt1_u_cim) <- qu1 * c_glom_u_cim + jmate_1_cim + ps_u_pt_difeff_cim / 3 * c_pt1_c_cim -
      ps_u_pt_difinf_cim / 3 * c_pt1_u_cim - qu2 * c_pt1_u_cim
    d/dt(pt2_u_cim) <- qu2 * c_pt1_u_cim + jmate_2_cim + ps_u_pt_difeff_cim / 3 * c_pt2_c_cim -
      ps_u_pt_difinf_cim / 3 * c_pt2_u_cim - qu3 * c_pt2_u_cim
    d/dt(pt3_u_cim) <- qu3 * c_pt2_u_cim + jmate_3_cim + ps_u_pt_difeff_cim / 3 * c_pt3_c_cim -
      ps_u_pt_difinf_cim / 3 * c_pt3_u_cim - qu4 * c_pt3_u_cim

    d/dt(dt_b_cim) <- qr4 * c_pt3_b_cim + ps_r_dt_difeff_cim * c_dt_c_cim -
      ps_r_dt_difinf_cim * fu_cim * c_dt_b_cim / rb_cim - qr5 * c_dt_b_cim
    d/dt(dt_c_cim) <- ps_r_dt_difinf_cim * fu_cim * c_dt_b_cim / rb_cim + ps_u_dt_difinf_cim * c_dt_u_cim -
      ps_r_dt_difeff_cim * c_dt_c_cim - ps_u_dt_difeff_cim * c_dt_c_cim
    d/dt(dt_u_cim) <- qu4 * c_pt3_u_cim + ps_u_dt_difeff_cim * c_dt_c_cim - ps_u_dt_difinf_cim * c_dt_u_cim - qu5 * c_dt_u_cim

    d/dt(cd_b_cim) <- qr5 * c_dt_b_cim + ps_r_cd_difeff_cim * c_cd_c_cim -
      ps_r_cd_difinf_cim * fu_cim * c_cd_b_cim / rb_cim - qr6 * c_cd_b_cim
    d/dt(cd_c_cim) <- ps_r_cd_difinf_cim * fu_cim * c_cd_b_cim / rb_cim + ps_u_cd_difinf_cim * c_cd_u_cim -
      ps_r_cd_difeff_cim * c_cd_c_cim - ps_u_cd_difeff_cim * c_cd_c_cim
    d/dt(cd_u_cim) <- qu5 * c_dt_u_cim + ps_u_cd_difeff_cim * c_cd_c_cim - ps_u_cd_difinf_cim * c_cd_u_cim - qu6 * c_cd_u_cim
    d/dt(urine_cim) <- qu6 * c_cd_u_cim



    # --- Competitive inhibition of metformin transport by cimetidine
    #     (Methods eq. 15). Eq. 15 divides the active clearance by
    #     (1 + I / Ki); because the carriers here are written as explicit
    #     Michaelis-Menten terms, the algebraically equivalent and
    #     mechanistically correct competitive form is used instead --
    #     Km is multiplied by (1 + I / Ki), which reproduces eq. 15
    #     exactly in the linear range where metformin sits (its Km values
    #     are three orders of magnitude above the concentrations reached).
    #     The inhibitor concentration at each site is the unbound
    #     cimetidine concentration on the blood side of the membrane for
    #     the basolateral carriers OCT1 and OCT2, and the tubular-cell
    #     concentration for the luminal carrier MATE, which MATE sees
    #     from its cis side. Both sides of each bidirectional carrier are
    #     inhibited, as competitive inhibition requires.
    ki_oct1 <- ki_oct1_um * mw_cim
    ki_oct2 <- ki_oct2_um * mw_cim
    ki_mate <- ki_mate_um * mw_cim
    inh_oct1 <- 1 + fu_cim * c_eh_cim / rb_cim / ki_oct1
    inh_oct2_1 <- 1 + fu_cim * c_pt1_b_cim / rb_cim / ki_oct2
    inh_oct2_2 <- 1 + fu_cim * c_pt2_b_cim / rb_cim / ki_oct2
    inh_oct2_3 <- 1 + fu_cim * c_pt3_b_cim / rb_cim / ki_oct2
    inh_mate_1 <- 1 + c_pt1_c_cim / ki_mate
    inh_mate_2 <- 1 + c_pt2_c_cim / ki_mate
    inh_mate_3 <- 1 + c_pt3_c_cim / ki_mate

    # ============ METFORMIN (victim) ============

    mw_met <- 129.16 # metformin_met free_met base_met, g_met/mol_met; converts_met in_met vitro_met Km_met from_met umol_met/L_met to_met ug_met/L_met

    ka_met <- exp(lka_met)
    ktrans_met <- exp(lktrans_met)
    r_mate_dif_met <- exp(lr_mate_dif_met)

    gamma_h_met <- 1 / enh

    ps_h_act_met <- clint_all_met / (beta_liver_met * (1 + rdif_met))
    ps_h_difinf_met <- ps_h_act_met * rdif_met
    ps_h_difeff_met <- ps_h_difinf_met / gamma_h_met
    cl_met_met <- clint_all_met / (1 - beta_liver_met) * rdif_met / ((1 + rdif_met) * gamma_h_met)
    km_oct1_met <- km_oct1_um_met * mw_met
    vmax_oct1_met <- ps_h_act_met * km_oct1_met

    ps_r_pt_difinf_met <- pd_met * sa_r_pt * 1000 * nvpt / (envpt - 1)
    ps_r_pt_difeff_met <- ps_r_pt_difinf_met * envpt
    ps_u_pt_difinf_met <- pd_met * sa_u_pt * 1000 * nupt / (enupt - 1)
    ps_u_pt_difeff_met <- ps_u_pt_difinf_met * enupt
    ps_r_dt_difinf_met <- pd_met * sa_r_dt * 1000 * nvdt / (envdt - 1)
    ps_r_dt_difeff_met <- ps_r_dt_difinf_met * envdt
    ps_u_dt_difinf_met <- pd_met * sa_u_dt * 1000 * nudt / (enudt - 1)
    ps_u_dt_difeff_met <- ps_u_dt_difinf_met * enudt
    ps_r_cd_difinf_met <- pd_met * sa_r_cd * 1000 * nvcd / (envcd - 1)
    ps_r_cd_difeff_met <- ps_r_cd_difinf_met * envcd
    ps_u_cd_difinf_met <- pd_met * sa_u_cd * 1000 * nucd / (enucd - 1)
    ps_u_cd_difeff_met <- ps_u_cd_difinf_met * enucd

    ps_mate_met <- ps_u_pt_difeff_met * r_mate_dif_met
    km_mate_met <- km_mate_um_met * mw_met
    vmax_mate_met <- ps_mate_met * km_mate_met
    ps_oct2_met <- ((r_mate_dif_met + 1) * ps_u_pt_difeff_met * (1 - beta_kidney_met) / beta_kidney_met -
      ps_r_pt_difeff_met) * r_oct2_met / envpt
    km_oct2_met <- km_oct2_um_met * mw_met
    vmax_oct2_met <- ps_oct2_met * km_oct2_met

    qh_p_met <- qh * (1 - ht)
    qh_e_met <- qh * ht
    qmuscle_p_met <- qmuscle * (1 - ht)
    qmuscle_e_met <- qmuscle * ht
    qskin_p_met <- qskin * (1 - ht)
    qskin_e_met <- qskin * ht
    qadipose_p_met <- qadipose * (1 - ht)
    qadipose_e_met <- qadipose * ht
    qr_p_met <- qr * (1 - ht)
    qr_e_met <- qr * ht
    qr1_p_met <- qr1 * (1 - ht)
    qr1_e_met <- qr1 * ht
    qr2_p_met <- qr2 * (1 - ht)
    qr2_e_met <- qr2 * ht
    qr3_p_met <- qr3 * (1 - ht)
    qr3_e_met <- qr3 * ht
    qr4_p_met <- qr4 * (1 - ht)
    qr4_e_met <- qr4 * ht
    qr5_p_met <- qr5 * (1 - ht)
    qr5_e_met <- qr5 * ht
    qr6_p_met <- qr6 * (1 - ht)
    qr6_e_met <- qr6 * ht

    v_plasma_met <- (1 - ht) * vblood
    v_rbc_met <- ht * vblood
    v_eh_p_met <- (1 - ht) * veh
    v_eh_e_met <- ht * veh
    v_muscle_p_met <- (1 - ht) * (0.026 + 0.12) * vmuscle + 0.854 * vmuscle
    v_muscle_e_met <- ht * (0.026 + 0.12) * vmuscle
    v_skin_p_met <- (1 - ht) * (0.019 + 0.302) * vskin + 0.679 * vskin
    v_skin_e_met <- ht * (0.019 + 0.302) * vskin
    v_adipose_p_met <- (1 - ht) * (0.01 + 0.135) * vadipose + 0.855 * vadipose
    v_adipose_e_met <- ht * (0.01 + 0.135) * vadipose
    vr_p_met <- 0.30 * (1 - ht) * vkidney
    vr_e_met <- 0.30 * ht * vkidney
    vr_cell_met <- 0.24 * vkidney
    vr_u_met <- 0.46 * vkidney
    v_glom_p_met <- 5 / 129 * vr_p_met
    v_pt_p_met <- 17 / 129 * vr_p_met
    v_dt_p_met <- 11 / 129 * vr_p_met
    v_cd_p_met <- 62 / 129 * vr_p_met
    v_glom_e_met <- 5 / 129 * vr_e_met
    v_pt_e_met <- 17 / 129 * vr_e_met
    v_dt_e_met <- 11 / 129 * vr_e_met
    v_cd_e_met <- 62 / 129 * vr_e_met
    v_glom_u_met <- 5 / 129 * vr_u_met
    v_pt_u_met <- 17 / 129 * vr_u_met
    v_dt_u_met <- 11 / 129 * vr_u_met
    v_cd_u_met <- 62 / 129 * vr_u_met
    v_pt_c_met <- 17 / 124 * vr_cell_met
    v_dt_c_met <- 11 / 124 * vr_cell_met
    v_cd_c_met <- 62 / 124 * vr_cell_met

    c_plasma_met <- plasma_met / v_plasma_met
    c_rbc_met <- rbc_met / v_rbc_met
    c_muscle_p_met <- muscle_p_met / v_muscle_p_met
    c_muscle_e_met <- muscle_e_met / v_muscle_e_met
    c_skin_p_met <- skin_p_met / v_skin_p_met
    c_skin_e_met <- skin_e_met / v_skin_e_met
    c_adipose_p_met <- adipose_p_met / v_adipose_p_met
    c_adipose_e_met <- adipose_e_met / v_adipose_e_met
    c_eh1_p_met <- is_liver1_p_met / (v_eh_p_met / 5)
    c_eh2_p_met <- is_liver2_p_met / (v_eh_p_met / 5)
    c_eh3_p_met <- is_liver3_p_met / (v_eh_p_met / 5)
    c_eh4_p_met <- is_liver4_p_met / (v_eh_p_met / 5)
    c_eh5_p_met <- is_liver5_p_met / (v_eh_p_met / 5)
    c_eh1_e_met <- is_liver1_e_met / (v_eh_e_met / 5)
    c_eh2_e_met <- is_liver2_e_met / (v_eh_e_met / 5)
    c_eh3_e_met <- is_liver3_e_met / (v_eh_e_met / 5)
    c_eh4_e_met <- is_liver4_e_met / (v_eh_e_met / 5)
    c_eh5_e_met <- is_liver5_e_met / (v_eh_e_met / 5)
    c_hc1_met <- int_liver1_met / (vhc / 5)
    c_hc2_met <- int_liver2_met / (vhc / 5)
    c_hc3_met <- int_liver3_met / (vhc / 5)
    c_hc4_met <- int_liver4_met / (vhc / 5)
    c_hc5_met <- int_liver5_met / (vhc / 5)
    c_glom_p_met <- glom_p_met / v_glom_p_met
    c_glom_e_met <- glom_e_met / v_glom_e_met
    c_glom_u_met <- glom_u_met / v_glom_u_met
    c_pt1_p_met <- pt1_p_met / (v_pt_p_met / 3)
    c_pt2_p_met <- pt2_p_met / (v_pt_p_met / 3)
    c_pt3_p_met <- pt3_p_met / (v_pt_p_met / 3)
    c_pt1_e_met <- pt1_e_met / (v_pt_e_met / 3)
    c_pt2_e_met <- pt2_e_met / (v_pt_e_met / 3)
    c_pt3_e_met <- pt3_e_met / (v_pt_e_met / 3)
    c_pt1_c_met <- pt1_c_met / (v_pt_c_met / 3)
    c_pt2_c_met <- pt2_c_met / (v_pt_c_met / 3)
    c_pt3_c_met <- pt3_c_met / (v_pt_c_met / 3)
    c_pt1_u_met <- pt1_u_met / (v_pt_u_met / 3)
    c_pt2_u_met <- pt2_u_met / (v_pt_u_met / 3)
    c_pt3_u_met <- pt3_u_met / (v_pt_u_met / 3)
    c_dt_p_met <- dt_p_met / v_dt_p_met
    c_dt_e_met <- dt_e_met / v_dt_e_met
    c_dt_c_met <- dt_c_met / v_dt_c_met
    c_dt_u_met <- dt_u_met / v_dt_u_met
    c_cd_p_met <- cd_p_met / v_cd_p_met
    c_cd_e_met <- cd_e_met / v_cd_e_met
    c_cd_c_met <- cd_c_met / v_cd_c_met
    c_cd_u_met <- cd_u_met / v_cd_u_met

    jo1_1_met <- vmax_oct1_met / 5 * (c_eh1_p_met / (km_oct1_met * inh_oct1 + c_eh1_p_met) - c_hc1_met * enh / r_oct1_met / (km_oct1_met * inh_oct1 + c_hc1_met))
    jo1_2_met <- vmax_oct1_met / 5 * (c_eh2_p_met / (km_oct1_met * inh_oct1 + c_eh2_p_met) - c_hc2_met * enh / r_oct1_met / (km_oct1_met * inh_oct1 + c_hc2_met))
    jo1_3_met <- vmax_oct1_met / 5 * (c_eh3_p_met / (km_oct1_met * inh_oct1 + c_eh3_p_met) - c_hc3_met * enh / r_oct1_met / (km_oct1_met * inh_oct1 + c_hc3_met))
    jo1_4_met <- vmax_oct1_met / 5 * (c_eh4_p_met / (km_oct1_met * inh_oct1 + c_eh4_p_met) - c_hc4_met * enh / r_oct1_met / (km_oct1_met * inh_oct1 + c_hc4_met))
    jo1_5_met <- vmax_oct1_met / 5 * (c_eh5_p_met / (km_oct1_met * inh_oct1 + c_eh5_p_met) - c_hc5_met * enh / r_oct1_met / (km_oct1_met * inh_oct1 + c_hc5_met))
    jo2_1_met <- vmax_oct2_met / 3 * (c_pt1_p_met / (km_oct2_met * inh_oct2_1 + c_pt1_p_met) - c_pt1_c_met * envpt / r_oct2_met / (km_oct2_met * inh_oct2_1 + c_pt1_c_met))
    jo2_2_met <- vmax_oct2_met / 3 * (c_pt2_p_met / (km_oct2_met * inh_oct2_2 + c_pt2_p_met) - c_pt2_c_met * envpt / r_oct2_met / (km_oct2_met * inh_oct2_2 + c_pt2_c_met))
    jo2_3_met <- vmax_oct2_met / 3 * (c_pt3_p_met / (km_oct2_met * inh_oct2_3 + c_pt3_p_met) - c_pt3_c_met * envpt / r_oct2_met / (km_oct2_met * inh_oct2_3 + c_pt3_c_met))
    jmate_1_met <- vmax_mate_met / 3 / (km_mate_met * inh_mate_1 + c_pt1_c_met) * c_pt1_c_met
    jmate_2_met <- vmax_mate_met / 3 / (km_mate_met * inh_mate_2 + c_pt2_c_met) * c_pt2_c_met
    jmate_3_met <- vmax_mate_met / 3 / (km_mate_met * inh_mate_3 + c_pt3_c_met) * c_pt3_c_met

    ktr_gi_met <- ka_met * (1 - fafg_met)^(1 / 3) / (1 - (1 - fafg_met)^(1 / 3))
    d/dt(transit1_met) <- -ktrans_met * transit1_met
    d/dt(intestine1_met) <- ktrans_met * transit1_met - ka_met * intestine1_met - ktr_gi_met * intestine1_met
    d/dt(intestine2_met) <- ktr_gi_met * intestine1_met - ka_met * intestine2_met - ktr_gi_met * intestine2_met
    d/dt(intestine3_met) <- ktr_gi_met * intestine2_met - ka_met * intestine3_met - ktr_gi_met * intestine3_met
    d/dt(a_feces_met) <- ktr_gi_met * intestine3_met

    d/dt(plasma_met) <- qh_p_met * c_eh5_p_met + qr6_p_met * c_cd_p_met + kout_rbc_met * rbc_met -
      qh_p_met * c_plasma_met - qr_p_met * c_plasma_met -
      qmuscle_p_met * (c_plasma_met - c_muscle_p_met / kp_muscle_met) -
      qskin_p_met * (c_plasma_met - c_skin_p_met / kp_skin_met) -
      qadipose_p_met * (c_plasma_met - c_adipose_p_met / kp_adipose_met) -
      kin_rbc_met * plasma_met
    d/dt(rbc_met) <- qh_e_met * c_eh5_e_met + qr6_e_met * c_cd_e_met + kin_rbc_met * plasma_met -
      qh_e_met * c_rbc_met - qr_e_met * c_rbc_met -
      qmuscle_e_met * (c_rbc_met - c_muscle_e_met) -
      qskin_e_met * (c_rbc_met - c_skin_e_met) -
      qadipose_e_met * (c_rbc_met - c_adipose_e_met) -
      kout_rbc_met * rbc_met

    d/dt(muscle_p_met) <- qmuscle_p_met * c_plasma_met + kout_rbc_met * muscle_e_met -
      kin_rbc_met * v_muscle_p_met * c_muscle_p_met / kp_muscle_met - qmuscle_p_met * c_muscle_p_met / kp_muscle_met
    d/dt(muscle_e_met) <- qmuscle_e_met * c_rbc_met + kin_rbc_met * v_muscle_p_met * c_muscle_p_met / kp_muscle_met -
      kout_rbc_met * muscle_e_met - qmuscle_e_met * c_muscle_e_met
    d/dt(skin_p_met) <- qskin_p_met * c_plasma_met + kout_rbc_met * skin_e_met -
      kin_rbc_met * v_skin_p_met * c_skin_p_met / kp_skin_met - qskin_p_met * c_skin_p_met / kp_skin_met
    d/dt(skin_e_met) <- qskin_e_met * c_rbc_met + kin_rbc_met * v_skin_p_met * c_skin_p_met / kp_skin_met -
      kout_rbc_met * skin_e_met - qskin_e_met * c_skin_e_met
    d/dt(adipose_p_met) <- qadipose_p_met * c_plasma_met + kout_rbc_met * adipose_e_met -
      kin_rbc_met * v_adipose_p_met * c_adipose_p_met / kp_adipose_met - qadipose_p_met * c_adipose_p_met / kp_adipose_met
    d/dt(adipose_e_met) <- qadipose_e_met * c_rbc_met + kin_rbc_met * v_adipose_p_met * c_adipose_p_met / kp_adipose_met -
      kout_rbc_met * adipose_e_met - qadipose_e_met * c_adipose_e_met

    d/dt(is_liver1_p_met) <- ka_met * (intestine1_met + intestine2_met + intestine3_met) +
      qh_p_met * (c_plasma_met - c_eh1_p_met) + kout_rbc_met * is_liver1_e_met +
      ps_h_difeff_met / 5 * c_hc1_met - jo1_1_met - ps_h_difinf_met / 5 * c_eh1_p_met - kin_rbc_met * is_liver1_p_met
    d/dt(is_liver2_p_met) <- qh_p_met * (c_eh1_p_met - c_eh2_p_met) + kout_rbc_met * is_liver2_e_met +
      ps_h_difeff_met / 5 * c_hc2_met - jo1_2_met - ps_h_difinf_met / 5 * c_eh2_p_met - kin_rbc_met * is_liver2_p_met
    d/dt(is_liver3_p_met) <- qh_p_met * (c_eh2_p_met - c_eh3_p_met) + kout_rbc_met * is_liver3_e_met +
      ps_h_difeff_met / 5 * c_hc3_met - jo1_3_met - ps_h_difinf_met / 5 * c_eh3_p_met - kin_rbc_met * is_liver3_p_met
    d/dt(is_liver4_p_met) <- qh_p_met * (c_eh3_p_met - c_eh4_p_met) + kout_rbc_met * is_liver4_e_met +
      ps_h_difeff_met / 5 * c_hc4_met - jo1_4_met - ps_h_difinf_met / 5 * c_eh4_p_met - kin_rbc_met * is_liver4_p_met
    d/dt(is_liver5_p_met) <- qh_p_met * (c_eh4_p_met - c_eh5_p_met) + kout_rbc_met * is_liver5_e_met +
      ps_h_difeff_met / 5 * c_hc5_met - jo1_5_met - ps_h_difinf_met / 5 * c_eh5_p_met - kin_rbc_met * is_liver5_p_met
    d/dt(is_liver1_e_met) <- qh_e_met * (c_rbc_met - c_eh1_e_met) + kin_rbc_met * is_liver1_p_met - kout_rbc_met * is_liver1_e_met
    d/dt(is_liver2_e_met) <- qh_e_met * (c_eh1_e_met - c_eh2_e_met) + kin_rbc_met * is_liver2_p_met - kout_rbc_met * is_liver2_e_met
    d/dt(is_liver3_e_met) <- qh_e_met * (c_eh2_e_met - c_eh3_e_met) + kin_rbc_met * is_liver3_p_met - kout_rbc_met * is_liver3_e_met
    d/dt(is_liver4_e_met) <- qh_e_met * (c_eh3_e_met - c_eh4_e_met) + kin_rbc_met * is_liver4_p_met - kout_rbc_met * is_liver4_e_met
    d/dt(is_liver5_e_met) <- qh_e_met * (c_eh4_e_met - c_eh5_e_met) + kin_rbc_met * is_liver5_p_met - kout_rbc_met * is_liver5_e_met
    d/dt(int_liver1_met) <- jo1_1_met + ps_h_difinf_met / 5 * c_eh1_p_met - ps_h_difeff_met / 5 * c_hc1_met - cl_met_met / 5 * c_hc1_met
    d/dt(int_liver2_met) <- jo1_2_met + ps_h_difinf_met / 5 * c_eh2_p_met - ps_h_difeff_met / 5 * c_hc2_met - cl_met_met / 5 * c_hc2_met
    d/dt(int_liver3_met) <- jo1_3_met + ps_h_difinf_met / 5 * c_eh3_p_met - ps_h_difeff_met / 5 * c_hc3_met - cl_met_met / 5 * c_hc3_met
    d/dt(int_liver4_met) <- jo1_4_met + ps_h_difinf_met / 5 * c_eh4_p_met - ps_h_difeff_met / 5 * c_hc4_met - cl_met_met / 5 * c_hc4_met
    d/dt(int_liver5_met) <- jo1_5_met + ps_h_difinf_met / 5 * c_eh5_p_met - ps_h_difeff_met / 5 * c_hc5_met - cl_met_met / 5 * c_hc5_met
    d/dt(a_metab_met) <- cl_met_met / 5 * (c_hc1_met + c_hc2_met + c_hc3_met + c_hc4_met + c_hc5_met)

    d/dt(glom_p_met) <- qr_p_met * c_plasma_met + kout_rbc_met * glom_e_met - qr1_p_met * c_glom_p_met -
      qu1 * c_glom_p_met - kin_rbc_met * glom_p_met
    d/dt(glom_e_met) <- qr_e_met * c_rbc_met + kin_rbc_met * glom_p_met - qr1_e_met * c_glom_e_met - kout_rbc_met * glom_e_met
    d/dt(glom_u_met) <- qgfr * (c_glom_p_met - c_glom_u_met)

    d/dt(pt1_p_met) <- qr1_p_met * c_glom_p_met + ps_r_pt_difeff_met / 3 * c_pt1_c_met + kout_rbc_met * pt1_e_met -
      jo2_1_met - ps_r_pt_difinf_met / 3 * c_pt1_p_met - qr2_p_met * c_pt1_p_met - kin_rbc_met * pt1_p_met
    d/dt(pt2_p_met) <- qr2_p_met * c_pt1_p_met + ps_r_pt_difeff_met / 3 * c_pt2_c_met + kout_rbc_met * pt2_e_met -
      jo2_2_met - ps_r_pt_difinf_met / 3 * c_pt2_p_met - qr3_p_met * c_pt2_p_met - kin_rbc_met * pt2_p_met
    d/dt(pt3_p_met) <- qr3_p_met * c_pt2_p_met + ps_r_pt_difeff_met / 3 * c_pt3_c_met + kout_rbc_met * pt3_e_met -
      jo2_3_met - ps_r_pt_difinf_met / 3 * c_pt3_p_met - qr4_p_met * c_pt3_p_met - kin_rbc_met * pt3_p_met
    d/dt(pt1_e_met) <- qr1_e_met * c_glom_e_met + kin_rbc_met * pt1_p_met - qr2_e_met * c_pt1_e_met - kout_rbc_met * pt1_e_met
    d/dt(pt2_e_met) <- qr2_e_met * c_pt1_e_met + kin_rbc_met * pt2_p_met - qr3_e_met * c_pt2_e_met - kout_rbc_met * pt2_e_met
    d/dt(pt3_e_met) <- qr3_e_met * c_pt2_e_met + kin_rbc_met * pt3_p_met - qr4_e_met * c_pt3_e_met - kout_rbc_met * pt3_e_met
    d/dt(pt1_c_met) <- jo2_1_met + ps_r_pt_difinf_met / 3 * c_pt1_p_met + ps_u_pt_difinf_met / 3 * c_pt1_u_met -
      ps_r_pt_difeff_met / 3 * c_pt1_c_met - jmate_1_met - ps_u_pt_difeff_met / 3 * c_pt1_c_met
    d/dt(pt2_c_met) <- jo2_2_met + ps_r_pt_difinf_met / 3 * c_pt2_p_met + ps_u_pt_difinf_met / 3 * c_pt2_u_met -
      ps_r_pt_difeff_met / 3 * c_pt2_c_met - jmate_2_met - ps_u_pt_difeff_met / 3 * c_pt2_c_met
    d/dt(pt3_c_met) <- jo2_3_met + ps_r_pt_difinf_met / 3 * c_pt3_p_met + ps_u_pt_difinf_met / 3 * c_pt3_u_met -
      ps_r_pt_difeff_met / 3 * c_pt3_c_met - jmate_3_met - ps_u_pt_difeff_met / 3 * c_pt3_c_met
    d/dt(pt1_u_met) <- qu1 * c_glom_u_met + jmate_1_met + ps_u_pt_difeff_met / 3 * c_pt1_c_met -
      ps_u_pt_difinf_met / 3 * c_pt1_u_met - qu2 * c_pt1_u_met
    d/dt(pt2_u_met) <- qu2 * c_pt1_u_met + jmate_2_met + ps_u_pt_difeff_met / 3 * c_pt2_c_met -
      ps_u_pt_difinf_met / 3 * c_pt2_u_met - qu3 * c_pt2_u_met
    d/dt(pt3_u_met) <- qu3 * c_pt2_u_met + jmate_3_met + ps_u_pt_difeff_met / 3 * c_pt3_c_met -
      ps_u_pt_difinf_met / 3 * c_pt3_u_met - qu4 * c_pt3_u_met

    d/dt(dt_p_met) <- qr4_p_met * c_pt3_p_met + ps_r_dt_difeff_met * c_dt_c_met + kout_rbc_met * dt_e_met -
      ps_r_dt_difinf_met * c_dt_p_met - qr5_p_met * c_dt_p_met - kin_rbc_met * dt_p_met
    d/dt(dt_e_met) <- qr4_e_met * c_pt3_e_met + kin_rbc_met * dt_p_met - qr5_e_met * c_dt_e_met - kout_rbc_met * dt_e_met
    d/dt(dt_c_met) <- ps_r_dt_difinf_met * c_dt_p_met + ps_u_dt_difinf_met * c_dt_u_met -
      ps_r_dt_difeff_met * c_dt_c_met - ps_u_dt_difeff_met * c_dt_c_met
    d/dt(dt_u_met) <- qu4 * c_pt3_u_met + ps_u_dt_difeff_met * c_dt_c_met - ps_u_dt_difinf_met * c_dt_u_met - qu5 * c_dt_u_met

    d/dt(cd_p_met) <- qr5_p_met * c_dt_p_met + ps_r_cd_difeff_met * c_cd_c_met + kout_rbc_met * cd_e_met -
      ps_r_cd_difinf_met * c_cd_p_met - qr6_p_met * c_cd_p_met - kin_rbc_met * cd_p_met
    d/dt(cd_e_met) <- qr5_e_met * c_dt_e_met + kin_rbc_met * cd_p_met - qr6_e_met * c_cd_e_met - kout_rbc_met * cd_e_met
    d/dt(cd_c_met) <- ps_r_cd_difinf_met * c_cd_p_met + ps_u_cd_difinf_met * c_cd_u_met -
      ps_r_cd_difeff_met * c_cd_c_met - ps_u_cd_difeff_met * c_cd_c_met
    d/dt(cd_u_met) <- qu5 * c_dt_u_met + ps_u_cd_difeff_met * c_cd_c_met - ps_u_cd_difinf_met * c_cd_u_met - qu6 * c_cd_u_met
    d/dt(urine_met) <- qu6 * c_cd_u_met



    # --- Observations. Cc is metformin in plasma (the DDI endpoint of
    #     Figure 3); Cb is metformin in whole blood and Ccim is
    #     cimetidine in plasma, the perpetrator exposure driving the
    #     inhibition.
    Cc <- c_plasma_met
    Cb <- (1 - ht) * c_plasma_met + ht * c_rbc_met
    Ccim <- c_blood_cim / rb_cim
    Cc ~ prop(propSd)
  })
}
