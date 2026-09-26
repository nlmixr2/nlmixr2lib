Nishiyama_2019_metformin_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, transporter-mediated, five-unit tandem dispersion",
    "liver, segmented nephron). Metformin after a single oral dose in",
    "healthy adults. Every systemic and capillary space is split into a",
    "plasma and an erythrocyte sub-compartment, because metformin",
    "distributes into red cells slowly enough that the blood-to-plasma",
    "ratio rises over the first hours after a dose. Absorption is a",
    "gastric transit compartment feeding three intestinal segments, each",
    "absorbing into the portal blood; the segment transit rate is set so",
    "that the three segments together deliver exactly the intestinal",
    "availability FaFg. The liver is five hepatic-extracellular units in",
    "series, each exchanging with its own hepatocyte unit through",
    "bidirectional OCT1 and passive diffusion. The kidney is a glomerulus",
    "plus a proximal tubule resolved into three segments, a distal tubule",
    "and a collecting duct, each carrying separate vascular, cellular and",
    "urinary-lumen spaces; active secretion is OCT2 at the basolateral",
    "membrane and MATE1/MATE2-K at the luminal membrane, and reabsorption",
    "is passive. The paper's central structural claim is that the",
    "electrogenic OCT1 and OCT2 steps are driven by a CONSTANT membrane",
    "potential, in contrast to the earlier 'electrochemical' model whose",
    "potential moved with metformin concentration; the constant-potential",
    "form is what lets the companion DDI model reproduce the cimetidine",
    "interaction with in vitro inhibition constants. Absorption is",
    "DOSE-SPECIFIC: the values below are the 1,500 mg fit (Table 1, upper",
    "panel). For the 250 mg dose use ktrans = 0.61 /h and fafg = 0.84; the",
    "vignette shows both. beta_kidney is not identifiable from the",
    "metformin data alone, so the authors fixed it at each of 0.1, 0.3,",
    "0.5 and 0.8 and refitted r_mate_dif; 0.1 is shipped here and the",
    "other three are tabulated in the vignette. The paper reports no",
    "between-subject variability and no residual-error model (the fit is",
    "a fixed-effects weighted least-squares fit in NAPP), so propSd is a",
    "placeholder and there are no etas."
  )
  reference <- paste(
    "Nishiyama K, Toshimoto K, Lee W, Ishiguro N, Bister B, Sugiyama Y.",
    "Physiologically-Based Pharmacokinetic Modeling Analysis for",
    "Quantitative Prediction of Renal Transporter-Mediated Interactions",
    "Between Metformin and Cimetidine.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8(6):396-406.",
    "doi:10.1002/psp4.12398.",
    "The ODE system and the hybrid-to-elementary parameter conversions are",
    "transcribed from Supplementary Material S2 ('Model equations for",
    "metformin', file PSP4-8-396-s008.pdf) and the Supplemental Text",
    "(PSP4-8-396-s007.pdf). Drug parameters are Table S1",
    "(PSP4-8-396-s003.pdf), body physiology is Table S2",
    "(PSP4-8-396-s004.pdf) and kidney physiology is Table S3",
    "(PSP4-8-396-s005.pdf). Fitted ka, ktrans and RMATE/dif are Table 1 of",
    "the article. See the vignette Errata for the transcription",
    "corrections applied to the published equation list.",
    sep = " "
  )
  vignette <- "Nishiyama_2019_metformin_cimetidine"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  # No covariates: the model is a single typical-subject physiology.
  covariateData <- list()

  # Compartment vocabulary. `transit1`, `plasma` and `urine` are canonical
  # (inst/references/compartment-names.md). `is_liver<n>` (hepatic
  # extracellular / sinusoidal space) and `int_liver<n>` (hepatocyte space)
  # follow the three registered models of this same Sugiyama-laboratory
  # framework -- Toshimoto_2017_irinotecan_pbpk.R,
  # Tsuchitani_2024_telmisartan_pbpk.R and Aoki_2024_bosentan_pbpk.R -- here
  # carrying a `_p` / `_e` suffix because this paper splits every vascular
  # space into plasma and erythrocyte. The nephron stems are paper-specific:
  # `glom_*` is the glomerulus, `pt1..3_*` the three proximal-tubule
  # segments, `dt_*` the distal tubule and `cd_*` the collecting duct, each
  # with `_p` (vascular plasma), `_e` (vascular erythrocyte), `_c`
  # (tubular cell) and `_u` (urinary lumen) spaces. `a_metab` and `a_feces`
  # are sinks, not biological matrices.
  paper_specific_compartments <- c(
    "intestine1",
    "intestine2",
    "intestine3",
    "a_feces",
    "a_metab",
    "rbc",
    "muscle_p",
    "muscle_e",
    "skin_p",
    "skin_e",
    "adipose_p",
    "adipose_e",
    "is_liver1_p",
    "is_liver2_p",
    "is_liver3_p",
    "is_liver4_p",
    "is_liver5_p",
    "is_liver1_e",
    "is_liver2_e",
    "is_liver3_e",
    "is_liver4_e",
    "is_liver5_e",
    "int_liver1",
    "int_liver2",
    "int_liver3",
    "int_liver4",
    "int_liver5",
    "glom_p",
    "glom_e",
    "glom_u",
    "pt1_p",
    "pt2_p",
    "pt3_p",
    "pt1_e",
    "pt2_e",
    "pt3_e",
    "pt1_c",
    "pt2_c",
    "pt3_c",
    "pt1_u",
    "pt2_u",
    "pt3_u",
    "dt_p",
    "dt_e",
    "dt_c",
    "dt_u",
    "cd_p",
    "cd_e",
    "cd_c",
    "cd_u"
  )

  compartmentData <- list(
    transit1 = list(analyte = "metformin", units = "ug", specimen = "administration site", verified = TRUE),
    intestine1 = list(analyte = "metformin", units = "ug", specimen = "administration site", verified = TRUE),
    intestine2 = list(analyte = "metformin", units = "ug", specimen = "administration site", verified = TRUE),
    intestine3 = list(analyte = "metformin", units = "ug", specimen = "administration site", verified = TRUE),
    a_feces = list(analyte = "metformin", units = "ug", specimen = "faeces", verified = TRUE),
    plasma = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    rbc = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    muscle_p = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    muscle_e = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    skin_p = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    skin_e = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    adipose_p = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    adipose_e = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    is_liver1_p = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    is_liver2_p = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    is_liver3_p = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    is_liver4_p = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    is_liver5_p = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    is_liver1_e = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    is_liver2_e = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    is_liver3_e = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    is_liver4_e = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    is_liver5_e = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    int_liver1 = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver2 = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver3 = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver4 = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver5 = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    a_metab = list(analyte = "metformin", units = "ug", specimen = "not applicable", verified = TRUE),
    glom_p = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    glom_e = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    glom_u = list(analyte = "metformin", units = "ug", specimen = "urine", verified = TRUE),
    pt1_p = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    pt2_p = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    pt3_p = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    pt1_e = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    pt2_e = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    pt3_e = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    pt1_c = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    pt2_c = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    pt3_c = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    pt1_u = list(analyte = "metformin", units = "ug", specimen = "urine", verified = TRUE),
    pt2_u = list(analyte = "metformin", units = "ug", specimen = "urine", verified = TRUE),
    pt3_u = list(analyte = "metformin", units = "ug", specimen = "urine", verified = TRUE),
    dt_p = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    dt_e = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    dt_c = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    dt_u = list(analyte = "metformin", units = "ug", specimen = "urine", verified = TRUE),
    cd_p = list(analyte = "metformin", units = "ug", specimen = "plasma", verified = TRUE),
    cd_e = list(analyte = "metformin", units = "ug", specimen = "blood cell", verified = TRUE),
    cd_c = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    cd_u = list(analyte = "metformin", units = "ug", specimen = "urine", verified = TRUE),
    urine = list(analyte = "metformin", units = "ug", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "healthy adults",
    dose_range = "250 mg and 1,500 mg metformin, single oral dose",
    notes = paste(
      "The model was fitted to two published single-dose data sets rather",
      "than to an individual-level population: the 1,500 mg plasma, blood",
      "and urinary-excretion profiles of Tucker et al. 1981",
      "(Br J Clin Pharmacol 12:235-246) and the 250 mg control arm of the",
      "cimetidine interaction study of Somogyi et al. 1987",
      "(Br J Clin Pharmacol 23:545-551). Body physiology is the standard",
      "70 kg adult of Davies & Morris 1993 (Table S2), so subject counts,",
      "demographics and covariate distributions do not enter the model."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Absorption -- fitted (Table 1, 1,500 mg panel). ka was estimated
    # from the 1,500 mg data and then reused for the 250 mg data;
    # ktrans and r_mate_dif were re-estimated per dose.
    # ---------------------------------------------------------------
    lka <- log(0.21)
    label("Intestinal absorption rate constant (1/h)") # Table 1, 1,500 mg: ka 0.21 +/- 0.013 /h
    lktrans <- log(2.4)
    label("Gastric transit rate constant into the first intestinal segment (1/h)") # Table 1, 1,500 mg: ktrans 2.4 +/- 0.36 /h
    fafg <- fixed(0.57)
    label("Intestinal availability Fa x Fg (fraction)") # Table S1: FaFg 0.57 at 1,500 mg; 0.84 at 250 mg; back-calculated from bioavailability

    # ---------------------------------------------------------------
    # Distribution -- Table S1. Kp values are the Rodgers-Leahy-Rowland
    # in silico tissue-to-plasma ratios, tabulated by the authors.
    # ---------------------------------------------------------------
    kp_muscle <- fixed(2.09)
    label("Muscle-to-plasma concentration ratio (unitless)") # Table S1: Kp,muscle 2.09
    kp_skin <- fixed(1.46)
    label("Skin-to-plasma concentration ratio (unitless)") # Table S1: Kp,skin 1.46
    kp_adipose <- fixed(0.27)
    label("Adipose-to-plasma concentration ratio (unitless)") # Table S1: Kp,adipose 0.27
    kin_rbc <- fixed(0.006)
    label("Plasma-to-erythrocyte partitioning rate constant (1/h)") # Table S1: kin,RBC 0.006 /h
    kout_rbc <- fixed(0.02)
    label("Erythrocyte-to-plasma partitioning rate constant (1/h)") # Table S1: kout,RBC 0.02 /h

    # ---------------------------------------------------------------
    # Liver -- Table S1 and the Supplemental Text eqs. (1)-(7).
    # ---------------------------------------------------------------
    clint_all <- fixed(10.7)
    label("Overall hepatic intrinsic clearance (L/h)") # Table S1: CLint,all 10.7 L/h, held constant; back-calculated from the 250 mg intravenous data of Tucker 1981
    rdif <- fixed(0.186)
    label("Hepatic passive-to-active uptake clearance ratio (unitless)") # Table S1: Rdif 0.186
    beta_liver <- fixed(0.5)
    label("Hepatic elimination fraction beta (unitless)") # Table S1: beta,liver 0.5, held constant; fitting was insensitive to 0.2 / 0.5 / 0.8
    r_oct1 <- fixed(1.32)
    label("OCT1 influx-to-efflux clearance ratio (unitless)") # Table S1: ROCT1,inf/eff 1.32
    km_oct1_um <- fixed(1470)
    label("OCT1 Michaelis constant for metformin (umol/L)") # Table S1 in vitro Km: OCT1 1,470 umol/L
    volt_h <- fixed(-40)
    label("Hepatocyte plasma-membrane potential (mV)") # Table S1: membrane potential -40 mV

    # ---------------------------------------------------------------
    # Kidney -- Table S1, Table S3 and Methods eqs. (5)-(14).
    # ---------------------------------------------------------------
    pd <- fixed(1.8e-5)
    label("Passive permeability from the PAMPA assay (m/h)") # Table S1: Pd 1.8e-5 m/h
    r_oct2 <- fixed(1.32)
    label("OCT2 influx-to-efflux clearance ratio (unitless)") # Table S1: ROCT2,inf/eff 1.32
    km_oct2_um <- fixed(1178)
    label("OCT2 Michaelis constant for metformin (umol/L)") # Table S1 in vitro Km: OCT2 geometric mean 1,178 umol/L (range 810-1,465)
    km_mate_um <- fixed(740)
    label("MATE Michaelis constant for metformin (umol/L)") # Table S1 in vitro Km: MATEs geometric mean 740 umol/L (range 283-1,980)
    beta_kidney <- fixed(0.1)
    label("Renal secretion rate-determining fraction beta (unitless)") # Methods: beta,kidney not identifiable, held at 0.1, 0.3, 0.5 or 0.8; 0.1 shipped
    lr_mate_dif <- log(153)
    label("Ratio of MATE intrinsic clearance to luminal passive efflux clearance (unitless)") # Table 1, 1,500 mg, beta,kidney 0.1: RMATE/dif 153 +/- 18.1

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
    propSd <- fixed(0.1)
    label("Proportional residual error (fraction)") # not reported; placeholder
  })

  model({
    # ===============================================================
    # Physical constants and fixed geometry. z = +1 because metformin
    # (pKa 12.3) is a monovalent cation at every physiological pH, so
    # the extracellular and intracellular ionised fractions are both 1
    # and every free-fraction multiplier in the published equations
    # collapses to unity (fu,met = 1, Table S1).
    # ===============================================================
    fara <- 96485 # Faraday constant, C/mol
    rgas <- 8.314 # gas constant, J/(mol K)
    tabs <- 310 # body temperature, K
    zval <- 1 # metformin valence at physiological pH
    mw <- 129.16 # metformin free base, g/mol; converts in vitro Km from umol/L to ug/L

    ka <- exp(lka)
    ktrans <- exp(lktrans)
    r_mate_dif <- exp(lr_mate_dif)

    # --- Nernst factors, one per membrane (Suppl. S2, 'Other equations') ---
    nh <- zval * volt_h / 1000 * fara / (rgas * tabs)
    enh <- exp(nh)
    gamma_h <- 1 / enh
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

    # --- Liver: hybrid CLint,all / Rdif / beta / gamma to elementary PS ---
    # Suppl. Text eqs. (1)-(5) and Suppl. S2 'Other equations'.
    ps_h_act <- clint_all / (beta_liver * (1 + rdif))
    ps_h_difinf <- ps_h_act * rdif
    ps_h_difeff <- ps_h_difinf / gamma_h
    cl_met <- clint_all / (1 - beta_liver) * rdif / ((1 + rdif) * gamma_h)
    km_oct1 <- km_oct1_um * mw
    vmax_oct1 <- ps_h_act * km_oct1

    # --- Kidney passive diffusion: Goldman flux with a constant potential ---
    ps_r_pt_difinf <- pd * sa_r_pt * 1000 * nvpt / (envpt - 1)
    ps_r_pt_difeff <- ps_r_pt_difinf * envpt
    ps_u_pt_difinf <- pd * sa_u_pt * 1000 * nupt / (enupt - 1)
    ps_u_pt_difeff <- ps_u_pt_difinf * enupt
    ps_r_dt_difinf <- pd * sa_r_dt * 1000 * nvdt / (envdt - 1)
    ps_r_dt_difeff <- ps_r_dt_difinf * envdt
    ps_u_dt_difinf <- pd * sa_u_dt * 1000 * nudt / (enudt - 1)
    ps_u_dt_difeff <- ps_u_dt_difinf * enudt
    ps_r_cd_difinf <- pd * sa_r_cd * 1000 * nvcd / (envcd - 1)
    ps_r_cd_difeff <- ps_r_cd_difinf * envcd
    ps_u_cd_difinf <- pd * sa_u_cd * 1000 * nucd / (enucd - 1)
    ps_u_cd_difeff <- ps_u_cd_difinf * enucd

    # --- Renal transporters, back-solved from RMATE/dif and beta,kidney ---
    ps_mate <- ps_u_pt_difeff * r_mate_dif
    km_mate <- km_mate_um * mw
    vmax_mate <- ps_mate * km_mate
    ps_oct2 <- ((r_mate_dif + 1) * ps_u_pt_difeff * (1 - beta_kidney) / beta_kidney -
      ps_r_pt_difeff) * r_oct2 / envpt
    km_oct2 <- km_oct2_um * mw
    vmax_oct2 <- ps_oct2 * km_oct2

    # --- Flows. Tubular water is reabsorbed in five equal steps between
    #     the glomerular filtrate and the final urine, and the reabsorbed
    #     water returns to the peritubular vasculature.
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
    qh_p <- qh * (1 - ht)
    qh_e <- qh * ht
    qmuscle_p <- qmuscle * (1 - ht)
    qmuscle_e <- qmuscle * ht
    qskin_p <- qskin * (1 - ht)
    qskin_e <- qskin * ht
    qadipose_p <- qadipose * (1 - ht)
    qadipose_e <- qadipose * ht
    qr_p <- qr * (1 - ht)
    qr_e <- qr * ht
    qr1_p <- qr1 * (1 - ht)
    qr1_e <- qr1 * ht
    qr2_p <- qr2 * (1 - ht)
    qr2_e <- qr2 * ht
    qr3_p <- qr3 * (1 - ht)
    qr3_e <- qr3 * ht
    qr4_p <- qr4 * (1 - ht)
    qr4_e <- qr4 * ht
    qr5_p <- qr5 * (1 - ht)
    qr5_e <- qr5 * ht
    qr6_p <- qr6 * (1 - ht)
    qr6_e <- qr6 * ht

    # --- Volumes. The tissue vascular fractions (0.026 + 0.12 for muscle,
    #     0.019 + 0.302 for skin, 0.01 + 0.135 for adipose) and the
    #     0.30 / 0.24 / 0.46 vascular / cellular / luminal split of the
    #     kidney are Suppl. S2 'Other equations'; the 5 / 17 / 11 / 62
    #     weights distribute the renal vascular and luminal space over
    #     the glomerulus, proximal tubule, distal tubule and collecting
    #     duct (denominator 129 for the vascular and luminal spaces, 124
    #     for the cellular space, which has no glomerular member).
    v_plasma <- (1 - ht) * vblood
    v_rbc <- ht * vblood
    v_eh_p <- (1 - ht) * veh
    v_eh_e <- ht * veh
    v_muscle_p <- (1 - ht) * (0.026 + 0.12) * vmuscle + 0.854 * vmuscle
    v_muscle_e <- ht * (0.026 + 0.12) * vmuscle
    v_skin_p <- (1 - ht) * (0.019 + 0.302) * vskin + 0.679 * vskin
    v_skin_e <- ht * (0.019 + 0.302) * vskin
    v_adipose_p <- (1 - ht) * (0.01 + 0.135) * vadipose + 0.855 * vadipose
    v_adipose_e <- ht * (0.01 + 0.135) * vadipose
    vr_p <- 0.30 * (1 - ht) * vkidney
    vr_e <- 0.30 * ht * vkidney
    vr_cell <- 0.24 * vkidney
    vr_u <- 0.46 * vkidney
    v_glom_p <- 5 / 129 * vr_p
    v_pt_p <- 17 / 129 * vr_p
    v_dt_p <- 11 / 129 * vr_p
    v_cd_p <- 62 / 129 * vr_p
    v_glom_e <- 5 / 129 * vr_e
    v_pt_e <- 17 / 129 * vr_e
    v_dt_e <- 11 / 129 * vr_e
    v_cd_e <- 62 / 129 * vr_e
    v_glom_u <- 5 / 129 * vr_u
    v_pt_u <- 17 / 129 * vr_u
    v_dt_u <- 11 / 129 * vr_u
    v_cd_u <- 62 / 129 * vr_u
    v_pt_c <- 17 / 124 * vr_cell
    v_dt_c <- 11 / 124 * vr_cell
    v_cd_c <- 62 / 124 * vr_cell

    # --- Concentrations (ug/L) from the amount states (ug) ---
    c_plasma <- plasma / v_plasma
    c_rbc <- rbc / v_rbc
    c_muscle_p <- muscle_p / v_muscle_p
    c_muscle_e <- muscle_e / v_muscle_e
    c_skin_p <- skin_p / v_skin_p
    c_skin_e <- skin_e / v_skin_e
    c_adipose_p <- adipose_p / v_adipose_p
    c_adipose_e <- adipose_e / v_adipose_e
    c_eh1_p <- is_liver1_p / (v_eh_p / 5)
    c_eh2_p <- is_liver2_p / (v_eh_p / 5)
    c_eh3_p <- is_liver3_p / (v_eh_p / 5)
    c_eh4_p <- is_liver4_p / (v_eh_p / 5)
    c_eh5_p <- is_liver5_p / (v_eh_p / 5)
    c_eh1_e <- is_liver1_e / (v_eh_e / 5)
    c_eh2_e <- is_liver2_e / (v_eh_e / 5)
    c_eh3_e <- is_liver3_e / (v_eh_e / 5)
    c_eh4_e <- is_liver4_e / (v_eh_e / 5)
    c_eh5_e <- is_liver5_e / (v_eh_e / 5)
    c_hc1 <- int_liver1 / (vhc / 5)
    c_hc2 <- int_liver2 / (vhc / 5)
    c_hc3 <- int_liver3 / (vhc / 5)
    c_hc4 <- int_liver4 / (vhc / 5)
    c_hc5 <- int_liver5 / (vhc / 5)
    c_glom_p <- glom_p / v_glom_p
    c_glom_e <- glom_e / v_glom_e
    c_glom_u <- glom_u / v_glom_u
    c_pt1_p <- pt1_p / (v_pt_p / 3)
    c_pt2_p <- pt2_p / (v_pt_p / 3)
    c_pt3_p <- pt3_p / (v_pt_p / 3)
    c_pt1_e <- pt1_e / (v_pt_e / 3)
    c_pt2_e <- pt2_e / (v_pt_e / 3)
    c_pt3_e <- pt3_e / (v_pt_e / 3)
    c_pt1_c <- pt1_c / (v_pt_c / 3)
    c_pt2_c <- pt2_c / (v_pt_c / 3)
    c_pt3_c <- pt3_c / (v_pt_c / 3)
    c_pt1_u <- pt1_u / (v_pt_u / 3)
    c_pt2_u <- pt2_u / (v_pt_u / 3)
    c_pt3_u <- pt3_u / (v_pt_u / 3)
    c_dt_p <- dt_p / v_dt_p
    c_dt_e <- dt_e / v_dt_e
    c_dt_c <- dt_c / v_dt_c
    c_dt_u <- dt_u / v_dt_u
    c_cd_p <- cd_p / v_cd_p
    c_cd_e <- cd_e / v_cd_e
    c_cd_c <- cd_c / v_cd_c
    c_cd_u <- cd_u / v_cd_u

    # --- Bidirectional, membrane-potential-driven transporter fluxes
    #     (Methods eqs. 3-6). Each is written once and referenced by both
    #     sides of the membrane so the balance closes by construction.
    jo1_1 <- vmax_oct1 / 5 * (c_eh1_p / (km_oct1 + c_eh1_p) - c_hc1 * enh / r_oct1 / (km_oct1 + c_hc1))
    jo1_2 <- vmax_oct1 / 5 * (c_eh2_p / (km_oct1 + c_eh2_p) - c_hc2 * enh / r_oct1 / (km_oct1 + c_hc2))
    jo1_3 <- vmax_oct1 / 5 * (c_eh3_p / (km_oct1 + c_eh3_p) - c_hc3 * enh / r_oct1 / (km_oct1 + c_hc3))
    jo1_4 <- vmax_oct1 / 5 * (c_eh4_p / (km_oct1 + c_eh4_p) - c_hc4 * enh / r_oct1 / (km_oct1 + c_hc4))
    jo1_5 <- vmax_oct1 / 5 * (c_eh5_p / (km_oct1 + c_eh5_p) - c_hc5 * enh / r_oct1 / (km_oct1 + c_hc5))
    jo2_1 <- vmax_oct2 / 3 * (c_pt1_p / (km_oct2 + c_pt1_p) - c_pt1_c * envpt / r_oct2 / (km_oct2 + c_pt1_c))
    jo2_2 <- vmax_oct2 / 3 * (c_pt2_p / (km_oct2 + c_pt2_p) - c_pt2_c * envpt / r_oct2 / (km_oct2 + c_pt2_c))
    jo2_3 <- vmax_oct2 / 3 * (c_pt3_p / (km_oct2 + c_pt3_p) - c_pt3_c * envpt / r_oct2 / (km_oct2 + c_pt3_c))
    jmate_1 <- vmax_mate / 3 / (km_mate + c_pt1_c) * c_pt1_c
    jmate_2 <- vmax_mate / 3 / (km_mate + c_pt2_c) * c_pt2_c
    jmate_3 <- vmax_mate / 3 / (km_mate + c_pt3_c) * c_pt3_c

    # ===============================================================
    # Absorption. ktr_gi is set so that each of the three intestinal
    # segments passes on (1 - FaFg)^(1/3) of what enters it; after three
    # segments exactly (1 - FaFg) of the dose is unabsorbed.
    # ===============================================================
    ktr_gi <- ka * (1 - fafg)^(1 / 3) / (1 - (1 - fafg)^(1 / 3))
    d/dt(transit1) <- -ktrans * transit1
    d/dt(intestine1) <- ktrans * transit1 - ka * intestine1 - ktr_gi * intestine1
    d/dt(intestine2) <- ktr_gi * intestine1 - ka * intestine2 - ktr_gi * intestine2
    d/dt(intestine3) <- ktr_gi * intestine2 - ka * intestine3 - ktr_gi * intestine3
    d/dt(a_feces) <- ktr_gi * intestine3

    # --- Systemic circulation, split into plasma and erythrocytes ---
    d/dt(plasma) <- qh_p * c_eh5_p + qr6_p * c_cd_p + kout_rbc * rbc -
      qh_p * c_plasma - qr_p * c_plasma -
      qmuscle_p * (c_plasma - c_muscle_p / kp_muscle) -
      qskin_p * (c_plasma - c_skin_p / kp_skin) -
      qadipose_p * (c_plasma - c_adipose_p / kp_adipose) -
      kin_rbc * plasma
    d/dt(rbc) <- qh_e * c_eh5_e + qr6_e * c_cd_e + kin_rbc * plasma -
      qh_e * c_rbc - qr_e * c_rbc -
      qmuscle_e * (c_rbc - c_muscle_e) -
      qskin_e * (c_rbc - c_skin_e) -
      qadipose_e * (c_rbc - c_adipose_e) -
      kout_rbc * rbc

    # --- Perfusion-limited tissues ---
    d/dt(muscle_p) <- qmuscle_p * c_plasma + kout_rbc * muscle_e -
      kin_rbc * v_muscle_p * c_muscle_p / kp_muscle - qmuscle_p * c_muscle_p / kp_muscle
    d/dt(muscle_e) <- qmuscle_e * c_rbc + kin_rbc * v_muscle_p * c_muscle_p / kp_muscle -
      kout_rbc * muscle_e - qmuscle_e * c_muscle_e
    d/dt(skin_p) <- qskin_p * c_plasma + kout_rbc * skin_e -
      kin_rbc * v_skin_p * c_skin_p / kp_skin - qskin_p * c_skin_p / kp_skin
    d/dt(skin_e) <- qskin_e * c_rbc + kin_rbc * v_skin_p * c_skin_p / kp_skin -
      kout_rbc * skin_e - qskin_e * c_skin_e
    d/dt(adipose_p) <- qadipose_p * c_plasma + kout_rbc * adipose_e -
      kin_rbc * v_adipose_p * c_adipose_p / kp_adipose - qadipose_p * c_adipose_p / kp_adipose
    d/dt(adipose_e) <- qadipose_e * c_rbc + kin_rbc * v_adipose_p * c_adipose_p / kp_adipose -
      kout_rbc * adipose_e - qadipose_e * c_adipose_e

    # --- Liver: five extracellular units in series, each with its own
    #     erythrocyte space and its own hepatocyte. Absorbed drug enters
    #     the first extracellular unit, i.e. the portal inlet.
    d/dt(is_liver1_p) <- ka * (intestine1 + intestine2 + intestine3) +
      qh_p * (c_plasma - c_eh1_p) + kout_rbc * is_liver1_e +
      ps_h_difeff / 5 * c_hc1 - jo1_1 - ps_h_difinf / 5 * c_eh1_p - kin_rbc * is_liver1_p
    d/dt(is_liver2_p) <- qh_p * (c_eh1_p - c_eh2_p) + kout_rbc * is_liver2_e +
      ps_h_difeff / 5 * c_hc2 - jo1_2 - ps_h_difinf / 5 * c_eh2_p - kin_rbc * is_liver2_p
    d/dt(is_liver3_p) <- qh_p * (c_eh2_p - c_eh3_p) + kout_rbc * is_liver3_e +
      ps_h_difeff / 5 * c_hc3 - jo1_3 - ps_h_difinf / 5 * c_eh3_p - kin_rbc * is_liver3_p
    d/dt(is_liver4_p) <- qh_p * (c_eh3_p - c_eh4_p) + kout_rbc * is_liver4_e +
      ps_h_difeff / 5 * c_hc4 - jo1_4 - ps_h_difinf / 5 * c_eh4_p - kin_rbc * is_liver4_p
    d/dt(is_liver5_p) <- qh_p * (c_eh4_p - c_eh5_p) + kout_rbc * is_liver5_e +
      ps_h_difeff / 5 * c_hc5 - jo1_5 - ps_h_difinf / 5 * c_eh5_p - kin_rbc * is_liver5_p
    d/dt(is_liver1_e) <- qh_e * (c_rbc - c_eh1_e) + kin_rbc * is_liver1_p - kout_rbc * is_liver1_e
    d/dt(is_liver2_e) <- qh_e * (c_eh1_e - c_eh2_e) + kin_rbc * is_liver2_p - kout_rbc * is_liver2_e
    d/dt(is_liver3_e) <- qh_e * (c_eh2_e - c_eh3_e) + kin_rbc * is_liver3_p - kout_rbc * is_liver3_e
    d/dt(is_liver4_e) <- qh_e * (c_eh3_e - c_eh4_e) + kin_rbc * is_liver4_p - kout_rbc * is_liver4_e
    d/dt(is_liver5_e) <- qh_e * (c_eh4_e - c_eh5_e) + kin_rbc * is_liver5_p - kout_rbc * is_liver5_e
    d/dt(int_liver1) <- jo1_1 + ps_h_difinf / 5 * c_eh1_p - ps_h_difeff / 5 * c_hc1 - cl_met / 5 * c_hc1
    d/dt(int_liver2) <- jo1_2 + ps_h_difinf / 5 * c_eh2_p - ps_h_difeff / 5 * c_hc2 - cl_met / 5 * c_hc2
    d/dt(int_liver3) <- jo1_3 + ps_h_difinf / 5 * c_eh3_p - ps_h_difeff / 5 * c_hc3 - cl_met / 5 * c_hc3
    d/dt(int_liver4) <- jo1_4 + ps_h_difinf / 5 * c_eh4_p - ps_h_difeff / 5 * c_hc4 - cl_met / 5 * c_hc4
    d/dt(int_liver5) <- jo1_5 + ps_h_difinf / 5 * c_eh5_p - ps_h_difeff / 5 * c_hc5 - cl_met / 5 * c_hc5
    d/dt(a_metab) <- cl_met / 5 * (c_hc1 + c_hc2 + c_hc3 + c_hc4 + c_hc5)

    # --- Kidney: glomerulus ---
    d/dt(glom_p) <- qr_p * c_plasma + kout_rbc * glom_e - qr1_p * c_glom_p -
      qu1 * c_glom_p - kin_rbc * glom_p
    d/dt(glom_e) <- qr_e * c_rbc + kin_rbc * glom_p - qr1_e * c_glom_e - kout_rbc * glom_e
    d/dt(glom_u) <- qgfr * (c_glom_p - c_glom_u)

    # --- Kidney: proximal tubule, three segments in series ---
    d/dt(pt1_p) <- qr1_p * c_glom_p + ps_r_pt_difeff / 3 * c_pt1_c + kout_rbc * pt1_e -
      jo2_1 - ps_r_pt_difinf / 3 * c_pt1_p - qr2_p * c_pt1_p - kin_rbc * pt1_p
    d/dt(pt2_p) <- qr2_p * c_pt1_p + ps_r_pt_difeff / 3 * c_pt2_c + kout_rbc * pt2_e -
      jo2_2 - ps_r_pt_difinf / 3 * c_pt2_p - qr3_p * c_pt2_p - kin_rbc * pt2_p
    d/dt(pt3_p) <- qr3_p * c_pt2_p + ps_r_pt_difeff / 3 * c_pt3_c + kout_rbc * pt3_e -
      jo2_3 - ps_r_pt_difinf / 3 * c_pt3_p - qr4_p * c_pt3_p - kin_rbc * pt3_p
    d/dt(pt1_e) <- qr1_e * c_glom_e + kin_rbc * pt1_p - qr2_e * c_pt1_e - kout_rbc * pt1_e
    d/dt(pt2_e) <- qr2_e * c_pt1_e + kin_rbc * pt2_p - qr3_e * c_pt2_e - kout_rbc * pt2_e
    d/dt(pt3_e) <- qr3_e * c_pt2_e + kin_rbc * pt3_p - qr4_e * c_pt3_e - kout_rbc * pt3_e
    d/dt(pt1_c) <- jo2_1 + ps_r_pt_difinf / 3 * c_pt1_p + ps_u_pt_difinf / 3 * c_pt1_u -
      ps_r_pt_difeff / 3 * c_pt1_c - jmate_1 - ps_u_pt_difeff / 3 * c_pt1_c
    d/dt(pt2_c) <- jo2_2 + ps_r_pt_difinf / 3 * c_pt2_p + ps_u_pt_difinf / 3 * c_pt2_u -
      ps_r_pt_difeff / 3 * c_pt2_c - jmate_2 - ps_u_pt_difeff / 3 * c_pt2_c
    d/dt(pt3_c) <- jo2_3 + ps_r_pt_difinf / 3 * c_pt3_p + ps_u_pt_difinf / 3 * c_pt3_u -
      ps_r_pt_difeff / 3 * c_pt3_c - jmate_3 - ps_u_pt_difeff / 3 * c_pt3_c
    d/dt(pt1_u) <- qu1 * c_glom_u + jmate_1 + ps_u_pt_difeff / 3 * c_pt1_c -
      ps_u_pt_difinf / 3 * c_pt1_u - qu2 * c_pt1_u
    d/dt(pt2_u) <- qu2 * c_pt1_u + jmate_2 + ps_u_pt_difeff / 3 * c_pt2_c -
      ps_u_pt_difinf / 3 * c_pt2_u - qu3 * c_pt2_u
    d/dt(pt3_u) <- qu3 * c_pt2_u + jmate_3 + ps_u_pt_difeff / 3 * c_pt3_c -
      ps_u_pt_difinf / 3 * c_pt3_u - qu4 * c_pt3_u

    # --- Kidney: distal tubule (passive reabsorption only) ---
    d/dt(dt_p) <- qr4_p * c_pt3_p + ps_r_dt_difeff * c_dt_c + kout_rbc * dt_e -
      ps_r_dt_difinf * c_dt_p - qr5_p * c_dt_p - kin_rbc * dt_p
    d/dt(dt_e) <- qr4_e * c_pt3_e + kin_rbc * dt_p - qr5_e * c_dt_e - kout_rbc * dt_e
    d/dt(dt_c) <- ps_r_dt_difinf * c_dt_p + ps_u_dt_difinf * c_dt_u -
      ps_r_dt_difeff * c_dt_c - ps_u_dt_difeff * c_dt_c
    d/dt(dt_u) <- qu4 * c_pt3_u + ps_u_dt_difeff * c_dt_c - ps_u_dt_difinf * c_dt_u - qu5 * c_dt_u

    # --- Kidney: collecting duct (passive reabsorption only) ---
    d/dt(cd_p) <- qr5_p * c_dt_p + ps_r_cd_difeff * c_cd_c + kout_rbc * cd_e -
      ps_r_cd_difinf * c_cd_p - qr6_p * c_cd_p - kin_rbc * cd_p
    d/dt(cd_e) <- qr5_e * c_dt_e + kin_rbc * cd_p - qr6_e * c_cd_e - kout_rbc * cd_e
    d/dt(cd_c) <- ps_r_cd_difinf * c_cd_p + ps_u_cd_difinf * c_cd_u -
      ps_r_cd_difeff * c_cd_c - ps_u_cd_difeff * c_cd_c
    d/dt(cd_u) <- qu5 * c_dt_u + ps_u_cd_difeff * c_cd_c - ps_u_cd_difinf * c_cd_u - qu6 * c_cd_u
    d/dt(urine) <- qu6 * c_cd_u

    # --- Observations. Cc is the plasma concentration the paper plots in
    #     Figure 2; Cb is the whole-blood concentration whose rise
    #     relative to plasma over the first hours is the erythrocyte
    #     distribution signature the model was built to capture.
    Cc <- c_plasma
    Cb <- (1 - ht) * c_plasma + ht * c_rbc
    Cc ~ prop(propSd)
  })
}
