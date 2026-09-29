Nishiyama_2019_cimetidine_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, transporter-mediated, segmented nephron).",
    "Cimetidine after a single oral dose in healthy adults, built as the",
    "perpetrator arm of a metformin drug-drug-interaction analysis. The",
    "structure is the companion metformin model of the same paper with",
    "three simplifications the authors state explicitly: the gastric",
    "transit compartment is removed and the three intestinal segments are",
    "collapsed into one with an absorption lag time; the five hepatic",
    "units are collapsed into a single extracellular / hepatocyte pair",
    "with LINEAR (not Michaelis-Menten) uptake; and the plasma and",
    "erythrocyte spaces are merged into one blood compartment with a",
    "blood-to-plasma ratio. The kidney keeps the full nephron -- a",
    "glomerulus, three proximal-tubule segments, a distal tubule and a",
    "collecting duct, each with vascular, cellular and urinary-lumen",
    "spaces -- with basolateral OCT2 and OAT3 uptake, luminal MATE",
    "efflux, and passive diffusion resolved into ionised and unionised",
    "components. Unlike metformin, cimetidine (pKa 6.9) is only partly",
    "ionised at physiological pH, so every passive-permeability and",
    "Nernst term is split by the Henderson-Hasselbalch fractions at the",
    "local pH, and the ionised species diffuses at a fraction lambda of",
    "the unionised rate. No parameter was fitted: all values are taken",
    "from published sources and used as reported, so the paper reports no",
    "between-subject variability and no residual-error model and propSd",
    "is a placeholder."
  )
  reference <- paste(
    "Nishiyama K, Toshimoto K, Lee W, Ishiguro N, Bister B, Sugiyama Y.",
    "Physiologically-Based Pharmacokinetic Modeling Analysis for",
    "Quantitative Prediction of Renal Transporter-Mediated Interactions",
    "Between Metformin and Cimetidine.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8(6):396-406.",
    "doi:10.1002/psp4.12398.",
    "The ODE system is transcribed from the 'Model equations for",
    "cimetidine' section of Supplementary Material S2",
    "(PSP4-8-396-s008.pdf) and the 'Development of the cimetidine PBPK",
    "model' section of the Supplemental Text (PSP4-8-396-s007.pdf).",
    "Drug parameters are Table S4 (PSP4-8-396-s006.pdf), body physiology",
    "is Table S2 (PSP4-8-396-s004.pdf) and kidney physiology is Table S3",
    "(PSP4-8-396-s005.pdf). The transporter Km and Vmax values are",
    "themselves taken by the authors from Burt HJ et al.,",
    "Eur J Pharm Sci. 2016;88:70-82, scaled by a relative activity factor.",
    "See the vignette Errata for the transcription corrections applied to",
    "the published equation list.",
    sep = " "
  )
  vignette <- "Nishiyama_2019_metformin_cimetidine"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  covariateData <- list()

  # Compartment vocabulary. `urine` is canonical
  # (inst/references/compartment-names.md). `is_liver` (hepatic
  # extracellular / sinusoidal space) and `int_liver` (hepatocyte space)
  # are the unsegmented forms of the `is_liver<n>` / `int_liver<n>` stems
  # used by the registered models of this same Sugiyama-laboratory
  # framework (Toshimoto_2017_irinotecan_pbpk.R,
  # Tsuchitani_2024_telmisartan_pbpk.R, Aoki_2024_bosentan_pbpk.R). The
  # nephron stems match the companion Nishiyama_2019_metformin_pbpk.R:
  # `glom_*` glomerulus, `pt1..3_*` proximal-tubule segments, `dt_*`
  # distal tubule, `cd_*` collecting duct, with `_b` (vascular blood),
  # `_c` (tubular cell) and `_u` (urinary lumen) spaces. Metformin's
  # `_p` / `_e` plasma / erythrocyte split does not appear here because
  # the cimetidine model merges them into one blood space.
  paper_specific_compartments <- c(
    "intestine",
    "a_feces",
    "a_metab",
    "blood",
    "muscle",
    "skin",
    "adipose",
    "is_liver",
    "int_liver",
    "glom_b",
    "glom_u",
    "pt1_b",
    "pt2_b",
    "pt3_b",
    "pt1_c",
    "pt2_c",
    "pt3_c",
    "pt1_u",
    "pt2_u",
    "pt3_u",
    "dt_b",
    "dt_c",
    "dt_u",
    "cd_b",
    "cd_c",
    "cd_u"
  )

  compartmentData <- list(
    intestine = list(analyte = "cimetidine", units = "ug", specimen = "administration site", verified = TRUE),
    a_feces = list(analyte = "cimetidine", units = "ug", specimen = "faeces", verified = TRUE),
    blood = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    muscle = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    is_liver = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    int_liver = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    a_metab = list(analyte = "cimetidine", units = "ug", specimen = "not applicable", verified = TRUE),
    glom_b = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    glom_u = list(analyte = "cimetidine", units = "ug", specimen = "urine", verified = TRUE),
    pt1_b = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    pt2_b = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    pt3_b = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    pt1_c = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    pt2_c = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    pt3_c = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    pt1_u = list(analyte = "cimetidine", units = "ug", specimen = "urine", verified = TRUE),
    pt2_u = list(analyte = "cimetidine", units = "ug", specimen = "urine", verified = TRUE),
    pt3_u = list(analyte = "cimetidine", units = "ug", specimen = "urine", verified = TRUE),
    dt_b = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    dt_c = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    dt_u = list(analyte = "cimetidine", units = "ug", specimen = "urine", verified = TRUE),
    cd_b = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    cd_c = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    cd_u = list(analyte = "cimetidine", units = "ug", specimen = "urine", verified = TRUE),
    urine = list(analyte = "cimetidine", units = "ug", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "healthy adults",
    dose_range = "400 mg cimetidine, single oral dose",
    notes = paste(
      "The simulation was compared against the published 400 mg oral",
      "single-dose profile of Grahnen et al. 1979",
      "(Eur J Clin Pharmacol 16:335-340). Body physiology is the standard",
      "70 kg adult of Davies & Morris 1993 (Table S2), so subject counts,",
      "demographics and covariate distributions do not enter the model."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Physicochemical and absorption -- Table S4. None of the cimetidine
    # parameters was fitted in this paper: Table S4 cites Burt 2016 for
    # the compound layer, so every value is held constant.
    # ---------------------------------------------------------------
    pka <- fixed(6.9)
    label("Acid dissociation constant (unitless)") # Table S4: pKa 6.9
    lam <- fixed(0.1)
    label("Passive diffusion of the ionised relative to the unionised species (unitless)") # Supplemental Text: lambda set to 0.1 after Yoshikado 2017
    lka <- fixed(log(0.70))
    label("Intestinal absorption rate constant (1/h)") # Table S4: ka 0.70 /h
    ltlag <- fixed(log(0.15))
    label("Intestinal absorption lag time (h)") # Table S4: Tlag 0.15 h
    fafg <- fixed(0.92)
    label("Intestinal availability Fa x Fg (fraction)") # Table S4: FaFg 0.92

    # ---------------------------------------------------------------
    # Distribution -- Table S4. The Kp values are referenced to BLOOD,
    # not plasma, because the cimetidine model carries a single blood
    # compartment.
    # ---------------------------------------------------------------
    fu <- fixed(0.80)
    label("Unbound fraction in plasma (fraction)") # Table S4: fu,cim 0.80
    rb <- fixed(0.97)
    label("Blood-to-plasma concentration ratio (unitless)") # Table S4: Rb 0.97
    kp_muscle <- fixed(0.86)
    label("Muscle-to-blood concentration ratio (unitless)") # Table S4: Kp,muscle 0.86
    kp_skin <- fixed(0.72)
    label("Skin-to-blood concentration ratio (unitless)") # Table S4: Kp,skin 0.72
    kp_adipose <- fixed(0.24)
    label("Adipose-to-blood concentration ratio (unitless)") # Table S4: Kp,adipose 0.24

    # ---------------------------------------------------------------
    # Liver -- Table S4. Uptake is linear, so there is no hepatic Km.
    # ---------------------------------------------------------------
    ps_act <- fixed(12.0)
    label("Hepatic active uptake clearance (L/h)") # Table S4: PSact 12.0 L/h
    rdif <- fixed(1.16)
    label("Hepatic passive-to-active uptake clearance ratio (unitless)") # Table S4: Rdif 1.16
    cl_met_h <- fixed(11.3)
    label("Hepatic metabolic intrinsic clearance (L/h)") # Table S4: CLmet 11.3 L/h

    # ---------------------------------------------------------------
    # Kidney -- Table S4. The Vmax values already carry the relative
    # activity factor of Burt 2016, so they are whole-organ maxima.
    # ---------------------------------------------------------------
    pd <- fixed(7.9e-5)
    label("Passive permeability from the PAMPA assay (m/h)") # Table S4: Pd 7.9e-5 m/h
    r_oct2 <- fixed(1.32)
    label("OCT2 influx-to-efflux clearance ratio (unitless)") # Table S4: ROCT2,inf/eff 1.32
    km_oct2_um <- fixed(72.6)
    label("OCT2 Michaelis constant for cimetidine (umol/L)") # Table S4: Km,OCT2 72.6 umol/L
    vmax_oct2_umol <- fixed(7265)
    label("OCT2 maximum transport rate (umol/h)") # Table S4: Vmax,OCT2 7,265 umol/h
    km_oat3_um <- fixed(161)
    label("OAT3 Michaelis constant for cimetidine (umol/L)") # Table S4: Km,OAT3 161 umol/L
    vmax_oat3_umol <- fixed(4124)
    label("OAT3 maximum transport rate (umol/h)") # Table S4: Vmax,OAT3 4,124 umol/h
    km_mate_um <- fixed(7.7)
    label("MATE Michaelis constant for cimetidine (umol/L)") # Table S4: Km,MATE 7.7 umol/L
    vmax_mate_umol <- fixed(453)
    label("MATE maximum transport rate (umol/h)") # Table S4: Vmax,MATE 453 umol/h

    # ---------------------------------------------------------------
    # Body physiology -- Table S2 (shared with the metformin model).
    # ---------------------------------------------------------------
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
    # Kidney physiology -- Table S3 (shared with the metformin model).
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
    volt_h <- fixed(-40)
    label("Hepatocyte plasma-membrane potential (mV)") # Table S1: membrane potential -40 mV
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
    ph_b <- fixed(7.4)
    label("Vascular pH (unitless)") # Table S3: vessel pH 7.4
    ph_c <- fixed(7.2)
    label("Tubular-cell pH (unitless)") # Table S3: cell pH 7.2
    ph_u_pt <- fixed(7.0)
    label("Proximal-tubule luminal pH (unitless)") # Table S3: proximal tubule lumen pH 7.0
    ph_u_dt <- fixed(6.7)
    label("Distal-tubule luminal pH (unitless)") # Table S3: distal tubule lumen pH 6.7
    ph_u_cd <- fixed(6.4)
    label("Collecting-duct luminal pH (unitless)") # Table S3: collecting duct lumen pH 6.4

    propSd <- fixed(0.1)
    label("Proportional residual error (fraction)") # not reported; placeholder
  })

  model({
    fara <- 96485 # Faraday constant, C/mol
    rgas <- 8.314 # gas constant, J/(mol K)
    tabs <- 310 # body temperature, K
    zval <- 1 # cimetidine valence when protonated
    mw <- 252.34 # cimetidine free base, g/mol; converts Km and Vmax from umol to ug

    ka <- exp(lka)
    tlag <- exp(ltlag)

    # --- Henderson-Hasselbalch ionised / unionised fractions at each pH
    #     (Suppl. S2, cimetidine 'Other calculations').
    fr_ion <- 1 / (1 + 10^(ph_b - pka))
    fr_union <- fr_ion * 10^(ph_b - pka)
    fc_ion <- 1 / (1 + 10^(ph_c - pka))
    fc_union <- fc_ion * 10^(ph_c - pka)
    fu_ion_pt <- 1 / (1 + 10^(ph_u_pt - pka))
    fu_union_pt <- fu_ion_pt * 10^(ph_u_pt - pka)
    fu_ion_dt <- 1 / (1 + 10^(ph_u_dt - pka))
    fu_union_dt <- fu_ion_dt * 10^(ph_u_dt - pka)
    fu_ion_cd <- 1 / (1 + 10^(ph_u_cd - pka))
    fu_union_cd <- fu_ion_cd * 10^(ph_u_cd - pka)

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

    # --- Liver: linear uptake, Suppl. Text eqs. (5) and (7) ---
    gamma_h <- (fr_union + lam * fr_ion) / (fc_union + enh * lam * fc_ion)
    ps_h_difinf <- rdif * ps_act
    ps_h_difeff <- ps_h_difinf / gamma_h

    # --- Kidney passive diffusion. Pd is measured on the mixture present
    #     at blood pH, so it is first resolved into a unionised
    #     permeability and an ionised permeability lambda times smaller;
    #     each membrane then recombines them at the local pH, with the
    #     ionised component carrying the Goldman potential factor.
    pd_union <- pd / (fr_ion * lam + fr_union)
    pd_ion <- pd_union * lam
    ps_r_pt_difinf <- pd_union * sa_r_pt * 1000 * fr_union +
      pd_ion * sa_r_pt * 1000 * nvpt / (envpt - 1) * fr_ion
    gamma_r_pt <- (lam * fr_ion + fr_union) / (envpt * lam * fc_ion + fc_union)
    ps_r_pt_difeff <- ps_r_pt_difinf / gamma_r_pt
    ps_u_pt_difinf <- pd_union * sa_u_pt * 1000 * fu_union_pt +
      pd_ion * sa_u_pt * 1000 * nupt / (enupt - 1) * fu_ion_pt
    gamma_u_pt <- (lam * fu_ion_pt + fu_union_pt) / (enupt * lam * fc_ion + fc_union)
    ps_u_pt_difeff <- ps_u_pt_difinf / gamma_u_pt
    ps_r_dt_difinf <- pd_union * sa_r_dt * 1000 * fr_union +
      pd_ion * sa_r_dt * 1000 * nvdt / (envdt - 1) * fr_ion
    gamma_r_dt <- (lam * fr_ion + fr_union) / (envdt * lam * fc_ion + fc_union)
    ps_r_dt_difeff <- ps_r_dt_difinf / gamma_r_dt
    ps_u_dt_difinf <- pd_union * sa_u_dt * 1000 * fu_union_dt +
      pd_ion * sa_u_dt * 1000 * nudt / (enudt - 1) * fu_ion_dt
    gamma_u_dt <- (lam * fu_ion_dt + fu_union_dt) / (enudt * lam * fc_ion + fc_union)
    ps_u_dt_difeff <- ps_u_dt_difinf / gamma_u_dt
    ps_r_cd_difinf <- pd_union * sa_r_cd * 1000 * fr_union +
      pd_ion * sa_r_cd * 1000 * nvcd / (envcd - 1) * fr_ion
    gamma_r_cd <- (lam * fr_ion + fr_union) / (envcd * lam * fc_ion + fc_union)
    ps_r_cd_difeff <- ps_r_cd_difinf / gamma_r_cd
    ps_u_cd_difinf <- pd_union * sa_u_cd * 1000 * fu_union_cd +
      pd_ion * sa_u_cd * 1000 * nucd / (enucd - 1) * fu_ion_cd
    gamma_u_cd <- (lam * fu_ion_cd + fu_union_cd) / (enucd * lam * fc_ion + fc_union)
    ps_u_cd_difeff <- ps_u_cd_difinf / gamma_u_cd

    # --- Transporter Km and Vmax converted from umol to ug. The
    #     published conversion line additionally multiplies both by the
    #     unbound and ionised fractions; because that factor cancels in
    #     Vmax / Km it changes only where the carrier saturates, and as
    #     printed it is dimensionally inconsistent with the unbound,
    #     ionised driving concentrations used below. The in vitro Km is
    #     therefore used on its own scale. See vignette Errata.
    km_oct2 <- km_oct2_um * mw
    vmax_oct2 <- vmax_oct2_umol * mw
    km_oat3 <- km_oat3_um * mw
    vmax_oat3 <- vmax_oat3_umol * mw
    km_mate <- km_mate_um * mw
    vmax_mate <- vmax_mate_umol * mw

    qu1 <- qgfr
    dqu <- (qgfr - qurine) / 5
    qu2 <- qu1 - dqu
    qu3 <- qu2 - dqu
    qu4 <- qu3 - dqu
    qu5 <- qu4 - dqu
    qu6 <- qurine
    qr1 <- qr - qu1
    qr2 <- qr1 + dqu
    qr3 <- qr2 + dqu
    qr4 <- qr3 + dqu
    qr5 <- qr4 + dqu
    qr6 <- qr - qu6

    vr_b <- 0.30 * vkidney
    vr_cell <- 0.24 * vkidney
    vr_u <- 0.46 * vkidney
    v_glom_b <- 5 / 129 * vr_b
    v_pt_b <- 17 / 129 * vr_b
    v_dt_b <- 11 / 129 * vr_b
    v_cd_b <- 62 / 129 * vr_b
    v_glom_u <- 5 / 129 * vr_u
    v_pt_u <- 17 / 129 * vr_u
    v_dt_u <- 11 / 129 * vr_u
    v_cd_u <- 62 / 129 * vr_u
    v_pt_c <- 17 / 124 * vr_cell
    v_dt_c <- 11 / 124 * vr_cell
    v_cd_c <- 62 / 124 * vr_cell

    c_blood <- blood / vblood
    c_muscle <- muscle / vmuscle
    c_skin <- skin / vskin
    c_adipose <- adipose / vadipose
    c_eh <- is_liver / veh
    c_hc <- int_liver / vhc
    c_glom_b <- glom_b / v_glom_b
    c_glom_u <- glom_u / v_glom_u
    c_pt1_b <- pt1_b / (v_pt_b / 3)
    c_pt2_b <- pt2_b / (v_pt_b / 3)
    c_pt3_b <- pt3_b / (v_pt_b / 3)
    c_pt1_c <- pt1_c / (v_pt_c / 3)
    c_pt2_c <- pt2_c / (v_pt_c / 3)
    c_pt3_c <- pt3_c / (v_pt_c / 3)
    c_pt1_u <- pt1_u / (v_pt_u / 3)
    c_pt2_u <- pt2_u / (v_pt_u / 3)
    c_pt3_u <- pt3_u / (v_pt_u / 3)
    c_dt_b <- dt_b / v_dt_b
    c_dt_c <- dt_c / v_dt_c
    c_dt_u <- dt_u / v_dt_u
    c_cd_b <- cd_b / v_cd_b
    c_cd_c <- cd_c / v_cd_c
    c_cd_u <- cd_u / v_cd_u

    # --- Unbound, ionised concentrations drive the renal carriers ---
    cub1 <- fu * fr_ion * c_pt1_b / rb
    cub2 <- fu * fr_ion * c_pt2_b / rb
    cub3 <- fu * fr_ion * c_pt3_b / rb
    cuc1 <- fc_ion * c_pt1_c
    cuc2 <- fc_ion * c_pt2_c
    cuc3 <- fc_ion * c_pt3_c
    jo2_1 <- vmax_oct2 / 3 * (cub1 / (km_oct2 + cub1) - cuc1 * envpt / r_oct2 / (km_oct2 + cuc1))
    jo2_2 <- vmax_oct2 / 3 * (cub2 / (km_oct2 + cub2) - cuc2 * envpt / r_oct2 / (km_oct2 + cuc2))
    jo2_3 <- vmax_oct2 / 3 * (cub3 / (km_oct2 + cub3) - cuc3 * envpt / r_oct2 / (km_oct2 + cuc3))
    joat_1 <- vmax_oat3 / 3 * cub1 / (km_oat3 + cub1)
    joat_2 <- vmax_oat3 / 3 * cub2 / (km_oat3 + cub2)
    joat_3 <- vmax_oat3 / 3 * cub3 / (km_oat3 + cub3)
    jmate_1 <- vmax_mate / 3 * cuc1 / (km_mate + cuc1)
    jmate_2 <- vmax_mate / 3 * cuc2 / (km_mate + cuc2)
    jmate_3 <- vmax_mate / 3 * cuc3 / (km_mate + cuc3)

    # --- Absorption: one intestinal compartment emptying at ka / FaFg,
    #     of which ka is absorbed and the remainder lost to faeces.
    d/dt(intestine) <- -ka / fafg * intestine
    alag(intestine) <- tlag
    d/dt(a_feces) <- ka * (1 - fafg) / fafg * intestine

    d/dt(blood) <- qh * c_eh + qr6 * c_cd_b - qh * c_blood - qr * c_blood -
      qmuscle * (c_blood - c_muscle / kp_muscle) -
      qskin * (c_blood - c_skin / kp_skin) -
      qadipose * (c_blood - c_adipose / kp_adipose)
    d/dt(muscle) <- qmuscle * (c_blood - c_muscle / kp_muscle)
    d/dt(skin) <- qskin * (c_blood - c_skin / kp_skin)
    d/dt(adipose) <- qadipose * (c_blood - c_adipose / kp_adipose)

    d/dt(is_liver) <- ka * intestine + qh * (c_blood - c_eh) + ps_h_difeff * c_hc -
      fu * (ps_h_difinf + ps_act) * c_eh / rb
    d/dt(int_liver) <- fu * (ps_h_difinf + ps_act) * c_eh / rb -
      ps_h_difeff * c_hc - cl_met_h * c_hc
    d/dt(a_metab) <- cl_met_h * c_hc

    d/dt(glom_b) <- qr * c_blood - qr1 * c_glom_b - fu * qu1 * c_glom_b / rb
    d/dt(glom_u) <- qgfr * (fu * c_glom_b / rb - c_glom_u)

    d/dt(pt1_b) <- qr1 * c_glom_b + ps_r_pt_difeff / 3 * c_pt1_c - jo2_1 - joat_1 -
      ps_r_pt_difinf / 3 * fu * c_pt1_b / rb - qr2 * c_pt1_b
    d/dt(pt2_b) <- qr2 * c_pt1_b + ps_r_pt_difeff / 3 * c_pt2_c - jo2_2 - joat_2 -
      ps_r_pt_difinf / 3 * fu * c_pt2_b / rb - qr3 * c_pt2_b
    d/dt(pt3_b) <- qr3 * c_pt2_b + ps_r_pt_difeff / 3 * c_pt3_c - jo2_3 - joat_3 -
      ps_r_pt_difinf / 3 * fu * c_pt3_b / rb - qr4 * c_pt3_b
    d/dt(pt1_c) <- jo2_1 + joat_1 + ps_r_pt_difinf / 3 * fu * c_pt1_b / rb +
      ps_u_pt_difinf / 3 * c_pt1_u - ps_r_pt_difeff / 3 * c_pt1_c -
      jmate_1 - ps_u_pt_difeff / 3 * c_pt1_c
    d/dt(pt2_c) <- jo2_2 + joat_2 + ps_r_pt_difinf / 3 * fu * c_pt2_b / rb +
      ps_u_pt_difinf / 3 * c_pt2_u - ps_r_pt_difeff / 3 * c_pt2_c -
      jmate_2 - ps_u_pt_difeff / 3 * c_pt2_c
    d/dt(pt3_c) <- jo2_3 + joat_3 + ps_r_pt_difinf / 3 * fu * c_pt3_b / rb +
      ps_u_pt_difinf / 3 * c_pt3_u - ps_r_pt_difeff / 3 * c_pt3_c -
      jmate_3 - ps_u_pt_difeff / 3 * c_pt3_c
    d/dt(pt1_u) <- qu1 * c_glom_u + jmate_1 + ps_u_pt_difeff / 3 * c_pt1_c -
      ps_u_pt_difinf / 3 * c_pt1_u - qu2 * c_pt1_u
    d/dt(pt2_u) <- qu2 * c_pt1_u + jmate_2 + ps_u_pt_difeff / 3 * c_pt2_c -
      ps_u_pt_difinf / 3 * c_pt2_u - qu3 * c_pt2_u
    d/dt(pt3_u) <- qu3 * c_pt2_u + jmate_3 + ps_u_pt_difeff / 3 * c_pt3_c -
      ps_u_pt_difinf / 3 * c_pt3_u - qu4 * c_pt3_u

    d/dt(dt_b) <- qr4 * c_pt3_b + ps_r_dt_difeff * c_dt_c -
      ps_r_dt_difinf * fu * c_dt_b / rb - qr5 * c_dt_b
    d/dt(dt_c) <- ps_r_dt_difinf * fu * c_dt_b / rb + ps_u_dt_difinf * c_dt_u -
      ps_r_dt_difeff * c_dt_c - ps_u_dt_difeff * c_dt_c
    d/dt(dt_u) <- qu4 * c_pt3_u + ps_u_dt_difeff * c_dt_c - ps_u_dt_difinf * c_dt_u - qu5 * c_dt_u

    d/dt(cd_b) <- qr5 * c_dt_b + ps_r_cd_difeff * c_cd_c -
      ps_r_cd_difinf * fu * c_cd_b / rb - qr6 * c_cd_b
    d/dt(cd_c) <- ps_r_cd_difinf * fu * c_cd_b / rb + ps_u_cd_difinf * c_cd_u -
      ps_r_cd_difeff * c_cd_c - ps_u_cd_difeff * c_cd_c
    d/dt(cd_u) <- qu5 * c_dt_u + ps_u_cd_difeff * c_cd_c - ps_u_cd_difinf * c_cd_u - qu6 * c_cd_u
    d/dt(urine) <- qu6 * c_cd_u

    # --- Observations. Cc is plasma, which is what Figure 2i plots;
    #     Cb is the whole-blood concentration the ODEs carry.
    Cc <- c_blood / rb
    Cb <- c_blood
    Cc ~ prop(propSd)
  })
}
