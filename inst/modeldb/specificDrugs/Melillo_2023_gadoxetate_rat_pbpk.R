Melillo_2023_gadoxetate_rat_pbpk <- function() {
  description <- paste(
    "Preclinical (rat). Reduced PBPK model of the hepatobiliary MRI contrast",
    "agent gadoxetate (seven compartments, permeability-limited liver) coupled",
    "with empirical compartmental PK models of five perpetrator drugs given",
    "intravenously to the 250 g male Wistar rat: rifampicin (one-compartment),",
    "asunaprevir (two-compartment), ciclosporin (three-compartment, whole",
    "blood, with saturable red-cell binding), pioglitazone (one-compartment",
    "with a first-order oral depot) and bosentan (two-compartment). The",
    "unbound plasma concentration of each perpetrator inhibits gadoxetate",
    "active hepatic uptake (rOatp1b2 / OATP1B1) and biliary excretion (rMrp2",
    "/ MRP2) through 1 / (1 + cu / IC50) (Equation 2), with in vitro IC50",
    "values; renal clearance is not inhibited. Used for the prospective",
    "prediction of gadoxetate DCE-MRI plasma (spleen) and liver AUC ratios.",
    "Gadoxetate transporter parameters are the control-phase top-down fit of",
    "Scotcher 2021. Deterministic: no between-animal variability or residual",
    "error."
  )
  reference <- paste(
    "Melillo N, Scotcher D, Kenna JG, Green C, Hines CDG, Laitinen I,",
    "Hockings PD, Ogungbenro K, Gunwhy ER, Sourbron S, Waterton JC, Schuetz",
    "G, Galetin A. Use of In Vivo Imaging and Physiologically-Based Kinetic",
    "Modelling to Predict Hepatic Transporter Mediated Drug-Drug",
    "Interactions in Rats. Pharmaceutics. 2023;15(3):896.",
    "doi:10.3390/pharmaceutics15030896 (PMC10057977).",
    "Perpetrator PK models: Supplementary Section S1 (Eqs S1-S11, Tables",
    "S2-S6). Gadoxetate PBPK structure and parameters: Scotcher D et al.",
    "Mol Pharm. 2021;18(8):2997-3009, doi:10.1021/acs.molpharmaceut.1c00206",
    "(equation system S1, Table S2 physiology, Table 3 control-phase fit)."
  )
  vignette <- "Melillo_2023_gadoxetate_rat_pbpk"

  # Gadoxetate state names follow Scotcher_2021_gadoxetate_rat_pbpk.R:
  # `blood` = systemic blood; `vp_remainder` / `is_remainder` = rest-of-body
  # vascular / interstitial space; `is_spleen` = spleen extracellular space;
  # `is_liver` / `int_liver` = liver extracellular space / hepatocytes;
  # `urine` and `a_bile` are cumulative excretion records. The splanchnic
  # lump (stomach + gut + pancreas) has no registered organ token.
  paper_specific_compartments <- c("is_splanchnic")

  # Every dose and every state is an AMOUNT in umol, so the perpetrator
  # concentrations come out in umol/L (= uM), the units of the in vitro IC50s.
  units <- list(
    time = "h",
    dosing = "umol (gadoxetate and every perpetrator; mg/kg x 0.25 kg / MW x 1000)",
    concentration = "umol/L"
  )

  compartmentData <- list(
    blood = list(analyte = "gadoxetate", units = "umol", specimen = "whole blood", verified = TRUE),
    vp_remainder = list(analyte = "gadoxetate", units = "umol", specimen = "whole blood", verified = TRUE),
    is_remainder = list(analyte = "gadoxetate", units = "umol", specimen = "tissue", verified = TRUE),
    is_spleen = list(analyte = "gadoxetate", units = "umol", specimen = "tissue", verified = TRUE),
    is_splanchnic = list(analyte = "gadoxetate", units = "umol", specimen = "tissue", verified = TRUE),
    is_liver = list(analyte = "gadoxetate", units = "umol", specimen = "tissue", verified = TRUE),
    int_liver = list(analyte = "gadoxetate", units = "umol", specimen = "tissue", verified = TRUE),
    urine = list(analyte = "gadoxetate", units = "umol", specimen = "urine", verified = TRUE),
    a_bile = list(analyte = "gadoxetate", units = "umol", specimen = "bile", verified = TRUE),
    central_rif = list(analyte = "rifampicin", units = "umol", specimen = "plasma", verified = TRUE),
    central_asv = list(analyte = "asunaprevir", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1_asv = list(analyte = "asunaprevir", units = "umol", specimen = "tissue", verified = TRUE),
    central_ciclosporin = list(analyte = "ciclosporin", units = "umol", specimen = "whole blood", verified = TRUE),
    peripheral1_ciclosporin = list(analyte = "ciclosporin", units = "umol", specimen = "tissue", verified = TRUE),
    peripheral2_ciclosporin = list(analyte = "ciclosporin", units = "umol", specimen = "tissue", verified = TRUE),
    depot_pioglitazone = list(
      analyte = "pioglitazone",
      units = "umol",
      specimen = "administration site",
      verified = TRUE
    ),
    central_pioglitazone = list(analyte = "pioglitazone", units = "umol", specimen = "plasma", verified = TRUE),
    central_bosentan = list(analyte = "bosentan", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1_bosentan = list(analyte = "bosentan", units = "umol", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rat (Wistar)",
    n_subjects = 34L,
    n_studies = 1L,
    weight_median = "approximately 250 g at ordering (Section 2.3); the model is for the 250 g reference rat of Scotcher 2021",
    sex_female_pct = 0,
    disease_state = "healthy",
    dose_range = paste(
      "gadoxetate 25 umol/kg IV over 30 s, 30 or 60 min after vehicle",
      "(control phase) or a single IV perpetrator dose (inhibitory phase, 48",
      "h later): rifampicin 2 mg/kg, asunaprevir 5 mg/kg, bosentan 2 and",
      "4-6 mg/kg, ciclosporin 5 mg/kg, ketoconazole 3 mg/kg, pioglitazone",
      "0.4 mg/kg (Table 1)"
    ),
    regions = "multicentre preclinical imaging study (sites D, E, G1, G2 at 4.7 T and 7 T)",
    notes = paste(
      "The model was not fitted to the imaging data: it is a prospective",
      "(bottom-up) DDI prediction. Imaging cohort: 44 animals scanned in",
      "control and inhibitory phases (Table 1; n = 34 analysable profiles",
      "excluding the earlier rifampicin 10 mg/kg data). The perpetrator PK",
      "models were fitted to digitised literature rat PK data (Table S1).",
      "Ketoconazole has no PK model (its rat PK was read from the literature),",
      "so it is not represented."
    )
  )

  ini({
    # Gadoxetate: control-phase top-down fit of Scotcher 2021 Table 3, the
    # values Melillo 2023 Table S11 (CLuptake 2.17, CLbiliary 0.07 L/h) and
    # Table S10 (PS 0.172 mL/s = 0.62 L/h) cite as 'Table 3 in [1]'.
    lps_act_inf <- log(2.17)
    label("Active (OATP) hepatic uptake clearance CLuptake, uninhibited (L/h)") # Table S11 CLuptake 2.17 L/h (Scotcher 2021 Table 3)
    lcl_bile <- log(0.07)
    label("Biliary clearance from hepatocytes CLbiliary, uninhibited (L/h)") # Table S11 CLbiliary 0.07 L/h (Scotcher 2021 Table 3)
    lps_remainder <- log(0.62)
    label("Rest-of-body vascular-interstitial permeability-surface product PS (L/h)") # Table S10 PS 0.172 mL/s = 0.62 L/h (Scotcher 2021 Table 3)
    lps_dif <- fixed(log(0.014))
    label("Bidirectional passive diffusion clearance across the hepatocyte membrane CLpassive (L/h)") # Scotcher 2021 Table 3 CLpassive 0.014 (in vitro)
    fu_liver_cell <- fixed(0.648)
    label("Gadoxetate fraction unbound in hepatocytes fu,liv,cell (fraction)") # Table S11 footnote: fu in hepatocytes 0.648 (Scotcher 2021 Table 3)
    lcl_renal <- fixed(log(0.17))
    label("Renal clearance of gadoxetate from blood CLr, not inhibited (L/h)") # Scotcher 2021 Table 3 CLr 0.17
    qc <- fixed(6.62)
    label("Cardiac output QCO (L/h)") # Scotcher 2021 Table S2 cardiac output 6.62 L/h
    hct <- fixed(0.4183)
    label("Haematocrit (fraction)") # Scotcher 2021 Supporting Information Section 5: 0.4183

    # Rifampicin one-compartment model (Eq S1-S2), Table S2.
    lkel_rif <- log(0.294)
    label("Rifampicin elimination rate constant k10 (1/h)") # Table S2 k10 0.294 1/h (CV 10%)
    lvc_rif <- log(0.1821)
    label("Rifampicin central volume V1 (L)") # Table S2 V1 182.1 mL (CV 4%)

    # Asunaprevir two-compartment model (Eq S3-S5), Table S3.
    lkel_asv <- log(6.79)
    label("Asunaprevir elimination rate constant k10 (1/h)") # Table S3 k10 6.79 1/h (CV 19%)
    lk12_asv <- log(3.71)
    label("Asunaprevir central-to-peripheral rate constant k12 (1/h)") # Table S3 k12 3.71 1/h (CV 27%)
    lk21_asv <- log(1.60)
    label("Asunaprevir peripheral-to-central rate constant k21 (1/h)") # Table S3 k21 1.60 1/h (CV 34%)
    lvc_asv <- log(0.0960)
    label("Asunaprevir central volume V1 (L)") # Table S3 V1 96.0 mL (CV 17%)

    # Ciclosporin three-compartment whole-blood model (Eq S6-S9), Table S4.
    lkel_ciclosporin <- log(0.881)
    label("Ciclosporin elimination rate constant k10 (1/h)") # Table S4 k10 0.881 1/h (CV 6%)
    lk12_ciclosporin <- log(13.8)
    label("Ciclosporin central-to-first-peripheral rate constant k12 (1/h)") # Table S4 k12 13.8 1/h (CV 17%)
    lk21_ciclosporin <- log(14.6)
    label("Ciclosporin first-peripheral-to-central rate constant k21 (1/h)") # Table S4 k21 14.6 1/h (CV 16%)
    lk13_ciclosporin <- log(4.64)
    label("Ciclosporin central-to-second-peripheral rate constant k13 (1/h)") # Table S4 k13 4.64 1/h (CV 7%)
    lk31_ciclosporin <- log(0.715)
    label("Ciclosporin second-peripheral-to-central rate constant k31 (1/h)") # Table S4 K31 0.715 1/h (CV 9%)
    lvc_ciclosporin <- log(0.0635)
    label("Ciclosporin central (blood) volume V1 (L)") # Table S4 V1 63.5 mL (CV 3%)

    # Pioglitazone one-compartment model with an oral depot (Eq S1-S2),
    # Table S5 (titled 'two-compartment' but listing only ka, F, k10, V1;
    # the text and Figures S5-S6 call it one-compartment).
    lka_pioglitazone <- log(1.2)
    label("Pioglitazone oral absorption rate constant ka (1/h)") # Table S5 ka 1.2 1/h (CV 39%)
    lfdepot_pioglitazone <- log(0.66)
    label("Pioglitazone oral bioavailability F (fraction)") # Table S5 F 0.66 (CV 19%)
    lkel_pioglitazone <- log(0.17)
    label("Pioglitazone elimination rate constant k10 (1/h)") # Table S5 k10 0.17 1/h (CV 26%)
    lvc_pioglitazone <- log(0.088)
    label("Pioglitazone central volume V1 (L)") # Table S5 V1 88 mL (CV 13%)

    # Bosentan two-compartment model (Eq S3-S5), Table S6; the authors
    # selected it over the three-compartment fit of Table S7.
    lkel_bosentan <- log(5.0)
    label("Bosentan elimination rate constant k10 (1/h)") # Table S6 k10 5.0 1/h (CV 11%)
    lk12_bosentan <- log(3.8)
    label("Bosentan central-to-peripheral rate constant k12 (1/h)") # Table S6 k12 3.8 1/h (CV 45%)
    lk21_bosentan <- log(3.5)
    label("Bosentan peripheral-to-central rate constant k21 (1/h)") # Table S6 k21 3.5 1/h (CV 48%)
    lvc_bosentan <- log(0.032)
    label("Bosentan central volume V1 (L)") # Table S6 V1 32 mL (CV 13%)

    # Plasma unbound fractions of the perpetrators (Table 3 'fu' column).
    fu_rif <- fixed(0.2)
    label("Rifampicin fraction unbound in plasma (fraction)") # Table 3 fu 0.2
    fu_asv <- fixed(0.012)
    label("Asunaprevir fraction unbound in plasma (fraction)") # Table 3 fu 0.012
    fu_ciclosporin <- fixed(0.1)
    label("Ciclosporin fraction unbound in plasma (fraction)") # Table 3 fu 0.1
    fu_pioglitazone <- fixed(0.01)
    label("Pioglitazone fraction unbound in plasma (upper bound, fraction)") # Table 3 fu '<0.01'
    fu_bosentan <- fixed(0.02)
    label("Bosentan fraction unbound in plasma (upper bound, fraction)") # Table 3 fu '<0.02'

    # Ciclosporin saturable red-cell binding (Eq S11 and its legend).
    kd_ciclosporin <- fixed(0.158)
    label("Ciclosporin erythrocyte dissociation constant KD (ug/mL)") # Eq S11 legend KD 0.158 ug/mL
    lbmax_ciclosporin <- fixed(log(4.64))
    label("Ciclosporin erythrocyte binding capacity nPT (ug/mL)") # Eq S11 legend nPT 4.64 ug Eq/mL
    mw_ciclosporin <- fixed(1202.6)
    label("Ciclosporin molar mass (g/mol), converts ug/mL to umol/L") # PubChem CID 5284373 C62H111N11O12, not printed in the paper

    # In vitro inhibition constants (umol/L), the lowest rat value (Table 2)
    # or, without rat data, the lowest human value (Table 3), per Section 2.5.
    ki_oatp_rif <- fixed(0.79)
    label("Rifampicin rOatp1b2 IC50 for gadoxetate active uptake (umol/L)") # Table 2 rOatp1b2 IC50 0.79 (single-source value; reproduces Table 4 predicted AUCR)
    ki_mrp2_rif <- fixed(20)
    label("Rifampicin rMrp2 IC50 for gadoxetate biliary excretion (umol/L)") # Table 2 rMrp2 IC50 20-53; lowest 20
    ki_oatp_asv <- fixed(0.3)
    label("Asunaprevir OATP1B1 IC50 (human surrogate) for active uptake (umol/L)") # Table 3 OATP1B1 IC50 0.55 (0.3-0.79); lowest 0.3
    ki_mrp2_asv <- fixed(11)
    label("Asunaprevir rMrp2 IC50 for biliary excretion (umol/L)") # Table 2 rMrp2 IC50 11
    # Ciclosporin: the prospective prediction used the human OATP1B1 Ki (the
    # lowest reported value) rather than the rat rOatp1b2 IC50 of 1.2 -- the
    # rat value gives ~35% uptake inhibition, whereas the paper reports 97-98%
    # inhibition of CLactive (Discussion) and a liver AUCR of 0.38 (Table 4),
    # which the human Ki reproduces. Section S1.3 motivates the low value by
    # ciclosporin's strong pre-incubation effect on OATP1B. See vignette Errata.
    ki_oatp_ciclosporin <- fixed(0.014)
    label("Ciclosporin OATP1B1 Ki for active uptake (umol/L)") # Table 3 OATP1B1 Ki 0.014 (lowest; reproduces Table 4 predicted AUCR and the 97-98% uptake inhibition)
    ki_mrp2_ciclosporin <- fixed(2.7)
    label("Ciclosporin MRP2 IC50 (median) for biliary excretion (umol/L)") # Table 3 MRP2 IC50 2.7 (median; Figure S11 '<25% inhibition of CLbiliary')
    ki_oatp_pioglitazone <- fixed(5.09)
    label("Pioglitazone OATP1B1 IC50 (human surrogate) for active uptake (umol/L)") # Table 3 OATP1B1 IC50 5.09
    ki_mrp2_pioglitazone <- fixed(133)
    label("Pioglitazone MRP2 IC50 lower bound (human surrogate) for biliary excretion (umol/L)") # Table 3 MRP2 IC50 '>133'
    ki_oatp_bosentan <- fixed(5.0)
    label("Bosentan OATP1B1 IC50 (human surrogate) for active uptake (umol/L)") # Table 3 OATP1B1 IC50 6.6 (5.0-8.2); lowest 5.0
  })

  model({
    # ------------------------------------------------------------------
    # 1. Physiology of the 250 g rat (Scotcher 2021 Table S2), identical to
    #    Scotcher_2021_gadoxetate_rat_pbpk.R. Organ volume (L) = mass (g) /
    #    density (kg/L) / 1000; vascular and interstitial volumes are the
    #    organ volume times the Table S2 fractions.
    # ------------------------------------------------------------------
    v_blood <- 15.77 / 1 / 1000
    v_vas_remainder <- (0.262 * 1.25 / 1.0505 + # lungs
      0.037 * 1.44 / 1.0355 + # brain
      0.262 * 0.84 / 1.03 + # heart
      0.105 * 1.84 / 1.05 + # kidneys
      0.041 * 15 / 1.4303 + # bone
      0.026 * 101 / 1.041 + # muscle
      0.019 * 47.58 / 1.183 + # skin
      0.01 * 15.63 / 0.916) / 1000 # fat
    v_is_remainder <- (0.188 * 1.25 / 1.0505 + # lungs
      0.004 * 1.44 / 1.0355 + # brain
      0.1 * 0.84 / 1.03 + # heart
      0.2 * 1.84 / 1.05 + # kidneys
      0.1 * 15 / 1.4303 + # bone
      0.12 * 101 / 1.041 + # muscle
      0.302 * 47.58 / 1.183 + # skin
      0.135 * 15.63 / 0.916) / 1000 # fat
    v_spleen <- 0.5 / 1.054 / 1000
    vb_spleen <- 0.282 * v_spleen
    vi_spleen <- 0.15 * v_spleen
    v_is_spleen <- vb_spleen + vi_spleen
    vb_splanchnic <- (0.032 * 1.15 / 1.05 + # stomach
      0.024 * 5.6 / 1.043 + # gut
      0.18 * 0.8 / 1.045) / 1000 # pancreas
    vi_splanchnic <- (0.1 * 1.15 / 1.05 + # stomach
      0.094 * 5.6 / 1.043 + # gut
      0.12 * 0.8 / 1.045) / 1000 # pancreas
    v_is_splanchnic <- vb_splanchnic + vi_splanchnic
    v_liver <- 9.15 / 1.08 / 1000
    vb_liver <- 0.115 * v_liver
    vi_liver <- 0.163 * v_liver
    v_is_liver <- vb_liver + vi_liver
    v_int_liver <- v_liver - v_is_liver

    # Extracellular-to-blood partition coefficients (Scotcher 2021 eq S3).
    kp_spleen <- (vi_spleen + vb_spleen * (1 - hct)) / ((vi_spleen + vb_spleen) * (1 - hct))
    kp_splanchnic <- (vi_splanchnic + vb_splanchnic * (1 - hct)) /
      ((vi_splanchnic + vb_splanchnic) * (1 - hct))
    kp_liver <- (vi_liver + vb_liver * (1 - hct)) / ((vi_liver + vb_liver) * (1 - hct))

    # Blood flows (L/h, Scotcher 2021 Table S2).
    q_spleen <- 0.053
    q_splanchnic <- 0.08 + 0.85 + 0.03
    q_hepatic_artery <- 0.14
    qh <- q_hepatic_artery + q_spleen + q_splanchnic
    q_remainder <- qc - qh

    # ------------------------------------------------------------------
    # 2. Perpetrator PK (Supplementary Eqs S1-S9). Amounts in umol, volumes
    #    in L, so concentrations are umol/L.
    # ------------------------------------------------------------------
    kel_rif <- exp(lkel_rif)
    vc_rif <- exp(lvc_rif)
    kel_asv <- exp(lkel_asv)
    k12_asv <- exp(lk12_asv)
    k21_asv <- exp(lk21_asv)
    vc_asv <- exp(lvc_asv)
    kel_ciclosporin <- exp(lkel_ciclosporin)
    k12_ciclosporin <- exp(lk12_ciclosporin)
    k21_ciclosporin <- exp(lk21_ciclosporin)
    k13_ciclosporin <- exp(lk13_ciclosporin)
    k31_ciclosporin <- exp(lk31_ciclosporin)
    vc_ciclosporin <- exp(lvc_ciclosporin)
    ka_pioglitazone <- exp(lka_pioglitazone)
    kel_pioglitazone <- exp(lkel_pioglitazone)
    vc_pioglitazone <- exp(lvc_pioglitazone)
    kel_bosentan <- exp(lkel_bosentan)
    k12_bosentan <- exp(lk12_bosentan)
    k21_bosentan <- exp(lk21_bosentan)
    vc_bosentan <- exp(lvc_bosentan)

    Cc_rif <- central_rif / vc_rif
    Cc_asv <- central_asv / vc_asv
    Cc_ciclosporin <- central_ciclosporin / vc_ciclosporin # whole blood
    Cc_pioglitazone <- central_pioglitazone / vc_pioglitazone
    Cc_bosentan <- central_bosentan / vc_bosentan

    # Unbound plasma concentrations (umol/L) driving inhibition (Eq S10).
    cu_rif <- fu_rif * Cc_rif
    cu_asv <- fu_asv * Cc_asv
    cu_pioglitazone <- fu_pioglitazone * Cc_pioglitazone
    cu_bosentan <- fu_bosentan * Cc_bosentan

    # Ciclosporin: unbound concentration from the whole-blood concentration
    # with saturable erythrocyte binding (Eq S11, Tanaka 2000), evaluated in
    # ug/mL as printed and converted to umol/L. The coefficient
    # (1 - Hct) * fu + Hct is transcribed as printed in Eq S11.
    bmax_ciclosporin <- exp(lbmax_ciclosporin)
    cb_csa_ugml <- Cc_ciclosporin * mw_ciclosporin / 1000
    a_csa <- (1 - hct) * fu_ciclosporin + hct
    b_csa <- cb_csa_ugml - a_csa * kd_ciclosporin - hct * bmax_ciclosporin
    cu_csa_ugml <- 0.5 * (b_csa + sqrt(b_csa^2 + 4 * a_csa * kd_ciclosporin * cb_csa_ugml)) / a_csa
    cu_ciclosporin <- cu_csa_ugml / mw_ciclosporin * 1000

    # ------------------------------------------------------------------
    # 3. Transporter inhibition (Equation 2). Each perpetrator was studied
    #    alone, where the sum below reduces exactly to 1 / (1 + cu / IC50).
    #    Bosentan does not inhibit rMrp2 (Table 2 'NI'), so it has no
    #    biliary term. Renal clearance is not inhibited (Section 2.5).
    # ------------------------------------------------------------------
    inh_uptake <- 1 + cu_rif / ki_oatp_rif + cu_asv / ki_oatp_asv +
      cu_ciclosporin / ki_oatp_ciclosporin + cu_pioglitazone / ki_oatp_pioglitazone +
      cu_bosentan / ki_oatp_bosentan
    inh_bile <- 1 + cu_rif / ki_mrp2_rif + cu_asv / ki_mrp2_asv +
      cu_ciclosporin / ki_mrp2_ciclosporin + cu_pioglitazone / ki_mrp2_pioglitazone
    f_uptake_remaining <- 1 / inh_uptake
    f_bile_remaining <- 1 / inh_bile

    ps_act_inf <- exp(lps_act_inf) / inh_uptake
    cl_bile <- exp(lcl_bile) / inh_bile
    ps_remainder <- exp(lps_remainder)
    ps_dif <- exp(lps_dif)
    cl_renal <- exp(lcl_renal)

    # ------------------------------------------------------------------
    # 4. Gadoxetate concentrations (umol/L).
    # ------------------------------------------------------------------
    c_blood <- blood / v_blood
    c_vas_remainder <- vp_remainder / v_vas_remainder
    c_is_remainder <- is_remainder / v_is_remainder
    c_is_spleen <- is_spleen / v_is_spleen
    c_is_splanchnic <- is_splanchnic / v_is_splanchnic
    c_is_liver <- is_liver / v_is_liver
    c_int_liver <- int_liver / v_int_liver

    # ------------------------------------------------------------------
    # 5. ODEs. Gadoxetate: Scotcher 2021 equation system S1 with the liver
    #    equations of Equation (1). Perpetrators: Eqs S1-S9.
    # ------------------------------------------------------------------
    d/dt(blood) <- -qc * c_blood - cl_renal * c_blood + q_remainder * c_vas_remainder +
      qh * c_is_liver / kp_liver
    d/dt(vp_remainder) <- q_remainder * (c_blood - c_vas_remainder) -
      ps_remainder * (c_vas_remainder - c_is_remainder)
    d/dt(is_remainder) <- ps_remainder * (c_vas_remainder - c_is_remainder)
    d/dt(is_spleen) <- q_spleen * (c_blood - c_is_spleen / kp_spleen)
    d/dt(is_splanchnic) <- q_splanchnic * (c_blood - c_is_splanchnic / kp_splanchnic)
    d/dt(is_liver) <- (qh - q_spleen - q_splanchnic) * c_blood +
      q_spleen * c_is_spleen / kp_spleen +
      q_splanchnic * c_is_splanchnic / kp_splanchnic -
      qh * c_is_liver / kp_liver -
      ps_act_inf * c_is_liver -
      ps_dif * (c_is_liver - fu_liver_cell * c_int_liver)
    d/dt(int_liver) <- ps_act_inf * c_is_liver +
      ps_dif * (c_is_liver - fu_liver_cell * c_int_liver) -
      cl_bile * fu_liver_cell * c_int_liver
    d/dt(urine) <- cl_renal * c_blood
    d/dt(a_bile) <- cl_bile * fu_liver_cell * c_int_liver

    d/dt(central_rif) <- -kel_rif * central_rif
    d/dt(central_asv) <- -(kel_asv + k12_asv) * central_asv + k21_asv * peripheral1_asv
    d/dt(peripheral1_asv) <- k12_asv * central_asv - k21_asv * peripheral1_asv
    d/dt(central_ciclosporin) <- -(kel_ciclosporin + k12_ciclosporin + k13_ciclosporin) * central_ciclosporin +
      k21_ciclosporin * peripheral1_ciclosporin + k31_ciclosporin * peripheral2_ciclosporin
    d/dt(peripheral1_ciclosporin) <- k12_ciclosporin * central_ciclosporin -
      k21_ciclosporin * peripheral1_ciclosporin
    # Eq S8 is printed 'k13 A3 - k31 A3'; the inflow is k13 * A1.
    d/dt(peripheral2_ciclosporin) <- k13_ciclosporin * central_ciclosporin -
      k31_ciclosporin * peripheral2_ciclosporin
    # Eq S1 is printed with 'ka * A1' in the central equation; the inflow is
    # ka * AG (the depot), as in Eq S3.
    d/dt(depot_pioglitazone) <- -ka_pioglitazone * depot_pioglitazone
    d/dt(central_pioglitazone) <- ka_pioglitazone * depot_pioglitazone -
      kel_pioglitazone * central_pioglitazone
    f(depot_pioglitazone) <- exp(lfdepot_pioglitazone)
    d/dt(central_bosentan) <- -(kel_bosentan + k12_bosentan) * central_bosentan +
      k21_bosentan * peripheral1_bosentan
    d/dt(peripheral1_bosentan) <- k12_bosentan * central_bosentan - k21_bosentan * peripheral1_bosentan

    # ------------------------------------------------------------------
    # 6. Gadoxetate observations. Cc is the blood concentration (umol/L).
    #    Delta R1 (1/s) uses the ex vivo relaxivities (1/s per mM) of
    #    Scotcher 2021 Table 1: blood 6.4 / 6.2 and hepatocytes 7.6 / 6 at
    #    4.7 / 7 T; spleen and liver extracellular space take the blood value.
    # ------------------------------------------------------------------
    Cc <- c_blood
    C_spleen <- c_is_spleen * v_is_spleen / v_spleen
    C_liver <- (is_liver + int_liver) / v_liver
    r1_blood_4p7t <- 6.4
    r1_hep_4p7t <- 7.6
    r1_blood_7t <- 6.2
    r1_hep_7t <- 6
    dR1_spleen_4p7t <- c_is_spleen / 1000 * v_is_spleen * r1_blood_4p7t / v_spleen
    dR1_liver_4p7t <- (c_is_liver * v_is_liver * r1_blood_4p7t + c_int_liver * v_int_liver * r1_hep_4p7t) /
      (1000 * v_liver)
    dR1_spleen_7t <- c_is_spleen / 1000 * v_is_spleen * r1_blood_7t / v_spleen
    dR1_liver_7t <- (c_is_liver * v_is_liver * r1_blood_7t + c_int_liver * v_int_liver * r1_hep_7t) /
      (1000 * v_liver)
  })
}
