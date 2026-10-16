Karrer_2018_bisphenolS_pbpk <- function() {
  description <- "PBPK (whole-body, hand-coded in R/deSolve). Human eight-compartment physiologically based pharmacokinetic model for bisphenol S (BPS) and its glucuronide and sulfate, after peroral and dermal exposure. Same structure as the BPA model (Karrer_2018_bisphenolA_pbpk): perfusion-limited gonads, liver, fat, skin, brain, richly- and slowly-perfused lumps around a serum pool, presystemic gut and hepatic glucuronidation, hepatic sulfation, conjugate serum submodules, an enterohepatic-recirculation loop (67% of hepatic glucuronide), and two parallel dermal skin-surface depots. BPS uses SUBSTRATE-INHIBITION kinetics (paper Equation 4) for GUT glucuronidation and Michaelis-Menten for hepatic glucuronidation. The model was calibrated against the Oh et al. (2018) BPS biomonitoring study (increased peroral uptake and EHR). Transcribed from the deSolve source deposited in the Supplemental Material. Deterministic typical-value (no IIV, no residual error). Companion models: Karrer_2018_bisphenolA_pbpk, Karrer_2018_bisphenolF_pbpk, Karrer_2018_bisphenolAF_pbpk."
  reference <- paste(
    "Karrer C, Roiss T, von Goetz N, Gramec Skledar D, Peterlin Masic L, Hungerbuhler K. (2018).",
    "Physiologically Based Pharmacokinetic (PBPK) Modeling of the Bisphenols BPA, BPS, BPF, and BPAF",
    "with New Experimental Metabolic Parameters: Comparing the Pharmacokinetic Behavior of BPA with",
    "Its Substitutes. Environmental Health Perspectives 126(7):077002.",
    "doi:10.1289/EHP2739. PMCID PMC6108829.",
    sep = " "
  )
  vignette <- "Karrer_2018_bisphenols_pbpk"

  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nmol/L",
    amount = "nmol",
    weight = "kg"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Table 4 fixes the reference adult female at 60 kg. Scales every tissue volume (linear),",
        "the conjugate volumes of distribution (linear, = plasma volume), the saturable-metabolism",
        "vmax terms (liver vmax via liver mass; gut vmax linear in body weight), the first-order",
        "rate constants (gastric emptying, oral uptake, enterocyte-to-serum transport, EHR) by",
        "bw^-0.25, and the urinary clearances by bw^0.75. Cardiac output is NOT body-weight scaled:",
        "it is supplied per age group (see qcc). To simulate another age group, override WT, qcc and",
        "the fv_* / fq_* fractions with that group's row from Table S1.",
        sep = " "
      ),
      source_name = "bw (kg)"
    )
  )

  compartmentData <- list(
    stomach = list(analyte = "Bisphenol S", units = "nmol", specimen = "administration site", verified = TRUE),
    a_small_intestine = list(analyte = "Bisphenol S", units = "nmol", specimen = "tissue", verified = TRUE),
    a_liver = list(analyte = "Bisphenol S", units = "nmol", specimen = "tissue", verified = TRUE),
    a_fat = list(analyte = "Bisphenol S", units = "nmol", specimen = "tissue", verified = TRUE),
    a_gonads = list(analyte = "Bisphenol S", units = "nmol", specimen = "tissue", verified = TRUE),
    a_skin = list(analyte = "Bisphenol S", units = "nmol", specimen = "tissue", verified = TRUE),
    a_brain = list(analyte = "Bisphenol S", units = "nmol", specimen = "tissue", verified = TRUE),
    a_rapidly_perfused = list(analyte = "Bisphenol S", units = "nmol", specimen = "tissue", verified = TRUE),
    a_slowly_perfused = list(analyte = "Bisphenol S", units = "nmol", specimen = "tissue", verified = TRUE),
    plasma = list(analyte = "Bisphenol S", units = "nmol", specimen = "serum", verified = TRUE),
    urine = list(analyte = "Bisphenol S", units = "nmol", specimen = "urine", verified = TRUE),
    a_small_intestine_gluc = list(analyte = "BPS glucuronide", units = "nmol", specimen = "tissue", verified = TRUE),
    gut_lumen_gluc = list(analyte = "BPS glucuronide", units = "nmol", specimen = "bile", verified = TRUE),
    serum_gluc = list(analyte = "BPS glucuronide", units = "nmol", specimen = "serum", verified = TRUE),
    urine_gluc = list(analyte = "BPS glucuronide", units = "nmol", specimen = "urine", verified = TRUE),
    a_small_intestine_sulf = list(analyte = "BPS sulfate", units = "nmol", specimen = "tissue", verified = TRUE),
    gut_lumen_sulf = list(analyte = "BPS sulfate", units = "nmol", specimen = "bile", verified = TRUE),
    serum_sulf = list(analyte = "BPS sulfate", units = "nmol", specimen = "serum", verified = TRUE),
    urine_sulf = list(analyte = "BPS sulfate", units = "nmol", specimen = "urine", verified = TRUE),
    depot_td_thermalpaper = list(
      analyte = "Bisphenol S",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_td_personalcare = list(
      analyte = "Bisphenol S",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 0L,
    n_studies = 0L,
    age_range = "reference adult female, 30 years (Table 4); structure also covers the other Table S1 age groups",
    weight_range = "60 kg (reference adult female, Table 4)",
    sex_female_pct = 100,
    race_ethnicity = NA_character_,
    disease_state = "healthy reference individual",
    dose_range = "single peroral or dermal 500 ng/kg BW in the age-comparison scenarios; the Oh et al. (2018) calibration study used 8.75 ug deuterated BPS/kg BW oral in 7 adults",
    regions = "Switzerland (ETH Zurich); physiology Edginton et al. (2006) / ICRP (2002)",
    notes = paste(
      "Physiologically based forward model (n_subjects = 0). BPS was calibrated against the Oh et al.",
      "(2018) biomonitoring data (7 adults, single oral 8.75 ug/kg): the peroral uptake rate and the",
      "glucuronide clearance were increased and an EHR share of 67% introduced to match the measured",
      "BPS and BPS-g serum profiles. Partition coefficients from the Zhang and Zhang (2006) QSAR",
      "(Table 6, 14% deprotonated); glucuronidation kinetics measured in human liver and intestinal",
      "microsomes (Table 7, substrate-inhibition form in the gut, Michaelis-Menten in the liver);",
      "hepatic sulfation approximated by the BPA values. Parent molecular weight (249.88 g/mol) is",
      "derived from the Discussion's BPS-glucuronide weight (426 g/mol) minus the glucuronyl increment.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Reference physiology -- adult female, 60 kg (Table 4; fractions Table
    # S1 / C3). Chemical-independent; identical across the four bisphenols.
    # =====================================================================
    qcc <- fixed(5.90); label("Table S1 / C3, adult-female cardiac output (L/min); QC = qcc*60 L/h")

    fq_gonads <- fixed(0.000178); label("Table S1 / C3, gonad blood flow (fraction of cardiac output)")
    fq_liver <- fixed(0.240); label("Table S1 / C3, liver blood flow (fraction of cardiac output)")
    fq_fat <- fixed(0.0745); label("Table S1 / C3, fat blood flow (fraction of cardiac output)")
    fq_brain <- fixed(0.105); label("Table S1 / C3, brain blood flow (fraction of cardiac output)")
    fq_skin <- fixed(0.0437); label("Table S1 / C3, skin blood flow (fraction of cardiac output)")
    fq_slowly_perfused <- fixed(0.142); label("Table S1 / C3, slowly-perfused blood flow (fraction of cardiac output)")

    fv_plasma <- fixed(0.040); label("Table S1 / C3, plasma volume (fraction of body weight)")
    fv_fat <- fixed(0.317); label("Table S1 / C3, fat volume (fraction of body weight)")
    fv_liver <- fixed(0.0233); label("Table S1 / C3, liver volume (fraction of body weight)")
    fv_brain <- fixed(0.0217); label("Table S1 / C3, brain volume (fraction of body weight)")
    fv_skin <- fixed(0.0383); label("Table S1 / C3, skin volume (fraction of body weight)")
    fv_gonads <- fixed(0.000183); label("Table S1 / C3, gonad volume (fraction of body weight)")
    fv_slowly_perfused <- fixed(0.430); label("Table S1 / C3, slowly-perfused volume (fraction of body weight)")
    fv_rapidly_perfused <- fixed(0.0594); label("Table S1 / C3, richly-perfused volume (fraction of body weight)")

    # =====================================================================
    # BPS physicochemistry. Molecular weight 249.88 g/mol (BPS-glucuronide
    # 426 g/mol minus the glucuronyl increment 176.12); used for dose
    # conversion in the vignette, not in the ODEs.
    #
    # Tissue:serum partition coefficients -- Table 6 (Zhang and Zhang 2006 QSAR).
    # =====================================================================
    pc_liver <- fixed(0.846); label("Table 6, BPS liver:serum")
    pc_slowly_perfused <- fixed(0.881); label("Table 6, BPS slowly-perfused:serum")
    pc_brain <- fixed(0.810); label("Table 6, BPS brain:serum")
    pc_rapidly_perfused <- fixed(0.810); label("Table 6, BPS richly-perfused:serum")
    pc_fat <- fixed(0.435); label("Table 6, BPS fat:serum")
    pc_skin <- fixed(0.915); label("Table 6, BPS skin:serum")
    pc_gonads <- fixed(0.843); label("Table 6, BPS gonads:serum")

    # =====================================================================
    # Metabolism -- Table 7 ('used in basic model'). BPS uses SUBSTRATE
    # INHIBITION (Equation 4) for gut glucuronidation and Michaelis-Menten
    # for hepatic glucuronidation (no hepatic Ksi reported).
    # =====================================================================
    km_liver_gluc <- fixed(446000); label("Table 7, BPS hepatic glucuronidation Km (nM, visually fitted)")
    vmax_liver_gluc <- fixed(7810); label("Table 7, BPS hepatic glucuronidation vmax (nmol/h/g liver, visually fitted)")
    km_gut_gluc <- fixed(555000); label("Table 7, BPS enterocyte glucuronidation Km (nM, visually fitted)")
    k_substrateinh_gut_gluc <- fixed(711000); label("Table 7, BPS enterocyte glucuronidation substrate-inhibition constant Ksi (nM)")
    vmax_gut_gluc <- fixed(563); label("Table 7, BPS enterocyte glucuronidation vmax (nmol/h/kg bw)")
    km_liver_sulf <- fixed(10100); label("Table 7, hepatic sulfation Km (nM), approximated by the BPA value")
    vmax_liver_sulf <- fixed(149); label("Table 7, hepatic sulfation vmax (nmol/h/g liver), approximated by the BPA value")

    # =====================================================================
    # Uptake, EHR and clearance -- Table 8. BPS calibration increased oral
    # uptake and urinary clearance and introduced a 67% EHR share.
    # =====================================================================
    ge_c <- fixed(3.50); label("Table 8, gastric emptying coefficient (1/h/kgBW^-0.25)")
    k1_c <- fixed(5.00); label("Table 8, BPS oral uptake small-intestine->liver coefficient (1/h/kgBW^-0.25, calibrated)")
    kgiing_c <- fixed(50.0); label("Table 8, glucuronide enterocyte->serum transport coefficient (1/h/kgBW^-0.25)")
    kgiins_c <- fixed(50.0); label("Table 8, sulfate enterocyte->serum transport coefficient (1/h/kgBW^-0.25)")
    ehrrate_c <- fixed(0.20); label("Table 8, glucuronide reuptake-from-gut (EHR) coefficient (1/h/kgBW^-0.25)")
    kenterobpag_c <- fixed(0.20); label("Table 8, glucuronide biliary-deconjugation-to-parent coefficient (1/h/kgBW^-0.25)")
    kenterobpas_c <- fixed(0.0); label("deposited code, sulfate biliary deconjugation coefficient (none)")
    met1g <- fixed(0.33); label("Table 8, fraction of hepatic BPS-g taken up directly to serum (67% EHR)")
    met1s <- fixed(1.0); label("deposited code, fraction of hepatic BPS-s taken up directly to serum (no sulfate EHR)")
    cl_urine_bpa <- fixed(0.30); label("Table 8, BPS urinary clearance coefficient (L/h/kgBW^0.75, calibrated)")
    cl_urine_gluc <- fixed(0.35); label("Table 8, glucuronide urinary clearance coefficient (L/h/kgBW^0.75)")
    cl_urine_sulf <- fixed(0.030); label("Table 8, sulfate urinary clearance coefficient (L/h/kgBW^0.75)")
    enterocytes <- fixed(0.1223); label("Table 8, enterocyte distribution volume (L) (Gertz 2011)")

    # =====================================================================
    # Route-specific uptake -- Table 5 / deposited code (oral period = 3 min).
    # =====================================================================
    eoa_tp <- fixed(0.20); label("Table 5, dermal extent of absorption, thermal paper (fraction)")
    eoa_pcp <- fixed(0.60); label("Table 5, dermal extent of absorption, personal care products (fraction)")
    ahl_tp <- fixed(6.0); label("Table 5, dermal penetration half-life, thermal paper (h)")
    ahl_pcp <- fixed(0.167); label("Table 5, dermal penetration half-life, personal care products (h)")
    period_oral <- fixed(0.05); label("deposited code, peroral zero-order uptake period (h) (= 3/60)")
    period_dermal <- fixed(24.0); label("Table 5, dermal zero-order uptake period (h)")
  })

  model({
    # -------------------------------------------------------------------
    # 1. Cardiac output and regional blood flows (L/h).
    # -------------------------------------------------------------------
    qc <- qcc * 60

    q_gonad <- fq_gonads * qc
    q_liver <- fq_liver * qc
    q_fat <- fq_fat * qc
    q_brain <- fq_brain * qc
    q_skin <- fq_skin * qc
    q_slow <- fq_slowly_perfused * qc
    q_rich <- qc - q_liver - q_brain - q_fat - q_gonad - q_skin - q_slow

    # -------------------------------------------------------------------
    # 2. Compartment volumes (L).
    # -------------------------------------------------------------------
    v_plasma <- fv_plasma * WT
    v_fat <- fv_fat * WT
    v_liver <- fv_liver * WT
    v_brain <- fv_brain * WT
    v_skin <- fv_skin * WT
    v_gonad <- fv_gonads * WT
    v_slow <- fv_slowly_perfused * WT
    v_rich <- fv_rapidly_perfused * WT
    v_bodyg <- v_plasma
    v_bodys <- v_plasma

    # -------------------------------------------------------------------
    # 3. Allometric scaling.
    # -------------------------------------------------------------------
    vmax_liver_gluc_tot <- vmax_liver_gluc * v_liver * 1000
    vmax_liver_sulf_tot <- vmax_liver_sulf * v_liver * 1000
    vmax_gut_gluc_tot <- vmax_gut_gluc * WT

    ge <- ge_c / WT^0.25
    k1 <- k1_c / WT^0.25
    kgiing <- kgiing_c / WT^0.25
    kgiins <- kgiins_c / WT^0.25
    ehrrate <- ehrrate_c / WT^0.25
    kenterobpag <- kenterobpag_c / WT^0.25
    kenterobpas <- kenterobpas_c / WT^0.25

    cl_bpa <- cl_urine_bpa * WT^0.75
    cl_gluc <- cl_urine_gluc * WT^0.75
    cl_sulf <- cl_urine_sulf * WT^0.75

    met2g <- 1.0 - met1g
    met2s <- 1.0 - met1s

    # -------------------------------------------------------------------
    # 4. Dermal penetration rate constants (1/h).
    # -------------------------------------------------------------------
    kpen_tp <- log(2) / ahl_tp
    kpen_pcp <- log(2) / ahl_pcp
    rate_tp <- kpen_tp * depot_td_thermalpaper
    rate_pcp <- kpen_pcp * depot_td_personalcare

    # -------------------------------------------------------------------
    # 5. Concentrations (nmol/L).
    # -------------------------------------------------------------------
    c_art <- plasma / v_plasma
    c_gut <- a_small_intestine / enterocytes

    cv_liver <- a_liver / (v_liver * pc_liver)
    cv_fat <- a_fat / (v_fat * pc_fat)
    cv_gonad <- a_gonads / (v_gonad * pc_gonads)
    cv_skin <- a_skin / (v_skin * pc_skin)
    cv_brain <- a_brain / (v_brain * pc_brain)
    cv_rich <- a_rapidly_perfused / (v_rich * pc_rapidly_perfused)
    cv_slow <- a_slowly_perfused / (v_slow * pc_slowly_perfused)

    c_ven <- (cv_liver * q_liver + cv_skin * q_skin + cv_fat * q_fat +
      cv_rich * q_rich + cv_slow * q_slow + cv_gonad * q_gonad +
      cv_brain * q_brain) / qc

    c_gluc <- serum_gluc / v_bodyg
    c_sulf <- serum_sulf / v_bodys

    # -------------------------------------------------------------------
    # 6. Metabolic rates (nmol/h). SUBSTRATE INHIBITION in the gut
    #    (Equation 4); Michaelis-Menten in the liver.
    # -------------------------------------------------------------------
    r_gut_gluc <- vmax_gut_gluc_tot * c_gut /
      (km_gut_gluc + c_gut * (1 + c_gut / k_substrateinh_gut_gluc))
    r_liver_gluc <- vmax_liver_gluc_tot * cv_liver / (km_liver_gluc + cv_liver)
    r_liver_sulf <- vmax_liver_sulf_tot * cv_liver / (km_liver_sulf + cv_liver)

    r_gi_gluc <- kgiing * a_small_intestine_gluc
    r_gi_sulf <- kgiins * a_small_intestine_sulf
    r_oral <- k1 * a_small_intestine

    # -------------------------------------------------------------------
    # 7. Parent BPS mass balance (nmol/h).
    # -------------------------------------------------------------------
    d/dt(stomach) <- -ge * stomach
    d/dt(a_small_intestine) <- ge * stomach - r_gut_gluc - r_oral

    d/dt(a_liver) <- q_liver * (c_art - cv_liver) + r_oral -
      r_liver_gluc - r_liver_sulf +
      gut_lumen_gluc * kenterobpag + gut_lumen_sulf * kenterobpas

    d/dt(a_fat) <- q_fat * (c_art - cv_fat)
    d/dt(a_gonads) <- q_gonad * (c_art - cv_gonad)
    d/dt(a_skin) <- rate_tp + rate_pcp + q_skin * (c_art - cv_skin)
    d/dt(a_brain) <- q_brain * (c_art - cv_brain)
    d/dt(a_rapidly_perfused) <- q_rich * (c_art - cv_rich)
    d/dt(a_slowly_perfused) <- q_slow * (c_art - cv_slow)

    d/dt(plasma) <- qc * (c_ven - c_art) - cl_bpa * c_ven
    d/dt(urine) <- cl_bpa * c_ven

    # -------------------------------------------------------------------
    # 8. BPS-glucuronide.
    # -------------------------------------------------------------------
    d/dt(a_small_intestine_gluc) <- r_gut_gluc - r_gi_gluc
    d/dt(gut_lumen_gluc) <- met2g * (r_liver_gluc + r_gi_gluc) -
      gut_lumen_gluc * (ehrrate + kenterobpag)
    d/dt(serum_gluc) <- met1g * (r_liver_gluc + r_gi_gluc) +
      gut_lumen_gluc * ehrrate - cl_gluc * c_gluc
    d/dt(urine_gluc) <- cl_gluc * c_gluc

    # -------------------------------------------------------------------
    # 9. BPS-sulfate (met1s = 1; gut sulfate pool carried inert).
    # -------------------------------------------------------------------
    d/dt(a_small_intestine_sulf) <- -r_gi_sulf
    d/dt(gut_lumen_sulf) <- met2s * (r_liver_sulf + r_gi_sulf) -
      gut_lumen_sulf * (ehrrate + kenterobpas)
    d/dt(serum_sulf) <- met1s * (r_liver_sulf + r_gi_sulf) +
      gut_lumen_sulf * ehrrate - cl_sulf * c_sulf
    d/dt(urine_sulf) <- cl_sulf * c_sulf

    # -------------------------------------------------------------------
    # 10. Skin-surface depots.
    # -------------------------------------------------------------------
    d/dt(depot_td_thermalpaper) <- -rate_tp
    d/dt(depot_td_personalcare) <- -rate_pcp
    f(depot_td_thermalpaper) <- eoa_tp
    f(depot_td_personalcare) <- eoa_pcp
    dur(depot_td_thermalpaper) <- period_dermal
    dur(depot_td_personalcare) <- period_dermal
    dur(stomach) <- period_oral

    # -------------------------------------------------------------------
    # 11. Observations (nmol/L = nM).
    # -------------------------------------------------------------------
    Cc <- c_art
    Cgluc <- c_gluc
    Csulf <- c_sulf
    Cgonad <- a_gonads / v_gonad
  })
}
