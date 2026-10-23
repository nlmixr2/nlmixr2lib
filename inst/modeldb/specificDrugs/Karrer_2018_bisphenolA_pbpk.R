Karrer_2018_bisphenolA_pbpk <- function() {
  description <- "PBPK (whole-body, hand-coded in R/deSolve). Human eight-compartment physiologically based pharmacokinetic model for bisphenol A (BPA) and its glucuronide (BPA-g) and sulfate (BPA-s), after peroral and dermal exposure. Perfusion-limited gonads, liver, fat, skin, brain, richly-perfused and slowly-perfused tissue lumps around a well-mixed serum pool; presystemic glucuronidation and sulfation in the gut (enterocytes) and liver (Michaelis-Menten for BPA); conjugate serum submodules with urinary clearance; and an enterohepatic-recirculation (EHR) loop (10% of hepatic BPA-g). Dermal uptake enters the skin organ through two parallel skin-surface depots (thermal paper, personal care products). Transcribed from the deSolve source code deposited in the paper's Supplemental Material; reproduces the paper's published BPA predictions. Deterministic typical-value (no IIV, no residual error). Companion models: Karrer_2018_bisphenolS_pbpk, Karrer_2018_bisphenolF_pbpk, Karrer_2018_bisphenolAF_pbpk."
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
        "the fv_* / fq_* fractions with that group's row from Table S1 (adult male, adolescents,",
        "children, toddlers, infants).",
        sep = " "
      ),
      source_name = "bw (kg)"
    )
  )

  compartmentData <- list(
    # --- BPA (parent), amounts in nmol ---
    stomach = list(analyte = "Bisphenol A", units = "nmol", specimen = "administration site", verified = TRUE),
    a_small_intestine = list(analyte = "Bisphenol A", units = "nmol", specimen = "tissue", verified = TRUE),
    a_liver = list(analyte = "Bisphenol A", units = "nmol", specimen = "tissue", verified = TRUE),
    a_fat = list(analyte = "Bisphenol A", units = "nmol", specimen = "tissue", verified = TRUE),
    a_gonads = list(analyte = "Bisphenol A", units = "nmol", specimen = "tissue", verified = TRUE),
    a_skin = list(analyte = "Bisphenol A", units = "nmol", specimen = "tissue", verified = TRUE),
    a_brain = list(analyte = "Bisphenol A", units = "nmol", specimen = "tissue", verified = TRUE),
    a_rapidly_perfused = list(analyte = "Bisphenol A", units = "nmol", specimen = "tissue", verified = TRUE),
    a_slowly_perfused = list(analyte = "Bisphenol A", units = "nmol", specimen = "tissue", verified = TRUE),
    plasma = list(analyte = "Bisphenol A", units = "nmol", specimen = "serum", verified = TRUE),
    urine = list(analyte = "Bisphenol A", units = "nmol", specimen = "urine", verified = TRUE),
    # --- BPA-glucuronide, amounts in nmol ---
    a_small_intestine_gluc = list(analyte = "BPA glucuronide", units = "nmol", specimen = "tissue", verified = TRUE),
    gut_lumen_gluc = list(analyte = "BPA glucuronide", units = "nmol", specimen = "bile", verified = TRUE),
    serum_gluc = list(analyte = "BPA glucuronide", units = "nmol", specimen = "serum", verified = TRUE),
    urine_gluc = list(analyte = "BPA glucuronide", units = "nmol", specimen = "urine", verified = TRUE),
    # --- BPA-sulfate, amounts in nmol ---
    a_small_intestine_sulf = list(analyte = "BPA sulfate", units = "nmol", specimen = "tissue", verified = TRUE),
    gut_lumen_sulf = list(analyte = "BPA sulfate", units = "nmol", specimen = "bile", verified = TRUE),
    serum_sulf = list(analyte = "BPA sulfate", units = "nmol", specimen = "serum", verified = TRUE),
    urine_sulf = list(analyte = "BPA sulfate", units = "nmol", specimen = "urine", verified = TRUE),
    # --- dermal skin-surface depots, amounts in nmol ---
    depot_td_thermalpaper = list(
      analyte = "Bisphenol A",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_td_personalcare = list(
      analyte = "Bisphenol A",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 0L,
    n_studies = 0L,
    age_range = "reference adult female, 30 years (Table 4); the model structure also covers infants (6 d - 3 months), toddlers (1-3 y), children (3-10 y), adolescents (10-18 y) and adult males via the age-group physiology of Table S1",
    weight_range = "60 kg (reference adult female, Table 4)",
    sex_female_pct = 100,
    race_ethnicity = NA_character_,
    disease_state = "healthy reference individual",
    dose_range = "single peroral or dermal 500 ng/kg BW in the age-comparison scenarios; the Teeguarden et al. (2015) evaluation used 30 ug deuterated BPA/kg BW oral",
    regions = "Switzerland (ETH Zurich); physiology Edginton et al. (2006) / ICRP (2002)",
    notes = paste(
      "This is a physiologically based forward model, not a population fit: n_subjects = 0.",
      "The BPA model was recalibrated from Yang et al. (2015) against the Thayer et al. (2015)",
      "biomonitoring data (a scaling bug in the published gut-glucuronidation vmax was fixed) and",
      "evaluated against Teeguarden et al. (2015) and Hormann et al. (2014). Reference physiology is",
      "the adult female (60 kg, Table 4 / Table S1); the fv_* volume fractions, fq_* flow fractions",
      "and absolute cardiac output qcc are held at that group's values and are meant to be overridden",
      "per age group. Metabolism parameters are from Table 7 (hepatic glucuronidation Coughlin et al.",
      "2012; intestinal glucuronidation Trdan Lusin et al. 2012; hepatic sulfation Kurebayashi et al.",
      "2010). The deposited code sets the conjugate volume of distribution equal to the plasma volume",
      "(VbodygC <- VplasmaC), which is followed here rather than Table 8's 0.0435 L/kg.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Reference physiology -- adult female, 60 kg (Table 4; fractions from
    # Table S1 / Table C3). Cardiac output qcc is an ABSOLUTE value per age
    # group (L/min), NOT body-weight scaled. The volume fractions fv_* are
    # fractions of body weight; the flow fractions fq_* are fractions of
    # cardiac output; the richly-perfused flow is the remainder.
    # =====================================================================
    qcc <- fixed(5.90); label("Table S1 / Table C3, adult-female cardiac output (L/min); QC = qcc*60 L/h")

    fq_gonads <- fixed(0.000178); label("Table S1 / C3, adult-female gonad blood flow (fraction of cardiac output)")
    fq_liver <- fixed(0.240); label("Table S1 / C3, liver blood flow (fraction of cardiac output)")
    fq_fat <- fixed(0.0745); label("Table S1 / C3, fat blood flow (fraction of cardiac output)")
    fq_brain <- fixed(0.105); label("Table S1 / C3, brain blood flow (fraction of cardiac output)")
    fq_skin <- fixed(0.0437); label("Table S1 / C3, skin blood flow (fraction of cardiac output)")
    fq_slowly_perfused <- fixed(0.142); label("Table S1 / C3, 'muscle' (slowly perfused) blood flow (fraction of cardiac output)")

    fv_plasma <- fixed(0.040); label("Table S1 / C3, plasma volume (fraction of body weight)")
    fv_fat <- fixed(0.317); label("Table S1 / C3, fat volume (fraction of body weight)")
    fv_liver <- fixed(0.0233); label("Table S1 / C3, liver volume (fraction of body weight)")
    fv_brain <- fixed(0.0217); label("Table S1 / C3, brain volume (fraction of body weight)")
    fv_skin <- fixed(0.0383); label("Table S1 / C3, skin volume (fraction of body weight)")
    fv_gonads <- fixed(0.000183); label("Table S1 / C3, gonad volume (fraction of body weight)")
    fv_slowly_perfused <- fixed(0.430); label("Table S1 / C3, 'muscle' (slowly perfused) volume (fraction of body weight)")
    fv_rapidly_perfused <- fixed(0.0594); label("Table S1 / C3, richly-perfused volume (fraction of body weight)")

    # =====================================================================
    # BPA physicochemistry. Molecular weight 228.28 g/mol (deposited code) is
    # used only to convert an external ng/kg dose to the nmol the model doses
    # in; it is applied in the vignette, not in the ODEs, so it is not an
    # ini() parameter.
    #
    # Tissue:serum partition coefficients -- Table 6 (Doerge et al. 2011 for
    # BPA; Zhang and Zhang 2006 QSAR for skin).
    pc_liver <- fixed(0.730); label("Table 6, BPA liver:serum")
    pc_slowly_perfused <- fixed(2.70); label("Table 6, BPA slowly-perfused:serum")
    pc_brain <- fixed(2.80); label("Table 6, BPA brain:serum")
    pc_rapidly_perfused <- fixed(2.80); label("Table 6, BPA richly-perfused:serum")
    pc_fat <- fixed(5.00); label("Table 6, BPA fat:serum")
    pc_skin <- fixed(2.15); label("Table 6, BPA skin:serum")
    pc_gonads <- fixed(2.60); label("Table 6, BPA gonads:serum")

    # =====================================================================
    # Metabolism -- Table 7 ('used in basic model' column). BPA uses plain
    # Michaelis-Menten in both gut and liver (no substrate-inhibition Ksi).
    # vmax_liver_* are per gram of liver; vmax_gut_gluc is per kg body weight.
    # =====================================================================
    km_liver_gluc <- fixed(45800); label("Table 7, BPA hepatic glucuronidation Km (nM) (Coughlin et al. 2012)")
    vmax_liver_gluc <- fixed(9040); label("Table 7, BPA hepatic glucuronidation vmax (nmol/h/g liver)")
    km_gut_gluc <- fixed(58400); label("Table 7, BPA enterocyte glucuronidation Km (nM) (Trdan Lusin et al. 2012)")
    vmax_gut_gluc <- fixed(361); label("Table 7, BPA enterocyte glucuronidation vmax (nmol/h/kg bw)")
    km_liver_sulf <- fixed(10100); label("Table 7, BPA hepatic sulfation Km (nM) (Kurebayashi et al. 2010)")
    vmax_liver_sulf <- fixed(149); label("Table 7, BPA hepatic sulfation vmax (nmol/h/g liver)")

    # =====================================================================
    # Uptake, EHR and clearance -- Table 8. '_c' marks an allometric
    # coefficient scaled inside model(): rate constants by bw^-0.25,
    # clearances by bw^0.75.
    # =====================================================================
    ge_c <- fixed(3.50); label("Table 8, gastric emptying coefficient (1/h/kgBW^-0.25)")
    k1_c <- fixed(2.10); label("Table 8, BPA oral uptake small-intestine->liver coefficient (1/h/kgBW^-0.25)")
    kgiing_c <- fixed(50.0); label("Table 8, BPA-g enterocyte->serum transport coefficient (1/h/kgBW^-0.25)")
    kgiins_c <- fixed(50.0); label("Table 8, BPA-s enterocyte->serum transport coefficient (1/h/kgBW^-0.25)")
    ehrrate_c <- fixed(0.20); label("Table 8, BPA-g reuptake-from-gut (EHR) coefficient (1/h/kgBW^-0.25)")
    kenterobpag_c <- fixed(0.20); label("Table 8, BPA-g biliary-deconjugation-to-parent coefficient (1/h/kgBW^-0.25)")
    kenterobpas_c <- fixed(0.0); label("deposited code, BPA-s biliary deconjugation coefficient (none)")
    met1g <- fixed(0.9); label("Table 8, fraction of hepatic BPA-g taken up directly to serum (10% EHR)")
    met1s <- fixed(1.0); label("deposited code, fraction of hepatic BPA-s taken up directly to serum (no sulfate EHR)")
    cl_urine_bpa <- fixed(0.060); label("Table 8, BPA urinary clearance coefficient (L/h/kgBW^0.75)")
    cl_urine_gluc <- fixed(0.35); label("Table 8, BPA-g urinary clearance coefficient (L/h/kgBW^0.75)")
    cl_urine_sulf <- fixed(0.030); label("Table 8, BPA-s urinary clearance coefficient (L/h/kgBW^0.75)")
    enterocytes <- fixed(0.1223); label("Table 8, enterocyte distribution volume (L) (Gertz 2011)")

    # =====================================================================
    # Route-specific uptake -- Table 5 and deposited code. Peroral extent of
    # absorption is 100% (applied as the dose amount, so no f() is needed on
    # the stomach). The oral uptake period is the deposited code value of
    # 3 min (= 0.05 h), NOT the 0.25 h printed in Table 5: 3 min reproduces
    # the published tmax, 0.25 h does not (see the vignette Errata).
    # =====================================================================
    eoa_tp <- fixed(0.20); label("Table 5, dermal extent of absorption, thermal paper (fraction) (Toner et al. 2016)")
    eoa_pcp <- fixed(0.60); label("Table 5, dermal extent of absorption, personal care products (fraction) (Biedermann et al. 2010)")
    ahl_tp <- fixed(6.0); label("Table 5, dermal penetration half-life, thermal paper (h) (Demierre et al. 2012)")
    ahl_pcp <- fixed(0.167); label("Table 5, dermal penetration half-life, personal care products (h) (Biedermann et al. 2010)")
    period_oral <- fixed(0.05); label("deposited code, peroral zero-order uptake period (h) (= 3/60)")
    period_dermal <- fixed(24.0); label("Table 5, dermal zero-order uptake period (h)")
  })

  model({
    # -------------------------------------------------------------------
    # 1. Cardiac output and regional blood flows (L/h). Cardiac output is
    #    absolute (qcc in L/min); only the tissue volumes scale with WT.
    # -------------------------------------------------------------------
    qc <- qcc * 60

    q_gonad <- fq_gonads * qc
    q_liver <- fq_liver * qc
    q_fat <- fq_fat * qc
    q_brain <- fq_brain * qc
    q_skin <- fq_skin * qc
    q_slow <- fq_slowly_perfused * qc
    # Richly-perfused flow is the remainder of cardiac output, matching the
    # deposited code (Qrich <- QC - Qliver - Qbrain - Qfat - Qgonad - Qskin - Qslow).
    q_rich <- qc - q_liver - q_brain - q_fat - q_gonad - q_skin - q_slow

    # -------------------------------------------------------------------
    # 2. Compartment volumes (L). Tissue fractions of body weight; the
    #    conjugate distribution volumes equal the plasma volume.
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
    # 3. Allometric scaling of metabolism, uptake and clearance.
    #    Liver vmax: nmol/h/g liver * liver mass (v_liver L * 1000 g/L).
    #    Gut vmax: linear in body weight. Rate constants / bw^0.25,
    #    clearances * bw^0.75 -- deposited code.
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
    # 4. Dermal penetration rate constants (1/h) from the skin-surface
    #    depots into the skin organ.
    # -------------------------------------------------------------------
    kpen_tp <- log(2) / ahl_tp
    kpen_pcp <- log(2) / ahl_pcp
    rate_tp <- kpen_tp * depot_td_thermalpaper
    rate_pcp <- kpen_pcp * depot_td_personalcare

    # -------------------------------------------------------------------
    # 5. Concentrations (nmol/L). c_* are tissue concentrations; cv_* are
    #    the venous-blood concentrations leaving each tissue (c / kp).
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

    # Mixed venous return concentration (flow-weighted over the tissues).
    c_ven <- (cv_liver * q_liver + cv_skin * q_skin + cv_fat * q_fat +
      cv_rich * q_rich + cv_slow * q_slow + cv_gonad * q_gonad +
      cv_brain * q_brain) / qc

    c_gluc <- serum_gluc / v_bodyg
    c_sulf <- serum_sulf / v_bodys

    # -------------------------------------------------------------------
    # 6. Metabolic rates (nmol/h). BPA: plain Michaelis-Menten in both gut
    #    and liver.
    # -------------------------------------------------------------------
    r_gut_gluc <- vmax_gut_gluc_tot * c_gut / (km_gut_gluc + c_gut)
    r_liver_gluc <- vmax_liver_gluc_tot * cv_liver / (km_liver_gluc + cv_liver)
    r_liver_sulf <- vmax_liver_sulf_tot * cv_liver / (km_liver_sulf + cv_liver)

    # Enterocyte-to-serum transport of the conjugates formed in the gut.
    r_gi_gluc <- kgiing * a_small_intestine_gluc
    r_gi_sulf <- kgiins * a_small_intestine_sulf

    # Oral uptake of unconjugated BPA from small intestine into liver.
    r_oral <- k1 * a_small_intestine

    # -------------------------------------------------------------------
    # 7. Parent BPA mass balance (nmol/h). Doses enter stomach (oral) and
    #    the two skin-surface depots (dermal).
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
    # 8. BPA-glucuronide. Formed in liver (r_liver_gluc) and enterocytes
    #    (r_gi_gluc); split met1g directly to serum vs met2g into the
    #    gut EHR pool, which is reabsorbed (ehrrate) or deconjugated back to
    #    parent BPA in the liver (kenterobpag).
    # -------------------------------------------------------------------
    d/dt(a_small_intestine_gluc) <- r_gut_gluc - r_gi_gluc
    d/dt(gut_lumen_gluc) <- met2g * (r_liver_gluc + r_gi_gluc) -
      gut_lumen_gluc * (ehrrate + kenterobpag)
    d/dt(serum_gluc) <- met1g * (r_liver_gluc + r_gi_gluc) +
      gut_lumen_gluc * ehrrate - cl_gluc * c_gluc
    d/dt(urine_gluc) <- cl_gluc * c_gluc

    # -------------------------------------------------------------------
    # 9. BPA-sulfate. met1s = 1, so all hepatic sulfate goes directly to
    #    serum and the gut-lumen sulfate pool stays empty; carried for
    #    structural symmetry with the deposited code (gut sulfation off).
    # -------------------------------------------------------------------
    d/dt(a_small_intestine_sulf) <- -r_gi_sulf
    d/dt(gut_lumen_sulf) <- met2s * (r_liver_sulf + r_gi_sulf) -
      gut_lumen_sulf * (ehrrate + kenterobpas)
    d/dt(serum_sulf) <- met1s * (r_liver_sulf + r_gi_sulf) +
      gut_lumen_sulf * ehrrate - cl_sulf * c_sulf
    d/dt(urine_sulf) <- cl_sulf * c_sulf

    # -------------------------------------------------------------------
    # 10. Skin-surface depots. Dose enters with a modeled 24-h zero-order
    #     duration and the source-specific extent of absorption as
    #     bioavailability; drains first-order into the skin organ.
    # -------------------------------------------------------------------
    d/dt(depot_td_thermalpaper) <- -rate_tp
    d/dt(depot_td_personalcare) <- -rate_pcp
    f(depot_td_thermalpaper) <- eoa_tp
    f(depot_td_personalcare) <- eoa_pcp
    dur(depot_td_thermalpaper) <- period_dermal
    dur(depot_td_personalcare) <- period_dermal
    dur(stomach) <- period_oral

    # -------------------------------------------------------------------
    # 11. Observations (nmol/L = nM). Cc is the serum BPA concentration
    #     (the paper's reported 'serum' output). The conjugate and gonadal
    #     concentrations are returned alongside.
    # -------------------------------------------------------------------
    Cc <- c_art
    Cgluc <- c_gluc
    Csulf <- c_sulf
    Cgonad <- a_gonads / v_gonad
  })
}
