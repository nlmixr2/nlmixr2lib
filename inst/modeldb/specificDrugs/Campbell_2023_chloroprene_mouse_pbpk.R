Campbell_2023_chloroprene_mouse_pbpk <- function() {
  description <- paste(
    "Preclinical (female mouse). PBPK (whole-body, flow-limited, inhalation;",
    "recoded from the acslX / Berkeley-Madonna-style listing in the",
    "supplementary materials, itself translated from Yang et al. 2012).",
    "Parent beta-chloroprene (CP; 2-chloro-1,3-butadiene) is distributed over",
    "lung, liver, kidney, rapidly- and slowly-perfused tissue and fat as",
    "flow-limited compartments, with the lung in series (portal of entry):",
    "the gas-exchange arterial concentration is an algebraic steady-state lung",
    "equation coupling cardiac output to alveolar ventilation through the",
    "blood:air partition coefficient, so inhalation enters as a square-wave air",
    "concentration (ppm) rather than as a dose record, and lung venous blood",
    "supplies every other tissue. Saturable Michaelis-Menten hepatic and",
    "pulmonary oxidation of CP splits into the stable oxirane 1-CEO",
    "((1-chloroethenyl)oxirane) and the unstable 2-CEO. A linked 1-CEO",
    "submodel distributes 1-CEO over the same tissues (with a lung- and",
    "liver-epithelium diffusion barrier) and clears it by epoxide-hydrolase",
    "hydrolysis and, in the mouse only, a second microsomal oxidation to a",
    "diepoxide. The unstable 2-CEO and the diepoxide are lumped into a reactive",
    "product (RP) pool that reacts with glutathione (GSH, second order) and",
    "with other cellular constituents (first order); GSH follows zero-order",
    "synthesis and first-order turnover, so rising RP production depletes tissue",
    "GSH disproportionately with exposure. Deterministic typical-value model",
    "(no IIV, no residual error): the authors merged literature PBPK",
    "submodels and evaluated the result against in vitro kinetics and bioassay",
    "dose-response rather than fitting individuals. Outputs include venous CP",
    "blood concentration plus lung and liver GSH, RP and 1-CEO concentrations",
    "and cumulative metabolised / inhaled / exhaled amounts for the dose",
    "metrics (total CP metabolised per gram lung, average RP, average 1-CEO)."
  )
  reference <- paste(
    "Campbell JL Jr, Clewell HJ III, Van Landingham C, Gentry PR, Andersen ME.",
    "(2023). Using available in vitro metabolite identification and time course",
    "kinetics for beta-chloroprene and its metabolite, (1-chloroethenyl)oxirane,",
    "to include reactive oxidative metabolites and glutathione depletion in a",
    "PBPK model for beta-chloroprene. Front Pharmacol 14:1223808.",
    "doi:10.3389/fphar.2023.1223808. PMCID PMC10472072.",
    "Parent CP distribution/metabolism structure and physiology inherited from",
    "Clewell et al. 2019 and Yang et al. 2012 (reproduced in Tables 1, S2, S3, S6).",
    "Model code: Supplementary Material, 'Chloroprene model code for the epoxy submodel'.",
    sep = " "
  )
  vignette <- "Campbell_2023_chloroprene"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L",
    amount = "mg",
    weight = "kg"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "The model is parameterised for a female mouse with BW = 0.04 kg",
        "(NTP 1998 time- and survival-weighted average control body weight,",
        "Table S2). Carried as a covariate so the model can be re-scaled:",
        "all tissue volumes are linear in body weight, and cardiac output,",
        "alveolar ventilation, every metabolic Vmax and the GSH turnover rate",
        "are allometric (cardiac output and ventilation at exponent 0.75, Vmax",
        "at 0.75, GSH turnover KP at -0.3).",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  compartmentData <- list(
    # Parent chloroprene (CP), flow-limited tissue amounts (mg)
    a_lung = list(analyte = "chloroprene", units = "mg", specimen = "tissue", verified = TRUE),
    a_liver = list(analyte = "chloroprene", units = "mg", specimen = "tissue", verified = TRUE),
    a_kidney = list(analyte = "chloroprene", units = "mg", specimen = "tissue", verified = TRUE),
    a_slow = list(analyte = "chloroprene", units = "mg", specimen = "tissue", verified = TRUE),
    a_rapid = list(analyte = "chloroprene", units = "mg", specimen = "tissue", verified = TRUE),
    a_fat = list(analyte = "chloroprene", units = "mg", specimen = "tissue", verified = TRUE),
    a_inhaled = list(analyte = "chloroprene", units = "mg", specimen = "not applicable", verified = TRUE),
    a_exhaled = list(analyte = "chloroprene", units = "mg", specimen = "not applicable", verified = TRUE),
    a_met_liver = list(analyte = "chloroprene", units = "mg", specimen = "not applicable", verified = TRUE),
    a_met_lung = list(analyte = "chloroprene", units = "mg", specimen = "not applicable", verified = TRUE),
    # 1-CEO ((1-chloroethenyl)oxirane) amounts (mg); lung and liver each have an
    # epithelial sub-compartment separated from the vascular sub-compartment by
    # a diffusion barrier (permeation-area cross-product).
    a_ceo1_lung_epi = list(analyte = "(1-chloroethenyl)oxirane", units = "mg", specimen = "tissue", verified = TRUE),
    a_ceo1_lung = list(analyte = "(1-chloroethenyl)oxirane", units = "mg", specimen = "tissue", verified = TRUE),
    a_ceo1_liver_epi = list(analyte = "(1-chloroethenyl)oxirane", units = "mg", specimen = "tissue", verified = TRUE),
    a_ceo1_liver = list(analyte = "(1-chloroethenyl)oxirane", units = "mg", specimen = "tissue", verified = TRUE),
    a_ceo1_slow = list(analyte = "(1-chloroethenyl)oxirane", units = "mg", specimen = "tissue", verified = TRUE),
    a_ceo1_rapid = list(analyte = "(1-chloroethenyl)oxirane", units = "mg", specimen = "tissue", verified = TRUE),
    a_ceo1_fat = list(analyte = "(1-chloroethenyl)oxirane", units = "mg", specimen = "tissue", verified = TRUE),
    a_ceo1_exhaled = list(
      analyte = "(1-chloroethenyl)oxirane",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    # Reactive product (RP) pool, tissue concentration (uM)
    c_rp_liver = list(
      analyte = "reactive products (CP oxidative metabolites)",
      units = "uM",
      specimen = "tissue",
      verified = TRUE
    ),
    c_rp_lung = list(
      analyte = "reactive products (CP oxidative metabolites)",
      units = "uM",
      specimen = "tissue",
      verified = TRUE
    ),
    # Glutathione, tissue concentration (uM)
    c_gsh_liver = list(analyte = "glutathione", units = "uM", specimen = "tissue", verified = TRUE),
    c_gsh_lung = list(analyte = "glutathione", units = "uM", specimen = "tissue", verified = TRUE),
    # AUC accumulators for the reported average-concentration dose metrics
    auc_rp_liver = list(
      analyte = "reactive products (CP oxidative metabolites)",
      units = "uM*h",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_rp_lung = list(
      analyte = "reactive products (CP oxidative metabolites)",
      units = "uM*h",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_ceo1_liver = list(
      analyte = "(1-chloroethenyl)oxirane",
      units = "uM*h",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_ceo1_lung = list(
      analyte = "(1-chloroethenyl)oxirane",
      units = "uM*h",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_gsh_liver = list(analyte = "glutathione", units = "uM*h", specimen = "not applicable", verified = TRUE),
    auc_gsh_lung = list(analyte = "glutathione", units = "uM*h", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "mouse (female B6C3F1)",
    n_subjects = 0L,
    n_studies = 0L,
    age_range = "adult (single reference female mouse; no age term in the model)",
    weight_range = "0.04 kg (NTP 1998 time- and survival-weighted average control body weight, Table S2)",
    sex_female_pct = 100,
    race_ethnicity = NA_character_,
    disease_state = "healthy reference animal (inhalation toxicology model)",
    dose_range = paste(
      "Inhalation of beta-chloroprene vapour. The NTP (1998) bioassay",
      "concentrations 0 (control), 12.8, 32 and 80 ppm, and single-6 h and",
      "14-day (6 h/day, 5 days/week) exposures across 0.3-several hundred ppm",
      "for the dose-metric and GSH-depletion analyses (Figures 3, 4; Table S7).",
      sep = " "
    ),
    regions = "United States (NTP bioassay; model developed by Ramboll / Andersen ToxConsulting)",
    notes = paste(
      "This is a forward toxicokinetic PBPK simulation, not a fit, so",
      "n_subjects = 0. Physiological parameters are from Brown et al. 1997 and",
      "Marino et al. 2006 (Table S2); parent-CP partition coefficients from",
      "Himmelstein et al. 2004 (Table S3); parent-CP oxidation kinetics from",
      "Clewell et al. 2019 (reproduced in Table S6); 1-CEO metabolite kinetics,",
      "partition coefficients and reactive-product / GSH constants from Table 1",
      "(Himmelstein et al. 2004, Campbell et al. 2015, Clewell et al. 2001,",
      "D'Souza et al. 1988, Csanady et al. 2003). The second microsomal",
      "oxidation of 1-CEO to a diepoxide is a mouse-only pathway. The companion",
      "female-rat parameterisation is not shipped because its parent-CP",
      "oxidation Vmax/Km are not reported in this paper (they live in the",
      "upstream Clewell et al. 2019 model).",
      sep = " "
    )
  )

  # Whole-body PBPK tissue, metabolite, reactive-product, glutathione and AUC
  # states use descriptive per-paper names rather than the canonical PK
  # compartment set; declare them so checkModelConventions() does not flag them.
  paper_specific_compartment_pattern <- "^(a_|c_rp_|c_gsh_|auc_)"

  ini({
    # =====================================================================
    # Every value is a recoded literature constant, so every parameter is
    # fixed(); nothing was estimated from data in this paper.
    # =====================================================================

    # --- Physiology, female mouse (Table S2) ---------------------------------
    qp_per_kg <- fixed(29.1); label("Alveolar ventilation (L/h/kg^0.75)") # Table S2, QPC, Brown et al. 1997
    qcc <- fixed(20.1); label("Cardiac output (L/h/kg^0.75)") # Table S2, QCC, Marino et al. 2006
    e_wt_flow <- fixed(0.75); label("Allometric exponent on flows and Vmax (unitless)") # supplement code, pow(BW,0.75)
    e_wt_kp <- fixed(-0.3); label("Allometric exponent on GSH turnover KP (unitless)") # Table 1, KPC units '/hr*kg BW^-0.3'

    qlc <- fixed(0.161); label("Flow to liver, fraction of cardiac output") # Table S2, QLC
    qfc <- fixed(0.07); label("Flow to fat, fraction of cardiac output") # Table S2, QFC
    qsc <- fixed(0.159); label("Flow to slowly-perfused tissue, fraction of cardiac output") # Table S2, QSC
    qkc <- fixed(0.09); label("Flow to kidney, fraction of cardiac output") # Table S2, QKC
    qrc <- fixed(0.52); label("Flow to rapidly-perfused tissue, fraction of cardiac output") # supplement code, QRC

    vlc <- fixed(0.0549); label("Liver volume, fraction of body weight") # Table S2, VLC
    vluc <- fixed(0.0073); label("Lung volume, fraction of body weight") # Table S2, VLUC
    vfc <- fixed(0.1); label("Fat volume, fraction of body weight") # Table S2, VFC
    vrc <- fixed(0.08098); label("Rapidly-perfused volume, fraction of body weight") # Table S2, VRC
    vsc <- fixed(0.384); label("Slowly-perfused volume, fraction of body weight") # Table S2, VSC
    vkc <- fixed(0.0167); label("Kidney volume, fraction of body weight") # Table S2, VKC
    vtot_b <- fixed(0.64398); label("Perfused-tissue fraction used to normalise tissue volumes") # supplement code, VTOT_B

    # --- Partition coefficients, parent CP, mouse (Table S3) -----------------
    pb <- fixed(7.83); label("CP blood:air partition coefficient (unitless)") # Table S3, Himmelstein et al. 2004
    plu <- fixed(2.38); label("CP lung:blood partition coefficient (unitless)") # Table S3
    pl <- fixed(1.26); label("CP liver:blood partition coefficient (unitless)") # Table S3
    pf <- fixed(17.35); label("CP fat:blood partition coefficient (unitless)") # Table S3
    ps <- fixed(0.59); label("CP slowly-perfused (muscle):blood partition coefficient (unitless)") # Table S3
    pr <- fixed(1.76); label("CP rapidly-perfused (kidney):blood partition coefficient (unitless)") # Table S3
    pk <- fixed(1.76); label("CP kidney:blood partition coefficient (unitless)") # Table S3, 'used for rapidly perfused'

    # --- Partition coefficients, 1-CEO (Table 1; IndusChemFate, logKow 1.22) --
    pb1 <- fixed(5.74); label("1-CEO blood:air partition coefficient (unitless)") # Table 1, PB1
    plu1 <- fixed(0.69); label("1-CEO lung:blood partition coefficient (unitless)") # Table 1, PLU1
    pl1 <- fixed(1.18); label("1-CEO liver:blood partition coefficient (unitless)") # Table 1, PL1
    pf1 <- fixed(5.15); label("1-CEO fat:blood partition coefficient (unitless)") # Table 1, PF1
    ps1 <- fixed(0.69); label("1-CEO slowly-perfused:blood partition coefficient (unitless)") # Table 1, PS1
    pr1 <- fixed(1.18); label("1-CEO rapidly-perfused:blood partition coefficient (unitless)") # Table 1, PR1

    # --- Molecular weights and ppm conversion -------------------------------
    mw <- fixed(88.5); label("Chloroprene molecular weight (g/mol)") # supplement code, MW
    mwceo <- fixed(104.5); label("1-CEO molecular weight (g/mol)") # supplement code, MWCEO
    molar_vol <- fixed(24450); label("Ideal-gas molar volume for ppm<->mg/L (mL/mol)") # supplement code, CIX=CONC*MW/24450

    # --- Parent CP oxidation, mouse (Table S6 central value; Clewell 2019) ----
    vmaxc <- fixed(7.99); label("CP oxidation Vmax, liver (mg/h/kg^0.75)") # Table S6, VMAXC
    km <- fixed(0.040); label("CP oxidation Km, liver (mg/L)") # Table S6, KM
    vmaxclu <- fixed(0.120); label("CP oxidation Vmax, lung (mg/h/kg^0.75)") # Table S6, VMAXCLU
    kmlu <- fixed(0.210); label("CP oxidation Km, lung (mg/L)") # Table S6, KMLU

    # --- Split of CP oxidation to 1-CEO vs 2-CEO (Table 1, mouse) ------------
    alphal <- fixed(0.02); label("Fraction of liver CP oxidation to 1-CEO (remainder to 2-CEO)") # Table 1, ALPHAL
    alphalu <- fixed(0.03); label("Fraction of lung CP oxidation to 1-CEO (remainder to 2-CEO)") # Table 1, ALPHALU
    # BETA is the fraction of formed 1-CEO that escapes intracellular first-pass
    # hydration to the diol and is available for distribution / metabolism.
    # Table 1 reports 0.33; equivalently 1 - 0.67, the 67% of butadiene-derived
    # epoxybutene further metabolised by co-localised EH (Campbell et al. 2015).
    beta <- fixed(0.33); label("Fraction of 1-CEO production available, not first-pass diol") # Table 1, BETA

    # --- 1-CEO metabolism, mouse (Table 1) ----------------------------------
    vmaxc1 <- fixed(10.65); label("1-CEO hydrolysis Vmax, liver (mg/h/kg^0.75)") # Table 1, VMAXC1
    km1 <- fixed(1.9); label("1-CEO hydrolysis Km, liver (mg/L)") # Table 1, KM1
    vmaxclu1 <- fixed(0.64); label("1-CEO hydrolysis Vmax, lung (mg/h/kg^0.75)") # Table 1, VMAXCLU1
    kmlu1 <- fixed(4.6); label("1-CEO hydrolysis Km, lung (mg/L)") # Table 1, KMLU1
    vmaxc1o <- fixed(2.25); label("1-CEO oxidation Vmax, liver, mouse only (mg/h/kg^0.75)") # Table 1, VMAXC1O
    km1o <- fixed(1.5); label("1-CEO oxidation Km, liver (mg/L)") # Table 1, KM1O
    lloxact <- fixed(0.42); label("Lung:liver ratio for 1-CEO oxidation (unitless)") # Table 1, LLOXACT
    pa1 <- fixed(1.0); label("1-CEO tissue<->blood permeation coefficient (fraction of blood flow)") # supplement code, PA1

    # --- Reactive product (RP) pool and GSH, mouse (Table 1) ----------------
    kgshlc <- fixed(0.13); label("2nd-order rate of RP reaction with GSH, liver (L/umol/h)") # Table 1, KGSHLC
    kgshluc <- fixed(0.13); label("2nd-order rate of RP reaction with GSH, lung (L/umol/h)") # Table 1, KGSHLUC
    # KFEE is a first-order disposal of RP by 'everything else' (text); Table 1
    # labels its units L/umol/h but that is a copy of the KGSH unit -- the rate
    # equation (Figure 2 Eq 2) and the text both use it first order, so the
    # value 35 is applied as 1/h. See vignette Errata.
    kfeec <- fixed(35); label("First-order RP disposal by non-GSH constituents, liver (1/h)") # Table 1, KFEEC
    lleeact <- fixed(0.14); label("Lung:liver ratio for RP disposal by non-GSH (unitless)") # Table 1, LLEEACT
    kpc <- fixed(0.06); label("GSH turnover rate constant (/h/kg^-0.3)") # Table 1, KPC
    gso <- fixed(7000); label("Initial / steady-state liver GSH concentration (uM)") # Table 1, GSO
    gsolu <- fixed(1500); label("Initial / steady-state lung GSH concentration (uM)") # Table 1, GSOLU

    # --- Exposure schedule (overridden per scenario via rxSolve(params=)) ----
    conc_air_ppm <- fixed(0); label("Inhaled chloroprene air concentration during an exposure window (ppm)") # supplement code, CONC
    exp_hours_per_day <- fixed(6); label("Duration of each daily exposure window (h)") # Methods, 6 h/day
    exp_days_per_week <- fixed(5); label("Exposure days per 7-day week") # Methods, 5 days/week

    # No IIV and no residual error: deterministic recoded literature PBPK model
    # evaluated against in vitro kinetics and bioassay dose-response, not fitted.
  })

  model({
    # -------------------------------------------------------------------
    # 1. Scaled volumes (L), blood flows (L/h) and metabolic capacities.
    #    Supplement code block 'Scaled parameters'. qc equals cardiac
    #    output because the fractional flows are renormalised by their own
    #    sum (qtot_frac) before being multiplied back up.
    # -------------------------------------------------------------------
    qp <- qp_per_kg * WT^e_wt_flow
    qci <- qcc * WT^e_wt_flow
    qtot_frac <- qrc + qlc + qkc + qfc + qsc
    ql <- qlc * qci / qtot_frac
    qf <- qfc * qci / qtot_frac
    qs <- qsc * qci / qtot_frac
    qk <- qkc * qci / qtot_frac
    qr <- qrc * qci / qtot_frac
    qc <- ql + qf + qs + qk + qr

    vtot <- vlc + vluc + vfc + vsc + vrc + vkc
    vl <- vlc * WT * (vtot_b / vtot)
    vlu <- vluc * WT * (vtot_b / vtot)
    vf <- vfc * WT * (vtot_b / vtot)
    vs <- vsc * WT * (vtot_b / vtot)
    vr <- vrc * WT * (vtot_b / vtot)
    vk <- vkc * WT * (vtot_b / vtot)

    vmax <- vmaxc * WT^e_wt_flow
    vmaxlu <- vmaxclu * WT^e_wt_flow
    vmax1 <- vmaxc1 * WT^e_wt_flow
    vmaxlu1 <- vmaxclu1 * WT^e_wt_flow
    vmax1o <- vmaxc1o * WT^e_wt_flow
    vmaxlu1o <- vmax1o * lloxact
    kfee <- kfeec
    kfeelu <- kfee * lleeact
    kp <- kpc * WT^e_wt_kp

    # -------------------------------------------------------------------
    # 2. Inhalation exposure as a repeated square wave: on for
    #    exp_hours_per_day each of the first exp_days_per_week days of every
    #    7-day (168 h) week. MOD(TIME,P) is written time - P*floor(time/P)
    #    and each IF/THEN switch as a comparison that rxode2 evaluates to 1/0.
    # -------------------------------------------------------------------
    t_week <- time - 168 * floor(time / 168)
    t_day <- time - 24 * floor(time / 24)
    air_on <- (t_week < 24 * exp_days_per_week) * (t_day < exp_hours_per_day)
    cix <- conc_air_ppm * mw / molar_vol # inhaled air concentration (mg/L)
    ci <- cix * air_on

    # -------------------------------------------------------------------
    # 3. Parent CP tissue venous-effluent concentrations (mg/L) and the
    #    algebraic gas-exchange (arterial) and mixed-venous concentrations.
    #    The lung is in series: cpu supplies the lung, lung venous (cvlu)
    #    supplies every other tissue, and cv returns to the lung.
    # -------------------------------------------------------------------
    cvlu <- a_lung / (vlu * plu)
    cvl <- a_liver / (vl * pl)
    cvk <- a_kidney / (vk * pk)
    cvs <- a_slow / (vs * ps)
    cvr <- a_rapid / (vr * pr)
    cvf <- a_fat / (vf * pf)
    cpu <- (qp * ci + (qf * cvf + ql * cvl + qs * cvs + qr * cvr + qk * cvk)) / (qp / pb + qc)
    cx <- cpu / pb
    cv <- (qf * cvf + ql * cvl + qs * cvs + qr * cvr + qk * cvk) / qc

    # CP oxidation (Michaelis-Menten) in liver and lung and its split to the
    # two oxiranes (mg/h of CP mass).
    ram <- vmax * cvl / (km + cvl)
    ramlu <- vmaxlu * cvlu / (kmlu + cvlu)
    ramcp_1ceo <- ram * alphal
    ramcp_2ceo <- ram * (1 - alphal)
    ramlucp_1ceo <- ramlu * alphalu
    ramlucp_2ceo <- ramlu * (1 - alphalu)

    # -------------------------------------------------------------------
    # 4. Parent CP mass balance (amounts, mg).
    # -------------------------------------------------------------------
    d/dt(a_lung) <- qc * (cpu - cvlu) - ramlu
    d/dt(a_liver) <- ql * (cvlu - cvl) - ram
    d/dt(a_kidney) <- qk * (cvlu - cvk)
    d/dt(a_slow) <- qs * (cvlu - cvs)
    d/dt(a_rapid) <- qr * (cvlu - cvr)
    d/dt(a_fat) <- qf * (cvlu - cvf)
    d/dt(a_inhaled) <- qp * ci
    d/dt(a_exhaled) <- qp * cx
    d/dt(a_met_liver) <- ram
    d/dt(a_met_lung) <- ramlu

    # -------------------------------------------------------------------
    # 5. 1-CEO submodel. 1-CEO amounts (mg). Lung and liver each split into
    #    an epithelial sub-compartment (where production and metabolism
    #    occur) and a vascular sub-compartment, coupled by a permeation-area
    #    cross-product (pa1 * flow). Kidney 1-CEO is lumped with the rapidly-
    #    perfused group (flow qr + qk), as in the supplement code.
    # -------------------------------------------------------------------
    clue1 <- a_ceo1_lung_epi / vlu
    cvlu1 <- a_ceo1_lung / (vlu * plu1)
    cle1 <- a_ceo1_liver_epi / vl
    cvl1 <- a_ceo1_liver / (vl * pl1)
    cvs1 <- a_ceo1_slow / (vs * ps1)
    cvr1 <- a_ceo1_rapid / (vr * pr1)
    cvf1 <- a_ceo1_fat / (vf * pf1)
    cv1 <- (qf * cvf1 + ql * cvl1 + qs * cvs1 + (qr + qk) * cvr1) / qc
    cpu1 <- (qc * cv1) / (qp / pb1 + qc)
    cx1 <- cpu1 / pb1
    # Total 1-CEO tissue concentration (both sub-compartments)
    cl1 <- (a_ceo1_liver_epi + a_ceo1_liver) / vl
    clu1 <- (a_ceo1_lung_epi + a_ceo1_lung) / vlu

    # 1-CEO hydrolysis (both species) and oxidation to diepoxide (mouse only)
    ram1 <- vmax1 * cvl1 / (km1 + cvl1)
    ramlu1 <- vmaxlu1 * cvlu1 / (kmlu1 + cvlu1)
    ram1o <- vmax1o * cvl1 / (km1o + cvl1)
    ramlu1o <- vmaxlu1o * cvlu1 / (km1o + cvlu1)

    # Conversion of CP-mass 1-CEO production to 1-CEO mass: * mwceo / mw.
    d/dt(a_ceo1_lung_epi) <- pa1 * qc * (cvlu1 - clue1 / plu1) + ramlucp_1ceo * beta * mwceo / mw - ramlu1 - ramlu1o
    d/dt(a_ceo1_lung) <- qc * (cpu1 - cvlu1) + pa1 * qc * (clue1 / plu1 - cvlu1)
    d/dt(a_ceo1_liver_epi) <- pa1 * ql * (cvl1 - cle1 / pl1) + ramcp_1ceo * beta * mwceo / mw - ram1 - ram1o
    d/dt(a_ceo1_liver) <- ql * (cvlu1 - cvl1) + pa1 * ql * (cle1 / pl1 - cvl1)
    d/dt(a_ceo1_slow) <- qs * (cvlu1 - cvs1)
    d/dt(a_ceo1_rapid) <- (qr + qk) * (cvlu1 - cvr1)
    d/dt(a_ceo1_fat) <- qf * (cvlu1 - cvf1)
    d/dt(a_ceo1_exhaled) <- qp * cx1

    # -------------------------------------------------------------------
    # 6. Reactive product (RP) pool and glutathione (GSH), tissue
    #    concentrations (uM). These are Figure 2 Eqs 2 and 3 written in
    #    concentration form (the supplement code integrates the equivalent
    #    amounts). RP is produced from 2-CEO (the (1-alpha) oxidation branch)
    #    and, in the mouse, from the 1-CEO diepoxide (second oxidation); it is
    #    consumed by second-order reaction with GSH and first-order reaction
    #    with other constituents. GSH has zero-order synthesis (kp * baseline)
    #    and first-order turnover, so RP production depletes GSH. Production
    #    terms convert metabolised mass to uM/h: (mg/h / (g/mol)) * 1000 / L.
    # -------------------------------------------------------------------
    rp_prod_liver <- ((ramcp_2ceo / mw) * 1000 + (ram1o / mwceo) * 1000) / vl
    rp_prod_lung <- ((ramlucp_2ceo / mw) * 1000 + (ramlu1o / mwceo) * 1000) / vlu

    d/dt(c_rp_liver) <- rp_prod_liver - kgshlc * c_gsh_liver * c_rp_liver - kfee * c_rp_liver
    d/dt(c_rp_lung) <- rp_prod_lung - kgshluc * c_gsh_lung * c_rp_lung - kfeelu * c_rp_lung
    d/dt(c_gsh_liver) <- kp * gso - kp * c_gsh_liver - kgshlc * c_rp_liver * c_gsh_liver
    d/dt(c_gsh_lung) <- kp * gsolu - kp * c_gsh_lung - kgshluc * c_rp_lung * c_gsh_lung
    c_gsh_liver(0) <- gso
    c_gsh_lung(0) <- gsolu

    # -------------------------------------------------------------------
    # 7. AUC accumulators for the reported average-concentration dose metrics
    #    (average over the simulated window = AUC / elapsed time).
    # -------------------------------------------------------------------
    d/dt(auc_rp_liver) <- c_rp_liver
    d/dt(auc_rp_lung) <- c_rp_lung
    d/dt(auc_ceo1_liver) <- cl1 / mwceo * 1000
    d/dt(auc_ceo1_lung) <- clu1 / mwceo * 1000
    d/dt(auc_gsh_liver) <- c_gsh_liver
    d/dt(auc_gsh_lung) <- c_gsh_lung

    # -------------------------------------------------------------------
    # 8. Outputs. Cc is the mixed venous chloroprene blood concentration
    #    (the systemic blood level); the mechanistic lung/liver outputs are
    #    the quantities the paper tracks for its dose metrics.
    # -------------------------------------------------------------------
    Cc <- cv
    Cexhaled <- cx
    Cgsh_lung <- c_gsh_lung
    Cgsh_liver <- c_gsh_liver
    Crp_lung <- c_rp_lung
    Crp_liver <- c_rp_liver
    Cceo1_lung <- clu1 / mwceo * 1000
    Cceo1_liver <- cl1 / mwceo * 1000
  })
}
