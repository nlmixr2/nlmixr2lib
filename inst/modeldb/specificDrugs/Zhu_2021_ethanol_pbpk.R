Zhu_2021_ethanol_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, hand-written MATLAB; dynamic flux-balance coupling",
    "to the Harvey genome-scale model reduced to its continuous limit).",
    "Ethanol and acetaldehyde disposition after an oral drink in adult men",
    "(Zhu 2021 PLoS Comput Biol). Fourteen flow-limited tissues (adipose,",
    "arterial blood, brain, small and large intestine, heart, kidney,",
    "liver, lung, muscle, pancreas, skin, spleen, stomach) per analyte plus",
    "stomach and small-intestine lumens for ethanol. Organ masses, blood",
    "flows and cardiac output from age / height / weight / body-fat",
    "correlations; partition coefficients from logP, fraction unbound and",
    "tissue lipid/water composition. Stomach and small-intestine",
    "absorption and stomach-to-SI transit are quadratic functions of the",
    "drink ethanol fraction (S2 Text, fitted to Mitchell 2014). Hepatic",
    "ADH and ALDH2 and gastric ADH follow Michaelis-Menten kinetics",
    "(Umulis 2005 constants) with an age factor on all three. ALDH2",
    "isoform activity (Table 4) and a steady blood disulfiram level (S3",
    "Text) scale hepatic ALDH2. Urine, sweat and breath excretion are the",
    "Harvey flux-balance bounds expressed as fixed fractions of the",
    "hepatic ADH rate (set them to 0 to recover the base PBPK used for",
    "Figs 6-8). Deterministic: no IIV, no residual error. Encoded from the",
    "authors' deposited code, which reproduces Figs 4-8; the printed kSI",
    "correlation (S2.3) has a sign typo. Male only. See the vignette."
  )
  reference <- paste(
    "Zhu L, Pei W, Thiele I, Mahadevan R (2021). Integration of a",
    "physiologically-based pharmacokinetic model with a whole-body,",
    "organ-resolved genome-scale model for characterization of ethanol and",
    "acetaldehyde metabolism. PLoS Comput Biol 17(8):e1009110.",
    "doi:10.1371/journal.pcbi.1009110.",
    "Model code: https://github.com/LMSE/HH-PBPK-Ethanol (commit b14305a).",
    "Michaelis-Menten constants from Umulis DM, Gurmen NM, Singh P,",
    "et al. (2005) Alcohol 35(1):3-12, doi:10.1016/j.alcohol.2004.11.004.",
    sep = " "
  )
  vignette <- "Zhu_2021_ethanol"
  units <- list(time = "min", dosing = "mmol", concentration = "mmol/L")

  # Organ tissue states and the acetaldehyde (`_acald`) mirror of each have
  # no canonical name in the register (a_pancreas, a_stomach, the lumen
  # of the small intestine, and the acetaldehyde metabolite suffix). The
  # whitelist below names them explicitly rather than registering new
  # canonicals for a single paper.
  paper_specific_compartments <- c(
    "a_pancreas",
    "a_stomach",
    "small_intestine_lumen"
  )
  paper_specific_compartment_pattern <- "_acald$"

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the organ-mass and blood-flow correlations (S1 Text Eq",
        "S1.16 for cardiac output; organVolM.m / organFlow.m in the",
        "deposited code, citing Stader 2019) and the ADH/ALDH age factor",
        "(maxrates.m): 1 at <= 25 years, 0.5 above 55 years, and",
        "-0.00102*AGE^2 + 0.067*AGE - 0.0663 between. Paper simulations use",
        "25.6 years (Umulis / Jones scenario, Figs 4-8) and 37.8 years",
        "(Mitchell scenario, Fig 3)."
      ),
      source_name = "age"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters DuBois body surface area (0.007184 * WT^0.425 * HT^0.725),",
        "the brain, intestine and lung mass correlations, and the heart and",
        "kidney blood-flow fractions. Paper simulations use 180 cm (Figs",
        "4-8) and 177.1 cm (Fig 3)."
      ),
      source_name = "height"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters body surface area, adipose mass (WT * BODYFAT_PCT / 100)",
        "and the residual mass (WT minus the sum of the correlated organ",
        "masses), which organVolM.m redistributes 10 / 10 / 80 percent to",
        "small intestine, large intestine and skin. Paper simulations use",
        "74.5 kg (Figs 4-8) and 82.66 kg (Fig 3). The dose is supplied in",
        "mmol, so a g/kg drink must be converted with WT in the event table."
      ),
      source_name = "bodyMass"
    ),
    BODYFAT_PCT = list(
      description = "Percent total body fat",
      units = "% (percent, 0-100)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Adipose mass = 0.01 * BODYFAT_PCT * WT (organVolM.m). All paper",
        "simulations use 20 percent."
      ),
      source_name = "fat"
    )
  )

  compartmentData <- list(
    stomach = list(analyte = "ethanol", units = "mmol", specimen = "administration site", verified = TRUE),
    small_intestine_lumen = list(
      analyte = "ethanol",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    a_fat = list(analyte = "ethanol", units = "mmol", specimen = "tissue", verified = TRUE),
    a_arterial = list(analyte = "ethanol", units = "mmol", specimen = "whole blood", verified = TRUE),
    a_brain = list(analyte = "ethanol", units = "mmol", specimen = "tissue", verified = TRUE),
    a_small_intestine = list(analyte = "ethanol", units = "mmol", specimen = "tissue", verified = TRUE),
    a_large_intestine = list(analyte = "ethanol", units = "mmol", specimen = "tissue", verified = TRUE),
    a_heart = list(analyte = "ethanol", units = "mmol", specimen = "tissue", verified = TRUE),
    a_kidney = list(analyte = "ethanol", units = "mmol", specimen = "tissue", verified = TRUE),
    a_liver = list(analyte = "ethanol", units = "mmol", specimen = "tissue", verified = TRUE),
    a_lung = list(analyte = "ethanol", units = "mmol", specimen = "tissue", verified = TRUE),
    a_muscle = list(analyte = "ethanol", units = "mmol", specimen = "tissue", verified = TRUE),
    a_pancreas = list(analyte = "ethanol", units = "mmol", specimen = "tissue", verified = TRUE),
    a_skin = list(analyte = "ethanol", units = "mmol", specimen = "tissue", verified = TRUE),
    a_spleen = list(analyte = "ethanol", units = "mmol", specimen = "tissue", verified = TRUE),
    a_stomach = list(analyte = "ethanol", units = "mmol", specimen = "tissue", verified = TRUE),
    a_fat_acald = list(analyte = "acetaldehyde", units = "mmol", specimen = "tissue", verified = TRUE),
    a_arterial_acald = list(analyte = "acetaldehyde", units = "mmol", specimen = "whole blood", verified = TRUE),
    a_brain_acald = list(analyte = "acetaldehyde", units = "mmol", specimen = "tissue", verified = TRUE),
    a_small_intestine_acald = list(analyte = "acetaldehyde", units = "mmol", specimen = "tissue", verified = TRUE),
    a_large_intestine_acald = list(analyte = "acetaldehyde", units = "mmol", specimen = "tissue", verified = TRUE),
    a_heart_acald = list(analyte = "acetaldehyde", units = "mmol", specimen = "tissue", verified = TRUE),
    a_kidney_acald = list(analyte = "acetaldehyde", units = "mmol", specimen = "tissue", verified = TRUE),
    a_liver_acald = list(analyte = "acetaldehyde", units = "mmol", specimen = "tissue", verified = TRUE),
    a_lung_acald = list(analyte = "acetaldehyde", units = "mmol", specimen = "tissue", verified = TRUE),
    a_muscle_acald = list(analyte = "acetaldehyde", units = "mmol", specimen = "tissue", verified = TRUE),
    a_pancreas_acald = list(analyte = "acetaldehyde", units = "mmol", specimen = "tissue", verified = TRUE),
    a_skin_acald = list(analyte = "acetaldehyde", units = "mmol", specimen = "tissue", verified = TRUE),
    a_spleen_acald = list(analyte = "acetaldehyde", units = "mmol", specimen = "tissue", verified = TRUE),
    a_stomach_acald = list(analyte = "acetaldehyde", units = "mmol", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 25L,
    n_studies = 2L,
    age_range = paste(
      "group means only: 37.8 years (Mitchell 2014, n = 15, absorption",
      "fit) and 25.6 years (Jones 1988, n = 10, acetaldehyde comparison)"
    ),
    weight_range = "group means only: 82.66 kg (Mitchell) and 74.5 kg (Jones / Umulis)",
    sex_female_pct = 0,
    disease_state = "Healthy adult men (literature data; no new subjects were studied).",
    dose_range = paste(
      "Single oral ethanol drinks of 0.25 and 0.5 g/kg at 5.1, 12.5, 20",
      "and 40 percent (w/w) ethanol; a repeated-drink scenario adds 14 g",
      "standard drinks at 60, 120 and 180 min."
    ),
    regions = "Not reported (published literature data).",
    notes = paste(
      "No subject-level data were fitted. The absorption correlations",
      "(S2 Text) were fitted to the mean blood-ethanol profiles of",
      "Mitchell et al. 2014 (15 men; mean age 37.8 y, 82.66 kg, 177.1 cm,",
      "20 percent body fat; 0.5 g/kg as beer, wine or spirits; Fig 3).",
      "Hepatic Michaelis-Menten constants come from Umulis 2005, whose",
      "model was fitted to Jones 1988 (10 men; the paper simulates a",
      "25.6 y, 74.5 kg, 180 cm, 20 percent body-fat man drinking",
      "0.25 g/kg; Fig 4). Every other paper scenario (Figs 5-8) uses that",
      "same simulated man. The authors' code does not support women",
      "(main.m: 'female (not yet supported)')."
    )
  )

  ini({
    # ---- Physicochemical properties (Table 1; organPartition.m) -------------
    logp_etoh <- fixed(-0.31); label("Ethanol lipophilicity log10 P (unitless)") # Table 1, 'Lipophilicity' ethanol
    fu_etoh <- fixed(0.99); label("Ethanol fraction unbound in blood (unitless)") # Table 1, 'Fraction Unbound' ethanol
    logp_acald <- fixed(-0.34); label("Acetaldehyde lipophilicity log10 P (unitless)") # Table 1, 'Lipophilicity' acetaldehyde
    fu_acald <- fixed(0.99); label("Acetaldehyde fraction unbound in blood (unitless)") # Table 1, 'Fraction Unbound' acetaldehyde

    # ---- Oral input (main.m) -----------------------------------------------
    # main.m: c0(15) = 0.8*gkg*bodyMass/molarMassEtOH*1000/volStomLumen,
    # '80% bioavailability'; volStomLumen = 1.10 'volume of gastric fluid in
    # L'. The small-intestine lumen shares this volume implicitly (the code
    # passes lumen concentration 1:1 from stomach to SI lumen).
    fstomach <- fixed(0.8); label("Fraction of the ingested ethanol dose entering the stomach lumen (unitless)") # main.m c0(15), '80% bioavailability'
    v_stomach_lumen <- fixed(1.10); label("Gastric (and small-intestine) luminal fluid volume (L)") # main.m volStomLumen = 1.10

    # ---- Absorption / transit vs drink ethanol fraction (S2 Text) ----------
    # k = c2*drink_frac^2 + c1*drink_frac + c0, drink_frac as a fraction
    # (0.2 = 20 percent w/w). Fitted by the authors to Mitchell 2014 (Fig 3).
    drink_frac <- fixed(0.2); label("Drink ethanol fraction, w/w (unitless; 0.2 = 20 percent)") # Figs 4-8 scenario; lumen decay of the deposited Fig 4 trajectory back-solves 0.200
    kstom_c2 <- fixed(0.7135); label("Stomach absorption rate, quadratic coefficient (1/min)") # S2 Text Eq S2.1
    kstom_c1 <- fixed(-0.0985); label("Stomach absorption rate, linear coefficient (1/min)") # S2 Text Eq S2.1
    kstom_c0 <- fixed(0.0112); label("Stomach absorption rate, intercept (1/min)") # S2 Text Eq S2.1
    kstomsi_c2 <- fixed(1.953); label("Stomach-to-small-intestine lumen transit rate, quadratic coefficient (1/min)") # S2 Text Eq S2.2
    kstomsi_c1 <- fixed(-0.168); label("Stomach-to-small-intestine lumen transit rate, linear coefficient (1/min)") # S2 Text Eq S2.2
    kstomsi_c0 <- fixed(0.0255); label("Stomach-to-small-intestine lumen transit rate, intercept (1/min)") # S2 Text Eq S2.2
    ksi_c2 <- fixed(-0.006); label("Small-intestine absorption rate, quadratic coefficient (1/min)") # S2 Text Eq S2.3
    # S2.3 prints '- 0.0686 * Drink%'; main.m has +0.0686. The code sign
    # reproduces the deposited small-intestine-lumen trajectory exactly.
    ksi_c1 <- fixed(0.0686); label("Small-intestine absorption rate, linear coefficient (1/min)") # main.m kSI (S2 Text Eq S2.3 prints the opposite sign)
    ksi_c0 <- fixed(0.0615); label("Small-intestine absorption rate, intercept (1/min)") # S2 Text Eq S2.3

    # ---- Michaelis-Menten kinetics (Table 2; maxrates.m) -------------------
    vmax_adh_liver <- fixed(2.2); label("Hepatic ADH Vmax (mmol/L/min)") # maxrates.m row 8 'liver, Umulis et al'; Table 2 source [11]
    km_adh_liver <- fixed(1); label("Hepatic ADH Km (mmol/L)") # maxrates.m row 8
    vmax_adh_stomach <- fixed(0.68); label("Gastric ADH Vmax (mmol/L/min)") # maxrates.m row 14 'stomach, Toroghi et al'
    km_adh_stomach <- fixed(41); label("Gastric ADH Km (mmol/L)") # maxrates.m row 14
    vmax_aldh_liver <- fixed(2.7); label("Hepatic ALDH2 Vmax (mmol/L/min)") # maxrates.m row 25 'liver, Umulis et al'
    km_aldh_liver <- fixed(1.6); label("Hepatic ALDH2 Km (mmol/L)") # maxrates.m row 25
    # Only 1/360 of the hepatic ADH ethanol flux enters the liver
    # acetaldehyde pool (ODE.m dC(25): '+ aLiver*r(8)/360').
    f_acald <- fixed(1 / 360); label("Fraction of the hepatic ADH flux appearing as free acetaldehyde (unitless)") # ODE.m dC(25), r(8)/360

    # ---- Scenario scalars (Figs 5-7) --------------------------------------
    fexpr_liver <- fixed(1); label("Hepatic enzyme expression relative to normal (unitless; Fig 5)") # ODE.m aLiver; Fig 5 varies 1, 0.75, 0.5, 0.25, 0.1
    f_aldh2 <- fixed(1); label("ALDH2 isoform activity relative to wild type (unitless; Table 4)") # Table 4 / aALDHtype.m; 1 = ALDH2.1 wild type
    c_disulfiram <- fixed(0); label("Steady blood disulfiram concentration (mg/L)") # main.m mglDisulfiram; Fig 7 uses 0, 2, 4, 6, 8 mg/L
    mw_disulfiram <- fixed(296.539); label("Disulfiram molar mass (g/mol)") # main.m Disulfiram = 1000*mglDisulfiram/296.539
    dsf_c2 <- fixed(0.0015); label("Disulfiram effect on ALDH2 activity, quadratic coefficient (1/uM^2)") # S3 Text Eq S3.1
    dsf_c1 <- fixed(-0.0734); label("Disulfiram effect on ALDH2 activity, linear coefficient (1/uM)") # S3 Text Eq S3.1

    # ---- Genome-scale-model excretion fluxes (ODE.m ModelVer 2) ------------
    # Harvey flux bounds are fractions of totalmetMax = MMMax/0.90, where
    # MMMax is the hepatic ADH rate in mmol/day; ODE.m subtracts each flux
    # /(1440*20) from the tissue concentration. The deposited Fig 4 / Fig 5
    # trajectories show the solver holding urine and sweat at their 10
    # percent upper bounds and breath at its fixed 0.5 percent, so each
    # term is (bound / 0.90 / 20) * hepatic ADH rate. Set all three to 0
    # for the base PBPK (ModelVer 1) that generated Figs 6-8.
    f_urine_wbm <- fixed(0.1 / (0.9 * 20)); label("Urinary ethanol excretion from kidney, fraction of hepatic ADH rate (unitless)") # ODE.m EX_etoh[u] ub = .1*totalmetMax, /(1440*20)
    f_sweat_wbm <- fixed(0.1 / (0.9 * 20)); label("Sweat ethanol excretion from skin, fraction of hepatic ADH rate (unitless)") # ODE.m EX_etoh[sw] ub = .1*totalmetMax, /(1440*20)
    f_breath_wbm <- fixed(0.005 / (0.9 * 20)); label("Breath ethanol excretion from lung, fraction of hepatic ADH rate (unitless)") # ODE.m EX_etoh[br] = .005*totalmetMax, /(1440*20)
  })

  model({
    # ---- 1. Body size, cardiac output (organFlow.m, S1 Text Eq S1.16) -----
    bsa <- 0.007184 * WT^0.425 * HT^0.725
    co <- (159 * bsa - 1.56 * AGE + 114) / 60 # L/min (Eq S1.16 gives L/h)

    # ---- 2. Organ masses (kg) and volumes (L), male (organVolM.m) ---------
    m_fat <- 0.01 * BODYFAT_PCT * WT
    m_blood <- exp(0.067 * bsa - 0.0025 * AGE + 1.7)
    m_brain <- exp(-0.0075 * AGE + 0.0078 * HT - 0.97)
    m_si0 <- 0.45 * 3e-6 * HT^2.49
    m_li0 <- 0.55 * 3e-6 * HT^2.49
    m_heart <- 0.34 * bsa + 0.0018 * AGE - 0.36
    m_kidney <- -0.00038 * AGE + 0.33
    m_liver <- exp(0.87 * bsa - 0.0014 * AGE - 1)
    m_lung <- exp(0.028 * HT + 0.0077 * AGE - 5.6)
    m_muscle <- 17.9 * bsa - 0.0667 * AGE - 1.22
    m_pancreas <- 0.103
    m_skin0 <- exp(-0.0058 * AGE + 1.13)
    m_spleen <- exp(1.13 * bsa - 3.93)
    m_stomach <- 1.05
    # Residual body mass redistributed 10 / 10 / 80 percent (organVolM.m)
    m_rest <- WT - (m_fat + m_blood + m_brain + m_si0 + m_li0 + m_heart +
      m_kidney + m_liver + m_lung + m_muscle + m_pancreas + m_skin0 +
      m_spleen + m_stomach)
    m_si <- m_si0 + 0.1 * m_rest
    m_li <- m_li0 + 0.1 * m_rest
    m_skin <- m_skin0 + 0.8 * m_rest

    # Tissue densities (kg/L), organVolM.m
    v_fat <- m_fat / 0.916
    v_blood <- m_blood / 1.060
    v_brain <- m_brain / 1.035
    v_si <- m_si / 1.044
    v_li <- m_li / 1.044
    v_heart <- m_heart / 1.030
    v_kidney <- m_kidney / 1.050
    v_liver <- m_liver / 1.080
    v_lung <- m_lung / 1.050
    v_muscle <- m_muscle / 1.041
    v_pancreas <- m_pancreas / 1.045
    v_skin <- m_skin / 1.116
    v_spleen <- m_spleen / 1.054
    v_stomach <- m_stomach / 1.050

    # ---- 3. Blood flows (L/min), fractions of cardiac output (organFlow.m) --
    fq_fat <- 0.01 * (0.044 * AGE + 3.9)
    fq_brain <- 0.01 * exp(-0.48 * bsa + 3.5)
    fq_si <- 0.01 * 14 / 2
    fq_li <- 0.01 * 14 / 2
    fq_heart <- 0.01 * (-0.72 * HT + 134)
    fq_kidney <- 0.01 * (-8.7 * bsa + 0.29 * HT - 0.081 * AGE - 13)
    fq_hepart <- 0.01 * 0.24 * (-0.108 * AGE + 27.9) # hepatic artery only
    fq_muscle <- 0.01 * 17.5
    fq_pancreas <- 0.01
    fq_skin <- 0.05
    fq_spleen <- 0.03
    # 'stomach fix': C(14) = 2 - sum(C) with the lung fraction = 1
    fq_stomach <- 1 - (fq_fat + fq_brain + fq_si + fq_li + fq_heart +
      fq_kidney + fq_hepart + fq_muscle + fq_pancreas + fq_skin + fq_spleen)
    q_fat <- fq_fat * co
    q_brain <- fq_brain * co
    q_si <- fq_si * co
    q_li <- fq_li * co
    q_heart <- fq_heart * co
    q_kidney <- fq_kidney * co
    q_hepart <- fq_hepart * co
    q_muscle <- fq_muscle * co
    q_pancreas <- fq_pancreas * co
    q_skin <- fq_skin * co
    q_spleen <- fq_spleen * co
    q_stomach <- fq_stomach * co
    q_liver_in <- q_si + q_li + q_hepart + q_pancreas + q_spleen + q_stomach # Eq S1.9

    # ---- 4. Tissue:blood partition coefficients (organPartition.m) --------
    # K = (P*(Vnl + 0.3*Vph) + Vw + 0.7*Vph) / (same for blood) * fu/fut,
    # tissue composition (Vnl, Vph, Vw) from the Peters textbook; adipose
    # uses fut = 1. Blood: Vnl 0.0035, Vph 0.00225, Vw 0.945.
    pow_e <- 10^logp_etoh
    fut_e <- 1 / (1 + (1 - fu_etoh) * 0.5 / fu_etoh)
    den_e <- pow_e * (0.0035 + 0.3 * 0.00225) + 0.945 + 0.7 * 0.00225
    kp_fat <- (pow_e * (0.79 + 0.3 * 0.002) + 0.180 + 0.7 * 0.002) / den_e * fu_etoh
    kp_brain <- (pow_e * (0.051 + 0.3 * 0.0565) + 0.770 + 0.7 * 0.0565) / den_e * fu_etoh / fut_e
    kp_si <- (pow_e * (0.0487 + 0.3 * 0.0163) + 0.718 + 0.7 * 0.0163) / den_e * fu_etoh / fut_e
    kp_li <- kp_si
    kp_heart <- (pow_e * (0.0115 + 0.3 * 0.0166) + 0.758 + 0.7 * 0.0166) / den_e * fu_etoh / fut_e
    kp_kidney <- (pow_e * (0.0207 + 0.3 * 0.0162) + 0.783 + 0.7 * 0.0162) / den_e * fu_etoh / fut_e
    kp_liver <- (pow_e * (0.0348 + 0.3 * 0.0252) + 0.751 + 0.7 * 0.0252) / den_e * fu_etoh / fut_e
    kp_lung <- (pow_e * (0.003 + 0.3 * 0.009) + 0.811 + 0.7 * 0.009) / den_e * fu_etoh / fut_e
    kp_muscle <- (pow_e * (0.0238 + 0.3 * 0.0072) + 0.760 + 0.7 * 0.0072) / den_e * fu_etoh / fut_e
    kp_pancreas <- (pow_e * (0.0723 + 0.3 * 0.0188) + 0.660 + 0.7 * 0.0188) / den_e * fu_etoh / fut_e
    kp_skin <- (pow_e * (0.0284 + 0.3 * 0.0111) + 0.718 + 0.7 * 0.0111) / den_e * fu_etoh / fut_e
    kp_spleen <- (pow_e * (0.0201 + 0.3 * 0.0198) + 0.788 + 0.7 * 0.0198) / den_e * fu_etoh / fut_e
    kp_stomach <- (pow_e * (0.0338 + 0.3 * 0.0182) + 0.784 + 0.7 * 0.0182) / den_e * fu_etoh / fut_e

    pow_a <- 10^logp_acald
    fut_a <- 1 / (1 + (1 - fu_acald) * 0.5 / fu_acald)
    den_a <- pow_a * (0.0035 + 0.3 * 0.00225) + 0.945 + 0.7 * 0.00225
    kpa_fat <- (pow_a * (0.79 + 0.3 * 0.002) + 0.180 + 0.7 * 0.002) / den_a * fu_acald
    kpa_brain <- (pow_a * (0.051 + 0.3 * 0.0565) + 0.770 + 0.7 * 0.0565) / den_a * fu_acald / fut_a
    kpa_si <- (pow_a * (0.0487 + 0.3 * 0.0163) + 0.718 + 0.7 * 0.0163) / den_a * fu_acald / fut_a
    kpa_li <- kpa_si
    kpa_heart <- (pow_a * (0.0115 + 0.3 * 0.0166) + 0.758 + 0.7 * 0.0166) / den_a * fu_acald / fut_a
    kpa_kidney <- (pow_a * (0.0207 + 0.3 * 0.0162) + 0.783 + 0.7 * 0.0162) / den_a * fu_acald / fut_a
    kpa_liver <- (pow_a * (0.0348 + 0.3 * 0.0252) + 0.751 + 0.7 * 0.0252) / den_a * fu_acald / fut_a
    kpa_lung <- (pow_a * (0.003 + 0.3 * 0.009) + 0.811 + 0.7 * 0.009) / den_a * fu_acald / fut_a
    kpa_muscle <- (pow_a * (0.0238 + 0.3 * 0.0072) + 0.760 + 0.7 * 0.0072) / den_a * fu_acald / fut_a
    kpa_pancreas <- (pow_a * (0.0723 + 0.3 * 0.0188) + 0.660 + 0.7 * 0.0188) / den_a * fu_acald / fut_a
    kpa_skin <- (pow_a * (0.0284 + 0.3 * 0.0111) + 0.718 + 0.7 * 0.0111) / den_a * fu_acald / fut_a
    kpa_spleen <- (pow_a * (0.0201 + 0.3 * 0.0198) + 0.788 + 0.7 * 0.0198) / den_a * fu_acald / fut_a
    kpa_stomach <- (pow_a * (0.0338 + 0.3 * 0.0182) + 0.784 + 0.7 * 0.0182) / den_a * fu_acald / fut_a

    # ---- 5. Absorption / transit rates (S2 Text; ksi sign from main.m) ----
    kstom <- kstom_c2 * drink_frac^2 + kstom_c1 * drink_frac + kstom_c0
    kstomsi <- kstomsi_c2 * drink_frac^2 + kstomsi_c1 * drink_frac + kstomsi_c0
    ksi <- ksi_c2 * drink_frac^2 + ksi_c1 * drink_frac + ksi_c0

    # ---- 6. Concentrations (mmol/L) ----------------------------------------
    c_stom_lumen <- stomach / v_stomach_lumen
    c_si_lumen <- small_intestine_lumen / v_stomach_lumen
    c_fat <- a_fat / v_fat
    c_art <- a_arterial / v_blood
    c_brain <- a_brain / v_brain
    c_si <- a_small_intestine / v_si
    c_li <- a_large_intestine / v_li
    c_heart <- a_heart / v_heart
    c_kidney <- a_kidney / v_kidney
    c_liver <- a_liver / v_liver
    c_lung <- a_lung / v_lung
    c_muscle <- a_muscle / v_muscle
    c_pancreas <- a_pancreas / v_pancreas
    c_skin <- a_skin / v_skin
    c_spleen <- a_spleen / v_spleen
    c_stom <- a_stomach / v_stomach

    ca_fat <- a_fat_acald / v_fat
    ca_art <- a_arterial_acald / v_blood
    ca_brain <- a_brain_acald / v_brain
    ca_si <- a_small_intestine_acald / v_si
    ca_li <- a_large_intestine_acald / v_li
    ca_heart <- a_heart_acald / v_heart
    ca_kidney <- a_kidney_acald / v_kidney
    ca_liver <- a_liver_acald / v_liver
    ca_lung <- a_lung_acald / v_lung
    ca_muscle <- a_muscle_acald / v_muscle
    ca_pancreas <- a_pancreas_acald / v_pancreas
    ca_skin <- a_skin_acald / v_skin
    ca_spleen <- a_spleen_acald / v_spleen
    ca_stom <- a_stomach_acald / v_stomach

    # Portal / hepatic-artery mixing into the liver (Eq S1.10) and venous
    # return to the lung (Eq S1.11). The hepatic artery carries c_art
    # without a partition coefficient.
    c_liver_in <- (q_si * c_si / kp_si + q_li * c_li / kp_li + q_hepart * c_art +
      q_pancreas * c_pancreas / kp_pancreas + q_spleen * c_spleen / kp_spleen +
      q_stomach * c_stom / kp_stomach) / q_liver_in
    c_ven <- (q_fat * c_fat / kp_fat + q_brain * c_brain / kp_brain +
      q_heart * c_heart / kp_heart + q_kidney * c_kidney / kp_kidney +
      q_liver_in * c_liver / kp_liver + q_muscle * c_muscle / kp_muscle +
      q_skin * c_skin / kp_skin) / co
    ca_liver_in <- (q_si * ca_si / kpa_si + q_li * ca_li / kpa_li + q_hepart * ca_art +
      q_pancreas * ca_pancreas / kpa_pancreas + q_spleen * ca_spleen / kpa_spleen +
      q_stomach * ca_stom / kpa_stomach) / q_liver_in
    ca_ven <- (q_fat * ca_fat / kpa_fat + q_brain * ca_brain / kpa_brain +
      q_heart * ca_heart / kpa_heart + q_kidney * ca_kidney / kpa_kidney +
      q_liver_in * ca_liver / kpa_liver + q_muscle * ca_muscle / kpa_muscle +
      q_skin * ca_skin / kpa_skin) / co

    # ---- 7. Enzyme kinetics (maxrates.m, ODE.m) ----------------------------
    # Age factor on ADH, applied to every Michaelis-Menten rate (male branch)
    if (AGE <= 25) {
      f_age <- 1
    } else if (AGE > 55) {
      f_age <- 0.5
    } else {
      f_age <- -0.00102 * AGE^2 + 0.067 * AGE - 0.0663
    }
    r_adh_liver <- f_age * vmax_adh_liver * c_liver / (km_adh_liver + c_liver)
    r_adh_stomach <- f_age * vmax_adh_stomach * c_stom / (km_adh_stomach + c_stom)
    r_aldh_liver <- f_age * vmax_aldh_liver * ca_liver / (km_aldh_liver + ca_liver)
    # ALDH2 activity: isoform (Table 4) times disulfiram inhibition (Eq S3.1)
    dsf_um <- 1000 * c_disulfiram / mw_disulfiram
    f_aldh <- f_aldh2 * (dsf_c2 * dsf_um^2 + dsf_c1 * dsf_um + 1)

    # ---- 8. ODEs: ethanol (amounts, mmol; Eqs S1.1-S1.15) -----------------
    d/dt(stomach) <- -(kstom + kstomsi) * stomach
    d/dt(small_intestine_lumen) <- kstomsi * stomach - ksi * small_intestine_lumen
    d/dt(a_fat) <- q_fat * (c_art - c_fat / kp_fat)
    # Eq S1.2 as printed and coded: arterial concentration relaxes at the
    # lung flow / lung volume rate, so blood volume only scales the amount.
    d/dt(a_arterial) <- v_blood * co / v_lung * (c_lung / kp_lung - c_art)
    d/dt(a_brain) <- q_brain * (c_art - c_brain / kp_brain)
    # Luminal concentration is added 1:1 to tissue concentration (Eqs S1.6,
    # S1.7), so the absorbed amount is scaled by tissue volume.
    d/dt(a_small_intestine) <- q_si * (c_art - c_si / kp_si) + v_si * ksi * c_si_lumen
    d/dt(a_large_intestine) <- q_li * (c_art - c_li / kp_li)
    d/dt(a_heart) <- q_heart * (c_art - c_heart / kp_heart)
    d/dt(a_kidney) <- q_kidney * (c_art - c_kidney / kp_kidney) - v_kidney * f_urine_wbm * r_adh_liver
    d/dt(a_liver) <- q_liver_in * (c_liver_in - c_liver / kp_liver) - v_liver * fexpr_liver^2 * r_adh_liver
    d/dt(a_lung) <- co * (c_ven - c_lung / kp_lung) - v_lung * f_breath_wbm * r_adh_liver
    d/dt(a_muscle) <- q_muscle * (c_art - c_muscle / kp_muscle)
    d/dt(a_pancreas) <- q_pancreas * (c_art - c_pancreas / kp_pancreas)
    d/dt(a_skin) <- q_skin * (c_art - c_skin / kp_skin) - v_skin * f_sweat_wbm * r_adh_liver
    d/dt(a_spleen) <- q_spleen * (c_art - c_spleen / kp_spleen)
    d/dt(a_stomach) <- q_stomach * (c_art - c_stom / kp_stomach) + v_stomach * kstom * c_stom_lumen -
      v_stomach * r_adh_stomach

    # ---- 9. ODEs: acetaldehyde ---------------------------------------------
    d/dt(a_fat_acald) <- q_fat * (ca_art - ca_fat / kpa_fat)
    d/dt(a_arterial_acald) <- v_blood * co / v_lung * (ca_lung / kpa_lung - ca_art)
    d/dt(a_brain_acald) <- q_brain * (ca_art - ca_brain / kpa_brain)
    d/dt(a_small_intestine_acald) <- q_si * (ca_art - ca_si / kpa_si)
    d/dt(a_large_intestine_acald) <- q_li * (ca_art - ca_li / kpa_li)
    d/dt(a_heart_acald) <- q_heart * (ca_art - ca_heart / kpa_heart)
    d/dt(a_kidney_acald) <- q_kidney * (ca_art - ca_kidney / kpa_kidney)
    d/dt(a_liver_acald) <- q_liver_in * (ca_liver_in - ca_liver / kpa_liver) +
      v_liver * fexpr_liver^2 * f_acald * r_adh_liver - v_liver * f_aldh * fexpr_liver * r_aldh_liver
    d/dt(a_lung_acald) <- co * (ca_ven - ca_lung / kpa_lung)
    d/dt(a_muscle_acald) <- q_muscle * (ca_art - ca_muscle / kpa_muscle)
    d/dt(a_pancreas_acald) <- q_pancreas * (ca_art - ca_pancreas / kpa_pancreas)
    d/dt(a_skin_acald) <- q_skin * (ca_art - ca_skin / kpa_skin)
    d/dt(a_spleen_acald) <- q_spleen * (ca_art - ca_spleen / kpa_spleen)
    d/dt(a_stomach_acald) <- q_stomach * (ca_art - ca_stom / kpa_stomach)

    # ---- 10. Dose and observation ------------------------------------------
    f(stomach) <- fstomach

    Cc <- c_art # blood ethanol (mmol/L)
    Cc_acald <- ca_art # blood acetaldehyde (mmol/L; x1000 for the paper's umol/L)
  })
}
