Clegg_2024_sodiumZirconiumCyclosilicate <- function() {
  description <- paste(
    "QSP. Four-compartment potassium-homeostasis model (gastrointestinal",
    "tract, blood, tissue extracellular fluid, tissue intracellular space)",
    "with reversible potassium binding to orally administered sodium",
    "zirconium cyclosilicate (SZC) in the gastrointestinal lumen, faecal",
    "loss of free and bound SZC, and eGFR-scaled renal potassium clearance",
    "with a sinusoidal diurnal rhythm, in adults with hyperkalaemia",
    "(Clegg 2024). An optional hemodialysis extension removes blood",
    "potassium into a fixed-concentration dialysate bath during sessions.",
    "Typical-value model: parameters are the medians of the per-subject",
    "fits; there is no IIV or residual error."
  )
  reference <- paste(
    "Clegg LE, Chu L, Nagard M, Boulton DW, Penland RC. Potassium",
    "homeostasis and therapeutic intervention with sodium zirconium",
    "cyclosilicate: A model-informed drug development case study.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(2):296-307.",
    "doi:10.1002/psp4.13084"
  )
  vignette <- "Clegg_2024_sodiumZirconiumCyclosilicate"

  # Potassium and SZC states of the Clegg 2024 Figure 1 four-compartment
  # model. None has a canonical nlmixr2lib role (the "drug" SZC never leaves
  # the gut, and the measured species is endogenous potassium), so the names
  # describe the species and the paper's compartment:
  #   gi_k         <- K_GI      (potassium in the lumped GI tract)
  #   blood_k      <- K_B       (potassium in plasma / serum)
  #   tissue_ecf_k <- K_E       (potassium in tissue extracellular fluid)
  #   tissue_icf_k <- K_I       (potassium in the tissue intracellular space)
  #   gi_szc       <- SZC_GI    (free SZC in the GI tract)
  #   gi_kszc      <- K.SZC_GI  (potassium-SZC complex in the GI tract)
  paper_specific_compartments <- c(
    "gi_k",
    "blood_k",
    "tissue_ecf_k",
    "tissue_icf_k",
    "gi_szc",
    "gi_kszc"
  )

  # SZC doses (g) go into gi_szc; meals (mmol K+) go into gi_k.
  dosing <- c("gi_szc", "gi_k")

  units <- list(
    time = "h",
    dosing = "g (SZC into gi_szc) and mmol (dietary potassium into gi_k; 1 mmol = 1 mEq for K+)",
    concentration = "mmol/L (Cc is serum potassium, sK+)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Scales the blood, tissue extracellular and tissue intracellular volumes proportionally (Clegg 2024 Table S1: total body water = 60% of weight at 1 kg/L). The GI-tract volume does not scale with weight.",
      source_name = "weight"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (eGFR)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Scales renal potassium clearance linearly as r_renal * eGFR / 60 (Clegg 2024 Model equations, dK_B/dt). The paper does not name the eGFR equation; Table S2 reports the virtual-population mean 44.11 (SD 25.26).",
      source_name = "eGFR"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description = "Hemodialysis session running (time-varying)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no dialysis running)",
      notes = "1 during a hemodialysis session, when blood potassium transfers to the fixed-concentration dialysate bath at rate r_D * (K_B - K_dialysate) (Clegg 2024 Figure 1 'Dialysis bath' and Figure 5). The Clegg 2024 hemodialysis application ran 4-h sessions on Monday, Wednesday and Friday from 8 AM to noon. 0 for the base model.",
      source_name = NA_character_
    ),
    RRT_HEMODIAL_STATUS = list(
      description = "Subject on chronic hemodialysis (hemodialysis application of the model)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no chronic hemodialysis; base model)",
      notes = "1 selects the Clegg 2024 hemodialysis application, which holds intracellular potassium at its steady-state value ('Fix intracellular K+ to steady-state value to prevent progressive K+ loss', Figure 5 flow chart). Set it for the whole simulation, not only during sessions.",
      source_name = NA_character_
    )
  )

  # The GI-lumen states hold what was swallowed (meals, SZC) and what is
  # left for faecal excretion; "administration site" is the closest member
  # of the specimen vocabulary for gastrointestinal luminal content.
  compartmentData <- list(
    gi_k = list(analyte = "potassium", units = "mmol", specimen = "administration site", verified = TRUE),
    blood_k = list(analyte = "potassium", units = "mmol", specimen = "serum", verified = TRUE),
    tissue_ecf_k = list(analyte = "potassium", units = "mmol", specimen = "tissue", verified = TRUE),
    tissue_icf_k = list(analyte = "potassium", units = "mmol", specimen = "tissue", verified = TRUE),
    gi_szc = list(
      analyte = "sodium zirconium cyclosilicate",
      units = "g",
      specimen = "administration site",
      verified = TRUE
    ),
    gi_kszc = list(
      analyte = "potassium-sodium zirconium cyclosilicate complex",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 820,
    n_studies = 2,
    age_range = "mean 65.27 (SD 12.24) years",
    age_median = "not reported (mean 65.27 years)",
    weight_range = "mean 89.19 (SD 22.02) kg",
    weight_median = "not reported (mean 89.19 kg)",
    sex_female_pct = 100 * (820 - 502) / 820,
    race_ethnicity = c(White = 100 * 690 / 820, Black = 100 * 103 / 820, Other = 100 * 27 / 820),
    disease_state = "Outpatients with hyperkalaemia (ZS-003: serum K+ 5-6.5 mmol/L; ZS-004/HARMONIZE: serum K+ >= 5.1 mmol/L). Type 2 diabetes 43%, congestive heart failure 10%, chronic kidney disease 55%; ACE inhibitor 42%, ARB 18%, MRA 6%, loop diuretic 28%. eGFR mean 44.11 (SD 25.26) mL/min/1.73 m^2.",
    dose_range = "Oral SZC 1.25-10 g three times daily for 48 h (acute phase), then 0-15 g once daily (maintenance: 12 days in ZS-003, 28 days in ZS-004), or placebo.",
    regions = "Multinational (phase III trials ZS-003, NCT01737697, and ZS-004/HARMONIZE, NCT02088073)",
    renal_function = "eGFR mean 44.11 (SD 25.26) mL/min/1.73 m^2",
    notes = "Virtual population of the 820 subjects with accepted per-subject fits (652 of 753 in ZS-003, 168 of 258 in ZS-004); Clegg 2024 Table S2. Each subject was fitted separately by least squares with physiological penalties; Table 1 reports the median (IQR) of the per-subject estimates. Applications: emergency treatment of severe hyperkalaemia (ENERGIZE design) and hemodialysis (DIALIZE design)."
  )

  ini({
    # Per-subject rate constants: Clegg 2024 Table 1, median of the 820
    # individual estimates (IQR in the trailing comment). The paper uses r
    # for rate constants to avoid confusion with potassium, K.
    lr_gb <- log(0.062); label("K+ transport GI tract to blood r_GB (1/h)") # Table 1 r_GB median 0.062 (IQR 0.041, 0.11)
    lr_bg <- log(0.23); label("K+ transport blood to GI tract r_BG (1/h)") # Table 1 r_BG median 0.23 (IQR 0.09, 3.8)
    lr_eb <- log(3.9); label("K+ transport tissue extracellular fluid to blood r_EB (1/h)") # Table 1 r_EB median 3.9 (IQR 0.14, 49.9)
    r_be_eb <- fixed(2.3); label("Ratio r_BE / r_EB, blood-to-tissue-ECF over tissue-ECF-to-blood rate (unitless; volume ratio V_E / V_B)") # Table 1 r_BE row 'r_BE = 2.3*r_EB (volume ratio)'; Methods
    lr_ei <- log(0.091); label("K+ uptake by cells in tissue r_EI (1/h)") # Table 1 r_EI median 0.091 (IQR 0.069, 0.12)
    lr_ie <- log(0.0011); label("K+ release by cells in tissue r_IE (1/h)") # Table 1 r_IE median 0.0011 (IQR 0.00083, 0.0015)
    lr_cl_k <- log(0.0065); label("Faecal K+ clearance from the GI tract r_CL,K (1/h)") # Table 1 r_CL,K median 0.0065 (IQR 0.0036, 0.011)
    lr_cl_szc <- log(0.42); label("Faecal clearance of free and K+-bound SZC r_CL,SZC (1/h)") # Table 1 r_CL,SZC median 0.42 (IQR 0.26, 0.79)
    lr_renal <- log(0.18); label("Renal K+ clearance from blood r_renal at eGFR 60 mL/min/1.73 m^2 (1/h)") # Table 1 r_renal median 0.18 (IQR 0.13, 0.29)
    p_diurnal <- 0.20; label("Magnitude p of the diurnal variation in renal K+ clearance (unitless, 0-1)") # Table 1 p median 0.20 (IQR 0.073, 0.49)

    # K+-SZC binding: the same for all subjects. Figure S2c small-intestine
    # (pH 4.5) row k_on = 69.6 /(s*M), k_off = 2.14 /s, then r_on was
    # increased six-fold in vivo (Methods 'Model evaluation and tuning';
    # Figure 2). Unit conversion: 69.6 /(M*s) * 3600 s/h / 1000 mM/M =
    # 250.56 /(mM*h); times 6 = 1503.36 /(mM*h). 2.14 /s * 3600 = 7704 /h.
    lr_on <- fixed(log(1503.36)); label("K+-SZC binding rate constant r_on (1/(mM*h))") # Figure S2c k_on 69.6 /(s*M) x 6 (Methods; Figure 2)
    lr_off <- fixed(log(7704)); label("K+-SZC unbinding rate constant r_off (1/h)") # Figure S2c k_off 2.14 /s
    mw_szc <- fixed(410); label("Molecular weight per K+ binding site used to convert SZC mass to moles (g/mol)") # Figure S2 legend 'assuming a ZS-9 molecular weight of 410 g/mol'

    # Physiological volumes: Clegg 2024 Table S1 (fractions of body weight at
    # a density of 1 kg/L) and Figure 1.
    f_vb <- fixed(0.048); label("Blood (plasma) volume per kg body weight (L/kg)") # Table S1 V_B = 20% of ECF = 0.2 x 0.4 x 0.6; Figure 1 'V_B = 4.8% of weight'
    f_ve <- fixed(0.110769); label("Tissue extracellular fluid volume per kg body weight (L/kg)") # Table S1 V_E = 8/13 of TIF, TIF = 75% of ECF, ECF = 40% of TBW, TBW = 60% of weight: 0.6 x 0.4 x 0.75 x 8/13 = 0.110769
    f_vi <- fixed(0.36); label("Tissue intracellular volume per kg body weight (L/kg)") # Table S1 V_I = 60% of TBW = 0.6 x 0.6; Figure 1 'V_I = 36% of weight'
    vgi <- fixed(3); label("GI tract volume, not weight-scaled (L)") # Table S1 V_GI '3 L', not weight-scaled; Figure 1

    # Dietary potassium. Meals are dosing events into gi_k (25% at 8 AM, 33%
    # at noon, 42% at 6 PM; Methods 'Information sources'); this constant
    # intake only sets the steady-state initial conditions (Figure S1).
    k_intake <- fixed(100); label("Daily dietary K+ intake used for the steady-state initial conditions (mmol/day)") # Methods 'K+ consumption was assumed to be 100 mEq/day for all subjects'

    # Hemodialysis extension (inactive unless RRT_HEMODIAL_ACTIVE = 1).
    lr_d <- fixed(log(1.05)); label("K+ transfer rate from blood to the dialysate bath r_D (1/h)") # Figure 5 flow chart 'r_D = 1.05 hr-1'
    k_dialysate <- fixed(2); label("Dialysate K+ concentration (mmol/L)") # Results 'dialysate K+ of 2 mmol/L'; Figure 5 flow chart
  })

  model({
    # Individual parameters
    r_gb <- exp(lr_gb)
    r_bg <- exp(lr_bg)
    r_eb <- exp(lr_eb)
    r_be <- r_be_eb * r_eb
    r_ei <- exp(lr_ei)
    r_ie <- exp(lr_ie)
    r_cl_k <- exp(lr_cl_k)
    r_cl_szc <- exp(lr_cl_szc)
    r_renal <- exp(lr_renal)
    r_on <- exp(lr_on)
    r_off <- exp(lr_off)
    r_d <- exp(lr_d)

    # Compartment volumes (L)
    vb <- f_vb * WT
    ve <- f_ve * WT
    vi <- f_vi * WT

    # Renal K+ clearance: r_f = r_renal * eGFR / 60 (Model equations). m(t)
    # = 1 - p * cos(pi * t_d / 12) is lowest at midnight and highest at noon;
    # cos() has a 24-h period, so model time t is used directly. Model time 0
    # must therefore be a midnight.
    r_f <- r_renal * CRCL / 60
    m_diurnal <- 1 - p_diurnal * cos(pi * t / 12)

    # Steady-state initial conditions under constant K+ intake and no SZC
    # (Figure S1), written as amounts (concentration x volume).
    r_in <- k_intake / 24
    kb_ss <- (r_in / vb * r_gb / (r_cl_k + r_gb)) /
      (r_bg + r_f - r_gb * r_bg / (r_cl_k + r_gb))
    gi_k(0) <- (r_in + r_bg * vb * kb_ss) / (r_cl_k + r_gb)
    blood_k(0) <- kb_ss * vb
    tissue_ecf_k(0) <- kb_ss * vb * r_be / r_eb
    tissue_icf_k(0) <- kb_ss * vb * r_be / r_eb * r_ei / r_ie

    # Concentrations (mmol/L). SZC mass (g) becomes moles of K+ binding
    # sites at mw_szc g/mol, one K+ per site.
    k_gi <- gi_k / vgi
    k_b <- blood_k / vb
    szc_gi <- gi_szc / mw_szc * 1000 / vgi
    kszc_gi <- gi_kszc / vgi

    # Net K+-SZC binding flux in the GI tract (mmol/h)
    bind <- (r_on * k_gi * szc_gi - r_off * kszc_gi) * vgi

    # Hemodialysis: one-way transfer from blood to the fixed-[K+] dialysate
    # bath during a session (mmol/h); intracellular K+ held at its steady
    # state in the hemodialysis application.
    dialysis <- RRT_HEMODIAL_ACTIVE * r_d * (k_b - k_dialysate) * vb
    icf_free <- 1 - RRT_HEMODIAL_STATUS

    # Clegg 2024 Model equations, each multiplied through by its compartment
    # volume so the states are amounts.
    d/dt(gi_k) <- -r_cl_k * gi_k - bind - r_gb * gi_k + r_bg * blood_k
    d/dt(blood_k) <- -r_f * m_diurnal * blood_k + r_gb * gi_k - r_bg * blood_k +
      r_eb * tissue_ecf_k - r_be * blood_k - dialysis
    d/dt(tissue_ecf_k) <- -r_eb * tissue_ecf_k + r_be * blood_k -
      r_ei * tissue_ecf_k + r_ie * tissue_icf_k
    d/dt(tissue_icf_k) <- icf_free * (r_ei * tissue_ecf_k - r_ie * tissue_icf_k)
    d/dt(gi_szc) <- -r_cl_szc * gi_szc - bind * mw_szc / 1000
    d/dt(gi_kszc) <- -r_cl_szc * gi_kszc + bind

    # Outputs (mmol/L)
    k_e <- tissue_ecf_k / ve
    k_i <- tissue_icf_k / vi
    Cc <- k_b
  })
}
