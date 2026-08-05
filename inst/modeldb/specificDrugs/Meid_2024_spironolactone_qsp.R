Meid_2024_spironolactone_qsp <- function() {
  description <- paste(
    "QSP. Whole-body plasma-potassium homeostasis with renal handling",
    "along the nephron and aldosterone feedback, pharmacologically",
    "modulated by the mineralocorticoid-receptor antagonist",
    "spironolactone. Sixteen ODE states: a three-transit oral",
    "spironolactone absorption chain plus a spironolactone central",
    "concentration, a two-compartment canrenone (active metabolite)",
    "disposition, an extracellular and an intracellular potassium pool,",
    "luminal potassium amounts and principal-cell potassium",
    "concentrations in the distal convoluted tubule (DCT), connecting",
    "tubule / cortical collecting duct (CNT/CCD) segments, cumulative",
    "urinary potassium excretion, and a mineralocorticoid-receptor",
    "occupancy state. Potassium is filtered at the glomerulus,",
    "reabsorbed in the proximal tubule / loop of Henle in proportion to",
    "sodium reabsorption, secreted into DCT / CNT / CCD lumen by",
    "Goldman-Hodgkin-Katz passive flux coupled to Na/K-ATPase active",
    "basolateral uptake, and partly reabsorbed in the medullary",
    "collecting duct. Aldosterone rises exponentially with plasma",
    "potassium and falls with sodium intake, and scales luminal",
    "potassium permeability; canrenone inhibits that pathway through an",
    "Emax mineralocorticoid-receptor occupancy model. Five parameters",
    "carry inter-individual variability and were the ones the source",
    "estimated per patient by Bayesian maximum-a-posteriori updating",
    "from electronic health records: potassium intake, sodium intake,",
    "extracellular fluid volume, mineralocorticoid-receptor abundance,",
    "and a hyperaldosteronism effect. The structural model is inherited",
    "unchanged from Maddah and Hallow 2022; Meid 2024 Appendix C prints",
    "the complete operationalized rxode2 source that this file",
    "transcribes.",
    sep = " "
  )
  reference <- paste(
    "Meid AD, Scherkl C, Metzner M, Czock D, Seidling HM (2024).",
    "Real-World Application of a Quantitative Systems Pharmacology (QSP)",
    "Model to Predict Potassium Concentrations from Electronic Health",
    "Records: A Pilot Case towards Prescribing Monitoring of",
    "Spironolactone. Pharmaceuticals 17(8):1041.",
    "doi:10.3390/ph17081041.",
    "Structural model and physiological constants inherited from Maddah",
    "E, Hallow KM (2022), A quantitative systems pharmacology model of",
    "plasma potassium regulation by the kidney and aldosterone,",
    "J Pharmacokinet Pharmacodyn 49:471-486,",
    "doi:10.1007/s10928-022-09815-x; every value used here is taken from",
    "the Meid 2024 Appendix C listing, not from the upstream paper",
    "(which is not open access).",
    "Spironolactone / canrenone disposition traces to Gardiner 1989",
    "(J Clin Pharmacol 29:342-347) as cited by Meid 2024 reference 15.",
    sep = " "
  )
  vignette <- "Meid_2024_spironolactone_qsp"

  # Mechanistic renal / potassium-homeostasis states that have no
  # canonical nlmixr2lib compartment role. Names are lower-cased
  # transcriptions of the Meid 2024 Appendix C state names so the file
  # can be audited line-by-line against the published listing:
  #   k_ecf                <- K
  #   intracellular_k      <- intracellular_K
  #   dct_luminal_k_amount <- DCT_luminal_K_amount   (CNT / CCD likewise)
  #   dct_cell_k_conc      <- DCT_cell_K_conc        (CNT / CCD likewise)
  #   potassium_excretion  <- potassium_excretion_rate
  #   mra_effect           <- MRA_effect
  paper_specific_compartments <- c(
    "k_ecf", "intracellular_k",
    "dct_luminal_k_amount", "cnt_luminal_k_amount", "ccd_luminal_k_amount",
    "dct_cell_k_conc", "cnt_cell_k_conc", "ccd_cell_k_conc",
    "potassium_excretion", "mra_effect"
  )

  units <- list(
    time          = "min",
    dosing        = "ug (oral spironolactone: enter a 25 mg tablet as amt = 25000. The Appendix C listing does not state its amount unit, but it is fixed by EC50_spiro = 1.8296 being on the same scale as the canrenone state: micrograms is the only scale on which the model both reproduces published canrenone exposure and produces a spironolactone effect at all -- see the DOSE UNIT note in model())",
    concentration = "mEq/mL (Cc is plasma potassium, exactly as Meid 2024 Appendix C sets `Cc = plasma_K`; multiply by 1000 for the clinical mmol/L scale, so the physiological normal 0.0042 mEq/mL is 4.2 mmol/L). The spironolactone and canrenone states are in ng/mL"
  )

  compartmentData <- list(
    depot = list(
      analyte = "spironolactone", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    transit1 = list(
      analyte = "spironolactone", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    transit2 = list(
      analyte = "spironolactone", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "spironolactone", units = "ng/mL (the source carries this state as a CONCENTRATION -- d/dt is divided by V1_spiro)",
      specimen = "plasma", verified = TRUE
    ),
    central_canrenone = list(
      analyte = "canrenone", units = "ng/mL (concentration state; d/dt divided by V_canrenone)",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_canrenone = list(
      analyte = "canrenone", units = "ng/mL (concentration state; d/dt divided by V2_canrenone)",
      specimen = "plasma", verified = TRUE
    ),
    k_ecf = list(
      analyte = "potassium", units = "mEq",
      specimen = "plasma", verified = TRUE
    ),
    intracellular_k = list(
      analyte = "potassium", units = "mEq",
      specimen = "tissue", verified = TRUE
    ),
    dct_luminal_k_amount = list(
      analyte = "potassium", units = "mEq (single-nephron distal convoluted tubule lumen)",
      specimen = "urine", verified = TRUE
    ),
    cnt_luminal_k_amount = list(
      analyte = "potassium", units = "mEq (single-nephron connecting tubule lumen)",
      specimen = "urine", verified = TRUE
    ),
    ccd_luminal_k_amount = list(
      analyte = "potassium", units = "mEq (single-nephron cortical collecting duct lumen)",
      specimen = "urine", verified = TRUE
    ),
    dct_cell_k_conc = list(
      analyte = "potassium", units = "mEq/mL (distal convoluted tubule principal-cell concentration)",
      specimen = "tissue", verified = TRUE
    ),
    cnt_cell_k_conc = list(
      analyte = "potassium", units = "mEq/mL (connecting tubule principal-cell concentration)",
      specimen = "tissue", verified = TRUE
    ),
    ccd_cell_k_conc = list(
      analyte = "potassium", units = "mEq/mL (cortical collecting duct principal-cell concentration)",
      specimen = "tissue", verified = TRUE
    ),
    potassium_excretion = list(
      analyte = "potassium", units = "mEq (cumulative whole-body urinary excretion; the source names this state potassium_excretion_rate but integrates a rate, so it holds an AMOUNT)",
      specimen = "urine", verified = TRUE
    ),
    mra_effect = list(
      analyte = "mineralocorticoid receptor", units = "fraction of receptor NOT inhibited (1 = no antagonist present)",
      specimen = "not applicable", verified = TRUE
    )
  )

  covariateData <- list(
    SOD = list(
      description        = "Plasma (serum) sodium concentration, used as the model's `norm_plasma_Na` -- the sodium concentration of glomerular filtrate. Meid 2024 Section 4.1.2: 'sodium plasma concentrations and GFR estimates were used as covariate information in the model'; Appendix C sets `norm_plasma_Na = actual_norm_plasma_Na`. Pilot cohort mean 139.22 (SD 3.67) mmol/L, range 133-143 (Table A2).",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters as `filtered_Na = SNGFR * (SOD / 1000)`, i.e. the model divides by 1000 to reach mEq/mL. Last-observation-carried-forward was used for missing values in the source (Section 4.1.2). Physiological reference is 140 mmol/L: that is the value at which the printed Appendix C initial conditions form an exact steady state (see the derivation note on nephrons_ref in ini()).",
      source_name        = "actual_norm_plasma_Na"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (CKD-EPI) as reported in the electronic health record. Meid 2024 Table A2 reports eGFR (CKD-EPI) in mL/min/1.73 m^2; pilot cohort mean 72 (SD 28), range 27-103. Appendix C sets `GFR = actual_GFR` and consumes it directly as an absolute filtration rate in mL/min without BSA de-normalization.",
      units              = "mL/min/1.73 m^2 (consumed by the model as mL/min)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "CRCL does two things. (1) It sets the glomerular filtration rate. (2) Per Meid 2024 Section 2 ('The associated nephron number must be adjusted proportionally if GFR is used as a covariate from patient data'), it scales the nephron count: number_of_nephrons = nephrons_ref * CRCL / gfr_ref. That proportional scaling holds single-nephron GFR constant, which is why the source found that 'changing the GFR in the model alone did not meaningfully influence potassium predictions' -- the renal-function signal reaches plasma potassium through the nephron count, not through SNGFR. Last-observation-carried-forward was used for missing values in the source.",
      source_name        = "actual_GFR (and, via the proportional rule, actual_number_of_nephron)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 9L,
    n_studies      = 1L,
    age_range      = "56-88 years (mean 68, SD 12)",
    weight_range   = "54-104 kg (mean 87, SD 14)",
    sex_female_pct = 22,
    race_ethnicity = NA_character_,
    disease_state  = "Hospitalised adults newly initiating spironolactone at a German tertiary-care hospital (Heidelberg University Hospital, 2500 beds, 17 wards), admitted 1 January - 31 March 2023. Median length of stay 10.2 days. Comorbidities (Elixhauser): congestive heart failure 67%, cardiac arrhythmias 56%, uncomplicated hypertension 44%, hypothyroidism 22%, liver disease 22%, fluid and electrolyte disorders 22%, alcohol abuse 22%, uncomplicated diabetes 11%, coagulopathy 11%. Baseline eGFR (CKD-EPI) mean 72 (SD 28) mL/min/1.73 m^2, range 27-103; baseline sodium 139.22 (SD 3.67) mmol/L; baseline potassium 4.12 (SD 0.34) mmol/L. Patients with potassium binders, potassium infusions, or incomplete spironolactone dosing schedules were excluded, as were patients with other immediate risk factors for potassium derailment (haemolytic anaemia, renal impairment, diarrhoea).",
    dose_range     = "Oral spironolactone 25-100 mg (median 25 mg): 25 mg in 6 patients, 50 mg in 1, 50/100 mg in 1, 100 mg in 1. Four of nine also received oral potassium supplementation. Co-medication: ACE inhibitors 44%, angiotensin receptor blockers 11%, high-ceiling diuretics 56%, low-ceiling diuretics 11%.",
    regions        = "Germany",
    notes          = "Inter-individual variability was estimated by FOCEi non-linear mixed-effects fitting in a SEPARATE, earlier sample of 20 randomly selected patients who initiated spironolactone in hospital (10 of whom also received oral potassium supplementation) -- Meid 2024 Section 4.4. The nine patients described here are the pilot validation sample used for the Bayesian a-posteriori prediction assessment (Table 1: average fold error 1.06, absolute average fold error 1.19, percent prediction error 7.3% [95% CI 5.6; 9.0] over all nine). Only sodium and eGFR entered the model as covariates; the demographic and comorbidity fields above describe the cohort but are not model inputs."
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight (Table A2, mean 87 SD 14 kg). Reported to characterise the pilot cohort; no allometric or body-size scaling appears anywhere in the Appendix C model listing.",
      units       = "kg",
      type        = "continuous",
      notes       = "Documented-but-unused: the QSP model's extracellular fluid volume is estimated per patient rather than scaled from weight."
    ),
    AGE = list(
      description = "Age (Table A2, mean 68 SD 12 years).",
      units       = "years",
      type        = "continuous",
      notes       = "Documented-but-unused: cohort description only; not referenced in model()."
    ),
    SEXF = list(
      description = "Sex (Table A2; 22% female).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Documented-but-unused: cohort description only; not referenced in model()."
    )
  )

  ini({
    # ===================================================================
    # ESTIMATED / INDIVIDUALISED PARAMETERS
    #
    # These are the five parameters that Meid 2024 selected by local and
    # global sensitivity analysis plus identifiability analysis (Section
    # 2, Figures 1 and A4-A5) and then updated per patient by Bayesian
    # maximum-a-posteriori estimation. The paper reports their IIV and
    # their physiologically plausible ranges but does not print the
    # population typical values, so each typical value below is derived
    # from the reference state that Appendix C itself prints (the
    # function-argument initial conditions and the "norm_" / "0"-suffixed
    # reference constants). Each derivation is verified against a printed
    # constant, and each derived value falls inside the paper's own
    # global-sensitivity-analysis range (Section 4.3.2).
    # ===================================================================

    # Derivation: the reference state that Appendix C prints is an exact
    # steady state of the potassium system, and pinning it recovers Kin.
    # With number_of_nephrons = 2e6 (see nephrons_ref below), the
    # printed CCD initial condition gives a cortical-collecting-duct
    # potassium delivery of
    #   N * CCD_luminal_K_amount(0) / CCD_volume * CCD_water_out
    #     = 0.0945867 mEq/min,
    # which reproduces the printed constant potassium_in_MCD0 =
    # 0.0945866666666666 to seven digits. Subtracting the printed
    # medullary reabsorption N * K_reabsorption_MCD_rate0 = 0.0145867
    # leaves urinary excretion = 0.08000 mEq/min, and at steady state
    # d/dt(K) = 0 forces Kin = that excretion. Equivalently
    # potassium_in_MCD0 * (1 - eta_MCD0) = 0.0945866666666666 *
    # (1 - 0.154214829433323) = 0.08 exactly. 0.08 mEq/min = 115 mEq/day,
    # and lies inside the paper's global-SA range [0.073; 0.084] mEq/min.
    lkin <- log(0.08); label("Dietary potassium intake Kin (mEq/min)") # derived from Meid 2024 Appendix C potassium_in_MCD0 and eta_MCD0 (= 0.08 exactly); inside the Section 4.3.2 range [0.073; 0.084]

    # Derivation: Appendix C prints norm_Na_intake = 0.0694444444444444
    # mEq/min (= 100 mEq/day) and builds the aldosterone equation as
    # max(0, exp(-m_Na_ALDO * (Nain - norm_Na_intake)) - 1), which is
    # exactly zero -- i.e. aldosterone sits at norm_Aldo -- only when
    # Nain = norm_Na_intake. Independently, Nain = 0.0694444444444444
    # with 2e6 nephrons reproduces two more printed constants exactly:
    # Na_out_loH0 = Nain / 0.05 = 1.38888888888889, and the medullary
    # sodium reabsorption N * (Na_out_LoH * 0.5 - Nain / N) = 0.625 =
    # norm_Na_reabsorption_MCD. Inside the Section 4.3.2 range
    # [0.01; 0.17] mEq/min.
    lnain <- log(0.0694444444444444); label("Dietary sodium intake Nain (mEq/min)") # derived from Meid 2024 Appendix C norm_Na_intake / Na_out_loH0 / norm_Na_reabsorption_MCD; inside the Section 4.3.2 range [0.01; 0.17]

    # Derivation: Appendix C prints K_init = 63 mEq and norm_plasma_K =
    # 0.0042 mEq/mL, and defines plasma_K = K / V_ecf. The printed
    # initial condition therefore sits at the physiological normal only
    # when V_ecf = 63 / 0.0042 = 15000 mL, the textbook adult
    # extracellular fluid volume. Inside the Section 4.3.2 range
    # [10,000; 25,000] mL.
    lv_ecf <- log(15000); label("Extracellular fluid volume V_ecf (mL)") # derived from Meid 2024 Appendix C K_init = 63 and norm_plasma_K = 0.0042; inside the Section 4.3.2 range [10000; 25000]

    # Appendix C prints MR_value_init = 1 and Section 4.3.2 gives the
    # global-SA range as "factors for mineral corticoid receptor
    # abundance [0.8; 1.2]" -- a multiplicative factor centred on 1.
    lmr <- log(1); label("Mineralocorticoid-receptor abundance factor MR (unitless; 1 = normal)") # Meid 2024 Appendix C MR_value_init = 1; Section 4.3.2 range [0.8; 1.2] confirms the factor is centred on 1

    # Figure A4 caption: "The nominal (starting) value of the latter is
    # set to zero in case of no hyperaldosteronism (A) or set to 0.1 in
    # case of little hyperaldosteronism (B)." Zero is also the only value
    # consistent with the printed Appendix C initial conditions: at
    # hyperaldo_effect = 0.1 the aldosterone term rises to 0.59 and
    # Aldo_effect_on_K_secretion jumps to 11.35, which is nowhere near
    # the steady state the printed state vector encodes. Set to 0
    # (FIXED) so the packaged typical individual is the normal-aldosterone
    # reference; simulate a hyperaldosteronism cohort by raising it (0.1
    # per Figure A4B, or anywhere in the Section 4.3.2 range [-0.5; 0.5]).
    # NOTE: the source applies IIV multiplicatively to this parameter
    # (hyperaldo_effect = Thyperaldo_effect * exp(eta)), so at a typical
    # value of exactly 0 the 36% IIV below has no effect. That is a
    # faithful transcription of the published parameterisation, not an
    # encoding slip.
    hyperaldo <- fixed(0); label("Hyperaldosteronism effect on aldosterone secretion (mEq-equivalent aldosterone offset; 0 = none)") # Meid 2024 Figure A4 caption (nominal 0; 0.1 = little hyperaldosteronism); Section 4.3.2 global-SA range [-0.5; 0.5]

    # ===================================================================
    # PHYSICAL CONSTANTS (Meid 2024 Appendix C)
    # Renamed from the source's single letters F / R / T, which collide
    # with rxode2 / R reserved names; values are unchanged.
    # ===================================================================
    faraday <- fixed(97); label("Faraday-type constant used in the Goldman-Hodgkin-Katz normalisation (source name F)")   # Meid 2024 Appendix C `F = 97`
    rgas    <- fixed(8.3145); label("Universal gas constant (J/(mol*K); source name R)")                                  # Meid 2024 Appendix C `R = 8.3145`
    tkelvin <- fixed(310.6); label("Body temperature (K; source name T)")                                                 # Meid 2024 Appendix C `T = 310.6`

    # ===================================================================
    # POTASSIUM SYSTEM REFERENCE VALUES (Meid 2024 Appendix C)
    # ===================================================================
    q_k_intracellular        <- fixed(465.87); label("Intracellular/extracellular potassium exchange rate constant (mL/min)")        # Meid 2024 Appendix C `Q_K_intracellular = 465.87`
    kinfusion                <- fixed(0); label("Intravenous potassium infusion rate (mEq/min; 0 = none, as in the source)")         # Meid 2024 Appendix C `Kinfusion = 0`
    norm_na_intake           <- fixed(0.0694444444444444); label("Reference sodium intake (mEq/min; 100 mEq/day)")                   # Meid 2024 Appendix C `norm_Na_intake = 0.0694444444444444`
    norm_aldo                <- fixed(0.49); label("Reference plasma aldosterone concentration (units not stated in the source; the value is 0.49)")                            # Meid 2024 Appendix C `norm_Aldo = 0.49`
    norm_plasma_k            <- fixed(0.0042); label("Reference plasma potassium concentration (mEq/mL; 4.2 mmol/L)")                # Meid 2024 Appendix C `norm_plasma_K = 0.0042`; Section 2 states this was FIXED to the population-typical value to reduce model complexity
    nom_intracellular_k_conc <- fixed(150); label("Reference intracellular potassium concentration (mEq/L)")                         # Meid 2024 Appendix C `nom_intracellular_K_conc = 150`
    v_ic                     <- fixed(25000); label("Intracellular fluid volume (mL)")                                               # Meid 2024 Appendix C `V_ic = 25000`
    principal_cell_intracellular_na <- fixed(0.008); label("Principal-cell intracellular sodium concentration (mEq/mL)")             # Meid 2024 Appendix C `principal_cell_intracellular_Na = 0.008`

    # Derivation: Appendix C prints no nephron count, but three printed
    # constants each pin it to 2,000,000 when Nain = 0.0694444444444444
    # and plasma sodium = 140 mmol/L:
    #   Na_out_loH0              = N * Nain / 0.05 / N       -> 1.38888888888889
    #   norm_Na_reabsorption_MCD = N * (Na_out_LoH/2 - Nain/N) -> 0.625
    #   potassium_in_MCD0        = N * CCD_K_out             -> 0.0945866666666666
    # and a fourth cross-check, eta_MCD0 = K_reabsorption_MCD_rate0 /
    # CCD_K_out = 0.154214829433323, also reproduces exactly. 2e6 is the
    # classical adult total nephron count (about 1e6 per kidney).
    nephrons_ref <- fixed(2000000); label("Reference total nephron count at gfr_ref (nephrons)")                                     # derived from Meid 2024 Appendix C Na_out_loH0 / norm_Na_reabsorption_MCD / potassium_in_MCD0 / eta_MCD0 (all reproduce exactly at 2e6)

    # Meid 2024 Section 2 mandates that nephron number be scaled
    # proportionally with the GFR covariate but never states the GFR that
    # corresponds to the reference nephron count. Set to the rounded
    # standard normal adult GFR; with nephrons_ref = 2e6 this gives a
    # single-nephron GFR of 62.5 nL/min, the textbook value. The model is
    # insensitive to gfr_ref at fixed single-nephron GFR (the paper's own
    # finding in Section 2), but gfr_ref does set how steeply nephron
    # count -- and therefore potassium excretion capacity -- falls with
    # declining renal function, so it is documented in the vignette.
    gfr_ref <- fixed(125); label("Reference glomerular filtration rate for the nephron-count scaling (mL/min)")                      # NOT reported in Meid 2024; rounded-standard normal adult GFR, chosen so SNGFR = gfr_ref / nephrons_ref = 62.5 nL/min. See vignette Errata.

    # ===================================================================
    # NEPHRON GEOMETRY (Meid 2024 Appendix C)
    # ===================================================================
    cnt_diameter <- fixed(0.0024); label("Connecting tubule diameter (cm)")                       # Meid 2024 Appendix C `CNT_diameter = 0.0024`
    cnt_length   <- fixed(0.4); label("Connecting tubule length (cm)")                            # Meid 2024 Appendix C `CNT_length = 0.4`
    dct_diameter <- fixed(0.0015); label("Distal convoluted tubule diameter (cm)")                # Meid 2024 Appendix C `DCT_diameter = 0.0015`
    dct_length   <- fixed(0.5); label("Distal convoluted tubule length (cm)")                     # Meid 2024 Appendix C `DCT_length = 0.5`
    ccd_diameter <- fixed(0.0025); label("Cortical collecting duct diameter (cm)")                # Meid 2024 Appendix C `CCD_diameter = 0.0025`
    ccd_length   <- fixed(0.2); label("Cortical collecting duct length (cm)")                     # Meid 2024 Appendix C `CCD_length = 0.2`
    sv_cnt       <- fixed(0.006); label("Connecting tubule cell surface-to-volume ratio (1/cm)")  # Meid 2024 Appendix C `SV_CNT = 0.006`
    sv_dct       <- fixed(0.0075); label("Distal convoluted tubule cell surface-to-volume ratio (1/cm)") # Meid 2024 Appendix C `SV_DCT = 0.0075`
    sv_ccd       <- fixed(0.004); label("Cortical collecting duct cell surface-to-volume ratio (1/cm)")  # Meid 2024 Appendix C `SV_CCD = 0.004`
    principal_fraction_cnt <- fixed(0.6); label("Principal-cell fraction of the connecting tubule (unitless)")        # Meid 2024 Appendix C `principal_fraction_CNT = 0.6`
    principal_fraction_ccd <- fixed(0.75); label("Principal-cell fraction of the cortical collecting duct (unitless)") # Meid 2024 Appendix C `principal_fraction_CCD = 0.75`
    dct_volume   <- fixed(8.83572933822129e-07); label("Single-nephron distal convoluted tubule luminal volume (mL)")  # Meid 2024 Appendix C `DCT_volume = 8.83572933822129e-07`
    cnt_volume   <- fixed(1.80955736846772e-06); label("Single-nephron connecting tubule luminal volume (mL)")         # Meid 2024 Appendix C `CNT_volume = 1.80955736846772e-06`
    ccd_volume   <- fixed(9.8174770424681e-07); label("Single-nephron cortical collecting duct luminal volume (mL)")   # Meid 2024 Appendix C `CCD_volume = 9.8174770424681e-07`
    cnt_water_reabs_fraction <- fixed(0.7); label("Fractional water reabsorption in the connecting tubule (unitless)")             # Meid 2024 Appendix C `CNT_water_reabs_fraction = 0.7`
    ccd_water_reabs_fraction <- fixed(0.75); label("Fractional water reabsorption in the cortical collecting duct (unitless)")     # Meid 2024 Appendix C `CCD_water_reabs_fraction = 0.75`

    # ===================================================================
    # TRANSPORT AND FEEDBACK PARAMETERS (Meid 2024 Appendix C)
    # ===================================================================
    norm_na_reabsorption_mcd <- fixed(0.624999999999999); label("Reference whole-body medullary collecting duct sodium reabsorption (mEq/min)") # Meid 2024 Appendix C `norm_Na_reabsorption_MCD = 0.624999999999999`
    k_reabsorption_mcd_rate0 <- fixed(7.2933333333333e-09); label("Single-nephron baseline medullary collecting duct potassium reabsorption rate (mEq/min)") # Meid 2024 Appendix C `K_reabsorption_MCD_rate0 = 7.2933333333333e-09`
    baseline_k_luminal_permeability <- fixed(2.4935e-05); label("Baseline luminal potassium permeability (cm/min)")             # Meid 2024 Appendix C `baseline_K_luminal_permeability = 2.4935e-05`
    k_basolateral_permeability      <- fixed(3.43e-05); label("Basolateral potassium permeability (cm/min)")                    # Meid 2024 Appendix C `K_basolateral_permeability = 3.43e-05`
    j_na_active_max                 <- fixed(0.0001466); label("Maximum Na/K-ATPase active sodium flux (mEq/(cm^2*min))")       # Meid 2024 Appendix C `J_Na_active_max = 0.0001466`
    luminal_potential_difference    <- fixed(-18.4); label("Luminal transmembrane potential difference (mV)")                   # Meid 2024 Appendix C `luminal_potential_difference = -18.4`
    basolateral_potential_difference <- fixed(-78.2); label("Basolateral transmembrane potential difference (mV)")              # Meid 2024 Appendix C `basolateral_potential_difference = -78.2`
    m_k_aldo         <- fixed(951.2); label("Plasma potassium effect on plasma aldosterone (mL/mEq; exponential slope)")        # Meid 2024 Appendix C `m_K_ALDO = 951.2`
    m_na_aldo        <- fixed(15.569); label("Sodium intake effect on plasma aldosterone (min/mEq; exponential slope)")         # Meid 2024 Appendix C `m_Na_ALDO = 15.569`
    aldo_ksec_scale  <- fixed(103.5); label("Aldosterone effect on luminal potassium permeability (reciprocal of the aldosterone unit; linear scale)")      # Meid 2024 Appendix C `Aldo_KSec_scale = 103.5`
    m_plasmak_mcd    <- fixed(0.000000883); label("Plasma potassium effect on medullary collecting duct potassium reabsorption (unitless exponential slope)") # Meid 2024 Appendix C `m_plasmaK_MCD = 0.000000883`
    m_na_mcd         <- fixed(0.69775); label("Sodium-coupled potassium secretion factor in the medullary collecting duct (unitless)") # Meid 2024 Appendix C `m_Na_MCD = 0.69775`

    # ===================================================================
    # SPIRONOLACTONE / CANRENONE PK AND MRA PHARMACODYNAMICS
    # (Meid 2024 Appendix C; traced by Meid 2024 to Gardiner 1989.)
    # All values are FIXED -- the source estimated none of them and
    # carried them unchanged from the upstream model.
    # ===================================================================
    emax_spiro   <- fixed(0.9978); label("Maximum fractional mineralocorticoid-receptor inhibition by canrenone (unitless)")  # Meid 2024 Appendix C `E_MAX_spiro = 0.9978`
    ec50_spiro   <- fixed(1.8296); label("Canrenone concentration at half-maximal receptor inhibition (ng/mL; = 5.4 nmol/L at MW 340.5)")                # Meid 2024 Appendix C `EC50_spiro = 1.8296`
    koff_mra     <- fixed(3.4035); label("Mineralocorticoid-receptor antagonist off-rate (1/min; kon is set equal to koff)")  # Meid 2024 Appendix C `Koff_MRA = 3.4035` and `Kon_MRA = Koff_MRA`
    ka_spiro     <- fixed(0.01524458); label("Spironolactone transit/absorption rate constant (1/min; shared by all three transit steps)") # Meid 2024 Appendix C `Ka_spiro = 0.01524458`
    v1_spiro     <- fixed(7.15696); label("Spironolactone central volume of distribution (L)")                                # Meid 2024 Appendix C `V1_spiro = 7.15696`
    cl_spiro     <- fixed(8.07626); label("Spironolactone total clearance (L/min)")                                          # Meid 2024 Appendix C `CL_spiro = 8.07626`
    lcl_canrenone <- fixed(log(0.222487)); label("Canrenone clearance (L/min); log-transformed per the `cl`-parameter naming convention") # Meid 2024 Appendix C `CL_canrenone = 0.222487`
    v_canrenone  <- fixed(70.47); label("Canrenone central volume of distribution (L)")                                       # Meid 2024 Appendix C `V_canrenone = 70.47`
    v2_canrenone <- fixed(8.021); label("Canrenone peripheral volume of distribution (L)")                                    # Meid 2024 Appendix C `V2_canrenone = 8.021`
    lq_canrenone <- fixed(log(0.110275)); label("Canrenone inter-compartmental clearance (L/min); log-transformed per the `q`-parameter naming convention") # Meid 2024 Appendix C `Q_canrenone = 0.110275`
    spiro_fmetabolized   <- fixed(0.19311); label("Fraction of spironolactone clearance forming canrenone (unitless)")        # Meid 2024 Appendix C `Spiro_Fmetabolized = 0.19311`
    # `Spiro_bioavailability = 0.91097` is declared in Meid 2024 Appendix C
    # but never applied to the depot compartment anywhere in the printed
    # listing, so it is recorded here rather than carried as a live
    # parameter. See the bioavailability note at the end of model().

    # ===================================================================
    # INTER-INDIVIDUAL VARIABILITY (Meid 2024 Section 2, fourth workflow
    # step). Reported as %CV in historical data from 20 separate
    # patients; converted with the log-normal identity
    # omega^2 = log(1 + CV^2). The source caps the two 131% values at 80%
    # "in the subsequent application", so the applied 80% is used here
    # and the uncapped 131% estimate is recorded in the comment.
    # ===================================================================
    etalkin  ~ 0.49469624    # Meid 2024 Section 2: Kin CV = 131%, capped to 80% in the application -> log(1 + 0.80^2) = 0.49469624
    etalnain ~ 0.00270035    # Meid 2024 Section 2: Nain CV = 5.2% -> log(1 + 0.052^2) = 0.00270035
    etalv_ecf ~ 0.03188566   # Meid 2024 Section 2: V_ecf CV = 18% -> log(1 + 0.18^2) = 0.03188566
    etalmr   ~ 0.49469624    # Meid 2024 Section 2: MR CV = 131%, capped to 80% in the application -> log(1 + 0.80^2) = 0.49469624
    etahyperaldo ~ 0.12183267 # Meid 2024 Section 2: hyperaldo_effect CV = 36% -> log(1 + 0.36^2) = 0.12183267; inert while the typical value is 0 (see the hyperaldo note above)

    # ===================================================================
    # RESIDUAL ERROR
    # The source is a Bayesian prediction application and never prints a
    # fitted residual-error model. Figure A6 states the magnitude the
    # authors themselves used for simulated potassium measurements:
    # "the five simulated potassium values (multiplied with mean 1 and
    # SD of 0.075)" -- a proportional error of 7.5% CV. Figure A7 repeats
    # the experiment at SD 0.1 as a larger-noise scenario.
    # ===================================================================
    propSd <- fixed(0.075); label("Proportional residual error on plasma potassium (fraction)") # Meid 2024 Figure A6 caption (simulated potassium values multiplied by mean 1, SD 0.075); Figure A7 uses 0.1 for the higher-noise scenario
  })

  model({
    # -------------------------------------------------------------------
    # 1. Individual parameters. The source writes these as
    #    `Kin = exp(log(TKin) + eta.Kin)` etc. (log-normal IIV on four
    #    parameters) and `hyperaldo_effect = Thyperaldo_effect *
    #    exp(eta.hyperaldo_effect)` (multiplicative IIV on the fifth).
    # -------------------------------------------------------------------
    kin   <- exp(lkin + etalkin)
    nain  <- exp(lnain + etalnain)
    v_ecf <- exp(lv_ecf + etalv_ecf)
    mr    <- exp(lmr + etalmr)
    hyperaldo_effect <- hyperaldo * exp(etahyperaldo)

    # -------------------------------------------------------------------
    # 2. Covariates. `SOD` is plasma sodium (mmol/L) and `CRCL` is
    #    CKD-EPI eGFR. Nephron count scales proportionally with renal
    #    function so that single-nephron GFR stays at its reference value
    #    (Meid 2024 Section 2).
    # -------------------------------------------------------------------
    norm_plasma_na     <- SOD
    gfr                <- CRCL
    number_of_nephrons <- nephrons_ref * (CRCL / gfr_ref)

    # -------------------------------------------------------------------
    # 3. Spironolactone -> canrenone -> mineralocorticoid-receptor
    #    occupancy. `mra_effect` is the fraction of receptor NOT blocked;
    #    it multiplies the aldosterone-driven potassium secretion below.
    #
    #    DOSE UNIT. Meid 2024 Appendix C never states the amount unit of
    #    the depot compartment, but `EC50_spiro = 1.8296` must be on the
    #    same scale as the `canrenone` state, and that pins it. Entering
    #    the dose in MICROGRAMS (25 mg -> amt = 25000) makes `central`
    #    and `central_canrenone` come out in ng/mL, which gives:
    #      * canrenone Cmax 36 / 72 / 144 ng/mL for 25 / 50 / 100 mg,
    #        matching published canrenone exposure;
    #      * EC50 = 1.8296 ng/mL = 5.4 nmol/L (canrenone MW 340.5),
    #        matching canrenone's mineralocorticoid-receptor affinity;
    #      * a plasma-potassium rise of +0.21 / +0.31 / +0.44 mmol/L at
    #        25 / 50 / 100 mg once daily, the clinically expected effect.
    #    Entering the dose in milligrams instead puts canrenone three
    #    orders of magnitude below EC50, so the receptor is never
    #    meaningfully occupied and spironolactone changes plasma
    #    potassium by less than 0.003 mmol/L at any dose -- which would
    #    contradict the source's own Discussion, where the model is
    #    described as capable of predicting "a clinically too high
    #    increase in potassium" in response to spironolactone.
    # -------------------------------------------------------------------
    cl_canrenone <- exp(lcl_canrenone)
    q_canrenone  <- exp(lq_canrenone)
    e_mra_spiro  <- emax_spiro * central_canrenone / (central_canrenone + ec50_spiro)
    kon_mra      <- koff_mra

    # -------------------------------------------------------------------
    # 4. Potassium pools and intracellular shift.
    # -------------------------------------------------------------------
    plasma_k              <- k_ecf / v_ecf
    intracellular_k_conc  <- intracellular_k / v_ic
    intracellular_potassium_flux <-
      q_k_intracellular *
      (plasma_k - norm_plasma_k - (intracellular_k_conc - nom_intracellular_k_conc))

    # -------------------------------------------------------------------
    # 5. Aldosterone feedback. Aldosterone rises exponentially with
    #    plasma potassium, is suppressed by sodium intake above the
    #    reference, and is offset by the hyperaldosteronism term. The
    #    resulting secretion driver is gated by receptor availability
    #    (mra_effect) and receptor abundance (mr).
    # -------------------------------------------------------------------
    aldo <- max(
      0,
      norm_aldo * exp(m_k_aldo * (plasma_k - norm_plasma_k)) +
        max(0, exp(-m_na_aldo * (nain - norm_na_intake)) - 1) +
        hyperaldo_effect
    )
    aldo_effect_on_k_secretion <- mra_effect * mr * max(0, 1 + aldo_ksec_scale * (aldo - norm_aldo))

    # -------------------------------------------------------------------
    # 6. Glomerulus and proximal tubule / loop of Henle. Potassium is
    #    reabsorbed in the same fraction as sodium.
    # -------------------------------------------------------------------
    sngfr                   <- gfr / number_of_nephrons
    filtered_k              <- max(0, sngfr * plasma_k)
    filtered_na             <- sngfr * (norm_plasma_na / 1000)
    max_pt_loh_reabsorption <- filtered_na * 0.96
    na_out_loh              <- max(nain / (1 - 0.95) / number_of_nephrons, filtered_na - max_pt_loh_reabsorption)
    eta_pt_loh              <- (filtered_na - na_out_loh) / filtered_na
    k_reabsorption_pt_loh   <- eta_pt_loh * filtered_k
    pt_loh_k_out            <- filtered_k - k_reabsorption_pt_loh

    # -------------------------------------------------------------------
    # 7. Luminal potassium concentrations in the three cortical segments.
    # -------------------------------------------------------------------
    dct_luminal_k_conc <- max(0, dct_luminal_k_amount) / dct_volume
    cnt_luminal_k_conc <- max(0, cnt_luminal_k_amount) / cnt_volume
    ccd_luminal_k_conc <- max(0, ccd_luminal_k_amount) / ccd_volume

    # -------------------------------------------------------------------
    # 8. Goldman-Hodgkin-Katz passive fluxes and Na/K-ATPase active
    #    basolateral uptake, per segment. `F`, `R` and `T` in the source
    #    listing are the constants renamed faraday / rgas / tkelvin here.
    # -------------------------------------------------------------------
    normalized_luminal_potential_difference     <- faraday * luminal_potential_difference / (rgas * tkelvin)
    normalized_basolateral_potential_difference <- faraday * basolateral_potential_difference / (rgas * tkelvin)

    dct_k_in     <- pt_loh_k_out
    dct_water_in <- sngfr * (1 - eta_pt_loh)

    dct_k_passive_flux_lumenal <-
      baseline_k_luminal_permeability * normalized_luminal_potential_difference *
      (-dct_luminal_k_conc + dct_cell_k_conc * exp(-normalized_luminal_potential_difference)) /
      (1 - exp(-normalized_luminal_potential_difference))
    dct_k_passive_flux_basolateral <-
      -k_basolateral_permeability * normalized_basolateral_potential_difference *
      (-plasma_k + dct_cell_k_conc * exp(-normalized_basolateral_potential_difference)) /
      (1 - exp(-normalized_basolateral_potential_difference))

    k_k       <- (0.1 * (1 + 140 / 18.5)) / 1000
    k_na_dct  <- (0.2 * (1 + dct_cell_k_conc * 1000 / 8.33)) / 1000
    j_na_active_max_eff <- j_na_active_max * max(0, aldo_effect_on_k_secretion)
    dct_k_active_flux_basolateral <-
      (2 / 3) * j_na_active_max_eff *
      ((principal_cell_intracellular_na / (principal_cell_intracellular_na + k_na_dct))^3) *
      ((plasma_k / (plasma_k + k_k))^2)

    cnt_k_passive_flux_lumenal <-
      baseline_k_luminal_permeability * normalized_luminal_potential_difference *
      (-cnt_luminal_k_conc + cnt_cell_k_conc * exp(-normalized_luminal_potential_difference)) /
      (1 - exp(-normalized_luminal_potential_difference))
    cnt_k_passive_flux_basolateral <-
      -k_basolateral_permeability * normalized_basolateral_potential_difference *
      (-plasma_k + cnt_cell_k_conc * exp(-normalized_basolateral_potential_difference)) /
      (1 - exp(-normalized_basolateral_potential_difference))
    k_na_cnt <- (0.2 * (1 + cnt_cell_k_conc * 1000 / 8.33)) / 1000
    cnt_k_active_flux_basolateral <-
      (2 / 3) * j_na_active_max_eff *
      ((principal_cell_intracellular_na / (principal_cell_intracellular_na + k_na_cnt))^3) *
      ((plasma_k / (plasma_k + k_k))^2)

    ccd_k_passive_flux_lumenal <-
      baseline_k_luminal_permeability * normalized_luminal_potential_difference *
      (-ccd_luminal_k_conc + ccd_cell_k_conc * exp(-normalized_luminal_potential_difference)) /
      (1 - exp(-normalized_luminal_potential_difference))
    ccd_k_passive_flux_basolateral <-
      -k_basolateral_permeability * normalized_basolateral_potential_difference *
      (-plasma_k + ccd_cell_k_conc * exp(-normalized_basolateral_potential_difference)) /
      (1 - exp(-normalized_basolateral_potential_difference))
    k_na_ccd <- (0.2 * (1 + ccd_cell_k_conc * 1000 / 8.33)) / 1000
    ccd_k_active_flux_basolateral <-
      (2 / 3) * j_na_active_max_eff *
      ((principal_cell_intracellular_na / (principal_cell_intracellular_na + k_na_ccd))^3) *
      ((plasma_k / (plasma_k + k_k))^2)

    # -------------------------------------------------------------------
    # 9. Segment surface areas, secretion rates, and water / potassium
    #    flow between segments.
    # -------------------------------------------------------------------
    dct_sa <- pi * dct_diameter * dct_length
    cnt_sa <- pi * cnt_diameter * cnt_length
    ccd_sa <- pi * ccd_diameter * ccd_length

    dct_k_secretion_rate <- dct_k_passive_flux_lumenal * dct_sa
    cnt_k_secretion_rate <- cnt_k_passive_flux_lumenal * cnt_sa * principal_fraction_cnt
    ccd_k_secretion_rate <- ccd_k_passive_flux_lumenal * ccd_sa * principal_fraction_ccd

    dct_water_out <- dct_water_in
    cnt_water_out <- dct_water_out * (1 - cnt_water_reabs_fraction)
    ccd_water_out <- dct_water_in * (1 - ccd_water_reabs_fraction)

    dct_k_out <- dct_luminal_k_conc * dct_water_in
    cnt_k_out <- cnt_luminal_k_conc * dct_water_in
    ccd_k_out <- ccd_luminal_k_conc * ccd_water_out

    # -------------------------------------------------------------------
    # 10. Medullary collecting duct: baseline reabsorption modulated by
    #     plasma potassium and by sodium delivery. `cd_k_out` is the
    #     whole-body urinary potassium excretion rate.
    # -------------------------------------------------------------------
    k_mcd_effect       <- max(0, exp(m_plasmak_mcd * ((norm_plasma_k - plasma_k) / norm_plasma_k)) - 1)
    na_reabsorption_mcd <- max(0, na_out_loh * 0.5 - nain / number_of_nephrons)
    k_secretion_mcd    <- max(0, (norm_na_reabsorption_mcd / number_of_nephrons - na_reabsorption_mcd) * m_na_mcd)
    k_reabsorption_mcd_rate <- k_reabsorption_mcd_rate0 + k_mcd_effect + k_secretion_mcd
    k_reabsorption_cd  <- min(k_reabsorption_mcd_rate, ccd_k_out)
    eta_mcd            <- max(0, k_reabsorption_mcd_rate / ccd_k_out)
    cd_k_out           <- number_of_nephrons * (ccd_k_out - k_reabsorption_cd)

    # -------------------------------------------------------------------
    # 11. ODE system, in the order the source declares it.
    #     Spironolactone absorption is a three-step chain sharing one
    #     rate constant; `central` and the two canrenone states are
    #     CONCENTRATION states (the source divides each derivative by the
    #     corresponding volume).
    # -------------------------------------------------------------------
    d/dt(depot)    <- -ka_spiro * depot
    d/dt(transit1) <-  ka_spiro * depot - ka_spiro * transit1
    d/dt(transit2) <-  ka_spiro * transit1 - ka_spiro * transit2

    d/dt(k_ecf)           <- kin + kinfusion - cd_k_out - intracellular_potassium_flux
    d/dt(intracellular_k) <- intracellular_potassium_flux

    d/dt(dct_luminal_k_amount) <- dct_k_in  + dct_k_secretion_rate - dct_k_out
    d/dt(cnt_luminal_k_amount) <- dct_k_out + cnt_k_secretion_rate - cnt_k_out
    d/dt(ccd_luminal_k_amount) <- cnt_k_out + ccd_k_secretion_rate - ccd_k_out

    d/dt(dct_cell_k_conc) <- (-dct_k_passive_flux_lumenal + dct_k_passive_flux_basolateral + dct_k_active_flux_basolateral) / sv_dct
    d/dt(cnt_cell_k_conc) <- (-cnt_k_passive_flux_lumenal + cnt_k_passive_flux_basolateral + cnt_k_active_flux_basolateral) / sv_cnt
    d/dt(ccd_cell_k_conc) <- (-ccd_k_passive_flux_lumenal + ccd_k_passive_flux_basolateral + ccd_k_active_flux_basolateral) / sv_ccd

    d/dt(potassium_excretion) <- cd_k_out
    d/dt(mra_effect)          <- kon_mra * (1 - e_mra_spiro) - koff_mra * mra_effect

    d/dt(central) <-
      (ka_spiro * transit2 -
         cl_spiro * (1 - spiro_fmetabolized) * central -
         cl_spiro * spiro_fmetabolized * central) / v1_spiro
    d/dt(central_canrenone) <-
      (cl_spiro * spiro_fmetabolized * central -
         cl_canrenone * central_canrenone -
         (q_canrenone * central_canrenone - q_canrenone * peripheral1_canrenone)) / v_canrenone
    d/dt(peripheral1_canrenone) <-
      (q_canrenone * central_canrenone - q_canrenone * peripheral1_canrenone) / v2_canrenone

    # -------------------------------------------------------------------
    # 12. Initial conditions -- the physiological equilibrium that the
    #     source passes in as function-argument defaults. Two are written
    #     as the product that generates the printed number rather than as
    #     the number itself, so that the baseline stays physiological
    #     when the corresponding volume carries IIV:
    #       k_ecf(0)           = 15000 * 0.0042 = 63       (printed K_init)
    #       intracellular_k(0) = 25000 * 150    = 3750000  (printed init)
    #     All spironolactone / canrenone states start at zero.
    # -------------------------------------------------------------------
    k_ecf(0)                <- v_ecf * norm_plasma_k
    intracellular_k(0)      <- v_ic * nom_intracellular_k_conc
    dct_luminal_k_amount(0) <- 0.000000005836802
    cnt_luminal_k_amount(0) <- 0.00000001527986
    ccd_luminal_k_amount(0) <- 0.00000003744125
    dct_cell_k_conc(0)      <- 0.150002
    cnt_cell_k_conc(0)      <- 0.1501202
    ccd_cell_k_conc(0)      <- 0.1520374
    potassium_excretion(0)  <- 3225.623
    mra_effect(0)           <- 1

    # -------------------------------------------------------------------
    # 13. Outputs. The source sets `Cc = plasma_K`, so the canonical
    #     observation variable of this model is plasma potassium in
    #     mEq/mL, NOT a drug concentration. `plasmaK` restates it on the
    #     clinical mmol/L scale and `Cc_canrenone` exposes the active
    #     metabolite that drives the receptor effect.
    #
    #     NOTE ON BIOAVAILABILITY: the source declares
    #     `Spiro_bioavailability = 0.91097` but never applies it to the
    #     depot compartment in the printed listing. This file reproduces
    #     that behaviour -- there is deliberately no `f(depot)` line. To
    #     apply it, add `f(depot) <- 0.91097` here.
    # -------------------------------------------------------------------
    plasmaK      <- plasma_k * 1000
    Cc_canrenone <- central_canrenone
    Cc <- plasma_k
    Cc ~ prop(propSd)
  })
}
