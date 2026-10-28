Augustin_2023_warfarin_qsp <- function() {
  description <- paste(
    "QSP. Warfarin clinical-trial (CT) model of Augustin 2023: the Wajima 2009 humoral",
    "coagulation-network model (51 coagulation / vitamin K states) driven by a one-compartment",
    "oral warfarin PK model, with the population model the authors assembled from Hartmann",
    "2016 / 2020 and Hamberg 2010. Warfarin inhibits the two VKORC-catalysed vitamin K",
    "conversions (Hill Imax), which lowers vitamin K hydroquinone and with it the synthesis of",
    "factors II, VII, IX, X and proteins C and S. Covariates: CYP2C9 *2 / *3 allele counts and",
    "age on the warfarin elimination rate, VKORC1 -1639A allele count on the warfarin EC50 and",
    "on both VKORC-catalysed vitamin K conversion rates. IIV: lognormal on warfarin volume,",
    "elimination rate, EC50 and the two conversion rates; normal (20% CV) on ten coagulation-",
    "factor production rates. Inter-occasion variability: the vitamin K input rate is",
    "redrawn every treatment day (OCC = treatment day 1-56, SD 10%). Simulations start from",
    "the Wajima initial concentrations and need a 100-day (2400 h) drug-free run-in before",
    "the first dose so each individual reaches its own steady state. The INR is not an",
    "algebraic output: it is read from a simulated prothrombin-time test on a blood sample",
    "(sample diluted 1:2, 300 nM tissue factor added, PT = time for the fibrin AUC state",
    "fibrin_auc to reach 1500 nM*s, INR = PT / 11.8 s), which the validation vignette",
    "implements as a second solve of this same model. Solve with method = \"lsoda\": under",
    "the default liblsoda, round-off near zero can set off a spurious coagulation burst.",
    "No parameter was estimated by Augustin 2023; values were adopted from the cited",
    "models or hand-calibrated."
  )
  reference <- paste(
    "Augustin D, Lambert B, Robinson M, Wang K, Gavaghan D.",
    "Simulating clinical trials for model-informed precision dosing: using warfarin",
    "treatment as a use case. Front Pharmacol. 2023;14:1270443. doi:10.3389/fphar.2023.1270443.",
    "Model structure: Section 2.2 Eq 6-9 and Supplementary Appendix S1 (Tables S1.1-S1.3,",
    "Eq S1.1-S1.9). Parameter values: the authors' code and data deposit named in Section 2.2,",
    "https://github.com/DavAug/mipd-warfarin (commit 4710c062), models/wajima_coagulation_model.xml",
    "(SBML) and models/hartmann_coagulation_model_parameters.csv. Underlying models:",
    "Wajima T, Isbister GK, Duffull SB. Clin Pharmacol Ther. 2009;86(3):290-298",
    "(doi:10.1038/clpt.2009.87); Hartmann S et al. CPT Pharmacometrics Syst Pharmacol.",
    "2016;5(10):554-564 (doi:10.1002/psp4.12111) and Front Pharmacol. 2020;11:1041",
    "(doi:10.3389/fphar.2020.01041); Hamberg AK et al. Clin Pharmacol Ther.",
    "2010;87(6):727-734 (doi:10.1038/clpt.2010.37)."
  )
  vignette <- "Augustin_2023_warfarin_qsp"

  # Wajima 2009 coagulation-network species keep their SBML / Table S1.2 names;
  # fibrin_auc is the auxiliary fibrin-integral state of Section S1.1.2.
  paper_specific_compartments <- c(
    "coagulation_factor_ii",
    "activated_coagulation_factor_ii",
    "coagulation_factor_v",
    "activated_coagulation_factor_v",
    "coagulation_factor_vii",
    "activated_coagulation_factor_vii",
    "coagulation_factor_viii",
    "activated_coagulation_factor_viii",
    "coagulation_factor_ix",
    "activated_coagulation_factor_ix",
    "coagulation_factor_x",
    "activated_coagulation_factor_x",
    "coagulation_factor_xi",
    "activated_coagulation_factor_xi",
    "coagulation_factor_xii",
    "activated_coagulation_factor_xii",
    "coagulation_factor_xiii",
    "activated_coagulation_factor_xiii",
    "contact_system_activator",
    "fibrin",
    "cross_linked_fibrin",
    "fibrin_degradation_product",
    "d_dimer",
    "fibrinogen",
    "kallikrein",
    "prekallikrein",
    "plasmin",
    "plasminogen",
    "protein_c",
    "activated_protein_c",
    "protein_s",
    "thrombomodulin",
    "tissue_factor",
    "tissue_factor_pathway_inhibitor",
    "vitamin_k",
    "vitamin_k_epoxide",
    "vitamin_k_hydroquinone",
    "vitamin_k_peripheral",
    "aii_atiii_heparin_complex",
    "aii_at_complex",
    "aii_tmod_complex",
    "av_ax_complex",
    "vii_tf_complex",
    "avii_tf_complex",
    "avii_tf_ax_tfpi_complex",
    "aix_aviii_complex",
    "aix_atiii_heparin_complex",
    "ax_atiii_heparin_complex_complex",
    "ax_tfpi_complex",
    "apc_ps_complex",
    "atiii_heparin_complex",
    "fibrin_auc"
  )

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "warfarin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "warfarin", units = "mg", specimen = "plasma", verified = TRUE),
    coagulation_factor_ii = list(analyte = "coagulation factor ii", units = "nM", specimen = "plasma", verified = TRUE),
    activated_coagulation_factor_ii = list(
      analyte = "activated coagulation factor ii",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    coagulation_factor_v = list(analyte = "coagulation factor v", units = "nM", specimen = "plasma", verified = TRUE),
    activated_coagulation_factor_v = list(
      analyte = "activated coagulation factor v",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    coagulation_factor_vii = list(
      analyte = "coagulation factor vii",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    activated_coagulation_factor_vii = list(
      analyte = "activated coagulation factor vii",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    coagulation_factor_viii = list(
      analyte = "coagulation factor viii",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    activated_coagulation_factor_viii = list(
      analyte = "activated coagulation factor viii",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    coagulation_factor_ix = list(analyte = "coagulation factor ix", units = "nM", specimen = "plasma", verified = TRUE),
    activated_coagulation_factor_ix = list(
      analyte = "activated coagulation factor ix",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    coagulation_factor_x = list(analyte = "coagulation factor x", units = "nM", specimen = "plasma", verified = TRUE),
    activated_coagulation_factor_x = list(
      analyte = "activated coagulation factor x",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    coagulation_factor_xi = list(analyte = "coagulation factor xi", units = "nM", specimen = "plasma", verified = TRUE),
    activated_coagulation_factor_xi = list(
      analyte = "activated coagulation factor xi",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    coagulation_factor_xii = list(
      analyte = "coagulation factor xii",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    activated_coagulation_factor_xii = list(
      analyte = "activated coagulation factor xii",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    coagulation_factor_xiii = list(
      analyte = "coagulation factor xiii",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    activated_coagulation_factor_xiii = list(
      analyte = "activated coagulation factor xiii",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    contact_system_activator = list(
      analyte = "contact system activator",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    fibrin = list(analyte = "fibrin", units = "nM", specimen = "plasma", verified = TRUE),
    cross_linked_fibrin = list(analyte = "cross linked fibrin", units = "nM", specimen = "plasma", verified = TRUE),
    fibrin_degradation_product = list(
      analyte = "fibrin degradation product",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    d_dimer = list(analyte = "d dimer", units = "nM", specimen = "plasma", verified = TRUE),
    fibrinogen = list(analyte = "fibrinogen", units = "nM", specimen = "plasma", verified = TRUE),
    kallikrein = list(analyte = "kallikrein", units = "nM", specimen = "plasma", verified = TRUE),
    prekallikrein = list(analyte = "prekallikrein", units = "nM", specimen = "plasma", verified = TRUE),
    plasmin = list(analyte = "plasmin", units = "nM", specimen = "plasma", verified = TRUE),
    plasminogen = list(analyte = "plasminogen", units = "nM", specimen = "plasma", verified = TRUE),
    protein_c = list(analyte = "protein c", units = "nM", specimen = "plasma", verified = TRUE),
    activated_protein_c = list(analyte = "activated protein c", units = "nM", specimen = "plasma", verified = TRUE),
    protein_s = list(analyte = "protein s", units = "nM", specimen = "plasma", verified = TRUE),
    thrombomodulin = list(analyte = "thrombomodulin", units = "nM", specimen = "plasma", verified = TRUE),
    tissue_factor = list(analyte = "tissue factor", units = "nM", specimen = "plasma", verified = TRUE),
    tissue_factor_pathway_inhibitor = list(
      analyte = "tissue factor pathway inhibitor",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    vitamin_k = list(analyte = "vitamin k", units = "nM", specimen = "plasma", verified = TRUE),
    vitamin_k_epoxide = list(analyte = "vitamin k epoxide", units = "nM", specimen = "plasma", verified = TRUE),
    vitamin_k_hydroquinone = list(
      analyte = "vitamin k hydroquinone",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    vitamin_k_peripheral = list(
      analyte = "vitamin K (peripheral pool)",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    aii_atiii_heparin_complex = list(
      analyte = "aii atiii heparin complex",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    aii_at_complex = list(analyte = "aii at complex", units = "nM", specimen = "plasma", verified = TRUE),
    aii_tmod_complex = list(analyte = "aii tmod complex", units = "nM", specimen = "plasma", verified = TRUE),
    av_ax_complex = list(analyte = "av ax complex", units = "nM", specimen = "plasma", verified = TRUE),
    vii_tf_complex = list(analyte = "vii tf complex", units = "nM", specimen = "plasma", verified = TRUE),
    avii_tf_complex = list(analyte = "avii tf complex", units = "nM", specimen = "plasma", verified = TRUE),
    avii_tf_ax_tfpi_complex = list(
      analyte = "avii tf ax tfpi complex",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    aix_aviii_complex = list(analyte = "aix aviii complex", units = "nM", specimen = "plasma", verified = TRUE),
    aix_atiii_heparin_complex = list(
      analyte = "aix atiii heparin complex",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    ax_atiii_heparin_complex_complex = list(
      analyte = "ax atiii heparin complex complex",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    ax_tfpi_complex = list(analyte = "ax tfpi complex", units = "nM", specimen = "plasma", verified = TRUE),
    apc_ps_complex = list(analyte = "apc ps complex", units = "nM", specimen = "plasma", verified = TRUE),
    atiii_heparin_complex = list(analyte = "atiii heparin complex", units = "nM", specimen = "plasma", verified = TRUE),
    fibrin_auc = list(analyte = "fibrin (time integral)", units = "nM*h", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    CYP2C9_S2_COUNT = list(
      description = "Number of CYP2C9*2 alleles (0, 1 or 2)",
      units = "(count, 0/1/2 alleles per subject)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Each *2 allele lowers the typical warfarin elimination rate by",
        "e_cyp2c9s2_kel / 2 (additive per-allele contributions, Hamberg 2010 structure,",
        "Eq S1.5). The CYP2C9*1 count is implied (2 - *2 count - *3 count) and is not a",
        "model input."
      ),
      source_name = "CYP2C9 genotype (*1*1, *1*2, *1*3, *2*2, *2*3, *3*3)"
    ),
    CYP2C9_S3_COUNT = list(
      description = "Number of CYP2C9*3 alleles (0, 1 or 2)",
      units = "(count, 0/1/2 alleles per subject)",
      type = "continuous",
      reference_category = NULL,
      notes = "Each *3 allele lowers the typical warfarin elimination rate by e_cyp2c9s3_kel / 2 (Eq S1.5).",
      source_name = "CYP2C9 genotype (*1*1, *1*2, *1*3, *2*2, *2*3, *3*3)"
    ),
    VKORC1_1639G_COUNT = list(
      description = "Number of VKORC1 G alleles (2 = GG, 1 = GA, 0 = AA)",
      units = "(count, 0/1/2 alleles per subject)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The paper writes the VKORC1 genotype as GG / GA / AA (the -1639 G>A promoter",
        "polymorphism used by Hamberg 2010). The model uses the A-allele count",
        "2 - VKORC1_1639G_COUNT: each A allele lowers the typical warfarin EC50 by",
        "e_vkorc1aa_ec50 / 2 and both VKORC-catalysed vitamin K conversion rates by",
        "e_vkorc1aa_vk_vkh2 / 2 and e_vkorc1aa_vko_vk / 2 (Eq S1.3, S1.8)."
      ),
      source_name = "VKORC1 genotype (GG, GA, AA)"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred at 71 years: kel is multiplied by 1 - tanh(e_age_kel * (AGE - 71)) (Eq S1.5).",
      source_name = "age"
    ),
    OCC = list(
      description = "Treatment day index for the daily vitamin K input inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "OCC = k on treatment day k (hours 24 * (k - 1) to 24 * k after the first dose),",
        "k = 1-56; each day draws its own vitamin K input multiplier (Section 2.2 Eq 7,",
        "Algorithms S2-S4). Use OCC = 0 for the drug-free run-in and to switch the",
        "inter-occasion variability off."
      ),
      source_name = "simulation day"
    )
  )

  population <- list(
    species = "human (virtual population; the model was assembled, not fitted)",
    n_subjects = 1000L,
    n_studies = 0L,
    disease_state = "Patients starting warfarin anticoagulation (target INR 2-3), simulated",
    dose_range = paste(
      "Oral warfarin once daily; simulated trials used 10 / 7.5 / 5 mg induction",
      "doses followed by INR-guided adjustment, and MIPD-chosen doses of 0-30 mg"
    ),
    age_range = paste(
      "Simulated cohorts: phase I median 69 (5th-95th percentile 58-79) years,",
      "phase II 65 (56-76), phase III 50 (40-64) (Supplementary Table S1)"
    ),
    cyp2c9_freq = paste(
      "Phase III cohort: *1*1 66.1%, *1*2 18.4%, *1*3 12.3%, *2*2 1.4%, *2*3 1.2%,",
      "*3*3 0.6% (Supplementary Table S1)"
    ),
    vkorc1_freq = "Phase III cohort: GG 36.5%, GA 48.5%, AA 15% (Supplementary Table S1)",
    regions = "Not applicable (virtual population modelled after Hamberg 2010 and IWPC 2009 data)",
    notes = paste(
      "The clinical-trial model is a simulator. Its 53-state mechanistic core is the",
      "Wajima 2009 coagulation network; its population model combines Hartmann 2016 /",
      "2020 (production-rate IIV, VKORC1 / CYP2C9 effects) with Hamberg 2010 (EC50,",
      "elimination-rate, volume and age effects). The authors hand-calibrated the VKORC1",
      "effect on the vitamin K conversion rates against a baseline-INR distribution",
      "(Figure S1.1) and the GG EC50 and its IIV against IWPC 2009 maintenance doses",
      "(Figure S1.2). n_subjects is the size of the simulated MIPD-trial cohort."
    )
  )

  ini({
    # Values are the population-model parameters of the authors' deposit file
    # models/hartmann_coagulation_model_parameters.csv; "csv line N" is the line
    # number in that file. Parameter names follow Supplementary Table S1.3.
    # Nothing here was estimated: values were adopted from Wajima 2009 /
    # Hartmann 2016, 2020 / Hamberg 2010 or hand-calibrated (Section S1.2), so all
    # are fixed.

    # ---- Warfarin PK (Section 2.2 Eq 6; Hamberg 2010 structure, Section S1.2)
    lka <- fixed(log(2)); label("Warfarin absorption rate constant ka (1/h)")  # csv line 267
    lvc <- fixed(3.353); label("Log typical warfarin volume of distribution (log L)")  # csv line 263 (exp = 28.6 L)
    lkel <- fixed(-3.716); label("Log typical warfarin elimination rate, CYP2C9 *1/*1 at age 71 (log 1/h)")  # csv line 331; Eq S1.5
    e_cyp2c9s2_kel <- fixed(0.495); label("Fractional reduction of the warfarin elimination rate for CYP2C9 *2/*2 vs *1/*1 (half per *2 allele) (unitless)")  # csv line 333; Eq S1.5
    e_cyp2c9s3_kel <- fixed(0.757); label("Fractional reduction of the warfarin elimination rate for CYP2C9 *3/*3 vs *1/*1 (half per *3 allele) (unitless)")  # csv line 334; Eq S1.5
    e_age_kel <- fixed(0.00571); label("Age coefficient r_age of the warfarin elimination rate, 1 - tanh(r_age * (AGE - 71)) (1/year)")  # csv line 335; Eq S1.5

    # ---- Warfarin pharmacodynamics: inhibition of the VKORC-catalysed vitamin K
    # conversions (Table S1.1 vitamin K equations; Eq S1.3-S1.4, S1.8)
    lec50 <- fixed(-2); label("Log typical warfarin EC50 for VKORC inhibition, VKORC1 GG (log mg/L)")  # csv line 372; hand-calibrated against IWPC maintenance doses (Section S1.2, Figure S1.2)
    e_vkorc1aa_ec50 <- fixed(0.532); label("Fractional reduction of the warfarin EC50 for VKORC1 AA vs GG (half per A allele) (unitless)")  # csv line 374; Eq S1.3
    maximal_inhibitory_effect_warfarin <- fixed(1); label("Table S1.3 'maximal inhibitory effect warfarin' (unitless)")  # csv line 411
    hill_warfarin <- fixed(1.15); label("Table S1.3 'gamma': Hill coefficient of warfarin inhibition (unitless)")  # csv line 336; Section S4.4.2 gives gamma = 1.15 following Hamberg 2010
    lconversion_rate_vk_vkh2 <- fixed(-3.781); label("Log typical VKORC-catalysed conversion rate, vitamin K to vitamin K hydroquinone, VKORC1 GG (log 1/h)")  # csv line 279; Eq S1.8
    lconversion_rate_vko_vk <- fixed(-1.478); label("Log typical VKORC-catalysed conversion rate, vitamin K epoxide to vitamin K, VKORC1 GG (log 1/h)")  # csv line 283; Eq S1.8
    e_vkorc1aa_vk_vkh2 <- fixed(0.2); label("Fractional reduction of conversion_rate_vk_vkh2 for VKORC1 AA vs GG (half per A allele) (unitless)")  # csv line 281; 20% total reduction, Section S1.2 / Figure S1.1
    e_vkorc1aa_vko_vk <- fixed(0.2); label("Fractional reduction of conversion_rate_vko_vk for VKORC1 AA vs GG (half per A allele) (unitless)")  # csv line 285; 20% total reduction, Section S1.2 / Figure S1.1

    # ---- Compartment sizes
    v_blood <- fixed(24); label("Size of the central (blood) compartment of the coagulation network (L)")  # csv line 262
    v_vk_peripheral <- fixed(1); label("Size of the peripheral vitamin K compartment (L)")  # csv line 440

    # ---- Coagulation-factor production rates with normal IIV (Eq S1.2), population means
    tvproduction_rate_ii <- fixed(139.4); label("Table S1.3 'production rate ii', population mean (1/h)")  # csv line 413
    tvproduction_rate_v <- fixed(1.15); label("Table S1.3 'production rate v', population mean (nM/h)")  # csv line 425
    tvproduction_rate_vii <- fixed(12); label("Table S1.3 'production rate vii', population mean (1/h)")  # csv line 427
    tvproduction_rate_ix <- fixed(26.1); label("Table S1.3 'production rate ix', population mean (1/h)")  # csv line 415
    tvproduction_rate_x <- fixed(31.4); label("Table S1.3 'production rate x', population mean (1/h)")  # csv line 430
    tvproduction_rate_xi <- fixed(2.9); label("Table S1.3 'production rate xi', population mean (nM/h)")  # csv line 432
    tvproduction_rate_xii <- fixed(4.5); label("Table S1.3 'production rate xii', population mean (nM/h)")  # csv line 434
    tvproduction_rate_xiii <- fixed(0.25); label("Table S1.3 'production rate xiii', population mean (nM/h)")  # csv line 436
    tvproduction_rate_pc <- fixed(30.1); label("Table S1.3 'production rate pc', population mean (1/h)")  # csv line 417
    tvproduction_rate_ps <- fixed(49.8); label("Table S1.3 'production rate ps', population mean (1/h)")  # csv line 421

    # ---- Remaining mechanistic rate constants of the coagulation network (no IIV)
    degradation_rate_ii <- fixed(0.00999713138267); label("Table S1.3 'degradation rate ii' (1/h)")  # csv line 297
    degradation_rate_ii_activated <- fixed(67.4); label("Table S1.3 'degradation rate ii activated' (1/h)")  # csv line 298
    half_maximal_effect_concentration_av_ax_complex_ii <- fixed(10); label("Table S1.3 'half maximal effect concentration av ax complex ii' (nM)")  # csv line 346
    half_maximal_effect_concentration_ax_ii <- fixed(500); label("Table S1.3 'half maximal effect concentration ax ii' (nM)")  # csv line 351
    maximal_activation_rate_av_ax_complex_ii <- fixed(100); label("Table S1.3 'maximal activation rate av ax complex ii' (1/h)")  # csv line 380
    maximal_activation_rate_ax_ii <- fixed(9); label("Table S1.3 'maximal activation rate ax ii' (1/h)")  # csv line 385
    degradation_rate_v <- fixed(0.0430711610487); label("Table S1.3 'degradation rate v' (1/h)")  # csv line 312
    degradation_rate_v_activated <- fixed(20); label("Table S1.3 'degradation rate v activated' (1/h)")  # csv line 313
    half_maximal_effect_concentration_thrombin_v <- fixed(10); label("Table S1.3 'half maximal effect concentration thrombin v' (nM)")  # csv line 367
    half_maximal_effect_concentration_apc_ps_complex_av <- fixed(1); label("Table S1.3 'half maximal effect concentration apc ps complex av' (nM)")  # csv line 341
    maximal_activation_rate_thrombin_v <- fixed(50000); label("Table S1.3 'maximal activation rate thrombin v' (1/h)")  # csv line 393
    maximal_degradation_rate_apc_ps_complex_av <- fixed(100); label("Table S1.3 'maximal degradation rate apc ps complex av' (1/h)")  # csv line 404
    degradation_rate_vii <- fixed(0.12); label("Table S1.3 'degradation rate vii' (1/h)")  # csv line 314
    degradation_rate_vii_activated <- fixed(20); label("Table S1.3 'degradation rate vii activated' (1/h)")  # csv line 315
    half_maximal_effect_concentration_thrombin_vii <- fixed(10); label("Table S1.3 'half maximal effect concentration thrombin vii' (nM)")  # csv line 368
    half_maximal_effect_concentration_aix_vii <- fixed(10); label("Table S1.3 'half maximal effect concentration aix vii' (nM)")  # csv line 339
    half_maximal_effect_concentration_ax_vii <- fixed(10); label("Table S1.3 'half maximal effect concentration ax vii' (nM)")  # csv line 352
    half_maximal_effect_concentration_avii_tf_complex_vii <- fixed(10); label("Table S1.3 'half maximal effect concentration avii tf complex vii' (nM)")  # csv line 348
    maximal_activation_rate_thrombin_vii <- fixed(0.1); label("Table S1.3 'maximal activation rate thrombin vii' (1/h)")  # csv line 394
    maximal_activation_rate_aix_vii <- fixed(0.2); label("Table S1.3 'maximal activation rate aix vii' (1/h)")  # csv line 378
    maximal_activation_rate_ax_vii <- fixed(1); label("Table S1.3 'maximal activation rate ax vii' (1/h)")  # csv line 386
    maximal_activation_rate_avii_tf_complex_vii <- fixed(1); label("Table S1.3 'maximal activation rate avii tf complex vii' (1/h)")  # csv line 382
    degradation_rate_viii <- fixed(0.0571428571429); label("Table S1.3 'degradation rate viii' (1/h)")  # csv line 317
    degradation_rate_viii_activated <- fixed(20); label("Table S1.3 'degradation rate viii activated' (1/h)")  # csv line 318
    half_maximal_effect_concentration_thrombin_viii <- fixed(0.0001); label("Table S1.3 'half maximal effect concentration thrombin viii' (nM)")  # csv line 369
    half_maximal_effect_concentration_apc_ps_complex_aviii <- fixed(1); label("Table S1.3 'half maximal effect concentration apc ps complex aviii' (nM)")  # csv line 343
    maximal_activation_rate_thrombin_viii <- fixed(5000); label("Table S1.3 'maximal activation rate thrombin viii' (1/h)")  # csv line 395
    maximal_degradation_rate_apc_ps_complex_aviii <- fixed(50); label("Table S1.3 'maximal degradation rate apc ps complex aviii' (1/h)")  # csv line 406
    production_rate_viii <- fixed(0.04); label("Table S1.3 'production rate viii' (nM/h)")  # csv line 429
    degradation_rate_ix <- fixed(0.0291294642857); label("Table S1.3 'degradation rate ix' (1/h)")  # csv line 299
    degradation_rate_ix_activated <- fixed(20); label("Table S1.3 'degradation rate ix activated' (1/h)")  # csv line 300
    half_maximal_effect_concentration_avii_tf_complex_ix <- fixed(1); label("Table S1.3 'half maximal effect concentration avii tf complex ix' (nM)")  # csv line 347
    half_maximal_effect_concentration_axi_ix <- fixed(10); label("Table S1.3 'half maximal effect concentration axi ix' (nM)")  # csv line 354
    maximal_activation_rate_avii_tf_complex_ix <- fixed(70); label("Table S1.3 'maximal activation rate avii tf complex ix' (1/h)")  # csv line 381
    maximal_activation_rate_axi_ix <- fixed(7); label("Table S1.3 'maximal activation rate axi ix' (1/h)")  # csv line 388
    degradation_rate_x <- fixed(0.0180149168101); label("Table S1.3 'degradation rate x' (1/h)")  # csv line 319
    degradation_rate_x_activated <- fixed(20); label("Table S1.3 'degradation rate x activated' (1/h)")  # csv line 320
    half_maximal_effect_concentration_avii_x <- fixed(10); label("Table S1.3 'half maximal effect concentration avii x' (nM)")  # csv line 350
    half_maximal_effect_concentration_avii_tf_complex_x <- fixed(200); label("Table S1.3 'half maximal effect concentration avii tf complex x' (nM)")  # csv line 349
    half_maximal_effect_concentration_aix_x <- fixed(10); label("Table S1.3 'half maximal effect concentration aix x' (nM)")  # csv line 340
    half_maximal_effect_concentration_aix_aviii_complex_x <- fixed(0.2); label("Table S1.3 'half maximal effect concentration aix aviii complex x' (nM)")  # csv line 338
    maximal_activation_rate_avii_x <- fixed(1e-08); label("Table S1.3 'maximal activation rate avii x' (1/h)")  # csv line 384
    maximal_activation_rate_avii_tf_complex_x <- fixed(900); label("Table S1.3 'maximal activation rate avii tf complex x' (1/h)")  # csv line 383
    maximal_activation_rate_aix_x <- fixed(0.02); label("Table S1.3 'maximal activation rate aix x' (1/h)")  # csv line 379
    maximal_activation_rate_aix_aviii_complex_x <- fixed(2); label("Table S1.3 'maximal activation rate aix aviii complex x' (1/h)")  # csv line 377
    degradation_rate_xi <- fixed(0.0947712418301); label("Table S1.3 'degradation rate xi' (1/h)")  # csv line 323
    degradation_rate_xi_activated <- fixed(20); label("Table S1.3 'degradation rate xi activated' (1/h)")  # csv line 324
    half_maximal_effect_concentration_axii_xi <- fixed(1); label("Table S1.3 'half maximal effect concentration axii xi' (nM)")  # csv line 356
    half_maximal_effect_concentration_thrombin_xi <- fixed(10); label("Table S1.3 'half maximal effect concentration thrombin xi' (nM)")  # csv line 370
    maximal_activation_rate_axii_xi <- fixed(7); label("Table S1.3 'maximal activation rate axii xi' (1/h)")  # csv line 389
    maximal_activation_rate_thrombin_xi <- fixed(10); label("Table S1.3 'maximal activation rate thrombin xi' (1/h)")  # csv line 396
    degradation_rate_xii <- fixed(0.012); label("Table S1.3 'degradation rate xii' (1/h)")  # csv line 325
    degradation_rate_xii_activated <- fixed(20); label("Table S1.3 'degradation rate xii activated' (1/h)")  # csv line 326
    half_maximal_effect_concentration_ca_xii <- fixed(1); label("Table S1.3 'half maximal effect concentration ca xii' (nM)")  # csv line 358
    half_maximal_effect_concentration_kallikrein_xii <- fixed(1); label("Table S1.3 'half maximal effect concentration kallikrein xii' (nM)")  # csv line 360
    maximal_activation_rate_ca_xii <- fixed(7); label("Table S1.3 'maximal activation rate ca xii' (1/h)")  # csv line 390
    maximal_activation_rate_kallikrein_xii <- fixed(70); label("Table S1.3 'maximal activation rate kallikrein xii' (1/h)")  # csv line 391
    degradation_rate_xiii <- fixed(0.00355618776671); label("Table S1.3 'degradation rate xiii' (1/h)")  # csv line 327
    degradation_rate_xiii_activated <- fixed(0.69); label("Table S1.3 'degradation rate xiii activated' (1/h)")  # csv line 328
    half_maximal_effect_concentration_thrombin_xiii <- fixed(1); label("Table S1.3 'half maximal effect concentration thrombin xiii' (nM)")  # csv line 371
    maximal_activation_rate_thrombin_xiii <- fixed(7); label("Table S1.3 'maximal activation rate thrombin xiii' (1/h)")  # csv line 397
    degradation_rate_ca <- fixed(0.05); label("Table S1.3 'degradation rate ca' (1/h)")  # csv line 293
    degradation_rate_fibrinogen <- fixed(0.0320004471522); label("Table S1.3 'degradation rate fibrinogen' (1/h)")  # csv line 296
    degradation_rate_fibrin <- fixed(0.05); label("Table S1.3 'degradation rate fibrin' (1/h)")  # csv line 295
    degradation_rate_xf <- fixed(0.05); label("Table S1.3 'degradation rate xf' (1/h)")  # csv line 321
    half_maximal_effect_concentration_axiii_fibrin <- fixed(10); label("Table S1.3 'half maximal effect concentration axiii fibrin' (nM)")  # csv line 357
    half_maximal_effect_concentration_thrombin_fibrinogen <- fixed(0.5); label("Table S1.3 'half maximal effect concentration thrombin fibrinogen' (nM)")  # csv line 365
    half_maximal_effect_concentration_plasmin_fibrinogen <- fixed(500); label("Table S1.3 'half maximal effect concentration plasmin fibrinogen' (nM)")  # csv line 362
    half_maximal_effect_concentration_plasmin_fibrin <- fixed(10); label("Table S1.3 'half maximal effect concentration plasmin fibrin' (nM)")  # csv line 361
    half_maximal_effect_concentration_plasmin_xf <- fixed(100); label("Table S1.3 'half maximal effect concentration plasmin xf' (nM)")  # csv line 363
    half_maximal_effect_concentration_apc_ps_complex_xf <- fixed(1); label("Table S1.3 'half maximal effect concentration apc ps complex xf' (nM)")  # csv line 345
    maximal_conversion_rate_axiii_fibrin <- fixed(7); label("Table S1.3 'maximal conversion rate axiii fibrin' (1/h)")  # csv line 400
    maximal_conversion_rate_thrombin_fibrinogen <- fixed(20000); label("Table S1.3 'maximal conversion rate thrombin fibrinogen' (1/h)")  # csv line 402
    maximal_degradation_rate_plasmin_fibrinogen <- fixed(500); label("Table S1.3 'maximal degradation rate plasmin fibrinogen' (1/h)")  # csv line 409
    maximal_degradation_rate_plasmin_fibrin <- fixed(7); label("Table S1.3 'maximal degradation rate plasmin fibrin' (1/h)")  # csv line 408
    maximal_degradation_rate_plasmin_xf <- fixed(7); label("Table S1.3 'maximal degradation rate plasmin xf' (1/h)")  # csv line 410
    maximal_degradation_rate_apc_ps_complex_xf <- fixed(1); label("Table S1.3 'maximal degradation rate apc ps complex xf' (1/h)")  # csv line 407
    production_rate_fibrinogen <- fixed(286.26); label("Table S1.3 'production rate fibrinogen' (nM/h)")  # csv line 412
    degradation_rate_fdp <- fixed(3.5); label("Table S1.3 'degradation rate fdp' (1/h)")  # csv line 294
    degradation_rate_d_dimer <- fixed(0.1); label("Table S1.3 'degradation rate d dimer' (1/h)")  # csv line 322
    degradation_rate_kallikrein <- fixed(20); label("Table S1.3 'degradation rate kallikrein' (1/h)")  # csv line 301
    degradation_rate_prekallikrein <- fixed(0.05); label("Table S1.3 'degradation rate prekallikrein' (1/h)")  # csv line 306
    half_maximal_effect_concentration_axii_prekallikrein <- fixed(1); label("Table S1.3 'half maximal effect concentration axii prekallikrein' (nM)")  # csv line 355
    maximal_conversion_rate_axii_prekallikrein <- fixed(7); label("Table S1.3 'maximal conversion rate axii prekallikrein' (1/h)")  # csv line 399
    production_rate_prekallikrein <- fixed(22.5); label("Table S1.3 'production rate prekallikrein' (nM/h)")  # csv line 420
    degradation_rate_plasminogen <- fixed(0.0518497887945); label("Table S1.3 'degradation rate plasminogen' (1/h)")  # csv line 305
    degradation_rate_plasmin <- fixed(20); label("Table S1.3 'degradation rate plasmin' (1/h)")  # csv line 304
    half_maximal_effect_concentration_apc_ps_complex_plasminogen <- fixed(1); label("Table S1.3 'half maximal effect concentration apc ps complex plasminogen' (nM)")  # csv line 344
    half_maximal_effect_concentration_fibrin_plasminogen <- fixed(10000); label("Table S1.3 'half maximal effect concentration fibrin plasminogen' (nM)")  # csv line 359
    half_maximal_effect_concentration_thrombin_plasminogen <- fixed(5000); label("Table S1.3 'half maximal effect concentration thrombin plasminogen' (nM)")  # csv line 366
    maximal_conversion_rate_apc_ps_complex_plasminogen <- fixed(2); label("Table S1.3 'maximal conversion rate apc ps complex plasminogen' (1/h)")  # csv line 398
    maximal_conversion_rate_fibrin_plasminogen <- fixed(5); label("Table S1.3 'maximal conversion rate fibrin plasminogen' (1/h)")  # csv line 401
    maximal_conversion_rate_thrombin_plasminogen <- fixed(7); label("Table S1.3 'maximal conversion rate thrombin plasminogen' (1/h)")  # csv line 403
    production_rate_plasminogen <- fixed(111.7); label("Table S1.3 'production rate plasminogen' (nM/h)")  # csv line 419
    degradation_rate_pc <- fixed(0.0501666666667); label("Table S1.3 'degradation rate pc' (1/h)")  # csv line 302
    degradation_rate_pc_activated <- fixed(20.4); label("Table S1.3 'degradation rate pc activated' (1/h)")  # csv line 303
    half_maximal_effect_concentration_aii_tmod_complex_pc <- fixed(0.0001); label("Table S1.3 'half maximal effect concentration aii tmod complex pc' (nM)")  # csv line 337
    maximal_activation_rate_aii_tmod_complex_pc <- fixed(1); label("Table S1.3 'maximal activation rate aii tmod complex pc' (1/h)")  # csv line 376
    degradation_rate_ps <- fixed(0.0166); label("Table S1.3 'degradation rate ps' (1/h)")  # csv line 307
    degradation_rate_tmod <- fixed(0.05); label("Table S1.3 'degradation rate tmod' (1/h)")  # csv line 311
    production_rate_tmod <- fixed(2.5); label("Table S1.3 'production rate tmod' (nM/h)")  # csv line 424
    degradation_rate_tissue_factor <- fixed(0.05); label("Table S1.3 'degradation rate tissue factor' (1/h)")  # csv line 310
    conversion_rate_vkh2_vko <- fixed(0.228); label("Table S1.3 'conversion rate vkh2 vko' (1/h)")  # csv line 282
    elimination_rate_vk <- fixed(0.2052); label("Table S1.3 'elimination rate vk' (1/h)")  # csv line 330
    transition_rate_vk_central_to_peripheral <- fixed(0.0587); label("Table S1.3 'transition rate vk central to peripheral' (1/h)")  # csv line 438
    transition_rate_vk_peripheral_to_central <- fixed(0.0122); label("Table S1.3 'transition rate vk peripheral to central' (1/h)")  # csv line 439
    input_rate_vk <- fixed(0.2052); label("Table S1.3 'input rate vk' (nM/h)")  # csv line 375
    degradation_rate_tfpi <- fixed(20); label("Table S1.3 'degradation rate tfpi' (1/h)")  # csv line 309
    production_rate_tfpi <- fixed(50); label("Table S1.3 'production rate tfpi' (nM/h)")  # csv line 423
    degradation_rate_tat <- fixed(0.2); label("Table S1.3 'degradation rate tat' (1/h)")  # csv line 308
    complex_formation_rate_aii_atiii_heparin_complex <- fixed(8.35294117647); label("Table S1.3 'complex formation rate aii atiii heparin complex' (1/(nM*h))")  # csv line 268
    complex_formation_rate_aii_tmod <- fixed(0); label("Table S1.3 'complex formation rate aii tmod' (1/(nM*h))")  # csv line 269
    degradation_rate_aii_tmod_complex <- fixed(20); label("Table S1.3 'degradation rate aii tmod complex' (1/h)")  # csv line 286
    complex_formation_rate_av_ax <- fixed(2); label("Table S1.3 'complex formation rate av ax' (1/(nM*h))")  # csv line 273
    degradation_rate_av_ax_complex <- fixed(20); label("Table S1.3 'degradation rate av ax complex' (1/h)")  # csv line 289
    half_maximal_effect_concentration_apc_ps_complex_av_ax_complex <- fixed(1); label("Table S1.3 'half maximal effect concentration apc ps complex av ax complex' (nM)")  # csv line 342
    maximal_degradation_rate_apc_ps_complex_av_ax_complex <- fixed(2); label("Table S1.3 'maximal degradation rate apc ps complex av ax complex' (1/h)")  # csv line 405
    complex_formation_rate_vii_tf <- fixed(2); label("Table S1.3 'complex formation rate vii tf' (1/(nM*h))")  # csv line 278
    degradation_rate_vii_tf_complex <- fixed(7); label("Table S1.3 'degradation rate vii tf complex' (1/h)")  # csv line 316
    half_maximal_effect_concentration_ax_vii_tf_complex <- fixed(1); label("Table S1.3 'half maximal effect concentration ax vii tf complex' (nM)")  # csv line 353
    half_maximal_effect_concentration_tf_vii_tf_complex <- fixed(1); label("Table S1.3 'half maximal effect concentration tf vii tf complex' (nM)")  # csv line 364
    maximal_activation_rate_ax_vii_tf_complex <- fixed(70); label("Table S1.3 'maximal activation rate ax vii tf complex' (1/h)")  # csv line 387
    maximal_activation_rate_tf_vii_tf_complex <- fixed(1000); label("Table S1.3 'maximal activation rate tf vii tf complex' (1/h)")  # csv line 392
    complex_formation_rate_avii_tf <- fixed(2); label("Table S1.3 'complex formation rate avii tf' (1/(nM*h))")  # csv line 274
    degradation_rate_avii_tf_complex <- fixed(20); label("Table S1.3 'degradation rate avii tf complex' (1/h)")  # csv line 291
    complex_formation_rate_aix_atiii_heparin_complex <- fixed(1.17647058824); label("Table S1.3 'complex formation rate aix atiii heparin complex' (1/(nM*h))")  # csv line 270
    complex_formation_rate_aix_aviii <- fixed(100); label("Table S1.3 'complex formation rate aix aviii' (1/(nM*h))")  # csv line 271
    degradation_rate_aix_aviii_complex <- fixed(20); label("Table S1.3 'degradation rate aix aviii complex' (1/h)")  # csv line 287
    complex_formation_rate_ax_atiii_heparin_complex <- fixed(1.17647058824); label("Table S1.3 'complex formation rate ax atiii heparin complex' (1/(nM*h))")  # csv line 276
    complex_formation_rate_ax_tfpi <- fixed(2); label("Table S1.3 'complex formation rate ax tfpi' (1/(nM*h))")  # csv line 277
    degradation_rate_ax_tfpi_complex <- fixed(20); label("Table S1.3 'degradation rate ax tfpi complex' (1/h)")  # csv line 292
    complex_formation_rate_apc_ps <- fixed(1); label("Table S1.3 'complex formation rate apc ps' (1/(nM*h))")  # csv line 272
    degradation_rate_apc_ps_complex <- fixed(20); label("Table S1.3 'degradation rate apc ps complex' (1/h)")  # csv line 288
    complex_formation_rate_avii_tf_ax_tfpi <- fixed(2); label("Table S1.3 'complex formation rate avii tf ax tfpi' (1/(nM*h))")  # csv line 275
    degradation_rate_avii_tf_ax_tfpi_complex <- fixed(20); label("Table S1.3 'degradation rate avii tf ax tfpi complex' (1/h)")  # csv line 290
    elimination_rate_atiii_heparin_complex <- fixed(0.693); label("Table S1.3 'elimination rate atiii heparin complex' (1/h)")  # csv line 329

    # ---- Initial concentrations of the network states (csv amount / compartment size)
    bl_coagulation_factor_ii <- fixed(1394.4); label("Initial concentration of coagulation factor ii at the start of the drug-free run-in (nM)")  # csv line 233: 33465.6 nmol / 24 L
    bl_coagulation_factor_v <- fixed(26.7); label("Initial concentration of coagulation factor v at the start of the drug-free run-in (nM)")  # csv line 235: 640.8 nmol / 24 L
    bl_coagulation_factor_vii <- fixed(10); label("Initial concentration of coagulation factor vii at the start of the drug-free run-in (nM)")  # csv line 236: 240 nmol / 24 L
    bl_coagulation_factor_viii <- fixed(0.7); label("Initial concentration of coagulation factor viii at the start of the drug-free run-in (nM)")  # csv line 237: 16.8 nmol / 24 L
    bl_coagulation_factor_ix <- fixed(89.6); label("Initial concentration of coagulation factor ix at the start of the drug-free run-in (nM)")  # csv line 234: 2150.4 nmol / 24 L
    bl_coagulation_factor_x <- fixed(174.3); label("Initial concentration of coagulation factor x at the start of the drug-free run-in (nM)")  # csv line 238: 4183.2 nmol / 24 L
    bl_coagulation_factor_xi <- fixed(30.6); label("Initial concentration of coagulation factor xi at the start of the drug-free run-in (nM)")  # csv line 239: 734.4 nmol / 24 L
    bl_coagulation_factor_xii <- fixed(375); label("Initial concentration of coagulation factor xii at the start of the drug-free run-in (nM)")  # csv line 240: 9000 nmol / 24 L
    bl_coagulation_factor_xiii <- fixed(70.3); label("Initial concentration of coagulation factor xiii at the start of the drug-free run-in (nM)")  # csv line 241: 1687.2 nmol / 24 L
    bl_fibrinogen <- fixed(8945.5); label("Initial concentration of fibrinogen at the start of the drug-free run-in (nM)")  # csv line 245: 214692 nmol / 24 L
    bl_prekallikrein <- fixed(450); label("Initial concentration of prekallikrein at the start of the drug-free run-in (nM)")  # csv line 249: 10800 nmol / 24 L
    bl_plasminogen <- fixed(2154.3); label("Initial concentration of plasminogen at the start of the drug-free run-in (nM)")  # csv line 248: 51703.2 nmol / 24 L
    bl_protein_c <- fixed(60); label("Initial concentration of protein c at the start of the drug-free run-in (nM)")  # csv line 250: 1440 nmol / 24 L
    bl_protein_s <- fixed(300); label("Initial concentration of protein s at the start of the drug-free run-in (nM)")  # csv line 251: 7200 nmol / 24 L
    bl_thrombomodulin <- fixed(50); label("Initial concentration of thrombomodulin at the start of the drug-free run-in (nM)")  # csv line 252: 1200 nmol / 24 L
    bl_tissue_factor_pathway_inhibitor <- fixed(2.5); label("Initial concentration of tissue factor pathway inhibitor at the start of the drug-free run-in (nM)")  # csv line 254: 60 nmol / 24 L
    bl_vitamin_k <- fixed(1); label("Initial concentration of vitamin k at the start of the drug-free run-in (nM)")  # csv line 256: 24 nmol / 24 L
    bl_vitamin_k_epoxide <- fixed(0.1); label("Initial concentration of vitamin k epoxide at the start of the drug-free run-in (nM)")  # csv line 257: 2.4 nmol / 24 L
    bl_vitamin_k_hydroquinone <- fixed(0.1); label("Initial concentration of vitamin k hydroquinone at the start of the drug-free run-in (nM)")  # csv line 258: 2.4 nmol / 24 L
    bl_vitamin_k_peripheral <- fixed(115.475409836); label("Initial concentration of vitamin k peripheral at the start of the drug-free run-in (nM)")  # csv line 261: 115.475409836 nmol / 1 L

    # ---- IIV (Section S1.2). Lognormal parameters: variance = (csv log std)^2.
    etalvc ~ fixed(0.002916)  # csv line 264: log std 0.054
    etalkel ~ fixed(0.005041)  # csv line 332: log std 0.071
    etalec50 ~ fixed(0.01)  # csv line 373: log std 0.1
    etalconversion_rate_vk_vkh2 ~ fixed(0.0576)  # csv line 280: log std 0.24
    etalconversion_rate_vko_vk ~ fixed(0.0576)  # csv line 284: log std 0.24
    # Normal IIV on production rates, written as a proportional deviation:
    # rate = mean * (1 + eta), variance = (csv std / csv mean)^2 (20% CV, Eq S1.2).
    etatvproduction_rate_ii ~ fixed(0.04)  # csv lines 413-414: std 27.88 / mean 139.4
    etatvproduction_rate_v ~ fixed(0.037602268431)  # csv lines 425-426: std 0.223 / mean 1.15
    etatvproduction_rate_vii ~ fixed(0.04)  # csv lines 427-428: std 2.4 / mean 12
    etatvproduction_rate_ix ~ fixed(0.04)  # csv lines 415-416: std 5.22 / mean 26.1
    etatvproduction_rate_x ~ fixed(0.04)  # csv lines 430-431: std 6.28 / mean 31.4
    etatvproduction_rate_xi ~ fixed(0.04)  # csv lines 432-433: std 0.58 / mean 2.9
    etatvproduction_rate_xii ~ fixed(0.04)  # csv lines 434-435: std 0.9 / mean 4.5
    etatvproduction_rate_xiii ~ fixed(0.04)  # csv lines 436-437: std 0.05 / mean 0.25
    etatvproduction_rate_pc ~ fixed(0.04)  # csv lines 417-418: std 6.02 / mean 30.1
    etatvproduction_rate_ps ~ fixed(0.04)  # csv lines 421-422: std 9.96 / mean 49.8

    # ---- Inter-occasion variability of the vitamin K input rate (Section 2.2 Eq 7):
    # multiplier ~ N(1, 0.1^2), redrawn every treatment day. One eta per day,
    # OCC = 1-56 (the longest simulated trial, phase III, ran 8 weeks).
    etaiov_vk_input_1 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_2 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_3 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_4 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_5 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_6 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_7 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_8 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_9 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_10 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_11 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_12 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_13 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_14 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_15 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_16 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_17 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_18 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_19 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_20 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_21 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_22 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_23 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_24 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_25 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_26 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_27 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_28 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_29 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_30 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_31 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_32 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_33 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_34 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_35 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_36 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_37 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_38 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_39 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_40 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_41 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_42 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_43 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_44 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_45 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_46 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_47 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_48 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_49 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_50 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_51 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_52 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_53 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_54 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_55 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1
    etaiov_vk_input_56 ~ fixed(0.01)  # Section 2.2: sigma_eta = 0.1

    # ---- Measurement model (Section 2.2 Eq 9): lognormal, scale 0.1
    expSd <- fixed(0.1); label("Lognormal measurement-error scale, applied to warfarin concentration and to INR (unitless)")  # csv line 441; Eq 9 sigma = 0.1
  })
  model({
    # VKORC1 -1639A allele count; GG is the reference genotype (Eq S1.3)
    n_vkorc1_a <- 2 - VKORC1_1639G_COUNT

    # Warfarin PK. CYP2C9 enters as per-allele fractional reductions and age as
    # 1 - tanh(r_age * (AGE - 71)), the authors' bounded form of Hamberg 2010 (Eq S1.5).
    ka <- exp(lka)
    vc <- exp(lvc + etalvc)
    kel <- exp(lkel + etalkel) *
      (1 - (e_cyp2c9s2_kel * CYP2C9_S2_COUNT + e_cyp2c9s3_kel * CYP2C9_S3_COUNT) / 2) *
      (1 - tanh(e_age_kel * (AGE - 71)))

    # VKORC1 genotype on the warfarin EC50 and on both VKORC-catalysed conversions (Eq S1.3, S1.8)
    ec50 <- exp(lec50 + etalec50) * (1 - e_vkorc1aa_ec50 * n_vkorc1_a / 2)
    conversion_rate_vk_vkh2 <- exp(lconversion_rate_vk_vkh2 + etalconversion_rate_vk_vkh2) *
      (1 - e_vkorc1aa_vk_vkh2 * n_vkorc1_a / 2)
    conversion_rate_vko_vk <- exp(lconversion_rate_vko_vk + etalconversion_rate_vko_vk) *
      (1 - e_vkorc1aa_vko_vk * n_vkorc1_a / 2)

    # Coagulation-factor production rates (Eq S1.2)
    production_rate_ii <- tvproduction_rate_ii * (1 + etatvproduction_rate_ii)
    production_rate_v <- tvproduction_rate_v * (1 + etatvproduction_rate_v)
    production_rate_vii <- tvproduction_rate_vii * (1 + etatvproduction_rate_vii)
    production_rate_ix <- tvproduction_rate_ix * (1 + etatvproduction_rate_ix)
    production_rate_x <- tvproduction_rate_x * (1 + etatvproduction_rate_x)
    production_rate_xi <- tvproduction_rate_xi * (1 + etatvproduction_rate_xi)
    production_rate_xii <- tvproduction_rate_xii * (1 + etatvproduction_rate_xii)
    production_rate_xiii <- tvproduction_rate_xiii * (1 + etatvproduction_rate_xiii)
    production_rate_pc <- tvproduction_rate_pc * (1 + etatvproduction_rate_pc)
    production_rate_ps <- tvproduction_rate_ps * (1 + etatvproduction_rate_ps)

    # Daily vitamin K input with inter-occasion variability; OCC = 0 (the run-in,
    # or any record outside treatment days 1-56) carries no IOV.
    iov_vk_input <- (OCC == 1) * etaiov_vk_input_1 +
      (OCC == 2) * etaiov_vk_input_2 +
      (OCC == 3) * etaiov_vk_input_3 +
      (OCC == 4) * etaiov_vk_input_4 +
      (OCC == 5) * etaiov_vk_input_5 +
      (OCC == 6) * etaiov_vk_input_6 +
      (OCC == 7) * etaiov_vk_input_7 +
      (OCC == 8) * etaiov_vk_input_8 +
      (OCC == 9) * etaiov_vk_input_9 +
      (OCC == 10) * etaiov_vk_input_10 +
      (OCC == 11) * etaiov_vk_input_11 +
      (OCC == 12) * etaiov_vk_input_12 +
      (OCC == 13) * etaiov_vk_input_13 +
      (OCC == 14) * etaiov_vk_input_14 +
      (OCC == 15) * etaiov_vk_input_15 +
      (OCC == 16) * etaiov_vk_input_16 +
      (OCC == 17) * etaiov_vk_input_17 +
      (OCC == 18) * etaiov_vk_input_18 +
      (OCC == 19) * etaiov_vk_input_19 +
      (OCC == 20) * etaiov_vk_input_20 +
      (OCC == 21) * etaiov_vk_input_21 +
      (OCC == 22) * etaiov_vk_input_22 +
      (OCC == 23) * etaiov_vk_input_23 +
      (OCC == 24) * etaiov_vk_input_24 +
      (OCC == 25) * etaiov_vk_input_25 +
      (OCC == 26) * etaiov_vk_input_26 +
      (OCC == 27) * etaiov_vk_input_27 +
      (OCC == 28) * etaiov_vk_input_28 +
      (OCC == 29) * etaiov_vk_input_29 +
      (OCC == 30) * etaiov_vk_input_30 +
      (OCC == 31) * etaiov_vk_input_31 +
      (OCC == 32) * etaiov_vk_input_32 +
      (OCC == 33) * etaiov_vk_input_33 +
      (OCC == 34) * etaiov_vk_input_34 +
      (OCC == 35) * etaiov_vk_input_35 +
      (OCC == 36) * etaiov_vk_input_36 +
      (OCC == 37) * etaiov_vk_input_37 +
      (OCC == 38) * etaiov_vk_input_38 +
      (OCC == 39) * etaiov_vk_input_39 +
      (OCC == 40) * etaiov_vk_input_40 +
      (OCC == 41) * etaiov_vk_input_41 +
      (OCC == 42) * etaiov_vk_input_42 +
      (OCC == 43) * etaiov_vk_input_43 +
      (OCC == 44) * etaiov_vk_input_44 +
      (OCC == 45) * etaiov_vk_input_45 +
      (OCC == 46) * etaiov_vk_input_46 +
      (OCC == 47) * etaiov_vk_input_47 +
      (OCC == 48) * etaiov_vk_input_48 +
      (OCC == 49) * etaiov_vk_input_49 +
      (OCC == 50) * etaiov_vk_input_50 +
      (OCC == 51) * etaiov_vk_input_51 +
      (OCC == 52) * etaiov_vk_input_52 +
      (OCC == 53) * etaiov_vk_input_53 +
      (OCC == 54) * etaiov_vk_input_54 +
      (OCC == 55) * etaiov_vk_input_55 +
      (OCC == 56) * etaiov_vk_input_56
    vk_input <- input_rate_vk * (1 + iov_vk_input)

    # Warfarin plasma concentration (mg/L)
    Cc <- central / vc

    # Reaction rates of the Wajima 2009 network as deposited in
    # models/wajima_coagulation_model.xml (reaction ids kept), expressed per litre
    # of blood (nM/h); species names denote concentrations (nM).
    r_atiii_heparin_complex_elimination <- elimination_rate_atiii_heparin_complex * atiii_heparin_complex
    r_prothrombin_production <- production_rate_ii * vitamin_k_hydroquinone
    r_prothrombin_degradation <- degradation_rate_ii * coagulation_factor_ii
    r_av_ax_complex_prothrombin_activation <- (maximal_activation_rate_av_ax_complex_ii * av_ax_complex * coagulation_factor_ii) / (half_maximal_effect_concentration_av_ax_complex_ii + av_ax_complex)
    r_ax_prothrombin_activation <- (maximal_activation_rate_ax_ii * activated_coagulation_factor_x * coagulation_factor_ii) / (half_maximal_effect_concentration_ax_ii + activated_coagulation_factor_x)
    r_thrombin_degradation <- degradation_rate_ii_activated * activated_coagulation_factor_ii
    r_aii_atiii_heparin_complex_formation <- complex_formation_rate_aii_atiii_heparin_complex * activated_coagulation_factor_ii * atiii_heparin_complex
    r_tat_degradation <- degradation_rate_tat * aii_at_complex
    r_aii_tmod_complex_formation <- complex_formation_rate_aii_tmod * activated_coagulation_factor_ii * thrombomodulin
    r_aii_tmod_complex_degradation <- degradation_rate_aii_tmod_complex * aii_tmod_complex
    r_v_production <- production_rate_v
    r_v_degradation <- degradation_rate_v * coagulation_factor_v
    r_thrombin_v_activation <- (maximal_activation_rate_thrombin_v * activated_coagulation_factor_ii * coagulation_factor_v) / (half_maximal_effect_concentration_thrombin_v + activated_coagulation_factor_ii)
    r_av_degradation <- degradation_rate_v_activated * activated_coagulation_factor_v
    r_apc_ps_complex_av_degradation <- (maximal_degradation_rate_apc_ps_complex_av * apc_ps_complex * activated_coagulation_factor_v) / (half_maximal_effect_concentration_apc_ps_complex_av + apc_ps_complex)
    r_av_ax_complex_formation <- complex_formation_rate_av_ax * activated_coagulation_factor_v * activated_coagulation_factor_x
    r_av_ax_complex_degradation <- degradation_rate_av_ax_complex * av_ax_complex
    r_apc_ps_complex_av_ax_complex_degradation <- (maximal_degradation_rate_apc_ps_complex_av_ax_complex * apc_ps_complex * av_ax_complex) / (half_maximal_effect_concentration_apc_ps_complex_av_ax_complex + apc_ps_complex)
    r_vii_production <- production_rate_vii * vitamin_k_hydroquinone
    r_vii_degradation <- degradation_rate_vii * coagulation_factor_vii
    r_thrombin_vii_activation <- (maximal_activation_rate_thrombin_vii * activated_coagulation_factor_ii * coagulation_factor_vii) / (half_maximal_effect_concentration_thrombin_vii + activated_coagulation_factor_ii)
    r_aix_vii_activation <- (maximal_activation_rate_aix_vii * activated_coagulation_factor_ix * coagulation_factor_vii) / (half_maximal_effect_concentration_aix_vii + activated_coagulation_factor_ix)
    r_ax_vii_activation <- (maximal_activation_rate_ax_vii * activated_coagulation_factor_x * coagulation_factor_vii) / (half_maximal_effect_concentration_ax_vii + activated_coagulation_factor_x)
    r_avii_tf_complex_vii_activation <- (maximal_activation_rate_avii_tf_complex_vii * avii_tf_complex * coagulation_factor_vii) / (half_maximal_effect_concentration_avii_tf_complex_vii + avii_tf_complex)
    r_vii_tf_complex_formation <- complex_formation_rate_vii_tf * coagulation_factor_vii * tissue_factor
    r_vii_tf_complex_degradation <- degradation_rate_vii_tf_complex * vii_tf_complex
    r_ax_vii_tf_complex_activation <- (maximal_activation_rate_ax_vii_tf_complex * activated_coagulation_factor_x * vii_tf_complex) / (half_maximal_effect_concentration_ax_vii_tf_complex + activated_coagulation_factor_x)
    r_tf_vii_tf_complex_activation <- (maximal_activation_rate_tf_vii_tf_complex * tissue_factor * vii_tf_complex) / (half_maximal_effect_concentration_tf_vii_tf_complex + tissue_factor)
    r_activated_vii_degradation <- degradation_rate_vii_activated * activated_coagulation_factor_vii
    r_avii_tf_complex_formation <- complex_formation_rate_avii_tf * activated_coagulation_factor_vii * tissue_factor
    r_avii_tf_complex_degradation <- degradation_rate_avii_tf_complex * avii_tf_complex
    r_avii_tf_ax_tfpi_complex_formation <- complex_formation_rate_avii_tf_ax_tfpi * avii_tf_complex * ax_tfpi_complex
    r_avii_tf_ax_tfpi_complex_degradation <- degradation_rate_avii_tf_ax_tfpi_complex * avii_tf_ax_tfpi_complex
    r_viii_production <- production_rate_viii
    r_viii_degradation <- degradation_rate_viii * coagulation_factor_viii
    r_thrombin_viii_activation <- (maximal_activation_rate_thrombin_viii * activated_coagulation_factor_ii * coagulation_factor_viii) / (half_maximal_effect_concentration_thrombin_viii + activated_coagulation_factor_ii)
    r_activated_viii_degradation <- degradation_rate_viii_activated * activated_coagulation_factor_viii
    r_apc_ps_complex_activated_viii_degradation <- (maximal_degradation_rate_apc_ps_complex_aviii * apc_ps_complex * activated_coagulation_factor_viii) / (half_maximal_effect_concentration_apc_ps_complex_aviii + apc_ps_complex)
    r_ix_production <- production_rate_ix * vitamin_k_hydroquinone
    r_ix_degradation <- degradation_rate_ix * coagulation_factor_ix
    r_axi_ix_activation <- (maximal_activation_rate_axi_ix * activated_coagulation_factor_xi * coagulation_factor_ix) / (half_maximal_effect_concentration_axi_ix + activated_coagulation_factor_xi)
    r_avii_tf_complex_ix_activation <- (maximal_activation_rate_avii_tf_complex_ix * avii_tf_complex * coagulation_factor_ix) / (half_maximal_effect_concentration_avii_tf_complex_ix + avii_tf_complex)
    r_activated_ix_degradation <- degradation_rate_ix_activated * activated_coagulation_factor_ix
    r_aix_aviii_complex_formation <- complex_formation_rate_aix_aviii * activated_coagulation_factor_ix * activated_coagulation_factor_viii
    r_aix_aviii_complex_degradation <- degradation_rate_aix_aviii_complex * aix_aviii_complex
    r_aix_atiii_heparin_complex_formation <- complex_formation_rate_aix_atiii_heparin_complex * activated_coagulation_factor_ix * atiii_heparin_complex
    r_x_production <- production_rate_x * vitamin_k_hydroquinone
    r_x_degradation <- degradation_rate_x * coagulation_factor_x
    r_avii_x_activation <- (maximal_activation_rate_avii_x * activated_coagulation_factor_vii * coagulation_factor_x) / (half_maximal_effect_concentration_avii_x + activated_coagulation_factor_vii)
    r_avii_tf_complex_x_activation <- (maximal_activation_rate_avii_tf_complex_x * avii_tf_complex * coagulation_factor_x) / (half_maximal_effect_concentration_avii_tf_complex_x + avii_tf_complex)
    r_aix_x_activation <- (maximal_activation_rate_aix_x * activated_coagulation_factor_ix * coagulation_factor_x) / (half_maximal_effect_concentration_aix_x + activated_coagulation_factor_ix)
    r_aix_aviii_complex_x_activation <- (maximal_activation_rate_aix_aviii_complex_x * aix_aviii_complex * coagulation_factor_x) / (half_maximal_effect_concentration_aix_aviii_complex_x + aix_aviii_complex)
    r_ax_degradation <- degradation_rate_x_activated * activated_coagulation_factor_x
    r_ax_atiii_heparin_complex_complex_formation <- complex_formation_rate_ax_atiii_heparin_complex * activated_coagulation_factor_x * atiii_heparin_complex
    r_ax_tfpi_complex_formation <- complex_formation_rate_ax_tfpi * activated_coagulation_factor_x * tissue_factor_pathway_inhibitor
    r_ax_tfpi_complex_degradation <- degradation_rate_ax_tfpi_complex * ax_tfpi_complex
    r_xi_production <- production_rate_xi
    r_xi_degradation <- degradation_rate_xi * coagulation_factor_xi
    r_thrombin_xi_activation <- (maximal_activation_rate_thrombin_xi * activated_coagulation_factor_ii * coagulation_factor_xi) / (half_maximal_effect_concentration_thrombin_xi + activated_coagulation_factor_ii)
    r_axii_xi_activation <- (maximal_activation_rate_axii_xi * activated_coagulation_factor_xii * coagulation_factor_xi) / (half_maximal_effect_concentration_axii_xi + activated_coagulation_factor_xii)
    r_activated_xi_degradation <- degradation_rate_xi_activated * activated_coagulation_factor_xi
    r_xii_production <- production_rate_xii
    r_xii_degradation <- degradation_rate_xii * coagulation_factor_xii
    r_ca_xii_activation <- (maximal_activation_rate_ca_xii * contact_system_activator * coagulation_factor_xii) / (half_maximal_effect_concentration_ca_xii + contact_system_activator)
    r_kallikrein_xii_activation <- (maximal_activation_rate_kallikrein_xii * kallikrein * coagulation_factor_xii) / (half_maximal_effect_concentration_kallikrein_xii + kallikrein)
    r_activated_xii_degradation <- degradation_rate_xii_activated * activated_coagulation_factor_xii
    r_xiii_production <- production_rate_xiii
    r_xiii_degradation <- degradation_rate_xiii * coagulation_factor_xiii
    r_thrombin_xiii_activation <- (maximal_activation_rate_thrombin_xiii * activated_coagulation_factor_ii * coagulation_factor_xiii) / (half_maximal_effect_concentration_thrombin_xiii + activated_coagulation_factor_ii)
    r_axiii_degradation <- degradation_rate_xiii_activated * activated_coagulation_factor_xiii
    r_ca_degradation <- degradation_rate_ca * contact_system_activator
    r_fibrinogen_production <- production_rate_fibrinogen
    r_fibrinogen_degradation <- degradation_rate_fibrinogen * fibrinogen
    r_thrombin_fibrinogen_conversion <- (maximal_conversion_rate_thrombin_fibrinogen * activated_coagulation_factor_ii * fibrinogen) / (half_maximal_effect_concentration_thrombin_fibrinogen + activated_coagulation_factor_ii)
    r_plasmin_fibrinogen_degradation <- (maximal_degradation_rate_plasmin_fibrinogen * plasmin * fibrinogen) / (half_maximal_effect_concentration_plasmin_fibrinogen + plasmin)
    r_fibrin_degradation <- degradation_rate_fibrin * fibrin
    r_axiii_fibrin_conversion <- (maximal_conversion_rate_axiii_fibrin * activated_coagulation_factor_xiii * fibrin) / (half_maximal_effect_concentration_axiii_fibrin + activated_coagulation_factor_xiii)
    r_plasmin_fibrin_degradation <- (maximal_degradation_rate_plasmin_fibrin * plasmin * fibrin) / (half_maximal_effect_concentration_plasmin_fibrin + plasmin)
    r_xf_degradation <- degradation_rate_xf * cross_linked_fibrin
    r_plasmin_xf_degradation <- (maximal_degradation_rate_plasmin_xf * plasmin * cross_linked_fibrin) / (half_maximal_effect_concentration_plasmin_xf + plasmin)
    r_apc_ps_complex_xf_degradation <- (maximal_degradation_rate_apc_ps_complex_xf * apc_ps_complex * cross_linked_fibrin) / (half_maximal_effect_concentration_apc_ps_complex_xf + apc_ps_complex)
    r_fdp_degradation <- degradation_rate_fdp * fibrin_degradation_product
    r_d_dimer_degradation <- degradation_rate_d_dimer * d_dimer
    r_prekallikrein_production <- production_rate_prekallikrein
    r_prekallikrein_degradation <- degradation_rate_prekallikrein * prekallikrein
    r_axii_prekallikrein_conversion <- (maximal_conversion_rate_axii_prekallikrein * activated_coagulation_factor_xii * prekallikrein) / (half_maximal_effect_concentration_axii_prekallikrein + activated_coagulation_factor_xii)
    r_kallikrein_degradation <- degradation_rate_kallikrein * kallikrein
    r_plasminogen_production <- production_rate_plasminogen
    r_plasminogen_degradation <- degradation_rate_plasminogen * plasminogen
    r_thrombin_plasminogen_conversion <- (maximal_conversion_rate_thrombin_plasminogen * activated_coagulation_factor_ii * plasminogen) / (half_maximal_effect_concentration_thrombin_plasminogen + activated_coagulation_factor_ii)
    r_fibrin_plasminogen_conversion <- (maximal_conversion_rate_fibrin_plasminogen * fibrin * plasminogen) / (half_maximal_effect_concentration_fibrin_plasminogen + fibrin)
    r_apc_ps_complex_plasminogen_conversion <- (maximal_conversion_rate_apc_ps_complex_plasminogen * apc_ps_complex * plasminogen) / (half_maximal_effect_concentration_apc_ps_complex_plasminogen + apc_ps_complex)
    r_plasmin_degradation <- degradation_rate_plasmin * plasmin
    r_protein_c_production <- production_rate_pc * vitamin_k_hydroquinone
    r_protein_c_degradation <- degradation_rate_pc * protein_c
    r_aii_tmod_complex_pc_activation <- (maximal_activation_rate_aii_tmod_complex_pc * aii_tmod_complex * protein_c) / (half_maximal_effect_concentration_aii_tmod_complex_pc + aii_tmod_complex)
    r_activated_protein_c_degradation <- degradation_rate_pc_activated * activated_protein_c
    r_apc_ps_complex_formation <- complex_formation_rate_apc_ps * activated_protein_c * protein_s
    r_apc_ps_complex_degradation <- degradation_rate_apc_ps_complex * apc_ps_complex
    r_protein_s_production <- production_rate_ps * vitamin_k_hydroquinone
    r_protein_s_degradation <- degradation_rate_ps * protein_s
    r_thrombomodulin_production <- production_rate_tmod
    r_thrombomodulin_degradation <- degradation_rate_tmod * thrombomodulin
    r_tissue_factor_degradation <- degradation_rate_tissue_factor * tissue_factor
    r_tfpi_production <- production_rate_tfpi
    r_tfpi_degradation <- degradation_rate_tfpi * tissue_factor_pathway_inhibitor
    r_vitamin_k_intake <- vk_input
    r_vk_degradation <- elimination_rate_vk * vitamin_k
    r_vitamin_k_conversion <- conversion_rate_vk_vkh2 * vitamin_k * (1 - (maximal_inhibitory_effect_warfarin * Cc^hill_warfarin) / (ec50^hill_warfarin + Cc^hill_warfarin))
    r_vkh2_conversion <- conversion_rate_vkh2_vko * vitamin_k_hydroquinone
    r_vitamin_k_epoxide_conversion <- conversion_rate_vko_vk * vitamin_k_epoxide * (1 - (maximal_inhibitory_effect_warfarin * Cc^hill_warfarin) / (ec50^hill_warfarin + Cc^hill_warfarin))
    r_vk_transition_to_peripheral <- transition_rate_vk_central_to_peripheral * vitamin_k
    r_vk_transition_to_central <- v_vk_peripheral / v_blood * transition_rate_vk_peripheral_to_central * vitamin_k_peripheral

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(coagulation_factor_ii) <- r_prothrombin_production - r_prothrombin_degradation - r_av_ax_complex_prothrombin_activation - r_ax_prothrombin_activation
    d/dt(activated_coagulation_factor_ii) <- r_av_ax_complex_prothrombin_activation + r_ax_prothrombin_activation - r_thrombin_degradation - r_aii_atiii_heparin_complex_formation - r_aii_tmod_complex_formation
    d/dt(coagulation_factor_v) <- r_v_production - r_v_degradation - r_thrombin_v_activation
    d/dt(activated_coagulation_factor_v) <- r_thrombin_v_activation - r_av_degradation - r_apc_ps_complex_av_degradation - r_av_ax_complex_formation
    d/dt(coagulation_factor_vii) <- r_vii_production - r_vii_degradation - r_thrombin_vii_activation - r_aix_vii_activation - r_ax_vii_activation - r_avii_tf_complex_vii_activation - r_vii_tf_complex_formation
    d/dt(activated_coagulation_factor_vii) <- r_thrombin_vii_activation + r_aix_vii_activation + r_ax_vii_activation + r_avii_tf_complex_vii_activation - r_activated_vii_degradation - r_avii_tf_complex_formation
    d/dt(coagulation_factor_viii) <- r_viii_production - r_viii_degradation - r_thrombin_viii_activation
    d/dt(activated_coagulation_factor_viii) <- r_thrombin_viii_activation - r_activated_viii_degradation - r_apc_ps_complex_activated_viii_degradation - r_aix_aviii_complex_formation
    d/dt(coagulation_factor_ix) <- r_ix_production - r_ix_degradation - r_axi_ix_activation - r_avii_tf_complex_ix_activation
    d/dt(activated_coagulation_factor_ix) <- r_axi_ix_activation + r_avii_tf_complex_ix_activation - r_activated_ix_degradation - r_aix_aviii_complex_formation - r_aix_atiii_heparin_complex_formation
    d/dt(coagulation_factor_x) <- r_x_production - r_x_degradation - r_avii_x_activation - r_avii_tf_complex_x_activation - r_aix_x_activation - r_aix_aviii_complex_x_activation
    d/dt(activated_coagulation_factor_x) <- r_avii_x_activation + r_avii_tf_complex_x_activation + r_aix_x_activation + r_aix_aviii_complex_x_activation - r_av_ax_complex_formation - r_ax_degradation - r_ax_atiii_heparin_complex_complex_formation - r_ax_tfpi_complex_formation
    d/dt(coagulation_factor_xi) <- r_xi_production - r_xi_degradation - r_thrombin_xi_activation - r_axii_xi_activation
    d/dt(activated_coagulation_factor_xi) <- r_thrombin_xi_activation + r_axii_xi_activation - r_activated_xi_degradation
    d/dt(coagulation_factor_xii) <- r_xii_production - r_xii_degradation - r_ca_xii_activation - r_kallikrein_xii_activation
    d/dt(activated_coagulation_factor_xii) <- r_ca_xii_activation + r_kallikrein_xii_activation - r_activated_xii_degradation
    d/dt(coagulation_factor_xiii) <- r_xiii_production - r_xiii_degradation - r_thrombin_xiii_activation
    d/dt(activated_coagulation_factor_xiii) <- r_thrombin_xiii_activation - r_axiii_degradation
    d/dt(contact_system_activator) <- -r_ca_degradation
    d/dt(fibrin) <- r_thrombin_fibrinogen_conversion - r_fibrin_degradation - r_axiii_fibrin_conversion - r_plasmin_fibrin_degradation
    d/dt(cross_linked_fibrin) <- r_axiii_fibrin_conversion - r_xf_degradation - r_plasmin_xf_degradation - r_apc_ps_complex_xf_degradation
    d/dt(fibrin_degradation_product) <- r_fibrinogen_degradation + r_plasmin_fibrinogen_degradation + r_fibrin_degradation + r_plasmin_fibrin_degradation - r_fdp_degradation
    d/dt(d_dimer) <- r_xf_degradation + r_plasmin_xf_degradation + r_apc_ps_complex_xf_degradation - r_d_dimer_degradation
    d/dt(fibrinogen) <- r_fibrinogen_production - r_fibrinogen_degradation - r_thrombin_fibrinogen_conversion - r_plasmin_fibrinogen_degradation
    d/dt(kallikrein) <- r_axii_prekallikrein_conversion - r_kallikrein_degradation
    d/dt(prekallikrein) <- r_prekallikrein_production - r_prekallikrein_degradation - r_axii_prekallikrein_conversion
    d/dt(plasmin) <- r_thrombin_plasminogen_conversion + r_fibrin_plasminogen_conversion + r_apc_ps_complex_plasminogen_conversion - r_plasmin_degradation
    d/dt(plasminogen) <- r_plasminogen_production - r_plasminogen_degradation - r_thrombin_plasminogen_conversion - r_fibrin_plasminogen_conversion - r_apc_ps_complex_plasminogen_conversion
    d/dt(protein_c) <- r_protein_c_production - r_protein_c_degradation - r_aii_tmod_complex_pc_activation
    d/dt(activated_protein_c) <- r_aii_tmod_complex_pc_activation - r_activated_protein_c_degradation - r_apc_ps_complex_formation
    d/dt(protein_s) <- r_protein_s_production - r_apc_ps_complex_formation - r_protein_s_degradation
    d/dt(thrombomodulin) <- r_thrombomodulin_production - r_aii_tmod_complex_formation - r_thrombomodulin_degradation
    d/dt(tissue_factor) <- -r_vii_tf_complex_formation - r_avii_tf_complex_formation - r_tissue_factor_degradation
    d/dt(tissue_factor_pathway_inhibitor) <- r_tfpi_production - r_ax_tfpi_complex_formation - r_tfpi_degradation
    d/dt(vitamin_k) <- r_vitamin_k_intake + r_vitamin_k_epoxide_conversion + r_vk_transition_to_central - r_vk_degradation - r_vitamin_k_conversion - r_vk_transition_to_peripheral
    d/dt(vitamin_k_epoxide) <- r_vkh2_conversion - r_vitamin_k_epoxide_conversion
    d/dt(vitamin_k_hydroquinone) <- r_vitamin_k_conversion - r_vkh2_conversion
    d/dt(vitamin_k_peripheral) <- v_blood / v_vk_peripheral * (r_vk_transition_to_peripheral - r_vk_transition_to_central)
    d/dt(aii_atiii_heparin_complex) <- r_aii_atiii_heparin_complex_formation
    d/dt(aii_at_complex) <- r_thrombin_degradation - r_tat_degradation
    d/dt(aii_tmod_complex) <- r_aii_tmod_complex_formation - r_aii_tmod_complex_degradation
    d/dt(av_ax_complex) <- r_av_ax_complex_formation - r_av_ax_complex_degradation - r_apc_ps_complex_av_ax_complex_degradation
    d/dt(vii_tf_complex) <- r_vii_tf_complex_formation - r_vii_tf_complex_degradation - r_ax_vii_tf_complex_activation - r_tf_vii_tf_complex_activation
    d/dt(avii_tf_complex) <- r_ax_vii_tf_complex_activation + r_tf_vii_tf_complex_activation + r_avii_tf_complex_formation - r_avii_tf_complex_degradation - r_avii_tf_ax_tfpi_complex_formation
    d/dt(avii_tf_ax_tfpi_complex) <- r_avii_tf_ax_tfpi_complex_formation - r_avii_tf_ax_tfpi_complex_degradation
    d/dt(aix_aviii_complex) <- r_aix_aviii_complex_formation - r_aix_aviii_complex_degradation
    d/dt(aix_atiii_heparin_complex) <- r_aix_atiii_heparin_complex_formation
    d/dt(ax_atiii_heparin_complex_complex) <- r_ax_atiii_heparin_complex_complex_formation
    d/dt(ax_tfpi_complex) <- r_ax_tfpi_complex_formation - r_avii_tf_ax_tfpi_complex_formation - r_ax_tfpi_complex_degradation
    d/dt(apc_ps_complex) <- r_apc_ps_complex_formation - r_apc_ps_complex_degradation
    d/dt(atiii_heparin_complex) <- -r_atiii_heparin_complex_elimination - r_aii_atiii_heparin_complex_formation - r_aix_atiii_heparin_complex_formation - r_ax_atiii_heparin_complex_complex_formation
    # Auxiliary fibrin integral used by the prothrombin-time test (Section S1.1.2)
    d/dt(fibrin_auc) <- fibrin

    coagulation_factor_ii(0) <- bl_coagulation_factor_ii
    coagulation_factor_v(0) <- bl_coagulation_factor_v
    coagulation_factor_vii(0) <- bl_coagulation_factor_vii
    coagulation_factor_viii(0) <- bl_coagulation_factor_viii
    coagulation_factor_ix(0) <- bl_coagulation_factor_ix
    coagulation_factor_x(0) <- bl_coagulation_factor_x
    coagulation_factor_xi(0) <- bl_coagulation_factor_xi
    coagulation_factor_xii(0) <- bl_coagulation_factor_xii
    coagulation_factor_xiii(0) <- bl_coagulation_factor_xiii
    fibrinogen(0) <- bl_fibrinogen
    prekallikrein(0) <- bl_prekallikrein
    plasminogen(0) <- bl_plasminogen
    protein_c(0) <- bl_protein_c
    protein_s(0) <- bl_protein_s
    thrombomodulin(0) <- bl_thrombomodulin
    tissue_factor_pathway_inhibitor(0) <- bl_tissue_factor_pathway_inhibitor
    vitamin_k(0) <- bl_vitamin_k
    vitamin_k_epoxide(0) <- bl_vitamin_k_epoxide
    vitamin_k_hydroquinone(0) <- bl_vitamin_k_hydroquinone
    vitamin_k_peripheral(0) <- bl_vitamin_k_peripheral

    Cc ~ lnorm(expSd)
  })
}
