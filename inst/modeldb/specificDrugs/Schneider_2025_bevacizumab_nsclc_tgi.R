Schneider_2025_bevacizumab_nsclc_tgi <- function() {
  description <- paste(
    "QSP (semi-mechanistic tumor growth inhibition). Joint model of",
    "non-small cell lung cancer lesion growth and response to bevacizumab plus",
    "cisplatin, carboplatin, pemetrexed, paclitaxel, docetaxel, gemcitabine,",
    "erlotinib and drozitumab (apomab), fitted to individual lesion longest",
    "diameters from 11 pooled clinical trials (Schneider 2025). Each drug has",
    "a fixed literature PK model without IIV. The proliferating tumor volume",
    "follows Gompertz growth scaled to its own baseline; cytotoxic drugs move",
    "cells irreversibly into a three-compartment death chain (Q_gamma), and",
    "microtubule agents plus the direct bevacizumab effect move cells into a",
    "reversibly injured pool (Q_rho) from which a fraction repairs. Every drug",
    "effect is a weighted concentration attenuated by acquired resistance",
    "exp(-lambda * AUC). Bevacizumab additionally amplifies both kill terms",
    "through a copy of its PK delayed by tau (transient vascular",
    "normalisation). Each lesion is a separate subject; the observation is the",
    "lesion longest diameter of a sphere of the summed tumor volume.",
    sep = " "
  )
  reference <- paste(
    "Schneider BK, Benzekry S, Mochel JP. Optimizing First-Line Therapeutics",
    "in Non-Small Cell Lung Cancer: Insights From Joint Modeling and",
    "Large-Scale Data Analysis. CPT Pharmacometrics Syst Pharmacol.",
    "2025;14(10):1637-1649. doi:10.1002/psp4.70079.",
    "Structural equations and parameter mapping from the deposited Monolix",
    "Mlxtran model (Supporting Information, PSP-2024-0266-s03.txt);",
    "pharmacodynamic estimates from Table 3; pharmacokinetic values from",
    "Table 2 and the deposited code (PK fixed from Lu 2008 bevacizumab,",
    "Urien 2005 cisplatin, Latz 2006 pemetrexed, Camidge 2010 drozitumab,",
    "Henningsson 2001 paclitaxel, Joerger 2007 carboplatin, Jiang 2008",
    "gemcitabine, Slaviero 2004 docetaxel, Lu 2006 / Frohna 2006 erlotinib).",
    sep = " "
  )
  vignette <- "Schneider_2025_bevacizumab_nsclc_tgi"

  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "mm (the observed variable tumor_size is a lesion longest diameter; drug concentrations Cc_<drug> are in mg/mL)"
  )
  dosing <- c(
    "central_bevacizumab",
    "depot_erlotinib",
    "central_cisplatin",
    "central_pemetrexed",
    "central_drozitumab",
    "central_paclitaxel",
    "central_carboplatin",
    "central_gemcitabine",
    "central_docetaxel",
    "central_bevacizumab_lagged"
  )

  # The tau-delayed copy of the bevacizumab PK drives vascular normalisation
  # and receives the same doses as central_bevacizumab; the reversibly
  # injured tumor pool has no registered TGI-state name.
  paper_specific_compartments <- c(
    "central_bevacizumab_lagged",
    "peripheral1_bevacizumab_lagged",
    "injured_cells"
  )

  covariateData <- list(
    TUM_VOL = list(
      description = "Observed baseline volume of the modelled lesion, computed from its longest diameter D0 (mm) at time 0 as the volume of a sphere, pi / 6 * D0^3.",
      units = "mm^3",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per lesion (each lesion is a separate subject). Sets the initial condition of the proliferating tumor volume (cycling_cells(0) = TUM_VOL) and is the reference volume of the baseline-scaled Gompertz growth term, so the predicted diameter at time 0 equals the observed one. Schneider 2025 fixed initial tumor volumes to the time-0 measurement of each lesion rather than estimating a baseline; the deposited code carries it as the input v0_err. Because every tumor term is proportional to volume and growth is scaled to TUM_VOL, the relative time course is independent of TUM_VOL.",
      source_name = "v0_err"
    )
  )

  compartmentData <- list(
    central_bevacizumab = list(analyte = "bevacizumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_bevacizumab = list(analyte = "bevacizumab", units = "mg", specimen = "tissue", verified = TRUE),
    central_bevacizumab_lagged = list(
      analyte = "bevacizumab (copy delayed by tau, drives vascular normalisation)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_bevacizumab_lagged = list(
      analyte = "bevacizumab (copy delayed by tau, drives vascular normalisation)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    central_cisplatin = list(analyte = "cisplatin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_cisplatin = list(analyte = "cisplatin", units = "mg", specimen = "tissue", verified = TRUE),
    central_pemetrexed = list(analyte = "pemetrexed", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_pemetrexed = list(analyte = "pemetrexed", units = "mg", specimen = "tissue", verified = TRUE),
    central_drozitumab = list(analyte = "drozitumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_drozitumab = list(analyte = "drozitumab", units = "mg", specimen = "tissue", verified = TRUE),
    central_paclitaxel = list(analyte = "paclitaxel", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_paclitaxel = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2_paclitaxel = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    central_carboplatin = list(analyte = "carboplatin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_carboplatin = list(analyte = "carboplatin", units = "mg", specimen = "tissue", verified = TRUE),
    central_gemcitabine = list(analyte = "gemcitabine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_gemcitabine = list(analyte = "gemcitabine", units = "mg", specimen = "tissue", verified = TRUE),
    central_docetaxel = list(analyte = "docetaxel", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_docetaxel = list(analyte = "docetaxel", units = "mg", specimen = "tissue", verified = TRUE),
    depot_erlotinib = list(analyte = "erlotinib", units = "mg", specimen = "administration site", verified = TRUE),
    central_erlotinib = list(analyte = "erlotinib", units = "mg", specimen = "plasma", verified = TRUE),
    auc_central_cisplatin = list(
      analyte = "cisplatin",
      units = "mg*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_central_carboplatin = list(
      analyte = "carboplatin",
      units = "mg*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_central_pemetrexed = list(
      analyte = "pemetrexed",
      units = "mg*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_central_drozitumab = list(
      analyte = "drozitumab",
      units = "mg*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_central_erlotinib = list(
      analyte = "erlotinib",
      units = "mg*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_central_gemcitabine = list(
      analyte = "gemcitabine",
      units = "mg*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_central_paclitaxel = list(
      analyte = "paclitaxel",
      units = "mg*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_central_docetaxel = list(
      analyte = "docetaxel",
      units = "mg*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_central_bevacizumab = list(
      analyte = "bevacizumab",
      units = "mg*day/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    cycling_cells = list(
      analyte = "proliferating (unperturbed) tumor cells",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    injured_cells = list(
      analyte = "reversibly injured tumor cells",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells1 = list(
      analyte = "irreversibly damaged (dying) tumor cells",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells2 = list(
      analyte = "irreversibly damaged (dying) tumor cells",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells3 = list(
      analyte = "irreversibly damaged (dying) tumor cells",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1963L,
    n_studies = 11L,
    disease_state = "non-small cell lung cancer treated with bevacizumab-containing first- and later-line regimens",
    dose_range = "bevacizumab with platinum doublets (cisplatin or carboplatin plus pemetrexed, paclitaxel, gemcitabine or docetaxel), erlotinib or drozitumab (PRO95780, apomab), at the doses of the contributing trials",
    regions = "multinational (Roche and Eli Lilly sponsored trials shared through CSDR / Vivli)",
    notes = "Eleven phase 2-4 trials (Table 1: NCT00800202, AVF0757G, NCT01004250, NCT00312728, NCT00480831, NCT00700180, NCT01364012, NCT00762034, NCT00806923, NCT02596958, NCT00388206). Of 3686 candidate patients, 2586 had complete data (6197 lesions, 29,885 samples). After requiring >= 3 measurements per lesion and removing monotonic non-responders, 4450 lesions from 1963 patients (26,515 samples) remained; about 85% of lesions were used for estimation and 15% reserved for external validation. Each lesion was treated as a separate subject (no distinction between within- and between-patient variability). Demographics were not reported. Vinorelbine appears in one trial but has no effect term in the model."
  )

  ini({
    # ---- Pharmacokinetics: literature models fixed without IIV (Table 2 and
    # the deposited Mlxtran PK block, which gives the values in mL and
    # mL/day; divided by 1000 here to L and L/day). The paper states
    # between-patient PK variability was not included.
    lcl_bevacizumab <- fixed(log(0.216291))
    label("Bevacizumab clearance (L/day)") # Supplement s03 cl_bev = 216.291 mL/day; Table 2 Cl = 0.216 L/day (Lu 2008)
    lvc_bevacizumab <- fixed(log(2.80412))
    label("Bevacizumab central volume (L)") # Supplement s03 v1_bev = 2804.12 mL; Table 2 V1 = 2.80 L
    lk12_bevacizumab <- fixed(log(0.223))
    label("Bevacizumab central-to-peripheral rate constant (1/day)") # Supplement s03 k12_bev = 0.223; Table 2 k12 = 0.223 /day
    lk21_bevacizumab <- fixed(log(0.215))
    label("Bevacizumab peripheral-to-central rate constant (1/day)") # Supplement s03 k21_bev = 0.215; Table 2 k21 = 0.215 /day

    lcl_cisplatin <- fixed(log(6.48))
    label("Cisplatin clearance (L/day)") # Supplement s03 cl_cis = 6480 mL/day; Table 2 Cl = 6.48 L/day (Urien 2005)
    lvc_cisplatin <- fixed(log(22.3))
    label("Cisplatin central volume (L)") # Supplement s03 v1_cis = 22300 mL; Table 2 V1 = 22.3 L
    lq_cisplatin <- fixed(log(456))
    label("Cisplatin intercompartmental clearance (L/day)") # Supplement s03 q_cis = 456000 mL/day; Table 2 Q1 = 456.0 L/day
    lvp_cisplatin <- fixed(log(77))
    label("Cisplatin peripheral volume (L)") # Supplement s03 v2_cis = 77000 mL; Table 2 V2 = 77.0 L

    lcl_pemetrexed <- fixed(log(131.904))
    label("Pemetrexed clearance (L/day)") # Supplement s03 cl_pem = 131904 mL/day; Table 2 Cl = 131.9 L/day (Latz 2006)
    lvc_pemetrexed <- fixed(log(12.9))
    label("Pemetrexed central volume (L)") # Supplement s03 v1_pem = 12900 mL; Table 2 V1 = 12.9 L
    lq_pemetrexed <- fixed(log(20.736))
    label("Pemetrexed intercompartmental clearance (L/day)") # Supplement s03 q_pem = 20736 mL/day; Table 2 Q1 = 20.7 L/day
    lvp_pemetrexed <- fixed(log(3.38))
    label("Pemetrexed peripheral volume (L)") # Supplement s03 v2_pem = 3380 mL; Table 2 V2 = 3.38 L

    lcl_drozitumab <- fixed(log(0.328))
    label("Drozitumab (apomab) clearance (L/day)") # Supplement s03 cl_apo = 328 mL/day; Table 2 Cl = 0.328 L/day (Camidge 2010)
    lvc_drozitumab <- fixed(log(3.97))
    label("Drozitumab (apomab) central volume (L)") # Supplement s03 v1_apo = 3970 mL; Table 2 V1 = 3.97 L
    lq_drozitumab <- fixed(log(0.793))
    label("Drozitumab (apomab) intercompartmental clearance (L/day)") # Supplement s03 q_apo = 793 mL/day; Table 2 Q1 = 0.793 L/day
    lvp_drozitumab <- fixed(log(3.84))
    label("Drozitumab (apomab) peripheral volume (L)") # Supplement s03 v2_apo = 3840 mL; Table 2 V2 = 3.84 L

    lcl_paclitaxel <- fixed(log(10296))
    label("Paclitaxel clearance (L/day)") # Supplement s03 cl_pac = 10296000 mL/day; Table 2 Cl = 10296.0 L/day (Henningsson 2001)
    lvc_paclitaxel <- fixed(log(229))
    label("Paclitaxel central volume (L)") # Supplement s03 v1_pac = 229000 mL; Table 2 V1 = 229.0 L
    lq_paclitaxel <- fixed(log(3216))
    label("Paclitaxel intercompartmental clearance to the first peripheral (L/day)") # Supplement s03 q1_pac = 3216000 mL/day; Table 2 Q1 = 3216.0 L/day
    lvp_paclitaxel <- fixed(log(856))
    label("Paclitaxel first peripheral volume (L)") # Supplement s03 v2_pac = 856000 mL; Table 2 V2 = 856.0 L
    lq2_paclitaxel <- fixed(log(5112))
    label("Paclitaxel intercompartmental clearance to the second peripheral (L/day)") # Supplement s03 q2_pac = 5112000 mL/day; Table 2 Q2 = 5112.0 L/day
    lvp2_paclitaxel <- fixed(log(30.3))
    label("Paclitaxel second peripheral volume (L)") # Supplement s03 v3_pac = 30300 mL; Table 2 V3 = 30.3 L

    lcl_carboplatin <- fixed(log(177.12))
    label("Carboplatin clearance (L/day)") # Supplement s03 cl_car = 177120 mL/day; Table 2 Cl = 177.12 L/day (Joerger 2007)
    lvc_carboplatin <- fixed(log(11.9))
    label("Carboplatin central volume (L)") # Supplement s03 v1_car = 11900 mL; Table 2 V1 = 11.9 L
    lq_carboplatin <- fixed(log(2172))
    label("Carboplatin intercompartmental clearance (L/day)") # Supplement s03 q_car = 2172000 mL/day; Table 2 Q1 = 2172.0 L/day
    lvp_carboplatin <- fixed(log(8.23))
    label("Carboplatin peripheral volume (L)") # Supplement s03 v2_car = 8230 mL; Table 2 V2 = 8.23 L

    lcl_gemcitabine <- fixed(log(3888))
    label("Gemcitabine clearance (L/day)") # Supplement s03 cl_gem = 3888000 mL/day; Table 2 Cl = 3888.0 L/day (Jiang 2008)
    lvc_gemcitabine <- fixed(log(15))
    label("Gemcitabine central volume (L)") # Supplement s03 v1_gem = 15000 mL; Table 2 V1 = 15.0 L
    lq_gemcitabine <- fixed(log(1008))
    label("Gemcitabine intercompartmental clearance (L/day)") # Supplement s03 q_gem = 1008000 mL/day; Table 2 Q1 = 1008.0 L/day
    lvp_gemcitabine <- fixed(log(15))
    label("Gemcitabine peripheral volume (L)") # Supplement s03 v2_gem = 15000 mL; Table 2 V2 = 15.0 L

    lcl_docetaxel <- fixed(log(723.12))
    label("Docetaxel clearance (L/day)") # Supplement s03 cl_doc = 723120 mL/day; Table 2 Cl = 723.12 L/day (Slaviero 2004)
    lvc_docetaxel <- fixed(log(7.9))
    label("Docetaxel central volume (L)") # Supplement s03 v1_doc = 7900 mL; Table 2 V1 = 7.9 L
    lk12_docetaxel <- fixed(log(27.12))
    label("Docetaxel central-to-peripheral rate constant (1/day)") # Supplement s03 k12_doc = 27.12; Table 2 k12 = 27.12 /day
    lk21_docetaxel <- fixed(log(3.6))
    label("Docetaxel peripheral-to-central rate constant (1/day)") # Supplement s03 k21_doc = 3.6; Table 2 k21 = 3.6 /day

    lka_erlotinib <- fixed(log(21.36))
    label("Erlotinib absorption rate constant (1/day)") # Supplement s03 ka_erl = 21.36; Table 2 footnote ka = 21.36 /day (Lu 2006, Frohna 2006)
    lcl_erlotinib <- fixed(log(102.96))
    label("Erlotinib clearance (L/day)") # Supplement s03 cl_erl = 102960 mL/day; Table 2 Cl = 102.96 L/day
    lvc_erlotinib <- fixed(log(210))
    label("Erlotinib volume (L)") # Supplement s03 v1_erl = 210000 mL; Table 2 V1 = 210.0 L
    lfdepot_erlotinib <- fixed(log(0.6))
    label("Erlotinib oral bioavailability (fraction)") # Supplement s03 F_erl = .6; Table 2 footnote bioavailability = 60%

    # ---- Tumor growth (Table 3; Equation 2)
    lalpha <- log(0.00304)
    label("Gompertz growth rate at baseline volume alpha (1/day)") # Table 3 alpha = 0.00304 day-1 (RSE 2.72%)
    lbeta <- log(0.0398)
    label("Gompertz growth deceleration beta (1/day)") # Table 3 beta = 0.0398 day-1 (RSE 1.37%)
    ltau <- log(0.254)
    label("Delay between bevacizumab administration and vascular normalisation tau (day)") # Table 3 tau = 0.254 day (RSE 0.723%)
    lkk <- log(28)
    label("Transit rate constant through the irreversible death chain kk (1/day)") # Table 3 kk = 28 day-1 (RSE 7.3%)

    # ---- Drug-effect weights w_d (Equation 4; slope per mg/mL, effects in
    # 1/day). Row order of Table 3 matches the input order of the deposited
    # code, which maps each row to its role in the effect equations.
    lw_cisplatin_rel <- log(0.917)
    label("Platinum weight w_cis_rel: cisplatin potency, also scales carboplatin (mL/mg/day)") # Table 3 'w cisplatin relative to carboplatin' = 0.917 (RSE 7.94%); code w_cis_rel
    lw_carboplatin <- log(0.257)
    label("Carboplatin potency relative to cisplatin within the platinum effect (unitless)") # Table 3 'w carboplatin' = 0.257 (RSE 1.89%); code w_car
    lw_pemetrexed <- log(23.8)
    label("Pemetrexed (antifolate) potency weight (mL/mg/day)") # Table 3 'w pemetrexed' = 23.8 (RSE 3.19%); code w_pem
    lw_drozitumab <- log(0.75)
    label("Drozitumab (DR5 agonist) potency weight (mL/mg/day)") # Table 3 'w apomab' = 0.75 (RSE 0.456%); code w_apo
    lw_erlotinib <- log(0.62)
    label("Erlotinib (EGFR) potency weight (mL/mg/day)") # Table 3 'w erlotinib' = 0.62 (RSE 0.968%); code w_erl
    lw_gemcitabine <- log(1)
    label("Gemcitabine (nucleoside analog) potency weight (mL/mg/day)") # Table 3 'w gemcitabine' = 1 (RSE 1.45%); code w_gem
    lw_paclitaxel_rel <- log(0.758)
    label("Microtubule weight w_pac_rel: paclitaxel potency, also scales docetaxel (mL/mg/day)") # Table 3 'w paclitaxel relative to docetaxel' = 0.758 (RSE 1.26%); code w_pac_rel
    lw_docetaxel <- log(0.939)
    label("Docetaxel potency relative to paclitaxel within the microtubule effect (unitless)") # Table 3 'w docetaxel' = 0.939 (RSE 0.248%); code w_doc
    lw_bevacizumab <- log(4.14)
    label("Direct (anti-VEGF) bevacizumab potency weight on reversible injury (mL/mg/day)") # Table 3 first 'w bevacizumab,rho' row = 4.14 (RSE 6.82%); code w_bev

    # ---- Acquired resistance lambda_d (Equation 6; per unit AUC, mL/mg/day)
    lres_cisplatin_rel <- log(0.00101)
    label("Platinum resistance rate lambda_cis_rel per unit platinum AUC (mL/mg/day)") # Table 3 'lambda cisplatin relative to carboplatin' = 0.00101 (RSE 0.907%); code L_cis_rel
    lres_carboplatin <- log(0.000995)
    label("Carboplatin AUC weight relative to cisplatin AUC in platinum resistance (unitless)") # Table 3 'lambda carboplatin' = 0.000995 (RSE 0.963%); code L_car
    lres_pemetrexed <- log(0.000796)
    label("Pemetrexed resistance rate per unit AUC (mL/mg/day)") # Table 3 'lambda pemetrexed' = 0.000796 (RSE 0.911%); code L_pem
    lres_drozitumab <- log(0.000656)
    label("Drozitumab resistance rate per unit AUC (mL/mg/day)") # Table 3 'lambda apomab' = 0.000656 (RSE 0.206%); code L_apo
    lres_erlotinib <- log(0.00119)
    label("Erlotinib resistance rate per unit AUC (mL/mg/day)") # Table 3 'lambda erlotinib' = 0.00119 (RSE 1.81%); code L_erl
    lres_gemcitabine <- log(0.00152)
    label("Gemcitabine resistance rate per unit AUC (mL/mg/day)") # Table 3 'lambda gemcitabine' = 0.00152 (RSE 1.09%); code L_gem
    lres_paclitaxel_rel <- log(0.00111)
    label("Microtubule resistance rate lambda_pac_rel per unit microtubule-agent AUC (mL/mg/day)") # Table 3 'lambda paclitaxel relative to docetaxel' = 0.00111 (RSE 0.918%); code L_pac_rel
    lres_docetaxel <- log(0.00134)
    label("Docetaxel AUC weight relative to paclitaxel AUC in microtubule resistance (unitless)") # Table 3 'lambda docetaxel' = 0.00134 (RSE 0.776%); code L_doc
    lres_bevacizumab <- log(0.0164)
    label("Bevacizumab direct-effect resistance rate per unit AUC (mL/mg/day)") # Table 3 'lambda bevacizumab,rho' = 0.0164 (RSE 3.79%); code L_bev

    # ---- Bevacizumab vascular-normalisation amplification (Equations 1-2)
    lw_bevacizumab_gamma <- log(0.727)
    label("Amplification of irreversible kill Q_gamma per delayed bevacizumab concentration (mL/mg)") # Table 3 'w bevacizumab,gamma' = 0.727 (RSE 0.878%); code w_bev_gamma
    lw_bevacizumab_rho <- log(1.39)
    label("Amplification of reversible injury Q_rho per delayed bevacizumab concentration (mL/mg)") # Table 3 second 'w bevacizumab,rho' row = 1.39 (RSE 2.23%); code w_bev_rho

    # ---- Reversible injury (Equation 2)
    lkinj <- log(0.0746)
    label("Exit rate constant from the reversibly injured pool kkinj (1/day)") # Table 3 kkinj = 0.0746 (RSE 4.89%); code kk_i
    logitfrepair <- logit(0.364)
    label("Fraction of injured cells that repair and return to proliferation p (fraction)") # Table 3 p = 0.364 (RSE 1.98%); code imv_r_perc; logit-normal

    # ---- IIV (Table 3 'IIV' column = SD of the random effect; variance = SD^2).
    # Log-normal for every parameter except p, which is logit-normal. No
    # correlations were estimated.
    etalalpha ~ 1.22^2 # Table 3 IIV alpha 1.22 (SD of the random effect)
    etalbeta ~ 0.62^2 # Table 3 IIV beta 0.62 (SD of the random effect)
    etaltau ~ 0.304^2 # Table 3 IIV tau 0.304 (SD of the random effect)
    etalkk ~ 3.41^2 # Table 3 IIV kk 3.41 (SD of the random effect)
    etalw_cisplatin_rel ~ 3.59^2 # Table 3 IIV w cisplatin rel 3.59 (SD of the random effect)
    etalw_carboplatin ~ 0.808^2 # Table 3 IIV w carboplatin 0.808 (SD of the random effect)
    etalw_pemetrexed ~ 1.41^2 # Table 3 IIV w pemetrexed 1.41 (SD of the random effect)
    etalw_drozitumab ~ 0.194^2 # Table 3 IIV w apomab 0.194 (SD of the random effect)
    etalw_erlotinib ~ 0.407^2 # Table 3 IIV w erlotinib 0.407 (SD of the random effect)
    etalw_gemcitabine ~ 0.616^2 # Table 3 IIV w gemcitabine 0.616 (SD of the random effect)
    etalw_paclitaxel_rel ~ 0.537^2 # Table 3 IIV w paclitaxel rel 0.537 (SD of the random effect)
    etalw_docetaxel ~ 0.105^2 # Table 3 IIV w docetaxel 0.105 (SD of the random effect)
    etalw_bevacizumab ~ 3.44^2 # Table 3 IIV first w bevacizumab,rho row 3.44 (SD of the random effect)
    etalres_cisplatin_rel ~ 0.386^2 # Table 3 IIV lambda cisplatin rel 0.386 (SD of the random effect)
    etalres_carboplatin ~ 0.409^2 # Table 3 IIV lambda carboplatin 0.409 (SD of the random effect)
    etalres_pemetrexed ~ 0.385^2 # Table 3 IIV lambda pemetrexed 0.385 (SD of the random effect)
    etalres_drozitumab ~ 0.0877^2 # Table 3 IIV lambda apomab 0.0877 (SD of the random effect)
    etalres_erlotinib ~ 0.777^2 # Table 3 IIV lambda erlotinib 0.777 (SD of the random effect)
    etalres_gemcitabine ~ 0.47^2 # Table 3 IIV lambda gemcitabine 0.47 (SD of the random effect)
    etalres_paclitaxel_rel ~ 0.386^2 # Table 3 IIV lambda paclitaxel rel 0.386 (SD of the random effect)
    etalres_docetaxel ~ 0.327^2 # Table 3 IIV lambda docetaxel 0.327 (SD of the random effect)
    etalres_bevacizumab ~ 1.69^2 # Table 3 IIV lambda bevacizumab,rho 1.69 (SD of the random effect)
    etalw_bevacizumab_gamma ~ 0.373^2 # Table 3 IIV w bevacizumab,gamma 0.373 (SD of the random effect)
    etalw_bevacizumab_rho ~ 0.957^2 # Table 3 IIV second w bevacizumab,rho row 0.957 (SD of the random effect)
    etalkinj ~ 2.31^2 # Table 3 IIV kkinj 2.31 (SD of the random effect)
    etalogitfrepair ~ 1.37^2 # Table 3 IIV p 1.37 (SD of the random effect) (logit scale)

    # ---- Residual error: Monolix 'combined1', y = f + (a + b * f) * e
    addSd <- 0.0457
    label("Additive residual SD a (mm)") # Table 3 a = 0.0457 (RSE 2.53%)
    propSd <- 0.105
    label("Proportional residual SD b (fraction)") # Table 3 b = 0.105 (RSE 1%)
  })
  model({
    # 1. Pharmacokinetic parameters (L, L/day, 1/day)
    cl_bevacizumab <- exp(lcl_bevacizumab)
    vc_bevacizumab <- exp(lvc_bevacizumab)
    k12_bevacizumab <- exp(lk12_bevacizumab)
    k21_bevacizumab <- exp(lk21_bevacizumab)
    cl_cisplatin <- exp(lcl_cisplatin)
    vc_cisplatin <- exp(lvc_cisplatin)
    q_cisplatin <- exp(lq_cisplatin)
    vp_cisplatin <- exp(lvp_cisplatin)
    cl_pemetrexed <- exp(lcl_pemetrexed)
    vc_pemetrexed <- exp(lvc_pemetrexed)
    q_pemetrexed <- exp(lq_pemetrexed)
    vp_pemetrexed <- exp(lvp_pemetrexed)
    cl_drozitumab <- exp(lcl_drozitumab)
    vc_drozitumab <- exp(lvc_drozitumab)
    q_drozitumab <- exp(lq_drozitumab)
    vp_drozitumab <- exp(lvp_drozitumab)
    cl_paclitaxel <- exp(lcl_paclitaxel)
    vc_paclitaxel <- exp(lvc_paclitaxel)
    q_paclitaxel <- exp(lq_paclitaxel)
    vp_paclitaxel <- exp(lvp_paclitaxel)
    q2_paclitaxel <- exp(lq2_paclitaxel)
    vp2_paclitaxel <- exp(lvp2_paclitaxel)
    cl_carboplatin <- exp(lcl_carboplatin)
    vc_carboplatin <- exp(lvc_carboplatin)
    q_carboplatin <- exp(lq_carboplatin)
    vp_carboplatin <- exp(lvp_carboplatin)
    cl_gemcitabine <- exp(lcl_gemcitabine)
    vc_gemcitabine <- exp(lvc_gemcitabine)
    q_gemcitabine <- exp(lq_gemcitabine)
    vp_gemcitabine <- exp(lvp_gemcitabine)
    cl_docetaxel <- exp(lcl_docetaxel)
    vc_docetaxel <- exp(lvc_docetaxel)
    k12_docetaxel <- exp(lk12_docetaxel)
    k21_docetaxel <- exp(lk21_docetaxel)
    ka_erlotinib <- exp(lka_erlotinib)
    cl_erlotinib <- exp(lcl_erlotinib)
    vc_erlotinib <- exp(lvc_erlotinib)
    fdepot_erlotinib <- exp(lfdepot_erlotinib)

    # 2. Individual pharmacodynamic parameters
    alpha <- exp(lalpha + etalalpha)
    beta <- exp(lbeta + etalbeta)
    tau <- exp(ltau + etaltau)
    kk <- exp(lkk + etalkk)
    w_cisplatin_rel <- exp(lw_cisplatin_rel + etalw_cisplatin_rel)
    w_carboplatin <- exp(lw_carboplatin + etalw_carboplatin)
    w_pemetrexed <- exp(lw_pemetrexed + etalw_pemetrexed)
    w_drozitumab <- exp(lw_drozitumab + etalw_drozitumab)
    w_erlotinib <- exp(lw_erlotinib + etalw_erlotinib)
    w_gemcitabine <- exp(lw_gemcitabine + etalw_gemcitabine)
    w_paclitaxel_rel <- exp(lw_paclitaxel_rel + etalw_paclitaxel_rel)
    w_docetaxel <- exp(lw_docetaxel + etalw_docetaxel)
    w_bevacizumab <- exp(lw_bevacizumab + etalw_bevacizumab)
    res_cisplatin_rel <- exp(lres_cisplatin_rel + etalres_cisplatin_rel)
    res_carboplatin <- exp(lres_carboplatin + etalres_carboplatin)
    res_pemetrexed <- exp(lres_pemetrexed + etalres_pemetrexed)
    res_drozitumab <- exp(lres_drozitumab + etalres_drozitumab)
    res_erlotinib <- exp(lres_erlotinib + etalres_erlotinib)
    res_gemcitabine <- exp(lres_gemcitabine + etalres_gemcitabine)
    res_paclitaxel_rel <- exp(lres_paclitaxel_rel + etalres_paclitaxel_rel)
    res_docetaxel <- exp(lres_docetaxel + etalres_docetaxel)
    res_bevacizumab <- exp(lres_bevacizumab + etalres_bevacizumab)
    w_bevacizumab_gamma <- exp(lw_bevacizumab_gamma + etalw_bevacizumab_gamma)
    w_bevacizumab_rho <- exp(lw_bevacizumab_rho + etalw_bevacizumab_rho)
    kinj <- exp(lkinj + etalkinj)
    frepair <- expit(logitfrepair + etalogitfrepair)

    # 3. Drug concentrations in mg/mL (amount in mg over volume in L * 1000),
    #    the concentration scale of the deposited code (volumes in mL).
    Cc_bevacizumab <- central_bevacizumab / (vc_bevacizumab * 1000)
    Cc_bevacizumab_lagged <- central_bevacizumab_lagged / (vc_bevacizumab * 1000)
    Cc_cisplatin <- central_cisplatin / (vc_cisplatin * 1000)
    Cc_pemetrexed <- central_pemetrexed / (vc_pemetrexed * 1000)
    Cc_drozitumab <- central_drozitumab / (vc_drozitumab * 1000)
    Cc_paclitaxel <- central_paclitaxel / (vc_paclitaxel * 1000)
    Cc_carboplatin <- central_carboplatin / (vc_carboplatin * 1000)
    Cc_gemcitabine <- central_gemcitabine / (vc_gemcitabine * 1000)
    Cc_docetaxel <- central_docetaxel / (vc_docetaxel * 1000)
    Cc_erlotinib <- central_erlotinib / (vc_erlotinib * 1000)

    # 4. Pharmacokinetic ODEs. Every bevacizumab dose must be given to both
    #    central_bevacizumab (direct effect, AUC) and
    #    central_bevacizumab_lagged (same kinetics, delayed by tau).
    d/dt(central_bevacizumab) <- -(cl_bevacizumab / vc_bevacizumab + k12_bevacizumab) * central_bevacizumab + k21_bevacizumab * peripheral1_bevacizumab
    d/dt(peripheral1_bevacizumab) <- k12_bevacizumab * central_bevacizumab - k21_bevacizumab * peripheral1_bevacizumab
    d/dt(central_bevacizumab_lagged) <- -(cl_bevacizumab / vc_bevacizumab + k12_bevacizumab) * central_bevacizumab_lagged + k21_bevacizumab * peripheral1_bevacizumab_lagged
    d/dt(peripheral1_bevacizumab_lagged) <- k12_bevacizumab * central_bevacizumab_lagged - k21_bevacizumab * peripheral1_bevacizumab_lagged
    d/dt(central_cisplatin) <- -(cl_cisplatin + q_cisplatin) / vc_cisplatin * central_cisplatin + q_cisplatin / vp_cisplatin * peripheral1_cisplatin
    d/dt(peripheral1_cisplatin) <- q_cisplatin / vc_cisplatin * central_cisplatin - q_cisplatin / vp_cisplatin * peripheral1_cisplatin
    d/dt(central_pemetrexed) <- -(cl_pemetrexed + q_pemetrexed) / vc_pemetrexed * central_pemetrexed + q_pemetrexed / vp_pemetrexed * peripheral1_pemetrexed
    d/dt(peripheral1_pemetrexed) <- q_pemetrexed / vc_pemetrexed * central_pemetrexed - q_pemetrexed / vp_pemetrexed * peripheral1_pemetrexed
    d/dt(central_drozitumab) <- -(cl_drozitumab + q_drozitumab) / vc_drozitumab * central_drozitumab + q_drozitumab / vp_drozitumab * peripheral1_drozitumab
    d/dt(peripheral1_drozitumab) <- q_drozitumab / vc_drozitumab * central_drozitumab - q_drozitumab / vp_drozitumab * peripheral1_drozitumab
    d/dt(central_paclitaxel) <- -(cl_paclitaxel + q_paclitaxel + q2_paclitaxel) / vc_paclitaxel * central_paclitaxel + q_paclitaxel / vp_paclitaxel * peripheral1_paclitaxel + q2_paclitaxel / vp2_paclitaxel * peripheral2_paclitaxel
    d/dt(peripheral1_paclitaxel) <- q_paclitaxel / vc_paclitaxel * central_paclitaxel - q_paclitaxel / vp_paclitaxel * peripheral1_paclitaxel
    d/dt(peripheral2_paclitaxel) <- q2_paclitaxel / vc_paclitaxel * central_paclitaxel - q2_paclitaxel / vp2_paclitaxel * peripheral2_paclitaxel
    d/dt(central_carboplatin) <- -(cl_carboplatin + q_carboplatin) / vc_carboplatin * central_carboplatin + q_carboplatin / vp_carboplatin * peripheral1_carboplatin
    d/dt(peripheral1_carboplatin) <- q_carboplatin / vc_carboplatin * central_carboplatin - q_carboplatin / vp_carboplatin * peripheral1_carboplatin
    d/dt(central_gemcitabine) <- -(cl_gemcitabine + q_gemcitabine) / vc_gemcitabine * central_gemcitabine + q_gemcitabine / vp_gemcitabine * peripheral1_gemcitabine
    d/dt(peripheral1_gemcitabine) <- q_gemcitabine / vc_gemcitabine * central_gemcitabine - q_gemcitabine / vp_gemcitabine * peripheral1_gemcitabine
    d/dt(central_docetaxel) <- -(cl_docetaxel / vc_docetaxel + k12_docetaxel) * central_docetaxel + k21_docetaxel * peripheral1_docetaxel
    d/dt(peripheral1_docetaxel) <- k12_docetaxel * central_docetaxel - k21_docetaxel * peripheral1_docetaxel
    d/dt(depot_erlotinib) <- -ka_erlotinib * depot_erlotinib
    d/dt(central_erlotinib) <- ka_erlotinib * depot_erlotinib - cl_erlotinib / vc_erlotinib * central_erlotinib

    # 5. Cumulative exposure (mg*day/mL) driving acquired resistance
    d/dt(auc_central_cisplatin) <- Cc_cisplatin
    d/dt(auc_central_carboplatin) <- Cc_carboplatin
    d/dt(auc_central_pemetrexed) <- Cc_pemetrexed
    d/dt(auc_central_drozitumab) <- Cc_drozitumab
    d/dt(auc_central_erlotinib) <- Cc_erlotinib
    d/dt(auc_central_gemcitabine) <- Cc_gemcitabine
    d/dt(auc_central_paclitaxel) <- Cc_paclitaxel
    d/dt(auc_central_docetaxel) <- Cc_docetaxel
    d/dt(auc_central_bevacizumab) <- Cc_bevacizumab

    # 6. Mechanism-of-action effects with resistance (1/day; Equations 4 and 6
    #    as written in the deposited code)
    eff_platinum <- w_cisplatin_rel * (Cc_cisplatin + w_carboplatin * Cc_carboplatin) * exp(-res_cisplatin_rel * (auc_central_cisplatin + res_carboplatin * auc_central_carboplatin))
    eff_antifolate <- w_pemetrexed * Cc_pemetrexed * exp(-res_pemetrexed * auc_central_pemetrexed)
    eff_dr5 <- w_drozitumab * Cc_drozitumab * exp(-res_drozitumab * auc_central_drozitumab)
    eff_egfr <- w_erlotinib * Cc_erlotinib * exp(-res_erlotinib * auc_central_erlotinib)
    eff_nucleoside <- w_gemcitabine * Cc_gemcitabine * exp(-res_gemcitabine * auc_central_gemcitabine)
    eff_microtubule <- w_paclitaxel_rel * (Cc_paclitaxel + w_docetaxel * Cc_docetaxel) * exp(-res_paclitaxel_rel * (auc_central_paclitaxel + res_docetaxel * auc_central_docetaxel))
    eff_antivegf <- w_bevacizumab * Cc_bevacizumab * exp(-res_bevacizumab * auc_central_bevacizumab)

    # Vascular normalisation by the tau-delayed bevacizumab amplifies both
    # the irreversible kill (Q_gamma) and the reversible injury (Q_rho).
    q_gamma <- (1 + w_bevacizumab_gamma * Cc_bevacizumab_lagged) * (eff_platinum + eff_antifolate + eff_dr5 + eff_egfr + eff_nucleoside)
    q_rho <- (1 + w_bevacizumab_rho * Cc_bevacizumab_lagged) * (eff_microtubule + eff_antivegf)

    # 7. Tumor dynamics (mm^3; Equations 2-3). Growth is Gompertzian scaled to
    #    the lesion's own baseline volume: the code's
    #    alpha + beta*log(v0/Vc) - beta*log(v/Vc) equals alpha - beta*log(v/v0).
    #    The positivity guards reproduce the if/else branches of the code.
    v0 <- max(TUM_VOL, 1e-50)
    v_pos <- max(0, cycling_cells)
    inj_pos <- max(0, injured_cells)
    if (cycling_cells > 0) {
      dv <- (alpha - beta * log(cycling_cells / v0)) * cycling_cells - q_gamma * cycling_cells - q_rho * cycling_cells + kinj * inj_pos * frepair
    } else {
      dv <- kinj * inj_pos * frepair
    }
    if (injured_cells > 0) {
      dinj <- q_rho * v_pos - kinj * injured_cells - q_gamma * injured_cells
    } else {
      dinj <- q_rho * v_pos
    }
    if (damaged_cells1 > 0) {
      dz1 <- q_gamma * v_pos + q_gamma * inj_pos + kinj * inj_pos * (1 - frepair) - kk * damaged_cells1
    } else {
      dz1 <- q_gamma * v_pos + q_gamma * inj_pos + kinj * inj_pos * (1 - frepair)
    }
    if (damaged_cells2 > 0) {
      dz2 <- kk * damaged_cells1 - kk * damaged_cells2
    } else {
      dz2 <- kk * max(damaged_cells1, 0)
    }
    if (damaged_cells3 > 0) {
      dz3 <- kk * damaged_cells2 - kk * damaged_cells3
    } else {
      dz3 <- kk * max(damaged_cells2, 0)
    }
    d/dt(cycling_cells) <- dv
    d/dt(injured_cells) <- dinj
    d/dt(damaged_cells1) <- dz1
    d/dt(damaged_cells2) <- dz2
    d/dt(damaged_cells3) <- dz3
    cycling_cells(0) <- TUM_VOL

    # 8. Bioavailability and the tau delay of the vascular-normalisation copy
    f(depot_erlotinib) <- fdepot_erlotinib
    alag(central_bevacizumab_lagged) <- tau

    # 9. Observation: longest diameter (mm) of a sphere holding the total
    #    tumor volume, bounded to [0, 1e9] as in the code. A slightly negative
    #    solver volume after eradication has a negative real cube root, which
    #    the code's max(D, 0) maps to 0; clamping the volume first gives the
    #    same 0 without a NaN from a fractional power of a negative number.
    tumor_volume <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3 + injured_cells
    tumor_size <- min(2 * (max(tumor_volume, 0) * 3 / (4 * pi))^(1 / 3), 1e9)
    tumor_size ~ add(addSd) + prop(propSd) + combined1()
  })
}
