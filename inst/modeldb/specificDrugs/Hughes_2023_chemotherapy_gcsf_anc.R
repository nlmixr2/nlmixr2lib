Hughes_2023_chemotherapy_gcsf_anc <- function() {
  description <- paste(
    "Semi-mechanistic Friberg-type kinetic-pharmacodynamic (K-PD) model of",
    "absolute neutrophil count (ANC) during the first three cycles of",
    "cytotoxic chemotherapy with or without subcutaneous G-CSF (filgrastim or",
    "pegfilgrastim) in adults treated for metastatic breast cancer, thoracic",
    "cancer or lymphoma at a single US cancer centre (real-world EHR data).",
    "Chemotherapy is a BSA-scaled, population-average-normalized dose into a",
    "one-compartment K-PD state with fixed CL = 5 L/h and V = 50 L; its",
    "concentration kills proliferating cells linearly (SLOPE). G-CSF is a",
    "fixed one-compartment SC model (ka, CL and F taken from Melhem 2018, V",
    "shared with the chemotherapy state) whose central amount additively",
    "stimulates proliferation (ST1) and multiplies the bone-marrow transit",
    "rate (ST2). Neutrophils follow prol -> transit1-3 -> circ with",
    "asymmetric feedback (CIRC0/circ)^gamma only below baseline. Covariates:",
    "baseline ANC (sets CIRC0), age on SLOPE and ST1, serum albumin on MTT,",
    "chemotherapy-regimen indicators on SLOPE and MTT, and pegfilgrastim vs",
    "filgrastim on the G-CSF PK. This is the PKPD layer of the hybrid",
    "PKPD/machine-learning model; the XGBoost layer is not included."
  )
  reference <- paste(
    "Hughes JH, Woo KH, Keizer RJ, Goswami S. (2023). Clinical decision",
    "support for chemotherapy-induced neutropenia using a hybrid",
    "pharmacodynamic/machine learning model. CPT Pharmacometrics Syst",
    "Pharmacol 12(11):1764-1776. doi:10.1002/psp4.13019."
  )
  vignette <- "Hughes_2023_chemotherapy_gcsf_anc"
  # CIRC0 is the observed baseline ANC times exp(etalcirc0); there is no
  # estimated typical CIRC0, so the eta has no lcirc0 partner.
  paper_specific_etas <- "etalcirc0"

  units <- list(
    time = "h",
    dosing = paste(
      "normalized dose (unitless): each chemotherapy and G-CSF dose was scaled",
      "to body surface area and divided by the training-population average",
      "dose of that drug, so a population-average dose has amount 1;",
      "chemotherapy doses within 24 h are combined into one dose whose amount",
      "is the mean of their normalized amounts"
    ),
    concentration = "ANC in 10^9/L (k/uL)",
    notes = paste(
      "Chemotherapy doses go into 'depot_kpd' (bolus) and G-CSF doses into",
      "'depot' (SC). The chemotherapy concentration driving SLOPE is",
      "depot_kpd / V in normalized-dose units per L. The G-CSF effects ST1 and",
      "ST2 act on the G-CSF central AMOUNT (normalized-dose units), not on a",
      "concentration, as in the published NONMEM code."
    )
  )

  covariateData <- list(
    NEUT = list(
      description = "Baseline absolute neutrophil count, used to set the per-subject baseline CIRC0",
      units = "10^9/L (k/uL)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Observed pre-treatment ANC (NONMEM column CIRC0LIN). Values below",
        "0.5 x 10^9/L are raised to 0.5 before use (Supplement 2 $PK:",
        "'IF(CIRCBL.LT.0.5) CIRCBL = 0.5'). CIRC0 = max(NEUT, 0.5) * exp(eta),",
        "the Dansirikul baseline method (Methods 'PKPD model development').",
        "Supply in 10^9/L, the paper's ANC unit, not the register's",
        "cells/mm^3 default."
      ),
      source_name = "CIRC0LIN"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects (AGE / 57.2)^theta on SLOPE and on ST1 (Table 2; Supplement 2 $PK).",
      source_name = "AGE"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect (ALB / 4.1 g/dL)^0.326 on MTT (Table 2). The paper",
        "reports albumin in g/dL; the model converts the canonical g/L column",
        "to g/dL internally (divide by 10)."
      ),
      source_name = "ALB"
    ),
    FORM_GCSF_PEG = list(
      description = "G-CSF formulation indicator (1 = pegfilgrastim, 0 = filgrastim)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (filgrastim, or no G-CSF)",
      notes = paste(
        "NONMEM PEG = 1 when PEGF = 1 or PEGF300 = 1 (Supplement 2 $PK).",
        "Multiplies the G-CSF absorption rate by 0.1528455, clearance by 0.258",
        "and bioavailability by 0.646 (Table 2, values from Melhem 2018).",
        "Has no effect on a subject who receives no G-CSF."
      ),
      source_name = "PEG"
    ),
    TRT_CYCLOPHOSPHAMIDE_DOXORUBICIN = list(
      description = "Doxorubicin + cyclophosphamide (AC) regimen indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all regimen indicators 0 = reference regimen group)",
      notes = paste(
        "1 for AC alone or AC with rituximab, hyaluronidase-rituximab,",
        "obinutuzumab or pembrolizumab (Supplement 2 $PK REG codes 3, 301,",
        "367, 368, 101). Multiplies SLOPE by 1.68 and MTT by 0.458."
      ),
      source_name = "REG"
    ),
    TRT_CARBOPLATIN_ETOPOSIDE = list(
      description = "Carboplatin + etoposide regimen indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all regimen indicators 0 = reference regimen group)",
      notes = paste(
        "1 for carboplatin + etoposide alone or with atezolizumab, durvalumab,",
        "nivolumab or pembrolizumab (REG 57, 292, 35, 453, 191). Multiplies",
        "SLOPE by 1.44."
      ),
      source_name = "REG"
    ),
    TRT_CARBOPLATIN_PACLITAXEL = list(
      description = "Carboplatin + paclitaxel regimen indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all regimen indicators 0 = reference regimen group)",
      notes = paste(
        "1 for carboplatin + paclitaxel alone or with any of atezolizumab,",
        "bevacizumab, cetuximab, ipilimumab, nivolumab, pembrolizumab,",
        "pertuzumab or trastuzumab (REG 108, 350, 316, 216, 408, 291, 361,",
        "300, 591, 20, 63, 71, 76, 223, 293, 412). Multiplies SLOPE by 0.683."
      ),
      source_name = "REG"
    ),
    TRT_CARBOPLATIN_PEMETREXED = list(
      description = "Carboplatin + pemetrexed regimen indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all regimen indicators 0 = reference regimen group)",
      notes = paste(
        "1 for carboplatin + pemetrexed alone or with bevacizumab,",
        "pembrolizumab or trastuzumab (REG 502, 27, 26, 14, 362). Multiplies",
        "SLOPE by 1.41 and MTT by 1.71."
      ),
      source_name = "REG"
    ),
    TRT_CISPLATIN_PEMETREXED = list(
      description = "Cisplatin + pemetrexed regimen indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all regimen indicators 0 = reference regimen group)",
      notes = paste(
        "1 for cisplatin + pemetrexed alone or with bevacizumab, pembrolizumab",
        "or trastuzumab (REG 52, 82, 28, 415). Multiplies SLOPE by 1.74 and",
        "MTT by 1.89."
      ),
      source_name = "REG"
    ),
    TRT_CYCLOPHOSPHAMIDE_FLUOROURACIL_METHOTREXATE = list(
      description = "Cyclophosphamide + methotrexate + fluorouracil (CMF) regimen indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all regimen indicators 0 = reference regimen group)",
      notes = paste(
        "1 for CMF alone or with trastuzumab with or without pertuzumab (REG",
        "8, 154, 94). Multiplies SLOPE by 0.392 and MTT by 1.57."
      ),
      source_name = "REG"
    ),
    TRT_DOCETAXEL = list(
      description = "Docetaxel (without another cytotoxic agent) regimen indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all regimen indicators 0 = reference regimen group)",
      notes = paste(
        "1 for docetaxel alone or with atezolizumab, bevacizumab,",
        "pembrolizumab, pertuzumab + trastuzumab, ramucirumab or trastuzumab",
        "(REG 384, 188, 25, 420, 119, 29, 421). Multiplies SLOPE by 0.269."
      ),
      source_name = "REG"
    ),
    TRT_GEMCITABINE = list(
      description = "Gemcitabine (without another cytotoxic agent) regimen indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all regimen indicators 0 = reference regimen group)",
      notes = paste(
        "1 for gemcitabine alone or with atezolizumab, bevacizumab,",
        "necitumumab, pembrolizumab, pertuzumab and/or trastuzumab,",
        "ramucirumab or rituximab (REG 443, 202, 572, 10, 241, 308, 226, 81,",
        "326, 640, 168). Multiplies SLOPE by 0.416."
      ),
      source_name = "REG"
    ),
    TRT_PACLITAXEL = list(
      description = "Paclitaxel (without another cytotoxic agent) regimen indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all regimen indicators 0 = reference regimen group)",
      notes = paste(
        "1 for paclitaxel alone or with atezolizumab, bevacizumab, cetuximab,",
        "nivolumab, pembrolizumab, pertuzumab, ramucirumab, rituximab and/or",
        "trastuzumab (REG 315, 496, 89, 163, 565, 335, 134, 448, 513, 477,",
        "309, 478, 644, 540, 1, 16, 489, 663, 83, 178, 183, 662, 9, 58, 124,",
        "177, 11, 54). Multiplies SLOPE by 0.417; its MTT effect was fixed to 1."
      ),
      source_name = "REG"
    )
  )

  compartmentData <- list(
    depot_kpd = list(analyte = "chemotherapy", units = "normalized dose", specimen = "not applicable", verified = TRUE),
    depot = list(analyte = "filgrastim", units = "normalized dose", specimen = "administration site", verified = TRUE),
    central = list(analyte = "filgrastim", units = "normalized dose", specimen = "plasma", verified = TRUE),
    prol = list(analyte = "neutrophil precursors", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    transit1 = list(analyte = "neutrophil precursors", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    transit2 = list(analyte = "neutrophil precursors", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    transit3 = list(analyte = "neutrophil precursors", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    circ = list(analyte = "neutrophils", units = "10^9/L", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 5952L,
    n_studies = 1L,
    age_range = "19.6-89.4 years (training set); 19.6-89.4 years overall",
    age_median = "57.7 years (training set)",
    weight_range = "28-187 kg (training set)",
    weight_median = "71 kg (training set)",
    sex_female_pct = 75.8,
    race_ethnicity = c(White = 74.1, Asian = 9.2, Black = 8.6, Other = 4.0, Unknown = 3.9, AIAN_NHPI = 0.2),
    disease_state = paste(
      "Metastatic breast cancer (52%), metastatic thoracic cancer (NSCLC, SCLC",
      "or mesothelioma; 35%) or lymphoma (DLBCL or follicular; 14%) treated",
      "with chemotherapy at Memorial Sloan Kettering Cancer Center, 2016-2020"
    ),
    dose_range = paste(
      "Real-world chemotherapy regimens (first three cycles per patient),",
      "with filgrastim or pegfilgrastim in 2923 of 5952 training patients"
    ),
    regions = "United States (single centre)",
    notes = paste(
      "Retrospective EHR cohort of 9121 patients randomly split into",
      "training (5952), cross-validation (1339) and testing (1830) sets",
      "(Table 1); the PKPD model was estimated on the training set with",
      "21,831 ANC observations (median 3 per patient). Serum albumin median",
      "4.2 g/dL (range 1.7-5.3). Percentages above are training-set values",
      "from Table 1."
    )
  )

  ini({
    # Chemotherapy K-PD and G-CSF PK -- all fixed (no drug concentrations in the data).
    lkel <- fixed(log(5 / 50)); label("Chemotherapy K-PD elimination rate CLchemo/V (1/h)") # Table 2: CLchemo = 5 L/h fixed, V = 50 L fixed; Supplement 2 $DES DADT(1) = -CL*A(1)/V
    lvc <- fixed(log(50)); label("Volume V shared by the chemotherapy K-PD state and G-CSF central (L)") # Table 2: V = 50 L fixed; Supplement 2 $PK 'GCSF V fixed to same as chemo V'
    lka <- fixed(log(0.123)); label("G-CSF SC absorption rate, filgrastim (1/h)") # Table 2: kA = 0.123 h-1 fixed (Melhem 2018)
    e_peg_ka <- fixed(log(0.1528455)); label("Log ratio of pegfilgrastim to filgrastim absorption rate (unitless)") # Supplement 2 $PK: KSC = 0.123 * (0.1528455 ** PEG); Table 2 rounds to 0.153
    lcl <- fixed(log(13.3)); label("G-CSF clearance, filgrastim (L/h)") # Table 2: CLGCSF = 13.3 L/h fixed; Supplement 2 THETA(19) = 13.3 FIX
    e_peg_cl <- fixed(log(0.258)); label("Log ratio of pegfilgrastim to filgrastim clearance (unitless)") # Table 2: CLGCSF = theta * 0.258^PEG
    lfdepot <- fixed(log(1)); label("G-CSF SC bioavailability, filgrastim (fraction)") # Supplement 2 $PK: F7 = 0.646 ** PEG (filgrastim F = 1)
    e_peg_fdepot <- fixed(log(0.646)); label("Log bioavailability of pegfilgrastim relative to filgrastim (unitless)") # Table 2: F = theta^PEG = 0.646

    # Neutrophil PD
    lslope <- log(106.682); label("Linear chemotherapy kill coefficient SLOPE (L per normalized dose)") # Supplement 2 THETA(3) = 106.682; Table 2: SLOPE = 107 (RSE 2.8%)
    lmtt <- log(126.246); label("Mean transit time through the maturation chain MTT (h)") # Supplement 2 THETA(4) = 126.246; Table 2: MTT = 126 h (RSE 1.4%)
    lgamma <- fixed(log(0.16)); label("Feedback exponent gamma (unitless)") # Table 2: FB exponent = 0.16 fixed (Kloft 2006 mean); Supplement 2 THETA(5) = 0.16 FIX
    lstm1 <- log(0.960093); label("G-CSF stimulation of proliferation ST1 (per normalized-dose unit of G-CSF central amount)") # Supplement 2 THETA(20) = 0.960093; Table 2: ST1 = 0.960 (RSE 4.9%)
    lstm2 <- log(0.701233); label("G-CSF stimulation of maturation ST2 (per normalized-dose unit of G-CSF central amount)") # Supplement 2 THETA(21) = 0.701233; Table 2: ST2 = 0.701 (RSE 6.2%)

    # Continuous covariate effects (power models)
    e_age_slope <- 0.751727; label("Power exponent of AGE/57.2 on SLOPE (unitless)") # Supplement 2 THETA(27) = 0.751727; Table 2: Sage = 0.752 (RSE 11%)
    e_age_stm1 <- 0.160139; label("Power exponent of AGE/57.2 on ST1 (unitless)") # Supplement 2 THETA(28) = 0.160139; Table 2: Page = 0.160 (RSE 17%)
    e_alb_mtt <- 0.325695; label("Power exponent of ALB/4.1 g/dL on MTT (unitless)") # Supplement 2 THETA(29) = 0.325695; Table 2: Malb = 0.326 (RSE 6.0%)

    # Regimen effects on SLOPE (multiplicative theta^indicator; Figure 2b)
    e_trt_cyclophosphamide_doxorubicin_slope <- log(1.68205); label("Log SLOPE ratio, AC regimen (unitless)") # Supplement 2 THETA(6) = 1.68205 (AC)
    e_trt_carboplatin_etoposide_slope <- log(1.43615); label("Log SLOPE ratio, carboplatin + etoposide (unitless)") # Supplement 2 THETA(7) = 1.43615
    e_trt_carboplatin_paclitaxel_slope <- log(0.682858); label("Log SLOPE ratio, carboplatin + paclitaxel (unitless)") # Supplement 2 THETA(9) = 0.682858
    e_trt_carboplatin_pemetrexed_slope <- log(1.41292); label("Log SLOPE ratio, carboplatin + pemetrexed (unitless)") # Supplement 2 THETA(10) = 1.41292
    e_trt_cisplatin_pemetrexed_slope <- log(1.74323); label("Log SLOPE ratio, cisplatin + pemetrexed (unitless)") # Supplement 2 THETA(11) = 1.74323; Table 2 Sdrug upper bound 1.74
    e_trt_cyclophosphamide_fluorouracil_methotrexate_slope <- log(0.392131); label("Log SLOPE ratio, CMF regimen (unitless)") # Supplement 2 THETA(12) = 0.392131
    e_trt_docetaxel_slope <- log(0.26865); label("Log SLOPE ratio, docetaxel (unitless)") # Supplement 2 THETA(13) = 0.26865; Table 2 Sdrug lower bound 0.269
    e_trt_gemcitabine_slope <- log(0.416426); label("Log SLOPE ratio, gemcitabine (unitless)") # Supplement 2 THETA(14) = 0.416426
    e_trt_paclitaxel_slope <- log(0.416792); label("Log SLOPE ratio, paclitaxel (unitless)") # Supplement 2 THETA(16) = 0.416792

    # Regimen effects on MTT (multiplicative theta^indicator; Figure 2b)
    e_trt_cyclophosphamide_doxorubicin_mtt <- log(0.457641); label("Log MTT ratio, AC regimen (unitless)") # Supplement 2 THETA(22) = 0.457641; Table 2 Mdrug lower bound 0.458
    e_trt_carboplatin_pemetrexed_mtt <- log(1.71121); label("Log MTT ratio, carboplatin + pemetrexed (unitless)") # Supplement 2 THETA(23) = 1.71121
    e_trt_cisplatin_pemetrexed_mtt <- log(1.88982); label("Log MTT ratio, cisplatin + pemetrexed (unitless)") # Supplement 2 THETA(24) = 1.88982; Table 2 Mdrug upper bound 1.89
    e_trt_cyclophosphamide_fluorouracil_methotrexate_mtt <- log(1.5689); label("Log MTT ratio, CMF regimen (unitless)") # Supplement 2 THETA(25) = 1.5689

    # IIV (exponential). Supplement 2 $OMEGA gives variances; Table 2 reports their square roots as %.
    etalslope ~ 0.608128 # Supplement 2 $OMEGA 1 = 0.608128; Table 2 omega SLOPE 78.0%
    etalmtt ~ 0.108427 # Supplement 2 $OMEGA 2 = 0.108427; Table 2 omega MTT 32.9%
    etalstm1 ~ 0.00997778 # Supplement 2 $OMEGA 4 = 0.00997778; Table 2 omega ST1 9.99%
    etalcirc0 ~ fixed(0.162848) # Supplement 2: CIRC0 = CIRCBL*EXP(ETA(3)*SIGMA(1,1)**0.5) with unit $OMEGA 3 and $SIGMA 0.162848; Table 2 omega CIRC0 40.4%, equal to the residual variability

    # Residual error: additive on log(ANC)
    expSd <- sqrt(0.162848); label("Residual SD of log(ANC) (log scale)") # Supplement 2 $SIGMA 0.162848 with Y = LOG(IPRED + FLAG) + EPS(1); Table 2 sigma 40.4%
  })

  model({
    # Covariate transformations
    alb_gdl <- ALB / 10
    circbl <- NEUT
    if (circbl < 0.5) circbl <- 0.5

    # Regimen multipliers (all indicators 0 = reference regimen group)
    drg_slope <- exp(
      e_trt_cyclophosphamide_doxorubicin_slope * TRT_CYCLOPHOSPHAMIDE_DOXORUBICIN +
        e_trt_carboplatin_etoposide_slope * TRT_CARBOPLATIN_ETOPOSIDE +
        e_trt_carboplatin_paclitaxel_slope * TRT_CARBOPLATIN_PACLITAXEL +
        e_trt_carboplatin_pemetrexed_slope * TRT_CARBOPLATIN_PEMETREXED +
        e_trt_cisplatin_pemetrexed_slope * TRT_CISPLATIN_PEMETREXED +
        e_trt_cyclophosphamide_fluorouracil_methotrexate_slope * TRT_CYCLOPHOSPHAMIDE_FLUOROURACIL_METHOTREXATE +
        e_trt_docetaxel_slope * TRT_DOCETAXEL +
        e_trt_gemcitabine_slope * TRT_GEMCITABINE +
        e_trt_paclitaxel_slope * TRT_PACLITAXEL
    )
    drg_mtt <- exp(
      e_trt_cyclophosphamide_doxorubicin_mtt * TRT_CYCLOPHOSPHAMIDE_DOXORUBICIN +
        e_trt_carboplatin_pemetrexed_mtt * TRT_CARBOPLATIN_PEMETREXED +
        e_trt_cisplatin_pemetrexed_mtt * TRT_CISPLATIN_PEMETREXED +
        e_trt_cyclophosphamide_fluorouracil_methotrexate_mtt * TRT_CYCLOPHOSPHAMIDE_FLUOROURACIL_METHOTREXATE
    )

    # Individual parameters
    kel <- exp(lkel)
    vc <- exp(lvc)
    ka <- exp(lka + e_peg_ka * FORM_GCSF_PEG)
    cl <- exp(lcl + e_peg_cl * FORM_GCSF_PEG)
    fdepot <- exp(lfdepot + e_peg_fdepot * FORM_GCSF_PEG)
    slope <- exp(lslope + etalslope) * drg_slope * (AGE / 57.2)^e_age_slope
    mtt <- exp(lmtt + etalmtt) * drg_mtt * (alb_gdl / 4.1)^e_alb_mtt
    gamma <- exp(lgamma)
    stm1 <- exp(lstm1 + etalstm1) * (AGE / 57.2)^e_age_stm1
    stm2 <- exp(lstm2)
    circ0 <- circbl * exp(etalcirc0)
    ktr <- 4 / mtt

    # Drug effects. G-CSF acts on its central AMOUNT (Supplement 2 $DES ST1/ST2 use A(8)).
    edrug <- slope * depot_kpd / vc
    st1 <- stm1 * central
    st2 <- 1 + stm2 * central

    # Feedback acts only while circulating neutrophils are below baseline.
    feedb <- 1
    if (circ < circ0) feedb <- (circ0 / circ)^gamma

    d/dt(depot_kpd) <- -kel * depot_kpd
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central
    d/dt(prol) <- ktr * prol * (1 - edrug) * feedb + ktr * st1 * prol - ktr * st2 * prol
    d/dt(transit1) <- ktr * st2 * prol - ktr * st2 * transit1
    d/dt(transit2) <- ktr * st2 * transit1 - ktr * st2 * transit2
    d/dt(transit3) <- ktr * st2 * transit2 - ktr * st2 * transit3
    d/dt(circ) <- ktr * st2 * transit3 - ktr * circ

    f(depot) <- fdepot

    prol(0) <- circ0
    transit1(0) <- circ0
    transit2(0) <- circ0
    transit3(0) <- circ0
    circ(0) <- circ0

    ANC <- circ
    ANC ~ lnorm(expSd)
  })
}
