Thakre_2025_venetoclax <- function() {
  description <- paste(
    "Semi-mechanistic population PK/PD model of venetoclax plus azacitidine",
    "in treatment-naive adults with higher-risk myelodysplastic syndromes",
    "(MDS). PK: two-compartment venetoclax model with first-order absorption,",
    "an absorption lag and first-order elimination (the legacy venetoclax",
    "model re-estimated in MDS), with food and a power function of dose on",
    "relative bioavailability, CYP3A / OATP1B3 inhibitors and rituximab on",
    "CL/F, and sex and patient status on V2/F. PD: normal and MDS",
    "hematopoietic stem cells in the bone marrow differentiate into three",
    "parallel Friberg-type lineages (neutrophils, platelets, red blood cells;",
    "proliferating precursor, three transit compartments, circulating pool)",
    "with feedback from the circulating counts relative to healthy reference",
    "values. MDS stem cells proliferate faster and differentiate 100-fold",
    "less than normal stem cells; marrow crowding raises apoptosis of every",
    "marrow cell type. Venetoclax concentration linearly increases apoptosis",
    "of normal and MDS stem cells and of neutrophil and platelet precursors;",
    "azacitidine (an on/off flag, days 1-7 of each 28-day cycle) lowers",
    "proliferation of stem cells and of neutrophil and platelet precursors.",
    "Outputs: venetoclax Cc, percent bone-marrow blasts, and blood",
    "neutrophil, platelet and red blood cell counts."
  )
  reference <- paste(
    "Thakre N, Maier C, Zha J, Engelhardt B, Wolff JE, Mensing S.",
    "Semi-mechanistic population PK/PD model to aid clinical understanding",
    "of myelodysplastic syndromes following treatment with Venetoclax and",
    "Azacitidine. CPT Pharmacometrics Syst Pharmacol. 2025;14(3):448-458.",
    "doi:10.1002/psp4.13284.",
    "PD parameter values from Table 2; PK parameter values from Supplementary",
    "Table S1; model structure from the NONMEM control stream in the",
    "Supporting Information (Data S1). The legacy PK model is Jones AK et al.,",
    "AAPS J 2016;18(5):1192-1202, doi:10.1208/s12248-016-9927-9."
  )
  vignette <- "Thakre_2025_venetoclax"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "mg/L (= ug/mL)"
  )

  covariateData <- list(
    BMBLAST_PCT = list(
      description = "Baseline percentage of blasts in the bone marrow (aspirate, or biopsy when no aspirate was taken).",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column PERBLA. Sets the initial MDS stem-cell pool,",
        "hsc_mds(0) = BMBLAST_PCT / 100 x 10^12 marrow mononuclear cells,",
        "and the initial normal stem-cell pool as 0.1 x hsc_mds(0). Must be",
        "> 0. Study range 1-23% (median 12%, Table 1)."
      ),
      source_name = "PERBLA"
    ),
    NEUT = list(
      description = "Observed baseline (pre-treatment) absolute neutrophil count.",
      units = "10^9 cells/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column NEUT, used only as the per-subject initial condition of",
        "the circulating neutrophil pool (circ_anc(0) = NEUT x 5 L blood",
        "volume) and of the steady-state precursor pools. Time-fixed. Must be",
        "> 0. Study range 0.11-6.58 x 10^9/L (median 0.86, Table 1)."
      ),
      source_name = "NEUT"
    ),
    PLT_BASE = list(
      description = "Observed baseline (pre-treatment) platelet count.",
      units = "10^9 cells/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column PLAT. Initial condition of the circulating platelet",
        "pool (circ_plt(0) = PLT_BASE x 5 L) and of the steady-state",
        "precursor pools. Time-fixed. Study range 5-430 x 10^9/L (median 57,",
        "Table 1)."
      ),
      source_name = "PLAT"
    ),
    RBC_BASE = list(
      description = "Observed baseline (pre-treatment) red blood cell count.",
      units = "10^12 cells/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column RBC, which the source data set carries in 10^9 cells/L",
        "(Table 1: median 2975, range 1584.89-4300 x 10^9/L). The canonical",
        "unit is 10^12 cells/L, so divide the source value by 1000 (median",
        "2.975 x 10^12/L); model() multiplies by 1000 to return to the",
        "10^9-cell scale of the ODEs. Initial condition of the circulating",
        "red-cell pool and the steady-state precursor pools. Time-fixed."
      ),
      source_name = "RBC"
    ),
    CONMED_AZACITIDINE = list(
      description = "Azacitidine administration flag: 1 on the days azacitidine is given, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no azacitidine on that day)",
      notes = paste(
        "Source column REAZA. TIME-VARYING: the source has no azacitidine PK",
        "and implements its administration as a 0/1 flag (Methods). In the",
        "study and in the source simulations azacitidine 75 mg/m^2 was given",
        "once daily on days 1-7 of each 28-day cycle, so set",
        "CONMED_AZACITIDINE = 1 from day 0 to day 7 of every cycle and 0 for",
        "the remaining 21 days. Enters as faza^CONMED_AZACITIDINE on stem-cell",
        "and on neutrophil / platelet precursor proliferation."
      ),
      source_name = "REAZA"
    ),
    DOSE = list(
      description = "Administered venetoclax dose (mg) of the current dose record.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Use case (a) of the DOSE register entry: drives the dose",
        "nonlinearity on relative bioavailability, F1 x (DOSE / 400)^-0.183",
        "(Table S1 'Nonlinearity with Dose'; control stream",
        "'F1 = F1ETA*(DOSE/400)**(-0.183)'). Supply it on dose and",
        "observation records. Must be > 0; for a subject who receives no",
        "venetoclax (azacitidine monotherapy) give any positive value, e.g.",
        "400, since no dose is ever absorbed. Study doses 100, 200, 400 and",
        "800 mg once daily."
      ),
      source_name = "DOSE"
    ),
    FED = list(
      description = "Fed-vs-fasted indicator for the dose record: 1 = dosed with food, 0 = dosed fasting.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 with FED_LOWFAT = 1 (low-fat meal) is the bioavailability reference",
      notes = paste(
        "Legacy-model food levels (Table S1): fasting, low fat (reference),",
        "moderate fat, high fat, and 'any meal' (fed, meal type not",
        "recorded). FED = 0 selects the fasting factor 0.334. With FED = 1,",
        "exactly one of FED_LOWFAT / FED_MEDFAT / FED_HIGHFAT selects that",
        "meal; FED = 1 with all three 0 selects 'any meal' (factor 1.22),",
        "the condition the source simulations assumed."
      ),
      source_name = "FOOD"
    ),
    FED_LOWFAT = list(
      description = "Low-fat-meal indicator for the dose record: 1 = dosed after a low-fat meal.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (low-fat meal; the bioavailability reference, F1 factor 1)",
      notes = "Table S1 expresses every food effect relative to the low-fat meal.",
      source_name = "FOOD"
    ),
    FED_MEDFAT = list(
      description = "Moderate-fat-meal indicator for the dose record: 1 = dosed after a moderate-fat meal.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "F1 x 1.32 against the low-fat reference (Table S1 'Food effect on F1 (Moderate Fat Meal vs Low Fat Meal)').",
      source_name = "FOOD"
    ),
    FED_HIGHFAT = list(
      description = "High-fat-meal indicator for the dose record: 1 = dosed after a high-fat meal.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "F1 x 1.43 against the low-fat reference (Table S1 'Food effect on F1 (High Fat Meal vs Low Fat Meal)').",
      source_name = "FOOD"
    ),
    CONMED_CYP3A4_INH_STRONG = list(
      description = "Strong CYP3A inhibitor coadministered: 1 = yes.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no strong CYP3A inhibitor)",
      notes = paste(
        "Source column C3A4_HIB = 3. CL/F x 0.175 (Table S1). Mutually",
        "exclusive with CONMED_CYP3A4_INH_MOD."
      ),
      source_name = "C3A4_HIB"
    ),
    CONMED_CYP3A4_INH_MOD = list(
      description = "Moderate (but no strong) CYP3A inhibitor coadministered: 1 = yes.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no moderate or strong CYP3A inhibitor)",
      notes = "Source column C3A4_HIB = 2. CL/F x 0.819 (Table S1).",
      source_name = "C3A4_HIB"
    ),
    CONMED_OATP1B_INH = list(
      description = "OATP1B3 inhibitor coadministered: 1 = yes.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no OATP1B3 inhibitor)",
      notes = paste(
        "Source column OATP_HIB >= 1. CL/F x 0.853, a value printed only in",
        "the supplementary control stream ('IF(OATP_HIB.GE.1) CL = CL*0.853');",
        "Table S1 does not list the OATP effect, although the Methods state",
        "that the legacy covariate components were retained."
      ),
      source_name = "OATP_HIB"
    ),
    CONMED_RITUX = list(
      description = "Rituximab coadministered: 1 = yes.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no rituximab)",
      notes = paste(
        "Legacy-model covariate (CLL/NHL combination therapy); CL/F x 1.22",
        "(Table S1 'Rituximab Coadministration on CL/F'). No MDS patient",
        "received rituximab, so set 0 for this population."
      ),
      source_name = "RITUXIMAB"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator: 1 = healthy volunteer, 0 = patient.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (healthy volunteer is the V2/F reference of the legacy model)",
      notes = paste(
        "Legacy-model covariate: Table S1 'Healthy vs Diseased on V2/F' = 1.70,",
        "applied as 1.70^(1 - DIS_HEALTHY) so that patients carry a 1.70-fold",
        "higher V2/F, the orientation of the same legacy parameter in",
        "Samineni 2022 (1.71 for CLL/SLL/NHL patients vs healthy). Every MDS",
        "patient has DIS_HEALTHY = 0."
      ),
      source_name = "HV"
    ),
    SEXF = list(
      description = "Sex: 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Table S1 'Sex on V2/F' = 0.689; control stream",
        "'IF(SEX.EQ.1) V2 = V2*0.689'. The female orientation follows the",
        "legacy model (Samineni 2022 'V2/F,SEX=female' = 0.68). 14 of 57",
        "patients (25%) were female."
      ),
      source_name = "SEX"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "venetoclax", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "venetoclax", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "venetoclax", units = "mg", specimen = "plasma", verified = TRUE),
    hsc = list(
      analyte = "normal hematopoietic stem cells",
      units = "10^9 cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    hsc_mds = list(
      analyte = "MDS stem cells (bone-marrow blasts)",
      units = "10^9 cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    prol_anc = list(
      analyte = "granulocyte precursors",
      units = "10^9 cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit1_anc = list(
      analyte = "maturing granulocytes",
      units = "10^9 cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2_anc = list(
      analyte = "maturing granulocytes",
      units = "10^9 cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3_anc = list(
      analyte = "maturing granulocytes",
      units = "10^9 cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    circ_anc = list(analyte = "neutrophils", units = "10^9 cells", specimen = "whole blood", verified = TRUE),
    prol_plt = list(
      analyte = "thrombocyte precursors",
      units = "10^9 cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit1_plt = list(
      analyte = "maturing thrombocytes",
      units = "10^9 cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2_plt = list(
      analyte = "maturing thrombocytes",
      units = "10^9 cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3_plt = list(
      analyte = "maturing thrombocytes",
      units = "10^9 cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    circ_plt = list(analyte = "platelets", units = "10^9 cells", specimen = "whole blood", verified = TRUE),
    prol_rbc = list(
      analyte = "erythroid precursors",
      units = "10^9 cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit1_rbc = list(
      analyte = "maturing erythroid cells",
      units = "10^9 cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2_rbc = list(
      analyte = "maturing erythroid cells",
      units = "10^9 cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3_rbc = list(
      analyte = "maturing erythroid cells",
      units = "10^9 cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    circ_rbc = list(analyte = "red blood cells", units = "10^9 cells", specimen = "whole blood", verified = TRUE)
  )

  # Paper-specific PD states: the normal and MDS stem-cell pools at the head
  # of the three lineages, and the red-cell lineage of the Friberg chain (the
  # register carries cell-type suffixes for anc / plt but not for rbc).
  paper_specific_compartments <- c(
    "hsc",
    "hsc_mds",
    "prol_rbc",
    "transit1_rbc",
    "transit2_rbc",
    "transit3_rbc",
    "circ_rbc"
  )

  population <- list(
    species = "human",
    n_subjects = 57L,
    n_studies = 1L,
    n_observations = paste(
      "342 venetoclax plasma concentrations, 242 bone-marrow blast",
      "(biopsy / aspirate) samples, 1940 neutrophil, 1951 platelet and 1780",
      "red blood cell counts"
    ),
    age_range = "26-85 years (median 71, mean 68.7)",
    weight_range = "51.5-174 kg (median 82.5, mean 84.7)",
    sex_female_pct = 24.6,
    race_ethnicity = c(White = 96.5, Asian = 3.5),
    disease_state = paste(
      "Treatment-naive higher-risk myelodysplastic syndromes. Baseline",
      "bone-marrow blasts 1-23% (median 12%); neutrophils 0.11-6.58 x 10^9/L",
      "(median 0.86); platelets 5-430 x 10^9/L (median 57); red blood cells",
      "1.58-4.30 x 10^12/L (median 2.98)."
    ),
    dose_range = paste(
      "Venetoclax 400 or 800 mg once daily for 28 of 28 days (original",
      "protocol, 10 patients, with a cycle-1 ramp-up) or 100, 200 or 400 mg",
      "once daily for 14 of 28 days (revised protocol); 2 patients received",
      "azacitidine alone. Azacitidine 75 mg/m^2 once daily on days 1-7 of",
      "each 28-day cycle."
    ),
    renal_function = "normal 51%, mild impairment 37%, moderate impairment 12%",
    hepatic_function = "normal 89%, mild impairment 5.3%, moderate impairment 5.3%",
    notes = paste(
      "Phase 1b study NCT02942290 (M15-531), data cut-off 21 August 2019;",
      "demographics from Table 1. The PK layer was re-estimated in MDS with",
      "the legacy venetoclax model (505 subjects with CLL/SLL, NHL and",
      "healthy volunteers) as an informative prior."
    )
  )

  ini({
    # ---- Venetoclax PK (Table S1; time in day, dose in mg, conc in mg/L) ----
    lka <- log(3.71); label("First-order absorption rate constant ka (1/day)") # Table S1 'Ka (1/day)' = 3.71
    lcl <- log(451); label("Apparent clearance CL/F (L/day)") # Table S1 'CL/F (L/day)' = 451
    lvc <- log(117); label("Apparent central volume of distribution V2/F in healthy volunteers (L)") # Table S1 'V2/F (L)' = 117
    lq <- log(98.0); label("Apparent intercompartmental clearance Q/F (L/day)") # Table S1 'Q/F (L/day)' = 98.0
    lvp <- log(120); label("Apparent peripheral volume of distribution V3/F (L)") # Table S1 'V3/F (L)' = 120
    ltlag <- log(0.0533); label("Absorption lag time (day)") # Table S1 'Lag time in absorption' = 0.0533
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1 at the reference condition: low-fat meal, 400 mg (unitless)") # reference level; control stream 'F1 = 1'

    e_fasted_fdepot <- 0.334; label("Multiplicative factor on F1 when fasting vs a low-fat meal (unitless)") # Table S1 'Food effect on F1 (Fasting vs Low Fat Meal)' = 0.334
    e_medfat_fdepot <- 1.32; label("Multiplicative factor on F1 with a moderate-fat meal vs a low-fat meal (unitless)") # Table S1 'Food effect on F1 (Moderate Fat Meal vs Low Fat Meal)' = 1.32
    e_highfat_fdepot <- 1.43; label("Multiplicative factor on F1 with a high-fat meal vs a low-fat meal (unitless)") # Table S1 'Food effect on F1 (High Fat Meal vs Low Fat Meal)' = 1.43
    e_fed_fdepot <- 1.22; label("Multiplicative factor on F1 with any meal (type unrecorded) vs a low-fat meal (unitless)") # Table S1 'Food effect on F1 (Any Meal vs Low Fat Meal)' = 1.22
    e_dose_fdepot <- -0.183; label("Power exponent of (DOSE / 400 mg) on F1 (unitless)") # Table S1 'Nonlinearity with Dose' = -0.183

    e_conmed_cyp3a4_inh_mod_cl <- 0.819; label("Multiplicative factor on CL/F with a moderate CYP3A inhibitor (unitless)") # Table S1 'Moderate CYP3A Inhibitors on CL/F' = 0.819
    e_conmed_cyp3a4_inh_strong_cl <- 0.175; label("Multiplicative factor on CL/F with a strong CYP3A inhibitor (unitless)") # Table S1 'Strong CYP3A Inhibitors on CL/F' = 0.175
    e_conmed_oatp1b_inh_cl <- 0.853; label("Multiplicative factor on CL/F with an OATP1B3 inhibitor (unitless)") # supplementary control stream 'IF(OATP_HIB.GE.1) CL = CL*0.853'; not in Table S1
    e_conmed_ritux_cl <- 1.22; label("Multiplicative factor on CL/F with rituximab coadministration (unitless)") # Table S1 'Rituximab Coadministration on CL/F' = 1.22
    e_sexf_vc <- 0.689; label("Multiplicative factor on V2/F for women (unitless)") # Table S1 'Sex on V2/F' = 0.689
    e_dis_patient_vc <- 1.70; label("Multiplicative factor on V2/F for patients vs healthy volunteers (unitless)") # Table S1 'Healthy vs Diseased on V2/F' = 1.70

    # ---- Hematopoiesis PD (Table 2; final estimates) ----
    # The source parameterises every PD parameter as 10^(THETA + ETA). The
    # typical values below are the Table 2 natural-scale estimates on the
    # natural-log scale; the IIV variances are converted to the natural-log
    # scale further down.
    lmtt_plt <- log(57.7); label("Mean transit time of the platelet lineage MTTP (h)") # Table 2 'MTTP' = 57.7 h
    lkdel_plt <- log(0.0220); label("Elimination rate of circulating platelets kdel,P (1/day)") # Table 2 'kdel,P' = 0.0220 1/day
    lmtt_anc <- log(287); label("Mean transit time of the neutrophil lineage MTTN (h)") # Table 2 'MTTN' = 287 h
    lkdel_anc <- log(0.536); label("Elimination rate of circulating neutrophils kdel,N (1/day)") # Table 2 'kdel,N' = 0.536 1/day
    lmtt_rbc <- log(36.7); label("Mean transit time of the red blood cell lineage MTTR (h)") # Table 2 'MTTR' = 36.7 h
    lkdel_rbc <- log(0.00183); label("Elimination rate of circulating red blood cells kdel,R (1/day)") # Table 2 'kdel,R' = 0.00183 1/day
    lkprol_stem <- log(0.781); label("Proliferation rate of normal stem cells kprol,S (1/day)") # Table 2 'kprol,S' = 0.781 1/day
    lkmat_plt <- log(2.02e-5); label("Stem-cell maturation rate to the platelet precursor kmat,P (1/day)") # Table 2 'kmat,P' = 2.02 x 10^-5 1/day
    lkmat_anc <- log(0.00292); label("Stem-cell maturation rate to the neutrophil precursor kmat,N (1/day)") # Table 2 'kmat,N' = 0.00292 1/day
    lkmat_rbc <- log(0.00776); label("Stem-cell maturation rate to the red blood cell precursor kmat,R (1/day)") # Table 2 'kmat,R' = 0.00776 1/day
    lfp_mds <- log(0.000113); label("Fractional increase in proliferation of MDS stem cells fp (unitless)") # Table 2 'fp' = 0.000113
    lslope_stem <- log(0.169); label("Venetoclax slope on normal stem-cell apoptosis SLPS (L/mg)") # Table 2 'SLPS' = 0.169
    lslope_mds <- log(0.140); label("Venetoclax slope on MDS stem-cell apoptosis SLPMDS (L/mg)") # Table 2 'SLPMDS' = 0.140
    lslope_prol_anc <- log(0.241); label("Venetoclax slope on neutrophil precursor apoptosis SLPN1 (L/mg)") # Table 2 'SLPN1' = 0.241
    lslope_prol_plt <- log(0.543); label("Venetoclax slope on platelet precursor apoptosis SLPP1 (L/mg)") # Table 2 'SLPP1' = 0.543
    laza_stem <- log(1 / 0.760 - 1); label("Azacitidine effect on stem cells, log of (1/alphaSMDS - 1) (unitless)") # Table 2 'alphaSMDS' = 0.760 = 1/(1 + 10^THETA15); inverted to the log scale
    laza_prol <- log(1 / (1 - 0.788) - 2); label("Azacitidine effect on neutrophil and platelet precursors, log of (1/(1 - alphaN1P1) - 2) (unitless)") # Table 2 'alphaN1P1' = 0.788 = 1 - 1/(2 + 10^THETA19); inverted to the log scale
    lkdel_prol_plt <- log(32.7); label("Elimination rate of platelet precursors kdel,P1 (1/day)") # Table 2 'kdel,P1' = 32.7 1/day
    lkdel_prol_anc <- log(0.690); label("Elimination rate of neutrophil precursors kdel,N1 (1/day)") # Table 2 'kdel,N1' = 0.690 1/day
    lkdel_prol_rbc <- log(0.460); label("Elimination rate of red blood cell precursors kdel,R1 (1/day)") # Table 2 'kdel,R1' = 0.460 1/day

    # Fixed structural constants of the hematopoiesis model.
    gamma <- fixed(0.125); label("Feedback and crowding exponent for platelets and red blood cells (unitless)") # Methods 'gamma = 0.125'; control stream 'GAM = 0.125 ; fixed'
    gamma_anc <- fixed(0.148); label("Feedback exponent for neutrophils (unitless)") # Methods 'gammaN = 0.148'; control stream 'GAMN = 0.148 ;fixed'
    fdiff <- fixed(0.01); label("Impaired-differentiation factor of MDS stem cells fdiff (unitless)") # control stream 'IMPDIFF = 0.01 ; impaired differentiation of MDS stem cells (fixed)'
    fnorm <- fixed(0.1); label("Initial normal stem cells as a fraction of initial MDS stem cells (unitless)") # control stream 'FNORM = 0.1'; 'A04 = FNORM*A05'
    n_bm <- fixed(1000); label("Mononuclear cells in the bone marrow (10^9 cells)") # control stream 'A05 = (PERBLA/100)*1E12/1E9 ; ... assuming 10^12 mononuclear cells in BM'
    vblood <- fixed(5); label("Blood volume (L)") # control stream 'VPB = 5 ;blood volume'
    hb_plt <- fixed(300); label("Healthy reference platelet count HBP (10^9/L)") # control stream 'HBP = 300*VPB'
    hb_anc <- fixed(5.15); label("Healthy reference neutrophil count HBN (10^9/L)") # control stream 'HBN = 5.15*VPB'
    hb_rbc <- fixed(5400); label("Healthy reference red blood cell count HBR (10^9/L)") # control stream 'HBR = 5400*VPB'

    # ---- IIV ----
    # PK: Table S1 variances (log-normal).
    etalcl ~ 0.160 # Table S1 'IIV on CL/F' = 0.160 (41.7 %CV)
    etalvc ~ 0.203 # Table S1 'IIV on V2/F' = 0.203 (47.4 %CV)
    etalfdepot ~ 0.0935 # Table S1 'IIV on F1' = 0.0935 (31.3 %CV)
    # PD: Table 2 prints IIV as %CV = sqrt(exp(omega^2) - 1) x 100 of the
    # base-10 ETA, so omega^2(log10) = log(1 + CV^2), and the natural-log
    # variance used here is omega^2(log10) x log(10)^2.
    etalmtt_plt ~ 2.4786 # Table 2 MTTP IIV 77.2 %CV
    etalkdel_plt ~ 2.8830 # Table 2 kdel,P IIV 85.0 %CV
    etalmtt_anc ~ 0.4281 # Table 2 MTTN IIV 29.0 %CV
    etalkdel_anc ~ 1.8782 # Table 2 kdel,N IIV 65.2 %CV
    etalmtt_rbc ~ 0.9579 # Table 2 MTTR IIV 44.5 %CV
    etalkdel_rbc ~ 2.9459 # Table 2 kdel,R IIV 86.2 %CV
    etalkprol_stem ~ 1.0787 # Table 2 kprol,S IIV 47.5 %CV
    etalkmat_plt ~ 1.3082 # Table 2 kmat,P IIV 52.9 %CV
    etalkmat_anc ~ 7.1557 # Table 2 kmat,N IIV 169 %CV
    etalkmat_rbc ~ 4.2571 # Table 2 kmat,R IIV 111 %CV
    etalfp_mds ~ 1.2172 # Table 2 fp IIV 50.8 %CV
    etalslope_stem ~ 1.3878 # Table 2 SLPS IIV 54.7 %CV
    etalslope_mds ~ 0.7906 # Table 2 SLPMDS IIV 40.1 %CV
    etalslope_prol_anc ~ 4.5725 # Table 2 SLPN1 IIV 117 %CV
    etalaza_stem ~ 0.8164 # Table 2 alphaSMDS IIV 40.8 %CV
    etalkdel_prol_plt ~ 3.4153 # Table 2 kdel,P1 IIV 95.1 %CV
    etalkdel_prol_anc ~ 1.6961 # Table 2 kdel,N1 IIV 61.4 %CV
    etalkdel_prol_rbc ~ 1.9073 # Table 2 kdel,R1 IIV 65.8 %CV
    etalaza_prol ~ 2.1692 # Table 2 alphaN1P1 IIV 71.1 %CV
    etalslope_prol_plt ~ 8.4480 # Table 2 SLPP1 IIV 198 %CV

    # ---- Residual error ----
    # PK: Table S1 values are variances (they match the legacy-model sigma^2
    # in Samineni 2022), so SD = sqrt(variance).
    propSd <- sqrt(0.222); label("Proportional residual SD for venetoclax Cc (fraction)") # Table S1 'Proportional Error' = 0.222 (variance)
    addSd <- sqrt(3.07e-7); label("Additive residual SD for venetoclax Cc (mg/L)") # Table S1 'Additive Error' = 3.07 x 10^-7 (variance)
    # PD: Table 2 RUV values are variances (cf. control stream $SIGMA). Blasts
    # are additive on the logit scale; the three counts are additive on the
    # log10 scale, converted to the natural-log SD of lnorm() by x log(10).
    addSd_bmblast <- sqrt(0.859); label("Additive residual SD of bone-marrow blasts on the logit scale (unitless)") # Table 2 'RUV BM blasts' = 0.859 (variance)
    expSd_ANC <- sqrt(0.0914) * log(10); label("Log-scale residual SD of the neutrophil count (natural-log units)") # Table 2 'RUV neutrophils' = 0.0914 (variance, log10 scale)
    expSd_PLT <- sqrt(0.0424) * log(10); label("Log-scale residual SD of the platelet count (natural-log units)") # Table 2 'RUV platelets' = 0.0424 (variance, log10 scale)
    expSd_RBC <- sqrt(0.00202) * log(10); label("Log-scale residual SD of the red blood cell count (natural-log units)") # Table 2 'RUV RBCs' = 0.00202 (variance, log10 scale)
  })

  model({
    # ---- Venetoclax PK ----
    cl <- exp(lcl + etalcl) *
      e_conmed_cyp3a4_inh_strong_cl^CONMED_CYP3A4_INH_STRONG *
      e_conmed_cyp3a4_inh_mod_cl^CONMED_CYP3A4_INH_MOD *
      e_conmed_oatp1b_inh_cl^CONMED_OATP1B_INH *
      e_conmed_ritux_cl^CONMED_RITUX
    vc <- exp(lvc + etalvc) * e_dis_patient_vc^(1 - DIS_HEALTHY) * e_sexf_vc^SEXF
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka)
    tlag <- exp(ltlag)

    # Food factor on F1: fasting 0.334; low fat 1 (reference); moderate fat
    # 1.32; high fat 1.43; fed with no meal type recorded ('any meal') 1.22.
    f_meal <- FED_LOWFAT + FED_MEDFAT * e_medfat_fdepot + FED_HIGHFAT * e_highfat_fdepot +
      (1 - FED_LOWFAT - FED_MEDFAT - FED_HIGHFAT) * e_fed_fdepot
    f_food <- (1 - FED) * e_fasted_fdepot + FED * f_meal
    fdepot <- exp(lfdepot + etalfdepot) * f_food * (DOSE / 400)^e_dose_fdepot

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Hematopoiesis PD parameters ----
    mtt_plt <- exp(lmtt_plt + etalmtt_plt)
    kdel_plt <- exp(lkdel_plt + etalkdel_plt)
    kdel_prol_plt <- exp(lkdel_prol_plt + etalkdel_prol_plt)
    mtt_anc <- exp(lmtt_anc + etalmtt_anc)
    kdel_anc <- exp(lkdel_anc + etalkdel_anc)
    kdel_prol_anc <- exp(lkdel_prol_anc + etalkdel_prol_anc)
    mtt_rbc <- exp(lmtt_rbc + etalmtt_rbc)
    kdel_rbc <- exp(lkdel_rbc + etalkdel_rbc)
    kdel_prol_rbc <- exp(lkdel_prol_rbc + etalkdel_prol_rbc)
    kprol_stem <- exp(lkprol_stem + etalkprol_stem)
    kmat_plt <- exp(lkmat_plt + etalkmat_plt)
    kmat_anc <- exp(lkmat_anc + etalkmat_anc)
    kmat_rbc <- exp(lkmat_rbc + etalkmat_rbc)
    fp_mds <- exp(lfp_mds + etalfp_mds)
    slope_stem <- exp(lslope_stem + etalslope_stem)
    slope_mds <- exp(lslope_mds + etalslope_mds)
    slope_prol_anc <- exp(lslope_prol_anc + etalslope_prol_anc)
    slope_prol_plt <- exp(lslope_prol_plt + etalslope_prol_plt)
    # Azacitidine multipliers on proliferation: faza_stem in (0, 1) for stem
    # cells ('DEAZA = 1/(1+10**(MU_15+ETA(15)))') and faza_prol in (0.5, 1)
    # for neutrophil and platelet precursors ('DEAZAB = 1-1/(2+10**(...))').
    aza_stem <- exp(laza_stem + etalaza_stem)
    aza_prol <- exp(laza_prol + etalaza_prol)
    faza_stem <- 1 / (1 + aza_stem)
    faza_prol <- 1 - 1 / (2 + aza_prol)

    # Transit rate over the four maturation stages (proliferating precursor
    # plus three transit compartments); MTT in hours, rates in 1/day.
    ktr_plt <- 4 / (mtt_plt / 24 + 1e-6)
    ktr_anc <- 4 / (mtt_anc / 24 + 1e-6)
    ktr_rbc <- 4 / (mtt_rbc / 24 + 1e-6)

    # ---- Pre-treatment (diseased) steady state, 10^9 cells ----
    a0_hsc_mds <- BMBLAST_PCT / 100 * n_bm
    a0_hsc <- fnorm * a0_hsc_mds
    a0_circ_anc <- NEUT * vblood
    a0_circ_plt <- PLT_BASE * vblood
    a0_circ_rbc <- RBC_BASE * 1000 * vblood
    fb0_anc <- (hb_anc / NEUT)^gamma_anc
    fb0_plt <- (hb_plt / PLT_BASE)^gamma
    fb0_rbc <- (hb_rbc / (RBC_BASE * 1000))^gamma
    # Stem-cell apoptosis rates that hold the normal stem cells at steady
    # state at baseline; the MDS rate omits the (1 + fp) proliferation gain,
    # so the MDS clone expands without treatment.
    kdel_stem <- kprol_stem - (kmat_plt * fb0_plt + kmat_anc * fb0_anc + kmat_rbc * fb0_rbc)
    kdel_mds <- kprol_stem - fdiff * (kmat_plt * fb0_plt + kmat_anc * fb0_anc + kmat_rbc * fb0_rbc)
    # Precursor proliferation rates that hold each lineage at its observed
    # baseline count (control stream KPROLP / KPROLN / KPROLR).
    kprol_plt <- (a0_circ_plt * kdel_prol_plt * kdel_plt * mtt_plt +
      (-96 * a0_hsc_mds * fdiff - 96 * a0_hsc) * fb0_plt * kmat_plt + 96 * a0_circ_plt * kdel_plt) /
      (a0_circ_plt * kdel_plt * fb0_plt * mtt_plt)
    kprol_anc <- (a0_circ_anc * kdel_prol_anc * kdel_anc * mtt_anc +
      (-96 * a0_hsc_mds * fdiff - 96 * a0_hsc) * fb0_anc * kmat_anc + 96 * a0_circ_anc * kdel_anc) /
      (a0_circ_anc * kdel_anc * fb0_anc * mtt_anc)
    kprol_rbc <- (a0_circ_rbc * kdel_prol_rbc * kdel_rbc * mtt_rbc +
      (-96 * a0_hsc_mds * fdiff - 96 * a0_hsc) * fb0_rbc * kmat_rbc + 96 * a0_circ_rbc * kdel_rbc) /
      (a0_circ_rbc * kdel_rbc * fb0_rbc * mtt_rbc)
    a0_prol_anc <- ((a0_hsc_mds * fdiff + a0_hsc) * fb0_anc * kmat_anc) / (ktr_anc - fb0_anc * kprol_anc + kdel_prol_anc)
    a0_prol_plt <- ((a0_hsc_mds * fdiff + a0_hsc) * fb0_plt * kmat_plt) / (ktr_plt - fb0_plt * kprol_plt + kdel_prol_plt)
    a0_prol_rbc <- ((a0_hsc_mds * fdiff + a0_hsc) * fb0_rbc * kmat_rbc) / (ktr_rbc - fb0_rbc * kprol_rbc + kdel_prol_rbc)

    hsc(0) <- a0_hsc
    hsc_mds(0) <- a0_hsc_mds
    prol_anc(0) <- a0_prol_anc
    transit1_anc(0) <- a0_prol_anc
    transit2_anc(0) <- a0_prol_anc
    transit3_anc(0) <- a0_prol_anc
    circ_anc(0) <- a0_circ_anc
    prol_plt(0) <- a0_prol_plt
    transit1_plt(0) <- a0_prol_plt
    transit2_plt(0) <- a0_prol_plt
    transit3_plt(0) <- a0_prol_plt
    circ_plt(0) <- a0_circ_plt
    prol_rbc(0) <- a0_prol_rbc
    transit1_rbc(0) <- a0_prol_rbc
    transit2_rbc(0) <- a0_prol_rbc
    transit3_rbc(0) <- a0_prol_rbc
    circ_rbc(0) <- a0_circ_rbc

    # ---- ODEs ----
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag
    f(depot) <- fdepot

    Cc <- central / vc
    # Venetoclax concentration driving the PD (control stream
    # 'IF(CP.LE.1E-6) CP=0').
    cp <- Cc * (Cc > 1e-6)

    # Marrow crowding: total marrow cells relative to baseline (= 1 at t = 0).
    crowd <- ((hsc + hsc_mds + prol_anc + prol_plt + prol_rbc) /
      (a0_hsc + a0_hsc_mds + a0_prol_anc + a0_prol_plt + a0_prol_rbc))^gamma
    # Feedback of low circulating counts relative to the healthy reference.
    fb_anc <- (hb_anc * vblood / circ_anc + 1e-6)^gamma_anc
    fb_plt <- (hb_plt * vblood / circ_plt + 1e-6)^gamma
    fb_rbc <- (hb_rbc * vblood / circ_rbc + 1e-6)^gamma
    kmat_tot <- kmat_plt * fb_plt + kmat_anc * fb_anc + kmat_rbc * fb_rbc

    d/dt(hsc) <- faza_stem^CONMED_AZACITIDINE * kprol_stem * hsc - kmat_tot * hsc -
      kdel_stem * (1 + slope_stem * cp) * hsc * crowd
    d/dt(hsc_mds) <- faza_stem^CONMED_AZACITIDINE * (1 + fp_mds) * kprol_stem * hsc_mds -
      fdiff * kmat_tot * hsc_mds - kdel_mds * (1 + slope_mds * cp) * hsc_mds * crowd

    d/dt(prol_anc) <- faza_prol^CONMED_AZACITIDINE * fb_anc * kprol_anc * prol_anc - ktr_anc * prol_anc +
      kmat_anc * fb_anc * hsc + fdiff * kmat_anc * fb_anc * hsc_mds -
      kdel_prol_anc * prol_anc * (1 + slope_prol_anc * cp) * crowd
    d/dt(transit1_anc) <- ktr_anc * prol_anc - ktr_anc * transit1_anc
    d/dt(transit2_anc) <- ktr_anc * transit1_anc - ktr_anc * transit2_anc
    d/dt(transit3_anc) <- ktr_anc * transit2_anc - ktr_anc * transit3_anc
    d/dt(circ_anc) <- ktr_anc * transit3_anc - kdel_anc * circ_anc

    d/dt(prol_plt) <- faza_prol^CONMED_AZACITIDINE * fb_plt * kprol_plt * prol_plt - ktr_plt * prol_plt +
      kmat_plt * fb_plt * hsc + fdiff * kmat_plt * fb_plt * hsc_mds -
      kdel_prol_plt * (1 + slope_prol_plt * cp) * prol_plt * crowd
    d/dt(transit1_plt) <- ktr_plt * prol_plt - ktr_plt * transit1_plt
    d/dt(transit2_plt) <- ktr_plt * transit1_plt - ktr_plt * transit2_plt
    d/dt(transit3_plt) <- ktr_plt * transit2_plt - ktr_plt * transit3_plt
    d/dt(circ_plt) <- ktr_plt * transit3_plt - kdel_plt * circ_plt

    d/dt(prol_rbc) <- fb_rbc * kprol_rbc * prol_rbc - ktr_rbc * prol_rbc +
      kmat_rbc * fb_rbc * hsc + fdiff * kmat_rbc * fb_rbc * hsc_mds -
      kdel_prol_rbc * prol_rbc * crowd
    d/dt(transit1_rbc) <- ktr_rbc * prol_rbc - ktr_rbc * transit1_rbc
    d/dt(transit2_rbc) <- ktr_rbc * transit1_rbc - ktr_rbc * transit2_rbc
    d/dt(transit3_rbc) <- ktr_rbc * transit2_rbc - ktr_rbc * transit3_rbc
    d/dt(circ_rbc) <- ktr_rbc * transit3_rbc - kdel_rbc * circ_rbc

    # ---- Observations ----
    bmblast <- 100 * hsc_mds / n_bm # percent blasts in the bone marrow
    ANC <- circ_anc / vblood # 10^9 cells/L
    PLT <- circ_plt / vblood # 10^9 cells/L
    RBC <- circ_rbc / 1000 / vblood # 10^12 cells/L

    Cc ~ add(addSd) + prop(propSd)
    # Source: LOGIT = LOG((PBLA+0.1)/(100.1-PBLA)), additive on that scale,
    # which is logitNorm() on the interval (-0.1, 100.1).
    bmblast ~ logitNorm(addSd_bmblast, -0.1, 100.1)
    ANC ~ lnorm(expSd_ANC)
    PLT ~ lnorm(expSd_PLT)
    RBC ~ lnorm(expSd_RBC)
  })
}
