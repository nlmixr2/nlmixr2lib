Gallo_2021_remdesivir_pbpk <- function() {
  description <- paste(
    "PBPK (hybrid whole-body, Magnolia/ACSL). Remdesivir (RDV) and its",
    "alanine metabolite GS-704277, nucleoside monophosphate, nucleoside",
    "GS-441524 and active nucleoside triphosphate GS-443902 in healthy adults",
    "after intravenous infusion (Gallo 2021). Venous plasma RDV, GS-704277 and",
    "GS-441524 are a fitted forcing function: a linear two-compartment RDV",
    "model with first-order conversion RDV -> GS-704277 -> GS-441524 and",
    "first-order elimination of each metabolite. Venous plasma feeds a lung",
    "extracellular space that drains to arterial plasma, which perfuses eleven",
    "further tissues (adipose, bone, brain, gut, heart, kidney, liver, muscle,",
    "rest of body, skin, spleen; gut and spleen drain through the liver).",
    "Every tissue has an extracellular and an intracellular space; only RDV",
    "and GS-441524 cross the cell membrane (unbound-concentration flux with an",
    "intracellular:plasma partition coefficient), and each intracellular space",
    "carries the metabolic scheme RDV -> GS-704277 -> monophosphate <->",
    "GS-441524, monophosphate -> GS-443902 -> elimination. A separate",
    "peripheral blood mononuclear cell (PBMC) module, driven by venous plasma",
    "and written in first-order rate constants, was calibrated to reported",
    "PBMC GS-443902 data; the tissue clearances are those rate constants",
    "scaled by each tissue's intracellular volume. Tissue outflows leave the",
    "system (venous plasma is prescribed, not a mass balance). Species masses",
    "are carried without molecular-weight conversion, as in the authors' code;",
    "the micromolar outputs use the code's conversion factors. The 20 percent",
    "CV Monte-Carlo variability the paper applies is encoded as fixed",
    "log-normal etas; no residual error was reported."
  )
  reference <- paste(
    "Gallo JM. Hybrid physiologically-based pharmacokinetic model for",
    "remdesivir: Application to SARS-CoV-2.",
    "Clin Transl Sci. 2021;14:1082-1091. doi:10.1111/cts.12975.",
    "Parameter values are Supplementary Tables S1-S3 (CTS-14-1082-s001.pdf);",
    "the ODE system is transcribed from the author's deposited Magnolia code",
    "ffRDV_v50.csl and ffRDV_v50_multipleDose.csl at",
    "https://github.com/jmgallo/PBPK-Model-for-Remdesivir (commit ed83bc3,",
    "2020-10-26), which the supplement cites.",
    sep = " "
  )
  vignette <- "Gallo_2021_remdesivir_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  # Species masses are carried without molecular-weight conversion (as in the
  # deposited code), so 'mg' of a metabolite is mg of the species formed
  # 1:1 by mass from its precursor.
  compartmentData <- list(
    central = list(analyte = "remdesivir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "remdesivir", units = "mg", specimen = "plasma", verified = TRUE),
    central_gs704277 = list(analyte = "GS-704277", units = "mg/L", specimen = "plasma", verified = TRUE),
    central_gs441524 = list(analyte = "GS-441524", units = "mg/L", specimen = "plasma", verified = TRUE),
    arterial = list(analyte = "remdesivir", units = "mg", specimen = "plasma", verified = TRUE),
    arterial_gs704277 = list(analyte = "GS-704277", units = "mg", specimen = "plasma", verified = TRUE),
    arterial_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "plasma", verified = TRUE),
    pbmc = list(analyte = "remdesivir", units = "mg/L", specimen = "blood cell", verified = TRUE),
    pbmc_gs704277 = list(analyte = "GS-704277", units = "mg/L", specimen = "blood cell", verified = TRUE),
    pbmc_gs441524mp = list(
      analyte = "GS-441524 monophosphate",
      units = "mg/L",
      specimen = "blood cell",
      verified = TRUE
    ),
    pbmc_gs441524 = list(analyte = "GS-441524", units = "mg/L", specimen = "blood cell", verified = TRUE),
    pbmc_gs443902 = list(analyte = "GS-443902", units = "mg/L", specimen = "blood cell", verified = TRUE),
    is_adipose = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_adipose = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_adipose_gs704277 = list(analyte = "GS-704277", units = "mg", specimen = "tissue", verified = TRUE),
    int_adipose_gs441524mp = list(
      analyte = "GS-441524 monophosphate",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    is_adipose_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_adipose_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_adipose_gs443902 = list(analyte = "GS-443902", units = "mg", specimen = "tissue", verified = TRUE),
    is_bone = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_bone = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_bone_gs704277 = list(analyte = "GS-704277", units = "mg", specimen = "tissue", verified = TRUE),
    int_bone_gs441524mp = list(analyte = "GS-441524 monophosphate", units = "mg", specimen = "tissue", verified = TRUE),
    is_bone_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_bone_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_bone_gs443902 = list(analyte = "GS-443902", units = "mg", specimen = "tissue", verified = TRUE),
    is_brain = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_brain = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_brain_gs704277 = list(analyte = "GS-704277", units = "mg", specimen = "tissue", verified = TRUE),
    int_brain_gs441524mp = list(
      analyte = "GS-441524 monophosphate",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    is_brain_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_brain_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_brain_gs443902 = list(analyte = "GS-443902", units = "mg", specimen = "tissue", verified = TRUE),
    is_gut = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_gut = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_gut_gs704277 = list(analyte = "GS-704277", units = "mg", specimen = "tissue", verified = TRUE),
    int_gut_gs441524mp = list(analyte = "GS-441524 monophosphate", units = "mg", specimen = "tissue", verified = TRUE),
    is_gut_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_gut_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_gut_gs443902 = list(analyte = "GS-443902", units = "mg", specimen = "tissue", verified = TRUE),
    is_heart = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_heart = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_heart_gs704277 = list(analyte = "GS-704277", units = "mg", specimen = "tissue", verified = TRUE),
    int_heart_gs441524mp = list(
      analyte = "GS-441524 monophosphate",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    is_heart_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_heart_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_heart_gs443902 = list(analyte = "GS-443902", units = "mg", specimen = "tissue", verified = TRUE),
    is_kidney = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_kidney = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_kidney_gs704277 = list(analyte = "GS-704277", units = "mg", specimen = "tissue", verified = TRUE),
    int_kidney_gs441524mp = list(
      analyte = "GS-441524 monophosphate",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    is_kidney_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_kidney_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_kidney_gs443902 = list(analyte = "GS-443902", units = "mg", specimen = "tissue", verified = TRUE),
    is_liver = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_liver = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_liver_gs704277 = list(analyte = "GS-704277", units = "mg", specimen = "tissue", verified = TRUE),
    int_liver_gs441524mp = list(
      analyte = "GS-441524 monophosphate",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    is_liver_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_liver_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_liver_gs443902 = list(analyte = "GS-443902", units = "mg", specimen = "tissue", verified = TRUE),
    is_lung = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_lung = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_lung_gs704277 = list(analyte = "GS-704277", units = "mg", specimen = "tissue", verified = TRUE),
    int_lung_gs441524mp = list(analyte = "GS-441524 monophosphate", units = "mg", specimen = "tissue", verified = TRUE),
    is_lung_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_lung_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_lung_gs443902 = list(analyte = "GS-443902", units = "mg", specimen = "tissue", verified = TRUE),
    is_lung_gs704277 = list(analyte = "GS-704277", units = "mg", specimen = "tissue", verified = TRUE),
    is_muscle = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_muscle = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_muscle_gs704277 = list(analyte = "GS-704277", units = "mg", specimen = "tissue", verified = TRUE),
    int_muscle_gs441524mp = list(
      analyte = "GS-441524 monophosphate",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    is_muscle_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_muscle_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_muscle_gs443902 = list(analyte = "GS-443902", units = "mg", specimen = "tissue", verified = TRUE),
    is_remainder = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_remainder = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_remainder_gs704277 = list(analyte = "GS-704277", units = "mg", specimen = "tissue", verified = TRUE),
    int_remainder_gs441524mp = list(
      analyte = "GS-441524 monophosphate",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    is_remainder_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_remainder_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_remainder_gs443902 = list(analyte = "GS-443902", units = "mg", specimen = "tissue", verified = TRUE),
    is_skin = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_skin = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_skin_gs704277 = list(analyte = "GS-704277", units = "mg", specimen = "tissue", verified = TRUE),
    int_skin_gs441524mp = list(analyte = "GS-441524 monophosphate", units = "mg", specimen = "tissue", verified = TRUE),
    is_skin_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_skin_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_skin_gs443902 = list(analyte = "GS-443902", units = "mg", specimen = "tissue", verified = TRUE),
    is_spleen = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_spleen = list(analyte = "remdesivir", units = "mg", specimen = "tissue", verified = TRUE),
    int_spleen_gs704277 = list(analyte = "GS-704277", units = "mg", specimen = "tissue", verified = TRUE),
    int_spleen_gs441524mp = list(
      analyte = "GS-441524 monophosphate",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    is_spleen_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_spleen_gs441524 = list(analyte = "GS-441524", units = "mg", specimen = "tissue", verified = TRUE),
    int_spleen_gs443902 = list(analyte = "GS-443902", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "healthy adults (phase 1)",
    dose_range = paste(
      "single IV infusions of 3-225 mg over 0.5 or 2 h and a 150 mg daily",
      "1-h infusion multiple-dose cohort; simulated clinical regimen 200 mg",
      "loading dose then 100 mg daily for 4 days (1- or 2-h infusions)"
    ),
    weight_range = "physiology fixed to one 80 kg White adult male (PK-Sim, NHANES 1997)",
    notes = paste(
      "The model was built from digitized mean plasma RDV, GS-704277 and",
      "GS-441524 concentration-time data and reported Cmax/AUC for the",
      "healthy-volunteer phase 1 cohorts of Humeniuk et al. 2020 (single",
      "doses 3-225 mg as 0.5- or 2-h infusions, plus a 150 mg multiple-dose",
      "cohort), and on reported PBMC GS-443902 Cmax, C24 and AUC from those",
      "cohorts and from the EU compassionate-use assessment. Individual data",
      "were not used; subject counts and demographics are not reported by",
      "Gallo 2021. Organ volumes and plasma flows describe a single typical",
      "80 kg White male (PK-Sim, NHANES 1997) with hematocrit 0.45."
    )
  )

  ini({
    # =================================================================
    # FORCING FUNCTION -- Table S1. Fitted by maximum likelihood with an
    # additive error model to the digitized mean plasma data (Methods,
    # 'Basic characteristics of the hybrid PBPK model'); the paper
    # reports no standard errors and no residual-error magnitudes.
    # Venous plasma RDV is a linear two-compartment model; GS-704277 and
    # GS-441524 are first-order conversion / elimination chains written
    # directly in concentration (deposited code: dApv = frdva*RDVpv -
    # kela*Apv; dNpv = fan*Apv - keln*Npv).
    # =================================================================
    lcl <- log(44.6)
    label("Total systemic clearance of remdesivir (L/h)") # Table S1 CL = 44.6
    lvc <- log(4.6)
    label("Central (venous plasma) volume of remdesivir (L)") # Table S1 Vc = 4.6
    lk12 <- log(2.1)
    label("Central-to-peripheral transfer rate constant of remdesivir (1/h)") # Table S1 k12 = 2.1
    lk21 <- log(1.86)
    label("Peripheral-to-central transfer rate constant of remdesivir (1/h)") # Table S1 k21 = 1.86
    lk_gs704277_form <- log(0.17)
    label("Conversion rate constant, plasma remdesivir to GS-704277 (1/h)") # Table S1 krdva = 0.17
    lkel_gs704277 <- log(1.1)
    label("Elimination rate constant of plasma GS-704277 (1/h)") # Table S1 kela = 1.1
    lk_gs441524_form <- log(0.32)
    label("Conversion rate constant, plasma GS-704277 to GS-441524 (1/h)") # Table S1 kan = 0.32
    lkel_gs441524 <- log(0.04)
    label("Elimination rate constant of plasma GS-441524 (1/h)") # Table S1 keln = 0.04 (row printed 'Elimination rate constant for A'; symbol keln)

    # Plasma unbound fractions (Table S2 footnote '^'; Methods 'The unbound
    # fraction of RDV in human plasma has been reported as 0.12 and N within
    # the range of 1'). The paper assumes the extracellular unbound fraction
    # equals the plasma unbound fraction.
    fu_p <- fixed(0.12)
    label("Unbound fraction of remdesivir in plasma and extracellular fluid (unitless)") # Table S2 footnote: fu RDV = 0.12
    fu_p_gs441524 <- fixed(1)
    label("Unbound fraction of GS-441524 in plasma and extracellular fluid (unitless)") # Table S2 footnote: fu N = 1

    # =================================================================
    # PBMC MODULE -- Table S2 'PBMC' row (partition coefficients and the
    # '*' first-order transport rate constants) and the Table S3 footnote
    # (first-order metabolic / elimination rate constants). The PBMC
    # states are concentrations; the rate constants were 'arrived at
    # iteratively' against reported PBMC GS-443902 Cmax, C24 and AUC
    # (Results), so they are calibrated, not fixed from the literature.
    # =================================================================
    lkin_pbmc <- log(9)
    label("Plasma-PBMC transport rate constant of remdesivir (1/h)") # Table S2 PBMC Transport RDV = 9.0 (footnote '*': 1/h)
    lkin_pbmc_gs441524 <- log(1)
    label("Plasma-PBMC transport rate constant of GS-441524 (1/h)") # Table S2 PBMC Transport N = 1.0 (footnote '*': 1/h)
    lkp_pbmc <- fixed(log(1))
    label("PBMC:plasma partition coefficient of remdesivir (unitless)") # Table S2 PBMC Partition Coefficient RDV = 1.0
    lkp_pbmc_gs441524 <- fixed(log(1))
    label("PBMC:plasma partition coefficient of GS-441524 (unitless)") # Table S2 PBMC Partition Coefficient N = 1.0
    lk_rdva_pbmc <- log(10)
    label("PBMC metabolic rate constant, remdesivir to GS-704277 (1/h)") # Table S3 footnote mKrdva = 10.
    lk_amp_pbmc <- log(1)
    label("PBMC metabolic rate constant, GS-704277 to nucleoside monophosphate (1/h)") # Table S3 footnote mKamp = 1.
    lk_mpn_pbmc <- log(2)
    label("PBMC metabolic rate constant, nucleoside monophosphate to GS-441524 (1/h)") # Table S3 footnote mKmpn = 2.
    lk_nmp_pbmc <- log(0.5)
    label("PBMC metabolic rate constant, GS-441524 to nucleoside monophosphate (1/h)") # Table S3 footnote mKnmp = 0.5
    lk_mptn_pbmc <- log(10)
    label("PBMC metabolic rate constant, nucleoside monophosphate to GS-443902 (1/h)") # Table S3 footnote mKmptn = 10.
    lk_tn_pbmc <- log(0.03)
    label("PBMC elimination rate constant of GS-443902 (1/h)") # Table S3 footnote eKtn = 0.03

    # =================================================================
    # TISSUES -- Table S2 (partition coefficients, transport clearances)
    # and Table S3 (intracellular metabolic and elimination clearances).
    # Each clearance is the corresponding PBMC rate constant times the
    # tissue intracellular volume (Table S3 footnote), rounded as printed;
    # the printed values are used. Partition coefficients are PK-Sim
    # in silico predictions (Methods 'Parameters for the PBPK model').
    # =================================================================

    # ---- Adipose (Table S2 and Table S3 row 'Adipose') ----
    lkp_adipose <- fixed(log(12.2))
    label("Adipose intracellular:plasma partition coefficient of remdesivir (unitless)") # Table S2 Partition Coefficient RDV = 12.2
    lkp_adipose_gs441524 <- fixed(log(0.17))
    label("Adipose intracellular:plasma partition coefficient of GS-441524 (unitless)") # Table S2 Partition Coefficient N = 0.17
    lclin_adipose <- fixed(log(156.4))
    label("Adipose extracellular-intracellular transport clearance of remdesivir (L/h)") # Table S2 Transport Clearance RDV = 156.4
    lclin_adipose_gs441524 <- fixed(log(17.38))
    label("Adipose extracellular-intracellular transport clearance of GS-441524 (L/h)") # Table S2 Transport Clearance N = 17.38
    lcl_rdva_adipose <- fixed(log(173.8))
    label("Adipose metabolic clearance, remdesivir to GS-704277 (L/h)") # Table S3 mCLrdva = 173.8
    lcl_amp_adipose <- fixed(log(17.38))
    label("Adipose metabolic clearance, GS-704277 to nucleoside monophosphate (L/h)") # Table S3 mCLamp = 17.38
    lcl_mpn_adipose <- fixed(log(34.76))
    label("Adipose metabolic clearance, nucleoside monophosphate to GS-441524 (L/h)") # Table S3 mCLmpn = 34.76
    lcl_nmp_adipose <- fixed(log(8.68))
    label("Adipose metabolic clearance, GS-441524 to nucleoside monophosphate (L/h)") # Table S3 mCLnmp = 8.68
    lcl_mptn_adipose <- fixed(log(173.8))
    label("Adipose metabolic clearance, nucleoside monophosphate to GS-443902 (L/h)") # Table S3 mCLmptn = 173.8
    lcl_tn_adipose <- fixed(log(0.521))
    label("Adipose elimination clearance of GS-443902 (L/h)") # Table S3 eCLtn = 0.521

    # ---- Bone (Table S2 and Table S3 row 'Bone') ----
    lkp_bone <- fixed(log(4.23))
    label("Bone intracellular:plasma partition coefficient of remdesivir (unitless)") # Table S2 Partition Coefficient RDV = 4.23
    lkp_bone_gs441524 <- fixed(log(0.51))
    label("Bone intracellular:plasma partition coefficient of GS-441524 (unitless)") # Table S2 Partition Coefficient N = 0.51
    lclin_bone <- fixed(log(93.1))
    label("Bone extracellular-intracellular transport clearance of remdesivir (L/h)") # Table S2 Transport Clearance RDV = 93.1
    lclin_bone_gs441524 <- fixed(log(10.3))
    label("Bone extracellular-intracellular transport clearance of GS-441524 (L/h)") # Table S2 Transport Clearance N = 10.3
    lcl_rdva_bone <- fixed(log(100.3))
    label("Bone metabolic clearance, remdesivir to GS-704277 (L/h)") # Table S3 mCLrdva = 100.3
    lcl_amp_bone <- fixed(log(10.3))
    label("Bone metabolic clearance, GS-704277 to nucleoside monophosphate (L/h)") # Table S3 mCLamp = 10.3
    lcl_mpn_bone <- fixed(log(20.6))
    label("Bone metabolic clearance, nucleoside monophosphate to GS-441524 (L/h)") # Table S3 mCLmpn = 20.6
    lcl_nmp_bone <- fixed(log(5.2))
    label("Bone metabolic clearance, GS-441524 to nucleoside monophosphate (L/h)") # Table S3 mCLnmp = 5.2
    lcl_mptn_bone <- fixed(log(100.3))
    label("Bone metabolic clearance, nucleoside monophosphate to GS-443902 (L/h)") # Table S3 mCLmptn = 100.3
    lcl_tn_bone <- fixed(log(0.31))
    label("Bone elimination clearance of GS-443902 (L/h)") # Table S3 eCLtn = 0.31

    # ---- Brain (Table S2 and Table S3 row 'Brain') ----
    lkp_brain <- fixed(log(1.8))
    label("Brain intracellular:plasma partition coefficient of remdesivir (unitless)") # Table S2 Partition Coefficient RDV = 1.8
    lkp_brain_gs441524 <- fixed(log(0.82))
    label("Brain intracellular:plasma partition coefficient of GS-441524 (unitless)") # Table S2 Partition Coefficient N = 0.82
    lclin_brain <- fixed(log(13.1))
    label("Brain extracellular-intracellular transport clearance of remdesivir (L/h)") # Table S2 Transport Clearance RDV = 13.1
    lclin_brain_gs441524 <- fixed(log(1.45))
    label("Brain extracellular-intracellular transport clearance of GS-441524 (L/h)") # Table S2 Transport Clearance N = 1.45
    lcl_rdva_brain <- fixed(log(14.5))
    label("Brain metabolic clearance, remdesivir to GS-704277 (L/h)") # Table S3 mCLrdva = 14.5
    lcl_amp_brain <- fixed(log(1.45))
    label("Brain metabolic clearance, GS-704277 to nucleoside monophosphate (L/h)") # Table S3 mCLamp = 1.45
    lcl_mpn_brain <- fixed(log(2.9))
    label("Brain metabolic clearance, nucleoside monophosphate to GS-441524 (L/h)") # Table S3 mCLmpn = 2.9
    lcl_nmp_brain <- fixed(log(0.725))
    label("Brain metabolic clearance, GS-441524 to nucleoside monophosphate (L/h)") # Table S3 mCLnmp = 0.725
    lcl_mptn_brain <- fixed(log(14.5))
    label("Brain metabolic clearance, nucleoside monophosphate to GS-443902 (L/h)") # Table S3 mCLmptn = 14.5
    lcl_tn_brain <- fixed(log(0.044))
    label("Brain elimination clearance of GS-443902 (L/h)") # Table S3 eCLtn = 0.044

    # ---- GI (Table S2 and Table S3 row 'GI') ----
    lkp_gut <- fixed(log(1.09))
    label("GI intracellular:plasma partition coefficient of remdesivir (unitless)") # Table S2 Partition Coefficient RDV = 1.09
    lkp_gut_gs441524 <- fixed(log(0.81))
    label("GI intracellular:plasma partition coefficient of GS-441524 (unitless)") # Table S2 Partition Coefficient N = 0.81
    lclin_gut <- fixed(log(10.3))
    label("GI extracellular-intracellular transport clearance of remdesivir (L/h)") # Table S2 Transport Clearance RDV = 10.3
    lclin_gut_gs441524 <- fixed(log(1.14))
    label("GI extracellular-intracellular transport clearance of GS-441524 (L/h)") # Table S2 Transport Clearance N = 1.14
    lcl_rdva_gut <- fixed(log(11.4))
    label("GI metabolic clearance, remdesivir to GS-704277 (L/h)") # Table S3 mCLrdva = 11.4
    lcl_amp_gut <- fixed(log(1.14))
    label("GI metabolic clearance, GS-704277 to nucleoside monophosphate (L/h)") # Table S3 mCLamp = 1.14
    lcl_mpn_gut <- fixed(log(2.28))
    label("GI metabolic clearance, nucleoside monophosphate to GS-441524 (L/h)") # Table S3 mCLmpn = 2.28
    lcl_nmp_gut <- fixed(log(0.57))
    label("GI metabolic clearance, GS-441524 to nucleoside monophosphate (L/h)") # Table S3 mCLnmp = 0.57
    lcl_mptn_gut <- fixed(log(11.4))
    label("GI metabolic clearance, nucleoside monophosphate to GS-443902 (L/h)") # Table S3 mCLmptn = 11.4
    lcl_tn_gut <- fixed(log(0.034))
    label("GI elimination clearance of GS-443902 (L/h)") # Table S3 eCLtn = 0.034

    # ---- Heart (Table S2 and Table S3 row 'Heart') ----
    lkp_heart <- fixed(log(2.88))
    label("Heart intracellular:plasma partition coefficient of remdesivir (unitless)") # Table S2 Partition Coefficient RDV = 2.88
    lkp_heart_gs441524 <- fixed(log(0.76))
    label("Heart intracellular:plasma partition coefficient of GS-441524 (unitless)") # Table S2 Partition Coefficient N = 0.76
    lclin_heart <- fixed(log(2.56))
    label("Heart extracellular-intracellular transport clearance of remdesivir (L/h)") # Table S2 Transport Clearance RDV = 2.56
    lclin_heart_gs441524 <- fixed(log(0.32))
    label("Heart extracellular-intracellular transport clearance of GS-441524 (L/h)") # Table S2 Transport Clearance N = 0.32
    lcl_rdva_heart <- fixed(log(3.2))
    label("Heart metabolic clearance, remdesivir to GS-704277 (L/h)") # Table S3 mCLrdva = 3.2
    lcl_amp_heart <- fixed(log(0.32))
    label("Heart metabolic clearance, GS-704277 to nucleoside monophosphate (L/h)") # Table S3 mCLamp = 0.32
    lcl_mpn_heart <- fixed(log(0.64))
    label("Heart metabolic clearance, nucleoside monophosphate to GS-441524 (L/h)") # Table S3 mCLmpn = 0.64
    lcl_nmp_heart <- fixed(log(0.16))
    label("Heart metabolic clearance, GS-441524 to nucleoside monophosphate (L/h)") # Table S3 mCLnmp = 0.16
    lcl_mptn_heart <- fixed(log(3.2))
    label("Heart metabolic clearance, nucleoside monophosphate to GS-443902 (L/h)") # Table S3 mCLmptn = 3.2
    lcl_tn_heart <- fixed(log(0.0096))
    label("Heart elimination clearance of GS-443902 (L/h)") # Table S3 eCLtn = 0.0096

    # ---- Kidney (Table S2 and Table S3 row 'Kidney') ----
    lkp_kidney <- fixed(log(0.95))
    label("Kidney intracellular:plasma partition coefficient of remdesivir (unitless)") # Table S2 Partition Coefficient RDV = 0.95
    lkp_kidney_gs441524 <- fixed(log(0.8))
    label("Kidney intracellular:plasma partition coefficient of GS-441524 (unitless)") # Table S2 Partition Coefficient N = 0.8
    lclin_kidney <- fixed(log(2.3))
    label("Kidney extracellular-intracellular transport clearance of remdesivir (L/h)") # Table S2 Transport Clearance RDV = 2.3
    lclin_kidney_gs441524 <- fixed(log(0.25))
    label("Kidney extracellular-intracellular transport clearance of GS-441524 (L/h)") # Table S2 Transport Clearance N = 0.25
    lcl_rdva_kidney <- fixed(log(2.5))
    label("Kidney metabolic clearance, remdesivir to GS-704277 (L/h)") # Table S3 mCLrdva = 2.5
    lcl_amp_kidney <- fixed(log(0.25))
    label("Kidney metabolic clearance, GS-704277 to nucleoside monophosphate (L/h)") # Table S3 mCLamp = 0.25
    lcl_mpn_kidney <- fixed(log(0.5))
    label("Kidney metabolic clearance, nucleoside monophosphate to GS-441524 (L/h)") # Table S3 mCLmpn = 0.5
    lcl_nmp_kidney <- fixed(log(0.125))
    label("Kidney metabolic clearance, GS-441524 to nucleoside monophosphate (L/h)") # Table S3 mCLnmp = 0.125
    lcl_mptn_kidney <- fixed(log(2.5))
    label("Kidney metabolic clearance, nucleoside monophosphate to GS-443902 (L/h)") # Table S3 mCLmptn = 2.5
    lcl_tn_kidney <- fixed(log(0.0075))
    label("Kidney elimination clearance of GS-443902 (L/h)") # Table S3 eCLtn = 0.0075

    # ---- Liver (Table S2 and Table S3 row 'Liver') ----
    lkp_liver <- fixed(log(1.21))
    label("Liver intracellular:plasma partition coefficient of remdesivir (unitless)") # Table S2 Partition Coefficient RDV = 1.21
    lkp_liver_gs441524 <- fixed(log(0.78))
    label("Liver intracellular:plasma partition coefficient of GS-441524 (unitless)") # Table S2 Partition Coefficient N = 0.78
    lclin_liver <- fixed(log(14.3))
    label("Liver extracellular-intracellular transport clearance of remdesivir (L/h)") # Table S2 Transport Clearance RDV = 14.3
    lclin_liver_gs441524 <- fixed(log(1.59))
    label("Liver extracellular-intracellular transport clearance of GS-441524 (L/h)") # Table S2 Transport Clearance N = 1.59
    lcl_rdva_liver <- fixed(log(15.9))
    label("Liver metabolic clearance, remdesivir to GS-704277 (L/h)") # Table S3 mCLrdva = 15.9
    lcl_amp_liver <- fixed(log(1.59))
    label("Liver metabolic clearance, GS-704277 to nucleoside monophosphate (L/h)") # Table S3 mCLamp = 1.59
    lcl_mpn_liver <- fixed(log(3.18))
    label("Liver metabolic clearance, nucleoside monophosphate to GS-441524 (L/h)") # Table S3 mCLmpn = 3.18
    lcl_nmp_liver <- fixed(log(0.795))
    label("Liver metabolic clearance, GS-441524 to nucleoside monophosphate (L/h)") # Table S3 mCLnmp = 0.795
    lcl_mptn_liver <- fixed(log(15.9))
    label("Liver metabolic clearance, nucleoside monophosphate to GS-443902 (L/h)") # Table S3 mCLmptn = 15.9
    lcl_tn_liver <- fixed(log(0.048))
    label("Liver elimination clearance of GS-443902 (L/h)") # Table S3 eCLtn = 0.048

    # ---- Lung (Table S2 and Table S3 row 'Lung') ----
    lkp_lung <- fixed(log(0.32))
    label("Lung intracellular:plasma partition coefficient of remdesivir (unitless)") # Table S2 Partition Coefficient RDV = 0.32
    lkp_lung_gs441524 <- fixed(log(0.84))
    label("Lung intracellular:plasma partition coefficient of GS-441524 (unitless)") # Table S2 Partition Coefficient N = 0.84
    lclin_lung <- fixed(log(2.5))
    label("Lung extracellular-intracellular transport clearance of remdesivir (L/h)") # Table S2 Transport Clearance RDV = 2.5
    lclin_lung_gs441524 <- fixed(log(0.28))
    label("Lung extracellular-intracellular transport clearance of GS-441524 (L/h)") # Table S2 Transport Clearance N = 0.28
    lcl_rdva_lung <- fixed(log(2.8))
    label("Lung metabolic clearance, remdesivir to GS-704277 (L/h)") # Table S3 mCLrdva = 2.8
    lcl_amp_lung <- fixed(log(0.28))
    label("Lung metabolic clearance, GS-704277 to nucleoside monophosphate (L/h)") # Table S3 mCLamp = 0.28
    lcl_mpn_lung <- fixed(log(0.56))
    label("Lung metabolic clearance, nucleoside monophosphate to GS-441524 (L/h)") # Table S3 mCLmpn = 0.56
    lcl_nmp_lung <- fixed(log(0.14))
    label("Lung metabolic clearance, GS-441524 to nucleoside monophosphate (L/h)") # Table S3 mCLnmp = 0.14
    lcl_mptn_lung <- fixed(log(2.8))
    label("Lung metabolic clearance, nucleoside monophosphate to GS-443902 (L/h)") # Table S3 mCLmptn = 2.8
    lcl_tn_lung <- fixed(log(0.0084))
    label("Lung elimination clearance of GS-443902 (L/h)") # Table S3 eCLtn = 0.0084

    # ---- Muscle (Table S2 and Table S3 row 'Muscle') ----
    lkp_muscle <- fixed(log(1.35))
    label("Muscle intracellular:plasma partition coefficient of remdesivir (unitless)") # Table S2 Partition Coefficient RDV = 1.35
    lkp_muscle_gs441524 <- fixed(log(0.84))
    label("Muscle intracellular:plasma partition coefficient of GS-441524 (unitless)") # Table S2 Partition Coefficient N = 0.84
    lclin_muscle <- fixed(log(246.2))
    label("Muscle extracellular-intracellular transport clearance of remdesivir (L/h)") # Table S2 Transport Clearance RDV = 246.2
    lclin_muscle_gs441524 <- fixed(log(27.36))
    label("Muscle extracellular-intracellular transport clearance of GS-441524 (L/h)") # Table S2 Transport Clearance N = 27.36
    lcl_rdva_muscle <- fixed(log(273.6))
    label("Muscle metabolic clearance, remdesivir to GS-704277 (L/h)") # Table S3 mCLrdva = 273.6
    lcl_amp_muscle <- fixed(log(27.36))
    label("Muscle metabolic clearance, GS-704277 to nucleoside monophosphate (L/h)") # Table S3 mCLamp = 27.36
    lcl_mpn_muscle <- fixed(log(54.72))
    label("Muscle metabolic clearance, nucleoside monophosphate to GS-441524 (L/h)") # Table S3 mCLmpn = 54.72
    lcl_nmp_muscle <- fixed(log(13.68))
    label("Muscle metabolic clearance, GS-441524 to nucleoside monophosphate (L/h)") # Table S3 mCLnmp = 13.68
    lcl_mptn_muscle <- fixed(log(273.6))
    label("Muscle metabolic clearance, nucleoside monophosphate to GS-443902 (L/h)") # Table S3 mCLmptn = 273.6
    lcl_tn_muscle <- fixed(log(0.82))
    label("Muscle elimination clearance of GS-443902 (L/h)") # Table S3 eCLtn = 0.82

    # ---- Rest of Body (Table S2 and Table S3 row 'Rest of Body') ----
    lkp_remainder <- fixed(log(1.9))
    label("Rest of Body intracellular:plasma partition coefficient of remdesivir (unitless)") # Table S2 Partition Coefficient RDV = 1.9
    lkp_remainder_gs441524 <- fixed(log(0.74))
    label("Rest of Body intracellular:plasma partition coefficient of GS-441524 (unitless)") # Table S2 Partition Coefficient N = 0.74
    lclin_remainder <- fixed(log(9))
    label("Rest of Body extracellular-intracellular transport clearance of remdesivir (L/h)") # Table S2 Transport Clearance RDV = 9
    lclin_remainder_gs441524 <- fixed(log(1))
    label("Rest of Body extracellular-intracellular transport clearance of GS-441524 (L/h)") # Table S2 Transport Clearance N = 1
    lcl_rdva_remainder <- fixed(log(10))
    label("Rest of Body metabolic clearance, remdesivir to GS-704277 (L/h)") # Table S3 mCLrdva = 10
    lcl_amp_remainder <- fixed(log(1))
    label("Rest of Body metabolic clearance, GS-704277 to nucleoside monophosphate (L/h)") # Table S3 mCLamp = 1
    lcl_mpn_remainder <- fixed(log(2))
    label("Rest of Body metabolic clearance, nucleoside monophosphate to GS-441524 (L/h)") # Table S3 mCLmpn = 2
    lcl_nmp_remainder <- fixed(log(0.5))
    label("Rest of Body metabolic clearance, GS-441524 to nucleoside monophosphate (L/h)") # Table S3 mCLnmp = 0.5
    lcl_mptn_remainder <- fixed(log(10))
    label("Rest of Body metabolic clearance, nucleoside monophosphate to GS-443902 (L/h)") # Table S3 mCLmptn = 10
    lcl_tn_remainder <- fixed(log(0.03))
    label("Rest of Body elimination clearance of GS-443902 (L/h)") # Table S3 eCLtn = 0.03

    # ---- Skin (Table S2 and Table S3 row 'Skin') ----
    lkp_skin <- fixed(log(1.7))
    label("Skin intracellular:plasma partition coefficient of remdesivir (unitless)") # Table S2 Partition Coefficient RDV = 1.7
    lkp_skin_gs441524 <- fixed(log(0.81))
    label("Skin intracellular:plasma partition coefficient of GS-441524 (unitless)") # Table S2 Partition Coefficient N = 0.81
    lclin_skin <- fixed(log(22.2))
    label("Skin extracellular-intracellular transport clearance of remdesivir (L/h)") # Table S2 Transport Clearance RDV = 22.2
    lclin_skin_gs441524 <- fixed(log(2.47))
    label("Skin extracellular-intracellular transport clearance of GS-441524 (L/h)") # Table S2 Transport Clearance N = 2.47
    lcl_rdva_skin <- fixed(log(24.7))
    label("Skin metabolic clearance, remdesivir to GS-704277 (L/h)") # Table S3 mCLrdva = 24.7
    lcl_amp_skin <- fixed(log(2.47))
    label("Skin metabolic clearance, GS-704277 to nucleoside monophosphate (L/h)") # Table S3 mCLamp = 2.47
    lcl_mpn_skin <- fixed(log(4.94))
    label("Skin metabolic clearance, nucleoside monophosphate to GS-441524 (L/h)") # Table S3 mCLmpn = 4.94
    lcl_nmp_skin <- fixed(log(1.24))
    label("Skin metabolic clearance, GS-441524 to nucleoside monophosphate (L/h)") # Table S3 mCLnmp = 1.24
    lcl_mptn_skin <- fixed(log(24.7))
    label("Skin metabolic clearance, nucleoside monophosphate to GS-443902 (L/h)") # Table S3 mCLmptn = 24.7
    lcl_tn_skin <- fixed(log(0.074))
    label("Skin elimination clearance of GS-443902 (L/h)") # Table S3 eCLtn = 0.074

    # ---- Spleen (Table S2 and Table S3 row 'Spleen') ----
    lkp_spleen <- fixed(log(0.41))
    label("Spleen intracellular:plasma partition coefficient of remdesivir (unitless)") # Table S2 Partition Coefficient RDV = 0.41
    lkp_spleen_gs441524 <- fixed(log(0.81))
    label("Spleen intracellular:plasma partition coefficient of GS-441524 (unitless)") # Table S2 Partition Coefficient N = 0.81
    lclin_spleen <- fixed(log(1))
    label("Spleen extracellular-intracellular transport clearance of remdesivir (L/h)") # Table S2 Transport Clearance RDV = 1
    lclin_spleen_gs441524 <- fixed(log(0.11))
    label("Spleen extracellular-intracellular transport clearance of GS-441524 (L/h)") # Table S2 Transport Clearance N = 0.11
    lcl_rdva_spleen <- fixed(log(1.1))
    label("Spleen metabolic clearance, remdesivir to GS-704277 (L/h)") # Table S3 mCLrdva = 1.1
    lcl_amp_spleen <- fixed(log(0.11))
    label("Spleen metabolic clearance, GS-704277 to nucleoside monophosphate (L/h)") # Table S3 mCLamp = 0.11
    lcl_mpn_spleen <- fixed(log(0.22))
    label("Spleen metabolic clearance, nucleoside monophosphate to GS-441524 (L/h)") # Table S3 mCLmpn = 0.22
    lcl_nmp_spleen <- fixed(log(0.055))
    label("Spleen metabolic clearance, GS-441524 to nucleoside monophosphate (L/h)") # Table S3 mCLnmp = 0.055
    lcl_mptn_spleen <- fixed(log(1.1))
    label("Spleen metabolic clearance, nucleoside monophosphate to GS-443902 (L/h)") # Table S3 mCLmptn = 1.1
    lcl_tn_spleen <- fixed(log(0.0033))
    label("Spleen elimination clearance of GS-443902 (L/h)") # Table S3 eCLtn = 0.0033

    # =================================================================
    # MONTE-CARLO VARIABILITY -- the paper does not estimate
    # between-subject variability; every figure is a 500-replicate
    # Monte-Carlo simulation with a 20% CV on a named parameter subset
    # (Methods 'Model performance'; Results; Figure 3-5 captions):
    #   Figure 3, S1, S2 (plasma): all eight forcing-function parameters;
    #   Figure 4 (PBMC): tKrdvpbmc and eKtnpbmc;
    #   Figure 5, S5 (lung): tCLrdvlu and eCLtnluic.
    # The distribution is not stated; it is encoded here as log-normal,
    # omega^2 = log(1 + 0.2^2) = 0.0392207. The subsets were applied one
    # figure at a time, so simulating all twelve etas together is broader
    # than any single published figure (see the vignette).
    # =================================================================
    etalcl ~ fixed(0.0392207)
    etalvc ~ fixed(0.0392207)
    etalk12 ~ fixed(0.0392207)
    etalk21 ~ fixed(0.0392207)
    etalk_gs704277_form ~ fixed(0.0392207)
    etalkel_gs704277 ~ fixed(0.0392207)
    etalk_gs441524_form ~ fixed(0.0392207)
    etalkel_gs441524 ~ fixed(0.0392207)
    etalkin_pbmc ~ fixed(0.0392207)
    etalk_tn_pbmc ~ fixed(0.0392207)
    etalclin_lung ~ fixed(0.0392207)
    etalcl_tn_lung ~ fixed(0.0392207)

    # Residual error: the forcing function was fitted with an additive
    # error model (Methods) but no magnitude is reported, so the additive
    # SDs are fixed to zero.
    addSd <- fixed(0)
    label("Additive residual error, plasma remdesivir (mg/L; not reported)")
    addSd_gs704277 <- fixed(0)
    label("Additive residual error, plasma GS-704277 (mg/L; not reported)")
    addSd_gs441524 <- fixed(0)
    label("Additive residual error, plasma GS-441524 (mg/L; not reported)")
  })
  model({
    # ---------------- Forcing function (Table S1) ----------------
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    kel <- cl / vc
    k_gs704277_form <- exp(lk_gs704277_form + etalk_gs704277_form)
    kel_gs704277 <- exp(lkel_gs704277 + etalkel_gs704277)
    k_gs441524_form <- exp(lk_gs441524_form + etalk_gs441524_form)
    kel_gs441524 <- exp(lkel_gs441524 + etalkel_gs441524)

    # ---------------- PBMC module (Tables S2, S3) ----------------
    kin_pbmc <- exp(lkin_pbmc + etalkin_pbmc)
    kin_pbmc_gs441524 <- exp(lkin_pbmc_gs441524)
    kp_pbmc <- exp(lkp_pbmc)
    kp_pbmc_gs441524 <- exp(lkp_pbmc_gs441524)
    k_rdva_pbmc <- exp(lk_rdva_pbmc)
    k_amp_pbmc <- exp(lk_amp_pbmc)
    k_mpn_pbmc <- exp(lk_mpn_pbmc)
    k_nmp_pbmc <- exp(lk_nmp_pbmc)
    k_mptn_pbmc <- exp(lk_mptn_pbmc)
    k_tn_pbmc <- exp(lk_tn_pbmc + etalk_tn_pbmc)

    # ---------------- Physiology and tissue parameters ----------------
    # Plasma flows (L/h), Table S2 'Plasma flow rate' column (PK-Sim, 80 kg
    # White male, blood flows converted to plasma with hematocrit 0.45). The
    # lung receives the whole cardiac output (Table S2 footnote '#').
    q_adipose <- 12.42 # Table S2 Adipose
    q_bone <- 8.91 # Table S2 Bone
    q_brain <- 21.06 # Table S2 Brain
    q_gut <- 26.19 # Table S2 GI
    q_heart <- 7.02 # Table S2 Heart
    q_kidney <- 35.64 # Table S2 Kidney
    q_liver <- 11.34 # Table S2 Liver
    q_muscle <- 30.78 # Table S2 Muscle
    q_remainder <- 1.89 # Table S2 Rest of Body
    q_skin <- 8.91 # Table S2 Skin
    q_spleen <- 4.59 # Table S2 Spleen
    qc <- 168.75 # Table S2 Lung / Arterial Plasma / Venous Plasma 168.75 (cardiac output)

    # Extracellular (is_) and intracellular (int_) volumes (L), Table S2.
    v_is_adipose <- 3.62 # Table S2 Adipose Extracellular
    v_int_adipose <- 17.38 # Table S2 Adipose Intracellular
    v_is_bone <- 1.38 # Table S2 Bone Extracellular
    v_int_bone <- 10.34 # Table S2 Bone Intracellular
    v_is_brain <- 0.039 # Table S2 Brain Extracellular
    v_int_brain <- 1.45 # Table S2 Brain Intracellular
    v_is_gut <- 0.144 # Table S2 GI Extracellular
    v_int_gut <- 1.14 # Table S2 GI Intracellular
    v_is_heart <- 0.075 # Table S2 Heart Extracellular
    v_int_heart <- 0.32 # Table S2 Heart Intracellular
    v_is_kidney <- 0.143 # Table S2 Kidney Extracellular
    v_int_kidney <- 0.25 # Table S2 Kidney Intracellular
    v_is_liver <- 0.6 # Table S2 Liver Extracellular
    v_int_liver <- 1.59 # Table S2 Liver Intracellular
    v_is_lung <- 0.615 # Table S2 Lung Extracellular
    v_int_lung <- 0.28 # Table S2 Lung Intracellular
    v_is_muscle <- 5.89 # Table S2 Muscle Extracellular
    v_int_muscle <- 27.36 # Table S2 Muscle Intracellular
    v_is_remainder <- 0.22 # Table S2 Rest of Body Extracellular
    v_int_remainder <- 1 # Table S2 Rest of Body Intracellular
    v_is_skin <- 1.24 # Table S2 Skin Extracellular
    v_int_skin <- 2.47 # Table S2 Skin Intracellular
    v_is_spleen <- 0.042 # Table S2 Spleen Extracellular 0.042 (deposited code: 0.07; see vignette)
    v_int_spleen <- 0.11 # Table S2 Spleen Intracellular
    v_arterial <- 0.23 # Table S2 Arterial Plasma volume (footnote '^')

    # Back-transformed tissue parameters
    kp_adipose <- exp(lkp_adipose)
    kp_adipose_gs441524 <- exp(lkp_adipose_gs441524)
    clin_adipose <- exp(lclin_adipose)
    clin_adipose_gs441524 <- exp(lclin_adipose_gs441524)
    cl_rdva_adipose <- exp(lcl_rdva_adipose)
    cl_amp_adipose <- exp(lcl_amp_adipose)
    cl_mpn_adipose <- exp(lcl_mpn_adipose)
    cl_nmp_adipose <- exp(lcl_nmp_adipose)
    cl_mptn_adipose <- exp(lcl_mptn_adipose)
    cl_tn_adipose <- exp(lcl_tn_adipose)
    kp_bone <- exp(lkp_bone)
    kp_bone_gs441524 <- exp(lkp_bone_gs441524)
    clin_bone <- exp(lclin_bone)
    clin_bone_gs441524 <- exp(lclin_bone_gs441524)
    cl_rdva_bone <- exp(lcl_rdva_bone)
    cl_amp_bone <- exp(lcl_amp_bone)
    cl_mpn_bone <- exp(lcl_mpn_bone)
    cl_nmp_bone <- exp(lcl_nmp_bone)
    cl_mptn_bone <- exp(lcl_mptn_bone)
    cl_tn_bone <- exp(lcl_tn_bone)
    kp_brain <- exp(lkp_brain)
    kp_brain_gs441524 <- exp(lkp_brain_gs441524)
    clin_brain <- exp(lclin_brain)
    clin_brain_gs441524 <- exp(lclin_brain_gs441524)
    cl_rdva_brain <- exp(lcl_rdva_brain)
    cl_amp_brain <- exp(lcl_amp_brain)
    cl_mpn_brain <- exp(lcl_mpn_brain)
    cl_nmp_brain <- exp(lcl_nmp_brain)
    cl_mptn_brain <- exp(lcl_mptn_brain)
    cl_tn_brain <- exp(lcl_tn_brain)
    kp_gut <- exp(lkp_gut)
    kp_gut_gs441524 <- exp(lkp_gut_gs441524)
    clin_gut <- exp(lclin_gut)
    clin_gut_gs441524 <- exp(lclin_gut_gs441524)
    cl_rdva_gut <- exp(lcl_rdva_gut)
    cl_amp_gut <- exp(lcl_amp_gut)
    cl_mpn_gut <- exp(lcl_mpn_gut)
    cl_nmp_gut <- exp(lcl_nmp_gut)
    cl_mptn_gut <- exp(lcl_mptn_gut)
    cl_tn_gut <- exp(lcl_tn_gut)
    kp_heart <- exp(lkp_heart)
    kp_heart_gs441524 <- exp(lkp_heart_gs441524)
    clin_heart <- exp(lclin_heart)
    clin_heart_gs441524 <- exp(lclin_heart_gs441524)
    cl_rdva_heart <- exp(lcl_rdva_heart)
    cl_amp_heart <- exp(lcl_amp_heart)
    cl_mpn_heart <- exp(lcl_mpn_heart)
    cl_nmp_heart <- exp(lcl_nmp_heart)
    cl_mptn_heart <- exp(lcl_mptn_heart)
    cl_tn_heart <- exp(lcl_tn_heart)
    kp_kidney <- exp(lkp_kidney)
    kp_kidney_gs441524 <- exp(lkp_kidney_gs441524)
    clin_kidney <- exp(lclin_kidney)
    clin_kidney_gs441524 <- exp(lclin_kidney_gs441524)
    cl_rdva_kidney <- exp(lcl_rdva_kidney)
    cl_amp_kidney <- exp(lcl_amp_kidney)
    cl_mpn_kidney <- exp(lcl_mpn_kidney)
    cl_nmp_kidney <- exp(lcl_nmp_kidney)
    cl_mptn_kidney <- exp(lcl_mptn_kidney)
    cl_tn_kidney <- exp(lcl_tn_kidney)
    kp_liver <- exp(lkp_liver)
    kp_liver_gs441524 <- exp(lkp_liver_gs441524)
    clin_liver <- exp(lclin_liver)
    clin_liver_gs441524 <- exp(lclin_liver_gs441524)
    cl_rdva_liver <- exp(lcl_rdva_liver)
    cl_amp_liver <- exp(lcl_amp_liver)
    cl_mpn_liver <- exp(lcl_mpn_liver)
    cl_nmp_liver <- exp(lcl_nmp_liver)
    cl_mptn_liver <- exp(lcl_mptn_liver)
    cl_tn_liver <- exp(lcl_tn_liver)
    kp_lung <- exp(lkp_lung)
    kp_lung_gs441524 <- exp(lkp_lung_gs441524)
    clin_lung <- exp(lclin_lung + etalclin_lung)
    clin_lung_gs441524 <- exp(lclin_lung_gs441524)
    cl_rdva_lung <- exp(lcl_rdva_lung)
    cl_amp_lung <- exp(lcl_amp_lung)
    cl_mpn_lung <- exp(lcl_mpn_lung)
    cl_nmp_lung <- exp(lcl_nmp_lung)
    cl_mptn_lung <- exp(lcl_mptn_lung)
    cl_tn_lung <- exp(lcl_tn_lung + etalcl_tn_lung)
    kp_muscle <- exp(lkp_muscle)
    kp_muscle_gs441524 <- exp(lkp_muscle_gs441524)
    clin_muscle <- exp(lclin_muscle)
    clin_muscle_gs441524 <- exp(lclin_muscle_gs441524)
    cl_rdva_muscle <- exp(lcl_rdva_muscle)
    cl_amp_muscle <- exp(lcl_amp_muscle)
    cl_mpn_muscle <- exp(lcl_mpn_muscle)
    cl_nmp_muscle <- exp(lcl_nmp_muscle)
    cl_mptn_muscle <- exp(lcl_mptn_muscle)
    cl_tn_muscle <- exp(lcl_tn_muscle)
    kp_remainder <- exp(lkp_remainder)
    kp_remainder_gs441524 <- exp(lkp_remainder_gs441524)
    clin_remainder <- exp(lclin_remainder)
    clin_remainder_gs441524 <- exp(lclin_remainder_gs441524)
    cl_rdva_remainder <- exp(lcl_rdva_remainder)
    cl_amp_remainder <- exp(lcl_amp_remainder)
    cl_mpn_remainder <- exp(lcl_mpn_remainder)
    cl_nmp_remainder <- exp(lcl_nmp_remainder)
    cl_mptn_remainder <- exp(lcl_mptn_remainder)
    cl_tn_remainder <- exp(lcl_tn_remainder)
    kp_skin <- exp(lkp_skin)
    kp_skin_gs441524 <- exp(lkp_skin_gs441524)
    clin_skin <- exp(lclin_skin)
    clin_skin_gs441524 <- exp(lclin_skin_gs441524)
    cl_rdva_skin <- exp(lcl_rdva_skin)
    cl_amp_skin <- exp(lcl_amp_skin)
    cl_mpn_skin <- exp(lcl_mpn_skin)
    cl_nmp_skin <- exp(lcl_nmp_skin)
    cl_mptn_skin <- exp(lcl_mptn_skin)
    cl_tn_skin <- exp(lcl_tn_skin)
    kp_spleen <- exp(lkp_spleen)
    kp_spleen_gs441524 <- exp(lkp_spleen_gs441524)
    clin_spleen <- exp(lclin_spleen)
    clin_spleen_gs441524 <- exp(lclin_spleen_gs441524)
    cl_rdva_spleen <- exp(lcl_rdva_spleen)
    cl_amp_spleen <- exp(lcl_amp_spleen)
    cl_mpn_spleen <- exp(lcl_mpn_spleen)
    cl_nmp_spleen <- exp(lcl_nmp_spleen)
    cl_mptn_spleen <- exp(lcl_mptn_spleen)
    cl_tn_spleen <- exp(lcl_tn_spleen)

    # Tissue concentrations (mg/L)
    Cis_adipose <- is_adipose / v_is_adipose
    Cint_adipose <- int_adipose / v_int_adipose
    Cint_adipose_gs704277 <- int_adipose_gs704277 / v_int_adipose
    Cint_adipose_gs441524mp <- int_adipose_gs441524mp / v_int_adipose
    Cis_adipose_gs441524 <- is_adipose_gs441524 / v_is_adipose
    Cint_adipose_gs441524 <- int_adipose_gs441524 / v_int_adipose
    Cint_adipose_gs443902 <- int_adipose_gs443902 / v_int_adipose
    Cis_bone <- is_bone / v_is_bone
    Cint_bone <- int_bone / v_int_bone
    Cint_bone_gs704277 <- int_bone_gs704277 / v_int_bone
    Cint_bone_gs441524mp <- int_bone_gs441524mp / v_int_bone
    Cis_bone_gs441524 <- is_bone_gs441524 / v_is_bone
    Cint_bone_gs441524 <- int_bone_gs441524 / v_int_bone
    Cint_bone_gs443902 <- int_bone_gs443902 / v_int_bone
    Cis_brain <- is_brain / v_is_brain
    Cint_brain <- int_brain / v_int_brain
    Cint_brain_gs704277 <- int_brain_gs704277 / v_int_brain
    Cint_brain_gs441524mp <- int_brain_gs441524mp / v_int_brain
    Cis_brain_gs441524 <- is_brain_gs441524 / v_is_brain
    Cint_brain_gs441524 <- int_brain_gs441524 / v_int_brain
    Cint_brain_gs443902 <- int_brain_gs443902 / v_int_brain
    Cis_gut <- is_gut / v_is_gut
    Cint_gut <- int_gut / v_int_gut
    Cint_gut_gs704277 <- int_gut_gs704277 / v_int_gut
    Cint_gut_gs441524mp <- int_gut_gs441524mp / v_int_gut
    Cis_gut_gs441524 <- is_gut_gs441524 / v_is_gut
    Cint_gut_gs441524 <- int_gut_gs441524 / v_int_gut
    Cint_gut_gs443902 <- int_gut_gs443902 / v_int_gut
    Cis_heart <- is_heart / v_is_heart
    Cint_heart <- int_heart / v_int_heart
    Cint_heart_gs704277 <- int_heart_gs704277 / v_int_heart
    Cint_heart_gs441524mp <- int_heart_gs441524mp / v_int_heart
    Cis_heart_gs441524 <- is_heart_gs441524 / v_is_heart
    Cint_heart_gs441524 <- int_heart_gs441524 / v_int_heart
    Cint_heart_gs443902 <- int_heart_gs443902 / v_int_heart
    Cis_kidney <- is_kidney / v_is_kidney
    Cint_kidney <- int_kidney / v_int_kidney
    Cint_kidney_gs704277 <- int_kidney_gs704277 / v_int_kidney
    Cint_kidney_gs441524mp <- int_kidney_gs441524mp / v_int_kidney
    Cis_kidney_gs441524 <- is_kidney_gs441524 / v_is_kidney
    Cint_kidney_gs441524 <- int_kidney_gs441524 / v_int_kidney
    Cint_kidney_gs443902 <- int_kidney_gs443902 / v_int_kidney
    Cis_liver <- is_liver / v_is_liver
    Cint_liver <- int_liver / v_int_liver
    Cint_liver_gs704277 <- int_liver_gs704277 / v_int_liver
    Cint_liver_gs441524mp <- int_liver_gs441524mp / v_int_liver
    Cis_liver_gs441524 <- is_liver_gs441524 / v_is_liver
    Cint_liver_gs441524 <- int_liver_gs441524 / v_int_liver
    Cint_liver_gs443902 <- int_liver_gs443902 / v_int_liver
    Cis_lung <- is_lung / v_is_lung
    Cint_lung <- int_lung / v_int_lung
    Cint_lung_gs704277 <- int_lung_gs704277 / v_int_lung
    Cint_lung_gs441524mp <- int_lung_gs441524mp / v_int_lung
    Cis_lung_gs441524 <- is_lung_gs441524 / v_is_lung
    Cint_lung_gs441524 <- int_lung_gs441524 / v_int_lung
    Cint_lung_gs443902 <- int_lung_gs443902 / v_int_lung
    Cis_muscle <- is_muscle / v_is_muscle
    Cint_muscle <- int_muscle / v_int_muscle
    Cint_muscle_gs704277 <- int_muscle_gs704277 / v_int_muscle
    Cint_muscle_gs441524mp <- int_muscle_gs441524mp / v_int_muscle
    Cis_muscle_gs441524 <- is_muscle_gs441524 / v_is_muscle
    Cint_muscle_gs441524 <- int_muscle_gs441524 / v_int_muscle
    Cint_muscle_gs443902 <- int_muscle_gs443902 / v_int_muscle
    Cis_remainder <- is_remainder / v_is_remainder
    Cint_remainder <- int_remainder / v_int_remainder
    Cint_remainder_gs704277 <- int_remainder_gs704277 / v_int_remainder
    Cint_remainder_gs441524mp <- int_remainder_gs441524mp / v_int_remainder
    Cis_remainder_gs441524 <- is_remainder_gs441524 / v_is_remainder
    Cint_remainder_gs441524 <- int_remainder_gs441524 / v_int_remainder
    Cint_remainder_gs443902 <- int_remainder_gs443902 / v_int_remainder
    Cis_skin <- is_skin / v_is_skin
    Cint_skin <- int_skin / v_int_skin
    Cint_skin_gs704277 <- int_skin_gs704277 / v_int_skin
    Cint_skin_gs441524mp <- int_skin_gs441524mp / v_int_skin
    Cis_skin_gs441524 <- is_skin_gs441524 / v_is_skin
    Cint_skin_gs441524 <- int_skin_gs441524 / v_int_skin
    Cint_skin_gs443902 <- int_skin_gs443902 / v_int_skin
    Cis_spleen <- is_spleen / v_is_spleen
    Cint_spleen <- int_spleen / v_int_spleen
    Cint_spleen_gs704277 <- int_spleen_gs704277 / v_int_spleen
    Cint_spleen_gs441524mp <- int_spleen_gs441524mp / v_int_spleen
    Cis_spleen_gs441524 <- is_spleen_gs441524 / v_is_spleen
    Cint_spleen_gs441524 <- int_spleen_gs441524 / v_int_spleen
    Cint_spleen_gs443902 <- int_spleen_gs443902 / v_int_spleen
    Cis_lung_gs704277 <- is_lung_gs704277 / v_is_lung

    # ---------------- Plasma and PBMC concentrations ----------------
    # Venous plasma remdesivir (mg/L). central_gs704277 and
    # central_gs441524 are concentration states (mg/L): the forcing
    # function writes both metabolites directly in concentration, with
    # no volume (deposited code, 'A model' and 'N model').
    Cc <- central / vc
    Cc_gs704277 <- central_gs704277
    Cc_gs441524 <- central_gs441524
    Carterial <- arterial / v_arterial
    Carterial_gs704277 <- arterial_gs704277 / v_arterial
    Carterial_gs441524 <- arterial_gs441524 / v_arterial

    # ---------------- ODEs: forcing function (venous plasma) ----------------
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_gs704277) <- k_gs704277_form * Cc - kel_gs704277 * central_gs704277
    d/dt(central_gs441524) <- k_gs441524_form * central_gs704277 - kel_gs441524 * central_gs441524

    # ---------------- ODEs: arterial plasma (amounts) ----------------
    # Arterial plasma is fed by the lung extracellular outflow at the
    # cardiac output.
    d/dt(arterial) <- qc * Cis_lung - qc * Carterial
    d/dt(arterial_gs704277) <- qc * Cis_lung_gs704277 - qc * Carterial_gs704277
    d/dt(arterial_gs441524) <- qc * Cis_lung_gs441524 - qc * Carterial_gs441524

    # ---------------- ODEs: PBMC (concentrations, mg/L) ----------------
    # Driven by venous plasma; uptake does not deplete the forcing
    # function. Only remdesivir and GS-441524 cross the PBMC membrane.
    d/dt(pbmc) <- kin_pbmc * fu_p * (Cc - pbmc / kp_pbmc) - k_rdva_pbmc * pbmc
    d/dt(pbmc_gs704277) <- k_rdva_pbmc * pbmc - k_amp_pbmc * pbmc_gs704277
    d/dt(pbmc_gs441524mp) <- k_amp_pbmc * pbmc_gs704277 + k_nmp_pbmc * pbmc_gs441524 - (k_mpn_pbmc + k_mptn_pbmc) * pbmc_gs441524mp
    d/dt(pbmc_gs441524) <- kin_pbmc_gs441524 * fu_p_gs441524 * (Cc_gs441524 - pbmc_gs441524 / kp_pbmc_gs441524) + k_mpn_pbmc * pbmc_gs441524mp - k_nmp_pbmc * pbmc_gs441524
    d/dt(pbmc_gs443902) <- k_mptn_pbmc * pbmc_gs441524mp - k_tn_pbmc * pbmc_gs443902

    # ---------------- ODEs: tissues ----------------
    # Extracellular (is_<tissue>) and intracellular (int_<tissue>)
    # amounts (mg). Membrane flux of RDV and GS-441524 is
    # clin * fu * (C_extracellular - C_intracellular / kp) (Methods,
    # transport equation). Metabolic fluxes are clearance times
    # intracellular concentration.

    # ---- Adipose ----
    d/dt(is_adipose) <- q_adipose * Carterial - q_adipose * Cis_adipose - clin_adipose * fu_p * (Cis_adipose - Cint_adipose / kp_adipose)
    d/dt(int_adipose) <- clin_adipose * fu_p * (Cis_adipose - Cint_adipose / kp_adipose) - cl_rdva_adipose * Cint_adipose
    d/dt(int_adipose_gs704277) <- cl_rdva_adipose * Cint_adipose - cl_amp_adipose * Cint_adipose_gs704277
    d/dt(int_adipose_gs441524mp) <- cl_amp_adipose * Cint_adipose_gs704277 + cl_nmp_adipose * Cint_adipose_gs441524 - (cl_mpn_adipose + cl_mptn_adipose) * Cint_adipose_gs441524mp
    d/dt(is_adipose_gs441524) <- q_adipose * Carterial_gs441524 - q_adipose * Cis_adipose_gs441524 - clin_adipose_gs441524 * fu_p_gs441524 * (Cis_adipose_gs441524 - Cint_adipose_gs441524 / kp_adipose_gs441524)
    d/dt(int_adipose_gs441524) <- clin_adipose_gs441524 * fu_p_gs441524 * (Cis_adipose_gs441524 - Cint_adipose_gs441524 / kp_adipose_gs441524) + cl_mpn_adipose * Cint_adipose_gs441524mp - cl_nmp_adipose * Cint_adipose_gs441524
    d/dt(int_adipose_gs443902) <- cl_mptn_adipose * Cint_adipose_gs441524mp - cl_tn_adipose * Cint_adipose_gs443902

    # ---- Bone ----
    d/dt(is_bone) <- q_bone * Carterial - q_bone * Cis_bone - clin_bone * fu_p * (Cis_bone - Cint_bone / kp_bone)
    d/dt(int_bone) <- clin_bone * fu_p * (Cis_bone - Cint_bone / kp_bone) - cl_rdva_bone * Cint_bone
    d/dt(int_bone_gs704277) <- cl_rdva_bone * Cint_bone - cl_amp_bone * Cint_bone_gs704277
    d/dt(int_bone_gs441524mp) <- cl_amp_bone * Cint_bone_gs704277 + cl_nmp_bone * Cint_bone_gs441524 - (cl_mpn_bone + cl_mptn_bone) * Cint_bone_gs441524mp
    d/dt(is_bone_gs441524) <- q_bone * Carterial_gs441524 - q_bone * Cis_bone_gs441524 - clin_bone_gs441524 * fu_p_gs441524 * (Cis_bone_gs441524 - Cint_bone_gs441524 / kp_bone_gs441524)
    d/dt(int_bone_gs441524) <- clin_bone_gs441524 * fu_p_gs441524 * (Cis_bone_gs441524 - Cint_bone_gs441524 / kp_bone_gs441524) + cl_mpn_bone * Cint_bone_gs441524mp - cl_nmp_bone * Cint_bone_gs441524
    d/dt(int_bone_gs443902) <- cl_mptn_bone * Cint_bone_gs441524mp - cl_tn_bone * Cint_bone_gs443902

    # ---- Brain ----
    d/dt(is_brain) <- q_brain * Carterial - q_brain * Cis_brain - clin_brain * fu_p * (Cis_brain - Cint_brain / kp_brain)
    d/dt(int_brain) <- clin_brain * fu_p * (Cis_brain - Cint_brain / kp_brain) - cl_rdva_brain * Cint_brain
    d/dt(int_brain_gs704277) <- cl_rdva_brain * Cint_brain - cl_amp_brain * Cint_brain_gs704277
    d/dt(int_brain_gs441524mp) <- cl_amp_brain * Cint_brain_gs704277 + cl_nmp_brain * Cint_brain_gs441524 - (cl_mpn_brain + cl_mptn_brain) * Cint_brain_gs441524mp
    d/dt(is_brain_gs441524) <- q_brain * Carterial_gs441524 - q_brain * Cis_brain_gs441524 - clin_brain_gs441524 * fu_p_gs441524 * (Cis_brain_gs441524 - Cint_brain_gs441524 / kp_brain_gs441524)
    d/dt(int_brain_gs441524) <- clin_brain_gs441524 * fu_p_gs441524 * (Cis_brain_gs441524 - Cint_brain_gs441524 / kp_brain_gs441524) + cl_mpn_brain * Cint_brain_gs441524mp - cl_nmp_brain * Cint_brain_gs441524
    d/dt(int_brain_gs443902) <- cl_mptn_brain * Cint_brain_gs441524mp - cl_tn_brain * Cint_brain_gs443902

    # ---- GI ----
    d/dt(is_gut) <- q_gut * Carterial - q_gut * Cis_gut - clin_gut * fu_p * (Cis_gut - Cint_gut / kp_gut)
    d/dt(int_gut) <- clin_gut * fu_p * (Cis_gut - Cint_gut / kp_gut) - cl_rdva_gut * Cint_gut
    d/dt(int_gut_gs704277) <- cl_rdva_gut * Cint_gut - cl_amp_gut * Cint_gut_gs704277
    d/dt(int_gut_gs441524mp) <- cl_amp_gut * Cint_gut_gs704277 + cl_nmp_gut * Cint_gut_gs441524 - (cl_mpn_gut + cl_mptn_gut) * Cint_gut_gs441524mp
    d/dt(is_gut_gs441524) <- q_gut * Carterial_gs441524 - q_gut * Cis_gut_gs441524 - clin_gut_gs441524 * fu_p_gs441524 * (Cis_gut_gs441524 - Cint_gut_gs441524 / kp_gut_gs441524)
    d/dt(int_gut_gs441524) <- clin_gut_gs441524 * fu_p_gs441524 * (Cis_gut_gs441524 - Cint_gut_gs441524 / kp_gut_gs441524) + cl_mpn_gut * Cint_gut_gs441524mp - cl_nmp_gut * Cint_gut_gs441524
    d/dt(int_gut_gs443902) <- cl_mptn_gut * Cint_gut_gs441524mp - cl_tn_gut * Cint_gut_gs443902

    # ---- Heart ----
    d/dt(is_heart) <- q_heart * Carterial - q_heart * Cis_heart - clin_heart * fu_p * (Cis_heart - Cint_heart / kp_heart)
    d/dt(int_heart) <- clin_heart * fu_p * (Cis_heart - Cint_heart / kp_heart) - cl_rdva_heart * Cint_heart
    d/dt(int_heart_gs704277) <- cl_rdva_heart * Cint_heart - cl_amp_heart * Cint_heart_gs704277
    d/dt(int_heart_gs441524mp) <- cl_amp_heart * Cint_heart_gs704277 + cl_nmp_heart * Cint_heart_gs441524 - (cl_mpn_heart + cl_mptn_heart) * Cint_heart_gs441524mp
    d/dt(is_heart_gs441524) <- q_heart * Carterial_gs441524 - q_heart * Cis_heart_gs441524 - clin_heart_gs441524 * fu_p_gs441524 * (Cis_heart_gs441524 - Cint_heart_gs441524 / kp_heart_gs441524)
    d/dt(int_heart_gs441524) <- clin_heart_gs441524 * fu_p_gs441524 * (Cis_heart_gs441524 - Cint_heart_gs441524 / kp_heart_gs441524) + cl_mpn_heart * Cint_heart_gs441524mp - cl_nmp_heart * Cint_heart_gs441524
    d/dt(int_heart_gs443902) <- cl_mptn_heart * Cint_heart_gs441524mp - cl_tn_heart * Cint_heart_gs443902

    # ---- Kidney ----
    d/dt(is_kidney) <- q_kidney * Carterial - q_kidney * Cis_kidney - clin_kidney * fu_p * (Cis_kidney - Cint_kidney / kp_kidney)
    d/dt(int_kidney) <- clin_kidney * fu_p * (Cis_kidney - Cint_kidney / kp_kidney) - cl_rdva_kidney * Cint_kidney
    d/dt(int_kidney_gs704277) <- cl_rdva_kidney * Cint_kidney - cl_amp_kidney * Cint_kidney_gs704277
    d/dt(int_kidney_gs441524mp) <- cl_amp_kidney * Cint_kidney_gs704277 + cl_nmp_kidney * Cint_kidney_gs441524 - (cl_mpn_kidney + cl_mptn_kidney) * Cint_kidney_gs441524mp
    d/dt(is_kidney_gs441524) <- q_kidney * Carterial_gs441524 - q_kidney * Cis_kidney_gs441524 - clin_kidney_gs441524 * fu_p_gs441524 * (Cis_kidney_gs441524 - Cint_kidney_gs441524 / kp_kidney_gs441524)
    d/dt(int_kidney_gs441524) <- clin_kidney_gs441524 * fu_p_gs441524 * (Cis_kidney_gs441524 - Cint_kidney_gs441524 / kp_kidney_gs441524) + cl_mpn_kidney * Cint_kidney_gs441524mp - cl_nmp_kidney * Cint_kidney_gs441524
    d/dt(int_kidney_gs443902) <- cl_mptn_kidney * Cint_kidney_gs441524mp - cl_tn_kidney * Cint_kidney_gs443902

    # ---- Liver ----
    # Hepatic artery plus the gut and spleen outflows; the liver drains
    # (q_liver + q_gut + q_spleen).
    d/dt(is_liver) <- q_liver * Carterial + q_gut * Cis_gut + q_spleen * Cis_spleen - (q_liver + q_gut + q_spleen) * Cis_liver - clin_liver * fu_p * (Cis_liver - Cint_liver / kp_liver)
    d/dt(int_liver) <- clin_liver * fu_p * (Cis_liver - Cint_liver / kp_liver) - cl_rdva_liver * Cint_liver
    d/dt(int_liver_gs704277) <- cl_rdva_liver * Cint_liver - cl_amp_liver * Cint_liver_gs704277
    d/dt(int_liver_gs441524mp) <- cl_amp_liver * Cint_liver_gs704277 + cl_nmp_liver * Cint_liver_gs441524 - (cl_mpn_liver + cl_mptn_liver) * Cint_liver_gs441524mp
    d/dt(is_liver_gs441524) <- q_liver * Carterial_gs441524 + q_gut * Cis_gut_gs441524 + q_spleen * Cis_spleen_gs441524 - (q_liver + q_gut + q_spleen) * Cis_liver_gs441524 - clin_liver_gs441524 * fu_p_gs441524 * (Cis_liver_gs441524 - Cint_liver_gs441524 / kp_liver_gs441524)
    d/dt(int_liver_gs441524) <- clin_liver_gs441524 * fu_p_gs441524 * (Cis_liver_gs441524 - Cint_liver_gs441524 / kp_liver_gs441524) + cl_mpn_liver * Cint_liver_gs441524mp - cl_nmp_liver * Cint_liver_gs441524
    d/dt(int_liver_gs443902) <- cl_mptn_liver * Cint_liver_gs441524mp - cl_tn_liver * Cint_liver_gs443902

    # ---- Lung ----
    # Venous plasma (the forcing function) enters the lung extracellular
    # space at the cardiac output and leaves it to arterial plasma.
    d/dt(is_lung) <- qc * Cc - qc * Cis_lung - clin_lung * fu_p * (Cis_lung - Cint_lung / kp_lung)
    d/dt(int_lung) <- clin_lung * fu_p * (Cis_lung - Cint_lung / kp_lung) - cl_rdva_lung * Cint_lung
    d/dt(int_lung_gs704277) <- cl_rdva_lung * Cint_lung - cl_amp_lung * Cint_lung_gs704277
    d/dt(int_lung_gs441524mp) <- cl_amp_lung * Cint_lung_gs704277 + cl_nmp_lung * Cint_lung_gs441524 - (cl_mpn_lung + cl_mptn_lung) * Cint_lung_gs441524mp
    d/dt(is_lung_gs441524) <- qc * Cc_gs441524 - qc * Cis_lung_gs441524 - clin_lung_gs441524 * fu_p_gs441524 * (Cis_lung_gs441524 - Cint_lung_gs441524 / kp_lung_gs441524)
    d/dt(int_lung_gs441524) <- clin_lung_gs441524 * fu_p_gs441524 * (Cis_lung_gs441524 - Cint_lung_gs441524 / kp_lung_gs441524) + cl_mpn_lung * Cint_lung_gs441524mp - cl_nmp_lung * Cint_lung_gs441524
    d/dt(int_lung_gs443902) <- cl_mptn_lung * Cint_lung_gs441524mp - cl_tn_lung * Cint_lung_gs443902
    # GS-704277 in the lung extracellular space only relays venous
    # GS-704277 to arterial plasma; it has no transport into lung cells.
    d/dt(is_lung_gs704277) <- qc * Cc_gs704277 - qc * Cis_lung_gs704277

    # ---- Muscle ----
    d/dt(is_muscle) <- q_muscle * Carterial - q_muscle * Cis_muscle - clin_muscle * fu_p * (Cis_muscle - Cint_muscle / kp_muscle)
    d/dt(int_muscle) <- clin_muscle * fu_p * (Cis_muscle - Cint_muscle / kp_muscle) - cl_rdva_muscle * Cint_muscle
    d/dt(int_muscle_gs704277) <- cl_rdva_muscle * Cint_muscle - cl_amp_muscle * Cint_muscle_gs704277
    d/dt(int_muscle_gs441524mp) <- cl_amp_muscle * Cint_muscle_gs704277 + cl_nmp_muscle * Cint_muscle_gs441524 - (cl_mpn_muscle + cl_mptn_muscle) * Cint_muscle_gs441524mp
    d/dt(is_muscle_gs441524) <- q_muscle * Carterial_gs441524 - q_muscle * Cis_muscle_gs441524 - clin_muscle_gs441524 * fu_p_gs441524 * (Cis_muscle_gs441524 - Cint_muscle_gs441524 / kp_muscle_gs441524)
    d/dt(int_muscle_gs441524) <- clin_muscle_gs441524 * fu_p_gs441524 * (Cis_muscle_gs441524 - Cint_muscle_gs441524 / kp_muscle_gs441524) + cl_mpn_muscle * Cint_muscle_gs441524mp - cl_nmp_muscle * Cint_muscle_gs441524
    d/dt(int_muscle_gs443902) <- cl_mptn_muscle * Cint_muscle_gs441524mp - cl_tn_muscle * Cint_muscle_gs443902

    # ---- Rest of Body ----
    d/dt(is_remainder) <- q_remainder * Carterial - q_remainder * Cis_remainder - clin_remainder * fu_p * (Cis_remainder - Cint_remainder / kp_remainder)
    d/dt(int_remainder) <- clin_remainder * fu_p * (Cis_remainder - Cint_remainder / kp_remainder) - cl_rdva_remainder * Cint_remainder
    d/dt(int_remainder_gs704277) <- cl_rdva_remainder * Cint_remainder - cl_amp_remainder * Cint_remainder_gs704277
    d/dt(int_remainder_gs441524mp) <- cl_amp_remainder * Cint_remainder_gs704277 + cl_nmp_remainder * Cint_remainder_gs441524 - (cl_mpn_remainder + cl_mptn_remainder) * Cint_remainder_gs441524mp
    d/dt(is_remainder_gs441524) <- q_remainder * Carterial_gs441524 - q_remainder * Cis_remainder_gs441524 - clin_remainder_gs441524 * fu_p_gs441524 * (Cis_remainder_gs441524 - Cint_remainder_gs441524 / kp_remainder_gs441524)
    d/dt(int_remainder_gs441524) <- clin_remainder_gs441524 * fu_p_gs441524 * (Cis_remainder_gs441524 - Cint_remainder_gs441524 / kp_remainder_gs441524) + cl_mpn_remainder * Cint_remainder_gs441524mp - cl_nmp_remainder * Cint_remainder_gs441524
    d/dt(int_remainder_gs443902) <- cl_mptn_remainder * Cint_remainder_gs441524mp - cl_tn_remainder * Cint_remainder_gs443902

    # ---- Skin ----
    d/dt(is_skin) <- q_skin * Carterial - q_skin * Cis_skin - clin_skin * fu_p * (Cis_skin - Cint_skin / kp_skin)
    d/dt(int_skin) <- clin_skin * fu_p * (Cis_skin - Cint_skin / kp_skin) - cl_rdva_skin * Cint_skin
    d/dt(int_skin_gs704277) <- cl_rdva_skin * Cint_skin - cl_amp_skin * Cint_skin_gs704277
    d/dt(int_skin_gs441524mp) <- cl_amp_skin * Cint_skin_gs704277 + cl_nmp_skin * Cint_skin_gs441524 - (cl_mpn_skin + cl_mptn_skin) * Cint_skin_gs441524mp
    d/dt(is_skin_gs441524) <- q_skin * Carterial_gs441524 - q_skin * Cis_skin_gs441524 - clin_skin_gs441524 * fu_p_gs441524 * (Cis_skin_gs441524 - Cint_skin_gs441524 / kp_skin_gs441524)
    d/dt(int_skin_gs441524) <- clin_skin_gs441524 * fu_p_gs441524 * (Cis_skin_gs441524 - Cint_skin_gs441524 / kp_skin_gs441524) + cl_mpn_skin * Cint_skin_gs441524mp - cl_nmp_skin * Cint_skin_gs441524
    d/dt(int_skin_gs443902) <- cl_mptn_skin * Cint_skin_gs441524mp - cl_tn_skin * Cint_skin_gs443902

    # ---- Spleen ----
    d/dt(is_spleen) <- q_spleen * Carterial - q_spleen * Cis_spleen - clin_spleen * fu_p * (Cis_spleen - Cint_spleen / kp_spleen)
    d/dt(int_spleen) <- clin_spleen * fu_p * (Cis_spleen - Cint_spleen / kp_spleen) - cl_rdva_spleen * Cint_spleen
    d/dt(int_spleen_gs704277) <- cl_rdva_spleen * Cint_spleen - cl_amp_spleen * Cint_spleen_gs704277
    d/dt(int_spleen_gs441524mp) <- cl_amp_spleen * Cint_spleen_gs704277 + cl_nmp_spleen * Cint_spleen_gs441524 - (cl_mpn_spleen + cl_mptn_spleen) * Cint_spleen_gs441524mp
    d/dt(is_spleen_gs441524) <- q_spleen * Carterial_gs441524 - q_spleen * Cis_spleen_gs441524 - clin_spleen_gs441524 * fu_p_gs441524 * (Cis_spleen_gs441524 - Cint_spleen_gs441524 / kp_spleen_gs441524)
    d/dt(int_spleen_gs441524) <- clin_spleen_gs441524 * fu_p_gs441524 * (Cis_spleen_gs441524 - Cint_spleen_gs441524 / kp_spleen_gs441524) + cl_mpn_spleen * Cint_spleen_gs441524mp - cl_nmp_spleen * Cint_spleen_gs441524
    d/dt(int_spleen_gs443902) <- cl_mptn_spleen * Cint_spleen_gs441524mp - cl_tn_spleen * Cint_spleen_gs443902

    # ---------------- Micromolar outputs ----------------
    # Conversion factors are those of the deposited code (uM per mg/L).
    # They are applied to species masses that the model carries without
    # molecular-weight conversion (see description); 1.66 matches
    # remdesivir (MW 602.6) and 3.4 GS-441524 (MW 291.3), while the
    # GS-443902 factor 2.16 is the code's value.
    Cis_lung_uM <- 1.66 * Cis_lung # ffRDV_v50_multipleDose.csl: RDVluec_uM = 1.66 * RDVluec (Figure 5a)
    Cpbmc_uM <- 1.66 * pbmc # ffRDV_v50.csl: RDVpbmc_uM = 1.66 * RDVpbmc
    Cpbmc_gs704277_uM <- 2.27 * pbmc_gs704277 # ffRDV_v50.csl: Apbmc_uM = Apbmc * 2.27
    Cpbmc_gs441524mp_uM <- 2.71 * pbmc_gs441524mp # ffRDV_v50.csl: MPpbmc_uM = MPpbmc * 2.71
    Cpbmc_gs441524_uM <- 3.4 * pbmc_gs441524 # ffRDV_v50.csl: Npbmc_uM = 3.4 * Npbmc
    Cpbmc_gs443902_uM <- 2.16 * pbmc_gs443902 # ffRDV_v50.csl: TNpbmc_uM = 2.16 * TNpbmc (Figure 4)
    Cint_adipose_gs443902_uM <- 2.16 * Cint_adipose_gs443902 # code: TN<tissue>ic_uM = 2.16 * TN<tissue>ic
    Cint_bone_gs443902_uM <- 2.16 * Cint_bone_gs443902 # code: TN<tissue>ic_uM = 2.16 * TN<tissue>ic
    Cint_brain_gs443902_uM <- 2.16 * Cint_brain_gs443902 # code: TN<tissue>ic_uM = 2.16 * TN<tissue>ic
    Cint_gut_gs443902_uM <- 2.16 * Cint_gut_gs443902 # code: TN<tissue>ic_uM = 2.16 * TN<tissue>ic
    Cint_heart_gs443902_uM <- 2.16 * Cint_heart_gs443902 # code: TN<tissue>ic_uM = 2.16 * TN<tissue>ic
    Cint_kidney_gs443902_uM <- 2.16 * Cint_kidney_gs443902 # code: TN<tissue>ic_uM = 2.16 * TN<tissue>ic
    Cint_liver_gs443902_uM <- 2.16 * Cint_liver_gs443902 # code: TN<tissue>ic_uM = 2.16 * TN<tissue>ic
    Cint_lung_gs443902_uM <- 2.16 * Cint_lung_gs443902 # code: TN<tissue>ic_uM = 2.16 * TN<tissue>ic
    Cint_muscle_gs443902_uM <- 2.16 * Cint_muscle_gs443902 # code: TN<tissue>ic_uM = 2.16 * TN<tissue>ic
    Cint_remainder_gs443902_uM <- 2.16 * Cint_remainder_gs443902 # code: TN<tissue>ic_uM = 2.16 * TN<tissue>ic
    Cint_skin_gs443902_uM <- 2.16 * Cint_skin_gs443902 # code: TN<tissue>ic_uM = 2.16 * TN<tissue>ic
    Cint_spleen_gs443902_uM <- 2.16 * Cint_spleen_gs443902 # code: TN<tissue>ic_uM = 2.16 * TN<tissue>ic

    Cc ~ add(addSd)
    Cc_gs704277 ~ add(addSd_gs704277)
    Cc_gs441524 ~ add(addSd_gs441524)
  })
}
