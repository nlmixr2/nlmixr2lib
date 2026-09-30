Stader_2021_bictegravir_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, Stader et al. Matlab 2017a framework, as deposited",
    "with the paper). Oral bictegravir (single-agent tablet) in adults aged",
    "20 to 99 years, used to predict the continuous effect of ageing on",
    "bictegravir pharmacokinetics in people living with HIV. Sixteen organs",
    "(lung, adipose, bone, brain, gonads, heart, kidney, muscle, skin,",
    "thymus, gut, spleen, pancreas, liver, lymph node and a remaining",
    "tissue) with vascular (_vas), interstitial (_ew) and intracellular",
    "(_iw) sub-compartments (the gut carries vascular and interstitial",
    "only), venous and arterial blood, and a compartmental absorption and",
    "transit (CAT) gut of stomach plus duodenum, jejunum, ileum and colon,",
    "each with luminal fluid, an uptake layer and enterocytes: 63 ODEs.",
    "Organ weights, blood flows, lymph flows, haematocrit, albumin, GFR,",
    "microsomal protein per gram liver, hepatic CYP3A4 / UGT1A1 and",
    "intestinal CYP3A4 abundances and GI transit times are age-, sex-,",
    "height- and weight-dependent regressions from the deposited",
    "virtual-population generator; the generator's per-subject random",
    "draws are carried as fixed etas so that rxSolve() with etas sampled",
    "reproduces the virtual population. Distribution by Rodgers and",
    "Rowland partitioning with permeability-limited cellular uptake;",
    "hepatic elimination by CYP3A4, UGT1A1 and an unassigned pathway,",
    "intestinal CYP3A4 metabolism in the enterocytes, and GFR-scaled renal",
    "clearance. Cc is the framework's reported plasma concentration, which",
    "the deposited code reads from the venous-blood state.",
    sep = " "
  )
  reference <- paste(
    "Stader F, Courlet P, Decosterd LA, Battegay M, Marzolini C.",
    "Physiologically-Based Pharmacokinetic Modeling Combined with Swiss HIV",
    "Cohort Study Data Supports No Dose Adjustment of Bictegravir in Elderly",
    "Individuals Living With HIV. Clin Pharmacol Ther. 2021;109(4):1025-1029.",
    "doi:10.1002/cpt.2178. Model code: Supplementary Material s002",
    "(CPT_Matlab_Code, Matlab source of the PBPK framework including",
    "Drug/DrugLibrary/bictegravir.m); drug parameters: Supplementary Table S1.",
    sep = " "
  )
  vignette <- "Stader_2021_bictegravir_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the age regressions of the deposited virtual-population",
        "generator: height and weight means (supplied here as covariates),",
        "lung, adipose, brain, gonad, heart, kidney, muscle and skin weights,",
        "blood weight, albumin, brain vascular fraction, cardiac output,",
        "adipose, kidney and liver blood-flow fractions, GFR and microsomal",
        "protein per gram liver. The generator draws integer ages (rounded",
        "Weibull draws, 61.73 scale and 1.55 shape) restricted to 20-99",
        "years; predictions outside 20-99 years are not supported by the",
        "framework.",
        sep = " "
      ),
      source_name = "BODY.Age"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "male (0)",
      notes = paste(
        "The deposited code stores sex as BODY.Sex with 1 = female and",
        "0 = male (PBPK_Population_Demographics.m, Generate_Sex: 'assign a",
        "1 to females and 0 to males'), so it maps onto SEXF without",
        "transformation.",
        sep = " "
      ),
      source_name = "BODY.Sex"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Individual height. The generator draws it as",
        "normal(mean = -0.0039*AGE^2 + 0.238*AGE - 12.5*SEXF + 176,",
        "CV 3.8%) and resets it to the mean when the resulting BMI falls",
        "outside 18.5-30 kg/m^2 (PBPK_Population_Demographics.m); the",
        "validation vignette reproduces that draw. Enters BSA, lung, bone,",
        "brain and gut weights, adipose weight, heart and kidney blood",
        "flows and intestinal segment lengths.",
        sep = " "
      ),
      source_name = "BODY.Height"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Individual weight. The generator draws it as",
        "normal(mean = -0.0039*AGE^2 + 1.12*HT + 0.611*AGE - 0.424*SEXF - 137,",
        "CV 15.2%) with the same BMI reset as height. Enters BSA, adipose,",
        "gonad and lymph-node weights, total lymph flow, and the remaining",
        "tissue as the body-weight balance.",
        sep = " "
      ),
      source_name = "BODY.Weight"
    )
  )

  # Organ sub-compartment states compose registered organ roots with the
  # registered `vas` / `ew` / `iw` suffixes. Thymus and gonads have no
  # registered bare organ root, and the CAT uptake layer and per-segment
  # enterocytes have no registered canonical, so those states are declared
  # here as paper-specific.
  paper_specific_compartments <- c(
    "thymus_vas",
    "thymus_ew",
    "thymus_iw",
    "gonads_vas",
    "gonads_ew",
    "gonads_iw",
    "duodenum_uptake",
    "jejunum_uptake",
    "ileum_uptake",
    "colon_uptake",
    "duodenum_enterocyte",
    "jejunum_enterocyte",
    "ileum_enterocyte",
    "colon_enterocyte"
  )

  # The etas are the per-subject random draws of the deposited
  # virtual-population generator, not estimated between-subject variances.
  paper_specific_etas <- c(
    "etahct",
    "etahsa",
    "etaw_adipose",
    "etaw_bone",
    "etaw_brain",
    "etaw_gonads",
    "etaw_heart",
    "etaw_kidney",
    "etaw_muscle",
    "etaw_skin",
    "etaw_thymus",
    "etaw_gut",
    "etaw_spleen",
    "etaw_pancreas",
    "etaw_liver",
    "etaw_lnode",
    "etaw_blood",
    "etaco",
    "etaltot",
    "etagfr",
    "etamppgl",
    "etamppgl_redraw",
    "etacyp3a4_liver",
    "etacyp3a4_liver_redraw",
    "etaugt1a1_liver",
    "etaugt1a1_liver_redraw",
    "etagastric",
    "etasitt",
    "etasitt_redraw",
    "etacolt",
    "etacolt_redraw",
    "etacyp3a4_gut",
    "etacyp3a4_gut_redraw"
  )

  compartmentData <- list(
    lung_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    lung_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    adipose_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    adipose_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    adipose_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    bone_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    bone_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    bone_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    brain_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    brain_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    brain_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    gonads_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    gonads_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    gonads_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    heart_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    heart_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    heart_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    kidney_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    muscle_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    muscle_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    muscle_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    skin_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    skin_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    skin_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    thymus_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    thymus_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    thymus_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    gut_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    gut_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    spleen_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    spleen_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    spleen_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    pancreas_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    pancreas_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    pancreas_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    liver_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    liver_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    liver_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    lnode_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    lnode_ew = list(analyte = "bictegravir", units = "mg", specimen = "lymph", verified = TRUE),
    lnode_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    other_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    other_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    other_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    venous = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    stomach = list(analyte = "bictegravir", units = "mg", specimen = "administration site", verified = TRUE),
    duodenum = list(analyte = "bictegravir", units = "mg", specimen = "administration site", verified = TRUE),
    jejunum = list(analyte = "bictegravir", units = "mg", specimen = "administration site", verified = TRUE),
    ileum = list(analyte = "bictegravir", units = "mg", specimen = "administration site", verified = TRUE),
    colon = list(analyte = "bictegravir", units = "mg", specimen = "administration site", verified = TRUE),
    duodenum_uptake = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    jejunum_uptake = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    ileum_uptake = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    colon_uptake = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    duodenum_enterocyte = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    jejunum_enterocyte = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    ileum_enterocyte = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    colon_enterocyte = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    a_feces = list(analyte = "bictegravir", units = "mg", specimen = "faeces", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 92L,
    n_studies = 2L,
    age_range = "20-99 years (virtual); 22.8-81.1 years (verification cohort)",
    sex_female_pct = 50,
    disease_state = "Healthy volunteers (model development) and people living with HIV (verification)",
    dose_range = "5-600 mg oral bictegravir, single dose and once daily",
    regions = "Switzerland (Swiss HIV Cohort Study verification data)",
    notes = paste(
      "A PBPK model, not a fit to individual data. The drug model was",
      "developed against published phase I data in healthy volunteers",
      "given 5-600 mg single and once-daily doses (Gallant 2017; Table S2)",
      "and verified against therapeutic drug monitoring samples from 60",
      "young (mean 42.2 years, range 22.8-54.7) and 32 elderly (mean 63.8",
      "years, range 55.0-81.1) people living with HIV in the Swiss HIV",
      "Cohort Study, all on 50 mg bictegravir with no CYP3A or UGT1A1",
      "inhibitor or inducer (Results; Table 1). Simulations use 10 trials",
      "of 10 virtual individuals (50% women) for the phase I and TDM",
      "comparisons, and 500 virtual individuals (50% women) per 5-year age",
      "band from 20 to 99 years for the ageing analysis (Methods).",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Drug parameters. Values are the deposited bictegravir drug file
    # (s002 CPT_Matlab_Code/Drug/DrugLibrary/bictegravir.m), which Table S1
    # reprints (rounded) with the same sources. All are fixed inputs of a
    # PBPK model, not estimates. The molecular weight (449.39 g/mol,
    # bictegravir.m DRUG.MolW; Table S1 'MW' 449.4) cancels out of this
    # linear system once the states are carried in mg, so it is not a
    # parameter here.
    logp <- fixed(1.28)
    label("Octanol:water partition coefficient logP (unitless)") # bictegravir.m DRUG.logP = 1.28; Table S1 'logP' 1.28
    pka <- fixed(9.81)
    label("pKa of the monoprotic acid (unitless)") # bictegravir.m DRUG.pka1 = 9.81 (DRUG.type = mono_acid); Table S1 'pKa 1' 9.8, 'drug type' ma
    bp <- fixed(0.64)
    label("Blood:plasma concentration ratio (unitless)") # bictegravir.m DRUG.BP = 0.64; Table S1 'BP' 0.64
    fu <- fixed(0.0025)
    label("Fraction unbound in plasma at the reference albumin concentration of 45.6 g/L (unitless)") # bictegravir.m DRUG.fu = 0.0025 (binding protein albumin); Table S1 'fup' 0.0025
    papp <- fixed(24.6)
    label("Caco-2 apparent permeability (1e-6 cm/s)") # bictegravir.m DRUG.Papp = 24.6; Table S1 'Papp' 24.6
    jin_all <- fixed(0.7)
    label("Cellular influx:efflux scalar, all tissues (unitless)") # bictegravir.m DRUG.JinScalarAll = 0.7; Table S1 'Tissue Scalar All' 0.7 (optimized)
    jin_liver <- fixed(2.0)
    label("Cellular influx:efflux scalar, liver-specific multiplier (unitless)") # bictegravir.m DRUG.JinScalar(liver) = 2.0; Table S1 'Tissue Scalar LI' 2.0 (optimized)
    clint_cyp3a4 <- fixed(0.114)
    label("CYP3A4 intrinsic clearance (uL/min/pmol CYP3A4)") # bictegravir.m DRUG.CLint_CYP_1(CYP3A4) = 0.114; Table S1 'CYP3A4 CLint' 0.114 (retrograde)
    clint_ugt1a1 <- fixed(0.292)
    label("UGT1A1 intrinsic clearance (uL/min/pmol UGT1A1)") # bictegravir.m DRUG.CLint_UGT_1(UGT1A1) = 0.292; Table S1 'UGT1A1 CLint' 0.292 (retrograde)
    clint_hep <- fixed(3.993)
    label("Hepatic intrinsic clearance not assigned to an enzyme (uL/min/mg microsomal protein)") # bictegravir.m DRUG.CLint = 3.993; Table S1 'Unspecified' 3.993 (retrograde)
    lcl_renal <- fixed(log(0.0043))
    label("Log renal clearance at GFR 130 mL/min (men) or 120 mL/min (women) and reference albumin (log L/h)") # bictegravir.m DRUG.CLrenal = 0.0043 L/h; Table S1 'CLrenal' 0.004

    # ---------------------------------------------------------------------
    # Virtual-population random draws. The deposited generator draws each
    # quantity as normrnd(Mean, (CV/100)*Mean), so the value is
    # Mean*(1 + eta) with eta ~ N(0, (CV/100)^2); the variances below are
    # (CV/100)^2 of the printed CVs. Source file noted per line
    # (s002 CPT_Matlab_Code/Population/...).
    etahct ~ fixed(0.020736) # HCT_CV = 14.4 (PBPK_Population_Tissue.m)
    etahsa ~ fixed(0.006241) # HSA_CV = 7.9 (PBPK_Population_Tissue.m)
    etaw_adipose ~ fixed(0.087616) # WAD_CV = 29.6 (PBPK_Population_Tissue.m)
    etaw_bone ~ fixed(0.017424) # WBO_CV = 13.2 (PBPK_Population_Tissue.m)
    etaw_brain ~ fixed(0.0081) # WBR_CV = 9.0 (PBPK_Population_Tissue.m)
    etaw_gonads ~ fixed(0.121104) # WGO_CV = 34.8 (PBPK_Population_Tissue.m)
    etaw_heart ~ fixed(0.041209) # WHE_CV = 20.3 (PBPK_Population_Tissue.m)
    etaw_kidney ~ fixed(0.045369) # WKI_CV = 21.3 (PBPK_Population_Tissue.m)
    etaw_muscle ~ fixed(0.013924) # WMU_CV = 11.8 (PBPK_Population_Tissue.m)
    etaw_skin ~ fixed(0.006889) # WSK_CV = 8.3 (PBPK_Population_Tissue.m)
    etaw_thymus ~ fixed(0.200704) # WTH_CV = 44.8 (PBPK_Population_Tissue.m)
    etaw_gut ~ fixed(0.005329) # WGU_CV = 7.3 (PBPK_Population_Tissue.m)
    etaw_spleen ~ fixed(0.267289) # WSP_CV = 51.7 (PBPK_Population_Tissue.m)
    etaw_pancreas ~ fixed(0.077284) # WPA_CV = 27.8 (PBPK_Population_Tissue.m)
    etaw_liver ~ fixed(0.056169) # WLI_CV = 23.7 (PBPK_Population_Tissue.m)
    etaw_lnode ~ fixed(0.04) # WLN_CV = 20 (PBPK_Population_Tissue.m)
    etaw_blood ~ fixed(0.010816) # WBL_CV = 10.4 (PBPK_Population_Tissue.m)
    etaco ~ fixed(0.044521) # CO_CV = 21.1 (PBPK_Population_Tissue.m)
    etaltot ~ fixed(0.09) # LT_CV = 30 (PBPK_Population_Tissue.m)
    etagfr ~ fixed(0.021609) # GFR_CV = 14.7 (PBPK_Population_Tissue.m)
    etamppgl ~ fixed(0.2116) # MPPGL_CV = 46 (PBPK_Population_Liver.m; 46 for every drug but ritonavir)
    etamppgl_redraw ~ fixed(0.002116) # Gen_LIPar redraw normrnd(Mean, (CV/1000)*Mean) = 4.6% (PBPK_Population_Liver.m)
    etacyp3a4_liver ~ fixed(0.6561) # hepatic CYP3A4_CV = 81 (PBPK_Population_Liver.m)
    etacyp3a4_liver_redraw ~ fixed(0.006561) # Gen_LIPar redraw at CV/10 = 8.1% (PBPK_Population_Liver.m)
    etaugt1a1_liver ~ fixed(0.36) # UGT1A1_CV = 60 (PBPK_Population_Liver.m)
    etaugt1a1_liver_redraw ~ fixed(0.0036) # Gen_LIPar redraw at CV/10 = 6.0% (PBPK_Population_Liver.m)
    etagastric ~ fixed(1) # GET = 0.25 + (1.00 - 0.25)*rand; a standard normal mapped to U(0,1) by phi() (PBPK_Population_GIT.m)
    etasitt ~ fixed(0.161604) # SIT_CV = 40.2 (PBPK_Population_GIT.m)
    etasitt_redraw ~ fixed(0.00161604) # Gen_GUPar redraw at CV/10 = 4.02% (PBPK_Population_GIT.m)
    etacolt ~ fixed(0.04) # CNT_CV = 20.0 (PBPK_Population_GIT.m)
    etacolt_redraw ~ fixed(0.0004) # Gen_GUPar redraw at CV/10 = 2.0% (PBPK_Population_GIT.m)
    etacyp3a4_gut ~ fixed(0.36) # intestinal CYP3A4_CV = 60 (PBPK_Population_GIT.m)
    etacyp3a4_gut_redraw ~ fixed(0.0036) # Gen_GUPar redraw at CV/10 = 6.0% (PBPK_Population_GIT.m)
  })

  model({
    # Paths below are files of the deposited s002 CPT_Matlab_Code folder.
    # ====================================================================
    # 1. Demographics (Population/PBPK_Population_Demographics.m)
    # ====================================================================
    bsa <- 0.007184 * HT^0.725 * WT^0.425 # BSA, DuBois & DuBois (m^2)

    # ====================================================================
    # 2. Blood (Population/PBPK_Population_Tissue.m, 'Blood parameters').
    # GenTisPar truncates at Min/Max only when length(Min) is 1 or 2,
    # i.e. for populations of one or two subjects; for the 100- and
    # 500-subject populations of the paper the draws are untruncated.
    # ====================================================================
    hct <- (0.443 - 0.033 * SEXF) * (1 + etahct) # HCT_Mean
    hsa <- (-0.0709 * AGE + 47.7) * (1 + etahsa) # HSA_Mean (g/L)

    # ====================================================================
    # 3. Organ weights in kg (PBPK_Population_Tissue.m, 'Organ weights').
    # A negative draw is replaced by the block's alternative (WXX_Alt);
    # a negative mean makes normrnd return NaN, which GenTisPar replaces
    # by the (negative) mean and then by the alternative.
    # ====================================================================
    w_lung <- exp(0.00771 * AGE + 0.0279 * HT - 5.58) # WLU_Mean; WLU_CV = 0

    w_adipose_tv <- 0.68 * WT - 0.56 * HT + 6.1 * SEXF + 65 # WAD_Mean
    w_adipose_alt <- ((0.1875 - 0.0796 * SEXF) * AGE + (14.296 + 17.305 * SEXF)) * (WT / 100) # WAD_Alt
    w_adipose_0 <- w_adipose_tv * (1 + etaw_adipose)
    if (w_adipose_tv < 0 || w_adipose_0 < 0) {
      w_adipose_0 <- w_adipose_alt
    }

    w_bone_tv <- exp(0.024 * HT - 1.9) # WBO_Mean
    w_bone <- w_bone_tv * (1 + etaw_bone)
    if (w_bone < 0) {
      w_bone <- w_bone_tv
    }

    w_brain_tv <- exp(-0.00075 * AGE + 0.00778 * HT - 0.97) # WBR_Mean
    w_brain <- w_brain_tv * (1 + etaw_brain)
    if (w_brain < 0) {
      w_brain <- w_brain_tv
    }

    w_gonads_tv <- -0.00022 * AGE - 0.00034 * WT - 0.030 * SEXF + 0.072 # WGO_Mean
    w_gonads_alt <- (0.00049 - 0.00037 * SEXF) * WT # WGO_Alt
    w_gonads <- w_gonads_tv * (1 + etaw_gonads)
    if (w_gonads_tv < 0 || w_gonads < 0) {
      w_gonads <- w_gonads_alt
    }

    w_heart_tv <- 0.34 * bsa + 0.0018 * AGE - 0.36 # WHE_Mean
    w_heart <- w_heart_tv * (1 + etaw_heart)
    if (w_heart < 0) {
      w_heart <- w_heart_tv
    }

    w_kidney_tv <- -0.00038 * AGE - 0.056 * SEXF + 0.33 # WKI_Mean
    w_kidney <- w_kidney_tv * (1 + etaw_kidney)
    if (w_kidney < 0) {
      w_kidney <- w_kidney_tv
    }

    w_muscle_tv <- 17.9 * bsa - 0.0667 * AGE - 5.68 * SEXF - 1.22 # WMU_Mean
    w_muscle_0 <- w_muscle_tv * (1 + etaw_muscle)
    if (w_muscle_0 < 0) {
      w_muscle_0 <- w_muscle_tv
    }

    w_skin_tv <- exp(-0.0058 * AGE - 0.37 * SEXF + 1.13) # WSK_Mean
    w_skin <- w_skin_tv * (1 + etaw_skin)
    if (w_skin < 0) {
      w_skin <- w_skin_tv
    }

    w_thymus_tv <- 0.0221 # WTH_Mean
    w_thymus <- w_thymus_tv * (1 + etaw_thymus)
    if (w_thymus < 0) {
      w_thymus <- w_thymus_tv
    }

    w_gut_tv <- 3e-06 * HT^2.49 # WGU_Mean
    w_gut <- w_gut_tv * (1 + etaw_gut)
    if (w_gut < 0) {
      w_gut <- w_gut_tv
    }

    w_spleen_tv <- exp(1.13 * bsa - 3.93) # WSP_Mean
    w_spleen <- w_spleen_tv * (1 + etaw_spleen)
    if (w_spleen < 0) {
      w_spleen <- 0.00228 * WT # WSP_Alt
    }

    w_pancreas_tv <- 0.103 # WPA_Mean
    w_pancreas <- w_pancreas_tv * (1 + etaw_pancreas)
    if (w_pancreas < 0) {
      w_pancreas <- w_pancreas_tv
    }

    w_liver_tv <- 0.52 * bsa^1.73 # WLI_Mean
    w_liver <- w_liver_tv * (1 + etaw_liver)
    if (w_liver < 0) {
      w_liver <- w_liver_tv
    }

    w_lnode_tv <- 0.00386 * WT # WLN_Mean (Gill et al. 2016)
    w_lnode <- w_lnode_tv * (1 + etaw_lnode)
    if (w_lnode < 0) {
      w_lnode <- w_lnode_tv
    }

    w_blood <- (-0.0098 * AGE - 1.89 * SEXF + 6.06) * (1 + etaw_blood) # WBL_Mean
    w_plasma <- w_blood * (1 - 0.91 * hct)
    w_rbc <- w_blood - w_plasma

    # Remaining tissue: body weight balance, held at >= 2% of body weight
    # by taking the shortfall 75:25 from adipose and muscle, and adipose
    # held at >= 2% of body weight by taking the shortfall from muscle.
    w_other_0 <- WT - w_blood - (w_lung + w_adipose_0 + w_bone + w_brain + w_gonads +
      w_heart + w_kidney + w_muscle_0 + w_skin + w_thymus + w_gut + w_spleen +
      w_pancreas + w_liver + w_lnode)
    mirema <- 0
    if (w_other_0 < 0.02 * WT) {
      mirema <- 0.02 * WT - w_other_0
    }
    w_adipose_1 <- w_adipose_0 - 0.75 * mirema
    w_muscle_1 <- w_muscle_0 - 0.25 * mirema
    miadma <- 0
    if (w_adipose_1 < 0.02 * WT) {
      miadma <- 0.02 * WT - w_adipose_1
    }
    w_adipose <- w_adipose_1 + miadma
    w_muscle <- w_muscle_1 - miadma
    w_other <- WT - w_blood - (w_lung + w_adipose + w_bone + w_brain + w_gonads +
      w_heart + w_kidney + w_muscle + w_skin + w_thymus + w_gut + w_spleen +
      w_pancreas + w_liver + w_lnode)
    w_all <- w_lung + w_adipose + w_bone + w_brain + w_gonads + w_heart + w_kidney +
      w_muscle + w_skin + w_thymus + w_gut + w_spleen + w_pancreas + w_liver +
      w_lnode + w_other + w_plasma + w_rbc

    # ====================================================================
    # 4. Organ volumes in L (PBPK_Population_Tissue.m, 'Organ density').
    # The remaining-tissue entries of the density and composition tables
    # are initialised to 1 and then replaced by the Worg-weighted mean
    # over all 18 compartments, a mean that includes the remaining tissue
    # itself at that initial value of 1 (as-run).
    # ====================================================================
    den_other <- (w_lung * 1.0 + w_adipose * 0.916 + w_bone * 1.9 + w_brain * 1.04 +
      w_gonads * 1.045 + w_heart * 1.03 + w_kidney * 1.05 + w_muscle * 1.041 +
      w_skin * 1.1 + w_thymus * 1.025 + w_gut * 1.042 + w_spleen * 1.06 +
      w_pancreas * 1.045 + w_liver * 1.08 + w_lnode * 1.0 + w_other * 1 +
      w_plasma * 1.027 + w_rbc * 1.09) / w_all
    v_lung <- w_lung / 1.0
    v_adipose <- w_adipose / 0.916
    v_bone <- w_bone / 1.9
    v_brain <- w_brain / 1.04
    v_gonads <- w_gonads / 1.045
    v_heart <- w_heart / 1.03
    v_kidney <- w_kidney / 1.05
    v_muscle <- w_muscle / 1.041
    v_skin <- w_skin / 1.1
    v_thymus <- w_thymus / 1.025
    v_gut <- w_gut / 1.042
    v_spleen <- w_spleen / 1.06
    v_pancreas <- w_pancreas / 1.045
    v_liver <- w_liver / 1.08
    v_lnode <- w_lnode / 1.0
    v_other <- w_other / den_other
    v_venous <- (2 / 3) * w_blood / 1.06 # Vvein = Wvein/1.06
    v_arterial <- (1 / 3) * w_blood / 1.06 # Vartery = Wartery/1.06

    # ====================================================================
    # 5. Tissue composition of the remaining tissue (PBPK_Population_Tissue.m,
    # 'Tissue composition'): Worg-weighted means of the tabulated FraEW,
    # FraIW, FraNL, FraNP and KpHSA, and of FraVas.
    # ====================================================================
    fe_other <- (w_lung * 0.348 + w_adipose * 0.141 + w_bone * 0.098 + w_brain * 0.092 +
      w_gonads * 0.239 + w_heart * 0.313 + w_kidney * 0.283 + w_muscle * 0.091 +
      w_skin * 0.623 + w_thymus * 0.150 + w_gut * 0.267 + w_spleen * 0.208 +
      w_pancreas * 0.12 + w_liver * 0.165 + w_lnode * 0.208 + w_other * 1 +
      w_plasma * 0.945 + w_rbc * 0) / w_all
    fi_other <- (w_lung * 0.463 + w_adipose * 0.039 + w_bone * 0.341 + w_brain * 0.678 +
      w_gonads * 0.561 + w_heart * 0.445 + w_kidney * 0.5 + w_muscle * 0.669 +
      w_skin * 0.0947 + w_thymus * 0.626 + w_gut * 0.451 + w_spleen * 0.58 +
      w_pancreas * 0.664 + w_liver * 0.586 + w_lnode * 0.58 + w_other * 1 +
      w_plasma * 0 + w_rbc * 0.666) / w_all
    fnl_other <- (w_lung * 0.003 + w_adipose * 0.79 + w_bone * 0.074 + w_brain * 0.051 +
      w_gonads * 0.007 + w_heart * 0.015 + w_kidney * 0.0207 + w_muscle * 0.0238 +
      w_skin * 0.0248 + w_thymus * 0.017 + w_gut * 0.0487 + w_spleen * 0.0201 +
      w_pancreas * 0.041 + w_liver * 0.0348 + w_lnode * 0.0201 + w_other * 1 +
      w_plasma * 0.35 + w_rbc * 0.17) / w_all
    fnp_other <- (w_lung * 0.009 + w_adipose * 0.002 + w_bone * 0.0011 + w_brain * 0.0565 +
      w_gonads * 0.0077 + w_heart * 0.0166 + w_kidney * 0.0162 + w_muscle * 0.0072 +
      w_skin * 0.0111 + w_thymus * 0.0092 + w_gut * 0.0163 + w_spleen * 0.0198 +
      w_pancreas * 0.0093 + w_liver * 0.0252 + w_lnode * 0.0198 + w_other * 1 +
      w_plasma * 0.225 + w_rbc * 0.29) / w_all
    kphsa_other <- (w_lung * 0.212 + w_adipose * 0.021 + w_bone * 0.1 + w_brain * 0.048 +
      w_gonads * 0.048 + w_heart * 0.157 + w_kidney * 0.13 + w_muscle * 0.025 +
      w_skin * 0.277 + w_thymus * 0.075 + w_gut * 0.158 + w_spleen * 0.097 +
      w_pancreas * 0.06 + w_liver * 0.086 + w_lnode * 0.097 + w_other * 1 +
      w_plasma * 1.0 + w_rbc * 0) / w_all
    fvas_brain <- -0.000545 * AGE + 0.056 # FraVas(brain)
    fvas_other <- (w_lung * 0.185 + w_adipose * 0.031 + w_bone * 0.05 + w_brain * fvas_brain +
      w_gonads * 0.069 + w_heart * 0.042 + w_kidney * 0.07 + w_muscle * 0.027 +
      w_skin * 0.05 + w_thymus * 0.05 + w_gut * 0.05 + w_spleen * 0.05 +
      w_pancreas * 0.05 + w_liver * 0.05 + w_lnode * 0.05 + w_other * 1 +
      w_plasma * 1 + w_rbc * 1) / w_all

    # ====================================================================
    # 6. Sub-compartment volumes (PBPK_Population_Tissue.m, 'Subcompartment
    # volume'): Vvas = FraVas*Vorg, Vint = Vorg*FraEW - Vvas*(1 - HCT),
    # Vcel = Vorg - Vvas - Vint. FraVas and FraEW per organ as tabulated.
    # ====================================================================
    vv_lung <- 0.185 * v_lung
    ve_lung <- v_lung * 0.348 - vv_lung * (1 - hct)
    vi_lung <- v_lung - vv_lung - ve_lung
    vv_adipose <- 0.031 * v_adipose
    ve_adipose <- v_adipose * 0.141 - vv_adipose * (1 - hct)
    vi_adipose <- v_adipose - vv_adipose - ve_adipose
    vv_bone <- 0.05 * v_bone
    ve_bone <- v_bone * 0.098 - vv_bone * (1 - hct)
    vi_bone <- v_bone - vv_bone - ve_bone
    vv_brain <- fvas_brain * v_brain
    ve_brain <- v_brain * 0.092 - vv_brain * (1 - hct)
    vi_brain <- v_brain - vv_brain - ve_brain
    vv_gonads <- 0.069 * v_gonads
    ve_gonads <- v_gonads * 0.239 - vv_gonads * (1 - hct)
    vi_gonads <- v_gonads - vv_gonads - ve_gonads
    vv_heart <- 0.042 * v_heart
    ve_heart <- v_heart * 0.313 - vv_heart * (1 - hct)
    vi_heart <- v_heart - vv_heart - ve_heart
    vv_kidney <- 0.07 * v_kidney
    ve_kidney <- v_kidney * 0.283 - vv_kidney * (1 - hct)
    vi_kidney <- v_kidney - vv_kidney - ve_kidney
    vv_muscle <- 0.027 * v_muscle
    ve_muscle <- v_muscle * 0.091 - vv_muscle * (1 - hct)
    vi_muscle <- v_muscle - vv_muscle - ve_muscle
    vv_skin <- 0.05 * v_skin
    ve_skin <- v_skin * 0.623 - vv_skin * (1 - hct)
    vi_skin <- v_skin - vv_skin - ve_skin
    vv_thymus <- 0.05 * v_thymus
    ve_thymus <- v_thymus * 0.150 - vv_thymus * (1 - hct)
    vi_thymus <- v_thymus - vv_thymus - ve_thymus
    # Gut: Vint uses the tissue-level Vvas (0.05*Vorg); PBPK_Population_GIT.m
    # then overwrites the gut Vvas with the sum of the four segment vascular
    # spaces, 0.05 of the duodenum, jejunum, ileum and colon shares.
    ve_gut <- v_gut * 0.267 - (0.05 * v_gut) * (1 - hct)
    vv_gut <- 0.05 * (0.060 + 0.279 + 0.316 + 0.192) * v_gut
    vv_spleen <- 0.05 * v_spleen
    ve_spleen <- v_spleen * 0.208 - vv_spleen * (1 - hct)
    vi_spleen <- v_spleen - vv_spleen - ve_spleen
    vv_pancreas <- 0.05 * v_pancreas
    ve_pancreas <- v_pancreas * 0.12 - vv_pancreas * (1 - hct)
    vi_pancreas <- v_pancreas - vv_pancreas - ve_pancreas
    vv_liver <- 0.05 * v_liver
    ve_liver <- v_liver * 0.165 - vv_liver * (1 - hct)
    vi_liver <- v_liver - vv_liver - ve_liver
    vv_lnode <- 0.05 * v_lnode
    ve_lnode <- v_lnode * 0.208 - vv_lnode * (1 - hct)
    vi_lnode <- v_lnode - vv_lnode - ve_lnode
    vv_other <- fvas_other * v_other
    ve_other <- v_other * fe_other - vv_other * (1 - hct)
    vi_other <- v_other - vv_other - ve_other

    # ====================================================================
    # 7. Blood flows (PBPK_Population_Tissue.m, 'blood flows'): percent of
    # cardiac output, then L/h. The remaining-tissue flow is the balance of
    # the listed organs; the brain flow is not subtracted from it (as-run).
    # ====================================================================
    co <- (159 * bsa - 1.56 * AGE + 114) * (1 + etaco) # CO_Mean (L/h)
    fq_adipose_0 <- (0.044 + 0.027 * SEXF) * AGE + 2.4 * SEXF + 3.9
    fq_bone <- 5
    fq_brain <- exp(-0.48 * bsa + 0.04 * SEXF + 3.5)
    fq_gonads <- -0.03 * SEXF + 0.05
    fq_heart_0 <- -0.72 * HT - 10 * SEXF + 134
    if (fq_heart_0 < 0 || fq_heart_0 > 12.0) {
      fq_heart_0 <- 0.04 + 0.01 * SEXF # 'values from Valentin (2002)', entered in percent units (as-run)
    }
    fq_kidney_0 <- -8.7 * bsa + 0.29 * HT - 0.081 * AGE - 13
    fq_muscle_0 <- -6.4 * SEXF + 17.5
    fq_skin <- 5
    fq_thymus <- 1.5
    fq_liver <- -0.108 * AGE + 1.04 * SEXF + 27.9
    fq_ha <- 6.5 # hepatic artery
    fq_pv <- fq_liver - fq_ha # portal vein
    fq_gut <- ((2 * SEXF + 14) * fq_pv) / (1.5 * SEXF + 19)
    fq_spleen <- (3 * fq_pv) / (1.5 * SEXF + 19)
    fq_pancreas <- (1 * fq_pv) / (1.5 * SEXF + 19)
    fq_by <- fq_pv - fq_gut - fq_spleen - fq_pancreas # arterial flow bypassing the portal organs
    fq_lnode <- 1.65
    fq_other_0 <- 100 - fq_adipose_0 - fq_bone - fq_gonads - fq_heart_0 - fq_kidney_0 -
      fq_muscle_0 - fq_skin - fq_thymus - fq_liver - fq_lnode
    fqmis <- 0
    if (fq_other_0 < 1) {
      fqmis <- 1 - fq_other_0
    }
    fq_adipose <- fq_adipose_0 - 0.25 * fqmis
    fq_heart <- fq_heart_0 - 0.25 * fqmis
    fq_kidney <- fq_kidney_0 - 0.25 * fqmis
    fq_muscle <- fq_muscle_0 - 0.25 * fqmis
    fq_other <- 100 - fq_adipose - fq_bone - fq_gonads - fq_heart - fq_kidney -
      fq_muscle - fq_skin - fq_thymus - fq_liver - fq_lnode

    q_lung <- co
    q_adipose <- fq_adipose / 100 * co
    q_bone <- fq_bone / 100 * co
    q_brain <- fq_brain / 100 * co
    q_gonads <- fq_gonads / 100 * co
    q_heart <- fq_heart / 100 * co
    q_kidney <- fq_kidney / 100 * co
    q_muscle <- fq_muscle / 100 * co
    q_skin <- fq_skin / 100 * co
    q_thymus <- fq_thymus / 100 * co
    q_gut <- fq_gut / 100 * co
    q_spleen <- fq_spleen / 100 * co
    q_pancreas <- fq_pancreas / 100 * co
    q_liver <- fq_liver / 100 * co
    q_lnode <- fq_lnode / 100 * co
    q_other <- fq_other / 100 * co
    q_ha <- fq_ha / 100 * co
    q_by <- fq_by / 100 * co

    # ====================================================================
    # 8. Lymph flows in L/h (PBPK_Population_Tissue.m, 'lymph flows'; Gill
    # et al. 2016); the remaining tissue takes the unlisted fraction.
    # ====================================================================
    ltot <- 0.00386 * WT * (1 + etaltot) # LT_Mean
    if (ltot < 0) {
      ltot <- 0.00386 * WT
    }
    l_lung <- 0.03 * ltot
    l_adipose <- 0.128 * ltot
    l_bone <- 0 * ltot
    l_brain <- 0.0105 * ltot
    l_gonads <- 0.013 * ltot
    l_heart <- 0.01 * ltot
    l_kidney <- 0.085 * ltot
    l_muscle <- 0.16 * ltot
    l_skin <- 0.073 * ltot
    l_thymus <- 0.011 * ltot
    l_gut <- 0.12 * ltot
    l_spleen <- 0 * ltot
    l_pancreas <- 0.003 * ltot
    l_liver <- 0.33 * ltot
    l_other <- (1 - (0.03 + 0.128 + 0 + 0.0105 + 0.013 + 0.01 + 0.085 + 0.16 + 0.073 +
      0.011 + 0.12 + 0 + 0.003 + 0.33)) * ltot

    # Plasma flow Porg = Qorg*(1 - HCT), QL = Qorg - Lorg, PL = Porg - Lorg
    # (PBPK_ODE_solution.m, 'Blood flows').
    ql_lung <- q_lung - l_lung
    ql_adipose <- q_adipose - l_adipose
    ql_bone <- q_bone - l_bone
    ql_brain <- q_brain - l_brain
    ql_gonads <- q_gonads - l_gonads
    ql_heart <- q_heart - l_heart
    ql_kidney <- q_kidney - l_kidney
    ql_muscle <- q_muscle - l_muscle
    ql_skin <- q_skin - l_skin
    ql_thymus <- q_thymus - l_thymus
    ql_gut <- q_gut - l_gut
    ql_spleen <- q_spleen - l_spleen
    ql_pancreas <- q_pancreas - l_pancreas
    ql_liver <- q_liver - l_liver
    ql_other <- q_other - l_other
    p_lung <- q_lung * (1 - hct)
    p_adipose <- q_adipose * (1 - hct)
    p_bone <- q_bone * (1 - hct)
    p_brain <- q_brain * (1 - hct)
    p_gonads <- q_gonads * (1 - hct)
    p_heart <- q_heart * (1 - hct)
    p_kidney <- q_kidney * (1 - hct)
    p_muscle <- q_muscle * (1 - hct)
    p_skin <- q_skin * (1 - hct)
    p_thymus <- q_thymus * (1 - hct)
    p_gut <- q_gut * (1 - hct)
    p_spleen <- q_spleen * (1 - hct)
    p_pancreas <- q_pancreas * (1 - hct)
    p_liver <- q_liver * (1 - hct)
    p_other <- q_other * (1 - hct)

    # ====================================================================
    # 9. Kidney and liver scalars
    # ====================================================================
    gfr <- exp(-0.0079 * AGE + 0.5 * bsa + 4.2) * (1 + etagfr) # GFR_Mean (mL/min), PBPK_Population_Tissue.m

    # Gen_LIPar (PBPK_Population_Liver.m): a draw outside [Min, Max] is
    # redrawn once from normrnd(Mean, (CV/1000)*Mean).
    mppgl_tv <- 10^(0.0000024 * AGE^3 - 0.00038 * AGE^2 + 0.0158 * AGE + 1.407) # MPPGL_Mean (mg/g liver), Barter 2008
    mppgl <- mppgl_tv * (1 + etamppgl)
    if (mppgl < 10 || mppgl > 110) {
      mppgl <- mppgl_tv * (1 + etamppgl_redraw)
    }
    ab_cyp3a4_liver <- 93.0 * (1 + etacyp3a4_liver) # CYP3A4_Mean (pmol/mg)
    if (ab_cyp3a4_liver < 18.6 || ab_cyp3a4_liver > 601) {
      ab_cyp3a4_liver <- 93.0 * (1 + etacyp3a4_liver_redraw)
    }
    ab_ugt1a1_liver <- 41.0 * (1 + etaugt1a1_liver) # UGT1A1_Mean (pmol/mg)
    if (ab_ugt1a1_liver < 4.0 || ab_ugt1a1_liver > 138) {
      ab_ugt1a1_liver <- 41.0 * (1 + etaugt1a1_liver_redraw)
    }

    # ====================================================================
    # 10. Gut transit and intestinal CYP3A4 (PBPK_Population_GIT.m).
    # Gen_GUPar redraws out-of-range values like Gen_LIPar.
    # ====================================================================
    t_stomach <- 0.25 + (1.00 - 0.25) * phi(etagastric) # GET ~ U(0.25, 1.00) h
    sitt <- 3.4 * (1 + etasitt) # SIT_Mean (h)
    if (sitt < 0.5 || sitt > 9.5) {
      sitt <- 3.4 * (1 + etasitt_redraw)
    }
    colt <- 17.2 * (1 + etacolt) # CNT_Mean (h)
    if (colt < 17.2 / 5 || colt > 17.2 * 5) {
      colt <- 17.2 * (1 + etacolt_redraw)
    }
    t_duodenum <- sitt * 0.091
    t_jejunum <- sitt * 0.426
    t_ileum <- sitt * 0.483
    t_colon <- colt
    ab_cyp3a4_gut <- 66.2 * (1 + etacyp3a4_gut) # intestinal CYP3A4_Mean (nmol)
    if (ab_cyp3a4_gut < 66.2 / 3 || ab_cyp3a4_gut > 66.2 * 3) {
      ab_cyp3a4_gut <- 66.2 * (1 + etacyp3a4_gut_redraw)
    }

    # Segment volumes (L): VsegCAT, luminal fluid VfluCAT, segment vascular
    # VvasCAT = 0.05*Vseg, segment interstitial VintCAT (share of the gut
    # Vint), uptake layer VlumCAT = Vseg - Vflu - Vvas - Vint.
    vseg_duodenum <- 0.060 * v_gut
    vseg_jejunum <- 0.279 * v_gut
    vseg_ileum <- 0.316 * v_gut
    vseg_colon <- 0.192 * v_gut
    vflu_duodenum <- 0.136 * vseg_duodenum
    vflu_jejunum <- 0.136 * vseg_jejunum
    vflu_ileum <- 0.136 * vseg_ileum
    vflu_colon <- 0.057 * vseg_colon
    vup_duodenum <- vseg_duodenum - vflu_duodenum - 0.05 * vseg_duodenum - 0.060 * ve_gut
    vup_jejunum <- vseg_jejunum - vflu_jejunum - 0.05 * vseg_jejunum - 0.279 * ve_gut
    vup_ileum <- vseg_ileum - vflu_ileum - 0.05 * vseg_ileum - 0.316 * ve_gut
    vup_colon <- vseg_colon - vflu_colon - 0.05 * vseg_colon - 0.192 * ve_gut

    # Lengths (cm; ICRP), cylinder radius and surface, enterocyte volume
    # VentCAT = 3040*(Surface/0.0001)*770*1e-15 (L), and the absorptive
    # surface SAB = Surface*FPC*FVI.
    len_duodenum <- 0.091 * (1.6 * HT)
    len_jejunum <- 0.426 * (1.6 * HT)
    len_ileum <- 0.483 * (1.6 * HT)
    len_colon <- 0.52 * HT + 18.5
    len_total <- len_duodenum + len_jejunum + len_ileum + len_colon
    surf_duodenum <- 2 * pi * sqrt((vseg_duodenum * 1000) / (pi * len_duodenum)) * len_duodenum
    surf_jejunum <- 2 * pi * sqrt((vseg_jejunum * 1000) / (pi * len_jejunum)) * len_jejunum
    surf_ileum <- 2 * pi * sqrt((vseg_ileum * 1000) / (pi * len_ileum)) * len_ileum
    surf_colon <- 2 * pi * sqrt((vseg_colon * 1000) / (pi * len_colon)) * len_colon
    vent_duodenum <- (3040 * (surf_duodenum / 0.0001)) * 770 * 1e-15
    vent_jejunum <- (3040 * (surf_jejunum / 0.0001)) * 770 * 1e-15
    vent_ileum <- (3040 * (surf_ileum / 0.0001)) * 770 * 1e-15
    vent_colon <- (3040 * (surf_colon / 0.0001)) * 770 * 1e-15

    # Absorption clearance (Drug/PBPK_Drug_absorption.m): Peff from Papp by
    # Sun et al. 2002, apportioned to the segments by length (Darwich 2010),
    # CLab = SAB*Peff*1e-4*3600*0.001 (L/h).
    peff <- 10^(0.6795 * log10(papp) - 0.3355)
    clab_duodenum <- surf_duodenum * 1.6 * 6.5 * peff * (len_duodenum / len_total) * 1e-4 * 3600 * 0.001
    clab_jejunum <- surf_jejunum * 1.6 * 8.6 * peff * (len_jejunum / len_total) * 1e-4 * 3600 * 0.001
    clab_ileum <- surf_ileum * 1.6 * 4.5 * peff * (len_ileum / len_total) * 1e-4 * 3600 * 0.001
    clab_colon <- surf_colon * 1.0 * 6.5 * peff * (len_colon / len_total) * 1e-4 * 3600 * 0.001 # FabsColon = 1
    kperup <- 1000 # DRUG.kPerUP = 0 in bictegravir.m, replaced by 1000 1/h (PBPK_Drug_PostProcessing.m)

    # ====================================================================
    # 11. Protein binding and partitioning (Drug/PBPK_Drug_distribution.m).
    # Monoprotic acid: KRio = 1 + 10^(pH - pKa). Albumin binder.
    # ====================================================================
    fup <- 1 / (1 + (((1 / fu) - 1) / 45.6) * hsa)
    krio_plasma <- 1 + 10^(7.4 - pka)
    krio_lung <- 1 + 10^(6.6 - pka)
    krio_710 <- 1 + 10^(7.1 - pka) # adipose, brain, heart
    krio_700 <- 1 + 10^(7.0 - pka) # bone, gonads, muscle, skin, thymus, gut, spleen, pancreas, lymph node, remaining
    krio_kidney <- 1 + 10^(7.22 - pka)
    krio_liver <- 1 + 10^(7.23 - pka)
    lipid_plasma <- (logp * 0.35 + (0.3 * logp + 0.7) * 0.225) / krio_plasma
    kapr <- ((1 / fup) - 1 - lipid_plasma) / hsa

    kpu_lung <- abs(krio_lung * 0.463 / krio_plasma + 0.348 + (logp * 0.003 + (0.3 * logp + 0.7) * 0.009) / krio_plasma + kapr * 0.212 * hsa)
    kpu_adipose <- abs(krio_710 * 0.039 / krio_plasma + 0.141 + (logp * 0.79 + (0.3 * logp + 0.7) * 0.002) / krio_plasma + kapr * 0.021 * hsa)
    kpu_bone <- abs(krio_700 * 0.341 / krio_plasma + 0.098 + (logp * 0.074 + (0.3 * logp + 0.7) * 0.0011) / krio_plasma + kapr * 0.1 * hsa)
    kpu_brain <- abs(krio_710 * 0.678 / krio_plasma + 0.092 + (logp * 0.051 + (0.3 * logp + 0.7) * 0.0565) / krio_plasma + kapr * 0.048 * hsa)
    kpu_gonads <- abs(krio_700 * 0.561 / krio_plasma + 0.239 + (logp * 0.007 + (0.3 * logp + 0.7) * 0.0077) / krio_plasma + kapr * 0.048 * hsa)
    kpu_heart <- abs(krio_710 * 0.445 / krio_plasma + 0.313 + (logp * 0.015 + (0.3 * logp + 0.7) * 0.0166) / krio_plasma + kapr * 0.157 * hsa)
    kpu_kidney <- abs(krio_kidney * 0.5 / krio_plasma + 0.283 + (logp * 0.0207 + (0.3 * logp + 0.7) * 0.0162) / krio_plasma + kapr * 0.13 * hsa)
    kpu_muscle <- abs(krio_700 * 0.669 / krio_plasma + 0.091 + (logp * 0.0238 + (0.3 * logp + 0.7) * 0.0072) / krio_plasma + kapr * 0.025 * hsa)
    kpu_skin <- abs(krio_700 * 0.0947 / krio_plasma + 0.623 + (logp * 0.0248 + (0.3 * logp + 0.7) * 0.0111) / krio_plasma + kapr * 0.277 * hsa)
    kpu_thymus <- abs(krio_700 * 0.626 / krio_plasma + 0.150 + (logp * 0.017 + (0.3 * logp + 0.7) * 0.0092) / krio_plasma + kapr * 0.075 * hsa)
    kpu_gut <- abs(krio_700 * 0.451 / krio_plasma + 0.267 + (logp * 0.0487 + (0.3 * logp + 0.7) * 0.0163) / krio_plasma + kapr * 0.158 * hsa)
    kpu_spleen <- abs(krio_700 * 0.58 / krio_plasma + 0.208 + (logp * 0.0201 + (0.3 * logp + 0.7) * 0.0198) / krio_plasma + kapr * 0.097 * hsa)
    kpu_pancreas <- abs(krio_700 * 0.664 / krio_plasma + 0.12 + (logp * 0.041 + (0.3 * logp + 0.7) * 0.0093) / krio_plasma + kapr * 0.06 * hsa)
    kpu_liver <- abs(krio_liver * 0.586 / krio_plasma + 0.165 + (logp * 0.0348 + (0.3 * logp + 0.7) * 0.0252) / krio_plasma + kapr * 0.086 * hsa)
    kpu_lnode <- abs(krio_700 * 0.58 / krio_plasma + 0.208 + (logp * 0.0201 + (0.3 * logp + 0.7) * 0.0198) / krio_plasma + kapr * 0.097 * hsa)
    kpu_other <- abs(krio_700 * fi_other / krio_plasma + fe_other + (logp * fnl_other + (0.3 * logp + 0.7) * fnp_other) / krio_plasma + kapr * kphsa_other * hsa)

    # Unbound fraction in interstitial water, fuint = 1/((KpHSA/FraEW)*(1/fup - 1) + 1).
    fuint_lung <- 1 / ((0.212 / 0.348) * ((1 / fup) - 1) + 1)
    fuint_adipose <- 1 / ((0.021 / 0.141) * ((1 / fup) - 1) + 1)
    fuint_bone <- 1 / ((0.1 / 0.098) * ((1 / fup) - 1) + 1)
    fuint_brain <- 1 / ((0.048 / 0.092) * ((1 / fup) - 1) + 1)
    fuint_gonads <- 1 / ((0.048 / 0.239) * ((1 / fup) - 1) + 1)
    fuint_heart <- 1 / ((0.157 / 0.313) * ((1 / fup) - 1) + 1)
    fuint_kidney <- 1 / ((0.13 / 0.283) * ((1 / fup) - 1) + 1)
    fuint_muscle <- 1 / ((0.025 / 0.091) * ((1 / fup) - 1) + 1)
    fuint_skin <- 1 / ((0.277 / 0.623) * ((1 / fup) - 1) + 1)
    fuint_thymus <- 1 / ((0.075 / 0.150) * ((1 / fup) - 1) + 1)
    fuint_spleen <- 1 / ((0.097 / 0.208) * ((1 / fup) - 1) + 1)
    fuint_pancreas <- 1 / ((0.06 / 0.12) * ((1 / fup) - 1) + 1)
    fuint_liver <- 1 / ((0.086 / 0.165) * ((1 / fup) - 1) + 1)
    fuint_lnode <- 1 / ((0.097 / 0.208) * ((1 / fup) - 1) + 1)
    fuint_other <- 1 / ((kphsa_other / fe_other) * ((1 / fup) - 1) + 1)

    # Unbound fraction in intracellular water,
    # fucel = 1/(1 + (logP*FraNL + (0.3*logP + 0.7)*FraNP)/KRio_plasma + KpHSA*HSA).
    fucel_lung <- 1 / (1 + (logp * 0.003 + (0.3 * logp + 0.7) * 0.009) / krio_plasma + 0.212 * hsa)
    fucel_adipose <- 1 / (1 + (logp * 0.79 + (0.3 * logp + 0.7) * 0.002) / krio_plasma + 0.021 * hsa)
    fucel_bone <- 1 / (1 + (logp * 0.074 + (0.3 * logp + 0.7) * 0.0011) / krio_plasma + 0.1 * hsa)
    fucel_brain <- 1 / (1 + (logp * 0.051 + (0.3 * logp + 0.7) * 0.0565) / krio_plasma + 0.048 * hsa)
    fucel_gonads <- 1 / (1 + (logp * 0.007 + (0.3 * logp + 0.7) * 0.0077) / krio_plasma + 0.048 * hsa)
    fucel_heart <- 1 / (1 + (logp * 0.015 + (0.3 * logp + 0.7) * 0.0166) / krio_plasma + 0.157 * hsa)
    fucel_kidney <- 1 / (1 + (logp * 0.0207 + (0.3 * logp + 0.7) * 0.0162) / krio_plasma + 0.13 * hsa)
    fucel_muscle <- 1 / (1 + (logp * 0.0238 + (0.3 * logp + 0.7) * 0.0072) / krio_plasma + 0.025 * hsa)
    fucel_skin <- 1 / (1 + (logp * 0.0248 + (0.3 * logp + 0.7) * 0.0111) / krio_plasma + 0.277 * hsa)
    fucel_thymus <- 1 / (1 + (logp * 0.017 + (0.3 * logp + 0.7) * 0.0092) / krio_plasma + 0.075 * hsa)
    fucel_gut <- 1 / (1 + (logp * 0.0487 + (0.3 * logp + 0.7) * 0.0163) / krio_plasma + 0.158 * hsa)
    fucel_spleen <- 1 / (1 + (logp * 0.0201 + (0.3 * logp + 0.7) * 0.0198) / krio_plasma + 0.097 * hsa)
    fucel_pancreas <- 1 / (1 + (logp * 0.041 + (0.3 * logp + 0.7) * 0.0093) / krio_plasma + 0.06 * hsa)
    fucel_liver <- 1 / (1 + (logp * 0.0348 + (0.3 * logp + 0.7) * 0.0252) / krio_plasma + 0.086 * hsa)
    fucel_lnode <- 1 / (1 + (logp * 0.0201 + (0.3 * logp + 0.7) * 0.0198) / krio_plasma + 0.097 * hsa)
    fucel_other <- 1 / (1 + (logp * fnl_other + (0.3 * logp + 0.7) * fnp_other) / krio_plasma + kphsa_other * hsa)

    # Membrane clearance CLin = CLout = ((Qorg - Lorg)/Kpu)/fup (L/h).
    clin_lung <- (ql_lung / kpu_lung) / fup
    clin_adipose <- (ql_adipose / kpu_adipose) / fup
    clin_bone <- (ql_bone / kpu_bone) / fup
    clin_brain <- (ql_brain / kpu_brain) / fup
    clin_gonads <- (ql_gonads / kpu_gonads) / fup
    clin_heart <- (ql_heart / kpu_heart) / fup
    clin_kidney <- (ql_kidney / kpu_kidney) / fup
    clin_muscle <- (ql_muscle / kpu_muscle) / fup
    clin_skin <- (ql_skin / kpu_skin) / fup
    clin_thymus <- (ql_thymus / kpu_thymus) / fup
    clin_gut <- (ql_gut / kpu_gut) / fup
    clin_spleen <- (ql_spleen / kpu_spleen) / fup
    clin_pancreas <- (ql_pancreas / kpu_pancreas) / fup
    clin_liver <- (ql_liver / kpu_liver) / fup
    clin_lnode <- (q_lnode / kpu_lnode) / fup
    clin_other <- (ql_other / kpu_other) / fup

    # Influx scalar Jtis = JinScalarAll*JinScalar (liver 2, others 1), and
    # vascular-to-interstitial scalar Fint = FinScalarAll*FinScalar with
    # FinScalarAll = 1. PBPK_Drug_PostProcessing.m reads FinScalar from rows
    # CompNo+2:3*CompNo of the stacked Kp/Jin/Fin array, one row too early,
    # so each tissue receives the next tissue's Jin scalar: the pancreas gets
    # the liver's 2.0 and every other tissue 1 (as-run).
    jtis_liver <- jin_all * jin_liver
    fint_pancreas <- jin_liver

    # ====================================================================
    # 12. Elimination (Drug/PBPK_Drug_elimination.m, PBPK_ODE_solution.m).
    # CLint uL/min/pmol (or /mg) * 60e-6 -> L/h/pmol (or /mg), scaled by
    # abundance, MPPGL and liver weight in g; applied to the unbound
    # intracellular liver concentration.
    # ====================================================================
    clint_liver <- (clint_cyp3a4 * ab_cyp3a4_liver + clint_ugt1a1 * ab_ugt1a1_liver + clint_hep) *
      60e-6 * mppgl * (w_liver * 1000)
    # Intestinal CYP3A4 split 0.136 / 0.544 / 0.320 over duodenum, jejunum,
    # ileum (Paine 1997); nmol -> pmol by 1e3.
    clmet_duodenum <- clint_cyp3a4 * 60e-6 * 0.136 * ab_cyp3a4_gut * 1000
    clmet_jejunum <- clint_cyp3a4 * 60e-6 * 0.544 * ab_cyp3a4_gut * 1000
    clmet_ileum <- clint_cyp3a4 * 60e-6 * 0.320 * ab_cyp3a4_gut * 1000
    # Renal clearance on the kidney vascular (blood) concentration, scaled
    # by GFR relative to 130 (men) or 120 (women) mL/min and by fup/fu.
    cl_renal <- exp(lcl_renal)
    cl_r <- cl_renal * (gfr / (130 - 10 * SEXF)) * (fup / fu)

    # ====================================================================
    # 13. Concentrations (mg/L). Vascular states hold blood concentrations
    # (the vascular-to-interstitial flux divides them by BP).
    # ====================================================================
    c_venous <- venous / v_venous
    c_arterial <- arterial / v_arterial
    cv_lung <- lung_vas / vv_lung
    ce_lung <- lung_ew / ve_lung
    ci_lung <- lung_iw / vi_lung
    cv_adipose <- adipose_vas / vv_adipose
    ce_adipose <- adipose_ew / ve_adipose
    ci_adipose <- adipose_iw / vi_adipose
    cv_bone <- bone_vas / vv_bone
    ce_bone <- bone_ew / ve_bone
    ci_bone <- bone_iw / vi_bone
    cv_brain <- brain_vas / vv_brain
    ce_brain <- brain_ew / ve_brain
    ci_brain <- brain_iw / vi_brain
    cv_gonads <- gonads_vas / vv_gonads
    ce_gonads <- gonads_ew / ve_gonads
    ci_gonads <- gonads_iw / vi_gonads
    cv_heart <- heart_vas / vv_heart
    ce_heart <- heart_ew / ve_heart
    ci_heart <- heart_iw / vi_heart
    cv_kidney <- kidney_vas / vv_kidney
    ce_kidney <- kidney_ew / ve_kidney
    ci_kidney <- kidney_iw / vi_kidney
    cv_muscle <- muscle_vas / vv_muscle
    ce_muscle <- muscle_ew / ve_muscle
    ci_muscle <- muscle_iw / vi_muscle
    cv_skin <- skin_vas / vv_skin
    ce_skin <- skin_ew / ve_skin
    ci_skin <- skin_iw / vi_skin
    cv_thymus <- thymus_vas / vv_thymus
    ce_thymus <- thymus_ew / ve_thymus
    ci_thymus <- thymus_iw / vi_thymus
    cv_gut <- gut_vas / vv_gut
    ce_gut <- gut_ew / ve_gut
    cv_spleen <- spleen_vas / vv_spleen
    ce_spleen <- spleen_ew / ve_spleen
    ci_spleen <- spleen_iw / vi_spleen
    cv_pancreas <- pancreas_vas / vv_pancreas
    ce_pancreas <- pancreas_ew / ve_pancreas
    ci_pancreas <- pancreas_iw / vi_pancreas
    cv_liver <- liver_vas / vv_liver
    ce_liver <- liver_ew / ve_liver
    ci_liver <- liver_iw / vi_liver
    cv_lnode <- lnode_vas / vv_lnode
    ce_lnode <- lnode_ew / ve_lnode
    ci_lnode <- lnode_iw / vi_lnode
    cv_other <- other_vas / vv_other
    ce_other <- other_ew / ve_other
    ci_other <- other_iw / vi_other
    cent_duodenum <- duodenum_enterocyte / vent_duodenum
    cent_jejunum <- jejunum_enterocyte / vent_jejunum
    cent_ileum <- ileum_enterocyte / vent_ileum
    cent_colon <- colon_enterocyte / vent_colon

    # ====================================================================
    # 14. ODEs (PBPK_ODE_solution.m, rhs_function) in amount form: each
    # Matlab dC/dt = (1/V)*(flux) becomes dA/dt = flux. Vascular (_vas),
    # interstitial (_ew) and cellular (_iw) sub-compartments.
    # ====================================================================
    d/dt(lung_vas) <- q_lung * c_venous - ql_lung * cv_lung - p_lung * cv_lung / bp + (p_lung - l_lung) * ce_lung
    d/dt(lung_ew) <- p_lung * cv_lung / bp - (p_lung - l_lung) * ce_lung - l_lung * ce_lung -
      clin_lung * jin_all * ce_lung * fuint_lung + clin_lung * ci_lung * fucel_lung
    d/dt(lung_iw) <- clin_lung * jin_all * ce_lung * fuint_lung - clin_lung * ci_lung * fucel_lung

    d/dt(adipose_vas) <- q_adipose * c_arterial - ql_adipose * cv_adipose - p_adipose * cv_adipose / bp + (p_adipose - l_adipose) * ce_adipose
    d/dt(adipose_ew) <- p_adipose * cv_adipose / bp - (p_adipose - l_adipose) * ce_adipose - l_adipose * ce_adipose -
      clin_adipose * jin_all * ce_adipose * fuint_adipose + clin_adipose * ci_adipose * fucel_adipose
    d/dt(adipose_iw) <- clin_adipose * jin_all * ce_adipose * fuint_adipose - clin_adipose * ci_adipose * fucel_adipose

    d/dt(bone_vas) <- q_bone * c_arterial - ql_bone * cv_bone - p_bone * cv_bone / bp + (p_bone - l_bone) * ce_bone
    d/dt(bone_ew) <- p_bone * cv_bone / bp - (p_bone - l_bone) * ce_bone - l_bone * ce_bone -
      clin_bone * jin_all * ce_bone * fuint_bone + clin_bone * ci_bone * fucel_bone
    d/dt(bone_iw) <- clin_bone * jin_all * ce_bone * fuint_bone - clin_bone * ci_bone * fucel_bone

    d/dt(brain_vas) <- q_brain * c_arterial - ql_brain * cv_brain - p_brain * cv_brain / bp + (p_brain - l_brain) * ce_brain
    d/dt(brain_ew) <- p_brain * cv_brain / bp - (p_brain - l_brain) * ce_brain - l_brain * ce_brain -
      clin_brain * jin_all * ce_brain * fuint_brain + clin_brain * ci_brain * fucel_brain
    d/dt(brain_iw) <- clin_brain * jin_all * ce_brain * fuint_brain - clin_brain * ci_brain * fucel_brain

    d/dt(gonads_vas) <- q_gonads * c_arterial - ql_gonads * cv_gonads - p_gonads * cv_gonads / bp + (p_gonads - l_gonads) * ce_gonads
    d/dt(gonads_ew) <- p_gonads * cv_gonads / bp - (p_gonads - l_gonads) * ce_gonads - l_gonads * ce_gonads -
      clin_gonads * jin_all * ce_gonads * fuint_gonads + clin_gonads * ci_gonads * fucel_gonads
    d/dt(gonads_iw) <- clin_gonads * jin_all * ce_gonads * fuint_gonads - clin_gonads * ci_gonads * fucel_gonads

    d/dt(heart_vas) <- q_heart * c_arterial - ql_heart * cv_heart - p_heart * cv_heart / bp + (p_heart - l_heart) * ce_heart
    d/dt(heart_ew) <- p_heart * cv_heart / bp - (p_heart - l_heart) * ce_heart - l_heart * ce_heart -
      clin_heart * jin_all * ce_heart * fuint_heart + clin_heart * ci_heart * fucel_heart
    d/dt(heart_iw) <- clin_heart * jin_all * ce_heart * fuint_heart - clin_heart * ci_heart * fucel_heart

    d/dt(kidney_vas) <- q_kidney * c_arterial - ql_kidney * cv_kidney - p_kidney * cv_kidney / bp + (p_kidney - l_kidney) * ce_kidney - cl_r * cv_kidney
    d/dt(kidney_ew) <- p_kidney * cv_kidney / bp - (p_kidney - l_kidney) * ce_kidney - l_kidney * ce_kidney -
      clin_kidney * jin_all * ce_kidney * fuint_kidney + clin_kidney * ci_kidney * fucel_kidney
    d/dt(kidney_iw) <- clin_kidney * jin_all * ce_kidney * fuint_kidney - clin_kidney * ci_kidney * fucel_kidney

    d/dt(muscle_vas) <- q_muscle * c_arterial - ql_muscle * cv_muscle - p_muscle * cv_muscle / bp + (p_muscle - l_muscle) * ce_muscle
    d/dt(muscle_ew) <- p_muscle * cv_muscle / bp - (p_muscle - l_muscle) * ce_muscle - l_muscle * ce_muscle -
      clin_muscle * jin_all * ce_muscle * fuint_muscle + clin_muscle * ci_muscle * fucel_muscle
    d/dt(muscle_iw) <- clin_muscle * jin_all * ce_muscle * fuint_muscle - clin_muscle * ci_muscle * fucel_muscle

    d/dt(skin_vas) <- q_skin * c_arterial - ql_skin * cv_skin - p_skin * cv_skin / bp + (p_skin - l_skin) * ce_skin
    d/dt(skin_ew) <- p_skin * cv_skin / bp - (p_skin - l_skin) * ce_skin - l_skin * ce_skin -
      clin_skin * jin_all * ce_skin * fuint_skin + clin_skin * ci_skin * fucel_skin
    d/dt(skin_iw) <- clin_skin * jin_all * ce_skin * fuint_skin - clin_skin * ci_skin * fucel_skin

    d/dt(thymus_vas) <- q_thymus * c_arterial - ql_thymus * cv_thymus - p_thymus * cv_thymus / bp + (p_thymus - l_thymus) * ce_thymus
    d/dt(thymus_ew) <- p_thymus * cv_thymus / bp - (p_thymus - l_thymus) * ce_thymus - l_thymus * ce_thymus -
      clin_thymus * jin_all * ce_thymus * fuint_thymus + clin_thymus * ci_thymus * fucel_thymus
    d/dt(thymus_iw) <- clin_thymus * jin_all * ce_thymus * fuint_thymus - clin_thymus * ci_thymus * fucel_thymus

    # Gut tissue: vascular and interstitial only; the enterocytes of the four
    # CAT segments release drug into the gut interstitium (GUPermScalar = 1).
    d/dt(gut_vas) <- q_gut * c_arterial - ql_gut * cv_gut - p_gut * cv_gut / bp + (p_gut - l_gut) * ce_gut
    d/dt(gut_ew) <- p_gut * cv_gut / bp - (p_gut - l_gut) * ce_gut - l_gut * ce_gut +
      clin_gut * fucel_gut * (cent_duodenum + cent_jejunum + cent_ileum + cent_colon)

    d/dt(spleen_vas) <- q_spleen * c_arterial - ql_spleen * cv_spleen - p_spleen * cv_spleen / bp + (p_spleen - l_spleen) * ce_spleen
    d/dt(spleen_ew) <- p_spleen * cv_spleen / bp - (p_spleen - l_spleen) * ce_spleen - l_spleen * ce_spleen -
      clin_spleen * jin_all * ce_spleen * fuint_spleen + clin_spleen * ci_spleen * fucel_spleen
    d/dt(spleen_iw) <- clin_spleen * jin_all * ce_spleen * fuint_spleen - clin_spleen * ci_spleen * fucel_spleen

    d/dt(pancreas_vas) <- q_pancreas * c_arterial - ql_pancreas * cv_pancreas - p_pancreas * fint_pancreas * cv_pancreas / bp + (p_pancreas - l_pancreas) * ce_pancreas
    d/dt(pancreas_ew) <- p_pancreas * fint_pancreas * cv_pancreas / bp - (p_pancreas - l_pancreas) * ce_pancreas - l_pancreas * ce_pancreas -
      clin_pancreas * jin_all * ce_pancreas * fuint_pancreas + clin_pancreas * ci_pancreas * fucel_pancreas
    d/dt(pancreas_iw) <- clin_pancreas * jin_all * ce_pancreas * fuint_pancreas - clin_pancreas * ci_pancreas * fucel_pancreas

    # Liver: portal outflow of gut, spleen and pancreas plus hepatic-artery
    # and bypass flow; metabolism on the unbound intracellular concentration
    # (hepatic transporters switched off, COMP.LItrans = OFF; no biliary
    # clearance).
    d/dt(liver_vas) <- ql_gut * cv_gut + ql_spleen * cv_spleen + ql_pancreas * cv_pancreas + (q_ha + q_by) * c_arterial -
      ql_liver * cv_liver - p_liver * cv_liver / bp + (p_liver - l_liver) * ce_liver
    d/dt(liver_ew) <- p_liver * cv_liver / bp - (p_liver - l_liver) * ce_liver - l_liver * ce_liver -
      clin_liver * jtis_liver * ce_liver * fuint_liver + clin_liver * ci_liver * fucel_liver
    d/dt(liver_iw) <- clin_liver * jtis_liver * ce_liver * fuint_liver - clin_liver * ci_liver * fucel_liver -
      clint_liver * ci_liver * fucel_liver

    # Lymph node: collects every tissue's interstitial lymph and returns it
    # through its vascular space to venous blood.
    d/dt(lnode_vas) <- q_lnode * c_arterial - q_lnode * cv_lnode + ltot * ce_lnode - ltot * cv_lnode
    d/dt(lnode_ew) <- l_lung * ce_lung + l_adipose * ce_adipose + l_bone * ce_bone + l_brain * ce_brain +
      l_gonads * ce_gonads + l_heart * ce_heart + l_kidney * ce_kidney + l_muscle * ce_muscle +
      l_skin * ce_skin + l_thymus * ce_thymus + l_gut * ce_gut + l_spleen * ce_spleen +
      l_pancreas * ce_pancreas + l_liver * ce_liver + l_other * ce_other - ltot * ce_lnode -
      clin_lnode * jin_all * ce_lnode * fuint_lnode + clin_lnode * ci_lnode * fucel_lnode
    d/dt(lnode_iw) <- clin_lnode * jin_all * ce_lnode * fuint_lnode - clin_lnode * ci_lnode * fucel_lnode

    d/dt(other_vas) <- q_other * c_arterial - ql_other * cv_other - p_other * cv_other / bp + (p_other - l_other) * ce_other
    d/dt(other_ew) <- p_other * cv_other / bp - (p_other - l_other) * ce_other - l_other * ce_other -
      clin_other * jin_all * ce_other * fuint_other + clin_other * ci_other * fucel_other
    d/dt(other_iw) <- clin_other * jin_all * ce_other * fuint_other - clin_other * ci_other * fucel_other

    # Venous and arterial blood pools.
    d/dt(venous) <- ql_adipose * cv_adipose + ql_bone * cv_bone + ql_brain * cv_brain +
      ql_gonads * cv_gonads + ql_heart * cv_heart + ql_kidney * cv_kidney +
      ql_muscle * cv_muscle + ql_skin * cv_skin + ql_thymus * cv_thymus +
      ql_liver * cv_liver + q_lnode * cv_lnode + ql_other * cv_other +
      ltot * cv_lnode - q_lung * c_venous
    d/dt(arterial) <- ql_lung * cv_lung - (q_adipose + q_bone + q_brain + q_gonads + q_heart +
      q_kidney + q_muscle + q_skin + q_thymus + q_gut + q_spleen + q_pancreas +
      q_ha + q_by + q_lnode + q_other) * c_arterial

    # CAT gut (hepatic and intestinal transporters off). The oral dose goes
    # into `stomach`; the 200 mL of water taken with it only dilutes the
    # stomach concentration and cancels in amount form.
    d/dt(stomach) <- -stomach / t_stomach
    d/dt(duodenum) <- stomach / t_stomach - duodenum / t_duodenum - clab_duodenum * duodenum / vflu_duodenum
    d/dt(jejunum) <- duodenum / t_duodenum - jejunum / t_jejunum - clab_jejunum * jejunum / vflu_jejunum
    d/dt(ileum) <- jejunum / t_jejunum - ileum / t_ileum - clab_ileum * ileum / vflu_ileum
    d/dt(colon) <- ileum / t_ileum - colon / t_colon - clab_colon * colon / vflu_colon
    d/dt(duodenum_uptake) <- clab_duodenum * duodenum / vflu_duodenum - kperup * duodenum_uptake
    d/dt(jejunum_uptake) <- clab_jejunum * jejunum / vflu_jejunum - kperup * jejunum_uptake
    d/dt(ileum_uptake) <- clab_ileum * ileum / vflu_ileum - kperup * ileum_uptake
    d/dt(colon_uptake) <- clab_colon * colon / vflu_colon - kperup * colon_uptake
    d/dt(duodenum_enterocyte) <- kperup * duodenum_uptake - clin_gut * cent_duodenum * fucel_gut -
      clmet_duodenum * cent_duodenum * fucel_gut
    d/dt(jejunum_enterocyte) <- kperup * jejunum_uptake - clin_gut * cent_jejunum * fucel_gut -
      clmet_jejunum * cent_jejunum * fucel_gut
    d/dt(ileum_enterocyte) <- kperup * ileum_uptake - clin_gut * cent_ileum * fucel_gut -
      clmet_ileum * cent_ileum * fucel_gut
    d/dt(colon_enterocyte) <- kperup * colon_uptake - clin_gut * cent_colon * fucel_gut
    d/dt(a_feces) <- colon / t_colon

    # Reported plasma concentration (ng/mL): PBPK_ExtractConcentration.m
    # takes the venous state (CONC.VB) times the molecular weight, with no
    # blood:plasma conversion.
    Cc <- 1000 * c_venous
  })
}
