Chu_2016_vegf165b_human_qsp <- function() {
  description <- paste(
    "QSP. Human three-compartment (normal tissue, blood, ischemic PAD calf",
    "muscle) model of VEGF isoform distribution in peripheral arterial",
    "disease, adding the anti-angiogenic splice variant VEGF165b to the",
    "Popel-lab whole-body VEGF model. Tracks free VEGF165a, VEGF165b and",
    "VEGF121; their binding to VEGFR1, VEGFR2 and neuropilin-1 on abluminal",
    "and luminal endothelial surfaces, to GAG sites in the ECM, endothelial",
    "and parenchymal basement membranes, to soluble VEGFR1 and to",
    "alpha-2-macroglobulin (native and fast form) in plasma; receptor",
    "insertion, coupling and internalization; secretion, degradation,",
    "lymphatic drainage, microvascular permeability and plasma clearance.",
    "152 ODE states and 410 reactions translated from the authors' SimBiology",
    "SBML export (Supplementary File 2). Endogenous steady-state model with",
    "no dosing, no between-subject variability and no residual error.",
    sep = " "
  )
  reference <- paste(
    "Chu LH, Ganta VC, Choi MH, Chen G, Finley SD, Annex BH, Popel AS. A",
    "multiscale computational model predicts distribution of anti-angiogenic",
    "isoform VEGF165b in peripheral arterial disease in human and mouse. Sci",
    "Rep. 2016;6:37030. doi:10.1038/srep37030.",
    sep = " "
  )
  vignette <- "Chu_2016_vegf165b"

  units <- list(
    time = "h",
    dosing = "none (endogenous steady-state model)",
    concentration = paste(
      "states in mol per cm3 of compartment (tissue or blood); outputs ending",
      "in _pM are free-fluid concentrations in pmol/L (interstitial fluid for",
      "the tissues, plasma for blood); outputs ending in _num are receptors",
      "per cell",
      sep = " "
    )
  )

  covariateData <- list()

  # One entry per ODE state, named after the SBML species (Supplementary File 2).
  compartmentData <- list(
    bld_V165 = list(analyte = "free VEGF165a", units = "mol/cm3 blood", specimen = "whole blood", verified = TRUE),
    bld_V121 = list(analyte = "free VEGF121", units = "mol/cm3 blood", specimen = "whole blood", verified = TRUE),
    nrm_V165 = list(
      analyte = "free VEGF165a",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V121 = list(
      analyte = "free VEGF121",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    bld_R1_L_norm_N1_L_norm = list(
      analyte = "luminal VEGFR1 (normal-tissue endothelium):luminal NRP1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R1_L_norm = list(
      analyte = "luminal VEGFR1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_N1_L_norm = list(
      analyte = "luminal NRP1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    nrm_R1 = list(
      analyte = "abluminal VEGFR1",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_N1 = list(
      analyte = "abluminal NRP1",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_R1_N1 = list(
      analyte = "VEGFR1:NRP1",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_R2 = list(
      analyte = "abluminal VEGFR2",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    bld_R2_L_norm = list(
      analyte = "luminal VEGFR2 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    dis_V165 = list(
      analyte = "free VEGF165a",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_MECM = list(
      analyte = "ECM GAG site",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V165_MECM = list(
      analyte = "VEGF165a:ECM GAG site",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_MEBM = list(
      analyte = "EBM GAG site",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V165_MEBM = list(
      analyte = "VEGF165a:EBM GAG site",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_MPBM = list(
      analyte = "PBM GAG site",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V165_MPBM = list(
      analyte = "VEGF165a:PBM GAG site",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_R1 = list(
      analyte = "abluminal VEGFR1",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V165_R1 = list(
      analyte = "VEGF165a:VEGFR1",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_R2 = list(
      analyte = "abluminal VEGFR2",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V165_R2 = list(
      analyte = "VEGF165a:VEGFR2",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_N1 = list(
      analyte = "abluminal NRP1",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V165_N1 = list(
      analyte = "VEGF165a:NRP1",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V121 = list(
      analyte = "free VEGF121",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V121_R1 = list(
      analyte = "VEGF121:VEGFR1",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V121_R2 = list(
      analyte = "VEGF121:VEGFR2",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V165_R2_N1 = list(
      analyte = "VEGF165a:VEGFR2:NRP1",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_R2_V121 = list(
      analyte = "VEGFR2:VEGF121",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V121_R1_N1 = list(
      analyte = "VEGF121:VEGFR1:NRP1",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_R1_N1 = list(
      analyte = "VEGFR1:NRP1",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    bld_R1_L_dis = list(
      analyte = "luminal VEGFR1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_N1_L_dis = list(
      analyte = "luminal NRP1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R1_L_dis_N1_L_dis = list(
      analyte = "luminal VEGFR1 (PAD calf endothelium):luminal NRP1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R2_L_dis = list(
      analyte = "luminal VEGFR2 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V165_R1_L_dis = list(
      analyte = "VEGF165a:luminal VEGFR1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V165_R2_L_dis = list(
      analyte = "VEGF165a:luminal VEGFR2 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V165_N1_L_dis = list(
      analyte = "VEGF165a:luminal NRP1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V121_R1_L_dis_N1_L_dis = list(
      analyte = "VEGF121:luminal VEGFR1 (PAD calf endothelium):luminal NRP1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R1_L_dis_V121 = list(
      analyte = "luminal VEGFR1 (PAD calf endothelium):VEGF121",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R2_L_dis_V121 = list(
      analyte = "luminal VEGFR2 (PAD calf endothelium):VEGF121",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R2_L_dis_V165_N1_L_dis = list(
      analyte = "luminal VEGFR2 (PAD calf endothelium):VEGF165a:luminal NRP1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V165_R1_L_norm = list(
      analyte = "VEGF165a:luminal VEGFR1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V165_R2_L_norm = list(
      analyte = "VEGF165a:luminal VEGFR2 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V165_N1_L_norm = list(
      analyte = "VEGF165a:luminal NRP1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V121_R1_L_norm_N1_L_norm = list(
      analyte = "VEGF121:luminal VEGFR1 (normal-tissue endothelium):luminal NRP1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R1_L_norm_V121 = list(
      analyte = "luminal VEGFR1 (normal-tissue endothelium):VEGF121",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R2_L_norm_V121 = list(
      analyte = "luminal VEGFR2 (normal-tissue endothelium):VEGF121",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R2_L_norm_V165_N1_L_norm = list(
      analyte = "luminal VEGFR2 (normal-tissue endothelium):VEGF165a:luminal NRP1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    nrm_MECM = list(
      analyte = "ECM GAG site",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V165_MECM = list(
      analyte = "VEGF165a:ECM GAG site",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_MEBM = list(
      analyte = "EBM GAG site",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V165_MEBM = list(
      analyte = "VEGF165a:EBM GAG site",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_MPBM = list(
      analyte = "PBM GAG site",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V165_MPBM = list(
      analyte = "VEGF165a:PBM GAG site",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V165_R1 = list(
      analyte = "VEGF165a:VEGFR1",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V165_R2 = list(
      analyte = "VEGF165a:VEGFR2",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V165_N1 = list(
      analyte = "VEGF165a:NRP1",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V121_R1 = list(
      analyte = "VEGF121:VEGFR1",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V121_R2 = list(
      analyte = "VEGF121:VEGFR2",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V121_R1_N1 = list(
      analyte = "VEGF121:VEGFR1:NRP1",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V165_R2_N1 = list(
      analyte = "VEGF165a:VEGFR2:NRP1",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_R1_N2 = list(
      analyte = "VEGFR1:NRP1 (orphan N2 species)",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V121_R1_N2 = list(
      analyte = "VEGF121:VEGFR1:NRP1 (orphan N2 species)",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_sR1 = list(
      analyte = "free sVEGFR1",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_sR1 = list(
      analyte = "free sVEGFR1",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_sR1_MECM = list(
      analyte = "sVEGFR1:ECM GAG site",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_sR1_MEBM = list(
      analyte = "sVEGFR1:EBM GAG site",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_sR1_MPBM = list(
      analyte = "sVEGFR1:PBM GAG site",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V165_sR1 = list(
      analyte = "VEGF165a:sVEGFR1",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V121_sR1 = list(
      analyte = "VEGF121:sVEGFR1",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_sR1_N1 = list(
      analyte = "sVEGFR1:NRP1",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V121_sR1_N1 = list(
      analyte = "VEGF121:sVEGFR1:NRP1",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    bld_sR1 = list(analyte = "free sVEGFR1", units = "mol/cm3 blood", specimen = "whole blood", verified = TRUE),
    bld_V165_sR1 = list(analyte = "VEGF165a:sVEGFR1", units = "mol/cm3 blood", specimen = "whole blood", verified = TRUE),
    bld_V121_sR1 = list(analyte = "VEGF121:sVEGFR1", units = "mol/cm3 blood", specimen = "whole blood", verified = TRUE),
    bld_sR1_N1_L_norm = list(
      analyte = "sVEGFR1:luminal NRP1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_sR1_N1_L_dis = list(
      analyte = "sVEGFR1:luminal NRP1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V121_sR1_N1_L_norm = list(
      analyte = "VEGF121:sVEGFR1:luminal NRP1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V121_sR1_N1_L_dis = list(
      analyte = "VEGF121:sVEGFR1:luminal NRP1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    dis_sR1_MECM = list(
      analyte = "sVEGFR1:ECM GAG site",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_sR1_MEBM = list(
      analyte = "sVEGFR1:EBM GAG site",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_sR1_MPBM = list(
      analyte = "sVEGFR1:PBM GAG site",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V165_sR1 = list(
      analyte = "VEGF165a:sVEGFR1",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V121_sR1 = list(
      analyte = "VEGF121:sVEGFR1",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_sR1_N1 = list(
      analyte = "sVEGFR1:NRP1",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V121_sR1_N1 = list(
      analyte = "VEGF121:sVEGFR1:NRP1",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    bld_a2M = list(analyte = "alpha-2-macroglobulin", units = "mol/cm3 blood", specimen = "whole blood", verified = TRUE),
    bld_V165_a2M = list(
      analyte = "VEGF165a:alpha-2-macroglobulin",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V121_a2M = list(
      analyte = "VEGF121:alpha-2-macroglobulin",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_V165_R1_L_norm = list(
      analyte = "alpha-2-macroglobulin:VEGF165a:luminal VEGFR1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_V165_R2_L_norm = list(
      analyte = "alpha-2-macroglobulin:VEGF165a:luminal VEGFR2 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_V165_N1_L_norm = list(
      analyte = "alpha-2-macroglobulin:VEGF165a:luminal NRP1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R2_L_norm_a2M_V165_N1_L_norm = list(
      analyte = "luminal VEGFR2 (normal-tissue endothelium):alpha-2-macroglobulin:VEGF165a:luminal NRP1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R1_L_norm_V121_a2M = list(
      analyte = "luminal VEGFR1 (normal-tissue endothelium):VEGF121:alpha-2-macroglobulin",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R2_L_norm_V121_a2M = list(
      analyte = "luminal VEGFR2 (normal-tissue endothelium):VEGF121:alpha-2-macroglobulin",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_V121_R1_L_norm_N1_L_norm = list(
      analyte = "alpha-2-macroglobulin:VEGF121:luminal VEGFR1 (normal-tissue endothelium):luminal NRP1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_V165_R1_L_dis = list(
      analyte = "alpha-2-macroglobulin:VEGF165a:luminal VEGFR1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_V165_R2_L_dis = list(
      analyte = "alpha-2-macroglobulin:VEGF165a:luminal VEGFR2 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_V165_N1_L_dis = list(
      analyte = "alpha-2-macroglobulin:VEGF165a:luminal NRP1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R2_L_dis_a2M_V165_N1_L_dis = list(
      analyte = "luminal VEGFR2 (PAD calf endothelium):alpha-2-macroglobulin:VEGF165a:luminal NRP1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R1_L_dis_V121_a2M = list(
      analyte = "luminal VEGFR1 (PAD calf endothelium):VEGF121:alpha-2-macroglobulin",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R2_L_dis_V121_a2M = list(
      analyte = "luminal VEGFR2 (PAD calf endothelium):VEGF121:alpha-2-macroglobulin",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_V121_R1_L_dis_N1_L_dis = list(
      analyte = "alpha-2-macroglobulin:VEGF121:luminal VEGFR1 (PAD calf endothelium):luminal NRP1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_fast = list(
      analyte = "activated alpha-2-macroglobulin",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V165_a2M_fast = list(
      analyte = "VEGF165a:activated alpha-2-macroglobulin",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V121_a2M_fast = list(
      analyte = "VEGF121:activated alpha-2-macroglobulin",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_fast_V165_R1_L_norm = list(
      analyte = "activated alpha-2-macroglobulin:VEGF165a:luminal VEGFR1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_fast_V165_R2_L_norm = list(
      analyte = "activated alpha-2-macroglobulin:VEGF165a:luminal VEGFR2 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_fast_V165_N1_L_norm = list(
      analyte = "activated alpha-2-macroglobulin:VEGF165a:luminal NRP1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R2_L_norm_a2M_fast_V165_N1_L_norm = list(
      analyte = "luminal VEGFR2 (normal-tissue endothelium):activated alpha-2-macroglobulin:VEGF165a:luminal NRP1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R1_L_norm_V121_a2M_fast = list(
      analyte = "luminal VEGFR1 (normal-tissue endothelium):VEGF121:activated alpha-2-macroglobulin",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R2_L_norm_V121_a2M_fast = list(
      analyte = "luminal VEGFR2 (normal-tissue endothelium):VEGF121:activated alpha-2-macroglobulin",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_fast_V121_R1_L_norm_N1_L_norm = list(
      analyte = "activated alpha-2-macroglobulin:VEGF121:luminal VEGFR1 (normal-tissue endothelium):luminal NRP1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_fast_V165_R1_L_dis = list(
      analyte = "activated alpha-2-macroglobulin:VEGF165a:luminal VEGFR1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_fast_V165_R2_L_dis = list(
      analyte = "activated alpha-2-macroglobulin:VEGF165a:luminal VEGFR2 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_fast_V165_N1_L_dis = list(
      analyte = "activated alpha-2-macroglobulin:VEGF165a:luminal NRP1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R2_L_dis_a2M_fast_V165_N1_L_dis = list(
      analyte = "luminal VEGFR2 (PAD calf endothelium):activated alpha-2-macroglobulin:VEGF165a:luminal NRP1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R1_L_dis_V121_a2M_fast = list(
      analyte = "luminal VEGFR1 (PAD calf endothelium):VEGF121:activated alpha-2-macroglobulin",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_R2_L_dis_V121_a2M_fast = list(
      analyte = "luminal VEGFR2 (PAD calf endothelium):VEGF121:activated alpha-2-macroglobulin",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_fast_V121_R1_L_dis_N1_L_dis = list(
      analyte = "activated alpha-2-macroglobulin:VEGF121:luminal VEGFR1 (PAD calf endothelium):luminal NRP1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    nrm_V165b = list(
      analyte = "free VEGF165b",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V165b_MPBM = list(
      analyte = "VEGF165b:PBM GAG site",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V165b_MECM = list(
      analyte = "VEGF165b:ECM GAG site",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V165b_MEBM = list(
      analyte = "VEGF165b:EBM GAG site",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V165b = list(
      analyte = "free VEGF165b",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V165b_MPBM = list(
      analyte = "VEGF165b:PBM GAG site",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V165b_MECM = list(
      analyte = "VEGF165b:ECM GAG site",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V165b_MEBM = list(
      analyte = "VEGF165b:EBM GAG site",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    bld_V165b = list(analyte = "free VEGF165b", units = "mol/cm3 blood", specimen = "whole blood", verified = TRUE),
    bld_V165b_R1_L_norm = list(
      analyte = "VEGF165b:luminal VEGFR1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V165b_R2_L_norm = list(
      analyte = "VEGF165b:luminal VEGFR2 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V165b_sR1 = list(analyte = "VEGF165b:sVEGFR1", units = "mol/cm3 blood", specimen = "whole blood", verified = TRUE),
    nrm_V165b_R1 = list(
      analyte = "VEGF165b:VEGFR1",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V165b_R2 = list(
      analyte = "VEGF165b:VEGFR2",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    nrm_V165b_sR1 = list(
      analyte = "VEGF165b:sVEGFR1",
      units = "mol/cm3 normal tissue",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V165b_R1 = list(
      analyte = "VEGF165b:VEGFR1",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V165b_R2 = list(
      analyte = "VEGF165b:VEGFR2",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    dis_V165b_sR1 = list(
      analyte = "VEGF165b:sVEGFR1",
      units = "mol/cm3 PAD calf muscle",
      specimen = "tissue",
      verified = TRUE
    ),
    bld_V165b_R1_L_dis = list(
      analyte = "VEGF165b:luminal VEGFR1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V165b_R2_L_dis = list(
      analyte = "VEGF165b:luminal VEGFR2 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V165b_a2M = list(
      analyte = "VEGF165b:alpha-2-macroglobulin",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_V165b_a2M_fast = list(
      analyte = "VEGF165b:activated alpha-2-macroglobulin",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_V165b_R1_L_dis = list(
      analyte = "alpha-2-macroglobulin:VEGF165b:luminal VEGFR1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_V165b_R2_L_dis = list(
      analyte = "alpha-2-macroglobulin:VEGF165b:luminal VEGFR2 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_fast_V165b_R1_L_dis = list(
      analyte = "activated alpha-2-macroglobulin:VEGF165b:luminal VEGFR1 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_fast_V165b_R2_L_dis = list(
      analyte = "activated alpha-2-macroglobulin:VEGF165b:luminal VEGFR2 (PAD calf endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_V165b_R1_L_norm = list(
      analyte = "alpha-2-macroglobulin:VEGF165b:luminal VEGFR1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_V165b_R2_L_norm = list(
      analyte = "alpha-2-macroglobulin:VEGF165b:luminal VEGFR2 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_fast_V165b_R1_L_norm = list(
      analyte = "activated alpha-2-macroglobulin:VEGF165b:luminal VEGFR1 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    ),
    bld_a2M_fast_V165b_R2_L_norm = list(
      analyte = "activated alpha-2-macroglobulin:VEGF165b:luminal VEGFR2 (normal-tissue endothelium)",
      units = "mol/cm3 blood",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "peripheral arterial disease (ischemic calf muscle)",
    notes = paste(
      "Deterministic systems model, not fitted to individual data. Geometry of",
      "the normal-tissue and blood compartments follows Finley and Popel 2012",
      "(ref. 17); the diseased compartment is the combined lateral and medial",
      "gastrocnemius and soleus muscles (738 cm3, Table 2). Secretion ratios",
      "of VEGF165b to total VEGF165 (0.1 normal, 0.1 blood/endothelium, 1.0",
      "PAD calf) were chosen by the authors so that steady-state VEGF stays",
      "within the 0-10 pM range observed in PAD muscle interstitial fluid",
      "(Hoier 2013) and so that VEGF165b dominates in PAD calf, as in",
      "patient biopsies and mouse hindlimb ischemia western blots.",
      sep = " "
    )
  )

  ini({
    VEGF165Secr_tis_dis <- fixed(0.018); label("Total VEGF165 (a+b) secretion by PAD calf myocytes [molecules/cell/s]") # Suppl. 2 SBML; Results (0.018 = 0.02 total minus 10% VEGF121)
    ECM_Vol_tis_norm <- fixed(0.061987); label("Normal tissue ECM fluid volume fraction [cm3/cm3 tissue]") # Suppl. 2 SBML (Finley 2012, ref 17)
    ECM_Vol_blood <- fixed(0.53876); label("Blood compartment plasma-ECM volume fraction used as plasma fraction [cm3/cm3 blood]") # Suppl. 2 SBML
    ECM_Vol_tis_dis <- fixed(0.3187); label("PAD calf ECM fluid volume fraction [cm3/cm3 tissue]") # Table 2 (31.87%)
    EBM_Vol_tis_norm <- fixed(8.7e-05); label("Normal tissue EBM fluid volume fraction [cm3/cm3 tissue]") # Suppl. 2 SBML
    EBM_Vol_blood <- fixed(0.00043); label("Blood compartment EBM volume fraction [cm3/cm3 blood]") # Suppl. 2 SBML
    EBM_Vol_tis_dis <- fixed(0.0045); label("PAD calf EBM fluid volume fraction [cm3/cm3 tissue]") # Table 2 (0.45%)
    PBM_Vol_blood <- fixed(0.00421); label("Blood compartment PBM volume fraction [cm3/cm3 blood]") # Suppl. 2 SBML
    PBM_Vol_tis_norm <- fixed(0.000307); label("Normal tissue PBM fluid volume fraction [cm3/cm3 tissue]") # Suppl. 2 SBML
    PBM_Vol_tis_dis <- fixed(0.0058); label("PAD calf PBM fluid volume fraction [cm3/cm3 tissue]") # Table 2 (0.58%)
    ECM_conc_tis_norm <- fixed(7.5e-07); label("Normal tissue ECM GAG binding-site concentration [M]") # Suppl. 2 SBML
    EBM_conc_tis_norm <- fixed(1.3e-05); label("Normal tissue EBM GAG binding-site concentration [M]") # Suppl. 2 SBML
    PBM_conc_tis_norm <- fixed(1.3e-05); label("Normal tissue PBM GAG binding-site concentration [M]") # Suppl. 2 SBML
    ECM_conc_tis_dis <- fixed(7.5e-07); label("PAD calf ECM GAG binding-site concentration [M]") # Suppl. 2 SBML
    EBM_conc_tis_dis <- fixed(1.3e-05); label("PAD calf EBM GAG binding-site concentration [M]") # Suppl. 2 SBML
    PBM_conc_tis_dis <- fixed(1.3e-05); label("PAD calf PBM GAG binding-site concentration [M]") # Suppl. 2 SBML
    VesselSA_Vol_tis_norm <- fixed(108); label("Normal tissue abluminal endothelial surface area per tissue volume [cm2/cm3]") # Suppl. 2 SBML
    VesselSA_Vol_blood <- fixed(108); label("Endothelial surface area per tissue volume used in luminal coupling rates [cm2/cm3]") # Suppl. 2 SBML
    VesselSA_Vol_tis_dis <- fixed(417); label("PAD calf receptor-bearing surface area per tissue volume [cm2/cm3]") # Table 2 (417 cm2/cm3 tissue)
    FiberSA_Vol_tis_norm <- fixed(664); label("Normal tissue myocyte surface area per tissue volume [cm2/cm3]") # Suppl. 2 SBML
    VesselCellSurfArea_tis_norm <- fixed(1e-05); label("Surface area of one normal-tissue endothelial cell [cm2]") # Suppl. 2 SBML
    VesselCellSurfArea_blood <- fixed(1e-05); label("Surface area of one luminal endothelial cell [cm2]") # Suppl. 2 SBML
    VesselCellSurfArea_tis_dis <- fixed(1e-05); label("Surface area of one PAD-calf endothelial cell [cm2]") # Suppl. 2 SBML
    FiberCellSurfArea_tis_norm <- fixed(1.85e-05); label("Surface area of one normal-tissue myocyte [cm2]") # Suppl. 2 SBML
    ECM_Kd_tis_norm <- fixed(2.38e-08); label("VEGF165-GAG dissociation constant, normal tissue [M]") # Suppl. 2 SBML; Table 4 (Kd 23.8, unit printed pM)
    ECM_Kd_tis_dis <- fixed(2.38e-08); label("VEGF165-GAG dissociation constant, PAD calf [M]") # Suppl. 2 SBML; Table 4 (Kd 23.8, unit printed pM)
    Avogadro <- fixed(6.0221367e+23); label("Avogadro constant [molecules/mol]") # Suppl. 2 SBML
    Volume_tis_norm <- fixed(61321); label("Normal tissue compartment volume [cm3]") # Suppl. 2 SBML
    Volume_blood <- fixed(5000); label("Blood volume [cm3]") # Suppl. 2 SBML
    Volume_plasma <- fixed(2717); label("Plasma volume [cm3]") # Suppl. 2 SBML
    Volume_tis_dis <- fixed(738); label("PAD calf compartment volume [cm3]") # Table 2 (738 cm3)
    kdissoc_R1_N1_tis_norm_A <- fixed(0.01); label("VEGFR1-NRP1 decoupling, normal, abluminal [1/s]") # Table 4 (1e-2)
    kdissoc_R1_N1_tis_norm_L <- fixed(0.01); label("VEGFR1-NRP1 decoupling, normal, luminal [1/s]") # Table 4 (1e-2)
    R1_recnumber_tis_norm_A <- fixed(550); label("VEGFR1 density per endothelial cell, normal, abluminal [receptors/cell]") # Suppl. 2 SBML (differs from Table 3; see vignette)
    R2_recnumber_tis_norm_A <- fixed(350); label("VEGFR2 density per endothelial cell, normal, abluminal [receptors/cell]") # Suppl. 2 SBML (differs from Table 3; see vignette)
    N1_recnumber_tis_norm_A <- fixed(17500); label("NRP1 density per endothelial cell, normal, abluminal [receptors/cell]") # Suppl. 2 SBML (differs from Table 3; see vignette)
    R1_recnumber_tis_norm_L <- fixed(550); label("VEGFR1 density per endothelial cell, normal, luminal [receptors/cell]") # Suppl. 2 SBML (differs from Table 3; see vignette)
    R2_recnumber_tis_norm_L <- fixed(350); label("VEGFR2 density per endothelial cell, normal, luminal [receptors/cell]") # Suppl. 2 SBML (differs from Table 3; see vignette)
    N1_recnumber_tis_norm_L <- fixed(17500); label("NRP1 density per endothelial cell, normal, luminal [receptors/cell]") # Suppl. 2 SBML (differs from Table 3; see vignette)
    R1_recnumber_tis_dis_A <- fixed(3750); label("VEGFR1 density per endothelial cell, PAD calf, abluminal [receptors/cell]") # Suppl. 2 SBML (differs from Table 3; see vignette)
    R2_recnumber_tis_dis_A <- fixed(300); label("VEGFR2 density per endothelial cell, PAD calf, abluminal [receptors/cell]") # Suppl. 2 SBML (differs from Table 3; see vignette)
    N1_recnumber_tis_dis_A <- fixed(20000); label("NRP1 density per endothelial cell, PAD calf, abluminal [receptors/cell]") # Suppl. 2 SBML (differs from Table 3; see vignette)
    R1_recnumber_tis_dis_L <- fixed(3750); label("VEGFR1 density per endothelial cell, PAD calf, luminal [receptors/cell]") # Suppl. 2 SBML (differs from Table 3; see vignette)
    R2_recnumber_tis_dis_L <- fixed(300); label("VEGFR2 density per endothelial cell, PAD calf, luminal [receptors/cell]") # Suppl. 2 SBML (differs from Table 3; see vignette)
    N1_recnumber_tis_dis_L <- fixed(20000); label("NRP1 density per endothelial cell, PAD calf, luminal [receptors/cell]") # Suppl. 2 SBML (differs from Table 3; see vignette)
    c_V121 <- fixed(0.00107966850554509); label("Plasma clearance of free VEGF121 [1/s]") # Table 5 (VEGF clearance 1.1e-3 1/s)
    c_V165 <- fixed(0.00107966850554509); label("Plasma clearance of free VEGF165 (also VEGF165b) [1/s]") # Table 5 (VEGF clearance 1.1e-3 1/s)
    kperm_blood_and_tis_dis <- fixed(4e-07); label("VEGF microvascular permeability, PAD calf-blood [cm/s]") # Table 5 (4.0e-7 cm/s)
    kperm_blood_and_tis_norm <- fixed(4e-08); label("VEGF microvascular permeability, normal-blood [cm/s]") # Table 5 (4.0e-8 cm/s)
    k_lymph <- fixed(0.033333333); label("Lymphatic flow [cm3/s]") # Suppl. 2 SBML
    kdeg <- fixed(0.000193); label("Interstitial VEGF degradation rate [1/s]") # Suppl. 2 SBML; Table 5 (degradation 1.9e-4 1/s)
    k_int <- fixed(0.00028); label("Abluminal receptor internalization rate [1/s]") # Suppl. 2 SBML
    k_int_L <- fixed(0.00028); label("Luminal receptor internalization rate [1/s]") # Suppl. 2 SBML
    sR <- fixed(0.00028); label("Abluminal receptor insertion rate scale [1/s]") # Suppl. 2 SBML
    sR_L <- fixed(0.00028); label("Luminal receptor insertion rate scale [1/s]") # Suppl. 2 SBML
    kc_R1_N1_A <- fixed(100000000000000); label("VEGFR1-NRP1 coupling, abluminal [1/(mol/cm2)/s]") # Table 4 (1e14)
    kc_R1_N1_L <- fixed(100000000000000); label("VEGFR1-NRP1 coupling, luminal [1/(mol/cm2)/s]") # Table 4 (1e14)
    kon_V165_R1_A <- fixed(30000000); label("V165-VEGFR1 association, abluminal [1/M/s]") # Table 4 (3e7)
    kon_V165_R2_A <- fixed(10000000); label("V165-VEGFR2 association, abluminal [1/M/s]") # Table 4 (1e7)
    kon_V165_N1_A <- fixed(3200000); label("V165-NRP1 association, abluminal [1/M/s]") # Table 4 (3.2e6)
    kon_V121_R1_A <- fixed(30000000); label("V121-VEGFR1 association, abluminal [1/M/s]") # Table 4 (3e7)
    kon_V121_R2_A <- fixed(10000000); label("V121-VEGFR2 association, abluminal [1/M/s]") # Table 4 (1e7)
    kc_V165N1_R2_A <- fixed(100000000000000); label("VEGF165:NRP1-VEGFR2 coupling, abluminal [1/(mol/cm2)/s]") # Table 4 (1e14)
    kc_V165R2_N1_A <- fixed(31000000000000); label("VEGF165:VEGFR2-NRP1 coupling, abluminal [1/(mol/cm2)/s]") # Table 4 (3.1e13)
    kon_V165_R1_L <- fixed(30000000); label("V165-VEGFR1 association, luminal [1/M/s]") # Table 4 (3e7)
    kon_V165_R2_L <- fixed(10000000); label("V165-VEGFR2 association, luminal [1/M/s]") # Table 4 (1e7)
    kon_V165_N1_L <- fixed(3200000); label("V165-NRP1 association, luminal [1/M/s]") # Table 4 (3.2e6)
    kon_V121_R1_L <- fixed(30000000); label("V121-VEGFR1 association, luminal [1/M/s]") # Table 4 (3e7)
    kon_V121_R2_L <- fixed(10000000); label("V121-VEGFR2 association, luminal [1/M/s]") # Table 4 (1e7)
    kc_V165N1_R2_L <- fixed(100000000000000); label("VEGF165:NRP1-VEGFR2 coupling, luminal [1/(mol/cm2)/s]") # Table 4 (1e14)
    kc_V165R2_N1_L <- fixed(31000000000000); label("VEGF165:VEGFR2-NRP1 coupling, luminal [1/(mol/cm2)/s]") # Table 4 (3.1e13)
    koff_V165_MECM <- fixed(0.00069); label("VEGF165 (and VEGF165b)-GAG (ECM) dissociation [1/s]") # Suppl. 2 SBML (Table 4 prints 1e-2)
    koff_V165_MPBM <- fixed(0.00069); label("VEGF165 (and VEGF165b)-GAG (PBM) dissociation [1/s]") # Suppl. 2 SBML (Table 4 prints 1e-2)
    koff_V165_MEBM <- fixed(0.00069); label("VEGF165 (and VEGF165b)-GAG (EBM) dissociation [1/s]") # Suppl. 2 SBML (Table 4 prints 1e-2)
    koff_V165_R1_tis_dis_A <- fixed(0.001); label("V165-VEGFR1 dissociation, PAD calf, abluminal [1/s]") # Table 4 (1e-3)
    koff_V165_R2_tis_dis_A <- fixed(0.001); label("V165-VEGFR2 dissociation, PAD calf, abluminal [1/s]") # Table 4 (1e-3)
    koff_V165_N1_tis_dis_A <- fixed(0.001); label("V165-NRP1 dissociation, PAD calf, abluminal [1/s]") # Table 4 (1e-3)
    koff_V165N1_R2_tis_dis_A <- fixed(0.001); label("VEGF165:VEGFR2:NRP1 decoupling (V165N1_R2), PAD calf, abluminal [1/s]") # Table 4 (1e-3)
    koff_V165R2_N1_tis_dis_A <- fixed(0.001); label("VEGF165:VEGFR2:NRP1 decoupling (V165R2_N1), PAD calf, abluminal [1/s]") # Table 4 (1e-3)
    koff_V121_R1_tis_dis_A <- fixed(0.001); label("V121-VEGFR1 dissociation, PAD calf, abluminal [1/s]") # Table 4 (1e-3)
    koff_V121_R2_tis_dis_A <- fixed(0.001); label("V121-VEGFR2 dissociation, PAD calf, abluminal [1/s]") # Table 4 (1e-3)
    kdissoc_R1_N1_tis_dis_A <- fixed(0.01); label("VEGFR1-NRP1 decoupling, PAD calf, abluminal [1/s]") # Table 4 (1e-2)
    kdissoc_R1_N1_tis_dis_L <- fixed(0.01); label("VEGFR1-NRP1 decoupling, PAD calf, luminal [1/s]") # Table 4 (1e-2)
    koff_V165_R1_tis_dis_L <- fixed(0.001); label("V165-VEGFR1 dissociation, PAD calf, luminal [1/s]") # Table 4 (1e-3)
    koff_V165_R2_tis_dis_L <- fixed(0.001); label("V165-VEGFR2 dissociation, PAD calf, luminal [1/s]") # Table 4 (1e-3)
    koff_V165_N1_tis_dis_L <- fixed(0.001); label("V165-NRP1 dissociation, PAD calf, luminal [1/s]") # Table 4 (1e-3)
    koff_V165R2_N1_tis_dis_L <- fixed(0.001); label("VEGF165:VEGFR2:NRP1 decoupling (V165R2_N1), PAD calf, luminal [1/s]") # Table 4 (1e-3)
    koff_V165N1_R2_tis_dis_L <- fixed(0.001); label("VEGF165:VEGFR2:NRP1 decoupling (V165N1_R2), PAD calf, luminal [1/s]") # Table 4 (1e-3)
    koff_V121_R1_tis_dis_L <- fixed(0.001); label("V121-VEGFR1 dissociation, PAD calf, luminal [1/s]") # Table 4 (1e-3)
    koff_V121_R2_tis_dis_L <- fixed(0.001); label("V121-VEGFR2 dissociation, PAD calf, luminal [1/s]") # Table 4 (1e-3)
    koff_V165_R1_tis_norm_A <- fixed(0.001); label("V165-VEGFR1 dissociation, normal, abluminal [1/s]") # Table 4 (1e-3)
    koff_V165_R2_tis_norm_A <- fixed(0.001); label("V165-VEGFR2 dissociation, normal, abluminal [1/s]") # Table 4 (1e-3)
    koff_V165_N1_tis_norm_A <- fixed(0.001); label("V165-NRP1 dissociation, normal, abluminal [1/s]") # Table 4 (1e-3)
    koff_V165N1_R2_tis_norm_A <- fixed(0.001); label("VEGF165:VEGFR2:NRP1 decoupling (V165N1_R2), normal, abluminal [1/s]") # Table 4 (1e-3)
    koff_V165R2_N1_tis_norm_A <- fixed(0.001); label("VEGF165:VEGFR2:NRP1 decoupling (V165R2_N1), normal, abluminal [1/s]") # Table 4 (1e-3)
    koff_V121_R1_tis_norm_A <- fixed(0.001); label("V121-VEGFR1 dissociation, normal, abluminal [1/s]") # Table 4 (1e-3)
    koff_V121_R2_tis_norm_A <- fixed(0.001); label("V121-VEGFR2 dissociation, normal, abluminal [1/s]") # Table 4 (1e-3)
    koff_V165_R1_tis_norm_L <- fixed(0.001); label("V165-VEGFR1 dissociation, normal, luminal [1/s]") # Table 4 (1e-3)
    koff_V165_R2_tis_norm_L <- fixed(0.001); label("V165-VEGFR2 dissociation, normal, luminal [1/s]") # Table 4 (1e-3)
    koff_V165_N1_tis_norm_L <- fixed(0.001); label("V165-NRP1 dissociation, normal, luminal [1/s]") # Table 4 (1e-3)
    koff_V165N1_R2_tis_norm_L <- fixed(0.001); label("VEGF165:VEGFR2:NRP1 decoupling (V165N1_R2), normal, luminal [1/s]") # Table 4 (1e-3)
    koff_V165R2_N1_tis_norm_L <- fixed(0.001); label("VEGF165:VEGFR2:NRP1 decoupling (V165R2_N1), normal, luminal [1/s]") # Table 4 (1e-3)
    koff_V121_R1_tis_norm_L <- fixed(0.001); label("V121-VEGFR1 dissociation, normal, luminal [1/s]") # Table 4 (1e-3)
    koff_V121_R2_tis_norm_L <- fixed(0.001); label("V121-VEGFR2 dissociation, normal, luminal [1/s]") # Table 4 (1e-3)
    ratioV165 <- fixed(0.111111); label("VEGF121-to-VEGF165 myocyte secretion ratio [-]") # Suppl. 2 SBML; Results (VEGF121 = 10% of total VEGF)
    Kd_VEGF_a2M <- fixed(4e-06); label("VEGF-a2M dissociation constant [M]") # Table 4 (4.0 uM)
    koff_VEGF_a2M <- fixed(0.0001); label("VEGF-a2M dissociation rate [1/s]") # Table 4 (1e-4 1/s)
    c_a2M <- fixed(3.85081766977747e-05); label("alpha-2-macroglobulin and VEGF:a2M clearance [1/s]") # Table 5 (3.9e-5 1/s)
    ksyn_a2M <- fixed(6.27760296527124e-14); label("alpha-2-macroglobulin synthesis [mol/cm3/s]") # Suppl. 2 SBML; Table 5 (3.5e10 molecules/cm3/s)
    Kd_VEGF_a2M_fast <- fixed(4.2e-07); label("VEGF-a2Mfast dissociation constant [M]") # Table 4 (0.42 uM)
    koff_VEGF_a2M_fast <- fixed(0.0001); label("VEGF-a2Mfast dissociation rate [1/s]") # Table 4 (1e-4 1/s)
    ksyn_a2M_fast <- fixed(3.13880148263562e-14); label("Activated alpha-2-macroglobulin synthesis [mol/cm3/s]") # Table 5 (1.9e10 molecules/cm3/s)
    c_a2M_fast <- fixed(0.003850817669777); label("Activated a2M and VEGF:a2Mfast clearance [1/s]") # Table 5 (3.9e-3 1/s)
    koff_sR1_MECM <- fixed(0.01); label("sVEGFR1-GAG (ECM) dissociation [1/s]") # Table 4 (1e-2)
    koff_sR1_MPBM <- fixed(0.01); label("sVEGFR1-GAG (PBM) dissociation [1/s]") # Table 4 (1e-2)
    koff_sR1_MEBM <- fixed(0.01); label("sVEGFR1-GAG (EBM) dissociation [1/s]") # Table 4 (1e-2)
    koff_V165_sR1 <- fixed(0.001); label("VEGF165 (and VEGF165b)-sVEGFR1 dissociation [1/s]") # Table 4 (1e-3)
    koff_sR1_N1_tis_norm <- fixed(0.01); label("sVEGFR1-NRP1 dissociation, normal abluminal [1/s]") # Table 4 (1e-2)
    koff_sR1_N1_L_norm <- fixed(0.01); label("sVEGFR1-NRP1 dissociation, normal luminal [1/s]") # Table 4 (1e-2)
    koff_sR1_N1_L_dis <- fixed(0.01); label("sVEGFR1-NRP1 dissociation, PAD luminal [1/s]") # Table 4 (1e-2)
    koff_sR1_N1_tis_dis <- fixed(0.01); label("sVEGFR1-NRP1 dissociation, PAD calf [1/s]") # Table 4 (1e-2)
    kon_sR1_MECM <- fixed(420000); label("sVEGFR1-GAG (ECM) association [1/M/s]") # Table 4 (4.2e5)
    kon_sR1_MEBM <- fixed(420000); label("sVEGFR1-GAG (EBM) association [1/M/s]") # Table 4 (4.2e5)
    kon_sR1_MPBM <- fixed(420000); label("sVEGFR1-GAG (PBM) association [1/M/s]") # Table 4 (4.2e5)
    kon_V165_sR1 <- fixed(30000000); label("VEGF165 (and VEGF165b)-sVEGFR1 association [1/M/s]") # Table 4 (3e7)
    kon_V121_sR1 <- fixed(30000000); label("VEGF121-sVEGFR1 association [1/M/s]") # Table 4 (3e7)
    kon_sR1_N1_norm <- fixed(5560000); label("sVEGFR1-NRP1 association, normal abluminal [1/M/s]") # Table 4 (5.6e6)
    kon_sR1_N1_luminal <- fixed(5560000); label("sVEGFR1-NRP1 association, luminal [1/M/s]") # Table 4 (5.6e6)
    kon_sR1_N1_dis <- fixed(5560000); label("sVEGFR1-NRP1 association, PAD calf [1/M/s]") # Table 4 (5.6e6)
    kdeg_sR1 <- fixed(0.000193); label("sVEGFR1 and VEGF:sVEGFR1 degradation rate [1/s]") # Table 5 (1.9e-4 1/s)
    c_sR1 <- fixed(5e-06); label("Plasma clearance of sVEGFR1 and VEGF:sVEGFR1 [1/s]") # Table 5 (sVEGFR1 5.0e-6 1/s); Suppl. 2 SBML uses it for the complexes too
    sR1Secr_normEC <- fixed(1e-12); label("sVEGFR1 secretion by normal endothelium [nmol/cm2/s scale]") # Suppl. 2 SBML
    sR1Secr_disEC <- fixed(1e-12); label("sVEGFR1 secretion by PAD endothelium [nmol/cm2/s scale]") # Suppl. 2 SBML
    kperm_blood_and_tis_norm_sR1 <- fixed(1.5e-08); label("sVEGFR1 permeability, normal-blood [cm/s]") # Table 5 (1.5e-8 cm/s)
    kperm_blood_and_tis_dis_sR1 <- fixed(1.5e-07); label("sVEGFR1 permeability, PAD calf-blood [cm/s]") # Suppl. 2 SBML (Table 5 prints 3.0e-7)
    kperm_blood_and_tis_norm_V165_sR1 <- fixed(1.5e-08); label("VEGF165:sVEGFR1 permeability, normal-blood [cm/s]") # Table 5 (1.5e-8 cm/s)
    kperm_blood_and_tis_dis_V165_sR1 <- fixed(1.5e-07); label("VEGF165:sVEGFR1 permeability, PAD calf-blood [cm/s]") # Table 5 (1.5e-7 cm/s)
    kperm_blood_and_tis_norm_V121_sR1 <- fixed(1.5e-08); label("VEGF121:sVEGFR1 permeability, normal-blood [cm/s]") # Table 5 (1.5e-8 cm/s)
    kperm_blood_and_tis_dis_V121_sR1 <- fixed(1.5e-07); label("VEGF121:sVEGFR1 permeability, PAD calf-blood [cm/s]") # Table 5 (1.5e-7 cm/s)
    factor_disEC_secr <- fixed(1); label("Multiplier on PAD-calf endothelial VEGF secretion [-]") # Suppl. 2 SBML
    koff_V121_sR1 <- fixed(0.001); label("VEGF121-sVEGFR1 dissociation [1/s]") # Table 4 (1e-3)
    VEGF165Secr_tis_norm <- fixed(0.018); label("Total VEGF165 (a+b) secretion by normal-tissue myocytes [molecules/cell/s]") # Suppl. 2 SBML; Results (0.018 = 0.02 total minus 10% VEGF121)
    diseaseCellSurfArea_tis_dis <- fixed(3.464e-05); label("PAD calf myocyte area used for secretion scaling [cm2]") # Table 2 (3,464 um2)
    diseasedSA_Vol_tis_dis <- fixed(417); label("PAD calf myocyte surface area per tissue volume [cm2/cm3]") # Table 2 (417 cm2/cm3 tissue)
    mL_to_L <- fixed(1000); label("Conversion factor cm3 per L [-]") # Suppl. 2 SBML
    M_to_nM <- fixed(1000000000); label("Conversion factor nM per M [-]") # Suppl. 2 SBML
    M_to_uM <- fixed(1000000); label("Conversion factor uM per M [-]") # Suppl. 2 SBML
    M_to_pM <- fixed(1000000000000); label("Conversion factor pM per M [-]") # Suppl. 2 SBML
    V165bSecretionRatio_Normal <- fixed(0.1); label("VEGF165b fraction of VEGF165 secreted in normal tissue [-]") # Suppl. 2 SBML; Conclusions (0.1 in normal)
    V165bSecretionRatio_Disease <- fixed(1); label("VEGF165b fraction of VEGF165 secreted in PAD calf [-]") # Suppl. 2 SBML; Conclusions (1 in the disease compartment)
    VEGF165Secr_EC <- fixed(0.031); label("Total VEGF165 (a+b) secretion by endothelial cells [molecules/cell/s]") # Suppl. 2 SBML; Table 1 header and Methods (0.031 molecules/cell/s)
    VEGF121Secr_EC <- fixed(0.003444444); label("VEGF121 secretion by endothelial cells [molecules/cell/s]") # Suppl. 2 SBML (0.031 x 1/9)
    V165bSecretionRatio_Blood <- fixed(0.1); label("VEGF165b fraction of VEGF165 secreted by endothelial cells (luminal and abluminal) [-]") # Suppl. 2 SBML; Conclusions (0.1 in blood)
  })

  model({
    # Rates are per second (SimBiology model); model time is in hours.
    sec_per_hr <- 3600

    # Derived constants and assignment rules from the SBML export
    Interstit_Vol_blood <- (EBM_Vol_blood + PBM_Vol_blood + ECM_Vol_blood)
    Interstit_Vol_tis_dis <- (EBM_Vol_tis_dis + PBM_Vol_tis_dis + ECM_Vol_tis_dis)
    Interstit_Vol_tis_norm <- (EBM_Vol_tis_norm + PBM_Vol_tis_norm + ECM_Vol_tis_norm)
    S_NB <- (VesselSA_Vol_tis_norm * Volume_tis_norm)
    S_NB_dis <- (VesselSA_Vol_tis_dis * Volume_tis_dis)
    VEGF121Secr_tis_dis <- (VEGF165Secr_tis_dis * ratioV165)
    VEGF121Secr_tis_norm <- (VEGF165Secr_tis_norm * ratioV165)
    VesselSA_Vol_blood_tis_dis <- ((VesselSA_Vol_tis_dis * Volume_tis_dis) / Volume_blood)
    VesselSA_Vol_blood_tis_norm <- ((VesselSA_Vol_tis_norm * Volume_tis_norm) / Volume_blood)
    c_V165b <- c_V165
    k_int_N1_tis_dis_A <- k_int
    k_int_N1_tis_dis_L <- k_int_L
    k_int_N1_tis_norm_A <- k_int
    k_int_N1_tis_norm_L <- k_int_L
    k_int_R1N1_tis_dis_A <- k_int
    k_int_R1N1_tis_dis_L <- k_int_L
    k_int_R1N1_tis_norm_A <- k_int
    k_int_R1N1_tis_norm_L <- k_int_L
    k_int_R1_tis_dis_A <- k_int
    k_int_R1_tis_dis_L <- k_int_L
    k_int_R1_tis_norm_A <- k_int
    k_int_R1_tis_norm_L <- k_int_L
    k_int_R2_tis_dis_A <- k_int
    k_int_R2_tis_dis_L <- k_int_L
    k_int_R2_tis_norm_A <- k_int
    k_int_R2_tis_norm_L <- k_int_L
    k_int_V121R1N1_tis_dis_A <- k_int
    k_int_V121R1N1_tis_dis_L <- k_int_L
    k_int_V121R1N1_tis_norm_A <- k_int
    k_int_V121R1N1_tis_norm_L <- k_int_L
    k_int_V121R1_tis_dis_A <- k_int
    k_int_V121R1_tis_dis_L <- k_int_L
    k_int_V121R1_tis_norm_A <- k_int
    k_int_V121R1_tis_norm_L <- k_int_L
    k_int_V121R2_tis_dis_A <- k_int
    k_int_V121R2_tis_dis_L <- k_int_L
    k_int_V121R2_tis_norm_A <- k_int
    k_int_V121R2_tis_norm_L <- k_int_L
    k_int_V121sR1N1_L_dis <- k_int_L
    k_int_V121sR1N1_L_norm <- k_int_L
    k_int_V121sR1N1_tis_dis <- k_int
    k_int_V121sR1N1_tis_dis_A <- k_int
    k_int_V121sR1N1_tis_norm <- k_int
    k_int_V121sR1N1_tis_norm_A <- k_int
    k_int_V165N1_tis_dis_A <- k_int
    k_int_V165N1_tis_dis_L <- k_int_L
    k_int_V165N1_tis_norm_A <- k_int
    k_int_V165N1_tis_norm_L <- k_int_L
    k_int_V165R1_tis_dis_A <- k_int
    k_int_V165R1_tis_dis_L <- k_int_L
    k_int_V165R1_tis_norm_A <- k_int
    k_int_V165R1_tis_norm_L <- k_int_L
    k_int_V165R2N1_tis_dis_A <- k_int
    k_int_V165R2N1_tis_dis_L <- k_int_L
    k_int_V165R2N1_tis_norm_A <- k_int
    k_int_V165R2N1_tis_norm_L <- k_int_L
    k_int_V165R2_tis_dis_A <- k_int
    k_int_V165R2_tis_dis_L <- k_int_L
    k_int_V165R2_tis_norm_A <- k_int
    k_int_V165R2_tis_norm_L <- k_int_L
    k_int_V165bR1_tis_dis_A <- k_int
    k_int_V165bR1_tis_norm_A <- k_int
    k_int_V165bR2_tis_dis_A <- k_int
    k_int_V165bR2_tis_norm_A <- k_int
    k_int_sR1N1_L_dis <- k_int_L
    k_int_sR1N1_L_norm <- k_int_L
    k_int_sR1N1_tis_dis <- k_int
    k_int_sR1N1_tis_dis_A <- k_int
    k_int_sR1N1_tis_norm <- k_int
    k_int_sR1N1_tis_norm_A <- k_int
    k_lymph_dis <- k_lymph
    k_lymph_from_tis_dis <- ((k_lymph_dis / Volume_tis_dis) / Interstit_Vol_tis_dis)
    k_lymph_from_tis_norm <- (k_lymph / (Volume_tis_norm * Interstit_Vol_tis_norm))
    k_lymph_to_blood <- ((k_lymph / Volume_blood) / Interstit_Vol_tis_norm)
    k_lymph_to_blood_dis <- ((k_lymph_dis / Volume_blood) / Interstit_Vol_tis_dis)
    kc_R1_N1_tis_dis_A <- (kc_R1_N1_A / VesselSA_Vol_tis_dis)
    kc_R1_N1_tis_dis_L <- (((kc_R1_N1_L / VesselSA_Vol_blood) * Volume_blood) / Volume_tis_dis)
    kc_R1_N1_tis_norm_A <- (kc_R1_N1_A / VesselSA_Vol_tis_norm)
    kc_R1_N1_tis_norm_L <- (((kc_R1_N1_L / VesselSA_Vol_blood) * Volume_blood) / Volume_tis_norm)
    kc_V121R1_N1_tis_dis_A <- (kc_R1_N1_A / VesselSA_Vol_tis_dis)
    kc_V121R1_N1_tis_norm_A <- (kc_R1_N1_A / VesselSA_Vol_tis_norm)
    kc_V165N1_R2_tis_dis_A <- (kc_V165N1_R2_A / VesselSA_Vol_tis_dis)
    kc_V165N1_R2_tis_dis_L <- (((kc_V165N1_R2_L / VesselSA_Vol_blood) * Volume_blood) / Volume_tis_dis)
    kc_V165N1_R2_tis_norm_A <- (kc_V165N1_R2_A / VesselSA_Vol_tis_norm)
    kc_V165N1_R2_tis_norm_L <- (((kc_V165N1_R2_L / VesselSA_Vol_blood) * Volume_blood) / Volume_tis_norm)
    kc_V165R2_N1_tis_dis_A <- (kc_V165R2_N1_A / VesselSA_Vol_tis_dis)
    kc_V165R2_N1_tis_dis_L <- (((kc_V165R2_N1_L / VesselSA_Vol_blood) * Volume_blood) / Volume_tis_dis)
    kc_V165R2_N1_tis_norm_A <- (kc_V165R2_N1_A / VesselSA_Vol_tis_norm)
    kc_V165R2_N1_tis_norm_L <- (((kc_V165R2_N1_L / VesselSA_Vol_blood) * Volume_blood) / Volume_tis_norm)
    kdissoc_V121_R1_N1_tis_dis_A <- kdissoc_R1_N1_tis_dis_A
    kdissoc_V121_R1_N1_tis_norm_A <- kdissoc_R1_N1_tis_norm_A
    koff_V165_MEBM_tis_dis <- koff_V165_MEBM
    koff_Diseased_V165b_MEBM <- koff_V165_MEBM_tis_dis
    koff_V165_MECM_tis_dis <- koff_V165_MECM
    koff_Diseased_V165b_MECM <- koff_V165_MECM_tis_dis
    koff_V165_MPBM_tis_dis <- koff_V165_MPBM
    koff_Diseased_V165b_MPBM <- koff_V165_MPBM_tis_dis
    koff_Diseased_V165b_R1 <- koff_V165_R1_tis_dis_A
    koff_Diseased_V165b_R2 <- koff_V165_R2_tis_dis_A
    koff_Diseased_V165b_sR1 <- koff_V165_sR1
    koff_V165_MEBM_tis_norm <- koff_V165_MEBM
    koff_Normal_V165b_MEBM <- koff_V165_MEBM_tis_norm
    koff_V165_MECM_tis_norm <- koff_V165_MECM
    koff_Normal_V165b_MECM <- koff_V165_MECM_tis_norm
    koff_V165_MPBM_tis_norm <- koff_V165_MPBM
    koff_Normal_V165b_MPBM <- koff_V165_MPBM_tis_norm
    koff_Normal_V165b_R1 <- koff_V165_R1_tis_norm_A
    koff_Normal_V165b_R2 <- koff_V165_R2_tis_norm_A
    koff_Normal_V165b_sR1 <- koff_V165_sR1
    koff_V165b_R1_tis_dis_L <- koff_V165_R1_tis_dis_L
    koff_V165b_R1_tis_norm_L <- koff_V165_R1_tis_norm_L
    koff_V165b_R2_tis_dis_L <- koff_V165_R2_tis_dis_L
    koff_V165b_R2_tis_norm_L <- koff_V165_R2_tis_norm_L
    koff_V165b_sR1 <- koff_V165_sR1
    kon_V165_MEBM_tis_dis <- (((koff_V165_MEBM / ECM_Kd_tis_dis) * mL_to_L) / Interstit_Vol_tis_dis)
    kon_Diseased_V165b_MEBM <- kon_V165_MEBM_tis_dis
    kon_V165_MECM_tis_dis <- (((koff_V165_MECM / ECM_Kd_tis_dis) * 1000) / Interstit_Vol_tis_dis)
    kon_Diseased_V165b_MECM <- kon_V165_MECM_tis_dis
    kon_V165_MPBM_tis_dis <- (((koff_V165_MPBM / ECM_Kd_tis_dis) * 1000) / Interstit_Vol_tis_dis)
    kon_Diseased_V165b_MPBM <- kon_V165_MPBM_tis_dis
    kon_Diseased_V165b_R1 <- ((kon_V165_R1_A * mL_to_L) / Interstit_Vol_tis_dis)
    kon_Diseased_V165b_R2 <- ((kon_V165_R2_A * mL_to_L) / Interstit_Vol_tis_dis)
    kon_Diseased_V165b_sR1 <- ((kon_V165_sR1 * mL_to_L) / Interstit_Vol_tis_dis)
    kon_V165_MEBM_tis_norm <- (((koff_V165_MEBM / ECM_Kd_tis_norm) * 1000) / Interstit_Vol_tis_norm)
    kon_Normal_V165b_MEBM <- kon_V165_MEBM_tis_norm
    kon_V165_MECM_tis_norm <- (((koff_V165_MECM / ECM_Kd_tis_norm) * 1000) / Interstit_Vol_tis_norm)
    kon_Normal_V165b_MECM <- kon_V165_MECM_tis_norm
    kon_V165_MPBM_tis_norm <- (((koff_V165_MPBM / ECM_Kd_tis_norm) * 1000) / Interstit_Vol_tis_norm)
    kon_Normal_V165b_MPBM <- kon_V165_MPBM_tis_norm
    kon_V165_R1_tis_norm_A <- ((kon_V165_R1_A * mL_to_L) / Interstit_Vol_tis_norm)
    kon_Normal_V165b_R1 <- kon_V165_R1_tis_norm_A
    kon_V165_R2_tis_norm_A <- ((kon_V165_R2_A * mL_to_L) / Interstit_Vol_tis_norm)
    kon_Normal_V165b_R2 <- kon_V165_R2_tis_norm_A
    kon_V165_sR1_tis_norm <- ((kon_V165_sR1 * mL_to_L) / Interstit_Vol_tis_norm)
    kon_Normal_V165b_sR1 <- kon_V165_sR1_tis_norm
    kon_V121_R1_tis_dis_A <- ((kon_V121_R1_A * mL_to_L) / Interstit_Vol_tis_dis)
    kon_V121_R1_tis_dis_L <- ((kon_V121_R1_L * mL_to_L) / Interstit_Vol_blood)
    kon_V121_R1_tis_norm_A <- ((kon_V121_R1_A * mL_to_L) / Interstit_Vol_tis_norm)
    kon_V121_R1_tis_norm_L <- ((kon_V121_R1_L * mL_to_L) / Interstit_Vol_blood)
    kon_V121_R2_tis_dis_A <- ((kon_V121_R2_A * mL_to_L) / Interstit_Vol_tis_dis)
    kon_V121_R2_tis_dis_L <- ((kon_V121_R2_L * mL_to_L) / Interstit_Vol_blood)
    kon_V121_R2_tis_norm_A <- ((kon_V121_R2_A * mL_to_L) / Interstit_Vol_tis_norm)
    kon_V121_R2_tis_norm_L <- ((kon_V121_R2_L * mL_to_L) / Interstit_Vol_blood)
    kon_V121_sR1_blood <- ((kon_V121_sR1 * mL_to_L) / Interstit_Vol_blood)
    kon_V121_sR1_tis_dis <- ((kon_V121_sR1 * mL_to_L) / Interstit_Vol_tis_dis)
    kon_V121_sR1_tis_norm <- ((kon_V121_sR1 * mL_to_L) / Interstit_Vol_tis_norm)
    kon_V165_N1_tis_dis_A <- ((kon_V165_N1_A * mL_to_L) / Interstit_Vol_tis_dis)
    kon_V165_N1_tis_dis_L <- ((kon_V165_N1_L * mL_to_L) / Interstit_Vol_blood)
    kon_V165_N1_tis_norm_A <- ((kon_V165_N1_A * mL_to_L) / Interstit_Vol_tis_norm)
    kon_V165_N1_tis_norm_L <- ((kon_V165_N1_L * mL_to_L) / Interstit_Vol_blood)
    kon_V165_R1_tis_dis_A <- ((kon_V165_R1_A * mL_to_L) / Interstit_Vol_tis_dis)
    kon_V165_R1_tis_dis_L <- ((kon_V165_R1_L * mL_to_L) / Interstit_Vol_blood)
    kon_V165_R1_tis_norm_L <- ((kon_V165_R1_L * mL_to_L) / Interstit_Vol_blood)
    kon_V165_R2_tis_dis_A <- ((kon_V165_R2_A * mL_to_L) / Interstit_Vol_tis_dis)
    kon_V165_R2_tis_dis_L <- ((kon_V165_R2_L * mL_to_L) / Interstit_Vol_blood)
    kon_V165_R2_tis_norm_L <- ((kon_V165_R2_L * mL_to_L) / Interstit_Vol_blood)
    kon_V165_sR1_blood <- ((kon_V165_sR1 * mL_to_L) / Interstit_Vol_blood)
    kon_V165_sR1_tis_dis <- ((kon_V165_sR1 * mL_to_L) / Interstit_Vol_tis_dis)
    kon_V165b_R1_tis_dis_L <- kon_V165_R1_tis_dis_L
    kon_V165b_R1_tis_norm_L <- kon_V165_R1_tis_norm_L
    kon_V165b_R2_tis_dis_L <- kon_V165_R2_tis_dis_L
    kon_V165b_R2_tis_norm_L <- kon_V165_R2_tis_norm_L
    kon_V165b_sR1_blood <- kon_V165_sR1_blood
    kon_VEGF_a2M <- (((koff_VEGF_a2M / Kd_VEGF_a2M) * 1000) / Interstit_Vol_blood)
    kon_VEGF_a2M_fast <- (((koff_VEGF_a2M_fast / Kd_VEGF_a2M_fast) * 1000) / Interstit_Vol_blood)
    kon_sR1_MEBM_tis_dis <- ((kon_sR1_MEBM * mL_to_L) / Interstit_Vol_tis_dis)
    kon_sR1_MEBM_tis_norm <- ((kon_sR1_MEBM * 1000) / Interstit_Vol_tis_norm)
    kon_sR1_MECM_tis_dis <- ((kon_sR1_MECM * 1000) / Interstit_Vol_tis_dis)
    kon_sR1_MECM_tis_norm <- ((kon_sR1_MECM * 1000) / Interstit_Vol_tis_norm)
    kon_sR1_MPBM_tis_dis <- ((kon_sR1_MPBM * 1000) / Interstit_Vol_tis_dis)
    kon_sR1_MPBM_tis_norm <- ((kon_sR1_MPBM * 1000) / Interstit_Vol_tis_norm)
    kon_sR1_N1_L_dis <- ((kon_sR1_N1_luminal * mL_to_L) / Interstit_Vol_blood)
    kon_sR1_N1_L_norm <- ((kon_sR1_N1_luminal * 1000) / Interstit_Vol_blood)
    kon_sR1_N1_tis_dis <- ((kon_sR1_N1_dis * mL_to_L) / Interstit_Vol_tis_dis)
    kon_sR1_N1_tis_norm <- ((kon_sR1_N1_norm * 1000) / Interstit_Vol_tis_norm)
    kpermeability_V121_blood_tis_dis <- ((kperm_blood_and_tis_dis * VesselSA_Vol_tis_dis) / Interstit_Vol_blood)
    kpermeability_V121_blood_tis_dis_prime <- ((kperm_blood_and_tis_dis * S_NB_dis) / Volume_plasma)
    kpermeability_V121_blood_tis_norm <- ((kperm_blood_and_tis_norm * VesselSA_Vol_tis_norm) / Interstit_Vol_blood)
    kpermeability_V121_blood_tis_norm_prime <- ((kperm_blood_and_tis_norm * S_NB) / Volume_plasma)
    kpermeability_V121_tis_dis_blood <- ((kperm_blood_and_tis_dis * VesselSA_Vol_tis_dis) / Interstit_Vol_tis_dis)
    kpermeability_V121_tis_dis_blood_prime <- (((kperm_blood_and_tis_dis * S_NB_dis) / Volume_blood) / Interstit_Vol_tis_dis)
    kpermeability_V121_tis_norm_blood <- ((kperm_blood_and_tis_norm * VesselSA_Vol_tis_norm) / Interstit_Vol_tis_norm)
    kpermeability_V121_tis_norm_blood_prime <- (((kperm_blood_and_tis_norm * S_NB) / Volume_blood) / Interstit_Vol_tis_norm)
    kpermeability_V165_blood_tis_dis <- ((kperm_blood_and_tis_dis * VesselSA_Vol_tis_dis) / Interstit_Vol_blood)
    kpermeability_V165_blood_tis_dis_prime <- ((kperm_blood_and_tis_dis * S_NB_dis) / Volume_plasma)
    kpermeability_V165_blood_tis_norm <- ((kperm_blood_and_tis_norm * VesselSA_Vol_tis_norm) / Interstit_Vol_blood)
    kpermeability_V165_blood_tis_norm_prime <- ((kperm_blood_and_tis_norm * S_NB) / Volume_plasma)
    kpermeability_V165_tis_dis_blood <- ((kperm_blood_and_tis_dis * VesselSA_Vol_tis_dis) / Interstit_Vol_tis_dis)
    kpermeability_V165_tis_dis_blood_prime <- (((kperm_blood_and_tis_dis * S_NB_dis) / Volume_blood) / Interstit_Vol_tis_dis)
    kpermeability_V165_tis_norm_blood <- ((kperm_blood_and_tis_norm * VesselSA_Vol_tis_norm) / Interstit_Vol_tis_norm)
    kpermeability_V165_tis_norm_blood_prime <- (((kperm_blood_and_tis_norm * S_NB) / Volume_blood) / Interstit_Vol_tis_norm)
    kpermeability_V165b_blood_tis_dis <- kpermeability_V165_blood_tis_dis
    kpermeability_V165b_blood_tis_dis_prime <- ((kperm_blood_and_tis_dis * S_NB_dis) / Volume_plasma)
    kpermeability_V165b_blood_tis_norm <- kpermeability_V165_blood_tis_norm
    kpermeability_V165b_blood_tis_norm_prime <- kpermeability_V165_blood_tis_norm_prime
    kpermeability_V165b_tis_dis_blood <- kpermeability_V165_tis_dis_blood
    kpermeability_V165b_tis_dis_blood_prime <- (((kperm_blood_and_tis_dis * S_NB_dis) / Volume_blood) / Interstit_Vol_tis_dis)
    kpermeability_V165b_tis_norm_blood <- kpermeability_V165_tis_norm_blood
    kpermeability_V165b_tis_norm_blood_prime <- kpermeability_V165_tis_norm_blood_prime
    kpermeability_dis_from_blood_V121_sR1 <- ((kperm_blood_and_tis_dis_V121_sR1 * VesselSA_Vol_tis_dis) / Interstit_Vol_blood)
    kpermeability_dis_from_blood_V165_sR1 <- ((kperm_blood_and_tis_dis_V165_sR1 * VesselSA_Vol_tis_dis) / Interstit_Vol_blood)
    kpermeability_dis_from_blood_V165b_sR1 <- kpermeability_dis_from_blood_V165_sR1
    kpermeability_dis_from_blood_sR1 <- ((kperm_blood_and_tis_dis_sR1 * VesselSA_Vol_tis_dis) / Interstit_Vol_blood)
    kpermeability_dis_to_blood_V121_sR1 <- ((kperm_blood_and_tis_dis_V121_sR1 * VesselSA_Vol_tis_dis) / Interstit_Vol_tis_dis)
    kpermeability_dis_to_blood_V165_sR1 <- ((kperm_blood_and_tis_dis_V165_sR1 * VesselSA_Vol_tis_dis) / Interstit_Vol_tis_dis)
    kpermeability_dis_to_blood_V165b_sR1 <- kpermeability_dis_to_blood_V165_sR1
    kpermeability_dis_to_blood_sR1 <- ((kperm_blood_and_tis_dis_sR1 * VesselSA_Vol_tis_dis) / Interstit_Vol_tis_dis)
    kpermeability_norm_from_blood_V121_sR1 <- ((kperm_blood_and_tis_norm_V121_sR1 * VesselSA_Vol_tis_norm) / Interstit_Vol_blood)
    kpermeability_norm_from_blood_V165_sR1 <- ((kperm_blood_and_tis_norm_V165_sR1 * VesselSA_Vol_tis_norm) / Interstit_Vol_blood)
    kpermeability_norm_from_blood_V165b_sR1 <- kpermeability_norm_from_blood_V165_sR1
    kpermeability_norm_from_blood_sR1 <- ((kperm_blood_and_tis_norm_sR1 * VesselSA_Vol_tis_norm) / Interstit_Vol_blood)
    kpermeability_norm_to_blood_V121_sR1 <- ((kperm_blood_and_tis_norm_V121_sR1 * VesselSA_Vol_tis_norm) / Interstit_Vol_tis_norm)
    kpermeability_norm_to_blood_V165_sR1 <- ((kperm_blood_and_tis_norm_V165_sR1 * VesselSA_Vol_tis_norm) / Interstit_Vol_tis_norm)
    kpermeability_norm_to_blood_V165b_sR1 <- kpermeability_norm_to_blood_V165_sR1
    kpermeability_norm_to_blood_sR1 <- ((kperm_blood_and_tis_norm_sR1 * VesselSA_Vol_tis_norm) / Interstit_Vol_tis_norm)
    qV121_abluminal_disEC <- ((((VEGF121Secr_EC * factor_disEC_secr) / Avogadro) / VesselCellSurfArea_tis_dis) * VesselSA_Vol_tis_dis)
    qV121_abluminal_normEC <- (((VEGF121Secr_EC / Avogadro) / VesselCellSurfArea_tis_norm) * VesselSA_Vol_tis_norm)
    qV121_luminal_disEC <- (((VEGF121Secr_EC / Avogadro) / VesselCellSurfArea_blood) * VesselSA_Vol_blood_tis_dis)
    qV121_luminal_normEC <- (((VEGF121Secr_EC / Avogadro) / VesselCellSurfArea_blood) * VesselSA_Vol_blood_tis_norm)
    qV121_tis_dis <- (((VEGF121Secr_tis_dis / Avogadro) / diseaseCellSurfArea_tis_dis) * diseasedSA_Vol_tis_dis)
    qV121_tis_norm <- (((VEGF121Secr_tis_norm / Avogadro) / FiberCellSurfArea_tis_norm) * FiberSA_Vol_tis_norm)
    qV165_abluminal_disEC <- (((((1 - V165bSecretionRatio_Blood) * VEGF165Secr_EC * factor_disEC_secr) / Avogadro) / VesselCellSurfArea_tis_dis) * VesselSA_Vol_tis_dis)
    qV165_abluminal_normEC <- (((((1 - V165bSecretionRatio_Blood) * VEGF165Secr_EC) / Avogadro) / VesselCellSurfArea_tis_norm) * VesselSA_Vol_tis_norm)
    qV165_luminal_disEC <- (((((1 - V165bSecretionRatio_Blood) * VEGF165Secr_EC) / Avogadro) / VesselCellSurfArea_blood) * VesselSA_Vol_blood_tis_dis)
    qV165_luminal_normEC <- (((((1 - V165bSecretionRatio_Blood) * VEGF165Secr_EC) / Avogadro) / VesselCellSurfArea_blood) * VesselSA_Vol_blood_tis_norm)
    qV165_tis_dis <- (((((1 - V165bSecretionRatio_Disease) * VEGF165Secr_tis_dis) / Avogadro) / diseaseCellSurfArea_tis_dis) * diseasedSA_Vol_tis_dis)
    qV165_tis_norm <- (((((1 - V165bSecretionRatio_Normal) * VEGF165Secr_tis_norm) / Avogadro) / FiberCellSurfArea_tis_norm) * FiberSA_Vol_tis_norm)
    qV165b_abluminal_disEC <- ((((V165bSecretionRatio_Blood * VEGF165Secr_EC * factor_disEC_secr) / Avogadro) / VesselCellSurfArea_tis_dis) * VesselSA_Vol_tis_dis)
    qV165b_abluminal_normEC <- ((((V165bSecretionRatio_Blood * VEGF165Secr_EC) / Avogadro) / VesselCellSurfArea_tis_norm) * VesselSA_Vol_tis_norm)
    qV165b_luminal_disEC <- ((((V165bSecretionRatio_Blood * VEGF165Secr_EC * factor_disEC_secr) / Avogadro) / VesselCellSurfArea_blood) * VesselSA_Vol_blood_tis_dis)
    qV165b_luminal_normEC <- ((((V165bSecretionRatio_Blood * VEGF165Secr_EC) / Avogadro) / VesselCellSurfArea_blood) * VesselSA_Vol_blood_tis_norm)
    qV165b_tis_dis <- ((((V165bSecretionRatio_Disease * VEGF165Secr_tis_dis) / Avogadro) / diseaseCellSurfArea_tis_dis) * diseasedSA_Vol_tis_dis)
    qV165b_tis_norm <- ((((V165bSecretionRatio_Normal * VEGF165Secr_tis_norm) / Avogadro) / FiberCellSurfArea_tis_norm) * FiberSA_Vol_tis_norm)
    qsR1_abluminal_disEC <- (sR1Secr_disEC * VesselSA_Vol_tis_dis * 1e-9)
    qsR1_abluminal_normEC <- (sR1Secr_normEC * VesselSA_Vol_tis_norm * 1e-9)
    qsR1_luminal_disEC <- (sR1Secr_disEC * VesselSA_Vol_blood_tis_dis * 1e-9)
    qsR1_luminal_normEC <- (sR1Secr_normEC * VesselSA_Vol_blood_tis_norm * 1e-9)
    sR_N1_tis_dis_A <- ((((sR * N1_recnumber_tis_dis_A) / VesselCellSurfArea_tis_dis) / Avogadro) * VesselSA_Vol_tis_dis)
    sR_N1_tis_dis_L <- ((((sR_L * N1_recnumber_tis_dis_L) / VesselCellSurfArea_blood) / Avogadro) * VesselSA_Vol_blood_tis_dis)
    sR_N1_tis_norm_A <- ((((sR * N1_recnumber_tis_norm_A) / VesselCellSurfArea_tis_norm) / Avogadro) * VesselSA_Vol_tis_norm)
    sR_N1_tis_norm_L <- ((((sR_L * N1_recnumber_tis_norm_L) / VesselCellSurfArea_blood) / Avogadro) * VesselSA_Vol_blood_tis_norm)
    sR_R1_tis_dis_A <- ((((sR * R1_recnumber_tis_dis_A) / VesselCellSurfArea_tis_dis) / Avogadro) * VesselSA_Vol_tis_dis)
    sR_R1_tis_dis_L <- ((((sR_L * R1_recnumber_tis_dis_L) / VesselCellSurfArea_blood) / Avogadro) * VesselSA_Vol_blood_tis_dis)
    sR_R1_tis_norm_A <- ((((sR * R1_recnumber_tis_norm_A) / VesselCellSurfArea_tis_norm) / Avogadro) * VesselSA_Vol_tis_norm)
    sR_R1_tis_norm_L <- ((((sR_L * R1_recnumber_tis_norm_L) / VesselCellSurfArea_blood) / Avogadro) * VesselSA_Vol_blood_tis_norm)
    sR_R2_tis_dis_A <- ((((sR * R2_recnumber_tis_dis_A) / VesselCellSurfArea_tis_dis) / Avogadro) * VesselSA_Vol_tis_dis)
    sR_R2_tis_dis_L <- ((((sR_L * R2_recnumber_tis_dis_L) / VesselCellSurfArea_blood) / Avogadro) * VesselSA_Vol_blood_tis_dis)
    sR_R2_tis_norm_A <- ((((sR * R2_recnumber_tis_norm_A) / VesselCellSurfArea_tis_norm) / Avogadro) * VesselSA_Vol_tis_norm)
    sR_R2_tis_norm_L <- ((((sR_L * R2_recnumber_tis_norm_L) / VesselCellSurfArea_blood) / Avogadro) * VesselSA_Vol_blood_tis_norm)

    # Initial receptor and GAG binding-site amounts (mol/cm3 of compartment)
    bld_R1_L_norm(0) <- (((R1_recnumber_tis_norm_L / VesselCellSurfArea_blood) / Avogadro) * VesselSA_Vol_blood_tis_norm)
    bld_N1_L_norm(0) <- (((N1_recnumber_tis_norm_L / VesselCellSurfArea_blood) / Avogadro) * VesselSA_Vol_blood_tis_norm)
    nrm_R1(0) <- (((R1_recnumber_tis_norm_A / VesselCellSurfArea_tis_norm) / Avogadro) * VesselSA_Vol_tis_norm)
    nrm_N1(0) <- (((N1_recnumber_tis_norm_A / VesselCellSurfArea_tis_norm) / Avogadro) * VesselSA_Vol_tis_norm)
    nrm_R2(0) <- (((R2_recnumber_tis_norm_A / VesselCellSurfArea_tis_norm) / Avogadro) * VesselSA_Vol_tis_norm)
    bld_R2_L_norm(0) <- (((R2_recnumber_tis_norm_L / VesselCellSurfArea_blood) / Avogadro) * VesselSA_Vol_blood_tis_norm)
    dis_MECM(0) <- ((ECM_conc_tis_dis * ECM_Vol_tis_dis) / 1000)
    dis_MEBM(0) <- ((EBM_conc_tis_dis * EBM_Vol_tis_dis) / 1000)
    dis_MPBM(0) <- ((PBM_conc_tis_dis * PBM_Vol_tis_dis) / 1000)
    dis_R1(0) <- (((R1_recnumber_tis_dis_A / VesselCellSurfArea_tis_dis) / Avogadro) * VesselSA_Vol_tis_dis)
    dis_R2(0) <- (((R2_recnumber_tis_dis_A / VesselCellSurfArea_tis_dis) / Avogadro) * VesselSA_Vol_tis_dis)
    dis_N1(0) <- (((N1_recnumber_tis_dis_A / VesselCellSurfArea_tis_dis) / Avogadro) * VesselSA_Vol_tis_dis)
    bld_R1_L_dis(0) <- (((R1_recnumber_tis_dis_L / VesselCellSurfArea_blood) / Avogadro) * VesselSA_Vol_blood_tis_dis)
    bld_N1_L_dis(0) <- (((N1_recnumber_tis_dis_L / VesselCellSurfArea_blood) / Avogadro) * VesselSA_Vol_blood_tis_dis)
    bld_R2_L_dis(0) <- (((R2_recnumber_tis_dis_L / VesselCellSurfArea_blood) / Avogadro) * VesselSA_Vol_blood_tis_dis)
    nrm_MECM(0) <- ((ECM_conc_tis_norm * ECM_Vol_tis_norm) / 1000)
    nrm_MEBM(0) <- ((EBM_conc_tis_norm * EBM_Vol_tis_norm) / 1000)
    nrm_MPBM(0) <- ((PBM_conc_tis_norm * PBM_Vol_tis_norm) / 1000)

    # Reaction fluxes (mol/cm3/s), one per SBML reaction
    J_Blood_Secretion_From_Normal_V165 <- qV165_luminal_normEC
    J_Blood_Secretion_From_Normal_V121 <- qV121_luminal_normEC
    J_Normal_Secretion_V165_1 <- qV165_tis_norm
    J_Normal_Secretion_V121_1 <- qV121_tis_norm
    J_Blood_Internalization_From_Normal_R1_N1 <- (k_int_R1N1_tis_norm_L * bld_R1_L_norm_N1_L_norm)
    J_Blood_Reaction_Normal_04 <- ((kc_R1_N1_tis_norm_L * bld_R1_L_norm * bld_N1_L_norm) - (kdissoc_R1_N1_tis_norm_L * bld_R1_L_norm_N1_L_norm))
    J_Normal_Reaction_06 <- ((kc_R1_N1_tis_norm_A * nrm_R1 * nrm_N1) - (kdissoc_R1_N1_tis_norm_A * nrm_R1_N1))
    J_Normal_Abluminal_R2_1 <- sR_R2_tis_norm_A
    J_Normal_Internalization_R1 <- (k_int_R1_tis_norm_A * nrm_R1)
    J_Normal_Internalization_R2 <- (k_int_R2_tis_norm_A * nrm_R2)
    J_Normal_Internalization_N1 <- (k_int_N1_tis_norm_A * nrm_N1)
    J_Normal_Internalization_R1_N1 <- (k_int_R1N1_tis_norm_A * nrm_R1_N1)
    J_Blood_Secretion_From_Normal_R1 <- sR_R1_tis_norm_L
    J_Blood_Secretion_From_Normal_R2 <- sR_R2_tis_norm_L
    J_Blood_Secretion_From_Normal_N1 <- sR_N1_tis_norm_L
    J_Blood_Internalization_From_Normal_R1 <- (k_int_R1_tis_norm_L * bld_R1_L_norm)
    J_Blood_Internalization_From_Normal_R2 <- (k_int_R2_tis_norm_L * bld_R2_L_norm)
    J_Blood_Internalization_From_Normal_N1 <- (k_int_N1_tis_norm_L * bld_N1_L_norm)
    J_Diseased_Reaction_16 <- ((kon_V165_MECM_tis_dis * dis_V165 * dis_MECM) - (koff_V165_MECM_tis_dis * dis_V165_MECM))
    J_Diseased_Reaction_15 <- ((kon_V165_MEBM_tis_dis * dis_V165 * dis_MEBM) - (koff_V165_MEBM_tis_dis * dis_V165_MEBM))
    J_Diseased_Reaction_17 <- ((kon_V165_MPBM_tis_dis * dis_V165 * dis_MPBM) - (koff_V165_MPBM_tis_dis * dis_V165_MPBM))
    J_Diseased_Reaction_19 <- ((kon_V165_R1_tis_dis_A * dis_V165 * dis_R1) - (koff_V165_R1_tis_dis_A * dis_V165_R1))
    J_Diseased_Reaction_20 <- ((kon_V165_R2_tis_dis_A * dis_V165 * dis_R2) - (koff_V165_R2_tis_dis_A * dis_V165_R2))
    J_Diseased_Reaction_18 <- ((kon_V165_N1_tis_dis_A * dis_V165 * dis_N1) - (koff_V165_N1_tis_dis_A * dis_V165_N1))
    J_Diseased_Reaction_04 <- ((kon_V121_R1_tis_dis_A * dis_R1 * dis_V121) - (koff_V121_R1_tis_dis_A * dis_V121_R1))
    J_Diseased_Reaction_12 <- ((kon_V121_R2_tis_dis_A * dis_R2 * dis_V121) - (koff_V121_R2_tis_dis_A * dis_V121_R2))
    J_Diseased_Lymphatic_Drainage_V165 <- (k_lymph_from_tis_dis * dis_V165)
    J_Diseased_Lymphatic_Drainage_V121 <- (k_lymph_from_tis_dis * dis_V121)
    J_Diseased_Secretion_V165_1 <- qV165_tis_dis
    J_Diseased_Secretion_V165_2 <- qV165_abluminal_disEC
    J_Diseased_Secretion_V121_1 <- qV121_tis_dis
    J_Diseased_Secretion_V121_2 <- qV121_abluminal_disEC
    J_Diseased_Reaction_29 <- ((kc_V165N1_R2_tis_dis_A * dis_R2 * dis_V165_N1) - (koff_V165N1_R2_tis_dis_A * dis_V165_R2_N1))
    J_Diseased_Reaction_02 <- ((kc_V165R2_N1_tis_dis_A * dis_V165_R2 * dis_N1) - (koff_V165R2_N1_tis_dis_A * dis_V165_R2_N1))
    J_Diseased_Internalization_R1_V121 <- (k_int_V121R1_tis_dis_A * dis_V121_R1)
    J_Diseased_Internalization_R2_V121 <- (k_int_V121R2_tis_dis_A * dis_R2_V121)
    J_Diseased_Internalization_V165_R1 <- (k_int_V165R1_tis_dis_A * dis_V165_R1)
    J_Diseased_Internalization_V165_R2 <- (k_int_V165R2_tis_dis_A * dis_V165_R2)
    J_Diseased_Internalization_V165_N1 <- (k_int_V165N1_tis_dis_A * dis_V165_N1)
    J_Diseased_Internalization_V121_R1_N1 <- (k_int_V121R1N1_tis_dis_A * dis_V121_R1_N1)
    J_Diseased_Internalization_R2_V165_N1 <- (k_int_V165R2N1_tis_dis_A * dis_V165_R2_N1)
    J_Diseased_Reaction_03 <- ((kc_R1_N1_tis_dis_A * dis_R1 * dis_N1) - (kdissoc_R1_N1_tis_dis_A * dis_R1_N1))
    J_Diseased_Reaction_01 <- ((kc_R1_N1_tis_dis_A * dis_N1 * dis_V121_R1) - (kdissoc_R1_N1_tis_dis_A * dis_V121_R1_N1))
    J_Diseased_Abluminal_R1_1 <- sR_R1_tis_dis_A
    J_Diseased_Internalization_R1 <- (k_int_R1_tis_dis_A * dis_R1)
    J_Diseased_Internalization_R2 <- (k_int_R2_tis_dis_A * dis_R2)
    J_Diseased_Internalization_N1 <- (k_int_N1_tis_dis_A * dis_N1)
    J_Blood_Reaction_Disease_04 <- ((kc_R1_N1_tis_dis_L * bld_R1_L_dis * bld_N1_L_dis) - (kdissoc_R1_N1_tis_dis_L * bld_R1_L_dis_N1_L_dis))
    J_Blood_Secretion_From_Diseased_R1 <- sR_R1_tis_dis_L
    J_Blood_Secretion_From_Diseased_R2 <- sR_R2_tis_dis_L
    J_Blood_Secretion_From_Diseased_N1 <- sR_N1_tis_dis_L
    J_Blood_Internalization_From_Diseased_R1 <- (k_int_R1_tis_dis_L * bld_R1_L_dis)
    J_Blood_Internalization_From_Diseased_R2 <- (k_int_R2_tis_dis_L * bld_R2_L_dis)
    J_Blood_Internalization_From_Diseased_N1 <- (k_int_N1_tis_dis_L * bld_N1_L_dis)
    J_Blood_Internalization_From_Diseased_R1_N1 <- (k_int_R1N1_tis_dis_L * bld_R1_L_dis_N1_L_dis)
    J_Diseased_Internalization_R1_N1 <- (k_int_R1N1_tis_dis_A * dis_R1_N1)
    J_Blood_Lymphatic_Drainage_From_Normal_V121 <- (k_lymph_to_blood * nrm_V121)
    J_Blood_Permeable_Transfer_To_Normal_V165 <- (kpermeability_V165_blood_tis_norm_prime * bld_V165)
    J_Blood_Permeable_Transfer_To_Normal_V121 <- (kpermeability_V121_blood_tis_norm_prime * bld_V121)
    J_Blood_Reaction_Disease_12 <- ((kon_V165_R1_tis_dis_L * bld_V165 * bld_R1_L_dis) - (koff_V165_R1_tis_dis_L * bld_V165_R1_L_dis))
    J_Blood_Reaction_Disease_13 <- ((kon_V165_R2_tis_dis_L * bld_V165 * bld_R2_L_dis) - (koff_V165_R2_tis_dis_L * bld_V165_R2_L_dis))
    J_Blood_Reaction_Disease_11 <- ((kon_V165_N1_tis_dis_L * bld_V165 * bld_N1_L_dis) - (koff_V165_N1_tis_dis_L * bld_V165_N1_L_dis))
    J_Blood_Reaction_Disease_09 <- ((kon_V121_R1_tis_dis_L * bld_V121 * bld_R1_L_dis_N1_L_dis) - (koff_V121_R1_tis_dis_L * bld_V121_R1_L_dis_N1_L_dis))
    J_Blood_Reaction_Disease_05 <- ((kon_V121_R1_tis_dis_L * bld_R1_L_dis * bld_V121) - (koff_V121_R1_tis_dis_L * bld_R1_L_dis_V121))
    J_Blood_Reaction_Disease_07 <- ((kon_V121_R2_tis_dis_L * bld_R2_L_dis * bld_V121) - (koff_V121_R2_tis_dis_L * bld_R2_L_dis_V121))
    J_Blood_Reaction_Disease_06 <- ((kc_V165N1_R2_tis_dis_L * bld_R2_L_dis * bld_V165_N1_L_dis) - (koff_V165N1_R2_tis_dis_L * bld_R2_L_dis_V165_N1_L_dis))
    J_Blood_Reaction_Disease_03 <- ((kc_V165R2_N1_tis_dis_L * bld_N1_L_dis * bld_V165_R2_L_dis) - (koff_V165R2_N1_tis_dis_L * bld_R2_L_dis_V165_N1_L_dis))
    J_Blood_Internalization_From_Diseased_R1_V121 <- (k_int_V121R1_tis_dis_L * bld_R1_L_dis_V121)
    J_Blood_Internalization_From_Diseased_R2_V121 <- (k_int_V121R2_tis_dis_L * bld_R2_L_dis_V121)
    J_Blood_Internalization_From_Diseased_V165_R1 <- (k_int_V165R1_tis_dis_L * bld_V165_R1_L_dis)
    J_Blood_Internalization_From_Diseased_V165_R2 <- (k_int_V165R2_tis_dis_L * bld_V165_R2_L_dis)
    J_Blood_Internalization_From_Diseased_V165_N1 <- (k_int_V165N1_tis_dis_L * bld_V165_N1_L_dis)
    J_Blood_Internalization_From_Diseased_R2_V165_N1 <- (k_int_V165R2N1_tis_dis_L * bld_R2_L_dis_V165_N1_L_dis)
    J_Blood_Internalization_From_Diseased_V121_R1_N1 <- (k_int_V121R1N1_tis_dis_L * bld_V121_R1_L_dis_N1_L_dis)
    J_Blood_Reaction_Disease_02 <- ((kc_R1_N1_tis_dis_L * bld_N1_L_dis * bld_R1_L_dis_V121) - (kdissoc_R1_N1_tis_dis_L * bld_V121_R1_L_dis_N1_L_dis))
    J_Diseased_Permeable_Transfer_To_Blood_V165 <- (kpermeability_V165_tis_dis_blood * dis_V165)
    J_Diseased_Permeable_Transfer_To_Blood_V121 <- (kpermeability_V121_tis_dis_blood * dis_V121)
    J_Diseased_Permeable_Transfer_From_Blood_V165 <- (kpermeability_V165_blood_tis_dis * bld_V165)
    J_Diseased_Permeable_Transfer_From_Blood_V121 <- (kpermeability_V121_blood_tis_dis * bld_V121)
    J_Blood_Lymphatic_Drainage_From_Diseased_V165 <- (k_lymph_to_blood_dis * dis_V165)
    J_Blood_Lymphatic_Drainage_From_Diseased_V121 <- (k_lymph_to_blood_dis * dis_V121)
    J_Blood_Permeable_Transfer_From_Normal_V165 <- (kpermeability_V165_tis_norm_blood_prime * nrm_V165)
    J_Blood_Permeable_Transfer_From_Normal_V121 <- (kpermeability_V121_tis_norm_blood_prime * nrm_V121)
    J_Blood_Permeable_Transfer_From_Diseased_V165 <- (kpermeability_V165_tis_dis_blood_prime * dis_V165)
    J_Blood_Permeable_Transfer_From_Diseased_V121 <- (kpermeability_V121_tis_dis_blood_prime * dis_V121)
    J_Blood_Permeable_Transfer_To_Diseased_V165 <- (kpermeability_V165_blood_tis_dis_prime * bld_V165)
    J_Blood_Permeable_Transfer_to_Diseased_V121 <- (kpermeability_V121_blood_tis_dis_prime * bld_V121)
    J_Blood_Clearance_V165 <- (c_V165 * bld_V165)
    J_Blood_Clearance_V121 <- (c_V121 * bld_V121)
    J_Normal_Lymphatic_Drainage_V165 <- (k_lymph_from_tis_norm * nrm_V165)
    J_Normal_Lymphatic_Drainage_V121 <- (k_lymph_from_tis_norm * nrm_V121)
    J_Normal_Permeable_Transfer_To_Blood_V165 <- (kpermeability_V165_tis_norm_blood * nrm_V165)
    J_Normal_Permeable_Transfer_To_Blood_V121 <- (kpermeability_V121_tis_norm_blood * nrm_V121)
    J_Normal_Permeable_Transfer_From_Blood_V165 <- (kpermeability_V165_blood_tis_norm * bld_V165)
    J_Normal_Permeable_Transfer_From_Blood_V121 <- (kpermeability_V121_blood_tis_norm * bld_V121)
    J_Blood_Reaction_Normal_12 <- ((kon_V165_R1_tis_norm_L * bld_V165 * bld_R1_L_norm) - (koff_V165_R1_tis_norm_L * bld_V165_R1_L_norm))
    J_Blood_Reaction_Normal_13 <- ((kon_V165_R2_tis_norm_L * bld_V165 * bld_R2_L_norm) - (koff_V165_R2_tis_norm_L * bld_V165_R2_L_norm))
    J_Blood_Reaction_Normal_11 <- ((kon_V165_N1_tis_norm_L * bld_V165 * bld_N1_L_norm) - (koff_V165_N1_tis_norm_L * bld_V165_N1_L_norm))
    J_Blood_Reaction_Normal_09 <- ((kon_V121_R1_tis_norm_L * bld_V121 * bld_R1_L_norm_N1_L_norm) - (koff_V121_R1_tis_norm_L * bld_V121_R1_L_norm_N1_L_norm))
    J_Blood_Reaction_Normal_05 <- ((kon_V121_R1_tis_norm_L * bld_R1_L_norm * bld_V121) - (koff_V121_R1_tis_norm_L * bld_R1_L_norm_V121))
    J_Blood_Reaction_Normal_07 <- ((kon_V121_R2_tis_norm_L * bld_R2_L_norm * bld_V121) - (koff_V121_R2_tis_norm_L * bld_R2_L_norm_V121))
    J_Blood_Internalization_From_Normal_R2_N1 <- (k_int_V165R2N1_tis_norm_L * bld_R2_L_norm_V165_N1_L_norm)
    J_Blood_Reaction_Normal_06 <- ((kc_V165N1_R2_tis_norm_L * bld_R2_L_norm * bld_V165_N1_L_norm) - (koff_V165N1_R2_tis_norm_L * bld_R2_L_norm_V165_N1_L_norm))
    J_Blood_Reaction_Normal_03 <- ((kc_V165R2_N1_tis_norm_L * bld_N1_L_norm * bld_V165_R2_L_norm) - (koff_V165R2_N1_tis_norm_L * bld_R2_L_norm_V165_N1_L_norm))
    J_Blood_Internalization_From_Normal_V165_R1 <- (k_int_V165R1_tis_norm_L * bld_V165_R1_L_norm)
    J_Blood_Internalization_From_Normal_V165_R2 <- (k_int_V165R2_tis_norm_L * bld_V165_R2_L_norm)
    J_Blood_Internalization_From_Normal_V165_N1 <- (k_int_V165N1_tis_norm_L * bld_V165_N1_L_norm)
    J_Blood_Internalization_From_Normal_V121_R1_N1 <- (k_int_V121R1N1_tis_norm_L * bld_V121_R1_L_norm_N1_L_norm)
    J_Blood_Internalization_From_Normal_R1_V121 <- (k_int_V121R1_tis_norm_L * bld_R1_L_norm_V121)
    J_Blood_Internalization_From_Normal_R2_V121 <- (k_int_V121R2_tis_norm_L * bld_R2_L_norm_V121)
    J_Blood_Reaction_Normal_02 <- ((kc_R1_N1_tis_norm_L * bld_N1_L_norm * bld_R1_L_norm_V121) - (kdissoc_R1_N1_tis_norm_L * bld_V121_R1_L_norm_N1_L_norm))
    J_Normal_Reaction_22 <- ((kon_V165_MECM_tis_norm * nrm_V165 * nrm_MECM) - (koff_V165_MECM_tis_norm * nrm_V165_MECM))
    J_Normal_Reaction_21 <- ((kon_V165_MEBM_tis_norm * nrm_V165 * nrm_MEBM) - (koff_V165_MEBM_tis_norm * nrm_V165_MEBM))
    J_Normal_Reaction_23 <- ((kon_V165_MPBM_tis_norm * nrm_V165 * nrm_MPBM) - (koff_V165_MPBM_tis_norm * nrm_V165_MPBM))
    J_Normal_Reaction_26 <- ((kon_V165_R1_tis_norm_A * nrm_V165 * nrm_R1) - (koff_V165_R1_tis_norm_A * nrm_V165_R1))
    J_Normal_Reaction_27 <- ((kon_V165_R2_tis_norm_A * nrm_V165 * nrm_R2) - (koff_V165_R2_tis_norm_A * nrm_V165_R2))
    J_Normal_Reaction_24 <- ((kon_V165_N1_tis_norm_A * nrm_V165 * nrm_N1) - (koff_V165_N1_tis_norm_A * nrm_V165_N1))
    J_Normal_Reaction_07 <- ((kon_V121_R1_tis_norm_A * nrm_R1 * nrm_V121) - (koff_V121_R1_tis_norm_A * nrm_V121_R1))
    J_Normal_Reaction_09 <- ((kon_V121_R2_tis_norm_A * nrm_R2 * nrm_V121) - (koff_V121_R2_tis_norm_A * nrm_V121_R2))
    J_Normal_Reaction_15 <- ((kon_V121_R1_tis_norm_A * nrm_V121 * nrm_R1_N1) - (koff_V121_R1_tis_norm_A * nrm_V121_R1_N1))
    J_Normal_Reaction_04 <- ((kc_R1_N1_tis_norm_A * nrm_N1 * nrm_V121_R1) - (kdissoc_R1_N1_tis_norm_A * nrm_V121_R1_N1))
    J_Normal_Reaction_08 <- ((kc_V165N1_R2_tis_norm_A * nrm_R2 * nrm_V165_N1) - (koff_V165N1_R2_tis_norm_A * nrm_V165_R2_N1))
    J_Normal_Reaction_05 <- ((kc_V165R2_N1_tis_norm_A * nrm_N1 * nrm_V165_R2) - (koff_V165R2_N1_tis_norm_A * nrm_V165_R2_N1))
    J_Normal_Internalization_R1_V121 <- (k_int_V121R1_tis_norm_A * nrm_V121_R1)
    J_Normal_Internalization_R2_V121 <- (k_int_V121R2_tis_norm_A * nrm_V121_R2)
    J_Normal_Internalization_V165_R1 <- (k_int_V165R1_tis_norm_A * nrm_V165_R1)
    J_Normal_Internalization_V165_R2 <- (k_int_V165R2_tis_norm_A * nrm_V165_R2)
    J_Normal_Internalization_V165_N1 <- (k_int_V165N1_tis_norm_A * nrm_V165_N1)
    J_Normal_Internalization_R2_V165_N1 <- (k_int_V165R2N1_tis_norm_A * nrm_V165_R2_N1)
    J_Normal_Internalization_V121_R1_N1 <- (k_int_V121R1N1_tis_norm_A * nrm_V121_R1_N1)
    J_Normal_Degradation_V165 <- (kdeg * nrm_V165)
    J_Diseased_Degradation_V165 <- (kdeg * dis_V165)
    J_Diseased_Degradation_V121 <- (kdeg * dis_V121)
    J_Diseased_Reaction_10 <- ((kon_V121_R1_tis_dis_A * dis_V121 * dis_R1_N2) - (koff_V121_R1_tis_dis_A * dis_V121_R1_N2))
    J_Normal_Secretion_sR1 <- qsR1_abluminal_normEC
    J_Diseased_Secretion_sR1 <- qsR1_abluminal_disEC
    J_Normal_Reaction_11 <- ((kon_sR1_MECM_tis_norm * nrm_sR1 * nrm_MECM) - (koff_sR1_MECM * nrm_sR1_MECM))
    J_Normal_Reaction_10 <- ((kon_sR1_MEBM_tis_norm * nrm_sR1 * nrm_MEBM) - (koff_sR1_MEBM * nrm_sR1_MEBM))
    J_Normal_Reaction_12 <- ((kon_sR1_MPBM_tis_norm * nrm_sR1 * nrm_MPBM) - (koff_sR1_MPBM * nrm_sR1_MPBM))
    J_Normal_Reaction_28 <- ((kon_V165_sR1_tis_norm * nrm_V165 * nrm_sR1) - (koff_V165_sR1 * nrm_V165_sR1))
    J_Normal_Reaction_19 <- ((kon_V121_sR1_tis_norm * nrm_V121 * nrm_sR1) - (koff_V121_sR1 * nrm_V121_sR1))
    J_Normal_Reaction_13 <- ((kon_sR1_N1_tis_norm * nrm_sR1 * nrm_N1) - (koff_sR1_N1_tis_norm * nrm_sR1_N1))
    J_Normal_Reaction_16 <- ((kon_V121_sR1_tis_norm * nrm_V121 * nrm_sR1_N1) - (koff_V121_sR1 * nrm_V121_sR1_N1))
    J_Blood_Reaction_Soluble_02 <- ((kon_V165_sR1_blood * bld_V165 * bld_sR1) - (koff_V165_sR1 * bld_V165_sR1))
    J_Blood_Reaction_Soluble_01 <- ((kon_V121_sR1_blood * bld_V121 * bld_sR1) - (koff_V121_sR1 * bld_V121_sR1))
    J_Blood_Reaction_Normal_08 <- ((kon_sR1_N1_L_norm * bld_sR1 * bld_N1_L_norm) - (koff_sR1_N1_L_norm * bld_sR1_N1_L_norm))
    J_Blood_Reaction_Disease_08 <- ((kon_sR1_N1_L_dis * bld_sR1 * bld_N1_L_dis) - (koff_sR1_N1_L_dis * bld_sR1_N1_L_dis))
    J_Blood_Reaction_Normal_10 <- ((kon_V121_sR1_blood * bld_V121 * bld_sR1_N1_L_norm) - (koff_V121_sR1 * bld_V121_sR1_N1_L_norm))
    J_Blood_Reaction_Disease_10 <- ((kon_V121_sR1_blood * bld_V121 * bld_sR1_N1_L_dis) - (koff_V121_sR1 * bld_V121_sR1_N1_L_dis))
    J_Diseased_Reaction_06 <- ((kon_sR1_MECM_tis_dis * dis_sR1 * dis_MECM) - (koff_sR1_MECM * dis_sR1_MECM))
    J_Diseased_Reaction_05 <- ((kon_sR1_MEBM_tis_dis * dis_sR1 * dis_MEBM) - (koff_sR1_MEBM * dis_sR1_MEBM))
    J_Diseased_Reaction_07 <- ((kon_sR1_MPBM_tis_dis * dis_sR1 * dis_MPBM) - (koff_sR1_MPBM * dis_sR1_MPBM))
    J_Diseased_Reaction_21 <- ((kon_V165_sR1_tis_dis * dis_V165 * dis_sR1) - (koff_V165_sR1 * dis_V165_sR1))
    J_Diseased_Reaction_13 <- ((kon_V121_sR1_tis_dis * dis_V121 * dis_sR1) - (koff_V121_sR1 * dis_V121_sR1))
    J_Diseased_Reaction_08 <- ((kon_sR1_N1_tis_dis * dis_sR1 * dis_N1) - (koff_sR1_N1_tis_dis * dis_sR1_N1))
    J_Diseased_Reaction_11 <- ((kon_V121_sR1_tis_dis * dis_V121 * dis_sR1_N1) - (koff_V121_sR1 * dis_V121_sR1_N1))
    J_Normal_Lymphatic_Drainage_sR1 <- (k_lymph_from_tis_norm * nrm_sR1)
    J_Diseased_Lymphatic_Drainage_sR1 <- (k_lymph_from_tis_dis * dis_sR1)
    J_Blood_Permeable_Transfer_To_Diseased_sR1 <- ((kpermeability_dis_from_blood_sR1 * bld_sR1 * Volume_tis_dis) / Volume_blood)
    J_Normal_Permeable_Transfer_From_Blood_sR1 <- (kpermeability_norm_from_blood_sR1 * bld_sR1)
    J_Diseased_Permeable_Transfer_From_Blood_sR1 <- (kpermeability_dis_from_blood_sR1 * bld_sR1)
    J_Blood_Permeable_Transfer_From_Diseased_sR1 <- ((kpermeability_dis_to_blood_sR1 * dis_sR1 * Volume_tis_dis) / Volume_blood)
    J_Blood_Lymphatic_Drainage_From_Diseased_sR1 <- (k_lymph_to_blood_dis * dis_sR1)
    J_Normal_Lymphatic_Drainage_V165_sR1 <- (k_lymph_from_tis_norm * nrm_V165_sR1)
    J_Diseased_Degradation_V165_sR1 <- (kdeg_sR1 * dis_V165_sR1)
    J_Blood_Permeable_Transfer_To_Diseased_V165_sR1 <- ((kpermeability_dis_from_blood_V165_sR1 * bld_V165_sR1 * Volume_tis_dis) / Volume_blood)
    J_Normal_Permeable_Transfer_From_Blood_V165_sR1 <- (kpermeability_norm_from_blood_V165_sR1 * bld_V165_sR1)
    J_Diseased_Permeable_Transfer_From_Blood_V165_sR1 <- (kpermeability_dis_from_blood_V165_sR1 * bld_V165_sR1)
    J_Blood_Permeable_Transfer_From_Diseased_V165_sR1 <- ((kpermeability_dis_to_blood_V165_sR1 * dis_V165_sR1 * Volume_tis_dis) / Volume_blood)
    J_Blood_Lymphatic_Drainage_From_Normal_V165_sR1 <- (k_lymph_to_blood * nrm_V165_sR1)
    J_Normal_Degradation_V121_sR1 <- (kdeg_sR1 * nrm_V121_sR1)
    J_Diseased_Lymphatic_Drainage_V121_sR1 <- (k_lymph_from_tis_dis * dis_V121_sR1)
    J_Blood_Permeable_Transfer_To_Diseased_V121_sR1 <- ((kpermeability_dis_from_blood_V121_sR1 * bld_V121_sR1 * Volume_tis_dis) / Volume_blood)
    J_Normal_Permeable_Transfer_From_Blood_V121_sR1 <- (kpermeability_norm_from_blood_V121_sR1 * bld_V121_sR1)
    J_Diseased_Permeable_Transfer_From_Blood_V121_sR1 <- (kpermeability_dis_from_blood_V121_sR1 * bld_V121_sR1)
    J_Blood_Permeable_Transfer_From_Diseased_V121_sR1 <- ((kpermeability_dis_to_blood_V121_sR1 * dis_V121_sR1 * Volume_tis_dis) / Volume_blood)
    J_Blood_Lymphatic_Drainage_From_Normal_V121_sR1 <- (k_lymph_to_blood * nrm_V121_sR1)
    J_Blood_Degradation_Clearance_sR1 <- ((c_sR1 + kdeg_sR1) * bld_sR1)
    J_Blood_Degradation_Clearance_V165_sR1 <- ((c_sR1 + kdeg_sR1) * bld_V165_sR1)
    J_Blood_Degradation_Clearance_V121_sR1 <- ((c_sR1 + kdeg_sR1) * bld_V121_sR1)
    J_Normal_Internalization_sR1_N1 <- (k_int_sR1N1_tis_norm * nrm_sR1_N1)
    J_Blood_Internalization_From_Normal_sR1_N1 <- (k_int_sR1N1_L_norm * bld_sR1_N1_L_norm)
    J_Blood_Internalization_From_Diseased_sR1_N1 <- (k_int_sR1N1_L_dis * bld_sR1_N1_L_dis)
    J_Diseased_Internalization_sR1_N1 <- (k_int_sR1N1_tis_dis * dis_sR1_N1)
    J_Normal_Internalization_V121_sR1_N1 <- (k_int_V121sR1N1_tis_norm * nrm_V121_sR1_N1)
    J_Blood_Internalization_From_Normal_V121_sR1_N1 <- (k_int_V121sR1N1_L_norm * bld_V121_sR1_N1_L_norm)
    J_Blood_Internalization_From_Diseased_V121_sR1_N1 <- (k_int_V121sR1N1_L_dis * bld_V121_sR1_N1_L_dis)
    J_Diseased_Internalization_V121_sR1_N1 <- (k_int_V121sR1N1_tis_dis * dis_V121_sR1_N1)
    J_Blood_Reaction_Disease_01 <- ((kon_sR1_N1_L_dis * bld_V121_sR1 * bld_N1_L_dis) - (koff_sR1_N1_L_dis * bld_V121_sR1_N1_L_dis))
    J_Blood_Reaction_Normal_01 <- ((kon_sR1_N1_L_norm * bld_V121_sR1 * bld_N1_L_norm) - (koff_sR1_N1_L_norm * bld_V121_sR1_N1_L_norm))
    J_Diseased_Reaction_30 <- ((kon_sR1_N1_tis_dis * dis_V121_sR1 * dis_N1) - (koff_sR1_N1_tis_dis * dis_V121_sR1_N1))
    J_Normal_Reaction_01 <- ((kon_sR1_N1_tis_norm * nrm_V121_sR1 * nrm_N1) - (koff_sR1_N1_tis_norm * nrm_V121_sR1_N1))
    J_Blood_Reaction_a2M_01 <- ((kon_VEGF_a2M * bld_V165 * bld_a2M) - (koff_VEGF_a2M * bld_V165_a2M))
    J_Blood_Clearance_a2M <- (c_a2M * bld_a2M)
    J_Blood_Clearance_V165_a2M <- (c_a2M * bld_V165_a2M)
    J_Blood_Clearance_V121_a2M <- (c_a2M * bld_V121_a2M)
    J_Blood_Synthesis_a2M <- ksyn_a2M
    J_Blood_Reaction_Normal_a2M_05 <- ((kon_V165_R1_tis_norm_L * bld_V165_a2M * bld_R1_L_norm) - (koff_V165_R1_tis_norm_L * bld_a2M_V165_R1_L_norm))
    J_Blood_Reaction_Normal_a2M_06 <- ((kon_V165_R2_tis_norm_L * bld_V165_a2M * bld_R2_L_norm) - (koff_V165_R2_tis_norm_L * bld_a2M_V165_R2_L_norm))
    J_Blood_Reaction_Normal_a2M_07 <- ((kon_V165_N1_tis_norm_L * bld_V165_a2M * bld_N1_L_norm) - (koff_V165_N1_tis_norm_L * bld_a2M_V165_N1_L_norm))
    J_Blood_Reaction_Normal_a2M_08 <- ((kc_V165N1_R2_tis_norm_L * bld_R2_L_norm * bld_a2M_V165_N1_L_norm) - (koff_V165N1_R2_tis_norm_L * bld_R2_L_norm_a2M_V165_N1_L_norm))
    J_Blood_Reaction_Normal_a2M_09 <- ((kc_V165R2_N1_tis_norm_L * bld_N1_L_norm * bld_a2M_V165_R2_L_norm) - (koff_V165R2_N1_tis_norm_L * bld_R2_L_norm_a2M_V165_N1_L_norm))
    J_Blood_Reaction_Normal_a2M_10 <- ((kon_V121_R1_tis_norm_L * bld_R1_L_norm * bld_V121_a2M) - (koff_V121_R1_tis_norm_L * bld_R1_L_norm_V121_a2M))
    J_Blood_Reaction_Normal_a2M_11 <- ((kon_V121_R2_tis_norm_L * bld_R2_L_norm * bld_V121_a2M) - (koff_V121_R2_tis_norm_L * bld_R2_L_norm_V121_a2M))
    J_Blood_Reaction_Normal_a2M_12 <- ((kon_V121_R1_tis_norm_L * bld_V121_a2M * bld_R1_L_norm_N1_L_norm) - (koff_V121_R1_tis_norm_L * bld_a2M_V121_R1_L_norm_N1_L_norm))
    J_Blood_Reaction_Normal_a2M_13 <- ((kc_R1_N1_tis_norm_L * bld_N1_L_norm * bld_R1_L_norm_V121_a2M) - (kdissoc_R1_N1_tis_norm_L * bld_a2M_V121_R1_L_norm_N1_L_norm))
    J_Blood_Reaction_Diseased_a2M_05 <- ((kon_V165_R1_tis_dis_L * bld_V165_a2M * bld_R1_L_dis) - (koff_V165_R1_tis_dis_L * bld_a2M_V165_R1_L_dis))
    J_Blood_Reaction_Diseased_a2M_06 <- ((kon_V165_R2_tis_dis_L * bld_V165_a2M * bld_R2_L_dis) - (koff_V165_R2_tis_dis_L * bld_a2M_V165_R2_L_dis))
    J_Blood_Reaction_Diseased_a2M_07 <- ((kon_V165_N1_tis_dis_L * bld_V165_a2M * bld_N1_L_dis) - (koff_V165_N1_tis_dis_L * bld_a2M_V165_N1_L_dis))
    J_Blood_Reaction_Diseased_a2M_08 <- ((kc_V165N1_R2_tis_dis_L * bld_R2_L_dis * bld_a2M_V165_N1_L_dis) - (koff_V165N1_R2_tis_dis_L * bld_R2_L_dis_a2M_V165_N1_L_dis))
    J_Blood_Reaction_Diseased_a2M_09 <- ((kc_V165R2_N1_tis_dis_L * bld_N1_L_dis * bld_a2M_V165_R2_L_dis) - (koff_V165R2_N1_tis_dis_L * bld_R2_L_dis_a2M_V165_N1_L_dis))
    J_Blood_Reaction_Diseased_a2M_10 <- ((kon_V121_R1_tis_dis_L * bld_R1_L_dis * bld_V121_a2M) - (koff_V121_R1_tis_dis_L * bld_R1_L_dis_V121_a2M))
    J_Blood_Reaction_Diseased_a2M_11 <- ((kon_V121_R2_tis_dis_L * bld_R2_L_dis * bld_V121_a2M) - (koff_V121_R2_tis_dis_L * bld_R2_L_dis_V121_a2M))
    J_Blood_Reaction_Diseased_a2M_12 <- ((kon_V121_R1_tis_dis_L * bld_V121_a2M * bld_R1_L_dis_N1_L_dis) - (koff_V121_R1_tis_dis_L * bld_a2M_V121_R1_L_dis_N1_L_dis))
    J_Blood_Reaction_Diseased_a2M_13 <- ((kc_R1_N1_tis_dis_L * bld_N1_L_dis * bld_R1_L_dis_V121_a2M) - (kdissoc_R1_N1_tis_dis_L * bld_a2M_V121_R1_L_dis_N1_L_dis))
    J_Blood_Reaction_a2M_fast_01 <- ((kon_VEGF_a2M_fast * bld_V165 * bld_a2M_fast) - (koff_VEGF_a2M_fast * bld_V165_a2M_fast))
    J_Blood_Reaction_a2M_fast_02 <- ((kon_VEGF_a2M_fast * bld_V121 * bld_a2M_fast) - (koff_VEGF_a2M_fast * bld_V121_a2M_fast))
    J_Blood_Synthesis_a2M_fast <- ksyn_a2M_fast
    J_Blood_Clearance_a2M_fast <- (c_a2M_fast * bld_a2M_fast)
    J_Blood_Clearance_V165_a2M_fast <- (c_a2M_fast * bld_V165_a2M_fast)
    J_Blood_Clearance_V121_a2M_fast <- (c_a2M_fast * bld_V121_a2M_fast)
    J_Blood_Reaction_Normal_a2M_fast_03 <- ((kon_V165_R1_tis_norm_L * bld_V165_a2M_fast * bld_R1_L_norm) - (koff_V165_R1_tis_norm_L * bld_a2M_fast_V165_R1_L_norm))
    J_Blood_Reaction_Normal_a2M_fast_04 <- ((kon_V165_R2_tis_norm_L * bld_V165_a2M_fast * bld_R2_L_norm) - (koff_V165_R2_tis_norm_L * bld_a2M_fast_V165_R2_L_norm))
    J_Blood_Reaction_Normal_a2M_fast_05 <- ((kon_V165_N1_tis_norm_L * bld_V165_a2M_fast * bld_N1_L_norm) - (koff_V165_N1_tis_norm_L * bld_a2M_fast_V165_N1_L_norm))
    J_Blood_Reaction_Normal_a2M_fast_06 <- ((kc_V165N1_R2_tis_norm_L * bld_R2_L_norm * bld_a2M_fast_V165_N1_L_norm) - (koff_V165N1_R2_tis_norm_L * bld_R2_L_norm_a2M_fast_V165_N1_L_norm))
    J_Blood_Reaction_Normal_a2M_fast_07 <- ((kc_V165R2_N1_tis_norm_L * bld_N1_L_norm * bld_a2M_fast_V165_R2_L_norm) - (koff_V165R2_N1_tis_norm_L * bld_R2_L_norm_a2M_fast_V165_N1_L_norm))
    J_Blood_Reaction_Normal_a2M_fast_08 <- ((kon_V121_R1_tis_norm_L * bld_R1_L_norm * bld_V121_a2M_fast) - (koff_V121_R1_tis_norm_L * bld_R1_L_norm_V121_a2M_fast))
    J_Blood_Reaction_Normal_a2M_fast_09 <- ((kon_V121_R2_tis_norm_L * bld_R2_L_norm * bld_V121_a2M_fast) - (koff_V121_R2_tis_norm_L * bld_R2_L_norm_V121_a2M_fast))
    J_Blood_Reaction_Normal_a2M_fast_10 <- ((kon_V121_R1_tis_norm_L * bld_V121_a2M_fast * bld_R1_L_norm_N1_L_norm) - (koff_V121_R1_tis_norm_L * bld_a2M_fast_V121_R1_L_norm_N1_L_norm))
    J_Blood_Reaction_Normal_a2M_fast_11 <- ((kc_R1_N1_tis_norm_L * bld_N1_L_norm * bld_R1_L_norm_V121_a2M_fast) - (kdissoc_R1_N1_tis_norm_L * bld_a2M_fast_V121_R1_L_norm_N1_L_norm))
    J_Blood_Reaction_Diseased_a2M_fast_03 <- ((kon_V165_R1_tis_dis_L * bld_V165_a2M_fast * bld_R1_L_dis) - (koff_V165_R1_tis_dis_L * bld_a2M_fast_V165_R1_L_dis))
    J_Blood_Reaction_Diseased_a2M_fast_04 <- ((kon_V165_R2_tis_dis_L * bld_V165_a2M_fast * bld_R2_L_dis) - (koff_V165_R2_tis_dis_L * bld_a2M_fast_V165_R2_L_dis))
    J_Blood_Reaction_Diseased_a2M_fast_05 <- ((kon_V165_N1_tis_dis_L * bld_V165_a2M_fast * bld_N1_L_dis) - (koff_V165_N1_tis_dis_L * bld_a2M_fast_V165_N1_L_dis))
    J_Blood_Reaction_Diseased_a2M_fast_06 <- ((kc_V165N1_R2_tis_dis_L * bld_R2_L_dis * bld_a2M_fast_V165_N1_L_dis) - (koff_V165N1_R2_tis_dis_L * bld_R2_L_dis_a2M_fast_V165_N1_L_dis))
    J_Blood_Reaction_Diseased_a2M_fast_07 <- ((kc_V165R2_N1_tis_dis_L * bld_N1_L_dis * bld_a2M_fast_V165_R2_L_dis) - (koff_V165R2_N1_tis_dis_L * bld_R2_L_dis_a2M_fast_V165_N1_L_dis))
    J_Blood_Reaction_Diseased_a2M_fast_08 <- ((kon_V121_R1_tis_dis_L * bld_R1_L_dis * bld_V121_a2M_fast) - (koff_V121_R1_tis_dis_L * bld_R1_L_dis_V121_a2M_fast))
    J_Blood_Reaction_Diseased_a2M_fast_09 <- ((kon_V121_R2_tis_dis_L * bld_R2_L_dis * bld_V121_a2M_fast) - (koff_V121_R2_tis_dis_L * bld_R2_L_dis_V121_a2M_fast))
    J_Blood_Reaction_Diseased_a2M_fast_10 <- ((kon_V121_R1_tis_dis_L * bld_V121_a2M_fast * bld_R1_L_dis_N1_L_dis) - (koff_V121_R1_tis_dis_L * bld_a2M_fast_V121_R1_L_dis_N1_L_dis))
    J_Blood_Reaction_Diseased_a2M_fast_11 <- ((kc_R1_N1_tis_dis_L * bld_N1_L_dis * bld_R1_L_dis_V121_a2M_fast) - (kdissoc_R1_N1_tis_dis_L * bld_a2M_fast_V121_R1_L_dis_N1_L_dis))
    J_Normal_Reaction_31 <- ((kon_Normal_V165b_MPBM * nrm_V165b * nrm_MPBM) - (koff_Normal_V165b_MPBM * nrm_V165b_MPBM))
    J_Normal_Reaction_30 <- ((kon_Normal_V165b_MECM * nrm_V165b * nrm_MECM) - (koff_Normal_V165b_MECM * nrm_V165b_MECM))
    J_Normal_Reaction_29 <- ((kon_Normal_V165b_MEBM * nrm_V165b * nrm_MEBM) - (koff_Normal_V165b_MEBM * nrm_V165b_MEBM))
    J_Diseased_Reaction_24 <- ((kon_Diseased_V165b_MPBM * dis_V165b * dis_MPBM) - (koff_Diseased_V165b_MPBM * dis_V165b_MPBM))
    J_Diseased_Reaction_23 <- ((kon_Diseased_V165b_MECM * dis_V165b * dis_MECM) - (koff_Diseased_V165b_MECM * dis_V165b_MECM))
    J_Diseased_Reaction_22 <- ((kon_Diseased_V165b_MEBM * dis_V165b * dis_MEBM) - (koff_Diseased_V165b_MEBM * dis_V165b_MEBM))
    J_Blood_Reaction_Normal_14 <- ((kon_V165b_R1_tis_norm_L * bld_V165b * bld_R1_L_norm) - (koff_V165b_R1_tis_norm_L * bld_V165b_R1_L_norm))
    J_Blood_Reaction_Normal_15 <- ((kon_V165b_R2_tis_norm_L * bld_V165b * bld_R2_L_norm) - (koff_V165b_R2_tis_norm_L * bld_V165b_R2_L_norm))
    J_Blood_Reaction_Soluble_03 <- ((kon_V165b_sR1_blood * bld_V165b * bld_sR1) - (koff_V165b_sR1 * bld_V165b_sR1))
    J_Normal_Reaction_32 <- ((kon_Normal_V165b_R1 * nrm_V165b * nrm_R1) - (koff_Normal_V165b_R1 * nrm_V165b_R1))
    J_Normal_Reaction_33 <- ((kon_Normal_V165b_R2 * nrm_V165b * nrm_R2) - (koff_Normal_V165b_R2 * nrm_V165b_R2))
    J_Normal_Reaction_34 <- ((kon_Normal_V165b_sR1 * nrm_V165b * nrm_sR1) - (koff_Normal_V165b_sR1 * nrm_V165b_sR1))
    J_Diseased_Reaction_25 <- ((kon_Diseased_V165b_R1 * dis_V165b * dis_R1) - (koff_Diseased_V165b_R1 * dis_V165b_R1))
    J_Diseased_Reaction_26 <- ((kon_Diseased_V165b_R2 * dis_V165b * dis_R2) - (koff_Diseased_V165b_R2 * dis_V165b_R2))
    J_Diseased_Reaction_27 <- ((kon_Diseased_V165b_sR1 * dis_V165b * dis_sR1) - (koff_Diseased_V165b_sR1 * dis_V165b_sR1))
    J_Normal_Abluminal_R1_1 <- sR_R1_tis_norm_A
    J_Diseased_Abluminal_R2_1 <- sR_R2_tis_dis_A
    J_Diseased_Abluminal_R2_2 <- ((-k_int_R2_tis_dis_A) * dis_R2)
    J_Diseased_Abluminal_R2_3 <- ((-kc_V165N1_R2_tis_dis_A) * dis_V165_N1 * dis_R2)
    J_Diseased_Abluminal_R1_2 <- (kdissoc_R1_N1_tis_dis_A * dis_R1_N1)
    J_Diseased_Abluminal_R1_3 <- ((-kc_R1_N1_tis_dis_A) * dis_R1 * dis_N1)
    J_Diseased_Abluminal_N1_1 <- sR_N1_tis_dis_A
    J_Diseased_Abluminal_N1_2 <- ((-k_int_N1_tis_dis_A) * dis_N1)
    J_Diseased_Abluminal_N1_3 <- ((-kc_V121R1_N1_tis_dis_A) * dis_V121_R1 * dis_N1)
    J_Diseased_Abluminal_N1_4 <- (kdissoc_R1_N1_tis_dis_A * dis_V121_R1_N1)
    J_Diseased_Abluminal_N1_5 <- ((-kc_V165R2_N1_tis_dis_A) * dis_V165_R2 * dis_N1)
    J_Diseased_Abluminal_N1_6 <- ((-kc_R1_N1_tis_dis_A) * dis_N1 * dis_R1)
    J_Diseased_Abluminal_N1_7 <- (kdissoc_R1_N1_tis_dis_A * dis_R1_N1)
    J_Diseased_Abluminal_V165_R1_1 <- ((-k_int_V165R1_tis_dis_A) * dis_V165_R1)
    J_Diseased_Abluminal_V165_R2_1 <- ((-k_int_V165R2_tis_dis_A) * dis_V165_R2)
    J_Diseased_Abluminal_V165_R2_2 <- ((-kc_V165R2_N1_tis_dis_A) * dis_V165_R2 * dis_N1)
    J_Diseased_Abluminal_V165b_R1_1 <- ((-k_int_V165bR1_tis_dis_A) * dis_V165b_R1)
    J_Diseased_Abluminal_V165b_R2_1 <- ((-k_int_V165bR2_tis_dis_A) * dis_V165b_R2)
    J_Diseased_Abluminal_V165_N1_1 <- ((-k_int_V165N1_tis_dis_A) * dis_V165_N1)
    J_Diseased_Abluminal_V165_N1_2 <- ((-kc_V165N1_R2_tis_dis_A) * dis_V165_N1 * dis_R2)
    J_Diseased_Abluminal_V165_R2_N1_1 <- ((-k_int_V165R2N1_tis_dis_A) * dis_V165_R2_N1)
    J_Diseased_Abluminal_V165_R2_N1_2 <- (kc_V165R2_N1_tis_dis_A * dis_V165_R2 * dis_N1)
    J_Diseased_Abluminal_V165_R2_N1_3 <- (kc_V165N1_R2_tis_dis_A * dis_V165_N1 * dis_R2)
    J_Diseased_Abluminal_V121_R1_1 <- ((-k_int_V121R1_tis_dis_A) * dis_V121_R1)
    J_Diseased_Abluminal_V121_R1_2 <- ((-kc_R1_N1_tis_dis_A) * dis_V121_R1 * dis_N1)
    J_Diseased_Abluminal_V121_R1_3 <- (kdissoc_R1_N1_tis_dis_A * dis_V121_R1_N1)
    J_Diseased_Abluminal_V121_R2_1 <- ((-k_int_V121R2_tis_dis_A) * dis_V121_R2)
    J_Diseased_Abluminal_V121_R1_N1_1 <- ((-k_int_V121R1N1_tis_dis_A) * dis_V121_R1_N1)
    J_Diseased_Abluminal_V121_R1_N1_2 <- (kc_V121R1_N1_tis_dis_A * dis_V121_R1 * dis_N1)
    J_Diseased_Abluminal_V121_R1_N1_3 <- ((-kdissoc_V121_R1_N1_tis_dis_A) * dis_V121_R1_N1)
    J_Diseased_Abluminal_R1_N1_1 <- ((-k_int_R1N1_tis_dis_A) * dis_R1_N1)
    J_Diseased_Abluminal_R1_N1_2 <- (kc_R1_N1_tis_dis_A * dis_R1 * dis_N1)
    J_Diseased_Abluminal_R1_N1_3 <- ((-kdissoc_R1_N1_tis_dis_A) * dis_R1_N1)
    J_Diseased_Abluminal_V121_sR1_N1_1 <- ((-k_int_V121sR1N1_tis_dis_A) * dis_V121_sR1_N1)
    J_Diseased_Abluminal_sR1_N1_1 <- ((-k_int_sR1N1_tis_dis_A) * dis_sR1_N1)
    J_Normal_Abluminal_R1_2 <- ((-kc_R1_N1_tis_norm_A) * nrm_R1 * nrm_N1)
    J_Normal_Abluminal_R1_3 <- (kdissoc_R1_N1_tis_norm_A * nrm_R1_N1)
    J_Normal_Abluminal_R2_2 <- ((-k_int_R2_tis_norm_A) * nrm_R2)
    J_Normal_Abluminal_R2_3 <- ((-kc_V165N1_R2_tis_norm_A) * nrm_V165_N1 * nrm_R2)
    J_Normal_Abluminal_N1_1 <- sR_N1_tis_norm_A
    J_Normal_Abluminal_N1_2 <- ((-k_int_N1_tis_norm_A) * nrm_N1)
    J_Normal_Abluminal_N1_3 <- ((-kc_V121R1_N1_tis_norm_A) * nrm_V121_R1 * nrm_N1)
    J_Normal_Abluminal_N1_4 <- (kdissoc_R1_N1_tis_norm_A * nrm_V121_R1_N1)
    J_Normal_Abluminal_N1_5 <- ((-kc_V165R2_N1_tis_norm_A) * nrm_V165_R2 * nrm_N1)
    J_Normal_Abluminal_N1_6 <- ((-kc_R1_N1_tis_norm_A) * nrm_N1 * nrm_R1)
    J_Normal_Abluminal_N1_7 <- (kdissoc_R1_N1_tis_norm_A * nrm_R1_N1)
    J_Normal_Abluminal_V165_R2_1 <- ((-kc_V165R2_N1_tis_norm_A) * nrm_V165_R2 * nrm_N1)
    J_Normal_Abluminal_V165_R2_2 <- ((-k_int_V165R2_tis_norm_A) * nrm_V165_R2)
    J_Normal_Abluminal_V165_R1_1 <- ((-k_int_V165R1_tis_norm_A) * nrm_V165_R1)
    J_Normal_Abluminal_V165b_R1_1 <- ((-k_int_V165bR1_tis_norm_A) * nrm_V165b_R1)
    J_Normal_Abluminal_V165b_R2_1 <- ((-k_int_V165bR2_tis_norm_A) * nrm_V165b_R2)
    J_Normal_Abluminal_V165_N1_1 <- ((-k_int_V165N1_tis_norm_A) * nrm_V165_N1)
    J_Normal_Abluminal_V165_N1_2 <- ((-kc_V165N1_R2_tis_norm_A) * nrm_V165_N1 * nrm_R2)
    J_Normal_Abluminal_V165_R2_N1_1 <- ((-k_int_V165R2N1_tis_norm_A) * nrm_V165_R2_N1)
    J_Normal_Abluminal_V165_R2_N1_2 <- (kc_V165R2_N1_tis_norm_A * nrm_V165_R2 * nrm_N1)
    J_Normal_Abluminal_V165_R2_N1_3 <- (kc_V165N1_R2_tis_norm_A * nrm_V165_N1 * nrm_R2)
    J_Normal_Abluminal_V121_R1_1 <- ((-k_int_V121R1_tis_norm_A) * nrm_V121_R1)
    J_Normal_Abluminal_V121_R1_2 <- ((-kc_R1_N1_tis_norm_A) * nrm_V121_R1 * nrm_N1)
    J_Normal_Abluminal_V121_R1_3 <- (kdissoc_R1_N1_tis_norm_A * nrm_V121_R1_N1)
    J_Normal_Abluminal_V121_R2_1 <- ((-k_int_V121R2_tis_norm_A) * nrm_V121_R2)
    J_Normal_Abluminal_V121_R1_N1_1 <- ((-k_int_V121R1N1_tis_norm_A) * nrm_V121_R1_N1)
    J_Normal_Abluminal_V121_R1_N1_2 <- (kc_V121R1_N1_tis_norm_A * nrm_V121_R1 * nrm_N1)
    J_Normal_Abluminal_V121_R1_N1_3 <- ((-kdissoc_V121_R1_N1_tis_norm_A) * nrm_V121_R1_N1)
    J_Normal_Abluminal_R1_N1_1 <- ((-k_int_R1N1_tis_norm_A) * nrm_R1_N1)
    J_Normal_Abluminal_R1_N1_2 <- (kc_R1_N1_tis_norm_A * nrm_R1 * nrm_N1)
    J_Normal_Abluminal_R1_N1_3 <- ((-kdissoc_R1_N1_tis_norm_A) * nrm_R1_N1)
    J_Normal_Abluminal_V121_sR1_N1_1 <- ((-k_int_V121sR1N1_tis_norm_A) * nrm_V121_sR1_N1)
    J_Normal_Abluminal_sR1_N1_1 <- ((-k_int_sR1N1_tis_norm_A) * nrm_sR1_N1)
    J_Blood_Secretion_From_Diseased_V165 <- qV165_luminal_disEC
    J_Blood_Secretion_From_Diseased_V121 <- qV121_luminal_disEC
    J_Diseased_Degradation_V121_sR1 <- (kdeg_sR1 * dis_V121_sR1)
    J_Diseased_Permeable_Transfer_To_Blood_V121_sR1 <- (kpermeability_dis_to_blood_V121_sR1 * dis_V121_sR1)
    J_Normal_Lymphatic_Drainage_V121_sR1 <- (k_lymph_from_tis_norm * nrm_V121_sR1)
    J_Normal_Permeable_Transfer_to_Blood_V121_sR1 <- (kpermeability_norm_to_blood_V121_sR1 * nrm_V121_sR1)
    J_Diseased_Lymphatic_Drainage_V165_sR1 <- (k_lymph_from_tis_dis * dis_V165_sR1)
    J_Diseased_Permeable_Transfer_To_Blood_V165_sR1 <- (kpermeability_dis_to_blood_V165_sR1 * dis_V165_sR1)
    J_Normal_Permeable_Transfer_to_Blood_V165_sR1 <- (kpermeability_norm_to_blood_V165_sR1 * nrm_V165_sR1)
    J_Normal_Degradation_V165_sR1 <- (kdeg_sR1 * nrm_V165_sR1)
    J_Diseased_Degradation_sR1 <- (kdeg_sR1 * dis_sR1)
    J_Diseased_Permeable_Transfer_To_Blood_sR1 <- (kpermeability_dis_to_blood_sR1 * dis_sR1)
    J_Normal_Degradation_sR1 <- (kdeg_sR1 * nrm_sR1)
    J_Normal_Permeable_Transfer_to_Blood_sR1 <- (kpermeability_norm_to_blood_sR1 * nrm_sR1)
    J_Normal_Secretion_V121_2 <- qV121_abluminal_normEC
    J_Normal_Secretion_V165_2 <- qV165_abluminal_normEC
    J_Blood_Lymphatic_Drainage_From_Diseased_V121_sR1 <- (k_lymph_to_blood_dis * dis_V121_sR1)
    J_Diseased_Degradation_V165b <- (kdeg * dis_V165b)
    J_Normal_Degradation_V165b <- (kdeg * nrm_V165b)
    J_Diseased_Secretion_V165b_1 <- qV165b_tis_dis
    J_Normal_Secretion_V165b_1 <- qV165b_tis_norm
    J_Normal_Lymphatic_Drainage_V165b <- (k_lymph_from_tis_norm * nrm_V165b)
    J_Diseased_Lymphatic_Drainage_V165b <- (k_lymph_from_tis_dis * dis_V165b)
    J_Normal_Permeable_Transfer_From_Blood_V165b <- (kpermeability_V165b_blood_tis_norm * bld_V165b)
    J_Diseased_Permeable_Transfer_From_Blood_V165b <- (kpermeability_V165b_blood_tis_dis * bld_V165b)
    J_Normal_Permeable_Transfer_To_Blood_V165b <- (kpermeability_V165b_tis_norm_blood * nrm_V165b)
    J_Diseased_Permeable_Transfer_To_Blood_V165b <- (kpermeability_V165b_tis_dis_blood * dis_V165b)
    J_Blood_Lymphatic_Drainage_From_Diseased_V165_sR1 <- (k_lymph_to_blood_dis * dis_V165_sR1)
    J_Blood_Lymphatic_Drainage_From_Normal_V165 <- (k_lymph_to_blood * nrm_V165)
    J_Normal_Degradation_V121 <- (kdeg * nrm_V121)
    J_Blood_Reaction_a2M_02 <- ((kon_VEGF_a2M * bld_V121 * bld_a2M) - (koff_VEGF_a2M * bld_V121_a2M))
    J_Blood_Reaction_Disease_14 <- ((kon_V165b_R1_tis_dis_L * bld_V165b * bld_R1_L_dis) - (koff_V165b_R1_tis_dis_L * bld_V165b_R1_L_dis))
    J_Blood_Reaction_Disease_15 <- ((kon_V165b_R2_tis_dis_L * bld_V165b * bld_R2_L_dis) - (koff_V165b_R2_tis_dis_L * bld_V165b_R2_L_dis))
    J_Blood_Secretion_From_Normal_V165b <- qV165b_luminal_normEC
    J_Blood_Clearance_V165b <- (c_V165b * bld_V165b)
    J_Blood_Secretion_From_Diseased_V165b <- qV165b_luminal_disEC
    J_Blood_Permeable_Transfer_To_Normal_V165b <- (kpermeability_V165b_blood_tis_norm_prime * bld_V165b)
    J_Blood_Permeable_Transfer_To_Diseased_V165b <- (kpermeability_V165b_blood_tis_dis_prime * bld_V165b)
    J_Blood_Lymphatic_Drainage_From_Normal_sR1 <- (k_lymph_to_blood * nrm_sR1)
    J_Blood_Secretion_From_Normal_sR1 <- qsR1_luminal_normEC
    J_Blood_Secretion_From_Diseased_sR1 <- qsR1_luminal_disEC
    J_Normal_Secretion_V165b_2 <- qV165b_abluminal_normEC
    J_Blood_Permeable_Transfer_From_Normal_sR1 <- ((kpermeability_norm_to_blood_sR1 * nrm_sR1 * Volume_tis_norm) / Volume_blood)
    J_Blood_Permeable_Transfer_To_Normal_sR1 <- ((kpermeability_norm_from_blood_sR1 * bld_sR1 * Volume_tis_norm) / Volume_blood)
    J_Blood_Permeable_Transfer_To_Normal_V165_sR1 <- ((kpermeability_norm_from_blood_V165_sR1 * bld_V165_sR1 * Volume_tis_norm) / Volume_blood)
    J_Blood_Permeable_Transfer_To_Normal_V121_sR1 <- ((kpermeability_norm_from_blood_V121_sR1 * bld_V121_sR1 * Volume_tis_norm) / Volume_blood)
    J_Blood_Permeable_Transfer_From_Normal_V165_sR1 <- ((kpermeability_norm_to_blood_V165_sR1 * nrm_V165_sR1 * Volume_tis_norm) / Volume_blood)
    J_Blood_Permeable_Transfer_From_Normal_V121_sR1 <- ((kpermeability_norm_to_blood_V121_sR1 * nrm_V121_sR1 * Volume_tis_norm) / Volume_blood)
    J_Blood_Lymphatic_Drainage_From_Normal_V165b <- (k_lymph_to_blood * nrm_V165b)
    J_Blood_Lymphatic_Drainage_From_Diseased_V165b <- (k_lymph_to_blood_dis * dis_V165b)
    J_Blood_Permeable_Transfer_From_Normal_V165b <- (kpermeability_V165b_tis_norm_blood_prime * nrm_V165b)
    J_Blood_Permeable_Transfer_From_Diseased_V165b <- (kpermeability_V165b_tis_dis_blood_prime * dis_V165b)
    J_Diseased_Secretion_V165b_2 <- qV165b_abluminal_disEC
    J_Blood_Degradation_Clearance_V165b_sR1 <- ((c_sR1 + kdeg_sR1) * bld_V165b_sR1)
    J_Blood_Lymphatic_Drainage_From_Diseased_V165b_sR1 <- (k_lymph_to_blood_dis * dis_V165b_sR1)
    J_Blood_Lymphatic_Drainage_From_Normal_V165b_sR1 <- (k_lymph_to_blood * nrm_V165b_sR1)
    J_Blood_Permeable_Transfer_From_Diseased_V165b_sR1 <- ((kpermeability_dis_to_blood_V165b_sR1 * dis_V165b_sR1 * Volume_tis_dis) / Volume_blood)
    J_Blood_Permeable_Transfer_From_Normal_V165b_sR1 <- ((kpermeability_norm_to_blood_V165b_sR1 * nrm_V165b_sR1 * Volume_tis_norm) / Volume_blood)
    J_Blood_Permeable_Transfer_To_Diseased_V165b_sR1 <- ((kpermeability_dis_from_blood_V165b_sR1 * bld_V165b_sR1 * Volume_tis_dis) / Volume_blood)
    J_Blood_Permeable_Transfer_To_Normal_V165b_sR1 <- ((kpermeability_norm_from_blood_V165b_sR1 * bld_V165b_sR1 * Volume_tis_norm) / Volume_blood)
    J_Diseased_Degradation_V165b_sR1 <- (kdeg_sR1 * dis_V165b_sR1)
    J_Diseased_Lymphatic_Drainage_V165b_sR1 <- (k_lymph_from_tis_dis * dis_V165b_sR1)
    J_Diseased_Permeable_Transfer_From_Blood_V165b_sR1 <- (kpermeability_dis_from_blood_V165b_sR1 * bld_V165b_sR1)
    J_Diseased_Permeable_Transfer_To_Blood_V165b_sR1 <- (kpermeability_dis_to_blood_V165b_sR1 * dis_V165b_sR1)
    J_Normal_Degradation_V165b_sR1 <- (kdeg_sR1 * nrm_V165b_sR1)
    J_Normal_Lymphatic_Drainage_V165b_sR1 <- (k_lymph_from_tis_norm * nrm_V165b_sR1)
    J_Normal_Permeable_Transfer_From_Blood_V165b_sR1 <- (kpermeability_norm_from_blood_V165b_sR1 * bld_V165b_sR1)
    J_Normal_Permeable_Transfer_to_Blood_V165b_sR1 <- (kpermeability_norm_to_blood_V165b_sR1 * nrm_V165b_sR1)
    J_Blood_Reaction_a2M_03 <- ((kon_VEGF_a2M * bld_V165b * bld_a2M) - (koff_VEGF_a2M * bld_V165b_a2M))
    J_Blood_Reaction_a2M_fast_03 <- ((kon_VEGF_a2M_fast * bld_V165b * bld_a2M_fast) - (koff_VEGF_a2M_fast * bld_V165b_a2M_fast))
    J_Blood_Reaction_Diseased_a2M_14 <- ((kon_V165b_R1_tis_dis_L * bld_V165b_a2M * bld_R1_L_dis) - (koff_V165b_R1_tis_dis_L * bld_a2M_V165b_R1_L_dis))
    J_Blood_Reaction_Diseased_a2M_15 <- ((kon_V165b_R2_tis_dis_L * bld_V165b_a2M * bld_R2_L_dis) - (koff_V165b_R2_tis_dis_L * bld_a2M_V165b_R2_L_dis))
    J_Blood_Reaction_Diseased_a2M_fast_12 <- ((kon_V165b_R1_tis_dis_L * bld_V165b_a2M_fast * bld_R1_L_dis) - (koff_V165b_R1_tis_dis_L * bld_a2M_fast_V165b_R1_L_dis))
    J_Blood_Reaction_Diseased_a2M_fast_13 <- ((kon_V165b_R2_tis_dis_L * bld_V165b_a2M_fast * bld_R2_L_dis) - (koff_V165b_R2_tis_dis_L * bld_a2M_fast_V165b_R2_L_dis))
    J_Blood_Reaction_Normal_a2M_14 <- ((kon_V165b_R1_tis_norm_L * bld_V165b_a2M * bld_R1_L_norm) - (koff_V165b_R1_tis_norm_L * bld_a2M_V165b_R1_L_norm))
    J_Blood_Reaction_Normal_a2M_15 <- ((kon_V165b_R2_tis_norm_L * bld_V165b_a2M * bld_R2_L_norm) - (koff_V165b_R2_tis_norm_L * bld_a2M_V165b_R2_L_norm))
    J_Blood_Reaction_Normal_a2M_fast_12 <- ((kon_V165b_R1_tis_norm_L * bld_V165b_a2M_fast * bld_R1_L_norm) - (koff_V165b_R1_tis_norm_L * bld_a2M_fast_V165b_R1_L_norm))
    J_Blood_Reaction_Normal_a2M_fast_13 <- ((kon_V165b_R2_tis_norm_L * bld_V165b_a2M_fast * bld_R2_L_norm) - (koff_V165b_R2_tis_norm_L * bld_a2M_fast_V165b_R2_L_norm))
    J_Diseased_Internalization_V165b_R1 <- (k_int_V165bR1_tis_dis_A * dis_V165b_R1)
    J_Diseased_Internalization_V165b_R2 <- (k_int_V165bR2_tis_dis_A * dis_V165b_R2)
    J_Normal_Internalization_V165b_R1 <- (k_int_V165bR1_tis_norm_A * nrm_V165b_R1)
    J_Normal_Internalization_V165b_R2 <- (k_int_V165bR2_tis_norm_A * nrm_V165b_R2)
    J_Blood_Clearance_V165b_a2M_fast <- (c_a2M * bld_V165b_a2M_fast)
    J_Blood_Clearance_V165b_a2M <- (c_a2M * bld_V165b_a2M)
    J_Blood_Internalization_From_Diseased_V165b_R1 <- (k_int_V165R1_tis_dis_L * bld_V165b_R1_L_dis)
    J_Blood_Internalization_From_Diseased_V165b_R2 <- (k_int_V165R2_tis_dis_L * bld_V165b_R2_L_dis)
    J_Blood_Internalization_From_Normal_V165b_R1 <- (k_int_V165R1_tis_norm_L * bld_V165b_R1_L_norm)
    J_Blood_Internalization_From_Normal_V165b_R2 <- (k_int_V165R2_tis_norm_L * bld_V165b_R2_L_norm)

    # Mass balances
    d/dt(bld_V165) <- sec_per_hr * (J_Blood_Secretion_From_Normal_V165 - J_Blood_Permeable_Transfer_To_Normal_V165 - J_Blood_Reaction_Disease_12 - J_Blood_Reaction_Disease_13 - J_Blood_Reaction_Disease_11 + J_Blood_Lymphatic_Drainage_From_Diseased_V165 + J_Blood_Permeable_Transfer_From_Normal_V165 + J_Blood_Permeable_Transfer_From_Diseased_V165 - J_Blood_Permeable_Transfer_To_Diseased_V165 - J_Blood_Clearance_V165 - J_Blood_Reaction_Normal_12 - J_Blood_Reaction_Normal_13 - J_Blood_Reaction_Normal_11 - J_Blood_Reaction_Soluble_02 - J_Blood_Reaction_a2M_01 - J_Blood_Reaction_a2M_fast_01 + J_Blood_Secretion_From_Diseased_V165 + J_Blood_Lymphatic_Drainage_From_Normal_V165)
    d/dt(bld_V121) <- sec_per_hr * (J_Blood_Secretion_From_Normal_V121 + J_Blood_Lymphatic_Drainage_From_Normal_V121 - J_Blood_Permeable_Transfer_To_Normal_V121 - J_Blood_Reaction_Disease_09 - J_Blood_Reaction_Disease_05 - J_Blood_Reaction_Disease_07 + J_Blood_Lymphatic_Drainage_From_Diseased_V121 + J_Blood_Permeable_Transfer_From_Normal_V121 + J_Blood_Permeable_Transfer_From_Diseased_V121 - J_Blood_Permeable_Transfer_to_Diseased_V121 - J_Blood_Clearance_V121 - J_Blood_Reaction_Normal_09 - J_Blood_Reaction_Normal_05 - J_Blood_Reaction_Normal_07 - J_Blood_Reaction_Soluble_01 - J_Blood_Reaction_Normal_10 - J_Blood_Reaction_Disease_10 - J_Blood_Reaction_a2M_fast_02 + J_Blood_Secretion_From_Diseased_V121 - J_Blood_Reaction_a2M_02)
    d/dt(nrm_V165) <- sec_per_hr * (J_Normal_Secretion_V165_1 - J_Normal_Lymphatic_Drainage_V165 - J_Normal_Permeable_Transfer_To_Blood_V165 + J_Normal_Permeable_Transfer_From_Blood_V165 - J_Normal_Reaction_22 - J_Normal_Reaction_21 - J_Normal_Reaction_23 - J_Normal_Reaction_26 - J_Normal_Reaction_27 - J_Normal_Reaction_24 - J_Normal_Degradation_V165 - J_Normal_Reaction_28 + J_Normal_Secretion_V165_2)
    d/dt(nrm_V121) <- sec_per_hr * (J_Normal_Secretion_V121_1 - J_Normal_Lymphatic_Drainage_V121 - J_Normal_Permeable_Transfer_To_Blood_V121 + J_Normal_Permeable_Transfer_From_Blood_V121 - J_Normal_Reaction_07 - J_Normal_Reaction_09 - J_Normal_Reaction_15 - J_Normal_Reaction_19 - J_Normal_Reaction_16 + J_Normal_Secretion_V121_2 - J_Normal_Degradation_V121)
    d/dt(bld_R1_L_norm_N1_L_norm) <- sec_per_hr * (-J_Blood_Internalization_From_Normal_R1_N1 + J_Blood_Reaction_Normal_04 - J_Blood_Reaction_Normal_09 - J_Blood_Reaction_Normal_a2M_12 - J_Blood_Reaction_Normal_a2M_fast_10)
    d/dt(bld_R1_L_norm) <- sec_per_hr * (-J_Blood_Reaction_Normal_04 + J_Blood_Secretion_From_Normal_R1 - J_Blood_Internalization_From_Normal_R1 - J_Blood_Reaction_Normal_12 - J_Blood_Reaction_Normal_05 - J_Blood_Reaction_Normal_a2M_05 - J_Blood_Reaction_Normal_a2M_10 - J_Blood_Reaction_Normal_a2M_fast_03 - J_Blood_Reaction_Normal_a2M_fast_08 - J_Blood_Reaction_Normal_14 - J_Blood_Reaction_Normal_a2M_14 - J_Blood_Reaction_Normal_a2M_fast_12)
    d/dt(bld_N1_L_norm) <- sec_per_hr * (-J_Blood_Reaction_Normal_04 + J_Blood_Secretion_From_Normal_N1 - J_Blood_Internalization_From_Normal_N1 - J_Blood_Reaction_Normal_11 - J_Blood_Reaction_Normal_03 - J_Blood_Reaction_Normal_02 - J_Blood_Reaction_Normal_08 - J_Blood_Reaction_Normal_01 - J_Blood_Reaction_Normal_a2M_07 - J_Blood_Reaction_Normal_a2M_09 - J_Blood_Reaction_Normal_a2M_13 - J_Blood_Reaction_Normal_a2M_fast_05 - J_Blood_Reaction_Normal_a2M_fast_07 - J_Blood_Reaction_Normal_a2M_fast_11)
    d/dt(nrm_R1) <- sec_per_hr * (-J_Normal_Reaction_06 - J_Normal_Internalization_R1 - J_Normal_Reaction_26 - J_Normal_Reaction_07 - J_Normal_Reaction_32 + J_Normal_Abluminal_R1_1 + J_Normal_Abluminal_R1_2 + J_Normal_Abluminal_R1_3)
    d/dt(nrm_N1) <- sec_per_hr * (-J_Normal_Reaction_06 - J_Normal_Internalization_N1 - J_Normal_Reaction_24 - J_Normal_Reaction_04 - J_Normal_Reaction_05 - J_Normal_Reaction_13 - J_Normal_Reaction_01 + J_Normal_Abluminal_N1_1 + J_Normal_Abluminal_N1_2 + J_Normal_Abluminal_N1_3 + J_Normal_Abluminal_N1_4 + J_Normal_Abluminal_N1_5 + J_Normal_Abluminal_N1_6 + J_Normal_Abluminal_N1_7)
    d/dt(nrm_R1_N1) <- sec_per_hr * (J_Normal_Reaction_06 - J_Normal_Internalization_R1_N1 - J_Normal_Reaction_15 + J_Normal_Abluminal_R1_N1_1 + J_Normal_Abluminal_R1_N1_2 + J_Normal_Abluminal_R1_N1_3)
    d/dt(nrm_R2) <- sec_per_hr * (J_Normal_Abluminal_R2_1 - J_Normal_Internalization_R2 - J_Normal_Reaction_27 - J_Normal_Reaction_09 - J_Normal_Reaction_08 - J_Normal_Reaction_33 + J_Normal_Abluminal_R2_2 + J_Normal_Abluminal_R2_3)
    d/dt(bld_R2_L_norm) <- sec_per_hr * (J_Blood_Secretion_From_Normal_R2 - J_Blood_Internalization_From_Normal_R2 - J_Blood_Reaction_Normal_13 - J_Blood_Reaction_Normal_07 - J_Blood_Reaction_Normal_06 - J_Blood_Reaction_Normal_a2M_06 - J_Blood_Reaction_Normal_a2M_08 - J_Blood_Reaction_Normal_a2M_11 - J_Blood_Reaction_Normal_a2M_fast_04 - J_Blood_Reaction_Normal_a2M_fast_06 - J_Blood_Reaction_Normal_a2M_fast_09 - J_Blood_Reaction_Normal_15 - J_Blood_Reaction_Normal_a2M_15 - J_Blood_Reaction_Normal_a2M_fast_13)
    d/dt(dis_V165) <- sec_per_hr * (-J_Diseased_Reaction_16 - J_Diseased_Reaction_15 - J_Diseased_Reaction_17 - J_Diseased_Reaction_19 - J_Diseased_Reaction_20 - J_Diseased_Reaction_18 - J_Diseased_Lymphatic_Drainage_V165 + J_Diseased_Secretion_V165_1 + J_Diseased_Secretion_V165_2 - J_Diseased_Permeable_Transfer_To_Blood_V165 + J_Diseased_Permeable_Transfer_From_Blood_V165 - J_Diseased_Degradation_V165 - J_Diseased_Reaction_21)
    d/dt(dis_MECM) <- sec_per_hr * (-J_Diseased_Reaction_16 - J_Diseased_Reaction_06 - J_Diseased_Reaction_23)
    d/dt(dis_V165_MECM) <- sec_per_hr * (J_Diseased_Reaction_16)
    d/dt(dis_MEBM) <- sec_per_hr * (-J_Diseased_Reaction_15 - J_Diseased_Reaction_05 - J_Diseased_Reaction_22)
    d/dt(dis_V165_MEBM) <- sec_per_hr * (J_Diseased_Reaction_15)
    d/dt(dis_MPBM) <- sec_per_hr * (-J_Diseased_Reaction_17 - J_Diseased_Reaction_07 - J_Diseased_Reaction_24)
    d/dt(dis_V165_MPBM) <- sec_per_hr * (J_Diseased_Reaction_17)
    d/dt(dis_R1) <- sec_per_hr * (-J_Diseased_Reaction_19 - J_Diseased_Reaction_04 - J_Diseased_Reaction_03 + J_Diseased_Abluminal_R1_1 - J_Diseased_Internalization_R1 - J_Diseased_Reaction_25 + J_Diseased_Abluminal_R1_2 + J_Diseased_Abluminal_R1_3)
    d/dt(dis_V165_R1) <- sec_per_hr * (J_Diseased_Reaction_19 - J_Diseased_Internalization_V165_R1 + J_Diseased_Abluminal_V165_R1_1)
    d/dt(dis_R2) <- sec_per_hr * (-J_Diseased_Reaction_20 - J_Diseased_Reaction_12 - J_Diseased_Reaction_29 - J_Diseased_Internalization_R2 - J_Diseased_Reaction_26 + J_Diseased_Abluminal_R2_1 + J_Diseased_Abluminal_R2_2 + J_Diseased_Abluminal_R2_3)
    d/dt(dis_V165_R2) <- sec_per_hr * (J_Diseased_Reaction_20 - J_Diseased_Reaction_02 - J_Diseased_Internalization_V165_R2 + J_Diseased_Abluminal_V165_R2_1 + J_Diseased_Abluminal_V165_R2_2)
    d/dt(dis_N1) <- sec_per_hr * (-J_Diseased_Reaction_18 - J_Diseased_Reaction_02 - J_Diseased_Reaction_03 - J_Diseased_Reaction_01 - J_Diseased_Internalization_N1 - J_Diseased_Reaction_08 - J_Diseased_Reaction_30 + J_Diseased_Abluminal_N1_1 + J_Diseased_Abluminal_N1_2 + J_Diseased_Abluminal_N1_3 + J_Diseased_Abluminal_N1_4 + J_Diseased_Abluminal_N1_5 + J_Diseased_Abluminal_N1_6 + J_Diseased_Abluminal_N1_7)
    d/dt(dis_V165_N1) <- sec_per_hr * (J_Diseased_Reaction_18 - J_Diseased_Reaction_29 - J_Diseased_Internalization_V165_N1 + J_Diseased_Abluminal_V165_N1_1 + J_Diseased_Abluminal_V165_N1_2)
    d/dt(dis_V121) <- sec_per_hr * (-J_Diseased_Reaction_04 - J_Diseased_Reaction_12 - J_Diseased_Lymphatic_Drainage_V121 + J_Diseased_Secretion_V121_1 + J_Diseased_Secretion_V121_2 - J_Diseased_Permeable_Transfer_To_Blood_V121 + J_Diseased_Permeable_Transfer_From_Blood_V121 - J_Diseased_Degradation_V121 - J_Diseased_Reaction_10 - J_Diseased_Reaction_13 - J_Diseased_Reaction_11)
    d/dt(dis_V121_R1) <- sec_per_hr * (J_Diseased_Reaction_04 - J_Diseased_Internalization_R1_V121 - J_Diseased_Reaction_01 + J_Diseased_Abluminal_V121_R1_1 + J_Diseased_Abluminal_V121_R1_2 + J_Diseased_Abluminal_V121_R1_3)
    d/dt(dis_V121_R2) <- sec_per_hr * (J_Diseased_Reaction_12 + J_Diseased_Abluminal_V121_R2_1)
    d/dt(dis_V165_R2_N1) <- sec_per_hr * (J_Diseased_Reaction_29 + J_Diseased_Reaction_02 - J_Diseased_Internalization_R2_V165_N1 + J_Diseased_Abluminal_V165_R2_N1_1 + J_Diseased_Abluminal_V165_R2_N1_2 + J_Diseased_Abluminal_V165_R2_N1_3)
    d/dt(dis_R2_V121) <- sec_per_hr * (-J_Diseased_Internalization_R2_V121)
    d/dt(dis_V121_R1_N1) <- sec_per_hr * (-J_Diseased_Internalization_V121_R1_N1 + J_Diseased_Reaction_01 + J_Diseased_Abluminal_V121_R1_N1_1 + J_Diseased_Abluminal_V121_R1_N1_2 + J_Diseased_Abluminal_V121_R1_N1_3)
    d/dt(dis_R1_N1) <- sec_per_hr * (J_Diseased_Reaction_03 - J_Diseased_Internalization_R1_N1 + J_Diseased_Abluminal_R1_N1_1 + J_Diseased_Abluminal_R1_N1_2 + J_Diseased_Abluminal_R1_N1_3)
    d/dt(bld_R1_L_dis) <- sec_per_hr * (-J_Blood_Reaction_Disease_04 + J_Blood_Secretion_From_Diseased_R1 - J_Blood_Internalization_From_Diseased_R1 - J_Blood_Reaction_Disease_12 - J_Blood_Reaction_Disease_05 - J_Blood_Reaction_Diseased_a2M_05 - J_Blood_Reaction_Diseased_a2M_10 - J_Blood_Reaction_Diseased_a2M_fast_03 - J_Blood_Reaction_Diseased_a2M_fast_08 - J_Blood_Reaction_Disease_14 - J_Blood_Reaction_Diseased_a2M_14 - J_Blood_Reaction_Diseased_a2M_fast_12)
    d/dt(bld_N1_L_dis) <- sec_per_hr * (-J_Blood_Reaction_Disease_04 + J_Blood_Secretion_From_Diseased_N1 - J_Blood_Internalization_From_Diseased_N1 - J_Blood_Reaction_Disease_11 - J_Blood_Reaction_Disease_03 - J_Blood_Reaction_Disease_02 - J_Blood_Reaction_Disease_08 - J_Blood_Reaction_Disease_01 - J_Blood_Reaction_Diseased_a2M_07 - J_Blood_Reaction_Diseased_a2M_09 - J_Blood_Reaction_Diseased_a2M_13 - J_Blood_Reaction_Diseased_a2M_fast_05 - J_Blood_Reaction_Diseased_a2M_fast_07 - J_Blood_Reaction_Diseased_a2M_fast_11)
    d/dt(bld_R1_L_dis_N1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Disease_04 - J_Blood_Internalization_From_Diseased_R1_N1 - J_Blood_Reaction_Disease_09 - J_Blood_Reaction_Diseased_a2M_12 - J_Blood_Reaction_Diseased_a2M_fast_10)
    d/dt(bld_R2_L_dis) <- sec_per_hr * (J_Blood_Secretion_From_Diseased_R2 - J_Blood_Internalization_From_Diseased_R2 - J_Blood_Reaction_Disease_13 - J_Blood_Reaction_Disease_07 - J_Blood_Reaction_Disease_06 - J_Blood_Reaction_Diseased_a2M_06 - J_Blood_Reaction_Diseased_a2M_08 - J_Blood_Reaction_Diseased_a2M_11 - J_Blood_Reaction_Diseased_a2M_fast_04 - J_Blood_Reaction_Diseased_a2M_fast_06 - J_Blood_Reaction_Diseased_a2M_fast_09 - J_Blood_Reaction_Disease_15 - J_Blood_Reaction_Diseased_a2M_15 - J_Blood_Reaction_Diseased_a2M_fast_13)
    d/dt(bld_V165_R1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Disease_12 - J_Blood_Internalization_From_Diseased_V165_R1)
    d/dt(bld_V165_R2_L_dis) <- sec_per_hr * (J_Blood_Reaction_Disease_13 - J_Blood_Reaction_Disease_03 - J_Blood_Internalization_From_Diseased_V165_R2)
    d/dt(bld_V165_N1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Disease_11 - J_Blood_Reaction_Disease_06 - J_Blood_Internalization_From_Diseased_V165_N1)
    d/dt(bld_V121_R1_L_dis_N1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Disease_09 - J_Blood_Internalization_From_Diseased_V121_R1_N1 + J_Blood_Reaction_Disease_02)
    d/dt(bld_R1_L_dis_V121) <- sec_per_hr * (J_Blood_Reaction_Disease_05 - J_Blood_Internalization_From_Diseased_R1_V121 - J_Blood_Reaction_Disease_02)
    d/dt(bld_R2_L_dis_V121) <- sec_per_hr * (J_Blood_Reaction_Disease_07 - J_Blood_Internalization_From_Diseased_R2_V121)
    d/dt(bld_R2_L_dis_V165_N1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Disease_06 + J_Blood_Reaction_Disease_03 - J_Blood_Internalization_From_Diseased_R2_V165_N1)
    d/dt(bld_V165_R1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_12 - J_Blood_Internalization_From_Normal_V165_R1)
    d/dt(bld_V165_R2_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_13 - J_Blood_Reaction_Normal_03 - J_Blood_Internalization_From_Normal_V165_R2)
    d/dt(bld_V165_N1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_11 - J_Blood_Reaction_Normal_06 - J_Blood_Internalization_From_Normal_V165_N1)
    d/dt(bld_V121_R1_L_norm_N1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_09 - J_Blood_Internalization_From_Normal_V121_R1_N1 + J_Blood_Reaction_Normal_02)
    d/dt(bld_R1_L_norm_V121) <- sec_per_hr * (J_Blood_Reaction_Normal_05 - J_Blood_Internalization_From_Normal_R1_V121 - J_Blood_Reaction_Normal_02)
    d/dt(bld_R2_L_norm_V121) <- sec_per_hr * (J_Blood_Reaction_Normal_07 - J_Blood_Internalization_From_Normal_R2_V121)
    d/dt(bld_R2_L_norm_V165_N1_L_norm) <- sec_per_hr * (-J_Blood_Internalization_From_Normal_R2_N1 + J_Blood_Reaction_Normal_06 + J_Blood_Reaction_Normal_03)
    d/dt(nrm_MECM) <- sec_per_hr * (-J_Normal_Reaction_22 - J_Normal_Reaction_11 - J_Normal_Reaction_30)
    d/dt(nrm_V165_MECM) <- sec_per_hr * (J_Normal_Reaction_22)
    d/dt(nrm_MEBM) <- sec_per_hr * (-J_Normal_Reaction_21 - J_Normal_Reaction_10 - J_Normal_Reaction_29)
    d/dt(nrm_V165_MEBM) <- sec_per_hr * (J_Normal_Reaction_21)
    d/dt(nrm_MPBM) <- sec_per_hr * (-J_Normal_Reaction_23 - J_Normal_Reaction_12 - J_Normal_Reaction_31)
    d/dt(nrm_V165_MPBM) <- sec_per_hr * (J_Normal_Reaction_23)
    d/dt(nrm_V165_R1) <- sec_per_hr * (J_Normal_Reaction_26 - J_Normal_Internalization_V165_R1 + J_Normal_Abluminal_V165_R1_1)
    d/dt(nrm_V165_R2) <- sec_per_hr * (J_Normal_Reaction_27 - J_Normal_Reaction_05 - J_Normal_Internalization_V165_R2 + J_Normal_Abluminal_V165_R2_1 + J_Normal_Abluminal_V165_R2_2)
    d/dt(nrm_V165_N1) <- sec_per_hr * (J_Normal_Reaction_24 - J_Normal_Reaction_08 - J_Normal_Internalization_V165_N1 + J_Normal_Abluminal_V165_N1_1 + J_Normal_Abluminal_V165_N1_2)
    d/dt(nrm_V121_R1) <- sec_per_hr * (J_Normal_Reaction_07 - J_Normal_Reaction_04 - J_Normal_Internalization_R1_V121 + J_Normal_Abluminal_V121_R1_1 + J_Normal_Abluminal_V121_R1_2 + J_Normal_Abluminal_V121_R1_3)
    d/dt(nrm_V121_R2) <- sec_per_hr * (J_Normal_Reaction_09 - J_Normal_Internalization_R2_V121 + J_Normal_Abluminal_V121_R2_1)
    d/dt(nrm_V121_R1_N1) <- sec_per_hr * (J_Normal_Reaction_15 + J_Normal_Reaction_04 - J_Normal_Internalization_V121_R1_N1 + J_Normal_Abluminal_V121_R1_N1_1 + J_Normal_Abluminal_V121_R1_N1_2 + J_Normal_Abluminal_V121_R1_N1_3)
    d/dt(nrm_V165_R2_N1) <- sec_per_hr * (J_Normal_Reaction_08 + J_Normal_Reaction_05 - J_Normal_Internalization_R2_V165_N1 + J_Normal_Abluminal_V165_R2_N1_1 + J_Normal_Abluminal_V165_R2_N1_2 + J_Normal_Abluminal_V165_R2_N1_3)
    d/dt(dis_R1_N2) <- sec_per_hr * (-J_Diseased_Reaction_10)
    d/dt(dis_V121_R1_N2) <- sec_per_hr * (J_Diseased_Reaction_10)
    d/dt(nrm_sR1) <- sec_per_hr * (J_Normal_Secretion_sR1 - J_Normal_Reaction_11 - J_Normal_Reaction_10 - J_Normal_Reaction_12 - J_Normal_Reaction_28 - J_Normal_Reaction_19 - J_Normal_Reaction_13 - J_Normal_Lymphatic_Drainage_sR1 + J_Normal_Permeable_Transfer_From_Blood_sR1 - J_Normal_Reaction_34 - J_Normal_Degradation_sR1 - J_Normal_Permeable_Transfer_to_Blood_sR1)
    d/dt(dis_sR1) <- sec_per_hr * (J_Diseased_Secretion_sR1 - J_Diseased_Reaction_06 - J_Diseased_Reaction_05 - J_Diseased_Reaction_07 - J_Diseased_Reaction_21 - J_Diseased_Reaction_13 - J_Diseased_Reaction_08 - J_Diseased_Lymphatic_Drainage_sR1 + J_Diseased_Permeable_Transfer_From_Blood_sR1 - J_Diseased_Reaction_27 - J_Diseased_Degradation_sR1 - J_Diseased_Permeable_Transfer_To_Blood_sR1)
    d/dt(nrm_sR1_MECM) <- sec_per_hr * (J_Normal_Reaction_11)
    d/dt(nrm_sR1_MEBM) <- sec_per_hr * (J_Normal_Reaction_10)
    d/dt(nrm_sR1_MPBM) <- sec_per_hr * (J_Normal_Reaction_12)
    d/dt(nrm_V165_sR1) <- sec_per_hr * (J_Normal_Reaction_28 - J_Normal_Lymphatic_Drainage_V165_sR1 + J_Normal_Permeable_Transfer_From_Blood_V165_sR1 - J_Normal_Permeable_Transfer_to_Blood_V165_sR1 - J_Normal_Degradation_V165_sR1)
    d/dt(nrm_V121_sR1) <- sec_per_hr * (J_Normal_Reaction_19 - J_Normal_Degradation_V121_sR1 + J_Normal_Permeable_Transfer_From_Blood_V121_sR1 - J_Normal_Reaction_01 - J_Normal_Lymphatic_Drainage_V121_sR1 - J_Normal_Permeable_Transfer_to_Blood_V121_sR1)
    d/dt(nrm_sR1_N1) <- sec_per_hr * (J_Normal_Reaction_13 - J_Normal_Reaction_16 - J_Normal_Internalization_sR1_N1 + J_Normal_Abluminal_sR1_N1_1)
    d/dt(nrm_V121_sR1_N1) <- sec_per_hr * (J_Normal_Reaction_16 - J_Normal_Internalization_V121_sR1_N1 + J_Normal_Reaction_01 + J_Normal_Abluminal_V121_sR1_N1_1)
    d/dt(bld_sR1) <- sec_per_hr * (-J_Blood_Reaction_Soluble_02 - J_Blood_Reaction_Soluble_01 - J_Blood_Reaction_Normal_08 - J_Blood_Reaction_Disease_08 - J_Blood_Permeable_Transfer_To_Diseased_sR1 + J_Blood_Permeable_Transfer_From_Diseased_sR1 + J_Blood_Lymphatic_Drainage_From_Diseased_sR1 - J_Blood_Degradation_Clearance_sR1 - J_Blood_Reaction_Soluble_03 + J_Blood_Lymphatic_Drainage_From_Normal_sR1 + J_Blood_Secretion_From_Normal_sR1 + J_Blood_Secretion_From_Diseased_sR1 + J_Blood_Permeable_Transfer_From_Normal_sR1 - J_Blood_Permeable_Transfer_To_Normal_sR1)
    d/dt(bld_V165_sR1) <- sec_per_hr * (J_Blood_Reaction_Soluble_02 - J_Blood_Permeable_Transfer_To_Diseased_V165_sR1 + J_Blood_Permeable_Transfer_From_Diseased_V165_sR1 + J_Blood_Lymphatic_Drainage_From_Normal_V165_sR1 - J_Blood_Degradation_Clearance_V165_sR1 + J_Blood_Lymphatic_Drainage_From_Diseased_V165_sR1 - J_Blood_Permeable_Transfer_To_Normal_V165_sR1 + J_Blood_Permeable_Transfer_From_Normal_V165_sR1)
    d/dt(bld_V121_sR1) <- sec_per_hr * (J_Blood_Reaction_Soluble_01 - J_Blood_Permeable_Transfer_To_Diseased_V121_sR1 + J_Blood_Permeable_Transfer_From_Diseased_V121_sR1 + J_Blood_Lymphatic_Drainage_From_Normal_V121_sR1 - J_Blood_Degradation_Clearance_V121_sR1 - J_Blood_Reaction_Disease_01 - J_Blood_Reaction_Normal_01 + J_Blood_Lymphatic_Drainage_From_Diseased_V121_sR1 - J_Blood_Permeable_Transfer_To_Normal_V121_sR1 + J_Blood_Permeable_Transfer_From_Normal_V121_sR1)
    d/dt(bld_sR1_N1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_08 - J_Blood_Reaction_Normal_10 - J_Blood_Internalization_From_Normal_sR1_N1)
    d/dt(bld_sR1_N1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Disease_08 - J_Blood_Reaction_Disease_10 - J_Blood_Internalization_From_Diseased_sR1_N1)
    d/dt(bld_V121_sR1_N1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_10 - J_Blood_Internalization_From_Normal_V121_sR1_N1 + J_Blood_Reaction_Normal_01)
    d/dt(bld_V121_sR1_N1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Disease_10 - J_Blood_Internalization_From_Diseased_V121_sR1_N1 + J_Blood_Reaction_Disease_01)
    d/dt(dis_sR1_MECM) <- sec_per_hr * (J_Diseased_Reaction_06)
    d/dt(dis_sR1_MEBM) <- sec_per_hr * (J_Diseased_Reaction_05)
    d/dt(dis_sR1_MPBM) <- sec_per_hr * (J_Diseased_Reaction_07)
    d/dt(dis_V165_sR1) <- sec_per_hr * (J_Diseased_Reaction_21 - J_Diseased_Degradation_V165_sR1 + J_Diseased_Permeable_Transfer_From_Blood_V165_sR1 - J_Diseased_Lymphatic_Drainage_V165_sR1 - J_Diseased_Permeable_Transfer_To_Blood_V165_sR1)
    d/dt(dis_V121_sR1) <- sec_per_hr * (J_Diseased_Reaction_13 - J_Diseased_Lymphatic_Drainage_V121_sR1 + J_Diseased_Permeable_Transfer_From_Blood_V121_sR1 - J_Diseased_Reaction_30 - J_Diseased_Degradation_V121_sR1 - J_Diseased_Permeable_Transfer_To_Blood_V121_sR1)
    d/dt(dis_sR1_N1) <- sec_per_hr * (J_Diseased_Reaction_08 - J_Diseased_Reaction_11 - J_Diseased_Internalization_sR1_N1 + J_Diseased_Abluminal_sR1_N1_1)
    d/dt(dis_V121_sR1_N1) <- sec_per_hr * (J_Diseased_Reaction_11 - J_Diseased_Internalization_V121_sR1_N1 + J_Diseased_Reaction_30 + J_Diseased_Abluminal_V121_sR1_N1_1)
    d/dt(bld_a2M) <- sec_per_hr * (-J_Blood_Reaction_a2M_01 - J_Blood_Clearance_a2M + J_Blood_Synthesis_a2M - J_Blood_Reaction_a2M_02 - J_Blood_Reaction_a2M_03)
    d/dt(bld_V165_a2M) <- sec_per_hr * (J_Blood_Reaction_a2M_01 - J_Blood_Clearance_V165_a2M - J_Blood_Reaction_Normal_a2M_05 - J_Blood_Reaction_Normal_a2M_06 - J_Blood_Reaction_Normal_a2M_07 - J_Blood_Reaction_Diseased_a2M_05 - J_Blood_Reaction_Diseased_a2M_06 - J_Blood_Reaction_Diseased_a2M_07)
    d/dt(bld_V121_a2M) <- sec_per_hr * (-J_Blood_Clearance_V121_a2M - J_Blood_Reaction_Normal_a2M_10 - J_Blood_Reaction_Normal_a2M_11 - J_Blood_Reaction_Normal_a2M_12 - J_Blood_Reaction_Diseased_a2M_10 - J_Blood_Reaction_Diseased_a2M_11 - J_Blood_Reaction_Diseased_a2M_12 + J_Blood_Reaction_a2M_02)
    d/dt(bld_a2M_V165_R1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_05)
    d/dt(bld_a2M_V165_R2_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_06 - J_Blood_Reaction_Normal_a2M_09)
    d/dt(bld_a2M_V165_N1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_07 - J_Blood_Reaction_Normal_a2M_08)
    d/dt(bld_R2_L_norm_a2M_V165_N1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_08 + J_Blood_Reaction_Normal_a2M_09)
    d/dt(bld_R1_L_norm_V121_a2M) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_10 - J_Blood_Reaction_Normal_a2M_13)
    d/dt(bld_R2_L_norm_V121_a2M) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_11)
    d/dt(bld_a2M_V121_R1_L_norm_N1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_12 + J_Blood_Reaction_Normal_a2M_13)
    d/dt(bld_a2M_V165_R1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_05)
    d/dt(bld_a2M_V165_R2_L_dis) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_06 - J_Blood_Reaction_Diseased_a2M_09)
    d/dt(bld_a2M_V165_N1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_07 - J_Blood_Reaction_Diseased_a2M_08)
    d/dt(bld_R2_L_dis_a2M_V165_N1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_08 + J_Blood_Reaction_Diseased_a2M_09)
    d/dt(bld_R1_L_dis_V121_a2M) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_10 - J_Blood_Reaction_Diseased_a2M_13)
    d/dt(bld_R2_L_dis_V121_a2M) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_11)
    d/dt(bld_a2M_V121_R1_L_dis_N1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_12 + J_Blood_Reaction_Diseased_a2M_13)
    d/dt(bld_a2M_fast) <- sec_per_hr * (-J_Blood_Reaction_a2M_fast_01 - J_Blood_Reaction_a2M_fast_02 + J_Blood_Synthesis_a2M_fast - J_Blood_Clearance_a2M_fast - J_Blood_Reaction_a2M_fast_03)
    d/dt(bld_V165_a2M_fast) <- sec_per_hr * (J_Blood_Reaction_a2M_fast_01 - J_Blood_Clearance_V165_a2M_fast - J_Blood_Reaction_Normal_a2M_fast_03 - J_Blood_Reaction_Normal_a2M_fast_04 - J_Blood_Reaction_Normal_a2M_fast_05 - J_Blood_Reaction_Diseased_a2M_fast_03 - J_Blood_Reaction_Diseased_a2M_fast_04 - J_Blood_Reaction_Diseased_a2M_fast_05)
    d/dt(bld_V121_a2M_fast) <- sec_per_hr * (J_Blood_Reaction_a2M_fast_02 - J_Blood_Clearance_V121_a2M_fast - J_Blood_Reaction_Normal_a2M_fast_08 - J_Blood_Reaction_Normal_a2M_fast_09 - J_Blood_Reaction_Normal_a2M_fast_10 - J_Blood_Reaction_Diseased_a2M_fast_08 - J_Blood_Reaction_Diseased_a2M_fast_09 - J_Blood_Reaction_Diseased_a2M_fast_10)
    d/dt(bld_a2M_fast_V165_R1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_fast_03)
    d/dt(bld_a2M_fast_V165_R2_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_fast_04 - J_Blood_Reaction_Normal_a2M_fast_07)
    d/dt(bld_a2M_fast_V165_N1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_fast_05 - J_Blood_Reaction_Normal_a2M_fast_06)
    d/dt(bld_R2_L_norm_a2M_fast_V165_N1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_fast_06 + J_Blood_Reaction_Normal_a2M_fast_07)
    d/dt(bld_R1_L_norm_V121_a2M_fast) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_fast_08 - J_Blood_Reaction_Normal_a2M_fast_11)
    d/dt(bld_R2_L_norm_V121_a2M_fast) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_fast_09)
    d/dt(bld_a2M_fast_V121_R1_L_norm_N1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_fast_10 + J_Blood_Reaction_Normal_a2M_fast_11)
    d/dt(bld_a2M_fast_V165_R1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_fast_03)
    d/dt(bld_a2M_fast_V165_R2_L_dis) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_fast_04 - J_Blood_Reaction_Diseased_a2M_fast_07)
    d/dt(bld_a2M_fast_V165_N1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_fast_05 - J_Blood_Reaction_Diseased_a2M_fast_06)
    d/dt(bld_R2_L_dis_a2M_fast_V165_N1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_fast_06 + J_Blood_Reaction_Diseased_a2M_fast_07)
    d/dt(bld_R1_L_dis_V121_a2M_fast) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_fast_08 - J_Blood_Reaction_Diseased_a2M_fast_11)
    d/dt(bld_R2_L_dis_V121_a2M_fast) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_fast_09)
    d/dt(bld_a2M_fast_V121_R1_L_dis_N1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_fast_10 + J_Blood_Reaction_Diseased_a2M_fast_11)
    d/dt(nrm_V165b) <- sec_per_hr * (-J_Normal_Reaction_31 - J_Normal_Reaction_30 - J_Normal_Reaction_29 - J_Normal_Reaction_32 - J_Normal_Reaction_33 - J_Normal_Reaction_34 - J_Normal_Degradation_V165b + J_Normal_Secretion_V165b_1 - J_Normal_Lymphatic_Drainage_V165b + J_Normal_Permeable_Transfer_From_Blood_V165b - J_Normal_Permeable_Transfer_To_Blood_V165b + J_Normal_Secretion_V165b_2)
    d/dt(nrm_V165b_MPBM) <- sec_per_hr * (J_Normal_Reaction_31)
    d/dt(nrm_V165b_MECM) <- sec_per_hr * (J_Normal_Reaction_30)
    d/dt(nrm_V165b_MEBM) <- sec_per_hr * (J_Normal_Reaction_29)
    d/dt(dis_V165b) <- sec_per_hr * (-J_Diseased_Reaction_24 - J_Diseased_Reaction_23 - J_Diseased_Reaction_22 - J_Diseased_Reaction_25 - J_Diseased_Reaction_26 - J_Diseased_Reaction_27 - J_Diseased_Degradation_V165b + J_Diseased_Secretion_V165b_1 - J_Diseased_Lymphatic_Drainage_V165b + J_Diseased_Permeable_Transfer_From_Blood_V165b - J_Diseased_Permeable_Transfer_To_Blood_V165b + J_Diseased_Secretion_V165b_2)
    d/dt(dis_V165b_MPBM) <- sec_per_hr * (J_Diseased_Reaction_24)
    d/dt(dis_V165b_MECM) <- sec_per_hr * (J_Diseased_Reaction_23)
    d/dt(dis_V165b_MEBM) <- sec_per_hr * (J_Diseased_Reaction_22)
    d/dt(bld_V165b) <- sec_per_hr * (-J_Blood_Reaction_Normal_14 - J_Blood_Reaction_Normal_15 - J_Blood_Reaction_Soluble_03 - J_Blood_Reaction_Disease_14 - J_Blood_Reaction_Disease_15 + J_Blood_Secretion_From_Normal_V165b - J_Blood_Clearance_V165b + J_Blood_Secretion_From_Diseased_V165b - J_Blood_Permeable_Transfer_To_Normal_V165b - J_Blood_Permeable_Transfer_To_Diseased_V165b + J_Blood_Lymphatic_Drainage_From_Normal_V165b + J_Blood_Lymphatic_Drainage_From_Diseased_V165b + J_Blood_Permeable_Transfer_From_Normal_V165b + J_Blood_Permeable_Transfer_From_Diseased_V165b - J_Blood_Reaction_a2M_03 - J_Blood_Reaction_a2M_fast_03)
    d/dt(bld_V165b_R1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_14 - J_Blood_Internalization_From_Normal_V165b_R1)
    d/dt(bld_V165b_R2_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_15 - J_Blood_Internalization_From_Normal_V165b_R2)
    d/dt(bld_V165b_sR1) <- sec_per_hr * (J_Blood_Reaction_Soluble_03 - J_Blood_Degradation_Clearance_V165b_sR1 + J_Blood_Lymphatic_Drainage_From_Diseased_V165b_sR1 + J_Blood_Lymphatic_Drainage_From_Normal_V165b_sR1 + J_Blood_Permeable_Transfer_From_Diseased_V165b_sR1 + J_Blood_Permeable_Transfer_From_Normal_V165b_sR1 - J_Blood_Permeable_Transfer_To_Diseased_V165b_sR1 - J_Blood_Permeable_Transfer_To_Normal_V165b_sR1)
    d/dt(nrm_V165b_R1) <- sec_per_hr * (J_Normal_Reaction_32 + J_Normal_Abluminal_V165b_R1_1 - J_Normal_Internalization_V165b_R1)
    d/dt(nrm_V165b_R2) <- sec_per_hr * (J_Normal_Reaction_33 + J_Normal_Abluminal_V165b_R2_1 - J_Normal_Internalization_V165b_R2)
    d/dt(nrm_V165b_sR1) <- sec_per_hr * (J_Normal_Reaction_34 - J_Normal_Degradation_V165b_sR1 - J_Normal_Lymphatic_Drainage_V165b_sR1 + J_Normal_Permeable_Transfer_From_Blood_V165b_sR1 - J_Normal_Permeable_Transfer_to_Blood_V165b_sR1)
    d/dt(dis_V165b_R1) <- sec_per_hr * (J_Diseased_Reaction_25 + J_Diseased_Abluminal_V165b_R1_1 - J_Diseased_Internalization_V165b_R1)
    d/dt(dis_V165b_R2) <- sec_per_hr * (J_Diseased_Reaction_26 + J_Diseased_Abluminal_V165b_R2_1 - J_Diseased_Internalization_V165b_R2)
    d/dt(dis_V165b_sR1) <- sec_per_hr * (J_Diseased_Reaction_27 - J_Diseased_Degradation_V165b_sR1 - J_Diseased_Lymphatic_Drainage_V165b_sR1 + J_Diseased_Permeable_Transfer_From_Blood_V165b_sR1 - J_Diseased_Permeable_Transfer_To_Blood_V165b_sR1)
    d/dt(bld_V165b_R1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Disease_14 - J_Blood_Internalization_From_Diseased_V165b_R1)
    d/dt(bld_V165b_R2_L_dis) <- sec_per_hr * (J_Blood_Reaction_Disease_15 - J_Blood_Internalization_From_Diseased_V165b_R2)
    d/dt(bld_V165b_a2M) <- sec_per_hr * (J_Blood_Reaction_a2M_03 - J_Blood_Reaction_Diseased_a2M_14 - J_Blood_Reaction_Diseased_a2M_15 - J_Blood_Reaction_Normal_a2M_14 - J_Blood_Reaction_Normal_a2M_15 - J_Blood_Clearance_V165b_a2M)
    d/dt(bld_V165b_a2M_fast) <- sec_per_hr * (J_Blood_Reaction_a2M_fast_03 - J_Blood_Reaction_Diseased_a2M_fast_12 - J_Blood_Reaction_Diseased_a2M_fast_13 - J_Blood_Reaction_Normal_a2M_fast_12 - J_Blood_Reaction_Normal_a2M_fast_13 - J_Blood_Clearance_V165b_a2M_fast)
    d/dt(bld_a2M_V165b_R1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_14)
    d/dt(bld_a2M_V165b_R2_L_dis) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_15)
    d/dt(bld_a2M_fast_V165b_R1_L_dis) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_fast_12)
    d/dt(bld_a2M_fast_V165b_R2_L_dis) <- sec_per_hr * (J_Blood_Reaction_Diseased_a2M_fast_13)
    d/dt(bld_a2M_V165b_R1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_14)
    d/dt(bld_a2M_V165b_R2_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_15)
    d/dt(bld_a2M_fast_V165b_R1_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_fast_12)
    d/dt(bld_a2M_fast_V165b_R2_L_norm) <- sec_per_hr * (J_Blood_Reaction_Normal_a2M_fast_13)

    # Output rules (pM in the interstitial fluid / plasma, receptors per cell)
    bld_R1_L_dis_num <- ((bld_R1_L_dis * VesselCellSurfArea_blood * Avogadro) / VesselSA_Vol_blood_tis_dis)
    bld_R2_L_dis_num <- ((bld_R2_L_dis * VesselCellSurfArea_blood * Avogadro) / VesselSA_Vol_blood_tis_dis)
    bld_N1_L_dis_num <- ((bld_N1_L_dis * VesselCellSurfArea_blood * Avogadro) / VesselSA_Vol_blood_tis_dis)
    dis_R1_num <- ((dis_R1 * VesselCellSurfArea_tis_dis * Avogadro) / VesselSA_Vol_tis_dis)
    dis_R2_num <- ((dis_R2 * VesselCellSurfArea_tis_dis * Avogadro) / VesselSA_Vol_tis_dis)
    dis_N1_num <- ((dis_N1 * VesselCellSurfArea_tis_dis * Avogadro) / VesselSA_Vol_tis_dis)
    nrm_V165_pM <- (((nrm_V165 / Interstit_Vol_tis_norm) * 1000) / 1e-12)
    bld_V165_pM <- (((bld_V165 / Interstit_Vol_blood) * 1000) / 1e-12)
    bld_V121_pM <- (((bld_V121 / Interstit_Vol_blood) * 1000) / 1e-12)
    dis_V165_pM <- (((dis_V165 / Interstit_Vol_tis_dis) * mL_to_L) / 1e-12)
    dis_V121_pM <- (((dis_V121 / Interstit_Vol_tis_dis) * 1000) / 1e-12)
    bld_a2M_fast_nM <- ((bld_a2M_fast / Interstit_Vol_blood) * mL_to_L * M_to_nM)
    bld_a2M_uM <- ((bld_a2M / Interstit_Vol_blood) * mL_to_L * M_to_uM)
    nrm_sR1_pM <- (((nrm_sR1 / Interstit_Vol_tis_norm) * mL_to_L) / 1e-12)
    bld_sR1_pM <- (((bld_sR1 / Interstit_Vol_blood) * mL_to_L) / 1e-12)
    dis_sR1_pM <- ((dis_sR1 / Interstit_Vol_tis_dis) * mL_to_L * M_to_pM)
    nrm_hVEGF_sR1_pM <- ((((nrm_V165_sR1 + nrm_V165b_sR1 + nrm_V121_sR1) / Interstit_Vol_tis_norm) * mL_to_L) / 1e-12)
    dis_hVEGF_sR1_pM <- (((dis_V165_sR1 + dis_V165b_sR1 + dis_V121_sR1) / Interstit_Vol_tis_dis) * mL_to_L * M_to_pM)
    bld_V165b_R1 <- (bld_V165b_R1_L_norm + bld_V165b_R1_L_dis)
    bld_V165b_R2 <- (bld_V165b_R2_L_norm + bld_V165b_R2_L_dis)
    nrm_V165b_pM <- ((nrm_V165b / Interstit_Vol_tis_norm) * mL_to_L * M_to_pM)
    dis_V165b_pM <- ((dis_V165b / Interstit_Vol_tis_dis) * mL_to_L * M_to_pM)
    bld_V165b_pM <- ((bld_V165b / Interstit_Vol_blood) * mL_to_L * M_to_pM)
    bld_R2 <- (bld_R2_L_norm + bld_R2_L_dis)
    bld_R1 <- (bld_R1_L_norm + bld_R1_L_dis)
    bld_hVEGF_sR1_pM <- (((bld_V165_sR1 + bld_V165b_sR1 + bld_V121_sR1) / Interstit_Vol_blood) * mL_to_L * M_to_pM)
    nrm_V121_pM <- (((nrm_V121 / Interstit_Vol_tis_norm) * 1000) / 1e-12)
    nrm_VEGF165_121_total_norm_pM <- (nrm_V165_pM + nrm_V165b_pM + nrm_V121_pM)
    bld_VEGF165_121_total_blood_pM <- (bld_V165_pM + bld_V165b_pM + bld_V121_pM)
    dis_VEGF165_121_total_tis_dis_pM <- (dis_V165_pM + dis_V165b_pM + dis_V121_pM)
  })
}
