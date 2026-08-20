Proctor_2026_durvalumab_pbpk <- function() {
  description <- "PBPK / QSP (whole-body two-pore, 15 organs, 219 ODE states). Durvalumab (anti-PD-L1) with cachexia-driven time-dependent clearance, jointly describing endogenous albumin, endogenous IgG and the dosed antibody through a shared, competable FcRn recycling pathway. The endosomal degradation rate constant of FcRn-unbound protein decays exponentially over time, reproducing the ~12% fall in durvalumab clearance seen in cancer patients whose cachexia improves on treatment. Extends the Liu 2024 translational two-pore PBPK model; only the three parameters of the kdeg(t) decay were fitted, to longitudinal serum albumin (no clinical PK was used in the fit)."
  reference <- "Proctor JR, Wong H. Albumin Levels Are Predictive of Cachexia-Induced Time-Dependent Clearance of Therapeutic Antibodies: A Physiologically Based Pharmacokinetic Model of Durvalumab. CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70185. doi:10.1002/psp4.70185"
  vignette <- "Proctor_2026_durvalumab_cachexia"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # Structural equations are Appendix S1 Eq S1-S30; every physiological constant
  # and initial condition below is transcribed from the complete SimBiology model
  # export shipped as Appendix S3 (Compartments / Species / Parameters / Rules /
  # Reactions sheets), which supplies the values at full precision where main-text
  # Table S1 / S2 round them. See the vignette source-trace table.
  #
  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: every state below was checked against the
  # Appendix S3 Species sheet, which names each species and its compartment scope.
  compartmentData <- list(
    vp_lung                = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_heart               = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_muscle              = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_skin                = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_adipose             = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_bone                = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_brain               = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_kidney              = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_liver               = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_small_intestine     = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_large_intestine     = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_pancreas            = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_thymus              = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_spleen              = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_other               = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    is_lung                = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    is_heart               = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    is_muscle              = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    is_skin                = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    is_adipose             = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    is_bone                = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    is_brain               = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    is_kidney              = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    is_liver               = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    is_small_intestine     = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    is_large_intestine     = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    is_pancreas            = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    is_thymus              = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    is_spleen              = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    is_other               = list(analyte = "durvalumab", units = "nmol", specimen = "tissue", verified = TRUE),
    eu_lung                = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_heart               = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_muscle              = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_skin                = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_adipose             = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_bone                = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_brain               = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_kidney              = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_liver               = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_small_intestine     = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_large_intestine     = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_pancreas            = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_thymus              = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_spleen              = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_other               = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_lung                = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_heart               = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_muscle              = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_skin                = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_adipose             = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_bone                = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_brain               = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_kidney              = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_liver               = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_small_intestine     = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_large_intestine     = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_pancreas            = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_thymus              = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_spleen              = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_other               = list(analyte = "durvalumab", units = "nmol", specimen = "endosome", verified = TRUE),
    plasma                 = list(analyte = "durvalumab", units = "nmol", specimen = "plasma", verified = TRUE),
    lnode                  = list(analyte = "durvalumab", units = "nmol", specimen = "lymph", verified = TRUE),
    urine                  = list(analyte = "durvalumab", units = "nmol", specimen = "urine", verified = TRUE),
    vp_lung_alb            = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_heart_alb           = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_muscle_alb          = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_skin_alb            = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_adipose_alb         = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_bone_alb            = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_brain_alb           = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_kidney_alb          = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_liver_alb           = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_small_intestine_alb = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_large_intestine_alb = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_pancreas_alb        = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_thymus_alb          = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_spleen_alb          = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_other_alb           = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    is_lung_alb            = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    is_heart_alb           = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    is_muscle_alb          = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    is_skin_alb            = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    is_adipose_alb         = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    is_bone_alb            = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    is_brain_alb           = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    is_kidney_alb          = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    is_liver_alb           = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    is_small_intestine_alb = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    is_large_intestine_alb = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    is_pancreas_alb        = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    is_thymus_alb          = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    is_spleen_alb          = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    is_other_alb           = list(analyte = "albumin", units = "nmol", specimen = "tissue", verified = TRUE),
    eu_lung_alb            = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_heart_alb           = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_muscle_alb          = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_skin_alb            = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_adipose_alb         = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_bone_alb            = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_brain_alb           = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_kidney_alb          = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_liver_alb           = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_small_intestine_alb = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_large_intestine_alb = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_pancreas_alb        = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_thymus_alb          = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_spleen_alb          = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_other_alb           = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_lung_alb            = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_heart_alb           = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_muscle_alb          = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_skin_alb            = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_adipose_alb         = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_bone_alb            = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_brain_alb           = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_kidney_alb          = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_liver_alb           = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_small_intestine_alb = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_large_intestine_alb = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_pancreas_alb        = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_thymus_alb          = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_spleen_alb          = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_other_alb           = list(analyte = "albumin", units = "nmol", specimen = "endosome", verified = TRUE),
    plasma_alb             = list(analyte = "albumin", units = "nmol", specimen = "plasma", verified = TRUE),
    lnode_alb              = list(analyte = "albumin", units = "nmol", specimen = "lymph", verified = TRUE),
    urine_alb              = list(analyte = "albumin", units = "nmol", specimen = "urine", verified = TRUE),
    vp_lung_igg            = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_heart_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_muscle_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_skin_igg            = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_adipose_igg         = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_bone_igg            = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_brain_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_kidney_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_liver_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_small_intestine_igg = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_large_intestine_igg = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_pancreas_igg        = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_thymus_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_spleen_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    vp_other_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    is_lung_igg            = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    is_heart_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    is_muscle_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    is_skin_igg            = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    is_adipose_igg         = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    is_bone_igg            = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    is_brain_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    is_kidney_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    is_liver_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    is_small_intestine_igg = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    is_large_intestine_igg = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    is_pancreas_igg        = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    is_thymus_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    is_spleen_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    is_other_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "tissue", verified = TRUE),
    eu_lung_igg            = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_heart_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_muscle_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_skin_igg            = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_adipose_igg         = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_bone_igg            = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_brain_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_kidney_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_liver_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_small_intestine_igg = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_large_intestine_igg = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_pancreas_igg        = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_thymus_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_spleen_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_other_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_lung_igg            = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_heart_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_muscle_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_skin_igg            = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_adipose_igg         = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_bone_igg            = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_brain_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_kidney_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_liver_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_small_intestine_igg = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_large_intestine_igg = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_pancreas_igg        = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_thymus_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_spleen_igg          = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_other_igg           = list(analyte = "endogenous IgG", units = "nmol", specimen = "endosome", verified = TRUE),
    plasma_igg             = list(analyte = "endogenous IgG", units = "nmol", specimen = "plasma", verified = TRUE),
    lnode_igg              = list(analyte = "endogenous IgG", units = "nmol", specimen = "lymph", verified = TRUE),
    urine_igg              = list(analyte = "endogenous IgG", units = "nmol", specimen = "urine", verified = TRUE),
    fr_lung_igg            = list(analyte = "unoccupied FcRn (IgG binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_heart_igg           = list(analyte = "unoccupied FcRn (IgG binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_muscle_igg          = list(analyte = "unoccupied FcRn (IgG binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_skin_igg            = list(analyte = "unoccupied FcRn (IgG binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_adipose_igg         = list(analyte = "unoccupied FcRn (IgG binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_bone_igg            = list(analyte = "unoccupied FcRn (IgG binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_brain_igg           = list(analyte = "unoccupied FcRn (IgG binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_kidney_igg          = list(analyte = "unoccupied FcRn (IgG binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_liver_igg           = list(analyte = "unoccupied FcRn (IgG binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_small_intestine_igg = list(analyte = "unoccupied FcRn (IgG binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_large_intestine_igg = list(analyte = "unoccupied FcRn (IgG binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_pancreas_igg        = list(analyte = "unoccupied FcRn (IgG binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_thymus_igg          = list(analyte = "unoccupied FcRn (IgG binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_spleen_igg          = list(analyte = "unoccupied FcRn (IgG binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_other_igg           = list(analyte = "unoccupied FcRn (IgG binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_lung_alb            = list(analyte = "unoccupied FcRn (albumin binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_heart_alb           = list(analyte = "unoccupied FcRn (albumin binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_muscle_alb          = list(analyte = "unoccupied FcRn (albumin binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_skin_alb            = list(analyte = "unoccupied FcRn (albumin binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_adipose_alb         = list(analyte = "unoccupied FcRn (albumin binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_bone_alb            = list(analyte = "unoccupied FcRn (albumin binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_brain_alb           = list(analyte = "unoccupied FcRn (albumin binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_kidney_alb          = list(analyte = "unoccupied FcRn (albumin binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_liver_alb           = list(analyte = "unoccupied FcRn (albumin binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_small_intestine_alb = list(analyte = "unoccupied FcRn (albumin binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_large_intestine_alb = list(analyte = "unoccupied FcRn (albumin binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_pancreas_alb        = list(analyte = "unoccupied FcRn (albumin binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_thymus_alb          = list(analyte = "unoccupied FcRn (albumin binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_spleen_alb          = list(analyte = "unoccupied FcRn (albumin binding site)", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_other_alb           = list(analyte = "unoccupied FcRn (albumin binding site)", units = "nmol", specimen = "endosome", verified = TRUE)
  )

  # No covariates: the model is a typical-subject mechanistic simulation. Serum
  # albumin is an OUTPUT of this model (the biomarker whose time course drives
  # kdeg), not an input covariate -- which is the paper's central point and the
  # reason ALB does not appear here as it does in empirical durvalumab popPK
  # models such as Ogasawara_2020_durvalumab.
  covariateData <- list()

  covariatesDataExcluded <- list(
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Serum albumin is the biomarker the paper is about, but it enters this",
        "model as a simulated STATE (the `albumin` output, in g/L) rather than as",
        "a covariate on a clearance parameter. Recorded here so that a reader",
        "searching the register for albumin-driven durvalumab clearance finds this",
        "model; contrast Ogasawara_2020_durvalumab, where ALB is a genuine",
        "power-model covariate on CL."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    disease_state  = "Adults with advanced cancer receiving durvalumab; the cohort exhibits cachexia that improves during treatment (rising serum albumin).",
    dose_range     = "Simulated at 10 mg/kg IV every 2 weeks for 60 weeks (paper Methods 2.2); the Appendix S3 dose object also stores a 1 mg/kg (473.33 nmol) weekly regimen.",
    weight_range   = "71 kg reference adult (inherited from the Shah & Betts 2012 / Liu 2024 human physiology; recovered exactly from the Appendix S3 dose conversion 1 mg/kg = 473.3333 nmol at MW 150000 g/mol).",
    notes          = paste(
      "The three kdeg(t) parameters were estimated by fminsearch against",
      "longitudinal mean serum albumin digitized (WebPlotDigitizer v4.0) from the",
      "durvalumab popPK analysis of Baverel 2018; NO clinical pharmacokinetic data",
      "entered the fit. All remaining parameters were fixed from Liu 2024. Because",
      "the fit used digitized summary means rather than individual data, the",
      "reported RSEs (all < 10%) do not carry subject-level variance.",
      "An alternative implementation with time-varying total endosomal FcRn was",
      "rejected by the authors on RMSE (Table S4); only the kdeg version is packaged."
    ),
    scope_note     = paste(
      "Mechanistic typical-value model: no IIV and no informative residual error",
      "(the paper reports a fit RMSE of 0.245 g/L for albumin, not a residual-error",
      "model). Intended for deterministic simulation. Dosing is in nmol into the",
      "`plasma` compartment: dose_nmol = dose_mg / 150000 * 1e6.",
      "Simulations must include a pre-treatment equilibration window of length",
      "`tequil` (default 2400 h = 100 days, the Appendix S3 `equil_time`); the",
      "kdeg decay begins at t = tequil, so the first dose belongs at t = tequil."
    )
  )

  ini({
    # ---- The only three estimated parameters (Table 1) ----------------------
    # Table 1 prints these with a "1/day" units column, but the model is
    # parameterised in 1/HOUR: the Appendix S3 Parameters sheet declares
    # kdeg_inf / kdeg_T in "1 / hour", and at kdeg = 18.318 /h the simulated
    # steady-state plasma albumin lands on Table S2's CssALB (6.0e-4 M = 40.2
    # g/L) to three digits, whereas a 1/day reading gives ~66 g/L, which is
    # physiologically impossible. Recorded as an erratum in the vignette.
    lkdeg_exp_inf       <- log(18.318)      ; label("Time-invariant endosomal degradation rate constant kdeg,inf (1/h)")        # Table 1 (SE 0.292, RSE 1.6%); Appendix S3 kdeg_inf = 18.318
    lkdeg_exp_component <- log(9.4793)      ; label("Time-varying component of the degradation rate constant kdeg,T (1/h)")     # Table 1 (SE 0.110, RSE 1.2%); Appendix S3 kdeg_T = 9.4793
    lkdeg_exp_kdes      <- log(0.010019/24) ; label("Rate constant of the decline in kdeg, k_des (1/h)")                        # Table 1 k_tv 0.0100 1/day (SE 0.000621, RSE 6.2%); Appendix S3 ktv = 1.0019e-2 /day, converted to /h

    # ---- Fixed system parameters (Table S2 / Appendix S3 Parameters) --------
    lclup    <- fixed(log(1.22))      ; label("Pinocytosis uptake clearance per unit endosomal volume (1/h)")      # Table S2 CL_up 1.22 1/hour
    lfcrn_igg <- fixed(log(62.67))    ; label("Endosomal FcRn concentration, IgG binding site (uM)")               # Table S2 FcRn_IGG 6.267e-5 M
    lfcrn_alb <- fixed(log(125.34))   ; label("Endosomal FcRn concentration, albumin binding site (uM)")           # Table S2 FcRn_ALB 1.2534e-4 M
    lcss_alb <- fixed(log(600))       ; label("Steady-state plasma albumin concentration (uM)")                    # Table S2 Css,ALB 6.0e-4 M = 40.2 g/L at MW 67000
    lcss_igg <- fixed(log(66))        ; label("Steady-state plasma endogenous IgG concentration (uM)")             # Table S2 Css,IGG 6.6e-5 M = 9.9 g/L at MW 150000

    # Pre-treatment equilibration window. kdeg is held at its baseline
    # kdeg(0) = kdeg,inf + kdeg,T until t = tequil, after which it decays; this
    # is the max(0, time - equil_time) in Appendix S3 Rule_1. Dose at t = tequil.
    tequil   <- fixed(2400)           ; label("Pre-treatment equilibration window before kdeg starts to decline (h)")  # Appendix S3 equil_time = 100 day

    # The paper fitted digitized summary means and reports no residual-error
    # model for either endpoint, so residual error is fixed at zero per the
    # standing policy for unreported RUV (documented in the vignette Errata).
    propSd   <- fixed(0)              ; label("Proportional residual error on durvalumab plasma concentration (fraction)")
  })

  model({
    # Time-varying endosomal degradation rate constant of FcRn-unbound
    # protein, Eq S16 / SimBiology Rule_1 (Appendix S3, Rules sheet):
    #   kdeg = kdeg_exp_inf + kdeg_exp_component * exp(-kdes * max(0, t - tequil))
    # tequil is the pre-treatment equilibration window during which kdeg is
    # held at its baseline kdeg(0) = kdegInf + kdegT (cachectic state).
    kdegInf <- exp(lkdeg_exp_inf)
    kdegT   <- exp(lkdeg_exp_component)
    kdes    <- exp(lkdeg_exp_kdes)
    tdecl   <- max(0, t - tequil)
    kdeg    <- kdegInf + kdegT * exp(-kdes * tdecl)

    # ================= Physiological constants =================
    # Organ volumes (mL) and flows (mL/h): Table S1 / Appendix S3 Compartments
    # and Parameters sheets. vv=vascular plasma, vi=interstitial, ve=endosomal.
    vv_lung <- 55.0; vi_lung <- 300.0; ve_lung <- 5.0
    vv_heart <- 13.1; vi_heart <- 48.8; ve_heart <- 1.71
    vv_muscle <- 662.0; vi_muscle <- 3910.0; ve_muscle <- 150.0
    vv_skin <- 127.0; vi_skin <- 1125.0; ve_skin <- 17.0
    vv_adipose <- 148.0; vi_adipose <- 2289.0; ve_adipose <- 67.3
    vv_bone <- 224.0; vi_bone <- 1891.0; ve_bone <- 50.8
    vv_brain <- 31.9; vi_brain <- 261.0; ve_brain <- 7.25
    vv_kidney <- 18.2; vi_kidney <- 49.8; ve_kidney <- 1.66
    vv_liver <- 183.0; vi_liver <- 429.0; ve_liver <- 10.7
    vv_small_intestine <- 6.15; vi_small_intestine <- 67.1; ve_small_intestine <- 1.93
    vv_large_intestine <- 8.74; vi_large_intestine <- 95.3; ve_large_intestine <- 2.74
    vv_pancreas <- 5.7; vi_pancreas <- 18.0; ve_pancreas <- 0.518
    vv_thymus <- 0.353; vi_thymus <- 1.09; ve_thymus <- 0.0321
    vv_spleen <- 26.8; vi_spleen <- 44.3; ve_spleen <- 1.11
    vv_other <- 204.0; vi_other <- 831.0; ve_other <- 24.3

    q_lung <- 181913.0; j_lung <- 363.826
    q_heart <- 7752.0; j_heart <- 15.504
    q_muscle <- 33469.0; j_muscle <- 66.938
    q_skin <- 11626.0; j_skin <- 23.252
    q_adipose <- 11233.0; j_adipose <- 22.466
    q_bone <- 2591.0; j_bone <- 5.182
    q_brain <- 21453.0; j_brain <- 42.906
    q_kidney <- 36402.0; j_kidney <- 72.804
    q_liver <- 13210.0; j_liver <- 95.549
    q_small_intestine <- 12368.0; j_small_intestine <- 24.736
    q_large_intestine <- 12867.0; j_large_intestine <- 25.734
    q_pancreas <- 3056.0; j_pancreas <- 6.112
    q_thymus <- 353.0; j_thymus <- 0.706
    q_spleen <- 6343.0; j_spleen <- 12.686
    q_other <- 5521.0; j_other <- 11.042

    vpl   <- 1412.0569999999998   # central plasma, residual of total plasma minus organ vascular volumes
    vln   <- 274.0     # lymph node
    vbl   <- 31.54   # renal-filtrate (bladder) holding volume
    lnlf  <- 363.826   # lymph-node outflow to plasma (mL/h) = q_lung/500
    sigis <- 0.2    # lymphatic reflection coefficient
    fr    <- 0.715      # fraction of FcRn-bound recycled to vascular side
    gfr   <- 7200.0    # glomerular filtration rate (mL/h)
    qurine <- 58.0  # urine flow (mL/h)

    clup  <- exp(lclup)   # pinocytosis uptake per unit endosomal volume (1/h)
    clup_lung <- clup * ve_lung
    clup_heart <- clup * ve_heart
    clup_muscle <- clup * ve_muscle
    clup_skin <- clup * ve_skin
    clup_adipose <- clup * ve_adipose
    clup_bone <- clup * ve_bone
    clup_brain <- clup * ve_brain
    clup_kidney <- clup * ve_kidney
    clup_liver <- clup * ve_liver
    clup_small_intestine <- clup * ve_small_intestine
    clup_large_intestine <- clup * ve_large_intestine
    clup_pancreas <- clup * ve_pancreas
    clup_thymus <- clup * ve_thymus
    clup_spleen <- clup * ve_spleen
    clup_other <- clup * ve_other

    # Total arterial plasma flow leaving the lung to all systemic organs.
    sumq_sys <- q_heart + q_muscle + q_skin + q_adipose + q_bone + q_brain + q_kidney + q_liver + q_small_intestine + q_large_intestine + q_pancreas + q_thymus + q_spleen + q_other
    qout_liver <- (q_liver - j_liver) + (q_spleen - j_spleen) + (q_pancreas - j_pancreas) + (q_small_intestine - j_small_intestine) + (q_large_intestine - j_large_intestine)

    # ================= Substance-specific constants =================
    # durvalumab (exogenous IgG1)
    spino_mab <- 1.0
    kon_mab   <- 560.0    # 1/(uM*h), converted from 560000000.0 1/(M*h)
    koff_mab  <- 23.9
    clren_mab <- gfr * 1.836703796979123e-08 * (1 - 0.0)

    # endogenous serum albumin
    spino_alb <- 0.185
    kon_alb   <- 27.0    # 1/(uM*h), converted from 27000000.0 1/(M*h)
    koff_alb  <- 30.24
    clren_alb <- gfr * 0.000400730683682927 * (1 - 0.97)

    # endogenous IgG
    spino_igg <- 1.0
    kon_igg   <- 560.0    # 1/(uM*h), converted from 560000000.0 1/(M*h)
    koff_igg  <- 23.9
    clren_igg <- gfr * 1.836703796979123e-08 * (1 - 0.0)

    # Two-pore extravasation clearance CL_TP,L + CL_TP,S per organ (mL/h),
    # Eq S21-S30. These are time-invariant (SimBiology initialAssignment), so
    # they are evaluated here as folded numeric constants; the generating
    # formula is documented in the vignette source-trace table.
    cldif_mab_lung <- 0.6231353179072838; clcnv_mab_lung <- 126.3374758731036
    cldif_mab_heart <- 0.026554149425369638; clcnv_mab_heart <- 5.383717012903415
    cldif_mab_muscle <- 0.11464664952498664; clcnv_mab_muscle <- 23.24401763478643
    cldif_mab_skin <- 0.03982437322230997; clcnv_mab_skin <- 8.074186531477698
    cldif_mab_adipose <- 0.03847816827853163; clcnv_mab_adipose <- 7.801250413563476
    cldif_mab_bone <- 0.00887536134689534; clcnv_mab_bone <- 1.7994337952054633
    cldif_mab_brain <- 0.07348634773251482; clcnv_mab_brain <- 14.89897846721065
    cldif_mab_kidney <- 0.1246935174641777; clcnv_mab_kidney <- 25.28096835703175
    cldif_mab_liver <- 0.16364953711588245; clcnv_mab_liver <- 33.179100675045696
    cldif_mab_small_intestine <- 0.04236606296348964; clcnv_mab_small_intestine <- 8.589501034002765
    cldif_mab_large_intestine <- 0.044075366441722254; clcnv_mab_large_intestine <- 8.936053509420569
    cldif_mab_pancreas <- 0.010468199257472862; clcnv_mab_pancreas <- 2.1223734767070224
    cldif_mab_thymus <- 0.001209186628890026; clcnv_mab_thymus <- 0.24515636036569988
    cldif_mab_spleen <- 0.021727679283426147; clcnv_mab_spleen <- 4.405175053256755
    cldif_mab_other <- 0.018911952912469773; clcnv_mab_other <- 3.8343010356346436

    cldif_alb_lung <- 57.93932566335267; clcnv_alb_lung <- 150.26571662501527
    cldif_alb_heart <- 2.4690134984432657; clcnv_alb_heart <- 6.403389726281894
    cldif_alb_muscle <- 10.659882969478543; clcnv_alb_muscle <- 27.646420375248802
    cldif_alb_skin <- 3.7028832472783026; clcnv_alb_skin <- 9.60343252808995
    cldif_alb_adipose <- 3.5777126713123337; clcnv_alb_adipose <- 9.278802476177052
    cldif_alb_bone <- 0.8252340008341721; clcnv_alb_bone <- 2.14024545675908
    cldif_alb_brain <- 6.8327846468141615; clcnv_alb_brain <- 17.720835887245286
    cldif_alb_kidney <- 11.594044036420508; clcnv_alb_kidney <- 30.069168319932082
    cldif_alb_liver <- 15.21618748469786; clcnv_alb_liver <- 39.463202074078225
    cldif_alb_small_intestine <- 3.9392103907051492; clcnv_alb_small_intestine <- 10.216347282592164
    cldif_alb_large_intestine <- 4.098141987160669; clcnv_alb_large_intestine <- 10.628536585148238
    cldif_alb_pancreas <- 0.9733365907175723; clcnv_alb_pancreas <- 2.5243497166560203
    cldif_alb_thymus <- 0.11243056823406514; clcnv_alb_thymus <- 0.29158882525509655
    cldif_alb_spleen <- 2.0202467260868984; clcnv_alb_spleen <- 5.2395125172608426
    cldif_alb_other <- 1.7584395671962432; clcnv_alb_other <- 4.560515309443026

    cldif_igg_lung <- 0.6231353179072838; clcnv_igg_lung <- 126.3374758731036
    cldif_igg_heart <- 0.026554149425369638; clcnv_igg_heart <- 5.383717012903415
    cldif_igg_muscle <- 0.11464664952498664; clcnv_igg_muscle <- 23.24401763478643
    cldif_igg_skin <- 0.03982437322230997; clcnv_igg_skin <- 8.074186531477698
    cldif_igg_adipose <- 0.03847816827853163; clcnv_igg_adipose <- 7.801250413563476
    cldif_igg_bone <- 0.00887536134689534; clcnv_igg_bone <- 1.7994337952054633
    cldif_igg_brain <- 0.07348634773251482; clcnv_igg_brain <- 14.89897846721065
    cldif_igg_kidney <- 0.1246935174641777; clcnv_igg_kidney <- 25.28096835703175
    cldif_igg_liver <- 0.16364953711588245; clcnv_igg_liver <- 33.179100675045696
    cldif_igg_small_intestine <- 0.04236606296348964; clcnv_igg_small_intestine <- 8.589501034002765
    cldif_igg_large_intestine <- 0.044075366441722254; clcnv_igg_large_intestine <- 8.936053509420569
    cldif_igg_pancreas <- 0.010468199257472862; clcnv_igg_pancreas <- 2.1223734767070224
    cldif_igg_thymus <- 0.001209186628890026; clcnv_igg_thymus <- 0.24515636036569988
    cldif_igg_spleen <- 0.021727679283426147; clcnv_igg_spleen <- 4.405175053256755
    cldif_igg_other <- 0.018911952912469773; clcnv_igg_other <- 3.8343010356346436

    # ================= Concentrations (uM) =================
    cv_mab_lung <- vp_lung / vv_lung; ci_mab_lung <- is_lung / vi_lung; ceu_mab_lung <- eu_lung / ve_lung; ceb_mab_lung <- eb_lung / ve_lung
    cv_mab_heart <- vp_heart / vv_heart; ci_mab_heart <- is_heart / vi_heart; ceu_mab_heart <- eu_heart / ve_heart; ceb_mab_heart <- eb_heart / ve_heart
    cv_mab_muscle <- vp_muscle / vv_muscle; ci_mab_muscle <- is_muscle / vi_muscle; ceu_mab_muscle <- eu_muscle / ve_muscle; ceb_mab_muscle <- eb_muscle / ve_muscle
    cv_mab_skin <- vp_skin / vv_skin; ci_mab_skin <- is_skin / vi_skin; ceu_mab_skin <- eu_skin / ve_skin; ceb_mab_skin <- eb_skin / ve_skin
    cv_mab_adipose <- vp_adipose / vv_adipose; ci_mab_adipose <- is_adipose / vi_adipose; ceu_mab_adipose <- eu_adipose / ve_adipose; ceb_mab_adipose <- eb_adipose / ve_adipose
    cv_mab_bone <- vp_bone / vv_bone; ci_mab_bone <- is_bone / vi_bone; ceu_mab_bone <- eu_bone / ve_bone; ceb_mab_bone <- eb_bone / ve_bone
    cv_mab_brain <- vp_brain / vv_brain; ci_mab_brain <- is_brain / vi_brain; ceu_mab_brain <- eu_brain / ve_brain; ceb_mab_brain <- eb_brain / ve_brain
    cv_mab_kidney <- vp_kidney / vv_kidney; ci_mab_kidney <- is_kidney / vi_kidney; ceu_mab_kidney <- eu_kidney / ve_kidney; ceb_mab_kidney <- eb_kidney / ve_kidney
    cv_mab_liver <- vp_liver / vv_liver; ci_mab_liver <- is_liver / vi_liver; ceu_mab_liver <- eu_liver / ve_liver; ceb_mab_liver <- eb_liver / ve_liver
    cv_mab_small_intestine <- vp_small_intestine / vv_small_intestine; ci_mab_small_intestine <- is_small_intestine / vi_small_intestine; ceu_mab_small_intestine <- eu_small_intestine / ve_small_intestine; ceb_mab_small_intestine <- eb_small_intestine / ve_small_intestine
    cv_mab_large_intestine <- vp_large_intestine / vv_large_intestine; ci_mab_large_intestine <- is_large_intestine / vi_large_intestine; ceu_mab_large_intestine <- eu_large_intestine / ve_large_intestine; ceb_mab_large_intestine <- eb_large_intestine / ve_large_intestine
    cv_mab_pancreas <- vp_pancreas / vv_pancreas; ci_mab_pancreas <- is_pancreas / vi_pancreas; ceu_mab_pancreas <- eu_pancreas / ve_pancreas; ceb_mab_pancreas <- eb_pancreas / ve_pancreas
    cv_mab_thymus <- vp_thymus / vv_thymus; ci_mab_thymus <- is_thymus / vi_thymus; ceu_mab_thymus <- eu_thymus / ve_thymus; ceb_mab_thymus <- eb_thymus / ve_thymus
    cv_mab_spleen <- vp_spleen / vv_spleen; ci_mab_spleen <- is_spleen / vi_spleen; ceu_mab_spleen <- eu_spleen / ve_spleen; ceb_mab_spleen <- eb_spleen / ve_spleen
    cv_mab_other <- vp_other / vv_other; ci_mab_other <- is_other / vi_other; ceu_mab_other <- eu_other / ve_other; ceb_mab_other <- eb_other / ve_other
    cpl_mab <- plasma / vpl; cln_mab <- lnode / vln; cbl_mab <- urine / vbl

    cv_alb_lung <- vp_lung_alb / vv_lung; ci_alb_lung <- is_lung_alb / vi_lung; ceu_alb_lung <- eu_lung_alb / ve_lung; ceb_alb_lung <- eb_lung_alb / ve_lung
    cv_alb_heart <- vp_heart_alb / vv_heart; ci_alb_heart <- is_heart_alb / vi_heart; ceu_alb_heart <- eu_heart_alb / ve_heart; ceb_alb_heart <- eb_heart_alb / ve_heart
    cv_alb_muscle <- vp_muscle_alb / vv_muscle; ci_alb_muscle <- is_muscle_alb / vi_muscle; ceu_alb_muscle <- eu_muscle_alb / ve_muscle; ceb_alb_muscle <- eb_muscle_alb / ve_muscle
    cv_alb_skin <- vp_skin_alb / vv_skin; ci_alb_skin <- is_skin_alb / vi_skin; ceu_alb_skin <- eu_skin_alb / ve_skin; ceb_alb_skin <- eb_skin_alb / ve_skin
    cv_alb_adipose <- vp_adipose_alb / vv_adipose; ci_alb_adipose <- is_adipose_alb / vi_adipose; ceu_alb_adipose <- eu_adipose_alb / ve_adipose; ceb_alb_adipose <- eb_adipose_alb / ve_adipose
    cv_alb_bone <- vp_bone_alb / vv_bone; ci_alb_bone <- is_bone_alb / vi_bone; ceu_alb_bone <- eu_bone_alb / ve_bone; ceb_alb_bone <- eb_bone_alb / ve_bone
    cv_alb_brain <- vp_brain_alb / vv_brain; ci_alb_brain <- is_brain_alb / vi_brain; ceu_alb_brain <- eu_brain_alb / ve_brain; ceb_alb_brain <- eb_brain_alb / ve_brain
    cv_alb_kidney <- vp_kidney_alb / vv_kidney; ci_alb_kidney <- is_kidney_alb / vi_kidney; ceu_alb_kidney <- eu_kidney_alb / ve_kidney; ceb_alb_kidney <- eb_kidney_alb / ve_kidney
    cv_alb_liver <- vp_liver_alb / vv_liver; ci_alb_liver <- is_liver_alb / vi_liver; ceu_alb_liver <- eu_liver_alb / ve_liver; ceb_alb_liver <- eb_liver_alb / ve_liver
    cv_alb_small_intestine <- vp_small_intestine_alb / vv_small_intestine; ci_alb_small_intestine <- is_small_intestine_alb / vi_small_intestine; ceu_alb_small_intestine <- eu_small_intestine_alb / ve_small_intestine; ceb_alb_small_intestine <- eb_small_intestine_alb / ve_small_intestine
    cv_alb_large_intestine <- vp_large_intestine_alb / vv_large_intestine; ci_alb_large_intestine <- is_large_intestine_alb / vi_large_intestine; ceu_alb_large_intestine <- eu_large_intestine_alb / ve_large_intestine; ceb_alb_large_intestine <- eb_large_intestine_alb / ve_large_intestine
    cv_alb_pancreas <- vp_pancreas_alb / vv_pancreas; ci_alb_pancreas <- is_pancreas_alb / vi_pancreas; ceu_alb_pancreas <- eu_pancreas_alb / ve_pancreas; ceb_alb_pancreas <- eb_pancreas_alb / ve_pancreas
    cv_alb_thymus <- vp_thymus_alb / vv_thymus; ci_alb_thymus <- is_thymus_alb / vi_thymus; ceu_alb_thymus <- eu_thymus_alb / ve_thymus; ceb_alb_thymus <- eb_thymus_alb / ve_thymus
    cv_alb_spleen <- vp_spleen_alb / vv_spleen; ci_alb_spleen <- is_spleen_alb / vi_spleen; ceu_alb_spleen <- eu_spleen_alb / ve_spleen; ceb_alb_spleen <- eb_spleen_alb / ve_spleen
    cv_alb_other <- vp_other_alb / vv_other; ci_alb_other <- is_other_alb / vi_other; ceu_alb_other <- eu_other_alb / ve_other; ceb_alb_other <- eb_other_alb / ve_other
    cpl_alb <- plasma_alb / vpl; cln_alb <- lnode_alb / vln; cbl_alb <- urine_alb / vbl

    cv_igg_lung <- vp_lung_igg / vv_lung; ci_igg_lung <- is_lung_igg / vi_lung; ceu_igg_lung <- eu_lung_igg / ve_lung; ceb_igg_lung <- eb_lung_igg / ve_lung
    cv_igg_heart <- vp_heart_igg / vv_heart; ci_igg_heart <- is_heart_igg / vi_heart; ceu_igg_heart <- eu_heart_igg / ve_heart; ceb_igg_heart <- eb_heart_igg / ve_heart
    cv_igg_muscle <- vp_muscle_igg / vv_muscle; ci_igg_muscle <- is_muscle_igg / vi_muscle; ceu_igg_muscle <- eu_muscle_igg / ve_muscle; ceb_igg_muscle <- eb_muscle_igg / ve_muscle
    cv_igg_skin <- vp_skin_igg / vv_skin; ci_igg_skin <- is_skin_igg / vi_skin; ceu_igg_skin <- eu_skin_igg / ve_skin; ceb_igg_skin <- eb_skin_igg / ve_skin
    cv_igg_adipose <- vp_adipose_igg / vv_adipose; ci_igg_adipose <- is_adipose_igg / vi_adipose; ceu_igg_adipose <- eu_adipose_igg / ve_adipose; ceb_igg_adipose <- eb_adipose_igg / ve_adipose
    cv_igg_bone <- vp_bone_igg / vv_bone; ci_igg_bone <- is_bone_igg / vi_bone; ceu_igg_bone <- eu_bone_igg / ve_bone; ceb_igg_bone <- eb_bone_igg / ve_bone
    cv_igg_brain <- vp_brain_igg / vv_brain; ci_igg_brain <- is_brain_igg / vi_brain; ceu_igg_brain <- eu_brain_igg / ve_brain; ceb_igg_brain <- eb_brain_igg / ve_brain
    cv_igg_kidney <- vp_kidney_igg / vv_kidney; ci_igg_kidney <- is_kidney_igg / vi_kidney; ceu_igg_kidney <- eu_kidney_igg / ve_kidney; ceb_igg_kidney <- eb_kidney_igg / ve_kidney
    cv_igg_liver <- vp_liver_igg / vv_liver; ci_igg_liver <- is_liver_igg / vi_liver; ceu_igg_liver <- eu_liver_igg / ve_liver; ceb_igg_liver <- eb_liver_igg / ve_liver
    cv_igg_small_intestine <- vp_small_intestine_igg / vv_small_intestine; ci_igg_small_intestine <- is_small_intestine_igg / vi_small_intestine; ceu_igg_small_intestine <- eu_small_intestine_igg / ve_small_intestine; ceb_igg_small_intestine <- eb_small_intestine_igg / ve_small_intestine
    cv_igg_large_intestine <- vp_large_intestine_igg / vv_large_intestine; ci_igg_large_intestine <- is_large_intestine_igg / vi_large_intestine; ceu_igg_large_intestine <- eu_large_intestine_igg / ve_large_intestine; ceb_igg_large_intestine <- eb_large_intestine_igg / ve_large_intestine
    cv_igg_pancreas <- vp_pancreas_igg / vv_pancreas; ci_igg_pancreas <- is_pancreas_igg / vi_pancreas; ceu_igg_pancreas <- eu_pancreas_igg / ve_pancreas; ceb_igg_pancreas <- eb_pancreas_igg / ve_pancreas
    cv_igg_thymus <- vp_thymus_igg / vv_thymus; ci_igg_thymus <- is_thymus_igg / vi_thymus; ceu_igg_thymus <- eu_thymus_igg / ve_thymus; ceb_igg_thymus <- eb_thymus_igg / ve_thymus
    cv_igg_spleen <- vp_spleen_igg / vv_spleen; ci_igg_spleen <- is_spleen_igg / vi_spleen; ceu_igg_spleen <- eu_spleen_igg / ve_spleen; ceb_igg_spleen <- eb_spleen_igg / ve_spleen
    cv_igg_other <- vp_other_igg / vv_other; ci_igg_other <- is_other_igg / vi_other; ceu_igg_other <- eu_other_igg / ve_other; ceb_igg_other <- eb_other_igg / ve_other
    cpl_igg <- plasma_igg / vpl; cln_igg <- lnode_igg / vln; cbl_igg <- urine_igg / vbl

    cfr_igg_lung <- fr_lung_igg / ve_lung
    cfr_igg_heart <- fr_heart_igg / ve_heart
    cfr_igg_muscle <- fr_muscle_igg / ve_muscle
    cfr_igg_skin <- fr_skin_igg / ve_skin
    cfr_igg_adipose <- fr_adipose_igg / ve_adipose
    cfr_igg_bone <- fr_bone_igg / ve_bone
    cfr_igg_brain <- fr_brain_igg / ve_brain
    cfr_igg_kidney <- fr_kidney_igg / ve_kidney
    cfr_igg_liver <- fr_liver_igg / ve_liver
    cfr_igg_small_intestine <- fr_small_intestine_igg / ve_small_intestine
    cfr_igg_large_intestine <- fr_large_intestine_igg / ve_large_intestine
    cfr_igg_pancreas <- fr_pancreas_igg / ve_pancreas
    cfr_igg_thymus <- fr_thymus_igg / ve_thymus
    cfr_igg_spleen <- fr_spleen_igg / ve_spleen
    cfr_igg_other <- fr_other_igg / ve_other

    cfr_alb_lung <- fr_lung_alb / ve_lung
    cfr_alb_heart <- fr_heart_alb / ve_heart
    cfr_alb_muscle <- fr_muscle_alb / ve_muscle
    cfr_alb_skin <- fr_skin_alb / ve_skin
    cfr_alb_adipose <- fr_adipose_alb / ve_adipose
    cfr_alb_bone <- fr_bone_alb / ve_bone
    cfr_alb_brain <- fr_brain_alb / ve_brain
    cfr_alb_kidney <- fr_kidney_alb / ve_kidney
    cfr_alb_liver <- fr_liver_alb / ve_liver
    cfr_alb_small_intestine <- fr_small_intestine_alb / ve_small_intestine
    cfr_alb_large_intestine <- fr_large_intestine_alb / ve_large_intestine
    cfr_alb_pancreas <- fr_pancreas_alb / ve_pancreas
    cfr_alb_thymus <- fr_thymus_alb / ve_thymus
    cfr_alb_spleen <- fr_spleen_alb / ve_spleen
    cfr_alb_other <- fr_other_alb / ve_other

    # ================= Two-pore extravasation flux (nmol/h) =================
    # Appendix S3 reactions <SUB>_<ORG>_is_convection_LP / _SP:
    #   rate = (C_vascular - C_interstitial) * CL_TP
    tp_mab_lung <- (cv_mab_lung - ci_mab_lung) * (cldif_mab_lung + clcnv_mab_lung)
    tp_mab_heart <- (cv_mab_heart - ci_mab_heart) * (cldif_mab_heart + clcnv_mab_heart)
    tp_mab_muscle <- (cv_mab_muscle - ci_mab_muscle) * (cldif_mab_muscle + clcnv_mab_muscle)
    tp_mab_skin <- (cv_mab_skin - ci_mab_skin) * (cldif_mab_skin + clcnv_mab_skin)
    tp_mab_adipose <- (cv_mab_adipose - ci_mab_adipose) * (cldif_mab_adipose + clcnv_mab_adipose)
    tp_mab_bone <- (cv_mab_bone - ci_mab_bone) * (cldif_mab_bone + clcnv_mab_bone)
    tp_mab_brain <- (cv_mab_brain - ci_mab_brain) * (cldif_mab_brain + clcnv_mab_brain)
    tp_mab_kidney <- (cv_mab_kidney - ci_mab_kidney) * (cldif_mab_kidney + clcnv_mab_kidney)
    tp_mab_liver <- (cv_mab_liver - ci_mab_liver) * (cldif_mab_liver + clcnv_mab_liver)
    tp_mab_small_intestine <- (cv_mab_small_intestine - ci_mab_small_intestine) * (cldif_mab_small_intestine + clcnv_mab_small_intestine)
    tp_mab_large_intestine <- (cv_mab_large_intestine - ci_mab_large_intestine) * (cldif_mab_large_intestine + clcnv_mab_large_intestine)
    tp_mab_pancreas <- (cv_mab_pancreas - ci_mab_pancreas) * (cldif_mab_pancreas + clcnv_mab_pancreas)
    tp_mab_thymus <- (cv_mab_thymus - ci_mab_thymus) * (cldif_mab_thymus + clcnv_mab_thymus)
    tp_mab_spleen <- (cv_mab_spleen - ci_mab_spleen) * (cldif_mab_spleen + clcnv_mab_spleen)
    tp_mab_other <- (cv_mab_other - ci_mab_other) * (cldif_mab_other + clcnv_mab_other)

    tp_alb_lung <- (cv_alb_lung - ci_alb_lung) * (cldif_alb_lung + clcnv_alb_lung)
    tp_alb_heart <- (cv_alb_heart - ci_alb_heart) * (cldif_alb_heart + clcnv_alb_heart)
    tp_alb_muscle <- (cv_alb_muscle - ci_alb_muscle) * (cldif_alb_muscle + clcnv_alb_muscle)
    tp_alb_skin <- (cv_alb_skin - ci_alb_skin) * (cldif_alb_skin + clcnv_alb_skin)
    tp_alb_adipose <- (cv_alb_adipose - ci_alb_adipose) * (cldif_alb_adipose + clcnv_alb_adipose)
    tp_alb_bone <- (cv_alb_bone - ci_alb_bone) * (cldif_alb_bone + clcnv_alb_bone)
    tp_alb_brain <- (cv_alb_brain - ci_alb_brain) * (cldif_alb_brain + clcnv_alb_brain)
    tp_alb_kidney <- (cv_alb_kidney - ci_alb_kidney) * (cldif_alb_kidney + clcnv_alb_kidney)
    tp_alb_liver <- (cv_alb_liver - ci_alb_liver) * (cldif_alb_liver + clcnv_alb_liver)
    tp_alb_small_intestine <- (cv_alb_small_intestine - ci_alb_small_intestine) * (cldif_alb_small_intestine + clcnv_alb_small_intestine)
    tp_alb_large_intestine <- (cv_alb_large_intestine - ci_alb_large_intestine) * (cldif_alb_large_intestine + clcnv_alb_large_intestine)
    tp_alb_pancreas <- (cv_alb_pancreas - ci_alb_pancreas) * (cldif_alb_pancreas + clcnv_alb_pancreas)
    tp_alb_thymus <- (cv_alb_thymus - ci_alb_thymus) * (cldif_alb_thymus + clcnv_alb_thymus)
    tp_alb_spleen <- (cv_alb_spleen - ci_alb_spleen) * (cldif_alb_spleen + clcnv_alb_spleen)
    tp_alb_other <- (cv_alb_other - ci_alb_other) * (cldif_alb_other + clcnv_alb_other)

    tp_igg_lung <- (cv_igg_lung - ci_igg_lung) * (cldif_igg_lung + clcnv_igg_lung)
    tp_igg_heart <- (cv_igg_heart - ci_igg_heart) * (cldif_igg_heart + clcnv_igg_heart)
    tp_igg_muscle <- (cv_igg_muscle - ci_igg_muscle) * (cldif_igg_muscle + clcnv_igg_muscle)
    tp_igg_skin <- (cv_igg_skin - ci_igg_skin) * (cldif_igg_skin + clcnv_igg_skin)
    tp_igg_adipose <- (cv_igg_adipose - ci_igg_adipose) * (cldif_igg_adipose + clcnv_igg_adipose)
    tp_igg_bone <- (cv_igg_bone - ci_igg_bone) * (cldif_igg_bone + clcnv_igg_bone)
    tp_igg_brain <- (cv_igg_brain - ci_igg_brain) * (cldif_igg_brain + clcnv_igg_brain)
    tp_igg_kidney <- (cv_igg_kidney - ci_igg_kidney) * (cldif_igg_kidney + clcnv_igg_kidney)
    tp_igg_liver <- (cv_igg_liver - ci_igg_liver) * (cldif_igg_liver + clcnv_igg_liver)
    tp_igg_small_intestine <- (cv_igg_small_intestine - ci_igg_small_intestine) * (cldif_igg_small_intestine + clcnv_igg_small_intestine)
    tp_igg_large_intestine <- (cv_igg_large_intestine - ci_igg_large_intestine) * (cldif_igg_large_intestine + clcnv_igg_large_intestine)
    tp_igg_pancreas <- (cv_igg_pancreas - ci_igg_pancreas) * (cldif_igg_pancreas + clcnv_igg_pancreas)
    tp_igg_thymus <- (cv_igg_thymus - ci_igg_thymus) * (cldif_igg_thymus + clcnv_igg_thymus)
    tp_igg_spleen <- (cv_igg_spleen - ci_igg_spleen) * (cldif_igg_spleen + clcnv_igg_spleen)
    tp_igg_other <- (cv_igg_other - ci_igg_other) * (cldif_igg_other + clcnv_igg_other)

    # ============================================================
    # durvalumab (exogenous IgG1)
    # ============================================================

    # --- Vascular plasma of each organ (Eq S3) ---
    d/dt(vp_lung) <- (q_lung + j_lung) * cpl_mab - sumq_sys * cv_mab_lung -
      spino_mab * clup_lung * cv_mab_lung + 2 * clup_lung * fr * ceb_mab_lung - tp_mab_lung
    d/dt(vp_heart) <- q_heart * cv_mab_lung - (q_heart - j_heart) * cv_mab_heart -
      spino_mab * clup_heart * cv_mab_heart + 2 * clup_heart * fr * ceb_mab_heart - tp_mab_heart
    d/dt(vp_muscle) <- q_muscle * cv_mab_lung - (q_muscle - j_muscle) * cv_mab_muscle -
      spino_mab * clup_muscle * cv_mab_muscle + 2 * clup_muscle * fr * ceb_mab_muscle - tp_mab_muscle
    d/dt(vp_skin) <- q_skin * cv_mab_lung - (q_skin - j_skin) * cv_mab_skin -
      spino_mab * clup_skin * cv_mab_skin + 2 * clup_skin * fr * ceb_mab_skin - tp_mab_skin
    d/dt(vp_adipose) <- q_adipose * cv_mab_lung - (q_adipose - j_adipose) * cv_mab_adipose -
      spino_mab * clup_adipose * cv_mab_adipose + 2 * clup_adipose * fr * ceb_mab_adipose - tp_mab_adipose
    d/dt(vp_bone) <- q_bone * cv_mab_lung - (q_bone - j_bone) * cv_mab_bone -
      spino_mab * clup_bone * cv_mab_bone + 2 * clup_bone * fr * ceb_mab_bone - tp_mab_bone
    d/dt(vp_brain) <- q_brain * cv_mab_lung - (q_brain - j_brain) * cv_mab_brain -
      spino_mab * clup_brain * cv_mab_brain + 2 * clup_brain * fr * ceb_mab_brain - tp_mab_brain
    d/dt(vp_kidney) <- q_kidney * cv_mab_lung - (q_kidney - j_kidney) * cv_mab_kidney - clren_mab * cv_mab_kidney -
      spino_mab * clup_kidney * cv_mab_kidney + 2 * clup_kidney * fr * ceb_mab_kidney - tp_mab_kidney
    d/dt(vp_liver) <- q_liver * cv_mab_lung + (q_spleen - j_spleen) * cv_mab_spleen + (q_pancreas - j_pancreas) * cv_mab_pancreas + (q_small_intestine - j_small_intestine) * cv_mab_small_intestine + (q_large_intestine - j_large_intestine) * cv_mab_large_intestine - qout_liver * cv_mab_liver -
      spino_mab * clup_liver * cv_mab_liver + 2 * clup_liver * fr * ceb_mab_liver - tp_mab_liver
    d/dt(vp_small_intestine) <- q_small_intestine * cv_mab_lung - (q_small_intestine - j_small_intestine) * cv_mab_small_intestine -
      spino_mab * clup_small_intestine * cv_mab_small_intestine + 2 * clup_small_intestine * fr * ceb_mab_small_intestine - tp_mab_small_intestine
    d/dt(vp_large_intestine) <- q_large_intestine * cv_mab_lung - (q_large_intestine - j_large_intestine) * cv_mab_large_intestine -
      spino_mab * clup_large_intestine * cv_mab_large_intestine + 2 * clup_large_intestine * fr * ceb_mab_large_intestine - tp_mab_large_intestine
    d/dt(vp_pancreas) <- q_pancreas * cv_mab_lung - (q_pancreas - j_pancreas) * cv_mab_pancreas -
      spino_mab * clup_pancreas * cv_mab_pancreas + 2 * clup_pancreas * fr * ceb_mab_pancreas - tp_mab_pancreas
    d/dt(vp_thymus) <- q_thymus * cv_mab_lung - (q_thymus - j_thymus) * cv_mab_thymus -
      spino_mab * clup_thymus * cv_mab_thymus + 2 * clup_thymus * fr * ceb_mab_thymus - tp_mab_thymus
    d/dt(vp_spleen) <- q_spleen * cv_mab_lung - (q_spleen - j_spleen) * cv_mab_spleen -
      spino_mab * clup_spleen * cv_mab_spleen + 2 * clup_spleen * fr * ceb_mab_spleen - tp_mab_spleen
    d/dt(vp_other) <- q_other * cv_mab_lung - (q_other - j_other) * cv_mab_other -
      spino_mab * clup_other * cv_mab_other + 2 * clup_other * fr * ceb_mab_other - tp_mab_other

    # --- Interstitial space of each organ (Eq S4) ---
    d/dt(is_lung) <- tp_mab_lung - (1 - sigis) * j_lung * ci_mab_lung -
      spino_mab * clup_lung * ci_mab_lung + 2 * clup_lung * (1 - fr) * ceb_mab_lung
    d/dt(is_heart) <- tp_mab_heart - (1 - sigis) * j_heart * ci_mab_heart -
      spino_mab * clup_heart * ci_mab_heart + 2 * clup_heart * (1 - fr) * ceb_mab_heart
    d/dt(is_muscle) <- tp_mab_muscle - (1 - sigis) * j_muscle * ci_mab_muscle -
      spino_mab * clup_muscle * ci_mab_muscle + 2 * clup_muscle * (1 - fr) * ceb_mab_muscle
    d/dt(is_skin) <- tp_mab_skin - (1 - sigis) * j_skin * ci_mab_skin -
      spino_mab * clup_skin * ci_mab_skin + 2 * clup_skin * (1 - fr) * ceb_mab_skin
    d/dt(is_adipose) <- tp_mab_adipose - (1 - sigis) * j_adipose * ci_mab_adipose -
      spino_mab * clup_adipose * ci_mab_adipose + 2 * clup_adipose * (1 - fr) * ceb_mab_adipose
    d/dt(is_bone) <- tp_mab_bone - (1 - sigis) * j_bone * ci_mab_bone -
      spino_mab * clup_bone * ci_mab_bone + 2 * clup_bone * (1 - fr) * ceb_mab_bone
    d/dt(is_brain) <- tp_mab_brain - (1 - sigis) * j_brain * ci_mab_brain -
      spino_mab * clup_brain * ci_mab_brain + 2 * clup_brain * (1 - fr) * ceb_mab_brain
    d/dt(is_kidney) <- tp_mab_kidney - (1 - sigis) * j_kidney * ci_mab_kidney -
      spino_mab * clup_kidney * ci_mab_kidney + 2 * clup_kidney * (1 - fr) * ceb_mab_kidney
    d/dt(is_liver) <- tp_mab_liver - (1 - sigis) * j_liver * ci_mab_liver -
      spino_mab * clup_liver * ci_mab_liver + 2 * clup_liver * (1 - fr) * ceb_mab_liver
    d/dt(is_small_intestine) <- tp_mab_small_intestine - (1 - sigis) * j_small_intestine * ci_mab_small_intestine -
      spino_mab * clup_small_intestine * ci_mab_small_intestine + 2 * clup_small_intestine * (1 - fr) * ceb_mab_small_intestine
    d/dt(is_large_intestine) <- tp_mab_large_intestine - (1 - sigis) * j_large_intestine * ci_mab_large_intestine -
      spino_mab * clup_large_intestine * ci_mab_large_intestine + 2 * clup_large_intestine * (1 - fr) * ceb_mab_large_intestine
    d/dt(is_pancreas) <- tp_mab_pancreas - (1 - sigis) * j_pancreas * ci_mab_pancreas -
      spino_mab * clup_pancreas * ci_mab_pancreas + 2 * clup_pancreas * (1 - fr) * ceb_mab_pancreas
    d/dt(is_thymus) <- tp_mab_thymus - (1 - sigis) * j_thymus * ci_mab_thymus -
      spino_mab * clup_thymus * ci_mab_thymus + 2 * clup_thymus * (1 - fr) * ceb_mab_thymus
    d/dt(is_spleen) <- tp_mab_spleen - (1 - sigis) * j_spleen * ci_mab_spleen -
      spino_mab * clup_spleen * ci_mab_spleen + 2 * clup_spleen * (1 - fr) * ceb_mab_spleen
    d/dt(is_other) <- tp_mab_other - (1 - sigis) * j_other * ci_mab_other -
      spino_mab * clup_other * ci_mab_other + 2 * clup_other * (1 - fr) * ceb_mab_other

    # --- Endosomal unbound / FcRn-bound (Eq S14, S15, S17, S18) ---
    bind_mab_lung <- kon_mab * ceu_mab_lung * cfr_igg_lung * ve_lung
    d/dt(eu_lung) <- spino_mab * clup_lung * (cv_mab_lung + ci_mab_lung) -
      kdeg * eu_lung - bind_mab_lung + koff_mab * eb_lung
    d/dt(eb_lung) <- bind_mab_lung - koff_mab * eb_lung - 2 * clup_lung * ceb_mab_lung
    bind_mab_heart <- kon_mab * ceu_mab_heart * cfr_igg_heart * ve_heart
    d/dt(eu_heart) <- spino_mab * clup_heart * (cv_mab_heart + ci_mab_heart) -
      kdeg * eu_heart - bind_mab_heart + koff_mab * eb_heart
    d/dt(eb_heart) <- bind_mab_heart - koff_mab * eb_heart - 2 * clup_heart * ceb_mab_heart
    bind_mab_muscle <- kon_mab * ceu_mab_muscle * cfr_igg_muscle * ve_muscle
    d/dt(eu_muscle) <- spino_mab * clup_muscle * (cv_mab_muscle + ci_mab_muscle) -
      kdeg * eu_muscle - bind_mab_muscle + koff_mab * eb_muscle
    d/dt(eb_muscle) <- bind_mab_muscle - koff_mab * eb_muscle - 2 * clup_muscle * ceb_mab_muscle
    bind_mab_skin <- kon_mab * ceu_mab_skin * cfr_igg_skin * ve_skin
    d/dt(eu_skin) <- spino_mab * clup_skin * (cv_mab_skin + ci_mab_skin) -
      kdeg * eu_skin - bind_mab_skin + koff_mab * eb_skin
    d/dt(eb_skin) <- bind_mab_skin - koff_mab * eb_skin - 2 * clup_skin * ceb_mab_skin
    bind_mab_adipose <- kon_mab * ceu_mab_adipose * cfr_igg_adipose * ve_adipose
    d/dt(eu_adipose) <- spino_mab * clup_adipose * (cv_mab_adipose + ci_mab_adipose) -
      kdeg * eu_adipose - bind_mab_adipose + koff_mab * eb_adipose
    d/dt(eb_adipose) <- bind_mab_adipose - koff_mab * eb_adipose - 2 * clup_adipose * ceb_mab_adipose
    bind_mab_bone <- kon_mab * ceu_mab_bone * cfr_igg_bone * ve_bone
    d/dt(eu_bone) <- spino_mab * clup_bone * (cv_mab_bone + ci_mab_bone) -
      kdeg * eu_bone - bind_mab_bone + koff_mab * eb_bone
    d/dt(eb_bone) <- bind_mab_bone - koff_mab * eb_bone - 2 * clup_bone * ceb_mab_bone
    bind_mab_brain <- kon_mab * ceu_mab_brain * cfr_igg_brain * ve_brain
    d/dt(eu_brain) <- spino_mab * clup_brain * (cv_mab_brain + ci_mab_brain) -
      kdeg * eu_brain - bind_mab_brain + koff_mab * eb_brain
    d/dt(eb_brain) <- bind_mab_brain - koff_mab * eb_brain - 2 * clup_brain * ceb_mab_brain
    bind_mab_kidney <- kon_mab * ceu_mab_kidney * cfr_igg_kidney * ve_kidney
    d/dt(eu_kidney) <- spino_mab * clup_kidney * (cv_mab_kidney + ci_mab_kidney) -
      kdeg * eu_kidney - bind_mab_kidney + koff_mab * eb_kidney
    d/dt(eb_kidney) <- bind_mab_kidney - koff_mab * eb_kidney - 2 * clup_kidney * ceb_mab_kidney
    bind_mab_liver <- kon_mab * ceu_mab_liver * cfr_igg_liver * ve_liver
    d/dt(eu_liver) <- spino_mab * clup_liver * (cv_mab_liver + ci_mab_liver) -
      kdeg * eu_liver - bind_mab_liver + koff_mab * eb_liver
    d/dt(eb_liver) <- bind_mab_liver - koff_mab * eb_liver - 2 * clup_liver * ceb_mab_liver
    bind_mab_small_intestine <- kon_mab * ceu_mab_small_intestine * cfr_igg_small_intestine * ve_small_intestine
    d/dt(eu_small_intestine) <- spino_mab * clup_small_intestine * (cv_mab_small_intestine + ci_mab_small_intestine) -
      kdeg * eu_small_intestine - bind_mab_small_intestine + koff_mab * eb_small_intestine
    d/dt(eb_small_intestine) <- bind_mab_small_intestine - koff_mab * eb_small_intestine - 2 * clup_small_intestine * ceb_mab_small_intestine
    bind_mab_large_intestine <- kon_mab * ceu_mab_large_intestine * cfr_igg_large_intestine * ve_large_intestine
    d/dt(eu_large_intestine) <- spino_mab * clup_large_intestine * (cv_mab_large_intestine + ci_mab_large_intestine) -
      kdeg * eu_large_intestine - bind_mab_large_intestine + koff_mab * eb_large_intestine
    d/dt(eb_large_intestine) <- bind_mab_large_intestine - koff_mab * eb_large_intestine - 2 * clup_large_intestine * ceb_mab_large_intestine
    bind_mab_pancreas <- kon_mab * ceu_mab_pancreas * cfr_igg_pancreas * ve_pancreas
    d/dt(eu_pancreas) <- spino_mab * clup_pancreas * (cv_mab_pancreas + ci_mab_pancreas) -
      kdeg * eu_pancreas - bind_mab_pancreas + koff_mab * eb_pancreas
    d/dt(eb_pancreas) <- bind_mab_pancreas - koff_mab * eb_pancreas - 2 * clup_pancreas * ceb_mab_pancreas
    bind_mab_thymus <- kon_mab * ceu_mab_thymus * cfr_igg_thymus * ve_thymus
    d/dt(eu_thymus) <- spino_mab * clup_thymus * (cv_mab_thymus + ci_mab_thymus) -
      kdeg * eu_thymus - bind_mab_thymus + koff_mab * eb_thymus
    d/dt(eb_thymus) <- bind_mab_thymus - koff_mab * eb_thymus - 2 * clup_thymus * ceb_mab_thymus
    bind_mab_spleen <- kon_mab * ceu_mab_spleen * cfr_igg_spleen * ve_spleen
    d/dt(eu_spleen) <- spino_mab * clup_spleen * (cv_mab_spleen + ci_mab_spleen) -
      kdeg * eu_spleen - bind_mab_spleen + koff_mab * eb_spleen
    d/dt(eb_spleen) <- bind_mab_spleen - koff_mab * eb_spleen - 2 * clup_spleen * ceb_mab_spleen
    bind_mab_other <- kon_mab * ceu_mab_other * cfr_igg_other * ve_other
    d/dt(eu_other) <- spino_mab * clup_other * (cv_mab_other + ci_mab_other) -
      kdeg * eu_other - bind_mab_other + koff_mab * eb_other
    d/dt(eb_other) <- bind_mab_other - koff_mab * eb_other - 2 * clup_other * ceb_mab_other

    # --- Central plasma (Eq S1) ---
    d/dt(plasma) <- (q_heart - j_heart) * cv_mab_heart +
      (q_muscle - j_muscle) * cv_mab_muscle +
      (q_skin - j_skin) * cv_mab_skin +
      (q_adipose - j_adipose) * cv_mab_adipose +
      (q_bone - j_bone) * cv_mab_bone +
      (q_brain - j_brain) * cv_mab_brain +
      (q_kidney - j_kidney) * cv_mab_kidney +
      (q_thymus - j_thymus) * cv_mab_thymus +
      (q_other - j_other) * cv_mab_other +
      qout_liver * cv_mab_liver +
      lnlf * cln_mab - (q_lung + j_lung) * cpl_mab

    # --- Lymph node (Eq S2) ---
    d/dt(lnode) <- (1 - sigis) * (
      j_lung * ci_mab_lung +
      j_heart * ci_mab_heart +
      j_muscle * ci_mab_muscle +
      j_skin * ci_mab_skin +
      j_adipose * ci_mab_adipose +
      j_bone * ci_mab_bone +
      j_brain * ci_mab_brain +
      j_kidney * ci_mab_kidney +
      j_liver * ci_mab_liver +
      j_small_intestine * ci_mab_small_intestine +
      j_large_intestine * ci_mab_large_intestine +
      j_pancreas * ci_mab_pancreas +
      j_thymus * ci_mab_thymus +
      j_spleen * ci_mab_spleen +
      j_other * ci_mab_other) - lnlf * cln_mab

    # --- Renal filtrate holding pool, cleared to urine (Eq S11) ---
    d/dt(urine) <- clren_mab * cv_mab_kidney - qurine * cbl_mab

    # ============================================================
    # endogenous serum albumin
    # ============================================================

    # --- Vascular plasma of each organ (Eq S3) ---
    d/dt(vp_lung_alb) <- (q_lung + j_lung) * cpl_alb - sumq_sys * cv_alb_lung -
      spino_alb * clup_lung * cv_alb_lung + 2 * clup_lung * fr * ceb_alb_lung - tp_alb_lung
    d/dt(vp_heart_alb) <- q_heart * cv_alb_lung - (q_heart - j_heart) * cv_alb_heart -
      spino_alb * clup_heart * cv_alb_heart + 2 * clup_heart * fr * ceb_alb_heart - tp_alb_heart
    d/dt(vp_muscle_alb) <- q_muscle * cv_alb_lung - (q_muscle - j_muscle) * cv_alb_muscle -
      spino_alb * clup_muscle * cv_alb_muscle + 2 * clup_muscle * fr * ceb_alb_muscle - tp_alb_muscle
    d/dt(vp_skin_alb) <- q_skin * cv_alb_lung - (q_skin - j_skin) * cv_alb_skin -
      spino_alb * clup_skin * cv_alb_skin + 2 * clup_skin * fr * ceb_alb_skin - tp_alb_skin
    d/dt(vp_adipose_alb) <- q_adipose * cv_alb_lung - (q_adipose - j_adipose) * cv_alb_adipose -
      spino_alb * clup_adipose * cv_alb_adipose + 2 * clup_adipose * fr * ceb_alb_adipose - tp_alb_adipose
    d/dt(vp_bone_alb) <- q_bone * cv_alb_lung - (q_bone - j_bone) * cv_alb_bone -
      spino_alb * clup_bone * cv_alb_bone + 2 * clup_bone * fr * ceb_alb_bone - tp_alb_bone
    d/dt(vp_brain_alb) <- q_brain * cv_alb_lung - (q_brain - j_brain) * cv_alb_brain -
      spino_alb * clup_brain * cv_alb_brain + 2 * clup_brain * fr * ceb_alb_brain - tp_alb_brain
    d/dt(vp_kidney_alb) <- q_kidney * cv_alb_lung - (q_kidney - j_kidney) * cv_alb_kidney - clren_alb * cv_alb_kidney -
      spino_alb * clup_kidney * cv_alb_kidney + 2 * clup_kidney * fr * ceb_alb_kidney - tp_alb_kidney
    d/dt(vp_liver_alb) <- q_liver * cv_alb_lung + (q_spleen - j_spleen) * cv_alb_spleen + (q_pancreas - j_pancreas) * cv_alb_pancreas + (q_small_intestine - j_small_intestine) * cv_alb_small_intestine + (q_large_intestine - j_large_intestine) * cv_alb_large_intestine - qout_liver * cv_alb_liver -
      spino_alb * clup_liver * cv_alb_liver + 2 * clup_liver * fr * ceb_alb_liver - tp_alb_liver
    d/dt(vp_small_intestine_alb) <- q_small_intestine * cv_alb_lung - (q_small_intestine - j_small_intestine) * cv_alb_small_intestine -
      spino_alb * clup_small_intestine * cv_alb_small_intestine + 2 * clup_small_intestine * fr * ceb_alb_small_intestine - tp_alb_small_intestine
    d/dt(vp_large_intestine_alb) <- q_large_intestine * cv_alb_lung - (q_large_intestine - j_large_intestine) * cv_alb_large_intestine -
      spino_alb * clup_large_intestine * cv_alb_large_intestine + 2 * clup_large_intestine * fr * ceb_alb_large_intestine - tp_alb_large_intestine
    d/dt(vp_pancreas_alb) <- q_pancreas * cv_alb_lung - (q_pancreas - j_pancreas) * cv_alb_pancreas -
      spino_alb * clup_pancreas * cv_alb_pancreas + 2 * clup_pancreas * fr * ceb_alb_pancreas - tp_alb_pancreas
    d/dt(vp_thymus_alb) <- q_thymus * cv_alb_lung - (q_thymus - j_thymus) * cv_alb_thymus -
      spino_alb * clup_thymus * cv_alb_thymus + 2 * clup_thymus * fr * ceb_alb_thymus - tp_alb_thymus
    d/dt(vp_spleen_alb) <- q_spleen * cv_alb_lung - (q_spleen - j_spleen) * cv_alb_spleen -
      spino_alb * clup_spleen * cv_alb_spleen + 2 * clup_spleen * fr * ceb_alb_spleen - tp_alb_spleen
    d/dt(vp_other_alb) <- q_other * cv_alb_lung - (q_other - j_other) * cv_alb_other -
      spino_alb * clup_other * cv_alb_other + 2 * clup_other * fr * ceb_alb_other - tp_alb_other

    # --- Interstitial space of each organ (Eq S4) ---
    d/dt(is_lung_alb) <- tp_alb_lung - (1 - sigis) * j_lung * ci_alb_lung -
      spino_alb * clup_lung * ci_alb_lung + 2 * clup_lung * (1 - fr) * ceb_alb_lung
    d/dt(is_heart_alb) <- tp_alb_heart - (1 - sigis) * j_heart * ci_alb_heart -
      spino_alb * clup_heart * ci_alb_heart + 2 * clup_heart * (1 - fr) * ceb_alb_heart
    d/dt(is_muscle_alb) <- tp_alb_muscle - (1 - sigis) * j_muscle * ci_alb_muscle -
      spino_alb * clup_muscle * ci_alb_muscle + 2 * clup_muscle * (1 - fr) * ceb_alb_muscle
    d/dt(is_skin_alb) <- tp_alb_skin - (1 - sigis) * j_skin * ci_alb_skin -
      spino_alb * clup_skin * ci_alb_skin + 2 * clup_skin * (1 - fr) * ceb_alb_skin
    d/dt(is_adipose_alb) <- tp_alb_adipose - (1 - sigis) * j_adipose * ci_alb_adipose -
      spino_alb * clup_adipose * ci_alb_adipose + 2 * clup_adipose * (1 - fr) * ceb_alb_adipose
    d/dt(is_bone_alb) <- tp_alb_bone - (1 - sigis) * j_bone * ci_alb_bone -
      spino_alb * clup_bone * ci_alb_bone + 2 * clup_bone * (1 - fr) * ceb_alb_bone
    d/dt(is_brain_alb) <- tp_alb_brain - (1 - sigis) * j_brain * ci_alb_brain -
      spino_alb * clup_brain * ci_alb_brain + 2 * clup_brain * (1 - fr) * ceb_alb_brain
    d/dt(is_kidney_alb) <- tp_alb_kidney - (1 - sigis) * j_kidney * ci_alb_kidney -
      spino_alb * clup_kidney * ci_alb_kidney + 2 * clup_kidney * (1 - fr) * ceb_alb_kidney
    d/dt(is_liver_alb) <- tp_alb_liver - (1 - sigis) * j_liver * ci_alb_liver -
      spino_alb * clup_liver * ci_alb_liver + 2 * clup_liver * (1 - fr) * ceb_alb_liver
    d/dt(is_small_intestine_alb) <- tp_alb_small_intestine - (1 - sigis) * j_small_intestine * ci_alb_small_intestine -
      spino_alb * clup_small_intestine * ci_alb_small_intestine + 2 * clup_small_intestine * (1 - fr) * ceb_alb_small_intestine
    d/dt(is_large_intestine_alb) <- tp_alb_large_intestine - (1 - sigis) * j_large_intestine * ci_alb_large_intestine -
      spino_alb * clup_large_intestine * ci_alb_large_intestine + 2 * clup_large_intestine * (1 - fr) * ceb_alb_large_intestine
    d/dt(is_pancreas_alb) <- tp_alb_pancreas - (1 - sigis) * j_pancreas * ci_alb_pancreas -
      spino_alb * clup_pancreas * ci_alb_pancreas + 2 * clup_pancreas * (1 - fr) * ceb_alb_pancreas
    d/dt(is_thymus_alb) <- tp_alb_thymus - (1 - sigis) * j_thymus * ci_alb_thymus -
      spino_alb * clup_thymus * ci_alb_thymus + 2 * clup_thymus * (1 - fr) * ceb_alb_thymus
    d/dt(is_spleen_alb) <- tp_alb_spleen - (1 - sigis) * j_spleen * ci_alb_spleen -
      spino_alb * clup_spleen * ci_alb_spleen + 2 * clup_spleen * (1 - fr) * ceb_alb_spleen
    d/dt(is_other_alb) <- tp_alb_other - (1 - sigis) * j_other * ci_alb_other -
      spino_alb * clup_other * ci_alb_other + 2 * clup_other * (1 - fr) * ceb_alb_other

    # --- Endosomal unbound / FcRn-bound (Eq S14, S15, S17, S18) ---
    bind_alb_lung <- kon_alb * ceu_alb_lung * cfr_alb_lung * ve_lung
    d/dt(eu_lung_alb) <- spino_alb * clup_lung * (cv_alb_lung + ci_alb_lung) -
      kdeg * eu_lung_alb - bind_alb_lung + koff_alb * eb_lung_alb
    d/dt(eb_lung_alb) <- bind_alb_lung - koff_alb * eb_lung_alb - 2 * clup_lung * ceb_alb_lung
    bind_alb_heart <- kon_alb * ceu_alb_heart * cfr_alb_heart * ve_heart
    d/dt(eu_heart_alb) <- spino_alb * clup_heart * (cv_alb_heart + ci_alb_heart) -
      kdeg * eu_heart_alb - bind_alb_heart + koff_alb * eb_heart_alb
    d/dt(eb_heart_alb) <- bind_alb_heart - koff_alb * eb_heart_alb - 2 * clup_heart * ceb_alb_heart
    bind_alb_muscle <- kon_alb * ceu_alb_muscle * cfr_alb_muscle * ve_muscle
    d/dt(eu_muscle_alb) <- spino_alb * clup_muscle * (cv_alb_muscle + ci_alb_muscle) -
      kdeg * eu_muscle_alb - bind_alb_muscle + koff_alb * eb_muscle_alb
    d/dt(eb_muscle_alb) <- bind_alb_muscle - koff_alb * eb_muscle_alb - 2 * clup_muscle * ceb_alb_muscle
    bind_alb_skin <- kon_alb * ceu_alb_skin * cfr_alb_skin * ve_skin
    d/dt(eu_skin_alb) <- spino_alb * clup_skin * (cv_alb_skin + ci_alb_skin) -
      kdeg * eu_skin_alb - bind_alb_skin + koff_alb * eb_skin_alb
    d/dt(eb_skin_alb) <- bind_alb_skin - koff_alb * eb_skin_alb - 2 * clup_skin * ceb_alb_skin
    bind_alb_adipose <- kon_alb * ceu_alb_adipose * cfr_alb_adipose * ve_adipose
    d/dt(eu_adipose_alb) <- spino_alb * clup_adipose * (cv_alb_adipose + ci_alb_adipose) -
      kdeg * eu_adipose_alb - bind_alb_adipose + koff_alb * eb_adipose_alb
    d/dt(eb_adipose_alb) <- bind_alb_adipose - koff_alb * eb_adipose_alb - 2 * clup_adipose * ceb_alb_adipose
    bind_alb_bone <- kon_alb * ceu_alb_bone * cfr_alb_bone * ve_bone
    d/dt(eu_bone_alb) <- spino_alb * clup_bone * (cv_alb_bone + ci_alb_bone) -
      kdeg * eu_bone_alb - bind_alb_bone + koff_alb * eb_bone_alb
    d/dt(eb_bone_alb) <- bind_alb_bone - koff_alb * eb_bone_alb - 2 * clup_bone * ceb_alb_bone
    bind_alb_brain <- kon_alb * ceu_alb_brain * cfr_alb_brain * ve_brain
    d/dt(eu_brain_alb) <- spino_alb * clup_brain * (cv_alb_brain + ci_alb_brain) -
      kdeg * eu_brain_alb - bind_alb_brain + koff_alb * eb_brain_alb
    d/dt(eb_brain_alb) <- bind_alb_brain - koff_alb * eb_brain_alb - 2 * clup_brain * ceb_alb_brain
    bind_alb_kidney <- kon_alb * ceu_alb_kidney * cfr_alb_kidney * ve_kidney
    d/dt(eu_kidney_alb) <- spino_alb * clup_kidney * (cv_alb_kidney + ci_alb_kidney) -
      kdeg * eu_kidney_alb - bind_alb_kidney + koff_alb * eb_kidney_alb
    d/dt(eb_kidney_alb) <- bind_alb_kidney - koff_alb * eb_kidney_alb - 2 * clup_kidney * ceb_alb_kidney
    bind_alb_liver <- kon_alb * ceu_alb_liver * cfr_alb_liver * ve_liver
    d/dt(eu_liver_alb) <- spino_alb * clup_liver * (cv_alb_liver + ci_alb_liver) -
      kdeg * eu_liver_alb - bind_alb_liver + koff_alb * eb_liver_alb
    d/dt(eb_liver_alb) <- bind_alb_liver - koff_alb * eb_liver_alb - 2 * clup_liver * ceb_alb_liver
    bind_alb_small_intestine <- kon_alb * ceu_alb_small_intestine * cfr_alb_small_intestine * ve_small_intestine
    d/dt(eu_small_intestine_alb) <- spino_alb * clup_small_intestine * (cv_alb_small_intestine + ci_alb_small_intestine) -
      kdeg * eu_small_intestine_alb - bind_alb_small_intestine + koff_alb * eb_small_intestine_alb
    d/dt(eb_small_intestine_alb) <- bind_alb_small_intestine - koff_alb * eb_small_intestine_alb - 2 * clup_small_intestine * ceb_alb_small_intestine
    bind_alb_large_intestine <- kon_alb * ceu_alb_large_intestine * cfr_alb_large_intestine * ve_large_intestine
    d/dt(eu_large_intestine_alb) <- spino_alb * clup_large_intestine * (cv_alb_large_intestine + ci_alb_large_intestine) -
      kdeg * eu_large_intestine_alb - bind_alb_large_intestine + koff_alb * eb_large_intestine_alb
    d/dt(eb_large_intestine_alb) <- bind_alb_large_intestine - koff_alb * eb_large_intestine_alb - 2 * clup_large_intestine * ceb_alb_large_intestine
    bind_alb_pancreas <- kon_alb * ceu_alb_pancreas * cfr_alb_pancreas * ve_pancreas
    d/dt(eu_pancreas_alb) <- spino_alb * clup_pancreas * (cv_alb_pancreas + ci_alb_pancreas) -
      kdeg * eu_pancreas_alb - bind_alb_pancreas + koff_alb * eb_pancreas_alb
    d/dt(eb_pancreas_alb) <- bind_alb_pancreas - koff_alb * eb_pancreas_alb - 2 * clup_pancreas * ceb_alb_pancreas
    bind_alb_thymus <- kon_alb * ceu_alb_thymus * cfr_alb_thymus * ve_thymus
    d/dt(eu_thymus_alb) <- spino_alb * clup_thymus * (cv_alb_thymus + ci_alb_thymus) -
      kdeg * eu_thymus_alb - bind_alb_thymus + koff_alb * eb_thymus_alb
    d/dt(eb_thymus_alb) <- bind_alb_thymus - koff_alb * eb_thymus_alb - 2 * clup_thymus * ceb_alb_thymus
    bind_alb_spleen <- kon_alb * ceu_alb_spleen * cfr_alb_spleen * ve_spleen
    d/dt(eu_spleen_alb) <- spino_alb * clup_spleen * (cv_alb_spleen + ci_alb_spleen) -
      kdeg * eu_spleen_alb - bind_alb_spleen + koff_alb * eb_spleen_alb
    d/dt(eb_spleen_alb) <- bind_alb_spleen - koff_alb * eb_spleen_alb - 2 * clup_spleen * ceb_alb_spleen
    bind_alb_other <- kon_alb * ceu_alb_other * cfr_alb_other * ve_other
    d/dt(eu_other_alb) <- spino_alb * clup_other * (cv_alb_other + ci_alb_other) -
      kdeg * eu_other_alb - bind_alb_other + koff_alb * eb_other_alb
    d/dt(eb_other_alb) <- bind_alb_other - koff_alb * eb_other_alb - 2 * clup_other * ceb_alb_other

    # --- Central plasma (Eq S1) ---
    d/dt(plasma_alb) <- (q_heart - j_heart) * cv_alb_heart +
      (q_muscle - j_muscle) * cv_alb_muscle +
      (q_skin - j_skin) * cv_alb_skin +
      (q_adipose - j_adipose) * cv_alb_adipose +
      (q_bone - j_bone) * cv_alb_bone +
      (q_brain - j_brain) * cv_alb_brain +
      (q_kidney - j_kidney) * cv_alb_kidney +
      (q_thymus - j_thymus) * cv_alb_thymus +
      (q_other - j_other) * cv_alb_other +
      qout_liver * cv_alb_liver +
      lnlf * cln_alb - (q_lung + j_lung) * cpl_alb +
      6.5156 * vpl   # zero-order synthesis ksyn * V_plasma

    # --- Lymph node (Eq S2) ---
    d/dt(lnode_alb) <- (1 - sigis) * (
      j_lung * ci_alb_lung +
      j_heart * ci_alb_heart +
      j_muscle * ci_alb_muscle +
      j_skin * ci_alb_skin +
      j_adipose * ci_alb_adipose +
      j_bone * ci_alb_bone +
      j_brain * ci_alb_brain +
      j_kidney * ci_alb_kidney +
      j_liver * ci_alb_liver +
      j_small_intestine * ci_alb_small_intestine +
      j_large_intestine * ci_alb_large_intestine +
      j_pancreas * ci_alb_pancreas +
      j_thymus * ci_alb_thymus +
      j_spleen * ci_alb_spleen +
      j_other * ci_alb_other) - lnlf * cln_alb

    # --- Renal filtrate holding pool, cleared to urine (Eq S11) ---
    d/dt(urine_alb) <- clren_alb * cv_alb_kidney - qurine * cbl_alb

    # ============================================================
    # endogenous IgG
    # ============================================================

    # --- Vascular plasma of each organ (Eq S3) ---
    d/dt(vp_lung_igg) <- (q_lung + j_lung) * cpl_igg - sumq_sys * cv_igg_lung -
      spino_igg * clup_lung * cv_igg_lung + 2 * clup_lung * fr * ceb_igg_lung - tp_igg_lung
    d/dt(vp_heart_igg) <- q_heart * cv_igg_lung - (q_heart - j_heart) * cv_igg_heart -
      spino_igg * clup_heart * cv_igg_heart + 2 * clup_heart * fr * ceb_igg_heart - tp_igg_heart
    d/dt(vp_muscle_igg) <- q_muscle * cv_igg_lung - (q_muscle - j_muscle) * cv_igg_muscle -
      spino_igg * clup_muscle * cv_igg_muscle + 2 * clup_muscle * fr * ceb_igg_muscle - tp_igg_muscle
    d/dt(vp_skin_igg) <- q_skin * cv_igg_lung - (q_skin - j_skin) * cv_igg_skin -
      spino_igg * clup_skin * cv_igg_skin + 2 * clup_skin * fr * ceb_igg_skin - tp_igg_skin
    d/dt(vp_adipose_igg) <- q_adipose * cv_igg_lung - (q_adipose - j_adipose) * cv_igg_adipose -
      spino_igg * clup_adipose * cv_igg_adipose + 2 * clup_adipose * fr * ceb_igg_adipose - tp_igg_adipose
    d/dt(vp_bone_igg) <- q_bone * cv_igg_lung - (q_bone - j_bone) * cv_igg_bone -
      spino_igg * clup_bone * cv_igg_bone + 2 * clup_bone * fr * ceb_igg_bone - tp_igg_bone
    d/dt(vp_brain_igg) <- q_brain * cv_igg_lung - (q_brain - j_brain) * cv_igg_brain -
      spino_igg * clup_brain * cv_igg_brain + 2 * clup_brain * fr * ceb_igg_brain - tp_igg_brain
    d/dt(vp_kidney_igg) <- q_kidney * cv_igg_lung - (q_kidney - j_kidney) * cv_igg_kidney - clren_igg * cv_igg_kidney -
      spino_igg * clup_kidney * cv_igg_kidney + 2 * clup_kidney * fr * ceb_igg_kidney - tp_igg_kidney
    d/dt(vp_liver_igg) <- q_liver * cv_igg_lung + (q_spleen - j_spleen) * cv_igg_spleen + (q_pancreas - j_pancreas) * cv_igg_pancreas + (q_small_intestine - j_small_intestine) * cv_igg_small_intestine + (q_large_intestine - j_large_intestine) * cv_igg_large_intestine - qout_liver * cv_igg_liver -
      spino_igg * clup_liver * cv_igg_liver + 2 * clup_liver * fr * ceb_igg_liver - tp_igg_liver
    d/dt(vp_small_intestine_igg) <- q_small_intestine * cv_igg_lung - (q_small_intestine - j_small_intestine) * cv_igg_small_intestine -
      spino_igg * clup_small_intestine * cv_igg_small_intestine + 2 * clup_small_intestine * fr * ceb_igg_small_intestine - tp_igg_small_intestine
    d/dt(vp_large_intestine_igg) <- q_large_intestine * cv_igg_lung - (q_large_intestine - j_large_intestine) * cv_igg_large_intestine -
      spino_igg * clup_large_intestine * cv_igg_large_intestine + 2 * clup_large_intestine * fr * ceb_igg_large_intestine - tp_igg_large_intestine
    d/dt(vp_pancreas_igg) <- q_pancreas * cv_igg_lung - (q_pancreas - j_pancreas) * cv_igg_pancreas -
      spino_igg * clup_pancreas * cv_igg_pancreas + 2 * clup_pancreas * fr * ceb_igg_pancreas - tp_igg_pancreas
    d/dt(vp_thymus_igg) <- q_thymus * cv_igg_lung - (q_thymus - j_thymus) * cv_igg_thymus -
      spino_igg * clup_thymus * cv_igg_thymus + 2 * clup_thymus * fr * ceb_igg_thymus - tp_igg_thymus
    d/dt(vp_spleen_igg) <- q_spleen * cv_igg_lung - (q_spleen - j_spleen) * cv_igg_spleen -
      spino_igg * clup_spleen * cv_igg_spleen + 2 * clup_spleen * fr * ceb_igg_spleen - tp_igg_spleen
    d/dt(vp_other_igg) <- q_other * cv_igg_lung - (q_other - j_other) * cv_igg_other -
      spino_igg * clup_other * cv_igg_other + 2 * clup_other * fr * ceb_igg_other - tp_igg_other

    # --- Interstitial space of each organ (Eq S4) ---
    d/dt(is_lung_igg) <- tp_igg_lung - (1 - sigis) * j_lung * ci_igg_lung -
      spino_igg * clup_lung * ci_igg_lung + 2 * clup_lung * (1 - fr) * ceb_igg_lung
    d/dt(is_heart_igg) <- tp_igg_heart - (1 - sigis) * j_heart * ci_igg_heart -
      spino_igg * clup_heart * ci_igg_heart + 2 * clup_heart * (1 - fr) * ceb_igg_heart
    d/dt(is_muscle_igg) <- tp_igg_muscle - (1 - sigis) * j_muscle * ci_igg_muscle -
      spino_igg * clup_muscle * ci_igg_muscle + 2 * clup_muscle * (1 - fr) * ceb_igg_muscle
    d/dt(is_skin_igg) <- tp_igg_skin - (1 - sigis) * j_skin * ci_igg_skin -
      spino_igg * clup_skin * ci_igg_skin + 2 * clup_skin * (1 - fr) * ceb_igg_skin
    d/dt(is_adipose_igg) <- tp_igg_adipose - (1 - sigis) * j_adipose * ci_igg_adipose -
      spino_igg * clup_adipose * ci_igg_adipose + 2 * clup_adipose * (1 - fr) * ceb_igg_adipose
    d/dt(is_bone_igg) <- tp_igg_bone - (1 - sigis) * j_bone * ci_igg_bone -
      spino_igg * clup_bone * ci_igg_bone + 2 * clup_bone * (1 - fr) * ceb_igg_bone
    d/dt(is_brain_igg) <- tp_igg_brain - (1 - sigis) * j_brain * ci_igg_brain -
      spino_igg * clup_brain * ci_igg_brain + 2 * clup_brain * (1 - fr) * ceb_igg_brain
    d/dt(is_kidney_igg) <- tp_igg_kidney - (1 - sigis) * j_kidney * ci_igg_kidney -
      spino_igg * clup_kidney * ci_igg_kidney + 2 * clup_kidney * (1 - fr) * ceb_igg_kidney
    d/dt(is_liver_igg) <- tp_igg_liver - (1 - sigis) * j_liver * ci_igg_liver -
      spino_igg * clup_liver * ci_igg_liver + 2 * clup_liver * (1 - fr) * ceb_igg_liver
    d/dt(is_small_intestine_igg) <- tp_igg_small_intestine - (1 - sigis) * j_small_intestine * ci_igg_small_intestine -
      spino_igg * clup_small_intestine * ci_igg_small_intestine + 2 * clup_small_intestine * (1 - fr) * ceb_igg_small_intestine
    d/dt(is_large_intestine_igg) <- tp_igg_large_intestine - (1 - sigis) * j_large_intestine * ci_igg_large_intestine -
      spino_igg * clup_large_intestine * ci_igg_large_intestine + 2 * clup_large_intestine * (1 - fr) * ceb_igg_large_intestine
    d/dt(is_pancreas_igg) <- tp_igg_pancreas - (1 - sigis) * j_pancreas * ci_igg_pancreas -
      spino_igg * clup_pancreas * ci_igg_pancreas + 2 * clup_pancreas * (1 - fr) * ceb_igg_pancreas
    d/dt(is_thymus_igg) <- tp_igg_thymus - (1 - sigis) * j_thymus * ci_igg_thymus -
      spino_igg * clup_thymus * ci_igg_thymus + 2 * clup_thymus * (1 - fr) * ceb_igg_thymus
    d/dt(is_spleen_igg) <- tp_igg_spleen - (1 - sigis) * j_spleen * ci_igg_spleen -
      spino_igg * clup_spleen * ci_igg_spleen + 2 * clup_spleen * (1 - fr) * ceb_igg_spleen
    d/dt(is_other_igg) <- tp_igg_other - (1 - sigis) * j_other * ci_igg_other -
      spino_igg * clup_other * ci_igg_other + 2 * clup_other * (1 - fr) * ceb_igg_other

    # --- Endosomal unbound / FcRn-bound (Eq S14, S15, S17, S18) ---
    bind_igg_lung <- kon_igg * ceu_igg_lung * cfr_igg_lung * ve_lung
    d/dt(eu_lung_igg) <- spino_igg * clup_lung * (cv_igg_lung + ci_igg_lung) -
      kdeg * eu_lung_igg - bind_igg_lung + koff_igg * eb_lung_igg
    d/dt(eb_lung_igg) <- bind_igg_lung - koff_igg * eb_lung_igg - 2 * clup_lung * ceb_igg_lung
    bind_igg_heart <- kon_igg * ceu_igg_heart * cfr_igg_heart * ve_heart
    d/dt(eu_heart_igg) <- spino_igg * clup_heart * (cv_igg_heart + ci_igg_heart) -
      kdeg * eu_heart_igg - bind_igg_heart + koff_igg * eb_heart_igg
    d/dt(eb_heart_igg) <- bind_igg_heart - koff_igg * eb_heart_igg - 2 * clup_heart * ceb_igg_heart
    bind_igg_muscle <- kon_igg * ceu_igg_muscle * cfr_igg_muscle * ve_muscle
    d/dt(eu_muscle_igg) <- spino_igg * clup_muscle * (cv_igg_muscle + ci_igg_muscle) -
      kdeg * eu_muscle_igg - bind_igg_muscle + koff_igg * eb_muscle_igg
    d/dt(eb_muscle_igg) <- bind_igg_muscle - koff_igg * eb_muscle_igg - 2 * clup_muscle * ceb_igg_muscle
    bind_igg_skin <- kon_igg * ceu_igg_skin * cfr_igg_skin * ve_skin
    d/dt(eu_skin_igg) <- spino_igg * clup_skin * (cv_igg_skin + ci_igg_skin) -
      kdeg * eu_skin_igg - bind_igg_skin + koff_igg * eb_skin_igg
    d/dt(eb_skin_igg) <- bind_igg_skin - koff_igg * eb_skin_igg - 2 * clup_skin * ceb_igg_skin
    bind_igg_adipose <- kon_igg * ceu_igg_adipose * cfr_igg_adipose * ve_adipose
    d/dt(eu_adipose_igg) <- spino_igg * clup_adipose * (cv_igg_adipose + ci_igg_adipose) -
      kdeg * eu_adipose_igg - bind_igg_adipose + koff_igg * eb_adipose_igg
    d/dt(eb_adipose_igg) <- bind_igg_adipose - koff_igg * eb_adipose_igg - 2 * clup_adipose * ceb_igg_adipose
    bind_igg_bone <- kon_igg * ceu_igg_bone * cfr_igg_bone * ve_bone
    d/dt(eu_bone_igg) <- spino_igg * clup_bone * (cv_igg_bone + ci_igg_bone) -
      kdeg * eu_bone_igg - bind_igg_bone + koff_igg * eb_bone_igg
    d/dt(eb_bone_igg) <- bind_igg_bone - koff_igg * eb_bone_igg - 2 * clup_bone * ceb_igg_bone
    bind_igg_brain <- kon_igg * ceu_igg_brain * cfr_igg_brain * ve_brain
    d/dt(eu_brain_igg) <- spino_igg * clup_brain * (cv_igg_brain + ci_igg_brain) -
      kdeg * eu_brain_igg - bind_igg_brain + koff_igg * eb_brain_igg
    d/dt(eb_brain_igg) <- bind_igg_brain - koff_igg * eb_brain_igg - 2 * clup_brain * ceb_igg_brain
    bind_igg_kidney <- kon_igg * ceu_igg_kidney * cfr_igg_kidney * ve_kidney
    d/dt(eu_kidney_igg) <- spino_igg * clup_kidney * (cv_igg_kidney + ci_igg_kidney) -
      kdeg * eu_kidney_igg - bind_igg_kidney + koff_igg * eb_kidney_igg
    d/dt(eb_kidney_igg) <- bind_igg_kidney - koff_igg * eb_kidney_igg - 2 * clup_kidney * ceb_igg_kidney
    bind_igg_liver <- kon_igg * ceu_igg_liver * cfr_igg_liver * ve_liver
    d/dt(eu_liver_igg) <- spino_igg * clup_liver * (cv_igg_liver + ci_igg_liver) -
      kdeg * eu_liver_igg - bind_igg_liver + koff_igg * eb_liver_igg
    d/dt(eb_liver_igg) <- bind_igg_liver - koff_igg * eb_liver_igg - 2 * clup_liver * ceb_igg_liver
    bind_igg_small_intestine <- kon_igg * ceu_igg_small_intestine * cfr_igg_small_intestine * ve_small_intestine
    d/dt(eu_small_intestine_igg) <- spino_igg * clup_small_intestine * (cv_igg_small_intestine + ci_igg_small_intestine) -
      kdeg * eu_small_intestine_igg - bind_igg_small_intestine + koff_igg * eb_small_intestine_igg
    d/dt(eb_small_intestine_igg) <- bind_igg_small_intestine - koff_igg * eb_small_intestine_igg - 2 * clup_small_intestine * ceb_igg_small_intestine
    bind_igg_large_intestine <- kon_igg * ceu_igg_large_intestine * cfr_igg_large_intestine * ve_large_intestine
    d/dt(eu_large_intestine_igg) <- spino_igg * clup_large_intestine * (cv_igg_large_intestine + ci_igg_large_intestine) -
      kdeg * eu_large_intestine_igg - bind_igg_large_intestine + koff_igg * eb_large_intestine_igg
    d/dt(eb_large_intestine_igg) <- bind_igg_large_intestine - koff_igg * eb_large_intestine_igg - 2 * clup_large_intestine * ceb_igg_large_intestine
    bind_igg_pancreas <- kon_igg * ceu_igg_pancreas * cfr_igg_pancreas * ve_pancreas
    d/dt(eu_pancreas_igg) <- spino_igg * clup_pancreas * (cv_igg_pancreas + ci_igg_pancreas) -
      kdeg * eu_pancreas_igg - bind_igg_pancreas + koff_igg * eb_pancreas_igg
    d/dt(eb_pancreas_igg) <- bind_igg_pancreas - koff_igg * eb_pancreas_igg - 2 * clup_pancreas * ceb_igg_pancreas
    bind_igg_thymus <- kon_igg * ceu_igg_thymus * cfr_igg_thymus * ve_thymus
    d/dt(eu_thymus_igg) <- spino_igg * clup_thymus * (cv_igg_thymus + ci_igg_thymus) -
      kdeg * eu_thymus_igg - bind_igg_thymus + koff_igg * eb_thymus_igg
    d/dt(eb_thymus_igg) <- bind_igg_thymus - koff_igg * eb_thymus_igg - 2 * clup_thymus * ceb_igg_thymus
    bind_igg_spleen <- kon_igg * ceu_igg_spleen * cfr_igg_spleen * ve_spleen
    d/dt(eu_spleen_igg) <- spino_igg * clup_spleen * (cv_igg_spleen + ci_igg_spleen) -
      kdeg * eu_spleen_igg - bind_igg_spleen + koff_igg * eb_spleen_igg
    d/dt(eb_spleen_igg) <- bind_igg_spleen - koff_igg * eb_spleen_igg - 2 * clup_spleen * ceb_igg_spleen
    bind_igg_other <- kon_igg * ceu_igg_other * cfr_igg_other * ve_other
    d/dt(eu_other_igg) <- spino_igg * clup_other * (cv_igg_other + ci_igg_other) -
      kdeg * eu_other_igg - bind_igg_other + koff_igg * eb_other_igg
    d/dt(eb_other_igg) <- bind_igg_other - koff_igg * eb_other_igg - 2 * clup_other * ceb_igg_other

    # --- Central plasma (Eq S1) ---
    d/dt(plasma_igg) <- (q_heart - j_heart) * cv_igg_heart +
      (q_muscle - j_muscle) * cv_igg_muscle +
      (q_skin - j_skin) * cv_igg_skin +
      (q_adipose - j_adipose) * cv_igg_adipose +
      (q_bone - j_bone) * cv_igg_bone +
      (q_brain - j_brain) * cv_igg_brain +
      (q_kidney - j_kidney) * cv_igg_kidney +
      (q_thymus - j_thymus) * cv_igg_thymus +
      (q_other - j_other) * cv_igg_other +
      qout_liver * cv_igg_liver +
      lnlf * cln_igg - (q_lung + j_lung) * cpl_igg +
      0.47 * vpl   # zero-order synthesis ksyn * V_plasma

    # --- Lymph node (Eq S2) ---
    d/dt(lnode_igg) <- (1 - sigis) * (
      j_lung * ci_igg_lung +
      j_heart * ci_igg_heart +
      j_muscle * ci_igg_muscle +
      j_skin * ci_igg_skin +
      j_adipose * ci_igg_adipose +
      j_bone * ci_igg_bone +
      j_brain * ci_igg_brain +
      j_kidney * ci_igg_kidney +
      j_liver * ci_igg_liver +
      j_small_intestine * ci_igg_small_intestine +
      j_large_intestine * ci_igg_large_intestine +
      j_pancreas * ci_igg_pancreas +
      j_thymus * ci_igg_thymus +
      j_spleen * ci_igg_spleen +
      j_other * ci_igg_other) - lnlf * cln_igg

    # --- Renal filtrate holding pool, cleared to urine (Eq S11) ---
    d/dt(urine_igg) <- clren_igg * cv_igg_kidney - qurine * cbl_igg

    # ================= Free FcRn binding sites (Eq S19, S20) =================
    # Endogenous IgG and durvalumab compete for the same FcRn site; albumin
    # binds a separate, non-competing site.
    d/dt(fr_lung_igg) <- 2 * clup_lung * (ceb_mab_lung + ceb_igg_lung) -
      bind_mab_lung - bind_igg_lung + koff_mab * eb_lung + koff_igg * eb_lung_igg
    d/dt(fr_heart_igg) <- 2 * clup_heart * (ceb_mab_heart + ceb_igg_heart) -
      bind_mab_heart - bind_igg_heart + koff_mab * eb_heart + koff_igg * eb_heart_igg
    d/dt(fr_muscle_igg) <- 2 * clup_muscle * (ceb_mab_muscle + ceb_igg_muscle) -
      bind_mab_muscle - bind_igg_muscle + koff_mab * eb_muscle + koff_igg * eb_muscle_igg
    d/dt(fr_skin_igg) <- 2 * clup_skin * (ceb_mab_skin + ceb_igg_skin) -
      bind_mab_skin - bind_igg_skin + koff_mab * eb_skin + koff_igg * eb_skin_igg
    d/dt(fr_adipose_igg) <- 2 * clup_adipose * (ceb_mab_adipose + ceb_igg_adipose) -
      bind_mab_adipose - bind_igg_adipose + koff_mab * eb_adipose + koff_igg * eb_adipose_igg
    d/dt(fr_bone_igg) <- 2 * clup_bone * (ceb_mab_bone + ceb_igg_bone) -
      bind_mab_bone - bind_igg_bone + koff_mab * eb_bone + koff_igg * eb_bone_igg
    d/dt(fr_brain_igg) <- 2 * clup_brain * (ceb_mab_brain + ceb_igg_brain) -
      bind_mab_brain - bind_igg_brain + koff_mab * eb_brain + koff_igg * eb_brain_igg
    d/dt(fr_kidney_igg) <- 2 * clup_kidney * (ceb_mab_kidney + ceb_igg_kidney) -
      bind_mab_kidney - bind_igg_kidney + koff_mab * eb_kidney + koff_igg * eb_kidney_igg
    d/dt(fr_liver_igg) <- 2 * clup_liver * (ceb_mab_liver + ceb_igg_liver) -
      bind_mab_liver - bind_igg_liver + koff_mab * eb_liver + koff_igg * eb_liver_igg
    d/dt(fr_small_intestine_igg) <- 2 * clup_small_intestine * (ceb_mab_small_intestine + ceb_igg_small_intestine) -
      bind_mab_small_intestine - bind_igg_small_intestine + koff_mab * eb_small_intestine + koff_igg * eb_small_intestine_igg
    d/dt(fr_large_intestine_igg) <- 2 * clup_large_intestine * (ceb_mab_large_intestine + ceb_igg_large_intestine) -
      bind_mab_large_intestine - bind_igg_large_intestine + koff_mab * eb_large_intestine + koff_igg * eb_large_intestine_igg
    d/dt(fr_pancreas_igg) <- 2 * clup_pancreas * (ceb_mab_pancreas + ceb_igg_pancreas) -
      bind_mab_pancreas - bind_igg_pancreas + koff_mab * eb_pancreas + koff_igg * eb_pancreas_igg
    d/dt(fr_thymus_igg) <- 2 * clup_thymus * (ceb_mab_thymus + ceb_igg_thymus) -
      bind_mab_thymus - bind_igg_thymus + koff_mab * eb_thymus + koff_igg * eb_thymus_igg
    d/dt(fr_spleen_igg) <- 2 * clup_spleen * (ceb_mab_spleen + ceb_igg_spleen) -
      bind_mab_spleen - bind_igg_spleen + koff_mab * eb_spleen + koff_igg * eb_spleen_igg
    d/dt(fr_other_igg) <- 2 * clup_other * (ceb_mab_other + ceb_igg_other) -
      bind_mab_other - bind_igg_other + koff_mab * eb_other + koff_igg * eb_other_igg

    d/dt(fr_lung_alb) <- 2 * clup_lung * ceb_alb_lung - bind_alb_lung +
      koff_alb * eb_lung_alb
    d/dt(fr_heart_alb) <- 2 * clup_heart * ceb_alb_heart - bind_alb_heart +
      koff_alb * eb_heart_alb
    d/dt(fr_muscle_alb) <- 2 * clup_muscle * ceb_alb_muscle - bind_alb_muscle +
      koff_alb * eb_muscle_alb
    d/dt(fr_skin_alb) <- 2 * clup_skin * ceb_alb_skin - bind_alb_skin +
      koff_alb * eb_skin_alb
    d/dt(fr_adipose_alb) <- 2 * clup_adipose * ceb_alb_adipose - bind_alb_adipose +
      koff_alb * eb_adipose_alb
    d/dt(fr_bone_alb) <- 2 * clup_bone * ceb_alb_bone - bind_alb_bone +
      koff_alb * eb_bone_alb
    d/dt(fr_brain_alb) <- 2 * clup_brain * ceb_alb_brain - bind_alb_brain +
      koff_alb * eb_brain_alb
    d/dt(fr_kidney_alb) <- 2 * clup_kidney * ceb_alb_kidney - bind_alb_kidney +
      koff_alb * eb_kidney_alb
    d/dt(fr_liver_alb) <- 2 * clup_liver * ceb_alb_liver - bind_alb_liver +
      koff_alb * eb_liver_alb
    d/dt(fr_small_intestine_alb) <- 2 * clup_small_intestine * ceb_alb_small_intestine - bind_alb_small_intestine +
      koff_alb * eb_small_intestine_alb
    d/dt(fr_large_intestine_alb) <- 2 * clup_large_intestine * ceb_alb_large_intestine - bind_alb_large_intestine +
      koff_alb * eb_large_intestine_alb
    d/dt(fr_pancreas_alb) <- 2 * clup_pancreas * ceb_alb_pancreas - bind_alb_pancreas +
      koff_alb * eb_pancreas_alb
    d/dt(fr_thymus_alb) <- 2 * clup_thymus * ceb_alb_thymus - bind_alb_thymus +
      koff_alb * eb_thymus_alb
    d/dt(fr_spleen_alb) <- 2 * clup_spleen * ceb_alb_spleen - bind_alb_spleen +
      koff_alb * eb_spleen_alb
    d/dt(fr_other_alb) <- 2 * clup_other * ceb_alb_other - bind_alb_other +
      koff_alb * eb_other_alb

    # ================= Initial conditions =================
    # Appendix S3 initialAssignment rules: endogenous albumin and IgG start at
    # their steady-state plasma concentration in the central plasma AND in every
    # organ vascular space; interstitial, endosomal, lymph and filtrate pools
    # start empty. Durvalumab starts at zero everywhere (dose supplied by data).
    cssalb <- exp(lcss_alb)
    cssigg <- exp(lcss_igg)
    plasma_alb(0) <- cssalb * vpl
    plasma_igg(0) <- cssigg * vpl
    vp_lung_alb(0) <- cssalb * vv_lung; vp_lung_igg(0) <- cssigg * vv_lung
    vp_heart_alb(0) <- cssalb * vv_heart; vp_heart_igg(0) <- cssigg * vv_heart
    vp_muscle_alb(0) <- cssalb * vv_muscle; vp_muscle_igg(0) <- cssigg * vv_muscle
    vp_skin_alb(0) <- cssalb * vv_skin; vp_skin_igg(0) <- cssigg * vv_skin
    vp_adipose_alb(0) <- cssalb * vv_adipose; vp_adipose_igg(0) <- cssigg * vv_adipose
    vp_bone_alb(0) <- cssalb * vv_bone; vp_bone_igg(0) <- cssigg * vv_bone
    vp_brain_alb(0) <- cssalb * vv_brain; vp_brain_igg(0) <- cssigg * vv_brain
    vp_kidney_alb(0) <- cssalb * vv_kidney; vp_kidney_igg(0) <- cssigg * vv_kidney
    vp_liver_alb(0) <- cssalb * vv_liver; vp_liver_igg(0) <- cssigg * vv_liver
    vp_small_intestine_alb(0) <- cssalb * vv_small_intestine; vp_small_intestine_igg(0) <- cssigg * vv_small_intestine
    vp_large_intestine_alb(0) <- cssalb * vv_large_intestine; vp_large_intestine_igg(0) <- cssigg * vv_large_intestine
    vp_pancreas_alb(0) <- cssalb * vv_pancreas; vp_pancreas_igg(0) <- cssigg * vv_pancreas
    vp_thymus_alb(0) <- cssalb * vv_thymus; vp_thymus_igg(0) <- cssigg * vv_thymus
    vp_spleen_alb(0) <- cssalb * vv_spleen; vp_spleen_igg(0) <- cssigg * vv_spleen
    vp_other_alb(0) <- cssalb * vv_other; vp_other_igg(0) <- cssigg * vv_other

    # Free FcRn pools start fully unoccupied at their total site concentration.
    fcrnigg <- exp(lfcrn_igg)
    fcrnalb <- exp(lfcrn_alb)
    fr_lung_igg(0) <- fcrnigg * ve_lung; fr_lung_alb(0) <- fcrnalb * ve_lung
    fr_heart_igg(0) <- fcrnigg * ve_heart; fr_heart_alb(0) <- fcrnalb * ve_heart
    fr_muscle_igg(0) <- fcrnigg * ve_muscle; fr_muscle_alb(0) <- fcrnalb * ve_muscle
    fr_skin_igg(0) <- fcrnigg * ve_skin; fr_skin_alb(0) <- fcrnalb * ve_skin
    fr_adipose_igg(0) <- fcrnigg * ve_adipose; fr_adipose_alb(0) <- fcrnalb * ve_adipose
    fr_bone_igg(0) <- fcrnigg * ve_bone; fr_bone_alb(0) <- fcrnalb * ve_bone
    fr_brain_igg(0) <- fcrnigg * ve_brain; fr_brain_alb(0) <- fcrnalb * ve_brain
    fr_kidney_igg(0) <- fcrnigg * ve_kidney; fr_kidney_alb(0) <- fcrnalb * ve_kidney
    fr_liver_igg(0) <- fcrnigg * ve_liver; fr_liver_alb(0) <- fcrnalb * ve_liver
    fr_small_intestine_igg(0) <- fcrnigg * ve_small_intestine; fr_small_intestine_alb(0) <- fcrnalb * ve_small_intestine
    fr_large_intestine_igg(0) <- fcrnigg * ve_large_intestine; fr_large_intestine_alb(0) <- fcrnalb * ve_large_intestine
    fr_pancreas_igg(0) <- fcrnigg * ve_pancreas; fr_pancreas_alb(0) <- fcrnalb * ve_pancreas
    fr_thymus_igg(0) <- fcrnigg * ve_thymus; fr_thymus_alb(0) <- fcrnalb * ve_thymus
    fr_spleen_igg(0) <- fcrnigg * ve_spleen; fr_spleen_alb(0) <- fcrnalb * ve_spleen
    fr_other_igg(0) <- fcrnigg * ve_other; fr_other_alb(0) <- fcrnalb * ve_other

    # ================= Observations =================
    # Durvalumab plasma concentration. States are nmol and volumes mL, so the
    # native concentration unit is uM; mwMab converts to mg/L for comparison
    # with the published popPK profile.
    # Molecular weights. ALB_MW = 67000 g/mol is given directly in the
    # Appendix S3 Parameters sheet. The IgG value is not printed anywhere in
    # the paper; 150000 g/mol is recovered EXACTLY from two independent
    # Appendix S3 quantities (a back-solve, not an assumed class value):
    #   ENIGG_theta_S = 1 - 0.8489*exp(-4e-5*MW) = 0.99789578727722794 -> MW = 150000
    #   RepeatDoses "1 mg/kg" = 473.3333 nmol = 71 kg * 1 mg/kg / 150000 g/mol
    mwMab <- 150000; mwAlb <- 67000; mwIgg <- 150000
    # Cc is reported in nmol/L to stay dimensionally consistent with the nmol
    # dosing unit. Multiply by mwMab/1e6 for mg/L (1 nmol/L = 1.5e-4 mg/L).
    Cc      <- cpl_mab * 1000
    albumin <- cpl_alb * mwAlb / 1e6   # g/L (Appendix S3 ALB_gL_assignment)
    igg     <- cpl_igg * mwIgg / 1e6   # g/L
    kdegOut <- kdeg

    Cc      ~ prop(propSd)
  })
}
