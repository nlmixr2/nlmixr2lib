Jones_2019_mab_human_pbpk <- function() {
  description <- "PBPK / QSP (whole-body, 15 organs, 620 ODE states). Platform model predicting human linear pharmacokinetics of an IgG1 monoclonal antibody a priori from two in vitro inputs: an AC-SINS polyspecificity score and the FcRn binding affinity at pH 6.0. Each organ carries vascular, vascular-side membrane, three pH-resolved endosomal transit compartments (early pH 7.4, sorting pH 6.0, recycling pH 7.4), interstitial-side membrane and interstitial space. Dosed antibody and endogenous IgG compete for a shared FcRn pool with explicit 1:1 and 2:1 stoichiometry, and non-specific charge-mediated binding to cell-membrane sites provides an FcRn-independent clearance route driven by the AC-SINS score. Extends the Shah and Betts 2012 platform topology."
  reference <- "Jones HM, Zhang Z, Jasper P, Luo H, Avery LB, King LE, Neubert H, Barton HA, Betts AM, Webster R. A Physiologically-Based Pharmacokinetic Model for the Prediction of Monoclonal Antibody Pharmacokinetics From In Vitro Data. CPT Pharmacometrics Syst Pharmacol. 2019;8(10):738-747. doi:10.1002/psp4.12461"
  vignette <- "Jones_2019_mab_pbpk"
  units <- list(time = "h", dosing = "umol", concentration = "ug/mL")

  # Structural equations: Supplementary "PBPK Model Equations" (s001.pdf).
  # Every numeric constant below is transcribed from the deposited Berkeley
  # Madonna listing (Supplementary Model Code, s002.txt), which is the
  # executable form of those equations and carries the values at full
  # precision where main-text Table 1 / Table 2 round them. Physiological
  # volumes and flows are Table S2 (human, 71 kg male), reproduced verbatim
  # in the Madonna listing. See the vignette source-trace table.
  #
  # Madonna states are CONCENTRATIONS (uM) and each equation divides its flux
  # sum by the compartment volume. The states below are AMOUNTS (umol), so
  # each d/dt() here is the Madonna numerator verbatim and every
  # concentration is recovered as amount / volume. This is an exact
  # re-arrangement, not a re-parameterisation.

  compartmentData <- list(
    vp_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    vp_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    vp_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    vp_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    vp_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    vp_adipose = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    vp_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    vp_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    vp_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    vp_small_intestine = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    vp_large_intestine = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    vp_pancreas = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    vp_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    vp_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    vp_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_adipose = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memvas_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memvas_pancreas = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    endoearly_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_adipose = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearly_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearly_pancreas = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_adipose = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endosort_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endosort_pancreas = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_adipose = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecyc_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecyc_pancreas = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    memint_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memint_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memint_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memint_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memint_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memint_adipose = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memint_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memint_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memint_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memint_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memint_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memint_pancreas = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memint_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memint_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memint_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    is_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    is_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    is_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    is_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    is_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    is_adipose = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    is_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    is_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    is_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    is_small_intestine = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    is_large_intestine = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    is_pancreas = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    is_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    is_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    is_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_adipose = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memvasfr1_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memvasfr1_pancreas = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_adipose = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memvasfr2_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memvasfr2_pancreas = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    endoearlyfr1_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_adipose = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearlyfr1_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearlyfr1_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearlyfr1_pancreas = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearlyfr1_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_adipose = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearlyfr2_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearlyfr2_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearlyfr2_pancreas = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearlyfr2_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_adipose = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endosortfr1_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endosortfr1_pancreas = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endosortfr1_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_adipose = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endosortfr2_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endosortfr2_pancreas = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endosortfr2_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_adipose = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecycfr1_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecycfr1_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecycfr1_pancreas = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecycfr1_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_adipose = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecycfr2_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecycfr2_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecycfr2_pancreas = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecycfr2_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "endosome", verified = TRUE),
    memintfr1_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_adipose = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memintfr1_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memintfr1_pancreas = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_adipose = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memintfr2_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memintfr2_pancreas = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasns_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasns_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasns_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasns_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasns_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasns_adipose = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasns_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasns_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasns_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasns_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memvasns_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memvasns_pancreas = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasns_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasns_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memvasns_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintns_lung = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintns_liver = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintns_heart = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintns_muscle = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintns_skin = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintns_adipose = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintns_bone = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintns_brain = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintns_kidney = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintns_small_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memintns_large_intestine = list(
      analyte = "monoclonal antibody",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memintns_pancreas = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintns_thymus = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintns_spleen = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    memintns_other = list(analyte = "monoclonal antibody", units = "umol", specimen = "tissue", verified = TRUE),
    vp_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    vp_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    vp_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    vp_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    vp_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    vp_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    vp_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    vp_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    vp_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    vp_small_intestine_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    vp_large_intestine_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    vp_pancreas_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    vp_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    vp_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    vp_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_small_intestine_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_large_intestine_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_pancreas_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvas_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    endoearly_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_small_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearly_large_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearly_pancreas_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearly_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_small_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endosort_large_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endosort_pancreas_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosort_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_small_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecyc_large_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecyc_pancreas_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecyc_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    memint_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memint_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memint_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memint_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memint_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memint_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memint_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memint_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memint_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memint_small_intestine_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memint_large_intestine_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memint_pancreas_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memint_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memint_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memint_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    is_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    is_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    is_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    is_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    is_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    is_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    is_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    is_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    is_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    is_small_intestine_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    is_large_intestine_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    is_pancreas_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    is_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    is_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    is_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_small_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memvasfr1_large_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memvasfr1_pancreas_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr1_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_small_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memvasfr2_large_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memvasfr2_pancreas_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memvasfr2_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    endoearlyfr1_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_small_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearlyfr1_large_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearlyfr1_pancreas_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearlyfr1_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr1_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_small_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearlyfr2_large_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearlyfr2_pancreas_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endoearlyfr2_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endoearlyfr2_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_small_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endosortfr1_large_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endosortfr1_pancreas_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr1_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_small_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endosortfr2_large_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endosortfr2_pancreas_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endosortfr2_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_small_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecycfr1_large_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecycfr1_pancreas_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecycfr1_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr1_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_small_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecycfr2_large_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecycfr2_pancreas_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "endosome",
      verified = TRUE
    ),
    endorecycfr2_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    endorecycfr2_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "endosome", verified = TRUE),
    memintfr1_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_small_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memintfr1_large_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memintfr1_pancreas_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr1_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_lung_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_liver_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_heart_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_muscle_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_skin_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_adipose_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_bone_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_brain_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_kidney_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_small_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memintfr2_large_intestine_igg = list(
      analyte = "endogenous IgG",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    memintfr2_pancreas_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_thymus_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_spleen_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    memintfr2_other_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemvas_lung = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemvas_liver = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemvas_heart = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemvas_muscle = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemvas_skin = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemvas_adipose = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemvas_bone = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemvas_brain = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemvas_kidney = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemvas_small_intestine = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemvas_large_intestine = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemvas_pancreas = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemvas_thymus = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemvas_spleen = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemvas_other = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnendoearly_lung = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendoearly_liver = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendoearly_heart = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendoearly_muscle = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendoearly_skin = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendoearly_adipose = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendoearly_bone = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendoearly_brain = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendoearly_kidney = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendoearly_small_intestine = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendoearly_large_intestine = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendoearly_pancreas = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendoearly_thymus = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendoearly_spleen = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendoearly_other = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendosort_lung = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendosort_liver = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendosort_heart = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendosort_muscle = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendosort_skin = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendosort_adipose = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendosort_bone = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendosort_brain = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendosort_kidney = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendosort_small_intestine = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendosort_large_intestine = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendosort_pancreas = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendosort_thymus = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendosort_spleen = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendosort_other = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendorecyc_lung = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendorecyc_liver = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendorecyc_heart = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendorecyc_muscle = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendorecyc_skin = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendorecyc_adipose = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendorecyc_bone = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendorecyc_brain = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendorecyc_kidney = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendorecyc_small_intestine = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendorecyc_large_intestine = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendorecyc_pancreas = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendorecyc_thymus = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendorecyc_spleen = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnendorecyc_other = list(analyte = "free FcRn", units = "umol", specimen = "endosome", verified = TRUE),
    fcrnmemint_lung = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemint_liver = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemint_heart = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemint_muscle = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemint_skin = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemint_adipose = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemint_bone = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemint_brain = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemint_kidney = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemint_small_intestine = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemint_large_intestine = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemint_pancreas = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemint_thymus = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemint_spleen = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    fcrnmemint_other = list(analyte = "free FcRn", units = "umol", specimen = "tissue", verified = TRUE),
    central = list(analyte = "monoclonal antibody", units = "umol", specimen = "plasma", verified = TRUE),
    lnode = list(analyte = "monoclonal antibody", units = "umol", specimen = "lymph", verified = TRUE),
    lnode_igg = list(analyte = "endogenous IgG", units = "umol", specimen = "lymph", verified = TRUE),
    depot_sc = list(analyte = "monoclonal antibody", units = "umol", specimen = "administration site", verified = TRUE),
    auc_central = list(analyte = "monoclonal antibody", units = "umol", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = "12 monoclonal antibodies (training set 7: mAbs 1-3, 5-8; test set 5: mAbs 11-14, 23); healthy volunteers or patients",
    n_studies = "single-dose intravenous clinical PK studies, in-house (Pfizer) or literature",
    disease_state = "healthy volunteers or patients; linear (saturating-dose) PK only",
    dose_range = "saturating doses where PK was linear",
    notes = paste(
      "Whole-body catabolic capacity was calibrated jointly against Tg32 mouse and human data,",
      "including plasma profiles from FcRn-knockout mice and from humans with",
      "beta-2-microglobulin-deficiency (familial hypercatabolic hypoproteinaemia).",
      "The model is deterministic: the source reports no between-subject variability",
      "and no residual-error model."
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight. The deposited listing carries BW = 71 kg only to convert a mg/kg dose into an amount; no model parameter is scaled by body weight.",
      units = "kg",
      type = "continuous",
      notes = "Not referenced in model(): doses are supplied directly in umol, so body weight enters only at the event table."
    )
  )

  ini({
    # ---- Calibrated parameters (Table 1 "Fitted" and Table 2) ----------
    # Left un-fixed because the source ESTIMATED them; no uncertainty is
    # reported for any of them, so they carry point estimates only.
    lendoscale <- log(603.7)      ; label("Human endothelial-cell scale factor relative to the Tg32 mouse (-)")   # Madonna human_endo_scale_factor; 1.422e9 * 603.7 = 0.86E12 = Table 1 human N_endo
    lkdratio   <- log(220.11)     ; label("Ratio of mAb-FcRn Kd at pH 7.4 to Kd at pH 6.0 (-)")                   # Table 1 Kd 7.4/Kd 6.0 ratio = 220; Madonna KD7_WT/KD6_WT = 154077/700 = 220.11
    lkonratio  <- log(83.7)       ; label("Ratio of first to second mAb-FcRn association rate constant (-)")      # Table 1 k on ratio (FcRn binding #1/#2) = 83.7; gives k_on_2nd = 8.06E7/83.7 = 9.63E5 1/M/h as reported in Methods
    lpsa       <- log(1.8051)     ; label("AC-SINS scaling parameter a (-)")                                      # Table 2 human PS_a = 1.81; Madonna PS_a = 1.8051
    lpsb       <- log(0.2624)     ; label("AC-SINS scaling parameter b (-)")                                      # Table 2 human PS_b = -0.262 (entered here as the magnitude; the equation subtracts it); Madonna PS_b = 0.2624
    lcmem      <- log(18.5)       ; label("Cell-membrane site density available for non-specific mAb binding (uM)")  # Table 2 C mem = 18.5 uM; Methods: C mem in human was set equal to C mem in Tg32 mouse
    lkintps    <- log(0.0380)     ; label("Internalisation rate of membrane-bound mAb from non-specific binding (1/h)")  # Table 2 human k int_PS = 0.0380 1/h

    # ---- Fixed parameters (Table 1 "Fixed" block and Methods) ----------
    lkd6       <- fixed(log(700)) ; label("mAb-FcRn equilibrium dissociation constant at pH 6.0 (nM)")            # Methods: binding affinity of 700 nM at pH 6.0 from Biacore across mAbs with no non-specific interactions; Madonna KD6_WT
    lclup      <- fixed(log(150)) ; label("Pinocytotic uptake rate of endothelial cells (nL/h per 1E6 cells)")    # Table 1 Uptake rate (CL up) = 150
    lpinotime  <- fixed(log(0.18)); label("Endosomal transit time (h)")                                           # Methods: T endo assumed 10.8 minutes = 0.18 h (Table 1 rounds to 11 min)
    lfcrntot   <- fixed(log(1022)); label("Total body FcRn amount (nmol)")                                        # Table 1 FcRn tissue levels 1022 nmol (human); Madonna total_FcRn_in_nmole_par
    lnendomouse <- fixed(log(1.422e9)) ; label("Tg32 mouse total endothelial-cell number, the scaling anchor (cells)")  # Table 1 Tg32 N endo = 14.2E8; Madonna Total_Endothelial_Cell
    lkdegfcrnab <- fixed(log(log(2) / 11.1)) ; label("Degradation rate constant of FcRn-bound mAb at pH 7.4 (1/h)")   # Table 1 kdeg_FcRn_Ab = 0.062 at the 11.1 h natural FcRn turnover half-life; Madonna LOGN(2)/11.1 = 0.06245 1/h. Table 1 labels the units min-1, which is a units typo: 0.062 is the per-HOUR value - see vignette Errata
    lkon6      <- fixed(log(80.6)); label("mAb-FcRn association rate constant at pH 6.0 (1/uM/h)")                # Methods: k on_1st = 8.06E+7 1/M/h = 80.6 1/uM/h; Madonna k_on_6_EXG
    lkon7      <- fixed(log(3.22)); label("mAb-FcRn association rate constant at pH 7.4 (1/uM/h)")                # Madonna k_on_7_EXG = 1.61E7/5 = 3.22E6 1/M/h = 3.22 1/uM/h
    probdeg    <- fixed(0.95)     ; label("Probability of mAb degradation in the absence of FcRn binding (fraction)")  # Madonna Prob_deg = 0.95. Table 1 prints 98%; the deposited executable value is used - see vignette Errata
    fr         <- fixed(0.715)    ; label("Apical (vascular-side) fraction of pinocytosis and recycling (-)")     # Methods: apical recycling fraction set at 0.715 (Shah and Betts); basolateral = 1 - FR
    frecycle   <- fixed(0.99)     ; label("Fraction of free FcRn returned directly to the endosomal pool (-)")    # Methods: FcRn_recycle_fraction set to 0.99
    sigis      <- fixed(0.2)      ; label("Lymphatic (interstitial) reflection coefficient, all organs (-)")      # Madonna sigma_IS[1..N_Organs] = 0.2
    e6apct     <- fixed(0.33)     ; label("Sorting-endosome fraction of total endosomal volume (-)")              # Madonna E6a_Vol_Pct = 0.33; early and recycling endosomes split the remainder equally
    tauvm      <- fixed(0.0166666666666667) ; label("Residence time of the vascular-side membrane compartment (h)")   # Madonna tau_VM = 1/60 h
    tauism     <- fixed(0.0166666666666667) ; label("Residence time of the interstitial-side membrane compartment (h)")  # Madonna tau_ISM = 1/60 h
    cigg0      <- fixed(66.6666666666667)   ; label("Endogenous IgG plasma concentration, held constant (uM)")    # Madonna EDG_mg_ml = 10 mg/mL at MW 150000 -> 66.67 uM; d/dt(C_EDG_Plasma) = 0
    mwmab      <- fixed(150000)   ; label("Monoclonal antibody molecular weight (g/mol)")                         # Madonna MW_EDG = 150000
    lka        <- fixed(log(0.0108333333333333)) ; label("First-order subcutaneous absorption rate constant (1/h)")   # Madonna ka = 0.26/24 1/h
    fsc        <- fixed(0.60)     ; label("Subcutaneous bioavailability (fraction)")                              # Madonna F = 0.60
    psscore    <- fixed(0)        ; label("AC-SINS polyspecificity score of the simulated antibody (-)")          # Madonna PS_Score = 0 (default: an antibody with no non-specific interaction). Assay range 0-25

    # The source is a deterministic PBPK platform fitted to mean plasma
    # profiles. It reports neither between-subject variability nor a
    # residual-error model, so residual error is fixed at zero per the
    # standing policy for unreported RUV (documented in the vignette Errata).
    propSd     <- fixed(0)        ; label("Proportional residual error on plasma mAb concentration (fraction)")
  })

  model({
    # ---------------- derived from ini() ----------------
    kd6      <- exp(lkd6)
    kd7      <- kd6 * exp(lkdratio)
    kon6     <- exp(lkon6)
    kon7     <- exp(lkon7)
    konratio <- exp(lkonratio)
    kon6b    <- kon6 / konratio      # second (2:1) association, pH 6.0
    kon7b    <- kon7 / konratio      # second (2:1) association, pH 7.4
    # Kd of the second binding event is scaled by the same ratio, so
    # koff_2nd = (Kd * ratio / 1000) * (kon / ratio) = koff_1st. Methods:
    # "k off_2nd was assumed to be the same as k off_1st".
    koff6    <- (kd6 / 1000) * kon6
    koff7    <- (kd7 / 1000) * kon7
    koff6b   <- (kd6 * konratio / 1000) * kon6b
    koff7b   <- (kd7 * konratio / 1000) * kon7b
    # Endogenous IgG is wild-type IgG: same affinities and rate constants
    # as the exogenous wild-type antibody (Madonna KD_6_EDG / KD_7_EDG
    # equal KD6_WT / KD7_WT and k_on_*_EDG equal k_on_*_EXG).
    kdegab   <- exp(lkdegfcrnab)
    cmem     <- exp(lcmem)
    kintps   <- exp(lkintps)
    konps    <- kon6                 # Methods: k onPS assumed equal to the mAb-FcRn k on
    # AC-SINS score -> non-specific binding affinity (Supplementary
    # "Binding affinity from polyspecificity"): log10 Kd = exp(PS_a - PS_b * PS_Score), Kd in nM.
    pskd     <- 10^(exp(exp(lpsa) - exp(lpsb) * psscore))
    koffps   <- (pskd / 1000) * konps
    ka       <- exp(lka)

    # ---------------- physiological constants ----------------
    # Table S2 (human, 71 kg male); volumes in L, plasma flows in L/h.
    vv_lung <-            0.055; vis_lung <-              0.3; vorg_lung <-                1; plq_lung <-          181.913; sigv_lung <-             0.95; ecf_lung <-           0.0834
    vv_liver <-            0.183; vis_liver <-            0.429; vorg_liver <-            2.143; plq_liver <-            13.21; sigv_liver <-             0.85; ecf_liver <-           0.1877
    vv_heart <-           0.0131; vis_heart <-           0.0488; vorg_heart <-            0.341; plq_heart <-            7.752; sigv_heart <-             0.95; ecf_heart <-           0.0011
    vv_muscle <-            0.662; vis_muscle <-             3.91; vorg_muscle <-           30.078; plq_muscle <-           33.469; sigv_muscle <-             0.95; ecf_muscle <-           0.1928
    vv_skin <-            0.127; vis_skin <-            1.125; vorg_skin <-            3.408; plq_skin <-           11.626; sigv_skin <-             0.95; ecf_skin <-           0.0819
    vv_adipose <-            0.148; vis_adipose <-            2.289; vorg_adipose <-           13.465; plq_adipose <-           11.233; sigv_adipose <-             0.95; ecf_adipose <-           0.0999
    vv_bone <-            0.224; vis_bone <-            1.891; vorg_bone <-           10.165; plq_bone <-            2.591; sigv_bone <-             0.85; ecf_bone <-           0.1478
    vv_brain <-           0.0319; vis_brain <-            0.261; vorg_brain <-             1.45; plq_brain <-           21.453; sigv_brain <-             0.99; ecf_brain <-           0.0115
    vv_kidney <-           0.0182; vis_kidney <-           0.0498; vorg_kidney <-            0.332; plq_kidney <-           36.402; sigv_kidney <-              0.9; ecf_kidney <-           0.0157
    vv_small_intestine <-          0.00615; vis_small_intestine <-           0.0671; vorg_small_intestine <-            0.385; plq_small_intestine <-           12.368; sigv_small_intestine <-              0.9; ecf_small_intestine <-           0.0121
    vv_large_intestine <-          0.00874; vis_large_intestine <-           0.0953; vorg_large_intestine <-            0.548; plq_large_intestine <-           12.867; sigv_large_intestine <-             0.95; ecf_large_intestine <-           0.0209
    vv_pancreas <-           0.0057; vis_pancreas <-            0.018; vorg_pancreas <-            0.104; plq_pancreas <-            3.056; sigv_pancreas <-              0.9; ecf_pancreas <-           0.0001
    vv_thymus <-         0.000353; vis_thymus <-          0.00109; vorg_thymus <-          0.00641; plq_thymus <-            0.353; sigv_thymus <-              0.9; ecf_thymus <-           0.0001
    vv_spleen <-           0.0268; vis_spleen <-           0.0443; vorg_spleen <-            0.221; plq_spleen <-            6.343; sigv_spleen <-             0.85; ecf_spleen <-           0.0499
    vv_other <-            0.204; vis_other <-            0.831; vorg_other <-            4.852; plq_other <-             9.19; sigv_other <-             0.95; ecf_other <-           0.0951
    vplasma <- 3.126   # Table S2 plasma volume 3126 mL
    vln     <- 0.274   # Table S2 lymph-node volume 274 mL

    # Lymph flow is 0.2% of plasma flow; the liver additionally receives the
    # lymph-corrected flow of its four drained organs (Madonna LF[LIVER]).
    # Table S2 lists "Other" plasma flow as 5521 mL/h and lymph node as
    # 3670 mL/h; the deposited listing lumps the two into PLQ[Other] = 9190,
    # the lymph node being represented instead by L_LymphNode below.
    lf_lung <- plq_lung * 0.002
    lf_heart <- plq_heart * 0.002
    lf_muscle <- plq_muscle * 0.002
    lf_skin <- plq_skin * 0.002
    lf_adipose <- plq_adipose * 0.002
    lf_bone <- plq_bone * 0.002
    lf_brain <- plq_brain * 0.002
    lf_kidney <- plq_kidney * 0.002
    lf_small_intestine <- plq_small_intestine * 0.002
    lf_large_intestine <- plq_large_intestine * 0.002
    lf_pancreas <- plq_pancreas * 0.002
    lf_thymus <- plq_thymus * 0.002
    lf_spleen <- plq_spleen * 0.002
    lf_other <- plq_other * 0.002
    lf_liver <- (plq_liver + plq_small_intestine - lf_small_intestine + plq_large_intestine - lf_large_intestine + plq_spleen - lf_spleen + plq_pancreas - lf_pancreas) * 0.002
    llymph <- lf_lung + lf_liver + lf_heart + lf_muscle + lf_skin + lf_adipose + lf_bone + lf_brain + lf_kidney + lf_small_intestine + lf_large_intestine + lf_pancreas + lf_thymus + lf_spleen + lf_other

    # ---------------- endothelial cells, endosomal volumes ----------------
    # Organ endothelial-cell number uses each organ's share of total body
    # FcRn as a proxy (Methods); CL_up and the membrane/endosomal volumes
    # follow from it. Madonna "SIMULATION EQUATIONS" block.
    clupcell <- exp(lclup)
    pinotime <- exp(lpinotime)
    nendo    <- exp(lnendomouse) * exp(lendoscale)
    oec_lung <- nendo * ecf_lung
    oec_liver <- nendo * ecf_liver
    oec_heart <- nendo * ecf_heart
    oec_muscle <- nendo * ecf_muscle
    oec_skin <- nendo * ecf_skin
    oec_adipose <- nendo * ecf_adipose
    oec_bone <- nendo * ecf_bone
    oec_brain <- nendo * ecf_brain
    oec_kidney <- nendo * ecf_kidney
    oec_small_intestine <- nendo * ecf_small_intestine
    oec_large_intestine <- nendo * ecf_large_intestine
    oec_pancreas <- nendo * ecf_pancreas
    oec_thymus <- nendo * ecf_thymus
    oec_spleen <- nendo * ecf_spleen
    oec_other <- nendo * ecf_other
    clup_lung <- clupcell * oec_lung * 1e-15
    clup_liver <- clupcell * oec_liver * 1e-15
    clup_heart <- clupcell * oec_heart * 1e-15
    clup_muscle <- clupcell * oec_muscle * 1e-15
    clup_skin <- clupcell * oec_skin * 1e-15
    clup_adipose <- clupcell * oec_adipose * 1e-15
    clup_bone <- clupcell * oec_bone * 1e-15
    clup_brain <- clupcell * oec_brain * 1e-15
    clup_kidney <- clupcell * oec_kidney * 1e-15
    clup_small_intestine <- clupcell * oec_small_intestine * 1e-15
    clup_large_intestine <- clupcell * oec_large_intestine * 1e-15
    clup_pancreas <- clupcell * oec_pancreas * 1e-15
    clup_thymus <- clupcell * oec_thymus * 1e-15
    clup_spleen <- clupcell * oec_spleen * 1e-15
    clup_other <- clupcell * oec_other * 1e-15
    vendo_lung <- clupcell * pinotime * oec_lung * 1e-6 * 1e-9
    vendo_liver <- clupcell * pinotime * oec_liver * 1e-6 * 1e-9
    vendo_heart <- clupcell * pinotime * oec_heart * 1e-6 * 1e-9
    vendo_muscle <- clupcell * pinotime * oec_muscle * 1e-6 * 1e-9
    vendo_skin <- clupcell * pinotime * oec_skin * 1e-6 * 1e-9
    vendo_adipose <- clupcell * pinotime * oec_adipose * 1e-6 * 1e-9
    vendo_bone <- clupcell * pinotime * oec_bone * 1e-6 * 1e-9
    vendo_brain <- clupcell * pinotime * oec_brain * 1e-6 * 1e-9
    vendo_kidney <- clupcell * pinotime * oec_kidney * 1e-6 * 1e-9
    vendo_small_intestine <- clupcell * pinotime * oec_small_intestine * 1e-6 * 1e-9
    vendo_large_intestine <- clupcell * pinotime * oec_large_intestine * 1e-6 * 1e-9
    vendo_pancreas <- clupcell * pinotime * oec_pancreas * 1e-6 * 1e-9
    vendo_thymus <- clupcell * pinotime * oec_thymus * 1e-6 * 1e-9
    vendo_spleen <- clupcell * pinotime * oec_spleen * 1e-6 * 1e-9
    vendo_other <- clupcell * pinotime * oec_other * 1e-6 * 1e-9
    vendotot <- vendo_lung + vendo_liver + vendo_heart + vendo_muscle + vendo_skin + vendo_adipose + vendo_bone + vendo_brain + vendo_kidney + vendo_small_intestine + vendo_large_intestine + vendo_pancreas + vendo_thymus + vendo_spleen + vendo_other
    # FcRn is assumed to sit only in endothelial cells at one concentration
    # throughout the body: total amount / total endosomal volume.
    fcrnconc <- exp(lfcrntot) * 1e-3 / vendotot
    e7pct  <- (1 - e6apct) / 2
    e7bpct <- (1 - e6apct) / 2
    ve7_lung <- vendo_lung * e7pct; ve6a_lung <- vendo_lung * e6apct; ve7b_lung <- vendo_lung * e7bpct
    ve7_liver <- vendo_liver * e7pct; ve6a_liver <- vendo_liver * e6apct; ve7b_liver <- vendo_liver * e7bpct
    ve7_heart <- vendo_heart * e7pct; ve6a_heart <- vendo_heart * e6apct; ve7b_heart <- vendo_heart * e7bpct
    ve7_muscle <- vendo_muscle * e7pct; ve6a_muscle <- vendo_muscle * e6apct; ve7b_muscle <- vendo_muscle * e7bpct
    ve7_skin <- vendo_skin * e7pct; ve6a_skin <- vendo_skin * e6apct; ve7b_skin <- vendo_skin * e7bpct
    ve7_adipose <- vendo_adipose * e7pct; ve6a_adipose <- vendo_adipose * e6apct; ve7b_adipose <- vendo_adipose * e7bpct
    ve7_bone <- vendo_bone * e7pct; ve6a_bone <- vendo_bone * e6apct; ve7b_bone <- vendo_bone * e7bpct
    ve7_brain <- vendo_brain * e7pct; ve6a_brain <- vendo_brain * e6apct; ve7b_brain <- vendo_brain * e7bpct
    ve7_kidney <- vendo_kidney * e7pct; ve6a_kidney <- vendo_kidney * e6apct; ve7b_kidney <- vendo_kidney * e7bpct
    ve7_small_intestine <- vendo_small_intestine * e7pct; ve6a_small_intestine <- vendo_small_intestine * e6apct; ve7b_small_intestine <- vendo_small_intestine * e7bpct
    ve7_large_intestine <- vendo_large_intestine * e7pct; ve6a_large_intestine <- vendo_large_intestine * e6apct; ve7b_large_intestine <- vendo_large_intestine * e7bpct
    ve7_pancreas <- vendo_pancreas * e7pct; ve6a_pancreas <- vendo_pancreas * e6apct; ve7b_pancreas <- vendo_pancreas * e7bpct
    ve7_thymus <- vendo_thymus * e7pct; ve6a_thymus <- vendo_thymus * e6apct; ve7b_thymus <- vendo_thymus * e7bpct
    ve7_spleen <- vendo_spleen * e7pct; ve6a_spleen <- vendo_spleen * e6apct; ve7b_spleen <- vendo_spleen * e7bpct
    ve7_other <- vendo_other * e7pct; ve6a_other <- vendo_other * e6apct; ve7b_other <- vendo_other * e7bpct
    vvm_lung <- clup_lung * tauvm; vism_lung <- clup_lung * tauism
    vvm_liver <- clup_liver * tauvm; vism_liver <- clup_liver * tauism
    vvm_heart <- clup_heart * tauvm; vism_heart <- clup_heart * tauism
    vvm_muscle <- clup_muscle * tauvm; vism_muscle <- clup_muscle * tauism
    vvm_skin <- clup_skin * tauvm; vism_skin <- clup_skin * tauism
    vvm_adipose <- clup_adipose * tauvm; vism_adipose <- clup_adipose * tauism
    vvm_bone <- clup_bone * tauvm; vism_bone <- clup_bone * tauism
    vvm_brain <- clup_brain * tauvm; vism_brain <- clup_brain * tauism
    vvm_kidney <- clup_kidney * tauvm; vism_kidney <- clup_kidney * tauism
    vvm_small_intestine <- clup_small_intestine * tauvm; vism_small_intestine <- clup_small_intestine * tauism
    vvm_large_intestine <- clup_large_intestine * tauvm; vism_large_intestine <- clup_large_intestine * tauism
    vvm_pancreas <- clup_pancreas * tauvm; vism_pancreas <- clup_pancreas * tauism
    vvm_thymus <- clup_thymus * tauvm; vism_thymus <- clup_thymus * tauism
    vvm_spleen <- clup_spleen * tauvm; vism_spleen <- clup_spleen * tauism
    vvm_other <- clup_other * tauvm; vism_other <- clup_other * tauism

    # ---------------- initial conditions ----------------
    # FcRn starts at its steady-state endosomal concentration in the three
    # endosomal compartments and at 1e-4 of it on the two membranes
    # (Madonna INIT block). Free FcRn is not synthesised or degraded.
    fcrnendoearly_lung(0) <- fcrnconc * ve7_lung; fcrnendosort_lung(0) <- fcrnconc * ve6a_lung; fcrnendorecyc_lung(0) <- fcrnconc * ve7b_lung
    fcrnendoearly_liver(0) <- fcrnconc * ve7_liver; fcrnendosort_liver(0) <- fcrnconc * ve6a_liver; fcrnendorecyc_liver(0) <- fcrnconc * ve7b_liver
    fcrnendoearly_heart(0) <- fcrnconc * ve7_heart; fcrnendosort_heart(0) <- fcrnconc * ve6a_heart; fcrnendorecyc_heart(0) <- fcrnconc * ve7b_heart
    fcrnendoearly_muscle(0) <- fcrnconc * ve7_muscle; fcrnendosort_muscle(0) <- fcrnconc * ve6a_muscle; fcrnendorecyc_muscle(0) <- fcrnconc * ve7b_muscle
    fcrnendoearly_skin(0) <- fcrnconc * ve7_skin; fcrnendosort_skin(0) <- fcrnconc * ve6a_skin; fcrnendorecyc_skin(0) <- fcrnconc * ve7b_skin
    fcrnendoearly_adipose(0) <- fcrnconc * ve7_adipose; fcrnendosort_adipose(0) <- fcrnconc * ve6a_adipose; fcrnendorecyc_adipose(0) <- fcrnconc * ve7b_adipose
    fcrnendoearly_bone(0) <- fcrnconc * ve7_bone; fcrnendosort_bone(0) <- fcrnconc * ve6a_bone; fcrnendorecyc_bone(0) <- fcrnconc * ve7b_bone
    fcrnendoearly_brain(0) <- fcrnconc * ve7_brain; fcrnendosort_brain(0) <- fcrnconc * ve6a_brain; fcrnendorecyc_brain(0) <- fcrnconc * ve7b_brain
    fcrnendoearly_kidney(0) <- fcrnconc * ve7_kidney; fcrnendosort_kidney(0) <- fcrnconc * ve6a_kidney; fcrnendorecyc_kidney(0) <- fcrnconc * ve7b_kidney
    fcrnendoearly_small_intestine(0) <- fcrnconc * ve7_small_intestine; fcrnendosort_small_intestine(0) <- fcrnconc * ve6a_small_intestine; fcrnendorecyc_small_intestine(0) <- fcrnconc * ve7b_small_intestine
    fcrnendoearly_large_intestine(0) <- fcrnconc * ve7_large_intestine; fcrnendosort_large_intestine(0) <- fcrnconc * ve6a_large_intestine; fcrnendorecyc_large_intestine(0) <- fcrnconc * ve7b_large_intestine
    fcrnendoearly_pancreas(0) <- fcrnconc * ve7_pancreas; fcrnendosort_pancreas(0) <- fcrnconc * ve6a_pancreas; fcrnendorecyc_pancreas(0) <- fcrnconc * ve7b_pancreas
    fcrnendoearly_thymus(0) <- fcrnconc * ve7_thymus; fcrnendosort_thymus(0) <- fcrnconc * ve6a_thymus; fcrnendorecyc_thymus(0) <- fcrnconc * ve7b_thymus
    fcrnendoearly_spleen(0) <- fcrnconc * ve7_spleen; fcrnendosort_spleen(0) <- fcrnconc * ve6a_spleen; fcrnendorecyc_spleen(0) <- fcrnconc * ve7b_spleen
    fcrnendoearly_other(0) <- fcrnconc * ve7_other; fcrnendosort_other(0) <- fcrnconc * ve6a_other; fcrnendorecyc_other(0) <- fcrnconc * ve7b_other
    fcrnmemvas_lung(0) <- fcrnconc * 1e-4 * vvm_lung; fcrnmemint_lung(0) <- fcrnconc * 1e-4 * vism_lung
    fcrnmemvas_liver(0) <- fcrnconc * 1e-4 * vvm_liver; fcrnmemint_liver(0) <- fcrnconc * 1e-4 * vism_liver
    fcrnmemvas_heart(0) <- fcrnconc * 1e-4 * vvm_heart; fcrnmemint_heart(0) <- fcrnconc * 1e-4 * vism_heart
    fcrnmemvas_muscle(0) <- fcrnconc * 1e-4 * vvm_muscle; fcrnmemint_muscle(0) <- fcrnconc * 1e-4 * vism_muscle
    fcrnmemvas_skin(0) <- fcrnconc * 1e-4 * vvm_skin; fcrnmemint_skin(0) <- fcrnconc * 1e-4 * vism_skin
    fcrnmemvas_adipose(0) <- fcrnconc * 1e-4 * vvm_adipose; fcrnmemint_adipose(0) <- fcrnconc * 1e-4 * vism_adipose
    fcrnmemvas_bone(0) <- fcrnconc * 1e-4 * vvm_bone; fcrnmemint_bone(0) <- fcrnconc * 1e-4 * vism_bone
    fcrnmemvas_brain(0) <- fcrnconc * 1e-4 * vvm_brain; fcrnmemint_brain(0) <- fcrnconc * 1e-4 * vism_brain
    fcrnmemvas_kidney(0) <- fcrnconc * 1e-4 * vvm_kidney; fcrnmemint_kidney(0) <- fcrnconc * 1e-4 * vism_kidney
    fcrnmemvas_small_intestine(0) <- fcrnconc * 1e-4 * vvm_small_intestine; fcrnmemint_small_intestine(0) <- fcrnconc * 1e-4 * vism_small_intestine
    fcrnmemvas_large_intestine(0) <- fcrnconc * 1e-4 * vvm_large_intestine; fcrnmemint_large_intestine(0) <- fcrnconc * 1e-4 * vism_large_intestine
    fcrnmemvas_pancreas(0) <- fcrnconc * 1e-4 * vvm_pancreas; fcrnmemint_pancreas(0) <- fcrnconc * 1e-4 * vism_pancreas
    fcrnmemvas_thymus(0) <- fcrnconc * 1e-4 * vvm_thymus; fcrnmemint_thymus(0) <- fcrnconc * 1e-4 * vism_thymus
    fcrnmemvas_spleen(0) <- fcrnconc * 1e-4 * vvm_spleen; fcrnmemint_spleen(0) <- fcrnconc * 1e-4 * vism_spleen
    fcrnmemvas_other(0) <- fcrnconc * 1e-4 * vvm_other; fcrnmemint_other(0) <- fcrnconc * 1e-4 * vism_other

    # ---------------- concentrations (umol / L = uM) ----------------
    c_central <- central / vplasma
    c_lnode  <- lnode / vln
    c_lnode_igg  <- lnode_igg / vln
    # Endogenous IgG in plasma is clamped (Madonna d/dt(C_EDG_Plasma) = 0).
    c_igg_plasma <- cigg0
    c_vp_lung <- vp_lung / vv_lung
    c_vp_liver <- vp_liver / vv_liver
    c_vp_heart <- vp_heart / vv_heart
    c_vp_muscle <- vp_muscle / vv_muscle
    c_vp_skin <- vp_skin / vv_skin
    c_vp_adipose <- vp_adipose / vv_adipose
    c_vp_bone <- vp_bone / vv_bone
    c_vp_brain <- vp_brain / vv_brain
    c_vp_kidney <- vp_kidney / vv_kidney
    c_vp_small_intestine <- vp_small_intestine / vv_small_intestine
    c_vp_large_intestine <- vp_large_intestine / vv_large_intestine
    c_vp_pancreas <- vp_pancreas / vv_pancreas
    c_vp_thymus <- vp_thymus / vv_thymus
    c_vp_spleen <- vp_spleen / vv_spleen
    c_vp_other <- vp_other / vv_other
    c_memvas_lung <- memvas_lung / vvm_lung
    c_memvas_liver <- memvas_liver / vvm_liver
    c_memvas_heart <- memvas_heart / vvm_heart
    c_memvas_muscle <- memvas_muscle / vvm_muscle
    c_memvas_skin <- memvas_skin / vvm_skin
    c_memvas_adipose <- memvas_adipose / vvm_adipose
    c_memvas_bone <- memvas_bone / vvm_bone
    c_memvas_brain <- memvas_brain / vvm_brain
    c_memvas_kidney <- memvas_kidney / vvm_kidney
    c_memvas_small_intestine <- memvas_small_intestine / vvm_small_intestine
    c_memvas_large_intestine <- memvas_large_intestine / vvm_large_intestine
    c_memvas_pancreas <- memvas_pancreas / vvm_pancreas
    c_memvas_thymus <- memvas_thymus / vvm_thymus
    c_memvas_spleen <- memvas_spleen / vvm_spleen
    c_memvas_other <- memvas_other / vvm_other
    c_endoearly_lung <- endoearly_lung / ve7_lung
    c_endoearly_liver <- endoearly_liver / ve7_liver
    c_endoearly_heart <- endoearly_heart / ve7_heart
    c_endoearly_muscle <- endoearly_muscle / ve7_muscle
    c_endoearly_skin <- endoearly_skin / ve7_skin
    c_endoearly_adipose <- endoearly_adipose / ve7_adipose
    c_endoearly_bone <- endoearly_bone / ve7_bone
    c_endoearly_brain <- endoearly_brain / ve7_brain
    c_endoearly_kidney <- endoearly_kidney / ve7_kidney
    c_endoearly_small_intestine <- endoearly_small_intestine / ve7_small_intestine
    c_endoearly_large_intestine <- endoearly_large_intestine / ve7_large_intestine
    c_endoearly_pancreas <- endoearly_pancreas / ve7_pancreas
    c_endoearly_thymus <- endoearly_thymus / ve7_thymus
    c_endoearly_spleen <- endoearly_spleen / ve7_spleen
    c_endoearly_other <- endoearly_other / ve7_other
    c_endosort_lung <- endosort_lung / ve6a_lung
    c_endosort_liver <- endosort_liver / ve6a_liver
    c_endosort_heart <- endosort_heart / ve6a_heart
    c_endosort_muscle <- endosort_muscle / ve6a_muscle
    c_endosort_skin <- endosort_skin / ve6a_skin
    c_endosort_adipose <- endosort_adipose / ve6a_adipose
    c_endosort_bone <- endosort_bone / ve6a_bone
    c_endosort_brain <- endosort_brain / ve6a_brain
    c_endosort_kidney <- endosort_kidney / ve6a_kidney
    c_endosort_small_intestine <- endosort_small_intestine / ve6a_small_intestine
    c_endosort_large_intestine <- endosort_large_intestine / ve6a_large_intestine
    c_endosort_pancreas <- endosort_pancreas / ve6a_pancreas
    c_endosort_thymus <- endosort_thymus / ve6a_thymus
    c_endosort_spleen <- endosort_spleen / ve6a_spleen
    c_endosort_other <- endosort_other / ve6a_other
    c_endorecyc_lung <- endorecyc_lung / ve7b_lung
    c_endorecyc_liver <- endorecyc_liver / ve7b_liver
    c_endorecyc_heart <- endorecyc_heart / ve7b_heart
    c_endorecyc_muscle <- endorecyc_muscle / ve7b_muscle
    c_endorecyc_skin <- endorecyc_skin / ve7b_skin
    c_endorecyc_adipose <- endorecyc_adipose / ve7b_adipose
    c_endorecyc_bone <- endorecyc_bone / ve7b_bone
    c_endorecyc_brain <- endorecyc_brain / ve7b_brain
    c_endorecyc_kidney <- endorecyc_kidney / ve7b_kidney
    c_endorecyc_small_intestine <- endorecyc_small_intestine / ve7b_small_intestine
    c_endorecyc_large_intestine <- endorecyc_large_intestine / ve7b_large_intestine
    c_endorecyc_pancreas <- endorecyc_pancreas / ve7b_pancreas
    c_endorecyc_thymus <- endorecyc_thymus / ve7b_thymus
    c_endorecyc_spleen <- endorecyc_spleen / ve7b_spleen
    c_endorecyc_other <- endorecyc_other / ve7b_other
    c_memint_lung <- memint_lung / vism_lung
    c_memint_liver <- memint_liver / vism_liver
    c_memint_heart <- memint_heart / vism_heart
    c_memint_muscle <- memint_muscle / vism_muscle
    c_memint_skin <- memint_skin / vism_skin
    c_memint_adipose <- memint_adipose / vism_adipose
    c_memint_bone <- memint_bone / vism_bone
    c_memint_brain <- memint_brain / vism_brain
    c_memint_kidney <- memint_kidney / vism_kidney
    c_memint_small_intestine <- memint_small_intestine / vism_small_intestine
    c_memint_large_intestine <- memint_large_intestine / vism_large_intestine
    c_memint_pancreas <- memint_pancreas / vism_pancreas
    c_memint_thymus <- memint_thymus / vism_thymus
    c_memint_spleen <- memint_spleen / vism_spleen
    c_memint_other <- memint_other / vism_other
    c_is_lung <- is_lung / vis_lung
    c_is_liver <- is_liver / vis_liver
    c_is_heart <- is_heart / vis_heart
    c_is_muscle <- is_muscle / vis_muscle
    c_is_skin <- is_skin / vis_skin
    c_is_adipose <- is_adipose / vis_adipose
    c_is_bone <- is_bone / vis_bone
    c_is_brain <- is_brain / vis_brain
    c_is_kidney <- is_kidney / vis_kidney
    c_is_small_intestine <- is_small_intestine / vis_small_intestine
    c_is_large_intestine <- is_large_intestine / vis_large_intestine
    c_is_pancreas <- is_pancreas / vis_pancreas
    c_is_thymus <- is_thymus / vis_thymus
    c_is_spleen <- is_spleen / vis_spleen
    c_is_other <- is_other / vis_other
    c_memvasfr1_lung <- memvasfr1_lung / vvm_lung
    c_memvasfr1_liver <- memvasfr1_liver / vvm_liver
    c_memvasfr1_heart <- memvasfr1_heart / vvm_heart
    c_memvasfr1_muscle <- memvasfr1_muscle / vvm_muscle
    c_memvasfr1_skin <- memvasfr1_skin / vvm_skin
    c_memvasfr1_adipose <- memvasfr1_adipose / vvm_adipose
    c_memvasfr1_bone <- memvasfr1_bone / vvm_bone
    c_memvasfr1_brain <- memvasfr1_brain / vvm_brain
    c_memvasfr1_kidney <- memvasfr1_kidney / vvm_kidney
    c_memvasfr1_small_intestine <- memvasfr1_small_intestine / vvm_small_intestine
    c_memvasfr1_large_intestine <- memvasfr1_large_intestine / vvm_large_intestine
    c_memvasfr1_pancreas <- memvasfr1_pancreas / vvm_pancreas
    c_memvasfr1_thymus <- memvasfr1_thymus / vvm_thymus
    c_memvasfr1_spleen <- memvasfr1_spleen / vvm_spleen
    c_memvasfr1_other <- memvasfr1_other / vvm_other
    c_memvasfr2_lung <- memvasfr2_lung / vvm_lung
    c_memvasfr2_liver <- memvasfr2_liver / vvm_liver
    c_memvasfr2_heart <- memvasfr2_heart / vvm_heart
    c_memvasfr2_muscle <- memvasfr2_muscle / vvm_muscle
    c_memvasfr2_skin <- memvasfr2_skin / vvm_skin
    c_memvasfr2_adipose <- memvasfr2_adipose / vvm_adipose
    c_memvasfr2_bone <- memvasfr2_bone / vvm_bone
    c_memvasfr2_brain <- memvasfr2_brain / vvm_brain
    c_memvasfr2_kidney <- memvasfr2_kidney / vvm_kidney
    c_memvasfr2_small_intestine <- memvasfr2_small_intestine / vvm_small_intestine
    c_memvasfr2_large_intestine <- memvasfr2_large_intestine / vvm_large_intestine
    c_memvasfr2_pancreas <- memvasfr2_pancreas / vvm_pancreas
    c_memvasfr2_thymus <- memvasfr2_thymus / vvm_thymus
    c_memvasfr2_spleen <- memvasfr2_spleen / vvm_spleen
    c_memvasfr2_other <- memvasfr2_other / vvm_other
    c_endoearlyfr1_lung <- endoearlyfr1_lung / ve7_lung
    c_endoearlyfr1_liver <- endoearlyfr1_liver / ve7_liver
    c_endoearlyfr1_heart <- endoearlyfr1_heart / ve7_heart
    c_endoearlyfr1_muscle <- endoearlyfr1_muscle / ve7_muscle
    c_endoearlyfr1_skin <- endoearlyfr1_skin / ve7_skin
    c_endoearlyfr1_adipose <- endoearlyfr1_adipose / ve7_adipose
    c_endoearlyfr1_bone <- endoearlyfr1_bone / ve7_bone
    c_endoearlyfr1_brain <- endoearlyfr1_brain / ve7_brain
    c_endoearlyfr1_kidney <- endoearlyfr1_kidney / ve7_kidney
    c_endoearlyfr1_small_intestine <- endoearlyfr1_small_intestine / ve7_small_intestine
    c_endoearlyfr1_large_intestine <- endoearlyfr1_large_intestine / ve7_large_intestine
    c_endoearlyfr1_pancreas <- endoearlyfr1_pancreas / ve7_pancreas
    c_endoearlyfr1_thymus <- endoearlyfr1_thymus / ve7_thymus
    c_endoearlyfr1_spleen <- endoearlyfr1_spleen / ve7_spleen
    c_endoearlyfr1_other <- endoearlyfr1_other / ve7_other
    c_endoearlyfr2_lung <- endoearlyfr2_lung / ve7_lung
    c_endoearlyfr2_liver <- endoearlyfr2_liver / ve7_liver
    c_endoearlyfr2_heart <- endoearlyfr2_heart / ve7_heart
    c_endoearlyfr2_muscle <- endoearlyfr2_muscle / ve7_muscle
    c_endoearlyfr2_skin <- endoearlyfr2_skin / ve7_skin
    c_endoearlyfr2_adipose <- endoearlyfr2_adipose / ve7_adipose
    c_endoearlyfr2_bone <- endoearlyfr2_bone / ve7_bone
    c_endoearlyfr2_brain <- endoearlyfr2_brain / ve7_brain
    c_endoearlyfr2_kidney <- endoearlyfr2_kidney / ve7_kidney
    c_endoearlyfr2_small_intestine <- endoearlyfr2_small_intestine / ve7_small_intestine
    c_endoearlyfr2_large_intestine <- endoearlyfr2_large_intestine / ve7_large_intestine
    c_endoearlyfr2_pancreas <- endoearlyfr2_pancreas / ve7_pancreas
    c_endoearlyfr2_thymus <- endoearlyfr2_thymus / ve7_thymus
    c_endoearlyfr2_spleen <- endoearlyfr2_spleen / ve7_spleen
    c_endoearlyfr2_other <- endoearlyfr2_other / ve7_other
    c_endosortfr1_lung <- endosortfr1_lung / ve6a_lung
    c_endosortfr1_liver <- endosortfr1_liver / ve6a_liver
    c_endosortfr1_heart <- endosortfr1_heart / ve6a_heart
    c_endosortfr1_muscle <- endosortfr1_muscle / ve6a_muscle
    c_endosortfr1_skin <- endosortfr1_skin / ve6a_skin
    c_endosortfr1_adipose <- endosortfr1_adipose / ve6a_adipose
    c_endosortfr1_bone <- endosortfr1_bone / ve6a_bone
    c_endosortfr1_brain <- endosortfr1_brain / ve6a_brain
    c_endosortfr1_kidney <- endosortfr1_kidney / ve6a_kidney
    c_endosortfr1_small_intestine <- endosortfr1_small_intestine / ve6a_small_intestine
    c_endosortfr1_large_intestine <- endosortfr1_large_intestine / ve6a_large_intestine
    c_endosortfr1_pancreas <- endosortfr1_pancreas / ve6a_pancreas
    c_endosortfr1_thymus <- endosortfr1_thymus / ve6a_thymus
    c_endosortfr1_spleen <- endosortfr1_spleen / ve6a_spleen
    c_endosortfr1_other <- endosortfr1_other / ve6a_other
    c_endosortfr2_lung <- endosortfr2_lung / ve6a_lung
    c_endosortfr2_liver <- endosortfr2_liver / ve6a_liver
    c_endosortfr2_heart <- endosortfr2_heart / ve6a_heart
    c_endosortfr2_muscle <- endosortfr2_muscle / ve6a_muscle
    c_endosortfr2_skin <- endosortfr2_skin / ve6a_skin
    c_endosortfr2_adipose <- endosortfr2_adipose / ve6a_adipose
    c_endosortfr2_bone <- endosortfr2_bone / ve6a_bone
    c_endosortfr2_brain <- endosortfr2_brain / ve6a_brain
    c_endosortfr2_kidney <- endosortfr2_kidney / ve6a_kidney
    c_endosortfr2_small_intestine <- endosortfr2_small_intestine / ve6a_small_intestine
    c_endosortfr2_large_intestine <- endosortfr2_large_intestine / ve6a_large_intestine
    c_endosortfr2_pancreas <- endosortfr2_pancreas / ve6a_pancreas
    c_endosortfr2_thymus <- endosortfr2_thymus / ve6a_thymus
    c_endosortfr2_spleen <- endosortfr2_spleen / ve6a_spleen
    c_endosortfr2_other <- endosortfr2_other / ve6a_other
    c_endorecycfr1_lung <- endorecycfr1_lung / ve7b_lung
    c_endorecycfr1_liver <- endorecycfr1_liver / ve7b_liver
    c_endorecycfr1_heart <- endorecycfr1_heart / ve7b_heart
    c_endorecycfr1_muscle <- endorecycfr1_muscle / ve7b_muscle
    c_endorecycfr1_skin <- endorecycfr1_skin / ve7b_skin
    c_endorecycfr1_adipose <- endorecycfr1_adipose / ve7b_adipose
    c_endorecycfr1_bone <- endorecycfr1_bone / ve7b_bone
    c_endorecycfr1_brain <- endorecycfr1_brain / ve7b_brain
    c_endorecycfr1_kidney <- endorecycfr1_kidney / ve7b_kidney
    c_endorecycfr1_small_intestine <- endorecycfr1_small_intestine / ve7b_small_intestine
    c_endorecycfr1_large_intestine <- endorecycfr1_large_intestine / ve7b_large_intestine
    c_endorecycfr1_pancreas <- endorecycfr1_pancreas / ve7b_pancreas
    c_endorecycfr1_thymus <- endorecycfr1_thymus / ve7b_thymus
    c_endorecycfr1_spleen <- endorecycfr1_spleen / ve7b_spleen
    c_endorecycfr1_other <- endorecycfr1_other / ve7b_other
    c_endorecycfr2_lung <- endorecycfr2_lung / ve7b_lung
    c_endorecycfr2_liver <- endorecycfr2_liver / ve7b_liver
    c_endorecycfr2_heart <- endorecycfr2_heart / ve7b_heart
    c_endorecycfr2_muscle <- endorecycfr2_muscle / ve7b_muscle
    c_endorecycfr2_skin <- endorecycfr2_skin / ve7b_skin
    c_endorecycfr2_adipose <- endorecycfr2_adipose / ve7b_adipose
    c_endorecycfr2_bone <- endorecycfr2_bone / ve7b_bone
    c_endorecycfr2_brain <- endorecycfr2_brain / ve7b_brain
    c_endorecycfr2_kidney <- endorecycfr2_kidney / ve7b_kidney
    c_endorecycfr2_small_intestine <- endorecycfr2_small_intestine / ve7b_small_intestine
    c_endorecycfr2_large_intestine <- endorecycfr2_large_intestine / ve7b_large_intestine
    c_endorecycfr2_pancreas <- endorecycfr2_pancreas / ve7b_pancreas
    c_endorecycfr2_thymus <- endorecycfr2_thymus / ve7b_thymus
    c_endorecycfr2_spleen <- endorecycfr2_spleen / ve7b_spleen
    c_endorecycfr2_other <- endorecycfr2_other / ve7b_other
    c_memintfr1_lung <- memintfr1_lung / vism_lung
    c_memintfr1_liver <- memintfr1_liver / vism_liver
    c_memintfr1_heart <- memintfr1_heart / vism_heart
    c_memintfr1_muscle <- memintfr1_muscle / vism_muscle
    c_memintfr1_skin <- memintfr1_skin / vism_skin
    c_memintfr1_adipose <- memintfr1_adipose / vism_adipose
    c_memintfr1_bone <- memintfr1_bone / vism_bone
    c_memintfr1_brain <- memintfr1_brain / vism_brain
    c_memintfr1_kidney <- memintfr1_kidney / vism_kidney
    c_memintfr1_small_intestine <- memintfr1_small_intestine / vism_small_intestine
    c_memintfr1_large_intestine <- memintfr1_large_intestine / vism_large_intestine
    c_memintfr1_pancreas <- memintfr1_pancreas / vism_pancreas
    c_memintfr1_thymus <- memintfr1_thymus / vism_thymus
    c_memintfr1_spleen <- memintfr1_spleen / vism_spleen
    c_memintfr1_other <- memintfr1_other / vism_other
    c_memintfr2_lung <- memintfr2_lung / vism_lung
    c_memintfr2_liver <- memintfr2_liver / vism_liver
    c_memintfr2_heart <- memintfr2_heart / vism_heart
    c_memintfr2_muscle <- memintfr2_muscle / vism_muscle
    c_memintfr2_skin <- memintfr2_skin / vism_skin
    c_memintfr2_adipose <- memintfr2_adipose / vism_adipose
    c_memintfr2_bone <- memintfr2_bone / vism_bone
    c_memintfr2_brain <- memintfr2_brain / vism_brain
    c_memintfr2_kidney <- memintfr2_kidney / vism_kidney
    c_memintfr2_small_intestine <- memintfr2_small_intestine / vism_small_intestine
    c_memintfr2_large_intestine <- memintfr2_large_intestine / vism_large_intestine
    c_memintfr2_pancreas <- memintfr2_pancreas / vism_pancreas
    c_memintfr2_thymus <- memintfr2_thymus / vism_thymus
    c_memintfr2_spleen <- memintfr2_spleen / vism_spleen
    c_memintfr2_other <- memintfr2_other / vism_other
    c_memvasns_lung <- memvasns_lung / vvm_lung
    c_memvasns_liver <- memvasns_liver / vvm_liver
    c_memvasns_heart <- memvasns_heart / vvm_heart
    c_memvasns_muscle <- memvasns_muscle / vvm_muscle
    c_memvasns_skin <- memvasns_skin / vvm_skin
    c_memvasns_adipose <- memvasns_adipose / vvm_adipose
    c_memvasns_bone <- memvasns_bone / vvm_bone
    c_memvasns_brain <- memvasns_brain / vvm_brain
    c_memvasns_kidney <- memvasns_kidney / vvm_kidney
    c_memvasns_small_intestine <- memvasns_small_intestine / vvm_small_intestine
    c_memvasns_large_intestine <- memvasns_large_intestine / vvm_large_intestine
    c_memvasns_pancreas <- memvasns_pancreas / vvm_pancreas
    c_memvasns_thymus <- memvasns_thymus / vvm_thymus
    c_memvasns_spleen <- memvasns_spleen / vvm_spleen
    c_memvasns_other <- memvasns_other / vvm_other
    c_memintns_lung <- memintns_lung / vism_lung
    c_memintns_liver <- memintns_liver / vism_liver
    c_memintns_heart <- memintns_heart / vism_heart
    c_memintns_muscle <- memintns_muscle / vism_muscle
    c_memintns_skin <- memintns_skin / vism_skin
    c_memintns_adipose <- memintns_adipose / vism_adipose
    c_memintns_bone <- memintns_bone / vism_bone
    c_memintns_brain <- memintns_brain / vism_brain
    c_memintns_kidney <- memintns_kidney / vism_kidney
    c_memintns_small_intestine <- memintns_small_intestine / vism_small_intestine
    c_memintns_large_intestine <- memintns_large_intestine / vism_large_intestine
    c_memintns_pancreas <- memintns_pancreas / vism_pancreas
    c_memintns_thymus <- memintns_thymus / vism_thymus
    c_memintns_spleen <- memintns_spleen / vism_spleen
    c_memintns_other <- memintns_other / vism_other
    c_vp_lung_igg <- vp_lung_igg / vv_lung
    c_vp_liver_igg <- vp_liver_igg / vv_liver
    c_vp_heart_igg <- vp_heart_igg / vv_heart
    c_vp_muscle_igg <- vp_muscle_igg / vv_muscle
    c_vp_skin_igg <- vp_skin_igg / vv_skin
    c_vp_adipose_igg <- vp_adipose_igg / vv_adipose
    c_vp_bone_igg <- vp_bone_igg / vv_bone
    c_vp_brain_igg <- vp_brain_igg / vv_brain
    c_vp_kidney_igg <- vp_kidney_igg / vv_kidney
    c_vp_small_intestine_igg <- vp_small_intestine_igg / vv_small_intestine
    c_vp_large_intestine_igg <- vp_large_intestine_igg / vv_large_intestine
    c_vp_pancreas_igg <- vp_pancreas_igg / vv_pancreas
    c_vp_thymus_igg <- vp_thymus_igg / vv_thymus
    c_vp_spleen_igg <- vp_spleen_igg / vv_spleen
    c_vp_other_igg <- vp_other_igg / vv_other
    c_memvas_lung_igg <- memvas_lung_igg / vvm_lung
    c_memvas_liver_igg <- memvas_liver_igg / vvm_liver
    c_memvas_heart_igg <- memvas_heart_igg / vvm_heart
    c_memvas_muscle_igg <- memvas_muscle_igg / vvm_muscle
    c_memvas_skin_igg <- memvas_skin_igg / vvm_skin
    c_memvas_adipose_igg <- memvas_adipose_igg / vvm_adipose
    c_memvas_bone_igg <- memvas_bone_igg / vvm_bone
    c_memvas_brain_igg <- memvas_brain_igg / vvm_brain
    c_memvas_kidney_igg <- memvas_kidney_igg / vvm_kidney
    c_memvas_small_intestine_igg <- memvas_small_intestine_igg / vvm_small_intestine
    c_memvas_large_intestine_igg <- memvas_large_intestine_igg / vvm_large_intestine
    c_memvas_pancreas_igg <- memvas_pancreas_igg / vvm_pancreas
    c_memvas_thymus_igg <- memvas_thymus_igg / vvm_thymus
    c_memvas_spleen_igg <- memvas_spleen_igg / vvm_spleen
    c_memvas_other_igg <- memvas_other_igg / vvm_other
    c_endoearly_lung_igg <- endoearly_lung_igg / ve7_lung
    c_endoearly_liver_igg <- endoearly_liver_igg / ve7_liver
    c_endoearly_heart_igg <- endoearly_heart_igg / ve7_heart
    c_endoearly_muscle_igg <- endoearly_muscle_igg / ve7_muscle
    c_endoearly_skin_igg <- endoearly_skin_igg / ve7_skin
    c_endoearly_adipose_igg <- endoearly_adipose_igg / ve7_adipose
    c_endoearly_bone_igg <- endoearly_bone_igg / ve7_bone
    c_endoearly_brain_igg <- endoearly_brain_igg / ve7_brain
    c_endoearly_kidney_igg <- endoearly_kidney_igg / ve7_kidney
    c_endoearly_small_intestine_igg <- endoearly_small_intestine_igg / ve7_small_intestine
    c_endoearly_large_intestine_igg <- endoearly_large_intestine_igg / ve7_large_intestine
    c_endoearly_pancreas_igg <- endoearly_pancreas_igg / ve7_pancreas
    c_endoearly_thymus_igg <- endoearly_thymus_igg / ve7_thymus
    c_endoearly_spleen_igg <- endoearly_spleen_igg / ve7_spleen
    c_endoearly_other_igg <- endoearly_other_igg / ve7_other
    c_endosort_lung_igg <- endosort_lung_igg / ve6a_lung
    c_endosort_liver_igg <- endosort_liver_igg / ve6a_liver
    c_endosort_heart_igg <- endosort_heart_igg / ve6a_heart
    c_endosort_muscle_igg <- endosort_muscle_igg / ve6a_muscle
    c_endosort_skin_igg <- endosort_skin_igg / ve6a_skin
    c_endosort_adipose_igg <- endosort_adipose_igg / ve6a_adipose
    c_endosort_bone_igg <- endosort_bone_igg / ve6a_bone
    c_endosort_brain_igg <- endosort_brain_igg / ve6a_brain
    c_endosort_kidney_igg <- endosort_kidney_igg / ve6a_kidney
    c_endosort_small_intestine_igg <- endosort_small_intestine_igg / ve6a_small_intestine
    c_endosort_large_intestine_igg <- endosort_large_intestine_igg / ve6a_large_intestine
    c_endosort_pancreas_igg <- endosort_pancreas_igg / ve6a_pancreas
    c_endosort_thymus_igg <- endosort_thymus_igg / ve6a_thymus
    c_endosort_spleen_igg <- endosort_spleen_igg / ve6a_spleen
    c_endosort_other_igg <- endosort_other_igg / ve6a_other
    c_endorecyc_lung_igg <- endorecyc_lung_igg / ve7b_lung
    c_endorecyc_liver_igg <- endorecyc_liver_igg / ve7b_liver
    c_endorecyc_heart_igg <- endorecyc_heart_igg / ve7b_heart
    c_endorecyc_muscle_igg <- endorecyc_muscle_igg / ve7b_muscle
    c_endorecyc_skin_igg <- endorecyc_skin_igg / ve7b_skin
    c_endorecyc_adipose_igg <- endorecyc_adipose_igg / ve7b_adipose
    c_endorecyc_bone_igg <- endorecyc_bone_igg / ve7b_bone
    c_endorecyc_brain_igg <- endorecyc_brain_igg / ve7b_brain
    c_endorecyc_kidney_igg <- endorecyc_kidney_igg / ve7b_kidney
    c_endorecyc_small_intestine_igg <- endorecyc_small_intestine_igg / ve7b_small_intestine
    c_endorecyc_large_intestine_igg <- endorecyc_large_intestine_igg / ve7b_large_intestine
    c_endorecyc_pancreas_igg <- endorecyc_pancreas_igg / ve7b_pancreas
    c_endorecyc_thymus_igg <- endorecyc_thymus_igg / ve7b_thymus
    c_endorecyc_spleen_igg <- endorecyc_spleen_igg / ve7b_spleen
    c_endorecyc_other_igg <- endorecyc_other_igg / ve7b_other
    c_memint_lung_igg <- memint_lung_igg / vism_lung
    c_memint_liver_igg <- memint_liver_igg / vism_liver
    c_memint_heart_igg <- memint_heart_igg / vism_heart
    c_memint_muscle_igg <- memint_muscle_igg / vism_muscle
    c_memint_skin_igg <- memint_skin_igg / vism_skin
    c_memint_adipose_igg <- memint_adipose_igg / vism_adipose
    c_memint_bone_igg <- memint_bone_igg / vism_bone
    c_memint_brain_igg <- memint_brain_igg / vism_brain
    c_memint_kidney_igg <- memint_kidney_igg / vism_kidney
    c_memint_small_intestine_igg <- memint_small_intestine_igg / vism_small_intestine
    c_memint_large_intestine_igg <- memint_large_intestine_igg / vism_large_intestine
    c_memint_pancreas_igg <- memint_pancreas_igg / vism_pancreas
    c_memint_thymus_igg <- memint_thymus_igg / vism_thymus
    c_memint_spleen_igg <- memint_spleen_igg / vism_spleen
    c_memint_other_igg <- memint_other_igg / vism_other
    c_is_lung_igg <- is_lung_igg / vis_lung
    c_is_liver_igg <- is_liver_igg / vis_liver
    c_is_heart_igg <- is_heart_igg / vis_heart
    c_is_muscle_igg <- is_muscle_igg / vis_muscle
    c_is_skin_igg <- is_skin_igg / vis_skin
    c_is_adipose_igg <- is_adipose_igg / vis_adipose
    c_is_bone_igg <- is_bone_igg / vis_bone
    c_is_brain_igg <- is_brain_igg / vis_brain
    c_is_kidney_igg <- is_kidney_igg / vis_kidney
    c_is_small_intestine_igg <- is_small_intestine_igg / vis_small_intestine
    c_is_large_intestine_igg <- is_large_intestine_igg / vis_large_intestine
    c_is_pancreas_igg <- is_pancreas_igg / vis_pancreas
    c_is_thymus_igg <- is_thymus_igg / vis_thymus
    c_is_spleen_igg <- is_spleen_igg / vis_spleen
    c_is_other_igg <- is_other_igg / vis_other
    c_memvasfr1_lung_igg <- memvasfr1_lung_igg / vvm_lung
    c_memvasfr1_liver_igg <- memvasfr1_liver_igg / vvm_liver
    c_memvasfr1_heart_igg <- memvasfr1_heart_igg / vvm_heart
    c_memvasfr1_muscle_igg <- memvasfr1_muscle_igg / vvm_muscle
    c_memvasfr1_skin_igg <- memvasfr1_skin_igg / vvm_skin
    c_memvasfr1_adipose_igg <- memvasfr1_adipose_igg / vvm_adipose
    c_memvasfr1_bone_igg <- memvasfr1_bone_igg / vvm_bone
    c_memvasfr1_brain_igg <- memvasfr1_brain_igg / vvm_brain
    c_memvasfr1_kidney_igg <- memvasfr1_kidney_igg / vvm_kidney
    c_memvasfr1_small_intestine_igg <- memvasfr1_small_intestine_igg / vvm_small_intestine
    c_memvasfr1_large_intestine_igg <- memvasfr1_large_intestine_igg / vvm_large_intestine
    c_memvasfr1_pancreas_igg <- memvasfr1_pancreas_igg / vvm_pancreas
    c_memvasfr1_thymus_igg <- memvasfr1_thymus_igg / vvm_thymus
    c_memvasfr1_spleen_igg <- memvasfr1_spleen_igg / vvm_spleen
    c_memvasfr1_other_igg <- memvasfr1_other_igg / vvm_other
    c_memvasfr2_lung_igg <- memvasfr2_lung_igg / vvm_lung
    c_memvasfr2_liver_igg <- memvasfr2_liver_igg / vvm_liver
    c_memvasfr2_heart_igg <- memvasfr2_heart_igg / vvm_heart
    c_memvasfr2_muscle_igg <- memvasfr2_muscle_igg / vvm_muscle
    c_memvasfr2_skin_igg <- memvasfr2_skin_igg / vvm_skin
    c_memvasfr2_adipose_igg <- memvasfr2_adipose_igg / vvm_adipose
    c_memvasfr2_bone_igg <- memvasfr2_bone_igg / vvm_bone
    c_memvasfr2_brain_igg <- memvasfr2_brain_igg / vvm_brain
    c_memvasfr2_kidney_igg <- memvasfr2_kidney_igg / vvm_kidney
    c_memvasfr2_small_intestine_igg <- memvasfr2_small_intestine_igg / vvm_small_intestine
    c_memvasfr2_large_intestine_igg <- memvasfr2_large_intestine_igg / vvm_large_intestine
    c_memvasfr2_pancreas_igg <- memvasfr2_pancreas_igg / vvm_pancreas
    c_memvasfr2_thymus_igg <- memvasfr2_thymus_igg / vvm_thymus
    c_memvasfr2_spleen_igg <- memvasfr2_spleen_igg / vvm_spleen
    c_memvasfr2_other_igg <- memvasfr2_other_igg / vvm_other
    c_endoearlyfr1_lung_igg <- endoearlyfr1_lung_igg / ve7_lung
    c_endoearlyfr1_liver_igg <- endoearlyfr1_liver_igg / ve7_liver
    c_endoearlyfr1_heart_igg <- endoearlyfr1_heart_igg / ve7_heart
    c_endoearlyfr1_muscle_igg <- endoearlyfr1_muscle_igg / ve7_muscle
    c_endoearlyfr1_skin_igg <- endoearlyfr1_skin_igg / ve7_skin
    c_endoearlyfr1_adipose_igg <- endoearlyfr1_adipose_igg / ve7_adipose
    c_endoearlyfr1_bone_igg <- endoearlyfr1_bone_igg / ve7_bone
    c_endoearlyfr1_brain_igg <- endoearlyfr1_brain_igg / ve7_brain
    c_endoearlyfr1_kidney_igg <- endoearlyfr1_kidney_igg / ve7_kidney
    c_endoearlyfr1_small_intestine_igg <- endoearlyfr1_small_intestine_igg / ve7_small_intestine
    c_endoearlyfr1_large_intestine_igg <- endoearlyfr1_large_intestine_igg / ve7_large_intestine
    c_endoearlyfr1_pancreas_igg <- endoearlyfr1_pancreas_igg / ve7_pancreas
    c_endoearlyfr1_thymus_igg <- endoearlyfr1_thymus_igg / ve7_thymus
    c_endoearlyfr1_spleen_igg <- endoearlyfr1_spleen_igg / ve7_spleen
    c_endoearlyfr1_other_igg <- endoearlyfr1_other_igg / ve7_other
    c_endoearlyfr2_lung_igg <- endoearlyfr2_lung_igg / ve7_lung
    c_endoearlyfr2_liver_igg <- endoearlyfr2_liver_igg / ve7_liver
    c_endoearlyfr2_heart_igg <- endoearlyfr2_heart_igg / ve7_heart
    c_endoearlyfr2_muscle_igg <- endoearlyfr2_muscle_igg / ve7_muscle
    c_endoearlyfr2_skin_igg <- endoearlyfr2_skin_igg / ve7_skin
    c_endoearlyfr2_adipose_igg <- endoearlyfr2_adipose_igg / ve7_adipose
    c_endoearlyfr2_bone_igg <- endoearlyfr2_bone_igg / ve7_bone
    c_endoearlyfr2_brain_igg <- endoearlyfr2_brain_igg / ve7_brain
    c_endoearlyfr2_kidney_igg <- endoearlyfr2_kidney_igg / ve7_kidney
    c_endoearlyfr2_small_intestine_igg <- endoearlyfr2_small_intestine_igg / ve7_small_intestine
    c_endoearlyfr2_large_intestine_igg <- endoearlyfr2_large_intestine_igg / ve7_large_intestine
    c_endoearlyfr2_pancreas_igg <- endoearlyfr2_pancreas_igg / ve7_pancreas
    c_endoearlyfr2_thymus_igg <- endoearlyfr2_thymus_igg / ve7_thymus
    c_endoearlyfr2_spleen_igg <- endoearlyfr2_spleen_igg / ve7_spleen
    c_endoearlyfr2_other_igg <- endoearlyfr2_other_igg / ve7_other
    c_endosortfr1_lung_igg <- endosortfr1_lung_igg / ve6a_lung
    c_endosortfr1_liver_igg <- endosortfr1_liver_igg / ve6a_liver
    c_endosortfr1_heart_igg <- endosortfr1_heart_igg / ve6a_heart
    c_endosortfr1_muscle_igg <- endosortfr1_muscle_igg / ve6a_muscle
    c_endosortfr1_skin_igg <- endosortfr1_skin_igg / ve6a_skin
    c_endosortfr1_adipose_igg <- endosortfr1_adipose_igg / ve6a_adipose
    c_endosortfr1_bone_igg <- endosortfr1_bone_igg / ve6a_bone
    c_endosortfr1_brain_igg <- endosortfr1_brain_igg / ve6a_brain
    c_endosortfr1_kidney_igg <- endosortfr1_kidney_igg / ve6a_kidney
    c_endosortfr1_small_intestine_igg <- endosortfr1_small_intestine_igg / ve6a_small_intestine
    c_endosortfr1_large_intestine_igg <- endosortfr1_large_intestine_igg / ve6a_large_intestine
    c_endosortfr1_pancreas_igg <- endosortfr1_pancreas_igg / ve6a_pancreas
    c_endosortfr1_thymus_igg <- endosortfr1_thymus_igg / ve6a_thymus
    c_endosortfr1_spleen_igg <- endosortfr1_spleen_igg / ve6a_spleen
    c_endosortfr1_other_igg <- endosortfr1_other_igg / ve6a_other
    c_endosortfr2_lung_igg <- endosortfr2_lung_igg / ve6a_lung
    c_endosortfr2_liver_igg <- endosortfr2_liver_igg / ve6a_liver
    c_endosortfr2_heart_igg <- endosortfr2_heart_igg / ve6a_heart
    c_endosortfr2_muscle_igg <- endosortfr2_muscle_igg / ve6a_muscle
    c_endosortfr2_skin_igg <- endosortfr2_skin_igg / ve6a_skin
    c_endosortfr2_adipose_igg <- endosortfr2_adipose_igg / ve6a_adipose
    c_endosortfr2_bone_igg <- endosortfr2_bone_igg / ve6a_bone
    c_endosortfr2_brain_igg <- endosortfr2_brain_igg / ve6a_brain
    c_endosortfr2_kidney_igg <- endosortfr2_kidney_igg / ve6a_kidney
    c_endosortfr2_small_intestine_igg <- endosortfr2_small_intestine_igg / ve6a_small_intestine
    c_endosortfr2_large_intestine_igg <- endosortfr2_large_intestine_igg / ve6a_large_intestine
    c_endosortfr2_pancreas_igg <- endosortfr2_pancreas_igg / ve6a_pancreas
    c_endosortfr2_thymus_igg <- endosortfr2_thymus_igg / ve6a_thymus
    c_endosortfr2_spleen_igg <- endosortfr2_spleen_igg / ve6a_spleen
    c_endosortfr2_other_igg <- endosortfr2_other_igg / ve6a_other
    c_endorecycfr1_lung_igg <- endorecycfr1_lung_igg / ve7b_lung
    c_endorecycfr1_liver_igg <- endorecycfr1_liver_igg / ve7b_liver
    c_endorecycfr1_heart_igg <- endorecycfr1_heart_igg / ve7b_heart
    c_endorecycfr1_muscle_igg <- endorecycfr1_muscle_igg / ve7b_muscle
    c_endorecycfr1_skin_igg <- endorecycfr1_skin_igg / ve7b_skin
    c_endorecycfr1_adipose_igg <- endorecycfr1_adipose_igg / ve7b_adipose
    c_endorecycfr1_bone_igg <- endorecycfr1_bone_igg / ve7b_bone
    c_endorecycfr1_brain_igg <- endorecycfr1_brain_igg / ve7b_brain
    c_endorecycfr1_kidney_igg <- endorecycfr1_kidney_igg / ve7b_kidney
    c_endorecycfr1_small_intestine_igg <- endorecycfr1_small_intestine_igg / ve7b_small_intestine
    c_endorecycfr1_large_intestine_igg <- endorecycfr1_large_intestine_igg / ve7b_large_intestine
    c_endorecycfr1_pancreas_igg <- endorecycfr1_pancreas_igg / ve7b_pancreas
    c_endorecycfr1_thymus_igg <- endorecycfr1_thymus_igg / ve7b_thymus
    c_endorecycfr1_spleen_igg <- endorecycfr1_spleen_igg / ve7b_spleen
    c_endorecycfr1_other_igg <- endorecycfr1_other_igg / ve7b_other
    c_endorecycfr2_lung_igg <- endorecycfr2_lung_igg / ve7b_lung
    c_endorecycfr2_liver_igg <- endorecycfr2_liver_igg / ve7b_liver
    c_endorecycfr2_heart_igg <- endorecycfr2_heart_igg / ve7b_heart
    c_endorecycfr2_muscle_igg <- endorecycfr2_muscle_igg / ve7b_muscle
    c_endorecycfr2_skin_igg <- endorecycfr2_skin_igg / ve7b_skin
    c_endorecycfr2_adipose_igg <- endorecycfr2_adipose_igg / ve7b_adipose
    c_endorecycfr2_bone_igg <- endorecycfr2_bone_igg / ve7b_bone
    c_endorecycfr2_brain_igg <- endorecycfr2_brain_igg / ve7b_brain
    c_endorecycfr2_kidney_igg <- endorecycfr2_kidney_igg / ve7b_kidney
    c_endorecycfr2_small_intestine_igg <- endorecycfr2_small_intestine_igg / ve7b_small_intestine
    c_endorecycfr2_large_intestine_igg <- endorecycfr2_large_intestine_igg / ve7b_large_intestine
    c_endorecycfr2_pancreas_igg <- endorecycfr2_pancreas_igg / ve7b_pancreas
    c_endorecycfr2_thymus_igg <- endorecycfr2_thymus_igg / ve7b_thymus
    c_endorecycfr2_spleen_igg <- endorecycfr2_spleen_igg / ve7b_spleen
    c_endorecycfr2_other_igg <- endorecycfr2_other_igg / ve7b_other
    c_memintfr1_lung_igg <- memintfr1_lung_igg / vism_lung
    c_memintfr1_liver_igg <- memintfr1_liver_igg / vism_liver
    c_memintfr1_heart_igg <- memintfr1_heart_igg / vism_heart
    c_memintfr1_muscle_igg <- memintfr1_muscle_igg / vism_muscle
    c_memintfr1_skin_igg <- memintfr1_skin_igg / vism_skin
    c_memintfr1_adipose_igg <- memintfr1_adipose_igg / vism_adipose
    c_memintfr1_bone_igg <- memintfr1_bone_igg / vism_bone
    c_memintfr1_brain_igg <- memintfr1_brain_igg / vism_brain
    c_memintfr1_kidney_igg <- memintfr1_kidney_igg / vism_kidney
    c_memintfr1_small_intestine_igg <- memintfr1_small_intestine_igg / vism_small_intestine
    c_memintfr1_large_intestine_igg <- memintfr1_large_intestine_igg / vism_large_intestine
    c_memintfr1_pancreas_igg <- memintfr1_pancreas_igg / vism_pancreas
    c_memintfr1_thymus_igg <- memintfr1_thymus_igg / vism_thymus
    c_memintfr1_spleen_igg <- memintfr1_spleen_igg / vism_spleen
    c_memintfr1_other_igg <- memintfr1_other_igg / vism_other
    c_memintfr2_lung_igg <- memintfr2_lung_igg / vism_lung
    c_memintfr2_liver_igg <- memintfr2_liver_igg / vism_liver
    c_memintfr2_heart_igg <- memintfr2_heart_igg / vism_heart
    c_memintfr2_muscle_igg <- memintfr2_muscle_igg / vism_muscle
    c_memintfr2_skin_igg <- memintfr2_skin_igg / vism_skin
    c_memintfr2_adipose_igg <- memintfr2_adipose_igg / vism_adipose
    c_memintfr2_bone_igg <- memintfr2_bone_igg / vism_bone
    c_memintfr2_brain_igg <- memintfr2_brain_igg / vism_brain
    c_memintfr2_kidney_igg <- memintfr2_kidney_igg / vism_kidney
    c_memintfr2_small_intestine_igg <- memintfr2_small_intestine_igg / vism_small_intestine
    c_memintfr2_large_intestine_igg <- memintfr2_large_intestine_igg / vism_large_intestine
    c_memintfr2_pancreas_igg <- memintfr2_pancreas_igg / vism_pancreas
    c_memintfr2_thymus_igg <- memintfr2_thymus_igg / vism_thymus
    c_memintfr2_spleen_igg <- memintfr2_spleen_igg / vism_spleen
    c_memintfr2_other_igg <- memintfr2_other_igg / vism_other
    c_fcrnmemvas_lung <- fcrnmemvas_lung / vvm_lung
    c_fcrnmemvas_liver <- fcrnmemvas_liver / vvm_liver
    c_fcrnmemvas_heart <- fcrnmemvas_heart / vvm_heart
    c_fcrnmemvas_muscle <- fcrnmemvas_muscle / vvm_muscle
    c_fcrnmemvas_skin <- fcrnmemvas_skin / vvm_skin
    c_fcrnmemvas_adipose <- fcrnmemvas_adipose / vvm_adipose
    c_fcrnmemvas_bone <- fcrnmemvas_bone / vvm_bone
    c_fcrnmemvas_brain <- fcrnmemvas_brain / vvm_brain
    c_fcrnmemvas_kidney <- fcrnmemvas_kidney / vvm_kidney
    c_fcrnmemvas_small_intestine <- fcrnmemvas_small_intestine / vvm_small_intestine
    c_fcrnmemvas_large_intestine <- fcrnmemvas_large_intestine / vvm_large_intestine
    c_fcrnmemvas_pancreas <- fcrnmemvas_pancreas / vvm_pancreas
    c_fcrnmemvas_thymus <- fcrnmemvas_thymus / vvm_thymus
    c_fcrnmemvas_spleen <- fcrnmemvas_spleen / vvm_spleen
    c_fcrnmemvas_other <- fcrnmemvas_other / vvm_other
    c_fcrnendoearly_lung <- fcrnendoearly_lung / ve7_lung
    c_fcrnendoearly_liver <- fcrnendoearly_liver / ve7_liver
    c_fcrnendoearly_heart <- fcrnendoearly_heart / ve7_heart
    c_fcrnendoearly_muscle <- fcrnendoearly_muscle / ve7_muscle
    c_fcrnendoearly_skin <- fcrnendoearly_skin / ve7_skin
    c_fcrnendoearly_adipose <- fcrnendoearly_adipose / ve7_adipose
    c_fcrnendoearly_bone <- fcrnendoearly_bone / ve7_bone
    c_fcrnendoearly_brain <- fcrnendoearly_brain / ve7_brain
    c_fcrnendoearly_kidney <- fcrnendoearly_kidney / ve7_kidney
    c_fcrnendoearly_small_intestine <- fcrnendoearly_small_intestine / ve7_small_intestine
    c_fcrnendoearly_large_intestine <- fcrnendoearly_large_intestine / ve7_large_intestine
    c_fcrnendoearly_pancreas <- fcrnendoearly_pancreas / ve7_pancreas
    c_fcrnendoearly_thymus <- fcrnendoearly_thymus / ve7_thymus
    c_fcrnendoearly_spleen <- fcrnendoearly_spleen / ve7_spleen
    c_fcrnendoearly_other <- fcrnendoearly_other / ve7_other
    c_fcrnendosort_lung <- fcrnendosort_lung / ve6a_lung
    c_fcrnendosort_liver <- fcrnendosort_liver / ve6a_liver
    c_fcrnendosort_heart <- fcrnendosort_heart / ve6a_heart
    c_fcrnendosort_muscle <- fcrnendosort_muscle / ve6a_muscle
    c_fcrnendosort_skin <- fcrnendosort_skin / ve6a_skin
    c_fcrnendosort_adipose <- fcrnendosort_adipose / ve6a_adipose
    c_fcrnendosort_bone <- fcrnendosort_bone / ve6a_bone
    c_fcrnendosort_brain <- fcrnendosort_brain / ve6a_brain
    c_fcrnendosort_kidney <- fcrnendosort_kidney / ve6a_kidney
    c_fcrnendosort_small_intestine <- fcrnendosort_small_intestine / ve6a_small_intestine
    c_fcrnendosort_large_intestine <- fcrnendosort_large_intestine / ve6a_large_intestine
    c_fcrnendosort_pancreas <- fcrnendosort_pancreas / ve6a_pancreas
    c_fcrnendosort_thymus <- fcrnendosort_thymus / ve6a_thymus
    c_fcrnendosort_spleen <- fcrnendosort_spleen / ve6a_spleen
    c_fcrnendosort_other <- fcrnendosort_other / ve6a_other
    c_fcrnendorecyc_lung <- fcrnendorecyc_lung / ve7b_lung
    c_fcrnendorecyc_liver <- fcrnendorecyc_liver / ve7b_liver
    c_fcrnendorecyc_heart <- fcrnendorecyc_heart / ve7b_heart
    c_fcrnendorecyc_muscle <- fcrnendorecyc_muscle / ve7b_muscle
    c_fcrnendorecyc_skin <- fcrnendorecyc_skin / ve7b_skin
    c_fcrnendorecyc_adipose <- fcrnendorecyc_adipose / ve7b_adipose
    c_fcrnendorecyc_bone <- fcrnendorecyc_bone / ve7b_bone
    c_fcrnendorecyc_brain <- fcrnendorecyc_brain / ve7b_brain
    c_fcrnendorecyc_kidney <- fcrnendorecyc_kidney / ve7b_kidney
    c_fcrnendorecyc_small_intestine <- fcrnendorecyc_small_intestine / ve7b_small_intestine
    c_fcrnendorecyc_large_intestine <- fcrnendorecyc_large_intestine / ve7b_large_intestine
    c_fcrnendorecyc_pancreas <- fcrnendorecyc_pancreas / ve7b_pancreas
    c_fcrnendorecyc_thymus <- fcrnendorecyc_thymus / ve7b_thymus
    c_fcrnendorecyc_spleen <- fcrnendorecyc_spleen / ve7b_spleen
    c_fcrnendorecyc_other <- fcrnendorecyc_other / ve7b_other
    c_fcrnmemint_lung <- fcrnmemint_lung / vism_lung
    c_fcrnmemint_liver <- fcrnmemint_liver / vism_liver
    c_fcrnmemint_heart <- fcrnmemint_heart / vism_heart
    c_fcrnmemint_muscle <- fcrnmemint_muscle / vism_muscle
    c_fcrnmemint_skin <- fcrnmemint_skin / vism_skin
    c_fcrnmemint_adipose <- fcrnmemint_adipose / vism_adipose
    c_fcrnmemint_bone <- fcrnmemint_bone / vism_bone
    c_fcrnmemint_brain <- fcrnmemint_brain / vism_brain
    c_fcrnmemint_kidney <- fcrnmemint_kidney / vism_kidney
    c_fcrnmemint_small_intestine <- fcrnmemint_small_intestine / vism_small_intestine
    c_fcrnmemint_large_intestine <- fcrnmemint_large_intestine / vism_large_intestine
    c_fcrnmemint_pancreas <- fcrnmemint_pancreas / vism_pancreas
    c_fcrnmemint_thymus <- fcrnmemint_thymus / vism_thymus
    c_fcrnmemint_spleen <- fcrnmemint_spleen / vism_spleen
    c_fcrnmemint_other <- fcrnmemint_other / vism_other

    # ---- Dosed antibody: organ vascular spaces ----
    d/dt(vp_lung) <- (plq_lung + lf_lung) * c_central - plq_lung * c_vp_lung - (1 - sigv_lung) * lf_lung * c_vp_lung - clup_lung * fr * c_vp_lung + clup_lung * fr * c_memvas_lung
    d/dt(vp_liver) <- plq_liver * c_vp_lung + (plq_spleen - lf_spleen) * c_vp_spleen + (plq_pancreas - lf_pancreas) * c_vp_pancreas + (plq_small_intestine - lf_small_intestine) * c_vp_small_intestine + (plq_large_intestine - lf_large_intestine) * c_vp_large_intestine - c_vp_liver * (plq_liver - lf_liver + plq_spleen - lf_spleen + plq_pancreas - lf_pancreas + plq_small_intestine - lf_small_intestine + plq_large_intestine - lf_large_intestine) - (1 - sigv_liver) * lf_liver * c_vp_liver - clup_liver * fr * c_vp_liver + clup_liver * fr * c_memvas_liver
    d/dt(vp_heart) <- plq_heart * c_vp_lung - (plq_heart - lf_heart) * c_vp_heart - (1 - sigv_heart) * lf_heart * c_vp_heart - clup_heart * fr * c_vp_heart + clup_heart * fr * c_memvas_heart
    d/dt(vp_muscle) <- plq_muscle * c_vp_lung - (plq_muscle - lf_muscle) * c_vp_muscle - (1 - sigv_muscle) * lf_muscle * c_vp_muscle - clup_muscle * fr * c_vp_muscle + clup_muscle * fr * c_memvas_muscle
    d/dt(vp_skin) <- plq_skin * c_vp_lung - (plq_skin - lf_skin) * c_vp_skin - (1 - sigv_skin) * lf_skin * c_vp_skin - clup_skin * fr * c_vp_skin + clup_skin * fr * c_memvas_skin
    d/dt(vp_adipose) <- plq_adipose * c_vp_lung - (plq_adipose - lf_adipose) * c_vp_adipose - (1 - sigv_adipose) * lf_adipose * c_vp_adipose - clup_adipose * fr * c_vp_adipose + clup_adipose * fr * c_memvas_adipose
    d/dt(vp_bone) <- plq_bone * c_vp_lung - (plq_bone - lf_bone) * c_vp_bone - (1 - sigv_bone) * lf_bone * c_vp_bone - clup_bone * fr * c_vp_bone + clup_bone * fr * c_memvas_bone
    d/dt(vp_brain) <- plq_brain * c_vp_lung - (plq_brain - lf_brain) * c_vp_brain - (1 - sigv_brain) * lf_brain * c_vp_brain - clup_brain * fr * c_vp_brain + clup_brain * fr * c_memvas_brain
    d/dt(vp_kidney) <- plq_kidney * c_vp_lung - (plq_kidney - lf_kidney) * c_vp_kidney - (1 - sigv_kidney) * lf_kidney * c_vp_kidney - clup_kidney * fr * c_vp_kidney + clup_kidney * fr * c_memvas_kidney
    d/dt(vp_small_intestine) <- plq_small_intestine * c_vp_lung - (plq_small_intestine - lf_small_intestine) * c_vp_small_intestine - (1 - sigv_small_intestine) * lf_small_intestine * c_vp_small_intestine - clup_small_intestine * fr * c_vp_small_intestine + clup_small_intestine * fr * c_memvas_small_intestine
    d/dt(vp_large_intestine) <- plq_large_intestine * c_vp_lung - (plq_large_intestine - lf_large_intestine) * c_vp_large_intestine - (1 - sigv_large_intestine) * lf_large_intestine * c_vp_large_intestine - clup_large_intestine * fr * c_vp_large_intestine + clup_large_intestine * fr * c_memvas_large_intestine
    d/dt(vp_pancreas) <- plq_pancreas * c_vp_lung - (plq_pancreas - lf_pancreas) * c_vp_pancreas - (1 - sigv_pancreas) * lf_pancreas * c_vp_pancreas - clup_pancreas * fr * c_vp_pancreas + clup_pancreas * fr * c_memvas_pancreas
    d/dt(vp_thymus) <- plq_thymus * c_vp_lung - (plq_thymus - lf_thymus) * c_vp_thymus - (1 - sigv_thymus) * lf_thymus * c_vp_thymus - clup_thymus * fr * c_vp_thymus + clup_thymus * fr * c_memvas_thymus
    d/dt(vp_spleen) <- plq_spleen * c_vp_lung - (plq_spleen - lf_spleen) * c_vp_spleen - (1 - sigv_spleen) * lf_spleen * c_vp_spleen - clup_spleen * fr * c_vp_spleen + clup_spleen * fr * c_memvas_spleen
    d/dt(vp_other) <- plq_other * c_vp_lung - (plq_other - lf_other) * c_vp_other - (1 - sigv_other) * lf_other * c_vp_other - clup_other * fr * c_vp_other + clup_other * fr * c_memvas_other

    # ---- mab: interstitial spaces ----
    d/dt(is_lung) <- (1 - sigv_lung) * lf_lung * c_vp_lung - clup_lung * (1 - fr) * c_is_lung + clup_lung * (1 - fr) * c_memint_lung - (1 - sigis) * lf_lung * c_is_lung
    d/dt(is_liver) <- (1 - sigv_liver) * lf_liver * c_vp_liver - clup_liver * (1 - fr) * c_is_liver + clup_liver * (1 - fr) * c_memint_liver - (1 - sigis) * lf_liver * c_is_liver
    d/dt(is_heart) <- (1 - sigv_heart) * lf_heart * c_vp_heart - clup_heart * (1 - fr) * c_is_heart + clup_heart * (1 - fr) * c_memint_heart - (1 - sigis) * lf_heart * c_is_heart
    d/dt(is_muscle) <- (1 - sigv_muscle) * lf_muscle * c_vp_muscle - clup_muscle * (1 - fr) * c_is_muscle + clup_muscle * (1 - fr) * c_memint_muscle - (1 - sigis) * lf_muscle * c_is_muscle
    d/dt(is_skin) <- (1 - sigv_skin) * lf_skin * c_vp_skin - clup_skin * (1 - fr) * c_is_skin + clup_skin * (1 - fr) * c_memint_skin - (1 - sigis) * lf_skin * c_is_skin
    d/dt(is_adipose) <- (1 - sigv_adipose) * lf_adipose * c_vp_adipose - clup_adipose * (1 - fr) * c_is_adipose + clup_adipose * (1 - fr) * c_memint_adipose - (1 - sigis) * lf_adipose * c_is_adipose
    d/dt(is_bone) <- (1 - sigv_bone) * lf_bone * c_vp_bone - clup_bone * (1 - fr) * c_is_bone + clup_bone * (1 - fr) * c_memint_bone - (1 - sigis) * lf_bone * c_is_bone
    d/dt(is_brain) <- (1 - sigv_brain) * lf_brain * c_vp_brain - clup_brain * (1 - fr) * c_is_brain + clup_brain * (1 - fr) * c_memint_brain - (1 - sigis) * lf_brain * c_is_brain
    d/dt(is_kidney) <- (1 - sigv_kidney) * lf_kidney * c_vp_kidney - clup_kidney * (1 - fr) * c_is_kidney + clup_kidney * (1 - fr) * c_memint_kidney - (1 - sigis) * lf_kidney * c_is_kidney
    d/dt(is_small_intestine) <- (1 - sigv_small_intestine) * lf_small_intestine * c_vp_small_intestine - clup_small_intestine * (1 - fr) * c_is_small_intestine + clup_small_intestine * (1 - fr) * c_memint_small_intestine - (1 - sigis) * lf_small_intestine * c_is_small_intestine
    d/dt(is_large_intestine) <- (1 - sigv_large_intestine) * lf_large_intestine * c_vp_large_intestine - clup_large_intestine * (1 - fr) * c_is_large_intestine + clup_large_intestine * (1 - fr) * c_memint_large_intestine - (1 - sigis) * lf_large_intestine * c_is_large_intestine
    d/dt(is_pancreas) <- (1 - sigv_pancreas) * lf_pancreas * c_vp_pancreas - clup_pancreas * (1 - fr) * c_is_pancreas + clup_pancreas * (1 - fr) * c_memint_pancreas - (1 - sigis) * lf_pancreas * c_is_pancreas
    d/dt(is_thymus) <- (1 - sigv_thymus) * lf_thymus * c_vp_thymus - clup_thymus * (1 - fr) * c_is_thymus + clup_thymus * (1 - fr) * c_memint_thymus - (1 - sigis) * lf_thymus * c_is_thymus
    d/dt(is_spleen) <- (1 - sigv_spleen) * lf_spleen * c_vp_spleen - clup_spleen * (1 - fr) * c_is_spleen + clup_spleen * (1 - fr) * c_memint_spleen - (1 - sigis) * lf_spleen * c_is_spleen
    d/dt(is_other) <- (1 - sigv_other) * lf_other * c_vp_other - clup_other * (1 - fr) * c_is_other + clup_other * (1 - fr) * c_memint_other - (1 - sigis) * lf_other * c_is_other

    # ---- mab: vascular-side membrane ----
    d/dt(memvas_lung) <- clup_lung * fr * c_vp_lung - clup_lung * fr * c_memvas_lung - clup_lung * fr * c_memvas_lung + clup_lung * fr * c_endorecyc_lung - kon7 * c_memvas_lung * c_fcrnmemvas_lung * vvm_lung + koff7 * c_memvasfr1_lung * vvm_lung - konps * c_memvas_lung * cmem * vvm_lung + koffps * c_memvasns_lung * vvm_lung
    d/dt(memvas_liver) <- clup_liver * fr * c_vp_liver - clup_liver * fr * c_memvas_liver - clup_liver * fr * c_memvas_liver + clup_liver * fr * c_endorecyc_liver - kon7 * c_memvas_liver * c_fcrnmemvas_liver * vvm_liver + koff7 * c_memvasfr1_liver * vvm_liver - konps * c_memvas_liver * cmem * vvm_liver + koffps * c_memvasns_liver * vvm_liver
    d/dt(memvas_heart) <- clup_heart * fr * c_vp_heart - clup_heart * fr * c_memvas_heart - clup_heart * fr * c_memvas_heart + clup_heart * fr * c_endorecyc_heart - kon7 * c_memvas_heart * c_fcrnmemvas_heart * vvm_heart + koff7 * c_memvasfr1_heart * vvm_heart - konps * c_memvas_heart * cmem * vvm_heart + koffps * c_memvasns_heart * vvm_heart
    d/dt(memvas_muscle) <- clup_muscle * fr * c_vp_muscle - clup_muscle * fr * c_memvas_muscle - clup_muscle * fr * c_memvas_muscle + clup_muscle * fr * c_endorecyc_muscle - kon7 * c_memvas_muscle * c_fcrnmemvas_muscle * vvm_muscle + koff7 * c_memvasfr1_muscle * vvm_muscle - konps * c_memvas_muscle * cmem * vvm_muscle + koffps * c_memvasns_muscle * vvm_muscle
    d/dt(memvas_skin) <- clup_skin * fr * c_vp_skin - clup_skin * fr * c_memvas_skin - clup_skin * fr * c_memvas_skin + clup_skin * fr * c_endorecyc_skin - kon7 * c_memvas_skin * c_fcrnmemvas_skin * vvm_skin + koff7 * c_memvasfr1_skin * vvm_skin - konps * c_memvas_skin * cmem * vvm_skin + koffps * c_memvasns_skin * vvm_skin
    d/dt(memvas_adipose) <- clup_adipose * fr * c_vp_adipose - clup_adipose * fr * c_memvas_adipose - clup_adipose * fr * c_memvas_adipose + clup_adipose * fr * c_endorecyc_adipose - kon7 * c_memvas_adipose * c_fcrnmemvas_adipose * vvm_adipose + koff7 * c_memvasfr1_adipose * vvm_adipose - konps * c_memvas_adipose * cmem * vvm_adipose + koffps * c_memvasns_adipose * vvm_adipose
    d/dt(memvas_bone) <- clup_bone * fr * c_vp_bone - clup_bone * fr * c_memvas_bone - clup_bone * fr * c_memvas_bone + clup_bone * fr * c_endorecyc_bone - kon7 * c_memvas_bone * c_fcrnmemvas_bone * vvm_bone + koff7 * c_memvasfr1_bone * vvm_bone - konps * c_memvas_bone * cmem * vvm_bone + koffps * c_memvasns_bone * vvm_bone
    d/dt(memvas_brain) <- clup_brain * fr * c_vp_brain - clup_brain * fr * c_memvas_brain - clup_brain * fr * c_memvas_brain + clup_brain * fr * c_endorecyc_brain - kon7 * c_memvas_brain * c_fcrnmemvas_brain * vvm_brain + koff7 * c_memvasfr1_brain * vvm_brain - konps * c_memvas_brain * cmem * vvm_brain + koffps * c_memvasns_brain * vvm_brain
    d/dt(memvas_kidney) <- clup_kidney * fr * c_vp_kidney - clup_kidney * fr * c_memvas_kidney - clup_kidney * fr * c_memvas_kidney + clup_kidney * fr * c_endorecyc_kidney - kon7 * c_memvas_kidney * c_fcrnmemvas_kidney * vvm_kidney + koff7 * c_memvasfr1_kidney * vvm_kidney - konps * c_memvas_kidney * cmem * vvm_kidney + koffps * c_memvasns_kidney * vvm_kidney
    d/dt(memvas_small_intestine) <- clup_small_intestine * fr * c_vp_small_intestine - clup_small_intestine * fr * c_memvas_small_intestine - clup_small_intestine * fr * c_memvas_small_intestine + clup_small_intestine * fr * c_endorecyc_small_intestine - kon7 * c_memvas_small_intestine * c_fcrnmemvas_small_intestine * vvm_small_intestine + koff7 * c_memvasfr1_small_intestine * vvm_small_intestine - konps * c_memvas_small_intestine * cmem * vvm_small_intestine + koffps * c_memvasns_small_intestine * vvm_small_intestine
    d/dt(memvas_large_intestine) <- clup_large_intestine * fr * c_vp_large_intestine - clup_large_intestine * fr * c_memvas_large_intestine - clup_large_intestine * fr * c_memvas_large_intestine + clup_large_intestine * fr * c_endorecyc_large_intestine - kon7 * c_memvas_large_intestine * c_fcrnmemvas_large_intestine * vvm_large_intestine + koff7 * c_memvasfr1_large_intestine * vvm_large_intestine - konps * c_memvas_large_intestine * cmem * vvm_large_intestine + koffps * c_memvasns_large_intestine * vvm_large_intestine
    d/dt(memvas_pancreas) <- clup_pancreas * fr * c_vp_pancreas - clup_pancreas * fr * c_memvas_pancreas - clup_pancreas * fr * c_memvas_pancreas + clup_pancreas * fr * c_endorecyc_pancreas - kon7 * c_memvas_pancreas * c_fcrnmemvas_pancreas * vvm_pancreas + koff7 * c_memvasfr1_pancreas * vvm_pancreas - konps * c_memvas_pancreas * cmem * vvm_pancreas + koffps * c_memvasns_pancreas * vvm_pancreas
    d/dt(memvas_thymus) <- clup_thymus * fr * c_vp_thymus - clup_thymus * fr * c_memvas_thymus - clup_thymus * fr * c_memvas_thymus + clup_thymus * fr * c_endorecyc_thymus - kon7 * c_memvas_thymus * c_fcrnmemvas_thymus * vvm_thymus + koff7 * c_memvasfr1_thymus * vvm_thymus - konps * c_memvas_thymus * cmem * vvm_thymus + koffps * c_memvasns_thymus * vvm_thymus
    d/dt(memvas_spleen) <- clup_spleen * fr * c_vp_spleen - clup_spleen * fr * c_memvas_spleen - clup_spleen * fr * c_memvas_spleen + clup_spleen * fr * c_endorecyc_spleen - kon7 * c_memvas_spleen * c_fcrnmemvas_spleen * vvm_spleen + koff7 * c_memvasfr1_spleen * vvm_spleen - konps * c_memvas_spleen * cmem * vvm_spleen + koffps * c_memvasns_spleen * vvm_spleen
    d/dt(memvas_other) <- clup_other * fr * c_vp_other - clup_other * fr * c_memvas_other - clup_other * fr * c_memvas_other + clup_other * fr * c_endorecyc_other - kon7 * c_memvas_other * c_fcrnmemvas_other * vvm_other + koff7 * c_memvasfr1_other * vvm_other - konps * c_memvas_other * cmem * vvm_other + koffps * c_memvasns_other * vvm_other

    # ---- mab: early endosome (pH 7.4) ----
    d/dt(endoearly_lung) <- clup_lung * fr * c_memvas_lung + clup_lung * (1 - fr) * c_memint_lung - clup_lung * c_endoearly_lung - kon7 * c_endoearly_lung * c_fcrnendoearly_lung * ve7_lung + koff7 * c_endoearlyfr1_lung * ve7_lung + kintps * c_memvasns_lung * vvm_lung + kintps * c_memintns_lung * vism_lung
    d/dt(endoearly_liver) <- clup_liver * fr * c_memvas_liver + clup_liver * (1 - fr) * c_memint_liver - clup_liver * c_endoearly_liver - kon7 * c_endoearly_liver * c_fcrnendoearly_liver * ve7_liver + koff7 * c_endoearlyfr1_liver * ve7_liver + kintps * c_memvasns_liver * vvm_liver + kintps * c_memintns_liver * vism_liver
    d/dt(endoearly_heart) <- clup_heart * fr * c_memvas_heart + clup_heart * (1 - fr) * c_memint_heart - clup_heart * c_endoearly_heart - kon7 * c_endoearly_heart * c_fcrnendoearly_heart * ve7_heart + koff7 * c_endoearlyfr1_heart * ve7_heart + kintps * c_memvasns_heart * vvm_heart + kintps * c_memintns_heart * vism_heart
    d/dt(endoearly_muscle) <- clup_muscle * fr * c_memvas_muscle + clup_muscle * (1 - fr) * c_memint_muscle - clup_muscle * c_endoearly_muscle - kon7 * c_endoearly_muscle * c_fcrnendoearly_muscle * ve7_muscle + koff7 * c_endoearlyfr1_muscle * ve7_muscle + kintps * c_memvasns_muscle * vvm_muscle + kintps * c_memintns_muscle * vism_muscle
    d/dt(endoearly_skin) <- clup_skin * fr * c_memvas_skin + clup_skin * (1 - fr) * c_memint_skin - clup_skin * c_endoearly_skin - kon7 * c_endoearly_skin * c_fcrnendoearly_skin * ve7_skin + koff7 * c_endoearlyfr1_skin * ve7_skin + kintps * c_memvasns_skin * vvm_skin + kintps * c_memintns_skin * vism_skin
    d/dt(endoearly_adipose) <- clup_adipose * fr * c_memvas_adipose + clup_adipose * (1 - fr) * c_memint_adipose - clup_adipose * c_endoearly_adipose - kon7 * c_endoearly_adipose * c_fcrnendoearly_adipose * ve7_adipose + koff7 * c_endoearlyfr1_adipose * ve7_adipose + kintps * c_memvasns_adipose * vvm_adipose + kintps * c_memintns_adipose * vism_adipose
    d/dt(endoearly_bone) <- clup_bone * fr * c_memvas_bone + clup_bone * (1 - fr) * c_memint_bone - clup_bone * c_endoearly_bone - kon7 * c_endoearly_bone * c_fcrnendoearly_bone * ve7_bone + koff7 * c_endoearlyfr1_bone * ve7_bone + kintps * c_memvasns_bone * vvm_bone + kintps * c_memintns_bone * vism_bone
    d/dt(endoearly_brain) <- clup_brain * fr * c_memvas_brain + clup_brain * (1 - fr) * c_memint_brain - clup_brain * c_endoearly_brain - kon7 * c_endoearly_brain * c_fcrnendoearly_brain * ve7_brain + koff7 * c_endoearlyfr1_brain * ve7_brain + kintps * c_memvasns_brain * vvm_brain + kintps * c_memintns_brain * vism_brain
    d/dt(endoearly_kidney) <- clup_kidney * fr * c_memvas_kidney + clup_kidney * (1 - fr) * c_memint_kidney - clup_kidney * c_endoearly_kidney - kon7 * c_endoearly_kidney * c_fcrnendoearly_kidney * ve7_kidney + koff7 * c_endoearlyfr1_kidney * ve7_kidney + kintps * c_memvasns_kidney * vvm_kidney + kintps * c_memintns_kidney * vism_kidney
    d/dt(endoearly_small_intestine) <- clup_small_intestine * fr * c_memvas_small_intestine + clup_small_intestine * (1 - fr) * c_memint_small_intestine - clup_small_intestine * c_endoearly_small_intestine - kon7 * c_endoearly_small_intestine * c_fcrnendoearly_small_intestine * ve7_small_intestine + koff7 * c_endoearlyfr1_small_intestine * ve7_small_intestine + kintps * c_memvasns_small_intestine * vvm_small_intestine + kintps * c_memintns_small_intestine * vism_small_intestine
    d/dt(endoearly_large_intestine) <- clup_large_intestine * fr * c_memvas_large_intestine + clup_large_intestine * (1 - fr) * c_memint_large_intestine - clup_large_intestine * c_endoearly_large_intestine - kon7 * c_endoearly_large_intestine * c_fcrnendoearly_large_intestine * ve7_large_intestine + koff7 * c_endoearlyfr1_large_intestine * ve7_large_intestine + kintps * c_memvasns_large_intestine * vvm_large_intestine + kintps * c_memintns_large_intestine * vism_large_intestine
    d/dt(endoearly_pancreas) <- clup_pancreas * fr * c_memvas_pancreas + clup_pancreas * (1 - fr) * c_memint_pancreas - clup_pancreas * c_endoearly_pancreas - kon7 * c_endoearly_pancreas * c_fcrnendoearly_pancreas * ve7_pancreas + koff7 * c_endoearlyfr1_pancreas * ve7_pancreas + kintps * c_memvasns_pancreas * vvm_pancreas + kintps * c_memintns_pancreas * vism_pancreas
    d/dt(endoearly_thymus) <- clup_thymus * fr * c_memvas_thymus + clup_thymus * (1 - fr) * c_memint_thymus - clup_thymus * c_endoearly_thymus - kon7 * c_endoearly_thymus * c_fcrnendoearly_thymus * ve7_thymus + koff7 * c_endoearlyfr1_thymus * ve7_thymus + kintps * c_memvasns_thymus * vvm_thymus + kintps * c_memintns_thymus * vism_thymus
    d/dt(endoearly_spleen) <- clup_spleen * fr * c_memvas_spleen + clup_spleen * (1 - fr) * c_memint_spleen - clup_spleen * c_endoearly_spleen - kon7 * c_endoearly_spleen * c_fcrnendoearly_spleen * ve7_spleen + koff7 * c_endoearlyfr1_spleen * ve7_spleen + kintps * c_memvasns_spleen * vvm_spleen + kintps * c_memintns_spleen * vism_spleen
    d/dt(endoearly_other) <- clup_other * fr * c_memvas_other + clup_other * (1 - fr) * c_memint_other - clup_other * c_endoearly_other - kon7 * c_endoearly_other * c_fcrnendoearly_other * ve7_other + koff7 * c_endoearlyfr1_other * ve7_other + kintps * c_memvasns_other * vvm_other + kintps * c_memintns_other * vism_other

    # ---- mab: sorting endosome (pH 6.0) ----
    d/dt(endosort_lung) <- clup_lung * c_endoearly_lung - clup_lung * c_endosort_lung - kon6 * c_endosort_lung * c_fcrnendosort_lung * ve6a_lung + koff6 * c_endosortfr1_lung * ve6a_lung
    d/dt(endosort_liver) <- clup_liver * c_endoearly_liver - clup_liver * c_endosort_liver - kon6 * c_endosort_liver * c_fcrnendosort_liver * ve6a_liver + koff6 * c_endosortfr1_liver * ve6a_liver
    d/dt(endosort_heart) <- clup_heart * c_endoearly_heart - clup_heart * c_endosort_heart - kon6 * c_endosort_heart * c_fcrnendosort_heart * ve6a_heart + koff6 * c_endosortfr1_heart * ve6a_heart
    d/dt(endosort_muscle) <- clup_muscle * c_endoearly_muscle - clup_muscle * c_endosort_muscle - kon6 * c_endosort_muscle * c_fcrnendosort_muscle * ve6a_muscle + koff6 * c_endosortfr1_muscle * ve6a_muscle
    d/dt(endosort_skin) <- clup_skin * c_endoearly_skin - clup_skin * c_endosort_skin - kon6 * c_endosort_skin * c_fcrnendosort_skin * ve6a_skin + koff6 * c_endosortfr1_skin * ve6a_skin
    d/dt(endosort_adipose) <- clup_adipose * c_endoearly_adipose - clup_adipose * c_endosort_adipose - kon6 * c_endosort_adipose * c_fcrnendosort_adipose * ve6a_adipose + koff6 * c_endosortfr1_adipose * ve6a_adipose
    d/dt(endosort_bone) <- clup_bone * c_endoearly_bone - clup_bone * c_endosort_bone - kon6 * c_endosort_bone * c_fcrnendosort_bone * ve6a_bone + koff6 * c_endosortfr1_bone * ve6a_bone
    d/dt(endosort_brain) <- clup_brain * c_endoearly_brain - clup_brain * c_endosort_brain - kon6 * c_endosort_brain * c_fcrnendosort_brain * ve6a_brain + koff6 * c_endosortfr1_brain * ve6a_brain
    d/dt(endosort_kidney) <- clup_kidney * c_endoearly_kidney - clup_kidney * c_endosort_kidney - kon6 * c_endosort_kidney * c_fcrnendosort_kidney * ve6a_kidney + koff6 * c_endosortfr1_kidney * ve6a_kidney
    d/dt(endosort_small_intestine) <- clup_small_intestine * c_endoearly_small_intestine - clup_small_intestine * c_endosort_small_intestine - kon6 * c_endosort_small_intestine * c_fcrnendosort_small_intestine * ve6a_small_intestine + koff6 * c_endosortfr1_small_intestine * ve6a_small_intestine
    d/dt(endosort_large_intestine) <- clup_large_intestine * c_endoearly_large_intestine - clup_large_intestine * c_endosort_large_intestine - kon6 * c_endosort_large_intestine * c_fcrnendosort_large_intestine * ve6a_large_intestine + koff6 * c_endosortfr1_large_intestine * ve6a_large_intestine
    d/dt(endosort_pancreas) <- clup_pancreas * c_endoearly_pancreas - clup_pancreas * c_endosort_pancreas - kon6 * c_endosort_pancreas * c_fcrnendosort_pancreas * ve6a_pancreas + koff6 * c_endosortfr1_pancreas * ve6a_pancreas
    d/dt(endosort_thymus) <- clup_thymus * c_endoearly_thymus - clup_thymus * c_endosort_thymus - kon6 * c_endosort_thymus * c_fcrnendosort_thymus * ve6a_thymus + koff6 * c_endosortfr1_thymus * ve6a_thymus
    d/dt(endosort_spleen) <- clup_spleen * c_endoearly_spleen - clup_spleen * c_endosort_spleen - kon6 * c_endosort_spleen * c_fcrnendosort_spleen * ve6a_spleen + koff6 * c_endosortfr1_spleen * ve6a_spleen
    d/dt(endosort_other) <- clup_other * c_endoearly_other - clup_other * c_endosort_other - kon6 * c_endosort_other * c_fcrnendosort_other * ve6a_other + koff6 * c_endosortfr1_other * ve6a_other

    # ---- mab: recycling endosome (pH 7.4) ----
    # Only the (1 - Prob_deg) fraction of unbound mAb escapes lysosomal routing.
    d/dt(endorecyc_lung) <- clup_lung * (1 - probdeg) * c_endosort_lung - clup_lung * c_endorecyc_lung - kon7 * c_endorecyc_lung * c_fcrnendorecyc_lung * ve7b_lung + koff7 * c_endorecycfr1_lung * ve7b_lung
    d/dt(endorecyc_liver) <- clup_liver * (1 - probdeg) * c_endosort_liver - clup_liver * c_endorecyc_liver - kon7 * c_endorecyc_liver * c_fcrnendorecyc_liver * ve7b_liver + koff7 * c_endorecycfr1_liver * ve7b_liver
    d/dt(endorecyc_heart) <- clup_heart * (1 - probdeg) * c_endosort_heart - clup_heart * c_endorecyc_heart - kon7 * c_endorecyc_heart * c_fcrnendorecyc_heart * ve7b_heart + koff7 * c_endorecycfr1_heart * ve7b_heart
    d/dt(endorecyc_muscle) <- clup_muscle * (1 - probdeg) * c_endosort_muscle - clup_muscle * c_endorecyc_muscle - kon7 * c_endorecyc_muscle * c_fcrnendorecyc_muscle * ve7b_muscle + koff7 * c_endorecycfr1_muscle * ve7b_muscle
    d/dt(endorecyc_skin) <- clup_skin * (1 - probdeg) * c_endosort_skin - clup_skin * c_endorecyc_skin - kon7 * c_endorecyc_skin * c_fcrnendorecyc_skin * ve7b_skin + koff7 * c_endorecycfr1_skin * ve7b_skin
    d/dt(endorecyc_adipose) <- clup_adipose * (1 - probdeg) * c_endosort_adipose - clup_adipose * c_endorecyc_adipose - kon7 * c_endorecyc_adipose * c_fcrnendorecyc_adipose * ve7b_adipose + koff7 * c_endorecycfr1_adipose * ve7b_adipose
    d/dt(endorecyc_bone) <- clup_bone * (1 - probdeg) * c_endosort_bone - clup_bone * c_endorecyc_bone - kon7 * c_endorecyc_bone * c_fcrnendorecyc_bone * ve7b_bone + koff7 * c_endorecycfr1_bone * ve7b_bone
    d/dt(endorecyc_brain) <- clup_brain * (1 - probdeg) * c_endosort_brain - clup_brain * c_endorecyc_brain - kon7 * c_endorecyc_brain * c_fcrnendorecyc_brain * ve7b_brain + koff7 * c_endorecycfr1_brain * ve7b_brain
    d/dt(endorecyc_kidney) <- clup_kidney * (1 - probdeg) * c_endosort_kidney - clup_kidney * c_endorecyc_kidney - kon7 * c_endorecyc_kidney * c_fcrnendorecyc_kidney * ve7b_kidney + koff7 * c_endorecycfr1_kidney * ve7b_kidney
    d/dt(endorecyc_small_intestine) <- clup_small_intestine * (1 - probdeg) * c_endosort_small_intestine - clup_small_intestine * c_endorecyc_small_intestine - kon7 * c_endorecyc_small_intestine * c_fcrnendorecyc_small_intestine * ve7b_small_intestine + koff7 * c_endorecycfr1_small_intestine * ve7b_small_intestine
    d/dt(endorecyc_large_intestine) <- clup_large_intestine * (1 - probdeg) * c_endosort_large_intestine - clup_large_intestine * c_endorecyc_large_intestine - kon7 * c_endorecyc_large_intestine * c_fcrnendorecyc_large_intestine * ve7b_large_intestine + koff7 * c_endorecycfr1_large_intestine * ve7b_large_intestine
    d/dt(endorecyc_pancreas) <- clup_pancreas * (1 - probdeg) * c_endosort_pancreas - clup_pancreas * c_endorecyc_pancreas - kon7 * c_endorecyc_pancreas * c_fcrnendorecyc_pancreas * ve7b_pancreas + koff7 * c_endorecycfr1_pancreas * ve7b_pancreas
    d/dt(endorecyc_thymus) <- clup_thymus * (1 - probdeg) * c_endosort_thymus - clup_thymus * c_endorecyc_thymus - kon7 * c_endorecyc_thymus * c_fcrnendorecyc_thymus * ve7b_thymus + koff7 * c_endorecycfr1_thymus * ve7b_thymus
    d/dt(endorecyc_spleen) <- clup_spleen * (1 - probdeg) * c_endosort_spleen - clup_spleen * c_endorecyc_spleen - kon7 * c_endorecyc_spleen * c_fcrnendorecyc_spleen * ve7b_spleen + koff7 * c_endorecycfr1_spleen * ve7b_spleen
    d/dt(endorecyc_other) <- clup_other * (1 - probdeg) * c_endosort_other - clup_other * c_endorecyc_other - kon7 * c_endorecyc_other * c_fcrnendorecyc_other * ve7b_other + koff7 * c_endorecycfr1_other * ve7b_other

    # ---- mab: interstitial-side membrane ----
    d/dt(memint_lung) <- clup_lung * (1 - fr) * c_is_lung - clup_lung * (1 - fr) * c_memint_lung + clup_lung * (1 - fr) * c_endorecyc_lung - clup_lung * (1 - fr) * c_memint_lung - kon7 * c_memint_lung * c_fcrnmemint_lung * vism_lung + koff7 * c_memintfr1_lung * vism_lung - konps * c_memint_lung * cmem * vism_lung + koffps * c_memintns_lung * vism_lung
    d/dt(memint_liver) <- clup_liver * (1 - fr) * c_is_liver - clup_liver * (1 - fr) * c_memint_liver + clup_liver * (1 - fr) * c_endorecyc_liver - clup_liver * (1 - fr) * c_memint_liver - kon7 * c_memint_liver * c_fcrnmemint_liver * vism_liver + koff7 * c_memintfr1_liver * vism_liver - konps * c_memint_liver * cmem * vism_liver + koffps * c_memintns_liver * vism_liver
    d/dt(memint_heart) <- clup_heart * (1 - fr) * c_is_heart - clup_heart * (1 - fr) * c_memint_heart + clup_heart * (1 - fr) * c_endorecyc_heart - clup_heart * (1 - fr) * c_memint_heart - kon7 * c_memint_heart * c_fcrnmemint_heart * vism_heart + koff7 * c_memintfr1_heart * vism_heart - konps * c_memint_heart * cmem * vism_heart + koffps * c_memintns_heart * vism_heart
    d/dt(memint_muscle) <- clup_muscle * (1 - fr) * c_is_muscle - clup_muscle * (1 - fr) * c_memint_muscle + clup_muscle * (1 - fr) * c_endorecyc_muscle - clup_muscle * (1 - fr) * c_memint_muscle - kon7 * c_memint_muscle * c_fcrnmemint_muscle * vism_muscle + koff7 * c_memintfr1_muscle * vism_muscle - konps * c_memint_muscle * cmem * vism_muscle + koffps * c_memintns_muscle * vism_muscle
    d/dt(memint_skin) <- clup_skin * (1 - fr) * c_is_skin - clup_skin * (1 - fr) * c_memint_skin + clup_skin * (1 - fr) * c_endorecyc_skin - clup_skin * (1 - fr) * c_memint_skin - kon7 * c_memint_skin * c_fcrnmemint_skin * vism_skin + koff7 * c_memintfr1_skin * vism_skin - konps * c_memint_skin * cmem * vism_skin + koffps * c_memintns_skin * vism_skin
    d/dt(memint_adipose) <- clup_adipose * (1 - fr) * c_is_adipose - clup_adipose * (1 - fr) * c_memint_adipose + clup_adipose * (1 - fr) * c_endorecyc_adipose - clup_adipose * (1 - fr) * c_memint_adipose - kon7 * c_memint_adipose * c_fcrnmemint_adipose * vism_adipose + koff7 * c_memintfr1_adipose * vism_adipose - konps * c_memint_adipose * cmem * vism_adipose + koffps * c_memintns_adipose * vism_adipose
    d/dt(memint_bone) <- clup_bone * (1 - fr) * c_is_bone - clup_bone * (1 - fr) * c_memint_bone + clup_bone * (1 - fr) * c_endorecyc_bone - clup_bone * (1 - fr) * c_memint_bone - kon7 * c_memint_bone * c_fcrnmemint_bone * vism_bone + koff7 * c_memintfr1_bone * vism_bone - konps * c_memint_bone * cmem * vism_bone + koffps * c_memintns_bone * vism_bone
    d/dt(memint_brain) <- clup_brain * (1 - fr) * c_is_brain - clup_brain * (1 - fr) * c_memint_brain + clup_brain * (1 - fr) * c_endorecyc_brain - clup_brain * (1 - fr) * c_memint_brain - kon7 * c_memint_brain * c_fcrnmemint_brain * vism_brain + koff7 * c_memintfr1_brain * vism_brain - konps * c_memint_brain * cmem * vism_brain + koffps * c_memintns_brain * vism_brain
    d/dt(memint_kidney) <- clup_kidney * (1 - fr) * c_is_kidney - clup_kidney * (1 - fr) * c_memint_kidney + clup_kidney * (1 - fr) * c_endorecyc_kidney - clup_kidney * (1 - fr) * c_memint_kidney - kon7 * c_memint_kidney * c_fcrnmemint_kidney * vism_kidney + koff7 * c_memintfr1_kidney * vism_kidney - konps * c_memint_kidney * cmem * vism_kidney + koffps * c_memintns_kidney * vism_kidney
    d/dt(memint_small_intestine) <- clup_small_intestine * (1 - fr) * c_is_small_intestine - clup_small_intestine * (1 - fr) * c_memint_small_intestine + clup_small_intestine * (1 - fr) * c_endorecyc_small_intestine - clup_small_intestine * (1 - fr) * c_memint_small_intestine - kon7 * c_memint_small_intestine * c_fcrnmemint_small_intestine * vism_small_intestine + koff7 * c_memintfr1_small_intestine * vism_small_intestine - konps * c_memint_small_intestine * cmem * vism_small_intestine + koffps * c_memintns_small_intestine * vism_small_intestine
    d/dt(memint_large_intestine) <- clup_large_intestine * (1 - fr) * c_is_large_intestine - clup_large_intestine * (1 - fr) * c_memint_large_intestine + clup_large_intestine * (1 - fr) * c_endorecyc_large_intestine - clup_large_intestine * (1 - fr) * c_memint_large_intestine - kon7 * c_memint_large_intestine * c_fcrnmemint_large_intestine * vism_large_intestine + koff7 * c_memintfr1_large_intestine * vism_large_intestine - konps * c_memint_large_intestine * cmem * vism_large_intestine + koffps * c_memintns_large_intestine * vism_large_intestine
    d/dt(memint_pancreas) <- clup_pancreas * (1 - fr) * c_is_pancreas - clup_pancreas * (1 - fr) * c_memint_pancreas + clup_pancreas * (1 - fr) * c_endorecyc_pancreas - clup_pancreas * (1 - fr) * c_memint_pancreas - kon7 * c_memint_pancreas * c_fcrnmemint_pancreas * vism_pancreas + koff7 * c_memintfr1_pancreas * vism_pancreas - konps * c_memint_pancreas * cmem * vism_pancreas + koffps * c_memintns_pancreas * vism_pancreas
    d/dt(memint_thymus) <- clup_thymus * (1 - fr) * c_is_thymus - clup_thymus * (1 - fr) * c_memint_thymus + clup_thymus * (1 - fr) * c_endorecyc_thymus - clup_thymus * (1 - fr) * c_memint_thymus - kon7 * c_memint_thymus * c_fcrnmemint_thymus * vism_thymus + koff7 * c_memintfr1_thymus * vism_thymus - konps * c_memint_thymus * cmem * vism_thymus + koffps * c_memintns_thymus * vism_thymus
    d/dt(memint_spleen) <- clup_spleen * (1 - fr) * c_is_spleen - clup_spleen * (1 - fr) * c_memint_spleen + clup_spleen * (1 - fr) * c_endorecyc_spleen - clup_spleen * (1 - fr) * c_memint_spleen - kon7 * c_memint_spleen * c_fcrnmemint_spleen * vism_spleen + koff7 * c_memintfr1_spleen * vism_spleen - konps * c_memint_spleen * cmem * vism_spleen + koffps * c_memintns_spleen * vism_spleen
    d/dt(memint_other) <- clup_other * (1 - fr) * c_is_other - clup_other * (1 - fr) * c_memint_other + clup_other * (1 - fr) * c_endorecyc_other - clup_other * (1 - fr) * c_memint_other - kon7 * c_memint_other * c_fcrnmemint_other * vism_other + koff7 * c_memintfr1_other * vism_other - konps * c_memint_other * cmem * vism_other + koffps * c_memintns_other * vism_other

    # ---- mab: 1:1 FcRn complexes ----
    d/dt(memvasfr1_lung) <- kon7 * c_memvas_lung * c_fcrnmemvas_lung * vvm_lung - koff7 * c_memvasfr1_lung * vvm_lung - kon7b * c_memvasfr1_lung * c_fcrnmemvas_lung * vvm_lung + koff7b * c_memvasfr2_lung * vvm_lung + clup_lung * fr * c_endorecycfr1_lung - clup_lung * fr * c_memvasfr1_lung - kdegab * c_memvasfr1_lung * vvm_lung
    d/dt(memvasfr1_liver) <- kon7 * c_memvas_liver * c_fcrnmemvas_liver * vvm_liver - koff7 * c_memvasfr1_liver * vvm_liver - kon7b * c_memvasfr1_liver * c_fcrnmemvas_liver * vvm_liver + koff7b * c_memvasfr2_liver * vvm_liver + clup_liver * fr * c_endorecycfr1_liver - clup_liver * fr * c_memvasfr1_liver - kdegab * c_memvasfr1_liver * vvm_liver
    d/dt(memvasfr1_heart) <- kon7 * c_memvas_heart * c_fcrnmemvas_heart * vvm_heart - koff7 * c_memvasfr1_heart * vvm_heart - kon7b * c_memvasfr1_heart * c_fcrnmemvas_heart * vvm_heart + koff7b * c_memvasfr2_heart * vvm_heart + clup_heart * fr * c_endorecycfr1_heart - clup_heart * fr * c_memvasfr1_heart - kdegab * c_memvasfr1_heart * vvm_heart
    d/dt(memvasfr1_muscle) <- kon7 * c_memvas_muscle * c_fcrnmemvas_muscle * vvm_muscle - koff7 * c_memvasfr1_muscle * vvm_muscle - kon7b * c_memvasfr1_muscle * c_fcrnmemvas_muscle * vvm_muscle + koff7b * c_memvasfr2_muscle * vvm_muscle + clup_muscle * fr * c_endorecycfr1_muscle - clup_muscle * fr * c_memvasfr1_muscle - kdegab * c_memvasfr1_muscle * vvm_muscle
    d/dt(memvasfr1_skin) <- kon7 * c_memvas_skin * c_fcrnmemvas_skin * vvm_skin - koff7 * c_memvasfr1_skin * vvm_skin - kon7b * c_memvasfr1_skin * c_fcrnmemvas_skin * vvm_skin + koff7b * c_memvasfr2_skin * vvm_skin + clup_skin * fr * c_endorecycfr1_skin - clup_skin * fr * c_memvasfr1_skin - kdegab * c_memvasfr1_skin * vvm_skin
    d/dt(memvasfr1_adipose) <- kon7 * c_memvas_adipose * c_fcrnmemvas_adipose * vvm_adipose - koff7 * c_memvasfr1_adipose * vvm_adipose - kon7b * c_memvasfr1_adipose * c_fcrnmemvas_adipose * vvm_adipose + koff7b * c_memvasfr2_adipose * vvm_adipose + clup_adipose * fr * c_endorecycfr1_adipose - clup_adipose * fr * c_memvasfr1_adipose - kdegab * c_memvasfr1_adipose * vvm_adipose
    d/dt(memvasfr1_bone) <- kon7 * c_memvas_bone * c_fcrnmemvas_bone * vvm_bone - koff7 * c_memvasfr1_bone * vvm_bone - kon7b * c_memvasfr1_bone * c_fcrnmemvas_bone * vvm_bone + koff7b * c_memvasfr2_bone * vvm_bone + clup_bone * fr * c_endorecycfr1_bone - clup_bone * fr * c_memvasfr1_bone - kdegab * c_memvasfr1_bone * vvm_bone
    d/dt(memvasfr1_brain) <- kon7 * c_memvas_brain * c_fcrnmemvas_brain * vvm_brain - koff7 * c_memvasfr1_brain * vvm_brain - kon7b * c_memvasfr1_brain * c_fcrnmemvas_brain * vvm_brain + koff7b * c_memvasfr2_brain * vvm_brain + clup_brain * fr * c_endorecycfr1_brain - clup_brain * fr * c_memvasfr1_brain - kdegab * c_memvasfr1_brain * vvm_brain
    d/dt(memvasfr1_kidney) <- kon7 * c_memvas_kidney * c_fcrnmemvas_kidney * vvm_kidney - koff7 * c_memvasfr1_kidney * vvm_kidney - kon7b * c_memvasfr1_kidney * c_fcrnmemvas_kidney * vvm_kidney + koff7b * c_memvasfr2_kidney * vvm_kidney + clup_kidney * fr * c_endorecycfr1_kidney - clup_kidney * fr * c_memvasfr1_kidney - kdegab * c_memvasfr1_kidney * vvm_kidney
    d/dt(memvasfr1_small_intestine) <- kon7 * c_memvas_small_intestine * c_fcrnmemvas_small_intestine * vvm_small_intestine - koff7 * c_memvasfr1_small_intestine * vvm_small_intestine - kon7b * c_memvasfr1_small_intestine * c_fcrnmemvas_small_intestine * vvm_small_intestine + koff7b * c_memvasfr2_small_intestine * vvm_small_intestine + clup_small_intestine * fr * c_endorecycfr1_small_intestine - clup_small_intestine * fr * c_memvasfr1_small_intestine - kdegab * c_memvasfr1_small_intestine * vvm_small_intestine
    d/dt(memvasfr1_large_intestine) <- kon7 * c_memvas_large_intestine * c_fcrnmemvas_large_intestine * vvm_large_intestine - koff7 * c_memvasfr1_large_intestine * vvm_large_intestine - kon7b * c_memvasfr1_large_intestine * c_fcrnmemvas_large_intestine * vvm_large_intestine + koff7b * c_memvasfr2_large_intestine * vvm_large_intestine + clup_large_intestine * fr * c_endorecycfr1_large_intestine - clup_large_intestine * fr * c_memvasfr1_large_intestine - kdegab * c_memvasfr1_large_intestine * vvm_large_intestine
    d/dt(memvasfr1_pancreas) <- kon7 * c_memvas_pancreas * c_fcrnmemvas_pancreas * vvm_pancreas - koff7 * c_memvasfr1_pancreas * vvm_pancreas - kon7b * c_memvasfr1_pancreas * c_fcrnmemvas_pancreas * vvm_pancreas + koff7b * c_memvasfr2_pancreas * vvm_pancreas + clup_pancreas * fr * c_endorecycfr1_pancreas - clup_pancreas * fr * c_memvasfr1_pancreas - kdegab * c_memvasfr1_pancreas * vvm_pancreas
    d/dt(memvasfr1_thymus) <- kon7 * c_memvas_thymus * c_fcrnmemvas_thymus * vvm_thymus - koff7 * c_memvasfr1_thymus * vvm_thymus - kon7b * c_memvasfr1_thymus * c_fcrnmemvas_thymus * vvm_thymus + koff7b * c_memvasfr2_thymus * vvm_thymus + clup_thymus * fr * c_endorecycfr1_thymus - clup_thymus * fr * c_memvasfr1_thymus - kdegab * c_memvasfr1_thymus * vvm_thymus
    d/dt(memvasfr1_spleen) <- kon7 * c_memvas_spleen * c_fcrnmemvas_spleen * vvm_spleen - koff7 * c_memvasfr1_spleen * vvm_spleen - kon7b * c_memvasfr1_spleen * c_fcrnmemvas_spleen * vvm_spleen + koff7b * c_memvasfr2_spleen * vvm_spleen + clup_spleen * fr * c_endorecycfr1_spleen - clup_spleen * fr * c_memvasfr1_spleen - kdegab * c_memvasfr1_spleen * vvm_spleen
    d/dt(memvasfr1_other) <- kon7 * c_memvas_other * c_fcrnmemvas_other * vvm_other - koff7 * c_memvasfr1_other * vvm_other - kon7b * c_memvasfr1_other * c_fcrnmemvas_other * vvm_other + koff7b * c_memvasfr2_other * vvm_other + clup_other * fr * c_endorecycfr1_other - clup_other * fr * c_memvasfr1_other - kdegab * c_memvasfr1_other * vvm_other
    d/dt(endoearlyfr1_lung) <- kon7 * c_endoearly_lung * c_fcrnendoearly_lung * ve7_lung - koff7 * c_endoearlyfr1_lung * ve7_lung - kon7b * c_endoearlyfr1_lung * c_fcrnendoearly_lung * ve7_lung + koff7b * c_endoearlyfr2_lung * ve7_lung - clup_lung * c_endoearlyfr1_lung + clup_lung * fr * c_memvasfr1_lung + clup_lung * (1 - fr) * c_memintfr1_lung
    d/dt(endoearlyfr1_liver) <- kon7 * c_endoearly_liver * c_fcrnendoearly_liver * ve7_liver - koff7 * c_endoearlyfr1_liver * ve7_liver - kon7b * c_endoearlyfr1_liver * c_fcrnendoearly_liver * ve7_liver + koff7b * c_endoearlyfr2_liver * ve7_liver - clup_liver * c_endoearlyfr1_liver + clup_liver * fr * c_memvasfr1_liver + clup_liver * (1 - fr) * c_memintfr1_liver
    d/dt(endoearlyfr1_heart) <- kon7 * c_endoearly_heart * c_fcrnendoearly_heart * ve7_heart - koff7 * c_endoearlyfr1_heart * ve7_heart - kon7b * c_endoearlyfr1_heart * c_fcrnendoearly_heart * ve7_heart + koff7b * c_endoearlyfr2_heart * ve7_heart - clup_heart * c_endoearlyfr1_heart + clup_heart * fr * c_memvasfr1_heart + clup_heart * (1 - fr) * c_memintfr1_heart
    d/dt(endoearlyfr1_muscle) <- kon7 * c_endoearly_muscle * c_fcrnendoearly_muscle * ve7_muscle - koff7 * c_endoearlyfr1_muscle * ve7_muscle - kon7b * c_endoearlyfr1_muscle * c_fcrnendoearly_muscle * ve7_muscle + koff7b * c_endoearlyfr2_muscle * ve7_muscle - clup_muscle * c_endoearlyfr1_muscle + clup_muscle * fr * c_memvasfr1_muscle + clup_muscle * (1 - fr) * c_memintfr1_muscle
    d/dt(endoearlyfr1_skin) <- kon7 * c_endoearly_skin * c_fcrnendoearly_skin * ve7_skin - koff7 * c_endoearlyfr1_skin * ve7_skin - kon7b * c_endoearlyfr1_skin * c_fcrnendoearly_skin * ve7_skin + koff7b * c_endoearlyfr2_skin * ve7_skin - clup_skin * c_endoearlyfr1_skin + clup_skin * fr * c_memvasfr1_skin + clup_skin * (1 - fr) * c_memintfr1_skin
    d/dt(endoearlyfr1_adipose) <- kon7 * c_endoearly_adipose * c_fcrnendoearly_adipose * ve7_adipose - koff7 * c_endoearlyfr1_adipose * ve7_adipose - kon7b * c_endoearlyfr1_adipose * c_fcrnendoearly_adipose * ve7_adipose + koff7b * c_endoearlyfr2_adipose * ve7_adipose - clup_adipose * c_endoearlyfr1_adipose + clup_adipose * fr * c_memvasfr1_adipose + clup_adipose * (1 - fr) * c_memintfr1_adipose
    d/dt(endoearlyfr1_bone) <- kon7 * c_endoearly_bone * c_fcrnendoearly_bone * ve7_bone - koff7 * c_endoearlyfr1_bone * ve7_bone - kon7b * c_endoearlyfr1_bone * c_fcrnendoearly_bone * ve7_bone + koff7b * c_endoearlyfr2_bone * ve7_bone - clup_bone * c_endoearlyfr1_bone + clup_bone * fr * c_memvasfr1_bone + clup_bone * (1 - fr) * c_memintfr1_bone
    d/dt(endoearlyfr1_brain) <- kon7 * c_endoearly_brain * c_fcrnendoearly_brain * ve7_brain - koff7 * c_endoearlyfr1_brain * ve7_brain - kon7b * c_endoearlyfr1_brain * c_fcrnendoearly_brain * ve7_brain + koff7b * c_endoearlyfr2_brain * ve7_brain - clup_brain * c_endoearlyfr1_brain + clup_brain * fr * c_memvasfr1_brain + clup_brain * (1 - fr) * c_memintfr1_brain
    d/dt(endoearlyfr1_kidney) <- kon7 * c_endoearly_kidney * c_fcrnendoearly_kidney * ve7_kidney - koff7 * c_endoearlyfr1_kidney * ve7_kidney - kon7b * c_endoearlyfr1_kidney * c_fcrnendoearly_kidney * ve7_kidney + koff7b * c_endoearlyfr2_kidney * ve7_kidney - clup_kidney * c_endoearlyfr1_kidney + clup_kidney * fr * c_memvasfr1_kidney + clup_kidney * (1 - fr) * c_memintfr1_kidney
    d/dt(endoearlyfr1_small_intestine) <- kon7 * c_endoearly_small_intestine * c_fcrnendoearly_small_intestine * ve7_small_intestine - koff7 * c_endoearlyfr1_small_intestine * ve7_small_intestine - kon7b * c_endoearlyfr1_small_intestine * c_fcrnendoearly_small_intestine * ve7_small_intestine + koff7b * c_endoearlyfr2_small_intestine * ve7_small_intestine - clup_small_intestine * c_endoearlyfr1_small_intestine + clup_small_intestine * fr * c_memvasfr1_small_intestine + clup_small_intestine * (1 - fr) * c_memintfr1_small_intestine
    d/dt(endoearlyfr1_large_intestine) <- kon7 * c_endoearly_large_intestine * c_fcrnendoearly_large_intestine * ve7_large_intestine - koff7 * c_endoearlyfr1_large_intestine * ve7_large_intestine - kon7b * c_endoearlyfr1_large_intestine * c_fcrnendoearly_large_intestine * ve7_large_intestine + koff7b * c_endoearlyfr2_large_intestine * ve7_large_intestine - clup_large_intestine * c_endoearlyfr1_large_intestine + clup_large_intestine * fr * c_memvasfr1_large_intestine + clup_large_intestine * (1 - fr) * c_memintfr1_large_intestine
    d/dt(endoearlyfr1_pancreas) <- kon7 * c_endoearly_pancreas * c_fcrnendoearly_pancreas * ve7_pancreas - koff7 * c_endoearlyfr1_pancreas * ve7_pancreas - kon7b * c_endoearlyfr1_pancreas * c_fcrnendoearly_pancreas * ve7_pancreas + koff7b * c_endoearlyfr2_pancreas * ve7_pancreas - clup_pancreas * c_endoearlyfr1_pancreas + clup_pancreas * fr * c_memvasfr1_pancreas + clup_pancreas * (1 - fr) * c_memintfr1_pancreas
    d/dt(endoearlyfr1_thymus) <- kon7 * c_endoearly_thymus * c_fcrnendoearly_thymus * ve7_thymus - koff7 * c_endoearlyfr1_thymus * ve7_thymus - kon7b * c_endoearlyfr1_thymus * c_fcrnendoearly_thymus * ve7_thymus + koff7b * c_endoearlyfr2_thymus * ve7_thymus - clup_thymus * c_endoearlyfr1_thymus + clup_thymus * fr * c_memvasfr1_thymus + clup_thymus * (1 - fr) * c_memintfr1_thymus
    d/dt(endoearlyfr1_spleen) <- kon7 * c_endoearly_spleen * c_fcrnendoearly_spleen * ve7_spleen - koff7 * c_endoearlyfr1_spleen * ve7_spleen - kon7b * c_endoearlyfr1_spleen * c_fcrnendoearly_spleen * ve7_spleen + koff7b * c_endoearlyfr2_spleen * ve7_spleen - clup_spleen * c_endoearlyfr1_spleen + clup_spleen * fr * c_memvasfr1_spleen + clup_spleen * (1 - fr) * c_memintfr1_spleen
    d/dt(endoearlyfr1_other) <- kon7 * c_endoearly_other * c_fcrnendoearly_other * ve7_other - koff7 * c_endoearlyfr1_other * ve7_other - kon7b * c_endoearlyfr1_other * c_fcrnendoearly_other * ve7_other + koff7b * c_endoearlyfr2_other * ve7_other - clup_other * c_endoearlyfr1_other + clup_other * fr * c_memvasfr1_other + clup_other * (1 - fr) * c_memintfr1_other
    d/dt(endosortfr1_lung) <- kon6 * c_endosort_lung * c_fcrnendosort_lung * ve6a_lung - koff6 * c_endosortfr1_lung * ve6a_lung - kon6b * c_endosortfr1_lung * c_fcrnendosort_lung * ve6a_lung + koff6b * c_endosortfr2_lung * ve6a_lung + clup_lung * c_endoearlyfr1_lung - clup_lung * c_endosortfr1_lung
    d/dt(endosortfr1_liver) <- kon6 * c_endosort_liver * c_fcrnendosort_liver * ve6a_liver - koff6 * c_endosortfr1_liver * ve6a_liver - kon6b * c_endosortfr1_liver * c_fcrnendosort_liver * ve6a_liver + koff6b * c_endosortfr2_liver * ve6a_liver + clup_liver * c_endoearlyfr1_liver - clup_liver * c_endosortfr1_liver
    d/dt(endosortfr1_heart) <- kon6 * c_endosort_heart * c_fcrnendosort_heart * ve6a_heart - koff6 * c_endosortfr1_heart * ve6a_heart - kon6b * c_endosortfr1_heart * c_fcrnendosort_heart * ve6a_heart + koff6b * c_endosortfr2_heart * ve6a_heart + clup_heart * c_endoearlyfr1_heart - clup_heart * c_endosortfr1_heart
    d/dt(endosortfr1_muscle) <- kon6 * c_endosort_muscle * c_fcrnendosort_muscle * ve6a_muscle - koff6 * c_endosortfr1_muscle * ve6a_muscle - kon6b * c_endosortfr1_muscle * c_fcrnendosort_muscle * ve6a_muscle + koff6b * c_endosortfr2_muscle * ve6a_muscle + clup_muscle * c_endoearlyfr1_muscle - clup_muscle * c_endosortfr1_muscle
    d/dt(endosortfr1_skin) <- kon6 * c_endosort_skin * c_fcrnendosort_skin * ve6a_skin - koff6 * c_endosortfr1_skin * ve6a_skin - kon6b * c_endosortfr1_skin * c_fcrnendosort_skin * ve6a_skin + koff6b * c_endosortfr2_skin * ve6a_skin + clup_skin * c_endoearlyfr1_skin - clup_skin * c_endosortfr1_skin
    d/dt(endosortfr1_adipose) <- kon6 * c_endosort_adipose * c_fcrnendosort_adipose * ve6a_adipose - koff6 * c_endosortfr1_adipose * ve6a_adipose - kon6b * c_endosortfr1_adipose * c_fcrnendosort_adipose * ve6a_adipose + koff6b * c_endosortfr2_adipose * ve6a_adipose + clup_adipose * c_endoearlyfr1_adipose - clup_adipose * c_endosortfr1_adipose
    d/dt(endosortfr1_bone) <- kon6 * c_endosort_bone * c_fcrnendosort_bone * ve6a_bone - koff6 * c_endosortfr1_bone * ve6a_bone - kon6b * c_endosortfr1_bone * c_fcrnendosort_bone * ve6a_bone + koff6b * c_endosortfr2_bone * ve6a_bone + clup_bone * c_endoearlyfr1_bone - clup_bone * c_endosortfr1_bone
    d/dt(endosortfr1_brain) <- kon6 * c_endosort_brain * c_fcrnendosort_brain * ve6a_brain - koff6 * c_endosortfr1_brain * ve6a_brain - kon6b * c_endosortfr1_brain * c_fcrnendosort_brain * ve6a_brain + koff6b * c_endosortfr2_brain * ve6a_brain + clup_brain * c_endoearlyfr1_brain - clup_brain * c_endosortfr1_brain
    d/dt(endosortfr1_kidney) <- kon6 * c_endosort_kidney * c_fcrnendosort_kidney * ve6a_kidney - koff6 * c_endosortfr1_kidney * ve6a_kidney - kon6b * c_endosortfr1_kidney * c_fcrnendosort_kidney * ve6a_kidney + koff6b * c_endosortfr2_kidney * ve6a_kidney + clup_kidney * c_endoearlyfr1_kidney - clup_kidney * c_endosortfr1_kidney
    d/dt(endosortfr1_small_intestine) <- kon6 * c_endosort_small_intestine * c_fcrnendosort_small_intestine * ve6a_small_intestine - koff6 * c_endosortfr1_small_intestine * ve6a_small_intestine - kon6b * c_endosortfr1_small_intestine * c_fcrnendosort_small_intestine * ve6a_small_intestine + koff6b * c_endosortfr2_small_intestine * ve6a_small_intestine + clup_small_intestine * c_endoearlyfr1_small_intestine - clup_small_intestine * c_endosortfr1_small_intestine
    d/dt(endosortfr1_large_intestine) <- kon6 * c_endosort_large_intestine * c_fcrnendosort_large_intestine * ve6a_large_intestine - koff6 * c_endosortfr1_large_intestine * ve6a_large_intestine - kon6b * c_endosortfr1_large_intestine * c_fcrnendosort_large_intestine * ve6a_large_intestine + koff6b * c_endosortfr2_large_intestine * ve6a_large_intestine + clup_large_intestine * c_endoearlyfr1_large_intestine - clup_large_intestine * c_endosortfr1_large_intestine
    d/dt(endosortfr1_pancreas) <- kon6 * c_endosort_pancreas * c_fcrnendosort_pancreas * ve6a_pancreas - koff6 * c_endosortfr1_pancreas * ve6a_pancreas - kon6b * c_endosortfr1_pancreas * c_fcrnendosort_pancreas * ve6a_pancreas + koff6b * c_endosortfr2_pancreas * ve6a_pancreas + clup_pancreas * c_endoearlyfr1_pancreas - clup_pancreas * c_endosortfr1_pancreas
    d/dt(endosortfr1_thymus) <- kon6 * c_endosort_thymus * c_fcrnendosort_thymus * ve6a_thymus - koff6 * c_endosortfr1_thymus * ve6a_thymus - kon6b * c_endosortfr1_thymus * c_fcrnendosort_thymus * ve6a_thymus + koff6b * c_endosortfr2_thymus * ve6a_thymus + clup_thymus * c_endoearlyfr1_thymus - clup_thymus * c_endosortfr1_thymus
    d/dt(endosortfr1_spleen) <- kon6 * c_endosort_spleen * c_fcrnendosort_spleen * ve6a_spleen - koff6 * c_endosortfr1_spleen * ve6a_spleen - kon6b * c_endosortfr1_spleen * c_fcrnendosort_spleen * ve6a_spleen + koff6b * c_endosortfr2_spleen * ve6a_spleen + clup_spleen * c_endoearlyfr1_spleen - clup_spleen * c_endosortfr1_spleen
    d/dt(endosortfr1_other) <- kon6 * c_endosort_other * c_fcrnendosort_other * ve6a_other - koff6 * c_endosortfr1_other * ve6a_other - kon6b * c_endosortfr1_other * c_fcrnendosort_other * ve6a_other + koff6b * c_endosortfr2_other * ve6a_other + clup_other * c_endoearlyfr1_other - clup_other * c_endosortfr1_other
    d/dt(endorecycfr1_lung) <- kon7 * c_endorecyc_lung * c_fcrnendorecyc_lung * ve7b_lung - koff7 * c_endorecycfr1_lung * ve7b_lung - kon7b * c_endorecycfr1_lung * c_fcrnendorecyc_lung * ve7b_lung + koff7b * c_endorecycfr2_lung * ve7b_lung + clup_lung * c_endosortfr1_lung - clup_lung * c_endorecycfr1_lung
    d/dt(endorecycfr1_liver) <- kon7 * c_endorecyc_liver * c_fcrnendorecyc_liver * ve7b_liver - koff7 * c_endorecycfr1_liver * ve7b_liver - kon7b * c_endorecycfr1_liver * c_fcrnendorecyc_liver * ve7b_liver + koff7b * c_endorecycfr2_liver * ve7b_liver + clup_liver * c_endosortfr1_liver - clup_liver * c_endorecycfr1_liver
    d/dt(endorecycfr1_heart) <- kon7 * c_endorecyc_heart * c_fcrnendorecyc_heart * ve7b_heart - koff7 * c_endorecycfr1_heart * ve7b_heart - kon7b * c_endorecycfr1_heart * c_fcrnendorecyc_heart * ve7b_heart + koff7b * c_endorecycfr2_heart * ve7b_heart + clup_heart * c_endosortfr1_heart - clup_heart * c_endorecycfr1_heart
    d/dt(endorecycfr1_muscle) <- kon7 * c_endorecyc_muscle * c_fcrnendorecyc_muscle * ve7b_muscle - koff7 * c_endorecycfr1_muscle * ve7b_muscle - kon7b * c_endorecycfr1_muscle * c_fcrnendorecyc_muscle * ve7b_muscle + koff7b * c_endorecycfr2_muscle * ve7b_muscle + clup_muscle * c_endosortfr1_muscle - clup_muscle * c_endorecycfr1_muscle
    d/dt(endorecycfr1_skin) <- kon7 * c_endorecyc_skin * c_fcrnendorecyc_skin * ve7b_skin - koff7 * c_endorecycfr1_skin * ve7b_skin - kon7b * c_endorecycfr1_skin * c_fcrnendorecyc_skin * ve7b_skin + koff7b * c_endorecycfr2_skin * ve7b_skin + clup_skin * c_endosortfr1_skin - clup_skin * c_endorecycfr1_skin
    d/dt(endorecycfr1_adipose) <- kon7 * c_endorecyc_adipose * c_fcrnendorecyc_adipose * ve7b_adipose - koff7 * c_endorecycfr1_adipose * ve7b_adipose - kon7b * c_endorecycfr1_adipose * c_fcrnendorecyc_adipose * ve7b_adipose + koff7b * c_endorecycfr2_adipose * ve7b_adipose + clup_adipose * c_endosortfr1_adipose - clup_adipose * c_endorecycfr1_adipose
    d/dt(endorecycfr1_bone) <- kon7 * c_endorecyc_bone * c_fcrnendorecyc_bone * ve7b_bone - koff7 * c_endorecycfr1_bone * ve7b_bone - kon7b * c_endorecycfr1_bone * c_fcrnendorecyc_bone * ve7b_bone + koff7b * c_endorecycfr2_bone * ve7b_bone + clup_bone * c_endosortfr1_bone - clup_bone * c_endorecycfr1_bone
    d/dt(endorecycfr1_brain) <- kon7 * c_endorecyc_brain * c_fcrnendorecyc_brain * ve7b_brain - koff7 * c_endorecycfr1_brain * ve7b_brain - kon7b * c_endorecycfr1_brain * c_fcrnendorecyc_brain * ve7b_brain + koff7b * c_endorecycfr2_brain * ve7b_brain + clup_brain * c_endosortfr1_brain - clup_brain * c_endorecycfr1_brain
    d/dt(endorecycfr1_kidney) <- kon7 * c_endorecyc_kidney * c_fcrnendorecyc_kidney * ve7b_kidney - koff7 * c_endorecycfr1_kidney * ve7b_kidney - kon7b * c_endorecycfr1_kidney * c_fcrnendorecyc_kidney * ve7b_kidney + koff7b * c_endorecycfr2_kidney * ve7b_kidney + clup_kidney * c_endosortfr1_kidney - clup_kidney * c_endorecycfr1_kidney
    d/dt(endorecycfr1_small_intestine) <- kon7 * c_endorecyc_small_intestine * c_fcrnendorecyc_small_intestine * ve7b_small_intestine - koff7 * c_endorecycfr1_small_intestine * ve7b_small_intestine - kon7b * c_endorecycfr1_small_intestine * c_fcrnendorecyc_small_intestine * ve7b_small_intestine + koff7b * c_endorecycfr2_small_intestine * ve7b_small_intestine + clup_small_intestine * c_endosortfr1_small_intestine - clup_small_intestine * c_endorecycfr1_small_intestine
    d/dt(endorecycfr1_large_intestine) <- kon7 * c_endorecyc_large_intestine * c_fcrnendorecyc_large_intestine * ve7b_large_intestine - koff7 * c_endorecycfr1_large_intestine * ve7b_large_intestine - kon7b * c_endorecycfr1_large_intestine * c_fcrnendorecyc_large_intestine * ve7b_large_intestine + koff7b * c_endorecycfr2_large_intestine * ve7b_large_intestine + clup_large_intestine * c_endosortfr1_large_intestine - clup_large_intestine * c_endorecycfr1_large_intestine
    d/dt(endorecycfr1_pancreas) <- kon7 * c_endorecyc_pancreas * c_fcrnendorecyc_pancreas * ve7b_pancreas - koff7 * c_endorecycfr1_pancreas * ve7b_pancreas - kon7b * c_endorecycfr1_pancreas * c_fcrnendorecyc_pancreas * ve7b_pancreas + koff7b * c_endorecycfr2_pancreas * ve7b_pancreas + clup_pancreas * c_endosortfr1_pancreas - clup_pancreas * c_endorecycfr1_pancreas
    d/dt(endorecycfr1_thymus) <- kon7 * c_endorecyc_thymus * c_fcrnendorecyc_thymus * ve7b_thymus - koff7 * c_endorecycfr1_thymus * ve7b_thymus - kon7b * c_endorecycfr1_thymus * c_fcrnendorecyc_thymus * ve7b_thymus + koff7b * c_endorecycfr2_thymus * ve7b_thymus + clup_thymus * c_endosortfr1_thymus - clup_thymus * c_endorecycfr1_thymus
    d/dt(endorecycfr1_spleen) <- kon7 * c_endorecyc_spleen * c_fcrnendorecyc_spleen * ve7b_spleen - koff7 * c_endorecycfr1_spleen * ve7b_spleen - kon7b * c_endorecycfr1_spleen * c_fcrnendorecyc_spleen * ve7b_spleen + koff7b * c_endorecycfr2_spleen * ve7b_spleen + clup_spleen * c_endosortfr1_spleen - clup_spleen * c_endorecycfr1_spleen
    d/dt(endorecycfr1_other) <- kon7 * c_endorecyc_other * c_fcrnendorecyc_other * ve7b_other - koff7 * c_endorecycfr1_other * ve7b_other - kon7b * c_endorecycfr1_other * c_fcrnendorecyc_other * ve7b_other + koff7b * c_endorecycfr2_other * ve7b_other + clup_other * c_endosortfr1_other - clup_other * c_endorecycfr1_other
    d/dt(memintfr1_lung) <- kon7 * c_memint_lung * c_fcrnmemint_lung * vism_lung - koff7 * c_memintfr1_lung * vism_lung - kon7b * c_memintfr1_lung * c_fcrnmemint_lung * vism_lung + koff7b * c_memintfr2_lung * vism_lung + clup_lung * (1 - fr) * c_endorecycfr1_lung - clup_lung * (1 - fr) * c_memintfr1_lung - kdegab * c_memintfr1_lung * vism_lung
    d/dt(memintfr1_liver) <- kon7 * c_memint_liver * c_fcrnmemint_liver * vism_liver - koff7 * c_memintfr1_liver * vism_liver - kon7b * c_memintfr1_liver * c_fcrnmemint_liver * vism_liver + koff7b * c_memintfr2_liver * vism_liver + clup_liver * (1 - fr) * c_endorecycfr1_liver - clup_liver * (1 - fr) * c_memintfr1_liver - kdegab * c_memintfr1_liver * vism_liver
    d/dt(memintfr1_heart) <- kon7 * c_memint_heart * c_fcrnmemint_heart * vism_heart - koff7 * c_memintfr1_heart * vism_heart - kon7b * c_memintfr1_heart * c_fcrnmemint_heart * vism_heart + koff7b * c_memintfr2_heart * vism_heart + clup_heart * (1 - fr) * c_endorecycfr1_heart - clup_heart * (1 - fr) * c_memintfr1_heart - kdegab * c_memintfr1_heart * vism_heart
    d/dt(memintfr1_muscle) <- kon7 * c_memint_muscle * c_fcrnmemint_muscle * vism_muscle - koff7 * c_memintfr1_muscle * vism_muscle - kon7b * c_memintfr1_muscle * c_fcrnmemint_muscle * vism_muscle + koff7b * c_memintfr2_muscle * vism_muscle + clup_muscle * (1 - fr) * c_endorecycfr1_muscle - clup_muscle * (1 - fr) * c_memintfr1_muscle - kdegab * c_memintfr1_muscle * vism_muscle
    d/dt(memintfr1_skin) <- kon7 * c_memint_skin * c_fcrnmemint_skin * vism_skin - koff7 * c_memintfr1_skin * vism_skin - kon7b * c_memintfr1_skin * c_fcrnmemint_skin * vism_skin + koff7b * c_memintfr2_skin * vism_skin + clup_skin * (1 - fr) * c_endorecycfr1_skin - clup_skin * (1 - fr) * c_memintfr1_skin - kdegab * c_memintfr1_skin * vism_skin
    d/dt(memintfr1_adipose) <- kon7 * c_memint_adipose * c_fcrnmemint_adipose * vism_adipose - koff7 * c_memintfr1_adipose * vism_adipose - kon7b * c_memintfr1_adipose * c_fcrnmemint_adipose * vism_adipose + koff7b * c_memintfr2_adipose * vism_adipose + clup_adipose * (1 - fr) * c_endorecycfr1_adipose - clup_adipose * (1 - fr) * c_memintfr1_adipose - kdegab * c_memintfr1_adipose * vism_adipose
    d/dt(memintfr1_bone) <- kon7 * c_memint_bone * c_fcrnmemint_bone * vism_bone - koff7 * c_memintfr1_bone * vism_bone - kon7b * c_memintfr1_bone * c_fcrnmemint_bone * vism_bone + koff7b * c_memintfr2_bone * vism_bone + clup_bone * (1 - fr) * c_endorecycfr1_bone - clup_bone * (1 - fr) * c_memintfr1_bone - kdegab * c_memintfr1_bone * vism_bone
    d/dt(memintfr1_brain) <- kon7 * c_memint_brain * c_fcrnmemint_brain * vism_brain - koff7 * c_memintfr1_brain * vism_brain - kon7b * c_memintfr1_brain * c_fcrnmemint_brain * vism_brain + koff7b * c_memintfr2_brain * vism_brain + clup_brain * (1 - fr) * c_endorecycfr1_brain - clup_brain * (1 - fr) * c_memintfr1_brain - kdegab * c_memintfr1_brain * vism_brain
    d/dt(memintfr1_kidney) <- kon7 * c_memint_kidney * c_fcrnmemint_kidney * vism_kidney - koff7 * c_memintfr1_kidney * vism_kidney - kon7b * c_memintfr1_kidney * c_fcrnmemint_kidney * vism_kidney + koff7b * c_memintfr2_kidney * vism_kidney + clup_kidney * (1 - fr) * c_endorecycfr1_kidney - clup_kidney * (1 - fr) * c_memintfr1_kidney - kdegab * c_memintfr1_kidney * vism_kidney
    d/dt(memintfr1_small_intestine) <- kon7 * c_memint_small_intestine * c_fcrnmemint_small_intestine * vism_small_intestine - koff7 * c_memintfr1_small_intestine * vism_small_intestine - kon7b * c_memintfr1_small_intestine * c_fcrnmemint_small_intestine * vism_small_intestine + koff7b * c_memintfr2_small_intestine * vism_small_intestine + clup_small_intestine * (1 - fr) * c_endorecycfr1_small_intestine - clup_small_intestine * (1 - fr) * c_memintfr1_small_intestine - kdegab * c_memintfr1_small_intestine * vism_small_intestine
    d/dt(memintfr1_large_intestine) <- kon7 * c_memint_large_intestine * c_fcrnmemint_large_intestine * vism_large_intestine - koff7 * c_memintfr1_large_intestine * vism_large_intestine - kon7b * c_memintfr1_large_intestine * c_fcrnmemint_large_intestine * vism_large_intestine + koff7b * c_memintfr2_large_intestine * vism_large_intestine + clup_large_intestine * (1 - fr) * c_endorecycfr1_large_intestine - clup_large_intestine * (1 - fr) * c_memintfr1_large_intestine - kdegab * c_memintfr1_large_intestine * vism_large_intestine
    d/dt(memintfr1_pancreas) <- kon7 * c_memint_pancreas * c_fcrnmemint_pancreas * vism_pancreas - koff7 * c_memintfr1_pancreas * vism_pancreas - kon7b * c_memintfr1_pancreas * c_fcrnmemint_pancreas * vism_pancreas + koff7b * c_memintfr2_pancreas * vism_pancreas + clup_pancreas * (1 - fr) * c_endorecycfr1_pancreas - clup_pancreas * (1 - fr) * c_memintfr1_pancreas - kdegab * c_memintfr1_pancreas * vism_pancreas
    d/dt(memintfr1_thymus) <- kon7 * c_memint_thymus * c_fcrnmemint_thymus * vism_thymus - koff7 * c_memintfr1_thymus * vism_thymus - kon7b * c_memintfr1_thymus * c_fcrnmemint_thymus * vism_thymus + koff7b * c_memintfr2_thymus * vism_thymus + clup_thymus * (1 - fr) * c_endorecycfr1_thymus - clup_thymus * (1 - fr) * c_memintfr1_thymus - kdegab * c_memintfr1_thymus * vism_thymus
    d/dt(memintfr1_spleen) <- kon7 * c_memint_spleen * c_fcrnmemint_spleen * vism_spleen - koff7 * c_memintfr1_spleen * vism_spleen - kon7b * c_memintfr1_spleen * c_fcrnmemint_spleen * vism_spleen + koff7b * c_memintfr2_spleen * vism_spleen + clup_spleen * (1 - fr) * c_endorecycfr1_spleen - clup_spleen * (1 - fr) * c_memintfr1_spleen - kdegab * c_memintfr1_spleen * vism_spleen
    d/dt(memintfr1_other) <- kon7 * c_memint_other * c_fcrnmemint_other * vism_other - koff7 * c_memintfr1_other * vism_other - kon7b * c_memintfr1_other * c_fcrnmemint_other * vism_other + koff7b * c_memintfr2_other * vism_other + clup_other * (1 - fr) * c_endorecycfr1_other - clup_other * (1 - fr) * c_memintfr1_other - kdegab * c_memintfr1_other * vism_other

    # ---- mab: 2:1 FcRn complexes ----
    d/dt(memvasfr2_lung) <- kon7b * c_memvasfr1_lung * c_fcrnmemvas_lung * vvm_lung - koff7b * c_memvasfr2_lung * vvm_lung + clup_lung * fr * c_endorecycfr2_lung - clup_lung * fr * c_memvasfr2_lung - kdegab * c_memvasfr2_lung * vvm_lung
    d/dt(memvasfr2_liver) <- kon7b * c_memvasfr1_liver * c_fcrnmemvas_liver * vvm_liver - koff7b * c_memvasfr2_liver * vvm_liver + clup_liver * fr * c_endorecycfr2_liver - clup_liver * fr * c_memvasfr2_liver - kdegab * c_memvasfr2_liver * vvm_liver
    d/dt(memvasfr2_heart) <- kon7b * c_memvasfr1_heart * c_fcrnmemvas_heart * vvm_heart - koff7b * c_memvasfr2_heart * vvm_heart + clup_heart * fr * c_endorecycfr2_heart - clup_heart * fr * c_memvasfr2_heart - kdegab * c_memvasfr2_heart * vvm_heart
    d/dt(memvasfr2_muscle) <- kon7b * c_memvasfr1_muscle * c_fcrnmemvas_muscle * vvm_muscle - koff7b * c_memvasfr2_muscle * vvm_muscle + clup_muscle * fr * c_endorecycfr2_muscle - clup_muscle * fr * c_memvasfr2_muscle - kdegab * c_memvasfr2_muscle * vvm_muscle
    d/dt(memvasfr2_skin) <- kon7b * c_memvasfr1_skin * c_fcrnmemvas_skin * vvm_skin - koff7b * c_memvasfr2_skin * vvm_skin + clup_skin * fr * c_endorecycfr2_skin - clup_skin * fr * c_memvasfr2_skin - kdegab * c_memvasfr2_skin * vvm_skin
    d/dt(memvasfr2_adipose) <- kon7b * c_memvasfr1_adipose * c_fcrnmemvas_adipose * vvm_adipose - koff7b * c_memvasfr2_adipose * vvm_adipose + clup_adipose * fr * c_endorecycfr2_adipose - clup_adipose * fr * c_memvasfr2_adipose - kdegab * c_memvasfr2_adipose * vvm_adipose
    d/dt(memvasfr2_bone) <- kon7b * c_memvasfr1_bone * c_fcrnmemvas_bone * vvm_bone - koff7b * c_memvasfr2_bone * vvm_bone + clup_bone * fr * c_endorecycfr2_bone - clup_bone * fr * c_memvasfr2_bone - kdegab * c_memvasfr2_bone * vvm_bone
    d/dt(memvasfr2_brain) <- kon7b * c_memvasfr1_brain * c_fcrnmemvas_brain * vvm_brain - koff7b * c_memvasfr2_brain * vvm_brain + clup_brain * fr * c_endorecycfr2_brain - clup_brain * fr * c_memvasfr2_brain - kdegab * c_memvasfr2_brain * vvm_brain
    d/dt(memvasfr2_kidney) <- kon7b * c_memvasfr1_kidney * c_fcrnmemvas_kidney * vvm_kidney - koff7b * c_memvasfr2_kidney * vvm_kidney + clup_kidney * fr * c_endorecycfr2_kidney - clup_kidney * fr * c_memvasfr2_kidney - kdegab * c_memvasfr2_kidney * vvm_kidney
    d/dt(memvasfr2_small_intestine) <- kon7b * c_memvasfr1_small_intestine * c_fcrnmemvas_small_intestine * vvm_small_intestine - koff7b * c_memvasfr2_small_intestine * vvm_small_intestine + clup_small_intestine * fr * c_endorecycfr2_small_intestine - clup_small_intestine * fr * c_memvasfr2_small_intestine - kdegab * c_memvasfr2_small_intestine * vvm_small_intestine
    d/dt(memvasfr2_large_intestine) <- kon7b * c_memvasfr1_large_intestine * c_fcrnmemvas_large_intestine * vvm_large_intestine - koff7b * c_memvasfr2_large_intestine * vvm_large_intestine + clup_large_intestine * fr * c_endorecycfr2_large_intestine - clup_large_intestine * fr * c_memvasfr2_large_intestine - kdegab * c_memvasfr2_large_intestine * vvm_large_intestine
    d/dt(memvasfr2_pancreas) <- kon7b * c_memvasfr1_pancreas * c_fcrnmemvas_pancreas * vvm_pancreas - koff7b * c_memvasfr2_pancreas * vvm_pancreas + clup_pancreas * fr * c_endorecycfr2_pancreas - clup_pancreas * fr * c_memvasfr2_pancreas - kdegab * c_memvasfr2_pancreas * vvm_pancreas
    d/dt(memvasfr2_thymus) <- kon7b * c_memvasfr1_thymus * c_fcrnmemvas_thymus * vvm_thymus - koff7b * c_memvasfr2_thymus * vvm_thymus + clup_thymus * fr * c_endorecycfr2_thymus - clup_thymus * fr * c_memvasfr2_thymus - kdegab * c_memvasfr2_thymus * vvm_thymus
    d/dt(memvasfr2_spleen) <- kon7b * c_memvasfr1_spleen * c_fcrnmemvas_spleen * vvm_spleen - koff7b * c_memvasfr2_spleen * vvm_spleen + clup_spleen * fr * c_endorecycfr2_spleen - clup_spleen * fr * c_memvasfr2_spleen - kdegab * c_memvasfr2_spleen * vvm_spleen
    d/dt(memvasfr2_other) <- kon7b * c_memvasfr1_other * c_fcrnmemvas_other * vvm_other - koff7b * c_memvasfr2_other * vvm_other + clup_other * fr * c_endorecycfr2_other - clup_other * fr * c_memvasfr2_other - kdegab * c_memvasfr2_other * vvm_other
    d/dt(endoearlyfr2_lung) <- kon7b * c_endoearlyfr1_lung * c_fcrnendoearly_lung * ve7_lung - koff7b * c_endoearlyfr2_lung * ve7_lung - clup_lung * c_endoearlyfr2_lung + clup_lung * fr * c_memvasfr2_lung + clup_lung * (1 - fr) * c_memintfr2_lung
    d/dt(endoearlyfr2_liver) <- kon7b * c_endoearlyfr1_liver * c_fcrnendoearly_liver * ve7_liver - koff7b * c_endoearlyfr2_liver * ve7_liver - clup_liver * c_endoearlyfr2_liver + clup_liver * fr * c_memvasfr2_liver + clup_liver * (1 - fr) * c_memintfr2_liver
    d/dt(endoearlyfr2_heart) <- kon7b * c_endoearlyfr1_heart * c_fcrnendoearly_heart * ve7_heart - koff7b * c_endoearlyfr2_heart * ve7_heart - clup_heart * c_endoearlyfr2_heart + clup_heart * fr * c_memvasfr2_heart + clup_heart * (1 - fr) * c_memintfr2_heart
    d/dt(endoearlyfr2_muscle) <- kon7b * c_endoearlyfr1_muscle * c_fcrnendoearly_muscle * ve7_muscle - koff7b * c_endoearlyfr2_muscle * ve7_muscle - clup_muscle * c_endoearlyfr2_muscle + clup_muscle * fr * c_memvasfr2_muscle + clup_muscle * (1 - fr) * c_memintfr2_muscle
    d/dt(endoearlyfr2_skin) <- kon7b * c_endoearlyfr1_skin * c_fcrnendoearly_skin * ve7_skin - koff7b * c_endoearlyfr2_skin * ve7_skin - clup_skin * c_endoearlyfr2_skin + clup_skin * fr * c_memvasfr2_skin + clup_skin * (1 - fr) * c_memintfr2_skin
    d/dt(endoearlyfr2_adipose) <- kon7b * c_endoearlyfr1_adipose * c_fcrnendoearly_adipose * ve7_adipose - koff7b * c_endoearlyfr2_adipose * ve7_adipose - clup_adipose * c_endoearlyfr2_adipose + clup_adipose * fr * c_memvasfr2_adipose + clup_adipose * (1 - fr) * c_memintfr2_adipose
    d/dt(endoearlyfr2_bone) <- kon7b * c_endoearlyfr1_bone * c_fcrnendoearly_bone * ve7_bone - koff7b * c_endoearlyfr2_bone * ve7_bone - clup_bone * c_endoearlyfr2_bone + clup_bone * fr * c_memvasfr2_bone + clup_bone * (1 - fr) * c_memintfr2_bone
    d/dt(endoearlyfr2_brain) <- kon7b * c_endoearlyfr1_brain * c_fcrnendoearly_brain * ve7_brain - koff7b * c_endoearlyfr2_brain * ve7_brain - clup_brain * c_endoearlyfr2_brain + clup_brain * fr * c_memvasfr2_brain + clup_brain * (1 - fr) * c_memintfr2_brain
    d/dt(endoearlyfr2_kidney) <- kon7b * c_endoearlyfr1_kidney * c_fcrnendoearly_kidney * ve7_kidney - koff7b * c_endoearlyfr2_kidney * ve7_kidney - clup_kidney * c_endoearlyfr2_kidney + clup_kidney * fr * c_memvasfr2_kidney + clup_kidney * (1 - fr) * c_memintfr2_kidney
    d/dt(endoearlyfr2_small_intestine) <- kon7b * c_endoearlyfr1_small_intestine * c_fcrnendoearly_small_intestine * ve7_small_intestine - koff7b * c_endoearlyfr2_small_intestine * ve7_small_intestine - clup_small_intestine * c_endoearlyfr2_small_intestine + clup_small_intestine * fr * c_memvasfr2_small_intestine + clup_small_intestine * (1 - fr) * c_memintfr2_small_intestine
    d/dt(endoearlyfr2_large_intestine) <- kon7b * c_endoearlyfr1_large_intestine * c_fcrnendoearly_large_intestine * ve7_large_intestine - koff7b * c_endoearlyfr2_large_intestine * ve7_large_intestine - clup_large_intestine * c_endoearlyfr2_large_intestine + clup_large_intestine * fr * c_memvasfr2_large_intestine + clup_large_intestine * (1 - fr) * c_memintfr2_large_intestine
    d/dt(endoearlyfr2_pancreas) <- kon7b * c_endoearlyfr1_pancreas * c_fcrnendoearly_pancreas * ve7_pancreas - koff7b * c_endoearlyfr2_pancreas * ve7_pancreas - clup_pancreas * c_endoearlyfr2_pancreas + clup_pancreas * fr * c_memvasfr2_pancreas + clup_pancreas * (1 - fr) * c_memintfr2_pancreas
    d/dt(endoearlyfr2_thymus) <- kon7b * c_endoearlyfr1_thymus * c_fcrnendoearly_thymus * ve7_thymus - koff7b * c_endoearlyfr2_thymus * ve7_thymus - clup_thymus * c_endoearlyfr2_thymus + clup_thymus * fr * c_memvasfr2_thymus + clup_thymus * (1 - fr) * c_memintfr2_thymus
    d/dt(endoearlyfr2_spleen) <- kon7b * c_endoearlyfr1_spleen * c_fcrnendoearly_spleen * ve7_spleen - koff7b * c_endoearlyfr2_spleen * ve7_spleen - clup_spleen * c_endoearlyfr2_spleen + clup_spleen * fr * c_memvasfr2_spleen + clup_spleen * (1 - fr) * c_memintfr2_spleen
    d/dt(endoearlyfr2_other) <- kon7b * c_endoearlyfr1_other * c_fcrnendoearly_other * ve7_other - koff7b * c_endoearlyfr2_other * ve7_other - clup_other * c_endoearlyfr2_other + clup_other * fr * c_memvasfr2_other + clup_other * (1 - fr) * c_memintfr2_other
    d/dt(endosortfr2_lung) <- kon6b * c_endosortfr1_lung * c_fcrnendosort_lung * ve6a_lung - koff6b * c_endosortfr2_lung * ve6a_lung + clup_lung * c_endoearlyfr2_lung - clup_lung * c_endosortfr2_lung
    d/dt(endosortfr2_liver) <- kon6b * c_endosortfr1_liver * c_fcrnendosort_liver * ve6a_liver - koff6b * c_endosortfr2_liver * ve6a_liver + clup_liver * c_endoearlyfr2_liver - clup_liver * c_endosortfr2_liver
    d/dt(endosortfr2_heart) <- kon6b * c_endosortfr1_heart * c_fcrnendosort_heart * ve6a_heart - koff6b * c_endosortfr2_heart * ve6a_heart + clup_heart * c_endoearlyfr2_heart - clup_heart * c_endosortfr2_heart
    d/dt(endosortfr2_muscle) <- kon6b * c_endosortfr1_muscle * c_fcrnendosort_muscle * ve6a_muscle - koff6b * c_endosortfr2_muscle * ve6a_muscle + clup_muscle * c_endoearlyfr2_muscle - clup_muscle * c_endosortfr2_muscle
    d/dt(endosortfr2_skin) <- kon6b * c_endosortfr1_skin * c_fcrnendosort_skin * ve6a_skin - koff6b * c_endosortfr2_skin * ve6a_skin + clup_skin * c_endoearlyfr2_skin - clup_skin * c_endosortfr2_skin
    d/dt(endosortfr2_adipose) <- kon6b * c_endosortfr1_adipose * c_fcrnendosort_adipose * ve6a_adipose - koff6b * c_endosortfr2_adipose * ve6a_adipose + clup_adipose * c_endoearlyfr2_adipose - clup_adipose * c_endosortfr2_adipose
    d/dt(endosortfr2_bone) <- kon6b * c_endosortfr1_bone * c_fcrnendosort_bone * ve6a_bone - koff6b * c_endosortfr2_bone * ve6a_bone + clup_bone * c_endoearlyfr2_bone - clup_bone * c_endosortfr2_bone
    d/dt(endosortfr2_brain) <- kon6b * c_endosortfr1_brain * c_fcrnendosort_brain * ve6a_brain - koff6b * c_endosortfr2_brain * ve6a_brain + clup_brain * c_endoearlyfr2_brain - clup_brain * c_endosortfr2_brain
    d/dt(endosortfr2_kidney) <- kon6b * c_endosortfr1_kidney * c_fcrnendosort_kidney * ve6a_kidney - koff6b * c_endosortfr2_kidney * ve6a_kidney + clup_kidney * c_endoearlyfr2_kidney - clup_kidney * c_endosortfr2_kidney
    d/dt(endosortfr2_small_intestine) <- kon6b * c_endosortfr1_small_intestine * c_fcrnendosort_small_intestine * ve6a_small_intestine - koff6b * c_endosortfr2_small_intestine * ve6a_small_intestine + clup_small_intestine * c_endoearlyfr2_small_intestine - clup_small_intestine * c_endosortfr2_small_intestine
    d/dt(endosortfr2_large_intestine) <- kon6b * c_endosortfr1_large_intestine * c_fcrnendosort_large_intestine * ve6a_large_intestine - koff6b * c_endosortfr2_large_intestine * ve6a_large_intestine + clup_large_intestine * c_endoearlyfr2_large_intestine - clup_large_intestine * c_endosortfr2_large_intestine
    d/dt(endosortfr2_pancreas) <- kon6b * c_endosortfr1_pancreas * c_fcrnendosort_pancreas * ve6a_pancreas - koff6b * c_endosortfr2_pancreas * ve6a_pancreas + clup_pancreas * c_endoearlyfr2_pancreas - clup_pancreas * c_endosortfr2_pancreas
    d/dt(endosortfr2_thymus) <- kon6b * c_endosortfr1_thymus * c_fcrnendosort_thymus * ve6a_thymus - koff6b * c_endosortfr2_thymus * ve6a_thymus + clup_thymus * c_endoearlyfr2_thymus - clup_thymus * c_endosortfr2_thymus
    d/dt(endosortfr2_spleen) <- kon6b * c_endosortfr1_spleen * c_fcrnendosort_spleen * ve6a_spleen - koff6b * c_endosortfr2_spleen * ve6a_spleen + clup_spleen * c_endoearlyfr2_spleen - clup_spleen * c_endosortfr2_spleen
    d/dt(endosortfr2_other) <- kon6b * c_endosortfr1_other * c_fcrnendosort_other * ve6a_other - koff6b * c_endosortfr2_other * ve6a_other + clup_other * c_endoearlyfr2_other - clup_other * c_endosortfr2_other
    d/dt(endorecycfr2_lung) <- kon7b * c_endorecycfr1_lung * c_fcrnendorecyc_lung * ve7b_lung - koff7b * c_endorecycfr2_lung * ve7b_lung + clup_lung * c_endosortfr2_lung - clup_lung * c_endorecycfr2_lung
    d/dt(endorecycfr2_liver) <- kon7b * c_endorecycfr1_liver * c_fcrnendorecyc_liver * ve7b_liver - koff7b * c_endorecycfr2_liver * ve7b_liver + clup_liver * c_endosortfr2_liver - clup_liver * c_endorecycfr2_liver
    d/dt(endorecycfr2_heart) <- kon7b * c_endorecycfr1_heart * c_fcrnendorecyc_heart * ve7b_heart - koff7b * c_endorecycfr2_heart * ve7b_heart + clup_heart * c_endosortfr2_heart - clup_heart * c_endorecycfr2_heart
    d/dt(endorecycfr2_muscle) <- kon7b * c_endorecycfr1_muscle * c_fcrnendorecyc_muscle * ve7b_muscle - koff7b * c_endorecycfr2_muscle * ve7b_muscle + clup_muscle * c_endosortfr2_muscle - clup_muscle * c_endorecycfr2_muscle
    d/dt(endorecycfr2_skin) <- kon7b * c_endorecycfr1_skin * c_fcrnendorecyc_skin * ve7b_skin - koff7b * c_endorecycfr2_skin * ve7b_skin + clup_skin * c_endosortfr2_skin - clup_skin * c_endorecycfr2_skin
    d/dt(endorecycfr2_adipose) <- kon7b * c_endorecycfr1_adipose * c_fcrnendorecyc_adipose * ve7b_adipose - koff7b * c_endorecycfr2_adipose * ve7b_adipose + clup_adipose * c_endosortfr2_adipose - clup_adipose * c_endorecycfr2_adipose
    d/dt(endorecycfr2_bone) <- kon7b * c_endorecycfr1_bone * c_fcrnendorecyc_bone * ve7b_bone - koff7b * c_endorecycfr2_bone * ve7b_bone + clup_bone * c_endosortfr2_bone - clup_bone * c_endorecycfr2_bone
    d/dt(endorecycfr2_brain) <- kon7b * c_endorecycfr1_brain * c_fcrnendorecyc_brain * ve7b_brain - koff7b * c_endorecycfr2_brain * ve7b_brain + clup_brain * c_endosortfr2_brain - clup_brain * c_endorecycfr2_brain
    d/dt(endorecycfr2_kidney) <- kon7b * c_endorecycfr1_kidney * c_fcrnendorecyc_kidney * ve7b_kidney - koff7b * c_endorecycfr2_kidney * ve7b_kidney + clup_kidney * c_endosortfr2_kidney - clup_kidney * c_endorecycfr2_kidney
    d/dt(endorecycfr2_small_intestine) <- kon7b * c_endorecycfr1_small_intestine * c_fcrnendorecyc_small_intestine * ve7b_small_intestine - koff7b * c_endorecycfr2_small_intestine * ve7b_small_intestine + clup_small_intestine * c_endosortfr2_small_intestine - clup_small_intestine * c_endorecycfr2_small_intestine
    d/dt(endorecycfr2_large_intestine) <- kon7b * c_endorecycfr1_large_intestine * c_fcrnendorecyc_large_intestine * ve7b_large_intestine - koff7b * c_endorecycfr2_large_intestine * ve7b_large_intestine + clup_large_intestine * c_endosortfr2_large_intestine - clup_large_intestine * c_endorecycfr2_large_intestine
    d/dt(endorecycfr2_pancreas) <- kon7b * c_endorecycfr1_pancreas * c_fcrnendorecyc_pancreas * ve7b_pancreas - koff7b * c_endorecycfr2_pancreas * ve7b_pancreas + clup_pancreas * c_endosortfr2_pancreas - clup_pancreas * c_endorecycfr2_pancreas
    d/dt(endorecycfr2_thymus) <- kon7b * c_endorecycfr1_thymus * c_fcrnendorecyc_thymus * ve7b_thymus - koff7b * c_endorecycfr2_thymus * ve7b_thymus + clup_thymus * c_endosortfr2_thymus - clup_thymus * c_endorecycfr2_thymus
    d/dt(endorecycfr2_spleen) <- kon7b * c_endorecycfr1_spleen * c_fcrnendorecyc_spleen * ve7b_spleen - koff7b * c_endorecycfr2_spleen * ve7b_spleen + clup_spleen * c_endosortfr2_spleen - clup_spleen * c_endorecycfr2_spleen
    d/dt(endorecycfr2_other) <- kon7b * c_endorecycfr1_other * c_fcrnendorecyc_other * ve7b_other - koff7b * c_endorecycfr2_other * ve7b_other + clup_other * c_endosortfr2_other - clup_other * c_endorecycfr2_other
    d/dt(memintfr2_lung) <- kon7b * c_memintfr1_lung * c_fcrnmemint_lung * vism_lung - koff7b * c_memintfr2_lung * vism_lung + clup_lung * (1 - fr) * c_endorecycfr2_lung - clup_lung * (1 - fr) * c_memintfr2_lung - kdegab * c_memintfr2_lung * vism_lung
    d/dt(memintfr2_liver) <- kon7b * c_memintfr1_liver * c_fcrnmemint_liver * vism_liver - koff7b * c_memintfr2_liver * vism_liver + clup_liver * (1 - fr) * c_endorecycfr2_liver - clup_liver * (1 - fr) * c_memintfr2_liver - kdegab * c_memintfr2_liver * vism_liver
    d/dt(memintfr2_heart) <- kon7b * c_memintfr1_heart * c_fcrnmemint_heart * vism_heart - koff7b * c_memintfr2_heart * vism_heart + clup_heart * (1 - fr) * c_endorecycfr2_heart - clup_heart * (1 - fr) * c_memintfr2_heart - kdegab * c_memintfr2_heart * vism_heart
    d/dt(memintfr2_muscle) <- kon7b * c_memintfr1_muscle * c_fcrnmemint_muscle * vism_muscle - koff7b * c_memintfr2_muscle * vism_muscle + clup_muscle * (1 - fr) * c_endorecycfr2_muscle - clup_muscle * (1 - fr) * c_memintfr2_muscle - kdegab * c_memintfr2_muscle * vism_muscle
    d/dt(memintfr2_skin) <- kon7b * c_memintfr1_skin * c_fcrnmemint_skin * vism_skin - koff7b * c_memintfr2_skin * vism_skin + clup_skin * (1 - fr) * c_endorecycfr2_skin - clup_skin * (1 - fr) * c_memintfr2_skin - kdegab * c_memintfr2_skin * vism_skin
    d/dt(memintfr2_adipose) <- kon7b * c_memintfr1_adipose * c_fcrnmemint_adipose * vism_adipose - koff7b * c_memintfr2_adipose * vism_adipose + clup_adipose * (1 - fr) * c_endorecycfr2_adipose - clup_adipose * (1 - fr) * c_memintfr2_adipose - kdegab * c_memintfr2_adipose * vism_adipose
    d/dt(memintfr2_bone) <- kon7b * c_memintfr1_bone * c_fcrnmemint_bone * vism_bone - koff7b * c_memintfr2_bone * vism_bone + clup_bone * (1 - fr) * c_endorecycfr2_bone - clup_bone * (1 - fr) * c_memintfr2_bone - kdegab * c_memintfr2_bone * vism_bone
    d/dt(memintfr2_brain) <- kon7b * c_memintfr1_brain * c_fcrnmemint_brain * vism_brain - koff7b * c_memintfr2_brain * vism_brain + clup_brain * (1 - fr) * c_endorecycfr2_brain - clup_brain * (1 - fr) * c_memintfr2_brain - kdegab * c_memintfr2_brain * vism_brain
    d/dt(memintfr2_kidney) <- kon7b * c_memintfr1_kidney * c_fcrnmemint_kidney * vism_kidney - koff7b * c_memintfr2_kidney * vism_kidney + clup_kidney * (1 - fr) * c_endorecycfr2_kidney - clup_kidney * (1 - fr) * c_memintfr2_kidney - kdegab * c_memintfr2_kidney * vism_kidney
    d/dt(memintfr2_small_intestine) <- kon7b * c_memintfr1_small_intestine * c_fcrnmemint_small_intestine * vism_small_intestine - koff7b * c_memintfr2_small_intestine * vism_small_intestine + clup_small_intestine * (1 - fr) * c_endorecycfr2_small_intestine - clup_small_intestine * (1 - fr) * c_memintfr2_small_intestine - kdegab * c_memintfr2_small_intestine * vism_small_intestine
    d/dt(memintfr2_large_intestine) <- kon7b * c_memintfr1_large_intestine * c_fcrnmemint_large_intestine * vism_large_intestine - koff7b * c_memintfr2_large_intestine * vism_large_intestine + clup_large_intestine * (1 - fr) * c_endorecycfr2_large_intestine - clup_large_intestine * (1 - fr) * c_memintfr2_large_intestine - kdegab * c_memintfr2_large_intestine * vism_large_intestine
    d/dt(memintfr2_pancreas) <- kon7b * c_memintfr1_pancreas * c_fcrnmemint_pancreas * vism_pancreas - koff7b * c_memintfr2_pancreas * vism_pancreas + clup_pancreas * (1 - fr) * c_endorecycfr2_pancreas - clup_pancreas * (1 - fr) * c_memintfr2_pancreas - kdegab * c_memintfr2_pancreas * vism_pancreas
    d/dt(memintfr2_thymus) <- kon7b * c_memintfr1_thymus * c_fcrnmemint_thymus * vism_thymus - koff7b * c_memintfr2_thymus * vism_thymus + clup_thymus * (1 - fr) * c_endorecycfr2_thymus - clup_thymus * (1 - fr) * c_memintfr2_thymus - kdegab * c_memintfr2_thymus * vism_thymus
    d/dt(memintfr2_spleen) <- kon7b * c_memintfr1_spleen * c_fcrnmemint_spleen * vism_spleen - koff7b * c_memintfr2_spleen * vism_spleen + clup_spleen * (1 - fr) * c_endorecycfr2_spleen - clup_spleen * (1 - fr) * c_memintfr2_spleen - kdegab * c_memintfr2_spleen * vism_spleen
    d/dt(memintfr2_other) <- kon7b * c_memintfr1_other * c_fcrnmemint_other * vism_other - koff7b * c_memintfr2_other * vism_other + clup_other * (1 - fr) * c_endorecycfr2_other - clup_other * (1 - fr) * c_memintfr2_other - kdegab * c_memintfr2_other * vism_other

    # ---- Dosed antibody: non-specific (AC-SINS driven) membrane binding ----
    d/dt(memvasns_lung) <- konps * c_memvas_lung * cmem * vvm_lung - koffps * c_memvasns_lung * vvm_lung - kintps * c_memvasns_lung * vvm_lung
    d/dt(memvasns_liver) <- konps * c_memvas_liver * cmem * vvm_liver - koffps * c_memvasns_liver * vvm_liver - kintps * c_memvasns_liver * vvm_liver
    d/dt(memvasns_heart) <- konps * c_memvas_heart * cmem * vvm_heart - koffps * c_memvasns_heart * vvm_heart - kintps * c_memvasns_heart * vvm_heart
    d/dt(memvasns_muscle) <- konps * c_memvas_muscle * cmem * vvm_muscle - koffps * c_memvasns_muscle * vvm_muscle - kintps * c_memvasns_muscle * vvm_muscle
    d/dt(memvasns_skin) <- konps * c_memvas_skin * cmem * vvm_skin - koffps * c_memvasns_skin * vvm_skin - kintps * c_memvasns_skin * vvm_skin
    d/dt(memvasns_adipose) <- konps * c_memvas_adipose * cmem * vvm_adipose - koffps * c_memvasns_adipose * vvm_adipose - kintps * c_memvasns_adipose * vvm_adipose
    d/dt(memvasns_bone) <- konps * c_memvas_bone * cmem * vvm_bone - koffps * c_memvasns_bone * vvm_bone - kintps * c_memvasns_bone * vvm_bone
    d/dt(memvasns_brain) <- konps * c_memvas_brain * cmem * vvm_brain - koffps * c_memvasns_brain * vvm_brain - kintps * c_memvasns_brain * vvm_brain
    d/dt(memvasns_kidney) <- konps * c_memvas_kidney * cmem * vvm_kidney - koffps * c_memvasns_kidney * vvm_kidney - kintps * c_memvasns_kidney * vvm_kidney
    d/dt(memvasns_small_intestine) <- konps * c_memvas_small_intestine * cmem * vvm_small_intestine - koffps * c_memvasns_small_intestine * vvm_small_intestine - kintps * c_memvasns_small_intestine * vvm_small_intestine
    d/dt(memvasns_large_intestine) <- konps * c_memvas_large_intestine * cmem * vvm_large_intestine - koffps * c_memvasns_large_intestine * vvm_large_intestine - kintps * c_memvasns_large_intestine * vvm_large_intestine
    d/dt(memvasns_pancreas) <- konps * c_memvas_pancreas * cmem * vvm_pancreas - koffps * c_memvasns_pancreas * vvm_pancreas - kintps * c_memvasns_pancreas * vvm_pancreas
    d/dt(memvasns_thymus) <- konps * c_memvas_thymus * cmem * vvm_thymus - koffps * c_memvasns_thymus * vvm_thymus - kintps * c_memvasns_thymus * vvm_thymus
    d/dt(memvasns_spleen) <- konps * c_memvas_spleen * cmem * vvm_spleen - koffps * c_memvasns_spleen * vvm_spleen - kintps * c_memvasns_spleen * vvm_spleen
    d/dt(memvasns_other) <- konps * c_memvas_other * cmem * vvm_other - koffps * c_memvasns_other * vvm_other - kintps * c_memvasns_other * vvm_other
    d/dt(memintns_lung) <- konps * c_memint_lung * cmem * vism_lung - koffps * c_memintns_lung * vism_lung - kintps * c_memintns_lung * vism_lung
    d/dt(memintns_liver) <- konps * c_memint_liver * cmem * vism_liver - koffps * c_memintns_liver * vism_liver - kintps * c_memintns_liver * vism_liver
    d/dt(memintns_heart) <- konps * c_memint_heart * cmem * vism_heart - koffps * c_memintns_heart * vism_heart - kintps * c_memintns_heart * vism_heart
    d/dt(memintns_muscle) <- konps * c_memint_muscle * cmem * vism_muscle - koffps * c_memintns_muscle * vism_muscle - kintps * c_memintns_muscle * vism_muscle
    d/dt(memintns_skin) <- konps * c_memint_skin * cmem * vism_skin - koffps * c_memintns_skin * vism_skin - kintps * c_memintns_skin * vism_skin
    d/dt(memintns_adipose) <- konps * c_memint_adipose * cmem * vism_adipose - koffps * c_memintns_adipose * vism_adipose - kintps * c_memintns_adipose * vism_adipose
    d/dt(memintns_bone) <- konps * c_memint_bone * cmem * vism_bone - koffps * c_memintns_bone * vism_bone - kintps * c_memintns_bone * vism_bone
    d/dt(memintns_brain) <- konps * c_memint_brain * cmem * vism_brain - koffps * c_memintns_brain * vism_brain - kintps * c_memintns_brain * vism_brain
    d/dt(memintns_kidney) <- konps * c_memint_kidney * cmem * vism_kidney - koffps * c_memintns_kidney * vism_kidney - kintps * c_memintns_kidney * vism_kidney
    d/dt(memintns_small_intestine) <- konps * c_memint_small_intestine * cmem * vism_small_intestine - koffps * c_memintns_small_intestine * vism_small_intestine - kintps * c_memintns_small_intestine * vism_small_intestine
    d/dt(memintns_large_intestine) <- konps * c_memint_large_intestine * cmem * vism_large_intestine - koffps * c_memintns_large_intestine * vism_large_intestine - kintps * c_memintns_large_intestine * vism_large_intestine
    d/dt(memintns_pancreas) <- konps * c_memint_pancreas * cmem * vism_pancreas - koffps * c_memintns_pancreas * vism_pancreas - kintps * c_memintns_pancreas * vism_pancreas
    d/dt(memintns_thymus) <- konps * c_memint_thymus * cmem * vism_thymus - koffps * c_memintns_thymus * vism_thymus - kintps * c_memintns_thymus * vism_thymus
    d/dt(memintns_spleen) <- konps * c_memint_spleen * cmem * vism_spleen - koffps * c_memintns_spleen * vism_spleen - kintps * c_memintns_spleen * vism_spleen
    d/dt(memintns_other) <- konps * c_memint_other * cmem * vism_other - koffps * c_memintns_other * vism_other - kintps * c_memintns_other * vism_other

    # ---- Endogenous IgG: organ vascular spaces ----
    d/dt(vp_lung_igg) <- (plq_lung + lf_lung) * c_igg_plasma - plq_lung * c_vp_lung_igg - (1 - sigv_lung) * lf_lung * c_vp_lung_igg - clup_lung * fr * c_vp_lung_igg + clup_lung * fr * c_memvas_lung_igg
    d/dt(vp_liver_igg) <- plq_liver * c_vp_lung_igg + (plq_spleen - lf_spleen) * c_vp_spleen_igg + (plq_pancreas - lf_pancreas) * c_vp_pancreas_igg + (plq_small_intestine - lf_small_intestine) * c_vp_small_intestine_igg + (plq_large_intestine - lf_large_intestine) * c_vp_large_intestine_igg - c_vp_liver_igg * (plq_liver - lf_liver + plq_spleen - lf_spleen + plq_pancreas - lf_pancreas + plq_small_intestine - lf_small_intestine + plq_large_intestine - lf_large_intestine) - (1 - sigv_liver) * lf_liver * c_vp_liver_igg - clup_liver * fr * c_vp_liver_igg + clup_liver * fr * c_memvas_liver_igg
    d/dt(vp_heart_igg) <- plq_heart * c_vp_lung_igg - (plq_heart - lf_heart) * c_vp_heart_igg - (1 - sigv_heart) * lf_heart * c_vp_heart_igg - clup_heart * fr * c_vp_heart_igg + clup_heart * fr * c_memvas_heart_igg
    d/dt(vp_muscle_igg) <- plq_muscle * c_vp_lung_igg - (plq_muscle - lf_muscle) * c_vp_muscle_igg - (1 - sigv_muscle) * lf_muscle * c_vp_muscle_igg - clup_muscle * fr * c_vp_muscle_igg + clup_muscle * fr * c_memvas_muscle_igg
    d/dt(vp_skin_igg) <- plq_skin * c_vp_lung_igg - (plq_skin - lf_skin) * c_vp_skin_igg - (1 - sigv_skin) * lf_skin * c_vp_skin_igg - clup_skin * fr * c_vp_skin_igg + clup_skin * fr * c_memvas_skin_igg
    d/dt(vp_adipose_igg) <- plq_adipose * c_vp_lung_igg - (plq_adipose - lf_adipose) * c_vp_adipose_igg - (1 - sigv_adipose) * lf_adipose * c_vp_adipose_igg - clup_adipose * fr * c_vp_adipose_igg + clup_adipose * fr * c_memvas_adipose_igg
    d/dt(vp_bone_igg) <- plq_bone * c_vp_lung_igg - (plq_bone - lf_bone) * c_vp_bone_igg - (1 - sigv_bone) * lf_bone * c_vp_bone_igg - clup_bone * fr * c_vp_bone_igg + clup_bone * fr * c_memvas_bone_igg
    d/dt(vp_brain_igg) <- plq_brain * c_vp_lung_igg - (plq_brain - lf_brain) * c_vp_brain_igg - (1 - sigv_brain) * lf_brain * c_vp_brain_igg - clup_brain * fr * c_vp_brain_igg + clup_brain * fr * c_memvas_brain_igg
    d/dt(vp_kidney_igg) <- plq_kidney * c_vp_lung_igg - (plq_kidney - lf_kidney) * c_vp_kidney_igg - (1 - sigv_kidney) * lf_kidney * c_vp_kidney_igg - clup_kidney * fr * c_vp_kidney_igg + clup_kidney * fr * c_memvas_kidney_igg
    d/dt(vp_small_intestine_igg) <- plq_small_intestine * c_vp_lung_igg - (plq_small_intestine - lf_small_intestine) * c_vp_small_intestine_igg - (1 - sigv_small_intestine) * lf_small_intestine * c_vp_small_intestine_igg - clup_small_intestine * fr * c_vp_small_intestine_igg + clup_small_intestine * fr * c_memvas_small_intestine_igg
    d/dt(vp_large_intestine_igg) <- plq_large_intestine * c_vp_lung_igg - (plq_large_intestine - lf_large_intestine) * c_vp_large_intestine_igg - (1 - sigv_large_intestine) * lf_large_intestine * c_vp_large_intestine_igg - clup_large_intestine * fr * c_vp_large_intestine_igg + clup_large_intestine * fr * c_memvas_large_intestine_igg
    d/dt(vp_pancreas_igg) <- plq_pancreas * c_vp_lung_igg - (plq_pancreas - lf_pancreas) * c_vp_pancreas_igg - (1 - sigv_pancreas) * lf_pancreas * c_vp_pancreas_igg - clup_pancreas * fr * c_vp_pancreas_igg + clup_pancreas * fr * c_memvas_pancreas_igg
    d/dt(vp_thymus_igg) <- plq_thymus * c_vp_lung_igg - (plq_thymus - lf_thymus) * c_vp_thymus_igg - (1 - sigv_thymus) * lf_thymus * c_vp_thymus_igg - clup_thymus * fr * c_vp_thymus_igg + clup_thymus * fr * c_memvas_thymus_igg
    d/dt(vp_spleen_igg) <- plq_spleen * c_vp_lung_igg - (plq_spleen - lf_spleen) * c_vp_spleen_igg - (1 - sigv_spleen) * lf_spleen * c_vp_spleen_igg - clup_spleen * fr * c_vp_spleen_igg + clup_spleen * fr * c_memvas_spleen_igg
    d/dt(vp_other_igg) <- plq_other * c_vp_lung_igg - (plq_other - lf_other) * c_vp_other_igg - (1 - sigv_other) * lf_other * c_vp_other_igg - clup_other * fr * c_vp_other_igg + clup_other * fr * c_memvas_other_igg

    # ---- igg: interstitial spaces ----
    d/dt(is_lung_igg) <- (1 - sigv_lung) * lf_lung * c_vp_lung_igg - clup_lung * (1 - fr) * c_is_lung_igg + clup_lung * (1 - fr) * c_memint_lung_igg - (1 - sigis) * lf_lung * c_is_lung_igg
    d/dt(is_liver_igg) <- (1 - sigv_liver) * lf_liver * c_vp_liver_igg - clup_liver * (1 - fr) * c_is_liver_igg + clup_liver * (1 - fr) * c_memint_liver_igg - (1 - sigis) * lf_liver * c_is_liver_igg
    d/dt(is_heart_igg) <- (1 - sigv_heart) * lf_heart * c_vp_heart_igg - clup_heart * (1 - fr) * c_is_heart_igg + clup_heart * (1 - fr) * c_memint_heart_igg - (1 - sigis) * lf_heart * c_is_heart_igg
    d/dt(is_muscle_igg) <- (1 - sigv_muscle) * lf_muscle * c_vp_muscle_igg - clup_muscle * (1 - fr) * c_is_muscle_igg + clup_muscle * (1 - fr) * c_memint_muscle_igg - (1 - sigis) * lf_muscle * c_is_muscle_igg
    d/dt(is_skin_igg) <- (1 - sigv_skin) * lf_skin * c_vp_skin_igg - clup_skin * (1 - fr) * c_is_skin_igg + clup_skin * (1 - fr) * c_memint_skin_igg - (1 - sigis) * lf_skin * c_is_skin_igg
    d/dt(is_adipose_igg) <- (1 - sigv_adipose) * lf_adipose * c_vp_adipose_igg - clup_adipose * (1 - fr) * c_is_adipose_igg + clup_adipose * (1 - fr) * c_memint_adipose_igg - (1 - sigis) * lf_adipose * c_is_adipose_igg
    d/dt(is_bone_igg) <- (1 - sigv_bone) * lf_bone * c_vp_bone_igg - clup_bone * (1 - fr) * c_is_bone_igg + clup_bone * (1 - fr) * c_memint_bone_igg - (1 - sigis) * lf_bone * c_is_bone_igg
    d/dt(is_brain_igg) <- (1 - sigv_brain) * lf_brain * c_vp_brain_igg - clup_brain * (1 - fr) * c_is_brain_igg + clup_brain * (1 - fr) * c_memint_brain_igg - (1 - sigis) * lf_brain * c_is_brain_igg
    d/dt(is_kidney_igg) <- (1 - sigv_kidney) * lf_kidney * c_vp_kidney_igg - clup_kidney * (1 - fr) * c_is_kidney_igg + clup_kidney * (1 - fr) * c_memint_kidney_igg - (1 - sigis) * lf_kidney * c_is_kidney_igg
    d/dt(is_small_intestine_igg) <- (1 - sigv_small_intestine) * lf_small_intestine * c_vp_small_intestine_igg - clup_small_intestine * (1 - fr) * c_is_small_intestine_igg + clup_small_intestine * (1 - fr) * c_memint_small_intestine_igg - (1 - sigis) * lf_small_intestine * c_is_small_intestine_igg
    d/dt(is_large_intestine_igg) <- (1 - sigv_large_intestine) * lf_large_intestine * c_vp_large_intestine_igg - clup_large_intestine * (1 - fr) * c_is_large_intestine_igg + clup_large_intestine * (1 - fr) * c_memint_large_intestine_igg - (1 - sigis) * lf_large_intestine * c_is_large_intestine_igg
    d/dt(is_pancreas_igg) <- (1 - sigv_pancreas) * lf_pancreas * c_vp_pancreas_igg - clup_pancreas * (1 - fr) * c_is_pancreas_igg + clup_pancreas * (1 - fr) * c_memint_pancreas_igg - (1 - sigis) * lf_pancreas * c_is_pancreas_igg
    d/dt(is_thymus_igg) <- (1 - sigv_thymus) * lf_thymus * c_vp_thymus_igg - clup_thymus * (1 - fr) * c_is_thymus_igg + clup_thymus * (1 - fr) * c_memint_thymus_igg - (1 - sigis) * lf_thymus * c_is_thymus_igg
    d/dt(is_spleen_igg) <- (1 - sigv_spleen) * lf_spleen * c_vp_spleen_igg - clup_spleen * (1 - fr) * c_is_spleen_igg + clup_spleen * (1 - fr) * c_memint_spleen_igg - (1 - sigis) * lf_spleen * c_is_spleen_igg
    d/dt(is_other_igg) <- (1 - sigv_other) * lf_other * c_vp_other_igg - clup_other * (1 - fr) * c_is_other_igg + clup_other * (1 - fr) * c_memint_other_igg - (1 - sigis) * lf_other * c_is_other_igg

    # ---- igg: vascular-side membrane ----
    d/dt(memvas_lung_igg) <- clup_lung * fr * c_vp_lung_igg - clup_lung * fr * c_memvas_lung_igg - clup_lung * fr * c_memvas_lung_igg + clup_lung * fr * c_endorecyc_lung_igg - kon7 * c_memvas_lung_igg * c_fcrnmemvas_lung * vvm_lung + koff7 * c_memvasfr1_lung_igg * vvm_lung
    d/dt(memvas_liver_igg) <- clup_liver * fr * c_vp_liver_igg - clup_liver * fr * c_memvas_liver_igg - clup_liver * fr * c_memvas_liver_igg + clup_liver * fr * c_endorecyc_liver_igg - kon7 * c_memvas_liver_igg * c_fcrnmemvas_liver * vvm_liver + koff7 * c_memvasfr1_liver_igg * vvm_liver
    d/dt(memvas_heart_igg) <- clup_heart * fr * c_vp_heart_igg - clup_heart * fr * c_memvas_heart_igg - clup_heart * fr * c_memvas_heart_igg + clup_heart * fr * c_endorecyc_heart_igg - kon7 * c_memvas_heart_igg * c_fcrnmemvas_heart * vvm_heart + koff7 * c_memvasfr1_heart_igg * vvm_heart
    d/dt(memvas_muscle_igg) <- clup_muscle * fr * c_vp_muscle_igg - clup_muscle * fr * c_memvas_muscle_igg - clup_muscle * fr * c_memvas_muscle_igg + clup_muscle * fr * c_endorecyc_muscle_igg - kon7 * c_memvas_muscle_igg * c_fcrnmemvas_muscle * vvm_muscle + koff7 * c_memvasfr1_muscle_igg * vvm_muscle
    d/dt(memvas_skin_igg) <- clup_skin * fr * c_vp_skin_igg - clup_skin * fr * c_memvas_skin_igg - clup_skin * fr * c_memvas_skin_igg + clup_skin * fr * c_endorecyc_skin_igg - kon7 * c_memvas_skin_igg * c_fcrnmemvas_skin * vvm_skin + koff7 * c_memvasfr1_skin_igg * vvm_skin
    d/dt(memvas_adipose_igg) <- clup_adipose * fr * c_vp_adipose_igg - clup_adipose * fr * c_memvas_adipose_igg - clup_adipose * fr * c_memvas_adipose_igg + clup_adipose * fr * c_endorecyc_adipose_igg - kon7 * c_memvas_adipose_igg * c_fcrnmemvas_adipose * vvm_adipose + koff7 * c_memvasfr1_adipose_igg * vvm_adipose
    d/dt(memvas_bone_igg) <- clup_bone * fr * c_vp_bone_igg - clup_bone * fr * c_memvas_bone_igg - clup_bone * fr * c_memvas_bone_igg + clup_bone * fr * c_endorecyc_bone_igg - kon7 * c_memvas_bone_igg * c_fcrnmemvas_bone * vvm_bone + koff7 * c_memvasfr1_bone_igg * vvm_bone
    d/dt(memvas_brain_igg) <- clup_brain * fr * c_vp_brain_igg - clup_brain * fr * c_memvas_brain_igg - clup_brain * fr * c_memvas_brain_igg + clup_brain * fr * c_endorecyc_brain_igg - kon7 * c_memvas_brain_igg * c_fcrnmemvas_brain * vvm_brain + koff7 * c_memvasfr1_brain_igg * vvm_brain
    d/dt(memvas_kidney_igg) <- clup_kidney * fr * c_vp_kidney_igg - clup_kidney * fr * c_memvas_kidney_igg - clup_kidney * fr * c_memvas_kidney_igg + clup_kidney * fr * c_endorecyc_kidney_igg - kon7 * c_memvas_kidney_igg * c_fcrnmemvas_kidney * vvm_kidney + koff7 * c_memvasfr1_kidney_igg * vvm_kidney
    d/dt(memvas_small_intestine_igg) <- clup_small_intestine * fr * c_vp_small_intestine_igg - clup_small_intestine * fr * c_memvas_small_intestine_igg - clup_small_intestine * fr * c_memvas_small_intestine_igg + clup_small_intestine * fr * c_endorecyc_small_intestine_igg - kon7 * c_memvas_small_intestine_igg * c_fcrnmemvas_small_intestine * vvm_small_intestine + koff7 * c_memvasfr1_small_intestine_igg * vvm_small_intestine
    d/dt(memvas_large_intestine_igg) <- clup_large_intestine * fr * c_vp_large_intestine_igg - clup_large_intestine * fr * c_memvas_large_intestine_igg - clup_large_intestine * fr * c_memvas_large_intestine_igg + clup_large_intestine * fr * c_endorecyc_large_intestine_igg - kon7 * c_memvas_large_intestine_igg * c_fcrnmemvas_large_intestine * vvm_large_intestine + koff7 * c_memvasfr1_large_intestine_igg * vvm_large_intestine
    d/dt(memvas_pancreas_igg) <- clup_pancreas * fr * c_vp_pancreas_igg - clup_pancreas * fr * c_memvas_pancreas_igg - clup_pancreas * fr * c_memvas_pancreas_igg + clup_pancreas * fr * c_endorecyc_pancreas_igg - kon7 * c_memvas_pancreas_igg * c_fcrnmemvas_pancreas * vvm_pancreas + koff7 * c_memvasfr1_pancreas_igg * vvm_pancreas
    d/dt(memvas_thymus_igg) <- clup_thymus * fr * c_vp_thymus_igg - clup_thymus * fr * c_memvas_thymus_igg - clup_thymus * fr * c_memvas_thymus_igg + clup_thymus * fr * c_endorecyc_thymus_igg - kon7 * c_memvas_thymus_igg * c_fcrnmemvas_thymus * vvm_thymus + koff7 * c_memvasfr1_thymus_igg * vvm_thymus
    d/dt(memvas_spleen_igg) <- clup_spleen * fr * c_vp_spleen_igg - clup_spleen * fr * c_memvas_spleen_igg - clup_spleen * fr * c_memvas_spleen_igg + clup_spleen * fr * c_endorecyc_spleen_igg - kon7 * c_memvas_spleen_igg * c_fcrnmemvas_spleen * vvm_spleen + koff7 * c_memvasfr1_spleen_igg * vvm_spleen
    d/dt(memvas_other_igg) <- clup_other * fr * c_vp_other_igg - clup_other * fr * c_memvas_other_igg - clup_other * fr * c_memvas_other_igg + clup_other * fr * c_endorecyc_other_igg - kon7 * c_memvas_other_igg * c_fcrnmemvas_other * vvm_other + koff7 * c_memvasfr1_other_igg * vvm_other

    # ---- igg: early endosome (pH 7.4) ----
    d/dt(endoearly_lung_igg) <- clup_lung * fr * c_memvas_lung_igg + clup_lung * (1 - fr) * c_memint_lung_igg - clup_lung * c_endoearly_lung_igg - kon7 * c_endoearly_lung_igg * c_fcrnendoearly_lung * ve7_lung + koff7 * c_endoearlyfr1_lung_igg * ve7_lung
    d/dt(endoearly_liver_igg) <- clup_liver * fr * c_memvas_liver_igg + clup_liver * (1 - fr) * c_memint_liver_igg - clup_liver * c_endoearly_liver_igg - kon7 * c_endoearly_liver_igg * c_fcrnendoearly_liver * ve7_liver + koff7 * c_endoearlyfr1_liver_igg * ve7_liver
    d/dt(endoearly_heart_igg) <- clup_heart * fr * c_memvas_heart_igg + clup_heart * (1 - fr) * c_memint_heart_igg - clup_heart * c_endoearly_heart_igg - kon7 * c_endoearly_heart_igg * c_fcrnendoearly_heart * ve7_heart + koff7 * c_endoearlyfr1_heart_igg * ve7_heart
    d/dt(endoearly_muscle_igg) <- clup_muscle * fr * c_memvas_muscle_igg + clup_muscle * (1 - fr) * c_memint_muscle_igg - clup_muscle * c_endoearly_muscle_igg - kon7 * c_endoearly_muscle_igg * c_fcrnendoearly_muscle * ve7_muscle + koff7 * c_endoearlyfr1_muscle_igg * ve7_muscle
    d/dt(endoearly_skin_igg) <- clup_skin * fr * c_memvas_skin_igg + clup_skin * (1 - fr) * c_memint_skin_igg - clup_skin * c_endoearly_skin_igg - kon7 * c_endoearly_skin_igg * c_fcrnendoearly_skin * ve7_skin + koff7 * c_endoearlyfr1_skin_igg * ve7_skin
    d/dt(endoearly_adipose_igg) <- clup_adipose * fr * c_memvas_adipose_igg + clup_adipose * (1 - fr) * c_memint_adipose_igg - clup_adipose * c_endoearly_adipose_igg - kon7 * c_endoearly_adipose_igg * c_fcrnendoearly_adipose * ve7_adipose + koff7 * c_endoearlyfr1_adipose_igg * ve7_adipose
    d/dt(endoearly_bone_igg) <- clup_bone * fr * c_memvas_bone_igg + clup_bone * (1 - fr) * c_memint_bone_igg - clup_bone * c_endoearly_bone_igg - kon7 * c_endoearly_bone_igg * c_fcrnendoearly_bone * ve7_bone + koff7 * c_endoearlyfr1_bone_igg * ve7_bone
    d/dt(endoearly_brain_igg) <- clup_brain * fr * c_memvas_brain_igg + clup_brain * (1 - fr) * c_memint_brain_igg - clup_brain * c_endoearly_brain_igg - kon7 * c_endoearly_brain_igg * c_fcrnendoearly_brain * ve7_brain + koff7 * c_endoearlyfr1_brain_igg * ve7_brain
    d/dt(endoearly_kidney_igg) <- clup_kidney * fr * c_memvas_kidney_igg + clup_kidney * (1 - fr) * c_memint_kidney_igg - clup_kidney * c_endoearly_kidney_igg - kon7 * c_endoearly_kidney_igg * c_fcrnendoearly_kidney * ve7_kidney + koff7 * c_endoearlyfr1_kidney_igg * ve7_kidney
    d/dt(endoearly_small_intestine_igg) <- clup_small_intestine * fr * c_memvas_small_intestine_igg + clup_small_intestine * (1 - fr) * c_memint_small_intestine_igg - clup_small_intestine * c_endoearly_small_intestine_igg - kon7 * c_endoearly_small_intestine_igg * c_fcrnendoearly_small_intestine * ve7_small_intestine + koff7 * c_endoearlyfr1_small_intestine_igg * ve7_small_intestine
    d/dt(endoearly_large_intestine_igg) <- clup_large_intestine * fr * c_memvas_large_intestine_igg + clup_large_intestine * (1 - fr) * c_memint_large_intestine_igg - clup_large_intestine * c_endoearly_large_intestine_igg - kon7 * c_endoearly_large_intestine_igg * c_fcrnendoearly_large_intestine * ve7_large_intestine + koff7 * c_endoearlyfr1_large_intestine_igg * ve7_large_intestine
    d/dt(endoearly_pancreas_igg) <- clup_pancreas * fr * c_memvas_pancreas_igg + clup_pancreas * (1 - fr) * c_memint_pancreas_igg - clup_pancreas * c_endoearly_pancreas_igg - kon7 * c_endoearly_pancreas_igg * c_fcrnendoearly_pancreas * ve7_pancreas + koff7 * c_endoearlyfr1_pancreas_igg * ve7_pancreas
    d/dt(endoearly_thymus_igg) <- clup_thymus * fr * c_memvas_thymus_igg + clup_thymus * (1 - fr) * c_memint_thymus_igg - clup_thymus * c_endoearly_thymus_igg - kon7 * c_endoearly_thymus_igg * c_fcrnendoearly_thymus * ve7_thymus + koff7 * c_endoearlyfr1_thymus_igg * ve7_thymus
    d/dt(endoearly_spleen_igg) <- clup_spleen * fr * c_memvas_spleen_igg + clup_spleen * (1 - fr) * c_memint_spleen_igg - clup_spleen * c_endoearly_spleen_igg - kon7 * c_endoearly_spleen_igg * c_fcrnendoearly_spleen * ve7_spleen + koff7 * c_endoearlyfr1_spleen_igg * ve7_spleen
    d/dt(endoearly_other_igg) <- clup_other * fr * c_memvas_other_igg + clup_other * (1 - fr) * c_memint_other_igg - clup_other * c_endoearly_other_igg - kon7 * c_endoearly_other_igg * c_fcrnendoearly_other * ve7_other + koff7 * c_endoearlyfr1_other_igg * ve7_other

    # ---- igg: sorting endosome (pH 6.0) ----
    d/dt(endosort_lung_igg) <- clup_lung * c_endoearly_lung_igg - clup_lung * c_endosort_lung_igg - kon6 * c_endosort_lung_igg * c_fcrnendosort_lung * ve6a_lung + koff6 * c_endosortfr1_lung_igg * ve6a_lung
    d/dt(endosort_liver_igg) <- clup_liver * c_endoearly_liver_igg - clup_liver * c_endosort_liver_igg - kon6 * c_endosort_liver_igg * c_fcrnendosort_liver * ve6a_liver + koff6 * c_endosortfr1_liver_igg * ve6a_liver
    d/dt(endosort_heart_igg) <- clup_heart * c_endoearly_heart_igg - clup_heart * c_endosort_heart_igg - kon6 * c_endosort_heart_igg * c_fcrnendosort_heart * ve6a_heart + koff6 * c_endosortfr1_heart_igg * ve6a_heart
    d/dt(endosort_muscle_igg) <- clup_muscle * c_endoearly_muscle_igg - clup_muscle * c_endosort_muscle_igg - kon6 * c_endosort_muscle_igg * c_fcrnendosort_muscle * ve6a_muscle + koff6 * c_endosortfr1_muscle_igg * ve6a_muscle
    d/dt(endosort_skin_igg) <- clup_skin * c_endoearly_skin_igg - clup_skin * c_endosort_skin_igg - kon6 * c_endosort_skin_igg * c_fcrnendosort_skin * ve6a_skin + koff6 * c_endosortfr1_skin_igg * ve6a_skin
    d/dt(endosort_adipose_igg) <- clup_adipose * c_endoearly_adipose_igg - clup_adipose * c_endosort_adipose_igg - kon6 * c_endosort_adipose_igg * c_fcrnendosort_adipose * ve6a_adipose + koff6 * c_endosortfr1_adipose_igg * ve6a_adipose
    d/dt(endosort_bone_igg) <- clup_bone * c_endoearly_bone_igg - clup_bone * c_endosort_bone_igg - kon6 * c_endosort_bone_igg * c_fcrnendosort_bone * ve6a_bone + koff6 * c_endosortfr1_bone_igg * ve6a_bone
    d/dt(endosort_brain_igg) <- clup_brain * c_endoearly_brain_igg - clup_brain * c_endosort_brain_igg - kon6 * c_endosort_brain_igg * c_fcrnendosort_brain * ve6a_brain + koff6 * c_endosortfr1_brain_igg * ve6a_brain
    d/dt(endosort_kidney_igg) <- clup_kidney * c_endoearly_kidney_igg - clup_kidney * c_endosort_kidney_igg - kon6 * c_endosort_kidney_igg * c_fcrnendosort_kidney * ve6a_kidney + koff6 * c_endosortfr1_kidney_igg * ve6a_kidney
    d/dt(endosort_small_intestine_igg) <- clup_small_intestine * c_endoearly_small_intestine_igg - clup_small_intestine * c_endosort_small_intestine_igg - kon6 * c_endosort_small_intestine_igg * c_fcrnendosort_small_intestine * ve6a_small_intestine + koff6 * c_endosortfr1_small_intestine_igg * ve6a_small_intestine
    d/dt(endosort_large_intestine_igg) <- clup_large_intestine * c_endoearly_large_intestine_igg - clup_large_intestine * c_endosort_large_intestine_igg - kon6 * c_endosort_large_intestine_igg * c_fcrnendosort_large_intestine * ve6a_large_intestine + koff6 * c_endosortfr1_large_intestine_igg * ve6a_large_intestine
    d/dt(endosort_pancreas_igg) <- clup_pancreas * c_endoearly_pancreas_igg - clup_pancreas * c_endosort_pancreas_igg - kon6 * c_endosort_pancreas_igg * c_fcrnendosort_pancreas * ve6a_pancreas + koff6 * c_endosortfr1_pancreas_igg * ve6a_pancreas
    d/dt(endosort_thymus_igg) <- clup_thymus * c_endoearly_thymus_igg - clup_thymus * c_endosort_thymus_igg - kon6 * c_endosort_thymus_igg * c_fcrnendosort_thymus * ve6a_thymus + koff6 * c_endosortfr1_thymus_igg * ve6a_thymus
    d/dt(endosort_spleen_igg) <- clup_spleen * c_endoearly_spleen_igg - clup_spleen * c_endosort_spleen_igg - kon6 * c_endosort_spleen_igg * c_fcrnendosort_spleen * ve6a_spleen + koff6 * c_endosortfr1_spleen_igg * ve6a_spleen
    d/dt(endosort_other_igg) <- clup_other * c_endoearly_other_igg - clup_other * c_endosort_other_igg - kon6 * c_endosort_other_igg * c_fcrnendosort_other * ve6a_other + koff6 * c_endosortfr1_other_igg * ve6a_other

    # ---- igg: recycling endosome (pH 7.4) ----
    # Only the (1 - Prob_deg) fraction of unbound mAb escapes lysosomal routing.
    d/dt(endorecyc_lung_igg) <- clup_lung * (1 - probdeg) * c_endosort_lung_igg - clup_lung * c_endorecyc_lung_igg - kon7 * c_endorecyc_lung_igg * c_fcrnendorecyc_lung * ve7b_lung + koff7 * c_endorecycfr1_lung_igg * ve7b_lung
    d/dt(endorecyc_liver_igg) <- clup_liver * (1 - probdeg) * c_endosort_liver_igg - clup_liver * c_endorecyc_liver_igg - kon7 * c_endorecyc_liver_igg * c_fcrnendorecyc_liver * ve7b_liver + koff7 * c_endorecycfr1_liver_igg * ve7b_liver
    d/dt(endorecyc_heart_igg) <- clup_heart * (1 - probdeg) * c_endosort_heart_igg - clup_heart * c_endorecyc_heart_igg - kon7 * c_endorecyc_heart_igg * c_fcrnendorecyc_heart * ve7b_heart + koff7 * c_endorecycfr1_heart_igg * ve7b_heart
    d/dt(endorecyc_muscle_igg) <- clup_muscle * (1 - probdeg) * c_endosort_muscle_igg - clup_muscle * c_endorecyc_muscle_igg - kon7 * c_endorecyc_muscle_igg * c_fcrnendorecyc_muscle * ve7b_muscle + koff7 * c_endorecycfr1_muscle_igg * ve7b_muscle
    d/dt(endorecyc_skin_igg) <- clup_skin * (1 - probdeg) * c_endosort_skin_igg - clup_skin * c_endorecyc_skin_igg - kon7 * c_endorecyc_skin_igg * c_fcrnendorecyc_skin * ve7b_skin + koff7 * c_endorecycfr1_skin_igg * ve7b_skin
    d/dt(endorecyc_adipose_igg) <- clup_adipose * (1 - probdeg) * c_endosort_adipose_igg - clup_adipose * c_endorecyc_adipose_igg - kon7 * c_endorecyc_adipose_igg * c_fcrnendorecyc_adipose * ve7b_adipose + koff7 * c_endorecycfr1_adipose_igg * ve7b_adipose
    d/dt(endorecyc_bone_igg) <- clup_bone * (1 - probdeg) * c_endosort_bone_igg - clup_bone * c_endorecyc_bone_igg - kon7 * c_endorecyc_bone_igg * c_fcrnendorecyc_bone * ve7b_bone + koff7 * c_endorecycfr1_bone_igg * ve7b_bone
    d/dt(endorecyc_brain_igg) <- clup_brain * (1 - probdeg) * c_endosort_brain_igg - clup_brain * c_endorecyc_brain_igg - kon7 * c_endorecyc_brain_igg * c_fcrnendorecyc_brain * ve7b_brain + koff7 * c_endorecycfr1_brain_igg * ve7b_brain
    d/dt(endorecyc_kidney_igg) <- clup_kidney * (1 - probdeg) * c_endosort_kidney_igg - clup_kidney * c_endorecyc_kidney_igg - kon7 * c_endorecyc_kidney_igg * c_fcrnendorecyc_kidney * ve7b_kidney + koff7 * c_endorecycfr1_kidney_igg * ve7b_kidney
    d/dt(endorecyc_small_intestine_igg) <- clup_small_intestine * (1 - probdeg) * c_endosort_small_intestine_igg - clup_small_intestine * c_endorecyc_small_intestine_igg - kon7 * c_endorecyc_small_intestine_igg * c_fcrnendorecyc_small_intestine * ve7b_small_intestine + koff7 * c_endorecycfr1_small_intestine_igg * ve7b_small_intestine
    d/dt(endorecyc_large_intestine_igg) <- clup_large_intestine * (1 - probdeg) * c_endosort_large_intestine_igg - clup_large_intestine * c_endorecyc_large_intestine_igg - kon7 * c_endorecyc_large_intestine_igg * c_fcrnendorecyc_large_intestine * ve7b_large_intestine + koff7 * c_endorecycfr1_large_intestine_igg * ve7b_large_intestine
    d/dt(endorecyc_pancreas_igg) <- clup_pancreas * (1 - probdeg) * c_endosort_pancreas_igg - clup_pancreas * c_endorecyc_pancreas_igg - kon7 * c_endorecyc_pancreas_igg * c_fcrnendorecyc_pancreas * ve7b_pancreas + koff7 * c_endorecycfr1_pancreas_igg * ve7b_pancreas
    d/dt(endorecyc_thymus_igg) <- clup_thymus * (1 - probdeg) * c_endosort_thymus_igg - clup_thymus * c_endorecyc_thymus_igg - kon7 * c_endorecyc_thymus_igg * c_fcrnendorecyc_thymus * ve7b_thymus + koff7 * c_endorecycfr1_thymus_igg * ve7b_thymus
    d/dt(endorecyc_spleen_igg) <- clup_spleen * (1 - probdeg) * c_endosort_spleen_igg - clup_spleen * c_endorecyc_spleen_igg - kon7 * c_endorecyc_spleen_igg * c_fcrnendorecyc_spleen * ve7b_spleen + koff7 * c_endorecycfr1_spleen_igg * ve7b_spleen
    d/dt(endorecyc_other_igg) <- clup_other * (1 - probdeg) * c_endosort_other_igg - clup_other * c_endorecyc_other_igg - kon7 * c_endorecyc_other_igg * c_fcrnendorecyc_other * ve7b_other + koff7 * c_endorecycfr1_other_igg * ve7b_other

    # ---- igg: interstitial-side membrane ----
    d/dt(memint_lung_igg) <- clup_lung * (1 - fr) * c_is_lung_igg - clup_lung * (1 - fr) * c_memint_lung_igg + clup_lung * (1 - fr) * c_endorecyc_lung_igg - clup_lung * (1 - fr) * c_memint_lung_igg - kon7 * c_memint_lung_igg * c_fcrnmemint_lung * vism_lung + koff7 * c_memintfr1_lung_igg * vism_lung
    d/dt(memint_liver_igg) <- clup_liver * (1 - fr) * c_is_liver_igg - clup_liver * (1 - fr) * c_memint_liver_igg + clup_liver * (1 - fr) * c_endorecyc_liver_igg - clup_liver * (1 - fr) * c_memint_liver_igg - kon7 * c_memint_liver_igg * c_fcrnmemint_liver * vism_liver + koff7 * c_memintfr1_liver_igg * vism_liver
    d/dt(memint_heart_igg) <- clup_heart * (1 - fr) * c_is_heart_igg - clup_heart * (1 - fr) * c_memint_heart_igg + clup_heart * (1 - fr) * c_endorecyc_heart_igg - clup_heart * (1 - fr) * c_memint_heart_igg - kon7 * c_memint_heart_igg * c_fcrnmemint_heart * vism_heart + koff7 * c_memintfr1_heart_igg * vism_heart
    d/dt(memint_muscle_igg) <- clup_muscle * (1 - fr) * c_is_muscle_igg - clup_muscle * (1 - fr) * c_memint_muscle_igg + clup_muscle * (1 - fr) * c_endorecyc_muscle_igg - clup_muscle * (1 - fr) * c_memint_muscle_igg - kon7 * c_memint_muscle_igg * c_fcrnmemint_muscle * vism_muscle + koff7 * c_memintfr1_muscle_igg * vism_muscle
    d/dt(memint_skin_igg) <- clup_skin * (1 - fr) * c_is_skin_igg - clup_skin * (1 - fr) * c_memint_skin_igg + clup_skin * (1 - fr) * c_endorecyc_skin_igg - clup_skin * (1 - fr) * c_memint_skin_igg - kon7 * c_memint_skin_igg * c_fcrnmemint_skin * vism_skin + koff7 * c_memintfr1_skin_igg * vism_skin
    d/dt(memint_adipose_igg) <- clup_adipose * (1 - fr) * c_is_adipose_igg - clup_adipose * (1 - fr) * c_memint_adipose_igg + clup_adipose * (1 - fr) * c_endorecyc_adipose_igg - clup_adipose * (1 - fr) * c_memint_adipose_igg - kon7 * c_memint_adipose_igg * c_fcrnmemint_adipose * vism_adipose + koff7 * c_memintfr1_adipose_igg * vism_adipose
    d/dt(memint_bone_igg) <- clup_bone * (1 - fr) * c_is_bone_igg - clup_bone * (1 - fr) * c_memint_bone_igg + clup_bone * (1 - fr) * c_endorecyc_bone_igg - clup_bone * (1 - fr) * c_memint_bone_igg - kon7 * c_memint_bone_igg * c_fcrnmemint_bone * vism_bone + koff7 * c_memintfr1_bone_igg * vism_bone
    d/dt(memint_brain_igg) <- clup_brain * (1 - fr) * c_is_brain_igg - clup_brain * (1 - fr) * c_memint_brain_igg + clup_brain * (1 - fr) * c_endorecyc_brain_igg - clup_brain * (1 - fr) * c_memint_brain_igg - kon7 * c_memint_brain_igg * c_fcrnmemint_brain * vism_brain + koff7 * c_memintfr1_brain_igg * vism_brain
    d/dt(memint_kidney_igg) <- clup_kidney * (1 - fr) * c_is_kidney_igg - clup_kidney * (1 - fr) * c_memint_kidney_igg + clup_kidney * (1 - fr) * c_endorecyc_kidney_igg - clup_kidney * (1 - fr) * c_memint_kidney_igg - kon7 * c_memint_kidney_igg * c_fcrnmemint_kidney * vism_kidney + koff7 * c_memintfr1_kidney_igg * vism_kidney
    d/dt(memint_small_intestine_igg) <- clup_small_intestine * (1 - fr) * c_is_small_intestine_igg - clup_small_intestine * (1 - fr) * c_memint_small_intestine_igg + clup_small_intestine * (1 - fr) * c_endorecyc_small_intestine_igg - clup_small_intestine * (1 - fr) * c_memint_small_intestine_igg - kon7 * c_memint_small_intestine_igg * c_fcrnmemint_small_intestine * vism_small_intestine + koff7 * c_memintfr1_small_intestine_igg * vism_small_intestine
    d/dt(memint_large_intestine_igg) <- clup_large_intestine * (1 - fr) * c_is_large_intestine_igg - clup_large_intestine * (1 - fr) * c_memint_large_intestine_igg + clup_large_intestine * (1 - fr) * c_endorecyc_large_intestine_igg - clup_large_intestine * (1 - fr) * c_memint_large_intestine_igg - kon7 * c_memint_large_intestine_igg * c_fcrnmemint_large_intestine * vism_large_intestine + koff7 * c_memintfr1_large_intestine_igg * vism_large_intestine
    d/dt(memint_pancreas_igg) <- clup_pancreas * (1 - fr) * c_is_pancreas_igg - clup_pancreas * (1 - fr) * c_memint_pancreas_igg + clup_pancreas * (1 - fr) * c_endorecyc_pancreas_igg - clup_pancreas * (1 - fr) * c_memint_pancreas_igg - kon7 * c_memint_pancreas_igg * c_fcrnmemint_pancreas * vism_pancreas + koff7 * c_memintfr1_pancreas_igg * vism_pancreas
    d/dt(memint_thymus_igg) <- clup_thymus * (1 - fr) * c_is_thymus_igg - clup_thymus * (1 - fr) * c_memint_thymus_igg + clup_thymus * (1 - fr) * c_endorecyc_thymus_igg - clup_thymus * (1 - fr) * c_memint_thymus_igg - kon7 * c_memint_thymus_igg * c_fcrnmemint_thymus * vism_thymus + koff7 * c_memintfr1_thymus_igg * vism_thymus
    d/dt(memint_spleen_igg) <- clup_spleen * (1 - fr) * c_is_spleen_igg - clup_spleen * (1 - fr) * c_memint_spleen_igg + clup_spleen * (1 - fr) * c_endorecyc_spleen_igg - clup_spleen * (1 - fr) * c_memint_spleen_igg - kon7 * c_memint_spleen_igg * c_fcrnmemint_spleen * vism_spleen + koff7 * c_memintfr1_spleen_igg * vism_spleen
    d/dt(memint_other_igg) <- clup_other * (1 - fr) * c_is_other_igg - clup_other * (1 - fr) * c_memint_other_igg + clup_other * (1 - fr) * c_endorecyc_other_igg - clup_other * (1 - fr) * c_memint_other_igg - kon7 * c_memint_other_igg * c_fcrnmemint_other * vism_other + koff7 * c_memintfr1_other_igg * vism_other

    # ---- igg: 1:1 FcRn complexes ----
    d/dt(memvasfr1_lung_igg) <- kon7 * c_memvas_lung_igg * c_fcrnmemvas_lung * vvm_lung - koff7 * c_memvasfr1_lung_igg * vvm_lung - kon7b * c_memvasfr1_lung_igg * c_fcrnmemvas_lung * vvm_lung + koff7b * c_memvasfr2_lung_igg * vvm_lung + clup_lung * fr * c_endorecycfr1_lung_igg - clup_lung * fr * c_memvasfr1_lung_igg - kdegab * c_memvasfr1_lung_igg * vvm_lung
    d/dt(memvasfr1_liver_igg) <- kon7 * c_memvas_liver_igg * c_fcrnmemvas_liver * vvm_liver - koff7 * c_memvasfr1_liver_igg * vvm_liver - kon7b * c_memvasfr1_liver_igg * c_fcrnmemvas_liver * vvm_liver + koff7b * c_memvasfr2_liver_igg * vvm_liver + clup_liver * fr * c_endorecycfr1_liver_igg - clup_liver * fr * c_memvasfr1_liver_igg - kdegab * c_memvasfr1_liver_igg * vvm_liver
    d/dt(memvasfr1_heart_igg) <- kon7 * c_memvas_heart_igg * c_fcrnmemvas_heart * vvm_heart - koff7 * c_memvasfr1_heart_igg * vvm_heart - kon7b * c_memvasfr1_heart_igg * c_fcrnmemvas_heart * vvm_heart + koff7b * c_memvasfr2_heart_igg * vvm_heart + clup_heart * fr * c_endorecycfr1_heart_igg - clup_heart * fr * c_memvasfr1_heart_igg - kdegab * c_memvasfr1_heart_igg * vvm_heart
    d/dt(memvasfr1_muscle_igg) <- kon7 * c_memvas_muscle_igg * c_fcrnmemvas_muscle * vvm_muscle - koff7 * c_memvasfr1_muscle_igg * vvm_muscle - kon7b * c_memvasfr1_muscle_igg * c_fcrnmemvas_muscle * vvm_muscle + koff7b * c_memvasfr2_muscle_igg * vvm_muscle + clup_muscle * fr * c_endorecycfr1_muscle_igg - clup_muscle * fr * c_memvasfr1_muscle_igg - kdegab * c_memvasfr1_muscle_igg * vvm_muscle
    d/dt(memvasfr1_skin_igg) <- kon7 * c_memvas_skin_igg * c_fcrnmemvas_skin * vvm_skin - koff7 * c_memvasfr1_skin_igg * vvm_skin - kon7b * c_memvasfr1_skin_igg * c_fcrnmemvas_skin * vvm_skin + koff7b * c_memvasfr2_skin_igg * vvm_skin + clup_skin * fr * c_endorecycfr1_skin_igg - clup_skin * fr * c_memvasfr1_skin_igg - kdegab * c_memvasfr1_skin_igg * vvm_skin
    d/dt(memvasfr1_adipose_igg) <- kon7 * c_memvas_adipose_igg * c_fcrnmemvas_adipose * vvm_adipose - koff7 * c_memvasfr1_adipose_igg * vvm_adipose - kon7b * c_memvasfr1_adipose_igg * c_fcrnmemvas_adipose * vvm_adipose + koff7b * c_memvasfr2_adipose_igg * vvm_adipose + clup_adipose * fr * c_endorecycfr1_adipose_igg - clup_adipose * fr * c_memvasfr1_adipose_igg - kdegab * c_memvasfr1_adipose_igg * vvm_adipose
    d/dt(memvasfr1_bone_igg) <- kon7 * c_memvas_bone_igg * c_fcrnmemvas_bone * vvm_bone - koff7 * c_memvasfr1_bone_igg * vvm_bone - kon7b * c_memvasfr1_bone_igg * c_fcrnmemvas_bone * vvm_bone + koff7b * c_memvasfr2_bone_igg * vvm_bone + clup_bone * fr * c_endorecycfr1_bone_igg - clup_bone * fr * c_memvasfr1_bone_igg - kdegab * c_memvasfr1_bone_igg * vvm_bone
    d/dt(memvasfr1_brain_igg) <- kon7 * c_memvas_brain_igg * c_fcrnmemvas_brain * vvm_brain - koff7 * c_memvasfr1_brain_igg * vvm_brain - kon7b * c_memvasfr1_brain_igg * c_fcrnmemvas_brain * vvm_brain + koff7b * c_memvasfr2_brain_igg * vvm_brain + clup_brain * fr * c_endorecycfr1_brain_igg - clup_brain * fr * c_memvasfr1_brain_igg - kdegab * c_memvasfr1_brain_igg * vvm_brain
    d/dt(memvasfr1_kidney_igg) <- kon7 * c_memvas_kidney_igg * c_fcrnmemvas_kidney * vvm_kidney - koff7 * c_memvasfr1_kidney_igg * vvm_kidney - kon7b * c_memvasfr1_kidney_igg * c_fcrnmemvas_kidney * vvm_kidney + koff7b * c_memvasfr2_kidney_igg * vvm_kidney + clup_kidney * fr * c_endorecycfr1_kidney_igg - clup_kidney * fr * c_memvasfr1_kidney_igg - kdegab * c_memvasfr1_kidney_igg * vvm_kidney
    d/dt(memvasfr1_small_intestine_igg) <- kon7 * c_memvas_small_intestine_igg * c_fcrnmemvas_small_intestine * vvm_small_intestine - koff7 * c_memvasfr1_small_intestine_igg * vvm_small_intestine - kon7b * c_memvasfr1_small_intestine_igg * c_fcrnmemvas_small_intestine * vvm_small_intestine + koff7b * c_memvasfr2_small_intestine_igg * vvm_small_intestine + clup_small_intestine * fr * c_endorecycfr1_small_intestine_igg - clup_small_intestine * fr * c_memvasfr1_small_intestine_igg - kdegab * c_memvasfr1_small_intestine_igg * vvm_small_intestine
    d/dt(memvasfr1_large_intestine_igg) <- kon7 * c_memvas_large_intestine_igg * c_fcrnmemvas_large_intestine * vvm_large_intestine - koff7 * c_memvasfr1_large_intestine_igg * vvm_large_intestine - kon7b * c_memvasfr1_large_intestine_igg * c_fcrnmemvas_large_intestine * vvm_large_intestine + koff7b * c_memvasfr2_large_intestine_igg * vvm_large_intestine + clup_large_intestine * fr * c_endorecycfr1_large_intestine_igg - clup_large_intestine * fr * c_memvasfr1_large_intestine_igg - kdegab * c_memvasfr1_large_intestine_igg * vvm_large_intestine
    d/dt(memvasfr1_pancreas_igg) <- kon7 * c_memvas_pancreas_igg * c_fcrnmemvas_pancreas * vvm_pancreas - koff7 * c_memvasfr1_pancreas_igg * vvm_pancreas - kon7b * c_memvasfr1_pancreas_igg * c_fcrnmemvas_pancreas * vvm_pancreas + koff7b * c_memvasfr2_pancreas_igg * vvm_pancreas + clup_pancreas * fr * c_endorecycfr1_pancreas_igg - clup_pancreas * fr * c_memvasfr1_pancreas_igg - kdegab * c_memvasfr1_pancreas_igg * vvm_pancreas
    d/dt(memvasfr1_thymus_igg) <- kon7 * c_memvas_thymus_igg * c_fcrnmemvas_thymus * vvm_thymus - koff7 * c_memvasfr1_thymus_igg * vvm_thymus - kon7b * c_memvasfr1_thymus_igg * c_fcrnmemvas_thymus * vvm_thymus + koff7b * c_memvasfr2_thymus_igg * vvm_thymus + clup_thymus * fr * c_endorecycfr1_thymus_igg - clup_thymus * fr * c_memvasfr1_thymus_igg - kdegab * c_memvasfr1_thymus_igg * vvm_thymus
    d/dt(memvasfr1_spleen_igg) <- kon7 * c_memvas_spleen_igg * c_fcrnmemvas_spleen * vvm_spleen - koff7 * c_memvasfr1_spleen_igg * vvm_spleen - kon7b * c_memvasfr1_spleen_igg * c_fcrnmemvas_spleen * vvm_spleen + koff7b * c_memvasfr2_spleen_igg * vvm_spleen + clup_spleen * fr * c_endorecycfr1_spleen_igg - clup_spleen * fr * c_memvasfr1_spleen_igg - kdegab * c_memvasfr1_spleen_igg * vvm_spleen
    d/dt(memvasfr1_other_igg) <- kon7 * c_memvas_other_igg * c_fcrnmemvas_other * vvm_other - koff7 * c_memvasfr1_other_igg * vvm_other - kon7b * c_memvasfr1_other_igg * c_fcrnmemvas_other * vvm_other + koff7b * c_memvasfr2_other_igg * vvm_other + clup_other * fr * c_endorecycfr1_other_igg - clup_other * fr * c_memvasfr1_other_igg - kdegab * c_memvasfr1_other_igg * vvm_other
    d/dt(endoearlyfr1_lung_igg) <- kon7 * c_endoearly_lung_igg * c_fcrnendoearly_lung * ve7_lung - koff7 * c_endoearlyfr1_lung_igg * ve7_lung - kon7b * c_endoearlyfr1_lung_igg * c_fcrnendoearly_lung * ve7_lung + koff7b * c_endoearlyfr2_lung_igg * ve7_lung - clup_lung * c_endoearlyfr1_lung_igg + clup_lung * fr * c_memvasfr1_lung_igg + clup_lung * (1 - fr) * c_memintfr1_lung_igg
    d/dt(endoearlyfr1_liver_igg) <- kon7 * c_endoearly_liver_igg * c_fcrnendoearly_liver * ve7_liver - koff7 * c_endoearlyfr1_liver_igg * ve7_liver - kon7b * c_endoearlyfr1_liver_igg * c_fcrnendoearly_liver * ve7_liver + koff7b * c_endoearlyfr2_liver_igg * ve7_liver - clup_liver * c_endoearlyfr1_liver_igg + clup_liver * fr * c_memvasfr1_liver_igg + clup_liver * (1 - fr) * c_memintfr1_liver_igg
    d/dt(endoearlyfr1_heart_igg) <- kon7 * c_endoearly_heart_igg * c_fcrnendoearly_heart * ve7_heart - koff7 * c_endoearlyfr1_heart_igg * ve7_heart - kon7b * c_endoearlyfr1_heart_igg * c_fcrnendoearly_heart * ve7_heart + koff7b * c_endoearlyfr2_heart_igg * ve7_heart - clup_heart * c_endoearlyfr1_heart_igg + clup_heart * fr * c_memvasfr1_heart_igg + clup_heart * (1 - fr) * c_memintfr1_heart_igg
    d/dt(endoearlyfr1_muscle_igg) <- kon7 * c_endoearly_muscle_igg * c_fcrnendoearly_muscle * ve7_muscle - koff7 * c_endoearlyfr1_muscle_igg * ve7_muscle - kon7b * c_endoearlyfr1_muscle_igg * c_fcrnendoearly_muscle * ve7_muscle + koff7b * c_endoearlyfr2_muscle_igg * ve7_muscle - clup_muscle * c_endoearlyfr1_muscle_igg + clup_muscle * fr * c_memvasfr1_muscle_igg + clup_muscle * (1 - fr) * c_memintfr1_muscle_igg
    d/dt(endoearlyfr1_skin_igg) <- kon7 * c_endoearly_skin_igg * c_fcrnendoearly_skin * ve7_skin - koff7 * c_endoearlyfr1_skin_igg * ve7_skin - kon7b * c_endoearlyfr1_skin_igg * c_fcrnendoearly_skin * ve7_skin + koff7b * c_endoearlyfr2_skin_igg * ve7_skin - clup_skin * c_endoearlyfr1_skin_igg + clup_skin * fr * c_memvasfr1_skin_igg + clup_skin * (1 - fr) * c_memintfr1_skin_igg
    d/dt(endoearlyfr1_adipose_igg) <- kon7 * c_endoearly_adipose_igg * c_fcrnendoearly_adipose * ve7_adipose - koff7 * c_endoearlyfr1_adipose_igg * ve7_adipose - kon7b * c_endoearlyfr1_adipose_igg * c_fcrnendoearly_adipose * ve7_adipose + koff7b * c_endoearlyfr2_adipose_igg * ve7_adipose - clup_adipose * c_endoearlyfr1_adipose_igg + clup_adipose * fr * c_memvasfr1_adipose_igg + clup_adipose * (1 - fr) * c_memintfr1_adipose_igg
    d/dt(endoearlyfr1_bone_igg) <- kon7 * c_endoearly_bone_igg * c_fcrnendoearly_bone * ve7_bone - koff7 * c_endoearlyfr1_bone_igg * ve7_bone - kon7b * c_endoearlyfr1_bone_igg * c_fcrnendoearly_bone * ve7_bone + koff7b * c_endoearlyfr2_bone_igg * ve7_bone - clup_bone * c_endoearlyfr1_bone_igg + clup_bone * fr * c_memvasfr1_bone_igg + clup_bone * (1 - fr) * c_memintfr1_bone_igg
    d/dt(endoearlyfr1_brain_igg) <- kon7 * c_endoearly_brain_igg * c_fcrnendoearly_brain * ve7_brain - koff7 * c_endoearlyfr1_brain_igg * ve7_brain - kon7b * c_endoearlyfr1_brain_igg * c_fcrnendoearly_brain * ve7_brain + koff7b * c_endoearlyfr2_brain_igg * ve7_brain - clup_brain * c_endoearlyfr1_brain_igg + clup_brain * fr * c_memvasfr1_brain_igg + clup_brain * (1 - fr) * c_memintfr1_brain_igg
    d/dt(endoearlyfr1_kidney_igg) <- kon7 * c_endoearly_kidney_igg * c_fcrnendoearly_kidney * ve7_kidney - koff7 * c_endoearlyfr1_kidney_igg * ve7_kidney - kon7b * c_endoearlyfr1_kidney_igg * c_fcrnendoearly_kidney * ve7_kidney + koff7b * c_endoearlyfr2_kidney_igg * ve7_kidney - clup_kidney * c_endoearlyfr1_kidney_igg + clup_kidney * fr * c_memvasfr1_kidney_igg + clup_kidney * (1 - fr) * c_memintfr1_kidney_igg
    d/dt(endoearlyfr1_small_intestine_igg) <- kon7 * c_endoearly_small_intestine_igg * c_fcrnendoearly_small_intestine * ve7_small_intestine - koff7 * c_endoearlyfr1_small_intestine_igg * ve7_small_intestine - kon7b * c_endoearlyfr1_small_intestine_igg * c_fcrnendoearly_small_intestine * ve7_small_intestine + koff7b * c_endoearlyfr2_small_intestine_igg * ve7_small_intestine - clup_small_intestine * c_endoearlyfr1_small_intestine_igg + clup_small_intestine * fr * c_memvasfr1_small_intestine_igg + clup_small_intestine * (1 - fr) * c_memintfr1_small_intestine_igg
    d/dt(endoearlyfr1_large_intestine_igg) <- kon7 * c_endoearly_large_intestine_igg * c_fcrnendoearly_large_intestine * ve7_large_intestine - koff7 * c_endoearlyfr1_large_intestine_igg * ve7_large_intestine - kon7b * c_endoearlyfr1_large_intestine_igg * c_fcrnendoearly_large_intestine * ve7_large_intestine + koff7b * c_endoearlyfr2_large_intestine_igg * ve7_large_intestine - clup_large_intestine * c_endoearlyfr1_large_intestine_igg + clup_large_intestine * fr * c_memvasfr1_large_intestine_igg + clup_large_intestine * (1 - fr) * c_memintfr1_large_intestine_igg
    d/dt(endoearlyfr1_pancreas_igg) <- kon7 * c_endoearly_pancreas_igg * c_fcrnendoearly_pancreas * ve7_pancreas - koff7 * c_endoearlyfr1_pancreas_igg * ve7_pancreas - kon7b * c_endoearlyfr1_pancreas_igg * c_fcrnendoearly_pancreas * ve7_pancreas + koff7b * c_endoearlyfr2_pancreas_igg * ve7_pancreas - clup_pancreas * c_endoearlyfr1_pancreas_igg + clup_pancreas * fr * c_memvasfr1_pancreas_igg + clup_pancreas * (1 - fr) * c_memintfr1_pancreas_igg
    d/dt(endoearlyfr1_thymus_igg) <- kon7 * c_endoearly_thymus_igg * c_fcrnendoearly_thymus * ve7_thymus - koff7 * c_endoearlyfr1_thymus_igg * ve7_thymus - kon7b * c_endoearlyfr1_thymus_igg * c_fcrnendoearly_thymus * ve7_thymus + koff7b * c_endoearlyfr2_thymus_igg * ve7_thymus - clup_thymus * c_endoearlyfr1_thymus_igg + clup_thymus * fr * c_memvasfr1_thymus_igg + clup_thymus * (1 - fr) * c_memintfr1_thymus_igg
    d/dt(endoearlyfr1_spleen_igg) <- kon7 * c_endoearly_spleen_igg * c_fcrnendoearly_spleen * ve7_spleen - koff7 * c_endoearlyfr1_spleen_igg * ve7_spleen - kon7b * c_endoearlyfr1_spleen_igg * c_fcrnendoearly_spleen * ve7_spleen + koff7b * c_endoearlyfr2_spleen_igg * ve7_spleen - clup_spleen * c_endoearlyfr1_spleen_igg + clup_spleen * fr * c_memvasfr1_spleen_igg + clup_spleen * (1 - fr) * c_memintfr1_spleen_igg
    d/dt(endoearlyfr1_other_igg) <- kon7 * c_endoearly_other_igg * c_fcrnendoearly_other * ve7_other - koff7 * c_endoearlyfr1_other_igg * ve7_other - kon7b * c_endoearlyfr1_other_igg * c_fcrnendoearly_other * ve7_other + koff7b * c_endoearlyfr2_other_igg * ve7_other - clup_other * c_endoearlyfr1_other_igg + clup_other * fr * c_memvasfr1_other_igg + clup_other * (1 - fr) * c_memintfr1_other_igg
    d/dt(endosortfr1_lung_igg) <- kon6 * c_endosort_lung_igg * c_fcrnendosort_lung * ve6a_lung - koff6 * c_endosortfr1_lung_igg * ve6a_lung - kon6b * c_endosortfr1_lung_igg * c_fcrnendosort_lung * ve6a_lung + koff6b * c_endosortfr2_lung_igg * ve6a_lung + clup_lung * c_endoearlyfr1_lung_igg - clup_lung * c_endosortfr1_lung_igg
    d/dt(endosortfr1_liver_igg) <- kon6 * c_endosort_liver_igg * c_fcrnendosort_liver * ve6a_liver - koff6 * c_endosortfr1_liver_igg * ve6a_liver - kon6b * c_endosortfr1_liver_igg * c_fcrnendosort_liver * ve6a_liver + koff6b * c_endosortfr2_liver_igg * ve6a_liver + clup_liver * c_endoearlyfr1_liver_igg - clup_liver * c_endosortfr1_liver_igg
    d/dt(endosortfr1_heart_igg) <- kon6 * c_endosort_heart_igg * c_fcrnendosort_heart * ve6a_heart - koff6 * c_endosortfr1_heart_igg * ve6a_heart - kon6b * c_endosortfr1_heart_igg * c_fcrnendosort_heart * ve6a_heart + koff6b * c_endosortfr2_heart_igg * ve6a_heart + clup_heart * c_endoearlyfr1_heart_igg - clup_heart * c_endosortfr1_heart_igg
    d/dt(endosortfr1_muscle_igg) <- kon6 * c_endosort_muscle_igg * c_fcrnendosort_muscle * ve6a_muscle - koff6 * c_endosortfr1_muscle_igg * ve6a_muscle - kon6b * c_endosortfr1_muscle_igg * c_fcrnendosort_muscle * ve6a_muscle + koff6b * c_endosortfr2_muscle_igg * ve6a_muscle + clup_muscle * c_endoearlyfr1_muscle_igg - clup_muscle * c_endosortfr1_muscle_igg
    d/dt(endosortfr1_skin_igg) <- kon6 * c_endosort_skin_igg * c_fcrnendosort_skin * ve6a_skin - koff6 * c_endosortfr1_skin_igg * ve6a_skin - kon6b * c_endosortfr1_skin_igg * c_fcrnendosort_skin * ve6a_skin + koff6b * c_endosortfr2_skin_igg * ve6a_skin + clup_skin * c_endoearlyfr1_skin_igg - clup_skin * c_endosortfr1_skin_igg
    d/dt(endosortfr1_adipose_igg) <- kon6 * c_endosort_adipose_igg * c_fcrnendosort_adipose * ve6a_adipose - koff6 * c_endosortfr1_adipose_igg * ve6a_adipose - kon6b * c_endosortfr1_adipose_igg * c_fcrnendosort_adipose * ve6a_adipose + koff6b * c_endosortfr2_adipose_igg * ve6a_adipose + clup_adipose * c_endoearlyfr1_adipose_igg - clup_adipose * c_endosortfr1_adipose_igg
    d/dt(endosortfr1_bone_igg) <- kon6 * c_endosort_bone_igg * c_fcrnendosort_bone * ve6a_bone - koff6 * c_endosortfr1_bone_igg * ve6a_bone - kon6b * c_endosortfr1_bone_igg * c_fcrnendosort_bone * ve6a_bone + koff6b * c_endosortfr2_bone_igg * ve6a_bone + clup_bone * c_endoearlyfr1_bone_igg - clup_bone * c_endosortfr1_bone_igg
    d/dt(endosortfr1_brain_igg) <- kon6 * c_endosort_brain_igg * c_fcrnendosort_brain * ve6a_brain - koff6 * c_endosortfr1_brain_igg * ve6a_brain - kon6b * c_endosortfr1_brain_igg * c_fcrnendosort_brain * ve6a_brain + koff6b * c_endosortfr2_brain_igg * ve6a_brain + clup_brain * c_endoearlyfr1_brain_igg - clup_brain * c_endosortfr1_brain_igg
    d/dt(endosortfr1_kidney_igg) <- kon6 * c_endosort_kidney_igg * c_fcrnendosort_kidney * ve6a_kidney - koff6 * c_endosortfr1_kidney_igg * ve6a_kidney - kon6b * c_endosortfr1_kidney_igg * c_fcrnendosort_kidney * ve6a_kidney + koff6b * c_endosortfr2_kidney_igg * ve6a_kidney + clup_kidney * c_endoearlyfr1_kidney_igg - clup_kidney * c_endosortfr1_kidney_igg
    d/dt(endosortfr1_small_intestine_igg) <- kon6 * c_endosort_small_intestine_igg * c_fcrnendosort_small_intestine * ve6a_small_intestine - koff6 * c_endosortfr1_small_intestine_igg * ve6a_small_intestine - kon6b * c_endosortfr1_small_intestine_igg * c_fcrnendosort_small_intestine * ve6a_small_intestine + koff6b * c_endosortfr2_small_intestine_igg * ve6a_small_intestine + clup_small_intestine * c_endoearlyfr1_small_intestine_igg - clup_small_intestine * c_endosortfr1_small_intestine_igg
    d/dt(endosortfr1_large_intestine_igg) <- kon6 * c_endosort_large_intestine_igg * c_fcrnendosort_large_intestine * ve6a_large_intestine - koff6 * c_endosortfr1_large_intestine_igg * ve6a_large_intestine - kon6b * c_endosortfr1_large_intestine_igg * c_fcrnendosort_large_intestine * ve6a_large_intestine + koff6b * c_endosortfr2_large_intestine_igg * ve6a_large_intestine + clup_large_intestine * c_endoearlyfr1_large_intestine_igg - clup_large_intestine * c_endosortfr1_large_intestine_igg
    d/dt(endosortfr1_pancreas_igg) <- kon6 * c_endosort_pancreas_igg * c_fcrnendosort_pancreas * ve6a_pancreas - koff6 * c_endosortfr1_pancreas_igg * ve6a_pancreas - kon6b * c_endosortfr1_pancreas_igg * c_fcrnendosort_pancreas * ve6a_pancreas + koff6b * c_endosortfr2_pancreas_igg * ve6a_pancreas + clup_pancreas * c_endoearlyfr1_pancreas_igg - clup_pancreas * c_endosortfr1_pancreas_igg
    d/dt(endosortfr1_thymus_igg) <- kon6 * c_endosort_thymus_igg * c_fcrnendosort_thymus * ve6a_thymus - koff6 * c_endosortfr1_thymus_igg * ve6a_thymus - kon6b * c_endosortfr1_thymus_igg * c_fcrnendosort_thymus * ve6a_thymus + koff6b * c_endosortfr2_thymus_igg * ve6a_thymus + clup_thymus * c_endoearlyfr1_thymus_igg - clup_thymus * c_endosortfr1_thymus_igg
    d/dt(endosortfr1_spleen_igg) <- kon6 * c_endosort_spleen_igg * c_fcrnendosort_spleen * ve6a_spleen - koff6 * c_endosortfr1_spleen_igg * ve6a_spleen - kon6b * c_endosortfr1_spleen_igg * c_fcrnendosort_spleen * ve6a_spleen + koff6b * c_endosortfr2_spleen_igg * ve6a_spleen + clup_spleen * c_endoearlyfr1_spleen_igg - clup_spleen * c_endosortfr1_spleen_igg
    d/dt(endosortfr1_other_igg) <- kon6 * c_endosort_other_igg * c_fcrnendosort_other * ve6a_other - koff6 * c_endosortfr1_other_igg * ve6a_other - kon6b * c_endosortfr1_other_igg * c_fcrnendosort_other * ve6a_other + koff6b * c_endosortfr2_other_igg * ve6a_other + clup_other * c_endoearlyfr1_other_igg - clup_other * c_endosortfr1_other_igg
    d/dt(endorecycfr1_lung_igg) <- kon7 * c_endorecyc_lung_igg * c_fcrnendorecyc_lung * ve7b_lung - koff7 * c_endorecycfr1_lung_igg * ve7b_lung - kon7b * c_endorecycfr1_lung_igg * c_fcrnendorecyc_lung * ve7b_lung + koff7b * c_endorecycfr2_lung_igg * ve7b_lung + clup_lung * c_endosortfr1_lung_igg - clup_lung * c_endorecycfr1_lung_igg
    d/dt(endorecycfr1_liver_igg) <- kon7 * c_endorecyc_liver_igg * c_fcrnendorecyc_liver * ve7b_liver - koff7 * c_endorecycfr1_liver_igg * ve7b_liver - kon7b * c_endorecycfr1_liver_igg * c_fcrnendorecyc_liver * ve7b_liver + koff7b * c_endorecycfr2_liver_igg * ve7b_liver + clup_liver * c_endosortfr1_liver_igg - clup_liver * c_endorecycfr1_liver_igg
    d/dt(endorecycfr1_heart_igg) <- kon7 * c_endorecyc_heart_igg * c_fcrnendorecyc_heart * ve7b_heart - koff7 * c_endorecycfr1_heart_igg * ve7b_heart - kon7b * c_endorecycfr1_heart_igg * c_fcrnendorecyc_heart * ve7b_heart + koff7b * c_endorecycfr2_heart_igg * ve7b_heart + clup_heart * c_endosortfr1_heart_igg - clup_heart * c_endorecycfr1_heart_igg
    d/dt(endorecycfr1_muscle_igg) <- kon7 * c_endorecyc_muscle_igg * c_fcrnendorecyc_muscle * ve7b_muscle - koff7 * c_endorecycfr1_muscle_igg * ve7b_muscle - kon7b * c_endorecycfr1_muscle_igg * c_fcrnendorecyc_muscle * ve7b_muscle + koff7b * c_endorecycfr2_muscle_igg * ve7b_muscle + clup_muscle * c_endosortfr1_muscle_igg - clup_muscle * c_endorecycfr1_muscle_igg
    d/dt(endorecycfr1_skin_igg) <- kon7 * c_endorecyc_skin_igg * c_fcrnendorecyc_skin * ve7b_skin - koff7 * c_endorecycfr1_skin_igg * ve7b_skin - kon7b * c_endorecycfr1_skin_igg * c_fcrnendorecyc_skin * ve7b_skin + koff7b * c_endorecycfr2_skin_igg * ve7b_skin + clup_skin * c_endosortfr1_skin_igg - clup_skin * c_endorecycfr1_skin_igg
    d/dt(endorecycfr1_adipose_igg) <- kon7 * c_endorecyc_adipose_igg * c_fcrnendorecyc_adipose * ve7b_adipose - koff7 * c_endorecycfr1_adipose_igg * ve7b_adipose - kon7b * c_endorecycfr1_adipose_igg * c_fcrnendorecyc_adipose * ve7b_adipose + koff7b * c_endorecycfr2_adipose_igg * ve7b_adipose + clup_adipose * c_endosortfr1_adipose_igg - clup_adipose * c_endorecycfr1_adipose_igg
    d/dt(endorecycfr1_bone_igg) <- kon7 * c_endorecyc_bone_igg * c_fcrnendorecyc_bone * ve7b_bone - koff7 * c_endorecycfr1_bone_igg * ve7b_bone - kon7b * c_endorecycfr1_bone_igg * c_fcrnendorecyc_bone * ve7b_bone + koff7b * c_endorecycfr2_bone_igg * ve7b_bone + clup_bone * c_endosortfr1_bone_igg - clup_bone * c_endorecycfr1_bone_igg
    d/dt(endorecycfr1_brain_igg) <- kon7 * c_endorecyc_brain_igg * c_fcrnendorecyc_brain * ve7b_brain - koff7 * c_endorecycfr1_brain_igg * ve7b_brain - kon7b * c_endorecycfr1_brain_igg * c_fcrnendorecyc_brain * ve7b_brain + koff7b * c_endorecycfr2_brain_igg * ve7b_brain + clup_brain * c_endosortfr1_brain_igg - clup_brain * c_endorecycfr1_brain_igg
    d/dt(endorecycfr1_kidney_igg) <- kon7 * c_endorecyc_kidney_igg * c_fcrnendorecyc_kidney * ve7b_kidney - koff7 * c_endorecycfr1_kidney_igg * ve7b_kidney - kon7b * c_endorecycfr1_kidney_igg * c_fcrnendorecyc_kidney * ve7b_kidney + koff7b * c_endorecycfr2_kidney_igg * ve7b_kidney + clup_kidney * c_endosortfr1_kidney_igg - clup_kidney * c_endorecycfr1_kidney_igg
    d/dt(endorecycfr1_small_intestine_igg) <- kon7 * c_endorecyc_small_intestine_igg * c_fcrnendorecyc_small_intestine * ve7b_small_intestine - koff7 * c_endorecycfr1_small_intestine_igg * ve7b_small_intestine - kon7b * c_endorecycfr1_small_intestine_igg * c_fcrnendorecyc_small_intestine * ve7b_small_intestine + koff7b * c_endorecycfr2_small_intestine_igg * ve7b_small_intestine + clup_small_intestine * c_endosortfr1_small_intestine_igg - clup_small_intestine * c_endorecycfr1_small_intestine_igg
    d/dt(endorecycfr1_large_intestine_igg) <- kon7 * c_endorecyc_large_intestine_igg * c_fcrnendorecyc_large_intestine * ve7b_large_intestine - koff7 * c_endorecycfr1_large_intestine_igg * ve7b_large_intestine - kon7b * c_endorecycfr1_large_intestine_igg * c_fcrnendorecyc_large_intestine * ve7b_large_intestine + koff7b * c_endorecycfr2_large_intestine_igg * ve7b_large_intestine + clup_large_intestine * c_endosortfr1_large_intestine_igg - clup_large_intestine * c_endorecycfr1_large_intestine_igg
    d/dt(endorecycfr1_pancreas_igg) <- kon7 * c_endorecyc_pancreas_igg * c_fcrnendorecyc_pancreas * ve7b_pancreas - koff7 * c_endorecycfr1_pancreas_igg * ve7b_pancreas - kon7b * c_endorecycfr1_pancreas_igg * c_fcrnendorecyc_pancreas * ve7b_pancreas + koff7b * c_endorecycfr2_pancreas_igg * ve7b_pancreas + clup_pancreas * c_endosortfr1_pancreas_igg - clup_pancreas * c_endorecycfr1_pancreas_igg
    d/dt(endorecycfr1_thymus_igg) <- kon7 * c_endorecyc_thymus_igg * c_fcrnendorecyc_thymus * ve7b_thymus - koff7 * c_endorecycfr1_thymus_igg * ve7b_thymus - kon7b * c_endorecycfr1_thymus_igg * c_fcrnendorecyc_thymus * ve7b_thymus + koff7b * c_endorecycfr2_thymus_igg * ve7b_thymus + clup_thymus * c_endosortfr1_thymus_igg - clup_thymus * c_endorecycfr1_thymus_igg
    d/dt(endorecycfr1_spleen_igg) <- kon7 * c_endorecyc_spleen_igg * c_fcrnendorecyc_spleen * ve7b_spleen - koff7 * c_endorecycfr1_spleen_igg * ve7b_spleen - kon7b * c_endorecycfr1_spleen_igg * c_fcrnendorecyc_spleen * ve7b_spleen + koff7b * c_endorecycfr2_spleen_igg * ve7b_spleen + clup_spleen * c_endosortfr1_spleen_igg - clup_spleen * c_endorecycfr1_spleen_igg
    d/dt(endorecycfr1_other_igg) <- kon7 * c_endorecyc_other_igg * c_fcrnendorecyc_other * ve7b_other - koff7 * c_endorecycfr1_other_igg * ve7b_other - kon7b * c_endorecycfr1_other_igg * c_fcrnendorecyc_other * ve7b_other + koff7b * c_endorecycfr2_other_igg * ve7b_other + clup_other * c_endosortfr1_other_igg - clup_other * c_endorecycfr1_other_igg
    d/dt(memintfr1_lung_igg) <- kon7 * c_memint_lung_igg * c_fcrnmemint_lung * vism_lung - koff7 * c_memintfr1_lung_igg * vism_lung - kon7b * c_memintfr1_lung_igg * c_fcrnmemint_lung * vism_lung + koff7b * c_memintfr2_lung_igg * vism_lung + clup_lung * (1 - fr) * c_endorecycfr1_lung_igg - clup_lung * (1 - fr) * c_memintfr1_lung_igg - kdegab * c_memintfr1_lung_igg * vism_lung
    d/dt(memintfr1_liver_igg) <- kon7 * c_memint_liver_igg * c_fcrnmemint_liver * vism_liver - koff7 * c_memintfr1_liver_igg * vism_liver - kon7b * c_memintfr1_liver_igg * c_fcrnmemint_liver * vism_liver + koff7b * c_memintfr2_liver_igg * vism_liver + clup_liver * (1 - fr) * c_endorecycfr1_liver_igg - clup_liver * (1 - fr) * c_memintfr1_liver_igg - kdegab * c_memintfr1_liver_igg * vism_liver
    d/dt(memintfr1_heart_igg) <- kon7 * c_memint_heart_igg * c_fcrnmemint_heart * vism_heart - koff7 * c_memintfr1_heart_igg * vism_heart - kon7b * c_memintfr1_heart_igg * c_fcrnmemint_heart * vism_heart + koff7b * c_memintfr2_heart_igg * vism_heart + clup_heart * (1 - fr) * c_endorecycfr1_heart_igg - clup_heart * (1 - fr) * c_memintfr1_heart_igg - kdegab * c_memintfr1_heart_igg * vism_heart
    d/dt(memintfr1_muscle_igg) <- kon7 * c_memint_muscle_igg * c_fcrnmemint_muscle * vism_muscle - koff7 * c_memintfr1_muscle_igg * vism_muscle - kon7b * c_memintfr1_muscle_igg * c_fcrnmemint_muscle * vism_muscle + koff7b * c_memintfr2_muscle_igg * vism_muscle + clup_muscle * (1 - fr) * c_endorecycfr1_muscle_igg - clup_muscle * (1 - fr) * c_memintfr1_muscle_igg - kdegab * c_memintfr1_muscle_igg * vism_muscle
    d/dt(memintfr1_skin_igg) <- kon7 * c_memint_skin_igg * c_fcrnmemint_skin * vism_skin - koff7 * c_memintfr1_skin_igg * vism_skin - kon7b * c_memintfr1_skin_igg * c_fcrnmemint_skin * vism_skin + koff7b * c_memintfr2_skin_igg * vism_skin + clup_skin * (1 - fr) * c_endorecycfr1_skin_igg - clup_skin * (1 - fr) * c_memintfr1_skin_igg - kdegab * c_memintfr1_skin_igg * vism_skin
    d/dt(memintfr1_adipose_igg) <- kon7 * c_memint_adipose_igg * c_fcrnmemint_adipose * vism_adipose - koff7 * c_memintfr1_adipose_igg * vism_adipose - kon7b * c_memintfr1_adipose_igg * c_fcrnmemint_adipose * vism_adipose + koff7b * c_memintfr2_adipose_igg * vism_adipose + clup_adipose * (1 - fr) * c_endorecycfr1_adipose_igg - clup_adipose * (1 - fr) * c_memintfr1_adipose_igg - kdegab * c_memintfr1_adipose_igg * vism_adipose
    d/dt(memintfr1_bone_igg) <- kon7 * c_memint_bone_igg * c_fcrnmemint_bone * vism_bone - koff7 * c_memintfr1_bone_igg * vism_bone - kon7b * c_memintfr1_bone_igg * c_fcrnmemint_bone * vism_bone + koff7b * c_memintfr2_bone_igg * vism_bone + clup_bone * (1 - fr) * c_endorecycfr1_bone_igg - clup_bone * (1 - fr) * c_memintfr1_bone_igg - kdegab * c_memintfr1_bone_igg * vism_bone
    d/dt(memintfr1_brain_igg) <- kon7 * c_memint_brain_igg * c_fcrnmemint_brain * vism_brain - koff7 * c_memintfr1_brain_igg * vism_brain - kon7b * c_memintfr1_brain_igg * c_fcrnmemint_brain * vism_brain + koff7b * c_memintfr2_brain_igg * vism_brain + clup_brain * (1 - fr) * c_endorecycfr1_brain_igg - clup_brain * (1 - fr) * c_memintfr1_brain_igg - kdegab * c_memintfr1_brain_igg * vism_brain
    d/dt(memintfr1_kidney_igg) <- kon7 * c_memint_kidney_igg * c_fcrnmemint_kidney * vism_kidney - koff7 * c_memintfr1_kidney_igg * vism_kidney - kon7b * c_memintfr1_kidney_igg * c_fcrnmemint_kidney * vism_kidney + koff7b * c_memintfr2_kidney_igg * vism_kidney + clup_kidney * (1 - fr) * c_endorecycfr1_kidney_igg - clup_kidney * (1 - fr) * c_memintfr1_kidney_igg - kdegab * c_memintfr1_kidney_igg * vism_kidney
    d/dt(memintfr1_small_intestine_igg) <- kon7 * c_memint_small_intestine_igg * c_fcrnmemint_small_intestine * vism_small_intestine - koff7 * c_memintfr1_small_intestine_igg * vism_small_intestine - kon7b * c_memintfr1_small_intestine_igg * c_fcrnmemint_small_intestine * vism_small_intestine + koff7b * c_memintfr2_small_intestine_igg * vism_small_intestine + clup_small_intestine * (1 - fr) * c_endorecycfr1_small_intestine_igg - clup_small_intestine * (1 - fr) * c_memintfr1_small_intestine_igg - kdegab * c_memintfr1_small_intestine_igg * vism_small_intestine
    d/dt(memintfr1_large_intestine_igg) <- kon7 * c_memint_large_intestine_igg * c_fcrnmemint_large_intestine * vism_large_intestine - koff7 * c_memintfr1_large_intestine_igg * vism_large_intestine - kon7b * c_memintfr1_large_intestine_igg * c_fcrnmemint_large_intestine * vism_large_intestine + koff7b * c_memintfr2_large_intestine_igg * vism_large_intestine + clup_large_intestine * (1 - fr) * c_endorecycfr1_large_intestine_igg - clup_large_intestine * (1 - fr) * c_memintfr1_large_intestine_igg - kdegab * c_memintfr1_large_intestine_igg * vism_large_intestine
    d/dt(memintfr1_pancreas_igg) <- kon7 * c_memint_pancreas_igg * c_fcrnmemint_pancreas * vism_pancreas - koff7 * c_memintfr1_pancreas_igg * vism_pancreas - kon7b * c_memintfr1_pancreas_igg * c_fcrnmemint_pancreas * vism_pancreas + koff7b * c_memintfr2_pancreas_igg * vism_pancreas + clup_pancreas * (1 - fr) * c_endorecycfr1_pancreas_igg - clup_pancreas * (1 - fr) * c_memintfr1_pancreas_igg - kdegab * c_memintfr1_pancreas_igg * vism_pancreas
    d/dt(memintfr1_thymus_igg) <- kon7 * c_memint_thymus_igg * c_fcrnmemint_thymus * vism_thymus - koff7 * c_memintfr1_thymus_igg * vism_thymus - kon7b * c_memintfr1_thymus_igg * c_fcrnmemint_thymus * vism_thymus + koff7b * c_memintfr2_thymus_igg * vism_thymus + clup_thymus * (1 - fr) * c_endorecycfr1_thymus_igg - clup_thymus * (1 - fr) * c_memintfr1_thymus_igg - kdegab * c_memintfr1_thymus_igg * vism_thymus
    d/dt(memintfr1_spleen_igg) <- kon7 * c_memint_spleen_igg * c_fcrnmemint_spleen * vism_spleen - koff7 * c_memintfr1_spleen_igg * vism_spleen - kon7b * c_memintfr1_spleen_igg * c_fcrnmemint_spleen * vism_spleen + koff7b * c_memintfr2_spleen_igg * vism_spleen + clup_spleen * (1 - fr) * c_endorecycfr1_spleen_igg - clup_spleen * (1 - fr) * c_memintfr1_spleen_igg - kdegab * c_memintfr1_spleen_igg * vism_spleen
    d/dt(memintfr1_other_igg) <- kon7 * c_memint_other_igg * c_fcrnmemint_other * vism_other - koff7 * c_memintfr1_other_igg * vism_other - kon7b * c_memintfr1_other_igg * c_fcrnmemint_other * vism_other + koff7b * c_memintfr2_other_igg * vism_other + clup_other * (1 - fr) * c_endorecycfr1_other_igg - clup_other * (1 - fr) * c_memintfr1_other_igg - kdegab * c_memintfr1_other_igg * vism_other

    # ---- igg: 2:1 FcRn complexes ----
    d/dt(memvasfr2_lung_igg) <- kon7b * c_memvasfr1_lung_igg * c_fcrnmemvas_lung * vvm_lung - koff7b * c_memvasfr2_lung_igg * vvm_lung + clup_lung * fr * c_endorecycfr2_lung_igg - clup_lung * fr * c_memvasfr2_lung_igg - kdegab * c_memvasfr2_lung_igg * vvm_lung
    d/dt(memvasfr2_liver_igg) <- kon7b * c_memvasfr1_liver_igg * c_fcrnmemvas_liver * vvm_liver - koff7b * c_memvasfr2_liver_igg * vvm_liver + clup_liver * fr * c_endorecycfr2_liver_igg - clup_liver * fr * c_memvasfr2_liver_igg - kdegab * c_memvasfr2_liver_igg * vvm_liver
    d/dt(memvasfr2_heart_igg) <- kon7b * c_memvasfr1_heart_igg * c_fcrnmemvas_heart * vvm_heart - koff7b * c_memvasfr2_heart_igg * vvm_heart + clup_heart * fr * c_endorecycfr2_heart_igg - clup_heart * fr * c_memvasfr2_heart_igg - kdegab * c_memvasfr2_heart_igg * vvm_heart
    d/dt(memvasfr2_muscle_igg) <- kon7b * c_memvasfr1_muscle_igg * c_fcrnmemvas_muscle * vvm_muscle - koff7b * c_memvasfr2_muscle_igg * vvm_muscle + clup_muscle * fr * c_endorecycfr2_muscle_igg - clup_muscle * fr * c_memvasfr2_muscle_igg - kdegab * c_memvasfr2_muscle_igg * vvm_muscle
    d/dt(memvasfr2_skin_igg) <- kon7b * c_memvasfr1_skin_igg * c_fcrnmemvas_skin * vvm_skin - koff7b * c_memvasfr2_skin_igg * vvm_skin + clup_skin * fr * c_endorecycfr2_skin_igg - clup_skin * fr * c_memvasfr2_skin_igg - kdegab * c_memvasfr2_skin_igg * vvm_skin
    d/dt(memvasfr2_adipose_igg) <- kon7b * c_memvasfr1_adipose_igg * c_fcrnmemvas_adipose * vvm_adipose - koff7b * c_memvasfr2_adipose_igg * vvm_adipose + clup_adipose * fr * c_endorecycfr2_adipose_igg - clup_adipose * fr * c_memvasfr2_adipose_igg - kdegab * c_memvasfr2_adipose_igg * vvm_adipose
    d/dt(memvasfr2_bone_igg) <- kon7b * c_memvasfr1_bone_igg * c_fcrnmemvas_bone * vvm_bone - koff7b * c_memvasfr2_bone_igg * vvm_bone + clup_bone * fr * c_endorecycfr2_bone_igg - clup_bone * fr * c_memvasfr2_bone_igg - kdegab * c_memvasfr2_bone_igg * vvm_bone
    d/dt(memvasfr2_brain_igg) <- kon7b * c_memvasfr1_brain_igg * c_fcrnmemvas_brain * vvm_brain - koff7b * c_memvasfr2_brain_igg * vvm_brain + clup_brain * fr * c_endorecycfr2_brain_igg - clup_brain * fr * c_memvasfr2_brain_igg - kdegab * c_memvasfr2_brain_igg * vvm_brain
    d/dt(memvasfr2_kidney_igg) <- kon7b * c_memvasfr1_kidney_igg * c_fcrnmemvas_kidney * vvm_kidney - koff7b * c_memvasfr2_kidney_igg * vvm_kidney + clup_kidney * fr * c_endorecycfr2_kidney_igg - clup_kidney * fr * c_memvasfr2_kidney_igg - kdegab * c_memvasfr2_kidney_igg * vvm_kidney
    d/dt(memvasfr2_small_intestine_igg) <- kon7b * c_memvasfr1_small_intestine_igg * c_fcrnmemvas_small_intestine * vvm_small_intestine - koff7b * c_memvasfr2_small_intestine_igg * vvm_small_intestine + clup_small_intestine * fr * c_endorecycfr2_small_intestine_igg - clup_small_intestine * fr * c_memvasfr2_small_intestine_igg - kdegab * c_memvasfr2_small_intestine_igg * vvm_small_intestine
    d/dt(memvasfr2_large_intestine_igg) <- kon7b * c_memvasfr1_large_intestine_igg * c_fcrnmemvas_large_intestine * vvm_large_intestine - koff7b * c_memvasfr2_large_intestine_igg * vvm_large_intestine + clup_large_intestine * fr * c_endorecycfr2_large_intestine_igg - clup_large_intestine * fr * c_memvasfr2_large_intestine_igg - kdegab * c_memvasfr2_large_intestine_igg * vvm_large_intestine
    d/dt(memvasfr2_pancreas_igg) <- kon7b * c_memvasfr1_pancreas_igg * c_fcrnmemvas_pancreas * vvm_pancreas - koff7b * c_memvasfr2_pancreas_igg * vvm_pancreas + clup_pancreas * fr * c_endorecycfr2_pancreas_igg - clup_pancreas * fr * c_memvasfr2_pancreas_igg - kdegab * c_memvasfr2_pancreas_igg * vvm_pancreas
    d/dt(memvasfr2_thymus_igg) <- kon7b * c_memvasfr1_thymus_igg * c_fcrnmemvas_thymus * vvm_thymus - koff7b * c_memvasfr2_thymus_igg * vvm_thymus + clup_thymus * fr * c_endorecycfr2_thymus_igg - clup_thymus * fr * c_memvasfr2_thymus_igg - kdegab * c_memvasfr2_thymus_igg * vvm_thymus
    d/dt(memvasfr2_spleen_igg) <- kon7b * c_memvasfr1_spleen_igg * c_fcrnmemvas_spleen * vvm_spleen - koff7b * c_memvasfr2_spleen_igg * vvm_spleen + clup_spleen * fr * c_endorecycfr2_spleen_igg - clup_spleen * fr * c_memvasfr2_spleen_igg - kdegab * c_memvasfr2_spleen_igg * vvm_spleen
    d/dt(memvasfr2_other_igg) <- kon7b * c_memvasfr1_other_igg * c_fcrnmemvas_other * vvm_other - koff7b * c_memvasfr2_other_igg * vvm_other + clup_other * fr * c_endorecycfr2_other_igg - clup_other * fr * c_memvasfr2_other_igg - kdegab * c_memvasfr2_other_igg * vvm_other
    d/dt(endoearlyfr2_lung_igg) <- kon7b * c_endoearlyfr1_lung_igg * c_fcrnendoearly_lung * ve7_lung - koff7b * c_endoearlyfr2_lung_igg * ve7_lung - clup_lung * c_endoearlyfr2_lung_igg + clup_lung * fr * c_memvasfr2_lung_igg + clup_lung * (1 - fr) * c_memintfr2_lung_igg
    d/dt(endoearlyfr2_liver_igg) <- kon7b * c_endoearlyfr1_liver_igg * c_fcrnendoearly_liver * ve7_liver - koff7b * c_endoearlyfr2_liver_igg * ve7_liver - clup_liver * c_endoearlyfr2_liver_igg + clup_liver * fr * c_memvasfr2_liver_igg + clup_liver * (1 - fr) * c_memintfr2_liver_igg
    d/dt(endoearlyfr2_heart_igg) <- kon7b * c_endoearlyfr1_heart_igg * c_fcrnendoearly_heart * ve7_heart - koff7b * c_endoearlyfr2_heart_igg * ve7_heart - clup_heart * c_endoearlyfr2_heart_igg + clup_heart * fr * c_memvasfr2_heart_igg + clup_heart * (1 - fr) * c_memintfr2_heart_igg
    d/dt(endoearlyfr2_muscle_igg) <- kon7b * c_endoearlyfr1_muscle_igg * c_fcrnendoearly_muscle * ve7_muscle - koff7b * c_endoearlyfr2_muscle_igg * ve7_muscle - clup_muscle * c_endoearlyfr2_muscle_igg + clup_muscle * fr * c_memvasfr2_muscle_igg + clup_muscle * (1 - fr) * c_memintfr2_muscle_igg
    d/dt(endoearlyfr2_skin_igg) <- kon7b * c_endoearlyfr1_skin_igg * c_fcrnendoearly_skin * ve7_skin - koff7b * c_endoearlyfr2_skin_igg * ve7_skin - clup_skin * c_endoearlyfr2_skin_igg + clup_skin * fr * c_memvasfr2_skin_igg + clup_skin * (1 - fr) * c_memintfr2_skin_igg
    d/dt(endoearlyfr2_adipose_igg) <- kon7b * c_endoearlyfr1_adipose_igg * c_fcrnendoearly_adipose * ve7_adipose - koff7b * c_endoearlyfr2_adipose_igg * ve7_adipose - clup_adipose * c_endoearlyfr2_adipose_igg + clup_adipose * fr * c_memvasfr2_adipose_igg + clup_adipose * (1 - fr) * c_memintfr2_adipose_igg
    d/dt(endoearlyfr2_bone_igg) <- kon7b * c_endoearlyfr1_bone_igg * c_fcrnendoearly_bone * ve7_bone - koff7b * c_endoearlyfr2_bone_igg * ve7_bone - clup_bone * c_endoearlyfr2_bone_igg + clup_bone * fr * c_memvasfr2_bone_igg + clup_bone * (1 - fr) * c_memintfr2_bone_igg
    d/dt(endoearlyfr2_brain_igg) <- kon7b * c_endoearlyfr1_brain_igg * c_fcrnendoearly_brain * ve7_brain - koff7b * c_endoearlyfr2_brain_igg * ve7_brain - clup_brain * c_endoearlyfr2_brain_igg + clup_brain * fr * c_memvasfr2_brain_igg + clup_brain * (1 - fr) * c_memintfr2_brain_igg
    d/dt(endoearlyfr2_kidney_igg) <- kon7b * c_endoearlyfr1_kidney_igg * c_fcrnendoearly_kidney * ve7_kidney - koff7b * c_endoearlyfr2_kidney_igg * ve7_kidney - clup_kidney * c_endoearlyfr2_kidney_igg + clup_kidney * fr * c_memvasfr2_kidney_igg + clup_kidney * (1 - fr) * c_memintfr2_kidney_igg
    d/dt(endoearlyfr2_small_intestine_igg) <- kon7b * c_endoearlyfr1_small_intestine_igg * c_fcrnendoearly_small_intestine * ve7_small_intestine - koff7b * c_endoearlyfr2_small_intestine_igg * ve7_small_intestine - clup_small_intestine * c_endoearlyfr2_small_intestine_igg + clup_small_intestine * fr * c_memvasfr2_small_intestine_igg + clup_small_intestine * (1 - fr) * c_memintfr2_small_intestine_igg
    d/dt(endoearlyfr2_large_intestine_igg) <- kon7b * c_endoearlyfr1_large_intestine_igg * c_fcrnendoearly_large_intestine * ve7_large_intestine - koff7b * c_endoearlyfr2_large_intestine_igg * ve7_large_intestine - clup_large_intestine * c_endoearlyfr2_large_intestine_igg + clup_large_intestine * fr * c_memvasfr2_large_intestine_igg + clup_large_intestine * (1 - fr) * c_memintfr2_large_intestine_igg
    d/dt(endoearlyfr2_pancreas_igg) <- kon7b * c_endoearlyfr1_pancreas_igg * c_fcrnendoearly_pancreas * ve7_pancreas - koff7b * c_endoearlyfr2_pancreas_igg * ve7_pancreas - clup_pancreas * c_endoearlyfr2_pancreas_igg + clup_pancreas * fr * c_memvasfr2_pancreas_igg + clup_pancreas * (1 - fr) * c_memintfr2_pancreas_igg
    d/dt(endoearlyfr2_thymus_igg) <- kon7b * c_endoearlyfr1_thymus_igg * c_fcrnendoearly_thymus * ve7_thymus - koff7b * c_endoearlyfr2_thymus_igg * ve7_thymus - clup_thymus * c_endoearlyfr2_thymus_igg + clup_thymus * fr * c_memvasfr2_thymus_igg + clup_thymus * (1 - fr) * c_memintfr2_thymus_igg
    d/dt(endoearlyfr2_spleen_igg) <- kon7b * c_endoearlyfr1_spleen_igg * c_fcrnendoearly_spleen * ve7_spleen - koff7b * c_endoearlyfr2_spleen_igg * ve7_spleen - clup_spleen * c_endoearlyfr2_spleen_igg + clup_spleen * fr * c_memvasfr2_spleen_igg + clup_spleen * (1 - fr) * c_memintfr2_spleen_igg
    d/dt(endoearlyfr2_other_igg) <- kon7b * c_endoearlyfr1_other_igg * c_fcrnendoearly_other * ve7_other - koff7b * c_endoearlyfr2_other_igg * ve7_other - clup_other * c_endoearlyfr2_other_igg + clup_other * fr * c_memvasfr2_other_igg + clup_other * (1 - fr) * c_memintfr2_other_igg
    d/dt(endosortfr2_lung_igg) <- kon6b * c_endosortfr1_lung_igg * c_fcrnendosort_lung * ve6a_lung - koff6b * c_endosortfr2_lung_igg * ve6a_lung + clup_lung * c_endoearlyfr2_lung_igg - clup_lung * c_endosortfr2_lung_igg
    d/dt(endosortfr2_liver_igg) <- kon6b * c_endosortfr1_liver_igg * c_fcrnendosort_liver * ve6a_liver - koff6b * c_endosortfr2_liver_igg * ve6a_liver + clup_liver * c_endoearlyfr2_liver_igg - clup_liver * c_endosortfr2_liver_igg
    d/dt(endosortfr2_heart_igg) <- kon6b * c_endosortfr1_heart_igg * c_fcrnendosort_heart * ve6a_heart - koff6b * c_endosortfr2_heart_igg * ve6a_heart + clup_heart * c_endoearlyfr2_heart_igg - clup_heart * c_endosortfr2_heart_igg
    d/dt(endosortfr2_muscle_igg) <- kon6b * c_endosortfr1_muscle_igg * c_fcrnendosort_muscle * ve6a_muscle - koff6b * c_endosortfr2_muscle_igg * ve6a_muscle + clup_muscle * c_endoearlyfr2_muscle_igg - clup_muscle * c_endosortfr2_muscle_igg
    d/dt(endosortfr2_skin_igg) <- kon6b * c_endosortfr1_skin_igg * c_fcrnendosort_skin * ve6a_skin - koff6b * c_endosortfr2_skin_igg * ve6a_skin + clup_skin * c_endoearlyfr2_skin_igg - clup_skin * c_endosortfr2_skin_igg
    d/dt(endosortfr2_adipose_igg) <- kon6b * c_endosortfr1_adipose_igg * c_fcrnendosort_adipose * ve6a_adipose - koff6b * c_endosortfr2_adipose_igg * ve6a_adipose + clup_adipose * c_endoearlyfr2_adipose_igg - clup_adipose * c_endosortfr2_adipose_igg
    d/dt(endosortfr2_bone_igg) <- kon6b * c_endosortfr1_bone_igg * c_fcrnendosort_bone * ve6a_bone - koff6b * c_endosortfr2_bone_igg * ve6a_bone + clup_bone * c_endoearlyfr2_bone_igg - clup_bone * c_endosortfr2_bone_igg
    d/dt(endosortfr2_brain_igg) <- kon6b * c_endosortfr1_brain_igg * c_fcrnendosort_brain * ve6a_brain - koff6b * c_endosortfr2_brain_igg * ve6a_brain + clup_brain * c_endoearlyfr2_brain_igg - clup_brain * c_endosortfr2_brain_igg
    d/dt(endosortfr2_kidney_igg) <- kon6b * c_endosortfr1_kidney_igg * c_fcrnendosort_kidney * ve6a_kidney - koff6b * c_endosortfr2_kidney_igg * ve6a_kidney + clup_kidney * c_endoearlyfr2_kidney_igg - clup_kidney * c_endosortfr2_kidney_igg
    d/dt(endosortfr2_small_intestine_igg) <- kon6b * c_endosortfr1_small_intestine_igg * c_fcrnendosort_small_intestine * ve6a_small_intestine - koff6b * c_endosortfr2_small_intestine_igg * ve6a_small_intestine + clup_small_intestine * c_endoearlyfr2_small_intestine_igg - clup_small_intestine * c_endosortfr2_small_intestine_igg
    d/dt(endosortfr2_large_intestine_igg) <- kon6b * c_endosortfr1_large_intestine_igg * c_fcrnendosort_large_intestine * ve6a_large_intestine - koff6b * c_endosortfr2_large_intestine_igg * ve6a_large_intestine + clup_large_intestine * c_endoearlyfr2_large_intestine_igg - clup_large_intestine * c_endosortfr2_large_intestine_igg
    d/dt(endosortfr2_pancreas_igg) <- kon6b * c_endosortfr1_pancreas_igg * c_fcrnendosort_pancreas * ve6a_pancreas - koff6b * c_endosortfr2_pancreas_igg * ve6a_pancreas + clup_pancreas * c_endoearlyfr2_pancreas_igg - clup_pancreas * c_endosortfr2_pancreas_igg
    d/dt(endosortfr2_thymus_igg) <- kon6b * c_endosortfr1_thymus_igg * c_fcrnendosort_thymus * ve6a_thymus - koff6b * c_endosortfr2_thymus_igg * ve6a_thymus + clup_thymus * c_endoearlyfr2_thymus_igg - clup_thymus * c_endosortfr2_thymus_igg
    d/dt(endosortfr2_spleen_igg) <- kon6b * c_endosortfr1_spleen_igg * c_fcrnendosort_spleen * ve6a_spleen - koff6b * c_endosortfr2_spleen_igg * ve6a_spleen + clup_spleen * c_endoearlyfr2_spleen_igg - clup_spleen * c_endosortfr2_spleen_igg
    d/dt(endosortfr2_other_igg) <- kon6b * c_endosortfr1_other_igg * c_fcrnendosort_other * ve6a_other - koff6b * c_endosortfr2_other_igg * ve6a_other + clup_other * c_endoearlyfr2_other_igg - clup_other * c_endosortfr2_other_igg
    d/dt(endorecycfr2_lung_igg) <- kon7b * c_endorecycfr1_lung_igg * c_fcrnendorecyc_lung * ve7b_lung - koff7b * c_endorecycfr2_lung_igg * ve7b_lung + clup_lung * c_endosortfr2_lung_igg - clup_lung * c_endorecycfr2_lung_igg
    d/dt(endorecycfr2_liver_igg) <- kon7b * c_endorecycfr1_liver_igg * c_fcrnendorecyc_liver * ve7b_liver - koff7b * c_endorecycfr2_liver_igg * ve7b_liver + clup_liver * c_endosortfr2_liver_igg - clup_liver * c_endorecycfr2_liver_igg
    d/dt(endorecycfr2_heart_igg) <- kon7b * c_endorecycfr1_heart_igg * c_fcrnendorecyc_heart * ve7b_heart - koff7b * c_endorecycfr2_heart_igg * ve7b_heart + clup_heart * c_endosortfr2_heart_igg - clup_heart * c_endorecycfr2_heart_igg
    d/dt(endorecycfr2_muscle_igg) <- kon7b * c_endorecycfr1_muscle_igg * c_fcrnendorecyc_muscle * ve7b_muscle - koff7b * c_endorecycfr2_muscle_igg * ve7b_muscle + clup_muscle * c_endosortfr2_muscle_igg - clup_muscle * c_endorecycfr2_muscle_igg
    d/dt(endorecycfr2_skin_igg) <- kon7b * c_endorecycfr1_skin_igg * c_fcrnendorecyc_skin * ve7b_skin - koff7b * c_endorecycfr2_skin_igg * ve7b_skin + clup_skin * c_endosortfr2_skin_igg - clup_skin * c_endorecycfr2_skin_igg
    d/dt(endorecycfr2_adipose_igg) <- kon7b * c_endorecycfr1_adipose_igg * c_fcrnendorecyc_adipose * ve7b_adipose - koff7b * c_endorecycfr2_adipose_igg * ve7b_adipose + clup_adipose * c_endosortfr2_adipose_igg - clup_adipose * c_endorecycfr2_adipose_igg
    d/dt(endorecycfr2_bone_igg) <- kon7b * c_endorecycfr1_bone_igg * c_fcrnendorecyc_bone * ve7b_bone - koff7b * c_endorecycfr2_bone_igg * ve7b_bone + clup_bone * c_endosortfr2_bone_igg - clup_bone * c_endorecycfr2_bone_igg
    d/dt(endorecycfr2_brain_igg) <- kon7b * c_endorecycfr1_brain_igg * c_fcrnendorecyc_brain * ve7b_brain - koff7b * c_endorecycfr2_brain_igg * ve7b_brain + clup_brain * c_endosortfr2_brain_igg - clup_brain * c_endorecycfr2_brain_igg
    d/dt(endorecycfr2_kidney_igg) <- kon7b * c_endorecycfr1_kidney_igg * c_fcrnendorecyc_kidney * ve7b_kidney - koff7b * c_endorecycfr2_kidney_igg * ve7b_kidney + clup_kidney * c_endosortfr2_kidney_igg - clup_kidney * c_endorecycfr2_kidney_igg
    d/dt(endorecycfr2_small_intestine_igg) <- kon7b * c_endorecycfr1_small_intestine_igg * c_fcrnendorecyc_small_intestine * ve7b_small_intestine - koff7b * c_endorecycfr2_small_intestine_igg * ve7b_small_intestine + clup_small_intestine * c_endosortfr2_small_intestine_igg - clup_small_intestine * c_endorecycfr2_small_intestine_igg
    d/dt(endorecycfr2_large_intestine_igg) <- kon7b * c_endorecycfr1_large_intestine_igg * c_fcrnendorecyc_large_intestine * ve7b_large_intestine - koff7b * c_endorecycfr2_large_intestine_igg * ve7b_large_intestine + clup_large_intestine * c_endosortfr2_large_intestine_igg - clup_large_intestine * c_endorecycfr2_large_intestine_igg
    d/dt(endorecycfr2_pancreas_igg) <- kon7b * c_endorecycfr1_pancreas_igg * c_fcrnendorecyc_pancreas * ve7b_pancreas - koff7b * c_endorecycfr2_pancreas_igg * ve7b_pancreas + clup_pancreas * c_endosortfr2_pancreas_igg - clup_pancreas * c_endorecycfr2_pancreas_igg
    d/dt(endorecycfr2_thymus_igg) <- kon7b * c_endorecycfr1_thymus_igg * c_fcrnendorecyc_thymus * ve7b_thymus - koff7b * c_endorecycfr2_thymus_igg * ve7b_thymus + clup_thymus * c_endosortfr2_thymus_igg - clup_thymus * c_endorecycfr2_thymus_igg
    d/dt(endorecycfr2_spleen_igg) <- kon7b * c_endorecycfr1_spleen_igg * c_fcrnendorecyc_spleen * ve7b_spleen - koff7b * c_endorecycfr2_spleen_igg * ve7b_spleen + clup_spleen * c_endosortfr2_spleen_igg - clup_spleen * c_endorecycfr2_spleen_igg
    d/dt(endorecycfr2_other_igg) <- kon7b * c_endorecycfr1_other_igg * c_fcrnendorecyc_other * ve7b_other - koff7b * c_endorecycfr2_other_igg * ve7b_other + clup_other * c_endosortfr2_other_igg - clup_other * c_endorecycfr2_other_igg
    d/dt(memintfr2_lung_igg) <- kon7b * c_memintfr1_lung_igg * c_fcrnmemint_lung * vism_lung - koff7b * c_memintfr2_lung_igg * vism_lung + clup_lung * (1 - fr) * c_endorecycfr2_lung_igg - clup_lung * (1 - fr) * c_memintfr2_lung_igg - kdegab * c_memintfr2_lung_igg * vism_lung
    d/dt(memintfr2_liver_igg) <- kon7b * c_memintfr1_liver_igg * c_fcrnmemint_liver * vism_liver - koff7b * c_memintfr2_liver_igg * vism_liver + clup_liver * (1 - fr) * c_endorecycfr2_liver_igg - clup_liver * (1 - fr) * c_memintfr2_liver_igg - kdegab * c_memintfr2_liver_igg * vism_liver
    d/dt(memintfr2_heart_igg) <- kon7b * c_memintfr1_heart_igg * c_fcrnmemint_heart * vism_heart - koff7b * c_memintfr2_heart_igg * vism_heart + clup_heart * (1 - fr) * c_endorecycfr2_heart_igg - clup_heart * (1 - fr) * c_memintfr2_heart_igg - kdegab * c_memintfr2_heart_igg * vism_heart
    d/dt(memintfr2_muscle_igg) <- kon7b * c_memintfr1_muscle_igg * c_fcrnmemint_muscle * vism_muscle - koff7b * c_memintfr2_muscle_igg * vism_muscle + clup_muscle * (1 - fr) * c_endorecycfr2_muscle_igg - clup_muscle * (1 - fr) * c_memintfr2_muscle_igg - kdegab * c_memintfr2_muscle_igg * vism_muscle
    d/dt(memintfr2_skin_igg) <- kon7b * c_memintfr1_skin_igg * c_fcrnmemint_skin * vism_skin - koff7b * c_memintfr2_skin_igg * vism_skin + clup_skin * (1 - fr) * c_endorecycfr2_skin_igg - clup_skin * (1 - fr) * c_memintfr2_skin_igg - kdegab * c_memintfr2_skin_igg * vism_skin
    d/dt(memintfr2_adipose_igg) <- kon7b * c_memintfr1_adipose_igg * c_fcrnmemint_adipose * vism_adipose - koff7b * c_memintfr2_adipose_igg * vism_adipose + clup_adipose * (1 - fr) * c_endorecycfr2_adipose_igg - clup_adipose * (1 - fr) * c_memintfr2_adipose_igg - kdegab * c_memintfr2_adipose_igg * vism_adipose
    d/dt(memintfr2_bone_igg) <- kon7b * c_memintfr1_bone_igg * c_fcrnmemint_bone * vism_bone - koff7b * c_memintfr2_bone_igg * vism_bone + clup_bone * (1 - fr) * c_endorecycfr2_bone_igg - clup_bone * (1 - fr) * c_memintfr2_bone_igg - kdegab * c_memintfr2_bone_igg * vism_bone
    d/dt(memintfr2_brain_igg) <- kon7b * c_memintfr1_brain_igg * c_fcrnmemint_brain * vism_brain - koff7b * c_memintfr2_brain_igg * vism_brain + clup_brain * (1 - fr) * c_endorecycfr2_brain_igg - clup_brain * (1 - fr) * c_memintfr2_brain_igg - kdegab * c_memintfr2_brain_igg * vism_brain
    d/dt(memintfr2_kidney_igg) <- kon7b * c_memintfr1_kidney_igg * c_fcrnmemint_kidney * vism_kidney - koff7b * c_memintfr2_kidney_igg * vism_kidney + clup_kidney * (1 - fr) * c_endorecycfr2_kidney_igg - clup_kidney * (1 - fr) * c_memintfr2_kidney_igg - kdegab * c_memintfr2_kidney_igg * vism_kidney
    d/dt(memintfr2_small_intestine_igg) <- kon7b * c_memintfr1_small_intestine_igg * c_fcrnmemint_small_intestine * vism_small_intestine - koff7b * c_memintfr2_small_intestine_igg * vism_small_intestine + clup_small_intestine * (1 - fr) * c_endorecycfr2_small_intestine_igg - clup_small_intestine * (1 - fr) * c_memintfr2_small_intestine_igg - kdegab * c_memintfr2_small_intestine_igg * vism_small_intestine
    d/dt(memintfr2_large_intestine_igg) <- kon7b * c_memintfr1_large_intestine_igg * c_fcrnmemint_large_intestine * vism_large_intestine - koff7b * c_memintfr2_large_intestine_igg * vism_large_intestine + clup_large_intestine * (1 - fr) * c_endorecycfr2_large_intestine_igg - clup_large_intestine * (1 - fr) * c_memintfr2_large_intestine_igg - kdegab * c_memintfr2_large_intestine_igg * vism_large_intestine
    d/dt(memintfr2_pancreas_igg) <- kon7b * c_memintfr1_pancreas_igg * c_fcrnmemint_pancreas * vism_pancreas - koff7b * c_memintfr2_pancreas_igg * vism_pancreas + clup_pancreas * (1 - fr) * c_endorecycfr2_pancreas_igg - clup_pancreas * (1 - fr) * c_memintfr2_pancreas_igg - kdegab * c_memintfr2_pancreas_igg * vism_pancreas
    d/dt(memintfr2_thymus_igg) <- kon7b * c_memintfr1_thymus_igg * c_fcrnmemint_thymus * vism_thymus - koff7b * c_memintfr2_thymus_igg * vism_thymus + clup_thymus * (1 - fr) * c_endorecycfr2_thymus_igg - clup_thymus * (1 - fr) * c_memintfr2_thymus_igg - kdegab * c_memintfr2_thymus_igg * vism_thymus
    d/dt(memintfr2_spleen_igg) <- kon7b * c_memintfr1_spleen_igg * c_fcrnmemint_spleen * vism_spleen - koff7b * c_memintfr2_spleen_igg * vism_spleen + clup_spleen * (1 - fr) * c_endorecycfr2_spleen_igg - clup_spleen * (1 - fr) * c_memintfr2_spleen_igg - kdegab * c_memintfr2_spleen_igg * vism_spleen
    d/dt(memintfr2_other_igg) <- kon7b * c_memintfr1_other_igg * c_fcrnmemint_other * vism_other - koff7b * c_memintfr2_other_igg * vism_other + clup_other * (1 - fr) * c_endorecycfr2_other_igg - clup_other * (1 - fr) * c_memintfr2_other_igg - kdegab * c_memintfr2_other_igg * vism_other

    # ---- Free FcRn ----
    # FcRn released by degradation of a bound complex is recovered (one
    # receptor per 1:1 complex, two per 2:1 complex). Only the two membrane
    # pools see kdeg_FcRn_Ab; the endosomal pools exchange by transit only.
    d/dt(fcrnmemvas_lung) <- clup_lung * fr * c_fcrnendorecyc_lung * (1 - frecycle) - clup_lung * fr * c_fcrnmemvas_lung + vvm_lung * (-kon7 * c_memvas_lung_igg * c_fcrnmemvas_lung + koff7 * c_memvasfr1_lung_igg + kdegab * c_memvasfr1_lung_igg - kon7b * c_memvasfr1_lung_igg * c_fcrnmemvas_lung + koff7b * c_memvasfr2_lung_igg + 2 * kdegab * c_memvasfr2_lung_igg - kon7 * c_memvas_lung * c_fcrnmemvas_lung + koff7 * c_memvasfr1_lung + kdegab * c_memvasfr1_lung - kon7b * c_memvasfr1_lung * c_fcrnmemvas_lung + koff7b * c_memvasfr2_lung + 2 * kdegab * c_memvasfr2_lung)
    d/dt(fcrnmemvas_liver) <- clup_liver * fr * c_fcrnendorecyc_liver * (1 - frecycle) - clup_liver * fr * c_fcrnmemvas_liver + vvm_liver * (-kon7 * c_memvas_liver_igg * c_fcrnmemvas_liver + koff7 * c_memvasfr1_liver_igg + kdegab * c_memvasfr1_liver_igg - kon7b * c_memvasfr1_liver_igg * c_fcrnmemvas_liver + koff7b * c_memvasfr2_liver_igg + 2 * kdegab * c_memvasfr2_liver_igg - kon7 * c_memvas_liver * c_fcrnmemvas_liver + koff7 * c_memvasfr1_liver + kdegab * c_memvasfr1_liver - kon7b * c_memvasfr1_liver * c_fcrnmemvas_liver + koff7b * c_memvasfr2_liver + 2 * kdegab * c_memvasfr2_liver)
    d/dt(fcrnmemvas_heart) <- clup_heart * fr * c_fcrnendorecyc_heart * (1 - frecycle) - clup_heart * fr * c_fcrnmemvas_heart + vvm_heart * (-kon7 * c_memvas_heart_igg * c_fcrnmemvas_heart + koff7 * c_memvasfr1_heart_igg + kdegab * c_memvasfr1_heart_igg - kon7b * c_memvasfr1_heart_igg * c_fcrnmemvas_heart + koff7b * c_memvasfr2_heart_igg + 2 * kdegab * c_memvasfr2_heart_igg - kon7 * c_memvas_heart * c_fcrnmemvas_heart + koff7 * c_memvasfr1_heart + kdegab * c_memvasfr1_heart - kon7b * c_memvasfr1_heart * c_fcrnmemvas_heart + koff7b * c_memvasfr2_heart + 2 * kdegab * c_memvasfr2_heart)
    d/dt(fcrnmemvas_muscle) <- clup_muscle * fr * c_fcrnendorecyc_muscle * (1 - frecycle) - clup_muscle * fr * c_fcrnmemvas_muscle + vvm_muscle * (-kon7 * c_memvas_muscle_igg * c_fcrnmemvas_muscle + koff7 * c_memvasfr1_muscle_igg + kdegab * c_memvasfr1_muscle_igg - kon7b * c_memvasfr1_muscle_igg * c_fcrnmemvas_muscle + koff7b * c_memvasfr2_muscle_igg + 2 * kdegab * c_memvasfr2_muscle_igg - kon7 * c_memvas_muscle * c_fcrnmemvas_muscle + koff7 * c_memvasfr1_muscle + kdegab * c_memvasfr1_muscle - kon7b * c_memvasfr1_muscle * c_fcrnmemvas_muscle + koff7b * c_memvasfr2_muscle + 2 * kdegab * c_memvasfr2_muscle)
    d/dt(fcrnmemvas_skin) <- clup_skin * fr * c_fcrnendorecyc_skin * (1 - frecycle) - clup_skin * fr * c_fcrnmemvas_skin + vvm_skin * (-kon7 * c_memvas_skin_igg * c_fcrnmemvas_skin + koff7 * c_memvasfr1_skin_igg + kdegab * c_memvasfr1_skin_igg - kon7b * c_memvasfr1_skin_igg * c_fcrnmemvas_skin + koff7b * c_memvasfr2_skin_igg + 2 * kdegab * c_memvasfr2_skin_igg - kon7 * c_memvas_skin * c_fcrnmemvas_skin + koff7 * c_memvasfr1_skin + kdegab * c_memvasfr1_skin - kon7b * c_memvasfr1_skin * c_fcrnmemvas_skin + koff7b * c_memvasfr2_skin + 2 * kdegab * c_memvasfr2_skin)
    d/dt(fcrnmemvas_adipose) <- clup_adipose * fr * c_fcrnendorecyc_adipose * (1 - frecycle) - clup_adipose * fr * c_fcrnmemvas_adipose + vvm_adipose * (-kon7 * c_memvas_adipose_igg * c_fcrnmemvas_adipose + koff7 * c_memvasfr1_adipose_igg + kdegab * c_memvasfr1_adipose_igg - kon7b * c_memvasfr1_adipose_igg * c_fcrnmemvas_adipose + koff7b * c_memvasfr2_adipose_igg + 2 * kdegab * c_memvasfr2_adipose_igg - kon7 * c_memvas_adipose * c_fcrnmemvas_adipose + koff7 * c_memvasfr1_adipose + kdegab * c_memvasfr1_adipose - kon7b * c_memvasfr1_adipose * c_fcrnmemvas_adipose + koff7b * c_memvasfr2_adipose + 2 * kdegab * c_memvasfr2_adipose)
    d/dt(fcrnmemvas_bone) <- clup_bone * fr * c_fcrnendorecyc_bone * (1 - frecycle) - clup_bone * fr * c_fcrnmemvas_bone + vvm_bone * (-kon7 * c_memvas_bone_igg * c_fcrnmemvas_bone + koff7 * c_memvasfr1_bone_igg + kdegab * c_memvasfr1_bone_igg - kon7b * c_memvasfr1_bone_igg * c_fcrnmemvas_bone + koff7b * c_memvasfr2_bone_igg + 2 * kdegab * c_memvasfr2_bone_igg - kon7 * c_memvas_bone * c_fcrnmemvas_bone + koff7 * c_memvasfr1_bone + kdegab * c_memvasfr1_bone - kon7b * c_memvasfr1_bone * c_fcrnmemvas_bone + koff7b * c_memvasfr2_bone + 2 * kdegab * c_memvasfr2_bone)
    d/dt(fcrnmemvas_brain) <- clup_brain * fr * c_fcrnendorecyc_brain * (1 - frecycle) - clup_brain * fr * c_fcrnmemvas_brain + vvm_brain * (-kon7 * c_memvas_brain_igg * c_fcrnmemvas_brain + koff7 * c_memvasfr1_brain_igg + kdegab * c_memvasfr1_brain_igg - kon7b * c_memvasfr1_brain_igg * c_fcrnmemvas_brain + koff7b * c_memvasfr2_brain_igg + 2 * kdegab * c_memvasfr2_brain_igg - kon7 * c_memvas_brain * c_fcrnmemvas_brain + koff7 * c_memvasfr1_brain + kdegab * c_memvasfr1_brain - kon7b * c_memvasfr1_brain * c_fcrnmemvas_brain + koff7b * c_memvasfr2_brain + 2 * kdegab * c_memvasfr2_brain)
    d/dt(fcrnmemvas_kidney) <- clup_kidney * fr * c_fcrnendorecyc_kidney * (1 - frecycle) - clup_kidney * fr * c_fcrnmemvas_kidney + vvm_kidney * (-kon7 * c_memvas_kidney_igg * c_fcrnmemvas_kidney + koff7 * c_memvasfr1_kidney_igg + kdegab * c_memvasfr1_kidney_igg - kon7b * c_memvasfr1_kidney_igg * c_fcrnmemvas_kidney + koff7b * c_memvasfr2_kidney_igg + 2 * kdegab * c_memvasfr2_kidney_igg - kon7 * c_memvas_kidney * c_fcrnmemvas_kidney + koff7 * c_memvasfr1_kidney + kdegab * c_memvasfr1_kidney - kon7b * c_memvasfr1_kidney * c_fcrnmemvas_kidney + koff7b * c_memvasfr2_kidney + 2 * kdegab * c_memvasfr2_kidney)
    d/dt(fcrnmemvas_small_intestine) <- clup_small_intestine * fr * c_fcrnendorecyc_small_intestine * (1 - frecycle) - clup_small_intestine * fr * c_fcrnmemvas_small_intestine + vvm_small_intestine * (-kon7 * c_memvas_small_intestine_igg * c_fcrnmemvas_small_intestine + koff7 * c_memvasfr1_small_intestine_igg + kdegab * c_memvasfr1_small_intestine_igg - kon7b * c_memvasfr1_small_intestine_igg * c_fcrnmemvas_small_intestine + koff7b * c_memvasfr2_small_intestine_igg + 2 * kdegab * c_memvasfr2_small_intestine_igg - kon7 * c_memvas_small_intestine * c_fcrnmemvas_small_intestine + koff7 * c_memvasfr1_small_intestine + kdegab * c_memvasfr1_small_intestine - kon7b * c_memvasfr1_small_intestine * c_fcrnmemvas_small_intestine + koff7b * c_memvasfr2_small_intestine + 2 * kdegab * c_memvasfr2_small_intestine)
    d/dt(fcrnmemvas_large_intestine) <- clup_large_intestine * fr * c_fcrnendorecyc_large_intestine * (1 - frecycle) - clup_large_intestine * fr * c_fcrnmemvas_large_intestine + vvm_large_intestine * (-kon7 * c_memvas_large_intestine_igg * c_fcrnmemvas_large_intestine + koff7 * c_memvasfr1_large_intestine_igg + kdegab * c_memvasfr1_large_intestine_igg - kon7b * c_memvasfr1_large_intestine_igg * c_fcrnmemvas_large_intestine + koff7b * c_memvasfr2_large_intestine_igg + 2 * kdegab * c_memvasfr2_large_intestine_igg - kon7 * c_memvas_large_intestine * c_fcrnmemvas_large_intestine + koff7 * c_memvasfr1_large_intestine + kdegab * c_memvasfr1_large_intestine - kon7b * c_memvasfr1_large_intestine * c_fcrnmemvas_large_intestine + koff7b * c_memvasfr2_large_intestine + 2 * kdegab * c_memvasfr2_large_intestine)
    d/dt(fcrnmemvas_pancreas) <- clup_pancreas * fr * c_fcrnendorecyc_pancreas * (1 - frecycle) - clup_pancreas * fr * c_fcrnmemvas_pancreas + vvm_pancreas * (-kon7 * c_memvas_pancreas_igg * c_fcrnmemvas_pancreas + koff7 * c_memvasfr1_pancreas_igg + kdegab * c_memvasfr1_pancreas_igg - kon7b * c_memvasfr1_pancreas_igg * c_fcrnmemvas_pancreas + koff7b * c_memvasfr2_pancreas_igg + 2 * kdegab * c_memvasfr2_pancreas_igg - kon7 * c_memvas_pancreas * c_fcrnmemvas_pancreas + koff7 * c_memvasfr1_pancreas + kdegab * c_memvasfr1_pancreas - kon7b * c_memvasfr1_pancreas * c_fcrnmemvas_pancreas + koff7b * c_memvasfr2_pancreas + 2 * kdegab * c_memvasfr2_pancreas)
    d/dt(fcrnmemvas_thymus) <- clup_thymus * fr * c_fcrnendorecyc_thymus * (1 - frecycle) - clup_thymus * fr * c_fcrnmemvas_thymus + vvm_thymus * (-kon7 * c_memvas_thymus_igg * c_fcrnmemvas_thymus + koff7 * c_memvasfr1_thymus_igg + kdegab * c_memvasfr1_thymus_igg - kon7b * c_memvasfr1_thymus_igg * c_fcrnmemvas_thymus + koff7b * c_memvasfr2_thymus_igg + 2 * kdegab * c_memvasfr2_thymus_igg - kon7 * c_memvas_thymus * c_fcrnmemvas_thymus + koff7 * c_memvasfr1_thymus + kdegab * c_memvasfr1_thymus - kon7b * c_memvasfr1_thymus * c_fcrnmemvas_thymus + koff7b * c_memvasfr2_thymus + 2 * kdegab * c_memvasfr2_thymus)
    d/dt(fcrnmemvas_spleen) <- clup_spleen * fr * c_fcrnendorecyc_spleen * (1 - frecycle) - clup_spleen * fr * c_fcrnmemvas_spleen + vvm_spleen * (-kon7 * c_memvas_spleen_igg * c_fcrnmemvas_spleen + koff7 * c_memvasfr1_spleen_igg + kdegab * c_memvasfr1_spleen_igg - kon7b * c_memvasfr1_spleen_igg * c_fcrnmemvas_spleen + koff7b * c_memvasfr2_spleen_igg + 2 * kdegab * c_memvasfr2_spleen_igg - kon7 * c_memvas_spleen * c_fcrnmemvas_spleen + koff7 * c_memvasfr1_spleen + kdegab * c_memvasfr1_spleen - kon7b * c_memvasfr1_spleen * c_fcrnmemvas_spleen + koff7b * c_memvasfr2_spleen + 2 * kdegab * c_memvasfr2_spleen)
    d/dt(fcrnmemvas_other) <- clup_other * fr * c_fcrnendorecyc_other * (1 - frecycle) - clup_other * fr * c_fcrnmemvas_other + vvm_other * (-kon7 * c_memvas_other_igg * c_fcrnmemvas_other + koff7 * c_memvasfr1_other_igg + kdegab * c_memvasfr1_other_igg - kon7b * c_memvasfr1_other_igg * c_fcrnmemvas_other + koff7b * c_memvasfr2_other_igg + 2 * kdegab * c_memvasfr2_other_igg - kon7 * c_memvas_other * c_fcrnmemvas_other + koff7 * c_memvasfr1_other + kdegab * c_memvasfr1_other - kon7b * c_memvasfr1_other * c_fcrnmemvas_other + koff7b * c_memvasfr2_other + 2 * kdegab * c_memvasfr2_other)
    d/dt(fcrnendoearly_lung) <- clup_lung * fr * c_fcrnmemvas_lung + clup_lung * (1 - fr) * c_fcrnmemint_lung - clup_lung * c_fcrnendoearly_lung + ve7_lung * (-kon7 * c_endoearly_lung_igg * c_fcrnendoearly_lung + koff7 * c_endoearlyfr1_lung_igg - kon7b * c_endoearlyfr1_lung_igg * c_fcrnendoearly_lung + koff7b * c_endoearlyfr2_lung_igg - kon7 * c_endoearly_lung * c_fcrnendoearly_lung + koff7 * c_endoearlyfr1_lung - kon7b * c_endoearlyfr1_lung * c_fcrnendoearly_lung + koff7b * c_endoearlyfr2_lung)
    d/dt(fcrnendoearly_liver) <- clup_liver * fr * c_fcrnmemvas_liver + clup_liver * (1 - fr) * c_fcrnmemint_liver - clup_liver * c_fcrnendoearly_liver + ve7_liver * (-kon7 * c_endoearly_liver_igg * c_fcrnendoearly_liver + koff7 * c_endoearlyfr1_liver_igg - kon7b * c_endoearlyfr1_liver_igg * c_fcrnendoearly_liver + koff7b * c_endoearlyfr2_liver_igg - kon7 * c_endoearly_liver * c_fcrnendoearly_liver + koff7 * c_endoearlyfr1_liver - kon7b * c_endoearlyfr1_liver * c_fcrnendoearly_liver + koff7b * c_endoearlyfr2_liver)
    d/dt(fcrnendoearly_heart) <- clup_heart * fr * c_fcrnmemvas_heart + clup_heart * (1 - fr) * c_fcrnmemint_heart - clup_heart * c_fcrnendoearly_heart + ve7_heart * (-kon7 * c_endoearly_heart_igg * c_fcrnendoearly_heart + koff7 * c_endoearlyfr1_heart_igg - kon7b * c_endoearlyfr1_heart_igg * c_fcrnendoearly_heart + koff7b * c_endoearlyfr2_heart_igg - kon7 * c_endoearly_heart * c_fcrnendoearly_heart + koff7 * c_endoearlyfr1_heart - kon7b * c_endoearlyfr1_heart * c_fcrnendoearly_heart + koff7b * c_endoearlyfr2_heart)
    d/dt(fcrnendoearly_muscle) <- clup_muscle * fr * c_fcrnmemvas_muscle + clup_muscle * (1 - fr) * c_fcrnmemint_muscle - clup_muscle * c_fcrnendoearly_muscle + ve7_muscle * (-kon7 * c_endoearly_muscle_igg * c_fcrnendoearly_muscle + koff7 * c_endoearlyfr1_muscle_igg - kon7b * c_endoearlyfr1_muscle_igg * c_fcrnendoearly_muscle + koff7b * c_endoearlyfr2_muscle_igg - kon7 * c_endoearly_muscle * c_fcrnendoearly_muscle + koff7 * c_endoearlyfr1_muscle - kon7b * c_endoearlyfr1_muscle * c_fcrnendoearly_muscle + koff7b * c_endoearlyfr2_muscle)
    d/dt(fcrnendoearly_skin) <- clup_skin * fr * c_fcrnmemvas_skin + clup_skin * (1 - fr) * c_fcrnmemint_skin - clup_skin * c_fcrnendoearly_skin + ve7_skin * (-kon7 * c_endoearly_skin_igg * c_fcrnendoearly_skin + koff7 * c_endoearlyfr1_skin_igg - kon7b * c_endoearlyfr1_skin_igg * c_fcrnendoearly_skin + koff7b * c_endoearlyfr2_skin_igg - kon7 * c_endoearly_skin * c_fcrnendoearly_skin + koff7 * c_endoearlyfr1_skin - kon7b * c_endoearlyfr1_skin * c_fcrnendoearly_skin + koff7b * c_endoearlyfr2_skin)
    d/dt(fcrnendoearly_adipose) <- clup_adipose * fr * c_fcrnmemvas_adipose + clup_adipose * (1 - fr) * c_fcrnmemint_adipose - clup_adipose * c_fcrnendoearly_adipose + ve7_adipose * (-kon7 * c_endoearly_adipose_igg * c_fcrnendoearly_adipose + koff7 * c_endoearlyfr1_adipose_igg - kon7b * c_endoearlyfr1_adipose_igg * c_fcrnendoearly_adipose + koff7b * c_endoearlyfr2_adipose_igg - kon7 * c_endoearly_adipose * c_fcrnendoearly_adipose + koff7 * c_endoearlyfr1_adipose - kon7b * c_endoearlyfr1_adipose * c_fcrnendoearly_adipose + koff7b * c_endoearlyfr2_adipose)
    d/dt(fcrnendoearly_bone) <- clup_bone * fr * c_fcrnmemvas_bone + clup_bone * (1 - fr) * c_fcrnmemint_bone - clup_bone * c_fcrnendoearly_bone + ve7_bone * (-kon7 * c_endoearly_bone_igg * c_fcrnendoearly_bone + koff7 * c_endoearlyfr1_bone_igg - kon7b * c_endoearlyfr1_bone_igg * c_fcrnendoearly_bone + koff7b * c_endoearlyfr2_bone_igg - kon7 * c_endoearly_bone * c_fcrnendoearly_bone + koff7 * c_endoearlyfr1_bone - kon7b * c_endoearlyfr1_bone * c_fcrnendoearly_bone + koff7b * c_endoearlyfr2_bone)
    d/dt(fcrnendoearly_brain) <- clup_brain * fr * c_fcrnmemvas_brain + clup_brain * (1 - fr) * c_fcrnmemint_brain - clup_brain * c_fcrnendoearly_brain + ve7_brain * (-kon7 * c_endoearly_brain_igg * c_fcrnendoearly_brain + koff7 * c_endoearlyfr1_brain_igg - kon7b * c_endoearlyfr1_brain_igg * c_fcrnendoearly_brain + koff7b * c_endoearlyfr2_brain_igg - kon7 * c_endoearly_brain * c_fcrnendoearly_brain + koff7 * c_endoearlyfr1_brain - kon7b * c_endoearlyfr1_brain * c_fcrnendoearly_brain + koff7b * c_endoearlyfr2_brain)
    d/dt(fcrnendoearly_kidney) <- clup_kidney * fr * c_fcrnmemvas_kidney + clup_kidney * (1 - fr) * c_fcrnmemint_kidney - clup_kidney * c_fcrnendoearly_kidney + ve7_kidney * (-kon7 * c_endoearly_kidney_igg * c_fcrnendoearly_kidney + koff7 * c_endoearlyfr1_kidney_igg - kon7b * c_endoearlyfr1_kidney_igg * c_fcrnendoearly_kidney + koff7b * c_endoearlyfr2_kidney_igg - kon7 * c_endoearly_kidney * c_fcrnendoearly_kidney + koff7 * c_endoearlyfr1_kidney - kon7b * c_endoearlyfr1_kidney * c_fcrnendoearly_kidney + koff7b * c_endoearlyfr2_kidney)
    d/dt(fcrnendoearly_small_intestine) <- clup_small_intestine * fr * c_fcrnmemvas_small_intestine + clup_small_intestine * (1 - fr) * c_fcrnmemint_small_intestine - clup_small_intestine * c_fcrnendoearly_small_intestine + ve7_small_intestine * (-kon7 * c_endoearly_small_intestine_igg * c_fcrnendoearly_small_intestine + koff7 * c_endoearlyfr1_small_intestine_igg - kon7b * c_endoearlyfr1_small_intestine_igg * c_fcrnendoearly_small_intestine + koff7b * c_endoearlyfr2_small_intestine_igg - kon7 * c_endoearly_small_intestine * c_fcrnendoearly_small_intestine + koff7 * c_endoearlyfr1_small_intestine - kon7b * c_endoearlyfr1_small_intestine * c_fcrnendoearly_small_intestine + koff7b * c_endoearlyfr2_small_intestine)
    d/dt(fcrnendoearly_large_intestine) <- clup_large_intestine * fr * c_fcrnmemvas_large_intestine + clup_large_intestine * (1 - fr) * c_fcrnmemint_large_intestine - clup_large_intestine * c_fcrnendoearly_large_intestine + ve7_large_intestine * (-kon7 * c_endoearly_large_intestine_igg * c_fcrnendoearly_large_intestine + koff7 * c_endoearlyfr1_large_intestine_igg - kon7b * c_endoearlyfr1_large_intestine_igg * c_fcrnendoearly_large_intestine + koff7b * c_endoearlyfr2_large_intestine_igg - kon7 * c_endoearly_large_intestine * c_fcrnendoearly_large_intestine + koff7 * c_endoearlyfr1_large_intestine - kon7b * c_endoearlyfr1_large_intestine * c_fcrnendoearly_large_intestine + koff7b * c_endoearlyfr2_large_intestine)
    d/dt(fcrnendoearly_pancreas) <- clup_pancreas * fr * c_fcrnmemvas_pancreas + clup_pancreas * (1 - fr) * c_fcrnmemint_pancreas - clup_pancreas * c_fcrnendoearly_pancreas + ve7_pancreas * (-kon7 * c_endoearly_pancreas_igg * c_fcrnendoearly_pancreas + koff7 * c_endoearlyfr1_pancreas_igg - kon7b * c_endoearlyfr1_pancreas_igg * c_fcrnendoearly_pancreas + koff7b * c_endoearlyfr2_pancreas_igg - kon7 * c_endoearly_pancreas * c_fcrnendoearly_pancreas + koff7 * c_endoearlyfr1_pancreas - kon7b * c_endoearlyfr1_pancreas * c_fcrnendoearly_pancreas + koff7b * c_endoearlyfr2_pancreas)
    d/dt(fcrnendoearly_thymus) <- clup_thymus * fr * c_fcrnmemvas_thymus + clup_thymus * (1 - fr) * c_fcrnmemint_thymus - clup_thymus * c_fcrnendoearly_thymus + ve7_thymus * (-kon7 * c_endoearly_thymus_igg * c_fcrnendoearly_thymus + koff7 * c_endoearlyfr1_thymus_igg - kon7b * c_endoearlyfr1_thymus_igg * c_fcrnendoearly_thymus + koff7b * c_endoearlyfr2_thymus_igg - kon7 * c_endoearly_thymus * c_fcrnendoearly_thymus + koff7 * c_endoearlyfr1_thymus - kon7b * c_endoearlyfr1_thymus * c_fcrnendoearly_thymus + koff7b * c_endoearlyfr2_thymus)
    d/dt(fcrnendoearly_spleen) <- clup_spleen * fr * c_fcrnmemvas_spleen + clup_spleen * (1 - fr) * c_fcrnmemint_spleen - clup_spleen * c_fcrnendoearly_spleen + ve7_spleen * (-kon7 * c_endoearly_spleen_igg * c_fcrnendoearly_spleen + koff7 * c_endoearlyfr1_spleen_igg - kon7b * c_endoearlyfr1_spleen_igg * c_fcrnendoearly_spleen + koff7b * c_endoearlyfr2_spleen_igg - kon7 * c_endoearly_spleen * c_fcrnendoearly_spleen + koff7 * c_endoearlyfr1_spleen - kon7b * c_endoearlyfr1_spleen * c_fcrnendoearly_spleen + koff7b * c_endoearlyfr2_spleen)
    d/dt(fcrnendoearly_other) <- clup_other * fr * c_fcrnmemvas_other + clup_other * (1 - fr) * c_fcrnmemint_other - clup_other * c_fcrnendoearly_other + ve7_other * (-kon7 * c_endoearly_other_igg * c_fcrnendoearly_other + koff7 * c_endoearlyfr1_other_igg - kon7b * c_endoearlyfr1_other_igg * c_fcrnendoearly_other + koff7b * c_endoearlyfr2_other_igg - kon7 * c_endoearly_other * c_fcrnendoearly_other + koff7 * c_endoearlyfr1_other - kon7b * c_endoearlyfr1_other * c_fcrnendoearly_other + koff7b * c_endoearlyfr2_other)
    d/dt(fcrnendosort_lung) <- clup_lung * c_fcrnendoearly_lung - clup_lung * c_fcrnendosort_lung + clup_lung * c_fcrnendorecyc_lung * frecycle + ve6a_lung * (-kon6 * c_endosort_lung_igg * c_fcrnendosort_lung + koff6 * c_endosortfr1_lung_igg - kon6b * c_endosortfr1_lung_igg * c_fcrnendosort_lung + koff6b * c_endosortfr2_lung_igg - kon6 * c_endosort_lung * c_fcrnendosort_lung + koff6 * c_endosortfr1_lung - kon6b * c_endosortfr1_lung * c_fcrnendosort_lung + koff6b * c_endosortfr2_lung)
    d/dt(fcrnendosort_liver) <- clup_liver * c_fcrnendoearly_liver - clup_liver * c_fcrnendosort_liver + clup_liver * c_fcrnendorecyc_liver * frecycle + ve6a_liver * (-kon6 * c_endosort_liver_igg * c_fcrnendosort_liver + koff6 * c_endosortfr1_liver_igg - kon6b * c_endosortfr1_liver_igg * c_fcrnendosort_liver + koff6b * c_endosortfr2_liver_igg - kon6 * c_endosort_liver * c_fcrnendosort_liver + koff6 * c_endosortfr1_liver - kon6b * c_endosortfr1_liver * c_fcrnendosort_liver + koff6b * c_endosortfr2_liver)
    d/dt(fcrnendosort_heart) <- clup_heart * c_fcrnendoearly_heart - clup_heart * c_fcrnendosort_heart + clup_heart * c_fcrnendorecyc_heart * frecycle + ve6a_heart * (-kon6 * c_endosort_heart_igg * c_fcrnendosort_heart + koff6 * c_endosortfr1_heart_igg - kon6b * c_endosortfr1_heart_igg * c_fcrnendosort_heart + koff6b * c_endosortfr2_heart_igg - kon6 * c_endosort_heart * c_fcrnendosort_heart + koff6 * c_endosortfr1_heart - kon6b * c_endosortfr1_heart * c_fcrnendosort_heart + koff6b * c_endosortfr2_heart)
    d/dt(fcrnendosort_muscle) <- clup_muscle * c_fcrnendoearly_muscle - clup_muscle * c_fcrnendosort_muscle + clup_muscle * c_fcrnendorecyc_muscle * frecycle + ve6a_muscle * (-kon6 * c_endosort_muscle_igg * c_fcrnendosort_muscle + koff6 * c_endosortfr1_muscle_igg - kon6b * c_endosortfr1_muscle_igg * c_fcrnendosort_muscle + koff6b * c_endosortfr2_muscle_igg - kon6 * c_endosort_muscle * c_fcrnendosort_muscle + koff6 * c_endosortfr1_muscle - kon6b * c_endosortfr1_muscle * c_fcrnendosort_muscle + koff6b * c_endosortfr2_muscle)
    d/dt(fcrnendosort_skin) <- clup_skin * c_fcrnendoearly_skin - clup_skin * c_fcrnendosort_skin + clup_skin * c_fcrnendorecyc_skin * frecycle + ve6a_skin * (-kon6 * c_endosort_skin_igg * c_fcrnendosort_skin + koff6 * c_endosortfr1_skin_igg - kon6b * c_endosortfr1_skin_igg * c_fcrnendosort_skin + koff6b * c_endosortfr2_skin_igg - kon6 * c_endosort_skin * c_fcrnendosort_skin + koff6 * c_endosortfr1_skin - kon6b * c_endosortfr1_skin * c_fcrnendosort_skin + koff6b * c_endosortfr2_skin)
    d/dt(fcrnendosort_adipose) <- clup_adipose * c_fcrnendoearly_adipose - clup_adipose * c_fcrnendosort_adipose + clup_adipose * c_fcrnendorecyc_adipose * frecycle + ve6a_adipose * (-kon6 * c_endosort_adipose_igg * c_fcrnendosort_adipose + koff6 * c_endosortfr1_adipose_igg - kon6b * c_endosortfr1_adipose_igg * c_fcrnendosort_adipose + koff6b * c_endosortfr2_adipose_igg - kon6 * c_endosort_adipose * c_fcrnendosort_adipose + koff6 * c_endosortfr1_adipose - kon6b * c_endosortfr1_adipose * c_fcrnendosort_adipose + koff6b * c_endosortfr2_adipose)
    d/dt(fcrnendosort_bone) <- clup_bone * c_fcrnendoearly_bone - clup_bone * c_fcrnendosort_bone + clup_bone * c_fcrnendorecyc_bone * frecycle + ve6a_bone * (-kon6 * c_endosort_bone_igg * c_fcrnendosort_bone + koff6 * c_endosortfr1_bone_igg - kon6b * c_endosortfr1_bone_igg * c_fcrnendosort_bone + koff6b * c_endosortfr2_bone_igg - kon6 * c_endosort_bone * c_fcrnendosort_bone + koff6 * c_endosortfr1_bone - kon6b * c_endosortfr1_bone * c_fcrnendosort_bone + koff6b * c_endosortfr2_bone)
    d/dt(fcrnendosort_brain) <- clup_brain * c_fcrnendoearly_brain - clup_brain * c_fcrnendosort_brain + clup_brain * c_fcrnendorecyc_brain * frecycle + ve6a_brain * (-kon6 * c_endosort_brain_igg * c_fcrnendosort_brain + koff6 * c_endosortfr1_brain_igg - kon6b * c_endosortfr1_brain_igg * c_fcrnendosort_brain + koff6b * c_endosortfr2_brain_igg - kon6 * c_endosort_brain * c_fcrnendosort_brain + koff6 * c_endosortfr1_brain - kon6b * c_endosortfr1_brain * c_fcrnendosort_brain + koff6b * c_endosortfr2_brain)
    d/dt(fcrnendosort_kidney) <- clup_kidney * c_fcrnendoearly_kidney - clup_kidney * c_fcrnendosort_kidney + clup_kidney * c_fcrnendorecyc_kidney * frecycle + ve6a_kidney * (-kon6 * c_endosort_kidney_igg * c_fcrnendosort_kidney + koff6 * c_endosortfr1_kidney_igg - kon6b * c_endosortfr1_kidney_igg * c_fcrnendosort_kidney + koff6b * c_endosortfr2_kidney_igg - kon6 * c_endosort_kidney * c_fcrnendosort_kidney + koff6 * c_endosortfr1_kidney - kon6b * c_endosortfr1_kidney * c_fcrnendosort_kidney + koff6b * c_endosortfr2_kidney)
    d/dt(fcrnendosort_small_intestine) <- clup_small_intestine * c_fcrnendoearly_small_intestine - clup_small_intestine * c_fcrnendosort_small_intestine + clup_small_intestine * c_fcrnendorecyc_small_intestine * frecycle + ve6a_small_intestine * (-kon6 * c_endosort_small_intestine_igg * c_fcrnendosort_small_intestine + koff6 * c_endosortfr1_small_intestine_igg - kon6b * c_endosortfr1_small_intestine_igg * c_fcrnendosort_small_intestine + koff6b * c_endosortfr2_small_intestine_igg - kon6 * c_endosort_small_intestine * c_fcrnendosort_small_intestine + koff6 * c_endosortfr1_small_intestine - kon6b * c_endosortfr1_small_intestine * c_fcrnendosort_small_intestine + koff6b * c_endosortfr2_small_intestine)
    d/dt(fcrnendosort_large_intestine) <- clup_large_intestine * c_fcrnendoearly_large_intestine - clup_large_intestine * c_fcrnendosort_large_intestine + clup_large_intestine * c_fcrnendorecyc_large_intestine * frecycle + ve6a_large_intestine * (-kon6 * c_endosort_large_intestine_igg * c_fcrnendosort_large_intestine + koff6 * c_endosortfr1_large_intestine_igg - kon6b * c_endosortfr1_large_intestine_igg * c_fcrnendosort_large_intestine + koff6b * c_endosortfr2_large_intestine_igg - kon6 * c_endosort_large_intestine * c_fcrnendosort_large_intestine + koff6 * c_endosortfr1_large_intestine - kon6b * c_endosortfr1_large_intestine * c_fcrnendosort_large_intestine + koff6b * c_endosortfr2_large_intestine)
    d/dt(fcrnendosort_pancreas) <- clup_pancreas * c_fcrnendoearly_pancreas - clup_pancreas * c_fcrnendosort_pancreas + clup_pancreas * c_fcrnendorecyc_pancreas * frecycle + ve6a_pancreas * (-kon6 * c_endosort_pancreas_igg * c_fcrnendosort_pancreas + koff6 * c_endosortfr1_pancreas_igg - kon6b * c_endosortfr1_pancreas_igg * c_fcrnendosort_pancreas + koff6b * c_endosortfr2_pancreas_igg - kon6 * c_endosort_pancreas * c_fcrnendosort_pancreas + koff6 * c_endosortfr1_pancreas - kon6b * c_endosortfr1_pancreas * c_fcrnendosort_pancreas + koff6b * c_endosortfr2_pancreas)
    d/dt(fcrnendosort_thymus) <- clup_thymus * c_fcrnendoearly_thymus - clup_thymus * c_fcrnendosort_thymus + clup_thymus * c_fcrnendorecyc_thymus * frecycle + ve6a_thymus * (-kon6 * c_endosort_thymus_igg * c_fcrnendosort_thymus + koff6 * c_endosortfr1_thymus_igg - kon6b * c_endosortfr1_thymus_igg * c_fcrnendosort_thymus + koff6b * c_endosortfr2_thymus_igg - kon6 * c_endosort_thymus * c_fcrnendosort_thymus + koff6 * c_endosortfr1_thymus - kon6b * c_endosortfr1_thymus * c_fcrnendosort_thymus + koff6b * c_endosortfr2_thymus)
    d/dt(fcrnendosort_spleen) <- clup_spleen * c_fcrnendoearly_spleen - clup_spleen * c_fcrnendosort_spleen + clup_spleen * c_fcrnendorecyc_spleen * frecycle + ve6a_spleen * (-kon6 * c_endosort_spleen_igg * c_fcrnendosort_spleen + koff6 * c_endosortfr1_spleen_igg - kon6b * c_endosortfr1_spleen_igg * c_fcrnendosort_spleen + koff6b * c_endosortfr2_spleen_igg - kon6 * c_endosort_spleen * c_fcrnendosort_spleen + koff6 * c_endosortfr1_spleen - kon6b * c_endosortfr1_spleen * c_fcrnendosort_spleen + koff6b * c_endosortfr2_spleen)
    d/dt(fcrnendosort_other) <- clup_other * c_fcrnendoearly_other - clup_other * c_fcrnendosort_other + clup_other * c_fcrnendorecyc_other * frecycle + ve6a_other * (-kon6 * c_endosort_other_igg * c_fcrnendosort_other + koff6 * c_endosortfr1_other_igg - kon6b * c_endosortfr1_other_igg * c_fcrnendosort_other + koff6b * c_endosortfr2_other_igg - kon6 * c_endosort_other * c_fcrnendosort_other + koff6 * c_endosortfr1_other - kon6b * c_endosortfr1_other * c_fcrnendosort_other + koff6b * c_endosortfr2_other)
    d/dt(fcrnendorecyc_lung) <- clup_lung * c_fcrnendosort_lung - clup_lung * c_fcrnendorecyc_lung + ve7b_lung * (-kon7 * c_endorecyc_lung_igg * c_fcrnendorecyc_lung + koff7 * c_endorecycfr1_lung_igg - kon7b * c_endorecycfr1_lung_igg * c_fcrnendorecyc_lung + koff7b * c_endorecycfr2_lung_igg - kon7 * c_endorecyc_lung * c_fcrnendorecyc_lung + koff7 * c_endorecycfr1_lung - kon7b * c_endorecycfr1_lung * c_fcrnendorecyc_lung + koff7b * c_endorecycfr2_lung)
    d/dt(fcrnendorecyc_liver) <- clup_liver * c_fcrnendosort_liver - clup_liver * c_fcrnendorecyc_liver + ve7b_liver * (-kon7 * c_endorecyc_liver_igg * c_fcrnendorecyc_liver + koff7 * c_endorecycfr1_liver_igg - kon7b * c_endorecycfr1_liver_igg * c_fcrnendorecyc_liver + koff7b * c_endorecycfr2_liver_igg - kon7 * c_endorecyc_liver * c_fcrnendorecyc_liver + koff7 * c_endorecycfr1_liver - kon7b * c_endorecycfr1_liver * c_fcrnendorecyc_liver + koff7b * c_endorecycfr2_liver)
    d/dt(fcrnendorecyc_heart) <- clup_heart * c_fcrnendosort_heart - clup_heart * c_fcrnendorecyc_heart + ve7b_heart * (-kon7 * c_endorecyc_heart_igg * c_fcrnendorecyc_heart + koff7 * c_endorecycfr1_heart_igg - kon7b * c_endorecycfr1_heart_igg * c_fcrnendorecyc_heart + koff7b * c_endorecycfr2_heart_igg - kon7 * c_endorecyc_heart * c_fcrnendorecyc_heart + koff7 * c_endorecycfr1_heart - kon7b * c_endorecycfr1_heart * c_fcrnendorecyc_heart + koff7b * c_endorecycfr2_heart)
    d/dt(fcrnendorecyc_muscle) <- clup_muscle * c_fcrnendosort_muscle - clup_muscle * c_fcrnendorecyc_muscle + ve7b_muscle * (-kon7 * c_endorecyc_muscle_igg * c_fcrnendorecyc_muscle + koff7 * c_endorecycfr1_muscle_igg - kon7b * c_endorecycfr1_muscle_igg * c_fcrnendorecyc_muscle + koff7b * c_endorecycfr2_muscle_igg - kon7 * c_endorecyc_muscle * c_fcrnendorecyc_muscle + koff7 * c_endorecycfr1_muscle - kon7b * c_endorecycfr1_muscle * c_fcrnendorecyc_muscle + koff7b * c_endorecycfr2_muscle)
    d/dt(fcrnendorecyc_skin) <- clup_skin * c_fcrnendosort_skin - clup_skin * c_fcrnendorecyc_skin + ve7b_skin * (-kon7 * c_endorecyc_skin_igg * c_fcrnendorecyc_skin + koff7 * c_endorecycfr1_skin_igg - kon7b * c_endorecycfr1_skin_igg * c_fcrnendorecyc_skin + koff7b * c_endorecycfr2_skin_igg - kon7 * c_endorecyc_skin * c_fcrnendorecyc_skin + koff7 * c_endorecycfr1_skin - kon7b * c_endorecycfr1_skin * c_fcrnendorecyc_skin + koff7b * c_endorecycfr2_skin)
    d/dt(fcrnendorecyc_adipose) <- clup_adipose * c_fcrnendosort_adipose - clup_adipose * c_fcrnendorecyc_adipose + ve7b_adipose * (-kon7 * c_endorecyc_adipose_igg * c_fcrnendorecyc_adipose + koff7 * c_endorecycfr1_adipose_igg - kon7b * c_endorecycfr1_adipose_igg * c_fcrnendorecyc_adipose + koff7b * c_endorecycfr2_adipose_igg - kon7 * c_endorecyc_adipose * c_fcrnendorecyc_adipose + koff7 * c_endorecycfr1_adipose - kon7b * c_endorecycfr1_adipose * c_fcrnendorecyc_adipose + koff7b * c_endorecycfr2_adipose)
    d/dt(fcrnendorecyc_bone) <- clup_bone * c_fcrnendosort_bone - clup_bone * c_fcrnendorecyc_bone + ve7b_bone * (-kon7 * c_endorecyc_bone_igg * c_fcrnendorecyc_bone + koff7 * c_endorecycfr1_bone_igg - kon7b * c_endorecycfr1_bone_igg * c_fcrnendorecyc_bone + koff7b * c_endorecycfr2_bone_igg - kon7 * c_endorecyc_bone * c_fcrnendorecyc_bone + koff7 * c_endorecycfr1_bone - kon7b * c_endorecycfr1_bone * c_fcrnendorecyc_bone + koff7b * c_endorecycfr2_bone)
    d/dt(fcrnendorecyc_brain) <- clup_brain * c_fcrnendosort_brain - clup_brain * c_fcrnendorecyc_brain + ve7b_brain * (-kon7 * c_endorecyc_brain_igg * c_fcrnendorecyc_brain + koff7 * c_endorecycfr1_brain_igg - kon7b * c_endorecycfr1_brain_igg * c_fcrnendorecyc_brain + koff7b * c_endorecycfr2_brain_igg - kon7 * c_endorecyc_brain * c_fcrnendorecyc_brain + koff7 * c_endorecycfr1_brain - kon7b * c_endorecycfr1_brain * c_fcrnendorecyc_brain + koff7b * c_endorecycfr2_brain)
    d/dt(fcrnendorecyc_kidney) <- clup_kidney * c_fcrnendosort_kidney - clup_kidney * c_fcrnendorecyc_kidney + ve7b_kidney * (-kon7 * c_endorecyc_kidney_igg * c_fcrnendorecyc_kidney + koff7 * c_endorecycfr1_kidney_igg - kon7b * c_endorecycfr1_kidney_igg * c_fcrnendorecyc_kidney + koff7b * c_endorecycfr2_kidney_igg - kon7 * c_endorecyc_kidney * c_fcrnendorecyc_kidney + koff7 * c_endorecycfr1_kidney - kon7b * c_endorecycfr1_kidney * c_fcrnendorecyc_kidney + koff7b * c_endorecycfr2_kidney)
    d/dt(fcrnendorecyc_small_intestine) <- clup_small_intestine * c_fcrnendosort_small_intestine - clup_small_intestine * c_fcrnendorecyc_small_intestine + ve7b_small_intestine * (-kon7 * c_endorecyc_small_intestine_igg * c_fcrnendorecyc_small_intestine + koff7 * c_endorecycfr1_small_intestine_igg - kon7b * c_endorecycfr1_small_intestine_igg * c_fcrnendorecyc_small_intestine + koff7b * c_endorecycfr2_small_intestine_igg - kon7 * c_endorecyc_small_intestine * c_fcrnendorecyc_small_intestine + koff7 * c_endorecycfr1_small_intestine - kon7b * c_endorecycfr1_small_intestine * c_fcrnendorecyc_small_intestine + koff7b * c_endorecycfr2_small_intestine)
    d/dt(fcrnendorecyc_large_intestine) <- clup_large_intestine * c_fcrnendosort_large_intestine - clup_large_intestine * c_fcrnendorecyc_large_intestine + ve7b_large_intestine * (-kon7 * c_endorecyc_large_intestine_igg * c_fcrnendorecyc_large_intestine + koff7 * c_endorecycfr1_large_intestine_igg - kon7b * c_endorecycfr1_large_intestine_igg * c_fcrnendorecyc_large_intestine + koff7b * c_endorecycfr2_large_intestine_igg - kon7 * c_endorecyc_large_intestine * c_fcrnendorecyc_large_intestine + koff7 * c_endorecycfr1_large_intestine - kon7b * c_endorecycfr1_large_intestine * c_fcrnendorecyc_large_intestine + koff7b * c_endorecycfr2_large_intestine)
    d/dt(fcrnendorecyc_pancreas) <- clup_pancreas * c_fcrnendosort_pancreas - clup_pancreas * c_fcrnendorecyc_pancreas + ve7b_pancreas * (-kon7 * c_endorecyc_pancreas_igg * c_fcrnendorecyc_pancreas + koff7 * c_endorecycfr1_pancreas_igg - kon7b * c_endorecycfr1_pancreas_igg * c_fcrnendorecyc_pancreas + koff7b * c_endorecycfr2_pancreas_igg - kon7 * c_endorecyc_pancreas * c_fcrnendorecyc_pancreas + koff7 * c_endorecycfr1_pancreas - kon7b * c_endorecycfr1_pancreas * c_fcrnendorecyc_pancreas + koff7b * c_endorecycfr2_pancreas)
    d/dt(fcrnendorecyc_thymus) <- clup_thymus * c_fcrnendosort_thymus - clup_thymus * c_fcrnendorecyc_thymus + ve7b_thymus * (-kon7 * c_endorecyc_thymus_igg * c_fcrnendorecyc_thymus + koff7 * c_endorecycfr1_thymus_igg - kon7b * c_endorecycfr1_thymus_igg * c_fcrnendorecyc_thymus + koff7b * c_endorecycfr2_thymus_igg - kon7 * c_endorecyc_thymus * c_fcrnendorecyc_thymus + koff7 * c_endorecycfr1_thymus - kon7b * c_endorecycfr1_thymus * c_fcrnendorecyc_thymus + koff7b * c_endorecycfr2_thymus)
    d/dt(fcrnendorecyc_spleen) <- clup_spleen * c_fcrnendosort_spleen - clup_spleen * c_fcrnendorecyc_spleen + ve7b_spleen * (-kon7 * c_endorecyc_spleen_igg * c_fcrnendorecyc_spleen + koff7 * c_endorecycfr1_spleen_igg - kon7b * c_endorecycfr1_spleen_igg * c_fcrnendorecyc_spleen + koff7b * c_endorecycfr2_spleen_igg - kon7 * c_endorecyc_spleen * c_fcrnendorecyc_spleen + koff7 * c_endorecycfr1_spleen - kon7b * c_endorecycfr1_spleen * c_fcrnendorecyc_spleen + koff7b * c_endorecycfr2_spleen)
    d/dt(fcrnendorecyc_other) <- clup_other * c_fcrnendosort_other - clup_other * c_fcrnendorecyc_other + ve7b_other * (-kon7 * c_endorecyc_other_igg * c_fcrnendorecyc_other + koff7 * c_endorecycfr1_other_igg - kon7b * c_endorecycfr1_other_igg * c_fcrnendorecyc_other + koff7b * c_endorecycfr2_other_igg - kon7 * c_endorecyc_other * c_fcrnendorecyc_other + koff7 * c_endorecycfr1_other - kon7b * c_endorecycfr1_other * c_fcrnendorecyc_other + koff7b * c_endorecycfr2_other)
    d/dt(fcrnmemint_lung) <- clup_lung * (1 - fr) * c_fcrnendorecyc_lung * (1 - frecycle) - clup_lung * (1 - fr) * c_fcrnmemint_lung + vism_lung * (-kon7 * c_memint_lung_igg * c_fcrnmemint_lung + koff7 * c_memintfr1_lung_igg + kdegab * c_memintfr1_lung_igg - kon7b * c_memintfr1_lung_igg * c_fcrnmemint_lung + koff7b * c_memintfr2_lung_igg + 2 * kdegab * c_memintfr2_lung_igg - kon7 * c_memint_lung * c_fcrnmemint_lung + koff7 * c_memintfr1_lung + kdegab * c_memintfr1_lung - kon7b * c_memintfr1_lung * c_fcrnmemint_lung + koff7b * c_memintfr2_lung + 2 * kdegab * c_memintfr2_lung)
    d/dt(fcrnmemint_liver) <- clup_liver * (1 - fr) * c_fcrnendorecyc_liver * (1 - frecycle) - clup_liver * (1 - fr) * c_fcrnmemint_liver + vism_liver * (-kon7 * c_memint_liver_igg * c_fcrnmemint_liver + koff7 * c_memintfr1_liver_igg + kdegab * c_memintfr1_liver_igg - kon7b * c_memintfr1_liver_igg * c_fcrnmemint_liver + koff7b * c_memintfr2_liver_igg + 2 * kdegab * c_memintfr2_liver_igg - kon7 * c_memint_liver * c_fcrnmemint_liver + koff7 * c_memintfr1_liver + kdegab * c_memintfr1_liver - kon7b * c_memintfr1_liver * c_fcrnmemint_liver + koff7b * c_memintfr2_liver + 2 * kdegab * c_memintfr2_liver)
    d/dt(fcrnmemint_heart) <- clup_heart * (1 - fr) * c_fcrnendorecyc_heart * (1 - frecycle) - clup_heart * (1 - fr) * c_fcrnmemint_heart + vism_heart * (-kon7 * c_memint_heart_igg * c_fcrnmemint_heart + koff7 * c_memintfr1_heart_igg + kdegab * c_memintfr1_heart_igg - kon7b * c_memintfr1_heart_igg * c_fcrnmemint_heart + koff7b * c_memintfr2_heart_igg + 2 * kdegab * c_memintfr2_heart_igg - kon7 * c_memint_heart * c_fcrnmemint_heart + koff7 * c_memintfr1_heart + kdegab * c_memintfr1_heart - kon7b * c_memintfr1_heart * c_fcrnmemint_heart + koff7b * c_memintfr2_heart + 2 * kdegab * c_memintfr2_heart)
    d/dt(fcrnmemint_muscle) <- clup_muscle * (1 - fr) * c_fcrnendorecyc_muscle * (1 - frecycle) - clup_muscle * (1 - fr) * c_fcrnmemint_muscle + vism_muscle * (-kon7 * c_memint_muscle_igg * c_fcrnmemint_muscle + koff7 * c_memintfr1_muscle_igg + kdegab * c_memintfr1_muscle_igg - kon7b * c_memintfr1_muscle_igg * c_fcrnmemint_muscle + koff7b * c_memintfr2_muscle_igg + 2 * kdegab * c_memintfr2_muscle_igg - kon7 * c_memint_muscle * c_fcrnmemint_muscle + koff7 * c_memintfr1_muscle + kdegab * c_memintfr1_muscle - kon7b * c_memintfr1_muscle * c_fcrnmemint_muscle + koff7b * c_memintfr2_muscle + 2 * kdegab * c_memintfr2_muscle)
    d/dt(fcrnmemint_skin) <- clup_skin * (1 - fr) * c_fcrnendorecyc_skin * (1 - frecycle) - clup_skin * (1 - fr) * c_fcrnmemint_skin + vism_skin * (-kon7 * c_memint_skin_igg * c_fcrnmemint_skin + koff7 * c_memintfr1_skin_igg + kdegab * c_memintfr1_skin_igg - kon7b * c_memintfr1_skin_igg * c_fcrnmemint_skin + koff7b * c_memintfr2_skin_igg + 2 * kdegab * c_memintfr2_skin_igg - kon7 * c_memint_skin * c_fcrnmemint_skin + koff7 * c_memintfr1_skin + kdegab * c_memintfr1_skin - kon7b * c_memintfr1_skin * c_fcrnmemint_skin + koff7b * c_memintfr2_skin + 2 * kdegab * c_memintfr2_skin)
    d/dt(fcrnmemint_adipose) <- clup_adipose * (1 - fr) * c_fcrnendorecyc_adipose * (1 - frecycle) - clup_adipose * (1 - fr) * c_fcrnmemint_adipose + vism_adipose * (-kon7 * c_memint_adipose_igg * c_fcrnmemint_adipose + koff7 * c_memintfr1_adipose_igg + kdegab * c_memintfr1_adipose_igg - kon7b * c_memintfr1_adipose_igg * c_fcrnmemint_adipose + koff7b * c_memintfr2_adipose_igg + 2 * kdegab * c_memintfr2_adipose_igg - kon7 * c_memint_adipose * c_fcrnmemint_adipose + koff7 * c_memintfr1_adipose + kdegab * c_memintfr1_adipose - kon7b * c_memintfr1_adipose * c_fcrnmemint_adipose + koff7b * c_memintfr2_adipose + 2 * kdegab * c_memintfr2_adipose)
    d/dt(fcrnmemint_bone) <- clup_bone * (1 - fr) * c_fcrnendorecyc_bone * (1 - frecycle) - clup_bone * (1 - fr) * c_fcrnmemint_bone + vism_bone * (-kon7 * c_memint_bone_igg * c_fcrnmemint_bone + koff7 * c_memintfr1_bone_igg + kdegab * c_memintfr1_bone_igg - kon7b * c_memintfr1_bone_igg * c_fcrnmemint_bone + koff7b * c_memintfr2_bone_igg + 2 * kdegab * c_memintfr2_bone_igg - kon7 * c_memint_bone * c_fcrnmemint_bone + koff7 * c_memintfr1_bone + kdegab * c_memintfr1_bone - kon7b * c_memintfr1_bone * c_fcrnmemint_bone + koff7b * c_memintfr2_bone + 2 * kdegab * c_memintfr2_bone)
    d/dt(fcrnmemint_brain) <- clup_brain * (1 - fr) * c_fcrnendorecyc_brain * (1 - frecycle) - clup_brain * (1 - fr) * c_fcrnmemint_brain + vism_brain * (-kon7 * c_memint_brain_igg * c_fcrnmemint_brain + koff7 * c_memintfr1_brain_igg + kdegab * c_memintfr1_brain_igg - kon7b * c_memintfr1_brain_igg * c_fcrnmemint_brain + koff7b * c_memintfr2_brain_igg + 2 * kdegab * c_memintfr2_brain_igg - kon7 * c_memint_brain * c_fcrnmemint_brain + koff7 * c_memintfr1_brain + kdegab * c_memintfr1_brain - kon7b * c_memintfr1_brain * c_fcrnmemint_brain + koff7b * c_memintfr2_brain + 2 * kdegab * c_memintfr2_brain)
    d/dt(fcrnmemint_kidney) <- clup_kidney * (1 - fr) * c_fcrnendorecyc_kidney * (1 - frecycle) - clup_kidney * (1 - fr) * c_fcrnmemint_kidney + vism_kidney * (-kon7 * c_memint_kidney_igg * c_fcrnmemint_kidney + koff7 * c_memintfr1_kidney_igg + kdegab * c_memintfr1_kidney_igg - kon7b * c_memintfr1_kidney_igg * c_fcrnmemint_kidney + koff7b * c_memintfr2_kidney_igg + 2 * kdegab * c_memintfr2_kidney_igg - kon7 * c_memint_kidney * c_fcrnmemint_kidney + koff7 * c_memintfr1_kidney + kdegab * c_memintfr1_kidney - kon7b * c_memintfr1_kidney * c_fcrnmemint_kidney + koff7b * c_memintfr2_kidney + 2 * kdegab * c_memintfr2_kidney)
    d/dt(fcrnmemint_small_intestine) <- clup_small_intestine * (1 - fr) * c_fcrnendorecyc_small_intestine * (1 - frecycle) - clup_small_intestine * (1 - fr) * c_fcrnmemint_small_intestine + vism_small_intestine * (-kon7 * c_memint_small_intestine_igg * c_fcrnmemint_small_intestine + koff7 * c_memintfr1_small_intestine_igg + kdegab * c_memintfr1_small_intestine_igg - kon7b * c_memintfr1_small_intestine_igg * c_fcrnmemint_small_intestine + koff7b * c_memintfr2_small_intestine_igg + 2 * kdegab * c_memintfr2_small_intestine_igg - kon7 * c_memint_small_intestine * c_fcrnmemint_small_intestine + koff7 * c_memintfr1_small_intestine + kdegab * c_memintfr1_small_intestine - kon7b * c_memintfr1_small_intestine * c_fcrnmemint_small_intestine + koff7b * c_memintfr2_small_intestine + 2 * kdegab * c_memintfr2_small_intestine)
    d/dt(fcrnmemint_large_intestine) <- clup_large_intestine * (1 - fr) * c_fcrnendorecyc_large_intestine * (1 - frecycle) - clup_large_intestine * (1 - fr) * c_fcrnmemint_large_intestine + vism_large_intestine * (-kon7 * c_memint_large_intestine_igg * c_fcrnmemint_large_intestine + koff7 * c_memintfr1_large_intestine_igg + kdegab * c_memintfr1_large_intestine_igg - kon7b * c_memintfr1_large_intestine_igg * c_fcrnmemint_large_intestine + koff7b * c_memintfr2_large_intestine_igg + 2 * kdegab * c_memintfr2_large_intestine_igg - kon7 * c_memint_large_intestine * c_fcrnmemint_large_intestine + koff7 * c_memintfr1_large_intestine + kdegab * c_memintfr1_large_intestine - kon7b * c_memintfr1_large_intestine * c_fcrnmemint_large_intestine + koff7b * c_memintfr2_large_intestine + 2 * kdegab * c_memintfr2_large_intestine)
    d/dt(fcrnmemint_pancreas) <- clup_pancreas * (1 - fr) * c_fcrnendorecyc_pancreas * (1 - frecycle) - clup_pancreas * (1 - fr) * c_fcrnmemint_pancreas + vism_pancreas * (-kon7 * c_memint_pancreas_igg * c_fcrnmemint_pancreas + koff7 * c_memintfr1_pancreas_igg + kdegab * c_memintfr1_pancreas_igg - kon7b * c_memintfr1_pancreas_igg * c_fcrnmemint_pancreas + koff7b * c_memintfr2_pancreas_igg + 2 * kdegab * c_memintfr2_pancreas_igg - kon7 * c_memint_pancreas * c_fcrnmemint_pancreas + koff7 * c_memintfr1_pancreas + kdegab * c_memintfr1_pancreas - kon7b * c_memintfr1_pancreas * c_fcrnmemint_pancreas + koff7b * c_memintfr2_pancreas + 2 * kdegab * c_memintfr2_pancreas)
    d/dt(fcrnmemint_thymus) <- clup_thymus * (1 - fr) * c_fcrnendorecyc_thymus * (1 - frecycle) - clup_thymus * (1 - fr) * c_fcrnmemint_thymus + vism_thymus * (-kon7 * c_memint_thymus_igg * c_fcrnmemint_thymus + koff7 * c_memintfr1_thymus_igg + kdegab * c_memintfr1_thymus_igg - kon7b * c_memintfr1_thymus_igg * c_fcrnmemint_thymus + koff7b * c_memintfr2_thymus_igg + 2 * kdegab * c_memintfr2_thymus_igg - kon7 * c_memint_thymus * c_fcrnmemint_thymus + koff7 * c_memintfr1_thymus + kdegab * c_memintfr1_thymus - kon7b * c_memintfr1_thymus * c_fcrnmemint_thymus + koff7b * c_memintfr2_thymus + 2 * kdegab * c_memintfr2_thymus)
    d/dt(fcrnmemint_spleen) <- clup_spleen * (1 - fr) * c_fcrnendorecyc_spleen * (1 - frecycle) - clup_spleen * (1 - fr) * c_fcrnmemint_spleen + vism_spleen * (-kon7 * c_memint_spleen_igg * c_fcrnmemint_spleen + koff7 * c_memintfr1_spleen_igg + kdegab * c_memintfr1_spleen_igg - kon7b * c_memintfr1_spleen_igg * c_fcrnmemint_spleen + koff7b * c_memintfr2_spleen_igg + 2 * kdegab * c_memintfr2_spleen_igg - kon7 * c_memint_spleen * c_fcrnmemint_spleen + koff7 * c_memintfr1_spleen + kdegab * c_memintfr1_spleen - kon7b * c_memintfr1_spleen * c_fcrnmemint_spleen + koff7b * c_memintfr2_spleen + 2 * kdegab * c_memintfr2_spleen)
    d/dt(fcrnmemint_other) <- clup_other * (1 - fr) * c_fcrnendorecyc_other * (1 - frecycle) - clup_other * (1 - fr) * c_fcrnmemint_other + vism_other * (-kon7 * c_memint_other_igg * c_fcrnmemint_other + koff7 * c_memintfr1_other_igg + kdegab * c_memintfr1_other_igg - kon7b * c_memintfr1_other_igg * c_fcrnmemint_other + koff7b * c_memintfr2_other_igg + 2 * kdegab * c_memintfr2_other_igg - kon7 * c_memint_other * c_fcrnmemint_other + koff7 * c_memintfr1_other + kdegab * c_memintfr1_other - kon7b * c_memintfr1_other * c_fcrnmemint_other + koff7b * c_memintfr2_other + 2 * kdegab * c_memintfr2_other)

    # ---- Central plasma and lymph node ----
    # Spleen, pancreas, small and large intestine drain through the liver, so
    # their flows leave the body pool carrying the LIVER vascular
    # concentration (Supplementary "Exogenous IgG in Central Plasma").
    d/dt(central) <- (plq_heart - lf_heart) * c_vp_heart + (plq_kidney - lf_kidney) * c_vp_kidney + (plq_muscle - lf_muscle) * c_vp_muscle + (plq_skin - lf_skin) * c_vp_skin + (plq_brain - lf_brain) * c_vp_brain + (plq_adipose - lf_adipose) * c_vp_adipose + (plq_thymus - lf_thymus) * c_vp_thymus + (plq_liver - lf_liver) * c_vp_liver + (plq_bone - lf_bone) * c_vp_bone + (plq_other - lf_other) * c_vp_other + (plq_spleen - lf_spleen) * c_vp_liver + (plq_pancreas - lf_pancreas) * c_vp_liver + (plq_small_intestine - lf_small_intestine) * c_vp_liver + (plq_large_intestine - lf_large_intestine) * c_vp_liver - (plq_lung + lf_lung) * c_central + llymph * c_lnode + ka * depot_sc
    d/dt(depot_sc) <- -ka * depot_sc
    # Subcutaneous bioavailability. The deposited listing carries F = 0.60 and
    # ka = 0.26/day as the assumed SC parameters used for mAb 11, the one
    # antibody with only subcutaneous clinical data (Methods).
    f(depot_sc) <- fsc
    d/dt(lnode) <- (1 - sigis) * lf_lung * c_is_lung + (1 - sigis) * lf_liver * c_is_liver + (1 - sigis) * lf_heart * c_is_heart + (1 - sigis) * lf_muscle * c_is_muscle + (1 - sigis) * lf_skin * c_is_skin + (1 - sigis) * lf_adipose * c_is_adipose + (1 - sigis) * lf_bone * c_is_bone + (1 - sigis) * lf_brain * c_is_brain + (1 - sigis) * lf_kidney * c_is_kidney + (1 - sigis) * lf_small_intestine * c_is_small_intestine + (1 - sigis) * lf_large_intestine * c_is_large_intestine + (1 - sigis) * lf_pancreas * c_is_pancreas + (1 - sigis) * lf_thymus * c_is_thymus + (1 - sigis) * lf_spleen * c_is_spleen + (1 - sigis) * lf_other * c_is_other - llymph * c_lnode
    d/dt(lnode_igg) <- (1 - sigis) * lf_lung * c_is_lung_igg + (1 - sigis) * lf_liver * c_is_liver_igg + (1 - sigis) * lf_heart * c_is_heart_igg + (1 - sigis) * lf_muscle * c_is_muscle_igg + (1 - sigis) * lf_skin * c_is_skin_igg + (1 - sigis) * lf_adipose * c_is_adipose_igg + (1 - sigis) * lf_bone * c_is_bone_igg + (1 - sigis) * lf_brain * c_is_brain_igg + (1 - sigis) * lf_kidney * c_is_kidney_igg + (1 - sigis) * lf_small_intestine * c_is_small_intestine_igg + (1 - sigis) * lf_large_intestine * c_is_large_intestine_igg + (1 - sigis) * lf_pancreas * c_is_pancreas_igg + (1 - sigis) * lf_thymus * c_is_thymus_igg + (1 - sigis) * lf_spleen * c_is_spleen_igg + (1 - sigis) * lf_other * c_is_other_igg - llymph * c_lnode_igg
    d/dt(auc_central) <- c_central

    # ---- Tissue biodistribution coefficients ----
    # Total antibody (free and bound, all subcompartments except the organ
    # vascular space) per unit organ volume, relative to plasma
    # (Madonna mAb_Biodistr_Coeff_in_tissue).
    kp_lung <- (vvm_lung * (c_memvas_lung + c_memvasfr1_lung + c_memvasfr2_lung + c_memvasns_lung) + ve7_lung * (c_endoearly_lung + c_endoearlyfr1_lung + c_endoearlyfr2_lung) + ve6a_lung * (c_endosort_lung + c_endosortfr1_lung + c_endosortfr2_lung) + ve7b_lung * (c_endorecyc_lung + c_endorecycfr1_lung + c_endorecycfr2_lung) + vism_lung * (c_memint_lung + c_memintfr1_lung + c_memintfr2_lung + c_memintns_lung) + vis_lung * c_is_lung) / vorg_lung / (1e-9 + c_central)
    kp_liver <- (vvm_liver * (c_memvas_liver + c_memvasfr1_liver + c_memvasfr2_liver + c_memvasns_liver) + ve7_liver * (c_endoearly_liver + c_endoearlyfr1_liver + c_endoearlyfr2_liver) + ve6a_liver * (c_endosort_liver + c_endosortfr1_liver + c_endosortfr2_liver) + ve7b_liver * (c_endorecyc_liver + c_endorecycfr1_liver + c_endorecycfr2_liver) + vism_liver * (c_memint_liver + c_memintfr1_liver + c_memintfr2_liver + c_memintns_liver) + vis_liver * c_is_liver) / vorg_liver / (1e-9 + c_central)
    kp_heart <- (vvm_heart * (c_memvas_heart + c_memvasfr1_heart + c_memvasfr2_heart + c_memvasns_heart) + ve7_heart * (c_endoearly_heart + c_endoearlyfr1_heart + c_endoearlyfr2_heart) + ve6a_heart * (c_endosort_heart + c_endosortfr1_heart + c_endosortfr2_heart) + ve7b_heart * (c_endorecyc_heart + c_endorecycfr1_heart + c_endorecycfr2_heart) + vism_heart * (c_memint_heart + c_memintfr1_heart + c_memintfr2_heart + c_memintns_heart) + vis_heart * c_is_heart) / vorg_heart / (1e-9 + c_central)
    kp_muscle <- (vvm_muscle * (c_memvas_muscle + c_memvasfr1_muscle + c_memvasfr2_muscle + c_memvasns_muscle) + ve7_muscle * (c_endoearly_muscle + c_endoearlyfr1_muscle + c_endoearlyfr2_muscle) + ve6a_muscle * (c_endosort_muscle + c_endosortfr1_muscle + c_endosortfr2_muscle) + ve7b_muscle * (c_endorecyc_muscle + c_endorecycfr1_muscle + c_endorecycfr2_muscle) + vism_muscle * (c_memint_muscle + c_memintfr1_muscle + c_memintfr2_muscle + c_memintns_muscle) + vis_muscle * c_is_muscle) / vorg_muscle / (1e-9 + c_central)
    kp_skin <- (vvm_skin * (c_memvas_skin + c_memvasfr1_skin + c_memvasfr2_skin + c_memvasns_skin) + ve7_skin * (c_endoearly_skin + c_endoearlyfr1_skin + c_endoearlyfr2_skin) + ve6a_skin * (c_endosort_skin + c_endosortfr1_skin + c_endosortfr2_skin) + ve7b_skin * (c_endorecyc_skin + c_endorecycfr1_skin + c_endorecycfr2_skin) + vism_skin * (c_memint_skin + c_memintfr1_skin + c_memintfr2_skin + c_memintns_skin) + vis_skin * c_is_skin) / vorg_skin / (1e-9 + c_central)
    kp_adipose <- (vvm_adipose * (c_memvas_adipose + c_memvasfr1_adipose + c_memvasfr2_adipose + c_memvasns_adipose) + ve7_adipose * (c_endoearly_adipose + c_endoearlyfr1_adipose + c_endoearlyfr2_adipose) + ve6a_adipose * (c_endosort_adipose + c_endosortfr1_adipose + c_endosortfr2_adipose) + ve7b_adipose * (c_endorecyc_adipose + c_endorecycfr1_adipose + c_endorecycfr2_adipose) + vism_adipose * (c_memint_adipose + c_memintfr1_adipose + c_memintfr2_adipose + c_memintns_adipose) + vis_adipose * c_is_adipose) / vorg_adipose / (1e-9 + c_central)
    kp_bone <- (vvm_bone * (c_memvas_bone + c_memvasfr1_bone + c_memvasfr2_bone + c_memvasns_bone) + ve7_bone * (c_endoearly_bone + c_endoearlyfr1_bone + c_endoearlyfr2_bone) + ve6a_bone * (c_endosort_bone + c_endosortfr1_bone + c_endosortfr2_bone) + ve7b_bone * (c_endorecyc_bone + c_endorecycfr1_bone + c_endorecycfr2_bone) + vism_bone * (c_memint_bone + c_memintfr1_bone + c_memintfr2_bone + c_memintns_bone) + vis_bone * c_is_bone) / vorg_bone / (1e-9 + c_central)
    kp_brain <- (vvm_brain * (c_memvas_brain + c_memvasfr1_brain + c_memvasfr2_brain + c_memvasns_brain) + ve7_brain * (c_endoearly_brain + c_endoearlyfr1_brain + c_endoearlyfr2_brain) + ve6a_brain * (c_endosort_brain + c_endosortfr1_brain + c_endosortfr2_brain) + ve7b_brain * (c_endorecyc_brain + c_endorecycfr1_brain + c_endorecycfr2_brain) + vism_brain * (c_memint_brain + c_memintfr1_brain + c_memintfr2_brain + c_memintns_brain) + vis_brain * c_is_brain) / vorg_brain / (1e-9 + c_central)
    kp_kidney <- (vvm_kidney * (c_memvas_kidney + c_memvasfr1_kidney + c_memvasfr2_kidney + c_memvasns_kidney) + ve7_kidney * (c_endoearly_kidney + c_endoearlyfr1_kidney + c_endoearlyfr2_kidney) + ve6a_kidney * (c_endosort_kidney + c_endosortfr1_kidney + c_endosortfr2_kidney) + ve7b_kidney * (c_endorecyc_kidney + c_endorecycfr1_kidney + c_endorecycfr2_kidney) + vism_kidney * (c_memint_kidney + c_memintfr1_kidney + c_memintfr2_kidney + c_memintns_kidney) + vis_kidney * c_is_kidney) / vorg_kidney / (1e-9 + c_central)
    kp_small_intestine <- (vvm_small_intestine * (c_memvas_small_intestine + c_memvasfr1_small_intestine + c_memvasfr2_small_intestine + c_memvasns_small_intestine) + ve7_small_intestine * (c_endoearly_small_intestine + c_endoearlyfr1_small_intestine + c_endoearlyfr2_small_intestine) + ve6a_small_intestine * (c_endosort_small_intestine + c_endosortfr1_small_intestine + c_endosortfr2_small_intestine) + ve7b_small_intestine * (c_endorecyc_small_intestine + c_endorecycfr1_small_intestine + c_endorecycfr2_small_intestine) + vism_small_intestine * (c_memint_small_intestine + c_memintfr1_small_intestine + c_memintfr2_small_intestine + c_memintns_small_intestine) + vis_small_intestine * c_is_small_intestine) / vorg_small_intestine / (1e-9 + c_central)
    kp_large_intestine <- (vvm_large_intestine * (c_memvas_large_intestine + c_memvasfr1_large_intestine + c_memvasfr2_large_intestine + c_memvasns_large_intestine) + ve7_large_intestine * (c_endoearly_large_intestine + c_endoearlyfr1_large_intestine + c_endoearlyfr2_large_intestine) + ve6a_large_intestine * (c_endosort_large_intestine + c_endosortfr1_large_intestine + c_endosortfr2_large_intestine) + ve7b_large_intestine * (c_endorecyc_large_intestine + c_endorecycfr1_large_intestine + c_endorecycfr2_large_intestine) + vism_large_intestine * (c_memint_large_intestine + c_memintfr1_large_intestine + c_memintfr2_large_intestine + c_memintns_large_intestine) + vis_large_intestine * c_is_large_intestine) / vorg_large_intestine / (1e-9 + c_central)
    kp_pancreas <- (vvm_pancreas * (c_memvas_pancreas + c_memvasfr1_pancreas + c_memvasfr2_pancreas + c_memvasns_pancreas) + ve7_pancreas * (c_endoearly_pancreas + c_endoearlyfr1_pancreas + c_endoearlyfr2_pancreas) + ve6a_pancreas * (c_endosort_pancreas + c_endosortfr1_pancreas + c_endosortfr2_pancreas) + ve7b_pancreas * (c_endorecyc_pancreas + c_endorecycfr1_pancreas + c_endorecycfr2_pancreas) + vism_pancreas * (c_memint_pancreas + c_memintfr1_pancreas + c_memintfr2_pancreas + c_memintns_pancreas) + vis_pancreas * c_is_pancreas) / vorg_pancreas / (1e-9 + c_central)
    kp_thymus <- (vvm_thymus * (c_memvas_thymus + c_memvasfr1_thymus + c_memvasfr2_thymus + c_memvasns_thymus) + ve7_thymus * (c_endoearly_thymus + c_endoearlyfr1_thymus + c_endoearlyfr2_thymus) + ve6a_thymus * (c_endosort_thymus + c_endosortfr1_thymus + c_endosortfr2_thymus) + ve7b_thymus * (c_endorecyc_thymus + c_endorecycfr1_thymus + c_endorecycfr2_thymus) + vism_thymus * (c_memint_thymus + c_memintfr1_thymus + c_memintfr2_thymus + c_memintns_thymus) + vis_thymus * c_is_thymus) / vorg_thymus / (1e-9 + c_central)
    kp_spleen <- (vvm_spleen * (c_memvas_spleen + c_memvasfr1_spleen + c_memvasfr2_spleen + c_memvasns_spleen) + ve7_spleen * (c_endoearly_spleen + c_endoearlyfr1_spleen + c_endoearlyfr2_spleen) + ve6a_spleen * (c_endosort_spleen + c_endosortfr1_spleen + c_endosortfr2_spleen) + ve7b_spleen * (c_endorecyc_spleen + c_endorecycfr1_spleen + c_endorecycfr2_spleen) + vism_spleen * (c_memint_spleen + c_memintfr1_spleen + c_memintfr2_spleen + c_memintns_spleen) + vis_spleen * c_is_spleen) / vorg_spleen / (1e-9 + c_central)
    kp_other <- (vvm_other * (c_memvas_other + c_memvasfr1_other + c_memvasfr2_other + c_memvasns_other) + ve7_other * (c_endoearly_other + c_endoearlyfr1_other + c_endoearlyfr2_other) + ve6a_other * (c_endosort_other + c_endosortfr1_other + c_endosortfr2_other) + ve7b_other * (c_endorecyc_other + c_endorecycfr1_other + c_endorecycfr2_other) + vism_other * (c_memint_other + c_memintfr1_other + c_memintfr2_other + c_memintns_other) + vis_other * c_is_other) / vorg_other / (1e-9 + c_central)

    # ---- Observation ----
    # Plasma concentration in ug/mL, the unit used in the published figures.
    Cc <- c_central * mwmab / 1000
    Cc ~ prop(propSd)
  })
}
