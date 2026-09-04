Yao_2023_sglt2_endpoints_mbma <- function() {
  description <- paste0(
    "MBMA. Class-level PK/PD/disease-endpoint model-based meta-analysis ",
    "linking sodium-glucose co-transporter-2 (SGLT2) inhibitor exposure to ",
    "fasting plasma glucose (FPG) and hemoglobin A1c (HbA1c) in type 2 ",
    "diabetes. Steady-state AUC(0-24 h) drives a shared-Emax, ",
    "drug-specific-EC50 model of the translational biomarker dUGEc (the ",
    "24 h change in urinary glucose excretion from baseline, normalised by ",
    "the baseline FPG); dUGEc lowers FPG linearly; FPG drives an ",
    "indirect-response turnover model of HbA1c. Both endpoints carry a ",
    "treatment-history-specific placebo response and a linear ",
    "disease-progression term. Fitted to 27 dUGEc, 848 FPG and 1219 HbA1c ",
    "summary-level data points digitised from 80 published dapagliflozin, ",
    "canagliflozin and empagliflozin trials, and externally validated ",
    "against ertugliflozin. Variability is BETWEEN-STUDY (inter-study), ",
    "encoded as study-level etas, so the model simulates study-arm mean ",
    "endpoint trajectories and is NOT suitable for individual-subject ",
    "simulation. Exposure enters through the AUC_DAPA / AUC_CANA / AUC_EMPA ",
    "covariates, which the companion PK models ",
    "modellib('Yao_2023_dapagliflozin_mbma'), ",
    "modellib('Yao_2023_canagliflozin_mbma') and ",
    "modellib('Yao_2023_empagliflozin_mbma') supply."
  )

  reference <- paste(
    "Yao X, Zhou J, Song L, Ren Y, Hu P, Liu D.",
    "A model-based meta analysis study of sodium glucose co-transporter-2",
    "inhibitors. CPT Pharmacometrics Syst Pharmacol. 2023;12(4):487-499.",
    "doi:10.1002/psp4.12934.",
    sep = " "
  )
  vignette <- "Yao_2023_sglt2_inhibitors"
  units <- list(
    time          = "week",
    dosing        = "n/a (exposure enters through the AUC_<drug> covariates)",
    concentration = paste0(
      "FPG in mg/dL; HbA1c in % (NGSP); dUGEc in g/(mg/dL). The AUC_<drug> ",
      "covariates are in ng*h/mL. Output Cc is unused."
    )
  )

  paper_specific_compartments <- c("hba1c_drug")

  covariateData <- list(
    AUC_DAPA = list(
      description        = paste0(
        "Dapagliflozin area under the plasma concentration-time curve over ",
        "the 24 h dosing interval at steady state. Study-arm-level, not ",
        "individual-level. Set to 0 for placebo arms and for arms treated ",
        "with a different SGLT2 inhibitor."
      ),
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Yao 2023 PK/PD model section: 'Drug exposure in patients with T2DM ",
        "expressed as area under the concentration-time curve (AUC) was ",
        "simulated using the established population PK model for each ",
        "agent. AUC at steady status was selected to drive PD changes.' ",
        "Compute it from modellib('Yao_2023_dapagliflozin_mbma') as total ",
        "daily dose (mg) divided by CL/F (L/h), times 1000 to convert ",
        "mg*h/L to ng*h/mL -- exact for a linear model at steady state. ",
        "At the 10 mg label dose this gives 1000 * 10 / 19.5 = 513 ng*h/mL."
      ),
      source_name        = "AUC0-24h"
    ),
    AUC_CANA = list(
      description        = paste0(
        "Canagliflozin area under the plasma concentration-time curve over ",
        "the 24 h dosing interval at steady state. Study-arm-level, not ",
        "individual-level. Set to 0 for placebo arms and for arms treated ",
        "with a different SGLT2 inhibitor."
      ),
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Compute it from modellib('Yao_2023_canagliflozin_mbma') as total ",
        "daily dose (mg) divided by CL/F (L/h), times 1000. At the 300 mg ",
        "label dose this gives 1000 * 300 / 12.0 = 25000 ng*h/mL."
      ),
      source_name        = "AUC0-24h"
    ),
    AUC_EMPA = list(
      description        = paste0(
        "Empagliflozin area under the plasma concentration-time curve over ",
        "the 24 h dosing interval at steady state. Study-arm-level, not ",
        "individual-level. Set to 0 for placebo arms and for arms treated ",
        "with a different SGLT2 inhibitor."
      ),
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Compute it from modellib('Yao_2023_empagliflozin_mbma') as total ",
        "daily dose (mg) divided by CL/F (L/h), times 1000. At the 25 mg ",
        "label dose this gives 1000 * 25 / 4.25 = 5882 ng*h/mL. Distinct ",
        "from DOSE_EMPA_MGD, which is the daily dose itself and is used by ",
        "the Baron 2016 and Riggs 2014 empagliflozin exposure-response ",
        "models."
      ),
      source_name        = "AUC0-24h"
    ),
    TRT_T2DM_NAIVE = list(
      description        = paste0(
        "1 = the study arm enrolled patients naive to oral hypoglycemic ",
        "agents; 0 = otherwise."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other treatment-history stratum)",
      notes              = paste0(
        "One of four mutually exclusive study-arm treatment-history ",
        "indicators; exactly one of TRT_T2DM_NAIVE, TRT_T2DM_NONNAIVE, ",
        "TRT_T2DM_ADDON and TRT_T2DM_MIXED is 1 on any record. Yao 2023 ",
        "Model development: 'According to treatment history and therapy ",
        "regimens, patients were divided into different treatment types: ",
        "naive therapy (patients naive to oral hypoglycemic agents), ",
        "non-naive therapy (patients had taken hypoglycemic agents but had ",
        "undergone a more than 2-week washout period), add-on therapy ",
        "(patients were in combination treatment in the study), and mixed ",
        "(contained both naive and non-naive therapy and unknown treatment ",
        "type).' The stratum selects the placebo-response magnitude on both ",
        "endpoints (Pfmax1-4 and Phmax1-4); every other parameter is shared ",
        "('Treatment type is not significant on other parameters here')."
      ),
      source_name        = "treatment type = naive"
    ),
    TRT_T2DM_NONNAIVE = list(
      description        = paste0(
        "1 = the study arm enrolled patients who had taken hypoglycemic ",
        "agents but had undergone a washout period of more than 2 weeks; ",
        "0 = otherwise."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other treatment-history stratum)",
      notes              = "See TRT_T2DM_NAIVE. Yao 2023 Model development.",
      source_name        = "treatment type = non-naive"
    ),
    TRT_T2DM_ADDON = list(
      description        = paste0(
        "1 = the study arm enrolled patients receiving the study drug as ",
        "add-on to ongoing antihyperglycemic combination treatment; ",
        "0 = otherwise."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other treatment-history stratum)",
      notes              = paste0(
        "See TRT_T2DM_NAIVE. This is the only stratum in which the placebo ",
        "FPG response is a decrease rather than an increase (Yao 2023 ",
        "Results: 'patients got their FPG under control with the help of ",
        "other agents in combination')."
      ),
      source_name        = "treatment type = add-on"
    ),
    TRT_T2DM_MIXED = list(
      description        = paste0(
        "1 = the study arm pooled naive and non-naive patients, or the ",
        "treatment history was not reported; 0 = otherwise."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other treatment-history stratum)",
      notes              = paste0(
        "See TRT_T2DM_NAIVE. This stratum is a meta-analysis artefact: it ",
        "exists because some published arms do not resolve treatment ",
        "history, so it cannot be reconstructed from prior-therapy and ",
        "concomitant-therapy indicators."
      ),
      source_name        = "treatment type = mixed"
    )
  )

  compartmentData <- list(
    hba1c_drug = list(
      analyte  = "Glycated hemoglobin (drug-attributable deviation from the placebo trajectory)",
      units    = "% (NGSP)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24817L,
    n_studies      = 80L,
    age_range      = paste0(
      "endpoint cohorts: mean 55.8-57.6 years across the placebo and three ",
      "drug groups (Yao 2023 Table S2)"
    ),
    weight_range   = "endpoint cohorts: mean 79.0-84.9 kg (Yao 2023 Table S2)",
    sex_female_pct = 43.3,
    disease_state  = paste0(
      "Type 2 diabetes mellitus with normal or mildly impaired renal ",
      "function (glomerular filtration rate above 60 mL/min/1.73 m2). ",
      "Studies in patients with moderate or severe kidney impairment or ",
      "hepatic insufficiency were excluded, as were trials using insulin. ",
      "Baseline FPG about 160 mg/dL and baseline HbA1c about 7.9%."
    ),
    dose_range     = paste0(
      "Placebo plus dapagliflozin 1-50 mg, canagliflozin 50-300 mg and ",
      "empagliflozin 1-100 mg once daily; treatment durations up to 104 ",
      "weeks (Yao 2023 Table S1)"
    ),
    regions        = "International (published clinical trials indexed on PubMed to July 2016)",
    notes          = paste0(
      "Model-based meta-analysis: the unit of observation is a published ",
      "study-arm mean, not an individual measurement. 27 summary-level ",
      "dUGEc, 848 FPG (195 placebo, 653 drug) and 1219 HbA1c (290 placebo, ",
      "929 drug) data points were extracted. n_subjects is the sum of the ",
      "FPG/HbA1c endpoint cohorts in Yao 2023 Table S1 (dapagliflozin ",
      "8324, canagliflozin 7004, empagliflozin 9489); patients contributing ",
      "PK and PD data are counted in the companion PK models instead. A ",
      "further 26 HbA1c data points from ertugliflozin were held out for ",
      "external validation and are not part of the fit."
    )
  )

  # ==========================================================================
  # ini(): Yao 2023 Table 2 (PK/PD and PD/end point model parameter estimates).
  #
  # Sign conventions follow the table, which reports SIGNED placebo effects
  # (negative = a fall from baseline) and a signed drug slope. See the
  # model() block and the vignette Assumptions and deviations section for the
  # two places where the printed equations had to be reconciled against the
  # paper's own Discussion text and Figure 3 visual predictive checks.
  #
  # Between-study variability (Yao 2023 PD/end point model section):
  #   - ADDITIVE for Pfmax and Phmax ("Additive inter-study variability was
  #     proposed for Pfmax and Phmax"), so the table's IIV(%) is read as a
  #     relative SD against the typical value and omega^2 = (IIV/100 *
  #     |value|)^2. This is the same convention used by the sibling T2DM MBMA
  #     modellib('Li_2015_taspoglutide_mbma').
  #   - EXPONENTIAL for every other parameter ("inter-study variability of
  #     other parameters in both of FPG and HbA1c models were expressed as
  #     exponential form as assumed to be log-normally distribution"), so
  #     omega^2 = log((IIV/100)^2 + 1).
  # Etas are named eta_study_* so they cannot be mistaken for the popPK
  # between-subject convention.
  # ==========================================================================
  ini({
    # ----- PK/PD layer: exposure to dUGEc (Yao 2023 Equation 2, Table 2) ----
    lemax <- log(0.606)
    label("Maximal drug effect on dUGEc, shared across the three drugs (g/(mg/dL))")  # Yao 2023 Table 2 (RSE 4.40%)

    lec50_dapagliflozin <- log(56.6)
    label("Steady-state AUC(0-24 h) giving half-maximal dUGEc for dapagliflozin (ng*h/mL)")  # Yao 2023 Table 2 (RSE 27.1%)

    lec50_canagliflozin <- log(2310)
    label("Steady-state AUC(0-24 h) giving half-maximal dUGEc for canagliflozin (ng*h/mL)")  # Yao 2023 Table 2 (RSE 23.3%)

    lec50_empagliflozin <- log(841)
    label("Steady-state AUC(0-24 h) giving half-maximal dUGEc for empagliflozin (ng*h/mL)")  # Yao 2023 Table 2 (RSE 30.4%)

    # ----- FPG endpoint (Yao 2023 Equations 3-4, Table 2) -------------------
    lrbase_fpg <- log(160)
    label("Population baseline fasting plasma glucose (mg/dL)")  # Yao 2023 Table 2 (RSE 1.00%)

    pfmax_naive <- 1.45
    label("Maximal placebo effect on FPG, naive stratum (mg/dL; signed, positive = rise)")  # Yao 2023 Table 2 Pfmax1 (RSE 26.0%)

    pfmax_nonnaive <- 1.90
    label("Maximal placebo effect on FPG, non-naive stratum (mg/dL; signed, positive = rise)")  # Yao 2023 Table 2 Pfmax2 (RSE 111%)

    pfmax_addon <- -1.37
    label("Maximal placebo effect on FPG, add-on stratum (mg/dL; signed, negative = fall)")  # Yao 2023 Table 2 Pfmax3 (RSE 74.0%)

    pfmax_mixed <- 4.30
    label("Maximal placebo effect on FPG, mixed stratum (mg/dL; signed, positive = rise)")  # Yao 2023 Table 2 Pfmax4 (RSE 55.0%)

    lkfp <- log(0.340)
    label("First-order onset rate constant of the FPG response Kfp (1/week)")  # Yao 2023 Table 2 (RSE 55.0%; half-life 2 weeks)

    disfp <- 3.13
    label("Disease-progression slope on FPG (mg/dL per 100 weeks)")  # Yao 2023 Table 2 DISfp (RSE 41.0%)

    slopefd <- -43.3
    label("Linear drug effect of dUGEc on FPG (mg/dL per g/(mg/dL); signed, negative = fall)")  # Yao 2023 Table 2 SLOPEfd (RSE 4.10%)

    # ----- HbA1c endpoint (Yao 2023 Equations 5-8, Table 2) -----------------
    lrbase_hba1c <- log(7.92)
    label("Population baseline HbA1c (%)")  # Yao 2023 Table 2 (RSE 1.00%)

    phmax_naive <- -0.200
    label("Maximal placebo effect on HbA1c, naive stratum (%; signed, negative = fall)")  # Yao 2023 Table 2 Phmax1 (RSE 17.0%)

    phmax_nonnaive <- 0.0510
    label("Maximal placebo effect on HbA1c, non-naive stratum (%; signed, positive = rise)")  # Yao 2023 Table 2 Phmax2 (RSE 109%)

    phmax_addon <- -0.230
    label("Maximal placebo effect on HbA1c, add-on stratum (%; signed, negative = fall)")  # Yao 2023 Table 2 Phmax3 (RSE 16.0%)

    phmax_mixed <- -0.06
    label("Maximal placebo effect on HbA1c, mixed stratum (%; signed, negative = fall)")  # Yao 2023 Table 2 Phmax4 (RSE 81.0%)

    lkhp <- log(0.240)
    label("First-order onset rate constant of the HbA1c placebo response Khp (1/week)")  # Yao 2023 Table 2 (RSE 21.0%)

    dishp <- 0.310
    label("Disease-progression slope on HbA1c (% per 100 weeks)")  # Yao 2023 Table 2 DIShp (RSE 33.0%)

    lkout <- log(0.200)
    label("First-order elimination rate constant of HbA1c Kout (1/week)")  # Yao 2023 Table 2 (RSE 3.00%; half-life 3.5 weeks)

    lkin2 <- log(0.500)
    label("FPG-independent zero-order production rate of HbA1c Kin2 (%/week)")  # Yao 2023 Table 2 (RSE 5.00%)

    # ----- Between-study variability; diagonal ------------------------------
    # Additive on the placebo maxima
    eta_study_pfmax_naive  ~ 0.00692723   # Yao 2023 Table 2; additive, (0.0574 * 1.45)^2
    eta_study_pfmax_addon  ~ 0.00875553   # Yao 2023 Table 2; additive, (0.0683 * 1.37)^2
    eta_study_pfmax_mixed  ~ 0.00139951   # Yao 2023 Table 2; additive, (0.0087 * 4.30)^2
    eta_study_phmax_naive  ~ 2.56e-08     # Yao 2023 Table 2; additive, (0.0008 * 0.200)^2
    eta_study_phmax_addon  ~ 8.9401e-08   # Yao 2023 Table 2; additive, (0.0013 * 0.230)^2
    eta_study_phmax_mixed  ~ 1.44e-08     # Yao 2023 Table 2; additive, (0.0020 * 0.06)^2
    # No inter-study variability is reported for Pfmax2 or Phmax2 (the
    # non-naive stratum; Table 2 IIV column "-").

    # Exponential on the remaining parameters
    eta_study_lrbase_fpg   ~ 0.00291176   # Yao 2023 Table 2; log(0.0540^2 + 1)
    eta_study_lkfp         ~ 0.167486     # Yao 2023 Table 2; log(0.427^2 + 1)
    eta_study_slopefd      ~ 0.0818233    # Yao 2023 Table 2; log(0.292^2 + 1)
    eta_study_lrbase_hba1c ~ 0.00151984   # Yao 2023 Table 2; log(0.0390^2 + 1)
    eta_study_lkhp         ~ 0.0829026    # Yao 2023 Table 2; log(0.294^2 + 1)
    eta_study_dishp        ~ 1.05386      # Yao 2023 Table 2; log(1.367^2 + 1)
    eta_study_lkout        ~ 0.0265407    # Yao 2023 Table 2; log(0.164^2 + 1)
    # No inter-study variability is reported for DISfp or Kin2 (Table 2 IIV
    # column "-").

    # ----- Residual error ---------------------------------------------------
    # All three residuals are reported as variances weighted by the square
    # root of the study sample size (Yao 2023 PD/end point model section;
    # Appendix S1 Equation S14). The values below are back-transformed to SD
    # at unit study weight; see the vignette Assumptions and deviations
    # section for the weighting.
    propSd_UGEc <- 0.471169
    label("Proportional residual SD on dUGEc at unit study weight (fraction)")  # Yao 2023 Table 2 PK/PD sigma^2_pro; sqrt(0.222)

    addSd_UGEc <- 0.254165
    label("Additive residual SD on dUGEc at unit study weight (g/(mg/dL))")  # Yao 2023 Table 2 PK/PD sigma^2_add; sqrt(0.0646)

    addSd_FPG <- 0.181659
    label("Additive residual SD on FPG at unit study weight (mg/dL)")  # Yao 2023 Table 2 sigma^2_add,FPG; sqrt(0.0330)

    addSd_HbA1c <- 0.0774597
    label("Additive residual SD on HbA1c at unit study weight (%)")  # Yao 2023 Table 2 sigma^2_add,HbA1c; sqrt(0.006)
  })

  model({
    # ======================================================================
    # 1. Study-level (inter-study-perturbed) parameters
    # ======================================================================
    emax      <- exp(lemax)
    ec50_dapa <- exp(lec50_dapagliflozin)
    ec50_cana <- exp(lec50_canagliflozin)
    ec50_empa <- exp(lec50_empagliflozin)

    fpgbase   <- exp(lrbase_fpg   + eta_study_lrbase_fpg)
    hba1cbase <- exp(lrbase_hba1c + eta_study_lrbase_hba1c)
    kfp       <- exp(lkfp  + eta_study_lkfp)
    khp       <- exp(lkhp  + eta_study_lkhp)
    kout      <- exp(lkout + eta_study_lkout)
    kin2      <- exp(lkin2)

    # Signed parameters carry their exponential inter-study variability as a
    # multiplicative factor so the sign is preserved.
    slopefd_i <- slopefd * exp(eta_study_slopefd)
    dishp_i   <- dishp   * exp(eta_study_dishp)

    # ======================================================================
    # 2. Treatment-history stratum selection
    #
    # Exactly one of the four indicators is 1 on any record, so each of these
    # sums picks out that stratum's placebo maximum. Yao 2023 estimated the
    # placebo model separately in each of the four groups and found
    # "Treatment type is not significant on other parameters here", so no
    # other parameter is stratified.
    # ======================================================================
    pfmax_naive_i <- pfmax_naive + eta_study_pfmax_naive
    pfmax_addon_i <- pfmax_addon + eta_study_pfmax_addon
    pfmax_mixed_i <- pfmax_mixed + eta_study_pfmax_mixed
    phmax_naive_i <- phmax_naive + eta_study_phmax_naive
    phmax_addon_i <- phmax_addon + eta_study_phmax_addon
    phmax_mixed_i <- phmax_mixed + eta_study_phmax_mixed

    pfmax <-
      pfmax_naive_i  * TRT_T2DM_NAIVE +
      pfmax_nonnaive * TRT_T2DM_NONNAIVE +
      pfmax_addon_i  * TRT_T2DM_ADDON +
      pfmax_mixed_i  * TRT_T2DM_MIXED

    phmax <-
      phmax_naive_i  * TRT_T2DM_NAIVE +
      phmax_nonnaive * TRT_T2DM_NONNAIVE +
      phmax_addon_i  * TRT_T2DM_ADDON +
      phmax_mixed_i  * TRT_T2DM_MIXED

    # ======================================================================
    # 3. Translational biomarker dUGEc (Yao 2023 Equations 1-2)
    #
    #   dUGEc = (UGE - UGE_baseline) / FPG_baseline          (Equation 1)
    #   dUGEc = Emax * AUC(0-24 h) / (EC50 + AUC(0-24 h))    (Equation 2)
    #
    # The three drugs share one Emax and differ only in EC50, so the sum
    # below is the same shape as the multi-drug MBMA precedent
    # modellib('Dodds_2013_psoriasis_biologics_mbma'): every study arm
    # received at most ONE SGLT2 inhibitor, so at most one term is non-zero.
    # A placebo arm sets all three AUC covariates to 0, every term collapses
    # to exactly 0, and the endpoints reduce to their placebo trajectories.
    # Setting two AUC columns non-zero simultaneously is outside the source's
    # calibration and would make the model ADD the two drug effects.
    # ======================================================================
    UGEc <- emax * (
      AUC_DAPA / (ec50_dapa + AUC_DAPA) +
      AUC_CANA / (ec50_cana + AUC_CANA) +
      AUC_EMPA / (ec50_empa + AUC_EMPA)
    )

    # ======================================================================
    # 4. FPG (Yao 2023 Equations 3-4)
    #
    #   FPG_placebo = FPG_baseline + Pfmax * (1 - exp(-Kfp * t)) + DISfp * t
    #   FPG         = FPG_placebo + SLOPEfd * dUGEc                (as printed)
    #
    # RECTIFICATION: the drug term carries the SAME (1 - exp(-Kfp * t)) onset
    # factor as the placebo term. Equation 4 as printed makes the drug effect
    # an instantaneous step at t = 0, which is contradicted by the paper's own
    # outputs in two independent places:
    #   (a) Discussion: "The FPG responses in drug effects with a Kfp of 0.34
    #       weeks-1 or a half-life of 2 weeks indicated a 2-week continuous
    #       treatment could show a significant decrease in FPG" -- Kfp is
    #       attributed to the DRUG effect, but Kfp appears only in Equation 3.
    #   (b) Figure 3 panels (c) and (d), the active-arm FPG visual predictive
    #       checks: the simulated median starts at the 160 mg/dL baseline at
    #       t = 0 and reaches its ~137 mg/dL plateau at about 9 weeks, which
    #       is 3 / Kfp = 3 / 0.340 = 8.8 weeks. An instantaneous step would
    #       put the simulated median at the plateau already at t = 0.
    # The rectification is the minimal one: it multiplies the drug term by the
    # onset factor and changes nothing else. See the vignette Assumptions and
    # deviations section.
    # ======================================================================
    onsetFpg   <- 1 - exp(-kfp * t)
    fpgPlacebo <- fpgbase + pfmax * onsetFpg + disfp / 100 * t
    fpg        <- fpgPlacebo + slopefd_i * UGEc * onsetFpg

    # ======================================================================
    # 5. HbA1c (Yao 2023 Equations 5-8)
    #
    #   HbA1c_placebo = HbA1c_baseline + Phmax * (1 - exp(-Khp * t))
    #                   + DIShp * t                              (Equation 5)
    #   Kin           = Kout * HbA1c_baseline - Kin2             (Equation 6)
    #   d(HbA1c_drug)/dt = FPG / FPG_baseline * Kin + Kin2
    #                      - Kout * HbA1c,  HbA1c_drug(0) = 0    (Equation 7)
    #   HbA1c         = HbA1c_placebo + HbA1c_drug               (Equation 8)
    #
    # RECTIFICATION: Equation 5 is printed with a MINUS sign on Phmax
    # ("HbA1c_baseline - Phmax * (1 - exp(-Khp * t))") while Table 2 reports
    # SIGNED Phmax values whose sign already encodes the direction. Taking
    # both literally double-negates: Phmax1 = -0.200 would make the naive
    # placebo arm RISE by 0.2%. The paper's own text says the opposite
    # ("the estimated Phmax value in HbA1c for naive groups (Phmax1), add-on
    # groups (Phmax3), and mixed groups (Phmax4) were -0.20%, -0.23%, and
    # -0.06% decrease from baseline while that of non-naive groups (Phmax2)
    # was 0.051% increase from baseline"), and so does Figure 3 panel (e):
    # the simulated placebo median dips about 0.15% below baseline by week 12
    # and only then climbs back above it under disease progression, which is
    # exactly what "+ Phmax" with Phmax1 = -0.200, Khp = 0.240 and
    # DIShp = 0.310%/100 weeks predicts (7.92 -> 7.77 at week 12 -> 8.03 at
    # week 100). The companion FPG Equation 3 uses "+ Pfmax" with the same
    # signed-value convention, so the plus sign is also the internally
    # consistent one. This model therefore uses "+ Phmax".
    #
    # READING OF Equation 7: the unsubscripted "HbA1c" in the elimination
    # term is the TOTAL of Equation 8 (placebo + drug), which is the literal
    # reading of the symbol and the only one under which the stated initial
    # condition HbA1c_drug(0) = 0 is stationary in the absence of drug and
    # placebo effects. Reading it as HbA1c_drug alone would drive that state
    # to HbA1c_baseline and double the predicted HbA1c.
    # ======================================================================
    hba1cPlacebo <- hba1cbase + phmax * (1 - exp(-khp * t)) + dishp_i / 100 * t
    kin          <- kout * hba1cbase - kin2
    hba1c        <- hba1cPlacebo + hba1c_drug

    hba1c_drug(0)    <- 0
    d/dt(hba1c_drug) <- fpg / fpgbase * kin + kin2 - kout * hba1c

    # ======================================================================
    # 6. Outputs and residual error
    # ======================================================================
    FPG   <- fpg
    HbA1c <- hba1c

    UGEc  ~ add(addSd_UGEc) + prop(propSd_UGEc)
    FPG   ~ add(addSd_FPG)
    HbA1c ~ add(addSd_HbA1c)
  })
}
