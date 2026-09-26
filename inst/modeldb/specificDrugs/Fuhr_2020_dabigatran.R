Fuhr_2020_dabigatran <- function() {
  description <- paste(
    "Direct concentration-effect pharmacodynamic model relating the unbound",
    "sum dabigatran plasma concentration (unbound dabigatran plus unbound",
    "dabigatran acyl-glucuronide, i.e. drug bound neither to plasma proteins",
    "nor to idarucizumab) to four coagulation times in adults: activated",
    "partial thromboplastin time (aPTT) and thrombin time (TT) as a combined",
    "linear plus Emax function, and diluted thrombin time (dTT) and ecarin",
    "clotting time (ECT) as linear functions (Fuhr 2020 ESM Equations 16-19).",
    "Only this PD layer is reproduced here. The idarucizumab and dabigatran",
    "whole-body PBPK layers, the idarucizumab-dabigatran binding model and the",
    "hemodialysis extension were built in the PK-Sim / MoBi platform and depend",
    "on platform-database organ volumes, blood flows and partition",
    "coefficients that the paper does not print, so they are not reproduced",
    "(see the vignette 'Assumptions and deviations'). The unbound sum",
    "dabigatran concentration is therefore supplied by the user as the",
    "time-varying covariate CU_DABIGATRAN_UM. No inter-individual variability",
    "and no residual-error model are reported (the PD data were digitised",
    "study-arm means), so both are encoded as zero.",
    sep = " "
  )
  reference <- paste(
    "Fuhr LM, Hanke N, Meibohm B, Lehr T. (2020). Effective Removal of",
    "Dabigatran by Idarucizumab or Hemodialysis: A Physiologically Based",
    "Pharmacokinetic Modeling Analysis. Clinical Pharmacokinetics",
    "59:809-825. doi:10.1007/s40262-019-00857-y.",
    "PD equations from Electronic Supplementary Material 1, section",
    "'Pharmacodynamic modeling', Equations 16-19.",
    sep = " "
  )
  vignette <- "Fuhr_2020_dabigatran"
  units <- list(
    time = "h",
    dosing = "(not applicable; the unbound sum dabigatran concentration is supplied as the CU_DABIGATRAN_UM covariate, not as a dose record)",
    concentration = "umol/L (unbound sum dabigatran via CU_DABIGATRAN_UM); the PD outputs aPTT, dTT, ECT and TT are in seconds"
  )

  covariateData <- list(
    CU_DABIGATRAN_UM = list(
      description = "Unbound sum dabigatran plasma concentration (unbound dabigatran plus unbound dabigatran acyl-glucuronide), umol/L, supplied as the time-varying driver of the four coagulation-time equations",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fuhr 2020 ESM Equation 3 defines the driver as",
        "'DABunb sum [umol/L] = [DAB] * fu DAB + [DABG] * fu DABG', the sum of",
        "the dabigatran and dabigatran-glucuronide plasma concentrations not",
        "bound to plasma proteins, with fu assumed equal for the two moieties;",
        "drug captured by idarucizumab is also excluded (main text section 2.3:",
        "'neither bound to plasma proteins nor captured by idarucizumab'). The",
        "ESM Equations 16-19 state '[DAB] = unbound sum dabigatran plasma",
        "concentrations'; the umol/L unit is confirmed by Figure 6 (left",
        "column, x-axis 'Unb sum DAB plasma [ng/mL]'): ECT = 209.70 * C + 34.85",
        "reaches the plotted 130 s at about 0.45 umol/L = 214 ng/mL, with the",
        "dabigatran molecular weight 471.5 g/mol (ESM: 'dabigatran (MW = 472",
        "g/mol)'). Divide an unbound concentration in ng/mL by 471.5 before",
        "supplying it here. The fitted data span 0 to about 220 ng/mL",
        "(0 to about 0.47 umol/L; Figure 6 left column). In the source the",
        "driver is the prediction of the authors' PK-Sim / MoBi PBPK model,",
        "which is not reproduced in nlmixr2lib; users supply it from their",
        "own dabigatran PK model (converted to unbound with plasma protein",
        "binding 35%, Fuhr 2020 Introduction, i.e. fu = 0.65) or from measured",
        "unbound concentrations. Set to 0 when no dabigatran is present or when",
        "all circulating dabigatran is captured by idarucizumab."
      ),
      source_name = "[DAB] / DABunb sum (Fuhr 2020 ESM Equations 3 and 16-19)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 3L,
    age_range = "study-arm means 24-72 years (Fuhr 2020 Table 1)",
    weight_range = "study-arm means 57-90 kg (Fuhr 2020 Table 1)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Caucasian and Japanese (Fuhr 2020 Table 1)",
    disease_state = paste(
      "Healthy middle-aged and elderly Caucasian adults, Caucasians with mild",
      "or moderate renal impairment, and healthy Japanese men, all on",
      "steady-state oral dabigatran etexilate, with or without subsequent",
      "intravenous idarucizumab"
    ),
    dose_range = paste(
      "Dabigatran etexilate 220 mg twice daily (healthy and elderly) or 150 mg",
      "twice daily (renal impairment) for 3.5 days, followed in most arms by",
      "idarucizumab 1000-7500 mg IV (bolus or 60-min infusion) (Fuhr 2020",
      "Table 1 and ESM Table S8)"
    ),
    regions = "Germany / Europe and Japan",
    renal_function = "creatinine clearance study-arm means 59-139 mL/min (Fuhr 2020 Table 1)",
    notes = paste(
      "The PD equations were fitted to digitised study-arm mean",
      "concentration-effect data pooled over 26 dabigatran-treated arms from",
      "Glund 2015, Glund 2016 / 2017 and Yasaka 2017 (Fuhr 2020 ESM Table S8),",
      "including coagulation measurements before and after idarucizumab.",
      "The number of individual subjects contributing to the PD fit is not",
      "reported (n_subjects = NA). Table 1 reports 0% female in the healthy",
      "Caucasian and Japanese studies and 25-50% female in the Glund 2017",
      "age / renal-impairment study; no pooled percentage is given. Mean",
      "relative deviations of the predicted coagulation times were 1.16",
      "(aPTT), 1.11 (dTT), 1.16 (ECT) and 1.27 (TT) (Fuhr 2020 section 3.4)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Fuhr 2020 ESM section 'Pharmacodynamic modeling', Equations 16-19,
    # with [DAB] = unbound sum dabigatran plasma concentration (umol/L):
    #
    #   aPTT = 17.94 * [DAB] / (0.04 + [DAB]) + 61.67 * [DAB] + 30.21  (16)
    #   dTT  = 88.30 * [DAB] + 31.59                                    (17)
    #   ECT  = 209.70 * [DAB] + 34.85                                   (18)
    #   TT   = 106.28 * [DAB] / (0.19 + [DAB]) + 144.02 * [DAB] + 12.93 (19)
    #
    # Main text section 3.4: 'a combined linear and Emax model was
    # selected ... to aPTT and TT coagulation times, while a linear model
    # was sufficient to describe the effect of unbound sum dabigatran on
    # dTT and ECT'. All coagulation times are in seconds (Figure 6 axes).
    # The paper reports no standard errors, IIV or residual-error model
    # for these fits.
    # ------------------------------------------------------------------

    # aPTT (ESM Equation 16)
    lrbase_aPTT <- log(30.21)
    label("Log of baseline aPTT at zero unbound sum dabigatran (s)") # ESM Eq. 16 intercept 30.21
    lemax_aPTT <- log(17.94)
    label("Log of Emax of the saturable aPTT component (s)") # ESM Eq. 16 Emax 17.94
    lec50_aPTT <- log(0.04)
    label("Log of EC50 of the saturable aPTT component (umol/L)") # ESM Eq. 16 EC50 0.04
    lslope_aPTT <- log(61.67)
    label("Log of the linear aPTT slope (s per umol/L)") # ESM Eq. 16 slope 61.67

    # dTT (ESM Equation 17)
    lrbase_dTT <- log(31.59)
    label("Log of baseline dTT at zero unbound sum dabigatran (s)") # ESM Eq. 17 intercept 31.59
    lslope_dTT <- log(88.30)
    label("Log of the linear dTT slope (s per umol/L)") # ESM Eq. 17 slope 88.30

    # ECT (ESM Equation 18)
    lrbase_ECT <- log(34.85)
    label("Log of baseline ECT at zero unbound sum dabigatran (s)") # ESM Eq. 18 intercept 34.85
    lslope_ECT <- log(209.70)
    label("Log of the linear ECT slope (s per umol/L)") # ESM Eq. 18 slope 209.70

    # TT (ESM Equation 19)
    lrbase_TT <- log(12.93)
    label("Log of baseline TT at zero unbound sum dabigatran (s)") # ESM Eq. 19 intercept 12.93
    lemax_TT <- log(106.28)
    label("Log of Emax of the saturable TT component (s)") # ESM Eq. 19 Emax 106.28
    lec50_TT <- log(0.19)
    label("Log of EC50 of the saturable TT component (umol/L)") # ESM Eq. 19 EC50 0.19
    lslope_TT <- log(144.02)
    label("Log of the linear TT slope (s per umol/L)") # ESM Eq. 19 slope 144.02

    # No residual-error model is reported (the fits are to digitised
    # study-arm means; model performance is summarised only as mean
    # relative deviation), so each additive SD is fixed at zero.
    addSd_aPTT <- fixed(0)
    label("Additive residual SD on aPTT (s; zero, not reported by the source)")
    addSd_dTT <- fixed(0)
    label("Additive residual SD on dTT (s; zero, not reported by the source)")
    addSd_ECT <- fixed(0)
    label("Additive residual SD on ECT (s; zero, not reported by the source)")
    addSd_TT <- fixed(0)
    label("Additive residual SD on TT (s; zero, not reported by the source)")
  })

  model({
    # 1. Unbound sum dabigatran plasma concentration (umol/L), supplied by
    #    the user. The effect is direct (no effect compartment): Fuhr 2020
    #    links the coagulation assays to the instantaneous plasma
    #    concentration ('The PD effect of dabigatran is directly correlated
    #    to its plasma concentration', Introduction).
    conc <- CU_DABIGATRAN_UM

    # 2. Typical-value parameters (no IIV reported).
    rbase_aPTT <- exp(lrbase_aPTT)
    emax_aPTT <- exp(lemax_aPTT)
    ec50_aPTT <- exp(lec50_aPTT)
    slope_aPTT <- exp(lslope_aPTT)
    rbase_dTT <- exp(lrbase_dTT)
    slope_dTT <- exp(lslope_dTT)
    rbase_ECT <- exp(lrbase_ECT)
    slope_ECT <- exp(lslope_ECT)
    rbase_TT <- exp(lrbase_TT)
    emax_TT <- exp(lemax_TT)
    ec50_TT <- exp(lec50_TT)
    slope_TT <- exp(lslope_TT)

    # 3. Coagulation times (s), Fuhr 2020 ESM Equations 16-19.
    aPTT <- emax_aPTT * conc / (ec50_aPTT + conc) + slope_aPTT * conc + rbase_aPTT
    dTT <- slope_dTT * conc + rbase_dTT
    ECT <- slope_ECT * conc + rbase_ECT
    TT <- emax_TT * conc / (ec50_TT + conc) + slope_TT * conc + rbase_TT

    # 4. Observations. Residual SDs are fixed at zero because none is
    #    reported; free them when refitting to individual data.
    aPTT ~ add(addSd_aPTT)
    dTT ~ add(addSd_dTT)
    ECT ~ add(addSd_ECT)
    TT ~ add(addSd_TT)
  })
}
