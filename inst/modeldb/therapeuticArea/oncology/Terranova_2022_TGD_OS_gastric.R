Terranova_2022_TGD_OS_gastric <- function() {
  description <- paste(
    "Joint tumor growth dynamics (TGD) + overall survival (OS)",
    "disease-progression model for advanced gastric cancer /",
    "gastroesophageal junction cancer (GC/GEJC) developed on the JAVELIN",
    "Gastric 100 phase III trial by Terranova et al. (Merck KGaA /",
    "Merck Institute of Pharmacometrics / EMD Serono / Metrum Research",
    "Group). N = 499 patients randomized 1:1 to avelumab 10 mg/kg every",
    "2 weeks (n = 249) or continued chemotherapy (n = 250) maintenance",
    "after 12 weeks of induction chemotherapy. The paper's primary",
    "finding is NEGATIVE: no treatment effect on OS or TGD was",
    "identified, consistent with the primary JAVELIN Gastric 100",
    "analysis. Machine learning (random forests + SIDEScreen + Boruta /",
    "permutation / random-splits variable-importance methods) was used",
    "to screen 89 (OS) / 52 (TGD) candidate covariates, then covariates",
    "surviving the ML screen were incorporated into two parametric",
    "sub-models: (1) a Gompertzian TGD model on tumor size measured",
    "by sum of longest diameters (SLD, mm) with the Vaghi 2020",
    "Kd = slope*Kg + intercept reduced-parameter reformulation, and",
    "(2) a log-logistic parametric time-to-event OS model with an",
    "accelerated-failure-time link to time-invariant covariates on the",
    "log-median-OS scale, a proportional-hazards link to time-varying",
    "laboratory covariates, and a per-arm shape parameter that captures",
    "the crossing of hazard curves between avelumab and chemotherapy.",
    "There is NO drug PK sub-model -- avelumab exposure is not carried",
    "explicitly; the treatment arm is encoded as the TRT categorical",
    "covariate (1 = chemotherapy = reference, 2 = avelumab).",
    "Extract-scope note: no drug PK; the two ODE outputs are `tumor`",
    "(SLD, mm) with proportional residual error and `sur` (survival",
    "probability), driven off `d/dt(cumhaz) = hazard`."
  )
  reference <- paste(
    "Terranova N, French J, Dai H, Wiens M, Khandelwal A, Ruiz-Garcia A,",
    "Manitz J, von Heydebreck A, Ruisi M, Chin K, Girard P,",
    "Venkatakrishnan K.",
    "Pharmacometric modeling and machine learning analyses of prognostic",
    "and predictive factors in the JAVELIN Gastric 100 phase III trial",
    "of avelumab.",
    "CPT Pharmacometrics Syst Pharmacol. 2022;11(3):333-347.",
    "doi:10.1002/psp4.12754.",
    sep = " "
  )
  vignette <- "Terranova_2022_TGD_OS_gastric"

  units <- list(
    time          = "days (both sub-models run on a single days time axis; the paper reports Kg / Kd in 1/year but converts to 1/day inside the model equations, and the OS TTE was fit with time in days)",
    dosing        = "n/a (no PK sub-model; avelumab dose 10 mg/kg every 2 weeks and chemotherapy regimen are encoded implicitly via the TRT covariate)",
    concentration = "mm (TGD observable = sum of longest diameters of RECIST 1.1 target lesions); probability (OS sub-model output `sur`)"
  )

  covariateData <- list(
    TRT = list(
      description = "Per-subject treatment-arm integer indicator selecting the per-arm shape parameter of the log-logistic OS hazard and enabling an ML-screened but ultimately not-detected treatment effect on TGD.",
      units       = "(categorical / integer-coded)",
      type        = "categorical",
      source_name = "treatment arm",
      notes       = paste(
        "Integer coding (Terranova 2022 Methods 'JAVELIN Gastric 100'):",
        "  1 = CHEMO (continued chemotherapy; oxaliplatin + fluoropyrimidine or best supportive care only if ineligible; n = 250) -- reference arm.",
        "  2 = AVEL (avelumab 10 mg/kg IV every 2 weeks; n = 249).",
        "Reference = 1 (chemotherapy). The per-arm log-logistic shape parameters (2.25 chemo, 1.63 avelumab per Table S2) are selected via the (TRT == 2) indicator inside model(); the shape difference was the only credible-interval-excluding-zero treatment effect. The TRT effect on log-median OS was estimated but was consistent with the null (median-OS HR 1.10, 95% CrI 0.935-1.32 per Table S2). No treatment effect on TGD (Kg or Kd) was retained after ML screening.",
        sep = "\n"
      )
    ),

    T_DIAG_CANCER = list(
      description = "Time since primary gastric-cancer / GEJC diagnosis at re-baseline (randomization).",
      units       = "days",
      type        = "continuous",
      source_name = "Tdiag",
      notes       = "Power effect on the Gompertz tumor-growth rate constant Kg: `Kg = tvKg * (T_DIAG_CANCER / 53)^-0.00291`. Reference 53 days = avelumab-arm median at randomization (Terranova 2022 Supplementary Methods). Additive-log effect on log-median OS (Table S2): `log_median_OS += 0.0436 * log(T_DIAG_CANCER)`. Broad range (both arms combined: 3 to 9080 days; median 48; Table S3). Supply at least a small floor value (>= 1 day) for newly-diagnosed subjects in simulation to keep the power form well-defined."
    ),

    LMET = list(
      description = "Baseline liver-metastasis indicator (radiologically or biopsy-documented) at re-baseline.",
      units       = "(binary)",
      type        = "binary",
      source_name = "LIVER",
      notes       = "Exponential effect on Kg via a multiplicative (1 + LIVER*0.0164) term in the Terranova 2022 Supplementary Methods equation. Prevalence 32.7% pooled (Table S6)."
    ),

    MET_GE3 = list(
      description = "Baseline number of metastatic sites >= 3 indicator (derived from the source NONMEM continuous count column NumMet).",
      units       = "(binary)",
      type        = "binary",
      source_name = "NumMet",
      notes       = paste(
        "Derived via `MET_GE3 = as.integer(NumMet >= 3)` at data-assembly time.",
        "Terranova 2022 encodes the source paper's continuous count `NumMet` as `exp(0.143 * (NumMet - 3))` on baseline tumor size BASE (Supplementary Methods).",
        "Per the nlmixr2lib count-covariate-decomposed-to-binary policy, the continuous count is binarised at the paper's reference value of 3.",
        "The per-metastasis coefficient 0.143 is reused as the single-step exponential effect: `BASE = ... * exp(0.143 * MET_GE3)`.",
        "Deviation from the paper's continuous form documented in vignette Errata.",
        "Cohort NumMet: mean 2.51, median 2, range 0-12 (both arms combined; Table S3).",
        sep = " "
      )
    ),

    RESP_SD = list(
      description = "Best RECIST 1.1 response at re-baseline is stable disease (SD).",
      units       = "(binary)",
      type        = "binary",
      source_name = "RES_SD",
      notes       = paste(
        "Coded 1 = SD, 0 = any other RECIST category (CR, PR, PD, NED, NE, non-CR/non-PD) at re-baseline.",
        "Enters BOTH sub-models: (a) TGD BASE covariate equation (Terranova 2022 Supplementary Methods): `BASE = ... * (1 + 0.644 * RESP_SD)`; (b) OS TTE (Table S2): `log_median_OS += -0.0919 * RESP_SD`.",
        "Prevalence 34.3% pooled at re-baseline (Table S6).",
        sep = " "
      )
    ),

    RESP_NONPDCR = list(
      description = "Best RECIST 1.1 response at re-baseline is NEITHER complete response (CR) NOR partial response (PR) -- i.e., non-responder by RECIST.",
      units       = "(binary)",
      type        = "binary",
      source_name = "RES_nonPR/nonCR",
      notes       = paste(
        "Coded 1 = not CR and not PR (any of SD, non-CR/non-PD, NED, PD, NE), 0 = responder (CR or PR).",
        "Enters TGD BASE only (Terranova 2022 Supplementary Methods): `BASE = ... * (1 - 0.0769 * RESP_NONPDCR)`.",
        "Distinct from `RESP_SD` (proper subset when both indicators = 1 would be a coding inconsistency) and `RESP_RESPONDER` (which uses the opposite reference-category framing).",
        sep = " "
      )
    ),

    RESP_RESPONDER = list(
      description = "Best RECIST 1.1 response at re-baseline is responder (CR or PR).",
      units       = "(binary)",
      type        = "binary",
      source_name = "Re-baseline responder vs other",
      notes       = paste(
        "Coded 1 = CR or PR at re-baseline, 0 = non-responder (any of SD, non-CR/non-PD, NED, PD, NE).",
        "Enters OS TTE only (Table S2): `log_median_OS += 0.146 * RESP_RESPONDER`.",
        "Complement of RESP_NONPDCR with a different reference category; both are needed because TGD and OS use different reference-category framings.",
        "Prevalence ~42.6% pooled (Table S6: 40.8% PR + 1.8% CR).",
        sep = " "
      )
    ),

    AGE = list(
      description = "Baseline age in years.",
      units       = "years",
      type        = "continuous",
      source_name = "AGE",
      notes       = "Additive-log effect on log-median OS (Table S2): `log_median_OS += 0.418 * log(AGE)`. Older age is prognostic for longer OS. Cohort AGE: mean 60.6, median 62, range 21-88 (both arms combined; Table S3)."
    ),

    HR = list(
      description = "Baseline heart rate (beats per minute).",
      units       = "beats/min",
      type        = "continuous",
      source_name = "HR",
      notes       = "Additive-log effect on log-median OS (Table S2): `log_median_OS += 0.0879 * log(HR)`. Time-invariant re-baseline value used. Cohort HR: mean 77.6, median 76, range 42-128 (both arms combined; Table S3)."
    ),

    SBP = list(
      description = "Baseline systolic blood pressure (mmHg).",
      units       = "mmHg",
      type        = "continuous",
      source_name = "SBP",
      notes       = "Additive-log effect on log-median OS (Table S2): `log_median_OS += 0.532 * log(SBP)`. Cohort SBP: mean 122, median 120, range 82-183 (both arms combined; Table S3). Time-invariant re-baseline value used."
    ),

    ALB = list(
      description = "Serum albumin concentration.",
      units       = "g/dL",
      type        = "continuous",
      source_name = "ALB",
      notes       = paste(
        "Enters OS TTE as a TIME-VARYING proportional-hazards multiplier (Table S2): `hazard *= exp(-1.26 * ALB)`.",
        "Terranova 2022 uses ALB in g/dL (Table S3 cohort mean 40.8 g/dL is g/L numerically but the units column of Table S3 reports 'g/dL'; treated here as g/dL per the paper's stated units).",
        "Cohort ALB (both arms): mean 40.8, median 41, range 27-49 (units per Table S3).",
        sep = " "
      )
    ),

    ALP = list(
      description = "Alkaline phosphatase activity.",
      units       = "IU/L",
      type        = "continuous",
      source_name = "ALP",
      notes       = "Enters OS TTE as a TIME-VARYING proportional-hazards multiplier (Table S2): `hazard *= exp(0.0364 * ALP)`. Cohort ALP: mean 115, median 87, range 31-1540 (both arms combined; Table S3)."
    ),

    AST = list(
      description = "Aspartate aminotransferase activity.",
      units       = "IU/L",
      type        = "continuous",
      source_name = "AST",
      notes       = "Enters OS TTE as a TIME-VARYING proportional-hazards multiplier (Table S2): `hazard *= exp(0.165 * AST)`. Time-varying re-baseline value used."
    ),

    CRP = list(
      description = "C-reactive protein.",
      units       = "mg/L",
      type        = "continuous",
      source_name = "CRP",
      notes       = "Enters OS TTE as a TIME-VARYING proportional-hazards multiplier (Table S2): `hazard *= exp(0.258 * CRP)`. Largest time-varying effect in the OS model: 5th-95th-percentile CRP change (0.5-34.9 mg/L) gives hazard ratio 0.699-2.01 at cohort median 2.0 mg/L."
    ),

    LDH = list(
      description = "Serum lactate dehydrogenase activity.",
      units       = "IU/L",
      type        = "continuous",
      source_name = "LDH",
      notes       = "Enters OS TTE as a TIME-VARYING proportional-hazards multiplier (Table S2): `hazard *= exp(0.618 * LDH)`. Cohort LDH: mean 210, median 175, range 83-1480 (both arms combined; Table S3)."
    ),

    NLR = list(
      description = "Neutrophil-to-lymphocyte ratio.",
      units       = "ratio (unitless)",
      type        = "continuous",
      source_name = "NLR",
      notes       = "Enters OS TTE as a TIME-VARYING proportional-hazards multiplier (Table S2): `hazard *= exp(0.339 * NLR)`. Cohort NLR: mean 4.13, median 3.30, range 0.667-32 (both arms combined; Table S3). Only NLR was retained among four highly correlated time-varying inflammation indices."
    ),

    CPK = list(
      description = "Serum creatine phosphokinase (creatine kinase).",
      units       = "IU/L",
      type        = "continuous",
      source_name = "CK",
      notes       = "Additive-log effect on log-median OS (Table S2): `log_median_OS += 0.0968 * log(CPK)`. Cohort CPK: mean 72.3, median 60.5, range 11-567 (both arms combined; Table S3). Time-invariant re-baseline value used."
    ),

    CRCL = list(
      description = "Estimated glomerular filtration rate (eGFR); reported in Table S3 as mL/min/1.73 m^2 per unit conventions of the JAVELIN Gastric 100 clinical database (values 0.412-3.08 in Table S3 suggest the raw units are mL/s/1.73 m^2 or scaled). Encoded here in the paper's reported units.",
      units       = "mL/min/1.73 m^2 (units as reported; document any conversion applied at data-assembly time)",
      type        = "continuous",
      source_name = "eGFR",
      notes       = "Additive-log effect on log-median OS (Table S2): `log_median_OS += 0.237 * log(CRCL)`. Time-invariant re-baseline value used. Table S3 reports eGFR values in an unusual numeric range (0.412-3.08), suggesting they may be pre-normalised; document the unit conversion applied at simulation time in the vignette Errata section."
    ),

    PRIOR_GAST = list(
      description = "Prior gastrectomy (surgical removal of part or all of the stomach) indicator at baseline.",
      units       = "(binary)",
      type        = "binary",
      source_name = "PRIOR_GAST",
      notes       = "Additive effect on log-median OS (Table S2): `log_median_OS += -0.0313 * PRIOR_GAST`. Prevalence 27.1% pooled (Table S6)."
    ),

    GGT = list(
      description = "Gamma-glutamyl transferase activity.",
      units       = "IU/L",
      type        = "continuous",
      source_name = "GGT",
      notes       = paste(
        "Additive-log effect on log-median OS (Table S2): `log_median_OS += 0.0968 * log(GGT)`.",
        "Higher GGT paradoxically prognostic for LONGER OS in this cohort (contrary to prior gastric-cancer literature); flagged as one of the largest time-invariant covariate effects (17% longer OS at 95th percentile vs median).",
        "Cohort GGT: mean 62.4, median 25, range 4-1070 (both arms combined; Table S3).",
        sep = " "
      )
    ),

    PERIT_CARC = list(
      description = "Peritoneal carcinomatosis (diffuse peritoneal-surface metastatic spread) at baseline / re-baseline.",
      units       = "(binary)",
      type        = "binary",
      source_name = "Peritoneal carcinomatosis",
      notes       = "Additive effect on log-median OS (Table S2): `log_median_OS += -0.181 * PERIT_CARC`. One of the meaningful effects (credible interval excludes zero and the posterior median exceeds the paper's +/-15% threshold). Prevalence 29.9% pooled (Table S6)."
    ),

    TUM_SLD = list(
      description = "Baseline tumor size measured as the sum of longest diameters (SLD) of RECIST 1.1 target lesions.",
      units       = "mm",
      type        = "continuous",
      source_name = "SLD baseline",
      notes       = "Additive effect on log-median OS (Table S2): `log_median_OS += 0.00102 * TUM_SLD`. Cohort baseline SLD: mean 43.5, median 32, range 0-310 mm (both arms combined; Table S3). Time-invariant baseline value used in the OS sub-model; the TGD sub-model separately estimates the individual BASE parameter as a random effect (not the same as this covariate)."
    ),

    TRIG = list(
      description = "Serum triglyceride concentration.",
      units       = "mg/dL",
      type        = "continuous",
      source_name = "TRIG",
      notes       = "Additive-log effect on log-median OS (Table S2): `log_median_OS += -0.0151 * log(TRIG)`. Time-invariant re-baseline value used."
    ),

    WHO_PS = list(
      description = "Eastern Cooperative Oncology Group (ECOG) performance status at re-baseline. In this Terranova 2022 model the covariate enters as a binary contrast (ECOG >= 1 vs ECOG 0) rather than as the raw ordinal 0/1/2 integer.",
      units       = "(integer score with per-model binarisation; document via covariateData[[WHO_PS]]$notes)",
      type        = "categorical",
      source_name = "Re-baseline ECOG PS: >=1 vs 0",
      notes       = paste(
        "The source model binarises the ordinal WHO_PS to ECOG >= 1 (vs ECOG 0 reference) for the OS TTE covariate list.",
        "Derive an intermediate `ECOG_GE1 = as.integer(WHO_PS >= 1)` at data-assembly time (paper reports 40.2% ECOG 0, 58.8% ECOG 1, 0.6% ECOG 2; Table S6).",
        "Additive effect on log-median OS (Table S2): `log_median_OS += -0.155 * ECOG_GE1`.",
        sep = " "
      )
    )
  )

  population <- list(
    species         = "human (adults with unresectable, HER2-negative, locally advanced or metastatic gastric cancer / gastroesophageal junction cancer who did not progress after 12 weeks of induction chemotherapy with oxaliplatin + a fluoropyrimidine)",
    n_subjects      = 499L,
    n_studies       = 1L,
    age_range       = "21-88 years (median 62, mean 60.6; Table S3 both-arms row)",
    weight_range    = "not reported at the pooled level in this paper (weight was included in the ML covariate screen per Table S1 but was not retained in either final model)",
    sex_female_pct  = NA_real_,
    race_ethnicity  = "Asian (Japan / Republic of Korea / Taiwan / Thailand) vs non-Asian; per-region distribution not tabulated at the pooled level. Asian-vs-non-Asian was a stratification factor in the JAVELIN Gastric 100 randomization; the paper explicitly reports Asian vs non-Asian region was NOT identified as a covariate in either OS or TGD models (Figures S3 and S8).",
    disease_state   = "advanced (unresectable, locally advanced or metastatic) HER2-negative gastric cancer / gastroesophageal junction cancer; all patients received 12 weeks of first-line induction chemotherapy with oxaliplatin + a fluoropyrimidine before randomization and had no progressive disease at re-baseline (complete response, partial response, or stable disease at re-baseline required for enrollment).",
    dose_range      = "n/a (no drug PK sub-model). Avelumab dose was 10 mg/kg IV every 2 weeks in the maintenance phase; chemotherapy dose per the induction regimen (or best supportive care only if ineligible for further chemotherapy). Encoded implicitly via the TRT covariate.",
    regions         = "multiregional phase III trial (JAVELIN Gastric 100), including Asian and non-Asian regions as a randomization stratification factor.",
    notes           = paste(
      "TGD sub-model: Gompertzian ODE `dy/dt = Kg*y - Kd*y*log(y)` with y(0) = BASE, using the Vaghi 2020 reformulation Kd = slope*Kg + intercept so a single random effect (etalKg) drives both Kg and Kd with perfect correlation.",
      "  Kg  (year^-1): typical value 0.271; power form `(T_DIAG_CANCER/53)^-0.00291 * (1 + LIVER*0.0164)`; IIV etalKg ~ N(0, 0.00163).",
      "  Kd  (1/(year*mm)): derived as `-18.8*Kg + 5.18`; IIV inherited fully from etalKg.",
      "  BASE (mm): typical value 27.0; multiplicative form `exp(0.143*(MET_GE3-0)) * (1 - 0.0769*RESP_NONPDCR) * (1 + 0.644*RESP_SD)`; IIV etalBase ~ N(0, 0.659); prior binarisation deviation documented in Errata.",
      "  Residual: proportional variance 0.0505 -> propSd = sqrt(0.0505) ~ 0.2247 (CV% = 22.5).",
      "",
      "OS TTE sub-model: log-logistic parametric hazard `h(t) = alpha * lam^alpha * t^(alpha-1) / (1 + (lam*t)^alpha)` with alpha = per-arm shape (2.25 chemo / 1.63 avelumab) and lam = 1/median_i (days^-1) where median_i on log scale accumulates the 12 time-invariant covariate effects; time-varying covariates enter multiplicatively on the hazard as `exp(sum(z_i(t) * alpha_i))`.",
      "  Chemotherapy-arm median OS = 444 days at reference covariate levels (Table S2).",
      "  Avelumab-arm median-OS HR vs chemo = 1.10 (Table S2; CI 0.935-1.32; consistent with null).",
      "  No IIV or residual error on the survival probability; the sub-model output `sur = exp(-cumhaz)` is a deterministic conditional-mean survival given the covariate trajectory. Bootstrap uncertainty on covariate coefficients is reproduced from Table S2 posterior CIs in the vignette.",
      "",
      "Both sub-models are drug-independent disease-progression models; treatment arm is a covariate, not a PK/PD input.",
      sep = "\n"
    )
  )

  ini({
    # ------------------------ TGD sub-model ------------------------
    # Structural parameters (Table 1 'Final model Estimate' column, doi:10.1002/psp4.12754)
    # All held FIXED because this file distributes the paper's published final estimates rather than re-estimating.
    lKg           <- fixed(log(0.271)); label("log typical Gompertz tumor-growth rate constant Kg (year^-1)")               # Table 1: theta1 = 0.271 year^-1
    lBase         <- fixed(log(27.0));  label("log typical baseline tumor size BASE (mm)")                                    # Table 1: theta2 = 27.0 mm
    slope_kd      <- fixed(-18.8);      label("Slope of the linear relationship Kd = slope*Kg + intercept (1/mm)")            # Table 1: theta3 = -18.8 1/mm
    intercept_kd  <- fixed(5.18);       label("Intercept of the linear relationship Kd = slope*Kg + intercept (1/(year*mm))") # Table 1: theta4 = 5.18 1/(year*mm)

    # Covariate effects on Kg (Terranova 2022 Supplementary Methods)
    e_tdiag_kg    <- fixed(-0.00291); label("Power exponent of (T_DIAG_CANCER / 53) on Kg")   # Table 1: theta5 = -0.00291
    e_liver_kg    <- fixed(0.0164);   label("Linear multiplicative effect of LMET on Kg")     # Table 1: theta6 = 0.0164

    # Covariate effects on BASE (Terranova 2022 Supplementary Methods; paper's continuous NumMet
    # binarised to MET_GE3 per the count-covariate policy -- deviation documented in vignette Errata)
    e_nmet_base   <- fixed(0.143);   label("Exponential effect of MET_GE3 on BASE (single-step; paper's per-metastasis form was exp(0.143*(NumMet-3)))") # Table 1: theta7 = 0.143
    e_rnpdcr_base <- fixed(-0.0769); label("Linear-multiplicative effect of RESP_NONPDCR on BASE")                            # Table 1: theta8 = -0.0769
    e_rsd_base    <- fixed(0.644);   label("Linear-multiplicative effect of RESP_SD on BASE")                                 # Table 1: theta9 = 0.644

    # TGD IIV (Table 1 'Interindividual variance parameters')
    etalKg   ~ 0.00163; # Table 1: Omega(1,1) = 0.00163 (CV% = 4.04)
    etalBase ~ 0.659;   # Table 1: Omega(2,2) = 0.659 (CV% = 96.6)

    # TGD residual error (Table 1 'Residual variance')
    propSd <- fixed(sqrt(0.0505)); label("TGD proportional residual SD (fraction on SLD in mm)") # Table 1: Sigma(1,1) = 0.0505 (CV% = 22.5) -> SD = sqrt(0.0505) ~ 0.2247

    # ------------------------ OS TTE sub-model ------------------------
    # Structural (Table S2 posterior means)
    log_med_chemo    <- fixed(log(444)); label("log typical median OS in the chemotherapy arm at reference covariate levels (log-days)") # Table S2: chemo median OS = 444 days
    e_trt_ave_median <- fixed(log(1.10)); label("Additive effect of TRT_AVE on log-median OS (Table S2 'HR for avelumab vs chemotherapy arms' = 1.10)")

    # Per-arm log-logistic shape (Table S2 'Shape parameter')
    lshape_chemo    <- fixed(log(2.25)); label("log log-logistic shape parameter in the chemotherapy arm (dimensionless)") # Table S2: chemo shape = 2.25
    lshape_ave      <- fixed(log(1.63)); label("log log-logistic shape parameter in the avelumab arm (dimensionless)")     # Table S2: avelumab shape = 1.63

    # Time-invariant covariate effects on log-median OS (Table S2 'Covariate effects on log median OS')
    e_log_age_median      <- fixed(0.418);    label("Additive effect of log(AGE) on log-median OS")            # Table S2: 0.418 (0.0720, 0.768)
    e_ecog_ge1_median     <- fixed(-0.155);   label("Additive effect of ECOG >= 1 (WHO_PS>=1) on log-median OS") # Table S2: -0.155 (-0.320, 0.00573)
    e_perit_carc_median   <- fixed(-0.181);   label("Additive effect of PERIT_CARC on log-median OS")           # Table S2: -0.181 (-0.349, -0.0190)
    e_log_cpk_median      <- fixed(0.0968);   label("Additive effect of log(CPK) on log-median OS")             # Table S2: 0.0968 (-0.0301, 0.229)
    e_log_tdiag_median    <- fixed(0.0436);   label("Additive effect of log(T_DIAG_CANCER) on log-median OS")   # Table S2: 0.0436 (-0.0261, 0.113)
    e_log_crcl_median     <- fixed(0.237);    label("Additive effect of log(CRCL) on log-median OS")            # Table S2: 0.237 (-0.0425, 0.508)
    e_prior_gast_median   <- fixed(-0.0313);  label("Additive effect of PRIOR_GAST on log-median OS")           # Table S2: -0.0313 (-0.292, 0.231)
    e_log_ggt_median      <- fixed(0.0968);   label("Additive effect of log(GGT) on log-median OS")             # Table S2: 0.0968 (-0.000139, 0.197)
    e_log_hr_median       <- fixed(0.0879);   label("Additive effect of log(HR) on log-median OS")              # Table S2: 0.0879 (-0.377, 0.544)
    e_responder_median    <- fixed(0.146);    label("Additive effect of RESP_RESPONDER on log-median OS")       # Table S2: 0.146 (-0.0644, 0.365)
    e_resp_sd_median      <- fixed(-0.0919);  label("Additive effect of RESP_SD on log-median OS")              # Table S2: -0.0919 (-0.307, 0.126)
    e_sld_median          <- fixed(0.00102);  label("Additive effect of TUM_SLD (baseline SLD, mm) on log-median OS") # Table S2: 0.00102 (-0.000500, 0.00251)
    e_log_sbp_median      <- fixed(0.532);    label("Additive effect of log(SBP) on log-median OS")             # Table S2: 0.532 (-0.0195, 1.15)
    e_log_trig_median     <- fixed(-0.0151);  label("Additive effect of log(TRIG) on log-median OS")            # Table S2: -0.0151 (-0.192, 0.170)

    # Time-varying covariate effects on the hazard (Table S2 'Time-varying covariate effects on log HR')
    e_alb_haz   <- fixed(-1.26);   label("Time-varying proportional-hazards coefficient on ALB (log HR per unit ALB)")   # Table S2: -1.26 (-2.09, -0.423)
    e_alp_haz   <- fixed(0.0364);  label("Time-varying proportional-hazards coefficient on ALP (log HR per unit ALP)")   # Table S2: 0.0364 (-0.209, 0.284)
    e_ast_haz   <- fixed(0.165);   label("Time-varying proportional-hazards coefficient on AST (log HR per unit AST)")   # Table S2: 0.165 (-0.0837, 0.418)
    e_crp_haz   <- fixed(0.258);   label("Time-varying proportional-hazards coefficient on CRP (log HR per unit CRP)")   # Table S2: 0.258 (0.170, 0.345)
    e_ldh_haz   <- fixed(0.618);   label("Time-varying proportional-hazards coefficient on LDH (log HR per unit LDH)")   # Table S2: 0.618 (0.331, 0.914)
    e_nlr_haz   <- fixed(0.339);   label("Time-varying proportional-hazards coefficient on NLR (log HR per unit NLR)")   # Table S2: 0.339 (0.184, 0.494)
  })

  model({
    # ------------------------ TGD sub-model ------------------------
    # Individual Kg and Kd on the annual scale (paper reports in year^-1);
    # rate constants are then converted to day^-1 before entering the ODE, since
    # the ODE runs on a days time axis.
    #
    # Terranova 2022 Supplementary Methods:
    #   Kg [/year] = tvKg * (Tdiag / 53)^(-0.00291) * (1 + LIVER * 0.0164)
    #   Kd [/year/mm] = slope * Kg + intercept
    # Both rate constants are then multiplied by (1/365) to give /day and /day/mm.

    Kg_annual <- exp(lKg + etalKg) *
                 (T_DIAG_CANCER / 53)^e_tdiag_kg *
                 (1 + LMET * e_liver_kg)
    Kg        <- Kg_annual / 365   # /day

    Kd_annual <- slope_kd * Kg_annual + intercept_kd   # 1/(year*mm)
    Kd        <- Kd_annual / 365                       # 1/(day*mm)

    # Individual baseline tumor size (mm).
    # Paper's continuous NumMet centered at 3 -> binarised MET_GE3 per policy.
    BASE <- exp(lBase + etalBase) *
            exp(e_nmet_base * MET_GE3) *
            (1 + e_rnpdcr_base * RESP_NONPDCR) *
            (1 + e_rsd_base * RESP_SD)

    # Gompertz ODE for tumor size (SLD in mm)
    # dy/dt = Kg*y - Kd*y*log(y)
    d/dt(tumor) <- Kg * tumor - Kd * tumor * log(tumor)
    tumor(0)    <- BASE

    # ------------------------ OS TTE sub-model ------------------------
    # Individual median OS (log scale) at the reference (chemotherapy) arm.
    # AGE, HR, SBP, T_DIAG_CANCER, CPK, GGT, CRCL, TRIG enter as additive log terms;
    # ECOG_GE1, PERIT_CARC, PRIOR_GAST, RESP_RESPONDER, RESP_SD, TRT_AVE, TUM_SLD
    # enter as additive linear terms.
    trt_ave  <- (TRT == 2)              # 1 if avelumab arm, 0 if chemotherapy (reference)
    ecog_ge1 <- (WHO_PS >= 1)           # binary indicator derived from ordinal WHO_PS

    log_median_os <- log_med_chemo +
                     e_trt_ave_median   * trt_ave +
                     e_log_age_median   * log(AGE)   +
                     e_ecog_ge1_median  * ecog_ge1  +
                     e_perit_carc_median* PERIT_CARC +
                     e_log_cpk_median   * log(CPK)   +
                     e_log_tdiag_median * log(T_DIAG_CANCER) +
                     e_log_crcl_median  * log(CRCL)  +
                     e_prior_gast_median* PRIOR_GAST +
                     e_log_ggt_median   * log(GGT)   +
                     e_log_hr_median    * log(HR)    +
                     e_responder_median * RESP_RESPONDER +
                     e_resp_sd_median   * RESP_SD    +
                     e_sld_median       * TUM_SLD    +
                     e_log_sbp_median   * log(SBP)   +
                     e_log_trig_median  * log(TRIG)

    median_os_i <- exp(log_median_os)   # days
    lam         <- 1 / median_os_i      # days^-1 (log-logistic scale)

    # Per-arm log-logistic shape
    shape <- exp(lshape_chemo + (lshape_ave - lshape_chemo) * trt_ave)

    # Log-logistic baseline hazard with a small time offset to avoid t^(shape-1)
    # singularity at t = 0 (shape < 1 case) -- paper's shape > 1 for both arms,
    # so h_0(0) = 0 analytically, but numerical rendering safer with del_t.
    del_t  <- 1.0e-6
    t_pos  <- t + del_t
    lam_t  <- lam * t_pos
    lam_t_alpha <- lam_t^shape
    h0     <- shape * (lam^shape) * t_pos^(shape - 1) / (1 + lam_t_alpha)

    # Time-varying proportional-hazards multiplier
    hr_tv <- exp(e_alb_haz * ALB +
                 e_alp_haz * ALP +
                 e_ast_haz * AST +
                 e_crp_haz * CRP +
                 e_ldh_haz * LDH +
                 e_nlr_haz * NLR)

    hazard <- h0 * hr_tv

    # Cumulative hazard + survival
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0

    sur <- exp(-cumhaz)

    # ------------------------ Observation model ------------------------
    # Only the TGD `tumor` state carries residual error; the OS `sur` output is
    # a deterministic conditional-mean survival probability given the covariate
    # trajectory (no per-observation residual model).
    tumor ~ prop(propSd)
  })
}
attr(Terranova_2022_TGD_OS_gastric, "message") <-
  "Joint TGD (Gompertzian SLD) + OS (log-logistic parametric TTE) disease-progression model for advanced gastric / GEJC (Terranova 2022, JAVELIN Gastric 100, n=499). No drug PK sub-model. TGD state `tumor` (SLD, mm) with proportional residual error propSd = sqrt(0.0505). OS state `sur` = exp(-cumhaz) driven by d/dt(cumhaz) = hazard, where hazard = log-logistic baseline * time-varying-covariate proportional multiplier. Per-arm shape (2.25 chemo / 1.63 avelumab) via the TRT integer covariate; treatment effect on median OS retained but consistent with null (HR 1.10, 95% CrI 0.935-1.32). Continuous NumMet count binarised to MET_GE3 per the count-covariate-decomposed-to-binary policy; deviation documented in vignette Errata."
Terranova_2022_TGD_OS_gastric
