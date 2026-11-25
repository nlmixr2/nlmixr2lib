BravoPadros_2026_odronextamab <- function() {
  description <- "Two-compartment population PK model for odronextamab (CD20xCD3 bispecific IgG4 T-cell engager) after IV infusion, with parallel linear and Michaelis-Menten elimination whose Vmax declines exponentially over time toward an asymptote, in adults with relapsed or refractory B-cell non-Hodgkin lymphoma (Bravo Padros 2026)"
  reference <- "Bravo Padros M, Conrado DJ, Srinivasan K, Harnisch LO, Davis JD, Zhu M. Pharmacokinetics of Odronextamab, A Bispecific T-Cell-Engaging Antibody, in Adult Patients With Relapsed or Refractory B-Cell Non-Hodgkin Lymphoma. CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70162. doi:10.1002/psp4.70162"
  vignette <- "BravoPadros_2026_odronextamab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL and (one shared exponent) on Vc and Vp, reference 72.6 kg (Appendix S1 control stream MED_WT = 72.6, the cohort median). Missing values were set to the median in the source analysis.",
      source_name = "WGTBL"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Fractional change (1 + theta) on Vc and Vp (Appendix S1: IF(SEXN.EQ.2) V_FEMALE = LOG(1 + THETA(25))). SEXN = 2 is female, so SEXF = as.integer(SEXN == 2).",
      source_name = "SEXN"
    ),
    ALB = list(
      description = "Serum albumin, time-varying",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying albumin; power effect on linear CL with reference 38.7 g/L (Appendix S1: CL_ALBTIM * (LOG(ALBTIM) - LOG(MED_ALB)), MED_ALB = 38.7). The only time-varying covariate of the model; supply ALB = ALB_BASE for a baseline-only simulation.",
      source_name = "ALB"
    ),
    ALB_BASE = list(
      description = "Baseline serum albumin, time-fixed",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (one shared exponent) on Vc and Vp with reference 38.7 g/L (Appendix S1: V_ALB * (LOG(ALBU) - LOG(MED_ALB)), ALBU = ALBBL). The paper keeps baseline albumin (on the volumes) and time-varying albumin (on CL) as separate data columns.",
      source_name = "ALBBL"
    ),
    IGG = list(
      description = "Baseline serum immunoglobulin G",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on linear CL with reference 38.2 umol/L (Appendix S1 MED_IGG = 38.2 umol/L). The paper reports IgG in umol/L; model() converts the canonical g/L column to umol/L assuming an IgG molar mass of 150 kDa (1 g/L = 6.667 umol/L), so the reference corresponds to 5.73 g/L.",
      source_name = "IGGBL"
    ),
    LYMPH_ABS = list(
      description = "Baseline absolute peripheral-blood lymphocyte count",
      units = "cells/uL",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on Vmax,0 and (one shared exponent) on Vc and Vp, reference 0.792 x 10^9/L = 792 cells/uL (Appendix S1 MED_LYM = 0.792 10^9/L). model() converts the canonical cells/uL column to 10^9/L.",
      source_name = "LYMBL"
    ),
    IL10 = list(
      description = "Baseline serum interleukin-10",
      units = "pg/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline only; power effect on linear CL with reference 3.01 ng/L (= pg/mL) (Appendix S1: CL_IL10 * (LOG(IL10) - LOG(MED_IL10)), MED_IL10 = 3.01). Missing for 5.3% of patients (Table 1), who were set to the median in the source analysis.",
      source_name = "IL10BL"
    ),
    IFNG = list(
      description = "Baseline serum interferon gamma",
      units = "pg/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline only; power effect on linear CL with reference 3.66 ng/L (= pg/mL) (Appendix S1 MED_IFNG = 3.66). Missing for 5.3% of patients (Table 1), who were set to the median in the source analysis.",
      source_name = "IFNGBL"
    ),
    TUMSZ = list(
      description = "Baseline tumour size as the sum of products of perpendicular diameters (SPPD)",
      units = "mm^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on Vmax,0 with reference 2930 mm^2 (Appendix S1 MED_TUMSIZE = 2930, the cohort median; Table 1).",
      source_name = "SPDTLBL"
    ),
    PRIOR_CART = list(
      description = "Prior CAR T cell therapy indicator, 1 = prior CAR T",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Fractional change (1 + theta) on Vmax,0 (Appendix S1: IF(PRCARTN.EQ.1) VMAX_PRIORCART = LOG(1 + THETA(16))). All prior-CAR-T patients were DLBCL CAR T failures in ELM-1 (Table 1).",
      source_name = "PRCARTN"
    ),
    TUMTP_FL = list(
      description = "Follicular lymphoma indicator, 1 = FL",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Fractional changes (1 + theta) on Vmax,0, Kv and Rv (Appendix S1: NHLTYPEN.EQ.3). Reference histology (TUMTP_FL = 0 and TUMTP_OTHER_NHL = 0) is DLBCL, CAR T naive or CAR T failure. Mutually exclusive with TUMTP_OTHER_NHL.",
      source_name = "NHLTYPEN"
    ),
    TUMTP_OTHER_NHL = list(
      description = "Other B-NHL histology indicator, 1 = mantle cell, marginal zone or other B-NHL",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Fractional changes (1 + theta) on Kv and Rv (Appendix S1: NHLTYPEN.GT.3, i.e. MCL, MZL and 'Other B-NHL' pooled; 80 of 507 patients, Table 1). The effect of this group on Vmax,0 was dropped from the final model. Reference histology is DLBCL.",
      source_name = "NHLTYPEN"
    ),
    ANN_ARBOR_LT4 = list(
      description = "Ann Arbor stage below IV at study entry, 1 = stage I, II or III",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Fractional change (1 + theta) on Vmax,0 (Appendix S1: IF(AASTGN.NE.4) VMAX_AASTG123 = LOG(1 + THETA(24))). Reference is stage IV (61.5% of patients, Table 1); ANN_ARBOR_LT4 = as.integer(AASTGN != 4).",
      source_name = "AASTGN"
    )
  )

  compartmentData <- list(
    central = list(analyte = "odronextamab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "odronextamab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 507L,
    n_studies = 2L,
    n_observations = 14618L,
    age_range = "22-89 years",
    age_median = "65 years",
    weight_range = "31.8-165 kg",
    weight_median = "72.6 kg",
    sex_female_pct = 38.3,
    race_ethnicity = c(White = 62),
    disease_state = "Relapsed or refractory CD20-positive B-cell non-Hodgkin lymphoma: DLBCL CAR T naive 39.8%, DLBCL after CAR T failure 11.6%, follicular lymphoma 32.7%, mantle cell lymphoma 4.9%, marginal zone lymphoma 4.9%, other B-NHL 5.9% (Table 1).",
    dose_range = "0.03-320 mg IV infusion over 1-4 h; step-up doses then QW treatment doses, followed in expansion and Phase II by Q2W (and optionally Q4W) maintenance, e.g. 0.7/4/20/160 mg QW then 320 mg Q2W (Table 1).",
    regions = "Global (ELM-1 and ELM-2 multinational trials)",
    notes = "Phase I ELM-1 (NCT02290951, n = 167) and Phase II ELM-2 (NCT03888105, n = 340). Ann Arbor stage IV 61.5%; prior CAR T 13.0%; median baseline albumin 38.7 g/L, IgG 38.2 umol/L, lymphocytes 0.792 x 10^9/L, IL-10 3.02 ng/L, IFN gamma 3.66 ng/L, tumour SPPD 2930 mm^2 (Table 1). Age and weight ranges from the quantile bins of Table 3."
  )

  ini({
    # Fit-for-purpose model, Table 2. Typical values are for the reference
    # patient: 72.6-kg male with DLBCL, CAR T naive, Ann Arbor stage IV,
    # tumour size 2930 mm^2, albumin 38.7 g/L, IgG 38.2 umol/L, IL-10
    # 3.01 ng/L, IFN gamma 3.66 ng/L, lymphocytes 0.792 x 10^9/L.
    lcl <- log(0.189); label("Linear clearance CL (L/day)") # Table 2: CL = 0.189 L/day
    lq <- log(1.21); label("Intercompartmental clearance Q (L/day)") # Table 2: Q = 1.21 L/day
    lvc <- log(4.99); label("Central volume of distribution Vc (L)") # Table 2: Vc = 4.99 L
    lvp <- log(4.42); label("Peripheral volume of distribution Vp (L)") # Table 2: Vp = 4.42 L
    lvmax <- log(2.93); label("Baseline Michaelis-Menten maximum elimination velocity Vmax,0 (mg/L/day)") # Table 2: Vmax,0 = 2.93 mg/L/day
    lkm <- log(1.95); label("Concentration at half of Vmax km (mg/L)") # Table 2: km = 1.95 mg/L
    lvmax_exp_kdes <- log(1.09); label("Rate constant of the decline in Vmax Kv (1/day)") # Table 2: Kv = 1.09 1/day
    lvmax_exp_finf <- log(0.262); label("Ratio of asymptotic to baseline Vmax Rv (unitless)") # Table 2: Rv = 0.262

    # Covariate effects (Table 2); equation forms from the Appendix S1
    # control stream: continuous covariates enter as power functions of
    # covariate / median, categorical ones as fractional changes (1 + theta).
    e_wt_cl <- 0.854; label("Power exponent of baseline weight on CL (unitless)") # Table 2: body weight effect on CL = 0.854
    e_alb_cl <- -0.933; label("Power exponent of time-varying albumin on CL (unitless)") # Table 2: albumin effect on CL = -0.933
    e_igg_cl <- 0.118; label("Power exponent of baseline IgG on CL (unitless)") # Table 2: IgG effect on CL = 0.118
    e_ifng_cl <- 0.0274; label("Power exponent of baseline IFN gamma on CL (unitless)") # Table 2: IFN gamma effect on CL = 0.0274
    e_il10_cl <- -0.0717; label("Power exponent of baseline IL-10 on CL (unitless)") # Table 2: IL-10 effect on CL = -0.0717
    e_wt_vc_vp <- 0.372; label("Power exponent of baseline weight on Vc and Vp (unitless)") # Table 2: body weight effect on Vc and Vp = 0.372
    e_alb_base_vc_vp <- -0.561; label("Power exponent of baseline albumin on Vc and Vp (unitless)") # Table 2: albumin effect on Vc and Vp = -0.561
    e_sexf_vc_vp <- -0.186; label("Fractional change in Vc and Vp for female patients (unitless)") # Table 2: female sex effect on Vc and Vp = -0.186
    e_lymph_abs_vc_vp <- 0.0895; label("Power exponent of baseline lymphocyte count on Vc and Vp (unitless)") # Table 2: lymphocyte count effect on Vc and Vp = 0.0895
    e_ann_arbor_lt4_vmax <- -0.285; label("Fractional change in Vmax,0 for Ann Arbor stage I-III (unitless)") # Table 2: Ann Arbor stage < 4 effect on Vmax,0 = -0.285
    e_lymph_abs_vmax <- 0.313; label("Power exponent of baseline lymphocyte count on Vmax,0 (unitless)") # Table 2: lymphocyte count effect on Vmax,0 = 0.313
    e_tumtp_fl_vmax <- -0.252; label("Fractional change in Vmax,0 for follicular lymphoma (unitless)") # Table 2: FL effect on Vmax,0 = -0.252
    e_prior_cart_vmax <- -0.201; label("Fractional change in Vmax,0 after prior CAR T therapy (unitless)") # Table 2: prior CAR T effect on Vmax,0 = -0.201
    e_tumsz_vmax <- 0.0601; label("Power exponent of baseline tumour size on Vmax,0 (unitless)") # Table 2: tumor size effect on Vmax,0 = 0.0601
    e_tumtp_fl_vmax_exp_kdes <- -0.891; label("Fractional change in Kv for follicular lymphoma (unitless)") # Table 2: FL effect on Kv = -0.891
    e_tumtp_other_nhl_vmax_exp_kdes <- -0.844; label("Fractional change in Kv for other B-NHL (unitless)") # Table 2: other B-NHL type effect on Kv = -0.844
    e_tumtp_fl_vmax_exp_finf <- -0.133; label("Fractional change in Rv for follicular lymphoma (unitless)") # Table 2: FL effect on Rv = -0.133
    e_tumtp_other_nhl_vmax_exp_finf <- -0.292; label("Fractional change in Rv for other B-NHL (unitless)") # Table 2: other B-NHL type effect on Rv = -0.292

    # Between-subject variability on Vc, Vmax,0 and Vp only (diagonal
    # $OMEGA, Appendix S1); none on CL (Results text).
    etalvc ~ 0.0983 # Table 2: BSV variance Vc = 0.0983
    etalvmax ~ 0.403 # Table 2: BSV variance Vmax,0 = 0.403
    etalvp ~ 0.709 # Table 2: BSV variance Vp = 0.709

    # Residual error: additive on log-transformed concentrations
    # (Appendix S1 $ERROR: Y = LOG(IPRED) + EPS(1)); Table 2 variance 0.136.
    expSd <- sqrt(0.136); label("Additive residual error on the log scale (SD)") # Table 2: residual variance = 0.136 (ln mg/L)
  })
  model({
    # Unit conversions of the canonical covariate columns to the units the
    # source medians use: lymphocytes cells/uL -> 10^9/L; IgG g/L -> umol/L
    # at an assumed IgG molar mass of 150 kDa.
    lymph_e9 <- LYMPH_ABS / 1000
    igg_umol <- IGG * 1e6 / 150000

    cl <- exp(lcl) *
      (WT / 72.6)^e_wt_cl *
      (ALB / 38.7)^e_alb_cl *
      (igg_umol / 38.2)^e_igg_cl *
      (IFNG / 3.66)^e_ifng_cl *
      (IL10 / 3.01)^e_il10_cl

    cov_v <- (WT / 72.6)^e_wt_vc_vp *
      (ALB_BASE / 38.7)^e_alb_base_vc_vp *
      (lymph_e9 / 0.792)^e_lymph_abs_vc_vp *
      (1 + e_sexf_vc_vp * SEXF)
    vc <- exp(lvc + etalvc) * cov_v
    vp <- exp(lvp + etalvp) * cov_v
    q <- exp(lq)

    vmax <- exp(lvmax + etalvmax) *
      (lymph_e9 / 0.792)^e_lymph_abs_vmax *
      (TUMSZ / 2930)^e_tumsz_vmax *
      (1 + e_prior_cart_vmax * PRIOR_CART) *
      (1 + e_tumtp_fl_vmax * TUMTP_FL) *
      (1 + e_ann_arbor_lt4_vmax * ANN_ARBOR_LT4)
    km <- exp(lkm)
    vmax_exp_kdes <- exp(lvmax_exp_kdes) *
      (1 + e_tumtp_fl_vmax_exp_kdes * TUMTP_FL) *
      (1 + e_tumtp_other_nhl_vmax_exp_kdes * TUMTP_OTHER_NHL)
    vmax_exp_finf <- exp(lvmax_exp_finf) *
      (1 + e_tumtp_fl_vmax_exp_finf * TUMTP_FL) *
      (1 + e_tumtp_other_nhl_vmax_exp_finf * TUMTP_OTHER_NHL)

    # Eq. 3, dVmax/dt = -Kv * Vmax * (1 - Vmax,asym / Vmax) with
    # Vmax(0) = Vmax,0 and Vmax,asym = Rv * Vmax,0, solved in closed form;
    # t is time since the first dose.
    vmax_t <- vmax * (vmax_exp_finf + (1 - vmax_exp_finf) * exp(-vmax_exp_kdes * t))

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    Cc <- central / vc
    # Eq. 1: KMM = Vmax / (km + Acentral / Vc)
    kel_mm <- vmax_t / (km + Cc)

    d/dt(central) <- -kel_mm * central - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> mg/L.
    Cc ~ lnorm(expSd)
  })
}
