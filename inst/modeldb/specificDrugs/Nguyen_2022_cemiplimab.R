Nguyen_2022_cemiplimab <- function() {
  description <- "Two-compartment population PK model for cemiplimab (anti-PD-1 IgG4) with zero-order IV infusion, first-order elimination and a time-varying sigmoid-Emax change in clearance, in adults with solid tumors (CSCC, BCC, NSCLC, other advanced malignancies), updated from the earlier cemiplimab model and externally validated in recurrent or metastatic cervical cancer (Nguyen 2022)"
  reference <- "Nguyen J-H, Epling D, Dolphin N, Paccaly A, Conrado D, Davis JD, Al-Huniti N. Population pharmacokinetics modeling and exposure-response analyses of cemiplimab in patients with recurrent or metastatic cervical cancer. CPT Pharmacometrics Syst Pharmacol. 2022;11(11):1458-1471. doi:10.1002/psp4.12855"
  vignette <- "Nguyen_2022_cemiplimab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL and Q (shared exponent CL_WGTBL) and on V1 and V2 (shared exponent VssWGTBL), reference 75 kg (Appendix S1 control stream, M_WGTBL = 75; the cohort median was 73.2 kg). Baseline value only.",
      source_name = "WGTBL"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Fractional change (1 + theta) on CL and on V1 for female patients (Appendix S1: IF (SEXN.EQ.2) CL_SEX = 1 + THETA(15)). SEXN = 2 is female, so SEXF = as.integer(SEXN == 2).",
      source_name = "SEXN"
    ),
    ALB = list(
      description = "Serum albumin, time-varying",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying albumin; power effect on CL with reference 38.6 g/L (Appendix S1: CL_ALB = (ALB/M_ALBBL)**THETA(20), M_ALBBL = 38.6, the cohort median baseline albumin). Supply the per-record albumin; a constant column equal to ALB_BASE reproduces a baseline-only analysis.",
      source_name = "ALB"
    ),
    ALB_BASE = list(
      description = "Baseline serum albumin, time-fixed",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on V1 with reference 38.6 g/L (Appendix S1: V1_ALBBL = (ALBBL/M_ALBBL)**THETA(17)). The paper keeps baseline albumin (on V1) and time-varying albumin (on CL) as separate data columns.",
      source_name = "ALBBL"
    ),
    ALT = list(
      description = "Baseline serum alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline ALT; power effect on CL with reference 19 IU/L (Appendix S1: CL_ALT = (ALTBL/M_ALTBL)**THETA(23), M_ALTBL = 19, the cohort median).",
      source_name = "ALTBL"
    ),
    TUMTP_CSCC = list(
      description = "Cutaneous squamous cell carcinoma tumor-type indicator, 1 = CSCC",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Fractional change (1 + theta) on baseline CL (Appendix S1: IF (TUMTYPN.EQ.1) CL_TUMOR = 1 + THETA(32)). Reference category (all tumor indicators 0) is NSCLC; the tumor indicators are mutually exclusive.",
      source_name = "TUMTYPN"
    ),
    TUMTP_BCC = list(
      description = "Basal cell carcinoma tumor-type indicator, 1 = BCC",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Fractional change (1 + theta) on baseline CL and on the maximal time-dependent change in CL Emax (Appendix S1: TUMTYPN.EQ.2 -> CL_TUMOR = 1 + THETA(33), EMAX_TUMOR = 1 + THETA(36)). Reference category is NSCLC.",
      source_name = "TUMTYPN"
    ),
    TUMTP_OTHER = list(
      description = "Other tumor types indicator, 1 = tumor type other than CSCC, BCC or NSCLC",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Fractional change (1 + theta) on T50, the time to 50% of the maximal change in CL (Appendix S1: IF (TUMTYPN.EQ.5) ET50_TUMOR = 1 + THETA(40), commented 'OTHERS'). The pool is the miscellaneous advanced solid tumors of the first-in-human Study 1423 that are not CSCC, BCC or NSCLC; the paper does not enumerate its composition. The CL and Emax effects of this group were fixed to 0 (same as the NSCLC reference).",
      source_name = "TUMTYPN"
    )
  )

  compartmentData <- list(
    central = list(analyte = "cemiplimab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "cemiplimab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1063L,
    n_studies = 4L,
    n_observations = 17312L,
    age_range = "27-96 years",
    age_median = "65 years",
    weight_range = "30.9-171.6 kg",
    weight_median = "73.2 kg",
    sex_female_pct = 29.4,
    race_ethnicity = c(
      White = 87.7,
      Black = 2.0,
      Asian = 4.4,
      `American Indian or Alaska Native` = 0.7,
      Other = 0.3,
      `Unknown, not reported or missing` = 5.0
    ),
    disease_state = "Advanced solid tumors: cutaneous squamous cell carcinoma (Study 1540, n = 188), basal cell carcinoma after hedgehog-inhibitor therapy (Study 1620, n = 132), PD-L1 >= 50% non-small cell lung cancer (Study 1624, n = 345) and advanced malignancies in the first-in-human study (Study 1423, n = 398).",
    dose_range = "1, 3 or 10 mg/kg Q2W; 3 mg/kg Q3W; 200 mg Q2W; 350 mg Q3W; all 30-min IV infusions",
    regions = "Multinational (not enumerated)",
    notes = "Model-building data pooled from NCT02383212, NCT02760498, NCT03132636 and NCT03088540 (Table S1, Table S2). Median baseline albumin 38.6 g/L (range 20-93), ALT 19 IU/L (4-234); ECOG 0/1 = 38.1/61.9% (Table S4). External validation (not used in estimation) in 292 women with recurrent or metastatic cervical cancer from EMPOWER-Cervical 1 (NCT03257267, Study 1676) on 350 mg Q3W: median weight 61.9 kg, albumin 39.1 g/L, ALT 14 IU/L, 29.1% Asian (Table S5)."
  )

  ini({
    # Typical values for a male NSCLC patient at the reference covariates
    # WT 75 kg, albumin 38.6 g/L, ALT 19 IU/L (Appendix S1 control stream $PK).
    lcl <- log(0.254); label("Baseline clearance at first dose TVCL0 (L/day)") # Table 1: TVCL0 = 0.254 L/day
    lq <- log(0.652); label("Intercompartmental clearance TVQ (L/day)") # Table 1: TVQ = 0.652 L/day
    lvc <- log(3.35); label("Central volume of distribution TVV1 (L)") # Table 1: TVV1 = 3.35 L
    lvp <- log(2.52); label("Peripheral volume of distribution TVV2 (L)") # Table 1: TVV2 = 2.52 L

    # Time-dependent CL: CL(t) = CL0 * exp(Emax * t^HILL / (T50^HILL + t^HILL))
    # (Table 1 footnote; Appendix S1). Emax is a log-scale change, so CL
    # approaches exp(-0.174) = 0.840 of baseline.
    # Emax is negative (CL falls over time), so it is kept on the linear scale.
    cl_time_max <- -0.174; label("Maximal log-scale change in CL with time TVEmax (unitless)") # Table 1: TVEmax = -0.174
    lcl_t50 <- log(73.7); label("Time to the half-maximal change in CL TVT50 (day)") # Table 1: TVT50 = 73.7 days
    lcl_time_hill <- fixed(log(2.50)); label("Hill exponent of the time-dependent change in CL (unitless)") # Table 1: HILL = 2.50, held constant per footnote c

    # Covariate effects (Table 1; equation forms from the Appendix S1 control stream)
    e_wt_cl_q <- 0.539; label("Power exponent of baseline weight on CL and Q (unitless)") # Table 1: CL_WGTBL = 0.539
    e_wt_vc_vp <- 0.499; label("Power exponent of baseline weight on V1 and V2 (unitless)") # Table 1: VssWGTBL = 0.499
    e_sexf_cl <- -0.137; label("Fractional change in CL for female patients (unitless)") # Table 1: CL_SEX = -0.137
    e_sexf_vc <- -0.0801; label("Fractional change in V1 for female patients (unitless)") # Table 1: V1_SEX = -0.0801
    e_alb_base_vc <- -0.217; label("Power exponent of baseline albumin on V1 (unitless)") # Table 1: V1_ALBBL = -0.217
    e_alb_cl <- -1.11; label("Power exponent of time-varying albumin on CL (unitless)") # Table 1: CL_ALB = -1.11
    e_alt_cl <- -0.0729; label("Power exponent of baseline ALT on CL (unitless)") # Table 1: CL_ALTBL = -0.0729
    e_tumtp_cscc_cl <- -0.216; label("Fractional change in CL for CSCC vs NSCLC (unitless)") # Table 1: CL_CSCC = -0.216
    e_tumtp_bcc_cl <- -0.211; label("Fractional change in CL for BCC vs NSCLC (unitless)") # Table 1: CL_BCC = -0.211
    e_tumtp_bcc_cl_time_max <- -0.872; label("Fractional change in Emax of time-dependent CL for BCC (unitless)") # Table 1: Emax_BCC = -0.872
    e_tumtp_other_cl_t50 <- -0.491; label("Fractional change in T50 of time-dependent CL for other tumor types (unitless)") # Table 1: T50_OTHER = -0.491

    # IIV: one eta shared by CL and Q, one shared by V1 and V2, no covariance
    # (Base model section; Appendix S1 $PK CL/Q use ETA(1), V1/V2 use ETA(2)).
    etalcl ~ 0.0892 # Table 1: ETA1 - CL_Q variance = 0.0892
    etalvc ~ 0.0709 # Table 1: ETA2 - V1_V2 variance = 0.0709

    # Residual error: additive on log-transformed concentrations
    # (Appendix S1 $ERROR: Y = LOG(CONC) + W*ERR(1), SIGMA 1 fixed, W = THETA(6)).
    expSd <- 0.241; label("Additive residual error on the log scale (SD)") # Table 1: RE = 0.241
  })
  model({
    # Covariate effects; reference values from the Appendix S1 $PK block.
    cov_cl <- (WT / 75)^e_wt_cl_q *
      (ALB / 38.6)^e_alb_cl *
      (ALT / 19)^e_alt_cl *
      (1 + e_sexf_cl * SEXF) *
      (1 + e_tumtp_cscc_cl * TUMTP_CSCC + e_tumtp_bcc_cl * TUMTP_BCC)
    cov_vc <- (WT / 75)^e_wt_vc_vp *
      (ALB_BASE / 38.6)^e_alb_base_vc *
      (1 + e_sexf_vc * SEXF)

    cl_base <- exp(lcl + etalcl) * cov_cl
    q <- exp(lq + etalcl) * (WT / 75)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * cov_vc
    vp <- exp(lvp + etalvc) * (WT / 75)^e_wt_vc_vp

    # Sigmoid-Emax change in CL over time since the first dose.
    cl_time_hill <- exp(lcl_time_hill)
    cl_time_max_i <- cl_time_max * (1 + e_tumtp_bcc_cl_time_max * TUMTP_BCC)
    cl_t50_i <- exp(lcl_t50) * (1 + e_tumtp_other_cl_t50 * TUMTP_OTHER)
    cl <- cl_base * exp(cl_time_max_i * t^cl_time_hill / (cl_t50_i^cl_time_hill + t^cl_time_hill))

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> mg/L.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
