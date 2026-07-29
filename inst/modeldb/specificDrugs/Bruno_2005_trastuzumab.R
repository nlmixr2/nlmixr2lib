Bruno_2005_trastuzumab <- function() {
  description <- "Two-compartment linear population PK model for intravenous trastuzumab in adults with HER2-positive metastatic breast cancer (MBC) or advanced solid tumors; covariate effects of number of metastatic sites (>= 4) and baseline HER2 shed extracellular domain (ECD) on clearance, and body weight and ECD on central volume (Bruno 2005, first published trastuzumab popPK)."
  reference <- "Bruno R, Washington CB, Lu JF, Lieberman G, Banken L, Klein P. Population pharmacokinetics of trastuzumab in patients with HER2+ metastatic breast cancer. Cancer Chemother Pharmacol. 2005;56(4):361-369. doi:10.1007/s00280-005-1026-z"
  vignette <- "Bruno_2005_trastuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on central volume V (exponent 0.556); reference 70 kg (median of the analysis population per Bruno 2005 Table 2). No retained WT effect on linear CL in the final model.",
      source_name        = "WT"
    ),
    HER2_ECD = list(
      description        = "Baseline serum concentration of HER2 shed extracellular domain",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effects on linear CL (exponent 0.041) and on V (exponent 0.105). Reference 8.23 ng/mL (median baseline ECD in the analysis population, per Bruno 2005 Results 'Covariate effects on pharmacokinetic parameters'). ECD is capped at 200 ng/mL inside the model() block before applying the power law, reflecting Bruno 2005's observation that 'the relatively few ECD levels above 200 ng/ml were not associated with further increases in clearance' (Results and Figure 2 panel A). Values below the assay LLOQ (5.5 or 3.4 ng/mL depending on which of the two ELISA assays was used) were set to half the LLOQ in the analysis dataset per Bruno 2005 Methods.",
      source_name        = "ECD"
    ),
    MET_GE4 = list(
      description        = "Indicator of baseline number of metastatic sites >= 4",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fewer than four metastatic sites at baseline)",
      notes              = "Time-fixed baseline indicator. Multiplicative additive effect on linear CL: typical CL is scaled by (1 + 0.221 * MET_GE4), i.e. +22.1% in MET_GE4 = 1 patients. Reference MET_GE4 = 0 per Bruno 2005 Table 3. Bruno 2005 Methods defines MET = 1 when the number of metastatic sites is 4 or greater and 0 otherwise; the canonical column MET_GE4 inherits the same definition. Bruno 2005 Table 2 reports 53/476 (11%) patients with MET_GE4 = 1.",
      source_name        = "MET"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 476L,
    n_studies        = 4L,
    n_observations   = 3249L,
    age_range        = "25-81 years",
    age_median       = "50 years",
    weight_range     = "42-119 kg",
    weight_median    = "70 kg",
    sex_female_pct   = NA_real_,
    disease_state    = "HER2-positive metastatic breast cancer (n = 460 across phase II and III studies) plus 16 patients with advanced HER2-overexpressing or non-overexpressing solid tumors from a phase I single-dose study. HER2 status determined by immunohistochemistry (IHC); the analysis population was predominantly IHC 3+ (76%) and 2+ (22%).",
    dose_range       = "Trastuzumab 4 mg/kg IV loading dose followed by 2 mg/kg weekly (phase II and pivotal phase II/III studies, 30-90 minute infusion) for up to 840 days. Phase I study: single-dose escalation 10-500 mg total dose (n = 16 evaluable). Concomitant chemotherapy in 211/476 patients: anthracycline plus cyclophosphamide (n = 133) or paclitaxel (n = 78). Concomitant chemotherapy did not influence CL or V (Bruno 2005 Results) and is not a retained model covariate.",
    regions          = "Multinational (the four pooled trials enrolled subjects under FDA-registration programs leading to trastuzumab approval).",
    reference_subject = "70 kg, baseline HER2_ECD 8.23 ng/mL, MET_GE4 = 0 (the typical patient Bruno 2005 uses for the final CL = 0.225 L/day and V = 2.95 L estimates).",
    notes            = "Patient characteristics from Bruno 2005 Table 2. The pooled dataset combined intensive single-dose profiles (phase I, n = 16) with peak-and-trough sampling from phase II (n = 46), pivotal phase II (n = 213, 4 mg/kg + 2 mg/kg weekly, 2nd / 3rd line MBC), and pivotal phase III (n = 235 randomized to trastuzumab plus chemotherapy in the H0648g trial, 1st-line MBC). Trastuzumab serum concentrations measured by ELISA (LLOQ 0.156 ug/mL); BLQ values excluded. ECD baseline values measured by ELISA (LLOQ 5.5 or 3.4 ng/mL depending on assay; sub-LLOQ values imputed as LLOQ/2). NONMEM V FOCE-NONE (first-order) was used; bootstrap (500 replicates) confirmed parameter stability. Sex distribution not tabulated in Bruno 2005 Table 2 but implied predominantly female (MBC cohort)."
  )

  ini({
    # Structural parameters (Bruno 2005 Table 3 final population estimates).
    # Reference subject: 70 kg adult, HER2_ECD 8.23 ng/mL, fewer than four
    # metastatic sites. Bruno 2005 parameterised the two-compartment model
    # with the microconstants K12 and K21 rather than the canonical Q and Vp;
    # the canonical values stored here are the deterministic identity
    #   Q  = K12 * V1
    #   Vp = K12 * V1 / K21
    # Bruno 2005 Table 3: K12 = 0.164 /day, K21 = 0.101 /day, V1 = 2.95 L,
    # giving Q = 0.4838 L/day and Vp = 4.7901 L. Concentration in the central
    # compartment is Cc = central / vc, with dose in mg and volumes in L ->
    # Cc in mg/L (= ug/mL).
    lcl <- log(0.225);   label("Linear CL for the reference subject (L/day)")               # Bruno 2005 Table 3, CL point estimate
    lvc <- log(2.95);    label("Central volume of distribution V1 for the reference subject (L)") # Bruno 2005 Table 3, V point estimate
    lq  <- log(0.4838);  label("Intercompartmental clearance Q (L/day) = K12 * V1")           # Bruno 2005 Table 3, K12 = 0.164 /day * V1 = 2.95 L
    lvp <- log(4.7901);  label("Peripheral volume of distribution Vp (L) = K12 * V1 / K21")   # Bruno 2005 Table 3, K12 * V1 / K21

    # Covariate effects (Bruno 2005 Table 3 and the final-model covariate
    # equation in Results). Bruno 2005 implements continuous covariates as
    # power functions on the median-centered ratio
    #   theta_typ = theta_pop * (cov / median)^h_cov
    # and dichotomous covariates as the multiplicative additive form
    #   theta_typ = theta_pop * (1 + h_cov * MET)
    # MET = 1 iff number of metastatic sites >= 4, i.e. the MET_GE4 indicator.
    e_met_ge4_cl <- 0.221; label("Multiplicative-additive coefficient of MET_GE4 on linear CL (unitless)") # Bruno 2005 Table 3, MET on CL
    e_her2_ecd_cl <- 0.041; label("Power exponent of HER2_ECD on linear CL (unitless; reference 8.23 ng/mL, ECD capped at 200 ng/mL)") # Bruno 2005 Table 3, ECD on CL
    e_wt_vc <- 0.556; label("Power exponent of WT on central volume V (unitless; reference 70 kg)")          # Bruno 2005 Table 3, WT on V
    e_her2_ecd_vc <- 0.105; label("Power exponent of HER2_ECD on central volume V (unitless; reference 8.23 ng/mL, ECD capped at 200 ng/mL)") # Bruno 2005 Table 3, ECD on V

    # Inter-individual variability (Bruno 2005 Table 3). Bruno 2005 reports
    # omega as %CV on log-normal parameters; convert via omega^2 = log(CV^2 + 1).
    # IIV is reported on CL (43% CV) and V1 (29% CV); Bruno 2005 also reports
    # IIV on the microconstants K12 (54% CV) and K21 (67% CV). The canonical
    # reparameterisation in this file expresses the typical-value structure
    # with Q and Vp; the K12/K21 random effects do not translate cleanly to
    # independent IIV on Q and Vp (Q and Vp both inherit variance from K12
    # and V1, with Vp inheriting additional variance from K21), and Bruno
    # 2005 does not report the joint covariance structure needed to
    # propagate them faithfully. The K12/K21 random effects are therefore
    # omitted here and documented in the vignette Assumptions and
    # deviations section; the retained CL and V1 IIVs dominate steady-state
    # exposure variability and are sufficient for the intended forward-
    # simulation use cases.
    etalcl ~ 0.169823   # CL 43% CV -- Bruno 2005 Table 3; omega^2 = log(0.43^2 + 1)
    etalvc ~ 0.080728   # V  29% CV -- Bruno 2005 Table 3; omega^2 = log(0.29^2 + 1)

    # Residual error: proportional only (Bruno 2005 Methods: "Residual error
    # was modeled as proportional, i.e., assuming a constant coefficient of
    # variation for error over the range of measured concentrations").
    # Table 3 reports sigma = 23% as a CV; stored as a fraction in nlmixr2's
    # linear-space proportional parameterisation.
    propSd <- 0.23; label("Proportional residual error (fraction)")  # Bruno 2005 Table 3, sigma 23%
  })

  model({
    # Cap HER2_ECD at 200 ng/mL before applying the power-law effects to
    # reproduce Bruno 2005's observation that "the relatively few ECD levels
    # above 200 ng/ml were not associated with further increases in
    # clearance" (Results 'Covariate effects on pharmacokinetic parameters'
    # and Figure 2 panel A). For ECD <= 200 ng/mL the cap is inactive.
    her2_ecd_eff <- min(HER2_ECD, 200)

    # Individual linear CL (Bruno 2005 final CL equation):
    #   CL_i = CL_pop * (HER2_ECD / 8.23)^e_her2_ecd_cl
    #               * (1 + e_met_ge4_cl * MET_GE4) * exp(eta_CL)
    cl <- exp(lcl + etalcl) *
      (her2_ecd_eff / 8.23)^e_her2_ecd_cl *
      (1 + e_met_ge4_cl * MET_GE4)

    # Individual central volume (Bruno 2005 final V equation):
    #   V_i = V_pop * (WT / 70)^e_wt_vc * (HER2_ECD / 8.23)^e_her2_ecd_vc
    #              * exp(eta_V)
    vc <- exp(lvc + etalvc) *
      (WT / 70)^e_wt_vc *
      (her2_ecd_eff / 8.23)^e_her2_ecd_vc

    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    Cc <- central / vc

    # Two-compartment IV model with linear elimination from the central
    # compartment (Bruno 2005 Results: "A two-compartment linear
    # pharmacokinetic model best described the data and accounted for the
    # long-term accumulation observed following weekly administration of
    # trastuzumab"). All trastuzumab doses enter the central compartment
    # directly via 30-90 minute IV infusion.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    Cc ~ prop(propSd)
  })
}
