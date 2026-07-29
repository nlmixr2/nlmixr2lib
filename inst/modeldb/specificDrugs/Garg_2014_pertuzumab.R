Garg_2014_pertuzumab <- function() {
  description <- "Two-compartment population PK model with first-order linear elimination from the central compartment for intravenous pertuzumab (PERJETA) in patients with a variety of HER2-targeted solid tumors (Garg 2014)"
  reference <- "Garg A, Quartino A, Li J, Jin J, Wada DR, Li H, Cortes J, McNally V, Ross G, Visich J, Lum B. Population pharmacokinetic and covariate analysis of pertuzumab, a HER2-targeted monoclonal antibody, and evaluation of a fixed, non-weight-based dose in patients with a variety of solid tumors. Cancer Chemother Pharmacol. 2014;74(4):819-829. doi:10.1007/s00280-014-2560-3"
  vignette <- "Garg_2014_pertuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    LBM = list(
      description        = "Lean body weight (canonical column LBM; source paper uses LBW)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on linear CL (exponent 0.516), Vc (exponent 0.747), and Vp (exponent 0.83), each centred at the cohort median LBW = 48 kg (Garg 2014 Table 1 footnote and CL/Vc/Vp covariate equations on page 823). Body-composition formula not stated in the paper; mAb popPK literature most commonly uses the Hume (1966) or James (1976) formula for LBW. Stored under canonical LBM (lean body mass) per inst/references/covariate-columns.md; LBW and LBM refer to the same quantity. Body weight, BSA, sex, age, race (Japanese vs non-Japanese), and other patient demographics were screened during covariate analysis but only LBW was retained on volumes.",
      source_name        = "LBW"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value (US convention, g/dL). Power effect on linear CL (exponent -1.06, i.e. approximately inversely proportional); reference 3.9 g/dL per Garg 2014 CL covariate equation (page 823, typical patient). Source column 'ALBU' maps to the canonical ALB covariate.",
      source_name        = "ALBU"
    )
  )

  population <- list(
    n_subjects       = 481L,
    n_studies        = 12L,
    n_observations   = 4525L,
    phase_mix        = "Pooled five phase I/Ib studies, six phase II studies, and the pivotal phase III trial CLEOPATRA. 392 patients received pertuzumab as a single agent; 356 received pertuzumab in combination with chemotherapy (gemcitabine, capecitabine, or docetaxel) and/or targeted therapy (trastuzumab or erlotinib). 20 patients were enrolled in CLEOPATRA.",
    age_median       = "60 years",
    weight_median    = "72 kg",
    sex_female_pct   = 62.0,
    race_ethnicity   = c(`non-Japanese` = 95.4, Japanese = 4.6),
    disease_state    = "Solid tumors including metastatic breast cancer (MBC), non-small cell lung cancer (NSCLC), ovarian cancer, prostate cancer, and other HER2-targeted indications. 81.5% (392/481) received pertuzumab without chemotherapy; the remainder received pertuzumab with chemotherapy and/or other targeted therapy.",
    dose_range       = "Body-weight-based: 0.5-25 mg/kg IV q3w (n = 39 in phase Ia dose-ranging; the three patients receiving 0.5 mg/kg were excluded from the population PK analysis owing to nonlinear PK at low concentrations and analyses were performed on doses >=2 mg/kg). Fixed dose: 1050 mg IV q3w (n = 96, 20%) or 840 mg IV loading dose followed by 420 mg IV q3w maintenance (n = 346, 72%; the labelled clinical regimen).",
    regions          = "Global (multinational phase I-III studies).",
    ecog_status      = "Almost all patients (479/481, 99.6%) had ECOG performance status 0 or 1.",
    reference_subject = "Lean body weight 48 kg, serum albumin 3.9 g/dL (the cohort median per Garg 2014 Table 1 footnote and CL/Vc/Vp covariate equations on page 823).",
    notes            = "Baseline demographics per Garg 2014 Results section 'Patient population' and Online Resource 2. The PK analysis dataset comprised 4,525 pertuzumab serum concentration time points; 115 (2.5%) data points from 65 patients were flagged as outliers (trough above peak, pharmacologically unexplained spike/drop, or below LLOQ) and omitted before fitting. Phase Ia patients receiving 0.5 mg/kg (n = 3) were excluded after preliminary nonlinear PK exploration was inconclusive. CRP, HER2 ECD, and HER2 expression were assessed by sensitivity analysis on subsets where the data were available, not as formal covariates."
  )

  ini({
    # Structural parameters (Garg 2014 Table 1, final FOCEI estimates).
    # Reference subject: LBW 48 kg, ALBU 3.9 g/dL (cohort median per Garg 2014
    # Table 1 footnote). Concentration in the central compartment is
    # Cc = central / vc, with dose in mg and volumes in L -> Cc in mg/L
    # (= ug/mL).
    lcl <- log(0.235); label("Linear CL for the reference subject (L/day)") # Garg 2014 Table 1, theta1
    lvc <- log(3.11);  label("Central volume of distribution for the reference subject (L)") # Garg 2014 Table 1, theta2
    lq  <- log(0.534); label("Intercompartmental clearance Q (L/day)") # Garg 2014 Table 1, theta3 (no covariate effect retained)
    lvp <- log(2.46);  label("Peripheral volume of distribution for the reference subject (L)") # Garg 2014 Table 1, theta4

    # Covariate effects (Garg 2014 Table 1 and CL / Vc / Vp covariate
    # equations on page 823). All four covariate effects are power-form,
    # centred at the cohort median LBW = 48 kg and ALBU = 3.9 g/dL.
    # CL_i = theta1 * (LBW / 48)^theta5 * (ALBU / 3.9)^theta7 * exp(eta_CL)
    # Vc_i = theta2 * (LBW / 48)^theta6                      * exp(eta_Vc)
    # Vp_i = theta4 * (LBW / 48)^theta8                      * exp(eta_Vp)
    e_lbw_cl <-  0.516; label("Power exponent of LBW on linear CL (unitless)") # Garg 2014 Table 1, theta5
    e_lbw_vc <-  0.747; label("Power exponent of LBW on Vc (unitless)") # Garg 2014 Table 1, theta6
    e_alb_cl <- -1.06;  label("Power exponent of ALB (ALBU) on linear CL (unitless)") # Garg 2014 Table 1, theta7
    e_lbw_vp <-  0.83;  label("Power exponent of LBW on Vp (unitless)") # Garg 2014 Table 1, theta8

    # Inter-individual variability. Garg 2014 Table 1 reports omega as %CV on
    # log-normal parameters; convert via omega^2 = log(CV^2 + 1):
    #   CL  34.1% CV -> 0.1100   (= log(1 + 0.341^2))
    #   Vc  18.5% CV -> 0.0337   (= log(1 + 0.185^2))
    #   Vp  45.9% CV -> 0.1912   (= log(1 + 0.459^2))
    # Garg 2014 Table 1 footnote also reports the off-diagonal covariances
    # directly (NONMEM OMEGA scale): omega_CL,Vc = 0.0239, omega_CL,Vp =
    # -0.0416, omega_Vc,Vp = 0.000179. Lower-triangle order is var(CL),
    # cov(CL,Vc), var(Vc), cov(CL,Vp), cov(Vc,Vp), var(Vp).
    etalcl + etalvc + etalvp ~ c(0.1100,
                                  0.0239,  0.0337,
                                 -0.0416,  0.000179, 0.1912) # Garg 2014 Table 1 + footnote

    # Residual error: Garg 2014 Table 1 reports an additive residual on
    # log-transformed concentrations (NONMEM "additive on log scale"), which
    # maps to nlmixr2's proportional residual on the linear scale. The
    # tabulated value is 18.1%.
    propSd <- 0.181; label("Proportional residual error (fraction)") # Garg 2014 Table 1, proportional residual error
  })

  model({
    # SI -> US-convention unit conversion (canonical ALB is in SI g/L per the
    # 2026-06-19 register standardization audit; the original calibration
    # used the g/dL reference value, so convert inline here).
    alb_gdL <- ALB * 0.1  # SI g/L -> US-convention g/dL (factor 0.1)

    # Individual PK parameters with LBW power scaling (reference 48 kg) on
    # CL / Vc / Vp and alb_gdL power scaling (reference 3.9 g/dL) on CL.
    cl <- exp(lcl + etalcl) * (LBM / 48)^e_lbw_cl * (alb_gdL / 3.9)^e_alb_cl
    vc <- exp(lvc + etalvc) * (LBM / 48)^e_lbw_vc
    q  <- exp(lq)
    vp <- exp(lvp + etalvp) * (LBM / 48)^e_lbw_vp

    # Two-compartment IV model with first-order linear elimination from the
    # central compartment (Garg 2014 Results: "linear two-compartment model
    # with first-order elimination from the central compartment"). Dose is
    # delivered directly into central; no depot or absorption process.
    Cc <- linCmt()
    Cc ~ prop(propSd)
  })
}
