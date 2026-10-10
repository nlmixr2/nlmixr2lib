vanEijk_2022_paclitaxel <- function() {
  description <- "Semi-physiological population PK/PD model for oral paclitaxel (drinking solution, ModraPac capsule and ModraPac tablet) co-administered with oral ritonavir in adult cancer patients: Weibull-type gut absorption into a 1 L well-stirred liver compartment whose intrinsic clearance is inhibited by the ritonavir plasma concentration (Imax model), two-compartment systemic disposition, an embedded two-compartment ritonavir PK model with inverse Gaussian absorption (Yu 2020), and a thrombospondin-1 (TSP-1) turnover model whose formation rate is stimulated by paclitaxel (Emax = 1)"
  reference <- paste(
    "van Eijk M, Yu H, Sawicki E, de Weger VA, Nuijen B, Dorlo TPC, Beijnen JH,",
    "Huitema ADR. Development of a population pharmacokinetic/pharmacodynamic",
    "model for various oral paclitaxel formulations co-administered with ritonavir",
    "and thrombospondin-1 based on data from early phase clinical studies.",
    "Cancer Chemother Pharmacol. 2022;90(1):71-82. doi:10.1007/s00280-022-04445-z.",
    "Embedded ritonavir PK model (structure and all ritonavir parameter values)",
    "from Yu H, Janssen JM, Sawicki E, van Hasselt JGC, de Weger VA, Nuijen B,",
    "Schellens JHM, Beijnen JH, Huitema ADR. A Population Pharmacokinetic Model",
    "of Oral Docetaxel Coadministered With Ritonavir to Support Early Clinical",
    "Development. J Clin Pharmacol. 2020;60(3):340-350. doi:10.1002/jcph.1532",
    "(van Eijk 2022 reference [29]; Table 2 and Equation 3).",
    sep = " "
  )
  vignette <- "vanEijk_2022_paclitaxel"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # TSP-1 is a turnover biomarker state of this paper only (thrombospondin-1
  # quantified relative to platelet count); it is not a reusable compartment.
  paper_specific_compartments <- c("tsp1")

  covariateData <- list(
    FORM_TABLET = list(
      description = "Paclitaxel formulation indicator: 1 = ModraPac tablet (spray-dried amorphous solid dispersion), 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 together with FORM_CAPSULE = 0 (paclitaxel drinking solution, the oral liquid reference with rF = 1).",
      notes = paste(
        "Three-level paclitaxel formulation carried as two indicators (FORM_TABLET,",
        "FORM_CAPSULE); the reference oral liquid is the drinking solution (registered",
        "IV paclitaxel formulation, 6 mg/mL in ethanol / Cremophor EL, given orally).",
        "Tablet: relative gut bioavailability rF = 0.97 and Weibull shape BETA = 3.57",
        "(shared with the capsule) versus BETA = 2.53 for the drinking solution",
        "(Table 1). FORM_TABLET and FORM_CAPSULE must not both be 1. Set on the",
        "paclitaxel dose records (and carry forward), because the Weibull shape is",
        "evaluated continuously during absorption.",
        sep = " "
      ),
      source_name = "formulation (drinking solution / ModraPac capsule / ModraPac tablet)"
    ),
    FORM_CAPSULE = list(
      description = "Paclitaxel formulation indicator: 1 = ModraPac capsule (freeze-dried amorphous solid dispersion), 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 together with FORM_TABLET = 0 (paclitaxel drinking solution, rF = 1).",
      notes = paste(
        "Capsule: relative gut bioavailability rF = 0.46 and Weibull shape BETA =",
        "3.57 (shared with the tablet) versus the drinking-solution reference",
        "(Table 1). The comparator (FORM_CAPSULE = 0, FORM_TABLET = 0) is the",
        "paclitaxel drinking solution, whose rF is fixed to 1.",
        sep = " "
      ),
      source_name = "formulation (drinking solution / ModraPac capsule / ModraPac tablet)"
    ),
    FORM_RTV_TABLET = list(
      description = "Ritonavir formulation indicator: 1 = Norvir tablet, 0 = Norvir soft-gel capsule",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ritonavir capsule, the Yu 2020 reference with relative F = 1).",
      notes = paste(
        "Covariate of the embedded Yu 2020 ritonavir model: relative ritonavir",
        "bioavailability of tablet versus capsule Ftablet/capsule = 1.06 with 30%",
        "CV between-subject variability (Yu 2020 Table 2). van Eijk 2022 does not",
        "report which ritonavir formulation its three studies used (all brand",
        "Norvir), so the value must be chosen by the user; the reference capsule",
        "(0) is the default used in the validation vignette.",
        sep = " "
      ),
      source_name = "Ftablet/capsule (Yu 2020)"
    ),
    OCC = list(
      description = "Dosing occasion index (1 or 2) for the between-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "van Eijk 2022 Methods: 'Each dose administration was considered an",
        "occasion'. In the analysis data a subject had at most two PK occasions",
        "(Study 2: two weekly doses; Study 3: the two doses of day 1 of the",
        "twice-daily schedule), so two occasion slots are encoded. OCC = 1 or 2",
        "selects the occasion's BOV eta on paclitaxel relative gut bioavailability",
        "(read at the dose record) and the occasion's within-subject eta on the",
        "ritonavir mean absorption time and absorption-time dispersion (Yu 2020",
        "Table 2). Any other value (e.g. 0) switches the occasion effects off.",
        "For long simulations draw a new eta per dose instead (see the vignette).",
        sep = " "
      ),
      source_name = "occasion (each dose administration)"
    )
  )

  compartmentData <- list(
    depot1 = list(analyte = "paclitaxel", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "paclitaxel", units = "mg", specimen = "administration site", verified = TRUE),
    liver = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    central = list(analyte = "paclitaxel", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "paclitaxel", units = "mg", specimen = "plasma", verified = TRUE),
    depot1_rtv = list(analyte = "ritonavir", units = "mg", specimen = "administration site", verified = TRUE),
    depot2_rtv = list(analyte = "ritonavir", units = "mg", specimen = "administration site", verified = TRUE),
    central_rtv = list(analyte = "ritonavir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_rtv = list(analyte = "ritonavir", units = "mg", specimen = "plasma", verified = TRUE),
    tsp1 = list(analyte = "thrombospondin-1", units = "ng/mL/10^6 platelets", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 58,
    n_subjects_pd = 36,
    n_studies = 3,
    disease_state = "adult cancer patients in three early-phase clinical studies of oral paclitaxel boosted with ritonavir",
    dose_range = paste(
      "Study 1: single 100 mg paclitaxel drinking solution with 100 or 200 mg ritonavir",
      "30 min before (17 patients); Study 2: 30 mg once weekly as drinking solution and",
      "ModraPac capsule in randomised order with 100 mg ritonavir 30 min before",
      "(4 patients); Study 3 (low-dose metronomic, NTR3632): ModraPac capsule 5-40 mg/day",
      "or tablet 40-60 mg/day given twice daily 7 h apart with ritonavir 200 mg/day (37 patients)"
    ),
    regions = "The Netherlands (Netherlands Cancer Institute)",
    notes = paste(
      "PK data from 58 patients in three studies (van Eijk 2022 Methods,",
      "'Pharmacokinetic data'); TSP-1 PD data from 36 of the 37 Study 3 patients",
      "('Pharmacodynamic data'). The paper reports no demographic table (age, weight,",
      "sex) and its Supplementary Table 1 (dose and sampling summary) is not part",
      "of the published supplementary file. PK sampling in Study 3 covered day 1",
      "only (Discussion).",
      sep = " "
    )
  )

  ini({
    # ---- Paclitaxel Weibull absorption (van Eijk 2022 Eq. 9, Table 1) -------
    # The paper parameterises the Weibull with a time scale ALPHA (h) and a
    # shape BETA. The registered rate-scaling canonical `ra` is the reciprocal
    # of the scale, so ra = 1 / ALPHA.
    lra_dose1 <- log(1 / 1.68); label("Weibull rate scale of the first daily paclitaxel dose, 1/ALPHA (1/h)") # Table 1 'ALPHA 1st daily dose' = 1.68 h (95% CI 1.52-1.87)
    lra_dose2 <- log(1 / 1.97); label("Weibull rate scale of the second daily paclitaxel dose, 1/ALPHA (1/h)") # Table 1 'ALPHA 2nd daily dose' = 1.97 h (95% CI 1.79-2.19)
    lgam1_sol <- log(2.53); label("Weibull shape BETA for the paclitaxel drinking solution (unitless)") # Table 1 'BETA drinking solution' = 2.53 (95% CI 2.34-2.76)
    lgam1_asd <- log(3.57); label("Weibull shape BETA for the ModraPac capsule and tablet (unitless)") # Table 1 'BETA tablet+capsule' = 3.57 (95% CI 2.99-4.52)

    # ---- Relative gut bioavailability (Table 1) ------------------------------
    lfdepot <- fixed(log(1)); label("Relative gut bioavailability rF of the paclitaxel drinking solution (fraction)") # Table 1 'rF drinking solution' = 1 FIX
    e_form_tablet_fdepot <- 0.97; label("Relative gut bioavailability of the ModraPac tablet vs the drinking solution (fraction)") # Table 1 'rF tablet' = 0.97 (95% CI 0.67-1.33)
    e_form_capsule_fdepot <- 0.46; label("Relative gut bioavailability of the ModraPac capsule vs the drinking solution (fraction)") # Table 1 'rF capsule' = 0.46 (95% CI 0.34-0.61)
    e_dose2_fdepot <- 0.59; label("Relative gut bioavailability of the second daily paclitaxel dose vs the first (fraction)") # Table 1 'rF 2nd/1st' = 0.59 (95% CI 0.48-0.74)

    # ---- Well-stirred liver and ritonavir inhibition (Eqs. 1-3, Table 1) -----
    lclint <- log(746); label("Uninhibited intrinsic clearance CLint0 of paclitaxel (L/h)") # Table 1 'CL int0' = 746 L/h (95% CI 585-937)
    limax <- log(570); label("Maximum ritonavir-induced reduction of paclitaxel intrinsic clearance Imax (L/h)") # Table 1 'I max' = 570 L/h (95% CI 400-776)
    lki <- log(375); label("Ritonavir plasma concentration giving half of Imax, KI (ng/mL)") # Table 1 'KI' = 375 ng/mL (95% CI 135-906)
    qh <- fixed(80); label("Hepatic blood flow QH (L/h)") # Methods 'Pharmacokinetic model': QH fixed at 80 L/h (ref [31])
    vh <- fixed(1); label("Liver compartment volume VH (L)") # Methods 'Pharmacokinetic model': VH fixed to 1 L (ref [32])
    fu <- fixed(0.13); label("Fraction of paclitaxel unbound in plasma, literature value (fraction)") # Methods 'Pharmacokinetic model': literature fu of 13% (ref [33])

    # ---- Paclitaxel systemic disposition (Table 1) ---------------------------
    lvc <- log(128); label("Central volume of distribution Vc (L)") # Table 1 'Vc' = 128 L (95% CI 105-151)
    lq <- log(33.4); label("Intercompartmental clearance Q (L/h)") # Table 1 'Q' = 33.4 L/h (95% CI 29.6-37.6)
    lvp <- log(375); label("Peripheral volume of distribution Vp (L)") # Table 1 'Vp' = 375 L (95% CI 311-465)

    # ---- TSP-1 turnover model (Eqs. 4-6, Table 2) ----------------------------
    lrbase <- log(43.8); label("Baseline TSP-1 concentration E_BASE (ng/mL/10^6 platelets)") # Table 2 'E BASE' = 43.8 (95% CI 39.7-48.5)
    lkout <- fixed(log(1 / 233)); label("TSP-1 first-order elimination rate kout = 1/Turnover, literature platelet survival (1/h)") # Table 2 'Turnover' = 233 h FIX (Results: platelet survival 9.7 days, ref [46]); Eq. 5 kout = 1/Turnover
    lec50 <- log(284); label("Paclitaxel plasma concentration giving half-maximal stimulation of TSP-1 formation EC50 (ng/mL)") # Table 2 'EC50' = 284 ng/mL (95% CI 122-724)

    # ---- Embedded ritonavir PK model (Yu 2020 Table 2 and Eq. 3) -------------
    # van Eijk 2022 applied the previously developed Yu 2020 ritonavir model to
    # its own data without re-estimating it ('Pharmacokinetic model' and
    # 'Oral paclitaxel pharmacokinetic model'), so every ritonavir value is
    # held at the published Yu 2020 estimate.
    lmat_rtv <- fixed(log(8.45)); label("Ritonavir mean absorption time MAT of the inverse Gaussian input, taken from Yu 2020 (h)") # Yu 2020 Table 2 'MAT' = 8.45 h (RSE 5%)
    lcvabs_rtv <- fixed(log(1.23)); label("Ritonavir relative dispersion CV of the inverse Gaussian absorption-time density, taken from Yu 2020 (unitless)") # Yu 2020 Table 2 'CV' = 123% (RSE 3%); Eq. 3 uses CV^2
    lcl_rtv <- fixed(log(7.72)); label("Ritonavir apparent clearance CL, taken from Yu 2020 (L/h)") # Yu 2020 Table 2 'CLRTV' = 7.72 L/h (RSE 9%)
    lvc_rtv <- fixed(log(23)); label("Ritonavir apparent central volume Vc, taken from Yu 2020 (L)") # Yu 2020 Table 2 'VcRTV' = 23 L (RSE 15%)
    lq_rtv <- fixed(log(3.99)); label("Ritonavir apparent intercompartmental clearance Q, taken from Yu 2020 (L/h)") # Yu 2020 Table 2 'QRTV' = 3.99 L/h (RSE 15%)
    lvp_rtv <- fixed(log(17.9)); label("Ritonavir apparent peripheral volume Vp, taken from Yu 2020 (L)") # Yu 2020 Table 2 'VpRTV' = 17.9 L (RSE 12%)
    lfdepot_rtv <- fixed(log(1)); label("Relative bioavailability of the first daily ritonavir capsule dose, taken from Yu 2020 (fraction)") # Yu 2020: F is the relative-bioavailability reference (typical value 1) carrying the 52.2% BSV of Table 2
    lfdepot2_rtv <- fixed(log(2.25)); label("Relative bioavailability of the second daily ritonavir dose vs the first, taken from Yu 2020 (fraction)") # Yu 2020 Table 2 'F2nd/1st,rtv' = 2.25 (RSE 7%)
    lftab_rtv <- fixed(log(1.06)); label("Relative bioavailability of the ritonavir tablet vs the capsule, taken from Yu 2020 (fraction)") # Yu 2020 Table 2 'Ftablet/capsule' = 1.06 (RSE 12%)

    # ---- Between-subject variability -----------------------------------------
    # CV% converted with omega^2 = log(1 + CV^2) (exponential model, Eq. 7).
    etalra ~ 0.116183 # Table 1 BSV 'ALPHA' 35.1 CV% (shared by both daily doses); shrinkage 7%
    etalclint ~ 0.061096 # Table 1 BSV 'CL int0' 25.1 CV%; shrinkage 24%
    etalvc ~ 0.254211 # Table 1 BSV 'Vc' 53.8 CV%; shrinkage 13%
    etalfdepot ~ 0.136211 # Table 1 BSV 'rF gut' 38.2 CV%; shrinkage 32%
    etalrbase ~ 0.076520 # Table 2 BSV 'E BASE' 28.2 CV%; shrinkage 4%
    etalcvabs_rtv ~ fixed(0.016251) # Yu 2020 Table 2 BSV 'CV' 12.8 CV%
    etalcl_rtv ~ fixed(0.197283) # Yu 2020 Table 2 BSV 'CLRTV' 46.7 CV%
    etalvc_rtv ~ fixed(0.628195) # Yu 2020 Table 2 BSV 'VcRTV' 93.5 CV%
    etalfdepot_rtv ~ fixed(0.240971) # Yu 2020 Table 2 BSV 'F' 52.2 CV%
    etalfdepot2_rtv ~ fixed(0.106363) # Yu 2020 Table 2 BSV 'F2nd/1st' 33.5 CV%
    etalftab_rtv ~ fixed(0.086178) # Yu 2020 Table 2 BSV 'Ftablet/capsule' 30 CV%

    # ---- Between-occasion / within-subject variability (OCC slots 1 and 2) ---
    etaiov_lfdepot_1 ~ 0.190425 # Table 1 BOV 'rF gut' 45.8 CV%, occasion 1
    etaiov_lfdepot_2 ~ fixed(0.190425) # same variance as occasion 1
    etaiov_lmat_rtv_1 ~ fixed(0.098071) # Yu 2020 Table 2 within-subject 'MAT' 32.1 CV%, occasion 1
    etaiov_lmat_rtv_2 ~ fixed(0.098071) # same variance as occasion 1
    etaiov_lcvabs_rtv_1 ~ fixed(0.048108) # Yu 2020 Table 2 within-subject 'CV' 22.2 CV%, occasion 1
    etaiov_lcvabs_rtv_2 ~ fixed(0.048108) # same variance as occasion 1

    # ---- Residual error (Eq. 8) ----------------------------------------------
    propSd <- 0.258; label("Proportional residual error, paclitaxel (fraction)") # Table 1 'sigma prop' = 25.8 CV%
    propSd_rtv <- fixed(0.352); label("Proportional residual error, ritonavir, taken from Yu 2020 (fraction)") # Yu 2020 Table 2 'Proportional residual error' = 35.2 CV%
    propSd_tsp1 <- 0.138; label("Proportional residual error, TSP-1 (fraction)") # Table 2 'sigma prop' = 13.8 CV%
  })

  model({
    # 1. Occasion indicators (Eq. 7: P_i = P * exp(eta_BSV + eta_BOV))
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_lfdepot_1 + oc2 * etaiov_lfdepot_2
    iov_mat_rtv <- oc1 * etaiov_lmat_rtv_1 + oc2 * etaiov_lmat_rtv_2
    iov_cvabs_rtv <- oc1 * etaiov_lcvabs_rtv_1 + oc2 * etaiov_lcvabs_rtv_2

    # 2. Individual paclitaxel parameters
    ra1 <- exp(lra_dose1 + etalra)
    ra2 <- exp(lra_dose2 + etalra)
    asd <- FORM_TABLET + FORM_CAPSULE
    gam1 <- exp(lgam1_sol) * (1 - asd) + exp(lgam1_asd) * asd
    clint0 <- exp(lclint + etalclint)
    imax <- exp(limax)
    ki <- exp(lki)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    fform <- (1 - asd) + e_form_tablet_fdepot * FORM_TABLET + e_form_capsule_fdepot * FORM_CAPSULE
    fgut <- exp(lfdepot + etalfdepot + iov_fdepot) * fform

    # 3. Individual ritonavir parameters (Yu 2020)
    mat_rtv <- exp(lmat_rtv + iov_mat_rtv)
    cvabs_rtv <- exp(lcvabs_rtv + etalcvabs_rtv + iov_cvabs_rtv)
    cl_rtv <- exp(lcl_rtv + etalcl_rtv)
    vc_rtv <- exp(lvc_rtv + etalvc_rtv)
    q_rtv <- exp(lq_rtv)
    vp_rtv <- exp(lvp_rtv)
    f_rtv <- exp(lfdepot_rtv + etalfdepot_rtv) *
      exp((lftab_rtv + etalftab_rtv) * FORM_RTV_TABLET)
    f2_rtv <- exp(lfdepot2_rtv + etalfdepot2_rtv)

    # 4. Paclitaxel Weibull absorption, Eq. 9 as printed:
    #    ka(t) = (BETA/ALPHA) * (t/ALPHA)^(BETA-1) * exp(-(t/ALPHA)^BETA),
    # applied as a time-varying first-order rate constant on the gut amount,
    # with t the time after the last dose into that gut compartment. The first
    # and second daily doses are dosed into depot1 and depot2, which carry the
    # dose-specific ALPHA and rF. tad0() is 0 before the first dose, and the
    # 1e-8 floor keeps t^(BETA-1) well defined at the dose instant (BETA > 1).
    t_d1 <- max(tad0(depot1), 1e-8)
    t_d2 <- max(tad0(depot2), 1e-8)
    ka1 <- gam1 * ra1 * (ra1 * t_d1)^(gam1 - 1) * exp(-(ra1 * t_d1)^gam1)
    ka2 <- gam1 * ra2 * (ra2 * t_d2)^(gam1 - 1) * exp(-(ra2 * t_d2)^gam1)

    # 5. Ritonavir inverse Gaussian input (Yu 2020 Eq. 3):
    #    Nin(t) = Dose * sqrt(MAT / (2 pi CV^2 t^3)) * exp(-(t - MAT)^2 / (2 CV^2 MAT t)),
    # an inverse Gaussian density with mean MAT and shape lambda = MAT / CV^2.
    # It is applied as its hazard f(t) / S(t) on the gut amount, which returns
    # exactly Dose * f(t) after a single dose and conserves mass when doses
    # overlap. The survivor S(t) is the closed-form inverse Gaussian survivor
    # (phi() = standard normal CDF); the hazard tends to lambda / (2 MAT^2).
    lam_rtv <- mat_rtv / cvabs_rtv^2
    hzinf_rtv <- lam_rtv / (2 * mat_rtv^2)
    t_r1 <- max(tad0(depot1_rtv), 1e-6)
    pdf_r1 <- sqrt(lam_rtv / (2 * pi * t_r1^3)) * exp(-lam_rtv * (t_r1 - mat_rtv)^2 / (2 * mat_rtv^2 * t_r1))
    surv_r1 <- phi(-sqrt(lam_rtv / t_r1) * (t_r1 / mat_rtv - 1)) -
      exp(2 * lam_rtv / mat_rtv) * phi(-sqrt(lam_rtv / t_r1) * (t_r1 / mat_rtv + 1))
    hz_r1 <- hzinf_rtv
    if (surv_r1 > 1e-10) {
      hz_r1 <- pdf_r1 / surv_r1
    }
    t_r2 <- max(tad0(depot2_rtv), 1e-6)
    pdf_r2 <- sqrt(lam_rtv / (2 * pi * t_r2^3)) * exp(-lam_rtv * (t_r2 - mat_rtv)^2 / (2 * mat_rtv^2 * t_r2))
    surv_r2 <- phi(-sqrt(lam_rtv / t_r2) * (t_r2 / mat_rtv - 1)) -
      exp(2 * lam_rtv / mat_rtv) * phi(-sqrt(lam_rtv / t_r2) * (t_r2 / mat_rtv + 1))
    hz_r2 <- hzinf_rtv
    if (surv_r2 > 1e-10) {
      hz_r2 <- pdf_r2 / surv_r2
    }

    # 6. Ritonavir-inhibited intrinsic clearance and well-stirred extraction
    #    (Eqs. 1-3); ritonavir plasma concentration in ng/mL (mg/L * 1000).
    Cc_rtv <- central_rtv / vc_rtv * 1000
    clint <- clint0 - imax * Cc_rtv / (ki + Cc_rtv)
    eh <- clint * fu / (qh + clint * fu)

    # 7. TSP-1 turnover (Eqs. 4-6)
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout)
    ec50 <- exp(lec50)
    kin0 <- kout * rbase

    # 8. ODE system. Liver: inflow from the gut and from plasma (QH * Cc),
    #    outflow QH * C_liver of which the fraction EH is eliminated and
    #    (1 - EH) returns to the central compartment (Fig. 1).
    d/dt(depot1) <- -ka1 * depot1
    d/dt(depot2) <- -ka2 * depot2
    d/dt(liver) <- ka1 * depot1 + ka2 * depot2 + qh * central / vc - qh * liver / vh
    d/dt(central) <- qh * (1 - eh) * liver / vh - qh * central / vc -
      q * central / vc + q * peripheral1 / vp
    d/dt(peripheral1) <- q * central / vc - q * peripheral1 / vp
    d/dt(depot1_rtv) <- -hz_r1 * depot1_rtv
    d/dt(depot2_rtv) <- -hz_r2 * depot2_rtv
    d/dt(central_rtv) <- hz_r1 * depot1_rtv + hz_r2 * depot2_rtv -
      cl_rtv / vc_rtv * central_rtv - q_rtv / vc_rtv * central_rtv +
      q_rtv / vp_rtv * peripheral1_rtv
    d/dt(peripheral1_rtv) <- q_rtv / vc_rtv * central_rtv - q_rtv / vp_rtv * peripheral1_rtv

    Cc <- central / vc * 1000
    kin <- kin0 * (1 + Cc / (ec50 + Cc))
    tsp1(0) <- rbase
    d/dt(tsp1) <- kin - kout * tsp1

    # 9. Bioavailability (relative gut bioavailability; the paclitaxel F_H is
    #    produced structurally by the liver compartment)
    f(depot1) <- fgut
    f(depot2) <- fgut * e_dose2_fdepot
    f(depot1_rtv) <- f_rtv
    f(depot2_rtv) <- f_rtv * f2_rtv

    # 10. Observations (Eq. 8)
    Cc ~ prop(propSd)
    Cc_rtv ~ prop(propSd_rtv)
    tsp1 ~ prop(propSd_tsp1)
  })
}
