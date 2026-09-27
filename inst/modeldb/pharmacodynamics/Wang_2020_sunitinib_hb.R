Wang_2020_sunitinib_hb <- function() {
  description <- "Semi-mechanistic transit-compartments-in-series-with-feedback-loop (Friberg-type) PK-PD model for hemoglobin in children and young adults (2-21 years) with refractory solid tumours receiving oral sunitinib: a proliferating pool, three transit compartments and a circulating pool with the feedback exponent fixed to 1, and a linear (first-order kPD) effect of plasma sunitinib concentration inhibiting proliferation (Wang 2020). The upstream sunitinib PK layer is the Wang 2020 two-compartment model with first-order absorption, lag time and BSA effects on CL/F and Vc/F, held fixed at its final estimates (sequential PK-PD)."
  reference <- paste(
    "Wang E, DuBois SG, Wetmore C, Khosravan R.",
    "Population pharmacokinetics-pharmacodynamics of sunitinib in pediatric",
    "patients with solid tumors.",
    "Cancer Chemother Pharmacol. 2020;86(2):181-192.",
    "doi:10.1007/s00280-020-04106-z.",
    "PD model structures (Figure 5) from Khosravan R, Motzer RJ, Fumagalli E,",
    "Rini BI. Clin Pharmacokinet. 2016;55(10):1251-1269.",
    "doi:10.1007/s40262-016-0404-5.",
    sep = " "
  )
  vignette <- "Wang_2020_sunitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BSA = list(
      description = "Baseline body surface area (DuBois and DuBois formula)",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters only through the upstream sunitinib PK layer (linear effect on CL/F, power effect on Vc/F; reference 1.47 m^2). No covariates were retained on the PD parameters (Wang 2020 Results).",
      source_name = "BSA"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "sunitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sunitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sunitinib", units = "mg", specimen = "plasma", verified = TRUE),
    prol = list(
      analyte = "hemoglobin precursor (proliferating pool)",
      units = "g/dL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "hemoglobin precursor (transit 1)",
      units = "g/dL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "hemoglobin precursor (transit 2)",
      units = "g/dL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "hemoglobin precursor (transit 3)",
      units = "g/dL",
      specimen = "not applicable",
      verified = TRUE
    ),
    circ = list(analyte = "hemoglobin", units = "g/dL", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 59L,
    n_studies = 2L,
    age_range = "2-21 years",
    weight_range = "16.2-100 kg (median 50.4 kg)",
    bsa_range = "0.66-2.14 m^2 (median 1.47 m^2)",
    sex_female_pct = 52.5,
    race_ethnicity = c(Asian = 5.1, NonAsian = 89.8, Unknown = 5.1),
    disease_state = "Children and young adults with refractory solid tumours, predominantly high-grade glioma, ependymoma, brain stem glioma, or sarcoma (Children's Oncology Group studies ADVL0612 and ACNS1021).",
    dose_range = "Sunitinib 15 or 20 mg/m^2 orally once daily on schedule 4/2 (4 weeks on, 2 weeks off).",
    regions = "United States and Canada (Children's Oncology Group).",
    notes = "Baseline demographics from Wang 2020 Table 2. The PK-PD analysis used the same 59 patients; the hemoglobin endpoint was modelled with the final PK model predictions of sunitinib concentration."
  )

  ini({
    # ---- Upstream sunitinib PK (Wang 2020 Table 3, sunitinib final model) ----
    # Held fixed: the PK-PD models were fitted sequentially on the final PK
    # model predictions (Wang 2020 Methods, Model development).
    lka <- fixed(log(0.38)); label("Sunitinib absorption rate constant ka (1/h)") # Table 3 sunitinib 'ka' = 0.38
    ltlag <- fixed(log(0.64)); label("Sunitinib absorption lag time tlag (h)") # Table 3 sunitinib 'tlag' = 0.64
    lcl <- fixed(log(24.1)); label("Sunitinib apparent clearance CL/F at BSA 1.47 m^2 (L/h)") # Table 3 sunitinib 'CL/F' = 24.1
    lvc <- fixed(log(1070)); label("Sunitinib apparent central volume Vc/F at BSA 1.47 m^2 (L)") # Table 3 sunitinib 'Vc/F' = 1070
    lvp <- fixed(log(63.8)); label("Sunitinib apparent peripheral volume Vp/F (L)") # Table 3 sunitinib 'Vp/F' = 63.8
    lq <- fixed(log(0.28)); label("Sunitinib apparent intercompartmental clearance Q/F (L/h)") # Table 3 sunitinib 'Q/F' = 0.28
    e_bsa_cl <- fixed(0.557); label("Linear slope of BSA on CL/F, per m^2 about 1.47 m^2 (unitless)") # Results text CL/F = 24.1 * [1 + 0.557 * (BSA - 1.47)]
    e_bsa_vc <- fixed(1.47); label("Power exponent of BSA/1.47 on Vc/F (unitless)") # Results text Vc/F = 1070 * (BSA/1.47)^1.47
    etalcl ~ fixed(0.1106) # Table 3 sunitinib omega (CL/F) = 34.2%, omega^2 = log(1 + CV^2)
    etalvc ~ fixed(0.05646) # Table 3 sunitinib omega (Vc/F) = 24.1%
    etalka ~ fixed(0.5705) # Table 3 sunitinib omega (ka) = 87.7%

    # ---- PD parameters (Wang 2020 Table 4) ----
    lrbase <- log(13); label("Baseline hemoglobin BASE (g/dL)") # Table 4 'Hemoglobin' block, row 'BASE' = 13.0 (RSE 1.8%); printed unit '10^9/L' is a typo, 13.0 is a hemoglobin concentration in g/dL
    lmtt <- log(1370); label("Mean transit time MTT from the proliferating pool to the circulation (h)") # Table 4 'Hemoglobin' block, row 'MTT' = 1370 (RSE 13%)
    lslope <- log(0.000317); label("Linear drug-effect coefficient kPD on proliferation (mL/ng)") # Table 4 'Hemoglobin' block, row 'kPD' = 0.000317 (RSE 59%)
    lgamma <- fixed(log(1)); label("Feedback exponent POW on (BASE / circ) (unitless)") # Table 4 'Hemoglobin' block, row 'POW' = Fixed to 1

    # IIV: Table 4 reports omega as CV%; omega^2 = log(1 + CV^2)
    etalrbase ~ 0.01453 # Table 4 'Hemoglobin' block, row 'omega (BASE)' = 12.1%
    etalslope ~ 2.062 # Table 4 'Hemoglobin' block, row 'omega (kPD)' = 262%

    propSd <- 0.058; label("Proportional residual error on hemoglobin (fraction)") # Table 4 'Hemoglobin' block, row 'sigma' = 5.8% (RSE 2.7%)
  })

  model({
    # ---- Upstream sunitinib PK ----
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)
    cl <- exp(lcl + etalcl) * (1 + e_bsa_cl * (BSA - 1.47))
    vc <- exp(lvc + etalvc) * (BSA / 1.47)^e_bsa_vc
    vp <- exp(lvp)
    q <- exp(lq)
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    # Sunitinib plasma concentration (mg / L * 1000 = ng/mL); drives the PD effect
    Cc <- central / vc * 1000

    # ---- PD: transit compartments in series with feedback loop ----
    # Khosravan 2016 Figure 5A: Kprol = Ktr = Kcirc, three transit compartments.
    # ktr = (n + 1) / MTT with n = 3 transit compartments (Friberg 2002).
    rbase <- exp(lrbase + etalrbase)
    mtt <- exp(lmtt)
    gamma <- exp(lgamma)
    ktr <- 4 / mtt
    slope <- exp(lslope + etalslope)
    edrug <- slope * Cc
    feed <- (rbase / circ)^gamma

    d/dt(prol) <- ktr * prol * (1 - edrug) * feed - ktr * prol
    d/dt(transit1) <- ktr * prol - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(circ) <- ktr * transit3 - ktr * circ

    prol(0) <- rbase
    transit1(0) <- rbase
    transit2(0) <- rbase
    transit3(0) <- rbase
    circ(0) <- rbase

    circ ~ prop(propSd)
  })
}
