Wang_2020_sunitinib_ast <- function() {
  description <- "Indirect-response PK-PD model for aspartate aminotransferase (AST) in children and young adults (2-21 years) with refractory solid tumours receiving oral sunitinib, with a linear (first-order kPD) effect of plasma sunitinib concentration inhibiting the first-order loss rate kout (Wang 2020). The upstream sunitinib PK layer is the Wang 2020 two-compartment model with first-order absorption, lag time and BSA effects on CL/F and Vc/F, held fixed at its final estimates (sequential PK-PD)."
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
    ast = list(analyte = "aspartate aminotransferase (AST)", units = "U/L", specimen = "serum", verified = TRUE)
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
    notes = "Baseline demographics from Wang 2020 Table 2. The PK-PD analysis used the same 59 patients; the AST endpoint was modelled with the final PK model predictions of sunitinib concentration."
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
    lrbase <- log(26.2); label("Baseline AST BASE (U/L)") # Table 4 'Aspartate transaminase' block, row 'BASE' = 26.2 (RSE 8.3%)
    lkout <- log(1.7); label("First-order loss rate constant kout (1/h)") # Table 4 'Aspartate transaminase' block, row 'kout' = 1.7 (RSE 703%)
    lslope <- log(0.00492); label("Linear drug-effect coefficient kPD inhibiting kout (mL/ng)") # Table 4 'Aspartate transaminase' block, row 'kPD' = 0.00492 (RSE 32.7%)

    # IIV: Table 4 reports omega as CV%; omega^2 = log(1 + CV^2)
    etalrbase ~ 0.08344 # Table 4 'Aspartate transaminase' block, row 'omega (BASE)' = 29.5%
    etalkout ~ 2.992 # Table 4 'Aspartate transaminase' block, row 'omega (kout)' = 435%
    etalslope ~ 0.0004409 # Table 4 'Aspartate transaminase' block, row 'omega (kPD)' = 2.1%

    propSd <- 0.312; label("Proportional residual error on AST (fraction)") # Table 4 'Aspartate transaminase' block, row 'sigma' = 31.2% (RSE 2.7%)
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

    # ---- PD: indirect response ----
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout + etalkout)
    kin <- rbase * kout
    slope <- exp(lslope + etalslope)

    # Drug inhibits the loss rate: kout * (1 - kPD * Cc)
    d/dt(ast) <- kin - kout * (1 - slope * Cc) * ast
    ast(0) <- rbase

    ast ~ prop(propSd)
  })
}
