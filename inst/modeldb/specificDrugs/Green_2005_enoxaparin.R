Green_2005_enoxaparin <- function() {
  description <- "Two-compartment, first-order absorption population PK model for subcutaneously administered enoxaparin (anti-Xa activity) in 38 adults with acute coronary syndromes and a wide range of renal function (Green 2005). Total clearance is the sum of a renal arm scaled linearly to estimated creatinine clearance (CRCL, Cockcroft-Gault with ideal body weight; reference 80 mL/min) and a covariate-free non-renal arm: CL = 0.681 * (CRCL / 80) + 0.229 L/h. Central volume of distribution scales linearly with total body weight (reference 80 kg): Vc = 5.22 * (WT / 80) L. A constant basal anti-Xa activity (49.9 IU/L) is added to the model prediction to represent endogenous and assay-baseline anti-Xa activity, per the Schoemaker parameterisation referenced in the paper. Inter-individual variability is log-normal on total CL, Vc, Q, and basal anti-Xa activity (paper Table 2 Covariate Model). Residual error is combined additive (52.4 IU/L) plus proportional (20.0 percent CV) on observed anti-Xa concentrations."
  reference <- paste(
    "Green B, Greenwood M, Saltissi D, Westhuyzen J, Kluver L, Rowell J,",
    "Atherton J. (2005). Dosing strategy for enoxaparin in patients with",
    "renal impairment presenting with acute coronary syndromes.",
    "British Journal of Clinical Pharmacology 59(3):281-290.",
    "doi:10.1111/j.1365-2125.2004.02253.x"
  )
  vignette <- "Green_2005_enoxaparin"
  units <- list(time = "hour", dosing = "IU", concentration = "IU/L")
  # Enoxaparin dose is expressed in anti-Xa International Units (IU) to
  # match the assay readout used to fit the model (anti-Xa IU/L). The
  # clinical dose (mg) is converted via the standard 100 IU per 1 mg
  # enoxaparin equivalence (Lovenox / Clexane prescribing information).

  covariateData <- list(
    WT = list(
      description        = "Total body weight at study enrolment",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear (not allometric) covariate on Vc per Green 2005 Table 2 Covariate Model: Vc = 5.22 * (WT / 80) L. Reference 80 kg approximates the cohort median (69 kg, range 32-95). Source column WT.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance, computed with ideal body weight (IBW) as the size descriptor (raw, NOT BSA-normalized).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The paper labels this quantity GFR but it is computed by the Cockcroft-Gault formula with IBW substituted for total body weight as the size descriptor (Green 2005 Methods 'Population analysis' and Discussion paragraph 5). IBW is computed as IBW_male = 50 + 2.3 * (HT_in - 60) and IBW_female = 45.5 + 2.3 * (HT_in - 60) per Green 2005 reference [19,20]. Stored under the canonical CRCL column per inst/references/covariate-columns.md, with the per-model description recording the IBW-based assay form (CRCL accepts raw mL/min when the source paper does not apply BSA normalization). Reference value 80 mL/min per Green 2005 Table 2 (renal arm normalisation) and Discussion paragraph 5 (US National Kidney Foundation DOQI threshold for 'normal' renal function). Linear additive effect on CL: CL = 0.681 * (CRCL / 80) + 0.229 L/h.",
      source_name        = "GFR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 38L,
    n_studies      = 1L,
    age_range      = "median 78 years (range 44-87)",
    weight_range   = "median 69 kg (range 32-95)",
    height_range   = "median 1.63 m (range 1.46-1.84)",
    bmi_range      = "median 25.3 kg/m^2 (range 14.1-34.1)",
    serum_creat_range = "median 0.13 mmol/L (range 0.05-0.25)",
    crcl_range     = "median 32 mL/min (range 16-117); 47.4 percent below 30 mL/min, 28.9 percent in 30-50 mL/min, 15.8 percent in 51-80 mL/min, 7.9 percent above 80 mL/min",
    sex_female_pct = 52.6,
    disease_state  = "Adults admitted with acute coronary syndrome to a tertiary Coronary Care Unit, recruited between January 2001 and July 2002. Final diagnoses: non-ST elevation MI 60.5 percent, unstable angina 26.3 percent, ST elevation MI 10.5 percent, non-cardiac pain 2.6 percent. Comorbidities included hypertension 76.3 percent, hypercholesterolaemia 65.8 percent, diabetes 21.1 percent.",
    dose_range     = "Subcutaneous enoxaparin (Aventis Pharma) every 12 hours in conjunction with 100-150 mg oral aspirin daily. Dosing was adapted prospectively: the first seven patients received 0.5 mg/kg every 12 h if CRCL 10-25 mL/min, 0.75 mg/kg every 12 h if CRCL 26-50 mL/min, and 1.0 mg/kg every 12 h if CRCL above 50 mL/min; subsequently all patients received 1.0 mg/kg every 12 h.",
    regions        = "Single-centre study at the Royal Brisbane and Women's Hospital Coronary Care Unit, Brisbane, Queensland, Australia.",
    n_observations = 313L,
    notes          = "Population analysis used NONMEM v5 (Globomax) with FOCE-I. 313 anti-Xa concentrations measured by automated chromogenic assay (ACL Futura Plus, IL Test Heparin kit). Sampling: pre-dose plus 1, 2, 3, 4, 8, 12 h after the first dose, then trough before each subsequent dose. Three of 41 enrolled patients were excluded (consent withdrawal, cardiac event, obstructive uropathy). Demographics in Table 1; covariates considered in the model build included total body weight, lean body weight, ideal body weight, adjusted body weight, allometric scaling of each, body surface area, body mass index, predicted normal weight, height, sex, and Cockcroft-Gault CRCL computed with each weight descriptor. Final retained covariates: total body weight on Vc (linear) and CRCL (with IBW size descriptor) on CL (linear additive); no covariate effects were retained on intercompartmental clearance Q or basal anti-Xa activity."
  )

  ini({
    # Structural parameters -- Green 2005 Table 2 'Covariate Model' column.
    # Reference subject: WT = 80 kg total body weight, CRCL = 80 mL/min
    # (Cockcroft-Gault using ideal body weight as the size descriptor).
    lcl_renal  <- log(0.681); label("Renal clearance at reference CRCL = 80 mL/min (CL_renal, L/h)")  # Green 2005 Table 2 Covariate Model: 'renal = 0.681 / 80 ml min-1 (GFR)' (SE 33.3 percent)
    lcl_nonren <- log(0.229); label("Non-renal clearance (CL_nonrenal, L/h)")                          # Green 2005 Table 2 Covariate Model: 'nonrenal = 0.229' (SE 49.8 percent)
    lvc        <- log(5.22);  label("Central volume of distribution at reference WT = 80 kg (Vc, L)")  # Green 2005 Table 2 Covariate Model: '5.22 / 80 kg (WT)' (SE 18.8 percent)
    lka        <- log(0.255); label("First-order absorption rate constant (Ka, 1/h)")                  # Green 2005 Table 2 Covariate Model: 'K_a = 0.255' (SE 16.2 percent)
    lvp        <- log(29.6);  label("Peripheral volume of distribution (Vp, L)")                       # Green 2005 Table 2 Covariate Model: 'V_p = 29.6' (SE 22.6 percent)
    lq         <- log(0.632); label("Inter-compartmental clearance (Q, L/h)")                          # Green 2005 Table 2 Covariate Model: 'Q = 0.632' (SE 21.2 percent)
    lrbase     <- log(49.9);  label("Basal anti-Xa activity (IU/L)")                                    # Green 2005 Table 2 Covariate Model: 'Basal anti-Xa activity = 49.9' (SE 30.1 percent)

    # Inter-individual variability. Paper Results paragraph 2: "log normal
    # between subject variability (BSV) on clearance (CL), central volume
    # compartment (Vc), intercompartmental clearance (Q) and basal anti-Xa
    # activity". Table 2 reports BSV magnitudes as percent CV. Without an
    # explicit conversion in the paper, the formal log-normal back-transform
    # omega^2 = log(1 + CV^2) is used (matches the convention adopted in
    # Pierre_2017_morphine.R for the same omega-as-CV reporting style).
    #   CL    : 32.7 percent CV -> omega^2 = log(1 + 0.327^2) = 0.10162
    #   Vc    : 34.4 percent CV -> omega^2 = log(1 + 0.344^2) = 0.11186
    #   Q     : 69.8 percent CV -> omega^2 = log(1 + 0.698^2) = 0.39690
    #   Basal : 76.6 percent CV -> omega^2 = log(1 + 0.766^2) = 0.46168
    # No correlations between etas are reported, so each is independent.
    etalcl    ~ log(1 + 0.327^2)  # Green 2005 Table 2 Covariate Model: omega_CL    = 32.7 percent CV (SE 62.3 percent)
    etalvc    ~ log(1 + 0.344^2)  # Green 2005 Table 2 Covariate Model: omega_Vc    = 34.4 percent CV (SE 72.2 percent)
    etalq     ~ log(1 + 0.698^2)  # Green 2005 Table 2 Covariate Model: omega_Q     = 69.8 percent CV (SE 41.9 percent)
    etalrbase ~ log(1 + 0.766^2)  # Green 2005 Table 2 Covariate Model: omega_Basal = 76.6 percent CV (SE 30.2 percent)

    # Residual error -- combined additive + proportional on observed anti-Xa
    # activity. Table 2 reports the additive component as an SD in IU/L and
    # the proportional component as percent CV; nlmixr2 expects SDs.
    addSd  <- 52.4; label("Additive residual SD on anti-Xa activity (IU/L)")     # Green 2005 Table 2 Covariate Model: 'sigma_1 = 52.4 IU/L' (SE 37.5 percent)
    propSd <- 0.20; label("Proportional residual SD on anti-Xa activity (fraction)")  # Green 2005 Table 2 Covariate Model: 'sigma_2 = 20.0 percent CV' (SE 35.6 percent)
  })
  model({
    # Individual structural parameters. Renal and non-renal clearance
    # arms combine additively, with the renal arm scaling linearly to
    # CRCL: CL = 0.681 * (CRCL / 80) + 0.229 L/h at typical level. The
    # single eta on CL applies to total CL per the paper's BSV
    # description (Results paragraph 2: 'log normal BSV on CL'),
    # matching the NONMEM coding TVCL = THETA_renal * (CRCL/80) +
    # THETA_nonren; CL = TVCL * EXP(ETA_CL).
    cl_renal  <- exp(lcl_renal) * (CRCL / 80)
    cl_nonren <- exp(lcl_nonren)
    cl        <- (cl_renal + cl_nonren) * exp(etalcl)

    # Central volume scales linearly with total body weight (NOT
    # allometric); the paper retained the linear form preferentially
    # over allometric scaling during covariate selection (Methods
    # 'Population analysis').
    vc <- exp(lvc + etalvc) * (WT / 80)
    vp <- exp(lvp)
    q  <- exp(lq + etalq)
    ka <- exp(lka)

    # Basal anti-Xa activity (constant offset added to the model
    # prediction): a measurement / endogenous baseline present in
    # plasma prior to dosing per the Schoemaker reference [29] cited
    # in the paper, with log-normal BSV across subjects.
    rbase <- exp(lrbase + etalrbase)

    # ODE system: first-order absorption from depot into central, with
    # two-compartment disposition.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka  * depot
    d/dt(central)     <-  ka  * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observed anti-Xa activity is the simulated enoxaparin concentration
    # plus the basal anti-Xa activity offset (IU/L). Dose in mg, vc in
    # L; the resulting model concentration is on the same scale as the
    # anti-Xa activity measurement used to fit the model (IU/L), with
    # the conversion absorbed into the structural parameters per the
    # paper's parameterisation.
    Cc <- central / vc + rbase

    Cc ~ add(addSd) + prop(propSd)
  })
}
