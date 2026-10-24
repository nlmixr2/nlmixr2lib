Fukaya_2023_glucarpidase_methotrexate <- function() {
  description <- "Population PK-PD model of glucarpidase (carboxypeptidase G2, CPG2) rescue after high-dose methotrexate (MTX) in Japanese patients with delayed MTX excretion (Fukaya 2023, phase 2 study): one-compartment CPG2 PK with power functions of BSA, coupled to two-compartment MTX PK with first-order renal elimination plus Michaelis-Menten hydrolysis whose Vmax is proportional to the plasma CPG2 concentration (Vmax = alpha * [CPG2])."
  reference <- "Fukaya Y, Kimura T, Hamada Y, Yoshimura K, Hiraga H, Yuza Y, Ogawa A, Hara J, Koh K, Kikuta A, Koga Y, Kawamoto H. Development of a population pharmacokinetics and pharmacodynamics model of glucarpidase rescue treatment after high-dose methotrexate therapy. Front Oncol. 2023;13:1003633. doi:10.3389/fonc.2023.1003633"
  vignette <- "Fukaya_2023_glucarpidase"
  units <- list(
    time = "h",
    dosing = "umol (methotrexate, into central); mg (glucarpidase protein mass, into central_cpg2)",
    concentration = "umol/L (methotrexate, Cc); mg/L (glucarpidase, Cc_cpg2)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters linearly: CLr = tvCLr * WT / CREAT / 60, Vc = tvVc * WT / BSA, Vp = tvVp * WT / 60 (Fukaya 2023 Section 3.2.2 and Table 4). The paper states 60 kg is the standardising weight (58.4 kg being a representative Japanese body weight).",
      source_name = "Body Weight"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Un-normalised power function on CPG2 CL and V (tvX * BSA^theta) and a divisor of MTX Vc (tvVc * WT / BSA). The BSA formula is not named in the paper.",
      source_name = "BSA"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Divisor of MTX renal clearance (tvCLr * WT / Scr / 60). The paper does not state whether serum creatinine was time-varying in the fit; it reports 1.23 +/- 0.89 mg/dL (mean +/- SD) immediately before CPG2 and 1.87 +/- 1.78 mg/dL on day 4 after CPG2 (Discussion).",
      source_name = "Scr"
    )
  )

  compartmentData <- list(
    central = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "umol", specimen = "tissue", verified = TRUE),
    central_cpg2 = list(analyte = "glucarpidase", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 15L,
    n_studies = 1L,
    age_range = "1-75 years",
    age_median = "15 years",
    weight_range = "10.7-78.1 kg",
    weight_median = "47.0 kg",
    height_range = "78.5-177.5 cm (median 156.0 cm)",
    sex_female_pct = 40,
    race_ethnicity = c(Asian = 100),
    disease_state = "Patients with delayed methotrexate excretion after high-dose MTX (>= 1 g/m^2): osteosarcoma (9), acute lymphocytic leukemia (3), non-Hodgkin lymphoma (2), medulloblastoma (1).",
    dose_range = "MTX 1.6-20.0 g (2.9-14.3 g/m^2); glucarpidase 50 U/kg as a 5-min IV infusion within 12 h of confirmed delayed excretion, with a second 50 U/kg dose when plasma MTX stayed > 1 umol/L more than 46 h after the first.",
    regions = "Japan (eight medical facilities)",
    renal_function = "Serum creatinine median 0.81 mg/dL (range 0.23-3.47) at baseline; MTX-related renal dysfunction in 8 of 14 patients.",
    notes = "Phase 2 multicentre single-arm study (JMA-IIA00097), December 2012 to February 2016; 88 plasma CPG2 (ELISA) and 115 plasma MTX (HPLC) concentrations. CPG2 popPK was fitted first and its post-hoc CL and V drove the CPG2 concentration in the MTX popPK-PD step (sequential fit; the two stages are combined here in one simulation model). Demographics from Fukaya 2023 Table 2."
  )

  ini({
    # Glucarpidase (CPG2) PK, phase 2 final model
    lcl_cpg2 <- log(0.238); label("Glucarpidase clearance at BSA = 1 m^2 (L/h)") # Table 4, phase 2, tvCLCPG2 = 0.238 L/h
    lvc_cpg2 <- log(1.200); label("Glucarpidase volume of distribution at BSA = 1 m^2 (L)") # Table 4, phase 2, tvVCPG2 = 1.200 L
    e_bsa_cl_cpg2 <- 1.440; label("Power exponent of BSA on glucarpidase clearance (unitless)") # Table 4, phase 2, theta1 = 1.440
    e_bsa_vc_cpg2 <- 1.561; label("Power exponent of BSA on glucarpidase volume (unitless)") # Table 4, phase 2, theta2 = 1.561

    # Methotrexate PK-PD, phase 2 final model
    lcl_renal <- log(3.248); label("MTX renal clearance per (WT/Scr/60) (L/h)") # Table 4, tvCLrMTX = 3.248 L/h
    lvc <- log(0.386); label("MTX central volume per (WT/BSA) (L)") # Table 4, tvVcMTX = 0.386 L
    lvp <- log(3.052); label("MTX peripheral volume per (WT/60) (L)") # Table 4, tvVpMTX = 3.052 L
    lq <- fixed(log(0.0778)); label("MTX intercompartmental clearance (L/h)") # Table 4, Q = 0.0778 L/h, held constant at the Fukuhara 2008 Japanese MTX value (Section 3.2.2)
    lalpha <- log(6.545e5); label("Conversion constant alpha, Vmax = alpha * [CPG2] (umol/h per mg/L of glucarpidase)") # Table 4, tv alpha = 6.545 x 10^5 (reported as L/h)
    km_cpg2 <- fixed(86); label("Michaelis-Menten constant of MTX hydrolysis by glucarpidase (umol/L)") # Table 4, Km = 86 umol, held constant; Section 3.2.2 equation, 86 umol/L (EMA Voraxaze assessment report)

    # IIV reported as CV%; omega^2 = log(CV^2 + 1)
    etalcl_cpg2 ~ 0.029827 # Table 4, phase 2, CL CV 17.4%
    etalvc_cpg2 ~ 0.047686 # Table 4, phase 2, V CV 22.1%
    etalcl_renal ~ 0.10636 # Table 4, CLrMTX CV 33.5%
    etalvc ~ 0.081286 # Table 4, VcMTX CV 29.1%
    etalvp ~ 0.59930 # Table 4, VpMTX CV 90.6%
    etalalpha ~ 0.49275 # Table 4, alpha CV 79.8%

    propSd <- 1.414; label("MTX proportional residual error (fraction)") # Table 4, popPK-PD residual variability 1.414; Section 2.6, Cij = Cij* x (1 + eps)
    addSd_cpg2 <- 0.100; label("Glucarpidase additive residual error (mg/L)") # Table 4, phase 2, residual variability 0.100 mg/L
  })

  model({
    # Glucarpidase disposition (one compartment, first-order elimination)
    cl_cpg2 <- exp(lcl_cpg2 + etalcl_cpg2) * BSA^e_bsa_cl_cpg2
    vc_cpg2 <- exp(lvc_cpg2 + etalvc_cpg2) * BSA^e_bsa_vc_cpg2
    kel_cpg2 <- cl_cpg2 / vc_cpg2

    # Methotrexate disposition
    cl_renal <- exp(lcl_renal + etalcl_renal) * WT / CREAT / 60
    vc <- exp(lvc + etalvc) * WT / BSA
    vp <- exp(lvp + etalvp) * WT / 60
    q <- exp(lq)
    alpha <- exp(lalpha + etalalpha)

    kr <- cl_renal / vc
    k12 <- q / vc
    k21 <- q / vp

    Cc_cpg2 <- central_cpg2 / vc_cpg2
    Cc <- central / vc

    # MTX hydrolysis by CPG2: Vmax = alpha * [CPG2] (Section 3.2.2 equation)
    vmax <- alpha * Cc_cpg2

    d/dt(central) <- -(kr + k12) * central + k21 * peripheral1 - vmax * Cc / (km_cpg2 + Cc)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_cpg2) <- -kel_cpg2 * central_cpg2

    Cc ~ prop(propSd)
    Cc_cpg2 ~ add(addSd_cpg2)
  })
}
