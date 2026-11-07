Zhang_2025_obesity_renal_clearance <- function() {
  description <- paste(
    "PBPK (renal clearance, hypothetical drug). Renal clearance (CLr) of a",
    "hypothetical small-molecule drug eliminated by glomerular filtration",
    "(GF) and active tubular secretion (ATS) in adults with normal weight,",
    "overweight or obesity (BMI 20 to 60 kg/m^2; Zhang 2025). CLr = CLGF +",
    "CLATS, with CLGF = fu * GFR and CLATS = (Qr - GFR) * fu * CLint_sec /",
    "(Qr + fu * CLint_sec / BP). Renal blood flow, GFR, kidney weight,",
    "albumin, alpha-1 acid glycoprotein and hematocrit follow the obese",
    "population regressions of Berton 2022 in BMI, height, age, sex and",
    "body surface area. fu scales with the binding protein (albumin or",
    "AAG) from its normal-weight value, BP from fu, hematocrit and the",
    "red-cell partition coefficient Kp, and CLint_sec = CLint_ATS * PTCPGK",
    "* kidney weight * relative transporter activity (rTA). The drug",
    "properties (normal-weight fu, Kp, CLint_ATS, binding protein) and rTA",
    "are fixed placeholders that the user replaces. The model has no",
    "compartments, ODEs, dosing, between-subject variability or residual",
    "error; each record gives the clearances implied by its covariates and",
    "the time column is not used.",
    sep = " "
  )
  reference <- paste(
    "Zhang T, Calvier EAM, Krekels EHJ, Knibbe CAJ. The Impact of Drug",
    "Properties and Severity of Obesity on Renal Drug Clearance Through",
    "Glomerular Filtration and Active Tubular Secretion: A Systematic",
    "Analysis Using PBPK Modeling. Pharm Res. 2025;42(7):1079-1088.",
    "doi:10.1007/s11095-025-03885-5. PMID 40571872. PMCID PMC12304071.",
    "Eqs. 1-3 of the main text; the subject-specific equations and the",
    "drug- and system-specific Eqs. 1-3 are from the Supplementary",
    "Material (sections 1 and 2). The blood-weight and heart-rate",
    "functions inside the printed cardiac-output equation are from",
    "Berton M, Bettonte S, Stader F, Battegay M, Marzolini C. Clin",
    "Pharmacokinet. 2022;61(9):1251-1270 (doi:10.1007/s40262-022-01132-3),",
    "Table 1.",
    sep = " "
  )
  vignette <- "Zhang_2025_obesity_renal_clearance"
  units <- list(
    time = paste(
      "n/a -- the model has no dynamics. Each record gives the renal",
      "clearance implied by its covariates, so the time column is ignored."
    ),
    dosing = "n/a (renal clearance model; no dosing)",
    concentration = paste(
      "n/a (no drug concentration). Outputs are clearances in mL/min (clr,",
      "cl_gf, cl_ats, clint_sec), the percent contributions pct_gf and",
      "pct_ats, and the system quantities gfr (mL/min), qr and co (L/min),",
      "kw and bw (kg), hsa and aag (g/L), hct, fu and bpr (fractions)."
    )
  )

  covariateData <- list(
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives body weight (BW = BMI * (HT/100)^2), fat mass, the ejection",
        "fraction in cardiac output, albumin and AAG (Supplementary section",
        "1). The paper studied typical subjects with BMI 20 (the",
        "normal-weight reference), 25, 30, 40, 50 and 60 kg/m^2; the Berton",
        "2022 regressions were developed for BMI 18.5 to 60 kg/m^2."
      ),
      source_name = "BMI"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Converts BMI to body weight and enters the fat-free mass index. The",
        "paper used 172 cm for every typical subject (the NHANES adult",
        "median)."
      ),
      source_name = "HT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the kidney-weight and GFR regressions. The paper set age to",
        "50 years (Supplementary section 1 note)."
      ),
      source_name = "Age"
    ),
    SEXF = list(
      description = "Sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Berton 2022 codes Sex as 0 = male, 1 = female, the same as SEXF, so",
        "no transformation is needed. Sex enters every regression linearly,",
        "and the paper used the average of the male and female values",
        "(Supplementary section 1 note); SEXF = 0.5 reproduces that",
        "sex-averaged typical subject exactly."
      ),
      source_name = "Sex"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the kidney-weight regression only. Neither paper prints the",
        "BSA formula used for the typical subjects; Berton 2022 states that",
        "it uses the Ashby-Thompson equation. BSA is therefore taken as a",
        "covariate so that the user's own formula is used; the vignette",
        "shows that the sex-averaged Ashby-Thompson BSA reproduces the",
        "printed kidney weights."
      ),
      source_name = "Body surface area"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    age_range = "50 years (fixed for all typical subjects)",
    weight_range = "59.2 to 177.5 kg (BMI 20 to 60 kg/m^2 at a height of 172 cm)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not applicable (Berton 2022 regressions were built from White adults)",
    disease_state = "Normal weight, overweight and (morbid) obesity; no renal impairment",
    dose_range = "Not applicable (no dosing)",
    regions = "Not applicable",
    notes = paste(
      "No patients were studied. Six typical subjects (BMI 20, 25, 30, 40,",
      "50 and 60 kg/m^2; height 172 cm; age 50 years; sex-averaged) were",
      "combined with 11,520 hypothetical drugs: normal-weight fu 5, 25, 50,",
      "75, 95 or 100 percent, binding to albumin or AAG, Kp 0.35, 1, 2, 3,",
      "4, 5, 10, 20, 30 or 40, and 32 equidistant CLint_ATS values from 2",
      "to 500 uL/min/mg protein, each under rTA 20, 50, 100, 150, 200 or",
      "250 percent (Methods)."
    )
  )

  ini({
    # ---- System constants ----
    # PTCPGK as printed (Supplement section 2: 'assumed to be 99.4,
    # according to the SimCYP Simulator (version 22)'). The units are not
    # printed. With CLint_ATS in uL/min/mg protein and kidney weight in kg,
    # CLint_ATS * PTCPGK * KW is in mL/min, which reproduces Table S2.
    ptcpgk <- fixed(99.4)
    label("Proximal tubule cells per gram kidney, PTCPGK (as printed; see the vignette)") # Supplement section 2: PTCPGK 99.4 (Simcyp V22), unaltered across BMI
    hr <- fixed(80)
    label("Heart rate (beats/min)") # Berton 2022 Table 1: HR truncated normal, mu = 80 (min 60, max 100); heart rate does not change with obesity
    bmi_nw <- fixed(20)
    label("BMI of the normal-weight reference subject for the fu scaling (kg/m^2)") # Methods: normal-weight subject BMI 20 kg/m2; Supplement Eq. 1 PP normal weight

    # ---- Hypothetical-drug properties (placeholders) ----
    # Methods list the values the paper swept; the defaults are one of
    # them. Replace them per record for another drug or scenario.
    fu_nw <- fixed(0.5)
    label("Fraction unbound in plasma in a normal-weight subject (fraction)") # Methods: fu normal weight 5, 25, 50, 75, 95 or 100 percent
    kp_rbc <- fixed(1)
    label("Red blood cell to plasma partition coefficient, Kp (unitless)") # Methods: Kp 0.35, 1, 2, 3, 4, 5, 10, 20, 30 or 40; Figs. 3 and 4 use Kp = 1
    clint_ats <- fixed(50)
    label("Intrinsic active tubular secretion clearance, CLint_ATS (uL/min/mg protein)") # Methods: 2 to 500 uL/min/mg protein, 32 equidistant values; 50 is the Fig. 3 intermediate value
    bind_aag <- fixed(0)
    label("Binding protein: 0 = human serum albumin, 1 = alpha-1 acid glycoprotein (unitless)") # Methods: drugs bind exclusively to either HSA or AAG
    rta <- fixed(1)
    label("Relative transporter activity versus normal weight, rTA (fraction)") # Methods: rTA 20, 50, 100, 150, 200 or 250 percent; Supplement section 2: integrated as fractions
  })

  model({
    # ---- Subject-specific physiology (Supplement section 1; Berton 2022) ----
    bw <- BMI * (HT / 100)^2
    # Fat-free mass index (kg/m^2) and fat mass (kg). The printed bracket
    # reads [FFMI - BMI]; the printed fat masses (14.1 kg at BMI 20, 83.4 kg
    # at BMI 60) need BMI - FFMI.
    ffmi <- 30.86 - 3.09 * SEXF - 0.13 * HT + 0.14 * bw
    fm <- (BMI - ffmi) * (HT / 100)^2
    ffat <- fm / bw
    # Cardiac output (L/min) = HR x EDV x EF / 1000, EDV = 43.13 + 12.96 x
    # blood weight. Blood weight is Berton 2022 Table 1 with the adipose
    # fraction as a fraction, as in the printed GFR equation; this gives
    # the printed cardiac outputs (5.86 to 8.96 L/min).
    blood_wt <- 2.58 + 0.057 * bw - 0.061 * ffat
    co <- hr * (43.13 + 12.96 * blood_wt) * (0.65 - 0.0018 * BMI) / 1000
    # Renal blood flow (L/min); the renal share of cardiac output is in percent.
    qr <- (20.57 - 1.76 * SEXF) * co / 100
    hsa <- 47.58 - 0.16 * BMI
    aag <- 0.48 + 0.012 * BMI
    hct <- 0.45 - 0.04 * SEXF
    # Kidney weight (kg), with the -0.048 x Sex coefficient printed in the
    # Supplement (Berton 2022 Table 1 prints -0.48).
    kw <- -0.11 - 0.048 * SEXF + 0.003 * AGE + 0.18 * BSA
    # GFR (mL/min); Qr enters in L/min (printed as Qr/1000 with Qr in mL/min).
    gfr <- 13.02 - 0.51 * AGE + 0.93 * ffat + 93.56 * qr

    # ---- Drug- and system-specific parameters (Supplement Eqs. 1-3) ----
    pp <- (1 - bind_aag) * hsa + bind_aag * aag
    pp_nw <- (1 - bind_aag) * (47.58 - 0.16 * bmi_nw) +
      bind_aag * (0.48 + 0.012 * bmi_nw)
    fu <- 1 / (1 + (1 - fu_nw) * pp / (fu_nw * pp_nw))
    bpr <- 1 + hct * (fu * kp_rbc - 1)
    # CLint_sec (mL/min): uL/min/mg protein x PTCPGK x kidney weight in kg.
    clint_sec <- clint_ats * ptcpgk * kw * rta

    # ---- Renal clearance (main text Eqs. 1-5), mL/min ----
    qr_ml <- qr * 1000
    cl_gf <- fu * gfr
    cl_ats <- (qr_ml - gfr) * fu * clint_sec / (qr_ml + fu * clint_sec / bpr)
    clr <- cl_gf + cl_ats
    pct_gf <- 100 * cl_gf / clr
    pct_ats <- 100 * cl_ats / clr
  })
}
