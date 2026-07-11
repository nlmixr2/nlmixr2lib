Freyer_2000_etoposide <- function() {
  description <- "One-compartment IV population PK model for etoposide in 24 small cell lung cancer patients on the AVI regimen (Freyer 2000), with a linear serum creatinine effect on clearance. Structural reduction from the paper's two-compartment model: only central CL and V were reported in Table 2, so the peripheral compartment is omitted."
  reference <- paste(
    "Freyer G, Tranchand B, Ligneau B, Ardiet C, Souquet P-J,",
    "Court-Fortune I, Riou R, Rebattu P, Boissel J-P,",
    "Trillet-Lenoir V, Girard P.",
    "Population pharmacokinetics of doxorubicin, etoposide and ifosfamide",
    "in small cell lung cancer patients: results of a multicentre study.",
    "Br J Clin Pharmacol. 2000;50(4):315-324.",
    "doi:10.1046/j.1365-2125.2000.00269.x.",
    sep = " "
  )
  vignette <- "Freyer_2000_AVI_regimen"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CREAT = list(
      description        = "Baseline serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear covariate effect on etoposide CL as reported in Freyer",
        "2000 Table 2 and Results Etoposide (Equation): CL = 3.34 -",
        "0.0083 * S_Cr with S_Cr in umol/L. The effect is additive on the",
        "linear-space CL, NOT multiplicative on log-CL. Cohort median",
        "SCR = 82 umol/L (Table 1) is used as the reference for the",
        "recentred parameterisation lcl = log(2.66) so that the typical",
        "CL at reference SCR equals exp(lcl). Cohort SCR range 44-144",
        "umol/L; do NOT extrapolate to SCR values that would drive the",
        "linear form 3.34 - 0.0083 * S_Cr non-positive (S_Cr > ~402",
        "umol/L). Source column named S_Cr in the paper."
      ),
      source_name        = "SCR"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 24L,
    n_studies        = 1L,
    n_courses        = 47L,
    age_range        = "45-70 years",
    age_median       = "57 years",
    weight_range     = "54-93 kg",
    weight_median    = "69 kg",
    height_range     = "159-179 cm",
    height_median    = "170 cm",
    bsa_range        = "1.54-2.10 m^2",
    bsa_median       = "1.79 m^2",
    creat_range      = "44-144 umol/L",
    creat_median     = "82 umol/L (Freyer 2000 Table 1; used as the linear-covariate recentring reference)",
    sex_female_pct   = "not reported",
    race_ethnicity   = "not reported (French multicentre study)",
    disease_state    = paste(
      "Small cell lung cancer (SCLC), either limited to the thorax or",
      "extensive; no severe hepatic or renal impairment (patients with",
      "ASAT/ALAT/alkaline phosphatase > 2x upper limit of normal and/or",
      "serum creatinine > 1.5x upper limit of normal were excluded from",
      "the primary therapeutic trials)."
    ),
    dose_range       = paste(
      "AVI regimen: etoposide 120 mg/m^2 IV over 30 min on days 1, 2, and",
      "3. Total etoposide dose per course = 360 mg/m^2. At the 1.79 m^2",
      "median BSA that is ~215 mg per day-1 infusion."
    ),
    regions          = "France (Lyon Saint-Etienne Thoracic Oncology Group, GLOT; multicentre)",
    covariates_tested = paste(
      "Screened: age, sex, stage of disease, body weight, height, body",
      "surface area, serum creatinine, ASAT, ALAT, alkaline phosphatase,",
      "gamma-GT, total bilirubin, LDH, total protein. Only serum",
      "creatinine (r = -0.37, decrease of 12 points in the NONMEM",
      "objective function on inclusion) was retained in the final model",
      "(Freyer 2000 Results Etoposide, Figure 6)."
    ),
    notes            = paste(
      "Baseline characteristics from Freyer 2000 Table 1. Etoposide",
      "assayed by HPLC-UV (LOQ 150 ng/mL, linear 50-10000 ng/mL, CV 3.5-",
      "10.2%). NONMEM V (FO method); extensive sampling (14 samples per",
      "course across days 1-3) or limited sampling (7 samples per course)."
    )
  )

  ini({
    # Structural PK (Freyer 2000 Table 2, etoposide row; Results Etoposide).
    # The paper reports CL as a LINEAR covariate model:
    #   CL_pop = 3.34 - 0.0083 * S_Cr           (Freyer 2000 Table 2 / Fig 6)
    # Recentred at cohort-median SCR = 82 umol/L (Freyer 2000 Table 1):
    #   CL_pop = 3.34 - 0.0083 * 82 = 2.66 L/h at reference SCR.
    lcl          <- log(2.66);   label("Typical CL at reference SCR = 82 umol/L (L/h)")   # Freyer 2000 Table 2: CL intercept 3.34 (s.d. 0.228); recentred by cohort median SCR (Table 1)
    e_creat_cl   <- -0.0083;     label("Linear SCR slope on CL (L/h per umol/L)")         # Freyer 2000 Table 2 & Results Equation: CL = 3.34 - 0.0083*S_Cr; slope treated as structural (no s.d. tabulated)
    creat_ref_cl <- fixed(82);   label("Reference SCR for the linear CL covariate (umol/L)")  # Freyer 2000 Table 1: cohort-median SCR
    lvc          <- log(6.4);    label("Typical central volume of distribution (L)")     # Freyer 2000 Table 2: V_d = 6.4 L (s.d. 0.863); Results Etoposide "V_d value was 6.38 l" (6.4 rounded)

    # IIV -- Freyer 2000 Table 2 reports NONMEM omega variances.
    etalcl ~ 0.0243  # Freyer 2000 Table 2: omega^2_CL = 0.0243 (s.d. 0.0107); sqrt(exp(0.0243)-1) = 15.7 pct CV, matches text 15.6 pct for CL
    etalvc ~ 0.0350  # Freyer 2000 Table 2: omega^2_Vd = 0.0350 (s.d. 0.0128); sqrt(exp(0.0350)-1) = 18.9 pct CV, matches text 18.7 pct for Vd

    # Residual error -- proportional only (Freyer 2000 Methods Data
    # analysis: "a proportional model was used for etoposide").
    propSd <- 0.2400; label("Proportional residual error (fraction)")  # Freyer 2000 Table 2: sigma^2_prop = 0.0576 (s.d. 0.0213) -> SD = sqrt(0.0576) = 0.2400 (24.0%)
  })

  model({
    # Individual PK parameters. The paper's covariate model is LINEAR in
    # SCR on linear-space CL (NOT multiplicative on log CL); log-normal IIV
    # is applied on top of the covariate-adjusted linear CL. Recentring at
    # SCR = 82 umol/L is algebraically identical to the paper's
    # CL = 3.34 - 0.0083 * S_Cr form.
    cl <- (exp(lcl) + e_creat_cl * (CREAT - creat_ref_cl)) * exp(etalcl)
    vc <- exp(lvc + etalvc)

    # Micro-constant
    kel <- cl / vc

    # ODE system -- 1-compartment IV (etoposide 30-min IV infusion, days 1,
    # 2, and 3). Dose is delivered directly into the central compartment.
    d/dt(central) <- -kel * central

    # Plasma concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
