Budha_2023_tislelizumab <- function() {
  description <- "Three-compartment population PK model for intravenous tislelizumab (anti-PD-1 IgG4) in patients with advanced tumors (Budha 2023)"
  reference <- "Budha N, Wu CY, Tang Z, et al. Model-based population pharmacokinetic analysis of tislelizumab in patients with advanced tumors. CPT Pharmacometrics Syst Pharmacol. 2023;12(1):95-109. doi:10.1002/psp4.12880"
  vignette <- "Budha_2023_tislelizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL and Vc with reference weight 65 kg (Budha 2023 Equations 5 and 6).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on Vc with reference age 60 years (Budha 2023 Equation 6).",
      source_name        = "AGE"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 41 g/L (Budha 2023 Equation 5). Budha 2023 reports albumin in g/L (SI), not g/dL.",
      source_name        = "ALB"
    ),
    TUMSZ = list(
      description        = "Baseline tumor size (sum of diameters for solid tumors; sum of products of perpendicular diameters for cHL)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 63 mm (Budha 2023 Equation 5; Table 2 median 63.3 mm).",
      source_name        = "TUMSZ"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Exponential effect on Vc for females relative to males (Budha 2023 Equation 6). Renamed from source column SEX (character 'Female'/'Male') to the canonical SEXF per covariate-columns.md.",
      source_name        = "SEX"
    ),
    ADA_POS = list(
      description        = "Antidrug-antibody status (treatment-emergent)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Exponential effect on CL for ADA-positive patients (Budha 2023 Equation 5). Renamed from source column ADA to the canonical ADA_POS per covariate-columns.md.",
      source_name        = "ADA"
    ),
    TUMTP_CHL = list(
      description        = "Tumor-type indicator for classical Hodgkin lymphoma",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all other tumor types: NSCLC, EC, HCC, UC, GC, CRC, NPC, OC, Other)",
      notes              = "Exponential effect on CL for cHL patients (Budha 2023 Equation 5). Derived from the source categorical column TUMTP as TUMTP_CHL = as.integer(TUMTP == 'cHL').",
      source_name        = "TUMTP"
    ),
    TUMTP_GC = list(
      description        = "Tumor-type indicator for gastric cancer",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all other tumor types)",
      notes              = "Exponential effect on CL for GC patients (Budha 2023 Equation 5). Derived from the source categorical column TUMTP as TUMTP_GC = as.integer(TUMTP == 'GC').",
      source_name        = "TUMTP"
    )
  )

  population <- list(
    n_subjects     = 2596,
    n_studies      = 12,
    age_range      = "18-90 years",
    age_median     = "60 years",
    weight_range   = "31.9-130 kg",
    weight_median  = "65 kg",
    sex_female_pct = 26.0,
    race_ethnicity = c(White = 20.3, Asian = 76.7, `Black/African American` = 0.4, Other = 1.7, Missing = 0.9),
    disease_state  = "Advanced / metastatic solid tumors or classical Hodgkin lymphoma (NSCLC 44.3%, EC 14.4%, HCC 12.2%, UC 5.8%, GC 3.9%, CRC 3.1%, cHL 2.7%, OC 2.0%, NPC 0.8%, Other 10.7%).",
    dose_range     = "0.5-10 mg/kg IV q2w or q3w, or 200 mg IV q3w (flat) across 12 studies",
    regions        = "Global (4 studies), China (7 studies), China/Korea (1 study)",
    notes          = "Baseline demographics and covariate distributions per Budha 2023 Tables 1 and 2. 14,473 serum concentration observations. ADA-positive 16.6%, ADA-negative 82.3%, missing 1.1%. ECOG PS 0: 31.5%, ECOG PS 1: 68.5%. Therapy: monotherapy 78.7%, combination 21.0%."
  )

  ini({
    # Structural parameters — typical values for a 65 kg male patient, 60 y, ALB 41 g/L,
    # TUMSZ 63 mm, ADA-negative, non-cHL/non-GC tumor (Budha 2023 Table 3).
    lcl  <- log(0.153); label("Clearance at reference covariates (CL, L/day)")            # Budha 2023 Table 3, exp(theta1)*24 = 0.153 L/day
    lvc  <- log(3.05);  label("Central volume of distribution at reference (Vc, L)")     # Budha 2023 Table 3, exp(theta2) = 3.05 L
    lq   <- log(0.740); label("Intercompartmental clearance to peripheral1 (Q2, L/day)") # Budha 2023 Table 3, exp(theta3)*24 = 0.740 L/day
    lvp  <- log(1.27);  label("Peripheral1 volume of distribution (V2, L)")              # Budha 2023 Table 3, exp(theta4) = 1.27 L
    lq2  <- log(0.092); label("Intercompartmental clearance to peripheral2 (Q3, L/day)") # Budha 2023 Table 3, exp(theta5)*24 = 0.092 L/day
    lvp2 <- log(2.10);  label("Peripheral2 volume of distribution (V3, L)")              # Budha 2023 Table 3, exp(theta6) = 2.10 L

    # Covariate effects on CL (Budha 2023 Equation 5)
    e_wt_cl        <-  0.565;  label("Power exponent of WT on CL (unitless)")                            # Budha 2023 Table 3, theta7
    e_alb_cl       <- -0.457;  label("Power exponent of ALB on CL (unitless)")                           # Budha 2023 Table 3, theta10
    e_tumsz_cl     <-  0.0735; label("Power exponent of TUMSZ on CL (unitless)")                         # Budha 2023 Table 3, theta11
    e_ada_cl       <-  0.111;  label("Exponential coefficient of ADA-positive on CL (unitless)")         # Budha 2023 Table 3, theta13
    e_tumtp_gc_cl  <-  0.069;  label("Exponential coefficient of gastric cancer tumor type on CL (unitless)") # Budha 2023 Table 3, theta14
    e_tumtp_chl_cl <- -0.216;  label("Exponential coefficient of cHL tumor type on CL (unitless)")       # Budha 2023 Table 3, theta15

    # Covariate effects on Vc (Budha 2023 Equation 6)
    e_wt_vc   <-  0.397;  label("Power exponent of WT on Vc (unitless)")                          # Budha 2023 Table 3, theta8
    e_sexf_vc <- -0.116;  label("Exponential coefficient of female sex on Vc (unitless)")         # Budha 2023 Table 3, theta9
    e_age_vc  <-  0.0966; label("Power exponent of AGE on Vc (unitless)")                         # Budha 2023 Table 3, theta12

    # IIV: omega^2 = log(CV^2 + 1). CL-Vc block is correlated (covariance 0.020 on log-scale).
    # V2 and V3 IIV are retained per the paper despite high shrinkage (55.8% and 44.4%).
    etalcl + etalvc ~ c(0.066881,
                        0.020, 0.027506)  # CL 26.3% CV, Vc 16.7% CV, cov 0.020 — Budha 2023 Table 3
    etalvp  ~ 0.443417                    # V2 74.7% CV — Budha 2023 Table 3 (shrinkage 55.8%)
    etalvp2 ~ 0.692148                    # V3 99.9% CV — Budha 2023 Table 3 (shrinkage 44.4%)

    # Combined additive + proportional residual error (Budha 2023 Equation 2, Table 3)
    propSd <- 0.126; label("Proportional residual error (fraction)")  # Budha 2023 Table 3, sigma_p = 12.6%
    addSd  <- 2.09;  label("Additive residual error (ug/mL)")         # Budha 2023 Table 3, sigma_a = 2.09 ug/mL
  })
  model({
    # Individual PK parameters with covariate adjustments (Budha 2023 Equations 5 and 6).
    # Reference covariate values: WT 65 kg, AGE 60 y, ALB 41 g/L, TUMSZ 63 mm, male, ADA-negative,
    # non-cHL / non-GC tumor. Power effects act on log-transformed continuous covariates in NONMEM,
    # which translate to (COV / ref)^exponent in linear space. Categorical effects act exponentially.
    cl <- exp(lcl + etalcl) *
      (WT    / 65)^e_wt_cl *
      (ALB   / 41)^e_alb_cl *
      (TUMSZ / 63)^e_tumsz_cl *
      exp(e_ada_cl       * ADA_POS) *
      exp(e_tumtp_gc_cl  * TUMTP_GC) *
      exp(e_tumtp_chl_cl * TUMTP_CHL)

    vc <- exp(lvc + etalvc) *
      (WT  / 65)^e_wt_vc *
      (AGE / 60)^e_age_vc *
      exp(e_sexf_vc * SEXF)

    q   <- exp(lq)
    vp  <- exp(lvp + etalvp)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2 + etalvp2)

    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-                    k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-                                                      k13 * central - k31 * peripheral2

    # Observation: dose in mg, volume in L -> mg/L = ug/mL
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
