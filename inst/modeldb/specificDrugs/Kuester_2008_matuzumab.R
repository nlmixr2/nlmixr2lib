Kuester_2008_matuzumab <- function() {
  description <- "Two-compartment population PK model for matuzumab (humanised anti-EGFR IgG1 monoclonal antibody) in adults with advanced carcinoma (Kuester 2008), with parallel first-order linear and Michaelis-Menten elimination from the central compartment; body weight on linear CL and central volume."
  reference <- "Kuester K, Kovar A, Lupfert C, Brockhaus B, Kloft C. Population pharmacokinetic data analysis of three phase I studies of matuzumab, a humanised anti-EGFR monoclonal antibody in clinical cancer development. Br J Cancer. 2008;98(5):900-906. doi:10.1038/sj.bjc.6604265"
  vignette <- "Kuester_2008_matuzumab"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear additive-offset effect on linear CL (CLL) and central volume (V1) per Kuester 2008 footnotes c and d (Table 3, p902): V1,individual = V1 * (1 + V1_WT * (WT - WTmedian)) * exp(eta_V1) and CLL,individual = CLL * (1 + CLL_WT * (WT - WTmedian)) * exp(eta_CLL + kappa_CLL). WTmedian = 71 kg (pooled cohort median, Kuester 2008 Table 2, p902). The paper does not state whether weight is time-varying; for a typical 90-patient phase I cohort with treatment durations up to ~1 year, baseline-fixed weight is the standard assumption.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 90L,
    n_observations = 1256L,
    n_studies      = 3L,
    age_range      = "29-82 years",
    age_median     = "60 years",
    weight_range   = "44-125 kg",
    weight_median  = "71 kg",
    sex_female_pct = 41.1,
    disease_state  = "Patients with advanced, non-resectable and/or metastatic histologically proven carcinoma (predominantly colorectal carcinoma in studies 2 and 3; advanced pancreatic cancer in study 1). Karnofsky performance status >= 60%, life expectancy >= 12 weeks (study 1: >= 8 weeks), no prior chemotherapy or radiotherapy within 3-4 weeks of study start, adequate renal, haematological, and hepatic function, age >= 18 years.",
    dose_range     = "Matuzumab 400-2000 mg by 1-hour intravenous infusion in 11 different regimens spanning q1w, q2w, and q3w schedules. Study 1: 400 mg q1w; 800 mg q2w; 800 mg q1w in combination with fixed-dose gemcitabine 1000 mg/m^2. Study 2: 1200 mg q1w/q2w/q3w; 400, 800, 1600 mg q3w. Study 3: 400, 800, 1200, 1600 mg q1w; 2000 mg q1w (reduced to 1600 mg from week 2). Some patients treated for ~1 year.",
    regions        = "Three phase I, open-label, non-randomised, uncontrolled, multicentre studies (4 study sites).",
    notes          = paste(
      "Baseline demographics from Kuester 2008 Table 2 (p902). Pooled cohort: 53 male / 37 female across 3 studies.",
      "Body surface area median 1.82 m^2 (range 1.34-2.59); body mass index median 24.9 kg/m^2 (range 15.9-37.0).",
      "Creatinine clearance median 91 mL/min (range 41-480); alkaline phosphatase median 171 U/L; lactate dehydrogenase median 426 U/L (Table 2).",
      "Concomitant chemotherapy (gemcitabine) administered to all 17 study-1 subjects only.",
      "Bioanalytical method: validated sandwich ELISA (Vanhoefer et al 2004). Estimation: NONMEM V level 1.1; ADVAN6 TRANS1 TOL5; FOCE with interaction (Methods, p901).",
      sep = " "
    ),
    iov_note = "Inter-occasion variability on CLL (22.8% CV) was estimated with every infusion treated as one occasion (Kuester 2008 Table 3 'Interoccasion variability' row, p902). IOV is not encoded structurally in this model file -- there is no standardized OCC indicator in the nlmixr2lib event-data schema. Downstream users who want to simulate IOV can add an OCC indicator column to the event dataset and a per-occasion eta in rxode2. See the vignette 'Assumptions and deviations' section."
  )

  ini({
    # Structural PK parameters - Kuester 2008 Table 3 final-model estimates (p902).
    # Reference covariate value for the typical subject: WT = 71 kg (pooled median).
    # Time unit is hour; CLL and Q are reported in mL/h (= 0.001 L/h); V1, V2 in L;
    # Vmax in mg/h; Km in mg/L (= ug/mL since 1 mg/L = 1 ug/mL).
    lcl   <- log(14.5 / 1000); label("Linear clearance CLL (L/h)")                                   # Kuester 2008 Table 3 CLL = 14.5 mL/h (RSE 4.1%)
    lvc   <- log(3.72);        label("Central volume of distribution V1 (L)")                        # Kuester 2008 Table 3 V1 = 3.72 L (RSE 3.0%)
    lq    <- log(38.3 / 1000); label("Inter-compartmental clearance Q (L/h)")                        # Kuester 2008 Table 3 Q = 38.3 mL/h (RSE 7.6%)
    lvp   <- log(1.84);        label("Peripheral volume of distribution V2 (L)")                     # Kuester 2008 Table 3 V2 = 1.84 L (RSE 9.0%)
    lvmax <- log(0.456);       label("Maximum Michaelis-Menten elimination rate Vmax (mg/h)")        # Kuester 2008 Table 3 Vmax = 0.456 mg/h (RSE 13.7%)
    lkm   <- log(4.0);         label("Michaelis-Menten constant Km (ug/mL)")                         # Kuester 2008 Table 3 Km = 4.0 mg/L = 4.0 ug/mL (RSE 29.8%)

    # Covariate effects - Kuester 2008 Table 3 final-model "Covariate influence" rows
    # (p902) parameterized per the footnotes c and d as:
    #   V1,individual  = V1  * (1 + V1_WT  * (WT - 71)) * exp(eta_V1)
    #   CLL,individual = CLL * (1 + CLL_WT * (WT - 71)) * exp(eta_CLL + kappa_CLL)
    # Reported coefficients: V1_WT = +0.44% per kg, CLL_WT = +0.87% per kg.
    e_wt_vc <- 0.0044; label("Linear coefficient of (WT - 71 kg) on V1 (1/kg)")                       # Kuester 2008 Table 3 V1_WT (RSE 35.2%)
    e_wt_cl <- 0.0087; label("Linear coefficient of (WT - 71 kg) on CLL (1/kg)")                      # Kuester 2008 Table 3 CLL_WT (RSE 28.2%)

    # Inter-individual variability - Kuester 2008 Table 3 ("Random effects" section, p902).
    # The paper reports CV% on the linear-parameter scale. Convert to log-normal variance:
    #   omega^2 = log(1 + CV^2)
    # IIV on CLL is independent (no reported correlations involving CLL). IIV on V1, V2,
    # and Vmax form a 3x3 correlated block with the reported correlation coefficients:
    #   r(V1, V2)  = 0.777
    #   r(V1, Vmax) = 0.875
    #   r(V2, Vmax) = 0.875
    #
    # Variances:
    #   CLL  CV 24.0% -> omega^2 = log(1 + 0.240^2) = 0.05598, omega = 0.2366
    #   V1   CV 21.9% -> omega^2 = log(1 + 0.219^2) = 0.04685, omega = 0.2164
    #   V2   CV 61.6% -> omega^2 = log(1 + 0.616^2) = 0.32149, omega = 0.5670
    #   Vmax CV 53.8% -> omega^2 = log(1 + 0.538^2) = 0.25416, omega = 0.5041
    #
    # Off-diagonal covariances cov_ij = r_ij * omega_i * omega_j:
    #   cov(V1, V2)   = 0.777 * 0.2164 * 0.5670 = 0.0954
    #   cov(V1, Vmax) = 0.875 * 0.2164 * 0.5041 = 0.0955
    #   cov(V2, Vmax) = 0.875 * 0.5670 * 0.5041 = 0.2501
    etalcl ~ 0.05598                                                                                  # Kuester 2008 Table 3 IIV CLL = 24.0% CV (RSE 20.5%)
    etalvc + etalvp + etalvmax ~ c(0.04685,
                                   0.0954, 0.32149,
                                   0.0955, 0.2501, 0.25416)                                           # Kuester 2008 Table 3 IIV V1/V2/Vmax + correlations (V1-V2 = 0.777, V1-Vmax = 0.875, V2-Vmax = 0.875)

    # Residual error - Kuester 2008 Table 3 ("Random effects" section, p902). Combined
    # proportional + additive: Cobs = Cpred * (1 + eps_prop) + eps_add. The additive
    # term was held fixed at 0.312 mg/L (= 0.312 ug/mL) "due to model stability (the
    # value was chosen from prior plausible successfully run models)" (Kuester 2008
    # Results, p904). 1 mg/L = 1 ug/mL so the units match this model's declaration.
    propSd <- 0.134;        label("Proportional residual error (fraction)")                          # Kuester 2008 Table 3 proportional = 13.4% CV (RSE 1.5%)
    addSd  <- fixed(0.312); label("Additive residual error (ug/mL)")                                 # Kuester 2008 Table 3 additive = 0.312 mg/L = 0.312 ug/mL (fixed)
  })

  model({
    # Individual PK parameters with the Kuester 2008 final-model covariate equations.
    # Reference subject: WT = 71 kg (pooled cohort median, Table 2). The linear
    # additive-offset form (1 + e_wt_param * (WT - 71)) is taken directly from the
    # Table 3 footnotes c and d (p902).
    cl   <- exp(lcl + etalcl) * (1 + e_wt_cl * (WT - 71))
    vc   <- exp(lvc + etalvc) * (1 + e_wt_vc * (WT - 71))
    vp   <- exp(lvp + etalvp)
    q    <- exp(lq)
    vmax <- exp(lvmax + etalvmax)
    km   <- exp(lkm)

    # Two-compartment IV PK with parallel linear and Michaelis-Menten elimination
    # from the central compartment (Kuester 2008 Figure 2, p902). Dose enters the
    # central compartment directly via 1-hour IV infusion; no depot compartment.
    # Concentration in ug/mL; central mass in mg, V1 in L (1 mg/L = 1 ug/mL).
    Cc <- central / vc

    d/dt(central)     <- -(cl / vc) * central -
                          vmax * Cc / (km + Cc) -
                          (q / vc) * central +
                          (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central -
                          (q / vp) * peripheral1

    Cc ~ add(addSd) + prop(propSd)
  })
}
