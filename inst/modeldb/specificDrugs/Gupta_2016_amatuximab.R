Gupta_2016_amatuximab <- function() {
  description <- "Two-compartment population PK model with parallel linear and Michaelis-Menten elimination for amatuximab in patients with advanced cancers / malignant pleural mesothelioma (Gupta 2016)"
  reference <- "Gupta A, Hussein Z, Hassan R, Wustner J, Maltzman JD, Wallin BA. Population pharmacokinetics and exposure-response relationship of amatuximab, an anti-mesothelin monoclonal antibody, in patients with malignant pleural mesothelioma and its application in dose selection. Cancer Chemother Pharmacol. 2016;77(4):733-743. doi:10.1007/s00280-016-2984-z"
  vignette <- "Gupta_2016_amatuximab"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on central volume (Vc); reference weight 70 kg from Gupta 2016 Table 2 covariate equation Vc = ThetaVc * (WGT/70)^ThetaWGT.",
      source_name        = "WGT"
    ),
    ADA_POS = list(
      description        = "Antidrug-antibody (ADA) indicator with study-specific threshold: 1 = ADA titer > 64, 0 = otherwise (ADA-negative or titer <= 64)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0L,
      notes              = "Multiplicative effect on linear CL (CL = ThetaCL * ThetaADA^ADA). Gupta 2016 categorized ADA positivity at five candidate titer thresholds (>1, >4, >64, >160, any positive); >64 was the retained form in the final model (Gupta 2016 p. 738 and Table 2). Time-varying in principle (titer changes during treatment) but entered as a single covariate value per record in NONMEM.",
      source_name        = "ADA"
    )
  )

  population <- list(
    n_subjects     = 199L,
    n_studies      = 4L,
    age_range      = "33-90 years",
    age_median     = "65 years",
    weight_range   = "35-134 kg",
    weight_median  = "74 kg",
    sex_female_pct = 35.7,
    race_ethnicity = "Caucasian 81.4% (162/199); Japanese 8.54% (17/199); Other 8.04% (16/199); missing 2.01% (4/199)",
    disease_state  = "Advanced mesothelin-expressing cancers (malignant pleural mesothelioma [MPM], pancreatic cancer, other solid tumors) pooled across four clinical studies.",
    dose_range     = "Amatuximab IV infusion: Phase I dose-escalation 12.5-400 mg/m^2 weekly (US) and 50-200 mg/m^2 weekly (Japan); Phase II pancreatic cancer 5 mg/kg weekly; Phase II MPM 5 mg/kg on Days 1 and 8 of each 21-day cycle (plus pemetrexed + cisplatin).",
    regions        = "North America, Japan, and international Phase II sites",
    baseline_albumin = "3.80 +/- 0.53 g/dL (median 3.80, range 2.38-5.46)",
    ecog_status    = "ECOG 0: 59.8%; ECOG 1: 39.2%; ECOG 2: 1%",
    notes          = paste(
      "Pooled analysis across four studies: US Phase I (NCT00325494, 24 patients),",
      "Japanese Phase I (NCT01018784, 17 patients), Phase II pancreatic cancer",
      "(NCT00570713, 71 patients contributing PK data), and Phase II MPM",
      "(NCT00738582, 87 patients). Baseline demographics per Gupta 2016 Table 1",
      "pharmacokinetic-analysis database column (N = 199). 3236 amatuximab serum",
      "concentrations (1.51% BLQ)."
    )
  )

  ini({
    # Structural parameters (Gupta 2016 Table 2, final model estimates).
    # Units: CL and Q in L/h, Vc and Vp in L, Vmax in mg/h, Km in mg/L
    # (converted from paper's 790 ng/mL so that mg dosing and L volumes give
    # concentrations consistent in mg/L = ug/mL).
    lcl   <- log(0.0299); label("Linear clearance CL at reference 70 kg, ADA-negative (L/h)")       # Gupta 2016 Table 2
    lvc   <- log(3.89);   label("Central volume of distribution Vc at reference 70 kg (L)")         # Gupta 2016 Table 2
    lq    <- log(0.0147); label("Intercompartmental clearance Q (L/h)")                             # Gupta 2016 Table 2
    lvp   <- log(2.62);   label("Peripheral volume of distribution Vp (L)")                         # Gupta 2016 Table 2
    lvmax <- log(0.173);  label("Maximum rate of saturable (target-mediated) elimination Vmax (mg/h)")  # Gupta 2016 Table 2
    lkm   <- log(0.790);  label("Michaelis-Menten constant Km (mg/L = ug/mL)")                      # Gupta 2016 Table 2 (790 ng/mL converted)

    # Covariate effects
    e_wt_vc  <- 0.597; label("Power exponent for body weight on Vc (unitless)")                     # Gupta 2016 Table 2 (ThetaWGT)
    e_ada_cl <- 1.49;  label("Multiplicative effect on CL for ADA titer > 64 (unitless factor)")    # Gupta 2016 Table 2 (ThetaADA)

    # Inter-individual variability. Paper reports CV% (Table 2 column header
    # "Inter-individual variability (CV %)"), converted to log-normal variance
    # via omega^2 = log(CV^2 + 1). CL and Vc are correlated (r = 0.411),
    # coded as a 2x2 block. Paper used a Manly transformation on eta-Vmax
    # (shape factor -0.181, 95% CI -0.450 to 0.0875, i.e., not statistically
    # distinguishable from log-normal); this implementation uses the standard
    # log-normal form — see the "Assumptions and deviations" section of the
    # validation vignette.
    etalcl + etalvc ~ c(0.056443,
                        0.023198, 0.056443)  # CV_CL = CV_Vc = 24.1%, r = 41.1% (Gupta 2016 Table 2)
    etalvmax ~ 0.931485                      # CV_Vmax = 124% (Gupta 2016 Table 2)

    # Residual variability (Gupta 2016 Table 2). Additive part fixed to 1/4 of
    # the LLOQ (98 ng/mL / 4 = 24.5 ng/mL = 0.0245 ug/mL).
    propSd <- 0.339;         label("Proportional residual error (CV fraction)")  # Gupta 2016 Table 2
    addSd  <- fixed(0.0245); label("Additive residual error (ug/mL)")            # Gupta 2016 Table 2 (24.5 ng/mL FIXED = 0.0245 ug/mL)
  })

  model({
    # Individual PK parameters
    cl   <- exp(lcl + etalcl) * e_ada_cl^ADA_POS          # CL = ThetaCL * ThetaADA^ADA (Gupta 2016 Table 2)
    vc   <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc         # Vc = ThetaVc * (WGT/70)^ThetaWGT (Gupta 2016 Table 2)
    q    <- exp(lq)
    vp   <- exp(lvp)
    vmax <- exp(lvmax + etalvmax)
    km   <- exp(lkm)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Concentration in central compartment (ug/mL = mg/L)
    Cc <- central / vc

    # Two-compartment IV model with parallel linear and Michaelis-Menten
    # elimination from central (Gupta 2016 Fig 1 / "Final PK model" section).
    d/dt(central)     <- -kel * central - vmax * Cc / (km + Cc) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                           k12 * central - k21 * peripheral1

    Cc ~ add(addSd) + prop(propSd)
  })
}
