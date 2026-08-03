Feillet_2008_sapropterin <- function() {
  description <- "Two-compartment population PK model with first-order oral absorption, an absorption lag, linear elimination, and an additive endogenous BH4 baseline for sapropterin dihydrochloride in adolescent and adult patients with BH4-responsive phenylketonuria (Feillet 2008)."
  reference <- "Feillet F, Clarke L, Meli C, et al. Pharmacokinetics of sapropterin in patients with phenylketonuria. Clinical Pharmacokinetics. 2008;47(12):817-825. doi:10.2165/0003088-200847120-00006"
  vignette <- "Feillet_2008_sapropterin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "sapropterin", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "sapropterin", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "sapropterin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline; adolescent and adult PKU cohort).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on apparent clearance and apparent central volume of distribution, normalized to a reference weight of 70 kg. Source: Feillet 2008 Equation 3 and Table III (CL/F reported as L/h/70 kg, V1/F as L/70 kg).",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 76L,
    n_studies      = 1L,
    age_range      = "9-50 years",
    age_mean       = "21.1 years (SD 9.64)",
    weight_range   = "28.2-144 kg",
    weight_mean    = "67.2 kg (SD 21.8)",
    height_range   = "126-191 cm",
    height_mean    = "165 cm (SD 13.3)",
    sex_female_pct = 42.3,
    race_ethnicity = c(White = 97, non_White = 3),
    disease_state  = "BH4-responsive phenylketonuria (PKU); 12-week fixed-dose phase of an open-label phase III extension study following the PKU-003 placebo-controlled trial.",
    dose_range     = "Oral once-daily sapropterin dihydrochloride 5, 10, or 20 mg/kg/day (100 mg tablets dissolved in water or juice; doses rounded up to the nearest 100 mg unit).",
    regions        = "26 centres in North America (Canada, USA) and Europe (France, Germany, Ireland, Italy, Poland, UK).",
    baseline_covariates = list(
      bsa_range_m2          = "1.05-2.65",
      bsa_mean_m2           = 1.72,
      alt_range_U_L         = "11-127",
      alt_mean_U_L          = 28.4,
      ast_range_U_L         = "14-43",
      ast_mean_U_L          = 25.7,
      bilirubin_range_mg_dL = "0.1-1.9",
      bilirubin_mean_mg_dL  = 0.55,
      scr_range_mg_dL       = "0.6-1.3",
      scr_mean_mg_dL        = 0.89,
      crcl_range_mL_min     = "48-231",
      crcl_mean_mL_min      = 114,
      baseline_phe_umol_L   = "53-2190 (mean 811, SD 393)"
    ),
    notes          = "Baseline demographics from Feillet 2008 Table II reflect all 78 enrolled patients (45 males, 33 females); the final model dataset used 265 plasma BH4 observations from 76 patients (315 observations - 38 below LLOQ - 4 unreported - 8 excluded as unreliable = 265). D-optimal sparse sampling (four samples per patient) obtained during weeks 6, 10, and 12 of the fixed-dose phase. BH4 measured indirectly via oxidation to L-biopterin and reversed-phase HPLC with tandem MS detection; nominal BH4 to L-biopterin conversion ratio was 47.3%. Bodyweight was the only covariate retained at p<0.001; sex, race, age, height, body surface area, serum creatinine, albumin, ALT, AST, total bilirubin and phenylalanine were screened but not retained."
  )

  ini({
    # Structural parameters -- Feillet 2008 Table III, final covariate model.
    # Reference weight for CL/F and V1/F is 70 kg (units L/h/70 kg and L/70 kg
    # per Table III header); V2/F and Q/F are population-wide (not weight-scaled).
    ltlag  <- log(0.275); label("Absorption lag time (tlag, h)")                              # Feillet 2008 Table III: tlag = 0.275 h
    lka    <- log(0.518); label("First-order absorption rate (ka, 1/h)")                      # Feillet 2008 Table III: ka = 0.518 h^-1
    lcl    <- log(2100);  label("Apparent clearance at WT = 70 kg (CL/F, L/h)")               # Feillet 2008 Table III: CL/F = 2100 L/h/70 kg
    lvc    <- log(8350);  label("Apparent central volume at WT = 70 kg (V1/F, L)")            # Feillet 2008 Table III: V1/F = 8350 L/70 kg
    lvp    <- log(4240);  label("Apparent peripheral volume (V2/F, L)")                       # Feillet 2008 Table III: V2/F = 4240 L (not weight-scaled)
    lq     <- log(862);   label("Apparent intercompartmental clearance (Q/F, L/h)")           # Feillet 2008 Table III: Q/F = 862 L/h (not weight-scaled)
    lrbase <- log(13.5);  label("Endogenous BH4 baseline plasma concentration (BASE, ng/mL)") # Feillet 2008 Table III: BASE = 13.5 ng/mL

    # Body-weight covariate effects (Equation 3, power form normalized to 70 kg).
    e_wt_cl <- 0.586; label("Power exponent of body weight on CL/F (unitless, reference 70 kg)") # Feillet 2008 Table III: power function on CL/F = 0.586
    e_wt_vc <- 1.13;  label("Power exponent of body weight on V1/F (unitless, reference 70 kg)") # Feillet 2008 Table III: power function on V1/F = 1.13

    # Inter-individual variability -- correlated CL/F and V1/F on the log scale.
    # Feillet 2008 Table III reports the SDs directly (log-normal ETA parameterisation,
    # Equation 1: P_j = TVP * exp(eta_j)), so omega_CL = 0.539, omega_V1 = 0.557 and
    # the CL-V1 correlation R = 0.336.
    #   omega^2(CL/F) = 0.539^2                     = 0.290521
    #   omega^2(V1/F) = 0.557^2                     = 0.310249
    #   cov(CL,V1)    = 0.336 * 0.539 * 0.557       = 0.100869
    etalcl + etalvc ~ c(0.290521, 0.100869, 0.310249)  # Feillet 2008 Table III: SD_CL = 0.539, SD_V1 = 0.557, R = 0.336

    # Residual error -- constant CV proportional error (Table III: 21.7% CV).
    propSd <- 0.217; label("Proportional residual error (constant CV, fraction)") # Feillet 2008 Table III: constant CV residual = 21.7%
  })

  model({
    # 1. Individual PK parameters (log-normal ETA on CL/F and V1/F only;
    #    Feillet 2008 did not include ETAs on ka, tlag, Q/F, V2/F, or BASE).
    ka    <- exp(lka)
    cl    <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc    <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp    <- exp(lvp)
    q     <- exp(lq)
    tlag  <- exp(ltlag)
    rbase <- exp(lrbase)

    # 2. Micro-constants for the two-compartment model.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 3. ODE system: oral dose into depot, two-compartment disposition.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 4. Absorption lag on the depot compartment (Feillet 2008 Fig. 1).
    alag(depot) <- tlag

    # 5. Observation: predicted plasma BH4 = drug-derived + endogenous baseline.
    #    Dose (amt) is in mg and vc is in L, so central/vc is in mg/L; multiply
    #    by 1000 to express Cc in ng/mL, matching the paper's reporting unit and
    #    the ng/mL scale of the estimated BASE parameter (Table III).
    Cc <- (central / vc) * 1000 + rbase
    Cc ~ prop(propSd)
  })
}
