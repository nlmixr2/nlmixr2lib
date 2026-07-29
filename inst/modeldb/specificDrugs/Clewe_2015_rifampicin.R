Clewe_2015_rifampicin <- function() {
  description <- "Pharmacometric pulmonary distribution model for rifampicin in adults without tuberculosis: a one-compartment plasma PK model with single-transit oral absorption coupled to a Smythe 2012 enzyme-pool autoinduction structure (MTT, N, EMAX, EC50, kENZ all fixed from the upstream Smythe 2012 model) plus two effect compartments capturing distribution from plasma to epithelial lining fluid (ELF) and alveolar cells (AC); CL/F and Vc/F are FFM-allometrically scaled to 70 kg, the ELF and AC equilibration rate constants kELF and kAC are fixed to an equivalent 1-min half-life (instantaneous distribution at the single 4-h post-dose BAL sampling time), and only the unbound steady-state ELF/plasma and AC/plasma concentration ratios are estimated (1.28 and 5.5 after correction for the 20% rifampicin plasma free fraction)."
  reference <- paste(
    "Clewe O., Goutelle S., Conte J. E. Jr., Simonsson U. S. H. (2015).",
    "A pharmacometric pulmonary model predicting the extent and rate of",
    "distribution from plasma to epithelial lining fluid and alveolar",
    "cells -- using rifampicin as an example.",
    "European Journal of Clinical Pharmacology 71(3):313-319.",
    "doi:10.1007/s00228-014-1798-3.",
    "PK structure (one-compartment + single transit + enzyme-pool",
    "autoinduction) and the fixed autoinduction parameters (MTT, N, EMAX,",
    "EC50, kENZ) are inherited from Smythe et al. (2012) Antimicrob",
    "Agents Chemother 56(4):2091-2098 doi:10.1128/AAC.05792-11;",
    "see modellib('Svensson_2016_rifampicin') for the same PK backbone",
    "applied to a TB sputum dataset.",
    sep = " "
  )
  vignette <- "Clewe_2015_rifampicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass at baseline, time-fixed per subject.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the allometric scaling of CL/F (exponent 0.75) and V/F (exponent 1.0), both standardized to a 70-kg patient. The Clewe 2015 final model uses FFM (the model with one fewer parameter than the equivalent NFM scaling that gave the same OFV; Clewe 2015 Results paragraph 2). FFM is computed from WT, HT, and SEXF via the Janmahasatian / Anderson-Holford WHSMAX/WHS50 formula (Clewe 2015 Eq. 3) -- men: WHSMAX = 42.92 kg/m^2, WHS50 = 30.93 kg/m^2; women: WHSMAX = 37.99 kg/m^2, WHS50 = 35.98 kg/m^2; FFM = WHSMAX * HT^2 * WT / (WHS50 * HT^2 + WT). The vignette shows the FFM derivation chunk; user data sets supply FFM directly as a covariate column.",
      source_name        = "FFM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 40L,
    n_studies      = 1L,
    age_range      = "adult (Conte 2004 source trial; specific range not reported in Clewe 2015)",
    weight_range   = "not reported in Clewe 2015 (Conte 2004 source cohort)",
    sex_female_pct = 50,
    race_ethnicity = "not reported in Clewe 2015",
    disease_state  = "Adults without active tuberculosis: 10 women without AIDS, 10 men without AIDS, 10 women with AIDS, and 10 men with AIDS (Clewe 2015 Methods, 'Data' paragraph 1). HIV / AIDS status is descriptive of the cohort; it was not tested as a covariate on the PK parameters in Clewe 2015 (Methods paragraph 2: 'The influence of potential subpopulation-specific properties or covariates were not explored in this analysis').",
    dose_range     = "Rifampicin 600 mg orally once daily for 5 days; bronchoalveolar lavage performed approximately 4 h after the last dose on day 5.",
    notes          = "Original concentration data source: Conte et al. 2004 (Antimicrob Agents Chemother). Plasma sampled twice (approximately 2 h and 4 h post-dose) on day 5; BAL fluid sampled once per subject approximately 4 h after the last dose, with ELF volume estimated by urea dilution and AC volume estimated from BAL cell count. Total observations: 76 plasma, 32 ELF, 36 AC. The model also assumes ELF / AC protein binding is negligible (Clewe 2015 Methods paragraph 5)."
  )

  ini({
    # =========================================================================
    # Rifampicin plasma PK (Clewe 2015 Table 1, 'Parameter estimates of the
    # final model'). Time unit is HOURS throughout (matching the paper's
    # parameter table). Plasma volume V/F gives concentration units of mg/L
    # (= ug/mL) given a 600-mg oral dose.
    # =========================================================================
    lcl     <- log(3.85)                ; label("Apparent oral clearance CL/F (L/h) at the preinduced state, FFM-standardized to 70 kg")  # Clewe 2015 Table 1: TV(CL/F)STD = 3.85 L/h (95% CI 2.26-8.68; RSE 3.1%)
    lvc     <- log(76.6)                ; label("Apparent central volume of distribution V/F (L), FFM-standardized to 70 kg")              # Clewe 2015 Table 1: TV(Vc/F)STD = 76.6 L (95% CI 60.85-88.83; RSE 2.7%)

    # =========================================================================
    # Smythe 2012 autoinduction parameters -- all fixed in Clewe 2015 because
    # the 5-day dosing window contained insufficient information about the
    # ~40-day time to autoinduction steady-state (Clewe 2015 Methods paragraph
    # 2). Values match Smythe 2012 Model 3 (single-transit absorption + enzyme
    # turnover); identical numbers also appear in Svensson 2016 Table 2.
    # =========================================================================
    lmtt    <- fixed(log(0.71))         ; label("Mean transit time (h, on log scale)")                                                      # Clewe 2015 Table 1: MTT = 0.71 h FIX; carried from Smythe 2012
    nn_fix  <- fixed(1)                 ; label("Number of transit compartments (integer, unitless)")                                       # Clewe 2015 Table 1: N = 1 FIX; carried from Smythe 2012
    lemax   <- fixed(log(1.04))         ; label("Maximal fractional increase in enzyme production rate (unitless)")                         # Clewe 2015 Table 1: EMAX = 1.04 FIX; carried from Smythe 2012
    lec50   <- fixed(log(0.0705))       ; label("Rifampicin plasma concentration producing half EMAX (mg/L)")                                # Clewe 2015 Table 1: EC50 = 0.0705 mg/L FIX; carried from Smythe 2012
    lkenz   <- fixed(log(0.0036))       ; label("First-order degradation rate of the enzyme pool (1/h)")                                    # Clewe 2015 Table 1: kENZ = 0.0036/h FIX; carried from Smythe 2012
    lfdepot <- fixed(log(1))            ; label("Oral bioavailability (fixed at 1 because CL and V are apparent F-relative)")               # Implicit anchor: CL/F and V/F are reported in Clewe 2015 Table 1, so F is structurally 1

    # =========================================================================
    # Plasma-to-ELF and plasma-to-AC effect-compartment distribution rate
    # constants. The 4-h post-dose BAL sampling design gave only one sample
    # per subject, so neither rate could be estimated. Both kELF and kAC
    # were therefore FIXED to the value 60 * log(2) = 41.59/h (~1-minute
    # equilibration half-life, effectively instantaneous; Clewe 2015 Results
    # paragraph 2 'an equilibration half-life of about 1 min').
    # =========================================================================
    lkelf   <- fixed(log(41.58))        ; label("Plasma-to-ELF effect-compartment rate constant (1/h)")                                      # Clewe 2015 Table 1: kELF = 41.58/h FIX (equivalent to a 1-min equilibration half-life)
    lkac    <- fixed(log(41.58))        ; label("Plasma-to-AC effect-compartment rate constant (1/h)")                                       # Clewe 2015 Table 1: kAC = 41.58/h FIX (equivalent to a 1-min equilibration half-life)

    # =========================================================================
    # Pseudo-steady-state ELF/plasma and AC/plasma concentration ratios
    # (extent of distribution). Only the total-plasma ratios R_ELF/plasma and
    # R_AC/plasma are model parameters; the unbound ratios reported in
    # Clewe 2015 Table 1 (R_ELF/unbound-plasma = 1.28, R_AC/unbound-plasma =
    # 5.5) are derived post-estimation by dividing by the 20% rifampicin
    # plasma free fraction (Clewe 2015 ref [23]). ELF / AC protein binding
    # is assumed negligible.
    # =========================================================================
    lrelf   <- log(0.26)                ; label("ELF / plasma concentration ratio at pseudo steady-state (unitless)")                       # Clewe 2015 Table 1: R_ELF/plasma = 0.26 (95% CI 0.21-0.31; RSE 4.3%); R_ELF/unbound-plasma = 1.28 after dividing by fu = 0.2
    lrac    <- log(1.1)                 ; label("AC / plasma concentration ratio at pseudo steady-state (unitless)")                         # Clewe 2015 Table 1: R_AC/plasma = 1.1 (95% CI 0.92-1.35; RSE 6.2%); R_AC/unbound-plasma = 5.5 after dividing by fu = 0.2

    # =========================================================================
    # IIV. Clewe 2015 reports IIV only on CL/F (88.8% CV); no IIV was
    # quantifiable for V/F or for the distribution / extent parameters
    # because of the sparse one-sample-per-subject BAL design (Clewe 2015
    # Discussion paragraph 2: 'In this analysis example, no IIV was
    # quantified for the parameters describing the rate or extent of
    # distribution. In order to allow for IIV to be quantified with good
    # precision, more than one sample per subject is needed.').
    #
    # omega^2 = log(1 + CV^2) for a log-normal random effect:
    # omega^2 = log(1 + 0.888^2) = log(1.7885) = 0.5814.
    # =========================================================================
    etalcl  ~ 0.5814                    # Clewe 2015 Table 1: IIV_CL/F = 88.8% CV (95% CI 9.43-106.77; RSE 24.2%)

    # =========================================================================
    # Residual error. Plasma, ELF, and AC each carry a proportional
    # residual SD on the linear concentration scale (Clewe 2015 Table 1
    # 'Plasma / ELF / AC proportional error').
    # =========================================================================
    propSd      <- 0.352                ; label("Plasma proportional residual error (fraction)")  # Clewe 2015 Table 1: Plasma proportional error = 35.2% (95% CI 25.11-45.42; RSE 3.6%)
    propSd_Celf <- 0.407                ; label("ELF proportional residual error (fraction)")     # Clewe 2015 Table 1: ELF proportional error = 40.7% (95% CI 30.26-54.76; RSE 2.9%)
    propSd_Cac  <- 0.371                ; label("AC proportional residual error (fraction)")      # Clewe 2015 Table 1: AC proportional error = 37.1% (95% CI 22.95-46.91; RSE 7.3%)
  })

  model({
    # --- 1. FFM allometric scaling (Clewe 2015 Eqs. 1 and 2):
    #        CL/F = TV(CL/F)_STD * (FFM/70)^0.75
    #        V/F  = TV(V/F)_STD  * (FFM/70)^1.0
    #        Standardized to a 70-kg subject.
    cl_base <- exp(lcl + etalcl) * (FFM / 70) ^ 0.75
    vc      <- exp(lvc)          * (FFM / 70) ^ 1.0

    # --- 2. Smythe 2012 autoinduction parameters (all log-fixed in ini()).
    mtt   <- exp(lmtt)
    kenz  <- exp(lkenz)
    emax  <- exp(lemax)
    ec50  <- exp(lec50)
    kelf  <- exp(lkelf)
    kac   <- exp(lkac)
    relf  <- exp(lrelf)
    rac   <- exp(lrac)

    # --- 3. Transit-compartment rate constant. Savic-Karlsson convention:
    #        ktr = (N + 1) / MTT. With N = 1, ktr = 2 / 0.71 = 2.817/h.
    nn  <- nn_fix
    ktr <- (nn + 1) / mtt

    # --- 4. Plasma concentration (mg/L) and effective clearance. The
    #        enz_pool state starts at 1 by construction (Smythe 2012 idiom:
    #        zero-order enzyme production = kENZ, so enz = 1 at steady state
    #        with no drug present); when rifampicin is present, the enzyme
    #        pool grows above 1, multiplying the preinduced clearance.
    Cc <- central / vc
    cl <- cl_base * enz_pool

    # --- 5. PK ODE system: depot -> transit1 -> central with enzyme pool
    #        autoinduction. The first-order elimination from central is
    #        written as cl * Cc rather than kel * central so the
    #        autoinduction multiplier is read directly off cl (= cl_base
    #        * enz_pool); the two forms are algebraically identical
    #        (cl * Cc == cl * central / vc == kel * central).
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot - ktr * transit1
    d/dt(central)  <-  ktr * transit1 - cl * Cc
    d/dt(enz_pool) <-  kenz * (1 + emax * Cc / (ec50 + Cc)) - kenz * enz_pool

    # Enzyme pool starts at unity by construction (Smythe 2012 normalization).
    enz_pool(0) <- 1.0

    # --- 6. ELF and AC effect compartments (Clewe 2015 Methods 'ELF and AC
    #        drug penetration', Eqs. 4-5):
    #          dC_ELF/dt = kELF * (R_ELF/plasma * Cp - C_ELF)
    #          dC_AC /dt = kAC  * (R_AC /plasma * Cp - C_AC )
    #        The state variables effect1 and effect2 hold the ELF and AC
    #        concentrations directly (mg/L), so the observation outputs
    #        Celf and Cac below are bare aliases.
    d/dt(effect1) <- kelf * (relf * Cc - effect1)
    d/dt(effect2) <- kac  * (rac  * Cc - effect2)

    # --- 7. Bioavailability (anchor at 1; CL/F and V/F are apparent
    #        F-relative).
    f(depot) <- exp(lfdepot)

    # --- 8. Observations.
    Celf <- effect1
    Cac  <- effect2

    Cc   ~ prop(propSd)
    Celf ~ prop(propSd_Celf)
    Cac  ~ prop(propSd_Cac)
  })
}
