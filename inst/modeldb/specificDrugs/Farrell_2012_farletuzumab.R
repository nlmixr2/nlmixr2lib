Farrell_2012_farletuzumab <- function() {
  description <- "Two-compartment population PK model for farletuzumab (humanized IgG1 anti-folate-receptor-alpha monoclonal antibody) with first-order linear elimination after IV infusion in women with advanced epithelial ovarian cancer (Farrell 2012)."
  reference <- "Farrell C, Schweizer C, Wustner J, Weil S, Namiki M, Nakano T, Nakai K, Phillips MD. Population pharmacokinetics of farletuzumab, a humanized monoclonal antibody against folate receptor alpha, in epithelial ovarian cancer. Cancer Chemother Pharmacol. 2012;70(5):727-734. doi:10.1007/s00280-012-1959-y"
  vignette <- "Farrell_2012_farletuzumab"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline) body weight used with power scaling on CL and Vc. Reference weight 66.2 kg is the overall median body weight in the analysis cohort (Farrell 2012 Table 1).",
      source_name        = "WT"
    ),
    PHASE2 = list(
      description        = "Indicator for the Phase II study (MORAb-003-002) of the Farrell 2012 pooled analysis",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Phase I study, MORAb-003-001)",
      notes              = "Switches the residual-error model between the Phase I proportional-only and the Phase II combined additive + proportional error reported in Farrell 2012 Table 3. Derived from the trial identifier.",
      source_name        = "STUDY"
    )
  )

  population <- list(
    n_subjects     = 79,
    n_studies      = 2,
    age_range      = "31-81 years",
    age_median     = "61 years",
    weight_range   = "44.5-118.2 kg",
    weight_median  = "66.2 kg",
    sex_female_pct = 100,
    race_ethnicity = c(Caucasian = 77, Asian = 8, Hispanic = 6, Other = 6, Black = 3),
    disease_state  = "Advanced epithelial ovarian, fallopian tube, or primary peritoneal cancer (Phase I relapsed after standard chemotherapy; Phase II relapsed platinum-sensitive).",
    dose_range     = "12.5-400 mg/m^2 weekly IV infusion",
    regions        = "Not reported in detail",
    notes          = "Demographics from Farrell 2012 Table 1. Two studies pooled (Phase I MORAb-003-001: 25 patients, 722 samples at doses 12.5-400 mg/m^2; Phase II MORAb-003-002: 54 patients, 1,750 samples predominantly at 100 mg/m^2). Final model development excluded 169 samples at 12.5 and 25 mg/m^2 doses due to apparent low-dose nonlinearity; all 79 patients contribute baseline demographics."
  )

  ini({
    # Structural parameters (final model) — Farrell 2012 Table 3.
    # Reference weight 66.2 kg = overall population median (Table 1; paper states
    # COV_ST is the median value of the covariate in the study population).
    lcl      <- log(0.00784); label("Clearance at reference body weight (CL, L/hr)")                     # Table 3: CL = 0.00784 L/h
    lvc      <- log(3.00);    label("Central volume of distribution at reference body weight (Vc, L)")  # Table 3: Vc = 3.00 L
    lq       <- log(0.0203);  label("Inter-compartmental clearance (Q, L/hr)")                          # Table 3: Q = 0.0203 L/h
    lvp      <- log(7.50);    label("Peripheral volume of distribution (Vp, L)")                        # Table 3: Vp = 7.50 L

    # Weight effect (power form on CL and Vc only; no weight effect on Q or Vp)
    e_wt_cl  <- 0.715;        label("Power exponent of body weight on CL (unitless)")                   # Table 3: CL ~ WT = 0.715
    e_wt_vc  <- 0.629;        label("Power exponent of body weight on Vc (unitless)")                   # Table 3: Vc ~ WT = 0.629

    # IIV — Farrell 2012 Table 3 reports diagonal omega^2 (variance on the log / exponential-error scale).
    # Narrative section "Population pharmacokinetic analysis" describes starting the analysis with a
    # CL-Vc correlation; Table 3 reports only diagonal omega^2 with no covariance term, so the final
    # model is implemented with diagonal IIV.
    etalcl   ~ 0.0616   # Table 3: omega^2 CL  = 0.0616 (24.8% CV)
    etalvc   ~ 0.0470   # Table 3: omega^2 Vc  = 0.0470 (21.7% CV)
    etalvp   ~ 1.180    # Table 3: omega^2 Vp  = 1.180  (109% CV)

    # Residual error — study-specific estimates (Farrell 2012 Table 3).
    # Phase I: proportional-only (sigma^2 prop = 0.0420 -> SD = 20.5%; sigma^2 add fixed = 0).
    # Phase II: combined additive + proportional (sigma^2 prop = 0.122 -> SD = 34.9%;
    # sigma^2 add = 63.0 ug^2/mL^2 -> SD = 7.94 ug/mL).
    propSdPh1 <- 0.205; label("Phase I proportional residual SD (fraction)")                            # Table 3: sqrt(0.0420) = 0.205
    propSdPh2 <- 0.349; label("Phase II proportional residual SD (fraction)")                           # Table 3: sqrt(0.122)  = 0.349
    addSdPh2  <- 7.94;  label("Phase II additive residual SD (ug/mL)")                                  # Table 3: sqrt(63.0)   = 7.94
  })

  model({
    # Individual PK parameters with power-weight scaling (reference 66.2 kg, Farrell 2012 Eq. for TVP_i)
    cl <- exp(lcl + etalcl) * (WT / 66.2)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 66.2)^e_wt_vc
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

    # Micro-constants for the two-compartment IV model
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # IV infusion into the central compartment; rate / duration supplied via the event table
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL
    Cc <- central / vc

    # Study-dependent residual error: PHASE2 = 1 selects Phase II (combined add + prop);
    # PHASE2 = 0 selects Phase I (proportional only; addSd collapses to 0).
    propSd <- propSdPh1 * (1 - PHASE2) + propSdPh2 * PHASE2
    addSd  <- addSdPh2 * PHASE2
    Cc ~ add(addSd) + prop(propSd)
  })
}
