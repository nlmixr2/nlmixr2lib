Chae_2012_metformin <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral metformin in healthy Korean adults, coupled to a three-transit Sun-Jusko signal-transduction PD model for the antihyperglycaemic effect (Chae 2012). Plasma drug concentration in the central compartment drives a Hill-type stimulation function DR = Emax * Cp^r / (EC50^r + Cp^r) that initiates a cascade of three secondary-messenger transit compartments (M1 -> M2 -> M3) with shared mean transit time tau. The third messenger M3 is the measured percent change in plasma glucose from baseline relative to a sugar-bolus control arm. Creatinine clearance enters CL/F as a power covariate with reference 106.5 mL/min and exponent 0.782."
  reference <- "Chae JW, Baek IH, Lee BY, Cho SK, Kwon KI. Population pharmacokinetic and pharmacodynamic analysis of metformin using the signal transduction model. Br J Clin Pharmacol. 2012;74(5):815-823. doi:10.1111/j.1365-2125.2012.04260.x"
  vignette <- "Chae_2012_metformin"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance as reported by the source paper (Chae 2012 Table 1, mL/min). The paper does not state the derivation method (measured vs. Cockcroft-Gault).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on CL/F with reference 106.5 mL/min (cohort median per Chae 2012 Table 1) and exponent 0.782 (Chae 2012 final-PK structural equation). Raw mL/min as reported by the paper; not BSA-normalised. Cohort range 90-123 mL/min in 42 healthy young Korean males.",
      source_name        = "CLCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 42L,
    n_studies      = 1L,
    age_range      = "21-31 years",
    age_median     = "27 years",
    weight_range   = "61-78 kg",
    weight_median  = "69 kg",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Healthy",
    dose_range     = "500 mg single oral dose (Diabex 500 mg tablet, Daewoong)",
    regions        = "Korea (Daejeon)",
    notes          = "All-male healthy Korean cohort (Chae 2012 Table 1). 504 metformin plasma concentrations + 504 glucose concentrations (1008 total observations). All volunteers consumed 12 g of sugar 20 min after metformin dosing; an identical no-metformin sugar-bolus control study was performed 1 week later to derive the percent-change-in-glucose PD endpoint. Baseline fasting plasma glucose 98 +/- 7 mg/dL; haemoglobin 16 +/- 0.8 g/dL; total bilirubin 1.1 +/- 0.3 mg/dL; height 1.74 +/- 0.06 m."
  )

  ini({
    # Structural PK parameters (Chae 2012 Table 2; estimated point estimates with %RSE).
    lka  <- log(0.41);  label("Absorption rate constant Ka (1/h)")             # Chae 2012 Table 2 (RSE 2.43%)
    lcl  <- log(52.6);  label("Apparent clearance CL/F (L/h)")                  # Chae 2012 Table 2 (RSE 4.18%)
    lvc  <- log(113);   label("Apparent central volume of distribution V/F (L)") # Chae 2012 Table 2 (RSE 56.6%)

    # Covariate effect on CL/F (Chae 2012 final structural equation, p.819:
    # CL/F = 52.6 * (CLcr / 106.5)^0.782). Reference 106.5 mL/min = cohort median.
    e_crcl_cl <- 0.782; label("Power exponent of (CRCL/106.5) on CL/F (unitless)") # Chae 2012 p.819 final-PK equation

    # Structural PD parameters (Chae 2012 Table 2; signal-transduction model).
    ltau   <- log(0.50); label("Mean transit time of each PD transit compartment tau (h)") # Chae 2012 Table 2 (RSE 2.97%)
    lemax  <- log(19.8); label("Maximum stimulation Emax driving DR (percent change in glucose)") # Chae 2012 Table 2 (RSE 3.17%)
    lec50  <- log(3.68); label("EC50 driving DR (ug/mL)")                       # Chae 2012 Table 2 (RSE 3.89%)
    lhill  <- log(0.55); label("Hill coefficient r on Cp in DR (unitless)")     # Chae 2012 Table 2 (RSE 9.05%)

    # IIV (Chae 2012 Table 2 'Interindividual variability (IIV) CV%' column).
    # Variance on the log scale: omega^2 = log(CV^2 + 1).
    # CL: 29.7% CV -> omega^2 = log(0.297^2 + 1) = 0.08457
    # V : 22.1% CV -> omega^2 = log(0.221^2 + 1) = 0.04769
    # r : 4.05% CV -> omega^2 = log(0.0405^2 + 1) = 0.001639
    # Ka, tau, Emax, EC50 IIV reported as 'Not estimated' (Table 2 daggers) and
    # are not given etas in this model. The paper's descriptive equation
    # 'Ka = 0.41 * exp(eta1)' on p.819 is inconsistent with Table 2's daggered
    # Ka IIV; Table 2 (with paired bootstrap CIs) is treated as canonical here
    # and the discrepancy is flagged in the vignette Assumptions section.
    etalcl  ~ 0.08457   # Chae 2012 Table 2 (CL IIV 29.7% CV)
    etalvc  ~ 0.04769   # Chae 2012 Table 2 (V IIV 22.1% CV)
    etalhill ~ 0.001639 # Chae 2012 Table 2 (r IIV 4.05% CV)

    # Residual error.
    # PK: additive 23.0 ng/mL = 0.023 ug/mL on metformin Cc (Chae 2012 Table 2; RSE 11.7%).
    # PD: proportional 40.4% on the percent-effect output (Chae 2012 Table 2; RSE 1.52%).
    addSd          <- 0.023; label("Additive residual SD on metformin Cc (ug/mL)")        # Chae 2012 Table 2 PK residual (23 ng/mL)
    propSd_pctEffect <- 0.404; label("Proportional residual SD on PD percent-effect output") # Chae 2012 Table 2 PD residual (40.4% CV)
  })

  model({
    # 1. Individual parameters
    ka   <- exp(lka)
    cl   <- exp(lcl + etalcl) * (CRCL / 106.5)^e_crcl_cl
    vc   <- exp(lvc + etalvc)
    tau  <- exp(ltau)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill + etalhill)

    # 2. Micro-constants
    kel <- cl / vc

    # 3. PK ODE system (one-compartment, first-order absorption)
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    Cc <- central / vc

    # 4. Signal-transduction PD chain. DR is the Hill stimulation function (algebraic);
    # M1 -> M2 -> M3 form a transit chain with shared mean transit time tau and the
    # third transit compartment M3 is the observed percent-effect output. M1/M2/M3
    # map to canonical transit1/transit2/transit3 compartment names per
    # compartment-names.md (the paper writes the chain as M1, M2, M3).
    DR <- (emax * Cc^hill) / (ec50^hill + Cc^hill)
    d/dt(transit1) <- (DR       - transit1) / tau
    d/dt(transit2) <- (transit1 - transit2) / tau
    d/dt(transit3) <- (transit2 - transit3) / tau
    pctEffect <- transit3

    # 5. Residual error
    Cc        ~ add(addSd)
    pctEffect ~ prop(propSd_pctEffect)
  })
}
