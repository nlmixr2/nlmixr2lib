Rambiritch_2016_glibenclamide <- function() {
  description <- "Two-compartment population PK model with first-order oral absorption for glibenclamide in poorly controlled South African adults with type 2 diabetes (Rambiritch 2016). All disposition parameters are apparent (CL/F, Vc/F, Vp/F, Q/F); F is not estimated. Concentration data were log-transformed prior to NONMEM fitting (LTBS), giving an effectively proportional residual error in linear space. No covariate effects were retained in the final model."
  reference <- "Rambiritch V, Naidoo P, Maharaj B, Pillai G. Population pharmacokinetic modeling of glibenclamide in poorly controlled South African type 2 diabetic subjects. Clin Pharmacol (Auckl). 2016 Jul 12;8:83-92. doi:10.2147/CPAA.S102676"
  vignette <- "Rambiritch_2016_glibenclamide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 24,
    n_studies      = 1,
    age_range      = "39-73 years (mean 54, SD 9)",
    age_median     = NULL,
    weight_range   = "42.0-107.8 kg (mean 71.1, SD 14.1)",
    weight_median  = NULL,
    sex_female_pct = 91,
    race_ethnicity = "South African (race not stratified in source); cohort recruited at the University of KwaZulu-Natal / RK Khan Regional Hospital, Chatsworth, South Africa.",
    disease_state  = "Poorly controlled type 2 diabetes mellitus (mean fasting blood glucose 15.4 mmol/L; mean HbA1c 13.9%); requiring oral antidiabetic therapy.",
    dose_range     = "Glibenclamide 0, 2.5, 5, 10, and 20 mg orally once daily; dose levels escalated at 2-week intervals; sampling performed after steady state at each dose level.",
    regions        = "South Africa (Durban / Chatsworth)",
    sampling_design = "Per dose level (days 14, 28, 42, 56, and 70 for 0, 2.5, 5, 10, and 20 mg respectively): post-breakfast samples at 0, 30, 60, 90, 120 min; post-lunch samples at 240, 270, 300, 330, 360, and 420 min. 841 observation records from 24 individuals (22 with complete dose-escalation data; subjects 14, 16, 20, and 24 had partial profiles per Methods).",
    notes          = "Baseline demographics in Table 1 of Rambiritch 2016. Twenty (20) of 24 enrolled subjects were female. Mean creatinine 64.1 umol/L (within normal range); cohort assumed to have normal renal function. No functional covariate effects (weight, age, sex, renal markers) entered the final model."
  )

  ini({
    # Structural PK parameters (apparent, /F) -- Rambiritch 2016 Table 3, "Final two-compartment model" column
    lka  <- log(0.53);  label("First-order oral absorption rate constant Ka (1/h)")        # Rambiritch 2016 Table 3 (Final two-compartment): Ka = 0.53 1/h, RSE 8.33%
    lcl  <- log(2.16);  label("Apparent clearance CL/F (L/h)")                              # Rambiritch 2016 Table 3 (Final two-compartment): CL/F = 2.16 L/h, RSE 7.41%
    lvc  <- log(11.70); label("Apparent central volume of distribution Vc/F (L)")           # Rambiritch 2016 Table 3 (Final two-compartment): V2/F = 11.70 L, RSE 9.49% (paper labels the central volume V2/F)
    lq   <- log(3.84);  label("Apparent intercompartmental clearance Q/F (L/h)")            # Rambiritch 2016 Table 3 (Final two-compartment): Q/F = 3.84 L/h, RSE 14.97%
    lvp  <- log(68.10); label("Apparent peripheral volume of distribution Vp/F (L)")        # Rambiritch 2016 Table 3 (Final two-compartment): V3/F = 68.10 L, RSE 8.81% (paper labels the peripheral volume V3/F)
    lfdepot <- fixed(log(1)); label("Bioavailability F (fixed at 1; absolute F not identifiable from oral-only data)")  # Apparent parameterisation (CL/F, V/F, Q/F): F is not separately estimated in Rambiritch 2016; log(1) FIXED is the standard apparent-parameter anchor.

    # Inter-individual variability (BSV) -- Rambiritch 2016 Methods Eq 1: P_j = TVP * exp(eta_j); paper reports
    # sqrt(omega^2) as an approximation of CV (Methods paragraph following Eq 1), so omega^2 = (BSV%/100)^2.
    etalka  ~ 0.08161  # BSV 28.57% CV on Ka (Rambiritch 2016 Table 3); omega^2 = 0.2857^2 = 0.08161
    etalcl  ~ 0.11499  # BSV 33.91% CV on CL/F (Rambiritch 2016 Table 3); omega^2 = 0.3391^2 = 0.11499
    etalvc  ~ 0.05308  # BSV 23.04% CV on Vc/F (Rambiritch 2016 Table 3); omega^2 = 0.2304^2 = 0.05308
    etalq   ~ 0.42706  # BSV 65.35% CV on Q/F (Rambiritch 2016 Table 3); omega^2 = 0.6535^2 = 0.42706
    # No eta on lvp: BSV on V3/F reported as 0.02% (Rambiritch 2016 Table 3), and Table 4 shows V3/F = 68.10 L
    # for every individual subject, indicating omega_V3 was effectively constrained to zero in the source NONMEM run.

    # Residual error -- Rambiritch 2016 Methods Eq 2 and Table 3 "Residual variability"
    # Concentration data were log-transformed prior to fitting (LTBS); the additive residual on log-scale
    # is equivalent to a proportional residual in linear space with sigma ~= CV for small sigma.
    # Reported variance = 0.189 (43.5% CV); sqrt(0.189) = 0.4347.
    propSd <- 0.4347; label("Proportional residual SD (fraction; LTBS-equivalent)")          # Rambiritch 2016 Table 3: residual variability variance 0.189 (43.5% CV); sqrt(0.189) = 0.4347.
  })

  model({
    # Individual PK parameters
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp)
    fdepot <- exp(lfdepot)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system -- Rambiritch 2016 Figure 1 model equations
    #   dA_central/dt    = Ka * depot - (CL/V_c + Q/V_c) * A_central + (Q/V_p) * A_peripheral
    #   dA_peripheral/dt = (Q/V_c) * A_central - (Q/V_p) * A_peripheral
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Bioavailability anchor on the depot compartment (apparent parameterisation: F = 1)
    f(depot) <- fdepot

    # Concentration: dose in mg, vc in L -> central/vc in mg/L = ng/uL.
    # Multiply by 1000 to convert to ng/mL so Cc matches the assay units used in Rambiritch 2016
    # (Cmax reported in ng/mL; Table 2 columns "Cmax (ng/mL)" and "AUCinf (ng h/mL)").
    Cc <- (central / vc) * 1000
    Cc ~ prop(propSd)
  })
}
