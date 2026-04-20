Masters_2022_avelumab <- function() {
  description <- "Two-compartment population PK model for avelumab (anti-PD-L1 IgG1) with time-dependent clearance in patients with advanced solid tumors (Masters 2022)"
  reference <- "Masters JC, Khandelwal A, di Pietro A, Dai H, Brar S. Model-informed drug development supporting the approval of the avelumab flat-dose regimen in patients with advanced renal cell carcinoma. CPT Pharmacometrics Syst Pharmacol. 2022;11(4):458-468. doi:10.1002/psp4.12771"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling on CL, V1, V2, and Q. Allometric exponent fixed at 1 for Q; exponents on CL, V1, V2 are estimated. Reference weight 80 kg (paper repeatedly cites ~80 kg as the median body weight used to select the 800 mg flat dose; reference weight inherited from the Wilkins 2019 base structural model is not explicitly stated in Masters 2022 but is consistent with 80 kg given the reported typical CL).",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 2315L,
    n_studies      = "5 pooled trials (488 aRCC patients from JAVELIN Renal 100 and 101 combination arm + 1827 solid-tumor monotherapy patients from the prior Wilkins 2019 dataset)",
    age_range      = "not tabulated in the main-text (pooled solid-tumor adult oncology population)",
    weight_range   = "aRCC sub-population 44.2-143.0 kg; prior solid-tumor monotherapy sub-population 30.4-204 kg",
    weight_median  = "aRCC sub-population 81.5 kg; prior solid-tumor monotherapy sub-population 70.6 kg",
    disease_state  = "Advanced / metastatic solid tumors (metastatic Merkel cell carcinoma, advanced / metastatic urothelial carcinoma, advanced renal cell carcinoma, and other solid tumors)",
    dose_range     = "10 mg/kg IV every 2 weeks (weight-based clinical regimen); 800 mg IV every 2 weeks simulated as the flat-dose alternative",
    regions        = "Global (multinational phase Ib / III oncology trials)",
    notes          = "Final-model pooled population per Masters 2022 Results 'Population PK analysis'. aRCC sub-population baseline demographics are in Tables S2 / S3 of the supplementary material."
  )

  ini({
    # Structural parameters referenced to an 80 kg adult. CL and Q are
    # reported in Masters 2022 Table 1 in L/h; converted to L/day (x 24)
    # because T50 is reported in days and this model keeps time in days.
    lcl    <- log(0.0269 * 24);  label("Baseline clearance CL for an 80 kg adult (L/day)")             # Masters 2022 Table 1: theta_CL = 0.0269 L/h
    lvc    <- log(3.196);        label("Central volume V1 for an 80 kg adult (L)")                     # Masters 2022 Table 1: theta_V1 = 3.196 L
    lvp    <- log(0.7278);       label("Peripheral volume V2 for an 80 kg adult (L)")                  # Masters 2022 Table 1: theta_V2 = 0.7278 L
    # Masters 2022 Table 1 prints theta_Q = 0.3352 L/h, but its printed RSE
    # (12.24%) and 95% CI (0.02548-0.04157 L/h) are only internally
    # consistent with theta_Q = 0.03352 L/h (a missing leading zero in the
    # printed point estimate). 0.03352 L/h also matches the Wilkins 2019
    # precedent from which the structural model was inherited.
    lq     <- log(0.03352 * 24); label("Intercompartmental clearance Q for an 80 kg adult (L/day)")    # Masters 2022 Table 1 theta_Q, corrected via its own 95% CI to 0.03352 L/h

    # Time-dependent clearance modifier (Hill-type):
    #   CL(t) = CL_baseline * (1 + Imax * t^gamma / (T50^gamma + t^gamma))
    # Imax is reported as -0.08533 (fractional decrease in CL at t >> T50).
    # Parameterized here via log|Imax| so the magnitude is log-normally
    # distributed with IIV and the negative sign is applied in model(),
    # keeping Imax < 0 across all simulated individuals.
    lImax  <- log(0.08533);      label("log|Imax| - magnitude of the fractional change in CL at t >> T50 (unitless)") # Masters 2022 Table 1: theta_Imax = -0.08533 (sign applied in model())
    lt50   <- log(99.24);        label("log T50 - time at which half of Imax is reached (days)")       # Masters 2022 Table 1: theta_T50 = 99.24 days
    lgamma <- log(2.086);        label("log gamma - Hill shape coefficient of the time-on-CL function (unitless)")    # Masters 2022 Table 1: theta_gamma = 2.086

    # Allometric exponents on body weight (reference 80 kg). Exponent on Q
    # is fixed to 1 per Masters 2022 Methods ("Study overview").
    e_wt_cl <- 0.4714; label("Allometric exponent on CL (unitless)")                                   # Masters 2022 Table 1: theta_weight_on_CL = 0.4714
    e_wt_v1 <- 0.4694; label("Allometric exponent on V1 (unitless)")                                   # Masters 2022 Table 1: theta_weight_on_V1 = 0.4694
    e_wt_v2 <- 0.5826; label("Allometric exponent on V2 (unitless)")                                   # Masters 2022 Table 1: theta_weight_on_V2 = 0.5826

    # Inter-individual variability. CL, V1, V2 form a 3x3 log-normal block;
    # Imax has an independent log-normal eta on its magnitude. Variance /
    # covariance values (lower-triangular, row-major order):
    #   row 1:  omega^2_CL                                = 0.09339
    #   row 2:  cov(CL,V1),  omega^2_V1                   = 0.03048, 0.03776
    #   row 3:  cov(CL,V2),  cov(V1,V2),  omega^2_V2      = 0.08418, 0.01799, 1.204
    # Source: Masters 2022 Table 1 (rows omega^2_CL, omega^2_V1, omega^2_V2,
    # cov_CL-V1, cov_CL-V2, cov_V1-V2).
    etalcl + etalvc + etalvp ~ c(
      0.09339,
      0.03048, 0.03776,
      0.08418, 0.01799, 1.204
    )
    etalImax ~ 0.1052  # Masters 2022 Table 1: omega^2_Imax = 0.1052

    # Residual error (combined proportional + additive). The table column
    # is "sigma" (standard deviation, not variance).
    propSd <- 0.1742; label("Proportional residual error SD (fraction)")                               # Masters 2022 Table 1: sigma_proportional_error = 0.1742
    addSd  <- 2.168;  label("Additive residual error SD (ug/mL)")                                      # Masters 2022 Table 1: sigma_additive_error = 2.168
  })
  model({
    # Allometric scaling (reference 80 kg); exponent on Q is fixed at 1.
    cl_base <- exp(lcl + etalcl) * (WT / 80)^e_wt_cl
    v1      <- exp(lvc + etalvc) * (WT / 80)^e_wt_v1
    v2      <- exp(lvp + etalvp) * (WT / 80)^e_wt_v2
    q       <- exp(lq)           * (WT / 80)

    # Time-dependent clearance modifier (Hill function of time since first dose).
    gamma  <- exp(lgamma)
    t50    <- exp(lt50)
    Imax_i <- -exp(lImax + etalImax)
    cl     <- cl_base * (1 + Imax_i * t^gamma / (t50^gamma + t^gamma))

    # Two-compartment micro-constants.
    kel <- cl / v1
    k12 <- q  / v1
    k21 <- q  / v2

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose is in mg and V1 in L, so central/V1 has units mg/L = ug/mL.
    Cc <- central / v1
    Cc ~ prop(propSd) + add(addSd)
  })
}
