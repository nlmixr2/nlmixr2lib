Huang_2017_vrc01 <- function() {
  description <- "Two-compartment population PK model for VRC01 (HIV-1 broadly neutralizing IgG1 monoclonal antibody) in healthy adults after IV or SC administration (Huang 2017)"
  reference <- "Huang Y, Zhang L, Ledgerwood J, et al. Population pharmacokinetics analysis of VRC01, an HIV-1 broadly neutralizing monoclonal antibody, in healthy adults. MAbs. 2017;9(5):792-800. doi:10.1080/19420862.2017.1311435"
  vignette <- "Huang_2017_vrc01"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Centered at 74.5 kg (median weight of IV groups in HVTN104). Exponential effect on CL and Vc; power-form (allometric) effect on Q and Vp.",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 84L,
    n_studies      = 1L,
    age_range      = "18-50 years",
    age_median     = "27 years",
    weight_range   = "53-114 kg",
    weight_median  = "72 kg overall; 74.5 kg in IV groups (used as the covariate-centering reference)",
    sex_female_pct = 50,
    race_ethnicity = "Not reported in Huang 2017 Table 1",
    disease_state  = "HIV-uninfected (low HIV-1 risk) healthy adults",
    dose_range     = "10-40 mg/kg IV every 4 or 8 weeks; 5 mg/kg SC every 2 weeks (with a 40 mg/kg IV loading dose in Group 3)",
    regions        = "United States (HVTN 104, NCT02165267)",
    notes          = "Baseline demographics from Huang 2017 Table 1. The HVTN 104 phase 1 study enrolled 42 men and 42 women across five dosing groups; the popPK model was fit to 1117 VRC01 serum concentrations. External validation used VRC602 data."
  )

  ini({
    # Structural parameters — population-typical values at the median IV-group
    # weight of 74.5 kg. CL and Vc are reference values at WT = 74.5 kg; Q and
    # Vp are reference values at the same allometric-scaling reference.
    lka     <- log(0.26); label("Absorption rate constant after SC administration (Ka, 1/day)")  # Huang 2017 Table 2 (final model, IV+SC)
    lcl     <- log(0.40); label("Clearance from the central compartment at WT = 74.5 kg (CL, L/day)")  # Huang 2017 Table 2 (final model, IV+SC)
    lvc     <- log(1.94); label("Volume of the central compartment at WT = 74.5 kg (Vc, L)")  # Huang 2017 Table 2 (final model, IV+SC)
    lq      <- log(0.84); label("Intercompartmental distribution clearance at WT = 74.5 kg (Q, L/day)")  # Huang 2017 Table 2 (final model, IV+SC)
    lvp     <- log(4.90); label("Volume of the peripheral compartment at WT = 74.5 kg (Vp, L)")  # Huang 2017 Table 2 (final model, IV+SC)
    lfdepot <- log(0.74); label("SC bioavailability relative to IV administration (F1, fraction)")  # Huang 2017 Table 2 (final model, IV+SC)

    # Body-weight covariate effects.
    # Form is given by Huang 2017 Methods ("Covariate model"):
    #   exponential: TV_theta = theta * exp(beta * (BW - 74.5))
    #   power:       TV_theta = theta * (BW / 74.5)^beta
    # Final-model selections (Table 2 footnote 2): exponential on CL and Vc;
    # power on Q and Vp.
    e_wt_cl <- 0.012; label("Body-weight effect on CL (exponential, fold/kg)")  # Huang 2017 Table 2
    e_wt_vc <- 0.010; label("Body-weight effect on Vc (exponential, fold/kg)")  # Huang 2017 Table 2
    e_wt_q  <- 0.69;  label("Body-weight effect on Q (power exponent, unitless)")  # Huang 2017 Table 2
    e_wt_vp <- 0.82;  label("Body-weight effect on Vp (power exponent, unitless)")  # Huang 2017 Table 2

    # Inter-individual variability (omega^2 reported on the log-normal internal
    # scale; Huang 2017 reports CV% as sqrt(omega^2) * 100). CL, Q, and Vp form
    # a correlated 3x3 block (off-diagonals = covariances from Table 2). Vc has
    # a separate diagonal element with no reported covariances. F1 and ka have
    # IIV variances fixed at 0 (Table 2 footnote 3).
    etalcl + etalq + etalvp ~ c(0.067,
                                0.034, 0.063,
                                0.050, 0.082, 0.120)  # Huang 2017 Table 2
    etalvc ~ 0.028  # Huang 2017 Table 2

    # Combination proportional + additive residual error (Huang 2017 Methods,
    # "Variability popPK model"): C_obs = C_pred * (1 + e1) + e2.
    # Reported variances are sigma1^2 = 0.042 and sigma2^2 = 0.456.
    propSd <- sqrt(0.042); label("Proportional residual error (fraction)")  # Huang 2017 Table 2
    addSd  <- sqrt(0.456); label("Additive residual error (ug/mL)")  # Huang 2017 Table 2
  })
  model({
    # Body-weight effects (exponential on CL/Vc, power on Q/Vp).
    wt_cl <- exp(e_wt_cl * (WT - 74.5))
    wt_vc <- exp(e_wt_vc * (WT - 74.5))
    wt_q  <- (WT / 74.5)^e_wt_q
    wt_vp <- (WT / 74.5)^e_wt_vp

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * wt_cl
    vc <- exp(lvc + etalvc) * wt_vc
    q  <- exp(lq  + etalq)  * wt_q
    vp <- exp(lvp + etalvp) * wt_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
