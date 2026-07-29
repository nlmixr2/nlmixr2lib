Jones_2011_PF04878691_lymphocyte <- function() {
  description <- paste(
    "Coupled PK + indirect-response pharmacodynamic model for absolute",
    "lymphocyte count during oral PF-04878691 (TLR7 agonist) administration",
    "in healthy adult volunteers (Jones 2011 BJCP). PK is the two-compartment",
    "time-varying clearance model from the companion file",
    "Jones_2011_PF04878691.R (Table 1; all PK structural parameters and",
    "IIVs fixed at the published Table 1 values so the PK forcing function",
    "is the published popPK profile). Drug stimulates the re-distribution",
    "(loss) of lymphocytes through a power function on kout:",
    "dLYMPH/dt = kin - kout * (1 + slope * Cc^gamma) * LYMPH, with baseline",
    "lymphocyte count rbase = kin / kout so that kin = rbase * kout",
    "(Jones 2011 Methods / Table 3 lymphocyte model). The typical Emax",
    "indirect-response model could not adequately identify the parameters",
    "given the limited number of dose levels studied, so the Emax * Cc^gamma",
    "/ (EC50^gamma + Cc^gamma) drug effect was replaced with the power",
    "function slope * Cc^gamma (Methods 'Population PK-OAS and",
    "PK-lymphocyte models')."
  )
  reference <- paste(
    "Jones HM, Chan PLS, van der Graaf PH, Webster R. Use of modelling and",
    "simulation techniques to support decision making on the progression of",
    "PF-04878691, a TLR7 agonist being developed for hepatitis C.",
    "Br J Clin Pharmacol. 2012;73(1):77-92.",
    "doi:10.1111/j.1365-2125.2011.04047.x.",
    "Companion PK file Jones_2011_PF04878691.R supplies the upstream",
    "two-compartment time-varying clearance PK structure."
  )
  vignette <- "Jones_2011_PF04878691_HCV"
  paper_specific_compartments <- c("lymph")
  paper_specific_etas <- c("etalrbase", "etalkout", "etalslope")
  paper_specific_residual_sds <- c("propSd_lymph")
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required scaling covariate for the inherited per-kg PK parameters; see Jones_2011_PF04878691.R. Population median 79 kg.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    age_range      = "21-55 years",
    age_median     = "34 years",
    weight_range   = "57-97 kg",
    weight_median  = "79 kg",
    sex_female_pct = 8,
    race_ethnicity = "Not tabulated in Jones 2011.",
    disease_state  = "Healthy adult volunteers (multiple-dose escalation Phase 1 study; ClinicalTrials.gov NCT00810758).",
    dose_range     = "PF-04878691 administered orally at 3, 6, or 9 mg twice weekly (days 1, 4, 8, 11) for 2 weeks; n = 6 active per dose cohort plus n = 2 placebo per cohort.",
    regions        = "Not specified.",
    notes          = "Absolute lymphocyte count measured by immunophenotyping. Sequential population fit: individual EBE PK parameters from the Table 1 PK model were used as the exposure driver for the lymphocyte fit. Transient dose- and time-dependent decreases in absolute lymphocyte count were observed in the 6 and 9 mg dose groups."
  )

  ini({
    # ----------------------------------------------------------------------
    # PK structural parameters from the companion Jones_2011_PF04878691.R
    # PK file (Table 1). All four structural PK theta entries are wrapped
    # in fixed() because the lymphocyte model was fit sequentially using
    # the individual EBE PK parameters from the Table 1 popPK fit -- the
    # PK parameters are inherited, not re-estimated. Reparameterised from
    # the paper's (CLF, CL0, DEG) into the canonical (cl_ss, cl_time, kdeg)
    # form per the companion PK file.
    # ----------------------------------------------------------------------
    lcl       <- fixed(log(1.7));    label("Steady-state apparent clearance per kg body weight (CL_SS = paper CLF, L/h/kg)")              # Table 1 (CLF = 1.7 L/h/kg)
    lcl_time  <- fixed(log(1.8));    label("Initial offset of the time-varying clearance component per kg (CL_TIME0 = CL0 - CLF, L/h/kg)") # Derived from Table 1 (CL0 = 3.5, CLF = 1.7)
    lkdeg     <- fixed(log(0.24));   label("Exponential decay rate of the time-varying clearance component (paper DEG, 1/h)")              # Table 1 (DEG = 0.24 1/h)
    lvc       <- fixed(log(3.3));    label("Apparent central volume of distribution per kg body weight (Vc, L/kg)")                        # Table 1 (Vc = 3.3 L/kg)
    lq        <- fixed(log(0.74));   label("Apparent intercompartmental clearance per kg body weight (Q, L/h/kg)")                         # Table 1 (Q = 0.74 L/h/kg)
    lvp       <- fixed(log(21));     label("Apparent peripheral volume of distribution per kg body weight (Vp, L/kg)")                     # Table 1 (Vp = 21 L/kg)
    lka       <- fixed(log(0.078));  label("First-order absorption rate constant (ka, 1/h)")                                               # Table 1 (ka = 0.078 1/h)

    # PK IIVs from Table 1 -- also fixed because the lymphocyte fit used
    # EBE PK parameters from the upstream PK fit.
    etalcl ~ fixed(0.067)                                                                                                                  # Table 1 (IIV CLF = 0.067)
    etalka ~ fixed(0.19)                                                                                                                   # Table 1 (IIV ka  = 0.19)

    # ----------------------------------------------------------------------
    # Lymphocyte indirect-response PD parameters (Jones 2011 Table 3). The
    # baseline lymphocyte count rbase together with kout fix the
    # production rate via kin = rbase * kout. Drug stimulates loss
    # (kout) through the power function slope * Cc^gamma.
    # Note on baseline units: Table 3 prints 'pg ml^-1' for the BASE row,
    # which is a paper typo (lymphocytes were measured by immunophenotyping
    # as an absolute count). The numerical value 1890 is the typical
    # absolute lymphocyte count in cells/uL; downstream users should treat
    # rbase as cells/uL (within the normal range 1000-4500 cells/uL).
    # ----------------------------------------------------------------------
    lkout    <- log(0.044);  label("Lymphocyte first-order elimination rate constant kout (1/h)")                                            # Table 3 (kout = 0.044 1/h, %CV 20)
    lslope   <- log(0.44);   label("Slope of the drug stimulation of lymphocyte loss per (ng/mL)^gamma (paper SLP, units 1 / (ng/mL)^gamma)") # Table 3 (SLP = 0.44, %CV 15)
    lrbase   <- log(1890);   label("Lymphocyte baseline absolute count (rbase = kin/kout, cells/uL; Table 3 unit label 'pg/mL' is a paper typo)") # Table 3 (BASE = 1890 cells/uL, %CV 4.7)
    lgamma   <- log(2.2);    label("Sigmoidicity exponent on Cc in the lymphocyte stimulation power function (gamma, unitless)")             # Table 3 (gamma = 2.2, %CV 17)

    # Lymphocyte IIVs from Table 3 (multiplicative exponential random
    # effects). Three IIVs: kout, slope, baseline. No IIV on gamma.
    etalkout  ~ 0.19                                                                                                                        # Table 3 (IIV kout  = 0.19,  %CV 83)
    etalslope ~ 0.20                                                                                                                        # Table 3 (IIV SLP   = 0.20,  %CV 53)
    etalrbase ~ 0.051                                                                                                                       # Table 3 (IIV BASE  = 0.051, %CV 29)

    # Lymphocyte residual error (proportional). Table 3 reports
    # residual error = 0.021 (%CV 21).
    propSd_lymph <- 0.021; label("Lymphocyte proportional residual SD (fraction of count)")                                                  # Table 3 (residual error = 0.021, %CV 21)
  })

  model({
    # ----------------------------------------------------------------------
    # Inherited PK (deterministic + per-subject log-normal IIV on cl_ss
    # and ka, exactly as Table 1 reports). See Jones_2011_PF04878691.R for
    # the standalone PK model rationale.
    # ----------------------------------------------------------------------
    cl_ss    <- exp(lcl + etalcl) * WT
    cl_time0 <- exp(lcl_time) * WT
    kdeg_cl  <- exp(lkdeg)
    vc       <- exp(lvc) * WT
    q        <- exp(lq)  * WT
    vp       <- exp(lvp) * WT
    ka       <- exp(lka + etalka)

    cl  <- cl_ss + cl_time0 * exp(-kdeg_cl * time)
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    Cc <- 1000 * central / vc   # PF-04878691 plasma concentration in ng/mL

    # ----------------------------------------------------------------------
    # Lymphocyte indirect-response PD (Jones 2011 Methods 'Population PKPD
    # analysis', lymphocyte model). Drug stimulates loss (kout) through a
    # power function; baseline kin / kout balance sets the lymphocyte
    # state at rbase.
    #     dLYMPH/dt = kin - kout * (1 + slope * Cc^gamma) * LYMPH
    #     kin       = rbase * kout
    #     LYMPH(0)  = rbase
    # ----------------------------------------------------------------------
    rbase       <- exp(lrbase + etalrbase)
    kout_lymph  <- exp(lkout  + etalkout)
    slope_lymph <- exp(lslope + etalslope)
    gamma_lymph <- exp(lgamma)

    kin_lymph <- rbase * kout_lymph

    lymph(0) <- rbase

    d/dt(lymph) <- kin_lymph - kout_lymph * (1 + slope_lymph * Cc^gamma_lymph) * lymph

    lymph ~ prop(propSd_lymph)
  })
}
