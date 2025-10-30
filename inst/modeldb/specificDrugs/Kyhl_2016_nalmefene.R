Kyhl_2016_nalmefene <- function() {
  reference <- "Kyhl LE, Li S, Faerch KU, Soegaard B, Larsen F, Areberg J. Population pharmacokinetics of nalmefene in healthy subjects and its relation to μ-opioid receptor occupancy. Br J Clin Pharmacol. 2016 Feb;81(2):290-300. doi: 10.1111/bcp.12805. Epub 2016 Jan 27. PMID: 26483076; PMCID: PMC4833148."
  covariates <-
    list(
      LBM = "Lean body mass, LBM (kg)",
      AGE = "Age (years)",
      TABLET = "1 = tablet formulation, 0 = solution formulation",
      FED = "1 = fed, 0 = fasted",
      RIA_ASSAY = "1 = radioimmunoassay (RIA) used for analysis; 0 = LC-MS/MS assay used for analysis"
    )
  units <-
    list(
      Dose = "mg",
      Concentration = "ng/mL"
    )
  # Notes:
  # * Assuming that IIV% = sqrt(exp(variance) - 1)*100 (the typical, correct
  #   method, not one of the shortcut methods)
  # * Values were switched to log-scale from linear scale to standardize the
  #   model presentation with other models. This does not affect the model
  #   results.
  # * V2 and V3 were renamed to vc and vp to align with modeling used in
  #   `nlmixr2lib`
  # * All values are from Table 3
  ini({
    lka_tablet <- log(0.751); label("Absorption rate for oral tablet")
    iiv_lka_tablet ~ sqrt(log((69.9/100)^2 + 1))
    lka_solution <- log(1.4)

    lcl <- log(60.4); label("Clearance (L/h)")
    iiv_lcl ~ sqrt(log((18.7/100)^2 + 1))
    allo_lbm_cl <- 0.626; label("Lean body mass on CL (power)")
    lvc <- log(266); label("Volume of distribution, central compartment (L)")
    iiv_lvc ~ sqrt(log((66.6/100)^2 + 1))
    e_age_vc <- -2.11; label("Effect of age on central volume of distribution (L/year)")
    lq <- log(109); label("Inter‐compartmental clearance (L/h)")
    iiv_q ~ sqrt(log((44.5/100)^2 + 1))
    lvp <- log(537); label("Volume of distribution, peripheral compartment (L)")
    iiv_lvp ~ sqrt(log((45.8/100)^2 + 1))

    lf_oral <- log(0.406); label("Absolute oral bioavailability")
    iiv_f_oral ~ sqrt(log((20.1/100)^2 + 1))
    e_fed_f_oral <- 0.294

    propSd <- 0.094; label("Proportional residual error (fraction)")
    addSd_ria <- 0.0255; label("Additive residual error for radioimmunoassay (RIA)")
    addSd_lcms <- fixed(0); label("Additive residual error for LC-MS/MS")
    
    emax_mu_opioid <- 99.4; label("Maximum effect on mu-opioid receptor occupancy (percent)")
    lec50_mu_opioid <- log(0.338); label("Concentration producing 50% of maximum effect on mu-opioid receptor occupancy (ng/mL)")
  })
  model({
    ka_tablet <- exp(lka_tablet + iiv_lka_tablet)
    ka_solution <- exp(lka_solution)
    ka <- ka_tablet*TABLET + ka_solution*(1 - TABLET)
    cl <- exp(lcl + iiv_lcl) * (LBM/56.28)^allo_lbm_cl
    vc <- exp(lvc + iiv_lvc) + e_age_vc*(AGE - 28)
    q <- exp(lq + iiv_q)
    vp <- exp(lvp + iiv_lvp)

    f_oral <- exp(lf_oral + iiv_f_oral) * (1 + e_fed_f_oral*FED)

    kel <- cl/vc
    k12 <- q/vc
    k21 <- q/vp

    d/dt(depot) <- -ka*depot
    d/dt(central) <-  ka*depot - kel*central - k12*central + k21*peripheral1
    d/dt(peripheral1) <- k12*central - k21*peripheral1
    
    f(depot) <- f_oral
    # Unit conversion from mg/L to ng/mL
    Cc <- central / vc * 1000

    addSd <- addSd_ria*RIA_ASSAY + addSd_lcms*(1 - RIA_ASSAY)
    Cc ~ add(addSd) + prop(propSd)

    # These were not estimated within the primary model    
    ec50_mu_opioid <- exp(lec50_mu_opioid)
    e_mu_opioid <- emax_mu_opioid*Cc/(ec50_mu_opioid + Cc)
  })
}
