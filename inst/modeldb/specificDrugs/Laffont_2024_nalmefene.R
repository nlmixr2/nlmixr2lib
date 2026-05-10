Laffont_2024_nalmefene <- function() {
  description <- "Population PK model for intranasal (IN) nalmefene HCl in healthy adult volunteers (Laffont 2024): two-compartment model with linear elimination, parallel zero-order plus lagged first-order absorption, and allometric body-weight scaling on apparent clearance."
  reference <- "Laffont CM, Purohit P, Delcamp N, Gonzalez-Garcia I, Skolnick P. Comparison of intranasal naloxone and intranasal nalmefene in a translational model assessing the impact of synthetic opioid overdose on respiratory depression and cardiac arrest. Front Psychiatry. 2024;15:1399803. doi:10.3389/fpsyt.2024.1399803"
  vignette <- "Laffont_2024_nalmefene"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power (allometric) scaling on apparent clearance with reference weight 74.7 kg (population median in nalmefene PK dataset).",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = "TODO: not reported in main text; population PK dataset pooled three healthy-volunteer studies (Crystal 2024 ref 22 [two PK studies, one nostril vs. two nostrils] and Ellison 2024 ref 23 [pharmacodynamic remifentanil-induced respiratory depression study]). Detail is in Supplementary Tables 1-2 (not on disk).",
    n_studies      = 3,
    age_range      = "Adult (healthy volunteers); detailed range in Supplementary Table 1 (not on disk).",
    weight_median  = "74.7 kg (median body weight reported as the allometric reference value in Table 1)",
    sex_female_pct = "TODO: not reported in main text",
    race_ethnicity = "TODO: not reported in main text",
    disease_state  = "Healthy adult volunteers",
    dose_range     = "Single 3 mg IN nalmefene HCl per nostril (1 or 2 nostrils, total 3 or 6 mg) and single 1 mg IM nalmefene HCl in PK studies; 3 mg IN nalmefene HCl in pharmacodynamic study under hypercapnic mask",
    regions        = "TODO: not reported in main text",
    notes          = "Pooled population PK analysis from Crystal et al. 2024 and Ellison et al. 2024. Body weight on apparent clearance has reference 74.7 kg. The pharmacodynamic study under a hypercapnic gas mixture mask reduced first-order intranasal absorption by ~35% (estimated parameter STDEFF = -0.349 on INKA); for opioid-overdose rescue simulations the authors used absorption parameters estimated outside the hypercapnic mask, i.e., STDEFF = 0. The Laffont_2024_nalmefene model file follows that rescue-setting convention."
  )

  ini({
    # Structural parameters -- Laffont 2024 Table 1 (IN nalmefene, reference WT = 74.7 kg)
    lcl     <- log(63.7);   label("Apparent clearance (CL/F, L/h)")                                  # Table 1: CL/F = 63.7 L/h
    e_wt_cl <- 0.572;       label("Allometric exponent of (WT/74.7) on CL/F (unitless)")             # Table 1: Exponent of (WT/74.7) for CL/F = 0.572
    lvc     <- log(15.2);   label("Apparent central volume of distribution (Vc/F, L)")               # Table 1: Vc/F = 15.2 L
    lq      <- log(81.3);   label("Apparent intercompartmental clearance (Q/F, L/h)")                # Table 1: Q/F = 81.3 L/h
    lvp     <- log(522);    label("Apparent peripheral volume of distribution (Vp/F, L)")            # Table 1: Vp/F = 522 L
    linka   <- log(0.497);  label("Intranasal first-order absorption rate constant (INKA, 1/h)")     # Table 1: INKA = 0.497 1/h
    ld2     <- log(0.302);  label("Zero-order absorption duration (D2, h)")                          # Table 1: D2 = 0.302 h
    linfk0  <- log(0.0485); label("Fraction of intranasal dose absorbed via zero-order route (INFK0)")# Table 1: INFK0 = 0.0485
    lalag1  <- log(0.0615); label("Lag time of first-order absorption (ALAG1, h)")                   # Table 1: ALAG1 = 0.0615 h

    # IIV -- Table 1 reports CV%; log-normal variance = log(1 + CV^2)
    etalcl   ~ log(1 + 0.154^2)                                                                     # Table 1: IIV CL/F = 15.4 %CV
    etalvc   ~ log(1 + 2.11^2)                                                                      # Table 1: IIV Vc/F = 211 %CV
    etalinka ~ log(1 + 0.398^2)                                                                     # Table 1: IIV INKA = 39.8 %CV

    # Residual error -- Table 1 reports sigma^2 = 0.111 (33.3 %CV) on the log-additive scale,
    # which is equivalent to a proportional residual error (NONMEM "additive on log-scale" ==
    # proportional in linear space); propSd is the SD = sqrt(sigma^2) = 0.333.
    propSd <- 0.333; label("Proportional residual error (fraction)")                                # Table 1: sigma^2 = 0.111 (33.3 %CV)
  })

  model({
    # Individual PK parameters -- IN nalmefene, rescue setting (no hypercapnic mask, so STDEFF = 0)
    cl    <- exp(lcl + etalcl) * (WT / 74.7)^e_wt_cl
    vc    <- exp(lvc + etalvc)
    q     <- exp(lq)
    vp    <- exp(lvp)
    inka  <- exp(linka + etalinka)
    d2    <- exp(ld2)
    infk0 <- exp(linfk0)
    alag  <- exp(lalag1)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Parallel absorption: a single intranasal dose is split between two simultaneous
    # input pathways, controlled by f() and dur(). The user supplies two dose records
    # at the same time with the same total amt -- one to depot for first-order absorption
    # (rate = 0, normal bolus) and one to central for zero-order infusion (rate = -2 to
    # invoke the modeled dur(central)). See the validation vignette for an event-table
    # example.
    d/dt(depot)       <- -inka * depot
    d/dt(central)     <-  inka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                  k12 * central - k21 * peripheral1

    # First-order route: bioavailable fraction of IN dose, with absorption lag
    f(depot)   <- 1 - infk0
    lag(depot) <- alag

    # Zero-order route: small fraction of IN dose infused directly into central over D2 hours
    f(central)   <- infk0
    dur(central) <- d2

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL; multiply by 1000 for ng/mL
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
