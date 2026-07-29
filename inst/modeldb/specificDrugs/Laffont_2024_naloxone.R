Laffont_2024_naloxone <- function() {
  description <- "Population PK model for intranasal (IN) naloxone HCl in healthy adult volunteers (Laffont 2024): two-compartment model with linear elimination and parallel zero-order plus lagged first-order absorption; Q/F and Vp/F fixed to literature values from Yassen 2007."
  reference <- "Laffont CM, Purohit P, Delcamp N, Gonzalez-Garcia I, Skolnick P. Comparison of intranasal naloxone and intranasal nalmefene in a translational model assessing the impact of synthetic opioid overdose on respiratory depression and cardiac arrest. Front Psychiatry. 2024;15:1399803. doi:10.3389/fpsyt.2024.1399803. Q/F and Vp/F fixed to values from Yassen A, Olofsen E, van Dorp E, Sarton E, Teppema L, Danhof M, Dahan A. Mechanism-based pharmacokinetic-pharmacodynamic modelling of the reversal of buprenorphine-induced respiratory depression by naloxone. Clin Pharmacokinet. 2007;46(11):965-980. doi:10.2165/00003088-200746110-00004"
  vignette <- "Laffont_2024_naloxone"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    # Laffont 2024 Table 2 reports no statistically significant covariates retained
    # in the IN naloxone final model (in particular, body weight is not on CL/F).
    # No covariate columns are required to simulate from this model.
  )

  population <- list(
    n_subjects     = "TODO: not reported in main text; the IN naloxone population PK dataset comes from the pharmacodynamic study by Ellison et al. 2024 (ref 23). Detail is in Supplementary Table 2 (not on disk).",
    n_studies      = 1,
    age_range      = "Adult (healthy volunteers); detailed range in Supplementary Table 2 (not on disk).",
    weight_median  = "TODO: not reported in main text",
    sex_female_pct = "TODO: not reported in main text",
    race_ethnicity = "TODO: not reported in main text",
    disease_state  = "Healthy adult volunteers",
    dose_range     = "Single 4 mg IN naloxone HCl (commercial 4 mg/0.1 mL nasal spray) in pharmacodynamic study",
    regions        = "TODO: not reported in main text",
    notes          = "IN naloxone PK data come from the pharmacodynamic study by Ellison et al. 2024 (ref 23) under hypercapnic-gas-mixture conditions. The paper notes that the hypercapnic-mask conditions did not affect IN naloxone absorption during the first 20 minutes post dose when compared with published data in healthy volunteers (refs 10, 23), so the population PK parameters in Table 2 were used as-is for opioid-overdose rescue simulations."
  )

  ini({
    # Structural parameters -- Laffont 2024 Table 2 (IN naloxone)
    lcl    <- log(396);          label("Apparent clearance (CL/F, L/h)")                          # Table 2: CL/F = 396 L/h
    lvc    <- log(65.7);         label("Apparent central volume of distribution (Vc/F, L)")       # Table 2: Vc/F = 65.7 L
    lq     <- fixed(log(284));   label("Apparent intercompartmental clearance (Q/F, L/h)")        # Table 2: Q/F = 284 (fixed to Yassen 2007 ref 30)
    lvp    <- fixed(log(102));   label("Apparent peripheral volume of distribution (Vp/F, L)")    # Table 2: Vp/F = 102 (fixed to Yassen 2007 ref 30)
    lka    <- log(0.998);        label("First-order absorption rate constant (KA, 1/h)")          # Table 2: KA = 0.998 1/h
    ld2    <- log(0.689);        label("Zero-order absorption duration (D2, h)")                  # Table 2: D2 = 0.689 h
    lfk0   <- log(0.183);        label("Fraction of intranasal dose absorbed via zero-order route (FK0)")  # Table 2: FK0 = 0.183
    lalag1 <- log(0.0717);       label("Lag time of first-order absorption (ALAG1, h)")            # Table 2: ALAG1 = 0.0717 h

    # IIV -- Table 2 reports CV%; log-normal variance = log(1 + CV^2)
    etalcl  ~ log(1 + 0.391^2)                                                                   # Table 2: IIV CL/F = 39.1 %CV
    etalvc  ~ log(1 + 2.40^2)                                                                    # Table 2: IIV Vc/F = 240 %CV
    etalfk0 ~ log(1 + 1.51^2)                                                                    # Table 2: IIV FK0 = 151 %CV

    # Residual error -- Table 2 reports sigma^2 = 0.104 (32.3 %CV) on the log-additive scale
    # (proportional in linear space); propSd = sqrt(sigma^2) = 0.323.
    propSd <- 0.323; label("Proportional residual error (fraction)")                             # Table 2: sigma^2 = 0.104 (32.3 %CV)
  })

  model({
    # Individual PK parameters
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    ka   <- exp(lka)
    d2   <- exp(ld2)
    fk0  <- exp(lfk0 + etalfk0)
    alag <- exp(lalag1)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Parallel absorption: a single intranasal dose is split between two simultaneous
    # input pathways, controlled by f() and dur(). The user supplies two dose records
    # at the same time with the same total amt -- one to depot for first-order absorption
    # (rate = 0, normal bolus) and one to central for zero-order infusion (rate = -2 to
    # invoke the modeled dur(central)). See the validation vignette for an event-table
    # example.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-               k12 * central - k21 * peripheral1

    # First-order route: bioavailable fraction of IN dose, with absorption lag
    f(depot)   <- 1 - fk0
    lag(depot) <- alag

    # Zero-order route: fraction of IN dose infused directly into central over D2 hours
    f(central)   <- fk0
    dur(central) <- d2

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL; multiply by 1000 for ng/mL
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
