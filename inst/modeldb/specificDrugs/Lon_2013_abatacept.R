Lon_2013_abatacept <- function() {
  description <- "Two-compartment population PK model with linear elimination and short-term zero-order SC absorption for abatacept (CTLA-4Ig Fc-fusion) in male Lewis rats with collagen-induced arthritis (Lon 2013)."
  reference <- "Lon HK, Liu D, DuBois DC, Almon RR, Jusko WJ. Modeling pharmacokinetics/pharmacodynamics of abatacept and disease progression in collagen-induced arthritic rats: a population approach. J Pharmacokinet Pharmacodyn. 2013;40(6):701-712. doi:10.1007/s10928-013-9341-1"
  vignette <- "Lon_2013_abatacept"
  units <- list(time = "day", dosing = "mg", concentration = "umol/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Lon 2013 reports PK parameters normalized to body weight (CL, CLD in mL/day/kg; V1, V2 in mL/kg). WT is used here to scale those to absolute units (mL/day and mL) for a dose given in mg. The study enrolled male Lewis rats weighing 150-175 g at arrival; no individual body-weight covariate effect was modeled in the source.",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 17L,
    n_studies      = 1L,
    species        = "Male Lewis rat (Rattus norvegicus), collagen-induced arthritis (CIA) model",
    age_range      = "6-9 weeks old at arrival; dosing initiated on day 21 post-collagen induction",
    weight_range   = "150-175 g at arrival",
    sex_female_pct = 0,
    disease_state  = "Collagen-induced arthritis (CIA) with paw swelling >= 50% increase in one or two hind paws by day 20 post-induction",
    dose_range     = "10 mg/kg IV single dose (n=5); 20 mg/kg SC single dose (n=6); 20 mg/kg SC on day 21 followed by 10 mg/kg SC on days 23, 25, 27, 29 (n=6)",
    regions        = "United States (University at Buffalo)",
    notes          = "Fifty male Lewis rats purchased; 23 arthritic rats with >=50% paw-swelling response were randomized to vehicle (single n=3, multiple n=3), IV 10 mg/kg (n=5), SC single 20 mg/kg (n=6), or SC multiple (n=6). Seventeen abatacept-treated rats contributed to the population PK fit. Plasma abatacept was quantified by sCTLA-4 ELISA (LLOQ 0.16 ng/mL, inter-day CV ~15%). Per Lon 2013 Methods (Animals) and Experimental Design."
  )

  ini({
    # Molecular weight used to convert mg/mL in the central compartment to the
    # umol/L reporting units of Lon 2013. The paper reports concentrations in
    # umol/L (Table 1 epsilon2 = 0.0365 umol/L; IC50 = 0.731 umol/L "circa
    # 67.3 ug/mL" in the Discussion), which corresponds to MW ~92000 g/mol for
    # abatacept and matches the Orencia label value of ~92,300 Da. The MW is
    # not stated explicitly in the paper; it is used here only for unit
    # conversion and is declared fixed.
    mw_abatacept <- fixed(92000); label("Molecular weight of abatacept (g/mol)") # Lon 2013 Discussion: "IC50 value (0.731 umol/L, circa 67.3 ug/mL)" implies MW ~92 kDa

    # Structural PK parameters - Lon 2013 Table 1 final estimates.
    # Body-weight-normalized: CL and CLD in mL/day/kg; V1 and V2 in mL/kg.
    # Scaled to absolute units inside model() by multiplying by WT (kg).
    lcl     <- log(21.8);  label("Clearance CL normalized to body weight (mL/day/kg)")               # Lon 2013 Table 1: CL = 21.8 mL/day/kg (%RSE 6.70)
    lvc     <- log(69.5);  label("Central compartment volume V1 normalized to body weight (mL/kg)")  # Lon 2013 Table 1: V1 = 69.5 mL/kg (%RSE 21.6)
    lq      <- log(27.5);  label("Distributional clearance CLD normalized to body weight (mL/day/kg)") # Lon 2013 Table 1: CLD = 27.5 mL/day/kg (%RSE 10.3)
    lvp     <- log(61.9);  label("Peripheral compartment volume V2 normalized to body weight (mL/kg)") # Lon 2013 Table 1: V2 = 61.9 mL/kg (%RSE 9.00)
    lfdepot <- log(0.592); label("SC bioavailability F (fraction)")                                    # Lon 2013 Table 1: F = 59.2% (%RSE 9.20)
    ltau    <- log(2.67);  label("SC zero-order input duration tau (day)")                            # Lon 2013 Table 1: tau = 2.67 day (%RSE 7.80)

    # Inter-individual variability. Lon 2013 Methods: "inter-animal variability
    # of the parameters was modeled with an exponential function" (theta_i =
    # theta * exp(eta_i)). Table 1 reports IIV as mean % with %RSE, taken as
    # the CV% of the log-normal distribution; omega^2 = log(CV^2 + 1).
    # CL IIV 9.67% -> omega^2 = log(0.0967^2 + 1) = 0.009307
    # V1 IIV 56.2% -> omega^2 = log(0.562^2 + 1)  = 0.274478
    etalcl ~ 0.009307     # Lon 2013 Table 1: IIV on CL = 9.67% (%RSE 47.2)
    etalvc ~ 0.274478     # Lon 2013 Table 1: IIV on V1 = 56.2% (%RSE 70.3)

    # Residual error - Lon 2013 Table 1 reports a combined proportional +
    # additive model. Proportional epsilon1 is reported as a percent (16.1%,
    # used as fraction 0.161). Additive epsilon2 is reported directly in
    # umol/L (0.0365 umol/L) - the paper's concentration unit.
    propSd  <- 0.161;  label("Proportional residual error (fraction)")                   # Lon 2013 Table 1: epsilon1 = 16.1% (%RSE 28.0)
    addSd   <- 0.0365; label("Additive residual error (umol/L)")                         # Lon 2013 Table 1: epsilon2 = 0.0365 umol/L (%RSE 45.0)
  })

  model({
    # Individual PK parameters. CL, CLD, V1, V2 are scaled from body-weight-
    # normalized (per-kg) values to absolute (mL/day and mL) using WT (kg).
    cl     <- exp(lcl + etalcl) * WT
    vc     <- exp(lvc + etalvc) * WT
    q      <- exp(lq)           * WT
    vp     <- exp(lvp)          * WT
    fdepot <- exp(lfdepot)
    tau    <- exp(ltau)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Short-term zero-order SC absorption. A dose to 'depot' releases drug
    # into 'central' at a constant bioavailable rate of F * Dose / tau for
    # a duration of tau (Lon 2013 Eq. 3). podo(depot) returns the prescribed
    # dose of the most recent depot dose; tad(depot) is time since that dose.
    # After tau, kzero drops to zero and any residual (1 - F) * Dose stays
    # in depot as the un-absorbed fraction. IV doses are given to 'central'
    # and bypass this mechanism.
    kzero <- fdepot * podo(depot) / tau
    if (tad(depot) > tau) kzero <- 0.0

    d/dt(depot)       <- -kzero
    d/dt(central)     <-  kzero - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-          k12 * central - k21 * peripheral1

    # Concentration: dose in mg and vc in mL give central/vc in mg/mL = g/L.
    # Convert to umol/L by dividing by MW (g/mol) and multiplying by 1e6
    # (mol -> umol). Matches the units of propSd / addSd in ini() and the
    # values in Lon 2013 Table 1.
    Cc <- (central / vc) * 1e6 / mw_abatacept
    Cc ~ add(addSd) + prop(propSd)
  })
}
