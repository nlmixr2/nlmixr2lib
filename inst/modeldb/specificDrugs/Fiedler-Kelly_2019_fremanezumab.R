`Fiedler-Kelly_2019_fremanezumab` <- function() {
  description <- "Two-compartment population PK model for fremanezumab (anti-CGRP IgG2 delta-a/kappa mAb) with first-order SC absorption and allometric weight scaling (Fiedler-Kelly 2019)"
  reference <- "Fiedler-Kelly JB, Cohen-Barak O, Morris DN, et al. Population pharmacokinetic modelling and simulation of fremanezumab in healthy subjects and patients with migraine. Br J Clin Pharmacol. 2019;85(12):2721-2733. doi:10.1111/bcp.14096"
  vignette <- "Fiedler-Kelly_2019_fremanezumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling on CL (exponent 1.05) and Vc (exponent 1.53) with reference weight 71 kg (population median). Treated as a time-invariant baseline covariate in the source paper.",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects        = 2546,
    n_studies         = 7,
    age_range         = "18-71 years",
    age_median        = "43 years",
    weight_range      = "43.5-131.8 kg",
    weight_median     = "70.8 kg",
    sex_female_pct    = 86.1,
    race_ethnicity    = c(Caucasian = 79.9, Other = 20.1),
    disease_state     = "Healthy adults (n = 74) and adults with chronic migraine or episodic migraine (n = 2474)",
    dose_range        = "225, 675, or 900 mg SC (Q4W or Q12W); 225 or 900 mg IV in 1 phase 1 study",
    regions           = "Global (phase 1 included Japanese and Caucasian cohorts)",
    ada_positive_pct  = 0.7,
    analysis_dataset  = "13745 fremanezumab concentrations pooled from 2 phase 1, 2 phase 2b, and 3 phase 3 studies",
    notes             = "Baseline demographics from Fiedler-Kelly 2019 Section 3.1 and Supporting Table S1."
  )

  ini({
    # Structural parameters - reference values for a 71 kg subject (population median)
    # Units: CL and Q in L/day; Vc and Vp in L; ka in 1/day; tlag in day
    lka     <- log(0.180);  label("Absorption rate constant (ka, 1/day)")                             # Fiedler-Kelly 2019 Table 2
    lcl     <- log(0.0902); label("Central clearance for a 71 kg subject (CL, L/day)")                # Fiedler-Kelly 2019 Table 2
    lvc     <- log(1.88);   label("Central volume of distribution (SC) for a 71 kg subject (Vc, L)")  # Fiedler-Kelly 2019 Table 2
    lvp     <- fixed(log(1.72));   label("Peripheral volume of distribution (Vp, L)")                 # Fiedler-Kelly 2019 Table 2 (FIXED in source)
    lq      <- fixed(log(0.262));  label("Intercompartmental clearance (Q, L/day)")                   # Fiedler-Kelly 2019 Table 2 (FIXED in source)
    lfdepot <- fixed(log(0.658));  label("SC absolute bioavailability (F, fraction)")                 # Fiedler-Kelly 2019 Table 2 (FIXED in source)
    ltlag   <- fixed(log(0.0803)); label("SC absorption lag time (tlag, day)")                        # Fiedler-Kelly 2019 Table 2 (FIXED in source)

    # Allometric exponents on weight (reference 71 kg)
    allo_cl <- 1.05; label("Allometric exponent on CL (unitless)")                                    # Fiedler-Kelly 2019 Table 2, footnote c
    allo_v  <- 1.53; label("Allometric exponent on Vc (unitless)")                                    # Fiedler-Kelly 2019 Table 2, footnotes d and e

    # IIV - diagonal omega matrix (no off-diagonal terms estimated per Fiedler-Kelly 2019 Section 3.2)
    # omega^2 = log(CV^2 + 1) since BSV was modelled using an exponential form
    # CL: 23.4%CV -> omega^2 = log(0.234^2 + 1) = 0.05336
    # Vc: 35.1%CV -> omega^2 = log(0.351^2 + 1) = 0.11620
    # ka: 59.0%CV -> omega^2 = log(0.590^2 + 1) = 0.29858
    etalcl ~ 0.05336                                                                                  # Fiedler-Kelly 2019 Table 2 (23.4%CV)
    etalvc ~ 0.11620                                                                                  # Fiedler-Kelly 2019 Table 2 (35.1%CV)
    etalka ~ 0.29858                                                                                  # Fiedler-Kelly 2019 Table 2 (59.0%CV)

    # Residual error - SC combined additive + proportional model (the primary therapeutic route)
    # Source reports variance-scale NONMEM $SIGMA entries; nlmixr2 expects SD.
    # Proportional variance 0.0531 -> propSd = sqrt(0.0531) = 0.2304
    # Additive variance 0.204 (ug/mL^2) -> addSd = sqrt(0.204) = 0.4517 ug/mL
    CcpropSd <- 0.2304; label("Proportional residual error (SC, fraction)")                           # Fiedler-Kelly 2019 Table 2 (variance 0.0531)
    CcaddSd  <- 0.4517; label("Additive residual error (SC, ug/mL)")                                  # Fiedler-Kelly 2019 Table 2 (variance 0.204)
  })

  model({
    # Individual PK parameters with allometric weight scaling (reference 71 kg)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 71)^allo_cl
    vc <- exp(lvc + etalvc) * (WT / 71)^allo_v
    vp <- exp(lvp)          * (WT / 71)^allo_v
    q  <- exp(lq)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system - 2-compartment with first-order SC absorption
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # SC bioavailability and lag time on the depot compartment
    f(depot)     <- exp(lfdepot)
    alag(depot)  <- exp(ltlag)

    # Observation: dose in mg, volume in L -> mg/L = ug/mL
    Cc <- central / vc
    Cc ~ add(CcaddSd) + prop(CcpropSd)
  })
}
