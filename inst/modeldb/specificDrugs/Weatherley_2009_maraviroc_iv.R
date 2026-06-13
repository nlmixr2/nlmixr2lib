Weatherley_2009_maraviroc_iv <- function() {
  description <- paste(
    "Four-compartment IV maraviroc population PK in 20 healthy young adult",
    "males receiving 3, 10, or 30 mg as a 1-hour IV infusion (Weatherley &",
    "McFadyen 2009 Br J Clin Pharmacol). NONMEM ADVAN7 + FOCE-I fit to log-",
    "transformed plasma concentrations from study A4001009 only. Exponential",
    "inter-subject variability on CL, V1, V2, Q3 and V3; proportional residual",
    "error (additive on the log scale). No dose effect on clearance over the",
    "3-30 mg range. This file extracts only the IV disposition analysis",
    "(paper Analysis 1, Table 3); the paper's two companion analyses (sigmoid",
    "Emax NAUC meta-regression of oral phase 1 data, and the S-PLUS closed-form",
    "mass balance model) are not ODE PK models and are described in the",
    "vignette but not implemented here."
  )
  reference <- paste(
    "Weatherley B, McFadyen L (2009).",
    "Maraviroc modelling strategy: use of early phase 1 data to support a",
    "semi-mechanistic population pharmacokinetic model.",
    "Br J Clin Pharmacol 68(5):648-657.",
    "doi:10.1111/j.1365-2125.2009.03455.x.",
    sep = " "
  )
  vignette <- "Weatherley_2009_maraviroc_iv"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 20L,
    n_studies      = 1L,
    age_range      = "22-42 years",
    weight_range   = "65-93 kg",
    sex_female_pct = 0,
    disease_state  = "Healthy adult male volunteers (Pfizer study A4001009, open ascending IV solution).",
    dose_range     = "3, 10, or 30 mg maraviroc as a 1-hour intravenous solution infusion (single dose).",
    notes          = paste(
      "473 plasma concentrations from 20 subjects across the three IV dose levels.",
      "Concentrations measured by SPE-LC-MS-MS (LLOQ 0.5 ng/mL).",
      "See Weatherley 2009 Table 1 (row A4001009, IV-only analysis)."
    )
  )

  ini({
    # Structural fixed effects (Table 3, Estimate column).
    # All values are population means from the final 4-compartment fit to
    # 473 log-transformed plasma concentrations from 20 subjects.
    lcl  <- log(48.0);  label("Clearance CL (L/h)")                        # Table 3 (48.0 L/h, 2.2% SE)
    lvc  <- log(13.5);  label("Central volume V1 (L)")                     # Table 3 (13.5 L, 14.4% SE)
    lvp  <- log(16.1);  label("First peripheral volume V2 (L)")            # Table 3 (16.1 L, 10.5% SE)
    lvp2 <- log(140);   label("Second peripheral volume V3 (L)")           # Table 3 (140 L, 5.8% SE)
    lvp3 <- log(36.4);  label("Third peripheral volume V4 (L)")            # Table 3 (36.4 L, 4.0% SE)
    lq   <- log(35.4);  label("Intercompartmental clearance Q2 (L/h)")     # Table 3 (35.4 L/h, 15.3% SE)
    lq2  <- log(9.13);  label("Intercompartmental clearance Q3 (L/h)")     # Table 3 (9.13 L/h, 4.8% SE)
    lq3  <- log(16.8);  label("Intercompartmental clearance Q4 (L/h)")     # Table 3 (16.8 L/h, 8.9% SE)

    # IIV (log-normal). NONMEM Table 3 reports the CV% on the exponential
    # eta (hCL, hV1, hV2, hQ3, hV3). Converted to log-scale variance via
    # omega^2 = log(CV^2 + 1):
    #   CL : 10.8% -> log(1 + 0.108^2) = 0.01160
    #   V1 : 25.6% -> log(1 + 0.256^2) = 0.06348
    #   V2 : 18.0% -> log(1 + 0.180^2) = 0.03189
    #   Q3 : 18.8% -> log(1 + 0.188^2) = 0.03473
    #   V3 : 16.9% -> log(1 + 0.169^2) = 0.02816
    # No IIV reported on V4, Q2, or Q4 in Table 3; etas on those parameters
    # are not declared (they take the population value).
    etalcl  ~ 0.01160                                                      # Table 3 hCL = 10.8% CV
    etalvc  ~ 0.06348                                                      # Table 3 hV1 = 25.6% CV
    etalvp  ~ 0.03189                                                      # Table 3 hV2 = 18.0% CV
    etalq2  ~ 0.03473                                                      # Table 3 hQ3 = 18.8% CV
    etalvp2 ~ 0.02816                                                      # Table 3 hV3 = 16.9% CV

    # Residual error. Methods (Model development for intravenous data,
    # paper page 650): "proportional residual variability was modelled
    # using additive form on log transformed data." Table 3 reports
    # 10.9% (7.1% SE). This is the standard NONMEM EPS(1) additive-on-
    # log-Y form that is equivalent to proportional error in linear
    # space; encode as propSd = 0.109.
    propSd <- 0.109; label("Proportional residual error (fraction)")       # Table 3 Residual = 10.9%
  })

  model({
    # Individual structural parameters
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc  + etalvc)
    vp  <- exp(lvp  + etalvp)
    vp2 <- exp(lvp2 + etalvp2)
    vp3 <- exp(lvp3)
    q   <- exp(lq)
    q2  <- exp(lq2  + etalq2)
    q3  <- exp(lq3)

    # Micro-constants for the linear 4-compartment IV ODE
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    k14 <- q3 / vc
    k41 <- q3 / vp3

    # Four-compartment IV disposition (NONMEM ADVAN7 equivalent).
    # IV-infusion doses enter `central` directly via the data set's
    # RATE / DUR column (1-hour infusion in the source study).
    d/dt(central)     <- -(kel + k12 + k13 + k14) * central +
                          k21 * peripheral1 +
                          k31 * peripheral2 +
                          k41 * peripheral3
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(peripheral3) <- k14 * central - k41 * peripheral3

    # Maraviroc plasma concentration in ng/mL.
    # central is in mg (matching the mg dosing convention); vc is in L.
    # Cc [ng/mL] = (central [mg] / vc [L]) * 1000.
    Cc <- 1000 * central / vc

    Cc ~ prop(propSd)
  })
}
