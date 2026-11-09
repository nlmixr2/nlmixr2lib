Weiss_2023_niclosamide_sheep <- function() {
  description <- "Preclinical (sheep). Three-compartment population PK model for plasma niclosamide in ewes after intravenous infusion and nebulized (inhaled) niclosamide ethanolamine, with first-order absorption from a lung absorption compartment, estimated nebulized bioavailability, linear elimination, and combined proportional and additive residual error (Weiss 2023)."
  reference <- "Weiss A, Bischof RJ, Landersdorfer CB, Nguyen TH, Davies A, Ibrahim J, Wynne P, Wright P, Ditzinger G, Montgomery AB, Meeusen E, McIntosh MP, Sommer MOA. Single-dose pharmacokinetics and lung function of nebulized niclosamide ethanolamine in sheep. Pharm Res. 2023;40(8):1915-1925. doi:10.1007/s11095-023-03559-0"
  vignette <- "Weiss_2023_niclosamide_sheep"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "niclosamide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "niclosamide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "niclosamide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "niclosamide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "sheep (Merino cross-bred ewes)",
    n_subjects = 10L,
    n_studies = 1L,
    age_range = "12 months (all animals)",
    weight_range = "25-35 kg (nebulized cohort, Methods 'Animals'); 40.6 +/- 2.1 kg (intravenous / oral cohort, Table I footnote)",
    sex_female_pct = 100,
    disease_state = "Healthy",
    dose_range = "Nebulized niclosamide ethanolamine (NEN) 10, 30 and 60 mg (6 mL of 0.167%, 0.5% and 1% NEN) delivered over 24-27 min; intravenous NEN 3 mg (3 mL of 0.1% NEN) infused over 30 min (Table I, Table S1).",
    regions = "Monash Animal Research Platform, Monash University, Australia.",
    notes = "Two cohorts of five ewes (Table I). Experiment 1 (sheep 1-5) received nebulized placebo and the three ascending nebulized doses with 7-day washouts; Experiment 2 (sheep 11-14, 18) received an oral 1132 mg niclosamide dose and the intravenous dose. The population PK model was fitted to the plasma concentrations after the intravenous and nebulized doses only; the oral arm and the epithelial lining fluid concentrations were analysed by non-compartmental methods and are not part of this model. Dose nebulized administration with the DELIVERED dose (Table S1: 0.271, 0.762 and 1.401 mg/kg NEN), not the nominal nebulizer charge: the paper's own Figure S3 simulations sit 0.6-0.75 times below a nominal-dose simulation, matching the delivery loss (see vignette). Doses are in mg NEN as the paper reports them; the paper does not state whether a salt correction to niclosamide base was applied before modelling."
  )

  ini({
    # Final-model estimates, Weiss 2023 Supplementary Table S2. Doses are
    # entered in mg of niclosamide ethanolamine (NEN) as the paper reports
    # them (Table I, Table S1); concentrations are ng/mL niclosamide.
    lka <- log(20.9); label("First-order absorption rate constant from the lung absorption compartment (1/h)")  # Table S2 ka = 20.9 h-1 (also Table II)
    lfdepot <- log(0.858); label("Bioavailability of nebulized dosing (fraction)")  # Table S2 Fneb = 85.8% (Table II: 86%)
    lcl <- log(170); label("Systemic elimination clearance (L/h)")  # Table S2 CL = 170 L/h
    lvc <- log(15.7); label("Central volume of distribution (L)")  # Table S2 V1 = 15.7 L
    lvp <- log(10.1); label("Shallow peripheral volume of distribution (L)")  # Table S2 V2 = 10.1 L
    lq <- log(57.4); label("Distribution clearance, central to shallow peripheral (L/h)")  # Table S2 CLd1 = 57.4 L/h
    lvp2 <- log(148); label("Deep peripheral volume of distribution (L)")  # Table S2 V3 = 148 L
    lq2 <- log(7.47); label("Distribution clearance, central to deep peripheral (L/h)")  # Table S2 CLd2 = 7.47 L/h

    # Log-normal IIV reported as %CV; omega^2 = log(CV^2 + 1). No covariance
    # is reported in Table S2, so the etas are independent.
    etalka ~ 0.455905  # Table S2 IIV ka = 76% CV
    etalcl ~ 0.109392  # Table S2 IIV CL = 34% CV
    etalvc ~ 0.408195  # Table S2 IIV V1 = 71% CV
    etalvp ~ 0.0754785  # Table S2 IIV V2 = 28% CV

    propSd <- 0.252; label("Proportional residual error (fraction)")  # Table S2 CV_CP = 25.2%
    addSd <- 0.26; label("Additive residual error (ng/mL)")  # Table S2 SD_CP = 0.26 ng/mL
  })
  model({
    ka <- exp(lka + etalka)
    fdepot <- exp(lfdepot)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq)
    vp2 <- exp(lvp2)
    q2 <- exp(lq2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Nebulized doses enter depot (the paper's absorption compartment Aabs)
    # as a zero-order input over the nebulization time; intravenous doses are
    # zero-order infusions into central. The paper writes Fneb on the
    # depot -> central flux (ka * Fneb * Aabs); applying it to the depot dose
    # gives the same central-compartment amounts.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    f(depot) <- fdepot

    # mg / L = ug/mL; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
