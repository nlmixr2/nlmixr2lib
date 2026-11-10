Kinsale_2024_islatravir_rat <- function() {
  description <- "Preclinical (rat). Two-compartment PK model with first-order absorption for plasma islatravir after a single subcutaneous bolus (0.1, 0.3 and 1 mg/kg) in adult female Wistar rats, fitted naive-pooled in the macro-constant (A, alpha, B, beta, Ka) parameterisation and used as the unit impulse response for deconvolving absorption from a reservoir-style poly(caprolactone) islatravir implant"
  reference <- paste(
    "Kinsale TS, Cottrell ML, Li L, Brand R, Gatto G, Luecke E, Norton C,",
    "Krovi A, Dumond JB, Rao G, Yeshwante S, Van Horne B, Van Der Straten A,",
    "Kashuba ADM, Johnson LM.",
    "Pharmacokinetic Modeling to Guide Preclinical Development of an",
    "Islatravir-Eluting Reservoir-Style Biodegradable Implant for Long-Acting",
    "HIV PrEP. Pharmaceutics. 2024;16(2):201.",
    "doi:10.3390/pharmaceutics16020201",
    sep = " "
  )
  vignette <- "Kinsale_2024_islatravir_rat"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # No covariates: the naive-pooled fit has a single set of structural
  # parameters for all three dose levels. Doses were given per kg of body
  # weight, so the per-animal dose in mg (dose level x body weight) is carried
  # by the event table rather than by a covariate.
  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "islatravir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "islatravir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "islatravir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "rat (Wistar)",
    n_subjects = 72,
    n_studies = 1,
    weight_range = "0.307 +/- 0.02 kg (mean +/- SD)",
    weight_median = "0.307 kg (mean)",
    sex_female_pct = 100,
    disease_state = "healthy adult female rats",
    dose_range = "0.1, 0.3 and 1 mg/kg single subcutaneous bolus (islatravir in DMSO, under the neck skin)",
    regions = "USA (Magee-Womens Research Institute, Pittsburgh)",
    notes = paste(
      "Kinsale 2024 Sections 2.2-2.3: 72 rats, 24 per dose level, sparse",
      "destructive design with 3 rats per time point (0.17, 0.5, 1, 2, 4, 6",
      "and 24 h) and one plasma sample per animal. 57 quantifiable plasma",
      "islatravir concentrations; 6 BLQ values (Beal M2) and 1 Grubbs",
      "outlier were excluded (Section 3.1). The mean administered dose over",
      "all animals was 145,803 ng (Section 2.6). The model was then used as",
      "the unit impulse response to deconvolve published plasma PK from 10",
      "Wistar rats carrying an islatravir-eluting implant (98.6 mg loaded)."
    )
  )

  ini({
    # Phoenix NLME naive-pooled two-compartment macro model (Supplementary
    # Table S1, model 6 'Base Model no Outlier' = the final model). The fit
    # reports the macro-constants A and B in ng/mL for the pooled data at
    # all three dose levels (a single predicted curve, Figure 2). Section
    # 2.6 converts them to the unit impulse response by dividing by the mean
    # administered dose Dref = 145,803 ng, so the dose-normalised disposition
    # is the biexponential (A/Dref) exp(-alpha t) + (B/Dref) exp(-beta t)
    # (ng/mL per ng). That impulse response is written here with the
    # canonical two-compartment identities so that the model is dose
    # proportional (Section 3.1 found Cmax and AUClast dose proportional):
    #
    #   vc  = Dref / (A + B)
    #   k21 = (A * beta + B * alpha) / (A + B)
    #   kel = alpha * beta / k21
    #   k12 = alpha + beta - k21 - kel
    #
    # with Table S1 Final A = 141.60 ng/mL and B = 1.39 ng/mL. The absorption
    # rate Ka feeds the impulse response through a first-order depot.
    lka <- log(12.23); label("First-order SC absorption rate constant Ka (1/h)") # Table S1 'Final' Ka = 12.23 1/h (CV 16.83%); Figure 2 inset Ka = 12 (CV 17%)
    lalpha <- log(1.21); label("Distribution-phase hybrid rate constant alpha (1/h)") # Table S1 'Final' alpha = 1.21 1/h (CV 4.6%); Figure 2 inset alpha = 1.2; Section 3.2 alpha1 = 1.2
    lbeta <- log(0.103); label("Terminal-phase hybrid rate constant beta (1/h)") # Table S1 'Final' beta = 0.103 1/h (CV 6.6%); Figure 2 inset beta = 0.103; Section 3.2 alpha2 = 0.103
    lvc <- log(1.01967); label("Central volume of distribution (L), Dref / (A + B)") # derived: Section 2.6 Dref = 145,803 ng; Table S1 'Final' A = 141.60 ng/mL, B = 1.39 ng/mL; 145803 / 142.99 = 1019.67 mL
    lk21 <- log(0.11376); label("Peripheral to central rate constant k21 (1/h), (A beta + B alpha) / (A + B)") # derived: Table S1 'Final' A = 141.60, B = 1.39 ng/mL, alpha = 1.21, beta = 0.103 1/h

    # Naive-pooled fit: there is no between-animal variability by design
    # (Section 2.6, 'naive pooled algorithm'), so no etas are declared.

    # Phoenix multiplicative error C * (1 + eps) is proportional error.
    propSd <- 0.76; label("Proportional residual error (fraction)") # Section 3.2 '76% residual error'; Table S1 model 6 RSE column = 0.76
  })

  model({
    # 1. Individual parameters (typical values; naive-pooled fit)
    ka <- exp(lka)
    alpha <- exp(lalpha)
    beta <- exp(lbeta)
    vc <- exp(lvc)
    k21 <- exp(lk21)

    # 2. Micro-constants from alpha * beta = kel * k21 and
    #    alpha + beta = kel + k12 + k21
    kel <- alpha * beta / k21
    k12 <- alpha + beta - k21 - kel

    # Clearance-parameterisation equivalents (derived outputs)
    cl <- kel * vc
    q <- k12 * vc
    vp <- q / k21

    # 3. ODE system: SC bolus into depot; first-order absorption to central
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Observation: central in mg, vc in L -> mg/L; x 1000 -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
