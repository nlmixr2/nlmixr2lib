Kratochwil_2018_antiIL17_mouse <- function() {
  description <- "Preclinical (mouse, C57BL/6J). Linear two-compartment population PK model for an anti-IL17 human IgG1 model antibody (RO5310074) after a single IV bolus, fit to 14C-nanotracer plasma concentrations measured by cavity ring-down spectroscopy and accelerator mass spectrometry; parameters are body-weight normalised (per kg)."
  reference <- "Kratochwil NA, Dueker SR, Muri D, Senn C, Yoon H, Yu BY, Lee GH, Dong F, Otteneder MB. Nanotracing and cavity-ring down spectroscopy: A new ultrasensitive approach in large molecule drug disposition studies. PLoS One. 2018;13(10):e0205435. doi:10.1371/journal.pone.0205435"
  vignette <- "Kratochwil_2018_antiIL17_mouse"
  # Doses and volumes are per kg body weight, as in Table 1. Concentration is
  # molar (the paper's Figure 4B/C axis), using MW 150000 g/mol.
  units <- list(time = "day", dosing = "nmol/kg", concentration = "nM")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "anti-IL17 IgG1 (RO5310074)", units = "nmol/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "anti-IL17 IgG1 (RO5310074)", units = "nmol/kg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse (C57BL/6J)",
    n_subjects = NA,
    n_studies = 1,
    sex_female_pct = 0,
    weight_range = "about 0.02 kg (adult; the paper scales its tracer dose to a 0.02 kg mouse)",
    disease_state = "Healthy male adult C57BL/6J mice (no IL-17 target-mediated disposition expected at this dose).",
    dose_range = "Single IV bolus (tail vein) of 75 nmol/kg anti-IL17 containing 2.1 nCi of 14C-propionate-labelled anti-IL17 (specific activity 1.24 nCi/nmol for the blended dose).",
    regions = "Switzerland (Roche, Basel); bioanalysis in South Korea",
    notes = "Composite (one blood sample per mouse per time) plasma sampling at 1, 3, 7, 24, 48, 72 and 96 h plus terminal samples at 7, 24, 168 and 336 h. The number of mice is not printed. The prediction-corrected VPC is drawn for 'all data', which the paper compares with the 125I-labelled anti-IL17 microtracer profiles (66 nmol/kg) of Eigenmann et al.; whether those gamma-counting data entered the fit is not stated. Estimation in Phoenix NLME 6.4 (FOCE-ELS)."
  )

  ini({
    # Table 1 of Kratochwil 2018 reports body-weight-normalised estimates in
    # mL/kg and mL/day/kg; they are kept in those units (central / vc in
    # nmol/mL is converted to nM in model()).
    lcl <- log(6.69); label("Clearance (mL/day/kg)") # Table 1: CL = 6.69 mL day-1 kg-1 (RSE 3.68%)
    lvc <- log(49.5); label("Central volume of distribution (mL/kg)") # Table 1: Vc = 49.5 mL kg-1 (RSE 6.05%)
    lq <- log(84.3); label("Intercompartmental clearance (mL/day/kg)") # Table 1: Q = 84.3 mL day-1 kg-1 (RSE 8.05%)
    lvp <- log(44.7); label("Peripheral volume of distribution (mL/kg)") # Table 1: Vp = 44.7 mL kg-1 (RSE 8.42%)

    # Only Vp carries inter-subject variability (Methods 'PK modeling', Eq 8):
    # no other eta could be estimated from the composite profiles.
    etalvp ~ 0.0871 # Table 1: omega2 (Vp) = 0.0871 (RSE 32.6%)

    propSd <- 0.153; label("Proportional residual error (fraction)") # Table 1: sigma proportional = 0.153 (RSE 11.6%); Eq 9 Cobs = Cpred*(1+eps)
  })
  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    # Eqs 5-7: dAc/dt = -CL*Cc - Q*(Cc - Cp); dAp/dt = Q*(Cc - Cp)
    d / dt(central) <- -cl / vc * central - q / vc * central + q / vp * peripheral1
    d / dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    # central (nmol/kg) / vc (mL/kg) is nmol/mL; x 1000 gives nmol/L (nM).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
