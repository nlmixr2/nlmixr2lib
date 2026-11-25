Badhan_2017_liposome_mouse_1cmt <- function() {
  description <- "Preclinical (mouse). Empirical 1-compartment model of DDA:TDB cationic-liposome adjuvant kinetics at the IM site of injection (quadriceps), from a cationic-liposome TB subunit-antigen vaccine (Badhan 2017 Table 2)"
  reference <- "Badhan RKS, Khadke S, Perrie Y. Application of Pharmacokinetics Modelling to Predict Human Exposure of a Cationic Liposomal Subunit Antigen Vaccine System. Pharmaceutics. 2017 Dec 6;9(4):57. doi:10.3390/pharmaceutics9040057. PMID 29215597; PMCID PMC5750663."
  vignette <- "Badhan_2017_liposomal_vaccine"
  units <- list(time = "day", dosing = "ug", concentration = "ug")

  # Empirical compartmental model (Equation 1). The single central compartment
  # is the IM site of injection (quadriceps); there is no reported volume, so
  # the apparent volume is 1 and Cc is the site-of-injection liposome amount (ug).
  compartmentData <- list(
    central = list(analyte = "DDA:TDB liposome", units = "ug", specimen = "administration site", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "mouse (female, 28 g)",
    n_subjects = NA_integer_,
    n_studies = 1,
    age_range = NA,
    weight_range = "28 g body weight",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state = "Healthy mice dosed IM with a dual-radiolabelled DDA:TDB cationic liposome adsorbing the Ag85B-ESAT-6 (H1) TB subunit antigen; site-of-injection liposome kinetics.",
    dose_range = "300 ug liposome IM",
    regions = NA,
    notes = "Empirical model selected by AIC: a 1-compartment model best described the liposome (Table 2; k10 only, terminal half-life 22.6 days). The paper also develops a minimal-PBPK model (Equations 7-10); it is described in the vignette but not packaged because it is not reproducible from the printed (dimensionally inconsistent) equations. Companion compartmental model: Badhan_2017_antigen_mouse_2cmt."
  )

  ini({
    lk10 <- log(0.0306); label("log site-of-injection elimination rate constant k10 (day^-1)")  # Table 2 (Liposome): K10 = 0.0306 day^-1 (t1/2 22.6 days)
  })

  model({
    k10 <- exp(lk10)
    # Equation 1: dC1/dt = -k10 * C1. Apparent volume 1; Cc is the central amount.
    d/dt(central) <- -k10 * central
    Cc <- central
  })
}
