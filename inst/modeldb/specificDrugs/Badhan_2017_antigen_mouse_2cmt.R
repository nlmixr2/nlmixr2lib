Badhan_2017_antigen_mouse_2cmt <- function() {
  description <- "Preclinical (mouse). Empirical 2-compartment model of Ag85B-ESAT-6 (H1) subunit-antigen kinetics at the IM site of injection (quadriceps), from a DDA:TDB cationic-liposome TB vaccine (Badhan 2017 Table 2)"
  reference <- "Badhan RKS, Khadke S, Perrie Y. Application of Pharmacokinetics Modelling to Predict Human Exposure of a Cationic Liposomal Subunit Antigen Vaccine System. Pharmaceutics. 2017 Dec 6;9(4):57. doi:10.3390/pharmaceutics9040057. PMID 29215597; PMCID PMC5750663."
  vignette <- "Badhan_2017_liposomal_vaccine"
  units <- list(time = "day", dosing = "ug", concentration = "ug")

  # Empirical compartmental model (Equations 2-3). The central compartment is
  # the IM site of injection (quadriceps); there is no reported volume, so the
  # apparent volume is 1 and Cc is the site-of-injection antigen amount (ug).
  compartmentData <- list(
    central = list(analyte = "Ag85B-ESAT-6 antigen", units = "ug", specimen = "administration site", verified = TRUE),
    peripheral1 = list(analyte = "Ag85B-ESAT-6 antigen", units = "ug", specimen = "not applicable", verified = TRUE)
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
    disease_state = "Healthy mice dosed IM with a dual-radiolabelled DDA:TDB cationic liposome adsorbing the Ag85B-ESAT-6 (H1) TB subunit antigen; site-of-injection antigen kinetics.",
    dose_range = "2 ug antigen IM",
    regions = NA,
    notes = "Empirical model selected by AIC: a 2-compartment model best described the antigen (Table 2; micro-constants k10/k12/k21, distribution half-life 0.0069 days, terminal half-life 2.62 days). The paper also develops a minimal-PBPK model (Equations 7-10); it is described in the vignette but not packaged because it is not reproducible from the printed (dimensionally inconsistent) equations. Companion compartmental model: Badhan_2017_liposome_mouse_1cmt."
  )

  ini({
    lk10 <- log(0.34); label("log central elimination rate constant k10 (day^-1)")  # Table 2 (Antigen): k10 = 0.34 day^-1
    lk12 <- log(22.26); label("log central-to-peripheral transfer rate constant k12 (day^-1)")  # Table 2 (Antigen): k12 = 22.26 day^-1
    lk21 <- log(77.58); label("log peripheral-to-central transfer rate constant k21 (day^-1)")  # Table 2 (Antigen): k21 = 77.58 day^-1
  })

  model({
    k10 <- exp(lk10)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    # Equations 2-3: dC1/dt = -(k10+k12)*C1 + k21*C2; dC2/dt = k12*C1 - k21*C2.
    # Apparent volume 1; Cc is the central (site-of-injection) amount.
    d/dt(central) <- -(k10 + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    Cc <- central
  })
}
