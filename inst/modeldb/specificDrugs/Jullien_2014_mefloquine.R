Jullien_2014_mefloquine <- function() {
  description <- paste(
    "Population PK model for oral mefloquine given as the artesunate-mefloquine",
    "fixed-dose combination (200 mg mefloquine base + 100 mg artesunate per",
    "tablet, two tablets once daily for 3 days) in Indian adults with acute",
    "uncomplicated Plasmodium falciparum malaria (Jullien 2014). Two-compartment",
    "disposition with first-order absorption and first-order elimination from",
    "the central compartment; body weight enters Vc/F and Vp/F as power",
    "functions centred on the cohort median weight. No covariate on CL/F. No",
    "IIV on Q/F. Proportional residual error on plasma concentrations.",
    sep = " "
  )
  reference <- paste(
    "Jullien V, Valecha N, Srivastava B, Sharma B, Kiechel JR (2014).",
    "Population pharmacokinetics of mefloquine, administered as a fixed-dose",
    "combination of artesunate-mefloquine in Indian patients for the treatment",
    "of acute uncomplicated Plasmodium falciparum malaria.",
    "Malaria Journal 13:187. doi:10.1186/1475-2875-13-187.",
    sep = " "
  )
  vignette <- "Jullien_2014_mefloquine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "mefloquine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mefloquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mefloquine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight at the first drug intake",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on Vc/F and Vp/F, normalised by the cohort median",
        "body weight (Methods, generalized allometric model equation). The",
        "median is not printed; 53 kg is the median of the individual body",
        "weights plotted in Figure 4 (cohort mean 53.3 kg, SD 7.3, range",
        "40-70 kg, Table 1)."
      ),
      source_name = "BW"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 77L,
    n_studies = 1L,
    n_observations = 355L,
    age_range = "18-45 years (mean 28.1, SD 9.0; Table 1)",
    weight_range = "40-70 kg (mean 53.3, SD 7.3; Table 1)",
    sex_female_pct = 3.9,
    race_ethnicity = "Indian",
    disease_state = "Acute uncomplicated Plasmodium falciparum malaria (asexual parasitaemia 1,000-100,000 /uL and fever >= 37.5 C)",
    dose_range = "Mefloquine 400 mg base once daily for 3 days (2 tablets of artesunate-mefloquine 100/200 mg FDC per day; 1200 mg total)",
    regions = "India (multicentre: National Institute of Malaria Research New Delhi, Goa, Mangalore)",
    parasitaemia = "Median 8356 /uL, range 1165-94693 (Table 1)",
    notes = paste(
      "74 men and 3 women (Results, Table 1). Single-arm open-label trial",
      "December 2007 - November 2008 (ISRCTN70618692). Plasma mefloquine by",
      "LC-MS/MS (LLOQ 9.9 ng/mL). Samples pre-dose D0, D3 (72 h), D7, and one",
      "random time on day 28, 35 or 42. NONMEM VI, FOCE-I."
    )
  )

  ini({
    lka <- log(0.163); label("Absorption rate constant ka (1/h)") # Table 2 final model Ka = 0.163 h-1 (SD 0.0263)
    lcl <- log(1.13); label("Apparent clearance CL/F (L/h)") # Table 2 final model CL/F = 1.13 L/h (SD 0.0501)
    lvc <- log(271); label("Apparent central volume Vc/F at the median body weight (L)") # Table 2 final model Vc/F = 271 L (SD 14.1)
    lq <- log(1.43); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 final model Q/F = 1.43 L/h (SD 0.152)
    lvp <- log(344); label("Apparent peripheral volume Vp/F at the median body weight (L)") # Table 2 final model Vp/F = 344 L (SD 40.6)

    e_wt_vc <- 0.87; label("Power exponent of body weight on Vc/F (unitless)") # Table 2 final model theta BW,Vc = 0.87 (SD 0.288)
    e_wt_vp <- 2.41; label("Power exponent of body weight on Vp/F (unitless)") # Table 2 final model theta BW,Vp = 2.41 (SD 1.08)

    etalka ~ 0.562 # Table 2 final model omega2 ka = 0.562 (SD 0.232)
    etalcl ~ 0.0894 # Table 2 final model omega2 CL/F = 0.0894 (SD 0.0218)
    etalvc ~ 0.0453 # Table 2 final model omega2 Vc/F = 0.0453 (SD 0.0163)
    etalvp ~ 0.295 # Table 2 final model omega2 Vp/F = 0.295 (SD 0.128)

    propSd <- 0.2631; label("Proportional residual error (SD, fraction)") # Table 2 final model sigma2 = 0.0692 -> SD = sqrt(0.0692) = 0.2631
  })

  model({
    # Generalized allometric model (Methods): P = TV(P) * (BW / median BW)^theta,
    # centred on the cohort median body weight of 53 kg.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (WT / 53)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (WT / 53)^e_wt_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = mg/L, the plasma concentration units of Figures 2-3.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
