vanderVeer_2024_gentamicin <- function() {
  description <- "Two-compartment population PK model of intravenous gentamicin in (near-)term neonates with perinatal asphyxia undergoing controlled therapeutic hypothermia, refitted by van der Veer 2024 to the merged model-building plus external-validation PharmaCool data set (original structure Bijleveld 2016): fixed allometric body-weight scaling (0.75 on CL and Q, 1 on Vc and Vp), a gestational-age power effect on CL, a study-day-5 (> 96 h postnatal age) multiplicative increase in CL, and inter-individual variability on the log-scale residual error."
  reference <- paste(
    "van der Veer MAA, de Haan TR, Franken LGW, Groenendaal F, Dijk PH,",
    "de Boode WP, Simons S, Dijkman KP, van Straaten HLM, Rijken M,",
    "Cools F, Nuytemans DHGM, van Kaam AH, Bijleveld YA, Mathot RAA.",
    "Predictive Performance of a Gentamicin Pharmacokinetic Model in Term",
    "Neonates with Perinatal Asphyxia Undergoing Controlled Therapeutic",
    "Hypothermia. Ther Drug Monit. 2024;46(3):376.",
    "doi:10.1097/FTD.0000000000001166.",
    sep = " "
  )
  vignette <- "vanderVeer_2024_gentamicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "gentamicin", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "gentamicin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight (birth weight in the source cohort)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Birth weight median 3.40 kg (range 2.09-5.07) in the model-building",
        "set and 3.17 kg (2.26-4.62) in the external set (Table 1).",
        "Allometric scaling to a 70 kg reference with fixed exponents 0.75",
        "on CL and Q and 1 on Vc and Vp (Table 2 footnote)."
      ),
      source_name = "BW"
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Median 40 weeks (range 36-42) in both the model-building and the",
        "external data sets (Table 1). Enters CL as (GA / GAmedian)^theta_CLGA",
        "(Table 2 footnote) with GAmedian = 40 weeks."
      ),
      source_name = "GA"
    ),
    PNA = list(
      description = "Postnatal age (chronological since birth)",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. Used only to derive the study-day-5 indicator",
        "theta_SD5 ('study day 5 (>96 h postnatal age)', Table 2 footnote):",
        "SD5 = 1 when PNA > 96 h, else 0. Canonical PNA is in months, so the",
        "model converts PNA_hours = PNA * 24 * 30.4375."
      ),
      source_name = "PNA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 86L,
    n_studies = 1L,
    n_observations = 935L,
    ga_range = "36-42 weeks (median 40 in both data sets)",
    weight_range = "2.09-5.07 kg birth weight (medians 3.40 and 3.17 kg in the two data sets)",
    pna_range = "1.5-5.2 days at end of study (medians 4.7 and 3.1 days in the two data sets)",
    sex_female_pct = 43,
    disease_state = paste(
      "(Near-)term neonates (GA > 36 weeks) with perinatal asphyxia meeting",
      "the criteria for controlled therapeutic hypothermia: 72 h hypothermia",
      "(days 2-3), rewarming (day 4) and normothermia (day 5). Congenital",
      "hepatic or renal pathology was an exclusion criterion."
    ),
    dose_range = paste(
      "Gentamicin 30-min IV infusions of 4 mg/kg every 24 h, with",
      "concentration-guided adjustment per local protocol; daily dose",
      "3.0-5.9 mg/kg across the merged cohort (Table 1)."
    ),
    regions = "10 Dutch and 2 Belgian neonatal intensive care units (PharmaCool Study, NTR2529).",
    notes = paste(
      "Merged data set = the 47-neonate / 612-sample model-building set of",
      "Bijleveld 2016 plus a 39-neonate / 323-sample external validation set",
      "(Table 1). Male: 27 of 47 and 22 of 39 (49 of 86); 43% female is",
      "derived. Refit and 1000-sample bootstrap in NONMEM 7.4.2 (FOCE-I)."
    )
  )

  ini({
    # Structural parameters: van der Veer 2024 Table 2, 'Refit: Merged Data
    # Set' column. Typical values are per 70 kg.
    lcl <- log(1.88); label("Clearance (L/h per 70 kg)") # Table 2 refit merged: CL = 1.88 L/h/70 kg (CV 3%)
    lvc <- log(33.2); label("Central volume of distribution (L per 70 kg)") # Table 2 refit merged: Vc = 33.2 L/70 kg (CV 6%)
    lq <- log(1.81); label("Intercompartmental clearance (L/h per 70 kg)") # Table 2 refit merged: Q = 1.81 L/h/70 kg (CV 7%)
    lvp <- log(32.2); label("Peripheral volume of distribution (L per 70 kg)") # Table 2 refit merged: Vp = 32.2 L/70 kg (CV 8%)

    # Allometric exponents, fixed (Table 2 footnote: TVCL = CL x (BW/70)^0.75,
    # TVVc = Vc x (BW/70)^1, TVQ = Q x (BW/70)^0.75, TVVp = Vp x (BW/70)^1).
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)") # Table 2 footnote: (BW/70)^0.75
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on Vc and Vp (unitless)") # Table 2 footnote: (BW/70)^1

    e_ga_cl <- 2.76; label("Power exponent of (GA/40 wk) on CL (unitless)") # Table 2 refit merged: theta_CLGA = 2.76 (CV 15%)
    e_sd5_cl <- 1.36; label("Multiplicative CL factor when PNA > 96 h (unitless)") # Table 2 refit merged: theta_SD5 = 1.36 (CV 4%)

    # IIV reported as %CV; omega^2 = log(1 + CV^2).
    #   CL 24.5% -> 0.0582925; Vc 30.1% -> 0.0867289; Vp 56.4% -> 0.276188.
    # The refit does not report the CL-Vc correlation; r = 0.81 of the
    # original fit (Bijleveld 2016 Results) is carried forward:
    #   cov = 0.81 * sqrt(0.0582925 * 0.0867289) = 0.0575935.
    etalcl + etalvc ~ c(
      0.0582925,
      0.0575935,
      0.0867289
    ) # Table 2 refit merged: IIV CL 24.5%, IIV Vc 30.1%; r = 0.81 from Bijleveld 2016
    etalvp ~ 0.276188 # Table 2 refit merged: IIV Vp 56.4%

    # Residual error additive on log-transformed concentrations (the original
    # Bijleveld 2016 model), i.e. lnorm() with expSd = sigma.
    expSd <- 0.16; label("Lognormal residual SD (additive on log-transformed concentration)") # Table 2 refit merged: additive error 0.16 (CV 7%)
    # omega^2 = log(1 + 0.523^2) = 0.241792
    etaexpSd ~ 0.241792 # Table 2 refit merged: IIV on additive residual error 52.3%
  })

  model({
    # Study-day-5 indicator: PNA > 96 h (Table 2 footnote). Canonical PNA
    # is in months (1 month = 30.4375 days).
    pna_hours <- PNA * 24 * 30.4375
    sd5_ind <- (pna_hours > 96) * 1.0

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q * e_sd5_cl^sd5_ind * (GA / 40)^e_ga_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    # Subject-level residual magnitude: sigma * exp(eta_RUV).
    expSdi <- expSd * exp(etaexpSd)
    Cc ~ lnorm(expSdi)
  })
}
