Panjasawatwong_2020_ethambutol <- function() {
  description <- "Two-compartment population PK model for oral ethambutol in plasma in Vietnamese children with tuberculous meningitis (Panjasawatwong 2020). Absorption through a depot plus two transit compartments sharing one rate constant ktr = 3/MTT. Clearance and Q/F carry fixed allometric weight scaling (exponent 0.75) and volumes exponent 1, referenced to 10.9 kg; clearance also carries a postmenstrual-age maturation function (Hill fixed to 1) normalised to a 3-year-old. IIV on relative bioavailability, Vc/F, MTT and Vp/F; additive-on-log-scale residual error. No CSF model (CSF ethambutol could not be quantified)."
  reference <- paste(
    "Panjasawatwong N, Wattanakul T, Hoglund RM, Bang ND, Pouplin T,",
    "Nosoongnoen W, Ngo VN, Day JN, Tarning J. (2020).",
    "Population pharmacokinetic properties of antituberculosis drugs in",
    "Vietnamese children with tuberculous meningitis.",
    "Antimicrob Agents Chemother 65(1):e00487-20.",
    "doi:10.1128/AAC.00487-20.",
    sep = " "
  )
  vignette <- "Panjasawatwong_2020_antituberculosis_tbm"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Fixed allometric scaling (WT/10.9)^0.75 on CL/F and Q/F and (WT/10.9)^1 on Vc/F and Vp/F (Methods 'Population pharmacokinetic analysis'); 10.9 kg is the cohort median (Table 1).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Postnatal age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Maturation of CL/F through postmenstrual age PMA (months) = AGE * 12 + 9.33 (Methods: 'assuming full-term gestation at 9.33 months'), MF = PMA^HILL / (PMA^HILL + MAT50^HILL) (Eq 4) with HILL fixed to 1, normalised to the 3-year-old typical patient (PMA 45.33 months; Table 5 footnote b).",
      source_name = "AGE"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "ethambutol", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "ethambutol", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "ethambutol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ethambutol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ethambutol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 100L,
    n_studies = 1L,
    age_range = "2 months to 15 years (0.167-15.0 years)",
    age_median = "3.0 years",
    weight_range = "4.0-43 kg",
    weight_median = "10.9 kg",
    sex_female_pct = 44,
    race_ethnicity = "Vietnamese",
    disease_state = "Suspected tuberculous meningitis (TBM); severity grade I 58%, II 24%, III 18%; HIV positive 4%, negative 92%, unknown 4%.",
    dose_range = "Ethambutol 15 mg/kg orally once daily for 8 months (WHO 2006 paediatric regimen), with isoniazid, rifampicin, pyrazinamide, streptomycin (first 2 months) and adjunctive dexamethasone.",
    regions = "Vietnam (Pham Ngoc Thach Hospital, Ho Chi Minh City), October 2009 to March 2011",
    notes = "Demographics from Panjasawatwong 2020 Table 1. 517 plasma ethambutol concentrations; 1 sample (0.193%) below the 8 ug/L LOQ was coded as missing. Plasma sampled on days 1, 14, 30 and 90."
  )

  ini({
    # Structural parameters: Table 5 'Population estimate', scaled to the
    # typical patient (10.9 kg, 3 years) per Table 5 footnote b.
    lcl <- log(28.2); label("Apparent oral clearance CL/F at 10.9 kg and 3 years (L/h)") # Table 5 CL/F = 28.2 L/h (RSE 4.68%)
    lvc <- log(98.6); label("Apparent central volume Vc/F at 10.9 kg (L)") # Table 5 Vc/F = 98.6 L (RSE 11.0%)
    lq <- log(16.9); label("Apparent intercompartmental clearance Q/F at 10.9 kg (L/h)") # Table 5 Q/F = 16.9 L/h (RSE 12.2%)
    lvp <- log(153); label("Apparent peripheral volume Vp/F at 10.9 kg (L)") # Table 5 Vp/F = 153 L (RSE 16.0%)
    lmtt <- log(1.8); label("Mean transit absorption time MTT (h)") # Table 5 MTT = 1.8 h (RSE 7.15%)
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability F (unitless)") # Table 5 F = 1 (fixed)

    # Postmenstrual-age maturation of clearance (Methods Eq 4)
    ltm50_cl <- log(3.99); label("Postmenstrual age at 50% maturation of clearance MAT50 (months)") # Table 5 MAT50 = 3.99 months (RSE 36.0%)
    lhill_cl <- fixed(log(1)); label("Hill coefficient of the clearance maturation function (unitless)") # Table 5 HILL = 1 (fixed); Results: estimate 76.1 was implausible

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F and Q/F (unitless)") # Methods: 'exponent fixed to 0.75'
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vc/F and Vp/F (unitless)") # Methods: 'exponent fixed to 1'

    # IIV: Table 5 reports %CV; omega^2 = log(1 + CV^2)
    etalfdepot ~ 0.039607 # Table 5 IIV F = 20.1% CV; log(1 + 0.201^2) = 0.039607
    etalvc ~ 0.232001 # Table 5 IIV Vc/F = 51.1% CV; log(1 + 0.511^2) = 0.232001
    etalmtt ~ 0.033296 # Table 5 IIV MTT = 18.4% CV; log(1 + 0.184^2) = 0.033296
    etalvp ~ 0.656160 # Table 5 IIV Vp/F = 96.3% CV; log(1 + 0.963^2) = 0.656160

    # Residual error: additive on log-transformed concentrations. Table 5
    # footnote a gives sigma as the VARIANCE; SD = sqrt(variance).
    expSd <- 0.44385; label("Residual error SD on the log scale, plasma (unitless)") # Table 5 sigma Plasma = 0.197 (variance); sqrt(0.197) = 0.44385
  })

  model({
    # 1. Derived covariate terms
    # Postmenstrual age in months (Methods: full-term gestation 9.33 months)
    pma <- AGE * 12 + 9.33
    hill_cl <- exp(lhill_cl)
    tm50_cl <- exp(ltm50_cl)
    mf <- pma^hill_cl / (pma^hill_cl + tm50_cl^hill_cl)
    # Typical patient of Table 5 footnote b is 3 years old (PMA 45.33 months)
    mf_ref <- 45.33^hill_cl / (45.33^hill_cl + tm50_cl^hill_cl)

    # 2. Individual parameters
    cl <- exp(lcl) * (WT / 10.9)^e_wt_cl * (mf / mf_ref)
    vc <- exp(lvc + etalvc) * (WT / 10.9)^e_wt_vc
    q <- exp(lq) * (WT / 10.9)^e_wt_cl
    vp <- exp(lvp + etalvp) * (WT / 10.9)^e_wt_vc
    mtt <- exp(lmtt + etalmtt)
    # Depot plus two transit compartments, ka = ktr
    ktr <- 3 / mtt

    # 3. ODE system
    Cc <- central / vc
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(central) <- ktr * transit2 - cl * Cc - q * Cc + q * peripheral1 / vp
    d/dt(peripheral1) <- q * Cc - q * peripheral1 / vp

    # 4. Bioavailability with IIV (population F fixed to 1)
    f(depot) <- exp(lfdepot + etalfdepot)

    # 5. Observation: additive error on log-transformed data
    Cc ~ lnorm(expSd)
  })
}
