Panjasawatwong_2020_rifampicin <- function() {
  description <- "One-compartment population PK model for oral rifampicin in plasma and cerebrospinal fluid (CSF) in Vietnamese children with tuberculous meningitis (Panjasawatwong 2020). Absorption through a depot plus two transit compartments sharing one rate constant ktr = 3/MTT; autoinduction via an enzyme-turnover pool (kENZ, Emax and EC50 fixed to the Smythe 2012 adult values) that multiplies the pre-induced clearance; a CSF compartment of age-dependent physiological volume exchanges with central through QCSF, with unbound (fu = 0.2 fixed) drug entering at QCSF * fu * PC and PC rising exponentially with CSF protein. Clearance and QCSF carry fixed allometric weight scaling (exponent 0.75), Vc exponent 1, referenced to 10.9 kg; clearance also carries a sigmoidal postmenstrual-age maturation function normalised to a 3-year-old. IIV on CL/F, Vc/F, MTT and PC; additive-on-log-scale residual error for plasma and CSF."
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
      notes = "Fixed allometric scaling (WT/10.9)^0.75 on CL/F and QCSF/F and (WT/10.9)^1 on Vc/F (Methods 'Population pharmacokinetic analysis'); 10.9 kg is the cohort median (Table 1). The CSF volume is NOT weight-scaled.",
      source_name = "WT"
    ),
    AGE = list(
      description = "Postnatal age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives two terms. (1) Maturation of CL/F through postmenstrual age PMA (months) = AGE * 12 + 9.33 (Methods: 'assuming full-term gestation at 9.33 months'), MF = PMA^HILL / (PMA^HILL + MAT50^HILL) (Eq 4), normalised to the 3-year-old typical patient (PMA 45.33 months; Table 3 footnote b). (2) The age-dependent CSF volume of Eq 1, age in years.",
      source_name = "AGE"
    ),
    CSF_TPRO = list(
      description = "CSF total protein concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Exponential effect on the blood-brain penetration multiplier: PC = 0.844 * exp(0.245 * CSF_TPRO) ('an increase of 1 g/liter in CSF protein concentration resulted in a 1.28-fold increase in PC'; Table 3 'CSF protein on PC (%)' = 24.5). The paper does not print the centring value; the uncentred form (reference 0 g/L) is the one that reproduces Figure 4's simulated CSF exposures at 0.2, 1.0 and 5.0 g/L (see the vignette). Cohort median 1.20 g/L, range 0.1-5 (Table 1).",
      source_name = "CSF protein"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rifampicin", units = "mg", specimen = "plasma", verified = TRUE),
    csf = list(analyte = "rifampicin", units = "mg", specimen = "CSF", verified = TRUE),
    enz_pool = list(
      analyte = "metabolising enzyme (relative amount, pre-induced = 1)",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
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
    disease_state = "Suspected tuberculous meningitis (TBM); severity grade I 58%, II 24%, III 18%; HIV positive 4%, negative 92%, unknown 4%. Baseline CSF protein median 1.20 g/L (range 0.1-5).",
    dose_range = "Rifampicin 10 mg/kg orally once daily for 8 months (WHO 2006 paediatric regimen), with isoniazid, pyrazinamide, ethambutol, streptomycin (first 2 months) and adjunctive dexamethasone.",
    regions = "Vietnam (Pham Ngoc Thach Hospital, Ho Chi Minh City), October 2009 to March 2011",
    notes = "Demographics from Panjasawatwong 2020 Table 1. 512 plasma and 155 CSF rifampicin concentrations; 26 plasma samples (5.08%) below the 8.0 ug/L LOQ were coded as missing. Plasma sampled on days 1, 14, 30 and 90; CSF on days 30 and 90 within 15 min of a plasma sample."
  )

  ini({
    # Structural parameters: Table 3 'Population estimate', scaled to the
    # typical patient (10.9 kg, 3 years) per Table 3 footnote b.
    lcl <- log(3.22); label("Apparent oral pre-induced clearance CL/F at 10.9 kg and 3 years (L/h)") # Table 3 CL/F = 3.22 L/h (RSE 8.94%)
    lvc <- log(12.3); label("Apparent central volume Vc/F at 10.9 kg (L)") # Table 3 Vc/F = 12.3 L (RSE 6.46%)
    lmtt <- log(1.25); label("Mean transit absorption time MTT (h)") # Table 3 MTT = 1.25 h (RSE 13.0%)
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability F (unitless)") # Table 3 F = 1 (fixed)

    # Postmenstrual-age maturation of clearance (Methods Eq 4)
    ltm50_cl <- log(6.81); label("Postmenstrual age at 50% maturation of clearance MAT50 (months)") # Table 3 MAT50 = 6.81 months (RSE 25.0%)
    lhill_cl <- log(1.38); label("Hill coefficient of the clearance maturation function (unitless)") # Table 3 HILL = 1.38 (RSE 33.0%)

    # Enzyme-turnover autoinduction (Methods Eqs 2-3), fixed to Smythe 2012
    lkenz <- fixed(log(0.00369)); label("First-order enzyme degradation rate constant kENZ (1/h)") # Table 3 Kenz = 0.00369 1/h (fixed; ref 6 Smythe 2012)
    lemax <- fixed(log(1.04)); label("Maximum fractional increase in enzyme production Emax (unitless)") # Table 3 Emax = 1.04 (fixed)
    lec50 <- fixed(log(0.0705)); label("Plasma rifampicin concentration giving half of Emax (mg/L)") # Table 3 EC50 = 70.5 ug/L (fixed) = 0.0705 mg/L

    # CSF distribution
    lqcsf <- log(0.00482); label("Apparent intercompartmental clearance between central and CSF QCSF/F at 10.9 kg (L/h)") # Table 3 QCSF/F = 0.00482 L/h (RSE 11.6%)
    fu <- fixed(0.2); label("Fraction unbound in plasma (unitless)") # Table 3 fu = 0.2 (fixed); Results cite refs 25-26
    lpc <- log(0.844); label("Blood-brain penetration multiplier PC on the unbound central-to-CSF flux at CSF protein 0 g/L (unitless)") # Table 3 PC = 0.844 (RSE 7.84%)
    e_csf_tpro_pc <- 0.245; label("Exponential coefficient of CSF protein on PC (per g/L)") # Table 3 'CSF protein on PC (%)' = 24.5 (RSE 17.6%); exp(0.245) = 1.28-fold per g/L (Results)

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F and QCSF/F (unitless)") # Methods: 'exponent fixed to 0.75'
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vc/F (unitless)") # Methods: 'exponent fixed to 1'

    # IIV: Table 3 reports %CV; omega^2 = log(1 + CV^2)
    etalcl ~ 0.036946 # Table 3 IIV CL/F = 19.4% CV; log(1 + 0.194^2) = 0.036946
    etalvc ~ 0.051549 # Table 3 IIV Vc/F = 23.0% CV; log(1 + 0.230^2) = 0.051549
    etalmtt ~ 0.543777 # Table 3 IIV MTT = 85.0% CV; log(1 + 0.850^2) = 0.543777
    etalpc ~ 0.047266 # Table 3 IIV PC = 22.0% CV; log(1 + 0.220^2) = 0.047266

    # Residual error: additive on log-transformed concentrations. Table 3
    # footnote a gives sigma as the VARIANCE; SD = sqrt(variance).
    expSd <- 0.71624; label("Residual error SD on the log scale, plasma (unitless)") # Table 3 sigma Plasma = 0.513 (variance); sqrt(0.513) = 0.71624
    expSd_Ccsf <- 0.55588; label("Residual error SD on the log scale, CSF (unitless)") # Table 3 sigma CSF = 0.309 (variance); sqrt(0.309) = 0.55588
  })

  model({
    # 1. Derived covariate terms
    # Postmenstrual age in months (Methods: full-term gestation 9.33 months)
    pma <- AGE * 12 + 9.33
    hill_cl <- exp(lhill_cl)
    tm50_cl <- exp(ltm50_cl)
    mf <- pma^hill_cl / (pma^hill_cl + tm50_cl^hill_cl)
    # Typical patient of Table 3 footnote b is 3 years old (PMA 45.33 months)
    mf_ref <- 45.33^hill_cl / (45.33^hill_cl + tm50_cl^hill_cl)
    # CSF volume (L) from Methods Eq 1, read as a sigmoid in percent of the
    # 150-mL adult volume (the equation as typeset has misplaced brackets and
    # no /100; see the vignette)
    agep <- AGE^1.071
    vcsf <- 0.150 * (38.78 + agep * (102.6 - 38.78) / (agep + 1.297^1.071)) / 100

    # 2. Individual parameters
    cl_pre <- exp(lcl + etalcl) * (WT / 10.9)^e_wt_cl * (mf / mf_ref)
    vc <- exp(lvc + etalvc) * (WT / 10.9)^e_wt_vc
    qcsf <- exp(lqcsf) * (WT / 10.9)^e_wt_cl
    pc <- exp(lpc + etalpc) * exp(e_csf_tpro_pc * CSF_TPRO)
    mtt <- exp(lmtt + etalmtt)
    # Depot plus two transit compartments, ka = ktr (Fig 2)
    ktr <- 3 / mtt
    kenz <- exp(lkenz)
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    # 3. ODE system (Fig 2). The enzyme pool starts at 1 and its zero-order
    # production rate equals kENZ, so it stays at 1 without drug (Eq 2);
    # clearance is CLpre times the relative enzyme amount (Eq 3).
    Cc <- central / vc
    Ccsf <- csf / vcsf
    cl <- cl_pre * enz_pool
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(central) <- ktr * transit2 - cl * Cc - qcsf * fu * pc * Cc + qcsf * Ccsf
    d/dt(csf) <- qcsf * fu * pc * Cc - qcsf * Ccsf
    d/dt(enz_pool) <- kenz * (1 + emax * Cc / (ec50 + Cc)) - kenz * enz_pool
    enz_pool(0) <- 1

    # 4. Bioavailability
    f(depot) <- exp(lfdepot)

    # 5. Observations: additive error on log-transformed data
    Cc ~ lnorm(expSd)
    Ccsf ~ lnorm(expSd_Ccsf)
  })
}
