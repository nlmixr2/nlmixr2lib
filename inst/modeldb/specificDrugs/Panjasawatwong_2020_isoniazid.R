Panjasawatwong_2020_isoniazid <- function() {
  description <- "Two-compartment population PK model for oral isoniazid in plasma and cerebrospinal fluid (CSF) in Vietnamese children with tuberculous meningitis (Panjasawatwong 2020). Absorption through a depot plus two transit compartments sharing one rate constant ktr = 3/MTT; a CSF compartment of age-dependent physiological volume exchanges with the central compartment through QCSF, with unbound (fu = 0.9 fixed) drug entering at QCSF * fu * PC. Clearance and intercompartmental flows carry fixed allometric weight scaling (exponent 0.75) and volumes exponent 1, referenced to 10.9 kg; clearance also carries a sigmoidal postmenstrual-age maturation function normalised to a 3-year-old, and is 56.4% lower in NAT2 slow acetylators. IIV on CL/F and Q/F; additive-on-log-scale residual error for plasma and CSF."
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
      notes = "Fixed allometric scaling (WT/10.9)^0.75 on CL/F, Q/F and QCSF/F and (WT/10.9)^1 on Vc/F and Vp/F (Methods 'Population pharmacokinetic analysis'); 10.9 kg is the cohort median (Table 1). The CSF volume is NOT weight-scaled.",
      source_name = "WT"
    ),
    AGE = list(
      description = "Postnatal age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives two terms. (1) Maturation of CL/F through postmenstrual age PMA (months) = AGE * 12 + 9.33 (Methods: 'assuming full-term gestation at 9.33 months'), MF = PMA^HILL / (PMA^HILL + MAT50^HILL) (Eq 4), normalised to the 3-year-old typical patient (PMA 45.33 months; Table 2 footnote b). (2) The age-dependent CSF volume of Eq 1, age in years.",
      source_name = "AGE"
    ),
    NAT2_SLOW = list(
      description = "NAT2 slow-acetylator phenotype indicator (1 = slow, 0 = fast or intermediate)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fast plus intermediate acetylators, pooled)",
      notes = "The final model uses two groups: 'fast' = genotype-predicted fast + intermediate + unknown (n = 72) and slow (n = 28). CL/F is multiplied by (1 - 0.564) in slow acetylators (Table 2 'Slow acetylators (%)' = 56.4). Children without a NAT2 genotype were assigned to the intermediate (NAT2_SLOW = 0) group.",
      source_name = "NAT2 acetylator status"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "isoniazid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "isoniazid", units = "mg", specimen = "plasma", verified = TRUE),
    csf = list(analyte = "isoniazid", units = "mg", specimen = "CSF", verified = TRUE)
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
    dose_range = "Isoniazid 5 mg/kg orally once daily for 8 months (WHO 2006 paediatric regimen), with rifampicin, pyrazinamide, ethambutol, streptomycin (first 2 months) and adjunctive dexamethasone.",
    regions = "Vietnam (Pham Ngoc Thach Hospital, Ho Chi Minh City), October 2009 to March 2011",
    nat2_phenotype = "Genotype-predicted fast 17, intermediate 47, slow 28, unknown 8",
    notes = "Demographics from Panjasawatwong 2020 Table 1. 523 plasma and 140 CSF isoniazid concentrations; 48 plasma samples (9.2%) below the 12 ug/L LOQ were handled with the M3 method. Plasma sampled on days 1, 14, 30 and 90; CSF on days 30 and 90 within 15 min of a plasma sample."
  )

  ini({
    # Structural parameters: Table 2 'Population estimate', scaled to the
    # typical patient (10.9 kg, 3 years) per Table 2 footnote b.
    lcl <- log(9.43); label("Apparent oral clearance CL/F in fast/intermediate acetylators at 10.9 kg and 3 years (L/h)") # Table 2 CL/F = 9.43 L/h (RSE 6.55%)
    lvc <- log(3.78); label("Apparent central volume Vc/F at 10.9 kg (L)") # Table 2 Vc/F = 3.78 L (RSE 33.6%)
    lq <- log(28.0); label("Apparent intercompartmental clearance Q/F at 10.9 kg (L/h)") # Table 2 Q/F = 28.0 L/h (RSE 21.5%)
    lvp <- log(15.3); label("Apparent peripheral volume Vp/F at 10.9 kg (L)") # Table 2 Vp/F = 15.3 L (RSE 9.60%)
    lmtt <- log(0.878); label("Mean transit absorption time MTT (h)") # Table 2 MTT = 0.878 h (RSE 10.9%)
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability F (unitless)") # Table 2 F = 1 (fixed)

    # Postmenstrual-age maturation of clearance (Methods Eq 4)
    ltm50_cl <- log(12.7); label("Postmenstrual age at 50% maturation of clearance MAT50 (months)") # Table 2 MAT50 = 12.7 months (RSE 6.85%)
    lhill_cl <- log(4.7); label("Hill coefficient of the clearance maturation function (unitless)") # Table 2 HILL = 4.7 (RSE 13.5%)

    # CSF distribution
    lqcsf <- log(13.7); label("Apparent intercompartmental clearance between central and CSF QCSF/F at 10.9 kg (L/h)") # Table 2 QCSF/F = 13.7 L/h (RSE 43.6%)
    fu <- fixed(0.9); label("Fraction unbound in plasma (unitless)") # Table 2 fu = 0.9 (fixed); Results cite literature value (ref 24)
    lpc <- log(1.65); label("Blood-brain penetration multiplier PC on the unbound central-to-CSF flux (unitless)") # Table 2 PC = 1.65 (RSE 5.24%)

    # Covariate effects
    e_nat2_slow_cl <- 0.564; label("Fractional reduction of CL/F in NAT2 slow acetylators (unitless)") # Table 2 'Slow acetylators (%)' = 56.4 (RSE 6.61%)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F, Q/F and QCSF/F (unitless)") # Methods: 'exponent fixed to 0.75'
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vc/F and Vp/F (unitless)") # Methods: 'exponent fixed to 1'

    # IIV: Table 2 reports %CV; omega^2 = log(1 + CV^2)
    etalcl ~ 0.12715 # Table 2 IIV CL/F = 36.8% CV; log(1 + 0.368^2) = 0.12715
    etalq ~ 0.70315 # Table 2 IIV Q/F = 101% CV; log(1 + 1.01^2) = 0.70315

    # Residual error: additive on log-transformed concentrations. Table 2
    # footnote a gives sigma as the VARIANCE; SD = sqrt(variance).
    expSd <- 0.68848; label("Residual error SD on the log scale, plasma (unitless)") # Table 2 sigma Plasma = 0.474 (variance); sqrt(0.474) = 0.68848
    expSd_Ccsf <- 0.41231; label("Residual error SD on the log scale, CSF (unitless)") # Table 2 sigma CSF = 0.170 (variance); sqrt(0.170) = 0.41231
  })

  model({
    # 1. Derived covariate terms
    # Postmenstrual age in months (Methods: full-term gestation 9.33 months)
    pma <- AGE * 12 + 9.33
    hill_cl <- exp(lhill_cl)
    tm50_cl <- exp(ltm50_cl)
    mf <- pma^hill_cl / (pma^hill_cl + tm50_cl^hill_cl)
    # Typical patient of Table 2 footnote b is 3 years old (PMA 45.33 months)
    mf_ref <- 45.33^hill_cl / (45.33^hill_cl + tm50_cl^hill_cl)
    # CSF volume (L) from Methods Eq 1, read as a sigmoid in percent of the
    # 150-mL adult volume (the equation as typeset has misplaced brackets and
    # no /100; see the vignette)
    agep <- AGE^1.071
    vcsf <- 0.150 * (38.78 + agep * (102.6 - 38.78) / (agep + 1.297^1.071)) / 100

    # 2. Individual parameters
    cl <- exp(lcl + etalcl) * (WT / 10.9)^e_wt_cl * (mf / mf_ref) *
      (1 - e_nat2_slow_cl * NAT2_SLOW)
    vc <- exp(lvc) * (WT / 10.9)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 10.9)^e_wt_cl
    vp <- exp(lvp) * (WT / 10.9)^e_wt_vc
    qcsf <- exp(lqcsf) * (WT / 10.9)^e_wt_cl
    pc <- exp(lpc)
    mtt <- exp(lmtt)
    # Depot plus two transit compartments, ka = ktr (Fig 2 structure)
    ktr <- 3 / mtt

    # 3. ODE system (Fig 2: central-to-CSF flux QCSF * fu * PC / Vc, return QCSF / VCSF)
    Cc <- central / vc
    Ccsf <- csf / vcsf
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(central) <- ktr * transit2 - cl * Cc - q * Cc + q * peripheral1 / vp -
      qcsf * fu * pc * Cc + qcsf * Ccsf
    d/dt(peripheral1) <- q * Cc - q * peripheral1 / vp
    d/dt(csf) <- qcsf * fu * pc * Cc - qcsf * Ccsf

    # 4. Bioavailability
    f(depot) <- exp(lfdepot)

    # 5. Observations: additive error on log-transformed data
    Cc ~ lnorm(expSd)
    Ccsf ~ lnorm(expSd_Ccsf)
  })
}
