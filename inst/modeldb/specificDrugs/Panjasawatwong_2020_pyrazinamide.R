Panjasawatwong_2020_pyrazinamide <- function() {
  description <- "One-compartment population PK model for oral pyrazinamide in plasma and cerebrospinal fluid (CSF) in Vietnamese children with tuberculous meningitis (Panjasawatwong 2020). Absorption through a depot plus three transit compartments sharing one rate constant ktr = 4/MTT; a CSF compartment of age-dependent physiological volume exchanges with the central compartment through QCSF, with unbound (fu = 0.9 fixed) drug entering at QCSF * fu * PC. Clearance and QCSF carry fixed allometric weight scaling (exponent 0.75) and Vc exponent 1, referenced to 10.9 kg; clearance also carries a sigmoidal postmenstrual-age maturation function normalised to a 3-year-old; weight-for-age z-score has linear effects on CL/F (+4.76% per unit) and Vc/F (-4.65% per unit). IIV on CL/F, Vc/F and MTT, interoccasion variability on F across the four sampling days; additive-on-log-scale residual error for plasma and CSF."
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
      notes = "Drives two terms. (1) Maturation of CL/F through postmenstrual age PMA (months) = AGE * 12 + 9.33 (Methods: 'assuming full-term gestation at 9.33 months'), MF = PMA^HILL / (PMA^HILL + MAT50^HILL) (Eq 4), normalised to the 3-year-old typical patient (PMA 45.33 months; Table 4 footnote b). (2) The age-dependent CSF volume of Eq 1, age in years.",
      source_name = "AGE"
    ),
    WAZ = list(
      description = "Weight-for-age z-score",
      units = "unitless (z-score)",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effects: CL/F * (1 + 0.0476 * WAZ) and Vc/F * (1 - 0.0465 * WAZ) (Table 4 'WAZ on CL/F (%)' = 4.76, 'WAZ on Vc/F (%)' = -4.65; Discussion: 'clearance decreased 4.76% per unit of WAZ decrease, while the central volume of distribution increased 4.65% per unit of WAZ decrease'). The paper does not print a centring value; the z-score's own reference (0, the growth-standard median) is used, so the Table 4 typical CL/F and Vc/F apply at WAZ = 0. Cohort median -1.93, range -5.52 to 1.99 (Table 1). The growth reference used to compute WAZ is not stated.",
      source_name = "WAZ"
    ),
    OCC = list(
      description = "Sampling occasion index (1 = day 1, 2 = day 14, 3 = day 30, 4 = day 90)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Selects the interoccasion eta on relative bioavailability (Methods: IOV 'between the four different sampling occasions (days 1, 14, 30, and 90)'). Time-varying within a subject. Values outside 1-4 switch IOV off (typical F).",
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "pyrazinamide", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "pyrazinamide", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "pyrazinamide", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "pyrazinamide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pyrazinamide", units = "mg", specimen = "plasma", verified = TRUE),
    csf = list(analyte = "pyrazinamide", units = "mg", specimen = "CSF", verified = TRUE)
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
    disease_state = "Suspected tuberculous meningitis (TBM); severity grade I 58%, II 24%, III 18%; HIV positive 4%, negative 92%, unknown 4%. Weight-for-age z-score median -1.93 (range -5.52 to 1.99).",
    dose_range = "Pyrazinamide 25 mg/kg orally once daily for the first 3 months (WHO 2006 paediatric regimen), with isoniazid, rifampicin, ethambutol, streptomycin (first 2 months) and adjunctive dexamethasone.",
    regions = "Vietnam (Pham Ngoc Thach Hospital, Ho Chi Minh City), October 2009 to March 2011",
    notes = "Demographics from Panjasawatwong 2020 Table 1. 519 plasma and 155 CSF pyrazinamide concentrations; 11 plasma samples (2.12%) below the 800 ug/L LOQ were coded as missing. Plasma sampled on days 1, 14, 30 and 90; CSF on days 30 and 90 within 15 min of a plasma sample."
  )

  ini({
    # Structural parameters: Table 4 'Population estimate', scaled to the
    # typical patient (10.9 kg, 3 years) per Table 4 footnote b.
    lcl <- log(1.07); label("Apparent oral clearance CL/F at 10.9 kg, 3 years and WAZ 0 (L/h)") # Table 4 CL/F = 1.07 L/h (RSE 4.02%)
    lvc <- log(7.38); label("Apparent central volume Vc/F at 10.9 kg and WAZ 0 (L)") # Table 4 Vc/F = 7.38 L (RSE 2.73%)
    lmtt <- log(0.457); label("Mean transit absorption time MTT (h)") # Table 4 MTT = 0.457 h (RSE 10.8%)
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability F (unitless)") # Table 4 F = 1 (fixed)

    # Postmenstrual-age maturation of clearance (Methods Eq 4)
    ltm50_cl <- log(12.1); label("Postmenstrual age at 50% maturation of clearance MAT50 (months)") # Table 4 MAT50 = 12.1 months (RSE 7.26%)
    lhill_cl <- log(2.73); label("Hill coefficient of the clearance maturation function (unitless)") # Table 4 HILL = 2.73 (RSE 29.0%)

    # CSF distribution
    lqcsf <- log(0.0964); label("Apparent intercompartmental clearance between central and CSF QCSF/F at 10.9 kg (L/h)") # Table 4 QCSF/F = 0.0964 L/h (RSE 13.6%)
    fu <- fixed(0.9); label("Fraction unbound in plasma (unitless)") # Table 4 fu = 0.9 (fixed); Results cite refs 25-26
    lpc <- log(1.02); label("Blood-brain penetration multiplier PC on the unbound central-to-CSF flux (unitless)") # Table 4 PC = 1.02 (RSE 1.70%)

    # Covariate effects
    e_waz_cl <- 0.0476; label("Fractional change in CL/F per unit weight-for-age z-score (unitless)") # Table 4 'WAZ on CL/F (%)' = 4.76 (RSE 18.3%)
    e_waz_vc <- -0.0465; label("Fractional change in Vc/F per unit weight-for-age z-score (unitless)") # Table 4 'WAZ on Vc/F (%)' = -4.65 (RSE 27.6%)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F and QCSF/F (unitless)") # Methods: 'exponent fixed to 0.75'
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vc/F (unitless)") # Methods: 'exponent fixed to 1'

    # IIV: Table 4 reports %CV; omega^2 = log(1 + CV^2)
    etalcl ~ 0.039221 # Table 4 IIV CL/F = 20.0% CV; log(1 + 0.200^2) = 0.039221
    etalvc ~ 0.032942 # Table 4 IIV Vc/F = 18.3% CV; log(1 + 0.183^2) = 0.032942
    etalmtt ~ 0.347846 # Table 4 IIV MTT = 64.5% CV; log(1 + 0.645^2) = 0.347846

    # IOV on F across the four sampling occasions, one shared variance
    # (nlmixr2 has no SAME keyword, so occasions 2-4 repeat it with fix())
    etaiov_fdepot_1 ~ 0.035464 # Table 4 IOV F = 19.0% CV; log(1 + 0.190^2) = 0.035464
    etaiov_fdepot_2 ~ fix(0.035464) # same IOV variance, occasion 2 (day 14)
    etaiov_fdepot_3 ~ fix(0.035464) # same IOV variance, occasion 3 (day 30)
    etaiov_fdepot_4 ~ fix(0.035464) # same IOV variance, occasion 4 (day 90)

    # Residual error: additive on log-transformed concentrations. Table 4
    # footnote a gives sigma as the VARIANCE; SD = sqrt(variance).
    expSd <- 0.16553; label("Residual error SD on the log scale, plasma (unitless)") # Table 4 sigma Plasma = 0.0274 (variance); sqrt(0.0274) = 0.16553
    expSd_Ccsf <- 0.10677; label("Residual error SD on the log scale, CSF (unitless)") # Table 4 sigma CSF = 0.0114 (variance); sqrt(0.0114) = 0.10677
  })

  model({
    # 1. Derived covariate terms
    iov_f <- (OCC == 1) * etaiov_fdepot_1 + (OCC == 2) * etaiov_fdepot_2 +
      (OCC == 3) * etaiov_fdepot_3 + (OCC == 4) * etaiov_fdepot_4
    # Postmenstrual age in months (Methods: full-term gestation 9.33 months)
    pma <- AGE * 12 + 9.33
    hill_cl <- exp(lhill_cl)
    tm50_cl <- exp(ltm50_cl)
    mf <- pma^hill_cl / (pma^hill_cl + tm50_cl^hill_cl)
    # Typical patient of Table 4 footnote b is 3 years old (PMA 45.33 months)
    mf_ref <- 45.33^hill_cl / (45.33^hill_cl + tm50_cl^hill_cl)
    # CSF volume (L) from Methods Eq 1, read as a sigmoid in percent of the
    # 150-mL adult volume (the equation as typeset has misplaced brackets and
    # no /100; see the vignette)
    agep <- AGE^1.071
    vcsf <- 0.150 * (38.78 + agep * (102.6 - 38.78) / (agep + 1.297^1.071)) / 100

    # 2. Individual parameters
    cl <- exp(lcl + etalcl) * (WT / 10.9)^e_wt_cl * (mf / mf_ref) * (1 + e_waz_cl * WAZ)
    vc <- exp(lvc + etalvc) * (WT / 10.9)^e_wt_vc * (1 + e_waz_vc * WAZ)
    qcsf <- exp(lqcsf) * (WT / 10.9)^e_wt_cl
    pc <- exp(lpc)
    mtt <- exp(lmtt + etalmtt)
    # Depot plus three transit compartments, ka = ktr
    ktr <- 4 / mtt

    # 3. ODE system (Fig 2 structure: central-to-CSF flux QCSF * fu * PC / Vc,
    # return QCSF / VCSF)
    Cc <- central / vc
    Ccsf <- csf / vcsf
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(central) <- ktr * transit3 - cl * Cc - qcsf * fu * pc * Cc + qcsf * Ccsf
    d/dt(csf) <- qcsf * fu * pc * Cc - qcsf * Ccsf

    # 4. Bioavailability with interoccasion variability
    f(depot) <- exp(lfdepot + iov_f)

    # 5. Observations: additive error on log-transformed data
    Cc ~ lnorm(expSd)
    Ccsf ~ lnorm(expSd_Ccsf)
  })
}
