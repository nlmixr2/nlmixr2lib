Maranchick_2026_pyrazinamide <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption, an absorption lag time and linear elimination for oral pyrazinamide in Ghanaian children with tuberculosis with or without HIV coinfection (Maranchick 2026); estimated allometric weight scaling on CL/F (exponent 0.70) and V/F (exponent 0.79) normalised to 15 kg, and an exponential HIV-positive effect raising CL/F by 18.5%."
  reference <- "Maranchick NF, Martyn-Dickens C, Enimil A, Yang H, Amissah AK, Dompreh A, Bosomtwe D, Sly-Moore E, Opoku T, Appiah AF, Asiedu P, Antwi S, Scheetz MH, Peloquin CA, Kwara A. Population pharmacokinetics of pyrazinamide and ethambutol in children with tuberculosis with or without HIV. Antimicrob Agents Chemother. 2026. doi:10.1128/aac.00909-25"
  vignette <- "Maranchick_2026_pyrazinamide_ethambutol"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Maranchick 2026 Materials and Methods
  # ('Study design': oral dosing; 'blood samples ... plasma aliquoted' and
  # quantified by LC-MS/MS).
  compartmentData <- list(
    depot   = list(analyte = "pyrazinamide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pyrazinamide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with estimated exponents on CL/F (0.70) and V/F (0.79), normalised to a 15 kg reference weight, per Maranchick 2026 Results 'Pyrazinamide' paragraph 1: 'Cl/F*(Weight/15)^0.7' and 'V/F*(Weight/15)^0.79'. The 15 kg normalisation constant is the value printed by the paper; it is close to but not identical to the cohort median weight of 16 kg (Table 1). Cohort weight range 4-60 kg.",
      source_name        = "Weight"
    ),
    HIV_POS = list(
      description        = "HIV coinfection status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HIV-negative; children with TB alone)",
      notes              = "1 = TB/HIV coinfected (n = 44), 0 = TB only (n = 41). Time-fixed per subject. Applied as an exponential effect on apparent clearance: cl = cl_typ * exp(e_hiv_pos_cl * HIV_POS) with e_hiv_pos_cl = 0.17 (Maranchick 2026 Table 2, row 'Exponent, HIV+ on Cl/F'). exp(0.17) - 1 = 0.185, matching the Results statement that children with TB/HIV had clearance 18.5% faster than children with TB alone. HIV medications (abacavir, dolutegravir, efavirenz, lopinavir/ritonavir) were tested separately as covariates but did not improve model fit.",
      source_name        = "HIV+"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 85L,
    n_studies        = 1L,
    age_range        = "0.3-14.5 years (median 5.0); 49.4% under 5 years, 18.8% under 2 years, 7.1% under 1 year",
    age_median       = "5.0 years",
    weight_range     = "4-60 kg (median 16)",
    weight_median    = "16 kg",
    sex_female_pct   = 38.8,
    hiv_positive_pct = 51.8,
    disease_state    = "Children with drug-susceptible tuberculosis, 41 (48.2%) with TB alone and 44 (51.8%) with TB/HIV coinfection. 24 (28.2%) malnourished (body-mass-index-for-age Z score below -2 SD). Of the TB/HIV participants, 29 (65.9%) received efavirenz-based antiretroviral therapy, 7 (15.9%) abacavir-based, 4 (9.1%) dolutegravir-based and 4 (9.1%) lopinavir/ritonavir-based.",
    dose_range       = "Pyrazinamide 35 mg/kg once daily (WHO target range 30-40 mg/kg); administered median 31.6 mg/kg (range 21.4-49.7) as part of the first-line HRZE regimen. Children under 25 kg received dispersible HRZ 50/75/150 mg tablets (1, 2, 3 or 4 tablets in the 4-<8, 8-<12, 12-<16 and 16-<25 kg weight bands); children 25 kg and over received adult HRZE 75/150/400/275 mg tablets.",
    regions          = "Ghana (Komfo Anokye Teaching Hospital, Kumasi).",
    notes            = "Two-arm PK study, enrolment February 2019 to June 2021, children 3 months to 14 years. PK sampling on one occasion after at least 4 weeks of HRZE therapy (steady state), with blood drawn at 0 (pre-dose), 1, 2, 4, 8 and 12 h post-dose after an overnight fast. 509 samples from 85 participants entered the final PZA model; 16 samples were below the limit of quantification (14 reported as 0 mg/L, 2 interval-censored). LC-MS/MS quantification range 0.5-100 mg/L. Fitted in Monolix2024R1 by SAEM; a lognormal distribution was assumed for all parameters. A maturation function on clearance, inter-occasion variability and HIV-medication covariates were each tested and did not improve model fit. Demographics from Table 1; parameter estimates from Table 2."
  )

  ini({
    # Structural PK parameters -- Maranchick 2026 Table 2, pyrazinamide rows.
    # Typical values are at the 15 kg reference weight in an HIV-negative child.
    ltlag <- log(0.26);  label("Absorption lag time (h)")                                      # Table 2 PZA 'tlag (h)' = 0.26 (RSE 46.86)
    lka   <- log(3.76);  label("First-order absorption rate constant (1/h)")                   # Table 2 PZA 'Ka (h-1)' = 3.76 (RSE 20.79)
    lvc   <- log(11.30); label("Apparent central volume V/F at WT = 15 kg (L)")                # Table 2 PZA 'V/F (L)' = 11.30 (RSE 2.87)
    lcl   <- log(1.27);  label("Apparent oral clearance CL/F at WT = 15 kg, HIV-negative (L/h)")  # Table 2 PZA 'Cl/F (L/h)' = 1.27 (RSE 4.46)

    # Allometric exponents on body weight -- estimated (the paper reports an RSE
    # for each), not fixed at the canonical 0.75 / 1. Results 'Pyrazinamide':
    # 'Weight was allometrically scaled using estimated coefficients on Cl/F
    # (Cl/F*(Weight/15)^0.7) and V/F (V/F*(Weight/15)^0.79)'.
    e_wt_cl <- 0.70; label("Allometric exponent on CL/F (unitless)")                           # Table 2 PZA 'Exponent, BWonCl/F' = 0.70 (RSE 9.74)
    e_wt_vc <- 0.79; label("Allometric exponent on V/F (unitless)")                            # Table 2 PZA 'Exponent, BWonV/F' = 0.79 (RSE 6.64)

    # HIV coinfection effect on clearance, exponential on the log-parameter scale:
    # cl = cl_typ * exp(e_hiv_pos_cl * HIV_POS). exp(0.17) - 1 = 0.185, reproducing
    # the Results statement of 18.5% faster clearance in children with TB/HIV.
    e_hiv_pos_cl <- 0.17; label("Exponential effect of HIV-positive status on CL/F (unitless)")  # Table 2 PZA 'Exponent, HIV+ on Cl/F' = 0.17 (RSE 31.29)

    # Inter-individual variability. Table 2 reports the SD of the random effects
    # first and the corresponding coefficient of variation in parentheses (per the
    # Table 2 footnote: 'IIV, interindividual variability (reported as standard
    # deviation of the random effects); CV, coefficient of variability'). nlmixr2
    # takes variance, so each entry below is the printed SD squared. The printed
    # CV values reproduce as sqrt(exp(omega^2) - 1) * 100, confirming the log scale.
    etaltlag ~ 0.5625   # Table 2 PZA tlag: SD 0.75 (87.16 CV) -- 0.75^2 = 0.5625
    etalka   ~ 0.3844   # Table 2 PZA Ka:   SD 0.62 (68.71 CV) -- 0.62^2 = 0.3844
    etalvc   ~ 0.0529   # Table 2 PZA V/F:  SD 0.23 (23.4 CV)  -- 0.23^2 = 0.0529
    etalcl   ~ 0.1024   # Table 2 PZA Cl/F: SD 0.32 (32.6 CV)  -- 0.32^2 = 0.1024

    # Combined additive-plus-proportional residual error (Monolix 'a' and 'b').
    addSd  <- 1.06; label("Additive residual SD (ug/mL)")            # Table 2 PZA 'Residual variability' a = 1.06 (RSE 11.62)
    propSd <- 0.06; label("Proportional residual SD (fraction)")     # Table 2 PZA 'Residual variability' b = 0.06 (RSE 13.33)
  })

  model({
    # Individual PK parameters. Allometric weight scaling normalised to 15 kg;
    # HIV coinfection accelerates apparent clearance.
    tlag <- exp(ltlag + etaltlag)
    ka   <- exp(lka + etalka)
    vc   <- exp(lvc + etalvc) * (WT / 15)^e_wt_vc
    cl   <- exp(lcl + etalcl) * (WT / 15)^e_wt_cl * exp(e_hiv_pos_cl * HIV_POS)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    alag(depot) <- tlag

    # Dose mg / (V/F in L) -> mg/L, which is numerically ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
