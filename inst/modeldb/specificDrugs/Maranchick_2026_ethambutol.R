Maranchick_2026_ethambutol <- function() {
  description <- "Two-compartment population pharmacokinetic model with first-order absorption, an absorption lag time and linear elimination for oral ethambutol in Ghanaian children with tuberculosis with or without HIV coinfection (Maranchick 2026); estimated allometric weight scaling on CL/F (exponent 0.70) and V1/F (exponent 0.62) normalised to 15.1 kg, and an exponential HIV-positive effect raising CL/F by 24.6%."
  reference <- "Maranchick NF, Martyn-Dickens C, Enimil A, Yang H, Amissah AK, Dompreh A, Bosomtwe D, Sly-Moore E, Opoku T, Appiah AF, Asiedu P, Antwi S, Scheetz MH, Peloquin CA, Kwara A. Population pharmacokinetics of pyrazinamide and ethambutol in children with tuberculosis with or without HIV. Antimicrob Agents Chemother. 2026. doi:10.1128/aac.00909-25"
  vignette <- "Maranchick_2026_pyrazinamide_ethambutol"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Maranchick 2026 Materials and Methods
  # ('Study design': oral dosing; 'blood samples ... plasma aliquoted' and
  # quantified by LC-MS/MS).
  compartmentData <- list(
    depot       = list(analyte = "ethambutol", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "ethambutol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ethambutol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with estimated exponents on CL/F (0.70) and V1/F (0.62), normalised to a 15.1 kg reference weight, per Maranchick 2026 Results 'Ethambutol' paragraph 1: 'Cl/F: Cl/F*(Weight/15.1)^0.7 and V1/F: V1/F*(Weight/15.1)^0.62'. The paper reports no weight exponent on Q/F or V2/F, so those two parameters are left unscaled here; see the vignette Errata. The 15.1 kg normalisation constant is the value printed by the paper; it is close to but not identical to the cohort median weight of 16 kg (Table 1). Cohort weight range 4-60 kg.",
      source_name        = "Weight"
    ),
    HIV_POS = list(
      description        = "HIV coinfection status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HIV-negative; children with TB alone)",
      notes              = "1 = TB/HIV coinfected, 0 = TB only. Time-fixed per subject. Applied as an exponential effect on apparent clearance: cl = cl_typ * exp(e_hiv_pos_cl * HIV_POS) with e_hiv_pos_cl = 0.22 (Maranchick 2026 Table 2, row 'Exponent, HIV+ on Cl/F'). exp(0.22) - 1 = 0.246, matching the Results statement that children with TB/HIV had clearance 25% faster than HIV-negative children. Adding HIV status dropped the objective function by 17 points. HIV medications were tested separately and did not enhance the model; lopinavir/ritonavir, which affected ethambutol in Tikiso et al., was used by only four participants here.",
      source_name        = "HIV+"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 84L,
    n_studies        = 1L,
    age_range        = "0.3-14.5 years (median 5.0); 49.4% under 5 years, 18.8% under 2 years",
    age_median       = "5.0 years",
    weight_range     = "4-60 kg (median 16)",
    weight_median    = "16 kg",
    sex_female_pct   = 38.8,
    hiv_positive_pct = 51.8,
    disease_state    = "Children with drug-susceptible tuberculosis, with TB alone or with TB/HIV coinfection. 24 (28.2%) of the enrolled 85 were malnourished (body-mass-index-for-age Z score below -2 SD). Of the TB/HIV participants, 65.9% received efavirenz-based antiretroviral therapy.",
    dose_range       = "Ethambutol 20 mg/kg once daily (WHO target range 15-25 mg/kg); administered median 21.4 mg/kg (range 14.3-34.2) as part of the first-line HRZE regimen. Children under 25 kg received single dispersible EMB 100 mg tablets alongside dispersible HRZ; children 25 kg and over received adult HRZE 75/150/400/275 mg tablets.",
    regions          = "Ghana (Komfo Anokye Teaching Hospital, Kumasi).",
    notes            = "Two-arm PK study, enrolment February 2019 to June 2021, children 3 months to 14 years. PK sampling on one occasion after at least 4 weeks of HRZE therapy (steady state), with blood drawn at 0 (pre-dose), 1, 2, 4, 8 and 12 h post-dose after an overnight fast. 501 samples from 84 participants entered the final EMB model: one TB/HIV participant was removed because all samples were below or near the limit of quantification (suggesting a missed dose or post-dose vomiting), and six further samples were removed as apparently mislabelled out of order. LC-MS/MS quantification range 0.05-10 mg/L. Fitted in Monolix2024R1 by SAEM; a lognormal distribution was assumed for all parameters. A maturation function on clearance, inter-occasion variability and HIV-medication covariates were each tested and did not improve model fit. Demographics from Table 1 (whole cohort, n = 85); parameter estimates from Table 2."
  )

  ini({
    # Structural PK parameters -- Maranchick 2026 Table 2, ethambutol rows.
    # Typical values are at the 15.1 kg reference weight in an HIV-negative child.
    ltlag <- log(0.67);   label("Absorption lag time (h)")                                     # Table 2 EMB 'tlag (h)' = 0.67 (RSE 7.04)
    lka   <- log(3.83);   label("First-order absorption rate constant (1/h)")                  # Table 2 EMB 'Ka (h-1)' = 3.83 (RSE 22.97)
    lvc   <- log(95.16);  label("Apparent central volume V1/F at WT = 15.1 kg (L)")            # Table 2 EMB 'V1/F (L)' = 95.16 (RSE 6.49)
    lcl   <- log(23.2);   label("Apparent oral clearance CL/F at WT = 15.1 kg, HIV-negative (L/h)")  # Table 2 EMB 'Cl/F (L/h)' = 23.2 (RSE 4.4)
    lq    <- log(11.25);  label("Apparent inter-compartmental clearance Q/F (L/h)")            # Table 2 EMB 'Q/F (L/h)' = 11.25 (RSE 6.69)
    lvp   <- log(162.41); label("Apparent peripheral volume V2/F (L)")                         # Table 2 EMB 'V2/F (L)' = 162.41 (RSE 16.45)

    # Allometric exponents on body weight -- estimated (the paper reports an RSE
    # for each), not fixed at the canonical 0.75 / 1. Results 'Ethambutol':
    # 'Weight was allometrically scaled on Cl/F and V1/F using estimated
    # coefficients (Cl/F: Cl/F*(Weight/15.1)^0.7 and V1/F: V1/F*(Weight/15.1)^0.62)'.
    # No weight exponent is reported for Q/F or V2/F, so neither is scaled.
    e_wt_cl <- 0.70; label("Allometric exponent on CL/F (unitless)")                           # Table 2 EMB 'Exponent, BWonCl/F' = 0.70 (RSE 10.31)
    e_wt_vc <- 0.62; label("Allometric exponent on V1/F (unitless)")                           # Table 2 EMB 'Exponent, BWonV1/F' = 0.62 (RSE 19.34)

    # HIV coinfection effect on clearance, exponential on the log-parameter scale:
    # cl = cl_typ * exp(e_hiv_pos_cl * HIV_POS). exp(0.22) - 1 = 0.246, reproducing
    # the Results statement of 25% faster clearance in children with TB/HIV.
    e_hiv_pos_cl <- 0.22; label("Exponential effect of HIV-positive status on CL/F (unitless)")  # Table 2 EMB 'Exponent, HIV+ on Cl/F' = 0.22 (RSE 21.03)

    # Inter-individual variability. Table 2 reports the SD of the random effects
    # first and the corresponding coefficient of variation in parentheses (per the
    # Table 2 footnote: 'IIV, interindividual variability (reported as standard
    # deviation of the random effects); CV, coefficient of variability'). nlmixr2
    # takes variance, so each entry below is the printed SD squared. The printed
    # CV values reproduce as sqrt(exp(omega^2) - 1) * 100, confirming the log scale.
    etaltlag ~ 0.1225   # Table 2 EMB tlag: SD 0.35 (36.35 CV)  -- 0.35^2 = 0.1225
    etalka   ~ 1.3924   # Table 2 EMB Ka:   SD 1.18 (172.61 CV) -- 1.18^2 = 1.3924
    etalvc   ~ 0.2601   # Table 2 EMB V1/F: SD 0.51 (54.63 CV)  -- 0.51^2 = 0.2601
    etalcl   ~ 0.1024   # Table 2 EMB Cl/F: SD 0.32 (33.26 CV)  -- 0.32^2 = 0.1024
    etalq    ~ 0.1156   # Table 2 EMB Q/F:  SD 0.34 (35.23 CV)  -- 0.34^2 = 0.1156
    etalvp   ~ 0.6561   # Table 2 EMB V2/F: SD 0.81 (96.75 CV)  -- 0.81^2 = 0.6561

    # Combined additive-plus-proportional residual error (Monolix 'a' and 'b'),
    # Results 'Ethambutol': 'A combination of additive and proportional error
    # models was utilized.'
    addSd  <- 0.03; label("Additive residual SD (ug/mL)")            # Table 2 EMB 'Residual variability' a = 0.03 (RSE 9.99)
    propSd <- 0.17; label("Proportional residual SD (fraction)")     # Table 2 EMB 'Residual variability' b = 0.17 (RSE 7.08)
  })

  model({
    # Individual PK parameters. Allometric weight scaling normalised to 15.1 kg on
    # CL/F and V1/F only; HIV coinfection accelerates apparent clearance.
    tlag <- exp(ltlag + etaltlag)
    ka   <- exp(lka + etalka)
    vc   <- exp(lvc + etalvc) * (WT / 15.1)^e_wt_vc
    cl   <- exp(lcl + etalcl) * (WT / 15.1)^e_wt_cl * exp(e_hiv_pos_cl * HIV_POS)
    q    <- exp(lq + etalq)
    vp   <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Dose mg / (V1/F in L) -> mg/L, which is numerically ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
