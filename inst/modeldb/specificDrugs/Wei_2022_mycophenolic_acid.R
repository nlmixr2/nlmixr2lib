Wei_2022_mycophenolic_acid <- function() {
  description <- "Two-compartment population PK model with first-order absorption and first-order elimination for mycophenolic acid (MPA) after oral mycophenolate mofetil dispersible tablets (MMFdt) in Chinese pediatric patients early after liver transplantation (Wei 2022). Body weight enters every structural parameter except Vp/F by fixed allometric scaling to a 7.5 kg reference (exponent 0.75 on CL/F and Q/F, 1 on Vc/F, -0.25 on Ka); the per-administration MMF dose in mg/kg enters CL/F as a power term normalised to 11.16 mg/kg. Vp/F is fixed at 269 L. Log-normal inter-individual variability on CL/F and Q/F only; exponential residual error."
  reference <- paste(
    "Wei Y, Wu D, Chen Y, Dong C, Qi J, Wu Y, Cai R, Zhou S, Li C, Niu L, Wu T, Xiao Y, Liu T.",
    "Population pharmacokinetics of mycophenolate mofetil in pediatric patients",
    "early after liver transplantation.",
    "Front Pharmacol. 2022;13:1002628.",
    "doi:10.3389/fphar.2022.1002628.",
    sep = " "
  )
  vignette <- "Wei_2022_mycophenolic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "mycophenolate mofetil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mycophenolic acid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mycophenolic acid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at the pharmacokinetic sampling occasion.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Fixed allometric scaling to a 7.5 kg reference (the cohort median, Wei 2022 Table 1: 7.5 kg, IQR 6.0-10.0, range 4.6-27.0): exponent 0.75 on CL/F and Q/F, 1 on Vc/F, -0.25 on Ka (Wei 2022 Equations 1-3 and 5). Vp/F is not weight-scaled (Equation 4).",
      source_name = "WT"
    ),
    DOSE_MMF_MGKG = list(
      description = "Mycophenolate mofetil dose per administration normalised to body weight.",
      units = "mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL/F as (DOSE_MMF_MGKG / 11.16)^0.452 (Wei 2022 Equation 1, Table 3 theta_DOSE). The paper calls the covariate 'DOSE, the MMFdt administered dose' without stating units; the 11.16 normaliser matches the Table 1 cohort median of 11.2 mg/kg/dose (IQR 10.0-15.0, range 8.9-61.5), so the covariate is read as the MMF (not MPA-equivalent) dose per administration in mg/kg. A total-mg reading is excluded: a 7.5 kg child on 11.2 mg/kg receives about 84 mg, an order of magnitude above 11.16. Set it to amt / WT on every record of a q12h regimen.",
      source_name = "DOSE"
    )
  )

  covariatesDataExcluded <- list(
    ALT = list(
      description = "Alanine aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Significant on CL/F in the univariate forward step (Wei 2022 Table 2, Model 3, dOFV 4.819) but not carried into the final model."
    ),
    GRWR = list(
      description = "Graft-to-recipient weight ratio.",
      units = "%",
      type = "continuous",
      notes = "Added on CL/F in forward steps (Wei 2022 Table 2, Models 6, 9, 10) and dropped at backward elimination (Model 12, dOFV 6.329, p > 0.05)."
    ),
    SNP_SLCO1B1_RS4149056 = list(
      description = "SLCO1B1 521T>C genotype.",
      units = "(genotype)",
      type = "categorical",
      notes = "Added on Q/F in forward steps (Wei 2022 Table 2, Models 5, 8, 10) and dropped at backward elimination (Model 13, dOFV 4.037, p > 0.05)."
    ),
    SNP_UGT1A8_RS1042597 = list(
      description = "UGT1A8 518C>G genotype.",
      units = "(genotype)",
      type = "categorical",
      notes = "Significant on Q/F in the univariate forward step (Wei 2022 Table 2, Model 4) but not retained (Model 7, p > 0.05). UGT1A9 -275T>A, UGT1A9 -2152C>T and UGT2B7 211G>T were genotyped and also showed no significant effect."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    n_observations = "115 MPA plasma concentrations (122 samples, 7 below the detection limit removed; the 21 remaining values below the 0.3 mg/L LLOQ were imputed as 0.15 mg/L by the M5 method).",
    age_range = "0.42-7.76 years",
    age_median = "0.74 years (IQR 0.61-1.69); 15 of 20 younger than 24 months",
    weight_range = "4.6-27.0 kg",
    weight_median = "7.5 kg (IQR 6.0-10.0)",
    height_median = "67.5 cm (IQR 62.2-80.0)",
    bsa_median = "0.39 m^2 (IQR 0.32-0.43)",
    sex_female_pct = 40,
    race_ethnicity = "Chinese (single centre, Nanning, Guangxi).",
    disease_state = "Pediatric first liver transplant recipients (16 liver cirrhosis after Kasai operation, 1 each biliary atresia, liver failure, hepatoblastoma and glycogen storage disease); 11 living and 9 deceased donors. Post-operative day at sampling 12 (IQR 10-14, range 4-39).",
    dose_range = "Oral or nasogastric MMF dispersible tablets q12h, starting 10-15 mg/kg per dose and adjusted clinically; 11.2 mg/kg per dose median (IQR 10.0-15.0, range 8.9-61.5) on the sampling day.",
    regions = "China",
    co_medication = "Tacrolimus and methylprednisolone in all patients (triple regimen). Meropenem 55 %, voriconazole 50 %, furosemide 45 %, lansoprazole 30 %, linezolid 25 %, fluconazole 20 %.",
    notes = "Steady-state sampling from day 4 of MMFdt before and 0.5, 1, 2, 4, 8 and 12 h after the morning dose (Wei 2022 Methods; Table 1 demographics)."
  )

  ini({
    lcl <- log(14.8); label("Apparent clearance CL/F at 7.5 kg and 11.16 mg/kg dose (L/h)") # Equation 1 and Table 3 CL/F = 14.8 L/h
    lka <- log(2.02); label("Absorption rate constant Ka at 7.5 kg (1/h)") # Equation 2 Ka = 2.02 1/h (Table 3 rounds to 2.0)
    lvc <- log(6.01); label("Apparent central volume Vc/F at 7.5 kg (L)") # Equation 3 Vc/F = 6.01 L (Table 3 rounds to 6.0)
    lvp <- fixed(log(269)); label("Apparent peripheral volume Vp/F (L)") # Equation 4 and Table 3 Vp/F = 269 L, fixed
    lq <- log(15.4); label("Apparent intercompartmental clearance Q/F at 7.5 kg (L/h)") # Equation 5 and Table 3 Q/F = 15.4 L/h

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)") # Equation 1 exponent 0.75
    e_wt_q <- fixed(0.75); label("Allometric exponent of body weight on Q/F (unitless)") # Equation 5 exponent 0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vc/F (unitless)") # Equation 3 (WT/7.5) linear
    e_wt_ka <- fixed(-0.25); label("Allometric exponent of body weight on Ka (unitless)") # Equation 2 exponent -0.25
    e_dose_cl <- 0.452; label("Power exponent of MMF dose per kg on CL/F (unitless)") # Equation 1 and Table 3 theta_DOSE = 0.452

    etalcl ~ 0.06 # Equation 1 exp(eta) variance 0.06; Table 3 IIV CL/F 24.5 percent = sqrt(0.06)
    etalq ~ 1.39 # Equation 5 exp(eta) variance 1.39; Table 3 IIV Q/F 117.9 percent = sqrt(1.39)

    expSd <- 0.503; label("Exponential residual error (log-scale SD)") # Table 3 RV = 50.3 percent; Results: RV represented as exponential
  })

  model({
    wt_ratio <- WT / 7.5

    cl <- exp(lcl + etalcl) * wt_ratio^e_wt_cl * (DOSE_MMF_MGKG / 11.16)^e_dose_cl
    ka <- exp(lka) * wt_ratio^e_wt_ka
    vc <- exp(lvc) * wt_ratio^e_wt_vc
    vp <- exp(lvp)
    q <- exp(lq + etalq) * wt_ratio^e_wt_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Doses are MMF mg as administered: the paper reports no MPA-equivalent
    # conversion, so the apparent CL/F, Vc/F, Q/F and Vp/F carry both the
    # bioavailability and the MMF-to-MPA mass ratio.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
