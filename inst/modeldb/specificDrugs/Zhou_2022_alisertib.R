Zhou_2022_alisertib <- function() {
  description <- paste0(
    "Paediatric two-compartment population PK model for the investigational ",
    "Aurora A kinase inhibitor alisertib (MLN8237) given orally to children ",
    "and adolescents aged 2 to 21 years with relapsed or refractory solid ",
    "tumours or acute leukaemias (Zhou 2022, n = 146 from the Children's ",
    "Oncology Group trials ADVL0812 and ADVL0921), with a THREE-TRANSIT-",
    "COMPARTMENT absorption chain sharing a single rate constant ktr and ",
    "linear elimination. Body surface area is a power covariate on the ",
    "apparent clearance CL/F (exponent 0.742) and the apparent central ",
    "volume V1/F (exponent 1.47); the reference BSA is not printed in the ",
    "source and is taken as the 1.25 m^2 cohort mean (see the vignette for ",
    "the back-solve that supports it). The enteric-coated tablet has a ",
    "relative bioavailability of 0.671 versus the powder-in-capsule ",
    "reference formulation. Doses are supplied in mg into transit1; the ",
    "observation Cc is in nmol/L, converted in model() with the paper's ",
    "molar mass of 518.92 g/mol. ",
    "Companion exposure-safety logistic models are ",
    "Zhou_2022_alisertib_stomatitis and Zhou_2022_alisertib_febrile_neutropenia."
  )
  reference <- paste(
    "Zhou X, Mould DR, Yuan Y, Fox E, Greengard E, Faller DV,",
    "Venkatakrishnan K. Population pharmacokinetics and exposure-safety",
    "relationships of alisertib in children and adolescents with advanced",
    "malignancies. J Clin Pharmacol. 2022;62(2):206-219.",
    "doi:10.1002/jcph.1958. Final parameter estimates from Table 2; model",
    "diagram in Supplemental Figure S1.",
    sep = " "
  )
  vignette <- "Zhou_2022_alisertib"
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")
  # The oral dose enters transit1, not a state named depot or central.
  dosing <- "transit1"

  covariateData <- list(
    BSA = list(
      description = "Baseline body surface area.",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on CL/F (exponent 0.742) and V1/F (exponent 1.47)",
        "(Table 2, BSACL and BSAV1). Neither the functional form nor the",
        "reference value is printed; the power form is the Mould-lab",
        "convention of the companion adult model (Zhou 2018) and the",
        "reference is taken as 1.25 m^2, the Table 1 cohort mean (SD 0.46).",
        "Body weight fitted marginally better by objective function, but BSA",
        "was retained because paediatric alisertib doses are BSA-based. The",
        "BSA formula is not stated in the source."
      ),
      source_name = "BSA"
    ),
    FORM_ALISERTIB_ECT = list(
      description = "Alisertib enteric-coated tablet formulation indicator (versus powder-in-capsule).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (powder-in-capsule)",
      notes = paste(
        "1 = enteric-coated tablet (all 100 ADVL0921 patients), 0 =",
        "powder-in-capsule (all 46 ADVL0812 patients). The tablet's relative",
        "bioavailability is 0.671 (Table 2, FECT 67.1 percent). Each",
        "formulation was used in only one study, so the effect is confounded",
        "with study design and sampling scheme, as the Discussion cautions."
      ),
      source_name = "Formulation"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight.",
      units = "kg",
      type = "continuous",
      notes = "Screened; the most significant body-size predictor by objective function, but BSA gave comparable fits and was carried forward instead (Results). Mean 41.9 kg (SD 24.55).",
      source_name = "Weight"
    ),
    AGE = list(
      description = "Baseline age.",
      units = "year",
      type = "continuous",
      notes = "Screened and not retained; after BSA normalisation, CL/F shows no residual age dependence (Figure 2). Mean 11.1 years (range 2 to 21).",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained. 67 of 146 (46 percent) female.",
      source_name = "Sex"
    ),
    RACE_BLACK = list(
      description = "Black race indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Race was screened and not retained (89 White, 24 Black, 6 Asian, 27 other or unknown).",
      source_name = "Race"
    ),
    RACE_HISPANIC = list(
      description = "Hispanic or Latino ethnicity indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Ethnicity was screened and not retained.",
      source_name = "Ethnicity"
    )
  )

  compartmentData <- list(
    transit1 = list(analyte = "alisertib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "alisertib", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "alisertib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "alisertib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "alisertib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 146L,
    n_studies = 2L,
    n_observations = 606L,
    age_range = "2 to 21 years (77 aged 2-11, 40 aged 12-16, 29 aged 17-21)",
    age_mean = "11.1 years",
    weight_mean = "41.9 kg (SD 24.55)",
    bsa_mean = "1.25 m^2 (SD 0.46)",
    sex_female_pct = 46,
    race_ethnicity = c(White = 61, Black = 16, Asian = 4, `Other/unknown` = 18),
    disease_state = "Relapsed or refractory solid tumours, neuroblastoma (ADVL0812) or solid tumours and acute leukaemias (ADVL0921; 19 haematological malignancies)",
    dose_range = "45 to 100 mg/m^2 once daily or 30 to 40 mg/m^2 twice daily as powder-in-capsule (ADVL0812, 25 to 150 mg); 80 mg/m^2 once daily as enteric-coated tablet (ADVL0921, 40 to 160 mg); days 1-7 of 21-day cycles",
    regions = "United States and Canada (Children's Oncology Group sites)",
    notes = paste0(
      "Baseline characteristics are Table 1 of Zhou 2022. ADVL0812 ",
      "(NCT02444884, n = 46) sampled richly to 24 h on cycle 1 day 1; ",
      "ADVL0921 (NCT01154816, n = 100) sampled sparsely to 8 h. Concentrations ",
      "were measured by LC-MS/MS (LLOQ 10 nmol/L). 80 mg/m^2 once daily as the ",
      "enteric-coated tablet is the paediatric MTD and recommended phase 2 dose."
    )
  )

  ini({
    # Zhou 2022 Table 2, 'Final Model Parameters for the Pediatric
    # Population PK Model'. NONMEM 7.3. Apparent parameters refer to the
    # powder-in-capsule formulation (Figure 2 caption; Table S1 footnote).
    lcl <- log(1.84); label("Apparent oral clearance CL/F at the reference BSA (L/h)") # Table 2, CL/F 1.84 L/h, SE 12.0%
    lvc <- log(24.1); label("Apparent central volume V1/F at the reference BSA (L)") # Table 2, V1/F 24.1 L, SE 14.4%
    lq <- log(2.66); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2, Q/F 2.66 L/h, SE 19.2%
    lvp <- log(32.3); label("Apparent peripheral volume V2/F (L)") # Table 2, V2/F 32.3 L, SE 20.2%
    lktr <- log(2.35); label("Transit-compartment rate constant Ktr (1/h)") # Table 2, Ktr 2.35 1/h, SE 9.7%

    e_form_ect_fdepot <- 0.671; label("Relative bioavailability of the enteric-coated tablet versus powder-in-capsule (fraction)") # Table 2, FECT 67.1%, SE 12.5%
    e_bsa_cl <- 0.742; label("Power exponent of BSA on CL/F (unitless)") # Table 2, BSACL 0.742, SE 19.5%
    e_bsa_vc <- 1.47; label("Power exponent of BSA on V1/F (unitless)") # Table 2, BSAV1 1.47, SE 22.2%

    # Table 2 'BSV (ratio)' is read as the log-scale SD, as in the
    # authors' companion adult model (Zhou 2018, which reads 0.518 back
    # as 'CV 51.8%'); variances are its square. CL/F-V1/F correlation
    # 0.583 (Table 2 footnote a): cov = 0.583 * 0.581 * 0.699 = 0.23677.
    # Q/F carries no BSV.
    etalcl + etalvc ~ c(0.337561, 0.236770, 0.488601) # Table 2, BSV 0.581 (CL/F) and 0.699 (V1/F), correlation 0.583
    etalvp ~ 0.868624 # Table 2, BSV 0.932 (V2/F)
    etalktr ~ 0.291600 # Table 2, BSV 0.540 (Ktr)

    propSd <- 0.59; label("Proportional residual error (fraction)") # Table 2, CCV 0.59, SE 6.8%
  })

  model({
    # BSA reference 1.25 m^2 (Table 1 cohort mean); not printed in the source.
    cl <- exp(lcl + etalcl) * (BSA / 1.25)^e_bsa_cl
    vc <- exp(lvc + etalvc) * (BSA / 1.25)^e_bsa_vc
    q <- exp(lq)
    vp <- exp(lvp + etalvp)
    ktr <- exp(lktr + etalktr)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Supplemental Figure S1: oral dose -> three transit compartments at
    # the common rate Ktr -> central; two-compartment disposition.
    d/dt(transit1) <- -ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(central) <- ktr * transit3 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Enteric-coated tablet relative bioavailability (Table 2, FECT).
    f(transit1) <- e_form_ect_fdepot^FORM_ALISERTIB_ECT

    # mg/L -> nmol/L with the molar mass 518.92 g/mol printed in Methods.
    Cc <- central / vc / 518.92 * 1e6
    Cc ~ prop(propSd)
  })
}
