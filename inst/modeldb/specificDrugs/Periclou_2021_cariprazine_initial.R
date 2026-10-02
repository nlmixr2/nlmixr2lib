Periclou_2021_cariprazine_initial <- function() {
  description <- "Population PK model (final model on the INITIAL dataset, superseded in the same paper by Periclou_2021_cariprazine) for oral cariprazine and its active metabolites desmethyl-cariprazine (DCAR) and didesmethyl-cariprazine (DDCAR) in adults with schizophrenia or bipolar mania: two-compartment cariprazine with zero-order input into a depot followed by first-order absorption; cariprazine elimination forms one-compartment DCAR and DCAR elimination forms one-compartment DDCAR; additive-linear and power ideal-body-weight, age, race (Black, Asian) and sex covariates and a first-dose shift on cariprazine Vc; fitted to samples within 25 h of dosing"
  reference <- "Periclou A, Phillips L, Ghahramani P, Kapas M, Carrothers T, Khariton T. Population Pharmacokinetics of Cariprazine and its Major Metabolites. Eur J Drug Metab Pharmacokinet. 2021;46(1):53-69. doi:10.1007/s13318-020-00650-4 (initial-dataset model: Supplemental Equation Set 2 and Supplemental Table 3)"
  vignette <- "Periclou_2021_cariprazine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    IBW = list(
      description = "Ideal body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred at 64.5 kg (the model-development cohort median, Supplemental Table 2). Enters additively-linearly on cariprazine CL/F and Vc/F, on DCAR CL/F and on DDCAR Vc/F, and as a power of IBW/64.5 on DCAR Vc/F and DDCAR CL/F (Supplemental Equations A1-A6). The paper does not state which ideal-body-weight formula was used. In the updated analysis (Periclou_2021_cariprazine) IBW was replaced by total body weight.",
      source_name = "IBW"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Additive-linear effect (AGE - 40) on DDCAR Vc/F only (Supplemental Equation A6).",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Additive shift of +2.03 L/h on DCAR CL/F only (Supplemental Equation A3). Same orientation as the source 'Female' indicator.",
      source_name = "Female"
    ),
    RACE_BLACK = list(
      description = "Black or African-American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Black)",
      notes = "Additive shifts on cariprazine CL/F (-1.76 L/h), DCAR CL/F (+24.4 L/h), DDCAR CL/F (+4.23 L/h) and DDCAR Vc/F (+1180 L) (Supplemental Equations A1, A3, A5, A6).",
      source_name = "Black"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Additive shift of -4.01 L/h on cariprazine CL/F only (Supplemental Equation A1). The initial dataset did not include the Japanese study A002-A11, so Asian patients here were essentially the Indian-site patients of the phase 2/3 studies.",
      source_name = "Asian"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "cariprazine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cariprazine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cariprazine", units = "mg", specimen = "plasma", verified = TRUE),
    central_dcar = list(analyte = "desmethyl-cariprazine", units = "mg", specimen = "plasma", verified = TRUE),
    central_ddcar = list(analyte = "didesmethyl-cariprazine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 2049L,
    n_studies = 10L,
    age_range = "18-65 years (updated dataset; the initial-dataset demographics are not reported separately)",
    disease_state = "Adults with schizophrenia or with manic/mixed episodes of bipolar I disorder",
    dose_range = "Oral cariprazine 0.5-18 mg once daily (doses of 21 mg and first doses < 1.5 mg excluded)",
    regions = "Multinational (US, Europe, India)",
    notes = "Initial dataset: three phase 1, three phase 2 and six phase 3 studies, of which the long-term open-label RGH-MD-11 and RGH-MD-17 were held out for validation. Samples collected more than 25 h after a dose were excluded, so the model describes the 0-25 h post-dose window. 11,412 cariprazine, 11,189 DCAR and 10,286 DDCAR samples from 2049, 2044 and 2002 patients (Periclou 2021 Section 3.2.1)."
  )

  ini({
    # Cariprazine (Supplemental Table 3; Supplemental Equations A1-A2)
    ld1 <- log(3.14)
    label("Duration of zero-order input of the dose into the depot, DUR (h)") # Suppl Table 3 'DUR 3.14 (2.93, 3.38)'
    lka <- log(0.578)
    label("First-order absorption rate constant, Ka (1/h)") # Suppl Table 3 'Ka 0.578 (0.501, 0.683)'
    lcl <- log(22.8)
    label("Cariprazine apparent clearance CL/F for a non-Black non-Asian patient with IBW 64.5 kg (L/h)") # Suppl Table 3 'CL/F 22.8 (22.3, 23.3)'; Suppl Eq A1
    e_ibw_cl <- 0.183
    label("Linear effect of (IBW - 64.5) on cariprazine CL/F (L/h/kg)") # Suppl Table 3 'Linear effect of IBW (L/h/kg) 0.183'
    e_black_cl <- -1.76
    label("Additive shift on cariprazine CL/F for Black race (L/h)") # Suppl Table 3 'Additional shift in black patients (L/h) -1.76'
    e_asian_cl <- -4.01
    label("Additive shift on cariprazine CL/F for Asian race (L/h)") # Suppl Table 3 'Additional shift in Asian patients (L/h) -4.01'
    lvc <- log(454)
    label("Cariprazine apparent central volume Vc/F for IBW 64.5 kg after the second and later doses (L)") # Suppl Table 3 'VC/F 454 (397, 515)'; Suppl Eq A2
    e_fd_vc <- fixed(1.47)
    label("Proportional shift on cariprazine Vc/F after the first dose (fraction)") # Suppl Table 3 'Proportional shift in Vc for first dose 1.47 FIXED'
    e_ibw_vc <- 8.55
    label("Linear effect of (IBW - 64.5) on cariprazine Vc/F (L/kg)") # Suppl Table 3 'Linear effect of IBW (L/kg) 8.55'
    lq <- log(92.3)
    label("Cariprazine apparent intercompartmental clearance Q/F (L/h)") # Suppl Table 3 'Q/F 92.3 (66.5, 126)'
    lvp <- log(415)
    label("Cariprazine apparent peripheral volume Vp/F (L)") # Suppl Table 3 'VP/F 415 (334, 490)'

    # DCAR (Supplemental Equations A3-A4)
    lcl_dcar <- log(70.9)
    label("DCAR apparent clearance DCL/F for a non-Black male with IBW 64.5 kg (L/h)") # Suppl Table 3 'DCL/F 70.9 (69.0, 72.7)'; Suppl Eq A3
    e_ibw_cl_dcar <- 1.33
    label("Linear effect of (IBW - 64.5) on DCAR CL/F (L/h/kg)") # Suppl Table 3 'Linear effect of IBW (L/h/kg) 1.33'
    e_black_cl_dcar <- 24.4
    label("Additive shift on DCAR CL/F for Black race (L/h)") # Suppl Table 3 'Additional shift in black patients (L/h) 24.4'
    e_sexf_cl_dcar <- 2.03
    label("Additive shift on DCAR CL/F for female sex (L/h)") # Suppl Table 3 'Additional shift in female patients (L/h) 2.03'
    lvc_dcar <- log(176)
    label("DCAR apparent central volume DVC/F for IBW 64.5 kg (L)") # Suppl Table 3 'DVC/F 176 (159, 198)'; Suppl Eq A4
    e_ibw_vc_dcar <- 3.16
    label("Power exponent of IBW/64.5 on DCAR Vc/F (unitless)") # Suppl Table 3 'Power effect of IBW 3.16'

    # DDCAR (Supplemental Equations A5-A6)
    lcl_ddcar <- log(6.74)
    label("DDCAR apparent clearance DDCL/F for a non-Black patient with IBW 64.5 kg (L/h)") # Suppl Table 3 'DDCL/F 6.74 (6.46, 7.04)'; Suppl Eq A5
    e_ibw_cl_ddcar <- 1.12
    label("Power exponent of IBW/64.5 on DDCAR CL/F (unitless)") # Suppl Table 3 'Power effect of IBW 1.12'
    e_black_cl_ddcar <- 4.23
    label("Additive shift on DDCAR CL/F for Black race (L/h)") # Suppl Table 3 'Additional shift in black patients (L/h) 4.23'
    lvc_ddcar <- log(2220)
    label("DDCAR apparent central volume DDVC/F for a non-Black 40-year-old with IBW 64.5 kg (L)") # Suppl Table 3 'DDVC/F 2220 (2120, 2316)'; Suppl Eq A6
    e_age_vc_ddcar <- 27.4
    label("Linear effect of (AGE - 40) on DDCAR Vc/F (L/year)") # Suppl Table 3 'Linear effect of Age (L/y) 27.4'
    e_ibw_vc_ddcar <- 39.7
    label("Linear effect of (IBW - 64.5) on DDCAR Vc/F (L/kg)") # Suppl Table 3 'Linear effect of IBW (L/kg) 39.7'
    e_black_vc_ddcar <- 1180
    label("Additive shift on DDCAR Vc/F for Black race (L)") # Suppl Table 3 'Additional shift in black patients (L) 1180'; Suppl Eq A6 prints it as a multiplier, see model()

    # IIV: Supplemental Table 3 reports %CV; converted to log-normal variance
    # omega^2 = log(1 + CV^2). Off-diagonal elements are not reported.
    etalka ~ 0.65815 # Suppl Table 3 'Ka IIV 96.5 %CV'
    etalcl ~ 0.11061 # Suppl Table 3 'CL/F IIV 34.2 %CV'
    etalvc ~ 0.17697 # Suppl Table 3 'VC/F IIV 44.0 %CV'
    etalcl_dcar ~ 0.18741 # Suppl Table 3 'DCL/F IIV 45.4 %CV'
    etalvc_dcar ~ 0.75311 # Suppl Table 3 'DVC/F IIV 106 %CV'
    etalcl_ddcar ~ 0.38847 # Suppl Table 3 'DDCL/F IIV 68.9 %CV'
    etalvc_ddcar ~ 0.40160 # Suppl Table 3 'DDVC/F IIV 70.3 %CV'

    # Residual error: not reported for the initial-dataset models.
    propSd <- fixed(0)
    label("Proportional residual error on cariprazine (fraction; ZERO - not reported in source)") # not reported
    propSd_dcar <- fixed(0)
    label("Proportional residual error on DCAR (fraction; ZERO - not reported in source)") # not reported
    propSd_ddcar <- fixed(0)
    label("Proportional residual error on DDCAR (fraction; ZERO - not reported in source)") # not reported
  })

  model({
    # First-dose indicator (the source's 'SD'): 1 after the first dose and before
    # the second, 0 from the second dose onward.
    fd <- 1 * (dosenum() <= 1)
    ibwc <- IBW - 64.5
    ibwn <- IBW / 64.5

    # Cariprazine (Suppl Eq A1-A2)
    d1 <- exp(ld1)
    ka <- exp(lka + etalka)
    cl <- (exp(lcl) + e_ibw_cl * ibwc + e_black_cl * RACE_BLACK + e_asian_cl * RACE_ASIAN) *
      exp(etalcl)
    vc <- (exp(lvc) + e_ibw_vc * ibwc) * (1 + e_fd_vc * fd) * exp(etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    # DCAR (Suppl Eq A3-A4)
    cl_dcar <- (exp(lcl_dcar) + e_ibw_cl_dcar * ibwc + e_black_cl_dcar * RACE_BLACK + e_sexf_cl_dcar * SEXF) *
      exp(etalcl_dcar)
    vc_dcar <- exp(lvc_dcar + etalvc_dcar) * ibwn^e_ibw_vc_dcar

    # DDCAR (Suppl Eq A5-A6). Eq A6 is printed as
    # '2220 + 27.4 x (Age - 40) + 39.7 x (IBW - 64.5) x (1180 x Black)', which
    # would make the volume collapse to 2220 L plus an age term in every
    # non-Black patient; Supplemental Table 3 lists 1180 L as an 'Additional
    # shift in black patients (L)', so the Black term is encoded additively.
    cl_ddcar <- (exp(lcl_ddcar) * ibwn^e_ibw_cl_ddcar + e_black_cl_ddcar * RACE_BLACK) * exp(etalcl_ddcar)
    vc_ddcar <- (exp(lvc_ddcar) + e_age_vc_ddcar * (AGE - 40) + e_ibw_vc_ddcar * ibwc +
      e_black_vc_ddcar * RACE_BLACK) *
      exp(etalvc_ddcar)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_dcar <- cl_dcar / vc_dcar
    kel_ddcar <- cl_ddcar / vc_ddcar

    # Figure 1 without the boxed compartments that were added only in the
    # updated models (second cariprazine peripheral, DCAR and DDCAR peripherals,
    # DDCAR transit): DCAR elimination feeds DDCAR directly.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_dcar) <- kel * central - kel_dcar * central_dcar
    d/dt(central_ddcar) <- kel_dcar * central_dcar - kel_ddcar * central_ddcar

    # Zero-order input into the depot; dose records must carry rate = -2.
    dur(depot) <- d1

    # mg/L -> ng/mL
    Cc <- 1000 * central / vc
    Cc_dcar <- 1000 * central_dcar / vc_dcar
    Cc_ddcar <- 1000 * central_ddcar / vc_ddcar

    Cc ~ prop(propSd)
    Cc_dcar ~ prop(propSd_dcar)
    Cc_ddcar ~ prop(propSd_ddcar)
  })
}
