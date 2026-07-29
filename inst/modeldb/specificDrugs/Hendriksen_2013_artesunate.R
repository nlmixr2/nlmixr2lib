Hendriksen_2013_artesunate <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model of intramuscular artesunate (ARS)",
    "and its active metabolite dihydroartemisinin (DHA) in 70 African children aged",
    "7 months to 11 years admitted with severe Plasmodium falciparum malaria",
    "(Hendriksen 2013). Each species has a one-compartment apparent-volume disposition;",
    "ARS is delivered by a zero-order input over a 1-min fixed duration (the IM",
    "absorption from the injection site, fixed because too few samples were collected",
    "during the absorption phase to identify the rate) and is converted mole-for-mole",
    "to DHA with no separate parent elimination. Body weight is the dominant covariate",
    "(allometric scaling with fixed exponents 0.75 on apparent clearance and 1.0 on",
    "apparent volume for both species; reference 10.9 kg), with hemoglobin additionally",
    "lowering DHA clearance by 10.2% per g/dL above the reference 7.1 g/dL."
  )
  reference <- paste(
    "Hendriksen IC, Mtove G, Kent A, Gesase S, Reyburn H, Lemnge MM, Lindegardh N,",
    "Day NP, von Seidlein L, White NJ, Dondorp AM, Tarning J. Population",
    "pharmacokinetics of intramuscular artesunate in African children with severe",
    "malaria: implications for a practical dosing regimen.",
    "Clin Pharmacol Ther. 2013;93(5):443-450. doi:10.1038/clpt.2013.26"
  )
  vignette  <- "Hendriksen_2013_artesunate"
  units     <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight at admission",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with fixed exponents 0.75 on apparent CL and 1.0 on apparent Vc applied to both artesunate (parent) and dihydroartemisinin (DHA); reference weight 10.9 kg is the typical-value patient body weight (Hendriksen 2013 Table 2 footnote a). Time-fixed at admission in the source data set.",
      source_name        = "WT"
    ),
    HGB = list(
      description        = "Blood hemoglobin concentration at admission",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on DHA apparent clearance, centred on the reference hemoglobin 7.1 g/dL (Hendriksen 2013 Table 2 footnote a): cl_dihydroart = cl_dihydroart_typ * (1 + e_hgb_cl_dihydroart * (HGB - 7.1)) with e_hgb_cl_dihydroart = -0.102, i.e. DHA clearance increases 10.2% per g/dL decrease in hemoglobin (Results, p.445). Source paper reports hemoglobin in g/dL; multiply g/L values by 0.1 if simulating from SI inputs.",
      source_name        = "HGB"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 70L,
    n_studies      = 1L,
    n_observations = 274L,
    age_range      = "7 months to 11 years (median 2.5 years; 84% under 5 years)",
    weight_range   = "median 10.8 kg (IQR 9 to 13.5)",
    weight_typical = "10.9 kg (population typical value for allometric reference)",
    sex_female_pct = NA_real_,
    race_ethnicity = "African (Tanzanian)",
    disease_state  = "Severe Plasmodium falciparum malaria meeting one or more WHO severity criteria (coma, prostration, convulsions, respiratory distress / acidotic breathing, shock, severe anaemia, hypoglycaemia, hemoglobinuria, severe jaundice, or anuria/oliguria); case fatality 13% during this substudy.",
    hgb_range      = "median 7.1 g/dL (IQR 5.1 to 9.2; modeling range 2.72 to 13.6 g/dL)",
    dose_range     = "2.4 mg/kg intramuscular artesunate at 0, 12, and 24 h (then daily) per AQUAMAT protocol; this PK analysis used post-first-dose samples only.",
    sampling       = "Sparse: pre-dose plus four post-first-dose samples per patient in random windows 0-1, 1-4, 4-12, and 12-24 h; 274 quantifiable ARS and DHA concentrations across 70 patients.",
    regions        = "Tanzania (Teule Hospital, Muheza); part of the multicentre AQUAMAT trial (ISRCTN50258054).",
    notes          = "Demographics from Hendriksen 2013 Table 1; AQUAMAT PK substudy May 2009 to July 2010. Three patients (4.3%) were HIV-positive (none on antiretroviral therapy). 31 patients had received any prior oral antimalarial and 12 patients had received i.m. quinine within 24 h before admission; the model does not distinguish pretreatment subgroups."
  )

  ini({
    # Structural parameters (Hendriksen 2013 Table 2, "Population estimate" column).
    # Typical values are reported for a patient with WT = 10.9 kg and HGB =
    # 7.1 g/dL (Table 2 footnote a). The allometric / hemoglobin scaling in
    # model() centres these covariates so that exp(lcl) / exp(lvc) / exp(lcl_dihydroart)
    # / exp(lvc_dihydroart) reproduce the published typical values exactly when WT =
    # 10.9 and HGB = 7.1.
    lcl     <- log(45.8)
    label("Apparent artesunate elimination clearance at WT = 10.9 kg, CL/F (L/h)")  # Hendriksen 2013 Table 2: CL/F ARS = 45.8 L/h (%RSE 8.10)
    lvc     <- log(28.2)
    label("Apparent artesunate central volume at WT = 10.9 kg, V/F (L)")  # Hendriksen 2013 Table 2: V/F ARS = 28.2 L (%RSE 11.4)
    lcl_dihydroart <- log(22.4)
    label("Apparent DHA elimination clearance at WT = 10.9 kg and HGB = 7.1 g/dL, CL_DHA/F (L/h)")  # Hendriksen 2013 Table 2: CL/F DHA = 22.4 L/h (%RSE 8.40)
    lvc_dihydroart <- log(13.5)
    label("Apparent DHA central volume at WT = 10.9 kg, V_DHA/F (L)")  # Hendriksen 2013 Table 2: V/F DHA = 13.5 L (%RSE 9.69)

    # Zero-order IM input duration, fixed at 1 min = 1/60 h in the source
    # paper because too few samples were collected during the absorption
    # phase to estimate it (Hendriksen 2013 Results p.444 / Table 2: DUR
    # 1.00 min FIXED).
    ldur    <- fixed(log(1 / 60))
    label("IM zero-order input duration (h); fixed at 1 min")

    # Allometric exponents, fixed in the source paper ("body weight was used
    # as a fixed allometric function on all elimination clearance (power of
    # 0.75) and apparent volume of distribution (power of 1) parameters",
    # Hendriksen 2013 Results p.444).
    e_wt_cl     <- fixed(0.75)
    label("Allometric WT exponent on artesunate CL/F (unitless)")  # Hendriksen 2013 Results p.444
    e_wt_vc     <- fixed(1.00)
    label("Allometric WT exponent on artesunate V/F (unitless)")  # Hendriksen 2013 Results p.444
    e_wt_cl_dihydroart <- fixed(0.75)
    label("Allometric WT exponent on DHA CL/F (unitless)")  # Hendriksen 2013 Results p.444
    e_wt_vc_dihydroart <- fixed(1.00)
    label("Allometric WT exponent on DHA V/F (unitless)")  # Hendriksen 2013 Results p.444

    # Hemoglobin effect on DHA clearance. Paper text (Results p.445):
    # "DHA clearance increased 10.2% per unit (g/dl) of decrease of
    # hemoglobin." Encoded as a linear-deviation effect centred on
    # HGB_ref = 7.1 g/dL (Table 2 footnote a) so a negative coefficient
    # corresponds to higher CL at lower HGB:
    #   cl_dihydroart = cl_dihydroart_typ * (1 + e_hgb_cl_dihydroart * (HGB - 7.1))
    # At HGB = 6.1, multiplier = 1 + (-0.102) * (-1) = 1.102 (10.2% increase).
    e_hgb_cl_dihydroart <- -0.102
    label("Effect of HGB on DHA CL/F (per g/dL above 7.1)")  # Hendriksen 2013 Table 2: covariate effect "Negative effect of hemoglobin on CL/F DHA (%)" = 10.2 (%RSE 14.9); paper sign convention: increase per g/dL decrease in HGB

    # Inter-individual variability. Hendriksen 2013 Table 2 reports
    # variances for the three estimated etas and a single off-diagonal
    # entry labelled "correlation of random effects on CL/F and V/F" in
    # the footnote. The off-diagonal value 0.497 is therefore the
    # correlation coefficient; the covariance for the nlmixr2 lower-
    # triangle is corr * sqrt(var_CL * var_V) = 0.497 * sqrt(0.415 * 0.680)
    # = 0.2640 (see vignette Errata for the alternative covariance-direct
    # reading; the published validation results in the paper assume the
    # correlation reading).
    etalcl + etalvc ~ c(0.415,
                        0.2640, 0.680)  # Hendriksen 2013 Table 2: var(eta_CL_ARS) = 0.415, var(eta_V_ARS) = 0.680, corr = 0.497
    etalcl_dihydroart     ~ 0.306  # Hendriksen 2013 Table 2: var(eta_CL_DHA) = 0.306 (%RSE 37.9)

    # Residual error. The source paper modelled log-natural plasma
    # concentrations with additive residual on the log scale (NONMEM
    # "additive-on-log" = nlmixr2 proportional on linear scale per
    # references/naming-conventions.md). Table 2 reports the variance on
    # the log scale (sigma_ARS = 0.0942, sigma_DHA = 0.211); the propSd
    # parameters here are the corresponding SDs (sqrt of variance) on the
    # log scale, which equal the proportional CV in linear space to first
    # order.
    propSd     <- sqrt(0.0942)
    label("Proportional residual SD for artesunate plasma concentration (SD on log scale)")  # Hendriksen 2013 Table 2: sigma ARS = 0.0942 (variance, %RSE 29.2); SD = sqrt(0.0942) = 0.307
    propSd_dihydroart <- sqrt(0.211)
    label("Proportional residual SD for DHA plasma concentration (SD on log scale)")  # Hendriksen 2013 Table 2: sigma DHA = 0.211 (variance, %RSE 12.5); SD = sqrt(0.211) = 0.459
  })

  model({
    # Molecular weights (g/mol). Source paper modelled both species on a
    # molar basis with complete 1:1 in-vivo conversion of ARS to DHA. This
    # implementation tracks amounts in mass units (mg) for ease of dosing
    # in clinical mg/kg units, applying the mole-for-mole conversion via
    # the MW ratio at the metabolite-formation step.
    mw_ars <- 384.42  # artesunate, C19H28O8
    mw_dihydroart <- 284.35  # dihydroartemisinin, C15H24O5

    # Individual PK parameters (Hendriksen 2013 Results p.444). Allometric
    # scaling on CL with exponent 0.75 and on V with exponent 1, centred
    # on WT = 10.9 kg. Hemoglobin enters as a linear-deviation effect on
    # DHA CL/F centred on HGB = 7.1 g/dL.
    cl     <- exp(lcl     + etalcl)     * (WT / 10.9)^e_wt_cl
    vc     <- exp(lvc     + etalvc)     * (WT / 10.9)^e_wt_vc
    cl_dihydroart <- exp(lcl_dihydroart + etalcl_dihydroart) * (WT / 10.9)^e_wt_cl_dihydroart *
              (1 + e_hgb_cl_dihydroart * (HGB - 7.1))
    vc_dihydroart <- exp(lvc_dihydroart)              * (WT / 10.9)^e_wt_vc_dihydroart

    # Micro-constants. Under the paper's assumption of complete in-vivo
    # conversion of ARS to DHA, all artesunate clearance is metabolic
    # conversion, so the rate from central -> central_dihydroart is kel_ars =
    # cl/vc. DHA is eliminated linearly at kel_dihydroart = cl_dihydroart/vc_dihydroart.
    kel     <- cl     / vc
    kel_dihydroart <- cl_dihydroart / vc_dihydroart

    # Zero-order IM input duration (fixed at 1 min).
    dur_im <- exp(ldur)

    # ODE system. Dose records target central with dur(central) <- dur_im
    # for zero-order delivery; ARS is then converted mole-for-mole to DHA
    # (mass-rate factor mw_dihydroart / mw_ars).
    d/dt(central)     <- -kel * central
    d/dt(central_dihydroart) <-  kel * central * (mw_dihydroart / mw_ars) -
                          kel_dihydroart * central_dihydroart

    dur(central) <- dur_im

    # Plasma concentrations. With dose in mg and volume in L, central / vc
    # is mg/L (= ug/mL); multiply by 1000 to convert to ng/mL when
    # comparing against Hendriksen 2013 Table 2 post-hoc estimates
    # (Cmax_ARS = 943 ng/mL, AUC_ARS = 570 ng*h/mL, Cmax_DHA = 547 ng/mL,
    # AUC_DHA = 890 ng*h/mL).
    Cc     <- central     / vc
    Cc_dihydroart <- central_dihydroart / vc_dihydroart

    Cc     ~ prop(propSd)
    Cc_dihydroart ~ prop(propSd_dihydroart)
  })
}
