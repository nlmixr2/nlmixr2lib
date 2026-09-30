Sun_2021_samidorphan <- function() {
  description <- paste(
    "Two-compartment population PK model for oral samidorphan given as the",
    "olanzapine/samidorphan (OLZ/SAM) combination (or samidorphan alone) in",
    "healthy adults and adults with schizophrenia (Sun 2021; 521 subjects,",
    "11 studies). First-order absorption with an absorption lag time.",
    "Allometric body-weight scaling (fixed exponents 0.75 on CL/F and 1 on",
    "Vc/F, 70 kg reference). Multiplicative categorical effects on CL/F for",
    "rifampin coadministration, moderate hepatic impairment and severe renal",
    "impairment; a fed-state effect on ka; and effects on the lag time for",
    "the samidorphan-alone tablet (vs the OLZ/SAM bilayer tablet) and for",
    "the phase 3 study ALK3831-A305 (imputed dose times)."
  )
  reference <- paste(
    "Sun L, Mills R, Sadler BM, Rege B (2021). Population Pharmacokinetics",
    "of Olanzapine and Samidorphan When Administered in Combination in",
    "Healthy Subjects and Patients With Schizophrenia. J Clin Pharmacol",
    "61(11):1430-1441. doi:10.1002/jcph.1911.",
    sep = " "
  )
  vignette <- "Sun_2021_olanzapine_samidorphan"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "samidorphan", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "samidorphan", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "samidorphan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying ('WT, time-changing body weight' in the Table 4",
        "legend). Allometric power effects on CL/F (exponent fixed at 0.75)",
        "and Vc/F (exponent fixed at 1.0), Table 4 footnote a 'Fixed at",
        "allometric exponent'. Centered at 70 kg: the Figure 3 caption",
        "reference subject weighs 70 kg and the Results quote the body-weight",
        "ratios 'relative to 70 kg'. The Figure 3 reference AUCtau of 284",
        "ng*h/mL for 10 mg once daily agrees with Dose / CL = 10 mg /",
        "35.4 L/h = 282 ng*h/mL (a 76.4 kg median centering would give",
        "301 ng*h/mL). Cohort median 76.4 kg (range 46.5-130.0) per Table 1."
      ),
      source_name = "WT"
    ),
    FED = list(
      description = "Fed-vs-fasted dose-record indicator (1 = fed, 0 = fasted)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Multiplicative power-form effect on ka, 0.107^FED (Table 4 'Food",
        "(fed vs fasted) on Ka'; Table 2 '-90%'). Informed mainly by the",
        "ALK3831-A107 crossover food-effect study (Table S1)."
      ),
      source_name = "FOOD"
    ),
    CONMED_RIFAMPICIN = list(
      description = "Rifampin (rifampicin) coadministration indicator (1 = in the presence of rifampin)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (absence of rifampin)",
      notes = paste(
        "Multiplicative power-form effect on CL/F, 2.70^CONMED_RIFAMPICIN",
        "(Table 4 'Rifampin inducer effect (in the presence vs absence of",
        "rifampin) on CL/F'; Table 2 '+170%'). Estimated from the",
        "ALK-3831-A103 drug-drug interaction study, in which OLZ/SAM was",
        "given with rifampin 600 mg once daily on day 22 after rifampin on",
        "days 15-21 (Table S1), i.e. at established induction. The paper",
        "attributes the effect to CYP3A4 induction."
      ),
      source_name = "RIF"
    ),
    HEPIMP_MOD = list(
      description = "Moderate hepatic impairment indicator (Child-Pugh class B)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal hepatic function)",
      notes = paste(
        "Multiplicative power-form effect on CL/F, 0.810^HEPIMP_MOD (Table 4",
        "'Moderate hepatic impairment ... on CL/F'; Table 2 '-19%').",
        "Moderate impairment is a Child-Pugh score of 7-9 (class B) at",
        "screening (Table 1 footnote c); the 10 subjects came from the",
        "ALK3831-A105 hepatic-impairment study (Table S1)."
      ),
      source_name = "HEPATIC"
    ),
    RENALIMP_SEV = list(
      description = "Severe renal impairment indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal renal function)",
      notes = paste(
        "Multiplicative power-form effect on CL/F, 0.570^RENALIMP_SEV",
        "(Table 4 'Severe renal impairment ... on CL/F'; Table 2 '-43%').",
        "Table 2 describes the contrast as 'severe renal impairment vs",
        "normal renal function in a clinical study', i.e. the severely",
        "impaired subjects of the ALK3831-A106 renal-impairment study",
        "(Table S1). Severe impairment is 15-29 mL/min (CrCl) or",
        "15-29 mL/min/1.73 m^2 (eGFR) per Table 1 footnote d."
      ),
      source_name = "RENAL"
    ),
    FORM_SAM_TAB = list(
      description = "Samidorphan-alone immediate-release tablet (1) vs OLZ/SAM bilayer tablet (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (OLZ/SAM bilayer tablet)",
      notes = paste(
        "Multiplicative power-form effect on ALAG, 1.41^FORM_SAM_TAB (Table",
        "4 'Formulation (samidorphan tablet vs OLZ/SAM bilayer tablet) on",
        "ALAG'; Table 2 'Nonbilayer tablet vs bilayer tablet +41%'). The",
        "samidorphan-alone tablet was given in ALK33-301 and ALK33-B109",
        "(Table S1); 117 of 521 subjects (22%) received it (Table 1)."
      ),
      source_name = "FORM"
    ),
    STUDY_ALK3831A305 = list(
      description = "Phase 3 study ALK3831-A305 record indicator (1 = record from ALK3831-A305)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any of the other 10 pooled studies)",
      notes = paste(
        "Multiplicative power-form effect on ALAG, 10.1^STUDY_ALK3831A305",
        "(Table 4 'Change in ALAG' 10.1; Table 2 '+10-fold'). Table 4",
        "footnote c: dose times were not recorded in ALK3831-A305, so imputed",
        "dose timing was used and the change in ALAG for this study was",
        "estimated. It absorbs dose-time uncertainty rather than describing",
        "a physiological difference; set to 0 for simulation of a patient",
        "with known dose times."
      ),
      source_name = "STUDY"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 521L,
    n_studies = 11L,
    n_observations = 9321L,
    age_range = "18-73 years",
    age_median = "34 years",
    weight_range = "46.5-130.0 kg",
    weight_median = "76.4 kg",
    sex_female_pct = 29,
    race_ethnicity = c(White = 47, Black = 50, `Native American` = 2, Asian = 1, Other = 1),
    disease_state = "Healthy adults (57%) and adults with schizophrenia (43%)",
    dose_range = "Samidorphan 5-30 mg orally (single dose or once daily), as OLZ/SAM bilayer tablet or samidorphan tablet",
    regions = "Not reported (Alkermes-sponsored phase 1 and phase 3 studies)",
    renal_function = "CrCl median 117 mL/min (23-229); 3 subjects with severe impairment by CrCl",
    hepatic_function = "10 subjects with moderate hepatic impairment (Child-Pugh B)",
    notes = paste(
      "Sun 2021 Table 1 and Table S1. The 10 OLZ/SAM studies of the",
      "olanzapine analysis plus the samidorphan-alone study ALK33-B109 in",
      "nondependent recreational opioid users. 11.5% of samidorphan",
      "concentrations were below the 0.250 ng/mL LLOQ and were handled with",
      "the M3 method."
    )
  )

  ini({
    # Structural parameters -- Sun 2021 Table 4 'Estimate' column. Typical values
    # are for the Figure 3 reference subject: 70 kg, normal hepatic and renal
    # function, fasted, no rifampin, OLZ/SAM bilayer tablet, study other than A305.
    lcl <- log(35.4); label("Apparent clearance CL/F (L/h)") # Table 4 'CL/F (L/h)' = 35.4 (RSE 1.65%)
    lvc <- log(297); label("Apparent central volume Vc/F (L)") # Table 4 'Vc/F (L)' = 297 (RSE 1.63%)
    lvp <- log(124); label("Apparent peripheral volume Vp/F (L)") # Table 4 'Vp/F (L)' = 124 (RSE 8.87%)
    lka <- log(6.61); label("First-order absorption rate constant ka (1/h)") # Table 4 'Ka (h)' = 6.61 (RSE 14.2%); unit printed as 'h', a rate constant in 1/h
    ltlag <- log(0.323); label("Absorption lag time ALAG (h)") # Table 4 'ALAG (h)' = 0.323 (RSE 5.57%)
    lq <- log(12.1); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 4 'Q/F (L/h)' = 12.1 (RSE 7.89%)

    # Continuous covariate effects: ln(TVP) = ln(theta_P) + theta_COV * ln(COV / TVCOV).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)") # Table 4 'WT on CL/F' = 0.75 (fixed), footnote a 'Fixed at allometric exponent'
    e_wt_vc <- fixed(1.0); label("Allometric exponent of body weight on Vc/F (unitless)") # Table 4 'WT on Vc/F' = 1.0 (fixed), footnote a

    # Categorical covariate effects: TVP = theta_P * theta_CAT^CAT (Methods).
    e_conmed_rifampicin_cl <- 2.70; label("Rifampin coadministration multiplicative factor on CL/F (power-form base)") # Table 4 'Rifampin inducer effect ... on CL/F' = 2.70 (RSE 3.06%)
    e_hepimp_mod_cl <- 0.810; label("Moderate hepatic impairment multiplicative factor on CL/F (power-form base)") # Table 4 'Moderate hepatic impairment ... on CL/F' = 0.810 (RSE 9.04%)
    e_renalimp_sev_cl <- 0.570; label("Severe renal impairment multiplicative factor on CL/F (power-form base)") # Table 4 'Severe renal impairment ... on CL/F' = 0.570 (RSE 5.96%)
    e_fed_ka <- 0.107; label("Fed-state multiplicative factor on ka (power-form base)") # Table 4 'Food (fed vs fasted) on Ka' = 0.107 (RSE 36.9%)
    # Table 4 marks this row with footnote b ('Fixed at estimate from previous stable
    # model') but also prints an RSE and 95% CI; Table 2 does not flag it as fixed
    # and the Results list only the Q/F IIV as fixed in the samidorphan model, so it
    # is encoded as estimated.
    e_study_alk3831a305_tlag <- 10.1; label("Study ALK3831-A305 multiplicative factor on ALAG (power-form base)") # Table 4 'Change in ALAG' = 10.1 (RSE 11.0%), footnotes b, c
    e_form_sam_tab_tlag <- 1.41; label("Samidorphan-alone tablet multiplicative factor on ALAG (power-form base)") # Table 4 'Formulation (samidorphan tablet vs OLZ/SAM bilayer tablet) on ALAG' = 1.41 (RSE 5.80%)

    # Inter-individual variability: Table 4 reports omega^2 (log-scale variances);
    # CV% = sqrt(exp(omega^2) - 1) for omega^2 > 0.15 (table footnote), e.g.
    # sqrt(exp(1.76) - 1) = 219%. No IIV covariances are reported.
    etalcl ~ 0.087 # Table 4 IIV 'CL/F' = 0.087 (RSE 11.9%), CV 29.4%
    etalvc ~ 0.054 # Table 4 IIV 'Vc/F' = 0.054 (RSE 19.0%), CV 23.3%
    etalka ~ 1.76 # Table 4 IIV 'Ka' = 1.76 (RSE 16.8%), CV 219%
    etaltlag ~ 0.131 # Table 4 IIV 'ALAG' = 0.131 (RSE 24.8%), CV 36.2%
    etalvp ~ 0.681 # Table 4 IIV 'Vp/F' = 0.681 (RSE 24.5%), CV 98.8%
    etalq ~ fixed(0.223) # Table 4 IIV 'Q/F' = 0.223, CV 50.0%; Results: held at 50% to reduce model instability

    # Residual error: Yobs = Ypred * (1 + eps1) (Methods), proportional only.
    propSd <- sqrt(0.061); label("Proportional residual error (fraction)") # Table 4 'Residual variability in sigma^2 prop' = 0.061 (RSE 6.87%) -> SD 0.247 (CV 24.7%)
  })

  model({
    # Individual parameters (Methods 'Model Development' covariate forms; Table 4).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl *
      e_conmed_rifampicin_cl^CONMED_RIFAMPICIN *
      e_hepimp_mod_cl^HEPIMP_MOD * e_renalimp_sev_cl^RENALIMP_SEV
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)
    ka <- exp(lka + etalka) * e_fed_ka^FED
    tlag <- exp(ltlag + etaltlag) * e_study_alk3831a305_tlag^STUDY_ALK3831A305 *
      e_form_sam_tab_tlag^FORM_SAM_TAB

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with lagged first-order absorption (Figure S1).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Dose in mg and volumes in L give mg/L; x 1000 converts to ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
