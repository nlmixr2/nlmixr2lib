Sun_2021_olanzapine <- function() {
  description <- paste(
    "Two-compartment population PK model for oral olanzapine given as the",
    "olanzapine/samidorphan (OLZ/SAM) combination (or olanzapine alone) in",
    "healthy adults and adults with schizophrenia (Sun 2021; 601 subjects,",
    "10 studies). First-order absorption with a fixed absorption lag time and",
    "inter-occasion variability on ka. Allometric body-weight scaling (fixed",
    "exponents 0.75 on CL/F and 1 on Vc/F, 70 kg reference) and a power",
    "age effect on Vc/F (36-year reference). Multiplicative categorical",
    "effects on CL/F for rifampin coadministration, smoking, female sex,",
    "Black race, moderate hepatic impairment (fixed) and severe renal",
    "impairment, and a fed-state effect on relative bioavailability."
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
    depot = list(analyte = "olanzapine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "olanzapine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "olanzapine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying ('WT, time-changing body weight' in the Table 3",
        "legend). Allometric power effects on CL/F (exponent fixed at 0.75)",
        "and Vc/F (exponent fixed at 1.0), Table 3 footnote b 'Fixed at",
        "allometric exponent'. Centered at 70 kg: the Figure 3 caption",
        "reference subject weighs 70 kg and the Discussion quotes the",
        "body-weight ratios 'relative to the reference of 70 kg'. The Figure",
        "3 reference AUCtau of 635 ng*h/mL for 10 mg once daily agrees with",
        "Dose / CL = 10 mg / 15.5 L/h = 645 ng*h/mL (a 76.9 kg median",
        "centering would give 694 ng*h/mL). Cohort median 76.9 kg (range",
        "44.0-141.0) per Table 1."
      ),
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on Vc/F, exponent 0.356 (Table 3 'Age on Vc/F'), in",
        "the Methods centered form ln(TVP) = ln(theta_P) + theta_COV *",
        "ln(COV / TVCOV). Centered at 36 years: the Figure 3 caption",
        "reference subject is aged 36 years, equal to the olanzapine",
        "data-set median age (Table 1; Results). Range 18-73 years."
      ),
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Multiplicative power-form effect on CL/F, 0.862^SEXF (Table 3 'Sex",
        "(women vs men) on CL/F'; Table 2 '-14%'). 190 of 601 subjects (32%)",
        "were women (Table 1)."
      ),
      source_name = "SEX"
    ),
    RACE_BLACK = list(
      description = "Black race indicator (1 = Black, 0 = non-Black)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all non-Black races pooled: White, Native American, Asian, Hawaiian, Other)",
      notes = paste(
        "Multiplicative power-form effect on CL/F, 1.10^RACE_BLACK (Table 3",
        "'Race (Black vs non-Black) on CL/F'; Table 2 '+10%'). 255 of 601",
        "subjects (42%) were Black (Table 1)."
      ),
      source_name = "RACE"
    ),
    SMOKE = list(
      description = "Current-smoker indicator (1 = smoker, 0 = nonsmoker or smoking status not recorded)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (nonsmokers pooled with subjects whose smoking status was not recorded)",
      notes = paste(
        "Multiplicative power-form effect on CL/F, 1.30^SMOKE (Table 3",
        "'Smoking (smokers vs nonsmokers) on CL/F'). The reference level",
        "pools nonsmokers with not-recorded status: Table 2 'Smokers vs",
        "nonsmokers/not recorded' and the Discussion '30% higher than",
        "nonsmokers, including 29% of subjects with missing smoking",
        "status'. Table 1: 286 smokers (48%), 142 nonsmokers (24%), 173 not",
        "reported (29%)."
      ),
      source_name = "SMOK"
    ),
    FED = list(
      description = "Fed-vs-fasted dose-record indicator (1 = fed, 0 = fasted)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Multiplicative power-form effect on relative bioavailability,",
        "0.943^FED (Table 3 'Food (fed vs fasted) on F'; Table 2 '-6%').",
        "Informed mainly by the ALK3831-A107 crossover food-effect study",
        "(Table S1)."
      ),
      source_name = "FOOD"
    ),
    CONMED_RIFAMPICIN = list(
      description = "Rifampin (rifampicin) coadministration indicator (1 = in the presence of rifampin)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (absence of rifampin)",
      notes = paste(
        "Multiplicative power-form effect on CL/F, 1.80^CONMED_RIFAMPICIN",
        "(Table 3 'Rifampin inducer effect (in the presence vs absence of",
        "rifampin) on CL/F'; Table 2 '+80%'). Estimated from the",
        "ALK-3831-A103 drug-drug interaction study, in which OLZ/SAM was",
        "given with rifampin 600 mg once daily on day 22 after rifampin on",
        "days 15-21 (Table S1), i.e. at established induction. The paper",
        "attributes the effect to induction of UGT-mediated glucuronidation",
        "and CYP-mediated oxidation of olanzapine."
      ),
      source_name = "RIF"
    ),
    HEPIMP_MOD = list(
      description = "Moderate hepatic impairment indicator (Child-Pugh class B)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal hepatic function)",
      notes = paste(
        "Multiplicative power-form effect on CL/F, 0.875^HEPIMP_MOD, FIXED",
        "(Table 3 'Moderate hepatic impairment ... on CL/F' 0.875 (fixed),",
        "footnote a 'Fixed at estimate from previous stable model where",
        "effect was reliably estimated'). Moderate impairment is a",
        "Child-Pugh score of 7-9 (class B) at screening (Table 1 footnote",
        "c); the 10 subjects came from the ALK3831-A105 hepatic-impairment",
        "study (Table S1)."
      ),
      source_name = "HEPATIC"
    ),
    RENALIMP_SEV = list(
      description = "Severe renal impairment indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal renal function)",
      notes = paste(
        "Multiplicative power-form effect on CL/F, 0.801^RENALIMP_SEV",
        "(Table 3 'Severe renal impairment ... on CL/F'; Table 2 '-20%').",
        "Table 2 describes the contrast as 'severe renal impairment vs",
        "normal renal function in a clinical study', i.e. the severely",
        "impaired subjects of the ALK3831-A106 renal-impairment study",
        "(Table S1). Severe impairment is 15-29 mL/min (CrCl) or",
        "15-29 mL/min/1.73 m^2 (eGFR) per Table 1 footnote d. Mild and",
        "moderate renal impairment were not retained as covariates."
      ),
      source_name = "RENAL"
    ),
    OCC = list(
      description = "Integer-valued occasion index for the inter-occasion variability on ka",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Sun 2021 reports a single inter-occasion variance on ka (Table 3",
        "'Interoccasion variability in Ka' 0.630) but does not define an",
        "occasion or state the number of occasions. Three occasions are",
        "encoded, the largest number of separate PK sampling occasions per",
        "subject in the pooled studies (Table S1: three treatment periods in",
        "ALK3831-A101, sampling days 4, 8 and 13 in ALK3831-A109, three",
        "samples in ALK3831-A305). Each occasion carries its own eta with the",
        "shared variance (NONMEM BLOCK(1) SAME equivalent). Records with OCC",
        "outside 1-3 receive no IOV. For simulation of a single occasion",
        "set OCC = 1 throughout."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    FORM_OLZ_TAB = list(
      description = "Olanzapine-alone immediate-release tablet (1) vs OLZ/SAM bilayer tablet (0)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Formulation was screened (Methods 'Model Development') but not",
        "retained for olanzapine: Table 2 lists no formulation effect for",
        "olanzapine, while the samidorphan model retains one on ALAG."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 601L,
    n_studies = 10L,
    n_observations = 9905L,
    age_range = "18-73 years",
    age_median = "36 years",
    weight_range = "44.0-141.0 kg",
    weight_median = "76.9 kg",
    sex_female_pct = 32,
    race_ethnicity = c(White = 55, Black = 42, `Native American` = 1, Asian = 0.5, Hawaiian = 0.2, Other = 1),
    disease_state = "Healthy adults (41%) and adults with schizophrenia (59%)",
    dose_range = "Olanzapine 5-30 mg orally (single dose or once daily), as OLZ/SAM bilayer tablet or olanzapine tablet",
    regions = "Not reported (Alkermes-sponsored phase 1 and phase 3 studies)",
    renal_function = "CrCl median 117 mL/min (23-260); 3 subjects with severe impairment by CrCl",
    hepatic_function = "10 subjects with moderate hepatic impairment (Child-Pugh B)",
    notes = paste(
      "Sun 2021 Table 1 and Table S1. Nine phase 1 studies with rich",
      "sampling and one phase 3 study (ALK3831-A305) with sparse sampling;",
      "5.6% of olanzapine concentrations were below the 0.250 ng/mL LLOQ",
      "and were handled with the M3 method. Smoking: 48% smokers, 24%",
      "nonsmokers, 29% not reported."
    )
  )

  ini({
    # Structural parameters -- Sun 2021 Table 3 'Estimate' column. Typical values
    # are for the Figure 3 reference subject: 70 kg, 36 years, non-Black,
    # nonsmoking man with normal hepatic and renal function, fasted, no rifampin.
    lcl <- log(15.5); label("Apparent clearance CL/F (L/h)") # Table 3 'CL/F (L/h)' = 15.5 (RSE 2.85%)
    lvc <- log(656); label("Apparent central volume Vc/F (L)") # Table 3 'Vc/F (L)' = 656 (RSE 2.23%)
    lka <- log(0.861); label("First-order absorption rate constant ka (1/h)") # Table 3 'Ka (h)' = 0.861 (RSE 5.70%); unit printed as 'h', a rate constant in 1/h
    ltlag <- fixed(log(0.782)); label("Absorption lag time ALAG (h)") # Table 3 'ALAG (h)' = 0.782 (fixed), footnote a; Results 'lag time for absorption was fixed at 0.782 hour'
    lvp <- log(225); label("Apparent peripheral volume Vp/F (L)") # Table 3 'Vp/F (L)' = 225 (RSE 9.42%)
    lq <- log(6.15); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 'Q/F (L/h)' = 6.15 (RSE 18.4%)

    # Continuous covariate effects: ln(TVP) = ln(theta_P) + theta_COV * ln(COV / TVCOV)
    # (Methods 'Model Development').
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)") # Table 3 'WT on CL/F' = 0.75 (fixed), footnote b 'Fixed at allometric exponent'
    e_wt_vc <- fixed(1.0); label("Allometric exponent of body weight on Vc/F (unitless)") # Table 3 'WT on Vc/F' = 1.0 (fixed), footnote b
    e_age_vc <- 0.356; label("Power exponent of age on Vc/F (unitless)") # Table 3 'Age on Vc/F' = 0.356 (RSE 11.8%)

    # Categorical covariate effects: ln(TVP) = ln(theta_P) + CAT * ln(theta_CAT), i.e.
    # TVP = theta_P * theta_CAT^CAT (Methods 'Model Development').
    e_conmed_rifampicin_cl <- 1.80; label("Rifampin coadministration multiplicative factor on CL/F (power-form base)") # Table 3 'Rifampin inducer effect ... on CL/F' = 1.80 (RSE 4.45%)
    e_smoke_cl <- 1.30; label("Smoker multiplicative factor on CL/F (power-form base)") # Table 3 'Smoking (smokers vs nonsmokers) on CL/F' = 1.30 (RSE 3.64%)
    e_fed_fdepot <- 0.943; label("Fed-state multiplicative factor on relative bioavailability (power-form base)") # Table 3 'Food (fed vs fasted) on F' = 0.943 (RSE 2.26%)
    e_hepimp_mod_cl <- fixed(0.875); label("Moderate hepatic impairment multiplicative factor on CL/F (power-form base)") # Table 3 'Moderate hepatic impairment ... on CL/F' = 0.875 (fixed), footnote a
    e_renalimp_sev_cl <- 0.801; label("Severe renal impairment multiplicative factor on CL/F (power-form base)") # Table 3 'Severe renal impairment ... on CL/F' = 0.801 (RSE 5.67%)
    e_race_black_cl <- 1.10; label("Black race multiplicative factor on CL/F (power-form base)") # Table 3 'Race (Black vs non-Black) on CL/F' = 1.10 (RSE 3.22%)
    e_sexf_cl <- 0.862; label("Female sex multiplicative factor on CL/F (power-form base)") # Table 3 'Sex (women vs men) on CL/F' = 0.862 (RSE 3.56%)

    # Inter-individual variability: Table 3 reports omega^2 (log-scale variances);
    # the CV% column is sqrt(exp(omega^2) - 1) for omega^2 > 0.15 (table footnote),
    # e.g. sqrt(exp(0.171) - 1) = 43.2%. No IIV covariances are reported.
    etalcl ~ 0.171 # Table 3 IIV 'CL/F' = 0.171 (RSE 10.0%), CV 43.2%
    etalvc ~ 0.127 # Table 3 IIV 'Vc/F' = 0.127 (RSE 17.5%), CV 35.6%
    etalka ~ 0.209 # Table 3 IIV 'Ka' = 0.209 (RSE 50.2%), CV 48.2%
    etaltlag ~ 0.319 # Table 3 IIV 'ALAG' = 0.319 (RSE 9.53%), CV 61.3%
    etalvp ~ fixed(0.223) # Table 3 IIV 'VpF' = 0.223, CV 50.0%; Results: held at 50% to reduce model instability
    etalq ~ fixed(0.223) # Table 3 IIV 'Q/F' = 0.223, CV 50.0%; Results: held at 50% to reduce model instability

    # Inter-occasion variability on log-ka, one shared variance across occasions.
    # nlmixr2 has no BLOCK(1) SAME shortcut, so occasions 2 and 3 carry their own
    # etas with the variance fixed equal to the occasion-1 estimate.
    etaiov_ka_1 ~ 0.630 # Table 3 'Interoccasion variability in Ka' = 0.630 (RSE 14.2%), CV 93.7%
    etaiov_ka_2 ~ fixed(0.630) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_ka_3 ~ fixed(0.630) # SAME-equivalent: equal to the occasion-1 IOV variance

    # Residual error: Yobs = Ypred * (1 + eps1) (Methods), proportional only.
    propSd <- sqrt(0.0462); label("Proportional residual error (fraction)") # Table 3 'Residual variability in sigma^2 prop' = 0.0462 (RSE 6.26%) -> SD 0.215 (CV 21.5%)
  })

  model({
    # Occasion indicators multiplexing the IOV etas on log-ka.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3

    # Individual parameters (Methods 'Model Development' covariate forms; Table 3).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl *
      e_conmed_rifampicin_cl^CONMED_RIFAMPICIN * e_smoke_cl^SMOKE *
      e_hepimp_mod_cl^HEPIMP_MOD * e_renalimp_sev_cl^RENALIMP_SEV *
      e_race_black_cl^RACE_BLACK * e_sexf_cl^SEXF
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * (AGE / 36)^e_age_vc
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)
    ka <- exp(lka + etalka) * exp(iov_ka)
    tlag <- exp(ltlag + etaltlag)
    fdepot <- e_fed_fdepot^FED

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with lagged first-order absorption (Figure S1).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # Dose in mg and volumes in L give mg/L; x 1000 converts to ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
