AbouAuda_2024_gentamicin <- function() {
  description <- "One-compartment intravenous population PK model for gentamicin in pediatric patients aged 1-14 years with and without acute lymphoblastic leukemia (ALL), fitted to routine therapeutic-drug-monitoring peak and trough concentrations (Abou-Auda 2024). Clearance scales allometrically with body weight (exponent 0.75, reference 20 kg) and as a power function of Schwartz-estimated creatinine clearance (exponent 0.437, reference 120 mL/min); volume of distribution scales linearly with body weight (exponent 1, reference 20 kg). ALL status was screened but NOT retained: the paper's central finding is that gentamicin volume of distribution and clearance do not differ between ALL and non-ALL pediatric patients, so dosing requirements are the same in both groups. Residual variability is proportional."
  reference <- "Abou-Auda HS, Alotaibi F, Alsanea S, Alwhaibi A, Almutairi MM, Alrabiah Z, Alsultan A, Al Jeraisy M. Population pharmacokinetics of gentamicin in acute lymphoblastic leukemia pediatric patients compared to non-oncology patients. Saudi Pharm J. 2024;32:102060. doi:10.1016/j.jsps.2024.102060"
  vignette <- "AbouAuda_2024_gentamicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Gentamicin was given as a 30-minute IV infusion, so
  # the dose enters `central` directly and there is no depot state. The
  # specimen is verified: Abou-Auda 2024 Methods (Subjects) states that
  # SERUM samples were assayed by latex inhibition immunoassay (TDX, Abbott
  # Laboratories).
  compartmentData <- list(
    central = list(analyte = "gentamicin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Abou-Auda 2024 Table 1: mean (SD) 22.1 (17.1) kg over all 115 patients; 22.23 (20.3) kg in the ALL group and 22.02 (12.2) kg in the non-ALL group (p = 0.949). Reference weight 20 kg is the normalizing constant printed in both final-model equations in the Results (`(weight/20)`). Note that the cohort mean of 22.1 kg is NOT the reference: the mean is strongly right-skewed (SD 17.1 kg over ages 1-14 years) and the paper chose a round 20 kg. Weight was significant on both CL and Vd in the univariate screen (Table 2, CL R^2 = 0.38) and was retained on both parameters in the final model.",
      source_name        = "weight"
    ),
    CRCL = list(
      description        = "Creatinine clearance estimated with the modified Schwartz equation (Schwartz 2009)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Abou-Auda 2024 Methods (Data collection) and Table 1 footnote: 'Creatinine clearance (Clcr) was calculated using the modified Schwartz formula (Schwartz et al., 2009).' Reference value 120 mL/min is the normalizing constant printed in the final-model clearance equation in the Results (`(CLcr/120)`); it is not tabulated anywhere and is not equal to any reported cohort summary. UNITS CAVEAT: Table 1 heads the column 'Cl_cr (ml/min)', but the modified/bedside Schwartz 2009 equation returns a BSA-normalized value in mL/min/1.73 m^2. This extraction records the paper's own printed unit (mL/min) because that is what the source states, but a user supplying real data should be aware that the values the model was fitted to are almost certainly BSA-normalized Schwartz estimates; supplying a raw (un-normalized) creatinine clearance would silently rescale the renal term. See the vignette Assumptions and deviations section. Cohort values (Table 1): 149 (56.5) over all patients, 153.73 (62.8) in ALL, 164.88 (43.5) in non-ALL (p = 0.006) -- note that the printed all-patient mean of 149 lies BELOW both subgroup means, which is arithmetically impossible; the n-weighted mean of the two subgroups is 158.8 mL/min. This affects only the descriptive population metadata, not the model, whose reference value is the independently printed 120 mL/min. CLcr was significant on CL in the univariate screen (Table 2, R^2 = 0.35) and was one of the two covariates retained in the forward multiple regression (Table 2, combined R^2 = 0.69). It was NOT significant on Vd.",
      source_name        = "CLcr"
    )
  )

  # Screened in the Abou-Auda 2024 covariate analysis (Methods, PK analysis:
  # "The tested covariates included age, sex, height, serum albumin,
  # diagnosis, Scr concentration, ALL status, and Clcr") but NOT retained in
  # the final model, so these are documentation only and are not referenced
  # in model(). Only weight (on CL and Vd) and CLcr (on CL) survived.
  covariatesDataExcluded <- list(
    DIS_ALL = list(
      description        = "Acute lymphoblastic leukemia disease-state indicator (1 = ALL, 0 = non-ALL)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-ALL; the complement group here is the paper's non-oncology pediatric comparison cohort)",
      notes              = "THE PAPER'S PRIMARY RESEARCH QUESTION, and a deliberate negative result. ALL status ('oncology status') was significant on CL in the univariate screen (Table 2, p < 0.05) but had by far the weakest explanatory power of the four univariate covariates (R^2 = 0.07 versus 0.46 for age, 0.38 for weight, 0.35 for CLcr) and did NOT survive the forward multiple linear regression, in which only weight and CLcr remained (Table 2). Results: 'The average Cl in ALL vs. non-ALL patients was 1.9 vs. 2.2 L/hour (p value > 0.05)' and 'The average Vd in ALL vs. non-ALL patients was 6.3 vs. 6.4 L (p value > 0.05)'. The authors conclude that dosing requirements are the same in both groups, which CONTRADICTS the earlier Llanos-Paez 2020 pediatric finding of a 15% lower central Vd and 32% lower CL in oncology patients (see the shipped LlanosPaez_2020_gentamicin.R for that model). n = 63 ALL (54.8%) versus 52 non-ALL (45.2%).",
      source_name        = "diagnosis / ALL status / oncology status"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL and Vd (Methods, PK analysis). Age had the HIGHEST univariate R^2 of any covariate on CL (Table 2, R^2 = 0.46, p < 0.05) but was not retained in the forward multiple regression, which kept weight and CLcr; age is strongly collinear with body weight across a 1-14 year cohort. Table 1: mean (SD) 6.1 (4.2) years overall; 5.36 (4.48) in ALL and 6.88 (3.8) in non-ALL (p = 0.055).",
      source_name        = "age"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male",
      notes              = "Screened (Methods, PK analysis, listed as 'sex'); not retained. Table 1 reports males: 63 (54.8%) overall, 34 (53.9%) in ALL, 29 (55.7%) in non-ALL, so the cohort is 45.2% female. The source reports male counts; canonical SEXF is the complement (SEXF = 1 - SEXM).",
      source_name        = "sex"
    ),
    HT = list(
      description        = "Body height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL and Vd (Methods, PK analysis and Data collection); not retained. No height summary statistic is reported in Table 1, although height is an input to both the Schwartz CLcr estimate and the BSA calculation that ARE tabulated.",
      source_name        = "height"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL and Vd (Methods, PK analysis and Data collection); not retained. No albumin summary statistic is reported in Table 1 and no unit is stated, so the canonical SI unit is recorded here without a source-confirmed value.",
      source_name        = "serum albumin"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL and Vd (Methods, PK analysis); not retained as a standalone covariate, being an input to the Schwartz CLcr estimate that WAS retained. Table 1: mean (SD) 38.8 (12.9) umol/L overall; 41.85 (15.6) in ALL and 35.13 (7.58) in non-ALL (p = 0.005), i.e. significantly higher serum creatinine in the ALL group, which the authors attribute to chemotherapy-associated renal compromise.",
      source_name        = "Scr"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 115L,
    n_studies       = 1L,
    n_sites         = 1L,
    age_range       = "1 to 14 years (inclusion criterion); Table 1 mean (SD) 6.1 (4.2) years. Age strata: toddlers 18 (15.7%), children 88 (76.5%), adolescents 9 (7.8%).",
    weight_range    = "Full range not reported; Table 1 mean (SD) 22.1 (17.1) kg. BMI 16.9 (5.9); body surface area 0.8 (0.4) m^2; ideal body weight 25.5 (13.8) kg.",
    sex_female_pct  = 45.2,
    race_ethnicity  = "Not reported. Single-center Saudi Arabian cohort.",
    disease_state   = "Pediatric inpatients receiving intravenous gentamicin for more than 72 h, split into 63 (54.8%) with acute lymphoblastic leukemia and 52 (45.2%) without. Gentamicin was given as empirical therapy under a febrile-neutropenia protocol or for any suspected gram-negative bacterial infection. Excluded: critically ill children with renal failure, liver dysfunction, or burns, and samples drawn at inappropriate (early or late) peak/trough times.",
    dose_range      = "2.5 mg/kg per dose as a 30-minute intravenous infusion every 8 h for more than 72 h. Table 1 actual dose: 54.3 (44.2) mg/dose, i.e. 2.5 (0.64) mg/kg. All patients were dosed every 8 h.",
    regions         = "Saudi Arabia (King Abdullah Specialist Children's Hospital, Ministry of National Guard, Riyadh)",
    renal_function  = "Schwartz-estimated creatinine clearance 149 (56.5) mL/min overall, skewing supranormal; 153.73 (62.8) in ALL versus 164.88 (43.5) in non-ALL (p = 0.006). Serum creatinine 38.8 (12.9) umol/L overall; blood urea nitrogen 2.8 (1.9) mmol/L overall, 3.18 (1.9) in ALL versus 2.38 (1.7) in non-ALL (p = 0.026). Patients with renal failure were excluded by design, so the model carries no information about renal impairment. NOTE: the printed all-patient CLcr mean of 149 lies below BOTH subgroup means (153.73 and 164.88), which is arithmetically impossible; the n-weighted mean of the subgroups is 158.8 mL/min. Recorded as printed; see the vignette Assumptions and deviations section.",
    notes           = "Retrospective single-center cross-sectional study (IRB approval NRC21R/127/04, informed consent waived) using routine therapeutic-drug-monitoring data. Serum samples were drawn 30 minutes before the next dose (trough) and 1 h after the end of an infusion (peak), around the third dose to ensure steady state, and assayed by latex inhibition immunoassay (TDX, Abbott Laboratories). Sampling is therefore SPARSE (nominally two concentrations per subject) and the total number of concentration records is not reported. The paper states that the sparse data made models with more than one compartment infeasible, and that the one-compartment structure follows the authors' prior publication (Alsultan 2019). Fitted in Monolix 2023R1 with log-normally distributed PK parameters; constant, proportional, and combined residual-error models were tested. Covariate effects were assessed by linear regression on the empirical Bayes estimates with forward selection at alpha = 0.05, rather than within the population model itself. Model evaluation was by goodness-of-fit plots (Figs. 1 and 2) and a visual predictive check (Fig. 3); no numeric NCA or predictive-check summary is tabulated."
  )

  ini({
    # Structural parameters. The reference subject weighs 20 kg and has a
    # Schwartz-estimated creatinine clearance of 120 mL/min. Both reference
    # constants and both intercepts come from the two final-model equations
    # printed in the Results section immediately before the sentence "The PK
    # parameter estimates are shown in Table 3":
    #
    #   Cl = 2.13 * (weight/20)^0.75 * (CLcr/120)^0.437
    #   Vd = 7.3  * (weight/20)
    #
    # NOTE ON THE CLEARANCE INTERCEPT (2.13 versus 2.3). The paper is
    # internally inconsistent: the printed equation uses 2.13 L/h, while
    # Table 3 lists "Cl (L/hr) 2.3" (RSE 2.7%). The volume intercept 7.3 L
    # is identical in both places, which establishes that Table 3 is
    # reporting the equations' reference-subject intercepts rather than some
    # other quantity. This extraction uses the EQUATION value 2.13, per the
    # standing convention that a printed equation outranks a conflicting
    # table/text value, and because three independent checks favour it:
    #
    #   (1) Internal consistency. The Vd equation reproduces the reported
    #       EBE mean Vd (6.3 L in ALL, 6.4 L in non-ALL; n-weighted 6.345 L)
    #       at a typical weight of 17.4 kg, a plausible median for a cohort
    #       whose mean weight of 22.1 kg carries an SD of 17.1 kg. At that
    #       same 17.4 kg, reproducing the reported EBE mean CL (1.9 and
    #       2.2 L/h; n-weighted 2.036 L/h) requires CLcr = 137.6 mL/min with
    #       an intercept of 2.13, versus 115.4 mL/min with 2.3. Only the
    #       former is consistent with the cohort's reported CLcr (mean 149,
    #       right-skewed); 115.4 would sit below the model's own 120 mL/min
    #       reference.
    #   (2) Allometric back-scaling. 2.13 * (70/20)^0.75 = 5.45 L/h =
    #       90.8 mL/min at an adult 70 kg, matching the paper's own
    #       Introduction ("In healthy individuals, the clearance (Cl) of
    #       gentamicin is approximately 90 ml/min"). An intercept of 2.3
    #       gives 98.1 mL/min.
    #   (3) 2.13 does not round to 2.3, so the discrepancy is a
    #       transcription slip in one of the two places rather than a
    #       rounding difference.
    #
    # See the vignette Assumptions and deviations section.
    lcl <- log(2.13); label("Clearance at WT = 20 kg and CLcr = 120 mL/min (L/h)")  # Abou-Auda 2024 Results, final-model Cl equation (Table 3 lists 2.3; see note above)
    lvc <- log(7.3);  label("Volume of distribution at WT = 20 kg (L)")             # Abou-Auda 2024 Results, final-model Vd equation; Table 3: 7.3 L (RSE 3.6%)

    # Body-size exponents. Both are printed inside the final-model equations
    # with no standard error and both take their canonical theoretical
    # allometric values (0.75 on clearance, 1 on volume), so they are
    # encoded as fixed. The Vd equation is printed as `7.3 * (weight/20)`
    # with NO superscript on the closing parenthesis, i.e. an exponent of
    # exactly 1. Both equations are display mathematics that the markdown
    # preprocessor drops entirely (`<!-- formula-not-decoded -->`); they are
    # recovered with `pdftotext -layout`, which typesets the superscripts
    # onto the line above:
    #
    #   (        )0.75 (      )0.437
    #     weight         CLcr
    #   Cl = 2.13 x               x
    #       20           120
    #   (            )
    #      weight
    #   Vd = 7.3 x
    #        20
    #
    # Note that the exponent line renders for the Cl equation but there is
    # NO corresponding superscript on the Vd term in the same extraction,
    # which is what establishes the Vd exponent as exactly 1 rather than as
    # a dropped superscript. Same treatment as the structurally identical
    # Jung_2024_vancomycin.R.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/20) for CL (unitless)")  # Abou-Auda 2024 Results, final-model Cl equation: (weight/20)^0.75
    e_wt_vc <- fixed(1);    label("Allometric exponent on (WT/20) for Vd (unitless)")  # Abou-Auda 2024 Results, final-model Vd equation: (weight/20), no superscript

    # Renal-function effect on clearance. Estimated (a three-decimal,
    # non-theoretical value), so it is NOT wrapped in fixed(), even though
    # Table 3 omits it and therefore reports no RSE for it. Table 3 lists
    # only five rows (V, IIV V, Cl, IIV Cl, residual B) and omits every
    # covariate coefficient and reference constant, so the absence of this
    # exponent from the table is a gap in the table rather than evidence
    # that the exponent was fixed.
    e_crcl_cl <- 0.437; label("Power exponent on (CRCL/120) for CL (unitless)")  # Abou-Auda 2024 Results, final-model Cl equation: (CLcr/120)^0.437

    # Between-subject variability. Table 3 footnote: "IIV is expressed as
    # the coefficient of variation". Monolix's log-normal random-effect
    # parameterization X_i = X_pop * exp(eta_i) reports %CV as
    # sqrt(exp(omega^2) - 1) * 100, so omega^2 = log(CV^2 + 1). nlmixr2
    # takes the eta line as a VARIANCE on the log scale.
    etalcl ~ 0.0246575  # log(0.158^2 + 1); Abou-Auda 2024 Table 3: IIV Cl 15.8% CV (RSE 13.7%)
    etalvc ~ 0.0194104  # log(0.140^2 + 1); Abou-Auda 2024 Table 3: IIV V 14% CV (RSE 28%)

    # Residual variability. Table 3 reports a single parameter "Residual
    # variability B" = 0.3 with the footnote "B is residual unexplained
    # variability expressed as a proportional error". In Monolix's error
    # model y = f + b*f*e with e ~ N(0,1), b IS the proportional standard
    # deviation on the fraction scale, which is nlmixr2's propSd. The
    # Methods state that constant, proportional, and combined error models
    # were tested; the reported final model carries the proportional term
    # only, so no additive term is included.
    propSd <- 0.3; label("Proportional residual error (fraction)")  # Abou-Auda 2024 Table 3: B = 0.3 (RSE 8.2%)
  })

  model({
    # Individual PK parameters, exactly as printed in Abou-Auda 2024 Results:
    #   Cl = 2.13 * (weight/20)^0.75 * (CLcr/120)^0.437
    #   Vd = 7.3  * (weight/20)
    cl <- exp(lcl + etalcl) * (WT / 20)^e_wt_cl * (CRCL / 120)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 20)^e_wt_vc

    # Derived first-order elimination rate constant.
    kel <- cl / vc

    # One-compartment IV disposition. Gentamicin was administered as a
    # 30-minute intravenous infusion, so the dose enters `central` via the
    # event record (rate or duration set on the dose row); there is no
    # absorption compartment.
    d/dt(central) <- -kel * central

    # Dose in mg, vc in L, so central/vc is mg/L (equivalently ug/mL, the
    # conventional unit for gentamicin therapeutic drug monitoring).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
