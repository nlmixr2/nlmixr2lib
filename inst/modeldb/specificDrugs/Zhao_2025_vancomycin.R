Zhao_2025_vancomycin <- function() {
  description <- "One-compartment IV-infusion population PK model for vancomycin in 126 non-extremely preterm neonates treated in a Chinese neonatal intensive care unit (Zhao 2025). Clearance is 0.14 L/h at the cohort mean covariate values and scales as a power function of body weight (reference 2.12 kg, exponent 1.13), serum creatinine (reference 30.52 umol/L, exponent -0.15) and daily fluid input (reference 367.18 mL/24h, exponent 0.14), and is multiplied by exp(-0.20) when a diuretic is coadministered. Central volume is 1.04 L scaling with body weight (reference 2.12 kg, exponent 1.07). Daily fluid input and diuretic use are the novel covariates this paper contributes; postmenstrual age, albumin, blood urea nitrogen, urine volume and respiratory support were screened but not retained. NOTE: the published equations 6 and 7 print the covariate ratios WITHOUT their superscript exponents, which were lost in typesetting; the exponents are taken from Table 3 and are confirmed by back-calculation from the Table 4 dosing grid (see the vignette)."
  reference <- "Zhao K, Zhao F, Ju K, Chen H, Zhai X, Chang Y, Liu Z. Population pharmacokinetics of vancomycin in non-extremely preterm neonates based on real-world studies: influence of daily fluid input and diuretics. Microbiol Spectr. 2025;13(6):e02274-24. doi:10.1128/spectrum.02274-24"
  vignette <- "Zhao_2025_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at the start of vancomycin treatment",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Zhao 2025 Table 1: median 1.98 kg (IQR 1.35-2.98) for all 126 patients, 1.96 kg (IQR 1.35-2.94) for the 112-patient modeling group. Birth weight (median 1.45 kg) is a separate column and is NOT the covariate used here; equations 6 and 7 use WT, the weight at the start of vancomycin medication. The normalising constant 2.12 kg printed in equations 6 and 7 is the cohort MEAN weight, not the median (the same mean-not-median convention applies to the Scr and DFI reference values). The strongest single covariate in the univariate screen for both CL (dOFV -150.148) and V (dOFV -73.992), Table 2. Enters CL as (WT/2.12)^1.13 and V as (WT/2.12)^1.07; the exponents come from Table 3 because the printed equations lost their superscripts.",
      source_name        = "WT"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Zhao 2025 Table 1: median 32.48 umol/L (IQR 24.54-42.31) for all patients, 32.31 (IQR 24.42-41.76) for the modeling group. Reported in umol/L throughout, so no mg/dL conversion applies. The normalising constant in equation 6 is 30.52 umol/L (the cohort mean). Exclusion criterion (ii) removed serum creatinine values measured within 7 days of birth, because neonatal creatinine over that window still reflects maternal creatinine rather than the neonate's own renal function; the retained values are therefore post-equilibration. Enters CL as (CREAT/30.52)^-0.15. The negative exponent is the expected direction for a renally eliminated drug: higher creatinine means lower clearance. Univariate dOFV on CL -7.215 (Table 2); entered last on forward inclusion (step 8, dOFV -7.869) and survived backward elimination.",
      source_name        = "Scr"
    ),
    FLUID_IN_24H = list(
      description        = "Total fluid administered to the neonate over 24 hours (enteral plus parenteral), recorded at the start of vancomycin treatment",
      units              = "mL/24h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Zhao 2025 Table 1 'Daily fluid input (mL)': median 364 mL (IQR 252.91-473.10) for all patients, 365.15 (IQR 255.64-478.10) for the modeling group. The normalising constant in equation 6 is 367.18 mL/24h (the cohort mean). This is a whole-body daily INTAKE volume, not weight-normalised and not a net balance; it is one of the two covariates the paper contributes as novel (the other is diuretic use). Zhao 2025 Discussion: 'DFI was used instead of infusion volume as a more accurate measure of fluid balance in NICU neonates', and 'An increase in fluid volume results in elevated urine output and enhanced glomerular filtration, subsequently increasing the clearance of the hydrophilic drug vancomycin', which is the mechanism behind the positive exponent. The authors address the obvious collinearity with body weight explicitly: the weight-DFI correlation R = 0.7274 is weaker than the weight-PMA correlation R = 0.8696, so DFI was judged to carry information independent of weight. Enters CL as (FLUID_IN_24H/367.18)^0.14. Univariate dOFV on CL -20.894 (Table 2). Note that daily URINE volume (median 200.5 mL) was screened separately and NOT retained; it is documented in covariatesDataExcluded as URINE_VOL_24H.",
      source_name        = "DFI"
    ),
    CONMED_DIURETIC = list(
      description        = "Concomitant diuretic indicator: 1 = a diuretic was coadministered at the start of vancomycin treatment, 0 = no concomitant diuretic",
      units              = "(binary)",
      type               = "categorical",
      reference_category = "0 (no concomitant diuretic)",
      notes              = "Class composition in Zhao 2025 Materials and Methods: furosemide, spironolactone and hydrochlorothiazide, i.e. loop plus thiazide plus potassium-sparing, pooled into a single yes/no indicator. Table 1: 42/126 (33.3%) of all patients, 40/112 (35.7%) of the modeling group. Enters CL multiplicatively as exp(A) with A = -0.20 when a diuretic is used and A = 0 otherwise (Zhao 2025 text immediately following equation 6, and Table 3 row 'DA on CL' = -0.20, RSE 20.03%), i.e. an 18.1% reduction in clearance. The paper abbreviates the covariate 'DA' throughout Table 2, Table 3, Table 4 and Fig. 1. Zhao 2025 Discussion attributes the reduction to diuretic-associated renal injury (acute interstitial nephritis) rather than to a direct interaction, and concludes that the vancomycin dose should be reduced when a diuretic is coadministered; the Table 4 dosing grid does exactly that. Univariate dOFV on CL -4.609 (Table 2), entered second on forward inclusion (step 2, dOFV -15.724) and survived backward elimination.",
      source_name        = "DA"
    )
  )

  # Screened during covariate model building but NOT retained in the final
  # model, so they are documentation only and are not referenced in model().
  # Sources: Zhao 2025 Materials and Methods (covariate screening list),
  # Table 1 (distributions) and Table 2 (univariate dOFV and the stepwise
  # forward-inclusion / backward-elimination sequence).
  covariatesDataExcluded <- list(
    PAGE = list(
      description = "Postmenstrual age at the start of vancomycin treatment",
      units       = "weeks",
      type        = "continuous",
      notes       = "Zhao 2025 Table 1: median 35.7 weeks (IQR 32.9-39.98). The second strongest univariate covariate on CL (dOFV -113.109) and on V (dOFV -39.657), and it entered the multivariate model on forward inclusion (step 6, dOFV -4.543), but it was the FIRST covariate removed on backward elimination (step 1, dOFV +2.804 < 7.88). This is the notable negative result of the paper: most neonatal vancomycin models retain a PMA maturation term, and Zhao 2025 does not, having conditioned on weight, creatinine, fluid input and diuretic use first. Registered canonical PAGE is documented in weeks here rather than the register's default months, matching the source.",
      source_name = "PMA"
    ),
    GA = list(
      description = "Gestational age at delivery",
      units       = "weeks",
      type        = "continuous",
      notes       = "Zhao 2025 Table 1: median 32 weeks (IQR 29.38-37.85). Collected as a demographic characteristic; not carried into the covariate screen as a separate term (prematurity was screened as the binary PM indicator instead)."
    ),
    PRETERM = list(
      description = "Preterm birth indicator: 1 = born before 37 weeks gestation",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Zhao 2025 Table 1 'Prematurity': 93/126 (73.8%), defined in footnote b as birth before 37 weeks gestation. Abbreviated PM in Table 2. Univariate dOFV -36.347 on CL and -20.705 on V; not carried into the multivariate model. Documentation only, so this name is deliberately NOT registered in inst/references/covariate-columns.md."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Zhao 2025 Table 1: 86/126 male (68.3%), so 40/126 female (31.7%). Screened as SEX in Table 2 with essentially no signal (dOFV -0.169 on CL, -1.965 on V); not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Zhao 2025 Table 1: median 28.5 g/L (IQR 25.07-31.55). Third strongest continuous univariate signal on CL (dOFV -18.17) and a signal on V (-8.144), but not carried into the multivariate model. Zhao 2025 Discussion contrasts this with Smits et al., who identified albumin as the most critical covariate for FREE vancomycin concentration; this cohort measured total concentration, and the Discussion records that clearance rose with albumin in the univariate screen but that albumin was not significant thereafter."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Zhao 2025 Table 1: median 3.57 mmol/L (IQR 2.31-5.86). No univariate signal on CL (dOFV -1.022) or V (-0.023), but it nonetheless entered the multivariate model on V at forward-inclusion step 4 (dOFV -6.347) before being removed at backward-elimination step 3 (dOFV +6.047 < 7.88)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Zhao 2025 Table 1: median 10.69 U/L (IQR 6.46-15.53). Univariate dOFV -5.607 on CL, -1.482 on V; not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Zhao 2025 Table 1: median 28.95 U/L (IQR 21.38-41.9). Essentially no univariate signal (dOFV -0.003 on CL, -0.943 on V); not retained."
    ),
    URINE_VOL_24H = list(
      description = "Daily urine volume",
      units       = "mL/24h",
      type        = "continuous",
      notes       = "Zhao 2025 Table 1 'Daily urine volume (mL)': median 200.5 mL (IQR 156.75-283.25). Abbreviated DUV. The second strongest univariate signal on CL (dOFV -48.043) and a strong one on V (-23.771), yet it did not enter the multivariate model at all once weight was in. Zhao 2025 Discussion records this as a deliberate replication of a prior negative result: 'our findings indicated that urine volume was not a significant covariate affecting clearance, in agreement with the results of the aforementioned reports'. This is the OUTPUT counterpart of the retained input covariate FLUID_IN_24H; the paper retains the intake and rejects the output."
    ),
    MECH_VENT = list(
      description = "Respiratory support indicator (oxygen therapy or mechanical ventilation)",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Zhao 2025 Table 1: 38/126 (30.2%). Abbreviated RS in Table 2. Univariate dOFV -5.495 on CL and -5.828 on V; not carried into the multivariate model. Named MECH_VENT here for documentation only; note that the source pools oxygen therapy with mechanical ventilation, which is broader than the registered MECH_VENT canonical (invasive mechanical ventilation only), which is a further reason it is not written into model()."
    ),
    NCIS = list(
      description = "Neonatal critical illness score at admission, three-level (>90 non-critical, 70-90 critical, <70 extremely critical)",
      units       = "(score)",
      type        = "categorical",
      notes       = "Zhao 2025 Table 1: >90 in 10/126 (7.9%), 70-90 in 88/126 (69.9%), <70 in 28/126 (22.2%); the only baseline characteristic differing significantly between the modeling and validation groups (P = 0.038). Univariate dOFV -2.626 on CL and -1.831 on V; not retained. Zhao 2025 lists as a limitation that NCIS was assessed at admission rather than at the time of study inclusion, 'which may have led to the omission of certain significant covariates'. Documentation only, so this name is deliberately NOT registered in inst/references/covariate-columns.md."
    ),
    CONMED_VASOACTIVE = list(
      description = "Concomitant vasoactive drug indicator (dopamine, dobutamine, epinephrine, norepinephrine)",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Zhao 2025 Table 1: 40/126 (31.7%). Abbreviated VAA in Table 2. Univariate dOFV -1.992 on CL and -0.348 on V; not retained. Zhao 2025 Discussion explicitly contrasts this null result with Tang et al., who reported that concomitant vasoactive drugs decrease vancomycin clearance, and notes that this cohort's 31.7% exposure rate was HIGHER than the 21.4% in that study, so the null is not an exposure-prevalence artifact. Documentation only, so this name is deliberately NOT registered in inst/references/covariate-columns.md."
    ),
    CONMED_ALBUMIN = list(
      description = "Concomitant human serum albumin administration indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Zhao 2025 Table 1: 25/126 (19.8%). Abbreviated HA in Table 2. Exactly zero univariate signal on CL (dOFV 0) and a modest one on V (-5.145); it nonetheless entered the multivariate model on V at forward-inclusion step 3 (dOFV -6.709) before being removed at backward-elimination step 4 (dOFV +5.651 < 7.88). Documentation only, so this name is deliberately NOT registered in inst/references/covariate-columns.md."
    ),
    CONMED_PIPTAZ = list(
      description = "Concomitant piperacillin-tazobactam indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Zhao 2025 Table 1: 12/126 (9.5%). Abbreviated PTZ in Table 2. Essentially no univariate signal (dOFV -0.001 on both CL and V); it nonetheless entered the multivariate model on CL at forward-inclusion step 5 (dOFV -5.639) before being removed at backward-elimination step 2 (dOFV +5.458 < 7.88). Documentation only, so this name is deliberately NOT registered in inst/references/covariate-columns.md."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 126L,
    n_studies        = 1L,
    n_centers        = 1L,
    n_concentrations = 276L,
    age_range        = "Postnatal age at admission median 0 days (IQR 0-2.8); inclusion required postnatal age at or below 28 days. Gestational age at delivery median 32 weeks (IQR 29.38-37.85). Postmenstrual age at the start of vancomycin median 35.7 weeks (IQR 32.9-39.98).",
    weight_range     = "Weight at the start of vancomycin median 1.98 kg (IQR 1.35-2.98); cohort mean 2.12 kg (the normalising constant in equations 6 and 7). Birth weight median 1.45 kg (IQR 1.14-2.88).",
    sex_female_pct   = 31.7,
    race_ethnicity   = "Not reported; single-center Chinese cohort (Xi'an, Shaanxi), so presumed predominantly Han Chinese.",
    disease_state    = "Neonates admitted to a neonatal intensive care unit and treated with intravenous vancomycin for suspected or confirmed gram-positive infection. Preterm (born before 37 weeks) 93/126 (73.8%); extremely low birth weight 23/126 (18.3%), very low birth weight 44/126 (34.9%), low birth weight 20/126 (15.9%). Infection site: bloodstream 95/126 (75.4%), neural 32/126 (25.4%), other 21/126 (16.7%). Neonatal critical illness score at admission: >90 (non-critical) 10/126, 70-90 (critical) 88/126, <70 (extremely critical) 28/126. Extremely preterm infants were NOT included, which is what the title's 'non-extremely preterm' refers to.",
    renal_function   = "Serum creatinine median 32.48 umol/L (IQR 24.54-42.31); blood urea nitrogen median 3.57 mmol/L (IQR 2.31-5.86); daily urine volume median 200.5 mL (IQR 156.75-283.25). Patients with congenital renal dysplasia or with chronic or acute renal insufficiency (with or without renal replacement therapy) were excluded, and creatinine values measured within 7 days of birth were discarded because they still reflect maternal creatinine.",
    co_medication    = "Diuretics (furosemide, spironolactone, hydrochlorothiazide) 42/126 (33.3%); vasoactive drugs (dopamine, dobutamine, epinephrine, norepinephrine) 40/126 (31.7%); human serum albumin 25/126 (19.8%); piperacillin-tazobactam 12/126 (9.5%); respiratory support 38/126 (30.2%). Mild hypothermia therapy (5/126), non-steroidal anti-inflammatory drugs (3/126) and cimetidine (3/126) were too infrequent to be screened as covariates.",
    dose_range       = "Intravenous vancomycin 10-15 mg/kg per dose every 8 to 12 hours, infused over 1 hour. Initial dose median 25 mg (IQR 18-36.75). Daily fluid input median 364 mL (IQR 252.91-473.10); cohort mean 367.18 mL/24h.",
    regions          = "China (Northwest Women's and Children's Hospital, Xi'an, Shaanxi). Retrospective real-world cohort, January 2019 to December 2023.",
    notes            = "126 patients contributing 276 vancomycin concentrations. Split by calendar time rather than at random: 112 patients (January 2019 to June 2023) built the model and 14 patients (July 2023 to December 2023, 24 concentrations) were held out for external validation. 143/276 (51.8%) of the concentrations were troughs; peaks and troughs were drawn 0.5 h after and 0.5 h before an infusion respectively, all at steady state after at least four doses. Assay: chemiluminescence immunoassay (VIVA, Siemens), calibration range 2-50 mg/L; concentrations outside that range were excluded. Estimation in Phoenix NLME 8.3.5.340 by first-order conditional estimation with extended least squares. One- and two-compartment structures were both tried and one compartment was retained. Model evaluation: 1000-iteration bootstrap (Table 3), visual predictive check with 1000 simulations (Fig. 3), and normalized prediction distribution errors (Fig. 4; t-test P = 0.0723, Fisher variance test P = 0.0711, Shapiro-Wilk P = 1). External validation on the 24 held-out concentrations: mean prediction error 2.74%, mean absolute prediction error 17.48%, F20% 75.00%, F30% 83.33%."
  )

  ini({
    # Structural parameters. Zhao 2025 Table 3, "Final model (RSE%)" column,
    # with the bootstrap median and 95% CI from the "Bootstrap (95% CI)"
    # column quoted alongside each value.
    #
    # IMPORTANT -- the printed equations lost their exponents. Zhao 2025
    # equations 6 and 7 are typeset as
    #   CL(L/h) = 0.14 x (WT/2.12) x (Scr/30.52) x (DFI/367.18) x e^A x exp(etaCL)
    #   V       = 1.04 x (WT/2.12)
    # with NO superscript on any of the parenthesised ratios. The omission is
    # in the publisher's own equation artwork (the EuropePMC image bundle for
    # PMC12131725 ships equations 1-7 as spectrum.02274-24.m001-m007.jpg, and
    # m006/m007 are missing the superscripts too), so it is a production error
    # rather than a text-extraction artifact. The exponents are the Table 3
    # covariate rows, and they are confirmed numerically by back-calculating
    # them out of the Table 4 dosing grid -- see the vignette section
    # "The missing exponents". Reading the equations literally (all exponents
    # equal to 1) is falsified twice over by Table 4: it inverts the sign of
    # the creatinine effect and inflates the fluid-input effect roughly
    # sevenfold.
    lcl <- log(0.14); label("Clearance at the reference covariate values (L/h)")       # Zhao 2025 Table 3: CL = 0.14 (RSE 3.15%); bootstrap 0.14 (95% CI 0.13-0.15)
    lvc <- log(1.04); label("Central volume at the reference body weight (L)")         # Zhao 2025 Table 3: V  = 1.04 (RSE 4.28%); bootstrap 1.03 (95% CI 0.97-1.10)

    # Covariate effects on clearance. All three continuous covariates enter as
    # power functions of the covariate divided by its cohort MEAN (2.12 kg,
    # 30.52 umol/L, 367.18 mL/24h), so cl reduces to exp(lcl) = 0.14 L/h for a
    # non-diuretic subject sitting at all three reference values -- which is
    # exactly the "typical CL value of 0.14 L/hour" the Results paragraph
    # after equation 7 reports.
    e_wt_cl <- 1.13; label("Body-weight power exponent on clearance (unitless)")                     # Zhao 2025 Table 3 "WTonCL": 1.13 (RSE 5.64%); bootstrap 1.15 (95% CI 1.01-1.36)
    e_creat_cl <- -0.15; label("Serum-creatinine power exponent on clearance (unitless)")            # Zhao 2025 Table 3 "Scr on CL": -0.15 (RSE printed as -31.34%); bootstrap -0.15 (95% CI -0.25 to -0.017)
    e_fluid_in_24h_cl <- 0.14; label("Daily-fluid-input power exponent on clearance (unitless)")     # Zhao 2025 Table 3 "DFI on CL": 0.14 (RSE 27.03%); bootstrap 0.13 (95% CI printed as "-0.087 to -0.21", a sign typo for 0.087 to 0.21 given the positive point estimate)
    e_conmed_diuretic_cl <- -0.20; label("Log-scale shift in clearance with a concomitant diuretic (unitless)")  # Zhao 2025 Table 3 "DA on CL": -0.20 (RSE printed as -20.03%); bootstrap -0.20 (95% CI -0.39 to -0.051); text after equation 6: "When diuretics were used concomitantly, A = -0.20, and when diuretics were not used, A = 0"

    # Covariate effect on central volume. Body weight is the only covariate
    # retained on V (backward-elimination step 5, "No effect chosen").
    e_wt_vc <- 1.07; label("Body-weight power exponent on central volume (unitless)")  # Zhao 2025 Table 3 "WTonV": 1.07 (RSE 7.61%); bootstrap 1.07 (95% CI 0.92-1.22)

    # Interindividual variability. Zhao 2025 Table 3 reports IIV on the "%CV"
    # scale (the Table 3 abbreviation footnote defines "%CV, coefficient of
    # variation"), with the exponential eta model of equation 1
    # (P_i = P_TV * exp(eta_i)). The final-model CL value is 4.97% CV, so
    #   omega^2 = log(1 + 0.0497^2) = 0.00246704  (omega = 0.04967).
    # For a CV this small the lognormal back-transform and the naive
    # omega = CV reading differ by 0.07%, so the choice is immaterial here;
    # the back-transform is used because it is the exact inverse of the
    # reported quantity.
    #
    # A 4.97% CV is strikingly small for a neonatal popPK model, and it is not
    # a transcription slip: the Results paragraph after equation 7 states it
    # in prose -- "the random inter-individual variability in CL in the final
    # model was significantly lower (24.84% versus 4.97%)" -- and the Table 3
    # bootstrap column corroborates it at 4.80% (95% CI 2.83-6.77). The base
    # model's 24.84% CV is the conventional-looking value; almost all of it is
    # absorbed by the four covariates. IIV on V was 0.90% CV with an RSE of
    # 215.55% in the base model and was dropped from the final model
    # altogether, so V carries no eta.
    etalcl ~ 0.00246704  # Zhao 2025 Table 3, final model "CL (%CV)" = 4.97 (RSE 15.70%); bootstrap 4.80 (95% CI 2.83-6.77); omega^2 = log(1 + 0.0497^2)

    # Proportional residual error. Zhao 2025 Results: "The examination of
    # residual variability supported the adoption of a proportional error
    # model", i.e. equation 3 Y = F * (1 + eps). Table 3 reports it on the
    # same %CV scale as the IIV rows, so the residual SD is 0.180.
    propSd <- 0.180; label("Proportional residual error (fraction)")  # Zhao 2025 Table 3, final model "Proportional (%CV)" = 18.0 (RSE 5.76%); bootstrap 17.6 (95% CI 14.81-20.69)
  })

  model({
    # Vancomycin clearance, Zhao 2025 equation 6 with the Table 3 exponents
    # restored (see the ini() comment and the vignette for why the printed
    # equation has none). WT in kg, CREAT in umol/L, FLUID_IN_24H in mL/24h,
    # CONMED_DIURETIC a 0/1 indicator. The diuretic term is the paper's e^A
    # with A = -0.20 * CONMED_DIURETIC, an 18.1% reduction in clearance.
    cl <- exp(lcl + etalcl) *
      (WT / 2.12)^e_wt_cl *
      (CREAT / 30.52)^e_creat_cl *
      (FLUID_IN_24H / 367.18)^e_fluid_in_24h_cl *
      exp(e_conmed_diuretic_cl * CONMED_DIURETIC)

    # Central volume, Zhao 2025 equation 7 with the Table 3 exponent restored.
    # No interindividual variability was retained on V.
    vc <- exp(lvc) * (WT / 2.12)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
