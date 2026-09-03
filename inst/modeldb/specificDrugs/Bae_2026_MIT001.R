Bae_2026_MIT001 <- function() {
  description <- "Three-compartment population PK model for MIT-001, a mitochondria-targeting ferroptosis inhibitor, with first-order subcutaneous absorption in healthy Korean adults (Bae 2026)"
  reference <- "Bae S, Kim E, Rhee SJ, Kim S, Yu KS, Lee S. Population Pharmacokinetic Analysis of MIT-001, a Novel Ferroptosis Inhibitor, for Dose Optimization. J Clin Pharmacol. 2026;66(4):e70189. doi:10.1002/jcph.70189"
  vignette <- "Bae_2026_MIT001"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Bae 2026 Figure 1 (model schematic:
  # depot -> V1 with Q12 to V2 and Q13 to V3, CL out of V1) and Table 2.
  compartmentData <- list(
    depot       = list(analyte = "MIT-001", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "MIT-001", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "MIT-001", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "MIT-001", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight. Scales both peripheral volumes with reference weight 68.4 kg (the pooled cohort median; Bae 2026 Table 2 footnotes b and c). V2 scales linearly (exponent structurally 1, not estimated); V3 scales with an estimated exponent of 0.742. No weight effect was retained on CL or V1.",
      source_name        = "Weight"
    ),
    ALT = list(
      description        = "Serum alanine aminotransferase activity",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline ALT. Power effect on CL with reference 16 IU/L (the pooled cohort median; Bae 2026 Table 2 footnote a). The exponent is negative, so higher ALT gives lower clearance, which the authors interpret as reduced hepatic metabolic capacity (MIT-001 urinary recovery was only 1%-3%). The cohort ALT range was narrow (6-53 IU/L in healthy participants), so extrapolation to hepatically impaired patients is not supported by these data.",
      source_name        = "ALT"
    )
  )

  # Screened in the covariate analysis (Bae 2026 Methods, "Population
  # Pharmacokinetic Model Development") but NOT retained in the final model.
  # Documented here so the provenance of the paper's covariate screen is
  # preserved without carrying unused-covariate warnings.
  covariatesDataExcluded <- list(
    AST = list(
      description = "Serum aspartate aminotransferase activity",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function candidate on CL; not retained. ALT was the transaminase selected."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened; not retained. Cohort range 4.0-5.1 g/dL (Bae 2026 Table 1). Reported here in the paper's g/dL; the covariate register's canonical unit is g/L (multiply by 10)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened; not retained. Consistent with negligible renal elimination (1%-3% of dose recovered in urine)."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened; not retained. Cohort range 20-44 years."
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened; not retained. Bae 2026 Table 1 reports height in metres (median 1.72-1.75 m, range 1.52-1.88 m); the covariate register's canonical unit is cm."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened; not retained. Body weight was the size descriptor selected, on the peripheral volumes only."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened using the paper's coding SEX = 1 male / 0 female (Bae 2026 Methods); not retained. Only 3 of 119 participants were female, so the cohort could not support a sex effect. The authors explicitly flag this as a limitation."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 119,
    n_studies      = 3,
    n_observations = 2859,
    age_range      = "20-44 years",
    age_median     = "26, 28 and 28.5 years in studies 1, 2 and 3",
    weight_range   = "55-87 kg",
    weight_median  = "67.6, 67.7 and 70.6 kg in studies 1, 2 and 3",
    sex_female_pct = 2.5,
    race_ethnicity = c(Asian = 100),
    disease_state  = "healthy adults (no disease); the model was developed to support dose selection for patients with oral mucositis",
    dose_range     = "IV 0.3-200 mg single dose and 3-30 mg once daily for 7 days; SC 10-40 mg single dose and 20-40 mg once daily for 7 days",
    regions        = "Republic of Korea",
    notes          = "Pooled from three phase 1 randomized, double-blind, placebo-controlled trials in healthy Korean adults: study 1 NCT01737424 (single ascending dose IV, 63 participants, 1185 observations), study 2 NCT03196804 (multiple ascending dose IV, 26 participants, 1038 observations), study 3 NCT05389696 (SC dose escalation with an IV bioavailability arm, 30 participants, 636 observations). Baseline demographics are Bae 2026 Table 1; per-study designs, sampling schedules and sample counts are Supplementary Table 1. The 63 / 26 split for studies 1 and 2 is taken from Supplementary Table 1; the main-text Table 1 'Male, n (%)' row reports these two counts the other way round, and three independent checks favour Supplementary Table 1 (see the vignette's 'Assumptions and deviations' section). The total of 119 participants is unaffected. 95 participants contributed IV data and 30 contributed SC data; 6 study-3 participants received both a 40 mg SC and a 40 mg IV dose after a two-week washout, which is what identifies bioavailability."
  )

  ini({
    # Structural parameters (Bae 2026 Table 2, "Structural model").
    # Reference covariate values: ALT 16 IU/L, body weight 68.4 kg.
    lcl     <- log(2.17);  label("Clearance at ALT 16 IU/L (CL, L/h)")                                  # Table 2, theta1
    lvc     <- log(6.89);  label("Central volume of distribution (V1, L)")                              # Table 2, V1
    lvp     <- log(23.6);  label("First peripheral volume of distribution at 68.4 kg (V2, L)")          # Table 2, theta3
    lvp2    <- log(68.3);  label("Second peripheral volume of distribution at 68.4 kg (V3, L)")         # Table 2, theta4
    lq      <- log(108);   label("Intercompartmental clearance between V1 and V2 (Q12, L/h)")           # Table 2, Q12
    lq2     <- log(32.9);  label("Intercompartmental clearance between V1 and V3 (Q13, L/h)")           # Table 2, Q13
    lka     <- log(0.757); label("First-order subcutaneous absorption rate constant (Ka, 1/h)")         # Table 2, Ka
    lfdepot <- log(0.82);  label("Subcutaneous bioavailability (F, fraction)")                          # Table 2, F

    # Covariate effects (Bae 2026 Table 2 footnotes a, b, c).
    e_alt_cl <- -0.176;      label("Power exponent on (ALT / 16) for CL (unitless)")                    # Table 2, theta2, footnote a
    e_wt_vp  <- fixed(1);    label("Power exponent on (WT / 68.4) for V2 (unitless)")                   # Footnote b prints no exponent: V2 = theta3 * (WT/68.4), i.e. linear in weight, not estimated
    e_wt_vp2 <- 0.742;       label("Power exponent on (WT / 68.4) for V3 (unitless)")                   # Table 2, theta5, footnote c

    # Inter-individual variability. Table 2 reports each IIV as a percent CV,
    # converted here with omega^2 = log(CV^2 + 1). The report scale is pinned
    # by the RSE column: an RSE of 7.8% on the CL IIV is below the theoretical
    # floor sqrt(2/(119-1)) = 13.0% for the RSE of a variance, so the reported
    # estimates and RSEs are both on the SD / CV scale (see the vignette's
    # "Assumptions and deviations" section).
    etalcl     ~ 0.0823622  # 29.3% CV, Table 2, IIV for CL
    etalvc     ~ 0.2183609  # 49.4% CV, Table 2, IIV for V1
    etalvp     ~ 0.0992379  # 32.3% CV, Table 2, IIV for V2
    etalvp2    ~ 0.0259055  # 16.2% CV, Table 2, IIV for V3
    etalq      ~ 0.1021811  # 32.8% CV, Table 2, IIV for Q12
    etalq2     ~ 0.0683649  # 26.6% CV, Table 2, IIV for Q13
    etalfdepot ~ 0.0080674  #  9.0% CV, Table 2, IIV for F

    # Residual error (Bae 2026 Table 2, Residual error row). Reported as
    # Proportional residual error 0.108, read as a proportional SD on the
    # same scale as the IIV block above.
    propSd <- 0.108; label("Proportional residual error (fraction)")                                    # Table 2, proportional residual error
  })
  model({
    # Individual parameters. Covariate forms are exactly Bae 2026 Table 2
    # footnotes a-c; no covariate was retained on V1, Q12, Q13, Ka or F.
    cl  <- exp(lcl  + etalcl)  * (ALT / 16)^e_alt_cl
    vc  <- exp(lvc  + etalvc)
    vp  <- exp(lvp  + etalvp)  * (WT / 68.4)^e_wt_vp
    vp2 <- exp(lvp2 + etalvp2) * (WT / 68.4)^e_wt_vp2
    q   <- exp(lq   + etalq)
    q2  <- exp(lq2  + etalq2)
    ka  <- exp(lka)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Bioavailability applies to the subcutaneous depot only; intravenous
    # doses go directly into central with F = 1.
    f(depot) <- exp(lfdepot + etalfdepot)

    # MIT-001 plasma concentrations are reported in ng/mL. With dose in mg and
    # vc in L, central/vc has units mg/L = ug/mL; multiply by 1000 for ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
