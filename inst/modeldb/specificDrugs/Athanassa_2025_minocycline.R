Athanassa_2025_minocycline <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "linear elimination for orally administered minocycline in critically ill",
    "adults with ventilator-associated pneumonia caused by multidrug-resistant",
    "Acinetobacter baumannii (Athanassa 2025). Minocycline was given by mouth",
    "(dispersed in water, via nasogastric or orogastric tube in mechanically",
    "ventilated patients) as a 200 mg loading dose followed by 100 mg every",
    "12 h. The model is parameterized in apparent oral terms: apparent",
    "clearance CL/F 6.55 L/h, apparent central volume V/F 183.3 L, and",
    "absorption rate constant ka 1.66 1/h, which together imply a terminal",
    "half-life of 19.4 h and an absorption half-life of 0.42 h (25 min), both",
    "as quoted in the paper's Discussion. Between-subject variability is",
    "exponential on CL/F (omega 0.59 on the log scale) and V/F (omega 1.16);",
    "variability on ka was tested but was not estimable from these data and is",
    "absent from the final model. Residual error is proportional (0.616).",
    "Nine covariates (age, sex, body weight, BMI, serum albumin, creatinine,",
    "creatinine clearance, haematocrit and SOFA score) were screened; none was",
    "retained, so the final model is covariate-free (see",
    "covariatesDataExcluded)."
  )
  reference <- paste(
    "Athanassa Z, Papakyriakopoulou P, Marquez Megias S, Saitani EM,",
    "Manioudaki S, Dimoula K, Petsa I, Valsami G, Sakagianni A, Koumaki V,",
    "Dokoumetzidis A, Tsakris A (2025).",
    "Population pharmacokinetic model of oral minocycline in critically ill",
    "adult patients with ventilator-associated pneumonia.",
    "J Antimicrob Chemother 80(6):1420-1426. doi:10.1093/jac/dkaf090",
    sep = " "
  )
  vignette <- "Athanassa_2025_minocycline"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482. Both states hold an amount of minocycline in mg. `depot` is the
  # gastrointestinal absorption site (drug dispersed in water and given orally
  # or down a nasogastric / orogastric tube); `central` is the sampled plasma
  # compartment, scaled by the apparent volume V/F to give the observation Cc.
  compartmentData <- list(
    depot = list(
      analyte = "minocycline", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "minocycline", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  # The final Athanassa 2025 model contains NO covariates. Every covariate the
  # paper recorded was screened on CL/F and V/F by stepwise forward inclusion /
  # backward elimination at P < 0.01 (Methods, "Covariate analysis"), and the
  # Results state that "none of these covariates significantly influenced an
  # apparent clearance or apparent volume of distribution (P > 0.01), apart
  # from serum albumin", which was itself then rejected on plausibility
  # grounds. They are documented here (and not in covariateData) so the
  # provenance of the covariate screen is preserved without declaring
  # covariates that model() never references.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age at study onset.",
      units       = "years",
      type        = "continuous",
      notes       = "Athanassa 2025 Table 1: median 69.5 years (IQR 63-76.5). Screened on CL/F and V/F; not significant at P < 0.01 and not retained (Results, 'Covariate analysis')."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Athanassa 2025 Table 1 reports 17 male of 24 patients, i.e. 7 female (29.2%). The table's printed male percentage of 70.1% appears to be a typographical slip for 70.8% (17/24). Screened on CL/F and V/F; not retained."
    ),
    WT = list(
      description = "Total body weight at study onset.",
      units       = "kg",
      type        = "continuous",
      notes       = "Athanassa 2025 Table 1: median 85.0 kg (IQR 70.0-106.5). Screened on CL/F and V/F; not retained. Note that no allometric scaling is present in the final model either -- the paper reports CL/F and V/F as unscaled population values."
    ),
    BMI = list(
      description = "Body mass index at study onset.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Athanassa 2025 Table 1: median 25 kg/m^2 (IQR 23.95-32.5). Screened on CL/F and V/F; not retained. Obesity was among the most frequent comorbidities in the cohort (Results, 'Patient characteristics')."
    ),
    ALB = list(
      description = "Serum albumin at study onset. Reported by the source in g/dL; the register canonical ALB is g/L (multiply by 10).",
      units       = "g/dL (as reported by the source; canonical ALB is g/L)",
      type        = "continuous",
      notes       = "Athanassa 2025 Table 1: median 2.45 g/dL (IQR 2.20-2.62), i.e. 24.5 g/L. This was the ONLY covariate to reach statistical significance, entering as a power function on the apparent volume of distribution with an estimated exponent of -3.57. The authors rejected it from the final model on plausibility grounds: 'the estimated exponent of -3.57 indicated an excessively strong reciprocal dependence of the volume of distribution on albumin levels. This extreme sensitivity, where a 20% difference in albumin resulted in a 2-fold change in volume, could not be justified' (Results, 'Covariate analysis'). The Discussion attributes the implausible magnitude to the narrow albumin range in this hypoalbuminaemic ICU cohort. The exponent is recorded here for provenance only; it is deliberately NOT implemented in model()."
    ),
    CREAT = list(
      description = "Serum creatinine at study onset.",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Athanassa 2025 Table 1: median 1.00 mg/dL (IQR 0.600-1.98). Screened on CL/F and V/F; not retained."
    ),
    CRCL = list(
      description = "Creatinine clearance, measured directly by 24-hour urine collection (raw mL/min, NOT BSA-normalized).",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Athanassa 2025 Table 1: median 57.7 mL/min (IQR 28.5-104); 11 of 24 patients had CLCR < 50 mL/min without renal replacement therapy. The Methods note the choice of measured rather than estimated CrCl 'for its accuracy in critically ill patients, where serum creatinine-based equations may be less reliable'. Screened on CL/F and V/F; not retained -- consistent with the Discussion's citation of Carney et al. that minocycline is cleared mainly by gastrointestinal excretion rather than renally."
    ),
    HCT = list(
      description = "Haematocrit at study onset.",
      units       = "%",
      type        = "continuous",
      notes       = "Athanassa 2025 Table 1: median 28.2% (IQR 25.6-32.2). Screened on CL/F and V/F; not retained."
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score at ICU admission, a 0-24 composite of six organ-system sub-scores.",
      units       = "(score)",
      type        = "continuous",
      notes       = "Athanassa 2025 Table 1: median 9 (IQR 8-10), indicating severe organ dysfunction. Screened on CL/F and V/F; not retained. SOFA is NOT currently a canonical entry in inst/references/covariate-columns.md (the registered ICU severity scores are SAPS_II and APACHE_II). No new canonical is proposed here because the covariate was rejected by the source authors and is never referenced in model(); covariatesDataExcluded is documentation only and is not validated against the register. A future extraction that RETAINS a SOFA effect should propose SOFA as a canonical at that time."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24,
    n_studies      = 1,
    n_observations = 172,
    age_median     = "69.5 years (IQR 63-76.5)",
    weight_median  = "85.0 kg (IQR 70.0-106.5)",
    bmi_median     = "25 kg/m^2 (IQR 23.95-32.5)",
    sex_female_pct = 29.2,
    disease_state  = paste(
      "Critically ill mechanically-ventilated ICU adults with",
      "ventilator-associated pneumonia caused by multidrug-resistant",
      "Acinetobacter baumannii; median SOFA score 9 (IQR 8-10) at admission;",
      "84% admitted for COVID-19 infection; VAP diagnosed after a median ICU",
      "stay of 16 days"
    ),
    renal_function = paste(
      "Measured creatinine clearance (24-h urine collection) median",
      "57.7 mL/min (IQR 28.5-104); 11 of 24 patients had CLCR < 50 mL/min,",
      "none on renal replacement therapy"
    ),
    dose_range     = "200 mg oral loading dose, then 100 mg every 12 h for at least 5 days",
    co_medication  = "Colistin 67%, meropenem 29%, ampicillin-sulbactam 29%",
    regions        = "Greece (single centre: Sismanogleio General Hospital, Athens; August 2021 to April 2023)",
    notes          = paste(
      "Baseline demographics from Athanassa 2025 Table 1. 182 plasma",
      "minocycline concentrations were collected; 10 (5.49%) below the limit",
      "of quantification were discarded, leaving 172 used for model",
      "development. Sampling was sparse: 20 patients were sampled at 12, 13,",
      "14, 18, 21.5, 48, 59.5, 61, 64, 68 and 71.5 h and 4 patients at 0.5, 1,",
      "2, 4, 6 and 12 h. Minocycline MICs were 4 mg/L or lower in 13 strains",
      "(susceptible) and 8 mg/L in 11 strains (intermediate). Serum albumin was",
      "low throughout (median 2.45 g/dL), which the authors cite as the reason",
      "the albumin-on-volume effect was poorly identified."
    )
  )

  ini({
    # Structural parameters -- Athanassa 2025 Table 2 ("PK parameter estimates
    # of the final PK model"), "Final model / Estimate (%RSE)" column. The
    # abstract restates all three: "183.3 L, 6.55 L/h and 1.66 h-1, for the
    # apparent volume of distribution (V/F), the apparent clearance (CL/F) and
    # the absorption rate constant (ka)".
    #
    # These are APPARENT oral parameters (relative to the unknown oral
    # bioavailability F). No separate bioavailability term is estimated, so no
    # lfdepot appears below and F is implicitly 1; the Discussion cites a
    # literature oral bioavailability of about 90% but does not fit it.
    lka <- log(1.66);  label("First-order absorption rate constant ka (1/h)")               # Athanassa 2025 Table 2: Ka = 1.66 (41.9% RSE); bootstrap median 1.61, 95% CI 0.47-3.55. Discussion cross-check: absorption half-life log(2)/1.66 = 0.42 h (25 min).
    lcl <- log(6.55);  label("Apparent oral clearance CL/F (L/h)")                          # Athanassa 2025 Table 2: CL/F = 6.55 (17.2% RSE); bootstrap median 6.48, 95% CI 2.90-9.22.
    lvc <- log(183.3); label("Apparent central volume of distribution V/F (L)")             # Athanassa 2025 Table 2: V/F = 183.3 (25.9% RSE); bootstrap median 183.4, 95% CI 102.0-374.7. Discussion cross-check: terminal half-life log(2)*183.3/6.55 = 19.4 h, matching the quoted "estimated average t1/2 of 19.4 h".
    #
    # No covariate-effect parameters: the final model is covariate-free. See
    # covariatesDataExcluded above for the nine screened covariates, including
    # the rejected albumin-on-volume power effect (exponent -3.57).

    # Between-subject variability -- Athanassa 2025 Table 2 rows "IIV-CL/F" and
    # "IIV-V/F", "Final model / Estimate (%RSE)" column. Exponential IIV on
    # CL/F and V/F only; the Results state that IIV "could not be estimated on
    # ka", so no etalka term exists in the final model and none is added here.
    # No off-diagonal covariance is reported, so the two etas are independent.
    #
    # SCALE. The Estimate column holds Monolix's omega, i.e. the STANDARD
    # DEVIATION of the random effect on the log scale, not a variance. The CL/F
    # row proves this: sqrt(exp(0.59^2) - 1) * 100 = 64.5%, which is the 64%
    # printed in the adjacent "%CV" column (and back-solving 64% returns
    # omega = 0.5859 ~ 0.59). Reading the column as a variance is falsified --
    # it would give a CL/F CV of 89.7%, not 64%.
    #
    # The V/F row's printed "%CV" of 116% is internally inconsistent with that
    # formula: omega = 1.16 implies a CV of 168.5%, whereas 116% is exactly
    # omega * 100, i.e. the naive approximation applied to the V/F row only.
    # The Estimate of 1.16 is corroborated twice independently and is used
    # here: its own RSE implies a 95% CI of 0.787-1.533, and the bootstrap
    # column gives a median of 1.13 with a 95% CI of 0.55-1.42 -- both centred
    # on 1.16, not on the 0.9233 that a literal 116% CV would require. See the
    # vignette "Assumptions and deviations" section.
    #
    # nlmixr2 ini() takes VARIANCES, so each omega is squared inline to keep
    # the published standard deviation visible in the source trace.
    etalcl ~ 0.59^2  # Athanassa 2025 Table 2: IIV-CL/F omega = 0.59 (23.8% RSE), %CV 64%, eta-shrinkage 16.9%; bootstrap median 0.57, 95% CI 0.15-1.22. Variance = 0.3481.
    etalvc ~ 1.16^2  # Athanassa 2025 Table 2: IIV-V/F omega = 1.16 (16.4% RSE), eta-shrinkage -2.7%; bootstrap median 1.13, 95% CI 0.55-1.42. Variance = 1.3456. Printed %CV of 116% is the naive omega*100; the lognormal CV is 168.5%.

    # Residual error -- Athanassa 2025 Table 2 row "Proportional error". The
    # Results state "A proportional error model was found to be most
    # appropriate for explaining the residual error", i.e. Monolix's
    # y = f * (1 + b * eps), which is nlmixr2's prop() error model with
    # propSd = b. No additive component was retained.
    propSd <- 0.616; label("Proportional residual error (fraction)")  # Athanassa 2025 Table 2: 0.616 (5.9% RSE); bootstrap median 0.612, 95% CI 0.523-0.692.
  })

  model({
    # Individual parameters. Exponential (log-normal) IIV on CL/F and V/F; ka
    # has no random effect because it was not estimable (Results, "popPK
    # model").
    ka  <- exp(lka)
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc + etalvc)

    # Micro-constant
    kel <- cl / vc

    # One-compartment model with first-order absorption and linear elimination
    # (Results, "popPK model"; abstract "A one-compartment model with
    # first-order absorption and linear elimination best described the data").
    # No lag time: the Methods tested models "with and without lag-time
    # (Tlag)" and the final model has none.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation. Apparent volume, so Cc is the plasma concentration
    # predicted for the administered oral dose without an explicit F term.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
