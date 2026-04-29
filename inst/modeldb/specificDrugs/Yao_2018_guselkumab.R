Yao_2018_guselkumab <- function() {
  description <- "One-compartment population PK model with first-order SC absorption and first-order elimination for guselkumab (anti-IL-23 p19 human IgG1-lambda mAb) in adults with moderate-to-severe plaque psoriasis (pooled phase 2 X-PLORE and phase 3 VOYAGE 1 / VOYAGE 2 trials; Yao 2018)"
  reference <- "Yao Z, Hu C, Zhu Y, Xu Z, Randazzo B, Wasfi Y, Chen Y, Sharma A, Zhou H. Population Pharmacokinetic Modeling of Guselkumab, a Human IgG1-lambda Monoclonal Antibody Targeting IL-23, in Patients with Moderate to Severe Plaque Psoriasis. J Clin Pharmacol. 2018;58(5):613-627. doi:10.1002/jcph.1063"
  vignette <- "Yao_2018_guselkumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F and V/F per Yao 2018 Table 4 footnotes b and c; reference value 87.1 kg (median of pooled X-PLORE + VOYAGE 1 + VOYAGE 2 PopPK analysis population, Table 2). Time-fixed at baseline.",
      source_name        = "BWT"
    ),
    DIAB = list(
      description        = "Past or current diabetes mellitus comorbidity",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no diabetes mellitus)",
      notes              = "Multiplicative effect on CL/F per Yao 2018 Table 4 footnote b: CL/F = 0.516 * (BWT/87.1)^0.998 * 1.12^DIAB * 1.11^RACE. With DIAB=1, CL/F is 12% higher than the non-diabetic reference (Results section, page 622, and Table 4). Diabetes prevalence in the analysis population was 8.9% (Table 3).",
      source_name        = "DIAB"
    ),
    RACE_WHITE = list(
      description        = "White race indicator (1 = White, 0 = non-White)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-White, pooling Asian, Black, and Other per Table 3 of Yao 2018)",
      notes              = "Yao 2018 Table 4 footnote b reports the race effect as 1.11^RACE where the source RACE indicator equals 1 for non-White and 0 for White (the 'Race (nonwhite)' row in Table 4). The canonical column in inst/references/covariate-columns.md is RACE_WHITE (1 = White), so values are inverted: RACE_WHITE = 1 - source_RACE. The model applies the 1.11 multiplier when RACE_WHITE = 0 via 1.11^(1 - RACE_WHITE). Non-White subjects (16.2% of the analysis population per Table 3: 12.2% Asian, 1.8% Black, 2.3% Other) had 11% higher CL/F than White subjects.",
      source_name        = "RACE"
    )
  )

  population <- list(
    n_subjects     = 1454L,
    n_studies      = 3L,
    age_range      = "18-82 years (median 44, IQR 34-53)",
    age_median     = 44,
    weight_range   = "45-198 kg (median 87.1, IQR 74.8-100)",
    weight_median  = 87.1,
    sex_female_pct = 29.2,
    race_ethnicity = "White 83.8%, Asian 12.2%, Black 1.8%, Other 2.3% (Yao 2018 Table 3)",
    disease_state  = "Adults with moderate-to-severe plaque psoriasis. Pooled data from X-PLORE (NCT01483599; phase 2 dose-ranging), VOYAGE 1 (NCT02207231; phase 3), and VOYAGE 2 (NCT02207244; phase 3).",
    dose_range     = "Phase 2 X-PLORE: 5, 15, 50, 100, or 200 mg SC at varying schedules through week 40. Phase 3 VOYAGE 1/2: 100 mg SC at weeks 0 and 4 then every 8 weeks. Sparse PK sampling.",
    regions        = "Multinational phase 2/3 trials (X-PLORE, VOYAGE 1, VOYAGE 2).",
    diabetes_pct   = 8.9,
    ada_positive_pct = 5.4,
    ethnicity_hispanic_pct = 6.3,
    notes          = "Demographics from Yao 2018 Tables 2 and 3. 13,014 PK records from 1454 patients in the final analysis dataset (Table 1). Median age 44 years, median weight 87.1 kg, 70.8% male. Diabetes prevalence 8.9%. ADA-positive prevalence 5.4%. Model-derived elimination half-life is approximately 18.1 days (Results, page 619). Body weight was the primary covariate, accounting for 28% (CL/F) and 32% (V/F) of the proportion of variance for IIV; diabetes and race had marginal but retained effects on CL/F (10% threshold for inclusion in the final reduced model)."
  )

  ini({
    # Structural parameters from Yao 2018 Table 4 (final reduced PopPK model).
    # Typical values are for the reference subject: 87.1 kg body weight (median),
    # no diabetes (DIAB = 0), White (RACE_WHITE = 1, equivalently source RACE = 0).
    lcl <- log(0.516); label("Apparent clearance at reference covariates (CL/F, L/day)")  # Yao 2018 Table 4 (CL/F = 0.516 L/day for 87.1 kg, non-diabetic, White)
    lvc <- log(13.5);  label("Apparent volume of distribution at reference covariates (V/F, L)")  # Yao 2018 Table 4 (V/F = 13.5 L for 87.1 kg)
    lka <- log(1.11);  label("First-order SC absorption rate constant (Ka, 1/day)")  # Yao 2018 Table 4 (Ka = 1.11 1/day)

    # Covariate effect exponents and multipliers (Yao 2018 Table 4 footnotes b and c):
    #   CL/F = 0.516 * (BWT/87.1)^0.998 * 1.12^DIAB * 1.11^RACE       (RACE = 1 for non-White)
    #   V/F  = 13.5  * (BWT/87.1)^0.829
    e_wt_cl   <- 0.998; label("Power exponent of body weight on CL/F (unitless)")  # Yao 2018 Table 4 footnote b
    e_wt_vc   <- 0.829; label("Power exponent of body weight on V/F (unitless)")   # Yao 2018 Table 4 footnote c
    e_diab_cl <- 1.12;  label("Multiplier on CL/F for diabetes (1.12^DIAB)")        # Yao 2018 Table 4 footnote b
    e_race_cl <- 1.11;  label("Multiplier on CL/F for non-White race (1.11^(1-RACE_WHITE))")  # Yao 2018 Table 4 footnote b

    # IIV. Yao 2018 Table 4 reports IIV as %CV with explicit footnote
    # "interindividual variability calculated as (variance)^(1/2) x 100%",
    # i.e., IIV/100 is the standard deviation in log-space (omega), so
    # omega^2 = (IIV/100)^2.
    #   IIV CL/F 35.6% -> omega^2 = 0.356^2 = 0.126736
    #   IIV V/F  28.0% -> omega^2 = 0.280^2 = 0.078400
    #   IIV ka  129%   -> omega^2 = 1.290^2 = 1.664100
    # Correlation between IIV of CL/F and V/F is 0.834 (Table 4):
    #   covariance = 0.834 * sqrt(0.126736 * 0.078400) = 0.083133
    etalcl + etalvc ~ c(0.126736,
                        0.083133, 0.078400)  # Yao 2018 Table 4 (IIV CL/F 35.6%, IIV V/F 28.0%, correlation 0.834)
    etalka ~ 1.664100                         # Yao 2018 Table 4 (IIV ka 129%; shrinkage 74.7%)

    # Residual error: combined additive and proportional (Yao 2018 Table 4 and
    # Methods "Base Model" section). Additive component fixed at 0.00289 ug/mL
    # based on a uniform-distribution probability characteristic associated
    # with the LLOQ of 0.01 ug/mL.
    propSd <- 0.200;   label("Proportional residual error (SD, fraction)")  # Yao 2018 Table 4 (Proportional residual error CV% = 20.0)
    addSd  <- 0.00289; label("Additive residual error (SD, ug/mL)")          # Yao 2018 Table 4 (Additive residual error fixed at 0.00289 ug/mL)
  })
  model({
    # Individual PK parameters. Reference subject: 87.1 kg, no diabetes, White.
    # Covariate forms per Yao 2018 Table 4 footnotes b (CL/F) and c (V/F):
    #   CL/F_i = 0.516 * (BWT/87.1)^0.998 * 1.12^DIAB * 1.11^(1 - RACE_WHITE) * exp(eta_CL)
    #   V/F_i  = 13.5  * (BWT/87.1)^0.829                                    * exp(eta_V)
    # The race factor uses (1 - RACE_WHITE) because the canonical RACE_WHITE
    # column inverts Yao's source RACE encoding (see covariateData notes).
    cl <- exp(lcl + etalcl) * (WT / 87.1)^e_wt_cl * e_diab_cl^DIAB * e_race_cl^(1 - RACE_WHITE)
    vc <- exp(lvc + etalvc) * (WT / 87.1)^e_wt_vc
    ka <- exp(lka + etalka)

    kel <- cl / vc

    # One-compartment linear PK with first-order SC absorption and first-order
    # elimination (Yao 2018 Methods "Base Model" and Results "Base Model"
    # section). Bioavailability is implicit in the apparent clearance and
    # volume (CL/F, V/F); SC dosing is administered into the depot compartment.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
