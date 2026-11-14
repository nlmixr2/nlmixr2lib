Tachet_2026_ruxolitinib <- function() {
  description <- "One-compartment population PK model for oral ruxolitinib in adults with myeloproliferative neoplasms, graft-versus-host disease and other haematologic or immune-mediated disorders treated in routine care (Tachet 2026), with first-order absorption (ka fixed), linear elimination, between-subject variability on apparent clearance only, a proportional residual error, body weight on V/F (linear, exponent fixed to 1) and a -39% effect of strong CYP3A or dual CYP2C9/CYP3A inhibitors on CL/F."
  reference <- "Tachet J, Meyer SC, Grandoni F, Haefliger E, Nagy-Hulliger M, Halter JP, Passweg JR, Haefliger D, Bardinet C, Decosterd LA, Guidi M, Girardin FR. Ruxolitinib Pharmacokinetics and Exposure-Toxicity Relationship in Hematologic Malignancies and Immune-Mediated Diseases: A Prospective Observational Study. Clin Pharmacol Ther. 2026;120:753-761. doi:10.1002/cpt.70367"
  vignette <- "Tachet_2026_ruxolitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  dosing <- "depot"

  compartmentData <- list(
    depot = list(analyte = "ruxolitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ruxolitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying in the source analysis (regressor updated at each visit). Enters V/F as (WT/73)^1 with the exponent fixed to 1 (Tachet 2026 Table 2 footnote, 'BWmedian, the population median value of 73 kg'). Table 1 prints the patient-level median as 72.0 kg (range 40.7-116.0); the model equation uses 73 kg, which is reproduced here.",
      source_name = "BW"
    ),
    CONMED_CYP3A4_INH_STRONG = list(
      description = "Concomitant strong CYP3A inhibitor or dual CYP2C9/CYP3A inhibitor at the time of sampling (1 = present, 0 = none, a moderate CYP3A inhibitor only, or status unknown)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no strong CYP3A or dual CYP2C9/CYP3A inhibitor)",
      notes = "Exponential effect on CL/F: exp(-0.50 * CONMED_CYP3A4_INH_STRONG), a 39% reduction (Tachet 2026 Table 2 footnote and Results). The 1 category pools strong CYP3A inhibitors (posaconazole in this cohort) with dual CYP2C9/CYP3A inhibitors (fluconazole, a moderate CYP3A4 inhibitor by the FDA classification); moderate CYP3A inhibitors (isavuconazole, letermovir, verapamil) and samples with missing inhibitor status were pooled into the 0 reference category because their CL estimates were similar to the no-inhibitor group (Tachet 2026 Table 1 footnotes and Supplementary Table S1). Time-varying per sample.",
      source_name = "CYP3A/2C9 inhibitors (strong CYP3A or dual CYP3A/2C9 inhibitor)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL and V/F (Supplementary Table S1: Delta2LL -0.43 on CL, -0.59 on V) and not retained."
    ),
    SEXF = list(
      description = "Biological sex, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Significant on V/F in univariate screening (Delta2LL -5.04, Supplementary Table S1) but dropped in favour of body weight at backward elimination; not significant on CL."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Significant on V/F in univariate screening (Delta2LL -6.63, Supplementary Table S1) but not retained in the final model."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Significant on V/F in univariate screening (Delta2LL -10.8, Supplementary Table S1) but not retained; body weight was the strongest size descriptor."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 77L,
    n_studies = 1L,
    n_observations = 221L,
    age_range = "median 64 years, range 17-89",
    weight_range = "median 72.0 kg, range 40.7-116.0",
    height_range = "median 172 cm, range 146-202",
    sex_female_pct = 40,
    disease_state = "Graft-versus-host disease (acute or chronic) 43%; myeloproliferative neoplasms (myelofibrosis 30%, polycythaemia vera 16%, essential thrombocythaemia 3%, MDS/MPN-RS-T 1%); other haematologic or immune-mediated disease (VEXAS 4%, APECED 1%, HLH 1%)",
    dose_range = "Ruxolitinib 5-25 mg orally once or twice daily, asymmetric twice-daily schedules, or every 48-72 h (Supplementary Table S2); all samples at steady state",
    regions = "Switzerland (University Hospitals of Lausanne, Bern and Basel, Hospital of Morges, private practices)",
    co_medication = "Strong CYP3A (posaconazole) or dual CYP2C9/CYP3A (fluconazole) inhibitors in 24 of 221 samples (11%); moderate CYP3A inhibitors (isavuconazole, letermovir, verapamil) in 8 (4%); inhibitor status missing in 8 (4%)",
    notes = "Prospective observational real-world study, August 2023 to November 2025. 221 steady-state concentrations (12 below the 0.5 ng/mL LLOQ, handled as censored), 87 of them from 11 patients in an intensive sub-study (pre-dose plus 7 samples over 8 h). Demographics from Tachet 2026 Table 1."
  )

  ini({
    lka <- fixed(log(3.24)); label("First-order absorption rate constant ka (1/h)") # Table 2, 'k a (h-1)' = 3.24 FIX
    lvc <- log(69); label("Apparent volume of distribution V/F at 73 kg (L)") # Table 2, 'V pop (L)' = 69 (RSE 5%)
    lcl <- log(14); label("Apparent clearance CL/F without a strong CYP3A or dual CYP2C9/CYP3A inhibitor (L/h)") # Table 2, 'CLpop (L)' = 14 (RSE 6%); unit printed as L, the parameter is a clearance in L/h

    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V/F (unitless)") # Table 2, 'beta V,BW' = 1 FIX
    e_conmed_cyp3a4_inh_strong_cl <- -0.50; label("Log-scale effect of a strong CYP3A or dual CYP2C9/CYP3A inhibitor on CL/F (unitless)") # Table 2, 'beta CL,CYP3A/2C9 inhibitors' = -0.50 (RSE 17%); exp(-0.50) = 0.61, the 39% reduction in Results

    etalcl ~ 0.12825 # Table 2, 'omega CL (CV, %)' = 37; omega^2 = log(1 + 0.37^2)

    propSd <- 0.42; label("Proportional residual error (fraction)") # Table 2, 'sigma prop (CV, %)' = 42 (RSE 6%)
  })

  model({
    ka <- exp(lka)
    # Table 2 footnote: log(Vi) = log(Vpop) + beta_V,BW * log(BWi / 73)
    vc <- exp(lvc) * (WT / 73)^e_wt_vc
    # Table 2 footnote: log(CLi) = log(CLpop) + beta * 1[CYP3A/2C9 inhibitors = 1] + eta_CL,i
    cl <- exp(lcl + e_conmed_cyp3a4_inh_strong_cl * CONMED_CYP3A4_INH_STRONG + etalcl)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # central in mg, vc in L -> mg/L; x 1000 -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
