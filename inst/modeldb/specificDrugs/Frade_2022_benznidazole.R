Frade_2022_benznidazole <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "linear elimination for oral benznidazole in Brazilian adults with",
    "chronic Chagas disease, fitted to whole-blood (dried blood spot)",
    "concentrations with Pmetrics NPAG (Frade 2022; n = 8). Clearance is",
    "proportional to body weight normalised to 72.8 kg. Typical values are",
    "the means of the non-parametric support-point distributions; IIV is a",
    "log-normal approximation of the published CV%."
  )
  reference <- paste(
    "Frade VP, Moreira CHV, Sabino EC, Bedor DCG, Ghilard FR, Oliveira CDL,",
    "Sanches C. Population pharmacokinetic modeling of benznidazole in",
    "Brazilian patients with chronic Chagas disease.",
    "Rev Inst Med Trop Sao Paulo. 2022;64:e4.",
    "doi:10.1590/S1678-9946202264004."
  )
  vignette <- "Frade_2022_benznidazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "benznidazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "benznidazole", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters clearance as CL = CL_pop * (WT / 72.8). Frade 2022 Results",
        "state only that 'weight normalized to 72.8 kg' was included as a",
        "covariate on clearance; no exponent appears in Table 2, so the",
        "Pmetrics 'normalized to' idiom is read as a linear-proportional",
        "ratio (exponent 1). Cohort mean weight 70.16 kg (SD 14.20)."
      ),
      source_name = "weight"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Tested on the structural parameters (Frade 2022 Methods, covariate analysis paragraph). Not retained.",
      source_name = "gender"
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on the structural parameters (Frade 2022 Methods). Not retained.",
      source_name = "age"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on the structural parameters (Frade 2022 Methods). Not retained.",
      source_name = "body mass index"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 8L,
    n_studies = 1L,
    n_observations = 8L,
    age_range = "40-60 years (4 aged 40-50, 4 aged 51-60)",
    age_mean = "50.25 years (SD 6.22)",
    weight_mean = "70.16 kg (SD 14.20)",
    sex_female_pct = 75,
    race_ethnicity = c(Mixed = 87.5, Black = 12.5),
    disease_state = paste(
      "Adults with chronic Chagas disease starting standard benznidazole",
      "treatment. HIV infection, renal or hepatic impairment, pregnancy",
      "and lactation were exclusion criteria; non-adherent patients were",
      "excluded."
    ),
    dose_range = "Oral benznidazole 5 mg/kg/day for 60 days.",
    sampling = paste(
      "One whole-blood dried blood spot sample per patient on treatment",
      "day 15 (steady state), 0.8-10.4 h after the dose (Figure 1 VPC),",
      "assayed by LC-MS/MS (Bedor et al.)."
    ),
    regions = "Brazil (Instituto de Infectologia Emilio Ribas and HCFMUSP, Sao Paulo).",
    notes = paste(
      "Prospective cohort interrupted by COVID-19. Demographics from Frade",
      "2022 Table 1 and Results; parameters from Table 2 (final covariate",
      "model)."
    )
  )

  ini({
    # Frade 2022 Table 2, final covariate model; Pmetrics NPAG mean column
    # (median column 1.67 / 6.32 / 37.66 in parentheses of each comment).
    lka <- log(1.66); label("Absorption rate constant Ka (1/h)") # Table 2 Ka mean 1.66 (SD 0.03, median 1.67); Abstract
    lcl <- log(6.27); label("Clearance for a 72.8 kg subject (L/h)") # Table 2 CL mean 6.27 (SD 0.09, median 6.32); Abstract
    lvc <- log(38.97); label("Central volume of distribution (L)") # Table 2 V mean 38.97 (SD 8.33, median 37.66); Abstract

    # Weight effect on clearance, Results: 'inclusion of weight normalized
    # to 72.8 kg as a covariate on clearance'. No exponent is estimated or
    # printed, so the linear-proportional exponent is a structural constant.
    e_wt_cl <- fixed(1); label("Exponent of (WT/72.8) on clearance (unitless)") # Results paragraph 2; exponent not printed, 1 assumed

    # NPAG estimates a non-parametric support-point distribution, not an
    # omega. The Table 2 CV% column is carried as a log-normal
    # approximation, omega^2 = log(CV^2 + 1).
    etalka ~ 0.00039992 # Table 2 Ka CV% 2.0 -> log(0.020^2 + 1)
    etalcl ~ 0.00019598 # Table 2 CL CV% 1.4 -> log(0.014^2 + 1)
    etalvc ~ 0.044778 # Table 2 V CV% 21.4 -> log(0.214^2 + 1)

    # Methods: 'Residual error was modelled as gamma * (1 + 0.1*concentration),
    # value = 5', i.e. Pmetrics SD = gamma * (C0 + C1 * C) with C0 = 1,
    # C1 = 0.1, gamma = 5. addSd = gamma * C0 and propSd = gamma * C1.
    addSd <- 5; label("Additive residual SD (mg/L)") # Methods: gamma 5 x C0 1
    propSd <- 0.5; label("Proportional residual SD (fraction)") # Methods: gamma 5 x C1 0.1
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 72.8)^e_wt_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg, volume in L -> whole-blood concentration in mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
