Dari_2023_ad26cov2s_selisa <- function() {
  description <- "QSP. Six-compartment mechanistic humoral-immune-response model of spike-protein-binding antibody (S-ELISA, EU/mL) after a single 5 x 10^10 viral-particle dose of the COVID-19 vaccine Ad26.COV2.S in SARS-CoV-2-seronegative adults aged 18-55 and >= 65 years (Dari 2023; final updated model). A K-PD virtual antigen compartment drives logistic memory-B-cell expansion and second-order differentiation into short-lived and long-lived plasma cells, which secrete antibody distributing between serum and a peripheral compartment. Age >= 60 years lowers the antigen-independent memory-B-cell growth rate; female sex raises the long-lived-plasma-cell antibody production rate. The antigen dose is entered in units of 10^8 viral particles (amt = 500 for the 5 x 10^10 vp dose); the paper does not state the dosing scale, and the maintainers selected this reading as the only power-of-ten scaling of the nominal dose that reproduces the Figure 2 VPC median and the Table 3 predicted persistence."
  reference <- paste(
    "Dari A, Boulton M, Neyens M, Le Gars M, Valenzuela B, Shukarev G,",
    "Cardenas V, Ruiz-Guinazu J, Sadoff J, Hoetelmans RMW, Perez Ruixo JJ.",
    "Quantifying Antibody Persistence After a Single Dose of COVID-19 Vaccine",
    "Ad26.COV2.S in Humans Using a Mechanistic Modeling and Simulation Approach.",
    "Clin Pharmacol Ther. 2023;113(2):380. doi:10.1002/cpt.2796.",
    sep = " "
  )
  vignette <- "Dari_2023_ad26cov2s"
  units <- list(
    time = "day",
    dosing = "1e8 viral particles",
    concentration = "EU/mL"
  )

  # Mechanistic immune-response states of the Dari 2023 Figure 1 schematic;
  # none has a canonical compartment name.
  paper_specific_compartments <- c(
    "antigen",
    "bcell_memory",
    "plasma_short",
    "plasma_long",
    "antibody",
    "antibody_periph"
  )

  compartmentData <- list(
    antigen = list(
      analyte = "SARS-CoV-2 spike protein (virtual K-PD antigen)",
      units = "1e8 viral particles",
      specimen = "not applicable",
      verified = TRUE
    ),
    bcell_memory = list(
      analyte = "spike-specific memory B cells",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    plasma_short = list(
      analyte = "spike-specific short-lived plasma cells",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    plasma_long = list(
      analyte = "spike-specific long-lived plasma cells",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    antibody = list(
      analyte = "anti-spike binding IgG (S-ELISA)",
      units = "EU/mL",
      specimen = "serum",
      verified = TRUE
    ),
    antibody_periph = list(
      analyte = "anti-spike binding IgG (S-ELISA)",
      units = "EU/mL",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    AGE_GE60 = list(
      description = "Age >= 60 years indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (< 60 years; the analysed younger cohort was 18-55 years)",
      notes = "Multiplicative fractional effect on the antigen-independent memory-B-cell growth rate kpM: kpM * (1 + e_age_ge60_kp_m * AGE_GE60). Table 1 and Table 2 dichotomise age at 60 years ('decrease for age >= 60 years').",
      source_name = "Age (dichotomic)"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Multiplicative fractional effect on the long-lived-plasma-cell antibody production rate kAbL: kAbL * (1 + e_sexf_kab_l * SEXF).",
      source_name = "Sex"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (power model on kpM, kAbS, kAbL) but not retained in the final model; Methods 'Covariate analysis' and Supplementary Materials section 1.7. Race (White non-Hispanic vs. other) and country were screened the same way and not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 220L,
    n_studies = 2L,
    n_observations = 1185L,
    age_range = "18-82 years (18-55 and >= 65 years enrolled; 56-64 years excluded)",
    age_median = "52.5 years",
    weight_range = "48.0-113 kg",
    weight_median = "72.1 kg",
    sex_female_pct = 46.8,
    race_ethnicity = c(
      `White non-Hispanic` = 86.4,
      `White Hispanic or Latino` = 10.9,
      Asian = 0.9,
      `American Indian or Alaskan Native` = 1.4,
      Other = 0.5
    ),
    disease_state = "Healthy SARS-CoV-2-seronegative adults (no prior infection or COVID-19 vaccination)",
    dose_range = "Ad26.COV2.S 5 x 10^10 viral particles, single dose",
    regions = "United States, Belgium, the Netherlands, Spain, Germany",
    notes = "Pooled COV1001 cohorts 1a and 3 (NCT04436276) and COV2001 group 5 (NCT04535453); sampling on days 1, 15, 29, 57, (64,) 71, 85, 169 and 239. Demographics from Table 1 (S-ELISA column). LLOQ 50.3 EU/mL, ULOQ 15,797.9 EU/mL."
  )

  ini({
    # Antigen (virtual K-PD compartment)
    k_ag <- fixed(0.065)
    label("Apparent decay rate of the vaccine-generated S-protein antigen (1/day)") # Table 2 'Ad26 antigen decay, kAg' = 0.065, fixed (Balelli 2020)

    # Memory B cells
    kp_ag <- fixed(1)
    label("Antigen-dependent memory-B-cell production rate (cells per antigen unit per day)") # Table 2 'kpAg' = 1, fixed for nonidentifiability
    lkp_m <- log(0.523)
    label("Antigen-independent memory-B-cell growth rate kpM, age < 60 years (1/day)") # Table 2 'kpM for age < 60 years' = 0.523 (95% CI 0.381-0.665)
    ke_m <- fixed(0)
    label("Memory-B-cell elimination rate (1/day)") # Table 2 'Mem B decay, keM' = 0 /year, fixed due to short follow-up
    n_mem <- fixed(3.019e4)
    label("Maximum memory-B-cell count for logistic growth (cells)") # Table 2 'Maximum mem B produced (N)' = 3.019 x 10^4, fixed (Goel 2021)

    # Plasma-cell differentiation and loss
    lkp_spl <- log(0.691e-5)
    label("Sum of short- and long-lived plasma-cell differentiation rates kpS + kpL (1/day per antigen unit)") # Table 2 'kpSpL = kpS + kpL' = 0.691 x 10^-5 (95% CI 0.517-0.865 x 10^-5)
    r_kp_s_l <- fixed(2.5 / 0.011)
    label("Ratio of short- to long-lived plasma-cell differentiation rates kpS/kpL (unitless)") # Table 2 'kpS/kpL' = 2.5/0.011, fixed (Balelli 2020)
    ke_s <- fixed(0.578)
    label("Short-lived plasma-cell elimination rate (1/day)") # Table 2 'S elimination, keS' = 0.578 /day, fixed (Pasin 2019)
    lke_l <- log(16.8e-4)
    label("Long-lived plasma-cell elimination rate (1/day)") # Table 2 'keL' S-ELISA final updated model = 16.8 x 10^-4 (95% CI 9.45-24.15 x 10^-4); 0 fixed in the day-85 final model

    # Antibody production and disposition
    lkab_s <- log(0.530)
    label("Antibody production rate per short-lived plasma cell kAbS ((EU/mL)/cell/day)") # Table 2 'kAbS' S-ELISA = 0.530 (95% CI 0.377-0.683)
    lkab_l <- log(4.43)
    label("Antibody production rate per long-lived plasma cell kAbL, males ((EU/mL)/cell/day)") # Table 2 'kAbL for males' S-ELISA = 4.43 (95% CI 3.70-5.16)
    kab_p <- fixed(0.2075)
    label("Serum-to-peripheral antibody distribution rate kAbP (1/day)") # Table 2 'kAbP' = 0.2075, fixed (Davda 2014)
    kp_ab <- fixed(0.2716)
    label("Peripheral-to-serum antibody return rate kPAb (1/day)") # Table 2 'kPAb' = 0.2716, fixed (Davda 2014)
    ke_ab <- fixed(0.0556)
    label("Serum antibody elimination rate keAb (1/day)") # Table 2 'Ab decay rate, keAb' = 0.0556, fixed (Davda 2014)

    # Covariate effects, form theta * (1 + estimate * indicator) (Supplementary Materials 1.7)
    e_age_ge60_kp_m <- -0.441
    label("Fractional change in kpM for age >= 60 years (unitless)") # Table 2 'kpM_AGE' S-ELISA = -0.441 (95% CI -0.614 to -0.268)
    e_sexf_kab_l <- 0.671
    label("Fractional change in kAbL for females (unitless)") # Table 2 'kAbL_SEX' S-ELISA = 0.671 (95% CI 0.140-1.20)

    # IIV, final updated model; omega^2 = log(1 + CV^2) from Table 2 CV% (footnote a)
    etalkp_m ~ 0.4801 # Table 2 IIV on kpM S-ELISA final updated CV 78.5% (100% in the day-85 final model)
    etalkab_s ~ 1.1509 # Table 2 IIV on kAbS S-ELISA final updated CV 147% (145% in the day-85 final model)
    etalkab_l ~ 1.0946 # Table 2 IIV on kAbL S-ELISA final updated CV 141% (162% in the day-85 final model)

    # Residual error: additive on log10 titers (Supplementary Materials 1.6). The
    # Table 2 'Proportional error' entry is the log10-scale SD x 100, so the
    # natural-log SD is 0.117 * log(10) = 0.2694. Reading it instead as a CV on
    # the original scale (natural-log SD 0.117) leaves the simulated 24-month
    # persistence of Table 3 too flat; see the vignette.
    expSd <- 0.2694
    label("Log-scale residual SD for S-ELISA titers (natural-log units)") # Table 2 'Proportional error' S-ELISA final updated 11.7 (log10 SD 0.117; 8.28 in the day-85 final model)
  })

  model({
    # Individual parameters
    kp_m <- exp(lkp_m + etalkp_m) * (1 + e_age_ge60_kp_m * AGE_GE60)
    kab_s <- exp(lkab_s + etalkab_s)
    kab_l <- exp(lkab_l + etalkab_l) * (1 + e_sexf_kab_l * SEXF)
    ke_l <- exp(lke_l)

    # kpS and kpL from their estimated sum and fixed ratio (Table 2 'Derived' rows)
    kp_spl <- exp(lkp_spl)
    kp_s <- kp_spl / (1 + 1 / r_kp_s_l)
    kp_l <- kp_spl - kp_s

    # Eq. 1-6. Seronegative participants: every state starts at zero and the
    # vaccine dose is a bolus into antigen.
    d/dt(antigen) <- -k_ag * antigen
    d/dt(bcell_memory) <- kp_ag * antigen + kp_m * bcell_memory * (1 - bcell_memory / n_mem) -
      ((kp_s + kp_l) * antigen + ke_m) * bcell_memory
    d/dt(plasma_short) <- kp_s * antigen * bcell_memory - ke_s * plasma_short
    d/dt(plasma_long) <- kp_l * antigen * bcell_memory - ke_l * plasma_long
    d/dt(antibody) <- kab_s * plasma_short + kab_l * plasma_long - ke_ab * antibody -
      kab_p * antibody + kp_ab * antibody_periph
    d/dt(antibody_periph) <- kab_p * antibody - kp_ab * antibody_periph

    # Observed serum binding-antibody concentration (EU/mL). Named Cc, as for the
    # anti-RBD IgG output of Dasti_2025_mrna1273_qsp, although it is an
    # antibody response rather than a concentration of the administered vaccine.
    Cc <- antibody
    Cc ~ lnorm(expSd)
  })
}
