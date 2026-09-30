Themans_2020_hydroxychloroquine <- function() {
  description <- "One-compartment first-order-absorption population PK model for oral hydroxychloroquine (HCQ) whole-blood concentrations in 48 hospitalised adult COVID-19 patients from two Belgian tertiary hospitals, with bioavailability fixed at 0.746 from Carmichael 2003 and a power (allometric-type) body-weight effect on clearance (Themans 2020)."
  reference <- "Themans P, Belkhir L, Dauby N, Yombi JC, De Greef J, Delongie KA, Vandeputte M, Nasreddine R, Wittebole X, Wuillaume F, Lescrainier C, Verlinden V, Kiridis S, Dogne JM, Hamdani J, Wallemacq P, Musuamba FT. Population Pharmacokinetics of Hydroxychloroquine in COVID-19 Patients: Implications for Dose Optimization. Eur J Drug Metab Pharmacokinet. 2020;45:703-713. doi:10.1007/s13318-020-00648-y. PMCID PMC7511144."
  vignette <- "Themans_2020_hydroxychloroquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "hydroxychloroquine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "hydroxychloroquine", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on clearance, CL_i = CL * (WT/80)^1.38 * exp(eta_CL) (Eq. 1). The reference weight COV_POP is not printed; 80 kg is the median weight of the model-building dataset (Table 1) and is the WT split used in Figures 2, 6 and S1. The maintainers confirmed it by reproducing Supplementary Table S1 (see vignette). Two patients with missing weight were imputed to the dataset median.",
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on clearance in the stepwise covariate search; did not reach the predefined significance level (Results paragraph 2)."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Significant on clearance in the forward step but dropped because of its correlation with body weight, which gave the larger OFV drop (Results paragraph 2, Fig. 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 48L,
    n_subjects_external_validation = 8L,
    n_observations = 84L,
    age_range = "21-93 years (median 58.5)",
    weight_range = "50-122 kg (median 80)",
    sex_female_pct = 45.8,
    race_ethnicity = "not reported",
    disease_state = "Hospitalised adults with COVID-19 (moderate to severe disease, some in ICU).",
    dose_range = "Belgian national protocol: 400 mg HCQ sulfate twice daily on day 1, then 200 mg twice daily on days 2-5 (Plaquenil 200 mg sulfate = 155 mg HCQ base per tablet).",
    regions = "Belgium (Cliniques Universitaires Saint-Luc and Centre Hospitalier Saint-Pierre, Brussels).",
    notes = "33 patients from an open-label single-arm study (Eudract 2020-001434-35) plus 15 standard-of-care patients formed the model-building dataset (n = 48; 26 male / 22 female); 8 further standard-of-care patients were used for external validation (Table 1). Sparse sampling: one opportunistic sample within 4 h post-dose plus one at end of treatment. Whole-blood HCQ by LC-MS/MS, linear range 100-3000 ng/mL; no sample was outside the quantification range. NONMEM 7.3, FOCE-I."
  )

  ini({
    lka <- log(9.3); label("First-order absorption rate constant (1/h)") # Table 2: Ka = 9.3 /h
    lcl <- log(15.7); label("Apparent clearance at 80 kg (L/h)") # Table 2: CL = 15.7 L/h
    lvc <- log(860.8); label("Apparent volume of distribution (L)") # Table 2: V = 860.8 L
    lfdepot <- fixed(log(0.746)); label("Oral bioavailability (from Carmichael 2003)") # Table 2: F = 0.746 (fixed)
    e_wt_cl <- 1.38; label("Power exponent of body weight on clearance (unitless)") # Table 2: WT effect on CL = 1.38

    etalcl ~ 0.15 # Table 2: IIV on CL (omega^2) = 0.15
    etalvc ~ 0.27 # Table 2: IIV on V (omega^2) = 0.27

    propSd <- 0.1703; label("Proportional residual error (fraction)") # Table 2: Eprop sigma^2 = 0.029; sqrt(0.029) = 0.1703
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 80)^e_wt_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- exp(lfdepot)

    # Dose in mg HCQ base, vc in L -> mg/L; x 1000 -> ng/mL whole blood
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
