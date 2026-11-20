Baranda_2024_sar439459 <- function() {
  description <- "Two-compartment population PK model with linear elimination for SAR439459 (anti-TGF-beta monoclonal antibody) given intravenously as monotherapy or with cemiplimab in adults with advanced solid tumours, with body weight on clearance and sex on central volume (Baranda 2024)"
  reference <- paste(
    "Baranda JC, Robbrecht D, Sullivan R, Doger B, Santoro A, Barve M,",
    "Grob J-J, Bechter O, Vieito M, de Miguel MJ, Schadendorf D, Johnson M,",
    "Pouzin C, Cantalloube C, Wang R, Lee J, Chen X, Demers B, Amrate A,",
    "Abbadessa G, Hodi FS. Safety, pharmacokinetics, pharmacodynamics, and",
    "antitumor activity of SAR439459, a TGFbeta inhibitor, as monotherapy and",
    "in combination with cemiplimab in patients with advanced solid tumors:",
    "Findings from a phase 1/1b study. Clin Transl Sci. 2024;17(6):e13854.",
    "doi:10.1111/cts.13854",
    sep = " "
  )
  vignette <- "Baranda_2024_sar439459"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL, (WT/70)^0.994 (Table S2 '+ betaWT'). The source does not state the reference weight; 70 kg is assumed (see the vignette).",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Exponential effect on Vc, Vc x exp(-0.321) for females (Table S2 '+ betaSEX_female'), i.e. 27.5% lower Vc in females.",
      source_name = "SEX"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_CEMIPLIMAB = list(
      description = "Concomitant cemiplimab (Part 1B combination arm)",
      units = "(binary)",
      type = "binary",
      notes = "Tested but not retained: 'Cemiplimab had no effect on SAR'459 PK exposure' (Results, Population pharmacokinetics). No point estimate is reported.",
      source_name = NA_character_
    )
  )

  compartmentData <- list(
    central = list(analyte = "SAR439459", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "SAR439459", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 52L,
    n_studies = 1L,
    age_range = "25-91 years",
    age_median = "60.5 years (Part 1A) and 63.0 years (Part 1B)",
    weight_range = "Not reported",
    weight_median = "Not reported",
    sex_female_pct = 48.1,
    race_ethnicity = "White 86.5%, Black or African American 7.7%, Asian 3.8%, unknown 1.9% (Table 1, Parts 1A and 1B)",
    disease_state = "Advanced unresectable or metastatic solid tumours without a suitable alternative therapy (dose-escalation Parts 1A and 1B of NCT03192345)",
    dose_range = "SAR439459 0.05-15 mg/kg IV Q2W or 22.5 mg/kg IV Q3W, alone (Part 1A) or followed by cemiplimab 3 mg/kg Q2W or 350 mg Q3W (Part 1B)",
    regions = "43 centres in Australia, Belgium, Germany, Spain, Estonia, France, United Kingdom, Italy, Netherlands, Republic of Korea, Taiwan, United States and Canada (whole study)",
    notes = "Model-building dataset: Part 1A (n = 28, monotherapy) and Part 1B (n = 24, with cemiplimab), rich sampling in cycle 1 (pre-dose, 0.5, 2.5, 4.5, 8.5, 24, 48, 72, 96, 168, 336 h, plus 504 h for Q3W) and sparse sampling thereafter. Serum ELISA, LLOQ 0.0078 ug/mL. Fitted by SAEM in Monolix 2019R1. The model was then used to derive exposures for 99 expansion-phase patients (Parts 2A and 2B)."
  )

  ini({
    # Structural parameters: Table S2 (supplement), final population PK model.
    lcl <- log(0.141); label("Clearance for a 70 kg patient (L/day)") # Table S2 'CL (L/day)' = 0.141 (RSE 7%)
    lvc <- log(4.06); label("Central volume of distribution for a male (L)") # Table S2 'Vc (L)' = 4.06 (RSE 5%)
    lq <- log(0.624); label("Intercompartmental clearance (L/day)") # Table S2 'Q (L/day)' = 0.624 (RSE 7%)
    lvp <- log(2.44); label("Peripheral volume of distribution (L)") # Table S2 'Vp (L)' = 2.44 (RSE 13%)

    # Covariate effects (Monolix log-scale coefficients).
    e_wt_cl <- 0.994; label("Power exponent of body weight on CL (unitless)") # Table S2 '+ betaWT' = 0.994 (RSE 30%)
    e_sexf_vc <- -0.321; label("Log-scale effect of female sex on Vc (unitless)") # Table S2 '+ betaSEX_female' = -0.321 (RSE 21%)

    # IIV: Table S2 'omega (%)' column is the Monolix omega, the SD of the
    # log-normal random effect; the variance is omega^2. No correlations reported.
    etalcl ~ 0.113569 # Table S2 omega CL = 33.7% -> 0.337^2
    etalvc ~ 0.050176 # Table S2 omega Vc = 22.4% -> 0.224^2
    etalvp ~ 0.534361 # Table S2 omega Vp = 73.1% -> 0.731^2

    # Residual error: proportional only.
    propSd <- 0.195; label("Proportional residual error (fraction)") # Table S2 'sigma prop (%)' = 19.5 (RSE 3%)
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc + e_sexf_vc * SEXF)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> concentration in mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
