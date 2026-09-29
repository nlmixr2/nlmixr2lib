Gao_2020_teicoplanin <- function() {
  description <- "Two-compartment IV-infusion population PK model for teicoplanin in 136 Chinese children aged 0.17-9.42 years with different renal functions (Gao 2020). Central volume scales with natural-log body weight, (ln(WT)/2.3)^0.14; peripheral volume with (WT/10)^0.19; clearance with (WT/10)^0.74 and modified-Schwartz eGFR (eGFR/118.99)^0.60. Inter-compartmental clearance has no covariate. Power residual error model with exponent 0.5."
  reference <- "Gao L, Xu H, Ye Q, Li S, Wang J, Mei Y, Niu C, Kang T, Chen C, Wang Y. Population Pharmacokinetics and Dosage Optimization of Teicoplanin in Children With Different Renal Functions. Front Pharmacol. 2020;11:552. doi:10.3389/fphar.2020.00552"
  vignette <- "Gao_2020_teicoplanin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "teicoplanin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "teicoplanin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Gao 2020 Table 1: mean 12.12 kg (SD 6.34), median 10 kg (range 3.5-38). Reference value 10 kg (population median) in Eqs. 16-17; Eq. 15 normalizes ln(WT) by 2.3 (approximately ln(10) = 2.303, i.e. ln of the median weight).",
      source_name = "WT"
    ),
    CRCL = list(
      description = "Modified-Schwartz estimated glomerular filtration rate (BSA-normalized)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Gao 2020 Methods: eGFR (mL/min/1.73 m^2) = 0.413 * height (cm) / serum creatinine (mg/dL) (bedside Schwartz 2009); serum creatinine by enzymatic assay (Roche cobas 8000 c702). Table 1: mean 116.92 (SD 38.45), median 118.99 (range 30.09-280). Reference value 118.99 mL/min/1.73 m^2 (population median) in Eq. 17. Stored under canonical CRCL per inst/references/covariate-columns.md, which accepts a creatinine-based BSA-normalized eGFR.",
      source_name = "eGFR"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Postnatal age",
      units = "years",
      type = "continuous",
      notes = "Gao 2020 Table 2: age-driven maturation (Model II) and age-dependent-exponent (Model IV) clearance models were tested but the simple allometric weight model (Model I) was retained; age is not in the final model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 136L,
    n_studies = 1L,
    age_range = "0.17-9.42 years (text: 0.09-9.42 years; youngest patient 2 months)",
    age_median = "1.25 years (mean 2.19, SD 2.25)",
    weight_range = "3.5-38 kg",
    weight_median = "10 kg (mean 12.12, SD 6.34)",
    sex_female_pct = 41.9,
    race_ethnicity = "Chinese (single center, Wuhan)",
    disease_state = "Children aged 0-10 years with Gram-positive bacterial infection receiving teicoplanin; no neonates.",
    dose_range = "Teicoplanin (Sanofi-Aventis) IV infusion: three loading doses of 10 mg/kg q12h followed by 10 mg/kg once daily maintenance, adjustable per clinical condition.",
    regions = "China (Wuhan Children's Hospital)",
    renal_function = "Modified-Schwartz eGFR median 118.99 mL/min/1.73 m^2 (range 30.09-280): augmented (>= 130) n = 42, normal (90-130) n = 63, mild insufficiency (60-90) n = 23, moderate insufficiency (30-60) n = 8.",
    n_concentrations = 155L,
    notes = "Baseline demographics from Gao 2020 Table 1 (prospective, February 2016 - January 2019). 155 serum concentrations (2.22-79.49 mg/L; 1-3 per patient; 150 at steady state) measured by HPLC-UV (linear range 2.0-180 mg/L). Model fitted in Phoenix NLME 8.1 (FOCE-ELS)."
  )

  ini({
    # Structural parameters (Gao 2020 Table 4 final-model estimates; Eqs. 15-18).
    # Reference subject: WT = 10 kg, eGFR = 118.99 mL/min/1.73 m^2.
    lvc <- log(2.31); label("Central volume of distribution at WT = 10 kg (L)") # Table 4 'theta V1 (L)' = 2.31 (SE 13.31%)
    lvp <- log(16.19); label("Peripheral volume of distribution at WT = 10 kg (L)") # Table 4 'theta V2 (L)' = 16.19 (SE 3.10%)
    lcl <- log(0.13); label("Clearance at WT = 10 kg and eGFR = 118.99 mL/min/1.73 m^2 (L/h)") # Table 4 'theta CL (L/h)' = 0.13 (SE 21.51%)
    lq <- log(0.23); label("Inter-compartmental clearance (L/h)") # Table 4 'theta Q (L/h)' = 0.23 (SE 13.13%)

    # Covariate effects (Eqs. 15-17)
    e_lnwt_vc <- 0.14; label("Power exponent on (ln(WT)/2.3) for central volume (unitless)") # Table 4 'theta 1' = 0.14 (SE 30.69%); Eq. 15
    e_wt_vp <- 0.19; label("Power exponent on (WT/10) for peripheral volume (unitless)") # Table 4 'theta 2' = 0.19 (SE 30.68%); Eq. 16
    e_wt_cl <- 0.74; label("Power exponent on (WT/10) for clearance (unitless)") # Table 4 'theta 3' = 0.74 (SE 29.67%); Eq. 17
    e_crcl_cl <- 0.60; label("Power exponent on (eGFR/118.99) for clearance (unitless)") # Table 4 'theta 4' = 0.60 (SE 30.49%); Eq. 17

    # IIV: Table 4 reports omega as the 'square root of inter-individual
    # variance' in percent (table footnote), so variance = (omega/100)^2.
    etalvc ~ 1.11155 # Table 4 'omega V1 (%)' = 105.43; 1.0543^2
    etalvp ~ 0.038338 # Table 4 'omega V2 (%)' = 19.58; 0.1958^2
    etalcl ~ 0.199541 # Table 4 'omega CL (%)' = 44.67; 0.4467^2
    etalq ~ 0.183698 # Table 4 'omega Q (%)' = 42.86; 0.4286^2

    # Power residual error Y = IPRED + IPRED^power * eps (Eq. 5)
    propSd <- 0.46; label("Power residual error coefficient (mg/L)^(1 - powExp)") # Table 4 'Residual variability sigma' = 0.46 (SE 30.20%)
    powExp <- fixed(0.5); label("Power residual error exponent (unitless)") # Results: 'the power model (Eq. 5) was the best fit with the power value of 0.5'; not in Table 4 so treated as fixed
  })
  model({
    vc <- exp(lvc + etalvc) * (log(WT) / 2.3)^e_lnwt_vc
    vp <- exp(lvp + etalvp) * (WT / 10)^e_wt_vp
    cl <- exp(lcl + etalcl) * (WT / 10)^e_wt_cl * (CRCL / 118.99)^e_crcl_cl
    q <- exp(lq + etalq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, vc in L -> mg/L
    Cc <- central / vc
    Cc ~ pow(propSd, powExp)
  })
}
