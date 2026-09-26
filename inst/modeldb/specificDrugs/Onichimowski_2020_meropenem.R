Onichimowski_2020_meropenem <- function() {
  description <- "Two-compartment IV population PK model for meropenem in 19 critically ill adults on continuous renal replacement therapy (CVVH or CVVHD) receiving 1 g as a 1-h infusion q8h (Onichimowski 2020). Elimination is a single first-order clearance attributed to CRRT (ClCRRT); central volume V1 scales with serum albumin by a power function with exponent -2.87 around the 24.6 g/L cohort median. Log-normal IIV on V1, ClCRRT and V2 (IIV on Q fixed to zero); combined additive + proportional residual error."
  reference <- "Onichimowski D, Bedzkowska A, Ziolkowski H, Jaroszewski J, Borys M, Czuczwar M, Wiczling P. Population pharmacokinetics of standard-dose meropenem in critically ill patients on continuous renal replacement therapy: a prospective observational trial. Pharmacol Rep. 2020;72(3):719-729. doi:10.1007/s43440-020-00104-3"
  vignette <- "Onichimowski_2020_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Onichimowski 2020 Table 1: median 24.6 g/L (range 15.6-31.8). Reference value 24.6 g/L (cohort median) in the power relationship V1,i = 27.9 * (ALB_i / 24.6)^-2.87 * exp(eta_V1,i) (Results, final-model equation).",
      source_name = "ALB"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "year",
      type = "continuous",
      notes = "Screened (eta-vs-covariate plots, Methods 'Covariance analysis'); not retained. Table 1 median 67 years (range 36-79)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened; not retained. Table 1 median 80 kg (range 60-100)."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. 5/19 female (Table 1)."
    ),
    DIS_SEPSIS = list(
      description = "Presence of sepsis per Surviving Sepsis Campaign (1 = septic)",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. 10/19 septic (Table 1). The paper notes septic patients had lower albumin (Figure 2, lower panel), so the ALB effect on V1 partly reflects sepsis. Renal function (serum creatinine, eGFR MDRD, Cockcroft-Gault eGFR, diuresis), APACHE II, SOFA, day of therapy and CRRT settings (anticoagulation type, filter day, blood flow, dialysate/substitute flow, UF net) were also screened and not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 19L,
    n_studies = 1L,
    n_concentrations = 256L,
    age_range = "36-79 years",
    age_median = "67 years",
    weight_range = "60-100 kg",
    weight_median = "80 kg",
    sex_female_pct = 26.3,
    race_ethnicity = "Not reported (single-centre Polish ICU cohort)",
    disease_state = "Critically ill adults with acute kidney injury and/or fluid overload on CRRT (9 CVVH with heparin anticoagulation, 10 CVVHD with regional citrate anticoagulation; AV 1000 polysulfone 1.8 m2 filter). 10/19 septic. Median APACHE II 31 (range 8-44), median SOFA 10 (range 4-17).",
    dose_range = "Meropenem 1 g IV as a 1-h infusion q8h (standard licensed dose).",
    regions = "Poland (tertiary medical/surgical ICU, Olsztyn).",
    renal_function = "All on CRRT; median diuresis 0 mL/h (range 0-90); median serum creatinine 1.55 mg/dL (range 0.6-3.7).",
    notes = "Demographics from Onichimowski 2020 Table 1. Arterial sampling at 0, 15, 30, 45, 60, 75, 90, 120, 180, 240 and 480 min after the start of infusion; 9 patients sampled after the first dose and 10 during subsequent days of therapy. Total meropenem by HPLC-UV, LLOQ 0.1 mg/L. NONMEM 7.3 FOCE-I, ADVAN3 TRANS4."
  )

  ini({
    # Structural parameters: Onichimowski 2020 Table 2 final-model 'Estimate' column.
    lvc <- log(27.9); label("Central volume V1 at ALB = 24.6 g/L (L)") # Table 2: theta_V1 = 27.9 L (RSE 17.9%)
    lcl <- log(15.1); label("CRRT clearance ClCRRT (L/h)") # Table 2: theta_ClCRRT = 15.1 L/h (RSE 10.1%)
    lq <- log(21.1); label("Inter-compartmental clearance Q (L/h)") # Table 2: theta_Q = 21.1 L/h (RSE 16.4%)
    lvp <- log(33.7); label("Peripheral volume V2 (L)") # Table 2: theta_V2 = 33.7 L (RSE 28.1%)

    # Covariate effect: V1,i = 27.9 * (ALB_i / 24.6)^-2.87 * exp(eta_V1,i) (Results equation)
    e_alb_vc <- -2.87; label("Power exponent of (ALB / 24.6) on V1 (unitless)") # Table 2: theta_beta,V1 (power function) = -2.87 (RSE 21.4%)

    # IIV: Table 2 reports %CV with footnote %CV = sqrt(exp(omega^2) - 1) * 100,
    # so omega^2 = log(CV^2 + 1). IIV on Q was fixed to 0 and is omitted.
    etalvc ~ 0.24839 # log(0.531^2 + 1); Table 2 omega2 V1 = 53.1 %CV
    etalcl ~ 0.17477 # log(0.437^2 + 1); Table 2 omega2 Cl = 43.7 %CV
    etalvp ~ 0.54971 # log(0.856^2 + 1); Table 2 omega2 V2 = 85.6 %CV

    # Residual error: Cobs = C1 + C1 * eps_prop + eps_add (Methods Eq. 4)
    addSd <- 0.881; label("Additive residual error (mg/L)") # Table 2: sigma add = 0.881 mg/L (RSE 28.4%)
    propSd <- 0.241; label("Proportional residual error (fraction)") # Table 2: sigma2 prop = 24.1 %CV (RSE 10.5%)
  })
  model({
    vc <- exp(lvc + etalvc) * (ALB / 24.6)^e_alb_vc
    cl <- exp(lcl + etalcl)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Methods Eq. 1-2: V1 dC1/dt = -Q C1 + Q C2 - ClCRRT C1; V2 dC2/dt = Q C1 - Q C2
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> mg/L
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
