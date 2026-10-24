Zylbersztajn_2021_meropenem <- function() {
  description <- "Two-compartment IV population PK model for unbound meropenem in 9 paediatric patients (2 months to 13 years) on extracorporeal membrane oxygenation (Zylbersztajn 2021). Pmetrics NPAG non-parametric fit parameterised by clearance, central volume and the rate constants Kcp / Kpc. Total body weight scales the central volume linearly (reference 20.4 kg); clearance carries no covariate. Typical values are the medians of the per-kg individual-parameter summary (Table 5) anchored at the 20.4 kg reference weight."
  reference <- "Zylbersztajn B, Parker S, Navea D, Izquierdo G, Ortiz P, Torres JP, Fajardo C, Diaz R, Valverde C, Roberts J. Population Pharmacokinetics of Vancomycin and Meropenem in Pediatric Extracorporeal Membrane Oxygenation Support. Front Pharmacol. 2021;12:709332. doi:10.3389/fphar.2021.709332"
  vignette <- "Zylbersztajn_2021_vancomycin_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 3 Pmetrics model: 'Vw = V*(WGT/20.4)', a linear (exponent 1) scaling of the central",
        "volume normalised to 20.4 kg, the meropenem cohort's mean weight (Table 2: mean 20.4,",
        "median 16, range 3.5-45 kg). Weight does not enter clearance ('Ke = CL/Vw').",
        "Source-paper alias: 'WGT'."
      ),
      source_name = "WGT"
    )
  )

  covariatesDataExcluded <- list(
    PAGE = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      notes = "Screened (Methods, Covariate Screening); not retained. Table 2: mean 68.3, median 48, range 2-165 months."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Screened; not retained. Table 2: mean 103, median 100, range 50-155 cm."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened as gender (male/female); not retained. Table 1: 5 female, 4 male."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Screened; not retained. Table 2: mean 90.4, median 57.5, range 1.2-423 mg/L."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened; not retained. Table 2 labels the unit g/L but the values (median 3.5, range 2.6-5.4) are g/dL."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened; not retained. Table 2: mean 3.7, median 0.95, range 0.41-30.1 mg/dL."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened; not retained for meropenem (it is retained on vancomycin clearance). Table 2: mean 0.41, median 0.3, range 0.08-1.27 mg/dL."
    )
  )

  compartmentData <- list(
    central = list(analyte = "meropenem (unbound)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meropenem (unbound)", units = "mg", specimen = "tissue", verified = FALSE)
  )

  population <- list(
    species = "human",
    n_subjects = 9L,
    n_studies = 1L,
    n_profiles = 11L,
    n_concentrations = 33L,
    age_range = "2-165 months",
    age_median = "48 months",
    weight_range = "3.5-45 kg",
    weight_median = "16 kg",
    sex_female_pct = 55.6,
    race_ethnicity = "Not reported (Chilean paediatric cohort)",
    disease_state = paste(
      "Critically ill children on ECMO (5 veno-arterial, 4 veno-venous; polymethylpentene oxygenator",
      "with phosphorylcholine coating) receiving meropenem for empiric (6) or definitive infection",
      "treatment: pneumonia 6, endocarditis 1, refractory septic shock 2. Four received CRRT."
    ),
    renal_function = "Serum creatinine median 0.3 mg/dL (range 0.08-1.27); estimated creatinine clearance median 178 (30-239) mL/min/1.73 m^2 for the whole study.",
    dose_range = "20-40 mg/kg IV as a 3-h infusion every 8 or 12 h (median 120 mg/kg/day in 2-3 doses), adjusted to renal function.",
    regions = "Chile (Clinica Las Condes and Roberto Del Rio Hospital, Santiago)",
    notes = paste(
      "Demographics from Zylbersztajn 2021 Tables 1-2. Prospective study, July 2017 - November 2019;",
      "neonates excluded. Steady-state sampling from 24 h, at 1 and 3 h after the end of infusion and",
      "pre-dose. Unbound plasma meropenem by HPLC-UV, LLOQ 0.16 mg/L. Fitted with Pmetrics 1.5.0",
      "(NPAG). 24-h fluid balance (median 4.5 mL, range -2615 to 1170) was also screened and not retained."
    )
  )

  ini({
    # Typical values. The Table 3 '#Primary variables' lines are the NPAG
    # search ranges, not estimates; the only numeric results are the
    # Table 5 per-kg summaries of the individual parameters
    # (Mean / SD / Median / CV% / Shrink%). The medians are used and
    # anchored at the Table 3 reference weight of 20.4 kg:
    #   Vc : 0.289 L/kg * 20.4 kg = 5.8956 L (exact under Vw = V*(WGT/20.4))
    #   CL : 0.144 L/h/kg * 20.4 kg = 2.9376 L/h, held flat because Table 3
    #        has no weight term on clearance (Ke = CL/Vw)
    lcl <- log(2.9376)
    label("Clearance CL (L/h)") # Table 5: CL median 0.144 L/h/kg x 20.4 kg (mean 0.139, SD 0.102, CV 73.3%)
    lvc <- log(5.8956)
    label("Central volume Vc at 20.4 kg (L)") # Table 5: Vc median 0.289 L/kg x 20.4 kg (mean 0.289, SD 0.295, CV 72.0%)
    lk12 <- log(1.30)
    label("Central-to-peripheral rate constant Kcp (1/h)") # Table 5: Kcp median 1.30 (mean 1.37, SD 0.992, CV 72.6%)
    lk21 <- log(1.34)
    label("Peripheral-to-central rate constant Kpc (1/h)") # Table 5: Kpc median 1.34 (mean 1.25, SD 0.984, CV 79.0%)

    # Interindividual variability. NPAG estimates a discrete joint density;
    # Table 5 gives only the marginal CV% of each parameter, carried here as
    # independent log-normal random effects with omega^2 = log(CV^2 + 1).
    #   CL  73.3% -> 0.430020
    #   Vc  72.0% -> 0.417657
    #   Kcp 72.6% -> 0.423355
    #   Kpc 79.0% -> 0.484954
    etalcl ~ 0.430020 # Table 5 CL, CV 73.3%
    etalvc ~ 0.417657 # Table 5 Vc, CV 72.0%
    etalk12 ~ 0.423355 # Table 5 Kcp, CV 72.6%
    etalk21 ~ 0.484954 # Table 5 Kpc, CV 79.0%

    # Residual error. Table 3 '#Error model': G = 1 with assay polynomial
    # coefficients 0.2, 0.1, 0, 0, i.e. SD = G * (C0 + C1 * Y). Pmetrics adds
    # the two terms linearly, which is nlmixr2's combined1() form.
    #   addSd  = 1 * 0.2 = 0.2 mg/L
    #   propSd = 1 * 0.1 = 0.1
    addSd <- 0.2
    label("Additive residual SD, gamma * C0 (mg/L)") # Table 3 error model: G = 1, C0 = 0.2
    propSd <- 0.1
    label("Proportional residual SD, gamma * C1 (fraction)") # Table 3 error model: G = 1, C1 = 0.1
  })
  model({
    # Table 3: Vw = V*(WGT/20.4); Ke = CL/Vw.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (WT / 20.4)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    kel <- cl / vc

    # Methods: 'Vp was calculated as (Vc*Kcp)/Kpc'. Declaring vp and q keeps
    # rxode2 from matching cl / vc against its one-compartment linear kernel
    # and dropping the peripheral compartment (see Braune_2018_meropenem.R).
    vp <- vc * k12 / k21
    q <- k12 * vc

    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Table 3: Y(1) = X(1)/Vw. Dose in mg, Vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
