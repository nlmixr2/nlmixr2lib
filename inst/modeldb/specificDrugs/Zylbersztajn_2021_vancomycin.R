Zylbersztajn_2021_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in 9 paediatric patients (2 months to 11 years) on extracorporeal membrane oxygenation (Zylbersztajn 2021). Pmetrics NPAG non-parametric fit parameterised by clearance, central volume and the rate constants Kcp / Kpc. Total body weight scales clearance and central volume linearly (reference 15.5 kg) and clearance is inversely proportional to serum creatinine (reference 0.34 mg/dL), as CLw = CL*(WT/15.5)*(0.56*(0.34/CREAT)). Typical values are the medians of the per-kg individual-parameter summary (Table 6) anchored at the reference patient."
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
        "Table 4 Pmetrics model: 'CLw = CL*(wgt/15.5)*...' and 'Vw = V*(wgt/15.5)', linear (exponent 1)",
        "scaling of clearance and central volume normalised to 15.5 kg, the vancomycin cohort's mean",
        "weight (Table 2: mean 15.5, median 10, range 3.5-37 kg). Source-paper alias: 'WGT' / 'wgt'."
      ),
      source_name = "WGT"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 4 Pmetrics model: 'CLw = CL*(wgt/15.5)*(0.56*(0.34/crea))', clearance inversely",
        "proportional to serum creatinine normalised to 0.34 mg/dL, the vancomycin cohort's mean",
        "(Table 2: mean 0.34, median 0.22, range 0.08-1.23 mg/dL). The Table 4 footnote gives the unit",
        "(mg/dl). Source-paper alias: 'CREA' / 'crea'."
      ),
      source_name = "CREA"
    )
  )

  covariatesDataExcluded <- list(
    PAGE = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      notes = "Screened (Methods, Covariate Screening); not retained. Table 2: mean 50.4, median 24, range 2-132 months."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Screened; not retained. Table 2: mean 92.7, median 82, range 50-140 cm."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened as gender (male/female); not retained. Table 1: 6 female, 3 male."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Screened; not retained. Table 2: mean 60.1, median 47, range 1.2-248 mg/L."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened; not retained. Table 2 labels the unit g/L but the values (median 3.55, range 2.8-5.4) are g/dL."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened; not retained. Table 2: mean 3.5, median 0.86, range 0.31-30.1 mg/dL."
    )
  )

  compartmentData <- list(
    central = list(analyte = "vancomycin (unbound)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin (unbound)", units = "mg", specimen = "tissue", verified = FALSE)
  )

  population <- list(
    species = "human",
    n_subjects = 9L,
    n_studies = 1L,
    n_profiles = 11L,
    n_concentrations = 44L,
    age_range = "2-132 months",
    age_median = "24 months",
    weight_range = "3.5-37 kg",
    weight_median = "10 kg",
    sex_female_pct = 66.7,
    race_ethnicity = "Not reported (Chilean paediatric cohort)",
    disease_state = paste(
      "Critically ill children on ECMO (5 veno-arterial, 4 veno-venous; polymethylpentene oxygenator",
      "with phosphorylcholine coating) receiving vancomycin for empiric (8) or definitive infection",
      "treatment: pneumonia 7, endocarditis 1, refractory septic shock 1. Two received CRRT."
    ),
    renal_function = "Serum creatinine median 0.22 mg/dL (range 0.08-1.23); estimated creatinine clearance median 178 (30-239) mL/min/1.73 m^2 for the whole study.",
    dose_range = "10-15 mg/kg IV as a 2-h infusion every 6 or 12 h (median 54 mg/kg/day in 2-4 doses), adjusted to renal function.",
    regions = "Chile (Clinica Las Condes and Roberto Del Rio Hospital, Santiago)",
    notes = paste(
      "Demographics from Zylbersztajn 2021 Tables 1-2. Prospective study, July 2017 - November 2019;",
      "neonates excluded. Steady-state sampling from 48 h, at 0.5, 1 and 2 h after the end of infusion",
      "and pre-dose. Vancomycin by immunoenzymatic assay (Cobas c 311), LLOQ 1.7 mg/L, reported by",
      "the authors as unbound. Fitted with Pmetrics 1.5.0 (NPAG). 24-h fluid balance (median 37.5 mL,",
      "range -1203 to 1220) was also screened and not retained."
    )
  )

  ini({
    # Typical values. The Table 4 '#Primary variables' lines are the NPAG
    # search ranges, not estimates; the only numeric results are the
    # Table 6 per-kg summaries of the individual parameters
    # (Mean / SD / Median / CV% / Shrink%). The medians are used and
    # anchored at the Table 4 reference patient (15.5 kg, creatinine
    # 0.34 mg/dL):
    #   Vc  : 0.397 L/kg * 15.5 kg = 6.1535 L (exact under Vw = V*(wgt/15.5))
    #   CLw : 0.061 L/h/kg * 15.5 kg = 0.9455 L/h at the reference patient.
    # The printed CLw equation keeps a literal 0.56 multiplier, so the
    # primary CL that reproduces CLw at the reference is 0.9455 / 0.56 =
    # 1.688393 L/h.
    lcl <- log(1.688393)
    label("Clearance primary CL before the 0.56 factor (L/h)") # Table 6: CL median 0.061 L/h/kg x 15.5 kg / 0.56 (mean 0.060, SD 0.055, CV 92.2%)
    lvc <- log(6.1535)
    label("Central volume Vc at 15.5 kg (L)") # Table 6: Vc median 0.397 L/kg x 15.5 kg (mean 0.419, SD 0.280, CV 61.0%)
    lk12 <- log(0.481)
    label("Central-to-peripheral rate constant Kcp (1/h)") # Table 6: Kcp median 0.481 (mean 0.483, SD 0.304, CV 63.0%)
    lk21 <- log(0.286)
    label("Peripheral-to-central rate constant Kpc (1/h)") # Table 6: Kpc median 0.286 (mean 0.248, SD 0.114, CV 45.8%)

    # Interindividual variability. NPAG estimates a discrete joint density;
    # Table 6 gives only the marginal CV% of each parameter, carried here as
    # independent log-normal random effects with omega^2 = log(CV^2 + 1).
    #   CL  92.2% -> 0.615231
    #   Vc  61.0% -> 0.316342
    #   Kcp 63.0% -> 0.334255
    #   Kpc 45.8% -> 0.190425
    etalcl ~ 0.615231 # Table 6 CL, CV 92.2%
    etalvc ~ 0.316342 # Table 6 Vc, CV 61.0%
    etalk12 ~ 0.334255 # Table 6 Kcp, CV 63.0%
    etalk21 ~ 0.190425 # Table 6 Kpc, CV 45.8%

    # Residual error. Table 4 '#Error model': G = 3 with assay polynomial
    # coefficients 2, 1, 0, 0, i.e. SD = G * (C0 + C1 * Y). Pmetrics adds the
    # two terms linearly, which is nlmixr2's combined1() form. Taken as
    # printed; the resulting 300% proportional term is very wide.
    #   addSd  = 3 * 2 = 6 mg/L
    #   propSd = 3 * 1 = 3
    addSd <- 6
    label("Additive residual SD, gamma * C0 (mg/L)") # Table 4 error model: G = 3, C0 = 2
    propSd <- 3
    label("Proportional residual SD, gamma * C1 (fraction)") # Table 4 error model: G = 3, C1 = 1
  })
  model({
    # Table 4: CLw = CL*(wgt/15.5)*(0.56*(0.34/crea)); Vw = V*(wgt/15.5);
    # Ke = CLw/Vw. The 0.56 is a printed constant of the secondary equation.
    cl <- exp(lcl + etalcl) * (WT / 15.5) * (0.56 * (0.34 / CREAT))
    vc <- exp(lvc + etalvc) * (WT / 15.5)
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

    # Table 4: Y(1) = X(1)/Vw. Dose in mg, Vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
