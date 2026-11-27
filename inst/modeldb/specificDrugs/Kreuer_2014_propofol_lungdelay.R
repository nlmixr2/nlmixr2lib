Kreuer_2014_propofol_lungdelay <- function() {
  description <- paste(
    "Exhaled-breath propofol model, time-delay approach (Kreuer 2014,",
    "approach 2; one adult surgical patient). The Marsh three-compartment",
    "intravenous propofol model (micro-constants and weight-proportional",
    "central volume fixed to the Marsh data set) drives a lung concentration",
    "through a first-order equilibration delay, dlung/dt = ke0lung * (Cc -",
    "lung), with ke0lung = 0.27 /min fitted in Excel to one patient's",
    "breath time series. Measured breath propofol x relates to the lung",
    "concentration y (ug/mL) by a zero-intercept second-order polynomial,",
    "y = q * x^2 + s * x, fitted separately to the breath concentration in",
    "ppb (Cbreath) and to the raw ion-mobility-spectrometer peak height in",
    "volt (Cbreath_volt); the outputs invert that polynomial (lower root,",
    "held at the vertex above the polynomial maximum). Deterministic",
    "typical-value model: one patient, no inter-individual variability and",
    "no residual-error model were reported."
  )
  reference <- paste(
    "Kreuer S, Hauschild A, Fink T, Baumbach JI, Maddula S, Volk T.",
    "Two different approaches for pharmacokinetic modeling of exhaled drug",
    "concentrations. Sci Rep. 2014;4:5423. doi:10.1038/srep05423.",
    "PMCID: PMC4067807.",
    "Plasma PK is the Marsh model as printed in Kreuer 2014 Table 1 (Marsh B,",
    "White M, Morton N, Kenny GN. Br J Anaesth. 1991;67:41-48.",
    "doi:10.1093/bja/67.1.41)."
  )
  vignette <- "Kreuer_2014_propofol"
  units <- list(
    time = "min",
    dosing = "mg",
    concentration = "ug/mL (plasma and lung); ppb (Cbreath); V (Cbreath_volt)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters only through the Marsh central volume, V1 = 0.228 L/kg * WT (Kreuer 2014 Table 1, '0,228*m [kg]'; linear, no reference weight). The studied patient's weight is not reported. Under target-controlled infusion the plasma, lung and breath predictions do not depend on WT; only the infused amount does.",
      source_name = "m"
    )
  )

  compartmentData <- list(
    central = list(analyte = "propofol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "propofol", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "propofol", units = "mg", specimen = "tissue", verified = TRUE),
    lung = list(analyte = "propofol", units = "ug/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1,
    n_studies = 1,
    age_range = "adult (age not reported)",
    weight_range = "not reported",
    sex_female_pct = NA,
    disease_state = "Adult patient scheduled for an abdominal surgical procedure under total intravenous anaesthesia",
    dose_range = "Propofol target-controlled infusion (Marsh, plasma target 3 ug/mL held constant) with remifentanil (Minto, plasma mode), stopped at the end of surgery; breath sampled every minute for 346 min",
    regions = "Germany (Saarland University Medical Center, Homburg)",
    notes = "Single-patient demonstration. Exhaled propofol measured by a multi-capillary-column ion mobility spectrometer (BioScout, 10 mL breath per sample, one value per minute), reported as peak height (V) and as ppb via a calibration curve. Parameters were fitted in Excel 2010 by ordinary least squares, maximising R2 with the Solver tool."
  )

  ini({
    # Marsh three-compartment PK, fixed (Kreuer 2014 Table 1, the 'Marsh data set'
    # used by the TCI pump; not estimated in this paper)
    lkel <- fixed(log(0.1190)); label("Elimination rate constant k10 (1/min)") # Table 1: k10 = 0,1190
    lk12 <- fixed(log(0.1140)); label("Rate constant central -> peripheral1 k12 (1/min)") # Table 1: k12 = 0,1140
    lk21 <- fixed(log(0.0550)); label("Rate constant peripheral1 -> central k21 (1/min)") # Table 1: k21 = 0,0550
    lk13 <- fixed(log(0.0419)); label("Rate constant central -> peripheral2 k13 (1/min)") # Table 1: k13 = 0,0419
    lk31 <- fixed(log(0.0033)); label("Rate constant peripheral2 -> central k31 (1/min)") # Table 1: k31 = 0,0033
    lvc <- fixed(log(0.228)); label("Central volume per kg body weight V1 (L/kg)") # Table 1: V1 = 0,228*m [kg]

    # Lung time delay (fitted to the ppb breath series)
    lke0 <- log(0.27); label("Plasma-lung equilibration rate constant ke0lung (1/min)") # Results 'Time delay of the propofol plasma concentration': ke0lung = 0.27 min-1 (also Abstract)

    # Breath -> lung polynomials, y = q * x^2 + s * x + 0 (y = lung ug/mL)
    cal_slope_breath <- 7.14; label("Linear coefficient s of the ppb breath polynomial (ug/mL per ppb)") # Results: 'y = -4.2*x2 + 7.14*x + 0' (ppb, R2 = 0.96, Fig 7b)
    cal_slope_breath_quad <- -4.2; label("Quadratic coefficient q of the ppb breath polynomial (ug/mL per ppb^2)") # Results: 'y = -4.2*x2 + 7.14*x + 0' (ppb, R2 = 0.96, Fig 7b)
    cal_slope_breathv <- 9.64; label("Linear coefficient s of the peak-height polynomial (ug/mL per V)") # Results: 'y = -7.7*x2 + 9.64*x + 0' (volt, R2 = 0.82, Fig 7a)
    cal_slope_breathv_quad <- -7.7; label("Quadratic coefficient q of the peak-height polynomial (ug/mL per V^2)") # Results: 'y = -7.7*x2 + 9.64*x + 0' (volt, R2 = 0.82, Fig 7a)
  })

  model({
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13)
    k31 <- exp(lk31)
    vc <- exp(lvc) * WT
    ke0 <- exp(lke0)

    d/dt(central) <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Calculated plasma concentration (mg/L = ug/mL)
    Cc <- central / vc

    # Lung concentration lags plasma by a first-order delay (Fig 3). The
    # printed update 'LCC = (PPC+ (LCCO-PPC))*exp(-ke0lung)' misplaces a
    # parenthesis (as printed it reduces to LCCO*exp(-ke0lung) and never
    # rises from zero); the exact one-minute step of this ODE is
    # PPC + (LCCO - PPC) * exp(-ke0lung), which is what the text describes.
    d/dt(lung) <- ke0 * (Cc - lung)

    # Breath outputs: invert lung = q * x^2 + s * x (q < 0) on its rising
    # branch, x = (sqrt(s^2 + 4 * q * lung) - s) / (2 * q). Above the vertex
    # lung = -s^2 / (4 * q) the polynomial has no real root; the breath value
    # is held at the vertex x = -s / (2 * q), as the fitted curves in Fig 6 and
    # Fig 7 are drawn flat there.
    disc_breath <- cal_slope_breath^2 + 4 * cal_slope_breath_quad * lung
    if (disc_breath < 0) {
      disc_breath <- 0
    }
    Cbreath <- (sqrt(disc_breath) - cal_slope_breath) / (2 * cal_slope_breath_quad)

    disc_breathv <- cal_slope_breathv^2 + 4 * cal_slope_breathv_quad * lung
    if (disc_breathv < 0) {
      disc_breathv <- 0
    }
    Cbreath_volt <- (sqrt(disc_breathv) - cal_slope_breathv) / (2 * cal_slope_breathv_quad)

    # No residual-error model: one patient fitted by ordinary least squares
    # in Excel; no error magnitude was reported.
  })
}
