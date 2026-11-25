Kreuer_2014_propofol_lungcmt <- function() {
  description <- paste(
    "Exhaled-breath propofol model, lung-compartment approach (Kreuer 2014,",
    "approach 1; one adult surgical patient). The Marsh three-compartment",
    "intravenous propofol model (micro-constants and weight-proportional",
    "central volume fixed to the Marsh data set) drives a lung concentration",
    "with separate in and out rate constants, dlung/dt = k1L * Cc - kL1 *",
    "lung (k1L = 0.209, kL1 = 0.138 /min), fitted in Excel to one patient's",
    "breath time series. The lung state is driven by plasma and does not",
    "deplete it. Measured breath propofol x relates to the lung",
    "concentration y (ug/mL) by a zero-intercept second-order polynomial,",
    "y = q * x^2 + s * x, fitted separately to the breath concentration in",
    "ppb (Cbreath) and to the raw ion-mobility-spectrometer peak height in",
    "volt (Cbreath_volt); the outputs invert that polynomial (lower root,",
    "held at the vertex above the polynomial maximum). The printed rate",
    "constants give a 4.54 ug/mL lung plateau at the 3 ug/mL plasma target,",
    "above the ppb polynomial maximum (3.61) and the ~3.6 ug/mL plateau of",
    "the paper's Fig 6b, so Cbreath sits at its 0.78 ppb cap through",
    "maintenance (see the vignette). Deterministic typical-value model: no",
    "inter-individual variability and no residual-error model were reported."
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
    notes = "Single-patient demonstration. Exhaled propofol measured by a multi-capillary-column ion mobility spectrometer (BioScout, 10 mL breath per sample, one value per minute), reported as peak height (V) and as ppb via a calibration curve. Parameters were fitted in Excel 2010 by ordinary least squares, maximising R2 with the Solver tool. The paper also reports the lung compartment as a volume VLung = 17 L and a clearance ClLung = 3.75 L/min; neither formula printed in the Fig 2 legend maps those onto the printed rate constants, so the model carries the rate constants, which are what the printed update equation uses."
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

    # Lung compartment (fitted to the ppb breath series)
    lkin_lung <- log(0.209); label("Plasma -> lung rate constant k1L (1/min)") # Results 'Propofol lung compartment concentration': 'k1L = 0,209 min-1'
    lkout_lung <- log(0.138); label("Lung loss rate constant kL1 (1/min)") # Results 'Propofol lung compartment concentration': 'kL1 = 0,138 min-1'

    # Breath -> lung polynomials, y = q * x^2 + s * x + 0 (y = lung ug/mL)
    cal_slope_breath <- 9.23; label("Linear coefficient s of the ppb breath polynomial (ug/mL per ppb)") # Results: 'y = -5.9*x2 + 9.23*x + 0' (ppb, R2 = 0.99, Fig 6b)
    cal_slope_breath_quad <- -5.9; label("Quadratic coefficient q of the ppb breath polynomial (ug/mL per ppb^2)") # Results: 'y = -5.9*x2 + 9.23*x + 0' (ppb, R2 = 0.99, Fig 6b)
    cal_slope_breathv <- 15.57; label("Linear coefficient s of the peak-height polynomial (ug/mL per V)") # Results: 'y = -12.6*x2 + 15.57*x + 0' (volt, R2 = 0.87, Fig 6a)
    cal_slope_breathv_quad <- -12.6; label("Quadratic coefficient q of the peak-height polynomial (ug/mL per V^2)") # Results: 'y = -12.6*x2 + 15.57*x + 0' (volt, R2 = 0.87, Fig 6a)
  })

  model({
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13)
    k31 <- exp(lk31)
    vc <- exp(lvc) * WT
    kin_lung <- exp(lkin_lung)
    kout_lung <- exp(lkout_lung)

    d/dt(central) <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Calculated plasma concentration (mg/L = ug/mL)
    Cc <- central / vc

    # Lung concentration (Fig 2 legend, 'LCC = (k1L*PPCO-kL1*LCCO)+LCCO', a
    # one-minute Euler step of this ODE). The lung runs in parallel with V1
    # (Discussion: 'The first compartment V1 is used parallel for both
    # models'), so it is driven by Cc and removes nothing from central.
    d/dt(lung) <- kin_lung * Cc - kout_lung * lung

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
