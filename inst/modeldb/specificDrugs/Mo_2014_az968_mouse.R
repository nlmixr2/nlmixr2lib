Mo_2014_az968_mouse <- function() {
  description <- "Preclinical (mouse, HCT116 human colon carcinoma xenograft). Lifespan-based tumor growth inhibition (LS TGI) model of the casein kinase 2 inhibitor AZ968 with a non-cycle-specific drug effect, in the linear-growth form (Mo 2014 Model Comparison case study, Tables 3 and 4). Tumor cells have a common lifespan td_div at the end of which they divide; because all observed tumors were above the size threshold, the division efficiency is p(w) = 1 + pwth / w (Eq. 16), where pwth = (p0 - 1) * wth is the identifiable product. The tumor production rate is the recursive lifespan product kin(t) = kin0 * prod over past lifespans of p(w(t - i td_div)), so the model is a delay differential equation with one delayed tumor volume per elapsed lifespan; these are encoded with rxode2 delay() for 40 lifespans, which makes the model valid for t < 40.5 * td_div (52 days at the published td_div). AZ968 removes proliferating cells at the second-order rate k2 * C(t) (Eq. 18) into a non-proliferating apoptotic pool that dies after td_apo (Eq. 24); the surviving fraction of the proliferating cells, exp(-integral of k2 * C), is carried as the cumulative kill hazard cumhaz. Following the authors' published MATLAB implementation, the number of past lifespans in each product is round(t / td_div), not the floor INT(t / td_div) printed in Equations 7 and 22 (see the vignette). Two-compartment AZ968 PK with first-order absorption from the intraperitoneal injection site and dose-specific elimination rate constant and central volume (non-linear PK; Table 3), selected by the DOSE_AZ968_MGKG covariate. Time is in days, doses in mg/kg into depot, plasma concentration in ng/mL and tumor volume in mm^3. Solve with method = 'dop853' and maxsteps = 1e7."
  reference <- paste(
    "Mo G, Gibbons F, Schroeder P, Krzyzanski W. Lifespan based",
    "pharmacokinetic-pharmacodynamic model of tumor growth inhibition by",
    "anticancer therapeutics. PLoS One. 2014;9(10):e109747.",
    "doi:10.1371/journal.pone.0109747"
  )
  vignette <- "Mo_2014_lifespan_tgi"
  units <- list(time = "day", dosing = "mg/kg", concentration = "ng/mL")

  covariateData <- list(
    DOSE_AZ968_MGKG = list(
      description = "AZ968 dose per administration of the animal's treatment group",
      units = "mg/kg",
      type = "continuous",
      reference_category = "10 mg/kg",
      notes = "Selects the dose-specific PK parameters of Mo 2014 Table 3: 20 mg/kg and 30 mg/kg use their own kel and V; every other value (including 10 mg/kg and the untreated controls) takes the 10 mg/kg values. The dose-specific values were fitted to single 10, 20 and 30 mg/kg intraperitoneal doses and are not an exposure-response; do not interpolate to other doses.",
      source_name = "dose group (10, 20 or 30 mg/kg)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "AZ968", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "AZ968", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "AZ968", units = "mg/kg", specimen = "tissue", verified = TRUE),
    cycling_cells = list(analyte = "proliferating tumor cells M", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(
      analyte = "apoptotic (non-proliferating) tumor cells A",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    cumhaz = list(
      analyte = "cumulative drug-induced kill hazard of proliferating cells",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (female NCr nude, HCT116 xenograft; PK in female CD1 mice)",
    n_subjects = NA_integer_,
    n_studies = 2L,
    sex_female_pct = 100,
    age_range = "5-6 weeks",
    disease_state = "subcutaneous HCT116 human colon carcinoma xenograft, treatment started at 150-200 mm^3",
    dose_range = "AZ968 10, 20 or 30 mg/kg intraperitoneal once weekly for three weeks (tumor study, Figure 8B); single 10, 20 or 30 mg/kg intraperitoneal doses (PK study, 3 mice per dose)",
    notes = "Raw data supplied by AstraZeneca (Mo 2014 Methods, Animal Data; Schroeder 2012 AACR abstract). The Methods give 10 to 15 mice per treatment group for the efficacy study and three mice per dose for the pharmacodynamic study. Mean data were fitted with MATLAB nlinfit; td_div and td_apo were found by grid search and then held fixed. Initial tumor volume 180 mm^3 (Figure 8 caption)."
  )

  ini({
    # Lifespan TGI parameters: Mo 2014 Table 4 (left, Lifespan TGI Model).
    ltd_div <- fixed(log(1.28)); label("Tumor cell lifespan (doubling time) T (day)") # Table 4: T = 1.28 day, grid search
    ltd_apo <- fixed(log(5.56)); label("Apoptosis duration TA (day)") # Table 4: TA = 5.56 day, grid search
    lpwth <- log(80.8); label("Linear-phase division-efficiency product pwth = (p0 - 1) wth (mm^3)") # Table 4: Pwth = 80.8 mm^3
    lkin0 <- log(49.9); label("Tumor production rate before the experiment kin0 (mm^3/day)") # Table 4: kin0 = 49.9 mm^3/day
    ldrugSlope <- log(2.30e-3); label("Second-order drug potency constant k2 (mL/ng/day)") # Table 4: k2 = 2.30e-3 mL/ng/day
    lrbase_tumor <- fixed(log(180)); label("Tumor volume at the first dose w0 (mm^3)") # Figure 8 caption: initial tumor volume fixed at 180 mm^3

    # AZ968 PK: Mo 2014 Table 3 (two-compartment, first-order absorption,
    # dose-specific kel and V).
    lka <- fixed(log(46.8)); label("AZ968 absorption rate constant ka (1/day)") # Table 3: ka = 46.8 1/day, held constant
    lkel <- log(33.9); label("AZ968 elimination rate constant at 10 mg/kg kel (1/day)") # Table 3: kel 10 mg/kg = 33.9 1/day
    lkel_dose20 <- log(26.0); label("AZ968 elimination rate constant at 20 mg/kg (1/day)") # Table 3: kel 20 mg/kg = 26.0 1/day
    lkel_dose30 <- log(24.7); label("AZ968 elimination rate constant at 30 mg/kg (1/day)") # Table 3: kel 30 mg/kg = 24.7 1/day
    lk12 <- log(8.74); label("AZ968 central-to-peripheral rate constant k12 (1/day)") # Table 3: k12 = 8.74 1/day
    lk21 <- log(11.2); label("AZ968 peripheral-to-central rate constant k21 (1/day)") # Table 3: k21 = 11.2 1/day
    lvc <- log(2.70); label("AZ968 central volume at 10 mg/kg V (L/kg)") # Table 3: V 10 mg/kg = 2.70 L/kg
    lvc_dose20 <- log(2.07); label("AZ968 central volume at 20 mg/kg (L/kg)") # Table 3: V 20 mg/kg = 2.07 L/kg
    lvc_dose30 <- log(1.97); label("AZ968 central volume at 30 mg/kg (L/kg)") # Table 3: V 30 mg/kg = 1.97 L/kg

    # Mo 2014 fitted the mean profiles by unweighted nonlinear least
    # squares (MATLAB nlinfit) and reports no residual-error estimate.
    addSd_tumor_vol <- fixed(0); label("Additive residual error on tumor volume (mm^3; not reported)") # not reported (nlinfit least squares)
  })

  model({
    td_div <- exp(ltd_div)
    td_apo <- exp(ltd_apo)
    pwth <- exp(lpwth)
    kin0 <- exp(lkin0)
    drugSlope <- exp(ldrugSlope)
    rbase <- exp(lrbase_tumor)

    # Dose-specific PK (Table 3); any dose other than 20 or 30 mg/kg takes
    # the 10 mg/kg values.
    kel <- exp(lkel)
    vc <- exp(lvc)
    if (DOSE_AZ968_MGKG == 20) {
      kel <- exp(lkel_dose20)
      vc <- exp(lvc_dose20)
    }
    if (DOSE_AZ968_MGKG == 30) {
      kel <- exp(lkel_dose30)
      vc <- exp(lvc_dose30)
    }
    ka <- exp(lka)
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    # AZ968 PK (dose in mg/kg into depot; mg/L x 1000 = ng/mL)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    Cc <- central / vc * 1000

    # Tumor volume w = M + A (Eq. 17); linear-phase efficiency p(w) (Eq. 16)
    tumor_vol <- cycling_cells + damaged_cells1
    pw <- 1 + pwth / tumor_vol
    # Non-cycle-specific linear drug effect (Eq. 18)
    eff <- drugSlope * Cc

    # Delayed tumor volume for i = 1..40 past lifespans: w_di = w(t - i T),
    # p_di = p(w_di). Before the experiment (t - i T < 0) the tumor is at
    # w0 with A = 0 (Eqs. 23, 25 and the supplement history function S6).
    # As in the supplement S5 code, p_di enters once round(t / T) >= i,
    # that is from t = (i - 0.5) T, and is then p(w0) until t = i T; the
    # code additionally holds p_d1 out until t = T.
    w_d1 <- rbase
    if (t >= td_div) {w_d1 <- delay(cycling_cells, td_div) + delay(damaged_cells1, td_div)}
    p_d1 <- 1; if (t >= td_div) {p_d1 <- 1 + pwth / w_d1}
    w_d2 <- rbase
    if (t >= 2 * td_div) {w_d2 <- delay(cycling_cells, 2 * td_div) + delay(damaged_cells1, 2 * td_div)}
    p_d2 <- 1; if (t >= 1.5 * td_div) {p_d2 <- 1 + pwth / w_d2}
    w_d3 <- rbase
    if (t >= 3 * td_div) {w_d3 <- delay(cycling_cells, 3 * td_div) + delay(damaged_cells1, 3 * td_div)}
    p_d3 <- 1; if (t >= 2.5 * td_div) {p_d3 <- 1 + pwth / w_d3}
    w_d4 <- rbase
    if (t >= 4 * td_div) {w_d4 <- delay(cycling_cells, 4 * td_div) + delay(damaged_cells1, 4 * td_div)}
    p_d4 <- 1; if (t >= 3.5 * td_div) {p_d4 <- 1 + pwth / w_d4}
    w_d5 <- rbase
    if (t >= 5 * td_div) {w_d5 <- delay(cycling_cells, 5 * td_div) + delay(damaged_cells1, 5 * td_div)}
    p_d5 <- 1; if (t >= 4.5 * td_div) {p_d5 <- 1 + pwth / w_d5}
    w_d6 <- rbase
    if (t >= 6 * td_div) {w_d6 <- delay(cycling_cells, 6 * td_div) + delay(damaged_cells1, 6 * td_div)}
    p_d6 <- 1; if (t >= 5.5 * td_div) {p_d6 <- 1 + pwth / w_d6}
    w_d7 <- rbase
    if (t >= 7 * td_div) {w_d7 <- delay(cycling_cells, 7 * td_div) + delay(damaged_cells1, 7 * td_div)}
    p_d7 <- 1; if (t >= 6.5 * td_div) {p_d7 <- 1 + pwth / w_d7}
    w_d8 <- rbase
    if (t >= 8 * td_div) {w_d8 <- delay(cycling_cells, 8 * td_div) + delay(damaged_cells1, 8 * td_div)}
    p_d8 <- 1; if (t >= 7.5 * td_div) {p_d8 <- 1 + pwth / w_d8}
    w_d9 <- rbase
    if (t >= 9 * td_div) {w_d9 <- delay(cycling_cells, 9 * td_div) + delay(damaged_cells1, 9 * td_div)}
    p_d9 <- 1; if (t >= 8.5 * td_div) {p_d9 <- 1 + pwth / w_d9}
    w_d10 <- rbase
    if (t >= 10 * td_div) {w_d10 <- delay(cycling_cells, 10 * td_div) + delay(damaged_cells1, 10 * td_div)}
    p_d10 <- 1; if (t >= 9.5 * td_div) {p_d10 <- 1 + pwth / w_d10}
    w_d11 <- rbase
    if (t >= 11 * td_div) {w_d11 <- delay(cycling_cells, 11 * td_div) + delay(damaged_cells1, 11 * td_div)}
    p_d11 <- 1; if (t >= 10.5 * td_div) {p_d11 <- 1 + pwth / w_d11}
    w_d12 <- rbase
    if (t >= 12 * td_div) {w_d12 <- delay(cycling_cells, 12 * td_div) + delay(damaged_cells1, 12 * td_div)}
    p_d12 <- 1; if (t >= 11.5 * td_div) {p_d12 <- 1 + pwth / w_d12}
    w_d13 <- rbase
    if (t >= 13 * td_div) {w_d13 <- delay(cycling_cells, 13 * td_div) + delay(damaged_cells1, 13 * td_div)}
    p_d13 <- 1; if (t >= 12.5 * td_div) {p_d13 <- 1 + pwth / w_d13}
    w_d14 <- rbase
    if (t >= 14 * td_div) {w_d14 <- delay(cycling_cells, 14 * td_div) + delay(damaged_cells1, 14 * td_div)}
    p_d14 <- 1; if (t >= 13.5 * td_div) {p_d14 <- 1 + pwth / w_d14}
    w_d15 <- rbase
    if (t >= 15 * td_div) {w_d15 <- delay(cycling_cells, 15 * td_div) + delay(damaged_cells1, 15 * td_div)}
    p_d15 <- 1; if (t >= 14.5 * td_div) {p_d15 <- 1 + pwth / w_d15}
    w_d16 <- rbase
    if (t >= 16 * td_div) {w_d16 <- delay(cycling_cells, 16 * td_div) + delay(damaged_cells1, 16 * td_div)}
    p_d16 <- 1; if (t >= 15.5 * td_div) {p_d16 <- 1 + pwth / w_d16}
    w_d17 <- rbase
    if (t >= 17 * td_div) {w_d17 <- delay(cycling_cells, 17 * td_div) + delay(damaged_cells1, 17 * td_div)}
    p_d17 <- 1; if (t >= 16.5 * td_div) {p_d17 <- 1 + pwth / w_d17}
    w_d18 <- rbase
    if (t >= 18 * td_div) {w_d18 <- delay(cycling_cells, 18 * td_div) + delay(damaged_cells1, 18 * td_div)}
    p_d18 <- 1; if (t >= 17.5 * td_div) {p_d18 <- 1 + pwth / w_d18}
    w_d19 <- rbase
    if (t >= 19 * td_div) {w_d19 <- delay(cycling_cells, 19 * td_div) + delay(damaged_cells1, 19 * td_div)}
    p_d19 <- 1; if (t >= 18.5 * td_div) {p_d19 <- 1 + pwth / w_d19}
    w_d20 <- rbase
    if (t >= 20 * td_div) {w_d20 <- delay(cycling_cells, 20 * td_div) + delay(damaged_cells1, 20 * td_div)}
    p_d20 <- 1; if (t >= 19.5 * td_div) {p_d20 <- 1 + pwth / w_d20}
    w_d21 <- rbase
    if (t >= 21 * td_div) {w_d21 <- delay(cycling_cells, 21 * td_div) + delay(damaged_cells1, 21 * td_div)}
    p_d21 <- 1; if (t >= 20.5 * td_div) {p_d21 <- 1 + pwth / w_d21}
    w_d22 <- rbase
    if (t >= 22 * td_div) {w_d22 <- delay(cycling_cells, 22 * td_div) + delay(damaged_cells1, 22 * td_div)}
    p_d22 <- 1; if (t >= 21.5 * td_div) {p_d22 <- 1 + pwth / w_d22}
    w_d23 <- rbase
    if (t >= 23 * td_div) {w_d23 <- delay(cycling_cells, 23 * td_div) + delay(damaged_cells1, 23 * td_div)}
    p_d23 <- 1; if (t >= 22.5 * td_div) {p_d23 <- 1 + pwth / w_d23}
    w_d24 <- rbase
    if (t >= 24 * td_div) {w_d24 <- delay(cycling_cells, 24 * td_div) + delay(damaged_cells1, 24 * td_div)}
    p_d24 <- 1; if (t >= 23.5 * td_div) {p_d24 <- 1 + pwth / w_d24}
    w_d25 <- rbase
    if (t >= 25 * td_div) {w_d25 <- delay(cycling_cells, 25 * td_div) + delay(damaged_cells1, 25 * td_div)}
    p_d25 <- 1; if (t >= 24.5 * td_div) {p_d25 <- 1 + pwth / w_d25}
    w_d26 <- rbase
    if (t >= 26 * td_div) {w_d26 <- delay(cycling_cells, 26 * td_div) + delay(damaged_cells1, 26 * td_div)}
    p_d26 <- 1; if (t >= 25.5 * td_div) {p_d26 <- 1 + pwth / w_d26}
    w_d27 <- rbase
    if (t >= 27 * td_div) {w_d27 <- delay(cycling_cells, 27 * td_div) + delay(damaged_cells1, 27 * td_div)}
    p_d27 <- 1; if (t >= 26.5 * td_div) {p_d27 <- 1 + pwth / w_d27}
    w_d28 <- rbase
    if (t >= 28 * td_div) {w_d28 <- delay(cycling_cells, 28 * td_div) + delay(damaged_cells1, 28 * td_div)}
    p_d28 <- 1; if (t >= 27.5 * td_div) {p_d28 <- 1 + pwth / w_d28}
    w_d29 <- rbase
    if (t >= 29 * td_div) {w_d29 <- delay(cycling_cells, 29 * td_div) + delay(damaged_cells1, 29 * td_div)}
    p_d29 <- 1; if (t >= 28.5 * td_div) {p_d29 <- 1 + pwth / w_d29}
    w_d30 <- rbase
    if (t >= 30 * td_div) {w_d30 <- delay(cycling_cells, 30 * td_div) + delay(damaged_cells1, 30 * td_div)}
    p_d30 <- 1; if (t >= 29.5 * td_div) {p_d30 <- 1 + pwth / w_d30}
    w_d31 <- rbase
    if (t >= 31 * td_div) {w_d31 <- delay(cycling_cells, 31 * td_div) + delay(damaged_cells1, 31 * td_div)}
    p_d31 <- 1; if (t >= 30.5 * td_div) {p_d31 <- 1 + pwth / w_d31}
    w_d32 <- rbase
    if (t >= 32 * td_div) {w_d32 <- delay(cycling_cells, 32 * td_div) + delay(damaged_cells1, 32 * td_div)}
    p_d32 <- 1; if (t >= 31.5 * td_div) {p_d32 <- 1 + pwth / w_d32}
    w_d33 <- rbase
    if (t >= 33 * td_div) {w_d33 <- delay(cycling_cells, 33 * td_div) + delay(damaged_cells1, 33 * td_div)}
    p_d33 <- 1; if (t >= 32.5 * td_div) {p_d33 <- 1 + pwth / w_d33}
    w_d34 <- rbase
    if (t >= 34 * td_div) {w_d34 <- delay(cycling_cells, 34 * td_div) + delay(damaged_cells1, 34 * td_div)}
    p_d34 <- 1; if (t >= 33.5 * td_div) {p_d34 <- 1 + pwth / w_d34}
    w_d35 <- rbase
    if (t >= 35 * td_div) {w_d35 <- delay(cycling_cells, 35 * td_div) + delay(damaged_cells1, 35 * td_div)}
    p_d35 <- 1; if (t >= 34.5 * td_div) {p_d35 <- 1 + pwth / w_d35}
    w_d36 <- rbase
    if (t >= 36 * td_div) {w_d36 <- delay(cycling_cells, 36 * td_div) + delay(damaged_cells1, 36 * td_div)}
    p_d36 <- 1; if (t >= 35.5 * td_div) {p_d36 <- 1 + pwth / w_d36}
    w_d37 <- rbase
    if (t >= 37 * td_div) {w_d37 <- delay(cycling_cells, 37 * td_div) + delay(damaged_cells1, 37 * td_div)}
    p_d37 <- 1; if (t >= 36.5 * td_div) {p_d37 <- 1 + pwth / w_d37}
    w_d38 <- rbase
    if (t >= 38 * td_div) {w_d38 <- delay(cycling_cells, 38 * td_div) + delay(damaged_cells1, 38 * td_div)}
    p_d38 <- 1; if (t >= 37.5 * td_div) {p_d38 <- 1 + pwth / w_d38}
    w_d39 <- rbase
    if (t >= 39 * td_div) {w_d39 <- delay(cycling_cells, 39 * td_div) + delay(damaged_cells1, 39 * td_div)}
    p_d39 <- 1; if (t >= 38.5 * td_div) {p_d39 <- 1 + pwth / w_d39}
    w_d40 <- rbase
    if (t >= 40 * td_div) {w_d40 <- delay(cycling_cells, 40 * td_div) + delay(damaged_cells1, 40 * td_div)}
    p_d40 <- 1; if (t >= 39.5 * td_div) {p_d40 <- 1 + pwth / w_d40}
    prod_div <- 1 *
      p_d1 * p_d2 * p_d3 * p_d4 * p_d5 * p_d6 * p_d7 * p_d8 * p_d9 * p_d10 *
      p_d11 * p_d12 * p_d13 * p_d14 * p_d15 * p_d16 * p_d17 * p_d18 * p_d19 * p_d20 *
      p_d21 * p_d22 * p_d23 * p_d24 * p_d25 * p_d26 * p_d27 * p_d28 * p_d29 * p_d30 *
      p_d31 * p_d32 * p_d33 * p_d34 * p_d35 * p_d36 * p_d37 * p_d38 * p_d39 * p_d40

    # Apoptotic outflow E(C(t - TA)) M(t - TA) (Eq. 24)
    c_apo <- 0
    m_apo <- 0
    if (t >= td_apo) {
      c_apo <- delay(central, td_apo) / vc * 1000
      m_apo <- delay(cycling_cells, td_apo)
    }
    eff_apo <- drugSlope * c_apo

    # Proliferating pool (Eq. 22 with round() for INT()), apoptotic pool
    # (Eq. 24) and the surviving fraction exp(-cumhaz) of Eq. 22
    cycling_cells(0) <- rbase
    d/dt(cycling_cells) <- kin0 * (pw - 1) * exp(-cumhaz) * prod_div - eff * cycling_cells
    d/dt(damaged_cells1) <- eff * cycling_cells - eff_apo * m_apo
    d/dt(cumhaz) <- eff

    tumor_vol ~ add(addSd_tumor_vol)
  })
}
