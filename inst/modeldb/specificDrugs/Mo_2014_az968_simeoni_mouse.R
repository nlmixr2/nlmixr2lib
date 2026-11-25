Mo_2014_az968_simeoni_mouse <- function() {
  description <- "Preclinical (mouse, HCT116 human colon carcinoma xenograft). Reference Simeoni tumor growth inhibition model of the casein kinase 2 inhibitor AZ968 in its linear-growth form, fitted by Mo 2014 to the same data as the lifespan model Mo_2014_az968_mouse for comparison (Eq. 33, Table 4 right). Proliferating cells grow at the zero-order-in-volume rate lambda1 * x1 / w (linear growth of the whole tumor); AZ968 moves proliferating cells at the second-order rate k2 * C(t) into a chain of three damaged-cell transit compartments with rate constant k1 (mean transit time 4 / k1 in Mo 2014). Two-compartment AZ968 PK with first-order absorption from the intraperitoneal injection site and dose-specific elimination rate constant and central volume (non-linear PK; Table 3), selected by the DOSE_AZ968_MGKG covariate; the PK is the same as in Mo_2014_az968_mouse. Time is in days, doses in mg/kg into depot, plasma concentration in ng/mL and tumor volume in mm^3."
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
    cycling_cells = list(analyte = "proliferating tumor cells x1", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(analyte = "damaged tumor cells x2", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells2 = list(analyte = "damaged tumor cells x3", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells3 = list(analyte = "damaged tumor cells x4", units = "mm^3", specimen = "tumor", verified = TRUE)
  )

  population <- list(
    species = "mouse (female NCr nude, HCT116 xenograft; PK in female CD1 mice)",
    n_subjects = NA_integer_,
    n_studies = 2L,
    sex_female_pct = 100,
    age_range = "5-6 weeks",
    disease_state = "subcutaneous HCT116 human colon carcinoma xenograft, treatment started at 150-200 mm^3",
    dose_range = "AZ968 10, 20 or 30 mg/kg intraperitoneal once weekly for three weeks (tumor study, Figure 8C); single 10, 20 or 30 mg/kg intraperitoneal doses (PK study, 3 mice per dose)",
    notes = "Same AstraZeneca data as Mo_2014_az968_mouse (Mo 2014 Methods, Animal Data). Mean data fitted with MATLAB nlinfit. Initial tumor volume 180 mm^3 (Figure 8 caption)."
  )

  ini({
    # Simeoni TGI parameters: Mo 2014 Table 4 (right, Simeoni TGI Model).
    ltumorLinGrowth <- log(43.8); label("Linear tumor growth rate lambda1 (mm^3/day)") # Table 4: lambda1 = 43.8 (units printed as 1/day)
    ldamageTransit <- log(1.12); label("Damaged-cell transit rate constant k1 (1/day)") # Table 4: k1 = 1.12 1/day
    ldrugSlope <- log(2.0e-3); label("Second-order drug potency constant k2 (mL/ng/day)") # Table 4: k2 = 2.0e-3 mL/ng/day
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
    tumorLinGrowth <- exp(ltumorLinGrowth)
    damageTransit <- exp(ldamageTransit)
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

    # Linearised Simeoni TGI model (Eq. 33)
    tumor_vol <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3
    cycling_cells(0) <- rbase
    d/dt(cycling_cells) <- tumorLinGrowth * cycling_cells / tumor_vol - drugSlope * Cc * cycling_cells
    d/dt(damaged_cells1) <- drugSlope * Cc * cycling_cells - damageTransit * damaged_cells1
    d/dt(damaged_cells2) <- damageTransit * (damaged_cells1 - damaged_cells2)
    d/dt(damaged_cells3) <- damageTransit * (damaged_cells2 - damaged_cells3)

    tumor_vol ~ add(addSd_tumor_vol)
  })
}
