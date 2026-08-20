Olafuyi_2025_propyleneGlycol_neonate <- function() {
  description <- paste(
    "One-compartment reduction of the Olafuyi 2025 Simcyp full-PBPK model for the",
    "excipient propylene glycol (PG) in term neonates: first-order oral absorption",
    "with parallel saturable (Michaelis-Menten, ADH-mediated hepatic) and linear",
    "renal elimination, where the maximum metabolic rate carries a Hill-type",
    "alcohol-dehydrogenase ontogeny function of age (18% of adult activity at",
    "birth). Reproduces the paper's developmental clearance saturation and the",
    "580 mg/L plasma level associated with toxicity (PLAT).",
    sep = " "
  )
  reference <- paste(
    "Olafuyi O, Michelet R, Garle M, Allegaert K (2025). Exploring the Impact of",
    "Developmental Clearance Saturation on Propylene Glycol Exposure in Adults and",
    "Term Neonates Using Physiologically Based Pharmacokinetic Modeling.",
    "The Journal of Clinical Pharmacology 65(3):272-284. doi:10.1002/jcph.6150.",
    "Adult counterpart: modellib('Olafuyi_2025_propyleneGlycol_adult').",
    sep = " "
  )
  vignette <- "Olafuyi_2025_propyleneGlycol"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Scales the central volume only. Olafuyi 2025 Table 1 reports the",
        "optimized pediatric steady-state volume of distribution per kilogram",
        "(Vss = 0.40 L/kg, obtained by applying an empirically optimized Kp",
        "scalar of 0.44 to the predicted tissue partition values), so",
        "vc = 0.40 * WT. Because the paper's dose-response work is expressed in",
        "mg/kg, this makes the initial concentration after an mg/kg dose",
        "weight-independent (C0 = dose_mg_per_kg / 0.40). The paper reports",
        "hepatic and renal clearance as absolute L/h for the simulated neonatal",
        "population and publishes no weight-scaling exponent for either, so",
        "neither clearance arm is weight-scaled here (see the vignette",
        "Assumptions and deviations section)."
      ),
      source_name        = "Body weight (Simcyp virtual pediatric population, neonatal age group)"
    ),
    AGE = list(
      description        = "Postnatal age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives the alcohol-dehydrogenase ontogeny function that scales vmax",
        "(Olafuyi 2025 Equation 1, a Hill function of age in years with",
        "Fbirth = 0.18, Age50 = 0.9 years and n = 1.4). At AGE = 0 the function",
        "returns exactly 0.18, reproducing the F_ADH activity,neonate value of",
        "18% in Table 1. The paper states its dose-exposure work covers term",
        "neonates only, because clinically observed PG concentration-time data",
        "are not available for older pediatric age groups; the distribution and",
        "renal parameters above are term-neonate values, so the model should be",
        "used at neonatal ages even though Equation 1 itself is continuous in",
        "age."
      ),
      source_name        = "Age (years; the 'Age' term in Equation 1)"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "propyleneGlycol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "propyleneGlycol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 200L,
    n_studies      = 1L,
    age_range      = "Term neonates (Simcyp built-in virtual pediatric population, neonatal age group)",
    weight_range   = "Not reported; Simcyp virtual neonatal body-weight distribution",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported (Simcyp built-in virtual pediatric population)",
    disease_state  = "Term neonates; propylene glycol administered as a pharmaceutical excipient (e.g. in intravenous paracetamol or phenobarbital formulations)",
    dose_range     = paste(
      "Simulated intravenous PG doses of 0.75-7500 mg/kg per event, given 6-, 8-,",
      "12-, and 24-hourly (Methods, 'Determination of Dose and Clearance",
      "Relationship'). Figures 4 and 5 simulate 0.75, 50, 100, and 200 mg/kg per",
      "event; the proposed safe neonatal total daily dose is 25-50 mg/kg/day."
    ),
    regions        = "Validation data from Belgium (University Hospitals Leuven; internal study number B-32220084836)",
    notes          = paste(
      "Simulation population: 200 virtual subjects per scenario from the Simcyp",
      "(Version 20) built-in pediatric population, neonatal age group. The",
      "compound layer is shared with the adult model; the neonatal model differs",
      "in its Kp scalar (0.44, empirically optimized), its resulting Vss",
      "(0.40 L/kg), its renal clearance, and the ADH ontogeny scaling of vmax.",
      "The fractional ADH activity at birth was optimized from 13% to 18% to",
      "recover the fraction metabolized reported by De Cock et al; see",
      "modellib('DeCock_2012_propyleneGlycol') for the empirical population-PK",
      "model fitted to the Leuven neonatal cohort that provides this model's",
      "observed comparator data. Not a fitted population-PK analysis: the",
      "variability terms below are Simcyp virtual-population output %CVs, not",
      "estimated omegas."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Compound layer. Every value is a fixed model INPUT or a fixed quantity
    # derived from a printed model OUTPUT -- nothing here was estimated by
    # nlmixr2, so all are wrapped in fixed().
    # ---------------------------------------------------------------------

    # Absorption (Table 1). Simcyp mechanistic predictions from PG
    # physicochemical properties (footnote a); shared with the adult model.
    lka     <- fixed(log(7.56)); label("First-order absorption rate constant (1/h)")  # Table 1: ka = 7.56 per h
    lfdepot <- fixed(log(0.99)); label("Fraction of an oral dose absorbed (unitless)")  # Table 1: fa = 0.99

    # Distribution (Table 1, "Optimized Pediatric Model" column).
    lvc <- fixed(log(0.40)); label("Steady-state volume of distribution per body weight (L/kg)")  # Table 1: pediatric Vss = 0.40 L/kg (Kp scalar 0.44)

    # Michaelis constant, shared with the adult model: Table 2 reports
    # Km = 25.1 mM (95% CI 15.7-39.2) in pooled human liver cytosol, converted
    # to the plasma mass scale with the Table 1 molecular weight
    # (25.1 mmol/L * 76.06 mg/mmol = 1909.1 mg/L). The paper applies no
    # developmental change to Km -- only ADH activity (vmax) is age-dependent.
    lkm <- fixed(log(1909.1)); label("Michaelis constant for ADH-mediated metabolism (mg/L)")  # Table 2: Km = 25.1 mM x MW 76.06 g/mol (Table 1)

    # Maximum metabolic rate at FULL (adult-equivalent) ADH activity, scaled to
    # the neonatal liver. As in the adult model the in vitro Vmax of
    # 1.57 nmole/min/mg cytosolic protein (Table 2) cannot be scaled to a
    # whole-body rate from on-disk sources, so it is back-solved from the
    # paper's own printed neonatal hepatic clearance and then divided by the
    # ontogeny factor at birth so that the ontogeny function below reproduces
    # it: vmax = CL_hepatic * km / Fbirth = 0.063 L/h * 1909.1 mg/L / 0.18.
    # After multiplication by the Equation 1 ontogeny factor at AGE = 0 this
    # returns 120.3 mg/h, i.e. a hepatic clearance of 0.063 L/h.
    lvmax <- fixed(log(668.19)); label("Maximum rate of ADH-mediated hepatic metabolism at full adult ADH activity (mg/h)")  # derived: Table 3 predicted neonatal hepatic CL 0.063 L/h x km 1909.1 mg/L / Fbirth 0.18

    # Alcohol-dehydrogenase ontogeny, Equation (1):
    #   Fractional change with age = Fbirth + (Fmax - Fbirth) * Age^n /
    #                                (Age50^n + Age^n)
    # with Fmax = 1, Fbirth = 0.18, Age50 = 0.9 years, n = 1.4. Fmax is 1 by
    # definition (maximal fraction of expression), so it is written as the
    # literal 1 in model() rather than carried as a parameter.
    beta_vmax <- fixed(0.18); label("Fraction of adult ADH activity at birth (unitless)")  # Equation 1: Fbirth = 0.18 (Table 1: F_ADH activity,neonate = 18%)
    t50_vmax  <- fixed(0.9);  label("Age at half-maximal ADH activity (years)")            # Equation 1: Age50 = 0.9
    lhill     <- fixed(log(1.4)); label("Hill exponent of the ADH ontogeny function (unitless)")  # Equation 1: n = 1.4

    # Linear renal excretion. Table 3 predicted neonatal renal CL = 0.011 L/h.
    # The paper derives the underlying neonatal glomerular filtration rate from
    # Equation (2), GFR (mL/min) = -19.8 + 89.04*BSA - 7.16*BSA^2, fitted to the
    # renal clearance of PG reported by De Cock et al; that equation
    # parameterises the Simcyp POPULATION renal function rather than PG renal
    # clearance itself (PG is extensively reabsorbed, so its renal clearance is
    # well below GFR), and it is numerically unstable at neonatal body surface
    # areas, so the printed clearance is used directly. See the vignette
    # Assumptions and deviations section.
    lcl_renal <- fixed(log(0.011)); label("Renal clearance (L/h)")  # Table 3: predicted neonatal renal CL = 0.011 L/h

    # ---------------------------------------------------------------------
    # Variability. Table 3 reports a %CV alongside each predicted geometric-mean
    # clearance across the 200-subject virtual population. These are Simcyp
    # virtual-population OUTPUT dispersions, not estimated random effects, so
    # they are fixed. No %CV is reported for the volume of distribution, so vc
    # carries no eta. omega^2 = log(1 + CV^2).
    # ---------------------------------------------------------------------
    etalvmax     ~ fixed(0.169646)  # Table 3: predicted neonatal hepatic CL %CV = 43;  log(1 + 0.43^2)
    etalcl_renal ~ fixed(0.989541)  # Table 3: predicted neonatal renal   CL %CV = 130; log(1 + 1.30^2)

    # Residual error is not reported: the source is a simulation study and no
    # residual-error model was fitted to observed concentrations.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported in source)")  # not reported: the source performs simulation only
  })

  model({
    # 1. ADH ontogeny (Equation 1). Fmax = 1 by definition, so
    #    (Fmax - Fbirth) is written as (1 - beta_vmax). At AGE = 0 this returns
    #    beta_vmax = 0.18 exactly, matching Table 1.
    hill        <- exp(lhill)
    matur_vmax  <- beta_vmax + (1 - beta_vmax) * AGE^hill / (t50_vmax^hill + AGE^hill)

    # 2. Individual parameters. vc is L (Vss per kg x body weight); vmax is
    #    mg/h; km and Cc are mg/L; cl_renal is L/h.
    vc       <- exp(lvc) * WT
    ka       <- exp(lka)
    vmax     <- exp(lvmax + etalvmax) * matur_vmax
    km       <- exp(lkm)
    cl_renal <- exp(lcl_renal + etalcl_renal)

    # 3. Plasma concentration, used both as the Michaelis-Menten driver and as
    #    the observation.
    Cc <- central / vc

    # 4. One-compartment disposition with parallel saturable hepatic and linear
    #    renal elimination. At AGE = 0 and Cc << km the metabolic term collapses
    #    to (vmax * 0.18 / km) * Cc = 0.063 L/h * Cc, recovering the paper's
    #    predicted term-neonate hepatic clearance.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - vmax * Cc / (km + Cc) - cl_renal * Cc

    # 5. Fraction absorbed on the oral route. Intravenous doses go directly to
    #    central and bypass this term, which is how every simulation reported in
    #    the paper (Figures 3-5) was run.
    f(depot) <- exp(lfdepot)

    # 6. Observation
    Cc ~ prop(propSd)
  })
}
