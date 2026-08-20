Olafuyi_2025_propyleneGlycol_adult <- function() {
  description <- paste(
    "One-compartment reduction of the Olafuyi 2025 Simcyp full-PBPK model for the",
    "excipient propylene glycol (PG) in healthy adults: first-order oral absorption",
    "with parallel saturable (Michaelis-Menten, ADH-mediated hepatic) and linear",
    "renal elimination. Reproduces the paper's dose-dependent clearance saturation",
    "and the 580 mg/L plasma level associated with toxicity (PLAT).",
    sep = " "
  )
  reference <- paste(
    "Olafuyi O, Michelet R, Garle M, Allegaert K (2025). Exploring the Impact of",
    "Developmental Clearance Saturation on Propylene Glycol Exposure in Adults and",
    "Term Neonates Using Physiologically Based Pharmacokinetic Modeling.",
    "The Journal of Clinical Pharmacology 65(3):272-284. doi:10.1002/jcph.6150.",
    "Term-neonate counterpart: modellib('Olafuyi_2025_propyleneGlycol_neonate').",
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
        "Scales the central volume only. Olafuyi 2025 Table 1 reports the adult",
        "steady-state volume of distribution per kilogram (Vss = 0.80 L/kg), so",
        "vc = 0.80 * WT. Because the paper's dose-response work is expressed in",
        "mg/kg, this makes the initial concentration after an mg/kg dose",
        "weight-independent (C0 = dose_mg_per_kg / 0.80), which is why the",
        "published saturation thresholds are quoted in mg/kg. The paper reports",
        "hepatic and renal clearance as absolute L/h for the simulated adult",
        "population and publishes no weight-scaling exponent for either, so",
        "neither clearance arm is weight-scaled here (see the vignette",
        "Assumptions and deviations section)."
      ),
      source_name        = "Body weight (Simcyp virtual healthy-volunteer population)"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "propyleneGlycol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "propyleneGlycol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 200L,
    n_studies      = 2L,
    age_range      = "Adults (Simcyp built-in virtual healthy-volunteer population)",
    weight_range   = "Not reported; Simcyp virtual healthy-volunteer body-weight distribution",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported (Simcyp built-in virtual healthy-volunteer population)",
    disease_state  = "Healthy adults; propylene glycol administered as a pharmaceutical excipient",
    dose_range     = paste(
      "Simulated intravenous PG doses of 0.75-7500 mg/kg per event, given 6-, 8-,",
      "12-, and 24-hourly (Methods, 'Determination of Dose and Clearance",
      "Relationship'). Figures 4 and 5 simulate 0.75, 200, 500, and 1000 mg/kg",
      "per event; the proposed safe adult total daily dose is 100-200 mg/kg/day."
    ),
    regions        = "Not applicable (simulation study; UK-based investigators)",
    notes          = paste(
      "Simulation population: 200 virtual subjects per scenario from the Simcyp",
      "(Version 20) built-in healthy-volunteer population. The compound layer was",
      "built from PG physicochemical properties (Table 1) plus alcohol-",
      "dehydrogenase kinetics measured in pooled human liver cytosol in this study",
      "(Table 2). Model qualification used digitised concentration-time profiles",
      "from two published adult clinical studies (Table S3, supplement not on",
      "disk) and the clearance comparison in Table 3. Not a fitted population-PK",
      "analysis: no observed individual-level data were fitted, so the variability",
      "terms below are Simcyp virtual-population output %CVs, not estimated omegas."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Compound layer. Every value is a fixed model INPUT or a fixed quantity
    # derived from a printed model OUTPUT -- nothing here was estimated by
    # nlmixr2, so all are wrapped in fixed().
    # ---------------------------------------------------------------------

    # Absorption (Table 1, "Optimized Adult Model" column). Both are Simcyp
    # mechanistic predictions from PG physicochemical properties (footnote a).
    lka     <- fixed(log(7.56)); label("First-order absorption rate constant (1/h)")  # Table 1: ka = 7.56 per h
    lfdepot <- fixed(log(0.99)); label("Fraction of an oral dose absorbed (unitless)")  # Table 1: fa = 0.99

    # Distribution (Table 1). Vss is reported per kilogram; footnote c states it
    # is the calculated average volume of distribution reported in the
    # literature (reference 6), i.e. a model INPUT, and the adult Kp scalar of
    # 1.0 means the predicted tissue partition values were used unscaled.
    lvc <- fixed(log(0.80)); label("Steady-state volume of distribution per body weight (L/kg)")  # Table 1: Vss = 0.80 L/kg

    # Saturable ADH-mediated hepatic metabolism.
    #
    # km: Table 2 reports Km = 25.1 mM (95% CI 15.7-39.2) measured in pooled
    # human liver cytosol. Converted to the plasma mass scale with the Table 1
    # molecular weight: 25.1 mmol/L * 76.06 mg/mmol = 1909.1 mg/L.
    lkm <- fixed(log(1909.1)); label("Michaelis constant for ADH-mediated metabolism (mg/L)")  # Table 2: Km = 25.1 mM x MW 76.06 g/mol (Table 1)
    #
    # vmax: the in vitro Vmax of 1.57 nmole/min/mg cytosolic protein (Table 2)
    # cannot be scaled to a whole-body rate from on-disk sources, because the
    # liver weight and the cytosolic-protein-per-gram-liver scalar behind
    # Simcyp's in-vitro-to-in-vivo extrapolation are platform system parameters
    # and are not printed anywhere in the paper. It is instead back-solved from
    # the paper's own printed hepatic clearance: in the linear (sub-saturating)
    # concentration range the Michaelis-Menten term reduces to
    # (vmax / km) * Cc, so vmax = CL_hepatic * km = 4.6 L/h * 1909.1 mg/L.
    # See the vignette source-trace and Assumptions sections.
    lvmax <- fixed(log(8781.9)); label("Maximum rate of ADH-mediated hepatic metabolism (mg/h)")  # derived: Table 3 predicted adult hepatic CL 4.6 L/h x km 1909.1 mg/L

    # Linear renal excretion. Table 3 predicted adult renal CL = 3.7 L/h; this
    # is the value consistent with the 55% hepatic / 45% renal split the paper
    # reports for adults in Figure 2 and the Results text (4.6 / (4.6 + 3.7) =
    # 55%). Table 1's CLR input of 3.3 L/h is the pre-simulation estimate that
    # this predicted geometric mean derives from (footnote d: 45% of a total CL
    # of 0.1 L/kg reported by Yu et al 1985).
    lcl_renal <- fixed(log(3.7)); label("Renal clearance (L/h)")  # Table 3: predicted adult renal CL = 3.7 L/h

    # ---------------------------------------------------------------------
    # Variability. Table 3 reports a %CV alongside each predicted geometric-mean
    # clearance across the 200-subject virtual population. These are Simcyp
    # virtual-population OUTPUT dispersions (driven by the platform's
    # physiological variability), not estimated random effects, so they are
    # fixed. No %CV is reported for the volume of distribution, so vc carries
    # no eta. omega^2 = log(1 + CV^2).
    # ---------------------------------------------------------------------
    etalvmax     ~ fixed(0.184417)  # Table 3: predicted adult hepatic CL %CV = 45; log(1 + 0.45^2)
    etalcl_renal ~ fixed(0.051540)  # Table 3: predicted adult renal   CL %CV = 23; log(1 + 0.23^2)

    # Residual error is not reported: the source is a simulation study and no
    # residual-error model was fitted to observed concentrations.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported in source)")  # not reported: the source performs simulation only
  })

  model({
    # 1. Individual parameters. vc is L (Vss per kg x body weight); vmax is
    #    mg/h; km and Cc are mg/L; cl_renal is L/h.
    vc       <- exp(lvc) * WT
    ka       <- exp(lka)
    vmax     <- exp(lvmax + etalvmax)
    km       <- exp(lkm)
    cl_renal <- exp(lcl_renal + etalcl_renal)

    # 2. Plasma concentration, used both as the Michaelis-Menten driver and as
    #    the observation.
    Cc <- central / vc

    # 3. One-compartment disposition with parallel saturable hepatic and linear
    #    renal elimination. At Cc << km the metabolic term collapses to
    #    (vmax / km) * Cc = 4.6 L/h * Cc, recovering the paper's predicted
    #    adult hepatic clearance; as Cc approaches km it saturates, which is the
    #    mechanism the paper set out to quantify.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - vmax * Cc / (km + Cc) - cl_renal * Cc

    # 4. Fraction absorbed on the oral route. Intravenous doses go directly to
    #    central and bypass this term, which is how every simulation reported in
    #    the paper (Figures 3-5) was run.
    f(depot) <- exp(lfdepot)

    # 5. Observation
    Cc ~ prop(propSd)
  })
}
