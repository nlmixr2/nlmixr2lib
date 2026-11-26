Oggianu_2020_atomoxetine <- function() {
  description <- paste(
    "One-compartment oral pharmacokinetic reduction of the Simcyp",
    "minimal-PBPK model for atomoxetine in healthy adult CYP2D6 extensive",
    "and poor metabolizers (Oggianu 2020), built as the perpetrator in a",
    "trazodone-atomoxetine interaction prediction. The source model was",
    "run in the Simcyp simulator (version 14) and its equations are not",
    "published, but its compound layer is fully reported in Table 2: oral",
    "clearance CL/F by CYP2D6 phenotype with its between-subject CV, an",
    "absorption rate constant, fraction absorbed 1, blood-to-plasma ratio,",
    "and a steady-state volume of 0.71 L/kg with no single adjusting",
    "compartment, so the distribution is mono-exponential. Systemic",
    "clearance and hepatic first-pass availability follow from CL/F, the",
    "blood-to-plasma ratio and the 90 L/h hepatic blood flow through the",
    "well-stirred liver model; the systemic volume is Vss times body",
    "weight less the 1648 g liver. No parameter is fitted here. The",
    "reduction reproduces the paper's own predicted mean steady-state",
    "Cmax for 20, 40 and 60 mg twice daily within about 6% in all five",
    "arms, and the mean AUC within about 8% in four of them; the fifth",
    "(20 mg extensive metabolizers) lies below a published mean that one",
    "outlying virtual trial pulls upward (see the vignette). The CYP3A4 inhibition",
    "constant used for the interaction is documented but not encoded,",
    "because the trazodone victim model is a Simcyp full-PBPK model that",
    "cannot be reduced from the published inputs.",
    sep = " "
  )
  reference <- paste(
    "Oggianu L, Ke AB, Chetty M, Picollo R, Petrucci V, Calisti F,",
    "Garofolo F, Tongiani S. Estimation of an Appropriate Dose of",
    "Trazodone for Pediatric Insomnia and the Potential for a",
    "Trazodone-Atomoxetine Interaction. CPT Pharmacometrics Syst",
    "Pharmacol. 2020;9(2):77-86. doi:10.1002/psp4.12480.",
    sep = " "
  )
  vignette <- "Oggianu_2020_atomoxetine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(
      analyte = "atomoxetine",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "atomoxetine",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Scales the systemic volume through the per-kilogram steady-state",
        "volume of Oggianu 2020 Table 2 (Vss = 0.71 L/kg); clearance is",
        "not weight-scaled because the paper enters CL/F as an absolute",
        "value in L/h. The median body weight of the paper's adult Simcyp",
        "healthy-volunteer virtual population is 73 kg (Table 4), the",
        "value used for the validation simulations."
      ),
      source_name = "BW"
    ),
    CYP2D6_PM = list(
      description = "1 = CYP2D6 poor metabolizer, 0 = CYP2D6 extensive metabolizer",
      units = "(binary)",
      type = "binary",
      reference_category = "CYP2D6 extensive metabolizer (CYP2D6_PM = 0)",
      notes = paste(
        "Oggianu 2020 Table 2 gives a separate CL/F and CV for CYP2D6",
        "extensive metabolizers (EM: 26.4 L/h, CV 55.7%) and poor",
        "metabolizers (PM: 2.55 L/h, CV 18%). Only these two phenotypes",
        "were simulated; intermediate and ultrarapid metabolizers are not",
        "represented, so CYP2D6_PM = 0 means extensive metabolizer here."
      ),
      source_name = "CYP2D6 EM / PM"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    age_range = "20-50 years (Simcyp healthy-volunteer virtual trials)",
    weight_median = "73 kg (adult virtual population, Table 4)",
    sex_female_pct = 50,
    disease_state = "Healthy adult volunteers, CYP2D6 extensive and poor metabolizers",
    dose_range = "20, 40 and 60 mg twice daily (verification); 60 mg twice daily for 12 days (interaction)",
    regions = "Simcyp Caucasian healthy-volunteer population",
    notes = paste(
      "This is a PBPK simulation analysis, not a population-PK fit.",
      "The atomoxetine verification trials were 10 trials of 7 subjects",
      "(20 mg twice daily, all male) and 10 trials of 8 subjects (40 mg",
      "twice daily, 50% female), compared against Sauer 2003 and the FDA",
      "Strattera clinical pharmacology review; the interaction trial was",
      "10 trials of 10 subjects aged 20-50 years, 50% female",
      "(Supplementary Information S3). CL/F and ka were taken from the",
      "population-PK analysis in the FDA Strattera review; the volume",
      "was optimized by the authors to the 20 mg twice-daily CYP2D6 PM",
      "profile of Sauer 2003."
    )
  )

  ini({
    # Every value is fixed: nothing was estimated for this reduction.
    # Values are Oggianu 2020 Table 2 inputs, or the physiological
    # scalars printed in the Methods 'Estimation of CYP3A4-mediated
    # metabolism for trazodone' (same Simcyp healthy-volunteer population).

    lka <- fixed(log(0.926))
    label("First-order absorption rate constant ka (1/h)")
    # Table 2 'Ka (hour-1)' = 0.926, 'Estimated by Pop-PK analysis' (ref 28,
    # FDA Strattera clinical pharmacology review). Table 2 'Fa' = 1.

    lcl <- fixed(log(26.4))
    label("Population mean apparent oral clearance CL/F in CYP2D6 extensive metabolizers (L/h)")
    # Table 2 'CL/F (L/hour)' CYP2D6 EM: 26.4 (CV%: 55.7). Simcyp treats
    # this input as the arithmetic mean of a log-normal distribution; the
    # median is recovered in model() by subtracting half the log-scale
    # variance.

    e_cyp2d6_pm_cl <- fixed(log(2.55 / 26.4))
    label("Log ratio of CL/F in CYP2D6 poor versus extensive metabolizers (unitless)")
    # Table 2 'CL/F (L/hour)' CYP2D6 PM: 2.55 (CV%: 18); log(2.55/26.4) = -2.337271.

    sdlcl_em <- fixed(0.519824)
    label("Log-scale SD of CL/F between subjects, CYP2D6 extensive metabolizers (unitless)")
    # Table 2 CV 55.7%: sqrt(log(1 + 0.557^2)) = 0.519824.

    sdlcl_pm <- fixed(0.178567)
    label("Log-scale SD of CL/F between subjects, CYP2D6 poor metabolizers (unitless)")
    # Table 2 CV 18%: sqrt(log(1 + 0.18^2)) = 0.178567.

    lvss <- fixed(log(0.71))
    label("Steady-state volume of distribution per kg body weight Vss (L/kg)")
    # Table 2 'Vss (L/kg)' = 0.71, 'Optimized; observed Vss is 0.85 L/kg'.
    # Table 2 'Model' = Minimal-PBPK with no SAC entries.

    lvh <- fixed(log(1.648))
    label("Liver volume subtracted from the systemic volume (L)")
    # Methods (trazodone CLint scaling): average liver weight 1,648 g,
    # taken as 1.648 L.

    lqh <- fixed(log(90))
    label("Hepatic blood flow QH (L/h)")
    # Methods, text below Equation 2: 'QH is the blood flow in the
    # hepatic vein (90 L/hour)'.

    bp <- fixed(0.605)
    label("Blood-to-plasma concentration ratio B:P (unitless)")
    # Table 2 'B/P' = 0.605, 'Predicted by Simcyp'.

    etalcl ~ fixed(1) # standard normal; scaled by sdlcl_em or sdlcl_pm in model()

    propSd <- fixed(0)
    label("Proportional residual error SD (fraction; zero, the source is a simulation model)")
  })

  model({
    ka <- exp(lka)
    vss <- exp(lvss)
    vh <- exp(lvh)
    qh <- exp(lqh)

    # Apparent oral clearance. Simcyp samples CL/F log-normally with the
    # tabulated value as the arithmetic mean, so the median is
    # mean * exp(-sd^2 / 2) and etalcl ~ N(0, 1) carries the
    # phenotype-specific spread.
    sdcl <- sdlcl_em * (1 - CYP2D6_PM) + sdlcl_pm * CYP2D6_PM
    clpo <- exp(lcl + e_cyp2d6_pm_cl * CYP2D6_PM - sdcl^2 / 2 + sdcl * etalcl)

    # Well-stirred liver with fa = fG = 1 and no renal clearance: the
    # oral blood clearance equals fub * CLint, so the plasma fraction
    # unbound cancels and only B:P and QH are needed.
    clint_b <- clpo / bp
    fh <- qh / (qh + clint_b)
    cl <- fh * clpo

    # Systemic volume of the minimal PBPK model: Vss * body weight less
    # the separately modelled liver.
    vc <- vss * WT - vh
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- fh

    # mg / L = ug/mL; x 1000 for ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
