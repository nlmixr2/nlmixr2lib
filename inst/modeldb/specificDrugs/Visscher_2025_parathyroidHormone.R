Visscher_2025_parathyroidHormone <- function() {
  description <- paste(
    "Single-patient pharmacokinetic-pharmacodynamic model for recombinant",
    "human parathyroid hormone rhPTH(1-84) (Natpar) in chronic",
    "postsurgical hypoparathyroidism, built in Edsim++ v2.0.5.167 and used",
    "to compare once-daily, multiple-daily and continuous subcutaneous",
    "administration. One-compartment PK with an Erlang-5 transit",
    "absorption chain (dose site plus four transit compartments, all",
    "transfers at ka) and a constant additive endogenous baseline",
    "concentration. Elimination is organ-function-scaled (Visscher 2025",
    "Eq. 3): CLi = CL * (fe * RENALFUNC_REL + (1 - fe) * HEPFUNC_REL),",
    "splitting total clearance into a renal fraction fe = 0.3 that tracks",
    "relative renal function and a hepatic remainder that tracks relative",
    "liver function. Two independent, analyte-specific effect arms hang",
    "off the plasma concentration, each a two-state concentration-driven",
    "delay chain (the paper's transfer compartment T0x plus a virtual",
    "central compartment C0x, sharing one equilibration rate constant)",
    "feeding a sigmoid Emax function whose Emax is an ASYMPTOTE rather",
    "than an increment. The readouts are the renal clearances of",
    "phosphate and of calcium expressed as a percentage of creatinine",
    "clearance (fractional excretion; Eq. 2 and Supplementary Eq. 5).",
    "Phosphate is the better-identified arm (goodness-of-fit R^2 = 0.66);",
    "the calcium arm is poorly identified (R^2 = 0.06) because concomitant",
    "alfacalcidol, calcium carbonate and chlortalidone perturb calcium",
    "homeostasis. Deterministic typical-value model: the source reports no",
    "inter-individual variability and no residual-error model, as it was",
    "fitted to the intensive data of a single patient."
  )
  reference <- paste(
    "Visscher M, Schuls-Fouchier M, Berends AMA, Muller Kobold AC,",
    "Punt NC, Touw DJ.",
    "Personalized parathyroid hormone therapy for hypoparathyroidism:",
    "Insights from pharmacokinetic-pharmacodynamic modelling.",
    "Br J Clin Pharmacol. 2025;91(4):1233-1240.",
    "doi:10.1111/bcp.16342."
  )
  vignette <- "Visscher_2025_parathyroidHormone"

  # T01/T02 (the paper's "transfer compartment") and C03/C02 (the paper's
  # "virtual central compartment") are paper-mechanistic effect-delay states,
  # one pair per analyte. Supporting Information Figure S1 caption: "The delay
  # on the effect parameters was modeled by using a transfer compartment (T02
  # and T01 object) connected to a virtual central compartment (C02 and C03
  # object). The virtual compartments do not participate in the mass balance
  # of the model." Suffixed per analyte following the precedent of
  # Heo_2016_amlodipine_valsartan.R (`c("effect_sbp", "effect_dbp")`).
  paper_specific_compartments <- c(
    "transfer_phosphate", "effect_phosphate",
    "transfer_calcium", "effect_calcium"
  )

  # The model is encoded natively in the source's own molar units, so no
  # molecular weight enters ini() or model() and every value below is
  # directly source-traceable. Table S1 states the clinical doses in ug; the
  # ug <-> pmol conversion (MW 9425 g/mol for rhPTH(1-84)) is done in the
  # validation vignette only.
  units <- list(time = "h", dosing = "pmol", concentration = "pmol/L")

  compartmentData <- list(
    depot              = list(analyte = "parathyroid hormone", units = "pmol", specimen = "administration site", verified = TRUE),
    transit1           = list(analyte = "parathyroid hormone", units = "pmol", specimen = "administration site", verified = TRUE),
    transit2           = list(analyte = "parathyroid hormone", units = "pmol", specimen = "administration site", verified = TRUE),
    transit3           = list(analyte = "parathyroid hormone", units = "pmol", specimen = "administration site", verified = TRUE),
    transit4           = list(analyte = "parathyroid hormone", units = "pmol", specimen = "administration site", verified = TRUE),
    central            = list(analyte = "parathyroid hormone", units = "pmol", specimen = "plasma", verified = TRUE),
    transfer_phosphate = list(analyte = "parathyroid hormone", units = "pmol/L", specimen = "not applicable", verified = TRUE),
    effect_phosphate   = list(analyte = "parathyroid hormone", units = "pmol/L", specimen = "not applicable", verified = TRUE),
    transfer_calcium   = list(analyte = "parathyroid hormone", units = "pmol/L", specimen = "not applicable", verified = TRUE),
    effect_calcium     = list(analyte = "parathyroid hormone", units = "pmol/L", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    RENALFUNC_REL = list(
      description        = "Relative renal function (1 = normal)",
      units             = "(dimensionless)",
      type              = "continuous",
      reference_category = NULL,
      notes = paste(
        "Dimensionless organ-function scalar of the Edsim++ / MwPharm ORG",
        "('P') object, consumed by Visscher 2025 Eq. 3 as the multiplier on",
        "the renal fraction fe of total clearance:",
        "CLi = CL * (fe * RENALFUNC_REL + (1 - fe) * HEPFUNC_REL).",
        "The paper states that the equation 'includes the variations in the",
        "GFR by correcting for the relative renal function (RF)' and shows",
        "an increasing GFR over the 13-day study in Figure S2, but it",
        "reports NO numeric RF values and NO reference GFR denominator from",
        "which RF could be recovered, so the value is left at 1 (normal) in",
        "the validation vignette. Table S3 gives urine creatinine but not",
        "the 24-h urine volume that Eq. 1 needs, so the per-interval GFR",
        "cannot be recomputed from the supplement either. A user with a",
        "measured GFR can scale renal elimination by supplying this column",
        "as RENALFUNC_REL = GFR / GFR_normal.",
        "Source name RF; the bare two-letter name is a rejected alias of the",
        "rheumatoid-factor covariate and is deliberately not reused."
      ),
      source_name        = "RF"
    ),
    HEPFUNC_REL = list(
      description        = "Relative liver function (1 = normal)",
      units             = "(dimensionless)",
      type              = "continuous",
      reference_category = NULL,
      notes = paste(
        "Dimensionless organ-function scalar multiplying the hepatic",
        "remainder (1 - fe) of total clearance in Visscher 2025 Eq. 3.",
        "Methods, section 2.3: 'LF is the relative liver function (equal to",
        "1 in this model)', so the source holds it at 1 throughout and no",
        "hepatic-impairment data inform it. Retained in the model rather",
        "than folded into CL so that the paper's own decomposition of",
        "clearance into a renal and a hepatic arm stays visible and",
        "user-perturbable."
      ),
      source_name        = "LF"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1,
    n_studies      = 1,
    age_median     = "43 years",
    sex_female_pct = 0,
    disease_state  = paste(
      "Chronic primary (iatrogenic postsurgical) hypoparathyroidism",
      "following total thyroidectomy for medullary thyroid carcinoma,",
      "complicated by severe hypercalciuria and recurrent nephrolithiasis"
    ),
    dose_range     = paste(
      "rhPTH(1-84) subcutaneous: 100 ug QD (days 1-2), 50 ug BID",
      "(days 3-4), 35 ug TID (days 5-6), 25 ug QID (days 7-8), then",
      "continuous subcutaneous infusion 100 ug/day (days 9-10),",
      "75 ug/day (day 11), 50 ug/day (day 12), 25 ug/day (day 13)"
    ),
    regions        = "the Netherlands (University Medical Center Groningen)",
    co_medication  = paste(
      "Alfacalcidol 1.25 ug QD tapered to 0.25 ug QD then stopped",
      "(days 0-2), calcium carbonate 500 mg QD (days 0-3),",
      "chlortalidone 12.5 mg BID throughout (Table S1)"
    ),
    notes = paste(
      "Single-patient (n = 1) intensive study. Baseline demographics and",
      "the full day-by-day medication scheme are in Supporting Information",
      "Table S1; serum sampling in Table S2 and timed urine collections in",
      "Table S3. Body weight is not reported, so no allometric scaling is",
      "possible or needed. The Abstract describes the patient as 42 years",
      "old and Methods section 2.1 as 43 years old; the Methods value is",
      "used here. Endogenous PTH before rhPTH(1-84) administration was",
      "0.58 pmol/L. The patient maintained a consistent diet throughout",
      "the study period."
    )
  )

  ini({
    # --- Fixed parameters (Table S4, "Fixed parameters" block) ---------------
    # I01.F bioavailability, fixed at the Natpar SmPC literature value.
    # Methods 2.3: "bioavailability was set at a literature value of 0.53".
    lfdepot <- fixed(log(0.53)); label("Subcutaneous bioavailability (fraction)")  # Table S4, I01.F = 0.53 (no SE)

    # C01.Cb endogenous baseline concentration. Methods 2.3: "Before the
    # administration of PTH, an endogenous PTH concentration of 0.58 pmol/L
    # was determined and this was used as the baseline level in the PKPD
    # model." Enters as a constant additive offset on the plasma
    # concentration, not as a decaying initial condition.
    bl_pth <- fixed(0.58); label("Endogenous baseline PTH concentration (pmol/L)")  # Table S4, C01.Cb = 0.58 (no SE)

    # ORG.nfe "Normal fraction excreted unchanged" -- the renal share of total
    # clearance in healthy subjects, against which RENALFUNC_REL then scales
    # (Eq. 3). Methods 2.3: "According to the study by Hruska et al, the renal
    # clearance of PTH in dogs was estimated at an average of 30% of the total
    # clearance, therefore for our base model an NFE of 0.3 was used."
    fe <- fixed(0.3); label("Fraction of clearance excreted unchanged in urine, healthy-subject reference (fraction)")  # Table S4, ORG.nfe = 0.3 (no SE)

    # --- PK fit (Table S4, "PK fit" block) ----------------------------------
    lka <- log(8.39);   label("Transit-chain first-order absorption rate constant (1/h)")  # Table S4, I02.k = 8.39 (SE 2.81)
    lcl <- log(39.93);  label("Population clearance (L/h)")                                # Table S4, O01.CL = 39.93 (SE 2.79)
    lvc <- log(120.71); label("Central volume of distribution (L)")                         # Table S4, C01.V = 120.71 (SE 11.71)

    # --- PD fit, phosphate (Table S4, "PD fit phosphate" block) -------------
    # Well-identified arm: goodness-of-fit R^2 = 0.66 (Discussion / Figure 3B).
    lke0_phosphate  <- log(0.45); label("Effect-delay equilibration rate constant, phosphate arm (1/h)")  # Table S4, T01.k = 0.45 (SE 0.08)
    e0_phosphate    <- 14.18;     label("Phosphate clearance in the absence of drug, as % of creatinine clearance")  # Table S4, E01.E0 = 14.18 (SE 4.72)
    emax_phosphate  <- 50.62;     label("Maximum (asymptotic) phosphate clearance, as % of creatinine clearance")    # Table S4, E01.Emax = 50.62 (SE 9.65)
    lec50_phosphate <- log(6.35); label("PTH concentration at half-maximal phosphate effect (pmol/L)")  # Table S4, E01.EC50 = 6.35 (SE 1.05)
    hill_phosphate  <- 3.82;      label("Hill coefficient, phosphate arm (unitless)")                    # Table S4, E01.j = 3.82 (SE 2.18)

    # --- PD fit, calcium (Table S4, "PD fit calcium" block) -----------------
    # POORLY IDENTIFIED arm. Every value below carries the Table S4 dagger
    # footnote "These parameters were hard to fit since they were very
    # sensitive for the initial value upon fitting", the SEs exceed (often far
    # exceed) the estimates, and the goodness-of-fit R^2 is 0.06 (Discussion /
    # Figure 3A). The Discussion attributes this to concomitant alfacalcidol,
    # calcium carbonate and chlortalidone perturbing calcium homeostasis.
    # Reproduced faithfully; do not treat as equal-confidence with phosphate.
    lke0_calcium  <- log(4.16); label("Effect-delay equilibration rate constant, calcium arm (1/h)")  # Table S4, T02.k = 4.16 (SE 63.19; dagger footnote)
    e0_calcium    <- 2.57;      label("Calcium clearance in the absence of drug, as % of creatinine clearance")  # Table S4, E02.E0 = 2.57 (SE 4.51; dagger footnote)
    emax_calcium  <- 2.05;      label("Maximum (asymptotic) calcium clearance, as % of creatinine clearance")    # Table S4, E02.Emax = 2.05 (SE 3.85; dagger footnote)
    lec50_calcium <- log(3.56); label("PTH concentration at half-maximal calcium effect (pmol/L)")  # Table S4, E02.EC50 = 3.56 (SE 34.88; dagger footnote)
    hill_calcium  <- 3.42;      label("Hill coefficient, calcium arm (unitless)")                    # Table S4, E02.j = 3.42 (SE 100.71; dagger footnote)
  })

  model({
    # 1. Individual parameters. No inter-individual variability: the source is
    #    a single-patient fit and reports no OMEGA block.
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)

    ke0_phosphate  <- exp(lke0_phosphate)
    ec50_phosphate <- exp(lec50_phosphate)
    ke0_calcium    <- exp(lke0_calcium)
    ec50_calcium   <- exp(lec50_calcium)

    # 2. Organ-function-scaled individual clearance (Visscher 2025 Eq. 3):
    #      CLi = CL * (NFE * RF + (1 - NFE) * LF)
    #    The renal arm carries the fraction excreted unchanged fe and scales
    #    with relative renal function; the hepatic remainder scales with
    #    relative liver function. Both covariates are 1 for normal organ
    #    function, which is the value the source used (LF explicitly, RF by
    #    absence of reported values).
    cli <- cl * (fe * RENALFUNC_REL + (1 - fe) * HEPFUNC_REL)
    kel <- cli / vc

    # 3. PK ODEs. Erlang-5 subcutaneous absorption: the dose site plus four
    #    transit compartments, five first-order transfers all at ka
    #    (Methods 2.3: "to match the delay for the absorption of PTH, five
    #    transit compartments were included in the absorption compartment
    #    based on the fitting of the model to the data"). Only one absorption
    #    rate constant (I02.k) is reported, so the chain is homogeneous.
    d/dt(depot)    <- -ka * depot
    d/dt(transit1) <-  ka * depot    - ka * transit1
    d/dt(transit2) <-  ka * transit1 - ka * transit2
    d/dt(transit3) <-  ka * transit2 - ka * transit3
    d/dt(transit4) <-  ka * transit3 - ka * transit4
    d/dt(central)  <-  ka * transit4 - kel * central

    f(depot) <- exp(lfdepot)

    # 4. Plasma concentration = exogenous drug plus the constant endogenous
    #    baseline (Edsim++ C01.Cb). This is the quantity Figure 2A and Figure
    #    4A plot as "PTH concentration in the central compartment (C01)".
    Cc <- central / vc + bl_pth

    # 5. Effect-delay chains. One transfer compartment (T01 / T02) feeding one
    #    virtual central compartment (C03 / C02) per analyte, both transfers
    #    driven by the single reported equilibration rate constant T0x.k. The
    #    virtual compartments hold a CONCENTRATION and are excluded from the
    #    mass balance (Figure S1 caption), so they are driven by Cc rather
    #    than exchanging amounts with `central`. Both states start at the
    #    endogenous baseline so the simulation begins at rest rather than
    #    ramping up from zero. Figure 4B plots the effect against the virtual
    #    compartment C03, confirming the second state is the effect driver.
    transfer_phosphate(0) <- bl_pth
    effect_phosphate(0)   <- bl_pth
    transfer_calcium(0)   <- bl_pth
    effect_calcium(0)     <- bl_pth

    d/dt(transfer_phosphate) <- ke0_phosphate * (Cc - transfer_phosphate)
    d/dt(effect_phosphate)   <- ke0_phosphate * (transfer_phosphate - effect_phosphate)
    d/dt(transfer_calcium)   <- ke0_calcium * (Cc - transfer_calcium)
    d/dt(effect_calcium)     <- ke0_calcium * (transfer_calcium - effect_calcium)

    # 6. Readouts: the renal clearance of phosphate and of calcium expressed
    #    as a percentage of creatinine clearance (Eq. 2, expanded in
    #    Supplementary Eq. 5 as
    #    CL(%) = (Curine_X * Cplasma_creat) / (Cplasma_X * Curine_creat)).
    #    Emax is an ASYMPTOTE, not an increment: the readout runs from e0 at
    #    zero drug to emax at saturating drug. The calcium arm has
    #    emax (2.05%) BELOW e0 (2.57%), so only the asymptote reading
    #    reproduces the Conclusion that "a permanent increase in phosphate
    #    clearance was observed which corresponds to a permanent reduction in
    #    the calcium clearance"; an increment reading would make calcium
    #    clearance rise instead.
    clrel_phosphate <-
      e0_phosphate + (emax_phosphate - e0_phosphate) *
      effect_phosphate^hill_phosphate /
      (ec50_phosphate^hill_phosphate + effect_phosphate^hill_phosphate)
    clrel_calcium <-
      e0_calcium + (emax_calcium - e0_calcium) *
      effect_calcium^hill_calcium /
      (ec50_calcium^hill_calcium + effect_calcium^hill_calcium)

    # No residual-error model: Table S4 reports no SIGMA and the source fit a
    # single patient's intensive data deterministically in Edsim++. The model
    # is intended for typical-value simulation.
  })
}
