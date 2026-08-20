Adeojo_2024_levonorgestrel <- function() {
  description <- paste(
    "MBMA. Two-compartment population PK model for levonorgestrel fit to pooled",
    "MEAN plasma concentration-time profiles digitised from four published trials",
    "(Adeojo 2024). The four trials span single 0.25 mg oral, single 0.25 mg",
    "intravenous, single 0.75 mg oral, two-dose 0.75 mg oral, and 150 mg subdermal",
    "implant administration in women, with two of the trials contributing arms",
    "co-administered with efavirenz 600 mg orally once daily. Absorption is",
    "first-order and applies to oral dosing only (ka fixed at 4.15 /h); the",
    "intravenous and subdermal-implant routes enter the central compartment",
    "directly, so the estimated bioavailability F is an oral bioavailability.",
    "Concomitant efavirenz enters as a single binary covariate with two effects:",
    "clearance rises from 5.86 to 10.1 L/h (+72.4%) through CYP3A4 induction and",
    "oral bioavailability falls from 0.837 to 0.533 (-36.3%) through gut-wall",
    "CYP3A4 induction. The model carries BETWEEN-STUDY variability on oral",
    "bioavailability only (21.2% CV) and no between-subject variability, because",
    "it was fit to study-level mean profiles rather than individual",
    "concentrations; it therefore simulates study-level mean concentration-time",
    "curves, NOT individual patient concentrations. This is the reference",
    "mixed-effects model of the source paper, used there to anchor a Simcyp",
    "physiologically-based model by retrograde determination; the Simcyp PBPK",
    "model itself is not reproducible from the published sources and is not",
    "packaged here (see the validation vignette Assumptions and deviations).",
    sep = " "
  )
  reference <- paste(
    "Adeojo LW, Patel RC, Sambol NC.",
    "A Physiologically-Based Pharmacokinetic Simulation to Evaluate Approaches",
    "to Mitigate Efavirenz-Induced Decrease in Levonorgestrel Exposure with a",
    "Contraceptive Implant.",
    "Pharmaceutics. 2024;16(8):1050.",
    "doi:10.3390/pharmaceutics16081050.",
    "Parameter estimates from Supplementary Table S2; study designs from",
    "Supplementary Table S1.",
    sep = " "
  )
  vignette <- "Adeojo_2024_levonorgestrel"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "levonorgestrel", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "levonorgestrel", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "levonorgestrel", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CONMED_EFV = list(
      description        = "Concomitant efavirenz coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no efavirenz; no antiretroviral therapy)",
      notes              = paste(
        "Time-fixed per study arm. In every efavirenz-exposed arm of the four pooled",
        "trials the efavirenz dose was 600 mg orally once daily (Carten 2012 oral",
        "levonorgestrel arm and Scarsi 2016 implant arm; Adeojo 2024 Supplementary",
        "Table S1), so the indicator encodes that regimen specifically and does not",
        "carry efavirenz dose or exposure information. Applied as TWO independent",
        "linear-deviation effects: on systemic clearance,",
        "cl *= (1 + 0.7235 * CONMED_EFV), raising CL from 5.86 to 10.1 L/h",
        "(Adeojo 2024 Supplementary Table S2 rows 'CL' and 'CLwith efavirenz'); and on",
        "oral bioavailability, fdepot *= (1 + (-0.3632) * CONMED_EFV), lowering F from",
        "0.837 to 0.533 (Supplementary Table S2 rows 'Foral' and",
        "'Foral, with efavirenz'). Both effects are attributed by the paper to CYP3A4",
        "induction, the larger relative effect on the oral route reflecting the",
        "additional contribution of gut-wall CYP3A4 (Adeojo 2024 Discussion: CL/F is",
        "predicted to increase ~2.7-fold with efavirenz vs ~1.7-fold for systemic CL).",
        "Because the bioavailability effect is carried on f(depot), it applies only to",
        "oral doses; intravenous and subdermal-implant input reach the central",
        "compartment directly and are affected by the clearance term alone.",
        sep = " "
      ),
      source_name        = "efavirenz"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 80L,
    n_studies      = 4L,
    sex_female_pct = 100,
    age_range      = "Not reported in the pooled analysis; the constituent trials enrolled adult women of reproductive potential.",
    disease_state  = paste(
      "Pooled across four trials in adult women (Adeojo 2024 Supplementary",
      "Table S1): Back 1987 (n = 5, healthy women, oral and intravenous",
      "levonorgestrel crossover, co-administered with ethinylestradiol),",
      "Kook 2002 (n = 16, single oral dose), Carten 2012 (n = 21, two oral",
      "doses 12 h apart, crossover with and without efavirenz), and",
      "Scarsi 2016 (n = 18 without efavirenz and n = 20 with efavirenz,",
      "150 mg subdermal implant, women living with HIV on efavirenz-based",
      "antiretroviral therapy).",
      sep = " "
    ),
    dose_range     = paste(
      "Levonorgestrel 0.25 mg single oral, 0.25 mg single intravenous,",
      "0.75 mg single oral, 0.75 mg x 2 oral (12 h apart, PK sampling from",
      "12 h), and 150 mg subdermal implant (two rods, sampling over 48",
      "weeks). Efavirenz 600 mg orally once daily in the interaction arms.",
      sep = " "
    ),
    notes          = paste(
      "This is a meta-analytic fit to AGGREGATE data. The observations are the",
      "MEAN concentration among subjects at each nominal time point in each",
      "trial, digitised from the published figures with GraphClick 3.0;",
      "individual concentration-time data were not available to the authors",
      "(Adeojo 2024 Section 2.2). Two levels of random effects were recognised:",
      "between-study variability of F and within-study variability of the mean",
      "plasma concentrations. The limited number of studies did not permit",
      "estimation of between-study variability on any other parameter. Fit in",
      "NONMEM 7.4.2. Subject count 80 is the sum of distinct participants",
      "across the four trials (5 + 16 + 21 + 18 + 20); the Back 1987 oral/IV",
      "and Carten 2012 with/without-efavirenz arms are crossovers within the",
      "same participants.",
      sep = " "
    )
  )

  ini({
    # ---- Structural fixed effects (Adeojo 2024 Supplementary Table S2, 'Point Estimate' column) ----
    lka  <- fixed(log(4.15)); label("First-order absorption rate constant for oral dosing (1/h)")  # Supp Table S2 'ka' = 4.15 /h; footnote *: fixed value, not estimated, due to limited data in the absorption phase; chosen to give tmax ~1-1.5 h
    lcl  <- log(5.86);        label("Systemic clearance in the absence of efavirenz (L/h)")        # Supp Table S2 'CL' = 5.86 L/h (95% CI 4.97, 6.75)
    lvc  <- log(49.5);        label("Central volume of distribution (L)")                          # Supp Table S2 'Vcentral' = 49.5 L (95% CI 36.3, 62.7)
    lq   <- log(10.2);        label("Inter-compartmental clearance (L/h)")                         # Supp Table S2 'Q' = 10.2 L/h (95% CI 6.9, 13.5); also main-text Table 1 'Q (L/h) 10.2', sourced 'MEM'
    lvp  <- log(105);         label("Peripheral volume of distribution (L)")                       # Supp Table S2 'Vperipheral' = 105 L (95% CI 85.6, 124.4)
    lfdepot <- log(0.837);    label("Oral bioavailability in the absence of efavirenz (unitless)") # Supp Table S2 'Foral' = 0.837 (95% CI 0.736, 0.938); identifiable from the Back 1987 oral/IV crossover arm

    # ---- Concomitant-efavirenz covariate effects (Adeojo 2024 Supplementary Table S2) ----
    # The source reports the efavirenz condition as a second absolute value for CL and for
    # F rather than as a coefficient, so each effect is recorded here in the registry's
    # linear-deviation form, parameter *= (1 + coefficient * CONMED_EFV), with the
    # coefficient computed from the two published point estimates. Both reproduce the
    # published efavirenz value exactly:
    #   cl:     5.86 * (1 + 0.7235)  = 10.1  L/h
    #   fdepot: 0.837 * (1 - 0.3632) = 0.533
    e_conmed_efv_cl      <-  0.7235; label("Factor change in CL with concomitant efavirenz (fraction; +72.4% with CONMED_EFV = 1)")               # Derived from Supp Table S2: (10.1 - 5.86) / 5.86 = 0.72355; 'CLwith efavirenz' = 10.1 L/h (95% CI 9.09, 11.1)
    e_conmed_efv_fdepot  <- -0.3632; label("Factor change in oral bioavailability with concomitant efavirenz (fraction; -36.3% with CONMED_EFV = 1)") # Derived from Supp Table S2: (0.533 - 0.837) / 0.837 = -0.36320; 'Foral, with efavirenz' = 0.533 (95% CI 0.427, 0.639)

    # ---- Between-STUDY variability (Adeojo 2024 Supplementary Table S2, 'Point Estimate
    #      Inter-Study Variability' column) ----
    # This is a STUDY-level random effect on a meta-analytic fit to trial mean profiles, not
    # a between-subject eta; see the pbpk-qsp-mbma convention for eta_study_<param>. Oral
    # bioavailability is the only parameter carrying it -- the paper states the limited
    # between-study data did not permit estimating inter-study variability of any other
    # parameter (Section 2.2). The 21.2% is read as a log-normal CV, so the internal
    # variance is omega^2 = log(1 + CV^2) = log(1 + 0.212^2) = 0.043963.
    eta_study_lfdepot ~ 0.043963  # Supp Table S2 'Foral' inter-study variability = 21.2% (95% CI 0, 30.1)

    # ---- Residual error (Adeojo 2024 Supplementary Table S2, 'Intra-study variability in
    #      mean concentrations' rows) ----
    # Because the observations are trial MEAN concentrations, this residual describes the
    # scatter of the observed study means about the model prediction, not the scatter of
    # individual patient concentrations.
    addSd  <- 0.214; label("Additive residual SD on study mean concentrations (ng/mL)")            # Supp Table S2 'Additive component' = 0.214 ng/mL (95% CI 0.157, 0.258)
    propSd <- 0.084; label("Proportional residual SD on study mean concentrations (fraction)")     # Supp Table S2 'Proportional component' = 8.4% (95% CI 0, 12.2)
  })

  model({
    # Amounts are carried in mg and volumes in L, so central / vc is mg/L. The source
    # reports levonorgestrel concentrations in ng/mL (implant concentrations in pg/mL),
    # and the residual-error magnitudes above are on the ng/mL scale, so convert:
    # 1 mg/L = 1000 ng/mL.
    mgPerLToNgPerMl <- 1000

    # ---- Derived multiplicative covariate factors (both evaluate to 1 without efavirenz) ----
    cl_efv_factor <- 1 + e_conmed_efv_cl     * CONMED_EFV
    f_efv_factor  <- 1 + e_conmed_efv_fdepot * CONMED_EFV

    # ---- Typical PK parameters ----
    # No between-subject variability: the fit is to study-level mean profiles.
    ka <- exp(lka)
    cl <- exp(lcl) * cl_efv_factor
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # ---- Two-compartment micro-constants ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- ODE system ----
    # Oral doses go to depot and are absorbed first-order. Intravenous doses and the
    # subdermal-implant zero-order release are administered directly into central (the
    # implant release rates in Adeojo 2024 Table 1 are stated as the amounts REACHING THE
    # SYSTEMIC CIRCULATION, footnote a, so no bioavailability term applies to that route).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1

    # ---- Bioavailability ----
    # Applied to the depot compartment only, so it acts on oral doses alone; this is the
    # paper's Foral. The between-study eta and the efavirenz effect both act here. The
    # theta + eta sum is kept on its own simple line so nlmixr2's mu-referencing can
    # detect it; embedding it directly in the f(depot) assignment defeats that.
    fdepot   <- exp(lfdepot + eta_study_lfdepot) * f_efv_factor
    f(depot) <- fdepot

    # ---- Observation and residual error ----
    Cc <- central / vc * mgPerLToNgPerMl
    Cc ~ add(addSd) + prop(propSd)
  })
}
