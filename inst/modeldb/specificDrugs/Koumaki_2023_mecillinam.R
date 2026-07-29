Koumaki_2023_mecillinam <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous mecillinam",
    "(amdinocillin) in healthy adults, used to determine systemic MIC",
    "breakpoints against Enterobacterales for intermittent, extended, and",
    "continuous infusion regimens (Koumaki 2023). CL, Q, Vc and Vp are",
    "parameterised per kilogram of body weight (linear exponent 1, not",
    "allometric). The unbound concentration Cu = fu * Cc (fu = 0.9, i.e.",
    "10 percent protein binding) is the driver of the PK/PD target",
    "fT>MIC >= 40 percent of the dosing interval used for probability of",
    "target attainment. Parameters were obtained by MCMC (Stan)",
    "reanalysis of published AGGREGATE concentration mean and SD data,",
    "not individual-level data; no residual error was estimated."
  )
  reference <- paste(
    "Koumaki V, Dokoumetzidis A, Angelerou MGF, Baka S, Balakrishnan I,",
    "Tsakris A. (2023). Pharmacokinetic/pharmacodynamic determination of",
    "systemic MIC breakpoints for intermittent, extended, and continuous",
    "infusion dosage regimens of mecillinam.",
    "Microbiology Spectrum 11(2):e03441-22.",
    "doi:10.1128/spectrum.03441-22.",
    "The underlying concentration-time data (aggregate mean and SD, 12",
    "healthy volunteers, single 10 mg/kg 15-min IV infusion) are from",
    "Gambertoglio JG, Barriere SL, Lin ET, Conte JE Jr. (1980).",
    "Pharmacokinetics of mecillinam in healthy subjects.",
    "Antimicrob Agents Chemother 18:952-956. doi:10.1128/AAC.18.6.952",
    "(Koumaki 2023 reference 8). The aggregate-data estimation",
    "methodology is Karakitsios E, Dokoumetzidis A (2019) PAGE 28 abstr",
    "8895 (Koumaki 2023 reference 9).",
    sep = " "
  )
  vignette <- "Koumaki_2023_mecillinam"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear per-kg scaling of CL, Q, Vc and Vp (weight exponent 1,",
        "NOT 0.75 allometry): Koumaki 2023 Table 1 reports every",
        "structural parameter already normalised per kilogram",
        "(mL/min/kg and mL/kg), and the Monte Carlo section states the",
        "model 'was parameterized per kilogram'. Do not retrofit",
        "allometric exponents unless re-fitting on the original data.",
        "Because dosing is in mg (not mg/kg), body weight is the only",
        "source of between-subject exposure differences beyond the IIV",
        "terms. Koumaki 2023 drew weights from lognormal distributions",
        "fitted to the U.S. NCHS 2015-2018 adult weight percentiles",
        "(reference 11), stratified by sex."
      ),
      source_name        = "WT"
    )
  )

  # Sex is documented by the source but is NOT a covariate in the PK model:
  # Koumaki 2023 used it only to stratify the body-weight distribution that
  # fed the Monte Carlo simulation (5,000 females + 5,000 males), and no
  # sex effect on any PK parameter was estimated.
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Not a PK covariate. Koumaki 2023 'Monte Carlo simulations'",
        "used sex only to select which NCHS body-weight distribution to",
        "sample from (5,000 females and 5,000 males); no sex effect on",
        "CL, Q, Vc or Vp was estimated or reported. Sex enters the",
        "simulated exposure exclusively through WT."
      ),
      source_name        = "SEX"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12,
    n_studies      = 1,
    age_range      = "Adults; not reported in Koumaki 2023 (see Gambertoglio 1980 for the original cohort).",
    weight_range   = "Not reported in Koumaki 2023; the original dose was 10 mg/kg, so body weight was recorded but the range was not carried into the reanalysis.",
    sex_female_pct = "Not reported in Koumaki 2023.",
    race_ethnicity = "Not reported in Koumaki 2023.",
    disease_state  = "Healthy volunteers. The model was subsequently applied by Monte Carlo simulation to patients with systemic infections (urosepsis, pyelonephritis, bacteraemia) caused by Enterobacterales, including ESBL- and carbapenemase-producing strains.",
    dose_range     = "Estimation data: single 10 mg/kg dose as a 15-min IV infusion. Simulated regimens: 1,000 mg TID, 1,000 mg QID and 1,200 mg QID as 20-min, 2-h and 4-h infusions, plus continuous infusion of 150-4,800 mg/day.",
    regions        = "United States (Gambertoglio 1980 source data, San Francisco); analysis performed in Greece and the United Kingdom.",
    notes          = paste(
      "IMPORTANT: parameters were NOT estimated from individual-level",
      "data. Koumaki 2023 'Population pharmacokinetic model' reanalysed",
      "the AGGREGATE concentration mean and SD values published by",
      "Gambertoglio 1980 (n = 12 healthy volunteers) using the",
      "Karakitsios & Dokoumetzidis (2019) aggregate-data method: at each",
      "MCMC step a large virtual population was simulated from lognormal",
      "parameter distributions and the simulated mean and SD profiles",
      "were fitted to the observed aggregate mean and SD. Estimation used",
      "Stan via RStan 2.19.2 in R 3.5.1. Consequently the reported",
      "'interindividual variability' absorbs BOTH true between-subject",
      "variability and residual/assay error, and no separate residual",
      "error term exists (see the propSd note in ini()).",
      "Koumaki 2023 Table 1 gives RSEs for every estimate. Mecillinam",
      "was assumed to follow linear PK (Monte Carlo section)."
    )
  )

  ini({
    # ---- Structural parameters (Koumaki 2023 Table 1) ----
    # Table 1 reports per-kilogram values in mL/min/kg (clearances) and
    # mL/kg (volumes). They are converted here to the model's L/h and L
    # units; the model() block then multiplies by WT (exponent 1).
    #   clearance: mL/min/kg * 60 / 1000 = L/h/kg
    #   volume:    mL/kg          / 1000 = L/kg
    lcl <- log(0.207);   label("Weight-normalized clearance (L/h/kg)")                   # Table 1: CL = 3.45 mL/min/kg (RSE 1.8%) -> 3.45*60/1000 = 0.207
    lvc <- log(0.12312); label("Weight-normalized central volume of distribution (L/kg)") # Table 1: V1 = 123.12 mL/kg (RSE 16.1%) -> 123.12/1000
    lq  <- log(0.4044);  label("Weight-normalized intercompartmental clearance (L/h/kg)") # Table 1: Q = 6.74 mL/min/kg (RSE 43.4%) -> 6.74*60/1000 = 0.4044
    lvp <- log(0.08056); label("Weight-normalized peripheral volume of distribution (L/kg)") # Table 1: V2 = 80.56 mL/kg (RSE 21.1%) -> 80.56/1000

    # ---- Plasma protein binding ----
    # Koumaki 2023 Monte Carlo section: "Protein binding was assumed to be
    # 10% (12)" (reference 12 = pivmecillinam SmPC). Assumed, not estimated,
    # hence fixed(). fu drives the free-drug PK/PD target %fT>MIC.
    fu <- fixed(0.9); label("Fraction unbound in plasma (assumed; 10% protein binding)")  # Monte Carlo section, ref 12

    # ---- Inter-individual variability (Koumaki 2023 Table 1) ----
    # Table 1 reports IIV as a percentage of a lognormal distribution
    # ("drawn from lognormal distributions with known mean and standard
    # deviation"), i.e. a CV%. Converted to the internal log-scale
    # variance via omega^2 = log(CV^2 + 1).
    etalcl ~ 0.0103503   # Table 1: IIV on CL = 10.2% (RSE 17.1%) -> log(1 + 0.102^2)
    etalvc ~ 0.0744431   # Table 1: IIV on V1 = 27.8% (RSE 35.6%) -> log(1 + 0.278^2)
    etalvp ~ 0.1218636   # Table 1: IIV on V2 = 36.0% (RSE 35.6%) -> log(1 + 0.360^2)
    # No IIV on Q: Koumaki 2023 Table 1 leaves the Q variability cell blank
    # and the Results state "intercompartmental clearance (Q), 6.74 mL/min/kg
    # where no IIV was estimated". Q is therefore a fixed-effect-only
    # parameter; an eta is deliberately NOT introduced.

    # ---- Residual error ----
    # NOT reported and structurally absent: the aggregate-data method fitted
    # the observed concentration mean AND SD, attributing the entire observed
    # spread to the lognormal parameter distributions above. No separate
    # residual/assay error term was estimated. Fixed at zero per the
    # standing policy for unreported RUV, so the model loads and simulates
    # deterministically given the etas. Users re-fitting on individual-level
    # data must replace this with an estimated proportional (or combined)
    # error term. See vignette Errata.
    propSd <- fixed(0); label("Proportional residual error (fraction; FIXED AT ZERO - not estimated in source)")  # Table 1 / Methods: no residual error term
  })
  model({
    # 1. Individual PK parameters. Every structural parameter scales
    #    LINEARLY with body weight (exponent 1) because Table 1 is reported
    #    per kilogram. Q carries no eta (see ini()).
    cl <- exp(lcl + etalcl) * WT
    vc <- exp(lvc + etalvc) * WT
    q  <- exp(lq)           * WT
    vp <- exp(lvp + etalvp) * WT

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 3. ODE system. Mecillinam is given intravenously only (bolus,
    #    intermittent/extended infusion, or continuous infusion), so doses
    #    go directly to `central`; there is no depot.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 4. Observations. Cc is total plasma mecillinam (mg/L; dose in mg,
    #    volume in L). Cu is the unbound concentration that drives the
    #    %fT>MIC PK/PD target; it is a deterministic transform of Cc and is
    #    not separately measured, so it carries no residual error.
    Cc <- central / vc
    Cu <- fu * Cc

    Cc ~ prop(propSd)
  })
}
