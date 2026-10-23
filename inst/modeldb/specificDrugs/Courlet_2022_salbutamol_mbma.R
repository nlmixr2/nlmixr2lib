Courlet_2022_salbutamol_mbma <- function() {
  description <- paste(
    "MBMA. One-compartment population PK model of salbutamol (albuterol) in",
    "plasma and urine, fitted jointly to individual and aggregate (study-mean)",
    "data from 13 published studies in healthy volunteers, trained athletes and",
    "asthmatics (Courlet 2022, WADA doping-control analysis). An inhaled dose",
    "splits into a fixed 20% lung fraction absorbed first-order from a lung",
    "depot and an 80% swallowed fraction absorbed first-order from the gut with",
    "the oral bioavailability; an oral dose enters the gut depot only. Urinary",
    "excretion of unchanged salbutamol is a first-order arm of the total",
    "clearance into a cumulative urine compartment, and the urine concentration",
    "is that amount divided by the urine volume produced at a constant rate",
    "since the last bladder voiding. Physical exercise lowers urine production",
    "by 15%; urine-production IIV differs between USG-corrected and uncorrected",
    "records; aggregate records scale every random effect by N^-0.25."
  )
  reference <- paste(
    "Courlet P, Buclin T, Biollaz J, Mazzoni I, Rabin O, Guidi M.",
    "Model-based meta-analysis of salbutamol pharmacokinetics and practical",
    "implications for doping control. CPT Pharmacometrics Syst Pharmacol.",
    "2022;11(4):469-481. doi:10.1002/psp4.12773.",
    sep = " "
  )
  vignette <- "Courlet_2022_salbutamol_mbma"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  compartmentData <- list(
    depot_lung = list(analyte = "salbutamol", units = "ug", specimen = "administration site", verified = TRUE),
    depot_oral = list(analyte = "salbutamol", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "salbutamol", units = "ug", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "salbutamol", units = "ug", specimen = "urine", verified = TRUE),
    urine_vol = list(analyte = "not applicable", units = "L", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    ROUTE_ORAL = list(
      description = "Oral (systemic) administration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (inhalation)",
      notes = paste(
        "Per-dose-record indicator. 1 = oral salbutamol (tablet or solution):",
        "the whole dose enters `depot_oral` with oral bioavailability F2.",
        "0 = inhalation: the dose is given as TWO records with the same `amt`,",
        "one into `depot_lung` (bioavailable fraction F1 = 0.2) and one into",
        "`depot_oral` (bioavailable fraction (1 - F1) x F2 for the swallowed",
        "part), as in Figure 1 and Supplementary Material S1 Equation 6 of",
        "Courlet 2022. The indicator only scales the gut-depot bioavailability."
      ),
      source_name = "administration mode (Table S3)"
    ),
    USG_CORRECTED = list(
      description = "Urine-specific-gravity-corrected urine concentration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (uncorrected urine concentration)",
      notes = paste(
        "Per-record. 1 = the urine concentration was adjusted to a urine",
        "specific gravity of 1.020 g/mL, corrected = (1.020 - 1) / (USG - 1) x",
        "raw (Supplementary Material S1 Equation 1); 0 = uncorrected. Selects",
        "which of the two separately estimated IIVs on urine production applies",
        "(Table 1: 74% CV uncorrected, 38% CV corrected). The paper's",
        "simulations (Figures 4 and 6) report USG-corrected concentrations."
      ),
      source_name = "UR_PROD uncorrected / corrected (Table 1)"
    ),
    EXERCISE = list(
      description = "PK session performed during physical exercise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (at rest, or activity state not reported)",
      notes = paste(
        "1 = PK sampling performed during physical activity, 0 = at rest",
        "(Supplementary Material S1 Equation 4). Missing activity data were",
        "pooled with rest (Results, Covariate analysis). Multiplies urine",
        "production by (1 + e_exercise_urprod), i.e. 0.0467 -> 0.0396 L/h."
      ),
      source_name = "physical / exercise"
    ),
    DTYPE_AGGREGATED = list(
      description = "Aggregate (study-mean) record indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (individual-subject record)",
      notes = paste(
        "1 = the record is an arithmetic-mean concentration of N_ARM subjects",
        "read from a published study (Methods, Data); 0 = an individual",
        "concentration. Selects the residual-error magnitudes: individual",
        "records use the plasma 40% proportional and urine 46% + 18 ng/mL",
        "combined errors; aggregate records use the plasma 77% and urine 131%",
        "proportional errors, each weighted by N_ARM^-0.25 (Table 1)."
      ),
      source_name = "individual vs aggregate data (Table S3)"
    ),
    N_ARM = list(
      description = "Number of individuals contributing to an aggregate PK profile",
      units = "participants",
      type = "count",
      reference_category = NULL,
      notes = paste(
        "Nobs of Courlet 2022 (8 to 30 in the aggregate studies, Table S3).",
        "Set to 1 on individual records. Every random effect of an aggregate",
        "record is weighted as omega^2 = omega_raw^2 / sqrt(Nobs) and",
        "sigma = sigma_raw / sqrt(Nobs) (Table 1 footnote), which with the",
        "Results figures (plasma 46% / 40% at Nobs = 8 / 13 from 77%; urine",
        "78% / 56% at Nobs = 8 / 30 from 131%) is an SD weight of Nobs^-0.25",
        "on both the etas and the residual SDs. See the vignette Errata for",
        "the 1 / sqrt(Nobs) form printed in the Methods equation."
      ),
      source_name = "Nobs"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 13,
    age_range = "Study means 21-29 years; individual data 22-45 years where reported (Table S3)",
    weight_range = "Study means 70-84 kg; individual data 65-89 kg where reported (Table S3)",
    sex_female_pct = NA_real_,
    disease_state = paste(
      "Healthy volunteers, trained / elite athletes and asthmatics (with and",
      "without elite-athlete status); all studies conducted in males or",
      "without sex information (Methods, Covariate exploration)."
    ),
    dose_range = paste(
      "Inhaled 200-1600 ug per administration (median 800 ug), daily",
      "200-1600 ug; oral 4000-12000 ug per administration (median 8000 ug)",
      "(Results, Literature data)."
    ),
    regions = "Not reported (13 studies selected by WADA from the literature and WADA-funded projects)",
    notes = paste(
      "13 studies (Table S3): Pichon 2005 (10, aggregate), Sporer 2007 (30,",
      "aggregate), Sporer 2008 (8, aggregate), Mareck 2011 (28, individual),",
      "Delbeke 2007 (individual, size not reported), Elers 2010 (10), Elers",
      "2011 (20), Elers 2012 (18, aggregate), Pillard 2018 (12, aggregate),",
      "Haase 2015 (13), Jacobson 2013 (27, individual), Backer 2013 (17,",
      "individual) and Hostrup 2014 (9, individual). The reported sizes sum",
      "to 202 excluding Delbeke 2007, so the total number of subjects is not",
      "stated. 121 plasma concentrations (43 individual, 78 mean-profile",
      "points) and 796 urine concentrations (747 individual, 49 mean-profile",
      "points). Urine measured as free or non-sulfated salbutamol, which the",
      "model does not distinguish."
    )
  )

  ini({
    # Absorption. Estimated on plasma data alone (Table S4) and then fixed
    # for the joint plasma + urine model (Methods, Base structural model).
    lfdepot_lung <- fixed(log(0.2)); label("Fraction of an inhaled dose deposited in the lungs, F1 (fraction)") # Table 1: F1 = 0.2 FIX
    logitfdepot_oral <- fixed(-0.328); label("Logit of the oral bioavailability, logit F2 (unitless)") # Table 1: logit F2 = -0.328 FIX; F2 = 0.419
    lka_lung <- fixed(log(31.6)); label("Absorption rate constant after inhalation, ka1 (1/h)") # Table 1: ka1 = 31.6 1/h FIX
    lka_oral <- fixed(log(1.47)); label("Absorption rate constant after oral administration, ka2 (1/h)") # Table 1: ka2 = 1.47 1/h FIX

    # Disposition and urinary excretion
    lvc <- log(205); label("Central volume of distribution, V3 (L)") # Table 1: V3 = 205 L (RSE 2%)
    lkurine <- log(0.0468); label("Urinary excretion rate constant from plasma to urine, k34 (1/h)") # Table 1: k34 = 0.0468 1/h (RSE 14%)
    lcl <- log(28); label("Total salbutamol clearance, CL = (k30 + k34) x V3 (L/h)") # Table 1: CL = 28 L/h (RSE 8%); Supplementary Material S1 Equation 2
    lurprod <- log(0.0467); label("Urine production rate at rest, UR_PROD (L/h)") # Table 1: UR_PROD = 0.0467 L/h

    # Covariate effect
    e_exercise_urprod <- -0.153; label("Fractional change in urine production during physical exercise (unitless)") # Table 1: theta physical,UR_PROD = -0.153 (RSE 43%)

    # IIV. Table 1 reports CV%; variances below are log(CV^2 + 1).
    etalcl ~ 0.14158 # Table 1: omega CL = 39% CV (RSE 12%)
    etalurprod ~ 0.43671 # Table 1: omega UR_PROD uncorrected = 74% CV (RSE 49%)
    etalurprod_usg ~ 0.13491 # Table 1: omega UR_PROD corrected = 38% CV (RSE 16%)

    # Residual error. The aggregate SDs are the unweighted values that are
    # multiplied by N_ARM^-0.25 in model().
    propSd <- 0.40; label("Proportional residual error, plasma, individual records (fraction)") # Table 1: 40% CV (RSE 7%)
    propSd_aggregate <- 0.77; label("Proportional residual error, plasma, aggregate records before N weighting (fraction)") # Table 1: 77% CV (RSE 12%)
    propSd_urineSal <- 0.46; label("Proportional residual error, urine, individual records (fraction)") # Table 1: 46% CV (RSE 2%)
    addSd_urineSal <- 18; label("Additive residual error, urine, individual records (ng/mL)") # Table 1: 18 ng/mL (RSE 29%)
    propSd_urineSal_aggregate <- 1.31; label("Proportional residual error, urine, aggregate records before N weighting (fraction)") # Table 1: 131% CV (RSE 12%)
  })

  model({
    # Random-effect weight for aggregate records: omega^2 / sqrt(N) and
    # sigma / N^0.25 (Table 1 footnote and Results, Model evaluation).
    # N_ARM = 1 on individual records gives a weight of 1.
    wN <- 1 / N_ARM^0.25

    fdepot_lung <- exp(lfdepot_lung)
    foral <- expit(logitfdepot_oral)
    ka_lung <- exp(lka_lung)
    ka_oral <- exp(lka_oral)

    vc <- exp(lvc)
    cl <- exp(lcl + etalcl * wN)
    kurine <- exp(lkurine)
    urprod <- exp(lurprod + (etalurprod * (1 - USG_CORRECTED) + etalurprod_usg * USG_CORRECTED) * wN) *
      (1 + e_exercise_urprod * EXERCISE)

    # CL = (k30 + k34) x V3 (Supplementary Material S1 Equation 2); k30 is the
    # non-urinary elimination of free salbutamol.
    kel <- cl / vc
    k30 <- kel - kurine

    d/dt(depot_lung) <- -ka_lung * depot_lung
    d/dt(depot_oral) <- -ka_oral * depot_oral
    d/dt(central) <- ka_lung * depot_lung + ka_oral * depot_oral - k30 * central - kurine * central
    d/dt(urine) <- kurine * central
    # Urine volume accumulates at the constant production rate; the bladder
    # is voided (urine and urine_vol reset to 0) before each administration
    # (Methods, Base structural model).
    d/dt(urine_vol) <- urprod

    # Inhaled doses go to both depots; the gut depot carries the swallowed
    # (1 - F1) fraction with oral bioavailability F2 (Figure 1).
    f(depot_lung) <- fdepot_lung
    f(depot_oral) <- foral * (ROUTE_ORAL + (1 - ROUTE_ORAL) * (1 - fdepot_lung))

    # ug / L = ng/mL
    Cc <- central / vc
    urineSal <- urine / urine_vol

    propCc <- (1 - DTYPE_AGGREGATED) * propSd + DTYPE_AGGREGATED * propSd_aggregate * wN
    propUrine <- (1 - DTYPE_AGGREGATED) * propSd_urineSal + DTYPE_AGGREGATED * propSd_urineSal_aggregate * wN
    addUrine <- (1 - DTYPE_AGGREGATED) * addSd_urineSal

    Cc ~ prop(propCc)
    urineSal ~ add(addUrine) + prop(propUrine)
  })
}
