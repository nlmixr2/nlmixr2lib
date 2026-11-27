Aruldhas_2019_isoniazid <- function() {
  description <- paste(
    "One-compartment population PK model for oral isoniazid in Indian",
    "children aged 2-16 years with pulmonary or lymph-node tuberculosis",
    "(Aruldhas 2019), sampled at steady state in the intensive phase of",
    "daily or thrice-weekly therapy. Absorption through five fixed transit",
    "compartments with the absorption rate set equal to the transit rate",
    "(ka = ktr = 6 / MTT). Fixed allometric scaling of CL/F (exponent 0.75)",
    "and V/F (exponent 1) on body weight centred at 19.4 kg. A two-class",
    "mixture model on apparent clearance separates slow (2.59 L/h) and fast",
    "(7.79 L/h) acetylators, gated by the MIX_FAST_ELIM class column.",
    "Relative bioavailability fixed to 1 with between-subject variability;",
    "BSV also on V/F and MTT, none on CL/F; additive error on log",
    "concentrations. Companion rifampicin model:",
    "modellib('Aruldhas_2019_rifampicin').",
    sep = " "
  )
  reference <- paste(
    "Aruldhas BW, Hoglund RM, Ranjalkar J, Tarning J, Mathew SK,",
    "Verghese VP, Bose A, Mathew BS. Optimization of dosing regimens of",
    "isoniazid and rifampicin in children with tuberculosis in India.",
    "Br J Clin Pharmacol. 2019;85(3):644-654. doi:10.1111/bcp.13846.",
    sep = " "
  )
  vignette <- "Aruldhas_2019_tuberculosis"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "isoniazid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fixed allometric function on CL/F (exponent 0.75) and V/F",
        "(exponent 1) (Methods, 'Population pharmacokinetic modelling';",
        "Results: estimating the exponent did not improve the fit). Centred",
        "at 19.4 kg, the 'typical child' of the Table 3 footnote. Table 2",
        "reports median (IQR) weight 19.5 (13.7-33.7) kg in the 41 enrolled",
        "children.",
        sep = " "
      ),
      source_name = "WT"
    ),
    MIX_FAST_ELIM = list(
      description = paste(
        "Mixture-model class indicator for isoniazid acetylator status:",
        "1 = fast acetylator (typical CL/F 7.79 L/h at 19.4 kg),",
        "0 = slow acetylator (typical CL/F 2.59 L/h at 19.4 kg)",
        sep = " "
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (slow acetylator)",
      notes = paste(
        "Not a measured covariate: acetylator genotype and phenotype were",
        "unknown, and the class is the latent subpopulation of a NONMEM",
        "mixture model on CL/F (Methods; Results, 'Pharmacokinetic",
        "properties of isoniazid'). The paper states that 31% of children",
        "were fast acetylators (Results; Table 3 footnote; Discussion,",
        "'30.8% were found to be fast acetylators and 69.2% to be slow",
        "acetylators'). The paper's own exposure summaries do not agree",
        "with that split -- see the validation vignette before drawing a",
        "population. For a virtual cohort, draw MIX_FAST_ELIM per subject",
        "from a Bernoulli distribution with the fast-acetylator",
        "probability; for typical-value work, simulate each class",
        "separately, as the paper's Table 5 does.",
        sep = " "
      ),
      source_name = "MIXTURE subpopulation"
    )
  )

  # Screened in the stepwise covariate search and not retained ('No other
  # covariates could be retained in the stepwise covariate approach').
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Methods); not retained."
    ),
    SEXF = list(
      description = "Biological sex, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the stepwise covariate search (Methods); not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Methods); not retained."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Methods); not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Methods); not retained."
    ),
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Evaluated as an allometric size descriptor alternative to body",
        "weight (Methods); body weight was retained.",
        sep = " "
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 39L,
    n_studies = 1L,
    n_observations = 290L,
    age_range = "2-16 years; median 7 (IQR 3.5-13) in the 41 enrolled children (Table 2)",
    weight_range = "median 19.5 kg (IQR 13.7-33.7) in the 41 enrolled children (Table 2)",
    sex_female_pct = 29.3,
    race_ethnicity = "Indian (single centre, Vellore, Tamil Nadu)",
    disease_state = paste(
      "Newly diagnosed pulmonary (n = 36) or lymph-node (n = 5)",
      "tuberculosis, haemoglobin >= 10 g/dL; relapse, multidrug-resistant",
      "tuberculosis and HIV co-infection excluded. Sampled during the last",
      "2 weeks of the 8-week intensive phase (steady state), after an",
      "overnight fast.",
      sep = " "
    ),
    dose_range = paste(
      "Isoniazid with rifampicin, pyrazinamide and ethambutol by the",
      "Revised National Tuberculosis Control Programme weight bands",
      "(Table 1, Regimen 2: 75-300 mg per dose), given daily (WHO,",
      "n = 14) or thrice weekly (RNTCP directly observed therapy,",
      "n = 27).",
      sep = " "
    ),
    regions = "India (Christian Medical College, Vellore, Tamil Nadu)",
    notes = paste(
      "41 children enrolled; isoniazid was not measured in 2 because of an",
      "LC-MS/MS breakdown, leaving 290 concentrations from 39 children",
      "(7% below the 0.01 mg/L LLOQ, imputed as LLOQ/2). Samples pre-dose",
      "and at 0.5, 1, 1.5, 2, 2.5, 4 and 6 h. Acetylator status was not",
      "measured.",
      sep = " "
    )
  )

  ini({
    # Table 3 population estimates are 'for a typical child weighing 19.4 kg
    # with tuberculosis' (Table 3 footnote).

    # Absorption: five fixed transit compartments with ka = ktr (Results:
    # 'kA and kTR were set to be equal since estimating them separately did
    # not improve the model fit'). The mean transit time spans the 5 + 1 = 6
    # first-order transfers from depot to central, so ktr = 6 / MTT. The
    # vignette confirms this reproduces the published median Tmax of 1.10 h.
    lmtt <- log(0.547) ; label("Mean absorption transit time MTT (h)") # Table 3 row 'MTT(h)' = 0.547 (RSE 12.6%)

    # Mixture model on CL/F. The slow-acetylator clearance is the reference
    # (MIX_FAST_ELIM = 0); the fast-acetylator clearance is reached through
    # the log-ratio e_mix_fast_elim_cl = log(7.79 / 2.59).
    lcl <- log(2.59) ; label("Apparent clearance CL/F of slow acetylators (L/h, 19.4-kg child)") # Table 3 row 'CL SLOW /F (l h-1)' = 2.59 (RSE 9.95%)
    e_mix_fast_elim_cl <- log(7.79 / 2.59) ; label("Log ratio of fast- to slow-acetylator CL/F (unitless)") # Table 3 rows 'CL FAST /F (l h-1)' = 7.79 (RSE 7.38%) and 'CL SLOW /F' = 2.59

    lvc <- log(29.7) ; label("Apparent volume of distribution V/F (L, 19.4-kg child)") # Table 3 row 'Vc/F (l)' = 29.7 (RSE 7.66%)

    # 'Relative bioavailability (F) was fixed to unity for the population but
    # with an estimated interindividual variability' (Methods).
    lfdepot <- fixed(log(1)) ; label("Relative oral bioavailability F (unitless)") # Table 3 row 'F (%)' = 100 fix

    # 'Effect of both body weight and fat-free mass were evaluated with
    # allometric scaling on all clearances and volumes parameters with an
    # exponent of 0.75 and 1, respectively' (Methods); estimating the
    # exponent did not improve the fit (Results).
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of body weight on CL/F (unitless)") # Methods, 'Population pharmacokinetic modelling': exponent 0.75
    e_wt_vc <- fixed(1) ; label("Allometric exponent of body weight on V/F (unitless)") # Methods, 'Population pharmacokinetic modelling': exponent 1

    # Between-subject variability. Table 3 reports '%CV for BSV'; converted
    # with omega^2 = log(1 + CV^2), the exact log-normal relation. No BSV was
    # estimated on CL/F (Table 3 BSV column is blank for both CL rows).
    etalvc ~ 0.05328 # Table 3 row 'Vc/F' BSV = 23.4% CV (RSE 14.1%); shrinkage 25.4%
    etalmtt ~ 0.3818 # Table 3 row 'MTT(h)' BSV = 68.2% CV (RSE 13.9%); shrinkage 13.8%
    etalfdepot ~ 0.1609 # Table 3 row 'F (%)' BSV = 41.8% CV (RSE 12.1%); shrinkage 23.5%

    # Residual error: 'additive error on the individually predicted
    # logarithmic concentrations' (Methods), i.e. lnorm() here. The Table 3
    # RUV of 0.0967 is read as the VARIANCE sigma^2: its bootstrap 95% CI
    # (0.0639-0.141) is twice as wide, relative to the estimate, as its
    # printed RSE of 9.79% implies -- the signature of a variance whose RSE
    # is reported on the SD scale. The BSV rows, which are on the CV scale,
    # show no such factor. See the vignette.
    expSd <- sqrt(0.0967) ; label("Additive residual error on log concentration (SD)") # Table 3 row 'RUV' = 0.0967 (RSE 9.79%), read as a variance
  })
  model({
    # Allometric size scaling centred at the 19.4-kg typical child.
    cl <- exp(lcl + e_mix_fast_elim_cl * MIX_FAST_ELIM) * (WT / 19.4)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 19.4)^e_wt_vc
    mtt <- exp(lmtt + etalmtt)

    # Five transit compartments; ka = ktr, so MTT covers 6 transfers.
    ktr <- 6 / mtt
    ka <- ktr
    kel <- cl / vc

    d / dt(depot) <- -ktr * depot
    d / dt(transit1) <- ktr * depot - ktr * transit1
    d / dt(transit2) <- ktr * transit1 - ktr * transit2
    d / dt(transit3) <- ktr * transit2 - ktr * transit3
    d / dt(transit4) <- ktr * transit3 - ktr * transit4
    d / dt(transit5) <- ktr * transit4 - ktr * transit5
    d / dt(central) <- ka * transit5 - kel * central

    f(depot) <- exp(lfdepot + etalfdepot)

    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
