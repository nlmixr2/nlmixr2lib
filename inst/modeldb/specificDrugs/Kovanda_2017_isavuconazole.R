Kovanda_2017_isavuconazole <- function() {
  description <- paste(
    "One-compartment population PK model for isavuconazole (dosed as the prodrug",
    "isavuconazonium sulfate, doses expressed as isavuconazole) in adults with invasive",
    "aspergillosis, mucormycosis or other invasive fungal disease from the phase 3 SECURE",
    "and VITAL trials, with first-order oral absorption after a lag time, i.v. infusion",
    "into the central compartment and linear elimination (Kovanda 2017). Oral",
    "bioavailability differs between records with and without mucositis. Fitted",
    "non-parametrically in Pmetrics; the Table 2 medians are encoded as lognormal medians",
    "and the tabulated CV percentages as independent lognormal marginal variances.",
    sep = " "
  )
  reference <- paste(
    "Kovanda LL, Marty FM, Maertens J, Desai AV, Lademacher C, Engelhardt M, Lu Q,",
    "Hope WW. Impact of Mucositis on Absorption and Systemic Drug Exposure of",
    "Isavuconazole. Antimicrobial Agents and Chemotherapy. 2017;61(6):e00101-17.",
    "doi:10.1128/AAC.00101-17. PMCID PMC5444187.",
    sep = " "
  )
  vignette <- "Kovanda_2017_isavuconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Kovanda 2017 Figure 3: compartment 1 is the gut
  # (oral 'Nbolus (1)') and compartment 2 the central compartment receiving
  # the i.v. infusion ('RATEIV (1)').
  compartmentData <- list(
    depot = list(
      analyte = "isavuconazole",
      units = "mg",
      specimen = "administration site",
      verified = TRUE,
      notes = "Oral doses only (Figure 3, compartment 1, carrying Tlag, Ka and F). Intravenous doses go directly into central."
    ),
    central = list(
      analyte = "isavuconazole",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = "Methods, 'Plasma PK sampling': isavuconazole plasma concentrations measured by a validated LC-MS/MS method."
    )
  )

  covariateData <- list(
    DIS_MUCOSITIS = list(
      description = "Mucositis indicator (1 = the record falls in an episode of mucositis or stomatitis, 0 = no mucositis)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no mucositis; the 194 nonmucositis patients of Table 1)",
      notes = paste(
        "Methods, 'PPK modeling': 'The presence of mucositis (yes = 1; no = 0) was used",
        "as a covariate on oral bioavailability (F) as a secondary equation, which took",
        "the following form: F = F1 x (1 - MUC) + F12 x MUC'. It acts on the oral",
        "bioavailability only; clearance, volume, Ka and the lag time are shared.",
        "Patients were counted as mucositis patients only when oral dosing and the",
        "plasma samples coincided with the episode of mucositis, so the indicator is",
        "time-varying within a subject; it matters only at oral dose records.",
        "Mucositis was identified from case-report-form medical-history and",
        "adverse-event MedDRA terms (mucosal inflammation, radiation mucositis,",
        "stomatitis, gastrointestinal inflammation); severity was not graded."
      ),
      source_name = "MUC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 250,
    n_studies = 2,
    age_range = "18-92 years (median 52; mucositis median 50, nonmucositis median 52; Table 1)",
    weight_range = "mean 70.0 +/- 18.3 kg (mucositis 71.7 +/- 18.1, nonmucositis 69.5 +/- 18.4; Table 1)",
    sex_female_pct = 38,
    race_ethnicity = "White 197 (79%), Asian 38 (15%), Black 10 (4%), Other 5 (2%) (Table 1).",
    disease_state = paste(
      "Adults with invasive aspergillosis or other filamentous-fungal disease (SECURE,",
      "NCT00412893) or with mucormycosis, other rare invasive fungal disease, or invasive",
      "aspergillosis with renal impairment (VITAL, NCT00634049) who received oral",
      "isavuconazonium sulfate with plasma samples available. 56 had mucositis during",
      "oral therapy (89.3% with a haematological malignancy, 78.2% neutropenic at",
      "baseline, 26.8% after allogeneic HSCT) and 194 did not (Table 1)."
    ),
    dose_range = paste(
      "Loading 372 mg isavuconazonium sulfate (200 mg isavuconazole) every 8 h for the",
      "first 48 h (i.v. in SECURE, i.v. or oral in VITAL), then 372 mg (200 mg",
      "isavuconazole) once daily i.v. or oral for up to 84 days (SECURE) or 180 days",
      "(VITAL). Route at the investigator's discretion. Median total duration 82 days",
      "(i.v. 7.5 days, oral 73 days; Table 1)."
    ),
    notes = paste(
      "Trough sampling 24 h after the previous dose on treatment days 7, 14 and 42 and",
      "at the end of therapy; full 24-h profiles in a subset of 43 patients (6 with",
      "mucositis). Nonparametric fit in Pmetrics v1.4.1; observed-versus-median",
      "posterior predicted r2 = 0.813 after the Bayesian step; bias 0.11, imprecision",
      "0.938 (Results, 'Population PK (PPK) model')."
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural parameters: the MEDIAN column of Kovanda 2017 Table 2
    # ('Median parameter estimates from the PPK model'), nonmucositis group,
    # which is the MUC = 0 reference of the published F equation. The
    # mucositis-group medians of the shared parameters are identical (Ka,
    # CL, lag) or within 0.6% (V 347.7 L), consistent with one Pmetrics
    # population density summarised by group. Table 2 labels CL and V as
    # 'CL/F' and 'V/F', but Figure 3 routes i.v. infusions straight into
    # the central compartment and F is the absolute oral bioavailability
    # (Abstract: 98.3% and 99.8%), so CL and V are systemic.
    # ------------------------------------------------------------------------
    lka <- log(7.9); label("First-order absorption rate constant (1/h)")
    # Table 2, 'Ka (h-1)' nonmucositis median = 7.9 (mean 6.5 +/- 3.0, range 0.0-8.0, CV 46%)
    lcl <- log(1.9); label("Clearance (L/h)")
    # Table 2, 'CL/F (liters/h)' nonmucositis median = 1.9 (mean 2.3 +/- 1.1, range 0.1-5.9, CV 47%)
    lvc <- log(349.8); label("Central volume of distribution (L)")
    # Table 2, 'V/F (liters)' nonmucositis median = 349.8 (mean 354.1 +/- 182.5, range 5.8-895.5, CV 52%)
    ltlag <- log(1.0); label("Oral absorption lag time (h)")
    # Table 2, 'Lag time (h)' nonmucositis median = 1.0 (mean 1.3 +/- 1.3, range 0.0-5.0, CV 103%)
    lfdepot <- log(0.998); label("Oral bioavailability without mucositis, F1 (fraction)")
    # Table 2, 'F (%)' nonmucositis median = 99.8 (mean 97.4 +/- 6.9, range 70.2-99.9, CV 7%); Abstract '99.8%'

    # Mucositis on F. The paper estimates F12 as a separate parameter
    # (F = F1 x (1 - MUC) + F12 x MUC); it is re-expressed here as the F1
    # baseline times (1 + e * DIS_MUCOSITIS), which is the same function:
    # (0.983 / 0.998) - 1 = -0.01503.
    e_dis_mucositis_fdepot <- -0.01503; label("Fractional effect of mucositis on oral bioavailability (unitless)")
    # Derived from Table 2, 'F (%)' mucositis median = 98.3 (mean 86.0 +/- 18.5, range 50.3-99.7, CV 21%) vs nonmucositis median 99.8; Abstract '98.3% and 99.8%'

    # ------------------------------------------------------------------------
    # Inter-individual variability. Pmetrics estimates a discrete
    # nonparametric joint density, not a parametric omega. Table 2 gives
    # each parameter's CV%, carried here as a lognormal marginal with
    # omega^2 = log(CV^2 + 1) from the nonmucositis (reference) column. No
    # covariances are reported, so the etas are independent. The F CVs are
    # NOT encoded: both medians sit within 2% of the upper bound of 1, so a
    # lognormal marginal around them is physically impossible (F > 1) and a
    # logit-normal one is not identified from a CV. Ka's median sits at the
    # top of its 0-8 1/h search range with the mean below it, so its
    # lognormal marginal is an approximation (see vignette).
    # ------------------------------------------------------------------------
    etalka ~ 0.191942 # Table 2 nonmucositis Ka CV 46% -> log(0.46^2 + 1)
    etalcl ~ 0.199588 # Table 2 nonmucositis CL/F CV 47% -> log(0.47^2 + 1)
    etalvc ~ 0.239332 # Table 2 nonmucositis V/F CV 52% -> log(0.52^2 + 1)
    etaltlag ~ 0.723143 # Table 2 nonmucositis lag time CV 103% -> log(1.03^2 + 1)

    # ------------------------------------------------------------------------
    # Residual error. Methods, 'PPK modeling': 'Data were weighted by the
    # inverse of the estimated assay variance'. Neither the assay-error
    # polynomial nor a gamma/lambda is printed, so the residual error is
    # fixed at zero per the policy for unreported residual error.
    # ------------------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported in source)")
  })
  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # F = F1 x (1 - MUC) + F12 x MUC, with F12 = F1 x (1 + e).
    f(depot) <- exp(lfdepot) * (1 + e_dis_mucositis_fdepot * DIS_MUCOSITIS)
    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
