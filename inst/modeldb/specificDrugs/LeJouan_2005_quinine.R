# Population pharmacokinetic model for oral quinine in Cameroonian children
# with uncomplicated Plasmodium falciparum malaria (Le Jouan 2005,
# Antimicrob Agents Chemother 49(9):3658-3662;
# doi:10.1128/aac.49.9.3658-3662.2005). One-compartment first-order
# absorption with linear-in-time evolution of the unbound fraction (fu),
# linear-in-weight typical CL/F and V/F, and a published Bayesian inverse-
# Hill PD model relating parasite-eradication time to average exposure;
# only the PK layer is encoded here (see vignette Assumptions and deviations
# for the PD posterior summary).

LeJouan_2005_quinine <- function() {
  description <- paste(
    "Population PK model for oral quinine in Cameroonian children",
    "(aged 0.55-6.7 years) with uncomplicated Plasmodium falciparum",
    "malaria (Le Jouan 2005). One-compartment with first-order",
    "absorption, time-varying free fraction fu = 0.15 + 0.001*(t - 36)",
    "anchored at the literature value fu(t=36 h) = 0.15 (Babalola 1989)",
    "and clamped to its t = 72 h value beyond the studied window, and",
    "linear-in-body-weight apparent clearance CL/F = fu * 0.53 * WT and",
    "apparent volume V/F = fu * (57 + 3.8 * WT). Doses are oral quinine",
    "base in mg.",
    sep = " "
  )
  reference <- paste(
    "Le Jouan M, Jullien V, Tetanye E, Tran A, Rey E, Treluyer J-M,",
    "Tod M, Pons G (2005). Quinine pharmacokinetics and",
    "pharmacodynamics in children with malaria caused by Plasmodium",
    "falciparum. Antimicrobial Agents and Chemotherapy 49(9):3658-3662.",
    "doi:10.1128/aac.49.9.3658-3662.2005.",
    sep = " "
  )
  vignette <- "LeJouan_2005_quinine"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject at admission. Le Jouan 2005 enrolled 30",
        "Cameroonian children with mean (SD) body weight 13.6 (3.8) kg",
        "(Results, Patient characteristics). The paper reports a linear-",
        "in-WT structural form for CL/F and V/F:",
        "CL/F = fu * (theta1 + theta5 * WT) with theta1 = 0 (fixed",
        "because the intercept was not significantly different from 0;",
        "Results) and theta5 = 0.53 L/h/kg (Table 2);",
        "V/F = fu * (theta2 + theta6 * WT) with theta2 = 57 L (Table 2)",
        "and theta6 = 3.8 L/kg (Table 2). The model encodes these forms",
        "directly inside model() with theta2, theta6, and the typical",
        "weight WT_ref = 15 kg (used by the paper for the typical-child",
        "reference values reported in the Pharmacokinetics-of-quinine",
        "narrative) baked into the scaling expression. The Le Jouan",
        "cohort weight range was 8-23 kg (approximated from age 0.55-6.7",
        "years and mean 13.6 +/- 3.8 kg); extrapolation outside this",
        "range is not supported by the source data.",
        sep = " "
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 30L,
    n_studies       = 1L,
    age_range       = "0.55-6.7 years (Results, Patient characteristics)",
    age_median      = "2.8 years (mean; SD 1.7)",
    weight_range    = "approx. 8-23 kg (cohort mean 13.6 +/- 3.8 kg)",
    weight_median   = "13.6 kg (cohort mean)",
    sex_female_pct  = 100 * 13 / 30,
    race_ethnicity  = c(Black = 100),
    disease_state   = paste(
      "Uncomplicated Plasmodium falciparum malaria; initial parasitaemia",
      "1404-176000 / uL (median 16500 / uL). Inclusion required age",
      "0.5-6 years, infrequent vomiting, first dose within 14 h of",
      "diagnosis, and expected hospital stay of at least 5 days.",
      "Exclusions included cerebral malaria, heart failure, plasma",
      "creatinine above the age-adjusted mean + 2 SD, transaminase >",
      "2x normal, recent antimalarials, and recent enzyme inducers or",
      "inhibitors."
    ),
    dose_range      = paste(
      "Oral quinine base 8.3 mg/kg every 8 h for 5 days (15 doses) as",
      "a 2% formiate-salt syrup (Methods, Study design). Doses for the",
      "model must be supplied in quinine-base mg, matching the paper's",
      "salt-conversion convention."
    ),
    regions         = "Cameroon (Pediatric Unit of Yaounde Central Hospital)",
    notes           = paste(
      "Demographics from Results, Patient characteristics. Sex was 17",
      "boys and 13 girls. The follow-up rate was 100% and no",
      "quinine-attributable side effects were recorded. Sampling: 0, 1,",
      "2, 3, 4, 8, 24, 48, and 56 h after onset of treatment (9 samples",
      "per patient, with 5-9 measured concentrations per patient",
      "available for the population PK analysis). The bioanalytical",
      "assay was liquid chromatography with fluorescence detection",
      "after liquid-liquid extraction; the limit of quantification was",
      "1 mg/L, interday CV < 10%, and bias < 5%. The paper also reports",
      "a Bayesian PD model for the time to a 4-log reduction in",
      "parasitaemia (T_er) as an inverse-Hill function of average",
      "quinine plasma concentration C_av = AUC(0-72 h) / 72: T_er =",
      "T_min * (1 + (C50 / C_av)^s), with WinBugs posterior means",
      "T_min = 35 h (CV 45%, 5-95th percentile 32-38), C50 = 6.6 mg/L",
      "(CV 67%, 5-95th percentile 4.3-9.2), and sigmoidicity s = 2",
      "(fixed; selected from {1, 2, 3, 4} as the best fit). The PD",
      "layer is documented in the vignette and not encoded in this",
      "model file -- see Assumptions and deviations for the rationale."
    )
  )

  ini({
    # Structural population parameters are anchored at the paper's typical
    # 15-kg child at reference time t_ref = 36 h, where the literature
    # value of the unbound fraction is fu(t=36 h) = 0.15 (Le Jouan 2005
    # Methods, citing Babalola 1989 = paper ref 4). At that reference:
    #   CL/F_ref = fu_ref * theta5 * WT_ref = 0.15 * 0.53 * 15 = 1.1925 L/h
    #   V/F_ref  = fu_ref * (theta2 + theta6 * WT_ref)
    #            = 0.15  * (57 + 3.8 * 15) = 0.15 * 114 = 17.10 L
    # The linear-in-WT and linear-in-fu scaling is applied inside model().
    lka  <- log(0.934)   ; label("Absorption rate constant ka (1/h)")                  # Le Jouan 2005 Table 2: Ka = 0.934 (SE 0.244)
    lcl  <- log(1.1925)  ; label("CL/F at typical 15-kg child, fu = 0.15 (L/h)")        # Derived from Table 2 theta5 = 0.53 L/h/kg with theta1 = 0 (fixed) and fu_ref = 0.15
    lvc  <- log(17.10)   ; label("V/F at typical 15-kg child, fu = 0.15 (L)")           # Derived from Table 2 theta2 = 57 L and theta6 = 3.8 L/kg with fu_ref = 0.15

    # Linear time-varying free fraction (Le Jouan 2005 Methods, eq. 1):
    # fu = b * (t - 36) + 0.15. The slope b is fixed at the typical value
    # 0.001/h because the assay did not measure the unbound fraction
    # directly (Results, paragraph 1: "the typical value of b ... could
    # not be estimated and was fixed to 0.001/h to be consistent with the
    # available data"). The 37% CV inter-individual variability on b is
    # retained via etalbfu.
    lbfu <- fixed(log(0.001)) ; label("Slope of fu vs time (1/h, fixed)")              # Le Jouan 2005 Table 2: b = 0.001/h (fixed)

    # Inter-individual variability. Le Jouan 2005 Table 2 reports CV%
    # ("Interindividual CV (%)" column). The internal log-normal
    # variances are recovered by omega^2 = log((CV/100)^2 + 1):
    #   CL  CV 33%  -> log(0.33^2 + 1)  = 0.103376
    #   V   CV 32%  -> log(0.32^2 + 1)  = 0.097488
    #   ka  CV 113% -> log(1.13^2 + 1)  = 0.822873
    #   b   CV 37%  -> log(0.37^2 + 1)  = 0.128335
    # CL/F-V/F correlation = 0.64 (Table 2 footnote b):
    #   cov(etalcl, etalvc) = 0.64 * sqrt(0.103376 * 0.097488) = 0.064249
    etalcl + etalvc ~ c(0.103376, 0.064249, 0.097488)                                  # Le Jouan 2005 Table 2: CL/F CV 33%, V/F CV 32%, corr = 0.64
    etalka          ~ 0.822873                                                          # Le Jouan 2005 Table 2: Ka CV 113%
    etalbfu         ~ 0.128335                                                          # Le Jouan 2005 Table 2: b CV 37% (typical b fixed; only IIV estimated)

    # Residual error. Le Jouan 2005 Methods describes the residual as
    # C_obs = C_pred * exp(epsilon) with var(epsilon) = 0.048 (Table 2:
    # "Variance(epsilon) 0.048 (0.011)"); this is the log-additive form
    # that maps to proportional error in nlmixr2's linear-concentration
    # space, with propSd = sqrt(0.048) = 0.219 (matching the paper's
    # narrative "CV of the residual error was 22%"). The convention
    # follows references/parameter-names.md 'Residual error' and is
    # consistent with the Kloprogge_2014_quinine model in this package.
    propSd <- sqrt(0.048) ; label("Proportional residual SD on linear concentration scale") # Le Jouan 2005 Table 2: Variance(epsilon) = 0.048 (SE 0.011)
  })

  model({
    # Time-varying free fraction (Le Jouan 2005 Methods, eq. 1):
    #   fu(t) = fu_ref + bfu * (t - t_ref)
    # with fu_ref = 0.15 the literature value at t_ref = 36 h and bfu
    # the (mostly-fixed) per-individual slope. The paper asserts linear
    # evolution only for t in [0, 72] h; beyond t = 72 h fu is held
    # constant at its t = 72 value to prevent unsupported extrapolation
    # of the protein-binding term during the multiple-dose tail.
    fu_ref <- 0.15
    t_ref  <- 36
    t_max  <- 72
    bfu    <- exp(lbfu + etalbfu)
    fu     <- fu_ref + bfu * (min(t, t_max) - t_ref)

    # Individual structural PK parameters with the linear-in-WT, linear-
    # in-fu scaling from Table 2. The numeric constants 57 (theta2) and
    # 3.8 (theta6) are the paper's V/F intercept and per-kg slope; the
    # ratio (57 + 3.8 * WT) / (57 + 3.8 * 15) re-centers V/F at the
    # 15-kg reference embedded in lvc. The CL/F scaling reduces to
    # (WT / 15) because theta1 = 0 (Results, Pharmacokinetics).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (fu / fu_ref) * (WT / 15)
    vc <- exp(lvc + etalvc) * (fu / fu_ref) * (57 + 3.8 * WT) / (57 + 3.8 * 15)

    # One-compartment disposition with first-order absorption.
    kel <- cl / vc
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Quinine plasma concentration in mg/L (== ug/mL). Dose units are
    # quinine-base mg, vc is L, so central/vc has units mg/L.
    Cc <- central / vc

    # Proportional residual error on the linear-concentration scale
    # (Le Jouan 2005's log-multiplicative residual maps to proportional;
    # see propSd comment in ini()).
    Cc ~ prop(propSd)
  })
}
