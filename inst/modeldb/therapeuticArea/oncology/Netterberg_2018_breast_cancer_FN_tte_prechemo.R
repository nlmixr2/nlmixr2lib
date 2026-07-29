Netterberg_2018_breast_cancer_FN_tte_prechemo <- function() {
  description <- "Parametric time-to-event (TTE) submodel for febrile neutropenia (FN) in early breast cancer patients (n=49) receiving adjuvant FEC and docetaxel chemotherapy: the 'prior-to-chemotherapy' variant whose hazard depends only on baseline-available covariates. Time-constant exponential baseline hazard h0 with a log-linear age effect h(t) = h0 * exp(beta1 * (AGE - 54)) (Netterberg 2018 Eq. 3a). The hazard is fixed to 0 for t < 84 h (3.5 days) per the paper: no patient experienced grade 3 neutropenia before that time."
  reference <- paste(
    "Netterberg I, Karlsson MO, Nielsen EI, Quartino AL, Lindman H, Friberg LE.",
    "The risk of febrile neutropenia in breast cancer patients following adjuvant",
    "chemotherapy is predicted by the time course of interleukin-6 and C-reactive",
    "protein by modelling.",
    "Br J Clin Pharmacol. 2018;84(3):490-500.",
    "doi:10.1111/bcp.13477.",
    "Companion biomarker model:",
    "modellib('Netterberg_2018_breast_cancer_FN_biomarkers').",
    sep = " "
  )
  vignette <- "Netterberg_2018_breast_cancer_FN"
  units <- list(
    time = "hour",
    dosing = "n/a (no drug doses; the cycle anchors t = 0 at the start of chemotherapy)",
    concentration = "probability (the model output `sur` is a survival probability for FN, not a drug concentration)"
  )

  covariateData <- list(
    AGE = list(
      description = "Subject age at enrolment.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per subject. Used as a centred deviation `(AGE - 54)` inside the hazard (Netterberg 2018 Eq. 3a). Reference value 54 years is the cohort median (Table 1; range 31-73). Effect coefficient beta1 = 0.0754 per year (Table 2). The Discussion notes that a 70-year-old patient has a 3.3x higher predicted FN risk than a 54-year-old (exp(0.0754 * (70 - 54)) = 3.34).",
      source_name = "Age"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 49L,
    n_studies = 1L,
    age_range = "31-73 years; median 54 years (Netterberg 2018 Table 1)",
    weight_range = "54-111 kg; median 70 kg (Netterberg 2018 Table 1)",
    sex_female_pct = 100,
    disease_state = "Early breast cancer; receiving adjuvant chemotherapy.",
    dose_range = "Most (n=39) received three cycles of FEC (epirubicin 75 mg/m^2 + 5-FU 600 mg/m^2 + cyclophosphamide 600 mg/m^2) followed by three cycles of docetaxel 80 mg/m^2 Q3W. Six patients received the reverse order; small numbers received variant regimens; trastuzumab was added per local routine care when applicable.",
    regions = "Sweden (Uppsala University Hospital, ethics approval Dnr 2006/353).",
    notes = "11 patients developed FN; 12 FN episodes overall (6 in cycle 1, 6 in cycle 4). One patient developed FN in both cycles 1 and 4. All cycle-1 FN events occurred during FEC; all cycle-4 FN events occurred during docetaxel. Only one event per cycle was allowed per patient (Netterberg 2018 Methods 'TTE model for development of FN'). Hazard delayed: h(t) = 0 for t < 84 h (3.5 days), reflecting that no patient had grade 3 neutropenia before that time."
  )

  ini({
    # ----- Baseline hazard (Eq. 3a; Table 2 prior-to-chemotherapy column) -----
    lh0   <- log(5.70e-3); label("Log baseline FN hazard h0 (1/h)")               # Netterberg 2018 Table 2 prior-to-chemo h0 = 5.70e-3 (RSE 31%)

    # ----- Log-linear age effect (Eq. 3a) -----
    beta1 <- 0.0754;        label("Log-linear effect of (AGE - 54) on FN hazard (1/year)")  # Netterberg 2018 Table 2 prior-to-chemo beta1 = 0.0754 (RSE 40%)

    # No IIV, no residual error: parametric TTE is fit by exact likelihood with
    # no $OMEGA / $SIGMA estimated in the source NONMEM run. Forward simulation
    # exposes `hazard` and `sur = exp(-cumhaz)` as derived outputs.
  })

  model({
    # Hazard delay: 0 until 84 h (3.5 days). The Heaviside-style switch keeps
    # the integral well-behaved at the boundary.
    h0    <- exp(lh0)
    haz_raw <- h0 * exp(beta1 * (AGE - 54))
    hazard  <- ifelse(t < 84, 0, haz_raw)

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)
  })
}
