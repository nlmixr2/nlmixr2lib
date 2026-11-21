Bienert_2025_sufentanil <- function() {
  description <- paste(
    "Two-compartment population PK model for epidural sufentanil in adults",
    "undergoing major abdominal surgery under general anaesthesia with",
    "epidural analgesia (epidural bolus plus continuous epidural infusion",
    "for up to about 3 days). Absorption from the epidural space is itself",
    "two-compartmental: the dose enters an epidural absorption compartment",
    "that drains first order into plasma (KA) and exchanges with a slowly",
    "returning epidural distribution compartment (epidural fat; KA14 out,",
    "KA41 back), which produces flip-flop kinetics. Linear elimination from",
    "the central compartment. Full Bayesian fit (Stan/Torsten); typical",
    "values are posterior means. Log-normal IIV on all seven parameters",
    "with a full correlation matrix; residual error additive on the log",
    "scale with a Student-t distribution."
  )
  reference <- paste(
    "Bienert A, Borsuk-De Moor A, Okunska P, Ber J, Siluk D, Waszczuk M,",
    "Kusza K, Bartkowiak T, Szrama J, Kluzik A, Koszel T, Wiczling P.",
    "Pharmacokinetics of Sufentanil After Epidural Administration During",
    "the Course of Extensive Abdominal Surgery.",
    "Clin Pharmacokinet 2025;64(9):1331-1340. doi:10.1007/s40262-025-01543-y.",
    "Correlations of the random effects are from the authors' public code",
    "and output deposit (https://github.com/wiczling/sufentanil-epidural,",
    "deliv/tables/torsten/3b/parameters-table.csv) for the same final model."
  )
  vignette <- "Bienert_2025_sufentanil"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "pg/mL"
  )

  # The epidural distribution compartment (the paper's A4, "a distribution
  # compartment within the absorption site", attributed to epidural fat in
  # the Discussion) exchanges only with the dosing compartment, never with
  # plasma. No registered compartment has that role: `abs_site` is
  # downstream of the depot and exchanges with central, and `depot2` is a
  # second dosing route. It is therefore whitelisted here.
  paper_specific_compartments <- c("depot_peripheral")

  dosing <- "depot"

  compartmentData <- list(
    depot = list(analyte = "sufentanil", units = "ug", specimen = "administration site", verified = TRUE),
    depot_peripheral = list(analyte = "sufentanil", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sufentanil", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sufentanil", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Explored graphically against the individual random effects (Figure S5); no covariate model was fitted."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Explored graphically against the individual random effects (Figure S6); no covariate model was fitted."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Explored graphically against the individual random effects (Figure S7); no covariate model was fitted."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Explored graphically against the individual random effects (Figure S8); no covariate model was fitted."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 23L,
    n_studies = 1L,
    age_range = "30-83 years",
    age_median = "60 years (IQR 42.5-65)",
    weight_median = "71 kg (IQR 59-80.5)",
    sex_female_pct = 52,
    disease_state = paste(
      "Adults (ASA physical status 2-3) undergoing major oncologic or",
      "non-oncologic abdominal surgery under general anaesthesia combined",
      "with epidural analgesia (epidural catheter at Th8/9 to Th12/L1)."
    ),
    dose_range = paste(
      "Epidural sufentanil: an initial dose with an amide local anaesthetic",
      "(ropivacaine or bupivacaine), followed by a continuous epidural",
      "infusion titrated to clinical response during surgery and for the",
      "following 3 days. Median (IQR) total intra-operative epidural",
      "sufentanil dose 15 (10.9-20.6) ug (Table S1)."
    ),
    regions = "Poland (single centre, Poznan University of Medical Sciences)",
    notes = paste(
      "Study NCT06069219. Demographics from Tables S1 and S2: 12 female /",
      "11 male; ASA 2 in 9 and ASA 3 in 14 patients; general anaesthesia",
      "with sevoflurane (20) or desflurane (3) after IV fentanyl and",
      "propofol induction. Plasma sampled at 5 min, 0.5, 1, 2, 4, 6 and",
      "12 h, then every 24 h, before the end of infusion, and at 3, 5, 20",
      "and 40 min and 1, 2, 6, 12, 36 and 72 h after the infusion ended;",
      "LC-MS/MS assay. Two very high concentrations just after the start",
      "of the infusion were flagged as outliers (Figure 2). Stan/Torsten",
      "full Bayesian estimation, 4 chains x 1000 post-warmup iterations,",
      "weakly informative priors on disposition taken from published IV",
      "sufentanil PK."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Bienert 2025 Table 1, posterior Mean column.
    # The paper's own typical-subject statements use the posterior means
    # (Discussion: 'mean posterior clearance ... 64.9 L/h', induction dose
    # 'Rss/KA (approx Rss/3.3)'), so the means are the point estimates.
    # Dose in ug, CL/Q in L/h, V in L; the Stan code reports
    # concentration as amount / V1 * 1000, i.e. pg/mL.
    # ------------------------------------------------------------------
    lcl <- log(64.9); label("Clearance (L/h)")                                       # Table 1: thetaCL mean = 64.9 L/h
    lq <- log(48.2); label("Inter-compartmental clearance (L/h)")                    # Table 1: thetaQ mean = 48.2 L/h
    lvc <- log(10.2); label("Central volume of distribution (L)")                    # Table 1: thetaV1 mean = 10.2 L
    lvp <- log(506); label("Peripheral volume of distribution (L)")                  # Table 1: thetaV2 mean = 506 L
    lka <- log(3.27); label("Absorption rate constant, epidural space to plasma (1/h)") # Table 1: thetaKA mean = 3.27 1/h
    lka14 <- log(14.9); label("Rate constant, epidural absorption to epidural distribution compartment (1/h)") # Table 1: thetaKA14 mean = 14.9 1/h
    lka41 <- log(0.14); label("Rate constant, epidural distribution compartment back to epidural absorption compartment (1/h)") # Table 1: thetaKA41 mean = 0.14 1/h

    # ------------------------------------------------------------------
    # IIV -- log(P_i) ~ N(log(theta_P), Omega) with
    # Omega = diag(omega) * rho * diag(omega) (Statistical Methods).
    # omega is the SD of the log-scale random effect (Table 1 mean column),
    # so the variance is omega^2 and the covariance is
    # rho_ij * omega_i * omega_j. The posterior-mean correlations rho_ij
    # are not printed in the paper; they are from the authors' deposited
    # output for the same fit (parameters-table.csv, rows rho[i,j]).
    # Order of the 7 x 7 block is the paper's P = [CL, Q, V1, V2, KA,
    # KA14, KA41].
    # ------------------------------------------------------------------
    etalcl + etalq + etalvc + etalvp + etalka + etalka14 + etalka41 ~ c(
      0.49^2,
      0.0315 * 0.49 * 0.45, 0.45^2,
      0.253 * 0.49 * 1.51, -0.0422 * 0.45 * 1.51, 1.51^2,
      -0.018 * 0.49 * 0.41, 0.00239 * 0.45 * 0.41, -0.0451 * 1.51 * 0.41, 0.41^2,
      -0.0125 * 0.49 * 0.36, 0.0184 * 0.45 * 0.36, 0.00317 * 1.51 * 0.36, 0.0258 * 0.41 * 0.36, 0.36^2,
      0.025 * 0.49 * 0.38, 0.0203 * 0.45 * 0.38, -0.0955 * 1.51 * 0.38, -0.024 * 0.41 * 0.38, 0.0758 * 0.36 * 0.38, 0.38^2,
      0.187 * 0.49 * 0.5, 0.106 * 0.45 * 0.5, 0.111 * 1.51 * 0.5, -0.0355 * 0.41 * 0.5, -0.0811 * 0.36 * 0.5, 0.112 * 0.38 * 0.5, 0.5^2
    ) # Table 1: omegaCL 0.49, omegaQ 0.45, omegaV1 1.51, omegaV2 0.41, omegaKA 0.36, omegaKA14 0.38, omegaKA41 0.5; rho from the deposited parameters-table.csv

    # ------------------------------------------------------------------
    # Residual error -- log(Cobs) ~ student_t(nu, log(f), sigma)
    # (Statistical Methods), i.e. additive on the log scale with a
    # Student-t distribution of scale sigma and nu degrees of freedom.
    # rxode2's lnorm(expSd) + dt(nu) uses the same scale parameterisation
    # (the simulated log residual has variance expSd^2 * nu / (nu - 2)).
    # ------------------------------------------------------------------
    expSd <- 0.31; label("Residual scale on the log concentration (Student-t)") # Table 1: sigma mean = 0.31
    nu <- 3.81; label("Degrees of freedom of the Student-t residual error")    # Table 1: nu mean = 3.81
  })
  model({
    # 1. Individual parameters (log-normal IIV on every parameter).
    cl <- exp(lcl + etalcl)
    q <- exp(lq + etalq)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka + etalka)
    ka14 <- exp(lka14 + etalka14)
    ka41 <- exp(lka41 + etalka41)

    # 2. Micro-constants (1/h).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODEs. The paper writes the disposition part in concentrations
    # (V1 dC1/dt = KA*A1 - CL*C1 - Q*C1 + Q*C2); with amounts
    # central = V1*C1 and peripheral1 = V2*C2 that is the system below,
    # identical to the K matrix of the deposited Stan code.
    # depot = A1 (epidural absorption compartment, receives the dose);
    # depot_peripheral = A4 (epidural distribution compartment).
    d/dt(depot) <- -(ka + ka14) * depot + ka41 * depot_peripheral
    d/dt(depot_peripheral) <- ka14 * depot - ka41 * depot_peripheral
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Observation: ug / L * 1000 = pg/mL (the Stan code's '/ V1 * 1000').
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd) + dt(nu)
  })
}
