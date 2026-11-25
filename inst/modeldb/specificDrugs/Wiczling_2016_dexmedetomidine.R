Wiczling_2016_dexmedetomidine <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous dexmedetomidine",
    "given as a long-term continuous infusion (45-229 h) to critically ill",
    "children in a paediatric intensive care unit (Wiczling 2016). Fully",
    "Bayesian fit (WinBUGS) with informative priors from a pooled paediatric",
    "analysis; typical values are posterior means. Allometric weight scaling",
    "on CL, Q, V1 and V2 (exponents 0.75 and 1, 70 kg reference) and a Hill",
    "maturation of CL on postmenstrual age. All four disposition parameters",
    "differ between the two sampling occasions (0-24 h after infusion start",
    "and 0-8 h after infusion end) by fractional changes that the authors'",
    "code switches in smoothly at 30 h after infusion start. Full 4 x 4",
    "log-normal IIV block; residual error additive on the log scale with a",
    "Student-t distribution."
  )
  reference <- paste(
    "Wiczling P, Bartkowska-Sniatkowska A, Szerkus O, Siluk D,",
    "Rosada-Kurasinska J, Warzybok J, Borsuk A, Kaliszan R, Grzeskowiak E,",
    "Bienert A. The pharmacokinetics of dexmedetomidine during long-term",
    "infusion in critically ill pediatric patients. A Bayesian approach with",
    "informative priors. J Pharmacokinet Pharmacodyn 2016;43(3):315-324.",
    "doi:10.1007/s10928-016-9474-0."
  )
  vignette <- "Wiczling_2016_dexmedetomidine"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ng/mL"
  )

  dosing <- "central"

  compartmentData <- list(
    central = list(analyte = "dexmedetomidine", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dexmedetomidine", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling (WT/70)^0.75 on CL and Q and (WT/70)^1 on V1 and",
        "V2 (Methods equation 5; WinBUGS code). The code uses the weight at",
        "the subject's first record, i.e. baseline weight."
      ),
      source_name = "BW / weight"
    ),
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the Hill maturation of CL on postmenstrual age (Methods",
        "equation 6). The authors' WinBUGS code computes postmenstrual age in",
        "weeks as 4 * age_months + 40, i.e. a 40-week (term) gestation and",
        "4 weeks per month; model() reproduces that conversion. The code uses",
        "the age at the subject's first record."
      ),
      source_name = "age (months; Table 2 'Age, months')"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Explored graphically against the individual random effects (Figure 7S); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 38L,
    n_studies = 1L,
    n_observations = 470L,
    age_range = "1.4-188.6 months",
    age_median = "70 months",
    weight_range = "4.7-60 kg",
    weight_median = "18.5 kg",
    sex_female_pct = 39.5,
    disease_state = paste(
      "Critically ill children in a paediatric intensive care unit receiving",
      "dexmedetomidine in addition to standard sufentanil + midazolam",
      "sedation: acute respiratory failure (18), severe sepsis or septic",
      "shock (10), multiple or brain trauma (8), acute cardiac insufficiency",
      "(2). PRISM score median 0.5 (0-11). Severe renal or hepatic",
      "insufficiency excluded."
    ),
    dose_range = paste(
      "Continuous IV infusion started at 0.8 ug/kg/h, titrated in steps of",
      "0.2 ug/kg/h to a Cook sedation score of 7-14 (maximum 1.4 ug/kg/h),",
      "then reduced by 0.2 ug/kg/h steps to the minimum before stopping.",
      "Median total dose 1153.8 ug (248.8-4732.2); median infusion duration",
      "97.3 h (45.0-229.2) (Table 2)."
    ),
    regions = "Poland (single centre, Poznan University of Medical Sciences)",
    notes = paste(
      "Sampling: before the infusion and at 1, 4, 8, 12, 16, 20 and 24 h",
      "(occasion 1); just before the infusion stopped and 5, 10, 20 min and",
      "1, 2, 4 and 6 h after (occasion 2). LC-MS/MS assay. WinBUGS 1.4.3 with",
      "informative priors on the typical values and IIV from Potts 2009",
      "(Paediatr Anaesth 19:1119-1129; pooled analysis of 95 children,",
      "Table 1); 3 chains x 6000 iterations, 3000 burn-in, thinning 3."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Typical values -- Table 3, posterior mean column, for a 70-kg
    # subject. Dose in ug, CL/Q in L/h, V in L -> concentration ug/L =
    # ng/mL (Figure 1 and 2 axes).
    # ------------------------------------------------------------------
    lcl <- log(41.6); label("Clearance for a 70-kg adult (L/h)")                          # Table 3: theta CL mean = 41.6 L/h/70 kg (90% HDI 39.0-44.3)
    lq <- log(56.8); label("Inter-compartmental clearance for a 70-kg adult (L/h)")      # Table 3: theta Q mean = 56.8 L/h/70 kg (90% HDI 43.5-73.5)
    lvc <- log(52.0); label("Central volume of distribution for a 70-kg adult (L)")      # Table 3: theta V1 mean = 52.0 L/70 kg (90% HDI 43.2-59.6)
    lvp <- log(70.4); label("Peripheral volume of distribution for a 70-kg adult (L)")   # Table 3: theta V2 mean = 70.4 L/70 kg (90% HDI 63.0-79.8)

    # Maturation of CL: PMA^Hill / (PMA^Hill + TE50^Hill), PMA in weeks
    # (Methods equation 6).
    ltm50_cl <- log(42.5); label("Postmenstrual age at 50% CL maturation, TE50 (weeks)") # Table 3: theta TE50 mean = 42.5 weeks
    lhill_cl <- log(2.45); label("Hill coefficient of the CL maturation (unitless)")      # Table 3: theta Hill mean = 2.45 (90% HDI 1.72-3.39)

    # Allometric exponents -- fixed by theory (Methods equation 5).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT on CL (unitless)") # Methods eq. 5: K = 0.75 for clearance
    e_wt_q <- fixed(0.75); label("Allometric exponent of WT on Q (unitless)")   # Methods eq. 5: K = 0.75 for clearance
    e_wt_vc <- fixed(1); label("Allometric exponent of WT on V1 (unitless)")    # Methods eq. 5: K = 1 for distribution volumes
    e_wt_vp <- fixed(1); label("Allometric exponent of WT on V2 (unitless)")    # Methods eq. 5: K = 1 for distribution volumes

    # ------------------------------------------------------------------
    # Occasion-2 fractional changes. The paper reports fP (occasion 2 /
    # occasion 1 ratio, fP = 1 for occasion 1); here e = fP - 1 so that
    # P(occasion 2) = P * (1 + e). See model() for the time switch.
    # ------------------------------------------------------------------
    e_occ2_cl <- 0.31; label("Fractional change in CL on occasion 2 (unitless)")  # Table 3: f CL mean = 1.31 (90% HDI 0.910-1.82); 1.31 - 1
    e_occ2_q <- 0.02; label("Fractional change in Q on occasion 2 (unitless)")    # Table 3: f Q mean = 1.02 (90% HDI 0.722-1.40); 1.02 - 1
    e_occ2_vc <- 0.50; label("Fractional change in V1 on occasion 2 (unitless)")  # Table 3: f V1 mean = 1.50 (90% HDI 1.33-1.65); 1.50 - 1
    e_occ2_vp <- -0.14; label("Fractional change in V2 on occasion 2 (unitless)") # Table 3: f V2 mean = 0.86 (90% HDI 0.630-1.17); 0.86 - 1

    # ------------------------------------------------------------------
    # IIV -- posterior-mean covariance matrix Omega (Table 2S, lower
    # triangle, order CL, Q, V1, V2 as in the WinBUGS code logtheta[i, 1:4]).
    # Check against Results: CV = sqrt(exp(omega^2) - 1) = 56, 83, 152,
    # 68 %; corr(Q, V1) = 0.543 / sqrt(0.527 * 1.200) = 0.68.
    # ------------------------------------------------------------------
    etalcl + etalq + etalvc + etalvp ~ c(
      0.272,
      -0.018, 0.527,
      -0.003, 0.543, 1.200,
      -0.047, 0.152, 0.211, 0.383
    ) # Table 2S Omega means: CL 0.272; Q,CL -0.018; Q 0.527; V1,CL -0.003; V1,Q 0.543; V1 1.200; V2,CL -0.047; V2,Q 0.152; V2,V1 0.211; V2 0.383

    # ------------------------------------------------------------------
    # Residual error -- log(Cobs) ~ t(log(C), sigmaC^-2, nu) (Methods
    # equations 4 and 7; WinBUGS dt(logCHat, tauC, nu), tauC = 1/sigmaC^2):
    # additive on the log scale, Student-t with scale sigmaC and nu
    # degrees of freedom. rxode2's lnorm(expSd) + dt(nu) uses the same
    # scale parameterisation.
    # ------------------------------------------------------------------
    expSd <- 0.275; label("Residual scale on the log concentration (Student-t)") # Table 2S: sigma C mean = 0.275 (90% HDI 0.233-0.317)
    nu <- 1.75; label("Degrees of freedom of the Student-t residual error")      # Table 2S: nu mean = 1.75 (90% HDI 1.42-2.14)
  })
  model({
    # 1. Covariates. Postmenstrual age in weeks as in the authors' WinBUGS
    # code: 4 * age_months + 40.
    pma_wk <- 4 * PNA + 40
    tm50_cl <- exp(ltm50_cl)
    hill_cl <- exp(lhill_cl)
    mat_cl <- pma_wk^hill_cl / (pma_wk^hill_cl + tm50_cl^hill_cl)

    # 2. Occasion-1 individual parameters (Methods equations 5 and 6).
    cl1 <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * mat_cl
    q1 <- exp(lq + etalq) * (WT / 70)^e_wt_q
    vc1 <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp1 <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp

    # 3. Occasion switch. The authors' ODE component
    # (DexTwoCptModelTCLRK45OCC) moves every parameter from its occasion-1
    # to its occasion-2 value with the logistic weight
    # 1 / (1 + exp(NNN * (t - TCL50))), NNN = 1 /h and TCL50 = 30 h hardcoded
    # (theta[i, 9] and theta[i, 10] in the WinBUGS code); the observation
    # volume Vc uses the same switch. Time is hours since the infusion
    # started: occasion-1 samples end at 24 h and the shortest infusion was
    # 45 h, so 30 h falls between the two occasions. The occasion-2 values
    # share the occasion-1 random effects (theta[i, 5:8] = f * theta[i, 1:4]).
    # The occasion-2 weight 1 - 1 / (1 + exp(t - 30)) is expit(t - 30).
    occ2 <- expit(t - 30)
    cl <- cl1 * (1 + e_occ2_cl * occ2)
    q <- q1 * (1 + e_occ2_q * occ2)
    vc <- vc1 * (1 + e_occ2_vc * occ2)
    vp <- vp1 * (1 + e_occ2_vp * occ2)

    # 4. Micro-constants (1/h) and ODEs in amounts (ug), as in the authors'
    # code; the infusion enters central.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation: ug / L = ng/mL.
    Cc <- central / vc
    Cc ~ lnorm(expSd) + dt(nu)
  })
}
