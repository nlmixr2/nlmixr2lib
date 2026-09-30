Jeong_2020_tiropramide <- function() {
  description <- "One-compartment population PK model for oral tiropramide in healthy Korean adult men, with an absorption lag time followed by two sequential first-order absorption steps (depot -> transit1 -> central) and linear centred total-protein effects on CL/F and V/F (Jeong 2020)."
  reference <- "Jeong SH, Jang JH, Cho HY, Lee YB. Population Pharmacokinetic Analysis of Tiropramide in Healthy Korean Subjects. Pharmaceutics. 2020;12(4):374. doi:10.3390/pharmaceutics12040374"
  vignette <- "Jeong_2020_tiropramide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "The paper reports total protein in g/dL (median 7.6 g/dL, range 6.7-8.3; Table 1). The model converts the canonical g/L column inline (tpro_gdL = TPRO / 10) and applies the paper's linear centred form (1 + (TP - 7.6) * slope) on CL/F and V/F (Equation 1). The linear form becomes non-positive above TP = 8.55 g/dL for V/F and 8.84 g/dL for CL/F, so the model should not be used outside the observed 6.7-8.3 g/dL range.",
      source_name = "Totalproteins"
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Tested on Ka1 in the stepwise search (Table 5, model 4; dOFV -1.95) but not retained."
    ),
    ABCB1_C1236T_HET = list(
      description = "ABCB1 1236C>T genotype (the paper screened CC / CT / TT as an index variable)",
      units = "(binary)",
      type = "binary",
      notes = "Tested on Ka1 in the stepwise search (Table 5, model 6; dOFV -0.93) but not retained. The exact index coding was not reported."
    ),
    CYP2D6_STAR10_HET = list(
      description = "CYP2D6 *1/*10 heterozygous intermediate-metabolizer indicator",
      units = "(binary)",
      type = "binary",
      notes = "CYP2D6 (*1/*1, *1/*10, *10/*10) tested on Ka2 in the stepwise search (Table 5, model 7; dOFV -5.42 for 2 added parameters) but not retained."
    ),
    CYP2D6_STAR10_HOM = list(
      description = "CYP2D6 *10/*10 homozygous intermediate-metabolizer indicator",
      units = "(binary)",
      type = "binary",
      notes = "See CYP2D6_STAR10_HET. PEPT1 1287G>C (Table 5, model 5), OCT2 808G>T, ABCB1 2677G>T/A and 3435C>T were also screened without significant effect (Section 3.5)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tiropramide", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "tiropramide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tiropramide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    n_observations = 288L,
    age_range = "19-29 years",
    age_median = "23 years",
    weight_range = "52.5-82.8 kg",
    weight_median = "66.8 kg",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state = "healthy volunteers",
    dose_range = "100 mg single oral tablet (reference formulation of a two-way crossover bioequivalence study), fasted",
    regions = "Republic of Korea",
    total_protein_range = "6.7-8.3 g/dL (median 7.6)",
    notes = "Healthy Korean men from a randomized single-dose crossover bioequivalence study; only the reference-formulation period was analysed. Sampling 0-12 h. Demographics from Table 1; genotypes (ABCB1, CYP2D6, OCT2, PEPT1) in Table 2. Estimated in Phoenix NLME 8.1 (FOCE-ELS)."
  )

  ini({
    lktr <- log(3.187); label("First absorption-phase rate constant, depot to transit1 (paper Ka1, 1/h)") # Table 6 final model: tvKa1 = 3.187 1/h
    lka <- log(3.183); label("Second absorption-phase rate constant, transit1 to central (paper Ka2, 1/h)") # Table 6 final model: tvKa2 = 3.183 1/h
    lvc <- log(1889.250002); label("Apparent central volume V/F at total protein 7.6 g/dL (L)") # Table 6 final model: tvV/F = 1,889,250.002 mL
    lcl <- log(466.711101); label("Apparent clearance CL/F at total protein 7.6 g/dL (L/h)") # Table 6 final model: tvCL/F = 466,711.101 mL/h
    ltlag <- log(0.196); label("Absorption lag time on the dosing depot (h)") # Table 6 final model: tvTlag = 0.196 h

    e_tpro_cl <- -0.804; label("Linear slope of total protein (g/dL, centred at 7.6) on CL/F (per g/dL)") # Table 6 final model: dCl/FdTotalproteins = -0.804; Equation 1
    e_tpro_vc <- -1.049; label("Linear slope of total protein (g/dL, centred at 7.6) on V/F (per g/dL)") # Table 6 final model: dV/FdTotalproteins = -1.049; Equation 1

    etalktr ~ 0.557 # Table 6 final model: omega2 Ka1 = 0.557
    etalka ~ 0.556 # Table 6 final model: omega2 Ka2 = 0.556
    etalvc ~ 0.326 # Table 6 final model: omega2 V/F = 0.326
    etalcl ~ 0.160 # Table 6 final model: omega2 Cl/F = 0.160
    etaltlag ~ 0.107 # Table 6 final model: omega2 Tlag = 0.107

    expSd <- 0.357; label("Additive residual error on log-transformed concentration (SD, log scale)") # Table 6 final model: sigma = 0.357; Section 2.6 Cobs = Cpred * exp(eps)
  })

  model({
    # Total protein is supplied in g/L (canonical); the paper's covariate
    # equation is written in g/dL and centred at the 7.6 g/dL median.
    tpro_gdL <- TPRO / 10

    ktr <- exp(lktr + etalktr)
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc) * (1 + (tpro_gdL - 7.6) * e_tpro_vc)
    cl <- exp(lcl + etalcl) * (1 + (tpro_gdL - 7.6) * e_tpro_cl)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc

    # Figure 2: dose enters depot 1, which empties into depot 2 (transit1)
    # at Ka1 after the lag time; depot 2 empties into central at Ka2.
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ka * transit1
    d/dt(central) <- ka * transit1 - kel * central

    alag(depot) <- tlag

    # Dose in mg, volume in L -> mg/L; x 1000 for ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
