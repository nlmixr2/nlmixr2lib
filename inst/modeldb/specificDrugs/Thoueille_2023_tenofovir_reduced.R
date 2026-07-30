Thoueille_2023_tenofovir_reduced <- function() {
  description <- paste(
    "One-compartment population PK model for tenofovir after oral tenofovir",
    "alafenamide in people living with HIV, with creatinine clearance on",
    "apparent clearance (reduced model, used by the authors for all",
    "chronic-kidney-disease simulations)."
  )
  reference <- paste(
    "Thoueille P, Alves Saldanha S, Desfontaine V, Kusejko K, Courlet P,",
    "Andre P, Cavassini M, Decosterd LA, Buclin T, Guidi M, and the Swiss HIV",
    "Cohort Study. Population pharmacokinetic modelling to characterize the",
    "effect of chronic kidney disease on tenofovir exposure after tenofovir",
    "alafenamide administration. J Antimicrob Chemother. 2023;78:1433-1443.",
    "doi:10.1093/jac/dkad103"
  )
  vignette <- "Thoueille_2023_tenofovir_ckd"
  # Doses were converted to nmol of tenofovir alafenamide and concentrations to
  # nmol/mL of tenofovir for the population analysis (Methods, "PopPK
  # analysis"). Complete, irreversible 1:1 molar conversion of tenofovir
  # alafenamide into tenofovir is assumed (Methods, "Structural model"), so a
  # dose expressed in nmol of tenofovir alafenamide enters the model as the
  # same number of nmol of tenofovir.
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/mL")

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance estimated with the Cockcroft-Gault equation.",
        "NOTE: raw Cockcroft-Gault in mL/min, NOT normalized to 1.73 m^2 body",
        "surface area, so the values are not interchangeable with the",
        "BSA-normalized eGFR variants that also map to this canonical."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference value CL_CR-ref = 100 mL/min (Table S2 footnote). Enters as",
        "a linear centered fractional effect on CL/F. Study range 33-203",
        "mL/min, median 93 mL/min (Table 1). The authors also screened CKD-EPI",
        "eGFR, which was significant but fit less well than Cockcroft-Gault",
        "CLCR (Results: dOFV -242 for eGFR vs -285 for CLCR); only CLCR was",
        "retained. Precedent for the raw, non-BSA-normalized Cockcroft-Gault",
        "form under this canonical: Delattre_2010_amikacin.R."
      ),
      source_name = "CLCR"
    ),
    CONMED_COBICISTAT = list(
      description = paste(
        "Concomitant cobicistat coadministration indicator: 1 = tenofovir",
        "alafenamide 10 mg once daily given with cobicistat, 0 = tenofovir",
        "alafenamide 25 mg once daily without cobicistat."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (tenofovir alafenamide 25 mg once daily, no cobicistat; F fixed to 1)",
      notes = paste(
        "Cobicistat was forced into the model as a covariate on relative",
        "bioavailability F to absorb the 25 mg vs 10 mg dose difference",
        "(Methods, 'Structural model'), so in this dataset the indicator is",
        "perfectly confounded with dose level: every cobicistat-treated",
        "subject received 10 mg and every other subject 25 mg. It therefore",
        "encodes the combined dose-normalization and boosting effect, not a",
        "pure drug-drug-interaction effect. 277 of 486 subjects (43%) received",
        "cobicistat (Table 1). Cobicistat was deliberately excluded from the",
        "P-glycoprotein-inhibitor covariate analysis (Methods, 'Covariate",
        "model')."
      ),
      source_name = "Cobicistat"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 486,
    n_observations = 793,
    n_studies = 1,
    age_range = "19-79 years",
    age_median = "51 years",
    weight_range = "42-142 kg",
    weight_median = "74 kg",
    height_median = "173 cm",
    bmi_median = "24 kg/m^2",
    sex_female_pct = 30,
    race_ethnicity = c(White = 71, Black = 23, `Hispanic American` = 4, Asian = 2),
    disease_state = "HIV-1 infection (people living with HIV) on tenofovir alafenamide-based antiretroviral therapy",
    renal_function = "serum creatinine median 85 (range 42-256) umol/L; Cockcroft-Gault CLCR median 93 (range 33-203) mL/min; CKD-EPI eGFR median 87 (range 23-153) mL/min/1.73 m^2",
    dose_range = "tenofovir alafenamide 25 mg once daily, or 10 mg once daily with cobicistat",
    regions = "Switzerland",
    notes = paste(
      "Baseline demographics from Thoueille 2023 Table 1, model-building",
      "dataset column (n = 486 subjects, 793 tenofovir concentrations).",
      "Subjects were enrolled in Swiss HIV Cohort Study project #815 or",
      "followed in the Lausanne routine therapeutic-drug-monitoring programme",
      "between January 2017 and January 2021. Mostly sparse sampling (1-2",
      "samples per subject) plus four subjects with detailed 0-24 h profiles.",
      "Steady state was assumed for all individuals. An independent external",
      "validation dataset of 83 subjects / 84 observations is described in the",
      "same table but was not used for estimation."
    )
  )

  ini({
    # Structural parameters -- Thoueille 2023 Table S2 ("Reduced model" column).
    # ka was fixed, not estimated: the small number of samples collected soon
    # after dose intake could not support it, and tenofovir shows flip-flop
    # kinetics after tenofovir alafenamide dosing (Methods, "Structural
    # model"; Table S2 reports "2 FIX").
    lka <- fixed(log(2)); label("First-order absorption rate constant ka (1/h)")                          # Table S2 (ka = 2 h^-1 FIX)
    lcl <- log(43.2); label("Apparent clearance CL/F of tenofovir at CLCR = 100 mL/min (L/h)")             # Table S2 (CL_TFV = 43.2 L/h)
    lvc <- log(2330); label("Apparent central volume of distribution V/F of tenofovir (L)")                # Table S2 (V_TFV = 2330 L)
    # F is a relative bioavailability anchored at 1 for the 25 mg
    # tenofovir alafenamide regimen without cobicistat; Table S2 reports
    # "1 FIX".
    lfdepot <- fixed(log(1)); label("Relative bioavailability F, tenofovir alafenamide 25 mg without cobicistat (unitless)")  # Table S2 (F = 1 FIX)

    # Covariate effects -- Table S2 "Reduced model:" equations.
    #   TVCL_TFV = CL_TFV * (1 + theta_CLCR * (CLCR - CLCR_ref) / CLCR_ref),
    #     with CLCR_ref = 100 mL/min
    #   TVF      = 1 + theta_Cobicistat
    e_crcl_cl <- 0.827; label("Fractional CLCR effect on CL/F per unit relative deviation from 100 mL/min")  # Table S2 (theta_CLCR = 0.827)
    e_cobi_fdepot <- 1.11; label("Fractional increase in F with cobicistat coadministration")                # Table S2 (theta_Cobicistat = 1.11)

    # BSV was supported only on F; adding BSV on ka, CL or V gave dOFV > -1
    # (Results, "Structural, statistical and covariate models"). Table S2
    # reports BSV on F as 22.5% CV. Parameters were assumed log-normally
    # distributed (Methods, "Structural model"), so the variance is
    # omega^2 = log(1 + CV^2) = log(1 + 0.225^2) = 0.0493852.
    etalfdepot ~ 0.0493852                                                                                   # Table S2 (BSV on F = 22.5% CV)

    # An additive error model best described tenofovir residual variability
    # (Results, "Structural, statistical and covariate models").
    addSd <- 0.0101; label("Additive residual error (nmol/mL)")                                               # Table S2 (sigma_add = 0.0101 nmol/mL)
  })

  model({
    ka <- exp(lka)
    # Linear centered CLCR effect on apparent clearance; reference 100 mL/min.
    cl <- exp(lcl) * (1 + e_crcl_cl * (CRCL - 100) / 100)
    vc <- exp(lvc)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Relative bioavailability carries the only between-subject random effect
    # and the cobicistat / dose-normalization effect.
    f(depot) <- exp(lfdepot + etalfdepot) * (1 + e_cobi_fdepot * CONMED_COBICISTAT)

    # `central` is in nmol and `vc` in L, so `central / vc` is nmol/L; the
    # extra factor of 1000 converts to the nmol/mL scale the authors used for
    # the concentrations and for sigma_add.
    Cc <- central / vc / 1000
    Cc ~ add(addSd)
  })
}
