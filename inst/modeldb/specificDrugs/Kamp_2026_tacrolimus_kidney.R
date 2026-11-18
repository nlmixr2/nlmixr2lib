Kamp_2026_tacrolimus_kidney <- function() {
  description <- paste(
    "Two-compartment population PK model for twice-daily immediate-release",
    "oral tacrolimus (whole blood) in adult kidney-only transplant (KT)",
    "recipients during the first post-transplant year (Kamp 2026 model 2).",
    "First-order absorption with an absorption lag time, linear",
    "elimination, and oral bioavailability fixed at 0.23 with log-normal",
    "between-subject variability on it, so CL, Vc, Q and Vp are systemic (not",
    "apparent) values. Correlated log-normal IIV on Ka, Vc and CL;",
    "proportional residual error. No covariates were retained.",
    sep = " "
  )

  reference <- paste(
    "Kamp J, Gelinck A, Lyousoufi M, Meziyerh S, Spijker HS, Huurman VAL,",
    "van der Boog PJM, de Vries APJ, Moes DJAR. Differential Oral Absorption",
    "Pharmacokinetics of Tacrolimus in Simultaneous Kidney-Pancreas Transplant",
    "Compared to Kidney Transplant Recipients. Clin Pharmacokinet.",
    "2026;65(10):1517-1525. doi:10.1007/s40262-026-01689-3. PMCID PMC13623778.",
    "Final estimates from Table 2 (kidney model); structure and the OMEGA",
    "covariances from the NONMEM control stream in Supplementary Data S2.",
    sep = " "
  )

  vignette <- "Kamp_2026_tacrolimus"

  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list()

  covariatesDataExcluded <- list(
    POD = list(
      description = "Time after transplantation (post-operative day)",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as a three-level categorical effect on Ka (0-6 weeks",
        "reference, 6-26 weeks, > 26 weeks) and separately on the lag time;",
        "neither was retained (Supplemental Table S2 models 4 and 5:",
        "dOFV +10.09 and +12.75, i.e. worse fits)."
      ),
      source_name = "time after transplantation"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "tacrolimus",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "tacrolimus",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "tacrolimus",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 94L,
    n_studies = 1L,
    n_observations = 1073L,
    age_range = "19-76 years",
    age_median = "44 years",
    weight_range = "40-95 kg",
    weight_median = "68 kg",
    height_median = "175 cm (range 140-189)",
    sex_female_pct = 42.6,
    disease_state = paste(
      "Kidney-only transplant recipients (4 of 94 with type 1 diabetes; one",
      "patient older than 55 years) on tacrolimus / mycophenolate / prednisolone triple therapy,",
      "first year after transplantation."
    ),
    renal_function = "Serum creatinine median 123 umol/L (range 43-1101).",
    dose_range = paste(
      "Twice-daily immediate-release oral tacrolimus; protocol start 5 mg BID,",
      "individual doses 0.5-15 mg (median 4 mg)."
    ),
    regions = "Netherlands (Leiden University Medical Center)",
    notes = paste(
      "Retrospective single-centre cohort transplanted 2005-2020 (Kamp 2026",
      "Table 1). Sampling: C0 from about day 4, a 7-point AUC (0-6 h) at week 1",
      "and abbreviated 4-point AUCs (0-3 h) at 6 weeks, 6 months and 1 year;",
      "sampling moments 3-535 days after transplantation. Each sampling",
      "occasion was entered as a separate ID (no explicit IOV)."
    )
  )

  ini({
    # Structural parameters: Table 2 (kidney model); the more precise values
    # are the $THETA final values of Supplementary Data S2, which agree
    # with Table 2 at its printed rounding.
    lka <- log(0.336); label("Absorption rate constant Ka (1/h)") # Table 2 KA = 0.336 1/h (RSE 17%); S2 THETA(2) = 0.336
    lcl <- log(5.35); label("Clearance CL (L/h)") # Table 2 CL = 5.35 L/h (RSE 5%); S2 THETA(4) = 5.35
    lvc <- log(11.8); label("Central volume of distribution Vc (L)") # Table 2 Vc = 11.8 L (RSE 50%); S2 THETA(3) = 11.8
    lq <- log(14.9); label("Intercompartmental clearance Q (L/h)") # Table 2 Q = 14.9 L/h (RSE 13%); S2 THETA(6) = 14.9
    lvp <- log(177); label("Peripheral volume of distribution Vp (L)") # Table 2 Vp = 177 L (RSE 33%); S2 THETA(5) = 177
    ltlag <- log(0.533); label("Absorption lag time (h)") # Table 2 lag time = 0.53 h (RSE 7%); S2 THETA(7) = 0.533
    lfdepot <- fixed(log(0.23)); label("Oral bioavailability F (fraction)") # Table 2 oral bioavailability = 23% (Fixed); S2 THETA(1) = 0.23 FIX

    # IIV: S2 $OMEGA (log-normal variances). CV = sqrt(exp(omega^2) - 1)
    # reproduces Table 2: F 53.5%, Ka 67.8%, Vc 143.7%, CL 23.7%.
    etalfdepot ~ 0.252 # S2 $OMEGA IIV F1 = 0.252; Table 2 IIV F = 53.5% CV
    etalka + etalvc + etalcl ~ c(
      0.378,
      0.2, 1.12,
      0.0994, 0.133, 0.0547
    ) # S2 $OMEGA BLOCK(3) KA, V2, CL; Table 2 IIV Ka 67.8%, Vc 143.7%, CL 23.7% CV

    # Residual error: S2 W = THETA(8) * IPRED with $SIGMA 1 FIX, so THETA(8)
    # is the proportional SD.
    propSd <- 0.162; label("Proportional residual error (fraction)") # Table 2 proportional error = 16.2% (RSE 11%); S2 THETA(8) = 0.162
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    tlag <- exp(ltlag)
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # Dose in mg and volume in L give mg/L; x 1000 converts to ug/L (the
    # units of the measured whole-blood concentrations).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
