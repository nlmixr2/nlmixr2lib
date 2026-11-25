Kamp_2026_tacrolimus_combined <- function() {
  description <- paste(
    "Two-compartment population PK model for twice-daily immediate-release",
    "oral tacrolimus (whole blood) fitted jointly to adult simultaneous",
    "pancreas-kidney (SPK) and kidney-only (KT) transplant recipients during",
    "the first post-transplant year (Kamp 2026 model 3). First-order",
    "absorption with an absorption lag time, linear elimination, and oral",
    "bioavailability fixed at 0.23 with log-normal between-subject variability",
    "on it, so CL, Vc, Q and Vp are systemic (not apparent) values. Transplant",
    "type is the only covariate: SPK recipients have a 42.3% lower absorption",
    "rate constant Ka than KT recipients. Correlated log-normal IIV on Ka, Vc",
    "and CL; proportional residual error.",
    sep = " "
  )

  reference <- paste(
    "Kamp J, Gelinck A, Lyousoufi M, Meziyerh S, Spijker HS, Huurman VAL,",
    "van der Boog PJM, de Vries APJ, Moes DJAR. Differential Oral Absorption",
    "Pharmacokinetics of Tacrolimus in Simultaneous Kidney-Pancreas Transplant",
    "Compared to Kidney Transplant Recipients. Clin Pharmacokinet.",
    "2026;65(10):1517-1525. doi:10.1007/s40262-026-01689-3. PMCID PMC13623778.",
    "Final estimates from Supplemental Table S3 (combined model); structure,",
    "covariate equation and the OMEGA covariances from the NONMEM control",
    "stream in Supplementary Data S3.",
    sep = " "
  )

  vignette <- "Kamp_2026_tacrolimus"

  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    TX_KIDNEY_PANCREAS = list(
      description = paste(
        "Transplant type indicator: 1 = simultaneous pancreas-kidney (SPK)",
        "transplant recipient, 0 = kidney-only transplant recipient."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (kidney-only transplant)",
      notes = paste(
        "Kamp 2026 Methods eq. 2: KA = TVKA * (1 + theta_SPK); Supplementary",
        "Data S3 codes it as IF (NPTX.EQ.1) TVKA = TVKA * (1 + THETA(9)).",
        "Time-fixed per subject."
      ),
      source_name = "NPTX"
    )
  )

  covariatesDataExcluded <- list(
    POD = list(
      description = "Time after transplantation (post-operative day)",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as a three-level categorical effect on Ka and on the lag",
        "time in the separate SPK and KT models (Supplemental Tables S1 and",
        "S2); not retained in either, and not carried into the combined model."
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
    n_subjects = 179L,
    n_studies = 1L,
    n_observations = 2348L,
    age_range = "19-76 years",
    age_median = "45 years (SPK) and 44 years (KT)",
    weight_range = "40-104 kg",
    weight_median = "73 kg (SPK) and 68 kg (KT)",
    sex_female_pct = 39.1,
    disease_state = paste(
      "85 simultaneous pancreas-kidney and 94 kidney-only transplant",
      "recipients (80 with type 1 diabetes: 76 SPK, 4 KT) on tacrolimus /",
      "mycophenolate / prednisolone triple therapy, first year after",
      "transplantation."
    ),
    dose_range = paste(
      "Twice-daily immediate-release oral tacrolimus; protocol start 5 mg BID,",
      "individual doses 0.5-15 mg (median 4 mg)."
    ),
    regions = "Netherlands (Leiden University Medical Center)",
    notes = paste(
      "Retrospective single-centre cohort transplanted 2005-2020 (Kamp 2026",
      "Table 1): 1275 SPK and 1073 KT samples. Each sampling occasion was",
      "entered as a separate ID (no explicit IOV)."
    )
  )

  ini({
    # Structural parameters: Supplemental Table S3. The more precise lag time
    # is the $THETA final value of Supplementary Data S3, which agrees with
    # Table S3 at its printed rounding.
    lka <- log(0.360); label("Absorption rate constant Ka in kidney-only transplant recipients (1/h)") # Table S3 Ka = 0.360 1/h (RSE 10%; bootstrap median 0.361); S3 THETA(2) prints 0.336
    lcl <- log(5.8); label("Clearance CL (L/h)") # Table S3 CL = 5.8 L/h (RSE 3%); S3 THETA(4) = 5.8
    lvc <- log(14.0); label("Central volume of distribution Vc (L)") # Table S3 Vc = 14.0 L (RSE 17%); S3 THETA(3) = 14
    lq <- log(16.0); label("Intercompartmental clearance Q (L/h)") # Table S3 Q = 16.0 L/h (RSE 8%); S3 THETA(6) = 16
    lvp <- log(154); label("Peripheral volume of distribution Vp (L)") # Table S3 Vp = 154 L (RSE 15%); S3 THETA(5) = 154
    ltlag <- log(0.546); label("Absorption lag time (h)") # Table S3 lag time = 0.55 h (RSE 4%); S3 THETA(7) = 0.546
    lfdepot <- fixed(log(0.23)); label("Oral bioavailability F (fraction)") # Table S3 oral bioavailability = 23% (Fixed); S3 THETA(1) = 0.23 FIX

    # Covariate effect: KA = TVKA * (1 + theta_SPK) for SPK recipients.
    e_tx_kidney_pancreas_ka <- -0.423; label("Fractional change in Ka for simultaneous pancreas-kidney recipients (unitless)") # Table S3 'Ka - SPK (% reduction when SPK)' = 42.3% (RSE 10%); S3 THETA(9) = -0.423

    # IIV: S3 $OMEGA (log-normal variances). CV = sqrt(exp(omega^2) - 1)
    # reproduces Table S3: F 51.0%, Ka 76.4%, Vc 122.2%, CL 20.6%.
    etalfdepot ~ 0.231 # S3 $OMEGA IIV F1 = 0.231; Table S3 IIV F1 = 51.0% CV
    etalka + etalvc + etalcl ~ c(
      0.46,
      0.163, 0.914,
      0.0738, 0.124, 0.0414
    ) # S3 $OMEGA BLOCK(3) KA, V2, CL; Table S3 IIV Ka 76.4%, Vc 122.2%, CL 20.6% CV

    # Residual error: S3 W = THETA(8) * IPRED with $SIGMA 1 FIX, so THETA(8)
    # is the proportional SD.
    propSd <- 0.156; label("Proportional residual error (fraction)") # Table S3 proportional error = 15.6%; S3 THETA(8) = 0.156
  })

  model({
    ka <- exp(lka + etalka) * (1 + e_tx_kidney_pancreas_ka * TX_KIDNEY_PANCREAS)
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
