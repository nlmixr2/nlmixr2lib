Maganda_2015_lumefantrine <- function() {
  description <- paste(
    "Two-compartment population PK model for oral lumefantrine in",
    "HIV-malaria co-infected Tanzanian adults treated with",
    "artemether-lumefantrine (six doses over 3 days) while on",
    "nevirapine-based or efavirenz-based antiretroviral therapy, or not yet",
    "on antiretroviral therapy (Maganda 2015). First-order absorption with",
    "an absorption lag time; body weight enters allometrically on all",
    "apparent clearances (exponent 0.75) and volumes (exponent 1) centred",
    "at 70 kg. Concomitant efavirenz and nevirapine scale the relative",
    "bioavailability (0.42- and 1.32-fold). Between-subject variability on",
    "KA, Vc/F, Q/F, Vp/F and F, between-occasion variability on CL/F, and",
    "a combined additive and proportional residual error."
  )
  reference <- paste(
    "Maganda BA, Ngaimisi E, Kamuhabwa AAR, Aklillu E, Minzi OMS (2015).",
    "The influence of nevirapine and efavirenz-based anti-retroviral",
    "therapy on the pharmacokinetics of lumefantrine and anti-malarial dose",
    "recommendation in HIV-malaria co-treatment. Malaria Journal 14:179.",
    "doi:10.1186/s12936-015-0695-2"
  )
  vignette <- "Maganda_2015_lumefantrine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "lumefantrine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lumefantrine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "lumefantrine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of CL/F, Q/F (exponent 3/4) and Vc/F, Vp/F",
        "(exponent 1) to a 70 kg reference; KA and the lag time are not",
        "weight-scaled (Maganda 2015 Methods, 'Pharmacokinetic analysis').",
        "Cohort median weights were 55-58 kg by arm (Table 1)."
      ),
      source_name = "weight"
    ),
    CONMED_EFV = list(
      description = "Concomitant efavirenz-based antiretroviral therapy (600 mg once at night for at least two months)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not on antiretroviral therapy; the control arm)",
      notes = paste(
        "Relative bioavailability multiplied by 0.42 when CONMED_EFV = 1",
        "(Maganda 2015 Table 2 'F1 EFV'). The full-model clearance effect",
        "(CLARM2, Figure 2) was not retained. Mutually exclusive with",
        "CONMED_NVP in the source design."
      ),
      source_name = "ARM (EFV-arm)"
    ),
    CONMED_NVP = list(
      description = "Concomitant nevirapine-based antiretroviral therapy (200 mg twice daily for at least two months)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not on antiretroviral therapy; the control arm)",
      notes = paste(
        "Relative bioavailability multiplied by 1.32 when CONMED_NVP = 1",
        "(Maganda 2015 Table 2 'F1 NEV'). The full-model clearance effect",
        "(CLARM1, Figure 2) was not retained. Mutually exclusive with",
        "CONMED_EFV in the source design."
      ),
      source_name = "ARM (NVP-arm)"
    ),
    OCC = list(
      description = "Dose occasion index, 1 to 6 across the six-dose artemether-lumefantrine regimen",
      units = "(count)",
      type = "categorical",
      reference_category = "n/a (selects the between-occasion random effect on CL/F)",
      notes = paste(
        "Maganda 2015 estimates between-occasion variability on CL/F",
        "(Table 2) but does not state how an occasion was defined. The",
        "maintainers use one occasion per artemether-lumefantrine dose, the",
        "convention of the related Lohy Das 2018 lumefantrine model: doses",
        "at 0, 8, 24, 36, 48 and 60 h give OCC = 1 for 0 <= t < 8 h, 2 for",
        "8 <= t < 24 h, 3 for 24 <= t < 36 h, 4 for 36 <= t < 48 h, 5 for",
        "48 <= t < 60 h and 6 for t >= 60 h. For regimens with more than six",
        "doses, hold OCC at 6 after the sixth dose."
      ),
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 269L,
    n_studies = 1L,
    n_observations = 1514L,
    age_range = "19-67 years (adults >= 18 years; arm medians 38-43 years)",
    weight_range = "36-92 kg (Table 1 bracketed values, labelled IQR but spanning the range; arm medians 55-58 kg; inclusion required >= 35 kg)",
    sex_female_pct = 68.8,
    race_ethnicity = "Black African (Tanzanian)",
    disease_state = paste(
      "HIV-1 infection with microscopically confirmed uncomplicated",
      "Plasmodium falciparum malaria"
    ),
    dose_range = paste(
      "Artemether-lumefantrine (Coartem, 80 mg artemether / 480 mg",
      "lumefantrine per dose) at 0, 8, 24, 36, 48 and 60 h, each dose",
      "with 250 mL milk (3.5% fat)"
    ),
    regions = "Tanzania (Bagamoyo District Hospital HIV clinic)",
    co_medication = paste(
      "Three arms: nevirapine-based ART (n = 128), efavirenz-based ART",
      "(n = 66), and ART-naive control (n = 75)"
    ),
    notes = paste(
      "Demographics from Maganda 2015 Table 1. Female percentage is the",
      "subject-weighted mean of the per-arm values (65%, 79.5%, 52.3%).",
      "Sparse sampling: pre-dose 5 trough, up to three samples between 49",
      "and 120 h, and day 7 and day 14 samples; 17 samples (1.1%) below the",
      "50 ng/mL LLOQ were excluded."
    )
  )

  ini({
    # Structural parameters, Maganda 2015 Table 2 'Estimates (95% CI)'.
    # Apparent (relative to F = 1) values for a 70 kg patient.
    lka <- log(0.032)
    label("Absorption rate constant KA (1/h)") # Table 2 'Ka (hr-1)' = 0.032 (0.029-0.034)
    lvc <- log(25.6)
    label("Apparent central volume of distribution Vc/F (L)") # Table 2 'Vc/F (liters)' = 25.6 (16.21-34.99)
    lcl <- log(4.54)
    label("Apparent clearance CL/F (L/h)") # Table 2 'CL/F (liters/hr)' = 4.54 (3.913-5.167)
    lq <- log(1.23)
    label("Apparent inter-compartmental clearance Q/F (L/h)") # Table 2 'Q/F (liters/hr)' = 1.23 (0.99-1.47)
    lvp <- log(203)
    label("Apparent peripheral volume of distribution Vp/F (L)") # Table 2 'Vp/F (liters)' = 203 (167.13-238.87)
    ltlag <- log(1.45)
    label("Absorption lag time (h)") # Table 2 't lag (hr)' = 1.45 (1.25-1.65)
    lfdepot <- fixed(log(1))
    label("Relative bioavailability F, control arm (unitless)") # Table 2 'F1' = 1, footnote '*' fixed to this value

    # Allometric exponents, Methods 'Pharmacokinetic analysis': 'n is
    # exponent of value 3/4 for clearances and 1 for volumes'. Fixed, not
    # estimated.
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F and Q/F (unitless)") # Methods: exponent 3/4 for clearances
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on Vc/F and Vp/F (unitless)") # Methods: exponent 1 for volumes

    # Concomitant antiretroviral therapy on relative bioavailability. Table 2
    # reports each as the relative bioavailability of the arm (control = 1),
    # so it enters as a multiplicative factor raised to the 0/1 indicator.
    e_efv_fdepot <- 0.42
    label("Relative bioavailability on efavirenz-based ART vs control (ratio)") # Table 2 'F1 EFV' = 0.42 (0.34-0.5)
    e_nvp_fdepot <- 1.32
    label("Relative bioavailability on nevirapine-based ART vs control (ratio)") # Table 2 'F1 NEV' = 1.32 (1.08-1.52)

    # Between-subject variability. Table 2 reports BSV as a percentage and
    # Methods state 'Exponential models were assumed for all the random
    # effects', so the log-scale variance is omega^2 = log(CV^2 + 1).
    # CL/F and the lag time have BSV fixed to 0 (Table 2 '0*'), so they
    # carry no between-subject eta.
    etalka ~ 0.0807501 # Table 2 'Ka' BSV 29% (RSE 18%, shrinkage 30%); log(0.29^2 + 1)
    etalvc ~ 0.51426 # Table 2 'Vc/F' BSV 82% (RSE 24%, shrinkage 51%); log(0.82^2 + 1)
    etalq ~ 0.0703653 # Table 2 'Q/F' BSV 27% (RSE 24%, shrinkage 23%); log(0.27^2 + 1)
    etalvp ~ 0.141586 # Table 2 'Vp/F' BSV 39% (RSE 22%, shrinkage 37%); log(0.39^2 + 1)
    etalfdepot ~ 0.199588 # Table 2 'F1' BSV 47% (RSE 11%, shrinkage 5%); log(0.47^2 + 1)

    # Between-occasion variability on CL/F, one eta per dose occasion with
    # a shared variance (the first slot carries the estimate, the rest are
    # fixed to it).
    etaiov_cl_1 ~ 0.0354637 # Table 2 'CL/F' BOV 19% (RSE 22%, shrinkage 44%); log(0.19^2 + 1)
    etaiov_cl_2 ~ fixed(0.0354637)
    etaiov_cl_3 ~ fixed(0.0354637)
    etaiov_cl_4 ~ fixed(0.0354637)
    etaiov_cl_5 ~ fixed(0.0354637)
    etaiov_cl_6 ~ fixed(0.0354637)

    # Combined additive and proportional residual error (Methods: 'A
    # combined (additive and proportional) error model was used'), both
    # read as standard deviations.
    addSd <- 26.30
    label("Additive residual error SD (ng/mL)") # Table 2 'ADD' = 26.30 (4.15-48.44)
    propSd <- 0.083
    label("Proportional residual error SD (fraction)") # Table 2 'PROP' = 0.083 (0.06-0.10)
  })

  model({
    # Dose-occasion indicators (see covariateData$OCC for the time mapping).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6

    # Individual parameters; allometric scaling to 70 kg on clearances and
    # volumes only (Methods: 'All the parameter values, except Ka and tlag
    # were standardized for total body weight').
    ka <- exp(lka + etalka)
    cl <- exp(lcl + iov_cl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc
    tlag <- exp(ltlag)

    # Relative bioavailability with the antiretroviral-arm factors and
    # exponential between-subject variability.
    fdepot <- exp(lfdepot + etalfdepot) *
      e_efv_fdepot^CONMED_EFV * e_nvp_fdepot^CONMED_NVP

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # Dose in mg and volume in L give mg/L; x 1000 gives ng/mL (the unit of
    # the paper's concentrations, LLOQ 50 ng/mL).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
