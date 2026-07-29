Karunajeewa_2009_pyrimethamine <- function() {
  description <- paste(
    "Population PK model for pyrimethamine (PYR) in 60 Papua New",
    "Guinean women (30 pregnant, second or third trimester; 30",
    "age-matched nonpregnant controls) given a single oral 1,500 mg",
    "sulfadoxine / 75 mg pyrimethamine dose for intermittent",
    "presumptive treatment of malaria in pregnancy (Karunajeewa",
    "2009). Two-compartment disposition with first-order absorption",
    "and no lag, fit as a separate NONMEM dataset from the parent",
    "SDOX/NASDOX dataset. Allometric scaling at reference WT = 70 kg",
    "is applied to all apparent volumes (exponent 1) and all",
    "apparent clearances (exponent 0.75). Pregnancy is the only",
    "retained covariate; it enters as additive terms on apparent",
    "CL/F (+0.439 L/h/70 kg), Vc/F (+76 L/70 kg) and Vp/F (+98",
    "L/70 kg). Between-subject variability on CL/F, Vc/F and Vp/F",
    "is correlated (3x3 block, correlations 0.797 / 0.756 / 0.731",
    "from Table 4); BSV on Q/F and ka is independent. The companion",
    "model for the co-administered sulfadoxine plus its NASDOX",
    "metabolite is shipped as 'Karunajeewa_2009_sulfadoxine'",
    "(separate NONMEM dataset, fit independently in the source",
    "publication)."
  )
  reference <- paste(
    "Karunajeewa HA, Salman S, Mueller I, Baiwog F, Gomorrai S, Law I,",
    "Page-Sharp M, Rogerson S, Siba P, Ilett KF, Davis TME.",
    "Pharmacokinetic properties of sulfadoxine-pyrimethamine in",
    "pregnant women. Antimicrob Agents Chemother.",
    "2009;53(10):4368-4376. doi:10.1128/AAC.00335-09.",
    sep = " "
  )
  vignette <- "Karunajeewa_2009_sulfadoxinePyrimethamine"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ug/L (= ng/mL) for pyrimethamine plasma"
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling applied at reference WT = 70 kg",
        "(Anderson & Holford 2008 convention, paper Methods).",
        "Exponent 1.0 on all apparent volumes (Vc/F, Vp/F);",
        "exponent 0.75 on all apparent clearances (CL/F, Q/F).",
        "Baseline body weight was 54.0 +/- 6.4 kg (pregnant) and",
        "51.8 +/- 5.5 kg (nonpregnant) per Table 1."
      ),
      source_name        = "WT"
    ),
    PREG = list(
      description        = "Pregnancy status indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = second- or third-trimester pregnant; 0 = age-matched",
        "nonpregnant control. Pregnancy enters the structural model",
        "via additive terms on CL/F, Vc/F, and Vp/F:",
        "CL_typical = exp(lcl) + e_preg_cl * PREG;",
        "Vc_typical = exp(lvc) + e_preg_vc * PREG;",
        "Vp_typical = exp(lvp) + e_preg_vp * PREG. The additive",
        "parameterisation matches the verbatim Table 4 published",
        "values (in L/h/70 kg or L/70 kg). None of age, gestational",
        "age, hemoglobin, parasitemia, or blood glucose was",
        "retained in the final PYR model (Results paragraph",
        "'Pharmacokinetics of pyrimethamine')."
      ),
      source_name        = "PREG"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 60,
    n_subjects_pregnant     = 30,
    n_subjects_nonpregnant  = 30,
    n_studies      = 1,
    age_range      = "Mean 26.0 +/- 5.9 years (pregnant); 25.5 +/- 8.9 years (nonpregnant) (Table 1).",
    weight_range   = "Mean 54.0 +/- 6.4 kg (pregnant); 51.8 +/- 5.5 kg (nonpregnant) (Table 1).",
    sex_female_pct = 100,
    race_ethnicity = "Melanesian (Papua New Guinean); the local population at the study site is described as almost exclusively Melanesian (Methods 'Study site and sample').",
    disease_state  = paste(
      "Asymptomatic first-time antenatal attendees and age-matched",
      "nonpregnant village controls. P. falciparum parasitemia at",
      "baseline in 43% (pregnant) and 23% (nonpregnant); P. vivax",
      "and P. malariae also present at lower rates (Table 1). No",
      "subject had axillary temperature > 37.5 degC at enrolment.",
      "Median gestational age 22 weeks [IQR 20-28] for the pregnant",
      "cohort."
    ),
    dose_range     = paste(
      "Single oral 1,500 mg sulfadoxine + 75 mg pyrimethamine",
      "(Fansidar, Roche) under direct supervision, with three daily",
      "doses of chloroquine 450 mg base co-administered per Papua",
      "New Guinea national IPTp guidelines (Methods 'Clinical",
      "procedures')."
    ),
    regions        = "Papua New Guinea (Alexishafen Health Centre, Madang Province, north coast).",
    notes          = paste(
      "Sampling: pre-dose and 1, 2, 4, 6, 12, 18, 24, 30, 48, 72 h",
      "then 7, 10, 14, 28, 42 days post-dose. PYR assayed by",
      "HPLC-UV (LOQ 2.5 ug/L). NONMEM v6.2.0 with FOCE-INTER",
      "estimation. Bootstrap n = 1000."
    )
  )

  ini({
    # ============================================================
    # PYR structural parameters
    # ----- Karunajeewa 2009 Table 4 'Final covariate model' column
    # All apparent volumes and clearances are reported per 70 kg.
    # ============================================================
    lka <- log(1.84)
    label("PYR first-order absorption rate, ka (1/h)")                                    # Table 4: ka = 1.84 /h
    lcl <- log(0.87)
    label("PYR apparent CL/F in nonpregnant women at WT = 70 kg (L/h)")                   # Table 4: CL/F = 0.87 L/h/70kg
    lvc <- log(146)
    label("PYR apparent central volume Vc/F in nonpregnant women at WT = 70 kg (L)")      # Table 4: Vc/F = 146 L/70kg
    lq <- log(0.51)
    label("PYR apparent inter-compartmental clearance Q/F at WT = 70 kg (L/h)")           # Table 4: Q/F = 0.51 L/h/70kg
    lvp <- log(76.8)
    label("PYR apparent peripheral volume Vp/F in nonpregnant women at WT = 70 kg (L)")   # Table 4: Vp/F = 76.8 L/70kg
    lfdepot <- fixed(log(1))
    label("PYR relative bioavailability F (unitless, FIXED at 1)")                        # Methods: all V and CL parameters expressed relative to F

    # ============================================================
    # Pregnancy effects (all additive) on CL/F, Vc/F, Vp/F.
    # Encoded as TVCL = exp(lcl) + e_preg_cl * PREG (etc.) to
    # preserve the verbatim published parameter values in their
    # original units.
    # ============================================================
    e_preg_cl <- 0.439
    label("Pregnancy additive effect on PYR CL/F at WT = 70 kg (L/h)")                    # Table 4: Pregnancy on CL/F = 0.439 L/h/70kg
    e_preg_vc <- 76
    label("Pregnancy additive effect on PYR Vc/F at WT = 70 kg (L)")                      # Table 4: Pregnancy on V_C/F = 76 L/70kg
    e_preg_vp <- 98
    label("Pregnancy additive effect on PYR Vp/F at WT = 70 kg (L)")                      # Table 4: Pregnancy on V_P/F = 98 L/70kg

    # ============================================================
    # Correlated 3x3 IIV block on CL/F, Vc/F, Vp/F.
    # omega^2 = log(1 + CV^2); off-diagonals = rho * sqrt(var_i *
    # var_j). CV% and correlations from Table 4 'Final covariate
    # model' column.
    # ============================================================
    etalcl + etalvc + etalvp ~ c(
      log(1 + 0.276^2),
      0.797 * sqrt(log(1 + 0.276^2) * log(1 + 0.246^2)), log(1 + 0.246^2),
      0.756 * sqrt(log(1 + 0.276^2) * log(1 + 0.425^2)),
      0.731 * sqrt(log(1 + 0.246^2) * log(1 + 0.425^2)), log(1 + 0.425^2)
    )
    # Table 4: BSV CL/F = 27.6% CV; BSV V_C/F = 24.6% CV; BSV V_P/F = 42.5% CV.
    # Correlations R(CL,Vc) = 0.797; R(CL,Vp) = 0.756; R(Vc,Vp) = 0.731.

    # Independent IIV on Q/F and ka.
    etalq  ~ log(1 + 0.575^2)
    # Table 4: BSV Q/F = 57.5% CV   ->  omega^2 = log(1 + 0.575^2)
    etalka ~ log(1 + 0.107^2)
    # Table 4: BSV ka = 10.7% CV   ->  omega^2 = log(1 + 0.107^2). The
    # bootstrap median for BSV ka (Table 4 'Bootstrap' column = 112
    # [67-187]) disagrees markedly with the Final-covariate-model
    # estimate of 10.7% CV; the published Final-covariate-model
    # value is reproduced verbatim here and the discrepancy is
    # documented in the vignette Errata.

    # ============================================================
    # Residual error (proportional only).
    # ============================================================
    propSd <- 0.192
    label("PYR proportional residual SD (fraction)")               # Table 4: proportional error PYR = 19.2% CV (column-header units "ug/liter" in Table 4 are a typo; the value is a unitless fractional CV per the SDOX/NASDOX Table 2 header)
  })

  model({
    # ------------------------------------------------------------
    # Allometric scaling at reference WT = 70 kg (Methods
    # 'Population pharmacokinetic analysis'). Exponents fixed at
    # the Anderson & Holford 2008 canonical values: 1.0 on
    # volumes, 0.75 on clearance and inter-compartmental flow.
    # ------------------------------------------------------------
    wt_ratio <- WT / 70

    # ------------------------------------------------------------
    # PYR typical-value structure: additive pregnancy effects on
    # CL/F, Vc/F, Vp/F. The eta multiplies the
    # pregnancy-conditional typical value (NONMEM-style TVCL =
    # THETA(1) + THETA(2) * PREG; CL = TVCL * EXP(ETA), and the
    # analogous form for Vc and Vp).
    # ------------------------------------------------------------
    cl_typ <- exp(lcl) + e_preg_cl * PREG
    vc_typ <- exp(lvc) + e_preg_vc * PREG
    vp_typ <- exp(lvp) + e_preg_vp * PREG

    cl <- cl_typ * exp(etalcl) * wt_ratio^0.75
    vc <- vc_typ * exp(etalvc) * wt_ratio
    vp <- vp_typ * exp(etalvp) * wt_ratio
    q  <- exp(lq + etalq)      * wt_ratio^0.75
    ka <- exp(lka + etalka)

    # ------------------------------------------------------------
    # PYR 2-compartment disposition with first-order absorption.
    # ------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          cl * central / vc -
                          q * central / vc +
                          q * peripheral1 / vp
    d/dt(peripheral1) <-  q * central / vc -
                          q * peripheral1 / vp

    # PYR bioavailability anchored at 1 (FIXED).
    f(depot) <- exp(lfdepot)

    # Concentration in ug/L (= ng/mL): central is in mg, vc in L,
    # so central/vc is in mg/L; multiply by 1000 to get ug/L.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
