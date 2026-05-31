Karunajeewa_2009_sulfadoxine <- function() {
  description <- paste(
    "Population PK model for sulfadoxine (SDOX) and its primary",
    "N-acetylsulfadoxine (NASDOX) metabolite in 60 Papua New Guinean",
    "women (30 pregnant, second or third trimester; 30 age-matched",
    "nonpregnant controls) given a single oral 1,500 mg sulfadoxine /",
    "75 mg pyrimethamine dose for intermittent presumptive treatment",
    "of malaria in pregnancy (Karunajeewa 2009). SDOX is described by",
    "first-order absorption (no lag) into a 2-compartment disposition",
    "with separate non-metabolic clearance CL/F (renal excretion) and",
    "metabolic formation clearance CLM/F that drains SDOX into a",
    "1-compartment NASDOX disposition. NASDOX elimination clearance",
    "is fixed at 10 times the structural SDOX non-metabolic CL/F",
    "(rapid formation-rate-limited renal excretion of the metabolite,",
    "Bell 1985). Allometric scaling is applied to all apparent",
    "volumes (exponent 1) and all apparent clearances (exponent",
    "0.75) at reference WT = 70 kg. Pregnancy is the only retained",
    "covariate, entering as an additive term on the structural SDOX",
    "non-metabolic CL/F (+0.0181 L/h/70 kg). The companion model for",
    "the co-administered pyrimethamine is shipped as",
    "'Karunajeewa_2009_pyrimethamine' (separate NONMEM dataset, fit",
    "independently in the source publication)."
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
    concentration = "mg/L for sulfadoxine; mg/L for N-acetylsulfadoxine"
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
        "Exponent 1.0 on all apparent volumes (Vc, Vp for SDOX;",
        "Vc for NASDOX); exponent 0.75 on all apparent clearances",
        "(CL/F, CLM/F, Q/F for SDOX; CL/F NASDOX). Baseline body",
        "weight was 54.0 +/- 6.4 kg (pregnant) and 51.8 +/- 5.5 kg",
        "(nonpregnant) per Table 1."
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
        "only via an additive term on SDOX non-metabolic clearance",
        "(CL/F): CL_typical = exp(lcl) + e_preg_cl * PREG. None of",
        "age, gestational age, hemoglobin, parasitemia, or blood",
        "glucose was retained in the final SDOX/NASDOX model",
        "(Results paragraph 'Pharmacokinetics of SDOX and NASDOX')."
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
      "then 7, 10, 14, 28, 42 days post-dose. SDOX, NASDOX, and",
      "pyrimethamine assayed by HPLC-UV (LOQ 0.1 mg/L, 0.02 mg/L,",
      "and 2.5 ug/L respectively). NONMEM v6.2.0 with FOCE-INTER",
      "estimation. Bootstrap n = 1000."
    )
  )

  ini({
    # ============================================================
    # SDOX structural parameters
    # ----- Karunajeewa 2009 Table 2 'Final covariate model' column
    # All apparent volumes and clearances are reported per 70 kg
    # (Methods 'Population pharmacokinetic analysis').
    # ============================================================
    lka <- log(0.769)
    label("SDOX first-order absorption rate, ka (1/h)")                                          # Table 2: ka = 0.769 /h
    lcl <- log(0.0379)
    label("SDOX non-metabolic apparent clearance CL/F in nonpregnant women at WT = 70 kg (L/h)") # Table 2: CL/F SDOX = 0.0379 L/h/70kg
    lvc <- log(15.8)
    label("SDOX apparent central volume Vc/F at WT = 70 kg (L)")                                 # Table 2: Vc/F SDOX = 15.8 L/70kg
    lq <- log(0.0052)
    label("SDOX apparent inter-compartmental clearance Q/F at WT = 70 kg (L/h)")                 # Table 2: CL_Q/F = 0.0052 L/h/70kg
    lvp <- log(1.11)
    label("SDOX apparent peripheral volume Vp/F at WT = 70 kg (L)")                              # Table 2: Vp/F SDOX = 1.11 L/70kg
    lclm <- log(0.0227)
    label("SDOX -> NASDOX formation (metabolic) apparent clearance CLM/F at WT = 70 kg (L/h)")   # Table 2: CL_M/F = 0.0227 L/h/70kg
    lvc_nasdox <- log(3.69)
    label("NASDOX apparent volume V/F at WT = 70 kg (L)")                                        # Table 2: V/F NASDOX = 3.69 L/70kg
    lfdepot <- fixed(log(1))
    label("SDOX relative bioavailability F (unitless, FIXED at 1)")                              # Methods 'Population pharmacokinetic analysis': all V and CL parameters expressed relative to F

    # ============================================================
    # Pregnancy effect on SDOX non-metabolic CL/F (additive)
    # Final model: CL/F_typical = exp(lcl) + e_preg_cl * PREG.
    # The additive parameterisation matches the published Table 2
    # row "Pregnancy on CL/F SDOX (liters/h/70 kg)" verbatim.
    # ============================================================
    e_preg_cl <- 0.0181
    label("Pregnancy additive effect on SDOX non-metabolic CL/F at WT = 70 kg (L/h)")            # Table 2: Pregnancy on CL/F SDOX = 0.0181 L/h/70kg

    # ============================================================
    # IIV (log-normal). omega^2 = log(1 + CV^2). CV% from Table 2
    # 'Final covariate model' column.
    # ============================================================
    etalcl  ~ log(1 + 0.281^2)
    # Table 2: BSV CL/F SDOX = 28.1% CV   ->  omega^2 = log(1 + 0.281^2)
    etalvc  ~ log(1 + 0.199^2)
    # Table 2: BSV V_C/F SDOX = 19.9% CV  ->  omega^2 = log(1 + 0.199^2)
    etalka  ~ log(1 + 0.996^2)
    # Table 2: BSV ka = 99.6% CV          ->  omega^2 = log(1 + 0.996^2)
    etalclm ~ log(1 + 0.273^2)
    # Table 2: BSV CL_M/F = 27.3% CV      ->  omega^2 = log(1 + 0.273^2)

    # ============================================================
    # Residual error.
    # SDOX: proportional only.
    # NASDOX: combined proportional + additive.
    # ============================================================
    propSd <- 0.241
    label("SDOX proportional residual SD (fraction)")              # Table 2: proportional error SDOX = 24.1%
    propSd_nasdox <- 0.160
    label("NASDOX proportional residual SD (fraction)")            # Table 2: proportional error NASDOX = 16.0%
    addSd_nasdox <- 0.2
    label("NASDOX additive residual SD (mg/L)")                    # Table 2: additive error NASDOX = 0.2 mg/L
  })

  model({
    # ------------------------------------------------------------
    # Allometric scaling at reference WT = 70 kg (Methods
    # 'Population pharmacokinetic analysis'). Exponents fixed at
    # the Anderson & Holford 2008 canonical values: 1.0 on
    # volumes, 0.75 on clearances and inter-compartmental flow.
    # ------------------------------------------------------------
    wt_ratio <- WT / 70

    # ------------------------------------------------------------
    # SDOX typical CL/F: additive pregnancy effect on the structural
    # non-metabolic clearance. PREG = 1 raises CL/F by 0.0181
    # L/h/70 kg in line with the additive Table 2 parameterisation.
    # The eta multiplies the pregnancy-conditional typical value
    # (NONMEM-style TVCL = THETA(1) + THETA(2) * PREG ; CL = TVCL *
    # EXP(ETA)).
    # ------------------------------------------------------------
    cl_typ <- exp(lcl) + e_preg_cl * PREG
    cl <- cl_typ * exp(etalcl) * wt_ratio^0.75
    vc <- exp(lvc + etalvc) * wt_ratio
    q  <- exp(lq)           * wt_ratio^0.75
    vp <- exp(lvp)          * wt_ratio
    ka <- exp(lka + etalka)
    clm <- exp(lclm + etalclm) * wt_ratio^0.75
    vc_nasdox <- exp(lvc_nasdox) * wt_ratio

    # NASDOX elimination clearance fixed at 10 x SDOX non-metabolic
    # CL/F per Figure 2 caption: "CL/F NASDOX (fixed at 10 x CL/F
    # SDOX)" (Bell 1985 i.v. comparison showed renal CL of NASDOX
    # is approximately 10-fold that of SDOX). The 10x relationship
    # is structural and carries the SDOX-CL individual eta through
    # for mass-balance consistency.
    cl_nasdox <- 10 * cl

    # ------------------------------------------------------------
    # SDOX 2-compartment disposition with first-order absorption.
    # Total SDOX elimination from central = renal CL + metabolic
    # formation of NASDOX, i.e., (cl + clm) * Cc.
    # ------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl + clm) * central / vc -
                          q * central / vc +
                          q * peripheral1 / vp
    d/dt(peripheral1) <-  q * central / vc -
                          q * peripheral1 / vp

    # NASDOX formation from SDOX central; 1-compartment metabolite
    # disposition with first-order elimination (Figure 2).
    d/dt(central_nasdox) <-  clm * central / vc -
                             cl_nasdox * central_nasdox / vc_nasdox

    # SDOX bioavailability anchored at 1 (FIXED).
    f(depot) <- exp(lfdepot)

    # ------------------------------------------------------------
    # Observations.
    # SDOX: plasma concentration in central / Vc (mg/L).
    # NASDOX: plasma concentration in central_nasdox / Vc_nasdox
    #   (mg/L).
    # ------------------------------------------------------------
    Cc        <- central        / vc
    Cc_nasdox <- central_nasdox / vc_nasdox

    Cc        ~ prop(propSd)
    Cc_nasdox ~ add(addSd_nasdox) + prop(propSd_nasdox)
  })
}
