Faraj_2022_marzeptacogAlfa <- function() {
  description <- paste(
    "Two-compartment population PK model for marzeptacog alfa (activated)",
    "(MarzAA), an activated recombinant human factor VII variant, after",
    "intravenous and subcutaneous dosing in adult males with haemophilia A",
    "or B with or without inhibitors (Faraj 2022). Subcutaneous doses enter",
    "an absorption depot by a zero-order input of duration D1 and are then",
    "absorbed first-order into the central compartment; bioavailability is",
    "logit-transformed. Elimination is linear. Allometric body-weight",
    "scaling (reference 70 kg) uses fixed exponents of 0.75 on CL and Q and",
    "1 on Vc and Vp. The observed FVIIa clotting activity (ng/mL) is the",
    "sum of an estimated endogenous FVIIa baseline and the drug",
    "concentration. IIV on CL and Vc (correlated), baseline, F and ka;",
    "inter-occasion variability on D1, F and Vc; proportional residual",
    "error."
  )
  reference <- paste(
    "Faraj A, Knudsen T, Desai S, Neuman L, Blouse GE, Simonsson USH.",
    "Phase III dose selection of marzeptacog alfa (activated) informed by",
    "population pharmacokinetic modeling: A novel hemostatic drug.",
    "CPT Pharmacometrics Syst Pharmacol. 2022;11(12):1628-1637.",
    "doi:10.1002/psp4.12872"
  )
  vignette <- "Faraj_2022_marzeptacogAlfa"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  # Dose records: subcutaneous doses go to `depot` with rate = -2 so that
  # rxode2 applies the modelled zero-order duration dur(depot); intravenous
  # doses go to `central`. Doses are in ug (ug/kg dose times body weight).
  compartmentData <- list(
    depot = list(
      analyte = "marzeptacog alfa (activated)",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "marzeptacog alfa (activated)", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "marzeptacog alfa (activated)", units = "ug", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of CL and Q (exponent fixed at 0.75) and Vc and Vp (exponent fixed at 1), reference 70 kg (Faraj 2022 Methods 'Model building' equations; Data S1 $PK). Body weight also converts the per-kg dose into the absolute dose in ug.",
      source_name = "WT"
    ),
    OCC = list(
      description = "Integer-valued occasion indicator for inter-occasion variability.",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Faraj 2022 Data S1 multiplexes the inter-occasion etas on the data column REGI: REGI = 1 carries no IOV, and REGI = 2..9 each carry their own IOV eta on D1, F and Vc with one shared variance per parameter ($OMEGA BLOCK(1) SAME). OCC here takes the REGI value directly: OCC = 1 switches all IOV off, OCC = 2..9 select the occasion slot. The paper does not define REGI further; the IOV is described as applying to subcutaneous data obtained on more than one occasion per individual (Methods).",
      source_name = "REGI"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Age on ka was statistically significant in the stepwise covariate search (p < 0.01) but changed the typical Cmax after 60 ug/kg by +18% / -22% at the 10th / 90th age percentiles (23 / 52 years) versus the median 32 years, below the 30% clinical-significance criterion, so it was not retained (Faraj 2022 Results 'Model development'). No coefficient is printed."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 46L,
    n_studies = 3L,
    n_observations = 1225L,
    age_range = "18-62 years (10th-90th percentile 23-52 years)",
    age_median = "32 years",
    weight_range = "43-120 kg (10th-90th percentile 60-94 kg)",
    weight_median = "75 kg",
    sex_female_pct = 0,
    race_ethnicity = "Not reported",
    disease_state = "Haemophilia A or B, with or without inhibitors (NCT01439971: severe haemophilia A/B with inhibitors)",
    dose_range = "Single i.v. 4.5-30 ug/kg; single s.c. 30-120 ug/kg; s.c. 60 ug/kg two or three times 3 h apart; daily s.c. 30 ug/kg (escalated to 60 ug/kg) for 50 days",
    regions = "Not reported",
    trials = "NCT01439971 (phase I, i.v.), NCT04072237 (phase I/IIa, i.v. and s.c.), NCT03407651 (phase II, i.v. and s.c. daily prophylaxis)",
    notes = "Faraj 2022 Methods 'Study design and patients' and Table 1. 19 BLOQ observations (1.6%) and 3 further samples were excluded. Observations are functional FVIIa:Clot activity calibrated against MarzAA (Data S2), so the endogenous FVIIa baseline is included in each observation."
  )

  ini({
    # Structural parameters: Faraj 2022 Table S1 (final estimates) for a
    # 70 kg subject. Table S1 prints CL and Q in mL/h and V in mL; they are
    # converted to L/h and L so that ug doses give ug/L = ng/mL.
    lcl <- log(1.085); label("Clearance at 70 kg (L/h)") # Table S1: CL = 1085 mL/h (RSE 7.5%)
    lvc <- log(3.393); label("Central volume of distribution at 70 kg (L)") # Table S1: Vc = 3393 mL (RSE 6.5%)
    lvp <- log(0.619); label("Peripheral volume of distribution at 70 kg (L)") # Table S1: Vp = 619 mL (RSE 20%)
    lq <- log(0.213); label("Inter-compartmental clearance at 70 kg (L/h)") # Table S1: Q = 213 mL/h (RSE 39%)
    lrbase <- log(0.96); label("Endogenous baseline FVIIa activity (ng/mL)") # Table S1: BASE = 0.96 ng/mL (RSE 12%)
    lka <- log(0.05); label("First-order absorption rate constant from the depot (1/h)") # Table S1: ka = 0.05 1/h (RSE 15%)
    logitfdepot <- log(0.26 / (1 - 0.26)); label("Logit of subcutaneous bioavailability (logit fraction)") # Table S1: F = 0.26 (RSE 16%); Data S1 PHI = LOG(TVF1/(1-TVF1))
    ld1 <- log(0.875); label("Duration of zero-order input into the depot (h)") # Table S1: D1 = 0.875 h (RSE 24%)

    # Allometric exponents fixed (Methods 'Model building'; Data S1 $PK).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL and Q (unitless)") # Methods: fixed exponent 0.75 on clearance terms
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vc and Vp (unitless)") # Methods: fixed exponent 1 on volume terms

    # Inter-individual variability. Table S1 prints random effects as
    # percentages; they are read as 100 * omega (SD of the eta), so the
    # variance is (value / 100)^2. This is the only reading that applies to
    # the logit-scale F eta, whose percentage cannot be a log-normal CV.
    # Table S1 'Covariance between CL and V (%)' = 92 is read as the
    # correlation 0.92: covariance = 0.92 * 0.28 * 0.26 = 0.066976.
    etalcl + etalvc ~ c(0.0784, 0.066976, 0.0676) # Table S1: IIV CL 28%, correlation 92%, IIV V 26%
    etalrbase ~ 0.3136 # Table S1: IIV baseline FVIIa 56% -> 0.56^2
    etalogitfdepot ~ 0.2916 # Table S1: IIV F 54% (logit scale) -> 0.54^2
    etalka ~ 0.09 # Table S1: IIV ka 30% -> 0.30^2

    # Inter-occasion variability, one shared variance per parameter across
    # occasions 2-9 (Data S1 $OMEGA BLOCK(1) SAME); occasions 3-9 are fixed
    # to the occasion-2 value.
    etaiov_ld1_2 ~ 0.3721 # Table S1: IOV D1 61% -> 0.61^2
    etaiov_ld1_3 ~ fixed(0.3721)
    etaiov_ld1_4 ~ fixed(0.3721)
    etaiov_ld1_5 ~ fixed(0.3721)
    etaiov_ld1_6 ~ fixed(0.3721)
    etaiov_ld1_7 ~ fixed(0.3721)
    etaiov_ld1_8 ~ fixed(0.3721)
    etaiov_ld1_9 ~ fixed(0.3721)
    etaiov_logitfdepot_2 ~ 0.0441 # Table S1: IOV F 21% (logit scale) -> 0.21^2
    etaiov_logitfdepot_3 ~ fixed(0.0441)
    etaiov_logitfdepot_4 ~ fixed(0.0441)
    etaiov_logitfdepot_5 ~ fixed(0.0441)
    etaiov_logitfdepot_6 ~ fixed(0.0441)
    etaiov_logitfdepot_7 ~ fixed(0.0441)
    etaiov_logitfdepot_8 ~ fixed(0.0441)
    etaiov_logitfdepot_9 ~ fixed(0.0441)
    etaiov_lvc_2 ~ 0.2025 # Table S1: IOV V 45% -> 0.45^2
    etaiov_lvc_3 ~ fixed(0.2025)
    etaiov_lvc_4 ~ fixed(0.2025)
    etaiov_lvc_5 ~ fixed(0.2025)
    etaiov_lvc_6 ~ fixed(0.2025)
    etaiov_lvc_7 ~ fixed(0.2025)
    etaiov_lvc_8 ~ fixed(0.2025)
    etaiov_lvc_9 ~ fixed(0.2025)

    # Residual error: Data S1 Y = IPRED*EXP(EPS(1)) without log-transformed
    # data, i.e. proportional (Results 'Model development').
    propSd <- 0.29; label("Proportional residual error (fraction)") # Table S1: proportional error 29% (RSE 3.1%)
  })

  model({
    # Occasion indicators (Data S1: REGI = 1 has no IOV; REGI = 2..9 each
    # select an IOV slot).
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    oc8 <- (OCC == 8)
    oc9 <- (OCC == 9)
    iov_ld1 <- oc2 * etaiov_ld1_2 + oc3 * etaiov_ld1_3 + oc4 * etaiov_ld1_4 +
      oc5 * etaiov_ld1_5 + oc6 * etaiov_ld1_6 + oc7 * etaiov_ld1_7 +
      oc8 * etaiov_ld1_8 + oc9 * etaiov_ld1_9
    iov_logitfdepot <- oc2 * etaiov_logitfdepot_2 + oc3 * etaiov_logitfdepot_3 +
      oc4 * etaiov_logitfdepot_4 + oc5 * etaiov_logitfdepot_5 +
      oc6 * etaiov_logitfdepot_6 + oc7 * etaiov_logitfdepot_7 +
      oc8 * etaiov_logitfdepot_8 + oc9 * etaiov_logitfdepot_9
    iov_lvc <- oc2 * etaiov_lvc_2 + oc3 * etaiov_lvc_3 + oc4 * etaiov_lvc_4 +
      oc5 * etaiov_lvc_5 + oc6 * etaiov_lvc_6 + oc7 * etaiov_lvc_7 +
      oc8 * etaiov_lvc_8 + oc9 * etaiov_lvc_9

    # Individual parameters (Data S1 $PK). Q and Vp carry no IIV.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc + iov_lvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    rbase <- exp(lrbase + etalrbase)
    ka <- exp(lka + etalka)
    fdepot <- expit(logitfdepot + etalogitfdepot + iov_logitfdepot)
    d1 <- exp(ld1 + iov_ld1)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Subcutaneous input: zero-order release of duration D1 into the depot
    # (dose records need rate = -2); bioavailability applies to the depot
    # only, as F1 does in Data S1.
    f(depot) <- fdepot
    dur(depot) <- d1

    # Observed FVIIa activity = endogenous baseline + drug concentration
    # (Data S1 IPRED = BASE + A(2)/S2); ug / L = ng/mL.
    Cc <- rbase + central / vc
    Cc ~ prop(propSd)
  })
}
