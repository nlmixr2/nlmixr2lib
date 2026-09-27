Khalil_2020_tapentadol_allometryFixed <- function() {
  description <- "One-compartment population PK model for tapentadol oral solution and 1-hour intravenous infusion in 148 children from birth (including preterm neonates) to <18 years with acute pain (Khalil 2020), the variant with the allometric exponents FIXED to the theoretical 0.75 (CL) and 1 (V). First-order absorption from a depot preceded by an absorption lag time, first-order elimination. CL and V are systemic (IV data identify F). CL carries a hyperbolic (Hill fixed to 1) postmenstrual-age maturation function; oral bioavailability decays exponentially with postnatal age from twice the adult value at birth to the adult value. Log-normal IIV on CL, V (correlated) and Ka; combined proportional plus additive residual error. The companion Khalil_2020_tapentadol_allometryEstimated model estimates the two exponents instead."
  reference <- paste(
    "Khalil F, Choi SL, Watson E, Tzschentke TM, Lefeber C, Eerdekens M, Freijer J.",
    "Population Pharmacokinetics of Tapentadol in Children from Birth to <18 Years Old.",
    "J Pain Res. 2020;13:3107-3123.",
    "doi:10.2147/JPR.S269549",
    sep = " "
  )
  vignette <- "Khalil_2020_tapentadol"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  compartmentData <- list(
    depot = list(analyte = "tapentadol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tapentadol", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric size scaling (WT/70)^0.75 on CL and (WT/70)^1 on V, reference weight 70 kg",
        "(Khalil 2020 Methods equations 1-2; Table 3 caption). Cohort range 1.6-80 kg",
        "(Table 1). Treated as time-fixed; every trial was single-dose."
      ),
      source_name = "WT (Khalil 2020 equations 1-2)"
    ),
    PAGE = list(
      description = "Postmenstrual age (gestational age at birth plus postnatal age)",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "WEEKS, not the register-default months: Khalil 2020 writes the CL maturation",
        "function PMA / (PMA + PMA50) with PMA and PMA50 in weeks (Methods equation 1;",
        "Table 3 'PMA50 (wks)'). PMA = gestational age + postnatal age, with gestational",
        "age 40 weeks for term-born children (Methods, 'Simulation of Virtual Pediatric",
        "Population'). For a term-born child PAGE = 40 + PNA_months * 30.4375 / 7."
      ),
      source_name = "PMA (Khalil 2020 equation 1)"
    ),
    PNA = list(
      description = "Postnatal age (chronological time since birth)",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the exponential decay of oral bioavailability F = Fadult * (1 + exp(-k * PNA))",
        "(Methods equation 3). The paper states PNA and k in weeks; the canonical PNA is in",
        "months, so model() converts with PNA_weeks = PNA * 30.4375 / 7. Only affects oral",
        "doses (depot)."
      ),
      source_name = "PNA (Khalil 2020 equation 3)"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator; 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Tested for a univariate influence on CL or V and not statistically significant (Methods, 'Model Development'; Results). No point estimate reported.",
      source_name = "sex"
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested for a univariate influence on CL or V and not statistically significant (Methods, 'Model Development'; Results). No point estimate or units reported; recorded as the register default.",
      source_name = "creatinine clearance"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 148L,
    n_studies = 4L,
    n_observations = "569 quantifiable tapentadol serum concentrations (Khalil 2020 Methods, 'Population Pharmacokinetic Modeling Data Set')",
    age_range = "Birth (including preterm neonates, gestational age >=24 weeks) to <18 years; group medians 15 y, 9 y, 3 y, 14 months, 3 months, 12 days (term neonates) and 12 days (preterm) (Table 1)",
    weight_range = "1.6-80 kg; group medians 59.7, 29.5, 16.4, 10, 5.9, 3.6 and 2.45 kg (Table 1)",
    sex_female_pct = 45.9,
    disease_state = "Acute pain (mostly postsurgical or procedural) severe enough to require an opioid",
    dose_range = "Single dose. Oral solution 1.0 mg/kg (>=2 years), 0.75 mg/kg (6 months to <2 years), 0.6 mg/kg (1 to <6 months), 0.5 mg/kg (birth to <1 month); IV 1-hour infusion 0.4 mg/kg (<2 years) or 0.3-0.4 mg/kg (preterm, by gestational and postnatal age) (Table 1)",
    regions = "Multinational (trial sites listed in the paper's supplement)",
    notes = paste(
      "Pooled from four single-dose open-label phase 2 PK trials: NCT01729728 (n=56) and",
      "NCT01134536 (n=36), oral solution, 2 to <18 years; NCT02221674 (n=18), oral solution,",
      "birth to <2 years (term); EudraCT 2014-002259-24 (n=38), IV, preterm to <2 years",
      "(Table 2). Preterm neonates were dosed IV only. Female percentage (45.9%) is the",
      "N-weighted sum of Table 1's per-group percentages (68/148; Table 2 gives the same 68). 38 oral samples from 17",
      "patients who vomited within 3 h or did not take the full dose, and one contaminated IV",
      "sample, were excluded. NONMEM 7.2."
    )
  )

  ini({
    # ========================================================================
    # Structural parameters (Khalil 2020 Table 3, 'Final Model with Fixed
    # Exponents'). One-compartment model with first-order absorption after a
    # lag time. CL and V are SYSTEMIC, typical values at 70 kg and full
    # maturation: the IV arm identifies F, and the Discussion recovers
    # CL/F = 94.6 / 0.349 = 271 L/h and V/F = 414 / 0.349 = 1186 L.
    # ========================================================================
    lcl <- log(94.6); label("Clearance CL at 70 kg, full maturation (L/h)") # Table 3 fixed-exponent model: CL = 94.6 L/h (RSE 10.5%; bootstrap 76.6-117.3)
    lvc <- log(414); label("Volume of distribution V at 70 kg (L)") # Table 3 fixed-exponent model: V = 414 L (RSE 9%; bootstrap 342-497)
    lka <- log(2.19); label("Absorption rate constant Ka (1/h)") # Table 3 fixed-exponent model: Ka = 2.19 1/h (RSE 14.7%; bootstrap 1.65-3.25)
    lfdepot <- log(0.349); label("Oral bioavailability at full maturation, Fadult (fraction)") # Table 3 fixed-exponent model: F = 0.349 (RSE 10.1%; bootstrap 0.283-0.432)
    ltlag <- log(0.266); label("Absorption lag time TLAG (h)") # Table 3 fixed-exponent model: TLAG = 0.266 h (RSE 0.7%; bootstrap 0.245-0.291)

    # ---- Maturation ---------------------------------------------------------
    # CL = CLTV * PMA / (PMA + PMA50) * (WT/70)^n  (Methods equation 1); the
    # Hill exponent on PMA is fixed to 1 in both final models.
    ltm50_cl <- log(34.8); label("Postmenstrual age at 50% CL maturation, PMA50 (weeks)") # Table 3 fixed-exponent model: PMA50 = 34.8 wks (RSE 20.8%; bootstrap 22.4-51.2)
    hill_mat <- fixed(1); label("Hill exponent of the CL maturation function (unitless)") # Table 3: HILL exponent = 1 FIXED
    # F = Fadult * (1 + exp(-k * PNA))  (Methods equation 3; the minus sign is
    # dropped by text extraction and confirmed on the rendered equation).
    e_pna_fdepot <- 0.122; label("Rate constant k of the postnatal-age decay of oral bioavailability (1/week)") # Table 3 fixed-exponent model: k = 0.122 wks^-1 (RSE 29.1%; bootstrap 0.068-1.657)

    # ---- Allometric exponents (fixed to theory) -----------------------------
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)") # Table 3: Exponent CL-WT = 0.75 FIXED
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V (unitless)") # Table 3: Exponent V-WT = 1 FIXED

    # ========================================================================
    # IIV: Pi = PTV * exp(eta_i) (Methods equation 4). Table 3 reports omega^2
    # (variances) directly. CL-V block correlation 0.0867/sqrt(0.0961*0.13) = 0.78.
    # ========================================================================
    etalcl + etalvc ~ c(
      0.0961, # Table 3 fixed-exponent model: IIV CL (omega^2) = 0.0961 (RSE 21.1%, 15.9% shrinkage)
      0.0867, 0.13 # Table 3 fixed-exponent model: Cov CL-V = 0.0867 (bootstrap 0.041-0.207); IIV V (omega^2) = 0.13 (RSE 30.5%, 22.4% shrinkage)
    )
    etalka ~ 2 # Table 3 fixed-exponent model: IIV Ka (omega^2) = 2 (RSE 27.9%, 26.8% shrinkage)

    # ========================================================================
    # Residual error: Co = Cp * (1 + eps_p) + eps_a (Methods equation 5).
    # Table 3's abbreviation list defines sigma as a standard deviation.
    # ========================================================================
    propSd <- 0.327; label("Proportional residual error SD (fraction)") # Table 3 fixed-exponent model: Proportional error (sigma) = 0.327 (RSE 14.9%)
    addSd <- 0.48; label("Additive residual error SD (ng/mL)") # Table 3 fixed-exponent model: Additive error = 0.48 ng/mL (RSE 47.1%)
  })

  model({
    # Postnatal age in weeks (the paper's unit) from the canonical months.
    pna_wk <- PNA * 30.4375 / 7

    # Methods equation 1: hyperbolic PMA maturation times allometric weight.
    tm50_cl <- exp(ltm50_cl)
    mat_cl <- PAGE^hill_mat / (PAGE^hill_mat + tm50_cl^hill_mat)
    cl <- exp(lcl + etalcl) * mat_cl * (WT / 70)^e_wt_cl
    # Methods equation 2.
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)
    # Methods equation 3: F is twice Fadult at birth and decays to Fadult.
    fdepot <- exp(lfdepot) * (1 + exp(-e_pna_fdepot * pna_wk))

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- fdepot
    alag(depot) <- tlag

    # Dose in mg, V in L -> mg/L = ug/mL; x1000 for ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}
