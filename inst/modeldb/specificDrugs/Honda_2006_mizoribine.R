Honda_2006_mizoribine <- function() {
  description <- "One-compartment oral PK model for mizoribine in healthy Caucasian male volunteers (Honda 2006); first-order absorption with a fixed absorption-lag time, apparent volume of distribution V/F linear in body weight, apparent oral clearance CL/F linear in Cockcroft-Gault creatinine clearance (CLcr), and additive residual error on the serum-concentration scale."
  reference <- "Honda M, Itoh H, Suzuki T, Hashimoto Y. Population pharmacokinetics of higher-dose mizoribine in healthy male volunteers. Biol Pharm Bull. 2006 Dec;29(12):2460-2464. doi:10.1248/bpb.29.2460"
  vignette <- "Honda_2006_mizoribine"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear scaling on apparent volume of distribution: V/F = theta_3 * WT (Honda 2006 Eq. 3). Cohort range 54.4-98.2 kg, mean 75.6 kg (Honda 2006 Materials and Methods, Pharmacokinetic Data).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age (baseline)",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the Cockcroft-Gault formula for creatinine clearance: CLcr (L/h) = ((140 - AGE) * WT / (72 * CREAT)) * (60 / 1000), with CREAT in mg/dL (Honda 2006 Eq. 5). Cohort range 18-45 years, mean 25.8.",
      source_name        = "AGE"
    ),
    CREAT = list(
      description        = "Serum creatinine (baseline)",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the Cockcroft-Gault formula for creatinine clearance (Honda 2006 Eq. 5). Individual serum-creatinine values are not tabulated; the derived cohort-mean CLcr was 7.54 (+/- 1.40) L/h, i.e. ~125.7 mL/min.",
      source_name        = "Scr"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 36L,
    n_studies      = 2L,
    age_range      = "18-45 years",
    age_median     = "mean 25.8 years",
    weight_range   = "54.4-98.2 kg",
    weight_median  = "mean 75.6 kg",
    sex_female_pct = 0,
    race_ethnicity = c(Caucasian = 100),
    disease_state  = "Healthy adult male volunteers with normal renal function.",
    dose_range     = "Oral mizoribine (50 mg tablets). Single-dose study: 3, 6, 9, or 12 mg/kg (24 subjects across 4 dose groups of 6). Multiple-dose study: 6 mg/kg once daily for 5 days, or 6 mg/kg every 12 h for 7 days (12 subjects across 2 groups of 6).",
    regions        = "Phase 1 clinical-pharmacology cohort. Underlying serum-concentration data were from a prior phase 1 study (Stypinski et al. 2005, Br J Clin Pharmacol, in press at the time of Honda 2006).",
    renal_function = "Cockcroft-Gault creatinine clearance mean 7.54 (+/- 1.40) L/h across all 36 subjects.",
    notes          = "Demographics summarised in the Pharmacokinetic Data section of Materials and Methods. The PK data were originally obtained in the previous phase 1 study and re-analysed here with NONMEM (first-order conditional estimation) using PREDPP subroutines ADVAN2 and TRANS2 (one-compartment with first-order absorption)."
  )

  ini({
    # Structural parameters - Honda 2006 Table 1 final-model point estimates.
    # The reference subject is a Caucasian male with the cohort's median weight
    # and Cockcroft-Gault CLcr; the V/F and CL/F coefficients below combine
    # with WT and CLcr inside model() to yield individual-specific values.
    ltlag <- log(0.349); label("Absorption lag time (h)")                                                   # Honda 2006 Table 1 theta_1 = 0.349 h (no IIV per text: "we did not introduce interindividual variability in ALAG for simplicity")
    lka   <- log(0.869); label("Absorption rate constant ka (1/h)")                                         # Honda 2006 Table 1 theta_2 = 0.869 1/h
    lvc   <- log(0.834); label("Apparent volume of distribution per kg body weight, (V/F)/WT (L/kg)")       # Honda 2006 Table 1 theta_3 = 0.834 L/kg; V/F = theta_3 * WT (Eq. 3)
    lcl   <- log(1.93);  label("Apparent oral clearance per L/h Cockcroft-Gault CLcr, (CL/F)/CLcr (L/h per L/h)") # Honda 2006 Table 1 theta_4 = 1.93; CL/F = theta_4 * CLcr (Eq. 4)

    # Inter-individual variability - Honda 2006 Table 1 reports omega_KA,
    # omega_V/F, and omega_CL/F (i.e. the SDs on the log scale), distinct from
    # the omega^2 (variance) notation used inside the model equations
    # (Eqs. 2-4: "variance of omega^2_KA / omega^2_V/F / omega^2_CL/F"). The
    # same authors' Honda 2005 carvedilol paper (Biol Pharm Bull 28(9): 1683-7)
    # tabulates omega^2 explicitly (e.g. "omega^2_CL/F = 0.130, CV = 36.1%",
    # CV = sqrt(omega^2) * 100); Honda 2006 drops the squared notation in
    # the results table and tabulates the SD instead. The variances used in
    # ini() therefore equal the squared Table 1 values, which also reproduces
    # the typical-CV magnitudes seen in the Bayesian individual estimates
    # (Table 2): KA ~ 11-25%, V/F ~ 23-39%.
    etalka ~ 0.0724  # = 0.269^2; Honda 2006 Table 1 omega_KA = 0.269 (log-scale SD; CV ~ 27%)
    etalvc ~ 0.1253  # = 0.354^2; Honda 2006 Table 1 omega_V/F = 0.354 (log-scale SD; CV ~ 36%)
    etalcl ~ 0.1102  # = 0.332^2; Honda 2006 Table 1 omega_CL/F = 0.332 (log-scale SD; CV ~ 34%)

    # Residual error - additive in linear concentration space; Honda 2006 Eq. 6:
    # C_ij = C*_ij + epsilon_ij with epsilon ~ N(0, sigma^2). Table 1 reports
    # sigma (the SD) directly in ug/mL, again paired with an explicit sigma^2
    # notation in the equation.
    addSd <- 0.352; label("Additive residual error SD (ug/mL)")                                              # Honda 2006 Table 1 sigma = 0.352 ug/mL
  })

  model({
    # Cockcroft-Gault creatinine clearance, expressed in L/h, exactly as
    # Honda 2006 Eq. 5: CLcr (L/h) = ((140 - AGE) * WT / (72 * CREAT)) *
    # (60 / 1000), with CREAT in mg/dL. The 60/1000 factor converts the raw
    # Cockcroft-Gault output from mL/min to L/h. The cohort is all male, so
    # the female correction factor (* 0.85) is not applied.
    crcl <- ((140 - AGE) * WT / (72 * CREAT)) * (60 / 1000)

    # Individual PK parameters
    tlag <- exp(ltlag)
    ka   <- exp(lka + etalka)
    vc   <- exp(lvc + etalvc) * WT        # V/F = theta_3 * WT * exp(eta_V/F); Honda 2006 Eq. 3
    cl   <- exp(lcl + etalcl) * crcl      # CL/F = theta_4 * CLcr * exp(eta_CL/F); Honda 2006 Eq. 4

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Absorption lag time on the depot compartment (no IIV).
    alag(depot) <- tlag

    # Dose in mg, vc in L -> Cc in mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
