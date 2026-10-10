Li_2022_gumokimab <- function() {
  description <- paste(
    "Sequential population PK/PD model for gumokimab (AK111, a humanized",
    "anti-IL-17A IgG1 monoclonal antibody) in Chinese adults with",
    "moderate-to-severe plaque psoriasis (Li 2022 phase 1b).",
    "PK: one-compartment disposition with first-order SC absorption and",
    "first-order elimination (apparent CL 0.182 L/day, V 6.65 L).",
    "PD: indirect-response model for the Psoriasis Area and Severity Index",
    "(PASI) in which serum gumokimab inhibits plaque formation (Imax fixed at",
    "1, IC50 0.52 ug/mL) and a placebo effect multiplies the plaque-loss rate",
    "by 1 + PLBmax (Kplb fixed at 0, so the placebo effect is constant from",
    "the first dose). Baseline PASI is the steady state kin / kout, and the",
    "percentage of body surface area affected by psoriasis enters kout as a",
    "power covariate normalized to 31%."
  )
  reference <- paste(
    "Li Q, Qiao J, Jin H, Chen B, He Z, Wang G, Ni X, Wang M, Xia M, Li B,",
    "Chen R, Hu P. Population pharmacokinetic/pharmacodynamic analysis of",
    "AK111, an IL-17A monoclonal antibody, in subjects with",
    "moderate-to-severe plaque psoriasis. Front Pharmacol. 2022;13:966176.",
    "doi:10.3389/fphar.2022.966176"
  )
  vignette <- "Li_2022_gumokimab"

  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(
      analyte = "gumokimab",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "gumokimab",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    ),
    pasi = list(
      analyte = "none",
      units = "PASI units (0-72 clinical score)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    BSA_AFFECTED_PCT = list(
      description = "Percentage of total body surface area affected by psoriasis at baseline",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on kout normalized to 31% (Li 2022 Eq 7,",
        "(BSA/31)^thetaBSA). The paper's 'BSA (%)' / 'body surface area (BSA)",
        "involvement' is the extent of psoriatic skin, NOT the body surface",
        "area in m^2 (canonical BSA). Cohort mean 33.3% (SD 14.7%), Table 1."
      ),
      source_name = "BSA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 47L,
    n_studies = 1L,
    n_observations = "516 serum PK samples (AK111 arms) and 344 PASI scores (AK111 and placebo arms)",
    age_range = "mean 38.2 years (SD 8.9)",
    weight_range = "mean 67.5 kg (SD 10.0)",
    sex_female_pct = 21.3,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Moderate-to-severe plaque psoriasis; mean baseline PASI 20 (SD 5.6),",
      "mean BSA affected 33.3% (SD 14.7%), mean disease duration 15.7 years;",
      "8.5% with prior biologic use."
    ),
    dose_range = paste(
      "75, 150, 300 or 450 mg SC (or placebo) at weeks 0, 1, 4 and 8;",
      "9 AK111 and 3 placebo subjects per dose cohort."
    ),
    regions = "China (single centre, Peking Union Medical College Hospital)",
    notes = paste(
      "Li 2022 Table 1 and Methods (phase 1b, randomized, double-blind,",
      "placebo-controlled). 48 enrolled; one 450 mg subject was excluded for",
      "delayed administration, leaving 47 for the PK/PD analysis. PK samples",
      "to day 141; PASI on days 8, 15, 29, 57, 85, 113 and 141. NONMEM 7.3,",
      "FOCE-I, sequential PK then PD (individual PK post hoc estimates fed",
      "the PD fit). No covariates were retained on PK."
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested on PK and PD parameters (Li 2022 Methods); not retained."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on PK and PD parameters (Li 2022 Methods); not retained."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Tested on PK and PD parameters (Li 2022 Methods); not retained."
    )
  )

  ini({
    # PK (Li 2022 Table 2, PK model rows)
    lka <- log(0.463); label("First-order SC absorption rate constant Ka (1/day)") # Table 2 Ka = 0.463 1/day (RSE 9.3%)
    lcl <- log(0.182); label("Apparent clearance CL (L/day)") # Table 2 CL = 0.182 L/day (RSE 7.2%)
    lvc <- log(6.65); label("Apparent volume of distribution V (L)") # Table 2 V = 6.65 L (RSE 7.8%); row text says 'peripheral', abstract says central

    # PD (Li 2022 Table 2, PD model rows; Eqs 6-8)
    lkin <- log(0.474); label("Zero-order rate constant of psoriatic plaque formation Kin (PASI/day)") # Table 2 Kin = 0.474 PASI/day (RSE 21.3%)
    lkout <- log(0.024); label("First-order rate constant of psoriatic plaque loss Kout at 31% BSA affected (1/day)") # Table 2 Kout = 0.024 1/day (RSE 21.5%)
    imax <- fixed(1); label("Maximum fractional inhibition of plaque formation Imax (unitless)") # Table 2 Imax = 1, Fixed
    lic50 <- log(0.52); label("Serum concentration at half-maximal inhibition IC50 (ug/mL)") # Table 2 IC50 = 0.52 ug/mL (RSE 66.4%)
    lpmax <- log(0.429); label("Maximum fractional increase of Kout by placebo PLBmax (unitless)") # Table 2 PLBmax = 0.429 (RSE 78.2%)
    kplb <- fixed(0); label("Decay rate constant of the placebo effect Kplb (1/day)") # Table 2 Kplb = 0, Fixed
    e_bsa_affected_pct_kout <- -0.572; label("Power exponent of (BSA affected / 31%) on Kout (unitless)") # Table 2 thetaBSA = -0.572

    # IIV: exponential (Eq 1). Table 2 reports IIV as % with no stated
    # convention; read as omega x 100 (the SD of eta), so omega^2 =
    # (IIV/100)^2. The Figures 10-13 simulated PASI90 plateaus support this
    # reading over a CV% read with omega^2 = log(1 + CV^2) (see vignette).
    etalka ~ 0.251001 # Table 2 IIV Ka 50.1% -> 0.501^2
    etalcl ~ 0.178084 # Table 2 IIV CL 42.2% -> 0.422^2
    etalvc ~ 0.215296 # Table 2 IIV V 46.4% -> 0.464^2
    etalkin ~ 0.0256 # Table 2 IIV Kin 16.0% -> 0.160^2 (shrinkage 98.8%)
    etalkout ~ 0.054289 # Table 2 IIV Kout 23.3% -> 0.233^2
    etalic50 ~ 2.598544 # Table 2 IIV IC50 161.2% -> 1.612^2
    etalpmax ~ 0.974169 # Table 2 IIV PLBmax 98.7% -> 0.987^2

    # Residual error: Table 2 sigma rows are NONMEM $SIGMA variances, so
    # each SD is the square root of the printed value.
    propSd <- 0.12; label("Proportional residual error for serum concentration (fraction)") # Table 2 sigma prop 0.0144 = 0.12^2
    addSd <- 1.2166; label("Additive residual error for serum concentration (ug/mL)") # Table 2 sigma addi 1.48 = 1.2166^2
    addSd_pasi <- 2.1284; label("Additive residual error for PASI score (PASI units)") # Table 2 sigma PASI 4.53 = 2.1284^2
  })

  model({
    # Individual PK parameters (Eq 1)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    # Individual PD parameters (Eqs 1, 7, 8)
    kin <- exp(lkin + etalkin)
    kout <- exp(lkout + etalkout) * (BSA_AFFECTED_PCT / 31)^e_bsa_affected_pct_kout
    ic50 <- exp(lic50 + etalic50)
    pmax <- exp(lpmax + etalpmax)

    # Placebo effect multiplies kout for every subject (Eq 8); t is time
    # since the first dose. With kplb = 0 the multiplier is 1 + pmax.
    plb <- 1 + pmax * exp(-kplb * t)

    Cc <- central / vc

    # Eqs 4-6
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(pasi) <- kin * (1 - imax * Cc / (ic50 + Cc)) - kout * plb * pasi

    # Baseline PASI at the drug- and placebo-free steady state
    pasi(0) <- kin / kout

    Cc ~ add(addSd) + prop(propSd)
    pasi ~ add(addSd_pasi)
  })
}
