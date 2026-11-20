Johnson_2025_levetiracetam <- function() {
  description <- paste(
    "One-compartment reduction of the Johnson 2025 neonatal",
    "physiologically based pharmacokinetic (PBPK) model of intravenous",
    "levetiracetam in term and preterm neonates (Simcyp V23 full PBPK).",
    "Volume is the paper's neonatal steady-state volume (0.96 L/kg,",
    "Kp scalar 1.56). Total clearance is the sum of the paper's two",
    "neonatal clearance equations: renal clearance scaled from the adult",
    "value by the Wu 2022 neonatal GFR model (birth weight, gestational",
    "age and postnatal age) times a maturation factor of 1.61, and",
    "metabolic (Type B esterase) clearance scaled allometrically",
    "(exponent 0.75) from the adult value times a quadratic ontogeny",
    "function of postmenstrual age, floored at zero. Typical-value",
    "model: no between-subject or residual variability."
  )
  reference <- paste(
    "Johnson A, Thomas N, Blumenthal M, Ikonomidou C, Lim SY.",
    "A Physiologically Based Modeling Approach to Evaluate Intravenous",
    "Levetiracetam Dosing in Term and Preterm Neonates.",
    "J Clin Pharmacol. 2025;65(10):1246.",
    "doi:10.1002/jcph.70037"
  )
  vignette <- "Johnson_2025_levetiracetam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "levetiracetam",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Current body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. Scales the volume (0.96 L/kg) and the allometric",
        "metabolic clearance (WT/70)^0.75 (Johnson 2025 Table 1",
        "footnote c, Equation 1). The paper's Simcyp simulations",
        "redefine body weight over time ('Redefine Subjects Over Time')."
      ),
      source_name = "WT"
    ),
    WT_BIRTH = list(
      description = "Birth weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed. Enters the Wu 2022 neonatal GFR model as bW / 1.75",
        "(Johnson 2025 Methods, 'PBPK Model Scaling to Term Neonates').",
        "The paper used the Simcyp preterm-population average birth",
        "weight for each gestational age, which is not printed."
      ),
      source_name = "bW"
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed. Enters the Wu 2022 GFR model as (GA/34)^9.99 and",
        "through the postnatal-age exponent (GA/34)^-1.71, and the",
        "ontogeny function through postmenstrual age = GA + PNA."
      ),
      source_name = "GA"
    ),
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying; supply it on every record (it changes",
        "clearance materially within the first two weeks of life).",
        "Canonical PNA is in months; the source equations use days",
        "(Wu 2022 GFR term (PNA/2)) and weeks (postmenstrual age), which",
        "model() recovers with PNA * 30.4375 days per month."
      ),
      source_name = "PNA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 4L,
    age_range = "Neonates, postnatal age 0-32 days (validation studies); simulations at postnatal age 1-14 days",
    ga_range = "28-41 weeks gestational age",
    weight_range = "0.9-4.7 kg (validation studies, Johnson 2025 Table 2)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported (Simcyp Sim-Pediatric and Sim-Preterm virtual populations)",
    disease_state = paste(
      "Term and preterm neonates, most with neonatal seizures in the",
      "validation studies (Sharpe 2012 and 2024, Merhar 2011,",
      "Lima-Rogel 2018); one cohort had 65% hypoxic-ischaemic",
      "encephalopathy (Lima-Rogel)."
    ),
    dose_range = "Intravenous 15-min infusions: loading doses 20-150 mg/kg, maintenance 5-50 mg/kg/day (validation) and 25 mg/kg every 8-12 h (simulated regimens)",
    regions = "USA, New Zealand, Mexico (validation studies)",
    notes = paste(
      "PBPK model: no subjects were fitted. The neonatal layer was",
      "calibrated to Sharpe 2012 (term, n = 18) and checked against",
      "Sharpe 2024 (term, n = 49), Merhar 2011 (term and late",
      "preterm, n = 18) and Lima-Rogel 2018 (term and preterm, GA",
      "28-40 weeks, n = 20) (Johnson 2025 Table 2). Virtual cohorts",
      "were n = 500 per simulation in Simcyp V23."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every value is fixed: Johnson 2025 is a Simcyp PBPK model and
    # nothing here was estimated against data. The only value not
    # printed in the paper is the adult GFR reference (lgfr_adult),
    # which is back-solved from the paper's own Table S6 (see below).
    # ------------------------------------------------------------------

    # Volume: neonatal Vss per kg of current body weight.
    lvc <- fixed(log(0.96))
    label("Neonatal steady-state volume of distribution per kg body weight (L/kg)") # Table 1 'Vss (L/kg)', Neonatal Model column = 0.96 (footnote b: Kp scalar 1.56)

    # Adult clearance components (Simcyp base model inputs).
    lcl_renal <- fixed(log(2.64))
    label("Adult renal clearance CL_R,adult (L/h)") # Table 1 'Renal CL (L/h)' = 2.64; Methods 'renal clearance ... for adults was set to 2.64 L/h'
    lcl_nonren <- fixed(log(1.36))
    label("Adult metabolic (Type B esterase) clearance CL_M,adult (L/h)") # Table 1 'Metabolic CL (L/h)' = 1.36; Equation 1 coefficient 1.36

    # Renal-clearance scaling (Table 1 footnote e, Equation 2):
    #   CL_R,neonate = CL_R,adult x GFR_neonate / GFR_adult x MF
    mf_renal <- fixed(1.61)
    label("Maturation factor on GFR-scaled neonatal renal clearance (unitless)") # Table 1 footnote e Equation 2 'x 1.61'; Results 'MF of 1.61'
    # GFR_adult is a Simcyp population value the paper does not print.
    # Back-solved by the maintainers: with every other value as printed,
    # the single value that best reproduces the eight weight-normalised
    # clearances of Table S6 (GA 28, 31, 35, 38 at the end of weeks 1
    # and 2; all within +/-9%) is 6.9 L/h (115 mL/min). The value is
    # insensitive to the assumed birth weights (6.8-7.1 L/h across
    # +/-0.15 kg).
    lgfr_adult <- fixed(log(6.9))
    label("Adult reference glomerular filtration rate GFR_adult (L/h)") # NOT printed; back-solved from Table S6 (see comment above)

    # Wu 2022 neonatal GFR model as printed in Johnson 2025 Methods:
    #   GFR (L/h) = 0.086 x [ bW/1.75
    #               + (GA/34)^9.99 x 0.154 x (PNA/2)^((GA/34)^-1.71) ]
    # with bW in kg, GA in weeks and PNA in days.
    lgfr_bw <- fixed(log(0.086))
    label("Wu 2022 GFR scale per unit of (birth weight / 1.75 kg) (L/h)") # Methods GFR equation coefficient 0.086
    gfr_pna_frac <- fixed(0.154)
    label("Wu 2022 postnatal GFR increment coefficient, relative to the birth-weight term (unitless)") # Methods GFR equation coefficient 0.154
    e_ga_gfr <- fixed(9.99)
    label("Wu 2022 exponent of (GA/34 weeks) on the postnatal GFR increment (unitless)") # Methods GFR equation exponent 9.99
    e_ga_gfr_pnaexp <- fixed(-1.71)
    label("Wu 2022 exponent of (GA/34 weeks) giving the power on (PNA/2 days) (unitless)") # Methods GFR equation exponent -1.71

    # Metabolic-clearance scaling (Table 1 footnote c, Equation 1):
    #   CL_M,neonate = 1.36 x (WT/70)^0.75 x OF
    e_wt_cl_nonren <- fixed(0.75)
    label("Allometric exponent of (WT/70 kg) on metabolic clearance (unitless)") # Table 1 'allometric scaling with exponent of 0.75'; Equation 1
    # Ontogeny function OF = a2 x PMA^2 + a1 x PMA + a0 (PMA in weeks),
    # as printed in Figure 2a. Table 1 footnote c and the Results text
    # print '-0.716' and '-177.1', which give OF < -900 at term; the
    # Figure 2a form gives OF = 0.20 at 38 weeks and 0.66 at 40 weeks,
    # matching the Results text ('about 20% ... at a PMA of 38 weeks and
    # about 60% at 40 weeks') and the zero crossing near 37.5 weeks.
    of_pma2 <- fixed(-0.0716)
    label("Ontogeny-function quadratic coefficient on PMA^2 (1/week^2)") # Figure 2a 'OF = -0.0716*PMA^2 + 5.808*PMA - 117.1'
    of_pma1 <- fixed(5.808)
    label("Ontogeny-function linear coefficient on PMA (1/week)") # Figure 2a; Table 1 footnote c; Results equation
    of_int <- fixed(-117.1)
    label("Ontogeny-function intercept (unitless)") # Figure 2a 'OF = -0.0716*PMA^2 + 5.808*PMA - 117.1'

    # No residual error: PBPK simulation output.
    propSd <- fixed(0)
    label("Proportional residual error (fraction; PBPK model, none reported)")
  })

  model({
    # Age axes. PNA is supplied in months (canonical); the source
    # equations use days (GFR) and weeks (PMA).
    pna_d <- PNA * 30.4375
    pma_wk <- GA + pna_d / 7

    # Wu 2022 neonatal GFR (L/h).
    ga_rel <- GA / 34
    gfr_neo <- exp(lgfr_bw) *
      (WT_BIRTH / 1.75 + ga_rel^e_ga_gfr * gfr_pna_frac * (pna_d / 2)^(ga_rel^e_ga_gfr_pnaexp))

    # Renal clearance, Equation 2.
    cl_renal <- exp(lcl_renal) * gfr_neo / exp(lgfr_adult) * mf_renal

    # Metabolic clearance, Equation 1. The quadratic is negative below
    # about 37.5 weeks and above about 43.6 weeks PMA; the paper states
    # metabolic clearance is 'low to nonexistent' at PMA < 37.5 weeks,
    # so negative values are floored at zero.
    ontog <- of_pma2 * pma_wk^2 + of_pma1 * pma_wk + of_int
    ontog_pos <- max(ontog, 0)
    cl_nonren <- exp(lcl_nonren) * (WT / 70)^e_wt_cl_nonren * ontog_pos

    cl <- cl_renal + cl_nonren
    vc <- exp(lvc) * WT
    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
