Shin_2017_levetiracetam <- function() {
  description <- "One-compartment population PK model for levetiracetam in Korean neonates with seizures (Shin 2017). Structural parameters (V, CL) reported on a per-kg-body-weight basis (linear scaling by body weight). Drug absorption was not modelled because trough-style sampling between 6 and 23 hours after dose did not capture the absorption phase; intravenous and oral doses are therefore modelled as bolus inputs directly into the central compartment with bioavailability fixed at 1."
  reference   <- "Shin JW, Jung YS, Park K, Lee SM, Eun HS, Park MS, Park KI, Namgung R. Experience and pharmacokinetics of Levetiracetam in Korean neonates with neonatal seizures. Korean J Pediatr. 2017 Feb;60(2):50-54. doi:10.3345/kjp.2017.60.2.50"
  vignette    <- "Shin_2017_levetiracetam"
  units       <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at the time of medication",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Shin 2017 reports V and CL on a per-kg basis (1.15 L/kg and 0.083 L/hr/kg respectively); the model applies a linear weight scaling (allometric exponent fixed to 1.0) so that vc = (V per kg) * WT and cl = (CL per kg) * WT. Reference weight is therefore 1 kg; the paper's typical individual is recovered at the cohort median WT = 4.3 kg (Results paragraph 2).",
      source_name        = "WT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 18,
    n_studies       = 1,
    n_observations  = 151,
    age_range       = "Postmenstrual age 22.3-66.0 weeks at medication; postnatal age 3.4-35.1 weeks (Shin 2017 Table 2). Gestational age at birth 24.3-39.9 weeks (Table 1).",
    age_median      = "Postmenstrual age 48.7 weeks (Shin 2017 Results paragraph 2); mean postmenstrual age 47.6 +/- 6.6 weeks (Table 2).",
    weight_range    = "0.535-10.45 kg at medication (Shin 2017 Table 2).",
    weight_median   = "4.3 kg at medication (Shin 2017 Results paragraph 2); mean 4.27 +/- 1.76 kg (Table 2). Birth weight: mean 2.39 +/- 1.28 kg (range 0.535-5.44 kg, Table 1).",
    sex_female_pct  = 39,
    race_ethnicity  = c(Asian_Korean = 100),
    disease_state   = "Neonates with electro-clinical or electrographic-only confirmed seizures. Etiologies: hypoxic-ischemic encephalopathy 67% (12/18), brain malformation 22% (4/18), meningoencephalitis 6% (1/18), intracerebral hemorrhage 6% (1/18) (Shin 2017 Table 1).",
    dose_range      = "Loading dose 4.9-59.5 mg/kg (mean 20.0 +/- 16.1); maintenance dose 4.5-99.5 mg/kg/day (mean 29.0 +/- 20.1). Mix of intravenous and oral administration; route switched per clinical status (Shin 2017 Table 2 and Methods section 2).",
    regions         = "South Korea (two neonatal intensive care units of Yonsei University College of Medicine, Seoul; June 2013 to June 2015).",
    co_medication   = "94% of patients received concomitant antiepileptic drugs in addition to levetiracetam: phenobarbital + LEV 61.1%, phenobarbital + phenytoin + LEV 33.3%, LEV monotherapy 5.5% (Shin 2017 Results paragraph 1). No statistically significant association between CL or V and concomitant AED use (Discussion paragraph 6).",
    notes           = "Retrospective single-centre population PK analysis. Population PK parameters in Shin 2017 Table 5. Therapeutic range cited as 6-20 ug/mL; observed concentrations 1-41 ug/mL with 55% of samples in range (Table 4)."
  )

  ini({
    # Structural parameters - final-model estimates from Shin 2017 Table 5.
    # Values are reported on a per-kg basis so that the typical individual at
    # the cohort median weight 4.3 kg has V = 1.15 * 4.3 = 4.945 L and
    # CL = 0.083 * 4.3 = 0.357 L/hr (matching Results paragraph 2).
    lvc <- log(1.15);  label("Central volume per kg body weight (V, L/kg)")        # Shin 2017 Table 5: V = 1.15 L/kg (RSE 29.7%)
    lcl <- log(0.083); label("Clearance per kg body weight (CL, L/hr/kg)")         # Shin 2017 Table 5: CL = 0.083 L/hr/kg (RSE 12.7%)

    # Inter-individual variability. Shin 2017 reports "exponential error model
    # for random inter-individual variability ... mean zero and variance
    # omega^2" (Methods section 4); log-normal IIV converts as
    # omega^2 = log(CV^2 + 1).
    etalvc ~ 0.4158                                                                # Shin 2017 Table 5: omega_V CV = 71.8% -> log(0.718^2 + 1) = 0.4158
    etalcl ~ 0.0884                                                                # Shin 2017 Table 5: omega_CL CV = 30.4% -> log(0.304^2 + 1) = 0.0884

    # Residual error. Shin 2017 does not report the residual error structure
    # or magnitude (Table 5 lists only V, CL, and the two IIV variances). A
    # proportional residual error with token initial value 0.20 (20% CV) is
    # used here as a reasonable starting point for re-fits; the value is an
    # extraction-time assumption, not a Shin 2017 estimate. See vignette
    # "Assumptions and deviations".
    propSd <- 0.20;    label("Proportional residual error (fraction)")             # Not reported by Shin 2017; assumed proportional with token initial value; see vignette
  })

  model({
    # Individual parameters. Per-kg structural parameters are linearly scaled
    # by weight (allometric exponent 1.0; reference weight 1 kg) so that
    # vc and cl recover the paper's typical individual values at the cohort
    # median WT = 4.3 kg.
    vc <- exp(lvc + etalvc) * WT
    cl <- exp(lcl + etalcl) * WT

    kel <- cl / vc

    # One-compartment model. Both intravenous and oral doses entered the
    # central compartment directly; Shin 2017 explicitly did not model the
    # absorption phase because sampling occurred 6-23 hours after dose
    # (Results paragraph 2 / Methods section 3). Bioavailability is fixed
    # at 1 for both routes.
    d/dt(central) <- -kel * central

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
