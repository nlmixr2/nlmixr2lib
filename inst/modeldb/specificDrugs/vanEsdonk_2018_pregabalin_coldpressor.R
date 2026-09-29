vanEsdonk_2018_pregabalin_coldpressor <- function() {
  description <- "Turnover (indirect-response) pharmacodynamic model for the cold pressor pain tolerance threshold after a single 300 mg oral dose of pregabalin in healthy adults, driven by a one-compartment pregabalin PK model fixed from the same paper; a linear decrease of the threshold over time is combined with a linear stimulatory pregabalin concentration effect on the production rate kin, with between-occasion variability on the baseline across the placebo and pregabalin visits (van Esdonk 2018)"
  reference <- paste(
    "van Esdonk MJ, Lindeman I, Okkerse P, de Kam ML, Groeneveld GJ, Stevens J. (2018).",
    "Population pharmacokinetic/pharmacodynamic analysis of nociceptive pain models",
    "following an oral pregabalin dose administration to healthy subjects.",
    "CPT Pharmacometrics Syst Pharmacol 7(9):573-580. doi:10.1002/psp4.12318."
  )
  vignette <- "vanEsdonk_2018_pregabalin_pain_models"

  paper_specific_compartments <- c("ptt_cp")

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the allometric scaling of the fixed pregabalin PK layer only (70 kg reference, exponents 0.75 on CL and 1 on Vd). No covariate effect was identified on any cold pressor PD parameter (van Esdonk 2018 Discussion: 'No other covariate relationships were identified in the developed PK/PD models').",
      source_name = "WT"
    ),
    OCC = list(
      description = "Study-visit occasion indicator for the between-occasion variability on the cold pressor baseline",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "1 = placebo visit, 2 = pregabalin visit. The paper estimated between-occasion variability (BOV) between these two visits within an individual because 'the baseline response to a pain model between two visits of the same individual could differ significantly'. Decomposed inside model() into binary indicators oc1 / oc2 that multiplex the two BOV etas onto the baseline. The placebo visit carries no dose, so Cc is zero throughout it and the drug term vanishes; the occasion indicator itself does not gate the drug effect.",
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "pregabalin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pregabalin", units = "mg", specimen = "plasma", verified = TRUE),
    ptt_cp = list(
      analyte = "cold pressor pain tolerance threshold",
      units = "s",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 15L,
    n_studies = 1L,
    age_range = "19-25 years",
    age_mean = "21.75 years (SD 1.61)",
    weight_range = "54.25-77.50 kg",
    weight_mean = "68.0 kg (SD 8.22)",
    sex_female_pct = 50,
    disease_state = "Healthy volunteers",
    dose_range = "Single 300 mg oral dose, and a placebo visit",
    regions = "The Netherlands (Centre for Human Drug Research, Leiden)",
    notes = "291 cold pressor measurements (148 placebo, 143 pregabalin) were available for model development. One of the 16 subjects who received pregabalin was excluded from the cold pressor analysis because of a continuously maximal pain tolerance threshold of 120 s on both the placebo and the pregabalin visit, so n_subjects is 15 for this endpoint while the parent PK model (modellib('vanEsdonk_2018_pregabalin')) used all 16. Baseline demographics from Table 1. The pain battery was performed 10 times per visit, including two predose measurements up to 1 h before dosing."
  )

  ini({
    # ==================================================================
    # PREGABALIN PK LAYER -- fixed, not re-estimated.
    # The paper used a SEQUENTIAL modelling approach: 'individual post hoc
    # Bayesian estimates of the developed PK model were added to the PD
    # dataset'. The PK layer is therefore a forcing function here and every
    # PK parameter (and its IIV) is wrapped in fixed(); only the PD
    # parameters below were estimated in this run. Values are van Esdonk
    # 2018 Table 2; see modellib('vanEsdonk_2018_pregabalin') for the
    # standalone PK model those values came from.
    # ==================================================================
    lka <- fixed(log(6.07))
    label("Log of apparent first-order absorption rate constant (1/h)")
    # Table 2: 'k a (/hour)' = 6.07

    ltlag <- fixed(log(0.495))
    label("Log of absorption lag time (h)")
    # Table 2: 'Lag time (hour)' = 0.495

    lvc <- fixed(log(31.1))
    label("Log of apparent central volume of distribution Vd/F at 70 kg (L)")
    # Table 2: 'V d /70 kg (L)' = 31.1

    lcl <- fixed(log(4.5))
    label("Log of apparent clearance CL/F at 70 kg (L/h)")
    # Table 2: 'CL/70 kg (L/hour)' = 4.5

    e_wt_cl <- fixed(0.75)
    label("Allometric body-weight exponent on apparent clearance (unitless)")
    # Covariate analysis: 'clearance (CL; exponent = 0.75)'

    e_wt_vc <- fixed(1)
    label("Allometric body-weight exponent on apparent volume of distribution (unitless)")
    # Covariate analysis: 'volume of distribution (V d ; exponent = 1)'

    etalka ~ fixed(2.6) # Table 2: 'omega 2 k a' = 2.6 (CV 353%)
    etaltlag ~ fixed(7.09e-5) # Table 2: 'omega 2 lag time' = 7.09E-5 (CV 0.842%)
    etalvc ~ fixed(0.0101) # Table 2: 'omega 2 V d / F' = 0.0101 (CV 10.1%)
    etalcl ~ fixed(0.00672) # Table 2: 'omega 2 CL/ F' = 0.00672 (CV 8.21%)

    # ==================================================================
    # COLD PRESSOR PD LAYER -- estimated.
    # All values are van Esdonk 2018 Table 3, 'Cold pressor PTT / Turnover
    # model population parameters (RSE)' column.
    # ==================================================================
    lrbase <- log(16.9)
    label("Log of baseline cold pressor pain tolerance threshold (s)")
    # Table 3 cold pressor: 'Baseline' = 16.9 seconds [RSE 16.8%]

    lkout <- log(0.39)
    label("Log of first-order loss rate constant of the pain tolerance threshold (1/h)")
    # Table 3 cold pressor: 'k out' = 0.39/hour [RSE 21%]

    # Signed, so NOT log-transformed (a negative typical value forbids a log
    # transform, per the slope_placebo / slope_drug register entry).
    slope_placebo <- -0.07
    label("Linear drift of the cold pressor threshold over time, absent drug (s/h)")
    # Table 3 cold pressor: 'Slope over time' = -0.07 seconds/hour [RSE 57.8%]

    slope_drug <- 0.135
    label("Linear fractional stimulation of kin per unit pregabalin concentration (1/(mg/L))")
    # Table 3 cold pressor: 'Slope pregabalin' = 0.135 1/mg/L [RSE 16.7%]

    # ------------------------------------------------------------------
    # Variability. Table 3 reports omega^2 with a %CV column; the cold
    # pressor column is internally consistent with log-normal etas under
    # CV = sqrt(exp(omega^2) - 1): 0.057 -> 24.2%, 0.283 -> 57.2%,
    # 0.738 -> 105%, matching all three printed CVs. The omega^2 values
    # are therefore log-scale variances.
    # ------------------------------------------------------------------
    etalrbase ~ 0.283 # Table 3 cold pressor: 'omega 2 baseline' = 0.283 (CV 57.2%, shrinkage 2%)
    etalkout ~ 0.738 # Table 3 cold pressor: 'omega 2 k out' = 0.738 (CV 105%, shrinkage 19%)

    # Between-occasion variability on the baseline across the two visits.
    # One shared magnitude is reported, so the second occasion's variance is
    # fixed equal to the first (the NONMEM `$OMEGA BLOCK(1) SAME` idiom, as
    # encoded by Vet_2016_midazolam.R and Chen_2023_nemonoxacin.R).
    etaiov_rbase_1 ~ 0.057 # Table 3 cold pressor: 'omega 2 BOV baseline' = 0.057 (CV 24.2%, shrinkage 27%)
    etaiov_rbase_2 ~ fixed(0.057) # shared BOV magnitude, second occasion

    # ------------------------------------------------------------------
    # Residual error. Table 3 reports sigma^2 = 0.041 as a VARIANCE (the
    # parenthesised 6% is shrinkage, per the row header 'Residual error
    # (shrinkage)'); nlmixr2 takes the SD. sqrt(0.041) = 0.202485.
    # ------------------------------------------------------------------
    propSd <- 0.202485
    label("Proportional residual error on the cold pressor threshold (fraction)")
    # Table 3 cold pressor: 'sigma 2 proportional' = 0.041 (shrinkage 6%)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Pregabalin PK forcing function (Figure 1a). Identical structure to
    #    modellib('vanEsdonk_2018_pregabalin'); carried here with its IIV so
    #    that each simulated subject gets an individual concentration-time
    #    profile, reproducing the paper's use of individual post hoc
    #    Bayesian PK estimates. Cc carries no residual error because it is a
    #    driver here, not an observed endpoint of this run.
    # ------------------------------------------------------------------
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    Cc <- central / vc

    # ------------------------------------------------------------------
    # 2. Between-occasion variability on the baseline. OCC = 1 is the
    #    placebo visit, OCC = 2 the pregabalin visit.
    # ------------------------------------------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_rbase <- oc1 * etaiov_rbase_1 + oc2 * etaiov_rbase_2

    # ------------------------------------------------------------------
    # 3. Cold pressor turnover model (Figure 1b: a pain tolerance threshold
    #    compartment with production kin and loss kout, with pregabalin
    #    acting on kin as 'CP: Linear (+)').
    #
    #    The threshold declines linearly over time in the absence of drug
    #    ('A linear decrease in the cold pressor PTT over time during the
    #    placebo occasion gave a significant improvement'), and pregabalin
    #    raises kin proportionally ('The drug effect was implemented in the
    #    structural model as ... a proportional effect on the k in'; 'A
    #    linear relationship between the pregabalin concentrations and the
    #    k in showed to be superior').
    #
    #    The time drift is carried on the TARGET level rather than added to
    #    kin directly, because Table 3 reports the slope in seconds/hour --
    #    the units of the threshold per unit time, not of kin per unit time.
    #    Writing kin = kout * (rbase + slope_placebo * t) makes the system
    #    asymptotically linear in t with exactly the reported slope, which is
    #    the 'linear decrease over time' the paper describes; adding the
    #    slope to kin instead would produce a constant offset, not a drift.
    #    The drug term is dimensionless (slope_drug has units 1/(mg/L)),
    #    confirming it multiplies rather than adds.
    #
    #    Time is measured from the dose, matching the paper's 'Time after
    #    dose (hours)' axis; the two predose measurements sit at negative
    #    times, where Cc is zero.
    # ------------------------------------------------------------------
    rbase <- exp(lrbase + etalrbase + iov_rbase)
    kout <- exp(lkout + etalkout)

    ptt_cp_target <- rbase + slope_placebo * time
    kin_cp <- kout * ptt_cp_target * (1 + slope_drug * Cc)

    d/dt(ptt_cp) <- kin_cp - kout * ptt_cp
    ptt_cp(0) <- rbase

    ptt_cp ~ prop(propSd)
  })
}
