vanEsdonk_2018_pregabalin_electricalstim <- function() {
  description <- "Turnover (indirect-response) pharmacodynamic model for the electrical stimulation pain tolerance threshold after a single 300 mg oral dose of pregabalin in healthy adults, driven by a one-compartment pregabalin PK model fixed from the same paper; the drug acts as an all-or-none (on/off) fractional stimulation of the production rate kin whenever pregabalin is present, with between-occasion variability on the baseline across the placebo and pregabalin visits (van Esdonk 2018)"
  reference <- paste(
    "van Esdonk MJ, Lindeman I, Okkerse P, de Kam ML, Groeneveld GJ, Stevens J. (2018).",
    "Population pharmacokinetic/pharmacodynamic analysis of nociceptive pain models",
    "following an oral pregabalin dose administration to healthy subjects.",
    "CPT Pharmacometrics Syst Pharmacol 7(9):573-580. doi:10.1002/psp4.12318."
  )
  vignette <- "vanEsdonk_2018_pregabalin_pain_models"

  paper_specific_compartments <- c("ptt_es")

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the allometric scaling of the fixed pregabalin PK layer only (70 kg reference, exponents 0.75 on CL and 1 on Vd). No covariate effect was identified on any electrical stimulation PD parameter (van Esdonk 2018 Discussion: 'No other covariate relationships were identified in the developed PK/PD models').",
      source_name = "WT"
    ),
    OCC = list(
      description = "Study-visit occasion indicator for the between-occasion variability on the electrical stimulation baseline",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "1 = placebo visit, 2 = pregabalin visit. BOV on the baseline is load-bearing for this endpoint rather than merely descriptive: the paper reports that 'no significant effect could be estimated when no variability in the baseline between the placebo and pregabalin treatment was included', and the drug effect became estimable only after the BOV was added. Decomposed inside model() into binary indicators oc1 / oc2 that multiplex the two BOV etas onto the baseline. The placebo visit carries no dose, so Cc is zero throughout it and the on/off switch stays off; the occasion indicator itself does not gate the drug effect.",
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "pregabalin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pregabalin", units = "mg", specimen = "plasma", verified = TRUE),
    ptt_es = list(
      analyte = "electrical stimulation pain tolerance threshold",
      units = "mA",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 16L,
    n_studies = 1L,
    age_range = "19-25 years",
    age_mean = "21.75 years (SD 1.61)",
    weight_range = "54.25-77.50 kg",
    weight_mean = "68.0 kg (SD 8.22)",
    sex_female_pct = 50,
    disease_state = "Healthy volunteers",
    dose_range = "Single 300 mg oral dose, and a placebo visit",
    regions = "The Netherlands (Centre for Human Drug Research, Leiden)",
    notes = "313 electrical stimulation (single stimulus) measurements (160 placebo, 153 pregabalin treated) were available for model development. Baseline demographics from Table 1. The pain battery was performed 10 times per visit, including two predose measurements up to 1 h before dosing. Observed thresholds ranged from 5 to 50 mA."
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
    # ELECTRICAL STIMULATION PD LAYER -- estimated.
    # All values are van Esdonk 2018 Table 3, 'Electrical stimulation PTT /
    # Turnover model population parameters (RSE)' column.
    # ==================================================================
    lrbase <- log(19.1)
    label("Log of baseline electrical stimulation pain tolerance threshold (mA)")
    # Table 3 electrical stimulation: 'Baseline' = 19.1 mA [RSE 7%]

    lkout <- log(0.494)
    label("Log of first-order loss rate constant of the pain tolerance threshold (1/h)")
    # Table 3 electrical stimulation: 'k out' = 0.494/hour [RSE 24%]

    # The paper first fitted a sigmoid Emax effect on kin but estimated an
    # EC50 below 1 ug/L -- i.e. below the 20 ug/L assay lower limit of
    # quantification -- so the maximal effect was already reached at the
    # lowest measurable pregabalin concentration. It therefore replaced the
    # Emax/EC50 pair with a single all-or-none 'on/off' effect of the same
    # size and better precision. What remains is the maximal FRACTIONAL
    # increment on kin, which is the canonical `lemax`; the potency term is
    # deliberately absent rather than missing.
    lemax <- log(0.322)
    label("Log of the on/off fractional stimulation of kin while pregabalin is present (unitless)")
    # Table 3 electrical stimulation: 'Effect pregabalin' = 0.322 [RSE 18%]

    # ------------------------------------------------------------------
    # Variability. The paper states IIV was drawn 'from a ln-normal
    # distribution', so the Table 3 omega^2 values are log-scale variances;
    # they are used as printed. NOTE that the two %CV cells in the
    # electrical stimulation column do not reconcile with them under
    # CV = sqrt(exp(omega^2) - 1), which gives 36.0% for 0.122 (printed
    # 22%) and 45.3% for 0.187 (printed 70%) -- unlike the cold pressor
    # column, where all three CVs reconcile exactly. The omega^2 values are
    # taken as authoritative because each is corroborated by its own
    # bootstrap confidence interval (0.070-0.172 and 0.02-0.60); the CV
    # discrepancy is recorded in the vignette Errata.
    # No IIV on the baseline was reported for this endpoint -- only BOV.
    # ------------------------------------------------------------------
    etalemax ~ 0.187 # Table 3 electrical stimulation: 'omega 2 effect' = 0.187 (shrinkage 25%)

    # Between-occasion variability on the baseline across the two visits.
    # One shared magnitude is reported, so the second occasion's variance is
    # fixed equal to the first (the NONMEM `$OMEGA BLOCK(1) SAME` idiom, as
    # encoded by Vet_2016_midazolam.R and Chen_2023_nemonoxacin.R).
    etaiov_rbase_1 ~ 0.122 # Table 3 electrical stimulation: 'omega 2 BOV baseline' = 0.122 (shrinkage 2%)
    etaiov_rbase_2 ~ fixed(0.122) # shared BOV magnitude, second occasion

    # ------------------------------------------------------------------
    # Residual error. Table 3 reports sigma^2 = 0.0143 as a VARIANCE (the
    # parenthesised 6% is shrinkage, per the row header 'Residual error
    # (shrinkage)'); nlmixr2 takes the SD. sqrt(0.0143) = 0.119583.
    # ------------------------------------------------------------------
    propSd <- 0.119583
    label("Proportional residual error on the electrical stimulation threshold (fraction)")
    # Table 3 electrical stimulation: 'sigma 2 proportional' = 0.0143 (shrinkage 6%)
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
    # 3. Electrical stimulation turnover model (Figure 1b: a pain tolerance
    #    threshold compartment with production kin and loss kout, with
    #    pregabalin acting on kin as 'ES: on/off (+)').
    #
    #    The on/off switch is driven by the presence of measurable drug
    #    rather than by the treatment arm, which reproduces both visits from
    #    one model: on the placebo visit no dose is given, so Cc stays zero
    #    and the switch never turns on, while on the pregabalin visit the
    #    switch turns on when absorption begins after the lag time and turns
    #    off again only once the concentration returns to zero.
    #
    #    No time trend is included: 'No effect over time in the placebo
    #    occasion was significant' for this endpoint, so there is no
    #    slope_placebo counterpart to the cold pressor model.
    # ------------------------------------------------------------------
    # Table 3 reports no IIV on the baseline (only BOV) and none on kout for
    # this endpoint, so neither carries an eta; the only PD eta is on the
    # drug effect.
    rbase <- exp(lrbase + iov_rbase)
    kout <- exp(lkout)
    emax <- exp(lemax + etalemax)

    drug_on <- (Cc > 0)
    kin_es <- kout * rbase * (1 + emax * drug_on)

    d/dt(ptt_es) <- kin_es - kout * ptt_es
    ptt_es(0) <- rbase

    ptt_es ~ prop(propSd)
  })
}
