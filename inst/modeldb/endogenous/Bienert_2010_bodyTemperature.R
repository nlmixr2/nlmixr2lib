Bienert_2010_bodyTemperature <- function() {
  description <- paste(
    "Population cosinor model of the 24-hour rhythm in body temperature of",
    "adult intensive care unit patients sedated with a prolonged constant",
    "propofol infusion (Bienert 2010). Body temperature oscillates about an",
    "individual mesor with a fixed 24-hour period: T = mesor * (1 + famp *",
    "cos(2 * pi / 24 * (clock time - acrophase))). The fractional amplitude",
    "is very small (0.326 percent of the 37.0 degC mesor, about 0.12 degC)",
    "and the peak falls at 21:24. IIV on the mesor only and proportional",
    "residual error. No drug term: the rhythm was estimated from the hourly",
    "vital-sign record and is not linked to the propofol concentration."
  )
  reference <- paste(
    "Bienert A, Kusza K, Wawrzyniak K, Grzeskowiak E, Kokot ZJ, Matysiak J,",
    "Grabowski T, Wolc A, Wiczling P, Regulski M (2010). Assessing circadian",
    "rhythms in propofol PK and PD during prolonged infusion in ICU patients.",
    "J Pharmacokinet Pharmacodyn 37(3):289-304. doi:10.1007/s10928-010-9161-5.",
    sep = " "
  )
  vignette <- "Bienert_2010_propofol_circadian"
  units <- list(time = "h", dosing = "n/a (endogenous rhythm; no dosing)", concentration = "degC")

  covariateData <- list(
    TCLOCK = list(
      description = "Wall-clock time of day at each record, in decimal hours (0-24)",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying per record. Bienert 2010 Eq. 7 writes the cosine",
        "argument as 2 * pi / 24 * (t - acrophase + T), with t the time since",
        "the start of infusion and T the clock time at the start of infusion,",
        "so t + T is the clock time of the record; because the period is",
        "24 h, (t + T) modulo 24 = TCLOCK gives the same value. Supply the",
        "local clock hour of each temperature record. The dark period of the",
        "unit was 20:00-08:00 (Figure 4 shading)."
      ),
      source_name = "t + T (Eq. 7; time since infusion start plus clock time at infusion start)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 10L,
    n_studies = 1L,
    age_range = "53-85 years",
    age_median = "63 years (mean, SD 11)",
    weight_range = "50-95 kg",
    weight_median = "72.6 kg (mean, SD 22.8)",
    sex_female_pct = 27.3,
    disease_state = paste(
      "Critically ill adult intensive care unit patients sedated with a",
      "constant propofol infusion for more than 24 h; clinically stable and",
      "unstimulated during recording; no vasopressors; neuromuscular",
      "blockade, head trauma and neurological disorders excluded."
    ),
    dose_range = "Propofol 1.5 mg/kg then 0.7-0.9 mg/kg/h (sedation, not a model input).",
    regions = "Poland (Poznan)",
    notes = paste(
      "Body temperature recorded every hour during the second and third day",
      "of infusion. Table 1: body temperature 37.0 +/- 0.5 degC (35.7-38.2).",
      "Figure 4 shows 10 patients (patient 3 excluded from the analysis)."
    )
  )

  ini({
    # Table 4 'Basic model with rhythm' column, body temperature (TEMP),
    # printed as 37(0.3) (1 + 0.00326(52) cos(2 pi / 24 (t - 21.4(9.5) + T))),
    # parenthesised values being RSEs.
    lrbase <- log(37)
    label("Body-temperature mesor (degC)") # Table 4: theta TEMP = 37 (RSE 0.3%)
    rbase_circ_famp <- 0.00326
    label("Fractional amplitude of the 24-hour body-temperature cosine (fraction of mesor)") # Table 4: amplitude 0.00326 (RSE 52%); Results '0.32% of the mesor'
    ltacro <- log(21.4)
    label("Acrophase, clock time of the body-temperature peak (h after midnight)") # Table 4: acrophase 21.4 (RSE 9.5%); Results '21.4 h (21:24 p.m.)'

    # Log-normal IIV on the mesor, printed as a CV in percent (Results 'low
    # IIV of 0.82%'), so omega^2 = CV^2.
    etalrbase ~ 6.77329e-05 # Table 4: 'omega TEMP2 (%)' = 0.823 (RSE 38%); 0.00823^2

    # Table 4 prints 0.0996 for the basic model and 0.966 for the model with
    # rhythm. The with-rhythm value is the one carried; see the vignette for
    # why 0.966 percent (about 0.36 degC) is the consistent reading.
    propSd <- 0.00966
    label("Proportional residual error on body temperature (fraction)") # Table 4: 'sigma prop,TEMP2 (%)' = 0.966 (RSE 19%)
  })

  model({
    rbase <- exp(lrbase + etalrbase)
    tacro <- exp(ltacro)

    # Eq. 7 with t + T replaced by the record's clock time.
    body_temperature <- rbase * (1 + rbase_circ_famp * cos(2 * pi / 24 * (TCLOCK - tacro)))

    body_temperature ~ prop(propSd)
  })
}
