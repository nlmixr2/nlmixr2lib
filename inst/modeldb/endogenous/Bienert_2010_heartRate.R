Bienert_2010_heartRate <- function() {
  description <- paste(
    "Population cosinor model of the 24-hour rhythm in heart rate of adult",
    "intensive care unit patients sedated with a prolonged constant propofol",
    "infusion (Bienert 2010). Heart rate oscillates about an individual mesor",
    "with a fixed 24-hour period: HR = mesor * (1 + famp * cos(2 * pi / 24 *",
    "(clock time - acrophase))). The fractional amplitude is small (4.8",
    "percent of the 90.5 beats/min mesor) and the peak falls at 14:00. IIV",
    "on the mesor only and proportional residual error. No drug term: the",
    "rhythm was estimated from the hourly vital-sign record and is not linked",
    "to the propofol concentration."
  )
  reference <- paste(
    "Bienert A, Kusza K, Wawrzyniak K, Grzeskowiak E, Kokot ZJ, Matysiak J,",
    "Grabowski T, Wolc A, Wiczling P, Regulski M (2010). Assessing circadian",
    "rhythms in propofol PK and PD during prolonged infusion in ICU patients.",
    "J Pharmacokinet Pharmacodyn 37(3):289-304. doi:10.1007/s10928-010-9161-5.",
    sep = " "
  )
  vignette <- "Bienert_2010_propofol_circadian"
  units <- list(time = "h", dosing = "n/a (endogenous rhythm; no dosing)", concentration = "beats/min")

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
        "local clock hour of each heart-rate record. The dark period of the",
        "unit was 20:00-08:00 (Figure 3 shading)."
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
      "Heart rate recorded every hour during the second and third day of",
      "infusion. Table 1: heart rate 89.8 +/- 15.7 beats/min (60-138).",
      "Figure 3 shows 10 patients (patient 3 excluded from the analysis)."
    )
  )

  ini({
    # Table 4 'Basic model with rhythm' column, heart rate (HR), printed
    # as 90.5(4.7) (1 + 0.048(36) cos(2 pi / 24 (t - 14(3.6) + T))),
    # parenthesised values being RSEs.
    lrbase <- log(90.5)
    label("Heart-rate mesor (beats/min)") # Table 4: theta HR = 90.5 (RSE 4.7%)
    rbase_circ_famp <- 0.048
    label("Fractional amplitude of the 24-hour heart-rate cosine (fraction of mesor)") # Table 4: amplitude 0.048 (RSE 36%); Results '4.8% of the messor'
    ltacro <- log(14)
    label("Acrophase, clock time of the heart-rate peak (h after midnight)") # Table 4: acrophase 14 (RSE 3.6%); Results 'maximum at 02:00 p.m.'

    # Log-normal IIV on the mesor, printed as a CV in percent (Results 'IIV
    # of 14.7%'), so omega^2 = CV^2.
    etalrbase ~ 0.021609 # Table 4: 'omega HR2 (%)' = 14.7 (RSE 34%); 0.147^2

    propSd <- 0.119
    label("Proportional residual error on heart rate (fraction)") # Table 4: 'sigma prop,HR2 (%)' = 11.9 (RSE 42%)
  })

  model({
    rbase <- exp(lrbase + etalrbase)
    tacro <- exp(ltacro)

    # Eq. 7 with t + T replaced by the record's clock time.
    hr <- rbase * (1 + rbase_circ_famp * cos(2 * pi / 24 * (TCLOCK - tacro)))

    hr ~ prop(propSd)
  })
}
