Bihorel_2021_BMS986166_nddhr <- function() {
  description <- "Exposure-response model for the day-1 nadir of time-matched, placebo-corrected heart rate (nDDHR) after the first oral dose of BMS-986166 in healthy adults (Bihorel 2021). Inhibitory sigmoid Emax function of the average BMS-986166-P (active phosphorylated metabolite) blood concentration on day 1: nDDHR = nDDHRplacebo + (MaxDelta - nDDHRplacebo) * Cavg^h / (IC50^h + Cavg^h), with both nDDHRplacebo (-9.08 bpm) and MaxDelta (-19.7 bpm) negative. One observation per subject, no inter-individual variability, proportional residual error. PD-only model: Cavg,Day1 is supplied as the CAV covariate, computed from the individual empirical-Bayes predictions of the companion population PK model modellib('Bihorel_2021_BMS986166')."
  reference <- "Bihorel S, Singhal S, Shevell D, Sun H, Xie J, Basdeo S, Liu A, Dutta S, Ludwig E, Huang H, Lin K, Fura A, Throup J, Girgis IG. Population Pharmacokinetic Analysis of BMS-986166, a Novel Selective Sphingosine-1-Phosphate-1 Receptor Modulator, and Exposure-Response Assessment of Lymphocyte Counts and Heart Rate in Healthy Participants. Clin Pharmacol Drug Dev. 2021;10(1):8-21. doi:10.1002/cpdd.878. PMCID: PMC7821288."
  vignette <- "Bihorel_2021_BMS986166"
  units <- list(time = "h", dosing = "n/a (PD-only model; no dose events)", concentration = "ng/mL", response = "bpm")

  covariateData <- list(
    CAV = list(
      description = "Individual average BMS-986166-P blood concentration on day 1 (Cavg,Day1)",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Cavg,Day1 = AUC(0-24 h) of BMS-986166-P after the first BMS-986166 dose divided by 24 h, predicted for each subject from the individual maximum a posteriori Bayesian estimates of the final population PK model (Bihorel 2021 Methods, Pharmacokinetic and Exposure-Response Model Development; companion model modellib('Bihorel_2021_BMS986166')). Model-predicted average rather than peak concentrations were used because the PK model underpredicts peaks while reproducing average exposure (Discussion). Subjects randomized to placebo carry CAV = 0, which returns nDDHR = nDDHRplacebo.",
      source_name = "Cavg,Day1"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 61L,
    n_studies = 2L,
    n_placebo = 13L,
    age_range = "19-52 years",
    weight_range = "58.3-104.1 kg",
    sex_female_pct = 3.2,
    disease_state = "Healthy adult volunteers",
    dose_range = "Placebo or BMS-986166 0.25-5 mg (first dose of the single-ascending dose Study IM018001 or of the multiple-ascending dose Study IM018003)",
    regions = "United States",
    notes = "61 nDDHR values (one per subject, 13 of them from placebo recipients) pooled from Studies IM018001 and IM018003 (Bihorel 2021 Results). One placebo subject with an nDDHR of -40.6 bpm was excluded as an outlier. Demographics for the 62 enrolled participants: 96.8% male, mean age 34.7 years, mean weight 84.2 kg. Heart rate came from continuous cardiac monitoring summarised as hourly averages; DDHR is the hourly average on day 1 minus the time-matched hourly average on day -1, when every subject received placebo."
  )

  ini({
    # Bihorel 2021 Table 2, nDDHR block, and the display equation in the
    # Results section 'Exposure-Response Analysis: HR'. Both nDDHRplacebo
    # and MaxDelta are negative (Results), so they are kept on the natural
    # scale rather than log-transformed.
    rbase <- -9.08; label("nDDHR in placebo recipients, nDDHRplacebo (bpm)") # Table 2 nDDHRplacebo -9.08 bpm (RSE 10.8%)
    nadir_max <- -19.7; label("nDDHR at maximal drug effect, MaxDelta (bpm)") # Table 2 MaxDelta -19.7 bpm (RSE 10.9%)
    lic50 <- log(1.26); label("Cavg,Day1 at half-maximal response IC50 (ng/mL)") # Table 2 IC50 1.26 ng/mL (RSE 40.1%)
    lhill <- log(1.84); label("Hill coefficient h (unitless)") # Table 2 h 1.84 (RSE 33.2%)

    # No IIV was estimated in the E-R models (Methods). Table 2 gives the
    # residual variance 0.167 with 40.9 %CV = 100 * sqrt(0.167): a
    # constant-CV (proportional) error.
    propSd <- 0.4087; label("Proportional residual error on nDDHR (fraction)") # Table 2 residual variability 0.167 (40.9 %CV) -> sqrt(0.167)
  })

  model({
    ic50 <- exp(lic50)
    hill <- exp(lhill)

    # Inhibitory sigmoid Emax (Bihorel 2021 Results display equation)
    nadir_hr_change <- rbase + (nadir_max - rbase) * CAV^hill / (ic50^hill + CAV^hill)

    nadir_hr_change ~ prop(propSd)
  })
}
