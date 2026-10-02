Bihorel_2021_BMS986166_nalc_md <- function() {
  description <- "Exposure-response model for the nadir absolute lymphocyte count (nALC) after 28 days of once-daily oral BMS-986166 in healthy adults (Bihorel 2021, Study IM018003). Inhibitory sigmoid Emax function of the average BMS-986166-P (active phosphorylated metabolite) blood concentration on day 28, applied on top of a fractional placebo reduction: nALC = ALC0 * (1 - Delta_placebo) * (1 - Imax * Cavg^h / (IC50^h + Cavg^h)). One observation per subject, no inter-individual variability, proportional residual error. PD-only model: Cavg,Day28 is supplied as the CAV covariate, computed from the individual empirical-Bayes predictions of the companion population PK model modellib('Bihorel_2021_BMS986166'). The single-dose counterpart is modellib('Bihorel_2021_BMS986166_nalc_sd')."
  reference <- "Bihorel S, Singhal S, Shevell D, Sun H, Xie J, Basdeo S, Liu A, Dutta S, Ludwig E, Huang H, Lin K, Fura A, Throup J, Girgis IG. Population Pharmacokinetic Analysis of BMS-986166, a Novel Selective Sphingosine-1-Phosphate-1 Receptor Modulator, and Exposure-Response Assessment of Lymphocyte Counts and Heart Rate in Healthy Participants. Clin Pharmacol Drug Dev. 2021;10(1):8-21. doi:10.1002/cpdd.878. PMCID: PMC7821288."
  vignette <- "Bihorel_2021_BMS986166"
  units <- list(
    time = "h",
    dosing = "n/a (PD-only model; no dose events)",
    concentration = "ng/mL",
    response = "10^3 cells/uL"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average BMS-986166-P blood concentration on day 28 (Cavg,Day28)",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Cavg,Day28 = AUC over the 28th once-daily dosing interval of BMS-986166-P divided by 24 h, predicted for each subject from the individual maximum a posteriori Bayesian estimates of the final population PK model (Bihorel 2021 Methods; companion model modellib('Bihorel_2021_BMS986166')). This is an average over the LAST dosing interval of a 28-day regimen, not a steady-state value: the parent half-life is roughly 10 days, so day 28 is about 85% of steady state. Placebo recipients carry CAV = 0, which returns nALC = ALC0 * (1 - Delta_placebo).",
      source_name = "Cavg,Day28"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 1L,
    n_placebo = 8L,
    age_range = "19-52 years",
    disease_state = "Healthy adult volunteers",
    dose_range = "Placebo or BMS-986166 0.25, 0.75 or 1.5 mg once daily for 28 days",
    regions = "United States",
    notes = "31 subjects of the multiple-ascending dose Study IM018003, 8 of them placebo recipients (Bihorel 2021 Results). nALC is the lowest ALC observed at any time after the first BMS-986166 or placebo dose; ALC was sampled through 1176 h after the 28th dose (Supplementary Table S1)."
  )

  ini({
    # Bihorel 2021 Table 2, 'nALC after repeated dosing' block, and the
    # 'Repeated dosing' display equation in the Results section
    # 'Exposure-Response Analysis: Lymphocyte Count'.
    lrbase <- log(1.95); label("Baseline ALC, ALC0 (10^3 cells/uL)") # Table 2 ALC0 1.95 (RSE 5.55%)
    limax <- log(0.768); label("Maximum fractional reduction in nALC relative to placebo, Imax (fraction)") # Table 2 Imax 0.768 (RSE 9.26%)
    lic50 <- log(1.72); label("Cavg,Day28 at half-maximal response IC50 (ng/mL)") # Table 2 IC50 1.72 ng/mL (RSE 23.9%)
    lhill <- log(1.65); label("Hill coefficient h (unitless)") # Table 2 h 1.65 (RSE 44.8%)
    placebo_effect <- 0.284; label("Fractional reduction in nALC from ALC0 in placebo recipients, Delta_placebo (fraction)") # Table 2 Delta_placebo 0.284 (RSE 29.4%)

    # No IIV was estimated in the E-R models (Methods). Table 2 gives the
    # residual variance 0.0880 with 29.7 %CV = 100 * sqrt(0.0880).
    propSd <- 0.2966; label("Proportional residual error on nALC (fraction)") # Table 2 residual variability 0.0880 (29.7 %CV) -> sqrt(0.0880)
  })

  model({
    rbase <- exp(lrbase)
    imax <- exp(limax)
    ic50 <- exp(lic50)
    hill <- exp(lhill)

    nadir_lymphocyte_count <- rbase * (1 - placebo_effect) * (1 - imax * CAV^hill / (ic50^hill + CAV^hill))

    nadir_lymphocyte_count ~ prop(propSd)
  })
}
