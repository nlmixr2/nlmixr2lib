Bihorel_2021_BMS986166_nalc_sd <- function() {
  description <- "Exposure-response model for the nadir absolute lymphocyte count (nALC) after a single oral dose of BMS-986166 in healthy adults (Bihorel 2021, Study IM018001). Inhibitory sigmoid Emax function of the average BMS-986166-P (active phosphorylated metabolite) blood concentration on day 1, applied on top of a fractional placebo reduction: nALC = ALC0 * (1 - Delta_placebo) * (1 - Imax * Cavg^h / (IC50^h + Cavg^h)). One observation per subject, no inter-individual variability, proportional residual error. PD-only model: Cavg,Day1 is supplied as the CAV covariate, computed from the individual empirical-Bayes predictions of the companion population PK model modellib('Bihorel_2021_BMS986166'). The repeated-dose counterpart is modellib('Bihorel_2021_BMS986166_nalc_md')."
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
      description = "Individual average BMS-986166-P blood concentration on day 1 (Cavg,Day1)",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Cavg,Day1 = AUC(0-24 h) of BMS-986166-P after the single BMS-986166 dose divided by 24 h, predicted for each subject from the individual maximum a posteriori Bayesian estimates of the final population PK model (Bihorel 2021 Methods; companion model modellib('Bihorel_2021_BMS986166')). Placebo recipients carry CAV = 0, which returns nALC = ALC0 * (1 - Delta_placebo).",
      source_name = "Cavg,Day1"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 2L,
    n_placebo = 6L,
    age_range = "19-52 years",
    disease_state = "Healthy adult volunteers",
    dose_range = "Placebo or single oral doses of BMS-986166 0.75, 2 or 5 mg",
    regions = "United States",
    notes = "31 subjects (Bihorel 2021 Results): the single-ascending dose Study IM018001 (6 of them placebo recipients) plus 1 subject from Study IM018003 who received only one BMS-986166 dose. nALC is the lowest ALC observed at any time after the first BMS-986166 or placebo dose, sampled up to 816 h after dosing (Supplementary Table S1)."
  )

  ini({
    # Bihorel 2021 Table 2, 'nALC after single dose' block, and the
    # 'Single dose' display equation in the Results section
    # 'Exposure-Response Analysis: Lymphocyte Count'. Table 2 footnote a
    # marks Imax and h as highly correlated (r^2 >= 0.8100).
    lrbase <- log(2.061); label("Baseline ALC, ALC0 (10^3 cells/uL)") # Table 2 ALC0 2.061 (RSE 6.257%)
    limax <- log(0.6747); label("Maximum fractional reduction in nALC relative to placebo, Imax (fraction)") # Table 2 Imax 0.6747 (RSE 20.59%)
    lic50 <- log(2.693); label("Cavg,Day1 at half-maximal response IC50 (ng/mL)") # Table 2 IC50 2.693 ng/mL (RSE 13.65%)
    lhill <- log(3.080); label("Hill coefficient h (unitless)") # Table 2 h 3.080 (RSE 76.98%)
    placebo_effect <- 0.1992; label("Fractional reduction in nALC from ALC0 in placebo recipients, Delta_placebo (fraction)") # Table 2 Delta_placebo 0.1992 (RSE 29.20%)

    # No IIV was estimated in the E-R models (Methods). Table 2 gives the
    # residual variance 0.08569 with 29.27 %CV = 100 * sqrt(0.08569).
    propSd <- 0.2927; label("Proportional residual error on nALC (fraction)") # Table 2 residual variability 0.08569 (29.27 %CV) -> sqrt(0.08569)
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
