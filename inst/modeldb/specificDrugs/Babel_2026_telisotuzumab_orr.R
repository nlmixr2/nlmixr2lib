Babel_2026_telisotuzumab_orr <- function() {
  description <- paste0(
    "Logistic-regression exposure-efficacy model for overall response ",
    "rate (ORR) by Independent Central Review in adults with c-Met ",
    "protein overexpressing, EGFR wild-type, non-squamous non-small ",
    "cell lung cancer treated with telisotuzumab vedotin 1.6 or ",
    "1.9 mg/kg every 2 weeks (Babel 2026, n = 193, LUMINOSITY phase 2 ",
    "study). The probability of response is ",
    "expit(a + b * log(CAV)), where CAV is the individual average ",
    "Teliso-V CONJUGATE serum concentration up to the time of event ",
    "or end of treatment, in ug/mL. There is no PK layer and no ODE: ",
    "exposure is supplied as data, and the source obtained it from ",
    "post hoc estimates of the companion population PK model ",
    "Babel_2026_telisotuzumab using the actual doses each patient ",
    "received. No covariate was retained: Babel 2026 screened ",
    "demographics, c-Met expression level, prior therapy, ECOG status ",
    "and ADA status by stepwise forward selection and backward ",
    "elimination and reports that 'no covariates were found to have a ",
    "significant effect on efficacy or safety'. The payload was NOT ",
    "predictive of efficacy - Babel 2026 reports flat ORR ",
    "relationships against MMAE exposure - which is why this model is ",
    "driven by conjugate and not by payload exposure. NOTE: Babel ",
    "2026 does not tabulate the regression coefficients anywhere, so ",
    "the two values here were recovered by digitising the fitted ",
    "curve of Figure 2; see the vignette for the digitisation, its ",
    "residuals and the Table 1 cross-check."
  )
  reference <- paste(
    "Babel H, Brunsdon P, Engelhardt B, Schmitt V, Ratajczak C, Mensing S,",
    "Menon RM, Parikh A. Population pharmacokinetics and exposure-response",
    "analyses for telisotuzumab vedotin in patients with c-Met protein",
    "overexpressing tumors.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70219.",
    "doi:10.1002/psp4.70219. PMCID PMC12945708.",
    "The regression coefficients are not tabulated by the source; they were",
    "digitised from the fitted line of Figure 2 and validated against the",
    "simulated response probabilities in Table 1.",
    sep = " "
  )
  vignette <- "Babel_2026_telisotuzumab"

  units <- list(
    time          = "n/a (static landmark exposure-response model; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_orr_central (probability of overall response by Independent Central Review, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Individual average serum concentration of the telisotuzumab vedotin CONJUGATE, computed up to the time of the response assessment or up to the end of treatment. Supplied as data: this model has no PK layer.",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Babel 2026 Methods, Exposure-Response Analyses: 'Individual exposure metrics for Teliso-V conjugate and unconjugated MMAE payload were estimated based on post hoc estimates from the population PK models, utilizing actual doses received by patients', with Cavg defined 'up to the time of the event or up to the end of treatment if no event occurred'. Reproduce the column with the companion popPK model modellib('Babel_2026_telisotuzumab'). Because actual rather than planned doses are used, this Cavg is systematically LOWER than the nominal steady-state Dose/(CL*tau): the median implied by the Table 1 simulation is about 6.2 ug/mL at 1.9 mg/kg Q2W against a nominal 7.0 ug/mL at the median 68.9 kg weight. Babel 2026 Figure 2 shows the analysis-set range spanning roughly 2.2 to 10.8 ug/mL, with binned quartile medians near 3.6, 5.1, 6.6 and 8.5 ug/mL. Enters on the NATURAL-LOG scale and UNCENTRED. The conjugate, not the payload, is the efficacy driver.",
      source_name        = "CavgADC"
    )
  )

  covariatesDataExcluded <- list(
    CMET_HIGH = list(
      description = "c-Met protein overexpression level indicator; 1 = high (at least 50% of tumour cells with strong 3+ immunohistochemistry staining), 0 = intermediate (at least 25% to below 50%).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Listed in Babel 2026 Table S2 as a covariate of interest for the exposure-efficacy regression, but not retained: 'No covariates were found to have a significant effect on efficacy or safety'. No point estimate exists on disk. Documented here to preserve the covariate screen without carrying a convention warning; the same applies to the remaining Table S2 efficacy covariates (age, sex, race, ethnicity, body weight, prior therapy, number of prior systemic therapies, baseline ECOG status, treatment-emergent ADA status and nAb status)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 193L,
    n_studies      = 1L,
    n_observations = "193 binary response records (one per patient; landmark analysis, no repeated measures)",
    disease_state  = "c-Met protein overexpressing, EGFR wild-type, non-squamous non-small cell lung cancer, previously treated; c-Met high defined as at least 50% of tumour cells with strong 3+ immunohistochemistry staining and c-Met intermediate as at least 25% to below 50%",
    dose_range     = "telisotuzumab vedotin 1.6 mg/kg (n = 25) or 1.9 mg/kg (n = 168) every 2 weeks as an intravenous infusion",
    regions        = "LUMINOSITY (NCT03539536) international sites",
    notes          = paste0(
      "Baseline demographics of this analysis set are in Babel 2026 ",
      "Table S5. Responses were adjudicated by Independent Central ",
      "Review; the phase 1 study was deliberately excluded from the ",
      "exposure-efficacy analysis because its c-Met overexpression ",
      "criteria differed from LUMINOSITY. The observed ORR in ",
      "LUMINOSITY was 29% overall and 35% in the c-Met high subgroup. ",
      "Babel 2026 notes that the small 1.6 mg/kg group and the ",
      "non-randomised assignment of dose groups leave the analysis ",
      "open to bias and confounding."
    )
  )

  ini({
    # ==================================================================
    # Babel 2026 evaluated 'linear and logarithmic logistic regression
    # analyses' (Methods, Exposure-Response Analyses) but does NOT
    # tabulate the coefficients of the selected model for any endpoint.
    # The two values below were recovered by digitising the fitted
    # solid line of Figure 2 at 300 dpi (278 usable pixel columns over
    # CavgADC 0.6 to 9.7 ug/mL) and fitting both candidate forms.
    #
    # Provenance: DIGITISED FROM FIGURE 2, not read from the text or a
    # table. This is a reporting gap - no competing printed value
    # exists - which is the only circumstance under which the library
    # admits a figure-derived parameter.
    #
    # Form selection is unambiguous: the logarithmic form fits the
    # digitised curve with a root-mean-square residual of 0.19
    # percentage points, the linear form with 0.70, and the linear form
    # additionally predicts a 3.8% response probability at CavgADC = 0
    # where the printed curve is visibly pinned to zero.
    #
    # Cross-check against printed values: Babel 2026 Table 1 gives
    # simulated median ORR probabilities of 31.4% at 1.9 mg/kg Q2W and
    # 22.2% at 1.6 mg/kg Q2W. For a log-logistic those two imply an
    # exposure ratio of exp((logit(0.314) - logit(0.222)) / b) = 1.213,
    # against a nominal dose ratio of 1.9/1.6 = 1.1875 - agreement to
    # 2.1% on a quantity that does not involve the intercept at all.
    # Solving for the intercept, the implied median CavgADC is
    # 6.30 ug/mL at 1.9 mg/kg, inside the Figure 2 observed range.
    # ==================================================================
    logit_ref  <- -5.292; label("Logit of the response probability at CavgADC = 1 ug/mL (unitless logit)")                           # digitised from Babel 2026 Figure 2; no printed value exists
    e_cav_orr  <- 2.450;  label("Log-odds of overall response per e-fold increase in conjugate CavgADC (unitless logit)")             # digitised from Babel 2026 Figure 2; no printed value exists

    # ==================================================================
    # The source fits a binomial logistic regression with an exact
    # Bernoulli likelihood: there is no between-subject random effect
    # and no residual error. rxode2 requires an observation
    # declaration, so the deterministic probability is emitted with a
    # tiny placeholder additive residual, mirroring
    # Oniki_2018_nafld_risk.R and the Fukae_2024_valemetostat_* family.
    # ==================================================================
    addSd_prob_orr_central <- fixed(0.001); label("Placeholder additive residual SD on the typical-value response probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Linear predictor on the natural-log exposure scale. CAV is
    # uncentred, so logit_ref is the logit at CavgADC = 1 ug/mL - a
    # value below the observed range, which is why it is large and
    # negative and should not be read as a reference-patient
    # probability.
    logit_orr_central <- logit_ref + e_cav_orr * log(CAV)
    prob_orr_central  <- expit(logit_orr_central)

    # Deterministic probability of overall response by Independent
    # Central Review. Downstream callers can sample binary outcomes
    # with rbinom(n, 1, prob_orr_central) on the rxSolve output.
    prob_orr_central ~ add(addSd_prob_orr_central)
  })
}
