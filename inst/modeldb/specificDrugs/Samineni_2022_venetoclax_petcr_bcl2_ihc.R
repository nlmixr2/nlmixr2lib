Samineni_2022_venetoclax_petcr_bcl2_ihc <- function() {
  description <- paste(
    "Binomial logistic-regression exposure-EFFICACY model for PET-CT",
    "complete response at the end of treatment as a function of the",
    "venetoclax steady-state nominal-dose AUC, fitted separately in BCL-2",
    "IHC-high (IHC 2+/3+, the BCL-2-IHC-positive subpopulation; n = 101,",
    "65 events) and BCL-2 IHC-low (IHC 0/1+; n = 72, 54 events) patients",
    "with previously untreated DLBCL in CAVALLI phase 2 (Samineni 2022",
    "Figure 4c). The two fitted curves are packaged in one file with",
    "TUM_BCL2_IHC_HIGH selecting the stratum. The paper prints no",
    "coefficients; they were digitised from the figure by the",
    "maintainers. Neither slope is significant (p = 0.622 high, 0.704",
    "low). No PK layer: the exposure is supplied as a data column derived",
    "from Samineni_2022_venetoclax. No random effects and no residual",
    "error (Bernoulli likelihood).",
    sep = " "
  )
  reference <- paste(
    "Samineni D, Huang W, Gibiansky L, Ding H, Zhang R, Li C, Sinha A,",
    "Rajwanshi R, Humphrey K, Bazeos A, Salem AH, Miles D.",
    "Population Pharmacokinetics and Exposure-Response Analyses for",
    "Venetoclax in Combination with R-CHOP in Relapsed/Refractory and",
    "Previously Untreated Patients with Diffuse Large B Cell Lymphoma.",
    "Adv Ther. 2022;39(1):598-618. doi:10.1007/s12325-021-01919-z.",
    "Fitted curves of Figure 4c, digitised (the source prints no coefficients for them).",
    sep = " "
  )
  vignette <- "Samineni_2022_venetoclax"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_VENETOCLAX_SS covariate column)",
    concentration = "prob_pet_cr (probability of PET-CT complete response at end of treatment, 0-1; also logit_pet_cr)"
  )

  covariateData <- list(
    AUC_VENETOCLAX_SS = list(
      description = "Individual venetoclax steady-state nominal-dose AUC over one day, AUCss nominal = Dnom * F1 / (CL/F) (Samineni 2022 Eq. 1), from empirical Bayes post hoc PK parameters, the protocol-assigned dose and day-1 covariates.",
      units = "ug*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Supplied as data: this model has no PK layer. The source derived it",
        "from the companion population PK model, packaged as",
        "Samineni_2022_venetoclax. Enters the logit uncentred, per 1",
        "ug*day/mL. The source exposure axis (Figures 4 and 5) runs from",
        "about 0.20 to 1.83 ug*day/mL with exposure quartile bounds near",
        "0.43, 0.55 and 0.68. Computing Eq. 1 from the printed Table 2 values",
        "gives exposures about 1.8-fold higher than that axis; see the",
        "vignette before chaining the PK model into this one.",
        sep = " "
      ),
      source_name = "AUCss nominal"
    ),
    TUM_BCL2_IHC_HIGH = list(
      description = "Tumour BCL-2 protein expression high by immunohistochemistry: 1 = IHC 2+ or 3+ (BCL-2 IHC positive), 0 = IHC 0 or 1+.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (BCL-2 IHC low, scores 0/1+)",
      notes = "Selects which of the two separately fitted Figure 4c curves applies. CAVALLI defined BCL-2 IHC positivity by a cut-off of 50% medium/high staining (Methods). Figure 4c header: high 65 events/101 patients, low 54 events/72 patients; the high group is the BCL-2-IHC-positive subpopulation (n = 101) of the Results.",
      source_name = "BCL-2.IHG1"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 173L,
    n_studies = 1L,
    n_observations = "173 binary PET-CT complete-response records, one per patient",
    age_range = "18-85 years (median 64) in the CAVALLI PK population (Table 1)",
    weight_range = "46.0-182 kg (median 77.0) in the CAVALLI PK population (Table 1)",
    sex_female_pct = 44.4,
    race_ethnicity = c(
      White = 71.1,
      Asian = 2.2,
      Black = 1.7,
      `Native Hawaiian or Pacific Islander` = 1.3,
      Unknown = 23.7
    ),
    disease_state = "Previously untreated (1L) diffuse large B-cell lymphoma, phase 2 of CAVALLI (GO27878); patients with a BCL-2 IHC result",
    dose_range = "Venetoclax 800 mg once daily on days 4-10 of cycle 1 and days 1-10 of cycles 2-8 of 21-day cycles, with R-CHOP",
    regions = "North America, Western and Eastern Europe and other regions (CAVALLI)",
    notes = "Exposure-efficacy analysis population (199 phase 2 patients). Demographic percentages are those of the 232-patient CAVALLI PK population in Table 1. PET-CR was assessed at the end-of-treatment visit, 6-8 weeks after day 1 of the last cycle; patients without an end-of-treatment assessment count as non-responders (Figure 4 legend)."
  )

  ini({
    # Samineni 2022 Figure 4c fitted line(s), digitised by the maintainers from the
    # 300-dpi figure image and fitted logit-linear in AUC (residual RMSE
    # under 0.006 on the logit scale). The paper prints no coefficients.
    logit_ref_bcl2_ihc_low <- 1.402; label("Logit of the PET-CR probability at zero venetoclax exposure, BCL-2 IHC low (unitless logit)") # Figure 4c blue (LOW) curve, digitised; p(0) = 0.803
    e_auc_logit_bcl2_ihc_low <- -0.515; label("Log-odds of PET-CR per 1 ug*day/mL of venetoclax AUCss nominal, BCL-2 IHC low (unitless logit)") # Figure 4c blue (LOW) curve, digitised; figure p = 0.704
    logit_ref_bcl2_ihc_high <- 0.837; label("Logit of the PET-CR probability at zero venetoclax exposure, BCL-2 IHC high (unitless logit)") # Figure 4c red (HIGH) curve, digitised; p(0) = 0.698
    e_auc_logit_bcl2_ihc_high <- -0.422; label("Log-odds of PET-CR per 1 ug*day/mL of venetoclax AUCss nominal, BCL-2 IHC high (unitless logit)") # Figure 4c red (HIGH) curve, digitised; figure p = 0.622

    # Logistic regression (glm, logit link): Bernoulli likelihood, no random
    # effects and no residual error. The tiny additive residual below exists
    # only so rxode2 has an error model to attach to the probability; it is
    # not a published quantity.
    addSd_prob_pet_cr <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Two separately fitted regressions; the indicator selects the stratum.
    logit_pet_cr <- (1 - TUM_BCL2_IHC_HIGH) * (logit_ref_bcl2_ihc_low + e_auc_logit_bcl2_ihc_low * AUC_VENETOCLAX_SS) +
      TUM_BCL2_IHC_HIGH * (logit_ref_bcl2_ihc_high + e_auc_logit_bcl2_ihc_high * AUC_VENETOCLAX_SS)
    prob_pet_cr <- expit(logit_pet_cr)

    prob_pet_cr ~ add(addSd_prob_pet_cr)
  })
}
