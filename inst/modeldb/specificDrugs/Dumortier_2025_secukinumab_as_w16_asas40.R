Dumortier_2025_secukinumab_as_w16_asas40 <- function() {
  description <- paste(
    "Landmark logistic sigmoid-Emax exposure-response model for the ",
    "probability of at least 40% improvement in the Assessment of ",
    "SpondyloArthritis International Society criteria (ASAS40) at Week 16 ",
    "in adults with ankylosing spondylitis treated with secukinumab ",
    "(Dumortier 2025; MEASURE 2, 4 and 5: NCT01649375, NCT02159053, ",
    "NCT02896127). logit(p) = E0 + (Emax - E0) * Cmin^2 / (EC50^2 + Cmin^2) ",
    "with the Hill coefficient fixed at 2, where Cmin is the ",
    "popPK-predicted pre-dose secukinumab serum concentration at the Week ",
    "16 efficacy assessment (Cmin,16W); 0 for placebo patients, supplied in ",
    "the CTROUGH covariate column in ug/mL. E0, Emax and EC50 are not ",
    "printed in the paper; they were recovered by pixel-digitising the ",
    "fitted curve in Figure S7 and refitting the same equation (RMSE 0.0005 ",
    "probability units over 198 traced points). No PK layer and no ODE: ",
    "drive CTROUGH from the companion popPK model ",
    "Dumortier_2025_secukinumab."
  )
  reference <- paste(
    "Dumortier T, Valenzuela G, Churchill M, Mijatovic J, Bruin G,",
    "Pricop L, Richards H, Renard D, Singhal A, Marathe A.",
    "Model-Informed Drug Development-Based Bridging from Subcutaneous to",
    "Intravenous Secukinumab Dosing: Approval in Psoriatic Arthritis and",
    "Axial Spondyloarthritis. Clin Pharmacol Ther. 2025;118(2):480-488.",
    "doi:10.1002/cpt.3716."
  )
  vignette <- "Dumortier_2025_secukinumab"
  units <- list(
    time = "n/a (static landmark model for the Week 16 assessment; no time dimension)",
    dosing = "n/a (no dose events; secukinumab exposure enters through the CTROUGH covariate)",
    concentration = "prob_asas40 (probability of response, 0-1); CTROUGH in ug/mL"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Model-predicted pre-dose secukinumab serum concentration at the Week 16 efficacy assessment",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "PopPK-predicted pre-dose secukinumab serum concentration at the Week ",
        "16 efficacy assessment (Cmin,16W); 0 for placebo patients. Dumortier ",
        "2025 used the individual prediction of the population PK model ",
        "(Dumortier_2025_secukinumab) at the assessment time; a user can obtain ",
        "it by simulating that model to Week 16. The analysis set is patients ",
        "randomised to placebo or SC secukinumab with Week 16 efficacy data (no ",
        "imputation); arms with a 10 mg/kg IV loading regimen are excluded."
      ),
      source_name = "Cmin,16W"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 3L,
    disease_state = "Ankylosing spondylitis (AS)",
    dose_range = "0 (placebo), 75 or 150 mg SC q4w with weekly SC loading",
    regions = "Multinational",
    notes = paste(
      "Exposure-efficacy analysis set of Dumortier 2025 Supplementary ",
      "Methods: patients randomised to placebo or SC secukinumab with Week 16 ",
      "efficacy data (no imputation); arms with a 10 mg/kg IV loading regimen ",
      "are excluded. The Week 16 analysis-set size is not reported. The ",
      "per-endpoint number of patients is not reported."
    )
  )

  ini({
    e0 <- -1.32; label("Logit of the response probability at zero secukinumab concentration, E0 (unitless logit)") # digitised from Figure S7 panel 'ASAS40' (fitted curve; not printed in the paper)
    emax <- 0.058; label("Logit of the response probability at infinite secukinumab concentration, Emax (unitless logit)") # digitised from Figure S7 panel 'ASAS40' (fitted curve; not printed in the paper)
    lec50 <- log(12.6); label("Pre-dose concentration giving half of the logit shift from E0 to Emax, EC50 (ug/mL)") # digitised from Figure S7 panel 'ASAS40' (fitted curve; not printed in the paper)
    lhill <- fixed(log(2)); label("Hill coefficient (unitless)") # Supplementary Methods 'this parameter was fixed at a value of 2'
    addSd_prob_asas40 <- fixed(0.001); label("Placeholder additive residual SD on the response probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Supplementary Methods: log(p / (1 - p)) = E0 + (Emax - E0) * Cmin^h / (EC50^h + Cmin^h)
    logit_asas40 <- e0 + (emax - e0) * CTROUGH^hill / (ec50^hill + CTROUGH^hill)
    prob_asas40 <- expit(logit_asas40)

    prob_asas40 ~ add(addSd_prob_asas40)
  })
}
