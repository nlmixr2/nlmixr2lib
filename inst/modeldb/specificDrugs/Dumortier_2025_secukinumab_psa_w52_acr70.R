Dumortier_2025_secukinumab_psa_w52_acr70 <- function() {
  description <- paste(
    "Landmark logistic sigmoid-Emax exposure-response model for the ",
    "probability of at least 70% improvement in the American College of ",
    "Rheumatology response criteria (ACR70) at Week 52 in adults with ",
    "psoriatic arthritis treated with secukinumab (Dumortier 2025; FUTURE ",
    "1-5: NCT01392326, NCT01752634, NCT01989468, NCT02294227, NCT02404350). ",
    "logit(p) = E0 + (Emax - E0) * Cmin^2 / (EC50^2 + Cmin^2) with the Hill ",
    "coefficient fixed at 2, where Cmin is the popPK-predicted pre-dose ",
    "secukinumab serum concentration at the Week 52 efficacy assessment ",
    "(Cmin,ss), supplied in the CTROUGH covariate column in ug/mL. E0, Emax ",
    "and EC50 are not printed in the paper; they were recovered by ",
    "pixel-digitising the fitted curve in Figure S3 and refitting the same ",
    "equation (RMSE 0.0014 probability units over 235 traced points). No PK ",
    "layer and no ODE: drive CTROUGH from the companion popPK model ",
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
    time = "n/a (static landmark model for the Week 52 assessment; no time dimension)",
    dosing = "n/a (no dose events; secukinumab exposure enters through the CTROUGH covariate)",
    concentration = "prob_acr70 (probability of response, 0-1); CTROUGH in ug/mL"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Model-predicted pre-dose secukinumab serum concentration at the Week 52 efficacy assessment",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "PopPK-predicted pre-dose secukinumab serum concentration at the Week ",
        "52 efficacy assessment (Cmin,ss). Dumortier 2025 used the individual ",
        "prediction of the population PK model (Dumortier_2025_secukinumab) at ",
        "the assessment time; a user can obtain it by simulating that model to ",
        "Week 52. The analysis set is patients randomised to SC secukinumab ",
        "with Week 52 efficacy data (no imputation); the 10 mg/kg IV-loaded ",
        "arms of FUTURE 1 are included at Week 52."
      ),
      source_name = "Cmin,ss"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1595L,
    n_studies = 5L,
    disease_state = "Psoriatic arthritis (PsA)",
    dose_range = "75, 150 or 300 mg SC q4w, with weekly SC loading or 10 mg/kg IV loading at Weeks 0, 2 and 4",
    regions = "Multinational",
    notes = paste(
      "Exposure-efficacy analysis set of Dumortier 2025 Supplementary ",
      "Methods: patients randomised to SC secukinumab with Week 52 efficacy ",
      "data (no imputation); the 10 mg/kg IV-loaded arms of FUTURE 1 are ",
      "included at Week 52. 1595 PsA patients contributed to the Week 52 ",
      "analysis (main text, Results). The per-endpoint number of patients is ",
      "not reported. The fitted curve is already at 0 at the left edge of the ",
      "plotted range (Cmin about 1 ug/mL), so E0 is only bounded from above: ",
      "E0 between -10 and -20 reproduces the curve to within 0.0017. The ",
      "least-squares optimum is reported."
    )
  )

  ini({
    e0 <- -14.683; label("Logit of the response probability at zero secukinumab concentration, E0 (unitless logit)") # digitised from Figure S3 panel 'ACR70' (fitted curve; not printed in the paper)
    emax <- -0.851; label("Logit of the response probability at infinite secukinumab concentration, Emax (unitless logit)") # digitised from Figure S3 panel 'ACR70' (fitted curve; not printed in the paper)
    lec50 <- log(1.49); label("Pre-dose concentration giving half of the logit shift from E0 to Emax, EC50 (ug/mL)") # digitised from Figure S3 panel 'ACR70' (fitted curve; not printed in the paper)
    lhill <- fixed(log(2)); label("Hill coefficient (unitless)") # Supplementary Methods 'this parameter was fixed at a value of 2'
    addSd_prob_acr70 <- fixed(0.001); label("Placeholder additive residual SD on the response probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Supplementary Methods: log(p / (1 - p)) = E0 + (Emax - E0) * Cmin^h / (EC50^h + Cmin^h)
    logit_acr70 <- e0 + (emax - e0) * CTROUGH^hill / (ec50^hill + CTROUGH^hill)
    prob_acr70 <- expit(logit_acr70)

    prob_acr70 ~ add(addSd_prob_acr70)
  })
}
