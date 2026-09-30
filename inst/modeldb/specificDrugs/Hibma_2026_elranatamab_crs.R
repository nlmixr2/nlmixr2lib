Hibma_2026_elranatamab_crs <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-safety model for any-grade ",
    "cytokine release syndrome (CRS) after the FIRST step-up priming dose ",
    "(12 mg SC on Cycle 1 Day 1) of the BCMA x CD3 bispecific antibody ",
    "elranatamab in adults with relapsed or refractory multiple myeloma ",
    "(Hibma 2026; MagnetisMM-3, N = 183, 79 events). The probability of ",
    "CRS is expit(-7.65 + 1.14 * TUM_BURDEN_HIGH + 1.35 * log(CTROUGH)), ",
    "where CTROUGH is the individual FREE elranatamab trough concentration ",
    "in ng/mL on Day 4, just before the second (32 mg) step-up dose, and ",
    "TUM_BURDEN_HIGH flags high baseline tumour burden (ESM Table S2). ",
    "The log is natural (odds ratio 3.86 = exp(1.35) per unit of ",
    "log(Ctrough)). There is no PK layer and no ODE: the exposure is ",
    "supplied as a data column, derived in the source from post-hoc ",
    "estimates of the companion population PK model Hibma_2026_elranatamab. ",
    "No random effect and no residual error are estimated (Bernoulli ",
    "likelihood). No exposure-response relationship was significant after ",
    "the second step-up dose, and none could be fitted after the first ",
    "full dose, so this is the paper's only CRS model."
  )
  reference <- paste(
    "Hibma JE, Irby D, Liu A, Elmeliegy M, King LE, Gifondorwa D, Jiang S,",
    "Poels KE, Soltantabar P, Lon H-K, Shtylla B, Wang D, Williams JH,",
    "Nicholas T. Elranatamab population pharmacokinetics and",
    "exposure-response for cytokine release syndrome in patients with",
    "relapsed or refractory multiple myeloma.",
    "Clin Pharmacokinet. 2026;65:1173-1192.",
    "doi:10.1007/s40262-026-01663-z.",
    "Exposure derives from the companion population PK model; see",
    "modellib('Hibma_2026_elranatamab').",
    sep = " "
  )
  vignette <- "Hibma_2026_elranatamab"
  units <- list(
    time = "n/a (static landmark exposure-safety regression for the first step-up dose interval; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_crs_stepup1 (probability of any-grade CRS after the first step-up priming dose, 0-1; also logit_crs_stepup1)"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Individual FREE elranatamab plasma trough concentration on Cycle 1 Day 4, after the 12 mg first step-up priming dose and before the 32 mg second step-up dose.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "(a) FREE (not total) elranatamab. (b) Landmark: the Day 4 trough",
        "of the first step-up dose interval (Days 1-4), which the paper",
        "notes 'also corresponds to the overall Cmax for the same interval'",
        "because SC absorption is still rising at Day 4. (c) Derived from",
        "individual post-hoc parameters of the final population PK model,",
        "packaged as Hibma_2026_elranatamab (free elranatamab `Cc` at",
        "time 3 days after a 12 mg SC dose). Enters on the NATURAL LOG",
        "scale with the concentration in ng/mL: ESM Table S2 prints the",
        "term as 'Log(Ctrough (ng/mL))' with estimate 1.35 and odds ratio",
        "3.86 = exp(1.35), so the unit is load-bearing -- supplying nM",
        "instead of ng/mL would shift the logit by 1.35 * log(148) = 6.7.",
        "Figure 6a spans log(Ctrough) of about 3.6-7 (roughly 35-1100",
        "ng/mL); the geometric mean peak free concentration after the first",
        "step-up dose was 212 ng/mL (90% CV) (Section 3.8)."
      ),
      source_name = "Log(Ctrough (ng/mL))"
    ),
    TUM_BURDEN_HIGH = list(
      description = "High baseline multiple-myeloma tumour burden indicator: 1 = high, 0 = low or intermediate (reference).",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Hibma 2026 Table 1 footnote: high tumour burden = bone-marrow",
        "plasma-cell infiltrate >= 80% OR serum M-spike >= 5 g/dL OR serum",
        "free light chain >= 5000 mg/L; low = infiltrate < 50% AND M-spike",
        "< 3 g/dL AND free light chain < 3000 mg/L; patients meeting",
        "neither definition are intermediate. The regression pools low and",
        "intermediate as the reference level ('reference tumor burden",
        "(i.e., low/intermediate)', Figure 6 caption). MagnetisMM-3",
        "exposure-response population (N = 183; printed under the Cohort A",
        "column of Table 1, but the counts sum to 183): 131 (72%) low or",
        "intermediate, 38 (21%) high, 14 (7.7%) missing. The paper does not say how the missing category entered",
        "the regression; ESM Table S2 prints a single 'TumorBurden High'",
        "coefficient, so missing is not a separate level here."
      ),
      source_name = "TumorBurden High"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      notes = "Prespecified in the full multivariable CRS model (Section 2.7) and removed by backward elimination (alpha = 0.01)."
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      notes = "Prespecified in the full multivariable CRS model (Section 2.7) and removed by backward elimination."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      notes = "Prespecified in the full multivariable CRS model (Section 2.7) and removed by backward elimination."
    ),
    SBCMA = list(
      description = "Baseline soluble BCMA concentration",
      units = "nM",
      type = "continuous",
      notes = "Prespecified in the full multivariable CRS model (Section 2.7) and removed by backward elimination. Creatinine clearance, platelet count, lactate dehydrogenase, race, myeloma type, immunogenicity, premedication, steroid route, prior CAR-T and prior BCMA-directed treatment were likewise screened and removed."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 183L,
    n_studies = 1L,
    n_observations = "183 binary any-grade CRS outcomes after the first step-up priming dose (79 events, 43%; 104 non-events)",
    age_range = "36-89 years (MagnetisMM-3 Cohorts A and B; medians 68 and 67; Table 1)",
    weight_range = "36.5-159.6 kg (MagnetisMM-3 Cohorts A and B; medians 72.0 and 69.3; Table 1)",
    disease_state = "Relapsed or refractory multiple myeloma",
    dose_range = "SC elranatamab 12 mg on C1D1, 32 mg on C1D4, 76 mg on C1D8 then 76 mg QW (two-step-up priming regimen)",
    regions = "Multinational (MagnetisMM-3, NCT04649359)",
    tumor_burden = "MagnetisMM-3 exposure-response population: low or intermediate 131 (72%), high 38 (21%), missing 14 (7.7%) (Table 1)",
    notes = paste(
      "Data cutoff March 2024. After the second step-up dose 35 CRS events",
      "(19%) occurred with no significant exposure-response relationship;",
      "after the first full dose 13 events (7%), too few for a regression",
      "(Section 3.6). Model discrimination AUROC 0.746 (ESM Figure S5)."
    )
  )

  ini({
    # ESM Table S2 ('Estimates from the ER model for CRS any grade after the
    # first step-up priming dose'), n = N = 183, change in deviance 37.6 on
    # 2 df, logLik -106.33, AIC 218.66. Fitted with R glm(family =
    # binomial). Covariates are neither centred nor scaled, so the intercept
    # is the logit at CTROUGH = 1 ng/mL and reference tumour burden.
    # Each printed odds ratio is exp() of its estimate: exp(1.35) = 3.857
    # (printed 3.86); exp(1.14) = 3.127 (printed 3.11 -- consistent under
    # rounding, e.g. b = 1.135 gives round(b, 2) = 1.14 and
    # round(exp(b), 2) = 3.11); these checks are re-run in the vignette.
    logit_ref <- -7.65; label("Logit of any-grade CRS after the first step-up dose at CTROUGH = 1 ng/mL and low/intermediate tumour burden (unitless logit)") # ESM Table S2 Intercept -7.65 (95% CI -10.5, -5.03), Z -5.46
    e_tum_burden_high_logit <- 1.14; label("Log-odds of CRS for high vs low/intermediate baseline tumour burden (unitless logit)") # ESM Table S2 'TumorBurden High' 1.14 (0.335, 1.97), Z 2.74, P 0.0062, OR 3.11 (1.4, 7.17)
    e_ctrough_logit <- 1.35; label("Log-odds of CRS per unit increase in the natural log of the Day 4 free elranatamab trough in ng/mL (unitless logit)") # ESM Table S2 'Log(Ctrough (ng/mL))' 1.35 (0.869, 1.88), Z 5.26, OR 3.86 (2.38, 6.55)

    # No between-subject variability and no residual error: Bernoulli
    # likelihood. The tiny fixed additive residual exists only so rxode2 has
    # an error model to attach to the probability; it is NOT a published
    # quantity (see vignette Assumptions and deviations).
    addSd_prob_crs_stepup1 <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Linear predictor (ESM Table S2). log() is the natural log, as in R's
    # glm(), confirmed by OR = exp(estimate).
    logit_crs_stepup1 <- logit_ref +
      e_tum_burden_high_logit * TUM_BURDEN_HIGH +
      e_ctrough_logit * log(CTROUGH)

    prob_crs_stepup1 <- expit(logit_crs_stepup1)

    prob_crs_stepup1 ~ add(addSd_prob_crs_stepup1)
  })
}
