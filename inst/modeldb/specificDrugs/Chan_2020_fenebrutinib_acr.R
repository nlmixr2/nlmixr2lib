Chan_2020_fenebrutinib_acr <- function() {
  description <- paste(
    "Longitudinal logistic exposure-response model for the probability of",
    "achieving ACR20, ACR50 and ACR70 responses in rheumatoid arthritis",
    "patients treated with fenebrutinib (GDC-0853) or placebo in the phase",
    "2 ANDES trial (Chan 2020). The three ACR thresholds are fitted",
    "simultaneously: each has its own baseline logit, ACR20 has its own",
    "placebo onset time and ACR50 / ACR70 share one, and all three share a",
    "sigmoidal Emax placebo time course, one Markov term on the patient's",
    "previous response for the same threshold, a single subject-level",
    "random effect on the logit, and a single Emax function of the",
    "individual steady-state daily fenebrutinib AUC whose maximum is scaled",
    "by enrollment region (Eastern Europe reference, US, Latin America).",
    "The AUC is a covariate column produced by Chan_2020_fenebrutinib."
  )
  reference <- paste(
    "Chan P, Yu J, Chinn L, Prohn M, Huisman J, Matzuka B, Hanley W,",
    "Tuckwell K, Quartino A. Population Pharmacokinetics, Efficacy",
    "Exposure-response Analysis, and Model-based Meta-analysis of",
    "Fenebrutinib in Subjects with Rheumatoid Arthritis. Pharm Res.",
    "2020;37(2):25. doi:10.1007/s11095-019-2752-y"
  )
  vignette <- "Chan_2020_fenebrutinib"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "prob_acr20 / prob_acr50 / prob_acr70 (probability of response, 0-1)"
  )

  covariateData <- list(
    AUC_FENEBRUTINIB = list(
      description = "Individual steady-state daily (0-24 h) AUC of fenebrutinib",
      units = "h*ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Empirical Bayes (post hoc) estimate from the Chan 2020 popPK model",
        "(Chan_2020_fenebrutinib) at the nominal dose (Methods: 'total daily",
        "area-under-the-concentration-time-curve (AUC) at steady state as",
        "predicted using nominal dose'). 0 for placebo. Time-fixed per",
        "subject; the drug effect is applied at every visit, not ramped in",
        "over time. Figure S3 medians: about 1200 (50 mg QD), 3400 (150 mg",
        "QD) and 7000 (200 mg BID) h*ng/mL."
      ),
      source_name = "AUC"
    ),
    PDV = list(
      description = "Previous observed responder status (0/1) for the same ACR threshold as the current record",
      units = "(count)",
      type = "count",
      reference_category = NULL,
      notes = paste(
        "Markov element of Model S2 ('MARKOV = THETA(4)*PDV'). A binary",
        "0/1 count here: 1 if the patient met the same ACR threshold at the",
        "previous visit. At the first post-baseline visit (day 7) PDV = 0,",
        "because every patient is a non-responder at baseline by",
        "definition. Visits were on days 7, 14, 28, 56 and 84, with",
        "non-responder imputation after discontinuation."
      ),
      source_name = "PDV"
    ),
    REGION_USA = list(
      description = "Enrolled at a United States study site: 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (enrolled outside the US)",
      notes = paste(
        "Model S2 REGION = 1 (GEOGUSEUKORLATAM = 1). Multiplies the drug",
        "Emax by 2.13. 22 of 467 patients (5%) (Table S3)."
      ),
      source_name = "GEOGUSEUKORLATAM"
    ),
    REGION_EASTEUROPE = list(
      description = "Enrolled at an Eastern European study site: 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not Eastern Europe)",
      notes = paste(
        "Model S2 REGION = 2 (GEOGUSEUKORLATAM = 2), the most common and",
        "reference region: 290 of 467 patients (62%) (Table S3). The",
        "third region, Latin America (Model S2 REGION = 3, source codes 3",
        "and 4; 155 patients, 33%), has no canonical indicator and is",
        "derived as 1 - REGION_USA - REGION_EASTEUROPE, so a Latin American",
        "patient is supplied with both indicators 0."
      ),
      source_name = "GEOGUSEUKORLATAM"
    )
  )

  covariatesDataExcluded <- list(
    RHEUMATOID_FACTOR = list(
      description = "Baseline rheumatoid factor",
      units = "IU/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Model S2 carries 'RFCOV = (RF/124)**THETA(15)' on the drug effect",
        "(124 IU/mL = Table S3 median), and the Results name baseline RF as",
        "a significant covariate, but Table S4 reports no estimate for",
        "THETA(15) and the Discussion states that 'region was the only",
        "covariate in the longitudinal E-R model'. THETA(15) is therefore",
        "treated as fixed at 0 (RFCOV = 1) and the term is omitted. The",
        "same pattern (a THETA or ETA coded in the control stream but absent",
        "from the estimate table) recurs for the popPK additive error and",
        "for five popPK etas."
      ),
      source_name = "RFIN"
    ),
    PRIOR_TNF = list(
      description = "Prior anti-TNF therapy with inadequate response (TNF-IR, cohort 2)",
      units = "(binary)",
      type = "binary",
      reference_category = NULL,
      notes = "Tested but not significant (Discussion: 'treatment history (MTX-IR vs. TNF-IR) was tested but was not a significant covariate in the E-R model'). 98 of 467 patients (Table S3).",
      source_name = "TNF-IR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 467L,
    n_studies = 1L,
    age_range = "19-75 years (median 52)",
    weight_range = "38-153 kg (median 71)",
    sex_female_pct = 80.5,
    disease_state = "Adults with moderate to severe active seropositive rheumatoid arthritis on stable methotrexate, with inadequate response to methotrexate (cohort 1) or to anti-TNF therapy (cohort 2)",
    dose_range = "Placebo or fenebrutinib 50 mg QD, 150 mg QD or 200 mg BID tablets for 12 weeks (adalimumab arm not included)",
    regions = "Eastern Europe 62%, Latin America 33%, US 5%",
    n_observations = "7005 binary observations (2335 each for ACR20, ACR50, ACR70) on days 7, 14, 28, 56 and 84",
    notes = "Chan 2020 Methods 'ACR20, ACR50, and ACR70', Table S3 (demographics) and Table S4 (estimates); Laplacian estimation in NONMEM."
  )

  ini({
    # Baseline logits, one per ACR threshold (Model S2 TVBASE). Table S4
    # prints them untransformed (they carry an RSE).
    logit_bl_acr20 <- -4.62; label("Baseline logit of ACR20 response (unitless)")  # Table S4, theta1 'ACR20 baseline' = -4.62 (RSE 10%)
    logit_bl_acr50 <- -6.47; label("Baseline logit of ACR50 response (unitless)")  # Table S4, theta2 'ACR50 baseline' = -6.47 (RSE 6.4%)
    logit_bl_acr70 <- -8.3; label("Baseline logit of ACR70 response (unitless)")  # Table S4, theta3 'ACR70 baseline' = -8.3 (RSE 5.7%)

    e_pdv_acr <- 0.934; label("Markov coefficient on the previous response for the same threshold (logit units)")  # Table S4, theta4 'Markov component' = 0.934 (RSE 17%)

    # Placebo (time) effect, a sigmoidal Emax in study day. THETA(8),
    # THETA(11) and THETA(12) enter Model S2 through EXP(), so Table S4
    # prints them back-transformed (no RSE, CI only).
    emax_pbo_acr <- 3.41; label("Maximum placebo effect over time (logit units)")  # Table S4, theta7 'Maximum placebo effect over time' = 3.41 (RSE 14.5%)
    lt50_pbo_acr20 <- log(21.5); label("Log time of 50% placebo effect for ACR20 (day)")  # Table S4, theta8 'Time of 50% placebo effect - ACR20 (d)' = 21.5 (95% CI 16.6-27.7)
    lt50_pbo_acr5070 <- log(32.8); label("Log time of 50% placebo effect for ACR50 and ACR70 (day)")  # Table S4, theta11 'Time of 50% placebo effect - ACR50 and ACR70 (d)' = 32.8 (95% CI 26.1-41.1)
    lhill_pbo_acr <- log(2.52); label("Log Hill coefficient of the placebo time course (unitless)")  # Table S4, theta12 'Hill coefficient on time course' = 2.52 (95% CI 1.5-4.23)

    # Drug effect: Emax in AUC, one term for all three thresholds.
    emax_acr <- 1.39; label("Maximum fenebrutinib effect, Eastern Europe (logit units)")  # Table S4, theta9 'Max drug effect over time - Eastern Europe' = 1.39 (RSE 29.8%)
    lec50_acr <- log(2650); label("Log daily steady-state AUC giving 50% of the maximum drug effect (h*ng/mL)")  # Table S4, theta10 'Exposure at which 50% drug effect (AUC ng.hr/mL)' = 2650 (95% CI 675-10400)
    # Model S2 comments THETA(13) / THETA(14) as a 'Fraction of EMAX'
    # multiplying THETA(9) (DEFF = THETA(9)*REGCOV*...); Table S4 prints them
    # with an RSE, i.e. untransformed, so they are the multipliers themselves.
    e_region_usa_emax_acr <- 2.13; label("Multiplier on the drug Emax for US patients (relative to Eastern Europe)")  # Table S4, theta13 'Max drug effect over time - US' = 2.13 (RSE 39.7%)
    e_region_latam_emax_acr <- 2.06; label("Multiplier on the drug Emax for Latin American patients (relative to Eastern Europe)")  # Table S4, theta14 'Max drug effect over time - Latin America' = 2.06 (RSE 25.8%)

    # Subject-level additive random effect on the logit (Model S2 ETA(1)).
    # Table S4 labels it 'IOV on baseline', but Model S2 has a single
    # subject-level ETA in a $PRED model and the Results describe 'an
    # interindividual variability term ... added on baseline'.
    etalogit_bl_acr ~ 4.85 # Table S4, omega16.1 'omega2 IOV on baseline' = 4.85 (RSE 13.3%; shrinkage 20.45%)

    # Placeholder residuals. The source likelihood is Bernoulli (Model S2
    # Y = P or 1 - P) and estimates no residual error.
    addSd_prob_acr20 <- fixed(0.001); label("Placeholder additive residual SD on prob_acr20 (not from source)")  # not from source; see vignette Assumptions and deviations
    addSd_prob_acr50 <- fixed(0.001); label("Placeholder additive residual SD on prob_acr50 (not from source)")  # not from source; see vignette Assumptions and deviations
    addSd_prob_acr70 <- fixed(0.001); label("Placeholder additive residual SD on prob_acr70 (not from source)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Region multiplier on the drug Emax (Model S2 REGCOV; Eastern Europe = 1).
    region_latam <- 1 - REGION_USA - REGION_EASTEUROPE
    regcov <- e_region_usa_emax_acr^REGION_USA * e_region_latam_emax_acr^region_latam

    # Placebo time course (Model S2 TMAX), with `time` = study day.
    hill <- exp(lhill_pbo_acr)
    t50_20 <- exp(lt50_pbo_acr20)
    t50_5070 <- exp(lt50_pbo_acr5070)
    pbo20 <- emax_pbo_acr * time^hill / (time^hill + t50_20^hill)
    pbo5070 <- emax_pbo_acr * time^hill / (time^hill + t50_5070^hill)

    # Drug effect (Model S2 DEFF; RFCOV = 1, see covariatesDataExcluded).
    deff <- emax_acr * regcov * AUC_FENEBRUTINIB / (AUC_FENEBRUTINIB + exp(lec50_acr))

    # Logit per threshold. PDV is the previous response for the threshold of
    # the current record, so on a data row only the output matching that
    # row's endpoint uses the right PDV.
    common <- deff + e_pdv_acr * PDV + etalogit_bl_acr
    logit_acr20 <- logit_bl_acr20 + pbo20 + common
    logit_acr50 <- logit_bl_acr50 + pbo5070 + common
    logit_acr70 <- logit_bl_acr70 + pbo5070 + common

    prob_acr20 <- expit(logit_acr20)
    prob_acr50 <- expit(logit_acr50)
    prob_acr70 <- expit(logit_acr70)

    prob_acr20 ~ add(addSd_prob_acr20)
    prob_acr50 ~ add(addSd_prob_acr50)
    prob_acr70 ~ add(addSd_prob_acr70)
  })
}
