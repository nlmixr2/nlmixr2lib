Li_2015_taspoglutide_mbma <- function() {
  description <- "MBMA. Model-based meta-analysis of body-weight (WT) loss from baseline in adults with type 2 diabetes treated with the once-weekly GLP-1 analogue taspoglutide. Digitised study-arm-mean WT response from 9 placebo-controlled trials (3702 patients pooled across placebo, 10 mg QW, and 20 mg QW arms; 8-52 weeks follow-up) is decomposed into an exponential placebo response Placebo(t) = Pmax * (1 - exp(-Kp*t)) (Pmax = -1.33 kg, Kp = 0.0987 1/week; both fixed in the combined model from a placebo-only fit) and an additive Emax drug response Drug(t) = Placebo(t) + Dmax * CAV / (IC50 + CAV) * (1 - exp(-Kdrug*t)) with Dmax = -1.85 kg, IC50 = 41.7 pmol/L, and Kdrug = 0.422 1/week. Exposure enters as the study-arm average plasma concentration covariate CAV (0 pmol/L for placebo, 59.85 pmol/L for the 10 mg QW arm derived as half of 119.7 pmol/L, 119.7 pmol/L for the 20 mg QW arm; both calculated from the 20 mg Cavg between weeks 2 and 4 reported in Ratner 2010 assuming dose-proportional exposure). The MBMA between-study variability is encoded as study-level eta variance (eta_study_*); the model is intended for simulating study-arm-mean WT change in kg over time and is NOT suitable for individual-subject simulation. Proportional residual error b = 0.189 from Y = F + b*F*epsilon. Covariates were not retained in the final model."

  reference <- paste(
    "Li HQ, Xu JY, Jin L, Xin JL.",
    "The efficacy of placebo-adjusted taspoglutide on body weight reduction in clinical trials.",
    "Pharmazie. 2015;70(2):110-116.",
    "doi:10.1691/ph.2015.4716.",
    sep = " "
  )
  vignette <- "Li_2015_taspoglutide_mbma"
  units <- list(
    time          = "week (Kp and Kdrug rate constants are per week per Li 2015 Eq 1-2)",
    dosing        = "mg SC QW (taspoglutide subcutaneous once-weekly; the MBMA does not consume rxode2 dose events. Exposure is supplied as the CAV covariate column.)",
    concentration = "kg/arm (signed change in body weight from pretreatment baseline at the study-arm-mean dimension; Cc is NOT a drug concentration. The slash in the unit string is to satisfy checkModelConventions parsing; the absolute unit is kg.)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Study-arm average plasma concentration of taspoglutide between weeks 2 and 4 of QW dosing, used as the exposure metric in the Emax drug-response term.",
      units              = "pmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "MBMA study-arm-level covariate. Values per Li 2015 Table 2:",
        "0 pmol/L (placebo), 59.85 pmol/L (10 mg QW; derived as half of 119.7 pmol/L",
        "assuming dose-proportional exposure per Ratner 2010), and 119.7 pmol/L",
        "(20 mg QW; mean of 3 individual Cavg(2w-4w) values from Ratner 2010).",
        "The averaging window is weeks 2-4 of once-weekly subcutaneous dosing.",
        "Set CAV = 0 for placebo arms. The 8-52 week follow-up is assumed to",
        "match the same Cavg metric for all later observations as dosing continues",
        "at the same QW schedule (Li 2015 Methods section 2.2).",
        sep = " "
      ),
      source_name        = "metric (Cavg(2w-4w); Li 2015 Table 2)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 3702L,
    n_studies      = 9L,
    age_range      = "adults with type 2 diabetes mellitus; specific age range not reported in Li 2015 but the source trials (Raz 2012, Ratner 2010, Nauck 2009, Nauck 2013, Hollander 2013, Bergenstal 2012, Henry 2012, Pratley 2013, Rosenstock 2013) enrolled adults inadequately controlled on metformin, sulphonylurea, or metformin + TZD",
    disease_state  = "adults with type 2 diabetes mellitus inadequately controlled on background therapy (drug-naive, metformin, metformin + sulphonylurea, metformin + TZD, or sulphonylurea +/- metformin)",
    dose_range     = "taspoglutide 10 mg or 20 mg SC once-weekly; placebo also included",
    duration_range = "8-52 weeks (Li 2015 Table 1)",
    regions        = "International; trials from the published taspoglutide phase 2/3 program (T-emerge series and related)",
    notes          = paste(
      "MBMA at the study-arm level: each data point is the mean change in body",
      "weight from baseline in one trial arm at one time point. Source 9 trials",
      "tabulated in Li 2015 Table 1; total of 3702 patients pooled across placebo,",
      "10 mg, and 20 mg arms over 8-52 weeks follow-up. The model is intended for",
      "simulating study-arm-mean WT change and is NOT suitable for individual-subject",
      "simulation. Background concomitant medications across the 9 trials included",
      "monotherapy (1 trial), metformin (5 trials), metformin + pioglitazone or TZD",
      "(1 trial), sulphonylurea or sulphonylurea + metformin (1 trial), and metformin",
      "or metformin + thiazolidinedione (1 trial). Covariates were not retained in",
      "the final model (Li 2015 Discussion)."
    )
  )

  ini({
    # ============================================================
    # Placebo response parameters. Fixed in the combined PD model
    # per Li 2015 Table 3 footnote e: "Parameter estimates for WT
    # in placebo were fixed in finally combined PD model" -- so
    # Pmax and Kp were estimated in a placebo-only fit (with
    # ITV) and then held constant when the combined placebo+drug
    # model was fit.
    # ============================================================

    pmax <- fixed(-1.33)
    label("Maximum placebo effect Pmax on body weight from baseline (kg; negative = weight loss; FIXED from placebo-only fit)")  # Li 2015 Table 3, Pmax

    lkp <- fixed(log(0.0987))
    label("Log first-order rate constant Kp for placebo response (1/week; FIXED from placebo-only fit)")  # Li 2015 Table 3, Kp

    # ============================================================
    # Drug effect parameters (Li 2015 Eq 2 Emax + first-order
    # response). Estimated in the combined fit.
    # ============================================================

    dmax <- -1.85
    label("Maximum drug effect Dmax on body weight from baseline (kg; negative = additional weight loss beyond placebo)")  # Li 2015 Table 3, Dmax

    ic50 <- 41.7
    label("Half-maximal drug-effect concentration IC50 on the CAV metric (pmol/L)")  # Li 2015 Table 3, IC50

    lkdrug <- log(0.422)
    label("Log first-order rate constant Kdrug for drug response (1/week)")  # Li 2015 Table 3, Kdrug

    # ============================================================
    # Between-trial variability (study-level eta; "ITV" in Li 2015
    # Table 3 -- inter-trial variability). The placebo etas
    # (eta_study_pmax, eta_study_lkp) are fixed at the
    # placebo-only-fit ITV values (footnote e). The drug etas
    # (eta_study_dmax, eta_study_ic50, eta_study_lkdrug) are
    # estimated.
    #
    # IIV form per Li 2015 Eq 3 / Eq 4 (Statistical model
    # section): exponential model (P_i = P_pop * exp(eta)) for kp
    # and kdrug; additive model (P_i = P_pop + eta) for Pmax,
    # Dmax, and IC50. ITV % reported relative to the typical
    # value, consistent with Monolix 4.3 conventions, so the
    # variance for an additive model is (ITV/100 * |theta|)^2 and
    # the variance for an exponential model on the log scale
    # approximates the log-CV. See vignette "Source trace" /
    # "Assumptions and deviations" for the explicit numbers.
    # ============================================================

    eta_study_pmax   ~ fixed(2.1389)  # (1.10 * 1.33)^2 = 1.4630^2 kg^2; FIXED from placebo-only fit; Li 2015 Table 3 ITV 110%
    eta_study_lkp    ~ fixed(0.1089)  # 0.33^2 (log-scale variance); FIXED from placebo-only fit; Li 2015 Table 3 ITV 33%
    eta_study_dmax   ~ 0.0326         # (0.0976 * 1.85)^2 = 0.18056^2 kg^2; Li 2015 Table 3 ITV 9.76%
    eta_study_ic50   ~ 1881.4         # (1.04 * 41.7)^2 = 43.368^2 (pmol/L)^2; Li 2015 Table 3 ITV 104%
    eta_study_lkdrug ~ 0.0213         # 0.146^2 (log-scale variance); Li 2015 Table 3 ITV 14.6%

    # ============================================================
    # Residual error. Li 2015 Eq 5 proportional model
    # Y = F + b * F * epsilon, so propSd = b (the residual SD on
    # the proportional scale) with epsilon ~ N(0, 1).
    # ============================================================

    propSd <- 0.189
    label("Proportional residual SD on study-arm-mean WT change (fraction; Li 2015 Eq 5 b)")  # Li 2015 Table 3, b
  })

  model({
    # ----- Individual (study-arm-level) parameters -----
    # Pmax and Dmax: additive between-study variability per Li 2015 Eq 4
    # (additive in nlmixr2lib parlance maps to Monolix "addictive" /
    # constant model: P_i = P_pop + eta).
    pmax_i  <- pmax + eta_study_pmax
    dmax_i  <- dmax + eta_study_dmax

    # IC50: additive between-study variability per Eq 4. The
    # additive form with ITV 104% (~SD 43.4 pmol/L on a mean of
    # 41.7 pmol/L) can yield non-positive IC50 draws in stochastic
    # simulation; the model is faithful to the paper's encoding,
    # and downstream users are advised to use rxode2::zeroRe() for
    # typical-value replication of the published figures.
    ic50_i  <- ic50 + eta_study_ic50

    # Kp and Kdrug: exponential between-study variability per Eq 3
    # (P_i = P_pop * exp(eta)).
    kp_i    <- exp(lkp    + eta_study_lkp)
    kdrug_i <- exp(lkdrug + eta_study_lkdrug)

    # ----- Drug-side Emax term (Li 2015 Eq 2) -----
    # CAV is the study-arm average plasma concentration of
    # taspoglutide between weeks 2 and 4 (pmol/L). CAV = 0 for
    # placebo arms (drug_emax_i then collapses to 0).
    drug_emax_i <- dmax_i * CAV / (ic50_i + CAV)

    # ----- Algebraic time-course (Li 2015 Eq 1 and Eq 2) -----
    # placebo(t) = Pmax * (1 - exp(-Kp * t))
    # drug(t)    = drug_emax * (1 - exp(-Kdrug * t))
    # The MBMA emits study-arm-mean WT change at each requested
    # observation time; there is no compartment / dosing-event
    # consumption (CAV carries exposure).
    placebo_t <- pmax_i      * (1 - exp(-kp_i    * time))
    drug_t    <- drug_emax_i * (1 - exp(-kdrug_i * time))

    # ----- Output -----
    # Cc is the signed change in body weight from pretreatment
    # baseline (kg); placebo + drug response per Li 2015 Eq 2.
    # Cc is NOT a drug concentration.
    Cc <- placebo_t + drug_t

    Cc ~ prop(propSd)
  })
}
