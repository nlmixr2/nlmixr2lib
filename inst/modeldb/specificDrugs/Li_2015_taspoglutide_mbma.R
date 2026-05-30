Li_2015_taspoglutide_mbma <- function() {
  description <- paste0(
    "MBMA. Coupled PD model-based meta-analysis of taspoglutide (long-acting ",
    "human glucagon-like peptide-1 analogue, once-weekly SC) net efficacy on ",
    "fasting plasma glucose (FPG) and glycosylated hemoglobin (HbA1c) in type ",
    "2 diabetes. Each endpoint is the sum of an exponential-to-asymptote ",
    "placebo response (Pmax, Kp) and a saturable Emax drug response (Dmax, ",
    "IC50, Kdrug) approached exponentially over time. The FPG drug effect is ",
    "driven by the study-arm-mean taspoglutide concentration between weeks 2 ",
    "and 4 (Cavg; supplied as the METRIC_TASPO_C covariate: 0 / 59.85 / 119.7 ",
    "pmol/L for placebo / 10 mg / 20 mg QW). The HbA1c drug effect is driven ",
    "by the model-predicted drug-induced FPG change (i.e. the placebo-adjusted ",
    "FPG response feeds the HbA1c Emax). Estimated on digitised study-arm-mean ",
    "PD data from 8 published clinical trials of taspoglutide monotherapy or ",
    "add-on therapy in type 2 diabetes (3,702 patients pooled, 8-52 week ",
    "treatment durations); a ninth trial (Rosenstock 2013) was held out for ",
    "external validation. Placebo Pmax and Kp were fitted on the placebo-only ",
    "subset first and held fixed in the final combined PD model. Between-trial ",
    "variability (ITV) is encoded as study-level etas (one eta per parameter); ",
    "the model is suitable for simulating study-arm-mean PD outcomes and is ",
    "NOT suitable for individual-subject simulation. Residual error is a ",
    "proportional/power model on each endpoint (the small power-correction ",
    "term is simplified to a plain proportional error in this implementation; ",
    "see vignette Assumptions and deviations)."
  )

  reference <- paste(
    "Li HQ, Xu JY, Jin L, Xin JL.",
    "Utilization of model-based meta-analysis to delineate the net efficacy",
    "of taspoglutide from the response of placebo in clinical trials.",
    "Saudi Pharm J. 2015 Jul;23(3):241-249.",
    "doi:10.1016/j.jsps.2014.11.008.",
    sep = " "
  )
  vignette <- "Li_2015_taspoglutide_mbma"
  units <- list(
    time          = "week",
    dosing        = "n/a (MBMA driver is the METRIC_TASPO_C covariate)",
    concentration = paste0(
      "FPG (mmol/L change from baseline; signed); HbA1c (% change from ",
      "baseline; signed). Output Cc is unused."
    )
  )

  covariateData <- list(
    METRIC_TASPO_C = list(
      description        = paste0(
        "Per-arm taspoglutide average plasma concentration between weeks 2 ",
        "and 4 (the drug-effect 'metric' in Li 2015). Set to 0 for placebo ",
        "arms; 59.85 pmol/L for taspoglutide 10 mg once-weekly SC arms; ",
        "119.7 pmol/L for taspoglutide 20 mg once-weekly SC arms."
      ),
      units              = "pmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "MBMA study-arm-level covariate. The canonical register in ",
        "inst/references/covariate-columns.md is for individual-level ",
        "pop-PK covariates and does not directly fit this MBMA study-arm ",
        "drug-exposure column; this column mirrors the inline-documented ",
        "precedent set by Vargo_2014_statins_ezetimibe_mbma. The values ",
        "are reported by Li 2015 Section 3.2: 'the average concentrations ",
        "(0, 59.85 and 119.7 pmol/l for placebo, taspoglutide 10 and 20 ",
        "mg) between 2nd and 4th weeks were chosen as metrics'. The ",
        "concentrations themselves were derived by Li et al. from Ratner ",
        "2010 (the source taspoglutide PK study)."
      ),
      source_name        = "Cavg.2-4w (Li 2015 Section 3.2)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 3702L,
    n_studies      = 8L,
    n_studies_validation = 1L,
    age_range      = paste0(
      "adults with type 2 diabetes (per-study demographics not pooled in ",
      "Li 2015 Table 1; underlying trials enrolled adults broadly per ",
      "Ratner 2010, Nauck 2009, Bergenstal 2012, Henry 2012, Raz 2012, ",
      "Hollander 2013, Pratley 2013, Nauck 2013)"
    ),
    age_median     = NA_character_,
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = paste0(
      "Type 2 diabetes mellitus, pooled across mono- and combination-therapy ",
      "regimens (drug-naive; inadequately controlled on metformin; failing ",
      "metformin + sulphonylurea; inadequately controlled on metformin + ",
      "TZD; etc., per Li 2015 Table 1)"
    ),
    dose_range     = paste0(
      "Placebo, taspoglutide 10 mg once weekly SC, taspoglutide 20 mg once ",
      "weekly SC; treatment durations 8-52 weeks (Li 2015 Table 1)"
    ),
    regions        = "International (publicly published phase 2/3 trials)",
    notes          = paste0(
      "Model fit by Monolix 4.2/4.3 (SAEM) to digitised mean PD data from 9 ",
      "publications: Nauck 2009 (T-emerge 1 setup, n=373), Ratner 2010 ",
      "(n=64), Nauck 2009 (n=148), Nauck 2013 (n=709), Hollander 2013 ",
      "T-emerge 7 (n=305), Bergenstal 2012 T-emerge 2 (n=481), Henry 2012 ",
      "T-emerge 3 (n=326), Pratley 2013 T-emerge 6 (n=499), Rosenstock 2013 ",
      "T-emerge 4 (n=797). Eight trials (n=3702) were used for model ",
      "development; Rosenstock 2013 (n=797) was the external validation ",
      "cohort. Mean PD response (change from baseline) is the unit of ",
      "observation, NOT individual subject concentration. See Li 2015 Table ",
      "1 for the per-study breakdown."
    )
  )

  # ============================================================================
  # ini(): structural parameters and ITV from Li 2015 Tables 2 (FPG) and 3
  # (HbA1c). Sign convention follows the paper: Pmax and Dmax are signed
  # negative numbers (reductions from baseline); IC50_F is positive pmol/L
  # taspoglutide concentration; IC50_Hb is signed negative mmol/L (the
  # half-maximal HbA1c efficacy point in FPG-reduction space). Kp and Kdrug
  # are positive rate constants (1/week) and are log-transformed.
  #
  # ITV interpretation (paper Section 2.6 Eq. 3-4; Monolix output):
  #   - Exponential model (Kp, Kdrug): P_i = P_pop * exp(eta), eta ~ N(0, omega^2).
  #     The table's ITV(%) is read as approximate CV%, so
  #     omega^2 = log((ITV/100)^2 + 1).
  #   - Additive model (Pmax, Dmax, IC50): P_i = P_pop + eta, eta ~ N(0, omega^2).
  #     The table's ITV(%) is read as relative SD against the typical value,
  #     so omega = (ITV/100) * |P_pop| and omega^2 = ((ITV/100)*|P_pop|)^2.
  # The paper does not state the units of the ITV column; this convention is
  # the most common Monolix reading. The interpretation is called out in the
  # vignette Assumptions and deviations section.
  #
  # Placebo parameters (Pmax_F, Kp_F, Pmax_Hb, Kp_Hb) were fitted on the
  # placebo-only data first and then fixed in the final combined PD model
  # (Tables 2 and 3 footnote: 'Parameter estimates for FPG/HbA1c in placebo
  # were fixed in finally combined PD model'). They are wrapped in fixed()
  # together with their ITV.
  # ============================================================================
  ini({
    # ----- FPG endpoint: placebo structural parameters (fixed in final model) -----
    pmax_f <- fixed(-0.371)
    label("Placebo Pmax on FPG (mmol/L; signed, negative = reduction)")  # Li 2015 Table 2 (fixed)

    lkp_f <- fixed(log(0.781))
    label("Log of placebo Kp on FPG (1/week)")  # Li 2015 Table 2 (Kp_F = 0.781 /week; fixed)

    # ----- FPG endpoint: drug structural parameters -----
    dmax_f <- -2.39
    label("Drug Dmax on FPG (mmol/L; signed, negative = reduction)")  # Li 2015 Table 2 (RSE 6%)

    ic50_f <- 25.3
    label("Drug IC50 on FPG (pmol/L taspoglutide; positive)")  # Li 2015 Table 2 (RSE reported as 0%, likely <0.5% rounded)

    lkdrug_f <- log(2.0)
    label("Log of drug Kdrug on FPG (1/week)")  # Li 2015 Table 2 (Kdrug_F = 2.0 /week; RSE 2%)

    # ----- HbA1c endpoint: placebo structural parameters (fixed in final model) -----
    pmax_hb <- fixed(-0.253)
    label("Placebo Pmax on HbA1c (%; signed, negative = reduction)")  # Li 2015 Table 3 (fixed)

    lkp_hb <- fixed(log(0.382))
    label("Log of placebo Kp on HbA1c (1/week)")  # Li 2015 Table 3 (Kp_Hb = 0.382 /week; fixed)

    # ----- HbA1c endpoint: drug structural parameters -----
    dmax_hb <- -1.74
    label("Drug Dmax on HbA1c (%; signed, negative = reduction)")  # Li 2015 Table 3 (RSE 8%)

    ic50_hb <- -1.81
    label("Drug IC50 on HbA1c (mmol/L FPG reduction; signed negative)")  # Li 2015 Table 3 (RSE 15%; IC50_Hb is the FPG-reduction value that produces 50% of Dmax_Hb)

    lkdrug_hb <- log(0.249)
    label("Log of drug Kdrug on HbA1c (1/week)")  # Li 2015 Table 3 (Kdrug_Hb = 0.249 /week; RSE 13%)

    # ----- Between-trial variability (ITV); diagonal -----
    # FPG endpoint
    eta_study_pmax_f   ~ fixed(0.012244)  # Li 2015 Table 2; additive, (0.298 * 0.371)^2; fixed with structural
    eta_study_lkp_f    ~ fixed(0.133593)  # Li 2015 Table 2; exponential, log(0.378^2 + 1); fixed with structural
    eta_study_dmax_f   ~ 0.351272         # Li 2015 Table 2; additive, (0.248 * 2.39)^2
    eta_study_ic50_f   ~ 1.887524         # Li 2015 Table 2; additive, (0.0543 * 25.3)^2
    eta_study_lkdrug_f ~ 0.015504         # Li 2015 Table 2; exponential, log(0.125^2 + 1)

    # HbA1c endpoint
    eta_study_pmax_hb   ~ fixed(0.001479)  # Li 2015 Table 3; additive, (0.152 * 0.253)^2; fixed with structural
    eta_study_lkp_hb    ~ fixed(0.021365)  # Li 2015 Table 3; exponential, log(0.147^2 + 1); fixed with structural
    eta_study_dmax_hb   ~ 0.022652         # Li 2015 Table 3; additive, (0.0865 * 1.74)^2
    eta_study_ic50_hb   ~ 0.235303         # Li 2015 Table 3; additive, (0.268 * 1.81)^2
    eta_study_lkdrug_hb ~ 0.038841         # Li 2015 Table 3; exponential, log(0.199^2 + 1)

    # ----- Residual error (Eq. 5 and 6; Tables 2 and 3) -----
    # Paper Eq. 5 (FPG): Y = F + b * F * eps, but Table 2 reports both b and
    # c, consistent with a Monolix 'power' error model Y = F + b * |F|^c * eps.
    # Paper Eq. 6 (HbA1c): Y = F + (a + b * F^c) * eps (Monolix 'combined-power').
    # nlmixr2 has no clean power-of-prediction syntax, so this implementation
    # uses prop() on FPG and add()+prop() on HbA1c with the 'b' coefficient as
    # propSd and the 'a' coefficient as addSd; the 'c' power exponent is not
    # encoded (see vignette Assumptions and deviations).
    propSd_FPG   <- 0.194
    label("Proportional residual SD on FPG (fraction; Li 2015 Table 2 'b')")

    addSd_HbA1c  <- 0.00532
    label("Additive residual SD on HbA1c (%; Li 2015 Table 3 'a')")
    propSd_HbA1c <- 0.0717
    label("Proportional residual SD on HbA1c (fraction; Li 2015 Table 3 'b')")
  })

  model({
    # ---------------------------------------------------------------------
    # MBMA study-arm-level driver: the study-arm-mean taspoglutide average
    # concentration between weeks 2 and 4 (Cavg) is supplied as a covariate
    # and is treated as a time-constant scalar over the simulation horizon
    # (the MBMA simplification used by Li 2015). Placebo arms supply 0; the
    # 10 mg QW arm supplies 59.85 pmol/L; the 20 mg QW arm supplies 119.7
    # pmol/L (Li 2015 Section 3.2).
    # ---------------------------------------------------------------------
    metric_c <- METRIC_TASPO_C

    # ---------------------------------------------------------------------
    # Trial-level (ITV-perturbed) structural parameters. Exponential model
    # for Kp and Kdrug; additive model for Pmax, Dmax, IC50.
    # ---------------------------------------------------------------------
    pmax_f_i    <- pmax_f    + eta_study_pmax_f
    kp_f_i      <- exp(lkp_f + eta_study_lkp_f)
    dmax_f_i    <- dmax_f    + eta_study_dmax_f
    ic50_f_i    <- ic50_f    + eta_study_ic50_f
    kdrug_f_i   <- exp(lkdrug_f + eta_study_lkdrug_f)

    pmax_hb_i   <- pmax_hb   + eta_study_pmax_hb
    kp_hb_i     <- exp(lkp_hb + eta_study_lkp_hb)
    dmax_hb_i   <- dmax_hb   + eta_study_dmax_hb
    ic50_hb_i   <- ic50_hb   + eta_study_ic50_hb
    kdrug_hb_i  <- exp(lkdrug_hb + eta_study_lkdrug_hb)

    # ---------------------------------------------------------------------
    # FPG endpoint (Li 2015 Eq. 1-2): placebo and drug compartments each
    # approach an Emax-defined asymptote exponentially. At t = 0 both states
    # start at zero (change-from-baseline scale).
    #
    #   Placebo asymptote = pmax_f
    #   Drug    asymptote = dmax_f * metric_c / (ic50_f + metric_c)
    # ---------------------------------------------------------------------
    fpg_drug_asymp <- dmax_f_i * metric_c / (ic50_f_i + metric_c)

    d/dt(fpg_placebo) <- kp_f_i    * (pmax_f_i        - fpg_placebo)
    d/dt(fpg_drug)    <- kdrug_f_i * (fpg_drug_asymp  - fpg_drug)

    # ---------------------------------------------------------------------
    # HbA1c endpoint (Li 2015 Eq. 1-2 applied to HbA1c with the drug-induced
    # FPG change as the metric -- Li 2015 Section 3.2: 'The metrics for HbA1c
    # data modeling were derived from the population prediction values for
    # taspoglutide 10 and 20 mg in FPG modeling. However, the metric for
    # placebo HbA1c data modeling were set to zero...').
    #
    # The metric is fpg_drug only (NOT fpg_placebo + fpg_drug); the placebo
    # HbA1c response is a separate first-order rise to its own asymptote.
    # When fpg_drug = 0 (placebo arm or t = 0) the Emax expression is
    # 0 / ic50_hb = 0 so no singularity arises in the physical operating
    # range.
    # ---------------------------------------------------------------------
    hba1c_drug_asymp <- dmax_hb_i * fpg_drug / (ic50_hb_i + fpg_drug)

    d/dt(hba1c_placebo) <- kp_hb_i    * (pmax_hb_i        - hba1c_placebo)
    d/dt(hba1c_drug)    <- kdrug_hb_i * (hba1c_drug_asymp - hba1c_drug)

    # ---------------------------------------------------------------------
    # Outputs: change-from-baseline on each endpoint = placebo response +
    # drug response (Li 2015 Eq. 2).
    # ---------------------------------------------------------------------
    FPG   <- fpg_placebo   + fpg_drug
    HbA1c <- hba1c_placebo + hba1c_drug

    FPG   ~ prop(propSd_FPG)
    HbA1c ~ add(addSd_HbA1c) + prop(propSd_HbA1c)
  })
}
