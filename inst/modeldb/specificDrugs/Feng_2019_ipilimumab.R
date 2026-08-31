Feng_2019_ipilimumab <- function() {
  description <- "Three-subpopulation mixture tumor-growth-dynamics (TGD) model for tumor burden (sum of longest diameters of target lesions, cm) in patients with advanced melanoma treated with ipilimumab 3 or 10 mg/kg every 3 weeks (phase III CA184-169 / NCT01515189). Tumor burden follows a modified Wang-type algebraic profile TB(t) = TB0 * exp(-TS * t) + TG * t + TBss, where the three latent subpopulations differ in which terms are active and in their parameter values: fast tumor growth (TBss = 0), no growth (TG structurally 0, non-zero steady-state plateau TBss), and intermediate tumor growth and shrinkage (TBss = 0). Subpopulation membership is supplied by the user through the paired binary indicators MIX_FAST_GROW and MIX_NO_GROW (both 0 = intermediate, the reference class); the estimated population class probabilities are 39.0% fast, 32.5% no-growth and 28.5% intermediate. Baseline tumor burden carries a log-linear ULN-normalized baseline LDH effect and the linear growth rate carries both an LDH effect and an ipilimumab exposure effect driven by Cavg1, the time-averaged serum concentration after the first dose, supplied as the covariate CAV from an upstream population PK analysis (see modellib('Feng_2014_ipilimumab')). Interindividual variability is log-normal, with a correlated growth / shrinkage block within each of the fast and intermediate subpopulations. Residual error is the paper's single-epsilon combined form, SD(TB) = addSd + propSd * TB, encoded as combined1(). This file does NOT carry the paper's overall-survival model: that is a semiparametric Cox proportional-hazards regression whose nonparametric baseline hazard lambda0(t) is not reported (and is not reportable in closed form), so it cannot be simulated; it is described in the paired vignette narrative instead."
  reference <- paste(
    "Feng Y, Wang X, Suryawanshi S, Bello A, Roy A. Linking tumor growth",
    "dynamics to survival in ipilimumab-treated patients with advanced",
    "melanoma using mixture tumor growth dynamic modeling.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8(11):825-834.",
    "doi:10.1002/psp4.12454",
    sep = " "
  )
  vignette <- "Feng_2019_ipilimumab"
  units    <- list(time = "week", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    MIX_FAST_GROW = list(
      description        = "Fast-tumor-growth mixture subpopulation indicator: 1 = subject belongs to the fast tumor-growth subpopulation, 0 = subject belongs to the no-growth or intermediate subpopulation. Paired with MIX_NO_GROW; both 0 denotes the intermediate tumor-growth-and-shrinkage reference class. Exactly one of MIX_FAST_GROW / MIX_NO_GROW may be 1.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not fast-growth). The all-zero pattern MIX_FAST_GROW = 0, MIX_NO_GROW = 0 is the intermediate tumor-growth-and-shrinkage subpopulation.",
      notes              = "Latent class, not a measured covariate: in the source it is the per-subject NONMEM $MIXTURE assignment. Estimated population probability P(fast) = TP1 / (1 + TP1 + TP2) = 1.20 / 3.078 = 0.390 (Feng 2019 Table 2 footnote e). For population simulation draw the class per subject from the categorical distribution (fast, no-growth, intermediate) = (0.390, 0.325, 0.285). Feng 2019 Table S2 instead reports the empirical per-arm posterior class assignments -- 3 mg/kg 46.9% fast / 24.2% no-growth / 28.9% intermediate; 10 mg/kg 40.9% / 28.7% / 30.4% -- which differ from the population probabilities by a few percentage points because they are per-subject maximum-posterior assignments rather than the prior class mixture; see the vignette Errata.",
      source_name        = "MIXTURE (component 1 = fast TG)"
    ),
    MIX_NO_GROW = list(
      description        = "No-growth mixture subpopulation indicator: 1 = subject belongs to the no-growth subpopulation, in which the linear tumor-growth term is structurally zero and tumor burden decays from baseline to the non-zero steady-state plateau TBss; 0 = subject belongs to the fast-growth or intermediate subpopulation. Paired with MIX_FAST_GROW; both 0 denotes the intermediate reference class.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not no-growth). The all-zero pattern MIX_FAST_GROW = 0, MIX_NO_GROW = 0 is the intermediate tumor-growth-and-shrinkage subpopulation.",
      notes              = "Latent class, not a measured covariate. Estimated population probability P(no growth) = 1 / (1 + TP1 + TP2) = 1 / 3.078 = 0.325 (Feng 2019 Table 2 footnote e). This is the class the paper added to the pre-immunotherapy Wang model specifically to capture the immunotherapy response pattern in which tumor burden asymptotically approaches a durable steady-state value rather than always eventually growing (Feng 2019 Methods, 'TGD model'). Feng 2019 Table S2 empirical per-arm assignment: 24.2% (3 mg/kg) and 28.7% (10 mg/kg).",
      source_name        = "MIXTURE (component 2 = no growth)"
    ),
    LDHR = list(
      description        = "Baseline serum lactate dehydrogenase normalized to the reporting laboratory's upper limit of normal (unitless ratio). Enters both baseline tumor burden and the linear tumor-growth rate through the log-linear form (1 + log(LDHR) * coefficient), which is centered at exactly the ULN: a patient at LDHR = 1 gets a multiplier of exactly 1.",
      units              = "(unitless ratio)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Feng 2019 Table 1: median 1.0, range 0.4-40.5 across the whole 688-patient analysis population (3 mg/kg 1.0, range 0.4-29.4; 10 mg/kg 1.0, range 0.5-40.5). Table 1 footnote a: 'LDH ratio indicates patient's actual value divided by the upper limit of normal. Log-transformed LDH ratio was used in tumor growth dynamics and overall survival model development due to skewed distribution.' Note this is the ULN-normalized ratio (canonical LDHR), NOT the absolute U/L activity (canonical LDH) used by the companion ipilimumab popPK models Feng_2014_ipilimumab.R and Sanghavi_2020_ipilimumab.R -- feeding a U/L value here would be off by roughly two orders of magnitude.",
      source_name        = "BLDHU"
    ),
    CAV = list(
      description        = "Cavg1, the time-averaged ipilimumab serum concentration over the interval after the first dose, in ug/mL. Time-fixed per subject (a single first-dose exposure summary, not a running average). Supplied by the user from an upstream population PK analysis; drives the exposure effect on the linear tumor-growth rate.",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Feng 2019 Methods: 'Patient-specific values of ipilimumab CL and Cavg1 were obtained from a population pharmacokinetics analysis' (reference 22, Vezina 2018, a JPKPD conference abstract not available as a full model). IMPORTANT -- Feng 2019 never states the units of Cavg1 anywhere in the text, tables, or supplement. ug/mL is inferred and is the only dimensionally feasible reading: the Table 2 coefficient is -0.00342 per unit and enters as the linear multiplier (1 + CAV * -0.00342), which stays positive only for CAV < 292. Simulating a typical 80 kg patient through modellib('Feng_2014_ipilimumab') and averaging the concentration over the first 3-week dosing interval gives Cavg1 = 20.6 ug/mL at 3 mg/kg and 68.6 ug/mL at 10 mg/kg, hence growth multipliers of 0.930 and 0.765 -- both positive, and lower at the higher dose, matching the paper's finding of a lower tumor growth rate with 10 mg/kg. In ng/mL the same exposures would give a multiplier of about -234 (a negative tumor-growth rate of absurd magnitude); in mg/mL the covariate effect would vanish to four significant figures, contradicting its retention in the final model. The paper does not define the averaging window precisely either, but the conclusion is insensitive to that: the steady-state-equivalent Dose / (CL * tau) reading gives 31.7 and 105.8 ug/mL, which is the same order of magnitude and the same qualitative result. ug/mL is also the concentration unit used by both registry ipilimumab popPK models from the same Bristol-Myers Squibb group. Set CAV = 0 to switch the exposure effect off. See the vignette Errata.",
      source_name        = "Cavg1"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 688L,
    n_studies       = 1L,
    age_mean        = "59.9 years (SD 14.0) overall; 60.9 (SD 13.3) in the 3 mg/kg arm and 59.0 (SD 14.6) in the 10 mg/kg arm",
    weight_mean     = "79.9 kg (SD 17.6) overall; 79.3 (SD 17.4) at 3 mg/kg and 80.6 (SD 17.7) at 10 mg/kg",
    sex_female_pct  = 37.4,
    disease_state   = "Previously treated or untreated unresectable or metastatic (advanced) melanoma. M stage: M0 or M1a 17.0%, M1b 20.9%, M1c 62.1%. ECOG performance status 0 in 70.6% and >= 1 in 29.5%.",
    dose_range      = "Ipilimumab 3 mg/kg (n = 343) or 10 mg/kg (n = 345) by 90-minute intravenous infusion every 3 weeks for four doses.",
    regions         = "Multicentre international phase III study CA184-169 (NCT01515189).",
    sampling_window = "Protocol-specified tumor assessments at weeks 12, 16 and 24; tumor burden measured as the sum of the longest diameters of target lesions under immune-related response criteria.",
    tumor_type      = "Cutaneous / advanced melanoma. Mean baseline tumor burden 9.6 cm (SD 8.7) overall; 10.0 cm (SD 8.5) at 3 mg/kg and 9.3 cm (SD 9.0) at 10 mg/kg.",
    assay           = "Tumor burden = sum of the longest diameters of target lesions (cm) by immune-related response criteria.",
    notes           = "Analysis population is the 688 of the CA184-169 randomised patients for whom tumor burden data were available. Baseline LDH ratio to ULN median 1.0 (range 0.4-40.5). The paper's final tumor-dynamics model is TGD-Model 11 (mixture with three subpopulations, Cavg1 effect on TG, LDH effect on TG and TB0; BIC 6032.321, the reference model in Table S1), and it is the model encoded here. The companion overall-survival analysis is a semiparametric Cox proportional-hazards model (OS-Model 1 / S3) that is NOT encoded in this file because its nonparametric baseline hazard is not reported."
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural typical values, Feng 2019 Table 2 "Estimate" column.
    #
    # The three subpopulations of the mixture are strata of a single joint
    # fit, so every stratum-specific quantity carries an explicit suffix
    # (_fast / _nogrow / _inter) and none keeps the bare canonical name --
    # see references/parameter-names.md "Stratum-suffixed parameters".
    #
    # All typical values are strictly positive and enter multiplicatively
    # with exponential IIV, so they are stored on the log scale.
    # ---------------------------------------------------------------------

    # Fast tumor-growth subpopulation (mixture component 1).
    ltb0_fast   <- log(10.2)    ; label("Typical baseline tumor burden TB0, fast-growth subpopulation (cm)")            # Feng 2019 Table 2, Fast: TB0 = 10.2 cm (95% CI 9.04-11.4)
    ltg_fast    <- log(0.328)   ; label("Typical linear tumor-growth rate TG, fast-growth subpopulation (cm/week)")     # Feng 2019 Table 2, Fast: TG = 0.328 cm/week (95% CI 0.252-0.405)
    lts_fast    <- log(0.00362) ; label("Typical tumor-shrinkage rate constant TS, fast-growth subpopulation (1/week)") # Feng 2019 Table 2, Fast: TS = 0.00362 /week (95% CI -5.95E-04 to 0.00783)

    # No-growth subpopulation (mixture component 2). TG is structurally zero
    # here -- Feng 2019 Methods: "Each of the mixture models included one
    # subpopulation to describe the no-growth subpopulation (TG fixed to
    # zero)" -- and the steady-state plateau TBss is estimated in its place.
    ltb0_nogrow  <- log(2.53)   ; label("Typical baseline tumor burden TB0, no-growth subpopulation (cm)")              # Feng 2019 Table 2, No growth: TB0 = 2.53 cm (95% CI 2.05-3.00)
    lts_nogrow   <- log(0.0458) ; label("Typical tumor-shrinkage rate constant TS, no-growth subpopulation (1/week)")   # Feng 2019 Table 2, No growth: TS = 0.0458 /week (95% CI 0.0340-0.0576)
    ltbss_nogrow <- log(1.71)   ; label("Typical steady-state tumor burden TBss, no-growth subpopulation (cm)")         # Feng 2019 Table 2, No growth: TBss = 1.71 cm (95% CI 1.15-2.26)

    # Intermediate tumor-growth-and-shrinkage subpopulation (mixture
    # component 3; the reference class when both MIX_ indicators are 0).
    # Feng 2019 labels these rows "TB0 P3", "TG P3" and "TS P3".
    ltb0_inter <- log(5.83)    ; label("Typical baseline tumor burden TB0, intermediate subpopulation (cm)")            # Feng 2019 Table 2, Intermediate: TB0 P3 = 5.83 cm (95% CI 4.69-6.98)
    ltg_inter  <- log(0.0236)  ; label("Typical linear tumor-growth rate TG, intermediate subpopulation (cm/week)")     # Feng 2019 Table 2, Intermediate: TG P3 = 0.0236 cm/week (95% CI 0.00746-0.0398)
    lts_inter  <- log(0.00299) ; label("Typical tumor-shrinkage rate constant TS, intermediate subpopulation (1/week)") # Feng 2019 Table 2, Intermediate: TS P3 = 0.00299 /week (95% CI 7.56E-04 to 0.00522)

    # Mixture-proportion parameters. Feng 2019 Table 2 footnote e: the
    # overall class probabilities are TP1 / (1 + TP1 + TP2) for mixture 1
    # (fast TG), 1 / (1 + TP1 + TP2) for mixture 2 (no growth), and
    # TP2 / (1 + TP1 + TP2) for mixture 3 (intermediate TG and TS).
    # These do not enter the tumor-burden prediction -- class membership is
    # supplied through MIX_FAST_GROW / MIX_NO_GROW -- but they are estimated
    # parameters of the published model and are carried here so a user can
    # draw the latent class. model() exposes them as p_fast / p_nogrow /
    # p_inter.
    tp1 <- 1.20  ; label("Mixture proportion parameter TP1 (fast-growth class weight, unitless)")   # Feng 2019 Table 2: TP1 = 1.20 (95% CI 0.719-1.67)
    tp2 <- 0.878 ; label("Mixture proportion parameter TP2 (intermediate class weight, unitless)")  # Feng 2019 Table 2: TP2 = 0.878 (95% CI 0.456-1.30)

    # ---------------------------------------------------------------------
    # Covariate effects, Feng 2019 Table 2 footnote f:
    #   TG_TV  = TG_REF  * (1 + log(BLDHU) * TG_BLDHU) * (1 + Cavg1 * TG_Cavg1)
    #   TB0_TV = TB0_REF * (1 + log(BLDHU) * TB0_BLDHU)
    # Note these are LINEAR-in-coefficient multipliers, not the exponential
    # exp(coef * log(cov)) form used elsewhere in the registry.
    #
    # ERRATUM: the printed footnote f writes the symbol "TGBLDHU" in BOTH
    # equations, i.e. it reuses the TG coefficient in the TB0 equation. That
    # is a typographical slip: Table 2 tabulates two distinct estimates with
    # distinct confidence intervals ("LDH effect on TB0" = 0.868 and
    # "LDH effect on TG" = 0.771), and the Results text describes them as
    # two separate findings ("Estimated LDH effects on TB0 and TG showed
    # that higher baseline LDH was associated with a higher TG rate and
    # higher baseline tumor size"). The TB0 equation therefore uses the TB0
    # coefficient. See the vignette Errata.
    # ---------------------------------------------------------------------
    e_ldhr_tb0 <- 0.868    ; label("Effect of log(LDHR) on baseline tumor burden TB0 (linear multiplier coefficient)")            # Feng 2019 Table 2: LDH effect on TB0 = 0.868 (95% CI 0.752-0.984)
    e_ldhr_tg  <- 0.771    ; label("Effect of log(LDHR) on linear tumor-growth rate TG (linear multiplier coefficient)")          # Feng 2019 Table 2: LDH effect on TG = 0.771 (95% CI 0.473-1.07)
    e_cav_tg   <- -0.00342 ; label("Effect of Cavg1 on linear tumor-growth rate TG (linear multiplier coefficient, per ug/mL)")   # Feng 2019 Table 2: Exposure (Cavg1) effect on TG = -0.00342 (95% CI -0.00690 to 6.69E-05)

    # ---------------------------------------------------------------------
    # Interindividual variability. Feng 2019 Methods: P_i = P_TV * exp(eta),
    # i.e. log-normal on every structural parameter. Table 2 reports the
    # random effects as variance with SD in parentheses on the diagonal and
    # as covariance with correlation in parentheses off-diagonal (footnote
    # c), so the numbers entered below are variances / covariances.
    #
    # Every reported SD and correlation back-computes exactly from the
    # variances, which confirms the transcription:
    #   sqrt(0.535) = 0.731  sqrt(0.360) = 0.600  sqrt(4.07)  = 2.02
    #   sqrt(0.385) = 0.621  sqrt(1.38)  = 1.17   sqrt(0.203) = 0.451
    #   sqrt(3.21)  = 1.79
    #   -1.06 / sqrt(0.360 * 4.07) = -0.876 (Table 2 reports -0.878)
    #    0.129 / sqrt(0.203 * 3.21) =  0.160 (Table 2 reports  0.159)
    #
    # TB0 variability is listed in Table 2 above the subpopulation
    # sub-headings and is therefore shared across all three classes; it
    # keeps the bare canonical name. The growth / shrinkage variabilities
    # are listed under their subpopulation sub-headings and are suffixed.
    # ---------------------------------------------------------------------
    etaltb0 ~ 0.535                                            # Feng 2019 Table 2: omega^2 TB0 = 0.535 (SD 0.731), 95% CI 0.451-0.619

    # Fast subpopulation: correlated TG / TS block.
    etaltg_fast + etalts_fast ~ c(0.360,
                                  -1.06, 4.07)                 # Feng 2019 Table 2, Fast: omega^2 TG = 0.360 (SD 0.600); cov(TG,TS) = -1.06 (corr -0.878); omega^2 TS = 4.07 (SD 2.02)

    # No-growth subpopulation: TS and TBss, no covariance reported.
    etalts_nogrow   ~ 0.385                                    # Feng 2019 Table 2, No growth: omega^2 TS   = 0.385 (SD 0.621), 95% CI 0.135-0.636
    etaltbss_nogrow ~ 1.38                                     # Feng 2019 Table 2, No growth: omega^2 TBss = 1.38  (SD 1.17),  95% CI 0.857-1.90

    # Intermediate subpopulation: correlated TG / TS block.
    etaltg_inter + etalts_inter ~ c(0.203,
                                    0.129, 3.21)               # Feng 2019 Table 2, Intermediate: omega^2 TG = 0.203 (SD 0.451); cov(TG,TS) = 0.129 (corr 0.159); omega^2 TS = 3.21 (SD 1.79)

    # ---------------------------------------------------------------------
    # Residual error. Feng 2019 Methods:
    #   TB_obs(t) = TB_pred(t) * (1 + theta_PROP * eps) + theta_ADD * eps
    # with a SINGLE eps ~ N(0, 1). Factoring the shared eps gives
    #   SD(TB_obs) = theta_ADD + theta_PROP * TB_pred,
    # i.e. the additive and proportional components add on the SD scale
    # rather than in quadrature. That is exactly rxode2's combined1()
    # variance model, so the error is declared as
    #   add(addSd) + prop(propSd) + combined1()
    # NOT the default add + prop (which would be combined2(), summing the
    # two components in quadrature and understating the SD).
    # ---------------------------------------------------------------------
    addSd  <- 0.125 ; label("Additive residual error SD (cm)")            # Feng 2019 Table 2: Additive error = 0.125 cm (95% CI 0.0927-0.158)
    propSd <- 0.167 ; label("Proportional residual error SD (fraction)")  # Feng 2019 Table 2: Proportional error = 0.167 (95% CI 0.159-0.175)
  })

  model({
    # Intermediate tumor-growth-and-shrinkage subpopulation is the reference
    # class: both mixture indicators zero.
    mix_inter <- 1 - MIX_FAST_GROW - MIX_NO_GROW

    # Estimated population class probabilities (Feng 2019 Table 2 footnote e).
    # Carried as derived quantities so the TP1 / TP2 estimates stay attached
    # to the model; they do not feed the tumor-burden prediction, which uses
    # the user-supplied class indicators instead.
    p_fast   <- tp1 / (1 + tp1 + tp2)
    p_nogrow <- 1   / (1 + tp1 + tp2)
    p_inter  <- tp2 / (1 + tp1 + tp2)

    # Subpopulation-specific typical values. Exactly one of the three gates
    # is 1 for any given subject, so each of these selects a single anchor.
    # TG has no no-growth term (structurally zero) and TBss exists only in
    # the no-growth class.
    tb0_tv  <- exp(ltb0_fast) * MIX_FAST_GROW + exp(ltb0_nogrow) * MIX_NO_GROW + exp(ltb0_inter) * mix_inter
    ts_tv   <- exp(lts_fast)  * MIX_FAST_GROW + exp(lts_nogrow)  * MIX_NO_GROW + exp(lts_inter)  * mix_inter
    tg_tv   <- exp(ltg_fast)  * MIX_FAST_GROW                                  + exp(ltg_inter)  * mix_inter
    tbss_tv <- exp(ltbss_nogrow) * MIX_NO_GROW

    # Subpopulation-specific interindividual variability. Feng 2019 estimates
    # separate TG / TS variance components per class (Table 2), so the eta
    # that is active for a subject is the one belonging to that subject's
    # class. The etas of the inactive classes are still drawn but are gated
    # out of the prediction.
    eta_ts   <- etalts_fast * MIX_FAST_GROW + etalts_nogrow * MIX_NO_GROW + etalts_inter * mix_inter
    eta_tg   <- etaltg_fast * MIX_FAST_GROW                               + etaltg_inter * mix_inter
    eta_tbss <- etaltbss_nogrow * MIX_NO_GROW

    # Covariate multipliers (Feng 2019 Table 2 footnote f). Both are linear
    # in the coefficient. LDHR = 1 (patient at the laboratory ULN) and
    # CAV = 0 (no exposure) each give a multiplier of exactly 1.
    ldh_eff_tb0 <- 1 + log(LDHR) * e_ldhr_tb0
    ldh_eff_tg  <- 1 + log(LDHR) * e_ldhr_tg
    cav_eff_tg  <- 1 + CAV * e_cav_tg

    # Individual parameters: typical value * covariate multipliers * exp(eta).
    tb0  <- tb0_tv  * ldh_eff_tb0 * exp(etaltb0)
    ts   <- ts_tv   * exp(eta_ts)
    tg   <- tg_tv   * ldh_eff_tg * cav_eff_tg * exp(eta_tg)
    tbss <- tbss_tv * exp(eta_tbss)

    # Algebraic tumor-burden profile (Feng 2019 Methods, "TGD model", the
    # three per-subpopulation structural equations written as one expression
    # -- for the fast and intermediate classes tbss is 0 and the equation is
    # TB0 * exp(-TS * t) + TG * t; for the no-growth class tg is 0 and the
    # equation is TB0 * exp(-TS * t) + TBss). `t` is time since the first
    # ipilimumab dose, in weeks.
    tumor_size <- tb0 * exp(-ts * t) + tg * t + tbss

    # Single-epsilon combined additive + proportional residual error; see the
    # derivation in ini().
    tumor_size ~ add(addSd) + prop(propSd) + combined1()
  })
}
