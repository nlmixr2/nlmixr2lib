Yin_2021_pexidartinib <- function() {
  description <- "Semi-mechanistic longitudinal tumor-size (RECIST) exposure-response PD model for pexidartinib in adult patients with tenosynovial giant cell tumor (TGCT), driven by pre-computed running average pexidartinib concentration Cavg (mg/L). Tumor size Y(t) declines from an individual baseline Y0 through a saturable drug effect gated by an onset first-order process: Y(t) = Y0 * (1 - Emax * (1 - exp(-kdrug * Cavg)) * (1 - exp(-konset * TAFD))) + growth * time, with Emax fixed at 0.999 and natural tumor growth fixed at 0 (Yin 2021 Results: the placebo-cohort growth rate estimate 0.227 cm/yr had a 95% CI that included zero and was fixed). Baseline Y0, drug-effect rate constant kdrug, and onset rate constant konset each carry three multiplicative covariate effects: joint extremity (upper vs lower reference), joint size (small vs large reference), and age centered at 44 years. Individual variability is a 3x3 log-normal block on baseline, kdrug, and konset. Residual error is a power-of-prediction form SD(Y) = 0.365 * Ypred^0.550. The pexidartinib PK backbone (Cavg input) is a separately extracted model, `modellib('Yin_2020_pexidartinib')`, corresponding to the same authors' Yin 2020 J Clin Pharmacol popPK publication (reference 5 in Yin 2021). This PD model does NOT carry the ORR proportional-odds logistic regression (Table S3) or the piecewise-exponential TTE hepatic-AE models (Table S5) that appear as parallel endpoints in Yin 2021; those are non-ODE statistical / survival regressions and are documented in the paired vignette narrative but not encoded here as separate model files."
  reference <- paste(
    "Yin O, Zahir H, French J, et al. Exposure-response analysis of efficacy",
    "and safety for pexidartinib in patients with tenosynovial giant cell",
    "tumor. CPT Pharmacometrics Syst Pharmacol. 2021;10(11):1422-1432.",
    "doi:10.1002/psp4.12712. PK backbone (Cavg input) adapted from Yin O,",
    "Wagner AJ, Kang J, et al. Population pharmacokinetic analysis of",
    "pexidartinib in healthy subjects and patients with tenosynovial giant",
    "cell tumor or other solid tumors. J Clin Pharmacol. 2020;61(4):480-492.",
    "doi:10.1002/jcph.1753; see modellib('Yin_2020_pexidartinib').",
    sep = " "
  )
  vignette <- "Yin_2021_pexidartinib"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CAV = list(
      description        = "Running average pexidartinib plasma concentration up to the current tumor-measurement time (Yin 2021 Methods: 'Cavg_ij is average concentration up to tumor measurement time j'). Time-varying per record; supplied by the user from an upstream individual-PK simulation (e.g., empirical-Bayes predictions or forward simulation of `modellib('Yin_2020_pexidartinib')` followed by cumulative time-averaging).",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Kdrug is reported in units of (mg/L)^-1 (Yin 2021 Table S2 parameter row 6), so CAV must be supplied in mg/L. When paired with the Yin 2020 popPK output which reports Cc in ng/mL, divide by 1000 before passing as CAV. Yin 2021 Results reports population Cavg percentiles as 3.8 mg/L (25th), 4.7 mg/L (50th), and 6 mg/L (75th) for pexidartinib-treated patients. Placebo periods have CAV = 0.",
      source_name        = "Cavg"
    ),
    AGE = list(
      description        = "Age at baseline; enters model with reference 44 years (Yin 2021 Figure 2 caption reference patient).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Yin 2021 Table S2 parameterizes the age effect as `exp(theta5 * (AGE - 44))` on baseline, `exp(theta11 * (AGE - 44))` on kdrug, and `exp(theta15 * (AGE - 44))` on konset; the reported exp(theta) values are per-year multipliers relative to the 44-year reference. Longitudinal RECIST-analysis population age summary: PLX108-01 43.2 (SD 13.0) years; ENLIVEN placebo 40.9 (SD 13.7); ENLIVEN placebo-crossover 47.7 (SD 13.2); ENLIVEN pexidartinib 44.5 (SD 13.7). Table S1.",
      source_name        = "AGE"
    ),
    JOINT_SMALL = list(
      description        = "Small-joint TGCT involvement indicator, 1 = tumor located in a small joint, 0 = large joint (reference category).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (large joint; the reference joint-size category in Yin 2021 Figure 2 caption reference patient)",
      notes              = "Time-fixed per subject. Yin 2021 does not enumerate which anatomic joints are grouped as 'small' vs 'large'; the paper treats it as an investigator-assigned binary stratifier of TGCT anatomic location. Small-joint TGCT is reported to have both a lower baseline tumor size and a larger drug-effect rate constant (kdrug), consistent with clinical observation that small-joint disease has greater treatment response. Longitudinal RECIST-analysis population: PLX108-01 large 74% / small 26%; ENLIVEN placebo 81% / 19%; ENLIVEN placebo-crossover 77% / 23%; ENLIVEN pexidartinib 64% / 36%. Table S1.",
      source_name        = "JOINT_SIZE (paper's small vs large indicator; encoded as small = 1)"
    ),
    TUMEXT_UPPER = list(
      description        = "Upper-extremity TGCT tumor location indicator, 1 = primary tumor in an upper extremity, 0 = lower extremity (reference category).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (lower extremity; the reference extremity category in Yin 2021 Figure 2 caption reference patient)",
      notes              = "Time-fixed per subject. Enters the model as multiplicative effects `exp(theta3)^TUMEXT_UPPER` on baseline, `exp(theta9)^TUMEXT_UPPER` on kdrug, and `exp(theta13)^TUMEXT_UPPER` on konset (Yin 2021 Table S2 rows 3, 9, 13). Longitudinal RECIST-analysis population is heavily lower-extremity: PLX108-01 lower 90% / upper 10%; ENLIVEN placebo 96% / 4%; ENLIVEN placebo-crossover 87% / 13%; ENLIVEN pexidartinib 91% / 9% (Table S1). The paper reports the upper-extremity 95% CIs as wide and non-significant for kdrug (0.101-2.75) and konset (0.741-82.4), consistent with the small upper-extremity cell counts.",
      source_name        = "Tumor extremity (paper's upper vs lower indicator; encoded as upper = 1)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 141L,
    n_studies       = 2L,
    n_observations  = 781L,
    age_median      = "~44 years (pooled longitudinal-RECIST analysis population; Yin 2021 Figure 2 caption reference age)",
    weight_median   = "~80 kg (pooled longitudinal-RECIST analysis population; Yin 2021 Figure 2 caption reference weight)",
    sex_female_pct  = NA_real_,
    race_ethnicity  = c(White = round(122 / 141 * 100, 1), `non-White` = round(19 / 141 * 100, 1)),
    disease_state   = "Adult patients with tenosynovial giant cell tumor (TGCT), pooled from Phase 1 dose-escalation study PLX108-01 (N = 31; 200-1200 mg/day pexidartinib) and Phase 3 randomised placebo-controlled study PLX108-10 ENLIVEN (N = 110 across placebo, placebo-crossover, and pexidartinib arms; 1000 mg/day x 2 weeks then 800 mg/day pexidartinib for the active-treatment arms, placebo comparator, and crossover from placebo to pexidartinib at Week 25). 18 patients from the pooled screening cohort were excluded (10 with no tumor observations, 7 excluded from PK modeling, 1 with only a single tumor observation), leaving 141 patients with 781 longitudinal RECIST and TVS observations.",
    dose_range      = "200-1200 mg/day pexidartinib (PLX108-01 dose escalation and extension). ENLIVEN active arms: 1000 mg/day for 2 weeks then 800 mg/day for 22 weeks (Part 1) or 800 mg/day continuous (Part 2 / placebo-crossover); placebo comparator during Part 1. The exposure metric fed into the PD model is running average pexidartinib concentration Cavg (mg/L), which is 0 for placebo periods.",
    regions         = "International (multi-regional PLX108-01 and ENLIVEN).",
    sampling_window = "Longitudinal RECIST measurements at Weeks 0, 6, 13, 25 (and later per protocol) by blinded independent central review; 781 observations across 141 patients.",
    tumor_type      = "Tenosynovial giant cell tumor (TGCT); pooled from PLX108-01 (Phase 1) and ENLIVEN (Phase 3), all with locally advanced disease.",
    assay           = "Blinded independent central review of RECIST 1.1 (sum of longest diameters of target lesions in cm; tumor volume score TVS also measured; the primary endpoint modelled here is RECIST).",
    notes           = "The published model in Yin 2021 Table S2 is the longitudinal RECIST final model; the paper also fits a parallel longitudinal TVS model with the same structural form but different parameter estimates (not tabulated in the supplement). The pexidartinib PK is not fit in this paper -- Cavg is derived from an upstream popPK model (Yin 2020, `modellib('Yin_2020_pexidartinib')`). Reference patient (Yin 2021 Figure 2 caption): 44-year-old patient presenting with a TGCT in the lower extremity in a large joint."
  )

  ini({
    # Structural typical values -- Yin 2021 Table S2 "Point estimate" column.
    # The paper reports "exp(theta_i)" values in the table for parameters
    # modelled on the log scale (baseline, kdrug, konset, and all covariate
    # multipliers), so the linear-scale values are taken directly from Table S2
    # and re-logged here for the standard `l` (log) prefix convention.
    lbase   <- log(9.36)     ; label("Typical baseline tumor size Y0 at reference (cm)")                              # Yin 2021 Table S2 theta1 = 9.36 (95% CI 8.04, 10.7), RSE 7.18%
    lkdrug  <- log(0.196)    ; label("Typical drug-effect rate constant kdrug at reference ((mg/L)^-1)")             # Yin 2021 Table S2 exp(theta6) = 0.196 (95% CI 0.101, 0.377), RSE 20.6%
    lkonset <- log(0.000225) ; label("Typical drug-effect onset rate constant konset at reference (1/h)")           # Yin 2021 Table S2 exp(theta8) = 0.000225 (95% CI 0.000119, 0.000425), RSE 3.86%

    # Fixed structural effects. Yin 2021 Results: Emax fixed at 0.999 (Table S2 theta7);
    # placebo-cohort tumor growth rate estimated at 0.227 cm/yr with a 95% CI that included
    # zero (-0.13, 0.583), so it was fixed to 0 for subsequent drug-effect model development.
    emax_pd   <- fixed(0.999) ; label("Maximal achievable drug effect Emax (unitless)")                        # Yin 2021 Table S2 theta7 FIXED = 0.999
    growth_pd <- fixed(0)     ; label("Natural tumor growth rate (cm/h)")                                # Yin 2021 Table S2 theta2 FIXED = 0.00 cm/yr; encoded as 0/h

    # Multiplicative covariate effects on baseline Y0, kdrug, and konset. The
    # values in Yin 2021 Table S2 are the already-exponentiated multipliers
    # `exp(theta_i)`; the parameterization is `Y0 = base_TV * exp(e_ext_base *
    # TUMEXT_UPPER) * exp(e_size_base * JOINT_SMALL) * exp(e_age_base * (AGE - 44))`,
    # so the stored values are `log(exp(theta_i))` = theta_i on the natural-log scale.
    e_ext_base   <- log(0.692) ; label("Effect of TUMEXT_UPPER on baseline Y0 (multiplicative on log scale)")         # Yin 2021 Table S2 exp(theta3)  = 0.692 (95% CI 0.418, 1.15),  RSE 70.0%
    e_size_base  <- log(0.581) ; label("Effect of JOINT_SMALL on baseline Y0 (multiplicative on log scale)")          # Yin 2021 Table S2 exp(theta4)  = 0.581 (95% CI 0.441, 0.767), RSE 26.1%
    e_age_base   <- log(1.01)  ; label("Effect of AGE per year (rel 44) on baseline Y0 (multiplicative on log scale)") # Yin 2021 Table S2 exp(theta5)  = 1.01  (95% CI 1.00, 1.02),   RSE 40.3%

    e_ext_kdrug  <- log(0.527) ; label("Effect of TUMEXT_UPPER on kdrug (multiplicative on log scale)")               # Yin 2021 Table S2 exp(theta9)  = 0.527 (95% CI 0.101, 2.75),  RSE 132%
    e_size_kdrug <- log(2.14)  ; label("Effect of JOINT_SMALL on kdrug (multiplicative on log scale)")                # Yin 2021 Table S2 exp(theta10) = 2.14  (95% CI 0.886, 5.19), RSE 59.1%
    e_age_kdrug  <- log(0.918) ; label("Effect of AGE per year (rel 44) on kdrug (multiplicative on log scale)")      # Yin 2021 Table S2 exp(theta11) = 0.918 (95% CI 0.880, 0.958), RSE 25.2%

    e_ext_konset  <- log(7.82) ; label("Effect of TUMEXT_UPPER on konset (multiplicative on log scale)")              # Yin 2021 Table S2 exp(theta13) = 7.82  (95% CI 0.741, 82.4), RSE 58.4%
    e_size_konset <- log(1.82) ; label("Effect of JOINT_SMALL on konset (multiplicative on log scale)")               # Yin 2021 Table S2 exp(theta14) = 1.82  (95% CI 0.772, 4.28), RSE 73.1%
    e_age_konset  <- log(1.01) ; label("Effect of AGE per year (rel 44) on konset (multiplicative on log scale)")     # Yin 2021 Table S2 exp(theta15) = 1.01  (95% CI 0.980, 1.05), RSE 122%

    # Inter-individual variability -- 3x3 log-normal block on (baseline, kdrug,
    # konset). Yin 2021 Table S2 reports Omega variances on the diagonal and
    # off-diagonal covariances (with correlations shown in the [] bracket
    # annotation). Values verified via correlation calc:
    #   rho(base,kdrug)  = -0.161 / sqrt(0.410 * 2.22) = -0.169 ~ -0.168 (paper)
    #   rho(base,konset) = -0.157 / sqrt(0.410 * 2.18) = -0.166 ~ -0.166 (paper)
    #   rho(kdrug,konset)=  0.925 / sqrt(2.22  * 2.18) =  0.421 (paper -0.166 appears to be a supplement typo;
    #                                                              the 0.925 numeric covariance is retained here)
    #   %CV verified: sqrt(exp(0.410) - 1) = 71.2% (paper); sqrt(exp(2.22) - 1) = 287% (paper);
    #                 sqrt(exp(2.18) - 1) = 280% (paper).
    etalbase + etalkdrug + etalkonset ~ c(0.410,
                                          -0.161,  2.22,
                                          -0.157,  0.925, 2.18)

    # Residual error -- power form: SD(tumor_size_pred) = propSd * tumor_size_pred^powExp.
    # Yin 2021 Table S2: proportional-residual SD = 0.365 (95% CI 0.595, 0.613 [note:
    # the reported CI in the supplement table brackets a value inconsistent with the
    # point estimate; the point estimate 0.365 is retained]) and power parameter
    # residuals = 0.550 (95% CI 0.518, 0.582), RSE 2.97%.
    propSd <- 0.365 ; label("Power-error SD coefficient (cm^(1 - powExp))")                                           # Yin 2021 Table S2 residual SD (proportional) = 0.365, RSE 4.20%
    powExp <- 0.550 ; label("Power-error exponent on tumor_size prediction (unitless)")                               # Yin 2021 Table S2 power parameter residuals = 0.550, RSE 2.97%
  })

  model({
    # Age enters as a per-year deviation from the 44-year reference (Yin 2021
    # Figure 2 caption reference patient). Joint size and tumor extremity are
    # binary indicators; the reference-patient values are JOINT_SMALL = 0
    # (large joint) and TUMEXT_UPPER = 0 (lower extremity).
    age_dev <- AGE - 44

    # Individual baseline tumor size Y0, drug-effect rate constant kdrug, and
    # drug-effect onset rate constant konset -- log-normal IIV with a 3x3
    # correlated block; multiplicative covariate effects (Yin 2021 Table S2).
    y0     <- exp(lbase   + etalbase   + e_ext_base   * TUMEXT_UPPER + e_size_base   * JOINT_SMALL + e_age_base   * age_dev)
    kdrug  <- exp(lkdrug  + etalkdrug  + e_ext_kdrug  * TUMEXT_UPPER + e_size_kdrug  * JOINT_SMALL + e_age_kdrug  * age_dev)
    konset <- exp(lkonset + etalkonset + e_ext_konset * TUMEXT_UPPER + e_size_konset * JOINT_SMALL + e_age_konset * age_dev)

    # Algebraic RECIST tumor size prediction (Yin 2021 Eq. 1 with the +growth
    # linear term collapsed because growth is FIXED to 0). CAV is supplied by
    # the user as a time-varying covariate representing the running-average
    # pexidartinib plasma concentration (mg/L) up to the current record time.
    # `t` is the rxode2 time variable; the paper's TAFD (time after first
    # dose) is `t` under the standard convention that the event table's
    # time origin is the first pexidartinib dose. Placebo periods carry
    # CAV = 0 and reduce tumor_size to y0.
    tumor_size <- y0 * (1 - emax_pd * (1 - exp(-kdrug * CAV)) * (1 - exp(-konset * t))) + growth_pd * t

    # Power-of-prediction residual error: SD(tumor_size) = propSd * tumor_size^powExp.
    tumor_size ~ pow(propSd, powExp)
  })
}
