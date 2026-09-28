Vong_2021_tofacitinib <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and an",
    "absorption lag time for oral immediate-release tofacitinib in 1096 adults",
    "with moderately to severely active ulcerative colitis, pooled from the",
    "phase 2 dose-ranging induction study A3921063, the phase 3 OCTAVE",
    "Induction 1 and 2 studies and the phase 3 OCTAVE Sustain maintenance",
    "study (Vong 2021). The model is parameterized in apparent oral clearance",
    "(CL/F) and apparent volume of distribution (V/F). CL/F varies with",
    "baseline creatinine clearance (power 0.354 on CRCL_BASE/108.7) and",
    "multiplicative factors for female sex (0.868) and Asian race (0.932);",
    "V/F varies with baseline body weight (power 0.585 on WT/72), age (power",
    "-0.116 on AGE/40) and female sex (0.845). Inter-individual variability is",
    "an exponential eta on CL/F only; V/F has no eta of its own, its",
    "individual deviation being the paper's 'scaling parameter' (0.392) times",
    "the CL/F eta. Ka carries inter-occasion variability (six occasions, no",
    "inter-individual variability). Residual error is proportional (additive",
    "on the log scale) with a magnitude that switches on time after dose at",
    "8 hours (41.6% at or before 8 h, 68.9% after), and the residual",
    "magnitude itself carries a 58.0% inter-individual variability (etaruv).",
    sep = " "
  )
  reference <- paste(
    "Vong C, Martin SW, Deng C, Xie R, Ito K, Su C, Sandborn WJ, Mukherjee A.",
    "Population Pharmacokinetics of Tofacitinib in Patients With Moderate to",
    "Severe Ulcerative Colitis. Clin Pharmacol Drug Dev. 2021; 10(3): 229-240.",
    "doi:10.1002/cpdd.899. PMCID: PMC7986169.",
    sep = " "
  )
  vignette <- "Vong_2021_tofacitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL_BASE = list(
      description = "Baseline creatinine clearance",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline, time-fixed, and not body-surface-area normalized (the paper",
        "reports BCCL in mL/min; the estimating equation is not stated, the",
        "sibling tofacitinib analyses from the same sponsor use Cockcroft-Gault).",
        "Enters CL/F as a power function normalized to 108.7 mL/min, the",
        "cohort median (Results, Table 2). Cohort mean 112.3 mL/min (SD 30.3),",
        "range 40.8-255.2; creatinine clearance below 40 mL/min was an",
        "exclusion criterion in every study, so the model should not be",
        "extrapolated to moderate or severe renal impairment."
      ),
      source_name = "BCCL"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = paste(
        "Table 3 reports the 'Female' effects as MULTIPLICATIVE fractions of the",
        "male reference ('estimates were computed as a fraction of the reference",
        "category', Table 3 footnote d): 0.868 on CL/F and 0.845 on V/F, i.e.",
        "13.2% lower CL/F and 15.5% lower V/F in females (Results). Cohort 41.5%",
        "female (Table 2)."
      ),
      source_name = "Female"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "non-Asian (RACE_ASIAN = 0)",
      notes = paste(
        "The paper evaluated each race category against the rest of the",
        "population, so the reference is all non-Asian patients (White, Black",
        "and Other pooled). Multiplicative factor 0.932 on CL/F (Table 3), i.e.",
        "6.8% lower CL/F in Asian patients. Cohort 10.9% Asian (Table 2)."
      ),
      source_name = "Asian"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline (time-fixed), reported as BWT. Enters V/F as a power function",
        "normalized to 72.0 kg, the cohort median. Cohort mean 73.6 kg",
        "(SD 16.6), range 37-154.5 (Table 2); 53 and 95 kg are the 10th and 90th",
        "percentiles used in Figure 1."
      ),
      source_name = "BWT"
    ),
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline (time-fixed). Enters V/F as a power function normalized to",
        "40.0 years, the cohort median. Cohort mean 41.3 years (SD 13.8), range",
        "18-80 (Table 2)."
      ),
      source_name = "Age"
    ),
    OCC = list(
      description = "Occasion index for inter-occasion variability on Ka",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Integer occasion 1-6; any other value (e.g. 0) switches the",
        "inter-occasion variability on Ka off. The paper estimated IOV on Ka",
        "over six occasions (eta-shrinkage 46.9-87.4% across 'the 6 occasions",
        "estimated', Results) and one shared magnitude (Table 3, 191.8%), so",
        "occasions 2-6 are fixed to the occasion-1 variance (the NONMEM",
        "$OMEGA BLOCK(1) SAME idiom). The paper does not list which study visit",
        "maps to which occasion; the protocol PK visits were baseline and weeks",
        "2, 4 and 8 of A3921063, baseline/week 2 and week 8 of OCTAVE Induction",
        "1 and 2, and baseline and weeks 8, 24 and 52 of OCTAVE Sustain",
        "(Table 1). Assign one occasion per PK visit in the order they occur.",
        "Ka has no inter-individual variability of its own, so a single-occasion",
        "simulation (OCC = 1 on every record) makes the occasion-1 eta act as a",
        "per-subject absorption deviation."
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "tofacitinib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "tofacitinib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1096,
    n_studies = 4,
    n_observations = 7231,
    age_median = "40 years",
    age_mean = "41.3 years (SD 13.8)",
    age_range = "18-80 years",
    weight_median = "72 kg",
    weight_mean = "73.6 kg (SD 16.6)",
    weight_range = "37-154.5 kg",
    sex_female_pct = 41.5,
    race_ethnicity = c(White = 81.3, Asian = 10.9, Other = 4.2, Black = 1.0),
    ethnicity = c(Hispanic = 3.4, NonHispanic = 80.9),
    disease_state = "moderately to severely active ulcerative colitis (baseline total Mayo score mean 8.9, range 3-12); 50.9% TNF-inhibitor naive, 46.6% prior TNF-inhibitor failure",
    dose_range = "0.5, 3, 5, 10 or 15 mg orally twice daily (immediate-release)",
    renal_function = "baseline creatinine clearance median 108.7 mL/min, range 40.8-255.2; creatinine clearance < 40 mL/min was an exclusion criterion",
    co_medication = "background 5-aminosalicylates (77.7%) and oral corticosteroids (42.8%) allowed; azathioprine, 6-mercaptopurine, methotrexate and TNF inhibitors prohibited",
    regions = "multinational phase 2 (A3921063, NCT00787202) and phase 3 (OCTAVE Induction 1 NCT01465763, OCTAVE Induction 2 NCT01458951, OCTAVE Sustain NCT01458574) studies",
    notes = paste(
      "Demographics from Table 2 of Vong 2021 (race does not sum to 100%; the",
      "remaining 2.6% is unreported). Patient numbers per study (Table 1):",
      "A3921063 143 (0.5 mg n = 31, 3 mg n = 33, 10 mg n = 31, 15 mg n = 48),",
      "OCTAVE Induction 1 483 (10 mg n = 468, 15 mg n = 15), OCTAVE Induction 2",
      "427 (10 mg n = 422, 15 mg n = 5), OCTAVE Sustain 380 (5 mg n = 191,",
      "10 mg n = 189); patients continuing from induction into maintenance are",
      "counted once in the 1096. Sampling was sparse in phase 3 (predose plus",
      "one post-dose sample, 0.5 h or 2 h after an in-clinic dose) and serial in",
      "phase 2 (predose, 0.25, 0.5, 1 and 2-3 h at baseline and week 8). About",
      "4% of observations were below the 0.100 ng/mL LLOQ and were discarded,",
      "and 519 week-2/week-4 predose samples from A3921063 were excluded for",
      "data-collection errors (Results)."
    )
  )

  ini({
    # ---- Structural parameters -------------------------------------------
    # Typical values are for the paper's reference patient: non-Asian male,
    # body weight 72.0 kg, age 40.0 years, baseline creatinine clearance
    # 108.7 mL/min (Results, Final Model Results; Figure 1 legend).
    lka <- log(9.85); label("Absorption rate constant (1/h)")                     # Results and Abstract, Ka = 9.85 h^-1 (bootstrap 95% CI 7.9-10.8); Table 3 prints the same estimate rounded to 9.9 (RSE 7.1%)
    lcl <- log(26.3); label("Apparent oral clearance CL/F (L/h)")                 # Table 3, 'CL/F, L/h' = 26.3 (RSE 1.2%), bootstrap 95% CI 25.5-27.2
    lvc <- log(115.8); label("Apparent volume of distribution V/F (L)")           # Table 3, 'V/F, L' = 115.8 (RSE 1.1%), bootstrap 95% CI 111.5-120.6
    ltlag <- log(0.236); label("Absorption lag time (h)")                         # Table 3, 'Lag time, h' = 0.236 (RSE 0.52%), bootstrap 95% CI 0.218-0.238

    # ---- Covariate effects on CL/F ---------------------------------------
    # Continuous covariates are power functions of (cov / median); categorical
    # covariates are MULTIPLICATIVE fractions of the reference category, so
    # their null value is 1 (Table 3 footnote d).
    e_crcl_base_cl <- 0.354; label("Power exponent on (CRCL_BASE/108.7) for CL/F (unitless)")   # Table 3, covariate CL/F ~ BCCL = 0.354 (RSE 7.57%), 95% CI 0.279-0.426
    e_sexf_cl <- 0.868; label("Multiplicative factor on CL/F for female vs male (unitless)")     # Table 3, covariate CL/F ~ Female = 0.868 (RSE 1.86%), 95% CI 0.829-0.909
    e_race_asian_cl <- 0.932; label("Multiplicative factor on CL/F for Asian vs non-Asian (unitless)")  # Table 3, covariate CL/F ~ Asian = 0.932 (RSE 2.16%), 95% CI 0.888-0.972

    # ---- Covariate effects on V/F ----------------------------------------
    e_wt_vc <- 0.585; label("Power exponent on (WT/72) for V/F (unitless)")                     # Table 3, covariate V/F ~ Body weight = 0.585 (RSE 5.03%), 95% CI 0.501-0.674
    e_age_vc <- -0.116; label("Power exponent on (AGE/40) for V/F (unitless)")                  # Table 3, covariate V/F ~ Age = -0.116 (RSE 15.05%), 95% CI -0.173 to -0.065
    e_sexf_vc <- 0.845; label("Multiplicative factor on V/F for female vs male (unitless)")      # Table 3, covariate V/F ~ Female = 0.845 (RSE 1.54%), 95% CI 0.802-0.887

    # ---- IIV / IOV -------------------------------------------------------
    # SCALE CONVENTION (assumption; see the vignette). Table 3 prints IIV and
    # IOV as percentages without stating whether they are 100*sqrt(omega^2)
    # or the exact log-normal CV 100*sqrt(exp(omega^2)-1). The variances below
    # use (percent/100)^2, the convention adopted for the sibling Xie 2019
    # tofacitinib analysis from the same sponsor. It matters little for CL/F
    # (22.2%) and a great deal for the Ka IOV (191.8%): under the exact-CV
    # reading that variance would be log(1 + 1.918^2) = 1.543 instead of 3.679.
    etalcl ~ 0.049284                                                            # Table 3, 'CL/F, L/h' IIV = 22.2% (RSE 7.7%): 0.222^2 = 0.049284
    # V/F has NO eta of its own (Table 3 IIV 'NA'): Equation 1b,
    # V_i = theta_TV_V * exp(eta_CL,i * theta_scale), implemented in model().
    vc_eta_scale <- 0.392; label("Scaling factor relating the V/F deviation to etalcl (unitless)")  # Table 3, 'Scaling parameter' = 0.392 (RSE 7.19%), bootstrap 95% CI 0.305-0.481
    # Ka inter-occasion variability, six occasions with one shared magnitude
    # (Results: 'eta-shrinkage for IOV, for the 6 occasions estimated'). Ka
    # carries no IIV ('only IOV was included as a random effect on Ka').
    etaiov_lka_1 ~ 3.678724                                                      # Table 3, 'Ka' IOV = 191.8% (RSE 7.0%): 1.918^2 = 3.678724
    etaiov_lka_2 ~ fixed(3.678724)                                               # Table 3, 'Ka' IOV = 191.8%; occasion 2 shares the occasion-1 variance
    etaiov_lka_3 ~ fixed(3.678724)                                               # Table 3, 'Ka' IOV = 191.8%; occasion 3 shares the occasion-1 variance
    etaiov_lka_4 ~ fixed(3.678724)                                               # Table 3, 'Ka' IOV = 191.8%; occasion 4 shares the occasion-1 variance
    etaiov_lka_5 ~ fixed(3.678724)                                               # Table 3, 'Ka' IOV = 191.8%; occasion 5 shares the occasion-1 variance
    etaiov_lka_6 ~ fixed(3.678724)                                               # Table 3, 'Ka' IOV = 191.8%; occasion 6 shares the occasion-1 variance
    etaruv ~ 0.3364                                                              # Table 3, IIV on both proportional-error rows = 58.0% (RSE 6.3%): 0.580^2 = 0.3364

    # ---- Residual error --------------------------------------------------
    # Log-transformed data with additive error on the log scale, i.e.
    # proportional on the linear scale (Methods), with two magnitudes split on
    # time after dose at 8 h (Results).
    propSd_early <- 0.416; label("Proportional residual error SD, time after dose <= 8 h (fraction)")  # Table 3, 'Proportional error, TAD <= 8 h, %' = 41.6 (RSE 2.6%)
    propSd_late <- 0.689; label("Proportional residual error SD, time after dose > 8 h (fraction)")    # Table 3, 'Proportional error, TAD > 8 h, %' = 68.9 (RSE 2.6%)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Covariate models (Results; Table 3 footnote d).
    #    Verified against every covariate effect the paper states in prose:
    #      (60/108.7)^0.354 = 0.810 -> 19% lower CL/F at BCCL 60 mL/min [paper: 19%]
    #      0.868                    -> 13.2% lower CL/F in females     [paper: 13.2%]
    #      0.932                    ->  6.8% lower CL/F in Asians      [paper: 6.8%]
    #      (53/72)^0.585   = 0.836  -> 16.4% lower V/F at 53 kg        [paper: 16.4%]
    #      (95/72)^0.585   = 1.176  -> 17.6% higher V/F at 95 kg       [paper: 17.6%]
    #      0.845                    -> 15.5% lower V/F in females      [paper: 15.5%]
    #      (80/40)^-0.116  = 0.923  ->  7.7% lower V/F at age 80       [paper: 7.7%]
    #      log(2) * 115.8 / 26.3    =  3.05 h elimination half-life    [paper: 3.05 h]
    # ------------------------------------------------------------------
    cov_cl <- (CRCL_BASE / 108.7)^e_crcl_base_cl *
      (1 + (e_sexf_cl - 1) * SEXF) *
      (1 + (e_race_asian_cl - 1) * RACE_ASIAN)
    cov_vc <- (WT / 72)^e_wt_vc *
      (AGE / 40)^e_age_vc *
      (1 + (e_sexf_vc - 1) * SEXF)

    # ------------------------------------------------------------------
    # 2. Inter-occasion variability on Ka, multiplexed by OCC (1-6).
    # ------------------------------------------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    iov_ka <- oc1 * etaiov_lka_1 + oc2 * etaiov_lka_2 + oc3 * etaiov_lka_3 +
      oc4 * etaiov_lka_4 + oc5 * etaiov_lka_5 + oc6 * etaiov_lka_6

    # ------------------------------------------------------------------
    # 3. Individual parameters (Equations 1a and 1b). V/F has no eta of its
    #    own: its log-scale deviation is vc_eta_scale * etalcl.
    # ------------------------------------------------------------------
    ka <- exp(lka + iov_ka)
    cl <- exp(lcl + etalcl) * cov_cl
    vc <- exp(lvc + vc_eta_scale * etalcl) * cov_vc
    tlag <- exp(ltlag)

    kel <- cl / vc

    # ------------------------------------------------------------------
    # 4. One-compartment disposition with lagged first-order absorption.
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # 5. Observation and error. Dose in mg and vc in L give mg/L; x1000 for
    #    ng/mL (assay LLOQ 0.100 ng/mL). The proportional magnitude switches
    #    on time after dose at 8 h and carries its own inter-individual
    #    variability, so it is assembled as a model variable. tad() is
    #    evaluated once on its own line.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc

    early_flag <- tad() <= 8
    propSdTad <- propSd_early * early_flag + propSd_late * (1 - early_flag)
    propSdCc <- propSdTad * exp(etaruv)
    Cc ~ prop(propSdCc)
  })
}
