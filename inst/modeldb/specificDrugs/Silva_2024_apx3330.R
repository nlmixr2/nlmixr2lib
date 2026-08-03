Silva_2024_apx3330 <- function() {
  description <- paste(
    "Population pharmacokinetic model for APX3330 (E3330), a selective",
    "APE1/Ref-1 redox inhibitor, measured as total quinone in plasma after",
    "oral dosing. Two-compartment disposition with first-order absorption and",
    "an absorption lag time, fitted to 1460 plasma concentrations pooled from",
    "four Eisai studies in 49 healthy Japanese male volunteers (single",
    "dose-escalation 10-600 mg, multiple-dose 120 mg once or twice daily, and",
    "a fasted-vs-fed crossover) and one Apexian phase I study in 19 patients",
    "with advanced solid tumors dosed 120-360 mg twice daily. All disposition",
    "parameters are apparent (oral) values because absolute bioavailability was",
    "never measured. Body weight is a power-model covariate on CL/F (exponent",
    "0.659) and V1/F (exponent 0.839) referenced to 70 kg; the oncology cohort",
    "carries a log-additive shift on CL/F, and dosing in the fed state carries a",
    "log-additive shift on the absorption lag time. Inter-individual variability",
    "is estimated on CL/F, V1/F and ka with a CL/F-V1/F correlation of 0.349,",
    "and residual error is proportional (16%). The paper's companion GastroPlus",
    "ACAT semi-mechanistic absorption model is a commercial-platform model and",
    "is not reproduced here; see the validation vignette."
  )
  reference <- paste(
    "Silva LL, Stratford RE, Messmann R, Kelley MR, Quinney SK.",
    "Bridging population pharmacokinetic and semimechanistic absorption",
    "modeling of APX3330.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(1):106-117.",
    "doi:10.1002/psp4.13061."
  )
  vignette <- "Silva_2024_apx3330"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-model covariate on CL/F and V1/F, normalized to a 70 kg",
        "reference. Silva 2024 Methods 'Population pharmacokinetic analysis':",
        "'Covariate effects were evaluated for body weight (WT; normalized to",
        "70 kg)'. The 70 kg reference is stated explicitly and is NOT the",
        "cohort median: healthy Japanese volunteers weighed 52-77 kg (mean",
        "62.4 kg per the semi-physiologic modeling section) and the oncology",
        "cohort averaged 88.3 kg (Results 'Population pharmacokinetic",
        "analysis'). Applied via Eq. 1, the continuous-covariate power model",
        "p_j = p_pop x (cov_j / cov_ref)^beta. Body weight was recorded on",
        "three to eight occasions in the oncology cohort with no trend over",
        "the study, so it is treated as time-fixed here.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    FED = list(
      description        = "Fed-vs-fasted dose-record indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted dose)",
      notes              = paste(
        "Per-dose-record indicator: 1 = APX3330 dose administered with food,",
        "0 = fasted. Identified in the 120 mg fasted-vs-fed single-dose",
        "crossover arm of the Japanese healthy-volunteer program (Silva 2024",
        "Methods 'Plasma concentration data'). Retained only on the absorption",
        "lag time (Results: delta -2LL = -34.06, 'explained more than 15% of",
        "IOV of this parameter'). Food also affected ka during forward",
        "inclusion but failed backward elimination and is therefore absent",
        "from the final model. The source does not report the meal",
        "composition, so the general FED indicator applies rather than",
        "FED_HIGHFAT or FED_LOWFAT; the Discussion only speculates that 'a",
        "high-fat meal ... is commonly used in phase I trials'.",
        sep = " "
      ),
      source_name        = "Food"
    ),
    DIS_CANCER = list(
      description        = "Advanced-solid-tumor cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy Japanese male volunteer)",
      notes              = paste(
        "Time-fixed per subject: 1 = patient with an advanced solid tumor in",
        "the Apexian phase I study, 0 = healthy Japanese male volunteer in the",
        "pooled Eisai studies. Silva 2024 calls this covariate 'subject",
        "source' and states explicitly that it was preferred over serum",
        "albumin 'because it represents a combination of all the differences",
        "between the two cohorts' (Results), including disease state, age and",
        "ethnicity as well as the albumin difference (3.9 +/- 0.28 g/dL in",
        "patients with cancer vs 4.8 +/- 0.16 g/dL in healthy volunteers).",
        "Adding subject source on CL/F gave delta -2LL = -39.52 versus -25.1",
        "for albumin. Serum albumin is therefore documented in",
        "covariatesDataExcluded rather than carried as a model covariate.",
        sep = " "
      ),
      source_name        = "Subject source"
    )
  )

  covariatesDataExcluded <- list(
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Screened as a candidate covariate on CL/F and significant on its own",
        "(delta -2LL = -25.1), but NOT retained in the final model: the",
        "categorical DIS_CANCER 'subject source' indicator gave a larger drop",
        "(delta -2LL = -39.52) and subsumes the albumin difference along with",
        "the other between-cohort differences. Silva 2024 Discussion: 'we",
        "opted to maintain it in the model instead of albumin concentration'.",
        "Reported cohort means are 3.9 +/- 0.28 g/dL (39 g/L) in patients with",
        "cancer and 4.8 +/- 0.16 g/dL (48 g/L) in healthy volunteers; both are",
        "within the normal range. No usable point estimate for an albumin",
        "coefficient is reported.",
        sep = " "
      ),
      source_name = "Albumin"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 68,
    n_studies      = 5,
    age_range      = "healthy Japanese volunteers: mean 25.4 years; oncology cohort age not reported in the main text (Table S3)",
    weight_range   = "healthy Japanese volunteers: 52-77 kg (mean 62.4 kg); oncology cohort: mean 88.3 kg",
    sex_female_pct = 8.8,
    race_ethnicity = c(Japanese = 72.1, White = 23.5, Hispanic = 2.9, Other = 1.5),
    disease_state  = "healthy Japanese male volunteers pooled with patients with advanced solid tumors",
    dose_range     = "10-600 mg oral single dose and 120 mg once or twice daily (healthy volunteers); 120-360 mg orally twice daily for 22 days (patients with cancer)",
    regions        = "Japan (Eisai healthy-volunteer studies); United States (Apexian oncology study)",
    notes          = paste(
      "Silva 2024 Results 'Population pharmacokinetic analysis'. The combined",
      "model was built on 1460 total-quinone plasma concentrations, of which",
      "211 (14.45%) came from the oncology study. The healthy-volunteer",
      "component pools four Eisai studies (single dose-escalation 10/30/60/",
      "120/180/240 mg; multiple dose 120 mg QD or BID; a 120 mg fasted-vs-fed",
      "single-dose crossover; and a high single-dose study of 300/420/600 mg)",
      "whose individual-level data were digitized from clinical trial reports",
      "with Engauge Digitizer; all 49 participants were Japanese men within",
      "+/- 20% of ideal body weight. The oncology component is a single",
      "Apexian phase I study in 19 patients: 13 (68.42%) men, 84.21% White",
      "and 10.53% Hispanic. The race_ethnicity percentages above are derived",
      "from those two cohort descriptions (49 Japanese volunteers plus 16",
      "White, 2 Hispanic and 1 other of the 19 oncology patients) because the",
      "source reports ethnicity per cohort, not pooled. Assay note: neither",
      "the Eisai HPLC/UV nor the Apexian HPLC-MS/MS method discriminated the",
      "quinone from the hydroquinone form; samples were oxidized before",
      "analysis, so every concentration the model was fitted to is TOTAL",
      "quinone. Neither dose level nor study period was a significant",
      "covariate, indicating linear PK across 10-600 mg.",
      sep = " "
    )
  )

  ini({
    # =========================================================================
    # Silva 2024 Table 2, 'Combined data' column (the paper's final model,
    # fitted to the pooled healthy-volunteer + oncology dataset). The
    # 'Healthy volunteers' column of the same table is the preceding stage of
    # the authors' staged model-development strategy ("Using the final
    # estimates from the model for healthy volunteers as initial estimates, a
    # final model was developed to describe data from both studies", Methods)
    # and is therefore not carried as a separate model file; setting
    # DIS_CANCER = 0 here reproduces the healthy-Japanese-volunteer reference
    # population (CL/F 193 vs 196 mL/h, t_lag 0.401 vs 0.413 h in the earlier
    # stage, i.e. within 2%).
    #
    # The covariate model is Eq. 1 (continuous, power) and Eq. 2 (categorical,
    # log-additive):
    #   Eq. 1   p_j = p_pop x (cov_j / cov_ref)^beta_cov x exp(eta_j + eta_OCC)
    #   Eq. 2   p_j = p_pop x exp(beta_cov x cov_j)      x exp(eta_j + eta_OCC)
    #
    # Applied to this model:
    #   t_lag = 0.401 x exp(0.815 x FED)                       h
    #   ka    = 0.894                                          1/h
    #   CL/F  = 193  x (WT/70)^0.659 x exp(0.409 x DIS_CANCER)  mL/h
    #   V1/F  = 4050 x (WT/70)^0.839                            mL
    #   Q/F   = 276                                             mL/h
    #   V2/F  = 4780                                            mL
    #
    # UNIT CONVERSION: Silva 2024 reports clearances in mL/h and volumes in
    # mL. This file works in L and L/h so that a dose in mg divided by a
    # volume in L yields mg/L == ug/mL, matching the concentration units of
    # the paper's own Table 3 (observed C_max 28.0 ug/mL at 120 mg fasted;
    # 120 mg / 4.05 L = 29.6 ug/mL). Every conversion is a factor of 1000 and
    # is noted on the parameter line.
    # =========================================================================

    ltlag <- log(0.401); label("Absorption lag time t_lag, fasted (h)")                       # Table 2 'Combined data', row 't lag (h)' = 0.401 (RSE 6.92%)
    lka   <- log(0.894); label("First-order absorption rate constant ka (1/h)")               # Table 2 'Combined data', row 'k a (1/h)' = 0.894 (RSE 9.47%)
    lcl   <- log(0.193); label("Apparent oral clearance CL/F at WT = 70 kg, healthy volunteer (L/h)")   # Table 2 'Combined data', row 'CL/F (mL/h)' = 193 mL/h (RSE 2.48%) = 0.193 L/h
    lvc   <- log(4.05);  label("Apparent central volume V1/F at WT = 70 kg (L)")              # Table 2 'Combined data', row 'V 1/F (mL)' = 4050 mL (RSE 2.51%) = 4.05 L
    lq    <- log(0.276); label("Apparent intercompartmental clearance Q/F (L/h)")             # Table 2 'Combined data', row 'Q/F (mL/h)' = 276 mL/h (RSE 4.02%) = 0.276 L/h
    lvp   <- log(4.78);  label("Apparent peripheral volume V2/F (L)")                         # Table 2 'Combined data', row 'V 2/F (mL)' = 4780 mL (RSE 2.00%) = 4.78 L

    # Covariate effects. All four were estimated (no FIX flag, each reported
    # with an RSE), so none is wrapped in fixed(). The two power exponents are
    # estimated rather than fixed at the allometric 0.75 / 1.
    e_wt_cl     <- 0.659; label("Power exponent of body weight on CL/F, referenced to 70 kg (unitless)")  # Table 2 'Combined data', row 'beta WT,CL'  = 0.659 (RSE 17.6%); Eq. 1
    e_wt_vc     <- 0.839; label("Power exponent of body weight on V1/F, referenced to 70 kg (unitless)")  # Table 2 'Combined data', row 'beta WT,V1'  = 0.839 (RSE 12.8%); Eq. 1
    e_cancer_cl <- 0.409; label("Log-additive shift on CL/F for the advanced-solid-tumor cohort (unitless)") # Table 2 'Combined data', row 'beta SubjectSource,CL' = 0.409 (RSE 13.0%); Eq. 2
    e_fed_tlag  <- 0.815; label("Log-additive shift on t_lag for a dose taken with food (unitless)")      # Table 2 'Combined data', row 'beta Food,tlag' = 0.815 (RSE 13.5%); Eq. 2

    # =========================================================================
    # Random effects, Silva 2024 Table 2 'Combined data', 'Random effects' and
    # 'Correlations' blocks. Monolix reports omega / gamma as STANDARD
    # DEVIATIONS on the log scale, whereas nlmixr2's ini() takes VARIANCES, so
    # each entry below is the square of the tabulated value. The SD reading is
    # confirmed twice by the paper's own prose:
    #   - "Despite the small IIV for CL/F and V1/F estimates" -- omega_CL/F =
    #     0.142 as an SD gives CV = sqrt(exp(0.142^2) - 1) = 14.3% ("small");
    #     read as a variance it would be a 39% CV, which is not small.
    #   - "variability in absorption parameters, ka and t lag, was over 50%" --
    #     omega_ka = 0.517 -> CV 55.4%, gamma_tlag = 0.515 -> CV 55.2%,
    #     gamma_ka = 0.564 -> CV 61.0%, all just over 50% as stated.
    #
    # INTER-OCCASION VARIABILITY. Table 2 reports IOV (gamma) on t_lag (0.515)
    # and ka (0.564) in addition to IIV (omega) on ka, CL/F and V1/F. This
    # file follows the established nlmixr2lib convention for between-occasion
    # variability (Svensson_2012_nevirapine.R, Bienczak_2016_nevirapine.R,
    # Svensson_2018_bedaquiline.R, Bukkems_2021_raltegravir.R): BOV is dropped
    # when a BSV term is reported on the same parameter, and folded in as a
    # BSV-equivalent when only BOV is reported. Concretely:
    #   - ka   has both omega (0.517) and gamma (0.564): keep omega, drop gamma.
    #   - t_lag has gamma (0.515) only:  fold it in as etaltlag ~ 0.515^2.
    # The occasion-indicator encoding used in Chen_2023_nemonoxacin.R is not
    # applicable here because Silva 2024 never states how many occasions the
    # IOV spans (the dose-escalation arm gave three single-dose occasions per
    # subject and the food-effect crossover gave two), so an occasion count
    # would have to be invented. See the vignette's Assumptions and deviations
    # section.
    # =========================================================================

    # Correlated IIV on CL/F and V1/F. Covariance = corr x SD_cl x SD_vc =
    # 0.349 x 0.142 x 0.141 = 0.0069877.
    etalcl + etalvc ~ c(0.020164,
                        0.0069877, 0.019881)  # Table 2 'omega CL/F' = 0.142 -> var 0.020164; 'omega V1/F' = 0.141 -> var 0.019881; 'Correlations CL/F ~ V1/F' = 0.349 (RSE 37.0%)

    etalka   ~ 0.267289  # Table 2 'omega ka'      = 0.517 (RSE 20.1%) -> var 0.517^2 = 0.267289 (CV 55.4%)
    etaltlag ~ 0.265225  # Table 2 'gamma t lag'   = 0.515 (RSE 9.44%) -> var 0.515^2 = 0.265225 (CV 55.2%); IOV folded in as a BSV-equivalent because no IIV was estimated on t_lag

    # Table 2 also reports a t_lag ~ ka correlation of -0.476 (RSE 26.4%).
    # That correlation is between the two IOV terms (gamma_tlag and gamma_ka),
    # not between the IIV terms: t_lag has no omega row at all. Because
    # gamma_ka is dropped above and etaltlag is a folded gamma_tlag, pairing
    # -0.476 with the IIV etalka would mix variability levels, so the two etas
    # are left independent here.

    # Residual error. Table 2 'Residual variability' row 'Proportional'; the
    # Methods state "The preferred error model was proportional". Monolix's
    # proportional model is y = f x (1 + b x eps) with eps ~ N(0, 1), so the
    # tabulated b is already on the SD scale and maps directly to propSd.
    propSd <- 0.16; label("Proportional residual error (fraction)")  # Table 2 'Combined data', row 'Proportional' = 0.16 (RSE 2.14%)
  })

  model({
    # Individual parameters. Eq. 1 power form for the two body-weight effects
    # (70 kg reference, per Methods) and Eq. 2 log-additive form for the two
    # categorical effects. Note that CL/F, V1/F, Q/F and V2/F are all APPARENT
    # (oral) quantities: absolute bioavailability of APX3330 was never
    # determined, so no separate F is estimated and none is anchored here.
    tlag <- exp(ltlag + e_fed_tlag * FED + etaltlag)
    ka   <- exp(lka + etalka)
    cl   <- exp(lcl + e_cancer_cl * DIS_CANCER + etalcl) * (WT / 70)^e_wt_cl
    vc   <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q    <- exp(lq)
    vp   <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with first-order absorption from a depot.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Total-quinone APX3330 plasma concentration in mg/L == ug/mL (dose in mg,
    # volumes in L). Proportional residual error.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
