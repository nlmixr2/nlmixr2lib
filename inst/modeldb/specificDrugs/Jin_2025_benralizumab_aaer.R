Jin_2025_benralizumab_aaer <- function() {
  description <- "Longitudinal Poisson exposure-response model for the asthma exacerbation rate under benralizumab in Chinese and non-Chinese patients with severe eosinophilic asthma (Jin 2025), coupling the Jin 2025 two-compartment subcutaneous popPK layer to an Emax drug effect on the log exacerbation rate, with prior-exacerbation count, maintenance oral corticosteroid use and Eastern-European region on the baseline rate and Chinese race on Emax"
  reference <- paste(
    "Jin Y, Guiastrennec B, Stuke M, Yao Y, Zhang Y, Barker P, Jison M,",
    "Penland RC, Ding J, Lukka PB.",
    "Population pharmacokinetics and exposure-response analysis of benralizumab in",
    "Chinese adults, adolescents, and pediatric participants with severe eosinophilic asthma.",
    "Clin Pharmacokinet. 2025;64:1233-1245. doi:10.1007/s40262-025-01538-9.",
    "Model equation from Resource 3 and parameter values from Resource 17 of the electronic",
    "supplementary material; the PK layer is reproduced from Resource 10 -- see",
    "modellib('Jin_2025_benralizumab').",
    "Updates the legacy exposure-response model of Chia YL, Yan L, Yu B, et al.",
    "Clin Pharmacol Ther. 2019;106:383-90; doi:10.1002/cpt.1371.",
    sep = " "
  )
  vignette <- "Jin_2025_benralizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot       = list(analyte = "benralizumab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "benralizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "benralizumab", units = "mg", specimen = "serum", verified = TRUE),
    cumhaz_exac = list(analyte = "expected cumulative asthma exacerbation count", units = "events", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "PK-layer covariate. Power effects normalised to 70 kg: (WT/70)^0.849 on CL, (WT/70)^0.799 on",
        "V2 and (WT/70)^0.639 on V3 (Jin 2025 Resource 10). Body weight was NOT tested on the",
        "exacerbation-rate parameters. Source column BWGT."
      ),
      source_name        = "BWGT"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator, 1 = Asian (including all Chinese participants), 0 = non-Asian",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = paste(
        "PK-layer covariate only: CL * exp(0.0952 * RACE_ASIAN), a 9.99% higher clearance",
        "(Jin 2025 Resource 10). Set RACE_ASIAN = 1 whenever RACE_CHINESE = 1, since Jin 2025",
        "Methods 2.1 classifies all Chinese subjects as Asian. Source column ASIAN."
      ),
      source_name        = "ASIAN"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody status, 1 = positive ADA titer, 0 = no ADA",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no anti-drug antibody detected)",
      notes              = paste(
        "PK-layer covariate only: CL * exp(0.762 * ADA_POS), a 114% higher clearance",
        "(Jin 2025 Resource 10). Source column ADA."
      ),
      source_name        = "ADA"
    ),
    NEXAC12M = list(
      description        = "Number of asthma exacerbations in the 12 months before study entry",
      units              = "count (events in the prior 12 months)",
      type               = "count",
      reference_category = NULL,
      notes              = paste(
        "Log-linear effect on the baseline exacerbation rate: + 0.17 * NEXAC12M inside the exponent",
        "(Jin 2025 Resource 17, 'beta_Base, PE -- Prior exacerbation on baseline' = 0.17, RSE 9.18%).",
        "Used as a RAW, UNCENTERED continuous count multiplied by a single fitted slope, exactly as",
        "the source fits it, so it must not be decomposed into band indicators; the corresponding",
        "empirical generalized linear model in Jin 2025 Resource 3 likewise enters COVAR13, 'the",
        "number of exacerbation in the last 12 months', as one linear term. This means the tabulated",
        "Base intercept (-6.64) is the log rate at NEXAC12M = 0, which is outside the enrolled range:",
        "SIROCCO, CALIMA and MIRACLE all required at least two exacerbations in the prior year.",
        "Source column PE."
      ),
      source_name        = "PE"
    ),
    CONMED_STEROID = list(
      description        = "Maintenance oral corticosteroid use at baseline, 1 = on maintenance OCS, 0 = not on maintenance OCS",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no maintenance oral corticosteroid use at baseline)",
      notes              = paste(
        "Log-linear effect on the baseline exacerbation rate: + 0.345 * CONMED_STEROID inside the",
        "exponent, i.e. exp(0.345) = 1.41-fold higher baseline rate (Jin 2025 Resource 17,",
        "'beta_Base, OCS -- Maintenance OCS on baseline' = 0.345, RSE 25%). Time-fixed per subject",
        "(baseline / chronic background use), which is the first of the two temporal grains the",
        "CONMED_STEROID register entry supports and the one it names severe asthma for. The",
        "corresponding empirical model in Jin 2025 Resource 3 uses COVAR04, 'the use of oral",
        "corticosteroids at baseline (1: use, 2: no use)'; note that the empirical model's coding is",
        "inverted relative to the canonical orientation used here (canonical 1 = use). Source column OCS."
      ),
      source_name        = "OCS"
    ),
    REGION_EASTEUROPE = list(
      description        = "Eastern European study-site indicator, 1 = enrolled at an Eastern European site, 0 = enrolled elsewhere",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-Eastern-European region: Asia, Europe, North America, or rest of world)",
      notes              = paste(
        "Log-linear effect on the baseline exacerbation rate: - 0.329 * REGION_EASTEUROPE inside the",
        "exponent, i.e. exp(-0.329) = 0.72-fold (28% lower) baseline rate at Eastern European sites",
        "(Jin 2025 Resource 17, 'beta_Base, EE -- Eastern Europe on baseline' = -0.329, RSE 22.6%).",
        "Jin 2025 Resource 3 defines the parent region variable REGION01 with five mutually exclusive",
        "levels -- '1: Asia, 2: Eastern Europe, 3: Europe, 4: North America, 5: Rest of the world' --",
        "so Eastern Europe is a level DISTINCT from (western) Europe and must not be folded into",
        "REGION_EUROPE. Only the Eastern Europe contrast was retained in the final longitudinal model.",
        "Time-fixed per subject. Source column EE."
      ),
      source_name        = "EE"
    ),
    RACE_CHINESE = list(
      description        = "Chinese-heritage race indicator, 1 = participant from mainland China or Taiwan, 0 = non-Chinese",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Chinese)",
      notes              = paste(
        "Scales the maximum treatment effect FRACTIONALLY rather than on a log scale:",
        "Emax * (1 + 1.27 * RACE_CHINESE), so Emax goes from -0.51 in non-Chinese participants to",
        "-0.51 * 2.27 = -1.158 in Chinese participants (Jin 2025 Resource 17,",
        "'beta_Emax, CHINESE -- CHINESE on Emax' = 1.27, RSE 27.3%). The (1 + beta) form is what",
        "reproduces the main text's statement that 'the maximal treatment effect significantly",
        "increased (+127%; p < 0.001)' in Chinese participants: the coefficient 1.27 IS the +127%.",
        "An exp(beta) reading would give +256% and an additive Emax + beta reading would flip the sign",
        "of the treatment effect, so both are excluded. Jin 2025 Methods 2.1: 'All subjects from",
        "mainland China and Taiwan were considered as Chinese.' Time-fixed per subject.",
        "Source column CHINESE."
      ),
      source_name        = "CHINESE"
    )
  )

  covariatesDataExcluded <- list(
    STUDY_MICP220 = list(
      description = "Study MI-CP220 indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Relevant only to the popPK layer, where it selects a study-specific absolute subcutaneous",
        "bioavailability -- see modellib('Jin_2025_benralizumab'). MI-CP220 is a phase II study and",
        "contributed no data to the exposure-response analysis, which used only the three phase III",
        "studies SIROCCO, CALIMA and MIRACLE (Jin 2025 Methods 2.1), so the stratum is omitted here",
        "and the reference bioavailability applies to every subject."
      )
    ),
    STUDY_AMES = list(
      description = "AMES study indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Relevant only to the popPK layer -- see modellib('Jin_2025_benralizumab'). AMES is a phase I",
        "healthy-volunteer study and contributed no data to the exposure-response analysis."
      )
    ),
    DOSE_HIGH = list(
      description = "200 mg subcutaneous dose indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Relevant only to the popPK layer -- see modellib('Jin_2025_benralizumab'). The three phase",
        "III studies in the exposure-response dataset used 30 mg subcutaneously only, so DOSE_HIGH is",
        "0 for every subject here and the effect is omitted."
      )
    ),
    RACE_ASIAN_ER = list(
      description = "Asian race indicator tested on the exposure-response parameters",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Jin 2025 Methods 2.5.1 evaluated the Asian covariate on both the exacerbation-rate baseline",
        "and Emax alongside the Chinese covariate, and Resource 6 / Resource 19 report the Asian",
        "steady-state AAER-ratio simulations (Asian 0.4878 [95% CI 0.4341, 0.5585] versus non-Asian",
        "0.6419 [0.5887, 0.7016]). The FINAL longitudinal model in Resource 17 retains only the",
        "CHINESE effect on Emax; no Asian exposure-response coefficient is tabulated anywhere in the",
        "paper or its supplement, so no Asian term can be encoded. Note that RACE_ASIAN IS used in",
        "this file, but as a PK-layer covariate on clearance only."
      )
    ),
    TRT_BENRALIZUMAB = list(
      description = "Benralizumab treatment-arm indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Not needed as an explicit covariate: the drug effect enters only through the",
        "concentration-driven Emax term, which is identically zero when no benralizumab is dosed",
        "(Cc = 0). The same model therefore serves the placebo and the benralizumab arms, which is",
        "how Jin 2025 fits all three treatment groups -- placebo, 30 mg Q4W and 30 mg Q8W -- in one",
        "longitudinal model."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 3205L,
    n_studies      = 3L,
    age_range      = "12-75 years",
    age_median     = "not published for the exposure-response subset",
    weight_range   = "not published for the exposure-response subset",
    weight_median  = "not published for the exposure-response subset",
    sex_female_pct = NA_real_,
    race_ethnicity = "Chinese and non-Chinese strata; the exposure-response dataset is the SIROCCO + CALIMA + MIRACLE pool, in which the Chinese participants come from MIRACLE.",
    disease_state  = "Severe, uncontrolled eosinophilic asthma on medium-to-high-dose inhaled corticosteroid plus a long-acting beta2-agonist, with benralizumab or placebo as add-on maintenance therapy.",
    dose_range     = "Placebo, benralizumab 30 mg subcutaneously every 4 weeks, or benralizumab 30 mg subcutaneously every 8 weeks (first three doses every 4 weeks).",
    regions        = "Multi-regional. Five geographic regions were defined (Asia, Eastern Europe, Europe, North America, rest of world); only the Eastern Europe contrast was retained.",
    notes          = paste(
      "Jin 2025 Methods 2.5: 'Data from three phase III studies (SIROCCO, CALIMA, and MIRACLE) were",
      "used in the ER analysis'. THE PAPER DOES NOT PRINT AN ANALYSIS-SET N for either",
      "exposure-response model, so n_subjects is the combined RANDOMISED total of those three",
      "studies (SIROCCO 1204 + CALIMA 1306 + MIRACLE 695 = 3205) and is an UPPER BOUND rather than",
      "an exact analysis-set size. The only exposure-response counts the paper enumerates are the",
      "per-quartile Q8W benralizumab recipients in Resource 15 (non-Chinese 196 + 198 + 201 + 193 =",
      "788 and Chinese 65 + 62 + 54 + 54 = 235, i.e. 1023 Q8W recipients), to which the placebo and",
      "30 mg Q4W arms must be added but are not tabulated. See the vignette Errata.",
      "The model was fitted to exacerbation events summarised in 8-WEEK intervals as a Poisson count",
      "(Jin 2025 Methods 2.5.1 and Resource 3).",
      "VALIDITY RANGE, stated verbatim in Jin 2025 Methods 2.5.1: 'the model is valid for doses",
      "leading to mean average concentrations above 500 ng/mL and mean Ctrough above 50 ng/mL.'",
      "Do not use it to extrapolate below that exposure.",
      "EC50 (1.76 ng/mL) is far below the assay LLOQ of 3.86 ng/mL and is estimated with 180% RSE,",
      "so the exposure-response curve is effectively FLAT and at its plateau across the whole 30 mg",
      "dose range; Jin 2025 attributes this to most subjects sitting in the efficacy-plateau region.",
      "Emax is nevertheless estimated precisely (12.3% RSE).",
      "The individual benralizumab concentrations driving the Emax term were the popPK model's",
      "individual (empirical Bayes) predictions; the PK layer is reproduced inline here so the model",
      "is self-contained."
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # PK layer, reproduced from Jin 2025 Resource 10 so this model is
    # self-contained. Jin 2025 fitted the exposure-response model conditional
    # on individual popPK predictions ("The model used individual parameters to
    # simulate benralizumab exposure over time", Methods 2.5.1), so every PK
    # value below is FIXED from that popPK step rather than estimated here.
    # See modellib('Jin_2025_benralizumab') for the full annotation of each
    # value, including the study- and dose-specific bioavailability strata that
    # are dropped here because the exposure-response dataset contains only the
    # three phase III studies at 30 mg SC.
    # Reference covariate values: 70 kg, non-Asian, ADA-negative.
    # ------------------------------------------------------------------------
    lcl     <- fixed(log(0.269));         label("Systemic elimination clearance CL (L/day)")                        # Jin 2025 Resource 10, "CL (L/day)" = 0.269
    lvc     <- fixed(log(3.02));          label("Central volume of distribution V2 (L)")                            # Jin 2025 Resource 10, "V2 (L)" = 3.02
    lq      <- fixed(log(1.05));          label("Intercompartmental clearance Q2 (L/day)")                          # Jin 2025 Resource 10, "Q2 (L/day)" = 1.05
    lvp     <- fixed(log(2.67));          label("Peripheral volume of distribution V3 (L)")                         # Jin 2025 Resource 10, "V3 (L)" = 2.67
    lka     <- fixed(log(log(2) / 3.02)); label("First-order subcutaneous absorption rate ka (1/day)")              # Jin 2025 Resource 10, "KAThalf (day)" = 3.02 -> ka = log(2)/3.02
    lfdepot <- fixed(log(0.539));         label("Absolute subcutaneous bioavailability Fa1 (fraction)")             # Jin 2025 Resource 10, "Fa1 (fraction)" = 0.539

    e_wt_cl    <- fixed(0.849);   label("Power exponent of WT/70 on CL (unitless)")                                 # Jin 2025 Resource 10, "Beta_CL, BWGT (kg)"
    e_wt_vc    <- fixed(0.799);   label("Power exponent of WT/70 on V2 (unitless)")                                 # Jin 2025 Resource 10, "Beta_V2, BWGT (kg)"
    e_wt_vp    <- fixed(0.639);   label("Power exponent of WT/70 on V3 (unitless)")                                 # Jin 2025 Resource 10, "Beta_V3, BWGT (kg)"
    e_asian_cl <- fixed(0.0952);  label("Log-scale effect of Asian race on CL (unitless)")                          # Jin 2025 Resource 10, "Beta_CL, ASIAN_1"
    e_ada_cl   <- fixed(0.762);   label("Log-scale effect of ADA positivity on CL (unitless)")                      # Jin 2025 Resource 10, "RCLADA"

    etalcl     ~ fixed(0.053824)  # Jin 2025 Resource 10: omega(CL) SD 0.232 -> 0.232^2
    etalvc     ~ fixed(0.077284)  # Jin 2025 Resource 10: omega(V2) SD 0.278 -> 0.278^2
    etalq      ~ fixed(0.007921)  # Jin 2025 Resource 10: omega(Q2) SD 0.089, held constant by the authors -> 0.089^2
    etalvp     ~ fixed(0.184041)  # Jin 2025 Resource 10: omega(V3) SD 0.429 -> 0.429^2
    etalka     ~ fixed(0.499849)  # Jin 2025 Resource 10: omega(KAThalf) SD 0.707 -> 0.707^2
    etalfdepot ~ fixed(0.075076)  # Jin 2025 Resource 10: omega(Fa1) SD 0.274 -> 0.274^2

    # ------------------------------------------------------------------------
    # Longitudinal asthma exacerbation rate layer. Jin 2025 Resource 3 prints
    # the model as
    #
    #   lambda_j(a,b) = integral_a^b exp( beta0 + beta * X_j
    #                                     + Emax_j * C_j(t) / (C_j(t) + EC50)
    #                                     + eta_j ) dt
    #
    # i.e. the expected count in an 8-week interval is the time integral of an
    # instantaneous rate whose log is linear in the covariates and carries a
    # concentration-driven Emax term. All values below are the final-model
    # estimates in Jin 2025 Resource 17.
    #
    # Emax is NEGATIVE: it reduces the log exacerbation rate, so exp(Emax) at
    # full effect is the achievable rate ratio versus placebo.
    # ------------------------------------------------------------------------
    lbase <- -6.64;  label("Baseline log asthma exacerbation rate at zero prior exacerbations (log events/day)")  # Jin 2025 Resource 17, "Base -- Asthma exacerbation rate" = -6.64 (RSE 1.14%)
    emax  <- -0.51;  label("Maximum benralizumab effect on the log exacerbation rate (unitless)")                 # Jin 2025 Resource 17, "Emax -- Maximum treatment effect" = -0.51 (RSE 12.3%)

    # EC50 is tabulated in ng/mL; this model carries concentration in mg/L
    # (units$concentration above), so the printed 1.76 ng/mL is divided by 1000.
    # 1.76 ng/mL == 0.00176 mg/L. This is a pure unit conversion, not a refit.
    lec50 <- log(1.76 / 1000);  label("Benralizumab concentration giving half-maximal effect EC50 (mg/L)")        # Jin 2025 Resource 17, "EC50" = 1.76 ng/mL (RSE 180%) -> 0.00176 mg/L

    # Covariate effects on the baseline log rate (additive inside the exponent,
    # matching the beta * X_j term of the Resource 3 equation).
    e_nexac12m_base <-  0.17;   label("Log-scale slope of the prior-exacerbation count on the baseline rate (per event)")  # Jin 2025 Resource 17, "beta_Base, PE" = 0.17 (RSE 9.18%)
    e_ocs_base      <-  0.345;  label("Log-scale effect of maintenance OCS use on the baseline rate (unitless)")          # Jin 2025 Resource 17, "beta_Base, OCS" = 0.345 (RSE 25%)
    e_ee_base       <- -0.329;  label("Log-scale effect of Eastern European region on the baseline rate (unitless)")       # Jin 2025 Resource 17, "beta_Base, EE" = -0.329 (RSE 22.6%)

    # Covariate effect on Emax. FRACTIONAL, not log-scale: Emax * (1 + beta).
    e_chinese_emax <- 1.27;  label("Fractional increase in the magnitude of Emax in Chinese participants (unitless; +127%)")  # Jin 2025 Resource 17, "beta_Emax, CHINESE" = 1.27 (RSE 27.3%)

    # Inter-individual variability on the baseline log rate. Resource 17 gives
    # the comment "Normally distributed", i.e. eta is additive on the log rate,
    # which is the eta_j term of the Resource 3 equation. The tabulated 1.04 is
    # read as the SD on that log scale (the value sits in the same column that
    # Resources 8, 10 and 11 state holds standard deviations), so the variance
    # is 1.04^2.
    etalbase ~ 1.0816  # Jin 2025 Resource 17: omega(Base) 1.04 (RSE 3.56%, shrinkage 37.7%) -> 1.04^2
  })
  model({
    # --- 1. PK layer (see modellib('Jin_2025_benralizumab')) ----------------
    cl <- exp(lcl + etalcl) *
          (WT / 70)^e_wt_cl *
          exp(e_asian_cl * RACE_ASIAN) *
          exp(e_ada_cl * ADA_POS)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp
    q  <- exp(lq  + etalq)
    ka <- exp(lka + etalka)
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    Cc <- central / vc

    # --- 2. Exacerbation-rate layer ----------------------------------------
    # Baseline log rate: intercept + covariates + additive individual random
    # effect (Jin 2025 Resource 3 beta0 + beta * X_j + eta_j).
    baselograte <- lbase + etalbase +
                   e_nexac12m_base * NEXAC12M +
                   e_ocs_base * CONMED_STEROID +
                   e_ee_base * REGION_EASTEUROPE

    # Maximum treatment effect, scaled fractionally in Chinese participants.
    emaxi <- emax * (1 + e_chinese_emax * RACE_CHINESE)

    # Concentration-driven Emax term. With no benralizumab on board Cc is 0 and
    # drugeff is 0, so the placebo arm falls out of the same expression.
    ec50 <- exp(lec50)
    drugeff <- emaxi * Cc / (Cc + ec50)

    # Instantaneous exacerbation rate, events/day.
    exacrate <- exp(baselograte + drugeff)

    # --- 3. Interval count via the cumulative rate integral -----------------
    # Resource 3 defines the Poisson mean over an interval (a, b) as the
    # integral of exacrate across it, so the integral is carried as a state.
    # `cumhaz_exac` is the expected cumulative number of exacerbations since t = 0;
    # the expected count in any interval is the increment of `cumhaz_exac` across
    # it, and the annualised rate is the increment over 365 days. Interval
    # counts over disjoint intervals are independent Poisson variables with
    # those increments as means.
    d/dt(cumhaz_exac) <- exacrate

    exac_count <- cumhaz_exac
    exac_count ~ pois(exac_count)
  })
}
