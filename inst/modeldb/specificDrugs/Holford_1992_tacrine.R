Holford_1992_tacrine <- function() {
  description <- "Population pharmacodynamic disease-progression model for the cognitive subscale of the Alzheimer's Disease Assessment Scale (ADAS-cog, 0-70 score) in patients with probable Alzheimer's disease treated with tacrine. Linear disease progression (baseline S0 + alpha*time) with a tacrine effect on the location of the progression curve (effect compartment driven by IBW-normalised daily dose rate, no estimable PK clearance because the response is slow relative to the 2-hour tacrine plasma half-life) and a placebo effect with asymmetric onset / elimination / tolerance dynamics (placebo response builds up during treatment, dissipates after treatment ends, and develops tolerance during continued treatment). Estimated by Holford and Peace 1992 on 909 patients (5253 ADAS-cog observations) pooled from two clinical trials of identical design: US protocol 970-01 (n = 632) and French protocol 970-04 (n = 277). The French cohort takes multiplicative scale factors on baseline status (FS04 = 1.08), placebo potency (Fpp4 = 1.76), and placebo elimination half-time (Ft1/2el-p4 = 2.78). Inter-individual variability is correlated across baseline S0, progression rate alpha, and tacrine potency beta_a (block of three) with diagonal IIV on placebo potency beta_p; the time constants of the effect compartments are typical-value only. Residual error is proportional. NOTE: the lead Holford 1992 PNAS 89:11471-11475 'Results and validation' paper supplies all parameter values but the exact ODE form of the placebo dynamics is described in the companion methodology paper (PNAS 89:11466-11470) which was not available on disk at extraction time; the ODE form here is the field-standard reconstruction (asymmetric on/off placebo compartment plus multiplicative tolerance) and is documented in the validation vignette's Assumptions and deviations section."
  reference <- paste(
    "Holford NHG, Peace KE. (1992).",
    "Results and validation of a population pharmacodynamic model for",
    "cognitive effects in Alzheimer patients treated with tacrine.",
    "Proc Natl Acad Sci USA 89(23):11471-11475.",
    "doi:10.1073/pnas.89.23.11471.",
    "Companion methodology paper:",
    "Holford NHG, Peace KE. (1992).",
    "Methodologic aspects of a population pharmacodynamic model for",
    "cognitive effects in Alzheimer patients treated with tacrine.",
    "Proc Natl Acad Sci USA 89(23):11466-11470.",
    "doi:10.1073/pnas.89.23.11466.",
    sep = " "
  )
  vignette <- "Holford_1992_tacrine"
  units <- list(
    time          = "day",
    dosing        = "mg/day (tacrine daily dose rate supplied as a time-varying covariate column DOSE; no PK ODE)",
    concentration = "(ADAS-cog total cognitive subscale score, 0-70, unitless)"
  )

  covariateData <- list(
    IBW = list(
      description        = "Ideal body weight (Holford-Peace 1992 Devine variant: men IBW (kg) = 52 + 0.75 * (height_cm - 152); women IBW (kg) = 49 + 0.67 * (height_cm - 152)). Used to size-normalise the tacrine dose-rate input to the effect compartment.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline per subject. Reference value 60 kg is the population mean IBW from the source paper Data section. Users without IBW pre-computed should derive it from HT and SEXF using the Holford-Peace 1992 formula above before passing the column to the model.",
      source_name        = "IBW"
    ),
    DOSE = list(
      description        = "Time-varying tacrine daily dose rate (mg/day). The source paper used dose rate, not a PK ODE, because tacrine clearance could not be estimated directly from the cognitive-endpoint data (the ADAS-cog response is slow relative to the ~2-hour plasma tacrine half-life). Set to 0 during placebo arms, run-in, and washout phases.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-record covariate updated step-wise across treatment phases (e.g., 0 -> 40 -> 80 -> 0 across a titration sequence). Reference value 80 mg/day is the highest active dose in the source trials; the tacrine potency parameter beta_a is reported per 80 mg/day so the model uses `DOSE / 80` as the relative dose multiplier. The use case matches the canonical DOSE entry's (b) interpretation (time-varying current administered dose feeding a derived exposure term without an explicit PK compartment).",
      source_name        = "DOSE"
    ),
    TRT_PHASE = list(
      description        = "Active double-blind treatment-phase indicator (1 = on placebo or active treatment, 0 = baseline / run-in / washout / off-treatment). Gates both the placebo build-up and the placebo tolerance development.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (off-treatment / baseline)",
      notes              = "Holford 1992 explicitly states the placebo response is the same magnitude during active treatment and during placebo treatment ('the placebo efficacy of active drug is the same as placebo alone'); both treatment phases drive the placebo effect compartment. TRT_PHASE = 1 during any time the subject is on either placebo or active tacrine; 0 during pre-randomisation baseline and post-washout follow-up.",
      source_name        = "TRT_PHASE"
    ),
    REGION_FRANCE = list(
      description        = "Protocol indicator: 1 = French protocol 970-04, 0 = US protocol 970-01. Drives multiplicative scale factors on baseline ADAS-cog (FS04 = 1.08), placebo potency (Fpp4 = 1.76), and placebo elimination half-time (Ft1/2el-p4 = 2.78).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (US protocol 970-01)",
      notes              = "Time-fixed per subject. The source-paper indicator was the protocol number (PROT in 1, 4); REGION_FRANCE = as.integer(PROT == 4). The paper attributes the France-vs-US placebo-response difference to 'cultural milieu, which modify the behavioural response to drugs but are independent of pharmacological activity' (Discussion).",
      source_name        = "PROT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 909L,
    n_studies      = 2L,
    n_observations = 5253L,
    age_range      = "Adults / elderly with probable Alzheimer's disease (specific age range not tabulated in the source 'Results and validation' paper; demographic detail lives in the companion methodology paper Holford and Peace 1992 PNAS 89:11466-11470 which was not on disk at extraction time).",
    age_median     = "(not reported in the on-disk source paper)",
    weight_range   = "(not tabulated; ideal body weight mean across the cohort was 60 kg per source paper Data section)",
    weight_median  = "(not reported; population mean IBW = 60 kg)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Probable Alzheimer's disease (Methods refers to the standard NINCDS-ADRDA 'probable AD' diagnostic criteria implied by the 970-01 / 970-04 protocol designs).",
    dose_range     = "0 mg/day (placebo), 40 mg/day, and 80 mg/day oral tacrine (the 970-01 and 970-04 trials used titration sequences across these three dose levels with placebo intervals).",
    regions        = "United States (protocol 970-01, n = 632) and France (protocol 970-04, n = 277).",
    notes          = "Pooled analysis of two clinical trials of identical design conducted under Parke-Davis sponsorship (US trial 970-01 led by Davis / Thal and the Tacrine Collaborative Study group; French trial 970-04 led by Forette and the French Tacrine Study group). Outcome was ADAS-cog (cognitive subscale of Alzheimer's Disease Assessment Scale, 0-70). The model accommodates titration sequences (placebo / 40 / 80 mg/day in three orderings labeled tacseq114, tacseq214, tacseq314) and does not require all patients to complete all phases; the full tacalll4 analysis dataset pools all phases and protocols. Observation period up to 5 months per subject."
  )

  ini({
    # Disease-progression structural parameters
    # Source: Holford 1992 Table 2 'Full data set parameter estimates'
    # S0 (baseline ADAS-cog) = 28.7 units, SE 0.44
    # alpha (progression rate) = 6.17 units/year, SE 1.27
    lS0     <- log(28.7);  label("Baseline ADAS-cog at study entry (units, 0-70 score)")  # Table 2, Disease class, S_o units
    lalpha  <- log(6.17);  label("Disease-progression rate (ADAS-cog units/year)")        # Table 2, Disease class, a units/year

    # Tacrine pharmacodynamic structural parameter
    # beta_a = -2.99 units per 80 mg/day at IBW = 60 kg, SE 0.67
    # Stored as the absolute magnitude (the negative sign is applied in model()).
    lbeta_a <- log(2.99);  label("Absolute tacrine potency at steady state on 80 mg/day at IBW 60 kg (ADAS-cog units, sign negative)")  # Table 2, Pharmacodynamic class, beta_a units/80 mg/day

    # Placebo pharmacodynamic structural parameter
    # beta_p = -1.42 units, SE 0.20 (US cohort; France gets the multiplicative Fpp4 scaling)
    lbeta_p <- log(1.42);  label("Absolute placebo potency at full effect (US cohort baseline; ADAS-cog units, sign negative)")  # Table 2, Pharmacodynamic class, beta_p units

    # Effect-compartment time constants (typical-value only; no IIV reported in Table 2)
    # t1/2,eq,a = 20.9 days (tacrine effect compartment equilibration half-time), SE 6.0
    # t1/2,eq,p = 1.58 days (placebo on-treatment build-up half-time), SE 0.56
    # t1/2,el,p = 61.0 days (placebo off-treatment elimination half-time, US baseline), SE 28.6
    # t1/2,tol,p = 13.5 days (placebo tolerance development half-time), SE 3.4
    lt12eqa  <- log(20.9); label("Tacrine effect-compartment equilibration half-time (days)")           # Table 2, Pharmacokinetic class, t1/2,eq,a days
    lt12eqp  <- log(1.58); label("Placebo on-treatment build-up (equilibration) half-time (days)")     # Table 2, Pharmacokinetic class, t1/2,eq,p days
    lt12elp  <- log(61.0); label("Placebo off-treatment elimination half-time (US cohort, days)")      # Table 2, Pharmacokinetic class, t1/2,el,p days
    lt12tolp <- log(13.5); label("Placebo tolerance-development half-time (days)")                      # Table 2, Pharmacokinetic class, t1/2,tol,p day

    # Protocol scale factors (France vs US reference)
    # Source: Holford 1992 Table 2 'Scale' class. Reported as multiplicative scalars
    # F_x4; encoded here as fractional shifts e_region_france_<param> = F_x4 - 1.
    # FS0_4         = 1.08, SE 0.03  -> e_region_france_s0     = 0.08
    # F_betap_4     = 1.76, SE 0.25  -> e_region_france_betap  = 0.76
    # F_t1/2,el,p,4 = 2.78, SE 1.09  -> e_region_france_t12elp = 1.78
    # All three are estimated (have SEs), not fixed.
    e_region_france_s0     <- 0.08; label("Fractional shift on baseline ADAS-cog S0 in France (REGION_FRANCE = 1) vs US reference")        # Table 2, Scale class, FS0,4
    e_region_france_betap  <- 0.76; label("Fractional shift on placebo potency beta_p in France vs US reference")                            # Table 2, Scale class, F_beta-p,4
    e_region_france_t12elp <- 1.78; label("Fractional shift on placebo elimination half-time t1/2,el,p in France vs US reference")           # Table 2, Scale class, F_t1/2,el,p,4

    # Inter-individual variability
    # Source: Holford 1992 Table 2 'Population CV, %' column and Results text on correlations.
    # Block on (lS0, lalpha, lbeta_a): correlations 0.31 / 0.37 / 0.52 from Results paragraph
    # ("the correlation of baseline ADASC score ... with the potency of tacrine (beta_a) was 0.31;
    # the correlation of S0 with the ADAS-cog progression rate (alpha) was 0.37, while the
    # correlation of beta_a with alpha was 0.52").
    # Variances on the log scale via omega^2 = log(1 + CV^2):
    #   S0     CV 37.7%  -> omega^2 = log(1 + 0.377^2) = 0.1331
    #   alpha  CV 208%   -> omega^2 = log(1 + 2.08^2)  = 1.6730
    #   beta_a CV 126%   -> omega^2 = log(1 + 1.26^2)  = 0.9508
    # Covariances cov(x, y) = corr(x, y) * sqrt(var(x) * var(y)):
    #   cov(S0,alpha)   = 0.37 * sqrt(0.1331 * 1.6730) = 0.1746
    #   cov(S0,beta_a)  = 0.31 * sqrt(0.1331 * 0.9508) = 0.1102
    #   cov(alpha,beta_a) = 0.52 * sqrt(1.6730 * 0.9508) = 0.6562
    etalS0 + etalalpha + etalbeta_a ~ c(0.1331,
                                        0.1746, 1.6730,
                                        0.1102, 0.6562, 0.9508)

    # Diagonal IIV on placebo potency
    #   beta_p CV 128%   -> omega^2 = log(1 + 1.28^2)  = 0.9702
    etalbeta_p ~ 0.9702

    # Residual error: proportional ("Error SD ADASC 3.14", SE 0.08).
    # The paper reports the residual SD on the ADAS-cog scale at typical observation;
    # interpreted as propSd ~= 3.14 / typical_ADASC (~30) = 0.105 since the prose
    # explicitly says proportional error is better than additive. See vignette
    # Assumptions and deviations for the interpretation rationale.
    propSd <- 0.105; label("Proportional residual error (fraction; implied from reported SD ADASC = 3.14 at typical ADAS-cog ~30)")  # Table 2, Error class, SD ADASC
  })

  model({
    # Time-conversion: rates published in /year (alpha) and /day (effect compartments).
    # Internal time variable t is in days (matches the source data's observation grid
    # of up to ~5 months = ~150 days).
    # alpha needs conversion: 6.17 units/year / 365.25 days/year = 0.01689 units/day.

    # 1. Individual structural parameters with log-normal IIV and protocol scaling
    S0_indiv     <- exp(lS0     + etalS0)     * (1 + e_region_france_s0     * REGION_FRANCE)
    alpha_yr     <- exp(lalpha  + etalalpha)
    alpha_day    <- alpha_yr / 365.25
    beta_a_indiv <- -exp(lbeta_a + etalbeta_a)
    beta_p_indiv <- -exp(lbeta_p + etalbeta_p) * (1 + e_region_france_betap  * REGION_FRANCE)
    t12elp_indiv <- exp(lt12elp)              * (1 + e_region_france_t12elp * REGION_FRANCE)

    # 2. Effect-compartment rate constants (1/day; typical-value, no IIV)
    keqa  <- log(2) / exp(lt12eqa)
    keqp  <- log(2) / exp(lt12eqp)
    kelp  <- log(2) / t12elp_indiv
    ktolp <- log(2) / exp(lt12tolp)

    # 3. Tacrine effect-compartment input
    # Source paper: tacrine response is proportional to the average steady-state
    # tacrine concentration, which is dose-rate / clearance. Clearance scales with
    # IBW (the source identified IBW as the best size descriptor over total body
    # weight or height). The unobservable absolute clearance is absorbed into
    # beta_a (defined at IBW = 60 kg and DOSE = 80 mg/day), leaving a relative
    # dose-rate input scaled by 60 / IBW.
    Dr_norm <- (DOSE / 80) * (60 / IBW)

    # 4. ODEs for the three effect / tolerance states. The compartments use
    # the canonical numbered `effect<n>` chain:
    #   effect1 -> tacrine effect (paper symbol Cea)
    #   effect2 -> placebo effect (paper symbol Cep)
    #   effect3 -> placebo tolerance (paper symbol Tolp; semantically a
    #              delayed component of the placebo system, encoded as a
    #              third numbered effect compartment for convention compliance).
    #
    # effect1 (tacrine) drives toward Dr_norm with rate keqa.
    d/dt(effect1) <- keqa * (Dr_norm - effect1)
    effect1(0)    <- 0

    # effect2 (placebo) has asymmetric on/off kinetics: when TRT_PHASE = 1 it
    # builds up toward 1 with rate keqp (t1/2 = 1.58 d); when TRT_PHASE = 0 it
    # dissipates with rate kelp (t1/2 = 61 d in the US cohort, 61 * 2.78 =
    # 170 d in the France cohort). See the vignette Assumptions and deviations
    # section for the ODE-form reconstruction rationale -- the companion paper
    # Holford and Peace 1992 PNAS 89:11466-11470 holds the exact form.
    d/dt(effect2) <- TRT_PHASE * keqp * (1 - effect2) - (1 - TRT_PHASE) * kelp * effect2
    effect2(0)    <- 0

    # effect3 (placebo tolerance) follows TRT_PHASE with first-order kinetics
    # (rate ktolp, t1/2 = 13.5 d): builds toward 1 during treatment, decays
    # toward 0 off-treatment. The effective placebo is attenuated
    # multiplicatively by (1 - effect3).
    d/dt(effect3) <- ktolp * (TRT_PHASE - effect3)
    effect3(0)    <- 0

    # 5. ADAS-cog observation (paper-named PD output, kept as ADAS_cog rather
    # than the canonical PK observation name Cc; see vignette Assumptions and
    # deviations). Linear disease progression (S0 + alpha*t) shifted by the
    # time-developing tacrine effect (beta_a * effect1) and the placebo
    # response with tolerance attenuation (beta_p * effect2 * (1 - effect3)).
    ADAS_cog <- S0_indiv + alpha_day * t + beta_a_indiv * effect1 + beta_p_indiv * effect2 * (1 - effect3)

    ADAS_cog ~ prop(propSd)
  })
}
