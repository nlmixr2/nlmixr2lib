Claret_2014_motesanib_tumorsize_OS <- function() {
  description <- paste(
    "Joint tumor-size (TS) / overall-survival (OS) exploratory model for",
    "first-line advanced nonsquamous non-small cell lung cancer (NSCLC),",
    "developed by Claret et al. on the phase III MONET1 study of",
    "carboplatin/paclitaxel (CP) plus motesanib 125 mg q.d. vs. CP plus",
    "placebo (934 TS-evaluable patients of 1,090 enrolled).",
    "The TS sub-model is the simplified, semi-mechanistic tumor growth",
    "inhibition (TGI) model TS(t) = BASE * exp(kge*t - (kdrug/lambda) *",
    "(1 - exp(-lambda*t))), encoded here as the equivalent ODE",
    "d/dt(tumor_size) = tumor_size * (kge - kdrug*exp(-lambda*t)); kge is",
    "the tumor growth rate, kdrug the tumor growth inhibition rate, and",
    "lambda the rate constant for treatment-effect decay (resistance and",
    "dose reductions). The OS sub-model is an accelerated failure time",
    "(AFT) log-normal survival regression whose treatment-effect predictor",
    "is time to tumor growth, TTG = (log(kdrug) - log(kge)) / lambda, the",
    "model-predicted time to nadir TS. TTG is computed inside model() from",
    "the individual TS parameters and enters as log(TTG + 100); baseline",
    "TS, Asian ethnicity and smoking history enter linearly on the",
    "log-survival-time scale. ERRATUM: Supplementary Table S1 prints the",
    "KDE and lambda ESTIMATES transposed; this file uses KDE = 0.0103 and",
    "lambda = 0.0193 per day rather than the printed 0.0193 and 0.0103.",
    "The correction is forced by the paper's own outputs -- because the TS",
    "fit is pooled across arms, its typical value must lie between the two",
    "published per-arm medians, which it does for both TTG (117 vs 93-125",
    "days) and week-8 TSratio (0.746 vs 0.69-0.77) only after the swap.",
    "Both candidate numbers are verbatim from Table S1; nothing is tuned.",
    "See the vignette for the full argument. IMPORTANT: the TS sub-model is a",
    "single pooled fit with NO treatment-arm term, so forward simulation",
    "with freshly drawn etas yields the same TS and TTG distribution for",
    "both arms; the paper's own trial simulations were conditioned on the",
    "observed per-patient TTG values rather than on fresh draws. See the",
    "vignette Assumptions and deviations section.",
    sep = " "
  )
  reference <- paste(
    "Claret L, Bruno R, Lu J-F, Sun Y-N, Hsu C-P.",
    "Exploratory modeling and simulation to support development of",
    "motesanib in Asian patients with non-small cell lung cancer based on",
    "MONET1 study results.",
    "Clin Pharmacol Ther. 2014;95(4):446-451.",
    "doi:10.1038/clpt.2014.11.",
    "TS model parameter estimates: Supplementary Table S1.",
    "OS model parameter estimates: Table 3;",
    "parameter covariance matrix: Supplementary Table S3.",
    "The companion model for TS-nonevaluable patients is",
    "modellib('Claret_2014_motesanib_OS_nonevaluable').",
    sep = " "
  )
  vignette <- "Claret_2014_motesanib_nsclc"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; the simplified TS model is not exposure-driven -- the motesanib 125 mg q.d. regimen enters only through the estimated kdrug and lambda)",
    concentration = "mm (tumor size, the sum of longest diameters of the target lesions; not a drug concentration)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Claret 2014 Methods ("TS denotes the
  # tumor size (sum of longest diameters of the target lesions, in mm)").
  compartmentData <- list(
    tumor_size = list(analyte = "tumor size (RECIST sum of longest diameters of target lesions)", units = "mm", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    RACE_ASIAN = list(
      description        = "1 = Asian patient, 0 = non-Asian. Independent baseline prognostic factor in the final OS model; Asian patients had longer OS that none of the other available baseline covariates explained.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = "Additive on the log-survival-time (AFT) scale: mu_os = ... + 0.319 * RACE_ASIAN, i.e. a exp(0.319) = 1.38-fold longer median OS in Asian patients at equal TTG, baseline TS and smoking history (Claret 2014 Table 3). 227 of the 1,090 MONET1 patients were Asian; 219 Asian patients were used for the phase III simulations, of whom 21 were TS-nonevaluable.",
      source_name        = "Asian ethnicity"
    ),
    SMOKE_NEVER = list(
      description        = "1 = never smoker at baseline, 0 = former or current smoker. Claret 2014 dichotomises smoking history the other way round (its indicator is 1 for former OR current smokers); the canonical register column is SMOKE_NEVER, so model() forms the paper's ever-smoker indicator as its complement, (1 - SMOKE_NEVER).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (former or current smoker) -- note this is the COMPLEMENT of the paper's reference category, which is the never-smoker group",
      notes              = "The published coefficient is -0.264 on an ever-smoker (former or current) indicator, relative to never smokers (Claret 2014 Table 3 and footnote d). Encoding it as e_smoke_ever_mu_os * (1 - SMOKE_NEVER) keeps BOTH the published intercept (1.079) and the published coefficient (-0.264) numerically intact, so they correspond directly to the Supplementary Table S3 covariance matrix. Re-parameterising onto SMOKE_NEVER directly would instead require shifting the intercept to 1.079 - 0.264 = 0.815 and would break that correspondence.",
      source_name        = "Smoking history"
    )
  )

  # Covariates screened in the Claret 2014 univariate Cox analysis (Table 2)
  # but NOT retained in the final multivariate OS model. Documented here to
  # preserve the covariate-screen provenance without carrying "declared but
  # not referenced" convention warnings. checkModelConventions() treats this
  # list as documentation only.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (female indicator).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Claret 2014 Table 2: significant in the univariate Cox analysis (P < 0.0001) but eliminated in the backward stepwise selection of the parametric multivariate model."
    ),
    HT = list(
      description = "Height.",
      units       = "cm",
      type        = "continuous",
      notes       = "Claret 2014 Table 2: univariate Cox P = 0.0005; not retained."
    ),
    BSA = list(
      description = "Body surface area.",
      units       = "m^2",
      type        = "continuous",
      notes       = "Claret 2014 Table 2: univariate Cox P = 0.1307 (not significant at P < 0.05); not carried into the full model."
    ),
    WHO_PS = list(
      description = "Eastern Cooperative Oncology Group (ECOG) performance status.",
      units       = "(integer score)",
      type        = "categorical",
      notes       = "Claret 2014 Table 2: univariate Cox P = 0.0022; not retained in the updated OS model, although it IS the sole prognostic factor in the historical first-line NSCLC OS model that this paper set out to improve on."
    ),
    ALB = list(
      description = "Baseline serum albumin.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Claret 2014 Table 2: univariate Cox P < 0.0001. Retained only in the laboratory-data version of the model (n = 848), which the authors did NOT adopt as final because it was developed on fewer patients (Results, 'Updated survival model development'). Baseline albumin IS the sole prognostic factor of the companion nonevaluable-patient model -- see modellib('Claret_2014_motesanib_OS_nonevaluable')."
    ),
    ALP = list(
      description = "Baseline alkaline phosphatase.",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Claret 2014 Table 2: univariate Cox P = 0.1083; not significant."
    ),
    TBILI = list(
      description = "Baseline total bilirubin.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Claret 2014 Table 2: univariate Cox P = 0.1153; not significant."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance.",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Claret 2014 Table 2: univariate Cox P = 0.2892; not significant."
    ),
    LDH = list(
      description = "Baseline lactate dehydrogenase.",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Claret 2014 Table 2: univariate Cox P < 0.0001. Retained only in the non-final laboratory-data version of the model (n = 848)."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase.",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Claret 2014 Table 2: univariate Cox P = 0.0276; not retained."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase.",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Claret 2014 Table 2: univariate Cox P = 0.0032; not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 934L,
    n_studies      = 1L,
    age_range      = "adults with advanced nonsquamous NSCLC; the age distribution is reported in the MONET1 primary publication (Scagliotti 2012, J Clin Oncol 30:2829-2836) rather than in Claret 2014",
    weight_range   = "not reported in Claret 2014",
    sex_female_pct = NA_real_,
    race_ethnicity = "227 of the 1,090 enrolled MONET1 patients (21 percent) were Asian; 219 Asian patients were resampled for the virtual phase III simulations, of whom 21 were TS-nonevaluable.",
    disease_state  = "Advanced (stage IIIB with pleural effusion or stage IV) or recurrent nonsquamous non-small cell lung cancer, previously untreated (first line).",
    dose_range     = "Carboplatin/paclitaxel (CP) plus motesanib 125 mg orally once daily, vs. CP plus placebo. The simplified TS model is not exposure-driven, so dose does not appear in the model equations.",
    regions        = "MONET1 was an international, randomized, placebo-controlled, double-blind phase III study (so described in the title of Claret 2014's reference 13, Scagliotti 2012, J Clin Oncol 30:2829-2836); the Asian subgroup is the subpopulation of interest in this paper. Claret 2014 does not print a trial-registry identifier, so none is recorded here.",
    biomarkers     = "Longitudinal tumor size (TS; RECIST sum of longest diameters of target lesions, mm) and overall survival (OS, days). Derived TS-response metrics: time to tumor growth (TTG, days) and the week-8 TS ratio (TSratio, unitless).",
    notes          = paste(
      "Of 1,090 enrolled patients, 934 (86 percent) were TS-evaluable: 453 of 541 in the motesanib arm and 481 of 549 in the placebo arm.",
      "The remaining 156 (14 percent) were nonevaluable and are covered by the companion model modellib('Claret_2014_motesanib_OS_nonevaluable').",
      "A median of four TS measurements were available per patient (range 2-27).",
      "",
      "Published TS-metric distributions, usable as validation targets for the TS sub-model (Claret 2014 Results, 'TS model'):",
      "  TSratio (TS at week 8 / baseline), median (2.5th-97.5th percentiles): 0.69 (0.24-0.99) motesanib arm; 0.77 (0.27-1.25) placebo arm.",
      "  TTG, median (2.5th-97.5th percentiles): 125 (28-399) days motesanib arm; 93 (-6 to 312) days placebo arm.",
      "Because the TS sub-model has no treatment-arm term, a forward simulation with fresh etas reproduces the POOLED distribution only; the per-arm split in the source arises from the individual empirical-Bayes estimates.",
      "",
      "Published OS simulation results, usable as validation targets for the OS sub-model:",
      "  Table 4 (MONET1 replication, updated OS model): predicted hazard ratio 0.88 (95 percent PI 0.78, 0.99) full population, observed 0.90; predicted 0.76 (95 percent PI 0.58, 1.028) Asian patients, observed 0.66.",
      "  Table 5 (virtual 36-month phase III in Asian patients): median OS 594 days motesanib vs. 450 days placebo, HR 0.74; power 72 / 81 / 88 percent for 400 / 500 / 600 patients.",
      sep = "\n"
    )
  )

  ini({
    # ==================================================================
    # TS sub-model: simplified tumor growth inhibition model.
    # Claret 2014 Supplementary Table S1 ("Simplified TGI Model
    # Parameter Estimates"). Columns are Estimate | RSE, % | SD, where
    # the table's abbreviation line defines SD as the "standard
    # deviation of inter-patient variability in diagonal covariance
    # matrix". The reported SD values are therefore log-scale standard
    # deviations and are SQUARED below to give the variances nlmixr2's
    # `eta ~ variance` syntax expects. The matrix is diagonal, so no
    # off-diagonal correlations are encoded.
    #
    # All three rate constants carry units of 1/day (Claret 2014
    # Methods: "The parameters KL, KDE, and lambda have inverse of time
    # units (e.g., day^-1)").
    # ==================================================================

    # Baseline tumor size BASE. The estimated typical value doubles as
    # the per-subject initial condition of the tumor_size state AND as
    # the baseline-TS covariate of the OS sub-model.
    lrbase   <- log(90)      ; label("Typical baseline tumor size BASE (mm)")                                            # Table S1: BASE (mm) = 90, RSE 2.37 %

    # Tumor growth rate KL. Mapped to the register's TGI-family generic
    # growth-rate constant `kge` (exponential growth law
    # d(TS)/dt = kge * TS), as used by tgi_no_sat_expo and
    # Struemper_2025_tumorsize_OS_nsclc.
    lkge      <- log(0.00107) ; label("Tumor growth rate KL (1/day)")                                                     # Table S1: KL (day^-1) = 0.00107, RSE 7.70 %

    # ------------------------------------------------------------------
    # ERRATUM: Supplementary Table S1 prints the KDE and lambda ESTIMATES
    # TRANSPOSED. They are exchanged here. The table as published reads
    # KDE = 0.0193 and lambda = 0.0103; this file uses KDE = 0.0103 and
    # lambda = 0.0193. Nothing is tuned -- the two candidate values are
    # both taken verbatim from Table S1, only their row assignment is
    # corrected, and the correction is forced by three of the paper's own
    # published outputs.
    #
    # The DECISIVE test is eta-free. The TS model is a single POOLED fit
    # with no treatment-arm term, so its TYPICAL value (all etas = 0) must
    # fall BETWEEN the two published per-arm medians (Results, "TS model").
    # It does so for BOTH metrics under the corrected assignment and for
    # NEITHER under the printed one:
    #
    #                        typical TTG    typical TSratio(week 8)
    #   published arm range   93 .. 125       0.69 .. 0.77
    #   as printed in S1        280.8           0.467      <- outside both
    #   corrected (this file)   117.3           0.746      <- inside both
    #
    # Because the typical value contains no etas, no covariance or
    # empirical-Bayes-shrinkage argument can account for the gap: it can
    # only be a wrong point estimate.
    #
    # Two further published outputs agree independently:
    #   (a) forward Monte Carlo over the Table S1 IIV gives median TTG
    #       104 d and median TSratio 0.79 (2.5th pct 0.269) corrected, vs
    #       257 d and 0.51 (2.5th pct 0.05) as printed; published pooled
    #       values are ~108 d and 0.69-0.77 (2.5th pct 0.24-0.27).
    #   (b) feeding the published TTG medians through the Table 3 OS
    #       regression reproduces the Table 5 simulated median OS for
    #       Asian patients (445 d computed vs 450 d published on placebo;
    #       517 d vs 594 d on motesanib). The printed Table S1 typical TTG
    #       of 281 d instead implies 864 d -- 1.8x the 479 d implied by the
    #       published POOLED TTG of ~108 d (the like-for-like comparison,
    #       both being pooled typical values), and still 1.5x Table 5's
    #       better arm.
    #
    # Only the ESTIMATE column is transposed here; RSE and SD stay on
    # their printed rows. The alternative reading -- that the whole rows
    # (label + all three values) are swapped -- is numerically almost
    # indistinguishable (median TTG 105 d vs 104 d) and so is not
    # load-bearing; the estimate-only reading is preferred because it
    # matches the published lower tails more closely (TSratio 2.5th pct
    # 0.269 vs 0.196 against a published 0.24-0.27). See the vignette
    # section "Supplementary Table S1 transposes KDE and lambda".
    # ------------------------------------------------------------------

    # Tumor growth inhibition rate KDE. In this SIMPLIFIED version of the
    # Claret TGI model the drug term is the bare constant KDE rather than
    # the exposure-driven product (kdrug * concentration) of the full
    # model, so the single estimate lumps the motesanib + CP effect
    # across both arms.
    lkdrug  <- log(0.0103)  ; label("Tumor growth inhibition rate KDE (1/day)")                                         # Table S1: estimate 0.0103 read off the lambda row (transposed, see ERRATUM); RSE 5.16 % and SD 0.679 stay on the printed KDE row

    # Rate constant for treatment-effect decay. Accounts for the
    # appearance of resistance and for dose reductions, if any.
    llambda <- log(0.0193)  ; label("Rate constant for treatment-effect decay lambda (1/day)")                          # Table S1: estimate 0.0193 read off the KDE row (transposed, see ERRATUM); RSE 3.72 % and SD 0.797 stay on the printed lambda row (the .doc-to-text conversion mangles the Greek lambda; the abbreviation line confirms the row)

    # Inter-individual variability: log-normal, diagonal covariance
    # matrix. variance = (Table S1 SD)^2. The KDE/lambda erratum above
    # transposes the ESTIMATE column only, so every SD below is read off
    # the row it is printed on, exactly as published.
    etalrbase   ~ 0.683^2    # Table S1: SD for BASE = 0.683  -> variance 0.4665
    etalkge      ~ 1.18^2     # Table S1: SD for KL   = 1.18   -> variance 1.3924
    etalkdrug  ~ 0.679^2    # Table S1: SD for KDE  = 0.679  -> variance 0.4610
    etallambda ~ 0.797^2    # Table S1: SD for lambda = 0.797 -> variance 0.6352

    # Residual error on TS. Claret 2014 Methods: "additive Gaussian
    # distribution was assumed for residual error variability"; the
    # Table S1 abbreviation line defines Sigma as the "standard-deviation
    # of the additive residual error", so 7.849 is already an SD in mm.
    addSd <- 7.849 ; label("Additive residual SD on tumor size (mm)")                                                   # Table S1: Sigma (mm) = 7.849, RSE 1.18 %

    # ==================================================================
    # OS sub-model: parametric accelerated failure time (AFT) regression
    # with a log-normal survival-time distribution, fitted with survreg()
    # in R. Claret 2014 Table 3, "Parameter estimates of the final OS
    # model for evaluable patients (n = 934)". Survival time is analysed
    # in DAYS (Table 3 footnote c).
    #
    # There is no IIV on the OS parameters: the AFT scale parameter IS
    # the model's only source of between-subject survival variability,
    # so a separate eta would be redundant (and is not reported).
    # The parameter covariance matrix is Supplementary Table S3.
    # ==================================================================
    mu_os_int_logdays   <- 1.079    ; label("Intercept of mu_OS on the log(days) AFT scale")                            # Table 3: Intercept = 1.079 (SE 0.436, Z 2.48, P 0.0133)
    e_tsb_mu_os     <- -0.00197 ; label("Slope of baseline tumor size on mu_OS (per mm)")                           # Table 3: Baseline TS (mm) = -0.00197 (SE 0.00035, Z -5.59, P < 0.00001)
    e_ttg_mu_os         <- 0.97695  ; label("Slope of log(TTG + 100) on mu_OS (per log-day)")                           # Table 3: Log(TTG*(day)) = 0.97695 (SE 0.0794, Z 12.3, P < 0.00001)
    e_race_asian_mu_os  <- 0.319    ; label("Additive effect of Asian ethnicity on mu_OS (log-days)")                   # Table 3: Asian ethnicity = 0.319 (SE 0.069, Z 4.62, P < 0.00001)

    # Table 3 footnote d defines this coefficient for "former and current
    # smokers as compared with those who had never smoked", i.e. on an
    # EVER-smoker indicator. model() therefore multiplies it by
    # (1 - SMOKE_NEVER); see covariateData[[SMOKE_NEVER]]$notes.
    e_smoke_ever_mu_os <- -0.264   ; label("Additive effect of ever-smoking (former or current vs never) on mu_OS (log-days)")  # Table 3: Smoking history = -0.264 (SE 0.063, Z -4.18, P 0.00003)

    sigma_os_log        <- -0.24    ; label("log(scale) of the log-normal AFT survival distribution (unitless)")        # Table 3: Log(scale) = -0.24 (SE 0.0282, Z -8.36, P < 0.00001); footnote e: "Scale is the standard deviation of log(OS)"
  })

  model({
    # ---- Individual TS parameters (log-normal IIV) ----
    rbase  <- exp(lrbase  + etalrbase)
    kge    <- exp(lkge    + etalkge)
    kdrug  <- exp(lkdrug  + etalkdrug)
    lambda <- exp(llambda + etallambda)

    # ---- TS sub-model (Claret 2014 Methods, TS model equation) ----
    #   TS(t) = BASE * exp[ KL*t - (KDE/lambda) * (1 - exp(-lambda*t)) ]
    # Differentiating that closed form gives the ODE encoded here:
    #   d(TS)/dt = TS * [ KL - KDE * exp(-lambda*t) ]
    # with TS(0) = BASE. The drug term decays from KDE at t = 0 towards 0
    # as resistance develops, so net growth resumes once
    # KDE*exp(-lambda*t) falls below KL -- which happens exactly at TTG.
    d/dt(tumor_size) <- tumor_size * (kge - kdrug * exp(-lambda * t))
    tumor_size(0)    <- rbase

    # ---- Time to tumor growth (TTG), the model-predicted time to nadir TS ----
    # Setting d(TS)/dt = 0 gives KDE*exp(-lambda*TTG) = KL, hence
    #   TTG = (log(KDE) - log(KL)) / lambda            [Claret 2014 Methods, metric 2]
    # TTG is in DAYS. The Methods sentence defining metric 2 says "the time to
    # nadir TS (in weeks)", but that word is contradicted by the same Methods
    # paragraph ("KL, KDE, and lambda have inverse of time units (e.g.,
    # day^-1)"), by the Table 3 row heading "Log(TTG*(day))", by Table 3
    # footnote c ("Survival time was analyzed in days") and by the Results
    # medians ("125 (28-399) vs. 93 (-6 to 312) days"). Treated as a
    # units-labelling error; see vignette Errata item 11.
    # TTG is NEGATIVE whenever KL > KDE (a tumor that never shrinks), so
    # the OS model uses the shifted log
    #   log(TTG*) = log(TTG + 100)                     [Claret 2014 Methods]
    # NOTE: the published model is undefined for TTG <= -100 days, and no
    # guard is added here because none is published -- the observed MONET1
    # TTG range was -6 to 399 days, comfortably inside the domain. Under
    # forward simulation with freshly drawn etas a small fraction of
    # subjects falls outside it and yields NaN; the vignette quantifies
    # that fraction rather than silently clipping it.
    ttg     <- (log(kdrug) - log(kge)) / lambda
    log_ttg <- log(ttg + 100)

    # ---- OS sub-model: log-normal AFT (Claret 2014 Table 3) ----
    # mu_OS is the mean of log(OS) in log(days); the exponentiated value
    # exp(mu_OS) is the MEDIAN survival time in days.
    # Baseline TS enters as the individual TS-model baseline `rbase`,
    # mirroring the joint-model convention of
    # modellib('Struemper_2025_tumorsize_OS_nsclc').
    smoke_ever <- 1 - SMOKE_NEVER

    mu_os_logdays <- mu_os_int_logdays +
      e_tsb_mu_os      * rbase +
      e_ttg_mu_os        * log_ttg +
      e_race_asian_mu_os * RACE_ASIAN +
      e_smoke_ever_mu_os * smoke_ever

    sigma_os <- exp(sigma_os_log)

    # Log-normal survivor function S(t) = 1 - Phi((log t - mu)/sigma) and
    # the corresponding density / hazard. del_t keeps log(t) finite at
    # t = 0; sur is floored inside the hazard quotient so the ratio stays
    # finite in the far tail.
    del_t     <- 1e-6
    t_pos     <- t + del_t
    z_os      <- (log(t_pos) - mu_os_logdays) / sigma_os
    sur       <- 1 - pnorm(z_os)
    pdf_os    <- exp(-0.5 * z_os * z_os) / (t_pos * sigma_os * sqrt(2 * pi))
    hazard    <- pdf_os / (sur + 1e-30)
    cumhazard <- -log(sur + 1e-30)

    # Median OS in days for the individual subject, exposed as a derived
    # output so typical-value simulations can be read off directly.
    median_os <- exp(mu_os_logdays)

    # ---- Observation: tumor size, additive residual error ----
    tumor_size ~ add(addSd)
  })
}
attr(Claret_2014_motesanib_tumorsize_OS, "message") <-
  "Joint TS/OS exploratory model for first-line nonsquamous NSCLC (Claret 2014 MONET1, n=934 TS-evaluable). TS = simplified Claret TGI with resistance decay (kge, kdrug, lambda, rbase); OS = log-normal AFT driven by log(TTG+100) plus baseline TS, RACE_ASIAN and smoking history. Outputs: tumor_size (observable, mm), ttg, log_ttg, mu_os_logdays, sigma_os, median_os, sur, hazard, cumhazard. ERRATUM: Supplementary Table S1 transposes the KDE and lambda estimates; this model uses the corrected KDE=0.0103, lambda=0.0193 per day, which is what reproduces the paper's own published TTG and TSratio medians. The TS sub-model has NO treatment-arm term, so fresh-eta simulation gives the same TTG distribution in both arms; the paper conditioned its trial simulations on observed per-patient TTG instead. log(TTG+100) is undefined for TTG <= -100 days and is deliberately left unguarded. Companion model for TS-nonevaluable patients: modellib('Claret_2014_motesanib_OS_nonevaluable')."
Claret_2014_motesanib_tumorsize_OS
