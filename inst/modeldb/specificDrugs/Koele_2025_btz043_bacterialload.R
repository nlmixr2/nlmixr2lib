Koele_2025_btz043_bacterialload <- function() {
  description <- "Joint bilinear exposure-response model for the decline in sputum mycobacterial load during 14 days of BTZ-043 monotherapy in adults with drug-susceptible pulmonary tuberculosis, fitted simultaneously to log10 colony-forming units on solid medium and log10 time to positivity in liquid MGIT culture. Both endpoints share one breakpoint (node) fixed at 48 h; the first-phase slope of each endpoint is an Emax function of the individual model-predicted BTZ-043total (BTZ-043 plus metabolite M2) AUC(0-24), with a single shared EC50, while the second-phase slope (Days 2-14) is exposure-independent. Baseline loads and slopes carry correlated between-subject variability and the residual error is additive on the log10 scale. This is a $PRED regression on time and exposure with no differential equations and no dosing events; the driving exposure is supplied per subject as a covariate column."
  reference <- paste(
    "Koele S. E., Heinrich N., De Jager V. R., Dreisbach J., Phillips P. P. J.,",
    "Gross-Demel P., Dawson R., Narunsky K., Wildner L. M., Mchugh T. D.,",
    "Te Brake L. H. M., Diacon A. H., Aarnoutse R. E., Hoelscher M.,",
    "Svensson E. M. (2025).",
    "Population pharmacokinetics and exposure-response relationship of the",
    "antituberculosis drug BTZ-043.",
    "Journal of Antimicrobial Chemotherapy 80(5):1319-1327.",
    "doi:10.1093/jac/dkaf076.",
    "Structural equations and random-effect variances transcribed from the",
    "final NONMEM control stream in the Supplementary data",
    "('Pharmacodynamic model code'); typical values from Table 3.",
    "The driving exposure is produced by the companion PK model",
    "modellib('Koele_2025_btz043').",
    sep = " "
  )
  vignette <- "Koele_2025_btz043"
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")

  # No d/dt() states: log_cfu and log_ttp are algebraic predictions. log_cfu
  # is the registered canonical log-CFU PD output; log_ttp is declared
  # paper-specific because time to positivity has not previously appeared as
  # a modelled endpoint in this library, and a compartment / output name is
  # promoted to canonical only once a second paper uses it.
  paper_specific_compartments <- c("log_ttp")

  covariateData <- list(
    AUC_BTZ043TOT = list(
      description        = "Individual model-predicted BTZ-043total (BTZ-043 plus metabolite M2, expressed as BTZ-043 equivalents) area under the plasma concentration-time curve over the 24 h dosing interval, used as the driver of the first-phase bacterial-load decline.",
      units              = "ng/mL*h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "This is the EXMET data item of the supplementary pharmacodynamic control stream. Koele 2025 Methods 'PD model development': 'The exposure-response relationship of BTZ-043 was investigated using individual model-derived estimates of the AUC0-24 and Cmax during the intensive PK sampling on Day 12 for Stage 1 and Day 14 for Stage 2', and Results: 'BTZ-043total AUC0-24 and the M2 AUC0-24 were identified as the most significant drivers of the decrease in bacterial load during the first 2 days of treatment (dOFV = -12.1 and -13.1, respectively) ... As the BTZ-043total AUC0-24 is easier to quantify in future studies, the exposure-response model using this parameter as driver was selected.' Time-fixed per subject: one steady-state AUC per participant, carried on every bacterial-load record. UNITS ARE LOAD-BEARING and are mass-based even though the companion PK model works in molar units: the estimated EC50 is 16 900 ng/mL*h. To generate this column from modellib('Koele_2025_btz043'), integrate Cc + Cc_m2 (nmol/L) over the 24 h interval and multiply by the BTZ-043 molecular weight in ng per nmol, 431.39/1000 = 0.43139; the bioanalytical method converts all M2 back to BTZ-043 before quantification, so BTZ-043total is the molar sum reported as BTZ-043 mass. Observed range is not tabulated in Koele 2025; the estimated EC50 is stated to be 'similar to the exposure obtained after a 500 mg BTZ-043 dose following a standard breakfast'.",
      source_name        = "EXMET"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 68L,
    n_studies      = 1L,
    age_range      = "18-57 years",
    age_median     = "27 years",
    weight_range   = "42-81 kg",
    weight_median  = "54 kg",
    sex_female_pct = 16.2,
    race_ethnicity = c(
      Black            = 64.7,
      `Cape-coloured`  = 33.8,
      White            = 1.5
    ),
    hiv_status     = "HIV-1 negative, 68/68 (100%)",
    disease_state  = "Adults aged 18-64 years with drug-susceptible pulmonary tuberculosis enrolled in the sequential Phase 1b/2a trial NCT04044001, receiving BTZ-043 monotherapy for 14 days.",
    dose_range     = "Oral BTZ-043 250-1750 mg once daily for 14 days.",
    sampling       = "Overnight sputum collected before treatment and on Days 2, 3, 4, 6, 8, 11 and 14; each sample cultured for colony-forming units on solid medium and for time to positivity in the MGIT liquid system, in duplicate.",
    regions        = "South Africa (TASK, Cape Town; University of Cape Town Lung Institute)",
    notes          = "Koele 2025 Results 'PD model': 921 cfu and 1113 TTP measurements were available, of which 7 cfu and 2 TTP observations were culture negative (cfu < 1.0 cfu/mL or TTP > 25 days). The proportion of culture-negative observations at baseline (defined as the first 2 days on treatment) was 0.9% for cfu and 0.0% for TTP. No effect of age, weight or study site was identified in the stepwise covariate search (forward inclusion P < 0.05, backward elimination P < 0.01). Demographics are those of the combined Stage 1 + 2 PK/PD analysis population, Koele 2025 Table 1."
  )

  ini({
    # -----------------------------------------------------------------------
    # Typical values from Koele 2025 Table 3 'Estimated PD model parameters'.
    # The supplementary $PRED block confirms every structural relationship;
    # where its $THETA initials are rounded relative to the published finals
    # (0.0039 vs 0.00386 for the first TTP slope, 0.0004 vs 0.000440 for the
    # second), Table 3 is used.
    #
    # Sign convention follows the paper: both slopes are POSITIVE numbers.
    # cfu DECLINES by beta * t and TTP INCREASES by beta * t, which is why
    # Table 3 labels the TTP slopes '-log10 h/h'.
    # -----------------------------------------------------------------------
    le0_cfu     <- log(6.20)      ; label("Baseline sputum bacterial load (log10 cfu/mL)")                              # Koele 2025 Table 3 'Baseline cfu (log10 cfu/mL) 6.20 (5.94-6.42)'; control stream THETA(1)
    le0_ttp     <- log(1.99)      ; label("Baseline time to positivity (log10 h)")                                      # Koele 2025 Table 3 'Baseline TTP (log10 h) 1.99 (1.97-2.02)'; control stream THETA(2)
    lemax_cfu   <- log(0.0270)    ; label("Maximum first-phase rate of cfu decline (log10 cfu/mL per h)")               # Koele 2025 Table 3 'Emax slope 1 cfu (log10 cfu/mL/h) 0.0270 (0.0172-0.0433)'; control stream THETA(3)
    lemax_ttp   <- log(0.00386)   ; label("Maximum first-phase rate of TTP increase (log10 h per h)")                   # Koele 2025 Table 3 'Emax slope 1 TTP (-log10 h/h) 0.00386 (0.00249-0.00626)'; control stream THETA(5)
    lec50       <- log(16900)     ; label("BTZ-043total AUC(0-24) producing half the maximum first-phase slope (ng/mL*h)") # Koele 2025 Table 3 'EC50 BTZ-043 total exposure (ng/mL*h) 16900 (5510-44300)'; control stream THETA(8); shared by both endpoints
    lslope2_cfu <- log(0.00254)   ; label("Second-phase (Days 2-14) rate of cfu decline (log10 cfu/mL per h)")          # Koele 2025 Table 3 'Slope 2 cfu (log10 cfu/mL/h) 0.00254 (0.00191-0.00311)'; control stream THETA(4)
    lslope2_ttp <- log(0.000440)  ; label("Second-phase (Days 2-14) rate of TTP increase (log10 h per h)")              # Koele 2025 Table 3 'Slope 2 TTP (-log10 h/h) 0.000440 (0.000351-0.000510)'; control stream THETA(6)
    lnode_time  <- fixed(log(48)) ; label("Node: time after the start of treatment at which the bilinear slope changes (h)") # Koele 2025 Table 3 'Node (h) 48 FIX'; control stream THETA(7) FIX. Results: 'Before the exposure-response analysis, the node parameter was fixed to the best estimate (48 h) to stabilize the estimation properties of the model.'

    # -----------------------------------------------------------------------
    # Between-subject variability. Koele 2025 Table 3 reports CV% with the
    # same footnote transform as the PK table, CV% = sqrt(e^OM2 - 1); the
    # variances below are the $OMEGA entries of the supplementary control
    # stream and were each checked against the printed CV% and correlations.
    #
    # BLOCK(3) on (baseline cfu, baseline TTP, first-phase TTP slope):
    #   0.0229                     -> CV 15.2%   (Table 3: 15.2)
    #  -0.0053  0.0026             -> CV  5.10%  (Table 3: 5.14)
    #  -0.0211 -0.0066  0.217      -> CV 49.2%
    # correlations -68.7%, -29.9%, -27.8% (Table 3: -68.3, -29.9, -27.8).
    #
    # NOTE ON THE TABLE 3 TYPO: Table 3 prints '75.9 (35.4-67.1)' for the
    # variability of the first TTP slope, a point estimate that lies OUTSIDE
    # its own confidence interval and that duplicates the value printed one
    # row below for the second cfu slope. The control stream's variance of
    # 0.217 gives CV = 49.2%, which sits inside the printed 35.4-67.1
    # interval, so 0.217 is used and 75.9 is treated as a copy-paste error.
    # See the vignette Errata section.
    #
    # BLOCK(2) on (second-phase cfu slope, second-phase TTP slope):
    #   0.455                      -> CV 75.9%   (Table 3: 75.9)
    #   0.247  0.221               -> CV 49.7%   (Table 3: 49.7), corr +77.9%
    # -----------------------------------------------------------------------
    etale0_cfu + etale0_ttp + etalemax_ttp ~ c(0.0229,
                                              -0.0053,  0.0026,
                                              -0.0211, -0.0066, 0.217)  # control stream $OMEGA BLOCK(3)
    etalslope2_cfu + etalslope2_ttp ~ c(0.455,
                                        0.247, 0.221)                    # control stream $OMEGA BLOCK(2)

    # -----------------------------------------------------------------------
    # Residual error, additive on the log10 scale. Koele 2025 Table 3
    # footnote a: 'The average additive errors and correlations between cfu
    # and TTP replicates is reported here. Parameter estimates for the
    # replicates are presented in the model code shown in the Supplementary
    # material.' The supplementary $SIGMA BLOCK(4) resolves the four
    # replicate-level SDs and their correlations; nlmixr2lib has no idiomatic
    # encoding for a replicate-level residual or for a residual correlation
    # ACROSS outputs, so the published averages are used here and the full
    # 4x4 matrix is reproduced in the vignette Errata.
    # -----------------------------------------------------------------------
    addSd_log_cfu <- 0.558  ; label("Additive residual error on cfu (log10 cfu/mL)") # Koele 2025 Table 3 'Additive error cfu (log10 cfu/mL) 0.558 (0.526-0.592)'; mean of sqrt(0.299) = 0.547 and sqrt(0.322) = 0.567 from the control stream $SIGMA BLOCK(4)
    addSd_log_ttp <- 0.0624 ; label("Additive residual error on TTP (log10 h)")      # Koele 2025 Table 3 'Additive error TTP (log10 h) 0.0624 (0.0597-0.0665)'; mean of sqrt(0.0043) = 0.0656 and sqrt(0.0035) = 0.0592 from the control stream $SIGMA BLOCK(4)
  })

  model({
    # 1. Baselines at the start of treatment. Exponential between-subject
    #    variability, matching the control stream's THETA * EXP(ETA).
    base_cfu <- exp(le0_cfu + etale0_cfu)
    base_ttp <- exp(le0_ttp + etale0_ttp)

    # 2. First-phase (0 to 48 h) slopes: an Emax function of the individual
    #    BTZ-043total AUC(0-24) with one EC50 shared by both endpoints
    #    (control stream BETA1BTZCFU = EXMET*THETA(3)/(THETA(8)+EXMET) and
    #    BETA1BTZTTP = EXMET*THETA(5)/(THETA(8)+EXMET)).
    #
    #    The between-subject eta on the TTP first-phase slope is applied here
    #    even though the supplementary $PRED omits it: it is declared as the
    #    third element of $OMEGA BLOCK(3) with a non-zero variance and two
    #    estimated correlations, Table 3 reports its variability, and the
    #    Results state 'IIV was identified on the baseline bacterial load,
    #    the first slope for TTP and the second slopes for cfu and TTP. The
    #    model was unable to determine IIV on the first cfu slope.' The
    #    omission from the printed $PRED is therefore a transcription slip,
    #    not a modelling choice; see the vignette Errata section. There is
    #    correspondingly NO eta on the cfu first-phase slope.
    emax_cfu <- exp(lemax_cfu)
    emax_ttp <- exp(lemax_ttp + etalemax_ttp)
    ec50     <- exp(lec50)
    beta1_cfu <- emax_cfu * AUC_BTZ043TOT / (ec50 + AUC_BTZ043TOT)
    beta1_ttp <- emax_ttp * AUC_BTZ043TOT / (ec50 + AUC_BTZ043TOT)

    # 3. Second-phase (48 h to end of treatment) slopes: exposure-independent.
    #    Koele 2025 Results: 'No exposure-response relationship could be
    #    identified over Days 2-14 after the start of treatment.'
    beta2_cfu <- exp(lslope2_cfu + etalslope2_cfu)
    beta2_ttp <- exp(lslope2_ttp + etalslope2_ttp)

    # 4. Bilinear time split about the node. tbefore is the time spent in the
    #    first phase and tafter the time spent in the second, so the two
    #    branches of the control stream's IF (TIME.GT.NODE) are written
    #    branch-free and remain continuous at the node.
    nodetime <- exp(lnode_time)
    tafter   <- (time > nodetime) * (time - nodetime)
    tbefore  <- time - tafter

    # 5. Predictions. time is time after the start of treatment, in hours.
    log_cfu <- base_cfu - beta1_cfu * tbefore - beta2_cfu * tafter
    log_ttp <- base_ttp + beta1_ttp * tbefore + beta2_ttp * tafter

    log_cfu ~ add(addSd_log_cfu)
    log_ttp ~ add(addSd_log_ttp)
  })
}
