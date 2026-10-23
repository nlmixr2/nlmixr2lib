Teng_2018_rrmm_pfs_orr_mbma <- function() {
  description <- paste0(
    "MBMA. Drug-independent study-level linear regression of MEDIAN ",
    "PROGRESSION-FREE SURVIVAL (mPFS, months) on the OVERALL RESPONSE RATE ",
    "(ORR = CR + VGPR + PR, percent) in relapsed/refractory multiple myeloma ",
    "(RRMM), fitted by ordinary least squares to 12 treatment arms of seven ",
    "randomised phase III trials published 2014-2017 (TOURMALINE-MM1, ASPIRE, ",
    "ENDEAVOR, PANORAMA-1, ELOQUENT-2, POLLUX, CASTOR; 4,924 patients): ",
    "mPFS = -23.46 + 0.56 * ORR, R^2 = 0.84. The arms span proteasome ",
    "inhibitors, an HDAC inhibitor, and anti-SLAMF7 and anti-CD38 antibodies ",
    "on lenalidomide- or bortezomib-dexamethasone backbones, which is the ",
    "basis of the paper's claim that the relationship is disease- rather than ",
    "treatment-specific. This is a STUDY-LEVEL model: one record is one trial ",
    "arm and the prediction is that arm's median PFS, never an individual ",
    "patient's survival time and never a drug concentration. It has no dose, ",
    "no exposure and no time axis. The model also carries step 1-2 of the ",
    "paper's probability-of-technical-success (PTS) framework: mpfs_se is the ",
    "standard error of the predicted median PFS that comes from the binomial ",
    "sampling error of an ORR observed in a phase II arm of N_ARM patients, ",
    "propagated through the slope, so P(mPFS > target) = 1 - pnorm(target, ",
    "mpfs, mpfs_se). The source recommends use only for ORR of 50-85%. The ",
    "residual SD is not printed and is recovered by refitting the paper's own ",
    "Table 1 dataset. A companion model on the >=VGPR rate is ",
    "Teng_2018_rrmm_pfs_vgpr_mbma."
  )

  reference <- paste(
    "Teng Z, Gupta N, Hua Z, Liu G, Samnotra V, Venkatakrishnan K, Labotka R.",
    "Model-Based Meta-Analysis for Multiple Myeloma: A Quantitative",
    "Drug-Independent Framework for Efficient Decisions in Oncology Drug",
    "Development.",
    "Clin Transl Sci. 2018;11(3):218-225.",
    "doi:10.1111/cts.12524.",
    "Regression: Results 'Relationship between response rates and PFS' and",
    "Figure 1a. Data: Table 1. PTS framework: Methods 'MBMA' steps 1-3 and",
    "Figure 2.",
    sep = " "
  )

  vignette <- "Teng_2018_rrmm_pfs_mbma"

  # Study-level MBMA with no drug input: nothing is dosed and the output is a
  # survival duration rather than a concentration, so the dosing and
  # concentration strings are declared n/a as in
  # Zierhut_2016_hcc_antiangiogenic_os_mbma.
  units <- list(
    time = "month",
    dosing = "n/a (no dosing; the predictor is the trial arm's response rate)",
    concentration = "n/a (no drug concentration is modelled)",
    response = "month (mpfs, the study arm's median progression-free survival, and its standard error mpfs_se)"
  )

  covariateData <- list(
    RESP_ORR_PCT = list(
      description = "Study-arm overall response rate: percentage (0-100) of the arm with a best response of PR or better (CR + VGPR + PR).",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "PERCENT scale, not a fraction: the slope is 0.56 months per ",
        "percentage point. Fitted range 54.6-87.1% (Table 1); the paper ",
        "recommends predicting only for an observed ORR of 50-85%, and notes ",
        "that the normal approximation behind mpfs_se needs n*p and n*(1-p) ",
        ">= 5. Five of the seven trials assessed response by IMWG uniform ",
        "criteria and two (PANORAMA-1, ELOQUENT-2) by the older EBMT ",
        "criteria; the paper pools them. The analysis-set mean ORR is ",
        "70.275% (the centring constant of the paper's Bayesian ",
        "re-parameterisation, PFS = 15.61 + 0.56 * ORRc)."
      ),
      source_name = "ORR (%) (Teng 2018 Table 1; Results regression and Figure 1a)"
    ),
    N_ARM = list(
      description = "Number of patients in the phase II arm in which RESP_ORR_PCT was observed; sets the binomial sampling error of that ORR.",
      units = "participants",
      type = "count",
      reference_category = NULL,
      notes = paste0(
        "Enters only mpfs_se, never mpfs. Teng 2018 Methods 'MBMA' step 1: ",
        "ORR ~ N(ORR_hat, ORR_hat * (1 - ORR_hat) / n2) with n2 the phase II ",
        "sample size; the worked example (Figure 2) uses n2 = 50. When the ",
        "model is used only to predict mPFS for a known arm, any positive ",
        "value may be supplied and mpfs_se ignored."
      ),
      source_name = "n2 (Teng 2018 Methods 'MBMA' step 1); N = 50 (Figure 2)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 4924L,
    n_studies = 7L,
    n_arms = 12L,
    age_range = "Not reported at the arm level; adult patients (Methods 'Data used in the MBMA framework').",
    weight_range = "Not reported; not a covariate.",
    sex_female_pct = NA_real_,
    disease_state = paste0(
      "Relapsed and/or refractory multiple myeloma after at least one prior ",
      "line of therapy, primarily 1-3 prior lines (Discussion)."
    ),
    dose_range = "Not modelled. Treatment enters only through the arm's response rate.",
    regions = "Multinational phase III trials.",
    treatments = paste0(
      "Table 1 arms: ixazomib-Rd and placebo-Rd (TOURMALINE-MM1), ",
      "carfilzomib-Rd and Rd (ASPIRE), carfilzomib-dexamethasone and Vd ",
      "(ENDEAVOR), panobinostat-Vd and placebo-Vd (PANORAMA-1), elotuzumab-Rd ",
      "and Rd (ELOQUENT-2), Rd (POLLUX control) and Vd (CASTOR control). The ",
      "daratumumab-Rd (POLLUX) and daratumumab-Vd (CASTOR) arms are listed in ",
      "Table 1 but excluded from the fit because their median PFS had not ",
      "been reached. Rd = lenalidomide-dexamethasone, Vd = ",
      "bortezomib-dexamethasone."
    ),
    trial_design = "Seven randomised controlled phase III trials identified from the van Beurden-Tan 2017 systematic review, reporting PFS between 2014 and 2017; arm sizes 247-465.",
    endpoint = "Median progression-free survival per arm, 7.2-26.3 months (mean 15.61) over the 12 fitted arms.",
    notes = paste0(
      "Study-level (aggregate) published data only; one record per trial ",
      "arm. Unweighted ordinary least squares in R ('Linear regression in ",
      "R'), not a mixed-effects model: there is no between-study random ",
      "effect and arms are not weighted by sample size. Refitting Table 1 ",
      "reproduces every printed coefficient, including the ",
      "leave-TOURMALINE-MM1-out fit -23.2609 + 0.5547 * ORR to all four ",
      "printed decimals."
    )
  )

  ini({
    # Teng 2018 Results 'Relationship between response rates and PFS':
    # 'PFS (months) = -23.46 + 0.56*ORR' (also printed in the Figure 1a panel).
    # Both are OLS point estimates, so neither is fixed.
    int_mpfs   <- -23.46 ; label("Intercept of median PFS at an ORR of 0% (months; an extrapolation, see the 50-85% supported range)") # Teng 2018 Results and Figure 1a, beta0 = -23.46
    e_orr_mpfs <- 0.56   ; label("Change in median PFS per percentage point of overall response rate (months/%)")                         # Teng 2018 Results and Figure 1a, beta1 = 0.56

    # Not printed: Teng 2018 reports only R^2 = 0.84. Derived as the OLS
    # residual standard error (10 degrees of freedom) of refitting
    # lm(PFS ~ ORR) to the 12 Table 1 arms with median PFS reached; the same
    # refit returns -23.4612, 0.5559 and R^2 0.8449, matching the printed
    # values, so the dataset and estimator are the paper's own.
    addSd <- 2.3472 ; label("Between-arm residual SD of median PFS (months)") # derived: OLS refit of Teng 2018 Table 1, not printed in the paper
  })

  model({
    # Teng 2018 Methods 'MBMA': PFS = beta0 + beta1 * ORR, ORR in percent.
    mpfs <- int_mpfs + e_orr_mpfs * RESP_ORR_PCT

    # Teng 2018 Methods 'MBMA' steps 1-2: the observed phase II ORR is taken
    # as the truth with a normal-approximation binomial error,
    # ORR ~ N(ORR_hat, ORR_hat * (1 - ORR_hat) / n2), and pushed through the
    # regression. The variance is on the proportion scale, hence the factor
    # of 100 to return to percentage points. Only this sampling error is
    # propagated; the regression's own parameter and residual uncertainty
    # are not part of this layer of the paper's framework.
    orrFrac <- RESP_ORR_PCT / 100
    mpfs_se <- e_orr_mpfs * 100 * sqrt(orrFrac * (1 - orrFrac) / N_ARM)

    mpfs ~ add(addSd)
  })
}
