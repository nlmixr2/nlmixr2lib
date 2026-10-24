Teng_2018_rrmm_pfs_vgpr_mbma <- function() {
  description <- paste0(
    "MBMA. Drug-independent study-level linear regression of MEDIAN ",
    "PROGRESSION-FREE SURVIVAL (mPFS, months) on the rate of VERY GOOD ",
    "PARTIAL RESPONSE OR BETTER (>=VGPR = CR + VGPR, percent) in ",
    "relapsed/refractory multiple myeloma (RRMM), fitted by ordinary least ",
    "squares to 12 treatment arms of seven randomised phase III trials ",
    "published 2014-2017 (4,924 patients): mPFS = 2.62 + 0.34 * >=VGPR, ",
    "R^2 = 0.75. The paper reports this as the weaker of its two ",
    "response-rate surrogates (the ORR model has R^2 = 0.84) and uses only ",
    "the ORR model for its probability-of-technical-success calculations. ",
    "This is a STUDY-LEVEL model: one record is one trial arm and the ",
    "prediction is that arm's median PFS, never an individual patient's ",
    "survival time. It has no dose, no exposure and no time axis. The ",
    "residual SD is not printed and is recovered by refitting the paper's ",
    "own Table 1 dataset. Companion model: Teng_2018_rrmm_pfs_orr_mbma."
  )

  reference <- paste(
    "Teng Z, Gupta N, Hua Z, Liu G, Samnotra V, Venkatakrishnan K, Labotka R.",
    "Model-Based Meta-Analysis for Multiple Myeloma: A Quantitative",
    "Drug-Independent Framework for Efficient Decisions in Oncology Drug",
    "Development.",
    "Clin Transl Sci. 2018;11(3):218-225.",
    "doi:10.1111/cts.12524.",
    "Regression: Results 'Relationship between response rates and PFS' and",
    "Figure 1b. Data: Table 1.",
    sep = " "
  )

  vignette <- "Teng_2018_rrmm_pfs_mbma"

  units <- list(
    time = "month",
    dosing = "n/a (no dosing; the predictor is the trial arm's response rate)",
    concentration = "n/a (no drug concentration is modelled)",
    response = "month (mpfs, the study arm's median progression-free survival)"
  )

  covariateData <- list(
    RESP_VGPR_PCT = list(
      description = "Study-arm rate of very good partial response or better: percentage (0-100) of the arm with a best response of VGPR or CR.",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "PERCENT scale, not a fraction: the slope is 0.34 months per ",
        "percentage point. Fitted range 15.7-69.9% (Table 1). Five of the ",
        "seven trials assessed response by IMWG uniform criteria and two ",
        "(PANORAMA-1, ELOQUENT-2) by the older EBMT criteria; the paper pools ",
        "them."
      ),
      source_name = ">=VGPR (%) (Teng 2018 Table 1; Results regression and Figure 1b)"
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
      "The same 12 arms as Teng_2018_rrmm_pfs_orr_mbma (Table 1, excluding ",
      "the daratumumab arms of POLLUX and CASTOR, whose median PFS had not ",
      "been reached)."
    ),
    trial_design = "Seven randomised controlled phase III trials identified from the van Beurden-Tan 2017 systematic review, reporting PFS between 2014 and 2017; arm sizes 247-465.",
    endpoint = "Median progression-free survival per arm, 7.2-26.3 months (mean 15.61) over the 12 fitted arms.",
    notes = paste0(
      "Study-level (aggregate) published data only; one record per trial ",
      "arm. Unweighted ordinary least squares in R, with no between-study ",
      "random effect. Refitting Table 1 reproduces the printed coefficients ",
      "(2.6193, 0.3405, R^2 0.7527)."
    )
  )

  ini({
    # Teng 2018 Results 'Relationship between response rates and PFS':
    # 'PFS (months) = 2.62 + 0.34*>=VGPR' (also printed in the Figure 1b
    # panel). Both are OLS point estimates, so neither is fixed.
    int_mpfs    <- 2.62 ; label("Intercept of median PFS at a >=VGPR rate of 0% (months)")                   # Teng 2018 Results and Figure 1b, beta0 = 2.62
    e_vgpr_mpfs <- 0.34 ; label("Change in median PFS per percentage point of >=VGPR rate (months/%)") # Teng 2018 Results and Figure 1b, beta1 = 0.34

    # Not printed: Teng 2018 reports only R^2 = 0.75. Derived as the OLS
    # residual standard error (10 degrees of freedom) of refitting
    # lm(PFS ~ VGPR) to the 12 Table 1 arms with median PFS reached.
    addSd <- 2.9638 ; label("Between-arm residual SD of median PFS (months)") # derived: OLS refit of Teng 2018 Table 1, not printed in the paper
  })

  model({
    # Teng 2018 Methods 'MBMA' and Results: PFS = beta0 + beta1 * >=VGPR,
    # >=VGPR in percent.
    mpfs <- int_mpfs + e_vgpr_mpfs * RESP_VGPR_PCT

    mpfs ~ add(addSd)
  })
}
