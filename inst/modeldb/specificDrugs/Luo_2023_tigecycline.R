Luo_2023_tigecycline <- function() {
  description <- "One-compartment linear IV population PK model for tigecycline in critically ill adult ICU patients. Clearance decreases linearly with APACHE II score above the cohort median of 22.5 points, and central volume decreases linearly with age above the cohort median of 72 years. Fitted by NONMEM 7.3.0 FOCE-I to 143 steady-state plasma concentrations from 54 ICU patients on the licensed 100 mg loading / 50 mg q12h maintenance regimen."
  reference   <- "Luo X, Wang S, Li D, Wen J, Sun N, Fan G. Population pharmacokinetics of tigecycline in critically ill patients. Front Pharmacol. 2023;14:1083464. doi:10.3389/fphar.2023.1083464"
  vignette    <- "Luo_2023_tigecycline"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "tigecycline", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    APACHE_II = list(
      description        = "Acute Physiology and Chronic Health Evaluation II score at ICU admission",
      units              = "points",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Cohort median 22.50 points (IQR 16.50-27.00; Luo 2023 Table 1).",
        "Reference / centring value = 22.50 points, the cohort median. The paper prints the",
        "clearance equation without a centring term (Results eq. 1: CL = (11.30 - 0.14 x APACHE II) x e^0.065),",
        "but that uncentred reading is falsified by three of the paper's own outputs -- see the",
        "ini() comment on e_apache_cl and the vignette 'Assumptions and deviations' section.",
        "The centred additive-linear form used here matches the published ICU precedent",
        "Swart_2004_midazolam.R (Q = 40.8 - (APACHE - 26) x 2.75)."
      ),
      source_name        = "APACHE II"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Cohort median 72.0 years (IQR 57.5-80.3; Luo 2023 Table 1);",
        "inclusion required age >= 18 years. Reference / centring value = 72.0 years, the cohort median.",
        "As for APACHE_II, the paper prints the volume equation uncentred (Results eq. 2:",
        "V = [105.00 x (1 - 0.0059 x AGE)] x e^0.160); the centred reading is adopted here.",
        "See the ini() comment on e_age_vc."
      ),
      source_name        = "Age"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as a candidate covariate (Methods 2.4) but not retained in the final model. Cohort median 68.0 kg (IQR 58.3-70.0; Table 1). The Discussion attributes the null weight effect to the small sample size and the narrow, elderly age distribution."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a candidate covariate (Methods 2.4) but not retained in the final model. Cohort median 33.37 U/L (IQR 21.80-77.15; Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened as a candidate covariate (Methods 2.4) but not retained in the final model. Cohort median 27.91 g/L (IQR 25.58-32.50; Table 1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 54L,
    n_studies      = 1L,
    age_range      = "Adults >= 18 years; median 72.0 (IQR 57.5-80.3)",
    age_median     = "72.0 years",
    weight_range   = "Median 68.0 kg (IQR 58.3-70.0)",
    weight_median  = "68.0 kg",
    sex_female_pct = 44.4,
    race_ethnicity = "Not reported (single-centre Chinese ICU cohort, Dalian)",
    disease_state  = paste(
      "Critically ill adult ICU patients with clinically confirmed or suspected Gram-positive",
      "and/or Gram-negative infection, mostly pulmonary. APACHE II median 22.50 (IQR 16.50-27.00).",
      "Baseline laboratory medians (Table 1): ALT 29.09 U/L, AST 33.37 U/L, ALP 89.30 U/L,",
      "total bilirubin 16.22 mmol/L, serum creatinine 77.49 umol/L, BUN 9.30 mmol/L,",
      "albumin 27.91 g/L, sodium 137.0 mmol/L. 28 of 54 received concomitant antifungal therapy."
    ),
    dose_range     = "Licensed regimen only: 100 mg IV loading dose then 50 mg IV q12h maintenance, for at least 3 days",
    regions        = "China (Second Affiliated Hospital of Dalian Medical University)",
    notes          = paste(
      "Retrospective single-centre study, samples collected December 2017 - July 2018.",
      "143 plasma concentrations from 54 patients (median 2.6 samples per patient), drawn at",
      "steady state pre-dose and 1, 2 and 4 h after administration. Median observed concentration",
      "444.0 ng/mL (IQR 222.3-716.6). Assay: LC-MS/MS, linear 1-2000 ng/mL (R^2 = 0.9907),",
      "LOD 10 ng/mL, intra- and inter-day RSD 4.15% and 2.74%. Estimation by NONMEM 7.3.0 FOCE-I;",
      "final model verified by 500 bootstrap replicates (442 successful, 88.4%), VPC and NPDE.",
      "Baseline demographics are Luo 2023 Table 1; parameter estimates are Luo 2023 Table 2."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Luo 2023 Table 2 reports four thetas, two omegas and one sigma. The
    # Results text (section 3.2) prints the two parameter equations as
    #
    #   CL (L/h) = (11.30 - 0.14 x APACHE II score) x e^0.065        (eq. 1)
    #   V  (L)   = [105.00 x (1 - 0.0059 x AGE)]   x e^0.160         (eq. 2)
    #
    # The trailing e^0.065 / e^0.160 are typesetting of exp(eta) with the
    # OMEGA printed in place of the eta symbol; 0.065 and 0.160 are Table 2
    # Omega1 and Omega2 (see the etalcl / etalvc lines below), not part of
    # the fixed-effect expression.
    # ------------------------------------------------------------------

    # Typical CL at the cohort-median APACHE II of 22.5 points.
    lcl <- log(11.30);  label("Clearance at APACHE_II = 22.5 points (L/h)")   # Luo 2023 Table 2 Theta1 = 11.30 (bootstrap median 10.90, 95% CI 8.75-13.40)

    # Typical central volume at the cohort-median age of 72 years.
    lvc <- log(105.00); label("Central volume at AGE = 72 years (L)")         # Luo 2023 Table 2 Theta2 = 105.00 (bootstrap median 100.30, 95% CI 54.7-160.0)

    # ------------------------------------------------------------------
    # Covariate slopes. Both are transcribed at their printed magnitudes;
    # what is NOT taken from the printed equations is the centring, which
    # eqs. 1-2 omit. The uncentred reading is falsified by three of the
    # paper's own outputs (all reproduced in the vignette):
    #
    #  1. Table 4 lists AUC0-24 = 8.85 / 13.27 / 17.70 mg.h/L for 24 h doses
    #     of 100 / 150 / 200 mg. Methods 2.6 defines AUC0-24 = dose24 / CL,
    #     so the paper's own typical CL is 100/8.85 = 11.30 L/h -- the value
    #     of eq. 1 at APACHE II = 22.5 under centring, not the 8.15 L/h that
    #     eq. 1 gives at the median APACHE II when read uncentred.
    #  2. Re-running the Monte Carlo PTA grid of Table 4 (72 cells) with a
    #     typical CL of 11.30 L/h reproduces it to an RMSE of 0.9 percentage
    #     points; with the uncentred 8.15 L/h the RMSE is 19.4 points, and
    #     cells reported as 40.96% come out at 89%.
    #  3. Simulating the study's own sampling design (steady state pre-dose and
    #     1, 2, 4 h on 50 mg q12h) with this model's IIV and RUV gives a median
    #     concentration 3.6% above the observed median of 444.0 ng/mL under
    #     centring, versus 56.3% above it when read uncentred (both measured in
    #     the validation vignette). A converged FOCE-I fit whose reported
    #     goodness-of-fit is symmetric about the line of identity cannot
    #     overpredict the median of its own data by half.
    #
    # The centred additive-linear form also matches the registered ICU
    # precedent Swart_2004_midazolam.R, whose published APACHE II effect is
    # written with the cohort mean subtracted: Q = 40.8 - (APACHE - 26) * 2.75.
    # ------------------------------------------------------------------
    e_apache_cl <- 0.14;   label("Additive linear slope of (APACHE_II - 22.5) on CL (L/h per point, subtractive)")  # Luo 2023 Table 2 Theta3 = 0.14 (bootstrap median 0.13, 95% CI 0.03-0.21); coefficient of eq. 1
    e_age_vc    <- 0.0059; label("Fractional linear slope of (AGE - 72) on Vc (per year, subtractive)")             # Luo 2023 Table 2 Theta4 = 0.0059 (bootstrap median 0.0049, 95% CI 0.0033-0.0079); coefficient of eq. 2

    # ------------------------------------------------------------------
    # Inter-individual variability. Results 3.2: "The index model was the
    # best for both the interindividual variation and residual variation
    # models" -- i.e. exponential IIV and exponential RUV. Table 2 Omega1 /
    # Omega2 are NONMEM OMEGA variances on the log scale; the corresponding
    # CVs are sqrt(exp(0.065) - 1) = 25.9% for CL and sqrt(exp(0.160) - 1)
    # = 41.7% for V.
    # ------------------------------------------------------------------
    etalcl ~ 0.065  # Luo 2023 Table 2 Omega1 = 0.065 (bootstrap median 0.060, 95% CI 0.021-0.106)
    # Table 2's printed 95% CI for Omega2, '0.018 - 0.0483', excludes both the
    # point estimate 0.160 and the bootstrap median 0.153; it is a typographical
    # error in the source. Recorded here, not used.
    etalvc ~ 0.160  # Luo 2023 Table 2 Omega2 = 0.160 (bootstrap median 0.153)

    # Residual error. Exponential ("index") residual model, so the Table 2
    # Sigma1 variance maps to lnorm() with SD sqrt(0.0316) = 0.17776 on the
    # log scale (17.8% CV). Read as an SD instead of a variance it would be
    # 3.16%, below the assay's own 4.15% intra-day RSD, which rules that
    # reading out.
    expSd <- 0.17776; label("Log-normal residual SD on Cc (log-scale SD)")  # Luo 2023 Table 2 Sigma1 = 0.0316 variance; sqrt(0.0316) = 0.17776
  })

  model({
    # Individual PK parameters. Both covariate effects are linear deviations
    # from the typical value at the cohort median: additive (L/h per APACHE II
    # point) on clearance, fractional (per year) on central volume. Exponential
    # IIV on each.
    cl  <- (exp(lcl) - (APACHE_II - 22.5) * e_apache_cl) * exp(etalcl)
    vc  <- exp(lvc + etalvc) * (1 - e_age_vc * (AGE - 72))

    kel <- cl / vc

    # One-compartment linear disposition with first-order elimination. Doses
    # are short IV infusions administered into the central compartment through
    # the event table.
    d/dt(central) <- -kel * central

    # Tigecycline plasma concentration (mg/L; doses in mg, volumes in L).
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
