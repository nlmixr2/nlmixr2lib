Yang_2010_rosuvastatin_mbma <- function() {
  description <- "MBMA. Literature-based meta-analysis simple Emax dose-response model for percentage reduction in low-density lipoprotein cholesterol (LDL-C) from baseline in adult hypercholesterolemia patients receiving rosuvastatin. Operates at the study-arm level over 14 dose-ranging trials (46 study-arm-mean effect samples; 9 Western trials and 5 Asian trials, total N substantially larger than 46 because each arm pools many patients). Output Cc is the study-arm mean percent LDL-C reduction from baseline (unsigned: Cc = 50 means a 50 percent reduction). The placebo intercept E0 (-0.802 percent, a small expected LDL-C increase under placebo) and the Hill / sigmoidicity exponent (1) are fixed at the values used by the source paper -- E0 from prior literature [Mandema 2005, ref 15] and gamma after the sigmoidal Emax fit produced unstable estimates. Race (Asian vs Western reference) is the only retained covariate and acts on ED50: ED50_Asian = ED50_Western * 0.564 (i.e. roughly twofold-lower ED50 in Asians). Between-trial variability is encoded as a single study-arm-level eta on the predicted output (SD 3.0 percent); residual error is additive (SD 3.1 percent). Baseline LDL-C was screened but not retained. Suitable simulation scope is study-arm-mean percent LDL-C reduction, NOT individual-subject LDL-C trajectories. The model also predicts only the steady-state effect (paper restricted to arms with at least 4 weeks of treatment)."

  reference <- paste(
    "Yang J, Li LJ, Wang K, He YC, Sheng YC, Xu L, Huang XH, Guo F, Zheng QS.",
    "Race differences: modeling the pharmacodynamics of rosuvastatin",
    "in Western and Asian hypercholesterolemia patients.",
    "Acta Pharmacologica Sinica. 2011;32(1):116-125",
    "(published online 13 Dec 2010).",
    "doi:10.1038/aps.2010.169.",
    "Placebo intercept E0 = -0.802 percent fixed from Mandema JW, Hermann D, Wang W et al.",
    "Model-based development of gemcabene, a new lipid-altering agent.",
    "AAPS J. 2005;7(3):E513-E522 (ref [15]).",
    sep = " "
  )
  vignette <- "Yang_2010_rosuvastatin_mbma"
  units <- list(
    time          = "week (placeholder; the model is a steady-state dose-response and time-independent -- the paper restricts to arms with at least 4 weeks of treatment)",
    dosing        = "mg/day (per-arm daily rosuvastatin dose supplied as the DOSE covariate column; the model is an MBMA dose-response and does not consume rxode2 dose events)",
    concentration = "%/arm (study-arm-mean unsigned percent LDL-C reduction from baseline; e.g. Cc = 50 means a 50 percent reduction. Output Cc is NOT a drug concentration; the slash in the unit string is to satisfy checkModelConventions parsing)"
  )

  covariateData <- list(
    DOSE = list(
      description        = "Per-arm daily rosuvastatin dose (mg/day; 0 for placebo arms).",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA study-arm-level covariate (study-arm mean daily dose, time-invariant within the dose-ranging arm). Dose-ranging arms used 0 (placebo), 1, 2, 2.5, 4, 5, 10, 20, 40, and 80 mg/day across the 14 trials (Yang 2010 Tables 1-2). Source paper uses 'Dose' as the independent variable in Eq 1.",
      source_name        = "Dose (Yang 2010 Eq 1, Tables 1-2)"
    ),
    RACE_ASIAN = list(
      description        = "Indicator that the study-arm population is Asian (1) versus Western (0).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = Western (predominantly White / Caucasian) -- the reference category in the source paper (RACE=0 in Yang 2010 covariate-modeling equation).",
      notes              = "MBMA study-arm-level indicator (a property of the trial cohort, not of an individual subject). Yang 2010 defines RACE=0 for Westerners (predominantly Whites / Caucasians) and RACE=1 for Asians (Chinese, Japanese, and South Asian subjects). The covariate enters as a multiplicative factor on ED50: ED50_arm = ED50_Western * e_asian_ed50^RACE_ASIAN, so RACE_ASIAN=1 multiplies ED50 by 0.564 (twofold-lower ED50, i.e. greater LDL-C sensitivity per mg of rosuvastatin) and RACE_ASIAN=0 leaves ED50 at the Western reference. This mirrors the canonical RACE_ASIAN entry in inst/references/covariate-columns.md; the MBMA application aggregates per arm rather than per subject.",
      source_name        = "RACE (Yang 2010 Eq for P = TVP * theta_RACE^RACE)"
    )
  )

  covariatesDataExcluded <- list(
    LDLC = list(
      description        = "Per-arm baseline LDL-C concentration. Tested in the Yang 2010 forward-inclusion covariate screen but NOT retained in the final model.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yang 2010 Methods describes screening baseline LDL-C as a continuous covariate alongside race. Only race-on-ED50 produced a significant OFV reduction (delta OFV = 7.095); baseline LDL-C did not reach significance and was dropped before the backward-elimination step. Recorded here for source-trace completeness; no effect parameter is included in ini() because the paper reports no point estimate for the dropped covariate.",
      source_name        = "Baseline LDL-C (Yang 2010 Methods / Tables 1-2)"
    )
  )

  population <- list(
    species         = "human",
    n_studies       = 14L,
    n_data_points   = 46L,
    age_range       = "adult hypercholesterolemia patients (per-trial age distributions not tabulated in the main paper; Yang 2010 Tables 1-2 list reference/year/duration/dose/N/baseline LDL-C and percent reduction per arm)",
    disease_state   = "adults with primary hypercholesterolemia (heterogeneous: includes some trials in patients with concomitant cardiovascular disease, switched-from-other-statin arms excluded). All arms restricted to at least 4 weeks of treatment to capture the steady-state LDL-C-lowering effect.",
    dose_range      = "rosuvastatin 0-80 mg/day (placebo, 1, 2, 2.5, 4, 5, 10, 20, 40, and 80 mg/day across the 14 dose-ranging trials per Yang 2010 Tables 1-2)",
    race_ethnicity  = "9 Western dose-ranging trials (predominantly Whites / Caucasians) and 5 Asian dose-ranging trials (Chinese, Japanese, and South Asian subjects, including 2 unpublished trials from the authors department). 22 additional one-dose trials (18 Western, 4 Asian) used for external visual-predictive-check validation rather than parameter estimation.",
    baseline_ldlc   = "Per-arm baseline LDL-C ranged about 153-219 mg/dL across the 14 dose-ranging trials (Yang 2010 Tables 1-2). Baseline LDL-C was screened as a covariate but not retained in the final model.",
    regions         = "International (Western: predominantly Whites / Caucasians; Asian: Chinese, Japanese, South Asian). Some trials multicenter, the majority double-blind, three placebo-controlled.",
    notes           = "MBMA at the study-arm level: each modeled data point is the mean percent LDL-C reduction in a group of patients at the steady-state timepoint in a single dose-ranging trial arm. The model is intended for simulating study-arm-mean percent LDL-C reductions and is NOT suitable for individual-subject simulation. Inter-trial variability eta on the predicted output is reported as a standard deviation (3.0 percent) -- variance is 9.0. Residual error is reported as a standard deviation (3.1 percent) on the same unsigned percent-reduction scale. The final model was developed against the 14 dose-ranging trials only; the 22 one-dose trials served as an independent visual-predictive-check holdout (Yang 2010 Figure 6)."
  )

  ini({
    # ============================================================
    # Simple Emax dose-response model (Yang 2010 Eq 1, simplified
    # after gamma fixed to 1):
    #   Y_obs_arm = E0 + Emax * DOSE^gamma / (ED50_arm^gamma + DOSE^gamma)
    #             + eta_study + epsilon
    # with ED50_arm = ED50_Western * theta_RACE^RACE_ASIAN.
    #
    # Parameter values are the Yang 2010 Table 4 "TVP" column.
    # Sign convention: Y is the unsigned percent LDL-C reduction
    # (Cc = 50 means a 50 percent reduction).
    # ============================================================

    # ----- Placebo intercept E0 (FIXED; paper Results paragraph 2) -----
    # E0 (-0.802 percent) was held FIXED at the literature value from
    # Mandema 2005 (ref [15]) because the unconstrained sigmoidal fit
    # produced large RSE on E0 (>900 percent per the paper Discussion).
    # The negative value implies a small expected LDL-C increase under
    # placebo in the literature pool; observed placebo arms (Yang 2010
    # Tables 1-2) actually show 0-3.6 percent reductions, so the
    # study-level eta absorbs the deviation.
    e0 <- fixed(-0.802)
    label("Placebo intercept E0 (percent LDL-C reduction at zero dose; FIXED from Mandema 2005)")  # Yang 2010 Table 4 footnote "E0 was presumed to be -0.802 percent based on the literature value [15]"

    # ----- Maximal LDL-C reduction Emax (Yang 2010 Table 4) -----
    lemax <- log(57.0)
    label("Maximal LDL-C reduction Emax (percent)")  # Yang 2010 Table 4 Emax TVP=57.0 (RSE 3.86%)

    # ----- ED50 in Western reference (Yang 2010 Table 4) -----
    # Log-transformed for positivity. exp(led50) recovers 1.74 mg/day.
    led50 <- log(1.74)
    label("ED50 for Western reference (mg/day)")  # Yang 2010 Table 4 ED50 TVP=1.74 mg (RSE 21.8%)

    # ----- Hill exponent / sigmoidicity gamma (FIXED at 1) -----
    # Yang 2010 Results: "gamma was fixed at 1; that is, a simple
    # Emax model was used and the estimated parameters were
    # acceptable". The fixed value drops the model from sigmoidal
    # Emax to simple Emax.
    lhill <- fixed(log(1))
    label("Hill sigmoidicity exponent (canonical lhill; FIXED to 1 per paper Results, reducing the sigmoidal Emax to a simple Emax)")  # Yang 2010 Table 4 "gamma" = 1 (FIXED)

    # ----- Asian race multiplicative effect on ED50 (Yang 2010 Table 4) -----
    # ED50_arm = exp(led50) * e_asian_ed50^RACE_ASIAN.
    # Value 0.564 means ED50 in Asians is about half that in Westerners.
    # RSE 28.55 percent so the estimate has appreciable uncertainty.
    e_asian_ed50 <- 0.564
    label("Asian-race multiplicative effect on ED50 (unitless; ED50_Asian / ED50_Western)")  # Yang 2010 Table 4 theta (race on ED50) TVP=0.564 (RSE 28.55%)

    # ============================================================
    # Inter-trial (study-arm-level) variability. Yang 2010 reports a
    # single inter-trial random effect that enters additively on the
    # predicted percent reduction; Table 4 lists it as a standard
    # deviation (SD 3.0 percent) per the footnote "Inter-trials
    # variability was calculated by taking the square root of eta".
    # nlmixr2 ini() takes the variance, so encode 3.0^2 = 9.0.
    #
    # Naming follows the Boucher_2018_naproxen_mbma precedent
    # (eta_study_<param> for MBMA study-arm-level etas, NOT individual
    # between-subject variability). Adding the eta to a constant E0 is
    # mathematically equivalent to adding it to the predicted Y, since
    # E0 is constant -- so the eta represents a trial-level intercept
    # shift on the percent-reduction scale.
    # ============================================================
    eta_study_e0 ~ 9.0  # Yang 2010 Table 4: inter-trial variability = SD 3.0% -> variance = 3.0^2 = 9.0

    # ============================================================
    # Residual error (Yang 2010 Table 4: "Residual error (SD) = 3.1").
    # Additive on the unsigned percent-reduction scale; arm-level
    # observation precision is treated as homoscedastic at this SD.
    # ============================================================
    addSd <- 3.1
    label("Additive residual SD on study-arm-mean percent LDL-C reduction")  # Yang 2010 Table 4 "Residual error (SD)" = 3.1
  })

  model({
    # ----- Derive structural quantities -----
    emax <- exp(lemax)
    hill <- exp(lhill)

    # ED50 per arm with Asian-race multiplier. RACE_ASIAN is binary
    # (0 = Western reference, 1 = Asian) so the power form collapses to
    # exp(led50) when RACE_ASIAN = 0 and exp(led50) * e_asian_ed50 when
    # RACE_ASIAN = 1, matching the paper covariate equation
    # P = TVP * theta_RACE^RACE.
    ed50 <- exp(led50) * e_asian_ed50^RACE_ASIAN

    # ----- Predicted study-arm mean percent LDL-C reduction -----
    # Sigmoidal Emax form retained explicitly (hill is fixed at 1) so
    # the structural equation matches Yang 2010 Eq 1 verbatim.
    # eta_study_e0 enters additively on the predicted output per the
    # paper's Y_obs = E_pred + eta + epsilon parameterization.
    Cc <- e0 + emax * DOSE^hill / (ed50^hill + DOSE^hill) + eta_study_e0

    Cc ~ add(addSd)
  })
}
