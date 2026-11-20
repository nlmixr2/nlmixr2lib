Perello_2024_hexasodiumFytate_cac <- function() {
  description <- paste0(
    "Simple Emax exposure-efficacy model relating the individual average ",
    "end-of-infusion plasma concentration (Cmax) of hexasodium fytate ",
    "(SNF472, measured as free acid IP6) to the change from baseline to ",
    "week 52 in natural-log coronary artery calcium (CAC) volume score, in ",
    "adults on maintenance hemodialysis with cardiovascular calcification ",
    "(Perello 2024, phase 2b CaLIPSO per-protocol subset; placebo, 300 mg or ",
    "600 mg thrice weekly for 52 weeks). The authors constrained Emax = -E0, ",
    "so the change in log CAC volume is E0 * EC50 / (EC50 + Cmax): E0 = ",
    "log(1.169) (16.9% progression without drug), EC50 = 12.2 uM, and ",
    "progression vanishes at saturating exposure (Table 5). The published ",
    "model also adjusted for baseline log CAC volume, but that coefficient ",
    "and its centring value are not reported; this file encodes the ",
    "typical-baseline curve (adjustment at its centring value), which is the ",
    "'Model Estimates' curve of Figure 4B. Additive residual on the log ",
    "scale; no between-subject random effect. No PK layer and no ODE: the ",
    "exposure metric is supplied as the CMAX data column. The companion ",
    "PK/PD model is modellib('Perello_2024_hexasodiumFytate_hydroxyapatite')."
  )
  reference <- paste(
    "Perello J, Alberti J, Torres JV, Ferrer MD, Perez MM, Bassissi F,",
    "Gold A, Raggi P, Chertow GM, Salcedo C.",
    "Hexasodium fytate exposure-response correlations in a randomized,",
    "placebo-controlled study of patients on dialysis with cardiovascular",
    "calcification. Front Pharmacol. 2024;15:1325186.",
    "doi:10.3389/fphar.2024.1325186.",
    sep = " "
  )
  vignette <- "Perello_2024_hexasodiumFytate"
  units <- list(
    time = "n/a (static landmark regression at week 52; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CMAX covariate column)",
    concentration = "lncacvolcfb (change from baseline to week 52 in natural-log CAC volume score, unitless log fold-change; 100 * (exp(lncacvolcfb) - 1) is the percent change)"
  )

  covariateData <- list(
    CMAX = list(
      description = "Individual average hexasodium fytate end-of-infusion plasma concentration (Cmax) over 52 weeks. Supplied as data: this model has no PK layer.",
      units = "uM",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TOTAL plasma IP6 (hexasodium fytate free acid, molecular weight",
        "660 g/mol) measured by LC-MS/MS within the last 10 min of a 2.5-h",
        "intradialytic infusion at weeks 1, 10, 22 and 52, averaged per",
        "subject -- an OBSERVED end-of-infusion concentration, not a model",
        "prediction. No accumulation over 52 weeks (Results 3.2), so",
        "effectively a steady-state Cmax. Values below the lower limit of",
        "quantification (0.758 uM) were set to LQ/2 = 0.379 uM (Table 5",
        "footnote); every placebo subject carries that value. Cmax quartile",
        "averages 0.4, 9.0, 17.4, 32.0 and 60.1 uM appear in Figure 4B."
      ),
      source_name = "Cmax (average Cmax over 52 weeks)"
    )
  )

  covariatesDataExcluded <- list(
    CAC_VOL_BL = list(
      description = "Baseline natural-log coronary artery calcium volume score.",
      units = "log(mm^3)",
      type = "continuous",
      notes = paste(
        "RETAINED in the published final model ('adjusting by baseline log",
        "CAC volume', Results 3.5; it is the third of the three Model DF in",
        "Table 5) but its coefficient and centring value are NOT reported",
        "anywhere in Perello 2024, which has no supplement. The term is",
        "therefore omitted here and the packaged model is the curve at the",
        "centring value, the 'Model Estimates' line of Figure 4B. Listed for",
        "provenance only; not a register canonical."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 55L,
    n_studies = 1L,
    n_observations = "one week-52 change score per subject; Table 5 reports 56 observations (Model DF 3 + Error DF 53, uncorrected total) against 55 per-protocol subjects in Results 3.1 -- the discrepancy is not explained in the paper",
    disease_state = "adults with kidney failure on maintenance hemodialysis and coronary artery calcification (baseline Agatston CAC score 100-3,500), per-protocol PK/PD and CAC subset of the phase 2b CaLIPSO study (NCT02966028)",
    dose_range = "placebo (15), hexasodium fytate 300 mg (20) or 600 mg (20) infused over 2.5 h into the dialysis line thrice weekly for 52 weeks",
    notes = paste0(
      "Demographics of this subset are not reported in Perello 2024. CAC ",
      "volume was scored centrally from non-contrast multi-detector CT at ",
      "baseline and week 52 (Methods 2.1); the primary CaLIPSO efficacy ",
      "endpoint was the change in log CAC volume score. Adjusted R2 of the ",
      "final model 0.1857 (Results 3.5)."
    )
  )

  ini({
    # Perello 2024 Table 5 ('Final simple model for the relation between
    # average Cmax over 52 Weeks and percentage change from baseline in
    # coronary artery calcium volume score') and Results 3.5.
    #
    # Table 5 prints BACK-TRANSFORMED values, 100 * (exp(theta) - 1):
    #   E0   = 16.9 %  (95% CI 6.6, 28.2)
    #   Emax = -14.5 % (95% CI -22.0, -6.19)
    # On the log scale the two CIs are mirror images:
    #   E0:   log(1.066) = 0.0639, log(1.169) = 0.1562, log(1.282) = 0.2484
    #   Emax: log(0.780) = -0.2485, log(0.855) = -0.1567, log(0.9381) = -0.0639
    # and both rows carry t = 3.39, which is the 'Emax = -E0' constraint
    # stated in Results 3.5. Emax is therefore not a separate parameter.
    e0 <- log(1.169)
    label("Change in natural-log CAC volume score from baseline to week 52 at zero exposure, E0 (log fold-change)") # Table 5, E0 = 16.9 % CAC change back-transformed (95% CI 6.6, 28.2); log(1 + 0.169)
    lec50 <- log(12.2)
    label("Average Cmax producing half of the maximal attenuation of CAC progression, EC50 (uM)") # Table 5, EC50 = 12.2 uM, SE 12.1, 95% CI -12.1 to 36.5 (imprecise; p = 0.319)

    # Additive residual on the log scale. SD = sqrt(Error SS / Error DF) =
    # sqrt(1.4 / 53) = 0.163. The Error SS is printed to one decimal, so
    # the SD is only known to about +/- 0.003.
    addSd <- 0.163
    label("Additive residual error SD on the change in natural-log CAC volume score (log fold-change)") # Table 5, Error SS 1.4 on 53 DF; sqrt(1.4/53) = 0.1625
  })

  model({
    ec50 <- exp(lec50)
    # Emax = -E0 (Results 3.5): E0 + Emax * CMAX / (EC50 + CMAX)
    emax <- -e0

    lncacvolcfb <- e0 + emax * CMAX / (ec50 + CMAX)

    lncacvolcfb ~ add(addSd)
  })
}
