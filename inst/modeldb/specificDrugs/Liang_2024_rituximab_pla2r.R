Liang_2024_rituximab_pla2r <- function() {
  description <- "Empirical mono-exponential decline of the serum anti-PLA2R autoantibody titer after rituximab in adults with primary membranous nephropathy; fitted independently of the CD20+ B cell count, so it is a separate model from the companion QSS TMDD PK/PD model (Liang 2024)"
  reference <- "Liang H, Deng Z, Niu S, Kong W, Liu Y, Wang S, Li H, Wang Y, Zheng D, Liu D. Dosing optimization of rituximab for primary membranous nephropathy by population pharmacokinetic and pharmacodynamic study. Front Pharmacol. 2024 Mar 26;15:1197651. doi:10.3389/fphar.2024.1197651"
  vignette <- "Liang_2024_rituximab"

  # The paper reports ke,PLA2R per DAY ("mean ke,PLA2R is 0.033 +/- 0.017,
  # corresponding to the mean half-life 21 days"; log(2)/0.033 = 21.0 days), so
  # this model runs on a day time base -- unlike the companion QSS TMDD model
  # Liang_2024_rituximab, which runs on hours because Table 2 is in 1/hr. There
  # is no dosing record: time is measured from the first rituximab infusion and
  # the decline is described empirically, not driven by a drug concentration.
  units <- list(time = "day", dosing = "none", concentration = "U/mL")

  compartmentData <- list(
    # Serum anti-PLA2R titer assayed by ELISA (Euroimmun, Lubeck, Germany);
    # Liang 2024 Section 2.2.
    antipla2r = list(
      analyte  = "anti-PLA2R IgG autoantibody titer",
      units    = "U/mL",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list()

  # Pearson correlation was run against the available physiological factors to
  # explain ke,PLA2R and none reached the paper's |r| > 0.7 screening threshold
  # (Liang 2024 Section 3.6), so the model carries no covariate terms. The
  # covariates the paper names elsewhere are recorded for provenance only.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      notes = "Among the physiological factors screened by Pearson correlation against ke,PLA2R; no factor reached |r| > 0.7 (Liang 2024 Section 3.6). Population range 19-76 years (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "binary",
      type = "categorical",
      notes = "Among the physiological factors screened against ke,PLA2R; no significant relationship identified (Liang 2024 Section 3.6). 10/41 (24.4%) female (Table 1)."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (BSA-normalized renal function)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Among the physiological factors screened against ke,PLA2R; no significant relationship identified (Liang 2024 Section 3.6). Baseline 92.0 +/- 24.6 mL/min/1.73 m^2 (Table 1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 36L,
    n_studies      = 1L,
    age_range      = "19-76 years (mean 52.8 +/- 14.9) in the parent cohort of 41",
    weight_range   = "50-101 kg (mean 75.4 +/- 11.3) in the parent cohort of 41",
    sex_female_pct = 24.4,
    disease_state  = "Adults with primary membranous nephropathy (PMN) treated off-label with rituximab. Baseline anti-PLA2R antibody titer 260.3 +/- 453.2 U/mL (range 5.4-2695) with 17/41 (41.5%) above the 150 U/mL high-titer threshold (Table 1).",
    renal_function = "eGFR 92.0 +/- 24.6 mL/min/1.73 m^2 (range 32-151) in the parent cohort of 41.",
    dose_range     = "Most patients received a monthly mini-dose of 100 mg IV rituximab; 10/41 received 200-500 mg for some months. The titer model itself is not dose-driven.",
    regions        = "Single centre: Department of Nephrology, Peking University Third Hospital, Beijing, China; retrospective cohort treated March 2019 to December 2021. Registered as ChiCTR2200057381.",
    notes          = paste(
      "Anti-PLA2R titers from 36 of the 41 patients were analysable (Liang 2024 Section 3.6); 276 titer",
      "observations were collected in total (Section 3.2). Only the DESCENDING stage was fitted: serologic",
      "relapse (defined as two consecutive ascending titers), titers separated by more than 3 months, and",
      "titers after 12 months from the first infusion were all excluded (Section 2.7). Patients who relapsed",
      "and then underwent a second round of titer reduction were fitted with two separate exponential",
      "functions. The fit gave R^2 > 0.8 in every patient but one (R^2 = 0.67)."
    )
  )

  # ---------------------------------------------------------------------------
  # Why this is a separate model file from Liang_2024_rituximab
  #
  # The authors tried and failed to link the anti-PLA2R titer to the CD20+ B cell
  # count of their PK/PD model: "we failed to establish the quantitatively
  # relationship between CD20+ B cell counts and anti-PLA2R titers" (Liang 2024
  # Section 3.6), because monthly sampling could not capture the depletion
  # transient and because titer re-emergence after B cell recovery was
  # inconsistent across patients. The titers "of 36 in all 41 patients were
  # analyzed without the integration of CD20+ B cell count", giving an empirical
  # decay whose rate constant is stated to be "independent of CD20+ B cell
  # counts". That makes it a second, non-hierarchical model on a different
  # endpoint, so it is packaged as its own file per the library's
  # replicate-the-author's-structure policy, sharing one vignette with its
  # companion.
  #
  # Distribution of the two parameters (the one derived step in this model)
  #
  # The paper reports both parameters as arithmetic mean +/- SD over individual
  # fits, not as a fitted population model. Encoding them as log-normal by
  # moment matching -- median = mean / sqrt(1 + CV^2) and omega^2 = log(1 + CV^2)
  # -- makes a simulated cohort reproduce the reported mean and SD by
  # construction, and is confirmed by two further statistics the paper reports
  # but that were NOT used to set the values:
  #   * ke,PLA2R: the resulting 95% interval is 0.0113-0.0759 /day against the
  #     reported observed range of 0.01-0.079 /day (Section 3.6).
  #   * initial titer: the resulting distribution puts 45.1% of patients above
  #     the 150 U/mL high-titer threshold against the reported 17/41 = 41.5%
  #     (Table 1). Reading the arithmetic mean 260.3 U/mL as the median instead
  #     gives 68%, which is decisively refuted by that count.
  # A normal (additive) reading of the "+/-" is refuted independently: mean 0.033
  # with SD 0.017 would put the lower 2.5% of ke,PLA2R at approximately zero and
  # admit negative rate constants, whereas the observed minimum is 0.01.
  # Because the mean is reproduced by construction, log(2) / mean(ke,PLA2R) still
  # returns the paper's headline 21-day mean half-life, computed the same way the
  # paper computed it.
  # ---------------------------------------------------------------------------

  ini({
    # Initial (pre-treatment) anti-PLA2R titer, the "A" of Liang 2024 Eq. 5. The
    # paper does not report a fitted population value, so the cohort baseline
    # titer of Table 1 is used; see the note above for the moment matching.
    le0   <- log(129.64);  label("Initial anti-PLA2R titer (A, U/mL)")                        # Liang 2024 Table 1: baseline anti-PLA2R titer 260.3 +/- 453.2 U/mL; log-normal median = 260.3 / sqrt(1 + (453.2/260.3)^2) = 129.64
    lkel  <- log(0.029336); label("Anti-PLA2R elimination rate constant (ke,PLA2R, 1/day)")   # Liang 2024 Section 3.6: mean ke,PLA2R = 0.033 +/- 0.017 /day; log-normal median = 0.033 / sqrt(1 + (0.017/0.033)^2) = 0.029336

    # Between-subject variability, log-normal, from the reported SDs.
    etale0  ~ 1.394093                                                                       # Liang 2024 Table 1: omega^2 = log(1 + (453.2/260.3)^2) = 1.394093
    etalkel ~ 0.235373                                                                       # Liang 2024 Section 3.6: omega^2 = log(1 + (0.017/0.033)^2) = 0.235373

    # The individual fits are reported only as R^2 (> 0.8 in every patient but
    # one, which gave 0.67; Liang 2024 Section 3.6). No residual standard
    # deviation is reported anywhere in the paper or its Supplementary Material,
    # so it is encoded as zero rather than invented. See the vignette Errata.
    propSd <- fixed(0); label("Proportional residual standard deviation for anti-PLA2R titer (fraction)")
  })

  model({
    e0  <- exp(le0 + etale0)
    kel <- exp(lkel + etalkel)

    # Liang 2024 Eq. 5: titer(t) = A * exp(-ke,PLA2R * t), with t measured in
    # days from the start of the descending stage. Written as the equivalent
    # first-order ODE so the titer is a proper compartment.
    antipla2r(0) <- e0
    d/dt(antipla2r) <- -kel * antipla2r

    antipla2r ~ prop(propSd)
  })
}
