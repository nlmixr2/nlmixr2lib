Ribba_2022_ctdna <- function() {
  description <- paste(
    "Stein bi-exponential population model of the circulating tumor DNA",
    "(ctDNA) time course in advanced non-small cell lung cancer (NSCLC),",
    "fit by Ribba et al. (Roche) to the atezolizumab arm of the phase III",
    "OAK study (n = 46 participants with serial ctDNA). ctDNA was measured",
    "with the Roche AVENIO next-generation-sequencing panel and reported as",
    "average mutant molecules per millilitre of plasma (MMPM) at four",
    "cycles (baseline, ~21, ~42 and ~63 days). The data were base-10 log",
    "transformed before fitting, so the modeled quantity is log10(MMPM) and",
    "the Stein baseline y0 is log10 of the observed baseline MMPM, supplied",
    "as a per-subject regressor rather than estimated. The model decomposes",
    "the ctDNA signal into an exponentially growing and an exponentially",
    "decaying fraction; the population nadir time log(kse/kge)/(kge+kse) =",
    "63.6 days = 9.1 weeks underpins the paper's conclusion that sampling",
    "at cycle 2 (21 days) is too early to capture the maximal ctDNA",
    "response. Inter-individual variability on both rate constants is close",
    "to 100% CV, which the authors highlight as the central challenge for",
    "ctDNA-informed decision making."
  )
  reference <- paste(
    "Ribba B, Roller A, Helms H-J, Stern M, Bleul C.",
    "Circulating tumor DNA: Opportunities and challenges for",
    "pharmacometric approaches.",
    "Front Pharmacol. 2022;13:1058220. doi:10.3389/fphar.2022.1058220.",
    "PMID: 36968790. PMCID: PMC10030934.",
    "Structural model from Stein WD et al., Clin Cancer Res.",
    "2011;17(4):907-917.",
    "Parameter values from Supplementary Data (Appendix), ctDNA table.",
    sep = " "
  )
  vignette <- "Ribba_2022_ctdna_tumor_size"

  units <- list(
    time          = "day",
    dosing        = "n/a (no PK input; treatment effect is absorbed into the empirical growth and decay rate constants of the atezolizumab arm)",
    concentration = "log10(MMPM) -- the observable `ctdna` is base-10 log-transformed average mutant molecules per mL of plasma, matching the scale on which Ribba 2022 fit the model; the residual-error parameter addSd is therefore in log10 units"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    growth_ctdna = list(analyte = "ctDNA", units = NA_character_, specimen = "plasma", verified = FALSE),
    shrink_ctdna = list(analyte = "ctDNA", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CTDNA = list(
      description = "Observed baseline (cycle 1 day 1) circulating tumor DNA burden, used as the Stein baseline regressor y0 after base-10 log transformation.",
      units       = "MMPM (mutant molecules per mL of plasma)",
      type        = "continuous",
      source_name = "y0",
      notes       = paste(
        "Ribba 2022 Eq. 1 and Supplementary Data: 'The baseline value was used as a regressor (not estimated).'",
        "Assayed with the Roche AVENIO panel in the OAK cohort.",
        "Enters model() as rbase_ctdna <- log10(CTDNA) and initialises both Stein sub-states, so ctdna(0) = log10(CTDNA) exactly.",
        "Ribba 2022 Figure 2F (top left) plots the model prediction normalised by this baseline, i.e. with rbase_ctdna set to 1, which is why every simulated curve in that panel starts at 1.",
        sep = " "
      )
    )
  )

  population <- list(
    species         = "human (adults with advanced non-small cell lung cancer)",
    n_subjects      = 46L,
    n_studies       = 1L,
    age_range       = "not reported in this paper (see Rittmeyer 2017 for the OAK cohort demographics)",
    weight_range    = "not reported in this paper",
    sex_female_pct  = NA_real_,
    race_ethnicity  = "not reported in this paper",
    disease_state   = "previously treated locally advanced or metastatic NSCLC (OAK study)",
    dose_range      = "atezolizumab 1200 mg IV every 3 weeks (OAK protocol dose; not a model input)",
    regions         = "multiregional (OAK was conducted across 31 countries; per-region counts not reported in this paper)",
    notes           = paste(
      "The ctDNA sub-cohort is the 46 atezolizumab-arm OAK participants with serial ctDNA measurements, out of 613 atezolizumab-arm participants and 1,225 randomised overall (Ribba 2022 Figure 1B).",
      "ctDNA sampling times: cycle 1 (baseline), cycle 2 (~21 days), cycle 3 (~42 days), cycle 4 (~63 days).",
      "Estimation was performed in Monolix 2021R1 (Lixoft) with a constant residual-error model on the log10 scale.",
      "The companion sum-of-longest-diameters fit on the same study is Ribba_2022_sld.R and the coupled fit is Ribba_2022_ctdna_sld_joint.R.",
      sep = "\n"
    )
  )

  ini({
    # ----- Stein bi-exponential rate constants for ctDNA (Ribba 2022 Supplementary Data, ctDNA table) -----
    lkge_ctdna <- log(0.0038)   ; label("log ctDNA growth rate kge_ctdna (1/day)")  # Appendix ctDNA table: kg_pop = 0.0038 1/day, RSE 30.1%
    lkse_ctdna <- log(0.0081)   ; label("log ctDNA decay rate kse_ctdna (1/day)")   # Appendix ctDNA table: ks_pop = 0.0081 1/day, RSE 27.4%

    # ----- Inter-individual variability (log-normal, uncorrelated) -----
    # Appendix: "Individual parameters, kg_i and ks_i were assumed to be log-normally
    # distributed" and "We assumed no correlation between the random effects associated
    # to the two parameters." The table reports the standard deviation of the random
    # effects (omega), so the nlmixr2 variance is omega^2.
    etalkge_ctdna ~ 0.6724   # Appendix ctDNA table: omega_kg = 0.82 (RSE 25.0%); 0.82^2 = 0.6724
    etalkse_ctdna ~ 1.0404   # Appendix ctDNA table: omega_ks = 1.02 (RSE 19.9%); 1.02^2 = 1.0404

    # ----- Residual error (constant model on the log10 scale) -----
    addSd <- 0.27   ; label("Constant residual SD on log10 ctDNA (log10 MMPM)")  # Appendix ctDNA table: a = 0.27 (RSE 4.8%), constant error model
  })

  model({
    # Individual Stein rate constants.
    kge_ctdna <- exp(lkge_ctdna + etalkge_ctdna)
    kse_ctdna <- exp(lkse_ctdna + etalkse_ctdna)

    # Stein baseline y0 on the modeling scale. Ribba 2022 log10-transformed the
    # MMPM data before fitting, so the regressor that multiplies the bracket is
    # log10 of the observed baseline MMPM, not the baseline MMPM itself.
    rbase_ctdna <- log10(CTDNA)

    # Ribba 2022 Eq. 1 (Stein 2011): y = y0 * (exp(-ks*t) + exp(kg*t) - 1).
    # Encoded as two exponential sub-states with a shared initial condition:
    #   growth_ctdna(t) = y0 * exp(kge_ctdna * t)
    #   shrink_ctdna(t) = y0 * exp(-kse_ctdna * t)
    #   ctdna           = growth_ctdna + shrink_ctdna - y0
    # At t = 0 this gives y0 + y0 - y0 = y0, the Stein boundary condition.
    d/dt(growth_ctdna) <-  kge_ctdna * growth_ctdna
    d/dt(shrink_ctdna) <- -kse_ctdna * shrink_ctdna
    growth_ctdna(0) <- rbase_ctdna
    shrink_ctdna(0) <- rbase_ctdna

    ctdna <- growth_ctdna + shrink_ctdna - rbase_ctdna

    ctdna ~ add(addSd)
  })
}
attr(Ribba_2022_ctdna, "message") <-
  "Stein bi-exponential ctDNA model for advanced NSCLC (Ribba 2022; OAK atezolizumab arm, n = 46). Observable `ctdna` is on the base-10 LOG scale: ctdna = log10(MMPM), and addSd = 0.27 is in log10 units. Requires the per-subject covariate CTDNA (baseline MMPM, untransformed); the model derives rbase_ctdna <- log10(CTDNA) internally. No dosing input -- treatment effect is absorbed into the empirical rate constants. Population nadir at log(kse/kge)/(kge+kse) = 63.6 days = 9.1 weeks. IIV is close to 100% CV on both rates, as emphasised by the authors."
Ribba_2022_ctdna
