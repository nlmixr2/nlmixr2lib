Ribba_2022_sld <- function() {
  description <- paste(
    "Stein bi-exponential population model of the tumor-size (sum of the",
    "longest diameters, SLD) time course in advanced non-small cell lung",
    "cancer (NSCLC), fit by Ribba et al. (Roche) to the phase III OAK study.",
    "This is the companion tumor-size fit to Ribba_2022_ctdna.R: the authors",
    "applied the same Stein structural model to SLD independently of the",
    "ctDNA data in order to test whether the two decay-rate constants",
    "correlate (they found r = 0.45, Figure 2F bottom left), which motivated",
    "the coupled model Ribba_2022_ctdna_sld_joint.R. Unlike the ctDNA fit,",
    "SLD was modeled on its natural millimetre scale with a combined",
    "additive-plus-proportional residual-error model, and the fit pooled",
    "both OAK arms (atezolizumab and docetaxel) with treatment arm as a",
    "categorical covariate so that individual parameters could subsequently",
    "be read out for the atezolizumab arm alone. The baseline SLD is",
    "supplied as a per-subject regressor rather than estimated."
  )
  reference <- paste(
    "Ribba B, Roller A, Helms H-J, Stern M, Bleul C.",
    "Circulating tumor DNA: Opportunities and challenges for",
    "pharmacometric approaches.",
    "Front Pharmacol. 2022;13:1058220. doi:10.3389/fphar.2022.1058220.",
    "PMID: 36968790. PMCID: PMC10030934.",
    "Structural model from Stein WD et al., Clin Cancer Res.",
    "2011;17(4):907-917.",
    "Parameter values from Supplementary Data (Appendix), tumor-size table.",
    sep = " "
  )
  vignette <- "Ribba_2022_ctdna_tumor_size"

  units <- list(
    time          = "day",
    dosing        = "n/a (no PK input; treatment effect is absorbed into the empirical growth and decay rate constants)",
    concentration = "mm (the observable `TS` is the RECIST 1.1 sum of the longest diameters of target lesions)"
  )

  covariateData <- list(
    TUM_SLD = list(
      description = "Observed baseline (cycle 1 day 1) sum of the longest diameters of target lesions per RECIST 1.1, used as the Stein baseline regressor SLD0.",
      units       = "mm",
      type        = "continuous",
      source_name = "SLD0",
      notes       = paste(
        "Ribba 2022 Eq. 2 defines SLD0 as the baseline value of the sum of the longest diameters; the Supplementary Data states the baseline was used as a regressor and not estimated, following the same process as the ctDNA fit.",
        "Initialises both Stein sub-states, so TS(0) = TUM_SLD exactly.",
        "The source paper does not report the OAK baseline-SLD distribution; the validation vignette samples a log-normal cohort centred on a typical advanced-NSCLC value and states that assumption explicitly.",
        sep = " "
      )
    )
  )

  # Treatment arm was carried in the source fit as a categorical covariate but no
  # arm-effect coefficient is reported anywhere in the paper or the Appendix, so
  # it cannot be encoded. Documented here rather than in covariateData so the
  # provenance survives without creating a declared-but-unreferenced covariate.
  covariatesDataExcluded <- list(
    TRT = list(
      description = "Treatment arm of the OAK study (atezolizumab vs docetaxel).",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste(
        "Supplementary Data: 'these results were obtained by analyzing simultaneously not only the patients treated with atezolizumab but also the patients treated with docetaxel (the other arm of the study). The arm was treated as a categorical covariate enabling the use of the individual parameters only for the atezolizumab arm.'",
        "No arm-effect coefficient on kge or kse is reported in the paper or the Appendix, so the covariate is documented but not implemented; the packaged parameters are the pooled-arm population values.",
        sep = " "
      )
    )
  )

  population <- list(
    species         = "human (adults with advanced non-small cell lung cancer)",
    n_subjects      = 1225L,
    n_studies       = 1L,
    age_range       = "not reported in this paper (see Rittmeyer 2017 for the OAK cohort demographics)",
    weight_range    = "not reported in this paper",
    sex_female_pct  = NA_real_,
    race_ethnicity  = "not reported in this paper",
    disease_state   = "previously treated locally advanced or metastatic NSCLC (OAK study)",
    dose_range      = "atezolizumab 1200 mg IV every 3 weeks or docetaxel 75 mg/m2 IV every 3 weeks (OAK protocol doses; not model inputs)",
    regions         = "multiregional (OAK was conducted across 31 countries; per-region counts not reported in this paper)",
    notes           = paste(
      "The SLD fit pooled both OAK arms; Ribba 2022 Figure 1B reports 1,225 randomised participants with 613 in the atezolizumab arm. The exact number of participants contributing longitudinal SLD to this fit is not stated in the paper or the Appendix.",
      "Estimation was performed in Monolix 2021R1 (Lixoft) with a combined additive-plus-proportional residual-error model, which the authors report was the best residual-error model for SLD.",
      "The companion ctDNA fit on the same study is Ribba_2022_ctdna.R (n = 46) and the coupled fit is Ribba_2022_ctdna_sld_joint.R.",
      sep = "\n"
    )
  )

  ini({
    # ----- Stein bi-exponential rate constants for SLD (Ribba 2022 Supplementary Data, tumor-size table) -----
    lkge <- log(0.0016)   ; label("log SLD growth rate kge (1/day)")  # Appendix tumor-size table: kgT_pop = 0.0016 1/day, RSE 14.3%
    lkse <- log(0.0014)   ; label("log SLD decay rate kse (1/day)")   # Appendix tumor-size table: ksT_pop = 0.0014 1/day, RSE 29.0%

    # ----- Inter-individual variability (log-normal, uncorrelated) -----
    # Appendix: the SLD data were modeled "using the same model and following the same
    # process as described above", i.e. log-normal individual parameters with no
    # correlation between the two random effects. The table reports the standard
    # deviation of the random effects (omega), so the nlmixr2 variance is omega^2.
    etalkge ~ 1.0609   # Appendix tumor-size table: omega_kgT = 1.03 (RSE 10.7%); 1.03^2 = 1.0609
    etalkse ~ 2.6896   # Appendix tumor-size table: omega_ksT = 1.64 (RSE 15.4%); 1.64^2 = 2.6896

    # ----- Residual error (combined additive + proportional) -----
    addSd  <- 0.65   ; label("Additive residual SD on SLD (mm)")               # Appendix tumor-size table: a = 0.65 mm (RSE 28.1%), constant part of the combined error model
    propSd <- 0.08   ; label("Proportional residual SD on SLD (fraction)")     # Appendix tumor-size table: b = 0.08 (RSE 7.1%), proportional part of the combined error model
  })

  model({
    # Individual Stein rate constants.
    kge <- exp(lkge + etalkge)
    kse <- exp(lkse + etalkse)

    # Ribba 2022 Eq. 2, SLD row (Stein 2011):
    #   SLD(t) = SLD0 * (exp(-ksT*t) + exp(kgT*t) - 1).
    # Encoded as two exponential sub-states with a shared initial condition:
    #   growth(t) = SLD0 * exp(kge * t)
    #   shrink(t) = SLD0 * exp(-kse * t)
    #   TS        = growth + shrink - SLD0
    # At t = 0 this gives SLD0 + SLD0 - SLD0 = SLD0, the Stein boundary condition.
    d/dt(growth) <-  kge * growth
    d/dt(shrink) <- -kse * shrink
    growth(0) <- TUM_SLD
    shrink(0) <- TUM_SLD

    TS <- growth + shrink - TUM_SLD

    TS ~ add(addSd) + prop(propSd)
  })
}
attr(Ribba_2022_sld, "message") <-
  "Stein bi-exponential tumor-size (SLD) model for advanced NSCLC (Ribba 2022; OAK, both arms pooled). Observable `TS` is the RECIST 1.1 sum of longest diameters in mm. Requires the per-subject covariate TUM_SLD (baseline SLD in mm), which initialises both Stein sub-states. No dosing input. Treatment arm was a categorical covariate in the source fit but no arm-effect coefficient is reported, so the packaged parameters are the pooled-arm population values (see covariatesDataExcluded). Companion models: Ribba_2022_ctdna.R and Ribba_2022_ctdna_sld_joint.R."
Ribba_2022_sld
