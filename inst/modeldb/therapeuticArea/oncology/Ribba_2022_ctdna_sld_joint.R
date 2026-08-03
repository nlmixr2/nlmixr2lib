Ribba_2022_ctdna_sld_joint <- function() {
  description <- paste(
    "Joint Stein bi-exponential model of circulating tumor DNA (ctDNA) and",
    "tumor size (sum of the longest diameters, SLD) in advanced non-small",
    "cell lung cancer, fit by Ribba et al. (Roche) to the phase III OAK",
    "study. Two Stein sub-systems run in parallel -- one on SLD in mm, one",
    "on base-10 log-transformed ctDNA in MMPM -- and are coupled through a",
    "single estimated link parameter zeta: the ctDNA decay rate is not free",
    "but is set to zeta times the SLD decay rate, encoding the hypothesis",
    "that both biomarkers are driven by one underlying treatment-response",
    "process. Every other population parameter (both growth rates, the SLD",
    "decay rate, and their between-subject variances) is FIXED to the value",
    "obtained from the two independent single-endpoint fits, Ribba_2022_sld.R",
    "and Ribba_2022_ctdna.R; only zeta, its variance, and the two constant",
    "residual-error terms were estimated. The estimated zeta of 1.94 means",
    "ctDNA falls about twice as fast as tumor size. The authors present the",
    "resulting individual fits as illustration only and state explicitly",
    "that the parameters were not estimated with sufficient precision for",
    "predictive use; this model is packaged for structural reproduction of",
    "Ribba 2022 Figure 2F (right), not for prediction."
  )
  reference <- paste(
    "Ribba B, Roller A, Helms H-J, Stern M, Bleul C.",
    "Circulating tumor DNA: Opportunities and challenges for",
    "pharmacometric approaches.",
    "Front Pharmacol. 2022;13:1058220. doi:10.3389/fphar.2022.1058220.",
    "PMID: 36968790. PMCID: PMC10030934.",
    "Structural model from Stein WD et al., Clin Cancer Res.",
    "2011;17(4):907-917.",
    "Parameter values from Supplementary Data (Appendix), joint-model table",
    "plus the two single-endpoint tables from which the fixed values are taken.",
    sep = " "
  )
  vignette <- "Ribba_2022_ctdna_tumor_size"

  units <- list(
    time          = "day",
    dosing        = "n/a (no PK input; treatment effect is absorbed into the empirical growth and decay rate constants)",
    concentration = "two outputs on different scales -- `TS` is the RECIST 1.1 sum of longest diameters in mm, and `ctdna` is base-10 log-transformed average mutant molecules per mL of plasma (log10 MMPM). The residual-error parameters follow their outputs: addSd_TS is in mm, addSd_ctdna is in log10 units."
  )

  covariateData <- list(
    TUM_SLD = list(
      description = "Observed baseline (cycle 1 day 1) sum of the longest diameters of target lesions per RECIST 1.1, used as the Stein baseline regressor SLD0.",
      units       = "mm",
      type        = "continuous",
      source_name = "SLD0",
      notes       = "Ribba 2022 Eq. 2: SLD0 is the baseline value of the sum of longest diameters, supplied as a regressor and not estimated. Initialises both SLD Stein sub-states, so TS(0) = TUM_SLD exactly."
    ),
    CTDNA = list(
      description = "Observed baseline (cycle 1 day 1) circulating tumor DNA burden, used as the Stein baseline regressor ctDNA0 after base-10 log transformation.",
      units       = "MMPM (mutant molecules per mL of plasma)",
      type        = "continuous",
      source_name = "ctDNA0",
      notes       = paste(
        "Ribba 2022 Eq. 2: ctDNA0 is the baseline ctDNA value, supplied as a regressor and not estimated.",
        "Assayed with the Roche AVENIO panel in the OAK cohort.",
        "Enters model() as rbase_ctdna <- log10(CTDNA) because the source paper log10-transformed the MMPM data before fitting, so ctdna(0) = log10(CTDNA) exactly.",
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
      "The joint fit requires both endpoints on the same participant, so it is limited to the 46 atezolizumab-arm OAK participants with serial ctDNA and longitudinal SLD (Ribba 2022 Figure 1B).",
      "Estimation was performed in Monolix 2021R1 (Lixoft). Only zeta, omega_zeta and the two constant residual-error terms were estimated; all remaining population parameters were fixed to the two independent single-endpoint fits.",
      "The authors state that the model 'can fit different types of profiles although parameters of the model were not estimated with sufficient precision to use this model for any predictive purpose' (Ribba 2022 Figure 2F caption).",
      sep = "\n"
    )
  )

  ini({
    # ----- SLD Stein rate constants: FIXED to the independent SLD fit (Ribba_2022_sld.R) -----
    # Appendix: "For the joint ctDNA and SLD model, all population parameters (fixed and
    # random effects) were fixed to the values reported above."
    lkge <- fixed(log(0.0016))   ; label("log SLD growth rate kge (1/day), from the independent SLD fit")  # Appendix tumor-size table: kgT_pop = 0.0016 1/day
    lkse <- fixed(log(0.0014))   ; label("log SLD decay rate kse (1/day), from the independent SLD fit")   # Appendix tumor-size table: ksT_pop = 0.0014 1/day

    # ----- ctDNA growth rate: FIXED to the independent ctDNA fit (Ribba_2022_ctdna.R) -----
    # Eq. 2 reuses the symbol kg from Eq. 1 for the ctDNA growth rate. The ctDNA DECAY
    # rate ks from Eq. 1 is NOT reused: in the joint model it is replaced by zeta * kse.
    lkge_ctdna <- fixed(log(0.0038)) ; label("log ctDNA growth rate kge_ctdna (1/day), from the independent ctDNA fit")  # Appendix ctDNA table: kg_pop = 0.0038 1/day

    # ----- Cross-endpoint link parameter (the only estimated structural parameter) -----
    lzeta <- log(1.94)   ; label("log zeta -- dimensionless multiplier linking the SLD decay rate to the ctDNA decay rate (kse_ctdna = zeta * kse)")  # Appendix joint-model table: zeta_pop = 1.94, RSE 37.3%

    # ----- Inter-individual variability -----
    # The three fixed etas carry the omegas of the two single-endpoint fits; only
    # omega_zeta was estimated. Table values are standard deviations of the random
    # effects, so the nlmixr2 variance is omega^2.
    etalkge       ~ fixed(1.0609)   # Appendix tumor-size table: omega_kgT = 1.03; 1.03^2 = 1.0609 (from the joint fit)
    etalkse       ~ fixed(2.6896)   # Appendix tumor-size table: omega_ksT = 1.64; 1.64^2 = 2.6896 (from the joint fit)
    etalkge_ctdna ~ fixed(0.6724)   # Appendix ctDNA table: omega_kg = 0.82; 0.82^2 = 0.6724 (from the joint fit)
    etalzeta      ~ 0.7396          # Appendix joint-model table: omega_zeta = 0.86, RSE 35.0%; 0.86^2 = 0.7396 (estimated)

    # ----- Residual error: constant on both outputs (Appendix: "the parameters of the two error models (assumed constant)") -----
    addSd_TS    <- 7.91    ; label("Constant residual SD on SLD (mm)")                     # Appendix joint-model table: a (SLD) = 7.91 mm, RSE 5.07%
    addSd_ctdna <- 0.024   ; label("Constant residual SD on log10 ctDNA (log10 MMPM)")     # Appendix joint-model table: a (ctDNA) = 0.024, RSE 6.25%
  })

  model({
    # ----- Individual parameters -----
    kge  <- exp(lkge  + etalkge)
    kse  <- exp(lkse  + etalkse)
    kge_ctdna <- exp(lkge_ctdna + etalkge_ctdna)
    zeta <- exp(lzeta + etalzeta)

    # The coupling. Ribba 2022 Eq. 2: the ctDNA decay exponent is zeta * ksT, i.e. the
    # SLD decay rate scaled by the link parameter -- "ksT is the decay rate which we
    # found again in the equation for ctDNA. The parameter zeta links the time course
    # dynamics of SLD and ctDNA data."
    kse_ctdna <- zeta * kse

    # Stein baselines on their respective modeling scales. SLD is modeled in mm;
    # ctDNA was log10-transformed before fitting, so the ctDNA regressor is log10(MMPM).
    rbase_ctdna <- log10(CTDNA)

    # ----- SLD sub-system: SLD(t) = SLD0 * (exp(-kse*t) + exp(kge*t) - 1) -----
    d/dt(growth) <-  kge * growth
    d/dt(shrink) <- -kse * shrink
    growth(0) <- TUM_SLD
    shrink(0) <- TUM_SLD

    # ----- ctDNA sub-system: ctDNA(t) = ctDNA0 * (exp(-zeta*kse*t) + exp(kge_ctdna*t) - 1) -----
    d/dt(growth_ctdna) <-  kge_ctdna * growth_ctdna
    d/dt(shrink_ctdna) <- -kse_ctdna * shrink_ctdna
    growth_ctdna(0) <- rbase_ctdna
    shrink_ctdna(0) <- rbase_ctdna

    TS    <- growth       + shrink       - TUM_SLD
    ctdna <- growth_ctdna + shrink_ctdna - rbase_ctdna

    TS    ~ add(addSd_TS)
    ctdna ~ add(addSd_ctdna)
  })
}
attr(Ribba_2022_ctdna_sld_joint, "message") <-
  "Joint ctDNA + tumor-size Stein model for advanced NSCLC (Ribba 2022; OAK atezolizumab arm, n = 46). Two outputs on DIFFERENT scales: `TS` in mm and `ctdna` in log10(MMPM). Coupled through kse_ctdna = zeta * kse with zeta = 1.94 (RSE 37.3%) -- ctDNA falls about twice as fast as tumor size. All other population parameters are fixed() to the two single-endpoint fits (Ribba_2022_sld.R, Ribba_2022_ctdna.R); only zeta, omega_zeta and the two constant residual SDs were estimated. Requires per-subject covariates TUM_SLD (mm) and CTDNA (baseline MMPM). The authors state the parameters were not estimated with sufficient precision for predictive use, so treat this as a structural reproduction of Figure 2F (right), not a predictive model."
Ribba_2022_ctdna_sld_joint
