Shin_2014_sevoflurane <- function() {
  description <- "Pharmacodynamic sigmoid Emax model for the probability of recovery of consciousness (ROC) vs end-tidal sevoflurane concentration (vol %) during emergence from general anesthesia in pediatric dental-surgery patients (Shin 2014). Mentality (intact vs severely mentally disabled, MENT_DISABLED) stratifies both the concentration at 50% probability of ROC (C50) and the Hill coefficient. NONMEM Bernoulli likelihood in the source paper; this implementation exposes the typical-value probability with a placeholder additive residual error (see vignette Assumptions and deviations)."
  reference <- paste(
    "Shin TJ, Noh GJ, Koo YS, Han DW. (2014).",
    "Modeling of Recovery Profiles in Mentally Disabled and Intact Patients",
    "after Sevoflurane Anesthesia; A Pharmacodynamic Analysis.",
    "Yonsei Med J 55(6):1624-1629.",
    "doi:10.3349/ymj.2014.55.6.1624.",
    sep = " "
  )
  vignette <- "Shin_2014_sevoflurane"
  units <- list(
    time          = "minute",
    dosing        = "(none; volatile anesthetic discontinuation at start of emergence)",
    concentration = "vol % (end-tidal sevoflurane, ETSEVO)"
  )

  covariateData <- list(
    MENT_DISABLED = list(
      description        = "Severe mental disability indicator (1 = mentally disabled per the pediatrics-department diagnosis described in Shin 2014 Methods 'Subjects'; 0 = mentally intact, ASA class 1).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (mentally intact)",
      notes              = "Source column MEN in Shin 2014 Appendix 1 ($INPUT MEN). Identical orientation: MEN = 1 mentally disabled, MEN = 0 mentally intact. Time-fixed per subject. New canonical entry registered in inst/references/covariate-columns.md alongside this model.",
      source_name        = "MEN"
    ),
    ETSEVO = list(
      description        = "End-tidal sevoflurane concentration recorded continuously during emergence from general anesthesia (vol %, monitored by GE Healthcare Datex-Ohmeda S5 Collect software per Shin 2014 Methods).",
      units              = "vol % (volume percent in the breathing circuit)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The DOSE column in Shin 2014 Appendix 1 ($PRED 'PROB = 1 - DOSE**GAM/(CE50**GAM + DOSE**GAM)') is the per-record Etsevo concentration -- not an administered dose. Decreases from approximately the 1 MAC maintenance value to 0 during emergence; sampled per observation record (per Shin 2014 Methods, ROC was observed every 15 seconds). New canonical entry registered in inst/references/covariate-columns.md alongside this model.",
      source_name        = "DOSE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20,
    n_studies      = 1,
    age_range      = "3-15 years (median 6 in both groups; Table 1 reports range 2-8 in normal and 2-12 in mentally disabled)",
    weight_range   = "9.5-35 kg (Table 1 median 22 kg normal, 15 kg mentally disabled)",
    sex_female_pct = 45,
    disease_state  = "Pediatric patients scheduled for dental treatment under general anesthesia; ten with severe mental disability (clinical diagnosis from the pediatrics department) and ten ASA class 1 mentally intact controls.",
    dose_range     = "Induction sevoflurane 8 vol % at 6 L/min; maintenance at age-adjusted 1 MAC with 50% N2O until end of dental procedure, then discontinuation. No drug dose modeled.",
    regions        = "Korea (Seoul National University Dental College)",
    notes          = "Total 476 binary ROC observations (Shin 2014 Results). Duration of sevoflurane anesthesia 70 min (60-155) in normal vs 120 min (90-230) in disabled (Table 1; p = 0.08). Hemodynamic variables not significantly different between groups. Excluded prior adverse reaction to anesthetics and cardiac / pulmonary / hepatic / renal disease."
  )

  ini({
    # -----------------------------------------------------------------
    # Fixed-effect estimates and IIV from Shin 2014 Table 2, Final
    # model (lowest OFV = 165.09 after stepwise inclusion of mentality
    # on both C50 and gamma). Source NONMEM control stream reproduced
    # verbatim in Appendix 1; CV(%) on C50 was 12.28 in the Final model
    # (Table 2). IIV on gamma was fixed to 0 ($OMEGA "0 FIX" on
    # ETA(2)). Bernoulli likelihood ($EST LIKELIHOOD LAPLACE METHOD=1
    # SIG=3); see vignette "Assumptions and deviations" for the
    # nlmixr2 typical-value placeholder used here.
    # -----------------------------------------------------------------

    lec50 <- log(0.37)
    label("Log of typical end-tidal sevoflurane concentration at 50% probability of ROC for mentally intact patients (vol %)")
    # Table 2 Final/Intact: C50 = 0.37 vol % (RSE 3.5%); bootstrap 0.34-0.39.
    # Maps to NONMEM Appendix 1 THETA(1).

    e_ment_disabled_ec50 <- log(0.19 / 0.37)
    label("Additive shift on log(C50) when MENT_DISABLED = 1 (log-vol % units); back-transforms to typical C50 = 0.19 vol % in disabled patients")
    # Table 2 Final/Disabled: C50 = 0.19 vol % (RSE 12.5%); bootstrap 0.14-0.24.
    # log(0.19/0.37) = -0.666. Maps to NONMEM Appendix 1 THETA(2)/THETA(1).

    lhill <- log(16.4)
    label("Log of Hill coefficient (steepness of the Etsevo-response curve) for mentally intact patients (unitless)")
    # Table 2 Final/Intact: gamma = 16.4 (RSE 19.0%); bootstrap 13.7-23.7.
    # Maps to NONMEM Appendix 1 THETA(3).

    e_ment_disabled_hill <- log(4.53 / 16.4)
    label("Additive shift on log(Hill) when MENT_DISABLED = 1 (unitless); back-transforms to typical gamma = 4.53 in disabled patients")
    # Table 2 Final/Disabled: gamma = 4.53 (RSE 22.3%); bootstrap median 4.56, 3.0-8.1.
    # log(4.53/16.4) = -1.286. Maps to NONMEM Appendix 1 THETA(4)/THETA(3).

    # IIV: exponential (log-normal) on C50; gamma ETA was FIX at 0
    # in the NONMEM Appendix 1 $OMEGA block (line "0 FIX"), so only
    # etalec50 is estimated. CV% on C50 reported as 12.28% in
    # Table 2 Final row; omega^2 = log(1 + CV^2) = log(1 + 0.1228^2) =
    # 0.01497 on the log-normal-variance scale used in nlmixr2 ini().
    etalec50 ~ 0.01497
    # Table 2 Final CV(%) = 12.28; NONMEM $OMEGA initial 0.1 [P] on ETA(1).

    # Residual-error placeholder. The source likelihood is Bernoulli
    # over the binary ROC observation (Shin 2014 Appendix 1:
    # "Y = R*PROB + (1-R)*(1-PROB)" under $EST LIKELIHOOD LAPLACE).
    # nlmixr2 / rxode2 do not natively express a Bernoulli observation
    # for a probability output in ini() / model(); following the
    # precedent of inst/modeldb/ddmore/Hansson_2013c_sunitinib.R and
    # inst/modeldb/specificDrugs/Sheng_2016_quinine_rat.R, the
    # observation here is the typical-value probability of ROC with a
    # placeholder additive residual the user can tune at fit time.
    # See the vignette "Assumptions and deviations" section.
    addSd_prob_roc <- 0.05
    label("Placeholder additive residual error on the typical-value probability of ROC (unitless probability); not from the source -- see vignette Assumptions")
  })

  model({
    # 1. Derived covariate terms.
    # Per-record end-tidal sevoflurane concentration (vol %). The NONMEM
    # control stream calls this column DOSE; in this implementation it is
    # the canonical ETSEVO covariate. The model uses it directly as the
    # concentration argument to the sigmoid Emax.
    conc <- ETSEVO

    # 2. Individual structural parameters. Mentality enters as a
    # categorical effect on log(C50) and log(Hill); the formulation
    #   exp(lec50 + e_ment_disabled_ec50 * MENT_DISABLED)
    # is mathematically equivalent to the NONMEM Appendix 1 piecewise
    #   CE50 = THETA(1) * (1 - MEN) + THETA(2) * MEN
    # for the two valid MEN values 0 and 1.
    ec50 <- exp(lec50 + e_ment_disabled_ec50 * MENT_DISABLED + etalec50)
    hill <- exp(lhill + e_ment_disabled_hill * MENT_DISABLED)

    # 3. Sigmoid Emax probability of ROC at the current Etsevo. The
    # source-paper form (Shin 2014 Methods, equation immediately above
    # Statistical analysis and Appendix 1 $PRED) is
    #   PROB = 1 - C^gamma / (C50^gamma + C^gamma)
    #        = C50^gamma / (C50^gamma + C^gamma).
    # P(ROC) decreases monotonically as the end-tidal concentration
    # increases (high anesthetic -> unconscious).
    prob_roc <- ec50^hill / (ec50^hill + conc^hill)

    # 4. Observation. Placeholder additive residual on the typical-value
    # probability; see ini() block comment and vignette Assumptions and
    # deviations for the deviation rationale relative to the source
    # Bernoulli likelihood.
    prob_roc ~ add(addSd_prob_roc)
  })
}
