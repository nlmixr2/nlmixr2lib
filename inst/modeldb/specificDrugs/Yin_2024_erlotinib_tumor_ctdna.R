Yin_2024_erlotinib_tumor_ctdna <- function() {
  description <- paste(
    "Tumor-dynamics model of the sum of longest diameters (SLD) of target",
    "lesions in adults with EGFR-mutant non-small cell lung cancer treated",
    "with erlotinib (Yin 2024, model with baseline ctDNA as a covariate).",
    "Same structure as Yin_2024_erlotinib_tumor: a drug-sensitive clone",
    "(`shrink`) seeded at the observed baseline SLD and a drug-resistant",
    "clone (`growth`) starting at zero, both growing at kge; on treatment",
    "the sensitive clone is killed at kse and converts to the resistant",
    "clone at kmut. The growth rate kge depends on the baseline plasma",
    "ctDNA variant allele frequency of the primary EGFR variant: reference",
    "VAF >= 1.74% (the cohort median), multiplied by 0.334 when VAF < 1.74%",
    "and by 0.281 when no baseline sample was available. No PK input;",
    "treatment enters through the on-treatment indicator TRT_PHASE.",
    sep = " "
  )
  reference <- paste(
    "Yin A, Veerman GDM, van Hasselt JGC, Steendam CMJ, Dubbink HJ,",
    "Guchelaar HJ, Friberg LE, Dingemans AMC, Mathijssen RHJ, Moes DJAR.",
    "Quantitative modeling of tumor dynamics and development of drug",
    "resistance in non-small cell lung cancer patients treated with",
    "erlotinib. CPT Pharmacometrics Syst Pharmacol. 2024;13(4):612-623. doi:10.1002/psp4.13105",
    sep = " "
  )
  vignette <- "Yin_2024_erlotinib"
  units <- list(
    time = "day",
    dosing = "n/a (no PK input; erlotinib treatment enters through TRT_PHASE)",
    concentration = "mm (the observable `TS` is the RECIST 1.1 sum of longest diameters)"
  )

  covariateData <- list(
    TUM_SLD = list(
      description = "Observed baseline sum of longest diameters of RECIST 1.1 target lesions.",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      source_name = "baseTS",
      notes = paste(
        "Per-subject regressor, not an estimated parameter: the initial condition of the drug-sensitive clone (Table 3 'TS_0 = Observed baseline'; Appendix S2 'A_0(4)=BASES').",
        "The baseline was fixed to the observed value 'considering the amount of available data' (Methods). Cohort median 76.6 mm, range 29-116 mm (Results).",
        sep = " "
      )
    ),
    TRT_PHASE = list(
      description = "On-treatment indicator: 1 from the start of erlotinib onwards, 0 before it.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (before erlotinib start)",
      source_name = "UNDERTREAT",
      notes = paste(
        "Time-varying data column (Appendix S2 $INPUT 'UNDERTREAT ; if treatment started: 1, yes').",
        "Switches on the kill rate kse, the sensitive-to-resistant transfer kmut and the growth of the resistant clone (Appendix S2 $DES). Before treatment the sensitive clone grows at kge only.",
        "Data were censored after erlotinib was stopped at progression (Results), so the model has no off-treatment regrowth phase after treatment.",
        sep = " "
      )
    ),
    CTDNA_VAF = list(
      description = "Baseline (start-of-treatment) variant allele frequency of the primary sensitising EGFR variant in plasma cell-free DNA.",
      units = "%",
      type = "continuous",
      reference_category = ">= 1.74% (the cohort median)",
      source_name = "baseVAF",
      notes = paste(
        "Primary EGFR variant = the activating variant found in the tumor biopsy (exon 19 deletion, p.L858R or p.K852R), measured by next-generation sequencing of baseline cfDNA (Methods; Results).",
        "Enters as a categorical effect on kge with the cut at 1.74%: the Table 3 row reads 'mutant EGFR VAF < 1.74%' and the commented Appendix S2 code tests 'baseVAF.GE.0.AND.baseVAF.LT.1.74', so a VAF of exactly 1.74% falls in the reference group.",
        "A patient with no primary EGFR variant detected at baseline has VAF = 0 (4 of the 12 with a sample) and falls in the < 1.74% group.",
        "Cohort: available for 12 of 18 patients, median 1.74%, range 0-62.74% (Results).",
        sep = " "
      )
    ),
    CTDNA_VAF_MISSING = list(
      description = "Baseline plasma ctDNA VAF unavailable (no baseline cfDNA sample): 1 = unavailable, 0 = available.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (baseline sample available)",
      source_name = "baseVAF < 0 (Appendix S2 'no sample group')",
      notes = paste(
        "6 of 18 patients had no baseline sample (Results). Set CTDNA_VAF to any finite non-negative placeholder for these patients; the VAF term is gated off when CTDNA_VAF_MISSING = 1.",
        "The missing group's kge factor (0.281) was not significantly different from the VAF < 1.74% group's (0.334) (Results).",
        sep = " "
      )
    )
  )

  covariatesDataExcluded <- list(
    TUM_TP53_MUT = list(
      description = "TP53 mutation detected at baseline (here in plasma ctDNA rather than in a tumor biopsy).",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Tested on kg, km and kd as a categorical covariate (detected / not detected / unavailable) and not significant (Results, Genetic biomarkers); detected in 4 of 18 patients.",
        "Baseline cfDNA concentration (median 1.44 ng/uL) and erlotinib exposure (linear or Emax on kd) were also tested and not retained; neither has a register entry, so they are documented in the vignette only.",
        sep = " "
      )
    )
  )

  compartmentData <- list(
    shrink = list(analyte = "tumor (drug-sensitive clone)", units = "mm", specimen = "tissue", verified = TRUE),
    growth = list(analyte = "tumor (drug-resistant clone)", units = "mm", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 18,
    n_studies = 1,
    age_range = "48-78 years",
    age_median = "66 years",
    weight_range = "46.1-109 kg",
    weight_median = "69.5 kg",
    sex_female_pct = 72.2,
    disease_state = "Non-small cell lung cancer with an activating EGFR mutation (exon 19 deletion, p.L858R or p.K852R), first-line erlotinib",
    dose_range = "75-150 mg orally once daily (all started at 150 mg)",
    regions = "The Netherlands (START-TKI multicentre real-world study, NCT05221372)",
    notes = paste(
      "Yin 2024 Table 1 (PK/PD data, N = 18). 155 SLD measurements, median 7 per patient (range 2-18), over a median 264 days of treatment (range 20-1168).",
      "All 18 patients progressed (16) or died (2); data were censored after erlotinib was stopped.",
      "Baseline SLD median 76.6 mm (range 29-116 mm).",
      sep = " "
    )
  )

  ini({
    lkge <- log(0.00204)
    label("Tumor growth rate constant, both clones, baseline EGFR VAF >= 1.74 (1/day)") # Table 3 'kg (/day)' model with baseline ctDNA = 0.00204 (RSE 25%)
    e_ctdna_vaf_kge <- 0.334
    label("Multiplicative factor on kge when baseline EGFR VAF < 1.74 (unitless)") # Table 3 'f1' = 0.334 (RSE 28%)
    e_ctdna_vaf_missing_kge <- 0.281
    label("Multiplicative factor on kge when the baseline ctDNA sample is unavailable (unitless)") # Table 3 'f2' = 0.281 (RSE 28%)
    lkse <- log(0.0123)
    label("Treatment kill rate constant of the sensitive clone (1/day)") # Table 3 'kd (/day)' model with baseline ctDNA = 0.0123 (RSE 18%)
    lkmut <- log(0.00824)
    label("Sensitive-to-resistant mutation rate constant (1/day)") # Table 3 'km (/day)' model with baseline ctDNA = 0.00824 (RSE 18%)

    # Table 3 reports IIV as CV%; converted to log-normal variances as
    # omega^2 = log(1 + CV^2).
    etalkge ~ 0.027183 # Table 3 IIV kg 16.6% CV -> log(1 + 0.166^2)
    etalkse ~ 0.363423 # Table 3 IIV kd 66.2% CV -> log(1 + 0.662^2)
    etalkmut ~ 0.289112 # Table 3 IIV km 57.9% CV -> log(1 + 0.579^2)

    propSd <- 0.0767
    label("Proportional residual error (fraction)") # Table 3 'Prop. err. (CV%)' model with baseline ctDNA = 7.67 (RSE 14%)
    addSd <- 1.13
    label("Additive residual error (mm)") # Table 3 'Add. err. (SD, mm)' model with baseline ctDNA = 1.13 (RSE 9%)
  })
  model({
    # Appendix S2 (commented covariate lines):
    #   IF(baseVAF.GE.0.AND.baseVAF.LT.1.74) KG1=THETA(6)*KG1
    #   IF(baseVAF.LT.0) KG1=THETA(7)*KG1  ; no sample group
    vaf_low <- (CTDNA_VAF < 1.74) * (1 - CTDNA_VAF_MISSING)
    kge <- exp(lkge + etalkge) *
      e_ctdna_vaf_kge^vaf_low *
      e_ctdna_vaf_missing_kge^CTDNA_VAF_MISSING
    kse <- exp(lkse + etalkse)
    kmut <- exp(lkmut + etalkmut)

    # Methods Eqs 1, 3, 4 with the on-treatment gating of Appendix S2 $DES:
    #   DADT(4) = KG1*A(4) - KD1*UNDERTREAT*A(4) - KM1*UNDERTREAT*A(4)
    #   DADT(5) = KM1*UNDERTREAT*A(4) + KG1*UNDERTREAT*A(5)
    d / dt(shrink) <- kge * shrink - kse * TRT_PHASE * shrink - kmut * TRT_PHASE * shrink
    d / dt(growth) <- kmut * TRT_PHASE * shrink + kge * TRT_PHASE * growth

    # Table 3: TS_0 = observed baseline, TR_0 = 0 fixed
    shrink(0) <- TUM_SLD
    growth(0) <- 0

    TS <- shrink + growth
    TS ~ add(addSd) + prop(propSd)
  })
}
