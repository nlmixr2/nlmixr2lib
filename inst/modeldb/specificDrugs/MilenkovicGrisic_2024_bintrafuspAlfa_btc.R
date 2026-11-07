MilenkovicGrisic_2024_bintrafuspAlfa_btc <- function() {
  description <- "Empirical tumor-size model (sum of longest diameters of target lesions, mm) for adults with locally advanced or metastatic biliary tract cancer receiving second-line bintrafusp alfa 1200 mg IV every 2 weeks (phase II INTR@PID BTC 047, NCT03833661). Closed-form profile TS(t) = BTS * exp(-ks * t) + kg * t with t in weeks from the first dose: exponential shrinkage of the baseline tumor plus linear (zero-order) growth. Log-normal IIV on BTS, ks and kg with a full 3 x 3 correlation block; residual error is additive on the log scale (log-transform-both-sides). Covariates: minor allele frequency, time since diagnosis, AST, CRP and Asian race on BTS; bintrafusp alfa dynamic AUC and IDH1 mutation on ks; minor allele frequency, IDH1 mutation, CRP and blood tumor mutational burden at or above the cohort median on kg. The paper prints the coefficients but not the covariate equations or reference values; the categorical effects use the fractional-change form that Figure 3 confirms, while the power form and reference values of the continuous covariates and the median-normalised linear AUC effect are maintainers' assumptions documented in the vignette."
  reference <- paste(
    "Milenkovic-Grisic A-M, Terranova N, Mould DR, Vugmeyster Y, Mrowiec T,",
    "Machl A, Girard P, Venkatakrishnan K, Khandelwal A. Tumor growth inhibition",
    "modeling in patients with second line biliary tract cancer and first line",
    "non-small cell lung cancer based on bintrafusp alfa trials.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(1):143-153.",
    "doi:10.1002/psp4.13068",
    sep = " "
  )
  vignette <- "MilenkovicGrisic_2024_tumor_dynamics"
  units <- list(time = "week", dosing = "mg", concentration = "mm (tumor size)")

  covariateData <- list(
    CTDNA_VAF = list(
      description = "Minor (variant) allele frequency of the somatic variants detected in the patient's genomic sample, in per cent, baseline (time-fixed). Paper wording: 'minor allele frequency for all detected somatic variants in the sample'.",
      units = "percent",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters BTS and kg as power effects (CTDNA_VAF / 1)^theta. Neither the covariate form, the units nor the reference value is printed in Milenkovic-Grisic 2024 or its supplement; the power form and the 1 percent reference are maintainers' assumptions. Must be > 0 (a power effect with a positive exponent sends BTS and kg to 0 at a frequency of 0). The paper does not say whether the sample was plasma ctDNA or tumor tissue; the companion genomic covariate on kg is a BLOOD tumor mutational burden, which suggests a plasma (liquid-biopsy) panel.",
      source_name = "Minor allele frequency (MAF)"
    ),
    T_DIAG_CANCER = list(
      description = "Time from the primary biliary tract cancer diagnosis to the first bintrafusp alfa dose, in days, baseline (time-fixed).",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters BTS as a power effect (T_DIAG_CANCER / 365)^-0.151. The paper prints neither the form, the units nor the reference value; the power form and the 365-day (1 year) reference are maintainers' assumptions. Must be > 0.",
      source_name = "Time since diagnosis"
    ),
    AST = list(
      description = "Baseline serum aspartate aminotransferase (time-fixed).",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters BTS as a power effect (AST / 30)^0.189. The paper prints neither the form nor the reference value; the power form and the 30 U/L reference are maintainers' assumptions.",
      source_name = "AST"
    ),
    CRP = list(
      description = "Baseline C-reactive protein, standard assay (time-fixed).",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters BTS as (CRP / 10)^0.0777 and kg as (CRP / 10)^0.152. The paper prints neither the form, the units nor the reference value; the power form, mg/L and the 10 mg/L reference are maintainers' assumptions. Must be > 0.",
      source_name = "CRP"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator: 1 = Asian, 0 = non-Asian.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian; Milenkovic-Grisic 2024 Figure 3 caption)",
      notes = "Fractional-change effect on BTS, BTS * (1 - 0.2 * RACE_ASIAN). The form reproduces Figure 3e (Asian curve starts at 0.8 x the reference curve).",
      source_name = "Asian race"
    ),
    TUM_IDH1_MUT = list(
      description = "Tumor IDH1 mutation indicator: 1 = IDH1-mutant tumor, 0 = IDH1 wild type.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no IDH1 mutation; Milenkovic-Grisic 2024 Figure 3 caption)",
      notes = "Fractional-change effects ks * (1 + 2.13 * TUM_IDH1_MUT) and kg * (1 + 1.17 * TUM_IDH1_MUT). The (1 + theta) form reproduces Figure 3f (mutant curve about 52 mm at week 20 and 91 mm at week 65 from a 68 mm start); an exp(theta) form does not.",
      source_name = "IDH1 mutation"
    ),
    TMB_BLOOD_HIGH = list(
      description = "Blood tumor mutational burden at or above the 50th percentile of the analysis population: 1 = at or above the median, 0 = below.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (blood TMB below the 50th percentile; Milenkovic-Grisic 2024 Figure 3 caption)",
      notes = "Fractional-change effect kg * (1 + 0.629 * TMB_BLOOD_HIGH); the form reproduces Figure 3g. The numeric median blood TMB that defines the threshold is not reported.",
      source_name = "Blood tumor mutational burden"
    ),
    AUC_BINTRAFUSP_REL = list(
      description = "Bintrafusp alfa dynamic AUC (individual area under the serum concentration-time curve predicted by the bintrafusp alfa popPK model from the actual dosing history, including interruptions) divided by the analysis-population median, dimensionless; 1 = the median exposure.",
      units = "(ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters ks as the linear effect ks * (1 + 0.339 * (AUC_BINTRAFUSP_REL - 1)). The paper calls the exposure a 'dynamic area under the concentration-time curve' from the popPK model of Vugmeyster 2022 (Cancer Chemother Pharmacol 90:369-379) but prints neither its units, its time window, the effect form nor a reference value. Normalising by the cohort median is the maintainers' choice: the dimensionless coefficient 0.339 is only plausible against a normalised AUC, and the Figure 3 median curve reaches its minimum where the typical ks applies. May be supplied time-varying if the dynamic AUC is updated over time.",
      source_name = "AUC (dynamic area under the curve)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 143L,
    n_studies = 1L,
    age_range = "median 65 years (range 39-83)",
    weight_range = "median 65 kg (range 38-100)",
    sex_female_pct = 40.6,
    disease_state = "Locally advanced or metastatic biliary tract cancer previously treated with first-line platinum-based chemotherapy (second line).",
    dose_range = "Bintrafusp alfa 1200 mg IV every 2 weeks.",
    regions = "Multicenter international phase II single-arm trial INTR@PID BTC 047 (NCT03833661).",
    sampling_window = "Tumor assessments every 6 weeks up to 18 months, then every 12 weeks (Table S1).",
    assay = "Tumor size = RECIST 1.1 sum of longest diameters of target lesions (mm), mainly CT and/or MRI.",
    notes = "143 of 157 treated patients (91.1%) had evaluable PK and tumor-size data; 85 men (59.4%). Hyperprogression (Matos criteria) in 43 patients. Missing covariates were imputed with the median (continuous) or mode (categorical)."
  )

  ini({
    # Typical values: Milenkovic-Grisic 2024 Table 1 (2L BTC, bintrafusp alfa),
    # cross-checked against the Results text for ks and kg.
    lrbase <- log(51.3)  ; label("Typical baseline tumor size BTS (mm)")                  # Table 1: Baseline TS (mm) = 51.3, RSE 11.4
    lkse   <- log(0.0159) ; label("Typical exponential tumor-shrinkage rate constant ks (1/week)") # Table 1: ks (week-1) = 0.0159, RSE 26.7; Results text 0.0159
    lp     <- log(0.624)  ; label("Typical linear tumor-growth rate constant kg (mm/week)")  # Table 1: kg (mm/week) = 0.624, RSE 20; Results text 0.624

    # Covariate effects, Table 1. Power exponents for the continuous covariates
    # and fractional changes for the categorical ones (see covariateData notes).
    e_vaf_rbase    <-  0.106  ; label("Power exponent of minor allele frequency on BTS (unitless)")          # Table 1: Minor allele frequency on baseline TS = 0.106, RSE 26.2
    e_tdiag_rbase  <- -0.151  ; label("Power exponent of time since diagnosis on BTS (unitless)")           # Table 1: Time since diagnosis on baseline TS = -0.151, RSE 35.9
    e_ast_rbase    <-  0.189  ; label("Power exponent of AST on BTS (unitless)")                            # Table 1: AST on baseline TS = 0.189, RSE 54
    e_crp_rbase    <-  0.0777 ; label("Power exponent of CRP on BTS (unitless)")                            # Table 1: CRP on baseline TS = 0.0777, RSE 43.8
    e_asian_rbase  <- -0.2    ; label("Fractional change in BTS for Asian race (unitless)")                 # Table 1: Asian race on baseline TS = -0.2, RSE 35.2
    e_auc_kse      <-  0.339  ; label("Linear slope of median-normalised dynamic AUC on ks (unitless)")     # Table 1: AUC effect = 0.339, RSE 15.9
    e_idh1_kse     <-  2.13   ; label("Fractional change in ks for IDH1 mutation (unitless)")               # Table 1: IDH1 mutation on ks = 2.13, RSE 55.4
    e_vaf_p        <-  0.262  ; label("Power exponent of minor allele frequency on kg (unitless)")          # Table 1: Minor allele frequency on kg = 0.262, RSE 33.2
    e_idh1_p       <-  1.17   ; label("Fractional change in kg for IDH1 mutation (unitless)")               # Table 1: IDH1 mutation on kg = 1.17, RSE 51.5
    e_crp_p        <-  0.152  ; label("Power exponent of CRP on kg (unitless)")                             # Table 1: CRP on kg = 0.152, RSE 55.1
    e_tmbblood_p   <-  0.629  ; label("Fractional change in kg for blood TMB at or above the median (unitless)") # Table 1: Blood tumor mutational burden on kg = 0.629, RSE 70.1

    # IIV: Table 1 reports %CV; omega^2 = log(1 + CV^2). Off-diagonals are
    # corr * sqrt(omega_i^2 * omega_j^2) from the Table 1 correlations.
    #   BTS 50 %CV -> 0.2231; ks 103.4 %CV -> 0.7271; kg 72.7 %CV -> 0.4243
    #   cov(BTS, ks) = -0.362 * sqrt(0.2231 * 0.7271) = -0.1458
    #   cov(BTS, kg) =  0.465 * sqrt(0.2231 * 0.4243) =  0.1431
    #   cov(ks, kg)  = -0.845 * sqrt(0.7271 * 0.4243) = -0.4694
    etalrbase + etalkse + etalp ~ c(0.2231,
                                    -0.1458, 0.7271,
                                    0.1431, -0.4694, 0.4243) # Table 1: IIV BTS 50, ks 103.4, kg 72.7 %CV; Corr(BTS,ks) -0.362, Corr(BTS,kg) 0.465, Corr(ks,kg) -0.845

    # Residual error: log-transform-both-sides additive error (Eq 2), reported
    # as a 15.6 %CV proportional RUV; the log-scale SD is 0.156.
    expSd <- 0.156 ; label("Additive residual error SD on the log scale (unitless)") # Table 1: Proportional residual unexplained variability = 15.6 %CV, RSE 6.2
  })

  model({
    # Covariate multipliers. Continuous covariates: power form with
    # maintainers' reference values (minor allele frequency 1 percent,
    # 365 days since diagnosis, AST 30 U/L, CRP 10 mg/L). Categorical
    # covariates: (1 + theta), as reproduced from Figure 3e-g. AUC: linear
    # in the median-normalised dynamic AUC.
    cov_rbase <- (CTDNA_VAF / 1)^e_vaf_rbase *
      (T_DIAG_CANCER / 365)^e_tdiag_rbase *
      (AST / 30)^e_ast_rbase *
      (CRP / 10)^e_crp_rbase *
      (1 + e_asian_rbase * RACE_ASIAN)
    cov_kse <- (1 + e_auc_kse * (AUC_BINTRAFUSP_REL - 1)) *
      (1 + e_idh1_kse * TUM_IDH1_MUT)
    cov_p <- (CTDNA_VAF / 1)^e_vaf_p *
      (CRP / 10)^e_crp_p *
      (1 + e_idh1_p * TUM_IDH1_MUT) *
      (1 + e_tmbblood_p * TMB_BLOOD_HIGH)

    rbase <- exp(lrbase + etalrbase) * cov_rbase
    kse <- exp(lkse + etalkse) * cov_kse
    p <- exp(lp + etalp) * cov_p

    # Milenkovic-Grisic 2024 Eq 1: TS(t) = BTS * exp(-ks * t) + kg * t, with t
    # in weeks from the first bintrafusp alfa dose.
    tumor_size <- rbase * exp(-kse * t) + p * t
    tumor_size ~ lnorm(expSd)
  })
}
