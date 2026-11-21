MilenkovicGrisic_2024_pembrolizumab_nsclc <- function() {
  description <- "Empirical tumor-size model (sum of longest diameters of target lesions, mm) for adults with previously untreated advanced non-small cell lung cancer with high PD-L1 expression (at least 80 percent PD-L1-positive tumor cells) receiving first-line pembrolizumab 200 mg IV every 3 weeks (pembrolizumab control arm of the randomised phase III INTR@PID Lung 037, NCT03631706). Closed-form profile TS(t) = BTS * exp(-ks * t) + kg * t with t in weeks from the first dose: exponential shrinkage of the baseline tumor plus linear (zero-order) growth. Log-normal IIV on BTS, ks and kg with a full 3 x 3 correlation block; residual error is additive on the log scale (log-transform-both-sides). No covariate, and no pembrolizumab exposure effect, was retained in the final model."
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

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 143L,
    n_studies = 1L,
    age_range = "median 67 years (range 36-88)",
    weight_range = "median 66.9 kg (range 40.0-141.5)",
    sex_female_pct = 22.4,
    disease_state = "Previously untreated advanced non-small cell lung cancer with high PD-L1 tumor expression (at least 80 percent PD-L1-positive tumor cells); first line.",
    dose_range = "Pembrolizumab 200 mg IV every 3 weeks.",
    regions = "Multicenter international adaptive phase III randomised controlled trial INTR@PID Lung 037 (NCT03631706), pembrolizumab arm.",
    sampling_window = "Tumor assessments every 6 weeks up to 18 months, then every 12 weeks (Table S1).",
    assay = "Tumor size = RECIST 1.1 sum of longest diameters of target lesions (mm), mainly CT or X-ray.",
    notes = "143 of 152 patients randomised to pembrolizumab (94.1%) had evaluable PK and tumor-size data; 111 men (77.6%). No hyperprogression was identified in this arm."
  )

  ini({
    # Typical values: Milenkovic-Grisic 2024 Table 2, pembrolizumab columns,
    # cross-checked against the Results text for ks and kg.
    lrbase <- log(80.4)   ; label("Typical baseline tumor size BTS (mm)")                         # Table 2 pembrolizumab: Baseline TS (mm) = 80.4, RSE 4.91
    lkse   <- log(0.0307) ; label("Typical exponential tumor-shrinkage rate constant ks (1/week)") # Table 2 pembrolizumab: ks (week-1) = 0.0307, RSE 9.06; Results text 0.0307
    lp     <- log(0.643)  ; label("Typical linear tumor-growth rate constant kg (mm/week)")        # Table 2 pembrolizumab: kg (mm/week) = 0.643, RSE 12.9; Results text 0.643

    # IIV: Table 2 reports %CV; omega^2 = log(1 + CV^2). Off-diagonals are
    # corr * sqrt(omega_i^2 * omega_j^2) from the Table 2 correlations.
    #   BTS 55.5 %CV -> 0.2685; ks 70.1 %CV -> 0.3997; kg 92.7 %CV -> 0.6202
    #   cov(BTS, ks) = -0.0228 * sqrt(0.2685 * 0.3997) = -0.00747
    #   cov(BTS, kg) =  0.288  * sqrt(0.2685 * 0.6202) =  0.1175
    #   cov(ks, kg)  = -0.0834 * sqrt(0.3997 * 0.6202) = -0.04153
    etalrbase + etalkse + etalp ~ c(0.2685,
                                    -0.00747, 0.3997,
                                    0.1175, -0.04153, 0.6202) # Table 2 pembrolizumab: IIV BTS 55.5, ks 70.1, kg 92.7 %CV; Corr(BTS,ks) -0.0228, Corr(BTS,kg) 0.288, Corr(ks,kg) -0.0834

    expSd <- 0.113 ; label("Additive residual error SD on the log scale (unitless)") # Table 2 pembrolizumab: Proportional residual unexplained variability = 11.3 %CV, RSE 2.42
  })

  model({
    rbase <- exp(lrbase + etalrbase)
    kse <- exp(lkse + etalkse)
    p <- exp(lp + etalp)

    # Milenkovic-Grisic 2024 Eq 1: TS(t) = BTS * exp(-ks * t) + kg * t, with t
    # in weeks from the first pembrolizumab dose.
    tumor_size <- rbase * exp(-kse * t) + p * t
    tumor_size ~ lnorm(expSd)
  })
}
