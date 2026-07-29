Olagunju_2018_efavirenz <- function() {
  description <- "One-compartment population PK model for oral efavirenz in HIV-positive pregnant women (Olagunju 2018), with composite CYP2B6 516G>T (rs3745274) and 983T>C (rs28399499) metaboliser status (slow / intermediate / fast) as a categorical covariate on CL/F and fixed-exponent allometric body-weight scaling on CL/F and V/F."
  reference <- paste(
    "Olagunju A, Schipani A, Bolaji O, Khoo S, Owen A.",
    "Evaluation of universal versus genotype-guided efavirenz dose reduction",
    "in pregnant women using population pharmacokinetic modeling.",
    "J Antimicrob Chemother. 2018;73(1):165-172. doi:10.1093/jac/dkx334."
  )
  vignette <- "Olagunju_2018_efavirenz"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline) per the Olagunju 2018 sparse + intensive PK design. Drives fixed-exponent allometric scaling on CL/F (exponent 0.75) and V/F (exponent 1.0) with standard reference weight 70 kg per Olagunju 2018 Methods 'Population Pharmacokinetic-Pharmacogenetic Model Development' paragraph 3.",
      source_name        = "WT"
    ),
    SNP_CYP2B6_RS3745274_T_COUNT = list(
      description        = "Count of CYP2B6 c.516G>T (rs3745274, p.Q172H) T-alleles per subject (0/1/2). 0 = GG homozygous wild-type, 1 = GT heterozygous, 2 = TT homozygous variant.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (germline genotype). Olagunju 2018 Methods 'Sample Collection, SNP Genotyping and Efavirenz Quantification' paragraph 1 -- this SNP is combined with rs28399499 to define a composite CYP2B6 metaboliser status (fast / intermediate / slow) by counting the total number of variant alleles across the two SNPs. Cohort allele-genotype frequencies (n = 77; Olagunju 2018 Table 1): GG 0.32, GT 0.54, TT 0.14.",
      source_name        = "CYP2B6 516G>T (rs3745274)"
    ),
    SNP_CYP2B6_RS28399499_C_COUNT = list(
      description        = "Count of CYP2B6 c.983T>C (rs28399499, p.I328T) C-alleles per subject (0/1/2). 0 = TT homozygous wild-type, 1 = TC heterozygous, 2 = CC homozygous variant.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (germline genotype). Olagunju 2018 Methods 'Sample Collection, SNP Genotyping and Efavirenz Quantification' paragraph 1 -- combined with rs3745274 to define a composite CYP2B6 metaboliser status. Cohort allele-genotype frequencies (n = 77; Olagunju 2018 Table 1): TT 0.75, TC 0.25, CC 0.00. The Olagunju 2018 cohort had no 983CC homozygotes, so the 'ultra-slow' metaboliser substratum is not represented in the fitted data; the packaged model still classifies a hypothetical 983CC subject as slow (each variant allele contributes to the composite count).",
      source_name        = "CYP2B6 983T>C (rs28399499)"
    )
  )

  population <- list(
    species               = "human",
    n_subjects            = 77,
    n_studies             = 1,
    age_range             = "18-39 years",
    age_median            = "27 years",
    weight_range          = "48-83 kg",
    weight_median         = "57 kg",
    sex_female_pct        = 100,
    race_ethnicity        = c(Nigerian = 100),
    disease_state         = "HIV-positive pregnant women receiving combination antiretroviral therapy containing 600 mg efavirenz daily plus two nucleoside reverse transcriptase inhibitors for >= 4 weeks. Patients on anti-tuberculosis drugs or other co-medications with known / uncertain interactions with antiretrovirals were excluded.",
    dose_range            = "600 mg orally once daily (evening dose). The packaged model and validation vignette also simulate 200 mg and 400 mg daily reduced-dose scenarios that the source paper investigates by Monte Carlo simulation in pregnant women stratified by CYP2B6 metaboliser status.",
    regions               = "Nigeria (3 hospitals in Benue State: Bishop Murray Medical Centre, Makurdi; St Monica's Hospital, Adikpo; St Mary's Hospital, Okpoga)",
    gestational_age_range = "11-36 weeks (median 28); trimester distribution 5% first / 25% second / 70% third (Olagunju 2018 Table 1)",
    cyp2b6_freq           = "516G>T (rs3745274): GG 0.32, GT 0.54, TT 0.14. 983T>C (rs28399499): TT 0.75, TC 0.25, CC 0.00 (Olagunju 2018 Table 1). The cohort had no 983CC homozygotes, so the predicted-very-slow / 'ultra-slow' metaboliser substratum is not represented.",
    notes                 = "ClinicalTrials.gov ID NCT02269462. Data set is the same one analysed in Olagunju et al. 2015 (Clin Pharmacol Ther 97:298-306, ref 11 of the present paper). 252 plasma efavirenz concentrations were available (77 sparse PK samples from 77 women + 175 intensive PK samples from 25 women drawn 0.5-24 h after dose, stratified by genotype). Age, body weight, gestational age, and CYP2B6 516G>T / 983T>C were tested as covariates under stepwise backward elimination; only the CYP2B6 metaboliser status was statistically significant (Methods paragraph 4 + Results 'Population Pharmacokinetic Analysis')."
  )

  ini({
    # ---- Structural PK parameters (Olagunju 2018 Table 2 final population model) ----
    lka <- log(0.61) ; label("First-order absorption rate constant ka (1/h)")              # Olagunju 2018 Table 2 final ka = 0.61 h^-1 (RSE 23%; 90% CI 0.3-0.9)
    lcl <- log(18.0) ; label("CL/F (L/h) -- CYP2B6 fast-metaboliser reference at WT = 70 kg") # Olagunju 2018 Table 2 CL/F_Fast = 18 L/h (RSE 9%; 90% CI 15-21.5); fast = no CYP2B6 variant alleles across rs3745274 + rs28399499 (i.e., 516GG and 983TT)
    lvc <- log(281)  ; label("V/F (L) at WT = 70 kg")                                       # Olagunju 2018 Table 2 V/F = 281 L (RSE 10%; 90% CI 241-320)

    # ---- Composite-metaboliser-status covariate effects on CL/F ----
    # Olagunju 2018 Table 2 reports CL/F directly for each of the three composite
    # metaboliser groups (fast / intermediate / slow). The fast group is the natural
    # reference (no variant alleles); the non-reference groups are encoded as
    # log-ratio multiplicative effects so that the single etalcl IIV applies
    # uniformly on the log-CL scale across all groups.
    e_intermed_cl <- log(16.1 / 18.0) ; label("Log-ratio of intermediate-metaboliser CL/F vs fast-metaboliser reference (unitless)") # Olagunju 2018 Table 2 CL/F_Intermediate = 16.1 L/h (RSE 7%; 90% CI 15.1-18); log(16.1/18) = -0.1117
    e_slow_cl     <- log(6.24 / 18.0) ; label("Log-ratio of slow-metaboliser CL/F vs fast-metaboliser reference (unitless)")         # Olagunju 2018 Table 2 CL/F_Slow = 6.24 L/h (RSE 11%; 90% CI 4.8-7.5); log(6.24/18) = -1.0593

    # ---- Fixed allometric exponents on body weight (paper Methods paragraph 3) ----
    # "An allometric weight model for clearance parameters is given by CLwt = (WT/WTstd)^0.75
    #  and for volume parameters is given Vwt = (WT/WTstd)^1, where WTstd = 70 kg."
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of (WT/70) on CL/F (unitless; fixed)") # Olagunju 2018 Methods paragraph 3
    e_wt_vc <- fixed(1.0)  ; label("Allometric exponent of (WT/70) on V/F (unitless; fixed)")  # Olagunju 2018 Methods paragraph 3

    # ---- IIV (diagonal omega; exponential errors, log-normal per equation 1) ----
    # Olagunju 2018 Methods equation 1: theta_i = theta_1 * exp(eta_i)
    # CV-to-variance conversion: omega^2 = log(CV^2 + 1)
    etalka ~ 0.56167   # 86.8% CV (Olagunju 2018 Table 2 RSE 21%, 90% CI 23-107); log(1 + 0.868^2) = 0.56167
    etalcl ~ 0.15614   # 41.1% CV (Olagunju 2018 Table 2 RSE 16%, 90% CI 36-43);  log(1 + 0.411^2) = 0.15614
    etalvc ~ 0.067379  # 26.4% CV (Olagunju 2018 Table 2 RSE 33%, 90% CI 10-40);  log(1 + 0.264^2) = 0.06738

    # ---- Residual error (proportional only; additive structure did not improve fit) ----
    # Olagunju 2018 Table 2: proportional residual variance = 0.085 (RSE 21%, 90% CI 0.06-0.11).
    # Per NONMEM convention $SIGMA reports the variance; the proportional SD on the linear
    # concentration scale is sqrt(0.085) = 0.292 (~29.2% CV). Consistent with the same
    # reporting style used by the co-author group in Schipani 2011 (variance 0.0085 -> SD 0.092).
    propSd <- sqrt(0.085) ; label("Proportional residual error (fraction)") # Olagunju 2018 Table 2: variance 0.085; propSd = sqrt(0.085) = 0.292
  })

  model({
    # 1. Derive composite CYP2B6 metaboliser status from the two SNP allele counts.
    #    Per Olagunju 2018 Methods 'Sample Collection ...' paragraph 1:
    #      Slow         = >= 2 variant alleles across rs3745274 + rs28399499
    #                     (e.g., 516GT + 983TC, 516TT + 983TT, 516TT + 983TC)
    #      Intermediate = exactly 1 variant allele
    #                     (e.g., 516GT + 983TT, 516GG + 983TC)
    #      Fast         = 0 variant alleles (516GG + 983TT)
    n_variant   <- SNP_CYP2B6_RS3745274_T_COUNT + SNP_CYP2B6_RS28399499_C_COUNT
    is_intermed <- (n_variant == 1)
    is_slow     <- (n_variant >= 2)

    # 2. Typical log-CL/F: fast-metaboliser reference + log-ratio shifts for the two
    #    non-reference groups (mutually exclusive indicators).
    ltvcl <- lcl + e_intermed_cl * is_intermed + e_slow_cl * is_slow

    # 3. Individual PK parameters with fixed-exponent allometric weight scaling
    #    (reference WT = 70 kg).
    ka <- exp(lka + etalka)
    cl <- exp(ltvcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc)   * (WT / 70)^e_wt_vc

    # 4. Micro-constants
    kel <- cl / vc

    # 5. ODE system: one-compartment model with first-order absorption (oral)
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 6. Observation and error
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
