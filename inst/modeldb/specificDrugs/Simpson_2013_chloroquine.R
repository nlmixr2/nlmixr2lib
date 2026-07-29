Simpson_2013_chloroquine <- function() {
  description <- "In vitro (P. falciparum). Sigmoid Emax inhibition model of chloroquine effect on hypoxanthine uptake by clinical Plasmodium falciparum isolates from the Thai-Myanmar border (Shoklo Malaria Research Unit, 1993-2005), with pfmdr1 genotype covariate effects on EC50. The 'subject' in the NLME framework is a parasite isolate (n=421 isolates with chloroquine data). STIM_CHLOROQUINE_NM is the per-well drug concentration in the in vitro hypoxanthine-uptake-inhibition assay; the model has no PK and no time evolution. E0 and Emax are fixed per Simpson 2013 Table 3 footnote."
  reference <- paste(
    "Simpson JA, Jamsen KM, Anderson TJC, Zaloumis S, Nair S, Woodrow C, White NJ, Nosten F, Price RN. (2013).",
    "Nonlinear Mixed-Effects Modelling of In Vitro Drug Susceptibility and Molecular Correlates of",
    "Multidrug Resistant Plasmodium falciparum.",
    "PLoS ONE 8(7):e69505.",
    "doi:10.1371/journal.pone.0069505.",
    sep = " "
  )
  vignette <- "Simpson_2013_chloroquine"
  units <- list(
    time          = "(unused; each record is one drug-well concentration in a 24+18 h hypoxanthine-uptake-inhibition assay)",
    dosing        = "(no PK dosing; STIM_CHLOROQUINE_NM is the in vitro applied well concentration in nM)",
    concentration = "(observation is normalised hypoxanthine uptake fraction E, range 0-1)"
  )

  covariateData <- list(
    STIM_CHLOROQUINE_NM = list(
      description        = "Applied chloroquine concentration in the in vitro hypoxanthine-uptake-inhibition assay well (nM). Per-record covariate. Doubling-dilution series 10.02 to 10255.9 nM, plus drug-free control well (CONC = 0). Drives the sigmoid Emax inhibition of normalised hypoxanthine uptake.",
      units              = "nM",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source paper used a doubling-dilution series from 10255.9 nM down to 10.02 nM (Methods, In vitro Drug Assay), plus drug-free controls. New canonical entry registered alongside this model in inst/references/covariate-columns.md.",
      source_name        = "C"
    ),
    PFMDR1_86Y = list(
      description        = "Plasmodium falciparum pfmdr1 codon-86 tyrosine mutant indicator (1 = single-copy pfmdr1 with 86Y mutation, Simpson 2013 Genotype 2; 0 = otherwise). Time-fixed per isolate.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = Simpson 2013 Genotype 1 single-copy WT 86N/1042N (when all four PFMDR1 indicators are 0)",
      notes              = "Parasite-genome indicator (not host pharmacogenetics). Mutually exclusive with PFMDR1_1042D, PFMDR1_CN2, and PFMDR1_CN3PLUS in the Simpson 2013 cohort (Thai pfmdr1 mutations occur almost exclusively on single-copy parasites; amplifications occur exclusively on WT 86N/1042N parasites). Genotype 2 prevalence in Simpson 2013: 5% of 490 isolates (20 chloroquine isolates per Table 3).",
      source_name        = "X1"
    ),
    PFMDR1_1042D = list(
      description        = "Plasmodium falciparum pfmdr1 codon-1042 aspartate mutant indicator (1 = single-copy pfmdr1 with 1042D mutation, Simpson 2013 Genotype 3; 0 = otherwise). Time-fixed per isolate.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = Simpson 2013 Genotype 1 single-copy WT (when all four PFMDR1 indicators are 0)",
      notes              = "Mutually exclusive with the other three PFMDR1 indicators. Genotype 3 prevalence in Simpson 2013: 5% of 490 isolates (19 chloroquine isolates per Table 3).",
      source_name        = "X2"
    ),
    PFMDR1_CN2 = list(
      description        = "Plasmodium falciparum pfmdr1 double-copy amplification indicator (1 = two copies of pfmdr1, all WT 86N/1042N, Simpson 2013 Genotype 4; 0 = otherwise). Time-fixed per isolate.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = Simpson 2013 Genotype 1 single-copy WT (when all four PFMDR1 indicators are 0)",
      notes              = "Mutually exclusive with the other three PFMDR1 indicators. Genotype 4 prevalence in Simpson 2013: 26% of 490 isolates (113 chloroquine isolates per Table 3).",
      source_name        = "X3"
    ),
    PFMDR1_CN3PLUS = list(
      description        = "Plasmodium falciparum pfmdr1 triple-or-more-copy amplification indicator (1 = three or more copies of pfmdr1, all WT 86N/1042N, Simpson 2013 Genotype 5; 0 = otherwise). Time-fixed per isolate.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = Simpson 2013 Genotype 1 single-copy WT (when all four PFMDR1 indicators are 0)",
      notes              = "Mutually exclusive with the other three PFMDR1 indicators. Genotype 5 prevalence in Simpson 2013: 15% of 490 isolates (57 chloroquine isolates per Table 3).",
      source_name        = "X4"
    )
  )

  population <- list(
    species        = "in vitro (Plasmodium falciparum clinical isolates from the Shoklo Malaria Research Unit, western Thai-Myanmar border, 1993-2005)",
    n_subjects     = 421,
    n_studies      = 1,
    age_range      = "(not applicable; parasite isolates)",
    weight_range   = "(not applicable; parasite isolates)",
    sex_female_pct = NA_real_,
    disease_state  = "acute P. falciparum malaria",
    dose_range     = "in vitro doubling-dilution series 10.02 to 10255.9 nM chloroquine, plus drug-free control well (Methods, In vitro Drug Assay)",
    regions        = "western Thai-Myanmar border (Shoklo Malaria Research Unit clinics, Mae Sot, Thailand)",
    n_isolates_chloroquine = 421,
    pfmdr1_distribution = "Genotype 1 (single-copy WT 86N/1042N): 212 of 421 (50.4%); Genotype 2 (single-copy 86Y): 20 of 421 (4.8%); Genotype 3 (single-copy 1042D): 19 of 421 (4.5%); Genotype 4 (double-copy WT): 113 of 421 (26.8%); Genotype 5 (triple+ copy WT): 57 of 421 (13.5%). Table 3.",
    observations_per_isolate = "median 22 observations per isolate across all drugs, range 7 to 44 (Results paragraph 1).",
    notes = "n=490 total isolates available across the four-drug study, with chloroquine data on 421 isolates (Results paragraph 1; Table 3). Drug-susceptibility assay = hypoxanthine uptake inhibition (Methods, In vitro Drug Assay): fresh isolates adjusted to 0.5-1.0% infected RBC at 1.5% haematocrit in RPMI-1640 + 10% heat-inactivated AB sera, dispensed into 96-well plates with duplicate serial dilutions, incubated 24 h, pulsed with [3H] hypoxanthine, incubated a further 18 h, then harvested. The Methods section reports the assay had been previously described in reference [13] of the paper."
  )

  ini({
    # ---------------------------------------------------------------------
    # Sigmoid Emax inhibition (Methods, Eq. 1):
    #   E = Emax - (Emax - E0) * C^gamma / (C^gamma + EC50^gamma)
    # Emax and E0 are fixed per Table 3 footnote ("#Emax fixed to 0.98;
    # #E0 fixed to 0.01"). The # symbol in Table 3 marks the chloroquine
    # row; both values are reused for the other three drugs in the
    # Simpson 2013 study because the paper's structural equation and the
    # normalised-uptake range are common across all four drugs.
    # ---------------------------------------------------------------------
    e0   <- fixed(0.01); label("Minimum normalised hypoxanthine uptake at maximal drug inhibition (fraction, unitless)") # Table 3 footnote (#E0 fixed to 0.01)
    emax <- fixed(0.98); label("Maximum normalised hypoxanthine uptake at zero drug (fraction, unitless)") # Table 3 footnote (#Emax fixed to 0.98)

    # ---------------------------------------------------------------------
    # Population EC50 and slope. Reference category (Genotype 1: single-copy
    # WT 86N/1042N) population estimate from Table 3 row "Chloroquine #"
    # column "Genotype 1 Single Copy WT": 242 nM (95% CI 223, 260). Slope
    # gamma population estimate from Table 1 NLME row chloroquine: 4.14
    # (95% reference range 1.85 to 9.25). The slope-covariate effects
    # (theta_5 to theta_8) are dropped (File S2 not on disk; main text:
    # "minimal changes in the slope of the concentration-effect profile
    # with the 1042D variant were observed" and "the slope of the
    # effect-concentration curve did not differ significantly for these
    # molecular comparisons"); documented in vignette Errata.
    # ---------------------------------------------------------------------
    lec50  <- log(242);  label("Population log-EC50 for the WT reference parasite (log of nM)")  # Table 3 Chloroquine Genotype 1 reference
    lgamma <- log(4.14); label("Population log-slope gamma of the sigmoid inhibition curve (log of unitless)") # Table 1 NLME row chloroquine

    # ---------------------------------------------------------------------
    # pfmdr1 genotype covariate effects on EC50 (Methods Eq. 2 modified;
    # values from Table 3 "Percent change" row, converted to fractions).
    # Mutually exclusive in the source cohort; at most one indicator = 1
    # per isolate. The multiplier `(1 + theta_1*X_1 + ... + theta_4*X_4)`
    # collapses to a per-genotype constant when applied.
    # ---------------------------------------------------------------------
    e_pfmdr1_86y_ec50     <-  0.44;  label("Proportional effect of single-copy 86Y mutant on EC50 (fraction)")  # Table 3 Chloroquine Genotype 2 percent change 44 (14, 73)
    e_pfmdr1_1042d_ec50   <-  0.48;  label("Proportional effect of single-copy 1042D mutant on EC50 (fraction)") # Table 3 Chloroquine Genotype 3 percent change 48 (-4, 100)
    e_pfmdr1_cn2_ec50     <- -0.10;  label("Proportional effect of double-copy WT on EC50 (fraction)")           # Table 3 Chloroquine Genotype 4 percent change -10 (-23, 3)
    e_pfmdr1_cn3plus_ec50 <- -0.10;  label("Proportional effect of triple-or-more-copy WT on EC50 (fraction)")   # Table 3 Chloroquine Genotype 5 percent change -10 (-28, 7)

    # ---------------------------------------------------------------------
    # Between-isolate (IIV) variances. EC50 variance is the Table 3
    # footnote value (after removing genotype-explained variance). Gamma
    # variance is from Table 1 NLME row chloroquine SD column (log_e
    # units) squared because Table 3 does not report a separate slope
    # variance (with the slope-covariate effects dropped, the gamma IIV
    # is the same as in the no-covariate base model).
    # ---------------------------------------------------------------------
    etalec50  ~ 0.39;     # Table 3 footnote: Between-isolate variance for EC50 chloroquine = 0.39 (SE 0.026)
    etalgamma ~ 0.1681;   # Table 1 NLME chloroquine slope SD = 0.41 log_e units; variance = 0.41^2 = 0.1681

    # ---------------------------------------------------------------------
    # Residual error (combined additive + proportional on the normalised
    # uptake fraction; Methods Eq. 3). Reported variances in Table 3
    # footnote are converted to SDs (sqrt(variance)) for the nlmixr2
    # propSd / addSd parameterisation.
    # ---------------------------------------------------------------------
    propSd <- sqrt(0.013); label("Proportional residual SD on normalised uptake fraction E (unitless)")  # Table 3 footnote: proportional variance 0.013 (SE 0.0018) chloroquine; SD = sqrt(0.013) = 0.114
    addSd  <- sqrt(0.001); label("Additive residual SD on normalised uptake fraction E (unitless)")       # Table 3 footnote: additive variance 0.001 (SE 0.0002) chloroquine; SD = sqrt(0.001) = 0.0316
  })

  model({
    # ---------------------------------------------------------------------
    # Per-record applied chloroquine concentration (nM). For the drug-free
    # control well, the data assembler sets STIM_CHLOROQUINE_NM = 0; with
    # gamma > 0 the term C^gamma evaluates to 0 and effect collapses to
    # Emax (no inhibition), matching the Methods Eq. 1 boundary.
    # ---------------------------------------------------------------------
    conc <- STIM_CHLOROQUINE_NM

    # ---------------------------------------------------------------------
    # Genotype multiplier on EC50 (Methods Eq. 2 with theta_1..theta_4).
    # Encoded as the proportional-shift form `1 + sum(theta * X)`. With
    # X_1..X_4 mutually exclusive in the Simpson 2013 cohort, this reduces
    # to a single per-genotype multiplicative factor per isolate.
    # ---------------------------------------------------------------------
    ec50_geno <- 1 + e_pfmdr1_86y_ec50     * PFMDR1_86Y +
                     e_pfmdr1_1042d_ec50   * PFMDR1_1042D +
                     e_pfmdr1_cn2_ec50     * PFMDR1_CN2 +
                     e_pfmdr1_cn3plus_ec50 * PFMDR1_CN3PLUS

    # ---------------------------------------------------------------------
    # Individual structural parameters (multiplicative log-normal IIV).
    # ---------------------------------------------------------------------
    ec50  <- exp(lec50 + etalec50) * ec50_geno
    gamma <- exp(lgamma + etalgamma)

    # ---------------------------------------------------------------------
    # Sigmoid Emax inhibition (Methods Eq. 1, paper's authoring orientation
    # E = Emax - (Emax - E0) * C^gamma / (C^gamma + EC50^gamma)). At C = 0
    # the second term is 0 and E = Emax; at C >> EC50 the second term
    # approaches Emax - E0 and E -> E0. Maximum inhibition occurs at E = 0
    # corresponding to no parasite growth.
    # ---------------------------------------------------------------------
    effect <- emax - (emax - e0) * conc^gamma / (conc^gamma + ec50^gamma)

    # ---------------------------------------------------------------------
    # Observation model: combined additive + proportional residual error
    # on the normalised hypoxanthine-uptake fraction (Methods Eq. 3).
    # ---------------------------------------------------------------------
    effect ~ add(addSd) + prop(propSd)
  })
}
