Simpson_2013_mefloquine <- function() {
  description <- "In vitro (P. falciparum). Sigmoid Emax inhibition model of mefloquine effect on hypoxanthine uptake by clinical Plasmodium falciparum isolates from the Thai-Myanmar border (Shoklo Malaria Research Unit, 1993-2005), with pfmdr1 genotype covariate effects on EC50. The 'subject' in the NLME framework is a parasite isolate (n=460 isolates with mefloquine data). STIM_MEFLOQUINE_NM is the per-well drug concentration in the in vitro hypoxanthine-uptake-inhibition assay; the model has no PK and no time evolution. E0 and Emax are fixed per Simpson 2013 Table 3 footnote."
  reference <- paste(
    "Simpson JA, Jamsen KM, Anderson TJC, Zaloumis S, Nair S, Woodrow C, White NJ, Nosten F, Price RN. (2013).",
    "Nonlinear Mixed-Effects Modelling of In Vitro Drug Susceptibility and Molecular Correlates of",
    "Multidrug Resistant Plasmodium falciparum.",
    "PLoS ONE 8(7):e69505.",
    "doi:10.1371/journal.pone.0069505.",
    sep = " "
  )
  vignette <- "Simpson_2013_mefloquine"
  units <- list(
    time          = "(unused; each record is one drug-well concentration in a 24+18 h hypoxanthine-uptake-inhibition assay)",
    dosing        = "(no PK dosing; STIM_MEFLOQUINE_NM is the in vitro applied well concentration in nM)",
    concentration = "(observation is normalised hypoxanthine uptake fraction E, range 0-1)"
  )

  covariateData <- list(
    STIM_MEFLOQUINE_NM = list(
      description        = "Applied mefloquine concentration in the in vitro hypoxanthine-uptake-inhibition assay well (nM). Per-record covariate. Doubling-dilution series 1.62 to 1646.6 nM, plus drug-free control well (CONC = 0). Drives the sigmoid Emax inhibition of normalised hypoxanthine uptake.",
      units              = "nM",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source paper used a doubling-dilution series from 1646.6 nM down to 1.62 nM (Methods, In vitro Drug Assay), plus drug-free controls. New canonical entry registered alongside this model in inst/references/covariate-columns.md.",
      source_name        = "C"
    ),
    PFMDR1_86Y = list(
      description        = "Plasmodium falciparum pfmdr1 codon-86 tyrosine mutant indicator (1 = single-copy pfmdr1 with 86Y mutation, Simpson 2013 Genotype 2; 0 = otherwise). Time-fixed per isolate.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = Simpson 2013 Genotype 1 single-copy WT 86N/1042N (when all four PFMDR1 indicators are 0)",
      notes              = "Parasite-genome indicator (not host pharmacogenetics). Mutually exclusive with the other three PFMDR1 indicators in the Simpson 2013 cohort. Genotype 2 prevalence for mefloquine: 25 of 460 isolates (5.4%) per Table 3.",
      source_name        = "X1"
    ),
    PFMDR1_1042D = list(
      description        = "Plasmodium falciparum pfmdr1 codon-1042 aspartate mutant indicator (1 = single-copy pfmdr1 with 1042D mutation, Simpson 2013 Genotype 3; 0 = otherwise). Time-fixed per isolate.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = Simpson 2013 Genotype 1 single-copy WT (when all four PFMDR1 indicators are 0)",
      notes              = "Mutually exclusive with the other three PFMDR1 indicators. Genotype 3 prevalence for mefloquine: 24 of 460 isolates (5.2%) per Table 3.",
      source_name        = "X2"
    ),
    PFMDR1_CN2 = list(
      description        = "Plasmodium falciparum pfmdr1 double-copy amplification indicator (1 = two copies of pfmdr1, all WT 86N/1042N, Simpson 2013 Genotype 4; 0 = otherwise). Time-fixed per isolate.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = Simpson 2013 Genotype 1 single-copy WT (when all four PFMDR1 indicators are 0)",
      notes              = "Mutually exclusive with the other three PFMDR1 indicators. Genotype 4 prevalence for mefloquine: 118 of 460 isolates (25.7%) per Table 3.",
      source_name        = "X3"
    ),
    PFMDR1_CN3PLUS = list(
      description        = "Plasmodium falciparum pfmdr1 triple-or-more-copy amplification indicator (1 = three or more copies of pfmdr1, all WT 86N/1042N, Simpson 2013 Genotype 5; 0 = otherwise). Time-fixed per isolate.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = Simpson 2013 Genotype 1 single-copy WT (when all four PFMDR1 indicators are 0)",
      notes              = "Mutually exclusive with the other three PFMDR1 indicators. Genotype 5 prevalence for mefloquine: 63 of 460 isolates (13.7%) per Table 3.",
      source_name        = "X4"
    )
  )

  population <- list(
    species        = "in vitro (Plasmodium falciparum clinical isolates from the Shoklo Malaria Research Unit, western Thai-Myanmar border, 1993-2005)",
    n_subjects     = 460,
    n_studies      = 1,
    age_range      = "(not applicable; parasite isolates)",
    weight_range   = "(not applicable; parasite isolates)",
    sex_female_pct = NA_real_,
    disease_state  = "acute P. falciparum malaria",
    dose_range     = "in vitro doubling-dilution series 1.62 to 1646.6 nM mefloquine, plus drug-free control well (Methods, In vitro Drug Assay)",
    regions        = "western Thai-Myanmar border (Shoklo Malaria Research Unit clinics, Mae Sot, Thailand)",
    n_isolates_mefloquine = 460,
    pfmdr1_distribution = "Genotype 1 (single-copy WT 86N/1042N): 230 of 460 (50.0%); Genotype 2 (single-copy 86Y): 25 of 460 (5.4%); Genotype 3 (single-copy 1042D): 24 of 460 (5.2%); Genotype 4 (double-copy WT): 118 of 460 (25.7%); Genotype 5 (triple+ copy WT): 63 of 460 (13.7%). Table 3.",
    observations_per_isolate = "median 22 observations per isolate across all drugs, range 7 to 44 (Results paragraph 1).",
    notes = "n=490 total isolates available across the four-drug study, with mefloquine data on 460 isolates (Results paragraph 1; Table 3). Mefloquine is the drug for which pfmdr1 amplification effects are most pronounced: double-copy parasites have 139% higher EC50 and triple-or-more-copy parasites have 188% higher EC50 than single-copy WT parasites (Table 3 percent-change row). Drug-susceptibility assay = hypoxanthine uptake inhibition (Methods, In vitro Drug Assay; previously described in reference [13])."
  )

  ini({
    # ---------------------------------------------------------------------
    # Sigmoid Emax inhibition (Methods, Eq. 1):
    #   E = Emax - (Emax - E0) * C^gamma / (C^gamma + EC50^gamma)
    # Emax and E0 fixed per Table 3 footnote (#Emax fixed to 0.98;
    # #E0 fixed to 0.01); the # symbol in Table 3 marks the chloroquine
    # row but the same fixed values apply to all four drugs because the
    # normalised-uptake response variable is bounded in [0, 1] across
    # the study.
    # ---------------------------------------------------------------------
    e0   <- fixed(0.01); label("Minimum normalised hypoxanthine uptake at maximal drug inhibition (fraction, unitless)") # Table 3 footnote (#E0 fixed to 0.01)
    emax <- fixed(0.98); label("Maximum normalised hypoxanthine uptake at zero drug (fraction, unitless)") # Table 3 footnote (#Emax fixed to 0.98)

    # ---------------------------------------------------------------------
    # Population EC50 and slope (reference category Genotype 1: single-copy
    # WT 86N/1042N). EC50 from Table 3 row "Mefloquine" Genotype 1 column:
    # 53.0 nM (95% CI 48.0, 58.1). Slope from Table 1 NLME row mefloquine:
    # 3.10 (95% reference range 1.39 to 6.92). Slope-covariate effects
    # (theta_5..theta_8) dropped per skill sidecar (File S2 not on disk;
    # main text describes slope effects as "minimal" / "not significant").
    # ---------------------------------------------------------------------
    lec50  <- log(53.0); label("Population log-EC50 for the WT reference parasite (log of nM)") # Table 3 Mefloquine Genotype 1 reference
    lgamma <- log(3.10); label("Population log-slope gamma of the sigmoid inhibition curve (log of unitless)") # Table 1 NLME row mefloquine

    # ---------------------------------------------------------------------
    # pfmdr1 genotype covariate effects on EC50 (Methods Eq. 2 modified;
    # values from Table 3 "Percent change" row, converted to fractions).
    # ---------------------------------------------------------------------
    e_pfmdr1_86y_ec50     <- -0.59;  label("Proportional effect of single-copy 86Y mutant on EC50 (fraction)")  # Table 3 Mefloquine Genotype 2 percent change -59 (-72, -46)
    e_pfmdr1_1042d_ec50   <- -0.42;  label("Proportional effect of single-copy 1042D mutant on EC50 (fraction)") # Table 3 Mefloquine Genotype 3 percent change -42 (-67, -17)
    e_pfmdr1_cn2_ec50     <-  1.39;  label("Proportional effect of double-copy WT on EC50 (fraction)")           # Table 3 Mefloquine Genotype 4 percent change 139 (102, 175)
    e_pfmdr1_cn3plus_ec50 <-  1.88;  label("Proportional effect of triple-or-more-copy WT on EC50 (fraction)")   # Table 3 Mefloquine Genotype 5 percent change 188 (126, 250)

    # ---------------------------------------------------------------------
    # Between-isolate variances. EC50 variance from Table 3 footnote (the
    # residual variance after the genotype effects are removed). Gamma
    # variance from Table 1 NLME row mefloquine SD column (no separate
    # slope variance reported in Table 3 because slope covariates were
    # dropped from this extraction).
    # ---------------------------------------------------------------------
    etalec50  ~ 0.56;    # Table 3 footnote: Between-isolate variance for EC50 mefloquine = 0.56 (SE 0.043)
    etalgamma ~ 0.1681;  # Table 1 NLME mefloquine slope SD = 0.41 log_e units; variance = 0.41^2 = 0.1681

    # ---------------------------------------------------------------------
    # Residual error (combined additive + proportional on the normalised
    # uptake fraction; Methods Eq. 3). Table 3 footnote variances
    # converted to SDs.
    # ---------------------------------------------------------------------
    propSd <- sqrt(0.010); label("Proportional residual SD on normalised uptake fraction E (unitless)") # Table 3 footnote: proportional variance 0.010 (SE 0.0011) mefloquine; SD = sqrt(0.010) = 0.100
    addSd  <- sqrt(0.001); label("Additive residual SD on normalised uptake fraction E (unitless)")     # Table 3 footnote: additive variance 0.001 (SE 0.0001) mefloquine; SD = sqrt(0.001) = 0.0316
  })

  model({
    # ---------------------------------------------------------------------
    # Per-record applied mefloquine concentration (nM). Set to 0 for the
    # drug-free control well; with gamma > 0 the C^gamma term is 0 and
    # effect collapses to Emax (no inhibition).
    # ---------------------------------------------------------------------
    conc <- STIM_MEFLOQUINE_NM

    # ---------------------------------------------------------------------
    # Genotype multiplier on EC50 (Methods Eq. 2 with theta_1..theta_4).
    # Mutually-exclusive indicator setup reduces to a per-genotype constant
    # per isolate.
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
    # Sigmoid Emax inhibition (Methods Eq. 1).
    # ---------------------------------------------------------------------
    effect <- emax - (emax - e0) * conc^gamma / (conc^gamma + ec50^gamma)

    # ---------------------------------------------------------------------
    # Observation model: combined additive + proportional residual error
    # on the normalised uptake fraction (Methods Eq. 3).
    # ---------------------------------------------------------------------
    effect ~ add(addSd) + prop(propSd)
  })
}
