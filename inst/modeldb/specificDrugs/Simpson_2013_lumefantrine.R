Simpson_2013_lumefantrine <- function() {
  description <- "In vitro (P. falciparum). Sigmoid Emax inhibition model of lumefantrine effect on hypoxanthine uptake by clinical Plasmodium falciparum isolates from the Thai-Myanmar border (Shoklo Malaria Research Unit, 1993-2005), with pfmdr1 genotype covariate effects on EC50. The 'subject' in the NLME framework is a parasite isolate (n=324 isolates with lumefantrine data). STIM_LUMEFANTRINE_NM is the per-well drug concentration in the in vitro hypoxanthine-uptake-inhibition assay; the model has no PK and no time evolution. E0 and Emax are fixed per Simpson 2013 Table 3 footnote."
  reference <- paste(
    "Simpson JA, Jamsen KM, Anderson TJC, Zaloumis S, Nair S, Woodrow C, White NJ, Nosten F, Price RN. (2013).",
    "Nonlinear Mixed-Effects Modelling of In Vitro Drug Susceptibility and Molecular Correlates of",
    "Multidrug Resistant Plasmodium falciparum.",
    "PLoS ONE 8(7):e69505.",
    "doi:10.1371/journal.pone.0069505.",
    sep = " "
  )
  vignette <- "Simpson_2013_lumefantrine"
  units <- list(
    time          = "(unused; each record is one drug-well concentration in a 24+18 h hypoxanthine-uptake-inhibition assay)",
    dosing        = "(no PK dosing; STIM_LUMEFANTRINE_NM is the in vitro applied well concentration in nM)",
    concentration = "(observation is normalised hypoxanthine uptake fraction E, range 0-1)"
  )

  covariateData <- list(
    STIM_LUMEFANTRINE_NM = list(
      description        = "Applied lumefantrine concentration in the in vitro hypoxanthine-uptake-inhibition assay well (nM). Per-record covariate. Doubling-dilution series 2.40 to 235.8 nM, plus drug-free control well (CONC = 0). Drives the sigmoid Emax inhibition of normalised hypoxanthine uptake.",
      units              = "nM",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source paper used a doubling-dilution series from 235.8 nM down to 2.40 nM (Methods, In vitro Drug Assay), plus drug-free controls. New canonical entry registered alongside this model in inst/references/covariate-columns.md.",
      source_name        = "C"
    ),
    PFMDR1_86Y = list(
      description        = "Plasmodium falciparum pfmdr1 codon-86 tyrosine mutant indicator (1 = single-copy pfmdr1 with 86Y mutation, Simpson 2013 Genotype 2; 0 = otherwise). Time-fixed per isolate.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = Simpson 2013 Genotype 1 single-copy WT 86N/1042N (when all four PFMDR1 indicators are 0)",
      notes              = "Parasite-genome indicator (not host pharmacogenetics). Mutually exclusive with the other three PFMDR1 indicators in the Simpson 2013 cohort. Genotype 2 prevalence for lumefantrine: 16 of 324 isolates (4.9%) per Table 3.",
      source_name        = "X1"
    ),
    PFMDR1_1042D = list(
      description        = "Plasmodium falciparum pfmdr1 codon-1042 aspartate mutant indicator (1 = single-copy pfmdr1 with 1042D mutation, Simpson 2013 Genotype 3; 0 = otherwise). Time-fixed per isolate.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = Simpson 2013 Genotype 1 single-copy WT (when all four PFMDR1 indicators are 0)",
      notes              = "Mutually exclusive with the other three PFMDR1 indicators. Genotype 3 prevalence for lumefantrine: 17 of 324 isolates (5.2%) per Table 3.",
      source_name        = "X2"
    ),
    PFMDR1_CN2 = list(
      description        = "Plasmodium falciparum pfmdr1 double-copy amplification indicator (1 = two copies of pfmdr1, all WT 86N/1042N, Simpson 2013 Genotype 4; 0 = otherwise). Time-fixed per isolate.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = Simpson 2013 Genotype 1 single-copy WT (when all four PFMDR1 indicators are 0)",
      notes              = "Mutually exclusive with the other three PFMDR1 indicators. Genotype 4 prevalence for lumefantrine: 83 of 324 isolates (25.6%) per Table 3.",
      source_name        = "X3"
    ),
    PFMDR1_CN3PLUS = list(
      description        = "Plasmodium falciparum pfmdr1 triple-or-more-copy amplification indicator (1 = three or more copies of pfmdr1, all WT 86N/1042N, Simpson 2013 Genotype 5; 0 = otherwise). Time-fixed per isolate.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = Simpson 2013 Genotype 1 single-copy WT (when all four PFMDR1 indicators are 0)",
      notes              = "Mutually exclusive with the other three PFMDR1 indicators. Genotype 5 prevalence for lumefantrine: 25 of 324 isolates (7.7%) per Table 3.",
      source_name        = "X4"
    )
  )

  population <- list(
    species        = "in vitro (Plasmodium falciparum clinical isolates from the Shoklo Malaria Research Unit, western Thai-Myanmar border, 1993-2005)",
    n_subjects     = 324,
    n_studies      = 1,
    age_range      = "(not applicable; parasite isolates)",
    weight_range   = "(not applicable; parasite isolates)",
    sex_female_pct = NA_real_,
    disease_state  = "acute P. falciparum malaria",
    dose_range     = "in vitro doubling-dilution series 2.40 to 235.8 nM lumefantrine, plus drug-free control well (Methods, In vitro Drug Assay)",
    regions        = "western Thai-Myanmar border (Shoklo Malaria Research Unit clinics, Mae Sot, Thailand)",
    n_isolates_lumefantrine = 324,
    pfmdr1_distribution = "Genotype 1 (single-copy WT 86N/1042N): 183 of 324 (56.5%); Genotype 2 (single-copy 86Y): 16 of 324 (4.9%); Genotype 3 (single-copy 1042D): 17 of 324 (5.2%); Genotype 4 (double-copy WT): 83 of 324 (25.6%); Genotype 5 (triple+ copy WT): 25 of 324 (7.7%). Table 3.",
    observations_per_isolate = "median 22 observations per isolate across all drugs, range 7 to 44 (Results paragraph 1).",
    notes = "n=490 total isolates available across the four-drug study, with lumefantrine data on 324 isolates (Results paragraph 1; Table 3). Lumefantrine had the highest STS-method exclusion rate among the four drugs (24.4% of isolates had CV > 15%, Table 2), reflecting its relatively narrow dynamic range across the doubling-dilution series. Drug-susceptibility assay = hypoxanthine uptake inhibition (Methods, In vitro Drug Assay; previously described in reference [13])."
  )

  ini({
    # ---------------------------------------------------------------------
    # Sigmoid Emax inhibition (Methods, Eq. 1):
    #   E = Emax - (Emax - E0) * C^gamma / (C^gamma + EC50^gamma)
    # Emax and E0 fixed per Table 3 footnote (#Emax fixed to 0.98;
    # #E0 fixed to 0.01).
    # ---------------------------------------------------------------------
    e0   <- fixed(0.01); label("Minimum normalised hypoxanthine uptake at maximal drug inhibition (fraction, unitless)") # Table 3 footnote (#E0 fixed to 0.01)
    emax <- fixed(0.98); label("Maximum normalised hypoxanthine uptake at zero drug (fraction, unitless)") # Table 3 footnote (#Emax fixed to 0.98)

    # ---------------------------------------------------------------------
    # Population EC50 and slope (Genotype 1 reference). EC50 from Table 3
    # row "Lumefantrine" Genotype 1 column: 35.7 nM (95% CI 31.4, 39.9).
    # Slope from Table 1 NLME row lumefantrine: 2.73 (95% reference range
    # 1.22, 6.10). Slope covariates dropped per skill sidecar.
    # ---------------------------------------------------------------------
    lec50  <- log(35.7); label("Population log-EC50 for the WT reference parasite (log of nM)") # Table 3 Lumefantrine Genotype 1 reference
    lgamma <- log(2.73); label("Population log-slope gamma of the sigmoid inhibition curve (log of unitless)") # Table 1 NLME row lumefantrine

    # ---------------------------------------------------------------------
    # pfmdr1 genotype covariate effects on EC50 (Table 3 percent-change row
    # converted to fractions).
    # ---------------------------------------------------------------------
    e_pfmdr1_86y_ec50     <- -0.31;  label("Proportional effect of single-copy 86Y mutant on EC50 (fraction)")  # Table 3 Lumefantrine Genotype 2 percent change -31 (-62, 0)
    e_pfmdr1_1042d_ec50   <- -0.57;  label("Proportional effect of single-copy 1042D mutant on EC50 (fraction)") # Table 3 Lumefantrine Genotype 3 percent change -57 (-76, -37)
    e_pfmdr1_cn2_ec50     <-  0.82;  label("Proportional effect of double-copy WT on EC50 (fraction)")           # Table 3 Lumefantrine Genotype 4 percent change 82 (46, 119)
    e_pfmdr1_cn3plus_ec50 <-  0.75;  label("Proportional effect of triple-or-more-copy WT on EC50 (fraction)")   # Table 3 Lumefantrine Genotype 5 percent change 75 (28, 122)

    # ---------------------------------------------------------------------
    # Between-isolate variances.
    # ---------------------------------------------------------------------
    etalec50  ~ 0.63;    # Table 3 footnote: Between-isolate variance for EC50 lumefantrine = 0.63 (SE 0.050)
    etalgamma ~ 0.1681;  # Table 1 NLME lumefantrine slope SD = 0.41 log_e units; variance = 0.41^2 = 0.1681

    # ---------------------------------------------------------------------
    # Residual error.
    # ---------------------------------------------------------------------
    propSd <- sqrt(0.019);  label("Proportional residual SD on normalised uptake fraction E (unitless)") # Table 3 footnote: proportional variance 0.019 (SE 0.0020) lumefantrine; SD = sqrt(0.019) = 0.1378
    addSd  <- sqrt(0.0009); label("Additive residual SD on normalised uptake fraction E (unitless)")     # Table 3 footnote: additive variance 0.0009 (SE 0.0002) lumefantrine; SD = sqrt(0.0009) = 0.0300
  })

  model({
    # ---------------------------------------------------------------------
    # Per-record applied lumefantrine concentration (nM); 0 for the
    # drug-free control well.
    # ---------------------------------------------------------------------
    conc <- STIM_LUMEFANTRINE_NM

    # ---------------------------------------------------------------------
    # Genotype multiplier on EC50.
    # ---------------------------------------------------------------------
    ec50_geno <- 1 + e_pfmdr1_86y_ec50     * PFMDR1_86Y +
                     e_pfmdr1_1042d_ec50   * PFMDR1_1042D +
                     e_pfmdr1_cn2_ec50     * PFMDR1_CN2 +
                     e_pfmdr1_cn3plus_ec50 * PFMDR1_CN3PLUS

    # ---------------------------------------------------------------------
    # Individual structural parameters.
    # ---------------------------------------------------------------------
    ec50  <- exp(lec50 + etalec50) * ec50_geno
    gamma <- exp(lgamma + etalgamma)

    # ---------------------------------------------------------------------
    # Sigmoid Emax inhibition (Methods Eq. 1).
    # ---------------------------------------------------------------------
    effect <- emax - (emax - e0) * conc^gamma / (conc^gamma + ec50^gamma)

    # ---------------------------------------------------------------------
    # Observation model.
    # ---------------------------------------------------------------------
    effect ~ add(addSd) + prop(propSd)
  })
}
