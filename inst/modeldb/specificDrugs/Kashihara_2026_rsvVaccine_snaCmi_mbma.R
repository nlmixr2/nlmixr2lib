Kashihara_2026_rsvVaccine_snaCmi_mbma <- function() {
  description <- paste0(
    "MBMA. Immune correlate of protection for respiratory syncytial virus ",
    "(RSV) vaccines in older adults, extended with cell-mediated immunity: ",
    "a study-level meta-regression linking placebo-corrected serum ",
    "neutralizing activity (SNA) against RSV subtype A AND the ",
    "placebo-corrected interferon-gamma response to vaccine efficacy (VE) at ",
    "three clinical severity levels -- RSV acute respiratory infection ",
    "(ARI), RSV lower respiratory tract disease with at least 2 clinical ",
    "symptoms (LRTD 2+), and RSV LRTD with at least 3 clinical symptoms ",
    "(LRTD 3+). VE is modelled on a negative-effect-size scale, ",
    "-log(1 - VE/100), as a severity-specific intercept plus a COMMON slope ",
    "on the log2 vaccine-to-placebo SNA ratio plus an interferon-gamma slope ",
    "that applies to LRTD 3+ ONLY, and is back-transformed to percent. ",
    "Fitted by weighted linear mixed effects (R nlme) with multiple ",
    "imputation of missing immunogenicity data across seven randomised ",
    "placebo-controlled phase 2b/3 trials of seven different RSV vaccines. ",
    "Variability is BETWEEN-TRIAL (inter-study), encoded as a single ",
    "study-level eta shared by all three severity levels, so the model ",
    "simulates trial-level efficacy and is NOT suitable for ",
    "individual-subject simulation. This is the paper's exploratory ",
    "SECONDARY analysis; the primary SNA-only model is ",
    "modellib('Kashihara_2026_rsvVaccine_sna_mbma')."
  )

  reference <- paste(
    "Kashihara Y, Qin L, Shimizu S, Diderichsen PM, Kotsuma M, Yoshihara K.",
    "Establishing Immune Correlates of Protection Against Respiratory",
    "Syncytial Virus Infection to Accelerate Vaccine Development:",
    "A Model-Based Meta-Analysis.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70133.",
    "doi:10.1002/psp4.70133.",
    sep = " "
  )

  vignette <- "Kashihara_2026_rsvVaccine_immune_correlates"

  units <- list(
    time          = paste0(
      "not applicable -- this is a static study-level meta-regression with ",
      "no time course. Immunogenicity is the value measured about 28 days ",
      "after a single dose; VE is accrued over one RSV season (fall to the ",
      "end of spring)."
    ),
    dosing        = paste0(
      "not applicable -- vaccine exposure enters through the ",
      "SNA_RSVA_RATIO and IFNG_RATIO covariates, not through dosing records."
    ),
    concentration = paste0(
      "ve_ari, ve_lrtd2 and ve_lrtd3 are vaccine efficacies in percent. The ",
      "output Cc is unused."
    )
  )

  covariateData <- list(
    SNA_RSVA_RATIO = list(
      description        = paste0(
        "Placebo-corrected serum neutralizing activity against RSV subtype ",
        "A: the ratio of the vaccine-arm SNA titer to the placebo-arm SNA ",
        "titer, on the original (untransformed) scale, measured about 28 ",
        "days after a single dose. Study-arm level, not individual level. A ",
        "value of 1 means equal titers in the vaccine and placebo arms."
      ),
      units              = "(fold; dimensionless ratio)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Kashihara 2026 Table 2 column 'RSV-A SNA', footnote b: ",
        "'Immunogenicity values were presented as the ratio to placebo on ",
        "the original scale.' The model uses the log2 of this ratio, which ",
        "is the paper's Equation 2. Observed values across the seven ",
        "analysis trials span 1.6 to 13.9."
      ),
      source_name        = "RSV-A SNA (ratio to placebo)"
    ),
    IFNG_RATIO = list(
      description        = paste0(
        "Placebo-corrected cell-mediated immunity: the ratio of the ",
        "vaccine-arm interferon-gamma response to the placebo-arm ",
        "interferon-gamma response, on the original (untransformed) scale, ",
        "measured about 28 days after a single dose. Study-arm level, not ",
        "individual level. A value of 1 means an equal response in the ",
        "vaccine and placebo arms, which zeroes the interferon-gamma term."
      ),
      units              = "(fold; dimensionless ratio)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Kashihara 2026 Table 2 column 'IFN-gamma', footnote b: ",
        "'Immunogenicity values were presented as the ratio to placebo on ",
        "the original scale.' The model uses the NATURAL log of this ratio ",
        "-- Equation 4, ",
        "CMI_i = ln(IFN-gamma_vaccine,i) - ln(IFN-gamma_placebo,i) -- ",
        "which is a DIFFERENT log base from the log2 used for SNA in ",
        "Equation 2. Only four trials reported or could be cross-matched ",
        "for interferon-gamma (VANIR 3.5, D4420C00005 9.0, CYPRESS 12.8, ",
        "RENOIR 20.6); the remaining three were multiply imputed. The ",
        "paper's Figure 3 simulates at 3.5 and 20.6 as the observed range, ",
        "so predictions outside roughly 3.5-20.6 extrapolate beyond the ",
        "calibration range. CD4+ T cell counts were the other ",
        "cell-mediated-immunity marker collected but were excluded from ",
        "the analysis because only two trials had data."
      ),
      source_name        = "IFN-gamma (ratio to placebo)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 132677L,
    n_studies      = 7L,
    age_range      = paste0(
      "Older adults. Six trials enrolled adults aged 60 years and older; ",
      "CYPRESS enrolled adults aged 65 years and older (Kashihara 2026 ",
      "Table 1)."
    ),
    disease_state  = paste0(
      "Community-dwelling older adults at risk of seasonal RSV infection; ",
      "efficacy endpoints are RSV-ARI, RSV-LRTD 2+ and RSV-LRTD 3+."
    ),
    dose_range     = paste0(
      "One dose of vaccine or placebo per participant. The seven trials ",
      "used seven different vaccines: RSVpreF3 120 ug with AS01E adjuvant ",
      "(AReSVi-006), mRNA-1345 50 ug (ConquerRSV), Ad26.RSV.preF 1e11 viral ",
      "particles with RSV preF 150 ug (CYPRESS), MEDI-7510 120 ug with ",
      "glucopyranosyl lipid adjuvant (D4420C00005), RSVpreF 120 ug ",
      "(RENOIR), RSV F 135 ug (Resolve), and MVA-BN-RSV 3e8 infectious ",
      "units (VANIR)."
    ),
    regions        = "International (published phase 2b and later trials)",
    notes          = paste0(
      "Model-based meta-analysis: the unit of observation is a published ",
      "trial-level vaccine efficacy, not an individual participant. This ",
      "secondary analysis was developed on the 11 VE values from the four ",
      "trials with complete RSV-A SNA and interferon-gamma data (CYPRESS, ",
      "D4420C00005, RENOIR and VANIR), then re-estimated on all seven ",
      "trials using multiple imputation (50 imputed datasets by chained ",
      "equations, R mice, method 'norm', pooled by Rubin's rules) without ",
      "altering the model structure. The multiply-imputed estimates are the ",
      "ones encoded here. n_subjects is the sum of the enrolled ",
      "participants in Kashihara 2026 Table 1 across all seven trials."
    )
  )

  # ==========================================================================
  # ini(): Kashihara 2026 Table 3, lower block ('SNA + CMI-VE model'),
  # columns headed 'RSV-A SNA / MI approach'.
  #
  # WHICH COLUMN: Table 3 reports this model twice, under a complete-case
  # (CCA) and a multiple-imputation (MI) approach. The MI column is the
  # final one and is what this file encodes, because Kashihara 2026 Results
  # 3.3 describes CCA as the development step and MI as the extension to
  # the full dataset: 'Initially, the model development process was started
  # using complete case data for both RSV-A SNA and IFN-gamma. ... In order
  # to account for the possible source of bias represented by the exclusion
  # of studies with incomplete data, the analysis was extended to all
  # available data using MI without altering the model structure.' The
  # paper's Figure 3, which is the published simulation of this model, is
  # likewise 'based on the MI model'; the CCA simulation is relegated to
  # supplementary Figure S9. Per
  # references/replicate-author-structure.md a base-plus-final pair is one
  # file holding the final. The CCA estimates are recorded in the vignette's
  # source-trace table for reference.
  #
  # STRUCTURE (Kashihara 2026 Results 3.3): 'The final model included the
  # three clinical severity-specific intercepts, a common SNA effect to
  # three severity levels, and the IFN-gamma effect on RSV-LRTD 3+.' So
  # there is ONE SNA slope rather than three, and the interferon-gamma slope
  # exists for LRTD 3+ ONLY -- confirmed by the Figure 3 caption: 'For
  # RSV-ARI and RSV-LRTD2+, no IFN-gamma effect on VE was estimated;
  # therefore, the simulated VE is based solely on RSV-A SNA.'
  #
  # Scale: every parameter below lives on the paper's negative-effect-size
  # scale, -log(1 - VE/100), and NONE is log-transformed.
  # ==========================================================================
  ini({
    # ----- Severity-specific intercepts -------------------------------------
    tve_ref_ari <- 0.051
    label("Transformed VE against RSV-ARI at equal vaccine and placebo SNA and interferon-gamma responses (negative effect size, unitless)")  # Kashihara 2026 Table 3, SNA + CMI-VE model, row 'RSV-ARI intercept', RSV-A SNA MI: 0.051, SE 0.247, 95% CI [-0.490, 0.591]

    tve_ref_lrtd2 <- -0.168
    label("Transformed VE against RSV-LRTD 2+ at equal vaccine and placebo SNA and interferon-gamma responses (negative effect size, unitless)")  # Kashihara 2026 Table 3, SNA + CMI-VE model, row 'RSV-LRTD2+ intercept', RSV-A SNA MI: -0.168, SE 0.247, 95% CI [-0.709, 0.373]

    tve_ref_lrtd3 <- -0.298
    label("Transformed VE against RSV-LRTD 3+ at equal vaccine and placebo SNA and interferon-gamma responses (negative effect size, unitless)")  # Kashihara 2026 Table 3, SNA + CMI-VE model, row 'RSV-LRTD3+ intercept', RSV-A SNA MI: -0.298, SE 0.309, 95% CI [-1.02, 0.425]

    # ----- Common SNA slope -------------------------------------------------
    # Shared by all three severity levels; contrast with the primary model,
    # which estimates one slope per level. Kashihara 2026 Results 3.3: 'The
    # common SNA effect in MI (0.321, SE = 0.091) was comparable to the SNA
    # effect on RSV-ARI (0.323, SE = 0.110) and RSV-LRTD2+ (0.285,
    # SE = 0.110) in the primary SNA-VE analysis but smaller than the SNA
    # effect on RSV-LRTD3+ (0.517, SE = 0.117) estimated in the primary
    # analysis. The smaller effect could be explained by part of the
    # contribution to VE for RSV-LRTD3+ being accounted for by the IFN-gamma
    # effect.'
    e_sna_tve <- 0.321
    label("Change in transformed VE per doubling of the vaccine-to-placebo RSV-A SNA ratio, common to all three severity levels (negative effect size per log2 ratio)")  # Kashihara 2026 Table 3, SNA + CMI-VE model, row 'SNA slope', RSV-A SNA MI: 0.321, SE 0.091, 95% CI [0.120, 0.521]

    # ----- Interferon-gamma slope, RSV-LRTD 3+ only -------------------------
    # The reported 95% CI is far wider than the standard error alone implies
    # because it is a Rubin's-rules interval on a t-distribution whose
    # degrees of freedom collapse when the fraction of missing information is
    # large -- only four of the seven trials have observed interferon-gamma
    # (Kashihara 2026 Discussion: 'The small number of observed IFN-gamma
    # data points (n = 4) limited characterization of its distribution').
    # The point estimate is the encoded quantity; the interval is recorded
    # for provenance. The paper is explicit that this analysis is
    # exploratory and that 'no definitive conclusions can be drawn regarding
    # the contribution of CMI to VE'.
    e_ifng_tve_lrtd3 <- 0.111
    label("Change in transformed VE against RSV-LRTD 3+ per unit increase in the natural-log vaccine-to-placebo interferon-gamma ratio (negative effect size per ln ratio)")  # Kashihara 2026 Table 3, SNA + CMI-VE model, row 'IFN-g slope on RSV-LRTD3+', RSV-A SNA MI: 0.111, SE 0.138, 95% CI [-2.00, 2.23]

    # ----- Between-trial variability ----------------------------------------
    # Kashihara 2026 Equation 3 carries the same eta_i as Equation 1: one
    # draw per trial, shared by all three severity levels, additive on the
    # transformed scale. Table 3 reports it as a standard deviation, so the
    # variance below is 0.318^2. Named eta_study_* because it is
    # BETWEEN-TRIAL, not between-subject.
    eta_study_tve ~ 0.101124   # Kashihara 2026 Table 3, SNA + CMI-VE model, row 'Between-trial variability (SD)', RSV-A SNA MI: 0.318; variance = 0.318^2

    # ----- Residual error ---------------------------------------------------
    # As in the primary model, the paper's SE-weighted residual scale is not
    # reported. These tiny fixed additive residuals are placeholders so
    # rxode2 has an error model to attach to each typical-value output; they
    # are NOT published quantities. See the vignette Assumptions and
    # deviations section.
    addSd_ve_ari <- fixed(0.001)
    label("Placeholder additive residual SD on VE against RSV-ARI (percent); NOT a published value")  # not from source -- the paper's SE-weighted residual scale is unreported; see vignette Assumptions and deviations

    addSd_ve_lrtd2 <- fixed(0.001)
    label("Placeholder additive residual SD on VE against RSV-LRTD 2+ (percent); NOT a published value")  # not from source -- the paper's SE-weighted residual scale is unreported; see vignette Assumptions and deviations

    addSd_ve_lrtd3 <- fixed(0.001)
    label("Placeholder additive residual SD on VE against RSV-LRTD 3+ (percent); NOT a published value")  # not from source -- the paper's SE-weighted residual scale is unreported; see vignette Assumptions and deviations
  })

  model({
    # ======================================================================
    # 1. Placebo-corrected immunogenicity
    #
    #   Equation 2: SNA_i = log2(SNA_vaccine,i) - log2(SNA_placebo,i)
    #   Equation 4: CMI_i = ln(IFNg_vaccine,i) - ln(IFNg_placebo,i)
    #
    # Note the DIFFERENT log bases: base 2 for SNA, natural log for CMI.
    # Both covariates are supplied as the ratio itself, exactly as
    # tabulated in Table 2, so each difference of logs collapses to the log
    # of the tabulated ratio. ln2 is spelled out because rxode2's parser
    # provides log() but not log2().
    # ======================================================================
    ln2 <- 0.693147180559945
    sna <- log(SNA_RSVA_RATIO) / ln2
    cmi <- log(IFNG_RATIO)

    # ======================================================================
    # 2. Transformed VE by severity level (Kashihara 2026 Equation 3)
    #
    #   -log(1 - VE_ik/100%) = Intercept_k + Slope_k * SNA_i
    #                          + Slope_IFN,k * CMI_i + eta_i + eps_ik
    #
    # with Slope_k common across k and Slope_IFN,k estimated for LRTD 3+
    # only, so the interferon-gamma term is ABSENT from the ARI and LRTD 2+
    # predictors rather than being present with a zero coefficient.
    # ======================================================================
    tve_ari   <- tve_ref_ari   + e_sna_tve * sna + eta_study_tve
    tve_lrtd2 <- tve_ref_lrtd2 + e_sna_tve * sna + eta_study_tve
    tve_lrtd3 <- tve_ref_lrtd3 + e_sna_tve * sna +
      e_ifng_tve_lrtd3 * cmi + eta_study_tve

    # ======================================================================
    # 3. Back-transformation to percent (Kashihara 2026 Methods 2.2.1)
    #
    #   VE = 100 * (1 - exp[-(Intercept + Slope * SNA)])
    # ======================================================================
    ve_ari   <- 100 * (1 - exp(-tve_ari))
    ve_lrtd2 <- 100 * (1 - exp(-tve_lrtd2))
    ve_lrtd3 <- 100 * (1 - exp(-tve_lrtd3))

    # ======================================================================
    # 4. Observations
    # ======================================================================
    ve_ari   ~ add(addSd_ve_ari)
    ve_lrtd2 ~ add(addSd_ve_lrtd2)
    ve_lrtd3 ~ add(addSd_ve_lrtd3)
  })
}
