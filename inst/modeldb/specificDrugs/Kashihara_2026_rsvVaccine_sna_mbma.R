Kashihara_2026_rsvVaccine_sna_mbma <- function() {
  description <- paste0(
    "MBMA. Immune correlate of protection for respiratory syncytial virus ",
    "(RSV) vaccines in older adults: a study-level meta-regression linking ",
    "placebo-corrected serum neutralizing activity (SNA) against RSV subtype ",
    "A to vaccine efficacy (VE) at three clinical severity levels -- RSV ",
    "acute respiratory infection (ARI), RSV lower respiratory tract disease ",
    "with at least 2 clinical symptoms (LRTD 2+), and RSV LRTD with at least ",
    "3 clinical symptoms (LRTD 3+). VE is modelled on a negative-effect-size ",
    "scale, -log(1 - VE/100), as a severity-specific intercept plus a ",
    "severity-specific slope on the log2 vaccine-to-placebo SNA ratio, and is ",
    "back-transformed to percent. Fitted by weighted linear mixed effects ",
    "(R nlme) to 19 published VE values from seven randomised ",
    "placebo-controlled phase 2b/3 trials of seven different RSV vaccines. ",
    "Variability is BETWEEN-TRIAL (inter-study), encoded as a single ",
    "study-level eta shared by all three severity levels, so the model ",
    "simulates trial-level efficacy and is NOT suitable for ",
    "individual-subject simulation. This is the paper's PRIMARY analysis ",
    "(complete-case, RSV-A SNA); the companion secondary model ",
    "modellib('Kashihara_2026_rsvVaccine_snaCmi_mbma') adds a ",
    "cell-mediated-immunity term."
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
      "SNA_RSVA_RATIO covariate, not through dosing records."
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
        "is the paper's Equation 2, ",
        "SNA_i = log2(SNA_vaccine,i) - log2(SNA_placebo,i). Observed values ",
        "across the seven analysis trials span 1.6 to 13.9 ",
        "(Resolve and VANIR 1.6; D4420C00005 2.5; ConquerRSV 7.7; ",
        "AReSVi-006 10.7; CYPRESS 13.0; RENOIR 13.9), so predictions ",
        "outside roughly 1.5-14 extrapolate beyond the calibration range."
      ),
      source_name        = "RSV-A SNA (ratio to placebo)"
    )
  )

  covariatesDataExcluded <- list(
    STUDY_DURATION = list(
      description = paste0(
        "Duration of the RSV season covered by the efficacy follow-up."
      ),
      units       = "(varies)",
      type        = "continuous",
      notes       = paste0(
        "Kashihara 2026 Discussion: 'The effects of study duration and the ",
        "collection time for SNA were evaluated in the exploratory ",
        "covariate analysis but were not statistically significant.' ",
        "Screened but not retained in the final model, and no point ",
        "estimate is reported, so it is documentation only."
      )
    ),
    SNA_COLLECTION_TIME = list(
      description = paste0(
        "Time after vaccination at which the SNA sample was drawn."
      ),
      units       = "day",
      type        = "continuous",
      notes       = paste0(
        "Kashihara 2026 Discussion; screened in the exploratory covariate ",
        "analysis and not statistically significant. Immunogenicity ",
        "responses assessed approximately 28 days after the single dose ",
        "were used throughout. No point estimate is reported."
      )
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
      "trial-level vaccine efficacy, not an individual participant. 19 VE ",
      "values across the seven trials entered the primary complete-case ",
      "analysis (three severity levels each for AReSVi-006, ConquerRSV, ",
      "CYPRESS, RENOIR and VANIR; two each for D4420C00005, which reports ",
      "no LRTD 3+ value, and Resolve, which reports no LRTD 2+ value). ",
      "n_subjects is the sum of the enrolled participants in Kashihara ",
      "2026 Table 1 (24966 + 35541 + 5782 + 1900 + 34284 + 11856 + 18348). ",
      "All trials were placebo-controlled and covered a single RSV season. ",
      "Immunogenicity endpoints were taken from the efficacy trial itself ",
      "or cross-matched from an earlier-phase trial of the same sponsor ",
      "with a similar design, population and intervention."
    )
  )

  # ==========================================================================
  # ini(): Kashihara 2026 Table 3, upper block ('SNA-VE model'), columns
  # headed 'Primary analysis based on CCA approach'. These are the RSV-A SNA
  # complete-case estimates, which the paper reports as its primary result
  # (Figure 1) and uses for its headline statement that a SNA ratio of 8
  # corresponds to about 70% VE against RSV-LRTD 3+.
  #
  # The columns headed 'Sensitivity analysis based on MI approach' with
  # RSV-B SNA are a sensitivity analysis, not a final model, and are
  # therefore NOT extracted (see the vignette's Assumptions and deviations
  # section, and references/replicate-author-structure.md).
  #
  # Scale: every parameter below lives on the paper's negative-effect-size
  # scale, -log(1 - VE/100), and NONE is log-transformed. The intercepts are
  # signed and small; the slopes are positive. Reported standard errors and
  # 95% confidence intervals are carried in the trailing comments.
  # ==========================================================================
  ini({
    # ----- Severity-specific intercepts -------------------------------------
    # The transformed VE when the vaccine and placebo arms have equal SNA
    # titers (SNA ratio = 1, so log2 ratio = 0). Kashihara 2026 Results 3.2:
    # 'Intercepts represented the VE when SNA levels were equal between
    # vaccinated and placebo groups. In this case, similar infection rates in
    # both groups would result in theoretical VE close to zero, which
    # explains the lack of statistical significance of the intercept
    # parameters.' All three 95% CIs span zero, as expected.
    tve_ref_ari <- 0.047
    label("Transformed VE against RSV-ARI at equal vaccine and placebo SNA titers (negative effect size, unitless)")  # Kashihara 2026 Table 3, row 'RSV-ARI intercept', primary CCA: 0.047, SE 0.296, 95% CI [-0.531, 0.626]

    tve_ref_lrtd2 <- -0.121
    label("Transformed VE against RSV-LRTD 2+ at equal vaccine and placebo SNA titers (negative effect size, unitless)")  # Kashihara 2026 Table 3, row 'RSV-LRTD2+ intercept', primary CCA: -0.121, SE 0.296, 95% CI [-0.700, 0.458]

    tve_ref_lrtd3 <- -0.293
    label("Transformed VE against RSV-LRTD 3+ at equal vaccine and placebo SNA titers (negative effect size, unitless)")  # Kashihara 2026 Table 3, row 'RSV-LRTD3+ intercept', primary CCA: -0.293, SE 0.297, 95% CI [-0.874, 0.287]

    # ----- Severity-specific SNA slopes -------------------------------------
    # Kashihara 2026 Results 3.2: 'The SNA slopes were statistically
    # significant.' All three 95% CIs exclude zero. The LRTD 3+ slope is the
    # steepest, which is the paper's central finding that protection against
    # the most severe endpoint rises fastest with neutralizing titer.
    e_sna_tve_ari <- 0.323
    label("Change in transformed VE against RSV-ARI per doubling of the vaccine-to-placebo RSV-A SNA ratio (negative effect size per log2 ratio)")  # Kashihara 2026 Table 3, row 'SNA slope on RSV-ARI', primary CCA: 0.323, SE 0.110, 95% CI [0.109, 0.537]

    e_sna_tve_lrtd2 <- 0.285
    label("Change in transformed VE against RSV-LRTD 2+ per doubling of the vaccine-to-placebo RSV-A SNA ratio (negative effect size per log2 ratio)")  # Kashihara 2026 Table 3, row 'SNA slope on RSV-LRTD2+', primary CCA: 0.285, SE 0.110, 95% CI [0.070, 0.501]

    e_sna_tve_lrtd3 <- 0.517
    label("Change in transformed VE against RSV-LRTD 3+ per doubling of the vaccine-to-placebo RSV-A SNA ratio (negative effect size per log2 ratio)")  # Kashihara 2026 Table 3, row 'SNA slope on RSV-LRTD3+', primary CCA: 0.517, SE 0.117, 95% CI [0.289, 0.745]

    # ----- Between-trial variability ----------------------------------------
    # Kashihara 2026 Equation 1: 'eta_i is between-trial variability in study
    # i estimated with an additive normal distribution.' The eta carries the
    # subscript i only, NOT ik, so ONE draw per trial is shared by all three
    # severity levels of that trial; it is added on the transformed scale.
    # Table 3 reports it as a standard deviation, so the variance below is
    # 0.315^2. Named eta_study_* rather than the popPK eta<param> convention
    # because it is BETWEEN-TRIAL, not between-subject.
    eta_study_tve ~ 0.099225   # Kashihara 2026 Table 3, row 'Between-trial variability (SD)', primary CCA: 0.315; variance = 0.315^2

    # ----- Residual error ---------------------------------------------------
    # Kashihara 2026 Equation 1 defines a residual eps_ik that was 'weighted
    # by the SE of VE at clinical severity level k in the active arm j of
    # study i', fitted in nlme. Neither the residual scale parameter nor the
    # per-observation weights are reported anywhere in the paper or its
    # Table 3, so no published value exists to encode. The tiny fixed
    # additive residuals below are placeholders that exist only so rxode2 has
    # an error model to attach to each typical-value output; they are NOT
    # published quantities. See the vignette Assumptions and deviations
    # section.
    addSd_ve_ari <- fixed(0.001)
    label("Placeholder additive residual SD on VE against RSV-ARI (percent); NOT a published value")  # not from source -- the paper's SE-weighted residual scale is unreported; see vignette Assumptions and deviations

    addSd_ve_lrtd2 <- fixed(0.001)
    label("Placeholder additive residual SD on VE against RSV-LRTD 2+ (percent); NOT a published value")  # not from source -- the paper's SE-weighted residual scale is unreported; see vignette Assumptions and deviations

    addSd_ve_lrtd3 <- fixed(0.001)
    label("Placeholder additive residual SD on VE against RSV-LRTD 3+ (percent); NOT a published value")  # not from source -- the paper's SE-weighted residual scale is unreported; see vignette Assumptions and deviations
  })

  model({
    # ======================================================================
    # 1. Placebo-corrected SNA (Kashihara 2026 Equation 2)
    #
    #   SNA_i = log2(SNA_vaccine,i) - log2(SNA_placebo,i) = log2(ratio)
    #
    # SNA_RSVA_RATIO is the ratio itself, exactly as tabulated in Table 2, so
    # the difference of logs collapses to the log of the tabulated ratio.
    # ln2 is spelled out because rxode2's parser provides log() but not
    # log2().
    # ======================================================================
    ln2 <- 0.693147180559945
    sna <- log(SNA_RSVA_RATIO) / ln2

    # ======================================================================
    # 2. Transformed VE by severity level (Kashihara 2026 Equation 1)
    #
    #   -log(1 - VE_ik/100%) = Intercept_k + Slope_k * SNA_i + eta_i + eps_ik
    #
    # The single between-trial eta is shared across all three severity
    # levels of a trial, per the subscript on eta_i in Equation 1.
    # ======================================================================
    tve_ari   <- tve_ref_ari   + e_sna_tve_ari   * sna + eta_study_tve
    tve_lrtd2 <- tve_ref_lrtd2 + e_sna_tve_lrtd2 * sna + eta_study_tve
    tve_lrtd3 <- tve_ref_lrtd3 + e_sna_tve_lrtd3 * sna + eta_study_tve

    # ======================================================================
    # 3. Back-transformation to percent (Kashihara 2026 Methods 2.2.1)
    #
    #   VE = 100 * (1 - exp[-(Intercept + Slope * SNA)])
    #
    # VE is bounded above by 100% and is unbounded below: a negative
    # transformed VE gives a negative efficacy, which is what the paper
    # observed for the two low-titer trials (D4420C00005 and Resolve).
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
