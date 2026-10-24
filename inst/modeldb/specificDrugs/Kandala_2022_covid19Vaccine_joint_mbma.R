Kandala_2022_covid19Vaccine_joint_mbma <- function() {
  description <- paste0(
    "MBMA (joint preclinical + clinical). COVID-19 vaccine immune correlate ",
    "of protection: the joint rhesus-macaque (RM) + clinical MBMA, the second ",
    "of three efficacy-scale models in Kandala 2022. RM viral load is ",
    "transformed to an efficacy scale (0-1) and modelled jointly with clinical ",
    "vaccine efficacy by a single sigmoid (Hill) function with SPECIES as a ",
    "covariate: the RM arm is driven by the log2 SN titre and the human arm ",
    "by the SN/HCS ratio, sharing the maximum efficacy (VEmax, fixed to 1) ",
    "while estimating species-specific Hill coefficients and 50%-protective ",
    "titres. The 50% protective SN titre in humans is 24% of the geometric-",
    "mean human convalescent titre (Figure 4). Fitted by generalized ",
    "non-linear least squares with trial-size weights; the species-specific ",
    "random effect was estimated near zero and dropped, so this is a ",
    "deterministic typical-value model and is NOT suitable for individual-",
    "level simulation. The two sibling efficacy models are ",
    "modellib('Kandala_2022_covid19Vaccine_calibrated_mbma') and ",
    "modellib('Kandala_2022_covid19Vaccine_clinical_mbma')."
  )

  reference <- paste(
    "Kandala B, Plock N, Chawla A, Largajolli A, Robey S, Watson K, et al.",
    "Accelerating model-informed decisions for COVID-19 vaccine candidates",
    "using a model-based meta-analysis approach.",
    "EBioMedicine. 2022;84:104264. doi:10.1016/j.ebiom.2022.104264.",
    sep = " "
  )

  vignette <- "Kandala_2022_covid19Vaccine"

  units <- list(
    time = paste0(
      "not applicable -- static study-level meta-analysis with no time course."
    ),
    dosing = paste0(
      "not applicable -- vaccine exposure enters through the SN_SARSCOV2_LOG2 ",
      "(RM arm) and SN_HCS_RATIO (human arm) covariates, selected by ",
      "SPECIES_MACAQUE, not through dosing records."
    ),
    concentration = paste0(
      "ve is the vaccine efficacy in percent (RM arm: viral-load-reduction ",
      "efficacy; human arm: clinical vaccine efficacy). The output Cc is ",
      "unused."
    )
  )

  covariateData <- list(
    SPECIES_MACAQUE = list(
      description = paste0(
        "Species indicator: 1 = rhesus-macaque (RM) datum, 0 = human datum. ",
        "Selects the RM viral-load-reduction arm versus the human clinical-",
        "efficacy arm of the joint sigmoid."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (human)",
      notes = paste0(
        "Kandala 2022 Equation 11 indicator I[Species = RM] / ",
        "I[Species = Human]. The register's SPECIES_MACAQUE canonical is ",
        "cynomolgus in its founding example; here it marks the rhesus-macaque ",
        "arm (both are the non-human-primate arm of a cross-species ",
        "translation). SPECIES_MACAQUE = 1 selects the log2 SN / IT50 RM ",
        "term; = 0 selects the SN/HCS human term."
      ),
      source_name = "Species"
    ),
    SN_SARSCOV2_LOG2 = list(
      description = paste0(
        "Study-arm log2 serum neutralizing titre against SARS-CoV-2 (RM arm ",
        "predictor), closest to and prior to viral challenge."
      ),
      units = "(log2 titre)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Kandala 2022 Equation 11, RM term log2[SN_RM]. Used only when ",
        "SPECIES_MACAQUE = 1."
      ),
      source_name = "log2[SN_RM]"
    ),
    SPECIMEN_NASAL = list(
      description = paste0(
        "Upper-vs-lower respiratory specimen indicator for the RM arm: 1 = ",
        "nasal / nasopharyngeal, 0 = BAL fluid / lung tissue. Selects the ",
        "specimen-specific RM 50%-protective titre."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (BAL fluid / lung tissue)",
      notes = paste0(
        "Kandala 2022 Figure 4 parameters log2[IT50]_RM,BAL vs ",
        "log2[IT50]_RM,nasal. SPECIMEN_NASAL = 1 selects it50_rm_nasal ",
        "(8.57); = 0 selects it50_rm_bal (4.92) log2 units. Applies only to ",
        "the RM arm (SPECIES_MACAQUE = 1)."
      ),
      source_name = "Specimen"
    ),
    SN_HCS_RATIO = list(
      description = paste0(
        "Study-arm ratio of the vaccine SN titre to the same-assay geometric-",
        "mean human convalescent serum (HCS) titre (human arm predictor)."
      ),
      units = "(fold; dimensionless ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Kandala 2022 Equation 11, human term SN_h/HCS. Used only when ",
        "SPECIES_MACAQUE = 0."
      ),
      source_name = "SN/HCS"
    )
  )

  population <- list(
    species = "human + rhesus macaque (Macaca mulatta)",
    n_studies = 19L,
    disease_state = paste0(
      "RM arm: SARS-CoV-2 challenge, protection = viral-load reduction. ",
      "Human arm: PCR-confirmed symptomatic COVID-19 of any severity, ",
      "efficacy = relative incidence-rate reduction versus placebo."
    ),
    dose_range = paste0(
      "13 RM challenge studies plus six clinical vaccine candidates ",
      "(BNT162b2, mRNA-1273, Ad26.COV2.S, 2019nCoV, Gam-COVID-Vac, ChAdOx1 ",
      "nCoV-19)."
    ),
    regions = "International (published non-clinical and clinical trials)",
    notes = paste0(
      "Model-based meta-analysis jointly fitting RM and clinical study-arm ",
      "summaries. RM viral loads were transformed to a 0-1 efficacy scale ",
      "before fitting. VE observations were weighted by trial size (clinical ",
      "weights 0.129-1, RM weights 0.1-1). The joint model's root-mean-",
      "squared prediction error was 6.6%. A species-specific random effect ",
      "was tested, estimated near zero, and dropped for a worse AIC."
    )
  )

  # ==========================================================================
  # ini(): Kandala 2022 Figure 4 parameter table (joint MBMA). VEmax and the
  # clinical Hill coefficient were fixed to 1 during model development (each
  # fixing reduced AIC by 2 units; see Results). Reported 95% confidence
  # intervals are carried in the trailing comments. No random effects or
  # residual SD are published (the species random effect was dropped), so the
  # residual below is a placeholder.
  # ==========================================================================
  ini({
    vemax <- fixed(1)
    label("Maximum vaccine efficacy, shared across species (fraction)")  # Kandala 2022 Figure 4 table, row 'VEmax': 1 (fixed)

    hill_rm <- 2.23
    label("Hill coefficient of the RM viral-load-reduction efficacy curve (unitless)")  # Kandala 2022 Figure 4 table, row 'gamma_RM': 2.23, 95% CI [1.50, 3.33]

    hill_clin <- fixed(1)
    label("Hill coefficient of the human clinical-efficacy curve (unitless)")  # Kandala 2022 Figure 4 table, row 'gamma_clin': 1 (fixed; estimated close to 1 but imprecise)

    it50_rm_bal <- 4.92
    label("Log2 SN titre for 50% RM efficacy in BAL fluid / lung tissue (log2 titre)")  # Kandala 2022 Figure 4 table, row 'log2[IT50]_RM,BAL': 4.92, 95% CI [4.06, 5.78]

    it50_rm_nasal <- 8.57
    label("Log2 SN titre for 50% RM efficacy in nasal tissue (log2 titre)")  # Kandala 2022 Figure 4 table, row 'log2[IT50]_RM,nasal': 8.57, 95% CI [6.86, 10.30]

    snhcs50 <- 0.24
    label("SN/HCS titre ratio for 50% clinical vaccine efficacy in humans (dimensionless)")  # Kandala 2022 Figure 4 table, row '[SN/HCS]_50': 0.24, 95% CI [0.07, 0.77]

    addSd <- fixed(0.001)
    label("Placeholder additive residual SD on VE (percent); NOT a published value")  # not from source -- joint-model residual scale is unreported; see vignette Assumptions and deviations
  })

  model({
    # ======================================================================
    # RM-arm specimen-specific 50% protective titre (Figure 4).
    # ======================================================================
    it50_rm <- it50_rm_bal * (1 - SPECIMEN_NASAL) +
      it50_rm_nasal * SPECIMEN_NASAL

    # ======================================================================
    # Species-covariate sigmoid (Kandala 2022 Equation 11). The RM term is
    # driven by the log2 SN titre and its specimen IT50; the human term by the
    # SN/HCS ratio and [SN/HCS]_50. SPECIES_MACAQUE selects the arm.
    # ======================================================================
    frac_rm <- SN_SARSCOV2_LOG2^hill_rm /
      (SN_SARSCOV2_LOG2^hill_rm + it50_rm^hill_rm)
    frac_h <- SN_HCS_RATIO^hill_clin /
      (SN_HCS_RATIO^hill_clin + snhcs50^hill_clin)

    frac <- frac_rm * SPECIES_MACAQUE + frac_h * (1 - SPECIES_MACAQUE)

    # Efficacy as a percentage (VEmax is a 0-1 fraction; the figure plots %).
    ve <- 100 * vemax * frac

    ve ~ add(addSd)
  })
}
