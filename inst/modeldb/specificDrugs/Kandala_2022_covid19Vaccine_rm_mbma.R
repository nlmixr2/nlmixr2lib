Kandala_2022_covid19Vaccine_rm_mbma <- function() {
  description <- paste0(
    "MBMA (preclinical, rhesus macaque). COVID-19 vaccine immune correlate ",
    "of protection: a study-level meta-analysis relating the log2 serum ",
    "neutralizing (SN) titre after SARS-CoV-2 vaccination to post-challenge ",
    "peak viral load (log10) in rhesus macaques (RM). Peak viral load falls ",
    "from a maximum (VLmax) toward a specimen-specific minimum (VLmin) as a ",
    "sigmoid (Hill) suppression function of log2 SN titre, with a lower 50%-",
    "protective titre (IT50) in bronchoalveolar-lavage/lung tissue than in ",
    "nasal tissue. Fitted by weighted non-linear mixed effects (R nlme) to ",
    "summary-level data from 13 RM challenge studies spanning 21 vaccine ",
    "candidates. Variability is BETWEEN-STUDY (inter-study / challenge-",
    "inoculum), encoded as a single study-level exponential random effect on ",
    "VLmax, so the model simulates study-arm-level viral load and is NOT ",
    "suitable for individual-animal simulation. This is the paper's first ",
    "model (Table 2); the three efficacy-scale models it feeds are ",
    "modellib('Kandala_2022_covid19Vaccine_calibrated_mbma'), ",
    "modellib('Kandala_2022_covid19Vaccine_joint_mbma') and ",
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
      "not applicable -- this is a static study-level meta-analysis with no ",
      "time course. SN titre is the value closest to and prior to viral ",
      "challenge; viral load is the post-challenge peak."
    ),
    dosing = paste0(
      "not applicable -- vaccine exposure enters through the ",
      "SN_SARSCOV2_LOG2 covariate, not through dosing records."
    ),
    concentration = paste0(
      "log10_viral_load is the log10 post-challenge peak viral load (sgRNA ",
      "copies). The output Cc is unused."
    )
  )

  covariateData <- list(
    SN_SARSCOV2_LOG2 = list(
      description = paste0(
        "Study-arm log2 serum neutralizing titre against SARS-CoV-2, ",
        "measured closest to and prior to viral challenge. The predictor ",
        "driving viral-load suppression."
      ),
      units = "(log2 titre)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Kandala 2022 Equation 3, term log2[SN_i,j]. Enters the sigmoid ",
        "suppression function directly as SN_SARSCOV2_LOG2^g, so the column ",
        "holds the already-log2-transformed titre (per the paper's data ",
        "derivation rules: SN titres were log-transformed base 2 prior to ",
        "analysis; control arms at or below the LLOQ/LLOD were set to 0)."
      ),
      source_name = "log2[SN]"
    ),
    SPECIMEN_NASAL = list(
      description = paste0(
        "Upper-vs-lower respiratory specimen indicator: 1 = nasal / ",
        "nasopharyngeal swab, 0 = bronchoalveolar-lavage fluid or lung ",
        "tissue. Selects the specimen-specific 50%-protective titre."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (BAL fluid / lung tissue)",
      notes = paste0(
        "Kandala 2022 Table 1 column 'Specimen', decomposed to a binary ",
        "upper-vs-lower indicator. The final RM MBMA estimated one IT50 for ",
        "BAL and lung tissue and, separately, IT50_nasal for nasal tissue; ",
        "SPECIMEN_NASAL = 1 selects it50_nasal (8.32) and = 0 selects it50 ",
        "(5.84) log2 units (Table 2)."
      ),
      source_name = "Specimen"
    )
  )

  population <- list(
    species = "rhesus macaque (Macaca mulatta)",
    n_studies = 13L,
    disease_state = paste0(
      "SARS-CoV-2 challenge in vaccinated (or monoclonal-antibody-treated, ",
      "or re-challenged) rhesus macaques; protection is reduction in ",
      "post-challenge peak viral load."
    ),
    dose_range = paste0(
      "13 challenge studies spanning 21 SARS-CoV-2 vaccine candidates across ",
      "RNA, DNA, adenovirus-vectored, subunit and VSV-recombinant platforms ",
      "(Kandala 2022 Table 1)."
    ),
    regions = "International (published non-clinical challenge studies)",
    notes = paste0(
      "Model-based meta-analysis: the unit of observation is a study-arm ",
      "summary (geometric-mean log2 SN titre and geometric-mean log10 peak ",
      "viral load), not an individual animal. Post-challenge peak ",
      "sub-genomic RNA (sgRNA) viral load in BAL fluid, lung and nasal ",
      "tissue was the protection endpoint; 103 data points entered the ",
      "exploratory analysis. Residual error was weighted by the inverse of ",
      "the number of animals in each trial arm."
    )
  )

  # ==========================================================================
  # ini(): Kandala 2022 Table 2 (RM MBMA parameter estimates). All structural
  # parameters are final maximum-likelihood estimates. The Hill coefficient,
  # VLmax and the two IT50 values are on the model's native scales (VLmax and
  # VLmin in log10 viral-load units; IT50 in log2 SN-titre units); none is
  # further log-transformed here. Reported 95% confidence intervals are
  # carried in the trailing comments.
  # ==========================================================================
  ini({
    vlmax <- 4.86
    label("Maximum log10 post-challenge peak viral load (VLmax, log10 copies)")  # Kandala 2022 Table 2, row 'VLmax': 4.86, 95% CI [4.13, 5.60]

    hill <- 4.14
    label("Hill coefficient of the SN-titre viral-load suppression curve (unitless)")  # Kandala 2022 Table 2, row 'g': 4.14, 95% CI [2.83, 5.45]

    it50 <- 5.84
    label("Log2 SN titre giving 50% viral-load suppression in BAL fluid / lung tissue (IT50, log2 titre)")  # Kandala 2022 Table 2, row 'IT50': 5.84, 95% CI [5.21, 6.47]

    it50_nasal <- 8.32
    label("Log2 SN titre giving 50% viral-load suppression in nasal tissue (IT50_nasal, log2 titre)")  # Kandala 2022 Table 2, row 'IT50,nasal': 8.32, 95% CI [7.11, 9.53]

    # VLmin is the minimum log10 peak viral load, fixed per study and specimen
    # to half the reported LLOQ/LLOD (Kandala 2022 Methods, Eq 2 definition).
    # A single across-study typical VLmin is required for a self-contained
    # model; the across-study mean is not printed as a number. It back-solves
    # to 1.6 log10 from the paper's own answer key: the headline 50% protective
    # titre of 24.3% of HCS (Figure 3) together with VLmax 4.86, Hill 4.14,
    # IT50 5.84 and the calibration factor Fc 143.9 force VLmin = 1.61 by
    # uniroot (at the 50% point VL equals (VLmin + VLmax)/2). Fixed, and not a
    # directly published value -- see the vignette Assumptions and deviations.
    vlmin <- fixed(1.6)
    label("Minimum log10 post-challenge peak viral load (VLmin, log10 copies); not a directly published value (back-solved)")  # not from source as a single number -- back-solved to 1.6 from Figure 3 + Table 2 + Fc; see vignette Assumptions and deviations

    # Study-level exponential random effect on VLmax (Kandala 2022 Eq 2:
    # VLmax_i,k,l = VLmax_i,k * exp(eta_i,l), eta ~ N(0, omega^2)). Table 2
    # reports it as a standard deviation (0.24), so the variance is 0.24^2.
    # Named eta_study_* because it is BETWEEN-STUDY (inter-study / challenge-
    # inoculum), not between-subject.
    eta_study_vlmax ~ 0.0576   # Kandala 2022 Table 2, row 'IIV on VLmax': 0.24 (SD), 95% CI [0.15, 0.39]; variance = 0.24^2

    # Additive residual on log10 viral load, weighted by 1/n_arm in the paper.
    addSd <- 1.91
    label("Additive residual SD on log10 peak viral load (log10 copies)")  # Kandala 2022 Table 2, row 'Residual standard error': 1.91, 95% CI [1.65, 2.22]
  })

  model({
    # ======================================================================
    # Specimen-specific 50% protective titre (Kandala 2022 Results: tissue
    # type as a covariate on IT50). SPECIMEN_NASAL = 1 selects nasal tissue.
    # ======================================================================
    it50_k <- it50 * (1 - SPECIMEN_NASAL) + it50_nasal * SPECIMEN_NASAL

    # ======================================================================
    # Study-level VLmax with its exponential random effect (Kandala 2022
    # Equation 2).
    # ======================================================================
    vlmax_i <- vlmax * exp(eta_study_vlmax)

    # ======================================================================
    # Sigmoid (Hill) suppression of peak viral load by SN titre (Kandala 2022
    # Equation 3):
    #   suppr = log2[SN]^g / (log2[SN]^g + log2[IT50]_k^g)
    # SN_SARSCOV2_LOG2 is already the log2 titre, and IT50 is in log2 units.
    # ======================================================================
    suppr <- SN_SARSCOV2_LOG2^hill /
      (SN_SARSCOV2_LOG2^hill + it50_k^hill)

    # ======================================================================
    # Peak viral load (Kandala 2022 Equation 1):
    #   VL = VLmax * exp( ln(VLmin / VLmax) * suppr )
    # suppr = 0 (no titre) gives VL = VLmax; suppr = 1 (saturating titre)
    # gives VL = VLmin; suppr = 0.5 gives the geometric mean of the two.
    # ======================================================================
    log10_viral_load <- vlmax_i * exp(log(vlmin / vlmax_i) * suppr)

    log10_viral_load ~ add(addSd)
  })
}
