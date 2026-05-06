Zecchin_2016_tumorovarian <- function() {
  description <- "Tumour-size dynamics model (sum of longest diameters, SLD) for advanced epithelial ovarian cancer with independent additive carboplatin and gemcitabine cytotoxic effects (Zecchin 2016 / DDMODEL00000217): exponential SLD growth with two drug-exposure-driven death-rate terms (no resistance, no synergy), additive residual error, and M3-method handling of below-LLOQ observations in the source NONMEM run."
  reference   <- paste(
    "Zecchin C, Gueorguieva I, Enas NH, Friberg LE.",
    "Models for change in tumour size, appearance of new lesions",
    "and survival probability in patients with advanced epithelial",
    "ovarian cancer.",
    "Br J Clin Pharmacol. 2016;82(3):717-727.",
    "doi:10.1111/bcp.12994.",
    "DDMORE Foundation Model Repository: DDMODEL00000217.",
    sep = " "
  )
  vignette    <- "Zecchin_2016_tumorovarian"
  units       <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; drug exposure enters as time-varying per-cycle AUC_CARBO and AUC_GEM covariates)",
    concentration = "mm (tumour-size endpoint, not a drug concentration)"
  )

  ddmore_id    <- "DDMODEL00000217"
  replicate_of <- NULL

  covariateData <- list(
    AUC_CARBO = list(
      description        = "Per-cycle average AUC of carboplatin (time-varying drug-exposure covariate driving the carboplatin cytotoxic-death term).",
      units              = "carboplatin AUC units (mg*min/mL)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Held step-wise constant within each chemotherapy cycle and reset at the start of the next cycle. Set to 0 in cycles where carboplatin is not administered. The DDMORE bundle's Simulated_SLD.csv encodes this column as AUC0; the source $INPUT NM-TRAN column is CB.",
      source_name        = "CB"
    ),
    AUC_GEM = list(
      description        = "Per-cycle average AUC of gemcitabine (sum of parent and active intracellular metabolite exposure, per Zecchin 2016 Methods); time-varying drug-exposure covariate driving the gemcitabine cytotoxic-death term.",
      units              = "gemcitabine AUC units (paper composite mol*day / 10^6 cells)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Held step-wise constant within each chemotherapy cycle and reset at the start of the next cycle. Set to 0 in carboplatin-monotherapy cycles. The DDMORE bundle's Simulated_SLD.csv encodes this column as AUC1; the source $INPUT NM-TRAN column is G.",
      source_name        = "G"
    )
  )

  population <- list(
    n_subjects     = 336L,
    n_studies      = 1L,
    age_range      = "advanced ovarian cancer adult population (specific range not available in this extraction)",
    weight_range   = "not available in this extraction",
    sex_female_pct = 100,
    disease_state  = "advanced (FIGO stage III/IV) epithelial ovarian cancer (recurrent / metastatic)",
    dose_range     = "Phase III chemotherapy: carboplatin monotherapy (target AUC 4-6 mg*min/mL Q3W) or carboplatin + gemcitabine combination per the trial protocol referenced by Zecchin 2016",
    notes          = "336 patients pooled from a randomised Phase III trial comparing carboplatin monotherapy vs carboplatin + gemcitabine combination chemotherapy in advanced ovarian cancer (per PubMed PMID 27136318 abstract for Zecchin 2016, Br J Clin Pharmacol 82(3):717-727). The linked publication PDF was not on disk for this extraction, so detailed baseline demographics (median age, weight, region) could not be transcribed; the population block reports only what the abstract and the DDMORE bundle's NONMEM run-time messages confirm. The model itself uses no patient-level covariates beyond the time-varying drug-AUC inputs; no allometric / age / sex effect was carried into the final SLD model."
  )

  ini({
    # Final estimates from Output_real_SLD.lst FINAL PARAMETER ESTIMATE block
    # (LAPLACIAN CONDITIONAL ESTIMATION WITH INTERACTION). The .lst reports
    # MINIMIZATION SUCCESSFUL with the caveat "HOWEVER, PROBLEMS OCCURRED
    # WITH THE MINIMIZATION" and very high THETA-THETA correlations (~0.98+
    # in the correlation matrix). The values are extracted as published; see
    # vignette Errata for context.
    lkg     <- log(0.611);  label("Tumour growth rate constant (1/day, internal /1000 scaling)")    # Output_real_SLD.lst FINAL TH1 = 6.11E-01
    lkd0    <- log(0.0497); label("Carboplatin-related tumour death rate constant (per AUC_CARBO unit, internal /1000 scaling)") # Output_real_SLD.lst FINAL TH2 = 4.97E-02
    lkd1    <- log(0.0164); label("Gemcitabine-related tumour death rate constant (per AUC_GEM unit, internal /100 scaling)")    # Output_real_SLD.lst FINAL TH3 = 1.64E-02
    libase  <- log(0.0713); label("Baseline SLD (m); converted to mm via internal *1000 scaling")    # Output_real_SLD.lst FINAL TH4 = 7.13E-02
    addSd_tumorSize <- 18.4 ; label("Additive residual error on SLD (mm)")  # Output_real_SLD.lst FINAL TH5 = 1.84E+01

    # IIV — diagonal $OMEGA from Output_real_SLD.lst FINAL OMEGA block. The
    # source NONMEM model (.mod $PK) wires the same ETA(2) onto BOTH KD0 and
    # KD1, so a single shared eta is added to both lkd0 and lkd1 in the
    # model() block via etalkd0_kd1.
    etalkg      ~ 1.72   # Output_real_SLD.lst FINAL OMEGA(1,1) = 1.72E+00 (large CV; ETAshrink ~55%)
    etalkd0_kd1 ~ 1.09   # Output_real_SLD.lst FINAL OMEGA(2,2) = 1.09E+00 (shared eta on KD0 and KD1; ETAshrink ~28%)
    etalibase   ~ 0.515  # Output_real_SLD.lst FINAL OMEGA(3,3) = 5.15E-01 (ETAshrink ~12%)
  })

  model({
    # Individual structural parameters (lognormal IIV). KD0 and KD1 share a
    # single random effect (etalkd0_kd1) per the source $PK block, which
    # writes KD0 = THETA(2)*EXP(ETA(2)) and KD1 = THETA(3)*EXP(ETA(2)).
    kg    <- exp(lkg + etalkg)
    kd0   <- exp(lkd0 + etalkd0_kd1)
    kd1   <- exp(lkd1 + etalkd0_kd1)
    ibase <- exp(libase + etalibase)

    # Initial SLD in mm. The source parameterizes IBASE in metres and
    # initializes A_0(1) = IBASE * 1000 to set the state in mm; the *1000
    # is preserved verbatim so the rxode2 simulation reproduces the source
    # NONMEM trajectory.
    tumorSize(0) <- ibase * 1000

    # Tumour size dynamics — direct port of $DES from Executable_SLD.mod:
    #   DADT(1) = KG/1000 * A(1)
    #             - (KD0/1000 * E0 + KD1/100 * E1) * A(1)
    # AUC_CARBO and AUC_GEM are time-varying covariates supplied per cycle
    # in the input data (per-cycle average drug AUCs). Carboplatin and
    # gemcitabine effects act independently and additively on the
    # cytotoxic death-rate term (no resistance, no synergy). The internal
    # /1000 (KG, KD0) and /100 (KD1) numerical scalings are inherited
    # verbatim from the source; they preserve the numerical equivalence
    # with the Output_real_SLD.lst final estimates and absorb the unit
    # conventions of the source dataset's exposure metrics.
    d/dt(tumorSize) <- kg / 1000 * tumorSize -
      (kd0 / 1000 * AUC_CARBO + kd1 / 100 * AUC_GEM) * tumorSize

    # Additive residual error. The source NONMEM model uses
    # SIGMA = 1 FIX with Y = IPRED + ERR(1)*W and W = FADD = 18.4 mm,
    # which is equivalent to a pure additive residual on SLD with
    # SD = FADD. The M3 method for below-LLOQ observations used by the
    # source $ERROR block is omitted from this nlmixr2 translation
    # because forward simulation does not exercise the censoring
    # likelihood; see vignette Assumptions and deviations.
    tumorSize ~ add(addSd_tumorSize)
  })
}
