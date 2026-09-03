Almond_2016_rifampicin_invitro <- function() {
  description <- paste(
    "In vitro (cryopreserved human hepatocytes, four donors). Three-parameter",
    "Emax concentration-response model for CYP3A4 induction by rifampicin,",
    "measured in parallel as the fold increase in CYP3A4 mRNA and as the fold",
    "increase in CYP3A4 catalytic activity (6-beta-hydroxytestosterone",
    "formation) over vehicle control. The model is",
    "fold = 1 + emax * C / (ec50 + C),",
    "which is the induction term of the Simcyp enzyme-turnover equation",
    "(Almond 2016 Eqs. 3 and 4) written on its own. The source reports the",
    "induction efficacy as Indmax, the maximum FOLD induction, which the paper",
    "defines explicitly as Emax + 1; this file therefore carries",
    "emax = Indmax - 1 and ec50 = IndC50, following the same translation used",
    "by Willemin_2024_interleukin6_cyp_pbpk.R and Chen_2024_interleukin6_cyp3a_pbpk.R.",
    "Inter-individual variability is BETWEEN HEPATOCYTE DONOR (n = 4), derived",
    "from the donor-level standard deviations in Table 2.",
    "This is the paper's own bench work and is the calibration standard for the",
    "rest of the study: the companion file Almond_2016_rifampicin_invivo.R",
    "carries the in vivo reference induction parameters for rifampicin against",
    "which the other inducers' in vitro values were calibrated (Eqs. 5 and 6).",
    "The Simcyp minimal-PBPK victim/perpetrator disposition models and the",
    "CYP3A4 enzyme-turnover ODEs themselves are NOT part of this file and are",
    "not reproducible from the published inputs: the paper never reports the",
    "enzyme degradation rate constants kdegH-3A4 and kdegG-3A4, the basal",
    "hepatic or intestinal CYP3A4 abundances, MPPGL, liver weight or organ",
    "blood flows. See the validation vignette for the full accounting.",
    sep = " "
  )

  reference <- paste(
    "Almond LM, Mukadam S, Gardner I, Okialda K, Wong S, Hatley O, Tay S,",
    "Rowland-Yeo K, Jamei M, Rostami-Hodjegan A, Kenny JR.",
    "Prediction of Drug-Drug Interactions Arising from CYP3A Induction Using a",
    "Physiologically Based Dynamic Model.",
    "Drug Metab Dispos. 2016;44(6):821-832. doi:10.1124/dmd.115.066845.",
    "Erratum: Drug Metab Dispos. 2016;44(6):877.",
    "doi:10.1124/dmd.115.066845err (corrects the units of the sixth and",
    "seventh columns of Table 3 from ng/L.h to ng/mL.h; Table 3 carries the",
    "clinical DDI exposures and is not a source of any value in this file).",
    "Induction parameters from Table 2 (mean and S.D. over four hepatocyte",
    "donors); incubation concentrations from Table 1; assay and curve-fitting",
    "methods from Materials and Methods, 'Generation of Induction Parameters",
    "In Vitro' and 'In Vitro Data Analysis'; the fold-induction functional form",
    "from Eqs. 3 and 4 (Materials and Methods, 'PBPK Modeling').",
    sep = " "
  )

  vignette <- "Almond_2016_cyp3a4_induction"

  units <- list(
    time          = "h",
    dosing        = "(none; static in vitro concentration-response model driven by an external rifampicin concentration covariate)",
    concentration = "(both observations are dimensionless fold-induction over vehicle control; the driving covariate CP_RIF_UM is the nominal rifampicin concentration in the culture medium in uM)"
  )

  covariateData <- list(
    CP_RIF_UM = list(
      description        = "Nominal rifampicin concentration in the hepatocyte culture medium, supplied as a covariate. Reused canonical: in this in-vitro model the column carries a culture-medium concentration rather than a plasma concentration, but the quantity and units (uM) are identical. Same reuse rationale as CP_TUVUSERTIB_NGML in Mukker_2026_tuvusertib_hERG.R, which carries a patch-clamp bath concentration.",
      units              = "umol/L (uM)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Almond 2016 Table 1: rifampicin was assayed at seven final concentrations in culture medium containing 0.1 percent dimethyl sulfoxide (v/v) -- 0.03, 0.1, 0.3, 1, 3, 10 and 30 uM. Concentration ranges were chosen per inducer 'with the aim of determining a robust Indmax and IndC50' (Materials and Methods).",
        "The top concentration is 100-fold above the fitted activity IndC50 of 0.30 uM, so the plateau of the curve is well determined; this is why rifampicin has the tightest donor-level S.D. on Indmax of the six inducers.",
        "Rifampicin MW = 823 g/mol (Supplemental Table 2), so 1 uM = 0.823 mg/L = 823 ng/mL.",
        "Set to 0 for the vehicle-control condition, at which the model returns fold = 1 by construction."
      ),
      source_name        = "concentration of the inducer"
    )
  )

  population <- list(
    species          = "in vitro (cryopreserved human hepatocytes, four donors: Hu1206, Hu1191, Hu1198, Hu4193)",
    n_subjects       = 4L,
    n_studies        = 1L,
    age_range        = NA_character_,
    weight_range     = NA_character_,
    sex_female_pct   = NA_real_,
    race_ethnicity   = NA_character_,
    disease_state    = "Not applicable -- primary human hepatocytes in culture.",
    dose_range       = "Nominal rifampicin concentrations 0.03, 0.1, 0.3, 1, 3, 10 and 30 uM in culture medium with 0.1 percent dimethyl sulfoxide (v/v) (Table 1).",
    regions          = NA_character_,
    notes            = paste(
      "Cryopreserved human hepatocytes from four donors (Hu1206, Hu1191, Hu1198, Hu4193; Life Technologies) were incubated with serial dilutions of the inducer prepared daily in dimethyl sulfoxide (Materials and Methods, 'Materials' and 'Generation of Induction Parameters In Vitro').",
      "TWO ENDPOINTS WERE MEASURED IN PARALLEL IN THE SAME INCUBATIONS: CYP3A4 catalytic activity, as 6-beta-hydroxytestosterone formation quantified by liquid chromatography-tandem mass spectrometry, and CYP3A4 mRNA, quantified with the QuantiGene Plex 2.0 assay (Affymetrix panel 11477). Because the two readouts come from the same cells in the same experiment, they are carried as two outputs of one model rather than as two model files.",
      "Cell toxicity and viability were monitored by lactate dehydrogenase leakage and AlamarBlue.",
      "Curve fitting was carried out on the data from EACH DONOR INDIVIDUALLY, and the mean Indmax and IndC50 were then calculated across donors (Materials and Methods, 'In Vitro Data Analysis'). Table 2 reports those means with the corresponding donor-level standard deviations, which is what this file's inter-individual variability encodes.",
      "Both a three-parameter model (Hill exponent constrained to 1) and a four-parameter sigmoidal model were fitted in GraphPad Prism version 5. The paper states that parameters derived from the two models 'were not significantly different; therefore, the values from the simpler model (three-parameter fit) were used for subsequent analysis'. This file therefore carries NO Hill exponent: the Hill exponent is 1 by the authors' own model-selection decision, not by an assumption made during extraction.",
      "Indmax is the maximum fold induction and is NOT corrected for baseline -- the paper states 'Indmax is the maximum fold induction and as such is not corrected for baseline (i.e., is equal to Emax + 1)'. This is the single most important reading trap in the source and is the reason this file carries emax = Indmax - 1.",
      "Rifampicin is the paper's reference inducer. Its in vitro values are the denominator of the calibration equations (Eqs. 5 and 6) applied to the other inducers, and its separately derived IN VIVO reference values are carried in Almond_2016_rifampicin_invivo.R."
    )
  )

  ini({
    # =====================================================================
    # CYP3A4 ACTIVITY endpoint (6-beta-hydroxytestosterone formation).
    # Almond 2016 Table 2, 'Activity' columns, rifampicin row:
    #   Indmax mean 22.7 (S.D. 7.8) fold; IndC50 mean 0.30 (S.D. 0.10) uM.
    # =====================================================================

    lemax_activity <- log(21.7)
    label("Log maximal FRACTIONAL increase in CYP3A4 activity over vehicle control (dimensionless; Emax = Indmax - 1)")
    # Table 2: Indmax (activity) = 22.7 fold. The paper defines Indmax as the
    # maximum fold induction, "equal to Emax + 1" (Materials and Methods,
    # 'In Vitro Data Analysis'), so emax = 22.7 - 1 = 21.7. The `emax` /
    # `ec50` naming (rather than a new `indmax` / `indc50` canonical) follows
    # the established translation in Willemin_2024_interleukin6_cyp_pbpk.R,
    # which converts each Simcyp Indmax to emax = Indmax - 1.

    lec50_activity <- log(0.30)
    label("Log rifampicin concentration producing half-maximal CYP3A4 activity induction (uM)")
    # Table 2: IndC50 (activity) = 0.30 uM. Note this is the IN VITRO
    # potency; the separately derived IN VIVO reference IndC50 for rifampicin
    # is 0.32 uM (Almond_2016_rifampicin_invivo.R) -- close agreement, but
    # the two are independent determinations and must not be interchanged.

    # =====================================================================
    # CYP3A4 mRNA endpoint (QuantiGene Plex 2.0).
    # Almond 2016 Table 2, 'mRNA' columns, rifampicin row:
    #   Indmax mean 29.9 (S.D. 7.0) fold; IndC50 mean 0.71 (S.D. 0.35) uM.
    # =====================================================================

    lemax_mrna <- log(28.9)
    label("Log maximal FRACTIONAL increase in CYP3A4 mRNA over vehicle control (dimensionless; Emax = Indmax - 1)")
    # Table 2: Indmax (mRNA) = 29.9 fold, so emax = 29.9 - 1 = 28.9.
    # Consistency check against the paper's own narrative (Results,
    # 'Induction Parameters Determined In Vitro'): efficacy measured by mRNA
    # was 1.3- to 2.0-fold higher than by activity across the six inducers.
    # Here 29.9 / 22.7 = 1.32, at the bottom of that stated range.

    lec50_mrna <- log(0.71)
    label("Log rifampicin concentration producing half-maximal CYP3A4 mRNA induction (uM)")
    # Table 2: IndC50 (mRNA) = 0.71 uM. Same narrative check: potency
    # measured by mRNA was 1.0- to 3.3-fold LOWER (i.e. IndC50 higher) than
    # by activity. Here 0.71 / 0.30 = 2.37, inside that stated range.

    # =====================================================================
    # BETWEEN-DONOR variability (n = 4 hepatocyte donors).
    #
    # Table 2 reports an arithmetic mean and a standard deviation over the
    # four donor-level fits. Each variance below is the log-normal variance
    # implied by that coefficient of variation, omega^2 = log(1 + CV^2).
    #
    # For the two efficacy parameters the CV is taken on emax, NOT on
    # Indmax: emax = Indmax - 1 is a SHIFT, so the reported S.D. carries over
    # unchanged while the mean drops by 1, and the CV must be recomputed on
    # the shifted mean. Using the unshifted CV would understate every
    # efficacy variance in this file.
    #
    # The source states no distributional form for the between-donor spread
    # (it reports mean +/- S.D. over n = 4 only), so the log-normal reading
    # is an extraction choice; it is flagged in the vignette Errata. No
    # correlation between parameters is reported, so the block is diagonal.
    # =====================================================================

    etalemax_activity ~ 0.121511
    # Table 2: S.D. 7.8 on Indmax 22.7 -> CV = 7.8 / 21.7 = 0.3594 on emax;
    # omega^2 = log(1 + 0.3594^2) = 0.121511.
    etalec50_activity ~ 0.105361
    # Table 2: S.D. 0.10 on IndC50 0.30 -> CV = 0.3333;
    # omega^2 = log(1 + 0.3333^2) = 0.105361.
    etalemax_mrna ~ 0.057011
    # Table 2: S.D. 7.0 on Indmax 29.9 -> CV = 7.0 / 28.9 = 0.2422 on emax;
    # omega^2 = log(1 + 0.2422^2) = 0.057011.
    etalec50_mrna ~ 0.217534
    # Table 2: S.D. 0.35 on IndC50 0.71 -> CV = 0.4930;
    # omega^2 = log(1 + 0.4930^2) = 0.217534.

    # =====================================================================
    # Residual error is NOT reported. The source fitted each donor's
    # concentration-response curve in GraphPad Prism and reports only the
    # resulting point estimates and their between-donor S.D.; there is no
    # residual-error model, no replicate-level variance and no goodness-of-
    # fit statistic. Per the standing policy on unreported residual error,
    # both are fixed at 0 so the model returns the deterministic curve for
    # each simulated donor. Flagged in the vignette Errata.
    # =====================================================================
    addSd_foldActivity <- fixed(0)
    label("Additive residual error on fold CYP3A4 activity induction (ZERO - not reported in source)")
    addSd_foldMrna <- fixed(0)
    label("Additive residual error on fold CYP3A4 mRNA induction (ZERO - not reported in source)")
  })

  model({
    emax_activity <- exp(lemax_activity + etalemax_activity)
    ec50_activity <- exp(lec50_activity + etalec50_activity)
    emax_mrna     <- exp(lemax_mrna + etalemax_mrna)
    ec50_mrna     <- exp(lec50_mrna + etalec50_mrna)

    # ===================================================================
    # Three-parameter Emax (Hill exponent 1) concentration-response for
    # CYP3A4 induction, expressed as fold over vehicle control:
    #
    #   fold(C) = 1 + emax * C / (ec50 + C)
    #
    # This is exactly the induction term of the paper's enzyme-turnover
    # equations (Eqs. 3 and 4), which read
    #   1 + (Indmax - 1) * I / (IndC50 + I),
    # with Indmax - 1 = emax and IndC50 = ec50.
    #
    # Structural checks asserted in the validation vignette:
    #   fold(0)    = 1               (vehicle control, by construction)
    #   fold(ec50) = 1 + emax / 2    (definition of the reported IndC50:
    #                                 "the concentration that yields half
    #                                 of the Emax")
    #   fold(Inf)  = 1 + emax = Indmax (the reported maximum fold induction)
    # ===================================================================
    foldActivity <- 1 + emax_activity * CP_RIF_UM / (ec50_activity + CP_RIF_UM)
    foldMrna     <- 1 + emax_mrna * CP_RIF_UM / (ec50_mrna + CP_RIF_UM)

    foldActivity ~ add(addSd_foldActivity)
    foldMrna ~ add(addSd_foldMrna)
  })
}
