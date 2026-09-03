Almond_2016_nifedipine_invitro <- function() {
  description <- paste(
    "In vitro (cryopreserved human hepatocytes, four donors). Three-parameter",
    "Emax concentration-response model for CYP3A4 induction by nifedipine,",
    "measured in parallel as the fold increase in CYP3A4 mRNA and as the fold",
    "increase in CYP3A4 catalytic activity (6-beta-hydroxytestosterone",
    "formation) over vehicle control. The model is",
    "fold = 1 + emax * C / (ec50 + C),",
    "the induction term of the Simcyp enzyme-turnover equation (Almond 2016",
    "Eqs. 3 and 4) written on its own. The source reports efficacy as Indmax,",
    "the maximum FOLD induction, defined by the paper as Emax + 1; this file",
    "carries emax = Indmax - 1 and ec50 = IndC50.",
    "Inter-individual variability is BETWEEN HEPATOCYTE DONOR (n = 4), derived",
    "from the donor-level standard deviations in Table 2, and is the WIDEST of",
    "the six inducers on both efficacy parameters (donor-level CV about 76",
    "percent on emax for each endpoint).",
    "Nifedipine occupies an unusual dual role in this paper: it is a CYP3A4",
    "SUBSTRATE (victim) in the clinical DDI simulations with rifampicin and",
    "phenobarbital, and separately a CYP3A4 INDUCER characterised in vitro",
    "here. The paper ran no DDI simulations with nifedipine as perpetrator, so",
    "these induction parameters carry no in vivo verification in this study.",
    "The Simcyp minimal-PBPK disposition models and the CYP3A4 enzyme-turnover",
    "ODEs are NOT part of this file and are not reproducible from the published",
    "inputs: the enzyme degradation rate constants kdegH-3A4 and kdegG-3A4,",
    "the basal CYP3A4 abundances, MPPGL, liver weight and organ blood flows",
    "are never reported. See the validation vignette for the full accounting.",
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
    "seventh columns of Table 3; Table 3 is not a source of any value here).",
    "Induction parameters from Table 2 (mean and S.D. over four hepatocyte",
    "donors); incubation concentrations from Table 1; assay and curve-fitting",
    "methods from Materials and Methods; fold-induction functional form from",
    "Eqs. 3 and 4.",
    sep = " "
  )

  vignette <- "Almond_2016_cyp3a4_induction"

  units <- list(
    time          = "h",
    dosing        = "(none; static in vitro concentration-response model driven by an external nifedipine concentration covariate)",
    concentration = "(both observations are dimensionless fold-induction over vehicle control; the driving covariate CP_NIFEDIPINE_UM is the nominal nifedipine concentration in the culture medium in uM)"
  )

  covariateData <- list(
    CP_NIFEDIPINE_UM = list(
      description        = "Nominal nifedipine concentration in the hepatocyte culture medium, supplied as a covariate. In this in-vitro model the column carries a culture-medium concentration rather than a plasma concentration, but the quantity and units (uM) are identical.",
      units              = "umol/L (uM)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Almond 2016 Table 1: nifedipine was assayed at seven final concentrations in culture medium containing 0.1 percent dimethyl sulfoxide (v/v) -- 0.03, 1, 2, 3, 10, 30 and 100 uM. Note the large gap between the lowest point (0.03 uM) and the next (1 uM); the remaining six points bracket the fitted activity IndC50 of 4.0 uM.",
        "Nifedipine MW = 346.3 g/mol (Supplemental Table 1, where nifedipine appears as a VICTIM drug), so 1 uM = 0.3463 mg/L.",
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
    dose_range       = "Nominal nifedipine concentrations 0.03, 1, 2, 3, 10, 30 and 100 uM in culture medium with 0.1 percent dimethyl sulfoxide (v/v) (Table 1).",
    regions          = NA_character_,
    notes            = paste(
      "Cryopreserved human hepatocytes from four donors (Hu1206, Hu1191, Hu1198, Hu4193; Life Technologies) were incubated with serial dilutions of the inducer prepared daily in dimethyl sulfoxide (Materials and Methods).",
      "Two endpoints were measured in parallel in the same incubations: CYP3A4 catalytic activity as 6-beta-hydroxytestosterone formation by liquid chromatography-tandem mass spectrometry, and CYP3A4 mRNA by the QuantiGene Plex 2.0 assay. Because both readouts come from the same cells in the same experiment, they are two outputs of one model rather than two model files.",
      "Curve fitting was carried out on each donor individually and the mean Indmax and IndC50 then calculated across donors; Table 2 reports those means with the donor-level standard deviations, which is what this file's inter-individual variability encodes.",
      "Both a three-parameter (Hill exponent constrained to 1) and a four-parameter sigmoidal model were fitted in GraphPad Prism version 5, and the paper states the two were 'not significantly different', so the three-parameter fit was used. This file therefore carries NO Hill exponent by the authors' own model-selection decision.",
      "Indmax is NOT baseline-corrected: the paper states it 'is equal to Emax + 1'. This is the reason the file carries emax = Indmax - 1.",
      "WIDEST BETWEEN-DONOR SPREAD IN THE STUDY. The activity Indmax is 15.6 with an S.D. of 11.3 and the mRNA Indmax is 30.0 with an S.D. of 22.0 -- donor-level CVs of about 76 percent on emax for both endpoints. Simulated nifedipine induction curves are correspondingly dispersed, and a four-donor mean with that spread should be treated as a weakly determined central estimate.",
      "SCOPE CAVEAT. Nifedipine is both a victim and a perpetrator in this paper, but only its VICTIM role reaches the clinical simulations: Table 3 (rifampicin, oral and i.v. nifedipine) and Table 4 (phenobarbital + nifedipine) treat it as a CYP3A4 substrate, and Supplemental Table 1 carries its victim-side compound file. It never appears as a perpetrator in Table 4, Fig. 5 or Table 6, so these induction parameters have no in vivo verification within this study.",
      "The calibration equations (Eqs. 5 and 6) can nonetheless be applied to these values, since calibration only requires the test compound's in vitro parameters and rifampicin's in vitro and in vivo parameters -- all of which are on disk. The vignette demonstrates that calculation for nifedipine alongside the three inducers the paper itself calibrated."
    )
  )

  ini({
    # =====================================================================
    # CYP3A4 ACTIVITY endpoint (6-beta-hydroxytestosterone formation).
    # Almond 2016 Table 2, 'Activity' columns, nifedipine row:
    #   Indmax mean 15.6 (S.D. 11.3) fold; IndC50 mean 4.0 (S.D. 1.9) uM.
    # =====================================================================

    lemax_activity <- log(14.6)
    label("Log maximal FRACTIONAL increase in CYP3A4 activity over vehicle control (dimensionless; Emax = Indmax - 1)")
    # Table 2: Indmax (activity) = 15.6 fold, so emax = 15.6 - 1 = 14.6.
    # The paper defines Indmax as the maximum fold induction, "equal to
    # Emax + 1" (Materials and Methods, 'In Vitro Data Analysis').

    lec50_activity <- log(4.0)
    label("Log nifedipine concentration producing half-maximal CYP3A4 activity induction (uM)")
    # Table 2: IndC50 (activity) = 4.0 uM -- the most potent of the six
    # inducers on the activity endpoint apart from rifampicin (0.30 uM).

    # =====================================================================
    # CYP3A4 mRNA endpoint (QuantiGene Plex 2.0).
    # Almond 2016 Table 2, 'mRNA' columns, nifedipine row:
    #   Indmax mean 30.0 (S.D. 22.0) fold; IndC50 mean 13.0 (S.D. 9.5) uM.
    # =====================================================================

    lemax_mrna <- log(29.0)
    label("Log maximal FRACTIONAL increase in CYP3A4 mRNA over vehicle control (dimensionless; Emax = Indmax - 1)")
    # Table 2: Indmax (mRNA) = 30.0 fold, so emax = 30.0 - 1 = 29.0.
    # Narrative check (Results, 'Induction Parameters Determined In Vitro'):
    # mRNA efficacy was 1.3- to 2.0-fold higher than activity efficacy; here
    # 30.0 / 15.6 = 1.92, near the top of that stated range.

    lec50_mrna <- log(13.0)
    label("Log nifedipine concentration producing half-maximal CYP3A4 mRNA induction (uM)")
    # Table 2: IndC50 (mRNA) = 13.0 uM. Narrative check: mRNA potency was
    # 1.0- to 3.3-fold lower (IndC50 higher) than activity potency; here
    # 13.0 / 4.0 = 3.25, which is the UPPER end of that stated range --
    # nifedipine defines the 3.3-fold bound the paper quotes.

    # =====================================================================
    # BETWEEN-DONOR variability (n = 4 hepatocyte donors).
    #
    # Table 2 reports an arithmetic mean and a standard deviation over the
    # four donor-level fits; each variance below is the log-normal variance
    # implied by that coefficient of variation, omega^2 = log(1 + CV^2).
    #
    # For the efficacy parameters the CV is taken on emax, NOT on Indmax:
    # emax = Indmax - 1 is a SHIFT, so the reported S.D. carries over
    # unchanged while the mean drops by 1, and the CV must be recomputed on
    # the shifted mean.
    #
    # The source states no distributional form for the between-donor spread,
    # so the log-normal reading is an extraction choice flagged in the
    # vignette Errata. No correlations are reported; the block is diagonal.
    # =====================================================================

    etalemax_activity ~ 0.469399
    # Table 2: S.D. 11.3 on Indmax 15.6 -> CV = 11.3 / 14.6 = 0.7740 on emax;
    # omega^2 = log(1 + 0.7740^2) = 0.469399. Largest efficacy variance in
    # the study on the activity endpoint.
    etalec50_activity ~ 0.203451
    # Table 2: S.D. 1.9 on IndC50 4.0 -> CV = 0.4750;
    # omega^2 = log(1 + 0.4750^2) = 0.203451.
    etalemax_mrna ~ 0.454576
    # Table 2: S.D. 22.0 on Indmax 30.0 -> CV = 22.0 / 29.0 = 0.7586 on emax;
    # omega^2 = log(1 + 0.7586^2) = 0.454576.
    etalec50_mrna ~ 0.427894
    # Table 2: S.D. 9.5 on IndC50 13.0 -> CV = 0.7308;
    # omega^2 = log(1 + 0.7308^2) = 0.427894.

    # =====================================================================
    # Residual error is NOT reported: the source fitted each donor's curve
    # in GraphPad Prism and reports only point estimates and between-donor
    # S.D. Per the standing policy on unreported residual error, both are
    # fixed at 0 so the model returns the deterministic curve per donor.
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
    # This is the induction term of the paper's enzyme-turnover equations
    # (Eqs. 3 and 4), which read 1 + (Indmax - 1) * I / (IndC50 + I).
    #
    # Structural checks asserted in the validation vignette:
    #   fold(0)    = 1                 (vehicle control)
    #   fold(ec50) = 1 + emax / 2      (definition of IndC50)
    #   fold(Inf)  = 1 + emax = Indmax (reported maximum fold induction)
    # ===================================================================
    foldActivity <- 1 + emax_activity * CP_NIFEDIPINE_UM /
      (ec50_activity + CP_NIFEDIPINE_UM)
    foldMrna <- 1 + emax_mrna * CP_NIFEDIPINE_UM /
      (ec50_mrna + CP_NIFEDIPINE_UM)

    foldActivity ~ add(addSd_foldActivity)
    foldMrna ~ add(addSd_foldMrna)
  })
}
