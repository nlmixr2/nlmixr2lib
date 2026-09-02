Almond_2016_carbamazepine_invitro <- function() {
  description <- paste(
    "In vitro (cryopreserved human hepatocytes, four donors). Three-parameter",
    "Emax concentration-response model for CYP3A4 induction by carbamazepine,",
    "measured in parallel as the fold increase in CYP3A4 mRNA and as the fold",
    "increase in CYP3A4 catalytic activity (6-beta-hydroxytestosterone",
    "formation) over vehicle control. The model is",
    "fold = 1 + emax * C / (ec50 + C),",
    "the induction term of the Simcyp enzyme-turnover equation (Almond 2016",
    "Eqs. 3 and 4) written on its own. The source reports efficacy as Indmax,",
    "the maximum FOLD induction, defined by the paper as Emax + 1; this file",
    "carries emax = Indmax - 1 and ec50 = IndC50.",
    "Inter-individual variability is BETWEEN HEPATOCYTE DONOR (n = 4), derived",
    "from the donor-level standard deviations in Table 2.",
    "Carbamazepine is one of the three non-rifampicin inducers the paper",
    "carried forward into clinical DDI predictions (Table 4: simvastatin,",
    "quinidine and zolpidem victims), using these values both uncalibrated and",
    "calibrated against the rifampicin reference by Eqs. 5 and 6.",
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
    dosing        = "(none; static in vitro concentration-response model driven by an external carbamazepine concentration covariate)",
    concentration = "(both observations are dimensionless fold-induction over vehicle control; the driving covariate CP_CARBAMAZEPINE_UM is the nominal carbamazepine concentration in the culture medium in uM)"
  )

  covariateData <- list(
    CP_CARBAMAZEPINE_UM = list(
      description        = "Nominal carbamazepine concentration in the hepatocyte culture medium, supplied as a covariate. In this in-vitro model the column carries a culture-medium concentration rather than a plasma concentration, but the quantity and units (uM) are identical.",
      units              = "umol/L (uM)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Almond 2016 Table 1: carbamazepine was assayed at seven final concentrations in culture medium containing 0.1 percent dimethyl sulfoxide (v/v) -- 1, 3, 10, 30, 100, 300 and 1000 uM.",
        "The top concentration is roughly 17-fold above the fitted IndC50 of about 59 uM on both endpoints, so the plateau is less well determined than for rifampicin; this shows up as the wide donor-level S.D. on the activity IndC50 (37.3 on a mean of 59.1).",
        "Carbamazepine MW = 236.3 g/mol (Supplemental Table 2), so 1 uM = 0.2363 mg/L.",
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
    dose_range       = "Nominal carbamazepine concentrations 1, 3, 10, 30, 100, 300 and 1000 uM in culture medium with 0.1 percent dimethyl sulfoxide (v/v) (Table 1).",
    regions          = NA_character_,
    notes            = paste(
      "Cryopreserved human hepatocytes from four donors (Hu1206, Hu1191, Hu1198, Hu4193; Life Technologies) were incubated with serial dilutions of the inducer prepared daily in dimethyl sulfoxide (Materials and Methods).",
      "Two endpoints were measured in parallel in the same incubations: CYP3A4 catalytic activity as 6-beta-hydroxytestosterone formation by liquid chromatography-tandem mass spectrometry, and CYP3A4 mRNA by the QuantiGene Plex 2.0 assay. Because both readouts come from the same cells in the same experiment, they are two outputs of one model rather than two model files.",
      "Curve fitting was carried out on each donor individually and the mean Indmax and IndC50 then calculated across donors; Table 2 reports those means with the donor-level standard deviations, which is what this file's inter-individual variability encodes.",
      "Both a three-parameter (Hill exponent constrained to 1) and a four-parameter sigmoidal model were fitted in GraphPad Prism version 5, and the paper states the two were 'not significantly different', so the three-parameter fit was used. This file therefore carries NO Hill exponent by the authors' own model-selection decision.",
      "Indmax is NOT baseline-corrected: the paper states it 'is equal to Emax + 1'. This is the reason the file carries emax = Indmax - 1.",
      "Carbamazepine autoinduces: the Simcyp simulations accounted for autoinduction because its metabolism was adequately defined (Materials and Methods, 'PBPK Modeling'). Supplemental Table 2 further assumes the active metabolite carbamazepine-10,11-epoxide has induction parameters 'Equal to parent', citing equipotency data (Oscarson 2006). Neither the autoinduction time course nor the metabolite model is part of this file -- both live in the unreproducible Simcyp layer.",
      "In the paper's DDI predictions these in vitro values were used in four ways: uncalibrated, and calibrated against the rifampicin in vivo reference with Indmax 8 or with the refined Indmax 16 (Eqs. 5 and 6; Table 6). The calibration is a transformation of these values, not a separate fit, and is demonstrated in the validation vignette rather than carried as extra model files."
    )
  )

  ini({
    # =====================================================================
    # CYP3A4 ACTIVITY endpoint (6-beta-hydroxytestosterone formation).
    # Almond 2016 Table 2, 'Activity' columns, carbamazepine row:
    #   Indmax mean 16.6 (S.D. 6.1) fold; IndC50 mean 59.1 (S.D. 37.3) uM.
    # =====================================================================

    lemax_activity <- log(15.6)
    label("Log maximal FRACTIONAL increase in CYP3A4 activity over vehicle control (dimensionless; Emax = Indmax - 1)")
    # Table 2: Indmax (activity) = 16.6 fold, so emax = 16.6 - 1 = 15.6.
    # The paper defines Indmax as the maximum fold induction, "equal to
    # Emax + 1" (Materials and Methods, 'In Vitro Data Analysis').

    lec50_activity <- log(59.1)
    label("Log carbamazepine concentration producing half-maximal CYP3A4 activity induction (uM)")
    # Table 2: IndC50 (activity) = 59.1 uM.

    # =====================================================================
    # CYP3A4 mRNA endpoint (QuantiGene Plex 2.0).
    # Almond 2016 Table 2, 'mRNA' columns, carbamazepine row:
    #   Indmax mean 21.9 (S.D. 12.4) fold; IndC50 mean 58.7 (S.D. 18.0) uM.
    # =====================================================================

    lemax_mrna <- log(20.9)
    label("Log maximal FRACTIONAL increase in CYP3A4 mRNA over vehicle control (dimensionless; Emax = Indmax - 1)")
    # Table 2: Indmax (mRNA) = 21.9 fold, so emax = 21.9 - 1 = 20.9.
    # Narrative check (Results, 'Induction Parameters Determined In Vitro'):
    # mRNA efficacy was 1.3- to 2.0-fold higher than activity efficacy across
    # the six inducers; here 21.9 / 16.6 = 1.32, inside that range.

    lec50_mrna <- log(58.7)
    label("Log carbamazepine concentration producing half-maximal CYP3A4 mRNA induction (uM)")
    # Table 2: IndC50 (mRNA) = 58.7 uM. Narrative check: mRNA potency was
    # 1.0- to 3.3-fold lower (IndC50 higher) than activity potency; here
    # 58.7 / 59.1 = 0.99, which is the 1.0-fold END of that stated range --
    # carbamazepine is the inducer that defines the lower bound of the
    # paper's own quoted IndC50 ratio range.

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

    etalemax_activity ~ 0.142281
    # Table 2: S.D. 6.1 on Indmax 16.6 -> CV = 6.1 / 15.6 = 0.3910 on emax;
    # omega^2 = log(1 + 0.3910^2) = 0.142281.
    etalec50_activity ~ 0.335278
    # Table 2: S.D. 37.3 on IndC50 59.1 -> CV = 0.6311;
    # omega^2 = log(1 + 0.6311^2) = 0.335278.
    etalemax_mrna ~ 0.301590
    # Table 2: S.D. 12.4 on Indmax 21.9 -> CV = 12.4 / 20.9 = 0.5933 on emax;
    # omega^2 = log(1 + 0.5933^2) = 0.301590.
    etalec50_mrna ~ 0.089869
    # Table 2: S.D. 18.0 on IndC50 58.7 -> CV = 0.3066;
    # omega^2 = log(1 + 0.3066^2) = 0.089869.

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
    foldActivity <- 1 + emax_activity * CP_CARBAMAZEPINE_UM /
      (ec50_activity + CP_CARBAMAZEPINE_UM)
    foldMrna <- 1 + emax_mrna * CP_CARBAMAZEPINE_UM /
      (ec50_mrna + CP_CARBAMAZEPINE_UM)

    foldActivity ~ add(addSd_foldActivity)
    foldMrna ~ add(addSd_foldMrna)
  })
}
