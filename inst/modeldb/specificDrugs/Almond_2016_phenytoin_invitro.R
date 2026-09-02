Almond_2016_phenytoin_invitro <- function() {
  description <- paste(
    "In vitro (cryopreserved human hepatocytes, four donors). Three-parameter",
    "Emax concentration-response model for CYP3A4 induction by phenytoin,",
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
    "Phenytoin is one of the three non-rifampicin inducers carried forward into",
    "clinical DDI predictions (Table 4: quinidine victim), using these values",
    "both uncalibrated and calibrated against the rifampicin reference.",
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
    dosing        = "(none; static in vitro concentration-response model driven by an external phenytoin concentration covariate)",
    concentration = "(both observations are dimensionless fold-induction over vehicle control; the driving covariate CP_PHENYTOIN_UM is the nominal phenytoin concentration in the culture medium in uM)"
  )

  covariateData <- list(
    CP_PHENYTOIN_UM = list(
      description        = "Nominal phenytoin concentration in the hepatocyte culture medium, supplied as a covariate. In this in-vitro model the column carries a culture-medium concentration rather than a plasma concentration, but the quantity and units (uM) are identical.",
      units              = "umol/L (uM)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Almond 2016 Table 1: phenytoin was assayed at seven final concentrations in culture medium containing 0.1 percent dimethyl sulfoxide (v/v) -- 1, 3, 10, 30, 100, 300 and 1000 uM.",
        "Phenytoin MW = 252.28 g/mol (Supplemental Table 2), so 1 uM = 0.2523 mg/L. The clinical DDI study of Table 4 adjusted the phenytoin dose to maintain plasma concentrations of 10-20 ug/mL, i.e. roughly 40-79 uM, which straddles the fitted activity IndC50 of 51.3 uM.",
        "The mRNA IndC50 of 123 uM carries a donor-level S.D. of 120 uM (CV about 98 percent), the largest relative between-donor spread of any parameter in Table 2; simulated mRNA potency for phenytoin is correspondingly wide.",
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
    dose_range       = "Nominal phenytoin concentrations 1, 3, 10, 30, 100, 300 and 1000 uM in culture medium with 0.1 percent dimethyl sulfoxide (v/v) (Table 1).",
    regions          = NA_character_,
    notes            = paste(
      "Cryopreserved human hepatocytes from four donors (Hu1206, Hu1191, Hu1198, Hu4193; Life Technologies) were incubated with serial dilutions of the inducer prepared daily in dimethyl sulfoxide (Materials and Methods).",
      "Two endpoints were measured in parallel in the same incubations: CYP3A4 catalytic activity as 6-beta-hydroxytestosterone formation by liquid chromatography-tandem mass spectrometry, and CYP3A4 mRNA by the QuantiGene Plex 2.0 assay. Because both readouts come from the same cells in the same experiment, they are two outputs of one model rather than two model files.",
      "Curve fitting was carried out on each donor individually and the mean Indmax and IndC50 then calculated across donors; Table 2 reports those means with the donor-level standard deviations, which is what this file's inter-individual variability encodes.",
      "Both a three-parameter (Hill exponent constrained to 1) and a four-parameter sigmoidal model were fitted in GraphPad Prism version 5, and the paper states the two were 'not significantly different', so the three-parameter fit was used. This file therefore carries NO Hill exponent by the authors' own model-selection decision.",
      "Indmax is NOT baseline-corrected: the paper states it 'is equal to Emax + 1'. This is the reason the file carries emax = Indmax - 1.",
      "Phenytoin was one of the two perpetrators (with carbamazepine) for which sufficient in vitro metabolism information was available to simulate the contribution of individual enzymes, so its Simcyp simulations accounted for autoinduction (Materials and Methods, 'PBPK Modeling'). That autoinduction time course belongs to the unreproducible Simcyp layer and does not affect the in vitro parameters carried here.",
      "Supplemental Table 2 additionally reports phenytoin induction parameters for two OTHER enzymes taken from the literature rather than generated in this study -- CYP2C9 Indmax 10.7 / IndC50 9.8 uM (Sahi 2009) and CYP2B6 Indmax 1.9 / IndC50 15.3 uM (Hariparsad 2008). Those are third-party values cited as Simcyp compound-file inputs, not this paper's bench work, and are deliberately NOT carried in this file; only the CYP3A4/5 parameters of Table 2 were generated as part of this study.",
      "In the paper's DDI predictions these in vitro values were used uncalibrated, and calibrated against the rifampicin in vivo reference with Indmax 8 or with the refined Indmax 16 (Eqs. 5 and 6; Table 6)."
    )
  )

  ini({
    # =====================================================================
    # CYP3A4 ACTIVITY endpoint (6-beta-hydroxytestosterone formation).
    # Almond 2016 Table 2, 'Activity' columns, phenytoin row:
    #   Indmax mean 13.6 (S.D. 3.7) fold; IndC50 mean 51.3 (S.D. 29.4) uM.
    # =====================================================================

    lemax_activity <- log(12.6)
    label("Log maximal FRACTIONAL increase in CYP3A4 activity over vehicle control (dimensionless; Emax = Indmax - 1)")
    # Table 2: Indmax (activity) = 13.6 fold, so emax = 13.6 - 1 = 12.6.
    # The paper defines Indmax as the maximum fold induction, "equal to
    # Emax + 1" (Materials and Methods, 'In Vitro Data Analysis').

    lec50_activity <- log(51.3)
    label("Log phenytoin concentration producing half-maximal CYP3A4 activity induction (uM)")
    # Table 2: IndC50 (activity) = 51.3 uM.

    # =====================================================================
    # CYP3A4 mRNA endpoint (QuantiGene Plex 2.0).
    # Almond 2016 Table 2, 'mRNA' columns, phenytoin row:
    #   Indmax mean 24.5 (S.D. 7.6) fold; IndC50 mean 123 (S.D. 120) uM.
    # =====================================================================

    lemax_mrna <- log(23.5)
    label("Log maximal FRACTIONAL increase in CYP3A4 mRNA over vehicle control (dimensionless; Emax = Indmax - 1)")
    # Table 2: Indmax (mRNA) = 24.5 fold, so emax = 24.5 - 1 = 23.5.
    # Narrative check (Results, 'Induction Parameters Determined In Vitro'):
    # mRNA efficacy was 1.3- to 2.0-fold higher than activity efficacy; here
    # 24.5 / 13.6 = 1.80, inside that stated range.

    lec50_mrna <- log(123)
    label("Log phenytoin concentration producing half-maximal CYP3A4 mRNA induction (uM)")
    # Table 2: IndC50 (mRNA) = 123 uM. Narrative check: mRNA potency was
    # 1.0- to 3.3-fold lower (IndC50 higher) than activity potency; here
    # 123 / 51.3 = 2.40, inside the range.

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

    etalemax_activity ~ 0.082714
    # Table 2: S.D. 3.7 on Indmax 13.6 -> CV = 3.7 / 12.6 = 0.2937 on emax;
    # omega^2 = log(1 + 0.2937^2) = 0.082714. This is the TIGHTEST efficacy
    # variance in the study.
    etalec50_activity ~ 0.284008
    # Table 2: S.D. 29.4 on IndC50 51.3 -> CV = 0.5731;
    # omega^2 = log(1 + 0.5731^2) = 0.284008.
    etalemax_mrna ~ 0.099475
    # Table 2: S.D. 7.6 on Indmax 24.5 -> CV = 7.6 / 23.5 = 0.3234 on emax;
    # omega^2 = log(1 + 0.3234^2) = 0.099475.
    etalec50_mrna ~ 0.668759
    # Table 2: S.D. 120 on IndC50 123 -> CV = 0.9756;
    # omega^2 = log(1 + 0.9756^2) = 0.668759. This is the LARGEST variance in
    # the study: the between-donor S.D. is nearly equal to the mean.

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
    foldActivity <- 1 + emax_activity * CP_PHENYTOIN_UM /
      (ec50_activity + CP_PHENYTOIN_UM)
    foldMrna <- 1 + emax_mrna * CP_PHENYTOIN_UM /
      (ec50_mrna + CP_PHENYTOIN_UM)

    foldActivity ~ add(addSd_foldActivity)
    foldMrna ~ add(addSd_foldMrna)
  })
}
