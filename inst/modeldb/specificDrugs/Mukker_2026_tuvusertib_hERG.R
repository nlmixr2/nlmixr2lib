Mukker_2026_tuvusertib_hERG <- function() {
  description <- paste(
    "In vitro (HEK-293). Sigmoidal Hill concentration-response model for",
    "inhibition of the hERG potassium channel by the investigational ATR",
    "kinase inhibitor tuvusertib (M1774), measured by GLP-compliant",
    "whole-cell patch clamp in HEK-293 cells stably expressing hERG.",
    "The model is:",
    "hergInh = imax * C^hill / (ic50^hill + C^hill),",
    "with ic50 = 1048 ng/mL (2.83 uM) and hill = 1.12, both fitted to",
    "tail-current block at nominal tuvusertib concentrations of 0.3, 1,",
    "3 and 10 uM (goodness of fit R^2 = 0.98). imax is fixed at 1",
    "(complete block at saturating concentration), the conventional",
    "three-parameter reduction of the Hill equation for a normalised",
    "fractional-block patch-clamp curve; the source reports only IC50",
    "and the Hill coefficient. The output is the fraction of hERG tail",
    "current blocked (0-1), so the model carries no inter-individual",
    "variability and no residual error -- it is a single fitted",
    "concentration-response curve, not a population model. This is the",
    "nonclinical arm of the paper's integrated QTc risk assessment: the",
    "IC50 sits ~2.5-fold above the clinical unbound steady-state Cmax",
    "of 415.95 ng/mL at the 180 mg QD recommended dose for expansion.",
    "Companion model files Mukker_2026_tuvusertib_QTcF.R and",
    "Mukker_2026_tuvusertib_HR.R carry the clinical concentration-QTc",
    "and concentration-heart-rate analyses.",
    sep = " "
  )

  reference <- paste(
    "Mukker JK, Yap TA, Tolcher AW, de Bono JS, Plummer R, Grosser G,",
    "van Amsterdam C, Schieferstein H, Witjes H, Diderichsen PM,",
    "Krebs-Brown A, Gao W, Strotmann R, Szucs Z, Gounaris I,",
    "Venkatakrishnan K. An Integrated Nonclinical and Clinical Risk",
    "Assessment of the Effects of Investigational ATRi Tuvusertib on QTc",
    "Interval in Patients With Solid Tumors.",
    "Clinical and Translational Science 2026;19(2):e70496.",
    "doi:10.1111/cts.70496.",
    "Parameter values are in Methods 2.1 and Results 3.1.",
    sep = " "
  )

  vignette <- "Mukker_2026_tuvusertib_QTc"

  units <- list(
    time          = "h",
    dosing        = "(none; static in vitro concentration-response model fed by an external tuvusertib concentration covariate)",
    concentration = "(observation hergInh is the fraction of hERG tail current blocked, dimensionless 0-1; driving covariate CP_TUVUSERTIB_NGML is the nominal bath tuvusertib concentration in ng/mL)"
  )

  covariateData <- list(
    CP_TUVUSERTIB_NGML = list(
      description        = "Nominal tuvusertib concentration in the patch-clamp bath, supplied as a covariate. Reused canonical: in this in-vitro model the column carries a bath concentration rather than a plasma concentration, but the quantity and units are identical.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The assay evaluated nominal tuvusertib concentrations of 0.3, 1, 3 and 10 uM (Mukker 2026 Methods 2.1 and Results 3.1). Converting at the molecular weight implied by the paper's own dual-unit IC50 report (1048 ng/mL / 2.83 uM = 370.3 g/mol) gives 111, 370, 1111 and 3703 ng/mL.",
        "The paper reports the IC50 in both units (2.83 uM and 1048 ng/mL); this file uses the ng/mL form so that the in-vitro potency is directly comparable with the clinical plasma concentrations driving the companion C-QTc models, which are in ng/mL.",
        "The clinically relevant comparator is the UNBOUND steady-state Cmax at 180 mg QD, 415.95 ng/mL (Mukker 2026 Results 3.1), giving the ~2.5-fold exposure margin quoted in the Abstract and Discussion. Comparing against a TOTAL plasma concentration would be a units error -- the patch-clamp bath concentration is a free concentration.",
        "The 10 uM top concentration blocked more than 30% of the tail current, which is what permitted a concentration-response curve to be fitted and an IC50 determined (Results 3.1)."
      ),
      source_name        = "nominal tuvusertib concentration"
    )
  )

  population <- list(
    species          = "in vitro (HEK-293 cells stably expressing the hERG potassium channel)",
    n_subjects       = NA_integer_,
    n_studies        = 1L,
    age_range        = NA_character_,
    weight_range     = NA_character_,
    sex_female_pct   = NA_real_,
    race_ethnicity   = NA_character_,
    disease_state    = "Not applicable -- recombinant cell line.",
    dose_range       = "Nominal tuvusertib bath concentrations 0.3, 1, 3 and 10 uM (approximately 111, 370, 1111 and 3703 ng/mL).",
    regions          = NA_character_,
    notes            = paste(
      "Good Laboratory Practice (GLP)-compliant whole-cell patch-clamp assay of the hERG potassium channel stably expressed in human embryonic kidney (HEK)-293 cells (Mukker 2026 Methods 2.1). Inhibition of hERG is the principal pharmacodynamic mechanism of ventricular action-potential / QT prolongation in humans.",
      "Assay validity: the test system was validated with the reference item E4031, a selective blocker of the cardiac rapid delayed-rectifier potassium current, which effectively blocked the hERG tail current (Results 3.1).",
      "Fit quality: R^2 = 0.98 for the concentration-response curve (Results 3.1).",
      "Study timing caveat carried by the authors: this hERG study and the companion dog telemetry study were designed and conducted BEFORE the 2022 ICH E14/S7B Q&A update and may not meet all of its current expectations, including its specific requirements for in-vitro study design and reference compounds (Mukker 2026 Introduction and Discussion).",
      "The paper's in vivo nonclinical arm (Beagle dog cardiovascular safety: a single 3 mg/kg oral dose in 4 males, and a 4-week repeat-dose study at 1 / 2.5 / 5 mg/kg/day) is DESCRIPTIVE only -- no compartmental or concentration-response model was fitted to the dog data, so it is not extractable and is not represented by a model file. Its quantitative results are reproduced as context in the validation vignette.",
      "Context for interpreting the margin: a ~2.5-fold free-Cmax / hERG-IC50 ratio is modest against the traditional 30-fold margin sought for high confidence of low proarrhythmic risk, which the authors acknowledge explicitly; the clinical C-QTc analysis is what carries the risk conclusion (Mukker 2026 Discussion)."
    )
  )

  ini({
    # ==================================================================
    # Sigmoidal Hill inhibition model for hERG tail-current block.
    # Mukker 2026 Results 3.1: 'The IC50 was determined to be 2.83 uM
    # (1048 ng/mL), with a Hill coefficient of 1.12. The goodness-of-fit
    # was evaluated using an R^2 value of 0.98 for the concentration-
    # response curve.'
    # ==================================================================

    lic50 <- log(1048)
    label("Log half-maximal inhibitory concentration for hERG tail current (ng/mL)")
    # Mukker 2026 Results 3.1: IC50 = 2.83 uM = 1048 ng/mL. The ng/mL
    # form is used so the potency is directly comparable with the
    # clinical plasma concentrations in the companion C-QTc models.
    # Internal consistency check on the paper's own dual-unit report:
    # 1048 / 2.83 = 370.3 g/mol, which is the molecular weight of
    # tuvusertib -- the two reported units are mutually consistent.

    lhill <- log(1.12)
    label("Log Hill coefficient of the hERG concentration-response curve (dimensionless)")
    # Mukker 2026 Results 3.1: Hill coefficient = 1.12. The paper calls
    # it 'Hill coefficient'; the canonical ini() name is `lhill` and the
    # bare in-model name `hill` (see parameter-names.md, 'Sigmoidal PD
    # shape parameters'). Close to 1, i.e. nearly hyperbolic binding.

    limax <- fixed(log(1))
    label("Log maximal fractional hERG block (1 = complete block; see notes)")
    # NOT reported in the source. Fixed at 1 (log(1) = 0), the standard
    # three-parameter reduction of the Hill equation used for a
    # normalised fractional-block patch-clamp curve, in which the IC50
    # is by definition the concentration producing half of COMPLETE
    # block. The paper's language supports this reading -- it reports a
    # 'half-maximal inhibitory concentration (IC50)' from a
    # concentration-response curve fitted once the top concentration
    # exceeded 30% block -- but the assumption is not stated
    # explicitly, so it is flagged in the vignette Errata. If a future
    # source reports a fitted Imax below 1, both this value and the
    # IC50 interpretation would need revisiting together.

    # ==================================================================
    # No inter-individual variability and no residual error. This is a
    # single fitted concentration-response curve from one GLP assay,
    # not a population model: the source reports point estimates and an
    # R^2 only, with no replicate-level variance components. Per the
    # standing policy on unreported residual error, addSd is fixed at
    # 0 so the model returns the deterministic curve.
    # ==================================================================
    addSd <- fixed(0)
    label("Additive residual error on fractional hERG block (ZERO - not reported in source)")
  })

  model({
    ic50 <- exp(lic50)
    hill <- exp(lhill)
    imax <- exp(limax)

    # ==================================================================
    # Fraction of hERG tail current blocked (0-1) at the bath
    # concentration CP_TUVUSERTIB_NGML. Standard sigmoidal Imax form
    # (parameter-names.md, 'Sigmoidal PD shape parameters'):
    #
    #   hergInh = imax * C^hill / (ic50^hill + C^hill)
    #
    # At C = ic50 this returns imax / 2 = 0.5, which is the definition
    # of the reported IC50 and is asserted in the validation vignette.
    # ==================================================================
    hergInh <- imax * CP_TUVUSERTIB_NGML^hill /
      (ic50^hill + CP_TUVUSERTIB_NGML^hill)

    hergInh ~ add(addSd)
  })
}
