Almond_2016_rifampicin_invivo <- function() {
  description <- paste(
    "Emax concentration-response model for CYP3A4 induction by rifampicin",
    "using the REFINED IN VIVO reference induction parameters that are this",
    "paper's headline result. The model is",
    "fold = 1 + emax * C / (ec50 + C),",
    "the induction term of the Simcyp enzyme-turnover equation (Almond 2016",
    "Eqs. 3 and 4) written on its own, with emax = Indmax - 1 = 15 and",
    "ec50 = IndC50 = 0.32 uM. The driving concentration C is the OPERATIONAL",
    "(unbound, site-of-interaction) rifampicin concentration, not total plasma:",
    "in the liver the paper uses fu,B * I_liver / (Kp / B:P), and in the gut it",
    "uses the enterocyte concentration directly.",
    "The base model (model A) inherited from Simcyp V12 used Indmax 8 with the",
    "same IndC50 0.32 uM, derived from the change in the 6-beta-hydroxycortisol",
    "to cortisol metabolic ratio during multiple-dose rifampicin. That base",
    "model systematically UNDERPREDICTED the interaction when the victim drug",
    "was given orally (geometric mean fold error 2.12 over all victims). The",
    "authors evaluated Indmax of 8-liver/16-gut (model B), 16 (model C), 12",
    "(model F) and 20 (model G), all at IndC50 0.32, and selected model C:",
    "raising Indmax to 16 in BOTH liver and gut gave the lowest overall fold",
    "error (GMFE 1.48 versus 2.12) and the highest proportion of predictions",
    "within the acceptance limits (79.3 percent versus 48.3 percent). This",
    "file carries model C. Both parameters are fixed(): Indmax was selected",
    "from a discrete grid by prediction accuracy rather than estimated by",
    "regression, and IndC50 was inherited unchanged from the base model.",
    "This model is also the CALIBRATION STANDARD of the study: Eqs. 5 and 6",
    "scale another inducer's in vitro Indmax and IndC50 by the ratio of these",
    "in vivo rifampicin values to the in vitro rifampicin values from the same",
    "assay (see Almond_2016_rifampicin_invitro.R).",
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
    "seventh columns of Table 3 from ng/L.h to ng/mL.h; Table 3 carries the",
    "clinical DDI exposures against which model C was selected, so the erratum",
    "is material to how those exposures are read even though no value in this",
    "file is taken from Table 3).",
    "Refined Indmax of 16 from Materials and Methods, 'Design of Virtual",
    "Studies' (model C) and Table 5; base Indmax 8 and IndC50 0.32 uM from",
    "'Design of Virtual Studies' and Supplemental Table 2; the derivation of",
    "the in vivo reference values from Materials and Methods, 'Derivation of",
    "Reference In Vivo Induction Parameters and their Role in Calibration';",
    "fold-induction functional form from Eqs. 3 and 4.",
    sep = " "
  )

  vignette <- "Almond_2016_cyp3a4_induction"

  units <- list(
    time          = "h",
    dosing        = "(none; static concentration-response model driven by an external rifampicin concentration covariate)",
    concentration = "(the observation is dimensionless fold-induction of CYP3A4 over the uninduced baseline; the driving covariate CP_RIF_UM is the operational unbound rifampicin concentration at the site of interaction, in uM)"
  )

  covariateData <- list(
    CP_RIF_UM = list(
      description        = "Operational unbound rifampicin concentration at the site of CYP3A4 interaction, supplied as a covariate. Reused canonical: the existing register entry describes an instantaneous rifampicin plasma concentration used as a perpetrator input, which is the same quantity and unit; here the value required is the UNBOUND concentration at the interaction site rather than total plasma.",
      units              = "umol/L (uM)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "WHICH CONCENTRATION TO SUPPLY. Almond 2016 Eq. 3 drives hepatic induction with fu,B-IN * I_t,Liv / (Kp_IN / B:P_IN) -- the unbound liver concentration corrected for tissue partitioning and blood-to-plasma ratio -- while Eq. 4 drives intestinal induction with the enterocyte concentration I_t,Gut directly. Supplying a TOTAL plasma concentration instead would overstate the driving concentration by roughly 1 / fu; rifampicin fu = 0.15 and B:P = 0.9 (Supplemental Table 2).",
        "SCALE. The fitted IndC50 is 0.32 uM. The Discussion reports that changing rifampicin fugut from 0.19 to 1 raised the simulated unbound portal-vein concentrations, but that 'in both cases the free concentrations exceeded the IndC50 for rifampicin (0.32 uM) across most of the dosing interval; hence, little effect on predictions was observed'. At clinical 600 mg daily dosing the driving concentration therefore sits ABOVE the IndC50 for most of the interval, i.e. on the plateau of this curve.",
        "Rifampicin MW = 823 g/mol (Supplemental Table 2), so 1 uM = 0.823 mg/L = 823 ng/mL.",
        "Set to 0 outside the rifampicin dosing window, at which the model returns fold = 1 (no induction) by construction."
      ),
      source_name        = "It (perpetrator concentration at time t in either the liver or the gut)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 357L,
    n_studies        = 29L,
    age_range        = "matched per simulated trial to the age range reported by each source clinical study",
    weight_range     = "not reported (Simcyp virtual populations generate weight from the built-in demographic distributions)",
    sex_female_pct   = NA_real_,
    race_ethnicity   = "White. The clinical DDI literature search was restricted to studies in white subjects; one otherwise eligible midazolam study was excluded because only one-third of its mixed-ethnicity cohort was white and the data were not stratified.",
    disease_state    = "healthy volunteers, except one simvastatin study in patients with cerebrotendinous xanthomatosis (Table 3)",
    dose_range       = "Rifampicin 300 mg twice daily or 450-600 mg once daily for 4-28 days. Victim drugs: midazolam 0.05 mg/kg or 1-2 mg i.v. and 2-15 mg orally; alfentanil 0.015 mg/kg or 1 mg i.v. and 0.06 mg/kg or 4 mg orally; nifedipine 0.02 mg/kg i.v. and 20 mg orally; alprazolam 1 mg, simvastatin 40 mg, zolpidem 20 mg and triazolam 0.5 mg orally (Table 3).",
    regions          = "not reported",
    notes            = paste(
      "DERIVATION OF THE REFERENCE VALUES. The in vivo concentration-induction response for rifampicin was derived from a study reporting the change in the metabolic ratio of 6-beta-hydroxycortisol to cortisol after multiple dosing of rifampicin 600 mg daily for 14 days, combined with separately published rifampicin concentration-time data (Materials and Methods, 'Derivation of Reference In Vivo Induction Parameters'). Because the endogenous cortisol metabolic ratio behaves like a ratio computed after intravenous administration, the paper notes it 'may not provide information on changes in gut metabolism' -- which is precisely the deficiency the refinement to Indmax 16 corrects.",
      "SELECTION OF THE REFINED VALUE. n_studies is the 29 rifampicin DDI studies of Table 3 (10 with the victim given intravenously, summing to 116 subjects, and 19 orally, summing to 241). n_subjects is the 357-subject sum of those reported study sizes and is an upper bound on distinct individuals: two rows report n = 52 for what appears to be the same large cohort contributing both an intravenous and an oral midazolam arm, so some subjects are counted twice. Each study was replicated as 10 virtual trials matched to the published age range, sex split, study size, dose, route, timing, frequency and duration for both perpetrator and victim (Materials and Methods, 'Design of Virtual Studies'). Model selection was by geometric mean fold error and by the proportion of predictions inside published acceptance limits (Table 5).",
      "PREDICTION ACCURACY OF THE SHIPPED MODEL C (Table 5, all victim drugs): GMFE 1.48 and 79.3 percent within acceptance limits, versus GMFE 2.12 and 48.3 percent for the base model A. The improvement is concentrated in the oral arm (GMFE 1.69 versus 2.81), because the base model's underprediction was an intestinal-induction deficiency; the intravenous arm was already accurate under model A (GMFE 1.24) and improves only slightly (1.15).",
      "INDEPENDENT SUPPORT FOR A HIGHER Indmax. The Discussion notes that other investigators have separately reported a need for rifampicin Indmax values of 11.5, 12.5 and 14.6, 'not dissimilar to the value of 16-fold used here', and that using those values in this model gave comparable prediction accuracy.",
      "SUPPLEMENT LABEL TYPO. Supplemental Table 2 lists the rifampicin base values on two adjacent rows labelled 'CYP3A4 Indmax (Fold, Emax +1) 8' and 'CYP3A5 IndC50 (uM) 0.32'. The second label is a typographical error for CYP3A4: the main text states the base reference values as 'Indmax = 8; IndC50 = 0.32' for CYP3A4, and Materials and Methods states that 'only CYP3A4 was considered as CYP3A5 induction is less well characterized'. The value 0.32 uM is used here as the CYP3A4 IndC50 on the authority of the main text.",
      "WHAT IS NOT CARRIED. Eqs. 3 and 4 also contain a mechanism-based-inhibition term, kinact * I / (KI + I), and Eqs. 1 and 2 contain competitive-inhibition terms. Neither applies to rifampicin as an inducer in this study: no kinact is reported for rifampicin anywhere on disk, and its reported CYP3A4 Ki of 10.5 uM (Supplemental Table 2) belongs to the competitive-inhibition term of Eqs. 1 and 2, which acts on the victim's clearance rather than on the enzyme pool. Carrying an inhibition term here would require constants the paper does not report."
    )
  )

  ini({
    # =====================================================================
    # Refined in vivo reference induction parameters for rifampicin
    # (model C of Almond 2016).
    #
    # Both values are fixed(), not estimated. Indmax was selected from the
    # discrete set {8, 12, 16, 20} by comparing DDI prediction accuracy
    # across 29 clinical studies (Materials and Methods, 'Design of Virtual
    # Studies'; Table 5); IndC50 was inherited unchanged from the base
    # model in every alternative the authors tested. Neither carries a
    # standard error or confidence interval in the source, which is itself
    # the signal that they are held constants rather than regression
    # estimates.
    # =====================================================================

    lemax <- fixed(log(15))
    label("Log maximal FRACTIONAL increase in CYP3A4 activity over the uninduced baseline (dimensionless; Emax = Indmax - 1)")
    # Materials and Methods, 'Design of Virtual Studies', model C: "Use of a
    # higher Indmax in both the gut and liver (16) but the same IndC50
    # (0.32)". Table 5 column 'Model C' confirms "Indmax 16, IndC50 0.32"
    # and reports it as the selected model. The paper defines Indmax as the
    # maximum fold induction, "equal to Emax + 1" (Materials and Methods,
    # 'In Vitro Data Analysis'), so emax = 16 - 1 = 15.
    #
    # The `emax` / `ec50` naming (rather than a new `indmax` / `indc50`
    # canonical) follows the established translation in
    # Willemin_2024_interleukin6_cyp_pbpk.R.
    #
    # Alternatives the authors evaluated, all at IndC50 0.32 uM, with the
    # all-victim GMFE from Table 5 (lower is better):
    #   model A  Indmax  8 liver and gut  GMFE 2.12  (Simcyp V12 base model)
    #   model B  Indmax  8 liver, 16 gut  GMFE 1.77
    #   model C  Indmax 16 liver and gut  GMFE 1.48  <- shipped here
    #   model F  Indmax 12 liver and gut  GMFE 1.63
    #   model G  Indmax 20 liver and gut  GMFE 1.51
    # Model C is not merely the lowest GMFE; it also ties model F for the
    # highest proportion within acceptance limits (79.3 percent).
    #
    # Model B is the one alternative this single-parameter file cannot
    # express, because it uses DIFFERENT Indmax values in liver (8) and gut
    # (16). Representing it would require two enzyme pools, which in turn
    # requires the tissue-specific kdeg values the paper never reports.

    lec50 <- fixed(log(0.32))
    label("Log operational unbound rifampicin concentration producing half-maximal CYP3A4 induction (uM)")
    # Materials and Methods, 'Design of Virtual Studies': the in vivo
    # reference values are "Indmax = 8; IndC50 = 0.32" for rifampicin, and
    # every refined model B-G retains "the same IndC50 (0.32)". Also in
    # Supplemental Table 2 (see the label-typo note in population$notes) and
    # quoted again in the Discussion as "the IndC50 for rifampicin (0.32
    # uM)".
    #
    # Derived from the 6-beta-hydroxycortisol / cortisol metabolic ratio
    # during rifampicin 600 mg daily for 14 days, combined with published
    # rifampicin concentration-time data (Materials and Methods,
    # 'Derivation of Reference In Vivo Induction Parameters'). The paper
    # flags that these two inputs come from separate studies and that
    # rifampicin pharmacokinetic variability means the plasma
    # concentrations may have differed between them -- the first of the
    # candidate explanations it offers for the base model's underprediction.
    #
    # Cross-check against the paper's own in vitro work: the in vitro
    # rifampicin activity IndC50 is 0.30 uM (Table 2), within 7 percent of
    # this in vivo value. The efficacy parameters do NOT agree so closely
    # (in vitro Indmax 22.7 activity / 29.9 mRNA against in vivo 16), which
    # is the entire motivation for the calibration approach of Eqs. 5 and 6.

    # =====================================================================
    # No inter-individual variability. This is a single pair of population
    # reference values selected against pooled clinical data; the source
    # reports no between-subject distribution for either parameter. The
    # between-subject variability in the paper's own simulations comes from
    # the Simcyp virtual population's demographic and physiological
    # sampling, which is not reproducible from the published inputs.
    #
    # No residual error is reported either -- the model was assessed by
    # geometric mean fold error against study-level AUC ratios, not by a
    # residual-error model. Per the standing policy on unreported residual
    # error, addSd is fixed at 0 so the model returns the deterministic
    # curve. Both omissions are flagged in the vignette Errata.
    # =====================================================================
    addSd <- fixed(0)
    label("Additive residual error on fold CYP3A4 induction (ZERO - not reported in source)")
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    # ===================================================================
    # Fold induction of CYP3A4 relative to the uninduced baseline:
    #
    #   fold(C) = 1 + emax * C / (ec50 + C)
    #
    # This is the induction term of the paper's enzyme-turnover equations
    # (Eqs. 3 and 4), which read 1 + (Indmax - 1) * I / (IndC50 + I).
    #
    # The full equations relax the active-enzyme pool toward this fold with
    # a time constant of 1 / kdeg:
    #
    #   d(Enz)/dt = kdeg * Enz0 * fold(C) - Enz * (kdeg + kinact * C / (KI + C))
    #
    # so fold(C) is the STEADY-STATE induction ratio reached under a
    # sustained concentration C in the absence of mechanism-based
    # inhibition. The dynamic form is NOT carried here because the paper
    # never reports kdegH-3A4 or kdegG-3A4; it reports only that a
    # sensitivity analysis was run on them. Users who need the time course
    # can supply their own kdeg and wrap this term, following the pattern
    # of Chen_2024_interleukin6_cyp3a_pbpk.R.
    #
    # Structural checks asserted in the validation vignette:
    #   fold(0)    = 1                (no rifampicin, uninduced baseline)
    #   fold(ec50) = 1 + emax / 2 = 8.5
    #   fold(Inf)  = 1 + emax = 16    (the refined Indmax of model C)
    # ===================================================================
    foldCyp3a4 <- 1 + emax * CP_RIF_UM / (ec50 + CP_RIF_UM)

    foldCyp3a4 ~ add(addSd)
  })
}
