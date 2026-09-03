Chetty_2014_efalizumab_cd11a_qsp <- function() {
  description <- paste(
    "QSP. CD11a target-turnover and target-engagement model for efalizumab, a",
    "humanised anti-CD11a monoclonal antibody once used for moderate-to-severe",
    "plaque psoriasis. CD11a (the alpha subunit of LFA-1, expressed on T",
    "lymphocytes) is synthesised at a zero-order rate and lost by constitutive",
    "degradation; efalizumab bound to CD11a forms a complex that is",
    "internalised and degraded faster than free CD11a, so drug exposure",
    "down-modulates the total CD11a pool and, far more steeply, the free",
    "(unoccupied) pool. The binding is written in the Michaelis-Menten /",
    "quasi-steady-state approximation of target-mediated drug disposition of",
    "Gibiansky et al. 2008, in which a single constant Km subsumes kon, koff",
    "and kint, so the drug-target complex is an algebraic function of the total",
    "target and the driving concentration rather than a separate ODE state.",
    "The efalizumab concentration is NOT modelled here. Chetty 2014 generated",
    "it with the Simcyp Simulator (Version 13 R1) Mechanistic FcRn whole-body",
    "PBPK model, a platform port whose tissue vascular / endothelial /",
    "interstitial volumes, plasma and lymph flows, reflection coefficients,",
    "endothelial uptake and recycling rate constants and FcRn abundances are",
    "Simcyp database outputs: the paper publishes no organ ODEs and no volume",
    "term of any kind (its only disposition parameter is a systemic clearance",
    "of 0.0227 L/h), so that layer is not reproducible from the on-disk sources",
    "and is deliberately NOT extracted. The concentration is instead supplied",
    "per record as the canonical time-varying covariate CP_EFALIZUMAB_UGML and",
    "converted to umol/L inside model(), following the Zhang_2026_ribociclib_qsp",
    "and Liang_2024_osimertinib_qsp precedents.",
    "KNOWN DEVIATION, quantified in the vignette: driven by a PLASMA efalizumab",
    "concentration this model predicts far deeper free-CD11a suppression than",
    "the paper's own Figure 3 (about 0.1 to 3 percent of baseline across the",
    "trough-to-peak range of the paper's Figure 2C, against 9 to 24 percent",
    "predicted in Figure 3). The paper never states which matrix drives target",
    "binding; inverting the model on the Figure 3 values implies a driver of",
    "roughly 0.05 to 0.16 ug/mL, one to two orders of magnitude below the",
    "plasma concentrations of the paper's own Figure 2C, i.e. an interstitial or",
    "tissue concentration internal to the un-extracted PBPK layer. The paper's",
    "Introduction advertises exactly that capability: 'input to the PD model can",
    "be done from a tissue interstitial compartment and not just from plasma,",
    "which is important when modeling membrane bound receptors'.",
    "Users must therefore choose the matrix their",
    "CP_EFALIZUMAB_UGML trajectory represents. The TOTAL CD11a pool is far less",
    "sensitive to that choice: it saturates at kdeg / kint = 18.5 percent of",
    "baseline for any driving concentration well above Km, against the roughly",
    "25 percent clinical down-modulation plateau the paper cites from Ng 2005.",
    sep = " "
  )
  reference <- paste(
    "Chetty M, Li L, Rose R, Machavaram K, Jamei M, Rostami-Hodjegan A,",
    "Gardner I. (2015). Prediction of the pharmacokinetics, pharmacodynamics,",
    "and efficacy of a monoclonal antibody, using a physiologically based",
    "pharmacokinetic FcRn model. Front Immunol 5:670.",
    "doi:10.3389/fimmu.2014.00670. PMCID PMC4283607.",
    "Every parameter value is from main-text Table 1. The Michaelis-Menten /",
    "quasi-steady-state target equations are the form Chetty 2014 attributes to",
    "its reference 23 (Gibiansky L, Gibiansky E, Kakkar T, Ma P. Approximations",
    "of the target-mediated drug disposition model and identifiability of model",
    "parameters. J Pharmacokinet Pharmacodyn 2008;35:573-591), and the paper",
    "enumerates their ingredients explicitly in the Methods bullet list (kon",
    "binding, koff dissociation, kint complex internalisation, saturation of",
    "binding sites, Ksyn synthesis and Kdeg degradation of CD11a, with Km",
    "subsuming kon, koff and kint). Chetty 2014 has no supplementary material",
    "(EuropePMC reports hasSuppl = N) and no erratum; the only linked item is a",
    "later 'Comment in' (AAPS J 2016;18:948-59).",
    "The companion efficacy model from the same paper is",
    "modellib('Chetty_2014_efalizumab_pasi').",
    sep = " "
  )
  vignette <- "Chetty_2014_efalizumab_cd11a_pasi"

  units <- list(
    time          = "h",
    dosing        = "(not applicable; the efalizumab concentration is supplied as the time-varying covariate CP_EFALIZUMAB_UGML)",
    concentration = "umol/L"
  )

  compartmentData <- list(
    # Gottlieb 2002 (Chetty 2014 reference 15) and Ng 2005 (reference 18)
    # measured CD11a expression on peripheral-blood T lymphocytes; Chetty 2014
    # Table 1 carries the pool as a molar concentration ("Rmax: CD11a abundance
    # 0.01 uM"), not as a per-cell receptor count.
    total_target = list(
      analyte  = "CD11a (total: free + efalizumab-bound)",
      units    = "umol/L",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  covariateData <- list(
    CP_EFALIZUMAB_UGML = list(
      description        = "Time-varying efalizumab concentration driving CD11a target engagement",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Converted to umol/L inside model() as CP_EFALIZUMAB_UGML / mw_efa *",
        "1000 (equivalently, divided by 148.841), because Km and the CD11a pool",
        "are reported in umol/L.",
        "Under the Michaelis-Menten / quasi-steady-state approximation the",
        "concentration that enters the binding term is the FREE (target-unbound)",
        "drug concentration. Efalizumab assays report total drug; the two differ",
        "by at most the bound-target concentration, which cannot exceed",
        "Rmax = 0.01 umol/L (1.49 ug/mL). The distinction is therefore",
        "immaterial above roughly 15 ug/mL and material at the sub-ug/mL",
        "concentrations that dominate free-CD11a suppression.",
        "WHICH MATRIX: Chetty 2014 never states whether its Simcyp",
        "implementation drove CD11a binding from plasma or from a tissue",
        "interstitial concentration; its Introduction notes only that 'input to",
        "the PD model can be done from a tissue interstitial compartment and not",
        "just from plasma, which is important when modeling membrane bound",
        "receptors'. Inverting this model on the paper's own Figure 3 predicted",
        "free-CD11a percentages implies a driver of roughly 0.05 to 0.16 ug/mL",
        "over the treatment period, one to two orders of magnitude below the",
        "plasma concentrations of its Figure 2C. See the vignette 'Assumptions",
        "and deviations' section.",
        "Set to 0 for untreated / pre-dose periods and for placebo arms, which",
        "returns total and free CD11a to Rmax.",
        "Reference values observed (Chetty 2014 Figure 2C, PREDICTED",
        "concentration-time profile for the Gottlieb escalating regimen, read",
        "off the panel): weekly peaks rising from about 5 ug/mL in week 1 to",
        "about 30 ug/mL in weeks 5 to 7, with pre-dose troughs rising from about",
        "0.5 ug/mL to about 7 ug/mL. Single-dose profiles for 1, 3 and 10 mg/kg",
        "are in Figure 2A."
      ),
      source_name        = "(none; computed by the Simcyp Mechanistic FcRn PBPK model, not a named data column)"
    )
  )

  population <- list(
    species        = "human (in silico; Simcyp virtual north European Caucasian healthy volunteers, verified against published psoriasis cohorts)",
    n_subjects     = 500L,
    n_studies      = 4L,
    age_range      = "25-50 years (simulated cohort)",
    sex_female_pct = 50,
    disease_state  = "Simulated in healthy volunteers; verified against clinical data from adults with moderate-to-severe plaque psoriasis",
    dose_range     = paste(
      "Single intravenous doses of 1, 3 and 10 mg/kg (Bauer 1999, Chetty 2014",
      "reference 12); multiple dosing escalating 0.3 mg/kg/week in week 1,",
      "0.4 in week 2, 0.6 in week 3 and 1 mg/kg/week for the following 4 weeks,",
      "each given as a 1 h infusion (Gottlieb 2002, reference 15); and",
      "1 mg/kg/week (Ng 2005, reference 18)"
    ),
    regions        = "north European Caucasian virtual population",
    notes          = paste(
      "Chetty 2014 Simulations section: 'Predictive studies used 5 trials with",
      "100 virtual north European Caucasian Healthy Volunteers each, aged",
      "between 25 and 50 years, with an equal proportion of males and females',",
      "hence n_subjects = 500 across 5 trials. Simulations replicating a",
      "clinical study instead used trial designs as close as possible to that",
      "study's design. n_studies = 4 counts the clinical datasets the PK and",
      "CD11a predictions were verified against: Bauer 1999 (single doses),",
      "Gottlieb 2000, Gottlieb 2002 (the escalating multiple-dose regimen and",
      "the CD11a time course of Figure 3) and Ng 2005.",
      "The paper cites its escalating-regimen study as reference 14 in the",
      "Methods and as reference 15 in the Results; the Figure 3 legend keys the",
      "observed series to 'Gottlieb 2002', which is reference 15",
      "(Arch Dermatol 138:591-600), so reference 15 is the escalating",
      "multiple-dose study and reference 14 (Gottlieb 2000, J Am Acad Dermatol",
      "42:428-435) is the single-dose study.",
      "No per-subject demographic table is published; the model carries no",
      "demographic covariates, so this metadata is descriptive only.",
      "Observed CD11a down-modulation reported for the class: 'typically to",
      "about 25% of baseline' (Chetty 2014 Discussion, citing Ng 2005), with",
      "'significant inter-individual variability'."
    )
  )

  ini({
    # ========================================================================
    # Molecular weight -- used only to convert the ug/mL covariate onto the
    # umol/L scale of Km and the CD11a pool.
    # ========================================================================
    mw_efa <- fixed(148841)
    label("Molecular weight of efalizumab (g/mol)")
    # Chetty 2014 Table 1, row "MW: molecular weight of efalizumab",
    # 148841 g/mol (source: reference 19, Boehncke 2007).

    # ========================================================================
    # CD11a pool size and turnover.
    # ========================================================================
    rbase <- fixed(0.01)
    label("Baseline (drug-free) total CD11a concentration, the paper's Rmax (umol/L)")
    # Chetty 2014 Table 1, row "Rmax: CD11a abundance", 0.01 uM, listed as
    # "Estimated". This is also the drug-free steady state of the ODE below,
    # because Table 1 defines Ksyn = Rmax * Kdeg, so ksyn / kdeg = rbase
    # identically and total_target(0) <- rbase starts the system at rest.

    lkdeg <- fixed(log(0.0185))
    label("First-order degradation rate constant of free CD11a (1/h)")
    # Chetty 2014 Table 1, row "Kdeg: degradation rate of the target ie CD11a",
    # 0.0185 1/h with CV% = 10 (source: reference 4, Chen and Balthasar 2012).
    # Table 1 glosses Kdeg as a "zero-order degradation rate" in the Figure 1
    # legend, but the value carries units of 1/h and Ksyn = Rmax * Kdeg only
    # balances if Kdeg multiplies a concentration, so it is first-order; the
    # Figure 1 legend also transposes the same way for Ksyn, which it calls
    # "first-order" while Table 1 gives it in uM/h. Table 1 is authoritative.

    lkint <- fixed(log(0.1))
    label("First-order internalisation / degradation rate constant of the efalizumab-CD11a complex (1/h)")
    # Chetty 2014 Table 1, row "Kint: internalization rate constant for the
    # complex", 0.1 (source: reference 18, Ng 2005). Table 1 prints the unit as
    # "l/h", which is dimensionally impossible for a rate constant acting on a
    # concentration; the Figure 1 legend gives the correct reading,
    # "kint - internalization rate constant of complex", i.e. 1/h. Read as
    # 0.1 1/h. This is the only unit correction applied anywhere in this file.

    lkm <- fixed(log(0.000573))
    label("Michaelis-Menten / quasi-steady-state constant for efalizumab-CD11a binding (umol/L)")
    # Chetty 2014 Table 1, row "Km: rate constant for receptor complex
    # internalization and degradation", 0.000573 uM (source: reference 20,
    # Peletier and Gabrielsson 2009). Table 1 calls it a "rate constant" but
    # prints it in uM and the Methods define it correctly: "Km (rate constant
    # for receptor complex internalization and degradation) is used in the MM
    # model and incorporates kon, koff, and kint", i.e. Km = (koff + kint)/kon,
    # a concentration. Read as 0.000573 umol/L (0.573 nmol/L), consistent with
    # efalizumab's reported low-nanomolar CD11a affinity.

    # ========================================================================
    # Between-subject variability. Chetty 2014 Table 1 attaches a CV% to
    # exactly two compound-file entries: CLiv (30%), which belongs to the
    # un-extracted PBPK layer, and Kdeg (10%), which belongs here. Encoded as
    # a log-normal variance, omega^2 = log(1 + CV^2), the library convention
    # for a reported CV% on a strictly positive rate constant. Chetty 2014
    # does not state which distribution family the Simcyp compound file used.
    # ========================================================================
    etalkdeg ~ log(1 + 0.10^2)

    # ========================================================================
    # Residual error. Chetty 2014 is a simulation study: it reports no
    # residual-error model and no assay variability for CD11a, so this is
    # encoded as zero rather than invented. See the vignette Errata.
    # ========================================================================
    propSd_total_target <- fixed(0)
    label("Proportional residual standard deviation for total CD11a (fraction)")
  })

  model({
    kdeg <- exp(lkdeg + etalkdeg)
    kint <- exp(lkint)
    km   <- exp(lkm)

    # Chetty 2014 Table 1, row "Ksyn: rate of synthesis of target", is given as
    # the derived quantity "Ksyn = Rmax * Kdeg" and evaluates to
    # 0.01 * 0.0185 = 0.000185 uM/h, exactly the printed value. Deriving it
    # here rather than hard-coding it keeps the drug-free steady state equal to
    # rbase for every subject, including under the Kdeg between-subject
    # variability above.
    ksyn <- rbase * kdeg

    # =====================================================================
    # 1. Driving efalizumab concentration, ug/mL -> umol/L.
    #        C[umol/L] = C[ug/mL] / MW[g/mol] * 1000
    #    With MW = 148841 g/mol this is C[ug/mL] / 148.841.
    # =====================================================================
    cefa <- CP_EFALIZUMAB_UGML / mw_efa * 1000

    # =====================================================================
    # 2. Quasi-steady-state partition of the total CD11a pool between free
    #    and efalizumab-bound receptor. Setting the binding flux equal to the
    #    complex-loss flux, kon * C * R = (koff + kint) * RC, gives
    #    RC / R = C / Km with Km = (koff + kint) / kon, hence
    #        RC = Rtot * C / (Km + C)
    #        R  = Rtot * Km / (Km + C)
    #    This is the Michaelis-Menten / quasi-steady-state approximation of
    #    Gibiansky et al. 2008 that Chetty 2014 Methods names (its reference
    #    23), with the paper's own list of ingredients: kon binding, koff
    #    dissociation, kint internalisation, saturation of binding sites, and
    #    Km subsuming all three.
    # =====================================================================
    foccCd11a  <- cefa / (km + cefa)
    boundCd11a <- total_target * foccCd11a
    freeCd11a  <- total_target - boundCd11a

    # =====================================================================
    # 3. Total CD11a turnover. Chetty 2014 Methods: the model accounted for
    #    "Changes in CD11a concentration due to synthesis (Ksyn) and
    #    degradation (Kdeg)" alongside "Internalization or degradation of the
    #    complex (kint)". Free receptor is lost at kdeg and bound receptor at
    #    kint, so the total pool obeys
    #        d(Rtot)/dt = ksyn - kdeg * R - kint * RC
    #    Because kint (0.1 /h) exceeds kdeg (0.0185 /h), drug exposure
    #    down-modulates the total pool; its saturating floor as C >> Km is
    #    ksyn / kint = rbase * kdeg / kint = 18.5% of baseline.
    # =====================================================================
    d/dt(total_target) <- ksyn - kdeg * freeCd11a - kint * boundCd11a

    total_target(0) <- rbase

    # =====================================================================
    # 4. Readouts. Chetty 2014 Figure 3 plots "%CD11a", the free (unoccupied)
    #    CD11a concentration as a percentage of its time-zero baseline;
    #    pctTotalCd11a is the companion for the total pool, which is what a
    #    non-competing-antibody flow-cytometry assay measures.
    # =====================================================================
    pctFreeCd11a  <- 100 * freeCd11a / rbase
    pctTotalCd11a <- 100 * total_target / rbase

    # Mass-balance check quantity: equals total_target at all times.
    sumCd11a <- freeCd11a + boundCd11a

    # The endpoint is the modelled ODE state, the total CD11a pool, following
    # the Liang_2024_rituximab precedent for a QSS-TMDD target state. The
    # paper's Figure 3 readout, free CD11a as a percentage of baseline, is the
    # derived pctFreeCd11a column above.
    total_target ~ prop(propSd_total_target)
  })
}
