Agarwal_2011_sorafenib_bcrp <- function() {
  description <- paste(
    "In vitro (MDCKII-Bcrp1). Inhibitory Emax model for the saturable,",
    "BCRP-mediated basolateral-to-apical (B-to-A) transcellular",
    "permeability of sorafenib across polarised MDCKII monolayers stably",
    "transfected with murine Bcrp1 (ABCG2). The model is:",
    "pappBa = emax * (1 - C / (km + C)),",
    "which is algebraically the hyperbolic decay pappBa = emax * km / (km + C).",
    "As the donor (basolateral) sorafenib concentration rises the",
    "transporter saturates, so the active efflux component of the B-to-A",
    "permeability falls from its tracer-concentration maximum toward zero.",
    "km = 3.6 ng/mL (5.5 nM) is the reported apparent affinity constant of",
    "BCRP for sorafenib and is the paper's headline in-vitro finding --",
    "sorafenib is a HIGH-affinity BCRP substrate, which is what makes BCRP",
    "rather than P-gp the dominant transporter restricting its entry into",
    "the brain. emax is NOT reported in the source and was recovered by",
    "digitising Figure 5 (see the parameter comment and the vignette",
    "Errata). The output is an apparent permeability in units of",
    "1e-6 cm/s, so the model carries no inter-individual variability and",
    "no residual error -- it is a single fitted concentration-response",
    "curve from one cell-culture experiment, not a population model.",
    "Companion model file Agarwal_2011_sorafenib_pgp.R carries the",
    "sigmoidal P-gp inhibition curve from the same paper.",
    sep = " "
  )

  reference <- paste(
    "Agarwal S, Sane R, Ohlfest JR, Elmquist WF.",
    "The role of the breast cancer resistance protein (ABCG2) in the",
    "distribution of sorafenib to the brain.",
    "J Pharmacol Exp Ther. 2011;336(1):223-233.",
    "doi:10.1124/jpet.110.175034.",
    "The model equation is Methods 'Calculation of Km value' Equation 2 and",
    "is reprinted inside the Figure 5 panel; the km estimate is in",
    "Results 'Directional Permeability of Sorafenib across MDCKII Cells'",
    "and in the Figure 5 panel.",
    sep = " "
  )

  vignette <- "Agarwal_2011_sorafenib_transporters"

  units <- list(
    time          = "h",
    dosing        = "(none; static in vitro concentration-response model driven by an external sorafenib concentration covariate)",
    concentration = "(observation pappBa is the apparent B-to-A permeability of sorafenib in units of 1e-6 cm/s; driving covariate CP_SORAFENIB_NGML is the sorafenib concentration applied to the donor/basolateral compartment in ng/mL)"
  )

  covariateData <- list(
    CP_SORAFENIB_NGML = list(
      description        = "Sorafenib concentration applied to the donor (basolateral) compartment of the Transwell, supplied as a covariate. Reused canonical: in this in-vitro model the column carries a donor-compartment buffer concentration rather than a plasma concentration, but the quantity and units are identical.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Methods 'Calculation of Km value' states that B-to-A permeability was measured with a range of concentrations from 2 ng/mL to 30 ug/mL applied to the donor side; the Figure 5 x-axis spans 1 to 100000 ng/mL and the plotted observations run from about 2 to about 32000 ng/mL, consistent with the Methods range.",
        "CAUTION -- the Figure 5 CAPTION instead states the range as '0.001 nM to 80 nM' (about 0.0006 to 52 ng/mL). That contradicts both the Methods text and the figure's own x-axis, and is one of three demonstrable errors in that caption (see the model file's km comment). The Methods range is the one this file follows.",
        "The paper reports km in both mass and molar units (3.6 ng/mL = 5.5 nM), which implies a molecular weight of about 655 g/mol. That is NOT the sorafenib free base (464.8 g/mol) but is close to sorafenib tosylate (637.0 g/mol), the salt form the authors purchased (Methods 'Chemicals and Reagents'). The companion P-gp model's dual-unit report (15.9 ug/mL = 25 uM) implies 636 g/mol, corroborating that both conversions in this paper use the tosylate molecular weight. Users converting between mass and molar units should apply the same convention the authors did, or the potencies will be inconsistent with the printed values."
      ),
      source_name        = "sorafenib concentration in the donor compartment (C)"
    )
  )

  population <- list(
    species          = "in vitro (MDCKII canine kidney epithelial cells stably transfected with murine Bcrp1)",
    n_subjects       = NA_integer_,
    n_studies        = 1L,
    age_range        = NA_character_,
    weight_range     = NA_character_,
    sex_female_pct   = NA_real_,
    race_ethnicity   = NA_character_,
    disease_state    = "Not applicable -- transfected cell line.",
    dose_range       = "Donor-side sorafenib concentrations from 2 ng/mL to 30 ug/mL (Methods 'Calculation of Km value').",
    regions          = NA_character_,
    notes            = paste(
      "Polarised MDCKII monolayers stably expressing murine Bcrp1, grown on six-well Transwells at 2e5 cells/well until confluent, with the receiver compartment sampled at 0, 10, 20, 30, 45, 60 and 90 min (Methods 'Directional Flux Studies in MDCKII Cells'). n = 3 per concentration for the Figure 5 affinity experiment.",
      "Apparent permeability was computed as Papp = (dQ/dt) / (A * C0) with monolayer area A = 4.67 cm^2 (Methods Equation 1). That equation is a data-reduction step applied BEFORE the model is fitted, not part of the model, so it is not encoded here; the model takes the resulting permeability as its observation.",
      "Assay validity: the prototypical BCRP substrate prazosin accumulated to about 10% of wild-type levels in the Bcrp1 transfects, and the BCRP inhibitor Ko143 (200 nM) abolished the directionality of sorafenib transport entirely, confirming the efflux signal is BCRP-mediated (Results, Figures 2 and 3).",
      "Scale check against an independent experiment: Figure 3 reports a B-to-A permeability of about 4.3e-6 cm/s in the Bcrp1 transfects at the 2 ng/mL tracer concentration versus about 0.29e-6 cm/s A-to-B, an efflux ratio of about 15 that matches the '15-fold greater' statement in the Results text. This model predicts 3.3e-6 cm/s at 2 ng/mL, the same order and consistent with between-experiment variability (Figure 3 used n = 6, Figure 5 n = 3).",
      "The paper's in vivo work (single 10 mg/kg IV dose and 48-h intraperitoneal infusions in FVB wild-type, Mdr1a/b(-/-), Bcrp1(-/-) and Mdr1a/b(-/-)Bcrp1(-/-) mice) was analysed by NONCOMPARTMENTAL analysis only -- no compartmental or physiological model was fitted to the mouse data, so it is not extractable and is not represented by a model file. Its quantitative results are reproduced as context in the validation vignette.",
      "Sorafenib is NOT a P-gp substrate in this system (no directional transport in MDR1 transfects, Figure 4), and it does not INHIBIT BCRP-mediated transport of prazosin or mitoxantrone (Figure 6B) despite being a high-affinity BCRP substrate -- the authors infer binding at a site distinct from those probes."
    )
  )

  ini({
    # ==================================================================
    # Inhibitory Emax model for saturable BCRP-mediated B-to-A
    # permeability. Agarwal 2011 Methods 'Calculation of Km value',
    # Equation 2, reprinted inside the Figure 5 panel:
    #
    #   E = Emax * (1 - C / (Km_app + C))
    #
    # 'where E is the B-to-A permeability of sorafenib, Emax is the
    # maximum permeability seen at the lowest sorafenib concentration,
    # Km_app is the affinity constant determined by the concentration of
    # sorafenib at which half-maximal inhibition is seen, and C is the
    # concentration of sorafenib in the donor compartment.'
    # ==================================================================

    lkm <- log(3.6)
    label("Log apparent affinity constant of BCRP for sorafenib (ng/mL)")
    # PRINTED. Agarwal 2011 Results: 'An inhibitory Emax model fitted to
    # the data estimated a Km_app of 5.5 +/- 1.2 nM'. The Figure 5 panel
    # prints the same estimate in both units: 'Km = 3.6 +/- 0.78 ng/ml
    # = 5.5 +/- 1.2 nM'. The +/- is the standard error of the estimate
    # (Figure 5 caption), not a variance component. The ng/mL form is
    # used here so the potency is on the same scale as the covariate
    # column, and because the molar conversion the authors used relies
    # on the tosylate molecular weight (see covariateData notes).
    #
    # INDEPENDENT CHECK: refitting Equation 2 to the sixteen observed
    # points digitised from Figure 5, with BOTH parameters free, returns
    # km = 3.13 (SE 0.41) ng/mL -- within one standard error of the
    # printed 3.6 +/- 0.78. The printed value is reproducible from the
    # published data and is the one encoded here.
    #
    # CAUTION on the Figure 5 caption, which is unreliable: it states
    # the affinity as '5 nM' (Results and the figure panel both say
    # 5.5), calls the model 'sigmoid' (Methods Equation 2 and the
    # equation printed in the panel are both the NON-sigmoid hyperbolic
    # form, with no Hill exponent), and gives the concentration range as
    # '0.001 nM to 80 nM' (Methods and the panel's own x-axis both say
    # 2 ng/mL to 30 ug/mL). This file follows Methods Equation 2 and the
    # Results/panel estimate; the caption discrepancies are recorded in
    # the vignette Errata.

    lemax <- log(5.07)
    label("Log maximum B-to-A permeability at tracer concentration (1e-6 cm/s)")
    # NOT REPORTED in the source in any table, text passage or figure
    # panel. DIGITISED from Figure 5. The published panel was rendered
    # at 400 dpi, the plot frame located from its axis rules (C = 1 and
    # C = 1e5 on the log x-axis; 0 and 5 x 1e-6 cm/s on the y-axis), and
    # the filled 'Observed' markers recovered by morphological erosion,
    # giving sixteen points spanning C = 1.9 to 3.2e4 ng/mL. Equation 2
    # was then refitted to them with km held at the PRINTED 3.6 ng/mL,
    # giving emax = 5.07 (SE 0.11) x 1e-6 cm/s, R^2 = 0.987. Per the
    # standing policy, figure-fitting fills gaps but never overrides a
    # printed value, so km stays at 3.6 rather than at the 3.13 the free
    # fit prefers.
    #
    # Bracketing / uncertainty: fitting Equation 2 to the plotted
    # 'Model Predicted' dashed curve instead of to the observed points
    # gives emax = 5.34 with km held at 3.6. The two routes bracket
    # emax at roughly 5.1-5.3 x 1e-6 cm/s. The observed-point route is
    # encoded because it fits the authors' own equation to the authors'
    # own data at the authors' own printed km; the dashed curve is NOT
    # internally consistent with Equation 2 at km = 3.6 (fitting the
    # curve with km free returns km = 2.27, well below the printed
    # value, and the curve's own half-maximal point implies emax = 4.94
    # while its C = 1 value implies emax = 6.15), which is a genuine
    # inconsistency in the published figure and is documented in the
    # vignette Errata.

    # ==================================================================
    # No inter-individual variability and no residual error. This is a
    # single fitted concentration-response curve from one cell-culture
    # experiment; the source reports a point estimate with its standard
    # error and replicate SDs on the plotted means, but no variance
    # components. Per the standing policy on unreported residual error,
    # addSd is fixed at 0 so the model returns the deterministic curve.
    # ==================================================================
    addSd <- fixed(0)
    label("Additive residual error on apparent permeability (ZERO - not reported in source)")
  })

  model({
    km <- exp(lkm)
    emax <- exp(lemax)

    # ==================================================================
    # Apparent B-to-A permeability of sorafenib (1e-6 cm/s) at the
    # donor-compartment concentration CP_SORAFENIB_NGML, written in the
    # exact form printed in Methods Equation 2 and the Figure 5 panel:
    #
    #   pappBa = emax * (1 - C / (km + C))
    #
    # This is algebraically emax * km / (km + C), so at C = km the
    # permeability is emax / 2 -- the 'concentration at which
    # half-maximal inhibition is seen' of the paper's own definition.
    # That identity is asserted in the validation vignette.
    # ==================================================================
    pappBa <- emax * (1 - CP_SORAFENIB_NGML / (km + CP_SORAFENIB_NGML))

    pappBa ~ add(addSd)
  })
}
