Agarwal_2011_sorafenib_pgp <- function() {
  description <- paste(
    "In vitro (MDCKII-MDR1). Sigmoidal Emax model for the inhibition of",
    "P-glycoprotein (P-gp, ABCB1) by sorafenib, measured as the fold",
    "increase in intracellular accumulation of the prototypical P-gp",
    "substrate vinblastine in polarised MDCKII monolayers stably",
    "transfected with human MDR1. The model is:",
    "vblAccum = e0 + (emax - e0) * C^hill / (ic50^hill + C^hill),",
    "with ic50 = 15900 ng/mL (15.9 ug/mL, 25 uM). Sorafenib is NOT itself",
    "a P-gp substrate in this system (Figure 4 shows no directional",
    "transport), so this curve describes sorafenib acting purely as a P-gp",
    "INHIBITOR -- a perpetrator, not a victim. The authors flag it as",
    "clinically relevant because sorafenib plasma concentrations on the",
    "accepted 100-400 mg b.i.d. regimen reach 1-15 uM, within range of the",
    "25 uM IC50, so sorafenib may alter the tissue distribution of",
    "co-administered P-gp substrates. Only ic50 is reported in the source;",
    "e0, emax and hill were recovered by digitising Figure 6A (see the",
    "parameter comments and the vignette Errata). The output is a",
    "dimensionless fold increase relative to untreated control, so the",
    "model carries no inter-individual variability and no residual error",
    "-- it is a single fitted concentration-response curve, not a",
    "population model. Companion model file Agarwal_2011_sorafenib_bcrp.R",
    "carries the BCRP affinity curve from the same paper.",
    sep = " "
  )

  reference <- paste(
    "Agarwal S, Sane R, Ohlfest JR, Elmquist WF.",
    "The role of the breast cancer resistance protein (ABCG2) in the",
    "distribution of sorafenib to the brain.",
    "J Pharmacol Exp Ther. 2011;336(1):223-233.",
    "doi:10.1124/jpet.110.175034.",
    "The model equation is Methods 'P-gp and BCRP Inhibition Assays'",
    "Equation 3 and is reprinted inside the Figure 6A panel; the ic50",
    "estimate is in Results 'Sorafenib as a P-gp or BCRP Inhibitor' and in",
    "the Figure 6A panel.",
    sep = " "
  )

  vignette <- "Agarwal_2011_sorafenib_transporters"

  units <- list(
    time          = "h",
    dosing        = "(none; static in vitro concentration-response model driven by an external sorafenib concentration covariate)",
    concentration = "(observation vblAccum is the fold increase in intracellular vinblastine accumulation relative to untreated control, dimensionless; driving covariate CP_SORAFENIB_NGML is the sorafenib concentration in the incubation medium in ng/mL)"
  )

  covariateData <- list(
    CP_SORAFENIB_NGML = list(
      description        = "Sorafenib concentration in the cell-assay incubation medium, supplied as a covariate. Reused canonical: in this in-vitro model the column carries an incubation-medium concentration rather than a plasma concentration, but the quantity and units are identical. The same column drives the companion BCRP model, where it carries a Transwell donor-compartment concentration.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Methods 'P-gp and BCRP Inhibition Assays' states that intracellular accumulation was determined in the presence of sorafenib concentrations spanning 20 ng/mL to 60 ug/mL. The Figure 6A x-axis is labelled in ug/mL and spans about 0.4 to 120 ug/mL, with plotted observations from about 0.64 to 64 ug/mL; the lower decades of the Methods range appear only in the Figure 6B (BCRP probe) panel, whose axis starts at 0.01 ug/mL.",
        "UNIT CONVERSION: the source prints this model's potency in ug/mL (15.9) while the covariate column is in ng/mL, so the ini() value is the exact, lossless conversion 15.9 ug/mL = 15900 ng/mL. A single ng/mL column is used for both Agarwal 2011 models so that one canonical covariate serves the whole paper.",
        "The paper reports ic50 in both mass and molar units (15.9 ug/mL = 25 uM), implying a molecular weight of 636 g/mol. That is NOT the sorafenib free base (464.8 g/mol) but matches sorafenib tosylate (637.0 g/mol), the salt the authors purchased (Methods 'Chemicals and Reagents'). The companion BCRP model's dual-unit report (3.6 ng/mL = 5.5 nM) implies about 655 g/mol, corroborating the same tosylate convention. Users converting between mass and molar units should apply the authors' convention or the potencies will disagree with the printed values.",
        "Clinical context supplied by the authors (Discussion): sorafenib plasma concentrations in humans on the accepted 100-400 mg b.i.d. regimen have been reported at 1-15 uM, i.e. the same order as this 25 uM IC50."
      ),
      source_name        = "sorafenib concentration (C)"
    )
  )

  population <- list(
    species          = "in vitro (MDCKII canine kidney epithelial cells stably transfected with human MDR1)",
    n_subjects       = NA_integer_,
    n_studies        = 1L,
    age_range        = NA_character_,
    weight_range     = NA_character_,
    sex_female_pct   = NA_real_,
    race_ethnicity   = NA_character_,
    disease_state    = "Not applicable -- transfected cell line.",
    dose_range       = "Sorafenib concentrations spanning 20 ng/mL to 60 ug/mL (Methods 'P-gp and BCRP Inhibition Assays'); the Figure 6A observations run from about 0.64 to 64 ug/mL.",
    regions          = NA_character_,
    notes            = paste(
      "Confluent MDCKII-MDR1 monolayers in 24-well plates at 2e5 cells/well; intracellular accumulation of tritiated vinblastine measured at 60 min in the presence of increasing sorafenib concentrations, normalised to the protein concentration of the solubilised cell fraction and expressed relative to the no-sorafenib control (Methods 'Intracellular Accumulation' and 'P-gp and BCRP Inhibition Assays'). n = 4 per concentration.",
      "Assay validity: accumulation of the prototypical P-gp substrate vinblastine in the MDR1 transfects was about 10% of wild-type levels, confirming functional P-gp expression (Results, Figure 1B).",
      "IMPORTANT ASYMMETRY. Sorafenib inhibits P-gp (this model) but is not transported by it -- accumulation of sorafenib itself in MDR1 transfects did not differ from wild-type (Figure 1B) and there was no directional flux (Figure 4). Conversely sorafenib IS a high-affinity BCRP substrate (companion model, km = 3.6 ng/mL) yet did NOT inhibit BCRP-mediated transport of either prazosin or mitoxantrone (Figure 6B), which is why no BCRP inhibition curve exists to extract. The authors infer that sorafenib binds BCRP at a site not overlapping those probes.",
      "The paper's in vivo work (single 10 mg/kg IV dose and 48-h intraperitoneal infusions in FVB wild-type, Mdr1a/b(-/-), Bcrp1(-/-) and Mdr1a/b(-/-)Bcrp1(-/-) mice, plus an elacridar interaction study) was analysed by NONCOMPARTMENTAL analysis only -- no compartmental or physiological model was fitted, so it is not extractable and is not represented by a model file. Its quantitative results are reproduced as context in the validation vignette.",
      "Note the in vitro / in vivo discordance the authors highlight: P-gp shows no effect on sorafenib in vitro, yet brain-to-plasma ratios rose about 10-fold in the double Mdr1a/b(-/-)Bcrp1(-/-) knockout versus about 4-fold in the Bcrp1(-/-) single knockout, so P-gp does restrict brain entry in vivo once BCRP is absent."
    )
  )

  ini({
    # ==================================================================
    # Sigmoid Emax model for P-gp inhibition by sorafenib. Agarwal 2011
    # Methods 'P-gp and BCRP Inhibition Assays', Equation 3, reprinted
    # inside the Figure 6A panel:
    #
    #   E = E0 + (Emax - E0) * C^gamma / (IC50^gamma + C^gamma)
    #
    # 'where E is the fold increase in the accumulation of substrate
    # seen in the presence of sorafenib relative to control (no
    # sorafenib), E0 is the accumulation of probe in the absence of
    # sorafenib normalized to unity, Emax is the maximum increase in
    # accumulation of probe, IC50 is the concentration of sorafenib at
    # which half-maximal effect is seen, C is the concentration of
    # sorafenib, and gamma is the shape factor that determines the slope
    # of the curve.'
    #
    # DIGITISATION NOTE applying to e0, emax and hill below. Only ic50
    # is reported numerically. The other three were recovered by pixel
    # analysis of the published Figure 6A 'Model Predicted' dashed
    # curve -- i.e. from the authors' own fitted curve, not from a
    # refit of the scatter. The panel was rendered at 400 dpi, the axes
    # calibrated from the '1', '10' and '100' x tick labels (one decade
    # = 478.1 px) and the 0-14 y frame, and the dashed curve traced in
    # the columns at least 30 px clear of every marker so that error
    # bars could not contaminate it. That trace is extremely well
    # determined (405 curve columns, R^2 = 0.99999 for the
    # four-parameter form) and SELF-VALIDATES: refitting with all four
    # parameters free returns ic50 = 15.48 ug/mL against the printed
    # 15.9 +/- 3.8, a 2.7% recovery error. The values below come from
    # refitting the traced curve with ic50 held at the PRINTED value,
    # per the standing policy that figure-fitting fills gaps but never
    # overrides a printed number. Quoted SEs are the curve-fit SEs and
    # describe how precisely the published curve was traced; they are
    # NOT estimation uncertainty from the original experiment.
    # ==================================================================

    lic50 <- log(15900)
    label("Log sorafenib concentration giving half-maximal P-gp inhibition (ng/mL)")
    # PRINTED. Agarwal 2011 Results: 'the IC50 of sorafenib for
    # inhibition of P-gp was estimated to be 25 +/- 6 uM'. The Figure 6A
    # panel prints the same estimate in both units: 'IC50 = 15.9 +/- 3.8
    # ug/ml = 25 +/- 6 uM'. The +/- is the standard error of the
    # estimate (Figure 6 caption), not a variance component.
    # Encoded as the exact conversion 15.9 ug/mL = 15900 ng/mL so that
    # the potency is on the same ng/mL scale as the covariate column
    # shared with the companion BCRP model.

    le0 <- log(1.130)
    label("Log baseline vinblastine accumulation with no sorafenib (fold of control)")
    # NOT REPORTED numerically. DIGITISED from the Figure 6A fitted
    # curve: e0 = 1.130 (curve-fit SE 0.002).
    #
    # DELIBERATE DEVIATION FROM A LITERAL READING OF THE TEXT, flagged
    # in the vignette Errata. Methods Equation 3 describes E0 as 'the
    # accumulation of probe in the absence of sorafenib normalized to
    # unity', which reads as E0 = 1 exactly, and the lowest observed
    # point in Figure 6A does sit at 1.00. But the authors' own plotted
    # curve has a left asymptote of 1.130, not 1.000, so E0 was
    # evidently ESTIMATED rather than held at 1: the quoted sentence
    # describes how the DATA were normalised, not a constraint on the
    # fitted parameter. Encoding 1.130 reproduces the published figure;
    # encoding 1.000 does not. Users who want the literal
    # normalised-to-unity form can set le0 <- log(1) and refit, which
    # shifts emax to about 10.25 and hill to about 2.02 (that variant
    # fits the observed scatter with R^2 = 0.95).

    lemax <- log(10.231)
    label("Log maximum vinblastine accumulation at saturating sorafenib (fold of control)")
    # NOT REPORTED numerically. DIGITISED from the Figure 6A fitted
    # curve with ic50 held at the printed 15.9 ug/mL: emax = 10.231
    # (curve-fit SE 0.004). Consistent with the visual plateau of the
    # dashed curve at about 10 fold and with the highest observed points
    # (9.9 fold at about 32 ug/mL, 9.5 fold at about 64 ug/mL, both
    # recovered by the same digitisation).

    lhill <- log(1.842)
    label("Log Hill shape factor of the P-gp inhibition curve (dimensionless)")
    # NOT REPORTED numerically. DIGITISED from the Figure 6A fitted
    # curve with ic50 held at the printed 15.9 ug/mL: hill = 1.842
    # (curve-fit SE 0.004). The source calls this parameter 'gamma',
    # 'the shape factor that determines the slope of the curve'; the
    # canonical ini() name is `lhill` with bare in-model name `hill`
    # (parameter-names.md, 'Sigmoidal PD shape parameters', which
    # directs that lhill be used regardless of whether the paper writes
    # gamma, hill or n).
    #
    # Stability check, since a Hill exponent is the parameter most
    # vulnerable to digitisation error. Three independent routes:
    #   (a) traced curve, ic50 ALSO free      -> hill = 1.905 (SE 0.002)
    #   (b) traced curve, ic50 fixed at 15.9  -> hill = 1.842 (SE 0.004)
    #   (c) the eleven digitised OBSERVED points, ic50 fixed at 15.9 and
    #       e0 fixed at 1                     -> hill = 2.02  (SE 0.46)
    # All three agree on hill near 1.8-2.0. Freeing e0 as well on the
    # eleven observed points alone gives hill = 4.0 with SE 1.4 (95% CI
    # roughly 1.2 to 6.8) -- that route simply cannot resolve the Hill
    # exponent from eleven noisy means, and its interval covers the
    # curve-derived value, so it neither supports nor contradicts it.
    # What every route does establish is that hill is comfortably above
    # 1, so the curve is genuinely sigmoidal rather than hyperbolic.

    # ==================================================================
    # No inter-individual variability and no residual error. This is a
    # single fitted concentration-response curve from one cell-culture
    # experiment; the source reports a point estimate with its standard
    # error and replicate SDs on the plotted means, but no variance
    # components. Per the standing policy on unreported residual error,
    # addSd is fixed at 0 so the model returns the deterministic curve.
    # ==================================================================
    addSd <- fixed(0)
    label("Additive residual error on fold increase in accumulation (ZERO - not reported in source)")
  })

  model({
    ic50 <- exp(lic50)
    e0 <- exp(le0)
    emax <- exp(lemax)
    hill <- exp(lhill)

    # ==================================================================
    # Fold increase in intracellular vinblastine accumulation relative
    # to untreated control, at sorafenib concentration
    # CP_SORAFENIB_NGML, in the exact form printed in Methods
    # Equation 3 and the Figure 6A panel:
    #
    #   vblAccum = e0 + (emax - e0) * C^hill / (ic50^hill + C^hill)
    #
    # At C = ic50 this returns (e0 + emax) / 2, the half-maximal effect
    # that defines the reported IC50. That identity is asserted in the
    # validation vignette.
    # ==================================================================
    vblAccum <-
      e0 + (emax - e0) * CP_SORAFENIB_NGML^hill /
        (ic50^hill + CP_SORAFENIB_NGML^hill)

    vblAccum ~ add(addSd)
  })
}
