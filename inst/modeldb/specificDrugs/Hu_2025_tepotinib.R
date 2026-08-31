Hu_2025_tepotinib <- function() {
  description <- paste(
    "Clinical PK submodule of a quantitative systems pharmacology (QSP) model.",
    "Two-compartment, first-order-absorption model for oral tepotinib in humans;",
    "the plasma-PK module that supplies drug exposure to the Hu 2025",
    "MET-aberrant NSCLC QSP framework. Only this PK module is packaged: the",
    "69-reaction MET / EGFR / ALK / ROS1 signalling and tumour-growth network of",
    "the same paper is not integrable from the published Supplementary Table S1",
    "(seven required quantities are reported nowhere) and is deferred.",
    "CAUTION -- the packaged parameters are transcribed verbatim from",
    "Supplementary Table S1, and that printed set does NOT reproduce the paper's",
    "own human tepotinib simulation in Supplementary Figure S4F: it gives a",
    "393 h terminal half-life where the figure shows about 29 h, and is roughly",
    "60-fold low on average over 13-199 h. The contradiction is internal to the",
    "publication and is documented in full in the validation vignette; no",
    "extractor-fitted value has been substituted. Deterministic: the source",
    "reports no between-subject or residual variability for any PK module."
  )
  reference <- paste(
    "Hu J, Zhao Y, Rao Q, Li G, Jiao Z, Wang H, Qu Y, Xu S, Gu Z, Wang T,",
    "Chen Z, Zhao C, Zhou G. Combining mechanistic quantitative systems",
    "pharmacology modeling and patient-derived organoid testing in MET-aberrant",
    "non-small cell lung cancer for high-throughput combination efficacy",
    "analysis and personalized treatment design.",
    "Front Pharmacol. 2025;16:1685468. doi:10.3389/fphar.2025.1685468.",
    "PMCID PMC12682884. The model definition is in Supplementary Table S1",
    "(Table1.xlsx; Species, Reactions and Parameters sheets) and Supplementary",
    "Figures S4/S9 (DataSheet1.pdf).",
    "The tepotinib PK parameters of Supplementary Table S1 are themselves cited",
    "to Johne A, Scheible H, Becker A, van Lier JJ, Wolna P, Meyring M.",
    "Invest New Drugs. 2020;38(5):1507-1519. doi:10.1007/s10637-020-00926-1",
    "(PMID 32221754), which is not on disk for this extraction."
  )
  vignette <- "Hu_2025_tepotinib"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Supplementary Table S1, Species sheet:
  # tepo_dose and tepo_c are assigned to the "central" compartment column and
  # tepo_p to "peripheral"; reactions v67 / v68 / v69 are the absorption,
  # elimination and inter-compartmental distribution fluxes respectively.
  compartmentData <- list(
    depot       = list(analyte = "tepotinib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "tepotinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tepotinib", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species       = "human",
    n_subjects    = NA_integer_,
    n_studies     = NA_integer_,
    disease_state = paste(
      "Advanced non-small cell lung cancer with MET exon 14 skipping mutations",
      "or MET amplification. The QSP framework applies this PK module to two",
      "virtual-patient cohorts of 5000 subjects each (Results, 'QSP model-based",
      "virtual clinical trials'); the PK module itself carries no",
      "between-subject variability, so it is a typical-value model."
    ),
    dose_range    = paste(
      "Tepotinib 500 mg once daily (containing 450 mg active moiety) is the",
      "clinically approved regimen simulated in Figure 7B; 600 mg QD was also",
      "simulated and gave no further gain in objective response rate.",
      "Supplementary Figure S4F, the human PK panel for tepotinib, is a 500 mg",
      "single oral dose."
    ),
    notes = paste(
      "PRINTED PARAMETERS DO NOT REPRODUCE THE PAPER'S OWN FIGURE. Every value",
      "in ini() is transcribed verbatim from Supplementary Table S1, Parameters",
      "sheet, where all five tepotinib PK entries are attributed to Johne 2020",
      "(PMID 32221754): ka_tepo = 0.2 1/h, kcl_tepo = 1.32 1/h, Q = 1.32 L/h,",
      "Vc = 34 L and vp_tepo = 726.8 L. Simulating that set at the 500 mg single",
      "oral dose of Supplementary Figure S4F gives a terminal half-life of",
      "393 h, whereas the simulation trace drawn in that figure declines",
      "log-linearly with a half-life of about 29 h; the geometric-mean ratio of",
      "predicted to figure concentration is 0.017 over 13-199 h (about 60-fold",
      "low), with a worst-case discrepancy near 460-fold. The mismatch is",
      "structural, not a scale factor: kcl_tepo = 1.32 1/h empties the central",
      "compartment with a 0.5 h half-life while Q/vp_tepo = 0.0018 1/h returns",
      "drug from the periphery far too slowly. No permutation of the five",
      "printed numbers reproduces the figure. Per operator ruling for this",
      "extraction, source fidelity was preferred over an extractor-side fit: the",
      "printed values ship as printed and the contradiction is documented rather",
      "than repaired. Users who need a profile matching Supplementary Figure S4F",
      "must supply their own disposition parameters; see the validation vignette",
      "for the quantified comparison.",
      "ADDITIONAL TRANSCRIPTION DEFECT. Reaction v69 as printed,",
      "'Q/Vc*tepo_c*Vc/Vp_tepo-Q/Vp_tepo*tepo_p*Vp_tepo/Vc', reduces to",
      "Q*tepo_c/Vp_tepo - Q*tepo_p/Vc: the two volume divisors are transposed",
      "relative to the standard inter-compartmental flux, and the expression is",
      "not mass-conserving. The mass-conserving form Q*(tepo_c/Vc -",
      "tepo_p/Vp_tepo) is encoded instead, since a distribution flux that",
      "creates or destroys drug cannot be what the authors integrated. This",
      "correction is not what causes the disagreement with Figure S4F: the",
      "literal printed flux fits that figure far worse still.",
      "STRUCTURE. Supplementary Table S1's footnote names tepotinib as the",
      "'two-compartment' clinical PK exemplar, reaction v69 is an",
      "inter-compartmental distribution flux, and the Species sheet carries a",
      "tepo_p peripheral species; the two-compartment structure encoded here",
      "therefore uses all five printed parameters. Supplementary Figure S9, the",
      "topology schematic, draws tepotinib with a central compartment only,",
      "which contradicts the reaction list; per the standing convention that a",
      "printed equation outranks a schematic, the reaction list is followed."
    )
  )

  ini({
    # ---- Absorption ----
    # All five tepotinib PK values are literature values re-used by the QSP
    # paper from Johne 2020 without re-fitting in this publication, so each is
    # encoded fixed().
    lka <- fixed(log(0.2)); label("First-order absorption rate constant ka (1/h)")   # Supplementary Table S1, Parameters sheet: ka_tepo = 0.2 (hour^-1), ref PMID 32221754; used by reaction v67 = ka_tepo*tepo_dose

    # ---- Disposition ----
    # Reaction v68 is kcl_tepo*tepo_c, i.e. a first-order elimination RATE
    # CONSTANT acting on the central amount, so the corresponding clearance is
    # CL = kcl_tepo * Vc = 1.32 1/h * 34 L = 44.88 L/h. The arithmetic is
    # written out here so the printed micro-constant stays in the source trace;
    # CL / Vc reproduces kcl_tepo exactly. Storing the flow (CL, Q) and volume
    # (Vc, Vp) rather than the micro-constants also prevents rxSolve()'s
    # default useLinCmt rewrite from silently collapsing the peripheral
    # compartment.
    lcl <- fixed(log(1.32 * 34)); label("Apparent oral clearance CL/F (L/h)")        # Supplementary Table S1: kcl_tepo = 1.32 (hour^-1) x Vc = 34 (L), both ref PMID 32221754
    lvc <- fixed(log(34));        label("Apparent central volume Vc/F (L)")          # Supplementary Table S1, Parameters sheet: Vc = 34 (L), ref PMID 32221754
    lq  <- fixed(log(1.32));      label("Apparent inter-compartmental clearance Q/F (L/h)")  # Supplementary Table S1, Parameters sheet: Q = 1.32 (L*hour^-1), ref PMID 32221754
    lvp <- fixed(log(726.8));     label("Apparent peripheral volume Vp/F (L)")       # Supplementary Table S1, Parameters sheet: vp_tepo = 726.8 (L), ref PMID 32221754

    # ---- Bioavailability anchor ----
    # Reaction v67 is ka_tepo*tepo_dose with no bioavailability term, so the QSP
    # module absorbs the full administered amount and the disposition
    # parameters above are apparent (F = 1) parameters. Johne 2020, the cited PK
    # source, is not on disk for this extraction, so no absolute bioavailability
    # is applied here.
    lfdepot <- fixed(log(1)); label("Bioavailability of the depot compartment (fraction)")  # Supplementary Table S1 reaction v67 carries no F term

    # ---- Residual error ----
    # The QSP PK modules are deterministic simulations; no residual-error and no
    # between-subject variance is reported anywhere in the paper or its
    # supplement. Encoded fixed(0) per the missing-value convention. Users
    # running stochastic simulations must supply their own magnitudes.
    propSd <- fixed(0); label("Proportional residual error on plasma tepotinib concentration (fraction); not published")  # no residual-error model reported in the paper or Supplementary Table S1
  })

  model({
    # 1. Individual parameters
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # 2. Two-compartment oral disposition, Supplementary Table S1 reactions
    #    v67 (absorption, ka_tepo*tepo_dose), v68 (elimination, kcl_tepo*tepo_c)
    #    and v69 (distribution). v69 is encoded in its mass-conserving form
    #    Q*(central/vc - peripheral1/vp); as printed its two volume divisors are
    #    transposed and the flux does not conserve drug. See population$notes.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    # 3. Bioavailability
    f(depot) <- exp(lfdepot)

    # 4. Observation. central is in mg and vc in L, so Cc is mg/L = ug/mL.
    #    Supplementary Figure S4F plots molar concentration; divide Cc by
    #    1000 * 492.58 (tepotinib free-base molecular weight, g/mol) to convert
    #    ug/mL to mol/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
attr(Hu_2025_tepotinib, "message") <-
  "Parameters are transcribed verbatim from Supplementary Table S1, but that printed set does NOT reproduce the paper's own human tepotinib simulation in Supplementary Figure S4F: terminal half-life 393 h vs about 29 h in the figure, and roughly 60-fold low over 13-199 h. The contradiction is internal to the publication and is documented, not repaired. See the validation vignette before using this model to predict tepotinib exposure."
Hu_2025_tepotinib
