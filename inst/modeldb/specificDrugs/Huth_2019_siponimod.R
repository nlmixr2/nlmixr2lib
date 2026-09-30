Huth_2019_siponimod <- function() {
  description <- paste(
    "Two-compartment oral population pharmacokinetic reduction of the",
    "Simcyp full-PBPK model for the sphingosine-1-phosphate receptor",
    "modulator siponimod, carrying the CYP2C9 diplotype and the",
    "CYP3A4 / CYP2C9 perpetrator drug-drug-interaction layer (Huth 2019).",
    "The source model was built in the Simcyp Population-based Simulator",
    "(V16) with a full perfusion-limited whole-body distribution model,",
    "so the per-organ partition coefficients, volumes and blood flows are",
    "platform database outputs and the PBPK model itself cannot be",
    "encoded here. What IS fully reported is the siponimod compound layer",
    "and, critically, both of its load-bearing disposition inputs are",
    "clinical observations rather than platform predictions: Table 5",
    "footnote f states the Vss of 1.45 L/kg is the observed mean Vss",
    "converted on the 80.08 kg mean study body weight, and footnote g",
    "states the retrograde clearance calculation used the observed",
    "geometric mean CL of 3.12 L/h.",
    "Two parameters could not be transcribed and were back-solved from",
    "the paper's own printed single-dose summary statistics, because a",
    "whole-body model publishes no central volume or distribution",
    "clearance: the central volume and the intercompartmental clearance",
    "were set to reproduce the Table 1 oral Tmax of 3.45 h and Cmax of",
    "7.76 ng/mL per mg. The peripheral volume is then fixed by the",
    "printed Vss. Those two numbers are corroborated out of sample by the",
    "Table 1 intravenous arm: a 2.92 h infusion of 0.25 mg gives a",
    "predicted Cmax of 2.99 ng/mL against the paper's 3.02 (-1.1 pct).",
    "The CYP2C9 diplotype enters as a relative CYP2C9 activity on the",
    "CYP2C9 clearance arm, derived from the six fraction-metabolised",
    "pairs in Table 5. That derivation is self-checking: the implied",
    "residual fraction closes each genotype row of Table 5 to 1.000, and",
    "the six genotype terminal half-lives reproduce Table 2 within 6 pct",
    "across a fivefold exposure range without any further adjustment.",
    "Coadministered CYP3A4 / CYP2C9 modulators act through relative",
    "enzyme activities that scale the matching clearance arm, gut-wall",
    "first-pass extraction and hepatic first-pass extraction together.",
    "Each modulator's activity was back-solved from its six published AUC",
    "ratios alone, with CYP2C9 activity held at 1 for the four",
    "perpetrators the paper classifies as CYP3A4-selective; the published",
    "Cmax ratios were held out entirely and are reproduced with a worst",
    "case of 7.2 pct and a median absolute error of 2.3 pct over all 42",
    "modulator-by-genotype cells.",
    "This is a typical-value simulation model: the source reports no",
    "inter-individual variance components and no residual-error model, so",
    "there are no etas and propSd is fixed at zero.",
    sep = " "
  )
  reference <- paste(
    "Huth F, Gardin A, Umehara K, He H. (2019). Prediction of the Impact",
    "of Cytochrome P450 2C9 Genotypes on the Drug-Drug Interaction",
    "Potential of Siponimod With Physiologically-Based Pharmacokinetic",
    "Modeling: A Comprehensive Approach for Drug Label Recommendations.",
    "Clin Pharmacol Ther 106(5):1113-1124. doi:10.1002/cpt.1547.",
    sep = " "
  )
  vignette <- "Huth_2019_siponimod"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The source is a whole-body PBPK, so these three
  # states are the compartmental reduction's own states rather than named
  # organs; `peripheral1` lumps every extravascular tissue the Simcyp
  # perfusion-limited model resolves separately.
  compartmentData <- list(
    depot = list(
      analyte = "siponimod",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "siponimod",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "siponimod",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    CYP2C9_S1_COUNT = list(
      description = "CYP2C9*1 (wild-type) allele count",
      units = "(count, 0/1/2 alleles per subject)",
      type = "continuous",
      reference_category = paste(
        "CYP2C9*1/*1 (CYP2C9_S1_COUNT == 2) is the reference diplotype;",
        "relative CYP2C9 activity 1 and total CL 3.12 L/h."
      ),
      notes = paste(
        "CYP2C9_S1_COUNT + CYP2C9_S2_COUNT + CYP2C9_S3_COUNT = 2 per",
        "subject. Huth 2019 simulates the six clinically relevant CYP2C9",
        "diplotypes in the white population; the Introduction gives their",
        "prevalences as 62-65 pct *1/*1, 20-24 pct *1/*2, 9-12 pct *1/*3,",
        "1-2 pct *2/*2, 1.4-1.7 pct *2/*3 and 0.3-0.4 pct *3/*3. The",
        "paper reports no CYP2C9 alleles beyond *1, *2 and *3, so the",
        "*5 / *6 / *8 count columns registered for the warfarin models",
        "are not carried here."
      ),
      source_name = "CYP2C9 genotype"
    ),
    CYP2C9_S2_COUNT = list(
      description = "CYP2C9*2 (rs1799853) reduced-function allele count",
      units = "(count, 0/1/2 alleles per subject)",
      type = "continuous",
      reference_category = "0 (no *2 allele); the *1/*1 diplotype is the model reference.",
      notes = paste(
        "Huth 2019 Introduction: subjects with the *1/*2 genotype behave",
        "as extensive metabolizers and *2/*2 as intermediate",
        "metabolizers. Enters only through the diplotype indicators."
      ),
      source_name = "CYP2C9 genotype"
    ),
    CYP2C9_S3_COUNT = list(
      description = "CYP2C9*3 (rs1057910) reduced-function allele count",
      units = "(count, 0/1/2 alleles per subject)",
      type = "continuous",
      reference_category = "0 (no *3 allele); the *1/*1 diplotype is the model reference.",
      notes = paste(
        "Huth 2019 Introduction: *1/*3 behaves as an intermediate",
        "metabolizer and *2/*3 and *3/*3 as poor metabolizers. The",
        "*3/*3 diplotype retains only 2.1 pct of the reference CYP2C9",
        "clearance arm, which is why CYP3A4 becomes the dominant",
        "elimination route in that subpopulation and why siponimod is",
        "contraindicated there."
      ),
      source_name = "CYP2C9 genotype"
    ),
    CONMED_ITRACONAZOLE = list(
      description = "Concomitant itraconazole (strong CYP3A4 inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP3A4 / CYP2C9 perpetrator coadministered",
      notes = paste(
        "Huth 2019 Table 3, steady-state arm. The paper classifies",
        "itraconazole as CYP3A4-selective, so its relative CYP2C9",
        "activity is held at 1 and only the CYP3A4 activity is",
        "back-solved. Huth 2019 also reports a single-dose itraconazole",
        "arm whose observed AUC ratios (0.90 for *1/*2 and 0.76 for",
        "*1/*3) were LOWER than 1 and were mispredicted by the source",
        "model itself; that anomaly is not encoded here."
      ),
      source_name = "itraconazole"
    ),
    CONMED_KETOCONAZOLE = list(
      description = "Concomitant ketoconazole (strong CYP3A4 inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP3A4 / CYP2C9 perpetrator coadministered",
      notes = paste(
        "Huth 2019 Table 3, steady-state arm, simulated only. The most",
        "potent CYP3A4 inhibitor in the paper (relative CYP3A4 activity",
        "0.084) and the arm that produces the largest predicted",
        "interaction of the whole analysis, an AUC ratio of 4.20 in the",
        "CYP2C9*3/*3 subpopulation."
      ),
      source_name = "ketoconazole"
    ),
    CONMED_ERYTHROMYCIN = list(
      description = "Concomitant erythromycin (moderate CYP3A4 inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP3A4 / CYP2C9 perpetrator coadministered",
      notes = "Huth 2019 Table 3, steady-state arm, simulated only; CYP3A4-selective.",
      source_name = "erythromycin"
    ),
    CONMED_FLUCONAZOLE = list(
      description = "Concomitant fluconazole (moderate CYP3A4 and CYP2C9 inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP3A4 / CYP2C9 perpetrator coadministered",
      notes = paste(
        "Huth 2019 Table 3. A dual perpetrator, so both relative",
        "activities are back-solved. This is the one inhibitor arm with",
        "clinical reference data: 200 mg twice daily on day 1 then",
        "200 mg once daily on days 2-19 with a single 4 mg siponimod",
        "dose on day 3 gave an observed AUC ratio of 1.98 against the",
        "source model's single-dose prediction of 2.15."
      ),
      source_name = "fluconazole"
    ),
    CONMED_FLUVOXAMINE = list(
      description = "Concomitant fluvoxamine (weak CYP3A4 and CYP2C9 inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP3A4 / CYP2C9 perpetrator coadministered",
      notes = paste(
        "Huth 2019 Table 3, steady-state arm, simulated only. A dual",
        "perpetrator, and the only one whose predicted AUC ratio",
        "DECREASES with worsening CYP2C9 function (1.42 in *1/*1 down to",
        "1.12 in *3/*3), because inhibiting CYP2C9 in a subpopulation",
        "that barely uses CYP2C9 buys nothing while the CYP3A4",
        "inhibition is only weak."
      ),
      source_name = "fluvoxamine"
    ),
    CONMED_RIFAMPICIN = list(
      description = "Concomitant rifampicin (strong CYP3A4 and moderate CYP2C9 inducer)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP3A4 / CYP2C9 perpetrator coadministered",
      notes = paste(
        "Huth 2019 Table 3. Dual inducer with clinical reference data:",
        "600 mg twice daily for 12 days against siponimod 2 mg once",
        "daily gave an observed AUC ratio of 0.43 and Cmax ratio of 0.55",
        "versus the source model's 0.32 and 0.50. Huth 2019 also reports",
        "that dropping the CYP2C9 induction parameters raises the",
        "predicted AUC ratio to 0.59, which is how the CYP2C9 arm of this",
        "perpetrator is identified."
      ),
      source_name = "rifampin"
    ),
    CONMED_EFV = list(
      description = "Concomitant efavirenz (moderate CYP3A4 inducer)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP3A4 / CYP2C9 perpetrator coadministered",
      notes = paste(
        "Huth 2019 Table 3, steady-state arm, simulated only. The",
        "Methods state the efavirenz model was built as a pure CYP3A4",
        "inducer because negligible CYP2C9 induction is expected from a",
        "moderate CYP3A4 inducer, so the relative CYP2C9 activity is held",
        "at 1. Here the reference category is simply no efavirenz rather",
        "than the alternative-antiretroviral-regimen reference used by",
        "the antiretroviral popPK models that share this column."
      ),
      source_name = "efavirenz"
    )
  )

  # Reported by the source model but not carried as covariates here.
  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. Huth 2019 Table 5 footnote f converts the observed",
        "Vss to 1.45 L/kg on the 80.08 kg mean body weight of the study",
        "subjects, so the Simcyp model does scale distribution volume",
        "with weight. This reduction fixes the 80.08 kg reference weight",
        "instead, because the clearance the paper anchors to is an",
        "absolute 3.12 L/h with no published weight relationship, and",
        "scaling volume alone would change the half-life with no basis",
        "in the source."
      ),
      units = "kg",
      type = "continuous",
      notes = "Simcyp healthy-volunteer population file; relationship not published."
    ),
    AGE = list(
      description = paste(
        "Age and female proportion. Huth 2019 Table S3 matched the",
        "simulation trial design to the actual clinical study's age",
        "ranges and female proportionality for the genotype",
        "verification, but no age or sex relationship is reported for any",
        "model parameter."
      ),
      units = "y",
      type = "continuous",
      notes = "Trial-design matching only; no published covariate relationship."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 15L,
    n_studies = 5L,
    age_range = paste(
      "Healthy adult volunteers; individual ages are not tabulated in",
      "the source. Simulation trial designs matched each clinical",
      "study's age range and female proportion (Table S3, not available when this model was built)."
    ),
    weight_median = "80.08 kg (mean body weight of the study subjects, Table 5 footnote f)",
    disease_state = paste(
      "Healthy volunteers throughout. The Methods state that population",
      "pharmacokinetic analysis found no clinically relevant PK",
      "difference between healthy volunteers and patients with multiple",
      "sclerosis, which is why the Simcyp healthy-volunteer population",
      "was used for every simulation."
    ),
    dose_range = paste(
      "Single intravenous 0.25 mg; single oral 0.1 to 75 mg; multiple",
      "oral 0.3 to 20 mg once daily for 28 days. The approved",
      "maintenance dose is 2 mg once daily, reduced to 1 mg once daily",
      "for CYP2C9*1/*3 and *2/*3 and contraindicated for *3/*3."
    ),
    regions = "Simcyp healthy-volunteer virtual population; genotype prevalences quoted for the white population.",
    genotypes = paste(
      "Six clinically relevant CYP2C9 diplotypes: *1/*1, *1/*2, *1/*3,",
      "*2/*2, *2/*3 and *3/*3. Observed genotype PK was available for",
      "*1/*1 (n = 12), *2/*3 (n = 6) and *3/*3 (n = 6); the other three",
      "diplotypes are model predictions only (Table 2)."
    ),
    studies = paste(
      "A single-ascending-dose study (0.1 to 75 mg, n = 6 to 8 per dose)",
      "and a multiple-ascending-dose study (0.3 to 20 mg once daily for",
      "28 days, n = 6 to 9 per dose) supplied the oral dose ladder of",
      "Table 1; an absolute-bioavailability study (n = 15) supplied the",
      "paired 0.25 mg intravenous and oral arms; a CYP2C9",
      "pharmacogenetic study supplied the genotype PK of Table 2 and the",
      "fluconazole interaction; and a dedicated rifampin interaction",
      "study supplied the induction reference data."
    ),
    notes = paste(
      "n_subjects records the 15 participants of the",
      "absolute-bioavailability study, which supplied both disposition",
      "anchors (the observed Vss and the observed clearance). This is a",
      "PBPK analysis rather than a population-PK fit, so there is no",
      "pooled analysis dataset and no estimated variance components. The",
      "simulated coefficients of variation in Table 1 (about 48 pct on",
      "AUC and 20 pct on Cmax) are Simcyp virtual-population output",
      "driven by demographic and enzyme-abundance distributions that are",
      "not published, so they are not encoded as omegas here."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter is fixed: this is a typical-value simulation model.
    # Values are either Huth 2019 Table 5 / Table 1 inputs, arithmetic
    # consequences of them, or (for the two distribution parameters a
    # whole-body PBPK never publishes) back-solved from the paper's own
    # printed summary statistics. Each case is labelled below.
    #
    # Shared quantities, all Huth 2019 Table 5 unless noted:
    #   MW    = 516.6 g/mol      B:P = 0.765      fu,plasma = 0.0002
    #   f_a   = 0.91             fraction available from the dosage form
    #   ka    = 0.687 1/h        lag = 1.5 h (optimized)
    #   Vss   = 1.45 L/kg on 80.08 kg = 116.116 L   (observed, footnote f)
    #   CL    = 3.12 L/h                            (observed, footnote g)
    #   renal CL = 0, biliary CL = 0: clearance is exclusively hepatic
    #                                 CYP-mediated metabolism
    #   fm CYP2C9 / CYP3A4, CYP2C9*1/*1 = 0.804 / 0.175
    #
    # One un-printed standard physiological constant is used: hepatic
    # blood flow QH = 90 L/h, and only to split the baseline oral
    # bioavailability into its gut and hepatic components. Siponimod is a
    # very low-extraction drug (eh = 3.12 / 90 = 0.0347), so the split is
    # nearly immaterial: varying QH by plus or minus 20 pct moves the
    # largest perpetrator-driven bioavailability change by under 1 pct.
    # ------------------------------------------------------------------

    # -- Absorption ----------------------------------------------------
    lka <- fixed(log(0.687))
    label("First-order absorption rate constant ka (1/h)")
    # Huth 2019 Table 5 row 'Absorption rate constant (1/hour)' = 0.687,
    # footnote d: estimated based on the population pharmacokinetic model.

    ltlag <- fixed(log(1.5))
    label("Absorption lag time (h)")
    # Huth 2019 Table 5 row 'Lag time (h)' = 1.5 (optimized), Methods:
    # optimized against the clinical Tmax at steady state of the
    # multiple-ascending-dose study.

    lfdepot <- fixed(log(0.91))
    label("Fraction of the oral dose available from the dosage form (unitless)")
    # Huth 2019 Table 5 row 'Fraction available from dosage form' = 0.91,
    # footnote c: based on the amount of siponimod excreted to feces in
    # the human mass-balance study. This is f_a, the ceiling on oral
    # bioavailability BEFORE first-pass loss; model() multiplies it by the
    # gut and hepatic availabilities to give the delivered f(depot).

    # -- Distribution --------------------------------------------------
    lvc <- fixed(log(43.72461))
    label("Central volume of distribution Vc (L)")
    # BACK-SOLVED, not transcribed. A full-PBPK model publishes no central
    # volume; Huth 2019 prints only the aggregate Vss. Vc and lq below
    # were solved jointly so the reduction reproduces the two Table 1
    # single-oral-dose shape statistics exactly: Tmax 3.45 h and Cmax
    # 7.76 ng/mL per mg (the ladder is strictly dose proportional from 0.1
    # to 75 mg). Held-out corroboration: the Table 1 intravenous arm, a
    # 2.92 h infusion of 0.25 mg, then gives a predicted Cmax of
    # 2.99 ng/mL against the paper's own predicted 3.02 (-1.1 pct).
    # A one-compartment reduction was tested first and rejected: with
    # V = Vss it over-predicts the printed Tmax by 86 pct.

    lvp <- fixed(log(72.39139))
    label("Peripheral volume of distribution Vp (L)")
    # NOT free. Fixed by the printed Vss once Vc is known:
    # Vp = 1.45 L/kg * 80.08 kg - 43.72461 = 116.116 - 43.725 = 72.391 L.

    lq <- fixed(log(30.02085))
    label("Intercompartmental clearance Q (L/h)")
    # BACK-SOLVED jointly with lvc; see the lvc comment.

    # -- Elimination, split by metabolic route ---------------------------
    # Huth 2019 Table 5 footnote g: the retrograde calculator was driven
    # by the observed geometric mean CL of 3.12 L/h, and Table 5 section 4
    # distributes it across the enzymes by their in-vitro fractions
    # metabolised. The renal and biliary clearance rows are both 0, so
    # every arm below is hepatic metabolism.
    lcl <- fixed(log(3.12))
    label("Total clearance CL in the CYP2C9*1/*1 reference genotype (L/h)")
    # Huth 2019 Table 1, intravenous 0.25 mg observed AUCinf of
    # 80.1 ng*h/mL gives 0.25 mg / 80.1 = 3.12 L/h, the value Table 5
    # footnote g names as the retrograde-calculator input.

    lcl_2c9 <- fixed(log(2.508480))
    label("CYP2C9-mediated clearance arm in the CYP2C9*1/*1 reference genotype (L/h)")
    # 3.12 * 0.804, where 0.804 is fm CYP2C9 for the CYP2C9*1/*1 row of
    # Huth 2019 Table 5 section 4.

    lcl_3a4 <- fixed(log(0.546000))
    label("CYP3A4-mediated clearance arm (L/h)")
    # 3.12 * 0.175, where 0.175 is fm CYP3A4 for the CYP2C9*1/*1 row of
    # Huth 2019 Table 5 section 4. Genotype-invariant by construction:
    # CYP2C9 genotype does not alter CYP3A4 activity, and it is exactly
    # this invariance that converts the paper's per-genotype fm pairs into
    # per-genotype total clearances (see the e_cyp2c9_* block).
    # model() forms the residual arm as lcl minus these two, i.e.
    # 3.12 * 0.021 = 0.0655 L/h, matching the CYP2B6 + CYP2C8 + CYP2C19
    # fractions of 0.004 + 0.017 + 0.001 = 0.022 printed in the same
    # table to within rounding.

    # -- Baseline first-pass extraction ----------------------------------
    # F = f_a * (1 - egut) * (1 - eh). Huth 2019 Table 1 gives the model's
    # own predicted absolute bioavailability directly, as the ratio of the
    # predicted oral and intravenous AUCinf at 0.25 mg:
    # F = 69.7 / 82.1 = 0.84897. With f_a = 0.91 and eh = 3.12 / 90 the
    # gut term is pinned with no freedom left.
    eh <- fixed(0.03466667)
    label("Hepatic extraction ratio in the CYP2C9*1/*1 reference genotype (unitless)")
    # 3.12 L/h / 90 L/h. model() rescales this with the current clearance,
    # so the hepatic blood flow cancels and only the ratio is load-bearing.

    egut <- fixed(0.03356880)
    label("Gut-wall extraction ratio (unitless)")
    # 1 - 0.84897 / (0.91 * (1 - 0.03466667)). Siponimod gut metabolism is
    # taken as entirely CYP3A4, the standard assumption for enterocyte
    # first pass; CYP2C9 abundance in the gut wall is negligible.

    # -- CYP2C9 diplotype effect on the CYP2C9 clearance arm --------------
    # Huth 2019 Table 5 section 4 prints an (fm CYP2C9 / fm CYP3A4) pair
    # for each of the six diplotypes. Because the CYP3A4 arm is
    # genotype-invariant, each pair fixes that genotype's total clearance
    # as CL_g = 0.546 / fm CYP3A4,g, and its CYP2C9 arm as the remainder
    # after the invariant CYP3A4 and residual arms are removed. Each
    # coefficient below is the log of that arm relative to CYP2C9*1/*1.
    #
    #  genotype  fm2C9  fm3A4   CL_g (L/h)  relative CYP2C9 activity
    #   *1/*1    0.804  0.175     3.120           1        (reference)
    #   *1/*2    0.788  0.189     2.889           0.907868
    #   *1/*3    0.678  0.287     1.902           0.514622
    #   *2/*2    0.727  0.244     2.238           0.648275
    #   *2/*3    0.616  0.343     1.592           0.390801
    #   *3/*3    0.074  0.822     0.664           0.021014
    #
    # Two independent checks on this derivation, neither used to fit it:
    # (i) the residual fraction implied by each CL_g closes that row of
    #     Table 5 to 1.000 (for example *3/*3: 0.074 + 0.822 + 0.099);
    # (ii) the resulting CL_g reproduce the six Table 2 predicted AUCinf
    #     ratios to within 5 pct and the six predicted terminal
    #     half-lives to within 6 pct, over a fivefold exposure range.
    e_cyp2c9_12_cl <- fixed(log(0.907868))
    label("log relative CYP2C9 clearance activity, CYP2C9*1/*2 versus *1/*1")
    e_cyp2c9_13_cl <- fixed(log(0.514622))
    label("log relative CYP2C9 clearance activity, CYP2C9*1/*3 versus *1/*1")
    e_cyp2c9_22_cl <- fixed(log(0.648275))
    label("log relative CYP2C9 clearance activity, CYP2C9*2/*2 versus *1/*1")
    e_cyp2c9_23_cl <- fixed(log(0.390801))
    label("log relative CYP2C9 clearance activity, CYP2C9*2/*3 versus *1/*1")
    e_cyp2c9_33_cl <- fixed(log(0.021014))
    label("log relative CYP2C9 clearance activity, CYP2C9*3/*3 versus *1/*1")

    # -- Perpetrator relative CYP3A4 activity -----------------------------
    # Each coefficient is the log of the relative CYP3A4 activity for that
    # perpetrator, back-solved from that arm's SIX published steady-state
    # AUC ratios in Huth 2019 Table 3 -- one free number against six
    # printed ones. The published Cmax ratios were held out entirely.
    # The recovered activities ladder monotonically with the paper's own
    # FDA potency classification, which is a free consistency check:
    #   ketoconazole 0.084 < itraconazole 0.209 (both strong CYP3A4)
    #     < erythromycin 0.312 (moderate) < fluconazole 0.339 (moderate)
    #     < fluvoxamine 0.909 (weak) < 1
    #     < efavirenz 2.870 (moderate inducer) < rifampicin 4.679 (strong).
    e_conmed_itraconazole_cyp3a4 <- fixed(log(0.208949))
    label("log relative CYP3A4 activity with itraconazole")
    e_conmed_ketoconazole_cyp3a4 <- fixed(log(0.084353))
    label("log relative CYP3A4 activity with ketoconazole")
    e_conmed_erythromycin_cyp3a4 <- fixed(log(0.312156))
    label("log relative CYP3A4 activity with erythromycin")
    e_conmed_fluconazole_cyp3a4 <- fixed(log(0.339165))
    label("log relative CYP3A4 activity with fluconazole")
    e_conmed_fluvoxamine_cyp3a4 <- fixed(log(0.908523))
    label("log relative CYP3A4 activity with fluvoxamine")
    e_conmed_rifampicin_cyp3a4 <- fixed(log(4.679147))
    label("log relative CYP3A4 activity with rifampicin")
    e_conmed_efv_cyp3a4 <- fixed(log(2.870164))
    label("log relative CYP3A4 activity with efavirenz")

    # -- Perpetrator relative CYP2C9 activity -----------------------------
    # Only the three perpetrators Huth 2019 classifies as dual
    # CYP3A4/CYP2C9 modulators carry a CYP2C9 term; for the other four the
    # relative CYP2C9 activity is held at exactly 1, matching the paper's
    # own arm labels in the Table 3 header and the Methods statement that
    # the efavirenz model was built as a pure CYP3A4 inducer. Holding
    # those four at 1 costs nothing: the worst AUC-ratio residual over all
    # 42 cells is 4.9 pct either way.
    e_conmed_fluconazole_cyp2c9 <- fixed(log(0.509023))
    label("log relative CYP2C9 activity with fluconazole")
    e_conmed_fluvoxamine_cyp2c9 <- fixed(log(0.665362))
    label("log relative CYP2C9 activity with fluvoxamine")
    e_conmed_rifampicin_cyp2c9 <- fixed(log(2.181306))
    label("log relative CYP2C9 activity with rifampicin")

    propSd <- fixed(0)
    label("Proportional residual error (none reported by the source)")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Typical-value structural parameters. No random effects.
    # ------------------------------------------------------------------
    ka <- exp(lka)
    tlag <- exp(ltlag)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)

    # ------------------------------------------------------------------
    # 2. CYP2C9 diplotype indicators built from the three allele counts,
    # which sum to 2 for every subject. All zero would be an ungenotyped
    # subject; CYP2C9*1/*1 needs no indicator because it is the reference.
    # ------------------------------------------------------------------
    is_cyp2c9_12 <- (CYP2C9_S1_COUNT == 1) * (CYP2C9_S2_COUNT == 1)
    is_cyp2c9_13 <- (CYP2C9_S1_COUNT == 1) * (CYP2C9_S3_COUNT == 1)
    is_cyp2c9_22 <- (CYP2C9_S2_COUNT == 2)
    is_cyp2c9_23 <- (CYP2C9_S2_COUNT == 1) * (CYP2C9_S3_COUNT == 1)
    is_cyp2c9_33 <- (CYP2C9_S3_COUNT == 2)

    a2c9_geno <- exp(
      e_cyp2c9_12_cl * is_cyp2c9_12 +
        e_cyp2c9_13_cl * is_cyp2c9_13 +
        e_cyp2c9_22_cl * is_cyp2c9_22 +
        e_cyp2c9_23_cl * is_cyp2c9_23 +
        e_cyp2c9_33_cl * is_cyp2c9_33
    )

    # ------------------------------------------------------------------
    # 3. Relative CYP3A4 and CYP2C9 activity contributed by a
    # coadministered perpetrator. All indicators zero (siponimod alone)
    # gives both activities equal to 1; the terms are additive on the log
    # scale so the flags are multiplicative on activity.
    # ------------------------------------------------------------------
    a3a4 <- exp(
      e_conmed_itraconazole_cyp3a4 * CONMED_ITRACONAZOLE +
        e_conmed_ketoconazole_cyp3a4 * CONMED_KETOCONAZOLE +
        e_conmed_erythromycin_cyp3a4 * CONMED_ERYTHROMYCIN +
        e_conmed_fluconazole_cyp3a4 * CONMED_FLUCONAZOLE +
        e_conmed_fluvoxamine_cyp3a4 * CONMED_FLUVOXAMINE +
        e_conmed_rifampicin_cyp3a4 * CONMED_RIFAMPICIN +
        e_conmed_efv_cyp3a4 * CONMED_EFV
    )
    a2c9 <- exp(
      e_conmed_fluconazole_cyp2c9 * CONMED_FLUCONAZOLE +
        e_conmed_fluvoxamine_cyp2c9 * CONMED_FLUVOXAMINE +
        e_conmed_rifampicin_cyp2c9 * CONMED_RIFAMPICIN
    )

    # ------------------------------------------------------------------
    # 4. Clearance. The residual arm (CYP2B6 + CYP2C8 + CYP2C19) is the
    # part of the reference clearance neither CYP2C9 nor CYP3A4 accounts
    # for; it is untouched by genotype and by every perpetrator the paper
    # simulates.
    # ------------------------------------------------------------------
    cl_ref <- exp(lcl)
    cl_2c9 <- exp(lcl_2c9)
    cl_3a4 <- exp(lcl_3a4)
    cl_other <- cl_ref - cl_2c9 - cl_3a4
    cl <- a2c9 * a2c9_geno * cl_2c9 + a3a4 * cl_3a4 + cl_other

    # ------------------------------------------------------------------
    # 5. Oral bioavailability. Multiplying an intrinsic clearance by a
    # factor A takes an extraction ratio E to E*A / (1 - E + E*A), so the
    # corresponding availability 1 - E is divided by (1 - E + E*A). The
    # gut wall is treated as entirely CYP3A4; the hepatic term simply
    # tracks the current total clearance, which is why the hepatic blood
    # flow cancels out of the ratio.
    # ------------------------------------------------------------------
    gfac <- 1 - egut + egut * a3a4
    fdepot <- exp(lfdepot) * ((1 - egut) / gfac) * (1 - eh * cl / cl_ref)

    # ------------------------------------------------------------------
    # 6. ODE system.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central +
      k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag
    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # 7. Observation. Doses are in mg and vc is in L, so central / vc is
    # in mg/L = ug/mL; multiply by 1000 to report ng/mL, the unit used
    # throughout Huth 2019 Tables 1 and 2.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
