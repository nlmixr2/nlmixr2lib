Zhang_2026_SYHA1805_human <- function() {
  description <- paste(
    "One-compartment oral pharmacokinetic reduction of the Simcyp",
    "(version 20) full-PBPK model for the non-steroidal farnesoid X",
    "receptor (FXR) agonist SYHA1805 in healthy adults (Zhang 2026), the",
    "model used to support first-in-human dose selection. The source",
    "model's whole-body mass-balance equations and its Rodgers-Rowland",
    "(method 2) tissue partition coefficients are platform outputs that the",
    "paper does not print, so the platform model itself cannot be encoded",
    "here. What IS reported is the drug layer of the human model (Table 3):",
    "the in vivo intravenous clearance input of 20.7 L/h (allometrically",
    "scaled from cynomolgus monkey), the predicted steady-state volume of",
    "0.21 L/kg, first-order absorption with ka 0.30 1/h and fraction",
    "absorbed 1.00, blood-to-plasma ratio 1.00 and fu 0.003. Oral",
    "bioavailability is fa x Fg x Fh with Fg = 1 and Fh from the",
    "well-stirred liver at the paper's human hepatic blood flow of",
    "1450 mL/min (Supplementary Equation S10), giving 0.762. No parameter",
    "is fitted here. The reduction reproduces the PBPK model's own",
    "predicted AUC at 30, 60 and 120 mg (Figure 4) to within 3 percent and",
    "its predicted Cmax to within 9 percent. The paper itself reports that",
    "this model over-predicts the observed single-ascending-dose exposure",
    "(1.8- to 2-fold on AUC). This is a typical-value simulation model: the",
    "source reports no inter-individual variance components and no",
    "residual-error model, so there are no etas and propSd is fixed at",
    "zero.",
    sep = " "
  )
  reference <- paste(
    "Zhang L, Zhang M, Zhou E, Zhang X, Qi H, Yao X, Liu D. (2026).",
    "First-in-Human Dose Selection,",
    "Pharmacokinetics Prediction, and Clinical Validation of SYHA1805, a",
    "Novel FXR Agonist, Using Allometric Scaling and PBPK Modeling.",
    "Pharmaceutics 18(7):862. doi:10.3390/pharmaceutics18070862.",
    sep = " "
  )
  vignette <- "Zhang_2026_SYHA1805"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # What each ODE state holds. Zhang 2026 section 2.6.1 and Table 3: 'Full
  # PBPK Model' distribution with Rodgers-Rowland 'Method 2' Kp and a Kp
  # scalar of 3; only the lumped predicted Vss is printed.
  compartmentData <- list(
    depot = list(
      analyte = "SYHA1805",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "SYHA1805",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  # Body weight is implicit in the L/kg volume input but is never used as
  # a covariate. Recorded as screened-but-not-carried.
  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. Zhang 2026 Table 3 expresses the human Vss in L/kg,",
        "while the clearance input (20.7 L/h) is an absolute value",
        "scaled to a 70 kg human (section 3.4: 'The predicted Vss,h",
        "(70 kg) was 15.1 L'). This reduction fixes the 70 kg reference",
        "weight rather than carrying a weight term, so that volume and",
        "clearance stay mutually consistent."
      ),
      units = "kg",
      type = "continuous",
      notes = paste(
        "Implicit in the L/kg volume input; not carried. The 70 kg",
        "reference weight is the paper's own allometric reference",
        "(section 3.4 and Supplementary Equations S7-S11)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 100L,
    n_studies = 1L,
    age_range = "20-50 years (virtual healthy volunteers; section 2.6.2)",
    weight_median = "70 kg (allometric reference weight; section 3.4)",
    sex_female_pct = 50,
    race_ethnicity = "Chinese (healthy Chinese volunteers; Conclusions)",
    disease_state = "Healthy adult volunteers.",
    dose_range = paste(
      "Single oral doses of 30, 60 and 120 mg (single-ascending-dose study",
      "CTR20202354; Figures 4 and 5). The proposed first-in-human range was",
      "30 mg starting dose to 3000 mg maximum, with a predicted effective",
      "dose of 1150 mg (section 3.4)."
    ),
    route = "oral",
    regions = "China (CTR20202354).",
    notes = paste(
      "This is a PBPK prediction rather than a population-PK fit. The",
      "human model was built bottom-up: clearance and Vss were",
      "allometrically scaled from cynomolgus monkey (Supplementary Table",
      "S7: monkey-based mean CLi.v 20.7 L/h, Vss 15.1 L; the rat-based",
      "values of 52.5 L/h and 125 L were not used), ka 0.30 1/h and the",
      "blood-to-plasma ratio 1.00 were carried over from the rat oral",
      "model, and fu 0.003 and the physicochemical inputs (log P 1.63,",
      "pKa 2.95, Caco-2 permeability 0.39e-6 cm/s) came from in vitro",
      "measurements (Table 3). Simulations used a Simcyp virtual population",
      "of 10 trials of 10 subjects (n = 100, 50% female, 20-50 years). The",
      "model was then checked against the single-ascending-dose study",
      "CTR20202354 at 30, 60 and 120 mg (Figures 4 and 5; sampling to",
      "72 h). The companion cynomolgus monkey reduction is",
      "Zhang_2026_SYHA1805_monkey; the rat PBPK model is not reproducible",
      "from the printed values (see the vignette)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter below is fixed: nothing was estimated in building
    # this reduction. Values are verbatim Zhang 2026 entries or arithmetic
    # consequences of them at the paper's own 70 kg reference weight.
    #
    # WHY THIS IS A ONE-COMPARTMENT MODEL. Table 3 prints a 'Full PBPK
    # Model' distribution with 'Method 2' Kp and a Kp scalar of 3.00, but
    # no tissue Kp, organ volume or blood flow, so the multi-tissue
    # structure cannot be rebuilt; its only lumped distribution quantity
    # is the predicted Vss. The reduction's agreement with the PBPK
    # model's own predicted AUC, Cmax and mean profile (Figures 4 and 5)
    # is shown in the vignette; nothing was tuned to reach it.
    # ------------------------------------------------------------------

    # Table 3, 'ka (1/h)', Oral (Human) = 0.30 (first-order absorption,
    # estimated in the rat oral model and carried to monkey and human;
    # section 2.6.1).
    lka <- fixed(log(0.3))
    label("Absorption rate constant ka (1/h)")

    # Table 3, 'CLi.v (mL/min)', Oral (Human) = 20.70. The human column is
    # in L/h, not the mL/min of the row label: section 3.4 and the
    # abstract give 'a human intravenous clearance of 20.7 L/h', which is
    # also the monkey-based mean of Supplementary Table S7.
    lcl <- fixed(log(20.7))
    label("Systemic clearance CL (L/h)")

    # Table 3, 'Vss input type', Oral (Human) = 'Predicted (0.21 L/kg)'.
    # At the 70 kg reference weight: 0.21 L/kg * 70 kg = 14.7 L (section
    # 3.4 gives the allometric Vss,h as 15.1 L at 70 kg, 2.7% higher).
    lvc <- fixed(log(14.7))
    label("Central volume vc (L), the lumped Vss at 70 kg")

    # Oral bioavailability F = fa * Fg * Fh, derived from printed values:
    #   fa = 1.00 (Table 3 'Fa', Oral (Human));
    #   Fg = 1 (no gut-wall clearance is specified; Table 3 'Fu(Gut)' = 1.00
    #        and elimination is the single in vivo clearance input);
    #   Fh = 1 - CL / (Qh * B/P) = 1 - 20.7 / (87.0 * 1.00) = 0.7621, with
    #        B/P = 1.00 (Table 3) and the human hepatic blood flow
    #        Qh = 1450 mL/min = 87.0 L/h (Supplementary Equation S10
    #        definitions), treating all of CL as hepatic (no renal
    #        clearance is printed).
    lfdepot <- fixed(log(0.7621))
    label("Oral bioavailability F = fa * Fg * Fh (fraction)")

    # Zhang 2026 is a PBPK simulation analysis, not a population-PK fit.
    # It reports no residual-error model and no inter-individual variance
    # components, so the residual error is fixed at zero and this is a
    # deterministic typical-value simulation model.
    propSd <- fixed(0)
    label("Proportional residual error SD (fraction; zero, no error model reported by the source)")
  })

  model({
    # 1. Individual parameters. No covariates and no random effects.
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)

    # 2. Micro-constant for elimination.
    kel <- cl / vc

    # 3. ODE system. First-order absorption from the gut lumen into a
    #    single plasma compartment; amounts in mg.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 4. Bioavailability on the oral dose.
    f(depot) <- exp(lfdepot)

    # 5. Observation. mg / L = ug/mL; multiply by 1000 for the ng/mL of
    #    Figures 4 and 5.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
