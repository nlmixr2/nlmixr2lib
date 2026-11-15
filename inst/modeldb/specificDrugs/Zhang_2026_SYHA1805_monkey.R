Zhang_2026_SYHA1805_monkey <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey). One-compartment intravenous and oral",
    "pharmacokinetic reduction of the Simcyp (version 20) full-PBPK model",
    "for the non-steroidal farnesoid X receptor (FXR) agonist SYHA1805 in",
    "cynomolgus monkeys (Zhang 2026), the species from which the human",
    "clearance and volume were scaled. The source model's whole-body",
    "mass-balance equations and Rodgers-Rowland (method 2) tissue",
    "partition coefficients are platform outputs that the paper does not",
    "print, so the platform model itself cannot be encoded here. What IS",
    "reported is the monkey drug layer (Table 3): the in vivo intravenous",
    "clearance input of 23.46 mL/min, the predicted steady-state volume of",
    "0.23 L/kg, first-order absorption with ka 0.30 1/h and fraction",
    "absorbed 0.02, blood-to-plasma ratio 1.00 and fu 0.002. The parameters",
    "are absolute values for a 2.81 kg monkey, the body weight implied both",
    "by the clearance input (23.46 / 8.35 mL/min/kg) and by the paper's",
    "monkey hepatic blood flow (122.5 / 43.6 mL/min/kg). Oral",
    "bioavailability is fa x Fg x Fh with Fg = 1 and Fh from the",
    "well-stirred liver at 122.5 mL/min, giving 0.0162. No parameter is",
    "fitted here. The reduction reproduces the PBPK model's own predicted",
    "oral AUC (Figure 3D) to within 4 percent and its predicted oral Cmax",
    "(Figure 3F) to within 15 percent. This is a typical-value simulation",
    "model: the source reports no inter-individual variance components and",
    "no residual-error model, so there are no etas and propSd is fixed at",
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

  # What each ODE state holds. Zhang 2026 section 2.6.1 and Table 3: the
  # monkey model adapts the rat 'Full PBPK Model' with species-specific Kp
  # (scalar 3.00) and clearance; only the lumped predicted Vss is printed.
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

  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. Zhang 2026 Table 3 gives the monkey Vss in L/kg and",
        "the clearance input as an absolute 23.46 mL/min. The reduction",
        "fixes the 2.81 kg weight that the paper's own numbers imply",
        "(23.46 mL/min / 8.35 mL/min/kg from Table 2; 122.5 mL/min /",
        "43.6 mL/min/kg from Supplementary Table S1 and Equation S10) and",
        "carries no weight term. Doses given in mg/kg must be multiplied by",
        "2.81 kg."
      ),
      units = "kg",
      type = "continuous",
      notes = paste(
        "Implicit in the L/kg volume input and in the absolute clearance;",
        "not carried. No monkey body weight is printed directly."
      )
    )
  )

  population <- list(
    species = "cynomolgus monkey",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "Not reported",
    weight_median = "2.81 kg (implied by the absolute clearance input and the per-kg clearance; not printed)",
    sex_female_pct = 50,
    disease_state = "Healthy cynomolgus monkeys.",
    dose_range = paste(
      "Single intravenous dose of 1 mg/kg, or single oral doses of 10, 30",
      "and 100 mg/kg (4 groups of 6, 3 per sex; section 2.3.2)."
    ),
    route = "intravenous; oral",
    regions = "China (Hainan Primate Research Center animals; WuXi AppTec).",
    notes = paste(
      "This is a PBPK model rather than a population-PK fit. The monkey",
      "model adapted the rat structure with species-specific Kp and the",
      "observed monkey intravenous clearance (Table 2: 8.35 mL/min/kg,",
      "Vdss 0.216 L/kg, t1/2 3.30 h after 1 mg/kg IV), and kept the rat",
      "ka (0.30 1/h) and blood-to-plasma ratio (1.00). It was checked",
      "against the intravenous and oral monkey data (Figure 3D-F). Observed",
      "oral bioavailability was very low (Table 2: 1.58%, 1.27% and 0.989%",
      "at 10, 30 and 100 mg/kg), which the model carries through its fa of",
      "0.02. Plasma protein binding was above 99.8% (Supplementary Table",
      "S5). The companion human reduction is Zhang_2026_SYHA1805_human."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter below is fixed: nothing was estimated in building
    # this reduction. Values are verbatim Zhang 2026 entries or arithmetic
    # consequences of them at the 2.81 kg body weight the paper's own
    # clearance numbers imply.
    #
    # WHY THIS IS A ONE-COMPARTMENT MODEL. Table 3 prints a 'Full PBPK
    # Model' distribution with 'Method 2' Kp and a Kp scalar of 3.00, but
    # no tissue Kp, organ volume or blood flow; the only lumped
    # distribution quantity is the predicted Vss.
    # ------------------------------------------------------------------

    # Table 3, 'ka (1/h)', Oral (Monkey) = 0.30 (carried over from the rat
    # oral model; section 3.5).
    lka <- fixed(log(0.3))
    label("Absorption rate constant ka (1/h)")

    # Table 3, 'CLi.v (mL/min)', I.V. and Oral (Monkey) = 23.46 mL/min
    # = 23.46 * 60 / 1000 = 1.4076 L/h (Table 2 observed CL 8.35
    # mL/min/kg x 2.81 kg).
    lcl <- fixed(log(1.4076))
    label("Systemic clearance CL (L/h)")

    # Table 3, 'Vss input type', I.V. and Oral (Monkey) = 'Predicted
    # (0.23 L/kg)'. At 2.81 kg: 0.23 * 2.81 = 0.6463 L.
    lvc <- fixed(log(0.6463))
    label("Central volume vc (L), the lumped Vss at 2.81 kg")

    # Oral bioavailability F = fa * Fg * Fh, derived from printed values:
    #   fa = 0.02 (Table 3 'Fa', Oral (Monkey));
    #   Fg = 1 (no gut-wall clearance is specified; Table 3 'Fu(Gut)' = 1.00);
    #   Fh = 1 - CL / (Qh * B/P) = 1 - 23.46 / (122.5 * 1.00) = 0.80849,
    #        with B/P = 1.00 (Table 3) and the monkey hepatic blood flow
    #        Qh = 122.5 mL/min (Supplementary Equation S10 definitions),
    #        treating all of CL as hepatic (no renal clearance is printed).
    #   F = 0.02 * 0.80849 = 0.01617.
    lfdepot <- fixed(log(0.01617))
    label("Oral bioavailability F = fa * Fg * Fh (fraction)")

    # No residual-error model and no inter-individual variance components
    # are reported, so the residual error is fixed at zero.
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

    # 3. ODE system. Oral doses go into depot; intravenous doses go
    #    straight into central. Amounts in mg.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 4. Bioavailability on the oral dose only.
    f(depot) <- exp(lfdepot)

    # 5. Observation in ng/mL (mg / L * 1000).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
