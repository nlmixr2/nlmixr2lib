Aston_2017_efalizumab_qsp <- function() {
  description <- paste(
    "QSP. Five-state PK / receptor / feedback model of efalizumab and CD11a",
    "down-modulation in moderate-to-severe psoriasis (Example 2, Sect. 7.2 of",
    "Aston et al. 2017, reproducing the model of Ng et al. 2005). A",
    "two-compartment PK backbone with first-order subcutaneous absorption and",
    "PARALLEL linear plus Michaelis-Menten elimination drives saturable",
    "down-modulation of the total CD11a on the T-cell surface; the CD11a",
    "production rate is itself a dynamic state that relaxes toward a",
    "hyperbolic, negative-feedback function of the CD11a level. That slow",
    "feedback is what makes total CD11a REBOUND to about 140% of baseline",
    "around day 50 after a single 3 mg/kg intravenous dose, which the same",
    "model without feedback does not do. The whole model is BODY-WEIGHT",
    "NORMALISED: every drug amount is per kg and Vc is 64.3 mL/kg. CD11a is",
    "carried as a percentage of its own baseline, so total_target starts at",
    "100. Deterministic illustration: no subjects were fitted in Aston et al.,",
    "so there is no inter-individual variability and no residual-error model.",
    sep = " "
  )
  reference <- paste(
    "Aston PJ, Derks G, Agoram BM, van der Graaf PH.",
    "A mathematical analysis of rebound in a target-mediated drug disposition",
    "model: II. With feedback.",
    "J Math Biol. 2017;75(1):39-73. doi:10.1007/s00285-016-1073-6.",
    "Companion model from the same paper: modellib('Aston_2017_omalizumab_qsp').",
    "The model structure and parameter values originate with",
    "Ng CM, Joshi A, Dedrick RL, Garovoy MR, Bauer RJ.",
    "Pharmacokinetic-pharmacodynamic-efficacy analysis of efalizumab in patients",
    "with moderate to severe psoriasis. Pharm Res. 2005;22(7):1088-1100.",
    "doi:10.1007/s11095-005-5642-4;",
    "Ng et al. (2005) was NOT available when this model was built, and every value below",
    "is transcribed from Sect. 7.2 of Aston et al. (2017), including Aston's own",
    "ten-fold correction to koff. See the vignette Errata.",
    sep = " "
  )
  vignette <- "Aston_2017_receptor_rebound"

  # Time is days. Every drug state is an AMOUNT PER KILOGRAM of body weight
  # (ug/kg), because the model is body-weight normalised throughout: Vc is
  # quoted as 64.3 mL/kg, Vm as a per-kg rate, and the illustrative dose as
  # 3 mg/kg. The `dosing` slot therefore records "ug" with the per-kg
  # normalisation carried in this comment and in the vignette, exactly as
  # Penney_2025_tce_tmdd_qsp.R records "nmol" for a mg/kg-dosed model.
  # Concentration is ug/mL, the units of Kmc and of the reported affinity.
  units <- list(time = "day", dosing = "ug", concentration = "ug/mL")

  # X3 is the TOTAL (free + drug-bound) CD11a on the T-cell surface, and the
  # bound species is not carried as a separate state, which is exactly the
  # canonical `total_target` role. X4 is the CD11a PRODUCTION RATE: a
  # first-order delay driven by a system state (total CD11a) with no mass
  # transfer, whose value is the production rate it modulates - the canonical
  # Gabrielsson-Hjorth `moderator1` role.
  compartmentData <- list(
    depot = list(analyte = "efalizumab", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "efalizumab", units = "ug", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "efalizumab", units = "ug", specimen = "serum", verified = TRUE),
    total_target = list(
      analyte = "CD11a on the T-cell surface, total (free + bound), as a percentage of baseline",
      units = "percent of baseline",
      specimen = "whole blood",
      verified = TRUE
    ),
    moderator1 = list(
      analyte = "CD11a production rate (the feedback moderator)",
      units = "percent of baseline per day",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    disease_state = "moderate-to-severe plaque psoriasis",
    dose_range = paste(
      "a single 3 mg/kg intravenous dose is the illustration of Fig. 11;",
      "the subcutaneous route is available through the depot compartment",
      "(Fa = 0.564)"
    ),
    notes = paste(
      "Aston et al. (2017) fitted no data: Sect. 7.2 re-simulates the published",
      "model of Ng et al. (2005) in order to demonstrate that its observed CD11a",
      "rebound is caused by the slow feedback on CD11a production. The subject",
      "count and demographics of the underlying psoriasis studies are reported in",
      "Ng et al. (2005), which was not available when this model was built, so they are",
      "left NA rather than guessed. No inter-individual variability and no",
      "residual-error model are reported anywhere in Aston et al., so every",
      "parameter is encoded with fixed(). The paper's own numerical anchors for",
      "this example are epsilon = koff / Vm2 = 7.13e-3, the slowest eigenvalue of",
      "the receptor-feedback block (-1.37e-2 in the non-dimensional time of",
      "tau = Vm2 * t), the slowest eigenvalue of the PK block (-8.87e-2 in the",
      "same non-dimensional time), and the Fig. 11 time courses; the vignette",
      "reproduces all of them."
    )
  )

  ini({
    # ---- PK disposition (the dX_sc, dX_1, dX_2 equations of Sect. 7.2) --------
    # All rates are first-order constants, as the source parameterises them; the
    # source reports no separate CL or Q.
    lka <- fixed(log(0.242)); label("First-order subcutaneous absorption rate constant ka (1/day)") # Sect. 7.2 (Ng 2005 Tables II-III): ka = 0.242 /day
    lkel <- fixed(log(0.114)); label("Linear central elimination rate constant k10 (1/day)") # Sect. 7.2 (Ng 2005 Tables II-III): k10 = 0.114 /day
    lk12 <- fixed(log(0.097)); label("Central-to-peripheral rate constant k12 (1/day)") # Sect. 7.2 (Ng 2005 Tables II-III): k12 = 0.097 /day
    lk21 <- fixed(log(0.193)); label("Peripheral-to-central rate constant k21 (1/day)") # Sect. 7.2 (Ng 2005 Tables II-III): k21 = 0.193 /day
    lfdepot <- fixed(log(0.564)); label("Subcutaneous bioavailability Fa (fraction)") # Sect. 7.2 (Ng 2005 Tables II-III): Fa = 0.564
    # Kept in mL/kg exactly as the source prints it. That is also the only
    # spelling that keeps the model dimensionally consistent: `central` is an
    # amount per kg (ug/kg), so `central / vc` is ug/mL, and `km * vc` is
    # 0.033 ug/mL * 64.3 mL/kg = 2.1219 ug/kg, which matches the units of
    # `central` in the Michaelis-Menten denominator. Writing Vc as 0.0643 L/kg
    # instead would make `km * vc` a thousand-fold too small.
    lvc <- fixed(log(64.3)); label("Central volume of distribution Vc (mL/kg)") # Sect. 7.2 (Ng 2005 Tables II-III): Vc = 64.3 mL/kg

    # ---- Saturable (Michaelis-Menten) elimination -----------------------------
    # Vm is printed in Sect. 7.2 as '26.9 ug/mL', which cannot be right
    # dimensionally: Vm is the saturating maximum of the term
    # Vm * X1 / (Kmc * Vc + X1) added to dX1/dt, so it must carry the units of an
    # amount per unit time, here ug/kg/day. The NUMBER 26.9 is confirmed by the
    # paper's own printed eigenvalue: the slowest PK eigenvalue quoted in
    # Sect. 7.2 as -8.87e-2 (non-dimensional, tau = Vm2 * t) is reproduced to
    # three significant figures only with Vm / (Kmc * Vc) = 26.9 / 2.1219 =
    # 12.68 /day. See the vignette Errata and its numerical demonstration.
    lvmax <- fixed(log(26.9)); label("Maximum rate of saturable central elimination Vm (ug/kg/day)") # Sect. 7.2 (Ng 2005 Tables II-III): Vm = 26.9, printed with the wrong units
    # Kmc enters BOTH the saturable elimination and the CD11a down-modulation as
    # Kmc * Vc, i.e. as an amount per kg. Aston notes it equals the affinity of
    # 0.033 ug/mL that Ng et al. report.
    km <- fixed(0.033); label("Michaelis constant of saturable elimination and of CD11a down-modulation Kmc (ug/mL)") # Sect. 7.2 (Ng 2005 Tables II-III): Kmc = 0.033 ug/mL

    # ---- CD11a turnover and drug-driven down-modulation -----------------------
    kdeg <- fixed(0.444); label("First-order CD11a loss rate constant k30 (1/day)") # Sect. 7.2 (Ng 2005 Tables II-III): k30 = 0.444 /day
    # Vm2 scales a FIRST-ORDER term in the CD11a equation
    # (-Vm2 * X3 * occupancy), so it is a maximal rate CONSTANT, not a
    # Michaelis-Menten Vmax; per the parameter register, such a quantity takes a
    # kmax-style name rather than vmax.
    kmax_total_target <- fixed(2.16); label("Maximum drug-driven first-order CD11a loss rate Vm2 (1/day)") # Sect. 7.2 (Ng 2005 Tables II-III): Vm2 = 2.16 /day
    kinmax <- fixed(334); label("Maximum zero-order CD11a production rate k03max (percent of baseline per day)") # Sect. 7.2 (Ng 2005 Tables II-III): k03max = 334 %CD11a/day

    # ---- Feedback moderator ---------------------------------------------------
    # Aston writes the moderator relaxation rate as koff, a name inherited from
    # Ng et al.; structurally it is the RESPONSE RATE of the CD11a-production
    # feedback, dX4/dt = koff * (k03max * H(X3) - X4), so it takes the canonical
    # moderator-chain turnover name ktol here. Sect. 7.2 states that Ng et al.'s
    # printed value of 0.00154 /day did NOT reproduce their own Fig. 3B and that
    # 0.0154 /day did; 0.0154 is the value Aston et al. use for Fig. 11 and for
    # the printed epsilon = koff / Vm2 = 7.13e-3, so it is the value encoded
    # here. Setting ktol to 0.00154 recovers the smaller rebound (about 110% of
    # baseline) that Aston reports for the uncorrected value.
    ktol <- fixed(0.0154); label("CD11a-production feedback response rate koff (1/day)") # Sect. 7.2: Aston's corrected value; Ng 2005 printed 0.00154 /day

    # The half-saturation constant Kmc03 of the feedback function
    # H(X3) = Kmc03 / (Kmc03 + X3) is the ONE parameter of this model that
    # Aston et al. never print. It is not guessed: CD11a is expressed as a
    # PERCENTAGE OF ITS OWN BASELINE, so the drug-free steady state of the CD11a
    # state is 100 by construction - which is what the Fig. 11 axis
    # ('%CD11a/baseline', both panels starting at and the no-feedback panel
    # returning to 100) shows. Given that baseline, Aston's Eq. (72) determines
    # Kmc03 exactly, and model() derives it below rather than hard-coding it.
    # The derived value, 15.33 %CD11a, is corroborated independently: it is the
    # only value that reproduces the paper's printed slowest feedback-block
    # eigenvalue of -1.37e-2. See the vignette Errata.
    bl_total_target <- fixed(100); label("Baseline total CD11a on the T-cell surface (percent of baseline, = 100 by definition)") # Fig. 11 axis '%CD11a/baseline'; CD11a is reported as a percentage of its own baseline

    # The source is a deterministic mathematical analysis and reports no
    # residual-error model and no inter-individual variability.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported in source)")
  })

  model({
    ka <- exp(lka)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    fdepot <- exp(lfdepot)
    vc <- exp(lvc)
    vmax <- exp(lvmax)

    # Drug-free steady state of the CD11a / feedback pair. Writing F = X4/k03max
    # as Aston does, Eq. (72) is k03max * F^2 + k30 * Kmc03 * F - k30 * Kmc03 = 0
    # and the baseline is X3 = k03max * F0 / k30. Inverting on the baseline:
    f0 <- kdeg * bl_total_target / kinmax
    km_moderator1 <- kinmax * f0 * f0 / (kdeg * (1 - f0))

    total_target(0) <- bl_total_target
    moderator1(0) <- kdeg * bl_total_target

    # Fractional saturation of the drug-binding site, X1 / (Kmc * Vc + X1).
    # Kmc is a concentration and Vc a volume per kg, so Kmc * Vc is an amount
    # per kg and matches the units of `central`.
    occ <- central / (km * vc + central)

    # dX_sc/dt = -ka * X_sc
    d/dt(depot) <- -ka * depot

    # dX_1/dt = -(k10 + k12) * X_1 + k21 * X_2 - Vm * X_1 / (Kmc * Vc + X_1)
    #           + Fa * ka * X_sc
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1 -
      vmax * occ + fdepot * ka * depot

    # dX_2/dt = k12 * X_1 - k21 * X_2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # dX_3/dt = X_4 - k30 * X_3 - Vm2 * X_3 * X_1 / (Kmc * Vc + X_1)
    d/dt(total_target) <- moderator1 - kdeg * total_target -
      kmax_total_target * total_target * occ

    # dX_4/dt = koff * (k03max * Kmc03 / (Kmc03 + X_3) - X_4)
    d/dt(moderator1) <- ktol *
      (kinmax * km_moderator1 / (km_moderator1 + total_target) - moderator1)

    # Efalizumab serum concentration (ug/mL) from the per-kg central amount.
    Cc <- central / vc

    # Total CD11a is the state itself; the FREE CD11a that Fig. 11 also plots is
    # Kmc * Vc * X_3 / (Kmc * Vc + X_1) (final paragraph of Sect. 7.2).
    totalCD11a <- total_target
    freeCD11a <- km * vc * total_target / (km * vc + central)

    Cc ~ prop(propSd)
  })
}
