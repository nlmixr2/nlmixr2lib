Jeong_2022_nafamostat <- function() {
  description <- paste(
    "Four-compartment intravenous pharmacokinetic reduction of the Simcyp",
    "full-PBPK model for the serine-protease inhibitor nafamostat in healthy",
    "adults (Jeong 2022), with algebraic unbound-plasma and lung-tissue",
    "outputs. The source model was built in the Simcyp simulator (version",
    "21) with a full-PBPK Rodgers and Rowland distribution model, an",
    "estimated Kp scalar and enzyme-kinetic esterase elimination; its",
    "whole-body equations, organ volumes, blood flows and tissue partition",
    "coefficients are platform outputs that are not published, so the",
    "platform model itself cannot be encoded here. This reduction is",
    "anchored to two paper values: the total clearance implied by the",
    "paper's own predicted AUC0-inf (Table 2; dose / AUC0-inf = 220.8 L/h,",
    "identical at 10, 20 and 40 mg) and the predicted steady-state volume of",
    "distribution of 11.66 L/kg (Table 1) at a 70 kg reference weight. The",
    "remaining distribution constants were recovered by fitting the paper's",
    "published predicted mean plasma profiles after 2-h infusions of 10, 20",
    "and 40 mg (Figure 1, extracted from the vector graphics of the article",
    "PDF); one- to three-compartment reductions cannot reproduce the",
    "sub-minute post-infusion drop of those curves. Held out from the fit,",
    "the reduction reproduces the predicted Cmax and AUC0-last of Table 2 to",
    "within 1 percent and the predicted 13-day continuous-infusion plasma",
    "profiles of Figure 2 to within 6 percent. Unbound plasma concentration",
    "is fu * Cc with the Table 1 fu of 0.46, and lung-tissue concentration",
    "is a constant 5.44 times plasma, the ratio of the Figure 2 lung and",
    "plasma curves (constant throughout both regimens because the",
    "perfusion-limited lung equilibrates within minutes). The global",
    "sensitivity analysis of physiological parameters cannot be reproduced",
    "by a compartmental reduction. This is a typical-value simulation",
    "model: the source reports no inter-individual variance components and",
    "no residual-error model, so there are no etas and propSd is fixed at",
    "zero.",
    sep = " "
  )
  reference <- paste(
    "Jeong HC, Chae YJ, Shin KH. (2022). Predicting the systemic exposure",
    "and lung concentration of nafamostat using physiologically-based",
    "pharmacokinetic modeling. Transl Clin Pharmacol 30(4):201-211.",
    "doi:10.12793/tcp.2022.30.e20.",
    sep = " "
  )
  vignette <- "Jeong_2022_nafamostat"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(
      analyte = "nafamostat",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "nafamostat",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    peripheral2 = list(
      analyte = "nafamostat",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    peripheral3 = list(
      analyte = "nafamostat",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  # Body weight is implicit in the L/kg volume input but is never printed.
  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. Jeong 2022 Table 1 expresses Vss in L/kg, and the",
        "second COVID-19 regimen is dosed per kilogram (4.8 mg/kg/24 h).",
        "This reduction fixes a 70 kg reference weight instead of carrying",
        "a weight term, because the clearance it is anchored to (Table 2",
        "dose / AUC0-inf) is an absolute value for the unpublished Simcyp",
        "virtual population, and scaling volume alone with weight would",
        "make the model internally inconsistent."
      ),
      units = "kg",
      type = "continuous",
      notes = paste(
        "Implicit in the L/kg volume input; not carried. No body weight is",
        "printed anywhere in the paper. The 70 kg value is the standing",
        "rounded-standard assumption. The ratio of the two Figure 2 plateaus",
        "(1.470) implies a mean simulated weight near 61 kg for the",
        "4.8 mg/kg regimen; the vignette uses that weight only to convert",
        "the per-kilogram regimen into an infusion rate."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1000L,
    n_studies = 1L,
    age_range = "20-40 years (virtual healthy volunteers; the clinical verification study enrolled healthy Chinese adults aged 20-30 years)",
    weight_median = "70 kg (assumed reference weight for the L/kg volume input; not reported)",
    sex_female_pct = 50,
    race_ethnicity = "Simcyp virtual healthy-volunteer population; verification data from healthy Chinese adults",
    disease_state = "Healthy adult volunteers (model applied to COVID-19 dosing regimens, but no patient physiology was used).",
    dose_range = paste(
      "Verification: single 2-h intravenous infusions of 10, 20 and 40 mg.",
      "Application: continuous intravenous infusion of 200 mg/24 h or",
      "4.8 mg/kg/24 h for 13 days."
    ),
    route = "intravenous",
    regions = "Simcyp virtual population (Methods, Model simulation and verification).",
    studies = paste(
      "Observed plasma profiles digitised by the authors from Cao et al. 2008",
      "(reference 22; n = 10 per dose, healthy Chinese volunteers, 10, 20 and",
      "40 mg 2-h infusions). Each simulation was 10 virtual trials of 100",
      "subjects aged 20-40 years with a 50:50 sex ratio."
    ),
    notes = paste(
      "This is a PBPK simulation analysis rather than a population-PK fit,",
      "so there is no pooled analysis dataset and no estimated variance",
      "components; the 5th-95th percentile bands of Figures 1 and 2 are the",
      "spread of the virtual population. Elimination as built in the source",
      "(Table 1): CES2 enzyme kinetics (Vmax 26,900 pmol/min/mg protein, Km",
      "1,790 uM), plasma esterase half-life 0.63 min (estimated), additional",
      "HLM and HLC CLint of 16.96 and 73 uL/min/mg protein, renal clearance",
      "0.56 L/h (estimated) and additional systemic clearance 0.02 L/h",
      "(estimated). Those are recorded for provenance; turning them into a",
      "systemic clearance needs unpublished Simcyp system values (microsomal",
      "and cytosolic protein per gram of liver, liver weight, plasma volume),",
      "so the reduction is instead anchored to the clearance implied by the",
      "paper's predicted AUC0-inf. Physicochemical inputs: MW 347.378 g/mol,",
      "logP 2.52, monoprotic base pKa 11.32, B/P 1.19, fu 0.46, Kp scalar",
      "2.18 (estimated)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter below is fixed: nothing was estimated against
    # observed data. Three kinds of value appear, marked on each line:
    #   * paper value -- a Table 1 entry, or a direct arithmetic
    #     consequence of Table 2 (no assumption);
    #   * paper value at 70 kg -- the L/kg Vss of Table 1 times the
    #     standing 70 kg reference weight (no body weight is printed);
    #   * figure-derived -- recovered by the maintainers from the
    #     predicted mean curves of Figure 1 (three doses, fitted jointly,
    #     log-scale least squares, clearance and Vss held at the two
    #     anchors above) or read off the ratio of the Figure 2 curves.
    #     The curves were read from the vector paths of the article PDF,
    #     not by eye; the Figure 1 apexes reproduce the Table 2 predicted
    #     Cmax to within 0.6 percent.
    #
    # WHY FOUR COMPARTMENTS. The Figure 1 curves fall about 60 percent in
    # the first three minutes after the infusion stops and then decline
    # with half-lives near 0.4 h and 3 h, while the 11.66 L/kg Vss needs a
    # deep, slowly equilibrating space. Two- and three-compartment
    # reductions with the same two anchors miss Figure 1 with log-scale
    # RMSE 0.37 and 0.085 (versus 0.050 here) and miss the Table 2 Cmax by
    # 14 and 10 percent; see the vignette.
    # ------------------------------------------------------------------

    # paper value: Table 2 predicted AUC0-inf = 45.30, 90.60 and 181.00
    # ng*h/mL after 10, 20 and 40 mg; dose / AUC0-inf = 220.8, 220.8 and
    # 221.0 L/h. Lumps esterase, hepatic, renal (0.56 L/h) and additional
    # systemic (0.02 L/h) clearance.
    lcl <- fixed(log(220.8))
    label("Total plasma clearance CL (L/h)")

    # figure-derived (Figure 1 fit)
    lvc <- fixed(log(7.61))
    label("Central compartment volume vc (L)")

    # figure-derived (Figure 1 fit)
    lq <- fixed(log(104))
    label("Intercompartmental clearance to peripheral1, q (L/h)")

    # figure-derived (Figure 1 fit)
    lvp <- fixed(log(42.2))
    label("Peripheral1 volume vp (L)")

    # figure-derived (Figure 1 fit)
    lq2 <- fixed(log(64.8))
    label("Intercompartmental clearance to peripheral2, q2 (L/h)")

    # figure-derived (Figure 1 fit)
    lvp2 <- fixed(log(232))
    label("Peripheral2 volume vp2 (L)")

    # figure-derived (Figure 1 fit)
    lq3 <- fixed(log(7.57))
    label("Intercompartmental clearance to peripheral3, q3 (L/h)")

    # paper value at 70 kg: Table 1 'Vss (L/kg)' = 11.66, so
    # Vss = 816.2 L; vp3 = 816.2 - 7.61 - 42.2 - 232 = 534.4 L.
    lvp3 <- fixed(log(534.4))
    label("Peripheral3 volume vp3 (L)")

    # paper value: Table 1 'fu' = 0.46 (predicted in Simcyp).
    fu <- fixed(0.46)
    label("Fraction unbound in plasma (unitless)")

    # figure-derived: ratio of Figure 2C (lung) to Figure 2A (plasma) =
    # 5.440 at every plotted time over the 13-day infusion; equals the
    # Simcyp-predicted lung Kp times the 2.18 Kp scalar of Table 1.
    lkp_lung <- fixed(log(5.44))
    label("Lung tissue to plasma total concentration ratio (unitless)")

    # Jeong 2022 is a PBPK simulation analysis, not a population-PK fit,
    # and reports no residual-error model. Rather than invent a variance,
    # the residual error is fixed at zero.
    propSd <- fixed(0)
    label("Proportional residual error SD (fraction; zero, no error model reported by the source)")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters. No covariates and no random effects.
    # ------------------------------------------------------------------
    cl  <- exp(lcl)
    vc  <- exp(lvc)
    q   <- exp(lq)
    vp  <- exp(lvp)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2)
    q3  <- exp(lq3)
    vp3 <- exp(lvp3)
    kp_lung <- exp(lkp_lung)

    # ------------------------------------------------------------------
    # 2. Micro-constants.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    k14 <- q3 / vc
    k41 <- q3 / vp3

    # ------------------------------------------------------------------
    # 3. ODE system (amounts in mg). Intravenous administration only;
    #    infusions are given through the event table's rate or dur.
    # ------------------------------------------------------------------
    d/dt(central) <- -(kel + k12 + k13 + k14) * central +
      k21 * peripheral1 + k31 * peripheral2 + k41 * peripheral3
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(peripheral3) <- k14 * central - k41 * peripheral3

    # ------------------------------------------------------------------
    # 4. Outputs. Doses are in mg and vc in L, so central / vc is mg/L;
    #    multiply by 1000 for ng/mL, the units of Figure 1 and Table 2.
    #    Cu is the unbound plasma concentration and Clung the total lung
    #    tissue concentration of Figure 2 (both ng/mL; divide by the
    #    molecular weight 347.378 g/mol and multiply by 1000 for nM).
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cu <- fu * Cc
    Clung <- kp_lung * Cc
    Cc ~ prop(propSd)
  })
}
