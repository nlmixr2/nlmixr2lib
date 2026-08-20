Sheiner_1979_tubocurarine <- function() {
  description <- "Two-compartment IV population PK linked to a hypothetical effect compartment and a sigmoid Emax pharmacodynamic model for d-tubocurarine (dTC) neuromuscular blockade in 20 adults undergoing elective surgery, 10 with normal renal function and 10 with chronic end-stage renal failure awaiting renal transplantation (Sheiner 1979, group 3). This is the founding effect-compartment paper: plasma concentration and effect (degree of thumb-adduction paralysis, 0 = none to 1 = complete) were fitted simultaneously by nonlinear mixed effects, and the effect compartment is parameterised so that its concentration is the equivalent steady-state plasma concentration Cpss. Renal failure is encoded exactly as the paper's baseline model does, as the absence of the renal elimination arm of clearance, with no other structural difference between the two groups."
  reference <- paste(
    "Sheiner LB, Stanski DR, Vozeh S, Miller RD, Ham J.",
    "Simultaneous modeling of pharmacokinetics and pharmacodynamics:",
    "application to d-tubocurarine.",
    "Clin Pharmacol Ther. 1979;25(3):358-371.",
    "doi:10.1002/cpt1979253358.",
    "Read from the full facsimile reprint bundled with",
    "Karlsson MO, In the cradle of pharmacometric methodology:",
    "introducing population PKPD modeling, simultaneous analysis and the",
    "effect-compartment model. Commentary on Sheiner et al.",
    "Clin Pharmacol Ther. 2025;117(6):1517-1532.",
    "doi:10.1002/cpt.3663 (PMC12087688), pages 1519-1532.",
    sep = " "
  )
  vignette <- "Sheiner_1979_tubocurarine"
  units    <- list(time = "min", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central     = list(analyte = "d-tubocurarine", units = "mg",    specimen = "plasma",         verified = TRUE),
    peripheral1 = list(analyte = "d-tubocurarine", units = "mg",    specimen = "plasma",         verified = TRUE),
    # The effect state holds a CONCENTRATION, not an amount: Sheiner 1979
    # Eq. 7 rescales the hypothetical effect-compartment amount Ae to the
    # equivalent steady-state plasma concentration Cpss = keo * Z / V1, so
    # d/dt(effect) = ke0 * (Cc - effect) and `effect` is in ug/mL. The
    # compartment is explicitly hypothetical (p. 360: it receives negligible
    # mass and its exponential does not enter the PK solution), so no
    # biological matrix applies.
    effect      = list(analyte = "d-tubocurarine", units = "ug/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    RENALIMP_SEV = list(
      description        = "Chronic end-stage renal failure indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal serum creatinine)",
      notes              = paste(
        "1 = chronic end-stage renal failure undergoing renal transplantation",
        "(anaemic, haemoglobin 6.8 +/- 1.5 g/100 mL; serum creatinine",
        "10.1 +/- 2.3 mg/100 mL in spite of recent dialysis; studied before",
        "insertion of the transplant kidney), 0 = normal serum creatinine.",
        "Sheiner 1979 p. 362 'The patients', group 3. The paper's severity",
        "scheme is dialysis-dependent end-stage disease, not a creatinine",
        "clearance cut-off. In the baseline model the ONLY distinction",
        "between the two groups is the absence of a renal elimination rate",
        "constant in the renal failure patients (p. 363)."
      ),
      source_name        = "RF"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Not reported by Sheiner 1979. Every dose in the paper is given in",
        "mg/kg and the disposition parameters recovered from Fig. 3 are",
        "therefore per kilogram; WT converts them to the absolute amounts",
        "rxode2 doses in. Linear (exponent 1) scaling is the encoding that",
        "reproduces the paper's per-kg parameterisation exactly, so the",
        "predicted concentration for a mg/kg dose is weight invariant and",
        "matches Fig. 3 for any WT. The paper never fitted a weight",
        "covariate; see the vignette Assumptions and deviations section."
      ),
      source_name        = NA_character_
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20L,                                            # Sheiner 1979 p. 362: group 3 = 10 normal + 10 renal failure
    n_studies      = 1L,                                             # Sheiner 1979 p. 362; part of these data previously reported as Miller 1977 J Pharmacol Exp Ther 202:1-7 (reference 21), refitted here
    age_range      = "not reported for group 3",                     # Sheiner 1979 p. 362 gives ages only for group 1 (21-49 years)
    weight_range   = "not reported",                                 # Sheiner 1979 reports doses in mg/kg but no weights
    disease_state  = paste(
      "Adults undergoing elective surgery. 10 patients with normal serum",
      "creatinine; 10 patients with chronic end-stage renal failure",
      "undergoing renal transplantation, anaemic (haemoglobin 6.8 +/- 1.5",
      "g/100 mL) with serum creatinine 10.1 +/- 2.3 mg/100 mL in spite of",
      "recent dialysis, all studied prior to insertion of the transplant",
      "kidney."
    ),
    dose_range     = "Single IV bolus d-tubocurarine 0.3 mg/kg (5 per group) or 0.5 mg/kg (5 per group)",  # Sheiner 1979 p. 362
    anesthesia     = paste(
      "Induction with nitrous oxide:oxygen and halothane, trachea intubated",
      "without drugs; maintenance with nitrous oxide (60%):oxygen and",
      "halothane at end-tidal 0.45-0.8%. dTC given 20 min after induction.",
      "Temperature, ventilation and dTC effect monitored as in group 1."
    ),                                                               # Sheiner 1979 p. 362
    n_observations = "7 to 9 plasma and effect data points per patient (samples at 3, 15, 30, 45, 60, 90 and 120 min, plus 150 and 180 min after the 0.5 mg/kg dose)",  # Sheiner 1979 p. 362
    notes          = paste(
      "Effect is the force of thumb adduction after supramaximal ulnar",
      "nerve stimulation (single 0.1 ms stimuli at 0.3 pulse/s), expressed",
      "as degree of paralysis with 0 = no paralysis and 1.0 = 100%",
      "paralysis (Sheiner 1979 p. 362). dTC plasma concentrations by",
      "radioimmunoassay (Horowitz and Spector), assay CV 8%, lower limit of",
      "sensitivity 0.05 ug/mL (p. 362). Because only 7 to 9 data points per",
      "patient were available, plasma and effect data for all 20 patients",
      "were fitted simultaneously by a maximum likelihood nonlinear mixed",
      "effects method estimating population means and interindividual",
      "variability (p. 363) - the first PKPD nonlinear mixed effects",
      "analysis published."
    )
  )

  ini({
    # --- Disposition -------------------------------------------------------
    # Sheiner 1979 states explicitly (p. 365, Results): "The pharmacokinetic
    # parameter estimates are not presented because they are not of central
    # concern to this paper." The biexponential population fit IS published,
    # as the heavy solid "estimated mean population response" lines of Fig. 3
    # (p. 364), one per dose (0.3 and 0.5 mg/kg) and per group (normal, renal
    # failure). The values below were recovered by digitising all four of
    # those curves at the nine sampling times and fitting a single
    # two-compartment IV bolus model under the paper's own structural
    # constraint that renal failure removes only the renal arm of clearance.
    # The four curves are reproduced with a median absolute deviation of 1.2%
    # and a maximum of 4.9%; the two dose arms within a group are recovered as
    # proportional to within 2-3%, confirming the linear disposition the
    # biexponential model assumes. Values are for a 70 kg adult and scale
    # linearly with WT (see covariateData$WT). FIGURE-DERIVED - see the
    # vignette Assumptions and deviations section.
    lvc        <- log(5.443);   label("Central volume of distribution Vc (L, 70 kg adult)")                        # Sheiner 1979 Fig. 3 (digitised); 77.8 mL/kg
    lvp        <- log(11.87);   label("Peripheral volume of distribution Vp (L, 70 kg adult)")                     # Sheiner 1979 Fig. 3 (digitised); 169.5 mL/kg
    lq         <- log(0.3654);  label("Inter-compartmental clearance Q (L/min, 70 kg adult)")                      # Sheiner 1979 Fig. 3 (digitised); 5.22 mL/min/kg
    lcl_nonren <- log(0.1037);  label("Non-renal clearance (L/min, 70 kg adult)")                                  # Sheiner 1979 Fig. 3 (digitised); 1.48 mL/min/kg - the only clearance retained in renal failure
    lcl_renal  <- log(0.06639); label("Renal clearance arm, absent in renal failure (L/min, 70 kg adult)")         # Sheiner 1979 Fig. 3 (digitised); 0.95 mL/min/kg = 39% of total clearance in normals

    # Linear weight scaling. Fixed structural constants, not estimates: the
    # paper dosed exclusively in mg/kg and reported no weights, so exponent 1
    # is the encoding that reproduces its per-kg parameterisation.
    e_wt_cl <- fixed(1); label("Weight exponent on clearance terms (unitless)")  # Sheiner 1979 p. 362 (all dosing in mg/kg); see covariateData$WT
    e_wt_vc  <- fixed(1); label("Weight exponent on volume terms (vc and vp) (unitless)")     # Sheiner 1979 p. 362 (all dosing in mg/kg); see covariateData$WT

    # --- Pharmacodynamics --------------------------------------------------
    # Table II, group 3 (N = 20) column - the simultaneous PK/PD nonlinear
    # mixed effects fit that also produced the Fig. 3 concentration curves, so
    # this is the parameter set that is self-consistent with the disposition
    # above. Table II reports mean +/- SEM. The group 1 column (N = 7,
    # individual nonlinear least-squares fits of the dual-infusion study) gives
    # ke0 0.13 +/- 0.015 /min, gamma 2.53 +/- 0.16 and Cpss(50) 0.37 +/- 0.02
    # ug/mL; the paper reports no statistically significant difference between
    # the two groups' dynamic parameters (p. 366).
    lke0   <- log(0.16);  label("Effect-compartment equilibration rate constant keo (1/min)")                       # Sheiner 1979 Table II group 3: keo = 0.16 +/- 0.14 /min; tabulated t1/2 = 4.33 min = log(2)/0.16
    lhill  <- log(2.30);  label("Hill / sigmoidicity coefficient gamma (unitless)")                                 # Sheiner 1979 Table II group 3: gamma = 2.30 +/- 0.01
    lec50  <- log(0.38);  label("Steady-state plasma concentration producing 50% paralysis, Cpss(50) (ug/mL)")      # Sheiner 1979 Table II group 3: Cpss(50) = 0.38 +/- 0.03 ug/mL
    lemax  <- fixed(log(1)); label("Maximum effect Emax (fraction of complete paralysis)")                          # Sheiner 1979 Eq. 1 and Eq. 8: E is the fraction of maximal response, 0 = no paralysis to 1.0 = 100% paralysis (p. 362)

    # No inter-individual variability or residual error terms are encoded.
    # The group 3 analysis did estimate interindividual variability
    # (p. 363: the method "directly estimates mean population values of the
    # pharmacokinetic and pharmacodynamic parameters along with their
    # interindividual variability"), but no omega or sigma value is tabulated
    # anywhere in the paper, and none may be invented. The model is therefore
    # the typical-value (population mean) prediction, which is exactly what
    # the heavy solid lines of Fig. 3 and Fig. 4 depict.
  })

  model({
    # Individual parameters. Renal failure removes the renal arm of clearance
    # and changes nothing else - Sheiner 1979 p. 363: "the baseline model
    # chosen was one in which the only distinction between N and RF patients
    # was the absence of a renal elimination rate constant in the RF patients."
    wt_cl <- (WT / 70)^e_wt_cl
    wt_v  <- (WT / 70)^e_wt_vc

    cl <- (exp(lcl_nonren) + exp(lcl_renal) * (1 - RENALIMP_SEV)) * wt_cl
    q  <- exp(lq)  * wt_cl
    vc <- exp(lvc) * wt_v
    vp <- exp(lvp) * wt_v

    ke0  <- exp(lke0)
    hill <- exp(lhill)
    ec50 <- exp(lec50)
    emax <- exp(lemax)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment mammillary disposition after IV bolus into `central`
    # (Sheiner 1979 Eq. 3, "a biexponential equation interpreted as a
    # 2-compartment mammillary model", p. 362).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    Cc <- central / vc                                              # ug/mL (= mg/L)

    # Hypothetical effect compartment, Sheiner 1979 Eq. 2 rescaled by Eq. 7.
    # Eq. 2 is dAe/dt = k1e * A1 - keo * Ae; Eq. 7 defines the equivalent
    # steady-state plasma concentration Cpss = keo * Ae / (k1e * V1), in which
    # k1e cancels ("Note that k1e cancels out of Equation 7 so that ... its
    # exact value is of no importance", p. 361). Substituting gives the form
    # below directly, so `effect` is Cpss in ug/mL and equals Cc at steady
    # state.
    d/dt(effect) <- ke0 * (Cc - effect)

    # Sigmoid Emax (Hill) relation between the effect-site concentration and
    # the degree of paralysis (Sheiner 1979 Eq. 1 and Eq. 8). With emax fixed
    # to 1, paralysis runs from 0 (no paralysis) to 1 (100% paralysis).
    paralysis <- emax * effect^hill / (ec50^hill + effect^hill)
  })
}
