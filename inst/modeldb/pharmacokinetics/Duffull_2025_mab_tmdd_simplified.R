Duffull_2025_mab_tmdd_simplified <- function() {
  description <- "Two-compartment TMDD model for an unnamed monoclonal antibody (mAb) simplified under the target-saturated limit, and the model the Duffull 2025 tutorial adopted to resolve the instability of its full counterpart (Equations 2, 4 and 5; Table 2 'Simplified model' column). Because the free-antibody amount in the central compartment exceeds Km*Vc by roughly two orders of magnitude over the 24 h observation window, the binding saturation fraction LC/(LC + Km*Vc) is set to 1: the binding constant Km drops out of the model entirely and target-mediated removal of antibody becomes kint*Rtot, independent of antibody amount. The total target still has zero-order synthesis and first-order internalisation, and kdeg survives only through the drug-free initial condition ksyn*Vc/kdeg. Cc is the free antibody concentration and Ctotal_target is the total target concentration, both nmol/L. All parameters were estimated with relative standard errors below 30% under the reduced clinical sampling design that made the full model unidentifiable. IMPORTANT: this approximation is valid only while the target remains saturated -- it makes target-mediated antibody removal a constant-rate sink, so the free-antibody state can be driven negative once antibody washes out. Restrict simulations to the 24 h window the source validated. The full model is packaged as Duffull_2025_mab_tmdd_qss."
  reference   <- "Duffull SB, Wright DFB, Zhu X, Liu X, Abulfathi A, Hishe H. A pharmacometric workflow for resolving model instability in model use-reuse settings. CPT Pharmacometrics Syst Pharmacol. 2025;14(10):1547-1556. doi:10.1002/psp4.70049"
  vignette    <- "Duffull_2025_tmdd_model_instability"
  units       <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # Issue #482: analyte and units are SOURCE-CONFIRMED from the Notes
  # paragraph under Duffull 2025 Equations 1-3 and from Figure 3, whose
  # axes are labelled in nmol. specimen is INFERRED: the source never
  # names the sampling matrix. verified = FALSE therefore records the
  # unverified specimen, not an unverified analyte.
  compartmentData <- list(
    central      = list(analyte = "monoclonal antibody (unnamed), free",       units = "nmol", specimen = "plasma", verified = FALSE),
    peripheral1  = list(analyte = "monoclonal antibody (unnamed), free",       units = "nmol", specimen = "plasma", verified = FALSE),
    total_target = list(analyte = "target receptor, free plus antibody-bound", units = "nmol", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list()

  population <- list(
    species        = "human (implied, not stated: the source frames the case example as a monoclonal antibody with sparse sampling 'common in mAb clinical trials')",
    n_subjects     = 80L,
    n_studies      = 1L,
    disease_state  = "Not reported. The source is a pharmacometric tutorial on resolving model instability; neither the antibody nor its target is named, and no demographics are given.",
    dose_range     = "Single intravenous bolus into the central compartment (Duffull 2025 Section 2.1.5). The dose amount is not stated in the text; Figure 3 shows the free-drug amount in the central compartment starting at 10000 nmol, so the simulated dose is 10 umol.",
    samples        = "This simplified model was estimated on the REDUCED sampling design: 400 free drug and 480 total target concentrations from 80 subjects (Duffull 2025 Section 2.1.1). The reduced design is the one under which the full model's Km was not deterministically identifiable.",
    designs        = "Reduced clinical design: 0, 0.5, 3, 5, 12, 24 h post-dose (6 times) (Duffull 2025 Section 2.1.3).",
    notes          = "The tutorial's datasets are SIMULATED, not observed: the reported observation counts are exactly the sampling grids times 80 subjects (80 x 6 = 480 total target, 80 x 5 = 400 free drug). Unlike the full model's nominal column, the values in Table 2's 'Simplified model' column ARE estimates: Section 2.1.5 reports that this model achieved 'full convergence, with both the estimation minimization (S) and covariance (C) steps completed' and that 'all parameters were estimated precisely with RSE for all parameters being well below 30%'."
  )

  ini({
    # ----- Antibody disposition (Duffull 2025 Table 2, 'Simplified model'
    # ----- column: final parameter estimates with their RSE%). -----
    lcl  <- log(4.15);  label("Linear (non-specific) clearance of free antibody from the central compartment (CL, L/h)")  # Duffull 2025 Table 2 row CL, simplified model 4.15 (RSE 6%)
    lvc  <- log(68.4);  label("Central compartment volume of distribution (Vc, L)")                                       # Duffull 2025 Table 2 row Vc, simplified model 68.4 (RSE 9%)
    lvp  <- log(75.8);  label("Peripheral compartment volume of distribution (Vp, L)")                                    # Duffull 2025 Table 2 row Vp, simplified model 75.8 (RSE 5%)
    lq   <- log(59.9);  label("Intercompartmental clearance between central and peripheral compartments (Q, L/h)")        # Duffull 2025 Table 2 row Q, simplified model 59.9 (RSE 11%)

    # ----- Target turnover (Duffull 2025 Table 2, 'Simplified model'
    # ----- column). -----
    # ksyn is reported per unit volume (nM/h); model() multiplies it by
    # Vc to obtain the zero-order synthesis rate of target AMOUNT
    # (nmol/h), matching Equation 5's "+ ksyn VC" term.
    lksyn <- log(1.16);  label("Zero-order synthesis rate of the target, per unit central volume (ksyn, nmol/L/h)")  # Duffull 2025 Table 2 row ksyn, simplified model 1.16 (RSE 6%)
    lkint <- log(11.5);  label("First-order internalisation (elimination) rate constant of the antibody-target complex (kint, 1/h)")  # Duffull 2025 Table 2 row kint, simplified model 11.5 (RSE 5%)

    # kdeg does NOT appear in Equations 4 or 5. It survives in this model
    # only through the drug-free initial condition Rtot(0) = ksyn*Vc/kdeg
    # that Duffull 2025 states alongside Equations 4-5, so it is
    # identified purely by the baseline total target concentration
    # (ksyn/kdeg = 1.16/0.373 = 3.11 nmol/L). Table 2 nonetheless reports
    # it as estimated in the simplified model with an RSE of 4%.
    lkdeg <- log(0.373); label("First-order degradation rate constant of free target; enters only via the baseline Rtot(0) = ksyn*Vc/kdeg (kdeg, 1/h)")  # Duffull 2025 Table 2 row kdeg, simplified model 0.373 (RSE 4%)

    # No binding constant: Km is absent from the simplified model. Table 2
    # shows an em-dash for both its estimate and its RSE in the
    # 'Simplified model' column, because setting the saturation fraction
    # LC/(LC + Km*Vc) to 1 removes Km from Equations 4 and 5 entirely.
    # This is the whole point of the simplification -- Km was the
    # parameter with RSE 1495.65% under this design.

    # ----- IIV (Duffull 2025 Table 2, 'Simplified model' column). -----
    # Reported as "IIV-<parameter> (CV%)" and encoded under the NONMEM
    # reporting convention omega^2 = (CV% / 100)^2, i.e. the reported
    # %CV taken as the log-scale standard deviation directly. This is
    # the same convention adopted for this author group elsewhere in the
    # library (Wright_2016_allopurinol.R). No IIV is reported on Vp, Q
    # or kdeg.
    etalcl   ~ 0.047961  # Duffull 2025 Table 2 row IIV-CL,   simplified model 21.9 CV% (RSE 15%) -> 0.219^2
    etalvc   ~ 0.096721  # Duffull 2025 Table 2 row IIV-Vc,   simplified model 31.1 CV% (RSE 17%) -> 0.311^2
    etalksyn ~ 0.106276  # Duffull 2025 Table 2 row IIV-ksyn, simplified model 32.6 CV% (RSE 14%) -> 0.326^2
    etalkint ~ 0.082369  # Duffull 2025 Table 2 row IIV-kint, simplified model 28.7 CV% (RSE 13%) -> 0.287^2

    # ----- Residual error: NOT REPORTED by the source. -----
    # Duffull 2025 Table 2 carries no sigma row for either output even in
    # the 'Simplified model' column, and no residual-error model is
    # reported for case example 1. The generic sigma values the paper
    # mentions (1E-3 in Section 1.2.1, 0.0225 in Section 1.2.2) are
    # illustrative $DESIGN set-up advice for the method, not estimates
    # for this model, and are not used. Both residuals are fixed at zero
    # so that no variance is invented. See the vignette's "Assumptions
    # and deviations" section.
    propSd               <- fixed(0); label("Proportional residual error on free antibody Cc (fraction; ZERO - not reported in source)")             # Duffull 2025: no residual-error model reported for case example 1
    propSd_Ctotal_target <- fixed(0); label("Proportional residual error on total target Ctotal_target (fraction; ZERO - not reported in source)")  # Duffull 2025: no residual-error model reported for case example 1
  })

  model({
    # ----- Individual parameters -----
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    vp   <- exp(lvp)
    q    <- exp(lq)
    ksyn <- exp(lksyn + etalksyn)
    kdeg <- exp(lkdeg)
    kint <- exp(lkint + etalkint)

    # ----- Micro-constants -----
    # Duffull 2025 parameterises Equations 2 and 4 in rate constants (ke,
    # k12, k21) but tabulates CL, Vc, Vp and Q, so the rate constants
    # are derived: ke = CL/Vc, k12 = Q/Vc, k21 = Q/Vp.
    ke  <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system (Duffull 2025 Equations 4, 2 and 5) -----
    # All three states are AMOUNTS in nmol.
    #
    # Equation 4: dLC/dt = -(ke + k12) LC + k21 LT - kint Rtot
    # The target-mediated loss term is the full model's
    #   kint * Rtot * LC / (LC + Km*Vc)
    # with the saturation fraction set to 1. It is therefore INDEPENDENT
    # of the antibody amount -- a constant-rate sink whose size is set by
    # the target only. That is what makes this form valid only while the
    # target is saturated (see the description and the vignette).
    d/dt(central) <- -(ke + k12) * central + k21 * peripheral1 - kint * total_target

    # Equation 2 (unchanged from the full model): dLT/dt = k12 LC - k21 LT
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Equation 5: dRtot/dt = -kint Rtot + ksyn VC
    # Applying the same saturation limit to Equation 3 cancels kdeg
    # exactly:
    #   -(kint - kdeg) Rtot - kdeg Rtot + ksyn Vc = -kint Rtot + ksyn Vc
    d/dt(total_target) <- -kint * total_target + ksyn * vc

    # Initial conditions (Duffull 2025, stated with Equations 4-5):
    #   LC(t = 0) = 0, LT(t = 0) = 0, Rtot(t = 0) = ksyn VC / kdeg
    # Note that this is the DRUG-FREE steady state of the FULL model
    # (governed by kdeg), which is not the steady state of Equation 5
    # (ksyn*Vc/kint). The total target therefore drops rapidly from
    # ksyn*Vc/kdeg toward ksyn*Vc/kint at the onset of dosing, as it
    # should: the target is saturated with antibody from the first
    # instant and is cleared by internalisation rather than by
    # degradation.
    total_target(0) <- ksyn * vc / kdeg

    # ----- Observations -----
    # Both source outputs are concentrations: "free drug concentrations"
    # and "total target concentrations (sum of free target and
    # drug-target complex)" (Duffull 2025 Section 2.1.1). Dividing the
    # nmol amounts by Vc (L) gives nmol/L.
    Cc            <- central      / vc
    Ctotal_target <- total_target / vc

    Cc            ~ prop(propSd)
    Ctotal_target ~ prop(propSd_Ctotal_target)
  })
}
