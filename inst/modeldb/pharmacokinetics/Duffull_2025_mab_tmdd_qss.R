Duffull_2025_mab_tmdd_qss <- function() {
  description <- "Two-compartment target-mediated drug disposition (TMDD) model with quasi-steady-state target binding confined to the central compartment, for an unnamed monoclonal antibody (mAb) and its soluble target; case example 1 of the Duffull 2025 model-instability tutorial (Equations 1-3; Table 2 'Nominal value' column). Free antibody exchanges with a peripheral compartment by first-order rate constants and is removed both by linear clearance and by saturable target-mediated internalisation of the antibody-target complex. The total target (free target plus antibody-target complex) is a dynamic state with zero-order synthesis and first-order degradation, initialised at its drug-free steady state ksyn*Vc/kdeg. Cc is the free antibody concentration and Ctotal_target is the total target concentration, both nmol/L. This is the FULL model that the tutorial shows to be structurally identifiable under a rich design but NOT deterministically identifiable under the reduced clinical design (Km relative standard error 1495.65%); the target-saturated simplification that the authors adopted instead is packaged separately as Duffull_2025_mab_tmdd_simplified. Parameter values are the nominal set used for the Fisher-information evaluations, not a fit to observed patient data: the tutorial's datasets were generated under the stated sampling designs."
  reference   <- "Duffull SB, Wright DFB, Zhu X, Liu X, Abulfathi A, Hishe H. A pharmacometric workflow for resolving model instability in model use-reuse settings. CPT Pharmacometrics Syst Pharmacol. 2025;14(10):1547-1556. doi:10.1002/psp4.70049"
  vignette    <- "Duffull_2025_tmdd_model_instability"
  units       <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte and units are SOURCE-CONFIRMED from the
  # Notes paragraph under Duffull 2025 Equations 1-3 ("LC represents the
  # ligand (drug) amount in central compartment; LT represents the
  # ligand (drug) amount in peripheral compartment; Rtot represents the
  # total amount of the receptor and receptor-drug complex") and from
  # Figure 3, whose axes are labelled in nmol. specimen is INFERRED: the
  # source never names the sampling matrix (it reports only "free drug
  # concentrations" and "total target concentrations"). verified = FALSE
  # therefore records the unverified specimen, not an unverified analyte.
  compartmentData <- list(
    central      = list(analyte = "monoclonal antibody (unnamed), free",       units = "nmol", specimen = "plasma", verified = FALSE),
    peripheral1  = list(analyte = "monoclonal antibody (unnamed), free",       units = "nmol", specimen = "plasma", verified = FALSE),
    total_target = list(analyte = "target receptor, free plus antibody-bound", units = "nmol", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  # Duffull 2025 turns covariate effects OFF for the identifiability
  # analyses ("Covariate effects will need to be turned off as the
  # covariate values will not be in your simulated data file as these are
  # seldom the cause of a structural identifiability issue", Section
  # 1.2.1) and reports no covariate model for this case example.
  covariatesDataExcluded <- list()

  population <- list(
    species        = "human (implied, not stated: the source frames the case example as a monoclonal antibody with sparse sampling 'common in mAb clinical trials')",
    n_subjects     = 80L,
    n_studies      = 1L,
    disease_state  = "Not reported. The source is a pharmacometric tutorial on resolving model instability; neither the antibody nor its target is named, and no demographics are given.",
    dose_range     = "Single intravenous bolus into the central compartment (Duffull 2025 Section 2.1.5). The dose amount is not stated in the text; Figure 3 shows the free-drug amount in the central compartment starting at 10000 nmol, so the simulated dose is 10 umol.",
    samples        = "Complete sampling design: 640 free drug and 720 total target concentrations from 80 subjects. Reduced sampling design: 400 free drug and 480 total target concentrations from the same 80 subjects (Duffull 2025 Section 2.1.1).",
    designs        = "Hypothetical rich design used for the structural-identifiability template data set: 0, 0.01, 0.2, 0.6, 3, 50, 64, 81, 91 h post-dose (9 times). Reduced clinical design used for the deterministic-identifiability analysis: 0, 0.5, 3, 5, 12, 24 h post-dose (6 times) (Duffull 2025 Section 2.1.3).",
    notes          = "The tutorial's datasets are SIMULATED, not observed: the reported observation counts are exactly the sampling grids times 80 subjects (80 x 9 = 720 total target, 80 x 6 = 480 total target), and the free-drug counts are one time point fewer in each case (80 x 8 = 640, 80 x 5 = 400). The model structure was 'adapted from a previously published TMDD example in the tutorial for $DESIGN in NONMEM' (Duffull 2025 Section 2.1.1, reference [24] = Bauer 2019), which is not on disk; the structure and every parameter value used here come from Duffull 2025 itself (Equations 1-3 and Table 2), so no value is inherited from the unavailable upstream tutorial."
  )

  ini({
    # ----- Antibody disposition (Duffull 2025 Table 2, 'Nominal value' column). -----
    # These are the nominal values at which the Fisher information matrix
    # was evaluated for the structural- and deterministic-identifiability
    # analyses (Section 2.1.3), i.e. the parameter set for the full TMDD
    # model that "performed well with a complete sampling dataset"
    # (Section 2.1.1). They are not accompanied by standard errors of
    # their own; Table 2's RSE% columns are the FIM-predicted precisions
    # of these values under the two designs.
    lcl  <- log(4.13);  label("Linear (non-specific) clearance of free antibody from the central compartment (CL, L/h)")  # Duffull 2025 Table 2 row CL, nominal value 4.13
    lvc  <- log(76.7);  label("Central compartment volume of distribution (Vc, L)")                                       # Duffull 2025 Table 2 row Vc, nominal value 76.7
    lvp  <- log(70);    label("Peripheral compartment volume of distribution (Vp, L)")                                    # Duffull 2025 Table 2 row Vp, nominal value 70
    lq   <- log(44.7);  label("Intercompartmental clearance between central and peripheral compartments (Q, L/h)")        # Duffull 2025 Table 2 row Q, nominal value 44.7

    # ----- Target turnover (Duffull 2025 Table 2, 'Nominal value' column). -----
    # ksyn is reported per unit volume (nM/h); model() multiplies it by
    # Vc to obtain the zero-order synthesis rate of target AMOUNT
    # (nmol/h), matching Equation 3's "+ ksyn Vc" term.
    lksyn <- log(1.09);  label("Zero-order synthesis rate of the target, per unit central volume (ksyn, nmol/L/h)")  # Duffull 2025 Table 2 row ksyn, nominal value 1.09 nM/h
    lkdeg <- log(0.349); label("First-order degradation rate constant of free target (kdeg, 1/h)")                   # Duffull 2025 Table 2 row kdeg, nominal value 0.349
    lkint <- log(11.1);  label("First-order internalisation (elimination) rate constant of the antibody-target complex (kint, 1/h)")  # Duffull 2025 Table 2 row kint, nominal value 11.1

    # ----- Quasi-steady-state binding constant. -----
    # Duffull 2025 calls this parameter "Km" and the Table 2 footnote
    # glosses it as the "Michaelis-Menten constant". Its ROLE in
    # Equations 1 and 3 is the quasi-steady-state (quasi-equilibrium)
    # binding constant of the Gibiansky 2008 QSS TMDD reduction -- it
    # scales the free-drug saturation of target-mediated internalisation
    # while the total target remains a dynamic state -- so it is
    # recorded under the canonical `lkss`, the same name the library's
    # PK_2cmt_tmdd_qss.R archetype uses for the same quantity. It is NOT
    # a metabolic Michaelis-Menten constant (canonical `lkm`), which
    # would describe saturable elimination with no target state.
    # This is the parameter the tutorial identifies as the source of the
    # instability: RSE 38.12% under the rich design but 1495.65% under
    # the reduced clinical design (Table 2).
    lkss <- log(0.13);  label("Quasi-steady-state antibody-target binding constant (Km in the source; Kss, nmol/L)")  # Duffull 2025 Table 2 row Km, nominal value 0.13 nM

    # ----- IIV (Duffull 2025 Table 2, 'Nominal value' column). -----
    # Table 2 reports between-subject variability as "IIV-<parameter>
    # (CV%)" with a nominal value of 25 for each of CL, Vc, ksyn and
    # kint. These are round assumed values supplied to $DESIGN, whose
    # $OMEGA inputs are variances. They are encoded here under the
    # NONMEM reporting convention
    #     omega^2 = (CV% / 100)^2
    # i.e. the reported %CV is taken as the log-scale standard deviation
    # directly. This is the same convention adopted for this author
    # group elsewhere in the library (see Wright_2016_allopurinol.R,
    # whose second author is this paper's second author). The
    # alternative lognormal reading omega^2 = log(1 + CV^2) would give
    # 0.06062 instead of 0.0625 -- a 3% difference in variance that the
    # source provides no reported omega^2 with which to discriminate.
    # No IIV is reported on Vp, Q, kdeg or Km. Section 2.1.2 records
    # that the between-subject variance of Km had to be FIXED: left
    # free it "became inflated, with an RSE of 10,367% and a shrinkage
    # of 100%".
    etalcl   ~ 0.0625  # Duffull 2025 Table 2 row IIV-CL,   nominal 25 CV% -> 0.25^2
    etalvc   ~ 0.0625  # Duffull 2025 Table 2 row IIV-Vc,   nominal 25 CV% -> 0.25^2
    etalksyn ~ 0.0625  # Duffull 2025 Table 2 row IIV-ksyn, nominal 25 CV% -> 0.25^2
    etalkint ~ 0.0625  # Duffull 2025 Table 2 row IIV-kint, nominal 25 CV% -> 0.25^2

    # ----- Residual error: NOT REPORTED by the source. -----
    # Duffull 2025 Table 2 carries no sigma row for either output, and
    # the paper reports no residual-error model for this case example.
    # The only sigma values it mentions are generic workflow advice for
    # setting up a $DESIGN run -- "fix $SIGMA initial values to be
    # small, for example, 1E-3" (Section 1.2.1) and "Often, a more
    # appropriate value might be 0.0225 (CV% = 15%)" (Section 1.2.2) --
    # which are illustrative defaults for the method, NOT estimates for
    # this model. They are therefore not used here. Both residuals are
    # fixed at zero so that no variance is invented; the packaged model
    # is a typical-value / IIV-only simulation model. See the vignette's
    # "Assumptions and deviations" section.
    propSd               <- fixed(0); label("Proportional residual error on free antibody Cc (fraction; ZERO - not reported in source)")               # Duffull 2025: no residual-error model reported for case example 1
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
    kss  <- exp(lkss)

    # ----- Micro-constants -----
    # Duffull 2025 parameterises Equations 1-2 in rate constants (ke,
    # k12, k21) but tabulates CL, Vc, Vp and Q, so the rate constants
    # are derived: ke = CL/Vc, k12 = Q/Vc, k21 = Q/Vp.
    ke  <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system (Duffull 2025 Equations 1-3) -----
    # All three states are AMOUNTS in nmol. kss is a concentration
    # (nmol/L), so kss * vc converts it to the amount scale on which
    # `central` is expressed -- this is the source's "Km Vc" term.
    #
    # The saturation fraction central/(central + kss*vc) is written
    # INLINE in both ODEs rather than via a named intermediate: in the
    # nlmixr2 model-function form a named intermediate that references
    # an ODE state can silently evaluate to zero inside d/dt(), which
    # would delete the entire target-mediated term with no error (see
    # the exaggeration test in the validation vignette, which gates
    # against exactly that failure).
    #
    # Equation 1: dLC/dt = -(ke + k12) LC + k21 LT - kint LC Rtot / (LC + Km Vc)
    d/dt(central) <- -(ke + k12) * central + k21 * peripheral1 -
      kint * total_target * (central / (central + kss * vc))

    # Equation 2: dLT/dt = k12 LC - k21 LT
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Equation 3: dRtot/dt = -(kint - kdeg) LC Rtot / (LC + Km Vc) - kdeg Rtot + ksyn Vc
    d/dt(total_target) <- -(kint - kdeg) * total_target * (central / (central + kss * vc)) -
      kdeg * total_target + ksyn * vc

    # Initial conditions (Duffull 2025, stated with Equations 4-5):
    #   LC(t = 0) = 0, LT(t = 0) = 0, Rtot(t = 0) = ksyn Vc / kdeg
    # The central and peripheral states start empty (the bolus arrives
    # from the event table); the total target starts at its drug-free
    # steady state, where no complex exists and the free target is
    # synthesised at ksyn*Vc and degraded at kdeg.
    total_target(0) <- ksyn * vc / kdeg

    # ----- Observations -----
    # Both source outputs are concentrations: "free drug concentrations"
    # and "total target concentrations (sum of free target and
    # drug-target complex)" (Duffull 2025 Section 2.1.1). Dividing the
    # nmol amounts by Vc (L) gives nmol/L.
    Cc             <- central      / vc
    Ctotal_target  <- total_target / vc

    Cc            ~ prop(propSd)
    Ctotal_target ~ prop(propSd_Ctotal_target)
  })
}
