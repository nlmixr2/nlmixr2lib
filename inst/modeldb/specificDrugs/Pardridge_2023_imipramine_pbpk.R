Pardridge_2023_imipramine_pbpk <- function() {
  description <- paste(
    "PBPK (brain, semi-mechanistic). Pardridge 2023 Pharmaceutical",
    "Research partly-flow / partly-compartmental model of imipramine",
    "delivery to brain that resolves plasma-protein binding as an",
    "explicit kinetic process inside the brain capillary. Seven brain",
    "states carry albumin-bound, alpha-1-acid-glycoprotein (AGP)-bound",
    "and free drug in the brain capillary, free AGP in the capillary,",
    "free drug in brain, and free / drug-occupied brain cytoplasmic",
    "binding protein. Drug crosses the blood-brain barrier by",
    "bi-directional first-order permeation and may be metabolised in",
    "brain. Imipramine shows a high degree of plasma-protein-mediated",
    "uptake from BOTH the AGP-bound and the albumin-bound pools",
    "(KG in vivo 90 uM vs KG in vitro 1.2 uM; KA in vivo >1,000 uM vs",
    "KA in vitro 42 uM), so free imipramine in brain is 18- to 31-fold",
    "higher than equilibrium dialysis of plasma predicts. Arterial input",
    "is a constant total plasma concentration (steady-state IV-infusion",
    "model, c_plasma_ss, Table II LT0 = 100 nM); the paper did not run a",
    "non-steady-state oral model for imipramine because detailed oral",
    "pharmacokinetic parameters were not available, so this model has no",
    "dosing compartment. All parameters are FIXED at the paper's",
    "literature-sourced values; the paper reports no IIV and no residual",
    "error because it is a deterministic simulation study rather than a",
    "fit to individual data."
  )
  reference <- paste(
    "Pardridge WM (2023). Physiologically Based Pharmacokinetic Model of",
    "Brain Delivery of Plasma Protein Bound Drugs. Pharmaceutical",
    "Research 40(3):661-674. doi:10.1007/s11095-023-03484-2.",
    sep = " "
  )
  vignette <- "Pardridge_2023_brain_plasma_protein_binding"
  units    <- list(time = "minute", dosing = "nmol/kg", concentration = "nM")

  # Paper-mechanistic brain states (Table I).  The `brain_vascular` /
  # `brain_extravascular` stems are the registered canonical brain
  # compartments; the trailing role tokens distinguish the bound and
  # free drug species and the binding proteins themselves, which have
  # no canonical equivalent.
  paper_specific_compartments <- c(
    "brain_vascular_drug_free",       # LF  - free (bioavailable) drug in brain capillary
    "brain_vascular_drug_alb",        # AL  - albumin-bound drug in brain capillary
    "brain_vascular_drug_agp",        # GL  - AGP(globulin)-bound drug in brain capillary
    "brain_vascular_agp_free",        # GF  - free (unbound) AGP in brain capillary
    "brain_extravascular_drug_free",  # LM  - free drug in brain
    "brain_extravascular_drug_prot",  # PL  - drug bound to brain cytoplasmic protein
    "brain_extravascular_prot_free"   # PF  - free brain cytoplasmic binding protein
  )

  # No covariates: every parameter is a fixed literature constant and the
  # paper simulates a single typical subject rather than a population.
  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    disease_state  = paste(
      "Not a fit to individual data. Deterministic simulation of a typical",
      "human parameterised from published human plasma-protein chemistry.",
      "The total plasma imipramine concentration of 100 nM used by the",
      "steady-state model is a pharmacologic plasma level in treated",
      "subjects (Table II, ref 11)."
    ),
    dose_range     = paste(
      "No dose events. The steady-state model holds the total plasma",
      "imipramine concentration constant at 100 nM, as during an IV",
      "infusion. The paper states that the non-steady-state model was",
      "examined only for propranolol, because detailed pharmacokinetic",
      "parameters are not available following oral administration of",
      "imipramine (Methods)."
    ),
    notes          = paste(
      "Species provenance is mixed and should be understood before reuse.",
      "The plasma-protein concentrations (albumin 800 uM, AGP 20 uM; and",
      "the metastatic-cancer values albumin 600 uM, AGP 70 uM) and the",
      "in vitro dissociation constants are human. The in vivo",
      "brain-capillary dissociation constants and the BBB permeation rate",
      "constants k3 / k4 derive from the rat blood-brain-barrier study of",
      "Riant et al. 1988 (ref 9). The observed Kp,brain values the paper",
      "compares against are mouse (23, ref 40) and rat (26-30, refs",
      "41-42). The author argues (Discussion) that rodent albumin and AGP",
      "differ substantially from the human proteins, so cross-species",
      "transfer of the binding constants is explicitly a limitation of",
      "the model rather than a validated assumption."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # All parameters are FIXED. Pardridge 2023 is a simulation study:
    # every value is taken from the cited literature or from the
    # author's own simulation-based calibration, and nothing is
    # estimated from data. Values are Table II, "Imipramine" column,
    # cross-checked against Supplementary Table S2 simulation 1
    # (= Table IV simulation 11).
    #
# Positive-constrained mechanistic constants are log-transformed
    # (`l` prefix) and back-transformed at the top of model(), which is
    # also what rxode2's mu-referencing parser requires: at most one
    # population parameter per expression. kmet_brain (k9) is
    # legitimately ZERO in the paper's basal simulation and stays on the
    # natural scale, because log(0) is undefined.
    # ---------------------------------------------------------------

    # ---- Drug binding to AGP (alpha-1-acid glycoprotein) in the brain
    # capillary. KG in vivo = k1/k2 = 90 uM vastly exceeds KG in vitro =
    # 1.2 uM: strongly enhanced dissociation, hence high PMU from AGP.
    lkoff_agp <- fixed(log(5400));      label("Drug dissociation rate from AGP in brain capillary, k1 (1/min)")             # Table II (ref 9); k1 = k2 * KG,in vivo = 0.06 * 90,000
    lkon_agp  <- fixed(log(0.06));      label("Drug association rate with AGP in brain capillary, k2 (1/(nM*min))")         # Table II (ref 16); 0.06 /nM/min = 1e6 /M/s

    # ---- Drug binding to albumin in the brain capillary. KA in vivo =
    # k7/k8 = 1,000 uM vs KA in vitro = 42 uM: unlike propranolol,
    # imipramine ALSO shows enhanced dissociation from albumin. Table II
    # reports k7 as ">= 6,000 /min"; the simulations use 6,000 exactly
    # (Supplementary Table S2), which is the value encoded here.
    lkoff_alb <- fixed(log(6000));      label("Drug dissociation rate from albumin in brain capillary, k7 (1/min)")         # Table II (ref 9), reported as ">= 6,000"; k7 = k8 * KA,in vivo = 0.006 * 1,000,000
    lkon_alb  <- fixed(log(0.006));     label("Drug association rate with albumin in brain capillary, k8 (1/(nM*min))")     # Table II (ref 19); 0.006 /nM/min = 1e5 /M/s

    # ---- Bi-directional blood-brain barrier permeation. k3 is derived
    # from the PS/F ratio of 2.5 with k10 = 60 /min; k4 = k3 * (VP/VT)
    # assumes symmetric transport (2.1 = 150 * 0.01 / 0.7).
    lkinflux_bbb <- fixed(log(150));    label("Drug influx rate constant plasma to brain across the BBB, k3 (1/min)")       # Table II (ref 9); k3 = k10 * PS/F = 60 * 2.5
    lkefflux_bbb <- fixed(log(2.1));    label("Drug efflux rate constant brain to plasma across the BBB, k4 (1/min)")       # Table II (ref 9); k4 = k3 * VP/VT = 150 * 0.01/0.7

    # ---- Drug binding to brain cytoplasmic protein. KP = k6/k5 = 83 nM.
    # Table II reports k6 as "<= 0.5 /min" because imipramine
    # sequestration in brain exceeds that of propranolol; the
    # simulations use 0.5 exactly (Supplementary Table S2).
    lkon_brainprot  <- fixed(log(0.006)); label("Drug association rate with brain binding protein, k5 (1/(nM*min))")        # Table II; 0.006 /nM/min = 1e5 /M/s, same value calibrated for propranolol
    lkoff_brainprot <- fixed(log(0.5));   label("Drug dissociation rate from brain binding protein, k6 (1/min)")            # Table II (ref 9), reported as "<= 0.5"

    # ---- Brain metabolism and brain capillary plasma flow.
    kmet_brain   <- fixed(0);           label("Drug metabolism rate constant in brain, k9 (1/min)")                         # Table II (ref 20): negligible brain metabolism of imipramine
    lkflow_brain <- fixed(log(60));     label("Brain capillary plasma flow rate constant, k10 (1/min)")                     # Table II (ref 15); 1-second capillary transit time

    # ---- Binding-protein concentrations.
    lc_alb       <- fixed(log(800000)); label("Total plasma / capillary albumin concentration, AF (nM)")                    # Table II (ref 12); 5.4 g/100 mL of 67 kDa albumin ~ 800 uM
    lc_agp       <- fixed(log(20000));  label("Total plasma / capillary AGP concentration, GT0 (nM)")                       # Table II (ref 12); 0.8 mg/mL of 42 kDa AGP ~ 20 uM
    lc_brainprot <- fixed(log(5000));   label("Total brain drug-binding protein concentration, PT (nM)")                    # Table II; calibrated by simulation against the observed Kp,brain (Results)

    # ---- In vitro dissociation constants governing the arterial pool.
    lkd_alb_invitro <- fixed(log(42000)); label("KD of drug binding to albumin in vitro, KA in vitro (nM)")                 # Table II (ref 13); 42 uM
    lkd_agp_invitro <- fixed(log(1200));  label("KD of drug binding to AGP in vitro, KG in vitro (nM)")                     # Table II (ref 14); 1.2 uM

    # ---- Brain physiologic volumes. Only the VP/VT ratio enters the
    # ODEs, so the L/kg normalisation cancels.
    lv_brain_vascular      <- fixed(log(0.01)); label("Brain capillary (vascular plasma) volume, VP (L/kg)")                # Table II (ref 15)
    lv_brain_extravascular <- fixed(log(0.7));  label("Brain extravascular volume, VT (L/kg)")                              # Table II (ref 15)

    # ---- Arterial input. Constant total plasma drug concentration as
    # during an IV infusion (the paper's steady-state model).
    c_plasma_ss <- fixed(100);    label("Constant total plasma drug concentration (IV-infusion steady state), LT0 (nM)")    # Table II (ref 11)
  })

  model({
    # ---- Back-transform the log-scale mechanistic constants ----------
    # One population parameter per expression, as rxode2's mu-referencing
    # parser requires. kmet_brain and c_plasma_ss are already on the
    # natural scale and are used directly below.
    koff_agp       <- exp(lkoff_agp)                # k1
    kon_agp        <- exp(lkon_agp)                 # k2
    kinflux_bbb    <- exp(lkinflux_bbb)             # k3
    kefflux_bbb    <- exp(lkefflux_bbb)             # k4
    kon_brainprot  <- exp(lkon_brainprot)           # k5
    koff_brainprot <- exp(lkoff_brainprot)          # k6
    koff_alb       <- exp(lkoff_alb)                # k7
    kon_alb        <- exp(lkon_alb)                 # k8
    kflow_brain    <- exp(lkflow_brain)             # k10
    c_alb          <- exp(lc_alb)                   # AF
    c_agp          <- exp(lc_agp)                   # GT0
    c_brainprot    <- exp(lc_brainprot)             # PT
    kd_alb_invitro <- exp(lkd_alb_invitro)          # KA in vitro
    kd_agp_invitro <- exp(lkd_agp_invitro)          # KG in vitro
    v_brain_vascular      <- exp(lv_brain_vascular)         # VP
    v_brain_extravascular <- exp(lv_brain_extravascular)    # VT

    # Total plasma drug concentration driving the brain model. Constant,
    # as in an IV infusion (Methods, "Steady State Model").
    LT0 <- c_plasma_ss

    # ---- Systemic arterial compartment (Eqs 1-4). Drug distributes
    # instantaneously between free, albumin-bound and AGP-bound pools at
    # the in vitro affinities. Free albumin is approximated by total
    # albumin because albumin is orders of magnitude more abundant than
    # drug (Methods).
    bterm <- 1 + c_alb / kd_alb_invitro + (c_agp + LT0) / kd_agp_invitro
    GL0   <- 0.5 * kd_agp_invitro *
      (bterm - sqrt(bterm^2 - 4 * c_agp * LT0 / kd_agp_invitro^2))     # Eq 1
    GF0   <- c_agp - GL0                                               # Eq 2
    fterm <- 1 + c_alb / kd_alb_invitro + GF0 / kd_agp_invitro
    LF0   <- LT0 / fterm                                               # Eq 3
    AL0   <- (c_alb * LT0 / kd_alb_invitro) / fterm                    # Eq 4

    # ---- Brain capillary and brain compartments (Eqs 18-24) ----------
    # The paper solved the imipramine case algebraically at steady state
    # (Eqs 5-11). The same system is encoded here as the differential
    # form so that the steady state is reached by integration; the two
    # agree to numerical precision (see the vignette).
    #
    # Volume ratio converting between the capillary plasma space and the
    # brain extravascular space across the BBB.
    vratio <- v_brain_extravascular / v_brain_vascular

    # Eq 18: AGP-bound drug in the brain capillary.
    d/dt(brain_vascular_drug_agp) <-
      kflow_brain * GL0 +
      kon_agp * brain_vascular_drug_free * brain_vascular_agp_free -
      koff_agp * brain_vascular_drug_agp -
      kflow_brain * brain_vascular_drug_agp

    # Eq 19: albumin-bound drug in the brain capillary.
    d/dt(brain_vascular_drug_alb) <-
      kflow_brain * AL0 +
      kon_alb * c_alb * brain_vascular_drug_free -
      koff_alb * brain_vascular_drug_alb -
      kflow_brain * brain_vascular_drug_alb

    # Eq 20: free (bioavailable) drug in the brain capillary.
    d/dt(brain_vascular_drug_free) <-
      kflow_brain * LF0 +
      koff_agp * brain_vascular_drug_agp +
      koff_alb * brain_vascular_drug_alb +
      kefflux_bbb * brain_extravascular_drug_free * vratio -
      kon_agp * brain_vascular_drug_free * brain_vascular_agp_free -
      kon_alb * c_alb * brain_vascular_drug_free -
      kinflux_bbb * brain_vascular_drug_free -
      kflow_brain * brain_vascular_drug_free

    # Eq 21: free drug in brain.
    d/dt(brain_extravascular_drug_free) <-
      kinflux_bbb * brain_vascular_drug_free / vratio -
      (kefflux_bbb + kmet_brain) * brain_extravascular_drug_free -
      kon_brainprot * brain_extravascular_drug_free *
        brain_extravascular_prot_free +
      koff_brainprot * brain_extravascular_drug_prot

    # Eq 22: free brain cytoplasmic binding protein.
    d/dt(brain_extravascular_prot_free) <-
      koff_brainprot * brain_extravascular_drug_prot -
      kon_brainprot * brain_extravascular_drug_free *
        brain_extravascular_prot_free

    # Eq 23: drug bound to brain cytoplasmic protein.
    d/dt(brain_extravascular_drug_prot) <-
      kon_brainprot * brain_extravascular_drug_free *
        brain_extravascular_prot_free -
      koff_brainprot * brain_extravascular_drug_prot

    # Eq 24: free AGP in the brain capillary.
    d/dt(brain_vascular_agp_free) <-
      kflow_brain * GF0 +
      koff_agp * brain_vascular_drug_agp -
      kon_agp * brain_vascular_drug_free * brain_vascular_agp_free -
      kflow_brain * brain_vascular_agp_free

    # ---- Initial conditions (Methods, non-steady-state model) --------
    # Binding proteins start fully unoccupied; all drug states start at 0.
    brain_vascular_agp_free(0)       <- c_agp
    brain_extravascular_prot_free(0) <- c_brainprot

    # ---- Observations -------------------------------------------------
    # Cc is the total drug concentration in the systemic (arterial)
    # plasma, the paper's LT0. The remaining outputs are the paper's
    # named variables (Table I).
    #
    # The ratio diagnostics Kp,brain = (LM + PL)/LT0 (Eq 12),
    # kp,uu,brain = LM/LF, and fu,plasma = LF/LT0 are computed in
    # post-processing rather than here, matching the sibling propranolol
    # model where LT0 can be 0 (see the vignette).
    Cc                <- LT0
    Cu_plasma_invitro <- LF0                              # C_u,plasma,in vitro
    Cu_plasma_invivo  <- brain_vascular_drug_free         # C_u,plasma,in vivo
    Cu_brain          <- brain_extravascular_drug_free    # C_u,brain,in vivo
    Cbrain_total      <- brain_extravascular_drug_free +
                         brain_extravascular_drug_prot    # LM + PL
  })
}
