Pardridge_2023_propranolol_pbpk <- function() {
  description <- paste(
    "PBPK (brain, semi-mechanistic). Pardridge 2023 Pharmaceutical",
    "Research partly-flow / partly-compartmental model of propranolol",
    "delivery to brain that resolves plasma-protein binding as an",
    "explicit kinetic process inside the brain capillary. Seven brain",
    "states carry albumin-bound, alpha-1-acid-glycoprotein (AGP)-bound",
    "and free drug in the brain capillary, free AGP in the capillary,",
    "free drug in brain, and free / drug-occupied brain cytoplasmic",
    "binding protein. Drug crosses the blood-brain barrier by",
    "bi-directional first-order permeation and may be metabolised in",
    "brain. The model quantifies plasma-protein-mediated uptake (PMU):",
    "because propranolol dissociates from AGP far faster in vivo in the",
    "capillary (KG in vivo 19 uM) than in vitro (KG in vitro 3.3 uM),",
    "the free drug in brain is roughly twofold higher than equilibrium",
    "dialysis of plasma predicts. Arterial input is either a constant",
    "total plasma concentration (steady-state IV-infusion model,",
    "c_plasma_ss, Table II LT0 = 100 nM) or a one-compartment",
    "first-order oral profile (non-steady-state model, Eq 13; dose the",
    "depot compartment and set c_plasma_ss to 0). All parameters are",
    "FIXED at the paper's literature-sourced values; the paper reports",
    "no IIV and no residual error because it is a deterministic",
    "simulation study rather than a fit to individual data."
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
  # no canonical equivalent.  `depot` and `central` are canonical.
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
      "human parameterised from published human plasma-protein chemistry",
      "and human propranolol pharmacokinetics. The total plasma",
      "propranolol concentration of 100 nM used by the steady-state model",
      "is the pharmacologic steady-state level during continuous IV",
      "infusion for supraventricular tachycardia (Table II, ref 10)."
    ),
    dose_range     = paste(
      "Non-steady-state model: single oral dose of 80 mg propranolol in a",
      "70 kg subject, i.e. 4,600 nmol/kg, simulated over 1,440 minutes",
      "(Table II, refs 21-22). Steady-state model: constant total plasma",
      "concentration of 100 nM as during IV infusion."
    ),
    notes          = paste(
      "Species provenance is mixed and should be understood before reuse.",
      "The plasma-protein concentrations (albumin 800 uM, AGP 20 uM; and",
      "the metastatic-cancer values albumin 600 uM, AGP 70 uM), the",
      "in vitro dissociation constants, and the oral propranolol",
      "pharmacokinetics are human. The in vivo brain-capillary",
      "dissociation constants and the BBB permeation rate constants k3 /",
      "k4 derive from rat carotid-arterial-injection studies (refs 8, 12)",
      "and the brain-tissue binding kinetics k5 / k6 from rat brain",
      "(ref 17). The observed Kp,brain of 9.7 that the paper matches is",
      "also rat (ref 39, Schneck 1977). The author argues (Discussion)",
      "that rodent albumin and AGP differ substantially from the human",
      "proteins, so cross-species transfer of the binding constants is",
      "explicitly a limitation of the model rather than a validated",
      "assumption."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # All parameters are FIXED. Pardridge 2023 is a simulation study:
    # every value is taken from the cited literature or from the
    # author's own simulation-based calibration, and nothing is
    # estimated from data. Values are Table II, "Propranolol" column,
    # cross-checked against Supplementary Table S1 simulation 1.
    #
# Positive-constrained mechanistic constants are log-transformed
    # (`l` prefix) and back-transformed at the top of model(), which is
    # also what rxode2's mu-referencing parser requires: at most one
    # population parameter per expression. The two parameters that are
    # legitimately ZERO in the paper's own simulations - kmet_brain (k9)
    # and c_plasma_ss (LT0) - stay on the natural scale, because log(0)
    # is undefined and both must be settable to 0.
    # ---------------------------------------------------------------

    # ---- Drug binding to AGP (alpha-1-acid glycoprotein) in the brain
    # capillary. The in vivo dissociation constant KG in vivo = k1/k2 =
    # 19 uM greatly exceeds KG in vitro = 3.3 uM; this enhanced
    # dissociation is the mechanism of plasma-protein-mediated uptake.
    lkoff_agp <- fixed(log(1140));      label("Drug dissociation rate from AGP in brain capillary, k1 (1/min)")             # Table II (ref 8); k1 = k2 * KG,in vivo = 0.06 * 19,000
    lkon_agp  <- fixed(log(0.06));      label("Drug association rate with AGP in brain capillary, k2 (1/(nM*min))")         # Table II (ref 16); 0.06 /nM/min = 1e6 /M/s

    # ---- Drug binding to albumin in the brain capillary. For
    # propranolol KA in vivo = k7/k8 = 290 uM equals KA in vitro, i.e.
    # there is NO enhanced dissociation from albumin (contrast
    # imipramine).
    lkoff_alb <- fixed(log(1740));      label("Drug dissociation rate from albumin in brain capillary, k7 (1/min)")         # Table II (ref 8); k7 = k8 * KA,in vivo = 0.006 * 290,000
    lkon_alb  <- fixed(log(0.006));     label("Drug association rate with albumin in brain capillary, k8 (1/(nM*min))")     # Table II (ref 18); 0.006 /nM/min = 1e5 /M/s

    # ---- Bi-directional blood-brain barrier permeation. k3 is derived
    # from the PS/F ratio of 1.1 with k10 = 60 /min; k4 = k3 * (VP/VT)
    # assumes symmetric transport (0.943 = 66 * 0.01 / 0.7).
    lkinflux_bbb <- fixed(log(66));     label("Drug influx rate constant plasma to brain across the BBB, k3 (1/min)")       # Table II (ref 8); k3 = k10 * PS/F = 60 * 1.1
    lkefflux_bbb <- fixed(log(0.943));  label("Drug efflux rate constant brain to plasma across the BBB, k4 (1/min)")       # Table II (ref 8); k4 = k3 * VP/VT = 66 * 0.01/0.7

    # ---- Drug binding to brain cytoplasmic protein. KP = k6/k5 = 87 nM.
    lkon_brainprot  <- fixed(log(0.006)); label("Drug association rate with brain binding protein, k5 (1/(nM*min))")        # Table II; 0.006 /nM/min = 1e5 /M/s, calibrated by simulation (Results; Table S1 sims 7/10 bracket it)
    lkoff_brainprot <- fixed(log(0.52));  label("Drug dissociation rate from brain binding protein, k6 (1/min)")            # Table II (ref 17), measured in vivo

    # ---- Brain metabolism and brain capillary plasma flow.
    kmet_brain   <- fixed(0);           label("Drug metabolism rate constant in brain, k9 (1/min)")                         # Table II (ref 17): negligible brain metabolism of propranolol
    lkflow_brain <- fixed(log(60));     label("Brain capillary plasma flow rate constant, k10 (1/min)")                     # Table II (ref 15); 1-second capillary transit time

    # ---- Binding-protein concentrations.
    lc_alb       <- fixed(log(800000)); label("Total plasma / capillary albumin concentration, AF (nM)")                    # Table II (ref 12); 5.4 g/100 mL of 67 kDa albumin ~ 800 uM
    lc_agp       <- fixed(log(20000));  label("Total plasma / capillary AGP concentration, GT0 (nM)")                       # Table II (ref 12); 0.8 mg/mL of 42 kDa AGP ~ 20 uM
    lc_brainprot <- fixed(log(5000));   label("Total brain drug-binding protein concentration, PT (nM)")                    # Table II; calibrated by simulation to the observed Kp,brain (Results; Fig 2)

    # ---- In vitro dissociation constants governing the arterial pool.
    lkd_alb_invitro <- fixed(log(290000)); label("KD of drug binding to albumin in vitro, KA in vitro (nM)")                # Table II (ref 8); 290 uM
    lkd_agp_invitro <- fixed(log(3300));   label("KD of drug binding to AGP in vitro, KG in vitro (nM)")                    # Table II (ref 8); 3.3 uM

    # ---- Brain physiologic volumes. Only the VP/VT ratio enters the
    # ODEs, so the L/kg normalisation cancels.
    lv_brain_vascular      <- fixed(log(0.01)); label("Brain capillary (vascular plasma) volume, VP (L/kg)")                # Table II (ref 15)
    lv_brain_extravascular <- fixed(log(0.7));  label("Brain extravascular volume, VT (L/kg)")                              # Table II (ref 15)

    # ---- Arterial input, steady-state model. Constant total plasma
    # drug concentration as during an IV infusion. Set to 0 to run the
    # non-steady-state oral model on its own (Eq 13).
    c_plasma_ss <- fixed(100);    label("Constant total plasma drug concentration (IV-infusion steady state), LT0 (nM)")    # Table II (ref 10)

    # ---- Arterial input, non-steady-state model (Eq 13). One-
    # compartment first-order oral absorption and elimination. Dose the
    # depot compartment with s = 4,600 nmol/kg (80 mg in a 70 kg
    # subject) to reproduce Fig 3.
    lfdepot <- fixed(log(0.3));    label("Fractional oral bioavailability, b (unitless)")                                   # Table II (ref 21)
    lka     <- fixed(log(0.023));  label("Oral absorption rate constant, k (1/min)")                                        # Table II (ref 22)
    lkel    <- fixed(log(0.0027)); label("Systemic elimination rate constant, d (1/min)")                                   # Table II (ref 22)
    lvc     <- fixed(log(5.0));    label("Systemic volume of distribution, V (L/kg)")                                       # Table II (ref 21)
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

    # ---- Systemic (arterial) pharmacokinetics ------------------------
    fdepot <- exp(lfdepot)
    ka     <- exp(lka)
    kel    <- exp(lkel)
    vc     <- exp(lvc)

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    f(depot)      <- fdepot

    # Total plasma drug concentration driving the brain model. The
    # central/vc term is exactly Eq 13,
    #   LT0[t] = b*s*k / (V*(k-d)) * (exp(-d*t) - exp(-k*t)),
    # for a single oral dose s to the depot; c_plasma_ss adds the
    # constant IV-infusion level of the steady-state model. Use one or
    # the other: keep c_plasma_ss = 100 with no dose for Table III, or
    # set c_plasma_ss = 0 and dose the depot for Fig 3 / Table V.
    LT0 <- c_plasma_ss + central / vc

    # ---- Systemic arterial compartment (Eqs 14-17; identical in form to
    # the steady-state Eqs 1-4). Drug distributes instantaneously between
    # free, albumin-bound and AGP-bound pools at the in vitro affinities.
    # Free albumin is approximated by total albumin because albumin is
    # orders of magnitude more abundant than drug (Methods).
    bterm <- 1 + c_alb / kd_alb_invitro + (c_agp + LT0) / kd_agp_invitro
    GL0   <- 0.5 * kd_agp_invitro *
      (bterm - sqrt(bterm^2 - 4 * c_agp * LT0 / kd_agp_invitro^2))     # Eq 14
    GF0   <- c_agp - GL0                                               # Eq 15
    fterm <- 1 + c_alb / kd_alb_invitro + GF0 / kd_agp_invitro
    LF0   <- LT0 / fterm                                               # Eq 16
    AL0   <- (c_alb * LT0 / kd_alb_invitro) / fterm                    # Eq 17

    # ---- Brain capillary and brain compartments (Eqs 18-24) ----------
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
    # kp,uu,brain = LM/LF, and fu,plasma = LF/LT0 are deliberately NOT
    # computed here: LT0 is exactly 0 at t = 0 in the oral simulation, so
    # they would emit 0/0. Compute them from these outputs in
    # post-processing (see the vignette).
    Cc                <- LT0
    Cu_plasma_invitro <- LF0                              # C_u,plasma,in vitro
    Cu_plasma_invivo  <- brain_vascular_drug_free         # C_u,plasma,in vivo
    Cu_brain          <- brain_extravascular_drug_free    # C_u,brain,in vivo
    Cbrain_total      <- brain_extravascular_drug_free +
                         brain_extravascular_drug_prot    # LM + PL
  })
}
