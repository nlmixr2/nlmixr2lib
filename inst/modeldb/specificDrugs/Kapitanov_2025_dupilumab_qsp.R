Kapitanov_2025_dupilumab_qsp <- function() {
  description <- paste(
    "QSP. Two-compartment physiologically-inspired PKRO (piPKRO) model of",
    "dupilumab (anti-IL4R IgG4) with full-binding TMDD in both the central",
    "and lumped-interstitial peripheral compartment (Kapitanov 2025 Case",
    "Study 2, calibrated to digitised first-in-human single-dose IV PK in",
    "healthy volunteers from Li 2015 at 1, 3, 8, and 12 mg/kg). The piPK",
    "framework replaces the classical (CL, Q, V1, V2) parameterisation",
    "with the physiological (V1, V2 = 13 L interstitial, t_half, t_dist,",
    "P_dist) macroparameter set of Kapitanov 2025 Eqs 3-7: drug",
    "elimination occurs at the same first-order rate in every compartment",
    "(k_el = log(2)/t_half), the peripheral volume is fixed at the total-",
    "body interstitial volume, and the peripheral partition coefficient",
    "P_dist is a first-class parameter. IL4R is expressed in both",
    "compartments with a bottom-up-derived concentration (Table S5;",
    "central 0.00605 nM, peripheral 0.127 nM); free drug + free receptor",
    "-> drug-receptor complex is described by mass-action forward binding",
    "with k_on = 1e-3 /(nM*s) and k_off = 3.3e-5 /s (K_D = 33 pM per",
    "Kuo 2010; Kapitanov 2025 ref 41). Complex degrades at the same rate",
    "as the free receptor (Section 3.1). Deterministic single-subject",
    "typical-value simulation; no IIV / residual error is reported by the",
    "paper for Case Study 2 and none is encoded here (fixed(0) alternative",
    "was not used because the paper reports no baseline IIV structure).",
    sep = " "
  )
  reference <- paste(
    "Kapitanov GI, Flowers D, Marcantonio DH, Lezon TR, Apgar JF, Hua F.",
    "A Tutorial on the Development of a Physiologically Inspired PKRO",
    "Model for Monoclonal Antibodies. CPT: Pharmacometrics & Systems",
    "Pharmacology. 2025. doi:10.1002/psp4.70160.",
    "Case Study 2 (2-compartment piPKRO fitted to digitised dupilumab",
    "single-dose IV PK from Li E et al., J Clin Pharmacol 2015 [Kapitanov",
    "2025 ref 35]).",
    "Bottom-up IL4R target burden derived per Marcantonio DH et al. 2020",
    "(Kapitanov 2025 ref 34). Dupilumab-IL4R affinity K_D = 33 pM from",
    "Kuo 2010 (Kapitanov 2025 ref 41). IL4R turnover rate from Andrews",
    "et al. (Kapitanov 2025 ref 40).",
    sep = " "
  )
  vignette <- "Kapitanov_2025_dupilumab"

  # Non-canonical compartments for the multi-compartment TMDD binding.
  # target_<location> is registered canonically only for csf / isf; the
  # piPKRO framework tracks free target and drug-target complex in each
  # PK compartment, so these are declared paper-specific.
  paper_specific_compartments <- c(
    "target_central", "complex_central",
    "target_peripheral1", "complex_peripheral1"
  )

  units <- list(
    time          = "h",
    dosing        = paste(
      "Dupilumab dose into the central compartment must be in mg.",
      "The bioavailability multiplier converts amt (mg) to nmol using",
      "the standard IgG4 antibody molecular weight (assumed 147 kDa;",
      "not stated in Kapitanov 2025). Case Study 2 simulates 1, 3, 8,",
      "and 12 mg/kg body-weight-normalised IV doses at 70 kg (i.e., 70,",
      "210, 560, and 840 mg absolute).",
      sep = " "
    ),
    concentration = paste(
      "Free dupilumab plasma concentration Cc is reported in mg/L to",
      "match Kapitanov 2025 Figure 4 (y-axis Concentration [ug/mL]).",
      "The underlying state concentration in nM (Cc_nM = central / vc)",
      "is also emitted for RO-related downstream comparisons. IL4R",
      "receptor and drug-receptor complex concentrations remain in",
      "nM. Receptor-occupancy outputs (RO_c, RO_p) are unitless",
      "fractions in [0, 1].",
      sep = " "
    )
  )

  # No data-carried covariates: Case Study 2 is a deterministic mean-fit
  # to digitised first-in-human PK in healthy volunteers with body
  # weight assumed at 70 kg (Kapitanov 2025 Section 3.1). Dose is a per-
  # simulation input, not a covariate.
  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    age_range      = "healthy adult volunteers",
    weight_range   = "assumed 70 kg for the deterministic simulation",
    weight_median  = "70 kg (Kapitanov 2025 Section 2.4 standard human)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Healthy adult volunteers from the dupilumab first-in-human",
      "single-ascending-dose IV study reported in Li E et al. 2015",
      "(Kapitanov 2025 ref 35). Mean PK profiles at 1, 3, 8, and 12",
      "mg/kg were digitised from Li 2015 for the Case Study 2 fit.",
      sep = " "
    ),
    dose_range     = paste(
      "Single IV bolus at 1, 3, 8, and 12 mg/kg (70, 210, 560, 840 mg",
      "at 70 kg).",
      sep = " "
    ),
    regions        = NA_character_,
    notes          = paste(
      "Kapitanov 2025 does not perform a per-subject fit for Case",
      "Study 2; the paper reports only typical values of the linear",
      "PK-related parameters after calibration to digitised mean PK",
      "profiles. No IIV / OMEGA and no residual-error / SIGMA are",
      "reported. The IL4R central and peripheral target concentrations",
      "(C_R,1 = 0.00605 nM, C_R,2 = 0.127 nM) were held fixed at the",
      "bottom-up estimates from Kapitanov 2025 Table S5 (Marcantonio",
      "2020 methodology).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Compartment volumes (Kapitanov 2025 Table S4, Case Study 2 fitted).
    # V_1 = central volume (L); V_2 = peripheral (lumped interstitial)
    # volume (L). V_2 was fitted rather than fixed at the 13 L
    # physiological interstitial default in this case study.
    # =====================================================================
    lvc <- fixed(log(2.29))   ; label("Central volume V_1 (L; Kapitanov 2025 Table S4 fitted 2-cpt)")     # Kapitanov 2025 Table S4, Fitted 2-compt model column: V_1 = 2.29 L
    lvp <- fixed(log(12.4))   ; label("Peripheral volume V_2 (L; Kapitanov 2025 Table S4 fitted 2-cpt)")  # Kapitanov 2025 Table S4, Fitted 2-compt model column: V_2 = 12.4 L

    # =====================================================================
    # piPK macroparameters (Kapitanov 2025 Eqs 5-7; Table S4 Case 2).
    # t_half = terminal-linear-elimination half-life (Eq 5): k_el =
    #   log(2)/t_half.
    # t_dist = central-peripheral distribution half-time (Eq 7):
    #   k_12 + k_21 = log(2)/t_dist.
    # P_dist = peripheral/central concentration ratio at quasi-steady-
    #   state in the absence of TMDD (Eq 6):
    #   P_dist = V_1 * k_12 / (V_2 * k_21).
    # Solving Eqs 6-7 jointly:
    #   k_21 = log(2) * V_1 / (t_dist * (V_1 + P_dist * V_2))
    #   k_12 = P_dist * V_2 * k_21 / V_1
    # =====================================================================
    lthalf <- fixed(log(32.8 * 24)) ; label("Linear elimination half-life t_half (h; 32.8 d converted; Kapitanov 2025 Table S4)")  # Kapitanov 2025 Table S4 Fitted 2-compt: t_half = 32.8 d; stored in h so k_el = log(2)/thalf gives 1/h
    ltdist <- fixed(log(54.3))      ; label("Central-peripheral distribution half-time t_dist (h; Kapitanov 2025 Table S4)")     # Kapitanov 2025 Table S4 Fitted 2-compt: t_dist,12 = 54.3 h
    lpdist <- fixed(log(0.352))     ; label("Peripheral partition coefficient P_dist (unitless; Kapitanov 2025 Table S4)")        # Kapitanov 2025 Table S4 Fitted 2-compt: P_dist,12 = 0.352

    # =====================================================================
    # IL4R receptor concentrations (Kapitanov 2025 Table S4 & Table S5).
    # C_R,1 = central compartment target concentration (nM).
    # C_R,2 = peripheral compartment target concentration (nM).
    # Derived from bottom-up cell-count sums per Marcantonio 2020
    # methodology (Kapitanov 2025 Table S5). Fixed structural inputs to
    # the piPKRO model; NOT fitted parameters.
    # =====================================================================
    lc_r_c <- fixed(log(0.00605))   ; label("Central compartment IL4R concentration C_R,1 (nM; bottom-up from Table S5)")    # Kapitanov 2025 Table S4 & Table S5 bottom-up sum for central compartment (lymphocyte + myeloid cells)
    lc_r_p <- fixed(log(0.127))     ; label("Peripheral compartment IL4R concentration C_R,2 (nM; bottom-up from Table S5)") # Kapitanov 2025 Table S4 & Table S5 bottom-up sum for peripheral compartment

    # =====================================================================
    # IL4R degradation and dupilumab-IL4R binding kinetics (Kapitanov
    # 2025 Table S4). Converted from the paper's mixed h / s / nM units
    # to the model's uniform h time unit and nM concentration unit.
    # =====================================================================
    lkdeg <- fixed(log(1))                 ; label("IL4R free-receptor and drug-receptor complex degradation rate k_deg (1/h; Kapitanov 2025 Table S4)")  # Kapitanov 2025 Table S4: k_deg = 1 /h (fixed from Kapitanov 2025 ref 40, Andrews et al.); complex degrades at the same rate per Section 3.1
    lkon  <- fixed(log(3.6))               ; label("Dupilumab-IL4R association rate k_on (1/(nM*h); converted from Table S4 1e-3 /(nM*s) * 3600)")        # Kapitanov 2025 Table S4: k_on = 1e-3 /(nM*s); converted to /(nM*h) by multiplying by 3600
    lkoff <- fixed(log(3.3e-5 * 3600))     ; label("Dupilumab-IL4R dissociation rate k_off (1/h; converted from Table S4 3.3e-5 /s * 3600)")               # Kapitanov 2025 Table S4: k_off = 3.3e-5 /s; K_D = k_off/k_on = 33 pM matches Kuo 2010 (ref 41)
  })

  model({
    # =====================================================================
    # 1. Fixed physical constants. Dupilumab molecular weight is not
    # stated in Kapitanov 2025; the standard IgG4 monoclonal-antibody
    # value of 147 kDa (147000 g/mol) is used for the mg/L <-> nM
    # conversion below. Kapitanov 2025 Figure 4 y-axis uses ug/mL (=
    # mg/L), implying the authors also assumed a standard IgG4 MW.
    # =====================================================================
    mw_dupi <- 147000  # g/mol   # standard human IgG4 mAb MW; not stated in Kapitanov 2025 but consistent with Figure 4 y-axis unit
    mgL_per_nM <- mw_dupi / 1e6  # 147/1000 mg/L per nM: c(mg/L) = c(nM) * MW / 1e6

    # =====================================================================
    # 2. Individual (typical-value) parameters. No IIV encoded because
    # the paper reports no per-subject fit or population variance.
    # =====================================================================
    vc <- exp(lvc)         # L
    vp <- exp(lvp)         # L
    thalf <- exp(lthalf)   # h
    tdist <- exp(ltdist)   # h
    pdist <- exp(lpdist)   # unitless
    c_r_c <- exp(lc_r_c)   # nM
    c_r_p <- exp(lc_r_p)   # nM
    kdeg  <- exp(lkdeg)    # 1/h
    kon   <- exp(lkon)     # 1/(nM*h)
    koff  <- exp(lkoff)    # 1/h

    # =====================================================================
    # 3. Derived rate constants. Kapitanov 2025 Eqs 5-7 (macroparameter
    # to rate-constant conversion). k_el is the same in central and
    # peripheral (piPK framework definition, Section 2.1).
    # =====================================================================
    kel <- log(2) / thalf
    # Solve Eqs 6-7 for k_12, k_21:
    #   k_12 + k_21 = log(2) / tdist         (Eq 7)
    #   V_1 * k_12  = pdist * V_2 * k_21     (Eq 6, rearranged)
    denom_tp <- vc + pdist * vp
    k21 <- log(2) * vc / (tdist * denom_tp)
    k12 <- pdist * vp * k21 / vc

    # Target zero-order synthesis rate constant (Kapitanov 2025 Eq 23):
    # at steady state before dosing, k_syn * V = k_deg * R_ss (amounts).
    # In concentration space, k_syn (nM/h) = k_deg * C_R.
    ksyn_c <- kdeg * c_r_c   # nM/h in central
    ksyn_p <- kdeg * c_r_p   # nM/h in peripheral

    # =====================================================================
    # 4. State-to-concentration conversions. State amounts are in nmol
    # (drug), so state / volume gives concentration in nM.
    # =====================================================================
    Cd_c   <- central              / vc  # free drug in central (nM)
    Cd_p   <- peripheral1          / vp  # free drug in peripheral (nM)
    Cr_c   <- target_central       / vc  # free receptor in central (nM)
    Cr_p   <- target_peripheral1   / vp  # free receptor in peripheral (nM)
    Cdr_c  <- complex_central      / vc  # drug-receptor complex in central (nM)
    Cdr_p  <- complex_peripheral1  / vp  # drug-receptor complex in peripheral (nM)

    # =====================================================================
    # 5. ODE system (Kapitanov 2025 Eqs 3-4 for drug transport, Eqs
    # 11-16 for the full binding ODE system). State amounts are in
    # nmol; concentration-form binding terms (k_on * C_D * C_R) are
    # multiplied by compartment volume V_i to give amount-form fluxes
    # (nmol/h). Elimination rate k_el is applied in both central and
    # peripheral for free drug per the piPK framework.
    # =====================================================================

    # Free drug (Eqs 11, 14 -- amount form via multiplication by V_i):
    d/dt(central)     <- -k12 * central - kel * central + k21 * peripheral1 -
                          vc * kon * Cd_c * Cr_c + koff * complex_central
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1 - kel * peripheral1 -
                          vp * kon * Cd_p * Cr_p + koff * complex_peripheral1

    # Free target (Eqs 12, 15). Synthesis in amount-form is ksyn * V.
    d/dt(target_central)     <-  ksyn_c * vc - kdeg * target_central -
                                  vc * kon * Cd_c * Cr_c + koff * complex_central
    d/dt(target_peripheral1) <-  ksyn_p * vp - kdeg * target_peripheral1 -
                                  vp * kon * Cd_p * Cr_p + koff * complex_peripheral1

    # Drug-receptor complex (Eqs 13, 16). Assumed to degrade at k_deg
    # per Section 3.1 ("The dupilumab:IL4R complex was assumed to
    # eliminate at the same rate as IL4R.").
    d/dt(complex_central)     <-  vc * kon * Cd_c * Cr_c - koff * complex_central -
                                   kdeg * complex_central
    d/dt(complex_peripheral1) <-  vp * kon * Cd_p * Cr_p - koff * complex_peripheral1 -
                                   kdeg * complex_peripheral1

    # =====================================================================
    # 6. Initial conditions. Drug and complex compartments start at 0.
    # Target starts at the steady-state amount k_syn * V / k_deg = C_R
    # * V.
    # =====================================================================
    target_central(0)     <- c_r_c * vc
    target_peripheral1(0) <- c_r_p * vp

    # =====================================================================
    # 7. Bioavailability-style unit conversion for IV dosing. Event-
    # table `amt` is in mg (clinical unit). State amount is in nmol.
    # Delivered amount (nmol) = amt (mg) * 1e-3 (g/mg) * (1 mol / MW g)
    # * 1e9 (nmol/mol) = amt * 1e6 / MW. For MW = 147000 g/mol this
    # gives amt * 6.803 nmol/mg.
    # =====================================================================
    f(central) <- 1e6 / mw_dupi

    # =====================================================================
    # 8. Observation variables. Cc is dupilumab plasma concentration in
    # mg/L to match Kapitanov 2025 Figure 4 axes. Cc_nM is the free-
    # drug concentration in nM (used inside the model for binding
    # kinetics). RO_c and RO_p are the receptor-occupancy fractions.
    # =====================================================================
    Cc      <- Cd_c * mgL_per_nM              # dupilumab in mg/L
    Cc_nM   <- Cd_c                            # dupilumab in nM (matches paper Table S4 units for binding)
    Cc_p    <- Cd_p * mgL_per_nM              # peripheral dupilumab in mg/L
    Cc_p_nM <- Cd_p                            # peripheral dupilumab in nM

    RO_c    <- Cdr_c / (Cr_c + Cdr_c)         # central receptor occupancy (unitless)
    RO_p    <- Cdr_p / (Cr_p + Cdr_p)         # peripheral receptor occupancy (unitless)

    # No residual error is reported by Kapitanov 2025 for Case Study 2
    # (deterministic mean-fit to digitised PK). The model is intended
    # for simulation rather than estimation.
  })
}
