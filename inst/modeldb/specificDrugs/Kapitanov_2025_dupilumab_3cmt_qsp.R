Kapitanov_2025_dupilumab_3cmt_qsp <- function() {
  description <- paste(
    "QSP. Three-compartment physiologically-inspired PKRO (piPKRO) model",
    "of dupilumab (anti-IL4R IgG4) with full-binding TMDD in central, a",
    "lumped-interstitial peripheral, and an inflamed-skin site-of-action",
    "(SoA) compartment (Kapitanov 2025 Case Study 3 Approach 2). Extends",
    "the 2-compartment Case Study 2 model by adding a peripheral2",
    "compartment representing 50 percent of the total interstitial skin",
    "volume (0.563 L), with dupilumab distribution parameterised by a",
    "fixed skin partition coefficient (P_dist,13 = 0.3) and a fixed",
    "central-to-skin distribution half-time (t_dist,13 = 30 h). IL4R",
    "concentration in the inflamed-skin SoA compartment is set to 2.02",
    "nM (Table S4) reflecting higher immune-cell infiltrate and higher",
    "per-cell IL4R expression in atopic-dermatitis skin (Table S5). All",
    "other parameters (t_half, P_dist,12, C_R,1, C_R,2, k_deg, k_on,",
    "k_off) are carried unchanged from the Case Study 2 fit. The lumped",
    "peripheral volume V_2 is kept at 12.4 L per Kapitanov 2025 Table S4",
    "Approach 2 (note: Section 3.3 narrates that V_2 was decreased by",
    "subtracting the SoA volume, but Table S4 Approach 2 reports V_2 =",
    "12.4 L identical to Case 2 -- this file follows Table S4). The",
    "central-peripheral distribution half-time t_dist,12 = 55.9 h is",
    "recomputed for the 3-cpt structure per Eq 22. Case",
    "Study 3 illustrates that a specific SoA compartment with high local",
    "target burden and slow distribution can lead to substantial local",
    "drug depletion and lower receptor occupancy at the SoA than in the",
    "central and lumped-peripheral compartments (Kapitanov 2025 Figure",
    "5). Deterministic single-subject typical-value simulation; no IIV /",
    "residual error is reported by the paper.",
    sep = " "
  )
  reference <- paste(
    "Kapitanov GI, Flowers D, Marcantonio DH, Lezon TR, Apgar JF, Hua F.",
    "A Tutorial on the Development of a Physiologically Inspired PKRO",
    "Model for Monoclonal Antibodies. CPT: Pharmacometrics & Systems",
    "Pharmacology. 2025. doi:10.1002/psp4.70160.",
    "Case Study 3 Approach 2 (3-compartment piPKRO with fixed P_dist and",
    "t_dist to a skin site-of-action compartment).",
    "See modellib('Kapitanov_2025_dupilumab_qsp') for the parent 2-",
    "compartment Case Study 2 model this extends.",
    sep = " "
  )
  vignette <- "Kapitanov_2025_dupilumab"

  # peripheral2 is the canonical name for the second peripheral (SoA)
  # compartment; target_peripheral2 and complex_peripheral2 are non-
  # canonical paper-mechanistic extensions of the multi-compartment TMDD
  # binding.
  paper_specific_compartments <- c(
    "target_central", "complex_central",
    "target_peripheral1", "complex_peripheral1",
    "target_peripheral2", "complex_peripheral2"
  )

  units <- list(
    time          = "h",
    dosing        = paste(
      "Dupilumab dose into the central compartment must be in mg.",
      "Bioavailability multiplier converts amt (mg) to nmol using MW",
      "147 kDa (standard IgG4 mAb; not stated in Kapitanov 2025). Case",
      "Study 3 simulates 1, 3, 8, 12 mg/kg body-weight-normalised IV",
      "doses at 70 kg (70, 210, 560, 840 mg absolute).",
      sep = " "
    ),
    concentration = paste(
      "Free dupilumab plasma concentration Cc reported in mg/L (paper",
      "Figure 5 y-axis). Cc_p and Cc_soa are the free peripheral and",
      "SoA dupilumab concentrations in mg/L; the underlying nM state",
      "concentrations are also emitted (Cc_nM, Cc_p_nM, Cc_soa_nM).",
      "Receptor occupancy outputs (RO_c, RO_p, RO_soa) are unitless",
      "fractions in [0, 1].",
      sep = " "
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    central             = list(analyte = "dupilumab", units = NA_character_, specimen = "plasma", verified = FALSE),
    peripheral1         = list(analyte = "dupilumab", units = NA_character_, specimen = "plasma", verified = FALSE),
    peripheral2         = list(analyte = "dupilumab", units = NA_character_, specimen = "plasma", verified = FALSE),
    target_central      = list(analyte = "IL4R", units = NA_character_, specimen = "plasma", verified = FALSE),
    target_peripheral1  = list(analyte = "IL4R", units = NA_character_, specimen = "plasma", verified = FALSE),
    target_peripheral2  = list(analyte = "IL4R", units = NA_character_, specimen = "plasma", verified = FALSE),
    complex_central     = list(analyte = "dupilumab-IL4R complex", units = NA_character_, specimen = "plasma", verified = FALSE),
    complex_peripheral1 = list(analyte = "dupilumab-IL4R complex", units = NA_character_, specimen = "plasma", verified = FALSE),
    complex_peripheral2 = list(analyte = "dupilumab-IL4R complex", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    age_range      = "healthy adult volunteers (upstream PK data)",
    weight_range   = "assumed 70 kg",
    weight_median  = "70 kg (Kapitanov 2025 Section 2.4 standard human)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "The dupilumab PK layer inherits from Case Study 2, i.e. the",
      "healthy-volunteer first-in-human study of Li E et al. 2015",
      "(Kapitanov 2025 ref 35). The IL4R concentration in the SoA",
      "compartment (C_R,3 = 2.02 nM) reflects atopic-dermatitis",
      "inflamed skin per the bottom-up Table S5 estimate: elevated",
      "immune-cell infiltrate (macrophages, T cells, mast cells, DCs,",
      "eosinophils, neutrophils, keratinocytes, fibroblasts) and",
      "elevated per-cell IL4R expression. The dupilumab clinical",
      "indication that motivates Case Study 3 is moderate-to-severe",
      "atopic dermatitis.",
      sep = " "
    ),
    dose_range     = paste(
      "Single IV bolus at 1, 3, 8, and 12 mg/kg (70, 210, 560, 840 mg",
      "at 70 kg) simulated to demonstrate SoA receptor-occupancy",
      "prediction.",
      sep = " "
    ),
    regions        = NA_character_,
    notes          = paste(
      "Kapitanov 2025 does not perform a per-subject fit for Case",
      "Study 3; the parameters are structural / mechanistic",
      "assumptions on top of the Case Study 2 fitted PK. No IIV / OMEGA",
      "and no residual error / SIGMA are reported. Section 3.3",
      "explicitly notes: 'to our knowledge, no clinical RO data are",
      "available to validate our model predictions' -- the SoA",
      "compartment is an illustrative extension.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Compartment volumes (Kapitanov 2025 Table S4, Case 3 Approach 2).
    # V_1 = central; V_2 = lumped peripheral (12.4 L, unchanged from
    # Case 2, i.e. the SoA volume is NOT subtracted for direct
    # comparability with Case Study 2 Table S4 row); V_3 = SoA (skin,
    # 50 percent of total interstitial skin volume = 0.5 * 1.125 =
    # 0.5625 L, reported as 0.563 L).
    # =====================================================================
    lvc  <- fixed(log(2.29))  ; label("Central volume V_1 (L; Kapitanov 2025 Table S4 Case 3)")            # Kapitanov 2025 Table S4, 3-compt model Approach 2 column: V_1 = 2.29 L
    lvp  <- fixed(log(12.4))  ; label("Peripheral volume V_2 (L; Kapitanov 2025 Table S4 Case 3)")         # Kapitanov 2025 Table S4, 3-compt model Approach 2 column: V_2 = 12.4 L
    lvp2 <- fixed(log(0.563)) ; label("SoA volume V_3 = 50% interstitial skin (L; Kapitanov 2025 Table S4)") # Kapitanov 2025 Table S4, 3-compt model Approach 2 column: V_3 = 0.563 L; Section 3.3: 50% of 1.125 L total interstitial skin

    # =====================================================================
    # piPK macroparameters. t_half is unchanged from Case 2; t_dist,12
    # was recomputed via Eq 22 for the 3-cpt structure; P_dist,12 was
    # kept unchanged from Case 2 per Kapitanov 2025 Approach 1/2
    # narrative. t_dist,13 and P_dist,13 are the Approach-2 fixed SoA
    # macroparameters. Approach 1 alternative values (t_dist,13 = 150 h,
    # P_dist,13 = 0.352) are noted in the source-trace comments.
    # =====================================================================
    lthalf <- fixed(log(32.8 * 24)) ; label("Linear elimination half-life t_half (h; 32.8 d converted; Kapitanov 2025 Table S4)")  # Kapitanov 2025 Table S4 3-compt Approach 2: t_half = 32.8 d (unchanged from Case 2)
    ltdist <- fixed(log(55.9))      ; label("Central-peripheral1 distribution half-time t_dist,12 (h; Kapitanov 2025 Table S4)")  # Kapitanov 2025 Table S4 3-compt Approach 2: t_dist,12 = 55.9 h (recomputed via Eq 22 vs Case 2's 54.3 h)
    lpdist <- fixed(log(0.352))     ; label("Peripheral1 partition coefficient P_dist,12 (unitless; Kapitanov 2025 Table S4)")     # Kapitanov 2025 Table S4 3-compt Approach 2: P_dist,12 = 0.352 (unchanged from Case 2)

    ltdist_p2 <- fixed(log(30))     ; label("Central-SoA distribution half-time t_dist,13 (h; Kapitanov 2025 Approach 2)")   # Kapitanov 2025 Table S4 3-compt Approach 2: t_dist,13 = 30 h (Approach 1 alternative: 150 h)
    lpdist_p2 <- fixed(log(0.3))    ; label("SoA partition coefficient P_dist,13 (unitless; Kapitanov 2025 Approach 2)")     # Kapitanov 2025 Table S4 3-compt Approach 2: P_dist,13 = 0.3 (Approach 1 alternative: 0.352)

    # =====================================================================
    # IL4R receptor concentrations (Kapitanov 2025 Table S4 and Table
    # S5). C_R,3 = 2.02 nM in the skin SoA compartment reflects
    # atopic-dermatitis pathology; higher than C_R,2 = 0.127 nM by a
    # factor of about 16, driven by higher immune-cell infiltration and
    # higher IL4R expression on inflamed keratinocytes and fibroblasts
    # (Kapitanov 2025 Section 3.3 and refs 43-48).
    # =====================================================================
    lc_r_c <- fixed(log(0.00605))   ; label("Central IL4R concentration C_R,1 (nM; Table S5 bottom-up)")   # Kapitanov 2025 Table S4 & Table S5 bottom-up sum
    lc_r_p <- fixed(log(0.127))     ; label("Peripheral IL4R concentration C_R,2 (nM; Table S5 bottom-up)") # Kapitanov 2025 Table S4 & Table S5 bottom-up sum
    lc_r_soa <- fixed(log(2.02))    ; label("SoA IL4R concentration C_R,3 (nM; inflamed skin; Table S5)")   # Kapitanov 2025 Table S4 & Table S5 skin bottom-up sum for atopic-dermatitis-inflamed skin

    # =====================================================================
    # IL4R kinetics (unchanged from Case 2).
    # =====================================================================
    lkdeg <- fixed(log(1))                 ; label("IL4R free-receptor and complex degradation rate k_deg (1/h)")  # Kapitanov 2025 Table S4: k_deg = 1 /h (Kapitanov 2025 ref 40, Andrews et al.)
    lkon  <- fixed(log(3.6))               ; label("Dupilumab-IL4R association rate k_on (1/(nM*h); 1e-3 /(nM*s) * 3600)")  # Kapitanov 2025 Table S4: k_on = 1e-3 /(nM*s)
    lkoff <- fixed(log(3.3e-5 * 3600))     ; label("Dupilumab-IL4R dissociation rate k_off (1/h; 3.3e-5 /s * 3600)")         # Kapitanov 2025 Table S4: k_off = 3.3e-5 /s
  })

  model({
    # =====================================================================
    # 1. Physical constants (see 2-compartment sibling for provenance
    # discussion). MW = 147 kDa standard IgG4.
    # =====================================================================
    mw_dupi   <- 147000                    # g/mol
    mgL_per_nM <- mw_dupi / 1e6            # 0.147 mg/L per nM

    # =====================================================================
    # 2. Individual (typical-value) parameters.
    # =====================================================================
    vc   <- exp(lvc)         # L
    vp   <- exp(lvp)         # L
    vp2  <- exp(lvp2)        # L (SoA skin)

    thalf <- exp(lthalf)                  # h
    tdist <- exp(ltdist)                  # h (central-peripheral)
    pdist <- exp(lpdist)                  # unitless (peripheral)
    tdist_p2 <- exp(ltdist_p2)            # h (central-SoA)
    pdist_p2 <- exp(lpdist_p2)            # unitless (SoA)

    c_r_c   <- exp(lc_r_c)   # nM
    c_r_p   <- exp(lc_r_p)   # nM
    c_r_soa <- exp(lc_r_soa) # nM

    kdeg <- exp(lkdeg)       # 1/h
    kon  <- exp(lkon)        # 1/(nM*h)
    koff <- exp(lkoff)       # 1/h

    # =====================================================================
    # 3. Derived rate constants. Kapitanov 2025 Eqs 5-7 applied to each
    # of the two central-peripheral distributions.
    # =====================================================================
    kel <- log(2) / thalf                              # 1/h; same in all compartments (piPK framework)

    # Central <-> peripheral1 (Kapitanov 2025 Eqs 6-7)
    denom_p1 <- vc + pdist * vp
    k21 <- log(2) * vc / (tdist * denom_p1)            # 1/h
    k12 <- pdist * vp * k21 / vc                        # 1/h

    # Central <-> peripheral2 / SoA (Kapitanov 2025 Eqs 6-7 extended)
    denom_p2 <- vc + pdist_p2 * vp2
    k31 <- log(2) * vc / (tdist_p2 * denom_p2)         # 1/h
    k13 <- pdist_p2 * vp2 * k31 / vc                   # 1/h

    # Target synthesis (Kapitanov 2025 Eq 23: k_syn = R_ss * k_deg).
    ksyn_c   <- kdeg * c_r_c   # nM/h
    ksyn_p   <- kdeg * c_r_p   # nM/h
    ksyn_soa <- kdeg * c_r_soa # nM/h

    # =====================================================================
    # 4. Concentrations from state amounts. State amounts in nmol;
    # concentrations in nM.
    # =====================================================================
    Cd_c    <- central             / vc
    Cd_p    <- peripheral1         / vp
    Cd_soa  <- peripheral2         / vp2
    Cr_c    <- target_central      / vc
    Cr_p    <- target_peripheral1  / vp
    Cr_soa  <- target_peripheral2  / vp2
    Cdr_c   <- complex_central     / vc
    Cdr_p   <- complex_peripheral1 / vp
    Cdr_soa <- complex_peripheral2 / vp2

    # =====================================================================
    # 5. ODE system (Kapitanov 2025 Eqs 9-16 extended). Amount-form ODEs
    # with concentration-form binding fluxes multiplied by V_i.
    # =====================================================================

    # Free drug (Eqs 9, 10, 14 -- piPK 3-cpt with elimination in every
    # compartment and TMDD in all three).
    d/dt(central)     <- -k12 * central - k13 * central - kel * central +
                          k21 * peripheral1 + k31 * peripheral2 -
                          vc * kon * Cd_c * Cr_c + koff * complex_central
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1 - kel * peripheral1 -
                          vp * kon * Cd_p * Cr_p + koff * complex_peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2 - kel * peripheral2 -
                          vp2 * kon * Cd_soa * Cr_soa + koff * complex_peripheral2

    # Free target (Eqs 12, 15 extended to peripheral2).
    d/dt(target_central)     <-  ksyn_c   * vc  - kdeg * target_central -
                                  vc  * kon * Cd_c   * Cr_c   + koff * complex_central
    d/dt(target_peripheral1) <-  ksyn_p   * vp  - kdeg * target_peripheral1 -
                                  vp  * kon * Cd_p   * Cr_p   + koff * complex_peripheral1
    d/dt(target_peripheral2) <-  ksyn_soa * vp2 - kdeg * target_peripheral2 -
                                  vp2 * kon * Cd_soa * Cr_soa + koff * complex_peripheral2

    # Drug-receptor complex (Eqs 13, 16 extended to peripheral2). Complex
    # eliminates at k_deg per Section 3.1.
    d/dt(complex_central)     <-  vc  * kon * Cd_c   * Cr_c   - koff * complex_central -
                                   kdeg * complex_central
    d/dt(complex_peripheral1) <-  vp  * kon * Cd_p   * Cr_p   - koff * complex_peripheral1 -
                                   kdeg * complex_peripheral1
    d/dt(complex_peripheral2) <-  vp2 * kon * Cd_soa * Cr_soa - koff * complex_peripheral2 -
                                   kdeg * complex_peripheral2

    # =====================================================================
    # 6. Initial conditions (steady-state target, no drug).
    # =====================================================================
    target_central(0)     <- c_r_c   * vc
    target_peripheral1(0) <- c_r_p   * vp
    target_peripheral2(0) <- c_r_soa * vp2

    # =====================================================================
    # 7. IV dosing (bioavailability converts mg -> nmol via MW).
    # =====================================================================
    f(central) <- 1e6 / mw_dupi

    # =====================================================================
    # 8. Observations. Cc, Cc_p, Cc_soa in mg/L (paper Figure 5). Cc_nM
    # variants in nM for cross-check against binding kinetics. RO in
    # each compartment.
    # =====================================================================
    Cc        <- Cd_c   * mgL_per_nM
    Cc_nM     <- Cd_c
    Cc_p      <- Cd_p   * mgL_per_nM
    Cc_p_nM   <- Cd_p
    Cc_soa    <- Cd_soa * mgL_per_nM
    Cc_soa_nM <- Cd_soa

    RO_c   <- Cdr_c   / (Cr_c   + Cdr_c)
    RO_p   <- Cdr_p   / (Cr_p   + Cdr_p)
    RO_soa <- Cdr_soa / (Cr_soa + Cdr_soa)
  })
}
