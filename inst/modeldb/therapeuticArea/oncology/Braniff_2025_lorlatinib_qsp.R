Braniff_2025_lorlatinib_qsp <- function() {
  description <- paste(
    "QSP. Minimal MAPK/PI3K signaling pathway coupled to a shell-and-core",
    "tumor growth model for the ALK inhibitor lorlatinib, used as the",
    "mechanistic engine of a virtual-population (Vpop) selection algorithm",
    "calibrated to a phase II non-small-cell lung cancer study (Braniff",
    "2025; clinical data from Solomon 2018, NCT01970865 cohort with at",
    "least two sum-of-longest-diameter measurements, truncated at day 500).",
    "Fourteen ODE states: a two-compartment lorlatinib PK model with",
    "sequential zero- and first-order absorption and time-varying",
    "auto-induced clearance CL(t) = CLmax * (1 - exp(-kc * t)) (covariates",
    "fixed at their median values, giving a median concentration-time",
    "course; inherited from the Chen 2021 population PK model); five",
    "normalised signaling species (pALK, RAS, pERK, pAKT, pS6) driven by",
    "sigmoid ALK inhibition and a basal growth-factor-receptor activation",
    "term gfr_act; and a six-state shell-and-core tumor model in which a",
    "proliferating shell (cycling_cells) transfers cells to a necrotic core",
    "(damaged_cells1) and on through four clearance-delay compartments",
    "(damaged_cells2-damaged_cells5). Drug effect enters through reduced",
    "proliferation (driven by pS6) and increased apoptosis (driven by",
    "pAKT), combined into a single kill term kkill. The observed sum of",
    "longest diameters is the diameter of the combined shell-and-core",
    "sphere at 1e5 cells per microlitre. The paper reports Vpop SAMPLING",
    "RANGES (Tables S1, S3) rather than point estimates for the eleven",
    "sampled quantities; the ini() values here are a single documented",
    "reference plausible patient and the ranges are recorded per parameter",
    "- see the vignette for the sampling-based reproduction of the",
    "published virtual population.",
    sep = " "
  )
  reference <- paste(
    "Braniff N, Joshi T, Cassidy T, Trogdon M, Kumar R, Poels K, Allen R,",
    "Musante CJ, Shtylla B. (2025).",
    "An integrated quantitative systems pharmacology virtual population",
    "approach for calibration with oncology efficacy endpoints.",
    "CPT Pharmacometrics Syst Pharmacol 14(2):268-278.",
    "doi:10.1002/psp4.13270",
    sep = " "
  )
  vignette <- "Braniff_2025_lorlatinib_qsp"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  # The five MAPK/PI3K signaling species are paper-mechanistic QSP states
  # with no canonical analogue in inst/references/compartment-names.md (see
  # that file's "Paper-specific compartments" section); the same treatment
  # is used by Minucci_2024_CART_qsp.R for its CAR T-cell phenotype states.
  # cycling_cells and damaged_cells1..damaged_cells5 are canonical.
  paper_specific_compartments <- c("palk", "ras", "perk", "pakt", "ps6")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The signaling species are normalised to their
  # untreated baseline (1 = fully active) and are therefore unitless.
  compartmentData <- list(
    depot          = list(analyte = "lorlatinib", units = "mg", specimen = "administration site", verified = TRUE),
    central        = list(analyte = "lorlatinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1    = list(analyte = "lorlatinib", units = "mg", specimen = "plasma", verified = TRUE),
    palk           = list(analyte = "phosphorylated ALK", units = "fraction of untreated baseline", specimen = "tumor", verified = TRUE),
    ras            = list(analyte = "active RAS", units = "fraction of untreated baseline", specimen = "tumor", verified = TRUE),
    perk           = list(analyte = "phosphorylated ERK", units = "fraction of untreated baseline", specimen = "tumor", verified = TRUE),
    pakt           = list(analyte = "phosphorylated AKT", units = "fraction of untreated baseline", specimen = "tumor", verified = TRUE),
    ps6            = list(analyte = "phosphorylated S6", units = "fraction of untreated baseline", specimen = "tumor", verified = TRUE),
    cycling_cells  = list(analyte = "proliferating tumor cells (shell)", units = "cells", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(analyte = "necrotic tumor cells (core entry)", units = "cells", specimen = "tumor", verified = TRUE),
    damaged_cells2 = list(analyte = "necrotic tumor cells (clearance delay 1)", units = "cells", specimen = "tumor", verified = TRUE),
    damaged_cells3 = list(analyte = "necrotic tumor cells (clearance delay 2)", units = "cells", specimen = "tumor", verified = TRUE),
    damaged_cells4 = list(analyte = "necrotic tumor cells (clearance delay 3)", units = "cells", specimen = "tumor", verified = TRUE),
    damaged_cells5 = list(analyte = "necrotic tumor cells (clearance delay 4)", units = "cells", specimen = "tumor", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 155L,
    n_studies      = 1L,
    disease_state  = paste(
      "Adults with ALK-positive advanced non-small cell lung cancer",
      "(NSCLC) enrolled in the global phase 2 lorlatinib study of",
      "Solomon 2018 (NCT01970865). Braniff 2025 Methods 'Dataset used",
      "for fitting' restricts the analysis to the subset of patients",
      "with at least two sum-of-longest-diameter (SLD) measurements and",
      "truncates follow-up at the 500th day of treatment, with all",
      "patients still on trial at that point right-censored. The",
      "virtual populations shown in Braniff 2025 Figures 4 and 5 are",
      "sized at 155 virtual patients to match this observed cohort."
    ),
    dose_range     = "Lorlatinib 100 mg orally once daily (the phase 2 dose of Solomon 2018).",
    endpoints      = paste(
      "Three scored outputs are matched by the Vpop selection algorithm",
      "(Braniff 2025 Methods 'Parameter sampling strategy'): baseline",
      "SLD (log-transformed before fitting the scoring distribution),",
      "best percentage change in SLD from baseline, and dropout time",
      "(time to progression or censoring, whichever came first). A",
      "two-component Gaussian mixture was fit jointly to these three",
      "endpoints and used as the acceptance score."
    ),
    notes          = paste(
      "Braniff 2025 is a virtual-population METHODOLOGY paper: the",
      "mechanistic model below is the engine, and the paper's",
      "contribution is a Metropolis-Hastings plausible-patient step",
      "extended to sample time-to-event endpoints, followed by the",
      "Allen 2016 acceptance-rejection virtual-population step. The",
      "paper reports Vpop sampling RANGES (Table S1 for the nine",
      "signaling / tumor-growth parameters, Table S3 for the two tumor",
      "physiology quantities) and NOT point estimates; ten thousand",
      "plausible patients were generated. Because the analysis fixes",
      "the PK covariates at their median values, the model carries no",
      "inter-individual variability on the PK parameters: all",
      "population heterogeneity is generated by sampling the eleven",
      "Vpop quantities and the necrotic initial-condition split",
      "phi_n0. Braniff 2025 Discussion states explicitly: 'in this work",
      "we do not incorporate pharmacokinetic variability, instead using",
      "a median model'."
    )
  )

  ini({
    # =====================================================================
    # NOTE ON PARAMETER PROVENANCE
    #
    # Braniff 2025 reports NO point estimates for the eleven quantities
    # that the virtual-population algorithm samples. Table S1 gives
    # minimum / maximum sampling bounds for nine signaling and tumor-growth
    # parameters; Table S3 gives bounds for the two tumor-physiology
    # quantities. Every such parameter below carries its published range in
    # the trailing comment, and the ini() value is a single documented
    # REFERENCE PLAUSIBLE PATIENT obtained by digitising the published
    # marginal distributions:
    #   * GFR, delta_shell, alpha and kg0 -- Braniff 2025 Figure 6a, which
    #     plots responder / non-responder quartiles on a "% value between
    #     feasible bounds" axis. The pooled median is approximated as
    #     0.65 * (responder median) + 0.35 * (non-responder median), the
    #     responder fraction being read from the Figure 4b waterfall (best
    #     percentage change crosses -30% at about the 35th percentile).
    #   * the remaining seven -- Braniff 2025 Figure S2, the pairwise and
    #     marginal plausible-parameter distributions, plotted on the same
    #     "% value between feasible bounds" axis.
    # Digitisation precision is about +/- 5 percentage points of each range.
    # The reference patient satisfies every published feasibility
    # constraint (see the model() header) and is a STRONG responder; the
    # population median best percentage change of about -45% (Figure 4b) is
    # reproduced in the vignette by SAMPLING the published ranges, not by
    # this single exemplar. Do not treat these ini() values as published
    # estimates.
    # =====================================================================

    # ---------------------------------------------------------------------
    # Lorlatinib PK. Braniff 2025 Supplement Section S1 states the model is
    # taken from Chen 2021 (reference [1] of the supplement) with all
    # covariates fixed at their median values. Braniff Table S2 tabulates
    # the resulting median values in day / L units; each is exactly 24x the
    # corresponding per-hour value in Chen 2021 Table 4, which confirms the
    # transcription (e.g. ka 3.113 1/h x 24 = 74.712 1/day).
    #
    # Braniff's clearance form CL(t) = CLmax * (1 - exp(-kc * t)) starts at
    # CL = 0 rather than at Chen's single-dose CLI = 9.035 L/h; Braniff
    # Table S2 tabulates no CLI. Encoded exactly as printed - see the
    # vignette "Assumptions and deviations" for the early-exposure
    # consequence. At steady state both forms give CL -> CLmax.
    # ---------------------------------------------------------------------
    lka          <- fixed(log(74.712))  ; label("First-order absorption rate constant ka (1/day)")                              # Braniff 2025 Table S2 ka = 74.712 1/day (= Chen 2021 Table 4 ka 3.113 1/h x 24)
    lcl_exp_inf  <- fixed(log(347.33))  ; label("Maximum auto-induced clearance CLmax (L/day)")                                 # Braniff 2025 Table S2 CLmax = 347.33 L/day (= Chen 2021 Table 4 CLMX 14.472 L/h x 24)
    lcl_exp_kdes <- fixed(log(0.48))    ; label("Clearance auto-induction rate constant kc (1/day)")                            # Braniff 2025 Table S2 kc = 0.48 1/day (= Chen 2021 Table 4 theta_IND 0.020 1/h x 24)
    lvc          <- fixed(log(120.51))  ; label("Central volume of distribution Vc (L)")                                        # Braniff 2025 Table S2 Vc = 120.51 L (= Chen 2021 Table 4 V2 120.511 L)
    lvp          <- fixed(log(154.91))  ; label("Peripheral volume of distribution Vp (L)")                                     # Braniff 2025 Table S2 Vp = 154.91 L (= Chen 2021 Table 4 V3 154.905 L)
    lq           <- fixed(log(528.05))  ; label("Inter-compartmental clearance Q (L/day)")                                      # Braniff 2025 Table S2 Q = 528.05 L/day (= Chen 2021 Table 4 Q 22.002 L/h x 24)
    lfdepot      <- fixed(log(0.759))   ; label("Oral bioavailability F (fraction, unitless)")                                  # Braniff 2025 Table S2 F = 0.759 (= Chen 2021 Table 4 theta_F 0.759)
    ld1          <- fixed(log(0.047833)); label("Zero-order absorption input duration D1 (day)")                                # NOT in Braniff Table S2; inherited from the cited PK source Chen 2021 Table 4 theta_D1 = 1.148 h = 0.047833 day, as extracted in inst/modeldb/specificDrugs/Chen_2021_lorlatinib.R. Braniff Supplement S1 describes "sequential zero and first order absorption" but tabulates no D1. See vignette Errata.

    # ---------------------------------------------------------------------
    # MAPK / PI3K signaling pathway. All five species are normalised so
    # that the untreated steady state is 1. Braniff 2025 Supplement S1.
    # ---------------------------------------------------------------------
    lkout        <- fixed(log(115.2))   ; label("Signal-transduction (de-phosphorylation) rate kout (1/day)")                   # Braniff 2025 Table S2 kout = 115.2 1/day, sourced there to Yamazaki 2014 (supplement reference [2])
    imax         <- fixed(1.0)          ; label("Maximal pALK inhibition by the ALK inhibitor Imax (fraction, unitless)")       # Braniff 2025 Table S2 Imax = 1.0, "Fitted*" (prior model fitting with preclinical data)
    lgfr_act     <- fixed(log(3.4149))  ; label("Basal growth-factor-receptor activation of RAS (unitless)")                    # Braniff 2025 Table S1 range 0.018 to 20; reference value = 17% of the range, digitised from Figure 6a (responder median ~12%, non-responder median ~27%)
    lic50_drug   <- fixed(log(16.870))  ; label("ALK inhibitor concentration for half-maximal pALK inhibition IC50 (nM)")       # Braniff 2025 Table S1 range 1.5 to 32.24 nM; reference value = 50% of the range (Figure S2 marginal is approximately uniform)
    lhill_drug   <- fixed(log(0.73))    ; label("Hill coefficient for ALK inhibitor effect on pALK (unitless)")                 # Braniff 2025 Table S1 range 0.70 to 0.76; reference value = 50% of the range (Figure S2 marginal is approximately uniform). The range is narrow, so the reference value is close to exact.

    # ---------------------------------------------------------------------
    # Anti-proliferative (pS6-driven) and pro-apoptotic (pAKT-driven)
    # effects. Braniff 2025 Supplement S1, the f_prolif / f_apop / k_kill
    # block. Note that both Hill exponents are named with the canonical
    # `hill` stem plus the effect they shape: the paper writes them as
    # alpha (proliferation) and beta (apoptosis).
    # ---------------------------------------------------------------------
    lhill_prolif <- fixed(log(3.4558))  ; label("Hill coefficient alpha for pS6 effect on proliferation (unitless)")            # Braniff 2025 Table S1 alpha range 0.28 to 5.02; reference value = 67% of the range, digitised from Figure 6a (responder median ~70%, non-responder median ~62%)
    lkc50_prolif <- fixed(log(21.936))  ; label("pS6 level giving half-maximal effect on proliferation KC50,prolif (unitless)") # Braniff 2025 Table S1 range 0.00025 to 36.56; reference value = 60% of the range, digitised from the Figure S2 marginal (rises monotonically toward the upper bound)
    lhill_apop   <- fixed(log(3.9486))  ; label("Hill coefficient beta for pAKT effect on apoptosis (unitless)")                # Braniff 2025 Table S1 beta range 0.21 to 6.24; reference value = 62% of the range, digitised from the Figure S2 marginal (rises monotonically toward the upper bound)
    lkc50_apop   <- fixed(log(19.364))  ; label("pAKT level giving half-maximal effect on apoptosis KC50,apop (unitless)")      # Braniff 2025 Table S1 range 0.023 to 33.37; reference value = 58% of the range, digitised from the Figure S2 marginal
    lkkill_scale <- fixed(log(0.135))   ; label("In vitro to in vivo cell-killing scaling factor theta_kkill (unitless)")       # Braniff 2025 Table S2 theta_kkill = 0.135, "Fitted*" (prior model fitting with preclinical data)

    # ---------------------------------------------------------------------
    # Shell-and-core tumor growth. Braniff 2025 Supplement S1 and Figure 1.
    # ---------------------------------------------------------------------
    lkg0            <- fixed(log(0.646))   ; label("Basal net growth rate of the proliferating shell kg0 (1/day)")              # Braniff 2025 Table S1 range 0.1 to 1.5 1/day; reference value = 39% of the range, digitised from Figure 6a (responder median ~41%, non-responder median ~36%)
    ltau_necrotic   <- fixed(log(0.212))   ; label("Per-compartment transit time of the necrotic clearance chain tau (day)")    # Braniff 2025 Table S1 range 0.04 to 0.9 day; reference value = 20% of the range, digitised from the Figure S2 marginal (sharply right-skewed, peaked near the lower bound). Figure 1 labels the Nn1->...->Nn4->0 transitions with rate 1/tau.
    lrho_core       <- fixed(log(0.582))   ; label("Core-shell ratio rho_core = sum(necrotic) / proliferating (unitless)")      # Braniff 2025 Table S3 range 0.15 to 0.75; reference value = 72% of the range, digitised from the Figure S2 marginal (strongly left-skewed, peaked near the upper bound)
    ldelta_shell    <- fixed(log(3.652))   ; label("Proliferating shell thickness delta_shell (mm)")                            # Braniff 2025 Table S3 range 0.07 to 1.3 cm = 0.7 to 13 mm; reference value = 24% of the range, digitised from Figure 6a (responder median ~27%, non-responder median ~19%). Expressed in mm here so it shares units with the SLD output.
    phi_n0          <- fixed(0.5294)       ; label("Fraction of the necrotic mass initially in the core-entry compartment (unitless)") # Braniff 2025 Supplement S2 randomises phi_n0 uniformly on [phi_n0_min, 1] where phi_n0_min = 1 - 4 * kg0 * tau * Np0 / NnTot; the reference value is the MEAN of that uniform distribution for this parameter set (phi_n0_min = 0.0587). Not a published point estimate.

    # ---------------------------------------------------------------------
    # Fixed structural constants.
    # ---------------------------------------------------------------------
    eps_bifurc  <- fixed(0.015)   ; label("Offset of the necrotic fraction below its stasis value (unitless)")                  # Braniff 2025 Supplement S2: "We thus require that (1 - epsilon) * phi_n^SS = phi_n to ensure a positive growth rate; here we use epsilon = 0.015."
    cell_vol    <- fixed(1e-5)    ; label("Volume of a single tumor cell (mm^3)")                                               # Braniff 2025 Supplement S1: "can be converted to a volume by assuming 1e5 cells per uL", i.e. 1 cell = 1e-5 uL = 1e-5 mm^3
    mw_drug     <- fixed(406.41)  ; label("Molar mass of lorlatinib (g/mol)")                                                   # Physical constant for lorlatinib (C21H19FN6O2); required because Braniff Table S1 reports IC50 in nM while the PK model returns mg/L. Not a Braniff parameter - see vignette Errata.

    # ---------------------------------------------------------------------
    # Residual error. Braniff 2025 fits no residual-error model: the QSP
    # model is deterministic and ALL population variability is generated by
    # the Vpop parameter sampling. Both endpoints therefore carry a fixed
    # zero residual SD, per the convention for unreported RUV. Simulate
    # with rxode2::rxSolve(mod, ..., omega = NA) for the deterministic
    # typical-value trajectories the paper works with.
    # ---------------------------------------------------------------------
    propSd     <- fixed(0) ; label("Proportional residual SD for plasma lorlatinib (fraction) - not reported by Braniff 2025")  # Braniff 2025 reports no residual-error model; fixed at 0 rather than invented
    addSd_sld  <- fixed(0) ; label("Additive residual SD for the sum of longest diameters (mm) - not reported by Braniff 2025") # Braniff 2025 reports no residual-error model; fixed at 0 rather than invented
  })

  model({
    # =====================================================================
    # 0. Back-transform.
    # =====================================================================
    ka           <- exp(lka)
    cl_exp_inf   <- exp(lcl_exp_inf)
    cl_exp_kdes  <- exp(lcl_exp_kdes)
    vc           <- exp(lvc)
    vp           <- exp(lvp)
    q            <- exp(lq)
    fdepot       <- exp(lfdepot)
    d1           <- exp(ld1)
    kout         <- exp(lkout)
    gfr_act      <- exp(lgfr_act)
    ic50_drug    <- exp(lic50_drug)
    hill_drug    <- exp(lhill_drug)
    hill_prolif  <- exp(lhill_prolif)
    kc50_prolif  <- exp(lkc50_prolif)
    hill_apop    <- exp(lhill_apop)
    kc50_apop    <- exp(lkc50_apop)
    kkill_scale  <- exp(lkkill_scale)
    kg0          <- exp(lkg0)
    tau_necrotic <- exp(ltau_necrotic)
    rho_core     <- exp(lrho_core)
    delta_shell  <- exp(ldelta_shell)

    # =====================================================================
    # 1. Derived tumor-geometry and initial-condition quantities.
    #
    #    Braniff 2025 Supplement S2 derives the initial state of the tumor
    #    from the two sampled physiology quantities rho_core and
    #    delta_shell. Everything in this block is algebra from that
    #    section, evaluated once from the parameters so that the model is
    #    self-initialising for any sampled parameter vector (this is what
    #    lets the vignette simulate a virtual population in one rxSolve
    #    call rather than looping subject by subject).
    #
    #    rad_const converts a cell COUNT to the radius of the equivalent
    #    sphere: V = n * cell_vol and r = (3V / 4pi)^(1/3), so
    #    r = rad_const * n^(1/3) with rad_const = (3 * cell_vol / 4pi)^(1/3).
    #    Braniff's printed expression for Np omits this dimensional
    #    constant (it equates a length cubed with a cell count); the
    #    constant is recovered from the supplement's own statement that
    #    there are 1e5 cells per microlitre. See vignette Errata.
    # =====================================================================
    rad_const <- (3 * cell_vol / (4 * pi)) ^ (1 / 3)

    # Proliferating-cell count implied by the shell thickness and the
    # core-shell ratio (Braniff Supplement S2, "Initial condition
    # randomization (a)"): delta_shell = r_total - r_core.
    np0    <- (delta_shell / rad_const) ^ 3 /
              ((1 + rho_core) ^ (1 / 3) - rho_core ^ (1 / 3)) ^ 3
    nn_tot <- rho_core * np0

    # Necrotic FRACTION phi_n and its stasis value. Braniff Supplement S2:
    # phi_n = rho_core / (1 + rho_core), and ke is assigned so that
    # (1 - eps) * phi_n^SS = phi_n, which places the tumor just on the
    # growing side of the growth/clearance bifurcation.
    phi_n     <- rho_core / (1 + rho_core)
    phi_n_ss  <- phi_n / (1 - eps_bifurc)
    ke_necrotic <- 1 / ((1 / kg0) * (phi_n_ss / (1 - phi_n_ss)) -
                        4 * tau_necrotic)

    # Necrotic clearance-chain transit rate (Braniff Figure 1 labels the
    # Nn1 -> ... -> Nn4 -> 0 transitions 1/tau).
    ktr <- 1 / tau_necrotic

    # =====================================================================
    # 2. Lorlatinib PK (Braniff Supplement S1, first equation block).
    #    Two compartments in CONCENTRATION form in the source; encoded here
    #    in the equivalent amount form so that rxode2 dosing works
    #    naturally. CL(t) rises from 0 to CLmax with rate kc.
    # =====================================================================
    cl <- cl_exp_inf * (1 - exp(-cl_exp_kdes * time))

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl / vc) * central -
                          (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    dur(depot) <- d1
    f(depot)   <- fdepot

    # Plasma concentration in mg/L, and the molar concentration the
    # signaling model needs (Braniff Table S1 reports IC50 in nM).
    conc_mgl <- central / vc
    conc_nm  <- conc_mgl * 1e6 / mw_drug

    # =====================================================================
    # 3. MAPK / PI3K signaling (Braniff Supplement S1, second block).
    #    Each species relaxes toward its input with the common rate kout.
    #    All are normalised so that the untreated steady state is 1.
    # =====================================================================
    inh <- imax * conc_nm ^ hill_drug /
           (ic50_drug ^ hill_drug + conc_nm ^ hill_drug)

    d/dt(palk) <- kout * ((1 - inh) - palk)
    d/dt(ras)  <- kout * ((palk + gfr_act) / (1 + gfr_act) - ras)
    d/dt(perk) <- kout * (ras - perk)
    d/dt(pakt) <- kout * (ras - pakt)
    d/dt(ps6)  <- kout * (perk - ps6)

    palk(0) <- 1
    ras(0)  <- 1
    perk(0) <- 1
    pakt(0) <- 1
    ps6(0)  <- 1

    # =====================================================================
    # 4. Anti-proliferative and pro-apoptotic drug effect (Braniff
    #    Supplement S1, the f_prolif / f_apop / k_kill equations).
    #
    #    ps6_pos and akt_comp are numerical guards only. Both pS6 and
    #    (1 - pAKT) are non-negative for every admissible trajectory
    #    (pAKT relaxes toward RAS <= 1 from pAKT(0) = 1), but a
    #    round-off-sized negative value raised to a non-integer Hill
    #    exponent returns NaN and aborts the solve. Clamping at zero
    #    changes no admissible trajectory.
    # =====================================================================
    ps6_pos  <- max(ps6, 0)
    akt_comp <- max(1 - pakt, 0)

    f_prolif <- (kc50_prolif ^ hill_prolif + 1) * ps6_pos ^ hill_prolif /
                (kc50_prolif ^ hill_prolif + ps6_pos ^ hill_prolif)
    f_apop   <- akt_comp ^ hill_apop /
                (kc50_apop ^ hill_apop + akt_comp ^ hill_apop)
    kkill    <- kkill_scale * (1 - (f_prolif - f_apop))

    # =====================================================================
    # 5. Shell-and-core tumor growth (Braniff Supplement S1, third block).
    #
    #    kn is the shell -> core transfer rate that holds the core-shell
    #    ratio constant. Braniff prints
    #      kn = [kt * Nn4 + rho_core * (kg0 * (1 - kkill) * Np - kt * Nn4)] / Np
    #    but with rho_core read as the core-shell RATIO that expression
    #    does NOT hold the ratio constant, and the resulting model cannot
    #    reach the untreated doubling-time window the paper's own Table S3
    #    enforces (the fastest untreated doubling time attainable anywhere
    #    in the Table S1 / S3 parameter space is about 332 days, against a
    #    required 26 to 165 days), nor the tumor shrinkage of Figure 4b.
    #    Reading rho_core in this one equation as the necrotic FRACTION
    #    phi_n = rho_core / (1 + rho_core) - a quantity the supplement
    #    itself defines and uses throughout Section S2 - makes the printed
    #    expression ALGEBRAICALLY IDENTICAL to the ratio-preserving rate
    #      kn = [rho_core * kg0 * (1 - kkill) * Np + kt * Nn4] / ((1 + rho_core) * Np)
    #    because 1 - phi_n = 1 / (1 + rho_core). That reading is adopted
    #    here; it satisfies Table S3 and reproduces Figure 4b. See the
    #    vignette "Assumptions and deviations" for the full falsification.
    # =====================================================================
    grow  <- kg0 * (1 - kkill) * cycling_cells
    clear <- ktr * damaged_cells5
    kn    <- (rho_core * grow + clear) / ((1 + rho_core) * cycling_cells)

    d/dt(cycling_cells)  <- (1 - kkill) * kg0 * cycling_cells -
                            kn * cycling_cells
    d/dt(damaged_cells1) <- kn * cycling_cells - ke_necrotic * damaged_cells1
    d/dt(damaged_cells2) <- ke_necrotic * damaged_cells1 - ktr * damaged_cells2
    d/dt(damaged_cells3) <- ktr * (damaged_cells2 - damaged_cells3)
    d/dt(damaged_cells4) <- ktr * (damaged_cells3 - damaged_cells4)
    d/dt(damaged_cells5) <- ktr * (damaged_cells4 - damaged_cells5)

    # Initial conditions (Braniff Supplement S2). The total necrotic mass
    # rho_core * Np0 is split between the core-entry compartment (fraction
    # phi_n0) and the four clearance-delay compartments (equal shares).
    cycling_cells(0)  <- np0
    damaged_cells1(0) <- phi_n0 * nn_tot
    damaged_cells2(0) <- (1 - phi_n0) * nn_tot / 4
    damaged_cells3(0) <- (1 - phi_n0) * nn_tot / 4
    damaged_cells4(0) <- (1 - phi_n0) * nn_tot / 4
    damaged_cells5(0) <- (1 - phi_n0) * nn_tot / 4

    # =====================================================================
    # 6. Observations.
    #    sld is the diameter of the combined shell-and-core sphere
    #    (Braniff Figure 1 caption: "The sum of longest diameter (SLD) for
    #    each simulated patient is defined as the diameter of the combined
    #    shell-and-core mass"), reported in mm to match the clinical
    #    endpoint axis of Figure 4a. Cc is plasma lorlatinib in ng/mL.
    # =====================================================================
    n_tot <- cycling_cells + damaged_cells1 + damaged_cells2 +
             damaged_cells3 + damaged_cells4 + damaged_cells5
    sld   <- 2 * rad_const * n_tot ^ (1 / 3)
    Cc    <- conc_mgl * 1000

    Cc  ~ prop(propSd)
    sld ~ add(addSd_sld)
  })
}
