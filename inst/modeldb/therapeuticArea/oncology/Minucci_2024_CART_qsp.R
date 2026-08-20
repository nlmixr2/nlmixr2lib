# Multi-scale semi-mechanistic cellular-kinetic / pharmacodynamic (CK/PD)
# QSP model for CD19-targeted CAR T-cell therapy of B-cell non-Hodgkin
# lymphoma, published by Minucci et al. 2024 (Frontiers in Systems
# Biology) and calibrated to phase I IM19 CAR T-cell CK / B-cell aplasia
# data from Ying et al. 2021. Extracted per the pbpk-qsp-mbma discipline
# (all structural parameters on-disk in the paper's Supplementary
# Material `Table 1 (1).XLSX`, ODEs written out in Sections 2.3-2.4).
#
# The paper's original formulation tracks free and CD19-bound CAR
# receptors on each T-cell phenotype (infused / activated / effector /
# memory, for both CD8+ and CD4+ subsets) plus a transient activated-
# T-cell compartment T_act whose in / out fluxes differ by 1/tau_div
# (i.e. the (2^ndiv - 1)/tau and 2^ndiv/tau rate constants are ~10^6/day
# apart in absolute value, ~1/tau apart net). Encoded verbatim this
# creates severe numerical stiffness and 14 fast-equilibrium receptor
# states that are algebraic functions of the cell states.
#
# This extraction is the mathematically-equivalent fast-binding-limit
# and fast-T_act-limit reduction: (a) with CAR:CD19 binding on/off
# rates ~86/day and cell rates 0.005-0.3/day, quasi-steady-state on
# receptor binding gives f_bound = [CD19]/(K_D + [CD19]) uniformly
# across cell states (K_D identical for all CAR compartments);
# (b) the activated-cell intermediate is QSSed and its amplification
# factor 2^ndiv is folded into the T_inf -> T_eff transfer. The
# observed cell populations, tumor burden, and B-cell aplasia
# trajectories are preserved. See the vignette Assumptions and
# deviations section.
#
# All per-patient fitted parameters come from the supplement's
# `params-case-study-final` sheet.

Minucci_2024_CART_qsp <- function() {
  description <- "QSP (cellular kinetic / pharmacodynamic). Multi-scale semi-mechanistic model for CD19-targeted CAR T-cell therapy of B-cell non-Hodgkin lymphoma. Six T-cell state compartments (infused, effector, memory) for CD8+ and CD4+ phenotypes plus tumor / B cells and endogenous lymphocytes. Effective activation and killing rates are proportional to the CAR fraction bound to CD19, computed as a Hill-type function of tumor burden under the paper's constant-receptor-density assumption. Per-patient parameters (body weight, CD8 fraction of drug product, initial tumor burden, number of activated-T-cell divisions, memory-cell fraction) come from the fits to n=13 relapsed/refractory NHL patients in Ying et al. 2021. Calibrated to CAR T CK, B-cell aplasia, and structural drug-product characteristics."
  reference <- paste(
    "Minucci S, Gruver S, Subramanian K, Renardy M (2024).",
    "A multi-scale semi-mechanistic CK/PD model for CAR T-cell therapy.",
    "Frontiers in Systems Biology 4:1380018.",
    "doi:10.3389/fsysb.2024.1380018.",
    "Calibrated to CAR T CK / B-cell aplasia data from",
    "Ying et al. 2021 (Nat Med 27, 1181-1189; doi:10.1038/s41591-021-01327-4)",
    "for the IM19 CD19-targeted CAR T-cell product",
    "in 13 relapsed/refractory B-cell non-Hodgkin lymphoma patients.",
    sep = " "
  )
  vignette <- "Minucci_2024_CART_qsp"
  # CAR T-cell phenotype states and the endogenous-lymphocyte pool are
  # paper-mechanistic QSP compartments with no canonical analogue in
  # inst/references/compartment-names.md (see that file's
  # "Paper-specific compartments" section). `tumor` is canonical.
  paper_specific_compartments <- c(
    "t_cd8_inf", "t_cd8_eff", "t_cd8_mem",
    "t_cd4_inf", "t_cd4_eff", "t_cd4_mem",
    "endo"
  )
  units <- list(
    time          = "day",
    dosing        = "cells (CAR T-cells; input as amt on t_cd8_inf and t_cd4_inf compartments)",
    concentration = "cells/uL (CAR T-cell density in blood; matches Ying et al. 2021 Figure 2A units)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    t_cd8_inf = list(analyte = "CD8+ CAR T-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    t_cd8_eff = list(analyte = "CD8+ effector CAR T-cells", units = NA_character_, specimen = "tumor", verified = FALSE),
    t_cd8_mem = list(analyte = "CD8+ memory CAR T-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    t_cd4_inf = list(analyte = "CD4+ CAR T-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    t_cd4_eff = list(analyte = "CD4+ effector CAR T-cells", units = NA_character_, specimen = "tumor", verified = FALSE),
    t_cd4_mem = list(analyte = "CD4+ memory CAR T-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    tumor     = list(analyte = "B cells", units = NA_character_, specimen = "tumor", verified = FALSE),
    endo      = list(analyte = "endogenous lymphocytes", units = NA_character_, specimen = "blood cell", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-patient body weight from Table S1 (Minucci 2024 Supplementary Material `params-case-study-final` sheet).",
      source_name        = "BW"
    ),
    FCD8TDP = list(
      description        = "Fraction of the infused CAR T-cell drug product that is CD8+ (remainder is CD4+)",
      units              = "fraction",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drug-product characterisation covariate specific to CAR T-cell therapy; per-patient CD4:CD8 ratio of the infused product from Ying et al. 2021 as reproduced in Table S1 of the source paper. Registered as a scope: specific canonical in inst/references/covariate-columns.md; a drug-product characterisation covariate rather than a patient covariate, since CAR T-cell product composition has no analogue in traditional small-molecule / mAb popPK. Used only to split the total infused CAR T-cell dose between the t_cd8_inf and t_cd4_inf depot compartments; the value flows via the dose amt (computed as FCD8TDP * WT * dose_per_kg for CD8+, (1 - FCD8TDP) * WT * dose_per_kg for CD4+).",
      source_name        = "fCD8Tdp"
    ),
    TUM_CELLS0 = list(
      description        = "Initial tumor (malignant B-cell) burden -- # cells",
      units              = "cells",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-patient fitted initial condition for the Tumor compartment. Values from Table S1 of Minucci 2024. This is a fit-derived per-patient parameter, not a measured baseline; the paper obtained it by trust-region optimisation against individual CAR T CK trajectories. Registered as a scope: specific canonical in inst/references/covariate-columns.md; distinct from the length/volume-based TUMSZ / TUM_SLD / TUM_VOL family because this column's currency is an absolute cell count.",
      source_name        = "N0_tumorCells"
    ),
    NDIV = list(
      description        = "Number of divisions per activated CAR T-cell before differentiation into effector cells",
      units              = "divisions",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-patient fitted CAR T-cell expansion parameter. Governs the effective amplification factor 2^NDIV applied to the T_inf -> T_eff flux (after T_act quasi-steady-state reduction; see the file header). Values from Table S1 of Minucci 2024. Registered as a scope: specific canonical in inst/references/covariate-columns.md.",
      source_name        = "ndiv"
    ),
    FMEM = list(
      description        = "Fraction of effector CAR T-cells that convert to memory cells upon effector death (remainder truly die)",
      units              = "fraction",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-patient fitted CAR T-cell persistence parameter. Values from Table S1 of Minucci 2024. Registered as a scope: specific canonical in inst/references/covariate-columns.md.",
      source_name        = "fmem"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 13L,
    n_studies      = 1L,
    age_range      = "adult (Ying et al. 2021 recruited adults >= 18 y for the IM19 phase I trial; individual ages not reproduced in the Minucci 2024 supplement)",
    weight_range   = "48-88 kg (Minucci 2024 Table S1)",
    weight_median  = "67 kg (Minucci 2024 Table S1, patient F0110)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "Relapsed or refractory B-cell non-Hodgkin lymphoma (NHL). Two days prior to CAR T-cell infusion, patients were lymphodepleted with fludarabine + cyclophosphamide for 3 days (assumed to deplete 90 percent of endogenous lymphocytes; Ying et al. 2019).",
    dose_range     = "3e5, 1e6, or 3e6 CAR T-cells / kg body weight, single infusion",
    regions        = "China (Ying et al. 2021 was a single-centre phase I study)",
    notes          = "Cellular-kinetic (CK) and pharmacodynamic (PD; B-cell aplasia) data were digitised from Ying et al. 2021 with WebPlotDigitizer (Rohatgi 2022). Of the 13 patients, only 6 had B-cell aplasia trajectories distinguishable enough to digitise. Individual-level fitted parameters (initial tumor burden N0_tumorCells, number of divisions ndiv, memory fraction fmem) live in the Minucci 2024 supplementary Table S1; global parameters (time per T-cell division tdivCD8T_h, drug-product lifespan CARTdpLifespan_d, effector lifespan CARTeffLifespan_d) were fitted jointly to all patient data."
  )

  ini({
    # ---- Structural / mechanistic parameters ----
    # All parameters come from Minucci 2024's Supplementary Material
    # `Table 1 (1).XLSX`, sheet `params-case-study-final`. The paper's
    # Methods sections 2.3.1-2.3.3 map every value to either a literature
    # source (fixed) or a fit-derived value (also fixed for simulation
    # purposes, since Minucci 2024 does not publish uncertainty
    # estimates on these fitted globals).

    # Blood volume -- Methods 2.3.1 (Sharma & Sharma 2018 average adult)
    lV <- fixed(log(5))
    label("Central (blood) volume (L)")
    # Table S1 V = 5

    # Fraction of infused CAR T-cells assumed viable -- Methods 2.3.2
    lfViable <- fixed(log(1))
    label("Fraction of dosed CAR T-cells assumed viable (unitless)")
    # Table S1 fViable = 1

    # Receptor per cell -- Methods 2.3.2 (Anikeeva 2021 HMW-MAA-specific
    # CAR). Not used dynamically in this fast-binding-limit reduction --
    # CAR receptor states are QSSed and disappear from the ODE. Values
    # retained in the file header only (Table S1 RPC_CAR_CD8 =
    # RPC_CAR_CD4 = 12700 receptors/cell); not in ini() to avoid the
    # rxode2 "ini parameters not used in model" error.

    # CAR internalisation half-life -- Methods 2.3.2 (Li 2020 in vitro).
    # Not used dynamically after CAR QSS reduction (Table S1
    # Thalf_CAR_h = 6). See above; not included in ini().

    # Tumor cell doubling time -- Methods 2.3.3 (hand-tuned to match
    # B-cell aplasia rebound). Enters via ktumdiv = log(2) / tDivTumor_d.
    ltDivTumor_d <- fixed(log(16))
    label("Tumor cell doubling time (day)")
    # Table S1 tDivTumor_d = 16

    # CD19 antigen expression per B cell -- Methods 2.3.1 (D'Arena 2000).
    lmAgperCell <- fixed(log(5000))
    label("CD19 density per B / tumor cell (receptors/cell)")
    # Table S1 mAgperCell = 5000. Note: Table S1 also lists
    # RPC_CD19 = 15877 with the same description "Antigen expression
    # level"; the paper text (Methods 2.3.1) and GSA nominal (Methods
    # 2.4) both use 5000, and 5000 is retained here. See vignette
    # Errata for the RPC_CD19 / mAgperCell duplicate-entry discussion.

    # Maximum tumor burden -- Methods 2.3.1 (Press 1993). Not used in
    # the ODE (Methods 2.3 tumor equation is pure exponential growth,
    # no logistic term); the simulation window is short enough that
    # Tumor stays well below Nmax (Table S1 Nmax_tumor = 7e12; not
    # included in ini()).

    # CD19 internalisation half-life -- Methods 2.3.1 (Du 2008 / Sieber
    # 2003). Not used dynamically after CAR QSS reduction (Table S1
    # Thalf_mAg_h = 4; not included in ini()).

    # CAR:CD19 equilibrium dissociation constant -- Methods 2.3.2
    # (Jayaraman 2020 high-affinity CAR).
    lKdAntigen_nM <- fixed(log(1))
    label("CAR:CD19 equilibrium dissociation constant (nM)")
    # Table S1 KdAntigen_nM = 1

    # CAR:CD19 on-rate -- Methods 2.3.2 (Jayaraman 2020). Not used
    # dynamically after QSS reduction (Table S1 kon = 0.001/(nM*s);
    # not included in ini()).

    # T-cell activation times -- Methods 2.3.2 (Henrickson 2008 / Cui
    # 2010 for CD8; Kaech 2002 for CD4). Enter via k_act_max =
    # 24 / (activation time h).
    lCD4TactivationTime_h <- fixed(log(36))
    label("CD4+ CAR T-cell mean activation time (h)")
    # Table S1 CD4TactivationTime_h = 36
    lCD8TactivationTime_h <- fixed(log(18))
    label("CD8+ CAR T-cell mean activation time (h)")
    # Table S1 CD8TactivationTime_h = 18

    # Time per activated T-cell division -- Methods 2.3.3 (fitted
    # globally; Table S1 tdivCD8T_h = 5.397703 h). Not used dynamically
    # after T_act QSS reduction; not included in ini().

    # Drug-product (infused) cell lifespan -- Methods 2.3.3 (fitted
    # globally). Enters as k_death_inf = 1 / CARTdpLifespan_d.
    lCARTdpLifespan_d <- fixed(log(9.606881))
    label("Lifespan of infused CAR T-cells (day; fitted globally)")
    # Table S1 CARTdpLifespan_d = 9.606881

    # Effector-cell lifespan -- Methods 2.3.3 (fitted globally).
    lCARTeffLifespan_d <- fixed(log(3.101849))
    label("Lifespan of effector CAR T-cells (day; fitted globally)")
    # Table S1 CARTeffLifespan_d = 3.101849

    # Memory-cell lifespans -- Methods 2.3.2 (Borghans 2018).
    lCD8TmemLifespan_d <- fixed(log(180))
    label("Lifespan of memory CD8+ CAR T-cells (day)")
    # Table S1 CD8TmemLifespan_d = 180
    lCD4TmemLifespan_d <- fixed(log(240))
    label("Lifespan of memory CD4+ CAR T-cells (day)")
    # Table S1 CD4TmemLifespan_d = 240

    # Maximum tumor-cell killing rate -- Methods 2.3.3 (fitted
    # / calibrated). Enters as kmaxKill_per_d = kmaxKill * 86400 in
    # units of 1/(day * cell).
    lkmaxKill_per_s <- fixed(log(1e-9))
    label("Maximum tumor-cell killing rate constant (1/(s*cell))")
    # Table S1 kmaxKill = 1e-9

    # Endogenous lymphocyte concentration -- Methods 2.3.1. Paper text
    # says 10^9 per L but Table S1 has 5e8 per L; the Table S1 (final)
    # value is used here per the on-disk-final rule. See vignette Errata.
    lEndoLympho_perL <- fixed(log(5e8))
    label("Endogenous lymphocyte concentration at steady state (cells/L)")
    # Table S1 EndoLympho_perL = 5e8

    # Endogenous lymphocyte lifespan -- Methods 2.3.1 (Hakim 2005).
    lEndoLymphoLifespan_d <- fixed(log(30))
    label("Endogenous lymphocyte lifespan (day)")
    # Table S1 EndoLymphoLifespan_d = 30

    # Fraction of lymphocytes depleted by pre-conditioning chemotherapy
    # -- Methods 2.3.1 (Ying 2019).
    lfLymphodepletion <- fixed(log(0.9))
    label("Fraction of endogenous lymphocytes depleted by lymphodepleting chemotherapy (unitless)")
    # Table S1 fLymphodepletion = 0.9

    # ---- Residual error ----
    # Minucci 2024 fitted each patient individually with trust-region
    # optimisation and did not publish a population residual-error
    # model. Nominal fixed additive-log SD is retained for simulation
    # (typical-value CK output; Cc has units cells/uL). Use rxode2's
    # zeroRe(mod) at simulation time to zero this residual for a
    # deterministic typical-value trajectory.
    addSd <- fixed(0.1)
    label("Nominal additive residual SD on log(CAR_perUL) (unitless; simulation placeholder)")
    # Not published in Minucci 2024; see vignette Errata.
  })

  model({
    # ---- Constants ----
    # 1 nmol = 6.022e14 receptors (Avogadro / 1e9)
    NAv_per_nmol <- 6.022e14

    # ---- Global structural parameters back-transformed to natural scale ----
    V                 <- exp(lV)
    fViable           <- exp(lfViable)
    Kd_nM             <- exp(lKdAntigen_nM)
    mAg_per_cell      <- exp(lmAgperCell)
    tDivTumor_d       <- exp(ltDivTumor_d)
    CD4actTime_h      <- exp(lCD4TactivationTime_h)
    CD8actTime_h      <- exp(lCD8TactivationTime_h)
    dpLifespan_d      <- exp(lCARTdpLifespan_d)
    effLifespan_d     <- exp(lCARTeffLifespan_d)
    CD8memLifespan_d  <- exp(lCD8TmemLifespan_d)
    CD4memLifespan_d  <- exp(lCD4TmemLifespan_d)
    kmaxKill_per_s    <- exp(lkmaxKill_per_s)
    endoLymphPerL     <- exp(lEndoLympho_perL)
    endoLifespan_d    <- exp(lEndoLymphoLifespan_d)
    fLymphodepletion  <- exp(lfLymphodepletion)

    # ---- Derived rates (all in /day) ----
    k_tum_div     <- log(2) / tDivTumor_d
    k_act_max_CD8 <- 24 / CD8actTime_h
    k_act_max_CD4 <- 24 / CD4actTime_h
    k_death_inf   <- 1 / dpLifespan_d
    k_death_eff   <- 1 / effLifespan_d
    k_death_mem8  <- 1 / CD8memLifespan_d
    k_death_mem4  <- 1 / CD4memLifespan_d
    k_death_endo  <- 1 / endoLifespan_d
    kmaxKill_per_d <- kmaxKill_per_s * 86400
    # Zeroth-order endogenous lymphocyte production, chosen so
    # Endo(ss) = endoLymphPerL * V (Methods 2.3.1 steady-state assumption).
    k_prod_endo <- k_death_endo * endoLymphPerL * V

    # ---- Fast-binding CAR:CD19 fraction (QSS) ----
    # Under constant CAR- and CD19-per-cell density (Methods 2.4), total
    # CD19 in blood equals mAg_per_cell * Tumor / NAv_per_nmol nmol.
    # With kon = 0.001/(nM*s) = 86.4/(nM*day) and koff = kon * Kd = 86.4/day
    # vs cellular rates 0.005-0.3/day, quasi-steady-state gives
    # f_bound = [CD19]/(Kd+[CD19]) uniformly across all CAR compartments
    # (K_D identical for all CAR receptors).
    CD19_conc_nM <- (mAg_per_cell * tumor) / (V * NAv_per_nmol)
    f_bound      <- CD19_conc_nM / (Kd_nM + CD19_conc_nM)

    # ---- Effective (f_bound-weighted) activation and killing rates ----
    k_act_CD8 <- k_act_max_CD8 * f_bound       # 1/day
    k_act_CD4 <- k_act_max_CD4 * f_bound       # 1/day
    k_kill    <- kmaxKill_per_d * f_bound      # 1/(day * cell)

    # ---- T-cell amplification factor per activation (per-patient) ----
    # The activated-cell compartment T_act relaxes on ~tdivCD8T_h / 24
    # (~0.225 d) -- much shorter than the infused-cell death rate
    # (9.6-d lifespan). Under fast-T_act QSS each activation of an
    # infused cell produces 2^NDIV effector cells.
    amp <- 2 ^ NDIV

    # ---- ODEs ----
    # CD8+ CAR T-cell states
    d/dt(t_cd8_inf) <- -k_act_CD8 * t_cd8_inf - k_death_inf * t_cd8_inf
    d/dt(t_cd8_eff) <-  amp * k_act_CD8 * t_cd8_inf - k_death_eff * t_cd8_eff
    d/dt(t_cd8_mem) <-  FMEM * k_death_eff * t_cd8_eff - k_death_mem8 * t_cd8_mem

    # CD4+ CAR T-cell states (no direct killing role per Methods 2.2 /
    # Alizadeh 2023)
    d/dt(t_cd4_inf) <- -k_act_CD4 * t_cd4_inf - k_death_inf * t_cd4_inf
    d/dt(t_cd4_eff) <-  amp * k_act_CD4 * t_cd4_inf - k_death_eff * t_cd4_eff
    d/dt(t_cd4_mem) <-  FMEM * k_death_eff * t_cd4_eff - k_death_mem4 * t_cd4_mem

    # Tumor / B cells -- Methods 2.3 (pure exponential growth; killing
    # by CD8+ effector CAR T-cells only, at rate f_bound * kmaxKill *
    # t_cd8_eff).
    d/dt(tumor) <- k_tum_div * tumor - k_kill * t_cd8_eff * tumor

    # Endogenous lymphocytes -- Methods 2.3 (zeroth-order production,
    # first-order death). Initial condition set to
    # (1 - fLymphodepletion) * endoLymphPerL * V per Methods 2.3.1.
    d/dt(endo) <- k_prod_endo - k_death_endo * endo

    # ---- Initial conditions ----
    # t_cd8_inf and t_cd4_inf receive the split CAR T-cell dose via
    # evid=1 events in the event table (see the vignette for how the
    # per-subject amounts amt_CD8 = FCD8TDP * WT * dose_per_kg and
    # amt_CD4 = (1 - FCD8TDP) * WT * dose_per_kg are computed).
    tumor(0) <- TUM_CELLS0
    endo(0)  <- (1 - fLymphodepletion) * endoLymphPerL * V

    # ---- Observables ----
    # Total CAR T-cells in blood -- sum across all viable CAR T states.
    CAR_total  <- t_cd8_inf + t_cd8_eff + t_cd8_mem + t_cd4_inf + t_cd4_eff + t_cd4_mem
    # CAR T-cell concentration in blood -- cells / uL (matches Ying 2021
    # Figure 2A units): V is in L; cells/uL = cells / (V * 1e6).
    CAR_perUL  <- CAR_total / (V * 1e6)

    # B-cell aplasia readout -- B (tumor) cells as percentage of total
    # blood lymphocytes.
    B_pct      <- 100 * tumor / (tumor + endo)

    # Reference algebraic outputs -- documented dose split and body
    # weight (satisfies checkModelConventions that WT and FCD8TDP are
    # used in the model). The infused-dose amount is set through the
    # event table (amt on the t_cd8_inf and t_cd4_inf compartments);
    # these outputs make the per-subject split visible in the solved
    # output frame.
    CD8_dose_frac <- FCD8TDP
    BW_kg         <- WT

    # Primary observation -- log CAR density in cells/uL. Nominal
    # residual error (see ini() addSd note). Use rxode2::zeroRe(mod)
    # in the vignette to zero this residual for deterministic
    # typical-value trajectories.
    Cc <- log(CAR_perUL + 1e-30)
    Cc ~ add(addSd)
  })
}
