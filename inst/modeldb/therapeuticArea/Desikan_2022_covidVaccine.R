Desikan_2022_covidVaccine <- function() {
  description <- paste(
    "QSP (immunodynamics, deterministic).",
    "Two-antigen, three-epitope, multi-B-cell / antibody model of the",
    "humoral recall response to a two-strain vaccine sequence (ancestral",
    "Wuhan strain WU and a novel variant OM = SARS-CoV-2 Omicron BA.1;",
    "the framework is generic and is also validated against H5N1",
    "influenza-vaccine data in Figure 5 of Desikan et al. 2022).",
    "Each antigen carries a conserved epitope class (C, shared across",
    "WU and OM) and a variant-unique epitope class (W on WU, O on OM).",
    "Fourteen ordinary differential equations track: two free antigens",
    "(agW, agV), six antibody-bound antigen intermediates (agWw / agWc",
    "/ agWwc / agVv / agVc / agVvc, where the subscripts denote which",
    "epitope class carries a bound antibody), three B-cell populations",
    "producing antibodies specific to each epitope class (Bc, Bw, Bv),",
    "and three corresponding antibody pools (Ac, Aw, Av).",
    "B-cell proliferation follows a saturating dose-response in the",
    "amount of cognate antigen exposed (free-epitope form); previously",
    "primed (memory) B cells can additionally be recruited by",
    "non-homotypic antigen at fraction f_b of the homotypic rate, which",
    "provides the original-antigenic-sin (OAS) mechanism.",
    "Reported observables (used to reproduce all figures) are the",
    "multiplicity-weighted antibody titres against each virus variant:",
    "Ab_WU = m_conserved*Ac + n_unique*Aw and",
    "Ab_OM = m_conserved*Ac + n_unique*Av,",
    "with m:n = 1:5 (Table 1).",
    "All ten model constants are FIXED at Table 1 values; there is no",
    "IIV and no residual error (deterministic simulation is the intended",
    "use). Vaccinations enter as bolus dose events on the appropriate",
    "free-antigen compartment (agW for WU-vaccinations; agV for",
    "OM-vaccinations); at the FIRST WU vaccination the naive WU-unique",
    "B-cell and antibody pools are seeded with n_unique units on Bw and",
    "Aw, and at the FIRST OM vaccination the naive OM-unique pools are",
    "seeded with n_unique units on Bv and Av. Conserved-epitope pools",
    "Bc and Ac start at their naive m_conserved units regardless of",
    "which vaccine is given first. See the validation vignette",
    "(Desikan_2022_covidVaccine) for reference vaccination schedules and",
    "reproductions of Figures 2A, 3, 4A, 4B and 5.",
    sep = " "
  )

  reference <- paste(
    "Desikan R, Linderman SL, Davis C, Zarnitsyna VI, Ahmed H, Antia R",
    "(2022). Vaccine models predict rules for updating vaccines against",
    "evolving pathogens such as SARS-CoV-2 and influenza in the context",
    "of pre-existing immunity.",
    "Frontiers in Immunology 13:985478.",
    "doi:10.3389/fimmu.2022.985478.",
    "Model equations, parameter values, dosing schedule and observation",
    "convention were transcribed from the paper's Materials and Methods",
    "(Table 1 and unnumbered equations for the WU response, with the",
    "OM response derived by symmetry) and from the authors' Zenodo /",
    "Frontiers supplementary MATLAB code (Data Sheet 1 archive:",
    "model.m, WW.m, OO.m, WWO.m, WWWWW.m, WWOOO.m, WWOW.m, WWWO.m and",
    "the MS_plots.m driver used to render Figures 2A, 3, 4, 5 and SI",
    "Figure 1).",
    sep = " "
  )

  vignette <- "Desikan_2022_covidVaccine"

  paper_specific_compartments <- c(
    "agW", "agV",
    "agWw", "agWc", "agWwc",
    "agVv", "agVc", "agVvc",
    "Bc", "Bw", "Bv",
    "Ac", "Aw", "Av"
  )

  units <- list(
    time = "day",
    dosing = paste(
      "Antigen dose amounts are on the model's dimensionless 'scaled",
      "concentration' (s) scale; Table 1 defines Ag_dose = 1e5 s for",
      "vaccinations #1 and #2 and 0.5e5 s for vaccinations #3, #4 and",
      "#5. Vaccinations enter as bolus dose events on the appropriate",
      "free-antigen compartment (agW for WU vaccinations, agV for OM",
      "vaccinations). At the first WU vaccination the vignette adds a",
      "one-time seeding dose of n_unique = 5 to each of Bw and Aw, and",
      "at the first OM vaccination it adds n_unique = 5 to each of Bv",
      "and Av, matching the naive-precursor initial conditions used in",
      "WW.m / OO.m and the state resets in WWO.m / WWOOO.m."
    ),
    concentration = paste(
      "Antibody titres and B-cell counts are reported on arbitrary",
      "units (AU); Table 1 states that the initial concentration of B",
      "cells is rescaled to 1 and 'Antibody = B cell at steady state'.",
      "For the H5N1 comparison in Figure 5 the paper anchors AU to a",
      "5e3 background titre floor."
    )
  )

  covariateData <- list()

  population <- list(
    species = paste(
      "Species-agnostic mechanistic model. Calibrated qualitatively",
      "against three data sources: (i) mice (BALB/c) two-dose WU vs OM",
      "primary vaccination (Ying et al. 2022, Figure 2A);",
      "(ii) non-human primates (rhesus macaques) WU-WU-WU vs WU-WU-OM",
      "three-dose Moderna mRNA-1273 vs mRNA-Omicron sequence (Gagne",
      "et al. 2022, Figure 3); (iii) humans, both healthy WU-vaccinated",
      "recipients plus post-Omicron-infection subjects (Khan et al.",
      "2022, Figure 2C), and H5N1-vaccinated adults (Ellebedy et al.",
      "2014 Figure 5B / 5D, reused in Figure 5)."
    ),
    n_subjects = "Framework model, not fitted to a single cohort; qualitative fits to the mouse (n approx 5-10 per group), primate (n = 5-10) and human (n = 30-40) cohorts cited above.",
    n_studies = 5L,
    age_range = "Mouse, primate and human adult cohorts (species-specific).",
    disease_state = "Healthy vaccine recipients (SARS-CoV-2 or H5N1 influenza) and, in the Khan 2022 dataset, community SARS-CoV-2 Omicron-infected individuals with or without prior WU vaccination.",
    dose_range = "Dose amounts are dimensionless (scaled concentration units s). Reference schedule (Gagne et al. simulation, Figure 3): vaccinations at days 0, 28, 287, 525 and 700 (weeks 0, 4, 41, 75, 100); vaccinations #1 and #2 each deliver 1e5 s, vaccinations #3-#5 each deliver 0.5e5 s. The mouse two-dose schedule (Figure 2A, WW / OO) uses Ag_dose = 1e3 with a 21-day boost interval; the human OM-infection schedule (Figure 2C) uses Ag_dose = 1e5 with a 28+139+2 day sequence.",
    regions = "Not applicable (mechanistic model).",
    notes = "The scale of the initial B-cell / antibody pool encodes epitope multiplicity (m_conserved = 1 for the single conserved epitope class C; n_unique = 5 for the five variant-unique epitopes on each of WU and OM; Table 1 states 'results are qualitatively similar for other values of m and n'). No individual-level covariates are used."
  )

  ini({
    # ================================================================
    # ALL VALUES ARE FIXED CONSTANTS (Desikan 2022 Table 1). No IIV, no
    # residual error. Reproducing the numerical fits in the paper (Figs
    # 2A, 3, 4A/B, 5, SI Fig 1) requires the exact values below and the
    # dosing schedule described in units$dosing.
    # ================================================================

    # ---- Rate constants (day^-1 or s^-1 day^-1) ----
    k_bind   <- fixed(0.0005)   ; label("Rate constant for antibody-antigen binding (s^-1 day^-1)")               # Table 1 row 'k'
    d_ag     <- fixed(1)        ; label("First-order decay rate of free and bound antigen (day^-1)")             # Table 1 row 'd_Ag'
    d_ab     <- fixed(0.1)      ; label("First-order decay rate of antibody (day^-1)")                            # Table 1 row 'd_Ab'
    lambda_b <- fixed(1)        ; label("Maximum proliferation rate of B cells (day^-1)")                         # Table 1 row 'lambda' (shown as l in the table)
    phi      <- fixed(100)      ; label("Antigen amount that half-maximally stimulates B-cell proliferation (s)") # Table 1 row 'phi' (shown as f_dose in the table)
    p_ab     <- fixed(0.1)      ; label("Antibody production rate per unit B cell (day^-1)")                      # Table 1 row 'p'
    d_b      <- fixed(log(2)/47); label("First-order decay rate of B cells (day^-1)")                             # Table 1 row 'd_B' = ln(2)/47 (47-day half-life; paper cites Cromer et al. 2021)

    # ---- Dimensionless mixing / structural fractions ----
    f_b      <- fixed(0.075)    ; label("Fractional B-cell stimulation in secondary responses by non-homotypic antigen (dimensionless)") # Table 1 row 'f'

    # ---- Epitope multiplicities (dimensionless) ----
    # Table 1 fixes m:n = 1:5. These multiplicities enter as (i) the
    # naive precursor levels for the corresponding B-cell and antibody
    # pools (Bc(0) = Ac(0) = m_conserved; Bw(0)/Aw(0) or Bv(0)/Av(0)
    # seeded to n_unique on the first respective vaccination), and (ii)
    # the coefficients that convert per-epitope-class antibody titres
    # to per-antigen titres in the observables (Ab_WU = m*Ac + n*Aw,
    # Ab_OM = m*Ac + n*Av).
    m_conserved <- fixed(1)     ; label("Number of conserved-epitope classes per antigen (m in Table 1)") # Table 1 row 'm:n = 1:5'
    n_unique    <- fixed(5)     ; label("Number of unique-epitope classes per antigen (n in Table 1)")    # Table 1 row 'm:n = 1:5'
  })

  model({
    # ================================================================
    # ODE system for the humoral recall response (Desikan 2022,
    # Materials and Methods, and MATLAB supplement model.m).
    #
    # Compartments (in ODE declaration order; the same order used in
    # the MATLAB state vector x(1..14)):
    #   agW    : free WU antigen                                (x1)
    #   agV    : free OM antigen                                (x2)
    #   agWw   : WU antigen with a W-epitope Ab bound            (x3)
    #   agWc   : WU antigen with a C-epitope Ab bound            (x4)
    #   agWwc  : WU antigen with both epitope classes bound      (x5)
    #   agVv   : OM antigen with an O-epitope Ab bound           (x6)
    #   agVc   : OM antigen with a C-epitope Ab bound            (x7)
    #   agVvc  : OM antigen with both epitope classes bound      (x8)
    #   Bc     : conserved-epitope-specific B cells              (x9)
    #   Bw     : WU-unique-epitope-specific B cells              (x10)
    #   Bv     : OM-unique-epitope-specific B cells              (x11)
    #   Ac     : conserved-epitope-specific antibody             (x12)
    #   Aw     : WU-unique-epitope-specific antibody             (x13)
    #   Av     : OM-unique-epitope-specific antibody             (x14)
    #
    # Vaccination enters as a bolus dose event on agW (WU vaccine)
    # or agV (OM vaccine). Bw / Aw are seeded on the first WU
    # vaccination and Bv / Av on the first OM vaccination via
    # additional dose events in the event table (see vignette). Bc
    # and Ac start at their naive m_conserved level at t = 0 via the
    # initial-condition assignments below.
    # ================================================================

    # -------- Antigen kinetics (WU side) --------
    # Free WU antigen is cleared by first-order decay (rate d_ag) and
    # by consumption when an antibody specific to any epitope class
    # present on WU (W or C) binds. Table 1 row 'd_Ag' and unnumbered
    # equation dW/dt in the Methods.
    d/dt(agW)   <- -k_bind * agW * (Aw + Ac) - d_ag * agW                                # model.m line 33: dx(1)
    d/dt(agWw)  <-  k_bind * agW * Aw    - k_bind * agWw * Ac - d_ag * agWw              # model.m line 35: dx(3)
    d/dt(agWc)  <-  k_bind * agW * Ac    - k_bind * agWc * Aw - d_ag * agWc              # model.m line 36: dx(4)
    d/dt(agWwc) <-  k_bind * agWw * Ac + k_bind * agWc * Aw - d_ag * agWwc               # model.m line 37: dx(5)

    # -------- Antigen kinetics (OM side, by symmetry) --------
    # Same structure as the WU side, with V/O in place of W and the
    # variant-unique antibody Av in place of Aw. Derived by symmetry
    # from the Methods; encoded explicitly in model.m lines 34, 38-40.
    d/dt(agV)   <- -k_bind * agV * (Av + Ac) - d_ag * agV                                # model.m line 34: dx(2)
    d/dt(agVv)  <-  k_bind * agV * Av    - k_bind * agVv * Ac - d_ag * agVv              # model.m line 38: dx(6)
    d/dt(agVc)  <-  k_bind * agV * Ac    - k_bind * agVc * Av - d_ag * agVc              # model.m line 39: dx(7)
    d/dt(agVvc) <-  k_bind * agVv * Ac + k_bind * agVc * Av - d_ag * agVvc               # model.m line 40: dx(8)

    # -------- B-cell dynamics --------
    # Each B-cell population is decayed first-order at d_b (47-day
    # half-life) and stimulated by a saturating dose-response in the
    # total amount of cognate antigen present with the target epitope
    # class free of bound antibody. Non-cognate antigen contributes at
    # a reduced fraction f_b (original antigenic sin: memory cells
    # cross-react with the antigenically altered epitope at low
    # efficiency).
    #
    # - Bc responds to any antigen with a free C epitope (both WU and
    #   OM antigens share the conserved epitope, so both contribute
    #   fully). The "C-epitope-free" pool is (agW + agWw + agV + agVv):
    #   free antigen and antigen with only the unique-epitope class
    #   bound still expose the C epitope.
    # - Bw responds to WU antigen with a free W epitope (agW + agWc,
    #   i.e. free WU antigen and WU antigen with only the C epitope
    #   bound; the W epitope is unique to WU so OM contributes only
    #   at fraction f_b via cross-reactivity of memory cells).
    # - Bv responds symmetrically to OM antigen with a free O epitope
    #   (agV + agVc), with WU antigen contributing at fraction f_b.
    d/dt(Bc) <- -d_b * Bc + lambda_b * Bc *
                (agW + agWw + agV + agVv) /
                (phi + agW + agWw + agV + agVv)                                          # model.m line 41: dx(9)

    d/dt(Bw) <- -d_b * Bw + lambda_b * Bw *
                (agW + agWc + f_b * agV + f_b * agVc) /
                (phi + agW + agWc + f_b * agV + f_b * agVc)                              # model.m line 42: dx(10)

    d/dt(Bv) <- -d_b * Bv + lambda_b * Bv *
                (agV + agVc + f_b * agW + f_b * agWc) /
                (phi + agV + agVc + f_b * agW + f_b * agWc)                              # model.m line 43: dx(11)

    # -------- Antibody dynamics --------
    # Each antibody pool is decayed first-order at d_ab (10-day half-
    # life) and produced at rate p_ab per unit B cell. Table 1 states
    # this is 'rescaled so that Antibody = B cell at steady state';
    # d_ab / p_ab = 0.1 / 0.1 = 1, i.e. Ac_ss = Bc, so the naive Ac(0)
    # matches Bc(0) = m_conserved.
    d/dt(Ac) <- -d_ab * Ac + p_ab * Bc                                                    # model.m line 44: dx(12)
    d/dt(Aw) <- -d_ab * Aw + p_ab * Bw                                                    # model.m line 45: dx(13)
    d/dt(Av) <- -d_ab * Av + p_ab * Bv                                                    # model.m line 46: dx(14)

    # -------- Naive-precursor initial conditions --------
    # Bc / Ac start at their m_conserved naive level regardless of
    # which vaccine is given first; the WU-unique and OM-unique pools
    # (Bw / Aw / Bv / Av) start at zero and are seeded on the first
    # respective vaccination via dose events in the event table (see
    # vignette; MATLAB precedent: WW.m Bw_0 / Aw_0 = 5 with Bv_0 /
    # Av_0 = 0, and the state resets Bv <- 5, Av <- 5 at the first OM
    # exposure in WWO.m / WWOOO.m).
    Bc(0) <- m_conserved
    Ac(0) <- m_conserved

    # -------- Observables (multiplicity-weighted variant titres) --------
    # The paper reports antibody titres against a variant as the
    # multiplicity-weighted sum over its epitope classes. All figures
    # in Desikan 2022 (2A, 3, 4A, 4B, 5, SI Fig 1) and their MATLAB
    # plotting counterparts (MS_plots.m lines 63-67, 191-195, 238-256,
    # 280-284, 358-361, 384-388, 425-431, 494-533) use this form:
    #   Ab_WU = m_conserved * Ac + n_unique * Aw
    #   Ab_OM = m_conserved * Ac + n_unique * Av
    # Separate epitope-class titres (Ab_uniqueOM and Ab_conserved)
    # are reported in Figure 4B and Figure 5.
    Ab_WU        <- m_conserved * Ac + n_unique * Aw
    Ab_OM        <- m_conserved * Ac + n_unique * Av
    Ab_uniqueWU  <- n_unique * Aw
    Ab_uniqueOM  <- n_unique * Av
    Ab_conserved <- m_conserved * Ac
  })
}
