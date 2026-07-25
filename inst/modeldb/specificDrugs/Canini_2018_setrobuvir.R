Canini_2018_setrobuvir <- function() {
  description <- "Integrated PK / viral-kinetic (VK) model for setrobuvir, a non-nucleoside HCV NS5B polymerase inhibitor, in adults chronically infected with HCV genotype 1a or 1b (Canini 2018). PK is a two-compartment model with first-order absorption and an absorption lag time, parameterised as rate constants (ka, ke, k12, k21) with central volume Vc (Table 3, source rates in /h). PD is a sigmoid Emax inhibition of virion production by central-compartment concentration, with genotype-specific EC50 and Hill coefficient. VK is the standard Neumann-style HCV model reduced to two ODE states (productively infected cells and free virus) under the short-treatment assumption that uninfected target cells remain at their pretreatment steady state; the reduction is dI/dt = d*c*V - d*I and dV/dt = (1 - e(t))*I - c*V (equivalent to the paper's target-fixed form with production rate p normalised to 1, so I represents c*V0 at baseline and its value is a virion-production-rate surrogate rather than an absolute infected-cell count). Genotype (HCV_GT1B binary) switches EC50, Hill, and viral clearance rate c between the GT1a and GT1b typical values and their independent IIVs. Data source: 77 subjects across 4 studies (three healthy-volunteer PK studies A/B/C and one Phase 1 HCV-infected patient PK+VK study D; treated patients received 200, 400, or 800 mg BID for 3 days). Suitable for simulating BID or ascending-dose setrobuvir regimens and the resulting biphasic HCV RNA decline over 3 to 14 days by genotype, and for reproducing the paper's Figure 3 14-day projections. NOT suitable beyond the on-treatment window studied (no drug-resistant-variant emergence submodel, no long-term target-cell repopulation)."

  reference <- paste(
    "Canini L, Lemenuel-Diot A, Brennan BJ, Smith PF, Perelson AS.",
    "(2018). A pharmacokinetic/viral kinetic model to evaluate treatment",
    "of chronic HCV infection with a non-nucleoside polymerase inhibitor.",
    "Antiviral Therapy 23(4):353-361.",
    "doi:10.3851/IMP3216.",
    sep = " "
  )
  vignette <- "Canini_2018_setrobuvir"

  # The PK rate constants (ka, ke, k12, k21) are reported in Table 3 in /h;
  # the VK rate constants (c_1a, c_1b, d) are reported in Table 3 in /day.
  # This packaged model normalises everything to /h inside ini() (VK rates
  # divided by 24 before log-transform), so time in the event table is in
  # hours - consistent with the BID dosing interval (12 h) and the
  # absorption lag time (1.68 h).
  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    HCV_GT1B = list(
      description        = "HCV genotype-1 subtype indicator. 1 = patient infected with HCV genotype 1B; 0 = patient infected with HCV genotype 1A (the source-paper reference subtype).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HCV GT1A; 12 of 27 treated patients in Canini 2018 study D, 44 percent).",
      notes              = "Time-fixed per subject (HCV subtype is determined at infection and does not change over the 3-day treatment window). Switches EC50, Hill coefficient h, and viral clearance rate c between the two independently-estimated typical values and their independent IIVs. Encoding inside model(): ec50_ind = exp(lec50_1a + etalec50_1a) * (1 - HCV_GT1B) + exp(lec50_1b + etalec50_1b) * HCV_GT1B, and analogous switches for hill_ind and c_ind. Registered canonical (see inst/references/covariate-columns.md).",
      source_name        = "Genotype 1a / 1b column in Table 2; the indicator encodes GT1B = 1 to match the naming of HCV_GT1B in the covariate register."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 77L,
    n_studies      = 4L,
    n_patients_vk  = 27L,
    pk_hv_a        = 12L,
    pk_hv_b        = 23L,
    pk_hv_c        = 15L,
    dosing_regimen = "PK studies A/B/C (healthy volunteers): single dose 800 mg or 2000 mg (A); 400/600/800 mg once daily for 10 days after a single D0 dose (B); 200 mg BID for 6 days (C, danoprevir interaction study). VK study D (HCV-infected patients): 200 mg BID, 400 mg BID, or 800 mg BID for 3 days (11, 8, 8 patients per arm; 8 additional patients received placebo and were excluded from the analysis).",
    disease_state  = "Chronic HCV genotype-1 infection (44 percent GT1A, 56 percent GT1B in the 27 treated patients from study D); healthy volunteers for PK studies A, B, and C.",
    dose_range     = "200 mg to 2000 mg single or repeated oral setrobuvir. Multiple-dose regimens all BID or QD as above.",
    regions        = "Canini 2018 does not report study region; study D references Antiviral Therapy dose-ranging trial (Roche).",
    notes          = "Population sourced from Canini 2018 Methods, Table 1 (PK study designs), and Table 2 (baseline characteristics of the 27 treated patients: 5 GT1A + 6 GT1B in arm 1 200 mg; 3 GT1A + 5 GT1B in arm 2 400 mg; 4 GT1A + 4 GT1B in arm 3 800 mg). Median initial viral load ranged from 5.71 to 6.89 log10 IU/mL across dose arms (Table 2). Baseline viral load was not significantly associated with dose group (P = 0.35, Spearman). Baseline demographics beyond genotype and viral load are not tabulated in the trimmed manuscript (see Additional file 1 supplement, not on disk)."
  )

  ini({
    # =========================================================================
    # Setrobuvir PK - two-compartment with first-order absorption and lag time.
    # Rates and volumes from Table 3 (Estimates row). IIVs from Table 3 IIV
    # row, converted from reported CV percent to omega squared via
    # omega^2 = log(1 + CV^2). The residual PK error model in Table 3 footnote
    # is combined additive + proportional on linear-space Cc (mg/L).
    # =========================================================================

    ltlag <- log(1.68);   label("Absorption lag time (h)")                                   # Canini 2018 Table 3: Tlag = 1.68 h (RSE 8%)
    lka   <- log(0.26);   label("First-order absorption rate constant (1/h)")                # Canini 2018 Table 3: ka = 0.26 /h (RSE 19%)
    lvc   <- log(5.11);   label("Central volume of distribution Vc (L)")                     # Canini 2018 Table 3: Vc = 5.11 L (RSE 17%; source uses Vd, single-compartment-volume notation shared with the paper's ke*Vc parameterisation)
    lkel  <- log(0.063);  label("First-order elimination rate constant ke (1/h)")            # Canini 2018 Table 3: ke = 0.063 /h (RSE 17%); half-life = log(2)/ke = 11.0 h
    lk12  <- log(0.12);   label("Central-to-peripheral distribution rate constant k12 (1/h)")# Canini 2018 Table 3: k12 = 0.12 /h (RSE 28%)
    lk21  <- log(0.07);   label("Peripheral-to-central distribution rate constant k21 (1/h)")# Canini 2018 Table 3: k21 = 0.07 /h (RSE 11%)

    # ---- PK IIVs (CV percent -> omega squared) ----
    etaltlag ~ log(1 + 0.52^2)  # Canini 2018 Table 3 IIV: Tlag = 52% -> omega^2 = log(1 + 0.52^2) = 0.2361
    etalka   ~ log(1 + 0.42^2)  # Canini 2018 Table 3 IIV: ka   = 42% -> omega^2 = log(1 + 0.42^2) = 0.1622
    etalvc   ~ log(1 + 0.21^2)  # Canini 2018 Table 3 IIV: Vc   = 21% -> omega^2 = log(1 + 0.21^2) = 0.0431
    etalkel  ~ log(1 + 0.21^2)  # Canini 2018 Table 3 IIV: ke   = 21% -> omega^2 = log(1 + 0.21^2) = 0.0431
    etalk12  ~ log(1 + 0.60^2)  # Canini 2018 Table 3 IIV: k12  = 60% -> omega^2 = log(1 + 0.60^2) = 0.3075
    etalk21  ~ log(1 + 0.42^2)  # Canini 2018 Table 3 IIV: k21  = 42% -> omega^2 = log(1 + 0.42^2) = 0.1622

    # =========================================================================
    # Setrobuvir PD - sigmoid Emax inhibition of virion production. Emax is
    # essentially at the upper bound (0.997, RSE ~0%, no IIV per paper text
    # "IIV on all parameters except Emax was found to best describe the
    # data"). EC50 and Hill coefficient are estimated separately per HCV
    # genotype (1a vs 1b). EC50 in Canini 2018 Table 3 is reported in ng/mL;
    # the packaged model expresses EC50 in mg/L (= ug/mL) to match
    # units$concentration and the Cc unit (mg/L) derived from dose(mg) / Vc(L).
    # 17.8 ng/mL = 0.0178 mg/L; 8.45 ng/mL = 0.00845 mg/L.
    # =========================================================================

    lemax     <- log(0.997)
    label("Maximum drug effectiveness Emax (unitless fraction, <= 1)")                        # Canini 2018 Table 3: Emax = 0.997 (RSE ~0%); no IIV per paper text

    lec50_1a  <- log(17.8e-3)
    label("EC50 for HCV genotype 1a (mg/L; source 17.8 ng/mL)")                              # Canini 2018 Table 3: EC50 GT1A = 17.8 ng/mL (RSE 19%)
    lec50_1b  <- log(8.45e-3)
    label("EC50 for HCV genotype 1b (mg/L; source 8.45 ng/mL)")                              # Canini 2018 Table 3: EC50 GT1B = 8.45 ng/mL (RSE 16%)

    lhill_1a  <- log(3.46)
    label("Hill coefficient for HCV genotype 1a (unitless)")                                 # Canini 2018 Table 3: h GT1A = 3.46 (RSE 12%)
    lhill_1b  <- log(15.6)
    label("Hill coefficient for HCV genotype 1b (unitless)")                                 # Canini 2018 Table 3: h GT1B = 15.6 (RSE 19%)

    etalec50_1a ~ log(1 + 0.45^2)  # Canini 2018 Table 3 IIV: EC50 GT1A = 45% -> omega^2 = log(1 + 0.45^2) = 0.1841
    etalec50_1b ~ log(1 + 0.10^2)  # Canini 2018 Table 3 IIV: EC50 GT1B = 10% -> omega^2 = log(1 + 0.10^2) = 0.00995
    etalhill_1a ~ log(1 + 0.09^2)  # Canini 2018 Table 3 IIV: h    GT1A =  9% -> omega^2 = log(1 + 0.09^2) = 0.00805
    etalhill_1b ~ log(1 + 0.14^2)  # Canini 2018 Table 3 IIV: h    GT1B = 14% -> omega^2 = log(1 + 0.14^2) = 0.01942

    # =========================================================================
    # HCV viral-kinetic parameters. The paper's Table 3 reports c in /day
    # (GT1A: 3.86, GT1B: 6.98) and d in /day (0.0933). The packaged model
    # converts to /h (divide by 24) so PK and VK are on a single time axis.
    # Table 3 does NOT tabulate an IIV for c_1a, c_1b, or d despite the
    # paper text stating IIV was on "all parameters except Emax"; the
    # Results text specifically says d "showed a larger IIV (167%)" which
    # by strict table-position reading corresponds to V0's IIV, and the
    # Discussion cites 163% for d. The packaged model follows the table:
    # V0 (via lrbase) carries the 167% IIV, and c and d are typical-value
    # only. The 163-167% text mention for d is documented in the vignette
    # Errata section. See parameter-names.md policy on unreported IIV: never
    # invent variances.
    # =========================================================================

    lc_1a <- log(3.86 / 24);  label("Viral clearance rate c for HCV genotype 1a (1/h; source 3.86 /day)")   # Canini 2018 Table 3: c GT1A = 3.86 /d (RSE 12%); no IIV tabulated
    lc_1b <- log(6.98 / 24);  label("Viral clearance rate c for HCV genotype 1b (1/h; source 6.98 /day)")   # Canini 2018 Table 3: c GT1B = 6.98 /d (RSE  6%); no IIV tabulated
    ld    <- log(0.0933 / 24); label("Loss rate of productively infected cells d (1/h; source 0.0933 /day)")# Canini 2018 Table 3: d = 0.0933 /d (RSE 52%); no IIV tabulated - text says 163-167% (see Errata)

    # =========================================================================
    # Baseline viral load V0. Paper reports V0 = 6.35 log10 IU/mL. The
    # packaged model carries V0 as the natural-scale IU/mL value (virus
    # compartment initial condition) via lrbase = log(10^6.35) = 14.62; the
    # exponential IIV model then acts multiplicatively on the linear-scale
    # V0, so an IIV of CV=167% translates to individual V0 spanning
    # ~[0.1x, 10x] of the typical value on the linear scale (equivalently
    # ~[-1, +1] log10 IU/mL around 6.35), which matches Table 2's observed
    # range 4.28-7.34 log10 IU/mL.
    # =========================================================================

    lrbase   <- log(10^6.35)
    label("Baseline viral load V0 (log-natural of IU/mL; source Table 3 V0 = 6.35 log10 IU/mL)")     # Canini 2018 Table 3: V0 = 6.35 log10 IU/mL (RSE 2%)
    etalrbase ~ log(1 + 1.67^2)  # Canini 2018 Table 3 IIV: V0 = 167% -> omega^2 = log(1 + 1.67^2) = 1.335

    # =========================================================================
    # Residual error. Table 3 footnote:
    #   PK  - additive a1 = 5.38 ng/mL (RSE 5%) + proportional b = 0.15 (RSE 4%)
    #   VK  - additive a2 = 0.20 log10 IU/mL (RSE 7%) on log10(viral load)
    # =========================================================================

    propSd       <- 0.15;      label("PK proportional residual SD (fraction)")           # Canini 2018 Table 3 footnote: b = 0.15
    addSd        <- 5.38e-3;   label("PK additive residual SD (mg/L; source 5.38 ng/mL)")# Canini 2018 Table 3 footnote: a1 = 5.38 ng/mL = 0.00538 mg/L
    addSd_Vlog10 <- 0.20;      label("VK additive residual SD on log10 viral load (log10 IU/mL)")  # Canini 2018 Table 3 footnote: a2 = 0.20 log10 IU/mL
  })

  model({
    # =========================================================================
    # 1. Individual PK parameters (log-normal exponential-IIV back-transform).
    # =========================================================================
    tlag <- exp(ltlag + etaltlag)
    ka   <- exp(lka   + etalka)
    vc   <- exp(lvc   + etalvc)
    kel  <- exp(lkel  + etalkel)
    k12  <- exp(lk12  + etalk12)
    k21  <- exp(lk21  + etalk21)

    # =========================================================================
    # 2. Individual PD parameters. Emax has no IIV. EC50 and Hill switch
    # between the two genotype-specific typical values and IIVs via
    # HCV_GT1B. Viral clearance c likewise switches (no IIV). c_ind and
    # d_ind are already in /h from ini().
    # =========================================================================
    emax    <- exp(lemax)

    # Genotype-specific PD parameters. Split into per-genotype intermediates
    # so rxode2's mu-reference detector sees each `<theta> + <eta>` on its
    # own line (avoiding the THETA+THETA false-positive when two such
    # sums appear as separate summands on the same RHS - see the
    # rxode2-mu-ref-false-positive-parens-workaround memory).
    ec50_1a_i <- exp(lec50_1a + etalec50_1a)
    ec50_1b_i <- exp(lec50_1b + etalec50_1b)
    ec50_i    <- ec50_1a_i * (1 - HCV_GT1B) + ec50_1b_i * HCV_GT1B

    hill_1a_i <- exp(lhill_1a + etalhill_1a)
    hill_1b_i <- exp(lhill_1b + etalhill_1b)
    hill_i    <- hill_1a_i * (1 - HCV_GT1B) + hill_1b_i * HCV_GT1B

    c_1a_typ  <- exp(lc_1a)
    c_1b_typ  <- exp(lc_1b)
    c_ind     <- c_1a_typ * (1 - HCV_GT1B) + c_1b_typ * HCV_GT1B
    d_ind     <- exp(ld)

    V0_i      <- exp(lrbase + etalrbase)

    # =========================================================================
    # 3. PK ODEs (2-cpt with first-order absorption and lag time). Micro-
    # constants are the primary parameterisation (matches Canini 2018 Table 3
    # which reports ke, k12, k21 rather than CL and Q).
    # =========================================================================
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central +
                          k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # PK observable in mg/L (matches units$concentration).
    Cc <- central / vc

    # =========================================================================
    # 4. Drug effectiveness e(t) - sigmoid Emax on viral production
    # (Canini 2018 Equation 3). Cc is in mg/L; ec50_i is in mg/L (both
    # converted from the paper's ng/mL). A tiny epsilon on Cc keeps the
    # Cc^hill^ term well-defined at t = 0 when Cc = 0; e_t evaluates to 0
    # in that limit because 0^hill / (ec50^hill + 0^hill) = 0.
    # =========================================================================
    Cc_pos <- Cc * (Cc > 0)
    num_e  <- Cc_pos^hill_i
    den_e  <- num_e + ec50_i^hill_i
    e_t    <- emax * num_e / den_e

    # =========================================================================
    # 5. Viral-kinetic ODEs. Standard HCV target-cell model reduced to two
    # states under the paper's short-treatment assumption T(t) = T_0 (target
    # cells constant). With p (virion production per infected cell)
    # normalised to 1, the reduced system is
    #
    #   dI/dt = d * c * V - d * I
    #   dV/dt = (1 - e(t)) * I - c * V
    #
    # which reproduces the pre-treatment steady state I_0 = c * V_0 and
    # V_0 as given (see Canini 2018 Equations 1-2 and the T_0 = c*d/(b*p)
    # constraint immediately following). c and d have units of /h. V is in
    # IU/mL and I is in IU/(mL) as a virion-production-rate surrogate.
    # =========================================================================
    d/dt(infected) <- d_ind * c_ind * virus - d_ind * infected
    d/dt(virus)    <- (1 - e_t) * infected - c_ind * virus

    infected(0) <- c_ind * V0_i
    virus(0)    <- V0_i

    # =========================================================================
    # 6. Observations and residual error. Cc uses combined add+prop on the
    # linear scale (mg/L). Viral-load observation is on log10(V) with an
    # additive SD (Canini 2018 fit the VK model to log10 viral load
    # directly, per Data fitting and statistical methods). A 1e-30 floor
    # inside log10() prevents -Inf if virus underflows to zero after full
    # effectiveness.
    # =========================================================================
    Vlog10 <- log10(virus + 1e-30)

    Cc     ~ add(addSd) + prop(propSd)
    Vlog10 ~ add(addSd_Vlog10)
  })
}
