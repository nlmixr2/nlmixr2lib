HillMcManus_2017_febuxostat_lesinurad <- function() {
  description <- "Semi-mechanistic dual-drug PKPD model for the impact of non-adherence to febuxostat (xanthine oxidase inhibitor) plus lesinurad (URAT1 uricosuric) urate-lowering therapy in gout; combines published 2-compartment first-order absorption PK for each drug with a 4-compartment xanthine / uric-acid PD system (Hill-McManus 2017)"
  reference <- "Hill-McManus D, Soto E, Marshall S, Lane S, Hughes D. Impact of non-adherence on the safety and efficacy of uric acid-lowering therapies in the treatment of gout. Br J Clin Pharmacol. 2018;84(1):142-152. doi:10.1111/bcp.13427"
  vignette <- "HillMcManus_2017_febuxostat_lesinurad"
  units <- list(
    time = "hour",
    dosing = "mg",
    concentration = "mg/dL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed. Additive modifier of febuxostat CL/F (b_WT = 0.155 (dL/h) per kg, no reference subtraction) and power modifier of lesinurad Vc/F (reference 70 kg).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Creatinine clearance (Cockcroft-Gault, actual body weight)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per simulation. Additive modifier of febuxostat CL/F (b_CrCl = 0.142 (dL/h) per mL/min), power modifier of lesinurad CL/F (reference 87 mL/min), and power modifier of lesinurad direct-Emax E_Dmax (reference 87 mL/min). Hill-McManus 2017 applied a 15 mL/min downward shift to Cockcroft-Gault estimates and excluded values < 30 mL/min when matching the CRYSTAL cohort (Results p. 147).",
      source_name        = "CrCl"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 154,
    n_studies      = NA_integer_,
    age_range      = "healthy adults (TMX-99-001, ref 24); gout-patient simulation parameterised from CRYSTAL trial (study 304, ref 25)",
    weight_range   = "simulation cohort modelled from CRYSTAL (log-normal); reference 70 kg (lesinurad Vc) and 99 kg (representative-male 15-dL daily urine assumption, p. 147)",
    sex_female_pct = 5,
    disease_state  = "Gout patients with hyperuricaemia (CRYSTAL trial-representative simulation); two phenotypes considered -- renal under-excreters and overproducers of uric acid",
    dose_range     = "Febuxostat 80 mg PO QD; lesinurad 200 mg or 400 mg PO QD; combinations as in CRYSTAL",
    regions        = NA_character_,
    notes          = "PK parameter values were reproduced by Hill-McManus 2017 Table 2 from upstream publications (febuxostat: Khosravan 2006, ref 22; lesinurad: FDA review, ref 23). The 154 healthy volunteers (TMX-99-001 dose escalation; ref 24) supplied weighted-average CLX and CLUA estimates used as literature-fixed PD system parameters. Gout-patient simulation (Methods 'Gout patient simulation model') uses CRYSTAL-trial baselines (pretreatment sUA log-normal mean 8.83 mg/dL, SD 1.53; 95 percent male) with phenotype-specific scaling of system parameters (Table 3): under-excreter scales CLUA by theta4 / BUA, overproducer scales BX by BUA / theta4. These phenotype scalings are applied at simulation time and are demonstrated in the validation vignette."
  )

  ini({
    # ============================================================
    # Febuxostat 2-compartment PK with first-order absorption + lag
    # Source: Hill-McManus 2017 Table 2 (reproduced from Khosravan
    # 2006 ref 22). CL/F is ADDITIVE in CrCl and WT (Table 2 footnote a):
    #   CL/F = CL/F_0 + b_CrCl * CrCl + b_WT * WT
    # ============================================================
    lka_febx       <- log(13.7);  label("Febuxostat absorption rate constant (Ka, 1/h)")                           # Table 2: Ka 13.7 h^-1, BSV CV% 176
    lcl_febx       <- log(49.3);  label("Febuxostat baseline apparent clearance CL/F_0 (dL/h)")                   # Table 2: CL/F_0 49.3 dL/h, BSV CV% 18.3
    lvc_febx       <- log(322);   label("Febuxostat apparent central volume Vc/F_0 (dL)")                          # Table 2: Vc/F_0 322 dL, BSV NE
    lvp_febx       <- log(222);   label("Febuxostat apparent peripheral volume Vp/F (dL)")                         # Table 2: Vp/F 222 dL, BSV NE
    lq_febx        <- log(55.7);  label("Febuxostat apparent intercompartmental clearance Q/F (dL/h)")             # Table 2: Q/F 55.7 dL/h, BSV NE
    ltlag_febx     <- log(0.23);  label("Febuxostat absorption lag time Tlag (h)")                                 # Table 2: Tlag 0.23 h, BSV NE
    e_crcl_cl_febx <- fixed(0.142); label("CrCl additive effect on febuxostat CL/F ((dL/h) per mL/min)")          # Table 2: b_CrCl 0.142, NA
    e_wt_cl_febx   <- fixed(0.155); label("WT additive effect on febuxostat CL/F ((dL/h) per kg)")                # Table 2: b_WT 0.155, NA

    # ============================================================
    # Lesinurad 2-compartment PK with first-order absorption + lag
    # Source: Hill-McManus 2017 Table 2 (reproduced from FDA review
    # ref 23). CL/F is power on CrCl (reference 87 mL/min); Vc/F is
    # power on WT (reference 70 kg) (Table 2 footnote a and b):
    #   CL/F = CL/F_0 * (CrCl / 87)^b_CrCl
    #   Vc/F = Vc/F_0 * (WT  / 70)^b_WT
    # ============================================================
    lka_lesn       <- log(0.69);  label("Lesinurad absorption rate constant (Ka, 1/h)")                            # Table 2: Ka 0.69 h^-1, BSV CV% 121.7
    lcl_lesn       <- log(69.9);  label("Lesinurad apparent clearance CL/F_0 at CrCl 87 mL/min (dL/h)")           # Table 2: CL/F_0 69.9 dL/h, BSV CV% 63.4
    lvc_lesn       <- log(241);   label("Lesinurad apparent central volume Vc/F_0 at WT 70 kg (dL)")               # Table 2: Vc/F_0 241 dL, BSV CV% 12.2
    lvp_lesn       <- log(83);    label("Lesinurad apparent peripheral volume Vp/F (dL)")                          # Table 2: Vp/F 83 dL, BSV CV% 20.5
    lq_lesn        <- log(4.48);  label("Lesinurad apparent intercompartmental clearance Q/F (dL/h)")              # Table 2: Q/F 4.48 dL/h, BSV NE
    ltlag_lesn     <- log(0.233); label("Lesinurad absorption lag time Tlag (h)")                                  # Table 2: Tlag 0.233 h, BSV CV% 38.9
    e_crcl_cl_lesn <- fixed(0.322); label("Lesinurad CrCl power exponent on CL/F (unitless, reference 87 mL/min)") # Table 2: b_CrCl 0.322, NA
    e_wt_vc_lesn   <- fixed(0.511); label("Lesinurad WT power exponent on Vc/F (unitless, reference 70 kg)")       # Table 2: b_WT 0.511, NA

    # ============================================================
    # System PD parameters (Figure 1 + Table 1, healthy-subject typical values)
    # Amounts in mg, volumes in dL, clearances in dL/h.
    # ============================================================
    bx     <- 8.94;          label("Baseline amount of xanthine BX in serum (mg)")                                 # Table 1 theta_1, Estimated
    vx     <- 333;           label("Volume of xanthine distribution VX (dL)")                                       # Table 1 theta_2, Estimated
    clx    <- fixed(10.57);  label("Renal clearance of xanthine CLX (dL/h)")                                        # Table 1 theta_3, Literature
    bua    <- 703;           label("Baseline amount of uric acid BUA in serum (mg)")                                # Table 1 theta_4, Estimated
    vua    <- 154;           label("Volume of uric acid distribution VUA (dL)")                                     # Table 1 theta_5, Estimated
    clua   <- fixed(4.11);   label("Renal clearance of uric acid CLUA (dL/h)")                                      # Table 1 theta_6, Literature

    # ============================================================
    # Febuxostat PD: indirect inhibition INH1 (on k0, xanthine production)
    # and INH2 (on k1, xanthine -> UA formation), plus stimulation STIM1
    # (on k2, xanthine renal clearance step-function up at very low doses).
    # ============================================================
    emax_1 <- fixed(3);      label("Febuxostat STIM1 Emax on xanthine renal clearance (unitless)")                 # Table 1 theta_7, Assumed (10 mg dose -> 3-5x baseline)
    ec50_1 <- fixed(0.001);  label("Febuxostat STIM1 EC50 on xanthine renal clearance (mg/dL)")                    # Table 1 theta_8, Assumed (low-dose surrogate)
    imax_1 <- fixed(1);      label("Febuxostat INH1 Imax on xanthine production (unitless)")                       # Table 1 theta_9, Assumed
    lic501 <- log(0.1320);   label("Febuxostat INH1 IC50_1 on xanthine production (mg/dL)")                        # Table 1 theta_10, Estimated; eta3 SD^2 = 0.2
    imax_2 <- fixed(1);      label("Febuxostat INH2 Imax on xanthine -> UA formation (unitless)")                  # Table 1 theta_11, Assumed
    lic502 <- log(0.00113);  label("Febuxostat INH2 IC50_2 on UA formation (mg/dL)")                               # Table 1 theta_12, Estimated; eta3 SD^2 = 0.2

    # ============================================================
    # Lesinurad PD: direct-effect Emax parameters (literature, ref 23).
    # The indirect-response STIM2 parameters acting on k3 (UA renal
    # clearance) are derived from these in model() using the steady-
    # state equations on p. 146:
    #   Emax_2 = E0 / (E0 - E_Dmax * (CrCl/87)^b_CrCl) - 1
    #   EC50_2 = (Emax_2 * EC_D50) / (E0 / (E0 - E_Dmax / 2) - 1) - EC_D50
    # ============================================================
    e0                 <- fixed(6.77);   label("Lesinurad reference baseline sUA E0 used in direct-Emax derivation (mg/dL)")  # Table 1 theta_13, Literature
    le_dmax            <- log(2.55);     label("Lesinurad direct-Emax maximum reduction in sUA E_Dmax (mg/dL)")               # Table 1 theta_14, Literature; eta4 SD^2 = 0.346
    e_crcl_e_dmax_lesn <- fixed(0.564);  label("Lesinurad CrCl power exponent on E_Dmax (unitless, reference 87 mL/min)")     # Table 1 theta_15, Literature
    ec_d50             <- fixed(0.0974); label("Lesinurad direct-Emax EC_D50 (mg/dL)")                                         # Table 1 theta_16, Literature

    # ============================================================
    # IIV
    #
    # PK between-subject variability (Table 2 CV%): omega^2 = log(1 + CV^2)
    #   etalka_febx   CV 176%   -> log(1 + 1.76^2)  = 1.4117
    #   etalcl_febx   CV 18.3%  -> log(1 + 0.183^2) = 0.03315
    #   etalka_lesn   CV 121.7% -> log(1 + 1.217^2) = 0.91858
    #   etalcl_lesn   CV 63.4%  -> log(1 + 0.634^2) = 0.33129
    #   etalvc_lesn   CV 12.2%  -> log(1 + 0.122^2) = 0.01481
    #   etalvp_lesn   CV 20.5%  -> log(1 + 0.205^2) = 0.04113
    #   etaltlag_lesn CV 38.9%  -> log(1 + 0.389^2) = 0.14133
    #
    # PD random effects (Table 1 reports SD^2 directly; "Error model used: theta_i = theta_u exp(eta_i)"):
    #   etalic501_lic502  shared eta_3 on IC50_1 and IC50_2 (Methods p. 147: "IC50 parameters were
    #                      assumed to vary according to eta_3 with a coefficient of variation of 20%").
    #                      The Table 1 SD^2 = 0.2 entries appear in BOTH theta_10 and theta_12 rows
    #                      because the same eta_3 acts on both. CV 20% -> log(1 + 0.20^2) = 0.0392.
    #   etale_dmax        eta_4 on E_Dmax (Table 1 theta_14 SD^2 = 0.346).
    # ============================================================
    etalka_febx       ~ 1.4117
    etalcl_febx       ~ 0.03315
    etalka_lesn       ~ 0.91858
    etalcl_lesn       ~ 0.33129
    etalvc_lesn       ~ 0.01481
    etalvp_lesn       ~ 0.04113
    etaltlag_lesn     ~ 0.14133
    etalic501_lic502  ~ 0.0392
    etale_dmax        ~ 0.346
  })

  model({
    # --- Effective PK parameters with covariate effects ---
    # Febuxostat: CL/F additive in CrCl and WT (Table 2 footnote a)
    cl_febx <- exp(lcl_febx + etalcl_febx) +
               e_crcl_cl_febx * CRCL +
               e_wt_cl_febx   * WT
    vc_febx <- exp(lvc_febx)
    vp_febx <- exp(lvp_febx)
    q_febx  <- exp(lq_febx)
    ka_febx <- exp(lka_febx + etalka_febx)
    tlag_febx_h <- exp(ltlag_febx)

    # Lesinurad: CL/F power on CrCl (ref 87), Vc/F power on WT (ref 70)
    cl_lesn <- exp(lcl_lesn + etalcl_lesn) * (CRCL / 87)^e_crcl_cl_lesn
    vc_lesn <- exp(lvc_lesn + etalvc_lesn) * (WT   / 70)^e_wt_vc_lesn
    vp_lesn <- exp(lvp_lesn + etalvp_lesn)
    q_lesn  <- exp(lq_lesn)
    ka_lesn <- exp(lka_lesn + etalka_lesn)
    tlag_lesn_h <- exp(ltlag_lesn + etaltlag_lesn)

    # Micro-constants
    kel_febx <- cl_febx / vc_febx
    k12_febx <- q_febx  / vc_febx
    k21_febx <- q_febx  / vp_febx
    kel_lesn <- cl_lesn / vc_lesn
    k12_lesn <- q_lesn  / vc_lesn
    k21_lesn <- q_lesn  / vp_lesn

    # --- PK ODEs (dose in mg into depot_*, central in mg, conc = central/Vc in mg/dL) ---
    d/dt(depot_febx)       <- -ka_febx * depot_febx
    d/dt(central_febx)     <-  ka_febx * depot_febx - kel_febx * central_febx - k12_febx * central_febx + k21_febx * peripheral1_febx
    d/dt(peripheral1_febx) <-  k12_febx * central_febx - k21_febx * peripheral1_febx
    lag(depot_febx) <- tlag_febx_h

    d/dt(depot_lesn)       <- -ka_lesn * depot_lesn
    d/dt(central_lesn)     <-  ka_lesn * depot_lesn - kel_lesn * central_lesn - k12_lesn * central_lesn + k21_lesn * peripheral1_lesn
    d/dt(peripheral1_lesn) <-  k12_lesn * central_lesn - k21_lesn * peripheral1_lesn
    lag(depot_lesn) <- tlag_lesn_h

    # Drug plasma concentrations (mg/dL) feeding PD drug functions
    cf_t <- central_febx / vc_febx
    cl_t <- central_lesn / vc_lesn

    # --- Individual PD parameters ---
    ic50_1 <- exp(lic501 + etalic501_lic502)
    ic50_2 <- exp(lic502 + etalic501_lic502)

    # Lesinurad direct-Emax parameters with CrCl covariate (Table 1 footnote)
    e_dmax_i <- exp(le_dmax + etale_dmax) * (CRCL / 87)^e_crcl_e_dmax_lesn

    # Indirect-response derivation from direct-Emax (Hill-McManus 2017 p. 146)
    emax_2 <- e0 / (e0 - e_dmax_i)        - 1
    ec50_2 <- (emax_2 * ec_d50) /
              (e0 / (e0 - e_dmax_i / 2) - 1) - ec_d50

    # --- PD system rate parameters (Figure 1 Eqs 5-8)
    # k2 = CLX/VX; k3 = CLUA/VUA;
    # k1 = k3 * BUA / (BX * (M_UA / M_X)) (xanthine -> UA, mass-balance stoichiometry)
    # k0 = (k1 + k2) * BX (no-drug steady state of d(A_X)/dt = 0)
    m_x  <- 152.11   # xanthine MW (g/mol)
    m_ua <- 168.11   # uric acid MW (g/mol)
    k2 <- clx  / vx
    k3 <- clua / vua
    k1 <- k3 * bua / (bx * (m_ua / m_x))
    k0 <- (k1 + k2) * bx

    # --- PD drug functions (Figure 1 Eqs 9-12)
    # Figure 1 Eqs 9-10 are written as INH = IC50 / (IC50 + C), which is the
    # special case of the general inhibitory Emax form 1 - Imax * C / (IC50 + C)
    # when Imax = 1 (both Imax_1 and Imax_2 are Assumed = 1 in Table 1).
    # The general form is used here so the Imax_1 and Imax_2 parameter values
    # remain auditable and overridable.
    inh1  <- 1 - imax_1 * cf_t / (ic50_1 + cf_t)
    inh2  <- 1 - imax_2 * cf_t / (ic50_2 + cf_t)
    stim1 <- 1 + emax_1 * cf_t / (ec50_1 + cf_t)
    stim2 <- 1 + emax_2 * cl_t / (ec50_2 + cl_t)

    # --- PD ODEs (Figure 1 Eqs 1-4)
    # The UA production term `k1 * inh2 * xanthine` is scaled by the
    # stoichiometric mass-conversion factor (m_ua / m_x) to keep the
    # no-drug steady-state of the urate compartment consistent with
    # Figure 1 Eq 7 (which defines k1 = k3 * BUA / (BX * (m_ua / m_x))).
    # Figure 1 Eq 2 as printed omits this factor; including it is required
    # for mass balance because the xanthine state is in mg xanthine and the
    # urate state is in mg uric acid (1 mol xanthine -> 1 mol UA).
    d/dt(xanthine)       <-  k0 * inh1 - k1 * inh2 * xanthine - k2 * stim1 * xanthine
    d/dt(urate)          <-  k1 * inh2 * xanthine * (m_ua / m_x) - k3 * stim2 * urate
    d/dt(xanthine_urine) <-  k2 * stim1 * xanthine
    d/dt(urate_urine)    <-  k3 * stim2 * urate

    # Drug-free steady-state initial conditions
    xanthine(0)       <- bx
    urate(0)          <- bua
    xanthine_urine(0) <- 0
    urate_urine(0)    <- 0

    # Derived clinical outputs
    sUA <- urate    / vua        # serum uric acid concentration (mg/dL)
    sX  <- xanthine / vx         # serum xanthine concentration  (mg/dL)
    uUA <- urate_urine           # cumulative urinary urate    amount (mg)
    uX  <- xanthine_urine        # cumulative urinary xanthine amount (mg)
  })
}
