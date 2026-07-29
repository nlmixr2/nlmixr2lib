Rymut_2023_anti_tryptase <- function() {
  description <- "Mechanistic population PK/PD model for the anti-tryptase IgG4 monoclonal antibody MTPS9579A in healthy adults and adults with moderate-to-severe asthma (Rymut 2023). Two-compartment serum disposition with first-order SC absorption and allometric weight scaling on linear CL and central volume; quasi-equilibrium (QE) TMDD describes saturable binding of MTPS9579A to total monomeric serum tryptase; a mechanistic airway interstitial-fluid (ISF) compartment receives free mAb and mAb-monomer complex from the systemic circulation through lymph flow with vascular reflection coefficients; in the ISF, tryptase is secreted as the active tetramer (target_isf), spontaneously dissociates into inactive monomers (monomer_isf), and is rapidly disrupted by bound mAb (kbreak); free MTPS9579A binds tetramer and monomer with the same KD. Estimated systemic TMDD parameters come from a NONMEM 7.4.3 SAEM fit to 106 healthy Phase 1 subjects (Table 1); fixed mechanistic ISF parameters come from in vitro / physiological literature and from a healthy-subject visual fit of the upper-airway biodistribution coefficient at 3% (Methods, Table S2, Figure S2)."
  reference <- "Rymut SM, Henderson LM, Poon V, Staton TL, Cai F, Sukumaran S, Rhee H, Owen R, Ramanujan S, Yoshida K. A mechanistic PK/PD model to enable dose selection of the potent anti-tryptase antibody (MTPS9579A) in patients with moderate-to-severe asthma. Clin Transl Sci. 2023;16(4):694-703. doi:10.1111/cts.13483"
  vignette <- "Rymut_2023_anti_tryptase"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL (Cc, serum MTPS9579A); nM (TotalSerumTryptase, serum total tryptase); unitless (ActiveAirwayTryptase, fraction of baseline)"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject in source analysis. Allometric power covariate on linear CL (estimated exponent 0.820) and central volume V2 (estimated exponent 0.808), normalised to a 70 kg adult (Rymut 2023 Table 1; reference weight set in the NONMEM control stream Text S1 WTTYP = 70).",
      source_name        = "BWT"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 106L,
    n_studies        = 1L,
    study            = "Phase I single + multiple ascending dose study of MTPS9579A in healthy adults (Methods; underlying Phase 1 study referenced as ref 6 in the paper)",
    age_range        = "Healthy adults (specific age range not tabulated in Rymut 2023 main text or supplement)",
    weight_range     = "Single-ascending-dose (SAD) cohort 49.9-113.3 kg; multiple-ascending-dose (MAD) cohort 47.1-97.7 kg (Rymut 2023 Table 1 footnote a)",
    weight_median    = NA_character_,
    sex_female_pct   = NA_real_,
    race_ethnicity   = NA_character_,
    disease_state    = "Healthy subjects (model-building cohort, n = 106). Prospective simulations were performed for adults with moderate-to-severe asthma using a separate observational cohort (n = 15, 100% White, 26.7% female, mean age 42.7 years, baseline serum total tryptase 10 ng/mL vs 7 ng/mL healthy, baseline nasal active tryptase 4 ng/mL vs 0.4 ng/mL healthy; Table S1).",
    dose_range       = "Single ascending dose 30-300 mg SC or 300-3600 mg IV; multiple ascending dose 150-750 mg SC or 1350-3600 mg IV once every 4 weeks (Methods)",
    regions          = "Phase I study; geography not reported in the paper",
    n_observations   = "835 serum PK + 937 serum total tryptase observations (Methods)",
    notes            = "Population TMDD PK/PD model fit by SAEM in NONMEM 7.4.3 (Text S1; OFV = 7593.39; stochastic + reduced stochastic portions completed). The asthma observational cohort (n = 15) supplied the elevated baseline-tryptase constants used in prospective simulations but did not contribute parameter estimates (Methods)."
  )

  ini({
    # === Systemic PK structural parameters (Rymut 2023 Table 1) ==============
    lka     <- log(0.239);  label("First-order SC absorption rate (ka, 1/day)")                                  # Table 1 row "First order absorption rate (ka)" 0.239 1/day
    lcl     <- log(0.128);  label("Linear clearance (CL, L/day) at the 70 kg reference weight")                  # Table 1 row "Clearance (CL)" 0.128 L/day
    lvc     <- log(3.33);   label("Central volume of distribution (V2, L) at the 70 kg reference weight")        # Table 1 row "Central volume of distribution (V2)" 3.33 L
    lq      <- log(0.408);  label("Inter-compartmental clearance (Q, L/day)")                                    # Table 1 row "Intercompartmental clearance (Q)" 0.408 L/day
    lvp     <- log(2.28);   label("Peripheral volume of distribution (V3, L)")                                   # Table 1 row "Peripheral volume of distribution (V3)" 2.28 L
    lfdepot <- log(0.661);  label("SC bioavailability (Fsc, fraction)")                                          # Table 1 row "Subcutaneous bioavailability (Fsc)" 0.661

    # === Serum QE-TMDD parameters (Rymut 2023 Table 1) =======================
    lkss   <- log(0.0448);  label("QE dissociation constant Kss (nM) for MTPS9579A binding to total monomeric serum tryptase (= KD)")  # Table 1 row "Equilibrium dissociation constant (KD)" 0.0448 nM
    lbase  <- log(0.223);   label("Baseline total monomeric serum tryptase concentration (Base, nM)")                                   # Table 1 row "Baseline total tryptase (Base)" 0.223 nM
    lkdeg  <- log(20.7);    label("First-order tryptase degradation rate constant (kdeg, 1/day)")                                       # Table 1 row "Total tryptase degradation rate constant (kdeg)" 20.7 1/day
    lclint <- log(0.398);   label("Internalisation clearance of the MTPS9579A-total-tryptase complex (CLint, L/day)")                   # Table 1 row "Clearance of MTPS9579A-total tryptase complex (CLint)" 0.398 L/day

    # === Allometric covariate effects on linear CL and central volume (Rymut 2023 Table 1) ===
    e_wt_cl <- 0.820;  label("Allometric exponent of WT/70 on linear CL (unitless)")           # Table 1 row "Effect of weight on CL (exponential model)" 0.82 (RSE 23 percent)
    e_wt_vc <- 0.808;  label("Allometric exponent of WT/70 on central volume V2 (unitless)")    # Table 1 row "Effect of weight on V2 (exponential model)" 0.808 (RSE 18 percent)

    # === Fixed tissue-distribution / mechanism parameters (Rymut 2023 Table S2 + Methods) =====
    # The mechanistic airway PK/PD layer is parameterised from in vitro / physiological
    # literature and from healthy-subject visual fits; only the systemic TMDD block above
    # is re-estimated. Every value below is a structural assumption fixed by the authors,
    # hence wrapped in fixed(). For prospective simulation of adults with moderate-to-severe
    # asthma, the operator may override bl_tryptase_isf (40 ng/mL healthy -> 140 ng/mL asthma
    # per Figure 5b) and the systemic baseline tryptase via lbase (0.223 nM healthy -> ~0.335 nM
    # asthma, i.e. 1.5x healthy) before simulation. See the vignette for worked examples.
    kp_free          <- fixed(0.03);                       label("Biodistribution coefficient of free MTPS9579A from serum to airway ISF (fraction). Visual best fit across SAD + MAD cohorts; supported by the observed nasal-lining-fluid:serum concentration ratio (Methods, Figure 4, Figure S2)")  # Methods "biodistribution coefficient of 3 percent" + Figure S2
    kp_bound         <- fixed(0.03);                       label("Biodistribution coefficient of the MTPS9579A-tryptase complex from serum to airway ISF (fraction); assumed equal to kp_free per Methods")  # Methods "biodistribution coefficient was assumed to be the same for MTPS9579A and the MTPS9579A-tryptase complex"
    refl_lymph       <- fixed(0.2);                        label("Lymphatic reflection coefficient (unitless)")  # Table S2 (Shah and Betts 2012, ref 17)
    klf_tissue       <- fixed(1.2 * 24);                   label("Lymph turnover rate of the airway ISF (1/day); from 1.2 1/h * 24")  # Table S2 "Lymph turnover rate 1.2 1/h"; Methods: 0.2 percent of 182 L/h lung plasma flow / 0.3 L ISF volume
    bl_tryptase_isf  <- fixed(40);                         label("Baseline total tryptase in airway ISF (ng/mL); 40 healthy; 140 in moderate-to-severe asthma (Figure 5b)")  # Table S2 "Baseline total tryptase 40 ng/mL"; Methods + Figure 5b for asthma
    kel_tryp_isf     <- fixed(log(2) / (2 / 24));          label("Tryptase elimination rate constant in ISF (kel, 1/day); from elimination half-life of 2 h")                # Methods "elimination rate constant (kel) of 8.31 day-1 (half-life of 2 h)" + Table S2
    kdiss_tet_isf    <- fixed(log(2) / (0.5 / 24));        label("Spontaneous tetramer dissociation rate constant in ISF (kdiss, 1/day); from dissociation half-life of 30 min")  # Methods "physiologic irreversible dissociation rate constant (kdiss) of 33 day-1 (terminal half-life of 30 min)"
    # Derived rate-balance constants. Computed in ini() (not model()) because the rxode2
    # parser flags <ini_param> + <ini_param> as a between-subject-variability block
    # declaration even on the RHS of an arithmetic assignment, so the sum kel + kdiss
    # cannot be assembled inside model() from the two pieces above. Storing the
    # rate-balance numerics here keeps the per-line provenance comment in one place.
    sum_kel_isf      <- fixed(log(2) / (2 / 24) + log(2) / (0.5 / 24));                 label("Tetramer-removal aggregate rate constant kel + kdiss in ISF (1/day)")  # Derived from kel_tryp_isf + kdiss_tet_isf (Methods)
    f_tet_isf        <- fixed((log(2) / (2 / 24)) / (log(2) / (2 / 24) + log(2) / (0.5 / 24)));  label("Baseline fraction of ISF total tryptase mass present as active tetramer (= kel/(kel+kdiss))")  # Methods rate-balance partition (Text S2)
    kbreak_tet       <- fixed(1000);                       label("MTPS9579A-induced rapid tetramer disruption rate constant (kbreak, 1/day); from disruption half-life of 1 min")  # Methods "kdiss of 1000 day-1 (half-life of 1 min)" + Table S2
    kel_mono_ab      <- fixed(log(2) / (100 / 24));        label("Elimination rate constant of the airway monomer-MTPS9579A complex (1/day); from assumed complex half-life of 100 h (negligible)")  # Table S2 "Degradation half-life of mAb-monomeric tryptase complex 100 h"
    # mAb-tryptase association rate: 7.62e5 1/M/s -> 1/nM/day. Convert: 7.62e5 [1/M/s] * (1 mol/L / 1e9 nmol/L) * (86400 s/day) = 65.83 [1/nM/day].
    kon_isf          <- fixed(7.62e5 / 1e9 * 86400);       label("MTPS9579A-tryptase association rate constant in ISF (1/nM/day); converted from 7.62e5 1/M/s")  # Table S2 "Binding constant (Kon) 7.62e5 1/Ms" + in vitro binding (ref 5)
    kd_isf           <- fixed(0.0448);                     label("Equilibrium dissociation constant in ISF (KD, nM); equal to the serum Kss")  # Table S2 "KD 4.88e-11 M ~= 0.0488 nM" (Table 1 reports 0.0448 nM; same value used in Text S2)
    # MTPS9579A molecular weight: needed to convert mg dose -> nmol (state) and to back-convert nM -> ug/mL for the observed Cc output.
    mw_ab            <- fixed(155);                        label("MTPS9579A molecular weight (kDa = ug/nmol); IgG4 monoclonal antibody (Text S2)")           # Text S2 MWab = 155.0 ug/nmol
    mw_mono          <- fixed(32);                         label("Tryptase monomer molecular weight (kDa) used in ISF total-tryptase mass balance")          # Text S2 MWmono = 32.0 ug/nmol
    mw_tet           <- fixed(128);                        label("Tryptase tetramer molecular weight (kDa = 4 * MWmono) used in ISF total-tryptase mass balance")  # Text S2 MWtet = 128.0 ug/nmol

    # === IIV (Rymut 2023 Table 1; exponential model on log-parameters) =====
    # NONMEM .lst Final Parameter Estimates (Text S1, BLOCK structure preserved):
    #   ETA(1) on KA   variance 0.214
    #   ETA(2)+ETA(3) on CL,V2 block: var_CL = 0.0809, cov(CL,V2) = 0.0342, var_V2 = 0.0387 (corr ~= 0.612)
    #   ETA(4) on Kss  variance 0   FIXED -> not encoded
    #   ETA(5) on Base variance 0.197
    #   ETA(6) on kdeg variance 0.103
    #   ETA(7) on CLint variance 0.135
    etalka                ~ 0.214                                           # Table 1 row ka omega-squared 0.214
    etalcl + etalvc       ~ c(0.0809,
                              0.0342, 0.0387)                              # Table 1 rows CL omega-squared 0.0809 and V2 omega-squared 0.0387; off-diagonal from BLOCK(2) Text S1
    etalbase              ~ 0.197                                           # Table 1 row Base omega-squared 0.197
    etalkdeg              ~ 0.103                                           # Table 1 row kdeg omega-squared 0.103
    etalclint             ~ 0.135                                           # Table 1 row CLint omega-squared 0.135

    # === Residual error (Rymut 2023 Table 1) ================================
    # NONMEM thetarized parameterisation with SIGMA fixed at 1:
    #   SD = sqrt(prop_term^2 + add_term^2), prop_term = IPRED * theta_prop, add_term = theta_add.
    # Internal NONMEM units are nM throughout (the QE math, Kss, Base are all in nM, and the
    # AMT column was pre-converted from mg to nmol). For Cc reported in ug/mL, propSd transfers
    # unchanged (unit-free fraction) and addSd transforms as 2.72 nM * MW_ab(kDa) / 1000 = 0.4216 ug/mL.
    propSd                    <- 0.103;            label("Proportional residual error on serum MTPS9579A (fraction)")  # Table 1 row "Thetarized proportional error for MTPS9579A" 0.103
    addSd                     <- 2.72 * 155 / 1000; label("Additive residual error on serum MTPS9579A (ug/mL); converted from 2.72 nM via mAb MW 155 kDa")  # Table 1 row "Thetarized additive error for MTPS9579A" 2.72 (nM, converted)
    propSd_TotalSerumTryptase <- 0.245;            label("Proportional residual error on serum total tryptase (fraction)")  # Table 1 row "Thetarized proportional error for total tryptase" 0.245
  })

  model({
    # ----- Individual systemic PK / TMDD parameters (allometric on linear CL and V2) -----
    wt_cl <- (WT / 70)^e_wt_cl
    wt_vc <- (WT / 70)^e_wt_vc
    ka    <- exp(lka    + etalka)
    cl    <- exp(lcl    + etalcl)    * wt_cl
    vc    <- exp(lvc    + etalvc)    * wt_vc
    q     <- exp(lq)
    vp    <- exp(lvp)
    fsc   <- exp(lfdepot)
    kss   <- exp(lkss)
    base  <- exp(lbase  + etalbase)
    kdeg  <- exp(lkdeg  + etalkdeg)
    clint <- exp(lclint + etalclint)
    ksyn  <- base * kdeg  # zero-order target synthesis (steady state with kdeg * base)

    # ----- Reflection coefficients back-calculated from Kp (Methods) -----
    # Solving Kp = (1 - refl_tissue) / (1 - refl_lymph) for refl_tissue, given equal influx
    # and efflux lymph flow rates, gives refl_tissue = 1 - Kp * (1 - refl_lymph).
    refl_tissue_free  <- 1 - kp_free  * (1 - refl_lymph)
    refl_tissue_bound <- 1 - kp_bound * (1 - refl_lymph)

    # ----- Initial ISF tryptase concentrations (rate-balance partition; Text S2) -----
    # Total tryptase in ISF (ng/mL) is distributed between active tetramer and inactive monomer
    # by rate balance: f_tet = kel / (kel + kdiss), f_mono = kdiss / (kel + kdiss). The assay
    # reads tetramer in monomer-equivalent mass; division by the species MW converts mass to molar.
    # sum_kel_isf (= kel + kdiss) and f_tet_isf (= kel / sum_kel) are precomputed in ini() to
    # avoid the rxode2 parser's <ini_param> + <ini_param> BSV-block ambiguity.
    tet_0    <- bl_tryptase_isf * f_tet_isf       / mw_tet   # nM
    mono_0   <- bl_tryptase_isf * (1 - f_tet_isf) / mw_mono  # nM
    # Steady-state zero-order tetramer secretion rate; d/dt(target_isf) = 0 at baseline with no
    # drug aboard.
    vin_tryp <- sum_kel_isf * tet_0  # nM/day

    # ----- ISF binding rate constants -----
    koff_isf <- kon_isf * kd_isf  # 1/day

    # ----- QE-TMDD free / bound mAb concentrations in serum (Text S1 / Gibiansky 2008) -----
    # central tracks TOTAL mAb amount (nmol); total_target tracks TOTAL serum tryptase
    # concentration (nM). The QE quadratic for free mAb concentration is the standard
    # numerically-stable root.
    ctot      <- central / vc                                           # total mAb conc (nM)
    qss_disc  <- (ctot - total_target - kss)^2 + 4 * kss * ctot
    cfree     <- 0.5 * (ctot - total_target - kss) + 0.5 * sqrt(qss_disc)  # free mAb (nM)
    cbound    <- ctot - cfree                                           # mAb-monomer complex (nM)

    # ----- Lymph fluxes to / from the airway ISF (Text S2) -----
    # All fluxes carry units of nM/day; ISF concentrations are tracked directly because the
    # lymph turnover term klf_tissue absorbs the ISF volume into the rate constant.
    inf_free_isf  <- klf_tissue * (1 - refl_tissue_free)  * cfree              # nM/day
    inf_bound_isf <- klf_tissue * (1 - refl_tissue_bound) * cbound             # nM/day
    out_free_isf  <- klf_tissue * (1 - refl_lymph) * mab_isf                   # nM/day
    out_bound_isf <- klf_tissue * (1 - refl_lymph) * complex_monomer_isf       # nM/day

    # ----- ODE system -----
    # Systemic compartments (amounts in nmol; total_target tracks concentration in nM).
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot + q * (peripheral1 / vp) - q * cfree - cl * cfree - clint * total_target * cfree / (kss + cfree)
    d/dt(peripheral1)  <-  q * cfree - q * (peripheral1 / vp)
    d/dt(total_target) <-  ksyn - kdeg * total_target - (clint / vc - kdeg) * cfree * total_target / (kss + cfree)
    total_target(0)    <- base

    # ISF compartments (all states in nM):
    #   mab_isf              -- free MTPS9579A in airway ISF
    #   target_isf           -- free active tetrameric tryptase in airway ISF
    #   complex_isf          -- MTPS9579A-tetramer complex in airway ISF (no first-order
    #                           elim; converted to monomer-mAb complex via kbreak)
    #   monomer_isf          -- free inactive monomeric tryptase in airway ISF
    #   complex_monomer_isf  -- MTPS9579A-monomer complex in airway ISF (also receives the
    #                           systemic mAb-monomer complex via lymph influx)
    d/dt(mab_isf) <- inf_free_isf - out_free_isf - kon_isf * mab_isf * target_isf + koff_isf * complex_isf - kon_isf * mab_isf * monomer_isf + koff_isf * complex_monomer_isf

    d/dt(target_isf) <- vin_tryp - kel_tryp_isf * target_isf - kdiss_tet_isf * target_isf - kon_isf * mab_isf * target_isf + koff_isf * complex_isf
    target_isf(0) <- tet_0

    d/dt(complex_isf) <- kon_isf * mab_isf * target_isf - koff_isf * complex_isf - kbreak_tet * complex_isf

    d/dt(monomer_isf) <- -kel_tryp_isf * monomer_isf - kon_isf * mab_isf * monomer_isf + koff_isf * complex_monomer_isf + 3 * kbreak_tet * complex_isf + 4 * kdiss_tet_isf * target_isf
    monomer_isf(0) <- mono_0

    d/dt(complex_monomer_isf) <- kon_isf * mab_isf * monomer_isf - koff_isf * complex_monomer_isf + kbreak_tet * complex_isf - kel_mono_ab * complex_monomer_isf + inf_bound_isf - out_bound_isf

    # ----- Bioavailability + mg -> nmol unit conversion -----
    # Dose AMT is in mg; the state `central` is in nmol so that ctot = central / vc is in nM
    # (matching the units of Kss and Base). 1 mg of MTPS9579A = 1000 / mw_ab nmol.
    f(depot)   <- fsc * 1000 / mw_ab
    f(central) <- 1000 / mw_ab

    # ----- Observation outputs -----
    # Cc: serum MTPS9579A in ug/mL (= mg/L). 1 nM = 1 nmol/L * MW(g/mol) * 1e-9 g/nmol
    #     = MW(kDa) * 1e-3 mg/L = MW(kDa) / 1000 ug/mL.
    Cc <- ctot * mw_ab / 1000

    # Serum total tryptase in nM (model native; multiply by mw_mono to obtain ng/mL).
    TotalSerumTryptase <- total_target

    # Active airway tryptase relative to baseline (Rymut 2023 Figure 3 / Table 2):
    # ratio of tetramer to its baseline value (1 = baseline, 0 = full inhibition).
    ActiveAirwayTryptase <- target_isf / tet_0

    Cc                 ~ add(addSd) + prop(propSd)
    TotalSerumTryptase ~ prop(propSd_TotalSerumTryptase)
  })
}
