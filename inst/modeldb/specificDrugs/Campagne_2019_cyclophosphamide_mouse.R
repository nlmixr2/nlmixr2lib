Campagne_2019_cyclophosphamide_mouse <- function() {
  description <- "Preclinical (mouse). Plasma and brain/tumor extracellular-fluid (ECF) population PK model for cyclophosphamide (CTX) and its sequential metabolites 4-hydroxy-cyclophosphamide (4OH-CTX) and carboxyethylphosphoramide mustard (CEPM) in female CD-1 nude mice (non-tumor-bearing and orthotopic Group 3 medulloblastoma G3MB), following a single 130 mg/kg intraperitoneal dose of cyclophosphamide (Campagne 2019). Three sequential two-compartment plasma sub-models are linked by full (Fm = 1) conversion CTX -> 4OH-CTX -> CEPM (so reported CL and V for the two metabolites are apparent CL/F and V/F); each compound additionally has a one-compartment brain/tumor ECF sub-model linked to its plasma central via influx (CLin) and efflux (CLef) clearances driven by the unbound plasma concentration FU x Cp. ECF volume fixed at 0.001 L/kg (Stewart 2010, ref 26 of source). No covariate effects retained; pooled fit across non-tumor-bearing and G3MB mice."
  reference <- "Campagne O, Davis A, Zhong B, Nair S, Haberman V, Patel YT, Janke L, Roussel MF, Stewart CF. CNS Penetration of Cyclophosphamide and Metabolites in Mice Bearing Group 3 Medulloblastoma and Non-Tumor Bearing Mice. J Pharm Pharm Sci. 2019;22(1):553-568. doi:10.18433/jpps30608"
  vignette <- "Campagne_2019_cyclophosphamide_mouse"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")
  # The model is body-weight-normalised: structural clearances/volumes
  # carry per-kg units (L/h/kg, L/kg), and dose `amt` is interpreted as
  # umol per kg body weight. Stating dosing = "umol" (rather than
  # "umol/kg") keeps the lint dose-concentration dimension check happy
  # while the dose-per-body-weight convention is documented in
  # population$notes and in the vignette Assumptions section.

  covariateData <- list(
    # Campagne 2019 Methods (Pharmacokinetic modeling) tested
    # cohort / tumor-status differences post hoc on individual
    # Kp,uu exposures (two-tailed Welch t-test), not as a model
    # covariate effect. Table 2 reports a single pooled parameter
    # set across NTB and G3MB mice; no covariates are carried
    # in the structural model.
  )

  population <- list(
    species        = "mouse (CD-1 nude, female)",
    n_subjects     = 41L,
    n_studies      = 5L,
    sex_female_pct = 100,
    disease_state  = "Non-tumor-bearing (NTB) and orthotopic Group 3 medulloblastoma (G3MB, 1e5 luciferase-transduced cells stereotactically implanted) mouse model of pediatric medulloblastoma.",
    dose_range     = "130 mg/kg cyclophosphamide IP single dose (498 umol/kg using MW = 261.09 g/mol).",
    regions        = "USA (St. Jude Children's Research Hospital).",
    n_observations = "Plasma: 30 samples in the plasma-only study + 14 + 23 + 27 + 21 = 115 samples across the four microdialysis studies. ECF: 25 + 38 + 45 + 35 = 143 dialysate intervals across the four microdialysis studies; only the M3 (NTB, n=9) and M4 (G3MB, n=7) ECF data contributed to the final ECF fit (M1 and M2 ECF data were excluded because the derivatizing solution in the perfusate caused brain hemorrhage and probe-site contamination; see Results 'Impact of different methods to collect 4OH-CTX samples in ECF').",
    cohorts        = "Plasma-only PK study: 12 NTB mice with sparse sampling at 5 min, 0.25, 0.5, 1, 1.5, 2, 4, 6 h. Microdialysis study M1: 5 NTB, derivatizing solution in perfusate (ECF excluded). M2: 8 G3MB, derivatizing solution in perfusate (ECF excluded). M3: 9 NTB, derivatizing solution in collection tubes. M4: 7 G3MB, derivatizing solution in collection tubes. Plasma sampling during microdialysis used a limited-sampling design at 0.25, 1, 2 h.",
    notes          = "Plasma model parameters (Campagne 2019 Table 2 upper block) were estimated from pooled plasma data from all five studies; ECF model parameters (Table 2 lower block) from M3 + M4 ECF data fit simultaneously with each mouse's plasma posterior. Modelling in Monolix v2018R1 with SAEM; log-normal IIV, proportional residual error. Plasma protein binding (median FU): CTX 0.26, 4OH-CTX 0.39, CEPM 0.31 (Results 'Plasma protein binding studies'). Fraction-of-metabolite-formed Fm fixed to 1.0 for CTX -> 4OH-CTX and 4OH-CTX -> CEPM for identifiability (Results 'Plasma pharmacokinetics'); so the metabolite CL and V are apparent (CL/Fm and V/Fm). MW values for the dose-unit conversion: cyclophosphamide 261.09 g/mol, 4-hydroxy-cyclophosphamide 277.09 g/mol, carboxyethylphosphoramide mustard 293.08 g/mol -- but molar (1:1) conversion is preserved through the metabolic cascade in the model, so the same umol/kg amount carries through CTX, 4OH-CTX, and CEPM compartments without MW scaling."
  )

  ini({
    # ============================================================
    # Plasma model -- Campagne 2019 Table 2 upper block (population
    # estimates). RSE% in parentheses for context. CL is L/h/kg,
    # V is L/kg.
    # ============================================================

    # Cyclophosphamide (parent)
    lcl    <- log(4.4);     label("Cyclophosphamide central clearance CL (L/h/kg)")           # Table 2: CL = 4.4 (3.7% RSE)
    lvc    <- log(0.65);    label("Cyclophosphamide central volume V (L/kg)")                 # Table 2: V = 0.65 (5.0% RSE)
    lq     <- log(0.18);    label("Cyclophosphamide inter-compartmental clearance Q (L/h/kg)") # Table 2: Q = 0.18 (4.3% RSE)
    lvp    <- log(0.062);   label("Cyclophosphamide peripheral volume Vp (L/kg)")             # Table 2: Vp = 0.062 (2.6% RSE)

    # 4-Hydroxy-cyclophosphamide -- apparent CL/F, V/F with Fm = 1
    lcl_4ohctx <- log(11);    label("Apparent 4OH-CTX central clearance CL/Fm (L/h/kg)")        # Table 2: CL = 11 (4.9% RSE); Fm = 1 fixed
    lvc_4ohctx <- log(2.4);   label("Apparent 4OH-CTX central volume V/Fm (L/kg)")              # Table 2: V = 2.4 (6.2% RSE)
    lq_4ohctx  <- log(0.075); label("Apparent 4OH-CTX inter-compartmental clearance Q/Fm (L/h/kg)") # Table 2: Q = 0.075 (17% RSE)
    lvp_4ohctx <- log(0.22);  label("Apparent 4OH-CTX peripheral volume Vp/Fm (L/kg)")          # Table 2: Vp = 0.22 (20% RSE)

    # Carboxyethylphosphoramide mustard -- apparent CL/F, V/F with Fm = 1
    lcl_cepm <- log(6.4);   label("Apparent CEPM central clearance CL/Fm (L/h/kg)")            # Table 2: CL = 6.4 (3.0% RSE); Fm = 1 fixed
    lvc_cepm <- log(1.1);   label("Apparent CEPM central volume V/Fm (L/kg)")                  # Table 2: V = 1.1 (5.3% RSE)
    lq_cepm  <- log(1.6);   label("Apparent CEPM inter-compartmental clearance Q/Fm (L/h/kg)") # Table 2: Q = 1.6 (9.0% RSE)
    lvp_cepm <- log(1.7);   label("Apparent CEPM peripheral volume Vp/Fm (L/kg)")              # Table 2: Vp = 1.7 (7.1% RSE)

    # ============================================================
    # ECF (brain / tumor extracellular fluid) model
    # -- Campagne 2019 Table 2 lower block. CLin / CLef in L/h/kg.
    # ============================================================
    lclin    <- log(4.3e-4);   label("Cyclophosphamide ECF influx clearance CLin (L/h/kg)")     # Table 2: CLin = 4.3e-4 (16% RSE)
    lclef    <- log(2.4e-3);   label("Cyclophosphamide ECF efflux clearance CLef (L/h/kg)")     # Table 2: CLef = 2.4e-3 (7.0% RSE)
    lclin_4ohctx <- log(2.3e-4); label("4OH-CTX ECF influx clearance CLin (L/h/kg)")               # Table 2: CLin = 2.3e-4 (25% RSE)
    lclef_4ohctx <- log(3.3e-3); label("4OH-CTX ECF efflux clearance CLef (L/h/kg)")               # Table 2: CLef = 3.3e-3 (4.9% RSE)
    lclin_cepm   <- log(8.8e-5); label("CEPM ECF influx clearance CLin (L/h/kg)")                  # Table 2: CLin = 8.8e-5 (17% RSE)
    lclef_cepm   <- log(1.5e-3); label("CEPM ECF efflux clearance CLef (L/h/kg)")                  # Table 2: CLef = 1.5e-3 (7.2% RSE)

    # ECF volume fixed -- shared across the three compounds; Stewart
    # 2010 reference 26 of the source paper.
    lvecf <- fixed(log(0.001)); label("ECF central volume Vecf (L/kg), fixed per Campagne 2019 ref 26 (Stewart 2010)") # Methods 'Tumor and brain ECF pharmacokinetic modeling'

    # ============================================================
    # Plasma protein binding (fraction unbound) -- fixed to median
    # values measured by rapid equilibrium dialysis in CD-1 mouse
    # plasma (Campagne 2019 Results 'Plasma protein binding
    # studies'). Only unbound drug crosses the BBB, so CLin acts on
    # FU x Cp.
    # ============================================================
    fu        <- fixed(0.26); label("Cyclophosphamide plasma fraction unbound FU (unitless), fixed") # Results: median FU = 0.26 (range 0.24-0.28)
    fu_4ohctx <- fixed(0.39); label("4OH-CTX plasma fraction unbound FU (unitless), fixed")          # Results: median FU = 0.39 (range 0.28-0.48)
    fu_cepm   <- fixed(0.31); label("CEPM plasma fraction unbound FU (unitless), fixed")             # Results: median FU = 0.31 (range 0.29-0.34)

    # ============================================================
    # Inter-individual variability -- Campagne 2019 Table 2.
    # IIV reported as standard deviations of eta on the log scale
    # (Table 2 footnote: "Inter-individual variabilities are
    # reported in terms of standard deviations"; log-normal
    # parameterization Pi = thetaP * exp(eta_P)). Variance stored
    # here is the squared SD. Only the parameters listed in
    # Table 2 carry IIV; the rest are typical-value only.
    # ============================================================
    etalcl         ~ 0.11^2     # Table 2: eta-CL_CTX SD = 0.11 (17% RSE)
    etalcl_4ohctx  ~ 0.13^2     # Table 2: eta-CL_4OH-CTX SD = 0.13 (19% RSE)
    etalcl_cepm    ~ 0.094^2    # Table 2: eta-CL_CEPM SD = 0.094 (25% RSE)
    etalvc_cepm    ~ 0.22^2     # Table 2: eta-V_CEPM SD = 0.22 (21% RSE); CTX and 4OH-CTX V have no reported IIV
    etalclin       ~ 0.31^2     # Table 2: eta-CLin_CTX SD = 0.31 (51% RSE)
    etalclef       ~ 0.23^2     # Table 2: eta-CLef_CTX SD = 0.23 (22% RSE)
    etalclin_4ohctx ~ 0.86^2    # Table 2: eta-CLin_4OH-CTX SD = 0.86 (19% RSE); CLef_4OH-CTX has no reported IIV
    etalclin_cepm  ~ 0.62^2     # Table 2: eta-CLin_CEPM SD = 0.62 (19% RSE)
    etalclef_cepm  ~ 0.20^2     # Table 2: eta-CLef_CEPM SD = 0.20 (38% RSE)

    # ============================================================
    # Residual error -- proportional only (Campagne 2019 Table 2).
    # Values are SD on the linear concentration scale (Monolix
    # default proportional parameterization for log-normal data).
    # Parent observations use the bare propSd / propSd_Cecf form;
    # metabolite outputs add a registered metabolite suffix.
    # ============================================================
    propSd              <- 0.24; label("Cyclophosphamide plasma proportional residual SD (fraction)") # Table 2: eps_prop = 0.24 (11% RSE)
    propSd_4ohctx       <- 0.29; label("4OH-CTX plasma proportional residual SD (fraction)")          # Table 2: eps_prop = 0.29 (11% RSE)
    propSd_cepm         <- 0.18; label("CEPM plasma proportional residual SD (fraction)")             # Table 2: eps_prop = 0.18 (12% RSE)
    propSd_Cecf         <- 0.45; label("Cyclophosphamide ECF proportional residual SD (fraction)")    # Table 2: eps_prop_ECF = 0.45 (12% RSE)
    propSd_Cecf_4ohctx  <- 0.33; label("4OH-CTX ECF proportional residual SD (fraction)")             # Table 2: eps_prop_ECF = 0.33 (16% RSE)
    propSd_Cecf_cepm    <- 0.22; label("CEPM ECF proportional residual SD (fraction)")                # Table 2: eps_prop_ECF = 0.22 (15% RSE)
  })

  model({
    # ---- Individual PK parameters --------------------------------
    cl    <- exp(lcl + etalcl)
    vc    <- exp(lvc)
    q     <- exp(lq)
    vp    <- exp(lvp)

    cl_4ohctx <- exp(lcl_4ohctx + etalcl_4ohctx)
    vc_4ohctx <- exp(lvc_4ohctx)
    q_4ohctx  <- exp(lq_4ohctx)
    vp_4ohctx <- exp(lvp_4ohctx)

    cl_cepm <- exp(lcl_cepm + etalcl_cepm)
    vc_cepm <- exp(lvc_cepm + etalvc_cepm)
    q_cepm  <- exp(lq_cepm)
    vp_cepm <- exp(lvp_cepm)

    # ECF influx / efflux clearances
    clin         <- exp(lclin + etalclin)
    clef         <- exp(lclef + etalclef)
    clin_4ohctx  <- exp(lclin_4ohctx + etalclin_4ohctx)
    clef_4ohctx  <- exp(lclef_4ohctx)
    clin_cepm    <- exp(lclin_cepm + etalclin_cepm)
    clef_cepm    <- exp(lclef_cepm + etalclef_cepm)

    vecf <- exp(lvecf)

    # ---- Observed concentrations ---------------------------------
    Cc                <- central            / vc
    Cc_4ohctx         <- central_4ohctx     / vc_4ohctx
    Cc_cepm           <- central_cepm       / vc_cepm
    Cecf              <- ecf                / vecf
    Cecf_4ohctx       <- ecf_4ohctx         / vecf
    Cecf_cepm         <- ecf_cepm           / vecf

    # ---- ODE system ----------------------------------------------
    # Compartments hold molar amounts per kg body weight (umol/kg);
    # volumes are L/kg. Conc Cc = amount/Vc has units umol/L = uM,
    # matching the units the source paper reports its concentration
    # data in. The molar 1:1 stoichiometry of CTX -> 4OH-CTX (loss
    # of 2H, gain of 1 OH; molecular formula change) and
    # 4OH-CTX -> CEPM (further oxidation; net +O, -2H) preserves
    # mole count through the cascade.
    #
    # Plasma central terms (left to right): elimination to next
    # metabolite, inter-compartmental flow to peripheral, BBB
    # influx loss to ECF, BBB efflux return from ECF. ECF terms:
    # BBB influx from plasma (driven by unbound plasma conc FU x Cc),
    # BBB efflux to plasma. Per Campagne 2019 Methods "only the
    # transfer of each compound between plasma and ECF compartments
    # were modeled" -- no brain metabolism, no ECF -> next-metabolite
    # conversion.

    d/dt(central) <-
      -cl * Cc -
      q * (Cc - peripheral1 / vp) -
      clin * fu * Cc + clef * Cecf
    d/dt(peripheral1) <-
      q * (Cc - peripheral1 / vp)
    d/dt(ecf) <-
      clin * fu * Cc - clef * Cecf

    d/dt(central_4ohctx) <-
      cl * Cc -
      cl_4ohctx * Cc_4ohctx -
      q_4ohctx * (Cc_4ohctx - peripheral1_4ohctx / vp_4ohctx) -
      clin_4ohctx * fu_4ohctx * Cc_4ohctx + clef_4ohctx * Cecf_4ohctx
    d/dt(peripheral1_4ohctx) <-
      q_4ohctx * (Cc_4ohctx - peripheral1_4ohctx / vp_4ohctx)
    d/dt(ecf_4ohctx) <-
      clin_4ohctx * fu_4ohctx * Cc_4ohctx - clef_4ohctx * Cecf_4ohctx

    d/dt(central_cepm) <-
      cl_4ohctx * Cc_4ohctx -
      cl_cepm * Cc_cepm -
      q_cepm * (Cc_cepm - peripheral1_cepm / vp_cepm) -
      clin_cepm * fu_cepm * Cc_cepm + clef_cepm * Cecf_cepm
    d/dt(peripheral1_cepm) <-
      q_cepm * (Cc_cepm - peripheral1_cepm / vp_cepm)
    d/dt(ecf_cepm) <-
      clin_cepm * fu_cepm * Cc_cepm - clef_cepm * Cecf_cepm

    # ---- Residual error ------------------------------------------
    Cc            ~ prop(propSd)
    Cc_4ohctx     ~ prop(propSd_4ohctx)
    Cc_cepm       ~ prop(propSd_cepm)
    Cecf          ~ prop(propSd_Cecf)
    Cecf_4ohctx   ~ prop(propSd_Cecf_4ohctx)
    Cecf_cepm     ~ prop(propSd_Cecf_cepm)
  })
}
