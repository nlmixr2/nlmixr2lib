Scheuher_2023_ADC_mouse_qsp <- function() {
  description <- "QSP. Mouse in vivo platform model for HER2-targeting antibody-drug conjugates (T-DM1 default, N87 tumor xenograft; T-DXd variant supported via parameter overrides). Combines: (i) mouse plasma PK for ADC / naked antibody / free payload in central + peripheral compartments; (ii) mechanistic tumor uptake via a Krogh cylinder + surface exchange model; (iii) intracellular ADC processing (HER2 binding, endocytosis, recycling, degradation, endosomal payload release, cytosol transport, tubulin binding); and (iv) tumor growth inhibition via a Simeoni-style 4-stage transit-chain cascade with Hill-type kill on the proliferating stage. Amounts in nmol; concentrations amount/volume."
  reference   <- "Scheuher B, Ghusinga KR, McGirr K, Nowak M, Panday S, Apgar J, Subramanian K, Betts A. Towards a platform quantitative systems pharmacology (QSP) model for preclinical to clinical translation of antibody drug conjugates (ADCs). J Pharmacokinet Pharmacodyn. 2023;51(1):5-30. doi:10.1007/s10928-023-09884-6. Mouse in vivo model = Tables S1b, S2b-c, S3c-d."
  vignette    <- "Scheuher_2023_ADC_platform_qsp"
  units       <- list(time = "h", dosing = "nmol", concentration = "nM")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    adc_central         = list(analyte = "ADC", units = "nmol", specimen = "plasma", verified = FALSE),
    ab_central          = list(analyte = "naked antibody", units = "nmol", specimen = "plasma", verified = FALSE),
    pl_central          = list(analyte = "free payload", units = "nmol", specimen = "plasma", verified = FALSE),
    adc_peripheral      = list(analyte = "ADC", units = "nmol", specimen = "plasma", verified = FALSE),
    ab_peripheral       = list(analyte = "naked antibody", units = "nmol", specimen = "plasma", verified = FALSE),
    pl_peripheral       = list(analyte = "free payload", units = "nmol", specimen = "plasma", verified = FALSE),
    adc_ext_tumor       = list(analyte = "ADC", units = "nmol", specimen = "tumor", verified = FALSE),
    ab_ext_tumor        = list(analyte = "naked antibody", units = "nmol", specimen = "tumor", verified = FALSE),
    pl_ext_tumor        = list(analyte = "free payload", units = "nmol", specimen = "tumor", verified = FALSE),
    n1                  = list(analyte = "not applicable", units = "nmol", specimen = "not applicable", verified = FALSE),
    n2                  = list(analyte = "not applicable", units = "nmol", specimen = "not applicable", verified = FALSE),
    n3                  = list(analyte = "not applicable", units = "nmol", specimen = "not applicable", verified = FALSE),
    n4                  = list(analyte = "not applicable", units = "nmol", specimen = "not applicable", verified = FALSE),
    her2_surf_tumor     = list(analyte = "HER2", units = "nmol", specimen = "tumor", verified = FALSE),
    her2_adc_surf_tumor = list(analyte = "ADC-HER2 complex", units = "nmol", specimen = "tumor", verified = FALSE),
    her2_ab_surf_tumor  = list(analyte = "naked antibody-HER2 complex", units = "nmol", specimen = "tumor", verified = FALSE),
    her2_endo_tumor     = list(analyte = "endosomal HER2", units = "nmol", specimen = "tumor", verified = FALSE),
    her2_adc_endo_tumor = list(analyte = "ADC-HER2 endosome complex", units = "nmol", specimen = "tumor", verified = FALSE),
    her2_ab_endo_tumor  = list(analyte = "naked antibody-HER2 endosome complex", units = "nmol", specimen = "tumor", verified = FALSE),
    pl_endo_tumor       = list(analyte = "free payload in endosomes", units = "nmol", specimen = "tumor", verified = FALSE),
    pl_cyto_tumor       = list(analyte = "free payload in cytosol", units = "nmol", specimen = "tumor", verified = FALSE),
    t_cyto_tumor        = list(analyte = "tubulin-bound payload", units = "nmol", specimen = "tumor", verified = FALSE),
    tpl_cyto_tumor      = list(analyte = "payload bound to tubulin", units = "nmol", specimen = "tumor", verified = FALSE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list()

  population <- list(
    species        = "mouse (BT-474EEI, BT-474, N87, KPL-4 breast cancer xenograft; non-tumor-bearing mice for PK-only)",
    n_subjects     = NA_integer_,
    n_studies      = 5L,
    disease_state  = "HER2-expressing human breast cancer cell-line xenograft tumors in immunocompromised mice",
    dose_range     = "3 mg/kg IV bolus (non-tumor-bearing mouse PK, Erickson 2012 / Okamoto 2020); 300 ug/kg DM1-equivalent IV (tumor xenograft, Erickson 2012); PD studies at multiple doses (Table S2c cell-line-specific fits)",
    regions        = "Preclinical (Scheuher et al. 2023 fits + literature-digitized data)",
    notes          = "Cell-line-specific PD parameters (Table S2c): BT-474 T-DM1, BT-474EEI T-DM1, N87 T-DM1 + T-DXd, KPL-4 T-DM1. Default parameterization: N87 with T-DM1 (kkill_max=0.31/day, tau=0.25 day, kc50=485 nM, n_Hill=1, tdouble=12.37 day, klin=189.56 mm3/day). Trastuzumab does NOT bind rodent HER2, so this model has no systemic HER2 sink or soluble HER2 term (see human model for those extensions). See vignette Errata for deviations from paper's per-stage tumor-cell tracking."
  )

  ini({
    # -------------------------------------------------------------------
    # Structural parameters (T-DM1 in N87 xenograft, mouse)
    # Source: Scheuher 2023 Supplement Tables S2b (mouse PK), S2c (mouse PD, N87 with T-DM1)
    # All rates converted from 1/s to 1/h by *3600 for internal consistency.
    # -------------------------------------------------------------------

    # ADC deconjugation rate (T-DM1 in mouse)
    lkdec <- log(8.5e-7 * 3600); label("Deconjugation rate constant kdec (1/h)")   # Table S2b: 8.5e-7 /s = 3.06e-3 /h (T-DM1); fit to Erickson 2012

    # DAR (drug-to-antibody ratio) for T-DM1
    lDAR <- fixed(log(3.5)); label("Drug-to-antibody ratio DAR (T-DM1)")               # Table S2b: 3.5 for T-DM1

    # HER2 binding (association rate + KD)
    lkonab <- fixed(log(1e-4 * 3600)); label("HER2 binding association rate k_on_Ab (1/nM/h)")  # Table S2b: 1e-4 /nM/s = 0.36 /nM/h
    lkdab  <- log(0.3); label("HER2:ADC equilibrium binding constant K_D_Ab (nM)")     # Table S2b: 0.3 nM (carried from in vitro)

    # HER2 receptor kinetics (from in vitro fit, carried to mouse per Table S2b header)
    lkendoher2 <- log(4.27e-5 * 3600); label("HER2 endocytosis rate constant (1/h)")  # Table S2a: 4.27e-5 /s (via S2b carry-over)
    lkrecher2  <- log(2.4e-5 * 3600);  label("HER2 recycling rate constant (1/h)")    # Table S2a: 2.4e-5 /s
    lkdegher2  <- log(1.27e-4 * 3600); label("HER2 endosomal degradation (1/h)")      # Table S2a: 1.27e-4 /s

    # Endosomal linker cleavage (0 for T-DM1 non-cleavable linker)
    lkcleave <- fixed(log(1e-12)); label("Endosomal linker cleavage rate (1/h)")    # Table S2a: 0 /s (T-DM1)

    # Payload target binding (DM1:tubulin)
    lkonpl <- fixed(log(1e-3 * 3600)); label("Payload:target association (1/nM/h)")  # Table S2b: 1e-3 /nM/s
    lkdpl  <- log(930); label("Payload:target K_D_PL (nM)")                             # Table S2b: 930 nM (DM1:tubulin)

    # Payload cellular flux (T-DM1)
    lkinpl  <- log(5.95e-5 * 3600); label("Payload influx k_in_PL (1/h)")           # Table S2b: 5.95e-5 /s (T-DM1)
    lkoutpl <- log(3.95e-5 * 3600); label("Payload efflux k_out_PL (1/h)")          # Table S2b: 3.95e-5 /s (T-DM1)

    # Payload partition coefficient in tumor tissue
    lPcpl <- log(0.51); label("Payload tumor partition coefficient Pc_PL (unitless)")   # Table S2b: 0.51 (T-DM1)

    # Cell / target constants
    lVcell    <- fixed(log(2e-12));  label("Volume of a cell (L)")                     # Table S2b: 2e-12 L typical
    lRPCher2  <- fixed(log(1e6));    label("HER2 receptors per tumor cell")            # Table S2b: 1e6 (N87 approx)
    lTpercell <- log(65); label("Intracellular target T_per_cell (nM) - tubulin for DM1")  # Table S2b: 65 nM

    # Molecular weights (g/mol)
    lMWadc <- fixed(log(148781)); label("Molecular weight of ADC (g/mol, T-DM1)")       # Table S2b: 148781 (T-DM1)
    lMWab  <- fixed(log(145167)); label("Molecular weight of naked antibody (g/mol)")   # Table S2b: 145167
    lMWpl  <- fixed(log(1100));   label("Molecular weight of payload (g/mol, DM1)")     # Table S2b: 1100 (Lys-MCC-DM1)

    # Mouse physiological volumes and body weight
    lVcentral    <- fixed(log(1e-3));  label("Central compartment volume V_c (L)")     # Table S2b: 1e-3 L
    lVperipheral <- fixed(log(5.7e-3)); label("Peripheral compartment volume V_p (L)") # Table S2b: 5.7e-3 L
    lBW          <- fixed(log(0.02));  label("Body weight (kg)")                       # Table S2b: 2e-2 kg

    # Systemic elimination half-lives (converted to hours)
    lthalf_ab <- log(11.6 * 24); label("Antibody / ADC elimination half-life t_(1/2),Ab (hour)")  # Table S2b: 11.6 day = 278.4 hour
    lthalf_pl <- log(3);         label("Payload elimination half-life t_(1/2),PL (hour, T-DM1)")  # Table S2b: 3 hour

    # Distribution parameters (central <-> peripheral)
    ltdist12_ab <- log(2.65); label("Antibody/ADC central-to-peripheral distribution time t_dist12,Ab (hour)")  # Table S2b: 2.65 hour
    lPdist12_ab <- log(0.34); label("Antibody/ADC partition coefficient P_dist12,Ab (unitless)")               # Table S2b: 0.34

    ltdist12_pl <- log(0.1);   label("Payload central-to-peripheral distribution time t_dist12,PL (hour)")     # Table S2b: 0.1 hour
    lPdist12_pl <- log(83.02); label("Payload partition coefficient P_dist12,PL (unitless)")                    # Table S2b: 83.02

    # Krogh cylinder / vascular tumor transport parameters
    lPab <- log(50.34);  label("Vascular permeability of Ab/ADC P_Ab (um/day)")        # Table S2b: 50.34
    lDab <- log(1.74e5); label("Diffusivity of Ab/ADC D_Ab (um^2/day)")                 # Table S2b: 1.74e5
    lPpl <- log(924);    label("Vascular permeability of payload P_PL (um/day)")        # Table S2b: 924
    lDpl <- log(1.1e7);  label("Diffusivity of payload D_PL (um^2/day)")                # Table S2b: 1.1e7

    lRcap   <- fixed(log(8));    label("Capillary radius R_cap (um)")                  # Table S2b: 8 um
    lRkrogh <- fixed(log(75));   label("Krogh cylinder radius R_Krogh (um)")           # Table S2b: 75 um
    leps_ab <- fixed(log(0.24)); label("Tumor void fraction for Ab/ADC epsilon_Ab")    # Table S2b: 0.24
    leps_pl <- fixed(log(0.44)); label("Tumor void fraction for payload epsilon_PL")   # Table S2b: 0.44

    # -------------------------------------------------------------------
    # PD parameters (N87 T-DM1 default, Table S2c)
    # -------------------------------------------------------------------
    lpsi <- fixed(log(20)); label("Growth-phase-switch parameter psi")                 # Table S2c: 20
    lVtumor_max <- fixed(log(5e-3)); label("Maximum tumor volume V_tumor_max (L)")     # Table S2c: 5e-3 L (assumed)
    lVtumor_i <- log(2e-4); label("Initial tumor volume V_tumor_i (L)")                # 200 mm^3 typical mouse xenograft start (2e-4 L = 200 mm^3)

    lkkillmax <- log(0.31 / 24); label("Maximum tumor kill rate kkill_max (1/h) - N87+T-DM1")  # Table S2c N87+T-DM1: 0.31 /day = 0.01292 /h
    ltau      <- log(0.25 * 24); label("Kill transduction delay tau (hour) - N87+T-DM1")           # Table S2c: 0.25 day = 6 hour
    lkc50     <- log(485); label("Kill half-max concentration kc50 (nM) - N87+T-DM1")               # Table S2c: 485 nM
    lnHill    <- fixed(log(1)); label("Hill coefficient n_Hill - N87+T-DM1")                        # Table S2c: 1

    ltdouble <- log(12.37 * 24); label("Exponential tumor doubling time t_double (hour) - N87")   # Table S2c: 12.37 day = 296.88 hour
    lklin    <- log(189.56 / 24); label("Linear tumor growth rate k_lin (mm^3/h) - N87")        # Table S2c: 189.56 mm^3/day = 7.898 mm^3/h

    # -------------------------------------------------------------------
    # Residual error (placeholder - paper reports deterministic fits)
    # -------------------------------------------------------------------
    propSd <- fixed(0.15); label("Proportional residual error (fraction) - not from paper; 15% for simulation")  # Assumption: simulation-only default
  })

  model({
    # -------------------------------------------------------------------
    # 1. Parameter transformations
    # -------------------------------------------------------------------
    kdec       <- exp(lkdec)
    DAR        <- exp(lDAR)
    konab      <- exp(lkonab)
    kdab       <- exp(lkdab)
    koffab     <- konab * kdab                       # k_off,Ab = k_on,Ab * K_D,Ab (Table S2a inferred)
    kendoher2  <- exp(lkendoher2)
    krecher2   <- exp(lkrecher2)
    kdegher2   <- exp(lkdegher2)
    kcleave    <- exp(lkcleave)
    konpl      <- exp(lkonpl)
    kdpl       <- exp(lkdpl)
    koffpl     <- konpl * kdpl
    kinpl      <- exp(lkinpl)
    koutpl     <- exp(lkoutpl)
    Pcpl       <- exp(lPcpl)
    Vcell      <- exp(lVcell)
    RPCher2    <- exp(lRPCher2)
    Tpercell   <- exp(lTpercell)
    MWadc      <- exp(lMWadc)
    MWab       <- exp(lMWab)
    MWpl       <- exp(lMWpl)
    Vcentral   <- exp(lVcentral)
    Vperi      <- exp(lVperipheral)
    BW         <- exp(lBW)
    thalf_ab   <- exp(lthalf_ab)
    thalf_pl   <- exp(lthalf_pl)
    tdist12_ab <- exp(ltdist12_ab)
    Pdist12_ab <- exp(lPdist12_ab)
    tdist12_pl <- exp(ltdist12_pl)
    Pdist12_pl <- exp(lPdist12_pl)
    Pab        <- exp(lPab)     # um/day
    Dab        <- exp(lDab)     # um^2/day
    Ppl        <- exp(lPpl)
    Dpl        <- exp(lDpl)
    Rcap       <- exp(lRcap)
    Rkrogh     <- exp(lRkrogh)
    eps_ab     <- exp(leps_ab)
    eps_pl     <- exp(leps_pl)

    psi        <- exp(lpsi)
    Vtumor_max <- exp(lVtumor_max)
    Vtumor_i   <- exp(lVtumor_i)

    kkillmax <- exp(lkkillmax)
    tau_kill <- exp(ltau)
    kc50     <- exp(lkc50)
    nHill    <- exp(lnHill)
    tdouble  <- exp(ltdouble)
    klin_h   <- exp(lklin)      # mm^3/h

    # -------------------------------------------------------------------
    # 2. Derived / inferred quantities (Table S2 inferred parameters)
    # -------------------------------------------------------------------
    # Systemic elimination rates
    kelim_ab <- log(2) / thalf_ab
    kelim_pl <- log(2) / thalf_pl

    # Distribution rates central <-> peripheral (Table S2b inferred k_12/k_21)
    k12_ab <- (log(2) / tdist12_ab) * Pdist12_ab / (Pdist12_ab + Vcentral / Vperi)
    k21_ab <- (log(2) / tdist12_ab) * (Vcentral / Vperi) / (Pdist12_ab + Vcentral / Vperi)
    k12_pl <- (log(2) / tdist12_pl) * Pdist12_pl / (Pdist12_pl + Vcentral / Vperi)
    k21_pl <- (log(2) / tdist12_pl) * (Vcentral / Vperi) / (Pdist12_pl + Vcentral / Vperi)

    # Convert diffusion / permeability from day to hour (Krogh cylinder)
    Pab_h <- Pab / 24
    Dab_h <- Dab / 24
    Ppl_h <- Ppl / 24
    Dpl_h <- Dpl / 24

    # HER2 amount per cell at steady state (Table S2 inferred)
    # RPCher2 receptors * (1 mol / 6.022e23) * 1e9 nmol/mol
    HER2_0_per_cell_nmol <- RPCher2 / 6.022e23 * 1e9
    HER2_endo_0_per_cell_nmol <- kendoher2 / (krecher2 + kdegher2) * HER2_0_per_cell_nmol

    # Zeroth-order HER2 synthesis rate per cell (nmol/h per cell)
    ksynth_per_cell <- HER2_0_per_cell_nmol * kendoher2 * kdegher2 / (kdegher2 + krecher2)

    # Total living tumor cells and tumor volumes (updated dynamically from N states)
    Ntot <- n1 + n2 + n3 + n4
    # Guard against Ntot=0 during warmup
    Ntot_safe <- Ntot + 1e-12
    Vtumor <- (Ntot_safe * Vcell) / (1 - eps_pl)
    Vext_tumor_ab <- eps_ab * Vtumor    # Volume accessible to Ab/ADC
    Vext_tumor_pl <- eps_pl * Vtumor    # Volume accessible to payload
    Vint_tumor    <- (1 - eps_pl) * Vtumor  # Intracellular tumor volume

    # Krogh cylinder / surface-exchange tumor uptake rate constants
    # k_13,X = (2*P_X*Rcap/Rkrogh^2 + 6*D_X*Rcap/Vtumor^(2/3)) * Vtumor/Vcentral
    # k_31,X = (2*P_X*Rcap/Rkrogh^2 + 6*D_X*Rcap/Vtumor^(2/3)) * 1/eps_X
    # Note Vtumor in L must be converted for Rcap/Rkrogh in um to align units.
    # Simplification: use paper's dimensionless factor and treat Vtumor in mm^3 for
    # geometric terms (1 L = 1e6 mm^3, so Vtumor_mm3 = Vtumor * 1e6).
    Vtumor_mm3 <- Vtumor * 1e6

    krogh_ab <- 2 * Pab_h * Rcap / (Rkrogh^2) + 6 * Dab_h * Rcap / (Vtumor_mm3^(2/3) * 1000)
    krogh_pl <- 2 * Ppl_h * Rcap / (Rkrogh^2) + 6 * Dpl_h * Rcap / (Vtumor_mm3^(2/3) * 1000)

    k13_ab <- krogh_ab * Vtumor / Vcentral
    k31_ab <- krogh_ab / eps_ab
    k13_pl <- krogh_pl * Vtumor / Vcentral
    k31_pl <- krogh_pl / eps_pl

    # -------------------------------------------------------------------
    # 3. Reactions and rate laws (Table S3c; aggregated tumor-interior form)
    # See vignette Errata for the aggregation vs paper per-stage tracking.
    # -------------------------------------------------------------------

    # Systemic PK (Table S3c reactions v14-v27)
    v_dec_c   <- kdec * adc_central                  # v14 central deconjugation
    v_elim_adc_c <- kelim_ab * adc_central           # v15 ADC clearance
    v_elim_ab_c  <- kelim_ab * ab_central            # v16 Ab clearance
    v_elim_pl_c  <- kelim_pl * pl_central            # v17 PL clearance
    v_c2p_adc <- k12_ab * adc_central; v_p2c_adc <- k21_ab * adc_peripheral   # v18
    v_c2p_ab  <- k12_ab * ab_central;  v_p2c_ab  <- k21_ab * ab_peripheral    # v19
    v_c2p_pl  <- k12_pl * pl_central;  v_p2c_pl  <- k21_pl * pl_peripheral    # v20

    v_dec_p   <- kdec * adc_peripheral               # v21
    v_elim_adc_p <- kelim_ab * adc_peripheral        # v22
    v_elim_ab_p  <- kelim_ab * ab_peripheral         # v23
    v_elim_pl_p  <- kelim_pl * pl_peripheral         # v24

    # Central <-> tumor extracellular (Krogh cylinder)
    v_c2t_adc <- k13_ab * adc_central; v_t2c_adc <- k31_ab * adc_ext_tumor    # v25
    v_c2t_ab  <- k13_ab * ab_central;  v_t2c_ab  <- k31_ab * ab_ext_tumor     # v26
    v_c2t_pl  <- k13_pl * pl_central;  v_t2c_pl  <- k31_pl * pl_ext_tumor     # v27

    v_dec_t <- kdec * adc_ext_tumor                  # v28 tumor extracellular deconjugation

    # Tumor cell interior HER2/ADC dynamics (aggregate across N1+N2+N3+N4).
    # v29: HER2 synthesis scaled by total living cells
    v_synth <- ksynth_per_cell * Ntot_safe

    # v30: HER2 (tumor cell surface) + ADC(ext) binding
    v30f <- konab / Vext_tumor_ab * her2_surf_tumor * adc_ext_tumor
    v30r <- koffab * her2_adc_surf_tumor

    # v31: HER2 (tumor cell surface) + Ab(ext) binding
    v31f <- konab / Vext_tumor_ab * her2_surf_tumor * ab_ext_tumor
    v31r <- koffab * her2_ab_surf_tumor              # v31r typo fix per Errata (paper printed k_on,PL K_D,PL)

    # v32: cell-surface HER2:ADC deconjugation (in tumor cells)
    v32 <- kdec * her2_adc_surf_tumor

    # v33: HER2 endocytosis / recycling
    v33f <- kendoher2 * her2_surf_tumor
    v33r <- krecher2  * her2_endo_tumor

    # v34: HER2:ADC endocytosis / recycling
    v34f <- kendoher2 * her2_adc_surf_tumor
    v34r <- krecher2  * her2_adc_endo_tumor

    # v35: HER2:Ab endocytosis / recycling
    v35f <- kendoher2 * her2_ab_surf_tumor
    v35r <- krecher2  * her2_ab_endo_tumor

    # v36: endosomal HER2 degradation
    v36 <- kdegher2 * her2_endo_tumor

    # v37: endosomal HER2:ADC degradation -> DAR * payload_endo
    v37 <- kdegher2 * her2_adc_endo_tumor

    # v38: endosomal HER2:Ab degradation
    v38 <- kdegher2 * her2_ab_endo_tumor

    # v39: endosomal linker cleavage HER2:ADC -> HER2:Ab + DAR*payload
    v39 <- kcleave * her2_adc_endo_tumor

    # v40: payload endosome -> cytosol
    v40 <- kinpl * pl_endo_tumor

    # v41: cytosolic target binding
    v41f <- konpl / (Pcpl * Vint_tumor) * t_cyto_tumor * pl_cyto_tumor
    v41r <- koffpl * tpl_cyto_tumor

    # v79: payload cytosol <-> tumor extracellular
    v79f <- kinpl  * pl_ext_tumor                    # into cell (extracellular concentration used)
    v79r <- koutpl * pl_cyto_tumor                   # out of cell

    # -------------------------------------------------------------------
    # 4. Tumor growth and kill dynamics (Simeoni transit chain)
    # -------------------------------------------------------------------
    # Free intracellular unconjugated payload concentration (nM)
    Cpl_intra <- (pl_cyto_tumor + tpl_cyto_tumor) / Vint_tumor

    # Hill kill rate on proliferating N1 cells
    kkill <- kkillmax * Cpl_intra^nHill / (kc50^nHill + Cpl_intra^nHill)

    # Simeoni-style growth rate: exp for small tumors, transitioning to linear
    Nmax <- (1 - eps_pl) * Vtumor_max / Vcell
    # w = exp(1/psi) * log(2 * Ntot * Vcell / (tdouble * (1-eps_pl) * klin_h * 1e6)) style
    # Use the paper's inferred k_growth formulation
    growth_num <- (log(2) / tdouble) * (1 - Ntot_safe / Nmax)
    growth_denom_arg <- log(2) * Ntot_safe * Vcell / (tdouble * (1 - eps_pl) * klin_h * 1e-6)
    # Ensure positive base for the psi-th root
    growth_denom_pos <- 1 + max(growth_denom_arg, 0)^psi
    k_growth <- growth_num / (growth_denom_pos^(1/psi))

    # Cell death outflow fraction (for aggregate species leaving tumor when N4 -> death)
    death_frac <- (n4 / Ntot_safe) / tau_kill

    # -------------------------------------------------------------------
    # 5. ODE system (Table S3d aggregate form)
    # -------------------------------------------------------------------

    # Systemic (Table S3d Eqns 14-19)
    d/dt(adc_central) <- -v_dec_c - v_elim_adc_c - v_c2p_adc + v_p2c_adc - v_c2t_adc + v_t2c_adc
    d/dt(ab_central)  <-  v_dec_c - v_elim_ab_c  - v_c2p_ab  + v_p2c_ab  - v_c2t_ab  + v_t2c_ab
    d/dt(pl_central)  <-  v_dec_c * DAR + v_elim_adc_c * DAR - v_elim_pl_c - v_c2p_pl + v_p2c_pl - v_c2t_pl + v_t2c_pl
    d/dt(adc_peripheral) <- v_c2p_adc - v_p2c_adc - v_dec_p - v_elim_adc_p
    d/dt(ab_peripheral)  <- v_c2p_ab  - v_p2c_ab  + v_dec_p - v_elim_ab_p
    d/dt(pl_peripheral)  <- v_c2p_pl  - v_p2c_pl  + v_dec_p * DAR + v_elim_adc_p * DAR - v_elim_pl_p

    # Tumor extracellular (Table S3d Eqns 20-22, aggregate)
    d/dt(adc_ext_tumor) <- v_c2t_adc - v_t2c_adc - v_dec_t - v30f + v30r
    d/dt(ab_ext_tumor)  <- v_c2t_ab  - v_t2c_ab  + v_dec_t - v31f + v31r
    d/dt(pl_ext_tumor)  <- v_c2t_pl  - v_t2c_pl  + v_dec_t * DAR + v32 * DAR + v79r * Ntot_safe / Ntot_safe - v79f
    # Note: v79f/v79r are aggregate; paper sums over N1-N4 for extracellular <-> cytosolic exchange

    # Tumor cell staging (Simeoni transit chain, N1-N4)
    # v42: N1 proliferation, v46: N1 -> N2 kill, v47: Ns/tau, v48: N4/tau -> death
    d/dt(n1) <- k_growth * n1 - kkill * n1                       # (23) Table S3d
    d/dt(n2) <- kkill * n1 - n2 / tau_kill                       # (24)
    d/dt(n3) <- n2 / tau_kill - n3 / tau_kill                    # (25)
    d/dt(n4) <- n3 / tau_kill - n4 / tau_kill                    # (26)

    # Tumor cell interior species (aggregate across all 4 stages).
    # Sources: v29 synthesis (only in growing N1) + v33-v41 processing;
    # sinks: reactions + cell-death outflow proportional to death_frac.
    #
    # HER2 cell surface (v29 synth into N1; v33 endo; +recycling; binding v30/v31 sinks; death loss)
    d/dt(her2_surf_tumor) <- v_synth - v33f + v33r - v30f + v30r - v31f + v31r - her2_surf_tumor * death_frac

    # HER2:ADC cell surface (v30 binding, v34 endo, v32 deconjugation to HER2:Ab, death loss)
    d/dt(her2_adc_surf_tumor) <- v30f - v30r - v34f + v34r - v32 - her2_adc_surf_tumor * death_frac

    # HER2:Ab cell surface (v31 binding, v32 deconjugation, v35 endo, death loss)
    d/dt(her2_ab_surf_tumor) <- v31f - v31r + v32 - v35f + v35r - her2_ab_surf_tumor * death_frac

    # HER2 endosomal
    d/dt(her2_endo_tumor) <- v33f - v33r - v36 - her2_endo_tumor * death_frac

    # HER2:ADC endosomal (v34 endo in, v37 degradation, v39 cleavage, death loss)
    d/dt(her2_adc_endo_tumor) <- v34f - v34r - v37 - v39 - her2_adc_endo_tumor * death_frac

    # HER2:Ab endosomal (v35 endo in, v38 degradation, +v39 cleavage produces HER2:Ab_endo, death loss)
    d/dt(her2_ab_endo_tumor) <- v35f - v35r - v38 + v39 - her2_ab_endo_tumor * death_frac

    # Payload endosomal (from HER2:ADC degradation v37 and cleavage v39, x DAR; -v40 out to cytosol; death loss returns to ext)
    d/dt(pl_endo_tumor) <- (v37 + v39) * DAR - v40 - pl_endo_tumor * death_frac

    # Payload cytosolic (v40 from endosome, v41 target binding, v79f/v79r extracellular exchange, death loss)
    d/dt(pl_cyto_tumor) <- v40 - v41f + v41r + v79f - v79r - pl_cyto_tumor * death_frac

    # Cytosolic target (tubulin/TOPO-1): steady maintenance in living cells + turnover balance
    # v45 synthesis proportional to k_growth*N1; v41 binding; death loss
    v_target_synth <- k_growth * n1 * Vcell * Tpercell
    d/dt(t_cyto_tumor) <- v_target_synth - v41f + v41r - t_cyto_tumor * death_frac

    # Target:payload complex
    d/dt(tpl_cyto_tumor) <- v41f - v41r - tpl_cyto_tumor * death_frac

    # -------------------------------------------------------------------
    # 6. Initial conditions
    # -------------------------------------------------------------------
    # Initial tumor cells and steady-state HER2
    N_i <- (1 - eps_pl) * Vtumor_i / Vcell
    n1(0) <- N_i
    n2(0) <- 0
    n3(0) <- 0
    n4(0) <- 0
    her2_surf_tumor(0) <- HER2_0_per_cell_nmol * N_i
    her2_endo_tumor(0) <- HER2_endo_0_per_cell_nmol * N_i
    t_cyto_tumor(0)    <- Tpercell * Vcell * N_i

    # -------------------------------------------------------------------
    # 7. Observations (plasma concentrations in nM)
    # -------------------------------------------------------------------
    Cc         <- adc_central / Vcentral                     # ADC plasma conc (nM)
    Cab        <- ab_central  / Vcentral                     # total antibody plasma conc (nM)
    Cpl_plasma <- pl_central  / Vcentral                     # free payload plasma conc (nM)
    Vtumor_out <- Vtumor_mm3                                  # tumor volume (mm^3)
    Cpl_intra_out <- Cpl_intra                                # intracellular free payload (nM)

    Cc ~ prop(propSd)
  })
}
