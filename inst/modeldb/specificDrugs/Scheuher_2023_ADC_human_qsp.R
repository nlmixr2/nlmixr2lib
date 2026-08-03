Scheuher_2023_ADC_human_qsp <- function() {
  description <- "QSP. Human platform model for HER2-targeting antibody-drug conjugates in HER2+ metastatic breast cancer (T-DM1 default; T-DXd variant via parameter overrides). Extends the mouse model with: (i) HER2 receptor sinks on normal cells in the central and peripheral compartments (with binding, endocytosis, recycling, degradation); (ii) soluble HER2 (sHER2) shed from cell-surface HER2 into central + peripheral + tumor compartments, with reversible binding to ADC and Ab and its own turnover; and (iii) larger physiologic volumes (3 L central, 13 L peripheral for a 70 kg adult). Mouse-derived TGI parameters (kkill_max, kc50, tau, n_Hill) are carried over from N87 xenograft fits. Amounts in nmol; concentrations amount/volume."
  reference   <- "Scheuher B, Ghusinga KR, McGirr K, Nowak M, Panday S, Apgar J, Subramanian K, Betts A. Towards a platform quantitative systems pharmacology (QSP) model for preclinical to clinical translation of antibody drug conjugates (ADCs). J Pharmacokinet Pharmacodyn. 2023;51(1):5-30. doi:10.1007/s10928-023-09884-6. Human model = Tables S1c, S2d-e, S3e-f."
  vignette    <- "Scheuher_2023_ADC_platform_qsp"
  units       <- list(time = "h", dosing = "nmol", concentration = "nM")

  covariateData <- list()

  covariatesDataExcluded <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 4L,
    age_range      = "adult HER2+ metastatic breast cancer patients (simulated virtual cohorts)",
    weight_range   = "70 kg reference adult",
    disease_state  = "HER2+ (HER2 1+ = 2e4 receptors/cell; HER2 3+ = 1e6 receptors/cell) metastatic breast cancer",
    dose_range     = "T-DM1: 0.3-4.8 mg/kg IV Q3W; T-DXd: 0.8-8.0 mg/kg IV Q3W; MTD dosing 3.6 mg/kg (T-DM1) or 5.4 mg/kg (T-DXd) Q3W over 14-29 months",
    regions        = "Global clinical trials (phase 1 dose escalation + phase 2 + phase 3)",
    notes          = "Human parameters (Table S2d) inherit systemic PK from mouse with volume rescaling to human physiology; PD parameters (kkill_max, kc50, tau, n_Hill, n_Hill from N87) carried directly from mouse T-DM1 or T-DXd fits. Soluble HER2 baseline (8 ng/mL healthy) and shedding calibrated to human data. Virtual patient CVs from Table S2e (kkill_max CV 10%, kc50 CV 2000%/10000% for T-DM1/T-DXd). See vignette Errata for aggregation of the tumor-cell 4-stage cascade."
  )

  ini({
    # -------------------------------------------------------------------
    # Systemic PK parameters (T-DM1 in humans, Table S2d)
    # -------------------------------------------------------------------
    lkdec <- log(8.5e-7 * 3600); label("Deconjugation rate constant kdec (1/h) - T-DM1")   # Table S2d: 8.5e-7 /s (carried from mouse)
    lDAR  <- fixed(log(3.5)); label("Drug-to-antibody ratio DAR (T-DM1)")                     # Table S2b: 3.5

    # HER2 binding (identical to in vitro / mouse)
    lkonab <- fixed(log(1e-4 * 3600)); label("HER2 binding association rate (1/nM/h)")
    lkdab  <- log(0.3); label("HER2:ADC K_D_Ab (nM)")

    # HER2 receptor kinetics carried from mouse
    lkendoher2 <- log(4.27e-5 * 3600); label("HER2 endocytosis (1/h)")
    lkrecher2  <- log(2.4e-5  * 3600); label("HER2 recycling (1/h)")
    lkdegher2  <- log(1.27e-4 * 3600); label("HER2 endosomal degradation (1/h)")

    lkcleave <- fixed(log(1e-12)); label("Endosomal linker cleavage (1/h, T-DM1)")

    # Payload binding
    lkonpl <- fixed(log(1e-3 * 3600)); label("Payload:target association (1/nM/h)")
    lkdpl  <- log(930); label("Payload:target K_D_PL (nM, DM1)")
    lkinpl <- log(5.95e-5 * 3600); label("Payload influx (1/h)")
    lkoutpl <- log(3.95e-5 * 3600); label("Payload efflux (1/h)")
    lPcpl <- log(0.51); label("Payload tumor partition Pc_PL")

    # Cell / target constants
    lVcell    <- fixed(log(2e-12));  label("Cell volume (L)")                     # Table S2d: 2e-12 L
    lRPCher2  <- fixed(log(1e6));    label("HER2 receptors per tumor cell")       # Table S2d: 1e6 (HER2 3+)
    lRPCnormal <- fixed(log(2e4));   label("HER2 receptors per normal cell")      # Table S2d: 2e4
    lTpercell <- log(65); label("Target T_per_cell (nM tubulin for DM1)")

    # Human physiologic volumes
    lVcentral    <- fixed(log(3));   label("Human central volume V_c (L)")        # Table S2d: 3 L
    lVperipheral <- fixed(log(13));  label("Human peripheral volume V_p (L)")     # Table S2d: 13 L
    lBW          <- fixed(log(70));  label("Body weight (kg)")                    # Table S2d: 70 kg

    # Normal cell counts (calibrated to human PK)
    lNcellCentral <- fixed(log(2.74e8)); label("HER2+ normal cells in central compartment") # Table S2d: 2.74e8
    lNcellPeri    <- fixed(log(2.08e11)); label("HER2+ normal cells in peripheral compartment") # Table S2d: 2.08e11

    # Soluble HER2 physiology
    lCsHER2_c <- log(8);    label("Baseline sHER2 concentration in central (ng/mL healthy)")  # Table S2d: 8 ng/mL (Jensen 2003)
    lchi_shed_ab <- fixed(log(1e-12)); label("Fraction shedding when HER2 bound to Ab chi_shed")  # Table S2d: 0 (Perrier 2018)
    lMWsHER2 <- fixed(log(1e5)); label("Soluble HER2 molecular weight (g/mol)")
    lthalf_sher2 <- log(5); label("sHER2 half-life (hour)")                                    # Table S2d: 5 hour
    lthalf_sher2_ab <- log(11.6 * 24); label("sHER2:ADC/Ab complex half-life (hour)")          # Table S2d: 11.6 day

    # sHER2 distribution
    ltdist12_sher2 <- log(8);  label("sHER2 central-peripheral distribution time (hour)")     # Table S2d: 8 hour
    lPdist12_sher2 <- fixed(log(1)); label("sHER2 partition coefficient (unitless)")          # Table S2d: 1
    ltdist13_sher2 <- log(8);  label("sHER2 central-tumor distribution time (hour)")          # Table S2d: 8 hour
    lPdist13_sher2 <- fixed(log(0.24)); label("sHER2 tumor partition coefficient")            # Table S2d: 0.24

    # Systemic elimination half-lives
    lthalf_ab <- log(11.6 * 24); label("ADC/Ab elimination half-life (hour)")
    lthalf_pl <- log(3);         label("Payload half-life (hour, DM1)")

    # Distribution
    ltdist12_ab <- log(10);    label("Ab/ADC distribution time (hour, human)")     # Table S2d: 10 hour (calibrated)
    lPdist12_ab <- log(0.187); label("Ab/ADC partition coefficient (human typical)") # Table S2d: 0.187
    ltdist12_pl <- log(0.1);   label("Payload distribution time (hour)")
    lPdist12_pl <- log(83.02); label("Payload partition coefficient")

    # Krogh / tumor uptake
    lPab <- log(50.34);  label("Vascular permeability of Ab/ADC (um/day)")
    lDab <- log(1.74e5); label("Diffusivity of Ab/ADC (um^2/day)")
    lPpl <- log(924);    label("Vascular permeability of payload (um/day)")
    lDpl <- log(1.1e7);  label("Diffusivity of payload (um^2/day)")
    lRcap   <- fixed(log(8));    label("Capillary radius (um)")
    lRkrogh <- fixed(log(75));   label("Krogh cylinder radius (um)")
    leps_ab <- fixed(log(0.24)); label("Tumor void fraction, Ab")
    leps_pl <- fixed(log(0.44)); label("Tumor void fraction, PL")

    # -------------------------------------------------------------------
    # PD parameters (Table S2d + S2e: human HER2+ MBC, using N87 T-DM1 defaults)
    # -------------------------------------------------------------------
    lpsi <- fixed(log(20)); label("Growth-phase switch psi")
    lVtumor_max <- fixed(log(0.5238));  label("Max human tumor volume V_tumor_max (L)")  # Table S2d: 0.5238 L
    lVtumor_i   <- log(0.0016); label("Initial human tumor volume V_tumor_i (L)")         # Table S2d: 0.0016 L

    ltdouble <- log(25 * 24); label("Human tumor doubling time (hour)")                   # Table S2d: 25 day
    lklin    <- log(621 / 24); label("Human linear tumor growth k_lin (mm^3/h)")       # Table S2d: 621 mm^3/day

    # Kill parameters from mouse N87 T-DM1 (Table S2c)
    lkkillmax <- log(0.139 / 24); label("Kill max k_kill_max (1/h, T-DM1 mean)")       # Table S2e: mean 1.39e-1 /day
    ltau      <- log(0.25 * 24); label("Kill delay tau (hour, N87)")                       # Table S2c: 0.25 day (N87 T-DM1)
    lkc50     <- log(23.8); label("Kill kc_50 (nM, T-DM1 mean)")                           # Table S2e: 23.8 nM (calibrated)
    lnHill    <- fixed(log(1)); label("Hill coefficient (N87 T-DM1)")                      # Table S2c: 1

    # Residual error placeholder
    propSd <- fixed(0.15); label("Proportional residual error (assumption; paper deterministic)")
  })

  model({
    # -------------------------------------------------------------------
    # Parameter transformations
    # -------------------------------------------------------------------
    kdec       <- exp(lkdec)
    DAR        <- exp(lDAR)
    konab      <- exp(lkonab); kdab <- exp(lkdab); koffab <- konab * kdab
    kendoher2  <- exp(lkendoher2); krecher2 <- exp(lkrecher2); kdegher2 <- exp(lkdegher2)
    kcleave    <- exp(lkcleave)
    konpl      <- exp(lkonpl); kdpl <- exp(lkdpl); koffpl <- konpl * kdpl
    kinpl      <- exp(lkinpl); koutpl <- exp(lkoutpl); Pcpl <- exp(lPcpl)
    Vcell      <- exp(lVcell)
    RPCher2    <- exp(lRPCher2); RPCnormal <- exp(lRPCnormal)
    Tpercell   <- exp(lTpercell)
    Vcentral   <- exp(lVcentral); Vperi <- exp(lVperipheral); BW <- exp(lBW)
    NcellCentral <- exp(lNcellCentral); NcellPeri <- exp(lNcellPeri)

    CsHER2_c    <- exp(lCsHER2_c)
    chi_shed_ab <- exp(lchi_shed_ab)
    MWsHER2     <- exp(lMWsHER2)
    thalf_sher2 <- exp(lthalf_sher2)
    thalf_sher2_ab <- exp(lthalf_sher2_ab)

    tdist12_sher2 <- exp(ltdist12_sher2); Pdist12_sher2 <- exp(lPdist12_sher2)
    tdist13_sher2 <- exp(ltdist13_sher2); Pdist13_sher2 <- exp(lPdist13_sher2)

    thalf_ab <- exp(lthalf_ab); thalf_pl <- exp(lthalf_pl)
    tdist12_ab <- exp(ltdist12_ab); Pdist12_ab <- exp(lPdist12_ab)
    tdist12_pl <- exp(ltdist12_pl); Pdist12_pl <- exp(lPdist12_pl)

    Pab <- exp(lPab); Dab <- exp(lDab); Ppl <- exp(lPpl); Dpl <- exp(lDpl)
    Rcap <- exp(lRcap); Rkrogh <- exp(lRkrogh); eps_ab <- exp(leps_ab); eps_pl <- exp(leps_pl)

    psi <- exp(lpsi); Vtumor_max <- exp(lVtumor_max); Vtumor_i <- exp(lVtumor_i)
    tdouble <- exp(ltdouble); klin_h <- exp(lklin)
    kkillmax <- exp(lkkillmax); tau_kill <- exp(ltau); kc50 <- exp(lkc50); nHill <- exp(lnHill)

    # -------------------------------------------------------------------
    # Derived quantities
    # -------------------------------------------------------------------
    kelim_ab <- log(2) / thalf_ab
    kelim_pl <- log(2) / thalf_pl
    kelim_sher2 <- log(2) / thalf_sher2
    kelim_sher2_ab <- log(2) / thalf_sher2_ab

    # Distribution rates
    k12_ab <- (log(2) / tdist12_ab) * Pdist12_ab / (Pdist12_ab + Vcentral / Vperi)
    k21_ab <- (log(2) / tdist12_ab) * (Vcentral / Vperi) / (Pdist12_ab + Vcentral / Vperi)
    k12_pl <- (log(2) / tdist12_pl) * Pdist12_pl / (Pdist12_pl + Vcentral / Vperi)
    k21_pl <- (log(2) / tdist12_pl) * (Vcentral / Vperi) / (Pdist12_pl + Vcentral / Vperi)

    # HER2 amounts per cell
    HER2_0_tumor_pc  <- RPCher2 / 6.022e23 * 1e9
    HER2_0_normal_pc <- RPCnormal / 6.022e23 * 1e9
    HER2_endo_0_tumor_pc  <- kendoher2 / (krecher2 + kdegher2) * HER2_0_tumor_pc
    HER2_endo_0_normal_pc <- kendoher2 / (krecher2 + kdegher2) * HER2_0_normal_pc

    ksynth_tumor_pc  <- HER2_0_tumor_pc  * kendoher2 * kdegher2 / (kdegher2 + krecher2)
    ksynth_normal_pc <- HER2_0_normal_pc * kendoher2 * kdegher2 / (kdegher2 + krecher2)

    # Soluble HER2 shedding rate (inferred from Table S2d equation)
    # Steady-state balance: k_shed calibrated to sHER2 central baseline
    # 8 ng/mL of 1e5 g/mol -> 8e-11 mol/L -> 0.08 nM (unit chain includes *1e3 for mL->L)
    C_sher2_baseline_nM <- CsHER2_c * 1e-9 / MWsHER2 * 1e9 * 1000   # ng/mL -> M -> nM
    # Simplified k_shed calibration: at steady state, shedding * total normal HER2 = kelim_sher2 * sher2_amt
    # k_shed_HER2 = (C_sHER2 * Vcentral) * kelim_sher2 / (RPC_normal * (Ncell_central + eff transfer))
    # Use approximation: k_shed = kelim_sher2 * C_sher2_baseline * Vcentral / (HER2_0_normal_pc * NcellCentral)
    kshed_her2 <- kelim_sher2 * C_sher2_baseline_nM * Vcentral / (HER2_0_normal_pc * NcellCentral + 1e-30)

    # sHER2 distribution rates
    k12_sher2 <- (log(2) / tdist12_sher2) * Pdist12_sher2 / (Pdist12_sher2 + Vcentral / Vperi)
    k21_sher2 <- (log(2) / tdist12_sher2) * (Vcentral / Vperi) / (Pdist12_sher2 + Vcentral / Vperi)

    # Krogh cylinder tumor uptake (uses dynamic tumor volume)
    Ntot <- n1 + n2 + n3 + n4
    Ntot_safe <- Ntot + 1e-12
    Vtumor <- (Ntot_safe * Vcell) / (1 - eps_pl)
    Vext_tumor_ab <- eps_ab * Vtumor
    Vext_tumor_pl <- eps_pl * Vtumor
    Vint_tumor    <- (1 - eps_pl) * Vtumor
    Vtumor_mm3 <- Vtumor * 1e6

    Pab_h <- Pab / 24; Dab_h <- Dab / 24; Ppl_h <- Ppl / 24; Dpl_h <- Dpl / 24
    krogh_ab <- 2 * Pab_h * Rcap / (Rkrogh^2) + 6 * Dab_h * Rcap / (Vtumor_mm3^(2/3) * 1000)
    krogh_pl <- 2 * Ppl_h * Rcap / (Rkrogh^2) + 6 * Dpl_h * Rcap / (Vtumor_mm3^(2/3) * 1000)
    k13_ab <- krogh_ab * Vtumor / Vcentral
    k31_ab <- krogh_ab / eps_ab
    k13_pl <- krogh_pl * Vtumor / Vcentral
    k31_pl <- krogh_pl / eps_pl

    # sHER2 tumor distribution
    k13_sher2 <- (log(2) / tdist13_sher2) * Pdist12_sher2 / (Pdist13_sher2 + Vcentral / Vext_tumor_ab)
    k31_sher2 <- (log(2) / tdist13_sher2) * (Vcentral / Vext_tumor_ab) / (Pdist13_sher2 + Vcentral / Vext_tumor_ab)

    # -------------------------------------------------------------------
    # Reaction rates (systemic PK, tumor uptake, tumor cell interior)
    # -------------------------------------------------------------------
    # Systemic (central)
    v_dec_c   <- kdec * adc_central
    v_elim_adc_c <- kelim_ab * adc_central
    v_elim_ab_c  <- kelim_ab * ab_central
    v_elim_pl_c  <- kelim_pl * pl_central
    v_c2p_adc <- k12_ab * adc_central; v_p2c_adc <- k21_ab * adc_peripheral
    v_c2p_ab  <- k12_ab * ab_central;  v_p2c_ab  <- k21_ab * ab_peripheral
    v_c2p_pl  <- k12_pl * pl_central;  v_p2c_pl  <- k21_pl * pl_peripheral
    v_c2t_adc <- k13_ab * adc_central; v_t2c_adc <- k31_ab * adc_ext_tumor
    v_c2t_ab  <- k13_ab * ab_central;  v_t2c_ab  <- k31_ab * ab_ext_tumor
    v_c2t_pl  <- k13_pl * pl_central;  v_t2c_pl  <- k31_pl * pl_ext_tumor

    v_dec_p   <- kdec * adc_peripheral
    v_elim_adc_p <- kelim_ab * adc_peripheral
    v_elim_ab_p  <- kelim_ab * ab_peripheral
    v_elim_pl_p  <- kelim_pl * pl_peripheral

    v_dec_t <- kdec * adc_ext_tumor

    # Normal cell HER2 sinks (central compartment)
    v_synth_nc_c <- ksynth_normal_pc * NcellCentral
    v_synth_nc_p <- ksynth_normal_pc * NcellPeri

    # v79-v83, v88-v92 in paper: normal-cell HER2 binding + endocytosis + recycling + degradation
    v79f <- konab / Vcentral * her2_nc_c * adc_central; v79r <- koffab * her2adc_nc_c
    v80f <- konab / Vcentral * her2_nc_c * ab_central;  v80r <- koffab * her2ab_nc_c
    v81f <- kendoher2 * her2_nc_c;         v81r <- krecher2 * her2endo_nc_c
    v82f <- kendoher2 * her2adc_nc_c;      v82r <- krecher2 * her2adcendo_nc_c
    v83f <- kendoher2 * her2ab_nc_c;       v83r <- krecher2 * her2abendo_nc_c
    v84  <- kdegher2 * her2endo_nc_c
    v85  <- kdegher2 * her2adcendo_nc_c
    v86  <- kdegher2 * her2abendo_nc_c

    v88f <- konab / Vperi * her2_nc_p * adc_peripheral; v88r <- koffab * her2adc_nc_p
    v89f <- konab / Vperi * her2_nc_p * ab_peripheral;  v89r <- koffab * her2ab_nc_p
    v90f <- kendoher2 * her2_nc_p;         v90r <- krecher2 * her2endo_nc_p
    v91f <- kendoher2 * her2adc_nc_p;      v91r <- krecher2 * her2adcendo_nc_p
    v92f <- kendoher2 * her2ab_nc_p;       v92r <- krecher2 * her2abendo_nc_p
    v93  <- kdegher2 * her2endo_nc_p
    v94  <- kdegher2 * her2adcendo_nc_p
    v95  <- kdegher2 * her2abendo_nc_p

    # sHER2 dynamics (shedding + binding + clearance + distribution)
    v96  <- kshed_her2 * her2_nc_c                     # shed from HER2 (central)
    v97  <- kelim_sher2 * sher2_c                       # sHER2 clearance central
    v98  <- kshed_her2 * chi_shed_ab * her2adc_nc_c    # shed from HER2:ADC
    v99  <- kshed_her2 * chi_shed_ab * her2ab_nc_c     # shed from HER2:Ab
    v100 <- kelim_sher2_ab * sher2adc_c                # sHER2:ADC clearance
    v101 <- kelim_sher2_ab * sher2ab_c                 # sHER2:Ab clearance
    v102f <- konab / Vcentral * sher2_c * adc_central; v102r <- koffab * sher2adc_c
    v103f <- konab / Vcentral * sher2_c * ab_central;  v103r <- koffab * sher2ab_c

    v104 <- kshed_her2 * her2_nc_p
    v105 <- kelim_sher2 * sher2_p
    v106 <- kshed_her2 * chi_shed_ab * her2adc_nc_p
    v107 <- kshed_her2 * chi_shed_ab * her2ab_nc_p
    v108 <- kelim_sher2_ab * sher2adc_p
    v109 <- kelim_sher2_ab * sher2ab_p
    v110f <- konab / Vperi * sher2_p * adc_peripheral; v110r <- koffab * sher2adc_p
    v111f <- konab / Vperi * sher2_p * ab_peripheral;  v111r <- koffab * sher2ab_p

    v112f <- k12_sher2 * sher2_c; v112r <- k21_sher2 * sher2_p
    v113f <- k13_sher2 * sher2_c; v113r <- k31_sher2 * sher2_tumor
    v114f <- k12_sher2 * sher2adc_c; v114r <- k21_sher2 * sher2adc_p
    v115f <- k13_sher2 * sher2adc_c; v115r <- k31_sher2 * sher2adc_tumor
    v116f <- k12_sher2 * sher2ab_c; v116r <- k21_sher2 * sher2ab_p
    v117f <- k13_sher2 * sher2ab_c; v117r <- k31_sher2 * sher2ab_tumor

    # Tumor cell HER2 shedding into sHER2_tumor
    v118 <- kshed_her2 * her2_surf_tumor
    v119 <- kshed_her2 * chi_shed_ab * her2_ab_surf_tumor
    v120 <- kshed_her2 * chi_shed_ab * her2_adc_surf_tumor

    v121f <- konab / Vext_tumor_ab * sher2_tumor * adc_ext_tumor; v121r <- koffab * sher2adc_tumor
    v122f <- konab / Vext_tumor_ab * sher2_tumor * ab_ext_tumor;  v122r <- koffab * sher2ab_tumor

    v123 <- kdec * sher2adc_c
    v124 <- kdec * sher2adc_p
    v125 <- kdec * sher2adc_tumor
    v126 <- kdec * her2adc_nc_c
    v127 <- kdec * her2adc_nc_p

    # Tumor cell interior (aggregate; same as mouse model)
    v_synth_tumor <- ksynth_tumor_pc * Ntot_safe
    v30f <- konab / Vext_tumor_ab * her2_surf_tumor * adc_ext_tumor
    v30r <- koffab * her2_adc_surf_tumor
    v31f <- konab / Vext_tumor_ab * her2_surf_tumor * ab_ext_tumor
    v31r <- koffab * her2_ab_surf_tumor
    v32 <- kdec * her2_adc_surf_tumor
    v33f <- kendoher2 * her2_surf_tumor;     v33r <- krecher2 * her2_endo_tumor
    v34f <- kendoher2 * her2_adc_surf_tumor; v34r <- krecher2 * her2_adc_endo_tumor
    v35f <- kendoher2 * her2_ab_surf_tumor;  v35r <- krecher2 * her2_ab_endo_tumor
    v36 <- kdegher2 * her2_endo_tumor
    v37 <- kdegher2 * her2_adc_endo_tumor
    v38 <- kdegher2 * her2_ab_endo_tumor
    v39 <- kcleave * her2_adc_endo_tumor
    v40 <- kinpl * pl_endo_tumor
    v41f <- konpl / (Pcpl * Vint_tumor) * t_cyto_tumor * pl_cyto_tumor
    v41r <- koffpl * tpl_cyto_tumor
    v_pl_in_tumor  <- kinpl  * pl_ext_tumor
    v_pl_out_tumor <- koutpl * pl_cyto_tumor

    # -------------------------------------------------------------------
    # Tumor growth kinetics (Simeoni transit chain)
    # -------------------------------------------------------------------
    Cpl_intra <- (pl_cyto_tumor + tpl_cyto_tumor) / Vint_tumor
    kkill <- kkillmax * Cpl_intra^nHill / (kc50^nHill + Cpl_intra^nHill)

    Nmax <- (1 - eps_pl) * Vtumor_max / Vcell
    growth_num <- (log(2) / tdouble) * (1 - Ntot_safe / Nmax)
    growth_denom_arg <- log(2) * Ntot_safe * Vcell / (tdouble * (1 - eps_pl) * klin_h * 1e-6)
    growth_denom_pos <- 1 + max(growth_denom_arg, 0)^psi
    k_growth <- growth_num / (growth_denom_pos^(1/psi))

    death_frac <- (n4 / Ntot_safe) / tau_kill

    # -------------------------------------------------------------------
    # ODEs
    # -------------------------------------------------------------------
    # Systemic (central) - includes normal HER2 + sHER2 sink contributions
    d/dt(adc_central) <- -v_dec_c - v_elim_adc_c - v_c2p_adc + v_p2c_adc - v_c2t_adc + v_t2c_adc -
                          v79f + v79r - v102f + v102r
    d/dt(ab_central)  <-  v_dec_c - v_elim_ab_c  - v_c2p_ab  + v_p2c_ab  - v_c2t_ab  + v_t2c_ab -
                          v80f + v80r - v103f + v103r
    d/dt(pl_central)  <-  v_dec_c * DAR + v_elim_adc_c * DAR - v_elim_pl_c -
                          v_c2p_pl + v_p2c_pl - v_c2t_pl + v_t2c_pl + v126 * DAR + v123 * DAR

    d/dt(adc_peripheral) <- v_c2p_adc - v_p2c_adc - v_dec_p - v_elim_adc_p - v88f + v88r - v110f + v110r
    d/dt(ab_peripheral)  <- v_c2p_ab  - v_p2c_ab  + v_dec_p - v_elim_ab_p  - v89f + v89r - v111f + v111r
    d/dt(pl_peripheral)  <- v_c2p_pl  - v_p2c_pl  + v_dec_p * DAR + v_elim_adc_p * DAR - v_elim_pl_p +
                            v127 * DAR + v124 * DAR

    # Tumor extracellular
    d/dt(adc_ext_tumor) <- v_c2t_adc - v_t2c_adc - v_dec_t - v30f + v30r - v121f + v121r
    d/dt(ab_ext_tumor)  <- v_c2t_ab  - v_t2c_ab  + v_dec_t - v31f + v31r - v122f + v122r
    d/dt(pl_ext_tumor)  <- v_c2t_pl  - v_t2c_pl  + v_dec_t * DAR + v32 * DAR + v125 * DAR +
                            v_pl_out_tumor - v_pl_in_tumor

    # Normal HER2 in central + endo + degrades
    d/dt(her2_nc_c)         <- v_synth_nc_c - v79f + v79r - v80f + v80r - v81f + v81r - v96
    d/dt(her2adc_nc_c)      <- v79f - v79r - v82f + v82r - v98 - v126
    d/dt(her2ab_nc_c)       <- v80f - v80r - v83f + v83r - v99 + v126
    d/dt(her2endo_nc_c)     <- v81f - v81r - v84
    d/dt(her2adcendo_nc_c)  <- v82f - v82r - v85
    d/dt(her2abendo_nc_c)   <- v83f - v83r - v86

    # Normal HER2 in peripheral
    d/dt(her2_nc_p)         <- v_synth_nc_p - v88f + v88r - v89f + v89r - v90f + v90r - v104
    d/dt(her2adc_nc_p)      <- v88f - v88r - v91f + v91r - v106 - v127
    d/dt(her2ab_nc_p)       <- v89f - v89r - v92f + v92r - v107 + v127
    d/dt(her2endo_nc_p)     <- v90f - v90r - v93
    d/dt(her2adcendo_nc_p)  <- v91f - v91r - v94
    d/dt(her2abendo_nc_p)   <- v92f - v92r - v95

    # Soluble HER2 in central
    d/dt(sher2_c)      <-  v96 - v97 - v102f + v102r - v103f + v103r - v112f + v112r - v113f + v113r
    d/dt(sher2adc_c)   <- -v100 + v102f - v102r - v114f + v114r - v123
    d/dt(sher2ab_c)    <- -v101 + v103f - v103r - v116f + v116r

    # Soluble HER2 in peripheral
    d/dt(sher2_p)      <-  v104 - v105 - v110f + v110r - v111f + v111r + v112f - v112r
    d/dt(sher2adc_p)   <- -v108 + v110f - v110r + v114f - v114r - v124
    d/dt(sher2ab_p)    <- -v109 + v111f - v111r + v116f - v116r

    # Soluble HER2 in tumor (shed from tumor cell HER2 + distribution from central)
    d/dt(sher2_tumor)     <- v118 - v121f + v121r + v113f - v113r
    d/dt(sher2adc_tumor)  <- v121f - v121r + v115f - v115r - v125
    d/dt(sher2ab_tumor)   <- v122f - v122r + v117f - v117r + v125

    # Tumor cell interior (aggregate)
    d/dt(her2_surf_tumor)     <- v_synth_tumor - v33f + v33r - v30f + v30r - v31f + v31r - v118 - her2_surf_tumor * death_frac
    d/dt(her2_adc_surf_tumor) <- v30f - v30r - v34f + v34r - v32 - v120 - her2_adc_surf_tumor * death_frac
    d/dt(her2_ab_surf_tumor)  <- v31f - v31r + v32 - v35f + v35r - v119 - her2_ab_surf_tumor * death_frac
    d/dt(her2_endo_tumor)     <- v33f - v33r - v36 - her2_endo_tumor * death_frac
    d/dt(her2_adc_endo_tumor) <- v34f - v34r - v37 - v39 - her2_adc_endo_tumor * death_frac
    d/dt(her2_ab_endo_tumor)  <- v35f - v35r - v38 + v39 - her2_ab_endo_tumor * death_frac
    d/dt(pl_endo_tumor)       <- (v37 + v39) * DAR - v40 - pl_endo_tumor * death_frac
    d/dt(pl_cyto_tumor)       <- v40 - v41f + v41r + v_pl_in_tumor - v_pl_out_tumor - pl_cyto_tumor * death_frac
    v_target_synth <- k_growth * n1 * Vcell * Tpercell
    d/dt(t_cyto_tumor)        <- v_target_synth - v41f + v41r - t_cyto_tumor * death_frac
    d/dt(tpl_cyto_tumor)      <- v41f - v41r - tpl_cyto_tumor * death_frac

    # Tumor cell staging (Simeoni)
    d/dt(n1) <- k_growth * n1 - kkill * n1
    d/dt(n2) <- kkill * n1 - n2 / tau_kill
    d/dt(n3) <- n2 / tau_kill - n3 / tau_kill
    d/dt(n4) <- n3 / tau_kill - n4 / tau_kill

    # -------------------------------------------------------------------
    # Initial conditions
    # -------------------------------------------------------------------
    N_i <- (1 - eps_pl) * Vtumor_i / Vcell
    n1(0) <- N_i
    her2_surf_tumor(0) <- HER2_0_tumor_pc * N_i
    her2_endo_tumor(0) <- HER2_endo_0_tumor_pc * N_i
    t_cyto_tumor(0)    <- Tpercell * Vcell * N_i

    # Normal cell HER2 initial (steady state, cells * per-cell HER2)
    her2_nc_c(0)         <- HER2_0_normal_pc * NcellCentral
    her2_nc_p(0)         <- HER2_0_normal_pc * NcellPeri
    her2endo_nc_c(0)     <- HER2_endo_0_normal_pc * NcellCentral
    her2endo_nc_p(0)     <- HER2_endo_0_normal_pc * NcellPeri

    # Soluble HER2 initial (baseline)
    sher2_c(0)   <- C_sher2_baseline_nM * Vcentral
    sher2_p(0)   <- C_sher2_baseline_nM * Vperi * Pdist12_sher2

    # -------------------------------------------------------------------
    # Observations (nM plasma / mm^3 tumor / nM intracellular)
    # -------------------------------------------------------------------
    Cc         <- adc_central / Vcentral
    Cab        <- ab_central  / Vcentral
    Cpl_plasma <- pl_central  / Vcentral
    Csher2_plasma <- sher2_c / Vcentral
    Vtumor_out <- Vtumor_mm3
    Cpl_intra_out <- Cpl_intra

    Cc ~ prop(propSd)
  })
}
