Yang_2024_axatilimab <- function() {
  description <- "Semimechanistic population PK/PD model for axatilimab (anti-CSF-1R IgG4 monoclonal antibody) in healthy adults, patients with advanced solid tumors, and patients with chronic graft-versus-host disease (Yang 2024). Two-compartment IV PK with parallel linear clearance and CSF-1R-mediated saturable elimination via competitive Hill binding with circulating CSF-1; CSF-1, NCMC, AST, and CPK pharmacodynamics integrated as turnover indirect-response biomarkers."
  reference <- "Yang Y, Sokolov V, Volkova A, et al. Semimechanistic Population PK/PD Modeling of Axatilimab in Healthy Participants and Patients With Solid Tumors or Chronic Graft-Versus-Host Disease. Clin Pharmacol Ther. 2025;117(3):704-714. doi:10.1002/cpt.3503"
  vignette <- "Yang_2024_axatilimab"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on volume of distribution (Vd) with reference 73.6 kg (overall pooled-cohort median, Yang 2024 Table S3); exponent 0.7. Time-fixed at baseline; the paper states continuous covariates were adjusted to the population median.",
      source_name        = "WT"
    ),
    CSF1 = list(
      description        = "Baseline plasma CSF-1 (M-CSF) concentration; time-fixed at the pre-dose value",
      units              = "pg/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on linear clearance (CL, exponent 0.912) and on the typical baseline CSF-1 model parameter BL_CSF1 (exponent 0.656); reference value 549 pg/mL (overall pooled-cohort median, Yang 2024 Table S3, reported there as 549 ng/L which equals 549 pg/mL). Distinct from the model state csf1 (time-course of CSF-1 in nM driven by axatilimab/CSF-1R competitive binding).",
      source_name        = "BLCSF1"
    ),
    CPK = list(
      description        = "Baseline serum creatine phosphokinase (CPK / creatine kinase)",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on baseline NCMC parameter BL_NCMC; exponent 0.376; reference 63 U/L (overall pooled-cohort median, Yang 2024 Table S3). Distinct from the model state cpk (time-course of CPK in U/L driven by NCMC-dependent indirect response).",
      source_name        = "BLCPK"
    ),
    DIS_CANCER = list(
      description        = "Advanced-solid-tumor cohort indicator (1 = solid-tumor patient)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (cGVHD or healthy volunteer)",
      notes              = "One of two orthogonal indicators (with DIS_HV) decomposing the three-level participant population categorical (cGVHD reference, advanced solid tumor, healthy volunteer). Effect on baseline NCMC: BL_NCMC * exp(1.22 * DIS_CANCER + 0.618 * DIS_HV); reference category cGVHD when both indicators are 0. Yang 2024 Table 1 row 'Population with cancer on BLNCMC'.",
      source_name        = "POPULATION (Population type = 'Patients with cancer')"
    ),
    DIS_HV = list(
      description        = "Healthy-volunteer cohort indicator (1 = healthy volunteer, 0 = patient)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (cGVHD or advanced solid tumor)",
      notes              = "Paired with DIS_CANCER. Effect on baseline NCMC: BL_NCMC * exp(1.22 * DIS_CANCER + 0.618 * DIS_HV); reference category cGVHD. Yang 2024 Table 1 row 'Healthy population on BLNCMC'.",
      source_name        = "POPULATION (Population type = 'Healthy participants')"
    ),
    ADA_POS = list(
      description        = "Time-varying antidrug-antibody (ADA) positivity (1 = positive at time t, 0 = negative)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative; typical patient)",
      notes              = "Time-varying covariate on linear clearance: CL * (1 + 0.489 * ADA_POS(t)). When ADA_POS(t) = 1 the linear clearance increases by ~50.6% relative to the ADA-negative reference (Yang 2024 final covariate-effect summary). Source MLXTRAN regressor name 'ADACN'.",
      source_name        = "ADA"
    )
  )

  population <- list(
    n_subjects     = 325L,
    n_studies      = 4L,
    age_range      = "7-81 years (7 children aged 7-15 with cGVHD)",
    age_median     = "55 years",
    weight_range   = "18.1-151 kg",
    weight_median  = "73.6 kg",
    sex_female_pct = 39.4,
    race_ethnicity = c(White = 82.5, Asian = 5.5, Black = 4.3, Other = 1.8, Unknown = 5.8),
    disease_state  = "Pooled cohort: 14 healthy adults (SNDX-6352-0001), 33 adults with advanced or metastatic solid tumors (SNDX-6352-0502), 271 adults with chronic graft-versus-host disease (cGVHD) and 7 children with cGVHD aged 7-15 (SNDX-6352-0503 phase 1/2 dose-finding; AGAVE-201 phase 2).",
    dose_range     = "0.15 to 6 mg/kg IV; single dose, every 2 weeks (Q2W), or every 4 weeks (Q4W); 0.3 mg/kg Q2W is the approved cGVHD regimen.",
    regions        = "Multi-regional (multi-center US-led trials; phase 1 healthy-volunteer study had EudraCT registration in Europe).",
    n_observations = list(PK = 5048L, CSF1 = 1659L, NCMC = 1683L, AST = 4966L, CPK = 4571L),
    n_observations_below_LOQ_pct = 35.5,
    notes          = "Baseline demographics from Yang 2024 Tables S3 and S4. ADA prevalence (ever-positive) 40% overall. Plasma axatilimab measured by ELISA with LLOQ 157 ng/mL (phase 1 studies) or 150 ng/mL (phase 1/2 + phase 2 cGVHD studies); CSF-1 measured by R&D Systems Quantikine human M-CSF ELISA; NCMC measured by qualified flow-cytometry assay. Below-quantification PK samples (35.5% of all PK observations, almost all post-first-dose tail) handled by the M4 method (Beal 2001); the structural / statistical model carries all BLOQ data."
  )

  ini({
    # === Structural PK parameters (all on log scale) ===
    # Yang 2024 Table 1 - reference covariate values: WT 73.6 kg, baseline CSF-1 549 pg/mL,
    # baseline CPK 63 U/L, ADA-negative, cGVHD population (DIS_CANCER = DIS_HV = 0).
    lvc       <- log(3.48);   label("Central volume of distribution Vd (L) for a typical 73.6 kg cGVHD patient")                                              # Yang 2024 Table 1 row 'Volume of distribution (Vd)' = 3.48 L
    lcl       <- log(0.007);  label("Linear clearance CL (L/h) for a typical ADA-negative cGVHD patient at median CSF-1")                                       # Yang 2024 Table 1 row 'Clearance (CL)' = 0.007 L/h
    lq        <- log(0.015);  label("Inter-compartmental clearance Q (L/h)")                                                                                    # Yang 2024 Table 1 row 'Intercompartmental CL (Q)' = 0.015 L/h
    lvp       <- log(2.64);   label("Peripheral volume of distribution Vp (L)")                                                                                 # Yang 2024 Table 1 row 'Volume of distribution in the peripheral compartment (Vp)' = 2.64 L

    # === Target-mediated elimination (CSF-1R-driven, Hill cooperativity) ===
    lvmax     <- log(0.37);   label("Maximum rate of CSF-1R-mediated elimination of axatilimab/CSF-1 complexes Vmax (nM/h)")                                    # Yang 2024 Table 1 row 'Elimination rates of the CSF-1 and axatilimab complexes with CSF-1R (Vmax)' = 0.37 nM/h
    lkd_pk    <- log(1.11);   label("Dissociation constant of axatilimab/CSF-1R complex Kd_PK (nM)")                                                            # Yang 2024 Table 1 row 'Dissociation constant of axatilimab/CSF-1R complex (KdPK)' = 1.11 nM
    lnh       <- log(2.5);    label("Hill coefficient for axatilimab cooperativity binding to CSF-1R (Nh, unitless)")                                           # Yang 2024 Table 1 row 'Hill coefficient (Nh)' = 2.5

    # === CSF-1 turnover ===
    lbl_csf1  <- log(0.01);   label("Baseline CSF-1 plasma concentration BL_CSF1 (nM) for typical median-CSF-1 patient")                                        # Yang 2024 Table 1 row 'Baseline CSF-1 concentration (BLCSF1)' = 0.01 nM
    lkdeg_csf1 <- log(0.002); label("CSF-1R-independent CSF-1 elimination rate kdeg_CSF1 (1/h)")                                                                 # Yang 2024 Table 1 row 'CSF-1R-independent CSF-1 elimination rate (kdegCSF1)' = 0.002 1/h

    # === NCMC (nonclassical CD14+/CD16++ monocytic cell) turnover ===
    lbl_ncmc  <- log(16);     label("Baseline NCMC concentration BL_NCMC (cells/uL) for typical cGVHD patient at median CPK")                                   # Yang 2024 Table 1 row 'Baseline NCMC concentration (BLNCMC)' = 16 cells/uL
    lkdeg_ncmc <- log(0.021); label("NCMC elimination rate kdeg_NCMC (1/h)")                                                                                    # Yang 2024 Table 1 row 'NCMC elimination rate (kdegNCMC)' = 0.021 1/h

    # === ADA effect coefficient (kada) ===
    # Source: Yang 2024 Table 1 row 'ADA effect coefficient (kada)' = 0.489.
    # NOTE on internal source contradiction: the page-708 narrative states 'kada = 2.1
    # ... ~3.1-fold increase in CL', but Table 1, the MLXTRAN supplement (Supplement 2),
    # and the formal covariate-effect summary ('Positive ADA status was associated with
    # a 50.6% (95% CI, 47.0-54.2) increase in CL') all give kada = 0.489 with the
    # CL * (1 + kada * ADA(t)) parameterization. The Errata section of the validation
    # vignette captures this discrepancy; the canonical value used here is 0.489.
    e_ada_cl  <- 0.489;       label("ADA effect coefficient on linear CL (kada, applied as CL * (1 + e_ada_cl * ADA_POS(t)))")                                   # Yang 2024 Table 1 row 'ADA effect coefficient (kada)' = 0.489

    # === AST (aspartate aminotransferase) - NCMC-driven indirect response ===
    lbl_ast            <- log(35.5);  label("Baseline AST concentration BL_AST (U/L)")                                                                          # Yang 2024 Table 1 row 'Baseline AST concentration (BLAST)' = 35.5 U/L
    lvmax_ast_ncmc     <- log(0.011); label("Max rate of NCMC-dependent AST elimination Vmax_AST_NCMC (1/h)")                                                   # Yang 2024 Table 1 row 'Maximum rate of NCMC-dependent elimination of AST (VmaxAST_NCMC)' = 0.011 1/h
    lkdeg_ast          <- log(0.002); label("NCMC-independent AST elimination rate kdeg_AST (1/h)")                                                             # Yang 2024 Table 1 row 'AST elimination rate (kdegAST)' = 0.002 1/h
    lec50_ast_ncmc     <- log(11.8);  label("NCMC concentration giving 50% of max NCMC-dependent AST elimination EC50_AST_NCMC (cells/uL)")                     # Yang 2024 Table 1 row 'NCMC concentration resulting in 50% of maximum rate of NCMC-dependent elimination of AST (EC50AST_NCMC)' = 11.8 cells/uL

    # === CPK (creatine phosphokinase) - NCMC-driven indirect response ===
    lbl_cpk            <- log(101);   label("Baseline CPK concentration BL_CPK (U/L) for typical median-CPK patient")                                           # Yang 2024 Table 1 row 'Baseline CPK concentration (BLCPK)' = 101 U/L
    lvmax_cpk_ncmc     <- log(0.02);  label("Max rate of NCMC-dependent CPK elimination Vmax_CPK_NCMC (1/h)")                                                   # Yang 2024 Table 1 row 'Maximum rate of NCMC-dependent elimination of CPK (VmaxCPK_NCMC)' = 0.02 1/h
    lkdeg_cpk          <- log(0.001); label("NCMC-independent CPK elimination rate kdeg_CPK (1/h)")                                                             # Yang 2024 Table 1 row 'CPK elimination rate (kdegCPK)' = 0.001 1/h
    lec50_cpk_ncmc     <- log(19.2);  label("NCMC concentration giving 50% of max NCMC-dependent CPK elimination EC50_CPK_NCMC (cells/uL)")                     # Yang 2024 Table 1 row 'NCMC concentration resulting in 50% of maximum rate of NCMC-dependent elimination of CPK (EC50CPK_NCMC)' = 19.2 cells/uL

    # === Continuous-covariate effects (power form on the scaled covariate ratio) ===
    e_wt_vc           <- 0.7;    label("Power exponent of (WT/73.6) on Vd (unitless)")                                                                          # Yang 2024 Table 1 row 'Bodyweight on Vd' = 0.7
    e_csf1_cl         <- 0.912;  label("Power exponent of (CSF1/549 pg/mL) on linear CL (unitless)")                                                            # Yang 2024 Table 1 row 'CSF-1 on CL' = 0.912
    e_csf1_blcsf1     <- 0.656;  label("Power exponent of (CSF1/549 pg/mL) on BL_CSF1 (unitless)")                                                              # Yang 2024 Table 1 row 'CSF-1 on BLCSF1' = 0.656
    e_cancer_blncmc   <- 1.22;   label("Exponential coefficient on DIS_CANCER for BL_NCMC (unitless; cGVHD reference)")                                         # Yang 2024 Table 1 row 'Population with cancer on BLNCMC' = 1.22
    e_hv_blncmc       <- 0.618;  label("Exponential coefficient on DIS_HV for BL_NCMC (unitless; cGVHD reference)")                                             # Yang 2024 Table 1 row 'Healthy population on BLNCMC' = 0.618
    e_cpk_blncmc      <- 0.376;  label("Power exponent of (CPK/63 U/L) on BL_NCMC (unitless)")                                                                  # Yang 2024 Table 1 row 'BLCPK on BLNCMC' = 0.376

    # === Inter-individual variability (BSV; Monolix omega = SD on log-scale; ===
    # === nlmixr2 expects variance, so use omega^2.) ===
    # Yang 2024 Table 1 reports omega; squared values used here.
    etalvc      ~ 0.055225  # Yang 2024 Table 1 omega(Vd) = 0.235 -> 0.235^2 = 0.055225
    etalcl      ~ 1.1881    # Yang 2024 Table 1 omega(CL) = 1.09 -> 1.09^2 = 1.1881
    etalvmax    ~ 0.083521  # Yang 2024 Table 1 omega(Vmax) = 0.289 -> 0.289^2 = 0.083521 (eta-shrinkage 51.2%)
    etalbl_csf1 ~ 0.034969  # Yang 2024 Table 1 omega(BLCSF1) = 0.187 -> 0.187^2 = 0.034969
    etalbl_ncmc ~ 0.753424  # Yang 2024 Table 1 omega(BLNCMC) = 0.868 -> 0.868^2 = 0.753424
    etalbl_ast  ~ 0.215296  # Yang 2024 Table 1 omega(BLAST)  = 0.464 -> 0.464^2 = 0.215296
    etalbl_cpk  ~ 0.570025  # Yang 2024 Table 1 omega(BLCPK)  = 0.755 -> 0.755^2 = 0.570025

    # === Residual error (Yang 2024 Table 1 'Residual error model' rows) ===
    CcpropSd        <- 0.375; label("Proportional residual error on axatilimab Cc (fraction)")                                                                   # Yang 2024 Table 1 'bPK' = 0.375
    propSd_csf1obs   <- 0.321; label("Proportional residual error on plasma CSF-1 observation (fraction)")                                                      # Yang 2024 Table 1 'bCSF1' = 0.321
    addSd_ncmcobs    <- 0.977; label("Additive residual error on NCMC observation (cells/uL)")                                                                  # Yang 2024 Table 1 'aNCMC' = 0.977
    propSd_ncmcobs   <- 0.676; label("Proportional residual error on NCMC observation (fraction)")                                                              # Yang 2024 Table 1 'bNCMC' = 0.676
    propSd_astobs    <- 0.324; label("Proportional residual error on AST observation (fraction)")                                                               # Yang 2024 Table 1 'bAST' = 0.324
    propSd_cpkobs    <- 0.462; label("Proportional residual error on CPK observation (fraction)")                                                               # Yang 2024 Table 1 'bCPK' = 0.462
  })

  model({
    # === Fixed mechanistic constants ===
    # Yang 2024 fixed Kd_CSF1 to 0.048 nM (48 pM) per Roussel et al. (1988) Cell, since
    # simultaneous estimation with Kd_PK was not practically identifiable in the absence
    # of CSF-1R measurements. (Yang 2024 Methods, paragraph after Eq. 7.)
    Kd_CSF1 <- 0.048

    # Molecular weights (kDa). Used to convert PK from mg/L to nM and CSF-1 from nM to ng/mL.
    MW_PK   <- 150
    MW_CSF1 <- 60

    # === Individual PK parameters with Yang 2024 covariate adjustments ===
    vc          <- exp(lvc       + etalvc)      * (WT / 73.6)^e_wt_vc
    cl          <- exp(lcl       + etalcl)      * (CSF1 / 549)^e_csf1_cl * (1 + e_ada_cl * ADA_POS)
    q           <- exp(lq)
    vp          <- exp(lvp)
    Vmax        <- exp(lvmax     + etalvmax)
    Kd_PK       <- exp(lkd_pk)
    Nh          <- exp(lnh)
    BL_CSF1     <- exp(lbl_csf1  + etalbl_csf1) * (CSF1 / 549)^e_csf1_blcsf1
    BL_NCMC     <- exp(lbl_ncmc  + etalbl_ncmc) *
                     exp(e_cancer_blncmc * DIS_CANCER + e_hv_blncmc * DIS_HV) *
                     (CPK / 63)^e_cpk_blncmc
    kdeg_csf1   <- exp(lkdeg_csf1)
    kdeg_ncmc   <- exp(lkdeg_ncmc)
    BL_AST            <- exp(lbl_ast  + etalbl_ast)
    Vmax_AST_NCMC     <- exp(lvmax_ast_ncmc)
    kdeg_ast          <- exp(lkdeg_ast)
    EC50_AST_NCMC     <- exp(lec50_ast_ncmc)
    BL_CPK            <- exp(lbl_cpk  + etalbl_cpk)
    Vmax_CPK_NCMC     <- exp(lvmax_cpk_ncmc)
    kdeg_cpk          <- exp(lkdeg_cpk)
    EC50_CPK_NCMC     <- exp(lec50_cpk_ncmc)

    # === Plasma concentrations (Yang 2024 Eq. 4 and the analogous Cp form) ===
    # central in mg, vc in L -> Cc in mg/L. Cp in mg/L by analogous division.
    Cc      <- central / vc
    Cp      <- peripheral1 / vp
    # Convert axatilimab plasma concentration (mg/L = ug/mL) to nM for the
    # Hill saturable-clearance term: 1 mg/L / (MW kDa / 1000) gives nmol/L.
    Cc_nM   <- Cc * 1000 / MW_PK

    # === Saturable-clearance fractional-occupancy expressions (Yang 2024 Eq. 3 and 7) ===
    occ_PK   <- (Cc_nM / Kd_PK)^Nh
    occ_CSF1 <- csf1   / Kd_CSF1
    denom    <- 1 + occ_PK + occ_CSF1
    C1 <- occ_PK   / denom    # Eq. 3
    C2 <- occ_CSF1 / denom    # Eq. 7

    # === Steady-state synthesis-rate expressions (Yang 2024 Eqs. 6, 9 and supplement Eq. 19) ===
    BL_C2     <- (BL_CSF1 / Kd_CSF1) / (1 + BL_CSF1 / Kd_CSF1)                  # baseline value of C2 with Cc_nM = 0
    ksyn_csf1 <- kdeg_csf1 * BL_CSF1 + Vmax * BL_C2                              # Eq. 6
    ksyn_ncmc <- kdeg_ncmc * BL_NCMC / BL_C2                                     # Eq. 9
    ksyn_ast  <- BL_AST * (Vmax_AST_NCMC * BL_NCMC / (EC50_AST_NCMC + BL_NCMC) + kdeg_ast)  # supplement Eq. 19
    ksyn_cpk  <- BL_CPK * (Vmax_CPK_NCMC * BL_NCMC / (EC50_CPK_NCMC + BL_NCMC) + kdeg_cpk)  # supplement Eq. 19

    # === ODE system (Yang 2024 Eqs. 1, 2, 5, 8, 10, 11) ===
    # Axatilimab amounts (central, peripheral1) in mg; CSF-1 in nM; NCMC in cells/uL;
    # AST and CPK in U/L. Saturable-elimination amount-term: Vmax [nM/h] x C1 x vc [L]
    # converted to mg/h via x MW_PK / 1000 (1 nmol of MW_PK kDa weighs MW_PK/1000 mg).
    d/dt(central)     <- -cl * Cc * (1 + e_ada_cl * ADA_POS) - q * (Cc - Cp) -
                          Vmax * C1 * vc * MW_PK / 1000
    d/dt(peripheral1) <-  q * (Cc - Cp)
    d/dt(csf1)        <-  ksyn_csf1 - kdeg_csf1 * csf1 - Vmax * C2
    d/dt(ncmc)        <-  ksyn_ncmc * C2 - kdeg_ncmc * ncmc
    d/dt(ast)         <-  ksyn_ast  - (Vmax_AST_NCMC * ncmc / (EC50_AST_NCMC + ncmc) + kdeg_ast) * ast
    d/dt(cpk)         <-  ksyn_cpk  - (Vmax_CPK_NCMC * ncmc / (EC50_CPK_NCMC + ncmc) + kdeg_cpk) * cpk

    # NOTE: the linear-CL ODE term carries the (1 + e_ada_cl * ADA_POS) factor explicitly
    # rather than re-using the cl variable computed above. This keeps the time-varying
    # ADA effect dynamic at every solver step (cl is computed once per dose interval in
    # rxode2; ADA_POS may change with each event row).

    # === Initial conditions: biomarker states start at their individual baselines ===
    csf1(0) <- BL_CSF1
    ncmc(0) <- BL_NCMC
    ast(0)  <- BL_AST
    cpk(0)  <- BL_CPK

    # === Observations (one per modeled output) ===
    # Cc - axatilimab in mg/L (= ug/mL).
    # csf1obs - plasma CSF-1 in ng/mL: 1 nM x MW_CSF1 (kDa) gives ng/mL.
    csf1obs <- csf1 * MW_CSF1
    # ncmcobs - NCMC in cells/uL (state already in observation units).
    ncmcobs <- ncmc
    # astobs / cpkobs - hepatic / muscle enzymes in U/L (state already in observation units).
    astobs  <- ast
    cpkobs  <- cpk

    Cc       ~ prop(CcpropSd)
    csf1obs  ~ prop(propSd_csf1obs)
    ncmcobs  ~ add(addSd_ncmcobs) + prop(propSd_ncmcobs)
    astobs   ~ prop(propSd_astobs)
    cpkobs   ~ prop(propSd_cpkobs)
  })
}
