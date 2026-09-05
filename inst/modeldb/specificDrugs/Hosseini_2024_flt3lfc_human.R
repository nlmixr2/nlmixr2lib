Hosseini_2024_flt3lfc_human <- function() {
  description <- "QSP. Human translation of the FLT3L-Fc (RO7497987) minimal PBPK/PD model with expansion-enhanced target-mediated drug disposition, used to select the 700 ug first-in-human dose. Couples (a) a Cao/Jusko two-tissue-group mPBPK model of FLT3L-Fc in which drug binds one or two monomeric FLT3 receptors in plasma and the double-bound homodimer expands the total FLT3 receptor pool, (b) an empirical one-compartment SC model with Michaelis-Menten elimination for the comparator recombinant FLT3 ligand CDX-301, and (c) a shared transit-delayed indirect-response model for conventional dendritic cell (cDC1 and cDC2) expansion in peripheral blood. PK parameters translated from cynomolgus monkey; FLT3 target parameters refined against GS-3583 clinical PK (human-derived scenario); PD calibrated to CDX-301 and recombinant FLT3L data in healthy volunteers (Hosseini 2024 Table S2, human column)."
  reference   <- "Hosseini I, Fleisher B, Getz J, Decalf J, Kwong M, Ovacik M, Bainbridge TW, Moussion C, Rao GK, Gadkar K, Kamath AV, Ramanujan S. A Minimal PBPK/PD Model with Expansion-Enhanced Target-Mediated Drug Disposition to Support a First-in-Human Clinical Study Design for a FLT3L-Fc Molecule. Pharmaceutics. 2024 May 15;16(5):660. doi:10.3390/pharmaceutics16050660. PMCID PMC11125320. Structural equations from the Supplementary Materials 'Supplemental ODEs and Repeated Assignments' (SimBiology export); parameter values from Supplementary Table S2, 'Human Value' column. PD source data: Anandasabapathy 2015 Bone Marrow Transplant 50:924-930 (CDX-301) and Maraskovsky 2000 Blood 96:878-884 (recombinant human FLT3L); GS-3583 PK from Rajakumaraswamy 2021 J Clin Oncol 39:2559."
  vignette    <- "Hosseini_2024_flt3lfc"

  # Two routes, two molecules: FLT3L-Fc is given IV straight into `plasma`
  # (absolute ug), and CDX-301 subcutaneously into `depot` (ug/kg). Declared
  # explicitly because the registry's automatic detection only recognises
  # `depot` / `central`, and `central` here is the CDX-301 compartment, which
  # is reached only through `depot` and is never dosed directly.
  dosing <- c("plasma", "depot")

  # complex_sb / complex_db are the single-bound (drug:FLT3) and double-bound
  # (FLT3:drug:FLT3 homodimer) receptor complexes of the expansion-enhanced
  # TMDD mechanism (the paper's CRp and CRpR). cdc1 / cdc2 are the circulating
  # conventional dendritic cell counts, each fed by its own two-compartment
  # signal-transit chain, so the transit states carry the cell-type suffix.
  paper_specific_compartments <- c(
    "complex_sb", "complex_db",
    "transit1_cdc1", "transit2_cdc1", "transit1_cdc2", "transit2_cdc2",
    "cdc1", "cdc2"
  )

  units <- list(time = "h", dosing = "ug", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Note the two dose routes use different amount
  # conventions, exactly as the published model does: FLT3L-Fc is dosed into
  # `plasma` as an absolute amount (ug), whereas the CDX-301 sub-model is
  # body-weight normalised and is dosed into `depot` as ug/kg.
  compartmentData <- list(
    plasma        = list(analyte = "FLT3L-Fc", units = "ug", specimen = "plasma", verified = TRUE),
    tight         = list(analyte = "FLT3L-Fc", units = "ug", specimen = "tissue", verified = TRUE),
    leaky         = list(analyte = "FLT3L-Fc", units = "ug", specimen = "tissue", verified = TRUE),
    lymph         = list(analyte = "FLT3L-Fc", units = "ug", specimen = "lymph", verified = TRUE),
    target        = list(analyte = "free FLT3 receptor", units = "nM", specimen = "plasma", verified = TRUE),
    complex_sb    = list(analyte = "single-bound FLT3L-Fc:FLT3 complex", units = "nM", specimen = "plasma", verified = TRUE),
    complex_db    = list(analyte = "double-bound FLT3:FLT3L-Fc:FLT3 homodimer complex", units = "nM", specimen = "plasma", verified = TRUE),
    depot         = list(analyte = "CDX-301 (recombinant human FLT3 ligand)", units = "ug/kg", specimen = "administration site", verified = TRUE),
    central       = list(analyte = "CDX-301 (recombinant human FLT3 ligand)", units = "ug/kg", specimen = "plasma", verified = TRUE),
    transit1_cdc1 = list(analyte = "cDC1 expansion signal, first transit state", units = "nM", specimen = "not applicable", verified = TRUE),
    transit2_cdc1 = list(analyte = "cDC1 expansion signal, second transit state", units = "nM", specimen = "not applicable", verified = TRUE),
    transit1_cdc2 = list(analyte = "cDC2 expansion signal, first transit state", units = "nM", specimen = "not applicable", verified = TRUE),
    transit2_cdc2 = list(analyte = "cDC2 expansion signal, second transit state", units = "nM", specimen = "not applicable", verified = TRUE),
    cdc1          = list(analyte = "conventional type 1 dendritic cells", units = "cells/mL", specimen = "whole blood", verified = TRUE),
    cdc2          = list(analyte = "conventional type 2 dendritic cells", units = "cells/mL", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species       = "human",
    n_studies     = 3,
    disease_state = "Healthy adult volunteers",
    dose_range    = paste0(
      "PD calibration -- CDX-301 daily SC 3, 10, 25 or 75 ug/kg for 5 days, ",
      "25 ug/kg for 7 or 10 days (clinical study 1) and recombinant human ",
      "FLT3L daily SC 10-100 ug/kg for 14 days (clinical study 2). ",
      "FLT3L-Fc target-parameter refinement -- GS-3583 single IV 225 and 675 ug ",
      "(clinical study 3). FLT3L-Fc projections -- IV 0.01-1 mg/kg q3w; the ",
      "selected first-in-human dose is 700 ug (0.01 mg/kg for a 70 kg subject)."
    ),
    notes         = paste0(
      "Typical-value mechanistic (QSP) simulator: the paper calibrated the model ",
      "by particle-swarm optimisation in SimBiology/gQSPSim and reports point ",
      "estimates with no uncertainty, so no inter-individual variability and no ",
      "residual-error model are encoded. The PD arm was fitted to digitized ",
      "group-mean time courses, not to individual data. FLT3L-Fc and CDX-301 are ",
      "assumed to elicit the same exposure-response on cDC expansion, so the ",
      "molar concentration of either can drive the PD block; `fcflag` selects ",
      "which (0 = CDX-301, the calibration configuration; 1 = FLT3L-Fc, the ",
      "projection configuration). Table S2's human column is the human-derived ",
      "scenario, in which the FLT3 target parameters were refined against ",
      "GS-3583 PK; the cyno-derived scenario is the same structure with ",
      "rbase_target = 11.3 nM, vmax_prolif = 1.90e-3 /h, km_prolif = 0.547 and ",
      "hill_prolif = 3.93 (Table S2 cyno column) substituted."
    ),
    model_class   = "QSP / minimal PBPK with expansion-enhanced TMDD, coupled to a transit-delayed indirect-response dendritic cell expansion model (16 states)",
    n_states      = 16
  )

  ini({
    # ---- Human physiology (Hosseini 2024 Table S2, Human column) ----
    vplasma     <- fixed(2.6);     label("Plasma volume (L)")                                          # Table S2 Vplasma, Cao 2013
    vlymph      <- fixed(5.2);     label("Lymph volume (L)")                                           # Table S2 Vlymph, Cao 2013
    visf        <- fixed(15.6);    label("Total interstitial fluid volume (L)")                        # Table S2 ISF, Cao 2014; table unit label 'mL/kg' is a typo -- the SimBiology project names the parameter ISF_L, i.e. litres
    kp          <- fixed(0.8);     label("Fraction of interstitial space available to IgG1 (unitless)") # Table S2 Kp, Cao 2013
    ltot        <- fixed(0.121);   label("Total lymph flow (L/h)")                                     # Table S2 L, Cao 2013
    clp         <- fixed(6.94e-3); label("Nonspecific clearance of FLT3L-Fc from plasma (L/h)")        # Table S2 CLp, allometrically scaled from the cynomolgus monkey fit
    sigma_tight <- fixed(0.985);   label("Vascular reflection coefficient, tight tissues (unitless)")  # Table S2 sigma_tight, fit from cynomolgus monkey data
    sigma_leaky <- fixed(0.656);   label("Vascular reflection coefficient, leaky tissues (unitless)")  # Table S2 sigma_leaky, fit from cynomolgus monkey data
    sigma_lymph <- fixed(0.2);     label("Lymphatic reflection coefficient (unitless)")                # Table S2 sigma_lymph, Cao 2013

    # ---- FLT3 target turnover, binding and expansion (Table S2, Human column;
    #      human-derived scenario refined against GS-3583 PK) ----
    rbase_target <- fixed(4);      label("Baseline free FLT3 receptor concentration in plasma (nM)")   # Table S2 CR,0 human; 2.8-fold lower than the cynomolgus monkey value (Discussion)
    kdeg         <- fixed(0.016);  label("Turnover rate constant of FLT3 receptors and of all drug-receptor complexes (1/h)") # Table S2 kdeg; supplement repeated assignment kint = kdeg_central
    kon1         <- fixed(0.1);    label("Association rate constant, drug to first FLT3 receptor (1/(nM*h))")  # Table S2 kon1, assumption
    kd1          <- fixed(9);      label("Equilibrium dissociation constant, drug to first FLT3 receptor (nM)") # Table S2 KD1, Verstraete 2011
    kon2         <- fixed(3.6);    label("Association rate constant, single-bound complex to second FLT3 receptor (1/(nM*h))") # Table S2 kon2, in-house Biacore
    kd2          <- fixed(0.2);    label("Equilibrium dissociation constant, single-bound complex to second FLT3 receptor (nM)") # Table S2 KD2, in-house Biacore
    vmax_prolif  <- fixed(3.30e-3); label("Maximum rate of FLT3 receptor pool expansion driven by the double-bound complex (1/h)") # Table S2 vmprolif human
    km_prolif    <- fixed(0.02);   label("Double-bound receptor occupancy giving half-maximal expansion (fraction of total receptor)") # Table S2 kmprolif human; the supplement ODE compares km_prolif with RO_DB/100, so it is on the fraction scale despite Table S2 wording it as a percentage
    hill_prolif  <- fixed(1);      label("Hill coefficient of the receptor-expansion function (unitless)") # Table S2 alpha human

    # ---- CDX-301 empirical PK (Table S2, Human column) ----
    ka_cdx   <- fixed(0.04);   label("First-order SC absorption rate constant of CDX-301 (1/h)")      # Table S2 kabsCDX, fit to Anandasabapathy 2015 data
    vd_cdx   <- fixed(213);    label("Central volume of distribution of CDX-301 (mL/kg)")             # Table S2 VdCDX, fit to Anandasabapathy 2015 data
    vmax_cdx <- fixed(3.72);   label("Maximum nonlinear elimination rate of CDX-301 (ug/h/kg)")       # Table S2 vmCDX, fit to Anandasabapathy 2015 data; table unit label 'mL/h/kg' is a typo -- the supplement ODE names the parameter Vm_CDX_ughkg and uses it as an amount rate
    km_cdx   <- fixed(0.394);  label("CDX-301 concentration at half-maximal elimination (ug/mL)")     # Table S2 kmCDX, fit to Anandasabapathy 2015 data

    # ---- Dendritic cell expansion PD (Table S2, Human column) ----
    rbase_cdc1 <- fixed(1180);    label("Baseline cDC1 count in peripheral blood (cells/mL)")          # Table S2 initDC1, estimated from Anandasabapathy 2015
    kdeg_cdc1  <- fixed(1.60e-2); label("First-order turnover rate constant of the cDC1 pool (1/h)")   # Table S2 kdegDC1
    vmax_cdc1  <- fixed(1.61);    label("Maximum fractional rate of drug-induced cDC1 expansion (unitless multiplier on kdeg_cdc1)") # Table S2 vm1DC1
    km_cdc1    <- fixed(0.072);   label("Transit-signal concentration giving half-maximal cDC1 expansion (nM)") # Table S2 km1DC1
    ktr_cdc1   <- fixed(0.114);   label("Transit rate constant of the cDC1 signal chain (1/h)")        # Table S2 delDC1
    hill_cdc1  <- fixed(5.85);    label("Hill coefficient of the drug effect on cDC1 expansion (unitless)") # Table S2 n1DC1

    rbase_cdc2 <- fixed(12700);   label("Baseline cDC2 count in peripheral blood (cells/mL)")          # Table S2 initDC2, estimated from Maraskovsky 2000
    kdeg_cdc2  <- fixed(1.70e-3); label("First-order turnover rate constant of the cDC2 pool (1/h)")   # Table S2 kdegDC2
    vmax_cdc2  <- fixed(16.1);    label("Maximum fractional rate of drug-induced cDC2 expansion (unitless multiplier on kdeg_cdc2)") # Table S2 vm1DC2
    km_cdc2    <- fixed(0.209);   label("Transit-signal concentration giving half-maximal cDC2 expansion (nM)") # Table S2 km1DC2
    ktr_cdc2   <- fixed(0.053);   label("Transit rate constant of the cDC2 signal chain (1/h)")        # Table S2 delDC2
    hill_cdc2  <- fixed(0.888);   label("Hill coefficient of the drug effect on cDC2 expansion (unitless)") # Table S2 n1DC2

    kdeg2_dc <- fixed(3.64e-4);   label("Second-order apoptosis rate constant of expanded cDC2 (1/h)") # Table S2 kdeg2DC
    f_cdc1   <- fixed(0.476);     label("cDC1 second-order apoptosis rate as a fraction of the cDC2 value (unitless)") # Table S2 fDC1

    # ---- Configuration switch (supplement repeated assignment C_CDX_FC) ----
    fcflag <- fixed(0);           label("PD driver selector: 0 drives cDC expansion from CDX-301, 1 from FLT3L-Fc (unitless)")
  })

  model({
    # Molecular weights: FLT3L-Fc 83000 Da = 83 ug/nmol (Table S2 MWFC,
    # in-house measurement); CDX-301 35000 Da = 35 ug/nmol (Table S2 MWCDX,
    # Satyamitra 2020).
    mw_ug_nmol     <- 83
    mw_cdx_ug_nmol <- 35

    # Off-rates from the reported affinities. The supplement writes the
    # dissociation flux as kon*KD*complex, so koff = kon * KD.
    koff1 <- kon1 * kd1
    koff2 <- kon2 * kd2

    # Complex internalisation equals free-receptor turnover
    # (supplement repeated assignment kint_1h = kdeg_central_1h).
    kint <- kdeg

    # Zero-order receptor synthesis holding the free receptor at baseline in
    # the absence of drug (ksyn = kdeg * CR,0).
    ksyn <- kdeg * rbase_target

    # Cao 2013 mPBPK tissue partition. The tight tissue group holds 65% of the
    # available interstitial volume but receives only 33% of lymph flow; the
    # leaky group is the complement. Table S2's "fleaky = 0.65 / ftight = 0.35"
    # rows carry the volume split under a flow-fraction caption; the
    # supplement's SimBiology project settles it by storing these four
    # assignments verbatim -- "Vtight_L = 0.65*ISF_L*Kp",
    # "Vleaky_L = 0.35*ISF_L*Kp", "L_tight_Lh = 0.33*L_Lh" and
    # "L_leaky_Lh = 0.67*L_Lh" -- which is the Cao 2013 convention already used
    # by Yuan_2019_concizumab.
    vtight <- 0.65 * visf * kp
    vleaky <- 0.35 * visf * kp
    ltight <- 0.33 * ltot
    lleaky <- 0.67 * ltot

    # Concentrations from amounts
    cp_ugL     <- plasma / vplasma
    ctight_ugL <- tight / vtight
    cleaky_ugL <- leaky / vleaky
    clymph_ugL <- lymph / vlymph
    cp_nM      <- cp_ugL / mw_ug_nmol

    # CDX-301 is body-weight normalised: amounts in ug/kg over a volume in
    # mL/kg give ug/mL directly.
    ccdx_ugmL <- central / vd_cdx
    ccdx_nM   <- ccdx_ugmL * 1000 / mw_cdx_ug_nmol

    # Mass-action binding fluxes (nM/h)
    bind1 <- kon1 * cp_nM * target - koff1 * complex_sb
    bind2 <- kon2 * complex_sb * target - koff2 * complex_db

    # Receptor occupancy. Each double-bound complex sequesters two receptors,
    # so it is counted twice in the total pool (Eq 6).
    target_tot <- target + complex_sb + 2 * complex_db
    ro_sb <- complex_sb / max(target_tot, 1e-18) * 100                    # Eq 7
    ro_db <- 2 * complex_db / max(target_tot, 1e-18) * 100                # Eq 5
    ro    <- (1 - target / max(target_tot, 1e-18)) * 100                  # Eq 8

    # Fractional double-bound occupancy driving expansion, floored at zero so
    # solver round-off cannot raise a negative base to a fractional power.
    rodb_frac <- max(ro_db / 100, 0)
    expansion <- target_tot * vmax_prolif * rodb_frac^hill_prolif /
      (rodb_frac^hill_prolif + km_prolif^hill_prolif)

    # ---- FLT3L-Fc disposition, minimal PBPK (Eq 1 and supplement ODEs) ----
    d/dt(plasma) <- -bind1 * vplasma * mw_ug_nmol -
      (1 - sigma_tight) * ltight * cp_ugL -
      (1 - sigma_leaky) * lleaky * cp_ugL +
      ltot * clymph_ugL -
      clp * cp_ugL

    d/dt(tight) <- (1 - sigma_tight) * ltight * cp_ugL -
      (1 - sigma_lymph) * ltight * ctight_ugL

    d/dt(leaky) <- (1 - sigma_leaky) * lleaky * cp_ugL -
      (1 - sigma_lymph) * lleaky * cleaky_ugL

    d/dt(lymph) <- (1 - sigma_lymph) * ltight * ctight_ugL +
      (1 - sigma_lymph) * lleaky * cleaky_ugL -
      ltot * clymph_ugL

    # ---- FLT3 target and complexes in plasma (Eq 2-4) ----
    d/dt(target)     <- ksyn - kdeg * target - bind1 - bind2 + expansion
    d/dt(complex_sb) <- bind1 - bind2 - kint * complex_sb
    d/dt(complex_db) <- bind2 - kint * complex_db

    # ---- CDX-301 one-compartment SC model with Michaelis-Menten elimination
    #      (Eq 9-11); amounts in ug/kg ----
    d/dt(depot)   <- -ka_cdx * depot
    d/dt(central) <- ka_cdx * depot - vmax_cdx * ccdx_ugmL / (km_cdx + ccdx_ugmL)

    # ---- Shared dendritic cell expansion PD (Eq 12) ----
    # Either molecule can drive the PD block; fcflag selects which.
    cdrive <- (1 - fcflag) * ccdx_nM + fcflag * cp_nM

    d/dt(transit1_cdc1) <- ktr_cdc1 * (cdrive - transit1_cdc1)
    d/dt(transit2_cdc1) <- ktr_cdc1 * (transit1_cdc1 - transit2_cdc1)
    d/dt(transit1_cdc2) <- ktr_cdc2 * (cdrive - transit1_cdc2)
    d/dt(transit2_cdc2) <- ktr_cdc2 * (transit1_cdc2 - transit2_cdc2)

    # The supplement floors the driving signal at 1e-6 nM before the power
    # (repeated assignments C3_DC1 / C3_DC2) so a zero base never meets a
    # fractional Hill exponent.
    c3_cdc1 <- max(transit2_cdc1, 1e-6)
    c3_cdc2 <- max(transit2_cdc2, 1e-6)

    stim_cdc1 <- vmax_cdc1 * c3_cdc1^hill_cdc1 / (km_cdc1^hill_cdc1 + c3_cdc1^hill_cdc1)
    stim_cdc2 <- vmax_cdc2 * c3_cdc2^hill_cdc2 / (km_cdc2^hill_cdc2 + c3_cdc2^hill_cdc2)

    d/dt(cdc1) <- kdeg_cdc1 * rbase_cdc1 -
      kdeg_cdc1 * cdc1 +
      cdc1 * kdeg_cdc1 * stim_cdc1 -
      f_cdc1 * kdeg2_dc * max(cdc1 - rbase_cdc1, 0)^2 / rbase_cdc1

    d/dt(cdc2) <- kdeg_cdc2 * rbase_cdc2 -
      kdeg_cdc2 * cdc2 +
      cdc2 * kdeg_cdc2 * stim_cdc2 -
      kdeg2_dc * max(cdc2 - rbase_cdc2, 0)^2 / rbase_cdc2

    # ---- Initial conditions ----
    target(0) <- rbase_target
    cdc1(0)   <- rbase_cdc1
    cdc2(0)   <- rbase_cdc2

    # ---- Outputs ----
    Cc          <- cp_ugL / 1000                              # FLT3L-Fc, ug/mL
    Ccdx        <- ccdx_ugmL                                  # CDX-301, ug/mL
    DCtotal     <- cdc1 + cdc2                                # cells/mL
    DC1fold     <- cdc1 / rbase_cdc1
    DCtotalfold <- (cdc1 + cdc2) / (rbase_cdc1 + rbase_cdc2)
    RO          <- ro
    ROsb        <- ro_sb
    ROdb        <- ro_db
  })
}
