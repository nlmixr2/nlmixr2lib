Chiu_2014_tce_mouse_pbpk <- function() {
  description <- paste(
    "Preclinical (mouse). PBPK (whole-body, multi-metabolite, flow-limited",
    "tissue groups) recoded from the authors' published MCSim model code for",
    "the oral-gavage scenario. Toxicokinetics of trichloroethylene (TCE) and",
    "its oxidative and glutathione-conjugation metabolites after a single oral",
    "gavage dose in the mouse. The parent TCE disposes through stomach and",
    "duodenum absorption into seven flow-limited tissues (rapidly perfused,",
    "slowly perfused, fat, gut, liver, kidney, blood); arterial blood is the",
    "steady-state lung equation, so volatile first-pass pulmonary elimination",
    "is captured without a respiratory-tract compartment. Hepatic TCE is split",
    "by saturable (Michaelis-Menten) CYP oxidation into TCA, trichloroethanol",
    "(TCOH) and a direct-dichloroacetic-acid (DCA) branch, plus a saturable",
    "glutathione(GSH)-conjugation branch producing S-(1,2-dichlorovinyl)-",
    "glutathione (DCVG); renal GSH conjugation adds a parallel kidney source.",
    "TCOH feeds a two-tissue submodel with onward oxidation to TCA and",
    "glucuronidation to trichloroethanol-glucuronide (TCOG, with biliary",
    "recirculation); TCA is a protein-binding three-tissue submodel; DCA is a",
    "two-compartment submodel; DCVG converts to S-(1,2-dichlorovinyl)-cysteine",
    "(DCVC). Deterministic typical-value replication of the B6C3F1/J index",
    "strain (posterior-median population parameters); the paper's hierarchical",
    "Bayesian interstrain variability and parameter uncertainty are documented",
    "but not carried as a random-effects layer.",
    sep = " "
  )
  reference <- paste(
    "Chiu WA, Campbell JL Jr, Clewell HJ 3rd, Zhou YH, Wright FA, Guyton KZ,",
    "Rusyn I. Physiologically based pharmacokinetic (PBPK) modeling of",
    "interstrain variability in trichloroethylene metabolism in the mouse.",
    "Environ Health Perspect. 2014;122(5):456-463. doi:10.1289/ehp.1307623.",
    "PMID:24518055. PMCID:PMC4014769.",
    "Model structure and ODEs from the authors' MCSim code (Supplemental",
    "Material -- PBPK Model Code, TCE.risk.2.0.gamma.IS.pop.model), itself a",
    "modification of Hack et al. 2006 and Evans et al. 2009. Physiologic",
    "baselines and scaling relationships Supplemental Table S3; posterior-median",
    "scaling parameters Supplemental Table S6; key B6C3F1/J metabolism estimates",
    "main-text Table 2.",
    sep = " "
  )
  vignette <- "Chiu_2014_tce_mouse_pbpk"

  units <- list(
    time = "h",
    dosing = "single oral gavage of TCE into the stomach compartment (amt in mg)",
    concentration = paste(
      "TCE (Cc) mg/L in venous blood; the four metabolite outputs match",
      "Figure 4 ('in blood'): TCA (Cc_tca), DCA (Cc_dca), DCVG (Cc_dcvg) and",
      "DCVC (Cc_dcvc) in mmol/L"
    ),
    amount = "mg (DCVG amount is carried in mmol)",
    weight = "kg"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Mouse body weight. Every tissue volume is linear in WT; cardiac output,",
        "alveolar ventilation and the TCOH/TCOG/TCA/DCA/DCVG/DCVC rate constants",
        "are allometric (WT^0.75 for capacities, WT^-0.25 for first-order rate",
        "constants, per Supplemental Table S3). The default 0.0364 kg reproduces",
        "the main-text Table 2 total dose of 76.4 mg at the 2,100 mg/kg gavage",
        "('Dose (mg) [fixed] 76.4'); set the gavage amt to 2100 * WT.",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  # Every ODE state is a paper-mechanistic PBPK amount (mg, except ADCVGmol in
  # mmol) outside the canonical register: per-tissue TCE amounts, per-tissue
  # metabolite amounts, and metabolite distribution compartments. Names follow
  # the authors' MCSim state variables.
  compartmentData <- list(
    AStom = list(analyte = "Trichloroethylene", units = "mg", specimen = "tissue", verified = TRUE),
    ADuod = list(analyte = "Trichloroethylene", units = "mg", specimen = "tissue", verified = TRUE),
    ARap = list(analyte = "Trichloroethylene", units = "mg", specimen = "tissue", verified = TRUE),
    ASlw = list(analyte = "Trichloroethylene", units = "mg", specimen = "tissue", verified = TRUE),
    AFat = list(analyte = "Trichloroethylene", units = "mg", specimen = "tissue", verified = TRUE),
    AGut = list(analyte = "Trichloroethylene", units = "mg", specimen = "tissue", verified = TRUE),
    ALiv = list(analyte = "Trichloroethylene", units = "mg", specimen = "tissue", verified = TRUE),
    AKid = list(analyte = "Trichloroethylene", units = "mg", specimen = "tissue", verified = TRUE),
    ABld = list(analyte = "Trichloroethylene", units = "mg", specimen = "whole blood", verified = TRUE),
    ABodTCOH = list(analyte = "Trichloroethanol (TCOH)", units = "mg", specimen = "tissue", verified = TRUE),
    ALivTCOH = list(analyte = "Trichloroethanol (TCOH)", units = "mg", specimen = "tissue", verified = TRUE),
    ABodTCOG = list(
      analyte = "Trichloroethanol glucuronide (TCOG)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    ALivTCOG = list(
      analyte = "Trichloroethanol glucuronide (TCOG)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    ABileTCOG = list(
      analyte = "Trichloroethanol glucuronide (TCOG)",
      units = "mg",
      specimen = "bile",
      verified = TRUE
    ),
    APlasTCA = list(
      analyte = "Trichloroacetic acid (TCA)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    ABodTCA = list(
      analyte = "Trichloroacetic acid (TCA)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    ALivTCA = list(
      analyte = "Trichloroacetic acid (TCA)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    ADCA = list(
      analyte = "Dichloroacetic acid (DCA)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    ADCAper = list(
      analyte = "Dichloroacetic acid (DCA)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    ADCVGmol = list(
      analyte = "S-(1,2-dichlorovinyl)glutathione (DCVG)",
      units = "mmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    ADCVC = list(
      analyte = "S-(1,2-dichlorovinyl)-L-cysteine (DCVC)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    AMetOx = list(
      analyte = "Trichloroethylene, cumulative oxidative (CYP450) metabolism",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    AMetGSH = list(
      analyte = "Trichloroethylene, cumulative glutathione conjugation",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    AMetTCAprod = list(
      analyte = "Trichloroacetic acid, cumulative production",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    AMetDCAprod = list(
      analyte = "Dichloroacetic acid, cumulative production",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  paper_specific_compartments <- c(
    "AStom",
    "ADuod",
    "ARap",
    "ASlw",
    "AFat",
    "AGut",
    "ALiv",
    "AKid",
    "ABld",
    "ABodTCOH",
    "ALivTCOH",
    "ABodTCOG",
    "ALivTCOG",
    "ABileTCOG",
    "APlasTCA",
    "ABodTCA",
    "ALivTCA",
    "ADCA",
    "ADCAper",
    "ADCVGmol",
    "ADCVC",
    "AMetOx",
    "AMetGSH",
    "AMetTCAprod",
    "AMetDCAprod"
  )

  population <- list(
    species = "mouse (B6C3F1/J; one hybrid plus 16 inbred strains in the full analysis)",
    n_subjects = 17L,
    n_studies = 2L,
    age_range = "7-9 weeks (male)",
    weight_range = "~0.03-0.04 kg (default 0.0364 kg reproduces the 76.4 mg Table 2 dose)",
    sex_female_pct = 0,
    race_ethnicity = NA_character_,
    disease_state = "healthy (TCE metabolism / risk-assessment study)",
    dose_range = "single oral gavage of TCE 2,100 mg/kg in corn oil (10 mL/kg); serum metabolites at 2, 8 and 24 h",
    regions = "United States (toxicokinetic / risk-assessment study)",
    notes = paste(
      "This file is a deterministic typical-value replication of the B6C3F1/J",
      "index strain (the strain shown in Figure 4 and summarised in Table 2),",
      "for the oral-gavage route the mouse data use. Parameter values are the",
      "posterior medians of the population means (Supplemental Table S6); all",
      "interstrain scaling parameters (lnIS*, which are 0 for the B6C3F1/J",
      "reference) are held at 0. The authors fitted a hierarchical Bayesian",
      "population model (8 MCMC chains, 160,000 iterations) estimating both",
      "parameter uncertainty and interstrain variability (Tables S4-S6); those",
      "random-effects and uncertainty layers are documented but not carried",
      "here. The inhalation / closed-chamber / intravenous routes and the",
      "respiratory-tract and lung-oxidation compartments present in the full",
      "MCSim model are omitted because they are inactive for oral gavage.",
      "Renal GSH-conjugation parameters (VMaxKidDCVG, KMKidDCVG, FracKidDCVC)",
      "were fixed, not estimated, in the paper and are kept at their coded",
      "mouse defaults.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Fractional blood flows (fraction of cardiac output) -- mouse baselines
    # from the MCSim code (Species == 3 branch) / Supplemental Table S3. The
    # per-flow scaling parameters QFatC..QKidC are fixed at 1 (Table S3 note
    # e), so only cardiac output and the ventilation-perfusion ratio carry an
    # estimated scaling parameter (below).
    # =====================================================================
    fqc_fat <- fixed(0.07)
    label("Blood flow to fat, fraction of cardiac output (unitless)") # code Species==3 QFatCtmp 0.07; Table S3 QFatC0
    fqc_gut <- fixed(0.141)
    label("Blood flow to gut, fraction of cardiac output (unitless)") # code Species==3 QGutCtmp 0.141; Table S3 QGutC0
    fqc_liv <- fixed(0.02)
    label("Hepatic-artery blood flow, fraction of cardiac output (unitless)") # code Species==3 QLivCtmp 0.02; Table S3 QLivC0
    fqc_slw <- fixed(0.217)
    label("Blood flow to slowly perfused tissue, fraction of cardiac output (unitless)") # code Species==3 QSlwCtmp 0.217; Table S3 QSlwC0
    fqc_kid <- fixed(0.091)
    label("Blood flow to kidney, fraction of cardiac output (unitless)") # code Species==3 QKidCtmp 0.091; Table S3 QKidC0

    # =====================================================================
    # Fractional tissue volumes (fraction of body weight) -- mouse baselines
    # (code Species == 3 branch / Table S3). All ViC scaling parameters are
    # fixed at 1 (Table S3 note g). The perfused-volume fraction (0.8897) sets
    # the slowly perfused tissue as the residual.
    # =====================================================================
    fvc_fat <- fixed(0.07)
    label("Fat volume, fraction of body weight (unitless)") # code Species==3 VFat 0.07; Table S3 VFatC0
    fvc_gut <- fixed(0.049)
    label("Gut volume, fraction of body weight (unitless)") # code Species==3 VGut 0.049; Table S3 VGutC0
    fvc_liv <- fixed(0.055)
    label("Liver volume, fraction of body weight (unitless)") # code Species==3 VLiv 0.055; Table S3 VLivC0
    fvc_rap <- fixed(0.100)
    label("Rapidly perfused volume, fraction of body weight (unitless)") # code Species==3 VRap 0.100; Table S3 VRapC0
    fvc_kid <- fixed(0.017)
    label("Kidney volume, fraction of body weight (unitless)") # code Species==3 VKid 0.017; Table S3 VKidC0
    fvc_bld <- fixed(0.049)
    label("Blood volume, fraction of body weight (unitless)") # code Species==3 VBld 0.049; Table S3 VBldC0
    fvc_perf <- fixed(0.8897)
    label("Total perfused volume, fraction of body weight (unitless)") # code Species==3 VSlw total 0.8897; Table S3 VperfC0
    frac_plas <- fixed(0.52)
    label("Fraction of blood that is plasma (1 - hematocrit) (unitless)") # code Species==3 FracPlas 0.52; Table S3 FracPlas0

    # =====================================================================
    # TCE partition coefficients (tissue:blood, blood:air) -- mouse baselines
    # (code Species == 3). All TCE partition scaling parameters are fixed at
    # their defaults (Table S3 notes k-r), so these equal the coded baselines.
    # =====================================================================
    pb_tce <- fixed(15.0)
    label("TCE blood/air partition coefficient (unitless)") # code Species==3 PB 15; Table S3 PB0
    pfat_tce <- fixed(36.0)
    label("TCE fat/blood partition coefficient (unitless)") # code Species==3 PFat 36; Table S3 PFatC0
    pgut_tce <- fixed(1.9)
    label("TCE gut/blood partition coefficient (unitless)") # code Species==3 PGut 1.9; Table S3 PGutC0
    pliv_tce <- fixed(1.7)
    label("TCE liver/blood partition coefficient (unitless)") # code Species==3 PLiv 1.7; Table S3 PLivC0
    prap_tce <- fixed(1.9)
    label("TCE rapidly perfused/blood partition coefficient (unitless)") # code Species==3 PRap 1.9; Table S3 PRapC0
    pkid_tce <- fixed(2.1)
    label("TCE kidney/blood partition coefficient (unitless)") # code Species==3 PKid 2.1; Table S3 PKidC0
    pslw_tce <- fixed(2.4)
    label("TCE slowly perfused/blood partition coefficient (unitless)") # code Species==3 PSlw 2.4; Table S3 PSlwC0

    # =====================================================================
    # Oral absorption rate constants (/h). Estimated log-scaling parameters;
    # posterior medians from Table S6. kTD (fecal transfer) is 0 in the code
    # (100% absorption), so it is not carried.
    # =====================================================================
    lkTSD <- fixed(1.36)
    label("ln scaling for TCE stomach->duodenum transfer; kTSD = exp(lkTSD) (/h)") # Table S6 lnkTSD posterior median 1.36
    lkAS <- fixed(-0.0762)
    label("ln scaling for TCE stomach absorption; kAS = exp(lkAS) (/h)") # Table S6 lnkAS posterior median -0.0762
    lkAD <- fixed(-1.34)
    label("ln scaling for TCE duodenum absorption; kAD = exp(lkAD) (/h)") # Table S6 lnkAD posterior median -1.34

    # =====================================================================
    # Cardiac output and ventilation-perfusion scaling (Table S6 posterior
    # medians). QC = exp(lnQCC) * 11.6 * WT^0.75; VPR = exp(lnVPRC) * 2.5.
    # =====================================================================
    lnQCC <- fixed(0.152)
    label("ln scaling for cardiac output; QC = exp(lnQCC) * 11.6 * WT^0.75 (L/h)") # Table S6 lnQCC posterior median 0.152 (baseline QCC0 11.6)
    lnVPRC <- fixed(-0.319)
    label("ln scaling for ventilation-perfusion ratio; VPR = exp(lnVPRC) * 2.5 (unitless)") # Table S6 lnVPRC posterior median -0.319 (baseline VPR0 2.5)

    # =====================================================================
    # TCE hepatic oxidation (Michaelis-Menten). VMax = exp(lnVMaxC) * 2700 *
    # VLiv (mouse VMAX0 2700 mg/h/kg liver); KM = exp(lnKMC) * 36.
    # =====================================================================
    lnVMaxC <- fixed(-0.633)
    label("ln scaling for hepatic TCE oxidation VMax; VMax = exp(lnVMaxC) * 2700 * VLiv (mg/h)") # Table S6 lnVMaxC posterior median -0.633 (baseline 2700 mg/h/kg liver)
    lnKMC <- fixed(-2.44)
    label("ln scaling for hepatic TCE oxidation KM; KM = exp(lnKMC) * 36 (mg/L)") # Table S6 lnKMC posterior median -2.44 (baseline KM0 36 mg/L)

    # Oxidation splits. FracOther = expit(lnFracOtherC) is the fraction of
    # oxidation going to the direct-DCA branch; FracTCA (below) is the fraction
    # to TCA; the remainder goes to TCOH.
    lnFracOtherC <- fixed(-4.16)
    label("logit of oxidation fraction to the direct-DCA branch; FracOther = expit(lnFracOtherC)") # Table S6 lnFracOtherC posterior median -4.16
    lnFracTCAC <- fixed(-1.43)
    label("ln scaling for the TCA:TCOH oxidation split (baseline ratio 0.32)") # Table S6 lnFracTCAC posterior median -1.43 (logitFracTCA0 0.32)

    # =====================================================================
    # TCE hepatic GSH conjugation (to DCVG). VMaxDCVG = exp(lnVMaxDCVGC) * 300
    # * VLiv; the clearance scaling lnClDCVGC sets KMDCVG via VMaxDCVG / (VLiv
    # * 1.53 * exp(lnClDCVGC)).
    # =====================================================================
    lnVMaxDCVGC <- fixed(-4.96)
    label("ln scaling for hepatic GSH-conjugation VMax; VMaxDCVG = exp(lnVMaxDCVGC) * 300 * VLiv (mg/h)") # Table S6 lnVMaxDCVGC posterior median -4.96 (baseline 300 mg/h/kg liver)
    lnClDCVGC <- fixed(-3.63)
    label("ln scaling for hepatic GSH-conjugation clearance (baseline 1.53 L/h/kg liver)") # Table S6 lnClDCVGC posterior median -3.63 (ClDCVG0 1.53)

    # =====================================================================
    # TCOH submodel (oxidation to TCA, glucuronidation to TCOG, other loss).
    # VMaxTCOH = exp(lnVMaxTCOHC) * WT^0.75; etc.
    # =====================================================================
    lnVMaxTCOHC <- fixed(1.14)
    label("ln scaling for TCOH->TCA VMax; VMaxTCOH = exp(lnVMaxTCOHC) * WT^0.75 (mg/h)") # Table S6 lnVMaxTCOHC posterior median 1.14
    lnKMTCOH <- fixed(-0.0643)
    label("ln of TCOH->TCA KM; KMTCOH = exp(lnKMTCOH) (mg/L)") # Table S6 lnKMTCOH posterior median -0.0643
    lnVMaxGlucC <- fixed(4.63)
    label("ln scaling for TCOH glucuronidation VMax; VMaxGluc = exp(lnVMaxGlucC) * WT^0.75 (mg/h)") # Table S6 lnVMaxGlucC posterior median 4.63
    lnKMGluc <- fixed(3.41)
    label("ln of TCOH glucuronidation KM; KMGluc = exp(lnKMGluc) (mg/L)") # Table S6 lnKMGluc posterior median 3.41
    lnkMetTCOHC <- fixed(1.98)
    label("ln scaling for TCOH other-clearance; kMetTCOH = exp(lnkMetTCOHC) * WT^-0.25 (/h)") # Table S6 lnkMetTCOHC posterior median 1.98

    # =====================================================================
    # TCA submodel (urinary excretion, hepatic metabolism, plasma binding).
    # =====================================================================
    lnkUrnTCAC <- fixed(-2.11)
    label("ln scaling for TCA urinary excretion; kUrnTCA = exp(lnkUrnTCAC) * 0.6 * WT / VPlas (/h)") # Table S6 lnkUrnTCAC posterior median -2.11 (GFR_BW 0.6 L/h/kg)
    lnkMetTCAC <- fixed(-0.698)
    label("ln scaling for TCA hepatic metabolism; kMetTCA = exp(lnkMetTCAC) * WT^-0.25 (/h)") # Table S6 lnkMetTCAC posterior median -0.698
    lnPRBCPlasTCAC <- fixed(0.889)
    label("ln of the TCA red-blood-cell:plasma partition (sets blood:plasma ratio TCAPlas)") # Table S6 lnPRBCPlasTCAC posterior median 0.889 (PRBCPlasTCA0 0.5)
    lnkDissocC <- fixed(0.467)
    label("ln scaling for the TCA-protein dissociation constant; kDissoc = exp(lnkDissocC) * 107 (micromol/L)") # Table S6 lnkDissocC posterior median 0.467 (kDissoc0 107)
    lnBMaxkDC <- fixed(0.421)
    label("ln scaling for the TCA binding-site:dissociation ratio; BMax = exp(lnBMaxkDC) * 0.88 * kDissoc (micromol/L)") # Table S6 lnBMaxkDC posterior median 0.421 (BMaxkD0 0.88)

    # =====================================================================
    # TCOG submodel (biliary excretion, enterohepatic recirculation, urinary
    # excretion).
    # =====================================================================
    lnkBileC <- fixed(0.863)
    label("ln scaling for TCOG liver->bile; kBile = exp(lnkBileC) * WT^-0.25 (/h)") # Table S6 lnkBileC posterior median 0.863
    lnkEHRC <- fixed(-3.25)
    label("ln scaling for TCOG bile->TCOH recirculation; kEHR = exp(lnkEHRC) * WT^-0.25 (/h)") # Table S6 lnkEHRC posterior median -3.25
    lnkUrnTCOGC <- fixed(2.33)
    label("ln scaling for TCOG urinary excretion; kUrnTCOG = exp(lnkUrnTCOGC) * 0.6 * WT / (VBodTCOH * PBodTCOG) (/h)") # Table S6 lnkUrnTCOGC posterior median 2.33 (GFR_BW 0.6)

    # =====================================================================
    # DCVG -> DCVC and DCVC elimination.
    # =====================================================================
    lnkDCVGC <- fixed(2.19)
    label("ln scaling for DCVG->DCVC; kDCVG = exp(lnkDCVGC) * WT^-0.25 (/h)") # Table S6 lnkDCVGC posterior median 2.19
    lnkElimDCVCC <- fixed(-1.14)
    label("ln scaling for DCVC elimination; kElimDCVC = exp(lnkElimDCVCC) * WT^-0.25 (/h)") # Table S6 lnkElimDCVCC posterior median -1.14

    # =====================================================================
    # DCA two-compartment submodel (clearance, central<->peripheral transfer).
    # =====================================================================
    lnkClearDCAC <- fixed(-2.57)
    label("ln scaling for DCA clearance; kClearDCA = exp(lnkClearDCAC) * WT^-0.25 (/h)") # Table S6 lnkClearDCAC posterior median -2.57
    lnkDCAcen_perC <- fixed(1.25)
    label("ln scaling for DCA central->peripheral; kDCAcen_per = exp(lnkDCAcen_perC) * WT^-0.25 (/h)") # Table S6 lnkDCAcen_perC posterior median 1.25
    lnkDCAper_cenC <- fixed(-2.14)
    label("ln scaling for DCA peripheral->central; kDCAper_cen = exp(lnkDCAper_cenC) * WT^-0.25 (/h)") # Table S6 lnkDCAper_cenC posterior median -2.14

    # =====================================================================
    # Effective-partition scaling for the DCVG, DCVC and DCA distribution
    # volumes: V_X = VBld + exp(lnPeff_X) * (VBod + VLiv) (Table S3 note x).
    # =====================================================================
    lnPeffDCVG <- fixed(0.0903)
    label("ln effective body:blood partition for the DCVG distribution volume (unitless)") # Table S6 lnPeffDCVG posterior median 0.0903
    lnPeffDCVC <- fixed(5.98)
    label("ln effective body:blood partition for the DCVC distribution volume (unitless)") # Table S6 lnPeffDCVC posterior median 5.98
    lnPeffDCA <- fixed(3.67)
    label("ln effective body:blood partition for the DCA distribution volume (unitless)") # Table S6 lnPeffDCA posterior median 3.67

    # =====================================================================
    # Partition / metabolism constants held FIXED (not estimated) in the
    # paper. TCOH/TCOG tissue partitions (code Species == 3; Table S3 notes
    # v,w: scaling fixed at 0). Renal GSH conjugation parameters (Table S3;
    # not in Tables S4/S6, so fixed at their coded mouse defaults). Free-TCA
    # tissue partitions (Table S3 note t: scaling fixed at 0).
    # =====================================================================
    pbod_tcoh <- fixed(1.11)
    label("TCOH body/blood partition coefficient (unitless)") # code Species==3 PBodTCOH 1.11; Table S3 PBodTCOH0
    pliv_tcoh <- fixed(1.3)
    label("TCOH liver/blood partition coefficient (unitless)") # code Species==3 PLivTCOH 1.3; Table S3 PLivTCOH0
    pbod_tcog <- fixed(1.11)
    label("TCOG body/blood partition coefficient (unitless)") # code Species==3 PBodTCOG 1.11; Table S3 PBodTCOG0
    pliv_tcog <- fixed(1.3)
    label("TCOG liver/blood partition coefficient (unitless)") # code Species==3 PLivTCOG 1.3; Table S3 PLivTCOG0
    pbod_tca0 <- fixed(0.88)
    label("Free-TCA body/blood baseline partition (multiplied by TCAPlas) (unitless)") # code Species==3 PBodTCA 0.88; Table S3 PBodTCAC0
    pliv_tca0 <- fixed(1.18)
    label("Free-TCA liver/blood baseline partition (multiplied by TCAPlas) (unitless)") # code Species==3 PLivTCA 1.18; Table S3 PLivTCAC0
    vmaxkid0 <- fixed(60.0)
    label("Renal GSH-conjugation VMax per kg kidney (mg/h/kg kidney)") # code Species==3 VMaxKidDCVG baseline 60 (scaling fixed; not in Table S6)
    clkid0 <- fixed(0.34)
    label("Renal GSH-conjugation clearance per kg kidney (L/h/kg kidney)") # code Species==3 KMKidDCVG via Cl 0.34 (scaling fixed; not in Table S6)
    frackiddcvc <- fixed(0.5)
    label("Fraction of renal GSH conjugation going directly (first-pass) to DCVC (unitless)") # code FracKidDCVC = expit(0) = 0.5 (lnFracKidDCVCC FIXED)
  })

  model({
    # -------------------------------------------------------------------
    # Allometric body-size terms and molecular weights (g/mol; code
    # "Molecular Weights" block). Stoichiometry converts a mass rate of the
    # precursor into a mass rate of the product.
    # -------------------------------------------------------------------
    BW75 <- WT^0.75
    BW25 <- WT^0.25

    mwTCE <- 131.39
    mwTCA <- 163.5
    mwTCOH <- 149.5
    mwTCOHGluc <- 325.53
    mwDCA <- 129.0
    mwDCVC <- 216.1

    stochDCATCE <- mwDCA / mwTCE
    stochTCATCE <- mwTCA / mwTCE
    stochTCATCOH <- mwTCA / mwTCOH
    stochTCOHTCE <- mwTCOH / mwTCE
    stochGlucTCOH <- mwTCOHGluc / mwTCOH
    stochTCOHGluc <- mwTCOH / mwTCOHGluc
    stochDCVCTCE <- mwDCVC / mwTCE

    # -------------------------------------------------------------------
    # Tissue volumes (L). VSlw is the perfused residual; VBod (TCA body) and
    # VBodTCOH (TCOH/TCOG body) are lumped sums (Table S3 notes g,i,j).
    # -------------------------------------------------------------------
    VFat <- fvc_fat * WT
    VGut <- fvc_gut * WT
    VLiv <- fvc_liv * WT
    VRap <- fvc_rap * WT
    VKid <- fvc_kid * WT
    VBld <- fvc_bld * WT
    VSlw <- fvc_perf * WT - (VFat + VGut + VLiv + VRap + VKid + VBld)
    VPlas <- frac_plas * VBld
    VBod <- VFat + VGut + VRap + VKid + VSlw
    VBodTCOH <- VBod + VBld

    # -------------------------------------------------------------------
    # Cardiac output (L/h), tissue blood flows, and plasma flows. QRap is the
    # flow residual; QBod is cardiac output minus the gut+liver (portal+
    # arterial) flow (code "Blood Flows to Tissues").
    # -------------------------------------------------------------------
    QC <- exp(lnQCC) * 11.6 * BW75
    VPR <- exp(lnVPRC) * 2.5
    QP <- QC * VPR
    QFat <- fqc_fat * QC
    QGut <- fqc_gut * QC
    QLiv <- fqc_liv * QC
    QSlw <- fqc_slw * QC
    QKid <- fqc_kid * QC
    QRap <- QC - (QFat + QGut + QLiv + QSlw + QKid)
    QGutLiv <- QGut + QLiv
    QBod <- QC - QGutLiv
    QCPlas <- frac_plas * QC
    QBodPlas <- frac_plas * QBod
    QGutLivPlas <- frac_plas * QGutLiv

    # -------------------------------------------------------------------
    # Partition coefficients for the metabolite submodels.
    # -------------------------------------------------------------------
    TCAPlas <- frac_plas + (1 - frac_plas) * 0.5 * exp(lnPRBCPlasTCAC)
    PBodTCA <- TCAPlas * pbod_tca0
    PLivTCA <- TCAPlas * pliv_tca0
    kDissoc <- exp(lnkDissocC) * 107
    BMax <- kDissoc * exp(lnBMaxkDC) * 0.88

    # Effective distribution volumes for DCVG, DCVC, DCA (Table S3 note x).
    VDCVG <- VBld + exp(lnPeffDCVG) * (VBod + VLiv)
    VDCVC <- VBld + exp(lnPeffDCVC) * (VBod + VLiv)
    VDCA <- VBld + exp(lnPeffDCA) * (VBod + VLiv)

    # -------------------------------------------------------------------
    # Metabolic capacities and rate constants (code "TCE ... Metabolism
    # Constants"). Interstrain scaling is 1 for the B6C3F1/J index strain.
    # -------------------------------------------------------------------
    VMax <- VLiv * exp(lnVMaxC) * 2700
    KM <- 36 * exp(lnKMC)
    FracOther <- expit(lnFracOtherC)
    FracTCA <- 0.32 * exp(lnFracTCAC) * (1 - FracOther) / (1 + 0.32 * exp(lnFracTCAC))

    VMaxDCVG <- VLiv * 300 * exp(lnVMaxDCVGC)
    KMDCVG <- VMaxDCVG / (VLiv * 1.53 * exp(lnClDCVGC))
    VMaxKidDCVG <- VKid * vmaxkid0
    KMKidDCVG <- VMaxKidDCVG / (VKid * clkid0)

    VMaxTCOH <- BW75 * exp(lnVMaxTCOHC)
    KMTCOH <- exp(lnKMTCOH)
    VMaxGluc <- BW75 * exp(lnVMaxGlucC)
    KMGluc <- exp(lnKMGluc)
    kMetTCOH <- exp(lnkMetTCOHC) / BW25

    kUrnTCA <- exp(lnkUrnTCAC) * 0.6 * WT / VPlas
    kMetTCA <- exp(lnkMetTCAC) / BW25
    kBile <- exp(lnkBileC) / BW25
    kEHR <- exp(lnkEHRC) / BW25
    kUrnTCOG <- exp(lnkUrnTCOGC) * 0.6 * WT / (VBodTCOH * pbod_tcog)
    kDCVG <- exp(lnkDCVGC) / BW25
    kElimDCVC <- exp(lnkElimDCVCC) / BW25
    kClearDCA <- exp(lnkClearDCAC) / BW25
    kDCAcen_per <- exp(lnkDCAcen_perC) / BW25
    kDCAper_cen <- exp(lnkDCAper_cenC) / BW25

    kTSD <- exp(lkTSD)
    kAS <- exp(lkAS)
    kAD <- exp(lkAD)

    # -------------------------------------------------------------------
    # TCE tissue concentrations and venous effluents (C / partition).
    # -------------------------------------------------------------------
    CRap <- ARap / VRap
    CSlw <- ASlw / VSlw
    CFat <- AFat / VFat
    CGut <- AGut / VGut
    CLiv <- ALiv / VLiv
    CKid <- AKid / VKid
    CVRap <- CRap / prap_tce
    CVSlw <- CSlw / pslw_tce
    CVFat <- CFat / pfat_tce
    CVGut <- CGut / pgut_tce
    CVLiv <- CLiv / pliv_tce
    CVKid <- CKid / pkid_tce
    CVen <- ABld / VBld

    # Steady-state lung: arterial TCE is reduced from mixed-venous by volatile
    # first-pass pulmonary exhalation (the QP/PB term); no inhalation source
    # for oral gavage.
    CArt <- QC * CVen / (QC + QP / pb_tce)

    # -------------------------------------------------------------------
    # Oral absorption of TCE (stomach -> duodenum -> gut compartment). kTD
    # (fecal) is 0 so 100% is absorbed.
    # -------------------------------------------------------------------
    RAO <- kAS * AStom + kAD * ADuod
    d/dt(AStom) <- -AStom * (kAS + kTSD)
    d/dt(ADuod) <- kTSD * AStom - kAD * ADuod

    # -------------------------------------------------------------------
    # TCE saturable hepatic oxidation and GSH conjugation, and renal GSH
    # conjugation.
    # -------------------------------------------------------------------
    RAMetLiv1 <- VMax * CVLiv / (KM + CVLiv)
    RAMetLiv2 <- VMaxDCVG * CVLiv / (KMDCVG + CVLiv)
    RAMetKid <- VMaxKidDCVG * CVKid / (KMKidDCVG + CVKid)

    # -------------------------------------------------------------------
    # TCE flow-limited tissue ODEs (mg). Blood collects the venous effluents;
    # gut receives the oral absorption; liver carries both metabolism routes;
    # kidney carries the renal GSH conjugation.
    # -------------------------------------------------------------------
    d/dt(ABld) <- QFat * CVFat + QGutLiv * CVLiv + QSlw * CVSlw +
      QRap * CVRap + QKid * CVKid - CVen * QC
    d/dt(ARap) <- QRap * (CArt - CVRap)
    d/dt(ASlw) <- QSlw * (CArt - CVSlw)
    d/dt(AFat) <- QFat * (CArt - CVFat)
    d/dt(AGut) <- QGut * (CArt - CVGut) + RAO
    d/dt(ALiv) <- QLiv * (CArt - CVLiv) + QGut * (CVGut - CVLiv) -
      RAMetLiv1 - RAMetLiv2
    d/dt(AKid) <- QKid * (CArt - CVKid) - RAMetKid

    # -------------------------------------------------------------------
    # TCOH submodel. Produced from the non-TCA, non-other fraction of hepatic
    # oxidation; cleared by onward oxidation to TCA, glucuronidation to TCOG,
    # and a lumped "other" route; receives recirculated TCOG from bile.
    # -------------------------------------------------------------------
    CVBodTCOH <- ABodTCOH / VBodTCOH / pbod_tcoh
    CVLivTCOH <- ALivTCOH / VLiv / pliv_tcoh
    CTCOH <- (QBod * CVBodTCOH + QGutLiv * CVLivTCOH) / QC
    RAMetTCOHTCA <- VMaxTCOH * CVLivTCOH / (KMTCOH + CVLivTCOH)
    RAMetTCOHGluc <- VMaxGluc * CVLivTCOH / (KMGluc + CVLivTCOH)
    RAMetTCOH <- kMetTCOH * ALivTCOH
    RARecircTCOG <- kEHR * ABileTCOG
    d/dt(ABodTCOH) <- QBod * (CTCOH - CVBodTCOH)
    d/dt(ALivTCOH) <- QGutLiv * (CTCOH - CVLivTCOH) - RAMetTCOH -
      RAMetTCOHTCA - RAMetTCOHGluc +
      (1 - FracOther - FracTCA) * stochTCOHTCE * RAMetLiv1 +
      stochTCOHGluc * RARecircTCOG

    # -------------------------------------------------------------------
    # DCA two-compartment submodel. Produced from the "other" oxidation
    # fraction; cleared from the central compartment; exchanges with a
    # peripheral compartment.
    # -------------------------------------------------------------------
    d/dt(ADCA) <- FracOther * stochDCATCE * RAMetLiv1 - kClearDCA * ADCA -
      kDCAcen_per * ADCA + kDCAper_cen * ADCAper
    d/dt(ADCAper) <- kDCAcen_per * ADCA - kDCAper_cen * ADCAper

    # -------------------------------------------------------------------
    # TCA protein-binding submodel. Free plasma TCA from the binding
    # quadratic; tissue venous concentrations are bound (in equilibrium with
    # plasma) plus free (in equilibrium with tissue).
    # -------------------------------------------------------------------
    CPlasTCA <- APlasTCA / VPlas
    CPlasTCAMole <- CPlasTCA / mwTCA * 1000
    aTCA <- kDissoc + BMax - CPlasTCAMole
    bTCA <- 4 * kDissoc * CPlasTCAMole
    CPlasTCAFreeMole <- 0.5 * (sqrt(aTCA * aTCA + bTCA) - aTCA)
    CPlasTCAFree <- CPlasTCAFreeMole * mwTCA / 1000
    APlasTCAFree <- CPlasTCAFree * VPlas
    CPlasTCABnd <- CPlasTCA - CPlasTCAFree
    CBodTCA <- ABodTCA / VBod
    CLivTCA <- ALivTCA / VLiv
    CVBodTCAFree <- CBodTCA / PBodTCA
    CVBodTCA <- CPlasTCABnd + CVBodTCAFree
    CVLivTCAFree <- CLivTCA / PLivTCA
    CVLivTCA <- CPlasTCABnd + CVLivTCAFree
    RUrnTCA <- kUrnTCA * APlasTCAFree
    RAMetTCA <- kMetTCA * ALivTCA
    d/dt(APlasTCA) <- QBodPlas * CVBodTCA + QGutLivPlas * CVLivTCA -
      QCPlas * CPlasTCA - RUrnTCA
    d/dt(ABodTCA) <- QBodPlas * (CPlasTCAFree - CVBodTCAFree)
    d/dt(ALivTCA) <- QGutLivPlas * (CPlasTCAFree - CVLivTCAFree) - RAMetTCA +
      FracTCA * stochTCATCE * RAMetLiv1 + stochTCATCOH * RAMetTCOHTCA

    # -------------------------------------------------------------------
    # TCOG submodel. Produced in liver by glucuronidation of TCOH; excreted
    # in urine from the body, in bile from the liver; bile recirculates to
    # hepatic TCOH.
    # -------------------------------------------------------------------
    CVBodTCOG <- ABodTCOG / VBodTCOH / pbod_tcog
    CVLivTCOG <- ALivTCOG / VLiv / pliv_tcog
    CTCOG <- (QBod * CVBodTCOG + QGutLiv * CVLivTCOG) / QC
    RUrnTCOG <- kUrnTCOG * ABodTCOG
    RBileTCOG <- kBile * ALivTCOG
    d/dt(ABodTCOG) <- QBod * (CTCOG - CVBodTCOG) - RUrnTCOG
    d/dt(ALivTCOG) <- QGutLiv * (CTCOG - CVLivTCOG) +
      stochGlucTCOH * RAMetTCOHGluc - RBileTCOG
    d/dt(ABileTCOG) <- RBileTCOG - RARecircTCOG

    # -------------------------------------------------------------------
    # DCVG submodel (amount in mmol). Produced from hepatic GSH conjugation
    # and the systemic part of renal GSH conjugation; converted to DCVC.
    # -------------------------------------------------------------------
    RAMetDCVGmol <- kDCVG * ADCVGmol
    d/dt(ADCVGmol) <- (RAMetLiv2 + RAMetKid * (1 - frackiddcvc)) / mwTCE -
      RAMetDCVGmol

    # -------------------------------------------------------------------
    # DCVC submodel (amount in mg). Produced from DCVG conversion and the
    # first-pass renal GSH fraction; cleared by a lumped elimination route.
    # -------------------------------------------------------------------
    d/dt(ADCVC) <- RAMetDCVGmol * mwDCVC +
      RAMetKid * frackiddcvc * stochDCVCTCE - kElimDCVC * ADCVC

    # -------------------------------------------------------------------
    # Cumulative metabolism dose-metrics (mg), the published code's AMetOX,
    # AMetGSH, TotTCAProd and TotDCAProd outputs (Table 2). Pure accumulators
    # (lung oxidation is zero for oral gavage, so AMetOX is hepatic only).
    # -------------------------------------------------------------------
    d/dt(AMetOx) <- RAMetLiv1
    d/dt(AMetGSH) <- RAMetLiv2 + RAMetKid
    d/dt(AMetTCAprod) <- FracTCA * stochTCATCE * RAMetLiv1 +
      stochTCATCOH * RAMetTCOHTCA
    d/dt(AMetDCAprod) <- FracOther * stochDCATCE * RAMetLiv1

    # -------------------------------------------------------------------
    # Observations. Cc is venous blood TCE (mg/L). The four metabolite
    # outputs replicate Figure 4, all "in blood (mmol/L)": TCA in blood is
    # the plasma concentration times the blood:plasma ratio (TCAPlas),
    # converted to mmol/L; DCA, DCVG and DCVC are blood molar concentrations.
    # -------------------------------------------------------------------
    CDCAmol <- ADCA / mwDCA / VDCA
    CDCVGmol <- ADCVGmol / VDCVG
    CDCVCmol <- ADCVC / mwDCVC / VDCVC

    Cc <- CVen
    Cc_tca <- CPlasTCA * TCAPlas / mwTCA
    Cc_dca <- CDCAmol
    Cc_dcvg <- CDCVGmol
    Cc_dcvc <- CDCVCmol
  })
}
