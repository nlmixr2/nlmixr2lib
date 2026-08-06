Zhao_2024_ciprofloxacin_colistin_plasma <- function() {
  description <- "Semi-mechanistic PK/PD (QSP) model predicting Escherichia coli killing in human plasma (bloodstream infection) by intravenous ciprofloxacin and colistimethate sodium, alone and in combination. The in vitro time-kill PD model of Zhao 2024 (five pre-existing subpopulations on the 2x2 ciprofloxacin/colistin susceptible-resistant grid plus a pre-existing resting subpopulation, additive killing with EC50-shifting interaction factors) is driven by unbound plasma concentrations from two published critically-ill population PK models: a two-compartment model for ciprofloxacin (Khachman 2011) and the CMS1/CMS2 prodrug plus formed-colistin model of Kristoffersson 2020. Colistin EC50 is rescaled to the simulated strain MIC. Sibling models: Zhao_2024_ciprofloxacin_colistin_invitro (the underlying time-kill fit) and Zhao_2024_ciprofloxacin_colistin_kidney (the same model driven by kidney interstitial concentrations)."
  reference <- "Zhao C, Kristoffersson AN, Khan DD, Lagerback P, Lustig U, Cao S, Annerstedt C, Cars O, Andersson DI, Hughes D, Nielsen EI, Friberg LE. Quantifying combined effects of colistin and ciprofloxacin against Escherichia coli in an in silico pharmacokinetic-pharmacodynamic model. Sci Rep. 2024 May 22;14(1):11603. doi:10.1038/s41598-024-61518-0. Pharmacodynamic structure and estimates: main text Eqs 1-18 and Table 1. Population PK parameters and the plasma-driven implementation: the authors' deposited mrgsolve model, Supplementary zip, Supplementary/PKPD_CIPCST_run422b.cpp, with the strain-specific parameter sets and dosing regimens in Supplementary/MrgSolve_CIPCST_PKPD_clean_0326.Rmd. Upstream population PK sources cited by that file: Khachman D et al, J Antimicrob Chemother 2011;66:1798-1809 (ciprofloxacin) and Kristoffersson AN et al, Clin Microbiol Infect 2020;26:1644-1650 (colistimethate/colistin). Reproduces main text Fig. 3 (plasma panels) and Fig. 4."
  vignette <- "Zhao_2024_ciprofloxacin_colistin"
  units <- list(time = "h", dosing = "mg (ciprofloxacin) or mol (colistimethate)", concentration = "mg/L")

  # Bacterial subpopulation indices follow Fig. 1 of the paper (see the sibling
  # in vitro model). cipCentral/cipPeripheral/cipKidney carry ciprofloxacin;
  # cms1Central/cms1Peripheral/cms2Central/cms2Peripheral/cstCentral carry the
  # two colistimethate species and the formed colistin; cmsKidney/cstKidney are
  # the lumped kidney tissue-plus-tubule compartments of Supplementary Eqs S4-S5
  # (retained here because the formed-colistin kidney state feeds back into
  # nothing but is the quantity the sibling kidney model uses for effect).
  paper_specific_compartments <- c(
    "cipCentral", "cipPeripheral", "cipKidney",
    "cms1Central", "cms1Peripheral", "cms2Central", "cms2Peripheral",
    "cstCentral", "cmsKidney", "cstKidney",
    "bact_s1", "bact_r1", "bact_nc1",
    "bact_s2", "bact_r2", "bact_nc2",
    "bact_s3", "bact_r3", "bact_nc3",
    "bact_s4", "bact_r4", "bact_nc4",
    "bact_pr"
  )

  # Ciprofloxacin states are in mg and the colistimethate/colistin states in
  # mol, matching the units the two upstream population PK models were built
  # in; colistin is converted to mg/L with 1163 g/mol where it drives effect.
  # The bacterial states carry specimen "not applicable": they live at the
  # site of infection, which the controlled specimen vocabulary does not name.
  compartmentData <- list(
    cipCentral = list(analyte = "ciprofloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    cipPeripheral = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    cipKidney = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    cms1Central = list(analyte = "colistimethate sodium, species CMS1", units = "mol", specimen = "plasma", verified = TRUE),
    cms1Peripheral = list(analyte = "colistimethate sodium, species CMS1", units = "mol", specimen = "tissue", verified = TRUE),
    cms2Central = list(analyte = "colistimethate sodium, species CMS2", units = "mol", specimen = "plasma", verified = TRUE),
    cms2Peripheral = list(analyte = "colistimethate sodium, species CMS2", units = "mol", specimen = "tissue", verified = TRUE),
    cstCentral = list(analyte = "colistin formed from colistimethate", units = "mol", specimen = "plasma", verified = TRUE),
    cmsKidney = list(analyte = "colistimethate sodium", units = "mol", specimen = "tissue", verified = TRUE),
    cstKidney = list(analyte = "colistin", units = "mol", specimen = "tissue", verified = TRUE),
    bact_s1 = list(analyte = "Escherichia coli, ciprofloxacin-susceptible and colistin-susceptible subpopulation, growing (S) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_r1 = list(analyte = "Escherichia coli, ciprofloxacin-susceptible and colistin-susceptible subpopulation, resting (R) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_nc1 = list(analyte = "Escherichia coli, ciprofloxacin-susceptible and colistin-susceptible subpopulation, non-colony-forming (Nc) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_s2 = list(analyte = "Escherichia coli, ciprofloxacin-susceptible and colistin-resistant subpopulation, growing (S) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_r2 = list(analyte = "Escherichia coli, ciprofloxacin-susceptible and colistin-resistant subpopulation, resting (R) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_nc2 = list(analyte = "Escherichia coli, ciprofloxacin-susceptible and colistin-resistant subpopulation, non-colony-forming (Nc) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_s3 = list(analyte = "Escherichia coli, ciprofloxacin-resistant and colistin-susceptible subpopulation, growing (S) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_r3 = list(analyte = "Escherichia coli, ciprofloxacin-resistant and colistin-susceptible subpopulation, resting (R) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_nc3 = list(analyte = "Escherichia coli, ciprofloxacin-resistant and colistin-susceptible subpopulation, non-colony-forming (Nc) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_s4 = list(analyte = "Escherichia coli, ciprofloxacin-resistant and colistin-resistant subpopulation, growing (S) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_r4 = list(analyte = "Escherichia coli, ciprofloxacin-resistant and colistin-resistant subpopulation, resting (R) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_nc4 = list(analyte = "Escherichia coli, ciprofloxacin-resistant and colistin-resistant subpopulation, non-colony-forming (Nc) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_pr = list(analyte = "Escherichia coli, pre-existing resting (persister) subpopulation, present only in strain C47", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    BACT = list(
      description = "Bacterial strain indicator selecting the strain-specific colistin PD parameter set, the strain's colistin MIC used to rescale EC50, and the fraction of pre-existing resting bacteria in the inoculum.",
      units = "(categorical)",
      type = "categorical",
      source_name = "STRAIN",
      notes = paste(
        "Levels follow the FSTRAIN coding of the authors' NONMEM dataset, reduced to the three groups that the simulation chapter of MrgSolve_CIPCST_PKPD_clean_0326.Rmd steps through (STRAIN = 'C47', 'WT', 'LM'):",
        "47 = C47, the clinical urinary isolate; the only group carrying a pre-existing resting subpopulation, and the only one for which the combination beats ciprofloxacin monotherapy at every dose because colistin also kills those bacteria.",
        "347 = LM347, the MG1655 wild type; the one group whose colistin potency is REDUCED in combination (I_CST = 1.19 > 1).",
        "421 = the LM378/LM421 laboratory-mutant group, which shares one colistin PD parameter set (Table 1 column 'LM378, LM421'). Level 378 is accepted as a synonym of 421.",
        "There is no reference level: each level selects its own parameter set rather than contributing a contrast against a baseline.",
        "The simulations deliberately use MICs far above those of the source strains (Table 1 was fitted at MIC_CIP 0.023-1 mg/L); set micCip and micCst to the extrapolated values and BACT only to pick the PD parameter set."
      )
    ),
    WT = list(
      description = "Total body weight, used for the cardiac-output and kidney-volume allometry that sets renal blood flow and kidney volume.",
      units = "kg",
      type = "continuous",
      source_name = "WT",
      notes = "Fixed to 80 kg for every simulated patient (Supplementary Methods: 'a male with a body weight of 80 kg'). At 80 kg the allometric expressions reproduce the Q_kidney = 73 L/h and V_kidney = 0.34 L quoted in the Supplementary Methods."
    ),
    CRCL = list(
      description = "Creatinine clearance, driving ciprofloxacin total clearance, colistimethate renal clearance and colistin clearance.",
      units = "mL/min",
      type = "continuous",
      source_name = "CLCR",
      notes = "RAW creatinine clearance in mL/min, NOT normalised to 1.73 m^2 body surface area, matching the parameterisation of both upstream population PK models (precedent for the un-normalised form under this canonical: Delattre_2010_amikacin.R). Fixed to 90 mL/min in every simulation; the authors note in the Discussion that using a single CrCL value makes the simulated variability narrower than a real patient population with a wide CrCL range."
    )
  )

  population <- list(
    species = "human (simulated)",
    n_subjects = 1000L,
    n_studies = 0L,
    disease_state = "Simulated critically ill adults with an Escherichia coli bloodstream infection. The bacterial strains are deliberately extrapolated beyond those used to fit the PD model: MIC_CIP 1-8 mg/L (only strains above the MIC_CIP 1 mg/L that 400 mg q8h already clears without regrowth were of interest) and MIC_CST 2 or 4 mg/L, 2 mg/L being the EUCAST epidemiological cut-off for resistance.",
    model_system = "Monte Carlo simulation of 1000 virtual patients over 32 h in mrgsolve, linking two published critically-ill population PK models to the in vitro time-kill PD model. Reported inter-individual and inter-occasion variability of the PK models were included; all patients share WT 80 kg and CrCL 90 mL/min.",
    initial_inoculum = "10^6 CFU/mL at the time of the first dose",
    dose_range = "Ciprofloxacin 400 mg as a 1 h intravenous infusion q8h in every scenario. Colistimethate sodium 9 + 4.5, 4.5 + 2 or 2 + 1 MU (loading + maintenance) as a 0.5 h intravenous infusion q12h, alone or combined with ciprofloxacin.",
    notes = paste(
      "1 MU of colistimethate sodium corresponds to 45.9 micromol; doses are therefore given in mol into cms1Central, while ciprofloxacin is dosed in mg into cipCentral (Supplementary Methods, last paragraph).",
      "The colistin plasma concentration driving the effect is the unbound formed-colistin concentration, converted from mol/L to mg/L with 1163 g/mol and multiplied by the unbound fraction 0.34. Ciprofloxacin uses an unbound fraction of 0.65.",
      "Ciprofloxacin population PK (two-compartment, CrCL on clearance) is from Khachman 2011; colistimethate/colistin population PK (two prodrug species CMS1 and CMS2 each with a peripheral compartment, feeding a single formed-colistin compartment) is from Kristoffersson 2020. Neither was re-estimated here.",
      "The delay in formation of colistin from colistimethate is what produces the transient rise in bacterial load at the start of colistin monotherapy that the paper describes.",
      "Simulated colistin EC50 is rescaled from the fitted value by the ratio of the simulated to the source-strain MIC (Eq. 18, with the power fixed to 1). The authors' deposited code does NOT apply that rescaling to the pre-existing-resting-subpopulation EC50; this file reproduces that behaviour and flags it in the vignette Errata.",
      "Sibling models: Zhao_2024_ciprofloxacin_colistin_invitro (the underlying time-kill fit) and Zhao_2024_ciprofloxacin_colistin_kidney (identical apart from the two concentrations that drive the drug effect)."
    )
  )

  ini({
    # ---- Simulation design inputs ---------------------------------------
    # Exposed as fixed ini() entries so a scenario can be set with
    # rxSolve(params = c(...)). Defaults reproduce the C47 panel of Fig. 4 at
    # the lowest simulated ciprofloxacin MIC.
    inoc   <- fixed(1e6); label("Bacterial density at the time of the first dose (CFU/mL)")                 # Fig. 4 caption: "Antibiotics were administered at a bacterial concentration of 6 log10 CFU/mL"
    micCip <- fixed(1);   label("Ciprofloxacin MIC of the simulated strain (mg/L)")                          # Predictions of clinical drug effect: strains with MIC_CIP 1-8 mg/L were simulated in plasma
    micCst <- fixed(2);   label("Colistin MIC of the simulated strain (mg/L)")                               # Predictions of clinical drug effect: MIC_CST 2 or 4 mg/L, 2 mg/L being the EUCAST epidemiological cut-off

    # ---- Ciprofloxacin population PK (Khachman 2011) --------------------
    vcCip     <- fixed(38);    label("Ciprofloxacin central volume (L)")                                     # PKPD_CIPCST_run422b.cpp $PARAM TVVcCIP : 38 : (L) cipro V1 Khachman D et al (2011)
    vpCip     <- fixed(73);    label("Ciprofloxacin peripheral volume (L)")                                  # PKPD_CIPCST_run422b.cpp $PARAM TVVpCIP : 73 : (L) cipro Vp Khachman D et al (2011)
    qCip      <- fixed(60);    label("Ciprofloxacin intercompartmental clearance (L/h)")                     # PKPD_CIPCST_run422b.cpp $PARAM TVQCIP : 60 : (L/h) cipro Q Khachman D et al (2011)
    clCip     <- fixed(18);    label("Ciprofloxacin total clearance at the reference creatinine clearance (L/h)") # PKPD_CIPCST_run422b.cpp $PARAM TVCLCIP : 18 : (L/h) cipro CLtot Khachman D et al (2011)
    crclRef   <- fixed(91.7);  label("Reference creatinine clearance for ciprofloxacin clearance (mL/min)")  # PKPD_CIPCST_run422b.cpp $MAIN: CLCIP = TVCLCIP*pow((CLCR/91.7),0.42)*...
    crclExp   <- fixed(0.42);  label("Exponent of creatinine clearance on ciprofloxacin clearance")          # PKPD_CIPCST_run422b.cpp $MAIN: CLCIP = TVCLCIP*pow((CLCR/91.7),0.42)*...
    fuCip     <- fixed(0.65);  label("Ciprofloxacin unbound fraction in plasma")                             # Predictions of clinical drug effect: "Adopted unbound fraction in plasma was 0.65 for CIP"; Supplementary Methods ref 3 (Sadiq 2017)
    rsec      <- fixed(0.674); label("Ciprofloxacin tubular secretion factor")                               # Supplementary Methods Eq. S2: R_SEC set to 0.674
    kpCipKid  <- fixed(8.09);  label("Ciprofloxacin kidney-to-plasma distribution coefficient")              # Supplementary Methods: "Kp_CIPkidney ... was set to the reported value of 8.09"

    # ---- Colistimethate / colistin population PK (Kristoffersson 2020) --
    vcCms    <- fixed(1.52);     label("Central volume of CMS1 and CMS2 (L)")                                 # PKPD_CIPCST_run422b.cpp $PARAM Vc : 1.52 : Kristoffersson et al. CMI 2020
    vpCms    <- fixed(13.0);     label("Peripheral volume of CMS1 and CMS2 (L)")                              # PKPD_CIPCST_run422b.cpp $PARAM Vp : 13.0 : Kristoffersson et al. CMI 2020
    vCol     <- fixed(81.2);     label("Volume of distribution of formed colistin (L)")                       # PKPD_CIPCST_run422b.cpp $PARAM Vcol : 81.2 : Kristoffersson et al. CMI 2020
    clNrCms  <- fixed(5.49);     label("Non-renal clearance of CMS1 and CMS2, i.e. conversion to colistin (L/h)") # PKPD_CIPCST_run422b.cpp $PARAM CLnrcms : 5.49; Supplementary Methods: CL_NR,CMS = 5.49 L/h
    clCol    <- fixed(3.03);     label("Non-renal clearance of formed colistin (L/h)")                        # PKPD_CIPCST_run422b.cpp $PARAM CLcol : 3.03; Supplementary Methods: CL_NR,CST = 3.03 L/h
    siCrcl   <- fixed(0.338);    label("Proportionality factor for renal clearance of CMS1 and CMS2")         # PKPD_CIPCST_run422b.cpp $PARAM SIcrcl : 0.338 : Kristoffersson et al. CMI 2020
    q1Cms    <- fixed(604);      label("Intercompartmental clearance of CMS1 (L/h)")                          # PKPD_CIPCST_run422b.cpp $PARAM Q1 : 604 : Kristoffersson et al. CMI 2020
    q2Cms    <- fixed(7.29);     label("Intercompartmental clearance of CMS2 (L/h)")                          # PKPD_CIPCST_run422b.cpp $PARAM Q2 : 7.29 : Kristoffersson et al. CMI 2020
    sCrcl    <- fixed(0.182);    label("Slope of creatinine clearance on formed-colistin clearance")          # PKPD_CIPCST_run422b.cpp $PARAM SCrCL : 0.182, used as CLcol*(1 + SCrCL*(CLCR-80)*60/1000)
    kpCmsKid <- fixed(12.9);     label("CMS kidney-to-plasma distribution coefficient")                       # Supplementary Methods: "Kp_CMSkidney and Kp_CSTkidney ... the adopted values were 12.9 and 20.8"
    kpCstKid <- fixed(20.8);     label("Colistin kidney-to-plasma distribution coefficient")                  # Supplementary Methods: "Kp_CMSkidney and Kp_CSTkidney ... the adopted values were 12.9 and 20.8"
    funbCst  <- fixed(0.34);     label("Colistin unbound fraction in plasma")                                 # Predictions of clinical drug effect: "0.34 for CST"; Supplementary Methods: f_u,p equal to 0.34
    convCol  <- fixed(1163000);  label("Conversion factor from mol/L to mg/L for colistin (mg per mol)")      # Supplementary Methods: "using the conversion factor of 1 mol/L equivalent to 1163 g/L"
    ivfKid   <- fixed(0.196);    label("Interstitial volume fraction of the kidney")                          # Supplementary Methods: "dividing interstitial volume fraction (IVf = 0.196 for kidney)"

    # ---- Physiology (Supplementary Methods) ------------------------------
    coScale   <- fixed(15);       label("Cardiac-output allometric coefficient (L/h per kg^exponent)")        # PKPD_CIPCST_run422b.cpp $MAIN: CO = 15*pow(WT, 0.74)
    coExp     <- fixed(0.74);     label("Cardiac-output allometric exponent")                                 # PKPD_CIPCST_run422b.cpp $MAIN: CO = 15*pow(WT, 0.74)
    fracQkid  <- fixed(0.19);     label("Fraction of cardiac output perfusing the kidney")                    # PKPD_CIPCST_run422b.cpp $MAIN: QKID = CO*0.19. At WT 80 kg this gives 73 L/h, the Q_kidney of the Supplementary Methods
    vkidPerKg <- fixed(0.004247); label("Kidney volume per kg body weight (L/kg)")                            # PKPD_CIPCST_run422b.cpp $MAIN: VKID = 0.004247*WT. At WT 80 kg this gives 0.34 L, the V_kidney of the Supplementary Methods
    ratioEP   <- fixed(0.5);      label("Ratio of albumin in interstitial fluid to plasma")                   # Supplementary Methods Eq. S3: "E/P ... was set to 0.5 as reported for rat kidney"

    # ---- Bacterial growth model (Table 1) --------------------------------
    kgs  <- fixed(1.69);  label("Growth rate constant of the ciprofloxacin-susceptible subpopulations 1 and 2 (1/h)") # Table 1: k_g1,2 = 1.69 FIX
    kgr  <- fixed(0.316); label("Growth rate constant of the ciprofloxacin-resistant subpopulations 3 and 4 (1/h)")   # Table 1: k_g3,4 = 0.316 FIX
    kd   <- fixed(0.179); label("Natural death rate constant (1/h)")                                                   # Table 1: k_d = 0.179 FIX
    lbmax <- fixed(log(1.56e9)); label("Log of the maximum bacterial density in the system (log CFU/mL)")             # Table 1: Bmax = 1.56e9 CFU/mL (11% RSE); PKPD_CIPCST_run422b.cpp $PARAM Bmax : 1.56*10^(9). Carried on the log scale per the library naming convention and back-transformed in model()

    # ---- Ciprofloxacin PD (Table 1; fixed to Khan 2015) ------------------
    emaxCip    <- fixed(6.75);    label("Emax of ciprofloxacin on k_CIP,1,2 and k_CIP,3,4 (1/h)")            # Table 1: E_maxCIP = 6.75 FIX
    ec50CipSlp <- fixed(1.38);    label("Multiplier of the ciprofloxacin EC50-MIC power relation (mg/L)")     # Table 1: EC50_CIP,1,2 = 1.38 x MIC^0.996 FIX; cpp: EC50CI = 1.38*pow(CIPMIC,0.996)
    ec50CipPwr <- fixed(0.996);   label("Exponent of the ciprofloxacin EC50-MIC power relation")              # Table 1: EC50_CIP,1,2 = 1.38 x MIC^0.996 FIX
    ec50CipR   <- fixed(1.53);    label("EC50 of ciprofloxacin on the ciprofloxacin-resistant subpopulations 3 and 4 (mg/L)") # Table 1: EC50_CIP,3,4 = 1.53 FIX
    hillCip    <- fixed(1.60);    label("Hill factor of ciprofloxacin")                                       # Table 1: gamma_CIP = 1.60 FIX
    fmutCip    <- fixed(1.23e-6); label("Fraction of ciprofloxacin-resistant bacteria in the inoculum")       # Table 1: f_MUTCIP = 1.23 per 10^6 FIX
    fmutPrC47  <- fixed(3.91e-3); label("Fraction of pre-existing resting bacteria in the inoculum, strain C47") # Table 1: f_MUTP = 3910 per 10^6 FIX for C47, 0 otherwise; cpp MUTR : 3910*10^(-6)
    ksncMax    <- fixed(4.69);    label("Emax of the S-to-Nc transfer rate constant (1/h)")                   # Table 1: k_SNc,max = 4.69 FIX
    ksnc50     <- fixed(0.183);   label("EC50 of the S-to-Nc transfer rate constant (tr50)")                  # Table 1: t_r50 = 0.183 FIX
    hillSnc    <- fixed(20);      label("Hill factor of the S-to-Nc transfer rate constant")                  # Table 1: gamma_SNc = 20 FIX
    kncs       <- fixed(0.449);   label("Scaling factor of the Nc-to-S transfer rate constant (1/h)")         # Table 1: k_NcS = 0.449 FIX
    tnc        <- fixed(3.2);     label("Time after which transfer from S to Nc is switched off (h)")         # Table 1: t_Nc = 3.2 FIX; cpp $MAIN: if (self.time<=3.2) FLAG = 1

    # ---- Colistin PD (Table 1) -------------------------------------------
    emaxCstS      <- fixed(50);     label("Emax of colistin on the colistin-susceptible subpopulations 1 and 3 (1/h)")  # Table 1: E_maxCST,1,3 = 50 FIX
    emaxCstRC47   <- fixed(1.44);   label("Emax of colistin on subpopulations 2 and 4, strain C47 (1/h)")               # Table 1: E_maxCST,2,4 = 1.44 for C47; driver Rmd: emaxcst2 = 1.44 for STRAIN 'C47'
    emaxCstRLab   <- fixed(2.44);   label("Emax of colistin on subpopulations 2 and 4, laboratory strains (1/h)")       # Table 1: E_maxCST,2,4 = 2.44; driver Rmd: emaxcst2 = 2.44 for STRAIN 'WT' and 'LM'
    ec50Cst       <- fixed(0.110);  label("EC50 of colistin at the source-strain MIC, before MIC rescaling (mg/L)")     # Table 1: EC50_CST = 0.110; driver Rmd passes EC50 = 0.11*factor with factor = MIC_sim/MIC_mod
    emaxCstPr     <- fixed(2.44);   label("Emax of colistin on the pre-existing resting subpopulation, strain C47 (1/h)") # Table 1: E_maxCST,PR = 2.44; driver Rmd: emaxcst3 = 2.44 for STRAIN 'C47', 0 otherwise
    ec50CstPr     <- fixed(0.2504); label("EC50 of colistin on the pre-existing resting subpopulation (mg/L); NOT MIC-rescaled") # cpp $PARAM EC50CST3 : 0.110*(1+1.28) = 0.2504, equal to Table 1 EC50_CST,PR = 0.251. The driver Rmd never overrides EC50CST3, so unlike ec50Cst it is not rescaled to the simulated MIC; see the vignette Errata
    hillCstSC47   <- fixed(5.26);   label("Hill factor of colistin on subpopulations 1 and 3, strain C47")              # Table 1: gamma_CST,1,3 = 5.26 for C47; driver Rmd: gamd = 5.26
    hillCstS347   <- fixed(3.55);   label("Hill factor of colistin on subpopulations 1 and 3, strain LM347")            # Table 1: gamma_CST,1,3 = 3.55 for LM347; driver Rmd: gamd = 3.55 for STRAIN 'WT'
    hillCstSLm    <- fixed(20);     label("Hill factor of colistin on subpopulations 1 and 3, strains LM378/LM421")     # Table 1: gamma_CST,1,3 = 20 FIX; driver Rmd: gamd = 20 for STRAIN 'LM'
    hillCstRLab   <- fixed(0.270);  label("Hill factor of colistin on subpopulations 2 and 4, laboratory strains")      # Table 1: gamma_CST,2,4 = 0.270; driver Rmd: gamd2 = 0.27 for STRAIN 'WT' and 'LM'
    fmutCst       <- fixed(3.80e-6); label("Fraction of colistin-resistant bacteria in the inoculum")                   # Table 1: f_MUTCST = 3.80 per 10^6; cpp MUTC : 3.8*10^(-6)
    micCstC47     <- fixed(0.75);   label("Colistin MIC of the source strain C47, the denominator of the EC50 rescaling (mg/L)") # Methods: MIC_CST 0.75 mg/L for C47; driver Rmd: CSTMIC = 0.75 with CSTMICfactor = 2/CSTMIC
    micCstLab     <- fixed(0.5);    label("Colistin MIC of the source laboratory strains, the denominator of the EC50 rescaling (mg/L)") # Methods: MIC_CST 0.5 mg/L for LM347/LM378/LM421; driver Rmd: CSTMIC = 0.5

    # ---- Combination model (Table 1) -------------------------------------
    intCstC47 <- fixed(0.541); label("Factor on the colistin EC50 in combination, strains C47/LM378/LM421") # Table 1: I_CST = 0.541; driver Rmd: ICST = 0.541 for STRAIN 'C47' and 'LM'
    intCst347 <- fixed(1.19);  label("Factor on the colistin EC50 in combination, strain LM347")            # Table 1: I_CST = 1.19 for LM347; driver Rmd: ICST = 1.19 for STRAIN 'WT'
    intCip    <- fixed(1.60);  label("Factor on the ciprofloxacin EC50 in combination, shared across strains") # Table 1: I_CIP = 1.60

    # ---- PK variability (mrgsolve $OMEGA of the deposited model) ---------
    # A diagonal matrix; the labels below are the authors' own. mrgsolve draws
    # the entries labelled IOV* once per simulated subject, so they act as
    # additional inter-individual variability here; see the vignette Assumptions.
    etavcCip   ~ 0.3342555 # PKPD_CIPCST_run422b.cpp $OMEGA label IIVVcCIP = 0.3342555
    etaqCip    ~ 0.4176571 # PKPD_CIPCST_run422b.cpp $OMEGA label IIVQCIP = 0.4176571
    etaclCip   ~ 0.2790508 # IIVCLCIP (0.1696584) + IOVCLCIP (0.1093924). The deposited model uses these only as the sum exp(IIVCLCIP + IOVCLCIP), so the two independent normals collapse exactly into one with the summed variance
    etaclNrCms ~ 0.0843750 # IIVCLcms (0.0166767) + IOVCLcms (0.0676983). Likewise only ever used as the sum exp(IIVCLcms + IOVCLcms), in BOTH the non-renal and the renal CMS clearance, so the pairing between those two clearances is preserved
    etaclCol   ~ 0.0591144 # PKPD_CIPCST_run422b.cpp $OMEGA label IIVCLcol = 0.0591144 (SCALE = 1)
    etavCol    ~ 0.0541416 # PKPD_CIPCST_run422b.cpp $OMEGA label IOVfm = 0.0541416. NOT collapsed: the deposited model applies this same draw to the formed-colistin volume AND to the formed-colistin clearance, so it must stay a shared eta
    etavpCms   ~ 0.0566257 # PKPD_CIPCST_run422b.cpp $OMEGA label IOVVp = 0.0566257

    # ---- Residual error ---------------------------------------------------
    addSd <- fixed(1.2983); label("Additive residual SD on the log10 bacterial count (log10 CFU/mL)") # run422b_clean.mod $ERROR: W = SQRT(SIGMA(1) + SIGMA(2)) = sqrt(1.66489 + 0.0205917). The deposited mrgsolve model simulates PREDCFU without residual error; the term is carried here so the model is usable for stochastic simulation
  })

  model({
    eps <- 1e-30
    bmax <- exp(lbmax)

    # ---- 1. Strain-specific parameter selection (BACT) ------------------
    if (BACT == 47) {
      emaxCstR <- emaxCstRC47
      hillCstS <- hillCstSC47
      hillCstR <- hillCstSC47
      hillCstPr <- hillCstSC47
      emaxCstPrStrain <- emaxCstPr
      fmutPr <- fmutPrC47
      intCst <- intCstC47
      micCstSrc <- micCstC47
    } else if (BACT == 347) {
      emaxCstR <- emaxCstRLab
      hillCstS <- hillCstS347
      hillCstR <- hillCstRLab
      hillCstPr <- hillCstS347
      emaxCstPrStrain <- 0
      fmutPr <- 0
      intCst <- intCst347
      micCstSrc <- micCstLab
    } else {
      emaxCstR <- emaxCstRLab
      hillCstS <- hillCstSLm
      hillCstR <- hillCstRLab
      hillCstPr <- hillCstSLm
      emaxCstPrStrain <- 0
      fmutPr <- 0
      intCst <- intCstC47
      micCstSrc <- micCstLab
    }

    # ---- 2. EC50 values for the simulated strain ------------------------
    # Eq. 18 with the power fixed to 1: EC50_sim = EC50_mod x MIC_sim/MIC_mod.
    ec50CstSim <- ec50Cst * micCst / micCstSrc
    ec50Cip <- ec50CipSlp * micCip^ec50CipPwr

    # ---- 3. Individual PK parameters ------------------------------------
    vcCipI <- vcCip * exp(etavcCip)
    qCipI <- qCip * exp(etaqCip)
    clCipI <- clCip * (CRCL / crclRef)^crclExp * exp(etaclCip)
    # Ciprofloxacin renal clearance (Supplementary Eq. S2).
    clrCip <- (CRCL * 60 / 1000) * fuCip * (1 + rsec)

    crclCms <- CRCL * 60 / 1000
    clNrCmsI <- clNrCms * exp(etaclNrCms)
    clColI <- clCol * (1 + (CRCL - 80) * 60 / 1000 * sCrcl) * exp(etaclCol + etavCol)
    clrCmsI <- siCrcl * exp(etaclNrCms) * crclCms
    vColI <- vCol * exp(etavCol)
    vpCmsI <- vpCms * exp(etavpCms)

    # Micro-rate constants, transcribed from $MAIN of the deposited model.
    # Note k1p2p uses the CENTRAL volume, as the authors wrote it.
    k1c1p <- q1Cms / vcCms
    k1p1c <- q1Cms / vpCms
    k1c2c <- clNrCmsI / vcCms
    k1p2p <- clNrCmsI / vcCms
    k2ccol <- clNrCmsI / vcCms
    k2c2p <- q2Cms / vcCms
    k2p2c <- q2Cms / vpCmsI
    kcol0 <- clColI / vColI
    k1c0 <- clrCmsI / vcCms
    k2c0 <- clrCmsI / vcCms

    # Population (non-individualised) clearances used by the kidney model.
    clrCst <- clCol * sCrcl * (CRCL - 80) * 60 / 1000
    clNrCst <- clCol
    clrCms <- siCrcl * crclCms

    # ---- 4. Physiology --------------------------------------------------
    co <- coScale * WT^coExp
    qKid <- co * fracQkid
    vKid <- vkidPerKg * WT

    # ---- 5. Concentrations ----------------------------------------------
    # Unbound fraction in interstitial fluid (Supplementary Eq. S3).
    futCip <- 1 / (1 + ratioEP * (1 - fuCip) / fuCip)
    futCst <- 1 / (1 + ratioEP * (1 - funbCst) / funbCst)

    cipPlasma <- cipCentral / vcCipI * fuCip
    cipKidneyConc <- cipKidney / vKid * futCip

    cmsPlasma <- cms1Central / vcCms + cms2Central / vcCms
    cstPlasmaMolar <- cstCentral / vColI
    cstPlasma <- cstPlasmaMolar * convCol * funbCst
    cstKidneyConc <- cstKidney / vKid * convCol * futCst / ivfKid

    # ---- 6. Concentrations driving the antibacterial effect -------------
    # THIS BLOCK is the only difference from the sibling kidney model.
    cipConc <- cipPlasma
    cstConc <- cstPlasma

    # ---- 7. Interaction factors (Eq. 17) --------------------------------
    if (cstConc * cipConc > 0) {
      iCst <- intCst
      iCip <- intCip
    } else {
      iCst <- 1
      iCip <- 1
    }

    # ---- 8. Drug killing (Eq. 12 and Eq. 16) ----------------------------
    if (cipConc > 0) {
      kCipS <- emaxCip * cipConc^hillCip / ((ec50Cip * iCip)^hillCip + cipConc^hillCip)
      kCipR <- emaxCip * cipConc^hillCip / (ec50CipR^hillCip + cipConc^hillCip)
      ksncS <- ksncMax * (cipConc / ec50Cip)^hillSnc /
        ((cipConc / ec50Cip)^hillSnc + ksnc50^hillSnc)
      ksncR <- ksncMax * (cipConc / ec50CipR)^hillSnc /
        ((cipConc / ec50CipR)^hillSnc + ksnc50^hillSnc)
      kncsS <- kncs * ec50Cip / cipConc
      kncsR <- kncs * ec50CipR / cipConc
    } else {
      kCipS <- 0
      kCipR <- 0
      ksncS <- 0
      ksncR <- 0
      kncsS <- 0
      kncsR <- 0
    }

    if (cstConc > 0) {
      kCstS <- emaxCstS * cstConc^hillCstS / ((ec50CstSim * iCst)^hillCstS + cstConc^hillCstS)
      kCstR <- emaxCstR * cstConc^hillCstR / (ec50CstSim^hillCstR + cstConc^hillCstR)
      kCstPr <- emaxCstPrStrain * cstConc^hillCstPr /
        ((ec50CstPr * iCst)^hillCstPr + cstConc^hillCstPr)
    } else {
      kCstS <- 0
      kCstR <- 0
      kCstPr <- 0
    }

    kTot1 <- kCipS + kCstS
    kTot2 <- kCipS + kCstR
    kTot3 <- kCipR + kCstS
    kTot4 <- kCipR + kCstR
    # The pre-existing resting subpopulation regains ciprofloxacin
    # susceptibility only when colistin is also present.
    if (cstConc * cipConc > 0) {
      kTotPr <- kCipS + kCstPr
    } else {
      kTotPr <- kCstPr
    }

    if (t <= tnc) {
      flagNc <- 1
    } else {
      flagNc <- 0
    }

    # ---- 9. Density-dependent S-to-R transfer (Eqs. 9 and 10) -----------
    btot <- bact_s1 + bact_r1 + bact_nc1 +
      bact_s2 + bact_r2 + bact_nc2 +
      bact_s3 + bact_r3 + bact_nc3 +
      bact_s4 + bact_r4 + bact_nc4 + bact_pr
    ksrS <- btot / bmax * (kgs - kd)
    ksrR <- btot / bmax * (kgr - kd)

    # ---- 10. Growth cut-off below one bacterium in 2 mL -----------------
    vol <- 2
    if (bact_s1 * vol < 0.5) {
      grow1 <- 0
    } else {
      grow1 <- 1
    }
    if (bact_s2 * vol < 0.5) {
      grow2 <- 0
    } else {
      grow2 <- 1
    }
    if (bact_s3 * vol < 0.5) {
      grow3 <- 0
    } else {
      grow3 <- 1
    }
    if (bact_s4 * vol < 0.5) {
      grow4 <- 0
    } else {
      grow4 <- 1
    }

    # ---- 11. Ciprofloxacin PK ODEs (Khachman 2011 plus Supplementary Eq. S1) ----
    d/dt(cipCentral) <- -clCipI / vcCipI * cipCentral -
      qCipI / vcCipI * cipCentral + qCipI / vpCip * cipPeripheral
    d/dt(cipPeripheral) <- qCipI / vcCipI * cipCentral - qCipI / vpCip * cipPeripheral
    d/dt(cipKidney) <- qKid * cipCentral / vcCipI -
      qKid * cipKidney / vKid / kpCipKid - cipCentral / vcCipI * clrCip

    # ---- 12. Colistimethate and colistin PK ODEs -------------------------
    #         (Kristoffersson 2020 plus Supplementary Eqs. S4-S5)
    d/dt(cms1Central) <- k1p1c * cms1Peripheral - k1c1p * cms1Central -
      k1c2c * cms1Central - k1c0 * cms1Central
    d/dt(cms2Central) <- k2p2c * cms2Peripheral - k2c2p * cms2Central -
      k2ccol * cms2Central - k2c0 * cms2Central + k1c2c * cms1Central
    d/dt(cstCentral) <- k2ccol * cms2Central - kcol0 * cstCentral
    d/dt(cms1Peripheral) <- k1c1p * cms1Central - k1p1c * cms1Peripheral - k1p2p * cms1Peripheral
    d/dt(cms2Peripheral) <- k2c2p * cms2Central - k2p2c * cms2Peripheral + k1p2p * cms1Peripheral
    d/dt(cmsKidney) <- qKid * cmsPlasma - cmsKidney / kpCmsKid * clNrCms / WT -
      qKid * cmsKidney / vKid / kpCmsKid - cmsPlasma * clrCms
    d/dt(cstKidney) <- qKid * cstPlasmaMolar + cmsKidney / kpCmsKid * clNrCms / WT -
      qKid * cstKidney / vKid / kpCstKid - cstKidney * clNrCst / WT / kpCstKid -
      cstPlasmaMolar * clrCst

    # ---- 13. Bacterial ODEs (Eqs. 6-8 and 11) ---------------------------
    d/dt(bact_s1) <- grow1 * kgs * bact_s1 - (kd + kTot1) * bact_s1 -
      ksrS * bact_s1 + kncsS * bact_nc1 - ksncS * bact_s1 * flagNc
    d/dt(bact_r1) <- -kd * bact_r1 + ksrS * bact_s1
    d/dt(bact_nc1) <- ksncS * bact_s1 * flagNc - kncsS * bact_nc1 - (kd + kTot1) * bact_nc1

    d/dt(bact_s2) <- grow2 * kgs * bact_s2 - (kd + kTot2) * bact_s2 -
      ksrS * bact_s2 + kncsS * bact_nc2 - ksncS * bact_s2 * flagNc
    d/dt(bact_r2) <- -kd * bact_r2 + ksrS * bact_s2
    d/dt(bact_nc2) <- ksncS * bact_s2 * flagNc - kncsS * bact_nc2 - (kd + kTot2) * bact_nc2

    d/dt(bact_s3) <- grow3 * kgr * bact_s3 - (kd + kTot3) * bact_s3 -
      ksrR * bact_s3 + kncsR * bact_nc3 - ksncR * bact_s3 * flagNc
    d/dt(bact_r3) <- -kd * bact_r3 + ksrR * bact_s3
    d/dt(bact_nc3) <- ksncR * bact_s3 * flagNc - kncsR * bact_nc3 - (kd + kTot3) * bact_nc3

    d/dt(bact_s4) <- grow4 * kgr * bact_s4 - (kd + kTot4) * bact_s4 -
      ksrR * bact_s4 + kncsR * bact_nc4 - ksncR * bact_s4 * flagNc
    d/dt(bact_r4) <- -kd * bact_r4 + ksrR * bact_s4
    d/dt(bact_nc4) <- ksncR * bact_s4 * flagNc - kncsR * bact_nc4 - (kd + kTot4) * bact_nc4

    # The deposited mrgsolve model splits the pre-existing resting bacteria
    # across the 2x2 resistance grid as PR1-PR4. All four share one rate
    # equation and one killing rate, and only their sum is ever used, so they
    # are collapsed here into the single Persister compartment that the fitted
    # NONMEM model (run422b_clean.mod, A(15)) actually used. The two forms are
    # algebraically identical.
    d/dt(bact_pr) <- -(kd + kTotPr) * bact_pr

    # ---- 14. Initial conditions (Eqs. 1-5) ------------------------------
    bact_s1(0) <- inoc * (1 - fmutCip) * (1 - fmutCst) * (1 - fmutPr)
    bact_s2(0) <- inoc * (1 - fmutCip) * fmutCst * (1 - fmutPr)
    bact_s3(0) <- inoc * fmutCip * (1 - fmutCst) * (1 - fmutPr)
    bact_s4(0) <- inoc * fmutCip * fmutCst * (1 - fmutPr)
    bact_pr(0) <- inoc * fmutPr

    # ---- 15. Outputs ----------------------------------------------------
    # Nc bacteria are viable but do not form colonies and are excluded from the
    # plate count, exactly as in $ERROR of the fitted model.
    cfuObserved <- bact_s1 + bact_r1 + bact_s2 + bact_r2 +
      bact_s3 + bact_r3 + bact_s4 + bact_r4 + bact_pr
    Cc <- log10(max(cfuObserved, eps))
    Cc ~ add(addSd)
  })
}
