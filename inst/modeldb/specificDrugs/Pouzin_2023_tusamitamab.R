Pouzin_2023_tusamitamab <- function() {
  description <- "In vitro-in vivo extrapolation PK (IVIVE-PK) model of hepatic CYP3A4, CYP1A2 and CYP2B6 downregulation after tusamitamab ravtansine (SAR408701, anti-CEACAM5 IgG1-SPDB-DM4 ADC) in adults with advanced solid tumors (Pouzin 2023): indirect-turnover CYP abundance whose synthesis is reduced by unbound DM4 + Me-DM4 (inverse Hill, in vitro hepatocyte constants) and, for CYP3A4 only, mechanism-based inactivation by unbound DM4, driven by the integrated multi-analyte popPK model of Pouzin 2022 (DAR1-DAR8 ADC chains, naked antibody, DM4 and Me-DM4)."
  reference <- paste(
    "Pouzin C, Teutonico D, Fagniez N, Ziti-Ljajic S, Perreard-Dumaine A,",
    "Pardon M, Klieber S, Nguyen L. Prediction of CYP Down Regulation after",
    "Tusamitamab Ravtansine Administration (a DM4-Conjugate), Based on an In",
    "Vitro-In Vivo Extrapolation Approach. Clin Pharmacol Ther.",
    "2024;115(2):278-287. doi:10.1002/cpt.3102.",
    "PK layer: Pouzin C, Gibiansky L, Fagniez N, Chadjaa M, Tod M, Nguyen L.",
    "Integrated multiple analytes and semi-mechanistic population",
    "pharmacokinetic model of tusamitamab ravtansine, a DM4 anti-CEACAM5",
    "antibody-drug conjugate. J Pharmacokinet Pharmacodyn. 2022;49(3):381-394.",
    "doi:10.1007/s10928-021-09799-0.",
    sep = " "
  )
  vignette <- "Pouzin_2023_tusamitamab"
  units <- list(
    time = "day",
    dosing = "umol of antibody (= mg of ADC / 150)",
    concentration = "uM (umol/L) in antibody-equivalent for ADC and NAB; uM for free DM4 and MeDM4; CYP abundance as a fraction of baseline (unitless)"
  )

  # What each ODE state holds. The PK states are those of
  # Pouzin_2022_tusamitamab (verified = FALSE there as here); the three hepatic
  # CYP pools are this paper's IVIVE layer.
  compartmentData <- list(
    dar8_central = list(
      analyte = "tusamitamab ravtansine DAR8",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    dar8_peripheral1 = list(
      analyte = "tusamitamab ravtansine DAR8",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    dar7_central = list(
      analyte = "tusamitamab ravtansine DAR7",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    dar7_peripheral1 = list(
      analyte = "tusamitamab ravtansine DAR7",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    dar6_central = list(
      analyte = "tusamitamab ravtansine DAR6",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    dar6_peripheral1 = list(
      analyte = "tusamitamab ravtansine DAR6",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    dar5_central = list(
      analyte = "tusamitamab ravtansine DAR5",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    dar5_peripheral1 = list(
      analyte = "tusamitamab ravtansine DAR5",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    dar4_central = list(
      analyte = "tusamitamab ravtansine DAR4",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    dar4_peripheral1 = list(
      analyte = "tusamitamab ravtansine DAR4",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    dar3_central = list(
      analyte = "tusamitamab ravtansine DAR3",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    dar3_peripheral1 = list(
      analyte = "tusamitamab ravtansine DAR3",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    dar2_central = list(
      analyte = "tusamitamab ravtansine DAR2",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    dar2_peripheral1 = list(
      analyte = "tusamitamab ravtansine DAR2",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    dar1_central = list(
      analyte = "tusamitamab ravtansine DAR1",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    dar1_peripheral1 = list(
      analyte = "tusamitamab ravtansine DAR1",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    central_nab = list(analyte = "naked antibody chain", units = NA_character_, specimen = "plasma", verified = FALSE),
    peripheral1_nab = list(
      analyte = "naked antibody chain",
      units = NA_character_,
      specimen = "plasma",
      verified = FALSE
    ),
    central_dm4 = list(analyte = "DM4 catabolite", units = NA_character_, specimen = "plasma", verified = FALSE),
    central_medm4 = list(analyte = "MeDM4 catabolite", units = NA_character_, specimen = "plasma", verified = FALSE),
    enzyme_3a4 = list(
      analyte = "cytochrome P450 3A4 (hepatic)",
      units = "fraction of baseline abundance",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_1a2 = list(
      analyte = "cytochrome P450 1A2 (hepatic)",
      units = "fraction of baseline abundance",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_2b6 = list(
      analyte = "cytochrome P450 2B6 (hepatic)",
      units = "fraction of baseline abundance",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 254,
    n_studies = 1,
    age_range = "adults (not tabulated in either paper)",
    weight_range = "not tabulated in either paper",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state = "Adults with advanced solid tumors expressing CEACAM5 (PK layer: TED13751 first-in-human study, NCT02187848).",
    dose_range = "Simulated regimens: 100 mg/m^2 IV Q2W (up to 6 administrations) and 170 mg/m^2 IV Q3W (up to 4 administrations), 1.5-h infusions (Pouzin 2023 Methods, 'IVIVE-PK simulations').",
    regions = "Not specified.",
    in_vitro_systems = "CYP downregulation: fresh human hepatocytes from 3 donors (2 evaluable for CYP1A2) treated 72 h with Me-DM4 4-600 nM; mechanism-based inhibition: pooled human liver microsomes (50 donors) (Pouzin 2023 Methods; Supplementary Table S1).",
    external_validation = "9 patients (8 evaluable) receiving SAR566658 (same SPDB-DM4 construct, anti-CA6) 150 or 190 mg/m^2 Q3W with 2 mg oral midazolam on day -1 and day 4; observed vs predicted geometric-mean midazolam AUC ratio 1.17 vs 1.15 (Supplementary Material 3).",
    notes = "The IVIVE layer was not fitted to clinical data: its constants are in vitro measurements (hepatocyte downregulation, HLM inactivation, unbound fractions) and Simcyp V22 library CYP turnover values (Pouzin 2023 Table 2). No inter-individual variability was estimated on the CYP layer; the paper simulated typical profiles only. The DM4 and Me-DM4 drivers come from the Pouzin 2022 integrated popPK model (254 patients), carried here unchanged from Pouzin_2022_tusamitamab including its IIV and residual error.",
    molecular_weights = list(
      note = "Pouzin 2022 Table 2; dose in umol of antibody is mg of ADC / 150.",
      ADC = 150000, # g/mol
      NAB = 144522, # g/mol
      DM4 = 780, # g/mol
      MeDM4 = 794 # g/mol
    )
  )

  ini({
    # PK layer: the Pouzin 2022 integrated popPK model, which Pouzin 2023
    # (Methods, 'Tusamitamab ravtansine, DM4, and Me-DM4 population PK model')
    # used to generate the typical DM4 and Me-DM4 profiles driving the IVIVE
    # layer. Values, IIV and residual error are carried unchanged from
    # Pouzin_2022_tusamitamab; the 'Pouzin 2022' citations below refer to that
    # paper.

    # ----- Distribution (shared across all DAR0..DAR8 ADC species) -----
    lvc       <- log(3.37);  label("Central volume of distribution Vc, shared across DAR0-DAR8 (L)")            # Pouzin 2022 Table 4: Vc = 3.37 L
    lvp       <- log(2.54);  label("Peripheral volume of distribution Vp, shared across DAR0-DAR8 (L)")          # Pouzin 2022 Table 4: Vp = 2.54 L
    lq        <- log(0.543); label("Inter-compartmental clearance Q, shared across DAR0-DAR8 (L/day)")           # Pouzin 2022 Table 4: Q = 0.543 L/day

    # ----- Proteolytic / catabolic clearances -----
    lcladc    <- log(0.392); label("Proteolytic clearance of conjugated ADC species (DAR1-DAR8) CL_ADC (L/day)") # Pouzin 2022 Table 4: CL_ADC = 0.392 L/day
    lclnab    <- log(0.408); label("Proteolytic clearance of naked antibody (DAR0/NAB) CL_NAB (L/day)")          # Pouzin 2022 Table 4: CL_NAB = 0.408 L/day

    # ----- DAR_n -> DAR_(n-1) deconjugation rates (1/day) -----
    # All deconjugation occurs in the central compartment as a first-order
    # process. Pouzin 2022 estimates kdec1..kdec6 individually and sets
    # kdec7 = kdec8 = kdec6 (DAR7/DAR8 plasma levels too low to identify).
    lkdec1    <- log(0.0565); label("Deconjugation rate DAR1 -> DAR0 / NAB, kdec1 (1/day)")                       # Pouzin 2022 Table 4: kdec1 = 0.0565/day
    lkdec2    <- log(0.181);  label("Deconjugation rate DAR2 -> DAR1, kdec2 (1/day)")                             # Pouzin 2022 Table 4: kdec2 = 0.181/day
    lkdec3    <- log(0.340);  label("Deconjugation rate DAR3 -> DAR2, kdec3 (1/day)")                             # Pouzin 2022 Table 4: kdec3 = 0.340/day
    lkdec4    <- log(0.525);  label("Deconjugation rate DAR4 -> DAR3, kdec4 (1/day)")                             # Pouzin 2022 Table 4: kdec4 = 0.525/day
    lkdec5    <- log(0.751);  label("Deconjugation rate DAR5 -> DAR4, kdec5 (1/day)")                             # Pouzin 2022 Table 4: kdec5 = 0.751/day
    lkdec6    <- log(0.938);  label("Deconjugation rate DAR6 -> DAR5, kdec6 (1/day); also reused for kdec7, kdec8") # Pouzin 2022 Table 4: kdec6 = 0.938/day; kdec7 = kdec8 = kdec6 per Methods (DAR7/DAR8 not separately identifiable)

    # ----- DM4 and MeDM4 disposition -----
    # Pouzin 2022: V_DM4 and V_MeDM4 are each fixed at 1 L because formation-
    # limited kinetics prevented simultaneous identification of V and FR_MeDM4.
    # CL_DM4 and CL_MeDM4 are therefore APPARENT clearances (CL/V = kel) with
    # V fixed to 1 L.
    lcldm4    <- log(240);    label("Apparent DM4 clearance CL_DM4 (L/day; V_DM4 fixed to 1 L)")                 # Pouzin 2022 Table 4: CL_DM4 = 240 L/day
    lclmedm4  <- log(0.256);  label("Apparent MeDM4 clearance CL_MeDM4 (L/day; V_MeDM4 fixed to 1 L)")           # Pouzin 2022 Table 4: CL_MeDM4 = 0.256 L/day
    lvdm4     <- fixed(log(1));   label("Apparent DM4 volume V_DM4 (L)")                                  # Pouzin 2022 Methods: V_DM4 fixed to 1 L
    lvmedm4   <- fixed(log(1));   label("Apparent MeDM4 volume V_MeDM4 (L)")                              # Pouzin 2022 Methods: V_MeDM4 fixed to 1 L
    lfrmedm4  <- log(0.0107); label("Apparent fraction of DM4 elimination producing MeDM4, FR_MeDM4 (unitless)") # Pouzin 2022 Table 4: FR_MeDM4 = 0.0107

    # ----- Administered DAR fractions (Pouzin 2022 Tables 3 & 4) -----
    # The Monolix supplement parameterises the dose composition as ratios FR_i
    # over an implicit FR_DAR8 = 1 reference, with SUM = 1 + sum(FR_i) and
    # F_DARi = FR_i / SUM, F_DAR8 = 1/SUM. Typical values below are derived
    # from the published F_DARi (%) by dividing by F_DAR8 = 0.9%; this
    # exactly reproduces Table 4 typical fractions through the SUM transform.
    lfrdar1   <- fixed(log(0.085 / 0.009)); label("log FR_DAR1 ratio, typical FR_DAR1 = F_DAR1/F_DAR8 (unitless)")        # Pouzin 2022 Table 4: F_DAR1 = 8.5% (fixed)
    lfrdar2   <- fixed(log(0.175 / 0.009)); label("log FR_DAR2 ratio, typical FR_DAR2 = F_DAR2/F_DAR8 (unitless)")        # Pouzin 2022 Table 4: F_DAR2 = 17.5% (fixed)
    lfrdar3   <- fixed(log(0.218 / 0.009)); label("log FR_DAR3 ratio, typical FR_DAR3 = F_DAR3/F_DAR8 (unitless)")        # Pouzin 2022 Table 4: F_DAR3 = 21.8% (fixed)
    lfrdar4   <- fixed(log(0.199 / 0.009)); label("log FR_DAR4 ratio, typical FR_DAR4 = F_DAR4/F_DAR8 (unitless)")        # Pouzin 2022 Table 4: F_DAR4 = 19.9% (fixed)
    lfrdar5   <- fixed(log(0.142 / 0.009)); label("log FR_DAR5 ratio, typical FR_DAR5 = F_DAR5/F_DAR8 (unitless)")        # Pouzin 2022 Table 4: F_DAR5 = 14.2% (fixed)
    lfrdar6   <- fixed(log(0.071 / 0.009)); label("log FR_DAR6 ratio, typical FR_DAR6 = F_DAR6/F_DAR8 (unitless)")        # Pouzin 2022 Table 4: F_DAR6 = 7.1% (fixed)
    lfrdar7   <- fixed(log(0.028 / 0.009)); label("log FR_DAR7 ratio, typical FR_DAR7 = F_DAR7/F_DAR8 (unitless)")        # Pouzin 2022 Table 4: F_DAR7 = 2.8% (fixed)
    lfrnab    <- log(0.071 / 0.009);        label("log FR_NAB ratio, typical FR_NAB = F_NAB/F_DAR8; estimated (unitless)") # Pouzin 2022 Table 4: F_NAB = 7.1% (estimated)

    # ----- IVIVE layer: hepatic CYP turnover (Pouzin 2023 Table 2) -----
    # Every value below is an in vitro measurement or a Simcyp library default
    # used as a fixed input; none was estimated against clinical data. Rates
    # are printed per hour and converted to the model's day time unit (x 24).
    kdeg_3a4 <- fixed(0.0193 * 24); label("Hepatic CYP3A4 degradation rate constant (1/day)") # Table 2 kdeg CYP3A4 = 0.0193 h-1 (Simcyp V22 library default)
    kdeg_1a2 <- fixed(0.0183 * 24); label("Hepatic CYP1A2 degradation rate constant (1/day)") # Table 2 kdeg CYP1A2 = 0.0183 h-1 (Simcyp V22 library default)
    kdeg_2b6 <- fixed(0.0217 * 24); label("Hepatic CYP2B6 degradation rate constant (1/day)") # Table 2 kdeg CYP2B6 = 0.0217 h-1 (Simcyp V22 library default)

    cyp0_3a4 <- fixed(137);  label("Baseline hepatic CYP3A4 abundance CYP0 (pmol/mg microsomal protein)") # Table 2 CYP0 CYP3A4 = 137 pmol/mg protein
    cyp0_1a2 <- fixed(52.0); label("Baseline hepatic CYP1A2 abundance CYP0 (pmol/mg microsomal protein)") # Table 2 CYP0 CYP1A2 = 52.0 pmol/mg protein
    cyp0_2b6 <- fixed(21.6); label("Baseline hepatic CYP2B6 abundance CYP0 (pmol/mg microsomal protein)") # Table 2 CYP0 CYP2B6 = 21.6 pmol/mg protein

    fu    <- fixed(0.01);   label("Fraction unbound of DM4 and Me-DM4 in plasma (unitless)")              # Table 2 Fu = 0.01 (in vitro binding assay)
    fuhep <- fixed(0.0909); label("Fraction unbound of DM4 and Me-DM4 in the hepatocyte incubation (unitless)") # Table 2 FuHEP = 0.0909; Table 1 condition 3 (C3 = 1.00 mg/mL)
    fumic <- fixed(0.286);  label("Fraction unbound of DM4 in the microsomal MBI incubation (unitless)")  # Table 2 FuMIC = 0.286; Table 1 condition 2 (C2 = 0.250 mg/mL)

    # Downregulation by DM4 + Me-DM4 (Me-DM4 hepatocyte constants; DM4 assumed equipotent)
    rmax_3a4 <- fixed(0.0859); label("Remaining CYP3A4 synthesis fold at infinite perpetrator concentration Rmax (unitless)") # Table 2 Rmax CYP3A4 = 0.0859 (Table S1 mean of 3 donors)
    rmax_1a2 <- fixed(0.0194); label("Remaining CYP1A2 synthesis fold at infinite perpetrator concentration Rmax (unitless)") # Table 2 Rmax CYP1A2 = 0.0194 (Table S1 mean of 2 donors)
    rmax_2b6 <- fixed(0.212);  label("Remaining CYP2B6 synthesis fold at infinite perpetrator concentration Rmax (unitless)") # Table 2 Rmax CYP2B6 = 0.212 (Table S1 mean of 3 donors)
    ec50_3a4 <- fixed(0.0118); label("Me-DM4 / DM4 concentration for half-maximal CYP3A4 downregulation, nominal in vitro (uM)") # Table 2 EC50 CYP3A4 = 0.0118 uM
    ec50_1a2 <- fixed(0.0159); label("Me-DM4 / DM4 concentration for half-maximal CYP1A2 downregulation, nominal in vitro (uM)") # Table 2 EC50 CYP1A2 = 0.0159 uM
    ec50_2b6 <- fixed(0.0169); label("Me-DM4 / DM4 concentration for half-maximal CYP2B6 downregulation, nominal in vitro (uM)") # Table 2 EC50 CYP2B6 = 0.0169 uM

    # Mechanism-based inhibition of CYP3A4 by DM4 (HLM)
    kinact_3a4   <- fixed(6.12 * 24); label("Maximal CYP3A4 inactivation rate constant by DM4 kinact (1/day)") # Table 2 kinact = 6.12 h-1
    ki_inact_3a4 <- fixed(2.3);       label("DM4 concentration supporting half-maximal CYP3A4 inactivation KI, nominal in vitro (uM)") # Table 2 KI = 2.3 uM

    # ----- Inter-individual variability (omega = SD on log scale; Pouzin 2022 Table 4) -----
    # Variance entered as omega^2; omega is the SD reported in the "standard
    # deviation of the random effect" column of Table 4. For fixed CVs the
    # corresponding etas are wrapped in fixed(); FR_DAR8 carries no eta; its
    # variability is implicit through the FR/SUM renormalisation.
    etalvc      ~ 0.060025                  # 0.245^2  (Vc CV 24.5%)
    etalvp      ~ 0.366025                  # 0.605^2  (Vp CV 60.5%)
    etalq       ~ 0.279841                  # 0.529^2  (Q  CV 52.9%)
    etalcladc   ~ 0.219961                  # 0.469^2  (CL_ADC  CV 46.9%)
    etalclnab   ~ 0.119025                  # 0.345^2  (CL_NAB  CV 34.5%)
    etalkdec1   ~ 0.040804                  # 0.202^2  shared across all kdec_i (per Pouzin 2022 Methods)
    etalcldm4   ~ 0.133225                  # 0.365^2  (CL_DM4   CV 36.5%)
    etalclmedm4 ~ 0.427716                  # 0.654^2  (CL_MeDM4 CV 65.4%)
    etalfrmedm4 ~ 0.522729                  # 0.723^2  (FR_MeDM4 CV 72.3%)
    etalfrdar1  ~ fixed(0.008281)           # 0.091^2  (FR_DAR1; CV from Table 3 administered batch)
    etalfrdar2  ~ fixed(0.008281)           # 0.091^2  (FR_DAR2)
    etalfrdar3  ~ fixed(0.004356)           # 0.066^2  (FR_DAR3)
    etalfrdar4  ~ fixed(0.003481)           # 0.059^2  (FR_DAR4)
    etalfrdar5  ~ fixed(0.010404)           # 0.102^2  (FR_DAR5)
    etalfrdar6  ~ fixed(0.017424)           # 0.132^2  (FR_DAR6)
    etalfrdar7  ~ fixed(0.073984)           # 0.272^2  (FR_DAR7)
    etalfrnab   ~ 0.174724                  # 0.418^2  (FR_NAB; estimated, not fixed)

    # ----- Residual error (Pouzin 2022 Table 4) -----
    # ADC: combined error model. The published additive SD of 1.03 ug/mL is
    # converted to uM (umol/L) by dividing by the SAR408701 MW = 150,000 g/mol
    # and multiplying by 1000 mL/L: 1.03 / 150,000 * 1000 = 0.006867 uM.
    addSd        <- 0.006867; label("Additive residual error on ADC concentration (uM)")                          # Pouzin 2022 Table 4: a_ADC = 1.03 ug/mL (= 0.006867 uM)
    propSd       <- 0.089;    label("Proportional residual error on ADC concentration (fraction)")                # Pouzin 2022 Table 4: b_ADC   = 8.9%
    propSd_nab     <- 0.260;    label("Proportional residual error on NAB concentration (fraction)")                # Pouzin 2022 Table 4: b_NAB   = 26.0%
    propSd_dm4     <- 0.335;    label("Proportional residual error on DM4 concentration (fraction)")                # Pouzin 2022 Table 4: b_DM4   = 33.5%
    propSd_medm4   <- 0.500;    label("Proportional residual error on MeDM4 concentration (fraction)")              # Pouzin 2022 Table 4: b_MeDM4 = 50.0%
    addSd_DARavg   <- 0.219;    label("Additive residual error SD on average DAR, drug-to-antibody ratio (unitless)")                         # Pouzin 2022 Table 4: a_DARaverage = 0.219
  })

  model({
    # ----- Individual structural parameters -----
    vc       <- exp(lvc + etalvc)
    vp       <- exp(lvp + etalvp)
    q        <- exp(lq  + etalq)
    cladc    <- exp(lcladc   + etalcladc)
    clnab    <- exp(lclnab   + etalclnab)
    cldm4    <- exp(lcldm4   + etalcldm4)
    clmedm4  <- exp(lclmedm4 + etalclmedm4)
    vdm4     <- exp(lvdm4)
    vmedm4   <- exp(lvmedm4)
    frmedm4  <- exp(lfrmedm4 + etalfrmedm4)

    # Shared deconjugation eta (etalkdec1) applied to every kdec_i; same
    # IIV across the chain per Pouzin 2022 Methods. kdec7 and kdec8 are
    # constrained equal to kdec6.
    kdec1 <- exp(lkdec1 + etalkdec1)
    kdec2 <- exp(lkdec2 + etalkdec1)
    kdec3 <- exp(lkdec3 + etalkdec1)
    kdec4 <- exp(lkdec4 + etalkdec1)
    kdec5 <- exp(lkdec5 + etalkdec1)
    kdec6 <- exp(lkdec6 + etalkdec1)
    kdec7 <- kdec6
    kdec8 <- kdec6

    # ----- Administered DAR fractions: FR/SUM transform (Pouzin 2022 supplement) -----
    frdar1 <- exp(lfrdar1 + etalfrdar1)
    frdar2 <- exp(lfrdar2 + etalfrdar2)
    frdar3 <- exp(lfrdar3 + etalfrdar3)
    frdar4 <- exp(lfrdar4 + etalfrdar4)
    frdar5 <- exp(lfrdar5 + etalfrdar5)
    frdar6 <- exp(lfrdar6 + etalfrdar6)
    frdar7 <- exp(lfrdar7 + etalfrdar7)
    frnab  <- exp(lfrnab  + etalfrnab)
    sumfr  <- 1 + frdar1 + frdar2 + frdar3 + frdar4 + frdar5 + frdar6 + frdar7 + frnab
    fdar8  <- 1      / sumfr
    fdar7  <- frdar7 / sumfr
    fdar6  <- frdar6 / sumfr
    fdar5  <- frdar5 / sumfr
    fdar4  <- frdar4 / sumfr
    fdar3  <- frdar3 / sumfr
    fdar2  <- frdar2 / sumfr
    fdar1  <- frdar1 / sumfr
    fnab   <- frnab  / sumfr

    # ----- Micro-constants (shared distribution kinetics) -----
    keladc   <- cladc   / vc
    kelnab   <- clnab   / vc
    keldm4   <- cldm4   / vdm4
    kelmedm4 <- clmedm4 / vmedm4
    k12      <- q / vc
    k21      <- q / vp

    # ----- ODE system: 9 parallel two-compartment chains + DM4 + MeDM4 -----
    # All ADC and NAB amounts are tracked in umol of antibody equivalent.
    # DAR_n -> DAR_(n-1) is an irreversible first-order central process.

    # DAR8 chain: receives no upstream input (top of chain).
    d/dt(dar8_central)     <- -keladc * dar8_central - kdec8 * dar8_central -
                               k12 * dar8_central + k21 * dar8_peripheral1
    d/dt(dar8_peripheral1) <-  k12 * dar8_central - k21 * dar8_peripheral1

    # DAR7..DAR1 chains: each receives kdec_(n+1) * upstream_central and loses
    # via proteolysis (keladc) and deconjugation (kdec_n) plus 2-compartment
    # exchange.
    d/dt(dar7_central)     <-  kdec8 * dar8_central - keladc * dar7_central -
                               kdec7 * dar7_central -
                               k12 * dar7_central + k21 * dar7_peripheral1
    d/dt(dar7_peripheral1) <-  k12 * dar7_central - k21 * dar7_peripheral1

    d/dt(dar6_central)     <-  kdec7 * dar7_central - keladc * dar6_central -
                               kdec6 * dar6_central -
                               k12 * dar6_central + k21 * dar6_peripheral1
    d/dt(dar6_peripheral1) <-  k12 * dar6_central - k21 * dar6_peripheral1

    d/dt(dar5_central)     <-  kdec6 * dar6_central - keladc * dar5_central -
                               kdec5 * dar5_central -
                               k12 * dar5_central + k21 * dar5_peripheral1
    d/dt(dar5_peripheral1) <-  k12 * dar5_central - k21 * dar5_peripheral1

    d/dt(dar4_central)     <-  kdec5 * dar5_central - keladc * dar4_central -
                               kdec4 * dar4_central -
                               k12 * dar4_central + k21 * dar4_peripheral1
    d/dt(dar4_peripheral1) <-  k12 * dar4_central - k21 * dar4_peripheral1

    d/dt(dar3_central)     <-  kdec4 * dar4_central - keladc * dar3_central -
                               kdec3 * dar3_central -
                               k12 * dar3_central + k21 * dar3_peripheral1
    d/dt(dar3_peripheral1) <-  k12 * dar3_central - k21 * dar3_peripheral1

    d/dt(dar2_central)     <-  kdec3 * dar3_central - keladc * dar2_central -
                               kdec2 * dar2_central -
                               k12 * dar2_central + k21 * dar2_peripheral1
    d/dt(dar2_peripheral1) <-  k12 * dar2_central - k21 * dar2_peripheral1

    d/dt(dar1_central)     <-  kdec2 * dar2_central - keladc * dar1_central -
                               kdec1 * dar1_central -
                               k12 * dar1_central + k21 * dar1_peripheral1
    d/dt(dar1_peripheral1) <-  k12 * dar1_central - k21 * dar1_peripheral1

    # NAB (DAR0) chain: receives from DAR1 deconjugation; NAB clearance only,
    # no further deconjugation.
    d/dt(central_nab)     <-  kdec1 * dar1_central - kelnab * central_nab -
                              k12 * central_nab + k21 * peripheral1_nab
    d/dt(peripheral1_nab) <-  k12 * central_nab - k21 * peripheral1_nab

    # DM4 catabolite: every DAR>=1 deconjugation event releases exactly one
    # DM4 molecule (1:1 molar). Pouzin 2022 supplement Online Resource 1
    # input_1 expression.
    dm4_input <- kdec1 * dar1_central + kdec2 * dar2_central +
                 kdec3 * dar3_central + kdec4 * dar4_central +
                 kdec5 * dar5_central + kdec6 * dar6_central +
                 kdec7 * dar7_central + kdec8 * dar8_central
    d/dt(central_dm4)   <- dm4_input - keldm4 * central_dm4

    # MeDM4 catabolite: formed sequentially from DM4 elimination, with
    # apparent fraction frmedm4 routed into MeDM4.
    d/dt(central_medm4) <- frmedm4 * keldm4 * central_dm4 - kelmedm4 * central_medm4

    # ----- Bioavailability: split a single dose across 9 ADC chains -----
    # rxode2 / nlmixr2 do not natively split a single dose row across
    # multiple compartments; the vignette therefore issues 9 dose events at
    # each administration (one per chain). Setting f(<chain>) = f<chain>
    # encodes the per-chain administered fraction so the user enters the
    # TOTAL antibody amount on each row and rxode2 scales by f<chain>.
    f(dar1_central) <- fdar1
    f(dar2_central) <- fdar2
    f(dar3_central) <- fdar3
    f(dar4_central) <- fdar4
    f(dar5_central) <- fdar5
    f(dar6_central) <- fdar6
    f(dar7_central) <- fdar7
    f(dar8_central) <- fdar8
    f(central_nab)  <- fnab

    # ----- Observed concentrations -----
    Cc       <- (dar1_central + dar2_central + dar3_central + dar4_central +
                 dar5_central + dar6_central + dar7_central + dar8_central) / vc  # SAR408701 (DAR>=1) by Gyrolab
    Cc_nab   <- central_nab / vc                                                  # NAB (DAR0) by depletion immunoassay
    Cc_dm4   <- central_dm4   / vdm4                                              # Free DM4 by LC-MS/MS
    Cc_medm4 <- central_medm4 / vmedm4                                            # Free MeDM4 by LC-MS/MS

    # Average DAR: total conjugated DM4 / total antibody (ADC + NAB).
    Cc_tab  <- Cc + Cc_nab
    DARavg  <- (1 * dar1_central + 2 * dar2_central + 3 * dar3_central +
                4 * dar4_central + 5 * dar5_central + 6 * dar6_central +
                7 * dar7_central + 8 * dar8_central) / vc / Cc_tab

    # ----- IVIVE layer: hepatic CYP abundance (Pouzin 2023 Methods, p. 280-281) -----
    # Indirect turnover dCYP/dt = kdeg * CYP0 * fold - kdeg * CYP, written on
    # the abundance relative to baseline, enzyme = CYP / CYP0, so enzyme(0) = 1
    # and CYP0 only rescales the absolute-abundance outputs. The ODE is linear
    # in CYP, so the relative form is exact.
    #
    # Downregulation acts on synthesis through the inverse Hill expression of
    # Machavaram 2013, driven by the SUM of DM4 and Me-DM4 (DM4 assumed
    # equipotent with Me-DM4). Unbound plasma concentration (fu * C) is compared
    # with the unbound in vitro EC50 (ec50 * fuhep):
    #   fold = 1 + (rmax - 1) * fu * C / (ec50 * fuhep + fu * C)
    cu_perp   <- fu * (Cc_dm4 + Cc_medm4)
    fold_3a4  <- 1 + (rmax_3a4 - 1) * cu_perp / (ec50_3a4 * fuhep + cu_perp)
    fold_1a2  <- 1 + (rmax_1a2 - 1) * cu_perp / (ec50_1a2 * fuhep + cu_perp)
    fold_2b6  <- 1 + (rmax_2b6 - 1) * cu_perp / (ec50_2b6 * fuhep + cu_perp)

    # Mechanism-based inactivation of CYP3A4 by unbound DM4 only (Yeo 2010
    # form as printed: linear in concentration, unbound plasma DM4 against the
    # unbound microsomal KI). Me-DM4 is a reversible CYP3A4 inhibitor only and
    # is not included.
    kmbi_3a4  <- kinact_3a4 * fu * Cc_dm4 / (ki_inact_3a4 * fumic)

    d/dt(enzyme_3a4) <- kdeg_3a4 * fold_3a4 - kdeg_3a4 * enzyme_3a4 - kmbi_3a4 * enzyme_3a4
    d/dt(enzyme_1a2) <- kdeg_1a2 * fold_1a2 - kdeg_1a2 * enzyme_1a2
    d/dt(enzyme_2b6) <- kdeg_2b6 * fold_2b6 - kdeg_2b6 * enzyme_2b6
    enzyme_3a4(0) <- 1
    enzyme_1a2(0) <- 1
    enzyme_2b6(0) <- 1

    # Absolute abundance (pmol/mg microsomal protein) and the paper's
    # percentage change from baseline, 100 * (CYP(t) - CYP0) / CYP0.
    abund_3a4  <- cyp0_3a4 * enzyme_3a4
    abund_1a2  <- cyp0_1a2 * enzyme_1a2
    abund_2b6  <- cyp0_2b6 * enzyme_2b6
    pctchg_3a4 <- 100 * (enzyme_3a4 - 1)
    pctchg_1a2 <- 100 * (enzyme_1a2 - 1)
    pctchg_2b6 <- 100 * (enzyme_2b6 - 1)

    Cc       ~ add(addSd) + prop(propSd)
    Cc_nab   ~ prop(propSd_nab)
    Cc_dm4   ~ prop(propSd_dm4)
    Cc_medm4 ~ prop(propSd_medm4)
    DARavg   ~ add(addSd_DARavg)
  })
}
