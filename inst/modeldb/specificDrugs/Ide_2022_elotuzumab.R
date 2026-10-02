Ide_2022_elotuzumab <- function() {
  description <- "Two-compartment population PK model for elotuzumab (anti-SLAMF7 humanized IgG1) in multiple myeloma given as monotherapy or with lenalidomide/dexamethasone (Ld) or pomalidomide/dexamethasone (Pd) (Ide 2022); parallel linear and Michaelis-Menten elimination from the central compartment plus second-order target-mediated elimination from the peripheral compartment driven by a non-renewable target pool, with time-varying serum M protein on Vmax and backbone (monotherapy / Ld / Pd) effects on linear CL and KINT."
  reference <- "Ide T, Osawa M, Sanghavi K, Vezina HE. Population pharmacokinetic and exposure-response analyses of elotuzumab plus pomalidomide and dexamethasone for relapsed and refractory multiple myeloma. Cancer Chemother Pharmacol. 2022;89(1):129-140. doi:10.1007/s00280-021-04365-4"
  vignette <- "Ide_2022_elotuzumab"
  paper_specific_etas <- c("etaruv")
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "elotuzumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "elotuzumab", units = "mg", specimen = "tissue", verified = TRUE),
    target = list(analyte = "SLAMF7 (peripheral target)", units = "ug/mL", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling on CL (1.32), VC (0.345), Q (0.75 fixed) and VP (0.696), each as (WT/75)^theta (Supplementary Methods control stream: VWT = WT/75).",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Exponential effect on VC only (VCSEX = 0.797). The source NONMEM column SEXF already follows the canonical 1 = female convention. The female effect on CL is carried in the control stream as THETA(24) = 0 FIXED, i.e. not part of the model.",
      source_name = "SEXF"
    ),
    RACE_ASIAN = list(
      description = "Indicator for Asian race",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Exponential effect on VC only (VCRACE = 0.885). Derived from the source NONMEM column RACEN via ASIAN = as.integer(RACEN == 3). The Asian effect on CL is carried as THETA(25) = 0 FIXED, i.e. not part of the model.",
      source_name = "RACEN"
    ),
    B2M = list(
      description = "Baseline serum beta-2-microglobulin",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only as the thresholded indicator B2M >= 3.5 mg/L on VC (VCB2MICG>0.35 = 1.11). The source column B2MICG_Z is in mg/dL (Table 1: median 0.32 mg/dL) and the control stream thresholds it at B2MICG.GE.0.35 mg/dL, which is 3.5 mg/L in the canonical unit (1 mg/dL = 10 mg/L); supply B2M in mg/L. The lower threshold (B2MICG >= 0.2 mg/dL) is carried on VC and CL as THETA(33)/THETA(35) = 0 FIXED and the >= 0.35 mg/dL effect on CL as THETA(34) = 0 FIXED, so none of those enter the model. Missing values were imputed to the dataset median 0.32 mg/dL (= 3.2 mg/L, below the threshold).",
      source_name = "B2MICG_Z"
    ),
    MCPROT = list(
      description = "Time-varying serum M (monoclonal) protein concentration",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters un-log-transformed as exp(0.27 * MCPROT) on Vmax of the central Michaelis-Menten elimination; the Vmax reference (12.2 ug/mL/day) is MCPROT = 0 g/dL (Supplementary Table S2 abbreviations). Time-varying: the paper linearly interpolated measured M-protein to assign a value at every record time (Methods, Population pharmacokinetic analysis). Source NONMEM column TMCPROT_Z; missing values (-99) imputed to the dataset median 2.05 g/dL.",
      source_name = "TMCPROT_Z"
    ),
    COMBO_LEN_DEX = list(
      description = "Lenalidomide + dexamethasone combination-therapy indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0; the paper's reference patient is on Ld (COMBO_LEN_DEX = 1, COMBO_POM_DEX = 0), at which the Ld/Pd/monotherapy factors on CL and KINT are all 1.",
      notes = "Together with COMBO_POM_DEX this selects the backbone: Ld (1, 0), Pd (0, 1) or elotuzumab monotherapy (0, 0). The source control stream uses LENDEX = as.integer(STUDY != 204011), which is 1 for BOTH the Ld and the Pd studies (so it is really an 'any IMiD backbone' flag), and POMDEX = as.integer(STUDY == 204125). The library columns map back as LENDEX = COMBO_LEN_DEX + COMBO_POM_DEX, POMDEX = COMBO_POM_DEX. Never set both columns to 1.",
      source_name = "LENDEX"
    ),
    COMBO_POM_DEX = list(
      description = "Pomalidomide + dexamethasone combination-therapy indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no Pd; the paper's reference is Ld co-administration)",
      notes = "Multiplicative factors relative to Ld: CLPomDex = 0.811 on linear CL and KINTPomDex = 0.487 on KINT (19% and 51% decreases; Abstract and Results). Only study CA204-125 (ELOQUENT-3) contributed Pd patients; source column POMDEX = as.integer(STUDY == 204125).",
      source_name = "POMDEX"
    ),
    STUDY_PHASE3 = list(
      description = "Phase III study indicator; 1 = ELOQUENT-2 (CA204-004), 0 = the phase 1 and phase 2 studies (CA204-005, CA204-007, CA204-011, CA204-125/ELOQUENT-3)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase 1/2 studies)",
      notes = "Selects the residual-error magnitude only: the phase 1/2 studies have the saturable residual SD multiplied by SDphase1,2 = 0.843. It touches no structural or covariate parameter, so it does not change typical-value predictions. Source control stream: STOTHER = 1 if STUDY != 204004, i.e. STOTHER = 1 - STUDY_PHASE3. For a new simulation, 1 reproduces the ELOQUENT-2 residual magnitude and 0 the (smaller) phase 1/2 magnitude.",
      source_name = "STUDY"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Carried in the Supplementary Methods control stream as LOG(AGE/65)*THETA(26) on CL with THETA(26) = 0 FIXED: a CL covariate of the predecessor model (Ide 2020, J Clin Pharmacol 2021;61:64-73) that this update switches off."
    ),
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Control stream LOG(GFR/100)*THETA(27) on CL with THETA(27) = 0 FIXED; not part of the model."
    ),
    LDH = list(
      description = "Baseline lactate dehydrogenase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Control stream LOG(LDH/200)*THETA(28) on CL with THETA(28) = 0 FIXED; not part of the model."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Control stream LOG(ALB/3.5)*THETA(30) on CL (ALB in g/dL) with THETA(30) = 0 FIXED; not part of the model."
    ),
    HEPIMP = list(
      description = "Hepatic impairment indicator (NCI ODWG mild or worse)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal hepatic function)",
      notes = "Control stream THETA(29)*HEPA on CL with THETA(29) = 0 FIXED; not part of the model."
    ),
    ECOG_GE1 = list(
      description = "ECOG performance status >= 1 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG 0)",
      notes = "Control stream THETA(31)*ECOG1 on CL with THETA(31) = 0 FIXED; not part of the model."
    ),
    ECOG_GE2 = list(
      description = "ECOG performance status >= 2 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG <= 1)",
      notes = "Control stream THETA(32)*ECOG2 on CL with THETA(32) = 0 FIXED; not part of the model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 440L,
    n_observations = 8180L,
    n_studies = 5L,
    age_range = "37-88 years (median 66; mean 65.5, SD 9.8)",
    weight_range = "40-150 kg (median 75; mean 75.6, SD 16.7)",
    sex_female_pct = 41,
    race_ethnicity = c(White = 80, `Black/African American` = 6, Asian = 13, `Other/Pacific Islander` = 2),
    disease_state = "Multiple myeloma: relapsed/refractory (CA204-004 ELOQUENT-2, CA204-005, CA204-007, CA204-125 ELOQUENT-3) and high-risk smoldering myeloma (CA204-011, monotherapy).",
    dose_range = "Elotuzumab 10 or 20 mg/kg IV. With Ld: 10 mg/kg QW for two 28-day cycles then 10 mg/kg Q2W. With Pd (ELOQUENT-3): 10 mg/kg QW for two 28-day cycles then 20 mg/kg Q4W. Monotherapy (CA204-011): 10 or 20 mg/kg.",
    co_medication = "Monotherapy 31 (7%), lenalidomide/dexamethasone 349 (79%), pomalidomide/dexamethasone 60 (14%)",
    renal_function = "eGFR mean 74.1 (SD 23.5) mL/min/1.73 m^2, range 4.58-124; normal 28%, mild 46%, moderate 21%, severe 3%, renal failure 2%",
    hepatic_function = "Normal 91%, mild impairment 9%, moderate <1%",
    mprotein_baseline = "Serum M-protein mean 2.25 (SD 1.58) g/dL, median 2.05 (0-7.7)",
    b2m_baseline = "Beta-2 microglobulin mean 0.43 (SD 0.371) mg/dL, median 0.32 (0.04-3.47)",
    regions = "Global (five BMS-sponsored trials; CA204-005 enrolled Japanese patients)",
    notes = "Ide 2022 Table 1 (population PK dataset, N = 440, 8180 serum concentrations). ECOG 0/1/2: 51/43/7%. Anti-drug antibodies detected at least once in 23%. ELISA LLOQ 190 ng/mL."
  )

  ini({
    # Structural parameters (Supplementary Table S2) at the reference patient:
    # WT = 75 kg, lenalidomide/dexamethasone co-administration, male,
    # non-Asian, B2M < 3.5 mg/L, MCPROT = 0 g/dL for Vmax.
    lcl <- log(0.0834); label("Linear (nonspecific) clearance CL_REF (L/day)") # Table S2: CLREF = exp(theta1) = 0.0834 L/day
    lvc <- log(4.06); label("Central volume of distribution VC_REF (L)") # Table S2: VCREF = exp(theta2) = 4.06 L
    lq <- log(0.512); label("Intercompartmental clearance Q_REF (L/day)") # Table S2: QREF = exp(theta3) = 0.512 L/day
    lvp <- log(1.93); label("Peripheral volume of distribution VP_REF (L)") # Table S2: VPREF = exp(theta4) = 1.93 L
    lrmax <- log(849); label("Initial target concentration in the peripheral compartment RMAX (ug/mL)") # Table S2: RMAX = exp(theta5) = 849 ug/mL
    lkint <- log(0.216e-3); label("Second-order target-mediated elimination rate constant KINT_REF (mL/ug/day)") # Table S2: KINTREF = exp(theta6) = 0.216 x 10^-3 /day/(ug/mL); control stream KINT = EXP(MU_6+ETA(6))/1000
    lvmax <- log(12.2); label("Maximum Michaelis-Menten elimination rate VMAX_REF at MCPROT = 0 (ug/mL/day)") # Table S2: VMAXREF = exp(theta7) = 12.2 ug/mL/day
    lkm <- log(281); label("Michaelis-Menten constant KM (ug/mL)") # Table S2: KM = exp(theta8) = 281 ug/mL

    # Covariate effects (Supplementary Table S2).
    e_wt_cl <- 1.32; label("Power exponent of WT/75 on CL (unitless)") # Table S2: CLWT = theta9 = 1.32
    e_wt_vc <- 0.345; label("Power exponent of WT/75 on VC (unitless)") # Table S2: VCWT = theta10 = 0.345
    e_wt_q <- fixed(0.75); label("Power exponent of WT/75 on Q (unitless)") # Table S2: QWT = theta11 = 0.75 Fixed
    e_wt_vp <- 0.696; label("Power exponent of WT/75 on VP (unitless)") # Table S2: VPWT = theta12 = 0.696
    e_mono_cl <- log(0.825); label("Log of CLMono; CL is multiplied by CLMono^-1 under elotuzumab monotherapy (unitless)") # Table S2: CLMono = exp(theta19) = 0.825, applied as (CLMono)^-Mono
    e_combo_pom_dex_cl <- log(0.811); label("Log factor of Pd (vs Ld) co-administration on CL (unitless)") # Table S2: CLPomDex = exp(theta20) = 0.811
    e_sexf_vc <- log(0.797); label("Log factor of female sex on VC (unitless)") # Table S2: VCSEX = exp(theta17) = 0.797
    e_race_asian_vc <- log(0.885); label("Log factor of Asian race on VC (unitless)") # Table S2: VCRACE = exp(theta18) = 0.885
    e_b2m_ge35_vc <- log(1.11); label("Log factor of B2M >= 3.5 mg/L (0.35 mg/dL) on VC (unitless)") # Table S2: VCB2MICG>0.35 = exp(theta36) = 1.11
    e_mcprot_vmax <- 0.27; label("Coefficient of MCPROT on log(VMAX) (per g/dL)") # Table S2: VMAXMCPROT = theta23 = 0.27 (g/dL)^-1
    e_mono_kint <- log(9.78); label("Log of KINTMono; KINT is multiplied by KINTMono^-1 under elotuzumab monotherapy (unitless)") # Table S2: KINTMono = exp(theta21) = 9.78, applied as (KINTMono)^-Mono
    e_combo_pom_dex_kint <- log(0.487); label("Log factor of Pd (vs Ld) co-administration on KINT (unitless)") # Table S2: KINTPomDex = exp(theta22) = 0.487

    # Inter-individual variability (Supplementary Table S2); all diagonal
    # (control stream $OMEGA has no BLOCK).
    etalcl ~ 0.158 # Table S2: omega2 CL = 0.158 (CV 39.7%)
    etalvc ~ 0.0361 # Table S2: omega2 VC = 0.0361 (CV 19%)
    etalq ~ 0.454 # Table S2: omega2 Q = 0.454 (CV 67.4%)
    etalvp ~ 0.133 # Table S2: omega2 VP = 0.133 (CV 36.5%)
    etalrmax ~ 0.189 # Table S2: omega2 RMAX = 0.189 (CV 43.4%)
    etalkint ~ 1.69 # Table S2: omega2 KINT = 1.69 (CV 130%)
    etalvmax ~ fixed(0.0001) # Table S2: omega2 VMAX = 0.0001, held near zero per footnote d (IMPMAP requirement)
    etalkm ~ 0.385 # Table S2: omega2 KM = 0.385 (CV 62%)
    etaruv ~ 0.183 # Table S2: omega2 epsilon = 0.183 (CV 42.7%), IIV on the residual-error magnitude

    # Residual error: saturable log-scale SD (Supplementary Table S2; control
    # stream $ERROR W = (SDL-(SDL-SDH)*TY/(SD50+TY))*THETA(16)**STOTHER*EXP(ETA(9)),
    # Y = LOG(TY) + W*EPS(1), $SIGMA 1 FIXED).
    sdL <- 2.46; label("Log-scale residual SD at low concentrations SDL (unitless)") # Table S2: SDL = theta13 = 2.46
    sdH <- 0.0976; label("Log-scale residual SD at high concentrations SDH (unitless)") # Table S2: SDH = theta14 = 0.0976
    sd50 <- 6.17; label("Concentration at which the residual SD is (SDL+SDH)/2, SD50 (ug/mL)") # Table S2: SD50 = theta15 = 6.17 ug/mL
    sdPhase12 <- 0.843; label("Multiplier on the residual SD for the phase 1/2 studies (unitless)") # Table S2: SDphase1,2 = theta16 = 0.843
  })

  model({
    # ---- Backbone indicators --------------------------------------------
    # Control stream: LENDEX = 1 for every non-monotherapy study (Ld AND Pd),
    # POMDEX = 1 for the Pd study. mono = 1 - LENDEX.
    mono <- 1 - COMBO_LEN_DEX - COMBO_POM_DEX
    B2M_GE35 <- B2M >= 3.5 # control stream B2MICG.GE.0.35 (mg/dL) = 3.5 mg/L

    # ---- Individual parameters (control stream MU_1..MU_8) ---------------
    cl <- exp(lcl + etalcl + e_wt_cl * log(WT / 75) - e_mono_cl * mono + e_combo_pom_dex_cl * COMBO_POM_DEX)
    vc <- exp(lvc + etalvc + e_wt_vc * log(WT / 75) + e_sexf_vc * SEXF + e_race_asian_vc * RACE_ASIAN + e_b2m_ge35_vc * B2M_GE35)
    q <- exp(lq + etalq + e_wt_q * log(WT / 75))
    vp <- exp(lvp + etalvp + e_wt_vp * log(WT / 75))
    rmax <- exp(lrmax + etalrmax)
    kint <- exp(lkint + etalkint - e_mono_kint * mono + e_combo_pom_dex_kint * COMBO_POM_DEX)
    vmax <- exp(lvmax + etalvmax) * exp(e_mcprot_vmax * MCPROT)
    km <- exp(lkm + etalkm)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    Cc <- central / vc

    # ---- ODEs (control stream $DES) --------------------------------------
    # central, peripheral1 are amounts (mg); target is the peripheral target
    # concentration (ug/mL), hence the division by vp in its equation.
    d/dt(central) <- -k12 * central + k21 * peripheral1 - kel * central - vmax * central / (Cc + km)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1 - kint * peripheral1 * target
    d/dt(target) <- -kint * peripheral1 / vp * target
    target(0) <- rmax # control stream A_0(3) = RMAX

    # ---- Residual error (control stream $ERROR) --------------------------
    # Y = log(Cc) + W * EPS with EPS ~ N(0, 1) is lnorm(W) on Cc.
    W <- (sdL - (sdL - sdH) * Cc / (sd50 + Cc)) * sdPhase12^(1 - STUDY_PHASE3) * exp(etaruv)
    Cc ~ lnorm(W)
  })
}
