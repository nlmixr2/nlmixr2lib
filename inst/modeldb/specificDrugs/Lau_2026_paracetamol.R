Lau_2026_paracetamol <- function() {
  description <- "Semi-physiological popPK. Joint parent-and-metabolites population PK model for oral and intravenous paracetamol (acetaminophen, PCM) and its glucuronide (PCM-GLU), sulphate (PCM-SUL), and combined cysteine + mercapturate (PCM-CYS & PCM-MER) metabolites in adults with and without obesity (Lau 2026). Extends the intravenous-only model of van Rongen 2016 by resolving first-pass loss into its anatomical sites: a well-stirred liver whose three parallel intrinsic clearances (glucuronidation, sulphation, CYP2E1 oxidation) set pathway-specific hepatic extraction ratios against a weight-driven hepatic blood flow, and a gut wall whose intrinsic CYP2E1 clearance (a fixed fraction of the hepatic oxidative intrinsic clearance) sets the gut extraction ratio against Q_gut. Gut wall, portal vein, and liver are quasi-steady-state algebraic pools rather than ODE states, so absorbed and recirculated drug is presented to the liver on every pass while systemic drug is presented at hepatic blood flow. Enterohepatic recirculation routes 10 percent of newly formed PCM-GLU into a gallbladder compartment that empties over two fixed 6-minute windows per study, releasing drug to a reuptake depot that is deglucuronidated and reabsorbed as parent PCM. Lean body mass scales parent volume, glucuronidation and oxidation intrinsic clearances, and glucuronide elimination clearance; total body weight scales glucuronide volume and drives cardiac output. Study-specific multipliers on the oxidative intrinsic clearance and on all three metabolite elimination clearances separate the Chen oral cohort from the other two studies."
  reference <- paste(
    "Lau C, van Kesteren C, Smeenk RM, Beex-Oosterhuis MM, Koch BCP,",
    "Chan LN, Lin YS, van Rongen A, Knibbe CAJ, Huitema ADR,",
    "Huisman-Siebinga H (2026).",
    "Semi-physiological population pharmacokinetic modeling of oral and",
    "intravenous paracetamol to quantify presystemic metabolism and",
    "enterohepatic recirculation.",
    "CPT Pharmacometrics Syst Pharmacol 15(1):e70168.",
    "doi:10.1002/psp4.70168.",
    sep = " "
  )
  vignette <- "Lau_2026_paracetamol"
  units <- list(time = "min", dosing = "umol", concentration = "umol/L")

  # `depot2` is the enterohepatic reuptake pool (control stream compartment 7,
  # "REUPTAKE PCM FROM GALLBLADDER"): PCM-GLU released from the gallbladder is
  # deglucuronidated in the gut and re-absorbed as parent PCM at ka2, entering
  # the gut wall on the same footing as a fresh oral dose. It is therefore a
  # second absorption depot in the blessed depot<n> family, not a route-specific
  # depot_<route>.
  compartmentData <- list(
    depot            = list(analyte = "paracetamol", units = "umol", specimen = "administration site", verified = TRUE),
    central          = list(analyte = "paracetamol", units = "umol", specimen = "plasma", verified = TRUE),
    transit1_gluc    = list(analyte = "paracetamol glucuronide", units = "umol", specimen = "administration site", verified = TRUE),
    central_gluc     = list(analyte = "paracetamol glucuronide", units = "umol", specimen = "plasma", verified = TRUE),
    central_sulf     = list(analyte = "paracetamol sulphate", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1_sulf = list(analyte = "paracetamol sulphate", units = "umol", specimen = "plasma", verified = TRUE),
    transit1_cysmer  = list(analyte = "paracetamol cysteine + mercapturate", units = "umol", specimen = "administration site", verified = TRUE),
    central_cysmer   = list(analyte = "paracetamol cysteine + mercapturate", units = "umol", specimen = "plasma", verified = TRUE),
    gallbladder      = list(analyte = "paracetamol glucuronide", units = "umol", specimen = "bile", verified = TRUE),
    depot2           = list(analyte = "paracetamol", units = "umol", specimen = "administration site", verified = TRUE)
  )

  # `etaltgb` is the single ETA(11) that the control stream applies to every
  # gallbladder-emptying MTIME. Its typical value is assembled inside model()
  # from the six study-specific fixed start times rather than from a single
  # `ltgb` ini parameter, so it is declared paper-specific.
  paper_specific_etas <- c("etaltgb")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Two roles. (1) Power-law scaling on the glucuronide volume of distribution with reference TBW = 130.9 kg (the pooled median inherited from van Rongen 2016). (2) Drives cardiac output through the Young 2009 polynomial CO = (9119 - exp(9.164 - 0.0291*WT + 0.000391*WT^2 - 0.00000191*WT^3))/1000 L/min (Table S1), from which hepatic, portal, intestinal, mucosal and villous blood flows are derived; every extraction ratio in the model therefore responds to total body weight. Time-fixed at baseline. Source column 'TBW' / 'WT' renamed to canonical 'WT' on input.",
      source_name        = "WT"
    ),
    LBM = list(
      description        = "Lean body mass at baseline (Janmahasatian et al. 2005 formula).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-law scaling with reference LBW = 65.2 kg (pooled median inherited from van Rongen 2016) on parent V, the glucuronidation and CYP2E1 intrinsic hepatic clearances, and the glucuronide elimination clearance. The LBW exponent on the sulphation intrinsic hepatic clearance is retained in the model but held at 0 (control stream THETA(20) 0 FIX; Results: 'We did not find a statistically significant effect of LBW on ... CL_H,int of PCM to PCM-SUL'). Methods 2.3.1 states LBW was calculated as described by Janmahasatian et al. (reference 20). Source column 'LBW' renamed to canonical 'LBM' on input (same biological quantity, no value transformation).",
      source_name        = "LBW"
    ),
    STUDY_PAPAYA = list(
      description        = "1 = subject enrolled in Study 1 (the PAPAYA study, single 1000 mg oral paracetamol suspension, Albert Schweitzer Hospital / Erasmus MC, the Netherlands); 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0",
      notes              = "Selects the Study 1 gallbladder-emptying schedule (60 and 240 min post-dose) and the Study 1 fold increase in the combined PCM residual error (2.25) and the shared oral fold on the PCM-GLU residual error (0.872). Study 1 carries no structural clearance multiplier - it shares the reference clearances with Study 3. STUDY_PAPAYA = STUDY_CHEN = 0 selects Study 3 (van Rongen intravenous), the reference cohort.",
      source_name        = "FLG3"
    ),
    STUDY_CHEN = list(
      description        = "1 = subject enrolled in Study 2 (the previously published Chen et al. oral-liquid study, single 1500 mg oral paracetamol, United States); 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0",
      notes              = "Selects the Study 2 gallbladder-emptying schedule (480 and 660 min post-dose), the Study 2 multipliers on the oxidative intrinsic hepatic clearance (1.94) and on the glucuronide (0.650), sulphate (0.826) and cysteine + mercapturate (1.59) elimination clearances, the additive-only residual error for PCM-CYS (0.288 umol/L), the Study 2 fold increase in the combined PCM residual error (1.87), and the shared oral fold on the PCM-GLU residual error (0.872). STUDY_PAPAYA = STUDY_CHEN = 0 selects Study 3 (van Rongen intravenous), the reference cohort.",
      source_name        = "FLG1"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 69L,
    n_studies      = 3L,
    age_range      = "18-65 years",
    age_median     = "41-49 years by cohort (Table 1)",
    weight_range   = "53-198 kg (TBW); 36-99 kg (LBW)",
    weight_median  = "130.9 kg TBW / 65.2 kg LBW (covariate reference values inherited from van Rongen 2016)",
    sex_female_pct = 71.0,
    race_ethnicity = NA_character_,
    disease_state  = "53 adults with obesity (BMI 30-77.5 kg/m^2; scheduled for laparoscopic Roux-en-Y gastric bypass or gastric sleeve, sampled before surgery) and 16 adults without obesity (healthy participants or patients undergoing oral and maxillofacial surgery). Exclusions in Study 1: pregnancy, liver disease, paracetamol within 24 h.",
    dose_range     = "Study 1: single 1000 mg oral suspension (n = 30). Study 2: single 1500 mg oral liquid (n = 11). Study 3: single 2000 mg intravenous, with optional standard postoperative 1000 mg QID from 8 h onward (n = 28).",
    regions        = "the Netherlands (Studies 1 and 3); United States (Study 2)",
    notes          = "Pooled from three independent clinical studies (Table 1): Study 1 = PAPAYA (this paper, oral); Study 2 = Chen et al. reference [9] (oral); Study 3 = van Rongen et al. reference [5] (intravenous, already packaged as vanRongen_2016_acetaminophen). 41 participants received oral and 28 intravenous paracetamol. Observation counts: 782 PCM, 784 PCM-GLU, 783 PCM-SUL, 766 PCM-CYS & PCM-MER. Concentrations and doses are molar (umol/L, umol); molecular weights used for the conversion are PCM 151.16, PCM-GLU 327.29, PCM-SUL 231.23, PCM-CYS 270.30 and PCM-MER 312.24 g/mol (Methods 2.5)."
  )

  ini({
    # ------------------------------------------------------------------
    # All final estimates are from Table 2 of Lau 2026 (final
    # semi-physiological population PK model); the model structure and
    # every fixed / FIX flag are from the deposited NONMEM control stream
    # in Supplementary Materials 3 (Data S1), and the drug- and
    # system-specific constants and their equations are from Table S1
    # (Supplementary Materials 1). Time is in minutes, amounts in umol,
    # concentrations in umol/L, clearances and blood flows in L/min.
    # ------------------------------------------------------------------

    # ---- Parent paracetamol disposition and absorption ----------------
    lvc            <- log(53.7)     ; label("Paracetamol central volume at LBW = 65.2 kg (L)")                          # Table 2 V_PCM,LBW 65.2 kg = 53.7 L (RSE 3.2%, 95% CI 50.2-57.0); control stream THETA(1)
    e_lbm_vc       <- 1.00          ; label("Power exponent of LBW on paracetamol central volume")                      # Table 2 exponent N = 1.00 (RSE 14.2%, 95% CI 0.73-1.28); control stream THETA(17)
    lka            <- log(0.0531)   ; label("First-order oral absorption rate constant ka (1/min)")                     # Table 2 ka = 0.0531 1/min (RSE 12.3%, 95% CI 0.0415-0.0667); control stream THETA(2). Table 2 prints the unit as min^-1 and the control stream is minute-based; Results 3.2.1 quotes the same number as 0.0531 h^-1, which is a units typo in the narrative (see vignette Errata)
    lfdepot        <- log(0.745)    ; label("Fraction of the oral dose absorbed Fa (bioavailability on depot)")         # Table 2 Fa = 0.745 (RSE 3.3%, 95% CI 0.699-0.792); control stream THETA(3) applied as F1 on the oral dosing compartment
    lka2           <- log(0.00312)  ; label("Reabsorption rate constant of enterohepatically recirculated PCM ka2 (1/min)") # Table 2 ka2 = 0.00312 1/min (RSE 8.2%, 95% CI 0.0025-0.0035); control stream THETA(24)

    # ---- Intrinsic hepatic clearances (well-stirred liver) ------------
    # CL_H,int of the three metabolic routes. The population estimates
    # follow the expected hierarchy glucuronidation > sulphation >
    # oxidation (Results 3.2.1).
    lclint_gluc          <- log(0.300)  ; label("Intrinsic hepatic clearance to PCM-GLU at LBW = 65.2 kg (L/min)")       # Table 2 CL_H GLU,LBW 65.2 kg = 0.300 L/min (RSE 4.2%, 95% CI 0.278-0.327); control stream THETA(4)
    e_lbm_clint_gluc     <- 0.888       ; label("Power exponent of LBW on the glucuronidation intrinsic hepatic clearance") # Table 2 exponent O = 0.888 (RSE 16.3%, 95% CI 0.61-1.18); control stream THETA(19)
    lclint_sulf          <- log(0.0667) ; label("Intrinsic hepatic clearance to PCM-SUL at LBW = 65.2 kg (L/min)")       # Table 2 CL_H SUL,LBW 65.2 kg = 0.0667 L/min (RSE 4.4%, 95% CI 0.0608-0.0725); control stream THETA(5)
    e_lbm_clint_sulf     <- fixed(0)    ; label("Power exponent of LBW on the sulphation intrinsic hepatic clearance") # Control stream THETA(20) '0 FIX'; Results 3.2.1: 'We did not find a statistically significant effect of LBW on ... CL_H,int of PCM to PCM-SUL, as was previously identified in Study 3'
    lclint_cysmer        <- log(0.0417) ; label("Intrinsic hepatic clearance to PCM-CYS & PCM-MER at LBW = 65.2 kg (L/min)") # Table 2 CL_H CYS,LBW 65.2 kg = 0.0417 L/min (RSE 7.9%, 95% CI 0.0359-0.0486); control stream THETA(6)
    e_lbm_clint_cysmer   <- 1.36        ; label("Power exponent of LBW on the oxidative intrinsic hepatic clearance")    # Table 2 exponent P = 1.36 (RSE 16.9%, 95% CI 0.94-1.82); control stream THETA(21)
    e_study_chen_clint_cysmer <- 1.94   ; label("Study 2 (Chen) multiplier on the oxidative intrinsic hepatic clearance") # Table 2 factor R = 1.94 (RSE 13.0%, 95% CI 1.50-2.48); control stream THETA(29), entering as LTVCHCYS = log(THETA(6)) + FLG1*log(THETA(29))

    # ---- Gut wall (intestinal CYP2E1 metabolism) ----------------------
    # CL_G,int is estimated as a fraction of the oxidative CL_H,int
    # (Table 2 footnote a: "CL_G,int coded as fraction of intrinsic
    # hepatic clearance"); only the oxidative route is given an
    # intestinal arm (control stream: "only via CYP2E1").
    fclint_gut     <- 0.00474         ; label("Intestinal CYP2E1 intrinsic clearance as a fraction of the oxidative intrinsic hepatic clearance") # Table 2 'Fraction CL_G,int' = 0.00474 (RSE 20.7%, 95% CI 0.0032-0.0069); control stream THETA(8)
    fub            <- fixed(0.855)    ; label("Fraction of paracetamol unbound in blood f_UB")                            # Table 2 / Table S1 f_UB = 0.855 fixed, citing Avdeef 2012 (Table S1 reference 1); control stream THETA(7) '0.855 FIX'
    fug            <- fixed(1)        ; label("Fraction of paracetamol unbound in the gut wall f_UG")                     # Table S1 f_UG = 1 fixed, citing Frechen 2013 / Yang 2007 (Table S1 reference 2); control stream 'FUG= 1'
    lclperm        <- fixed(log(3.79)); label("Permeability clearance through the enterocytes CL_perm (L/min)")           # Table 2 / Table S1 CL_perm = 3.79 L/min fixed = P_app * A with P_app 31.6e-6 cm/s and A 200 m^2, citing Faassen 2003 (Table S1 reference 3); control stream THETA(9) '3.79 FIX'

    # ---- Fixed semi-physiological volumes (Table S1, Frechen 2013) ----
    lvh            <- fixed(log(1))   ; label("Liver volume V_H (L)")                                                     # Table S1 V_H = 1 L fixed, citing Frechen 2013; control stream 'VH= 1'
    lvpv           <- fixed(log(1))   ; label("Portal vein volume V_PV (L)")                                              # Table S1 V_PV = 1 L fixed, citing Frechen 2013; control stream 'VPV= 1'
    lvgw           <- fixed(log(1))   ; label("Gut wall volume V_GW (L)")                                                 # Table S1 V_GW = 1 L fixed, citing Frechen 2013; control stream 'VGW= 1'

    # ---- Metabolite disposition ---------------------------------------
    lvc_gluc       <- log(26.8)       ; label("PCM-GLU volume of distribution at TBW = 130.9 kg (L)")                     # Table 2 V_PCM-GLU,TBW 130.9 kg = 26.8 L (RSE 4.3%, 95% CI 24.6-29.1); control stream THETA(10)
    e_wt_vc_gluc   <- 0.485           ; label("Power exponent of TBW on the PCM-GLU volume of distribution")              # Table 2 exponent S = 0.485 (RSE 18.1%, 95% CI 0.301-0.646); control stream THETA(18)
    lvc_sulf       <- fixed(log(5.66)); label("PCM-SUL central volume of distribution (L)")                               # Table 2 V_PCM-SUL,central = V_PCM-SUL,peripheral = 5.66 L fixed on the basis of van Rongen 2016; Table S1 same; control stream THETA(11) '5.66 FIX'
    lvp_sulf       <- fixed(log(5.66)); label("PCM-SUL peripheral volume of distribution (L, equal to central)")          # Table S1 V_central,sulfate = V_peripheral,sulfate; control stream 'V8=V4'
    lq_sulf        <- log(0.695)      ; label("PCM-SUL inter-compartmental clearance Q (L/min)")                          # Table 2 Q = 0.695 L/min (RSE 22.5%, 95% CI 0.472-1.11); control stream THETA(13)
    lvc_cysmer     <- fixed(log(15.6)); label("PCM-CYS & PCM-MER combined volume of distribution (L)")                    # Table 2 V_PCM-CYS/MER = 15.6 L fixed on the basis of van Rongen 2016; Table S1 same; control stream THETA(12) '15.6 FIX'

    # ---- Metabolite elimination clearances ----------------------------
    lcle_gluc               <- log(0.193) ; label("PCM-GLU elimination clearance at LBW = 65.2 kg (L/min)")               # Table 2 CL_E GLU = 0.193 L/min (RSE 3.6%, 95% CI 0.181-0.208); control stream THETA(14)
    e_lbm_cle_gluc          <- 0.537      ; label("Power exponent of LBW on the PCM-GLU elimination clearance")           # Table 2 exponent T = 0.537 (RSE 23.0%, 95% CI 0.287-0.791); control stream THETA(25)
    e_study_chen_cle_gluc   <- 0.650      ; label("Study 2 (Chen) multiplier on the PCM-GLU elimination clearance")       # Table 2 factor U = 0.650 (RSE 8.8%, 95% CI 0.555-0.780); control stream THETA(26)
    lcle_sulf               <- log(0.100) ; label("PCM-SUL elimination clearance (L/min)")                                # Table 2 CL_E,SUL = 0.100 L/min (RSE 2.6%, 95% CI 0.096-0.106); control stream THETA(15)
    e_study_chen_cle_sulf   <- 0.826      ; label("Study 2 (Chen) multiplier on the PCM-SUL elimination clearance")       # Table 2 factor W = 0.826 (RSE 5.1%, 95% CI 0.747-0.912); control stream THETA(27)
    lcle_cysmer             <- log(0.397) ; label("PCM-CYS & PCM-MER elimination clearance (L/min)")                      # Table 2 CL_E CYS = 0.397 L/min (RSE 7.8%, 95% CI 0.344-0.462); control stream THETA(16)
    e_study_chen_cle_cysmer <- 1.59       ; label("Study 2 (Chen) multiplier on the PCM-CYS & PCM-MER elimination clearance") # Table 2 factor X = 1.59 (RSE 14.9%, 95% CI 1.20-2.11); control stream THETA(28)

    # ---- Metabolite formation transit rate constants ------------------
    # One transit compartment delays the appearance of each of the two
    # hepatically formed metabolites that showed a formation lag; the
    # sulphate arm has no transit compartment.
    lktr_cysmer    <- log(0.00332) ; label("PCM-CYS & PCM-MER formation transit rate constant ktr (1/min)")               # Table 2 ktr = 0.00332 1/min (RSE 2.5%, 95% CI 0.0032-0.0035); control stream THETA(30), used as K105 (transit -> central PCM-CYS)
    lktr_gluc      <- log(0.0687)  ; label("PCM-GLU formation transit rate constant ktr2 (1/min)")                        # Table 2 ktr2 = 0.0687 1/min (RSE 5.9%, 95% CI 0.061-0.076); control stream THETA(31); the deposited stream omits the 'K93=KTR2' assignment line but K93 is the compartment-9 -> compartment-3 transit rate and KTR2 is the only rate constant defined for it (see vignette Errata)

    # ---- Enterohepatic recirculation (gallbladder) --------------------
    # Gallbladder emptying start times are fixed per study to the values
    # that were found to fit each fasting state; each study has two
    # emptying events 180 min apart, each lasting 6 min. All six start
    # times carry the same single log-normal IIV (etaltgb).
    fgb            <- fixed(0.1)      ; label("Fraction of newly formed PCM-GLU routed into the gallbladder F_GB")        # Table 2 F_GB = 0.1 fixed; Results 3.2.1: 'The fraction of PCM-GLU routed into the gallbladder was fixed at 0.1, representing 10% of the total PCM-GLU formed'; control stream THETA(23)
    lkgg           <- fixed(log(10))  ; label("Gallbladder-to-reuptake transfer rate constant k_GG during emptying (1/min)") # Table 2 / Table S1 k_GG = 10 1/min fixed; control stream THETA(22)
    ldgb           <- fixed(log(6))   ; label("Duration of each gallbladder emptying window (min)")                       # Table 2 'Duration of gallbladder emptying' = 6 min fixed; Table S1 same; control stream THETA(33)
    ltgb1_papaya    <- fixed(log(60))  ; label("Study 1 (PAPAYA) first gallbladder emptying start time (min)")            # Table 2 'Start gallbladder emptying, Study 1' = 60 fixed; control stream THETA(32)
    ltgb2_papaya    <- fixed(log(240)) ; label("Study 1 (PAPAYA) second gallbladder emptying start time (min)")           # Table 2 'Start gallbladder emptying, Study 1' = 240 fixed; control stream THETA(36)
    ltgb1_chen      <- fixed(log(480)) ; label("Study 2 (Chen) first gallbladder emptying start time (min)")              # Table 2 'Start gallbladder emptying, Study 2' = 480 fixed; control stream THETA(35)
    ltgb2_chen      <- fixed(log(660)) ; label("Study 2 (Chen) second gallbladder emptying start time (min)")             # Table 2 'Start gallbladder emptying, Study 2' = 660 fixed; control stream THETA(38)
    ltgb1_vanrongen <- fixed(log(360)) ; label("Study 3 (van Rongen, intravenous) first gallbladder emptying start time (min)")  # Table 2 'Start gallbladder emptying, Study 3' = 360 fixed; control stream THETA(34)
    ltgb2_vanrongen <- fixed(log(540)) ; label("Study 3 (van Rongen, intravenous) second gallbladder emptying start time (min)") # Table 2 'Start gallbladder emptying, Study 3' = 540 fixed; control stream THETA(37)

    # ---- Interindividual variability ----------------------------------
    # Table 2 reports IIV as %CV of a log-normal (exponential) random
    # effect, so the internal variance is omega^2 = log(1 + CV^2). The
    # control stream confirms the mapping: the fixed 131.1% CV on the
    # gallbladder start time is $OMEGA "1 FIX", and log(1 + 1.311^2) =
    # 1.0002. IIV on the intestinal clearance was $OMEGA "0 FIX"
    # (control stream ETA(7)) and is therefore omitted entirely.
    etalvc            ~ 0.073414   # Table 2 IIV V_PCM        = 27.6% CV (RSE 21.4%) -> log(1 + 0.276^2)
    etalka            ~ 0.592332   # Table 2 IIV ka           = 89.9% CV (RSE 25.6%) -> log(1 + 0.899^2)
    etalfdepot        ~ 0.019963   # Table 2 IIV Fa           = 14.2% CV (RSE 29.0%) -> log(1 + 0.142^2)
    etalclint_gluc    ~ 0.090630   # Table 2 IIV CL_H GLU     = 30.8% CV (RSE 18.8%) -> log(1 + 0.308^2)
    etalclint_sulf    ~ 0.121227   # Table 2 IIV CL_H SUL     = 35.9% CV (RSE 19.7%) -> log(1 + 0.359^2)
    etalclint_cysmer  ~ 0.158904   # Table 2 IIV CL_H CYS     = 41.5% CV (RSE 23.4%) -> log(1 + 0.415^2)
    etalvc_gluc       ~ 0.036945   # Table 2 IIV V_GLU        = 19.4% CV (RSE 27.8%) -> log(1 + 0.194^2)
    etalcle_gluc      ~ 0.052867   # Table 2 IIV CL_E,GLU     = 23.3% CV (RSE 20.4%) -> log(1 + 0.233^2)
    etalcle_cysmer    ~ 0.179932   # Table 2 IIV CL_E,CYS     = 44.4% CV (RSE 23.5%) -> log(1 + 0.444^2)
    etaltgb           ~ fixed(1)   # Table 2 IIV start of gallbladder emptying = 131.1% CV, control stream $OMEGA 11 '1 FIX' (log(1 + 1.311^2) = 1.0002); Results 3.2.1 carries it over from the model of de Winter et al. to avoid estimation problems

    # ---- Residual variability -----------------------------------------
    # Combined proportional + additive error on the linear scale for PCM,
    # PCM-GLU and PCM-SUL (Methods Eq. 2). The control stream $ERROR block
    # multiplies the whole combined error term for PCM and PCM-GLU by a
    # study-specific fold factor. For PCM-CYS & PCM-MER the error is
    # proportional for Studies 1 and 3 (the additive $SIGMA is 0 FIX) and
    # purely additive for Study 2 (Methods 2.3.2).
    propSd        <- 0.139   ; label("Paracetamol proportional residual SD, reference study (fraction)")                  # Table 2 Proportional error PCM = 13.9% CV (RSE 8.1%, 95% CI 13.0-15.5); control stream $SIGMA 1
    addSd         <- 0.00002 ; label("Paracetamol additive residual SD, reference study (umol/L)")                        # Table 2 Additive error PCM = 0.00002 umol/L (RSE 49.0%); control stream $SIGMA 5. The reported 95% CI (0.00005-0.00028) does not bracket the point estimate; the Estimate column is used here (see vignette Errata). Numerically negligible against PCM concentrations of order 10-100 umol/L
    propSd_gluc   <- 0.137   ; label("PCM-GLU proportional residual SD, reference study (fraction)")                      # Table 2 Proportional error PCM-GLU = 13.7% CV (RSE 9.6%, 95% CI 12.6-15.2); control stream $SIGMA 2
    addSd_gluc    <- 6.08    ; label("PCM-GLU additive residual SD, reference study (umol/L)")                            # Table 2 Additive error PCM-GLU = 6.08 umol/L (RSE 22.6%, 95% CI 3.89-9.28); control stream $SIGMA 6
    propSd_sulf   <- 0.175   ; label("PCM-SUL proportional residual SD (fraction)")                                       # Table 2 Proportional error PCM-SUL = 17.5% CV (RSE 10.2%, 95% CI 15.9-19.3); control stream $SIGMA 3
    addSd_sulf    <- 2.99    ; label("PCM-SUL additive residual SD (umol/L)")                                             # Table 2 Additive error PCM-SUL = 2.99 umol/L (RSE 29.5%, 95% CI 1.45-4.89); control stream $SIGMA 7
    propSd_cysmer <- 0.285   ; label("PCM-CYS & PCM-MER proportional residual SD, Studies 1 and 3 (fraction)")            # Table 2 Proportional error PCM-CYS = 28.5% CV (RSE 6.5%, 95% CI 26.8-30.5); control stream $SIGMA 4, applied only when FLG1 = 0
    addSd_cysmer  <- 0.288   ; label("PCM-CYS & PCM-MER additive residual SD, Study 2 only (umol/L)")                     # Table 2 Additive error PCM-CYS Study 2 = 0.288 umol/L (RSE 12.3%, 95% CI 0.234-0.370); control stream $SIGMA 9, applied only when FLG1 = 1. The Studies 1 and 3 additive component is $SIGMA 8 '0 FIX'

    e_study_chen_resSd        <- 1.87  ; label("Study 2 (Chen) fold increase in the combined paracetamol residual error")   # Table 2 'Fold increase in combined residual variability (additive + proportional) for PCM, Study 2' = 1.87 (RSE 7.7%, 95% CI 1.59-2.16); control stream THETA(39)
    e_study_papaya_resSd      <- 2.25  ; label("Study 1 (PAPAYA) fold increase in the combined paracetamol residual error") # Table 2 'Fold increase in combined residual variability (additive + proportional) for PCM, Study 1' = 2.25 (RSE 6.5%, 95% CI 1.97-2.55); control stream THETA(40)
    e_study_oral_resSd_gluc   <- 0.872 ; label("Oral studies (1 and 2) fold change in the combined PCM-GLU residual error") # Table 2 'Fold increase in combined residual variability (additive + proportional) for PCM-GLU, Study 1 and 2' = 0.872 (RSE 6.2%, 95% CI 0.77-0.98); control stream THETA(41)
  })

  model({
    # ------------------------------------------------------------------
    # Study indicators. Study 3 (van Rongen, intravenous) is the
    # reference cohort and is selected when both indicators are 0.
    # ------------------------------------------------------------------
    studyOral <- STUDY_PAPAYA + STUDY_CHEN
    studyIv   <- 1 - studyOral

    # ------------------------------------------------------------------
    # Individual parameters. LBW (canonical LBM) enters as a power-law
    # multiplier with reference 65.2 kg; TBW (canonical WT) enters on the
    # PCM-GLU volume with reference 130.9 kg. Study 2 (Chen) multipliers
    # are applied as X^STUDY_CHEN, reproducing the control stream's
    # log(THETA)*FLG1 additive-on-log-scale form exactly.
    # ------------------------------------------------------------------
    vc           <- exp(lvc + etalvc) * (LBM / 65.2)^e_lbm_vc
    ka           <- exp(lka + etalka)
    fa           <- exp(lfdepot + etalfdepot)
    ka2          <- exp(lka2)

    clint_gluc   <- exp(lclint_gluc   + etalclint_gluc)   * (LBM / 65.2)^e_lbm_clint_gluc
    clint_sulf   <- exp(lclint_sulf   + etalclint_sulf)   * (LBM / 65.2)^e_lbm_clint_sulf
    clint_cysmer <- exp(lclint_cysmer + etalclint_cysmer) * (LBM / 65.2)^e_lbm_clint_cysmer *
                      e_study_chen_clint_cysmer^STUDY_CHEN

    vc_gluc      <- exp(lvc_gluc + etalvc_gluc) * (WT / 130.9)^e_wt_vc_gluc
    vc_sulf      <- exp(lvc_sulf)
    vp_sulf      <- exp(lvp_sulf)
    q_sulf       <- exp(lq_sulf)
    vc_cysmer    <- exp(lvc_cysmer)

    cle_gluc     <- exp(lcle_gluc + etalcle_gluc) * (LBM / 65.2)^e_lbm_cle_gluc *
                      e_study_chen_cle_gluc^STUDY_CHEN
    cle_sulf     <- exp(lcle_sulf) * e_study_chen_cle_sulf^STUDY_CHEN
    cle_cysmer   <- exp(lcle_cysmer + etalcle_cysmer) * e_study_chen_cle_cysmer^STUDY_CHEN

    ktr_gluc     <- exp(lktr_gluc)
    ktr_cysmer   <- exp(lktr_cysmer)

    # ------------------------------------------------------------------
    # System-specific physiology (Table S1). Cardiac output follows the
    # Young 2009 polynomial in total body weight; every downstream blood
    # flow is a fixed fraction of it, so all extraction ratios scale with
    # body size even where no explicit covariate was estimated.
    # ------------------------------------------------------------------
    vh       <- exp(lvh)
    vpv      <- exp(lvpv)
    vgw      <- exp(lvgw)
    clperm   <- exp(lclperm)

    co       <- (9119 - exp(9.164 - 0.0291 * WT + 0.000391 * WT^2 -
                              0.00000191 * WT^3)) / 1000
    qh       <- 0.25 * co
    qpv      <- 0.75 * qh
    qha      <- 0.25 * qh
    qint     <- 0.40 * qh
    qmucosa  <- 0.80 * qint
    qvilli   <- 0.60 * qmucosa
    qgut     <- (qvilli * clperm) / (qvilli + clperm)

    # ------------------------------------------------------------------
    # Extraction ratios. Each hepatic route contributes its own
    # extraction ratio and the fraction escaping hepatic metabolism is
    # 1 minus their sum, exactly as coded in the deposited control stream
    # ("FH=1-EHGLU-EHSUL-EHCYS"). Only the oxidative route has an
    # intestinal arm.
    # ------------------------------------------------------------------
    eh_gluc   <- (clint_gluc   * fub) / (qh + clint_gluc   * fub)
    eh_sulf   <- (clint_sulf   * fub) / (qh + clint_sulf   * fub)
    eh_cysmer <- (clint_cysmer * fub) / (qh + clint_cysmer * fub)
    fh        <- 1 - eh_gluc - eh_sulf - eh_cysmer

    clint_gut <- fclint_gut * clint_cysmer
    egut      <- (clint_gut * fug) / (qgut + clint_gut * fug)
    fg        <- 1 - egut

    # ------------------------------------------------------------------
    # Gallbladder emptying schedule. Two 6-minute windows per study at
    # study-specific fixed start times, both shifted by the same
    # log-normal random effect (control stream ETA(11) on every MTIME).
    # mtime() is rxode2's equivalent of NONMEM's MTIME: it inserts a
    # solver break point so the discontinuous gate is integrated exactly,
    # and the [start, end) window reproduces MPAST(n) - MPAST(n+1).
    # ------------------------------------------------------------------
    dgb  <- exp(ldgb)
    kgg  <- exp(lkgg)
    tgb1 <- (exp(ltgb1_papaya) * STUDY_PAPAYA + exp(ltgb1_chen) * STUDY_CHEN +
               exp(ltgb1_vanrongen) * studyIv) * exp(etaltgb)
    tgb2 <- (exp(ltgb2_papaya) * STUDY_PAPAYA + exp(ltgb2_chen) * STUDY_CHEN +
               exp(ltgb2_vanrongen) * studyIv) * exp(etaltgb)
    mtime(tgbStart1) <- tgb1
    mtime(tgbEnd1)   <- tgb1 + dgb
    mtime(tgbStart2) <- tgb2
    mtime(tgbEnd2)   <- tgb2 + dgb
    gbOpen <- ((t >= tgbStart1) & (t < tgbEnd1)) + ((t >= tgbStart2) & (t < tgbEnd2))

    # ------------------------------------------------------------------
    # Micro-constants for the metabolite compartments.
    # ------------------------------------------------------------------
    kel_gluc   <- cle_gluc   / vc_gluc
    kel_sulf   <- cle_sulf   / vc_sulf
    kel_cysmer <- cle_cysmer / vc_cysmer
    k_sulf_cp  <- q_sulf / vc_sulf
    k_sulf_pc  <- q_sulf / vp_sulf

    # ------------------------------------------------------------------
    # Quasi-steady-state gut wall, portal vein and liver pools. These are
    # algebraic (not ODE) states in the deposited control stream: each
    # pool's amount is set so that its outflow equals its inflow at every
    # instant. Absorbed drug (fresh oral dose at ka, enterohepatically
    # recirculated drug at ka2) enters the gut wall, escapes gut
    # metabolism with probability fg, joins systemic drug carried in the
    # portal vein, and is presented to the liver together with the
    # hepatic-artery flow. `liverin` is the total molar rate presented to
    # the liver (umol/min); fh of it survives into the systemic
    # circulation and each eh_* fraction becomes the corresponding
    # metabolite.
    # ------------------------------------------------------------------
    absrate <- ka * depot + ka2 * depot2
    agutw   <- absrate / (qvilli / vgw)
    apv     <- ((qvilli / vgw) * agutw * fg + (qpv / vc) * central) / (qpv / vpv)
    ah      <- ((qha / vc) * central + (qpv / vpv) * apv) / (qh / vh)
    liverin <- (qh / vh) * ah
    agutwm  <- egut * agutw

    # ------------------------------------------------------------------
    # ODE system - Figure 2 of Lau 2026 and the $DES block of the
    # deposited control stream. Compartment order matches the control
    # stream where it matters for dosing: depot (1) is the oral dosing
    # compartment, central (2) receives intravenous doses.
    # ------------------------------------------------------------------
    d/dt(depot)            <- -ka * depot
    d/dt(central)          <-  fh * liverin - (qha / vc) * central - (qpv / vc) * central
    d/dt(transit1_gluc)    <-  (1 - fgb) * liverin * eh_gluc - ktr_gluc * transit1_gluc
    d/dt(central_gluc)     <-  ktr_gluc * transit1_gluc - kel_gluc * central_gluc
    d/dt(central_sulf)     <-  liverin * eh_sulf - k_sulf_cp * central_sulf +
                               k_sulf_pc * peripheral1_sulf - kel_sulf * central_sulf
    d/dt(peripheral1_sulf) <-  k_sulf_cp * central_sulf - k_sulf_pc * peripheral1_sulf
    d/dt(transit1_cysmer)  <-  liverin * eh_cysmer - ktr_cysmer * transit1_cysmer
    d/dt(central_cysmer)   <-  ktr_cysmer * transit1_cysmer +
                               agutwm * (qvilli / vgw) - kel_cysmer * central_cysmer
    d/dt(gallbladder)      <-  fgb * liverin * eh_gluc - kgg * gallbladder * gbOpen
    d/dt(depot2)           <-  kgg * gallbladder * gbOpen - ka2 * depot2

    # Fraction absorbed applies to the oral dosing compartment (control
    # stream F1 = Fa on compartment 1). Intravenous doses go directly to
    # `central` and are not multiplied by Fa.
    f(depot) <- fa

    # ------------------------------------------------------------------
    # Plasma concentrations. Doses are in umol and volumes in L, so
    # amount / volume is umol/L, the unit the paper modelled in
    # (Methods 2.5). Users dosing in mg convert with the paracetamol
    # molecular weight: dose_umol = dose_mg / 0.15116.
    # ------------------------------------------------------------------
    Cc        <- central        / vc
    Cc_gluc   <- central_gluc   / vc_gluc
    Cc_sulf   <- central_sulf   / vc_sulf
    Cc_cysmer <- central_cysmer / vc_cysmer

    # ------------------------------------------------------------------
    # Residual error, reproducing the $ERROR block. The PCM and PCM-GLU
    # combined errors are multiplied by a study-specific fold factor; the
    # PCM-CYS & PCM-MER error switches between proportional (Studies 1
    # and 3) and additive (Study 2).
    # ------------------------------------------------------------------
    resFoldPcm     <- e_study_chen_resSd * STUDY_CHEN +
                      e_study_papaya_resSd * STUDY_PAPAYA + studyIv
    resFoldGluc    <- e_study_oral_resSd_gluc * studyOral + studyIv
    propSdCc       <- propSd        * resFoldPcm
    addSdCc        <- addSd         * resFoldPcm
    propSdCcGluc   <- propSd_gluc   * resFoldGluc
    addSdCcGluc    <- addSd_gluc    * resFoldGluc
    propSdCcCysmer <- propSd_cysmer * (1 - STUDY_CHEN)
    addSdCcCysmer  <- addSd_cysmer  * STUDY_CHEN

    Cc        ~ prop(propSdCc)       + add(addSdCc)
    Cc_gluc   ~ prop(propSdCcGluc)   + add(addSdCcGluc)
    Cc_sulf   ~ prop(propSd_sulf)    + add(addSd_sulf)
    Cc_cysmer ~ prop(propSdCcCysmer) + add(addSdCcCysmer)
  })
}
