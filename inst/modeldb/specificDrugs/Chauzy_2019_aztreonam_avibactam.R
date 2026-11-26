Chauzy_2019_aztreonam_avibactam <- function() {
  description <- "In vitro (four multidrug-resistant, metallo-beta-lactamase-producing Enterobacteriaceae: Escherichia coli 1266865, Citrobacter freundii 974673, Enterobacter cloacae 1285905 and 1318536). Semi-mechanistic PD model of static time-kill experiments with aztreonam and avibactam, alone and in combination, fitted jointly to bacterial counts and measured aztreonam concentrations: a susceptible (P1) and a less-susceptible (P2) bacterial subpopulation with shared logistic growth, sigmoidal Emax kill by aztreonam on both (P2 EC50 raised by a resistance factor) and by avibactam on P1 only, avibactam enhancement of aztreonam as a bi-exponential fall in the aztreonam EC50 with the avibactam concentration, and first-order aztreonam degradation by beta-lactamases that rises exponentially with log10 bacterial density and is inhibited by avibactam (fractional inhibitory Emax). Aztreonam is a dosed state; avibactam is a static covariate. Strain selected by BACT."
  reference <- "Chauzy A, Gaelzer Silva Torres B, Buyck J, de Jonge B, Adier C, Marchand S, Couet W, Gregoire N. Semimechanistic Pharmacodynamic Modeling of Aztreonam-Avibactam Combination to Understand Its Antimicrobial Activity Against Multidrug-Resistant Gram-Negative Bacteria. CPT Pharmacometrics Syst Pharmacol. 2019;8(11):815-824. doi:10.1002/psp4.12452. Structural equations: Supplementary Material S1 (Pharmacodynamic model, and the NONMEM 7.4 model code $PK / $DES / $ERROR). Parameter estimates: main-text Table 2 (identical to the $THETA values of the published code)."
  vignette <- "Chauzy_2019_aztreonam_avibactam"
  units <- list(time = "h", dosing = "mg/L", concentration = "mg/L")

  # Aztreonam is added to the flask at time 0 as a concentration (AMT in mg/L
  # into NONMEM compartment 3, COMP=(DRUG)); `central` therefore holds the
  # aztreonam concentration in the medium directly (implicit unit volume).
  dosing <- c("central")

  compartmentData <- list(
    bact_susceptible_susceptible = list(
      analyte = "P1, the susceptible bacterial subpopulation (killed by aztreonam and by avibactam)",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact_intermediate_resistant = list(
      analyte = "P2, the less-susceptible bacterial subpopulation (aztreonam EC50 raised by the resistance factor; not killed by avibactam)",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    central = list(
      analyte = "aztreonam in the time-kill medium (cation-adjusted Mueller-Hinton broth)",
      units = "mg/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    BACT = list(
      description = "Bacterial strain indicator. Selects the strain-specific growth, kill, enhancement, degradation and residual-error parameters that Table 2 reports per strain, and the strain's initial P2 fraction.",
      units = "(categorical)",
      type = "categorical",
      reference_category = "none -- every level has its own parameter set; there is no reference strain",
      notes = paste(
        "Levels are the STRAIN codes of the authors' NONMEM dataset (Supplementary Material S1 model code, $PK comments):",
        "1 = Escherichia coli 1266865 (NDM-5, TEM-OSBL(b), CMY-42; no aztreonam degradation observed, so no degradation sub-model);",
        "12 = Citrobacter freundii 974673 (NDM-1, SHV-12(2be), TEM-OSBL(2b), CTX-M-3, CMY-34);",
        "15 = Enterobacter cloacae 1285905 (NDM-1, CTX-M-15);",
        "16 = Enterobacter cloacae 1318536 (NDM-1, CTX-M-15).",
        "Any other value yields NA parameters and outputs."
      ),
      source_name = "STRAIN"
    ),
    CONC_AVI_MGL = list(
      description = "Static (time-invariant) avibactam concentration in the time-kill flask. Drives the avibactam sigmoidal Emax kill of P1, the bi-exponential decrease of the aztreonam EC50 (enhancement) on both subpopulations, and the fractional inhibition of aztreonam degradation.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- 0 mg/L is the avibactam-free arm (growth control and aztreonam alone); the time-kill grid spanned 0.004-64 mg/L (Methods, Static time-kill studies)",
      notes = "Held constant over the 27-h experiment: no relevant avibactam degradation was observed (Results), and the NONMEM dataset carries avibactam as the data item AVI rather than as a compartment.",
      source_name = "AVI"
    )
  )

  population <- list(
    species = "in vitro (Escherichia coli 1266865, Citrobacter freundii 974673, Enterobacter cloacae 1285905 and 1318536)",
    n_subjects = 4L,
    n_studies = 1L,
    disease_state = "Four multidrug-resistant Enterobacteriaceae clinical isolates from International Health Management Associates (Schaumburg, IL), each producing several beta-lactamases including the metallo-beta-lactamase NDM-1 or NDM-5 (Table 1). Aztreonam MICs 32-512 mg/L alone and 0.125-4 mg/L with 4 mg/L avibactam; avibactam MICs 16-32 mg/L (Table 1, EUCAST broth microdilution).",
    model_system = "Static (constant-concentration) time-kill experiments over 27 h at 37 degrees C with shaking (150 rpm), 20 mL cation-adjusted Mueller-Hinton broth, inoculum about 1.5 x 10^6 CFU/mL. Per strain 32 conditions: growth control, 6 avibactam-alone levels (0.004-64 mg/L), 5 aztreonam-alone levels (0.25-4 x MIC) and 20 combinations (aztreonam 0.25-4 x the MIC in the presence of 0-8 mg/L avibactam), each at least in duplicate. Bacterial counts at 0, 2, 4, 6, 8, 10, 15, 21, 24 and 27 h by spiral plating (lower limit of quantification 200 CFU/mL); aztreonam measured by LC-MS/MS at 0, 8, 15 and 24 h (every sampling time for aztreonam alone).",
    initial_inoculum = "Estimated 10^6.01 CFU/mL, shared by the four strains (Table 2)",
    dose_range = "Static concentrations: aztreonam 0.25-4 x MIC (MICs 32-512 mg/L alone, 0.125-4 mg/L with avibactam), avibactam 0.004-64 mg/L, alone and in combination",
    notes = paste(
      "Estimated in NONMEM 7.4 (FOCE with the Laplacian method, ADVAN9, TOL=6); bacterial counts below the limit of quantification were handled with the M3 method.",
      "One joint fit to all four strains: the maximum population size and the inoculum are shared, every other parameter is strain-specific.",
      "No random effects: $OMEGA on the growth rate is fixed to 0, so the model is a typical-value model with strain-specific residual errors (additive on log10 CFU/mL, proportional on aztreonam).",
      "Fitted in two steps (bacteria first with no degradation, then aztreonam data with the bacterial parameters fixed), then all parameters re-estimated simultaneously for the final model."
    )
  )

  ini({
    # --- Parameters shared by the four strains (Table 2; Code $THETA 5-6) ---
    log10bmax <- 8.44; label("Maximum bacterial population size supported by the system (log10 CFU/mL)") # Table 2: Nmax 8.44 (RSE 0.27%), all strains; Code THETA(5) LGNMAX
    log10inoc <- 6.01; label("Bacterial count at t = 0 (log10 CFU/mL)")                                  # Table 2: INOC 6.01 (RSE 0.39%), all strains; Code THETA(6) LGINOC

    # --- Fraction of P2 in the inoculum, fixed from plating experiments (Table 2; Code $PK MUT) ---
    fmut_s1  <- fixed(0.00001);   label("Fraction of the less-susceptible subpopulation P2 at t = 0, E. coli 1266865 (fraction)")      # Table 2: P2 fraction 0.00001 fix; Code IF(STRAIN.EQ.1) MUT = 0.00001
    fmut_s12 <- fixed(0.00001);   label("Fraction of the less-susceptible subpopulation P2 at t = 0, C. freundii 974673 (fraction)")   # Table 2: P2 fraction 0.00001 fix; Code IF(STRAIN.EQ.12) MUT = 0.00001
    fmut_s15 <- fixed(0.00002);   label("Fraction of the less-susceptible subpopulation P2 at t = 0, E. cloacae 1285905 (fraction)")   # Table 2: P2 fraction 0.00002 fix; Code IF(STRAIN.EQ.15) MUT = 0.00002
    fmut_s16 <- fixed(0.0000001); label("Fraction of the less-susceptible subpopulation P2 at t = 0, E. cloacae 1318536 (fraction)")   # Table 2: P2 fraction 0.0000001 fix; Code IF(STRAIN.EQ.16) MUT = 0.0000001

    # --- Growth rate constant, shared by P1 and P2 (Table 2 Kgrowth; Code THETA 1-4) ---
    kg_s1  <- 0.907; label("Bacterial growth rate constant, E. coli 1266865 (1/h)")      # Table 2: Kgrowth 0.907 (RSE 5.43%); Code THETA(1)
    kg_s12 <- 1.1;   label("Bacterial growth rate constant, C. freundii 974673 (1/h)")   # Table 2: Kgrowth 1.1 (RSE 3.64%); Code THETA(2)
    kg_s15 <- 1.15;  label("Bacterial growth rate constant, E. cloacae 1285905 (1/h)")   # Table 2: Kgrowth 1.15 (RSE 2.86%); Code THETA(3)
    kg_s16 <- 1.24;  label("Bacterial growth rate constant, E. cloacae 1318536 (1/h)")   # Table 2: Kgrowth 1.24 (RSE 1.23%); Code THETA(4)

    # --- Aztreonam maximum kill rate constant (Table 2 Emax,ATM; Code THETA 7-10) ---
    emax_atm_s1  <- 2.62; label("Maximum kill rate constant of aztreonam, E. coli 1266865 (1/h)")      # Table 2: Emax,ATM 2.62 (RSE 5.1%); Code THETA(7)
    emax_atm_s12 <- 1.76; label("Maximum kill rate constant of aztreonam, C. freundii 974673 (1/h)")   # Table 2: Emax,ATM 1.76 (RSE 3.1%); Code THETA(8)
    emax_atm_s15 <- 1.6;  label("Maximum kill rate constant of aztreonam, E. cloacae 1285905 (1/h)")   # Table 2: Emax,ATM 1.6 (RSE 2.3%); Code THETA(9)
    emax_atm_s16 <- 1.86; label("Maximum kill rate constant of aztreonam, E. cloacae 1318536 (1/h)")   # Table 2: Emax,ATM 1.86 (RSE 1.7%); Code THETA(10)

    # --- Bi-exponential aztreonam EC50: EC50 = A*exp(-alpha*AVI) + B*exp(-beta*AVI) (Table 2; Code THETA 11-26) ---
    ec50a_atm_s1  <- 31.8;     label("First term A of the aztreonam EC50 bi-exponential, E. coli 1266865 (mg/L)")      # Table 2: A 31.8 (RSE 11.5%); Code THETA(11)
    ec50a_atm_s12 <- 99.1;     label("First term A of the aztreonam EC50 bi-exponential, C. freundii 974673 (mg/L)")   # Table 2: A 99.1 (RSE 6.8%); Code THETA(12)
    ec50a_atm_s15 <- 2.94;     label("First term A of the aztreonam EC50 bi-exponential, E. cloacae 1285905 (mg/L)")   # Table 2: A 2.94 (RSE 7.7%); Code THETA(13)
    ec50a_atm_s16 <- 67.9;     label("First term A of the aztreonam EC50 bi-exponential, E. cloacae 1318536 (mg/L)")   # Table 2: A 67.9 (RSE 3.3%); Code THETA(14)
    ec50b_atm_s1  <- fixed(0); label("Second term B of the aztreonam EC50 bi-exponential, E. coli 1266865 (mg/L)")     # Table 2: B 0 fix (mono-exponential for this strain); Code THETA(15) FIX
    ec50b_atm_s12 <- 0.494;    label("Second term B of the aztreonam EC50 bi-exponential, C. freundii 974673 (mg/L)")  # Table 2: B 0.494 (RSE 6.5%); Code THETA(16)
    ec50b_atm_s15 <- 0.241;    label("Second term B of the aztreonam EC50 bi-exponential, E. cloacae 1285905 (mg/L)")  # Table 2: B 0.241 (RSE 6.7%); Code THETA(17)
    ec50b_atm_s16 <- 0.402;    label("Second term B of the aztreonam EC50 bi-exponential, E. cloacae 1318536 (mg/L)")  # Table 2: B 0.402 (RSE 5.7%); Code THETA(18)
    alpha_s1  <- 0.42;     label("Avibactam exponential constant on the A term of the aztreonam EC50, E. coli 1266865 (L/mg)")      # Table 2: alpha 0.42 (RSE 1.5%); Code THETA(19)
    alpha_s12 <- 37.2;     label("Avibactam exponential constant on the A term of the aztreonam EC50, C. freundii 974673 (L/mg)")   # Table 2: alpha 37.2 (RSE 1.2%); Code THETA(20)
    alpha_s15 <- 25.9;     label("Avibactam exponential constant on the A term of the aztreonam EC50, E. cloacae 1285905 (L/mg)")   # Table 2: alpha 25.9 (RSE 3.1%); Code THETA(21)
    alpha_s16 <- 42.3;     label("Avibactam exponential constant on the A term of the aztreonam EC50, E. cloacae 1318536 (L/mg)")   # Table 2: alpha 42.3 (RSE 2.8%); Code THETA(22)
    beta_s1   <- fixed(0); label("Avibactam exponential constant on the B term of the aztreonam EC50, E. coli 1266865 (L/mg)")      # Table 2: beta 0 fix; Code THETA(23) FIX
    beta_s12  <- 0.548;    label("Avibactam exponential constant on the B term of the aztreonam EC50, C. freundii 974673 (L/mg)")   # Table 2: beta 0.548 (RSE 1.8%); Code THETA(24)
    beta_s15  <- 0.602;    label("Avibactam exponential constant on the B term of the aztreonam EC50, E. cloacae 1285905 (L/mg)")   # Table 2: beta 0.602 (RSE 3.4%); Code THETA(25)
    beta_s16  <- 0.412;    label("Avibactam exponential constant on the B term of the aztreonam EC50, E. cloacae 1318536 (L/mg)")   # Table 2: beta 0.412 (RSE 5.8%); Code THETA(26)

    # --- Aztreonam Hill coefficient gamma (Table 2; Code THETA 27-30) ---
    hill_atm_s1  <- 1.48; label("Hill coefficient of the aztreonam kill, E. coli 1266865 (unitless)")      # Table 2: gamma 1.48 (RSE 8.1%); Code THETA(27)
    hill_atm_s12 <- 1.42; label("Hill coefficient of the aztreonam kill, C. freundii 974673 (unitless)")   # Table 2: gamma 1.42 (RSE 5.6%); Code THETA(28)
    hill_atm_s15 <- 2.65; label("Hill coefficient of the aztreonam kill, E. cloacae 1285905 (unitless)")   # Table 2: gamma 2.65 (RSE 7.0%); Code THETA(29)
    hill_atm_s16 <- 1.52; label("Hill coefficient of the aztreonam kill, E. cloacae 1318536 (unitless)")   # Table 2: gamma 1.52 (RSE 3.7%); Code THETA(30)

    # --- Resistance factor: aztreonam EC50 multiplier for P2 (Table 2; Code THETA 31-34) ---
    resfac_s1  <- 8.09; label("Aztreonam EC50 increasing factor for P2, E. coli 1266865 (unitless)")      # Table 2: Resistance factor 8.09 (RSE 5.7%); Code THETA(31)
    resfac_s12 <- 9.39; label("Aztreonam EC50 increasing factor for P2, C. freundii 974673 (unitless)")   # Table 2: Resistance factor 9.39 (RSE 4.2%); Code THETA(32)
    resfac_s15 <- 18;   label("Aztreonam EC50 increasing factor for P2, E. cloacae 1285905 (unitless)")   # Table 2: Resistance factor 18 (RSE 7.2%); Code THETA(33)
    resfac_s16 <- 29.1; label("Aztreonam EC50 increasing factor for P2, E. cloacae 1318536 (unitless)")   # Table 2: Resistance factor 29.1 (RSE 4.2%); Code THETA(34)

    # --- Avibactam bactericidal effect on P1 (Table 2; Code THETA 35-46) ---
    emax_avi_s1  <- 2.28; label("Maximum kill rate constant of avibactam on P1, E. coli 1266865 (1/h)")      # Table 2: Emax,AVI 2.28 (RSE 2.0%); Code THETA(35)
    emax_avi_s12 <- 1.29; label("Maximum kill rate constant of avibactam on P1, C. freundii 974673 (1/h)")   # Table 2: Emax,AVI 1.29 (RSE 7.1%); Code THETA(36)
    emax_avi_s15 <- 1.69; label("Maximum kill rate constant of avibactam on P1, E. cloacae 1285905 (1/h)")   # Table 2: Emax,AVI 1.69 (RSE 7.0%); Code THETA(37)
    emax_avi_s16 <- 2.14; label("Maximum kill rate constant of avibactam on P1, E. cloacae 1318536 (1/h)")   # Table 2: Emax,AVI 2.14 (RSE 4.3%); Code THETA(38)
    ec50_avi_s1  <- 14.5; label("Avibactam concentration giving half the maximum kill rate on P1, E. coli 1266865 (mg/L)")      # Table 2: EC50,AVI 14.5 (RSE 2.8%); Code THETA(39)
    ec50_avi_s12 <- 17.3; label("Avibactam concentration giving half the maximum kill rate on P1, C. freundii 974673 (mg/L)")   # Table 2: EC50,AVI 17.3 (RSE 6.8%); Code THETA(40)
    ec50_avi_s15 <- 14.8; label("Avibactam concentration giving half the maximum kill rate on P1, E. cloacae 1285905 (mg/L)")   # Table 2: EC50,AVI 14.8 (RSE 14.7%); Code THETA(41)
    ec50_avi_s16 <- 12;   label("Avibactam concentration giving half the maximum kill rate on P1, E. cloacae 1318536 (mg/L)")   # Table 2: EC50,AVI 12 (RSE 5.1%); Code THETA(42)
    hill_avi_s1  <- 12.9; label("Hill coefficient of the avibactam kill, E. coli 1266865 (unitless)")      # Table 2: phi 12.9 (RSE 29.5%); Code THETA(43)
    hill_avi_s12 <- 4.63; label("Hill coefficient of the avibactam kill, C. freundii 974673 (unitless)")   # Table 2: phi 4.63 (RSE 15.1%); Code THETA(44)
    hill_avi_s15 <- 2.1;  label("Hill coefficient of the avibactam kill, E. cloacae 1285905 (unitless)")   # Table 2: phi 2.1 (RSE 15.3%); Code THETA(45)
    hill_avi_s16 <- 8.91; label("Hill coefficient of the avibactam kill, E. cloacae 1318536 (unitless)")   # Table 2: phi 8.91 (RSE 22.9%); Code THETA(46)

    # --- Aztreonam degradation by beta-lactamases (Table 2; Code THETA 47-52) ---
    # E. coli 1266865 has no degradation sub-model (Table 2 '-'): the code sets
    # TVDEGMIN = 0, TVPSI = 1, TVIC50 = 1 and TVPSI2 = 1 for STRAIN 1. With
    # degmin_s1 = 0 the other three values have no effect.
    degmin_s1  <- fixed(0);  label("Minimum degradation rate constant of aztreonam, E. coli 1266865 (1/h)")       # Table 2: '-' (no degradation observed); Code IF(STRAIN.EQ.1) TVDEGMIN = 0
    degmin_s12 <- 0.00462;   label("Minimum degradation rate constant of aztreonam, C. freundii 974673 (1/h)")    # Table 2: Degmin 0.00462 (RSE 8.9%); Code THETA(47)
    degmin_s15 <- 0.000129;  label("Minimum degradation rate constant of aztreonam, E. cloacae 1285905 (1/h)")    # Table 2: Degmin 0.000129 (RSE 2.9%); Code THETA(48)
    degmin_s16 <- 0.000166;  label("Minimum degradation rate constant of aztreonam, E. cloacae 1318536 (1/h)")    # Table 2: Degmin 0.000166 (RSE 2.7%); Code THETA(49)
    psi_s1  <- fixed(1);     label("Exponential constant linking log10 bacterial density to aztreonam degradation, E. coli 1266865 (per log10 CFU/mL)")      # Code IF(STRAIN.EQ.1) TVPSI = 1 (inert, degmin_s1 = 0)
    psi_s12 <- 0.563;        label("Exponential constant linking log10 bacterial density to aztreonam degradation, C. freundii 974673 (per log10 CFU/mL)")   # Table 2: psi 0.563 (RSE 2.2%); Code THETA(50)
    psi_s15 <- fixed(1);     label("Exponential constant linking log10 bacterial density to aztreonam degradation, E. cloacae 1285905 (per log10 CFU/mL)")   # Table 2: psi 1 fix; Code THETA(51) FIX
    psi_s16 <- fixed(1);     label("Exponential constant linking log10 bacterial density to aztreonam degradation, E. cloacae 1318536 (per log10 CFU/mL)")   # Table 2: psi 1 fix; Code THETA(52) FIX

    # --- Inhibition of aztreonam degradation by avibactam (Table 2; Code THETA 53-58) ---
    ic50_avi_s1  <- fixed(1);  label("Avibactam concentration halving the aztreonam degradation rate, E. coli 1266865 (mg/L)")      # Code IF(STRAIN.EQ.1) TVIC50 = 1 (inert, degmin_s1 = 0)
    ic50_avi_s12 <- 0.0000213; label("Avibactam concentration halving the aztreonam degradation rate, C. freundii 974673 (mg/L)")   # Table 2: IC50 0.0000213 (RSE 48.4%); Code THETA(53)
    ic50_avi_s15 <- 0.000185;  label("Avibactam concentration halving the aztreonam degradation rate, E. cloacae 1285905 (mg/L)")   # Table 2: IC50 0.000185 (RSE 22.0%); Code THETA(54)
    ic50_avi_s16 <- 0.0128;    label("Avibactam concentration halving the aztreonam degradation rate, E. cloacae 1318536 (mg/L)")   # Table 2: IC50 0.0128 (RSE 8.3%); Code THETA(55)
    hill_inh_s1  <- fixed(1);  label("Hill coefficient of the avibactam inhibition of aztreonam degradation, E. coli 1266865 (unitless)")      # Code IF(STRAIN.EQ.1) TVPSI2 = 1 (inert, degmin_s1 = 0)
    hill_inh_s12 <- 0.207;     label("Hill coefficient of the avibactam inhibition of aztreonam degradation, C. freundii 974673 (unitless)")   # Table 2: varphi 0.207 (RSE 6.5%); Code THETA(56)
    hill_inh_s15 <- 0.522;     label("Hill coefficient of the avibactam inhibition of aztreonam degradation, E. cloacae 1285905 (unitless)")   # Table 2: varphi 0.522 (RSE 6.9%); Code THETA(57)
    hill_inh_s16 <- 1.48;      label("Hill coefficient of the avibactam inhibition of aztreonam degradation, E. cloacae 1318536 (unitless)")   # Table 2: varphi 1.48 (RSE 5.9%); Code THETA(58)

    # --- Residual error (Table 2; Code THETA 59-66 used as SDs with $SIGMA 1 FIX) ---
    addSd_log_cfu_s1  <- 0.97;  label("Additive residual SD on the log10 bacterial count, E. coli 1266865 (log10 CFU/mL)")      # Table 2: 0.97 (RSE 2.8%); Code THETA(59), SIG = SQRT(THETA(59)**2)
    addSd_log_cfu_s12 <- 0.708; label("Additive residual SD on the log10 bacterial count, C. freundii 974673 (log10 CFU/mL)")   # Table 2: 0.708 (RSE 2.2%); Code THETA(60)
    addSd_log_cfu_s15 <- 0.845; label("Additive residual SD on the log10 bacterial count, E. cloacae 1285905 (log10 CFU/mL)")   # Table 2: 0.845 (RSE 2.8%); Code THETA(61)
    addSd_log_cfu_s16 <- 0.782; label("Additive residual SD on the log10 bacterial count, E. cloacae 1318536 (log10 CFU/mL)")   # Table 2: 0.782 (RSE 2.1%); Code THETA(62)
    propSd_s1  <- 0.101; label("Proportional residual SD on the aztreonam concentration, E. coli 1266865 (fraction)")      # Table 2: 10.1 percent (RSE 3.8%); Code THETA(63) = 0.101, W = IPRD2*SQRT(THETA(63)**2)
    propSd_s12 <- 0.494; label("Proportional residual SD on the aztreonam concentration, C. freundii 974673 (fraction)")   # Table 2: 49.4 percent (RSE 3.6%); Code THETA(64) = 0.494
    propSd_s15 <- 0.244; label("Proportional residual SD on the aztreonam concentration, E. cloacae 1285905 (fraction)")   # Table 2: 24.4 percent (RSE 4.1%); Code THETA(65) = 0.244
    propSd_s16 <- 0.333; label("Proportional residual SD on the aztreonam concentration, E. cloacae 1318536 (fraction)")   # Table 2: 33.3 percent (RSE 3.4%); Code THETA(66) = 0.333
  })

  model({
    # ---- 1. Strain-specific parameter selection (Code $PK IF(STRAIN.EQ.nn)) ----
    if (BACT == 1) {
      fmut <- fmut_s1
      kg <- kg_s1
      emax_atm <- emax_atm_s1
      ec50a_atm <- ec50a_atm_s1
      ec50b_atm <- ec50b_atm_s1
      alpha <- alpha_s1
      beta_avi <- beta_s1
      hill_atm <- hill_atm_s1
      resfac <- resfac_s1
      emax_avi <- emax_avi_s1
      ec50_avi <- ec50_avi_s1
      hill_avi <- hill_avi_s1
      degmin <- degmin_s1
      psi <- psi_s1
      ic50_avi <- ic50_avi_s1
      hill_inh <- hill_inh_s1
      sd_cfu <- addSd_log_cfu_s1
      sd_atm <- propSd_s1
    } else if (BACT == 12) {
      fmut <- fmut_s12
      kg <- kg_s12
      emax_atm <- emax_atm_s12
      ec50a_atm <- ec50a_atm_s12
      ec50b_atm <- ec50b_atm_s12
      alpha <- alpha_s12
      beta_avi <- beta_s12
      hill_atm <- hill_atm_s12
      resfac <- resfac_s12
      emax_avi <- emax_avi_s12
      ec50_avi <- ec50_avi_s12
      hill_avi <- hill_avi_s12
      degmin <- degmin_s12
      psi <- psi_s12
      ic50_avi <- ic50_avi_s12
      hill_inh <- hill_inh_s12
      sd_cfu <- addSd_log_cfu_s12
      sd_atm <- propSd_s12
    } else if (BACT == 15) {
      fmut <- fmut_s15
      kg <- kg_s15
      emax_atm <- emax_atm_s15
      ec50a_atm <- ec50a_atm_s15
      ec50b_atm <- ec50b_atm_s15
      alpha <- alpha_s15
      beta_avi <- beta_s15
      hill_atm <- hill_atm_s15
      resfac <- resfac_s15
      emax_avi <- emax_avi_s15
      ec50_avi <- ec50_avi_s15
      hill_avi <- hill_avi_s15
      degmin <- degmin_s15
      psi <- psi_s15
      ic50_avi <- ic50_avi_s15
      hill_inh <- hill_inh_s15
      sd_cfu <- addSd_log_cfu_s15
      sd_atm <- propSd_s15
    } else if (BACT == 16) {
      fmut <- fmut_s16
      kg <- kg_s16
      emax_atm <- emax_atm_s16
      ec50a_atm <- ec50a_atm_s16
      ec50b_atm <- ec50b_atm_s16
      alpha <- alpha_s16
      beta_avi <- beta_s16
      hill_atm <- hill_atm_s16
      resfac <- resfac_s16
      emax_avi <- emax_avi_s16
      ec50_avi <- ec50_avi_s16
      hill_avi <- hill_avi_s16
      degmin <- degmin_s16
      psi <- psi_s16
      ic50_avi <- ic50_avi_s16
      hill_inh <- hill_inh_s16
      sd_cfu <- addSd_log_cfu_s16
      sd_atm <- propSd_s16
    } else {
      # An unrecognised strain code propagates NA through every output
      # rather than silently borrowing another strain's parameters.
      fmut <- NA
      kg <- NA
      emax_atm <- NA
      ec50a_atm <- NA
      ec50b_atm <- NA
      alpha <- NA
      beta_avi <- NA
      hill_atm <- NA
      resfac <- NA
      emax_avi <- NA
      ec50_avi <- NA
      hill_avi <- NA
      degmin <- NA
      psi <- NA
      ic50_avi <- NA
      hill_inh <- NA
      sd_cfu <- NA
      sd_atm <- NA
    }

    bmax <- 10^log10bmax
    inoc <- 10^log10inoc

    # ---- 2. Avibactam enhancement: bi-exponential aztreonam EC50 (Code $DES EC501) ----
    ec50_atm <- ec50a_atm * exp(-alpha * CONC_AVI_MGL) +
      ec50b_atm * exp(-beta_avi * CONC_AVI_MGL)

    # ---- 3. Kill rates (Code $DES ATMS, ATMR, AVIS) ----
    # Aztreonam kills both subpopulations; the P2 EC50 is multiplied by the
    # resistance factor. Avibactam kills P1 only.
    kill_atm_p1 <- emax_atm * central^hill_atm /
      (ec50_atm^hill_atm + central^hill_atm)
    kill_atm_p2 <- emax_atm * central^hill_atm /
      ((ec50_atm * resfac)^hill_atm + central^hill_atm)
    kill_avi_p1 <- emax_avi * CONC_AVI_MGL^hill_avi /
      (ec50_avi^hill_avi + CONC_AVI_MGL^hill_avi)

    # ---- 4. Logistic growth on the total count (Code $DES PLAT) ----
    ntot <- bact_susceptible_susceptible + bact_intermediate_resistant
    plat <- 1 - ntot / bmax

    # ---- 5. Aztreonam degradation (Code $DES NB, DEG, INH) ----
    # The degradation rate constant is Degmin * exp(psi * log10(P1 + P2)),
    # inhibited by avibactam through a fractional inhibitory Emax term.
    # The count is floored at a tiny positive value inside log10() so that a
    # negative solver excursion of a fully killed culture cannot produce NaN;
    # for any positive count the expression is the code's exactly.
    deg <- degmin * exp(log10(max(ntot, 1e-30)) * psi)
    inh <- 1 - CONC_AVI_MGL^hill_inh / (ic50_avi^hill_inh + CONC_AVI_MGL^hill_inh)

    # ---- 6. ODEs (Code $DES DADT(1)-DADT(3)) ----
    d/dt(bact_susceptible_susceptible) <- kg * plat * bact_susceptible_susceptible -
      (kill_atm_p1 + kill_avi_p1) * bact_susceptible_susceptible
    d/dt(bact_intermediate_resistant) <- kg * plat * bact_intermediate_resistant -
      kill_atm_p2 * bact_intermediate_resistant
    d/dt(central) <- -deg * inh * central

    # Initial conditions (Code $PK A_0(1), A_0(2)); aztreonam enters as a dose.
    bact_susceptible_susceptible(0) <- inoc * (1 - fmut)
    bact_intermediate_resistant(0) <- inoc * fmut

    # ---- 7. Observations (Code $ERROR IPRD1, IPRD2) ----
    # Total count on log10 CFU/mL with the 1e-6 CFU/mL offset of the code
    # (negative solver excursions clamped to 0 first). The fitted counts were
    # censored below 200 CFU/mL (M3 method).
    log_cfu <- log10(max(ntot, 0) + 0.000001)
    Cc <- central

    log_cfu ~ add(sd_cfu)
    Cc ~ prop(sd_atm)
  })
}
