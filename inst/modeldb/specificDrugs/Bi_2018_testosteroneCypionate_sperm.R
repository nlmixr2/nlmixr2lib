Bi_2018_testosteroneCypionate_sperm <- function() {
  description <- "Indirect response model for suppression of spermatogenesis during and after 14 weekly intramuscular injections of depot testosterone cypionate in 29 healthy men. Sperm count is produced at a zero-order rate that is inhibited by a sigmoid Emax function of the average total testosterone concentration over the preceding 18 weeks, supplied as the CAV exposure covariate from the companion PK model. The maximum fractional inhibition is estimated on the logit scale so it stays within (0, 1), and the inhibitory potency carries a baseline body-weight power effect."
  reference <- paste(
    "Bi Y, Perry PJ, Ellerby M, Murry DJ.",
    "Population Pharmacokinetic/Pharmacodynamic Modeling of Depot Testosterone Cypionate in Healthy Male Subjects.",
    "CPT Pharmacometrics Syst Pharmacol 2018;7(4):259-268.",
    "doi:10.1002/psp4.12287.",
    "Structural equations from the sperm-count NONMEM control stream in Supplementary Material S6;",
    "final parameter values from Table 3 (Sperm count block).",
    "The CAV exposure covariate is produced by modellib('Bi_2018_testosteroneCypionate').",
    sep = " "
  )
  vignette <- "Bi_2018_testosteroneCypionate"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  # The sperm pool is the only ODE state; the drug enters only through the CAV
  # exposure covariate, exactly as in the source paper's sequential fit.
  paper_specific_compartments <- c("sperm")

  compartmentData <- list(
    sperm = list(analyte = "spermatozoa", units = "million/mL", specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    CAV = list(
      description = "Average total testosterone concentration over the 18 weeks preceding the current time.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Bi 2018 tested rolling averages from 3 to 40 weeks (Supplementary S6 $INPUT columns CAVGP3 ... CAVGP40) and retained the 18-week window (CAVGP18) as the one giving the lowest objective function value. Averaging convention, from Bi 2018 Results: for observations earlier than 18 weeks the window runs from the start of the study to the time of measurement, i.e. CAV(t) = (AUC(t) - AUC(max(0, t - 126 days))) / min(t, 126 days). Values are NOT zero before dosing - they are the pre-treatment endogenous testosterone average (the Supplementary S9 example subject carries 4.868 ng/mL at time 0). In the source analysis the values came from the individual post hoc predictions of the companion population PK model; reproduce them with modellib('Bi_2018_testosteroneCypionate').",
      source_name = "CAVGP18"
    ),
    WT_BASE = list(
      description = "Per-subject baseline body weight, time-fixed.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters the inhibitory potency as a power function. Bi 2018 Results states the effect against a typical subject of median weight 85 kg: 'The Cavg50 of a heavier subject (95th percentile, 110 kg) was estimated to be 0.72 (95% CI 0.55-0.91) fold of that of a typical subject with median weight (85 kg)', and (110/85)^-1.27 = 0.721. The deposited Supplementary S6 sperm-count stream prints (WT/70); see the vignette Errata for why 85 kg is used here. Table 1 group medians 82.2 / 88.8 / 84.7 kg.",
      source_name = "WT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 29L,
    n_studies = 1L,
    age_range = "21-39 years",
    weight_range = "60.7-115 kg",
    sex_female_pct = 0,
    disease_state = "healthy men",
    dose_range = "100, 250 or 500 mg testosterone cypionate intramuscularly once weekly for 14 consecutive weeks (study weeks 2-15)",
    n_observations = 159L,
    follow_up = "40 weeks",
    regions = "United States",
    notes = "Semen samples were obtained at baseline (week 0) and at weeks 2, 16, 21, 28 and, in some subjects, week 40. Two of the 31 enrolled subjects were excluded from this endpoint because prior vasectomy left them unable to provide a sample, giving 29 subjects and 159 samples. 11 subjects missed 15 of 174 planned collections (8.6%) and 7 subjects (24.1%) missed the final collection, which Bi 2018 flags as a possible source of bias in the end-of-study recovery estimate. Sperm count was used in preference to sperm motility because the two were highly correlated (adjusted R^2 = 0.87); the motility fit is a sensitivity analysis reported in Supplementary Table S1 and is not extracted here."
  )

  ini({
    # --- Sperm turnover (Bi 2018 Table 3, Sperm count block) ---
    lkin_sperm <- log(6.17); label("Zero-order sperm production rate (million/mL/day)") # Table 3, Kin 6.17 /day
    lkout_sperm <- log(0.0696); label("First-order sperm loss rate constant (1/day)") # Table 3, Kout 0.0696 /day

    # --- Sigmoid Emax inhibition of sperm production by testosterone exposure ---
    logitemax_sperm <- 4.65; label("Logit of the maximum fractional inhibition of sperm production (unitless)") # Table 3, 'Phi (phi)' 4.65; Bi 2018 Results: Emax = exp(PHI)/(1 + exp(PHI)), so Emax = 0.9905
    lic50_sperm <- log(8.68); label("18-week average testosterone concentration giving half-maximal inhibition of sperm production (ng/mL)") # Table 3, 'Cavg50' 8.68 ng/mL
    lhill_sperm <- log(11.3); label("Hill coefficient for inhibition of sperm production (unitless)") # Table 3, 'Lambda; k' 11.3

    # --- Covariate effect ---
    e_bwt_ic50_sperm <- -1.27; label("Power exponent of (WT_BASE / 85) on the inhibitory potency (unitless)") # Table 3, 'h Bwt-Cavg50' -1.27

    # --- Inter-individual variability ---
    # Table 3 reports IIV as 100 * sqrt(omega), so omega is the squared percentage.
    etalkin_sperm ~ 0.219024 # Table 3 'IIV_Kin' 46.8%; 0.468^2
    etalogitemax_sperm ~ 5.4756 # Table 3 'IIV_Phi (phi)' 234%; 2.34^2. Enters ADDITIVELY on the logit scale: S6 codes PHI = TVPHI + ETA(2)
    etalic50_sperm ~ 0.053824 # Table 3 'IIV_Cavg50' 23.2%; 0.232^2
    etalhill_sperm ~ 1 # Table 3 'IIV_k' 100%; 1.00^2

    # --- Residual error ---
    # S6 sets W = THETA(6) directly on log(sperm + 1) with $SIGMA 1 FIX, so Table 3's
    # 'sigma^2 additive' row is a standard deviation despite its label.
    addSd_logSperm <- 0.429; label("Additive residual error on log(sperm count + 1) (log million/mL)") # Table 3, 'sigma^2 additive' 0.429
  })

  model({
    # --- Individual parameters ------------------------------------------------
    kin_sperm <- exp(lkin_sperm + etalkin_sperm)
    kout_sperm <- exp(lkout_sperm)
    ic50_sperm <- exp(lic50_sperm + etalic50_sperm) * (WT_BASE / 85)^e_bwt_ic50_sperm
    hill_sperm <- exp(lhill_sperm + etalhill_sperm)
    # Logistic transformation constraining the maximum inhibition to (0, 1).
    phi_sperm <- logitemax_sperm + etalogitemax_sperm
    emax_sperm <- expit(phi_sperm)

    # --- Initial condition ----------------------------------------------------
    # S6 sperm-count stream: A_0(1) = KIN / KOUT, the unsuppressed steady state.
    sperm(0) <- kin_sperm / kout_sperm

    # --- Indirect response ----------------------------------------------------
    cav <- max(CAV, 0)
    # S6: INH = 1 - EMAX * CAVGP18**LAM / (CAVGP18**LAM + LH50**LAM), then
    # DADT(1) = KIN * INH - KOUT * A(1). Note that INH there is the fraction of
    # production REMAINING, not the fraction inhibited, so it multiplies KIN
    # directly rather than entering as (1 - INH).
    remaining <- 1 - emax_sperm * cav^hill_sperm / (cav^hill_sperm + ic50_sperm^hill_sperm)
    d/dt(sperm) <- kin_sperm * remaining - kout_sperm * sperm

    # --- Observation ----------------------------------------------------------
    # S6: IPRED = LOG(A(1) + 1); the +1 offset matters because sperm counts reach
    # zero under complete suppression in the 250 and 500 mg groups.
    logSperm <- log(sperm + 1)
    logSperm ~ add(addSd_logSperm)
  })
}
