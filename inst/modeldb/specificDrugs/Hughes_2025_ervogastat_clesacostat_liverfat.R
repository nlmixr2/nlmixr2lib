Hughes_2025_ervogastat_clesacostat_liverfat <- function() {
  description <- "Joint exposure-response model for two liver-fat imaging endpoints, magnetic resonance imaging proton density fat fraction (MRI-PDFF, %) and FibroScan controlled attenuation parameter (CAP, dB/m), in adults with non-alcoholic fatty liver disease or presumed non-alcoholic steatohepatitis given ervogastat (DGAT2 inhibitor) alone or with clesacostat (ACC inhibitor). Each bounded endpoint is mapped to an unbounded scale with an Aranda-Ordaz dynamic transform-both-sides (CAP with flipped boundaries), and both transformed endpoints load on one shared latent liver-fat variable (common factor model) with an endpoint-specific baseline and scaling and an additional unique baseline random effect for CAP. The latent variable carries an exponential-onset placebo effect and a shared exponential-onset drug effect: an on/off ervogastat treatment effect that co-administered clesacostat reduces by 73.5%, plus a hyperbolic Emax effect of the individual average steady-state clesacostat concentration. Covariates: baseline ALT and age on the shared baseline, body weight on the CAP baseline, and study design and PNPLA3 genotype on the placebo effect. Observations are the transformed endpoints with additive residual error; back-transformed pdff and cap are also output."
  reference <- paste(
    "Hughes JH, Amin NB, Wojciechowski J, Vourvahis M.",
    "Exposure-response modeling of liver fat imaging endpoints in",
    "non-alcoholic fatty liver disease populations administered ervogastat",
    "alone and co-administered with clesacostat.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(2):317-330.",
    "doi:10.1002/psp4.13275",
    sep = " "
  )
  vignette <- "Hughes_2025_ervogastat_clesacostat"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CAV_ERVOGASTAT = list(
      description = "Individual average steady-state ervogastat plasma concentration",
      units = "ng/mL",
      type = "continuous",
      reference_category = "0 (no ervogastat; placebo or clesacostat monotherapy)",
      notes = paste(
        "Source CAVSS1, the individual predicted average steady-state",
        "ervogastat concentration from the companion population PK model",
        "(modellib('Hughes_2025_ervogastat'); Cav,ss = F * dose / (CL/F *",
        "tau)). Enters only as an on/off switch, CAVSS1 > 0, because the",
        "ervogastat EC50 was not identifiable and the Emax model was reduced",
        "to a treatment effect (Results; Table S6 run 9). Constant within a",
        "subject's treatment; 0 for placebo and clesacostat-alone arms."
      ),
      source_name = "CAVSS1"
    ),
    CAV_CLESACOSTAT = list(
      description = "Individual average steady-state clesacostat plasma concentration",
      units = "ng/mL",
      type = "continuous",
      reference_category = "0 (no clesacostat; placebo or ervogastat monotherapy)",
      notes = paste(
        "Source CAVSS2, the individual predicted average steady-state",
        "clesacostat concentration from the companion population PK model",
        "(modellib('Hughes_2025_clesacostat')). Drives a hyperbolic Emax",
        "term, Emax * CAV / (CAV + EC50) with EC50 = 389.306 ng/mL, and as an",
        "on/off switch (CAVSS2 > 0) the 73.5% inhibition of the ervogastat",
        "treatment effect (Code S3 DEFF2, DINT2). 0 for placebo and",
        "ervogastat-alone arms."
      ),
      source_name = "CAVSS2"
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source BALT. Log-linear additive effect on the shared (latent)",
        "baseline, (log(ALT) - log(40)) * 0.575256, applied to both endpoints",
        "through their scaling parameters (Code S3 BALTBSLN; Table 2 'Shared",
        "effect of baseline ALT on h(z(BSLpsi)); ref. 40 IU'). Skipped when",
        "BALT is 0 (missing values were imputed with the median, Methods)."
      ),
      source_name = "BALT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Log-linear additive effect on the shared baseline, (log(AGE) -",
        "log(51)) * -0.563598 (Code S3 AGEBSLN; Table 2 'Shared effect of",
        "age on h(z(BSLpsi)); ref. 51 years')."
      ),
      source_name = "AGE"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source BWT. Log-linear additive effect on the CAP baseline only,",
        "(log(WT) - log(92)) * 3.80204 (Code S3 BWTBSLN, applied when FLAG =",
        "169, the CAP records; Table 2 'Effect of baseline weight on",
        "h(z(BSLCAP)); ref. 92 kg')."
      ),
      source_name = "BWT"
    ),
    STUDY_NCT03513588 = list(
      description = "Phase 1b ervogastat study NCT03513588 (protocol C2541005) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the three outpatient phase 2a studies)",
      notes = paste(
        "1 = participant of NCT03513588, the 2-week inpatient phase 1b study",
        "of ervogastat in adults with NAFLD (Table S1), 0 otherwise. Code S3",
        "keys the effect on the protocol (PROT = 2541005); the paper",
        "describes it as the 'effect of inpatient study design on Pmax'",
        "(-2.74058, Table 2) and attributes the larger placebo response to",
        "the inpatient design and milder disease (Discussion)."
      ),
      source_name = "PROT"
    ),
    STUDY_PHASE2 = list(
      description = "Phase 2a study indicator (studies with a 2-week placebo run-in)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (NCT03513588, the only non-phase-2 study in the exposure-response data)",
      notes = paste(
        "1 = record from NCT03248882, NCT03776175 or NCT04399538 (Code S3",
        "PROT = 1171002, 3711001 or 3711005), the three phase 2a studies,",
        "each of which had a 2-week placebo run-in (Table S1). Selects the",
        "placebo-effect clock: time after the baseline observation for these",
        "studies, time after first dose otherwise (Code S3 PTIME)."
      ),
      source_name = "PROT"
    ),
    T_SCAN_TO_DOSE = list(
      description = "Time from the baseline liver-fat assessment to the first study dose",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Model time is time after the baseline observation (Code S3",
        "TIME = TABO). T_SCAN_TO_DOSE = TABO - TAFD, so the drug-effect clock",
        "is max(0, time - T_SCAN_TO_DOSE) (Code S3 DTIME = TAFD). Baseline",
        "imaging was on Day -14 in NCT03248882 and NCT04399538, Day -4 in",
        "NCT03776175 and Day -2 in NCT03513588 (Table S1), i.e. about 336,",
        "96 and 48 h."
      ),
      source_name = "TABO - TAFD"
    ),
    PNPLA3_CC = list(
      description = "PNPLA3 rs738409 C/C genotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (with PNPLA3_CG and PNPLA3_GG also 0: genotype not determined)",
      notes = paste(
        "Code S3 PNPLA3 = 1 (wild type C/C); -1.00284 on Pmax (Table 2). In",
        "this model the reference is the 35.1% of participants without a",
        "PNPLA3 genotype, not C/C, because the placebo response of the",
        "undetermined group differed substantially from C/C (Methods,",
        "Covariate model)."
      ),
      source_name = "PNPLA3"
    ),
    PNPLA3_CG = list(
      description = "PNPLA3 rs738409 C/G heterozygote indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (with PNPLA3_CC and PNPLA3_GG also 0: genotype not determined)",
      notes = "Code S3 PNPLA3 = 2; -0.694264 on Pmax (Table 2). Reference is the undetermined genotype (see PNPLA3_CC).",
      source_name = "PNPLA3"
    ),
    PNPLA3_GG = list(
      description = "PNPLA3 rs738409 G/G homozygote indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (with PNPLA3_CC and PNPLA3_CG also 0: genotype not determined)",
      notes = "Code S3 PNPLA3 = 3; -0.0834581 on Pmax (Table 2). Reference is the undetermined genotype (see PNPLA3_CC).",
      source_name = "PNPLA3"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 524L,
    n_studies = 4L,
    n_observations = 3100L,
    age_range = "19-72 years",
    age_median = "55 years",
    weight_range = "55-170 kg",
    weight_median = "93.8 kg",
    sex_female_pct = 50,
    race_ethnicity = c(White = 82.8, Black = 4, Asian = 10.7, Other = 2.5),
    disease_state = "Adults with NAFLD (46.6%) or presumed NASH (53.4%)",
    dose_range = paste(
      "Ervogastat 25-300 mg twice daily or 300 mg once daily and clesacostat",
      "2-50 mg once daily, 2.5-15 mg twice daily or 20 mg once daily, alone or",
      "combined, for 2-16 weeks, plus placebo (Table S1)"
    ),
    regions = "United States",
    notes = paste(
      "Demographics from Hughes 2025 Table 1 (exposure-response column):",
      "baseline ALT median 40 IU/L (6-341), BMI median 33.6 kg/m^2. PNPLA3",
      "C/C 22.7%, C/G 27.1%, G/G 15.1%, not determined 35.1%. 1398 MRI-PDFF",
      "(none at a boundary) and 1702 CAP measurements (3 at 100 and 66 at",
      "400 dB/m); the 3 lower-boundary CAP values were removed as outliers",
      "before the final estimation. Studies: NCT03513588 (phase 1b, 2",
      "weeks), NCT03248882 (16 weeks), NCT03776175 and NCT04399538 (6",
      "weeks). Estimated with Laplacian conditional estimation and a",
      "user-defined -2 log-likelihood in NONMEM 7.5 (Code S3, run46)."
    )
  )

  ini({
    # Transformed-scale baselines and transform coefficients, fixed from the
    # unstructured longitudinal model (Methods; Table 2 'Fixed')
    bsl_hpdff <- fixed(-1.16); label("Typical baseline MRI-PDFF on the transformed scale, h(z(BSL_MRI)) (unitless)") # Code S3 THETA(1) TVMRITBSLN = -1.16 FIX; Table 2
    bsl_hcap <- fixed(5.05); label("Typical baseline CAP on the flipped transformed scale, h(z(BSL_CAP)) (unitless)") # Code S3 THETA(2) TVCAPTBSLN = 5.05 FIX; Table 2
    alpha_pdff <- fixed(4.92); label("Aranda-Ordaz transform coefficient for MRI-PDFF (unitless)") # Code S3 THETA(10) TVMRIALPHA = 4.92 FIX; Table 2
    alpha_cap <- fixed(32.1); label("Aranda-Ordaz transform coefficient for CAP (unitless)") # Code S3 THETA(11) TVCAPALPHA = 32.1 FIX; Table 2

    # Scaling (factor loading) of the latent variable onto each endpoint
    s_pdff <- 0.698341; label("Scaling of the latent variable onto transformed MRI-PDFF (unitless)") # Code S3 THETA(14) TVMRILOAD = 0.698341; Table 2 's MRI' = 0.698
    s_cap_neg <- 2.76926; label("Negative of the scaling of the latent variable onto flipped transformed CAP (unitless)") # Code S3 THETA(15) TVCAPLOADNEG = 2.76926; Table 2 '-s CAP' = 2.77

    # Longitudinal placebo and drug effects on the latent variable
    lthalf_pbo <- log(49.0785); label("Half-life of the placebo-effect onset (day)") # Code S3 THETA(3) TVPBOKHALF = 49.0785; Table 2 't1/2,kpbo' = 49.1
    lthalf_drug <- log(18.8237); label("Shared half-life of the drug-effect onset (day)") # Code S3 THETA(4) TVDRGKHALF = 18.8237; Table 2 't1/2,ktrt' = 18.8
    pmax <- 0.393304; label("Maximum placebo effect on the latent variable (unitless)") # Code S3 THETA(5) TVPMAX = 0.393304; Table 2 'Pmax' = 0.393
    etrt_ervogastat <- 1.27219; label("Ervogastat treatment effect on the latent variable (unitless)") # Code S3 THETA(6) TVDGATEMAX = 1.27219; Table 2 'Etrt,er' = 1.27
    emax_clesacostat <- 2.32486; label("Clesacostat maximum effect on the latent variable (unitless)") # Code S3 THETA(7) TVACCIEMAX = 2.32486; Table 2 'Eemax,cl' = 2.32
    lec50_clesacostat <- log(389.306); label("Clesacostat average concentration giving half the maximum effect (ng/mL)") # Code S3 THETA(8) TVACCIEC50 = 389.306; Table 2 'EC50,cl' = 389
    imax_clesacostat <- -0.735388; label("Fractional change in the ervogastat treatment effect when clesacostat is co-administered (fraction)") # Code S3 THETA(9) TVACCIMAX = -0.735388; Table 2 'Itrt,cl-er' = -0.735

    # Covariate effects on the latent variable (additive, transformed scale
    # before scaling; Equations 11 and 12)
    e_alt_bsl <- 0.575256; label("Effect of log(ALT/40) on the shared baseline (unitless)") # Code S3 THETA(16) TVBALTBSLN = 0.575256; Table 2 = 0.575
    e_wt_bsl_cap <- 3.80204; label("Effect of log(WT/92) on the CAP baseline (unitless)") # Code S3 THETA(17) TVBWTBSLN = 3.80204; Table 2 = 3.80
    e_inpatient_pmax <- -2.74058; label("Effect of the inpatient NCT03513588 study design on Pmax (unitless)") # Code S3 THETA(18) TVINPATPMAX = -2.74058; Table 2 = -2.74
    e_pnpla3cc_pmax <- -1.00284; label("Effect of PNPLA3 C/C on Pmax (unitless)") # Code S3 THETA(19) TVPNPLA3CCPMAX = -1.00284; Table 2 = -1.00
    e_pnpla3cg_pmax <- -0.694264; label("Effect of PNPLA3 C/G on Pmax (unitless)") # Code S3 THETA(20) TVPNPLA3CGPMAX = -0.694264; Table 2 = -0.694
    e_pnpla3gg_pmax <- -0.0834581; label("Effect of PNPLA3 G/G on Pmax (unitless)") # Code S3 THETA(21) TVPNPLA3GGPMAX = -0.0834581; Table 2 = -0.0835
    e_age_bsl <- -0.563598; label("Effect of log(AGE/51) on the shared baseline (unitless)") # Code S3 THETA(22) TVAGEBSLN = -0.563598; Table 2 = -0.564

    # Random effects on the latent scale (Code S3 $OMEGA variances; Table 2
    # prints the SDs 1.00 FIX, 1.33, 0.959)
    etabsl ~ fixed(1) # Code S3 $OMEGA PPVBSLN = 1; Table 2 'Shared baseline liver fat' SD 1.00
    etabsl_hcap ~ 1.7795 # Code S3 $OMEGA PPVCAPBSLN = 1.7795; Table 2 'Unique baseline CAP' SD 1.33
    etapmax ~ 0.918996 # Code S3 $OMEGA PPVPMAX = 0.918996; Table 2 'Maximum placebo effect' SD 0.959

    # Residual SDs on the transformed scale (EPS variance 1 FIX)
    addSd_hpdff <- 0.204681; label("Residual SD of transformed MRI-PDFF (unitless)") # Code S3 THETA(12) TVMRISD = 0.204681; Table 2 'sigma MRI' = 0.205
    addSd_hcap <- 5.35106; label("Residual SD of flipped transformed CAP (unitless)") # Code S3 THETA(13) TVCAPSD = 5.35106; Table 2 'sigma CAP' = 5.35
  })

  model({
    # Clocks (h). Model time is time after the baseline observation.
    tafd <- t - T_SCAN_TO_DOSE
    dtime <- tafd
    if (dtime < 0) dtime <- 0
    # Placebo clock: time after baseline observation in the three phase 2a
    # studies with a placebo run-in, time after first dose otherwise
    ptime <- tafd
    if (STUDY_PHASE2 == 1) ptime <- t
    if (ptime < 0) ptime <- 0

    # Endpoint scalings; CAP uses flipped boundaries, so its scaling is the
    # negative of the estimated parameter (Code S3 LOAD = -LOAD)
    s_cap <- -s_cap_neg

    # Covariates on the latent scale (Code S3 COVBSLN, COVPMAX)
    cov_alt <- 0
    if (ALT > 0) cov_alt <- (log(ALT) - log(40)) * e_alt_bsl
    cov_age <- 0
    if (AGE > 0) cov_age <- (log(AGE) - log(51)) * e_age_bsl
    cov_wt <- 0
    if (WT > 0) cov_wt <- (log(WT) - log(92)) * e_wt_bsl_cap
    cov_pmax <- e_inpatient_pmax * STUDY_NCT03513588 +
      e_pnpla3cc_pmax * PNPLA3_CC + e_pnpla3cg_pmax * PNPLA3_CG + e_pnpla3gg_pmax * PNPLA3_GG

    # Individual baselines and placebo maxima on each transformed scale
    tbsl_pdff <- bsl_hpdff + s_pdff * (etabsl + cov_alt + cov_age)
    tbsl_cap <- bsl_hcap + s_cap * (etabsl + etabsl_hcap + cov_alt + cov_age + cov_wt)
    pmax_i <- pmax + etapmax + cov_pmax

    # Exponential onset of the placebo and drug effects (half-lives in days)
    kpbo <- log(2) / (exp(lthalf_pbo) * 24)
    kdrug <- log(2) / (exp(lthalf_drug) * 24)
    pbo_prof <- 1 - exp(-kpbo * ptime)
    drug_prof <- 1 - exp(-kdrug * dtime)

    # Drug effects on the latent variable (Code S3 DEFF1, DEFF2, DINT2)
    ec50_clesacostat <- exp(lec50_clesacostat)
    eff_cles <- emax_clesacostat * CAV_CLESACOSTAT / (CAV_CLESACOSTAT + ec50_clesacostat)
    on_cles <- 0
    if (CAV_CLESACOSTAT > 0) on_cles <- 1
    on_ervo <- 0
    if (CAV_ERVOGASTAT > 0) on_ervo <- 1
    eff_ervo <- etrt_ervogastat * on_ervo * (1 + imax_clesacostat * on_cles)
    drug_latent <- (eff_ervo + eff_cles) * drug_prof

    # The placebo effect (whose sign depends on the study and PNPLA3
    # covariates) adds to, and the drug effect subtracts from, the latent
    # variable; each is mapped to an endpoint through its scaling (Code S3
    # PMAX = LOAD * (...), EMAX = -LOAD * THETA)
    hpdff <- tbsl_pdff + s_pdff * (pmax_i * pbo_prof - drug_latent)
    hcap <- tbsl_cap + s_cap * (pmax_i * pbo_prof - drug_latent)

    # Back-transform to the original scales (Equation 1 inverted; MRI-PDFF
    # bounds 0-100%, CAP bounds 100-400 dB/m, CAP flipped as b - (y - a))
    pdff <- 100 * (1 - exp(log(alpha_pdff * exp(hpdff) + 1) / (-alpha_pdff)))
    cap_flipped <- 300 * (1 - exp(log(alpha_cap * exp(hcap) + 1) / (-alpha_cap))) + 100
    cap <- 400 - (cap_flipped - 100)

    hpdff ~ add(addSd_hpdff)
    hcap ~ add(addSd_hcap)
  })
}
