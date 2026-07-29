Netterberg_2017_docetaxel <- function() {
  description <- "Friberg-style semi-mechanistic myelosuppression PD model for docetaxel-induced neutropenia in adult cancer patients (DDMODEL00000224, Netterberg 2017 / Kloft 2006). The bundle's NM-TRAN .mod fixes parameter values at the Kloft 2006 docetaxel myelosuppression analysis (per the .mod's `; Parameter estimates as according to Kloft et al., 2006` block; MAXEVALS=0 in the bundle's $ESTIMATION confirms no re-fit) and Netterberg 2017 reuses the model unchanged to evaluate frequent-monitoring ANC prediction methodology. Structurally, the model has a self-renewing proliferation pool plus three transit compartments and a circulating compartment; docetaxel exposure is supplied as a time-varying plasma-concentration covariate (CP_MGL, mg/L) that drives a linear drug effect (1 - SL * CP_MGL) on proliferation, with a feedback term (BA / circ)^PO. Covariate effects on baseline ANC (sex, ECOG performance status, prior anticancer therapy, alpha-1 acid glycoprotein) and on the drug-effect slope (alpha-1 acid glycoprotein) follow Kloft 2006. Output is the absolute neutrophil count ANC in 10^9 cells/L."
  reference <- paste(
    "Netterberg I, Nielsen E I, Friberg L E, Karlsson M O. (2017).",
    "Model-based prediction of myelosuppression and recovery based on",
    "frequent neutrophil monitoring.",
    "Cancer Chemother Pharmacol 80(2):343-353.",
    "doi:10.1007/s00280-017-3366-x.",
    "Parameter values inherited unchanged from the docetaxel arm of the",
    "Kloft et al. 2006 cross-drug myelosuppression analysis",
    "(per the bundle's NM-TRAN .mod $PK / $THETA `; ... according to",
    "Kloft et al., 2006` block; the original Kloft 2006 publication is",
    "not on disk in this worktree).",
    "DDMORE Foundation Model Repository: DDMODEL00000224.",
    sep = " "
  )
  vignette <- "Netterberg_2017_docetaxel"
  # `dosing = "mg"` reflects the docetaxel mass-based dosing used by the upstream popPK model
  # that supplies CP_MGL (Kloft 2006 / Netterberg 2017 reference: 100 mg/m^2 IV over 1 h).
  # The PD-only model itself has no dosing compartment; CP_MGL is read directly from the
  # time-varying covariate column.
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L", anc = "10^9 cells/L", aag = "g/L")
  ddmore_id    <- "DDMODEL00000224"
  replicate_of <- NULL

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplicative effect on baseline ANC: BACOV *= (1 + e_sexf_ba * SEXF). Female subjects have ~12.1% lower baseline ANC than males (Kloft 2006 / Netterberg 2017 .mod THETA(11) = -0.121451). Source NM-TRAN column `SEX` uses Kloft 2006's 1 = male, 2 = female encoding; decompose to canonical SEXF via `SEXF = as.integer(SEX == 2)`. Time-fixed per subject.",
      source_name        = "SEX"
    ),
    ECOG_GE1 = list(
      description        = "Eastern Cooperative Oncology Group performance status indicator, 1 = ECOG >= 1, 0 = ECOG = 0.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG = 0; fully active / asymptomatic)",
      notes              = "Multiplicative effect on baseline ANC: BACOV *= (1 + e_ecogge1_ba * ECOG_GE1). Subjects with ECOG >= 1 have ~13.0% higher baseline ANC than ECOG = 0 subjects (Kloft 2006 / Netterberg 2017 .mod THETA(8) = 0.130406). Source NM-TRAN column `PERF` is the ordinal ECOG / WHO performance score; the .mod's $PK block binarizes via `IF(PERF.EQ.0.OR.PERF.EQ.-99) BAPERF = 0; IF(PERF.GE.1) BAPERF = THETA(8)`. Decompose to canonical via `ECOG_GE1 = as.integer(PERF >= 1)`. Time-fixed per subject.",
      source_name        = "PERF"
    ),
    PRIOR_ANTICANCER = list(
      description        = "Prior anticancer therapy of any modality (chemotherapy, radiotherapy, hormonal, targeted, immunotherapy, surgery), 1 = had prior anticancer therapy, 0 = treatment-naive.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (treatment-naive)",
      notes              = "Multiplicative effect on baseline ANC: BACOV *= (1 + e_pc_ba * PRIOR_ANTICANCER). Subjects with prior anticancer therapy have ~14.7% lower baseline ANC than treatment-naive subjects (Kloft 2006 / Netterberg 2017 .mod THETA(7) = -0.146837). Source NM-TRAN column `PC` (0 = no prior anticancer therapy, 1 = had prior anticancer therapy) maps directly. Time-fixed per subject.",
      source_name        = "PC"
    ),
    AAG = list(
      description        = "Serum alpha-1 acid glycoprotein concentration (acute-phase protein that binds basic and lipophilic drugs including docetaxel).",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Two distinct effects per Kloft 2006: (a) a piecewise-linear effect on baseline ANC with breakpoint at the cohort median 1.34 g/L (separate slopes e_aag_low_ba = 0.175 for AAG <= 1.34 and e_aag_high_ba = 0.495 for AAG > 1.34, applied as `BACOV *= (1 + slope * (AAG - 1.34))`); (b) a linear effect on the drug-effect slope SL (e_aag_sl = -0.351 applied as `SL *= (1 + e_aag_sl * (AAG - 1.34))`). Reference value 1.34 g/L is the Kloft 2006 cohort median. Source NM-TRAN column `AAG` (g/L) maps directly. Time-fixed per subject.",
      source_name        = "AAG"
    ),
    CP_MGL = list(
      description        = "Time-varying instantaneous docetaxel plasma concentration supplied per event row as the PD driver.",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drug-effect input: `drug = SL * CP_MGL` with SL in 1/(mg/L) units. CP_MGL is supplied as a time-varying covariate column rather than computed from a coupled PK model, reflecting the Kloft 2006 / Netterberg 2017 sequential-PK-then-PD design. The Kloft 2006 source paper feeds CP_MGL from a Bruno-style docetaxel popPK simulation (typical Cmax ~3 mg/L after a 100 mg/m^2 1-hour IV infusion, biexponential decay over ~72 hours). Set CP_MGL = 0 outside the drug-exposure window. Source NM-TRAN column `CP` (Kloft 2006 / Netterberg 2017 'predicted docetaxel concentration', mg/L) maps directly.",
      source_name        = "CP"
    )
  )

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Adult cancer patients receiving docetaxel chemotherapy. The Kloft 2006 source analysis pools data from multiple anticancer drugs (docetaxel, paclitaxel, etoposide, CPT-11, vinflunine) into a single Friberg-family myelosuppression analysis and reports drug-specific parameter sets; the DDMORE bundle for DDMODEL00000224 implements only the docetaxel parameter set, used by Netterberg 2017 as a fixed model for an ANC-prediction-methodology study (frequent-monitoring evaluation of nadir, time-to-baseline-recovery, and time-to-different-neutropenic-grade prediction).",
    dose_range     = "Intravenous docetaxel, typical 100 mg/m^2 over 1-hour infusion every 3 weeks (Kloft 2006 / Netterberg 2017 simulated trajectory). The model itself does not encode docetaxel dosing — exposure is consumed via the CP_MGL covariate column.",
    regions        = NA_character_,
    notes          = "Population demographic detail (n_subjects, age, weight, sex, race) is not reproduced in the DDMORE bundle for DDMODEL00000224 and neither the Netterberg 2017 nor the Kloft 2006 publication PDF is on disk in this worktree. The bundle's `Simulated_myelosuppression_dailyANC.csv` contains a single virtual subject (54 records, one docetaxel cycle, daily ANC monitoring over 21 days) used as a regression-style smoke test; it is not representative of the source-paper clinical cohort. The Kloft 2006 paper develops the myelosuppression model on a pooled multi-drug cancer-patient cohort; the Netterberg 2017 paper uses the docetaxel arm of that model unchanged to evaluate prediction-methodology under frequent-monitoring scenarios. See the validation vignette's Errata section for the full list of bundle-versus-publication caveats."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PD parameters - Kloft 2006 docetaxel arm.
    # Final parameter values are read from the bundle's NM-TRAN .mod
    # $THETA / $OMEGA / $SIGMA blocks (which the .mod's `; Parameter
    # estimates as according to Kloft et al., 2006` block declares as
    # publication-fixed, not initial-only). The bundle's
    # `Output_simulated_*.lst` runs the .mod with MAXEVALS=0 (model
    # evaluation only, no re-fit) and reports the same point values to
    # 3 sig figs in the FINAL PARAMETER ESTIMATE block, confirming no
    # parameter movement. The original Kloft 2006 paper is not on disk
    # in this worktree, so a side-by-side parameter-table comparison
    # against the publication was not performed; this is documented in
    # the vignette Errata.
    # ------------------------------------------------------------------
    lba <- log(5.21965);          label("Baseline ANC BA (10^9 cells/L)")                                                            # .mod $THETA(1); .lst FINAL TH 1 = 5.22
    lmt <- log(84.1791);          label("Mean transit time MT through the proliferation -> circulation chain (h)")                   # .mod $THETA(2); .lst FINAL TH 2 = 84.2

    # SL slope: Kloft 2006 reports the slope in L/umol; the .mod converts to 1/(mg/L) using the docetaxel
    # MW (807.88 g/mol, rounded to 808): `SL = THETA(3) / 808 * 1000` where THETA(3) carries L/umol units.
    # Pre-converting the value here keeps the model body free of the MW factor and lets the input
    # CP_MGL covariate be in mg/L directly (the bundle's NM-TRAN data file ships CP in mg/L).
    lsl <- log(15.5711 / 808 * 1000); label("Linear drug-effect slope SL on docetaxel concentration (1/(mg/L))")                     # .mod $THETA(3) = 15.5711 L/umol; pre-converted via docetaxel MW 808 -> 19.27 1/(mg/L)
    lpo <- log(0.144543);         label("Feedback exponent PO on the (BA / circ) feedback term (unitless)")                          # .mod $THETA(4); .lst FINAL TH 4 = 0.145

    # ------------------------------------------------------------------
    # Covariate effects - Kloft 2006 docetaxel arm.
    # All effects enter as multiplicative factors `(1 + theta * <indicator>)`
    # on either BACOV (baseline ANC) or SLCOV (drug-effect slope). The
    # piecewise-linear AAG effect on BACOV uses two separate slopes around
    # the cohort median 1.34 g/L.
    # ------------------------------------------------------------------
    e_sexf_ba       <- -0.121451;  label("Multiplicative effect of female sex on baseline ANC (unitless)")                          # .mod $THETA(11); .lst FINAL TH11 = -0.121
    e_ecogge1_ba    <-  0.130406;  label("Multiplicative effect of ECOG performance status >= 1 on baseline ANC (unitless)")        # .mod $THETA(8);  .lst FINAL TH 8 =  0.130
    e_pc_ba         <- -0.146837;  label("Multiplicative effect of prior anticancer therapy on baseline ANC (unitless)")            # .mod $THETA(7);  .lst FINAL TH 7 = -0.147
    e_aag_low_ba    <-  0.174677;  label("Piecewise-linear AAG slope on baseline ANC for AAG <= 1.34 g/L (1/(g/L))")                 # .mod $THETA(9);  .lst FINAL TH 9 =  0.175
    e_aag_high_ba   <-  0.494618;  label("Piecewise-linear AAG slope on baseline ANC for AAG  > 1.34 g/L (1/(g/L))")                 # .mod $THETA(10); .lst FINAL TH10 =  0.495
    e_aag_sl        <- -0.350693;  label("Linear AAG slope on the drug-effect slope SL, centered at AAG = 1.34 g/L (1/(g/L))")       # .mod $THETA(6);  .lst FINAL TH 6 = -0.351

    # ------------------------------------------------------------------
    # Inter-individual variability - Kloft 2006 docetaxel arm.
    # The .mod $OMEGA block is diagonal (3 ETAs on BA, MT, SL); no inter-
    # parameter correlations are reported. ETA on PO is fixed to 0 in the
    # bundle (Model_Accommodations.txt: "the OMEGA related to the gamma
    # parameter was set to 0 here (estimated to 0.0216452 in the original
    # publication)"). This deviation is documented in vignette Errata.
    # ------------------------------------------------------------------
    etalba ~ 0.0639703   # .mod $OMEGA(1,1); .lst FINAL ETA1 = 6.40e-2 (variance on log-eta scale)
    etalmt ~ 0.0191785   # .mod $OMEGA(2,2); .lst FINAL ETA2 = 1.92e-2
    etalsl ~ 0.128412    # .mod $OMEGA(3,3); .lst FINAL ETA3 = 1.28e-1

    # ------------------------------------------------------------------
    # Residual error - Kloft 2006 docetaxel arm.
    # The .mod $ERROR block fits log-additive residual error: `Y = LOG(F) + W * EPS(1)`
    # with `W = THETA(5)` and `$SIGMA 1 FIX`. This is exactly nlmixr2's lnorm()
    # residual error: ANC ~ lnorm(addSd) with addSd = W on the log scale.
    # ------------------------------------------------------------------
    addSd <- 0.424093;    label("Lognormal residual SD on log(ANC) (unitless on log scale)")                                          # .mod $THETA(5); .lst FINAL TH 5 = 0.424
  })

  model({
    # ------------------------------------------------------------------
    # Covariate-modified individual parameters
    # ------------------------------------------------------------------

    # Composite baseline-ANC covariate factor - Kloft 2006 .mod $PK BACOV
    # block: BACOV = (1 + BASEX) * (1 + BAPERF) * (1 + BAPC) * (1 + BAAAG)
    aag_dev <- AAG - 1.34
    ba_aag_eff <- ifelse(aag_dev <= 0, e_aag_low_ba * aag_dev, e_aag_high_ba * aag_dev)
    ba_cov <-
      (1 + e_sexf_ba    * SEXF) *
      (1 + e_ecogge1_ba * ECOG_GE1) *
      (1 + e_pc_ba      * PRIOR_ANTICANCER) *
      (1 + ba_aag_eff)

    # Drug-effect-slope covariate factor - Kloft 2006 .mod $PK SLCOV block
    sl_cov <- 1 + e_aag_sl * aag_dev

    # Individual parameters (PD layer - the only layer in this PD-only model)
    ba <- exp(lba + etalba) * ba_cov
    mt <- exp(lmt + etalmt)
    sl <- exp(lsl + etalsl) * sl_cov
    po <- exp(lpo)

    # Transit-rate constant: the chain has 3 transit compartments plus the
    # proliferation and circulating compartments, so KTR = (NN + 1) / MTT = 4 / MTT
    # per the .mod $PK block (`K = 4/MT`).
    ktr <- 4 / mt

    # ------------------------------------------------------------------
    # Drug effect (linear)
    # ------------------------------------------------------------------
    # CP_MGL is the time-varying docetaxel plasma concentration (mg/L),
    # supplied per event row from an upstream popPK simulation. SL has
    # been pre-converted to 1/(mg/L) so the product is unitless.
    drug <- sl * CP_MGL

    # ------------------------------------------------------------------
    # Friberg myelosuppression chain (Kloft 2006 docetaxel arm).
    # Compartment ordering follows the canonical Friberg precursor1..4 + circ
    # convention used elsewhere in nlmixr2lib (see Friberg_2002_paclitaxel.R):
    #   precursor1 = self-renewing proliferating pool (Kloft .mod A(2) "STEM")
    #   precursor2 = transit 1                          (Kloft .mod A(3) "TRANSIT1")
    #   precursor3 = transit 2                          (Kloft .mod A(4) "TRANSIT2")
    #   precursor4 = transit 3                          (Kloft .mod A(5) "TRANSIT3")
    #   circ       = circulating neutrophils            (Kloft .mod A(1) "CIRC", DEFOBS)
    # The .mod's auxiliary compartments 6-12 (return-to-baseline / grade-0 /
    # grade-4 / nadir bookkeeping integrators) are intentionally not implemented
    # here - they have no biological role and are computed in the validation
    # vignette as post-hoc summaries of the simulated ANC trajectory.
    # ------------------------------------------------------------------
    d/dt(circ)       <-  ktr * precursor4 - ktr * circ
    d/dt(precursor1) <-  ktr * precursor1 * (1 - drug) * (ba / circ)^po - ktr * precursor1
    d/dt(precursor2) <-  ktr * precursor1 - ktr * precursor2
    d/dt(precursor3) <-  ktr * precursor2 - ktr * precursor3
    d/dt(precursor4) <-  ktr * precursor3 - ktr * precursor4

    # Initial conditions: all five myelosuppression compartments start at the
    # individual covariate-modified baseline, reproducing the .mod's bioavailability
    # trick (`F1=BA; F2=BA; ...; F5=BA` combined with TIME=0 AMT=1 records on
    # CMT=1..5 in the bundle's NM-TRAN dataset, which initializes each compartment
    # to BA via F_i * 1 = BA at TIME=0).
    circ(0)       <- ba
    precursor1(0) <- ba
    precursor2(0) <- ba
    precursor3(0) <- ba
    precursor4(0) <- ba

    # ------------------------------------------------------------------
    # Observation: absolute neutrophil count ANC (10^9 cells/L) on linear scale.
    # The .mod fits log-additive residual error (`Y = LOG(F) + W*EPS(1)`,
    # SIGMA fixed at 1, W = THETA(5)); this is exactly nlmixr2's lnorm()
    # residual where the `addSd` argument is the SD on the log scale.
    # ------------------------------------------------------------------
    ANC <- circ
    ANC ~ lnorm(addSd)
  })
}
