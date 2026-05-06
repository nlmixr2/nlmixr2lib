Hansson_2013a_sunitinib <- function() {
  description <- "Population PD biomarker model for sunitinib in adults with imatinib-resistant gastrointestinal stromal tumours (GIST). Four indirect-response compartments for the soluble biomarkers VEGF, sVEGFR-2, sVEGFR-3, and sKIT, each driven by a per-cycle drug-exposure summary AUC = DOSE_MG / CLI. Sigmoid Imax inhibition with a Hill coefficient applies to VEGF (Kout) and sVEGFR-2 (Kin); simple Imax inhibition applies to sVEGFR-3 (Kin) and sKIT (Kin). A linear disease-progression term increases the baseline of VEGF and sKIT over time. The PD model has no PK ODE: the user supplies DOSE_MG (current daily sunitinib dose, mg, time-varying with on/off cycling) and CLI (subject-specific posthoc total plasma clearance, L/h, from an upstream popPK fit) as data columns. No covariates other than the two exposure inputs."
  reference <- paste(
    "Hansson EK, Amantea MA, Westwood P, Milligan PA, Houk BE,",
    "French J, Karlsson MO, Friberg LE.",
    "PKPD modeling of VEGF, sVEGFR-2, sVEGFR-3, and sKIT as predictors of",
    "tumor dynamics and overall survival following sunitinib treatment in GIST.",
    "CPT Pharmacometrics Syst Pharmacol. 2013;2(11):e84.",
    "doi:10.1038/psp.2013.61.",
    "DDMORE Foundation Model Repository: DDMODEL00000197.",
    sep = " "
  )
  vignette <- "Hansson_2013a_sunitinib"
  units <- list(time = "hour", dosing = "mg", concentration = "pg/mL")
  ddmore_id <- "DDMODEL00000197"
  replicate_of <- NULL

  covariateData <- list(
    DOSE_MG = list(
      description        = "Current administered sunitinib daily dose (mg) carried as a time-varying data column. Set to 0 during off-cycles (4 weeks on / 2 weeks off in the Hansson 2013a GIST cohort) or for placebo subjects so the derived AUC = DOSE_MG / CLI becomes 0.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "DDMORE-bundle simulated dataset reports DOSE_MG = 50 (treated) or 0 (placebo) for every record of every subject. The .mod feeds DOSE_MG into AUC = DOSE_MG / CLI in $PK at every event call, producing a per-cycle daily-AUC equivalent (mg*h/L). For typical-cohort vignette simulations the value is held at 50 mg during the 4-week on-cycles.",
      source_name        = "DOS"
    ),
    CLI = list(
      description        = "Individual posthoc total plasma clearance (L/h) of sunitinib from the paper's upstream 2-compartment popPK fit. Per-subject, time-fixed.",
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The DDMORE-bundle simulated dataset carries CLI = 32.819 L/h for subject 1; this value (which is broadly consistent with the typical sunitinib CL reported by Houk et al. 2010) is used as the typical-value reference for the validation vignette's virtual cohort. For a re-fit or new-population simulation the user must supply individual CL drawn from a sunitinib popPK model (the Hansson 2013 paper text describes the upstream PK as a 'previously developed 2-compartment model'; that popPK is not extracted into nlmixr2lib).",
      source_name        = "CL"
    )
  )

  population <- list(
    n_subjects     = 303L,
    n_studies      = 1L,
    age_range      = "adults with imatinib-resistant GIST (paper not on disk; Hansson 2013 main-text demographics not available in the bundle)",
    weight_range   = "not reported in the DDMORE bundle",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Imatinib-resistant gastrointestinal stromal tumours (GIST). Phase III trial included a placebo-controlled run-in (PLA = 1) and an active sunitinib arm (PLA = 0).",
    dose_range     = "Sunitinib 50 mg PO QD on a 4-weeks-on / 2-weeks-off schedule (standard GIST regimen at the time of the source study). Placebo arm: no sunitinib.",
    regions        = "Phase III multinational trial; specific regions not reported in the DDMORE bundle.",
    biomarkers     = "VEGF, sVEGFR-2 (soluble VEGF receptor 2), sVEGFR-3 (soluble VEGF receptor 3), sKIT (soluble KIT receptor); plasma concentrations measured serially across treatment cycles. Source dataset reports concentrations on linear (DVX) and log-transformed (DV) scales; the model fits log(observation).",
    notes          = "n_subjects = 303 and n_observations = 5394 are reported in the .lst header (`TOT. NO. OF INDIVIDUALS: 303`, `TOT. NO. OF OBS RECS: 5394`). The detailed baseline-demographics table (age, weight, sex, race, prior-imatinib-duration distributions) is in the linked publication (CPT Pharmacometrics Syst Pharmacol 2013;2:e84), which is not on disk; populate the missing population fields if the paper PDF becomes available. The bundle's simulated dataset is intentionally minimal (single subject) and is not representative of the published cohort."
  )

  ini({
    # ----------------------------------------------------------------------
    # Final estimates from Output_real_Biomarker_GIST.lst FINAL PARAMETER
    # ESTIMATE block (the .lst is an evaluation run with MAXEVALS=0 POSTHOC,
    # so .mod $THETA initial values equal the published final estimates).
    # ----------------------------------------------------------------------

    # ----- VEGF baseline + MRT + DP slope -----
    lbl_vegf      <- log(59.7);   label("VEGF baseline (pg/mL)")                                # Output_real .lst TH 1
    lmrt_vegf     <- log(91);     label("VEGF mean residence time (h)")                         # Output_real .lst TH 2
    hill_vegf     <- 3.31;        label("Hill coefficient for VEGF sigmoid Imax (unitless)")    # Output_real .lst TH 5

    # ----- sVEGFR-2 baseline + MRT + Hill -----
    lbl_svegfr2   <- log(8670);   label("sVEGFR-2 baseline (pg/mL)")                            # Output_real .lst TH 7
    lmrt_svegfr2  <- log(554);    label("sVEGFR-2 mean residence time (h)")                     # Output_real .lst TH 8
    hill_svegfr2  <- 1.54;        label("Hill coefficient for sVEGFR-2 sigmoid Imax (unitless)") # Output_real .lst TH 9

    # ----- sVEGFR-3 baseline + MRT (no Hill — simple Imax) -----
    lbl_svegfr3   <- log(63900);  label("sVEGFR-3 baseline (pg/mL)")                            # Output_real .lst TH10
    lmrt_svegfr3  <- log(401);    label("sVEGFR-3 mean residence time (h)")                     # Output_real .lst TH11

    # ----- sKIT baseline + MRT (no Hill — simple Imax) -----
    lbl_skit      <- log(39200);  label("sKIT baseline (pg/mL)")                                # Output_real .lst TH12
    lmrt_skit     <- log(2430);   label("sKIT mean residence time (h)")                         # Output_real .lst TH13

    # ----- Common IC50 (typical value shared across all 4 biomarkers) -----
    # The .mod uses THETA(4) as a single typical IC50 for VEGF, sVEGFR-2,
    # sVEGFR-3 and sKIT; per-biomarker variation enters via the BLOCK(4)
    # IIV on (etalic50_vegf, etalic50_svegfr2, etalic50_svegfr3, etalic50_skit).
    lic50         <- log(1.0);    label("Typical IC50 for biomarker drug effect (mg*h/L AUC)") # Output_real .lst TH 4

    # ----- Linear disease-progression slope (shared typical between VEGF and sKIT) -----
    # The .mod parameterizes TVSLO = THETA(6)/1000 and applies the same
    # typical value to both VEGF (DPSLO = TVSLO * exp(eta11)) and sKIT
    # (DPSLOS = TVSLO * exp(eta12)). Final TH 6 = 0.035 hence typical slope
    # = 3.5e-5 per hour; log-transform is safe at the typical value (positive),
    # though the .mod $THETA bound of (-0.06, 0.035) admits negative slopes.
    ldp_slope     <- log(3.5e-5); label("Typical linear disease-progression slope on VEGF and sKIT baselines (1/h)") # Output_real .lst TH 6 / 1000

    # ----------------------------------------------------------------------
    # Inter-individual variability — Output_real .lst FINAL OMEGA block.
    # All omega^2 values on the internal log scale; CV% conversion is the
    # standard log-normal: CV% = sqrt(exp(omega^2) - 1).
    # ----------------------------------------------------------------------

    # Diagonal IIV on the four biomarker baselines.
    etalbl_vegf     ~ 0.252;    label("IIV variance on log VEGF baseline (~57% CV)")        # Output_real .lst OMEGA(1,1)
    etalbl_svegfr2  ~ 0.0369;   label("IIV variance on log sVEGFR-2 baseline (~19% CV)")    # Output_real .lst OMEGA(2,2)
    etalbl_svegfr3  ~ 0.186;    label("IIV variance on log sVEGFR-3 baseline (~46% CV)")    # Output_real .lst OMEGA(3,3)
    etalbl_skit     ~ 0.254;    label("IIV variance on log sKIT baseline (~57% CV)")        # Output_real .lst OMEGA(4,4)

    # Shared MRT eta on VEGF, sVEGFR-2, sVEGFR-3 (the .mod assigns ETA(5) to
    # MRT, MRT2, and MRT3 — a single eta couples the three biomarkers' MRTs);
    # separate eta on sKIT MRT.
    etalmrt_pooled  ~ 0.0600;   label("IIV variance on log MRT shared by VEGF, sVEGFR-2, sVEGFR-3 (~25% CV)") # Output_real .lst OMEGA(5,5)
    etalmrt_skit    ~ 0.0753;   label("IIV variance on log sKIT MRT (~28% CV)")                              # Output_real .lst OMEGA(6,6)

    # BLOCK(4) IIV on the four biomarkers' IC50 — high inter-biomarker
    # correlation (per the FINAL OMEGA correlation matrix in the .lst,
    # off-diagonal correlations 0.62 to 0.92). Lower-triangle entries
    # in NONMEM order: (1,1) (2,1)(2,2) (3,1)(3,2)(3,3) (4,1)(4,2)(4,3)(4,4).
    etalic50_vegf + etalic50_svegfr2 + etalic50_svegfr3 + etalic50_skit ~
      c(0.253,
        0.198, 0.189,
        0.238, 0.252, 0.398,
        0.218, 0.297, 0.936, 5.77)
    # Output_real .lst OMEGA BLOCK(4) rows 7–10

    # Independent IIV on the disease-progression slope (separate per biomarker).
    etaldp_slope_vegf ~ 2.95;   label("IIV variance on log VEGF disease-progression slope (large; see vignette)")   # Output_real .lst OMEGA(11,11)
    etaldp_slope_skit ~ 3.01;   label("IIV variance on log sKIT disease-progression slope (large; see vignette)")   # Output_real .lst OMEGA(12,12)

    # ----------------------------------------------------------------------
    # Residual error — Output_real .lst FINAL THETA block, slots 14–18.
    # The .mod uses log-transform-both-sides residuals (Y = LOG(A) + W*EPS),
    # which maps to proportional residual error on the linear scale in
    # nlmixr2 (per references/naming-conventions.md NONMEM->nlmixr2 table).
    # ----------------------------------------------------------------------
    propSd_vegf       <- 0.445;  label("Proportional residual error on VEGF (fraction)")                       # Output_real .lst TH14
    propSd_svegfr2    <- 0.12;   label("Proportional residual error on sVEGFR-2 (fraction)")                   # Output_real .lst TH15
    # The .mod's W = sqrt(THETA(15)^2 + (THETA(16)/A(2))^2) is mathematically
    # equivalent to nlmixr2 `prop(p) + add(a)` with p = THETA(15), a = THETA(16);
    # see vignette Source trace for the algebraic identity.
    addSd_svegfr2     <- 583;    label("Additive residual error on sVEGFR-2 (pg/mL)")                          # Output_real .lst TH16
    propSd_svegfr3    <- 0.22;   label("Proportional residual error on sVEGFR-3 (fraction)")                   # Output_real .lst TH17
    propSd_skit       <- 0.224;  label("Proportional residual error on sKIT (fraction)")                       # Output_real .lst TH18
  })

  model({
    # IMAX is fixed at 1 in the .mod (THETA(3) = 1 FIX) and shared across all
    # four biomarkers. Hardcoded here rather than carried through ini().
    imax <- 1.0

    # Typical-value parameters with IIV
    bl_vegf      <- exp(lbl_vegf      + etalbl_vegf)
    bl_svegfr2   <- exp(lbl_svegfr2   + etalbl_svegfr2)
    bl_svegfr3   <- exp(lbl_svegfr3   + etalbl_svegfr3)
    bl_skit      <- exp(lbl_skit      + etalbl_skit)

    mrt_vegf     <- exp(lmrt_vegf     + etalmrt_pooled)   # ETA(5) shared
    mrt_svegfr2  <- exp(lmrt_svegfr2  + etalmrt_pooled)
    mrt_svegfr3  <- exp(lmrt_svegfr3  + etalmrt_pooled)
    mrt_skit     <- exp(lmrt_skit     + etalmrt_skit)     # ETA(6) separate

    ic50_vegf    <- exp(lic50 + etalic50_vegf)
    ic50_svegfr2 <- exp(lic50 + etalic50_svegfr2)
    ic50_svegfr3 <- exp(lic50 + etalic50_svegfr3)
    ic50_skit    <- exp(lic50 + etalic50_skit)

    dp_slope_vegf <- exp(ldp_slope + etaldp_slope_vegf)
    dp_slope_skit <- exp(ldp_slope + etaldp_slope_skit)

    # Per-cycle drug-exposure summary fed by the two data covariates.
    # mg / (L/h) = mg*h/L; this is the daily-AUC equivalent for the current
    # dose level. CLI is per-subject (time-fixed, posthoc upstream PK);
    # DOSE_MG is the current daily dose (time-varying with on/off cycling).
    auc <- DOSE_MG / CLI

    # Drug-effect functions: sigmoid Imax (with Hill) for VEGF and sVEGFR-2,
    # plain Imax (Hill = 1) for sVEGFR-3 and sKIT.
    eff_vegf    <- imax * auc^hill_vegf    / (ic50_vegf^hill_vegf       + auc^hill_vegf)
    eff_svegfr2 <- imax * auc^hill_svegfr2 / (ic50_svegfr2^hill_svegfr2 + auc^hill_svegfr2)
    eff_svegfr3 <- imax * auc / (ic50_svegfr3 + auc)
    eff_skit    <- imax * auc / (ic50_skit    + auc)

    # Linear disease progression on VEGF and sKIT baselines.
    # rxode2 evaluates the model at every integration step, so dp_* updates
    # continuously in t (whereas the NONMEM .mod re-evaluates DP1/DPS only
    # at $PK calls). At observation-dense schedules the two are equivalent;
    # for the typical-cohort vignette this is exact.
    dp_vegf <- bl_vegf * (1 + dp_slope_vegf * t)
    dp_skit <- bl_skit * (1 + dp_slope_skit * t)

    # Indirect-response rate constants.
    kout_vegf    <- 1 / mrt_vegf
    kout_svegfr2 <- 1 / mrt_svegfr2
    kout_svegfr3 <- 1 / mrt_svegfr3
    kout_skit    <- 1 / mrt_skit

    # Production rates: VEGF and sKIT carry disease-progression-modulated
    # baselines; sVEGFR-2 and sVEGFR-3 carry constant baselines.
    kin_vegf    <- dp_vegf      * kout_vegf
    kin_svegfr2 <- bl_svegfr2   * kout_svegfr2
    kin_svegfr3 <- bl_svegfr3   * kout_svegfr3
    kin_skit    <- dp_skit      * kout_skit

    # Initial conditions at t = 0: each biomarker at its individual baseline.
    vegf(0)    <- bl_vegf
    svegfr2(0) <- bl_svegfr2
    svegfr3(0) <- bl_svegfr3
    skit(0)    <- bl_skit

    # ODEs: drug effect inhibits Kout for VEGF (build-up via reduced removal)
    # and inhibits Kin for sVEGFR-2, sVEGFR-3, and sKIT (depletion via
    # reduced production), per the .mod $DES block.
    d/dt(vegf)    <- kin_vegf                       - kout_vegf    * (1 - eff_vegf)    * vegf
    d/dt(svegfr2) <- kin_svegfr2 * (1 - eff_svegfr2) - kout_svegfr2 * svegfr2
    d/dt(svegfr3) <- kin_svegfr3 * (1 - eff_svegfr3) - kout_svegfr3 * svegfr3
    d/dt(skit)    <- kin_skit    * (1 - eff_skit)    - kout_skit    * skit

    # Residual error — log-transform-both-sides on the .mod side translates
    # to proportional on the linear scale in nlmixr2; sVEGFR-2 has an
    # additional additive component (see ini() comment for the algebra).
    vegf    ~ prop(propSd_vegf)
    svegfr2 ~ add(addSd_svegfr2) + prop(propSd_svegfr2)
    svegfr3 ~ prop(propSd_svegfr3)
    skit    ~ prop(propSd_skit)
  })
}
