Wang_2018_daclatasvir_asunaprevir <- function() {
  description <- "MBMA PK + mechanistic HCV viral-dynamic (VD) model for daclatasvir (DCV, NS5A inhibitor) and asunaprevir (ASV, NS3/4A protease inhibitor) combination therapy in adults with genotype-1 chronic hepatitis C (Wang 2018). PK was developed by a model-based meta-analysis of arm-mean concentration data pooled from 26 trials (1067 concentration records, DCV 198 subjects in 30 arms / 7 trials, ASV 290 subjects in 35 arms / 11 trials); each drug uses an independent 2-compartment disposition model with inter-arm variability (IAV) encoded as study-arm-level etas (not between-subject variability). DCV absorption is first-order; ASV absorption is simultaneous zero- plus first-order with formulation-dependent fraction FK absorbed via the zero-order route (FK=0.184 for capsule/tablet, 0.334 for suspension/solution). The shared viral dynamics is a Neumann-style three-state target-cell model (uninfected target cells `target` T, productively infected cells `infected` I, free virions `virus` V) with most system constants (Tmax, d, R0, delta) FIXED to literature values from Neumann et al 1998; virion clearance c and production p are estimated. Each drug acts via its own effect compartment with a sigmoid-Emax inhibition of virion production (Emax=1) and an empirical exponentially time-increasing IC50 capturing the emergence of drug-resistant variants (Kr coefficient). Genotype subtype (GT1A vs GT1B) modifies IC50 by a fixed scaling factor (SCL_DCV=0.18, SCL_ASV=0.30, both GT1B/GT1A ratio) and modifies the DCV resistance rate (Kr_DCV=0.43 /day for GT1A vs 0.13 /day for GT1B). Combination efficacy follows the Bliss-additive form ECOMB/(1-ECOMB) = EDCV/(1-EDCV) + EASV/(1-EASV) (Eq 13). The model is intended for simulating arm-mean PK and viral-load trajectories under DCV monotherapy, ASV monotherapy, or DCV+ASV combination regimens; downstream NCA-style summaries (Cmax, Tmax, AUC) reproduce the published per-dose-group PK profiles, and viral-load trajectories reproduce the published biphasic decline and resistance-driven rebound shapes."

  reference <- paste(
    "Wang HC, Ren YP, Qiu Y, Zheng J, Li GL, Hu CP, Zhou TY, Lu W, Li L.",
    "(2018). Integrated pharmacokinetic/viral dynamic model for",
    "daclatasvir/asunaprevir in treatment of patients with genotype 1",
    "chronic hepatitis C. Acta Pharmacologica Sinica 39(1):140-153.",
    "doi:10.1038/aps.2017.84.",
    sep = " "
  )
  vignette <- "Wang_2018_daclatasvir_asunaprevir"

  # The PK model is reported in /h units (Table 3) and the VD model in /day
  # units (Table 4). The packaged model converts the PK rate constants to
  # /day inside model() (multiply L/h by 24 to get L/day; divide h-units
  # of D_ASV by 24 to get day-units) so the integrated PK + VD system is
  # time-consistent on a single day scale. The `units$time = "day"` field
  # therefore declares the time axis the consumer should use.
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/L"
  )

  covariateData <- list(
    HCV_GT1B = list(
      description        = "HCV genotype-1 subtype indicator. 1 = patient infected with HCV genotype 1B; 0 = patient infected with HCV genotype 1A (the source-paper reference subtype for the IC50 estimates).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HCV GT1A; 77 percent of the Wang 2018 VD-cohort, 55 of 72 patients).",
      notes              = "Time-fixed per subject (HCV subtype is determined at the time of infection and does not change over the modelled treatment window). Switches IC50,DCV from 0.041 to 0.0074 ug/L via the fixed scaling factor SCL_IC50_DCV = 0.18; switches IC50,ASV from 2.45 to 0.74 ug/L via SCL_IC50_ASV = 0.30; switches the DCV resistance coefficient Kr_DCV from 0.43 to 0.13 per day. Kr_ASV is the same for both subtypes. Encoding inside model(): ic50_dcv_t0 = exp(lic50_dcv_gt1a + etalic50_dcv) * scl_ic50_dcv^HCV_GT1B; kr_dcv = exp(lkr_dcv_gt1a + etalkr_dcv) * (1 - HCV_GT1B) + exp(lkr_dcv_gt1b + etalkr_dcv) * HCV_GT1B.",
      source_name        = "Genotype 1A (%) in Table 2; the indicator encodes the complement (GT1B = 1 - GT1A) so the canonical reference category matches the source-paper IC50,GT1A estimates."
    ),
    FORM_ASV_LIQUID = list(
      description        = "Asunaprevir formulation indicator. 1 = ASV given as a suspension or oral solution (high-fraction zero-order absorption route, FK = 0.334); 0 = ASV given as a capsule or tablet (low-fraction zero-order absorption, FK = 0.184). The covariate has no effect when no ASV dose is administered.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (capsule or tablet; the reference formulation in the Wang 2018 ASV PK fit).",
      notes              = "Per-dose-occasion covariate in principle (a participant could in principle receive both formulations across study occasions), but in the Wang 2018 trials each subject received a single ASV formulation. Switches the structural absorption fraction-via-zero-order route FK between fk_cap_asv (0.184) and fk_sol_asv (0.334). Both values were estimated with a SHARED IAV CV of 65.0 percent (Table 3); the shared variance is encoded as a single eta `eta_study_lfk_asv` on the logit of FK so the same study-arm random effect applies regardless of formulation. Reference: Wang 2018 Table 3 footnote (FK Cap/Tab and FK Sus/Sol IAV both 65.0 percent).",
      source_name        = "Formulation column in Table 1 (values 'Suspension' / 'Solution' map to FORM_ASV_LIQUID = 1; 'Capsule' / 'Tablet' map to 0). Registered as a specific-scope canonical in inst/references/covariate-columns.md alongside this extraction."
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 1730L,
    n_studies       = 26L,
    n_arms          = 72L,
    pk_dcv_subjects = 198L,
    pk_dcv_arms     = 30L,
    pk_dcv_trials   = 7L,
    pk_dcv_records  = 465L,
    pk_asv_subjects = 290L,
    pk_asv_arms     = 35L,
    pk_asv_trials   = 11L,
    pk_asv_records  = 602L,
    vd_subjects     = 72L,
    vd_trials       = 4L,
    vd_records      = 952L,
    age_range       = "20-83 years across the meta-database (per-arm medians in Wang 2018 Table 2 ranged 30 to 65 years; range 19 to 83)",
    weight_range    = "BMI 19-35 kg/m^2 across the meta-database (Wang 2018 Table 2; body weight not tabulated separately and not used as a covariate in the final model)",
    sex_female_pct  = NA_real_,
    race_ethnicity  = "Caucasian arm percentages ranged 0 to 90 percent across the included trials (Wang 2018 Table 2); Japanese / Asian sub-cohorts in studies AI444007, AI447005, AI447015, AI447028, AI447031, AI447036. Race not retained as a covariate in the final model.",
    disease_state   = "Adults with chronic HCV genotype-1 infection (77 percent GT1A, 23 percent GT1B in the VD cohort) plus healthy volunteer cohorts from the Phase 1 PK studies. The VD-modelling cohort was 82.6 percent treatment-naive for prior peg-interferon-alpha plus ribavirin (PR) therapy. Median baseline viral load 6.76 x 10^6 IU/mL.",
    dose_range      = "DCV: single-ascending and multiple-ascending dose ranges 1-200 mg in Phase 1; 30-60 mg QD in Phase 2/3. ASV: 10-1200 mg single dose, 10-600 mg BID in Phase 1, 100-600 mg BID in Phase 2/3. Combination Phase 3 regimens: DCV 60 mg QD plus ASV 100-200 mg BID (capsule), or DCV 60 mg QD plus ASV 200-600 mg BID (tablet).",
    regimens        = "Once-daily for DCV; twice-daily (BID) for ASV in nearly all Phase 2/3 dual-therapy regimens. Treatment durations 1 day (single-dose), 14 days (multiple-dose Phase 1), to 24 weeks (Phase 3 combination therapy).",
    regions         = "Global; trials included North American, European, Japanese, and Asian centres (Wang 2018 Table 1 indicates Japanese ethnicity sub-cohorts in 4 studies and Asian sub-cohort in 1 Phase 3 study).",
    notes           = "Population sourced from a model-based meta-analysis (MBMA) of 26 published or registered clinical trials covering DCV monotherapy, ASV monotherapy, and DCV+ASV combination therapy. The packaged model is parameterised at the study-arm level (inter-arm variability replaces between-subject variability for PK); the VD layer has genuine inter-individual variability because the underlying viral-load fits used individual patient data from 4 trials (AI444002, AI444004, AI447002, AI447004). Suitable simulation scope: study-arm-mean PK profiles and individual viral-load trajectories. Per Wang 2018 Materials and methods, all viral-load measurements were either supplied directly or back-derived by adding the per-arm mean baseline viral load to a published change-from-baseline value. See Wang 2018 Table 1 for the per-trial breakdown and Table 2 for baseline demographics."
  )

  ini({
    # =========================================================================
    # DCV PK (Wang 2018 Table 3 left columns; rates reported in /h, converted
    # to /day inside model() so the integrated PK + VD system is time-
    # consistent on a single day scale). All structural parameters log-
    # transformed; IAV is encoded as eta_study_* per the SKILL Step 3a MBMA
    # guidance to mark these as study-arm-level random effects, not subject-
    # level. IAV is reported as CV percent in Table 3; the eta variance is the
    # log-normal conversion omega^2 = log(CV^2 + 1).
    # =========================================================================

    # ---- Structural typical values ----
    lka     <- log(1.17);  label("DCV first-order absorption rate constant (1/h; converted to /day in model() by multiplying by 24)")  # Wang 2018 Table 3 DCV: Ka = 1.17 /h
    lcl     <- log(5.24);  label("DCV elimination clearance (L/h; converted to L/day in model() by multiplying by 24)")              # Wang 2018 Table 3 DCV: CL = 5.24 L/h
    lvc     <- log(42.9);  label("DCV central compartment volume (L)")                                                                # Wang 2018 Table 3 DCV: Vc = 42.9 L
    lq      <- log(2.62);  label("DCV inter-compartmental clearance (L/h; converted to L/day in model() by multiplying by 24)")     # Wang 2018 Table 3 DCV: Q  = 2.62 L/h
    lvp     <- log(25.0);  label("DCV peripheral compartment volume (L)")                                                             # Wang 2018 Table 3 DCV: Vp = 25.0 L

    # ---- DCV inter-arm variability (study-arm-level random effects, NOT BSV) ----
    eta_study_lka  ~ log(1 + 0.283^2)  # Wang 2018 Table 3 DCV Ka IAV CV = 28.3 percent -> omega^2 = log(1 + 0.283^2) = 0.0776
    eta_study_lcl  ~ log(1 + 0.232^2)  # Wang 2018 Table 3 DCV CL IAV CV = 23.2 percent -> omega^2 = log(1 + 0.232^2) = 0.0526
    eta_study_lvc  ~ log(1 + 0.226^2)  # Wang 2018 Table 3 DCV Vc IAV CV = 22.6 percent -> omega^2 = log(1 + 0.226^2) = 0.0499
    eta_study_lvp  ~ log(1 + 0.324^2)  # Wang 2018 Table 3 DCV Vp IAV CV = 32.4 percent -> omega^2 = log(1 + 0.324^2) = 0.0998
    # Wang 2018 Table 3 DCV Q has no reported IAV (no value), so no eta_study_lq.

    # =========================================================================
    # ASV PK (Wang 2018 Table 3 right columns; rates in /h converted to /day
    # in model() except for D, which is in h and divided by 24 to get day-
    # units). Simultaneous zero- plus first-order absorption: FK fraction of
    # the dose goes via the zero-order route over duration D, (1-FK) via
    # first-order at rate Ka. The Wang 2018 manuscript contains a
    # contradiction about the FK definition (Figure 1 caption: "fraction
    # absorbed by the first-order mechanism"; Table 3 description column:
    # "Fraction of dose absorbed by the zero-order mechanism"). This packaged
    # model follows the Table 3 description (FK = zero-order fraction); the
    # operator confirmed this choice in the sidecar request-001 / response-001
    # exchange. The Figure 1 caption is treated as a transcription error and
    # the contradiction is documented in the validation vignette's Errata
    # section.
    # =========================================================================

    # ---- Structural typical values ----
    lka_asv  <- log(0.0352);  label("ASV first-order absorption rate constant (1/h; * 24 -> /day in model())")  # Wang 2018 Table 3 ASV: Ka = 0.0352 /h
    lcl_asv  <- log(432);     label("ASV elimination clearance (L/h; * 24 -> L/day in model())")               # Wang 2018 Table 3 ASV: CL = 432 L/h
    lvc_asv  <- log(1720);    label("ASV central compartment volume (L)")                                       # Wang 2018 Table 3 ASV: Vc = 1720 L
    lq_asv   <- log(237);     label("ASV inter-compartmental clearance (L/h; * 24 -> L/day in model())")        # Wang 2018 Table 3 ASV: Q  = 237 L/h
    lvp_asv  <- log(20.5);    label("ASV peripheral compartment volume (L)")                                    # Wang 2018 Table 3 ASV: Vp = 20.5 L
    ld_asv   <- log(2.58);    label("ASV zero-order absorption duration (h; / 24 -> day in model())")           # Wang 2018 Table 3 ASV: D  = 2.58 h

    # FK = fraction of dose absorbed by the zero-order mechanism (Table 3
    # description; see Errata note above). Logit-transformed because FK is
    # bounded in (0, 1). The two formulation-specific estimates share a
    # single IAV variance (Table 3 reports 65 percent CV for both, with the
    # same RSE, indicating a shared estimated random-effect on a common
    # transformation). The packaged model encodes a single eta on the
    # logit-scale FK and lets the formulation choice select between two
    # typical-value logit anchors.
    logit_fk_cap_asv <- logit(0.184)  # Wang 2018 Table 3 ASV: FK_Cap/Tab = 0.184 (zero-order fraction in capsule and tablet formulations)
    label("Logit of zero-order absorption fraction FK for ASV capsule/tablet formulation (unitless)")

    logit_fk_sol_asv <- logit(0.334)  # Wang 2018 Table 3 ASV: FK_Sus/Sol = 0.334 (zero-order fraction in suspension and solution formulations)
    label("Logit of zero-order absorption fraction FK for ASV suspension/solution formulation (unitless)")

    # ---- ASV inter-arm variability (study-arm-level random effects, NOT BSV) ----
    eta_study_lcl_asv      ~ log(1 + 0.438^2)  # Wang 2018 Table 3 ASV CL IAV CV = 43.8 percent -> omega^2 = log(1 + 0.438^2) = 0.1737
    eta_study_lvc_asv      ~ log(1 + 0.643^2)  # Wang 2018 Table 3 ASV Vc IAV CV = 64.3 percent -> omega^2 = log(1 + 0.643^2) = 0.3504
    eta_study_lvp_asv      ~ log(1 + 0.226^2)  # Wang 2018 Table 3 ASV Vp IAV CV = 22.6 percent -> omega^2 = log(1 + 0.226^2) = 0.0499
    eta_study_ld_asv       ~ log(1 + 0.392^2)  # Wang 2018 Table 3 ASV D  IAV CV = 39.2 percent -> omega^2 = log(1 + 0.392^2) = 0.1449
    eta_study_logit_fk_asv ~ 0.65^2            # Wang 2018 Table 3 ASV FK IAV CV = 65.0 percent reported on the linear FK scale; we apply that magnitude as the SD of the logit-transformed FK random effect (no closed-form CV<->omega conversion because FK is bounded). See vignette Assumptions and deviations.
    # Wang 2018 Table 3 ASV Ka has no IAV reported, no eta_study_lka_asv.
    # Wang 2018 Table 3 ASV Q  has no IAV reported, no eta_study_lq_asv.

    # =========================================================================
    # VD shared parameters (Wang 2018 Table 4). Most system constants are
    # FIXED to literature values from Neumann et al 1998 [reference 15 in
    # Wang 2018]; only virion clearance c and production rate p are
    # estimated. The basic reproductive ratio R0 is FIXED; the infection
    # rate constant beta is derived as beta = R0 * delta * c / (p * Tmax),
    # so beta inherits IIV from c and p. The target-cell production rate
    # s = d * Tmax is also derived. All rate constants are in /day.
    # =========================================================================

    # ---- FIXED literature constants ----
    lTmax    <- fixed(log(18.5e6)); label("Maximum number of hepatocytes Tmax (cells/mL); FIXED")    # Wang 2018 Table 4: Tmax = 18.5e6 cells/mL (FIX)
    ld       <- fixed(log(0.003));  label("Death rate constant of uninfected target cells d (1/day); FIXED")  # Wang 2018 Table 4: d = 0.003 /day (FIX)
    lR0      <- fixed(log(7.15));   label("Basic reproductive ratio R0 (unitless); FIXED")           # Wang 2018 Table 4: R0 = 7.15 (FIX)
    ldelta   <- fixed(log(0.139));  label("Loss rate constant of infected cells delta (1/day); FIXED")  # Wang 2018 Table 4: delta = 0.139 /day (FIX)

    # ---- Estimated VD system parameters ----
    lc       <- log(20.4); label("Virion clearance rate constant c (1/day)")                          # Wang 2018 Table 4: c = 20.4 /day
    lp       <- log(148);  label("Virion production rate constant p (virions/cells/day)")             # Wang 2018 Table 4: p = 148 virions/cells/day
    etalc    ~ log(1 + 0.221^2)  # Wang 2018 Table 4: c IIV CV = 22.1 percent -> omega^2 = log(1 + 0.221^2) = 0.0479
    etalp    ~ log(1 + 1.411^2)  # Wang 2018 Table 4: p IIV CV = 141.1 percent -> omega^2 = log(1 + 1.411^2) = 1.0972

    # ---- Effect-compartment rate constants (drug-specific) ----
    lkce_dcv <- log(0.0041)
    label("DCV effect-compartment equilibration rate constant Kce,DCV (1/day)")  # Wang 2018 Table 4: Kce,DCV = 0.0041 /day (no IIV reported)
    lkce_asv <- log(1.19)
    label("ASV effect-compartment equilibration rate constant Kce,ASV (1/day)")  # Wang 2018 Table 4: Kce,ASV = 1.19 /day
    etalkce_asv ~ log(1 + 0.432^2)  # Wang 2018 Table 4: Kce,ASV IIV CV = 43.2 percent -> omega^2 = log(1 + 0.432^2) = 0.1683

    # ---- DCV antiviral effect (Wang 2018 Table 4) ----
    # Primary IC50 estimate is for GT1A; the GT1B value is obtained via a
    # FIXED scaling factor SCL_IC50_DCV = 0.18 (= IC50_GT1B / IC50_GT1A) so
    # ic50_dcv_gt1b = 0.041 * 0.18 = 0.00738 ~ 0.0074 ug/L (paper text).
    lic50_dcv_gt1a <- log(0.041)
    label("IC50 of DCV for GT1A virion production (ug/L)")                    # Wang 2018 Table 4: IC50,DCV,GT1A = 0.041 ug/L
    scl_ic50_dcv   <- fixed(0.18)
    label("Fixed scaling factor IC50_DCV_GT1B / IC50_DCV_GT1A (unitless); FIXED")  # Wang 2018 Table 4: SCL_IC50,DCV = 0.18 (FIX)
    etalic50_dcv   ~ log(1 + 2.194^2)  # Wang 2018 Table 4: IC50,DCV IIV CV = 219.4 percent -> omega^2 = log(1 + 2.194^2) = 1.7574

    # DCV sigmoid Emax shape factor (gamma); log-transformed because positive.
    lgamma_dcv     <- log(2.25)
    label("DCV sigmoid Emax shape factor gamma_DCV (unitless)")                # Wang 2018 Table 4: gamma_DCV = 2.25
    etalgamma_dcv  ~ log(1 + 0.293^2)  # Wang 2018 Table 4: gamma_DCV IIV CV = 29.3 percent -> omega^2 = log(1 + 0.293^2) = 0.0826

    # DCV resistance coefficients - one per genotype (the source paper holds
    # the IIV CV for the two genotypes equal at 68.2 percent, so a single eta
    # is used in model() with subtype-specific typical values).
    lkr_dcv_gt1a   <- log(0.43)
    label("Coefficient for DCV IC50 exponential time-increase in GT1A (1/day)")  # Wang 2018 Table 4: Kr,DCV,GT1A = 0.43 /day
    lkr_dcv_gt1b   <- log(0.13)
    label("Coefficient for DCV IC50 exponential time-increase in GT1B (1/day)")  # Wang 2018 Table 4: Kr,DCV,GT1B = 0.13 /day
    etalkr_dcv     ~ log(1 + 0.682^2)  # Wang 2018 Table 4: Kr,DCV IIV CV = 68.2 percent (same value for GT1A and GT1B) -> omega^2 = log(1 + 0.682^2) = 0.3895

    # ---- ASV antiviral effect (Wang 2018 Table 4) ----
    lic50_asv_gt1a <- log(2.45)
    label("IC50 of ASV for GT1A virion production (ug/L)")                     # Wang 2018 Table 4: IC50,ASV,GT1A = 2.45 ug/L
    scl_ic50_asv   <- fixed(0.30)
    label("Fixed scaling factor IC50_ASV_GT1B / IC50_ASV_GT1A (unitless); FIXED")  # Wang 2018 Table 4: SCL_IC50,ASV = 0.30 (FIX)
    etalic50_asv   ~ log(1 + 0.964^2)  # Wang 2018 Table 4: IC50,ASV IIV CV = 96.4 percent -> omega^2 = log(1 + 0.964^2) = 0.6516

    # ASV sigmoid Emax shape factor and resistance coefficient have no IIV
    # reported in Wang 2018 Table 4.
    lgamma_asv     <- log(2.01)
    label("ASV sigmoid Emax shape factor gamma_ASV (unitless)")                # Wang 2018 Table 4: gamma_ASV = 2.01 (no IIV reported)
    lkr_asv        <- log(0.007)
    label("Coefficient for ASV IC50 exponential time-increase (single value for both GT1A and GT1B; 1/day)")  # Wang 2018 Table 4: Kr,ASV = 0.007 /day (single value)

    # =========================================================================
    # Residual error (Wang 2018 Table 3 PK, Table 4 VD). Table 3 reports
    # sigma^2 for the residuals; we interpret the table values as the
    # estimated VARIANCE of the underlying epsilon on the log-transformed-
    # observation scale (the model's Equation 2 carries the residual as
    # additive on ln(Cobs)). Linear-space SD therefore equals sqrt(variance)
    # and is used as the proportional-residual / additive-residual magnitude
    # in nlmixr2's ~ prop(...) + add(...) syntax. The Table 3 ASV row carries
    # both a sigma^2_prop AND a sigma^2_add component; we map them to a
    # combined error model on the linear-space ASV concentration. VD
    # residuals (Table 4) are additive on log10(viral load) per Equation 12.
    # =========================================================================
    propSd     <- sqrt(0.422)
    label("DCV proportional residual SD on plasma concentration (sqrt of sigma^2_Prop, fraction)")  # Wang 2018 Table 3 DCV: sigma^2_Prop = 0.422
    propSd_asv <- sqrt(0.495)
    label("ASV proportional residual SD on plasma concentration (sqrt of sigma^2_Prop, fraction)")  # Wang 2018 Table 3 ASV: sigma^2_Prop = 0.495
    addSd_asv  <- sqrt(0.217)
    label("ASV additive residual SD on plasma concentration (sqrt of sigma^2_Add, ug/L)")           # Wang 2018 Table 3 ASV: sigma^2_Add = 0.217

    # VD residual error - Table 4 reports separate sigma^2 values for the
    # DCV monotherapy fit (0.27) and the ASV monotherapy fit (0.29). Both
    # are additive on log10(viral load). The packaged model uses the
    # arithmetic-mean variance (0.28) for a single integrated residual,
    # equivalent to averaging the per-fit values; the disagreement is small
    # (sqrt(0.27) = 0.520 vs sqrt(0.29) = 0.539, about 4 percent of the SD)
    # and the choice is documented in the vignette Assumptions section.
    addSd_Vlog10  <- sqrt(0.28)
    label("Viral-load additive residual SD on log10 scale (averaged across DCV and ASV monotherapy fits)")  # Wang 2018 Table 4: sigma^2_DCV = 0.27, sigma^2_ASV = 0.29; averaged
  })

  model({
    # =========================================================================
    # 1. PK time-unit conversion (h -> day). The Wang 2018 PK rates and Q
    # are reported in /h; rxode2 integrates the system on the unit declared
    # in the events table (day, per units$time). Multiply rate-units of
    # /h by 24 to obtain /day; divide h-units (D, duration of zero-order
    # absorption for ASV) by 24 to obtain day. The volumes Vc and Vp have
    # no time-unit conversion (litres are independent of /h vs /day).
    # =========================================================================
    h_per_day <- 24

    # =========================================================================
    # 2. DCV individual PK parameters (with study-arm random effects).
    # =========================================================================
    ka <- exp(lka + eta_study_lka) * h_per_day
    cl <- exp(lcl + eta_study_lcl) * h_per_day
    vc <- exp(lvc + eta_study_lvc)
    q  <- exp(lq)                  * h_per_day
    vp <- exp(lvp + eta_study_lvp)

    # =========================================================================
    # 3. ASV individual PK parameters (with study-arm random effects).
    # =========================================================================
    ka_asv <- exp(lka_asv) * h_per_day
    cl_asv <- exp(lcl_asv + eta_study_lcl_asv) * h_per_day
    vc_asv <- exp(lvc_asv + eta_study_lvc_asv)
    q_asv  <- exp(lq_asv)  * h_per_day
    vp_asv <- exp(lvp_asv + eta_study_lvp_asv)
    d_asv  <- exp(ld_asv  + eta_study_ld_asv) / h_per_day

    # ASV FK (zero-order absorption fraction) - formulation-dependent via
    # FORM_ASV_LIQUID (1 = suspension/solution, 0 = capsule/tablet). The
    # logit transform keeps FK strictly in (0, 1) even with the IAV eta.
    logit_fk_asv_typ <- logit_fk_cap_asv +
                        (logit_fk_sol_asv - logit_fk_cap_asv) * FORM_ASV_LIQUID
    fk_asv <- expit(logit_fk_asv_typ + eta_study_logit_fk_asv)

    # =========================================================================
    # 4. VD shared parameters - back-transform fixed and estimated VD params
    # to linear scale (all in /day after step 1's conversion).
    # =========================================================================
    Tmax_vd   <- exp(lTmax)
    d_vd      <- exp(ld)
    R0_vd     <- exp(lR0)
    delta_vd  <- exp(ldelta)
    c_vd      <- exp(lc + etalc)
    p_vd      <- exp(lp + etalp)

    # Derived VD constants. s = d * Tmax (Wang 2018 sets s/d = Tmax = 18.5e6
    # at the pre-infection steady state). beta is the infection rate
    # constant; derived from the basic reproductive ratio R0 = beta * p *
    # Tmax / (delta * c) so beta = R0 * delta * c / (p * Tmax). beta
    # inherits the IIV from c and p.
    s_vd    <- d_vd * Tmax_vd
    beta_vd <- R0_vd * delta_vd * c_vd / (p_vd * Tmax_vd)

    # =========================================================================
    # 5. Effect-compartment rate constants (already in /day per Table 4).
    # =========================================================================
    kce_dcv <- exp(lkce_dcv)
    kce_asv <- exp(lkce_asv + etalkce_asv)

    # =========================================================================
    # 6. Drug-effect parameters (sigmoid Emax + time-varying IC50).
    # IC50 scaling by genotype: ic50_gt1b = ic50_gt1a * SCL_IC50, so the
    # general expression is ic50 = ic50_gt1a * SCL_IC50^HCV_GT1B (= ic50_gt1a
    # when HCV_GT1B=0 and = ic50_gt1a * SCL_IC50 when HCV_GT1B=1).
    # =========================================================================
    ic50_dcv_t0   <- exp(lic50_dcv_gt1a + etalic50_dcv) * scl_ic50_dcv^HCV_GT1B
    ic50_asv_t0   <- exp(lic50_asv_gt1a + etalic50_asv) * scl_ic50_asv^HCV_GT1B

    # DCV resistance coefficient varies between GT1A and GT1B (Wang 2018
    # Table 4); ASV has a single Kr for both subtypes. The single etalkr_dcv
    # is shared across genotypes per the source paper's reported common
    # IIV magnitude.
    kr_dcv_indiv  <- exp(lkr_dcv_gt1a + etalkr_dcv) * (1 - HCV_GT1B) +
                     exp(lkr_dcv_gt1b + etalkr_dcv) * HCV_GT1B
    kr_asv_indiv  <- exp(lkr_asv)

    gamma_dcv     <- exp(lgamma_dcv + etalgamma_dcv)
    gamma_asv     <- exp(lgamma_asv)

    # Time-varying IC50 (Wang 2018 Equation 11): IC50(t) = IC50(0) *
    # exp(Kr * t). `t` is the rxode2 simulation time in days (matches
    # units$time).
    ic50_dcv      <- ic50_dcv_t0 * exp(kr_dcv_indiv * t)
    ic50_asv      <- ic50_asv_t0 * exp(kr_asv_indiv * t)

    # =========================================================================
    # 7. DCV PK ODE system (depot, central, peripheral1; first-order
    # absorption from depot at rate ka, linear elimination with rate
    # cl/vc, distribution with rate q to peripheral1).
    # =========================================================================
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl + q) * central / vc +
                          q * peripheral1 / vp
    d/dt(peripheral1) <-  q * central / vc - q * peripheral1 / vp

    # =========================================================================
    # 8. ASV PK ODE system. The simultaneous zero- plus first-order
    # absorption model routes the dose into TWO compartments in parallel:
    # a first-order depot (depot_asv, rate ka_asv) receiving fraction
    # (1 - fk_asv) of the dose, and a zero-order direct-to-central infusion
    # into central_asv receiving fraction fk_asv over duration d_asv. The
    # dose-split is achieved by the f() bioavailability anchors below.
    # =========================================================================
    d/dt(depot_asv)       <- -ka_asv * depot_asv
    d/dt(central_asv)     <-  ka_asv * depot_asv -
                              (cl_asv + q_asv) * central_asv / vc_asv +
                              q_asv * peripheral1_asv / vp_asv
    d/dt(peripheral1_asv) <-  q_asv * central_asv / vc_asv -
                              q_asv * peripheral1_asv / vp_asv

    # Dose-split bioavailability. The event-table consumer dosing this
    # ASV PK must supply TWO simultaneous dose records per administered
    # dose: amt = mg_dose, cmt = "depot_asv", evid = 1 (first-order arm);
    # amt = mg_dose, cmt = "central_asv", evid = 1, rate = -2 (zero-order
    # arm with duration controlled by dur()). The f() values split the
    # nominal dose between the two routes so the total absorbed amount
    # equals 1 * mg_dose regardless of fk_asv.
    f(depot_asv)     <- 1 - fk_asv
    f(central_asv)   <- fk_asv
    dur(central_asv) <- d_asv

    # =========================================================================
    # 9. Effect compartments (Wang 2018 Equation 9). Plasma concentrations
    # (mg/L = ug/mL) are converted to ug/L by multiplying by 1000 so they
    # match the IC50 units (ug/L) in the sigmoid Emax equation. The
    # `effect` and `effect_asv` compartments hold the lagged plasma-
    # concentration surrogate that drives inhibition.
    # =========================================================================
    cc_dcv_ugL <- (central     / vc)     * 1000  # mg/L * 1000 = ug/L
    cc_asv_ugL <- (central_asv / vc_asv) * 1000

    d/dt(effect)     <- kce_dcv * (cc_dcv_ugL - effect)
    d/dt(effect_asv) <- kce_asv * (cc_asv_ugL - effect_asv)

    # =========================================================================
    # 10. Sigmoid Emax inhibition of virion production (Wang 2018 Eq 10).
    # Emax is fixed at 1 (maximum inhibition = 100 percent). The minimum
    # of effect is 0 to avoid negative-base power issues; the effect
    # compartment never goes negative for non-negative dosing histories.
    # =========================================================================
    e_dcv <- (effect^gamma_dcv) /
             ((effect^gamma_dcv) + (ic50_dcv^gamma_dcv))
    e_asv <- (effect_asv^gamma_asv) /
             ((effect_asv^gamma_asv) + (ic50_asv^gamma_asv))

    # =========================================================================
    # 11. Combination efficacy (Wang 2018 Eq 13, Bliss-additive form):
    #   ECOMB / (1 - ECOMB) = EDCV / (1 - EDCV) + EASV / (1 - EASV)
    # Let r_dcv = EDCV / (1 - EDCV) and r_asv = EASV / (1 - EASV). Then
    #   r_comb = r_dcv + r_asv,  E_COMB = r_comb / (1 + r_comb).
    # The 1e-12 floor on (1 - E) prevents singularity when an individual
    # E approaches 1 at high effect-compartment concentration.
    # =========================================================================
    r_dcv  <- e_dcv / (1 - e_dcv + 1e-12)
    r_asv  <- e_asv / (1 - e_asv + 1e-12)
    e_total <- (r_dcv + r_asv) / (1 + r_dcv + r_asv)

    # =========================================================================
    # 12. Shared viral-dynamics ODE system (Wang 2018 Eqs 3-5).
    # =========================================================================
    d/dt(target)   <- s_vd - d_vd * target - beta_vd * virus * target
    d/dt(infected) <- beta_vd * virus * target - delta_vd * infected
    d/dt(virus)    <- (1 - e_total) * p_vd * infected - c_vd * virus

    # =========================================================================
    # 13. Steady-state initial conditions (Wang 2018 Eqs 6-8). At pre-
    # treatment steady state with V > 0:
    #   T0 = delta*c / (beta*p) = Tmax / R0
    #   V0 = s*p / (delta*c) - d/beta = (R0 - 1) * d / beta
    #   I0 = (c / p) * V0
    # All three expressions inherit individual variability from c and p
    # (and beta, which is derived from c, p).
    # =========================================================================
    target(0)   <- Tmax_vd / R0_vd
    virus(0)    <- (R0_vd - 1) * d_vd / beta_vd
    infected(0) <- (c_vd / p_vd) * ((R0_vd - 1) * d_vd / beta_vd)

    # =========================================================================
    # 14. Observation outputs and residual error.
    # - Cc: DCV plasma concentration in ug/L
    # - Cc_asv: ASV plasma concentration in ug/L
    # - Vlog10: log10 of free virus concentration (log10(IU/mL)); the
    #   1e-12 floor inside log10 prevents -Inf when virus -> 0 after
    #   complete viral eradication.
    # =========================================================================
    Cc     <- cc_dcv_ugL
    Cc_asv <- cc_asv_ugL
    Vlog10 <- log10(virus + 1e-12)

    Cc     ~ prop(propSd)
    Cc_asv ~ add(addSd_asv) + prop(propSd_asv)
    Vlog10 ~ add(addSd_Vlog10)
  })
}
