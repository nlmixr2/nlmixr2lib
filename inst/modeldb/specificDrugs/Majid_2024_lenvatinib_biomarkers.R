Majid_2024_lenvatinib_biomarkers <- function() {
  description <- "Simultaneous four-biomarker population PK/PD indirect-response model for serum VEGF, Tie-2, Ang-2 and FGF-23 during lenvatinib treatment of radioiodine-refractory differentiated thyroid cancer (Majid 2024). Each biomarker is a turnover pool with its own baseline, mean residence time and Hill coefficient; a single common Emax, a single common EC50 and a single common linear disease-progression slope are shared across all four. Lenvatinib exposure enters as the steady-state daily AUC covariate AUC_LEN and acts through a sigmoid Emax function that inhibits Kout for VEGF and FGF-23 (levels rise) and inhibits Kin for Tie-2 and Ang-2 (levels fall). There is no PK ODE: AUC_LEN is supplied as a data column, as in the source sequential fit. Companion models: modellib('Majid_2024_lenvatinib') for the PK layer and modellib('Majid_2024_lenvatinib_tumor') for the tumor-growth-inhibition layer."
  reference <- paste(
    "Majid O, Hayato S, Sreerama Reddy SH, Lalovic B, Hihara T, Hoshi T,",
    "Funahashi Y, Aluri J, Takenaka O, Yasuda S, Hussein Z.",
    "Population pharmacokinetic-pharmacodynamic modeling of serum biomarkers as",
    "predictors of tumor dynamics following lenvatinib treatment in patients with",
    "radioiodine-refractory differentiated thyroid cancer (RR-DTC).",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(6):954-969.",
    "doi:10.1002/psp4.13130.",
    "Equations 1-3 and Table 2; NONMEM control stream in Supplementary Text S2.",
    "Model form follows Hansson et al. CPT Pharmacometrics Syst Pharmacol.",
    "2013;2(11):e84; see modellib('Hansson_2013a_sunitinib').",
    sep = " "
  )
  vignette <- "Majid_2024_lenvatinib"
  paper_specific_compartments <- c("tie2", "ang2", "fgf23")
  units <- list(
    time          = "h",
    dosing        = "n/a (no drug-dosing events; lenvatinib exposure enters as the AUC_LEN covariate, not via a PK ODE)",
    concentration = "ng/mL (serum biomarker concentration)"
  )

  covariateData <- list(
    AUC_LEN = list(
      description        = "Lenvatinib steady-state daily AUC at the time of the biomarker assessment, driving the sigmoid Emax drug effect on every biomarker turnover pool.",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying, step-wise: held constant between assessments and updated when the dose level changes. In the source sequential fit the value is carried as the NONMEM data column LENAUC, computed in the upstream population PK run (Text S1 $PK) as AUC = 1000 * F1 * DGRP / CL, i.e. 1000 x relative bioavailability x current daily dose (mg) divided by the individual apparent clearance (L/h). Set to 0 for placebo subjects and off-treatment periods, which makes the drug-effect term vanish exactly and leaves only the linear disease-progression model -- the behaviour the paper describes for the placebo arm. Records with LENAUC > 10000 ng*h/mL were excluded from the source fit (Text S2 $DATA IGNORE).",
      source_name        = "LENAUC"
    )
  )

  compartmentData <- list(
    vegf = list(
      analyte  = "vascular endothelial growth factor (VEGF)",
      units    = "ng/mL",
      specimen = "serum",
      verified = TRUE
    ),
    tie2 = list(
      analyte  = "soluble TEK tyrosine kinase 2 (Tie-2)",
      units    = "ng/mL",
      specimen = "serum",
      verified = TRUE
    ),
    ang2 = list(
      analyte  = "angiopoietin 2 (Ang-2)",
      units    = "ng/mL",
      specimen = "serum",
      verified = TRUE
    ),
    fgf23 = list(
      analyte  = "fibroblast growth factor 23 (FGF-23)",
      units    = "ng/mL",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 560L,
    n_studies      = 2L,
    age_range      = "not reported as a range; 300 of 560 subjects (54 percent) were < 65 years and 260 (46 percent) were >= 65 years (Table S2)",
    weight_range   = "31.0-190.0 kg",
    weight_median  = "75.1 kg",
    sex_female_pct = 48.8,
    race_ethnicity = c(
      White                          = 74.1,
      Other_Asian                    = 9.3,
      Japanese                       = 7.3,
      Missing                        = 4.8,
      Black_African_American         = 2.1,
      Others                         = 2.0,
      Chinese                        = 0.2,
      Asian_not_Japanese_or_Chinese  = 0.2
    ),
    disease_state  = "Radioiodine-refractory differentiated thyroid cancer (RR-DTC).",
    dose_range     = "Lenvatinib 18 or 24 mg orally once daily starting dose, with protocol-driven reductions to 20, 14, 10, 8 or 4 mg for grade 2-4 adverse events; 131 of 560 subjects received placebo.",
    regions        = "Multicenter: phase 3 study E7080-G000-303 and phase 2 post-marketing study E7080-G000-211.",
    n_observations = "5132 biomarker observations from 560 RR-DTC patients.",
    biomarkers     = "Serum VEGF, soluble Tie-2, Ang-2 and FGF-23 by ELISA (R&D Systems DVE00, DTE200, DANG20; Kainos CY-4000 for FGF-23; an in-house Eisai ELISA for Ang-2 in study 211). Baseline medians: VEGF 0.42, Ang-2 3.21, Tie-2 15.11, FGF-23 0.100 ng/mL (Table S2). Thyroglobulin and thyroid stimulating hormone were measured but not modelled because of extreme variability.",
    notes          = "Demographics from Majid 2024 Table S2 (biomarkers PK/PD population, N = 560; baseline FGF-23 available for 541). Sex 287 male / 273 female gives 48.8 percent female. Race percentages computed from the Table S2 overall counts over N = 560. The paper's Results text quotes 560 patients for the biomarker dataset while the abstract quotes 558 patients for the tumor dataset."
  )

  ini({
    # ---- Per-biomarker baselines. Majid 2024 Table 2 prints the unit as
    # 'ng/L', and Table S2 additionally prints Tie-2 as 'ug/mL'. Both are
    # unit-label errors: the Methods 'Model-based PK/PD simulations' section
    # states the same four baseline values explicitly as ng/mL, and the
    # Methods bioanalysis detection ranges (Ang-2 0.188-64 ng/mL, Tie-2
    # 3.10-200 ng/mL, VEGF 31-4000 pg/mL = 0.031-4 ng/mL) bracket the
    # reported baselines only on the ng/mL scale. See the vignette Errata.
    lrbase_vegf  <- log(0.370);   label("VEGF baseline BM0 (ng/mL)")                    # Majid 2024 Table 2: VEGF Baseline = 0.370 (%RSE 2.01; bootstrap 0.366 [0.346-0.386])
    lrbase_tie2  <- log(14.6);    label("Tie-2 baseline BM0 (ng/mL)")                   # Majid 2024 Table 2: Tie-2 Baseline = 14.6 (%RSE 0.829)
    lrbase_ang2  <- log(3.36);    label("Ang-2 baseline BM0 (ng/mL)")                   # Majid 2024 Table 2: Ang-2 Baseline = 3.36 (%RSE 1.55)
    lrbase_fgf23 <- log(0.0990);  label("FGF-23 baseline BM0 (ng/mL)")                  # Majid 2024 Table 2: FGF-23 Baseline = 0.0990 (%RSE 2.09)

    # ---- Per-biomarker mean residence times; Kout = 1 / MRT (Text S2 $PK).
    lmrt_vegf    <- log(58.3);    label("VEGF mean residence time MRT = 1/Kout (h)")    # Majid 2024 Table 2: VEGF MRT = 58.3 h (%RSE 30.5)
    lmrt_tie2    <- log(354);     label("Tie-2 mean residence time MRT = 1/Kout (h)")   # Majid 2024 Table 2: Tie-2 MRT = 354 h (%RSE 5.82)
    lmrt_ang2    <- log(173);     label("Ang-2 mean residence time MRT = 1/Kout (h)")   # Majid 2024 Table 2: Ang-2 MRT = 173 h (%RSE 11.5)
    lmrt_fgf23   <- log(265);     label("FGF-23 mean residence time MRT = 1/Kout (h)")  # Majid 2024 Table 2: FGF-23 MRT = 265 h (%RSE 15.5)

    # ---- Hill coefficients. Fixed to 1 for VEGF and FGF-23 (plain Emax) and
    # estimated for Tie-2 and Ang-2 (sigmoid Emax), per Majid 2024 Results
    # and Text S2 $THETA slots 10 and 13, both flagged '1 FIX'.
    hill_vegf    <- fixed(1);     label("Hill coefficient for the VEGF drug effect (unitless)")    # Majid 2024 Text S2 $THETA 10: '1 FIX ; [HILL_VEGF]'; Results: 'Hill coefficient was fixed to 1 for VEGF and FGF-23'
    hill_tie2    <- 0.313;        label("Hill coefficient for the Tie-2 drug effect (unitless)")   # Majid 2024 Table 2: Tie-2 Hill coefficient = 0.313 (%RSE 11.6)
    hill_ang2    <- 4.27;         label("Hill coefficient for the Ang-2 drug effect (unitless)")   # Majid 2024 Table 2: Ang-2 Hill coefficient = 4.27 (%RSE 16.2)
    hill_fgf23   <- fixed(1);     label("Hill coefficient for the FGF-23 drug effect (unitless)")  # Majid 2024 Text S2 $THETA 13: '1 FIX ; [HILL_FGF23]'

    # ---- Parameters common to all four biomarkers. The paper collapsed the
    # per-biomarker estimates because they were similar and the OFV did not
    # deteriorate (Results, 'PK/PD analysis for biomarkers').
    lemax        <- log(0.344);   label("Common maximum drug effect Emax (unitless fraction)")     # Majid 2024 Table 2: All common E max = 0.344 (%RSE 2.33)
    lec50        <- log(930);     label("Common EC50 on lenvatinib daily AUC (ng*h/mL)")           # Majid 2024 Table 2: All common EC50 = 930 ng*h/mL (%RSE 7.10)
    ldp_slope    <- log(2.93e-6); label("Common linear disease-progression slope (1/h)")           # Majid 2024 Table 2: Disease progression slope = 2.93e-6 1/h (%RSE 8.05); Text S2 parameterises TVSLO = THETA(15)/1000

    # ---- Inter-individual variability. As in the PK model, Majid 2024
    # Table 2's abbreviation footnote defines the printed 'CV %' as
    # 'square root of variance x 100', so the variances below are
    # (CV%/100)^2 with no log-normal back-conversion. Corroborated by the
    # Text S2 $OMEGA initial estimates (variances on the same scale):
    # baselines 0.25 vs 0.399/0.066/0.144/0.251, MRT 2.5/1/1/1 vs
    # 2.372/0.036/0.607/0.994, Emax 0.25 vs 0.214, EC50 1.25/0.5/0.25/0.5
    # vs 0.845/0.372/0.104/0.560, DP slope 5/3.3/3/5 vs 5.52/3.72/4.88/16.
    etalrbase_vegf     ~ 0.399424;  label("IIV variance on log VEGF baseline")     # Majid 2024 Table 2: VEGF Baseline 63.2 CV% -> (0.632)^2
    etalrbase_tie2     ~ 0.066049;  label("IIV variance on log Tie-2 baseline")    # Majid 2024 Table 2: Tie-2 Baseline 25.7 CV% -> (0.257)^2
    etalrbase_ang2     ~ 0.143641;  label("IIV variance on log Ang-2 baseline")    # Majid 2024 Table 2: Ang-2 Baseline 37.9 CV% -> (0.379)^2
    etalrbase_fgf23    ~ 0.251001;  label("IIV variance on log FGF-23 baseline")   # Majid 2024 Table 2: FGF-23 Baseline 50.1 CV% -> (0.501)^2

    etalmrt_vegf       ~ 2.371600;  label("IIV variance on log VEGF MRT")          # Majid 2024 Table 2: VEGF MRT 154 CV% -> (1.54)^2
    etalmrt_tie2       ~ 0.036481;  label("IIV variance on log Tie-2 MRT")         # Majid 2024 Table 2: Tie-2 MRT 19.1 CV% -> (0.191)^2
    etalmrt_ang2       ~ 0.606841;  label("IIV variance on log Ang-2 MRT")         # Majid 2024 Table 2: Ang-2 MRT 77.9 CV% -> (0.779)^2
    etalmrt_fgf23      ~ 0.994009;  label("IIV variance on log FGF-23 MRT")        # Majid 2024 Table 2: FGF-23 MRT 99.7 CV% -> (0.997)^2

    # EC50 has one common typical value but a separate eta per biomarker
    # (Text S2: IC50 = THETA(14)*EXP(ETA(14)), IC502 = THETA(14)*EXP(ETA(15)),
    # IC503 = THETA(14)*EXP(ETA(16)), IC504 = THETA(14)*EXP(ETA(17))).
    etalec50_vegf      ~ 0.844561;  label("IIV variance on log EC50, VEGF")        # Majid 2024 Table 2: VEGF EC50 91.9 CV% -> (0.919)^2
    etalec50_tie2      ~ 0.372100;  label("IIV variance on log EC50, Tie-2")       # Majid 2024 Table 2: Tie-2 EC50 61.0 CV% -> (0.610)^2
    etalec50_ang2      ~ 0.103684;  label("IIV variance on log EC50, Ang-2")       # Majid 2024 Table 2: Ang-2 EC50 32.2 CV% -> (0.322)^2
    etalec50_fgf23     ~ 0.559504;  label("IIV variance on log EC50, FGF-23")      # Majid 2024 Table 2: FGF-23 EC50 74.8 CV% -> (0.748)^2

    # The disease-progression slope likewise has one common typical value and
    # a separate eta per biomarker (Text S2 DPSLO..DPSLO4 all use TVSLO).
    etaldp_slope_vegf  ~  5.522500; label("IIV variance on log disease-progression slope, VEGF")    # Majid 2024 Table 2: VEGF DP slope 235 CV% -> (2.35)^2
    etaldp_slope_tie2  ~  3.724900; label("IIV variance on log disease-progression slope, Tie-2")   # Majid 2024 Table 2: Tie-2 DP slope 193 CV% -> (1.93)^2
    etaldp_slope_ang2  ~  4.884100; label("IIV variance on log disease-progression slope, Ang-2")   # Majid 2024 Table 2: Ang-2 DP slope 221 CV% -> (2.21)^2
    etaldp_slope_fgf23 ~ 16.000000; label("IIV variance on log disease-progression slope, FGF-23")  # Majid 2024 Table 2: FGF-23 DP slope 400 CV% -> (4.00)^2

    # Emax carries a single eta shared by all four biomarkers (Text S2 sets
    # IMAX1..IMAX4 all to THETA(9)*EXP(ETA(9))).
    etalemax           ~ 0.214369;  label("IIV variance on log Emax, common to all four biomarkers") # Majid 2024 Table 2: All common E max 46.3 CV% -> (0.463)^2

    # ---- Residual error. Text S2 fits every biomarker on the log scale as
    # Y = LOG(A) + W*EPS(1) with $SIGMA 1 FIXED, so the THETA values are
    # already standard deviations. A log-scale additive error is
    # proportional error in nlmixr2's linear space. For VEGF the source uses
    # W = sqrt(prop^2 + (add/A)^2), which multiplied through by A is
    # sqrt((prop*A)^2 + add^2) -- exactly nlmixr2's default combined
    # add + prop form. The FGF-23 additive term is fixed to zero in the
    # source (Text S2 $THETA 21: '0 FIX'), so FGF-23 is purely proportional.
    propSd_vegf  <- 0.261;   label("Proportional residual error on VEGF (fraction)")   # Majid 2024 Table 2: VEGF Proportional error 26.1 CV% (%RSE 0.877)
    addSd_vegf   <- 0.0857;  label("Additive residual error on VEGF (ng/mL)")          # Majid 2024 Table 2: VEGF Additive = 0.0857 (%RSE 1.75)
    propSd_tie2  <- 0.116;   label("Proportional residual error on Tie-2 (fraction)")  # Majid 2024 Table 2: Tie-2 Proportional error 11.6 CV% (%RSE 0.656)
    propSd_ang2  <- 0.208;   label("Proportional residual error on Ang-2 (fraction)")  # Majid 2024 Table 2: Ang-2 Proportional error 20.8 CV% (%RSE 0.625)
    propSd_fgf23 <- 0.249;   label("Proportional residual error on FGF-23 (fraction)") # Majid 2024 Table 2: FGF-23 Proportional error 24.9 CV% (%RSE 1.08)
  })

  model({
    # ---- Individual parameters.
    rbase_vegf  <- exp(lrbase_vegf  + etalrbase_vegf)
    rbase_tie2  <- exp(lrbase_tie2  + etalrbase_tie2)
    rbase_ang2  <- exp(lrbase_ang2  + etalrbase_ang2)
    rbase_fgf23 <- exp(lrbase_fgf23 + etalrbase_fgf23)

    mrt_vegf    <- exp(lmrt_vegf    + etalmrt_vegf)
    mrt_tie2    <- exp(lmrt_tie2    + etalmrt_tie2)
    mrt_ang2    <- exp(lmrt_ang2    + etalmrt_ang2)
    mrt_fgf23   <- exp(lmrt_fgf23   + etalmrt_fgf23)

    ec50_vegf   <- exp(lec50 + etalec50_vegf)
    ec50_tie2   <- exp(lec50 + etalec50_tie2)
    ec50_ang2   <- exp(lec50 + etalec50_ang2)
    ec50_fgf23  <- exp(lec50 + etalec50_fgf23)

    dp_slope_vegf  <- exp(ldp_slope + etaldp_slope_vegf)
    dp_slope_tie2  <- exp(ldp_slope + etaldp_slope_tie2)
    dp_slope_ang2  <- exp(ldp_slope + etaldp_slope_ang2)
    dp_slope_fgf23 <- exp(ldp_slope + etaldp_slope_fgf23)

    emax <- exp(lemax + etalemax)

    # ---- Sigmoid Emax drug effect on the daily-AUC exposure metric
    # (Majid 2024 Equations 1 and 2). The deposited control stream adds a
    # 0.1 numerical guard to AUC inside this expression; it is dropped here
    # because it is a unit-scale-dependent epsilon rather than a model
    # feature (see the vignette Errata), and because dropping it makes the
    # placebo arm (AUC_LEN = 0) reduce exactly to the disease-progression
    # model the paper describes.
    eff_vegf  <- emax * AUC_LEN^hill_vegf  / (ec50_vegf^hill_vegf   + AUC_LEN^hill_vegf)
    eff_tie2  <- emax * AUC_LEN^hill_tie2  / (ec50_tie2^hill_tie2   + AUC_LEN^hill_tie2)
    eff_ang2  <- emax * AUC_LEN^hill_ang2  / (ec50_ang2^hill_ang2   + AUC_LEN^hill_ang2)
    eff_fgf23 <- emax * AUC_LEN^hill_fgf23 / (ec50_fgf23^hill_fgf23 + AUC_LEN^hill_fgf23)

    # ---- Linear disease progression on each baseline (Majid 2024
    # Equation 3: DP(t) = BM0 * (1 + DPslope * t), Kin = DP(t) * Kout).
    # rxode2 re-evaluates this at every integration step, so DP is a
    # continuous function of t, matching the printed equation; the NONMEM
    # stream re-evaluates it only at record times.
    dp_vegf  <- rbase_vegf  * (1 + dp_slope_vegf  * t)
    dp_tie2  <- rbase_tie2  * (1 + dp_slope_tie2  * t)
    dp_ang2  <- rbase_ang2  * (1 + dp_slope_ang2  * t)
    dp_fgf23 <- rbase_fgf23 * (1 + dp_slope_fgf23 * t)

    kout_vegf  <- 1 / mrt_vegf
    kout_tie2  <- 1 / mrt_tie2
    kout_ang2  <- 1 / mrt_ang2
    kout_fgf23 <- 1 / mrt_fgf23

    kin_vegf  <- dp_vegf  * kout_vegf
    kin_tie2  <- dp_tie2  * kout_tie2
    kin_ang2  <- dp_ang2  * kout_ang2
    kin_fgf23 <- dp_fgf23 * kout_fgf23

    # ---- Initial conditions: each pool starts at its individual baseline
    # (Majid 2024 Equations 1-2, 'BM(0): BM0'; Text S2 A_0(i) = BM0i).
    vegf(0)  <- rbase_vegf
    tie2(0)  <- rbase_tie2
    ang2(0)  <- rbase_ang2
    fgf23(0) <- rbase_fgf23

    # ---- Turnover ODEs. Lenvatinib inhibits Kout for VEGF and FGF-23, so
    # those pools accumulate (Equation 1); it inhibits Kin for Tie-2 and
    # Ang-2, so those pools deplete (Equation 2).
    d/dt(vegf)  <- kin_vegf  - kout_vegf * (1 - eff_vegf) * vegf
    d/dt(tie2)  <- kin_tie2  * (1 - eff_tie2)  - kout_tie2  * tie2
    d/dt(ang2)  <- kin_ang2  * (1 - eff_ang2)  - kout_ang2  * ang2
    d/dt(fgf23) <- kin_fgf23 - kout_fgf23 * (1 - eff_fgf23) * fgf23

    vegf  ~ add(addSd_vegf) + prop(propSd_vegf)
    tie2  ~ prop(propSd_tie2)
    ang2  ~ prop(propSd_ang2)
    fgf23 ~ prop(propSd_fgf23)
  })
}
