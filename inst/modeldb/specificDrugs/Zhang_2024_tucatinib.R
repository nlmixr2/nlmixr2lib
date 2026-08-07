Zhang_2024_tucatinib <- function() {
  description <- "Two-compartment population PK model for oral tucatinib (HER2-selective tyrosine kinase inhibitor) in healthy participants and patients with HER2-positive metastatic breast cancer or HER2-positive metastatic colorectal cancer (Zhang 2024, n=283 across seven phase I/II studies). First-order absorption into a depot preceded by an absorption lag time, linear elimination from central, and two-compartment disposition (central + peripheral1). Tumor type is the only covariate retained in the final model: metastatic breast cancer and metastatic colorectal cancer each lower CL relative to the healthy-participant reference, and metastatic colorectal cancer additionally lowers relative bioavailability, so that the net apparent clearance CL/F is 51.9% lower in mBC and 18.7% lower in mCRC. Inter-individual variability on CL and Vc is correlated; Vc carries a very large 135% CV."
  reference <- paste(
    "Zhang D, Taylor A, Zhao JJ, Endres CJ, Topletz-Erickson A.",
    "Population Pharmacokinetic Analysis of Tucatinib in Healthy Participants",
    "and Patients with Breast Cancer or Colorectal Cancer.",
    "Clin Pharmacokinet. 2024;63(10):1477-1487.",
    "doi:10.1007/s40262-024-01412-0.",
    "Correction: Clin Pharmacokinet. 2025;64(1):171-172.",
    "doi:10.1007/s40262-024-01458-0",
    "(corrects the Table 2 model-predicted Ctrough,ss geometric mean for mCRC;",
    "no final-model parameter in Table 1 is affected).",
    sep = " "
  )
  vignette <- "Zhang_2024_tucatinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the Zhang 2024 supplementary NONMEM
  # control stream ($SUBROUTINE ADVAN4 TRANS4: compartment 1 = depot,
  # compartment 2 = central with S2 = V2/1000 for a mg dose read out in ng/mL,
  # compartment 3 = peripheral) and against the paper's Methods description of
  # tucatinib plasma PK.
  compartmentData <- list(
    depot       = list(analyte = "tucatinib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "tucatinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tucatinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    TUMTP_BREAST = list(
      description        = "HER2-positive metastatic breast cancer indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not mBC; a healthy participant is the typical-value reference when paired with TUMTP_CRC = 0)",
      notes              = paste(
        "Time-fixed. Multiplicative fractional effect on CL only:",
        "-0.519, i.e. 51.9% lower clearance than a healthy participant",
        "(Zhang 2024 Table 1). No effect on relative bioavailability was",
        "retained: the paper reports that the mBC effect on Frel was removed",
        "because it was poorly estimated (RSE 50.4%) and increased the",
        "objective function value (Results section 3.3). Because Frel is",
        "unchanged, the net apparent clearance CL/F is also 51.9% lower,",
        "giving the published 2.1-fold AUCss increase versus healthy",
        "participants. HER2+ mBC patients were 63 / 283 = 22.3% of the",
        "pooled cohort (studies ONT-380-004 and ONT-380-005; Table S2).",
        sep = " "
      ),
      source_name        = "TUM"
    ),
    TUMTP_CRC = list(
      description        = "HER2-positive metastatic colorectal cancer indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not mCRC; a healthy participant is the typical-value reference when paired with TUMTP_BREAST = 0)",
      notes              = paste(
        "Time-fixed. Multiplicative fractional effects on BOTH CL (-0.705,",
        "i.e. 70.5% lower clearance) and relative bioavailability",
        "(-0.637, i.e. 63.7% lower Frel) versus a healthy participant",
        "(Zhang 2024 Table 1). The two effects partly offset: the net",
        "apparent clearance CL/F = CL / Frel is 112 * 0.295 / 0.363 =",
        "91.0 L/h, i.e. 18.7% lower than the healthy-participant 112 L/h,",
        "which is the figure quoted in Results section 3.4. HER2+ mCRC",
        "patients were 69 / 283 = 24.4% of the pooled cohort (study",
        "SGNTUC-017 / MOUNTAINEER; Table S2).",
        sep = " "
      ),
      source_name        = "TUM"
    )
  )

  # Covariates screened by Zhang 2024 but NOT retained in the final model.
  # Documented here (rather than in covariateData) so the provenance of the
  # paper's covariate search is preserved without declaring covariates that
  # model() never references. Every one of these entered the stepwise search
  # and dropped out during backward elimination, except where noted.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested on CL/F and Vc/F in the univariable screen; not retained (Zhang 2024 section 3.3). Cohort median 75.9 kg [40.7, 146] (Table S3)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL/F, Frel and Ka. Sex on CL/F was added during primary forward addition but dropped during backward elimination (Zhang 2024 section 3.3). Cohort 39.9% female (Table S2)."
    ),
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = "Tested on CL/F and Vc/F; not retained (Zhang 2024 section 3.3). Cohort median 48 years [18, 77] (Table S3)."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Tested on CL/F and Vc/F. Albumin on Vc/F was added during primary forward addition but dropped during backward elimination (Zhang 2024 section 3.3). Cohort median 4.30 g/dL [1.87, 5.20] (Table S3), reported in US convention rather than the canonical g/L."
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL/F and Vc/F. Race on CL/F was added during primary forward addition but dropped during backward elimination (Zhang 2024 section 3.3). Cohort 13.1% Black (Table S2)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL/F and Vc/F; not retained (Zhang 2024 section 3.3). Cohort 8.8% Asian (Table S2). Study SGNTUC-015 specifically compared Japanese and Caucasian participants and found similar tucatinib PK."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance (Cockcroft-Gault)",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Tested on CL/F; not retained (Zhang 2024 section 3.3). Only BASELINE creatinine was used, because tucatinib inhibits OCT2/MATE-mediated tubular creatinine secretion and therefore raises on-treatment serum creatinine without renal impairment (Methods section 2.2.2). Cohort median 111 mL/min [45.5, 280] (Table S3). Mild renal impairment had no effect; moderate/severe could not be evaluated."
    ),
    NCI_HEPATIC = list(
      description = "National Cancer Institute liver dysfunction category",
      units       = "(categorical: normal / mild)",
      type        = "categorical",
      notes       = "Tested on CL/F. Added during primary forward addition but dropped during backward elimination (Zhang 2024 section 3.3). Only normal (84.5%) and mild (14.8%) categories were represented (Table S2); moderate/severe hepatic impairment could not be evaluated. Not a canonical register column - this model does not reference it in model()."
    ),
    ECOG = list(
      description = "Eastern Cooperative Oncology Group performance status",
      units       = "(categorical: 0 / 1 / 2)",
      type        = "categorical",
      notes       = "Tested on CL/F. Added during primary forward addition but dropped during backward elimination (Zhang 2024 section 3.3). Cohort 80.6% ECOG 0, 0.7% ECOG 1, 18.7% ECOG 2 (Table S2)."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 283L,
    n_observations   = 3942L,
    n_studies        = 7L,
    age_range        = "18-77 years",
    age_median       = "48 years (pooled across the seven studies)",
    weight_range     = "40.7-146 kg",
    weight_median    = "75.9 kg (pooled across the seven studies)",
    sex_female_pct   = 39.9,
    race_ethnicity   = "White 72.4%, Black 13.1%, Asian 8.8%, American Indian 0.4%, Native Hawaiian 0.4%, missing 4.9% (Zhang 2024 Table S2)",
    disease_state    = "Pooled population: healthy participants (151 / 283 = 53.4%); HER2+ metastatic breast cancer (63 / 283 = 22.3%); HER2+ metastatic colorectal cancer (69 / 283 = 24.4%).",
    dose_range       = "Oral tucatinib tablet: 300 mg single dose and 300 mg BID in healthy participants; 50, 150 and 300 mg BID in the Japanese/Caucasian study SGNTUC-015; 300 mg BID (a minority at 350 mg BID) in the patient studies.",
    regions          = "Multinational; SGNTUC-015 enrolled matched Japanese and Caucasian cohorts. Regional breakdown otherwise not reported.",
    renal_function   = "Mild renal impairment represented; moderate and severe renal impairment not evaluable (insufficient data).",
    hepatic_function = "Normal (84.5%) and NCI mild (14.8%) liver dysfunction only; moderate and severe hepatic impairment not evaluable.",
    co_medication    = "Patients received tucatinib in combination with trastuzumab (SGNTUC-017), T-DM1 (ONT-380-004), or capecitabine +/- trastuzumab (ONT-380-005). CYP2C8-modifying drugs were prohibited across the tucatinib studies and were therefore not tested as covariates.",
    studies          = c("ARRAY-380-103 (healthy, 300 mg single dose, n=12)",
                         "ONT-380-012 / NCT03723395 (healthy DDI, 300 mg SD and BID, n=86)",
                         "SGNTUC-015 / NCT03914755 (healthy Japanese vs Caucasian, 50/150/300 mg BID, n=36)",
                         "SGNTUC-020 / NCT03826602 (healthy metformin DDI, 300 mg BID, n=17)",
                         "SGNTUC-017 / NCT03043313 MOUNTAINEER (HER2+ mCRC, 300 mg BID, n=69)",
                         "ONT-380-004 / NCT01983501 (HER2+ mBC + T-DM1, 300-350 mg BID, n=39)",
                         "ONT-380-005 / NCT02025192 (HER2+ mBC + capecitabine, 300-350 mg BID, n=24)"),
    notes            = paste(
      "Baseline demographics from Zhang 2024 Tables S2 and S3. Only fasted",
      "marketed-tablet data were included from the healthy-participant",
      "studies; the fed arm of ARRAY-380-103 (n = 11) was too small to support",
      "a food-effect covariate. HER2CLIMB (ONT-380-206) was excluded because",
      "time since the previous dose was not recorded. Patient PK in",
      "ONT-380-004 and ONT-380-005 was collected only to 6 h post-dose, and",
      "SGNTUC-017 contributed sparse PK with diary-reported dosing - the",
      "authors flag both as limitations on the patient-arm exposure estimates.",
      sep = " "
    )
  )

  ini({
    # ---- Structural population parameters (Zhang 2024 Table 1, "Final
    # population PK model estimate" column) ----
    # Reference participant is a HEALTHY participant (Table 1 footnote), i.e.
    # TUMTP_BREAST = 0 and TUMTP_CRC = 0. Apparent oral parameters are reported
    # as CL/F (L/h), Vc/F (L), Q/F (L/h), Vp/F (L); Ka in 1/h; Tlag in h.
    # The supplementary NONMEM control stream fixes F1 = 1 at the reference
    # (TF1 = 1), so the reported CL/F is numerically the model's `cl`.
    lcl     <- log(112)   ; label("Apparent oral clearance in healthy participants (CL/F, L/h)")        # Zhang 2024 Table 1 CL/F healthy = 112 (4.3% RSE); bootstrap 112 (105, 119)
    lvc     <- log(125)   ; label("Apparent central volume of distribution (Vc/F, L)")                  # Zhang 2024 Table 1 Vc/F = 125 (11.3% RSE); bootstrap 104 (62.4, 145)
    lq      <- log(89.9)  ; label("Apparent inter-compartmental clearance (Q/F, L/h)")                  # Zhang 2024 Table 1 Q/F = 89.9 (5.4% RSE); bootstrap 83.4 (72.4, 97)
    lvp     <- log(635)   ; label("Apparent peripheral volume of distribution (Vp/F, L)")               # Zhang 2024 Table 1 Vp/F = 635 (4.0% RSE); bootstrap 659 (594, 740)
    lka     <- log(0.424) ; label("First-order absorption rate constant (Ka, 1/h)")                     # Zhang 2024 Table 1 Ka = 0.424 (5.3% RSE); bootstrap 0.385 (0.333, 0.453)
    ltlag   <- log(0.392) ; label("Absorption lag time (Tlag, h)")                                      # Zhang 2024 Table 1 Tlag = 0.392 (0.56% RSE); bootstrap 0.395 (0.385, 0.445)

    # Relative bioavailability at the reference (healthy participant) is the
    # structural anchor F1 = 1: the supplementary control stream sets
    # `TF1 = 1` and estimates only the mCRC fractional change on it. Tucatinib
    # was given orally in every included study, so no absolute-F information
    # is present in the data.
    lfdepot <- fixed(log(1)) ; label("Relative oral bioavailability in healthy participants (Frel)")    # Zhang 2024 supplementary NONMEM control stream, $PK block: TF1 = 1

    # ---- Tumor-type covariate effects (multiplicative fractional change) ----
    # Control-stream form, $PK block:
    #   CLTUM = 1 + THETA(9)  when TUM == 0  (HER2+ mBC)
    #   CLTUM = 1 + THETA(10) when TUM == 1  (HER2+ mCRC)
    #   CLTUM = 1             when TUM <  0  (healthy participant, reference)
    #   F1TUM = 1 + THETA(11) when TUM == 1  (HER2+ mCRC only)
    # i.e. TV_with_cov = TV_ref * (1 + e_<cov>_<param> * IND), so the printed
    # coefficient IS the fractional change at IND = 1. Reference = healthy.
    e_tumtp_breast_cl  <- -0.519 ; label("HER2+ mBC (vs healthy participant) fractional change on CL (unitless)")     # Zhang 2024 Table 1 fractional change in CL/F, HER2+ mBC = -0.519 (6.3% RSE); bootstrap -0.523 (-0.578, -0.456)
    e_tumtp_crc_cl     <- -0.705 ; label("HER2+ mCRC (vs healthy participant) fractional change on CL (unitless)")    # Zhang 2024 Table 1 fractional change in CL/F, HER2+ mCRC = -0.705 (3.6% RSE); bootstrap -0.730 (-0.819, -0.618)
    e_tumtp_crc_fdepot <- -0.637 ; label("HER2+ mCRC (vs healthy participant) fractional change on Frel (unitless)")  # Zhang 2024 Table 1 fractional change in Frel, HER2+ mCRC = -0.637 (3.8% RSE); bootstrap -0.673 (-0.789, -0.512)

    # ---- Inter-individual variability ----
    # Table 1 reports BSV on the CV% scale, with the footnote
    #   "BSV CV% is calculated as sqrt(exp(variance) - 1) x 100%",
    # so omega^2 = log(CV^2 + 1):
    #   CL/F 45%    -> log(1 + 0.45^2)  = 0.184403
    #   Vc/F 135%   -> log(1 + 1.35^2)  = 1.037623
    #   Q/F  43.9%  -> log(1 + 0.439^2) = 0.176237
    #   Vp/F 42.2%  -> log(1 + 0.422^2) = 0.163889
    # The footnote formula is corroborated by the supplementary control
    # stream's $OMEGA initial estimates (0.15 -> 40.2%, 0.99 -> 130.0%),
    # which sit next to the final 45% and 135%.
    #
    # Table 1's off-diagonal row is labelled "Cov (CL/F, Vc/F) 0.470" but is
    # encoded here as a CORRELATION, not a covariance. Treated as a
    # covariance it implies a correlation of 0.470 / sqrt(0.184403 * 1.037623)
    # = 1.074, so Omega would not be positive semi-definite; the largest
    # covariance the reported CVs admit is 0.4429 even at the extremes of
    # their rounding intervals. A correlation of 0.470 is feasible and is
    # consistent with the rest of the Table 1 BSV block, which is reported
    # entirely on the normalised (CV%) scale rather than as raw variances.
    # Encoded covariance = 0.470 * sqrt(0.184403 * 1.037623) = 0.205590.
    # See the vignette "Assumptions and deviations" section.
    etalcl + etalvc ~ c(0.184403,
                        0.205590, 1.037623)  # Zhang 2024 Table 1: BSV CV CL/F 45% (4.6% RSE), BSV CV Vc/F 135% (5.7% RSE), Cov(CL/F, Vc/F) 0.470 (6.8% RSE) read as a correlation
    etalq  ~ 0.176237                        # Zhang 2024 Table 1 BSV CV Q/F  = 43.9% (7.6% RSE) -> log(1 + 0.439^2)
    etalvp ~ 0.163889                        # Zhang 2024 Table 1 BSV CV Vp/F = 42.2% (7.4% RSE) -> log(1 + 0.422^2)

    # ---- Residual error ----
    # Combined additive + proportional, matching the control stream's
    # W = SQRT(WA*WA + WP*WP) with WA = THETA(7) and WP = THETA(8)*IPRED and
    # $SIGMA 1 FIX, i.e. THETA(7) is an additive SD in ng/mL and THETA(8) is a
    # proportional CV. Concentrations are in ng/mL (control stream S2 =
    # V2/1000 with a mg dose).
    addSd  <- 1.07  ; label("Additive residual error SD (ng/mL)")            # Zhang 2024 Table 1 additive residual error SD = 1.07 ng/mL (2.14% RSE); bootstrap 1.06 (0.0693, 1.57)
    propSd <- 0.329 ; label("Proportional residual error (fraction)")        # Zhang 2024 Table 1 proportional residual error CV = 0.329 (0.73% RSE); bootstrap 0.326 (0.309, 0.343)
  })

  model({
    # ---- Tumor-type multipliers (multiplicative fractional change; the
    # healthy-participant reference has both indicators 0 and multiplier 1) ----
    cl_tumtp <- 1 +
                e_tumtp_breast_cl * TUMTP_BREAST +
                e_tumtp_crc_cl    * TUMTP_CRC
    f_tumtp  <- 1 +
                e_tumtp_crc_fdepot * TUMTP_CRC

    # ---- Individual PK parameters ----
    # The control stream places ETAs on CL, V2, Q and V3 only; KA and ALAG1
    # carry no inter-individual variability.
    ka     <- exp(lka)
    cl     <- exp(lcl + etalcl) * cl_tumtp
    vc     <- exp(lvc + etalvc)
    q      <- exp(lq  + etalq)
    vp     <- exp(lvp + etalvp)
    tlag   <- exp(ltlag)
    fdepot <- exp(lfdepot) * f_tumtp

    # ---- Micro-constants ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- ODE system: first-order absorption + 2-compartment disposition
    # (NONMEM ADVAN4 TRANS4) ----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- fdepot
    alag(depot) <- tlag

    # ---- Plasma concentration ----
    # central is in mg and vc in L -> mg/L; multiply by 1000 to report ng/mL,
    # matching the control stream's S2 = V2/1000 ("dose mg conc ng/mL").
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
