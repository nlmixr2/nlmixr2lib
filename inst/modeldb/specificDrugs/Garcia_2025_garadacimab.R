Garcia_2025_garadacimab <- function() {
  description <- "Two-compartment population PK model with first-order absorption and first-order elimination, coupled to a direct sigmoidal Emax inhibition PD model for FXIIa-mediated kallikrein activity, for the anti-activated-factor-XII (anti-FXIIa) monoclonal antibody garadacimab in healthy volunteers and patients with hereditary angioedema (HAE). Pooled across two phase I, one phase II, and two phase III studies. Absolute bioavailability of the subcutaneous depot is estimated on the logit scale with logit-scale IIV; intravenous doses go directly to the central compartment and bypass it. Body weight enters CL, Q, Vc, and Vp as an estimated power exponent (not fixed allometry), and Japanese heritage, Chinese heritage, HAE disease status, baseline serum creatinine, baseline alanine aminotransferase, and baseline total bilirubin enter CL. PD is a direct effect on garadacimab plasma concentration with no effect compartment. Sister model file from the same paper (HAE attack repeated-time-to-event exposure-response): modellib('Garcia_2025_garadacimab_hae_attack')."
  reference <- paste(
    "Garcia R, Cheng S, Glassman F, Sharma A, De Miguel-Lillo B, Wiens M,",
    "Johnston C, Lawo JP, Pragst I, French J, Polhamus D, Nandy P.",
    "Population pharmacokinetic/pharmacodynamic and exposure-response modeling",
    "of garadacimab in healthy volunteers and patients with hereditary angioedema.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(5):954-963.",
    "doi:10.1002/psp4.70009.",
    "Open Access under CC BY-NC 4.0.",
    "Structural equations and the NONMEM control streams are in Methods S1;",
    "final parameter estimates are in Tables S4-S7 of the supplement",
    "(Wiley file PSP4-14-954-s001.docx).",
    "Sister model file from the same paper:",
    "modellib('Garcia_2025_garadacimab_hae_attack').",
    sep = " "
  )
  vignette <- "Garcia_2025_garadacimab"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (the source column is baseline weight, BLWT). Enters CL and Q with one shared power exponent and Vc and Vp with a second shared power exponent, both centred at 70 kg: `(WT / 70)^e_wt_cl` and `(WT / 70)^e_wt_vc` (Methods S1 PopPK code: LCLWT = THETA(7)*LOG(BLWT/70), LQWT = THETA(7)*LOG(BLWT/70), LV2WT = THETA(8)*LOG(BLWT/70), LV3WT = THETA(8)*LOG(BLWT/70)). Note that the FINAL model ESTIMATES both exponents (1.16 and 0.843); the fixed-allometry values 0.75 and 1 belong to the BASE model only (Table S4). Observed range 43.3-153 kg, median 79.2 kg.",
      source_name        = "BLWT"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese-heritage indicator: 1 = Japanese, 0 = non-Japanese.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese)",
      notes              = "Time-fixed per subject. Multiplicative effect on CL (Methods S1 PopPK code: IF (JPN.EQ.1) LCLJPN = THETA(9), added in the log domain). Japanese subjects enter the pooled analysis from the phase I ethnobridging study NCT04580654. The paper concludes this effect is NOT clinically meaningful: although the point estimate falls marginally outside the 80-125% reference range, EBE simulations showed no difference between Japanese and non-Japanese patients with HAE (Results Section 3.2).",
      source_name        = "JPN"
    ),
    RACE_CHINESE = list(
      description        = "Chinese-heritage indicator: 1 = Chinese, 0 = non-Chinese.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Chinese)",
      notes              = "Time-fixed per subject. Multiplicative effect on CL (Methods S1 PopPK code: IF (CHN.EQ.1) LCLCHN = THETA(10), added in the log domain). The estimated effect is 1.02-fold, i.e. essentially null (Table S4).",
      source_name        = "CHN"
    ),
    DIS_HAE = list(
      description        = "Hereditary angioedema patient indicator: 1 = patient with HAE (HAE-C1INH-Type1, HAE-C1INH-Type2, or HAE-nC1INH), 0 = healthy volunteer.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer)",
      notes              = "Time-fixed per subject. Multiplicative effect on CL (Methods S1 PopPK code: IF (PAT.EQ.2) LCLPAT = THETA(11), added in the log domain; the source column PAT codes 2 = patient with HAE). The estimated effect is 1.05-fold; the paper reports model-predicted CL, Vc, and AUCtau,ss to be generally similar in patients with HAE and healthy volunteers (Results Section 3.3).",
      source_name        = "PAT"
    ),
    CREAT = list(
      description        = "Baseline serum creatinine.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Power effect on CL centred at 0.75 mg/dL (Methods S1 PopPK code: LCLCREAT = THETA(12)*LOG(BLCREAT/0.75)). US-convention units, NOT the SI umol/L: the reference subject is defined in Section 2.4 as having 'baseline serum creatinine of 0.75 mg/dL', and Figure 1 varies sCr over 0.6-1.0 mg/dL. The estimated exponent is -0.0343 with a 95% CI spanning zero, and the paper concludes there is no clinically meaningful effect on CL or AUCtau,ss (Results Section 3.2).",
      source_name        = "BLCREAT"
    ),
    ALT = list(
      description        = "Baseline serum alanine aminotransferase activity.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Power effect on CL centred at 25 U/L (Methods S1 PopPK code: LCLALT = THETA(13)*LOG(BLALT/25)). Reference subject baseline ALT 25 U/L (Section 2.4); Figure 1 varies ALT over 10-40 U/L. The paper concludes there is no clinically meaningful effect (Results Section 3.2).",
      source_name        = "BLALT"
    ),
    TBILI = list(
      description        = "Baseline total serum bilirubin.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Power effect on CL centred at 8 umol/L (Methods S1 PopPK code: LCLBILI = THETA(14)*LOG(BLBILI/8)). SI units: the reference subject is defined in Section 2.4 as having 'baseline bilirubin of 8 umol/L', and Figure 1 varies bilirubin over 4-13 umol/L, a range consistent with TOTAL rather than direct bilirubin. The paper concludes there is no clinically meaningful effect (Results Section 3.2).",
      source_name        = "BLBILI"
    )
  )

  # Covariates that the source screened, and in one case RETAINED in the final
  # model, but for which no usable point estimate is reported anywhere on disk.
  # Documentation only -- checkModelConventions() does not require these to be
  # referenced in model().
  covariatesDataExcluded <- list(
    PKK_ACTIVITY_BL = list(
      description        = "Baseline FXIIa-mediated kallikrein activity (the PD assay readout at baseline), reference value 0.134 POB.",
      units              = "POB (proportion of baseline, chromogenic S-2302 assay)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "RETAINED in the source's final PopPK/PD model as a power covariate on BOTH E0 and EC50, centred at 0.134 (Methods S1 PopPK/PD code: LE0PD = THETA(5)*LOG(BLPD2/0.134); LEC50PD = THETA(6)*LOG(BLPD2/0.134); missing values imputed to the 0.134 median). However, THETA(5) and THETA(6) are NOT reported: Table S6 tabulates only the four structural PD parameters, and the only rendering of the effect is the Figure S6 forest plot (EC50 only, no E0 panel), which is an image with no accompanying data table and whose x-axis anchors -- the 10th and 90th percentiles of baseline kallikrein activity -- are not tabulated in Table S3 or anywhere else on disk. The coefficients are therefore not digitisable. This model file is accordingly the REFERENCE-SUBJECT PopPK/PD model at baseline kallikrein activity = 0.134 POB, where both covariate terms are exactly zero by construction and the Table S6 structural estimates apply directly. The source's own conclusion is that this covariate is not clinically meaningful: the EC50 point estimates and 95% CIs normalised to the reference subject were almost fully contained within the 80-125% reference range (Results Section 3.4; Discussion). No canonical covariate column is proposed because the coefficient cannot be populated; note that the existing register entry PKK_BL is a DIFFERENT concept (baseline plasma prekallikrein CONCENTRATION in mg/L), not an enzymatic activity readout."
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "garadacimab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "garadacimab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "garadacimab", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 242L,
    n_studies      = 5L,
    age_range      = "12-73 years (median 41)",
    weight_range   = "43.3-153 kg (median 79.2)",
    sex_female_pct = NULL,
    race_ethnicity = NULL,
    disease_state  = "Pooled healthy volunteers and patients with hereditary angioedema due to C1-inhibitor deficiency or dysfunction (HAE-C1INH-Type1/Type2) or with normal C1 inhibitor (HAE-nC1INH, including HAE-FXII and HAE-PLG subtypes).",
    dose_range     = "Phase I single intravenous doses up to 10 mg/kg and subcutaneous doses; phase II subcutaneous 75, 200, or 600 mg once monthly (28 +/- 2 days); phase III subcutaneous 200 mg once monthly (30 +/- 4 days) after a loading dose of two 200 mg subcutaneous injections.",
    regions        = "multinational",
    biomarkers     = "FXIIa-mediated kallikrein activity, expressed as percent of baseline (POB), measured with an in-house chromogenic substrate (S-2302) enzymatic assay. Garadacimab plasma concentration measured by a validated clinical ELISA.",
    notes          = "PK analysis dataset: 242 unique participants who received at least one dose of garadacimab prior to one evaluable PK sample. The PopPK/PD (FXIIa-mediated kallikrein activity) dataset adds 20 unique placebo recipients to the same 242 garadacimab recipients. Pooled from five studies (Table S1): two phase I studies in healthy volunteers (ACTRN 12616001438448 and NCT04580654, the Japanese/White ethnobridging study), one phase II study in patients with HAE including its randomized placebo-controlled and open-label periods (NCT03712228), and two phase III studies in patients with HAE-C1INH (pivotal VANGUARD NCT04656418 and the open-label extension NCT04739059). Sex ratios were comparable across studies EXCEPT the phase I study ACTRN 12616001438448, which enrolled only healthy male volunteers per protocol, so no single pooled female percentage is reported by the source (Results Section 3.1; per-study counts in Table S3). Observations below the limit of quantification were excluded; covariates with >30% missing values were excluded and those with 10-30% missing were imputed by regression (Section 2.3)."
  )

  ini({
    # ------------------------------------------------------------------
    # PK structural parameters -- Table S4 FINAL model column. Reference
    # subject (Section 2.4 and the Figure 1 caption): non-Japanese,
    # non-Chinese healthy volunteer, baseline weight 70 kg, baseline serum
    # creatinine 0.75 mg/dL, baseline ALT 25 U/L, baseline bilirubin
    # 8 umol/L.
    #
    # CL and Q are TRUE (not apparent) clearances and Vc and Vp are TRUE
    # volumes: the pooled dataset contains intravenous as well as
    # subcutaneous data, so absolute bioavailability F1 is identifiable and
    # is estimated separately. Table 1 of the main text reports the derived
    # APPARENT quantities CL/F and Vc/F for patients, which is why those
    # values are roughly 1/0.387 times larger.
    # ------------------------------------------------------------------
    lcl <- log(0.00664); label("Clearance at reference covariates (CL, L/h)")                    # Table S4 final model CL = 0.00664 L/h (95% CI 0.00573-0.00769)
    lvc <- log(2.37);    label("Central volume of distribution at reference weight (V2, L)")     # Table S4 final model V2 = 2.37 L (95% CI 1.77-3.17)
    lq  <- log(0.00685); label("Intercompartmental clearance at reference weight (Q, L/h)")      # Table S4 final model Q = 0.00685 L/h (95% CI 0.00588-0.00799)
    lvp <- log(1.41);    label("Peripheral volume of distribution at reference weight (V3, L)")  # Table S4 final model V3 = 1.41 L (95% CI 1.34-1.50)
    lka <- log(0.00824); label("First-order subcutaneous absorption rate constant (ka, 1/h)")    # Table S4 final model ka = 0.00824 (95% CI 0.00763-0.00889); the Table S4 unit label "(L/h)" is a typo -- ka is a first-order rate constant, 1/h

    # Absolute bioavailability of the SC depot, estimated on the logit
    # scale (Methods S1 PopPK code: LOGITTVF1 = LOG(TVF1/(1-TVF1)),
    # F1 = EXPLOGITF1/(1+EXPLOGITF1), with ETA(3) added on the logit
    # scale). F1 applies to compartment 1 (the depot) only, so IV doses
    # administered into the central compartment are unaffected.
    logitfdepot <- log(0.387 / (1 - 0.387)); label("Logit of absolute subcutaneous bioavailability (unitless; F1 = 0.387)")  # Table S4 final model F1 = 0.387 (95% CI 0.344-0.431)

    # ------------------------------------------------------------------
    # PK covariate effects -- Table S4 FINAL model column.
    #
    # SCALE WARNING: Table S4 reports the two kinds of covariate effect on
    # DIFFERENT scales, which is confirmed both by the initial estimates in
    # the Methods S1 control stream and by reproducing all 15 rows of the
    # Figure 1 forest plot (see the vignette's source-trace section):
    #   * CONTINUOUS covariates (weight, sCr, ALT, bilirubin) enter as
    #     THETA*LOG(cov/ref), so the tabulated number IS the power
    #     EXPONENT and is used here as printed.
    #   * CATEGORICAL covariates (Japanese, Chinese, patient) enter as a
    #     THETA added in the log domain, and the tabulated number is the
    #     back-transformed FOLD CHANGE exp(THETA). These are therefore
    #     entered here as log(fold change).
    # ------------------------------------------------------------------
    e_wt_cl <- 1.16;  label("Power exponent of body weight on CL and Q, centred at 70 kg (unitless)")        # Table S4 final model "Weight effect on CL and Q" = 1.16 (95% CI 0.963-1.360); ESTIMATED, replacing the base model's fixed 0.750
    e_wt_vc <- 0.843; label("Power exponent of body weight on Vc and Vp, centred at 70 kg (unitless)")       # Table S4 final model "Weight effect on V2 and V3" = 0.843 (95% CI 0.657-1.03); ESTIMATED, replacing the base model's fixed 1.00

    e_race_japanese_cl <- log(1.27); label("Log fold change in CL for Japanese vs non-Japanese subjects")    # Table S4 final model "Japanese effect on CL" = 1.27-fold (95% CI 1.11-1.44)
    e_race_chinese_cl  <- log(1.02); label("Log fold change in CL for Chinese vs non-Chinese subjects")      # Table S4 final model "Chinese effect on CL" = 1.02-fold (95% CI 0.861-1.21)
    e_dis_hae_cl       <- log(1.05); label("Log fold change in CL for patients with HAE vs healthy volunteers") # Table S4 final model "Patient effect on CL" = 1.05-fold (95% CI 0.952-1.17)

    e_creat_cl <- -0.0343; label("Power exponent of baseline serum creatinine on CL, centred at 0.75 mg/dL (unitless)") # Table S4 final model "sCR effect on CL" = -0.0343 (95% CI -0.214 to 0.146)
    e_alt_cl   <- -0.0773; label("Power exponent of baseline ALT on CL, centred at 25 U/L (unitless)")                  # Table S4 final model "ALT effect on CL" = -0.0773 (95% CI -0.148 to -0.0067)
    e_tbili_cl <- -0.136;  label("Power exponent of baseline total bilirubin on CL, centred at 8 umol/L (unitless)")    # Table S4 final model "Bilirubin effect on CL" = -0.136 (95% CI -0.203 to -0.069)

    # ------------------------------------------------------------------
    # PK interindividual variability -- Table S5 FINAL model column, a full
    # 3x3 OMEGA BLOCK on CL, Vc, and the LOGIT of F1, in that order
    # (Methods S1 PopPK code: $OMEGA BLOCK(3) over ETA(1) CL, ETA(2) V2,
    # ETA(3) logit-F1). Variances are on the log scale for CL and Vc and on
    # the LOGIT scale for F1; the "SD 0.135" annotation Table S5 attaches to
    # IIV-F1 is the delta-method SD of F1 on its natural 0-1 scale, not the
    # value entered here. Q, V3, and ka carry no IIV ($OMEGA 0 FIX).
    #
    # Reported correlations round-trip exactly against these variances and
    # covariances: 0.263/sqrt(0.175*0.717) = 0.743, 0.154/sqrt(0.359*0.175)
    # = 0.614, and 0.174/sqrt(0.359*0.717) = 0.343.
    # ------------------------------------------------------------------
    # Table S5 final model, in lower-triangle order:
    #   IIV-CL          = 0.175  (CV% 43.8; 95% CI 0.098-0.253)
    #   V2-CL covariance= 0.263  (Corr 0.743)
    #   IIV-V2          = 0.717  (CV% 102)
    #   F1-CL covariance= 0.154  (Corr 0.614)
    #   F1-V2 covariance= 0.174  (Corr 0.342)
    #   IIV-F1          = 0.359  (logit scale)
    # (comments must sit OUTSIDE the c() -- rxode2 fails to parse a trailing
    # comment inside a multi-line omega c() block.)
    etalcl + etalvc + etalogitfdepot ~ c(
      0.175,
      0.263, 0.717,
      0.154, 0.174, 0.359
    )

    # ------------------------------------------------------------------
    # PD structural parameters -- Table S6 FINAL model column. Direct
    # sigmoidal Emax INHIBITION of FXIIa-mediated kallikrein activity by
    # garadacimab plasma concentration, with no effect compartment
    # (Methods S1 PopPK/PD code $ERROR block).
    # ------------------------------------------------------------------
    lemax <- log(0.988);  label("Maximum fractional inhibition of FXIIa-mediated kallikrein activity (Emax, unitless)") # Table S6 final model Emax = 0.988 (95% CI 0.981-0.996)
    lec50 <- log(17600);  label("Garadacimab concentration producing half-maximal inhibition (EC50, ng/mL)")            # Table S6 final model EC50 = 17600 ng/mL (95% CI 16400-18800)
    le0   <- log(98.8);   label("Baseline FXIIa-mediated kallikrein activity (E0, percent of baseline)")                # Table S6 final model E0 = 98.8 % of baseline (95% CI 94.9-103)
    lhill <- log(2.05);   label("Hill coefficient of the inhibitory sigmoid (gamma, unitless)")                         # Table S6 final model Hill coefficient = 2.05 (95% CI 1.88-2.23)

    # ------------------------------------------------------------------
    # PD interindividual variability -- Table S7 FINAL model column, a full
    # 3x3 OMEGA BLOCK on EC50, E0, and gamma, in that order (Methods S1
    # PopPK/PD code: $OMEGA BLOCK(3) over ETA(2) EC50, ETA(3) E0,
    # ETA(4) GAMMA). Emax carries no IIV ($OMEGA 0 FIX).
    #
    # ERRATUM -- Table S7 COVARIANCE ROW LABELS ARE SWAPPED. As printed,
    # the row "gamma-E0 = 0.0846" would give a correlation of
    # 0.0846/sqrt(0.150*0.0268) = 1.334, which is mathematically impossible.
    # The table's own correlation column identifies the correct pairing in
    # BOTH the base-model and final-model columns:
    #   final: 0.0846/sqrt(0.150*0.168 ) = 0.533 = the printed "gamma-E0"  Corr
    #          0.0254/sqrt(0.150*0.0268) = 0.401 = the printed "gamma-E50" Corr
    #   base : 0.0716/sqrt(0.137*0.175 ) = 0.462 = the printed "gamma-E0"  Corr
    #          0.0299/sqrt(0.137*0.0415) = 0.397 = the printed "gamma-E50" Corr
    # So the value printed against "gamma-E0" is cov(gamma, EC50) and the
    # value printed against "gamma-E50" is cov(gamma, E0). The corrected
    # assignment is used below and yields a positive-definite matrix
    # (eigenvalues 0.2476, 0.0754, 0.0217). The uncorrected reading does
    # not. The analogous round-trip on the PK block (Table S5) reproduces
    # all three printed correlations with the labels AS PRINTED, so the
    # swap is specific to Table S7.
    # ------------------------------------------------------------------
    # Table S7 final model, in lower-triangle order:
    #   IIV-EC50           = 0.168  (CV% 42.8; 95% CI 0.119-0.217)
    #   EC50-E0 covariance = 0.0148 (Corr 0.221)
    #   IIV-E0             = 0.0268 (CV% 16.5)
    #   cov(gamma, EC50)   = 0.0846 -- printed against "gamma-E0",  see erratum above
    #   cov(gamma, E0)     = 0.0254 -- printed against "gamma-E50", see erratum above
    #   IIV-gamma          = 0.150  (CV% 40.2)
    etalec50 + etale0 + etalhill ~ c(
      0.168,
      0.0148, 0.0268,
      0.0846, 0.0254, 0.150
    )

    # ------------------------------------------------------------------
    # Residual error. Table S5 and Table S7 report VARIANCES; nlmixr2 takes
    # standard deviations, so each is entered as sqrt(variance). The
    # reported CV% / SD annotations confirm the scale:
    #   sqrt(0.0395) = 0.1987 -> CV% 19.9 (Table S5)
    #   sqrt(0.0577) = 0.2402 -> CV% 24.0 (Table S7)
    #   sqrt(6.01)   = 2.4515 -> SD  2.45 (Table S7)
    # ------------------------------------------------------------------
    propSd            <- 0.198746; label("Proportional residual error SD on garadacimab plasma concentration (unitless)")                   # Table S5 final model proportional residual variance = 0.0395 (CV% 19.9)
    propSd_kallikrein <- 0.240208; label("Proportional residual error SD on FXIIa-mediated kallikrein activity (unitless)")                 # Table S7 final model proportional residual variance = 0.0577 (CV% 24.0)
    addSd_kallikrein  <- 2.451530; label("Additive residual error SD on FXIIa-mediated kallikrein activity (percent of baseline)")          # Table S7 final model additive residual variance = 6.01 (SD 2.45)
  })

  model({
    # --- PK individual parameters ------------------------------------
    # Body weight enters CL and Q with one shared exponent and Vc and Vp
    # with a second shared exponent. The categorical covariates and the
    # remaining continuous covariates act on CL only.
    cl <- exp(lcl + etalcl) *
      (WT / 70)^e_wt_cl *
      exp(e_race_japanese_cl * RACE_JAPANESE +
            e_race_chinese_cl * RACE_CHINESE +
            e_dis_hae_cl * DIS_HAE) *
      (CREAT / 0.75)^e_creat_cl *
      (ALT / 25)^e_alt_cl *
      (TBILI / 8)^e_tbili_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    ka <- exp(lka)

    # Absolute bioavailability of the SC depot (logit scale, logit-scale IIV).
    fdepot <- expit(logitfdepot + etalogitfdepot)

    # --- micro-constants ---------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # --- disposition --------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # F1 applies to the SC depot only; IV doses go straight into `central`
    # and bypass it (NONMEM ADVAN4 F1 on compartment 1).
    f(depot) <- fdepot

    # Amounts are in mg and volumes in L, so central/vc is mg/L; the factor
    # 1000 converts to ng/mL, matching the source's S2 = V2/1000 scaling
    # ("Dose in mg, conc in ng/mL").
    Cc <- central / vc * 1000

    # --- PD: direct sigmoidal Emax inhibition -------------------------
    emax <- exp(lemax)
    ec50 <- exp(lec50 + etalec50)
    e0   <- exp(le0 + etale0)
    hill <- exp(lhill + etalhill)

    # Numerical guard only: the Hill exponent is non-integer (2.05), so a
    # transiently negative solver value for Cc would make Cc^hill NaN and
    # abort the solve. Clamping at zero is a no-op for all physically
    # meaningful states.
    cpd <- max(Cc, 0)

    kallikrein <- e0 * (1 - emax * cpd^hill / (ec50^hill + cpd^hill))

    Cc         ~ prop(propSd)
    kallikrein ~ prop(propSd_kallikrein) + add(addSd_kallikrein)
  })
}
