Pohl_2022_linzagolix_e2 <- function() {
  description <- paste0(
    "Pharmacokinetic-oestradiol (PK-E2) exposure-response model for linzagolix, ",
    "an oral gonadotropin-releasing hormone (GnRH) receptor antagonist developed ",
    "for endometriosis, in healthy women and women with endometriosis (Pohl 2022). ",
    "A sigmoid inhibitory Emax function of the individual linzagolix daily AUC ",
    "(AUC_LZGX, supplied as a per-subject covariate column from the companion ",
    "population PK model Pohl_2022_linzagolix) suppresses an estimated baseline ",
    "oestradiol that differs between patients and healthy volunteers and carries ",
    "weight, age and race effects. An apparent placebo-arm E2 rise -- attributed ",
    "by the authors to loss of menstrual-cycle synchronisation after a ",
    "menses-anchored treatment start -- is carried as a saturating drift ",
    "compartment, and an additive term accounts for hormonal add-back therapy. ",
    "Log-transform-both-sides residual error with separate magnitudes for ",
    "patients and healthy volunteers."
  )
  reference <- paste(
    "Pohl O, Baron K, Riggs M, French J, Garcia R, Gotteland J-P.",
    "A model-based analysis to guide gonadotropin-releasing hormone receptor",
    "antagonist use for management of endometriosis.",
    "British Journal of Clinical Pharmacology. 2022;88(5):2359-2371.",
    "doi:10.1111/bcp.15171.",
    "Model equations are in the Supporting Information (sections 3.2 and 4.2);",
    "parameter values are in Table 5 of the main article."
  )
  vignette <- "Pohl_2022_linzagolix"

  # Non-canonical residual-SD names carried by this paper: the log-scale
  # residual error magnitude differs between endometriosis patients and healthy
  # volunteers (Table 5, "Residual variability (additive on log-scale)").
  paper_specific_residual_sds <- c("expSdPatient", "expSdHealthy")

  units <- list(
    time          = "h",
    dosing        = paste(
      "n/a (no drug dosing events; linzagolix exposure enters as the",
      "per-subject AUC_LZGX covariate from the companion population PK model)"
    ),
    concentration = "pg/mL (serum oestradiol; observation output -- not a drug concentration)",
    AUC_LZGX      = "ng*h/mL (linzagolix AUC over the 24 h once-daily dosing interval)"
  )

  covariateData <- list(
    AUC_LZGX = list(
      description        = "Linzagolix AUC over the 24 h once-daily dosing interval, supplied as a per-subject drug-exposure covariate from the companion population PK model.",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The source authors fit the PK-E2 model sequentially: individual daily",
        "AUC is generated from the population PK model of Table 4",
        "(Pohl_2022_linzagolix.R) and passed to this model as the exposure",
        "driver (Supporting Information section 3.2, 'linzagolix daily AUC was",
        "assumed to drive a proportional inhibitory change from estimated",
        "baseline over time'). For linear PK under once-daily dosing the",
        "steady-state value is Dose / (CL/F); Figure 2 of the source shows the",
        "resulting AUCss histograms by dose group. AUC_LZGX = 0 reproduces the",
        "placebo arm (the inhibitory term vanishes). NOTE the unit: Table 5",
        "reports AUC50 in ng*h/mL, so AUC_LZGX must be supplied in ng*h/mL",
        "(= 1000 x the ug*h/mL that the companion PK model's ug/mL",
        "concentrations integrate to)."
      ),
      source_name        = "linzagolix daily AUC"
    ),
    BL_E2 = list(
      description        = "Observed baseline serum oestradiol concentration for the subject.",
      units              = "pg/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The OBSERVED per-subject baseline E2, normalised to 52 pg/mL in the",
        "model equation, and distinct from the model-estimated baseline",
        "(exp(lbl_e2_*)). Supporting Information section 3.2: 'The linzagolix",
        "drug effect was assumed to be related to observed baseline E2 by an",
        "estimated exponent.' Table 3 gives study medians (53.0 pg/mL in",
        "EDELWEISS, 44.0 in KLH1204, 25.0 and 20.5 in the healthy-volunteer",
        "studies). The 52 pg/mL normalisation constant is printed in the",
        "Supporting Information equation itself, not in any table."
      ),
      source_name        = "baseline E2"
    ),
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on the estimated baseline oestradiol, normalised to a 58 kg reference subject (Supporting Information section 3.2 baseline-E2 equation; Table 5 row 'Baseline E2 ~ (weight 58 kg)' = -0.699).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on the estimated baseline oestradiol, normalised to a 35 year reference subject (Supporting Information section 3.2 baseline-E2 equation; Table 5 row 'Baseline E2 ~ (age 35 y)' = 0.0829). The source reports no statistically significant age effect (95% CI -0.157 to 0.323 spans zero) but retains the term under its full-covariate-model approach, so it is carried here as estimated.",
      source_name        = "AGE"
    ),
    RACE_WHITE = list(
      description        = "Caucasian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (Caucasian; the typical-value reference in this model)",
      notes              = paste(
        "Source paper dichotomises race as Caucasian vs non-Caucasian (Table 3",
        "'Percent Caucasian'). The Caucasian subgroup is the typical-value",
        "reference, so the effect is implemented on (1 - RACE_WHITE): the",
        "Supporting Information baseline-E2 equation writes the term as",
        "theta_7^ASIAN, a multiplicative factor raised to the indicator power.",
        "In this pooled data set the non-Caucasian stratum is entirely the",
        "Japanese KLH1204 cohort (Table 3: 0% Caucasian), which is why the",
        "source's equation labels the indicator ASIAN while its tables label",
        "the same effect non-Caucasian. Same reference-category orientation as",
        "the companion Pohl_2022_linzagolix.R and Hu_2014_bapineuzumab.R."
      ),
      source_name        = "race group (Caucasian / non-Caucasian)"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-volunteer indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (endometriosis patient)",
      notes              = paste(
        "1 = healthy volunteer, 0 = endometriosis patient. Selects BOTH the",
        "estimated baseline oestradiol (Supporting Information section 3.2:",
        "theta_BASE = theta_1 if patient, theta_10 if healthy volunteer) and",
        "the log-scale residual error magnitude (Table 5 'Patients' vs",
        "'Healthy'). The source found approximately 2-fold higher baseline E2",
        "in Caucasian patients than in Caucasian healthy volunteers",
        "(59.1 vs 26.6 pg/mL)."
      ),
      source_name        = "health status (healthy volunteer / patient)"
    ),
    T_ABT = list(
      description        = "Time elapsed on hormonal add-back therapy",
      units              = "weeks",
      type               = "continuous",
      reference_category = "0 (no add-back therapy -- the EDELWEISS phase 2b population and the entire dose-selection analysis)",
      notes              = paste(
        "Weeks since the start of oestrogen/progestogen add-back therapy; 0 for",
        "subjects never receiving it, which zeroes the term. Add-back therapy",
        "(1 mg oestradiol / 0.5 mg norethisterone acetate) was administered only",
        "in the phase 1 studies 16-OBE2109-011 and 17-OBE2109-008 (Table 1).",
        "Table 5 reports 'E2 increase rate on add-back therapy' = 1.58",
        "pg/mL/wk. CAUTION: the functional form of this term is NOT printed",
        "anywhere in the source -- the Supporting Information equation omits it",
        "entirely. It is encoded here as an additive linear-in-time increment,",
        "the only reading consistent with the published pg/mL/wk unit. See the",
        "vignette Assumptions and deviations section."
      ),
      source_name        = "add-back therapy"
    )
  )

  compartmentData <- list(
    e2_placebo = list(analyte = "apparent placebo-arm oestradiol drift (fraction of baseline)", units = "(fraction)", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 724,
    n_studies      = 4,
    age_range      = "18-48 years",
    weight_range   = "median 53.9-65.5 kg across the contributing studies (Table 3)",
    weight_median  = "58 kg (model reference weight)",
    age_median     = "35 years (model reference age); 32-37 years across the contributing studies (Table 3)",
    sex_female_pct = 100,
    race_ethnicity = "0-100% Caucasian by study (Table 3); the non-Caucasian stratum is the Japanese KLH1204 cohort",
    disease_state  = "endometriosis, or endometriosis with co-existing uterine fibroids, plus healthy premenopausal volunteers (approximately 15% of subjects and of observations were from healthy volunteers)",
    dose_range     = "placebo and 25-200 mg once daily for 42 days to 24 weeks",
    regions        = "Europe and Japan",
    notes          = paste(
      "4674 oestradiol observations from 724 subjects (Supporting Information",
      "section 1.2). The PK-E2 analysis set INCLUDED subjects receiving placebo",
      "but EXCLUDED the healthy volunteers of study KLH1101, so it differs from",
      "the population PK analysis set of the companion Pohl_2022_linzagolix.R",
      "(4250 concentrations, 756 subjects, no placebo subjects). 192 E2",
      "observations (4%) were below the limit of quantitation, mostly from",
      "KLH1204 whose quantitation limit was 10 pg/mL against 1 pg/mL in",
      "EDELWEISS; these were imputed at one-half the quantitation limit rather",
      "than dropped. Baseline demographics are in Table 3."
    )
  )

  ini({
    # ---- Estimated baseline oestradiol (Table 5 structural parameters) -------
    # Supporting Information 3.2: theta_BASE,i = theta_1 if patient,
    #                                           theta_10 if healthy volunteer.
    lbl_e2_patient <- log(59.1); label("Typical baseline oestradiol, endometriosis patients (pg/mL)")  # Table 5: 59.1 (95% CI 52.5-65.6)
    lbl_e2_healthy <- log(26.6); label("Typical baseline oestradiol, healthy volunteers (pg/mL)")      # Table 5: 26.6 (95% CI 23.3-29.8)

    # ---- Sigmoid inhibitory Emax on linzagolix daily AUC --------------------
    lauc50 <- log(1.68e5); label("Linzagolix daily AUC producing half-maximal oestradiol suppression (ng*h/mL)")  # Table 5: 1.68e5 (95% CI 1.44e5-1.91e5); Supporting Information 4.2
    lhill  <- log(1.78);   label("Sigmoidicity (Hill) exponent of the AUC-oestradiol relationship (unitless)")    # Table 5: 1.78 (95% CI 1.49-2.08); Supporting Information 4.2 "the sigmoidicity parameter was estimated to be 1.78"

    # Imax is NEVER reported by the source: Table 5 lists every other structural
    # parameter -- including one explicitly marked FIXED -- but no Imax row, and
    # the Supporting Information equation carries a theta_IMAX symbol with no
    # accompanying estimate. Fixed to 1 (complete suppression of ovarian
    # oestradiol production), which is both the mechanistically expected value
    # for a GnRH receptor antagonist and the ONLY value that reproduces the
    # source's own stated result. See the vignette Assumptions and deviations
    # section for the falsification of Imax < 1 and of the no-Hill (gamma = 1)
    # alternative against the paper's published 75-150 mg target-window claim.
    limax <- fixed(log(1)); label("Maximum fractional suppression of oestradiol, Imax (unitless)")

    # ---- Covariate effects on baseline oestradiol (Table 5) ------------------
    e_wt_bl_e2       <- -0.699;  label("Exponent of body weight on baseline oestradiol (58 kg reference, unitless)")   # Table 5: "Baseline E2 ~ (weight 58 kg)" -0.699 (95% CI -0.958 to -0.441)
    e_age_bl_e2      <-  0.0829; label("Exponent of age on baseline oestradiol (35 y reference, unitless)")            # Table 5: "Baseline E2 ~ (age 35 y)" 0.0829 (95% CI -0.157 to 0.323)
    e_nonwhite_bl_e2 <-  0.804;  label("Multiplicative factor on baseline oestradiol for non-Caucasian subjects (unitless)")  # Table 5: "Baseline E2 ~ non-Caucasian" 0.804 (95% CI 0.702-0.907); enters as theta_7^ASIAN
    e_bl_e2_drug     <- -0.120;  label("Exponent relating the observed baseline oestradiol to the linzagolix drug effect (unitless)")  # Table 5: "Baseline E2 ~ linzagolix drug effect" -0.120 (95% CI -0.212 to -0.0279)

    # ---- Apparent placebo-arm oestradiol rise (Table 5) ---------------------
    # Supporting Information 3.2: treatment started during menses, at the nadir
    # of the menstrual E2 cycle, so E2 appears to rise on placebo as cycle
    # synchronisation is lost across the trial.
    lemax_pbo <- log(0.65);          label("Maximum apparent fractional oestradiol increase on placebo (unitless)")  # Table 5: "Placebo increase factor" 0.65 (95% CI 0.465-0.834); Supporting Information 4.2 states the maximum increase is 1.650-fold, i.e. 1 + 0.65, confirming the (1 + PB) form of the equation
    lkp_e2 <- fixed(log(0.231));     label("Rate constant of the apparent placebo oestradiol process (1/week)")      # Table 5: "Placebo effect rate constant" 0.231 FIXED; Supporting Information 3.2 "fixed so that the maximum increase happened at week 12" -- 1-exp(-0.231*12) = 0.94

    # ---- Hormonal add-back therapy (Table 5) --------------------------------
    e_abt_e2 <- 1.58; label("Rate of oestradiol increase on hormonal add-back therapy (pg/mL/week)")  # Table 5: "E2 increase rate on add-back therapy" 1.58 (95% CI 0.990-2.16). Functional form NOT printed by the source; see covariateData$T_ABT

    # ---- Interindividual variability ----------------------------------------
    # Table 5 heading "Interindividual variability (additive on log-scale)";
    # NONMEM $OMEGA convention, so the tabulated number is a VARIANCE.
    etalbl_e2 ~ 0.310  # Table 5: IIV-baseline E2 0.310 (95% CI 0.262-0.358), shrinkage 11.9%

    # ---- Residual error -----------------------------------------------------
    # Supporting Information 3.2: "estimated using a log-transform both sides
    # approach with additive residual error with respect to natural-logarithm
    # transformed E2" -- i.e. lnorm() in nlmixr2. Table 5 heading "Residual
    # variability (additive on log-scale)" follows the NONMEM $SIGMA
    # convention, so each tabulated number is a VARIANCE and the SD below is
    # its square root.
    expSdPatient <- 0.78102; label("Log-scale residual SD, endometriosis patients")  # Table 5: "Patients" variance 0.610 (95% CI 0.571-0.649); sqrt(0.610) = 0.78102
    expSdHealthy <- 0.49092; label("Log-scale residual SD, healthy volunteers")      # Table 5: "Healthy" variance 0.241 (95% CI 0.179-0.303); sqrt(0.241) = 0.49092
  })

  model({
    # 1. Individual baseline oestradiol. The patient / healthy-volunteer switch
    #    is applied on the LOG scale so the model stays mu-referenced in
    #    etalbl_e2 (Supporting Information 3.2 theta_BASE switch, combined with
    #    the baseline-E2 covariate equation).
    lbl_e2 <- lbl_e2_patient * (1 - DIS_HEALTHY) + lbl_e2_healthy * DIS_HEALTHY

    bl_e2 <- exp(lbl_e2 + etalbl_e2) *
      (WT  / 58)^e_wt_bl_e2 *
      (AGE / 35)^e_age_bl_e2 *
      e_nonwhite_bl_e2^(1 - RACE_WHITE)

    # 2. Sigmoid inhibitory Emax drug effect on the individual daily AUC.
    #    The Supporting Information equation prints the inhibition as
    #    theta_IMAX * AUC / (theta_AUC50 + AUC) with NO exponent, yet Table 5
    #    reports a sigmoidicity parameter of 1.78 and both the main text and
    #    Supporting Information 3.2 describe "a sigmoid inhibitory Emax model".
    #    The Hill exponent is therefore restored to its standard position on
    #    both AUC terms; see the vignette for the falsification of the
    #    alternatives against the paper's own published outputs.
    imax   <- exp(limax)
    auc50  <- exp(lauc50)
    hill   <- exp(lhill)
    inhib  <- imax * AUC_LZGX^hill / (auc50^hill + AUC_LZGX^hill)

    # 3. Apparent placebo-arm drift, PB(t) = emax_pbo * (1 - exp(-kp * t)),
    #    encoded as the equivalent first-order approach-to-plateau so the
    #    trajectory is a proper state rather than an explicit function of time.
    #    lkp_e2 is published per WEEK; 168 h per week converts it to the
    #    model's hour time base.
    emax_pbo <- exp(lemax_pbo)
    kp_e2    <- exp(lkp_e2) / 168

    d/dt(e2_placebo) <- kp_e2 * (emax_pbo - e2_placebo)
    e2_placebo(0)    <- 0

    # 4. Observation. Supporting Information 3.2 core equation:
    #      E2_ij = E2_0,i * (1 - Imax * AUC_ij^hill / (AUC50^hill + AUC_ij^hill))
    #                     * (BL_E2_i / 52 pg/mL)^theta_E2BASE
    #                     * (1 + PB_ij)
    #    plus the add-back therapy increment, whose form is inferred from its
    #    published pg/mL/week unit (the source equation omits the term).
    e2 <- bl_e2 * (1 - inhib) *
      (BL_E2 / 52)^e_bl_e2_drug *
      (1 + e2_placebo) +
      e_abt_e2 * T_ABT

    # 5. Log-transform-both-sides residual error, magnitude switched by cohort.
    expSd <- expSdPatient * (1 - DIS_HEALTHY) + expSdHealthy * DIS_HEALTHY

    e2 ~ lnorm(expSd)
  })
}
