Moein_2026_crohns_sescd <- function() {
  description <- paste0(
    "Drug-free (placebo) longitudinal disease-progression model for the ",
    "Simple Endoscopic Score for Crohn's Disease (SES-CD) in adults with ",
    "moderately-to-severely active Crohn's disease, fitted to the placebo ",
    "arm of the phase 3 BERGAMOT trial (Moein 2026, n = 155). Item-level ",
    "item-response-theory bounded-integer (IRT-BI) model: each of the 20 ",
    "SES-CD subscores (size of ulcers, ulcerated surface, affected surface ",
    "and presence of narrowings in the ileum, right colon, transverse ",
    "colon, left colon and rectum; each 0-3) is a bounded-integer ordinal ",
    "endpoint with an item-specific baseline location and dispersion on ",
    "the probit scale, cut at qnorm(1/4), qnorm(2/4) and qnorm(3/4). One ",
    "baseline random effect is shared by the four items of each segment ",
    "and one dispersion random effect by all 20 items. A single latent ",
    "placebo trajectory, (pmax + eta) * (1 - exp(-log(2) / thalf * t)), ",
    "shifts every item's location equally. Includes the authors' ",
    "end-of-induction logistic dropout submodel (probability of not ",
    "entering maintenance as a function of the CDAI change from baseline ",
    "at end of induction), which governs which simulated subjects enter ",
    "the maintenance phase."
  )
  reference <- paste(
    "Moein A, Langenhorst J, Kassir N. A Longitudinal Disease Progression",
    "Model Characterizing Endoscopic Response and Remission in Crohn's",
    "Disease. J Pharmacokinet Pharmacodyn. 2026. doi:10.1007/s10928-026-10042-x.",
    "Final estimates from Table 1; model structure from Supplementary",
    "Material 1 (NONMEM control stream)."
  )
  vignette <- "Moein_2026_crohns_sescd"
  units <- list(
    time = "day (time since first dose; the control stream TSFDD column)",
    dosing = "n/a (no dose events; placebo-arm disease-progression model)",
    concentration = "n/a (SES-CD subscores are ordinal 0-3; escore_* are expected subscores on the 0-3 scale and escore_sescd on the 0-60 scale)"
  )

  covariateData <- list(
    SCORE_CDAI_DELTA = list(
      description = "Change in Crohn's Disease Activity Index (CDAI) from baseline to the end of induction (week 14), CDAI(week 14) - CDAI(baseline). Negative values are improvement.",
      units = "(CDAI points)",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only by the end-of-induction dropout submodel (control stream column CEOICDAI, 'CDAI change at the end of induction'). Time-fixed per subject. Not centred. BERGAMOT patients entered maintenance only with a CDAI decrease of at least 70 points, so the dropout probability is about 0.5 at a change of about -72 points.",
      source_name = "CEOICDAI"
    )
  )

  covariatesDataExcluded <- list(
    CRP = list(
      description = "Baseline C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "FREM covariate on the maximum placebo response. Figure S5 forest plot: relative pmax 0.53 [90% CI -0.19, 0.95] at 75.86 mg/L and 1.5 [1.1, 2.4] at 0.59 mg/L. FREM reports only conditional relative effects, not a parametric covariate relationship, so it is not coded in model()."
    ),
    SCORE_CALPRO = list(
      description = "Baseline fecal calprotectin",
      units = "ug/g",
      type = "continuous",
      notes = "FREM covariate on the maximum placebo response. Figure S5: relative pmax 0.40 [-0.15, 0.89] at 8511 ug/g and 1.5 [1.1, 2.1] at 66.07 ug/g. Not coded in model() (FREM forest plot only)."
    ),
    PRIOR_TNF = list(
      description = "Prior anti-TNF exposure (1 = TNF-experienced, 0 = TNF-naive)",
      units = "(binary)",
      type = "binary",
      notes = "FREM covariate on the maximum placebo response. Figure S5: relative pmax 0.69 [0.41, 0.83] for prior TNF and 1.6 [1.3, 2.1] for TNF-naive. Not coded in model() (FREM forest plot only)."
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened in the FREM; no clinically relevant effect on the maximum placebo response (Figure S5: 1.1 at 102 kg, 0.91 at 49 kg)."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened in the FREM; no clinically relevant effect (Figure S5: 1.0 at 48.98 g/L, 0.96 at 33.11 g/L)."
    ),
    WBC = list(
      description = "Baseline white blood cell count",
      units = "10^9/L",
      type = "continuous",
      notes = "Screened in the FREM; no clinically relevant effect (Figure S5: 1.0 at 14.13, 0.96 at 4.27)."
    ),
    NEUT = list(
      description = "Baseline neutrophil count",
      units = "10^9/L",
      type = "continuous",
      notes = "Screened in the FREM; no clinically relevant effect (Figure S5: 1.1 at 11.48, 0.93 at 2.75)."
    ),
    SCORE_SESCD = list(
      description = "Baseline SES-CD total score",
      units = "(score, 0-56)",
      type = "continuous",
      notes = "Screened in the FREM; no clinically relevant effect (Figure S5: 0.85 at 25.70, 1.2 at 4.68)."
    ),
    SCORE_CDAI = list(
      description = "Baseline CDAI",
      units = "(score, 0-600)",
      type = "continuous",
      notes = "Screened in the FREM; no clinically relevant effect (Figure S5: 0.96 at 457, 1.0 at 234)."
    ),
    SMOKE_CURRENT = list(
      description = "Current smoker (paired with SMOKE_NEVER)",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the FREM; not retained as clinically relevant (Figure S5: current 1.2 [0.84, 1.6], previous 1.3 [1.0, 1.9], never 0.78 [0.50, 0.95])."
    ),
    SMOKE_NEVER = list(
      description = "Never smoker (paired with SMOKE_CURRENT)",
      units = "(binary)",
      type = "binary",
      notes = "See SMOKE_CURRENT."
    ),
    DISLOC_ILEUM = list(
      description = "Ileum-only disease location (reference ileocolonic)",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the FREM; no clinically relevant effect (Figure S5: ileum only 0.93 [0.42, 1.5], ileum and colon 0.99 [0.76, 1.2])."
    ),
    DISLOC_COLON = list(
      description = "Colon-only disease location (reference ileocolonic)",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the FREM; no clinically relevant effect (Figure S5: colon only 1.1 [0.64, 1.5])."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 155L,
    n_studies = 1L,
    n_observations = "6,202 endoscopic subscore observations; endoscopies at baseline, end of induction (week 14) and end of maintenance (week 66)",
    age_range = "18-79 years (mean 36.9, SD 13.2; Table S2)",
    weight_range = "39.7-118.4 kg (mean 70.9, SD 16.2; Table S2)",
    sex_female_pct = 44.5,
    disease_state = "moderately-to-severely active Crohn's disease; baseline CDAI mean 324.9 (range 211-479), baseline SES-CD mean 13.6 (range 4-36); 65.8% TNF-experienced",
    dose_range = "none -- placebo in both the 14-week induction and 52-week maintenance phases",
    regions = "multinational (BERGAMOT, NCT02394028)",
    notes = paste0(
      "Placebo patients of the phase 3 BERGAMOT etrolizumab trial who ",
      "received placebo in induction and, if they achieved a CDAI-70 ",
      "response at week 14, in maintenance (155 in induction, 117 in ",
      "maintenance; Tables S2-S4). Disease location: ileum and colon 62.6%, ",
      "ileum only 18.7%, colon only 18.7%. Fitted in NONMEM 7.4.1 with ",
      "conditional Laplace estimation; covariates on the maximum placebo ",
      "response were evaluated by FREM (Figure S5) but not built into the ",
      "structural model."
    )
  )

  ini({
    # ==================================================================
    # All values are the final estimates of Moein 2026 Table 1 ('Value'
    # column, bounded-integer probit scale). The Supplementary Material 1
    # control stream $THETA / $OMEGA blocks carry the same estimates to
    # within the last printed digit; where they differ the comment names
    # the control-stream value.
    # ==================================================================

    # ----- Ileum (Table 1 'Ileum' block; control stream THETA(1)-THETA(8)) -----
    probitbase_ileum_ulcersize <- -0.050; label("Baseline location of the size of ulcers subscore, ileum, on the bounded-integer probit scale (unitless)")  # Table 1 'Ileum' 'BASE Size of ulcers' = -0.050 (score scale 1.44); THETA(1); the control stream $THETA prints -0.0503
    probitbase_ileum_ulcerated <- -0.185; label("Baseline location of the ulcerated surface subscore, ileum, on the bounded-integer probit scale (unitless)")  # Table 1 'Ileum' 'BASE Ulcerated surface' = -0.185 (score scale 1.28); THETA(3)
    probitbase_ileum_affected <- -0.0765; label("Baseline location of the affected surface subscore, ileum, on the bounded-integer probit scale (unitless)")  # Table 1 'Ileum' 'BASE Affected surface' = -0.0765 (score scale 1.409); THETA(5); the control stream $THETA prints -0.0767
    probitbase_ileum_narrowing <- -1.26; label("Baseline location of the presence of narrowings subscore, ileum, on the bounded-integer probit scale (unitless)")  # Table 1 'Ileum' 'BASE Presence of narrowings' = -1.26 (score scale 0.312); THETA(7)
    lsd_ileum_ulcersize <- log(0.400); label("Dispersion (SD) of the size of ulcers subscore, ileum, on the bounded-integer probit scale (unitless)")  # Table 1 'Ileum' 'SD Size of ulcers' = 0.400; THETA(2)
    lsd_ileum_ulcerated <- log(0.290); label("Dispersion (SD) of the ulcerated surface subscore, ileum, on the bounded-integer probit scale (unitless)")  # Table 1 'Ileum' 'SD Ulcerated surface' = 0.290; THETA(4)
    lsd_ileum_affected <- log(0.388); label("Dispersion (SD) of the affected surface subscore, ileum, on the bounded-integer probit scale (unitless)")  # Table 1 'Ileum' 'SD Affected surface' = 0.388; THETA(6)
    lsd_ileum_narrowing <- log(1.260); label("Dispersion (SD) of the presence of narrowings subscore, ileum, on the bounded-integer probit scale (unitless)")  # Table 1 'Ileum' 'SD Presence of narrowings' = 1.260; THETA(8)

    # ----- Right colon (Table 1 'Right colon' block; control stream THETA(9)-THETA(16)) -----
    probitbase_rcolon_ulcersize <- -0.545; label("Baseline location of the size of ulcers subscore, right colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Right colon' 'BASE Size of ulcers' = -0.545 (score scale 0.879); THETA(9)
    probitbase_rcolon_ulcerated <- -0.689; label("Baseline location of the ulcerated surface subscore, right colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Right colon' 'BASE Ulcerated surface' = -0.689 (score scale 0.736); THETA(11)
    probitbase_rcolon_affected <- -0.602; label("Baseline location of the affected surface subscore, right colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Right colon' 'BASE Affected surface' = -0.602 (score scale 0.821); THETA(13); the control stream $THETA prints -0.603
    probitbase_rcolon_narrowing <- -2.53; label("Baseline location of the presence of narrowings subscore, right colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Right colon' 'BASE Presence of narrowings' = -2.53 (score scale 0.017); THETA(15)
    lsd_rcolon_ulcersize <- log(0.372); label("Dispersion (SD) of the size of ulcers subscore, right colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Right colon' 'SD Size of ulcers' = 0.372; THETA(10)
    lsd_rcolon_ulcerated <- log(0.206); label("Dispersion (SD) of the ulcerated surface subscore, right colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Right colon' 'SD Ulcerated surface' = 0.206; THETA(12)
    lsd_rcolon_affected <- log(0.316); label("Dispersion (SD) of the affected surface subscore, right colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Right colon' 'SD Affected surface' = 0.316; THETA(14)
    lsd_rcolon_narrowing <- log(1.47); label("Dispersion (SD) of the presence of narrowings subscore, right colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Right colon' 'SD Presence of narrowings' = 1.47; THETA(16)

    # ----- Transverse colon (Table 1 'Transverse colon' block; control stream THETA(17)-THETA(24)) -----
    probitbase_tcolon_ulcersize <- -0.965; label("Baseline location of the size of ulcers subscore, transverse colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Transverse colon' 'BASE Size of ulcers' = -0.965 (score scale 0.503); THETA(17); the control stream $THETA prints -0.964
    probitbase_tcolon_ulcerated <- -1.09; label("Baseline location of the ulcerated surface subscore, transverse colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Transverse colon' 'BASE Ulcerated surface' = -1.09 (score scale 0.414); THETA(19)
    probitbase_tcolon_affected <- -0.931; label("Baseline location of the affected surface subscore, transverse colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Transverse colon' 'BASE Affected surface' = -0.931 (score scale 0.529); THETA(21); the control stream $THETA prints -0.93
    probitbase_tcolon_narrowing <- -3.22; label("Baseline location of the presence of narrowings subscore, transverse colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Transverse colon' 'BASE Presence of narrowings' = -3.22 (score scale 0.002); THETA(23)
    lsd_tcolon_ulcersize <- log(0.349); label("Dispersion (SD) of the size of ulcers subscore, transverse colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Transverse colon' 'SD Size of ulcers' = 0.349; THETA(18)
    lsd_tcolon_ulcerated <- log(0.237); label("Dispersion (SD) of the ulcerated surface subscore, transverse colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Transverse colon' 'SD Ulcerated surface' = 0.237; THETA(20)
    lsd_tcolon_affected <- log(0.481); label("Dispersion (SD) of the affected surface subscore, transverse colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Transverse colon' 'SD Affected surface' = 0.481; THETA(22)
    lsd_tcolon_narrowing <- log(1.10); label("Dispersion (SD) of the presence of narrowings subscore, transverse colon, on the bounded-integer probit scale (unitless)")  # Table 1 'Transverse colon' 'SD Presence of narrowings' = 1.10; THETA(24)

    # ----- Left colon (Table 1 'Left colon' block; control stream THETA(25)-THETA(32)) -----
    probitbase_lcolon_ulcersize <- -0.493; label("Baseline location of the size of ulcers subscore, left colon and sigmoid, on the bounded-integer probit scale (unitless)")  # Table 1 'Left colon' 'BASE Size of ulcers' = -0.493 (score scale 0.933); THETA(25)
    probitbase_lcolon_ulcerated <- -0.672; label("Baseline location of the ulcerated surface subscore, left colon and sigmoid, on the bounded-integer probit scale (unitless)")  # Table 1 'Left colon' 'BASE Ulcerated surface' = -0.672 (score scale 0.752); THETA(27)
    probitbase_lcolon_affected <- -0.515; label("Baseline location of the affected surface subscore, left colon and sigmoid, on the bounded-integer probit scale (unitless)")  # Table 1 'Left colon' 'BASE Affected surface' = -0.515 (score scale 0.910); THETA(29)
    probitbase_lcolon_narrowing <- -2.39; label("Baseline location of the presence of narrowings subscore, left colon and sigmoid, on the bounded-integer probit scale (unitless)")  # Table 1 'Left colon' 'BASE Presence of narrowings' = -2.39 (score scale 0.025); THETA(31)
    lsd_lcolon_ulcersize <- log(0.331); label("Dispersion (SD) of the size of ulcers subscore, left colon and sigmoid, on the bounded-integer probit scale (unitless)")  # Table 1 'Left colon' 'SD Size of ulcers' = 0.331; THETA(26)
    lsd_lcolon_ulcerated <- log(0.252); label("Dispersion (SD) of the ulcerated surface subscore, left colon and sigmoid, on the bounded-integer probit scale (unitless)")  # Table 1 'Left colon' 'SD Ulcerated surface' = 0.252; THETA(28)
    lsd_lcolon_affected <- log(0.375); label("Dispersion (SD) of the affected surface subscore, left colon and sigmoid, on the bounded-integer probit scale (unitless)")  # Table 1 'Left colon' 'SD Affected surface' = 0.375; THETA(30)
    lsd_lcolon_narrowing <- log(0.92); label("Dispersion (SD) of the presence of narrowings subscore, left colon and sigmoid, on the bounded-integer probit scale (unitless)")  # Table 1 'Left colon' 'SD Presence of narrowings' = 0.92; THETA(32)

    # ----- Rectum (Table 1 'Rectum' block; control stream THETA(33)-THETA(40)) -----
    probitbase_rectum_ulcersize <- -0.847; label("Baseline location of the size of ulcers subscore, rectum, on the bounded-integer probit scale (unitless)")  # Table 1 'Rectum' 'BASE Size of ulcers' = -0.847 (score scale 0.595); THETA(33)
    probitbase_rectum_ulcerated <- -0.932; label("Baseline location of the ulcerated surface subscore, rectum, on the bounded-integer probit scale (unitless)")  # Table 1 'Rectum' 'BASE Ulcerated surface' = -0.932 (score scale 0.527); THETA(35)
    probitbase_rectum_affected <- -0.759; label("Baseline location of the affected surface subscore, rectum, on the bounded-integer probit scale (unitless)")  # Table 1 'Rectum' 'BASE Affected surface' = -0.759 (score scale 0.672); THETA(37)
    probitbase_rectum_narrowing <- -2.85; label("Baseline location of the presence of narrowings subscore, rectum, on the bounded-integer probit scale (unitless)")  # Table 1 'Rectum' 'BASE Presence of narrowings' = -2.85 (score scale 0.007); THETA(39)
    lsd_rectum_ulcersize <- log(0.38); label("Dispersion (SD) of the size of ulcers subscore, rectum, on the bounded-integer probit scale (unitless)")  # Table 1 'Rectum' 'SD Size of ulcers' = 0.38; THETA(34)
    lsd_rectum_ulcerated <- log(0.257); label("Dispersion (SD) of the ulcerated surface subscore, rectum, on the bounded-integer probit scale (unitless)")  # Table 1 'Rectum' 'SD Ulcerated surface' = 0.257; THETA(36)
    lsd_rectum_affected <- log(0.398); label("Dispersion (SD) of the affected surface subscore, rectum, on the bounded-integer probit scale (unitless)")  # Table 1 'Rectum' 'SD Affected surface' = 0.398; THETA(38)
    lsd_rectum_narrowing <- log(0.919); label("Dispersion (SD) of the presence of narrowings subscore, rectum, on the bounded-integer probit scale (unitless)")  # Table 1 'Rectum' 'SD Presence of narrowings' = 0.919; THETA(40)
    # ----- Latent placebo trajectory (Equation 1) -----
    pmax <- -0.329; label("Maximum placebo response RMaxPL, a shift of every item location on the probit scale (unitless)")  # Table 1 'R MaxPL' = -0.329 (90% CI -0.487, -0.179; score scale -0.220); THETA(41). Linear scale: the value is negative and its random effect is additive
    lthalf_pbo <- log(13.8); label("Time to half of the maximum placebo response T50PL (week)")  # Table 1 'T50 PL' = 13.8 weeks (90% CI 11.32, 16.21); THETA(42)

    # ----- End-of-induction dropout (logistic regression) -----
    logitp_dropout <- 4.36; label("Logit of the probability of not entering maintenance at a CDAI change of 0 from baseline (unitless)")  # Table 1 'Intercept Dropout End of Induction' = 4.36 (90% CI 3.661, 5.251); THETA(43)
    e_score_cdai_delta_p_dropout <- 0.0607; label("Change in the dropout logit per CDAI point of change from baseline at end of induction (unitless per CDAI point)")  # Table 1 'Slope Dropout End of Induction' = 0.0607 (90% CI 0.050, 0.073); THETA(44)

    # ----- Interindividual variability (variances; NONMEM $OMEGA) -----
    # One baseline random effect per segment, shared by its four items.
    etaprobitbase_ileum ~ 0.712  # Table 1 'IIVBASE, Ileum' = 0.712 (shrinkage 16.5%); $OMEGA 0.713
    etaprobitbase_rcolon ~ 0.521 # Table 1 'IIVBASE, Right colon' = 0.521 (shrinkage 20.7%); $OMEGA 0.522
    etaprobitbase_tcolon ~ 1.07  # Table 1 'IIVBASE, Transverse colon' = 1.07 (shrinkage 25.2%); $OMEGA 1.07
    etaprobitbase_lcolon ~ 0.923 # Table 1 'IIVBASE, Left colon' = 0.923 (shrinkage 18.3%); $OMEGA 0.924
    etaprobitbase_rectum ~ 1.11  # Table 1 'IIVBASE, Rectum' = 1.11 (shrinkage 22.2%); $OMEGA 1.11
    # One exponential dispersion random effect shared by all 20 items.
    etalsd ~ 0.0923              # Table 1 'IIVSD IRT' = 0.0923 (shrinkage 14.6%); $OMEGA 0.0922, ETA(6)
    # Additive random effect on the maximum placebo response.
    etapmax ~ 0.648              # Table 1 'IIV RMaxPL' = 0.648 (shrinkage 18.9%); $OMEGA 0.648, ETA(7)
  })

  model({
    # Latent placebo response (Equation 1). The control stream evaluates
    # it in weeks, (TSFDD / 7); time here is in days.
    thalf_pbo <- exp(lthalf_pbo)
    pbo <- (pmax + etapmax) * (1 - exp(-log(2) / thalf_pbo * (time / 7)))

    # Bounded-integer cut points: qnorm(1:3 / 4), the control stream CO1-CO3.
    cut1 <- -0.6744897502
    cut2 <- 0
    cut3 <- 0.6744897502

    # Each item (Equation 2): location = baseline + segment eta + placebo
    # response, dispersion = item SD * exp(shared eta); per-score
    # probabilities from the standard normal CDF; escore_* is the
    # expected subscore (0-3).
    # Ileum: all four items share etaprobitbase_ileum (control stream ETA(1))
    mean_ileum_ulcersize <- probitbase_ileum_ulcersize + etaprobitbase_ileum + pbo
    sd_ileum_ulcersize <- exp(lsd_ileum_ulcersize + etalsd)
    p0_ileum_ulcersize <- phi((cut1 - mean_ileum_ulcersize) / sd_ileum_ulcersize)
    p1_ileum_ulcersize <- phi((cut2 - mean_ileum_ulcersize) / sd_ileum_ulcersize) - phi((cut1 - mean_ileum_ulcersize) / sd_ileum_ulcersize)
    p2_ileum_ulcersize <- phi((cut3 - mean_ileum_ulcersize) / sd_ileum_ulcersize) - phi((cut2 - mean_ileum_ulcersize) / sd_ileum_ulcersize)
    p3_ileum_ulcersize <- 1 - phi((cut3 - mean_ileum_ulcersize) / sd_ileum_ulcersize)
    escore_ileum_ulcersize <- p1_ileum_ulcersize + 2 * p2_ileum_ulcersize + 3 * p3_ileum_ulcersize

    mean_ileum_ulcerated <- probitbase_ileum_ulcerated + etaprobitbase_ileum + pbo
    sd_ileum_ulcerated <- exp(lsd_ileum_ulcerated + etalsd)
    p0_ileum_ulcerated <- phi((cut1 - mean_ileum_ulcerated) / sd_ileum_ulcerated)
    p1_ileum_ulcerated <- phi((cut2 - mean_ileum_ulcerated) / sd_ileum_ulcerated) - phi((cut1 - mean_ileum_ulcerated) / sd_ileum_ulcerated)
    p2_ileum_ulcerated <- phi((cut3 - mean_ileum_ulcerated) / sd_ileum_ulcerated) - phi((cut2 - mean_ileum_ulcerated) / sd_ileum_ulcerated)
    p3_ileum_ulcerated <- 1 - phi((cut3 - mean_ileum_ulcerated) / sd_ileum_ulcerated)
    escore_ileum_ulcerated <- p1_ileum_ulcerated + 2 * p2_ileum_ulcerated + 3 * p3_ileum_ulcerated

    mean_ileum_affected <- probitbase_ileum_affected + etaprobitbase_ileum + pbo
    sd_ileum_affected <- exp(lsd_ileum_affected + etalsd)
    p0_ileum_affected <- phi((cut1 - mean_ileum_affected) / sd_ileum_affected)
    p1_ileum_affected <- phi((cut2 - mean_ileum_affected) / sd_ileum_affected) - phi((cut1 - mean_ileum_affected) / sd_ileum_affected)
    p2_ileum_affected <- phi((cut3 - mean_ileum_affected) / sd_ileum_affected) - phi((cut2 - mean_ileum_affected) / sd_ileum_affected)
    p3_ileum_affected <- 1 - phi((cut3 - mean_ileum_affected) / sd_ileum_affected)
    escore_ileum_affected <- p1_ileum_affected + 2 * p2_ileum_affected + 3 * p3_ileum_affected

    mean_ileum_narrowing <- probitbase_ileum_narrowing + etaprobitbase_ileum + pbo
    sd_ileum_narrowing <- exp(lsd_ileum_narrowing + etalsd)
    p0_ileum_narrowing <- phi((cut1 - mean_ileum_narrowing) / sd_ileum_narrowing)
    p1_ileum_narrowing <- phi((cut2 - mean_ileum_narrowing) / sd_ileum_narrowing) - phi((cut1 - mean_ileum_narrowing) / sd_ileum_narrowing)
    p2_ileum_narrowing <- phi((cut3 - mean_ileum_narrowing) / sd_ileum_narrowing) - phi((cut2 - mean_ileum_narrowing) / sd_ileum_narrowing)
    p3_ileum_narrowing <- 1 - phi((cut3 - mean_ileum_narrowing) / sd_ileum_narrowing)
    escore_ileum_narrowing <- p1_ileum_narrowing + 2 * p2_ileum_narrowing + 3 * p3_ileum_narrowing

    # Right colon: all four items share etaprobitbase_rcolon (control stream ETA(2))
    mean_rcolon_ulcersize <- probitbase_rcolon_ulcersize + etaprobitbase_rcolon + pbo
    sd_rcolon_ulcersize <- exp(lsd_rcolon_ulcersize + etalsd)
    p0_rcolon_ulcersize <- phi((cut1 - mean_rcolon_ulcersize) / sd_rcolon_ulcersize)
    p1_rcolon_ulcersize <- phi((cut2 - mean_rcolon_ulcersize) / sd_rcolon_ulcersize) - phi((cut1 - mean_rcolon_ulcersize) / sd_rcolon_ulcersize)
    p2_rcolon_ulcersize <- phi((cut3 - mean_rcolon_ulcersize) / sd_rcolon_ulcersize) - phi((cut2 - mean_rcolon_ulcersize) / sd_rcolon_ulcersize)
    p3_rcolon_ulcersize <- 1 - phi((cut3 - mean_rcolon_ulcersize) / sd_rcolon_ulcersize)
    escore_rcolon_ulcersize <- p1_rcolon_ulcersize + 2 * p2_rcolon_ulcersize + 3 * p3_rcolon_ulcersize

    mean_rcolon_ulcerated <- probitbase_rcolon_ulcerated + etaprobitbase_rcolon + pbo
    sd_rcolon_ulcerated <- exp(lsd_rcolon_ulcerated + etalsd)
    p0_rcolon_ulcerated <- phi((cut1 - mean_rcolon_ulcerated) / sd_rcolon_ulcerated)
    p1_rcolon_ulcerated <- phi((cut2 - mean_rcolon_ulcerated) / sd_rcolon_ulcerated) - phi((cut1 - mean_rcolon_ulcerated) / sd_rcolon_ulcerated)
    p2_rcolon_ulcerated <- phi((cut3 - mean_rcolon_ulcerated) / sd_rcolon_ulcerated) - phi((cut2 - mean_rcolon_ulcerated) / sd_rcolon_ulcerated)
    p3_rcolon_ulcerated <- 1 - phi((cut3 - mean_rcolon_ulcerated) / sd_rcolon_ulcerated)
    escore_rcolon_ulcerated <- p1_rcolon_ulcerated + 2 * p2_rcolon_ulcerated + 3 * p3_rcolon_ulcerated

    mean_rcolon_affected <- probitbase_rcolon_affected + etaprobitbase_rcolon + pbo
    sd_rcolon_affected <- exp(lsd_rcolon_affected + etalsd)
    p0_rcolon_affected <- phi((cut1 - mean_rcolon_affected) / sd_rcolon_affected)
    p1_rcolon_affected <- phi((cut2 - mean_rcolon_affected) / sd_rcolon_affected) - phi((cut1 - mean_rcolon_affected) / sd_rcolon_affected)
    p2_rcolon_affected <- phi((cut3 - mean_rcolon_affected) / sd_rcolon_affected) - phi((cut2 - mean_rcolon_affected) / sd_rcolon_affected)
    p3_rcolon_affected <- 1 - phi((cut3 - mean_rcolon_affected) / sd_rcolon_affected)
    escore_rcolon_affected <- p1_rcolon_affected + 2 * p2_rcolon_affected + 3 * p3_rcolon_affected

    mean_rcolon_narrowing <- probitbase_rcolon_narrowing + etaprobitbase_rcolon + pbo
    sd_rcolon_narrowing <- exp(lsd_rcolon_narrowing + etalsd)
    p0_rcolon_narrowing <- phi((cut1 - mean_rcolon_narrowing) / sd_rcolon_narrowing)
    p1_rcolon_narrowing <- phi((cut2 - mean_rcolon_narrowing) / sd_rcolon_narrowing) - phi((cut1 - mean_rcolon_narrowing) / sd_rcolon_narrowing)
    p2_rcolon_narrowing <- phi((cut3 - mean_rcolon_narrowing) / sd_rcolon_narrowing) - phi((cut2 - mean_rcolon_narrowing) / sd_rcolon_narrowing)
    p3_rcolon_narrowing <- 1 - phi((cut3 - mean_rcolon_narrowing) / sd_rcolon_narrowing)
    escore_rcolon_narrowing <- p1_rcolon_narrowing + 2 * p2_rcolon_narrowing + 3 * p3_rcolon_narrowing

    # Transverse colon: all four items share etaprobitbase_tcolon (control stream ETA(3))
    mean_tcolon_ulcersize <- probitbase_tcolon_ulcersize + etaprobitbase_tcolon + pbo
    sd_tcolon_ulcersize <- exp(lsd_tcolon_ulcersize + etalsd)
    p0_tcolon_ulcersize <- phi((cut1 - mean_tcolon_ulcersize) / sd_tcolon_ulcersize)
    p1_tcolon_ulcersize <- phi((cut2 - mean_tcolon_ulcersize) / sd_tcolon_ulcersize) - phi((cut1 - mean_tcolon_ulcersize) / sd_tcolon_ulcersize)
    p2_tcolon_ulcersize <- phi((cut3 - mean_tcolon_ulcersize) / sd_tcolon_ulcersize) - phi((cut2 - mean_tcolon_ulcersize) / sd_tcolon_ulcersize)
    p3_tcolon_ulcersize <- 1 - phi((cut3 - mean_tcolon_ulcersize) / sd_tcolon_ulcersize)
    escore_tcolon_ulcersize <- p1_tcolon_ulcersize + 2 * p2_tcolon_ulcersize + 3 * p3_tcolon_ulcersize

    mean_tcolon_ulcerated <- probitbase_tcolon_ulcerated + etaprobitbase_tcolon + pbo
    sd_tcolon_ulcerated <- exp(lsd_tcolon_ulcerated + etalsd)
    p0_tcolon_ulcerated <- phi((cut1 - mean_tcolon_ulcerated) / sd_tcolon_ulcerated)
    p1_tcolon_ulcerated <- phi((cut2 - mean_tcolon_ulcerated) / sd_tcolon_ulcerated) - phi((cut1 - mean_tcolon_ulcerated) / sd_tcolon_ulcerated)
    p2_tcolon_ulcerated <- phi((cut3 - mean_tcolon_ulcerated) / sd_tcolon_ulcerated) - phi((cut2 - mean_tcolon_ulcerated) / sd_tcolon_ulcerated)
    p3_tcolon_ulcerated <- 1 - phi((cut3 - mean_tcolon_ulcerated) / sd_tcolon_ulcerated)
    escore_tcolon_ulcerated <- p1_tcolon_ulcerated + 2 * p2_tcolon_ulcerated + 3 * p3_tcolon_ulcerated

    mean_tcolon_affected <- probitbase_tcolon_affected + etaprobitbase_tcolon + pbo
    sd_tcolon_affected <- exp(lsd_tcolon_affected + etalsd)
    p0_tcolon_affected <- phi((cut1 - mean_tcolon_affected) / sd_tcolon_affected)
    p1_tcolon_affected <- phi((cut2 - mean_tcolon_affected) / sd_tcolon_affected) - phi((cut1 - mean_tcolon_affected) / sd_tcolon_affected)
    p2_tcolon_affected <- phi((cut3 - mean_tcolon_affected) / sd_tcolon_affected) - phi((cut2 - mean_tcolon_affected) / sd_tcolon_affected)
    p3_tcolon_affected <- 1 - phi((cut3 - mean_tcolon_affected) / sd_tcolon_affected)
    escore_tcolon_affected <- p1_tcolon_affected + 2 * p2_tcolon_affected + 3 * p3_tcolon_affected

    mean_tcolon_narrowing <- probitbase_tcolon_narrowing + etaprobitbase_tcolon + pbo
    sd_tcolon_narrowing <- exp(lsd_tcolon_narrowing + etalsd)
    p0_tcolon_narrowing <- phi((cut1 - mean_tcolon_narrowing) / sd_tcolon_narrowing)
    p1_tcolon_narrowing <- phi((cut2 - mean_tcolon_narrowing) / sd_tcolon_narrowing) - phi((cut1 - mean_tcolon_narrowing) / sd_tcolon_narrowing)
    p2_tcolon_narrowing <- phi((cut3 - mean_tcolon_narrowing) / sd_tcolon_narrowing) - phi((cut2 - mean_tcolon_narrowing) / sd_tcolon_narrowing)
    p3_tcolon_narrowing <- 1 - phi((cut3 - mean_tcolon_narrowing) / sd_tcolon_narrowing)
    escore_tcolon_narrowing <- p1_tcolon_narrowing + 2 * p2_tcolon_narrowing + 3 * p3_tcolon_narrowing

    # Left colon: all four items share etaprobitbase_lcolon (control stream ETA(4))
    mean_lcolon_ulcersize <- probitbase_lcolon_ulcersize + etaprobitbase_lcolon + pbo
    sd_lcolon_ulcersize <- exp(lsd_lcolon_ulcersize + etalsd)
    p0_lcolon_ulcersize <- phi((cut1 - mean_lcolon_ulcersize) / sd_lcolon_ulcersize)
    p1_lcolon_ulcersize <- phi((cut2 - mean_lcolon_ulcersize) / sd_lcolon_ulcersize) - phi((cut1 - mean_lcolon_ulcersize) / sd_lcolon_ulcersize)
    p2_lcolon_ulcersize <- phi((cut3 - mean_lcolon_ulcersize) / sd_lcolon_ulcersize) - phi((cut2 - mean_lcolon_ulcersize) / sd_lcolon_ulcersize)
    p3_lcolon_ulcersize <- 1 - phi((cut3 - mean_lcolon_ulcersize) / sd_lcolon_ulcersize)
    escore_lcolon_ulcersize <- p1_lcolon_ulcersize + 2 * p2_lcolon_ulcersize + 3 * p3_lcolon_ulcersize

    mean_lcolon_ulcerated <- probitbase_lcolon_ulcerated + etaprobitbase_lcolon + pbo
    sd_lcolon_ulcerated <- exp(lsd_lcolon_ulcerated + etalsd)
    p0_lcolon_ulcerated <- phi((cut1 - mean_lcolon_ulcerated) / sd_lcolon_ulcerated)
    p1_lcolon_ulcerated <- phi((cut2 - mean_lcolon_ulcerated) / sd_lcolon_ulcerated) - phi((cut1 - mean_lcolon_ulcerated) / sd_lcolon_ulcerated)
    p2_lcolon_ulcerated <- phi((cut3 - mean_lcolon_ulcerated) / sd_lcolon_ulcerated) - phi((cut2 - mean_lcolon_ulcerated) / sd_lcolon_ulcerated)
    p3_lcolon_ulcerated <- 1 - phi((cut3 - mean_lcolon_ulcerated) / sd_lcolon_ulcerated)
    escore_lcolon_ulcerated <- p1_lcolon_ulcerated + 2 * p2_lcolon_ulcerated + 3 * p3_lcolon_ulcerated

    mean_lcolon_affected <- probitbase_lcolon_affected + etaprobitbase_lcolon + pbo
    sd_lcolon_affected <- exp(lsd_lcolon_affected + etalsd)
    p0_lcolon_affected <- phi((cut1 - mean_lcolon_affected) / sd_lcolon_affected)
    p1_lcolon_affected <- phi((cut2 - mean_lcolon_affected) / sd_lcolon_affected) - phi((cut1 - mean_lcolon_affected) / sd_lcolon_affected)
    p2_lcolon_affected <- phi((cut3 - mean_lcolon_affected) / sd_lcolon_affected) - phi((cut2 - mean_lcolon_affected) / sd_lcolon_affected)
    p3_lcolon_affected <- 1 - phi((cut3 - mean_lcolon_affected) / sd_lcolon_affected)
    escore_lcolon_affected <- p1_lcolon_affected + 2 * p2_lcolon_affected + 3 * p3_lcolon_affected

    mean_lcolon_narrowing <- probitbase_lcolon_narrowing + etaprobitbase_lcolon + pbo
    sd_lcolon_narrowing <- exp(lsd_lcolon_narrowing + etalsd)
    p0_lcolon_narrowing <- phi((cut1 - mean_lcolon_narrowing) / sd_lcolon_narrowing)
    p1_lcolon_narrowing <- phi((cut2 - mean_lcolon_narrowing) / sd_lcolon_narrowing) - phi((cut1 - mean_lcolon_narrowing) / sd_lcolon_narrowing)
    p2_lcolon_narrowing <- phi((cut3 - mean_lcolon_narrowing) / sd_lcolon_narrowing) - phi((cut2 - mean_lcolon_narrowing) / sd_lcolon_narrowing)
    p3_lcolon_narrowing <- 1 - phi((cut3 - mean_lcolon_narrowing) / sd_lcolon_narrowing)
    escore_lcolon_narrowing <- p1_lcolon_narrowing + 2 * p2_lcolon_narrowing + 3 * p3_lcolon_narrowing

    # Rectum: all four items share etaprobitbase_rectum (control stream ETA(5))
    mean_rectum_ulcersize <- probitbase_rectum_ulcersize + etaprobitbase_rectum + pbo
    sd_rectum_ulcersize <- exp(lsd_rectum_ulcersize + etalsd)
    p0_rectum_ulcersize <- phi((cut1 - mean_rectum_ulcersize) / sd_rectum_ulcersize)
    p1_rectum_ulcersize <- phi((cut2 - mean_rectum_ulcersize) / sd_rectum_ulcersize) - phi((cut1 - mean_rectum_ulcersize) / sd_rectum_ulcersize)
    p2_rectum_ulcersize <- phi((cut3 - mean_rectum_ulcersize) / sd_rectum_ulcersize) - phi((cut2 - mean_rectum_ulcersize) / sd_rectum_ulcersize)
    p3_rectum_ulcersize <- 1 - phi((cut3 - mean_rectum_ulcersize) / sd_rectum_ulcersize)
    escore_rectum_ulcersize <- p1_rectum_ulcersize + 2 * p2_rectum_ulcersize + 3 * p3_rectum_ulcersize

    mean_rectum_ulcerated <- probitbase_rectum_ulcerated + etaprobitbase_rectum + pbo
    sd_rectum_ulcerated <- exp(lsd_rectum_ulcerated + etalsd)
    p0_rectum_ulcerated <- phi((cut1 - mean_rectum_ulcerated) / sd_rectum_ulcerated)
    p1_rectum_ulcerated <- phi((cut2 - mean_rectum_ulcerated) / sd_rectum_ulcerated) - phi((cut1 - mean_rectum_ulcerated) / sd_rectum_ulcerated)
    p2_rectum_ulcerated <- phi((cut3 - mean_rectum_ulcerated) / sd_rectum_ulcerated) - phi((cut2 - mean_rectum_ulcerated) / sd_rectum_ulcerated)
    p3_rectum_ulcerated <- 1 - phi((cut3 - mean_rectum_ulcerated) / sd_rectum_ulcerated)
    escore_rectum_ulcerated <- p1_rectum_ulcerated + 2 * p2_rectum_ulcerated + 3 * p3_rectum_ulcerated

    mean_rectum_affected <- probitbase_rectum_affected + etaprobitbase_rectum + pbo
    sd_rectum_affected <- exp(lsd_rectum_affected + etalsd)
    p0_rectum_affected <- phi((cut1 - mean_rectum_affected) / sd_rectum_affected)
    p1_rectum_affected <- phi((cut2 - mean_rectum_affected) / sd_rectum_affected) - phi((cut1 - mean_rectum_affected) / sd_rectum_affected)
    p2_rectum_affected <- phi((cut3 - mean_rectum_affected) / sd_rectum_affected) - phi((cut2 - mean_rectum_affected) / sd_rectum_affected)
    p3_rectum_affected <- 1 - phi((cut3 - mean_rectum_affected) / sd_rectum_affected)
    escore_rectum_affected <- p1_rectum_affected + 2 * p2_rectum_affected + 3 * p3_rectum_affected

    mean_rectum_narrowing <- probitbase_rectum_narrowing + etaprobitbase_rectum + pbo
    sd_rectum_narrowing <- exp(lsd_rectum_narrowing + etalsd)
    p0_rectum_narrowing <- phi((cut1 - mean_rectum_narrowing) / sd_rectum_narrowing)
    p1_rectum_narrowing <- phi((cut2 - mean_rectum_narrowing) / sd_rectum_narrowing) - phi((cut1 - mean_rectum_narrowing) / sd_rectum_narrowing)
    p2_rectum_narrowing <- phi((cut3 - mean_rectum_narrowing) / sd_rectum_narrowing) - phi((cut2 - mean_rectum_narrowing) / sd_rectum_narrowing)
    p3_rectum_narrowing <- 1 - phi((cut3 - mean_rectum_narrowing) / sd_rectum_narrowing)
    escore_rectum_narrowing <- p1_rectum_narrowing + 2 * p2_rectum_narrowing + 3 * p3_rectum_narrowing
    # Expected total SES-CD (0-60 by this 20-item construction).
    escore_sescd <- escore_ileum_ulcersize +
      escore_ileum_ulcerated +
      escore_ileum_affected +
      escore_ileum_narrowing +
      escore_rcolon_ulcersize +
      escore_rcolon_ulcerated +
      escore_rcolon_affected +
      escore_rcolon_narrowing +
      escore_tcolon_ulcersize +
      escore_tcolon_ulcerated +
      escore_tcolon_affected +
      escore_tcolon_narrowing +
      escore_lcolon_ulcersize +
      escore_lcolon_ulcerated +
      escore_lcolon_affected +
      escore_lcolon_narrowing +
      escore_rectum_ulcersize +
      escore_rectum_ulcerated +
      escore_rectum_affected +
      escore_rectum_narrowing

    # End-of-induction dropout: probability that a subject does not
    # enter the maintenance phase, from the CDAI change at week 14.
    prob_dropout <- expit(logitp_dropout + e_score_cdai_delta_p_dropout * SCORE_CDAI_DELTA)
    prob_continue <- 1 - prob_dropout

    # Ordinal endpoints. rxode2 numbers the categories 1-4; category k is
    # subscore k - 1 and carries p0, p1, p2 and 1 - p0 - p1 - p2.
    sescd_ileum_ulcersize ~ c(p0_ileum_ulcersize, p1_ileum_ulcersize, p2_ileum_ulcersize)
    sescd_ileum_ulcerated ~ c(p0_ileum_ulcerated, p1_ileum_ulcerated, p2_ileum_ulcerated)
    sescd_ileum_affected ~ c(p0_ileum_affected, p1_ileum_affected, p2_ileum_affected)
    sescd_ileum_narrowing ~ c(p0_ileum_narrowing, p1_ileum_narrowing, p2_ileum_narrowing)
    sescd_rcolon_ulcersize ~ c(p0_rcolon_ulcersize, p1_rcolon_ulcersize, p2_rcolon_ulcersize)
    sescd_rcolon_ulcerated ~ c(p0_rcolon_ulcerated, p1_rcolon_ulcerated, p2_rcolon_ulcerated)
    sescd_rcolon_affected ~ c(p0_rcolon_affected, p1_rcolon_affected, p2_rcolon_affected)
    sescd_rcolon_narrowing ~ c(p0_rcolon_narrowing, p1_rcolon_narrowing, p2_rcolon_narrowing)
    sescd_tcolon_ulcersize ~ c(p0_tcolon_ulcersize, p1_tcolon_ulcersize, p2_tcolon_ulcersize)
    sescd_tcolon_ulcerated ~ c(p0_tcolon_ulcerated, p1_tcolon_ulcerated, p2_tcolon_ulcerated)
    sescd_tcolon_affected ~ c(p0_tcolon_affected, p1_tcolon_affected, p2_tcolon_affected)
    sescd_tcolon_narrowing ~ c(p0_tcolon_narrowing, p1_tcolon_narrowing, p2_tcolon_narrowing)
    sescd_lcolon_ulcersize ~ c(p0_lcolon_ulcersize, p1_lcolon_ulcersize, p2_lcolon_ulcersize)
    sescd_lcolon_ulcerated ~ c(p0_lcolon_ulcerated, p1_lcolon_ulcerated, p2_lcolon_ulcerated)
    sescd_lcolon_affected ~ c(p0_lcolon_affected, p1_lcolon_affected, p2_lcolon_affected)
    sescd_lcolon_narrowing ~ c(p0_lcolon_narrowing, p1_lcolon_narrowing, p2_lcolon_narrowing)
    sescd_rectum_ulcersize ~ c(p0_rectum_ulcersize, p1_rectum_ulcersize, p2_rectum_ulcersize)
    sescd_rectum_ulcerated ~ c(p0_rectum_ulcerated, p1_rectum_ulcerated, p2_rectum_ulcerated)
    sescd_rectum_affected ~ c(p0_rectum_affected, p1_rectum_affected, p2_rectum_affected)
    sescd_rectum_narrowing ~ c(p0_rectum_narrowing, p1_rectum_narrowing, p2_rectum_narrowing)
    # Binary dropout endpoint (observed only at the end-of-induction
    # record): category 1 = continues into maintenance (DV 0), category
    # 2 = drops out (DV 1).
    dropout ~ c(prob_continue)
  })
}
