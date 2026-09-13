Moein_2025_etrolizumab <- function() {
  description <- paste0(
    "Two-compartment population PK model for etrolizumab, an IgG1 ",
    "humanized anti-beta7 integrin monoclonal antibody, with first-order ",
    "SC absorption and a clearance that decreases exponentially with time ",
    "since the SECOND dose, in adults with moderately-to-severely active ",
    "Crohn's disease or ulcerative colitis (Moein 2025, n = 2312 with ",
    "9200 PK observations pooled over eight trials, of whom 864 CD ",
    "patients contributing 4317 observations came from the phase 3 ",
    "BERGAMOT study, NCT02394028). This model UPDATES the ",
    "ulcerative-colitis-only predecessor modellib('Moein_2022_etrolizumab') ",
    "in two ways: the time-dependent clearance is now a continuous ",
    "function of time since the second dose rather than the predecessor's ",
    "stepwise function of the most-recent dose time, and Crohn's disease ",
    "data are added so that bioavailability is estimated per indication. ",
    "Typical maximum clearance reduction is 22.0% with an onset half-life ",
    "of 3.45 weeks; baseline body weight, albumin and C-reactive protein ",
    "are the most influential covariates on exposure. Six companion ",
    "landmark logistic exposure-response models in the ",
    "Moein_2025_etrolizumab_* family consume this model's predicted ",
    "single-dose week-4 trough as their exposure metric."
  )
  reference <- paste(
    "Moein A, Ribbing J, Ibrahim MMA, Zhang W, Kassir N.",
    "Population pharmacokinetics and exposure-response relationships of",
    "etrolizumab in patients with moderately-to-severely active Crohn's",
    "disease.",
    "J Clin Pharmacol. 2025;65(10):1208-1219.",
    "doi:10.1002/jcph.70043.",
    "Parameter estimates from Table 1; the time-dependent clearance is",
    "Equation 1 and Table 1 footnote c; covariate forms are Table 1",
    "footnotes a, b, d and g; the covariate-effect magnitudes were",
    "independently confirmed against the fourteen printed ratios of the",
    "Figure 1 forest plot. Predecessor UC-only model: Moein A, Lu T,",
    "Jonsson S, et al. CPT Pharmacometrics Syst Pharmacol.",
    "2022;11(9):1244-1255. doi:10.1002/psp4.12846; see",
    "modellib('Moein_2022_etrolizumab').",
    sep = " "
  )
  vignette <- "Moein_2025_etrolizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # `tssd` is not a biological compartment: it is a zero-initialised
  # elapsed-time integrator whose value reproduces Equation 1's TSSD (time
  # in days since the second dose). It is declared here so
  # checkModelConventions() does not flag it as an unregistered
  # compartment role. See the model() block for why an ODE state is the
  # right vehicle.
  paper_specific_compartments <- c("tssd")

  compartmentData <- list(
    depot       = list(analyte = "etrolizumab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "etrolizumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "etrolizumab", units = "mg", specimen = "plasma", verified = TRUE),
    tssd        = list(analyte = "not applicable", units = "day", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Allometric power scaling (WT / 70)^exponent, with DIFFERENT ",
        "exponents on the elimination/disposition-clearance pair (CL, Q; ",
        "exponent 0.819) and on the volume pair (Vc, Vp; exponent 0.752). ",
        "Mechanistic covariate carried over from the predecessor model, so ",
        "it is part of the base model rather than a selected covariate. ",
        "NOTE the exponent assignment: Moein 2025 Table 1 rows and their ",
        "footnotes a and b both put 0.819 on CL/Q and 0.752 on Vc/Vp, but ",
        "the Results prose transposes them ('exponents of 0.75 and 0.82, ",
        "respectively' for clearance and volume). The table assignment is ",
        "the correct one and is what is encoded here; see the vignette ",
        "Errata for the forest-plot arithmetic that settles it."
      ),
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Exponential effect on CL: exp(theta * (ALB - 41)). The reference ",
        "41 g/L is the reference-patient value in the Table 1 footnote ",
        "block and equals the pooled UC+CD baseline median (Table S4). ",
        "Clearance falls as albumin rises (theta < 0)."
      ),
      source_name        = "ALB"
    ),
    CRP = list(
      description        = "Baseline C-reactive protein",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Exponential effect on the NATURAL LOG of CRP: ",
        "exp(theta * (log(CRP) - log(5.47))). The log transform is ",
        "load-bearing -- Table 1 names the row 'Log(CRP) on CL' and ",
        "footnote f reads 'per log unit CRP', and the Results quote a ",
        "7.8% change per natural-log unit, which is exp(0.0748). The ",
        "reference 5.47 mg/L is the Table 1 reference-patient value ",
        "(the Figure 1 caption rounds it to 5.48 mg/L); it equals the ",
        "pooled baseline median of 5.46 mg/L in Table S4. Standard (not ",
        "high-sensitivity) CRP assay, as is typical of ",
        "moderate-to-severe IBD cohorts. CRP must be strictly positive."
      ),
      source_name        = "CRP"
    ),
    CRCL = list(
      description        = "Baseline BSA-normalized estimated glomerular filtration rate",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Exponential effect on CL: exp(theta * (CRCL - 94.9)). Reported ",
        "by Moein 2025 as 'GFR', a creatinine-based BSA-normalized ",
        "estimate, which is the creatinine-estimate branch of the ",
        "canonical CRCL definition. The reference 94.9 mL/min/1.73 m^2 ",
        "is the Table 1 reference-patient value and equals the pooled ",
        "baseline median (Table S4). The effect is small and positive ",
        "(0.2% higher CL per unit) and reaches only a 0.87-1.10 fold ",
        "range over the 2.5th-97.5th covariate percentiles (Figure 1)."
      ),
      source_name        = "GFR"
    ),
    SCORE_SESCD = list(
      description        = "Baseline Simple Endoscopic Score for Crohn's Disease (SES-CD)",
      units              = "(score, 0-56)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Exponential effect on CL: exp(theta * (SCORE_SESCD - 12)). The ",
        "reference score 12 is the Table 1 reference-patient value and ",
        "equals the CD baseline median (Table S4). Recorded only in CD ",
        "patients -- Moein 2025 scores UC patients at the reference ",
        "value so the effect cancels, and the Figure 1 percentiles for ",
        "this row are computed over CD patients alone. Set ",
        "SCORE_SESCD = 12 for UC subjects to reproduce the published ",
        "model. Higher endoscopic severity raises clearance."
      ),
      source_name        = "SES-CD"
    ),
    ADA_TITER = list(
      description        = "Time-varying antidrug antibody titer, cumulative maximum carried forward",
      units              = "(titer units)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Exponential effect on CL: exp(theta * ADA_TITER), i.e. centered ",
        "at zero titer, so ADA-negative samples contribute no effect. ",
        "Table 1 footnote d is explicit that the covariate is the ",
        "subject's CUMULATIVE MAXIMUM titer with the maximum carried ",
        "forward, so the column is non-decreasing in time within a ",
        "subject -- it is not the instantaneous titer. ADA-negative is ",
        "encoded as 0 on the linear-titer convention (the same ",
        "convention as the predecessor Moein_2022_etrolizumab, and ",
        "distinct from the reciprocal-dilution convention where ",
        "negatives are coded 1). The effect is minimal: 3.5% higher CL ",
        "per titer unit."
      ),
      source_name        = "ADAT"
    ),
    PRIOR_TNF = list(
      description        = "Prior anti-TNF biologic therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no prior anti-TNF, i.e. TNF-naive)",
      notes              = paste0(
        "FRACTIONAL multiplicative effect on CL: ",
        "(1 + theta * PRIOR_TNF), giving 5.86% higher clearance in ",
        "TNF-experienced patients. The fractional (1 + theta) form -- ",
        "rather than exp(theta) -- is the one carried by Table 1 ",
        "footnote g, and it is what reproduces the printed Figure 1 ",
        "ratio of 0.925 (exp(theta) would give 0.923). The predecessor ",
        "Moein_2022_etrolizumab uses the same fractional form for its ",
        "categorical CL covariates."
      ),
      source_name        = "Prior anti-TNF"
    ),
    IBD_CD = list(
      description        = "Crohn's disease indicator within the pooled UC + CD analysis",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ulcerative colitis with left-sided colitis, which also implies DISEXT_EP = 0 and DISEXT_OTHER = 0)",
      notes              = paste0(
        "ADDITIVE shift on the LOGIT of bioavailability, not on CL and ",
        "not multiplicative: logit(F) = logit(0.743) + theta * IBD_CD. ",
        "Table 1 footnote g states that the covariates on F are on the ",
        "logit scale, and the arithmetic confirms it exactly -- ",
        "expit(logit(0.743) - 0.314) = 0.678, the paper's printed CD ",
        "point estimate. Mutually exclusive with the UC disease-extent ",
        "indicators: a CD subject has DISEXT_EP = DISEXT_OTHER = 0."
      ),
      source_name        = "indication"
    ),
    DISEXT_EP = list(
      description        = "Ulcerative colitis disease extent: extensive colitis / pancolitis (vs left-sided colitis)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (left-sided colitis; implies DISEXT_OTHER = 0 as well)",
      notes              = paste0(
        "ADDITIVE shift on the LOGIT of bioavailability. Moein 2025 ",
        "pools extensive/pancolitis with the small 'other' UC extent ",
        "group into a single 'UC not left-sided colitis' category ",
        "carrying ONE coefficient, so DISEXT_EP and DISEXT_OTHER share ",
        "the estimate e_ucother_fdepot; expit(logit(0.743) - 0.244) = ",
        "0.693, the paper's printed value. This differs from the ",
        "predecessor Moein_2022_etrolizumab, which estimated the two ",
        "extents separately (and on CL); here they are deliberately ",
        "collapsed. Zero for CD subjects."
      ),
      source_name        = "DISSPR"
    ),
    DISEXT_OTHER = list(
      description        = "Ulcerative colitis disease extent: other (neither left-sided colitis nor extensive/pancolitis)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (left-sided colitis)",
      notes              = paste0(
        "ADDITIVE shift on the LOGIT of bioavailability, sharing the ",
        "single 'UC not left-sided colitis' coefficient with DISEXT_EP ",
        "(see that entry). Mutually exclusive with DISEXT_EP. Only 20 of ",
        "2312 subjects (1%) fall in this group (Table S5). Zero for CD ",
        "subjects."
      ),
      source_name        = "DISSPR"
    )
  )

  # Documented but NOT referenced in model(). Moein 2025 retained one
  # covariate that rxode2's error-model DSL cannot express: a covariate
  # on the magnitude of the residual error itself. `add()` and `prop()`
  # accept only bare estimated parameters, not expressions, so
  # `add(addSd * ruv_scale)` is a parse error. The model below therefore
  # encodes the PHASE 3 residual error, which is both the paper's
  # reference stratum and the only stratum any Crohn's-disease
  # simulation can occupy (every phase I/II subject is a UC patient).
  # The value is preserved here so a downstream user can reconstruct the
  # phase I/II stratum by hand. See the vignette Errata.
  covariatesDataExcluded <- list(
    STUDY_ETRO_PHASE12 = list(
      description        = "Phase I or phase II study indicator, a residual-error stratum",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (phase 3 study)",
      notes              = paste0(
        "RETAINED IN THE PUBLISHED MODEL but not expressible in ",
        "rxode2; see the block comment above. Table 1 reports ",
        "'Clinical study Phase I/II on RUV' = -0.230 (RSE 24.4%) as a ",
        "fractional change on the OVERALL residual error, scaling BOTH ",
        "the additive and the proportional component by ",
        "(1 - 0.230) = 0.770 -- a 23.0% smaller residual error in the ",
        "phase I/II studies than in phase 3 (footnotes e and g). ",
        "Scaling both components is equivalent to scaling the combined ",
        "residual SD, because sqrt((s*add)^2 + ((s*prop)*f)^2) = ",
        "s * sqrt(add^2 + (prop*f)^2); so the phase I/II stratum is ",
        "recovered by setting addSd = 0.426 * 0.770 = 0.328 ug/mL and ",
        "propSd = 0.201 * 0.770 = 0.155 via an ini() override. ",
        "Table S2 footnote 2 records that the analysis plan had ",
        "pre-specified this covariate on CL, F and ka but that it was ",
        "instead investigated on the residual-error magnitude, to ",
        "capture manufacturing-process and PK-assay changes between the ",
        "early studies and phase 3. Only 119 of 2312 subjects (5%) are ",
        "phase I/II and all of them are UC patients (Table S5)."
      ),
      source_name        = "Clinical study phase"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2312L,
    n_studies      = 8L,
    n_observations = "9200 etrolizumab serum concentrations (4317 of them from the 864 Crohn's disease patients)",
    age_range      = "18.0-79.0 years",
    age_median     = "37.0 years (UC 38.0, CD 36.0)",
    weight_range   = "35.0-216 kg",
    weight_median  = "72.0 kg (UC 72.0, CD 71.0)",
    sex_female_pct = 44,
    race_ethnicity = c(
      White = 83, Asian = 8, `Black or African American` = 2,
      Other = 3, `American Indian or Alaska Native` = 0,
      `Native Hawaiian or Pacific Islander` = 0, Multiple = 0, Unknown = 4
    ),
    disease_state  = "Moderately-to-severely active Crohn's disease (n = 864) or ulcerative colitis (n = 1448); TNF-naive and TNF-experienced",
    dose_range     = "SC: 105 mg Q4W, 210 mg Q4W with an extra 210 mg loading dose at week 2 (BERGAMOT), 315-420 mg (phase II), 0.5-3 mg/kg Q4W and 1-3 mg/kg single dose (phase I); IV: 0.3-10 mg/kg single dose and 4 mg/kg Q4W (phase I)",
    regions        = "Multinational: Eastern Europe 33%, Western Europe 25%, USA/Canada 23%, Asia 9%, Latin/Middle America 5%, Australia/New Zealand 4%, (South) Africa 1%",
    renal_function = "Baseline GFR 30.7-269 mL/min/1.73 m^2, median 94.9",
    notes          = paste0(
      "Baseline characteristics from Moein 2025 Tables S4 (continuous) ",
      "and S5 (categorical), 'All' column. Contributing studies ",
      "(Table S1): ABS4262g (phase I, UC, 38 subjects), EUCALYPTUS ",
      "(phase II, UC, 81), HIBISCUS I + HIBISCUS II (phase III, ",
      "TNF-naive UC, 285 combined), HICKORY (phase III, TNF inadequate ",
      "responder UC, 509), LAUREL (phase III, TNF-naive UC, 350), ",
      "GARDENIA (phase III, TNF-naive UC, 184) and BERGAMOT ",
      "(NCT02394028, phase III, CD, 864). n_studies counts HIBISCUS I ",
      "and HIBISCUS II separately; Table S1 lists them on one row, so ",
      "the table shows seven rows. GARDENIA was an external-validation ",
      "hold-out for the predecessor UC model but contributes to the ",
      "estimation data set here. Other pooled baseline medians: ",
      "albumin 41 g/L, CRP 5.46 mg/L, fecal calprotectin 1270 ug/g, ",
      "BMI 24.3 kg/m^2, disease duration 5.98 years; CD-only medians: ",
      "CDAI 314, SES-CD 12.0. 48% had prior anti-TNF therapy."
    )
  )

  ini({
    # ==================================================================
    # Structural PK parameters. Moein 2025 Table 1. All typical values
    # are for the reference patient defined in the Table 1 footnote
    # block: a phase 3, 70 kg ulcerative-colitis patient with
    # left-sided colitis, albumin 41 g/L, CRP 5.47 mg/L, GFR
    # 94.9 mL/min/1.73 m^2, SES-CD 12, no positive ADA titer and no
    # prior anti-TNF therapy.
    #
    # Every value below was cross-checked against the fourteen printed
    # ratios of the Figure 1 forest plot: reproducing the reference
    # patient's single-dose week-4 trough and each covariate
    # perturbation returns all fourteen ratios to within 0.3%. See the
    # vignette Source trace and Errata sections.
    # ==================================================================
    lka         <- log(0.213);   label("First-order SC absorption rate constant (1/day)")                       # Table 1: ka = 0.213 per day (RSE 5.55%)
    lcl         <- log(0.282);   label("Clearance before the second dose, for a 70 kg reference patient (L/day)") # Table 1: CL = 0.282 L/day (RSE 7.35%)
    lvc         <- log(2.78);    label("Central volume of distribution for a 70 kg reference patient (L)")        # Table 1: Vc = 2.78 L (RSE 5.00%)
    lvp         <- log(1.73);    label("Peripheral volume of distribution for a 70 kg reference patient (L)")     # Table 1: Vp = 1.73 L (RSE 9.85%)
    lq          <- log(0.489);   label("Intercompartmental clearance for a 70 kg reference patient (L/day)")       # Table 1: Q = 0.489 L/day (RSE 10.4%)
    logitfdepot <- logit(0.743); label("SC bioavailability of the reference UC left-sided-colitis patient (unitless logit)") # Table 1: F = 0.743 (RSE 7.48%); covariates on F act on this logit scale per footnote g

    # ------------------------------------------------------------------
    # Time-dependent clearance (Equation 1 / Table 1 footnote c):
    #   CL(TSSD) = CL_0 * (1 - Maxred * (1 - exp(-log(2)/(Onset*7) * TSSD)))
    # Maxred is a FRACTION in (0, 1) and is therefore held on the logit
    # scale, matching the predecessor model and letting its IIV (also
    # reported on the logit scale) attach directly.
    # ------------------------------------------------------------------
    logitmaxred <- logit(0.220); label("Maximum fractional reduction of CL over time (unitless logit)")          # Table 1: Maxred = 0.220 (RSE 3.52%); Results 95% CI 0.205-0.235
    lonset      <- log(3.45);    label("Half-life of the time-dependent CL change (weeks)")                    # Table 1: Onset = 3.45 weeks (RSE 9.03%); Results 95% CI 2.84-4.04 weeks

    # ------------------------------------------------------------------
    # Allometric exponents on body weight, reference 70 kg
    # (Table 1 footnotes a and b, TVP = Ppop * (WT/70)^theta).
    #
    # The exponents are NOT the conventional 0.75 / 1.0 pair and they
    # are NOT interchangeable: 0.819 goes on CL and Q, 0.752 on Vc and
    # Vp. The Results prose transposes them; the table does not. See the
    # WT covariateData note and the vignette Errata.
    # ------------------------------------------------------------------
    e_wt_cl_q  <- 0.819; label("Allometric exponent of WT on CL and Q (unitless)")                                               # Table 1 row 'Body weight on CL and Q' + footnote a: 0.819 (RSE 4.16%)
    e_wt_vc_vp <- 0.752; label("Allometric exponent of WT on Vc and Vp (unitless)")                                              # Table 1 row 'Body weight on Vc and Vp' + footnote b: 0.752 (RSE 5.90%)

    # ------------------------------------------------------------------
    # Continuous covariate effects on CL. Table 1 footnote d gives the
    # single exponential form for all of them:
    #   CovEff = exp(theta * (Cov - Cov_ref))
    # ------------------------------------------------------------------
    e_alb_cl   <- -0.0260;  label("Albumin effect on CL, exp(theta * (ALB - 41)) (per g/L)")                            # Table 1: Albumin on CL = -0.0260 (RSE 7.46%); Results '2.6% per g/L'
    e_crp_cl   <-  0.0748;  label("log-CRP effect on CL, exp(theta * (log(CRP) - log(5.47))) (per natural-log unit CRP)")   # Table 1: Log(CRP) on CL = 0.0748 (RSE 7.23%); Results '7.8% per log unit', = exp(0.0748)
    e_crcl_cl  <-  0.00202; label("GFR effect on CL, exp(theta * (CRCL - 94.9)) (per mL/min/1.73 m^2)")                 # Table 1: GFR on CL = 0.00202 (RSE 15.2%); Results '0.2% per mL/min/1.73 m2'
    e_sescd_cl <-  0.00656; label("SES-CD effect on CL, exp(theta * (SCORE_SESCD - 12)) (per score point)")             # Table 1: SES-CD on CL = 0.00656 (RSE 22.8%); Results '0.7% per unit score'
    e_adat_cl  <-  0.0342;  label("ADA titer effect on CL, exp(theta * ADA_TITER) (per titer unit)")                    # Table 1: ADAT on CL = 0.0342 (RSE 16.9%); Results '3.5% per titer unit'

    # ------------------------------------------------------------------
    # Categorical covariate effect on CL. Table 1 footnote g reports a
    # relative change versus the reference category, and the FRACTIONAL
    # (1 + theta) reading -- not exp(theta) -- is what reproduces the
    # printed Figure 1 ratio of 0.925 for a prior-anti-TNF patient.
    # ------------------------------------------------------------------
    e_priortnf_cl <- 0.0586; label("Prior anti-TNF fractional change in CL vs no prior anti-TNF (fraction)")                      # Table 1: Prior anti-TNF on CL = 0.0586 (RSE 24.2%); Results 'was 5.9% higher in TNF-experienced patients'

    # ------------------------------------------------------------------
    # Bioavailability covariates. Table 1 footnote g: "Since covariates
    # on F are on the logit scale". These are ADDITIVE SHIFTS on
    # logit(F), which the paper's own printed back-transforms confirm:
    #   expit(logit(0.743) - 0.244) = 0.693  (UC, not left-sided)
    #   expit(logit(0.743) - 0.314) = 0.678  (CD)
    # ------------------------------------------------------------------
    e_ucother_fdepot <- -0.244; label("UC-not-left-sided-colitis shift in logit(F) vs UC left-sided colitis (unitless logit)")          # Table 1: UC not left-sided colitis on F = -0.244 (RSE 33.8%); yields F = 0.693
    e_cd_fdepot      <- -0.314; label("Crohn's disease shift in logit(F) vs UC left-sided colitis (unitless logit)")                    # Table 1: CD indication on F = -0.314 (RSE 24.2%); yields F = 0.678

    # NOTE: Table 1's fourth covariate family -- the phase I/II
    # fractional change of -0.230 on the OVERALL residual error -- has
    # no ini() entry, because rxode2's error-model DSL cannot carry a
    # covariate on the residual magnitude. It is recorded in full in
    # `covariatesDataExcluded` above and in the vignette Errata.

    # ==================================================================
    # Between-subject variability. Table 1 reports each term on the
    # scale of the transform its parameter carries, and the note block
    # under Table 1 is explicit about which scale that is:
    #   - CL, Vtot, ka  -> approximate CV of a log-normal parameter, so
    #     the variance on the estimation (log) scale is log(1 + CV^2).
    #   - F, Maxred     -> SD on the LOGIT scale (footnote h), so the
    #     variance is the reported SD squared, with no back-transform.
    #
    # Footnote h's own consistency check confirms the logit reading: for
    # Maxred, p*(1-p)*SD / p = 0.220*0.780*0.784/0.220 = 0.611, exactly
    # the "approximate CV for Maxred" of 0.611 that the footnote quotes;
    # the same delta-method identity returns the footnote's F CVs of
    # 0.175, 0.209 and 0.219 for the three indication groups.
    # ==================================================================

    # CL and Maxred are correlated (Table 1: cor(CL,Maxred) = 0.263), so
    # they share a 2x2 block. The off-diagonal is the covariance implied
    # by that correlation and the two SDs:
    #   0.263 * sqrt(log(1 + 0.232^2)) * 0.784
    etalcl + etalogitmaxred ~ c(
      log(1 + 0.232^2),
      0.263 * sqrt(log(1 + 0.232^2)) * 0.784, 0.784^2
    )

    # A SINGLE eta on total volume. Table 1 footnote i: "IIV for Vtot,
    # that is, same magnitude and full correlation for the random effect
    # on Vc and Vp." Full correlation with equal magnitude IS one shared
    # random effect, so `etalvc` is deliberately reused on both volumes
    # in model() rather than declaring a second, perfectly-correlated
    # eta (which would be singular and would fail Cholesky).
    etalvc ~ log(1 + 0.132^2)                # Table 1: IIV Vtot CV = 0.132 (RSE 13.3%, shrinkage 71.2%)

    etalka ~ log(1 + 0.336^2)                # Table 1: IIV ka CV = 0.336 (RSE 14.4%, shrinkage 65.7%)

    etalogitfdepot ~ 0.681^2                 # Table 1: IIV F = 0.681, an SD on the logit scale per footnote h

    # ==================================================================
    # Residual error, combined additive + proportional. These are the
    # PHASE 3 (reference-stratum) magnitudes; the phase I/II stratum
    # scales both by (1 + e_phase12_ruv), applied in model().
    # ==================================================================
    propSd <- 0.201; label("Proportional residual error in the phase 3 stratum, reported as a CV (fraction)")                        # Table 1: Proportional residual error CV = 0.201 (RSE 2.79%)
    addSd  <- 0.426; label("Additive residual error in the phase 3 stratum (ug/mL)")                                   # Table 1: Additive residual error SD = 0.426 ug/mL (RSE 10.6%)
  })

  model({
    # ==================================================================
    # Time since the second dose (TSSD of Equation 1).
    #
    # Moein 2025 replaced the predecessor model's stepwise function of
    # the most-recent dose time with a CONTINUOUS function of time since
    # the second dose, and states that clearance is constant until that
    # second dose is given. TSSD therefore cannot be built from `tad()`
    # (which restarts at every dose) and it is not a simple lag off the
    # first dose either, because the second dose falls at week 2 in the
    # 210 mg loading arm and at week 4 on plain Q4W dosing.
    #
    # The integrator below is exact for ANY regimen and needs no data
    # column: `dosenum()` is 0 or 1 until the second dose is
    # administered, so the derivative is 0 and the state holds at its
    # zero initial condition; from the second dose onward the derivative
    # is 1 and the state accumulates elapsed time. For a single dose
    # `dosenum()` never reaches 2, so tssd stays 0 and clearance stays
    # constant -- which is exactly the behaviour the paper's single-dose
    # exposure metric (Ctrough,W4,adjusted) requires.
    #
    # Declared last among the ODEs so it does not disturb the
    # depot / central / peripheral1 compartment ordering.
    # ==================================================================
    maxred       <- expit(logitmaxred + etalogitmaxred)
    onset        <- exp(lonset)
    td_cl_factor <- 1 - maxred * (1 - exp(-log(2) / (onset * 7) * tssd))

    # Covariate effects on CL: exponential for the continuous set
    # (Table 1 footnote d), fractional for the categorical one
    # (footnote g).
    alb_cl      <- exp(e_alb_cl   * (ALB - 41))
    crp_cl      <- exp(e_crp_cl   * (log(CRP) - log(5.47)))
    crcl_cl     <- exp(e_crcl_cl  * (CRCL - 94.9))
    sescd_cl    <- exp(e_sescd_cl * (SCORE_SESCD - 12))
    adat_cl     <- exp(e_adat_cl  * ADA_TITER)
    priortnf_cl <- 1 + e_priortnf_cl * PRIOR_TNF

    # Individual PK parameters. Note the two distinct allometric
    # exponents and the single shared volume eta.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q *
      td_cl_factor * alb_cl * crp_cl * crcl_cl * sescd_cl *
      adat_cl * priortnf_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp <- exp(lvp + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(tssd)        <- (dosenum() >= 2)

    # Bioavailability. Additive shifts on the logit scale; the three
    # indication groups are mutually exclusive, so at most one shift
    # applies to any subject.
    f(depot) <- expit(
      logitfdepot + etalogitfdepot +
        e_ucother_fdepot * (DISEXT_EP + DISEXT_OTHER) +
        e_cd_fdepot * IBD_CD
    )

    # Concentration: dose in mg / volume in L = mg/L = ug/mL.
    Cc <- central / vc

    # Combined residual error, PHASE 3 stratum (see the
    # covariatesDataExcluded block for the phase I/II scaling).
    Cc ~ add(addSd) + prop(propSd)
  })
}
