Robarge_2017_efavirenz <- function() {
  description <- paste(
    "Two-compartment population PK model for a single 600 mg oral dose of",
    "efavirenz in 73 HIV-seronegative adult volunteers (Robarge 2017), with",
    "parallel zero- and first-order absorption (independent lag times t_lag1",
    "and t_lag2, zero-order duration D2), allometric fat-free-mass scaling on",
    "CL/F (exponent 3/4, fixed), allometric fat-mass scaling on Vp/F",
    "(exponent 1, fixed), and CYP2B6 metaboliser status (normal /",
    "intermediate / slow) reducing CL/F by 0%, 25% and 51% respectively.",
    "Bioavailability is fixed to F = 1 (no IV reference formulation);",
    "the first-order absorption fraction F1 = 0.414 was estimated and the",
    "zero-order fraction F2 = 1 - F1 = 0.586 was assigned by mass balance.",
    "All absorption-related typical values (F1, t_lag1, K_a, t_lag2, D2)",
    "were estimated in an interim model and then fixed prior to covariate",
    "evaluation; the IIVs on those absorption parameters were re-estimated",
    "in the final model. Block-structured between-subject variability is",
    "estimated on (CL/F, Q/F, V_p/F) with correlations rho(CL/F, V_p/F) =",
    "0.196 and rho(Q/F, V_p/F) = 0.849; rho(CL/F, Q/F) was fixed at 0."
  )
  reference <- paste(
    "Robarge JD, Metzger IF, Lu J, Thong N, Skaar TC, Desta Z, Bies RR (2017).",
    "Population pharmacokinetic modeling to estimate the contributions of",
    "genetic and nongenetic factors to efavirenz disposition.",
    "Antimicrob Agents Chemother 61(1):e01813-16.",
    "doi:10.1128/AAC.01813-16.",
    sep = " "
  )
  vignette <- "Robarge_2017_efavirenz"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass derived from total body weight, height and sex via the Janmahasatian et al. (2005) semimechanistic formula; FFM = TBW - FM.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline). Drives allometric scaling on CL/F with exponent 3/4 (fixed) and reference value 56 kg = the cohort median calculated FFM per Robarge 2017 Table 1. Janmahasatian formula citation: Clin Pharmacokinet 2005;44(10):1051-1065.",
      source_name        = "FFM"
    ),
    FM = list(
      description        = "Fat mass derived as the difference between total body weight and calculated fat-free mass (FM = TBW - FFM), where FFM is computed by the Janmahasatian et al. (2005) formula.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline). Drives allometric scaling on V_p/F with exponent 1 (fixed) and reference value 19 kg = the cohort median calculated FM per Robarge 2017 Table 1. Robarge 2017 Materials and Methods 'Covariate model development' paragraph 2: 'FM = TBW - FFM'.",
      source_name        = "FM"
    ),
    CYP2B6_IM = list(
      description        = "1 = CYP2B6 intermediate-metabolizer phenotype (one CYP2B6 reduced-function star allele), 0 = otherwise. Reference category (both CYP2B6_IM and CYP2B6_SM equal to 0) is the CYP2B6 normal-metabolizer (extensive) phenotype.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Normal / extensive metabolizer (CYP2B6_IM = 0 and CYP2B6_SM = 0)",
      notes              = "Time-fixed (germline genotype). Star alleles assigned per the Human Cytochrome P450 Allele Nomenclature Database; functional consequence per Robarge 2017 Materials and Methods 'CYP nomenclature and predicted metabolizer status' and Table S1 in the supplemental material. The intermediate phenotype includes for example a CYP2B6 *1/*18 genotype (Robarge 2017 Discussion paragraph 2). The Robarge 2017 cohort frequency for intermediate metabolizers is summarised in Table S1; the canonical reference category is the unmodified CL/F (reduction factor 1.0).",
      source_name        = "CYP2B6 metabolizer status (intermediate level)"
    ),
    CYP2B6_SM = list(
      description        = "1 = CYP2B6 slow-metabolizer phenotype (Robarge's poor-metabolizer level, e.g. CYP2B6 *6/*6), 0 = otherwise. Reference category is the CYP2B6 normal-metabolizer (extensive) phenotype.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Normal / extensive metabolizer (CYP2B6_IM = 0 and CYP2B6_SM = 0)",
      notes              = "Time-fixed (germline genotype). Robarge 2017 uses 'slow' for what other authors call 'poor' (CYP2B6 *6/*6 homozygotes etc.). Star allele assignment per Robarge 2017 Materials and Methods 'CYP nomenclature and predicted metabolizer status' and Table S1.",
      source_name        = "CYP2B6 metabolizer status (slow / poor level)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 73L,
    n_studies      = 1L,
    age_range      = "18-50 years (median 24)",
    age_median     = "24 years",
    weight_range   = "53.0-103.6 kg (median 72.7)",
    weight_median  = "72.7 kg",
    height_range   = "1.55-1.98 m (median 1.76)",
    bmi_range      = "17.8-32.2 kg/m^2 (median 24.0); 9 of 73 subjects (12%) clinically obese (BMI > 30)",
    ffm_range      = "35.6-75.1 kg (median 56.4)",
    fm_range       = "6.3-42.7 kg (median 18.8)",
    sex_female_pct = 100 * 27 / 73,
    race_ethnicity = c(
      Caucasian        = 100 * 52 / 73,
      AfricanAmerican  = 100 * 16 / 73,
      Asian            = 100 * 3  / 73,
      Indian           = 100 * 1  / 73,
      AmericanIndian   = 100 * 1  / 73
    ),
    disease_state  = "HIV-seronegative healthy volunteers free of significant medical conditions; nonsmokers or willing to refrain from tobacco / marijuana for >= 1 month prior to and during the study; female subjects confirmed non-pregnant.",
    dose_range     = "Single oral 600 mg efavirenz dose (Sustiva tablet; Bristol-Myers Squibb) administered in the morning after overnight fasting, followed 1 hour later by an unrelated CYP drug cocktail (250 mg tolbutamide + 20 mg omeprazole + 150 mg caffeine + 1 mg midazolam) used for separate CYP-activity phenotyping.",
    regions        = "United States (Indiana University School of Medicine Clinical Research Center, Indianapolis, IN)",
    sampling       = "1,132 plasma efavirenz concentrations from 73 subjects: nominal 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 10, 12, 16, 24, 48, 72, and 144 h post-dose; 14-16 samples per subject in 69 of 73 subjects (3 early-withdrawal subjects with 8, 8 and 14 samples; one subject with a single 145 h sample due to misplaced samples). 2 of 1,134 concentrations were below the LLOQ of 1 ng/mL and were excluded.",
    notes          = "ClinicalTrials.gov NCT00668395; subjects enrolled between August 2007 and April 2010. Body composition (FFM, FM) calculated by the Janmahasatian et al. (2005) semimechanistic formula from total body weight, height and sex. Estimation method FOCE-I in NONMEM 7.3."
  )

  ini({
    # ============================================================================
    # Absorption (parallel zero- and first-order routes).
    # Robarge 2017 Table 2: all five absorption typical values (F1, t_lag1, K_a,
    # t_lag2, D2) were estimated in an interim model and then fixed prior to
    # covariate evaluation; the bootstrap median / %RSE columns of Table 2 are
    # blank for these five rows, confirming their fixed status in the final fit.
    # Wrap each in fixed() per parameter-names.md "Fixed parameters" rules.
    # ============================================================================

    logitfdepot <- fixed(log(0.414 / (1 - 0.414)))
    label("Logit of first-order absorption fraction F1 (unitless; F1 = 0.414 fixed, F2 = 1 - F1 = 0.586 by mass balance)")  # Robarge 2017 Table 2: F1 = 0.414 (fixed); F2 = 0.586 (fixed); Results 'Population PK model' paragraph 2

    lka     <- fixed(log(0.504))
    label("First-order absorption rate constant k_a (1/h, fixed)")  # Robarge 2017 Table 2: K_a = 0.504 1/h (fixed)

    ltlag1  <- fixed(log(1.97))
    label("Absorption lag time prior to first-order absorption t_lag1 (h, fixed)")  # Robarge 2017 Table 2: t_lag1 = 1.97 h (fixed)

    ltlag2  <- fixed(log(0.445))
    label("Absorption lag time prior to zero-order absorption t_lag2 (h, fixed)")  # Robarge 2017 Table 2: t_lag2 = 0.445 h (fixed)

    ldur2   <- fixed(log(0.675))
    label("Duration of zero-order absorption into central D2 (h, fixed)")  # Robarge 2017 Table 2: D2 (duration of zero-order absorption) = 0.675 h (fixed)

    # ============================================================================
    # Disposition (estimated; reference covariate values FFM = 56 kg, FM = 19 kg,
    # CYP2B6 normal metabolizer).
    # Robarge 2017 Table 2 "Final model estimate" column.
    # ============================================================================

    lcl <- log(7.52)
    label("CL/F (L/h) -- typical value at FFM = 56 kg, CYP2B6 normal metabolizer")  # Robarge 2017 Table 2: CL/F = 7.52 L/h (%RSE 1.70; bootstrap median 7.56)

    lvc <- log(125)
    label("Vc/F (L) -- apparent central volume of distribution")  # Robarge 2017 Table 2: V_c/F = 125 L (%RSE 6.03; bootstrap median 127)

    lvp <- log(374)
    label("Vp/F (L) -- typical value at FM = 19 kg")  # Robarge 2017 Table 2: V_p/F = 374 L (%RSE 3.53; bootstrap median 386)

    lq  <- log(32.3)
    label("Q/F (L/h) -- apparent inter-compartmental clearance")  # Robarge 2017 Table 2: Q/F = 32.3 L/h (%RSE 5.21; bootstrap median 33.3)

    # ============================================================================
    # Covariate effects.
    # ============================================================================

    # --- Allometric exponents on body composition (fixed). ---
    # Robarge 2017 Table 2 / Results 'Covariate model development' paragraphs 1-2:
    # "An allometric exponent of 3/4 was fixed in the final model as estimation
    # of the exponent did not significantly lower the overall objective function
    # ... bootstrap estimation did not support an alternative parameterization
    # (mean of 100 bootstrap simulations = 0.72, 95% CI 0.60-0.83)."
    # "An allometric exponent of 1 was fixed in the final model as estimation of
    # the exponent resulted in a small improvement in overall model fit (Delta
    # OFV = -4.686, P < 0.05) and was not strongly supported by bootstrap
    # estimation (mean of 100 bootstrap simulations = 1.13, 95% CI 0.79-1.46)."
    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent of (FFM/56) on CL/F (unitless; fixed at 3/4)")  # Robarge 2017 Table 2 / Results 'Covariate model development' para 1
    e_fm_vp  <- fixed(1.0)
    label("Allometric exponent of (FM/19) on V_p/F (unitless; fixed at 1)")  # Robarge 2017 Table 2 / Results 'Covariate model development' para 2

    # --- CYP2B6 metaboliser-status effects on CL/F (estimated). ---
    # Robarge 2017 Table 2: reduction factor 0.752 (intermediate; %RSE 5.54,
    # bootstrap median 0.740) and 0.490 (slow / poor; %RSE 7.22, bootstrap
    # median 0.494). Encoded as multiplicative log-ratios so the etalcl IIV
    # applies uniformly on the log-CL scale across all three metaboliser groups.
    e_cyp2b6_im_cl <- log(0.752)
    label("Log-ratio of CL/F for CYP2B6 intermediate vs normal metabolizer (unitless)")  # Robarge 2017 Table 2: factor 0.752 (-25% vs normal)
    e_cyp2b6_sm_cl <- log(0.490)
    label("Log-ratio of CL/F for CYP2B6 slow / poor vs normal metabolizer (unitless)")  # Robarge 2017 Table 2: factor 0.490 (-51% vs normal)

    # ============================================================================
    # Between-subject variability (BSV).
    # Robarge 2017 Table 2 reports BSV as standard deviations on the log scale,
    # so omega^2 = SD^2 directly (no CV->variance conversion needed). Correlations
    # rho(CL/F, V_p/F) = 0.196 and rho(Q/F, V_p/F) = 0.849 are also reported;
    # rho(CL/F, Q/F) was not estimated and is fixed to 0.
    # Covariances: cov = rho * sd_i * sd_j.
    #   cov(lcl, lq)  = 0
    #   cov(lcl, lvp) = 0.196 * 0.257 * 0.374 = 0.018841
    #   cov(lq,  lvp) = 0.849 * 0.671 * 0.374 = 0.213099
    # ============================================================================

    # Block omega lower triangle, row-major:
    #   row 1: var(etalcl)                                  = 0.257^2  = 0.066049  (Table 2: SD 0.257; %RSE 6.85)
    #   row 2: cov(etalcl, etalq) FIXED 0; var(etalq)       = 0.671^2  = 0.450241  (Table 2: SD 0.671; %RSE 15.77)
    #   row 3: cov(etalcl, etalvp) = 0.196 * 0.257 * 0.374  = 0.018841 (Table 2 footnote c: rho = 0.196)
    #          cov(etalq,  etalvp) = 0.849 * 0.671 * 0.374  = 0.213099 (Table 2 footnote c: rho = 0.849)
    #          var(etalvp)                                   = 0.374^2  = 0.139876 (Table 2: SD 0.374; %RSE 7.01)
    etalcl + etalq + etalvp ~ c(
      0.066049,
      fixed(0), 0.450241,
      0.018841, 0.213099, 0.139876
    )

    etalvc    ~ 0.101124  # var = 0.318^2 (Robarge 2017 Table 2: SD 0.318; %RSE 7.16)
    etalka    ~ 0.931225  # var = 0.965^2 (Robarge 2017 Table 2: SD 0.965; %RSE 7.72)
    etaltlag1 ~ 0.073441  # var = 0.271^2 (Robarge 2017 Table 2: SD 0.271; %RSE 7.85)
    etaltlag2 ~ 0.223729  # var = 0.473^2 (Robarge 2017 Table 2: SD 0.473; %RSE 7.53)
    etaldur2  ~ 0.494209  # var = 0.703^2 (Robarge 2017 Table 2: SD 0.703; %RSE 6.83)

    # ============================================================================
    # Residual error.
    # Robarge 2017 Table 2 reports the proportional and additive residual error
    # in the "Final model estimate" column as NONMEM $SIGMA variances:
    #   proportional sigma^2 = 0.016     => proportional SD = sqrt(0.016) = 0.1265
    #   additive     sigma^2 = 4270      => additive    SD = sqrt(4270)   = 65.35 nmol/L
    # The bootstrap median for the proportional row (0.130) confirms the
    # variance reading of the point estimate column (sqrt(0.016) = 0.1265 is
    # consistent with the bootstrap SD 0.130). The table footnote's
    # "expressed as standard deviations" applies cleanly to the BSV rows
    # (where point estimate and bootstrap median agree to 3 sig figs across
    # every parameter); for the two residual-error rows the point-estimate
    # column is the underlying NONMEM $SIGMA variance, which is the standard
    # NONMEM output convention.
    # See vignette section "Assumptions and deviations" for the full
    # provenance reading.
    #
    # The model below uses concentrations in mg/L (= ug/mL); converting the
    # variance reading of the additive residual:
    #   65.35 nmol/L * (315.6745 g/mol) / 1e9 (mg/(nmol*L)) per (g/L per (mg/L))
    #     = 65.35 * 3.15675e-4 mg/L
    #     = 0.020628 mg/L
    # (efavirenz molecular weight 315.6745 g/mol; PubChem CID 64139.)
    # ============================================================================

    propSd <- 0.1265
    label("Proportional residual SD (unitless; ~12.65% CV)")  # Robarge 2017 Table 2: sigma^2_prop = 0.016 (variance, %RSE 8.90); sqrt(0.016) = 0.1265; bootstrap median SD = 0.130
    addSd  <- 0.020628
    label("Additive residual SD (mg/L)")  # Robarge 2017 Table 2: sigma^2_add = 4270 (nmol/L)^2 (variance, %RSE 10.78); sqrt(4270) = 65.35 nmol/L = 65.35 * 315.6745 / 1e6 mg/L = 0.020628 mg/L
  })

  model({
    # ---- 1. Derived covariate / typical-value terms ---------------------------
    # CYP2B6 metabolizer factor: reference (normal) has both indicators = 0.
    # The exponent (e_cyp2b6_im_cl * CYP2B6_IM + e_cyp2b6_sm_cl * CYP2B6_SM) is
    # 0 for normal subjects, log(0.752) for intermediate, log(0.490) for slow.
    ltvcl <- lcl +
      e_cyp2b6_im_cl * CYP2B6_IM +
      e_cyp2b6_sm_cl * CYP2B6_SM

    # ---- 2. Individual PK parameters ------------------------------------------
    # Robarge 2017 final-model equations (Table 2 footnote a):
    #   CL/F = CL_TV * (FFM/56)^(3/4) * factor_CYP2B6  for normal metabolizer factor = 1
    #   V_p/F = V_p,TV * (FM/19)^1
    f1     <- 1 / (1 + exp(-logitfdepot))           # F1 ~= 0.414 (fixed)
    ka     <- exp(lka     + etalka)
    tlag1  <- exp(ltlag1  + etaltlag1)
    tlag2  <- exp(ltlag2  + etaltlag2)
    dur2   <- exp(ldur2   + etaldur2)

    cl <- exp(ltvcl + etalcl) * (FFM / 56)^e_ffm_cl
    vc <- exp(lvc   + etalvc)
    vp <- exp(lvp   + etalvp) * (FM  / 19)^e_fm_vp
    q  <- exp(lq    + etalq)

    # ---- 3. Micro-constants ---------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- 4. ODE system: two-compartment with parallel zero / first-order
    #        absorption. The user-supplied dosing dataset is expected to provide
    #        TWO dose records per administration -- one targeting depot (for the
    #        first-order route, cmt = "depot", rate = 0) and one targeting
    #        central (for the modeled zero-order route, cmt = "central",
    #        rate = -2). Both records carry the same amt; the f() statements
    #        below split the dose into the F1 / 1-F1 fractions. ----------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-               k12 * central - k21 * peripheral1

    # ---- 5. Bioavailability and lag-time / duration assignments ---------------
    f(depot)     <- f1
    alag(depot)  <- tlag1
    f(central)   <- 1 - f1
    alag(central) <- tlag2
    dur(central) <- dur2

    # ---- 6. Observation and error model ---------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
