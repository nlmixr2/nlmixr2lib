Sano_2023_fesoterodine <- function() {
  description <- "Population pharmacokinetic model for 5-hydroxymethyl tolterodine (5-HMT), the active metabolite of fesoterodine, in pediatric patients aged 6-17 years with overactive bladder or neurogenic detrusor overactivity. One-compartment disposition with first-order absorption, a first-order absorption lag and first-order elimination, fitted to 428 5-HMT concentrations from 142 patients pooled across the phase II study 1066 (NCT00857896) and the phase III study 1047 (NCT01557244). Body weight is allometrically scaled onto CL/F (exponent 0.75 FIXED) and Vd/F (exponent 1 FIXED) referenced to 35 kg; CYP2D6 poor metabolizers carry 0.546-fold CL/F, female patients 0.862-fold CL/F and 0.634-fold Vd/F, and the beads-in-capsule (BIC) formulation carries 0.648-fold bioavailability relative to the tablet reference. Absorption is flip-flop (ka 0.0897 1/h is far slower than kel = CL/Vd = 1.05 1/h), so the reported 7.73 h terminal half-life is absorption-rate-limited. Inter-individual variability on CL/F, Vd/F and ka is a full 3x3 covariance block and residual error is additive on the log-transformed concentration scale (i.e. log-normal). Doses are supplied in micrograms so that amount/volume yields ng/mL, matching the source NM-TRAN dataset."
  reference <- "Sano Y, Shoji S, Shahin M, Sweeney K, Darekar A, Malhotra BK. Population Pharmacokinetic and Pharmacodynamic Modeling of Fesoterodine in Pediatric Patients with Neurogenic Detrusor Overactivity. Eur J Drug Metab Pharmacokinet. 2023 May;48(3):257-269. doi:10.1007/s13318-023-00818-8. PMID: 36892805. PMCID: PMC10175358."
  vignette <- "Sano_2023_fesoterodine"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometrically scaled onto CL/F and Vd/F as (WT/35)^exponent, referenced to 35 kg (Sano 2023 Methods section 2.3 and Online Resource 8a). Both exponents were FIXED a priori (0.75 on CL/F, 1 on Vd/F) rather than estimated, per Sano 2023 Table 2. Observed range in the pharmacokinetic analysis population is 11.7-85.0 kg (median 33.6 kg).",
      source_name        = "BWT"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) -- the most common category in the analysis population (52.1% male) and the typical-value reference in Sano 2023 Eqs. 1-2",
      notes              = "Multiplicative effect on both CL/F (0.862-fold) and Vd/F (0.634-fold) for female relative to the male reference (Sano 2023 Table 2). The source NM-TRAN column SEX is coded 0 = male, 1 = female (Online Resource 4 variable list), which matches the SEXF canonical orientation directly with no value inversion. Note that the original clinical datasets used 1 = male, 2 = female and were recoded before modeling. Neither 95% CI excluded the null value of 1.0 (CL/F 0.716-1.02; Vd/F 0.161-1.26), so both sex effects are imprecise; they are retained here because the authors carried them into the final model.",
      source_name        = "SEX"
    ),
    CYP2D6_PM = list(
      description        = "CYP2D6 poor-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP2D6 extensive metabolizer) -- the reference scenario in Sano 2023 Table 2",
      notes              = "Multiplicative 0.546-fold effect on CL/F for poor metabolizers relative to the extensive-metabolizer reference (Sano 2023 Table 2; 95% CI 0.343-0.721 excludes the null value of 1.0). Sano 2023 Table 2 footnote a states that intermediate metabolizers (IM) and ultra-rapid metabolizers (UM) were pooled with extensive metabolizers, so this column is a two-level PM-vs-everyone-else indicator rather than a four-level phenotype. The source NM-TRAN column CYPbi is coded 1 = extensive metabolizer, 0 = poor metabolizer (Online Resource 4), i.e. the INVERSE of this canonical: derive CYP2D6_PM = 1 - CYPbi. Only 3 of 142 patients (2.1%) were poor metabolizers.",
      source_name        = "CYPbi (inverted: CYP2D6_PM = 1 - CYPbi)"
    ),
    FORM_CAPSULE = list(
      description        = "Fesoterodine beads-in-capsule (BIC) formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fesoterodine tablet) -- the tablet is the structural bioavailability anchor with F fixed at 1",
      notes              = "1 = beads-in-capsule (BIC), 0 = tablet. Multiplicative 0.648-fold effect on the extent of absorption (relative bioavailability) for BIC relative to the tablet reference (Sano 2023 Table 2; 95% CI 0.546-0.787 excludes the null value of 1.0). Same orientation as Gupta 2016 lenvatinib: F is fixed at 1 on the non-capsule (tablet) arm and the capsule arm carries the estimated relative bioavailability. In the source studies BIC was given only to study 1047 cohort 2 (patients weighing 25 kg or less, 2 and 4 mg QD) and tablets only to study 1047 cohort 1 and study 1066 (4 and 8 mg QD), so formulation is strongly confounded with body weight and dose in this dataset. The authors also tested a formulation effect on ka instead of on F; it worsened the objective function value and was not retained (Sano 2023 Results section 3.1). The model-predicted relative bioavailability of 64.8% is lower than the 79.9-87.9% observed by non-compartmental analysis in healthy adults (NCT02160158); Sano 2023 Discussion attributes the discrepancy to design differences and to unknown fed/fasted state in study 1047.",
      source_name        = "BIC"
    )
  )

  population <- list(
    species           = "human",
    n_subjects        = 142L,
    n_observations    = 428L,
    n_studies         = 2L,
    age_range         = "6-17 years (median 10, mean 10.1, SD 3.00)",
    weight_range      = "11.7-85.0 kg (median 33.6, mean 36.2, SD 16.1)",
    sex_female_pct    = 47.9,
    race_ethnicity    = c(White = 50.7, Black = 3.5, Asian = 44.4, Other = 1.4),
    cyp2d6_pm_pct     = 2.1,
    formulation_pct   = c(Tablet = 64.8, BIC = 35.2),
    disease_state     = "Pediatric patients with overactive bladder (OAB) or neurogenic detrusor overactivity (NDO). Study 1066 enrolled OAB patients aged 8-17 years with approximately half having NDO; study 1047 enrolled patients aged 6-17 years with symptoms of NDO.",
    dose_range        = "Fesoterodine 4 mg tablet QD then 8 mg tablet QD (study 1066, 4-week periods); fesoterodine 4 or 8 mg tablet QD (study 1047 cohort 1, patients over 25 kg); fesoterodine 2 or 4 mg beads-in-capsule QD (study 1047 cohort 2, patients 25 kg or less). Study 1047 also had an oxybutynin extended-release comparator arm that contributed no 5-HMT data.",
    regions           = "Multinational (NCT00857896 phase II study 1066; NCT01557244 phase III study 1047).",
    notes             = "Demographics from Sano 2023 Table 1 (pharmacokinetic analysis population). Sampling was sparse: study 1066 collected pre-dose and 0.5-2, 2-4, 4-6 h post-dose samples at week 4 plus 8-10, 10-14, 14-16 and 16-20 h post-dose samples at week 8; study 1047 collected up to three samples per patient at week 4. The sparse design and lack of information on the absorption phase produced moderate eta shrinkage on Vd/F (42.1%) and ka (45.0%), which the authors flag as a caveat for empirical-Bayes-based diagnostics and exposure metrics."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. All values are the final population
    # estimates in Sano 2023 Table 2 ("Value" column), for a typical
    # 35 kg MALE CYP2D6 EXTENSIVE METABOLIZER receiving the TABLET
    # formulation (Table 2 legend).
    #
    # Absorption is flip-flop in this model: ka (0.0897 1/h) is much
    # slower than kel = CL/F / (Vd/F) = 71.6 / 68.1 = 1.05 1/h, so the
    # terminal phase is absorption-rate-limited. This reproduces the
    # paper's reported typical half-life of 7.73 h (= log(2)/0.0897)
    # and Tmax,ss of 2.55 h, neither of which is recoverable from
    # log(2)/kel = 0.66 h.
    # ------------------------------------------------------------------
    lcl   <- log(71.6);   label("Apparent clearance CL/F (L/h)")                     # Sano 2023 Table 2 (RSE 6.7%; bootstrap median 71.4, 95% CI 63.0-81.1)
    lvc   <- log(68.1);   label("Apparent central volume Vd/F (L)")                  # Sano 2023 Table 2 (RSE 29.7%; bootstrap median 73.7, 95% CI 10.8-127)
    lka   <- log(0.0897); label("First-order absorption rate constant ka (1/h)")     # Sano 2023 Table 2 (RSE 5.99%; bootstrap median 0.0902, 95% CI 0.0804-0.105)
    ltlag <- log(0.285);  label("Absorption lag time ALAG (h)")                      # Sano 2023 Table 2 (RSE 45.1%; bootstrap median 0.279, 95% CI 0.00263-0.479)

    # ------------------------------------------------------------------
    # Covariate effects. Sano 2023 Eqs. 1-2 define categorical effects
    # multiplicatively (theta_i = theta_TV * theta_cov for the minor
    # category, theta_i = theta_TV for the reference category), and
    # continuous covariates as a power function centered on a reference
    # value. Online Resource 8a writes the assembled forms:
    #   (CL/F)_i = theta_CL/F * (BWT_i/35)^0.75 * theta_CYP * theta_SEX * exp(eta)
    #   (Vd/F)_i = theta_Vd/F * (BWT_i/35)^1    * theta_SEX * exp(eta)
    # Multiplicative factors are stored on the LOG scale here so that
    # model() can assemble every parameter in one mu-referenced
    # exponential; the linear-scale factor from Table 2 is given in each
    # trailing comment.
    # ------------------------------------------------------------------
    e_wt_cl        <- fixed(0.75);        label("Allometric exponent of body weight on CL/F (referenced to 35 kg)") # Sano 2023 Table 2 ("0.75, Fixed" -- fixed a priori, not estimated; $THETA(4) 0.75 FIX in Online Resource 4)
    e_wt_vc        <- fixed(1);           label("Allometric exponent of body weight on Vd/F (referenced to 35 kg)") # Sano 2023 Table 2 ("1, Fixed" -- fixed a priori, not estimated; $THETA(5) 1 FIX in Online Resource 4)
    e_cyp2d6_pm_cl <- log(0.546);         label("Log multiplicative effect of CYP2D6 poor metabolizer on CL/F")     # Sano 2023 Table 2: 0.546 (RSE 12.5%; bootstrap median 0.534, 95% CI 0.343-0.721)
    e_sexf_cl      <- log(0.862);         label("Log multiplicative effect of female sex on CL/F")                  # Sano 2023 Table 2: 0.862 (RSE 8.52%; bootstrap median 0.862, 95% CI 0.716-1.02)
    e_sexf_vc      <- log(0.634);         label("Log multiplicative effect of female sex on Vd/F")                  # Sano 2023 Table 2: 0.634 (RSE 30.6%; bootstrap median 0.635, 95% CI 0.161-1.26)
    lfdepot        <- log(0.648);         label("Log relative bioavailability of the BIC formulation vs the tablet reference (log fraction, dimensionless)") # Sano 2023 Table 2 "Effect of fesoterodine formulations on F": 0.648 (RSE 9.38%; bootstrap median 0.65, 95% CI 0.546-0.787)

    # ------------------------------------------------------------------
    # Inter-individual variability -- full OMEGA BLOCK(3) on CL/F, Vd/F
    # and ka (Online Resource 4 $OMEGA BLOCK(3)).
    #
    # Sano 2023 Table 2 reports the DIAGONAL entries as "% CV" and the
    # OFF-DIAGONAL entries as bare covariances. The %CV column is the
    # sqrt-of-variance convention (omega^2 = (%CV/100)^2), NOT the
    # log-normal sqrt(exp(omega^2)-1) convention. Two independent checks
    # confirm this reading:
    #   (1) Under the log-normal convention the Vd/F-ka covariance 0.435
    #       would imply a correlation of 1.15, which is impossible; under
    #       the sqrt-of-variance convention it is 0.88, which is valid.
    #   (2) The resulting variances (0.214, 1.300, 0.187) sit right on
    #       top of the Online Resource 4 $OMEGA initial estimates
    #       (0.2, 0.9, 0.5), as expected for a converged fit.
    # The assembled 3x3 matrix is positive definite (leading minors
    # 0.214, 0.190, 0.00735).
    #
    # Block order is the NONMEM lower-triangular row order, which is the
    # same order Table 2 lists and the same order rxode2 expects:
    #   var(CL), cov(CL,Vd), var(Vd), cov(CL,ka), cov(Vd,ka), var(ka).
    #
    # Element-by-element source trace (all from Sano 2023 Table 2). The
    # comments are kept OUTSIDE the c(...) call below: a trailing `#`
    # comment inside an omega block's c(...) parses under source() but
    # breaks readModelDb()'s comment-to-label rewriter.
    #   0.214369 = 0.463^2  omega^2 CL/F   (46.3 %CV, RSE 7.98%, shrinkage 14.8%)
    #   0.298               cov(CL/F,Vd/F) (RSE 38.8%)
    #   1.2996   = 1.14^2   omega^2 Vd/F   (114 %CV, RSE 19.2%, shrinkage 42.1%)
    #   0.0815              cov(CL/F,ka)   (RSE 34.3%)
    #   0.435               cov(Vd/F,ka)   (RSE 26.7%)
    #   0.186624 = 0.432^2  omega^2 ka     (43.2 %CV, RSE 16.2%, shrinkage 45.0%)
    # ------------------------------------------------------------------
    etalcl + etalvc + etalka ~ c(
      0.214369,
      0.298,     1.2996,
      0.0815,    0.435,    0.186624
    )

    # ------------------------------------------------------------------
    # Residual unexplained variability. Sano 2023 Methods section 2.5 and
    # Online Resource 8c specify RUV in LOG-TRANSFORMED 5-HMT
    # concentration as purely additive:
    #     ln(Y_ij) = ln(F_ij) + eps_ij,   eps ~ N(0, sigma^2)
    # which is exactly a log-normal residual in linear space. The
    # NM-TRAN $ERROR block (Online Resource 4) implements it as
    # Y = IPRED + W*EPS(1) with IPRED = LOG(F), W = THETA(8) and $SIGMA
    # 1 FIX, so THETA(8) carries the whole residual magnitude.
    #
    # Encoded here as `lnorm(expSd)` rather than `prop(propSd)` because
    # the source form is exactly log-normal; the two differ materially at
    # this magnitude (a 0.381 log-scale SD is a 39.5% CV, not 38.1%).
    # ------------------------------------------------------------------
    expSd <- 0.381; label("Additive residual SD on the log-transformed concentration scale (log-normal)") # Sano 2023 Table 2 sigma: 38.1 %CV = 0.381 on the log scale (RSE 8.46%; bootstrap median 37.2, 95% CI 31.6-42.7; epsilon shrinkage 19.8%)
  })

  model({
    # Reference body weight for the allometric terms (Sano 2023 Methods
    # section 2.3: "referenced to 35 kg").
    wtRef <- 35

    # Individual parameters, assembled per Online Resource 8a. The
    # multiplicative covariate factors of the paper become additive terms
    # on the log scale: theta_TV * theta_cov^I  ==  exp(l_theta + e_cov*I).
    cl <- exp(lcl + e_wt_cl * log(WT / wtRef) + e_cyp2d6_pm_cl * CYP2D6_PM + e_sexf_cl * SEXF + etalcl)
    vc <- exp(lvc + e_wt_vc * log(WT / wtRef) + e_sexf_vc * SEXF + etalvc)
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)

    kel <- cl / vc

    # Bioavailability: the tablet arm is the structural F = 1 anchor and
    # the BIC arm carries the estimated relative bioavailability
    # (Online Resource 4: F1 = 1; IF(BIC.EQ.1) F1 = THETA(7)).
    f(depot) <- (1 - FORM_CAPSULE) + FORM_CAPSULE * exp(lfdepot)

    # First-order absorption lag on the depot (NM-TRAN ALAG1).
    alag(depot) <- tlag

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Doses are in micrograms and vc is in litres, so central/vc is in
    # ug/L == ng/mL, matching the source NM-TRAN scaling (S2 = V with
    # AMT in ug per the Online Resource 4 variable list).
    Cc <- central / vc

    Cc ~ lnorm(expSd)
  })
}
