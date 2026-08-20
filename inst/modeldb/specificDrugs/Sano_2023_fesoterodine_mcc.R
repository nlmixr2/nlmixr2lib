Sano_2023_fesoterodine_mcc <- function() {
  description <- "Population pharmacokinetic/pharmacodynamic exposure-response model relating 5-hydroxymethyl tolterodine (5-HMT, the active metabolite of fesoterodine) average steady-state plasma concentration to maximum cystometric capacity (MCC) in pediatric patients aged 6-16 years with neurogenic detrusor overactivity. An Emax model in which the maximum attainable MCC is the age-based expected bladder capacity (EBC) rather than an estimated parameter: MCC = BASE + (Emax - BASE) * Cavg,ss / (EC50 + Cavg,ss), with Emax = 30 * (AGE + 1) mL up to age 12 and a 390 mL plateau thereafter (FIXED, not estimated), and typical baseline MCC scaling with age by the same (AGE + 1)/13 factor from a 190 mL plateau. Fitted to 242 MCC observations (baseline and week 12) from 121 patients in the phase III study 1047 (NCT01557244); EC50 is 6.22 ng/mL. This is the pharmacodynamic half of a sequential PK-then-PD analysis: Cavg,ss is supplied as the CAV covariate from individual empirical-Bayes estimates of the companion population PK model, available in this library as modellib('Sano_2023_fesoterodine'). This is a PD-only model with no dose events and no ODE states; correlated inter-individual variability is carried on baseline MCC and Emax, and residual error is combined proportional plus additive."
  reference <- "Sano Y, Shoji S, Shahin M, Sweeney K, Darekar A, Malhotra BK. Population Pharmacokinetic and Pharmacodynamic Modeling of Fesoterodine in Pediatric Patients with Neurogenic Detrusor Overactivity. Eur J Drug Metab Pharmacokinet. 2023 May;48(3):257-269. doi:10.1007/s13318-023-00818-8. PMID: 36892805. PMCID: PMC10175358."
  vignette <- "Sano_2023_fesoterodine"
  units <- list(time = "week", dosing = "n/a (PD-only model; no dose events)", concentration = "ng/mL", response = "mL")

  covariateData <- list(
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives BOTH the maximum attainable MCC and the typical baseline MCC through the same piecewise-linear factor (AGE + 1)/13, capped at 1 for AGE > 12 (Sano 2023 Eqs. 5-8 and Online Resource 8b). The factor comes from the pediatric expected-bladder-capacity (EBC) rule quoted in Sano 2023 Methods section 2.3.3: EBC = [30 + AGE * 30] mL = 30 * (AGE + 1) mL up to age 12, constant at 390 mL from age 12 onwards. Age is used at its baseline value and is not advanced over the 12-week observation window. Observed range in the pharmacokinetic/pharmacodynamic analysis population is 6-16 years (median 9).",
      source_name        = "AGE"
    ),
    CAV = list(
      description        = "Individual 5-HMT average plasma concentration at steady state (Cavg,ss) over the once-daily dosing interval",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cavg,ss = F * DOSE / (CL/F * tau) (Sano 2023 Eq. 4), computed per patient from the individual empirical-Bayes CL/F of the companion population PK model (modellib('Sano_2023_fesoterodine')) with tau = 24 h. In the source NM-TRAN dataset this is the CSS column, supplied as data rather than computed inside the PD model, because Sano 2023 fitted the PK and PD models SEQUENTIALLY (Online Resource 5 is a $PRED model with no PK). SET CAV = 0 FOR THE BASELINE OCCASION: Online Resource 8b states explicitly that Cavg,ss at baseline was set to zero, which collapses the Emax term and makes the baseline prediction exactly BASE. Observed group medians in study 1047 were 1.21, 2.64, 2.49 and 4.22 ng/mL for 2 mg BIC, 4 mg BIC, 4 mg tablet and 8 mg tablet QD respectively (Online Resource 10).",
      source_name        = "CSS"
    )
  )

  population <- list(
    species           = "human",
    n_subjects        = 121L,
    n_observations    = 242L,
    n_studies         = 1L,
    age_range         = "6-16 years (median 9, mean 9.58, SD 2.72)",
    weight_range      = "11.7-85.0 kg (median 28, mean 33.5, SD 14.4)",
    sex_female_pct    = 49.6,
    race_ethnicity    = c(White = 44.6, Black = 1.7, Asian = 51.2, Other = 2.5),
    cyp2d6_pm_pct     = 1.7,
    baseline_mcc      = "Median 152 mL (range 16-451); mean 165 mL (SD 94.6)",
    disease_state     = "Pediatric patients with symptoms of neurogenic detrusor overactivity (NDO) enrolled in the phase III study 1047.",
    dose_range        = "Fesoterodine 4 or 8 mg tablet QD (cohort 1, patients over 25 kg); fesoterodine 2 or 4 mg beads-in-capsule QD (cohort 2, patients 25 kg or less). Patients randomized to the oxybutynin extended-release comparator arm are not part of this analysis population.",
    regions           = "Multinational (NCT01557244 phase III study 1047).",
    notes             = "Demographics from Sano 2023 Online Resource 11 (pharmacokinetic/pharmacodynamic analysis population). Each patient contributes exactly two MCC observations: baseline and week 12 (2 x 121 = 242). MCC was measured by multichannel cystometry through a dual-lumen urodynamic catheter, filling until voiding or leaking began, until MCC reached the expected bladder capacity, or at a detrusor pressure of at least 40 cm H2O (Sano 2023 Methods section 2.2.3). Because few observations lay above EC50, the exposure-response relationship is close to linear across the observed Cavg,ss range and the authors caution that EC50 carries wide uncertainty (95% CI 4.11-10.1 ng/mL)."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. Values are the final population estimates in
    # Sano 2023 Table 3 ("Value" column). BASE and Emax as tabulated are
    # the values for a patient OLDER THAN 12 years; the age scaling for
    # younger patients is applied inside model() per Table 3 footnotes a
    # and b, Eqs. 5-8, and Online Resource 8b.
    # ------------------------------------------------------------------
    lrbase <- log(190);        label("Typical baseline MCC plateau for age > 12 years (mL)")  # Sano 2023 Table 3 BASE: 190 mL (RSE 5.49%; bootstrap median 190, 95% CI 170-211)
    lec50  <- log(6.22);       label("5-HMT Cavg,ss producing half the maximum response (ng/mL)") # Sano 2023 Table 3 EC50: 6.22 ng/mL (RSE 21.8%; bootstrap median 6.21, 95% CI 4.11-10.1)

    # Emax was NOT estimated. Sano 2023 Methods section 2.4 assumes the
    # maximum attainable MCC equals the age-based expected bladder
    # capacity, "because Emax can be interpreted as the bladder capacity
    # when individuals completely improve their symptom". Table 3 reports
    # it as "390, Fixed" for age > 12; the NM-TRAN $PRED block hardcodes
    # EMAX = (AGE+1)*30 for AGE <= 12 and 13*30 = 390 otherwise, with no
    # $THETA. It nevertheless carries an ESTIMATED inter-individual
    # variance (etalemax below), so the typical value is fixed while the
    # between-patient spread around it is not.
    lemax  <- fixed(log(390)); label("Maximum attainable MCC (expected bladder capacity) plateau for age > 12 years (mL), from source") # Sano 2023 Table 3 Emax: "390, Fixed"; = 30 * (12 + 1) from the EBC rule in Methods section 2.2.3

    # ------------------------------------------------------------------
    # Inter-individual variability -- OMEGA BLOCK(2) on BASE and Emax
    # (Online Resource 5 $OMEGA BLOCK(2)).
    #
    # As in the companion PK model, Sano 2023 Table 3 reports the
    # DIAGONAL entries as "% CV" under the sqrt-of-variance convention
    # (omega^2 = (%CV/100)^2) and the OFF-DIAGONAL entry as a bare
    # covariance. Checks: the implied variances (0.238, 0.222) sit on the
    # Online Resource 5 $OMEGA initial estimates (0.3, 0.3), and the
    # implied correlation 0.122/sqrt(0.238144*0.221841) = 0.53 is valid
    # (it would also be valid, but less consistent with the initials,
    # under the log-normal convention).
    #
    # Block order is the NONMEM lower-triangular row order:
    #   var(BASE), cov(BASE,Emax), var(Emax).
    #
    # Element-by-element source trace (all from Sano 2023 Table 3). The
    # comments are kept OUTSIDE the c(...) call below: a trailing `#`
    # comment inside an omega block's c(...) parses under source() but
    # breaks readModelDb()'s comment-to-label rewriter.
    #   0.238144 = 0.488^2  omega^2 BASE   (48.8 %CV, RSE 7.89%, shrinkage 10.4%)
    #   0.122               cov(BASE,Emax) (RSE 28.3%)
    #   0.221841 = 0.471^2  omega^2 Emax   (47.1 %CV, RSE 14.3%, shrinkage 22.7%)
    # ------------------------------------------------------------------
    etalrbase + etalemax ~ c(
      0.238144,
      0.122,    0.221841
    )

    # ------------------------------------------------------------------
    # Residual unexplained variability -- combined proportional plus
    # additive on MCC (Sano 2023 Methods section 2.4 and Online Resource
    # 8d): Y_ij = F_ij * (1 + eps_PRP,ij) + eps_ADD,ij. The NM-TRAN
    # $ERROR equivalent (Online Resource 5) is
    # Y = E*(1 + THETA(3)*ERR(1)) + THETA(4)*ERR(2) with $SIGMA 1 FIX on
    # both, so THETA(3) and THETA(4) carry the residual magnitudes
    # directly as SDs.
    # ------------------------------------------------------------------
    propSd_MCC <- 0.0741; label("Proportional residual SD on MCC (fraction)")  # Sano 2023 Table 3 sigma_PRP: 7.41 %CV = 0.0741 (RSE 36.4%; bootstrap median 7.25, 95% CI 0.0729-12.5; epsilon shrinkage 35.3%)
    addSd_MCC  <- 34.6;   label("Additive residual SD on MCC (mL)")            # Sano 2023 Table 3 sigma_ADD: 34.6 mL (RSE 11.1%; bootstrap median 34, 95% CI 26.5-42.1; epsilon shrinkage 35.3%)
  })

  model({
    # Piecewise-linear age scaling shared by BASE and Emax, from the
    # pediatric expected-bladder-capacity rule (Sano 2023 Eqs. 5-8):
    #   Emax = 30 * (AGE + 1) mL   if AGE <= 12,   390 mL if AGE > 12
    #   BASE = theta_base * (AGE + 1)/13 if AGE <= 12, theta_base if AGE > 12
    # Both collapse to <plateau> * ageScale with the same factor, since
    # 390 * (AGE + 1)/13 == 30 * (AGE + 1).
    if (AGE <= 12) {
      ageScale <- (AGE + 1) / 13
    } else {
      ageScale <- 1
    }

    rbase <- exp(lrbase + etalrbase) * ageScale
    emax  <- exp(lemax  + etalemax)  * ageScale
    ec50  <- exp(lec50)

    # Emax exposure-response (Sano 2023 Eq. 3 / Online Resource 8b).
    # Emax here is the maximum ATTAINABLE MCC (the bladder-capacity
    # ceiling), not the maximum increment, so the drug effect is the
    # fraction (emax - rbase) of the remaining headroom above baseline.
    # CAV = 0 on the baseline occasion => MCC = rbase exactly.
    MCC <- rbase + (emax - rbase) * CAV / (ec50 + CAV)

    MCC ~ prop(propSd_MCC) + add(addSd_MCC)
  })
}
