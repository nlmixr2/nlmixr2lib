Lin_2026_sc1cmt_coxTte <- function() {
  description <- "Methodology reference. Data-generating time-to-event model from Lin 2026, the paper that shows how to compute a time-varying Cox partial likelihood inside NONMEM. There is no molecule and there are no patients: a hypothetical subcutaneous drug follows a one-compartment model with first-order absorption and elimination, and the event hazard is the product of a two-exponential bathtub-like baseline hazard and a hazard ratio driven by a time-invariant covariate, a time-varying covariate transform, and an Imax inhibitory effect of the simulated concentration. The Imax term sits on the LOG-hazard scale, so Imax = 2 is not bounded by 1 and saturating drug effect multiplies the hazard by exp(-2) = 0.135. Cumulative event hazard and cumulative censoring hazard are carried as ODE states and the survivor functions are exposed as derived outputs. All parameter values are author-chosen simulation constants, not estimates, and are therefore fixed()."
  reference <- paste(
    "Lin CW, Chen PW, Doshi S, Dutta S.",
    "Integration of Time-Varying Pharmacometric Modeling With Cox Regression",
    "for Time-to-Event Analysis in NONMEM.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(6):e70253.",
    "doi:10.1002/psp4.70253. PMCID PMC13106229. Open Access under CC BY-NC-ND.",
    "Equations 1-8 and the simulation parameter values are in Methods section 2.1",
    "(pages 2-3); the data-generating model is deposited verbatim as the",
    "$DES/$THETA blocks of Data S2 (psp470253-sup-0002-DataS2.mod), and an N = 20",
    "example dataset is deposited as Data S3 (psp470253-sup-0003-DataS3.csv).",
    "NOTE - PAPER ERRATUM: page 3 prints the label 'kh1' TWICE (0.03 and 0.01) and",
    "never prints 'kh2'; the second value is kh2. See the vignette for the proof.",
    sep = " "
  )
  vignette <- "Lin_2026_sc1cmt_coxTte"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = FALSE throughout, and deliberately so: the
  # analyte is confirmed (a hypothetical subcutaneous drug; Data S2 $DES gives
  # A(1) = depot, A(2) = central with C2 = A(2)/VC, A(3) = cumulative hazard),
  # but the SPECIMEN is not stated anywhere in the source. Lin 2026 invents the
  # drug and speaks only of "drug concentration", so "plasma" below is the
  # conventional reading of a one-compartment model, not something the paper
  # confirms. verified = TRUE requires analyte AND specimen from the source.
  compartmentData <- list(
    depot = list(analyte = "hypothetical subcutaneous drug", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "hypothetical subcutaneous drug", units = "mg", specimen = "plasma", verified = FALSE),
    cumhaz = list(analyte = "cumulative event hazard", units = NA_character_, specimen = "not applicable", verified = FALSE),
    cumhaz_cens = list(analyte = "cumulative censoring hazard", units = NA_character_, specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    SIMCOV_TI = list(
      description        = "Meaning-free simulation covariate entering the log hazard directly (Lin 2026 'COV1'). A per-subject standard-normal draw with no clinical interpretation, used to demonstrate a TIME-INVARIANT covariate effect on the event hazard.",
      units              = "(z-score; standard normal, mean 0, SD 1)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters Eq. 6 as beta1 * SIMCOV_TI with beta1 = 0.3. Lin 2026 Methods 2.1: 'COV1 and COV2 are generated from a standard normal distribution with mean of 0 and a standard deviation of 1.' Deposited as the constant-per-subject column COV1 in Data S3.",
      source_name        = "COV1"
    ),
    SIMCOV_TV = list(
      description        = "Meaning-free simulation covariate entering the log hazard through a time-varying transform (Lin 2026 'COV2'). A per-subject standard-normal draw with no clinical interpretation, used to demonstrate a TIME-VARYING covariate effect on the event hazard.",
      units              = "(z-score; standard normal, mean 0, SD 1)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The DATA column is time-fixed per subject; only the model-side transform COVT2(t) = SIMCOV_TV * log(t + 20) is time-varying (Lin 2026 Eq. 5). That is precisely the paper's demonstration -- a time-invariant draw given a time-varying effect. Enters Eq. 6 as beta2 * COVT2(t) with beta2 = 0.05. Lin 2026 Methods 2.1 calls the time function 'an arbitrary logarithmic time function for demonstration purposes' and states that 'the additive constant of 20 was included to ensure stability at early time points and has no clinical interpretation.' Deposited as the constant-per-subject column COV2 in Data S3.",
      source_name        = "COV2"
    )
  )

  population <- list(
    species        = "None (methodology paper; simulation-only data-generating model with no drug, no patients, and no fitted estimates).",
    n_subjects     = 1500L,
    n_studies      = 1L,
    disease_state  = "N/A (Monte Carlo clinical-trial simulation study; not a fit of any real molecule).",
    dose_range     = "Five equally sized dose groups given 1, 3, 10, 30 and 100 mg subcutaneously every 4 weeks for 24 weeks (six doses at days 0, 28, 56, 84, 112 and 140), with follow-up to 96 weeks (672 days). Events after 96 weeks were treated as censored (Lin 2026 Methods 2.1).",
    regions        = "N/A",
    scope_note     = paste(
      "Filed under inst/modeldb/pharmacodynamics/ (not specificDrugs/) because there is no drug;",
      "the registry's time-to-event / hazard models live in this category. The file stem follows the",
      "established hypothetical-drug family <route><n>cmt_<topic> set by",
      "inst/modeldb/pharmacokinetics/Beal_2001_iv1cmt_bql.R and the Schoning_2026_oral1cmt_* family.",
      "The paper is a methodology reference introducing a way to evaluate the Cox partial likelihood",
      "inside NONMEM; every value below is a simulation constant the authors chose, so all are fixed()."
    ),
    notes          = paste(
      "The paper deposits TWO NONMEM control streams, but they are two ESTIMATION METHODS applied to",
      "one data-generating model, not two models. Data S1 is the semi-parametric Cox partial-likelihood",
      "stream; it deliberately carries no baseline hazard and therefore cannot be simulated forward.",
      "Data S2 is the parametric two-exponential bathtub stream, which IS the data-generating model and",
      "is what this file encodes. In both streams KA, KE and VC are $INPUT data columns with $OMEGA 0 FIX,",
      "so the individual PK parameters are supplied per subject rather than estimated and the dummy ETA",
      "exists only to keep NONMEM happy. The inter-individual variability therefore belongs to the",
      "data-generating model of Methods 2.1, which is what is encoded here.",
      "Simulations were run at N = 20, 50, 150, 500 and 1500 subjects with 100 replicates each;",
      "n_subjects records the largest.",
      "The true parameter values beta1 = 0.3, beta2 = 0.05, Imax = 2 and IC50 = 8 mg/L are also printed",
      "as the 'true parameter values' footnotes to Tables 1 and 2."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # PK. Lin 2026 Methods 2.1 page 3:
    #   "Ka (1/day) = 0.3, Ke (1/day) ~ log normal(log(0.03), 0.5),
    #    V2 (L) ~ log normal(log(3), 0.5)"
    # Author-chosen simulation constants, hence fixed().
    # ------------------------------------------------------------------
    lka <- fixed(log(0.3));  label("Absorption rate constant Ka (1/day)")                   # Lin 2026 p. 3 (Ka = 0.3 1/day); Data S3 column KA is constant across all 20 subjects, confirming no IIV
    lkel <- fixed(log(0.03)); label("Elimination rate constant Ke (1/day)")                 # Lin 2026 p. 3 (median Ke = 0.03 1/day)
    lvc <- fixed(log(3));    label("Central volume of distribution V2 (L)")                 # Lin 2026 p. 3 (median V2 = 3 L)

    # The second argument of "log normal(log(m), 0.5)" is the SD on the log
    # scale (R rlnorm sdlog), NOT a variance, so the variance is 0.5^2 = 0.25.
    # Settled empirically against the deposited per-subject KE / VC values in
    # Data S3: pooled sd(log) = 0.4278 on 38 df gives a lower-tail chi-square
    # p of 0.112 under sdlog = 0.5 versus 0.00012 under variance = 0.5. See
    # the vignette's Assumptions and deviations section.
    etalkel ~ fixed(0.25)  # variance = 0.5^2; Lin 2026 p. 3 Ke ~ log normal(log(0.03), 0.5) with 0.5 read as the log-scale SD
    etalvc ~ fixed(0.25)   # variance = 0.5^2; Lin 2026 p. 3 V2 ~ log normal(log(3), 0.5) with 0.5 read as the log-scale SD

    # ------------------------------------------------------------------
    # Two-exponential bathtub-like baseline hazard (Lin 2026 Eq. 4):
    #   Hbaseline(t) = h1 * exp(-kh1 * t) + h2 * exp(+kh2 * t)
    # Note the PLUS sign on the second exponent -- that rising second term
    # is what makes the sum bathtub-shaped. Confirmed verbatim in the
    # Data S2 $DES block: BSHAZ = (AMP1*EXP(-HK1*T) + AMP2*EXP(HK2*T)).
    # Log-transformed here because all four are strictly positive; the
    # Data S2 $THETA block likewise bounds each below at 0.
    # ------------------------------------------------------------------
    lh1_haz <- fixed(log(0.02));    label("Baseline-hazard amplitude of the decaying (infant-mortality) term h1 (1/day)")  # Lin 2026 p. 3 h1 = 0.02 1/day; Data S2 $THETA 5 AMP1 = 0.02
    lkh1_haz <- fixed(log(0.03));   label("Baseline-hazard decay rate constant of the first term kh1 (1/day)")             # Lin 2026 p. 3 kh1 = 0.03 1/day; Data S2 $THETA 6 HK1 = 0.03
    lh2_haz <- fixed(log(0.0005));  label("Baseline-hazard amplitude of the rising (wear-out) term h2 (1/day)")            # Lin 2026 p. 3 h2 = 0.0005 1/day; Data S2 $THETA 7 AMP2 = 0.0005
    lkh2_haz <- fixed(log(0.01));   label("Baseline-hazard rise rate constant of the second term kh2 (1/day)")             # Data S2 $THETA 8 HK2 = 0.01 -- p. 3 mislabels this value 'kh1' a second time and never prints kh2; see the vignette Errata

    # Constant censoring hazard (Lin 2026 Eq. 8).
    lhcens_haz <- fixed(log(0.001)); label("Constant censoring hazard hcensor (1/day)")                                    # Lin 2026 p. 3 hcensor = 0.001 1/day

    # ------------------------------------------------------------------
    # Hazard-ratio terms (Lin 2026 Eq. 6):
    #   HR(t) = exp(beta1*COV1 + beta2*COVT2(t) - Imax*C2(t)/(IC50 + C2(t)))
    # ------------------------------------------------------------------
    e_simcov_ti_haz <- fixed(0.3);  label("Hazard log-linear coefficient beta1 on the time-invariant simulation covariate (per z-score unit)")  # Lin 2026 p. 3 beta1 = 0.3; Data S1/S2 $THETA 1 EFF_COV1 = 0.3; also the "true parameter values" footnote to Tables 1 and 2
    e_simcov_tv_haz <- fixed(0.05); label("Hazard log-linear coefficient beta2 on the time-varying covariate transform SIMCOV_TV * log(t + 20)") # Lin 2026 p. 3 beta2 = 0.05; Data S1/S2 $THETA 2 EFF_COV2 = 0.05; also the Tables 1-2 footnote

    # Imax is NOT bounded by 1: the inhibition term is subtracted on the
    # LOG-hazard scale, so a saturating effect multiplies the hazard by
    # exp(-Imax) = exp(-2) = 0.135.
    limax <- fixed(log(2));  label("Maximum drug inhibition of the log hazard, Imax (unitless; multiplies the hazard by exp(-Imax) at saturation)")  # Lin 2026 p. 3 Imax = 2 (the Data S2 $THETA 3 LN_IMAX entry is an estimation START value, not the truth)
    lic50 <- fixed(log(8));  label("Concentration producing half of Imax, IC50 (mg/L)")                                                             # Lin 2026 p. 3 IC50 = 8 mg/L (the Data S2 $THETA 4 LN_IC50 entry is an estimation START value, not the truth)

    # ------------------------------------------------------------------
    # Residual error placeholder. The source NONMEM runs use -2LL on the
    # survival / partial likelihood, so there is no observation-error
    # model. This nlmixr2 translation is intended for forward simulation:
    # `hazard`, `cumhaz`, `sur`, `hazard_cens` and `sur_cens` are exposed
    # as derived outputs and a tiny additive residual is attached to `sur`
    # so the nlmixr2 likelihood machinery accepts the model.
    # ------------------------------------------------------------------
    addSd <- fixed(0.001); label("Placeholder additive residual error on the survivor-function output `sur` (unitless); NOT from the source -- see the vignette Assumptions")
  })

  model({
    # ---- PK: one-compartment, first-order absorption (Eqs. 1-3) -------
    ka <- exp(lka)
    kel <- exp(lkel + etalkel)
    vc <- exp(lvc + etalvc)

    d/dt(depot) <- -ka * depot                  # Eq. 1
    d/dt(central) <- ka * depot - kel * central # Eq. 2

    Cc <- central / vc                          # Eq. 3

    # ---- Two-exponential bathtub baseline hazard (Eq. 4) -------------
    h1 <- exp(lh1_haz)
    kh1 <- exp(lkh1_haz)
    h2 <- exp(lh2_haz)
    kh2 <- exp(lkh2_haz)

    hazbase <- h1 * exp(-kh1 * t) + h2 * exp(kh2 * t)

    # ---- Time-varying covariate transform (Eq. 5) --------------------
    # The +20 offset is the paper's stability constant and carries no
    # clinical meaning.
    covt2 <- SIMCOV_TV * log(t + 20)

    # ---- Hazard ratio (Eq. 6) ----------------------------------------
    imax <- exp(limax)
    ic50 <- exp(lic50)

    hr <- exp(e_simcov_ti_haz * SIMCOV_TI +
                e_simcov_tv_haz * covt2 -
                imax * Cc / (ic50 + Cc))

    # ---- Event hazard, cumulative hazard, survivor function (Eq. 7) --
    hazard <- hazbase * hr
    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)

    # ---- Constant censoring hazard (Eq. 8) ---------------------------
    hazard_cens <- exp(lhcens_haz)
    d/dt(cumhaz_cens) <- hazard_cens
    cumhaz_cens(0) <- 0
    sur_cens <- exp(-cumhaz_cens)

    sur ~ add(addSd)
  })
}
