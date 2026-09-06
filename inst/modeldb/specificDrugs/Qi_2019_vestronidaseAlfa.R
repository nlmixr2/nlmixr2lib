Qi_2019_vestronidaseAlfa <- function() {
  description <- "Joint population PK + inhibitory-Imax-on-AUC exposure-response model for vestronidase alfa (recombinant human beta-glucuronidase, UX003) enzyme replacement therapy in patients with mucopolysaccharidosis type VII (Sly syndrome). The PK layer (Table 2) is a two-compartment model with zero-order intravenous input and linear elimination from the central compartment, carrying baseline body weight as the only retained covariate via two estimated allometric exponents (a shared exponent on CL and Q, and a shared exponent on Vc and Vp), both centered on a 20 kg reference weight. The PD layer (Table 4) is a static inhibitory Emax model linking the individual vestronidase alfa AUC over a dosing interval to the percent change from pretreatment baseline in urinary chondroitin sulfate (uCS) and urinary dermatan sulfate (uDS), fit sequentially on post hoc individual exposures. AUC is integrated inside an extra rxode2 state, so both uGAG observables are meaningful at the end of a dosing interval (the vignette documents the time-window discipline)."
  reference <- paste(
    "Qi Y, Mc Namara MP, Haller C, Song W, Gutierrez F, Kolodny E, Ma J.",
    "Pharmacokinetic and pharmacodynamic modeling to optimize the dose of",
    "vestronidase alfa, an enzyme replacement therapy for treatment of patients",
    "with mucopolysaccharidosis type VII: results from three trials.",
    "Clin Pharmacokinet. 2019;58(5):673-683. doi:10.1007/s40262-018-0721-y.",
    "Erratum: Clin Pharmacokinet. 2019;58(5):685.",
    "doi:10.1007/s40262-018-0726-6 (corrects a ClinicalTrials.gov identifier in",
    "the Introduction; no parameter value is affected).",
    sep = " "
  )
  vignette <- "Qi_2019_vestronidaseAlfa"
  paper_specific_compartments <- c("auc_tau")
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final population PK model (Qi 2019 Results;",
        "neutralizing antibody status, ADA titer, other demographics and laboratory",
        "tests were screened in a full covariate model with backward deletion at",
        "p = 0.001 and none was retained). Entered as a power (allometric) term",
        "(WT / 20)^exponent on all four disposition parameters, centered on a 20 kg",
        "reference weight (Table 2 footnote b). Two exponents were ESTIMATED rather",
        "than fixed at the canonical 0.75 / 1 allometric values (Table 2 footnote c):",
        "0.587 shared by CL and Q, and 0.483 shared by Vc and Vp. Baseline (time-fixed)",
        "body weight; the cohort median of 24.8 kg is the value the paper used for its",
        "dosing-regimen simulations."
      ),
      source_name        = "BWT"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "vestronidase alfa", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "vestronidase alfa", units = "mg", specimen = "serum", verified = TRUE),
    auc_tau     = list(analyte = "vestronidase alfa AUC over a dosing interval", units = "mg*h/L", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 23L,
    n_studies      = 3L,
    age_range      = "1.7-25.3 years",
    age_median     = "10.1 years (mean 11.0, SD 7.7)",
    weight_range   = "9.41-104.0 kg",
    weight_median  = "24.8 kg (mean 37.7, SD 27.7)",
    sex_female_pct = 52,
    race_ethnicity = c(White = 61, Asian = 13, Black = 4, Other = 22),
    disease_state  = "Mucopolysaccharidosis type VII (MPS VII, Sly syndrome), an ultra-rare lysosomal storage disorder caused by beta-glucuronidase deficiency.",
    dose_range     = "Vestronidase alfa 1, 2, or 4 mg/kg by intravenous infusion every other week (QOW), infused over approximately 4 h with the first 2.5% of the total volume given over the first hour.",
    regions        = "Multinational (three clinical trials: NCT01856218, NCT02418455, NCT02230566).",
    sampling       = "15 subjects contributed intensively sampled PK profiles and 8 contributed sparse profiles (Qi 2019 Results). Urinary CS and DS were assayed at representative steady-state visits by LC-MS/MS and creatinine-corrected.",
    notes          = paste(
      "Demographics from Qi 2019 Table 1 (overall column across the three studies:",
      "study 1 phase I/II, n = 3, ages 5.5-25.1 y; study 2 phase II, n = 8, ages",
      "1.7-5.0 y; study 3 phase III, n = 12, ages 8.5-25.3 y). Race percentages are",
      "the Table 1 overall column and sum to 100. The same 23 subjects support both",
      "the PK and the exposure-response layers; PK and PD were fit SEQUENTIALLY",
      "(individual post hoc AUC from the PK model was carried into the PD fit as the",
      "exposure predictor), not simultaneously."
    )
  )

  ini({
    # =========================================================================
    # PK structural parameters (Table 2; typical values for a 20 kg subject)
    # =========================================================================
    # Table 2 footnote a: CL, Vc, Q and Vp were estimated on a log scale and
    # subsequently back-transformed, so the printed point estimates are already
    # on the linear scale and log() is applied here to recover the estimated
    # parameter.
    lcl <- log(1.97);  label("Serum clearance for a 20 kg subject (CL, L/h)")                          # Table 2: CL = 1.97 L/h (RSE 9.61%; 95% CI 1.63, 2.38; bootstrap median 1.99)
    lvc <- log(1.52);  label("Central volume of distribution for a 20 kg subject (Vc, L)")             # Table 2: Vc = 1.52 L (RSE 9.64%; 95% CI 1.26, 1.83; bootstrap median 1.49)
    lq  <- log(0.931); label("Intercompartmental clearance for a 20 kg subject (Q, L/h)")              # Table 2: Q = 0.931 L/h (RSE 16.3%; 95% CI 0.676, 1.28; bootstrap median 0.940)
    lvp <- log(3.11);  label("Peripheral volume of distribution for a 20 kg subject (Vp, L)")          # Table 2: Vp = 3.11 L (RSE 4.23%; 95% CI 2.86, 3.38; bootstrap median 3.03)

    # ---- Allometric body-weight exponents (Table 2, "Covariates") ----
    # Both are ESTIMATED, not fixed: Table 2 footnote c states explicitly that
    # the reported confidence intervals "do not contain allometric exponents of
    # 0.75 and 1 for CL and V, respectively", i.e. the estimated exponents are
    # statistically distinct from the canonical allometric values. They are
    # therefore NOT wrapped in fixed(). Each exponent is SHARED by two
    # parameters, which is what the e_<cov>_<param1>_<param2> naming records.
    e_wt_cl_q  <- 0.587; label("Power exponent of (WT / 20 kg) shared by CL and Q (unitless)")         # Table 2: BWT on CL and Q = 0.587 (RSE 12.4%; 95% CI 0.444, 0.730; bootstrap median 0.574)
    e_wt_vc_vp <- 0.483; label("Power exponent of (WT / 20 kg) shared by Vc and Vp (unitless)")        # Table 2: BWT on Vc and Vp = 0.483 (RSE 17.7%; 95% CI 0.315, 0.651; bootstrap median 0.492)

    # ---- PK interindividual variability (Table 2, "Interindividual variability") ----
    # Exponential IIV, theta_i = theta_TV * exp(eta_i) (Qi 2019 Methods, Eq. 1).
    # Table 2 reports omega^2 DIRECTLY, so these are transcribed as-is and no
    # CV% back-transformation is performed. This matters here: the paper's CV%
    # column switches formula row by row (footnote d gives
    # CV = sqrt(exp(omega^2) - 1) only when omega^2 > 0.15, and a plain
    # sqrt(omega^2) otherwise), so inverting that column with a single formula
    # would misstate the Vc and Vp variances. Reported CV%: CL 43.0, Vc 30.7,
    # Q 86.9, Vp 15.4. The paper reports no off-diagonal covariances, so the
    # block is encoded as diagonal.
    etalcl ~ 0.170                                                                                     # Table 2: omega^2 for CL = 0.170 (RSE 50.2%; 95% CI 0.00281, 0.337)
    etalvc ~ 0.0944                                                                                    # Table 2: omega^2 for Vc = 0.0944 (RSE 54.0%)
    etalq  ~ 0.563                                                                                     # Table 2: omega^2 for Q = 0.563 (RSE 52.4%)
    etalvp ~ 0.0236                                                                                    # Table 2: omega^2 for Vp = 0.0236 (RSE 39.9%; 95% CI 0.00516, 0.0420)

    # ---- PK residual error (Table 2, "Residual variability") ----
    # Proportional error selected over combined proportional-plus-additive
    # (Qi 2019 Methods): DV = IPRED * (1 + eps), var(eps) = sigma^2 = 0.109.
    # propSd is a standard deviation, so sqrt(0.109) = 0.330; the paper's own
    # CV column prints 33.0%, which confirms the square-root convention.
    propSd <- sqrt(0.109); label("Proportional residual error on serum vestronidase alfa (fraction)")   # Table 2: sigma^2 proportional = 0.109 (RSE 18.7%; 95% CI 0.0690, 0.149; CV 33.0%)

    # =========================================================================
    # PD / exposure-response parameters (Table 4, inhibitory Imax on AUC)
    # =========================================================================
    # Qi 2019 Exposure-Response Analysis, displayed equation:
    #     I = I0 - Imax * AUC / (AUC + IC50)
    # I is the percent CHANGE FROM PRETREATMENT BASELINE in urinary GAG, so a
    # maximally suppressed subject approaches -Imax percent. The leading minus
    # sign is the paper's own (verified against the JATS MathML of the display
    # equation and confirmed by forward-simulating the curve against Figure 4:
    # at AUC = 110 ug*h/mL the equation gives -76.0% for uCS and -72.1% for
    # uDS, matching the plotted population-prediction curves).
    #
    # I0 is the no-drug response and was FIXED to 0 by the authors. It is kept
    # here as an explicit fixed(0) so the model form remains auditable rather
    # than silently dropping a term the paper's equation carries. It is NOT log
    # transformed: the value is exactly zero and the quantity is a percent
    # change that is not sign constrained, so the usual `l` prefix does not apply.
    rbase_ucs <- fixed(0); label("No-drug response for uCS, percent change from baseline (%)")         # Qi 2019 Methods (Exposure-Response Analysis): "I0 is the no-drug response and fixed to 0 in the model"
    rbase_uds <- fixed(0); label("No-drug response for uDS, percent change from baseline (%)")         # Qi 2019 Methods (Exposure-Response Analysis): "I0 is the no-drug response and fixed to 0 in the model"

    limax_ucs <- log(82.0); label("Maximal percent reduction from baseline in urinary chondroitin sulfate (Imax, %)") # Table 4: Imax uCS = 82.0% (RSE 5.3%)
    limax_uds <- log(76.9); label("Maximal percent reduction from baseline in urinary dermatan sulfate (Imax, %)")    # Table 4: Imax uDS = 76.9% (RSE 4.1%)
    lic50_ucs <- log(8.6);  label("Vestronidase alfa AUC producing half-maximal uCS reduction (IC50, ug*h/mL)")       # Table 4: IC50 uCS = 8.6 ug*h/mL (RSE 31.0%)
    lic50_uds <- log(7.3);  label("Vestronidase alfa AUC producing half-maximal uDS reduction (IC50, ug*h/mL)")       # Table 4: IC50 uDS = 7.3 ug*h/mL (RSE 24.5%)

    # ---- PD interindividual variability (Table 4) ----
    # The paper states IIV "was estimable for Imax only"; no IIV is reported on
    # IC50, and no off-diagonal covariance is reported between the two
    # endpoints, so both etas are diagonal. Table 4 reports %CV (not omega^2),
    # so a back-transformation is required. Applying the paper's own threshold
    # convention from Table 2 footnote d -- log-normal CV only when
    # omega^2 > 0.15, plain sqrt(omega^2) below it -- these variances are far
    # below the threshold, so omega^2 = CV^2:
    #   uCS: 0.079^2 = 0.006241  (log-normal alternative log(0.079^2 + 1) = 0.006222)
    #   uDS: 0.087^2 = 0.007569  (log-normal alternative log(0.087^2 + 1) = 0.007540)
    # The two conventions agree to within 0.3% of the variance at this
    # magnitude, so the choice is immaterial to any simulation.
    etalimax_ucs ~ 0.006241                                                                            # Table 4: IIV for Imax uCS = 7.9 %CV (RSE 35.9%)
    etalimax_uds ~ 0.007569                                                                            # Table 4: IIV for Imax uDS = 8.7 %CV (RSE 18.9%)

    # ---- PD residual error (Table 4) ----
    # Additive residual on the percent-change scale (the dependent variable of
    # the exposure-response model is a percent change from baseline), so the
    # units are percentage points, NOT ug*h/mL.
    addSd_uCSchange <- 6.1; label("Additive residual error on uCS percent change from baseline (%)")    # Table 4: additive residual error uCS = 6.1 (RSE 13.8%)
    addSd_uDSchange <- 6.4; label("Additive residual error on uDS percent change from baseline (%)")    # Table 4: additive residual error uDS = 6.4 (RSE 9.4%)
  })

  model({
    # =========================================================================
    # 1. Individual PK parameters
    # =========================================================================
    # Allometric scaling on baseline body weight, centered on the 20 kg
    # reference (Table 2 footnote b). One exponent is shared by the two
    # clearances and a second by the two volumes, exactly as the paper reports
    # them ("BWT on CL and Q" / "BWT on Vc and Vp").
    allom_cl <- (WT / 20)^e_wt_cl_q
    allom_v  <- (WT / 20)^e_wt_vc_vp

    cl <- exp(lcl + etalcl) * allom_cl
    q  <- exp(lq  + etalq)  * allom_cl
    vc <- exp(lvc + etalvc) * allom_v
    vp <- exp(lvp + etalvp) * allom_v

    # =========================================================================
    # 2. Two-compartment disposition with zero-order (infusion) input
    # =========================================================================
    # Dose is administered into `central` as an infusion (the vignette event
    # table supplies the infusion duration / rate); there is no absorption
    # compartment. Elimination is linear from the central compartment -- the
    # paper tested combined linear plus Michaelis-Menten elimination and found
    # no statistical improvement over the nested linear-only model.
    Cc <- central / vc

    d/dt(central)     <- -cl * Cc - q * Cc + q * (peripheral1 / vp)
    d/dt(peripheral1) <-            q * Cc - q * (peripheral1 / vp)

    # =========================================================================
    # 3. Exposure metric driving the PD layer
    # =========================================================================
    # The exposure-response model is driven by the individual AUC of serum
    # vestronidase alfa. This state integrates Cc from the start of the solve,
    # so at the end of one dosing interval auc_tau holds AUC(0-tau); for a
    # single dose carried to completion it converges to AUC(0-inf) = dose / CL.
    # Because the uGAG observables below are algebraic functions of an
    # accumulating integral, they are only interpretable at the END of a dosing
    # interval -- observing them mid-interval reports the effect of the
    # partial AUC accrued so far. The vignette documents this discipline.
    d/dt(auc_tau) <- Cc

    # =========================================================================
    # 4. Inhibitory Imax exposure-response on urinary GAG (Table 4)
    # =========================================================================
    #     I = I0 - Imax * AUC / (AUC + IC50)
    # with I0 fixed to 0, so each observable is the percent CHANGE from
    # pretreatment baseline and is therefore negative (a reduction).
    imax_ucs <- exp(limax_ucs + etalimax_ucs)
    imax_uds <- exp(limax_uds + etalimax_uds)
    ic50_ucs <- exp(lic50_ucs)
    ic50_uds <- exp(lic50_uds)

    uCSchange <- rbase_ucs - imax_ucs * auc_tau / (auc_tau + ic50_ucs)
    uDSchange <- rbase_uds - imax_uds * auc_tau / (auc_tau + ic50_uds)

    # =========================================================================
    # 5. Observations and residual error
    # =========================================================================
    # Three observables: serum concentration (ug/mL) and the two urinary-GAG
    # percent changes from baseline. dvid assignment follows the order of the
    # residual lines below, so observation rows in user event tables carry
    # dvid = 1L for Cc, 2L for uCSchange and 3L for uDSchange.
    Cc        ~ prop(propSd)
    uCSchange ~ add(addSd_uCSchange)
    uDSchange ~ add(addSd_uDSchange)
  })
}
