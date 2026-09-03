Dodds_2013_psoriasis_biologics_mbma <- function() {
  description <- paste0(
    "MBMA. Dose-response meta-analysis of monoclonal-antibody efficacy in ",
    "moderate-to-severe plaque psoriasis, used as the truth generator for a ",
    "clinical trial simulation comparing concentrated and distributed ",
    "first-in-patient trial designs. The meta-analysis pooled 27 randomized ",
    "controlled trials of anti-tumour-necrosis-factor agents, ustekinumab ",
    "and methotrexate, jointly modelling the probability of a PASI 50, 75 ",
    "and 90 response together with the mean percent improvement in Psoriasis ",
    "Area Severity Index (PASI); the percent-improvement arm of that ",
    "meta-analysis is what this model encodes. The log ratio of on-treatment ",
    "to baseline PASI is a placebo response plus an Emax function of dose, ",
    "with one Emax shared by all three marketed antibodies and an ED50 ",
    "estimated per antibody (adalimumab 16.9 mg, golimumab 45.5 mg, ",
    "ustekinumab 13.9 mg). Dose enters through one CONMED_<drug>_DOSE column ",
    "per antibody; there are no rxode2 dose events, no PK layer and no time ",
    "axis, because the model predicts a single per-subject PASI response at ",
    "the end of a 12-week window after one dose. Unlike a per-arm MBMA this ",
    "model DOES describe individual subjects: the source adds a subject-level ",
    "residual whose standard deviation is linear in the predicted response, ",
    "calibrated to the between-subject spread of the PHOENIX 1 and 2 ",
    "ustekinumab trials. Outputs are the log PASI ratio, the typical percent ",
    "PASI improvement, and the absolute percentage-point difference from ",
    "placebo that the source's Go/No-Go criterion is applied to. The five ",
    "hypothetical compounds of the source's Table 1 are deliberately NOT ",
    "included; see the validation vignette."
  )

  reference <- paste(
    "Dodds MG, Salinger DH, Mandema J, Gibbs JP, Gibbs MA.",
    "Clinical trial simulation to inform phase 2: comparison of concentrated",
    "vs. distributed first-in-patient study designs in psoriasis.",
    "CPT Pharmacometrics Syst Pharmacol. 2013;2(7):e58.",
    "doi:10.1038/psp.2013.32.",
    "The model equations are Eqs. 1-3 of the Methods section 'Clinical trial",
    "simulation'; the parameter values are in Table 1 and its footnote.",
    sep = " "
  )

  vignette <- "Dodds_2013_psoriasis_biologics_mbma"

  # Algebraic dose-response model: no rxode2 dose events are consumed (the dose
  # enters through the CONMED_<drug>_DOSE covariate columns), there is no time
  # axis, and the outputs are dimensionless PASI ratios and percentages rather
  # than drug concentrations. The `units` entries below follow the placeholder
  # convention already used by Mandema_2011_biologicDMARDs_mbma and
  # Vargo_2014_statins_ezetimibe_mbma so that checkModelConventions() sees a
  # parseable dosing / concentration pair.
  units <- list(
    time          = "week (placeholder; the model is a time-independent dose-response evaluated once per subject at the end of the source's 12-week post-dose observation window, not a time course)",
    dosing        = "mg/administration (single-dose milligrams supplied through the CONMED_<drug>_DOSE covariate columns, NOT as rxode2 dose events; the source's simulated dose levels are 0, 21, 70, 210 and 700 mg)",
    concentration = "%/subject (percent improvement in PASI score for one subject; the lpasi output is a dimensionless log ratio. Output is NOT a drug concentration; the slash satisfies checkModelConventions parsing)"
  )

  covariateData <- list(
    CONMED_ADALIMUMAB_DOSE = list(
      description        = "Single-dose adalimumab dose assigned to the subject in the simulated first-in-patient trial.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an adalimumab arm; 0 in every arm makes this a placebo subject. IMPORTANT UNITS WARNING: this is a SINGLE-DOSE milligram amount read against an ED50 of 16.9 mg, which is NOT the same dose metric as the identically-named column in Mandema_2011_biologicDMARDs_mbma (mg per administration on a once-every-2-weeks maintenance regimen, ED50 17.8 mg q2w). The two models must not share a dose column value. The source's simulated dose levels are 0, 21, 70, 210 and 700 mg, chosen as the deliverable volumes of a 70 mg/mL presentation (0.3, 1, 3 and 10 mL); 700 mg is the maximum feasible dose assumed tolerated in healthy volunteers.",
      source_name        = "Dose (Dodds 2013 Eq. 1; Table 1 adalimumab ED50 row; Table 2 dose-level assignments)"
    ),
    CONMED_GOLIMUMAB_DOSE = list(
      description        = "Single-dose golimumab dose assigned to the subject in the simulated first-in-patient trial.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a golimumab arm. Same single-dose milligram metric as CONMED_ADALIMUMAB_DOSE above, read against an ED50 of 45.5 mg; see that entry's units warning about the Mandema_2011_biologicDMARDs_mbma column of the same name. Golimumab is the least potent of the three marketed antibodies in this analysis, which is why the source's Table 4 shows its ED50 is the one most often recovered within twofold by a 16-subject trial.",
      source_name        = "Dose (Dodds 2013 Eq. 1; Table 1 golimumab ED50 row; Table 2 dose-level assignments)"
    ),
    CONMED_USTEKINUMAB_DOSE = list(
      description        = "Single-dose ustekinumab dose assigned to the subject in the simulated first-in-patient trial.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a ustekinumab arm. Same single-dose milligram metric as CONMED_ADALIMUMAB_DOSE above, read against an ED50 of 13.9 mg. Ustekinumab is the most potent of the three marketed antibodies here, and it is also the compound the source used to calibrate the subject-level variance model: Figure 3 compares simulated mean and standard deviation of percent PASI improvement at 45 and 90 mg against the observed PHOENIX 2 data. Those two doses are the licensed ustekinumab doses and are inside, but at the low end of, the 0-700 mg range this model was exercised over.",
      source_name        = "Dose (Dodds 2013 Eq. 1; Table 1 ustekinumab ED50 row; Figure 3 PHOENIX 2 comparison)"
    )
  )

  covariatesDataExcluded <- list(
    SCORE_PASI_BASE = list(
      description        = "Baseline PASI score before treatment. Screened as a meta-analysis covariate but NOT retained in the final dose-response model.",
      units              = "(score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Dodds 2013 Methods 'MBMA of biologic treatments for psoriasis' lists baseline PASI score among the characteristics evaluated for impact on the model parameters; no effect was retained and no point estimate is reported. Baseline PASI is nevertheless structurally present in this model as the DENOMINATOR of the response: the modelled quantity is log(PASI_treatment / PASI_base), so a subject's own baseline cancels out of every prediction and only the ratio is predicted. A user who wants absolute on-treatment PASI must multiply exp(lpasi) by their own baseline score.",
      source_name        = "baseline PASI score (Dodds 2013 Methods covariate list)"
    ),
    AGE = list(
      description        = "Patient age. Screened as a meta-analysis covariate but NOT retained in the final dose-response model.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Listed in Dodds 2013 Methods among 'endpoint, drug, drug class, regimen, indication (psoriasis vs. psoriatic rheumatoid arthritis), failure of prior treatment, baseline PASI score, disease duration, age, weight, and gender' evaluated for impact on model parameters. No effect retained and no point estimate reported.",
      source_name        = "age (Dodds 2013 Methods covariate list)"
    ),
    WT = list(
      description        = "Patient body weight. Screened as a meta-analysis covariate but NOT retained in the final dose-response model.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Listed in the Dodds 2013 Methods covariate screen; no effect retained and no point estimate reported. Note that all three antibodies enter this model as flat milligram doses, not mg/kg, so body weight does not enter through the dose units either.",
      source_name        = "weight (Dodds 2013 Methods covariate list)"
    ),
    SEXF = list(
      description        = "Sex. Screened as a meta-analysis covariate but NOT retained in the final dose-response model.",
      units              = "(binary)",
      type               = "binary",
      reference_category = NULL,
      notes              = "Listed as 'gender' in the Dodds 2013 Methods covariate screen; no effect retained and no point estimate reported.",
      source_name        = "gender (Dodds 2013 Methods covariate list)"
    ),
    DIS_DURATION_PSORIASIS = list(
      description        = "Duration of psoriasis at trial entry. Screened as a meta-analysis covariate but NOT retained in the final dose-response model.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Listed as 'disease duration' in the Dodds 2013 Methods covariate screen; no effect retained and no point estimate reported.",
      source_name        = "disease duration (Dodds 2013 Methods covariate list)"
    ),
    PRIOR_TX_FAILURE = list(
      description        = "Failure of prior treatment. Screened as a meta-analysis covariate but NOT retained in the final dose-response model.",
      units              = "(binary)",
      type               = "binary",
      reference_category = NULL,
      notes              = "Listed as 'failure of prior treatment' in the Dodds 2013 Methods covariate screen; no effect retained and no point estimate reported.",
      source_name        = "failure of prior treatment (Dodds 2013 Methods covariate list)"
    )
  )

  population <- list(
    species        = "human",
    n_studies      = 27L,
    disease_state  = "Moderate-to-severe plaque psoriasis. The meta-analysis pooled randomized controlled trials in psoriasis AND in psoriatic rheumatoid arthritis, and tested indication as a covariate; the dose-response encoded here is the psoriasis percent-PASI-improvement arm. The endpoint is percent improvement in Psoriasis Area Severity Index over baseline, assessed in a 12-week window after a single dose.",
    dose_range     = "Marketed-antibody doses as studied in the pooled trials; the clinical trial simulation exercised the model at single doses of 0, 21, 70, 210 and 700 mg, corresponding to 0, 0.3, 1, 3 and 10 mL of an assumed 70 mg/mL presentation. 700 mg intravenous is the maximum dose assumed to have been shown safe and well tolerated in healthy volunteers.",
    notes          = "Two DIFFERENT populations are involved and should not be confused. (i) The dose-response parameters come from a model-based meta-analysis of 27 published randomized controlled trials of anti-tumour-necrosis-factor agents, ustekinumab and methotrexate, identified by PubMed search plus review of previous meta-analyses, clinicaltrials.gov, conference abstracts and corporate websites; the source reports no subject counts, no demographic table and no per-trial listing for this database, so n_subjects and the demographic fields are unavailable. (ii) The subject-level variance model (Eqs. 2-3) was calibrated separately, to reproduce the between-subject variability of percent PASI improvement observed in the PHOENIX 1 and PHOENIX 2 ustekinumab trials. The simulated first-in-patient trials themselves were 16 subjects each, 9,999 replicates per design per compound, run in NONMEM 7.2. The Go/No-Go criterion the source applies to this model is an absolute improvement over placebo of more than 50 percentage points at the maximum feasible dose of 700 mg; the dose-response criterion is recovery of ED50 within twofold. Simulation scope: this model predicts ONE PASI observation per subject, with no time course and no PK; it is not a repeated-measures or longitudinal model."
  )

  ini({
    # ========================================================================
    # PROVENANCE FOR THE WHOLE ini() BLOCK
    #
    # Dodds 2013 reports its dose-response on the REAL scale (Table 1 gives the
    # maximal absolute difference from placebo in PASI percentage points, and
    # the table footnote gives the placebo response) while the model itself is
    # written in the LOG domain (Eq. 1: "the parameters E0, Emax, and ED50
    # reflect the response in the log domain"). The two log-domain values below
    # are therefore obtained by inverting the paper's own definitions, not by
    # reading a printed log-scale number:
    #
    #   exp(e0)                  = 1 - 0.095            (placebo response)
    #   exp(e0) - exp(e0 + emax) = 0.823                (maximal difference)
    #
    # The inversion is confirmed by three independent printed numbers that were
    # NOT used to derive it: the absolute difference from placebo at the maximum
    # feasible 700 mg dose, which the source states is 81.8% for adalimumab,
    # 81.0% for golimumab and 81.9% for ustekinumab. The parameters below
    # reproduce 81.82, 81.01 and 81.91. The vignette shows the arithmetic.
    #
    # Every parameter is fixed(): none was estimated here, and re-estimating
    # them against subject-level data would be a category error - they are
    # meta-analytic summaries of published trial results.
    # ========================================================================

    e0 <- fixed(-0.0998203353)
    label("Placebo response, i.e. log ratio of on-treatment to baseline PASI at zero dose (log ratio)")
    # Dodds 2013 Table 1 footnote: "Placebo response (E0) is 9.5% for all
    # compounds." A 9.5% improvement is a PASI ratio of 1 - 0.095 = 0.905, so
    # e0 = log(0.905) = -0.0998203353. The source calls this response "small but
    # non-negligible" and holds it common to all compounds.

    emax <- fixed(-2.4012156964)
    label("Maximal drug effect on the log PASI ratio, shared by all three marketed antibodies (log ratio)")
    # Dodds 2013 Table 1: "All of the marketed compounds were found to have
    # equal maximum change in PASI % improvement (absolute maximal difference
    # from placebo = 82.3%)." Inverting on the real scale,
    # exp(e0) - exp(e0 + emax) = 0.823, so emax = log(1 - 0.823/0.905)
    # = log(0.0906077) = -2.4012156964. NEGATIVE because an improvement in PASI
    # is a REDUCTION in the score and therefore a reduction in the ratio.

    # ------------------------------------------------------------------------
    # ED50 by drug, log-transformed. "The marketed compounds were found to have
    # different doses that produce half of the maximum effect (ED50 = 16.9,
    # 45.5, and 13.9 mg for adalimumab, golimumab, and ustekinumab,
    # respectively)." Single-dose milligrams - see the covariateData units
    # warning. Each value is independently confirmed by the twofold acceptance
    # window printed in the same table row.
    # ------------------------------------------------------------------------
    led50_adalimumab <- fixed(log(16.9))
    label("Log adalimumab ED50 (mg, single dose); back-transform ED50 = 16.9 mg")
    # Dodds 2013 Table 1 adalimumab row. The same row's "estimated ED50 within"
    # window of 8.50-33.8 mg is twofold either side of 16.9 (8.45-33.8), which
    # confirms the transcription.

    led50_golimumab <- fixed(log(45.5))
    label("Log golimumab ED50 (mg, single dose); back-transform ED50 = 45.5 mg")
    # Dodds 2013 Table 1 golimumab row; twofold window 22.8-91.0 mg confirms
    # 45.5 (22.75-91.0).

    led50_ustekinumab <- fixed(log(13.9))
    label("Log ustekinumab ED50 (mg, single dose); back-transform ED50 = 13.9 mg")
    # Dodds 2013 Table 1 ustekinumab row; twofold window 6.90-27.7 mg confirms
    # 13.9 (6.95-27.8, printed rounded outward).

    # ========================================================================
    # SUBJECT-LEVEL RESIDUAL VARIABILITY (Eqs. 2-3)
    #
    #   Y     = LPASI + (sERR*LPASI + iERR) * eps,  eps ~ N(0, 1)
    #   sERR  = -0.32,  iERR = 0.30
    #
    # so the residual standard deviation is LINEAR in the predicted response
    # rather than constant, and it GROWS as the response gets larger: the
    # source notes "sERR is given with a negative sign, as LPASI is negative".
    #
    # ENCODING. rxode2's combined1() variance model computes
    #   SD = addSd + propSd * prediction
    # on the SIGNED prediction, which is exactly the shape of Eq. 2 - but
    # rxode2 requires propSd > 0, and Eq. 2's slope is negative. The model
    # block therefore predicts pasiLogDrop = -LPASI, the log PASI REDUCTION, which
    # is positive over the whole dose range, and the slope becomes +0.32:
    #
    #   SD = 0.30 + 0.32 * (-LPASI) = 0.30 - 0.32 * LPASI = iERR + sERR * LPASI
    #
    # Because eps is symmetric about zero, negating the endpoint leaves the
    # distribution of the response unchanged; the encoded model is Eq. 2
    # exactly, not an approximation. The vignette checks the simulated standard
    # deviation against iERR + sERR*LPASI dose by dose.
    #
    # This is a BETWEEN-SUBJECT variance, not the between-study variance that
    # an MBMA usually reports: the source calibrated it against the observed
    # subject-level spread in PHOENIX 1 and 2, and simulated one PASI
    # observation per subject. There is no separate eta and no between-trial
    # random effect - the source's trials differ only through which subjects
    # they happen to draw.
    # ========================================================================
    addSd <- fixed(0.30)
    label("Intercept of the subject-level residual standard deviation, the source's iERR (log ratio)")
    # Dodds 2013 Eq. 3: "iERR = 0.30".

    propSd <- fixed(0.32)
    label("Slope of the subject-level residual standard deviation on the log PASI reduction, the magnitude of the source's sERR (unitless)")
    # Dodds 2013 Eq. 3: "sERR = -0.32". The sign is carried by the negated
    # endpoint pasiLogDrop = -lpasi; see the encoding note above.
  })

  model({
    # ---- Per-drug ED50, back-transformed -----------------------------------
    ed50_adalimumab   <- exp(led50_adalimumab)
    ed50_golimumab    <- exp(led50_golimumab)
    ed50_ustekinumab  <- exp(led50_ustekinumab)

    # ---- Saturating dose fraction Dose/(Dose + ED50), per drug -------------
    # Every simulated subject received at most ONE antibody, so at most one of
    # these fractions is non-zero. A placebo subject sets all three
    # CONMED_<drug>_DOSE columns to 0, every fraction collapses to 0, and the
    # response reduces to the placebo response e0. Setting two columns
    # non-zero simultaneously is outside the source's calibration and would
    # make the model ADD the two drug effects.
    f_adalimumab   <- CONMED_ADALIMUMAB_DOSE  / (CONMED_ADALIMUMAB_DOSE  + ed50_adalimumab)
    f_golimumab    <- CONMED_GOLIMUMAB_DOSE   / (CONMED_GOLIMUMAB_DOSE   + ed50_golimumab)
    f_ustekinumab  <- CONMED_USTEKINUMAB_DOSE / (CONMED_USTEKINUMAB_DOSE + ed50_ustekinumab)

    # ---- Drug effect on the log PASI ratio ---------------------------------
    # One Emax shared by all three antibodies; the drugs differ ONLY in
    # potency. This is the source's central structural finding for the
    # marketed compounds and is what makes ED50 the quantity a first-in-patient
    # trial has to recover.
    edrug <- emax * (f_adalimumab + f_golimumab + f_ustekinumab)

    # ---- Eq. 1: the log PASI ratio -----------------------------------------
    # LPASI = log(PASI_treatment / PASI_base) = E0 + Emax*Dose/(ED50 + Dose).
    # Negative throughout: PASI falls under treatment. The log transform is
    # what bounds percent improvement to [-Inf%, 100%] on the real scale.
    lpasi <- e0 + edrug

    # ---- Reporting quantities on the real scale ----------------------------
    # Percent improvement from baseline INCLUDING the placebo response, which
    # is what the source's Figures 1 and 3 plot.
    pasiPctImprove <- 100 * (1 - exp(lpasi))

    # Absolute difference from the placebo response, in PASI percentage
    # points. This is the quantity the Go/No-Go criterion is applied to: the
    # source requires more than 50 percentage points at the maximum feasible
    # 700 mg dose for a compound to advance.
    pasiPctImproveDrug <- 100 * (exp(e0) - exp(lpasi))

    # ---- Observation and error (Eqs. 2-3) ----------------------------------
    # The endpoint is the log PASI REDUCTION -lpasi, which is positive, so
    # that rxode2's combined1() variance model reproduces Eq. 2's negative
    # slope with a positive propSd. See the ini() encoding note.
    pasiLogDrop <- -lpasi
    pasiLogDrop ~ add(addSd) + prop(propSd) + combined1()
  })
}
