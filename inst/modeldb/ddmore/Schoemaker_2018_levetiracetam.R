Schoemaker_2018_levetiracetam <- function() {
  description <- "Combined adult / pediatric population PK-PD count model for levetiracetam (LEV) used in Schoemaker 2018 to scaffold a pediatric brivaracetam (BRV) extrapolation. The fitted compound is LEV; CAV is LEV plasma concentration (mg/L). Negative-binomial seizure-count likelihood with two-component mixture (responder vs non-responder), Box-Cox-transformed inter-individual variability on log baseline rate, and a Markovian dependence on the previous-period seizure count entering as a per-record covariate (PDV). The DDMORE bundle is a PD-only $PRED model: drug exposure enters as the data column CAV, so there is no PK ODE in this file. Adults contribute monthly (~28-day) counts (PDV unused, sentinel -99 in the source dataset, CHILD = 0); pediatrics contribute daily counts (NDAYS = 1, PDV = previous-day observed count, CHILD = 1)."
  reference <- paste(
    "Schoemaker R, Wade JR, Stockis A. (2018).",
    "Extrapolation of a Brivaracetam Exposure-Response Model from Adults to",
    "Children with Focal Seizures.",
    "Clin Pharmacokinet 57(7):843-854.",
    "doi:10.1007/s40262-017-0597-2.",
    "DDMORE Foundation Model Repository: DDMODEL00000239",
    "(adult+pediatric LEV PK/PD count model used in the publication to support",
    "the BRV pediatric scaling).",
    sep = " "
  )
  vignette <- "Schoemaker_2018_levetiracetam"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "mg/L"
  )
  ddmore_id <- "DDMODEL00000239"
  replicate_of <- NULL

  covariateData <- list(
    CHILD = list(
      description        = "Pediatric-vs-adult age-group indicator. 1 = pediatric (4-16 years, daily seizure-count records, NDAYS = 1, PDV = observed previous-day count, IIV active on overdispersion alpha and Markov amplitude); 0 = adult (>=18 years, monthly seizure-count records, NDAYS approx 28, PDV unused / sentinel -99, no IIV on overdispersion).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult)",
      notes              = "Source column PED (renamed to canonical CHILD). Drives four FIXED-zero peds offsets in the published model (LS00P on log baseline rate is the only non-zero peds offset, +0.420; the peds offsets on placebo, Emax, EC50, and the mixture logit are FIXED to zero). Also gates the Markov-amplitude term and the overdispersion IIV.",
      source_name        = "PED"
    ),
    TRT_PHASE = list(
      description        = "Double-blind treatment-phase indicator. 1 = record falls within the active treatment phase (placebo and Emax drug-effect terms are switched on); 0 = baseline / run-in / off-treatment (no placebo or drug effect; only the LS0 baseline-rate term contributes).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (baseline / off-treatment)",
      notes              = "Source column Q2 (renamed to canonical TRT_PHASE; new entry registered alongside this model in inst/references/covariate-columns.md). The .ctl computes LE = LS0 + Q2 * LTRTE, so TRT_PHASE multiplies the entire treatment effect (placebo + drug). New canonical because Q2 the column name collides with the canonical PK parameter q2 (inter-compartmental clearance to peripheral2).",
      source_name        = "Q2"
    ),
    PDV = list(
      description        = "Previous-period observed seizure count, supplied per-record as a covariate input. For pediatric subjects (CHILD = 1, daily counts) PDV is the observed seizure count on the immediately preceding day. For adults (CHILD = 0, monthly counts) PDV is unused; the source dataset uses the sentinel -99 to flag this.",
      units              = "(seizures per record interval)",
      type               = "count",
      reference_category = NULL,
      notes              = "Source column PDV (canonical name same as source; new entry registered alongside this model). Per operator decision (sidecar response-001 Q2 free-text answer) the Markov dependence on the previous count is preserved in this nlmixr2lib port by representing PDV as a per-record covariate the user supplies, rather than dropping the term or attempting a state-based approximation. Markov amplitude only acts when CHILD = 1 (.ctl `LS0 = LS00 + PED*LSMAX*PDV/(ES50+PDV)`); for CHILD = 0 the value is multiplied by 0 in model() so the sentinel -99 is harmless. New canonical because the name (previous-day-value seizure count) is intrinsically tied to count / Markov-feedback PD models and unlikely to be reused outside that family.",
      source_name        = "PDV"
    ),
    NDAYS = list(
      description        = "Length of the count-interval window (in days) over which the observed seizure count was tabulated. Multiplies the per-day rate to give the expected count for the interval (LAMB = exp(LE) * NDAYS in the .ctl).",
      units              = "days",
      type               = "count",
      reference_category = NULL,
      notes              = "Source column NDAYS (canonical name same as source; new entry registered alongside this model). Adult records: NDAYS = 28 (approximately 4 weeks; some bundle records use 30 or 31 reflecting the actual visit interval). Pediatric records: NDAYS = 1 (daily count). New canonical because the count-interval-length concept is specific to count / TTE PD models.",
      source_name        = "NDAYS"
    ),
    CAV = list(
      description        = "Average levetiracetam plasma concentration (mg/L) over the count-interval window. Drives the Emax drug effect via leff = lemax * CAV / (exp(lec50) + CAV). The CAV column is set to 0 for placebo records and for baseline/run-in records.",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CAV. The published model treats LEV exposure as a per-period summary; the bundle does not ship the upstream LEV popPK that produced these CAV values. For new simulations the user must supply CAV either from a LEV popPK fit they run separately or as a fixed dose-response scan. EC50 is on the LEV mg/L scale (exp(lec50) = exp(3.45) approximately 31.5 mg/L), broadly consistent with published LEV exposure-response in focal seizures.",
      source_name        = "CAV"
    )
  )

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    age_range      = "Pooled adult and pediatric (4-16 years) levetiracetam focal-seizure trial cohorts. Exact subject counts and demographic distributions are not in the DDMORE bundle and the publication PDF was not on disk for cross-check.",
    weight_range   = "(not extracted; bundle does not ship demographics)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Focal-onset seizures (uncontrolled on stable background antiepileptic therapy; LEV was added on as monotherapy or adjunctive therapy depending on the contributing trial).",
    dose_range     = "Levetiracetam plasma concentrations summarised as CAV per count-interval; the bundle does not document the dose levels the LEV cohorts received. Adult LEV add-on therapy is typically 1000-3000 mg/day in two divided doses.",
    regions        = "(not reported in the DDMORE bundle)",
    notes          = "Population field detail is intentionally sparse: the DDMORE bundle for DDMODEL00000239 ships the .ctl, the .res, and the simulated dataset (39065 rows; 6107 adult monthly-count records and 32958 pediatric daily-count records by the bundle's own PED column), but not a baseline-demographics table. The Schoemaker 2018 publication PDF was not on disk under /home/bill/github/mab_human_consensus/literature at extraction time, so subject counts, study counts, and demographic distributions are not populated. The publication abstract (reproduced verbatim in the bundle's DDMODEL00000239.rdf model-has-description block) confirms the model fit is to a combined adult + pediatric (4-16 years) LEV cohort and reports 33.5% as the mixture-responder fraction. Update population fields when the publication PDF becomes available."
  )

  ini({
    # ------------------------------------------------------------------
    # Final estimates from Output_real_P241.res FINAL PARAMETER ESTIMATE
    # block (lines 376-407 of the .res after MINIMIZATION SUCCESSFUL at
    # line 303). The .ctl $THETA / $OMEGA blocks carry the supplied
    # initial values, which differ from the .res final estimates; values
    # below are the .res finals.
    #
    # NONMEM TH index -> nlmixr2 parameter:
    #   TH 1 = lbase           (typical log baseline seizure rate, 1/day)
    #   TH 2 = es50            (Markov ES50, seizures-per-record)
    #   TH 3 = lsmax           (max Markov amplitude on log-rate scale)
    #   TH 4 = lplac           (typical log placebo effect)
    #   TH 5 = lemax            (typical log Emax -- note the .ctl applies
    #                           ETA(4) multiplicatively as TH(5)*EXP(ETA4),
    #                           so TH5 is on log-rate scale and is itself
    #                           negative; see model() block)
    #   TH 6 = lec50           (typical log EC50, mg/L LEV)
    #   TH 7 = lovdp           (typical log overdispersion alpha)
    #   TH 8 = bc_shape        (Box-Cox shape on baseline-rate eta)
    #   TH 9 = p_responder     (mixture probability of "responders")
    #   TH10 = p_responder_ped (peds offset on mixture logit; FIXED 0)
    #   TH11 = lbase_ped       (peds offset on log baseline rate)
    #   TH12 = lplac_ped       (peds offset on log placebo; FIXED 0)
    #   TH13 = lemax_ped       (peds offset on log Emax; FIXED 0)
    #   TH14 = lec50_ped       (peds offset on log EC50; FIXED 0)
    # ------------------------------------------------------------------

    lbase   <- -1.09;        label("Typical log baseline seizure rate (log(seizures/day))")             # Output_real_P241.res TH 1
    es50    <- 2.75;         label("Markov ES50 (seizures-per-record at which Markov amplitude is half-max)") # Output_real_P241.res TH 2
    lsmax   <- 1.28;         label("Maximum Markov amplitude on log-rate scale (added to log baseline rate as PDV -> Inf, only when CHILD = 1)") # Output_real_P241.res TH 3
    lplac   <- -0.16;        label("Typical log placebo effect during active treatment phase (added to log rate when TRT_PHASE = 1)") # Output_real_P241.res TH 4
    lemax   <- -3.13;        label("Typical log Emax (Emax drug effect on log rate; negative because increased exposure suppresses count)") # Output_real_P241.res TH 5
    lec50   <- 3.45;         label("Typical log EC50 (log(mg/L) of LEV; exp(lec50) approximately 31.5 mg/L)") # Output_real_P241.res TH 6
    lovdp   <- -2.24;        label("Typical log overdispersion alpha (log of unitless NB dispersion) for the negative-binomial likelihood")  # Output_real_P241.res TH 7
    bc_shape <- 0.442;       label("Box-Cox shape parameter for the IIV transform on the baseline-rate eta") # Output_real_P241.res TH 8
    p_responder <- 0.335;    label("Mixture probability of 'responder' subpopulation (1 - p_responder is the 'placebo-like' subpopulation)") # Output_real_P241.res TH 9

    # FIXED pediatric offsets -- the .ctl declares four 'Peds on ...' THETAs
    # of which only TH11 (peds offset on log baseline rate) is freely
    # estimated. The other three are FIXED 0 in the source; we keep them
    # FIXED 0 here so the structural form is preserved verbatim and a
    # downstream user can free them if they re-fit with a different cohort.
    p_responder_ped <- fixed(0); label("Peds offset on mixture-logit (FIXED 0 in source)") # Output_real_P241.res TH10 (FIXED)
    lbase_ped       <- 0.420;    label("Peds offset on log baseline seizure rate (added when CHILD = 1)") # Output_real_P241.res TH11
    lplac_ped       <- fixed(0); label("Peds offset on log placebo effect (FIXED 0 in source)")           # Output_real_P241.res TH12 (FIXED)
    lemax_ped       <- fixed(0); label("Peds offset on log Emax (FIXED 0 in source)")                     # Output_real_P241.res TH13 (FIXED)
    lec50_ped       <- fixed(0); label("Peds offset on log EC50 (FIXED 0 in source)")                     # Output_real_P241.res TH14 (FIXED)

    # ------------------------------------------------------------------
    # Inter-individual variability -- Output_real_P241.res FINAL OMEGA
    # block (lines 389-407). All entries are diagonal (OMEGA off-diagonals
    # are 0.00E+00 in the FINAL block).
    #
    # Each eta enters the model in a non-standard NONMEM way (Box-Cox,
    # multiplicative on TH(5), conditional on PED). The transformations
    # are reproduced verbatim in model() below; ini() carries only the
    # OMEGA variance values.
    # ------------------------------------------------------------------
    etalbase ~ 0.755    # ETA1 -- Box-Cox-transformed iid eta on baseline   # Output_real_P241.res OMEGA(1,1)
    etalsmax ~ 1.44     # ETA2 -- additive eta on log Markov amplitude      # Output_real_P241.res OMEGA(2,2)
    etalplac ~ 0.166    # ETA3 -- additive eta on log placebo               # Output_real_P241.res OMEGA(3,3)
    etalemax ~ 0.640    # ETA4 -- multiplicative eta on log Emax            # Output_real_P241.res OMEGA(4,4)
    etalovdp ~ 8.47     # ETA5 -- eta on log overdispersion (CHILD = 1 only) # Output_real_P241.res OMEGA(5,5)
  })

  model({
    # ------------------------------------------------------------------
    # Per-subject parameters -- translate the .ctl $PRED block.
    # ------------------------------------------------------------------

    # Box-Cox transform on the baseline-rate eta (.ctl line:
    # `TETA1 = (EXP(ETA(1))**SHP1 - 1) / SHP1`). At bc_shape = 0.442 this
    # is a moderate non-Gaussian shape on the IIV; at bc_shape -> 0 it
    # would degenerate to ETA1 itself (log-normal), at bc_shape = 1 it
    # would be exp(ETA1) - 1. nlmixr2 evaluates the expression element-
    # wise, so the same form transfers verbatim.
    bc_eta_base <- (exp(etalbase)^bc_shape - 1) / bc_shape

    # Subject-level log baseline seizure rate, with the pediatric offset.
    ls0_subject <- lbase + bc_eta_base + CHILD * lbase_ped

    # Subject-level log Markov amplitude (only used when CHILD = 1, see
    # log_rate_record below).
    lsmax_subject <- lsmax + etalsmax

    # Subject-level log placebo effect, with the pediatric offset.
    lplac_subject <- lplac + etalplac + CHILD * lplac_ped

    # Subject-level log Emax. The .ctl uses a multiplicative eta on the
    # negative typical value (`LEMAX = TVLEMAX * EXP(ETA(4)) + PED*LEMAXP`),
    # so etalemax scales the magnitude of the (negative) typical log Emax.
    # This is reproduced verbatim -- it is unusual but not a translation
    # error.
    lemax_subject <- lemax * exp(etalemax) + CHILD * lemax_ped

    # Subject-level log EC50, with the pediatric offset.
    lec50_subject <- lec50 + CHILD * lec50_ped

    # Subject-level log overdispersion. The .ctl multiplies ETA(5) by PED,
    # so adults (CHILD = 0) have no IIV on overdispersion in the source
    # (typical-value alpha = exp(-2.24) approximately 0.106 for adults).
    lovdp_subject <- lovdp + CHILD * etalovdp
    ovdp <- exp(lovdp_subject)

    # ------------------------------------------------------------------
    # Per-record terms.
    # ------------------------------------------------------------------

    # Markov amplitude term. PDV is supplied per-record as a covariate
    # (operator decision; sidecar response-001 Q2). For adult records
    # (CHILD = 0) the source dataset uses PDV = -99 as a sentinel; the
    # CHILD multiplier zeroes the term so the sentinel value is harmless.
    markov_addend <- CHILD * exp(lsmax_subject) * PDV / (es50 + PDV)
    log_rate_baseline <- ls0_subject + markov_addend

    # Drug-effect Hill term (Emax on CAV). CAV is in mg/L; exp(lec50) is
    # also in mg/L.
    leff_subject <- lemax_subject * CAV / (exp(lec50_subject) + CAV)

    # Treatment-phase contribution. The mixture-responder branch adds
    # placebo + Emax drug effect; the non-responder branch adds only
    # placebo. TRT_PHASE = 0 zeroes both branches' contributions, so
    # baseline records reduce to log_rate_baseline.
    log_rate_responder    <- log_rate_baseline + TRT_PHASE * (lplac_subject + leff_subject)
    log_rate_nonresponder <- log_rate_baseline + TRT_PHASE * lplac_subject

    # Per-record expected seizure count = exp(log rate per day) * NDAYS.
    expected_count_responder    <- exp(log_rate_responder)    * NDAYS
    expected_count_nonresponder <- exp(log_rate_nonresponder) * NDAYS

    # Mixture-responder probability -- a parameter exposed for downstream
    # use (the user computes the mixture-weighted expected count
    # externally; see vignette). Encoded as logit-additive in CHILD so
    # the FIXED-0 peds offset can be freed by a re-fitter without
    # restructuring the model.
    p_responder_subject <-
      1 / (1 + exp(-(log(p_responder / (1 - p_responder)) + CHILD * p_responder_ped)))

    # ------------------------------------------------------------------
    # Observation models -- Plan_2012_pain precedent.
    #
    # The source likelihood is a negative-binomial with overdispersion
    # alpha = exp(lovdp_subject) and per-record mean expected_count_*.
    # rxode2 / nlmixr2 do not natively express the source's joint NB
    # likelihood within the nlmixr2 ini()/model() syntax in this batch.
    # The simplification used here, matching the ddmore/Plan_2012_pain.R
    # precedent, is to declare a Poisson observation likelihood per
    # output branch with the same per-record mean. The deterministic
    # typical-value rate trajectory -- which is what the F.3
    # mechanistic-sanity vignette validates -- is unaffected; the
    # difference is only in the dispersion of stochastic VPC samples.
    # The vignette's "Assumptions and deviations" section calls this
    # out and exposes the source's overdispersion alpha as the model
    # variable `ovdp` so a downstream user can post-process Poisson
    # samples through a NB-correction step if they need it.
    # ------------------------------------------------------------------
    count_responder    <- expected_count_responder
    count_nonresponder <- expected_count_nonresponder

    count_responder    ~ pois(expected_count_responder)
    count_nonresponder ~ pois(expected_count_nonresponder)
  })
}
