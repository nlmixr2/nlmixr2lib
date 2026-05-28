Khan_2015_ciprofloxacin <- function() {
  description <- "Mechanism-based six-compartment PK/PD model for ciprofloxacin against Escherichia coli K-12 wild-type and quinolone-resistant single-step mutants in static in vitro time-kill experiments (Khan et al. 2015 J Antimicrob Chemother). Two co-existing bacterial subpopulations (susceptible-growing, plus a small pre-existing less-susceptible subpopulation seeded at MUT*1e-6 of the inoculum) each carry growing, resting, and drug-induced non-colony-forming states; drug effect is a Hill-Emax term on bacterial kill with strain-specific EC50 (LM202 estimated; LM347, LM378, LM534, LM625, LM693, LM707 fixed at the published per-strain values), a separate EC502 for the resistant subpopulation, density-dependent active->resting flux, and a model-time gate that switches off the active->non-colony-forming transition after THETA(19) hours."
  reference <- paste(
    "Khan DD, Lagerback P, Cao S, Lustig U, Nielsen EI, Cars O, Hughes D, Andersson DI, Friberg LE. (2015).",
    "A mechanism-based pharmacokinetic/pharmacodynamic model allows prediction of antibiotic killing from",
    "MIC values for WT and mutants. J Antimicrob Chemother 70(11):3051-3060.",
    "doi:10.1093/jac/dkv233.",
    "DDMORE Foundation Model Repository: DDMODEL00000225.",
    sep = " "
  )
  vignette <- "Khan_2015_ciprofloxacin"
  units <- list(time = "hour", dosing = "mg/L", concentration = "log CFU/mL")
  ddmore_id    <- "DDMODEL00000225"
  replicate_of <- NULL

  covariateData <- list(
    STR = list(
      description        = "Bacterial strain identifier (E. coli K-12 derivative). Values 347, 202, 378, 534, 625, 693, 707 select the matching strain-specific EC50 in `model()`. Strain LM202 is the wild-type reference (estimated EC50); the other six are quinolone-resistant single-step mutants with fixed published EC50s.",
      units              = "(integer code)",
      type               = "categorical",
      reference_category = "202",
      notes              = "Per-tube fixed covariate. The DDMORE bundle's cipro_simulated.csv ships only STR=202; cipro202.csv and cipro378.csv carry the real-data fits for those two strains. Any STR value outside {347, 202, 378, 534, 625, 693, 707} drives the strain cascade in `model()` to ec50 = 0, which is intentionally unsafe (drug effect would diverge); guard the input data accordingly. In vitro experimental-design column, not a population-PK covariate; not added to the canonical inst/references/covariate-columns.md register.",
      source_name        = "STR"
    ),
    CAB = list(
      description        = "Static (constant-during-tube) ciprofloxacin concentration in the in vitro time-kill experiment. Drives the Hill-Emax DRUGS / DRUGS2 kill terms and the active <-> non-colony-forming transition rates.",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-tube fixed covariate. CAB = 0 represents the no-drug control arm (drug effect is identically 0). The original Khan 2015 experiment static-dosed CAB at 0, 0.0625, 0.125, 0.25, 1, 2, 4, 8, and 16 x MIC for each strain; the bundle's cipro_simulated.csv replays a similar grid for STR = 202 (CAB in {0, 0.0029, 0.0059, 0.0118, 0.0235, 0.047, 0.094, 0.188, 0.376} mg/L). In vitro experimental-design column, not a population-PK covariate.",
      source_name        = "CAB"
    ),
    BASE = list(
      description        = "Per-tube logarithm (natural-ln) of the baseline bacterial inoculum (log CFU/mL). Sets sbase = exp(BASE), which seeds the susceptible-growing compartment at sbase * (1 - mut*1e-6) and the pre-existing resistant-growing compartment at sbase * mut * 1e-6.",
      units              = "log CFU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-tube fixed covariate; the .mod's $PK B2-method baseline-IIV term `exp(ETA(1)*sqrt(SIGMA(1)))` is dropped here (see vignette Errata). The DDMORE bundle's cipro_simulated.csv carries BASE values around 12.9-13 (ln of ~4-5e5 CFU/mL initial inoculum). In vitro experimental-design column, not a population-PK covariate.",
      source_name        = "BASE"
    )
  )

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "Static in vitro time-kill experiments on Escherichia coli K-12 MG1655 wild-type and six gyrA/gyrB/parC/parE single-step mutants (LM347, LM378, LM534, LM625, LM693, LM707) raised on the same K-12 MG1655 background. Bacteria were grown in Mueller-Hinton broth and exposed to a static ciprofloxacin concentration; viable counts were obtained by serial-dilution colony counting at 0, 1, 2, 4, 6, 9, 12, and 24 hours. NOT a human population-PK study.",
    dose_range     = "Static (constant-during-tube) ciprofloxacin concentrations: 0 (no-drug control) plus a per-strain grid of MIC multipliers (0.0625x, 0.125x, 0.25x, 1x, 2x, 4x, 8x, 16x of each strain's measured MIC). LM202 wild-type MIC = 0.047 mg/L. The DDMORE bundle's simulated dataset uses the same grid for STR = 202 only.",
    regions        = NA_character_,
    notes          = "DDMODEL00000225's bundle does NOT include the original Khan 2015 publication, the published Table 1 (per-strain MICs and EC50s) or Table 2 (model parameter estimates) -- the doi:10.1093/jac/dkv233 article is paywalled and is not on disk in this worktree. Population-style demographic fields (n_subjects, age, weight, sex) do not apply to an in vitro bacterial culture experiment and are intentionally NA. The simulated dataset cipro_simulated.csv ships 9 tubes (one per CAB level for STR = 202) with 4 sample-position replicates per observation time, totalling 207 rows. See the validation vignette's Assumptions and deviations / Errata section for the full list of bundle-versus-publication caveats."
  )

  ini({
    # Parameter values come from the .mod $THETA initial estimates in
    # executeable_cipro_pkpd.mod (lines 167-185). The bundle does not ship
    # an Output_real_*.lst; output_simulated.lst is a refit on a simulated
    # dataset and round-trip-recovers the .mod initials within ~3% for
    # THETA but with SIGMA(1) materially different (2.42 -> 1.78). Per
    # operator decision (sidecar 025 q1 = mod_initials), the .mod $THETA
    # initials are treated as the publication-derived point values that
    # seeded the bundle's simulated dataset; addSd uses sqrt(SIGMA(1) =
    # 2.41597) to match.
    #
    # Six of the seven strain-specific EC50s are FIX in the .mod
    # ($THETA(3), (5)..(9)), reflecting the publication's reported
    # per-mutant MIC values; only the wild-type LM202 EC50 ($THETA(4))
    # was estimated. The Hill exponent on the active->NC transition
    # ($THETA(17) = 20) is also FIX. We mirror those FIX flags here.

    # Sensitive subpopulation: growth and drug-effect Emax / EC50
    lkgs       <- log(1.70)              ; label("Growth rate of susceptible (sensitive) subpopulation kgs (1/h)")                           # .mod $THETA(1) = (1, 1.70); init 1.70
    lemax      <- log(5.24)              ; label("Maximum drug-induced bacterial-kill rate Emax (1/h)")                                       # .mod $THETA(2) = (0, 5.24); init 5.24
    lec50_347  <- fixed(log(0.037))      ; label("EC50 for E. coli LM347 (mg/L)")                                                              # .mod $THETA(3) = (0, 0.037) FIX
    lec50_202  <- log(0.057)             ; label("EC50 for E. coli LM202 wild-type (mg/L)")                                                    # .mod $THETA(4) = (0, 0.057); init 0.057 (only LM202 was estimated; ciprotest.csv carries LM202 only)
    lec50_378  <- fixed(log(0.65))       ; label("EC50 for E. coli LM378 (mg/L)")                                                              # .mod $THETA(5) = (0, 0.65) FIX
    lec50_534  <- fixed(log(0.30))       ; label("EC50 for E. coli LM534 (mg/L)")                                                              # .mod $THETA(6) = (0, 0.30) FIX
    lec50_625  <- fixed(log(1.0))        ; label("EC50 for E. coli LM625 (mg/L)")                                                              # .mod $THETA(7) = (0, 1.0) FIX
    lec50_693  <- fixed(log(31.0))       ; label("EC50 for E. coli LM693 (mg/L)")                                                              # .mod $THETA(8) = (0, 31) FIX
    lec50_707  <- fixed(log(92.0))       ; label("EC50 for E. coli LM707 (mg/L)")                                                              # .mod $THETA(9) = (0, 92) FIX
    lgam       <- log(1.98)              ; label("Hill coefficient for drug-effect Emax curve (unitless)")                                     # .mod $THETA(10) = (0.5, 1.98); init 1.98

    # Density-dependent active -> resting transition (proportionality
    # constant). PC is the .mod's THETA(11) raw value; the .mod multiplies
    # by 1e-7 inside $PK before use.
    lpc_raw    <- log(0.0186)            ; label("Raw active->resting proportionality constant PC (multiplied by 1e-7 in model() to give 1/CFU)") # .mod $THETA(11) = (0, 0.0186); init 0.0186

    # Pre-existing resistant subpopulation (subpop 2): independent growth
    # rate kgs2 and EC50 ec502; seeded at MUT * 1e-6 of the baseline.
    lkgs2      <- log(0.344)             ; label("Growth rate of pre-existing resistant subpopulation kgs2 (1/h)")                              # .mod $THETA(12) = (0.18, 0.344); init 0.344
    lec502     <- log(1.25)              ; label("EC50 of pre-existing resistant subpopulation EC502 (mg/L)")                                   # .mod $THETA(13) = (0, 1.25); init 1.25
    lmut_raw   <- log(0.81)              ; label("Raw pre-existing resistant fraction MUT (multiplied by 1e-6 in model() to give a unitless seed fraction)") # .mod $THETA(14) = (0, 0.81); init 0.81

    # Active <-> non-colony-forming (NC) transition (drug-driven).
    # ksnc_max gates active -> NC; knc_factor gates NC -> active back to
    # the active pool when drug pressure relaxes.
    lksnc_max  <- log(5.83)              ; label("Maximum rate of active -> non-colony-forming transition kSNc,max (1/h)")                      # .mod $THETA(15) = (0, 5.83); init 5.83
    lknc_factor <- log(0.17)             ; label("NC -> active transition rate factor sfNcS (unitless; multiplies EC50/CAB to give 1/h)")        # .mod $THETA(16) = (0, 0.17); init 0.17
    lgam_nc    <- fixed(log(20))         ; label("Hill exponent for active -> NC transition (unitless)")                                        # .mod $THETA(17) = (0, 20) FIX
    ltr50      <- log(0.24)              ; label("Fold-EC50 ratio at which active -> NC transition is 50% saturated (unitless)")                # .mod $THETA(18) = (0, 0.24); init 0.24

    # Model-time gate: active -> NC transition switches off after tmtime
    # hours (NONMEM's MTIME(2) mechanism translated to a step on `t`).
    ltmtime    <- log(5.3158)            ; label("Time after which the active -> NC transition shuts off MTIME(2) (h)")                         # .mod $THETA(19) = (2, 5.3158); init 5.3158

    # Across-tube residual error on log-bacterial-count scale. addSd =
    # sqrt(SIGMA(1)). The .mod additionally has four BLOCK(1) SAME SIGMA
    # slots gated on FLG2 (sample position 1..4); per operator decision
    # (sidecar 025 q3 = across_tube_only) those per-position random
    # effects are documented in the vignette Errata and not encoded here.
    addSd <- 1.5544                      ; label("Additive residual SD on log-bacterial-count scale (log CFU/mL); = sqrt(.mod $SIGMA(1) = 2.41597)") # .mod $SIGMA  2.41597
  })
  model({
    # Bacterial death rate is a hardcoded constant in the .mod $PK
    # (line 32: "KK = 0.179"); not a $THETA, not a fitted parameter.
    kk <- 0.179

    # Back-transformed structural parameters.
    kgs        <- exp(lkgs)
    kgs2       <- exp(lkgs2)
    emax       <- exp(lemax)
    gam        <- exp(lgam)
    ec502      <- exp(lec502)
    pc         <- exp(lpc_raw)   * 1e-7
    mut        <- exp(lmut_raw)
    ksnc_max   <- exp(lksnc_max)
    knc_factor <- exp(lknc_factor)
    gam_nc     <- exp(lgam_nc)
    tr50       <- exp(ltr50)
    tmtime     <- exp(ltmtime)

    # Strain-specific EC50 cascade (mg/L). The seven STR codes are
    # mutually exclusive; (STR == X) evaluates to 1 for the matching
    # strain and 0 otherwise, so the sum picks exactly one EC50 per
    # tube. Mirrors the .mod $PK lines 39-45 IF(STR.EQ.X) cascade.
    ec50 <- exp(lec50_347) * (STR == 347) +
            exp(lec50_202) * (STR == 202) +
            exp(lec50_378) * (STR == 378) +
            exp(lec50_534) * (STR == 534) +
            exp(lec50_625) * (STR == 625) +
            exp(lec50_693) * (STR == 693) +
            exp(lec50_707) * (STR == 707)

    # Hill-Emax drug effect on bacterial kill, separately for each
    # subpopulation. The .mod guards the division with IF(CAB > 1e-11);
    # since gam > 0 always (lower bound 0.5 in the .mod), CAB^gam is
    # exactly 0 at CAB = 0 and the +1e-30 below avoids 0/0 only at the
    # double-zero limit.
    drugs  <- emax * CAB^gam / (ec50^gam  + CAB^gam + 1e-30)
    drugs2 <- emax * CAB^gam / (ec502^gam + CAB^gam + 1e-30)

    # Active -> non-colony-forming (NC) transition rates, drug-driven
    # via a Hill curve in fold-EC50 with exponent gam_nc and half-max
    # at fold-ratio tr50.
    ksnc_1 <- ksnc_max * (CAB / ec50 )^gam_nc /
              ((CAB / ec50 )^gam_nc + tr50^gam_nc + 1e-30)
    ksnc_2 <- ksnc_max * (CAB / ec502)^gam_nc /
              ((CAB / ec502)^gam_nc + tr50^gam_nc + 1e-30)

    # NC -> active transition rates: fast when drug is gone, asymptote
    # to 0 when CAB >> EC50. The .mod adds 1e-10 to the denominator to
    # avoid division by zero when CAB == 0; reproduced verbatim.
    knc_s_1 <- knc_factor * ec50  / (CAB + 1e-10)
    knc_s_2 <- knc_factor * ec502 / (CAB + 1e-10)

    # MTIME(2) gate: 1 while t < tmtime, 0 thereafter. Closes off
    # the active -> NC transition once the experiment has run past
    # tmtime hours. The .mod implements this via NONMEM's
    # MPAST(1) - MPAST(2) construct against MTIME(1)=0 / MTIME(2)=tmtime.
    flag <- 1.0 * (t < tmtime)

    # Density-dependent active -> resting flux. The .mod uses
    # PC * (A1+A2+A3+A4+A5+A6); SR2 = SR (same rate for both
    # subpopulations); RS / RS2 = 0 (no resting -> active flux).
    total_bact <- bact_s + bact_r + bact_spe + bact_np + bact_rpe + bact_nppe
    sr  <- pc * total_bact
    sr2 <- sr

    # Per-tube baseline bacterial count (log scale -> linear via exp).
    # The .mod's B2 method seeds bact_s and bact_spe with sbase scaled
    # by the pre-existing resistant fraction; the .mod's
    # exp(ETA(1)*sqrt(SIGMA(1))) IIV term on the baseline is dropped
    # here (see vignette Errata).
    sbase <- exp(BASE)
    bact_s(0)   <- sbase * (1 - mut * 1e-6)
    bact_spe(0) <- sbase *      mut * 1e-6
    # bact_r, bact_np, bact_rpe, bact_nppe default to 0.

    # ODE system; line-for-line image of the .mod $DES block (lines
    # 117-123). NC compartments share the same drug-induced kill
    # term (kk + drugs) as the active subpopulation they came from.
    d/dt(bact_s)    <- kgs  * bact_s   - (kk + drugs)  * bact_s   - sr  * bact_s   + knc_s_1 * bact_np   - ksnc_1 * bact_s   * flag
    d/dt(bact_r)    <- -kk  * bact_r   + sr  * bact_s
    d/dt(bact_spe)  <- kgs2 * bact_spe - (kk + drugs2) * bact_spe - sr2 * bact_spe + knc_s_2 * bact_nppe - ksnc_2 * bact_spe * flag
    d/dt(bact_np)   <- ksnc_1 * bact_s   * flag - knc_s_1 * bact_np   - (kk + drugs)  * bact_np
    d/dt(bact_rpe)  <- -kk  * bact_rpe + sr2 * bact_spe
    d/dt(bact_nppe) <- ksnc_2 * bact_spe * flag - knc_s_2 * bact_nppe - (kk + drugs2) * bact_nppe

    # Observation: log of the colony-forming bacterial count. The .mod's
    # ATOT = A1 + A2 + A3 + A5 explicitly excludes A4 (NP) and A6 (NPPE)
    # because those non-colony-forming states do not contribute to a
    # plate count; reproduced verbatim. The +1e-8 floor matches the .mod
    # $ERROR line 132: IPRED = LOG(ATOT + 1E-8).
    counted_bact <- bact_s + bact_r + bact_spe + bact_rpe
    log_cfu <- log(counted_bact + 1e-8)
    log_cfu ~ add(addSd)
  })
}
