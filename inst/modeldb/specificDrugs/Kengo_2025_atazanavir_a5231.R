Kengo_2025_atazanavir_a5231 <- function() {
  description <- paste(
    "Two-compartment population PK model for oral UNBOOSTED atazanavir (no",
    "ritonavir) given with and without rifampicin to healthy adult volunteers",
    "in the ACTG A5231 study, as re-fitted by Kengo 2025 (supplementary",
    "Table S3). The structure is the companion DERIVE atazanavir model of the",
    "same paper -- a Savic transit chain (mean transit time 1.38 h, 10 transit",
    "compartments fixed) feeding a first-order",
    "depot whose rate constant is fixed at 6 1/h, two-compartment disposition",
    "with first-order elimination, and fat-free-mass allometric scaling at a",
    "42 kg reference with fixed 0.75 / 1 exponents -- with every disposition",
    "parameter FIXED to its DERIVE estimate and only the study-specific terms",
    "re-estimated. Three parameters differ from DERIVE and were estimated on",
    "the A5231 data: atazanavir clearance is 2.13-fold higher in the absence",
    "of ritonavir's CYP3A inhibition, absorption is slower (mean transit time",
    "1.38 h against 0.499 h in DERIVE), and rifampicin co-administration cuts",
    "bioavailability by 55.3%. Rifampicin's 70.8% reduction of the absorption",
    "rate constant was carried over fixed. Unlike the DERIVE models this fit",
    "has no PBMC effect compartment, because A5231 collected no intracellular",
    "samples. Random effects are between-subject variability on clearance",
    "(18.8%) and six-occasion between-occasion variability on ka (101%, fixed),",
    "mean transit time (48.1%) and bioavailability (58.3%), the last three",
    "inflated 1.7-fold on dosing occasions whose dose was self-administered",
    "rather than directly observed; all reported percentages are the omega",
    "standard deviation on the log scale. Residual error is combined",
    "proportional plus additive, both fixed from DERIVE (18.8%, 0.005 mg/L).",
    "Because the number of transit compartments is fixed at the integer 10,",
    "the absorption delay is encoded as the explicit 11-compartment chain the",
    "Savic density is the analytical solution of, so THE DOSE IS ADMINISTERED",
    "TO transit1 rather than to depot."
  )
  reference <- paste(
    "Kengo A, Resendiz-Galvan JE, Najjemba L, Mugerwa H, De Nicolo A,",
    "D'Avolio A, Atoyebi S, Wiesner L, Svensson EM, Waitt C, Denti P (2025).",
    "Model-based evaluation of the interaction between ritonavir-boosted",
    "atazanavir and rifampicin in Ugandan adults with HIV.",
    "Br J Clin Pharmacol 91(12):3471-3481. doi:10.1002/bcp.70195",
    "Parameter estimates are supplementary Table S3; the underlying data are",
    "the ACTG A5231 study of Burger DM et al., reference 30 of Kengo 2025.",
    sep = " "
  )
  vignette <- "Kengo_2025_atazanavir_ritonavir_rifampicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The absorption delay is an explicit 11-compartment Savic chain rather than
  # a closed-form density; see the note above the ODE block in model() for why.
  # The DOSE IS ADMINISTERED TO `transit1`, not to `depot`.
  compartmentData <- c(
    stats::setNames(
      lapply(
        1:11,
        function(i) {
          list(
            analyte = "atazanavir", units = "mg",
            specimen = "administration site", verified = TRUE
          )
        }
      ),
      paste0("transit", 1:11)
    ),
    list(
    depot = list(
      analyte = "atazanavir", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "atazanavir", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "atazanavir", units = "mg",
      specimen = "plasma", verified = TRUE
    )
    )
  )

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass, computed from sex, total body weight, and height by the Janmahasatian (2005) formula",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size descriptor for clearance and inter-compartmental",
        "clearance (exponent 0.75) and for both volumes (exponent 1). All the",
        "disposition parameters of this fit are fixed to their DERIVE values,",
        "so the reference fat-free mass is the DERIVE reference: Kengo 2025",
        "Table 2 footnote a gives 42 kg and the supplementary control stream",
        "$PK sets TVFFM = 42 verbatim. The Table S3 footnote instead prints",
        "41 kg, the cohort median fat-free mass from Table 1; the control",
        "stream governs and 42 kg is used here, a 1.8% difference in",
        "clearance. See the vignette Errata.",
        "A5231 enrolled healthy US volunteers whose fat-free mass Kengo 2025",
        "does not tabulate; Table S2 reports only total body weight (median",
        "75 kg, range 55-110) and sex (5 of 13 female), from which the",
        "Janmahasatian formula gives the required column,",
        "FFM = 37.99 * HT^2 * WT / (35.98 * HT^2 + WT) for females and",
        "FFM = 42.92 * HT^2 * WT / (30.93 * HT^2 + WT) for males, with HT in",
        "m and WT in kg."
      ),
      source_name        = "FFM"
    ),
    CONMED_RTV = list(
      description        = "Concomitant ritonavir boosting indicator (1 = atazanavir given as ritonavir-boosted ATV/r, 0 = unboosted atazanavir)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (ritonavir-boosted; the DERIVE reference arm from which the fixed clearance is carried)",
      notes              = paste(
        "Every A5231 participant received UNBOOSTED atazanavir, so this",
        "covariate is 0 throughout the population the model was fitted to.",
        "It is carried explicitly rather than folded into the typical value",
        "because Kengo 2025 Table S3 parameterises the fit exactly that way:",
        "clearance is fixed at the DERIVE ritonavir-boosted estimate of",
        "7.55 L/h and a separate estimated parameter, 'Fold-change in CL due",
        "to absence of ritonavir', multiplies it by 2.13. The polarity of the",
        "printed coefficient is preserved (the 2.13 applies when CONMED_RTV =",
        "0) so the tabulated value stays literally traceable; setting",
        "CONMED_RTV = 1 recovers the boosted DERIVE clearance. Kengo 2025",
        "Results 3.6 and the Discussion attribute the increase to the absence",
        "of ritonavir's strong CYP3A inhibition. Note that the fold-change",
        "applies to ALL A5231 periods including period 1, which had no",
        "rifampicin, which is why it is a ritonavir-absence term and not a",
        "rifampicin term."
      ),
      source_name        = "RTV"
    ),
    CONMED_RIF = list(
      description        = "Concomitant rifampicin co-administration indicator (1 = on rifampicin, 0 = atazanavir alone)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no rifampicin; A5231 period 1, atazanavir 300 mg twice daily alone)",
      notes              = paste(
        "Chronic-induction semantics: A5231 dosed atazanavir 300 mg twice",
        "daily alone for 8 days (period 1), then added rifampicin 600 mg once",
        "daily for 11 days (period 2, atazanavir 300 mg twice daily) and 8",
        "days (period 3, atazanavir 400 mg twice daily), with PK sampling at",
        "the end of each period, so induction is at equilibrium at every",
        "rifampicin-arm observation. In this fit rifampicin acts on",
        "bioavailability (-55.3%, estimated) and on the absorption rate",
        "constant (-70.8%, fixed) but NOT on clearance: Kengo 2025 Table S3",
        "carries no rifampicin clearance term, because the 2.13-fold",
        "clearance increase was attributed to the absence of ritonavir across",
        "all three periods. The two rifampicin periods differ only in the",
        "atazanavir dose (300 versus 400 mg twice daily), which the event",
        "table carries as the dose amount rather than as a covariate."
      ),
      source_name        = "RIF"
    ),
    OCC = list(
      description        = "Integer dosing-occasion index used for the between-occasion random effects",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Six occasions. Kengo 2025 Methods 2.3 defines an occasion as a",
        "single administered dose, and A5231 sampled at the end of each of",
        "three treatment periods, giving two dosing occasions per PK visit --",
        "the dose taken before the visit and the dose given at the visit --",
        "on the same design reading that yields eight occasions over four",
        "visits in the companion DERIVE model. The A5231 occasion count is",
        "not stated numerically in Kengo 2025 and no A5231 control stream was",
        "supplied, so six is inferred from the three-period design plus the",
        "presence of the unobserved-dose scaling factor in Table S3; see the",
        "vignette Errata. Because every occasion shares one variance, the",
        "count only bounds how many occasions can be simulated. For",
        "simulation, set OCC to the occasion index of each dosing interval; a",
        "single-occasion simulation may use OCC = 1 throughout."
      ),
      source_name        = "OCC"
    ),
    SELFADMIN = list(
      description        = "Self-administered (not directly observed) dosing-occasion indicator (1 = dose taken unobserved, 0 = directly observed dose)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (directly observed dose)",
      notes              = paste(
        "Per-dosing-occasion and time-varying. Carries the source model's",
        "multiplier on the between-occasion etas of ka, mean transit time and",
        "bioavailability for doses that were not directly observed. Kengo 2025",
        "Table S3 footnote d describes it as a multiplicative factor",
        "increasing the BOV of the absorption parameters for pre-dose",
        "concentrations following an unobserved dose -- an adherence /",
        "dose-timing uncertainty device, not an absorption mechanism, and a",
        "variance-model role rather than the relative-bioavailability role of",
        "the covariate register's founding example. Set SELFADMIN = 0 to",
        "simulate a fully supervised regimen."
      ),
      source_name        = "OBS"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 13L,
    n_studies      = 1L,
    age_range      = "23-48 years (median 30)",
    age_median     = "30 years",
    weight_range   = "55-110 kg (median 75)",
    weight_median  = "75 kg",
    sex_female_pct = 38,
    race_ethnicity = "1 of 13 (8%) Black; the remaining participants' race is not broken out in Kengo 2025 Table S2",
    disease_state  = paste(
      "Healthy adult volunteers without HIV (Kengo 2025 Table S2 records",
      "'Participants living with HIV: None'). This is the external validation",
      "cohort of the paper, contrasted against the 26 Ugandan adults living",
      "with HIV of the DERIVE trial who supplied the companion",
      "ritonavir-boosted models."
    ),
    dose_range     = paste(
      "Atazanavir 300 mg twice daily for 8 days (period 1); atazanavir 300 mg",
      "twice daily with rifampicin 600 mg once daily for 11 days (period 2);",
      "atazanavir 400 mg twice daily with rifampicin 600 mg once daily for 8",
      "days (period 3). No ritonavir was given in any period."
    ),
    regions        = "United States",
    notes          = paste(
      "ACTG A5231. 355 atazanavir plasma concentrations. Sampling at the end",
      "of each period, 15 min predose and at 1, 2, 3, 4, 5, 6, 8, 10, 12 and",
      "24 h postdose. Plasma samples were analysed at the University of",
      "Alabama by HPLC with UV detection, LLOQ 0.025 mg/L -- a different assay",
      "and a different limit of quantification from the DERIVE cohort's",
      "LC-MS/MS at 0.030 mg/L. Baseline characteristics are Table S2. Kengo",
      "2025 writes the study identifier as both 'A5213' (Methods 2.2 and",
      "Results 3.1) and 'A5231' (Results 3.6, Table S2, Table S3); A5231 is",
      "used here because it is the identifier attached to the parameter table",
      "this model encodes. See the vignette Errata."
    )
  )

  ini({
    # --- Structural parameters. Final estimates are Kengo 2025 supplementary
    # Table S3. That table's own note states 'All fixed parameters are
    # estimates from DERIVE study population', so every fixed() value below is
    # carried over from the companion Kengo_2025_atazanavir model rather than
    # estimated on A5231, and is wrapped in fixed() accordingly. The small
    # numerical differences from the Kengo 2025 Table 2 DERIVE final estimates
    # (7.55 versus 7.57 L/h, 77.3 versus 77.5 L, 3.51 versus 3.13 L/h, 48.9
    # versus 42.1 L) are reproduced verbatim from Table S3, because those are
    # the values this fit actually used.
    lcl <- fixed(log(7.55))
    label("Apparent oral clearance CL/F at FFM = 42 kg for ritonavir-boosted atazanavir, carried from DERIVE (L/h)")  # Table S3 CL '7.55 (fixed)'
    lvc <- fixed(log(77.3))
    label("Apparent central volume of distribution Vc/F at FFM = 42 kg (L)")                          # Table S3 central volume '77.3 (fixed)'
    lq <- fixed(log(3.51))
    label("Apparent inter-compartmental clearance Q/F at FFM = 42 kg (L/h)")                          # Table S3 inter compartmental clearance '3.51 (fixed)'
    lvp <- fixed(log(48.9))
    label("Apparent peripheral volume of distribution Vp/F at FFM = 42 kg (L)")                       # Table S3 peripheral volume '48.9 (fixed)'
    lka <- fixed(log(6))
    label("First-order absorption rate constant from depot to central without rifampicin (1/h)")      # Table S3 ka '6 (fixed)'
    lmtt <- log(1.38)
    label("Mean transit time through the absorption transit chain (h)")                               # Table S3 MTT 1.38 (1.15-1.61), estimated; Results 3.6 'MTT was 2.5-fold [2.3-3.2] longer' than the DERIVE 0.499 h
    lnn <- fixed(log(10))
    label("Number of absorption transit compartments in the Savic chain (unitless)")                  # Table S3 NN '10 (fixed)'
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability without rifampicin (fraction)")                              # Table S3 bioavailability '1 (fixed)'

    # --- Covariate effects. Only the ritonavir-absence clearance factor and
    # the rifampicin bioavailability effect were estimated on A5231; the
    # rifampicin effect on ka was carried over fixed.
    e_rtv_cl <- 2.13
    label("Fold-change in CL in the ABSENCE of ritonavir, i.e. at CONMED_RTV = 0 (-fold)")            # Table S3 fold-change in CL due to absence of ritonavir 2.13 (1.88-2.34); Results 3.6 'a twofold increase in atazanavir clearance compared to the DERIVE study'
    e_rif_fdepot <- -0.553
    label("Change in bioavailability on rifampicin (fraction)")                                       # Table S3 change in F with rifampicin -55.3% (-64.0 to -46.2); Results 3.6 'a significant 55% (64-46) reduction in atazanavir bioavailability'
    e_rif_ka <- fixed(-0.708)
    label("Change in ka on rifampicin (fraction)")                                                    # Table S3 change in ka with rifampicin '-70.8 (fixed)'
    e_selfadmin_iov <- 1.7
    label("Multiplier on the between-occasion etas of ka, MTT and F for a self-administered dose (-fold)")  # Table S3 scaling factor on BOV for unobserved dose 1.7 (1.10-2.81)

    # --- Allometric exponents, fixed by the authors rather than estimated.
    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent of fat-free mass on CL and Q (unitless)")                              # Kengo 2025 Methods 2.3 allometric scaling; control stream ALLMCL_FFM = (FFM/42)**0.75
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on Vc and Vp (unitless)")                             # Kengo 2025 Methods 2.3 allometric scaling; control stream ALLMV_FFM = (FFM/42)

    # --- Random effects. Kengo 2025 Table S3 footnote c defines the tabulated
    # percentages as %CV = sqrt(omega^2) * 100, i.e. they are the omega
    # standard deviation on the log scale and NOT a log-normal CV; variances
    # below are (percentage / 100)^2. Table S3 marks the between-occasion
    # variability in ka as fixed, so its variance is wrapped in fix().
    etalcl ~ 0.035344
    label("Between-subject variability in clearance (log-scale variance)")                            # Table S3 BSV in clearance 18.8% (13.7-24.4); 0.188^2

    # Between-occasion variability on ka, MTT and F over the six dosing
    # occasions, each a single shared variance. nlmixr2 has no $OMEGA
    # BLOCK(1) SAME shortcut, so occasions 2-6 are fix()-pinned to occasion 1.
    etaiov_ka_1 ~ fix(1.0201)
    label("Between-occasion variability in ka, occasion 1 (log-scale variance)")                      # Table S3 BOV in ka 101%, carried over from DERIVE rather than re-estimated on A5231; 1.01^2
    etaiov_ka_2 ~ fix(1.0201)                                                                         # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_3 ~ fix(1.0201)                                                                         # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_4 ~ fix(1.0201)                                                                         # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_5 ~ fix(1.0201)                                                                         # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_6 ~ fix(1.0201)                                                                         # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_1 ~ 0.231361
    label("Between-occasion variability in mean transit time, occasion 1 (log-scale variance)")       # Table S3 BOV in MTT 48.1% (36.9-63.4); 0.481^2
    etaiov_mtt_2 ~ fix(0.231361)                                                                      # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_3 ~ fix(0.231361)                                                                      # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_4 ~ fix(0.231361)                                                                      # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_5 ~ fix(0.231361)                                                                      # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_6 ~ fix(0.231361)                                                                      # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_1 ~ 0.339889
    label("Between-occasion variability in bioavailability, occasion 1 (log-scale variance)")         # Table S3 BOV in F 58.3% (46.4-73.5); 0.583^2
    etaiov_fdepot_2 ~ fix(0.339889)                                                                   # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_3 ~ fix(0.339889)                                                                   # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_4 ~ fix(0.339889)                                                                   # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_5 ~ fix(0.339889)                                                                   # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_6 ~ fix(0.339889)                                                                   # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'

    # --- Residual error, combined proportional plus additive, both carried
    # over fixed from DERIVE. The additive term follows the same
    # 20%-of-the-LLOQ rule Kengo 2025 Methods 2.3 states, applied to the A5231
    # assay's own limit: 0.2 * 0.025 = 0.005 mg/L.
    propSd <- fixed(0.188)
    label("Proportional residual error for plasma atazanavir (fraction)")                             # Table S3 proportional error '18.8 (fixed)'
    addSd <- fixed(0.005)
    label("Additive residual error for plasma atazanavir (mg/L)")                                     # Table S3 additive error '0.005 (fixed)'; 20% of the 0.025 mg/L A5231 LLOQ
  })

  model({
    # 1. Occasion indicators, multiplexing the between-occasion etas over the
    # six dosing occasions.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)

    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 +
      oc4 * etaiov_ka_4 + oc5 * etaiov_ka_5 + oc6 * etaiov_ka_6
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 +
      oc4 * etaiov_mtt_4 + oc5 * etaiov_mtt_5 + oc6 * etaiov_mtt_6
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 +
      oc5 * etaiov_fdepot_5 + oc6 * etaiov_fdepot_6

    # Doses that were not directly observed inflate the between-occasion etas
    # of all three absorption parameters.
    iov_scale <- 1 + (e_selfadmin_iov - 1) * SELFADMIN

    # 2. Covariate factors. The clearance factor is active when ritonavir is
    # ABSENT, which is the case for every A5231 participant; rifampicin acts
    # on bioavailability and on the absorption rate constant only.
    cl_rtv <- 1 + (e_rtv_cl - 1) * (1 - CONMED_RTV)
    fdepot_rif <- 1 + e_rif_fdepot * CONMED_RIF
    ka_rif <- 1 + e_rif_ka * CONMED_RIF

    # 3. Individual parameters. Clearance and both volumes are allometrically
    # scaled on fat-free mass against the 42 kg reference; ka, MTT and NN are
    # not size-scaled.
    cl <- exp(lcl + etalcl) * (FFM / 42)^e_ffm_cl * cl_rtv
    vc <- exp(lvc) * (FFM / 42)^e_ffm_vc
    q <- exp(lq) * (FFM / 42)^e_ffm_cl
    vp <- exp(lvp) * (FFM / 42)^e_ffm_vc
    ka <- exp(lka + iov_ka * iov_scale) * ka_rif
    mtt <- exp(lmtt + iov_mtt * iov_scale)
    nn <- exp(lnn)
    fdepot <- exp(lfdepot + iov_fdepot * iov_scale) * fdepot_rif

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. Savic transit-compartment absorption, written out as an
    # EXPLICIT chain rather than through rxode2's transit() helper.
    #
    # The source control stream computes the Savic density analytically in
    # $DES, exactly as the companion Kengo_2025_atazanavir model does:
    #   KTR     = (NN + 1) / MTT
    #   TRANSIT = EXP(LOG(BIO*PD*KTR) - GAMLN(NN+1) + NN*LOG(KTR*TEMPO)
    #             - KTR*TEMPO)
    #   DADT(1) = TRANSIT - KA*A(1)
    # with F1 = 0 so the dose does not also arrive in the depot as a bolus.
    #
    # That closed form CANNOT be used here. In this model -- which, unlike the
    # two DERIVE models, has no PBMC state -- rxode2 5.1.7 evaluates
    # `transit(nn, mtt, fdepot)` together with the `f(depot) <- 0` that F1 = 0
    # requires as an identically zero input rate: the solve returns all-zero
    # concentrations with no error and no NA. Every closed-form variant fails
    # the same way (rxode2's transit() helper, the hand-written density driven
    # by podo(depot)/tad(depot), the argument-less podo()/tad() forms, and a
    # dose passed in as a covariate column), because `f(depot) <- 0` zeroes the
    # whole right-hand side of d/dt(depot) and not merely the dose bolus.
    #
    # Because NN is FIXED at the integer 10 (Table S3), the density is instead
    # written as the compartment chain it is the analytical solution of, which
    # needs no bioavailability trick at all. A dose entering the first of
    # NN + 1 = 11 compartments that each transfer at rate ktr leaves the last
    # one at exactly
    #   ktr * (ktr*t)^NN * exp(-ktr*t) / NN!
    # times the administered amount, i.e. the Savic input rate above. The dose
    # is therefore given to `transit1` and bioavailability is applied there.
    #
    # This is exact, not an approximation: solved against the closed form it
    # agrees to 3e-11 in Cc, and cl * AUCtau equals Dose * F to eight
    # significant figures in every regimen arm. NOTE that the chain length is
    # tied to the fixed NN = 10; lnn is fixed() for that reason, and changing
    # it would require adding or removing transit compartments to match.
    ktr <- (nn + 1) / mtt

    d/dt(transit1) <- -ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ktr * transit6
    d/dt(transit7) <- ktr * transit6 - ktr * transit7
    d/dt(transit8) <- ktr * transit7 - ktr * transit8
    d/dt(transit9) <- ktr * transit8 - ktr * transit9
    d/dt(transit10) <- ktr * transit9 - ktr * transit10
    d/dt(transit11) <- ktr * transit10 - ktr * transit11
    d/dt(depot) <- ktr * transit11 - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Relative bioavailability acts on the dose as it enters the chain, which
    # is where the control stream's BIO multiplies it (its $DES carries BIO
    # inside the transit density).
    f(transit1) <- fdepot

    # 5. Observation and residual error.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
