Kengo_2025_ritonavir <- function() {
  description <- paste(
    "Two-compartment population PK model for oral ritonavir given as the",
    "booster of ritonavir-boosted atazanavir (ATV/r 300/100 mg) to Ugandan",
    "adults living with HIV, with and without co-administered rifampicin",
    "(DERIVE trial, Kengo 2025). Absorption is a Savic analytical transit",
    "chain (mean transit time 0.483 h, 12.3 estimated transit compartments)",
    "feeding a first-order depot (ka 1.02 1/h); disposition is",
    "two-compartment with first-order elimination. Clearance and both volumes",
    "are allometrically scaled on fat-free mass with a 42 kg reference and",
    "fixed 0.75 / 1 exponents. Rifampicin raises ritonavir clearance",
    "2.12-fold irrespective of the ATV/r dosing frequency, and cuts",
    "bioavailability by 68.8% on once-daily ATV/r; doubling ATV/r to twice",
    "daily partially restores bioavailability, leaving it 33.3% below",
    "reference. Rifampicin had no retained effect on ritonavir absorption",
    "rate. Intracellular ritonavir in peripheral blood mononuclear cells is a",
    "Sheiner-style effect compartment holding a concentration, equilibrating",
    "with plasma at an equilibration half-life of 1.40 h toward a",
    "pseudo-partition coefficient of 1.68, so ritonavir accumulates in PBMCs",
    "relative to plasma. Random effects are between-subject variability on",
    "clearance (16.4%) and eight-occasion between-occasion variability on ka",
    "(82.3%), mean transit time (43.6%) and bioavailability (55.5%); all",
    "reported percentages are the omega standard deviation on the log scale.",
    "Residual error is combined proportional plus additive, separately for",
    "plasma (25.6%, 0.001 mg/L) and PBMC (51.4%, 0.003 mg/L)."
  )
  reference <- paste(
    "Kengo A, Resendiz-Galvan JE, Najjemba L, Mugerwa H, De Nicolo A,",
    "D'Avolio A, Atoyebi S, Wiesner L, Svensson EM, Waitt C, Denti P (2025).",
    "Model-based evaluation of the interaction between ritonavir-boosted",
    "atazanavir and rifampicin in Ugandan adults with HIV.",
    "Br J Clin Pharmacol 91(12):3471-3481. doi:10.1002/bcp.70195",
    sep = " "
  )
  vignette <- "Kengo_2025_atazanavir_ritonavir_rifampicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482. The `effect` state is an exception to the usual "states hold an
  # amount" rule: the atazanavir control stream supplied with Kengo 2025 (the
  # ritonavir model shares its structure) integrates the PBMC state against a
  # concentration, DADT(4) = KE0*(PPC*C2 - A(4)) with C2 = A(2)/V, so the
  # state is a concentration in mg/L and not an amount in mg.
  compartmentData <- list(
    depot = list(
      analyte = "ritonavir", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "ritonavir", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "ritonavir", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    effect = list(
      analyte = "ritonavir", units = "mg/L",
      specimen = "blood cell", verified = TRUE
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
        "clearance (exponent 0.75) and for both volumes (exponent 1), with a",
        "fixed 42 kg reference. Kengo 2025 Table 2 footnote a states that all",
        "clearance and volume parameters for atazanavir AND ritonavir were",
        "allometrically scaled using fat-free mass and that the tabulated",
        "values refer to a typical participant of 67 kg total body weight and",
        "42 kg fat-free mass; the supplementary atazanavir control stream $PK",
        "sets TVFFM = 42 verbatim. The supplementary Table S1 and Table S3",
        "footnotes instead print 41 kg, the cohort median fat-free mass from",
        "Table 1 (41.0 kg, range 37.9-41.9); the control stream governs and",
        "42 kg is used here. See the vignette Errata.",
        "The control stream computes FFM by the Janmahasatian formula,",
        "FFM = 37.99 * HT^2 * WT / (35.98 * HT^2 + WT) for females and",
        "FFM = 42.92 * HT^2 * WT / (30.93 * HT^2 + WT) for males, with HT in",
        "m and WT in kg. Users should supply a measured or",
        "Janmahasatian-derived FFM column directly."
      ),
      source_name        = "FFM"
    ),
    CONMED_RIF = list(
      description        = "Concomitant rifampicin co-administration indicator (1 = on rifampicin, 0 = ATV/r alone)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no rifampicin; DERIVE visit 1, ATV/r 300/100 mg once daily alone)",
      notes              = paste(
        "Chronic-induction semantics: rifampicin 600 mg once daily was added",
        "to the ATV/r regimen after DERIVE visit 1 (day 7) and sampling for",
        "visit 2 was carried out on day 21, so the indicator switches on only",
        "after two weeks of daily dosing and induction is at equilibrium at",
        "every rifampicin-arm observation. Unlike the companion atazanavir",
        "model, ritonavir clearance carries a SINGLE rifampicin fold-change",
        "that applies to both the once- and twice-daily ATV/r arms: Kengo 2025",
        "Table 2 prints no 'Fold-change in CL for ATV/r BID + RIF' entry for",
        "ritonavir, and Results 3.4 reports only that rifampicin increased",
        "ritonavir clearance 2-fold. The rifampicin dose was raised from",
        "600 mg to 1200 mg once daily before visit 4, but Results 3.4 reports",
        "that increasing the rifampicin dose had no further effect on",
        "ritonavir pharmacokinetics, so no rifampicin-dose covariate enters."
      ),
      source_name        = "RIF"
    ),
    REGI_BID = list(
      description        = "ATV/r dosing-regimen indicator (1 = atazanavir/ritonavir 300/100 mg twice daily, 0 = 300/100 mg once daily)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ATV/r 300/100 mg once daily)",
      notes              = paste(
        "Modifies bioavailability only, and only in the presence of",
        "rifampicin: DERIVE never sampled twice-daily ATV/r without",
        "rifampicin, so setting REGI_BID = 1 with CONMED_RIF = 0 is",
        "extrapolation beyond the studied design and selects the unmodified",
        "reference bioavailability. With rifampicin, doubling the dosing",
        "frequency partially restores ritonavir bioavailability from 31.2% of",
        "reference to 66.7% of reference (Kengo 2025 Table 2 changes in F of",
        "-68.8% and -33.3%, and Results 3.4)."
      ),
      source_name        = "VISIT"
    ),
    OCC = list(
      description        = "Integer dosing-occasion index used for the between-occasion random effects",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Eight occasions. Kengo 2025 Methods 2.3 defines an occasion as a",
        "single administered dose, and the supplementary atazanavir control",
        "stream -- whose structure the ritonavir model shares -- multiplexes",
        "the between-occasion etas with IF (OCC==1) ... IF (OCC==8) blocks.",
        "The eight occasions span the four PK visits at two doses each, the",
        "dose taken at home before the visit and the dose given at the visit.",
        "For simulation, set OCC to the occasion index of each dosing",
        "interval; a single-occasion simulation may use OCC = 1 throughout.",
        "Unlike the companion atazanavir model, the ritonavir model carries no",
        "between-visit eta on clearance and no unobserved-dose inflation of",
        "the between-occasion etas -- Kengo 2025 Table 2 leaves both rows",
        "blank in the ritonavir column."
      ),
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    DOSE_RIF_MG = list(
      description = "Rifampicin daily dose",
      units       = "mg",
      type        = "continuous",
      notes       = paste(
        "The rifampicin dose was raised from 600 mg to 1200 mg once daily",
        "before DERIVE visit 4, but Kengo 2025 Results 3.4 reports that",
        "increasing the rifampicin dose had no further effect on ritonavir",
        "pharmacokinetics, attributed in the Discussion to near-maximal enzyme",
        "induction at the standard dose."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 26L,
    n_studies      = 1L,
    age_range      = "23-61 years (median 44)",
    age_median     = "44 years",
    weight_range   = "50-75 kg (median 67)",
    weight_median  = "67 kg",
    ffm_range      = "37.9-41.9 kg (median 41.0)",
    height_range   = "1.48-1.86 m (median 1.59)",
    bmi_range      = "19.9-31.6 kg/m2 (median 26.1)",
    sex_female_pct = 88,
    race_ethnicity = "Black African (23 of 23 with race recorded, 100%)",
    disease_state  = paste(
      "Adults living with HIV with undetectable viral load (<50 copies/mL) on",
      "ritonavir-boosted-atazanavir-based second-line antiretroviral therapy",
      "for at least 6 months. Participants with tuberculosis, hepatitis or",
      "other coinfections were excluded, as were pregnant and breastfeeding",
      "women, so rifampicin was given to participants without tuberculosis.",
      "All were taking lamivudine; 17 (65%) were also on tenofovir disoproxil",
      "fumarate, 8 (31%) on zidovudine and 1 (4%) on abacavir."
    ),
    dose_range     = paste(
      "Ritonavir 100 mg once daily as the booster of atazanavir/ritonavir",
      "300/100 mg at visit 1 (day 7); then rifampicin 600 mg once daily and",
      "dolutegravir 50 mg twice daily added with sampling at visit 2",
      "(day 21); then ATV/r increased to 300/100 mg twice daily, i.e.",
      "ritonavir 100 mg twice daily, with sampling at visit 3 (day 28); then",
      "rifampicin increased to 1200 mg once daily with sampling at visit 4",
      "(day 35)."
    ),
    regions        = "Uganda (Joint Clinical Research Centre, Kampala)",
    notes          = paste(
      "DERIVE (NCT04121195), an open-label, single-arm, dose-escalation",
      "study. Plasma sampling at predose and 0.5, 1, 2, 4, 6, 8 and 12 h",
      "postdose at every visit, with an extra 24 h sample at visit 1.",
      "857 plasma concentrations entered the pooled analysis of both drugs,",
      "of which 20 (2%) ritonavir samples were below the limit of",
      "quantification. Separate PBMC trough samples were taken at visits 1, 3",
      "and 4 and at 12 h postdose at visit 2. Lower limits of quantification",
      "were 0.005 mg/L for plasma ritonavir and 0.015 mg/L for the",
      "intracellular assay. Baseline characteristics are Table 1."
    )
  )

  ini({
    # --- Structural parameters. Final estimates are Kengo 2025 Table 2,
    # ritonavir column, with 95% confidence intervals from sampling importance
    # resampling. Values are quoted for a typical participant of 42 kg
    # fat-free mass on ATV/r 300/100 mg once daily WITHOUT rifampicin (DERIVE
    # visit 1), which is the model's reference level. Kengo 2025 Results 3.3
    # states the structure: two-compartment disposition (dOFV -271 versus one
    # compartment) with transit-compartment absorption (dOFV -52 versus a lag
    # time), the same shape as the companion atazanavir model.
    lcl <- log(9.67)
    label("Apparent oral clearance CL/F at FFM = 42 kg, ATV/r once daily without rifampicin (L/h)")  # Table 2 CL 9.67 (8.51-11.6)
    lvc <- log(55.4)
    label("Apparent central volume of distribution Vc/F at FFM = 42 kg (L)")                         # Table 2 central volume 55.4 (46.3-69.1)
    lq <- log(1.56)
    label("Apparent inter-compartmental clearance Q/F at FFM = 42 kg (L/h)")                         # Table 2 inter compartmental clearance 1.56 (1.17-2.15)
    lvp <- log(70.1)
    label("Apparent peripheral volume of distribution Vp/F at FFM = 42 kg (L)")                      # Table 2 peripheral volume 70.1 (39.7-125)
    lka <- log(1.02)
    label("First-order absorption rate constant from depot to central (1/h)")                        # Table 2 ka 1.02 (0.864-1.27); the column header prints '/L', a unit typo for the 1/h of a first-order rate constant
    lmtt <- log(0.483)
    label("Mean transit time through the absorption transit chain (h)")                              # Table 2 MTT 0.483 (0.428-0.545)
    lnn <- log(12.3)
    label("Number of absorption transit compartments in the Savic chain (unitless)")                 # Table 2 NN 12.3 (6.64-17.7)
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability at the reference regimen (fraction)")                       # Table 2 bioavailability '1 fixed'

    # --- PBMC effect-compartment parameters. Kengo 2025 reports the
    # equilibration half-life rather than the rate constant; the companion
    # atazanavir control stream parameterises and estimates the rate constant
    # directly (TVKE0 = THETA(17)), so ke0 = log(2) / t_half is the exact
    # inversion of the reported quantity: log(2) / 1.40 = 0.49511 1/h.
    lke0 <- log(log(2) / 1.40)
    label("Plasma-to-PBMC equilibration rate constant (1/h)")                                        # Table 2 PBMC equilibration half-life 1.40 h (1.38-1.63)
    lppc <- log(1.68)
    label("Plasma-to-PBMC pseudo-partition coefficient (fraction)")                                  # Table 2 PBMC pseudo-partition coefficient 1.68 (0.643-1.75); above 1, i.e. ritonavir accumulates in PBMCs

    # --- Regimen covariate effects. Ritonavir carries ONE rifampicin
    # fold-change on clearance covering both ATV/r dosing frequencies, and TWO
    # bioavailability levels. Table 2 prints no BID + RIF clearance entry and
    # no ka entry for ritonavir.
    e_rif_cl <- 2.12
    label("Fold-change in CL on rifampicin, both once- and twice-daily ATV/r (-fold)")               # Table 2 fold-change in CL for ATV/r QD + RIF 2.12 (1.95-2.31)
    e_rif_qd_fdepot <- -0.688
    label("Change in bioavailability when rifampicin is added to once-daily ATV/r (fraction)")       # Table 2 change in F for ATV/r QD + RIF -68.8% (-75.2 to -58.5); Results 3.4 'decreased its bioavailability to 31% (25-39)', and 1 - 0.688 = 0.312
    e_rif_bid_fdepot <- -0.333
    label("Change in bioavailability when rifampicin is given with twice-daily ATV/r (fraction)")    # Table 2 change in F for ATV/r BID + RIF -33.3% (-46.6 to -10.9); Results 3.4 'partially restored ritonavir bioavailability to 66% (54-89)', and 1 - 0.333 = 0.667

    # --- Allometric exponents, fixed by the authors rather than estimated.
    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent of fat-free mass on CL and Q (unitless)")                             # Kengo 2025 Methods 2.3 and Results 3.3 allometric scaling by fat-free mass
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on Vc and Vp (unitless)")                            # Kengo 2025 Methods 2.3 and Results 3.3 allometric scaling by fat-free mass

    # --- Random effects. Kengo 2025 Table 2 footnote c defines the tabulated
    # percentages as %CV = sqrt(omega^2) * 100, i.e. they are the omega
    # standard deviation on the log scale and NOT a log-normal CV. The
    # convention is confirmed against the supplementary atazanavir control
    # stream $OMEGA block, where sqrt(0.0734386) = 0.271 reproduces that
    # model's reported 27.6% BSV in clearance. Variances below are therefore
    # (percentage / 100)^2 using the Table 2 ritonavir final estimates.
    #
    # Table 2 leaves the ritonavir 'Between visit variability in clearance'
    # and 'Scaling factor on BOV for unobserved dose' rows blank, so this
    # model carries neither.
    etalcl ~ 0.026896
    label("Between-subject variability in clearance (log-scale variance)")                           # Table 2 BSV in clearance 16.4% (12.3-21.8); 0.164^2

    # Between-occasion variability on ka, MTT and F over the eight dosing
    # occasions, each a single shared variance in the source. nlmixr2 has no
    # $OMEGA BLOCK(1) SAME shortcut, so occasions 2-8 are fix()-pinned to
    # occasion 1 (the Abdelgawad_2024_linezolid pattern).
    etaiov_ka_1 ~ 0.677329
    label("Between-occasion variability in ka, occasion 1 (log-scale variance)")                     # Table 2 BOV in ka 82.3% (67.3-99.1); 0.823^2
    etaiov_ka_2 ~ fix(0.677329)                                                                      # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_3 ~ fix(0.677329)                                                                      # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_4 ~ fix(0.677329)                                                                      # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_5 ~ fix(0.677329)                                                                      # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_6 ~ fix(0.677329)                                                                      # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_7 ~ fix(0.677329)                                                                      # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_8 ~ fix(0.677329)                                                                      # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_1 ~ 0.190096
    label("Between-occasion variability in mean transit time, occasion 1 (log-scale variance)")      # Table 2 BOV in MTT 43.6% (37.7-51.6); 0.436^2
    etaiov_mtt_2 ~ fix(0.190096)                                                                     # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_3 ~ fix(0.190096)                                                                     # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_4 ~ fix(0.190096)                                                                     # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_5 ~ fix(0.190096)                                                                     # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_6 ~ fix(0.190096)                                                                     # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_7 ~ fix(0.190096)                                                                     # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_8 ~ fix(0.190096)                                                                     # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_1 ~ 0.308025
    label("Between-occasion variability in bioavailability, occasion 1 (log-scale variance)")        # Table 2 BOV in F 55.5% (48.1-63.8); 0.555^2
    etaiov_fdepot_2 ~ fix(0.308025)                                                                  # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_3 ~ fix(0.308025)                                                                  # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_4 ~ fix(0.308025)                                                                  # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_5 ~ fix(0.308025)                                                                  # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_6 ~ fix(0.308025)                                                                  # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_7 ~ fix(0.308025)                                                                  # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_8 ~ fix(0.308025)                                                                  # equal to occasion 1 per '$OMEGA BLOCK(1) SAME'

    # --- Residual error, one combined proportional-plus-additive model per
    # matrix. Kengo 2025 Methods 2.3 fixed the additive error to at least 20%
    # of the corresponding LLOQ, and the atazanavir control stream $ERROR
    # implements that as ADD = THETA + 0.2 * LLOQ with the THETA fixed to
    # zero. The tabulated ritonavir additive errors are exactly that:
    # 0.2 * 0.005 = 0.001 mg/L for plasma and 0.2 * 0.015 = 0.003 mg/L for
    # PBMC.
    propSd <- 0.256
    label("Proportional residual error for plasma ritonavir (fraction)")                             # Table 2 proportional error (plasma) 25.6% (24.1-27.8)
    addSd <- fixed(0.001)
    label("Additive residual error for plasma ritonavir (mg/L)")                                     # Table 2 additive error (plasma) '0.001 Fixed'; 20% of the 0.005 mg/L plasma LLOQ
    propSd_Cpbmc <- 0.514
    label("Proportional residual error for intracellular PBMC ritonavir (fraction)")                 # Table 2 proportional error PBMC 51.4% (42.9-61.1)
    addSd_Cpbmc <- fixed(0.003)
    label("Additive residual error for intracellular PBMC ritonavir (mg/L)")                         # Table 2 additive error PBMC '0.003 fixed'; 20% of the 0.015 mg/L intracellular LLOQ
  })

  model({
    # 1. Occasion indicators, multiplexing the between-occasion etas over the
    # eight dosing occasions.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    oc8 <- (OCC == 8)

    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 +
      oc4 * etaiov_ka_4 + oc5 * etaiov_ka_5 + oc6 * etaiov_ka_6 +
      oc7 * etaiov_ka_7 + oc8 * etaiov_ka_8
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 +
      oc4 * etaiov_mtt_4 + oc5 * etaiov_mtt_5 + oc6 * etaiov_mtt_6 +
      oc7 * etaiov_mtt_7 + oc8 * etaiov_mtt_8
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 +
      oc5 * etaiov_fdepot_5 + oc6 * etaiov_fdepot_6 +
      oc7 * etaiov_fdepot_7 + oc8 * etaiov_fdepot_8

    # 2. Regimen covariate factors. CONMED_RIF and REGI_BID together select
    # one of the three regimen levels the source fitted; with CONMED_RIF = 0
    # every factor collapses to 1 and the reference regimen is recovered.
    rif_qd <- CONMED_RIF * (1 - REGI_BID)
    rif_bid <- CONMED_RIF * REGI_BID

    # One clearance fold-change covering both rifampicin arms.
    cl_rif <- 1 + (e_rif_cl - 1) * CONMED_RIF
    # Two bioavailability levels; exactly one branch is active.
    fdepot_rif <- 1 + e_rif_qd_fdepot * rif_qd + e_rif_bid_fdepot * rif_bid

    # 3. Individual parameters. Clearance and both volumes are allometrically
    # scaled on fat-free mass against the 42 kg reference; ka, MTT and NN are
    # not size-scaled.
    cl <- exp(lcl + etalcl) * (FFM / 42)^e_ffm_cl * cl_rif
    vc <- exp(lvc) * (FFM / 42)^e_ffm_vc
    q <- exp(lq) * (FFM / 42)^e_ffm_cl
    vp <- exp(lvp) * (FFM / 42)^e_ffm_vc
    ka <- exp(lka + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    nn <- exp(lnn)
    fdepot <- exp(lfdepot + iov_fdepot) * fdepot_rif
    ke0 <- exp(lke0)
    ppc <- exp(lppc)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. Savic transit-compartment absorption, implemented by
    # rxode2's transit() closed form; the source control stream sets F1 = 0 so
    # the dose does not also arrive in the depot as a bolus, which
    # f(depot) <- 0 preserves.
    Cc <- central / vc

    d/dt(depot) <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # The PBMC state holds a concentration, not an amount.
    d/dt(effect) <- ke0 * (ppc * Cc - effect)

    f(depot) <- 0

    # 5. Observations and residual error.
    Cpbmc <- effect
    Cc ~ add(addSd) + prop(propSd)
    Cpbmc ~ add(addSd_Cpbmc) + prop(propSd_Cpbmc)
  })
}
