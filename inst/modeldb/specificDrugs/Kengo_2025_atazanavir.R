Kengo_2025_atazanavir <- function() {
  description <- paste(
    "Two-compartment population PK model for oral atazanavir given as",
    "ritonavir-boosted atazanavir (ATV/r 300/100 mg) to Ugandan adults living",
    "with HIV, with and without co-administered rifampicin (DERIVE trial,",
    "Kengo 2025). Absorption is a Savic analytical transit chain (mean transit",
    "time 0.499 h, 10 transit compartments fixed) feeding a first-order depot",
    "whose rate constant is fixed at 6 1/h; disposition is two-compartment with",
    "first-order elimination. Clearance and both volumes are allometrically",
    "scaled on fat-free mass with a 42 kg reference and fixed 0.75 / 1",
    "exponents. Rifampicin co-administration is carried as a three-level",
    "regimen categorical: relative to ATV/r once daily without rifampicin,",
    "adding rifampicin to once-daily ATV/r raises clearance 3.05-fold and cuts",
    "bioavailability by 52.5%, whereas doubling ATV/r to twice daily with",
    "rifampicin raises clearance only 2.03-fold and fully restores",
    "bioavailability; rifampicin slows absorption by 67.3% on ka in both",
    "rifampicin arms. Intracellular atazanavir in peripheral blood mononuclear",
    "cells is a Sheiner-style effect compartment holding a concentration,",
    "equilibrating with plasma at an equilibration half-life of 0.963 h toward",
    "a pseudo-partition coefficient of 0.653, neither of which was affected by",
    "rifampicin. Random effects are between-subject variability on clearance",
    "(27.6%), between-visit variability on clearance across the four study",
    "visits (17.5%), and eight-occasion between-occasion variability on ka",
    "(97.9%), mean transit time (59.4%) and bioavailability (48.2%), the last",
    "three inflated 1.63-fold on dosing occasions whose dose was",
    "self-administered rather than directly observed; all reported percentages",
    "are the omega standard deviation on the log scale. Residual error is",
    "combined proportional plus additive, separately for plasma (19.8%,",
    "0.006 mg/L) and PBMC (74.9%, 0.003 mg/L)."
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
  # amount" rule: the supplementary control stream's $DES integrates the PBMC
  # state against a concentration, DADT(4) = KE0*(PPC*C2 - A(4)) with
  # C2 = A(2)/V, so the state is a concentration in mg/L and not an amount in
  # mg. Its $ERROR block reads it back directly as CC = A(4).
  compartmentData <- list(
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
    ),
    effect = list(
      analyte = "atazanavir", units = "mg/L",
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
        "clearance and volume parameters were allometrically scaled using",
        "fat-free mass and that the tabulated values refer to a typical",
        "participant of 67 kg total body weight and 42 kg fat-free mass; the",
        "supplementary control stream $PK sets TVFFM = 42 verbatim. Note that",
        "the supplementary Table S1 and Table S3 footnotes instead print",
        "41 kg, which is the cohort median fat-free mass reported in Table 1",
        "(41.0 kg, range 37.9-41.9); the control stream governs and 42 kg is",
        "used here. The difference is a 1.8% shift in clearance",
        "((42/41)^0.75) and is documented in the vignette Errata.",
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
        "after two weeks of daily dosing and the CYP3A / P-glycoprotein",
        "induction is at equilibrium at every rifampicin-arm observation.",
        "Together with REGI_BID this reproduces the source model's three-level",
        "regimen categorical: the control stream's CL_VISIT selects THETA(12)",
        "at visit 2 and THETA(13) at visits 3 and 4, and its BIO_VISIT selects",
        "THETA(14) at visit 2 only. The rifampicin dose was raised from 600 mg",
        "to 1200 mg once daily before visit 4, but Kengo 2025 Results 3.4",
        "reports no significant effect of the higher rifampicin dose on",
        "atazanavir clearance or bioavailability, so this model carries no",
        "rifampicin-dose covariate and visits 3 and 4 share one level."
      ),
      source_name        = "RIF"
    ),
    REGI_BID = list(
      description        = "ATV/r dosing-regimen indicator (1 = atazanavir/ritonavir 300/100 mg twice daily, 0 = 300/100 mg once daily)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ATV/r 300/100 mg once daily)",
      notes              = paste(
        "Interacts with CONMED_RIF only: the model has no BID level in the",
        "absence of rifampicin because DERIVE never sampled twice-daily ATV/r",
        "without rifampicin. Setting REGI_BID = 1 with CONMED_RIF = 0 is",
        "therefore extrapolation beyond the studied design and selects the",
        "unmodified reference clearance and bioavailability. With rifampicin,",
        "doubling the dosing frequency lowers the rifampicin induction of",
        "clearance from 3.05-fold to 2.03-fold and fully restores",
        "bioavailability (Kengo 2025 Table 2 and Results 3.4)."
      ),
      source_name        = "VISIT"
    ),
    OCC = list(
      description        = "Integer dosing-occasion index used for the between-occasion and between-visit random effects",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Eight occasions, matching the supplementary control stream's",
        "IF (OCC==1) ... IF (OCC==8) multiplexers over ETA(18)-ETA(25) for",
        "bioavailability, ETA(26)-ETA(33) for ka and ETA(34)-ETA(41) for mean",
        "transit time. Kengo 2025 Methods 2.3 defines an occasion as a single",
        "administered dose. The eight occasions span the four PK visits at two",
        "doses each -- the dose self-administered at home before the visit and",
        "the dose given at the visit -- which is the only reading consistent",
        "with the stream's four VISIT levels, its eight OCC levels, and its",
        "SCALE_BOV device for pre-dose concentrations that follow an",
        "unobserved dose. The model derives the visit level for the",
        "between-visit clearance eta from OCC on that map (occasions 1-2 =",
        "visit 1, 3-4 = visit 2, 5-6 = visit 3, 7-8 = visit 4), so no separate",
        "visit column is needed. For simulation, set OCC to the occasion index",
        "of each dosing interval; a single-occasion simulation may use",
        "OCC = 1 throughout."
      ),
      source_name        = "OCC"
    ),
    SELFADMIN = list(
      description        = "Self-administered (not directly observed) dosing-occasion indicator (1 = dose taken unobserved, 0 = directly observed dose)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (directly observed dose)",
      notes              = paste(
        "Per-dosing-occasion and time-varying. Carries the control stream's",
        "IF (OBS.EQ.0) branch, which multiplies the between-occasion etas of",
        "ka, mean transit time and bioavailability by SCALE_BOV = THETA(10)",
        "for doses that were not directly observed. Kengo 2025 Table 2",
        "footnote d describes it as a multiplicative factor increasing the BOV",
        "of the absorption parameters for pre-dose concentrations following an",
        "unobserved dose -- an adherence / dose-timing uncertainty device, not",
        "an absorption mechanism. Note that this is a variance-model role",
        "rather than the relative-bioavailability role that the covariate",
        "register's founding example (Wallender 2021) uses; the sign",
        "convention is the same (1 = unobserved). Set SELFADMIN = 0 to",
        "simulate a fully supervised regimen."
      ),
      source_name        = "OBS"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_RTV_AUC = list(
      description = "Ritonavir AUC over the 0-24 h dosing interval",
      units       = "mg*h/L",
      type        = "continuous",
      notes       = paste(
        "Screened as a continuous covariate on atazanavir clearance",
        "(control stream $PK 'RTVAUC_CL = 1 + THETA(11)*(AUC_RTV2 -",
        "MEDIAN_AUC24_RTV)' centred on a median ritonavir AUC24 of",
        "7.0923 mg*h/L) but NOT retained: THETA(11) is '(-0.0616, 0, 1502)",
        "FIX', i.e. fixed at zero in the final run. Kengo 2025 Results 3.4",
        "reports that attempts to fit a ritonavir-based inhibition model for",
        "atazanavir led to limited improvement in fit and reduced parameter",
        "precision, because atazanavir and ritonavir share a clearance pathway",
        "and were given at a fixed dose ratio, making causality and",
        "correlation inseparable in this dataset."
      )
    ),
    DOSE_RIF_MG = list(
      description = "Rifampicin daily dose",
      units       = "mg",
      type        = "continuous",
      notes       = paste(
        "The rifampicin dose was raised from 600 mg to 1200 mg once daily",
        "before DERIVE visit 4, but Kengo 2025 Results 3.4 reports no",
        "significant effect of increasing the rifampicin dose on atazanavir",
        "clearance or bioavailability, attributed in the Discussion to",
        "near-maximal enzyme induction at the standard dose. Visits 3 and 4",
        "therefore share one regimen level and no rifampicin-dose term enters",
        "the model."
      )
    ),
    CONMED_TDF = list(
      description = "Concomitant tenofovir disoproxil fumarate indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on both atazanavir clearance (dOFV = 0.64, P = 0.424) and",
        "bioavailability (dOFV = 2.66, P = 0.103) and retained on neither",
        "(Kengo 2025 Results 3.4). 17 of 26 DERIVE participants (65%) were on",
        "TDF, 8 (31%) on zidovudine and 1 (4%) on abacavir."
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
      "Atazanavir/ritonavir 300/100 mg once daily at visit 1 (day 7); then",
      "rifampicin 600 mg once daily and dolutegravir 50 mg twice daily added",
      "with sampling at visit 2 (day 21); then ATV/r increased to 300/100 mg",
      "twice daily with sampling at visit 3 (day 28); then rifampicin",
      "increased to 1200 mg once daily with sampling at visit 4 (day 35)."
    ),
    regions        = "Uganda (Joint Clinical Research Centre, Kampala)",
    notes          = paste(
      "DERIVE (NCT04121195), an open-label, single-arm, dose-escalation",
      "study. Plasma sampling at predose and 0.5, 1, 2, 4, 6, 8 and 12 h",
      "postdose at every visit, with an extra 24 h sample at visit 1.",
      "857 plasma concentrations entered the pooled analysis of both drugs,",
      "of which 28 (3%) atazanavir samples were below the limit of",
      "quantification. Separate PBMC trough samples were taken at visits 1, 3",
      "and 4 and at 12 h postdose at visit 2. Lower limits of quantification",
      "were 0.030 mg/L for plasma atazanavir and 0.015 mg/L for the",
      "intracellular assay. Baseline characteristics are Table 1."
    )
  )

  ini({
    # --- Structural parameters. Final estimates are Kengo 2025 Table 2,
    # atazanavir column, with 95% confidence intervals from sampling
    # importance resampling. Values are quoted for a typical participant of
    # 42 kg fat-free mass on ATV/r 300/100 mg once daily WITHOUT rifampicin
    # (DERIVE visit 1), which is the model's reference level. Where a
    # parameter also appears in the supplementary NONMEM control stream it is
    # cross-referenced by THETA number, but that stream lists initial values
    # carried from a previous run, so Table 2 governs.
    lcl <- log(7.57)
    label("Apparent oral clearance CL/F at FFM = 42 kg, ATV/r once daily without rifampicin (L/h)")  # Table 2 CL 7.57 (6.42-9.06); THETA 1
    lvc <- log(77.5)
    label("Apparent central volume of distribution Vc/F at FFM = 42 kg (L)")                         # Table 2 central volume 77.5 (69.3-88.7); THETA 2
    lq <- log(3.13)
    label("Apparent inter-compartmental clearance Q/F at FFM = 42 kg (L/h)")                         # Table 2 inter compartmental clearance 3.13 (2.34-4.31); THETA 9
    lvp <- log(42.1)
    label("Apparent peripheral volume of distribution Vp/F at FFM = 42 kg (L)")                      # Table 2 peripheral volume 42.1 (26.6-79.0); THETA 8
    lka <- fixed(log(6))
    label("First-order absorption rate constant from depot to central without rifampicin (1/h)")     # Table 2 ka '6 Fixed'; THETA 3 FIX
    lmtt <- log(0.499)
    label("Mean transit time through the absorption transit chain (h)")                              # Table 2 MTT 0.499 (0.429-0.574); THETA 7
    lnn <- fixed(log(10))
    label("Number of absorption transit compartments in the Savic chain (unitless)")                 # Table 2 NN '10 Fixed'; THETA 15 FIX
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability at the reference regimen (fraction)")                       # Table 2 bioavailability '1 fixed'; THETA 4 FIX

    # --- PBMC effect-compartment parameters. Kengo 2025 reports the
    # equilibration half-life rather than the rate constant, but the control
    # stream parameterises and estimates the rate constant directly
    # (TVKE0 = THETA(17), initial 0.700687), so ke0 = log(2) / t_half is the
    # exact inversion of the reported quantity: log(2) / 0.963 = 0.71979 1/h.
    lke0 <- log(log(2) / 0.963)
    label("Plasma-to-PBMC equilibration rate constant (1/h)")                                        # Table 2 PBMC equilibration half-life 0.963 h (0.546-1.51); THETA 17
    lppc <- log(0.653)
    label("Plasma-to-PBMC pseudo-partition coefficient (fraction)")                                  # Table 2 PBMC pseudo-partition coefficient 0.653 (0.538-0.797); THETA 18

    # --- Regimen covariate effects. The source encodes rifampicin through a
    # three-level regimen categorical built in $PK as CL_VISIT (1 at visit 1,
    # 1 + THETA(12) at visit 2, 1 + THETA(13) at visits 3 and 4) and
    # BIO_VISIT (THETA(4) = 1 fixed except THETA(14) at visit 2), plus a
    # RIF_KA factor of 1 + THETA(16) whenever rifampicin is present. Table 2
    # prints the resulting fold-changes and percentage changes, which are the
    # values used here; the control stream initials reproduce them
    # (1 + 2.08781 = 3.09, 1 + 1.01864 = 2.02, 1 - 0.784495 = 0.216).
    e_rif_qd_cl <- 3.05
    label("Fold-change in CL when rifampicin is added to once-daily ATV/r (-fold)")                  # Table 2 fold-change in CL for ATV/r QD + RIF 3.05 (2.67-3.45); THETA 12
    e_rif_bid_cl <- 2.03
    label("Fold-change in CL when rifampicin is given with twice-daily ATV/r (-fold)")               # Table 2 fold-change in CL for ATV/r BID + RIF 2.03 (1.82-2.25); THETA 13
    e_rif_qd_fdepot <- -0.525
    label("Change in bioavailability when rifampicin is added to once-daily ATV/r (fraction)")       # Table 2 change in F for ATV/r QD + RIF -52.5% (-62.5 to -41.4); THETA 14
    e_rif_ka <- -0.673
    label("Change in ka on rifampicin, both once- and twice-daily ATV/r (fraction)")                 # Table 2 change in ka due to RIF -67.3% (-76.2 to -53.2); THETA 16
    e_selfadmin_iov <- 1.63
    label("Multiplier on the between-occasion etas of ka, MTT and F for a self-administered dose (-fold)")  # Table 2 scaling factor on BOV for unobserved dose 1.63 (1.26-2.10); THETA 10

    # --- Allometric exponents, fixed by the authors rather than estimated.
    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent of fat-free mass on CL and Q (unitless)")                             # Kengo 2025 Methods 2.3 allometric scaling; control stream ALLMCL_FFM = (FFM/42)**0.75
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on Vc and Vp (unitless)")                            # Kengo 2025 Methods 2.3 allometric scaling; control stream ALLMV_FFM = (FFM/42)

    # --- Random effects. Kengo 2025 Table 2 footnote c defines the tabulated
    # percentages as %CV = sqrt(omega^2) * 100, i.e. they are the omega
    # standard deviation on the log scale and NOT a log-normal CV. This is
    # confirmed against the control stream $OMEGA block: sqrt(0.0734386) =
    # 0.271 reproduces the 27.6% BSV in clearance, sqrt(0.861212) = 0.928 the
    # 97.9% BOV in ka, sqrt(0.357417) = 0.598 the 59.4% BOV in MTT,
    # sqrt(0.237034) = 0.487 the 48.2% BOV in F, and sqrt(0.0339448) = 0.184
    # the 17.5% BVV in clearance. Variances below are therefore
    # (percentage / 100)^2 using the Table 2 final estimates.
    #
    # The control stream fixes to zero every between-subject eta except the
    # one on clearance (OMEGA 2-9 are 'BLOCK(1) FIX 0'), and also fixes the
    # between-occasion eta on clearance (OMEGA 10) to zero, so clearance
    # carries between-subject and between-visit variability only. The
    # effect-compartment etas BSVKE0 and BVVPPC are likewise FIX 0.
    etalcl ~ 0.076176
    label("Between-subject variability in clearance (log-scale variance)")                           # Table 2 BSV in clearance 27.6% (19.7-36.7); 0.276^2

    # Between-visit variability on clearance across the four PK visits, a
    # single shared variance in the source ($OMEGA BLOCK(1) 0.0339448 followed
    # by three BLOCK(1) SAME repeats for ETA(42)-ETA(45)). nlmixr2 has no SAME
    # shortcut, so visits 2-4 are fix()-pinned to visit 1 (the
    # Abdelgawad_2024_linezolid / Svensson_2018_rifampicin pattern).
    etabvv_cl_1 ~ 0.030625
    label("Between-visit variability in clearance, visit 1 (log-scale variance)")                    # Table 2 BVV in clearance 17.5% (14.4-24.5); 0.175^2
    etabvv_cl_2 ~ fix(0.030625)                                                                      # OMEGA 43 equal to OMEGA 42 per '$OMEGA BLOCK(1) SAME'
    etabvv_cl_3 ~ fix(0.030625)                                                                      # OMEGA 44 equal to OMEGA 42 per '$OMEGA BLOCK(1) SAME'
    etabvv_cl_4 ~ fix(0.030625)                                                                      # OMEGA 45 equal to OMEGA 42 per '$OMEGA BLOCK(1) SAME'

    # Between-occasion variability on ka, MTT and F over the eight dosing
    # occasions, each a single shared variance ($OMEGA BLOCK(1) plus seven
    # SAME repeats).
    etaiov_ka_1 ~ 0.958441
    label("Between-occasion variability in ka, occasion 1 (log-scale variance)")                     # Table 2 BOV in ka 97.9% (80.8-125); 0.979^2
    etaiov_ka_2 ~ fix(0.958441)                                                                      # OMEGA 27 equal to OMEGA 26 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_3 ~ fix(0.958441)                                                                      # OMEGA 28 equal to OMEGA 26 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_4 ~ fix(0.958441)                                                                      # OMEGA 29 equal to OMEGA 26 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_5 ~ fix(0.958441)                                                                      # OMEGA 30 equal to OMEGA 26 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_6 ~ fix(0.958441)                                                                      # OMEGA 31 equal to OMEGA 26 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_7 ~ fix(0.958441)                                                                      # OMEGA 32 equal to OMEGA 26 per '$OMEGA BLOCK(1) SAME'
    etaiov_ka_8 ~ fix(0.958441)                                                                      # OMEGA 33 equal to OMEGA 26 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_1 ~ 0.352836
    label("Between-occasion variability in mean transit time, occasion 1 (log-scale variance)")      # Table 2 BOV in MTT 59.4% (49.3-73.6); 0.594^2
    etaiov_mtt_2 ~ fix(0.352836)                                                                     # OMEGA 35 equal to OMEGA 34 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_3 ~ fix(0.352836)                                                                     # OMEGA 36 equal to OMEGA 34 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_4 ~ fix(0.352836)                                                                     # OMEGA 37 equal to OMEGA 34 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_5 ~ fix(0.352836)                                                                     # OMEGA 38 equal to OMEGA 34 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_6 ~ fix(0.352836)                                                                     # OMEGA 39 equal to OMEGA 34 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_7 ~ fix(0.352836)                                                                     # OMEGA 40 equal to OMEGA 34 per '$OMEGA BLOCK(1) SAME'
    etaiov_mtt_8 ~ fix(0.352836)                                                                     # OMEGA 41 equal to OMEGA 34 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_1 ~ 0.232324
    label("Between-occasion variability in bioavailability, occasion 1 (log-scale variance)")        # Table 2 BOV in F 48.2% (40.7-53.6); 0.482^2
    etaiov_fdepot_2 ~ fix(0.232324)                                                                  # OMEGA 19 equal to OMEGA 18 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_3 ~ fix(0.232324)                                                                  # OMEGA 20 equal to OMEGA 18 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_4 ~ fix(0.232324)                                                                  # OMEGA 21 equal to OMEGA 18 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_5 ~ fix(0.232324)                                                                  # OMEGA 22 equal to OMEGA 18 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_6 ~ fix(0.232324)                                                                  # OMEGA 23 equal to OMEGA 18 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_7 ~ fix(0.232324)                                                                  # OMEGA 24 equal to OMEGA 18 per '$OMEGA BLOCK(1) SAME'
    etaiov_fdepot_8 ~ fix(0.232324)                                                                  # OMEGA 25 equal to OMEGA 18 per '$OMEGA BLOCK(1) SAME'

    # --- Residual error, one combined proportional-plus-additive model per
    # matrix. The control stream $ERROR builds each additive term as
    # THETA + 0.2 * LLOQ with THETA(6) and THETA(20) both FIX 0, so the
    # tabulated additive errors are exactly 20% of the respective limits of
    # quantification: 0.2 * 0.030 = 0.006 mg/L for plasma and
    # 0.2 * 0.015 = 0.003 mg/L for PBMC, which is what Table 2 prints. This
    # is the 'fixing the additive error to at least 20% of the corresponding
    # LLOQ' rule stated in Kengo 2025 Methods 2.3.
    propSd <- 0.198
    label("Proportional residual error for plasma atazanavir (fraction)")                            # Table 2 proportional error (plasma) 19.8% (18.2-21.1)
    addSd <- fixed(0.006)
    label("Additive residual error for plasma atazanavir (mg/L)")                                    # Table 2 additive error (plasma) '0.006 Fixed'; 20% of the 0.030 mg/L plasma LLOQ
    propSd_Cpbmc <- 0.749
    label("Proportional residual error for intracellular PBMC atazanavir (fraction)")                # Table 2 proportional error PBMC 74.9% (62.4-92.0)
    addSd_Cpbmc <- fixed(0.003)
    label("Additive residual error for intracellular PBMC atazanavir (mg/L)")                        # Table 2 additive error PBMC '0.003 fixed'; 20% of the 0.015 mg/L intracellular LLOQ
  })

  model({
    # 1. Occasion indicators. The source control stream multiplexes the
    # between-occasion etas with IF (OCC==k) blocks over eight occasions.
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

    # Doses that were not directly observed inflate the between-occasion etas
    # of all three absorption parameters, per $PK
    #   IF (OBS.EQ.0) THEN BOVKA = SCALE_BOV*BOVKA ... ENDIF
    iov_scale <- 1 + (e_selfadmin_iov - 1) * SELFADMIN

    # Between-visit variability on clearance. The control stream selects
    # ETA(42)-ETA(45) on a separate VISIT column; the DERIVE design pairs two
    # dosing occasions with each of the four PK visits, so the visit level is
    # derived here from OCC and no second column is needed.
    bvv_cl <- (oc1 + oc2) * etabvv_cl_1 + (oc3 + oc4) * etabvv_cl_2 +
      (oc5 + oc6) * etabvv_cl_3 + (oc7 + oc8) * etabvv_cl_4

    # 2. Regimen covariate factors. CONMED_RIF and REGI_BID together select
    # one of the three regimen levels the source fitted; with CONMED_RIF = 0
    # every factor collapses to 1 and the reference regimen is recovered.
    rif_qd <- CONMED_RIF * (1 - REGI_BID)
    rif_bid <- CONMED_RIF * REGI_BID

    # Clearance rises 3.05-fold on once-daily and 2.03-fold on twice-daily
    # ATV/r with rifampicin; exactly one of the two branches is active.
    cl_rif <- 1 + (e_rif_qd_cl - 1) * rif_qd + (e_rif_bid_cl - 1) * rif_bid
    # Bioavailability falls 52.5% on once-daily ATV/r with rifampicin and is
    # fully restored on twice-daily ATV/r (Table 2 prints no BID + RIF entry
    # and Results 3.4 states that doubling the dosing frequency restored
    # atazanavir bioavailability), so only the QD + RIF branch appears.
    fdepot_rif <- 1 + e_rif_qd_fdepot * rif_qd
    # Absorption slows by the same 67.3% in both rifampicin arms.
    ka_rif <- 1 + e_rif_ka * CONMED_RIF

    # 3. Individual parameters. Clearance and both volumes are allometrically
    # scaled on fat-free mass against the 42 kg reference; ka, MTT and NN are
    # not size-scaled.
    cl <- exp(lcl + etalcl + bvv_cl) * (FFM / 42)^e_ffm_cl * cl_rif
    vc <- exp(lvc) * (FFM / 42)^e_ffm_vc
    q <- exp(lq) * (FFM / 42)^e_ffm_cl
    vp <- exp(lvp) * (FFM / 42)^e_ffm_vc
    ka <- exp(lka + iov_ka * iov_scale) * ka_rif
    mtt <- exp(lmtt + iov_mtt * iov_scale)
    nn <- exp(lnn)
    fdepot <- exp(lfdepot + iov_fdepot * iov_scale) * fdepot_rif
    ke0 <- exp(lke0)
    ppc <- exp(lppc)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. Savic transit-compartment absorption: the control stream
    # computes the gamma-density input analytically in $DES,
    #   KTR     = (NN + 1) / MTT
    #   TRANSIT = EXP(LOG(BIO*PD*KTR) - GAMLN(NN+1) + NN*LOG(KTR*TEMPO)
    #             - KTR*TEMPO)
    #   DADT(1) = TRANSIT - KA*A(1)
    # with PD the most recent dose amount and TEMPO the time after that dose.
    # rxode2's transit() is exactly this closed form, using podo() and tad()
    # internally. The control stream sets F1 = 0 so the dose does not also
    # arrive in the depot as a bolus; f(depot) <- 0 preserves that.
    Cc <- central / vc

    d/dt(depot) <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # The PBMC state holds a concentration, not an amount: $DES integrates
    # DADT(4) = KE0*(PPC*C2 - A(4)) against C2 = A(2)/V, and $ERROR reads the
    # state back directly as CC = A(4).
    d/dt(effect) <- ke0 * (ppc * Cc - effect)

    f(depot) <- 0

    # 5. Observations and residual error.
    Cpbmc <- effect
    Cc ~ add(addSd) + prop(propSd)
    Cpbmc ~ add(addSd_Cpbmc) + prop(propSd_Cpbmc)
  })
}
