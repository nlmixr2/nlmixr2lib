Abdelgawad_2024_linezolid <- function() {
  description <- paste(
    "One-compartment population PK model for oral linezolid in plasma and",
    "lumbar cerebrospinal fluid (CSF) in adults with HIV-associated",
    "tuberculous meningitis receiving high-dose rifampicin (35 mg/kg)",
    "co-treatment (LASER-TBM PK substudy, Abdelgawad 2024). Absorption is a",
    "Savic transit chain (mean transit time 0.211 h, estimated chain length",
    "5.68) feeding a first-order depot (ka 1.21 1/h); elimination from the",
    "central compartment is saturable Michaelis-Menten with maximal clearance",
    "7.25 L/h and Michaelis-Menten constant 27.2 mg/L, so Vmax = CLmax * km =",
    "197 mg/h. CLmax and central volume (40.8 L) are allometrically scaled on",
    "fat-free mass with a 45 kg reference and fixed 0.75 / 1 exponents (FFM",
    "beat total body weight by dOFV -30 vs -7.7). CSF is a Sheiner-style",
    "effect compartment holding a concentration, equilibrating with plasma at",
    "0.198 1/h (equilibration half-life 3.5 h) toward a pseudo-partition",
    "coefficient PPC; PPC rises linearly with CSF total protein and plateaus",
    "at PPCmax = 0.365 once CSF protein reaches an estimated 1.18 g/L",
    "breakpoint (a broken-stick whose intercept and amplitude were fixed to 0",
    "and 1, so the slope 0.847 = 1/1.18 is fully determined by the",
    "breakpoint). Random effects are between-subject variability on CLmax",
    "(9.60%), between-visit variability on CLmax across the day-3 and day-28",
    "PK visits (20.3%), and five-occasion between-occasion variability on ka",
    "(87.9%) and mean transit time (110%); all reported percentages are the",
    "omega standard deviation on the log scale. Residual error is combined",
    "proportional plus additive, separately for plasma (21.5%, 0.173 mg/L)",
    "and CSF (91.5%, 0.02 mg/L)."
  )
  reference <- paste(
    "Abdelgawad N, Wasserman S, Abdelwahab MT, Davis A, Stek C, Wiesner L,",
    "Black J, Meintjes G, Wilkinson RJ, Denti P (2024).",
    "Linezolid population pharmacokinetic model in plasma and cerebrospinal",
    "fluid among patients with tuberculosis meningitis.",
    "J Infect Dis 229(4):1200-1208. doi:10.1093/infdis/jiad413",
    sep = " "
  )
  vignette <- "Abdelgawad_2024_linezolid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482. The `csf` state is an exception to the usual "states hold an
  # amount" rule: the source control stream's $DES carries the verbatim
  # comment "A(3) IS ACTUALLY CONC IN EFFECT CMT", i.e. the Sheiner effect
  # compartment integrates a concentration directly
  # (dA(3)/dt = KE0*(PPC*C2 - A(3))), so its units are mg/L and not mg.
  compartmentData <- list(
    depot = list(
      analyte = "linezolid", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "linezolid", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    csf = list(
      analyte = "linezolid", units = "mg/L",
      specimen = "CSF", verified = TRUE
    )
  )

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass, computed from sex, total body weight, and height by the Janmahasatian (2005) formula",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size descriptor for both maximal clearance CLmax",
        "(exponent 0.75) and central volume V (exponent 1), with a fixed",
        "reference of 45 kg -- the cohort median fat-free mass (Abdelgawad",
        "2024 Table 1 'Fat-free mass 45 (30-59) kg'; control-stream $PK",
        "'TVFFM = 45'). FFM was selected over total body weight on model",
        "fit (Results, PK Modeling: dOFV = -30 for FFM vs -7.7 for total",
        "body weight). Height was missing for 18 of 30 participants (Table 1",
        "footnote a); the authors imputed it from sex and weight before",
        "computing FFM. The control stream's $PK imputation model is",
        "HT = (0.00133 * WT + 1.51) * exp(eta) for females and",
        "HT = (0.00281 * WT + 1.53) * exp(eta) for males, with the imputation",
        "etas carrying fixed variances 0.00215 (female) and 0.00170 (male);",
        "FFM then follows Janmahasatian,",
        "FFM = 37.99 * HT^2 * WT / (35.98 * HT^2 + WT) for females and",
        "FFM = 42.92 * HT^2 * WT / (30.93 * HT^2 + WT) for males, with HT in",
        "m and WT in kg. That imputation is a missing-data device for the",
        "original fit, not part of the structural model; users should supply",
        "a measured or Janmahasatian-derived FFM column directly."
      ),
      source_name        = "FFMNEW"
    ),
    CSF_TPRO = list(
      description        = "Total protein concentration in lumbar cerebrospinal fluid",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives the pseudo-partition coefficient PPC through a piecewise-",
        "linear (broken-stick) function: PPC rises linearly from 0 and",
        "plateaus at PPCmax once CSF total protein reaches the estimated",
        "1.18 g/L breakpoint (Abdelgawad 2024 Table 2 and footnote e;",
        "Figure 2). The source paper writes this covariate in mg/mL, which",
        "is numerically identical to the canonical g/L (1 mg/mL = 1 g/L), so",
        "no conversion is applied. Time-varying: the cohort median fell from",
        "1.46 g/L at the day-3 visit to 0.750 g/L at the day-28 visit",
        "(Table 1). The control stream imputes missing values to the",
        "population median 0.995 g/L ($PK 'COV_MED = 0.995'), which is also",
        "the value used for the typical participant in the paper's Monte",
        "Carlo simulations (Methods, Simulations). Observed cohort range was",
        "0.22-54.7 g/L; the paper explicitly marks PPC values outside the",
        "observed range as extrapolation (Figure 2 dashed line)."
      ),
      source_name        = "CSF_PROTEIN"
    ),
    OCC = list(
      description        = "Sampling-occasion index used for the between-occasion and between-visit random effects",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Five occasions, matching the control stream's",
        "IF (OCC==1) ... IF (OCC==5) multiplexers over ETA(6)-ETA(10) for ka",
        "and ETA(11)-ETA(15) for mean transit time. The occasion-to-visit map",
        "is deterministic in the source: occasions 1, 2, and 5 belong to the",
        "day-3 PK visit and occasions 3 and 4 to the day-28 PK visit (the",
        "control stream's commented per-occasion assignments of BVVCLMAX to",
        "ETA(23) vs ETA(24), which the active code selects equivalently via",
        "IF (PK_VISIT == 3) / IF (PK_VISIT == 28)). Because the map is",
        "deterministic, no separate visit column is needed: the model derives",
        "the visit level from OCC. For simulation, set OCC to the occasion",
        "index of each dosing interval; a single-occasion simulation may use",
        "OCC = 1 throughout."
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 30L,
    n_studies      = 1L,
    age_range      = "27-56 years (median 40 at the day-3 visit)",
    age_median     = "40 years",
    weight_range   = "30-96 kg (median 58 at the day-3 visit)",
    weight_median  = "58 kg",
    ffm_range      = "30-59 kg (median 45 at the day-3 visit)",
    sex_female_pct = 40,
    race_ethnicity = "Not reported; the cohort was enrolled at four public hospitals in Cape Town and Gqeberha, South Africa",
    disease_state  = paste(
      "HIV-associated tuberculous meningitis (TBM). All participants were",
      "living with HIV; median CD4 count 137 cells/mm3 (range 2-890). Median",
      "CSF total protein 1.46 g/L at day 3 and 0.750 g/L at day 28; median",
      "CSF glucose 2.9 and 3.2 mmol/L. All received adjunctive dexamethasone."
    ),
    dose_range     = paste(
      "Linezolid 1200 mg orally once daily for the first 28 days, then",
      "600 mg once daily to day 56. At the day-3 visit all 30 participants",
      "were on 1200 mg; at the day-28 visit 13 were on 1200 mg and 4 on",
      "600 mg. Co-administered with high-dose rifampicin (35 mg/kg) plus",
      "isoniazid, pyrazinamide, and ethambutol, with or without aspirin."
    ),
    regions        = "South Africa (Cape Town and Gqeberha)",
    notes          = paste(
      "LASER-TBM phase IIb open-label trial PK substudy",
      "(ClinicalTrials.gov NCT03927313). Thirty participants underwent PK",
      "sampling on day 3 and 17 on day 28, contributing 247 plasma",
      "observations (6 below the limit of quantification, 2.43%) and 28",
      "lumbar CSF observations (7 BLQ, 25%). Day-3 sampling was intensive",
      "(predose and 0.5, 1, 2, 3, 6, 8-10, and 24 h postdose); day-28",
      "sampling was sparse (predose, 2 and 4 h postdose). One lumbar CSF",
      "sample was taken per visit, with the sampling time randomised across",
      "the 1-3, 3-6, 6-10, and 24 h postdose windows. Baseline",
      "characteristics are Table 1."
    )
  )

  ini({
    # --- Structural parameters. All final estimates are Abdelgawad 2024
    # Table 2 (95% CI and RSE from sampling/importance resampling). Where a
    # value also appears in the supplementary NONMEM control stream it is
    # cross-referenced by THETA number, but the control stream lists
    # initial values carried from a previous run, so Table 2 governs.
    #
    # The model is parameterised on maximal clearance CLmax and the
    # Michaelis-Menten constant km; the control stream forms the ODE rate
    # from VMAX = CLMAX * KM ($PK 'VMAX = CLMAX*KM'). The canonical
    # nlmixr2lib name for a Michaelis-Menten maximum rate is lvmax, so the
    # product is taken here: 7.25 L/h * 27.2 mg/L = 197.2 mg/h. Because km
    # carries no random effect in the source (OMEGA 25 BSVKM is FIX 0), an
    # eta on CLmax and an eta on Vmax are mathematically identical.
    lvmax <- log(7.25 * 27.2)
    label("Michaelis-Menten maximum elimination rate at FFM = 45 kg, Vmax = CLmax * km (mg/h)")  # Table 2 CLmax 7.25 (6.09-8.86) [RSE 9.93%] and km 27.2; THETA 1 x THETA 9
    lkm <- log(27.2)
    label("Michaelis-Menten constant, linezolid concentration at half of CLmax (mg/L)")          # Table 2 km 27.2 (16.0-46.4) [RSE 29.1%]; THETA 9
    lvc <- log(40.8)
    label("Central volume of distribution at FFM = 45 kg (L)")                                   # Table 2 V 40.8 (37.9-43.6) [RSE 3.65%]; THETA 2
    lka <- log(1.21)
    label("First-order absorption rate constant from depot to central (1/h)")                    # Table 2 ka 1.21 (.831-1.76) [RSE 19.6%]; THETA 3
    lmtt <- log(0.211)
    label("Mean transit time through the absorption transit chain (h)")                          # Table 2 MTT 0.211 (.112-.342) [RSE 28.6%]; THETA 4
    lntr <- log(5.68)
    label("Number of absorption transit compartments in the Savic chain (unitless)")             # Table 2 NN 5.68 (2.36-11.8) [RSE 43.5%]; THETA 6
    lfdepot <- fixed(log(1))
    label("Oral bioavailability (fraction)")                                                     # Table 2 F 1 fixed; THETA 5 FIX

    # --- CSF effect-compartment parameters.
    lke0 <- log(0.198)
    label("Plasma-to-CSF equilibration rate constant (1/h)")                                     # Table 2 k plasma-CSF 0.198 (.0849-.340) [RSE 33.7%]; THETA 10
    lppc <- log(0.365)
    label("Maximal pseudo-partition coefficient of CSF relative to plasma (fraction)")           # Table 2 PPCmax 0.365 (.238-.566) [RSE 23.2%]; THETA 11

    # --- Covariate effects.
    # The CSF-protein effect on PPC is a broken stick whose intercept and
    # amplitude were FIXED to 0 and 1 respectively, leaving the breakpoint
    # as the single estimated parameter. The slope 0.847 reported in Table 2
    # footnote e is then fully determined: slope = (amplitude - intercept) /
    # (breakpoint - 0) = 1 / 1.18 = 0.847. See the model() block for the
    # resulting expression and the vignette Errata for the printed-equation
    # discrepancy.
    e_csftpro_ppc <- 1.18
    label("CSF total protein breakpoint at which PPC reaches PPCmax (g/L)")                      # Table 2 CSF protein max 1.18 (.730-1.90) [RSE 24.4%]; THETA 14
    e_ffm_vmax <- fixed(0.75)
    label("Allometric exponent of fat-free mass on Vmax (unitless)")                             # Results PK Modeling allometric scaling; control stream $PK ALLMCL_FFM = (FFM/45)**0.75
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on central volume (unitless)")                   # Results PK Modeling allometric scaling; control stream $PK ALLMV_FFM = (FFM/45)

    # --- Random effects. Table 2 reports each variability as a percentage
    # that is the omega standard deviation on the log scale, NOT a CV%.
    # Verified against the control stream $OMEGA block: sqrt(0.772478) =
    # 0.879 reproduces the reported BOV in ka of 87.9%, and
    # sqrt(0.0411842) = 0.203 reproduces the reported BVV of 20.3%.
    # Variances below are therefore (percentage / 100)^2 using the Table 2
    # final estimates.
    etalvmax ~ 0.009216
    label("Between-subject variability in CLmax (log-scale variance)")                           # Table 2 BSV in CLmax 9.60% (3.44-13.9) [RSE 51.9%]; 0.0960^2

    # Between-visit variability on CLmax. Two visits (day 3 and day 28) with
    # a single shared variance: control stream $OMEGA BLOCK(1) 0.0411842
    # followed by $OMEGA BLOCK(1) SAME for ETA(23) / ETA(24). nlmixr2 has no
    # SAME shortcut, so the second visit's variance is fix()-pinned to the
    # first (the Svensson_2018_rifampicin / Jonsson_2011_ethambutol pattern).
    etabvv_vmax_1 ~ 0.041209
    label("Between-visit variability in CLmax at the day-3 visit (log-scale variance)")          # Table 2 BVV in CLmax 20.3% (15.3-26.9) [RSE 30.7%]; 0.203^2
    etabvv_vmax_2 ~ fix(0.041209)                                                                # OMEGA 24 equal to OMEGA 23 per $OMEGA BLOCK(1) SAME

    # Between-occasion variability on ka and MTT over five occasions, each a
    # single shared variance ($OMEGA BLOCK(1) + four SAME repeats).
    etaiov_ka_1 ~ 0.772641
    label("Between-occasion variability in ka, occasion 1 (log-scale variance)")                 # Table 2 BOV in ka 87.9% (66.4-110) [RSE 25.9%]; 0.879^2
    etaiov_ka_2 ~ fix(0.772641)                                                                  # OMEGA 7 equal to OMEGA 6 per $OMEGA BLOCK(1) SAME
    etaiov_ka_3 ~ fix(0.772641)                                                                  # OMEGA 8 equal to OMEGA 6 per $OMEGA BLOCK(1) SAME
    etaiov_ka_4 ~ fix(0.772641)                                                                  # OMEGA 9 equal to OMEGA 6 per $OMEGA BLOCK(1) SAME
    etaiov_ka_5 ~ fix(0.772641)                                                                  # OMEGA 10 equal to OMEGA 6 per $OMEGA BLOCK(1) SAME
    etaiov_mtt_1 ~ 1.21
    label("Between-occasion variability in mean transit time, occasion 1 (log-scale variance)")  # Table 2 BOV in MTT 110% (75.8-144) [RSE 32.8%]; 1.10^2
    etaiov_mtt_2 ~ fix(1.21)                                                                     # OMEGA 12 equal to OMEGA 11 per $OMEGA BLOCK(1) SAME
    etaiov_mtt_3 ~ fix(1.21)                                                                     # OMEGA 13 equal to OMEGA 11 per $OMEGA BLOCK(1) SAME
    etaiov_mtt_4 ~ fix(1.21)                                                                     # OMEGA 14 equal to OMEGA 11 per $OMEGA BLOCK(1) SAME
    etaiov_mtt_5 ~ fix(1.21)                                                                     # OMEGA 15 equal to OMEGA 11 per $OMEGA BLOCK(1) SAME

    # --- Residual error, one combined proportional-plus-additive model per
    # matrix. The control stream $ERROR builds each additive term as
    # THETA + 0.2 * LLOQ with LLOQ = 0.1 mg/L; Table 2 reports the resulting
    # total, which is confirmed by the CSF row: THETA(13) is FIX 0, so the
    # only way Table 2 can print 0.02 mg/L for CSF is if the tabulated
    # quantity is THETA + 0.2 * 0.1 = 0.02. The plasma row is read the same
    # way and taken directly as 0.173 mg/L.
    propSd <- 0.215
    label("Proportional residual error for plasma linezolid (fraction)")                         # Table 2 Plasma proportional error 21.5% (18.8-24.7) [RSE 7.06%]
    addSd <- 0.173
    label("Additive residual error for plasma linezolid (mg/L)")                                 # Table 2 Plasma additive error 0.173 (.0379-.355) [RSE 47.1%]
    propSd_Ccsf <- 0.915
    label("Proportional residual error for CSF linezolid (fraction)")                            # Table 2 CSF proportional error 91.5% (63.3-151) [RSE 23.4%]
    addSd_Ccsf <- fixed(0.02)
    label("Additive residual error for CSF linezolid (mg/L)")                                    # Table 2 CSF additive error 0.02 fixed; footnote c sets it to 20% of the 0.1 mg/L LLOQ
  })

  model({
    # 1. Occasion indicators. The source control stream multiplexes the
    # between-occasion etas with IF (OCC==k) blocks over five occasions.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)

    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 +
      oc4 * etaiov_ka_4 + oc5 * etaiov_ka_5
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 +
      oc4 * etaiov_mtt_4 + oc5 * etaiov_mtt_5

    # Between-visit variability on CLmax. The control stream selects
    # ETA(23) when PK_VISIT == 3 and ETA(24) when PK_VISIT == 28; its
    # commented per-occasion block gives the deterministic map
    # occasions 1, 2, 5 -> day 3 and occasions 3, 4 -> day 28, so the visit
    # level is derived here from OCC rather than requiring a second column.
    bvv_vmax <- (oc1 + oc2 + oc5) * etabvv_vmax_1 + (oc3 + oc4) * etabvv_vmax_2

    # 2. Individual parameters. CLmax and V are allometrically scaled on
    # fat-free mass against the 45 kg cohort median; because Vmax = CLmax *
    # km and km is not scaled, Vmax inherits the CLmax exponent.
    vmax <- exp(lvmax + etalvmax + bvv_vmax) * (FFM / 45)^e_ffm_vmax
    km <- exp(lkm)
    vc <- exp(lvc) * (FFM / 45)^e_ffm_vc
    ka <- exp(lka + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    ntr <- exp(lntr)
    ke0 <- exp(lke0)
    fdepot <- exp(lfdepot)

    # 3. Pseudo-partition coefficient. The control stream $PK builds this as
    #   SLP     = (AMP - INT) / (BRK - 0)  with AMP = 1 FIX, INT = 0 FIX
    #   COV_EFF = SLP * (COV - BRK)   for COV <= BRK, else 0
    #   TVPPC   = THETA(11) * (1 + COV_EFF)
    # which simplifies to PPC = PPCmax * min(CSF_TPRO / breakpoint, 1): the
    # coefficient rises linearly from 0 and plateaus at PPCmax. Table 2
    # footnote e prints this without the leading 1, which would make PPC
    # zero exactly at the breakpoint instead of maximal; the control stream
    # is definitive and the paper's own text (a 3% rise in PPC per
    # 0.1 mg/mL of CSF protein: 0.365 / 1.18 * 0.1 = 0.0309) corroborates
    # the form used here. See the vignette Errata.
    ppc <- exp(lppc) * min(CSF_TPRO / e_csftpro_ppc, 1)

    # 4. ODE system. Plasma concentration drives the saturable elimination
    # and the CSF effect compartment.
    Cc <- central / vc

    # Savic transit-compartment absorption. The control stream computes the
    # gamma-density input analytically in $DES,
    #   KTR     = (NN + 1) / MTT
    #   TRANSIT = EXP(LOG(BIO*PD*KTR) - GAMLN(NN+1) + NN*LOG(KTR*TEMPO)
    #             - KTR*TEMPO)
    #   DADT(1) = TRANSIT - KA*A(1)
    # with PD the most recent dose amount and TEMPO the time after that
    # dose. rxode2's transit() is exactly this closed form, using podo() and
    # tad() internally. The control stream sets F1 = 0 so the dose does not
    # also accumulate in the depot as a bolus; f(depot) <- 0 preserves that.
    d/dt(depot) <- transit(ntr, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - vmax * Cc / (km + Cc)
    # The CSF state holds a concentration, not an amount: $DES carries the
    # verbatim comment "A(3) IS ACTUALLY CONC IN EFFECT CMT".
    d/dt(csf) <- ke0 * (ppc * Cc - csf)

    f(depot) <- 0

    # 5. Observations and residual error.
    Ccsf <- csf
    Cc ~ add(addSd) + prop(propSd)
    Ccsf ~ add(addSd_Ccsf) + prop(propSd_Ccsf)
  })
}
