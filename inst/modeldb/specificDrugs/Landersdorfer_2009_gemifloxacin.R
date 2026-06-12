Landersdorfer_2009_gemifloxacin <- function() {
  description <- "Two-compartment population PK model for gemifloxacin in healthy adults with first-order absorption + lag time, additive renal (filtration + saturable Michaelis-Menten tubular secretion with competitive probenecid inhibition) and non-renal clearance, and treatment-arm-static probenecid effects on absorption rate, absorption lag, and non-renal clearance (Landersdorfer 2009 gemifloxacin / probenecid DDI study)."
  reference <- paste(
    "Landersdorfer CB, Kirkpatrick CMJ, Kinzig M, Bulitta JB, Holzgrabe U, Drusano GL, Sorgel F.",
    "Competitive inhibition of renal tubular secretion of gemifloxacin by probenecid.",
    "Antimicrob Agents Chemother. 2009;53(9):3902-3907.",
    "doi:10.1128/AAC.01200-08.",
    sep = " "
  )
  vignette <- "Landersdorfer_2009_gemifloxacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CONMED_PROBENECID = list(
      description        = "1 = probenecid co-administered (oral 4.5 g total split into 8 doses across the gemifloxacin sampling window); 0 = gemifloxacin alone (reference). Time-fixed at the treatment-arm level in the source crossover study.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no probenecid; gemifloxacin alone)",
      notes              = "Encodes the treatment-arm switch for the static (non-mechanism-specific) probenecid effects on absorption rate (Ka 0.839 -> 0.897 1/h), absorption lag (Tlag 0.223 -> 0.129 h), and non-renal clearance (CL_NR 25.2 -> 21.0 L/h) reported in Landersdorfer 2009 Table 3. The mechanism-specific competitive inhibition at the renal tubular transporter (apparent Km = Km * (1 + [P] / Kic)) is encoded separately via the time-varying CP_PRB_MGL covariate so that the on / off transition between probenecid pharmacokinetics drives the renal-clearance term smoothly. Source dosing schedule (Landersdorfer 2009 Methods): probenecid 1000 mg at -10 h and -2 h relative to gemifloxacin, 250 mg at +6 h and +14 h, 500 mg at +24 h, +36 h, +48 h, and +60 h.",
      source_name        = "(treatment arm indicator)"
    ),
    CP_PRB_MGL = list(
      description        = "Instantaneous plasma probenecid concentration (mg/L) supplied as a time-varying covariate driving competitive inhibition of the saturable renal tubular secretion of gemifloxacin: apparent Km = Km * (1 + CP_PRB_MGL / Ki).",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Set to 0 throughout the gemifloxacin-alone arm (CONMED_PROBENECID = 0) so the apparent-Km term collapses to Km. Set to the instantaneous plasma probenecid concentration during the gemifloxacin + probenecid arm. The Landersdorfer 2009 paper does NOT report a probenecid PK sub-model parameter table on disk -- the Methods describes the structural form (1-compartment + lag + parallel first-order + mixed-order elimination) but the corresponding parameter estimates are not in any on-disk table -- so users simulating the with-probenecid arm must supply CP_PRB_MGL from an external probenecid popPK source (e.g., literature digitisation of Landersdorfer 2009 Fig. 1C, or an unrelated published probenecid popPK model). Reference peak observed: ~30-60 mg/L during the gemifloxacin sampling window (Landersdorfer 2009 Fig. 1C).",
      source_name        = "[P]"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 17L,
    n_studies        = 1L,
    age_range        = "healthy adults; per-subject ages not tabulated in Landersdorfer 2009.",
    weight_range     = "mean 69.1 +/- 13 kg (Landersdorfer 2009 Results paragraph 1).",
    height_range     = "mean 173 +/- 10 cm (Landersdorfer 2009 Results paragraph 1).",
    sex_female_pct   = 47.1,
    race_ethnicity   = c(White = 100),
    disease_state    = "healthy Caucasian volunteers",
    dose_range       = "single 320 mg oral gemifloxacin dose (Factive 320 mg tablet); 4.5 g total oral probenecid in the co-administration period split into 8 doses (1000 mg at -10 h and -2 h, 250 mg at +6 h and +14 h, 500 mg at +24 h / +36 h / +48 h / +60 h).",
    regions          = "Germany (IBMP-Institute, Nuernberg-Heroldsberg).",
    notes            = "Randomized, two-way crossover with washout >= 7 days between periods. Each subject received gemifloxacin alone in one period and gemifloxacin + probenecid in the other. Plasma sampled at 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 36, 48 h post-gemifloxacin-dose (14 samples per period) plus a predose; urine collected in 10 intervals through 72 h post-dose. WinNonlin (v4.0.1) was used for compartmental modeling and NCA; NONMEM (v V) was used only for visual predictive checks. Baseline demographics from Landersdorfer 2009 Results paragraph 1."
  )

  ini({
    # Structural PK parameters -- reference (without-probenecid) typical values from
    # Landersdorfer 2009 Table 3 ("Without probenecid" column for the absorption /
    # non-renal-clearance parameters that carry an arm-specific point estimate, and
    # the shared row otherwise). All values for gemifloxacin in the final model
    # (Model 1: competitive inhibition of renal tubular secretion + static
    # inhibition of non-renal clearance).
    lka       <- log(0.839); label("Absorption rate constant, gemifloxacin alone (1/h)")             # Table 3, "Without probenecid" row, Ka column
    ltlag     <- log(0.223); label("Absorption lag time, gemifloxacin alone (h)")                    # Table 3, "Without probenecid" row, Tlag column
    lvc       <- log(89.5);  label("Central volume of distribution V1 (L)")                          # Table 3, V1 column (shared between arms)
    lvp       <- log(160);   label("Peripheral volume of distribution V2 (L)")                       # Table 3, V2 column (shared between arms)
    lq        <- log(39.8);  label("Inter-compartmental clearance CL_IC (L/h)")                      # Table 3, CL_IC column (shared between arms)
    lcl_nonren <- log(25.2); label("Non-renal clearance, gemifloxacin alone (L/h)")                  # Table 3, "Without probenecid" row, CL_NR column
    lvmax     <- log(113);   label("Maximum rate of renal tubular secretion Vmax_R (mg/h)")          # Table 3, V_max_R column (shared between arms)
    lkm       <- log(9.16);  label("Half-saturation gemifloxacin concentration for renal tubular secretion Km_R (mg/L)")  # Table 3, Km_R column (shared between arms)
    lki       <- log(69.3);  label("Competitive inhibition constant of probenecid at the renal tubular transporter Kic (mg/L)")  # Table 3, K_ic column (shared between arms; only meaningful when probenecid is co-administered)

    # Probenecid co-administration multiplicative effects on Ka, Tlag, and CL_NR
    # (treatment-arm-static "static inhibition" components from Landersdorfer 2009
    # Methods and Table 3). Encoded as fixed multiplicative offsets computed
    # from the With- and Without-probenecid columns of Table 3; the model is
    # not refit in this extraction so the values are treated as fixed.
    # Coefficient = (with / without) - 1; the model() block applies it as
    # param = param_baseline * (1 + e_prob_<param> * CONMED_PROBENECID).
    e_prob_ka        <- fixed( 0.06913); label("Probenecid effect on absorption rate Ka (fold-change, unitless)")            # 0.897 / 0.839 - 1; Table 3
    e_prob_ltlag     <- fixed(-0.42152); label("Probenecid effect on absorption lag Tlag (fold-change, unitless)")           # 0.129 / 0.223 - 1; Table 3
    e_prob_cl_nonren <- fixed(-0.16667); label("Probenecid effect on non-renal clearance CL_NR (fold-change, unitless)")     # 21.0 / 25.2  - 1; Table 3

    # IIV -- between-subject log-normal variances computed from the
    # Landersdorfer 2009 Table 3 "%coefficient of variation" entries via the
    # standard log-normal conversion omega^2 = log(1 + CV^2). For the
    # treatment-arm-dependent parameters (Ka, Tlag, CL_NR), the
    # without-probenecid arm's CV is used as the population reference
    # (with-probenecid CVs differ slightly -- Ka 104%, Tlag 104%, CL_NR 23% --
    # but only one population CV is carried in the model; see
    # the "Assumptions and deviations" section of the validation vignette).
    etalka       ~ 0.32527  # Ka  CV = 62% (Without arm);   log(1 + 0.62^2) = 0.32527
    etaltlag     ~ 0.73314  # Tlag CV = 104% (both arms);   log(1 + 1.04^2) = 0.73314
    etalcl_nonren ~ 0.06541 # CL_NR CV = 26% (Without arm); log(1 + 0.26^2) = 0.06541
    etalvc       ~ 0.32527  # V1  CV = 62%;                 log(1 + 0.62^2) = 0.32527
    etalvp       ~ 0.09176  # V2  CV = 31%;                 log(1 + 0.31^2) = 0.09176
    etalq        ~ 0.28998  # CL_IC CV = 58%;               log(1 + 0.58^2) = 0.28998
    etalvmax     ~ 0.04316  # Vmax_R CV = 21%;              log(1 + 0.21^2) = 0.04316
    etalkm       ~ 0.03922  # Km_R CV = 20%;                log(1 + 0.20^2) = 0.03922
    etalki       ~ 0.06541  # Kic CV = 26%;                 log(1 + 0.26^2) = 0.06541

    # Residual error -- Landersdorfer 2009 does not tabulate a residual-error
    # magnitude in Table 3 or anywhere else on disk. The bioanalytical assay
    # precision is 3.7 to 7.2% interday CV for gemifloxacin in plasma
    # (Landersdorfer 2009 Methods, "Determination of plasma and urine drug
    # concentrations"), and the compartmental model in WinNonlin would have
    # produced some residual estimate that the paper did not report. A
    # 10% proportional residual is used here as a placeholder consistent
    # with the assay floor plus a modest structural-misspecification buffer;
    # see the "Assumptions and deviations" section of the validation
    # vignette for the unit-interpretation rationale.
    propSd <- 0.10; label("Proportional residual error (fraction; placeholder, see vignette Errata)")  # Not reported; see vignette Errata
  })
  model({
    # Filtration clearance fixed to fu * GFR ~ 2 L/h (Landersdorfer 2009 Methods,
    # "Compartmental modeling": "the range of fu is between 0.3 and 0.4 and all
    # subjects had normal renal function, the renal filtration clearance of
    # gemifloxacin, fu * GFR, is about 2 liters/h and accounts for only about 6%
    # of the total body clearance. Therefore, the first-order component of renal
    # clearance was fixed to 2 liters/h"). Hardcoded model() constant rather
    # than an ini() parameter because it is not estimable and not a primary
    # canonical PK parameter.
    cl_filt <- 2

    # Individual PK parameters. The treatment-arm-static probenecid effects
    # on Ka, Tlag, and CL_NR are applied multiplicatively via the binary
    # CONMED_PROBENECID indicator; the mechanism-specific competitive Kic
    # inhibition on the saturable renal tubular secretion enters the ODE
    # below via the time-varying CP_PRB_MGL covariate.
    ka         <- exp(lka        + etalka)        * (1 + e_prob_ka        * CONMED_PROBENECID)
    tlag       <- exp(ltlag      + etaltlag)      * (1 + e_prob_ltlag     * CONMED_PROBENECID)
    cl_nonren  <- exp(lcl_nonren + etalcl_nonren) * (1 + e_prob_cl_nonren * CONMED_PROBENECID)
    vc         <- exp(lvc        + etalvc)
    vp         <- exp(lvp        + etalvp)
    q          <- exp(lq         + etalq)
    vmax       <- exp(lvmax      + etalvmax)
    km         <- exp(lkm        + etalkm)
    ki         <- exp(lki        + etalki)

    # Plasma gemifloxacin concentration (mg/L); dose in mg, V1 in L => central / vc is mg/L.
    Cc <- central / vc

    # ODE system: 2-compartment disposition with first-order absorption + lag,
    # additive renal (filtration + saturable Michaelis-Menten tubular secretion
    # with competitive probenecid inhibition) + linear non-renal elimination.
    #
    # CL_R(C, [P]) = (fu * GFR) + Vmax / (Km * (1 + [P] / Kic) + C)
    # CL_NR        = CL_NR (treatment-arm-static via CONMED_PROBENECID)
    # Total mass elimination rate from central = (CL_R + CL_NR) * C plus
    # distribution to / from peripheral1 via Q.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - cl_filt * Cc - vmax * Cc / (km * (1 + CP_PRB_MGL / ki) + Cc) - cl_nonren * Cc - q * Cc + q * peripheral1 / vp
    d/dt(peripheral1) <-  q * Cc - q * peripheral1 / vp

    alag(depot) <- tlag

    Cc ~ prop(propSd)
  })
}
