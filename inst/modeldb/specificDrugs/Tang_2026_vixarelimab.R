Tang_2026_vixarelimab <- function() {
  description <- "Two-compartment target-mediated drug-disposition (TMDD) model with a quasi-steady-state (QSS) approximation for subcutaneous and intravenous vixarelimab (anti-OSMR-beta monoclonal antibody) in healthy volunteers and patients with chronic pruritic conditions, with first-order absorption, logit-scale bioavailability and power body-weight effects on CL, Vc and Vp (Tang 2026, Table 3)"
  reference <- "Tang F, Dang S, Arjomandi A, Pan L, Hsu J, Kshirsagar S, Kassir N. Population Pharmacokinetic Analysis of Vixarelimab in Healthy Volunteers and Patients With Chronic Pruritic Conditions. CPT Pharmacometrics Syst Pharmacol. 2026;15(3):e70230. doi:10.1002/psp4.70230"
  vignette <- "Tang_2026_vixarelimab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The authors' NONMEM model works in MOLAR units (see the
  # unit-conversion note in model() below): the three drug states carry nmol and
  # `total_target` carries a CONCENTRATION in nM, exactly as NONMEM compartment 3
  # does (Supporting Information $DES: DADT(3) = KSYN - KDEG*A(3) - ...,
  # with KSYN = KDEG * R0 in nM/h and A_0(3) = R0 in nM).
  compartmentData <- list(
    depot        = list(analyte = "vixarelimab", units = "nmol", specimen = "administration site", verified = TRUE),
    central      = list(analyte = "vixarelimab", units = "nmol", specimen = "serum", verified = TRUE),
    peripheral1  = list(analyte = "vixarelimab", units = "nmol", specimen = "serum", verified = TRUE),
    total_target = list(analyte = "oncostatin M receptor beta (OSMR-beta), free + drug-bound", units = "nM", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline (time-fixed) weight. Entered as a power model centred on 81.7 kg on CL, Vc and Vp",
        "(Supporting Information $PK: CLWT = (WT/81.7)**THETA(11), and likewise for Vc and Vp).",
        "The 81.7 kg centring constant appears only in the deposited control stream; the paper text",
        "reports a cohort median weight of 82 kg (Table 2) and describes a Figure 4 reference subject",
        "of 82.1 kg. Body weight was the only statistically significant covariate retained after",
        "backward elimination (Section 3.2)."
      ),
      source_name        = "WT"
    )
  )

  # Screened in the pre-specified covariate analysis (Table S2) but NOT retained
  # in the final model, so they are documentation only and are never referenced
  # in model(). Age on Vc entered in forward selection (dOFV = -10.2) and was
  # removed again in backward elimination (Table S3, run 6).
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = "Tested on CL and Vc (Table S2). Entered forward selection on Vc but was dropped in backward elimination (Table S3 run 6); no point estimate is reported for the final model."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL and Vc (Table S2, coded SEX). Not statistically significant (Section 4); no point estimate reported."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = all other races)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL and Vc as 'Asian versus others' (Section 2.4, Table S2). Not statistically significant (Section 4); no point estimate reported."
    ),
    ADA_POS = list(
      description = "Anti-drug-antibody-positive indicator (1 = at least one positive post-baseline result, 0 = never positive)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL (Table S2, coded ATA; 18% positive per Table 2). No statistically significant effect on CL (Section 4); no point estimate reported."
    ),
    DISEASE_STATUS = list(
      description = "Patient population indicator (healthy volunteer versus patient with a chronic pruritic skin disorder)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL, F1, Ka and R0 (Table S2, coded DSSTAT). Skin disease was not a significant covariate (Section 4); no point estimate reported."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 274,
    n_studies      = 3,
    n_observations = "4032 serum concentrations collected; 3328 quantifiable post-dose samples used for model fitting (Section 3.1)",
    age_range      = "18-79 years",
    age_median     = "52 years",
    weight_range   = "48-158 kg",
    weight_median  = "82 kg",
    sex_female_pct = 58,
    race_ethnicity = c(White = 68, Black = 18, Asian = 9.9, Other = 4.3),
    disease_state  = "Healthy volunteers (14%) and patients with chronic pruritic skin conditions (86%): moderate-to-severe atopic dermatitis, prurigo nodularis, and chronic idiopathic pruritus, chronic idiopathic urticaria, lichen planus, lichen simplex chronicus or plaque psoriasis",
    dose_range     = "IV single doses 0.3, 1.5, 5, 7.5, 10 and 20 mg/kg; SC single doses 1.5 mg/kg and 360 mg; SC multiple doses 360 mg QW, 720 mg SC loading followed by 360 mg QW, and 120/360/540 mg Q4W (Table 1)",
    regions        = "Not reported",
    immunogenicity = "50 of 274 participants (18%) had at least one positive post-baseline anti-drug-antibody result (Table 2)",
    notes          = "Pooled analysis of the Phase 1 study KPL-716-C001 (Parts 1, 3 and 4) and the Phase 2 studies KPL-716-C201 (NCT03816891, Phases 2a and 2b) and KPL-716-C202 (NCT03858634). Baseline demographics from Table 2; study designs and dose regimens from Table 1."
  )

  ini({
    # ---- Structural disposition (typical 81.7 kg participant) ----
    # All final estimates are taken from Table 3 (the published final model, run 7).
    # The deposited run-7 control stream lists $THETA records that are run-7 INITIAL
    # estimates carried over from run 6: they agree with Table 3 for nine of eleven
    # thetas but differ for Vc (2.9 vs 3.04) and Kint (0.000146168 vs 0.000141), so
    # Table 3 is used throughout.
    lvc <- log(3.04);       label("Central volume of distribution (Vc, L)")                      # Table 3 row "Central volume of distribution, Vc (L)"
    lcl <- log(0.00649);    label("Linear (non-specific) clearance (CL, L/h)")                   # Table 3 row "Clearance, CL (L/h)"
    lq  <- log(0.0207);     label("Inter-compartmental clearance (Q, L/h)")                      # Table 3 row "Inter-compartmental clearance, Q (L/h)"
    lvp <- log(1.74);       label("Peripheral volume of distribution (Vp, L)")                   # Table 3 row "Peripheral volume of distribution, Vp (L)"
    lka <- log(0.0126);     label("First-order subcutaneous absorption rate constant (Ka, 1/h)") # Table 3 row "Absorption rate constant, Ka (h-1)"

    # Subcutaneous bioavailability, estimated on the LOGIT scale by the authors
    # ($PK: TVPHI = LOG(TVBIO/(1-TVBIO)); PHI = TVPHI + ETA(10); BIO = EXP(PHI)/(1+EXP(PHI))),
    # so the IIV below is additive on the logit scale rather than multiplicative.
    logitfdepot <- log(0.566 / (1 - 0.566)); label("Subcutaneous bioavailability on the logit scale (logit(F1), unitless)")  # Table 3 row "Subcutaneous bioavailability, F (%)" = 56.6%

    # ---- Target (OSMR-beta) turnover and QSS binding ----
    lrbase <- log(1.84);      label("Baseline total target (OSMR-beta) concentration (R0, nM)")             # Table 3 row "Baseline receptor level, R0 (nM)"
    lkdeg  <- log(0.097);     label("Target degradation rate constant (Kdeg, 1/h)")                         # Table 3 row "Receptor degradation rate constant, Kdeg (h-1)"
    lkss   <- log(0.0296);    label("Quasi-steady-state constant for drug-target binding (Kss, nM)")        # Table 3 row "Quasi steady-state rate constant, Kss (nM)"
    lkint  <- log(0.000141);  label("Drug-target complex internalisation rate constant (Kint, 1/h)")        # Table 3 row "Complex elimination rate constant Kint (h-1)"

    # ---- Body-weight covariate effects (power model centred on 81.7 kg) ----
    e_wt_cl <- 0.943; label("Power exponent of WT/81.7 on CL (unitless)")  # Table 3 row "Weight on CL"
    e_wt_vc <- 0.762; label("Power exponent of WT/81.7 on Vc (unitless)")  # Table 3 row "Weight on Vc"
    e_wt_vp <- 0.646; label("Power exponent of WT/81.7 on Vp (unitless)")  # Table 3 row "Weight on Vp"

    # ---- Interindividual variability ----
    # Table 3 reports IIV as CV%, defined in its footnote b as
    # CV% = sqrt(exp(omega^2) - 1) * 100%, so omega^2 = log(1 + (CV/100)^2).
    # The diagonal $OMEGA of the deposited run-7 stream (0.0321914, 0.0599205,
    # 0.609013, 0.167648, 0.232025, 0.162611) reproduces these CV% to the three
    # significant figures printed in Table 3; the Table 3 back-transforms are used
    # here for consistency with the other final estimates.
    # (Trailing comments on these lines are deliberately free of double quotes:
    # rxode2 promotes an unlabelled ini() trailing comment into label(...), and an
    # embedded quote breaks the generated label.)
    etalvc          ~ 0.0322358  # Table 3 between subject variability for Vc, CV% = 18.1%
    etalcl          ~ 0.0601549  # Table 3 between subject variability for CL, CV% = 24.9%
    etalq           ~ 0.6092524  # Table 3 between subject variability for Q, CV% = 91.6%
    etalrbase       ~ 0.1682091  # Table 3 between subject variability for R0, CV% = 42.8%
    etalka          ~ 0.2320010  # Table 3 between subject variability for Ka, CV% = 51.1%
    etalogitfdepot  ~ 0.1631736  # Table 3 between subject variability for F1, CV% = 42.1%; additive on the LOGIT scale

    # IIV on Vp approached zero, was poorly estimated, and was fixed to 0 to give
    # the final model (Section 3.2; $OMEGA record "0 FIX ; ETA_VP"). Note that
    # Table S3 describes run 7 as "Removing IIV on Q", which contradicts Section 3.2,
    # Table 3 (which reports an IIV for Q but none for Vp) and the deposited control
    # stream; the three agreeing sources are followed. See the vignette Errata.
    etalvp ~ fixed(0)  # Section 3.2; Supporting Information $OMEGA record 0 FIX on ETA_VP

    # Kdeg, Kss and Kint carry no IIV in the final model
    # (Supporting Information $OMEGA: "0 FIX" on ETA_KDEG, ETA_KSS and ETA_KINT).

    # ---- Residual error ----
    # The authors fitted log-transformed concentrations with an additive error on
    # the log scale ($ERROR: IPRED = LOG(CFREE + DEL); Y = IPRED + EPS(1)), which is
    # a log-normal residual on the linear scale, i.e. `~ lnorm()` in nlmixr2.
    # Table 3 reports it as CV% = sqrt(exp(sigma^2) - 1) * 100% = 19.6%, so
    # sigma^2 = log(1 + 0.196^2) = 0.0376965 and sigma = 0.1941558.
    expSd <- 0.1941558; label("Log-scale (log-additive) residual standard deviation")  # Table 3 row "Log-additive residual variability (CV%)" = 19.6%
  })
  model({
    # ---- Unit conversion (molar model, mass-unit dose and observation) ----
    # The authors' NONMEM model is entirely MOLAR: R0 and Kss are estimated in nM
    # (Table 3), the QSS algebra subtracts them from the drug concentration
    # (A(2)/VC - A(3) - KSS), and the deposited $INPUT carries both the mass columns
    # (AMTMG, DVUGML) and the molar columns (AMT, DVNM) that the authors' SAS
    # derivation produced. The mg <-> nmol factor itself is NOT reported in the
    # paper or in the Supporting Information.
    #
    # NON-PAPER VALUE: mw is set to the nominal molecular weight of a human IgG
    # monoclonal antibody, 150000 g/mol. Vixarelimab's molecular weight is not
    # stated anywhere in the paper or supplement. This follows the same convention
    # as Papachristos_2020_bevacizumab_qss.R and is documented in the vignette's
    # "Assumptions and deviations" section. The choice is benign for the exposures
    # that matter clinically: once the target is saturated the model reduces to a
    # linear two-compartment system in which the mass-unit predictions are exactly
    # MW-invariant (mw cancels between the dose conversion and the observation
    # conversion), so mw only shifts where the target-mediated arm saturates.
    mw <- 150000
    nmol_per_mg <- 1e6 / mw

    # ---- Individual parameters ----
    # Power body-weight model centred on 81.7 kg (Supporting Information $PK:
    # CLWT = (WT/81.7)**THETA(11); VCWT and VPWT identically with THETA(12)/THETA(13)).
    cl <- exp(lcl + e_wt_cl * log(WT / 81.7) + etalcl)
    vc <- exp(lvc + e_wt_vc * log(WT / 81.7) + etalvc)
    vp <- exp(lvp + e_wt_vp * log(WT / 81.7) + etalvp)
    q  <- exp(lq + etalq)
    ka <- exp(lka + etalka)

    rbase <- exp(lrbase + etalrbase)
    kdeg  <- exp(lkdeg)
    kss   <- exp(lkss)
    kint  <- exp(lkint)

    # Logit-scale bioavailability; the fixed effect and eta are collected on their
    # own line so the term stays in a mu-referenced position.
    logitfdepot_ind <- logitfdepot + etalogitfdepot
    fdepot <- exp(logitfdepot_ind) / (1 + exp(logitfdepot_ind))

    # ---- Micro-constants (Section 2.4 definitions) ----
    kel  <- cl / vc   # elimination rate constant from the central compartment
    kpt  <- q / vc    # central -> peripheral transfer rate constant
    ktp  <- q / vp    # peripheral -> central transfer rate constant
    ksyn <- kdeg * rbase  # zero-order target synthesis rate (nM/h)

    # ---- QSS algebra (Equation 6) ----
    # Ctot = Ac/Vc is the TOTAL drug concentration in the central compartment;
    # Cfree is the free drug concentration obtained from the positive root of the
    # quasi-steady-state binding quadratic.
    ctot  <- central / vc
    qss_d <- ctot - total_target - kss
    cfree <- 0.5 * (qss_d + sqrt(qss_d * qss_d + 4 * kss * ctot))

    # ---- Dosing ----
    # Doses are supplied in mg and converted to nmol. Subcutaneous doses go into
    # `depot` and additionally carry the estimated bioavailability; intravenous
    # doses go straight into `central` and are fully bioavailable.
    f(depot)   <- fdepot * nmol_per_mg
    f(central) <- nmol_per_mg

    # ---- ODEs (Equations 2-5) ----
    d/dt(depot)        <- -ka * depot                                                   # Equation 2
    d/dt(central)      <- ka * depot - (kel + kpt) * cfree * vc + ktp * peripheral1 -
      total_target * kint * cfree * vc / (kss + cfree)                                  # Equation 3
    d/dt(peripheral1)  <- kpt * cfree * vc - ktp * peripheral1                          # Equation 4
    d/dt(total_target) <- ksyn - kdeg * total_target -
      (kint - kdeg) * total_target * cfree / (kss + cfree)                              # Equation 5
    total_target(0)    <- rbase  # Supporting Information $PK: A_0(3) = BLR

    # ---- Observation ----
    # The assay does not detect vixarelimab when both binding sites are occupied
    # (Section 2.2), and the authors' $ERROR block predicts the FREE concentration
    # (IPRED = LOG(CFREE + DEL)), so the observed quantity is Cfree, not Ctot.
    # nM -> ug/mL: 1 nM = MW * 1e-9 g/L = MW / 1e6 mg/L = MW / 1e6 ug/mL.
    Cc <- cfree * mw / 1e6
    Cc ~ lnorm(expSd)
  })
}
