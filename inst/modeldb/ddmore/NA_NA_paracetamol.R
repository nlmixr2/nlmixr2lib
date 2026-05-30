# Mechanistic oral-glucose-tolerance-test (OGTT) model extracted from the
# DDMORE Foundation Model Repository entry DDMODEL00000228. Paracetamol is
# co-administered as a tracer for gastric emptying; its two-compartment
# disposition coupled with a saturable first-pass loss step determines the
# stomach -> small-intestine transit rate KS, which in turn drives glucose
# absorption from a four-segment GI tract (stomach -> duodenum, jejunum,
# ileum) into central and peripheral glucose compartments. Plasma insulin
# (INS) enters as a time-varying regressor through a one-compartment delay
# to drive insulin-dependent glucose clearance, and the incretin hormones
# GLP-1 (responding to duodenal + ileal glucose) and GIP (responding to
# duodenal glucose) are modelled as a pair of indirect-response loops with
# baseline-anchored states.
#
# No publication is linked in the DDMORE bundle (the .rdf carries no
# `model-described-in-literature` URI and the bundle contains no
# Model_Accommodations.txt). The on-disk listing
# `Output_real_run126c.lst` (NONMEM 7.3, FOCEI with eps-eta interaction,
# importance-sampling SE step, OBJ = -7902.844) is the authoritative source
# for the final estimates carried into this model file. Equations come from
# the `Executable_run126h.mod` $PK / $DES / $ERROR blocks. The FOCEI step
# terminated due to rounding errors (NSIG = 0.5), so the reported estimates
# are point values from the gradient-stopped step rather than a strictly
# converged minimum; see vignette Errata for the implication and a
# self-consistency cross-check against `Simulated_ddmoremockdata2.txt`.

NA_NA_paracetamol <- function() {
  description <- paste(
    "Mechanistic OGTT model from the DDMORE Foundation Model Repository",
    "(DDMODEL00000228) that uses paracetamol as a gastric-emptying tracer to",
    "drive a coupled paracetamol PK + glucose + GLP-1 + GIP system. Fifteen",
    "compartments span paracetamol stomach / intestine / central / peripheral",
    "(with saturable first-pass loss), glucose stomach / duodenum / jejunum /",
    "ileum / central / peripheral, two effect compartments for glucose-on-",
    "production and insulin-on-elimination delays, a cumulative first-pass-",
    "loss tally for paracetamol, and indirect-response states for the",
    "incretin hormones GLP-1 and GIP. The gastric-emptying rate KS is",
    "modulated downwards by duodenal glucose via a Hill function (IGD50 /",
    "GAM) and gated by a logistic lag(T-T50) profile; glucose absorption",
    "from each small-intestine segment is Michaelis-Menten in segment amount",
    "(KMG, RAMAXD / RAMAXJ / RAMAXI). Plasma insulin (INS) is a time-varying",
    "regressor that enters the central glucose compartment through a one-",
    "compartment effect delay (KIE). Type 2 diabetes mellitus (T2DM) is",
    "encoded as a binary indicator switching the glucose baseline (GSSH /",
    "GSSD), glucose clearance (CLGH / CLGD), insulin-dependent glucose",
    "clearance (CLGIH / CLGID), glucose bioavailability (FPGH / FPGD), and",
    "the empirical glucose-on-production exponent (GPRG = -2.79 healthy, 0",
    "T2DM). Body weight (WT) scales the central glucose volume linearly",
    "(VG * WT / 70). Outputs are observed on the linear scale: paracetamol",
    "concentration plus baseline noise (Cc, uM), glucose concentration",
    "(Cglu, mM), GLP-1 concentration (CGLP1), and GIP concentration (CGIP).",
    "Source listing reports the FOCEI step terminated due to rounding",
    "errors (NSIG = 0.5); see vignette Errata for the implication on",
    "parameter-precision claims."
  )
  reference <- paste(
    "DDMORE Foundation Model Repository: DDMODEL00000228.",
    "No linked publication identified in the bundle (the `.rdf` carries no",
    "model-described-in-literature URI, and the bundle contains no",
    "Model_Accommodations.txt). Source $PROBLEM line `Benjamins estimation",
    "and covariance`; License Registered to Uppsala University; NONMEM 7.3",
    "FOCEI run dated 3 April 2015 (run126c) followed by an importance-",
    "sampling-only standard-error step (run126h).",
    sep = " "
  )
  vignette <- "NA_NA_paracetamol"
  ddmore_id <- "DDMODEL00000228"
  replicate_of <- NULL
  paper_specific_compartments <- c("effect_glu_prod", "effect_ins", "cumloss_apap")

  units <- list(
    time          = "min",
    dosing        = "mg",                  # paracetamol oral dose; glucose dose enters in g, see description
    concentration = "umol/L"               # paracetamol; secondary outputs in mmol/L (glucose) and pmol/L (GLP-1, GIP), see description
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight, time-fixed at baseline. Used as a linear scaling factor on the central glucose volume (VG_i = THETA(9) * WT / 70). Reference weight 70 kg.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The bundle's `Simulated_ddmoremockdata2.txt` records per-subject WT in the BW column at approximately 87.9 kg (single representative subject in the simulated dataset; the real-data fit used 16 subjects with WT not echoed in the bundle).",
      source_name        = "BW"
    ),
    T2DM = list(
      description        = "Type-2-diabetes-mellitus indicator (1 = T2DM patient, 0 = normal-glucose-tolerance control). Switches the glucose baseline (GSSH vs GSSD), glucose clearance (CLGH vs CLGD), insulin-dependent glucose clearance (CLGIH vs CLGID), glucose bioavailability into central (FPGH = 0.909 vs FPGD = 1 FIXED), and the empirical glucose-on-production exponent (GPRG = -2.79 healthy, 0 T2DM).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "T2DM = 0 (normal-glucose-tolerance control)",
      notes              = "Distinct from the existing `DIAB` canonical (which does not distinguish Type 1 vs Type 2) because this model is specifically a Type-2-versus-healthy stratification of OGTT response. Carried per subject (time-fixed in the bundle's simulated dataset).",
      source_name        = "T2DM"
    ),
    INS = list(
      description        = "Plasma insulin concentration as a time-varying regressor input. Drives the insulin-on-glucose-elimination effect compartment through a first-order delay (KIE); the conversion factor 1/6.945 turns the bundle's pmol/L scale into the uU/mL scale used inside the model.",
      units              = "pmol/L (raw); the model internally rescales by 1/6.945 to uU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Same canonical name as `Bizzotto_2016_glucose.R`'s plasma-insulin regressor. The model uses linear interpolation between dataset rows via the `linear(INS)` declaration in `model()`. The bundle's source column is `INSU`; rename `INSU -> INS` before passing the dataset to `rxSolve`.",
      source_name        = "INSU"
    ),
    INS_BL = list(
      description        = "Baseline (fasting) plasma insulin concentration, time-fixed per subject. Used to initialise the insulin-on-glucose-elimination effect compartment so the OGTT starts at the fasting state (A_0(10) = INS_BL / 6.945) and to set the typical-value baseline glucose-production rate GPRO via ISS = INS_BL / 6.945.",
      units              = "pmol/L (raw); the model rescales by 1/6.945 to uU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "New canonical column registered alongside this extraction (specific scope). Distinct from `INS` (time-varying regressor): `INS_BL` is the time-fixed baseline-state anchor used in initial conditions and the baseline-glucose-production calculation. Bundle source column is `BASI`; rename `BASI -> INS_BL` before passing to `rxSolve`.",
      source_name        = "BASI"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 16L,
    n_studies      = 3L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "Mixed normal-glucose-tolerance and T2DM adults pooled from three studies of an OGTT challenge. Paracetamol (1500 mg / 25 mg infusion over 5 min at t = 15 min in the bundle's simulated dataset) is co-administered as a gastric-emptying tracer alongside an oral 75 g glucose load. Subject-level T2DM status is carried in the T2DM column (0/1); the simulated representative subject ships as T2DM = 0 (healthy control).",
    dose_range     = "Paracetamol (1500 mg oral, infused as zero-order over ~5 min in the bundle's simulated dataset). Oral glucose load (~25 g into stomach compartment 5 in the bundle's simulated dataset, though clinical OGTT protocols typically use 75 g).",
    regions        = NA_character_,
    notes          = "Subject count and study count from the `Output_real_run126c.lst` listing's `TOT. NO. OF INDIVIDUALS: 16` line and the three STUDY levels (1, 2, 3) referenced by the `IF(STUDY.EQ.n)` switches on APAPBL and T50 in `Executable_run126h.mod`. Demographic detail (age, sex, region, race, exact study identities) is not recoverable from the on-disk bundle; the linked publication is not on disk to consult."
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural paracetamol PK. Parameter values are the
    # FINAL ESTIMATES from `Output_real_run126c.lst` after MINIMIZATION
    # TERMINATED (FOCEI, OBJ = -7913.528) and before the importance-sampling
    # standard-error step. Equation forms are from
    # `Executable_run126h.mod` $PK and $DES blocks (lines 24-37, 158-163).
    # Time unit: minutes (per inline `.mod` comments on KA, KW, etc.).
    # ---------------------------------------------------------------------
    lcl       <- log(0.344)          ; label("Paracetamol apparent clearance CL (L/min)")             # `.lst` FINAL TH 1 = 3.44E-01
    lvc       <- log(27.0)           ; label("Paracetamol central volume V1 (L)")                     # `.lst` FINAL TH 2 = 2.70E+01
    lq        <- log(0.675)          ; label("Paracetamol inter-compartmental clearance Q2 (L/min)")  # `.lst` FINAL TH 3 = 6.75E-01
    lvp       <- log(22.0)           ; label("Paracetamol peripheral volume V2 (L)")                  # `.lst` FINAL TH 4 = 2.20E+01
    lka       <- fixed(log(0.140))   ; label("Paracetamol absorption rate KA (1/min) - FIXED")        # `.lst` FINAL TH 5 = 1.40E-01 (FIX)
    lapapbl   <- log(6.54)           ; label("Paracetamol baseline-noise concentration APAPBL (uM)")  # `.lst` FINAL TH 6 = 6.54E+00
    lvmax     <- log(17.0)           ; label("Paracetamol intestinal saturable first-pass Vmax (mg/min equivalent)")  # `.lst` FINAL TH 29 = 1.70E+01
    lkm_apap  <- log(168)            ; label("Paracetamol intestinal saturable first-pass KM (mg)")   # `.lst` FINAL TH 30 = 1.68E+02

    # ---------------------------------------------------------------------
    # Glucose disposition + intestinal-absorption structural parameters.
    # Parameter values from `Output_real_run126c.lst` FINAL TH 7-23; FIXED
    # status follows the `.lst` $THETA boundary table (LB == INIT == UB)
    # and the `Executable_run126h.mod` `FIX` flags. Equations from
    # `Executable_run126h.mod` $PK lines 40-58 and $DES lines 146-149,
    # 165-173.
    # ---------------------------------------------------------------------
    lgss_h    <- log(5.27)           ; label("Glucose baseline (healthy) GSSH (mM)")                  # `.lst` FINAL TH 7 = 5.27E+00
    lgss_d    <- log(7.48)           ; label("Glucose baseline (T2DM) GSSD (mM)")                     # `.lst` FINAL TH 8 = 7.48E+00
    lvg       <- fixed(log(9.33))    ; label("Glucose central volume VG (L) - FIXED")                 # `.lst` FINAL TH 9 = 9.33E+00 (FIX)
    lclg_h    <- fixed(log(0.0894))  ; label("Glucose insulin-independent clearance (healthy) CLGH (L/min) - FIXED")  # `.lst` FINAL TH 10 = 8.94E-02 (FIX)
    lclgi_h   <- log(0.00663)        ; label("Glucose insulin-dependent clearance (healthy) CLGIH (L/min/(uU/mL))")    # `.lst` FINAL TH 11 = 6.63E-03 (bounded at upper 8.3E-03)
    lq_g      <- fixed(log(0.442))   ; label("Glucose inter-compartmental clearance Q (L/min) - FIXED")  # `.lst` FINAL TH 12 = 4.42E-01 (FIX)
    lvp_g     <- fixed(log(8.56))    ; label("Glucose peripheral volume VP (L) - FIXED")              # `.lst` FINAL TH 13 = 8.56E+00 (FIX)
    lkge1     <- fixed(log(0.0573))  ; label("Glucose-on-production effect rate KGE1 (1/min) - FIXED") # `.lst` FINAL TH 14 = 5.73E-02 (FIX)
    lkie      <- fixed(log(0.0213))  ; label("Insulin-on-elimination effect rate KIE (1/min) - FIXED") # `.lst` FINAL TH 15 = 2.13E-02 (FIX)
    lclg_d    <- fixed(log(0.0287))  ; label("Glucose insulin-independent clearance (T2DM) CLGD (L/min) - FIXED")  # `.lst` FINAL TH 16 = 2.87E-02 (FIX)
    lclgi_d   <- log(0.00547)        ; label("Glucose insulin-dependent clearance (T2DM) CLGID (L/min/(uU/mL))")    # `.lst` FINAL TH 17 = 5.47E-03 (bounded at upper 8.3E-03)
    fpg_h     <- 0.909               ; label("Glucose bioavailability into central (healthy) FPGH (fraction, bounded [0,1])")  # `.lst` FINAL TH 18 = 9.09E-01 (bounded at upper 1)
    fpg_d     <- fixed(1)            ; label("Glucose bioavailability into central (T2DM) FPGD (fraction) - FIXED at 1")  # `.lst` FINAL TH 19 = 1.00E+00 (FIX)
    lramax_d  <- log(0.576)          ; label("Glucose duodenal-segment max absorption rate RAMAXD (g/min)")  # `.lst` FINAL TH 20 = 5.76E-01
    lramax_j  <- log(2.06)           ; label("Glucose jejunal-segment max absorption rate RAMAXJ (g/min)")   # `.lst` FINAL TH 21 = 2.06E+00
    lramax_i  <- log(1.33)           ; label("Glucose ileal-segment max absorption rate RAMAXI (g/min)")     # `.lst` FINAL TH 22 = 1.33E+00
    lkm_g     <- log(6.32)           ; label("Glucose segment-absorption KM (g)")                            # `.lst` FINAL TH 23 = 6.32E+00

    # ---------------------------------------------------------------------
    # Gastric-emptying and gut-glucose-inhibition parameters.
    # ---------------------------------------------------------------------
    lkw       <- fixed(log(0.140))   ; label("Water gastric-emptying rate KW (1/min) - FIXED")        # `.lst` FINAL TH 24 = 1.40E-01 (FIX)
    glumax    <- fixed(1)            ; label("Maximum fractional inhibition of gastric emptying by duodenal glucose GLUMAX - FIXED at 1")  # `.lst` FINAL TH 25 = 1.00E+00 (FIX)
    ligd50    <- log(7.42)           ; label("Duodenal-glucose amount giving 50 percent gastric-emptying inhibition IGD50 (g)")  # `.lst` FINAL TH 26 = 7.42E+00
    lgam      <- log(14.0)           ; label("Hill exponent on duodenal-glucose-mediated gastric-emptying inhibition GAM (unitless)")  # `.lst` FINAL TH 27 = 1.40E+01
    lt50      <- log(20.0)           ; label("Mid-point of the logistic lag profile gating gastric emptying T50 (min, bounded 15-30)")  # `.lst` FINAL TH 28 = 2.00E+01

    # ---------------------------------------------------------------------
    # Incretin (GLP-1, GIP) indirect-response structural parameters.
    # KEGIP = LOG(2)/7 and KEGLP1 = LOG(2)/1.5 are computed as constants
    # inside `model()` rather than declared as `ini()` parameters (the
    # `.mod` hard-codes them as half-lives 7 min and 1.5 min respectively;
    # treating them as constants matches the source's intent).
    # ---------------------------------------------------------------------
    lemglp1_d <- log(0.0279)         ; label("GLP-1 indirect-response stimulation by duodenal glucose EMGLP1D (1/g)")  # `.lst` FINAL TH 31 = 2.79E-02
    lemglp1_i <- log(1.98)           ; label("GLP-1 indirect-response stimulation by ileal glucose EMGLP1I (1/g)")     # `.lst` FINAL TH 32 = 1.98E+00
    le50glp1  <- log(0.0525)         ; label("Ileal-glucose-amount giving 50 percent GLP-1 stimulation E50GLP1I (g)")  # `.lst` FINAL TH 33 = 5.25E-02
    lemgip    <- log(14.3)           ; label("GIP indirect-response stimulation by duodenal glucose EMGIP (1/g)")       # `.lst` FINAL TH 34 = 1.43E+01
    le50gip   <- log(2.57)           ; label("Duodenal-glucose-amount giving 50 percent GIP stimulation E50GIP (g)")  # `.lst` FINAL TH 35 = 2.57E+00
    lbglp1    <- log(17.1)           ; label("Baseline GLP-1 concentration BGLP1 (pmol/L-equivalent units)")                              # `.lst` FINAL TH 36 = 1.71E+01
    lbgip     <- log(8.03)           ; label("Baseline GIP concentration BGIP (pmol/L-equivalent units)")                                 # `.lst` FINAL TH 37 = 8.03E+00

    # ---------------------------------------------------------------------
    # Inter-individual variability. Source `$OMEGA` block structure
    # (`Executable_run126h.mod` lines 258-278):
    #   ETA1 -> CL          (var 0.109)
    #   ETA2 -> V1          (var 0.143)
    #   ETA3,4,5 -> APAPBL  (BLOCK(1) SAME, var 0.0338) - one slot per STUDY
    #   ETA6 -> GSSH        (var 0.00177)
    #   ETA7 -> GSSD        (var 0.00762)
    #   ETA8 -> CLGH        (var 0.348 FIX)
    #   ETA9 -> CLGIH       (var 0.226)
    #   ETA10 -> CLGD       (var 0.348 FIX)
    #   ETA11 -> CLGID      (var 0.209)
    #   ETA12 -> KW         (var 0.25 FIX)
    #   ETA13 -> IGD50      (var 0.114)
    #   ETA14,15,16 -> T50  (BLOCK(1) SAME, var 0.0538) - one slot per STUDY
    #   ETA17 -> BGLP1      (var 0.0649)
    #   ETA18 -> BGIP       (var 0.177)
    # The 3-slot `BLOCK(1) SAME` structures for APAPBL and T50 are a NONMEM
    # idiom for partitioning a single eta across mutually exclusive study
    # cohorts (each subject is in one and only one STUDY level). Since
    # each subject contributes to exactly one of the three slots, the
    # effective variability structure is a single eta per parameter with
    # the BLOCK(1) variance value (0.0338 for APAPBL, 0.0538 for T50).
    # The per-study eta partitioning carries no additional information
    # at the typical-value / per-subject level and is collapsed here.
    # ---------------------------------------------------------------------
    etalcl     ~ 0.109               # `.lst` FINAL OMEGA ETA1 (diagonal) = 1.09E-01
    etalvc     ~ 0.143               # `.lst` FINAL OMEGA ETA2 (diagonal) = 1.43E-01
    etalapapbl ~ 0.0338              # `.lst` FINAL OMEGA ETA3..5 (BLOCK(1) SAME) = 3.38E-02
    etalgss_h  ~ 0.00177             # `.lst` FINAL OMEGA ETA6 (diagonal) = 1.77E-03
    etalgss_d  ~ 0.00762             # `.lst` FINAL OMEGA ETA7 (diagonal) = 7.62E-03
    etalclg_h  ~ fixed(0.348)        # `.lst` FINAL OMEGA ETA8 (diagonal) = 3.48E-01 (FIX)
    etalclgi_h ~ 0.226               # `.lst` FINAL OMEGA ETA9 (diagonal) = 2.26E-01
    etalclg_d  ~ fixed(0.348)        # `.lst` FINAL OMEGA ETA10 (diagonal) = 3.48E-01 (FIX)
    etalclgi_d ~ 0.209               # `.lst` FINAL OMEGA ETA11 (diagonal) = 2.09E-01
    etalkw     ~ fixed(0.25)         # `.lst` FINAL OMEGA ETA12 (diagonal) = 2.50E-01 (FIX)
    etaligd50  ~ 0.114               # `.lst` FINAL OMEGA ETA13 (diagonal) = 1.14E-01
    etalt50    ~ 0.0538              # `.lst` FINAL OMEGA ETA14..16 (BLOCK(1) SAME) = 5.38E-02
    etalbglp1  ~ 0.0649              # `.lst` FINAL OMEGA ETA17 (diagonal) = 6.49E-02
    etalbgip   ~ 0.177               # `.lst` FINAL OMEGA ETA18 (diagonal) = 1.77E-01

    # ---------------------------------------------------------------------
    # Residual error. Source `$ERROR` block uses log-additive epsilon on
    # each output (`Y = LOG(F + small) + EPS(i)`), which is exactly
    # nlmixr2's `lnorm(expSd)` residual structure where `expSd` is the SD
    # on the log scale. SD values are sqrt of the `.lst` FINAL SIGMA
    # diagonals: sqrt(0.0185) = 0.136 (paracetamol); sqrt(0.0123) = 0.111
    # (glucose); sqrt(0.129) = 0.359 (GLP-1); sqrt(0.314) = 0.561 (GIP).
    # ---------------------------------------------------------------------
    expSd        <- sqrt(0.0185)     ; label("Log-scale residual SD for paracetamol concentration (uM)")  # `.lst` FINAL SIGMA EPS1 (variance) = 1.85E-02
    expSd_Cglu   <- sqrt(0.0123)     ; label("Log-scale residual SD for glucose concentration (mM)")      # `.lst` FINAL SIGMA EPS2 (variance) = 1.23E-02
    expSd_Cglp1  <- sqrt(0.129)      ; label("Log-scale residual SD for GLP-1 concentration")             # `.lst` FINAL SIGMA EPS3 (variance) = 1.29E-01
    expSd_Cgip   <- sqrt(0.314)      ; label("Log-scale residual SD for GIP concentration")               # `.lst` FINAL SIGMA EPS4 (variance) = 3.14E-01
  })

  model({
    # The plasma-insulin time-course INS is supplied as a data-driven
    # regressor at every observation/dosing row; the `linear(INS)`
    # declaration tells rxode2 to linearly interpolate INS between rows
    # rather than carrying the row-time value forward. See the Bizzotto
    # 2016 glucose-kinetics model for the precedent and the canonical-
    # covariate entry `INS`.
    linear(INS)

    # -------------------------------------------------------------------
    # Individual structural parameters (log-normal eta around typical
    # value where IIV is declared).
    # -------------------------------------------------------------------
    cl       <- exp(lcl       + etalcl)
    vc       <- exp(lvc       + etalvc)
    q        <- exp(lq)
    vp       <- exp(lvp)
    ka       <- exp(lka)
    apapbl   <- exp(lapapbl   + etalapapbl)
    vmax     <- exp(lvmax)
    km_apap  <- exp(lkm_apap)

    gss_h    <- exp(lgss_h    + etalgss_h)
    gss_d    <- exp(lgss_d    + etalgss_d)
    vg       <- exp(lvg) * WT / 70           # BW scaling per `.mod` $PK line 43 (VG = THETA(9)*BW/70)
    clg_h    <- exp(lclg_h    + etalclg_h)
    clgi_h   <- exp(lclgi_h   + etalclgi_h)
    q_g      <- exp(lq_g)
    vp_g     <- exp(lvp_g)
    kge1     <- exp(lkge1)
    kie      <- exp(lkie)
    clg_d    <- exp(lclg_d    + etalclg_d)
    clgi_d   <- exp(lclgi_d   + etalclgi_d)
    ramax_d  <- exp(lramax_d)
    ramax_j  <- exp(lramax_j)
    ramax_i  <- exp(lramax_i)
    km_g     <- exp(lkm_g)

    kw       <- exp(lkw       + etalkw)
    igd50    <- exp(ligd50    + etaligd50)
    gam      <- exp(lgam)
    t50      <- exp(lt50      + etalt50)

    emglp1_d <- exp(lemglp1_d)
    emglp1_i <- exp(lemglp1_i)
    e50glp1  <- exp(le50glp1)
    emgip    <- exp(lemgip)
    e50gip   <- exp(le50gip)
    bglp1    <- exp(lbglp1    + etalbglp1)
    bgip     <- exp(lbgip     + etalbgip)

    # Incretin half-life-derived elimination rates (constants from
    # `Executable_run126h.mod` $PK lines 105-106).
    kegip    <- log(2) / 7        # GIP   half-life 7 min
    keglp1   <- log(2) / 1.5      # GLP-1 half-life 1.5 min

    # -------------------------------------------------------------------
    # T2DM stratification of glucose-arm parameters
    # (`Executable_run126h.mod` $PK lines 61-70). T2DM = 1 selects the
    # diabetic baseline / clearance / bioavailability / production-
    # exponent; T2DM = 0 selects the healthy values.
    # -------------------------------------------------------------------
    gss      <- gss_h   * (1 - T2DM) + gss_d   * T2DM
    clg      <- clg_h   * (1 - T2DM) + clg_d   * T2DM
    clgi     <- clgi_h  * (1 - T2DM) + clgi_d  * T2DM
    fpg      <- fpg_h   * (1 - T2DM) + fpg_d   * T2DM
    gprg     <- -2.79   * (1 - T2DM) + 0       * T2DM      # `.mod` $PK lines 67-68

    # -------------------------------------------------------------------
    # Baseline-state derivations
    # (`Executable_run126h.mod` $PK lines 72, 74-79).
    # ISS = INS_BL / 6.945 is the baseline insulin in uU/mL.
    # GPRO is the baseline glucose production rate (g/min equivalent).
    # -------------------------------------------------------------------
    iss      <- INS_BL / 6.945
    k12g     <- q_g / vg
    k21g     <- q_g / vp_g
    kg       <- clg  / vg
    kgi      <- clgi / vg
    gpro     <- gss * (kg + kgi * iss) * vg * 180 / 1000   # `.mod` $PK line 79

    # Paracetamol micro-constants for the central / peripheral block.
    k12a     <- q / vc
    k21a     <- q / vp
    kel      <- cl / vc

    # -------------------------------------------------------------------
    # GI transit-time geometry (`Executable_run126h.mod` $PK lines 109-122).
    # TT = total small-intestine transit time (min); FLD / FLJ / FLI are
    # the fractional lengths of duodenum / jejunum / ileum.
    # -------------------------------------------------------------------
    tt       <- 240
    fld      <- 0.08
    flj      <- 0.37
    fli      <- 0.55
    dtt      <- tt * fld
    kdj      <- 1 / dtt
    jtt      <- tt * flj
    kji      <- 1 / jtt
    itt      <- tt * fli
    kic      <- 1 / itt   # ileum -> colon (terminal drain on intestinal glucose)

    # -------------------------------------------------------------------
    # Derived dynamics (`Executable_run126h.mod` $DES lines 142-159).
    # RAD / RAJ / RAI are Michaelis-Menten absorption rates from each
    # small-intestine segment into central glucose. GLUINHI is the Hill
    # function modulating gastric emptying by duodenal glucose. KS is the
    # current gastric-emptying rate after glucose inhibition; LAG is a
    # logistic gate that prevents gastric emptying before time T50. FPL is
    # the saturable first-pass paracetamol loss from intestine.
    # -------------------------------------------------------------------
    glurat   <- effect_glu_prod / gss
    ggpr     <- (glurat + 0.0001)^gprg
    gprt     <- gpro * ggpr

    rad      <- ramax_d * duodenum_glu / (km_g + duodenum_glu)
    raj      <- ramax_j * jejunum_glu / (km_g + jejunum_glu)
    rai      <- ramax_i * ileum_glu   / (km_g + ileum_glu)
    ra       <- rad + raj + rai

    duoglu   <- duodenum_glu * (duodenum_glu >= 0)
    gluinhi  <- glumax * (duoglu + 0.00001)^gam /
                (igd50^gam + (duoglu + 0.00001)^gam)
    ks       <- kw * (1 - gluinhi)
    ks       <- ks * (ks >= 0)

    lag      <- 1 / (1 + exp(-10 * (t - t50)))
    fpl      <- vmax * intestine_apap / (km_apap + intestine_apap)

    # -------------------------------------------------------------------
    # ODE system (`Executable_run126h.mod` $DES lines 160-177). State
    # naming maps to the source compartment list:
    #   stomach_apap     <- COMP=STOMACH      (paracetamol stomach)
    #   intestine_apap   <- COMP=INTESTI      (paracetamol intestine)
    #   central          <- COMP=CENTP        (paracetamol central)
    #   peripheral1      <- COMP=PERIPHE      (paracetamol peripheral)
    #   stomach_glu      <- COMP=SG           (glucose stomach)
    #   duodenum_glu     <- COMP=DG           (glucose duodenum)
    #   central_glu      <- COMP=CENTG        (glucose central)
    #   peripheral1_glu  <- COMP=PERIG        (glucose peripheral)
    #   effect_glu_prod  <- COMP=EFF_G        (effect for glucose-on-production)
    #   effect_ins       <- COMP=EFF_I        (effect for insulin-on-elimination)
    #   jejunum_glu      <- COMP=JEJUN        (glucose jejunum)
    #   ileum_glu        <- COMP=ILEUM        (glucose ileum)
    #   cumloss_apap     <- COMP=KUML         (cumulative paracetamol first-pass loss)
    #   glp1             <- COMP=GLP1         (GLP-1 indirect-response state)
    #   gip              <- COMP=GIP          (GIP  indirect-response state)
    # -------------------------------------------------------------------
    d/dt(stomach_apap)    <- -ks * stomach_apap * lag
    d/dt(intestine_apap)  <-  ks * stomach_apap * lag - ka * intestine_apap - fpl
    d/dt(central)         <-  ka * intestine_apap - (kel + k12a) * central + k21a * peripheral1
    d/dt(peripheral1)     <-  k12a * central - k21a * peripheral1

    d/dt(stomach_glu)     <- -ks * stomach_glu * lag
    d/dt(duodenum_glu)    <-  ks * stomach_glu * lag - rad - kdj * duodenum_glu
    d/dt(central_glu)     <-  ra * fpg + gprt + k21g * peripheral1_glu -
                              (k12g + kg) * central_glu - kgi * central_glu * effect_ins
    d/dt(peripheral1_glu) <-  k12g * central_glu - k21g * peripheral1_glu
    d/dt(effect_glu_prod) <-  kge1 * (central_glu / vg / 180 * 1000 - effect_glu_prod)
    d/dt(effect_ins)      <-  kie  * (INS / 6.945 - effect_ins)
    d/dt(jejunum_glu)     <-  kdj * duodenum_glu - kji * jejunum_glu - raj
    d/dt(ileum_glu)       <-  kji * jejunum_glu - rai
    d/dt(cumloss_apap)    <-  fpl
    d/dt(glp1)            <- -keglp1 * (glp1 - bglp1) +
                              bglp1 * keglp1 * (emglp1_d * duodenum_glu +
                                                emglp1_i * ileum_glu /
                                                (e50glp1 + ileum_glu))
    d/dt(gip)             <- -kegip  * (gip - bgip) +
                              bgip  * kegip *
                              emgip * duodenum_glu / (duodenum_glu + e50gip)

    # -------------------------------------------------------------------
    # Initial conditions (`Executable_run126h.mod` $PK lines 125-139).
    # The paracetamol GI tract / central / peripheral states all start at
    # zero (default); the glucose central / peripheral / effect-
    # compartment states start at their steady-state values to anchor the
    # OGTT at a fasting baseline.
    # -------------------------------------------------------------------
    central_glu(0)     <- gss * vg * 180 / 1000
    peripheral1_glu(0) <- (k12g / k21g) * gss * vg * 180 / 1000
    effect_glu_prod(0) <- gss
    effect_ins(0)      <- iss
    glp1(0)            <- bglp1
    gip(0)             <- bgip

    # -------------------------------------------------------------------
    # Observation variables. Paracetamol central amount converted to a
    # uM concentration with the 151.17 g/mol molecular weight and central
    # volume vc (in L); glucose central amount converted from g to mM
    # via 1000 / (180 * vg). GLP-1 / GIP states are reported directly.
    # The baseline-noise APAPBL is added to the paracetamol prediction
    # before the log-normal residual, matching the `.mod` $ERROR line
    # `IF(CMT.EQ.3) IPRED = LOG(CAC + APAPBL + 0.00001)`.
    # -------------------------------------------------------------------
    Cc      <- central     / 151.17 / vc * 1000 + apapbl
    Cglu    <- central_glu / vg * 1000 / 180
    Cglp1   <- glp1
    Cgip    <- gip

    Cc      ~ lnorm(expSd)
    Cglu    ~ lnorm(expSd_Cglu)
    Cglp1   ~ lnorm(expSd_Cglp1)
    Cgip    ~ lnorm(expSd_Cgip)
  })
}
