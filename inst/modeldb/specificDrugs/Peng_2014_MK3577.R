Peng_2014_MK3577 <- function() {
  description <- paste(
    "Semi-mechanistic PD model for the glucagon receptor antagonist MK-3577",
    "(Merck) in healthy male subjects undergoing a glucagon challenge (Peng",
    "2014). Glucose / glucagon / insulin homeostasis is described with",
    "constant endogenous baselines Gss, Iss, GNss; glucose has central +",
    "peripheral distribution plus an effect compartment for delayed",
    "regulation of glucose production. MK-3577 drives an inhibitory Imax",
    "model on the glucagon-stimulation of glucose production (Imax,MK =",
    "0.961, IC50,MK = 13.9 nM) and a stimulatory Emax model on glucagon",
    "secretion (Emax,MK = 0.788 FIX, EC50,MK = 575 nM FIX) for the",
    "prechallenge compensatory feedback. Sandostatin (octreotide) PK is",
    "modeled with literature CL / V (0.121 L/kg/h, 0.194 L/kg) and inhibits",
    "endogenous insulin (IC50,S2 = 0.921 ng/mL) and glucagon (IC50,S1 = 5.50",
    "ng/mL) secretion at a fixed Imax of 1. The MK-3577 PK layer is NOT",
    "modeled here because the absorption rate ka, apparent volume V/F, and",
    "molecular weight for the mg-to-nM conversion are not reported in the",
    "on-disk paper or its tables; users supply the MK-3577 plasma",
    "concentration as a time-varying covariate column CP_MK3577_NM (nM)",
    "per the standing operator decision (extract PD layer only; PK supplied",
    "externally). See vignette Assumptions and deviations for the gap.",
    sep = " "
  )
  reference <- paste(
    "Peng JZ, Denney WS, Musser BJ, Liu R, Tsai K, Fang L, Reitman ML,",
    "Troyer MD, Engel SS, Xu L, Stoch A, Stone JA, Kowalski KG (2014).",
    "A Semi-mechanistic Model for the Effects of a Novel Glucagon Receptor",
    "Antagonist on Glucagon and the Interaction Between Glucose, Glucagon,",
    "and Insulin Applied to Adaptive Phase II Design.",
    "The AAPS Journal 16(6):1259-1270.",
    "doi:10.1208/s12248-014-9648-x.",
    sep = " "
  )
  vignette <- "Peng_2014_MK3577"

  paper_specific_compartments <- c(
    "glucose_peripheral", "glucagon", "sandostatin"
  )

  units <- list(
    time          = "h",
    dosing        = paste(
      "Time is in hours; supply infusions via cmt + amt + dur (2 hour",
      "infusions). Glucagon challenge: cmt = glucagon, amt = 360 ng/kg,",
      "dur = 2 (3 ng/kg/min for 2 h converts to 360 ng/kg total). Sandostatin",
      "(octreotide): cmt = sandostatin, amt = 3600 ng/kg, dur = 2 (30 ng/kg/min",
      "for 2 h). Basal insulin: cmt = insulin, amt = 12000 uU/kg, dur = 2",
      "(0.1 mIU/kg/min = 100 uU/kg/min for 2 h). MK-3577 is NOT dosed:",
      "supply CP_MK3577_NM (nM) as a time-varying covariate column.",
      sep = " "
    ),
    concentration = paste(
      "Glucose Gc in mg/dL; insulin Ic in uU/mL; glucagon GNc in pg/mL.",
      "Sandostatin internal state is amount per kg in ng/kg with concentration",
      "Cs = sandostatin / vsand in ng/mL.",
      sep = " "
    )
  )

  covariateData <- list(
    CP_MK3577_NM = list(
      description        = paste(
        "Instantaneous MK-3577 plasma concentration in nM supplied as a",
        "time-varying regressor that drives the Imax inhibition of glucagon",
        "stimulation on glucose production (IC50,MK = 13.9 nM) and the Emax",
        "stimulation of glucagon secretion (EC50,MK = 575 nM FIX). Set to 0",
        "for placebo subjects and for time points before any MK-3577 dose.",
        sep = " "
      ),
      units              = "nM",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Supplied externally because the on-disk Peng 2014 PDF does not",
        "report the MK-3577 absorption rate ka, apparent volume V/F, or",
        "molecular weight needed to derive a mg-dose-to-nM-plasma profile",
        "internally. Users with access to a published or internal MK-3577",
        "PK model should simulate that PK separately and feed the resulting",
        "Cp profile in nM through this column. For demonstration the",
        "validation vignette uses a stylised one-compartment first-order",
        "absorption profile with placeholder ka / V / MW values that are",
        "documented in the vignette Errata section as illustrative only.",
        sep = " "
      ),
      source_name        = "CMK"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 36L,
    n_studies      = 1L,
    age_range      = "FIM cohort 35 +/- 7.1 years (Table I)",
    weight_range   = NA_character_,
    sex_female_pct = 0,
    disease_state  = paste(
      "Healthy male subjects undergoing a glucagon challenge study (Peng 2014",
      "Methods, FIM Study). Randomized, double-blind, placebo-controlled,",
      "crossover with balanced incomplete block design (Table II); 36 healthy",
      "males received single oral doses of MK-3577 1-900 mg either in the AM",
      "or PM. Body mass index 24.7 +/- 2.5 kg/m^2 (Table I).",
      sep = " "
    ),
    dose_range     = paste(
      "MK-3577 single oral dose 0 (placebo) / 1 / 3 / 10 / 20 / 30 / 40 /",
      "100 / 300 / 600 / 900 mg, AM or PM (Table II). Glucagon challenge:",
      "2-h IV infusions of glucagon (3 ng/kg/min), Sandostatin / octreotide",
      "(30 ng/kg/min), and basal insulin (0.1 mIU/kg/min) starting 3, 12, or",
      "24 h post MK-3577 dose. Sandostatin blocks endogenous insulin and",
      "glucagon secretion; basal insulin reduces excessive glycemia.",
      sep = " "
    ),
    regions        = NA_character_,
    notes          = paste(
      "Baseline demographics from Peng 2014 Table I (FIM cohort). PK / PD",
      "samples collected predose and up to 14 h postdose in parts II and",
      "III and up to 32 h postdose in part IV. Glucose, glucagon, and",
      "insulin were measured at each sampling time. The companion T2DM",
      "patient model with GPRG1 = 0, CLGI scaled to 11% of healthy, and a",
      "fold-increase theta in baseline glucose (typical theta fixed at 1,",
      "IIV 51% CV fixed) is provided as modellib('Peng_2014_MK3577_t2dm').",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Glucose disposition (Peng 2014 Table III, "Glucose" rows; reference
    # weight is implicit in the dL/kg/h and dL/kg units because the paper
    # normalizes all PK parameters by body weight). Footnote a of Table
    # III notes that the V_GC and CL terms are not independently
    # identifiable from observed glucose alone (only the ratios are); the
    # paper-reported point estimates are used here for forward simulation.
    # =====================================================================
    lclg      <- log(0.463)    ; label("Insulin-independent glucose clearance CL_G (dL/kg/h)")            # Peng 2014 Table III: CL_G = 0.463 (%RSE 36)
    lclgi     <- log(0.102)    ; label("Insulin-dependent glucose clearance CL_GI (dL/kg/h per uU/mL)")    # Peng 2014 Table III: CL_GI = 0.102 (%RSE 30)
    lqg       <- log(0.180)    ; label("Glucose intercompartmental clearance Q_G (dL/kg/h)")              # Peng 2014 Table III: Q_G = 0.180 (%RSE 14)
    lvgc      <- log(0.845)    ; label("Glucose central volume V_GC (dL/kg)")                             # Peng 2014 Table III: V_GC = 0.845 (%RSE 23)
    lvgp      <- log(0.301)    ; label("Glucose peripheral volume V_GP (dL/kg)")                          # Peng 2014 Table III: V_GP = 0.301 (%RSE 15)
    lgss      <- log(91.9)     ; label("Glucose steady-state concentration G_ss (mg/dL)")                 # Peng 2014 Table III: G_ss = 91.9 (%RSE 1.3)
    lke0      <- log(0.084)    ; label("Glucose effect-compartment rate constant K_GE (1/h)")             # Peng 2014 Table III: K_GE = 0.084 (%RSE 43)

    # Glucose-feedback shape exponents. GPRG1 is negative in the paper so
    # it is kept on the linear scale (paper reports -2.08 for GPRG1 and
    # 4.05 for GPRG3); GPRG3 is positive but the linear form is preserved
    # for parallel readability with GPRG1.
    gprg1     <- -2.08         ; label("Glucose negative-feedback exponent GPRG1 (unitless; negative)")    # Peng 2014 Table III: GPRG1 = -2.08 (%RSE 26)
    gprg3     <- 4.05          ; label("Glucagon stimulatory-effect exponent GPRG3 (unitless; positive)")  # Peng 2014 Table III: GPRG3 = 4.05 (%RSE 10)

    # =====================================================================
    # MK-3577 effects on glucose production and glucagon secretion (Peng
    # 2014 Table III rows under Glucose and Glucagon). E_max,MK and EC_50,MK
    # were fixed during estimation to values from initial runs (see Table
    # III footer "Some parameter estimates were fixed... E_max,MK and
    # EC_50,MK").
    # =====================================================================
    limax_mk_gp <- log(0.961)             ; label("MK-3577 Imax for inhibition of glucagon-driven glucose production (unitless)")  # Peng 2014 Table III: Imax,MK = 0.961 (%RSE 1.7)
    lic50_mk_gp <- log(13.9)              ; label("MK-3577 IC50 for inhibition of glucagon-driven glucose production (nM)")        # Peng 2014 Table III: IC50,MK = 13.9 (%RSE 14)
    lemax_mk_gn <- fixed(log(0.788))      ; label("MK-3577 Emax for stimulation of glucagon secretion (unitless; FIX)")            # Peng 2014 Table III: Emax,MK = 0.788 FIX
    lec50_mk_gn <- fixed(log(575))        ; label("MK-3577 EC50 for stimulation of glucagon secretion (nM; FIX)")                   # Peng 2014 Table III: EC50,MK = 575 FIX

    # =====================================================================
    # Glucagon disposition + Sandostatin inhibition (Peng 2014 Table III
    # "Glucagon" rows). The product GN_ss * CL_GN is the steady-state
    # endogenous glucagon secretion rate.
    # =====================================================================
    lclgn     <- fixed(log(3.19))         ; label("Glucagon clearance CL_GN (L/kg/h; FIX-IIV but mean estimated)")  # Peng 2014 Table III: CL_GN = 3.19 (%RSE 3.4; IIV 18.4% FIX)
    lvgn      <- log(1.39)                ; label("Glucagon volume V_GN (L/kg)")                                     # Peng 2014 Table III: V_GN = 1.39 (%RSE 7.7)
    lgnss     <- log(58.3)                ; label("Glucagon steady-state concentration GN_ss (pg/mL)")               # Peng 2014 Table III: GN_ss = 58.3 (%RSE 3.3)
    lic50_s1  <- log(5.50)                ; label("Sandostatin IC50 for inhibition of glucagon secretion (ng/mL)")    # Peng 2014 Table III: IC50,S1 = 5.50 (%RSE 10)

    # =====================================================================
    # Insulin disposition + Sandostatin inhibition (Peng 2014 Table III
    # "Insulin" rows). I_ss * CL_I is the steady-state endogenous insulin
    # secretion rate.
    # =====================================================================
    lcli      <- log(1.40)                ; label("Insulin clearance CL_I (L/kg/h)")                                  # Peng 2014 Table III: CL_I = 1.40 (%RSE 13)
    lvi       <- log(0.320)               ; label("Insulin volume V_I (L/kg)")                                        # Peng 2014 Table III: V_I = 0.320 (%RSE 33)
    liss      <- log(4.14)                ; label("Insulin steady-state concentration I_ss (uU/mL)")                  # Peng 2014 Table III: I_ss = 4.14 (%RSE 4.8)
    liprg     <- log(2.30)                ; label("Glucose stimulatory-effect exponent on insulin secretion IPRG (unitless)")  # Peng 2014 Table III: IPRG = 2.30 (%RSE 13)
    lic50_s2  <- log(0.921)               ; label("Sandostatin IC50 for inhibition of insulin secretion (ng/mL)")     # Peng 2014 Table III: IC50,S2 = 0.921 (%RSE 23)

    # =====================================================================
    # Sandostatin (octreotide) PK: one-compartment with CL / V from the
    # Peng 2014 Results paragraph "A one-compartment model was used for
    # Sandostatin pharmacokinetics (0.121 L/kg/h and 0.194 L/kg as
    # clearance and volume of distribution, respectively)" (literature
    # values, not estimated within the FIM cohort).
    # =====================================================================
    lclsand   <- fixed(log(0.121))        ; label("Sandostatin clearance CL_S (L/kg/h; literature)")                 # Peng 2014 Results: CL_S = 0.121
    lvsand    <- fixed(log(0.194))        ; label("Sandostatin volume V_S (L/kg; literature)")                       # Peng 2014 Results: V_S = 0.194

    # =====================================================================
    # Inter-individual variability (Peng 2014 Table III "IIV (%CV)" column).
    # The paper fixes the IIV on Gss, Vgc, Iss, Cli, Gnss, Clgn to lead-
    # compound values; only IC50,MK has an estimated IIV (77.7% CV). %CV
    # is converted to a log-scale variance via omega^2 = log(CV^2 + 1).
    # =====================================================================
    etalgss      ~ fixed(log(0.061^2 + 1))    # Gss     IIV 6.1% CV  FIX -> omega^2 = log(1.003721) = 0.003714
    etalvgc      ~ fixed(log(0.288^2 + 1))    # Vgc     IIV 28.8% CV FIX -> omega^2 = log(1.082944) = 0.07969
    etaliss      ~ fixed(log(0.333^2 + 1))    # Iss     IIV 33.3% CV FIX -> omega^2 = log(1.110889) = 0.10512
    etalcli      ~ fixed(log(0.263^2 + 1))    # CLI     IIV 26.3% CV FIX -> omega^2 = log(1.069169) = 0.06690
    etalgnss     ~ fixed(log(0.106^2 + 1))    # GNss    IIV 10.6% CV FIX -> omega^2 = log(1.011236) = 0.01117
    etalclgn     ~ fixed(log(0.184^2 + 1))    # CLGN    IIV 18.4% CV FIX -> omega^2 = log(1.033856) = 0.03330
    etalic50_mk_gp ~ log(0.777^2 + 1)         # IC50,MK IIV 77.7% CV     -> omega^2 = log(1.603729) = 0.47225

    # =====================================================================
    # Residual error (Peng 2014 Table III "Description" column for the
    # residual rows). Glucose and glucagon were log-transformed during
    # estimation so an additive residual on log-scale corresponds to a
    # proportional residual in linear space (RESG and RESGN reported as
    # %CV). Insulin was not log-transformed; the combined error model
    # collapsed to additive (RESI reported as uU/mL SD).
    # =====================================================================
    propSd_Gc  <- 0.0754              ; label("Proportional residual SD on glucose concentration (fraction)")   # Peng 2014 Table III: RESG = 7.54% (%RSE 8.1)
    addSd_Ic   <- 1.38                ; label("Additive residual SD on insulin concentration (uU/mL)")           # Peng 2014 Table III: RESI = 1.38 uU/mL (%RSE 7.4)
    propSd_GNc <- 0.303               ; label("Proportional residual SD on glucagon concentration (fraction)")   # Peng 2014 Table III: RESGN = 30.3% (%RSE 5.0)
  })

  model({
    # Linearly interpolate the user-supplied MK-3577 plasma concentration
    # between event rows so a coarse profile (e.g., one row per hour) still
    # drives the Imax / Emax terms smoothly. Declared at the top of model()
    # before any mu-referenced parameter expression to avoid parser
    # confusion.
    linear(CP_MK3577_NM)

    # ---------------------------------------------------------------------
    # Individual structural parameters. Only the IIV-bearing entries take
    # an eta term (the paper's Table III "IIV (%CV)" column lists FIX or
    # an estimated value for the entries below; all other parameters are
    # typical-value only).
    # ---------------------------------------------------------------------
    clg       <- exp(lclg)
    clgi      <- exp(lclgi)
    qg        <- exp(lqg)
    vgc       <- exp(lvgc + etalvgc)
    vgp       <- exp(lvgp)
    gss       <- exp(lgss + etalgss)
    ke0       <- exp(lke0)

    imax_mk_gp <- exp(limax_mk_gp)
    ic50_mk_gp <- exp(lic50_mk_gp + etalic50_mk_gp)
    emax_mk_gn <- exp(lemax_mk_gn)
    ec50_mk_gn <- exp(lec50_mk_gn)

    clgn      <- exp(lclgn + etalclgn)
    vgn       <- exp(lvgn)
    gnss      <- exp(lgnss + etalgnss)
    ic50_s1   <- exp(lic50_s1)

    cli       <- exp(lcli + etalcli)
    vi        <- exp(lvi)
    iss       <- exp(liss + etaliss)
    iprg      <- exp(liprg)
    ic50_s2   <- exp(lic50_s2)

    clsand    <- exp(lclsand)
    vsand     <- exp(lvsand)

    # ---------------------------------------------------------------------
    # Derived endogenous quantities. The MK-3577 drive into the inhibitory
    # and stimulatory Emax forms is the user-supplied time-varying
    # covariate CP_MK3577_NM (nM); 0 reproduces the placebo / no-drug
    # condition exactly (1 - 0 = 1 in the inhibition factor, 1 + 0 = 1 in
    # the stimulation factor).
    # ---------------------------------------------------------------------
    cmk       <- CP_MK3577_NM                   # nM, supplied as covariate
    # Sandostatin concentration: state in ng/kg, vsand in L/kg gives ng/L;
    # divide by 1000 to get ng/mL so it matches IC50_S1 / IC50_S2 units.
    cs        <- sandostatin / vsand / 1000     # ng/mL
    ge        <- effect                         # glucose effect-compartment concentration (mg/dL)
    cgc       <- glucose / vgc                  # mg/kg / (dL/kg) = mg/dL
    # Glucagon and insulin volumes are L/kg but the units 1 ng/L = 1 pg/mL
    # and 1 mU/L = 1 uU/mL make the unscaled division already give the
    # paper-reported concentration unit.
    cgn       <- glucagon / vgn                 # ng/kg / (L/kg) = ng/L = pg/mL
    ci        <- insulin / vi                   # mU/kg / (L/kg) = mU/L = uU/mL

    # GPROD0 is the steady-state glucose production rate (mg/kg/h) at
    # baseline glucose Gss and baseline insulin Iss, derived from total
    # insulin-independent + insulin-dependent clearance contributions
    # (Peng 2014 Eq. 3).
    gprod0    <- gss * (clg + clgi * iss)

    # GPROD modulates GPROD0 through (1) glucose-self negative feedback
    # via the effect compartment GE, (2) MK-3577 inhibition of the
    # glucagon-driven amplification, and (3) glucagon-driven
    # amplification with a positive exponent (Peng 2014 Eq. 2). The
    # leading 0.5 averages the two parallel contributions.
    drug_inh  <- 1 - imax_mk_gp * cmk / (ic50_mk_gp + cmk)
    gprod     <- 0.5 * gprod0 * ((ge / gss)^gprg1 +
                                 drug_inh * (cgn / gnss)^gprg3)

    # ---------------------------------------------------------------------
    # ODE system. Per-kg amounts in the canonical "amount per kg" basis;
    # concentrations are derived as state / V. Glucose central +
    # peripheral + effect form a three-compartment chain on the glucose
    # axis; insulin and glucagon each have a single central compartment
    # with secretion (driven by glucose / drug) and first-order
    # elimination; Sandostatin has a literature one-compartment PK with
    # first-order elimination.
    # ---------------------------------------------------------------------
    d/dt(glucose)            <- gprod +
                                (qg / vgp) * glucose_peripheral -
                                (qg / vgc + clg / vgc + (clgi / vgc) * ci) * glucose
    d/dt(glucose_peripheral) <- (qg / vgc) * glucose - (qg / vgp) * glucose_peripheral
    d/dt(effect)             <- ke0 * (cgc - effect)
    d/dt(insulin)            <- iss * cli * (cgc / gss)^iprg *
                                (1 - cs / (ic50_s2 + cs)) -
                                (cli / vi) * insulin
    d/dt(glucagon)           <- gnss * clgn *
                                (1 + emax_mk_gn * cmk / (ec50_mk_gn + cmk)) *
                                (1 - cs / (ic50_s1 + cs)) -
                                (clgn / vgn) * glucagon
    d/dt(sandostatin)        <- -(clsand / vsand) * sandostatin

    # Initial conditions at baseline (steady state pre-MK-3577).
    glucose(0)            <- gss * vgc
    glucose_peripheral(0) <- gss * vgp
    effect(0)             <- gss
    insulin(0)            <- iss * vi
    glucagon(0)           <- gnss * vgn
    sandostatin(0)        <- 0

    # Observation variables (multi-output): glucose mg/dL, insulin uU/mL,
    # glucagon pg/mL. Bare names match the Hong 2013 glucose-insulin
    # precedent (Gc, Ic); GNc is the paper-named glucagon output.
    Gc  <- cgc
    Ic  <- ci
    GNc <- cgn

    Gc  ~ prop(propSd_Gc)
    Ic  ~ add(addSd_Ic)
    GNc ~ prop(propSd_GNc)
  })
}
