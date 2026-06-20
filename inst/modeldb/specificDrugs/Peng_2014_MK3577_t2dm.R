Peng_2014_MK3577_t2dm <- function() {
  description <- paste(
    "T2DM-patient adaptation of the Peng 2014 semi-mechanistic glucose /",
    "glucagon / insulin model with MK-3577 as the glucagon receptor",
    "antagonist (Peng 2014 Fig. 1b and 'CTS Method for T2DM Patients').",
    "Three structural changes vs. the healthy model: (1) GPRG1 = 0 (glucose",
    "self-regulation of GPROD is fully compromised in T2DM, matching",
    "Silber 2007); (2) CL_GI is scaled to 11% of the healthy value (lead-",
    "compound finding); (3) baseline glucose is elevated by a fold factor",
    "theta with typical value fixed at 1 (i.e., baseline FPG is 2 x healthy",
    "G_SS) and IIV fixed at 51% CV based on lead-compound data. The new",
    "T2DM baselines for insulin (I_SSP) and glucagon (GN_SSP) are derived",
    "from theta via Peng 2014 Eqs. 7-10 (insulin: I_SSP = I_SS * (1+theta)^IPRG;",
    "glucagon: closed-form rearrangement of Eqs. 9 and 10 conditioned on",
    "the T2DM-scaled CL_GI). The effect compartment for glucose negative",
    "feedback is omitted because GPRG1 = 0 makes its contribution",
    "vanish; the Sandostatin compartment is also omitted because the T2DM",
    "phase IIa CTS used in Peng 2014 had no glucagon challenge -- only",
    "endogenous homeostasis under multi-day MK-3577 dosing. The MK-3577 PK",
    "layer is NOT modeled here for the same reason as the healthy model",
    "(ka, V/F, MW not in the on-disk PDF); users supply the time-varying",
    "MK-3577 plasma concentration via CP_MK3577_NM (nM). See vignette",
    "Assumptions and deviations.",
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

  paper_specific_compartments <- c("glucose_peripheral", "glucagon")

  units <- list(
    time          = "h",
    dosing        = paste(
      "MK-3577 is NOT dosed: supply CP_MK3577_NM (nM) as a time-varying",
      "covariate column. No glucagon, Sandostatin, or basal insulin",
      "infusions are required for the T2DM CTS use case (the phase IIa",
      "study had no glucagon challenge).",
      sep = " "
    ),
    concentration = paste(
      "Glucose Gc in mg/dL; insulin Ic in uU/mL; glucagon GNc in pg/mL.",
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
        "internally. See modellib('Peng_2014_MK3577') for the healthy",
        "model and the shared vignette Assumptions and deviations section.",
        sep = " "
      ),
      source_name        = "CMK"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 118L,
    n_studies      = 1L,
    age_range      = "phase IIa cohort 54 +/- 10 years (Table I)",
    weight_range   = NA_character_,
    sex_female_pct = 47,
    disease_state  = paste(
      "Type 2 diabetes mellitus (Peng 2014 Methods, Phase IIa Study in T2DM",
      "Patients). Baseline FPG 152 +/- 35 mg/dL; A1c 7.6 +/- 0.8 %; 2-h",
      "PMG 223 +/- 63 mg/dL (Table I). Inclusion criterion in lead-compound",
      "internal studies used to derive the typical theta: baseline FPG",
      ">= 140 and <= 240 mg/dL.",
      sep = " "
    ),
    dose_range     = paste(
      "Phase IIa treatments: MK-3577 10 mg QD AM, 6 mg QD PM, 25 mg BID;",
      "metformin 1000 mg BID (active control); placebo. Prespecified interim",
      "analysis at N = 118 patients used the CTS output for adaptive dose",
      "adjustment (Tables V and VI).",
      sep = " "
    ),
    regions        = NA_character_,
    notes          = paste(
      "Baseline demographics from Peng 2014 Table I (Phase IIa cohort).",
      "The CTS used 82 patients per dose cohort and 1000 simulated trials.",
      "The actual baseline FPG was unavailable prior to the interim analysis",
      "due to blinding; the typical theta was fixed at 1 (twofold increase",
      "in baseline FPG vs. healthy) based on four internal T2DM studies",
      "applying the same FPG inclusion criterion, and IIV on theta was",
      "fixed at 51% CV from lead-compound data (Methods, p. 1264). The",
      "companion healthy-subjects model with glucagon challenge is",
      "modellib('Peng_2014_MK3577').",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Glucose disposition. Parameter values are inherited from the healthy
    # model (Peng 2014 Table III "Glucose" rows) -- Methods state "The
    # pharmacokinetic parameters for MK-3577 were kept the same between
    # healthy and T2DM subjects". CL_GI is scaled to 11% of the healthy
    # value (lead-compound finding); we encode this by storing CL_GI on
    # the healthy reference scale and applying the 0.11 scale factor in
    # the model() body. GPRG1 is forced to 0 in the model() body (not
    # estimated here) per the T2DM assumption.
    # =====================================================================
    lclg      <- log(0.463)    ; label("Insulin-independent glucose clearance CL_G (dL/kg/h)")             # Peng 2014 Table III: CL_G = 0.463 (%RSE 36)
    lclgi     <- log(0.102)    ; label("Healthy insulin-dependent glucose clearance CL_GI (dL/kg/h per uU/mL); scaled to 11% in model()")  # Peng 2014 Table III: CL_GI = 0.102 (%RSE 30); T2DM uses 0.11 * CL_GI
    lqg       <- log(0.180)    ; label("Glucose intercompartmental clearance Q_G (dL/kg/h)")               # Peng 2014 Table III: Q_G = 0.180 (%RSE 14)
    lvgc      <- log(0.845)    ; label("Glucose central volume V_GC (dL/kg)")                              # Peng 2014 Table III: V_GC = 0.845 (%RSE 23)
    lvgp      <- log(0.301)    ; label("Glucose peripheral volume V_GP (dL/kg)")                           # Peng 2014 Table III: V_GP = 0.301 (%RSE 15)
    lgss      <- log(91.9)     ; label("Healthy glucose steady-state G_ss (mg/dL); reference for T2DM")    # Peng 2014 Table III: G_ss = 91.9 (%RSE 1.3)

    gprg3     <- 4.05          ; label("Glucagon stimulatory-effect exponent GPRG3 (unitless; positive)")  # Peng 2014 Table III: GPRG3 = 4.05 (%RSE 10)

    # =====================================================================
    # MK-3577 drug-effect parameters (Peng 2014 Methods: "the model for
    # healthy subjects was adopted (including typical values, IIV, and
    # residual error) with modifications"). E_max,MK and EC_50,MK were
    # fixed during estimation in the healthy fit.
    # =====================================================================
    limax_mk_gp <- log(0.961)             ; label("MK-3577 Imax for inhibition of glucagon-driven glucose production (unitless)")  # Peng 2014 Table III: Imax,MK = 0.961 (%RSE 1.7)
    lic50_mk_gp <- log(13.9)              ; label("MK-3577 IC50 for inhibition of glucagon-driven glucose production (nM)")        # Peng 2014 Table III: IC50,MK = 13.9 (%RSE 14)
    lemax_mk_gn <- fixed(log(0.788))      ; label("MK-3577 Emax for stimulation of glucagon secretion (unitless; FIX)")            # Peng 2014 Table III: Emax,MK = 0.788 FIX
    lec50_mk_gn <- fixed(log(575))        ; label("MK-3577 EC50 for stimulation of glucagon secretion (nM; FIX)")                   # Peng 2014 Table III: EC50,MK = 575 FIX

    # =====================================================================
    # Glucagon disposition (no Sandostatin in T2DM CTS).
    # =====================================================================
    lclgn     <- fixed(log(3.19))         ; label("Glucagon clearance CL_GN (L/kg/h; FIX-IIV but mean estimated)")  # Peng 2014 Table III: CL_GN = 3.19 (%RSE 3.4; IIV 18.4% FIX)
    lvgn      <- log(1.39)                ; label("Glucagon volume V_GN (L/kg)")                                     # Peng 2014 Table III: V_GN = 1.39 (%RSE 7.7)
    lgnss     <- log(58.3)                ; label("Healthy glucagon steady-state GN_ss (pg/mL); reference for T2DM") # Peng 2014 Table III: GN_ss = 58.3 (%RSE 3.3)

    # =====================================================================
    # Insulin disposition (no Sandostatin in T2DM CTS).
    # =====================================================================
    lcli      <- log(1.40)                ; label("Insulin clearance CL_I (L/kg/h)")                                  # Peng 2014 Table III: CL_I = 1.40 (%RSE 13)
    lvi       <- log(0.320)               ; label("Insulin volume V_I (L/kg)")                                        # Peng 2014 Table III: V_I = 0.320 (%RSE 33)
    liss      <- log(4.14)                ; label("Healthy insulin steady-state I_ss (uU/mL); reference for T2DM")    # Peng 2014 Table III: I_ss = 4.14 (%RSE 4.8)
    liprg     <- log(2.30)                ; label("Glucose stimulatory-effect exponent on insulin secretion IPRG (unitless)")    # Peng 2014 Table III: IPRG = 2.30 (%RSE 13)

    # =====================================================================
    # Fold-increase in baseline glucose for T2DM (Peng 2014 Eq. 7). Typical
    # theta fixed at 1 (so baseline FPG is 2 x healthy G_ss); IIV fixed at
    # 51% CV from lead-compound data. Stored as ltheta = log(theta) so the
    # lognormal IIV is on the log-scale.
    # =====================================================================
    ltheta    <- fixed(log(1.0))          ; label("T2DM fold-increase in steady-state glucose theta (unitless; FIX)") # Peng 2014 Methods p. 1264: typical theta FIXED at 1 (twofold increase)

    # =====================================================================
    # IIV (Peng 2014 Table III + T2DM-specific note). Healthy IIVs on Gss,
    # Vgc, Iss, Cli, Gnss, Clgn are inherited from Table III. The T2DM
    # theta IIV is fixed at 51% CV (Peng 2014 Methods p. 1264: "The IIV
    # was fixed at 51% coefficient of variation (CV) based on the lead
    # compound data.").
    # =====================================================================
    etalgss       ~ fixed(log(0.061^2 + 1))    # Gss IIV 6.1% CV  FIX -> omega^2 = log(1.003721) = 0.003714
    etalvgc       ~ fixed(log(0.288^2 + 1))    # Vgc IIV 28.8% CV FIX -> omega^2 = log(1.082944) = 0.07969
    etaliss       ~ fixed(log(0.333^2 + 1))    # Iss IIV 33.3% CV FIX -> omega^2 = log(1.110889) = 0.10512
    etalcli       ~ fixed(log(0.263^2 + 1))    # CLI IIV 26.3% CV FIX -> omega^2 = log(1.069169) = 0.06690
    etalgnss      ~ fixed(log(0.106^2 + 1))    # GNss IIV 10.6% CV FIX -> omega^2 = log(1.011236) = 0.01117
    etalclgn      ~ fixed(log(0.184^2 + 1))    # CLGN IIV 18.4% CV FIX -> omega^2 = log(1.033856) = 0.03330
    etalic50_mk_gp ~ log(0.777^2 + 1)          # IC50,MK IIV 77.7% CV -> omega^2 = log(1.603729) = 0.47225
    etaltheta     ~ fixed(log(0.51^2 + 1))     # theta IIV 51% CV FIX -> omega^2 = log(1.2601) = 0.23120

    # =====================================================================
    # Residual error (Peng 2014: "model for healthy subjects was adopted
    # (including ... residual error)"). Same forms as healthy.
    # =====================================================================
    propSd_Gc  <- 0.0754              ; label("Proportional residual SD on glucose concentration (fraction)")   # Peng 2014 Table III: RESG = 7.54% (%RSE 8.1)
    addSd_Ic   <- 1.38                ; label("Additive residual SD on insulin concentration (uU/mL)")           # Peng 2014 Table III: RESI = 1.38 uU/mL (%RSE 7.4)
    propSd_GNc <- 0.303               ; label("Proportional residual SD on glucagon concentration (fraction)")   # Peng 2014 Table III: RESGN = 30.3% (%RSE 5.0)
  })

  model({
    # Linearly interpolate the user-supplied MK-3577 plasma concentration
    # between event rows so a coarse profile still drives the Imax / Emax
    # terms smoothly.
    linear(CP_MK3577_NM)

    # ---------------------------------------------------------------------
    # Individual structural parameters.
    # ---------------------------------------------------------------------
    clg       <- exp(lclg)
    clgi_h    <- exp(lclgi)            # healthy insulin-dependent clearance
    clgi_t    <- 0.11 * clgi_h          # T2DM scaled value (Peng 2014 Methods p. 1264)
    qg        <- exp(lqg)
    vgc       <- exp(lvgc + etalvgc)
    vgp       <- exp(lvgp)
    gss       <- exp(lgss + etalgss)    # healthy reference G_ss (subject-typical)

    imax_mk_gp <- exp(limax_mk_gp)
    ic50_mk_gp <- exp(lic50_mk_gp + etalic50_mk_gp)
    emax_mk_gn <- exp(lemax_mk_gn)
    ec50_mk_gn <- exp(lec50_mk_gn)

    clgn      <- exp(lclgn + etalclgn)
    vgn       <- exp(lvgn)
    gnss      <- exp(lgnss + etalgnss)  # healthy reference GN_ss (subject-typical)

    cli       <- exp(lcli + etalcli)
    vi        <- exp(lvi)
    iss       <- exp(liss + etaliss)    # healthy reference I_ss (subject-typical)
    iprg      <- exp(liprg)

    theta     <- exp(ltheta + etaltheta)

    # ---------------------------------------------------------------------
    # T2DM baselines derived from theta (Peng 2014 Eqs. 7, 8, 10). Eq. 7:
    # G_SSP = (1 + theta) * G_SS. Eq. 8: I_SSP = I_SS * (1 + theta)^IPRG.
    # Eq. 10: GN_SSP = GN_SS * ((2*(1+theta) * (CL_G + CL_GI*I_SSP) /
    # (CL_G + CL_GI*I_SS)) - 1)^(1/GPRG3). The CL_GI in both numerator
    # and denominator is the operative T2DM CL_GI -- the paper uses one
    # unsubscripted "CL_GI" symbol after stating the 11% scaling rule, so
    # both occurrences inherit the T2DM value. At theta = 0 the inner
    # ratio collapses to (CL_G + CL_GI_T2DM*I_SS)/(CL_G + CL_GI_T2DM*I_SS)
    # = 1, giving 2 - 1 = 1 and so GN_SSP = GN_SS (matches healthy
    # reference); for the typical theta = 1 the formula gives GN_SSP
    # approximately 84 pg/mL, consistent with hyperglucagonemia.
    # ---------------------------------------------------------------------
    gssp      <- (1 + theta) * gss
    issp      <- iss * (1 + theta)^iprg
    gnss_ratio_arg <- 2 * (1 + theta) * (clg + clgi_t * issp) /
                     (clg + clgi_t * iss) - 1
    gnssp     <- gnss * gnss_ratio_arg^(1 / gprg3)

    # ---------------------------------------------------------------------
    # Derived endogenous quantities. GPRG1 = 0 omits the GE term entirely;
    # the effect compartment is therefore not carried in this T2DM
    # variant (Peng 2014 Fig. 1b). The Sandostatin term is also removed
    # (no glucagon challenge in the phase IIa CTS).
    # ---------------------------------------------------------------------
    cmk       <- CP_MK3577_NM                    # nM, supplied as covariate
    cgc       <- glucose / vgc                   # mg/kg / (dL/kg) = mg/dL
    cgn       <- glucagon / vgn                  # ng/kg / (L/kg) = ng/L = pg/mL
    ci        <- insulin / vi                    # mU/kg / (L/kg) = mU/L = uU/mL

    # Reference GPROD_0 in the T2DM equation system (Peng 2014 Eq. 3 with
    # the operative T2DM CL_GI replacing the healthy CL_GI). The
    # combination of T2DM CL_GI and the GN_SSP-conditioned glucagon term
    # in Eq. 2 reproduces GPROD_0P = G_SSP * (CL_G + CL_GI_T2DM * I_SSP)
    # at baseline.
    gprod0    <- gss * (clg + clgi_t * iss)
    drug_inh  <- 1 - imax_mk_gp * cmk / (ic50_mk_gp + cmk)
    gprod     <- 0.5 * gprod0 * (1 + drug_inh * (cgn / gnss)^gprg3)

    # ---------------------------------------------------------------------
    # ODE system. The paper specifies that the same Eqs. 2-6 of the
    # healthy model are reused in T2DM, except that Eq. 6 substitutes
    # GN_SSP for GN_SS in the glucagon secretion rate. Eq. 5 keeps the
    # healthy reference G_SS and I_SS in (C_GC / G_SS)^IPRG and
    # I_SS * CL_I, because at baseline C_GC = G_SSP makes the secretion
    # numerically equal to I_SSP * CL_I via Eq. 8. CL_GI throughout the
    # glucose and GPROD_0 expressions uses the T2DM-scaled value. The
    # effect compartment is omitted because GPRG1 = 0.
    # ---------------------------------------------------------------------
    d/dt(glucose)            <- gprod +
                                (qg / vgp) * glucose_peripheral -
                                (qg / vgc + clg / vgc + (clgi_t / vgc) * ci) * glucose
    d/dt(glucose_peripheral) <- (qg / vgc) * glucose - (qg / vgp) * glucose_peripheral
    d/dt(insulin)            <- iss * cli * (cgc / gss)^iprg -
                                (cli / vi) * insulin
    d/dt(glucagon)           <- gnssp * clgn *
                                (1 + emax_mk_gn * cmk / (ec50_mk_gn + cmk)) -
                                (clgn / vgn) * glucagon

    # Initial conditions at T2DM steady state (derived from theta).
    glucose(0)            <- gssp * vgc
    glucose_peripheral(0) <- gssp * vgp
    insulin(0)            <- issp * vi
    glucagon(0)           <- gnssp * vgn

    # Observation variables: glucose mg/dL, insulin uU/mL, glucagon pg/mL.
    Gc  <- cgc
    Ic  <- ci
    GNc <- cgn

    Gc  ~ prop(propSd_Gc)
    Ic  ~ add(addSd_Ic)
    GNc ~ prop(propSd_GNc)
  })
}
