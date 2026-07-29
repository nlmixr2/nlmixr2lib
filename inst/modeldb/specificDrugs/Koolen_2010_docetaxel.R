Koolen_2010_docetaxel <- function() {
  description <- paste(
    "Five-compartment population PK model for intravenous and oral docetaxel",
    "with concomitant oral ritonavir in 36 adults with advanced cancer.",
    "Docetaxel: depot + single transit (Savic-style; ktr = 2/MAT) + three",
    "disposition compartments (central, peripheral1, peripheral2) with",
    "Bruno-style 3-compartment IV disposition (V_central, V_peripheral1,",
    "V_peripheral2; Q1 central-peripheral1, Q2 central-peripheral2).",
    "Clearance is parameterised via a well-stirred hepatic-extraction model",
    "(Wilkinson 1975) so that elimination is driven by intrinsic clearance",
    "CLi modulated by ritonavir plasma concentration via competitive",
    "inhibition (Ki = 0.028 ug/mL); CL = Q_hep * CLi / (CLi + Q_hep) with",
    "Q_hep fixed at 80 L/h. Hepatic bioavailability F_hep = Q_hep / (CLi +",
    "Q_hep) multiplies the depot -> transit transition to encode oral",
    "first-pass extraction. Gut bioavailability F_gut switches between F_doc",
    "= 0.19 (no ritonavir) and F_RTV = 0.39 (concomitant ritonavir, gated by",
    "the binary covariate CONMED_RTV). Polysorbate-80-driven micelle",
    "sequestration after IV docetaxel is encoded by a route-dependent central",
    "volume (V_central_iv = 9.8 L vs V_central_po = 44 L; gated by the",
    "binary per-dose-record covariate ROUTE_IV). Embedded one-compartment",
    "first-order-absorption ritonavir PK (depot_rtv + central_rtv) carries",
    "fixed typical-value parameters from Kappelhoff 2005 (CL/F = 10.5 L/h,",
    "V/F = 96.6 L, ka = 0.871 1/h, Tlag = 0.778 h) so that the ritonavir",
    "concentration that drives docetaxel CLi-inhibition is simulated within",
    "this single model file (modellib('Kappelhoff_2005_ritonavir') is the",
    "upstream source). Inter-individual variability on CLi,0, V_central_iv,",
    "V_central_po, MAT, F_depot (shared between F_doc and F_RTV), and Ki;",
    "correlated etas for CLi,0 ~ V_central_iv (rho = 0.446). Proportional",
    "residual error is encoded at the final-model typical value (32%); the",
    "source paper's separately-estimated higher 63% proportional RUV for the",
    "first 4 hours after oral administration is documented in the validation",
    "vignette's Assumptions and deviations section. Inter-occasion",
    "variability on CLi,0 (22%), MAT (52%), and F_RTV (44%) reported in the",
    "source is not propagated -- see vignette Assumptions and deviations."
  )
  reference <- paste(
    "Koolen SLW, Oostendorp RL, Beijnen JH, Schellens JHM, Huitema ADR.",
    "Population pharmacokinetics of intravenously and orally administered",
    "docetaxel with or without co-administration of ritonavir in patients",
    "with advanced cancer.",
    "Br J Clin Pharmacol. 2010;69(5):465-474.",
    "doi:10.1111/j.1365-2125.2010.03621.x.",
    "Embedded ritonavir PK parameters (CL/F, V/F, ka, Tlag) are inherited",
    "as fixed typical values from",
    "Kappelhoff BS, Huitema ADR, Crommentuyn KML, Mulder JW, Meenhorst PL,",
    "van Gorp ECM, Mairuhu ATA, Beijnen JH.",
    "Br J Clin Pharmacol. 2005;59(2):174-182.",
    "doi:10.1111/j.1365-2125.2004.02241.x;",
    "see modellib('Kappelhoff_2005_ritonavir')."
  )
  vignette <- "Koolen_2010_docetaxel"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CONMED_RTV = list(
      description        = paste(
        "1 = the docetaxel dose record is administered with concomitant oral",
        "ritonavir 100 mg (either simultaneously or 1 h prior to docetaxel),",
        "0 = docetaxel administered alone (no ritonavir). Time-fixed within",
        "an evaluated occasion; the binary indicator gates the switch",
        "between F_gut = F_doc (no RTV) and F_gut = F_RTV (concomitant RTV)",
        "inside model() and is independent of the per-time ritonavir plasma",
        "concentration that drives the CLi inhibition through the embedded",
        "ritonavir PK compartments."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (docetaxel alone, no ritonavir).",
      notes              = paste(
        "Koolen 2010 Methods Equation 4 (the non-concentration-dependent",
        "ritonavir effect on gut bioavailability) and Results Model",
        "development section: F1 = 0.19 (no RTV) vs F_RTV = 0.39 (RTV",
        "co-administered). For simulation set CONMED_RTV = 1 on every",
        "docetaxel dose record that is paired with a ritonavir dose on the",
        "same occasion, and 0 otherwise. Distinct from the concomitant",
        "ritonavir dose itself, which must be supplied as a separate dose",
        "event into the depot_rtv compartment so that the time-varying",
        "ritonavir plasma concentration (central_rtv / vc_rtv) can drive",
        "the competitive CLi inhibition."
      ),
      source_name        = "RTV (indicator in Equation 4)"
    ),
    ROUTE_IV = list(
      description        = paste(
        "1 = the dose record is administered intravenously (docetaxel IV",
        "infusion; dose into the central compartment), 0 = administered",
        "orally (docetaxel by mouth; dose into the depot compartment).",
        "Per-dose-record indicator that switches the docetaxel central",
        "volume of distribution between the IV value (V_central_iv = 9.8 L)",
        "and the oral value (V_central_po = 44.0 L). Time-varying within a",
        "subject when the same patient receives both IV and oral docetaxel",
        "on different occasions."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral docetaxel; reference V_central_po = 44.0 L).",
      notes              = paste(
        "Koolen 2010 Results 'Volume of distribution' section: 'Separate",
        "analyses of the volume of distribution of the central compartment",
        "for orally and intravenously administered docetaxel resulted in",
        "different estimates of 44 L (CV 17%) and 9.8 L (CV 10%),",
        "respectively.' The author attributes the IV-vs-oral discrepancy",
        "to polysorbate-80 micelle sequestration of docetaxel immediately",
        "after IV infusion (Discussion). For simulation set ROUTE_IV = 1 on",
        "IV dose records (cmt = central) and ROUTE_IV = 0 on oral dose",
        "records (cmt = depot); propagate the value forward through the",
        "observations associated with each dose. ROUTE_IV here flags an IV",
        "vs oral contrast (the canonical ROUTE_IV in nlmixr2lib is also",
        "used for IV vs SC, IV vs IP, etc. -- the reference category is",
        "paper-specific; see covariate-columns.md entry for the registered",
        "examples)."
      ),
      source_name        = "(NONMEM data column flagging IV vs oral administration; not explicitly named in the paper text)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 36L,
    n_studies       = 2L,
    n_observations  = 1025L,
    n_treatment_courses = 72L,
    age_range       = "31-73 years",
    age_median      = "54 years",
    sex_female_pct  = 44.4,
    disease_state   = paste(
      "Adults with histologically or cytologically confirmed advanced",
      "cancer for which no standard therapy was available (Koolen 2010",
      "Methods Patients section). Performance status <= 2 (WHO),",
      "life expectancy >= 3 months, adequate haematologic, hepatic",
      "(bilirubin < 20 umol/L; AST/ALT <= 1.5 x ULN or <= 3 x ULN with",
      "liver metastases), and renal (creatinine <= 160 umol/L or",
      "CrCl >= 50 mL/min) function. Patients on chronic H2-receptor",
      "antagonists or proton-pump inhibitors and patients with active",
      "bacterial / viral infections, CNS metastases, alcoholism, drug",
      "addiction, or pregnancy were excluded."
    ),
    dose_range      = paste(
      "Docetaxel intravenously 100 mg or 100 mg/m^2 (32 patients across",
      "the two cohorts) and orally 10 mg, 75 mg/m^2, or 100 mg with or",
      "without oral ritonavir 100 mg (administered simultaneously or 1 h",
      "prior to docetaxel). Per-arm patient counts per Koolen 2010",
      "Table 1: docetaxel IV 100 mg/m^2 = 15, docetaxel IV 100 mg = 17,",
      "oral docetaxel 75 mg/m^2 alone = 3, RTV 100 mg + oral docetaxel",
      "10 mg (1 h apart) = 6, RTV 100 mg + oral docetaxel 100 mg (1 h",
      "apart) = 15, RTV 100 mg + oral docetaxel 10 mg simultaneous = 5,",
      "RTV 100 mg + oral docetaxel 100 mg simultaneous = 11."
    ),
    regions         = "Netherlands (Netherlands Cancer Institute / Slotervaart Hospital, Amsterdam).",
    notes           = paste(
      "Pooled data from two phase I studies in the same institution (Koolen",
      "2010 Methods refs [4, 5]). 1025 docetaxel and 276 ritonavir plasma",
      "concentrations across 72 treatment courses. Standard dexamethasone +",
      "granisetron pre-treatment. Docetaxel assayed by validated LC-MS/MS",
      "(LLOQ 0.25 ng/mL, first study) or HPLC (LLOQ 10 ng/mL, second",
      "study); ritonavir by HPLC-UV (LLOQ 50 ng/mL). NONMEM VI ADVAN6",
      "TRANS1 with FO initialisation followed by FOCE final estimates.",
      "Patient-level ritonavir PK parameters were derived as MAP Bayesian",
      "posthoc estimates from the previously developed Kappelhoff 2005",
      "ritonavir popPK model (Koolen 2010 Methods Data analyses; reference",
      "[8] in the source). Sampling at predose, 0.25, 0.5, 0.75, 1, 1.25,",
      "1.5, 2, 3, 4, 6, 8, 10, 24, 36, and 48 h after administration."
    )
  )

  ini({
    # ============================================================
    # Docetaxel structural PK -- final model column of Koolen 2010
    # Table 2. Three-compartment disposition extended from Bruno
    # 1996 (Methods refs [17, 18]) with a depot + single transit
    # absorption chain (ktr = 2 / MAT, Methods Equation 2 and
    # Figure 1 caption) and a well-stirred hepatic-extraction model
    # (Wilkinson 1975, Methods Equations 1-3) that ties elimination
    # to intrinsic clearance modulated by ritonavir.
    # ============================================================
    lcli0   <- log(113)                ; label("Uninhibited intrinsic clearance, CLi,0 (L/h)")            # Table 2: CL i,0 = 113 L/h, CV 19%
    lvc_iv  <- log(9.8)                ; label("Docetaxel central volume after IV dosing, V_central_iv (L)") # Table 2: V2 (iv) = 9.8 L, CV 10%
    lvc_po  <- log(44.0)               ; label("Docetaxel central volume after oral dosing, V_central_po (L)") # Table 2: V2 (po) = 44.0 L, CV 17%
    lq      <- log(6.9)                ; label("Inter-compartmental clearance Q1 between central and peripheral1 (L/h)") # Table 2: Q1 = 6.9 L/h, CV 12%
    lvp     <- log(7.5)                ; label("Peripheral1 volume of distribution, V_peripheral1 (L)")   # Table 2: V3 = 7.5 L, CV 16%
    lq2     <- log(15.7)               ; label("Inter-compartmental clearance Q2 between central and peripheral2 (L/h)") # Table 2: Q2 = 15.7 L/h, CV 14%
    lvp2    <- log(376)                ; label("Peripheral2 volume of distribution, V_peripheral2 (L)")   # Table 2: V4 = 376 L, CV 15%
    lmat    <- log(1.3)                ; label("Docetaxel mean absorption time, MAT (h); ktr = 2 / MAT")  # Table 2: MAT = 1.3 h, CV 21%
    qhep    <- fixed(80)               ; label("Hepatic blood flow, Q_hep (L/h, FIXED)")                  # Table 2: Q = 80 L/h (fixed); Methods: 'The hepatic blood flow was fixed at 80 L h-1 ... estimation of this parameter resulted in a poor precision'.

    # Gut bioavailability of docetaxel. The canonical lfdepot is the
    # log-scale F when CONMED_RTV = 0; the additive log-scale
    # covariate effect e_rtv_fdepot shifts lfdepot to the with-RTV
    # value when CONMED_RTV = 1. The reported typical values give
    #   F_doc          = exp(lfdepot)                  = 0.19
    #   F_RTV          = exp(lfdepot + e_rtv_fdepot)   = 0.39
    # so e_rtv_fdepot = log(0.39 / 0.19) ~= 0.719. The shared 72%
    # CV on F (Koolen 2010 Table 2 reports the same 72% for the
    # F1 and F_RTV rows) is carried by a single eta on lfdepot.
    lfdepot      <- log(0.19)          ; label("Gut bioavailability of docetaxel, no ritonavir (F_doc, fraction)") # Table 2: F1 = 19%, CV 21%
    e_rtv_fdepot <- log(0.39 / 0.19)   ; label("Additive log-scale shift in F_gut for concomitant ritonavir (unitless)") # Table 2: F_RTV = 39%, CV 13% (paired with F1 = 19%)

    # Ritonavir-driven competitive inhibition of CLi:
    #   CLi(t) = CLi,0 / (1 + [RTV](t) / Ki)   (Methods Equation 3)
    # Ki is expressed in the same concentration unit as the embedded
    # ritonavir model output (mg/L = ug/mL).
    lki     <- log(0.028)              ; label("Ritonavir-on-docetaxel inhibition constant, Ki (ug/mL = mg/L)") # Table 2: K i = 0.028 ug/mL, CV 36%

    # ============================================================
    # Ritonavir structural PK -- inherited as FIXED typical-value
    # parameters from Kappelhoff 2005 (Koolen 2010 Methods ref [8]).
    # Koolen 2010 did not refit the ritonavir popPK; instead they
    # used MAP Bayesian posthoc estimates per subject from the
    # Kappelhoff 2005 NONMEM run. Holding the typical values here
    # makes the joint model self-contained for simulation; users
    # who need per-subject Bayesian posthocs should switch to
    # modellib('Kappelhoff_2005_ritonavir') and pass the resulting
    # ritonavir concentration profile in as a covariate.
    # ============================================================
    lcl_rtv   <- fixed(log(10.5))      ; label("Ritonavir apparent oral clearance, CL/F (L/h, FIXED at Kappelhoff 2005 Table 2 typical value)") # Kappelhoff 2005 Table 2: CL/F = 10.5 L/h
    lvc_rtv   <- fixed(log(96.6))      ; label("Ritonavir apparent volume of distribution, V/F (L, FIXED at Kappelhoff 2005 Table 2 typical value)") # Kappelhoff 2005 Table 2: V/F = 96.6 L
    lka_rtv   <- fixed(log(0.871))     ; label("Ritonavir first-order absorption rate, ka (1/h, FIXED at Kappelhoff 2005 Table 2 typical value)") # Kappelhoff 2005 Table 2: ka = 0.871 1/h
    ltlag_rtv <- fixed(log(0.778))     ; label("Ritonavir absorption lag time, Tlag (h, FIXED at Kappelhoff 2005 Table 2 typical value)") # Kappelhoff 2005 Table 2: Lag-time = 0.778 h

    # ============================================================
    # Inter-individual variability. Log-normal IIV on the natural
    # parameter ('exponential model' in Koolen 2010 Methods Data
    # analyses), so omega^2 = log(1 + CV^2) with CV reported as a
    # fraction from Table 2. The CLi,0 ~ V_central_iv correlation
    # 0.446 is inherited from the basic IV-only Bruno-style model
    # and is applied here to the (CLi,0, V_central_iv) eta block.
    # V_central_po, MAT, F_gut and Ki carry independent etas.
    # ============================================================
    etalcli0 + etalvc_iv ~ c(
      log(1 + 0.60^2),
      0.446 * sqrt(log(1 + 0.60^2) * log(1 + 0.45^2)),
      log(1 + 0.45^2)
    )
    # Table 2: IIV CLi,0 = 60% (RSE 26%), IIV V2 (iv) = 45% (RSE 27%),
    # 'CL ~ V2 Correlation' = 44.6%.

    etalvc_po ~ log(1 + 0.35^2)
    # Table 2: IIV V2 (po) = 35% (RSE 27%).

    etalmat   ~ log(1 + 0.87^2)
    # Table 2: IIV MAT = 87% (RSE 35%).

    # F_doc and F_RTV share a single inter-individual eta in the
    # source (Table 2 reports the same 72% CV for both rows),
    # implemented here by carrying one eta on lfdepot that is
    # applied regardless of the CONMED_RTV-driven typical-value
    # shift.
    etalfdepot ~ log(1 + 0.72^2)
    # Table 2: IIV F1 = IIV F_RTV = 72% (RSE 45%, shared estimate).

    etalki    ~ log(1 + 1.22^2)
    # Table 2: IIV K i = 122% (RSE 33%).

    # ============================================================
    # Residual error. Koolen 2010 Table 2 reports a 32%
    # proportional error in the final model and an additional 63%
    # proportional error for observations within the first 4 hours
    # after oral docetaxel administration. Only the structural
    # 32% term is propagated here; the elevated early-oral RUV
    # is documented in the validation vignette's Assumptions and
    # deviations section.
    # ============================================================
    propSd  <- 0.32                    ; label("Proportional residual error (fraction)")                 # Table 2: P = 32%, CV 14% (final model)
  })

  model({
    # ------------------------------------------------------------
    # Individual ritonavir PK parameters. Held at the Kappelhoff
    # 2005 typical values (no IIV in this model); used only to
    # compute the time-varying ritonavir plasma concentration that
    # drives docetaxel CLi inhibition.
    # ------------------------------------------------------------
    cl_rtv   <- exp(lcl_rtv)
    vc_rtv   <- exp(lvc_rtv)
    ka_rtv   <- exp(lka_rtv)
    tlag_rtv <- exp(ltlag_rtv)
    kel_rtv  <- cl_rtv / vc_rtv

    # ------------------------------------------------------------
    # Ritonavir plasma concentration (mg/L = ug/mL with central_rtv
    # in mg and vc_rtv in L). Drives the competitive CLi inhibition
    # via Koolen 2010 Methods Equation 3.
    # ------------------------------------------------------------
    crtv     <- central_rtv / vc_rtv

    # ------------------------------------------------------------
    # Individual docetaxel structural parameters.
    # ------------------------------------------------------------
    cli0     <- exp(lcli0 + etalcli0)
    vc_iv    <- exp(lvc_iv + etalvc_iv)
    vc_po    <- exp(lvc_po + etalvc_po)
    q        <- exp(lq)
    vp       <- exp(lvp)
    q2       <- exp(lq2)
    vp2      <- exp(lvp2)
    mat      <- exp(lmat + etalmat)
    ki       <- exp(lki + etalki)

    # Route-dependent central volume: V_central_iv when ROUTE_IV =
    # 1 (IV dose), V_central_po when ROUTE_IV = 0 (oral dose).
    # Koolen 2010 Discussion attributes the contrast to transient
    # polysorbate-80 micelle sequestration after IV infusion.
    vc       <- vc_iv * ROUTE_IV + vc_po * (1 - ROUTE_IV)

    # Gut bioavailability shifts between F_doc and F_RTV based on
    # the binary indicator CONMED_RTV. The shared eta is applied on
    # the log scale so that subjects with high F_doc also have high
    # F_RTV (Koolen 2010 Table 2 reports a single 72% CV that the
    # paper attributes to both F1 and F_RTV).
    lfdepot_eff <- lfdepot + e_rtv_fdepot * CONMED_RTV
    fdepot      <- exp(lfdepot_eff + etalfdepot)

    # ------------------------------------------------------------
    # Competitive ritonavir-on-docetaxel inhibition and the well-
    # stirred hepatic extraction model (Koolen 2010 Methods
    # Equations 1-3). Hepatic blood flow Q_hep is fixed at 80 L/h;
    # F_hep multiplies the depot -> transit transition to encode
    # the oral first-pass extraction, while the systemic elimination
    # term CL/Vc * central uses the same CL value because the
    # well-stirred model treats hepatic extraction and elimination
    # symmetrically.
    # ------------------------------------------------------------
    cli      <- cli0 / (1 + crtv / ki)
    fhep     <- qhep / (cli + qhep)
    cl       <- qhep * cli / (cli + qhep)

    # ------------------------------------------------------------
    # Absorption-chain rate constant (Koolen 2010 Figure 1 caption:
    # ktr = 2 / MAT for a single-transit Savic-style chain).
    # ------------------------------------------------------------
    ktr      <- 2 / mat

    # ------------------------------------------------------------
    # Micro-constants for the docetaxel disposition system.
    # ------------------------------------------------------------
    kel      <- cl / vc
    k12      <- q  / vc
    k21      <- q  / vp
    k13      <- q2 / vc
    k31      <- q2 / vp2

    # ------------------------------------------------------------
    # ODE system. Docetaxel: depot -> transit1 (with F_hep first-
    # pass extraction at the transition) -> central with parallel
    # exchanges to peripheral1 and peripheral2. Ritonavir: depot_rtv
    # -> central_rtv with a first-order absorption lag time.
    # Compartment numbering in Koolen 2010 Figure 1: A(1) = depot,
    # A(2) = transit, A(3) = central docetaxel, A(4) = peripheral1
    # docetaxel, A(5) = peripheral2 docetaxel, A(6) = depot_rtv,
    # A(7) = central_rtv.
    # ------------------------------------------------------------
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * fhep * depot - ktr * transit1
    d/dt(central)     <-  ktr * transit1 -
                          kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    d/dt(depot_rtv)   <- -ka_rtv  * depot_rtv
    d/dt(central_rtv) <-  ka_rtv  * depot_rtv - kel_rtv * central_rtv

    # ------------------------------------------------------------
    # Bioavailability and lag-time adjustments. F_gut multiplies
    # the oral docetaxel dose (Koolen 2010 Methods Equation 4 and
    # initial condition A(1) | t=0 = F_gut * Dose). Oral ritonavir
    # has no F multiplier (CL/F and V/F are apparent parameters in
    # the embedded Kappelhoff 2005 model, so F is absorbed into the
    # apparent terms).
    # ------------------------------------------------------------
    f(depot)        <- fdepot
    alag(depot_rtv) <- tlag_rtv

    # ------------------------------------------------------------
    # Observation. Docetaxel plasma concentration. Doses in mg and
    # volumes in L give central / vc in mg/L (= ug/mL). Koolen 2010
    # reports concentrations in ng/mL on the y-axis of Figures 3-4;
    # the validation vignette scales by 1000 when overlaying paper
    # figures.
    # ------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
