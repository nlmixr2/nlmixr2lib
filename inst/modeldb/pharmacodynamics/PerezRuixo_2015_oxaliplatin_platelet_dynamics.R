PerezRuixo_2015_oxaliplatin_platelet_dynamics <- function() {
  description <- "Semi-mechanistic PD-only platelet-dynamics model for peritoneal carcinomatosis (PC) patients treated with cytoreductive surgery (CRS) alone or CRS followed by hyperthermic intraperitoneal oxaliplatin (HIO). Extends the Harker (2000) cytokinetic model with (i) a self-renewing megakaryocyte progenitor pool (prol) coupled to a (PLT0/PLT)^gamma feedback on the proliferation rate, (ii) a five-compartment Friberg-Bulitta transit chain for platelet aging (transit1 -> transit2 -> ... -> transit5) with a random-destruction rate ks on each transit compartment, (iii) a fixed release factor eta = 4000 platelets per maturing megakaryocyte, (iv) a transient post-surgery stimulation of prol proliferation with maximum effect SPmax that attenuates first-order with rate kp starting at surgery time t = 0, and (v) a power-function drug effect E_drug = alpha * CP_OXA_MGL^beta that inhibits prol proliferation. The paper does NOT develop an oxaliplatin PK model of its own; CP_OXA_MGL is obtained from empirical Bayes estimates of an upstream popPK model (Perez-Ruixo 2013 Cancer Chemother Pharmacol 71:693-704, cited as reference 15 -- not packaged in nlmixr2lib), and is consumed here as an exogenous time-varying covariate column (oxaliplatin plasma concentration in mg/L). Prior splenectomy multiplies ks by Phi = 0.475, prolonging platelet lifespan from 3.23 to 7.78 days. Age, body surface area, sex, total proteins and HIO carrier solution were tested and dropped. The paper additionally models platelet transfusions as a bolus of TRF0 = 255 x 10^9/L platelet-count-equivalent decaying first-order with rate kt = 0.104/h (t1/2 = 6.66 h); this transfusion sub-model is NOT retained in the packaged model file (only 5 / 80 patients received transfusions in the source cohort, kt has the highest RSE 53.8% in the model, and the paper itself notes 'this finding should be interpreted with caution given the limited number of patients who received transfusions'). Downstream users who need transfusion simulation can layer an exogenous decay on the observation."
  reference <- paste(
    "Perez-Ruixo C, Valenzuela B, Peris JE, Bretcha-Boix P,",
    "Escudero-Ortiz V, Farre-Alegre J, Perez-Ruixo JJ.",
    "Platelet Dynamics in Peritoneal Carcinomatosis Patients Treated with",
    "Cytoreductive Surgery and Hyperthermic Intraperitoneal Oxaliplatin.",
    "AAPS J. 2016;18(1):245-257 (published online 17 Nov 2015).",
    "doi:10.1208/s12248-015-9839-0.",
    "CP_OXA_MGL driver comes from the upstream oxaliplatin popPK (Perez-Ruixo C,",
    "Valenzuela B, Peris JE, et al. Cancer Chemother Pharmacol.",
    "2013;71:693-704) not packaged in nlmixr2lib -- supplied here as an",
    "exogenous time-varying covariate.",
    sep = " "
  )
  vignette <- "PerezRuixo_2015_oxaliplatin_platelet_dynamics"
  units <- list(
    time          = "h",
    dosing        = "not applicable (PD-only model; no drug dose enters the ODE system -- the oxaliplatin drug effect is driven by the exogenous CP_OXA_MGL covariate column)",
    concentration = "10^9 cells/L (observed circulating platelet count, PLT); oxaliplatin plasma concentration covariate CP_OXA_MGL in mg/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    prol     = list(analyte = "megakaryocyte progenitor", units = NA_character_, specimen = "blood cell", verified = FALSE),
    transit1 = list(analyte = "young platelet", units = NA_character_, specimen = "administration site", verified = FALSE),
    transit2 = list(analyte = "young platelet", units = NA_character_, specimen = "administration site", verified = FALSE),
    transit3 = list(analyte = "young platelet", units = NA_character_, specimen = "administration site", verified = FALSE),
    transit4 = list(analyte = "young platelet", units = NA_character_, specimen = "administration site", verified = FALSE),
    transit5 = list(analyte = "young platelet", units = NA_character_, specimen = "administration site", verified = FALSE)
  )

  covariateData <- list(
    CP_OXA_MGL = list(
      description        = "Oxaliplatin plasma concentration driving the platelet-dynamics drug effect E_drug = alpha * CP_OXA_MGL^beta on megakaryocyte-progenitor proliferation. Time-varying within subject.",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Required input. The Perez-Ruixo 2015 PD model uses the empirical",
        "Bayes individual PK parameter estimates from an upstream popPK paper",
        "(Perez-Ruixo 2013 Cancer Chemother Pharmacol 71:693-704, reference",
        "15 in the source; NOT packaged in nlmixr2lib at extraction time) to",
        "predict oxaliplatin plasma concentrations, and then treats the",
        "predicted CP_OXA_MGL(t) as an input to the PD model. In this port CP_OXA_MGL is",
        "supplied as a time-varying covariate column in the event table; the",
        "vignette shows two options for populating it: (i) an offline PK",
        "simulation with a simple bolus + monoexponential decline as a",
        "placeholder driver tuned to typical HIO Cmax ~ 1 mg/L and half-life",
        "~ 14 h (per Perez-Ruixo 2013 and Ferron 2008 Cancer Chemother",
        "Pharmacol 62:679-683), or (ii) values digitised from Figure 5 of",
        "the source paper. Set CP_OXA_MGL = 0 for pre-HIO or CRS-alone (cohort C)",
        "subjects so E_drug becomes 0.",
        sep = " "),
      source_name        = "Cp"
    ),
    PRIOR_SPLEN = list(
      description        = "Prior splenectomy indicator, 1 = prior splenectomy (spleen surgically removed), 0 = spleen intact. Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (spleen intact)",
      notes              = paste(
        "Perez-Ruixo 2015 Table I: multiplicative factor Phi = 0.475 on the",
        "random platelet destruction rate constant ks for splenectomized",
        "patients (ks_i = ks * Phi^PRIOR_SPLEN). Splenectomy prolongs",
        "platelet lifespan from 3.23 days (spleen intact) to 7.78 days",
        "(post-splenectomy) per the paper's Discussion. 30 of 80 patients",
        "(37.5%) had prior splenectomy driven by the ovarian / GI",
        "carcinomatosis debulking that includes splenectomy when tumour",
        "invades the splenic hilum.",
        sep = " "),
      source_name        = "splenectomy"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened in the forward-inclusion / backward-elimination covariate",
        "analysis but did NOT reach the p<0.05 forward-inclusion threshold",
        "on any PD model parameter and was excluded from the final model",
        "(Perez-Ruixo 2015 Results: 'The exploratory graphical analysis of",
        "the correlation between age, BSA, sex, total proteins, and carrier",
        "solution with PD model parameters did not suggest any statistically",
        "significant association').",
        sep = " ")
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened but not retained (see AGE above)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (see AGE above)."
    ),
    ALB = list(
      description = "Total proteins (serum)",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Total serum proteins (the paper's covariate label) is functionally",
        "analogous to serum albumin ALB in downstream models; screened but",
        "not retained (see AGE above).",
        sep = " ")
    ),
    FORM_HIO_CARRIER = list(
      description = "HIO carrier solution: 1 = icodextrin 4%, 0 = dextrose 5%",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Perez-Ruixo 2015 Methods and Results: the two HIO carrier",
        "solutions (isotonic 4% icodextrin in cohort A, isotonic 5% dextrose",
        "in cohort B) were tested on PD model parameters and did not reach",
        "significance. Its effect is captured entirely upstream in the",
        "oxaliplatin PK layer (higher molecular weight icodextrin reduces",
        "peritoneum-to-plasma absorption and lowers CP_OXA_MGL), which is not",
        "packaged in this PD-only extraction; downstream CP_OXA_MGL already reflects",
        "the carrier's effect.",
        sep = " ")
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 80L,
    n_studies      = 3L,
    n_observations = 1386L,
    disease_state  = "Peritoneal carcinomatosis (PC), heterogeneous origin (ovarian, colorectal, gastric, pseudomyxoma peritonei, mesothelioma).",
    cohorts        = paste(
      "Cohort A: 41 patients (51.2%) treated with CRS followed by HIO in",
      "isotonic 4% icodextrin. Cohort B: 21 patients (26.3%) CRS + HIO in",
      "isotonic 5% dextrose. Cohort C: 18 patients (22.5%) CRS alone (no",
      "HIO).",
      sep = " "),
    splenectomized_pct = 37.5,
    dose_range     = paste(
      "HIO administered at 200 mg/L initial peritoneal concentration for",
      "30 or 60 min (dose-finding simulations explored 100-500 mg/L for",
      "30 and 60 min). Platelet transfusions administered in 5/80 (6.25%)",
      "patients.",
      sep = " "),
    regions        = "Spain (Hospital Quiron Torrevieja, Alicante; grant GE-079/11 Conselleria de Sanidad de Comunidad Valenciana).",
    notes          = paste(
      "Age, BSA, sex, total proteins and HIO carrier solution were tested",
      "as covariates and did not reach significance. Only prior splenectomy",
      "was retained (on ks). The paper does not tabulate baseline",
      "demographic ranges (age, weight, sex distribution) in the trimmed",
      "text on disk; cross-cohort demographics are described in prior",
      "publications by the same group (Perez-Ruixo 2013 Cancer Chemother",
      "Pharmacol; Perez-Ruixo 2013 Clin Pharmacokinet; Valenzuela 2011",
      "AAPS J). Estimation was NONMEM 7.1.2 FOCE. The covariance step",
      "failed so parameter RSE values are from a 76/100 successful",
      "bootstrap (Table I).",
      sep = " ")
  )

  ini({
    # ------------------------------------------------------------------
    # SYSTEM (PD) PARAMETERS - Perez-Ruixo 2015 Table I 'System-related
    # parameters'. All log-transformed on the ini() scale except gamma
    # (unitless feedback exponent, kept on the natural scale).
    # ------------------------------------------------------------------

    lplt0 <- log(237)
    label("Baseline circulating platelet count PLT0 (10^9 cells/L)")           # Perez-Ruixo 2015 Table I: PLT0 = 237 x 10^9/L (RSE 3.72%)

    lktr <- log(7.09e-3)
    label("Platelet-progenitor maturation / first-order transit rate constant ktr (1/h). At steady state kprol equals ktr per the paper's balance for dProl/dt = 0.")   # Perez-Ruixo 2015 Table I: ktr = 7.09e-3 /h (RSE 35.5%)

    lks <- log(8.86e-3)
    label("First-order random destruction rate constant ks (1/h) for circulating platelets in patients with an intact spleen")                                   # Perez-Ruixo 2015 Table I: ks = 8.86e-3 /h (RSE 16.8%)

    lgamma <- log(0.621)
    label("Feedback exponent gamma on (PLT0 / PLT_total) driving compensatory megakaryocyte-progenitor proliferation (unitless)")                                # Perez-Ruixo 2015 Table I: gamma = 0.621 (RSE 45.9%)

    # ------------------------------------------------------------------
    # SURGICAL-STRESS PARAMETERS - Perez-Ruixo 2015 Table I 'Surgical-related parameters'.
    # kprol is multiplied by (1 + SPmax * exp(-kp*(t - t_surg))) starting
    # at the time of surgery (encoded as t_surg = 0 in this model file so
    # the surgery term reduces to 1 + SPmax * exp(-kp * t) for t >= 0).
    # ------------------------------------------------------------------

    lspmax <- log(2.09)
    label("Maximum surgical-stress stimulation of megakaryocyte-progenitor proliferation SPmax (unitless additive multiplier; kprol_i = ktr * (1 + SPmax * exp(-kp * t)))")   # Perez-Ruixo 2015 Table I: SPmax = 2.09 (RSE 49.6%)

    lkp <- log(3.43e-3)
    label("First-order attenuation rate constant kp (1/h) for the surgical-stress effect on kprol; corresponds to a half-life of 8.42 days per the paper's Results")           # Perez-Ruixo 2015 Table I: kp = 3.43e-3 /h (RSE 44.9%)

    # ------------------------------------------------------------------
    # DRUG (oxaliplatin) EFFECT PARAMETERS - Perez-Ruixo 2015 Table I 'Drug-related parameter'.
    # CP_OXA_MGL is the plasma oxaliplatin concentration in mg/L, supplied as an
    # exogenous time-varying covariate column (see covariateData).
    # ------------------------------------------------------------------

    lalpha <- log(0.881)
    label("Coefficient alpha of the drug-effect power function E_drug = alpha * CP_OXA_MGL^beta on kprol (L/mg)")                                        # Perez-Ruixo 2015 Table I: alpha = 0.881 L/mg (RSE 17.2%)

    lbeta <- log(2.63)
    label("Exponent beta of the drug-effect power function E_drug = alpha * CP_OXA_MGL^beta (unitless)")                                                  # Perez-Ruixo 2015 Table I: beta = 2.63 (RSE 17.4%)

    # ------------------------------------------------------------------
    # SPLENECTOMY COVARIATE - Perez-Ruixo 2015 Table I 'Splenectomy factor (Phi)'.
    # Multiplicative effect on ks: ks_i = ks * Phi^PRIOR_SPLEN.
    # PRIOR_SPLEN = 0 -> ks_i = ks; PRIOR_SPLEN = 1 -> ks_i = 0.475 * ks.
    # ------------------------------------------------------------------

    lphi_splen <- log(0.475)
    label("Log-multiplicative splenectomy factor Phi on ks; ks_i = ks * exp(lphi_splen * PRIOR_SPLEN) so splenectomized patients (PRIOR_SPLEN = 1) have ks_i = 0.475 * ks")                                                          # Perez-Ruixo 2015 Table I: Phi = 0.475 (RSE 26.2%)

    # ------------------------------------------------------------------
    # TRANSFUSION-RELATED PARAMETERS - Perez-Ruixo 2015 Table I 'Transfusion-related parameters'.
    # The source paper models platelet transfusions as a bolus of TRF0 = 255
    # x 10^9 cells/L (RSE 29.2%; IIV n/a) decaying first-order with rate kt =
    # 0.104 /h (RSE 53.8%; IIV 62.8% CV). This transfusion sub-model is
    # intentionally NOT retained in the packaged model file: only 5 / 80
    # (6.25%) patients received transfusions in the source cohort, kt has
    # the highest RSE in the model, and the paper itself flagged the
    # limited-data caveat ('This finding should be interpreted with caution
    # given the limited number of patients who received transfusions'). See
    # description and vignette Assumptions for the deviation rationale.
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY - Perez-Ruixo 2015 Table I 'Interindividual variability (CV%)'.
    # Reported as CV%; converted to log-normal omega^2 = log((CV/100)^2 + 1).
    # IIV on kt (62.8% CV) is dropped alongside the transfusion sub-model.
    # No IIV was reported on gamma, SPmax, beta or on the splenectomy factor.
    # No between-eta correlations reported (diagonal Omega).
    # ------------------------------------------------------------------

    etalplt0  ~ 0.1027     # PLT0 CV 32.9%: omega2 = log(0.329^2 + 1)
    etalktr   ~ 0.2005     # ktr CV 47.1%: omega2 = log(0.471^2 + 1)
    etalks    ~ 0.4947     # ks CV 80.0%: omega2 = log(0.80^2 + 1)
    etalkp    ~ 0.5497     # kp CV 85.6%: omega2 = log(0.856^2 + 1)
    etalalpha ~ 0.2806     # alpha CV 56.9%: omega2 = log(0.569^2 + 1)

    # ------------------------------------------------------------------
    # RESIDUAL ERROR - Perez-Ruixo 2015 Methods 'Statistical Model'.
    # 'Residual variability in platelet counts was evaluated using an
    # additive error model after natural logarithmic transformation of the
    # observations and model predictions'. Per the skill's Phase 4
    # verification checklist ('NONMEM "additive on log-scale" is
    # proportional in nlmixr2's linear space'), encode as prop(propSd).
    # Reported epsilon CV = 25.5% (Table I row 'Residual variability').
    # ------------------------------------------------------------------

    propSd <- 0.255
    label("Proportional residual error on observed platelet count (fraction; equivalent to the paper's 'additive on log-scale' 25.5% CV epsilon)")                                                                       # Perez-Ruixo 2015 Table I: residual epsilon = 25.5% (RSE 9.98%)
  })

  model({
    # ---- 1. Fixed structural constants ----
    # Each maturing megakaryocyte releases approximately eta = 4000
    # platelets to the circulation (Perez-Ruixo 2015 Methods,
    # Kaushansky & Roth Wintrobe's Clinical Hematology 11th ed. p. 605
    # cited as reference 33). N_PLT = 5 transit compartments is the
    # paper's chosen chain length (Methods: 'N PLT was arbitrarily set
    # to 5 in order to have a large enough number of compartments that
    # result in a smoothed distribution of platelet lifespan').
    eta_plt <- 4000

    # ---- 2. Individual parameters ----
    plt0  <- exp(lplt0  + etalplt0)
    ktr   <- exp(lktr   + etalktr)
    ks_ref <- exp(lks   + etalks)                       # reference (spleen-intact) ks
    ks    <- ks_ref * exp(lphi_splen * PRIOR_SPLEN)     # splenectomy-adjusted ks
    gamma <- exp(lgamma)
    spmax <- exp(lspmax)
    kp    <- exp(lkp    + etalkp)
    alpha <- exp(lalpha + etalalpha)
    beta  <- exp(lbeta)

    # ---- 3. Effective proliferation rate kprol_eff(t) ----
    # kprol_eff = ktr * (surgery stimulation) * (1 - drug inhibition)
    # At steady state before surgery kprol = ktr (Perez-Ruixo 2015 Methods
    # paragraph following Eq. 3). Surgery time t_surg is encoded here as
    # t = 0 so the surgery term reduces to (1 + SPmax * exp(-kp * t))
    # for t >= 0; the drug term is (1 - alpha * CP_OXA_MGL^beta).
    surgery_effect <- 1 + spmax * exp(-kp * t)
    e_drug         <- alpha * CP_OXA_MGL^beta
    kprol_eff      <- ktr * surgery_effect * (1 - e_drug)

    # ---- 4. Feedback on the megakaryocyte progenitor pool ----
    # Total observed platelet count is the sum of the five transit
    # compartments (the paper's transfusion sub-compartment is not
    # retained in this packaged model file; see description and vignette
    # Assumptions for the deviation).
    plt_total <- transit1 + transit2 + transit3 + transit4 + transit5
    feedback  <- (plt0 / plt_total)^gamma

    # ---- 5. Steady-state initial conditions ----
    # Friberg-Bulitta transit chain with per-compartment random-destruction
    # rate ks gives a geometric-decay steady state: transit_i = transit_1 *
    # r^(i-1) with r = ktr / (ktr + ks); sum(transit_i, i=1..5) = plt0. The
    # progenitor-pool initial value follows from dProl/dt = 0 at t = 0
    # (pre-surgery equilibrium, kprol = ktr, CP_OXA_MGL = 0, feedback = 1) and the
    # first transit-compartment SS balance eta * ktr * prol0 = (ktr + ks) *
    # transit1(0). The paper's Eq. 4 states an approximate transit_i =
    # plt0 / N_PLT initialization; the exact geometric-decay form is used
    # here so PLT_total(0) = plt0 holds without transient drift. See
    # vignette Assumptions and deviations.
    r_ss <- ktr / (ktr + ks)
    t1_ss <- plt0 * (1 - r_ss) / (1 - r_ss^5)
    transit1(0) <- t1_ss
    transit2(0) <- t1_ss * r_ss
    transit3(0) <- t1_ss * r_ss^2
    transit4(0) <- t1_ss * r_ss^3
    transit5(0) <- t1_ss * r_ss^4
    prol(0)     <- t1_ss * (ktr + ks) / (eta_plt * ktr)

    # ---- 6. ODE system ----
    # Progenitor pool (Perez-Ruixo 2015 Eq. 1):
    d/dt(prol) <- kprol_eff * feedback * prol - ktr * prol

    # Five-compartment Friberg-Bulitta transit chain with per-compartment
    # random destruction ks (Perez-Ruixo 2015 Eqs. 2-3):
    d/dt(transit1) <- eta_plt * ktr * prol     - ktr * transit1 - ks * transit1
    d/dt(transit2) <- ktr * transit1           - ktr * transit2 - ks * transit2
    d/dt(transit3) <- ktr * transit2           - ktr * transit3 - ks * transit3
    d/dt(transit4) <- ktr * transit3           - ktr * transit4 - ks * transit4
    d/dt(transit5) <- ktr * transit4           - ktr * transit5 - ks * transit5

    # ---- 7. Observation ----
    # Observed PLT is the sum of the five transit compartments (Perez-Ruixo
    # 2015 Eq. 6, transfusion contribution omitted in this port; see
    # description). Proportional
    # residual encodes the paper's 'additive on log-scale' error model
    # (Methods 'Statistical Model' and skill Phase 4 verification
    # checklist 'NONMEM additive on log-scale is proportional in
    # nlmixr2s linear space').
    PLT <- plt_total
    PLT ~ prop(propSd)
  })
}
