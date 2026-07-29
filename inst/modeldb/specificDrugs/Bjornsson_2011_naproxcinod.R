Bjornsson_2011_naproxcinod <- function() {
  description <- paste(
    "Joint population PK / pain intensity (PI) / informative-dropout model",
    "for naproxen following oral administration of naproxcinod (a naproxen",
    "nitrate ester prodrug), naproxen, or placebo after wisdom-tooth",
    "extraction (Bjornsson 2011, 242 patients with moderate-to-severe",
    "post-surgical dental pain).",
    "PK: one-compartment disposition of unbound naproxen with parallel",
    "Savic transit-compartment absorption chains for naproxcinod",
    "(MTT 1.77 h, NN 3.58) and naproxen (MTT 0.500 h, NN 4.23) feeding a",
    "shared central compartment. Total naproxen is computed from unbound",
    "via a saturable albumin-binding equation Ctot = Cu + Bmax * Cu /",
    "(Km + Cu) (Bmax = 643 umol/L, Km = 0.549 umol/L). Relative",
    "bioavailability of naproxen via naproxcinod vs naproxen is 59.7%.",
    "PD: pain intensity on a 100-mm visual analogue scale modeled as",
    "PI(t) = PI_baseline * (1 - placebo(t)) * (1 - drug(t)),",
    "where placebo(t) = Pmax * (1 - exp(-kpl * t)) (Pmax 20.2%,",
    "kpl 0.237 /h; additive IIV on Pmax allows individual PI to either",
    "decrease or increase from baseline) and drug(t) is a sigmoid Emax",
    "function of unbound naproxen with Emax fixed at 1, EC50 0.135 umol/L,",
    "and Hill exponent 1.61. TTE: rescue-medication request modeled as a",
    "Weibull hazard (lambda 0.00999, alpha 0.729) with log-linear",
    "covariate effects of PI(t) and (PI_baseline - 55) on the slope of",
    "PI(t); the hazard is set to zero for t < 1.5 h to reflect the",
    "protocol's rescue-medication abstention window."
  )
  reference <- paste(
    "Bjornsson MA, Simonsson USH.",
    "Modelling of pain intensity and informative dropout in a dental",
    "pain model after naproxcinod, naproxen and placebo administration.",
    "Br J Clin Pharmacol. 2011;71(6):899-906.",
    "doi:10.1111/j.1365-2125.2011.03924.x.",
    sep = " "
  )
  vignette <- "Bjornsson_2011_naproxcinod"
  paper_specific_compartments <- c("depot_naproxcinod", "depot_naproxen", "cumhaz")

  units <- list(
    time          = "h",
    dosing        = paste(
      "umol of the parent compound dosed (naproxcinod or naproxen).",
      "Convert mg to umol using MW(naproxcinod) = 317.30 g/mol or",
      "MW(naproxen) = 230.26 g/mol. Naproxcinod doses must target",
      "compartment depot_naproxcinod and naproxen doses must target",
      "compartment depot_naproxen so that the matching transit chain",
      "and bioavailability are applied."
    ),
    concentration = "umol/L (both total and unbound naproxen)"
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age at enrolment (years).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline demographic captured in the trial (range 19-38 years, Bjornsson 2011 Table 1) but not retained in the final PK / PI / TTE model.",
      source_name        = "AGE"
    ),
    BMI = list(
      description        = "Body mass index (kg/m^2).",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline demographic captured in the trial (range 18-31 kg/m^2, Bjornsson 2011 Table 1) but not retained in the final PK / PI / TTE model. Body weight itself was not transcribed in the publication.",
      source_name        = "BMI"
    ),
    SEXF = list(
      description        = "Sex indicator: 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Trial cohort was 48% male / 52% female (Bjornsson 2011 Table 1). Sex was not retained in the final PK / PI / TTE model.",
      source_name        = "SEX"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 242L,
    n_studies      = 1L,
    age_range      = "19-38 years",
    age_median     = "~25 years (per-arm means in Bjornsson 2011 Table 1 ranged 24.0-25.6 years)",
    weight_range   = "not transcribed in the publication (BMI 18-31 kg/m^2 reported in Table 1)",
    sex_female_pct = 52,
    race_ethnicity = NULL,
    disease_state  = "Healthy adults undergoing surgical removal of mandibular wisdom teeth under local anaesthesia, with moderate-to-severe post-surgical dental pain (entry criterion PI >= 40 mm on 100-mm VAS within 6 h after local anaesthetic).",
    dose_range     = paste(
      "Single oral dose: naproxcinod 375 mg (1182 umol), 750 mg (2364 umol),",
      "1500 mg (4728 umol), or 2250 mg (7092 umol); naproxen 500 mg (2172 umol);",
      "or placebo (0 umol). Rescue medication (ibuprofen 400 mg) was available",
      "from 1.5 h onward."
    ),
    regions        = "United Kingdom (Eastman International Centre for Excellence in Dentistry, London)",
    notes          = paste(
      "PI measured on a 100-mm visual analogue scale (VAS) immediately before",
      "drug administration (baseline) and at 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6,",
      "7, and 8 h after drug administration (or until rescue-medication request,",
      "after which no further PI measurements were taken). Total naproxen",
      "plasma concentrations were measured in ~one-third of patients (15, 12,",
      "18, 15, 16 patients in the naproxcinod 375 / 750 / 1500 / 2250 mg and",
      "naproxen 500 mg arms respectively) at the same time points; unbound",
      "naproxen was measured by ultrafiltration at 1, 3, and 8 h post-dose.",
      "LOQ: 0.5 umol/L (total) and 5 nmol/L (unbound). Three subjects who took",
      "rescue medication before 1.5 h were excluded from the dropout analysis.",
      "Sponsor: AstraZeneca. NONMEM VI with Laplace estimation; PK and PD",
      "fit sequentially (PK final estimates fixed in the PD fit) and then",
      "jointly refit. Source: Bjornsson 2011 Methods 'Study design',",
      "'Pharmacokinetic model', and Table 1."
    )
  )

  ini({
    # ============================================================
    # Pharmacokinetic model (Bjornsson 2011 Table 2 final estimates).
    # CL_u/F and V_u/F are apparent UNBOUND clearance and volume.
    # Concentration of unbound naproxen in central is in umol/L.
    # ============================================================
    lcl <- log(515)
    label("Apparent unbound clearance CL_u/F (L/h)")
    # Table 2 CL_u/F = 515 L/h (RSE 12.1%); IIV 25% (RSE 37%).

    lvc <- log(4290)
    label("Apparent unbound volume of distribution V_u/F (L)")
    # Table 2 V_u/F = 4290 L (RSE 13.6%); IIV 44% (RSE 29%).

    # Savic transit-compartment absorption: parameters per formulation.
    # The depot drains at the canonical Savic rate ktr = (NN + 1) / MTT.
    lmtt_naproxcinod <- log(1.77)
    label("Mean transit time for naproxcinod absorption (h)")
    # Table 2 MTT_naproxcinod = 1.77 h (RSE 10.8%); IIV 58% (RSE 24%).

    lnn_naproxcinod <- log(3.58)
    label("Log of number of transit compartments for naproxcinod absorption (unitless)")
    # Table 2 NN_naproxcinod = 3.58 (RSE 9.9%); IIV 58% (RSE 26%).
    # Correlated with MTT_naproxcinod (rho = -52%).

    lmtt_naproxen <- log(0.500)
    label("Mean transit time for naproxen absorption (h)")
    # Table 2 MTT_naproxen = 0.500 h (RSE 23.8%); IIV 100% (RSE 60%).

    lnn_naproxen <- log(4.23)
    label("Log of number of transit compartments for naproxen absorption (unitless)")
    # Table 2 NN_naproxen = 4.23 (RSE 24.8%); IIV 64% (RSE 68%).

    # Saturable plasma-protein binding model (paper Methods 'Pharmacokinetic
    # model'): Ctot = Cu + Bmax * Cu / (Km + Cu). Linear PK is assumed for
    # the unbound state; total concentration is derived from unbound for the
    # residual-error model. Bmax ~643 umol/L corresponds to one binding site
    # per albumin molecule (albumin ~600-700 umol/L in plasma).
    lbmax <- log(643)
    label("Maximum saturable binding of naproxen to plasma proteins, Bmax (umol/L)")
    # Table 2 Bmax = 643 umol/L (RSE 7.1%); IIV 17% (RSE 44%).

    lkm <- log(0.549)
    label("Unbound naproxen concentration at half-maximum protein binding, Km (umol/L)")
    # Table 2 Km = 0.549 umol/L (RSE 10.2%); no IIV.

    # Relative bioavailability of naproxen via naproxcinod vs direct naproxen.
    # F = F_rel is applied as the bio argument to transit() for the
    # naproxcinod depot; the naproxen depot gets bio = 1 (reference).
    lfdepot_naproxcinod <- log(0.597)
    label("Relative bioavailability of naproxen from naproxcinod vs naproxen (fraction)")
    # Table 2 F_rel = 59.7% (RSE 14.6%); no IIV.

    # PK IIV. Lognormal: omega^2 = log(1 + CV^2) per the convention noted
    # in references/parameter-names.md (Section "IIV").
    etalcl ~ log(1 + 0.25^2)
    # Table 2 IIV CL_u/F = 25% CV (RSE 37%).

    etalvc ~ log(1 + 0.44^2)
    # Table 2 IIV V_u/F = 44% CV (RSE 29%).

    etalbmax ~ log(1 + 0.17^2)
    # Table 2 IIV Bmax = 17% CV (RSE 44%).

    # Correlated IIV between MTT_naproxcinod and NN_naproxcinod (rho = -52%).
    # Off-diagonal covariance = rho * sqrt(var_mtt * var_nn).
    etalmtt_naproxcinod + etalnn_naproxcinod ~ c(
      log(1 + 0.58^2),
      -0.52 * sqrt(log(1 + 0.58^2) * log(1 + 0.58^2)),
      log(1 + 0.58^2)
    )
    # Table 2 IIV MTT_NC = 58% (RSE 24), IIV NN_NC = 58% (RSE 26),
    # Corr(MTT_NC, NN_NC) = -52% (RSE 38).

    etalmtt_naproxen ~ log(1 + 1.00^2)
    # Table 2 IIV MTT_naproxen = 100% CV (RSE 60%).

    etalnn_naproxen ~ log(1 + 0.64^2)
    # Table 2 IIV NN_naproxen = 64% CV (RSE 68%).

    # Residual error: total naproxen has combined additive + proportional;
    # unbound naproxen has proportional only (paper Methods).
    addSd <- 6.19
    label("Additive residual SD for total naproxen (umol/L)")
    # Table 2 sigma_T,add = 6.19 umol/L (RSE 22.3%).

    propSd <- 0.0843
    label("Proportional residual SD for total naproxen (fraction)")
    # Table 2 sigma_T,prop = 8.43% (RSE 8.0%).

    propSd_Cc_unbound <- 0.186
    label("Proportional residual SD for unbound naproxen (fraction)")
    # Table 2 sigma_U,prop = 18.6% (RSE 11.0%).

    # ============================================================
    # Pain-intensity (PI) model (Bjornsson 2011 Table 3 final estimates).
    # PI(t) = PI_baseline * (1 - placebo(t)) * (1 - drug(t)), where
    # placebo(t) = Pmax * (1 - exp(-kpl * t)) and drug(t) is a sigmoid
    # Emax function of unbound naproxen with Emax fixed at 1. PI is
    # clamped to [0, 100] mm per paper Methods.
    # ============================================================
    lrbase <- log(52.7)
    label("Log of typical baseline PI on 100-mm VAS (mm)")
    # Table 3 PI_baseline = 52.7 mm (RSE 13.4%); IIV 32% (RSE 27%).

    pmax <- 0.202
    label("Maximum placebo effect (fractional decrease in PI)")
    # Table 3 P_max = 20.2% (RSE 12.2%); ADDITIVE IIV 120% (RSE 16%).
    # No log-transform because the additive eta allows Pmax to cross
    # zero (paper Methods: 'for Pmax where an additive model was used,
    # allowing for PI to either increase or decrease from the baseline
    # value').

    lkpl <- log(0.237)
    label("Placebo onset rate constant k_pl (1/h)")
    # Table 3 k_pl = 0.237 /h (RSE 68.8%); IIV 43% (RSE 39%).

    lemax <- fixed(log(1))
    label("Maximum drug effect on PI (fraction, fixed at 1)")
    # Paper Methods 'Pain intensity model' + Table 3: 'Emax was
    # estimated close to 1, the upper boundary, and was therefore
    # fixed to 1.'

    lec50 <- log(0.135)
    label("Unbound naproxen concentration for half-maximum drug effect, EC50 (umol/L)")
    # Table 3 EC50 = 0.135 umol/L (RSE 10.4%); IIV 120% (RSE 21%).

    lhill <- log(1.61)
    label("Sigmoid Emax shape factor (paper symbol gamma; unitless)")
    # Table 3 gamma = 1.61 (RSE 12.4%); no IIV.

    etalrbase ~ log(1 + 0.32^2)
    # Table 3 IIV PI_baseline = 32% CV (RSE 27%).

    # ADDITIVE IIV on Pmax. The paper's footnote reports IIV "in % of
    # the parameter estimate", so SD(eta_pmax) = 1.20 * P_max = 1.20 *
    # 0.202 = 0.2424; variance = SD^2 = 0.05876.
    etapmax ~ (1.20 * 0.202)^2
    # Table 3 IIV P_max = 120% (RSE 16%, additive).

    etalkpl ~ log(1 + 0.43^2)
    # Table 3 IIV k_pl = 43% CV (RSE 39%).

    etalec50 ~ log(1 + 1.20^2)
    # Table 3 IIV EC50 = 120% CV (RSE 21%).

    addSd_PI <- 7.82
    label("Additive residual SD for pain intensity (mm)")
    # Table 3 sigma_PI = 7.82 mm (RSE 13.3%).

    # ============================================================
    # Time-to-event (rescue medication) hazard (Bjornsson 2011 Table 3).
    # Weibull baseline hazard with log-linear effects of PI(t) and an
    # additional (PI_baseline - 55) modulator of the PI(t) slope:
    # h(t) = [lambda * alpha * (lambda * (t - 1.5))^(alpha - 1)] *
    #         exp([e_pi_haz + e_pibase_haz * (PI_baseline - 55)] * PI(t))
    # for t >= 1.5 h; h(t) = 0 for t < 1.5 h (paper Methods 'Model for
    # request of rescue medication').
    # ============================================================
    llambda_haz <- log(0.00999)
    label("Weibull scale parameter lambda for the rescue-medication hazard (1/h)")
    # Table 3 lambda = 0.00999 (RSE 15.6%).

    lalpha_haz <- log(0.729)
    label("Weibull shape parameter alpha for the rescue-medication hazard (unitless)")
    # Table 3 alpha = 0.729 (RSE 9.9%). alpha < 1 implies the hazard at
    # a given PI level decreases with time (paper Discussion).

    e_pi_haz <- 0.0782
    label("Log-linear coefficient on PI(t) in the hazard at the cohort-median baseline PI = 55 mm (1/mm)")
    # Table 3 q_PI = 0.0782 (RSE 9.2%).

    e_pibase_haz <- -0.00261
    label("Modulator of e_pi_haz by (PI_baseline - 55 mm) (1/mm^2)")
    # Table 3 q_baseline = -0.00261 (RSE 19.2%). The total slope on
    # PI(t) is (e_pi_haz + e_pibase_haz * (PI_baseline_i - 55)); higher
    # baseline PI gives a smaller slope and therefore a lower hazard
    # at a given PI level (paper Results: 'patients with a high
    # baseline PI had a lower hazard at a given PI ...').
  })

  model({
    # ------------------------------------------------------------
    # Reference constants (paper Methods).
    # ------------------------------------------------------------
    pibase_median <- 55  # cohort-median baseline PI used in the hazard

    # ------------------------------------------------------------
    # Individual PK parameters.
    # ------------------------------------------------------------
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    mtt_naproxcinod <- exp(lmtt_naproxcinod + etalmtt_naproxcinod)
    nn_naproxcinod  <- exp(lnn_naproxcinod  + etalnn_naproxcinod)
    mtt_naproxen    <- exp(lmtt_naproxen    + etalmtt_naproxen)
    nn_naproxen     <- exp(lnn_naproxen     + etalnn_naproxen)

    bmax <- exp(lbmax + etalbmax)
    km   <- exp(lkm)

    fdepot_naproxcinod <- exp(lfdepot_naproxcinod)
    fdepot_naproxen    <- 1

    # Savic depot-drain rate ktr = (NN + 1) / MTT for each chain.
    ktr_naproxcinod <- (nn_naproxcinod + 1) / mtt_naproxcinod
    ktr_naproxen    <- (nn_naproxen    + 1) / mtt_naproxen

    kel <- cl / vc

    # ------------------------------------------------------------
    # ODE system. Two parallel transit-absorption chains feed a
    # shared naproxen central compartment. The user routes a dose to
    # exactly one of depot_naproxcinod or depot_naproxen by setting
    # cmt on the event record; the corresponding transit() call reads
    # podo() / tad() for that compartment and produces the gamma-PDF
    # input rate, while f(depot_*) <- 0 suppresses the bolus.
    # ------------------------------------------------------------
    d/dt(depot_naproxcinod) <- transit(nn_naproxcinod, mtt_naproxcinod, fdepot_naproxcinod) -
                                 ktr_naproxcinod * depot_naproxcinod
    d/dt(depot_naproxen)    <- transit(nn_naproxen,    mtt_naproxen,    fdepot_naproxen) -
                                 ktr_naproxen    * depot_naproxen

    d/dt(central) <- ktr_naproxcinod * depot_naproxcinod +
                     ktr_naproxen    * depot_naproxen -
                     kel * central

    f(depot_naproxcinod) <- 0
    f(depot_naproxen)    <- 0

    # ------------------------------------------------------------
    # Concentrations: central holds unbound naproxen (umol). Total
    # naproxen is derived from unbound via the saturable albumin-
    # binding equation (paper Methods 'Pharmacokinetic model').
    # ------------------------------------------------------------
    Cc_unbound <- central / vc
    Cc         <- Cc_unbound + bmax * Cc_unbound / (km + Cc_unbound)

    # ------------------------------------------------------------
    # Pain-intensity model. Multiplicative combination of placebo
    # and drug effects on the individual baseline PI (paper Methods
    # 'Model for PI'). PI is clamped to [0, 100] mm per paper.
    # ------------------------------------------------------------
    pibase_i <- exp(lrbase + etalrbase)
    pmax_i   <- pmax + etapmax
    kpl      <- exp(lkpl + etalkpl)
    emax     <- exp(lemax)
    ec50     <- exp(lec50 + etalec50)
    hill     <- exp(lhill)

    # Clamp Cc_unbound to non-negative before raising to the Hill power.
    # rxode2's stiff solver can return a tiny negative numerical noise
    # value (~1e-18) on the placebo arm where central is identically
    # zero; raising a slightly negative base to a fractional power
    # returns NaN and poisons the entire downstream PI / hazard /
    # survival chain.
    cu_pos   <- ifelse(Cc_unbound < 0, 0, Cc_unbound)
    placebo  <- pmax_i * (1 - exp(-kpl * t))
    drug_eff <- emax * (cu_pos^hill) / (ec50^hill + cu_pos^hill)

    pi_raw <- pibase_i * (1 - placebo) * (1 - drug_eff)
    PI     <- ifelse(pi_raw > 100, 100, ifelse(pi_raw < 0, 0, pi_raw))

    # ------------------------------------------------------------
    # Time-to-event hazard for rescue-medication request. The Weibull
    # scale-shape baseline is modulated by exp(slope_pi * PI(t)),
    # where slope_pi = e_pi_haz + e_pibase_haz * (PI_baseline - 55).
    # Time is shifted by 1.5 h to match the protocol's rescue-
    # medication abstention window; the hazard is zero for t < 1.5 h.
    # A small offset DEL = 1e-6 keeps (lambda * t_haz)^(alpha - 1)
    # finite at t_haz = 0 without affecting the integrated hazard.
    # ------------------------------------------------------------
    t_haz <- ifelse(t < 1.5, 0, t - 1.5)

    lambda_haz <- exp(llambda_haz)
    alpha_haz  <- exp(lalpha_haz)
    del        <- 1e-6
    h0 <- lambda_haz * alpha_haz *
            (lambda_haz * (t_haz + del))^(alpha_haz - 1)

    slope_pi <- e_pi_haz + e_pibase_haz * (pibase_i - pibase_median)

    hazard <- ifelse(t < 1.5, 0, h0 * exp(slope_pi * PI))

    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur <- exp(-cumhaz)

    # ------------------------------------------------------------
    # Observation models. Total naproxen has combined add+prop
    # residual error; unbound naproxen and PI each have a single
    # residual term. The hazard and survival outputs do not carry
    # a residual error -- they are deterministic forward-simulation
    # outputs.
    # ------------------------------------------------------------
    Cc         ~ add(addSd) + prop(propSd)
    Cc_unbound ~ prop(propSd_Cc_unbound)
    PI         ~ add(addSd_PI)
  })
}
