Andersson_2017_nicotinicAcid_lean_rat <- function() {
  description <- paste(
    "Preclinical (rat, lean Sprague-Dawley). Nicotinic acid (NiAc) - insulin -",
    "free fatty acid (FFA) model for acute and chronic NiAc exposure.",
    "NiAc PK: two-compartment disposition with two parallel Michaelis-Menten",
    "elimination pathways and endogenous synthesis, i.v. input or s.c.",
    "mini-pump input whose delivered rate rises as an error function while",
    "the diluted catheter contents are flushed. Insulin: turnover model whose",
    "production is inhibited by NiAc (Imax-Hill), with a two-state moderator",
    "feedback (first state on production, second on loss) and an",
    "insulin-driven integral controller (regulator) on production. FFA:",
    "turnover model whose production is inhibited by NiAc with a dynamic",
    "efficacy that is down-regulated by a slow NiAc action compartment, a",
    "single moderator feedback on production, and an insulin-driven integral",
    "controller (regulator) on production. Separate lean-rat fit; the obese",
    "Zucker fit is Andersson_2017_nicotinicAcid_obese_rat."
  )
  reference <- paste(
    "Andersson R, Kroon T, Almquist J, Jirstrand M, Oakes ND, Evans ND,",
    "Chappell MJ, Gabrielsson J. Modeling of free fatty acid dynamics:",
    "insulin and nicotinic acid resistance under acute and chronic",
    "treatments. J Pharmacokinet Pharmacodyn. 2017;44(3):203-222.",
    "doi:10.1007/s10928-017-9512-6"
  )
  vignette <- "Andersson_2017_nicotinicAcid_rat"
  units <- list(
    time = "h",
    dosing = "umol/kg",
    concentration = "umol/L",
    insulin = "nmol/L",
    nefa = "mmol/L"
  )

  # s.c. mini-pump infusions are dosed into `pump`, whose first-dose time
  # tafd(pump) drives the catheter-dilution erf input; i.v. infusions are
  # dosed into `central`.
  dosing <- c("pump", "central")

  covariateData <- list()

  paper_specific_compartments <- c(
    "pump",
    "regulator_insulin",
    "regulator_nefa",
    "niac_action_nefa"
  )

  compartmentData <- list(
    pump = list(
      analyte = "nicotinic acid (NiAc)",
      units = "umol/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot = list(
      analyte = "nicotinic acid (NiAc)",
      units = "umol/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "nicotinic acid (NiAc)",
      units = "umol/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "nicotinic acid (NiAc)",
      units = "umol/kg",
      specimen = "tissue",
      verified = TRUE
    ),
    insulin = list(
      analyte = "insulin",
      units = "nmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    moderator1 = list(
      analyte = "insulin moderator M1I",
      units = "nmol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    moderator2 = list(
      analyte = "insulin moderator M2I",
      units = "nmol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    regulator_insulin = list(
      analyte = "insulin integral-controller output RI",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    nefa = list(
      analyte = "free fatty acids (FFA)",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    moderator3 = list(
      analyte = "FFA moderator MF",
      units = "mmol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    regulator_nefa = list(
      analyte = "insulin-driven FFA integral-controller output RF",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    niac_action_nefa = list(
      analyte = "nicotinic acid (NiAc) action level NF",
      units = "umol/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "rat (male lean Sprague-Dawley)",
    n_subjects = 38L,
    n_studies = 7L,
    age_range = "not reported",
    weight_range = "not reported",
    sex_female_pct = 0,
    disease_state = "healthy (lean), fasted before the acute experiment",
    dose_range = paste(
      "NiAc 170 nmol/min/kg (10.2 umol/kg/h) as s.c. mini-pump infusion",
      "(acute 5 h; after 5 days continuous; after 5 days 12 h on/off;",
      "12 h then off or 3.5 h step-down) or i.v. infusion (1 h then off or",
      "3.5 h step-down)."
    ),
    regions = "Sweden (AstraZeneca, Gothenburg)",
    notes = paste(
      "Table 1 lists 7, 6, 6, 4, 5, 5 and 5 lean rats in the NiAc naive, Cont.",
      "NiAc, Inter. NiAc, NiAc Off 1 h, NiAc Stp-Dwn 1 h, NiAc Off 12 h and",
      "NiAc Stp-Dwn 12 h protocols, with 2, 2, 2, 3, 2, 2 and 3 saline-infused",
      "controls in parentheses. Conscious groups were dosed s.c. by an",
      "implanted programmable pump (iPrecio SMP200); anesthetized groups",
      "(Inactin) were dosed i.v. (1 h) or s.c. (12 h). Lean and obese rats",
      "were fitted separately; each sub-model (NiAc, insulin, FFA) was fitted",
      "sequentially in Mathematica with FOCE, and inter-study variability on",
      "insulin and FFA baselines and moderator rates was handled as fixed",
      "study effects whose weighted summaries are the Table 2 estimates."
    )
  )

  ini({
    # NiAc PK (Table 2, lean Sprague-Dawley). Methods: only ka, delta and
    # Vmax1 were estimated; the remaining PK parameters were taken from
    # Tapani et al. (ref 23), so they are fixed() here.
    lka <- log(4.27); label("First-order absorption rate from the s.c. depot ka (1/h)") # Table 2: ka = 4.27 (RSE 13%)
    ldelta_pump <- log(77.4); label("Lumped catheter diffusion parameter delta of the pump erf input (min^0.5/h)") # Table 2: delta = 77.4 (RSE 15%)
    tidle_pump <- fixed(10080); label("Pump inactivation time t0 before infusion start (7 days, entered in minutes)") # Eq. 5 text: t0 = 7 days; minutes per Fig. 7a, see vignette
    lvmax <- log(2.64); label("Maximal elimination rate, pathway 1 Vmax1 (umol/kg/h)") # Table 2: Vmax1 = 2.64 (RSE 12%)
    lkm <- fixed(log(0.235)); label("Michaelis constant, pathway 1 Km1 (umol/L)") # Table 2: Km1 = 0.235 (literature)
    lvmax2 <- fixed(log(425)); label("Maximal elimination rate, pathway 2 Vmax2 (umol/kg/h)") # Table 2: Vmax2 = 425 (literature)
    lkm2 <- fixed(log(74.5)); label("Michaelis constant, pathway 2 Km2 (umol/L)") # Table 2: Km2 = 74.5 (literature)
    lvc <- fixed(log(0.393)); label("Plasma volume of distribution Vp (L/kg)") # Table 2: Vp = 0.393 (literature)
    lvp <- fixed(log(0.172)); label("Tissue volume of distribution Vt (L/kg)") # Table 2: Vt = 0.172 (literature)
    lq <- fixed(log(0.0511)); label("Inter-compartmental distribution Cld (L/kg/h)") # Table 2: Cld = 0.0511 (literature)
    lsynt <- fixed(log(0.213)); label("Endogenous NiAc synthesis Synt (umol/kg/h)") # Table 2: Synt = 0.213 (literature)

    # Insulin model (Table 2, lean)
    lrbase_insulin <- log(0.188); label("Baseline insulin I0 (nmol/L)") # Table 2: I0 = 0.188 (RSE 9.7%)
    lkout_insulin <- log(6.58); label("Fractional turnover rate of insulin koutI (1/h)") # Table 2: koutI = 6.58 (RSE 14%)
    imax_insulin <- 0.793; label("Efficacy of NiAc on insulin production ImaxNI (fraction)") # Table 2: ImaxNI = 0.793 (RSE 11%)
    lic50_insulin <- log(0.338); label("Potency of NiAc on insulin IC50NI (umol/L)") # Table 2: IC50NI = 0.338 (RSE 15%)
    lhill_insulin <- log(3.54); label("Hill coefficient of NiAc on insulin n (unitless)") # Table 2: n = 3.54 (RSE 6.6%)
    lktol_insulin <- log(0.646); label("Turnover rate of the insulin moderators ktolI (1/h)") # Table 2: ktolI = 0.646 (RSE 28%)
    lkoutr_insulin <- log(3.94); label("Integral gain of the insulin regulator koutRI (L/nmol/h)") # Table 2: koutRI = 3.94 (RSE 17%)

    # FFA model (Table 2, lean)
    lrbase_nefa <- log(0.707); label("Baseline FFA F0 (mmol/L)") # Table 2: F0 = 0.707 (RSE 5.0%)
    lkout_nefa <- log(428); label("Fractional turnover rate of FFA koutF (1/h)") # Table 2: koutF = 428 (RSE 140%)
    imax_nefa <- fixed(1); label("Initial efficacy of NiAc on FFA production ImaxNF (fraction)") # Table 2: ImaxNF = 1 (fixed)
    lic50_nefa <- log(0.436); label("Potency of NiAc on FFA IC50NF (umol/L)") # Table 2: IC50NF = 0.436 (RSE 12%)
    lhill_nefa <- log(1.24); label("Hill coefficient of NiAc on FFA m (unitless)") # Table 2: m = 1.24 (RSE 11%)
    lktol_nefa <- log(1.21); label("Turnover rate of the FFA moderator ktolF (1/h)") # Table 2: ktolF = 1.21 (RSE 67%)
    lkoutr_nefa <- log(0.965); label("Integral gain of the insulin-driven FFA regulator koutRF (L/nmol/h)") # Table 2: koutRF = 0.965 (RSE 29%)
    lkn_nefa <- log(0.00654); label("Turnover rate of the NiAc action compartment for FFA kNF (1/h)") # Table 2: kNF = 0.00654 (RSE 65%)
    ln50_nefa <- log(3.05); label("Potency of the NiAc action compartment for FFA N50F (umol/L)") # Table 2: N50F = 3.05 (RSE 160%)
    lgam_nefa <- fixed(log(1)); label("Hill coefficient of the FFA efficacy loss phi (unitless)") # Table 2: phi = 1 (fixed)
    sloss_nefa <- 0.807; label("Long-term NiAc efficacy loss on FFA SNF (fraction)") # Table 2: SNF = 0.807 (RSE 190%)

    # Between-subject variability: Table 2 footnote a gives CV% = 100*sqrt(omega^2),
    # so omega^2 = (CV/100)^2.
    etalka ~ 0.641601 # Table 2: BSV ka 80.1%
    etalvmax ~ 0.874225 # Table 2: BSV Vmax1 93.5%
    etalsynt ~ 0.444889 # Table 2: BSV Synt 66.7%
    etalrbase_insulin ~ 0.243049 # Table 2: BSV I0 49.3%
    etalic50_insulin ~ 1.2321 # Table 2: BSV IC50NI 111%
    etalktol_insulin ~ 0.881721 # Table 2: BSV ktolI 93.9%
    etalrbase_nefa ~ 0.031684 # Table 2: BSV F0 17.8%
    etalic50_nefa ~ 0.174724 # Table 2: BSV IC50NF 41.8%
    etalktol_nefa ~ 0.341056 # Table 2: BSV ktolF 58.4%

    # Residual error (Table 2)
    propSd <- 0.313; label("Proportional residual error on NiAc (fraction)") # Table 2: sigma propN = 0.313 (RSE 5.1%)
    addSd_insulin <- 0.0699; label("Additive residual error on insulin (nmol/L)") # Table 2: sigma addI = 0.0699 (RSE 3.3%)
    addSd_nefa <- 0.130; label("Additive residual error on FFA (mmol/L)") # Table 2: sigma addF = 0.130 (RSE 3.5%)
  })

  model({
    # Individual parameters (log-normal BSV)
    ka <- exp(lka + etalka)
    delta_pump <- exp(ldelta_pump)
    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm)
    vmax2 <- exp(lvmax2)
    km2 <- exp(lkm2)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    synt <- exp(lsynt + etalsynt)

    rbase_insulin <- exp(lrbase_insulin + etalrbase_insulin)
    kout_insulin <- exp(lkout_insulin)
    ic50_insulin <- exp(lic50_insulin + etalic50_insulin)
    hill_insulin <- exp(lhill_insulin)
    ktol_insulin <- exp(lktol_insulin + etalktol_insulin)
    koutr_insulin <- exp(lkoutr_insulin)

    rbase_nefa <- exp(lrbase_nefa + etalrbase_nefa)
    kout_nefa <- exp(lkout_nefa)
    ic50_nefa <- exp(lic50_nefa + etalic50_nefa)
    hill_nefa <- exp(lhill_nefa)
    ktol_nefa <- exp(lktol_nefa + etalktol_nefa)
    koutr_nefa <- exp(lkoutr_nefa)
    kn_nefa <- exp(lkn_nefa)
    n50_nefa <- exp(ln50_nefa)
    gam_nefa <- exp(lgam_nefa)

    # Steady-state relations (Appendix 2, Eqs. 39 and 41; regulator set
    # points from Eqs. 13 and 24)
    kin_insulin <- kout_insulin * rbase_insulin
    kinr_insulin <- koutr_insulin * rbase_insulin
    kin_nefa <- kout_nefa * rbase_nefa
    kinr_nefa <- koutr_nefa * rbase_insulin

    # s.c. pump input (Eqs. 3-5): the programmed rate is dosed into `pump`,
    # which empties into the s.c. depot within seconds (kpump is a numerical
    # device, not a paper parameter) scaled by erf(t * delta / sqrt(t0)),
    # where t is the time since the pump was first switched on (h) and t0 is
    # the 7-day pump inactivation time entered in minutes (see vignette).
    kpump <- 1000
    tpump <- tafd(pump)
    fpump <- erf(tpump * delta_pump / sqrt(tidle_pump))
    if (is.na(tpump)) fpump <- 0

    # NiAc disposition (Eqs. 1-2); i.v. doses go straight into `central`.
    Cc <- central / vc
    Ct <- peripheral1 / vp
    d/dt(pump) <- -kpump * pump
    d/dt(depot) <- fpump * kpump * pump - ka * depot
    d/dt(central) <- ka * depot + synt -
      vmax * Cc / (km + Cc) - vmax2 * Cc / (km2 + Cc) -
      q * Cc + q * Ct
    d/dt(peripheral1) <- q * Cc - q * Ct

    # Endogenous NiAc steady state (Appendix 2, Eqs. 33-37): positive root of
    # a*C^2 + b*C + c = 0 with alpha = Synt. Because a < 0 and c > 0 the
    # positive root is (-b - sqrt(b^2 - 4ac)) / (2a).
    aq <- synt - vmax - vmax2
    bq <- synt * (km + km2) - vmax * km2 - vmax2 * km
    cq <- synt * km * km2
    cp0 <- (-bq - sqrt(bq * bq - 4 * aq * cq)) / (2 * aq)
    central(0) <- cp0 * vc
    peripheral1(0) <- cp0 * vp

    # Insulin (Eqs. 7-15; lean rats have the fixed efficacy ImaxNI, Eq. 16)
    inh_insulin <- 1 - imax_insulin * Cc^hill_insulin /
      (ic50_insulin^hill_insulin + Cc^hill_insulin)
    d/dt(insulin) <- kin_insulin * regulator_insulin * inh_insulin *
      rbase_insulin / moderator1 -
      kout_insulin * moderator2 / rbase_insulin * insulin
    d/dt(moderator1) <- ktol_insulin * (insulin - moderator1)
    d/dt(moderator2) <- ktol_insulin * (moderator1 - moderator2)
    d/dt(regulator_insulin) <- kinr_insulin - koutr_insulin * insulin
    insulin(0) <- rbase_insulin
    moderator1(0) <- rbase_insulin
    moderator2(0) <- rbase_insulin
    regulator_insulin(0) <- 1

    # FFA (Eqs. 20-27) with dynamic NiAc efficacy
    eff_nefa <- imax_nefa * (1 - sloss_nefa * niac_action_nefa^gam_nefa /
      (n50_nefa^gam_nefa + niac_action_nefa^gam_nefa))
    inh_nefa <- 1 - eff_nefa * Cc^hill_nefa /
      (ic50_nefa^hill_nefa + Cc^hill_nefa)
    d/dt(nefa) <- kin_nefa * regulator_nefa * inh_nefa *
      rbase_nefa / moderator3 - kout_nefa * nefa
    d/dt(moderator3) <- ktol_nefa * (nefa - moderator3)
    d/dt(regulator_nefa) <- kinr_nefa - koutr_nefa * insulin
    d/dt(niac_action_nefa) <- kn_nefa * (Cc - niac_action_nefa)
    nefa(0) <- rbase_nefa
    moderator3(0) <- rbase_nefa
    regulator_nefa(0) <- 1
    niac_action_nefa(0) <- cp0

    Cc ~ prop(propSd)
    insulin ~ add(addSd_insulin)
    nefa ~ add(addSd_nefa)
  })
}
