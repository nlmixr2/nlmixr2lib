Xie_2000_m3g_rat <- function() {
  description <- "Preclinical (rat, male Sprague-Dawley). Blood-brain barrier (BBB) distributional model for morphine-3-glucuronide (M3G) in rat as published by Xie et al. (2000, Br J Pharmacol): a one-compartment plasma PK driven by an unbound systemic clearance CL_u = 3.8 mL/min from the paper's Model A, coupled to a two-compartment brain model (brain 1 = sampled brain extracellular fluid via striatal microdialysis, brain 2 = deeper redistribution compartment) with asymmetric BBB exchange (separate unbound influx CL_u,in and efflux CL_u,out across the BBB) and a symmetric intercompartmental clearance Q_br between the two brain compartments. The model captures a probenecid-sensitive organic-anion transport contribution to BBB influx: CL_u,in is 1.55-fold higher under co-administered probenecid (CONMED_PROBENECID = 1) while CL_u,out, Q_br, and the two brain volumes are unchanged."
  reference <- paste(
    "Xie R, Bouw MR, Hammarlund-Udenaes M.",
    "Modelling of the blood-brain barrier transport of",
    "morphine-3-glucuronide studied using microdialysis in the rat:",
    "involvement of probenecid-sensitive transport.",
    "Br J Pharmacol. 2000;131(8):1784-1792.",
    "doi:10.1038/sj.bjp.0703759.",
    sep = " "
  )
  vignette <- "Xie_2000_m3g_rat"
  units <- list(
    time          = "min",
    dosing        = "umol",
    concentration = "umol/L (= uM; reported as unbound M3G in arterial plasma and brain ECF)"
  )

  covariateData <- list(
    CONMED_PROBENECID = list(
      description        = "Indicator for probenecid co-administration: 1 during the day-2 probenecid co-infusion in the probenecid-arm rats, 0 otherwise (control arm both days, plus day 1 of the probenecid arm before the loading dose).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no probenecid co-administration)",
      notes              = "Xie 2000 dose / regimen: 70 umol/kg intravenous probenecid loading dose followed by a 70 umol/kg/h constant infusion (flow rate 1.3 mL/kg/h of a 5% sodium bicarbonate / saline vehicle), continued for the full ~8 h experimental window starting at the blank-period start on day 2. The canonical column is treated as time-varying per subject: 0 before the loading dose, 1 from the loading dose onward through the end of the post-infusion sampling period. The day-2 control arm and day 1 of both arms run with CONMED_PROBENECID = 0; the 5% bicarbonate vehicle alone had no detectable effect on M3G PK or BBB transport (paper Discussion, p1786).",
      source_name        = "(derived from study-day + treatment-arm indicators in the original NONMEM dataset; the paper itself reports parameter estimates separately for the 'without probenecid' pool (control arm both days plus probenecid arm day 1) and the 'with probenecid' arm (probenecid arm day 2), as in paper Table 3)"
    )
  )

  population <- list(
    species        = "rat (male Sprague-Dawley, Charles River Sweden)",
    n_subjects     = 14L,
    n_studies      = 1L,
    age_range      = "not reported (animals housed for at least 1 week pre-experiment after arrival)",
    weight_range   = "280-320 g body weight at experiment start",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = "healthy rats prepared with chronic left-femoral-artery / -vein cannulae for blood sampling and IV M3G / probenecid dosing, a right-jugular-vein CMA/20 microdialysis probe for blood, and an indwelling CMA/12 striatal microdialysis probe for brain ECF sampling; animals were allowed 24 h post-surgery recovery before the first experimental day",
    dose_range     = "exponential intravenous infusion of M3G over 4 h targeting a steady-state plasma concentration of 65 uM (= 30 ug/mL) using the Stanpump CCI system with rat M3G PK parameters from Ekblom et al. 1993; for the probenecid arm on day 2 additionally 70 umol/kg IV probenecid loading dose plus 70 umol/kg/h constant infusion in 5% sodium-bicarbonate / saline vehicle (flow rate 1.3 mL/kg/h, total ~8 h infusion duration)",
    regions        = "preclinical (in-vivo rat); Uppsala University, Sweden",
    notes          = "Two experimental groups of n = 7 rats each: a control group receiving the 4 h M3G infusion on both days, and a probenecid group receiving the same M3G regimen on both days plus a continuous probenecid infusion on day 2 only. The two daily 4 h M3G infusions per rat allow each animal to serve as its own within-subject control for the probenecid effect (paper Tables 1 - 3). Brain tissue was harvested at the end of day 2 in n = 4 animals per group for total-brain assays. The BBB transport parameters (Model B in Xie 2000) are reported per gram of brain tissue; brain volumes are in uL/g-brain and brain clearances in uL/min/g-brain (paper Table 3, p1788), so the rat-brain compartments below carry per-g-brain semantics while the central plasma compartment uses whole-rat semantics. Brain accounted for about 0.02% of total body amount of M3G at steady state, so plasma kinetics are essentially independent of brain redistribution (paper Discussion p1789)."
  )

  ini({
    # ------------------------------------------------------------------
    # Plasma PK (collapsed 1-compartment central representation of the
    # paper's Model A; see Assumptions and deviations in the validation
    # vignette). Final values from Xie 2000 Results section paragraph 4:
    #   "an unbound clearance of 3.8 mL/min" (population unbound CL_u
    #    from the joint arterial + venous Model A fit, NONMEM VI, OF).
    # The corresponding central volume V_c is not tabulated; it is
    # derived from CL_u = 3.8 mL/min and the paper's reported terminal
    # half-life in venous blood t_1/2,bl = 22 min (Table 2, control
    # arm + probenecid arm pooled, days 1 and 2 -- values 22.2, 21.9,
    # 21.4, 20.0 min round to a pooled t_1/2,bl = 22 min) as
    # V_c = CL_u * t_1/2,bl / ln(2) = 3.8 * 22 / 0.6931 = 120.6 mL,
    # rounded to 0.121 L (paper-derived).
    lcl     <- log(0.0038);   label("Unbound systemic clearance CL_u (L/min)")                                 # Xie 2000 Results p1787: CL_u = 3.8 mL/min from joint arterial-venous Model A
    lvc     <- log(0.1206);   label("Central (arterial plasma) volume of distribution V_c (L)")                # Derived: V_c = CL_u * t_1/2,bl / ln(2) = 3.8 * 22 / 0.6931 = 120.6 mL (Table 2 pooled t_1/2,bl ~= 22 min)

    # ------------------------------------------------------------------
    # BBB transport parameters (Model B). Table 3 of Xie 2000, p1788,
    # reports population estimates per g-brain. The "without probenecid"
    # column pools the control arm (days 1 + 2) and the probenecid arm
    # on day 1; the "with probenecid" column corresponds to the
    # probenecid arm on day 2. CL_u,out, Q_br, V_u,br1, V_u,br2 are
    # reported as a single value for both arms (no probenecid effect
    # demonstrated for these parameters; see paper Results p1787).
    lcluin       <- log(0.00011);  label("Unbound BBB influx clearance CL_u,in (L/min/g-brain) without probenecid")  # Table 3: 0.11 uL/min/g-brain (RSE 25%) -- without-probenecid pool
    lcluout      <- log(0.00115);  label("Unbound BBB efflux clearance CL_u,out (L/min/g-brain)")                    # Table 3: 1.15 uL/min/g-brain (RSE 34%) -- common to both arms
    lqbr         <- log(0.00095);  label("Intercompartmental clearance Q_br between brain 1 and brain 2 (L/min/g-brain)") # Xie 2000 Results p1787: Q_br = 0.95 uL/min/g-brain (preclinical fit; no per-arm split)
    lvubr1       <- log(0.000028); label("Brain 1 (brain ECF) unbound volume of distribution V_u,br1 (L/g-brain)")    # Xie 2000 Results p1787: V_u,br1 = 28 uL/g-brain
    lvubr2       <- log(0.000205); label("Brain 2 (deep brain) unbound volume of distribution V_u,br2 (L/g-brain)")   # Xie 2000 Results p1787: V_u,br2 = 205 uL/g-brain, derived as V_u,app - V_u,br1 (V_u,app = 250 uL/g-brain fixed in regression per Eq 3)

    # ------------------------------------------------------------------
    # Probenecid covariate effect on CL_u,in. Xie 2000 fits separate
    # population CL_u,in values for the with-probenecid (0.17) and
    # without-probenecid (0.11) data subsets (Table 3); the canonical
    # form here is a multiplicative exponential effect on the
    # reference (without-probenecid) value so that CONMED_PROBENECID = 0
    # reproduces the reference 0.11 uL/min/g-brain and CONMED_PROBENECID = 1
    # gives 0.11 * exp(log(0.17 / 0.11)) = 0.17 uL/min/g-brain.
    e_conmed_probenecid_cluin <- log(0.17 / 0.11); label("Exponential effect of CONMED_PROBENECID on CL_u,in (log(0.17/0.11) = 0.4353; the 1.55-fold paper-reported increase)") # Xie 2000 Table 3: with-probenecid CL_u,in = 0.17 uL/min/g-brain vs reference 0.11 uL/min/g-brain

    # ------------------------------------------------------------------
    # Inter-animal variability (Xie 2000 "IAV", exponential model per
    # equation 5: P_i = P_pop * exp(Z_i), Z_i ~ N(0, omega^2); the IAV
    # column in Table 3 is reported as a percentage, so the variance
    # entered here is (IAV / 100)^2 on the log scale).
    #
    # Table 3 separately reports IAV(CL_u,in) = 49% in the
    # without-probenecid pool and 20% in the with-probenecid arm;
    # the encoding below uses the without-probenecid IAV (49%) as a
    # single across-arm IIV because nlmixr2lib does not support
    # arm-dependent variance directly. Discussed in vignette
    # Assumptions and deviations.
    etalcluin    ~ 0.2401  # (0.49)^2; without-probenecid IAV(CL_u,in) = 49% (RSE 25%), Xie 2000 Table 3
    etalcluout   ~ 0.2500  # (0.50)^2; IAV(CL_u,out) = 50% (RSE 65%), Xie 2000 Table 3, common to both arms

    # ------------------------------------------------------------------
    # Residual error. Xie 2000 reports residual variability = 0.14
    # (RSE 20%) under the proportional model Y_obs = Y_pred * (1 + eps)
    # with eps ~ N(0, sigma^2) and explicitly states sigma is the
    # standard deviation (paper Methods, equation 6 description). Mapped
    # directly to propSd on the brain ECF observation.
    propSd <- 0.14; label("Proportional residual error on brain ECF observations (fraction; SD on the multiplicative noise term)") # Xie 2000 Table 3: residual variability = 0.14 (RSE 20%); paper Methods equation 6 confirms sigma is the SD
  })

  model({
    # ---- Individual structural parameters ----
    cl      <- exp(lcl)
    vc      <- exp(lvc)
    cluin   <- exp(lcluin + e_conmed_probenecid_cluin * CONMED_PROBENECID + etalcluin)
    cluout  <- exp(lcluout + etalcluout)
    qbr     <- exp(lqbr)
    vubr1   <- exp(lvubr1)
    vubr2   <- exp(lvubr2)

    # ---- Plasma elimination rate ----
    kel <- cl / vc

    # ---- Unbound concentrations driving the BBB transport ODEs ----
    # The plasma concentration is the only state-derived "concentration"
    # in the system (no protein-binding correction applied; M3G plasma
    # binding is about 7%, so unbound ~ total, paper Methods p1786).
    # The brain compartment "concentrations" are per-g-brain ratios of
    # amount to V_u,br (the latter has units of L/g-brain), giving
    # uM = umol/L on the same scale as Cu_pl.
    cu_pl     <- central     / vc
    cu_br_csf <- brain_csf   / vubr1
    cu_br_dp  <- brain_deep  / vubr2

    # ---- ODE system (Xie 2000 Model B, equation 2, plus a one-cmt
    # plasma upstream representing Model A; paper Discussion p1789
    # confirms the brain -> plasma feedback is negligible, ~0.02% of
    # body amount at steady state, so brain efflux does not feed back
    # into the central ODE).
    d/dt(central)    <- -kel * central
    d/dt(brain_csf)  <-  cluin  * cu_pl -
                         cluout * cu_br_csf -
                         qbr    * (cu_br_csf - cu_br_dp)
    d/dt(brain_deep) <-  qbr    * (cu_br_csf - cu_br_dp)

    # ---- Observation: unbound M3G in brain ECF (paper-named output) ----
    # The plasma concentration cu_pl is computed above for diagnostic /
    # simulation use but is not declared as an observation here because
    # Xie 2000 reports a residual-error estimate only for the brain ECF
    # fit (Model B, Table 3). The plasma error model belonging to Model A
    # is not tabulated in the published paper, so a plasma observation
    # with an invented residual-error magnitude would not be source-traceable.
    Cbrain_ecf <- cu_br_csf
    Cbrain_ecf ~ prop(propSd)
  })
}
