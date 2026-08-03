Ollier_2015_ropivacaine <- function() {
  description <- paste(
    "Semi-mechanistic population PK model for free and total ropivacaine",
    "after transversus abdominis plane (TAP) nerve block in adult patients",
    "undergoing liver resection surgery (Ollier 2015). One-compartment",
    "first-order absorption disposition for free ropivacaine (Rfree) with",
    "apparent clearance Cl/F and apparent volume V/F. Free ropivacaine binds",
    "reversibly to a latent unbound-binding-site pool (target) that",
    "approximates alpha-1 acid glycoprotein (AAG) via 1:1 mass-action",
    "kinetics: binding rate kbind (fixed at 100 uM^-1 h^-1 per the paper's",
    "reported insensitivity of kb over 10^2 - 10^15 uM^-1 h^-1) and",
    "dissociation constant Kd (fixed at 0.557 uM, prior-pinned in the paper",
    "and reported without RSE). The bound species (complex) plus free",
    "species give total ropivacaine. Binding-site production rate kin",
    "switches on at 12 h post-incision (the 2nd TAP bolus timepoint,",
    "representing the postoperative-inflammation-driven onset of AAG",
    "acute-phase response). Two retained covariates: allometric power on Vc",
    "for body weight (reference 70 kg, exponent 1.28); multiplicative",
    "exponential effect on Cl for major hepatic resection",
    "LIVER_RESECT_MAJOR = 1 (>=3 segments) vs 0 (2 segments), reducing free",
    "ropivacaine clearance from 1310 L/h to 620 L/h (53% drop). The paper's",
    "third retained covariate, a per-subject postoperative fibrinogen-slope",
    "effect on kin (beta = 0.422), was dropped in this packaged model",
    "because the population mean fibrinogen slope used to centre the",
    "covariate was not reported in any source on disk; see vignette",
    "Assumptions and deviations. Ropivacaine is dosed as 5 TAP boluses of 3",
    "mg/kg (10.9 umol/kg) at 0, 12, 24, 36, 48 h post-incision by protocol.",
    "Concentrations throughout are in molar units (uM) matching the paper."
  )
  reference <- paste(
    "Ollier E, Heritier F, Bonnet C, Hodin S, Beauchesne B, Molliex S,",
    "Delavenne X.",
    "Population pharmacokinetic model of free and total ropivacaine after",
    "transversus abdominis plane nerve block in patients undergoing liver",
    "resection.",
    "Br J Clin Pharmacol. 2015;80(1):67-74.",
    "doi:10.1111/bcp.12582.",
    sep = " "
  )
  vignette <- "Ollier_2015_ropivacaine"
  units <- list(time = "h", dosing = "umol", concentration = "uM")

  covariateData <- list(
    WT = list(
      description        = "Body weight at study entry.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline weight, used in the allometric power scaling of Vc with reference 70 kg and exponent 1.28 (Table 2, Weight on V = 1.28 with %RSE 37; equation on p. 69 'ln(theta_ij) = ln(theta_pop_j) + beta_W * (ln(W_i) - ln(70)) + eta_ij' -- the log-centred exponential form is algebraically the power form (W/70)^beta_W).",
      source_name        = "W"
    ),
    LIVER_RESECT_MAJOR = list(
      description        = "Binary indicator for major hepatic resection: 1 = three or more liver segments resected, 0 = two segments resected.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (2 segments)",
      notes              = "Time-fixed per subject. Ollier 2015 reports the raw counts n=7 with 2 segments, n=2 with 3, n=6 with 4, n=1 with 5 (Table 1); dichotomised into NSH=0 (2 segments) and NSH=1 (3, 4 or 5 segments) per Methods 'Covariate model'. Used as a multiplicative exponential effect on free ropivacaine clearance: `Cl = exp(lcl + etalcl) * exp(beta * LIVER_RESECT_MAJOR)` with `beta = log(620/1310) = -0.7480`, i.e., Cl drops from 1310 L/h in the reference cohort to 620 L/h in the major-resection cohort (Table 2 shows the two values directly rather than the coefficient).",
      source_name        = "NSH"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 16L,
    n_studies      = 1L,
    age_range      = "29-77 years",
    age_median     = "65 years",
    weight_range   = "44-130 kg",
    weight_median  = "75.5 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported in source paper.",
    disease_state  = "Adult patients undergoing hepatic resection for hepatocellular carcinoma (n=4), metastasis (n=10), or hepatocellular adenoma (n=2). Randomised to ropivacaine (n=19 enrolled; 3 excluded for pharmacokinetic profiles incompatible with a correctly placed TAP catheter, leaving n=16 for the analysis).",
    dose_range     = "Fixed regimen of five 3 mg/kg (10.9 umol/kg) ropivacaine boluses via a transversus abdominis plane (TAP) catheter at 0, 12, 24, 36, and 48 h post-surgical incision (Methods 'Ropivacaine administration'). The vignette packages a 210 mg (764.7 umol) bolus per dose for a 70 kg reference weight.",
    regions        = "France (single-centre; CHU Saint-Etienne).",
    n_excluded     = 3L,
    notes          = "Baseline preoperative biological data (Table 1) include plasma creatinine 71 (45-135) uM, prothrombin time 100 (87-100) %, AST 35 (25-61) IU/L, fibrinogen 3.0 (1.1-5.6) g/L. Postoperative day-0/1/2 fibrinogen medians 3.6 / 4.5 / 5.8 g/L (see vignette Assumptions and deviations for why the fibrinogen-slope covariate on kin was dropped)."
  )

  ini({
    # Structural PK parameters. Values from Ollier 2015 Table 2 ('Values of
    # the pharmacokinetic parameters'). Ropivacaine is administered by TAP
    # block (extravascular); apparent V and Cl correspond to V/F and Cl/F
    # (Methods 'Base model': 'As the bioavailability (F) could not be
    # quantified, the parameter values correspond to the ratios V/F and
    # Cl/F, respectively').
    lka   <- log(3.58);  label("Absorption rate ka (1/h)")                                # Table 2 ka = 3.58 1/h (%RSE 33)
    lvc   <- log(103);   label("Apparent central volume of distribution V/F (L)")         # Table 2 V  = 103  L   (%RSE 12), reference weight 70 kg
    lcl   <- log(1310);  label("Apparent free-ropivacaine clearance Cl/F (L/h), reference LIVER_RESECT_MAJOR = 0")  # Table 2 Cl NSH=0 = 1310 L/h (%RSE 16)

    # Binding kinetics. kbind (paper kb) is fixed at 100 uM^-1 h^-1 -- the
    # paper reports the value 10^12 but explicitly notes 'different values
    # for kb from 100 to 10^15 uM^-1 h^-1 have been tested, and no impact
    # on ropivacaine pharmacokinetics has been observed' (Methods 'Base
    # model'). The lower bound of the paper-authorised range is used here
    # for numerical stability of the LSODA solver; the binding half-life
    # under kbind = 100 is ~45 seconds, ~4x faster than the fastest PK
    # process (kel = 1310/103 = 12.7 1/h, half-life ~200 s), which is
    # enough to preserve the quasi-equilibrium approximation used by the
    # authors. See vignette Assumptions and deviations.
    #
    # kdiss (paper Kd) is fixed at 0.557 uM -- the paper's fitted value is
    # essentially the posterior mean of a strong normal prior N(0.58, 0.05)
    # from Aarons et al. and Table 2 reports it as '0.557 ( - )', i.e.,
    # without %RSE, indicating a prior-pinned quantity. Encoded as fixed()
    # per the parameter-names 'Fixed parameters' rule.
    lkbind <- fixed(log(100));    label("Free-to-bound ropivacaine binding rate constant kbind (1/(h*uM))")  # Methods 'Base model' (paper-authorised kb range 100 - 10^15 uM^-1 h^-1)
    lkdiss <- fixed(log(0.557));  label("Ropivacaine-AAG dissociation constant Kd (uM) (prior-pinned)")     # Table 2 Kd = 0.557 uM (no RSE, prior-informed)

    # Binding-site pool (latent variable BS approximating unbound AAG).
    # kin is the zero-order production rate of unbinding sites in uM/h.
    # Table 2 labels its units as '1/h', which is inconsistent with the
    # dBS/dt equation on p. 69 (dBS/dt = kin * I(t>=t2) + kub*Rbound -
    # kb*Rfree*BS, where every other term is uM/h -- the Table 2 unit label
    # is a typo). See vignette Assumptions and deviations. Value 1.57
    # (with %RSE 12) matches the physiological range for AAG acute-phase
    # increases.
    #
    # rbase is the initial binding-site concentration BS_12 = 71.2 uM per
    # Table 2. The paper's approximation that BS is constant from t=0 to
    # t=t_2 = 12 h (Methods 'Base model': 'We made the approximation that
    # binding site concentrations do not change between the first two
    # doses') is implemented by setting target(0) = rbase * V and
    # activating the kin production term only for t >= t_kin_start = 12 h.
    lkin   <- log(1.57);   label("Binding-site production rate kin (uM/h), active for t >= 12 h post-incision")    # Table 2 kin = 1.57 uM/h (%RSE 12); units corrected from Table 2 label '1/h'
    lrbase <- log(71.2);   label("Binding-site pool baseline concentration BS_12 (uM) during 0-12 h")             # Table 2 BS_12 = 71.2 uM (%RSE 12)

    # Covariate effects
    # WT on Vc: allometric power form, reference 70 kg, exponent 1.28 per
    # Table 2 Weight on V = 1.28 (%RSE 37). Paper equation: `ln(theta_ij) =
    # ln(theta_pop_j) + beta_W * (ln(W_i) - ln(70)) + eta_ij`, which is
    # algebraically `theta_i = theta_pop * (W_i / 70)^beta_W * exp(eta)`.
    e_wt_vc <- 1.28;                          label("Allometric exponent of WT on Vc (unitless; reference 70 kg)")  # Table 2 Weight on V = 1.28 (%RSE 37)

    # LIVER_RESECT_MAJOR on Cl: multiplicative exponential form. Paper
    # equation: `ln(theta_ij) = ln(theta_pop_j) + beta_NSH * NSH_i +
    # eta_ij`. Table 2 reports the two cohort clearances directly (Cl(NSH=0)
    # = 1310 L/h, Cl(NSH=1) = 620 L/h) rather than the beta coefficient;
    # beta = log(620/1310) = -0.7480 exactly reproduces the two reported
    # values.
    e_liverresectmajor_cl <- log(620 / 1310); label("Effect of LIVER_RESECT_MAJOR on log-Cl (unitless)")           # Table 2 Cl NSH=0 = 1310, Cl NSH=1 = 620 -> log(620/1310) = -0.7480

    # Inter-subject variability. Values in Table 2 are 'Intersubject SD'
    # (log-scale SD of the individual random effect eta), per the MONOLIX
    # reporting convention -- omega^2 = SD^2 directly. IIV was reported
    # only for the five parameters V, ka, Cl, kin, and BS_12 (Results
    # 'Population pharmacokinetic model': 'Intersubject variability was
    # detected for V, ka, Cl, kin and BS_H0'); Kd and kbind carry no IIV
    # (they are fixed structural constants).
    etalka    ~ 0.729                                                                                              # Table 2 IIV SD ka  = 0.854 (%RSE 22) -> 0.854^2
    etalvc    ~ 0.204                                                                                              # Table 2 IIV SD V   = 0.452 (%RSE 19) -> 0.452^2
    etalcl    ~ 0.168                                                                                              # Table 2 IIV SD Cl  = 0.410 (%RSE 22) -> 0.410^2
    etalkin   ~ 0.0225                                                                                             # Table 2 IIV SD kin = 0.150 (%RSE 94) -> 0.150^2 (%RSE flagged as poorly estimated in the source)
    etalrbase ~ 0.196                                                                                              # Table 2 IIV SD BS_12 = 0.443 (%RSE 19) -> 0.443^2

    # Residual variability. Table 2 reports separate additive-plus-
    # proportional residuals for free ropivacaine (Cc) and total
    # ropivacaine (Rtotal). The additive component units are uM (matching
    # the concentration units of the observations).
    propSd        <- 0.159;   label("Proportional residual error on free ropivacaine (fraction)")                  # Table 2 sigma_Free,proportional = 0.159 (%RSE 16)
    addSd         <- 0.00169; label("Additive residual error on free ropivacaine (uM)")                            # Table 2 sigma_Free,additive   = 0.00169 (%RSE 45)
    propSd_Rtotal <- 0.0966;  label("Proportional residual error on total ropivacaine (fraction)")                 # Table 2 sigma_Total,proportional = 0.0966 (%RSE 22)
    addSd_Rtotal  <- 0.318;   label("Additive residual error on total ropivacaine (uM)")                           # Table 2 sigma_Total,additive     = 0.318 (%RSE 36)
  })

  model({
    # ------------------------------------------------------------------
    # Individual parameters. Weight covariate on Vc is the allometric
    # power form derived from the paper's log-centred exponential
    # equation. LIVER_RESECT_MAJOR covariate on Cl is the multiplicative
    # exponential form.
    # ------------------------------------------------------------------
    ka     <- exp(lka + etalka)
    vc     <- exp(lvc + etalvc) * (WT / 70) ^ e_wt_vc
    cl     <- exp(lcl + etalcl) * exp(e_liverresectmajor_cl * LIVER_RESECT_MAJOR)
    kin    <- exp(lkin + etalkin)
    rbase  <- exp(lrbase + etalrbase)
    kbind  <- exp(lkbind)
    kdiss  <- exp(lkdiss)
    kunbind <- kbind * kdiss

    # ------------------------------------------------------------------
    # Postoperative onset of AAG acute-phase production. The paper models
    # binding-site production as switching on at the second dose (t = 12
    # h post-incision) because the postoperative inflammatory process
    # induces a ~12 h delay before AAG concentrations begin to rise
    # (Methods 'Base model': 'The postoperative inflammatory process
    # induced a delayed (~12 h) increase in binding site concentrations
    # [...] We made the approximation that binding site concentrations do
    # not change between the first two doses'). The onset time is a
    # biological feature of the postoperative recovery, not a dosing
    # artefact; hard-coded at 12 h.
    # ------------------------------------------------------------------
    t_kin_start <- 12
    kin_active  <- kin * (t >= t_kin_start)

    # ------------------------------------------------------------------
    # Species concentrations. All three species (free ropivacaine, bound
    # ropivacaine, unbound binding sites) share the apparent volume V.
    # ------------------------------------------------------------------
    Rfree  <- central / vc
    Rbound <- complex / vc
    BS     <- target  / vc

    # ------------------------------------------------------------------
    # Binding fluxes (concentration/time, uM/h). Multiplied by vc inside
    # each amount-space d/dt(...) below.
    # ------------------------------------------------------------------
    rate_bind_conc   <- kbind   * Rfree  * BS
    rate_unbind_conc <- kunbind * Rbound

    # ------------------------------------------------------------------
    # ODE system in amount space (rxode2 convention). State variables and
    # units:
    #   depot   (umol)   -- TAP-block ropivacaine dose site
    #   central (umol)   -- free ropivacaine
    #   complex (umol)   -- bound ropivacaine
    #   target  (umol)   -- unbound binding site (latent AAG proxy)
    # Paper concentration-space equations (p. 69) are premultiplied by vc
    # to convert to amount-space rates:
    #   d(Rfree)/dt   = Input + kub*Rbound - kb*Rfree*BS - (Cl/V)*Rfree
    #   d(Rbound)/dt  = kb*Rfree*BS - kub*Rbound
    #   d(BS)/dt      = kin*I(t>=t2) + kub*Rbound - kb*Rfree*BS
    # ------------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - (cl / vc) * central - rate_bind_conc * vc + rate_unbind_conc * vc
    d/dt(complex) <-  rate_bind_conc * vc - rate_unbind_conc * vc
    d/dt(target)  <-  kin_active * vc + rate_unbind_conc * vc - rate_bind_conc * vc

    # ------------------------------------------------------------------
    # Initial condition for the binding-site pool: BS(0) = BS_12
    # (Methods 'Base model' approximation that BS is constant over 0-12
    # h, so the value fitted for BS during the first dosing interval,
    # 71.2 uM, is used as the t=0 amount).
    # ------------------------------------------------------------------
    target(0) <- rbase * vc

    # ------------------------------------------------------------------
    # Observations. Cc = free ropivacaine (uM), the paper's primary
    # analyte and the toxicity-relevant species (published free-
    # ropivacaine toxic threshold 0.55 uM per Discussion). Rtotal = total
    # ropivacaine = Rfree + Rbound (uM), the sum used for the total-drug
    # residual-error model.
    # ------------------------------------------------------------------
    Cc     <- Rfree
    Rtotal <- Rfree + Rbound

    Cc     ~ add(addSd)        + prop(propSd)
    Rtotal ~ add(addSd_Rtotal) + prop(propSd_Rtotal)
  })
}
