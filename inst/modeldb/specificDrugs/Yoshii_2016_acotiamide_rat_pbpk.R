Yoshii_2016_acotiamide_rat_pbpk <- function() {
  description <- paste(
    "Preclinical (rat).",
    "PBPK/PD (semi-mechanistic, single-organ).",
    "Distribution of the gastroprokinetic agent acotiamide (Z-338 / YM443)",
    "from blood into rat stomach, and the resulting inhibition of gastric",
    "acetylcholinesterase (AChE), after a single 1.85 umol/kg intravenous",
    "bolus. Blood disposition is biexponential and is carried here as a",
    "two-compartment system whose analytic bolus solution is the source",
    "paper's Eq. 2. The stomach is resolved as three serial spaces: a",
    "vascular space perfused at gastric blood flow, a precursor pool",
    "reached from the vascular space by a blood-flow-independent,",
    "carrier-mediated permeation clearance, and a deep pool in which drug",
    "interacts non-specifically with cellular components under linear",
    "conditions. The stomach draws from the blood but does not feed back",
    "into it: the source paper fitted blood concentrations independently",
    "and used the resulting profile as a forcing function for the tissue",
    "model. Acetylcholine in the stomach follows a Dayneka Model II",
    "indirect response in which the unbound acotiamide concentration in",
    "the precursor pool inhibits the AChE-mediated hydrolysis rate kout.",
    "Gastric blood flow, the vascular-space volume, the influx permeation",
    "clearance and the AChE IC50 were fixed from independent sources; the",
    "remaining parameters were fitted in Phoenix WinNonlin 6.1.",
    "No between-subject variability or residual-error magnitude is",
    "reported (the model was fitted to group mean profiles), so all",
    "variance terms are fixed to zero."
  )
  reference <- paste(
    "Yoshii K, Iikura M, Hirayama M, Toda R, Kawabata Y.",
    "Physiologically-based pharmacokinetic and pharmacodynamic modeling",
    "for the inhibition of acetylcholinesterase by acotiamide, a novel",
    "gastroprokinetic agent for the treatment of functional dyspepsia, in",
    "rat stomach.",
    "Pharmaceutical Research. 2016;33(2):292-298.",
    "doi:10.1007/s11095-015-1787-y.",
    sep = " "
  )
  vignette <- "Yoshii_2016_acotiamide_rat_pbpk"

  # The three stomach spaces are the source paper's own constructs (Fig. 2 and
  # Eqs. 3-5): "vascular space", "precursor pool" and "deep pool". They are
  # declared paper-specific rather than mapped onto the canonical
  # `is_stomach` / `int_stomach` / `bound_stomach` PBPK sub-compartment
  # prefixes because the authors explicitly decline the anatomical reading --
  # "Acotiamide in the precursor pool of stomach plays an important role in
  # the elevation of ACh in the stomach though the anatomical meanings of
  # these compartments remain unclear" (Discussion). The Discussion does
  # offer candidate identities (Ve resembles the inulin / extracellular
  # space; VT matches the 12.0% cytosolic fraction; the deep pool is
  # suggested to be organelle) but presents them as hypotheses, and
  # `bound_<organ>` is in any case reserved for SATURABLE pools whereas this
  # deep pool is stated to "operate under linear conditions".
  #
  # `ach` carries a CONCENTRATION (nmol / g of stomach tissue), not an
  # amount -- Eq. 9 is written directly as dR/dt with R the observed ACh
  # concentration, and kin is reported in nmol/g of tissue/min.
  paper_specific_compartments <- c(
    "stomach_vascular", "stomach_precursor", "stomach_deep", "ach"
  )

  units <- list(time = "min", dosing = "nmol/kg", concentration = "uM")

  compartmentData <- list(
    # Concentrations here are BLOOD, not plasma: the source paper converted its
    # measured plasma values with Rbp = 0.84 (Eq. 1) before fitting.
    central           = list(analyte = "acotiamide", units = "nmol/kg", specimen = "whole blood", verified = TRUE),
    peripheral1       = list(analyte = "acotiamide", units = "nmol/kg", specimen = "whole blood", verified = TRUE),
    stomach_vascular  = list(analyte = "acotiamide", units = "nmol", specimen = "tissue", verified = TRUE),
    stomach_precursor = list(analyte = "acotiamide", units = "nmol", specimen = "tissue", verified = TRUE),
    stomach_deep      = list(analyte = "acotiamide", units = "nmol", specimen = "tissue", verified = TRUE),
    ach               = list(analyte = "acetylcholine", units = "nmol/g", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species        = "rat (male Sprague-Dawley)",
    n_subjects     = 66L,
    n_studies      = 1L,
    age_range      = "6 to 7 weeks",
    sex_female_pct = 0,
    disease_state  = paste(
      "Healthy male Sprague-Dawley rats (Charles River Japan), housed at",
      "23 +/- 3 degrees C and 55 +/- 20% humidity on a 12 h light/dark",
      "cycle with food and water ad libitum, acclimated for at least one",
      "week before the experiments."
    ),
    dose_range     = paste(
      "Acotiamide 1.85 umol/kg dissolved in 5% glucose solution as a single",
      "intravenous bolus into the femoral vein under isoflurane",
      "anaesthesia. Blood sampled from the abdominal aorta at 5, 10, 15 and",
      "30 min and 1, 2, 4, 6, 8, 24 and 48 h; the stomach was excised,",
      "rinsed and homogenised at each of those times for the acotiamide",
      "assay, and over 5 min to 4 h for the acetylcholine assay."
    ),
    regions        = "Japan (Zeria Pharmaceutical Co., Ltd, Saitama)",
    notes          = paste(
      "n_subjects is the 66 rats dosed across six separate experiments;",
      "each reported time point is the mean +/- S.E. of six rats, and the",
      "model was fitted to those group mean profiles rather than to",
      "individual animals, which is why no between-subject variability is",
      "reported. Blood concentrations were obtained from measured plasma",
      "concentrations via Ca = Cp * Rbp with a blood-to-plasma ratio",
      "Rbp = 0.84 (source paper Eq. 1); that conversion is a data-handling",
      "step applied before fitting and is therefore NOT part of the model",
      "below, whose predicted concentrations are BLOOD concentrations.",
      "Acotiamide and acetylcholine were measured by LC-MS/MS. Model",
      "fitting used the Phoenix model of WinNonlin version 6.1 (Pharsight).",
      "Fit quality for the PBPK part: Loglik 36.0, AIC -61.9 on 26",
      "observations and 5 parameters (source paper Table II); for the",
      "PBPK/PD part with all three PD parameters free, Loglik 43.9 and AIC",
      "-75.8. Animal protocols approved by the Animal Care and Use",
      "Committee of the Central Research Laboratories, Zeria",
      "Pharmaceutical Co., Ltd."
    )
  )

  ini({
    # ================================================================
    # Blood disposition (source paper Table I; source paper Eq. 2).
    #
    # Eq. 2 is the analytic single-bolus solution of a two-compartment
    # model, C1 = D(alpha-k2)/(V1(alpha-beta)) * exp(-alpha*t) +
    # D(k2-beta)/(V1(alpha-beta)) * exp(-beta*t), whose micro-constants
    # are exactly the four Table I rows below: k1 is the
    # central-to-peripheral rate constant, k2 the peripheral-to-central
    # rate constant, and CLtot/V1 the elimination rate constant. The ODE
    # form in model() is written instead of the closed form so that the
    # model accepts arbitrary dosing; the vignette checks the two against
    # each other to machine precision.
    #
    # V1 and CLtot are the only Table I parameters reported PER KILOGRAM;
    # every stomach volume and flow is an absolute per-rat value. Body
    # weight cancels exactly out of the blood concentration (dose per kg
    # divided by volume per kg), so the central and peripheral
    # compartments carry amount PER KILOGRAM and no body-weight covariate
    # is needed -- see the vignette Assumptions section.
    # ================================================================
    lvc  <- log(302)       ; label("Central blood volume V1 (mL/kg)")                              # Table I: V1 = 302 mL/kg, by fitting
    lcl  <- log(56.9)      ; label("Total blood clearance CLtot (mL/min/kg)")                      # Table I: CLtot = 56.9 mL/min/kg, by fitting
    lk12 <- log(0.126)     ; label("Central-to-peripheral rate constant k1 (1/min)")               # Table I: k1 = 0.126 1/min, by fitting
    lk21 <- log(0.0313)    ; label("Peripheral-to-central rate constant k2 (1/min)")               # Table I: k2 = 0.0313 1/min, by fitting

    # ================================================================
    # Stomach distribution (source paper Table I; source paper Eqs. 3-7).
    # Table I marks the source of each row; the three rows taken from
    # independent work (Ve, Qt and fb*PSinf) are encoded fixed(), as is
    # the in vitro IC50 below.
    # ================================================================
    lv_stomach_vascular  <- fixed(log(0.441))    ; label("Stomach vascular-space volume Ve (mL)")                                        # Table I: Ve = 0.441 mL, 'Calculated from Vi [10] and tissue weight [12]' (0.401 mL/g * 1.1 g)
    lv_stomach_precursor <- log(0.133)           ; label("Stomach precursor-pool volume VT (mL)")                                        # Table I: VT = 0.133 mL, by fitting
    lq_stomach           <- fixed(log(1.1))      ; label("Gastric blood flow Qt (mL/min)")                                               # Table I: Qt = 1.1 mL/min, Hosseini-Yeganeh and McLachlan [12]
    lclin_stomach        <- fixed(log(0.174))    ; label("Vascular-space to precursor-pool influx permeation clearance fb*PSinf (mL/min)")  # Table I: fb*PSinf = 0.174 mL/min, Yoshii et al., 2011 [10]
    lclef_stomach        <- log(0.00600)         ; label("Precursor-pool to vascular-space efflux permeation clearance fu*PSeff (mL/min)")  # Table I: fu*PSeff = 0.00600 mL/min, by fitting
    lkin_stomach_deep    <- log(0.0000320)       ; label("Precursor-pool to deep-pool association rate constant kass (1/min)")           # Table I: kass = 0.0000320 1/min, by fitting
    lkout_stomach_deep   <- log(0.00000485)      ; label("Deep-pool to precursor-pool dissociation rate constant kdis (1/min)")          # Table I: kdis = 0.00000485 1/min, by fitting

    # ================================================================
    # Acetylcholine indirect response (source paper Table I and Eq. 9).
    #
    # Dayneka Model II (inhibition of the loss rate): acotiamide inhibits
    # the AChE-mediated hydrolysis of ACh, so ACh rises. The maximum
    # achievable inhibition is 1 (complete), which is why Eq. 9 carries a
    # bare (1 - CT/(IC50 + CT)) factor with no Imax term and no Imax row
    # appears in Table I.
    #
    # Table I's kin / kout are the estimates obtained with IC50 FIXED to
    # the in vitro value of 1.79 uM; this is the final model and the one
    # used for the source paper's Fig. 5. Refitting with all three PD
    # parameters free gave kin = 0.00337, kout = 0.00446 and
    # IC50 = 2.10 uM (Results, 'Analysis of Stomach Concentration of
    # ACh') -- reported as a supporting analysis, not as the final model.
    # ================================================================
    lkin  <- log(0.00314)     ; label("Acetylcholine zero-order production rate kin (nmol/g of tissue/min)")             # Table I: kin = 0.00314 nmol/g of tissue/min, by fitting
    lkout <- log(0.00415)     ; label("Acetylcholine first-order AChE hydrolysis rate constant kout (1/min)")            # Table I: kout = 0.00415 1/min, by fitting
    lic50 <- fixed(log(1.79)) ; label("Acotiamide concentration giving 50 percent AChE inhibition IC50 (uM)")            # Table I: IC50 = 1.79 uM, in vitro study (Fig. 4); fixed when fitting the ACh data

    # ================================================================
    # Residual unexplained variability. The source paper fitted group
    # mean profiles in Phoenix WinNonlin and reports only Loglik and AIC
    # (Table II); no residual-error model or magnitude is given for any
    # of the three observed quantities. Fixed to zero rather than
    # invented -- see the vignette Errata.
    # ================================================================
    propSd          <- fixed(0) ; label("Proportional residual error on blood acotiamide (fraction; ZERO - not reported in source)")
    propSd_Cstomach <- fixed(0) ; label("Proportional residual error on stomach acotiamide (fraction; ZERO - not reported in source)")
    propSd_ach      <- fixed(0) ; label("Proportional residual error on stomach acetylcholine (fraction; ZERO - not reported in source)")
  })

  model({
    # ================================================================
    # 1. Individual parameters
    # ================================================================
    vc  <- exp(lvc)
    cl  <- exp(lcl)
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    v_stomach_vascular  <- exp(lv_stomach_vascular)
    v_stomach_precursor <- exp(lv_stomach_precursor)
    q_stomach           <- exp(lq_stomach)
    clin_stomach        <- exp(lclin_stomach)
    clef_stomach        <- exp(lclef_stomach)
    kin_stomach_deep    <- exp(lkin_stomach_deep)
    kout_stomach_deep   <- exp(lkout_stomach_deep)

    kin  <- exp(lkin)
    kout <- exp(lkout)
    ic50 <- exp(lic50)

    # ================================================================
    # 2. Derived stomach volumes (source paper Eqs. 6 and 7).
    #
    # Eq. 6, Vd = VT * kass / kdis, is the volume the deep pool must have
    # for the kass / kdis pair to be an equilibrium partition: at
    # equilibrium VT*CT*kass = Vd*Cd*kdis with CT = Cd. Vd is therefore
    # derived, not estimated, and correspondingly has no Table I row.
    # With the Table I values Vd = 0.133 * 0.0000320 / 0.00000485 =
    # 0.878 mL and Vstomach = 0.133 + 0.878 = 1.01 mL, against the 1.1 mL
    # anatomical stomach volume that Table I lists from the literature --
    # see the vignette Errata for that 8% gap.
    # ================================================================
    v_stomach_deep <- v_stomach_precursor * kin_stomach_deep / kout_stomach_deep   # Eq. 6
    v_stomach      <- v_stomach_precursor + v_stomach_deep                         # Eq. 7

    # ================================================================
    # 3. Concentrations. Amounts are nmol (nmol/kg in blood) and volumes
    # are mL, so every concentration below is nmol/mL = uM.
    # ================================================================
    Cc <- central           / vc                   # C1, arterial blood
    Ce <- stomach_vascular  / v_stomach_vascular    # Ce, stomach vascular space
    CT <- stomach_precursor / v_stomach_precursor   # CT, stomach precursor pool
    Cd <- stomach_deep      / v_stomach_deep        # Cd, stomach deep pool

    # ================================================================
    # 4. Blood disposition. Two-compartment ODE form of source paper
    # Eq. 2; the elimination rate constant is CLtot / V1.
    # ================================================================
    d/dt(central)     <- -(cl / vc + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ================================================================
    # 5. Stomach mass balance (source paper Eqs. 3-5), written on
    # amounts: AT = VT*CT is `stomach_precursor` and Ad = Vd*Cd is
    # `stomach_deep`, so VT*dCT/dt and Vd*dCd/dt are d/dt of those
    # amounts directly.
    #
    # The stomach draws from `central` through q_stomach but returns
    # nothing to it, exactly as published: Eq. 2 fits the blood data on
    # its own and serves as a forcing function for Eqs. 3-5, with no
    # stomach term appearing in the blood equation. CLtot already
    # absorbs whatever the stomach removes.
    #
    # NOTE on Eq. 5: the source paper prints
    #   Vd * dCd/dt = VT*CT*kass - Vd*Cd*kass,
    # with kass in BOTH terms. The second one is a typographical error
    # for kdis, and this is mechanical rather than a judgement call:
    # (a) Eq. 4's deep-pool influx term is +Vd*Cd*kdis, and the two
    # equations must be equal and opposite for the precursor/deep
    # exchange to conserve mass; (b) as printed the deep pool could
    # never reach the Vd = VT*kass/kdis equilibrium that the paper's own
    # Eq. 6 defines. kdis is used below.
    # ================================================================
    d/dt(stomach_vascular) <-
      q_stomach * (Cc - Ce) - clin_stomach * Ce + clef_stomach * CT                # Eq. 3

    d/dt(stomach_precursor) <-
      clin_stomach * Ce + kout_stomach_deep * stomach_deep -
      clef_stomach * CT - kin_stomach_deep * stomach_precursor                     # Eq. 4

    d/dt(stomach_deep) <-
      kin_stomach_deep * stomach_precursor - kout_stomach_deep * stomach_deep      # Eq. 5 (kdis; see note above)

    # ================================================================
    # 6. Whole-stomach acotiamide concentration (source paper Eq. 8).
    # The vascular space is excluded from both the numerator and the
    # denominator: Fig. 2 brackets Vstomach / Cstomach around the
    # precursor and deep pools only, and the excised stomachs were
    # rinsed and had the blood removed before assay.
    # ================================================================
    Cstomach <- (stomach_precursor + stomach_deep) / v_stomach                     # Eq. 8

    # ================================================================
    # 7. Acetylcholine indirect response (source paper Eq. 9). The pool
    # starts at its undisturbed baseline kin / kout = 0.757 nmol/g, and
    # the driver is the unbound acotiamide concentration in the
    # precursor pool, CT -- not blood and not whole stomach. That choice
    # is the paper's central claim (Results, 'PD Modeling').
    # ================================================================
    ach(0)    <- kin / kout
    d/dt(ach) <- kin - kout * (1 - CT / (ic50 + CT)) * ach                         # Eq. 9

    Cc       ~ prop(propSd)
    Cstomach ~ prop(propSd_Cstomach)
    ach      ~ prop(propSd_ach)
  })
}
