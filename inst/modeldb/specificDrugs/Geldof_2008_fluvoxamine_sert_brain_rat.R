Geldof_2008_fluvoxamine_sert_brain_rat <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). Pharmacokinetic-pharmacodynamic model",
    "relating ex vivo 5-HT transporter (SERT) occupancy in rat frontal cortex",
    "to fluvoxamine TOTAL BRAIN TISSUE concentration, fit by",
    "Geldof et al. (2008, Br J Pharmacol 154:1369-1378). This is the third of",
    "the three PK-PD models reported side by side in Table 3 of that paper; the",
    "plasma- and ECF-driven variants are",
    "Geldof_2008_fluvoxamine_sert_plasma_rat and",
    "Geldof_2008_fluvoxamine_sert_ecf_rat. Occupancy is related to the total",
    "brain concentration by a hyperbolic Bmax (Emax) function,",
    "B = Bmax * C / (EC50 + C), with no hysteresis. The driving concentration",
    "comes from a three-compartment plasma disposition coupled to the non-linear",
    "brain distribution model of the companion paper (Geldof 2008, Pharm Res",
    "25:792-804; see modellib('Geldof_2008_fluvoxamine_rat')): a lumped total-",
    "brain state whose shallow perfusion-limited (CSP) and deep-brain (CDB, the",
    "measured ECF) concentrations are recovered algebraically at each time step",
    "from a rapid-equilibrium saturable-efflux quadratic. All PK parameters are",
    "FIXED (not re-estimated in the present paper). Inter-individual variability",
    "is on EC50 only; residual variability on occupancy is additive.",
    sep = " "
  )
  reference <- paste(
    "Geldof M, Freijer JI, van Beijsterveldt L, Langlois X, Danhof M.",
    "Pharmacokinetic-pharmacodynamic modelling of fluvoxamine 5-HT transporter",
    "occupancy in rat frontal cortex.",
    "Br J Pharmacol. 2008;154(6):1369-1378.",
    "doi:10.1038/bjp.2008.179. PMCID: PMC2483389.",
    "Plasma PK parameters fixed at the mean post-hoc estimates of the upstream",
    "Geldof 2007 rat population three-compartment PK model",
    "(Eur J Pharm Sci 2007;30(1):45-55; doi:10.1016/j.ejps.2006.10.001) per",
    "Table 1 of the 2008 paper. Brain-distribution structure and the saturable",
    "efflux capacity N***max are inherited from the companion paper",
    "Geldof M, Freijer J, van Beijsterveldt L, Danhof M.",
    "Pharmacokinetic modeling of non-linear brain distribution of fluvoxamine",
    "in the rat. Pharm Res. 2008;25(4):792-804;",
    "doi:10.1007/s11095-007-9390-5; see modellib('Geldof_2008_fluvoxamine_rat').",
    sep = " "
  )
  vignette <- "Geldof_2008_fluvoxamine_sert_occupancy"

  # brain_total is a paper-specific lumped total-brain state (companion paper
  # Geldof 2008 Pharm Res Eq 10); same declaration as the already-merged
  # upstream model inst/modeldb/specificDrugs/Geldof_2008_fluvoxamine_rat.R.
  paper_specific_compartments <- c("brain_total")

  units <- list(
    time          = "min",
    dosing        = "ng",
    concentration = "ng/mL"
  )

  # Issue #482. NOTE: brain_total is a CONCENTRATION state (the lumped total
  # brain concentration CT of paper Eq 10 of the companion Pharm Res paper),
  # not an amount -- its units are therefore ng/mL rather than ng.
  compartmentData <- list(
    central     = list(analyte = "fluvoxamine", units = "ng",    specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fluvoxamine", units = "ng",    specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "fluvoxamine", units = "ng",    specimen = "plasma", verified = TRUE),
    brain_total = list(analyte = "fluvoxamine", units = "ng/mL", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    DOSE_FLUVOXAMINE_MGKG = list(
      description = "Administered fluvoxamine dose level (1, 3.7 or 7.3 mg/kg, 30 min IV infusion)",
      units       = "mg/kg",
      type        = "continuous",
      notes       = paste(
        "Screened against the PD parameter estimates and against the plasma PK",
        "parameters; no significant correlation with either was found, so dose",
        "was not retained as a covariate (paper Results p1373-p1374)."
      )
    ),
    STUDY = list(
      description = "Source protocol (brain-sampling study vs microdialysis study)",
      units       = NA,
      type        = "categorical",
      notes       = paste(
        "Screened as a covariate on the plasma PK and brain distribution",
        "parameters; no difference between the two protocols could be detected",
        "(paper Results p1373). Tables 1 and 2 nevertheless report the",
        "per-protocol mean post-hoc estimates separately; the pooled",
        "'Brain sampling + microdialysis' rows are used here."
      )
    )
  )

  population <- list(
    species        = "rat (male Wistar, Charles River Wiga GmbH, Sulzfeld, Germany)",
    n_subjects     = 47L,
    n_studies      = 2L,
    age_range      = "adult (specific age not reported; group-housed 1 week on arrival, then 2 days individually after cannulation surgery for the brain-sampling animals and 7 days for the microdialysis animals)",
    weight_range   = "226-250 g body weight at the start of the experiments",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = paste(
      "Healthy rats prepared with a permanent right-jugular-vein cannula for",
      "fluvoxamine administration and a permanent left-femoral-artery cannula",
      "for serial arterial blood sampling. The 26 microdialysis animals also",
      "carried a CMA/12 microdialysis guide cannula in the frontal cortex",
      "(AP +3.2, L -3.0, V -1.5 mm from bregma, Paxinos and Watson atlas) for",
      "sampling brain ECF. SERT occupancy in the frontal cortex was determined",
      "ex vivo by [3H]citalopram autoradiography on 20 um cryosections,",
      "expressed as a percentage of the labelling in the corresponding brain",
      "area of untreated control animals."
    ),
    dose_range     = paste(
      "Single 30 min IV infusion of fluvoxamine free base into the right",
      "jugular vein at 1 mg/kg (24 rats) or 7.3 mg/kg (23 rats) in the",
      "brain-sampling study that provided the SERT occupancy observations",
      "(flow rate 20 uL/min, BAS BeeHive pump). The companion microdialysis",
      "study (26 rats: 8 / 8 / 10 at 1 / 3.7 / 7.3 mg/kg) contributed the",
      "brain ECF concentration data."
    ),
    regions        = "preclinical (in-vivo rat); Leiden University, The Netherlands",
    notes          = paste(
      "n_subjects = 47 counts the brain-sampling animals that contributed the",
      "ex vivo SERT occupancy observations used to fit the PD model (paper",
      "Figure 4 caption: 'the total population of 47 rats'). A further 26",
      "animals in the companion microdialysis study contributed brain ECF and",
      "plasma concentration data; n_studies = 2 counts both protocols.",
      "Only ONE SERT occupancy observation could be obtained per animal, so",
      "inter- and intra-individual variability could not be separated",
      "(Discussion p1377).",
      "Relating occupancy to brain tissue rather than plasma concentration",
      "lowered the objective function by 12.2 points and reduced both the",
      "residual error and the IIV on EC50, which the authors attribute to the",
      "non-linear brain distribution (Results p1374). The ECF- and brain-tissue-",
      "driven fits share the same MVOF (370.7), Bmax, omega^2 and (to rounding)",
      "sigma^2 -- brain tissue and ECF concentrations differ only by a scaling",
      "factor in this model, so the two are reparameterisations of one fit and",
      "differ only in EC50 (14.8 vs 0.22 ng/mL).",
      "The brain tissue EC50 (14.8 ng/mL) is the ONLY one of the three reported",
      "EC50 values that lies above the 1 ng/mL bioanalytical LOQ, so this is",
      "the variant whose concentration-effect curve is anchored throughout by",
      "directly measured concentrations (paper Discussion p1377, Figure 3c)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Plasma PK -- three-compartment linear disposition, FIXED at the mean
    # post-hoc estimates of Geldof 2008 BJP Table 1, row
    # "Brain sampling + microdialysis". Not re-estimated in this paper.
    # ------------------------------------------------------------------
    lcl  <- fixed(log(29.6));  label("Plasma systemic clearance CL (mL/min)")                              # Geldof 2008 BJP Table 1, 'Brain sampling + microdialysis' row: CL = 29.6 mL/min
    lvc  <- fixed(log(294.4)); label("Plasma central volume of distribution V1 (mL)")                      # Geldof 2008 BJP Table 1, row 1: V1 = 294.4 mL
    lvp  <- fixed(log(858.1)); label("Plasma peripheral 1 volume of distribution V2 (mL)")                 # Geldof 2008 BJP Table 1, row 1: V2 = 858.1 mL
    lq   <- fixed(log(31.8));  label("Inter-compartmental clearance central <-> peripheral1 Q2 (mL/min)")  # Geldof 2008 BJP Table 1, row 1: Q2 = 31.8 mL/min
    lvp2 <- fixed(log(136.3)); label("Plasma peripheral 2 volume of distribution V3 (mL)")                 # Geldof 2008 BJP Table 1, row 1: V3 = 136.3 mL
    lq2  <- fixed(log(1.0));   label("Inter-compartmental clearance central <-> peripheral2 Q3 (mL/min)")  # Geldof 2008 BJP Table 1, row 1: Q3 = 1.0 mL/min

    # ------------------------------------------------------------------
    # Brain distribution -- FIXED at the mean post-hoc estimates of Geldof
    # 2008 BJP Table 2, row "Brain sampling + microdialysis". These are
    # inherited from the companion non-linear brain distribution paper
    # (Geldof 2008, Pharm Res 25:792-804) and were not re-estimated here;
    # the paper used them to predict each animal's ECF and brain tissue
    # concentration at its occupancy sampling time (Data analysis, p1373).
    #
    # NOTE ON PROVENANCE: the present BJP paper tabulates only kin, kout and
    # C50. The saturable-efflux capacity N***max is REQUIRED by the brain
    # distribution structure but is NOT reported in this paper; it is taken
    # from Table II of the companion Pharm Res paper, which is on disk and
    # is already extracted as modellib('Geldof_2008_fluvoxamine_rat').
    # ------------------------------------------------------------------
    lkin      <- fixed(log(0.2031)); label("Brain influx rate constant kin (1/min)")                         # Geldof 2008 BJP Table 2, 'Brain sampling + microdialysis' row: kin = 0.2031 /min
    lkout     <- fixed(log(0.0183)); label("Brain efflux rate constant kout (1/min)")                        # Geldof 2008 BJP Table 2, row 1: kout = 0.0183 /min
    lc50      <- fixed(log(710));    label("Deep-brain fluvoxamine concentration at 50% saturation of the active removal flux, C50 (ng/mL)")  # Geldof 2008 BJP Table 2, row 1: C50 = 710 ng/mL
    lNstarMax <- fixed(log(30700));  label("Lumped saturable active-efflux capacity N***max (ng/mL)")        # NOT reported in the present BJP paper. Carried from the companion paper Geldof 2008 Pharm Res 25:792-804 Table II: N***max = 30,700 (CV 92.5%), and from the already-extracted upstream model file inst/modeldb/specificDrugs/Geldof_2008_fluvoxamine_rat.R. Pharm Res Table II prints the unit as ng.h^-1, which is dimensionally inconsistent with N***max appearing as an additive term against CT and C50 in the partition quadratic; treated there and here as a typo for ng/mL. See vignette Errata.

    # ------------------------------------------------------------------
    # PD -- hyperbolic Bmax (Emax) model, estimated in THIS paper
    # (Geldof 2008 BJP Table 3, 'ECF' column):
    #   B = Bmax * C / (EC50 + C),  C = brain ECF concentration.
    # ------------------------------------------------------------------
    lemax <- log(94.9); label("Maximum SERT occupancy Bmax (% of control [3H]citalopram labelling)")  # Geldof 2008 BJP Table 3, Brain column: Bmax = 94.9 % (CV 1.1%)
    lec50 <- log(14.8); label("Total brain tissue fluvoxamine concentration at half-maximal SERT occupancy EC50 (ng/mL)")  # Geldof 2008 BJP Table 3, Brain column: EC50 = 14.8 ng/mL (CV 10.6%)

    etalec50 ~ 0.25  # Geldof 2008 BJP Table 3, Brain column: omega^2(eta_EC50) = 0.25 (CV 38.0%)

    addSd <- 5.5498; label("Additive residual SD on SERT occupancy (percentage points)")  # Geldof 2008 BJP Table 3, Brain column: sigma^2(eps1ij) = 30.8 (CV 27.1%); SD = sqrt(30.8) = 5.5498
  })

  model({
    # 1. Individual parameters
    cl   <- exp(lcl)
    vc   <- exp(lvc)
    vp   <- exp(lvp)
    q    <- exp(lq)
    vp2  <- exp(lvp2)
    q2   <- exp(lq2)

    kin      <- exp(lkin)
    kout     <- exp(lkout)
    c50      <- exp(lc50)
    NstarMax <- exp(lNstarMax)

    emax <- exp(lemax)
    ec50 <- exp(lec50 + etalec50)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 3. Brain distribution algebra (companion paper Geldof 2008 Pharm Res
    #    Eqs 8-10, 16 and Appendix Eqs 36, 47, 55-57; encoding carried from
    #    the already-extracted model Geldof_2008_fluvoxamine_rat).
    #
    #    After eliminating the (numerically unbounded) inter-compartmental
    #    diffusion constant kdiff under the rapid-equilibrium assumption, the
    #    shallow and deep brain compartments collapse to a single ODE on the
    #    lumped total brain concentration CT:
    #        dCT/dt = kin*Cp - kout*CSP                             (Eq 10)
    #    CSP and CDB are recovered from CT at each step via the quadratic
    #    obtained from the saturable Michaelis-Menten efflux out of the deep
    #    brain (Appendix Eq 47, in concentration units with VSP = VDB):
    #        CDB^2 + (C50 + N***max - CT)*CDB - CT*C50 = 0
    #    whose physical (non-negative) root is
    #        CDB = ((CT - C50 - N***max) + sqrt((C50 + N***max - CT)^2
    #               + 4*CT*C50)) / 2
    #    and the brain mass balance VT*CT = VSP*CSP + VDB*CDB with VSP = VDB
    #    gives CSP = 2*CT - CDB.
    #
    #    The published Eq 57 prints +N***max in the numerator, which sends
    #    CDB/CT to infinity as CT -> 0 and contradicts the low-CT limit
    #    CDB/CT -> C50/(C50 + N***max) implied by Eq 47. That sign is treated
    #    as a typesetting error, consistently with the upstream extraction.
    bt    <- brain_total
    diff_ <- bt - c50 - NstarMax
    disc  <- diff_ * diff_ + 4 * bt * c50
    cdb   <- (diff_ + sqrt(disc)) / 2
    csp   <- 2 * bt - cdb

    cp <- central / vc

    # 4. ODE system
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2
    d/dt(brain_total) <-  kin * cp - kout * csp

    # 5. Observations
    Cc     <- cp   # plasma fluvoxamine concentration (ng/mL); diagnostic output
    Cecf   <- cdb  # brain ECF (deep brain) fluvoxamine concentration (ng/mL); diagnostic output
    Cbrain <- bt   # total brain tissue fluvoxamine concentration (ng/mL); the PD driver

    # Hyperbolic Bmax model (paper Methods p1373) driven by the total brain
    # tissue concentration.
    sertOccupancy <- emax * Cbrain / (ec50 + Cbrain)

    sertOccupancy ~ add(addSd)
  })
}
