VelezdeMendizabal_2012_lumiracoxib_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Semi-mechanistic PD model of the formalin-induced antinociceptive",
    "response to lumiracoxib in adult female Wistar rats (Velez de",
    "Mendizabal 2012). No PK measurements were made: lumiracoxib was",
    "tracked through two virtual compartments -- intraplantar local",
    "(lumxLocal) and intrathecal central (lumxCns) -- each decaying",
    "monoexponentially at first-order rates K_D_Local and K_D_CNS",
    "from a bolus equal to the administered dose (10, 30, 100, or 300",
    "ug per route). The biphasic formalin-induced nociceptive response",
    "(flinch count per 1-min window) is modeled as the sum of an early",
    "phase PN1, a monoexponential decay from an initial pain load PN1_0",
    "with rate K_PN1 (insensitive to lumiracoxib), and a delayed phase",
    "PN2 built from upregulated COX-2 in the local and CNS compartments.",
    "Both COX-2 species are taken proportional to a pain-mediator signal",
    "MED whose time course is the analytical Erlang-transit kernel of",
    "Savic 2007 (MED0 = 1; chain length NC = 6.5; transit rate K_TR =",
    "0.233 min^-1), and the proportionality constants theta_COX2_L /",
    "theta_COX2_CNS scale MED to flinch units in the local and CNS arms",
    "respectively. Lumiracoxib inhibits upregulated COX-2 in each arm",
    "via E = 1 / (1 + LUMX) with an implicit IC50 of one dose unit (an",
    "IC50 parameter was tested and found not significant). Model is the",
    "second-pass selection (Table I of Velez de Mendizabal 2012); the",
    "IC50, delayed-COX-2, and Emax variants were rejected during model",
    "development."
  )
  reference <- paste(
    "Velez de Mendizabal N, Vasquez-Bahena D, Jimenez-Andrade JM,",
    "Ortiz MI, Castaneda-Hernandez G, Troconiz IF. (2012).",
    "Semi-mechanistic modeling of the interaction between the central",
    "and peripheral effects in the antinociceptive response to",
    "lumiracoxib in rats.",
    "AAPS J 14(4):904-914.",
    "doi:10.1208/s12248-012-9405-y.",
    sep = " "
  )
  vignette <- "VelezdeMendizabal_2012_lumiracoxib_rat"
  dosing <- c("lumxLocal", "lumxCns")
  paper_specific_compartments <- c("lumxLocal", "lumxCns")

  units <- list(
    time          = "min",
    dosing        = "ug (lumiracoxib bolus into the intraplantar and / or intrathecal virtual compartment)",
    concentration = "flinches per 1-min window (number of paw flinches; the only observed quantity, no lumiracoxib concentration was measured)"
  )

  covariateData <- list()

  population <- list(
    species        = "rat (female Wistar)",
    n_subjects     = 86L,
    n_studies      = 1L,
    age_range      = "6-7 weeks",
    weight_range   = "180-220 g",
    sex_female_pct = 100,
    disease_state  = paste(
      "Healthy rats subjected to the formalin test: 50 uL of 1% formalin",
      "injected subcutaneously into the dorsal surface of the right",
      "hind paw to elicit a biphasic nociceptive response measured as",
      "the number of paw flinches per 1-min window every 5 min for 60",
      "min after formalin injection."
    ),
    dose_range     = paste(
      "Lumiracoxib administered as 10, 30, 100, or 300 ug by",
      "intraplantar (i.pl., 20 min before formalin) or intrathecal",
      "(i.th., 10 min before formalin) route. Combination arm gave",
      "i.pl. + i.th. simultaneously in a fixed-ratio (13 + 13.52,",
      "26 + 27, 52 + 54, 104 + 108, 208 + 216 ug) based on per-route",
      "ED30. One control group received saline."
    ),
    regions        = "Mexico (CINVESTAV, Mexico City)",
    notes          = paste(
      "Animals were obtained from the in-house Wistar colony at CINVESTAV",
      "(Mexico City). 86 rats randomized into 14 groups of 6:",
      "4 i.pl.-dose groups, 4 i.th.-dose groups, 5 combined-route",
      "groups, 1 saline group. Intrathecal catheterization was performed",
      "under ketamine-xylazine anesthesia at least 5 days prior to dosing.",
      "Pain-related behavior was recorded as the number of flinches /",
      "shakings of the injected paw during 1-min windows every 5 min up",
      "to 60 min after formalin injection (Methods, Measurement of",
      "Antinociceptive Activity). Estimation was via NONMEM 7 FOCE",
      "INTERACTION."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Lumiracoxib virtual-compartment decay rates. Drug concentrations
    # were NOT measured; LUMX is in dose-administered units (ug) and
    # decays monoexponentially from the bolus value at the time of
    # injection (Eqs. 7-8). The IC50 in the inhibitory function E =
    # 1 / (1 + LUMX) (Eqs. 9-10) is structurally fixed to 1 ug-equivalent
    # because an explicit IC50 parameter was tested and reported as not
    # significant (Results paragraph 1).
    # ---------------------------------------------------------------------
    lkdLocal      <- log(0.129)   ; label("Local virtual-lumiracoxib decay rate K_D_Local (1/min)")    # Table I: K_D_Local = 0.129 (CV 5.34%)
    lkdCns        <- log(0.073)   ; label("Intrathecal virtual-lumiracoxib decay rate K_D_CNS (1/min)") # Table I: K_D_CNS = 0.073 (CV 22.16%)

    # ---------------------------------------------------------------------
    # COX-2 to flinch-count scaling factors (Eqs. 3-4). theta_COX2_L
    # and theta_COX2_CNS map the unit-less pain-mediator signal MED to
    # flinch counts in the local and CNS arms respectively. They are
    # the proportionality constants for upregulated COX-2 in the
    # injured paw and in the spinal cord.
    # ---------------------------------------------------------------------
    lthetaCox2Local <- log(94.0)  ; label("Local-COX-2 to flinch scaling theta_COX-2_L (flinches)")     # Table I: theta_COX-2_L = 94.0 (CV 7.79%); IAV 15.87% (CV 42.85%)
    lthetaCox2Cns   <- log(28.5)  ; label("CNS-COX-2 to flinch scaling theta_COX-2_CNS (flinches)")     # Table I: theta_COX-2_CNS = 28.5 (CV 24.91%)

    # ---------------------------------------------------------------------
    # First-phase pain. PN1 is a monoexponential decay (Eq. 1) from
    # an initial pain load PN1_0 at the time of formalin injection
    # (t = 0). The first phase is reported as insensitive to lumiracoxib.
    # ---------------------------------------------------------------------
    lpn10           <- log(18.7)  ; label("Initial first-phase pain load PN1_0 (flinches)")              # Table I: PN1,0 = 18.7 (CV 3.14%); IAV 22.4% (CV 20.51%)
    lkpn1           <- log(0.279) ; label("First-phase pain decay rate K_PN1 (1/min)")                   # Table I: K_PN1 = 0.279 (CV 6.37%)

    # ---------------------------------------------------------------------
    # Second-phase pain-mediator transit kernel (Eq. 2). MED is the
    # analytical Erlang-transit solution with MED0 = 1, chain length
    # NC, and transit rate K_TR. NC is estimated as a continuous value
    # (6.5) and the Erlang factorial NC! is implemented as the gamma
    # function Gamma(NC+1) via lgammafn() in the model() block so that
    # non-integer NC is handled smoothly. K_TR is the canonical
    # transit-rate name from parameter-names.md (ktr).
    # ---------------------------------------------------------------------
    lnc             <- log(6.5)   ; label("Erlang transit-chain length NC (count)")                      # Table I: NC = 6.5 (CV 4.92%); IAV 11.09% (CV 21.13%)
    lktr            <- log(0.233) ; label("Transit-chain rate K_TR (1/min)")                             # Table I: K_TR = 0.233 (CV 4.20%)

    # ---------------------------------------------------------------------
    # Inter-animal variability (paper Methods: "exponential" IIV
    # model; only the three parameters that Table I reports IAV for
    # carry an eta). omega^2 = log(CV^2 + 1) is the log-scale variance
    # for an exponential IIV with the reported IAV CV%.
    # ---------------------------------------------------------------------
    etalpn10           ~ log(0.224^2 + 1)   # Table I: PN1,0 IAV = 22.4% CV -> var = log(1.05018) = 0.04895
    etalthetaCox2Local ~ log(0.1587^2 + 1)  # Table I: theta_COX-2_L IAV = 15.87% CV -> var = log(1.02518) = 0.02486
    etalnc             ~ log(0.1109^2 + 1)  # Table I: NC IAV = 11.09% CV -> var = log(1.01230) = 0.01223

    # ---------------------------------------------------------------------
    # Residual error. "Residual error was modeled using an additive
    # model" (Methods, Data Analysis). The estimate is 2.93 flinches
    # (CV 19.11%, Table I).
    # ---------------------------------------------------------------------
    addSd_pain      <- 2.93        ; label("Additive residual SD on flinch counts (flinches)")           # Table I: Residual error = 2.93 (CV 19.11%)
  })

  model({
    # Individual structural parameters (exponential IIV on PN1_0,
    # theta_COX-2_L, and NC; the remaining parameters have no IIV per
    # Table I).
    kdLocal         <- exp(lkdLocal)
    kdCns           <- exp(lkdCns)
    thetaCox2Local  <- exp(lthetaCox2Local + etalthetaCox2Local)
    thetaCox2Cns    <- exp(lthetaCox2Cns)
    pn10            <- exp(lpn10 + etalpn10)
    kpn1            <- exp(lkpn1)
    nc              <- exp(lnc + etalnc)
    ktr             <- exp(lktr)

    # ---------------------------------------------------------------------
    # Time-since-formalin gating. By convention in this model the
    # simulation variable t is referenced to the formalin injection
    # (t = 0). Lumiracoxib bolus events are scheduled at t = -20 min
    # (intraplantar) and / or t = -10 min (intrathecal) relative to
    # formalin. PN1 and MED are 0 for t <= 0 (no formalin -> no
    # response); the (t > 0) Heaviside gate plus a 1e-12 floor on the
    # log-argument keeps the analytical Erlang kernel finite and
    # non-NaN before formalin.
    # ---------------------------------------------------------------------
    post            <- (t > 0)
    tSafe           <- post * t + (1 - post) * 1e-12

    # Eq. 1 (analytical form): PN1(t) = PN1_0 * exp(-K_PN1 * t), t > 0.
    pn1             <- post * pn10 * exp(-kpn1 * t)

    # Eq. 2: MED(t) = (K_TR * t)^NC * exp(-K_TR * t) / Gamma(NC+1),
    # MED0 = 1, for t > 0. Implemented in log space via lgammafn()
    # so that non-integer NC works smoothly.
    med             <- post * exp(nc * log(ktr * tSafe) - ktr * tSafe - lgammafn(nc + 1))

    # Eqs. 7-8: virtual-lumiracoxib monoexponential decay. LUMXLocal
    # and LUMXCNS are dosed in ug via bolus events at the time of
    # intraplantar / intrathecal injection (handled by the event
    # dataset, not inside model()).
    d/dt(lumxLocal) <- -kdLocal * lumxLocal
    d/dt(lumxCns)   <- -kdCns   * lumxCns

    # Eqs. 9-10: inhibitory effect of lumiracoxib on upregulated
    # COX-2 in each arm, with structurally fixed IC50 = 1 dose unit.
    eLocal          <- 1 / (1 + lumxLocal)
    eCns            <- 1 / (1 + lumxCns)

    # Eqs. 3-4 with the Eq. 11 drug-effect substitution: COX-2 levels
    # are proportional to MED, multiplied by the inhibitory effect.
    cox2Local       <- thetaCox2Local * med * eLocal
    cox2Cns         <- thetaCox2Cns   * med * eCns

    # Eqs. 5-6: second-phase pain (sum of local + CNS upregulated
    # COX-2 contributions) and total pain response (PN1 + PN2).
    pn2             <- cox2Local + cox2Cns
    pain            <- pn1 + pn2

    pain ~ add(addSd_pain)
  })
}
