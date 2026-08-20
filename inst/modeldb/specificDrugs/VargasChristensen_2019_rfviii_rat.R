VargasChristensen_2019_rfviii_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Two-compartment quasi-steady-state target-mediated drug",
    "disposition (QSS-TMDD) population PK model of recombinant Factor",
    "VIII (rFVIII) and its carrier protein von Willebrand factor (VWF)",
    "in hemophilia A rats (Vargas Christensen 2019). Jointly describes",
    "total rFVIII plasma concentration, total VWF plasma concentration,",
    "and the rFVIII:VWF complex luminescent oxygen channeling",
    "immunoassay (LOCI) signal in cps following single IV bolus rFVIII",
    "over a 285-fold dose range (17.5 to 5000 IU/kg). Unbound rFVIII",
    "distributes to a peripheral compartment via first-order rate",
    "constants k12 / k21 and is eliminated by linear clearance CL of",
    "unbound rFVIII. VWF follows zero-order synthesis (ksyn = kdeg * R0",
    "at steady state) and first-order degradation of the unbound",
    "species (kdeg). The rFVIII:VWF complex is eliminated with a",
    "Michaelis-Menten-shape apparent rate constant driven by the",
    "unbound rFVIII concentration",
    "(k_e,comp = Vmax,comp * Cfree / (KM,comp + Cfree)); at Cfree = 0",
    "the complex has no direct elimination and its effective half-life",
    "reduces to that of free VWF (ln(2)/kdeg = ~22 h; Figure 2 of the",
    "paper). The rFVIII:VWF complex plasma concentration maps to the",
    "observed LOCI signal via a 4-parametric logistic",
    "(base_cps + (max_cps - base_cps) * complex^gamma /",
    "(cps50^gamma + complex^gamma)). CL scales allometrically with body",
    "weight at fixed exponent 0.75 (reference weight 0.3 kg = the",
    "cohort median; the paper reports no explicit reference weight).",
    "No covariate on V. Gender and age were tested by the paper but not",
    "retained in the final model."
  )
  reference <- paste(
    "Vargas Christensen I, Loftager M, Rode F, Morck Nielsen H,",
    "Kreilgaard M, Larsen MS.",
    "Impact of capacity-limited binding on recombinant factor VIII",
    "and von Willebrand factor pharmacokinetics in hemophilia A rats.",
    "J Thromb Haemost. 2019 Jun;17(6):964-974.",
    "doi:10.1111/jth.14441. PMID:30924607."
  )
  vignette <- "VargasChristensen_2019_rfviii_rat"
  units <- list(
    time          = "h",
    dosing        = "nmol",
    concentration = "nmol/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    central      = list(analyte = "rFVIII (unbound)", units = "nmol", specimen = "plasma", verified = FALSE),
    peripheral1  = list(analyte = "rFVIII (unbound)", units = "nmol", specimen = "plasma", verified = FALSE),
    total_target = list(analyte = "rFVIII:VWF complex", units = "nmol", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric power scaling on CL of unbound rFVIII with fixed",
        "exponent 0.75 (Vargas Christensen 2019 Methods and Results,",
        "and Table 1). No effect on V. Reference weight is not",
        "explicitly reported in the paper; the median cohort weight",
        "(300 g = 0.3 kg; Methods) is assumed here as the WT/0.3 kg",
        "reference, so the typical CL value of 0.0287 L/h corresponds",
        "to a 0.3 kg rat. Cohort weight range 200-550 g (Methods).",
        "The NONMEM code in Data S1 (not on disk) would nail down the",
        "reference; if it differs, simulated CL scales by",
        "(0.3 kg / actual_ref_kg)^0.75."
      ),
      source_name        = "body weight"
    )
  )

  population <- list(
    species        = "rat (hemophilia A)",
    n_subjects     = 30L,
    n_studies      = 1L,
    age_range      = "8-14 weeks; mean 11 weeks",
    weight_range   = "200-550 g; median 300 g",
    sex_female_pct = 63.3,
    disease_state  = "hemophilia A",
    dose_range     = paste(
      "Single IV bolus rFVIII into tail vein at four dose levels:",
      "17.5 IU/kg (n=8: 5 F / 3 M), 100 IU/kg (n=6: 4 F / 2 M),",
      "1000 IU/kg (n=8: 5 F / 3 M), and 5000 IU/kg (n=8: 5 F / 3 M).",
      "Plasma samples collected over 48 h per a sampling matrix",
      "(Data S1). Samples drawn 2 min post-dose were excluded from",
      "the analysis (Methods)."
    ),
    regions        = "Denmark (Novo Nordisk A/S, Maaloev)",
    notes          = paste(
      "30 hemophilia A rats bred at Novo Nordisk A/S (Maaloev, Denmark)",
      "housed under standard conditions (20-23 C, 30-60% RH,",
      "12/12 light-dark cycle). All animal procedures approved by the",
      "Danish Animal Experiments Council and the Novo Nordisk Ethical",
      "Review Council. Total rFVIII was measured by FVIII:C chromogenic",
      "assay (Chromogenix, Turoctocog alfa SRM 7008 calibrator, LLOQ",
      "120 mU/mL); total VWF was measured by a novel LOCI (calibrated",
      "against literature-derived 14.94 nmol/L human VWF standard as",
      "no rat VWF calibrator was available; LLOQ 0.3 nmol/L); and the",
      "rFVIII:VWF complex was measured by a novel LOCI with no",
      "calibrator available (dependent variable expressed in cps).",
      "The estimated typical R0 (11.5 nmol/L) is lower than the",
      "reported human endogenous VWF (~40 nmol/L); the paper attributes",
      "the discrepancy to the use of a human VWF preparation as",
      "calibrator in a rat assay (Discussion). All parameter estimates",
      "come from Table 1; the NONMEM code lives in Data S1 (not on",
      "disk)."
    )
  )

  ini({
    # Structural PK. Table 1 reports CL = 28.7 mL/h and V = 9.9 mL for the
    # unbound rFVIII disposition; these are converted to L internally so
    # concentrations come out directly in nmol/L (Kss / R0 / KM,comp / cps50
    # are all reported in nmol/L, and all rate constants are 1/h).
    lcl <- log(0.0287);   label("Clearance of unbound rFVIII (L/h)")           # Table 1: CL = 28.7 mL/h -> 0.0287 L/h
    lvc <- log(0.0099);   label("Central volume of distribution (L)")          # Table 1: V  = 9.9  mL -> 0.0099 L

    # Two-compartment micro-rate constants describing UNBOUND rFVIII exchange
    # between central plasma and peripheral tissue. Reported in Table 1 as
    # 1/h rate constants (not as Q + Vp), so we keep the paper's rate-constant
    # parameterisation rather than back-computing Q / Vp.
    lk12 <- log(0.032);   label("Central-to-peripheral rate constant for unbound rFVIII (1/h)")  # Table 1: k12 = 0.032 1/h
    lk21 <- log(0.0488);  label("Peripheral-to-central rate constant for unbound rFVIII (1/h)")  # Table 1: k21 = 0.0488 1/h

    # QSS binding constant of the rFVIII:VWF complex; Kss = k_off / k_on +
    # k_e,FVIII:VWF / k_on in the Mager & Krzyzanski 2005 QSS approximation.
    lkss <- log(0.142);   label("QSS binding constant Kss (nmol/L)")           # Table 1: Kss = 0.142 nmol/L

    # VWF turnover parameters. R0 is the baseline VWF concentration; kdeg is
    # the first-order degradation rate of unbound VWF; ksyn is derived from
    # the steady-state balance ksyn = kdeg * R0 (Methods).
    lrbase <- log(11.5);  label("Baseline total VWF concentration R0 (nmol/L)")           # Table 1: R0   = 11.5  nmol/L
    lkdeg  <- log(0.031); label("Free-VWF first-order degradation rate constant kdeg (1/h)") # Table 1: kdeg = 0.031 1/h

    # rFVIII:VWF complex apparent elimination rate constant is a
    # Michaelis-Menten function of the unbound rFVIII concentration:
    #   k_e,comp = Vmax,comp * Cfree / (KM,comp + Cfree).
    # Vmax,comp is the maximal rate CONSTANT (1/h) at saturating Cfree, not
    # a mass-flux Vmax; KM,comp is the Cfree at which k_e,comp is half of
    # Vmax,comp. At Cfree = 0 (no drug) k_e,comp = 0, so the complex has no
    # direct elimination pathway and its effective loss reduces to the free-
    # VWF turnover kdeg (matches the ~22 h complex half-life at Cfree ~ 0
    # in Figure 2).
    lvmax_comp <- log(0.653); label("Maximal complex-elimination rate constant Vmax,comp (1/h)")          # Table 1: Vmax,comp = 0.653 1/h
    lkm_comp   <- log(0.871); label("Unbound rFVIII conc. at half-max complex elimination KM,comp (nmol/L)") # Table 1: KM,comp   = 0.871 nmol/L

    # Allometric CL exponent held fixed at the theory-based 0.75. Reference
    # weight = 0.3 kg (median of the reported cohort) applied in model()
    # below.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/0.3 kg) on CL")   # Table 1 + Results: Application of a fixed allometric exponent of 0.75 for CL

    # 4-parametric logistic mapping rFVIII:VWF complex plasma concentration
    # to the observed LOCI signal (Equation for LOCI signal in Methods;
    # standard 4-PL with baseline offset):
    #   signal = base_cps + (max_cps - base_cps) *
    #            complex^gamma / (cps50^gamma + complex^gamma).
    # base_cps is the plateau at zero complex; max_cps is the plateau at
    # saturating complex; cps50 is the complex conc at half-max; gamma is
    # the Hill / slope factor.
    lbase_cps <- log(298);    label("Baseline LOCI signal at zero rFVIII:VWF complex, base_cps (cps)")             # Table 1: base_cps = 298 cps
    lmax_cps  <- log(127000); label("Maximal LOCI signal at saturating rFVIII:VWF complex, max_cps (cps)")         # Table 1: max_cps  = 127000 cps
    lcps50    <- log(0.228);  label("rFVIII:VWF complex conc. at half-maximal LOCI signal cps50 (nmol/L)")         # Table 1: cps50    = 0.228 nmol/L
    lhill_cps <- log(2.24);   label("Hill coefficient of the LOCI 4-PL signal function, gamma (unitless)")         # Table 1: gamma    = 2.24

    # Inter-individual variability (log-normal). Table 1 reports IIV as %CV
    # computed as sqrt(exp(omega^2) - 1) * 100; the internal variance is
    # omega^2 = log(1 + CV^2).
    #   CL         IIV 57.1% CV -> omega^2 = log(1 + 0.571^2) = 0.282198
    #   R0         IIV 25.6% CV -> omega^2 = log(1 + 0.256^2) = 0.063478
    #   Vmax,comp  IIV 36.3% CV -> omega^2 = log(1 + 0.363^2) = 0.123782
    etalcl        ~ 0.282198   # Table 1: IIV CL        = 57.1% CV
    etalrbase     ~ 0.063478   # Table 1: IIV R0        = 25.6% CV
    etalvmax_comp ~ 0.123782   # Table 1: IIV Vmax,comp = 36.3% CV

    # Proportional residual error per output. Table 1 reports each error as
    # %CV of a normally-distributed proportional model (Methods); the
    # nlmixr2 prop() form takes the SD as a fraction, dimensionally the
    # same as %CV / 100.
    propSd             <- 0.276; label("Proportional residual SD for total rFVIII (fraction)")                # Table 1: sigma_rFVIII      = 27.6%
    propSd_totalVWF    <- 0.325; label("Proportional residual SD for total VWF (fraction)")                    # Table 1: sigma_VWF         = 32.5%
    propSd_lociComplex <- 0.355; label("Proportional residual SD for rFVIII:VWF complex LOCI signal (fraction)") # Table 1: sigma_rFVIII:VWF  = 35.5%
  })

  model({
    # Individual PK parameters. CL scales allometrically with WT at the
    # fixed 0.75 exponent (reference 0.3 kg = cohort median); no covariate
    # on V.
    cl        <- exp(lcl + etalcl) * (WT / 0.3) ^ e_wt_cl
    vc        <- exp(lvc)
    k12       <- exp(lk12)
    k21       <- exp(lk21)

    # QSS binding + VWF turnover (steady-state derivation ksyn = kdeg * R0).
    kss       <- exp(lkss)
    rbase     <- exp(lrbase + etalrbase)
    kdeg      <- exp(lkdeg)
    ksyn      <- kdeg * rbase

    # rFVIII:VWF complex Michaelis-Menten elimination parameters.
    vmax_comp <- exp(lvmax_comp + etalvmax_comp)
    km_comp   <- exp(lkm_comp)

    # 4-PL LOCI-signal parameters.
    base_cps  <- exp(lbase_cps)
    max_cps   <- exp(lmax_cps)
    cps50     <- exp(lcps50)
    hill_cps  <- exp(lhill_cps)

    # Compartment interpretation:
    #   central       = total rFVIII AMOUNT (nmol) in central plasma
    #                   (unbound + bound as rFVIII:VWF complex).
    #   peripheral1   = unbound rFVIII AMOUNT (nmol) in peripheral tissue.
    #   total_target  = total VWF CONCENTRATION (nmol/L) in central plasma
    #                   (unbound + bound as complex); state IS a concentration
    #                   per the Gibiansky et al. 2008 QSS-TMDD convention.
    # Baseline VWF total = R0 concentration (no drug, no complex).
    total_target(0) <- rbase

    # QSS algebra (Mager & Krzyzanski 2005 / Gibiansky et al. 2008 Eq. 7).
    ctot     <- central / vc
    rtot     <- total_target
    sum_ck   <- ctot + rtot + kss
    complex  <- (sum_ck - sqrt(sum_ck * sum_ck - 4 * ctot * rtot)) / 2
    cfree    <- ctot - complex
    rfree    <- rtot - complex

    # rFVIII:VWF complex apparent elimination rate constant, driven by
    # unbound rFVIII (Michaelis-Menten shape).
    k_e_comp <- vmax_comp * cfree / (km_comp + cfree)

    # ODEs.
    # Central total rFVIII: linear CL on Cfree, complex loss at rate
    # k_e_comp * complex_amount = k_e_comp * complex * vc, and exchange of
    # unbound rFVIII with peripheral (unbound_amount_in_central = cfree * vc).
    d/dt(central)      <- -cl * cfree - k_e_comp * complex * vc -
                            k12 * cfree * vc + k21 * peripheral1
    d/dt(peripheral1)  <-  k12 * cfree * vc - k21 * peripheral1
    # Total VWF: zero-order synthesis - first-order degradation of free VWF -
    # complex loss (the complex takes VWF with it).
    d/dt(total_target) <-  ksyn - kdeg * rfree - k_e_comp * complex

    # Observations.
    # Total rFVIII concentration (nmol/L) as measured by chromogenic assay.
    Cc          <- ctot
    # Total VWF concentration (nmol/L) as measured by VWF LOCI (paper-named
    # observation).
    totalVWF    <- rtot
    # rFVIII:VWF complex LOCI signal (cps) via the 4-parametric logistic;
    # no absolute complex concentration is measured because no rat VWF-FVIII
    # calibrator was available (Methods / Discussion).
    lociComplex <- base_cps + (max_cps - base_cps) *
                   complex ^ hill_cps / (cps50 ^ hill_cps + complex ^ hill_cps)

    # Residual error (all proportional, per Table 1).
    Cc          ~ prop(propSd)
    totalVWF    ~ prop(propSd_totalVWF)
    lociComplex ~ prop(propSd_lociComplex)
  })
}
