Zurlinden_2018_chlorpyrifos <- function() {
  description <- paste(
    "Preclinical (rat; applied cross-species to human). Direct-effect",
    "sigmoidal Emax (Hill) dose-response model relating the peak",
    "chlorpyrifos (CPF) concentration in the brain to a fractional spatial-",
    "memory cognitive deficit (Morris water maze latency-to-platform),",
    "Zurlinden et al. (2018) Equation 1. There is NO pharmacokinetic layer",
    "in this paper: the internal-dose driver (peak brain CPF, umol/L) is",
    "generated externally by the Poet et al. (2014) CPF PBPK model and is",
    "supplied here as the dosed concentration of a single bookkeeping state",
    "(implied volume of 1, no elimination), exactly as the response is a",
    "pure algebraic function of that concentration. Emax is a structural",
    "100% maximum deficit; E50 = 0.035 umol/L and the Hill coefficient",
    "gamma = 0.31 were fit by nonlinear least squares to three once-daily",
    "oral rat doses (Yan et al. 2012). The maximum cognitive deficit is",
    "100%. The 15% deficit point of departure (BMD15) corresponds to a peak",
    "brain CPF of 1.23e-4 umol/L. There is no between-subject variability",
    "and no residual-error model; the reported parameter CVs are the least-",
    "squares estimate uncertainties, not IIV. The paper's biomarker-of-",
    "exposure correlations (Equation 2) and exposure-space boundary",
    "polynomials (Equation 4) are empirical risk-assessment regressions and",
    "are described in the validation vignette rather than shipped as models.",
    sep = " "
  )
  reference <- paste(
    "Zurlinden TJ, Reisfeld B.",
    "A Novel Method for the Development of Environmental Public Health",
    "Indicators and Benchmark Dose Estimation Using a Health-Based End Point",
    "for Chlorpyrifos.",
    "Environ Health Perspect. 2018;126(4):047013.",
    "doi:10.1289/EHP1743.",
    "Internal-dose driver from the Poet et al. (2014) CPF PBPK model;",
    "dose-response calibrated on the Yan et al. (2012) rat Morris water maze",
    "study.",
    sep = " "
  )
  vignette <- "Zurlinden_2018_chlorpyrifos"
  units <- list(
    time = "h",
    dosing = "umol/L",
    concentration = "umol/L"
  )

  # The single ODE state carries the externally-supplied peak brain CPF
  # concentration (Poet et al. 2014 PBPK output), not a fitted PK state. The
  # implied volume is 1 (lvc fixed below) so the dosed amount equals the
  # concentration, mirroring Schulthess_2017_directEffect_stim. d/dt = 0 holds
  # the imposed peak for the whole solve. specimen "tissue" = whole brain.
  compartmentData <- list(
    central = list(
      analyte = "chlorpyrifos",
      units = "umol/L",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "rat (dose-response calibrated on Yan et al. 2012 Morris water maze; the same Equation 1 parameters are applied cross-species to human with human PBPK internal doses)",
    n_subjects = 0L,
    n_studies = 1L,
    disease_state = "not applicable (environmental organophosphate-insecticide neurotoxicity; CPF-induced spatial-memory deficit)",
    dose_range = paste(
      "calibration doses: 0, 1, 5, 10 mg/kg/day once-daily oral gavage CPF",
      "(Yan et al. 2012), giving simulated peak brain CPF of 0, 9.13e-4,",
      "4.62e-3 and 1.02e-2 umol/L (Zurlinden 2018 Table 2); the Emax curve",
      "spans peak brain CPF from ~1e-6 to ~1e-2 umol/L"
    ),
    notes = paste(
      "Zurlinden et al. (2018) is a proof-of-concept framework for deriving",
      "environmental public health indicators (EPHIs) and benchmark doses",
      "(BMDs) for chlorpyrifos using a neurobehavioral (cognitive-deficit)",
      "end point instead of acetylcholinesterase inhibition. Only the",
      "direct-effect dose-response (Equation 1) is encoded here. The",
      "fractional cognitive deficit is the control-normalized reduction in",
      "Morris water maze latency to platform (0.26, 0.34, 0.40 observed for",
      "1, 5, 10 mg/kg/day; Yan et al. 2012). Measurable biomarkers of",
      "exposure in the paper (peak plasma CPF, RBC AChE inhibition, plasma",
      "BuChE inhibition) are related to the deficit through the empirical",
      "Equation 2 polynomial and are not part of this model. No IIV and no",
      "residual error were estimated (deterministic nonlinear least-squares",
      "fit)."
    )
  )

  ini({
    # ================================================================
    # Zurlinden 2018 Equation 1 (sigmoidal Emax / Hill direct effect):
    #   CD = 100 * CB^gamma / (E50^gamma + CB^gamma)
    # CB = peak brain CPF concentration (umol/L); CD = percent cognitive
    # deficit. Parameters: Results 'Dose-response curve' paragraph,
    # E50 = 0.035 (CV 0.47) and gamma = 0.31 (CV 0.18), fit by nonlinear
    # least squares to Table 2 (peak brain CPF) against the Yan et al.
    # (2012) normalized latency-to-platform deficits.
    # ================================================================

    # Implied central volume. The paper has no PK and no volume; the brain
    # CPF concentration is dosed directly into 'central', so vc is fixed at
    # 1 to make Cc == the imposed peak brain CPF concentration.
    lvc <- fixed(log(1)); label("Implied volume of the driver state Vc (umol/L per umol; equals 1, see notes)") # no volume in source; concentration dosed directly

    # Structural maximum effect: Equation 1 is written with a literal
    # '100 x', i.e. a 100% maximum cognitive deficit, not an estimated Emax.
    lemax <- fixed(log(100)); label("Maximum cognitive deficit Emax (percent; structural 100% ceiling)") # Eq 1: the '100 x' factor; 'the maximum effect level represents a 100% deficit in cognitive ability'

    # Estimated Hill parameters (nonlinear least squares; CVs are estimate
    # uncertainty, NOT between-subject variability -- no IIV in this model).
    lec50 <- log(0.035); label("Peak brain CPF producing half-maximal cognitive deficit E50 (umol/L)") # Results 'Dose-response curve': E50 = 0.035 (CV 0.47)
    lhill <- log(0.31); label("Hill / sigmoidicity coefficient gamma (unitless)") # Results 'Dose-response curve': gamma = 0.31 (CV 0.18)
  })

  model({
    vc <- exp(lvc)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # The driver state holds the externally-supplied peak brain CPF
    # concentration (Poet 2014 PBPK). It has no kinetics in this paper, so it
    # is held constant at the dosed level for the whole solve.
    d/dt(central) <- 0
    Cc <- central / vc

    # Equation 1: sigmoidal Emax direct effect. CD is a percent (0-100)
    # cognitive deficit; deterministic (no residual-error endpoint), matching
    # the paper's nonlinear least-squares point model.
    CD <- emax * Cc^hill / (ec50^hill + Cc^hill)
  })
}
